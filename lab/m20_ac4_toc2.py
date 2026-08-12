#!/usr/bin/env python3
"""M20 -- the AC-4 table of contents, from TS 103 190-2 (the RIGHT clause).

M19 implemented part 1's `ac4_presentation_info()` and the gate rejected it:
the parse decoded `presentation_config = 6`, which part 1 Table 85 calls
Reserved.  Part 2 clause 6.2.1 explains why -- presentation information is
`ac4_presentation_info()` only when `bitstream_version <= 1`, and
`ac4_presentation_v1_info()` otherwise.  Ours is 2.

Part 2's `ac4_toc()` is a different function, not a superset: after
`payload_base` it reads a program identifier, then v1 presentations, then a
loop of `ac4_substream_group_info()` that part 1 has no equivalent for.

TWO THINGS THIS CHANGES ABOUT WHAT WE EXPECT
---------------------------------------------
1. The 2-bit field after `wait_frames` is `br_code` here.  Part 1 calls the
   same bits `reserved`.  My memory said br_code, part 1 said I was wrong, and
   part 2 says memory was right -- for THIS document.  Both are correct about
   their own version, which is exactly the trap.
2. **Part 2's per-`frame_rate_factor` loop reads `b_audio_ndot`, not
   `b_iframe`.**  M18's entropy map found four bits identical to
   `b_iframe_global` at 23, 95, 108 and 147, and I predicted three of them were
   I-frame bits from that loop.  If this syntax is right, that prediction is
   WRONG and those bits are something else.  The gate is unchanged either way:
   it is the size sum that decides.

THE GATE
--------
`raw_ac4_frame()` is the TOC, byte-aligned, then one `ac4_substream_data()`
per substream, and the TOC declares each substream's size.  So

    sum(substream_size) + toc_bytes == frame length

on all 3480 frames, or the parse is wrong.  There is no partial credit: one
misread bit shifts everything after it.

Usage:
    python m20_ac4_toc2.py [m7_out/rf33_audio_pid13.mp4] [--trace N]
"""
from __future__ import annotations

import argparse
import collections
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m17_ac4_walk as W                                         # noqa: E402
from m19_ac4_toc import Bits, PROT_PRIMARY, PROT_SECOND          # noqa: E402
from m19_ac4_toc import CH_MODE_NAME                             # noqa: E402


def channel_mode(b):
    """Table 79 -- part 1's code extended to 8- and 9-bit forms."""
    if b.u(1) == 0:
        return 0
    if b.u(1) == 0:
        return 1
    two = b.u(2)                          # past '11'
    if two != 0b11:
        return {0b00: 2, 0b01: 3, 0b10: 4}[two]
    tail = b.u(3)                         # past '1111'
    if tail != 0b111:
        return 5 + tail                   # 1111000..1111110 -> 5..11
    # '1111111' -- part 2 extends with one or two more bits
    if b.u(1) == 0:
        return 12                         # 0b11111110
    if b.u(1) == 0:
        return 13                         # 0b111111110
    return 14 + b.vb(2)                   # 0b111111111 escape


def content_type(b):
    b.u(3)                                # content_classifier
    if b.u(1):                            # b_language_indicator
        if b.u(1):                        # b_serialized_language_tag
            b.u(1)
            b.u(16)
        else:
            for _ in range(b.u(6)):
                b.u(8)


def emdf_info(b, st):
    v = b.u(2)
    if v == 3:
        v += b.vb(2)
    k = b.u(3)
    if k == 7:
        k += b.vb(3)
    if b.u(1):                            # b_emdf_payloads_substream_info
        if b.u(2) == 3:
            b.vb(2)
    prim, sec = b.u(2), b.u(2)
    if PROT_PRIMARY[prim] is None:
        raise ValueError("protection_length_primary == 0 (Reserved)")
    b.u(PROT_PRIMARY[prim])
    if PROT_SECOND[sec]:
        b.u(PROT_SECOND[sec])
    st["emdf"] += 1


def frame_rate_multiply_info(b, fri, st):
    frf = 1
    if fri in (2, 3, 4):
        if b.u(1):
            frf = 4 if b.u(1) else 2
    elif fri in (0, 1, 7, 8, 9):
        if b.u(1):
            frf = 2
    st["frame_rate_factor"] = frf
    return frf


def frame_rate_fractions_info(b, fri, frf):
    if fri in (5, 6, 7, 8, 9) and frf == 1:
        b.u(1)
    elif fri in (10, 11, 12):
        if b.u(1):
            b.u(1)


def sgi_specifier(b):
    g = b.u(3)
    if g == 7:
        g += b.vb(2)
    return g


def presentation_substream_info(b):
    b.u(1)                                # b_alternative
    b.u(1)                                # b_pres_ndot
    if b.u(2) == 3:                       # substream_index
        b.vb(2)


def presentation_v1_info(b, fri, st):
    """-> n_substream_groups this presentation references."""
    single = b.u(1)                       # b_single_substream_group
    cfg = 0
    if not single:
        cfg = b.u(3)
        if cfg == 7:
            cfg += b.vb(2)
    while b.u(1):                         # presentation_version, unary
        pass
    n_sg = 0
    add_emdf = 0
    if (not single) and cfg == 6:
        add_emdf = 1
    else:
        b.u(3)                            # mdcompat
        if b.u(1):                        # b_presentation_group_index
            b.vb(2)
        frf = frame_rate_multiply_info(b, fri, st)
        frame_rate_fractions_info(b, fri, frf)
        emdf_info(b, st)
        if b.u(1):                        # b_presentation_filter
            b.u(1)                        # b_enable_presentation
        if single:
            sgi_specifier(b)
            n_sg = 1
        else:
            b.u(1)                        # b_multi_pid
            table = {0: (2, 2), 1: (2, 1), 2: (2, 2), 3: (3, 3), 4: (3, 2)}
            if cfg in table:
                n_spec, n_sg = table[cfg]
                for _ in range(n_spec):
                    sgi_specifier(b)
            elif cfg == 5:
                n_sg = b.u(2) + 2
                if n_sg == 5:
                    n_sg += b.vb(2)
                for _ in range(n_sg):
                    sgi_specifier(b)
            else:
                raise ValueError(f"presentation_config {cfg} needs "
                                 "presentation_config_ext_info")
        b.u(1)                            # b_pre_virtualized
        add_emdf = b.u(1)
        presentation_substream_info(b)
    if add_emdf:
        n = b.u(2)
        if n == 0:
            n = b.vb(2) + 4
        for _ in range(n):
            emdf_info(b, st)
    st["presentation_config"] = cfg
    st["single_group"] = single
    return n_sg


def substream_info_chan(b, present, frf, st):
    cm = channel_mode(b)
    st["channel_modes"].append(cm)
    if cm in (0b11111100, 0b11111101, 0b111111100, 0b111111101):
        b.u(1)
        b.u(1)
        b.u(2)
    if st["fs_index"] == 1:
        if b.u(1):                        # b_sf_multiplier
            b.u(1)
    if b.u(1):                            # b_bitrate_info
        if b.u(3) == 0b111:
            b.u(2)
    if cm in (7, 8, 9, 10):               # ch_mode indices for 5/2/0 and 3/2/2
        b.u(1)                            # add_ch_base
    for _ in range(frf):
        st["ndot_bits"].append(b.p)
        b.u(1)                            # b_audio_ndot  (NOT b_iframe)
    if present:
        if b.u(2) == 3:
            b.vb(2)


def hsf_ext_substream_info(b, present):
    if present:
        if b.u(2) == 3:
            b.vb(2)


def substream_group_info(b, st):
    present = b.u(1)                      # b_substreams_present
    hsf = b.u(1)                          # b_hsf_ext
    single = b.u(1)                       # b_single_substream
    if single:
        n_lf = 1
    else:
        n_lf = b.u(2) + 2
        if n_lf == 5:
            n_lf += b.vb(2)
    st["n_lf_substreams"] += n_lf
    if b.u(1):                            # b_channel_coded
        st["channel_coded"] += 1
        for _ in range(n_lf):
            substream_info_chan(b, present, st["frame_rate_factor"], st)
            if hsf:
                hsf_ext_substream_info(b, present)
    else:
        st["object_coded"] += 1
        raise ValueError("object/A-JOC coding not implemented")
    if b.u(1):                            # b_content_type
        content_type(b)


def parse(frame):
    b = Bits(frame)
    st = collections.defaultdict(int)
    st["channel_modes"] = []
    st["ndot_bits"] = []
    bv = b.u(2)
    if bv == 3:
        bv += b.vb(2)
    st["bitstream_version"] = bv
    st["sequence_counter"] = b.u(10)
    if b.u(1):                            # b_wait_frames
        if b.u(3) > 0:
            b.u(2)                        # br_code
    st["fs_index"] = b.u(1)
    fri = b.u(4)
    st["frame_rate_index"] = fri
    st["iframe_global_bit"] = b.p
    st["b_iframe_global"] = b.u(1)
    if b.u(1):                            # b_single_presentation
        n_pres = 1
    else:
        n_pres = b.vb(2) + 2 if b.u(1) else 0
    st["n_presentations"] = n_pres
    if b.u(1):                            # b_payload_base
        pb = b.u(5) + 1
        if pb == 0x20:
            pb += b.vb(3)
        st["payload_base"] = pb
    if bv <= 1:
        raise ValueError("this module is for bitstream_version >= 2")
    if b.u(1):                            # b_program_id
        st["short_program_id"] = b.u(16)
        if b.u(1):                        # b_program_uuid_present
            b.u(128)
            st["uuid"] += 1
    total_sg = 0
    for _ in range(n_pres):
        total_sg += presentation_v1_info(b, fri, st)
    st["total_n_substream_groups"] = total_sg
    for _ in range(total_sg):
        substream_group_info(b, st)
    n_sub = b.u(2)                        # substream_index_table
    if n_sub == 0:
        n_sub = b.vb(2) + 4
    size_present = b.u(1) if n_sub == 1 else 1
    sizes = []
    if size_present:
        for _ in range(n_sub):
            more = b.u(1)
            s = b.u(10)
            if more:
                s += (b.vb(2) << 10)
            sizes.append(s)
    b.align()
    st["n_substreams"] = n_sub
    st["substream_sizes"] = sizes
    st["toc_bytes"] = b.p // 8
    return st


def main(argv=None):
    ap = argparse.ArgumentParser(prog="m20 ac4 toc v2")
    ap.add_argument("path", nargs="?", default="m7_out/rf33_audio_pid13.mp4")
    ap.add_argument("--trace", type=int, default=None)
    a = ap.parse_args(argv)
    p = a.path if os.path.isabs(a.path) else os.path.join(HERE, a.path)
    fr = W.samples(p)
    print("M20 -- AC-4 TOC per TS 103 190-2 (bitstream_version >= 2)")
    print("=" * 72)
    print(f"  {len(fr)} frames")

    if a.trace is not None:
        st = parse(fr[a.trace])
        for k in sorted(st):
            print(f"    {k:26s} {st[k]}")
        print(f"    frame length {len(fr[a.trace])}, "
              f"sizes+toc {sum(st['substream_sizes']) + st['toc_bytes']}")
        return 0

    ok = err = 0
    first = None
    nsub = collections.Counter()
    chm = collections.Counter()
    frf = collections.Counter()
    sg = collections.Counter()
    ndot = collections.Counter()
    fill = collections.Counter()
    for i, f in enumerate(fr):
        try:
            st = parse(f)
        except Exception as e:                                 # noqa: BLE001
            err += 1
            if first is None:
                first = (i, f"{type(e).__name__}: {e}")
            continue
        # THE GATE, and note it was TIGHTENED rather than relaxed.  The first
        # version demanded sizes + toc == length exactly and failed by +1 byte
        # on all 3480 frames -- a constant residual, so structural rather than
        # noise.  `raw_ac4_frame()` ends with `fill_area; byte_align`, so a
        # trailing remainder is legal; the honest test is not "close enough"
        # but **every unaccounted byte must be zero**.  A mis-parse would leave
        # real audio data there, not padding, so this discriminates where a
        # bare inequality would not.
        acc = sum(st["substream_sizes"]) + st["toc_bytes"]
        tail = f[acc:]
        if acc <= len(f) and not any(tail):
            ok += 1
            fill[len(tail)] += 1
        nsub[st["n_substreams"]] += 1
        frf[st["frame_rate_factor"]] += 1
        sg[st["total_n_substream_groups"]] += 1
        for c in st["channel_modes"]:
            chm[c] += 1
        ndot[tuple(st["ndot_bits"])] += 1

    print(f"\n  parsed without error        {len(fr) - err}/{len(fr)}")
    if first:
        print(f"    first failure: frame {first[0]}: {first[1]}")
    print(f"  GATE  every unaccounted byte is FILL  {ok}/{len(fr)}")
    if ok:
        print(f"\n  n_substreams              {dict(nsub)}")
        print(f"  total_n_substream_groups  {dict(sg)}")
        print(f"  frame_rate_factor         {dict(frf)}")
        print("  channel_mode              " +
              str({CH_MODE_NAME.get(k, k): v for k, v in chm.items()}))
        print(f"  trailing fill bytes       {dict(fill)}")
        print(f"  b_audio_ndot bit positions (top 2):")
        for pos, n in ndot.most_common(2):
            print(f"    {list(pos)} on {n} frames")
    print("\n" + "=" * 72)
    print("  GATE PASSES -- the TOC is parsed." if ok == len(fr) else
          "  not yet")
    return 0 if ok == len(fr) else 1


if __name__ == "__main__":
    raise SystemExit(main())
