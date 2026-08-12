#!/usr/bin/env python3
"""M19 -- the complete AC-4 table of contents, parsed, with a gate that decides.

Everything before this read the AC-4 header as far as it could be checked
without the standard.  The standard is now on disk, the three unknowns M19's
notes named are resolved, and this is the full `ac4_toc()` implemented from
clause 4.2.3.

THE GATE IS NOT AN OPINION
---------------------------
`raw_ac4_frame()` is `ac4_toc()`, byte-aligned, then one
`ac4_substream_data()` per substream.  The TOC declares each substream's SIZE.
The container declares the FRAME's size.  So:

    sum(substream_size) + len(toc) + padding  ==  frame length

on every frame, or the parse is wrong.  There is no way to be nearly right:
one misread bit shifts every subsequent field and the sizes become nonsense.

A SECOND, INDEPENDENT GATE, WRITTEN BEFORE THE SPEC ARRIVED
------------------------------------------------------------
M18's entropy map found FOUR bits that are identical to `b_iframe_global` on
all 3480 frames, at positions 23, 95, 108 and 147 -- with no spec at all.  The
parse must put `b_iframe` bits exactly there.  That map was built when it could
not be influenced by any reading of the standard, which is what makes it a real
check rather than a restatement.

THE THREE UNKNOWNS, RESOLVED
-----------------------------
  * protection_length_primary   (Table 175) 00=Reserved 01=8 10=32 11=128
    protection_length_secondary (Table 176) 00=0        01=8 10=32 11=128
  * channel_mode (Table 88) is a prefix code: 0=Mono, 10=Stereo,
    1100/1101/1110 = 3.0/5.0/5.1, 1111000..1111101 = the 7.x layouts,
    1111110 reserved, 1111111 escapes to variable_bits(2).
  * frame_rate_factor (Table 87): for frame_rate_index 2/3/4 it is 1 when
    b_multiplier is 0, else 2 or 4 by multiplier_bit.  It sets how many
    `b_iframe` bits each ac4_substream_info() carries -- which is where the
    extra I-frame bits in the entropy map should come from.

Usage:
    python m19_ac4_toc.py [m7_out/rf33_audio_pid13.mp4]
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

PROT_PRIMARY = {0: None, 1: 8, 2: 32, 3: 128}     # 0 is Reserved
PROT_SECOND = {0: 0, 1: 8, 2: 32, 3: 128}
CH_MODE_NAME = {0: "Mono", 1: "Stereo", 2: "3.0", 3: "5.0", 4: "5.1",
                5: "7.0 (3/4/0)", 6: "7.1 (3/4/0.1)", 7: "7.0 (5/2/0)",
                8: "7.1 (5/2/0.1)", 9: "7.0 (3/2/2)", 10: "7.1 (3/2/2.1)",
                11: "reserved"}


class Bits:
    """MSB-first reader that RAISES at the end rather than returning silence."""

    def __init__(self, data):
        self.d, self.p = data, 0

    def u(self, n):
        if (self.p + n) > len(self.d) * 8:
            raise EOFError(f"ran out at bit {self.p} wanting {n}")
        v = 0
        for _ in range(n):
            v = (v << 1) | ((self.d[self.p >> 3] >> (7 - (self.p & 7))) & 1)
            self.p += 1
        return v

    def vb(self, n):
        """variable_bits(n), Table 3 -- offset added on each continuation."""
        val = 0
        while True:
            val += self.u(n)
            if not self.u(1):
                return val
            val <<= n
            val += (1 << n)

    def align(self):
        while self.p & 7:
            self.p += 1


def channel_mode(b):
    """Table 88's prefix code -> (ch_mode index, bits consumed)."""
    if b.u(1) == 0:
        return 0, 1
    if b.u(1) == 0:
        return 1, 2
    two = b.u(2)                       # we are past '11'
    if two != 0b11:
        return {0b00: 2, 0b01: 3, 0b10: 4}[two], 4
    tail = b.u(3)                      # '1111' + 3
    if tail == 0b111:
        return 12 + b.vb(2), 7
    return 5 + tail, 7


def emdf_info(b, out):
    ver = b.u(2)
    if ver == 3:
        ver += b.vb(2)
    key = b.u(3)
    if key == 7:
        key += b.vb(3)
    if b.u(1):                         # b_emdf_payloads_substream_info
        if b.u(2) == 3:                # substream_index
            b.vb(2)
    prim = b.u(2)
    sec = b.u(2)
    if PROT_PRIMARY[prim] is None:
        raise ValueError("protection_length_primary == 0 is Reserved")
    b.u(PROT_PRIMARY[prim])
    if PROT_SECOND[sec]:
        b.u(PROT_SECOND[sec])
    out["emdf"] += 1


def substream_info(b, fs_index, frf, out):
    cm, _ = channel_mode(b)
    out["channel_modes"].append(cm)
    if fs_index == 1 and b.u(1):       # b_sf_multiplier
        b.u(1)
    if b.u(1):                         # b_bitrate_info
        # 3 bits, extended to 5 when the low 3 are all ones
        if b.u(3) == 0b111:
            b.u(2)
    if cm in (7, 8, 9, 10):
        b.u(1)                         # add_ch_base
    if b.u(1):                         # b_content_type
        content_type(b)
    for _ in range(frf):
        out["iframe_bits"].append(b.p)
        b.u(1)
    if b.u(2) == 3:                    # substream_index
        b.vb(2)


def content_type(b):
    b.u(3)                             # content_classifier
    if b.u(1):                         # b_language_indicator
        if b.u(1):                     # b_serialized_language_tag
            b.u(1)
            b.u(16)
        else:
            n = b.u(6)
            for _ in range(n):
                b.u(8)


def presentation_info(b, frame_rate_index, fs_index, out):
    single = b.u(1)
    cfg = 0
    if not single:
        cfg = b.u(3)
        if cfg == 7:
            cfg += b.vb(2)
    while b.u(1):                      # presentation_version, unary
        pass
    add_emdf = 0
    if (not single) and cfg == 6:
        add_emdf = 1
    else:
        b.u(3)                         # mdcompat
        if b.u(1):                     # b_belongs_to_presentation_id
            b.vb(2)
        # frame_rate_multiply_info
        frf = 1
        if frame_rate_index in (2, 3, 4):
            if b.u(1):
                frf = 4 if b.u(1) else 2
        elif frame_rate_index in (0, 1, 7, 8, 9):
            if b.u(1):
                frf = 2
        out["frame_rate_factor"] = frf
        emdf_info(b, out)
        if single:
            substream_info(b, fs_index, frf, out)
        else:
            hsf = b.u(1)
            n_main = {0: 2, 1: 2, 2: 2, 3: 3, 4: 3, 5: 1}.get(cfg)
            if n_main is None:
                raise ValueError(f"presentation_config {cfg} not implemented")
            substream_info(b, fs_index, frf, out)
            if hsf:
                if b.u(2) == 3:        # ac4_hsf_ext_substream_info
                    b.vb(2)
            for _ in range(n_main - 1):
                substream_info(b, fs_index, frf, out)
        b.u(1)                         # b_pre_virtualized
        add_emdf = b.u(1)
    if add_emdf:
        n = b.u(2)
        if n == 0:
            n = b.vb(2) + 4
        for _ in range(n):
            emdf_info(b, out)
    out["presentation_config"] = cfg
    out["single_substream"] = single


def parse_toc(frame):
    b = Bits(frame)
    out = collections.defaultdict(int)
    out["channel_modes"] = []
    out["iframe_bits"] = []
    bv = b.u(2)
    if bv == 3:
        bv += b.vb(2)
    out["bitstream_version"] = bv
    out["sequence_counter"] = b.u(10)
    if b.u(1):                         # b_wait_frames
        if b.u(3) > 0:
            b.u(2)                     # reserved
    fs = b.u(1)
    fri = b.u(4)
    out["fs_index"], out["frame_rate_index"] = fs, fri
    out["iframe_global_bit"] = b.p
    out["b_iframe_global"] = b.u(1)
    n_pres = 1 if b.u(1) else 0
    if not n_pres:
        n_pres = b.vb(2) + 2 if b.u(1) else 0
    out["n_presentations"] = n_pres
    payload_base = 0
    if b.u(1):
        payload_base = b.u(5) + 1
        if payload_base == 0x20:
            payload_base += b.vb(3)
    out["payload_base"] = payload_base
    for _ in range(n_pres):
        presentation_info(b, fri, fs, out)
    # substream_index_table
    n_sub = b.u(2)
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
    out["n_substreams"] = n_sub
    out["substream_sizes"] = sizes
    out["toc_bytes"] = b.p // 8
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(prog="m19 ac4 toc")
    ap.add_argument("path", nargs="?",
                    default="m7_out/rf33_audio_pid13.mp4")
    a = ap.parse_args(argv)
    p = a.path if os.path.isabs(a.path) else os.path.join(HERE, a.path)
    fr = W.samples(p)
    print("M19 -- the complete AC-4 TOC")
    print("=" * 72)
    print(f"  {len(fr)} frames")

    ok_sum = err = 0
    iframe_pos = collections.Counter()
    cfg = collections.Counter()
    chm = collections.Counter()
    frf = collections.Counter()
    nsub = collections.Counter()
    first_err = None
    for i, f in enumerate(fr):
        try:
            t = parse_toc(f)
        except Exception as e:                                 # noqa: BLE001
            err += 1
            if first_err is None:
                first_err = (i, f"{type(e).__name__}: {e}")
            continue
        if sum(t["substream_sizes"]) + t["toc_bytes"] == len(f):
            ok_sum += 1
        cfg[(t["single_substream"], t["presentation_config"])] += 1
        frf[t["frame_rate_factor"]] += 1
        nsub[t["n_substreams"]] += 1
        for c in t["channel_modes"]:
            chm[c] += 1
        iframe_pos[tuple([t["iframe_global_bit"]] + t["iframe_bits"])] += 1

    print(f"\n  parsed without error   {len(fr) - err}/{len(fr)}")
    if first_err:
        print(f"    first failure: frame {first_err[0]}: {first_err[1]}")
    print(f"  GATE 1  sizes + toc == frame length: {ok_sum}/{len(fr)}")
    print(f"\n  n_substreams        {dict(nsub)}")
    print(f"  frame_rate_factor   {dict(frf)}")
    print(f"  (single, config)    {dict(cfg)}")
    print(f"  channel_mode        " +
          str({CH_MODE_NAME.get(k, k): v for k, v in chm.items()}))
    print(f"\n  GATE 2  I-frame bit positions (global + per substream):")
    for pos, n in iframe_pos.most_common(3):
        print(f"    {list(pos)}  on {n} frames")
    want = [23, 95, 108, 147]
    hit = any(list(p) == want for p in iframe_pos)
    print(f"    M18's entropy map (built with NO spec) said {want}"
          f"  -> {'MATCH' if hit else 'does not match'}")
    print("\n" + "=" * 72)
    good = ok_sum == len(fr) and hit
    print("  BOTH GATES PASS -- the TOC is fully parsed." if good else
          "  NOT YET -- see which gate failed above.")
    return 0 if good else 1


if __name__ == "__main__":
    raise SystemExit(main())
