#!/usr/bin/env python3
"""M31 -- A-SPX subband group tables: the frequency skeleton of the high band.

M28 reads `aspx_config` and M28's gate 4 already showed the A-SPX range starts
at exactly the QMF subband where the MDCT stops (12000 Hz).  This builds the
tables that every A-SPX structure is indexed by -- the signal envelope groups
at high and low resolution, and the noise envelope groups -- because
`aspx_framing`, `aspx_hfgen_iwc_*` and `aspx_ec_data` all loop over
`num_sbg_sig_highres`, `num_sbg_sig_lowres` and `num_sbg_noise`.  Without these
counts the A-SPX payload cannot even be *skipped*, let alone decoded.

WHAT IS DERIVED, AND FROM WHAT (clause 5.7.6.3.1)
--------------------------------------------------
Everything comes from three header fields plus one field at the head of
`aspx_data`:

    aspx_master_freq_scale  picks the static template          (aspx_config)
    aspx_start_freq         cuts the low end of the master     (aspx_config)
    aspx_stop_freq          cuts the high end of the master    (aspx_config)
    aspx_noise_sbg          sets the noise group count         (aspx_config)
    aspx_xover_subband_offset  splits master -> signal range   (aspx_data, 3 b)

    Pseudocode 67  sbg_master   = template[2*start_freq + sbg]
    Pseudocode 68  sbg_sig_highres = master truncated at the crossover offset
                   sbx = first subband, num_sb_aspx = span
    Pseudocode 69  sbg_sig_lowres  = an exact decimation of highres by 2,
                                     anchored differently for odd vs even counts
    Pseudocode 70  sbg_noise       = a further thinning of lowres, with
                                     num_sbg_noise from a log2 span formula

THE GATES
----------
These tables are pure arithmetic on values read from the stream, so the useful
checks are structural properties the spec states independently:

  1. `num_sbg_noise <= 5` -- stated as a shall in clause 5.7.6.3.1.3.
  2. every table strictly increasing, and nested inside the master range.
  3. the master table starts AND ends on an even QMF subband -- clause
     5.7.6.3.1.2 says so explicitly, and it is not something the arithmetic
     forces, so it tests the template transcription and the index arithmetic.
  4. sbg_sig_lowres is a genuine subsequence of sbg_sig_highres, and its group
     count follows the floor rule.
  5. the A-SPX range starts at the MDCT's coded edge (M28 gate 4, repeated here
     against the derived sbx rather than the raw start frequency).

Usage:
    python m31_aspx_bands.py
"""
from __future__ import annotations

import argparse
import collections
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m17_ac4_walk as W                                          # noqa: E402
import m20_ac4_toc2 as M                                          # noqa: E402
import m27_sfb as B                                               # noqa: E402
import m28_channels as C                                          # noqa: E402
from m19_ac4_toc import Bits                                      # noqa: E402


# Table 191, frame_length 1536 at 48 kHz.  A 64-band QMF turns a 1536-sample
# frame into 24 QMF timeslots, and two QMF slots make one A-SPX slot.
NUM_TS_IN_ATS = {2048: 2, 1920: 2, 1536: 2, 1024: 1, 960: 1, 768: 1,
                 512: 1, 384: 1}
TS_OFFSET_HFGEN = {2048: 6, 1920: 6, 1536: 6, 1024: 3, 960: 3, 768: 3,
                   512: 3, 384: 3}


def num_aspx_timeslots(frame_length=1536, qmf_bands=64):
    """Pseudocode 75a.  Decides the 2-vs-1 bit widths in aspx_framing.

    aspx_framing's Note 1 says several border fields shrink to 1 bit when
    num_aspx_timeslots <= 8.  For this stream it is 12, so they stay at 2 --
    a width that would otherwise have to be guessed.
    """
    return (frame_length // qmf_bands) // NUM_TS_IN_ATS[frame_length]


def sbg_master(cfg):
    """Pseudocode 67."""
    tpl = (C.SBG_TEMPLATE_HIGHRES if cfg["master_freq_scale"] == 1
           else C.SBG_TEMPLATE_LOWRES)
    n = (22 if cfg["master_freq_scale"] == 1 else 20) \
        - 2 * cfg["start_freq"] - 2 * cfg["stop_freq"]
    base = 2 * cfg["start_freq"]
    return [tpl[base + i] for i in range(n + 1)]


def sbg_tables(cfg, xover):
    """Pseudocodes 68-70.  -> dict of every derived table and count."""
    master = sbg_master(cfg)
    n_master = len(master) - 1

    n_hi = n_master - xover                              # Pseudocode 68
    hi = [master[s + xover] for s in range(n_hi + 1)]
    sbx = hi[0]
    num_sb_aspx = hi[n_hi] - sbx

    n_lo = n_hi - n_hi // 2                              # Pseudocode 69
    lo = [hi[0]]
    for s in range(1, n_lo + 1):
        lo.append(hi[2 * s] if n_hi % 2 == 0 else hi[2 * s - 1])

    sba, sbz = master[0], master[n_master]               # Pseudocode 70
    n_noise = max(1, math.floor(cfg["noise_sbg"] * math.log2(sbz / sba) + 0.5))
    idx = [0]
    noise = [lo[0]]
    for s in range(1, n_noise + 1):
        idx.append(idx[s - 1]
                   + (n_lo - idx[s - 1]) // (n_noise + 1 - s))
        noise.append(lo[idx[s]])

    t = dict(master=master, n_master=n_master, sba=sba, sbz=sbz,
             hi=hi, n_hi=n_hi, lo=lo, n_lo=n_lo,
             noise=noise, n_noise=n_noise,
             sbx=sbx, num_sb_aspx=num_sb_aspx)
    t.update(patches(t, cfg))
    return t


def patches(t, cfg, base_samp_freq=48):
    """Pseudocode 71 -- the HF patch table.

    A-SPX does not transmit the high band; it COPIES low subbands upward and
    then reshapes them.  This decides which source block feeds which part of
    the A-SPX range, and it is the last purely structural piece before the
    high band can be synthesised.

    The loop is transcribed literally, including the `odd` parity term and the
    `sbg_master[sbg] - sb < 3` escape that jumps the search back to the top of
    the master table.  Both look arbitrary and both change the answer.
    """
    master, n_master = t["master"], t["n_master"]
    sba, sbx, span = t["sba"], t["sbx"], t["num_sb_aspx"]
    msb, usb = sba, sbx
    goal_sb = 43 if base_samp_freq == 48 else 46
    source_band_low = 4 if cfg["master_freq_scale"] == 1 else 2

    if goal_sb < sbx + span:
        sbg = 0
        for i in range(len(master)):
            if master[i] >= goal_sb:
                break
            sbg = i + 1
    else:
        sbg = n_master

    num_sb, start_sb = [], []
    guard = 0
    while True:
        guard += 1
        if guard > 64:                       # the spec's do/while has no bound
            raise RuntimeError("patch loop did not terminate")
        j = sbg
        sb = master[j]
        odd = (sb - 2 + sba) % 2
        while sb > (sba - source_band_low + msb - odd):
            j -= 1
            sb = master[j]
            odd = (sb - 2 + sba) % 2
        n = max(sb - usb, 0)
        num_sb.append(n)
        start_sb.append(sba - odd - n)
        if n > 0:
            usb = msb = sb
        else:
            msb = sbx
            num_sb.pop()
            start_sb.pop()
        if master[sbg] - sb < 3:
            sbg = n_master
        if sb == sbx + span:
            break

    if len(num_sb) > 1 and num_sb[-1] < 3:
        num_sb.pop()
        start_sb.pop()

    borders = [sbx]
    for n in num_sb:
        borders.append(borders[-1] + n)
    return dict(patch_num_sb=num_sb, patch_start_sb=start_sb,
                patches=borders, n_patches=len(num_sb))


def read_xover(f, T):
    """The 3-bit aspx_xover_subband_offset at the head of aspx_data_2ch.

    Only present in i-frames.  It sits immediately after the six channel
    elements, so decoding the element gives its exact bit position -- no
    search, no guess.
    """
    st = M.parse(f)
    if not st["b_iframe_global"]:
        return None, None
    o = st["toc_bytes"] + st["substream_sizes"][0]
    sub = f[o:o + st["substream_sizes"][1]]
    r = C.decode_element(sub, T, b_iframe=True)
    if r["bits"] > r["nbits"]:
        return None, None
    b = Bits(sub)
    b.p = r["bits"]
    return b.u(3), r["aspx"]


def main(argv=None):
    ap = argparse.ArgumentParser(prog="m31 aspx bands")
    ap.add_argument("path", nargs="?",
                    default="m7_out/rf33_audio_pid13.mp4")
    a = ap.parse_args(argv)
    p = a.path if os.path.isabs(a.path) else os.path.join(HERE, a.path)

    print("M31 -- A-SPX subband group tables")
    print("=" * 74)
    T = C.load_tables()
    fr = W.samples(p)

    seen = collections.Counter()
    cfg0 = None
    for f in fr:
        try:
            x, cfg = read_xover(f, T)
        except Exception:                                      # noqa: BLE001
            continue
        if x is None:
            continue
        seen[x] += 1
        cfg0 = cfg
    print(f"\n  aspx_xover_subband_offset over {sum(seen.values())} i-frames: "
          f"{dict(seen)}")
    if not cfg0:
        print("  no i-frames decoded")
        return 1
    xover = seen.most_common(1)[0][0]
    print(f"  aspx_config: start_freq {cfg0['start_freq']} stop_freq "
          f"{cfg0['stop_freq']} master_freq_scale {cfg0['master_freq_scale']} "
          f"noise_sbg {cfg0['noise_sbg']}")

    t = sbg_tables(cfg0, xover)
    hz = lambda sb: sb * C.QMF_HZ                              # noqa: E731
    print(f"\n  master  ({t['n_master']} groups)  {t['master']}")
    print(f"          {hz(t['sba']):.0f} .. {hz(t['sbz']):.0f} Hz")
    print(f"  sig hi  ({t['n_hi']} groups)  {t['hi']}")
    print(f"  sig lo  ({t['n_lo']} groups)  {t['lo']}")
    print(f"  noise   ({t['n_noise']} groups)  {t['noise']}")
    print(f"  sbx {t['sbx']} = {hz(t['sbx']):.0f} Hz, "
          f"num_sb_aspx {t['num_sb_aspx']} subbands "
          f"= {t['num_sb_aspx'] * C.QMF_HZ:.0f} Hz wide")
    nats = num_aspx_timeslots()
    print(f"\n  num_qmf_timeslots {1536 // 64}, num_ts_in_ats "
          f"{NUM_TS_IN_ATS[1536]} -> num_aspx_timeslots {nats}")
    print(f"  {nats} > 8, so aspx_framing's relative-border fields are "
          f"2 bits, not 1 (Note 1)")

    print(f"\n  HF patches ({t['n_patches']}): borders {t['patches']}")
    for i, (n, st) in enumerate(zip(t["patch_num_sb"], t["patch_start_sb"])):
        d0, d1 = t["patches"][i] * C.QMF_HZ, t["patches"][i + 1] * C.QMF_HZ
        print(f"    patch {i}: {n} subbands sourced from {st}..{st + n} "
              f"({st * C.QMF_HZ:.0f}..{(st + n) * C.QMF_HZ:.0f} Hz)"
              f"  ->  {d0:.0f}..{d1:.0f} Hz")

    print("\n  GATES")
    mono = all(all(v[i] < v[i + 1] for i in range(len(v) - 1))
               for v in (t["master"], t["hi"], t["lo"], t["noise"]))
    nested = (t["hi"][0] >= t["master"][0]
              and t["hi"][-1] <= t["master"][-1]
              and t["noise"][0] >= t["lo"][0]
              and t["noise"][-1] <= t["lo"][-1])
    even = t["master"][0] % 2 == 0 and t["master"][-1] % 2 == 0
    sub = set(t["lo"]).issubset(set(t["hi"]))
    lo_rule = t["n_lo"] == t["n_hi"] - t["n_hi"] // 2
    n_ok = t["n_noise"] <= 5
    mdct = B.SFB_OFFSET_1536[43] / 1536.0 * 24000.0
    edge = abs(hz(t["sbx"]) - mdct) < 1.0
    for label, ok in (
            (f"num_sbg_noise = {t['n_noise']} <= 5  (clause 5.7.6.3.1.3 shall)",
             n_ok),
            ("every table strictly increasing", mono),
            ("signal range nested in master, noise nested in lowres", nested),
            ("master starts AND ends on an even QMF subband", even),
            ("sbg_sig_lowres is a subsequence of sbg_sig_highres", sub),
            (f"num_sbg_sig_lowres follows the floor rule", lo_rule),
            (f"A-SPX starts at the MDCT's coded edge ({mdct:.0f} Hz)", edge)):
        print(f"    {'PASS' if ok else 'FAIL'}  {label}")

    p_ok = t["n_patches"] <= 5
    p_tile = t["patches"][0] == t["sbx"] and         t["patches"][-1] == t["sbx"] + t["num_sb_aspx"]
    p_src = all(st >= 0 and st + n <= t["sbx"]
                for st, n in zip(t["patch_start_sb"], t["patch_num_sb"]))
    for label, ok in (
            (f"num_sbg_patches = {t['n_patches']} <= 5  (5.7.6.3.1.4 shall)",
             p_ok),
            ("patches tile the A-SPX range exactly", p_tile),
            ("every patch sources from BELOW the crossover", p_src)):
        print(f"    {'PASS' if ok else 'FAIL'}  {label}")

    good = all((n_ok, mono, nested, even, sub, lo_rule, edge,
                p_ok, p_tile, p_src))
    print("\n" + "=" * 74)
    print("  A-SPX FREQUENCY SKELETON DERIVED.  These counts are what "
          "aspx_framing,\n  aspx_hfgen_iwc and aspx_ec_data are indexed by."
          if good else "  NOT established -- see the failing gate")
    return 0 if good else 1


if __name__ == "__main__":
    raise SystemExit(main())
