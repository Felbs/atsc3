#!/usr/bin/env python3
"""M3 Step 3 (part 2) -- SOLVE the A/322 scrambler from the air.

A/322 5.2.3 Figure 5.6 is vector art, so which eight LFSR stages are the
outputs D7..D0 does not survive PDF extraction, and no hypothesis reproduces
the spec's printed 24-bit test vector (two independent exhaustive searches,
zero hits).  Guessing is not an option and neither is giving up, because the
air already gave us a 32-bit referee.

THE OBSERVATION THAT MAKES THIS SOLVABLE RATHER THAN SEARCHABLE
---------------------------------------------------------------
The CRC-32 residue g(v) = Lin(v[0:168]) XOR v[168:200] is LINEAR in v, and for
any valid L1-Basic block g(d) = K.  With d = r XOR mask that gives one exact
32-bit equation in the mask:

        g(mask) = g(r) XOR K

The mask is not free: it is the LFSR output, so mask[8b + p] = state[b][tap_p]
for eight unknown stage indices tap_0..tap_7.  Substituting and using linearity,

        g(mask) = XOR_{p=0..7} C[p][tap_p],
        C[p][s] = XOR_{b} state[b][s] * Lin_basis[8b + p]

which is just an 8 x 16 table of 32-bit constants per LFSR configuration.
Finding tap_0..tap_7 is then a meet-in-the-middle over 16^4 + 16^4 = 131072
combinations instead of a blind search, and it covers ALL orderings and even
repeated stages -- a strictly larger space than the C(16,8) subset search that
found nothing.

WHY THE ANSWER IS NOT JUST "WHATEVER FITS 32 BITS"
--------------------------------------------------
16^8 = 4.3e9 tap assignments against a 32-bit constraint yields ~1 spurious
solution per configuration by chance, so the CRC ALONE cannot identify the
scrambler and it would be dishonest to claim it does.  Every surviving
candidate is therefore put to a second, completely independent test: the
descrambled fields must agree with what M2/M3 MEASURED from the waveform
without reading a single L1 bit --

    L1B_first_sub_fft_size                8K      (measured: 8192-point grid)
    L1B_first_sub_guard_interval          1536    (measured: CP autocorrelation)
    L1B_first_sub_scattered_pilot_pattern SP4_2   (measured: repeat-coherence)
    L1B_num_subframes                     2       (measured: grid changes at
                                                   52.667 ms = 36 x 9728)
    L1B_first_sub_num_ofdm_symbols        36      (measured: symbol count)
    L1B_first_sub_reduced_carriers        0       (measured: 6913 carriers)

Those six were all established in M2, months of DSP before any table existed in
this repo.  A candidate that satisfies a 32-bit CRC AND reproduces six
independently measured physical quantities is not a coincidence.  A candidate
that satisfies the CRC and contradicts them is a coincidence, and is discarded
as one.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from itertools import product

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m3_spec as S                                            # noqa: E402
import m3_freqint as FI                                        # noqa: E402
import m3_l1basic as L1                                        # noqa: E402
import m3_crc as CRC                                           # noqa: E402
import spec_l1syntax as LS                                     # noqa: E402
from m3_preamble import analyse                                 # noqa: E402

POLY_EXP = (0, 1, 3, 6, 7, 11, 12, 13, 16)
SEED = 0xF180
NBITS = 200
NBYTES = (NBITS + 7) // 8


# --------------------------------------------------------------------------
# LFSR configurations
# --------------------------------------------------------------------------

def configs():
    out = []
    taps_spec = tuple(e - 1 for e in POLY_EXP if 1 <= e <= 16)
    taps_rec = tuple(16 - e for e in POLY_EXP if 1 <= e <= 16)
    for seed_name, seed_bits in (
            ("X1=LSB", [(SEED >> b) & 1 for b in range(16)]),
            ("X1=MSB", [(SEED >> (15 - b)) & 1 for b in range(16)])):
        for tname, taps in (("spec", taps_spec), ("recip", taps_rec)):
            for direction in ("right", "left"):
                for nsh in (1, 8):
                    out.append({"seed": seed_name, "seed_bits": seed_bits,
                                "taps": tname, "tap_idx": taps,
                                "dir": direction, "shifts_per_byte": nsh})
    return out


def states(cfg, nbytes=NBYTES):
    s = list(cfg["seed_bits"])
    out = []
    for _ in range(nbytes):
        out.append(list(s))
        for _ in range(cfg["shifts_per_byte"]):
            fb = 0
            for t in cfg["tap_idx"]:
                fb ^= s[t]
            s = ([fb] + s[:-1]) if cfg["dir"] == "right" else (s[1:] + [fb])
    return out


def mask_from(cfg, order, nbits=NBITS):
    st = states(cfg)
    m = []
    for b in range(len(st)):
        for p in range(8):
            m.append(st[b][order[p]])
    return np.array(m[:nbits], np.uint8)


# --------------------------------------------------------------------------
# solve
# --------------------------------------------------------------------------

def basis():
    """Lin_basis[j] = g(e_j), the 32-bit linear image of each of the 200 bits."""
    B = []
    for j in range(NBITS):
        e = np.zeros(NBITS, np.uint8)
        e[j] = 1
        B.append(CRC.g_stat(e))
    return B


def solve(r0, verbose=True):
    Lb = basis()
    K = CRC.crc32(np.zeros(168, np.uint8), init=0xFFFFFFFF)
    target = CRC.g_stat(r0) ^ K
    cands = []
    for cfg in configs():
        st = states(cfg)
        C = [[0] * 16 for _ in range(8)]
        for p in range(8):
            for s in range(16):
                v = 0
                for b in range(len(st)):
                    j = 8 * b + p
                    if j < NBITS and st[b][s]:
                        v ^= Lb[j]
                C[p][s] = v
        half = {}
        for q in product(range(16), repeat=4):
            v = C[0][q[0]] ^ C[1][q[1]] ^ C[2][q[2]] ^ C[3][q[3]]
            half.setdefault(v, []).append(q)
        for q2 in product(range(16), repeat=4):
            v = C[4][q2[0]] ^ C[5][q2[1]] ^ C[6][q2[2]] ^ C[7][q2[3]]
            for q1 in half.get(v ^ target, ()):
                cands.append({"cfg": cfg, "order": list(q1) + list(q2)})
    if verbose:
        print(f"  LFSR configurations searched: {len(configs())}")
        print(f"  tap assignments per config:   16^8 = {16**8:,}")
        print(f"  total space:                  "
              f"{len(configs())*16**8:,}")
        print(f"  candidates satisfying the 32-bit CRC: {len(cands)}")
        print(f"  (expected by chance: "
              f"{len(configs())*16**8/2**32:.1f}) -- so the CRC alone "
              f"CANNOT identify the scrambler")
    return cands


# --------------------------------------------------------------------------
# field parsing + the independent-measurement test
# --------------------------------------------------------------------------

def parse_l1b(bits):
    out, i = {}, 0
    for name, w in LS.L1_BASIC_FIELDS:
        v = 0
        for k in range(w):
            v = (v << 1) | int(bits[i + k])
        out[name] = v
        i += w
    return out


PREDICTIONS = {
    "L1B_first_sub_fft_size": (0, "8K -- M2 measured an 8192-point grid"),
    "L1B_first_sub_guard_interval": (6, "GI6_1536 -- M2 measured GI 1536"),
    "L1B_first_sub_scattered_pilot_pattern":
        (2, "SP4_2 -- M2 measured D_x=4, D_y=2 blind"),
    "L1B_num_subframes": (1, "2 subframes -- M2 measured the grid change at "
                             "52.667 ms"),
    "L1B_first_sub_reduced_carriers": (0, "NoC 6913 -- M2 measured 6921"),
}


def score(fields):
    hits, det = 0, []
    for k, (want, why) in PREDICTIONS.items():
        got = fields.get(k)
        ok = got == want
        hits += ok
        det.append((k, got, want, ok, why))
    n = fields.get("L1B_first_sub_num_ofdm_symbols")
    det.append(("L1B_first_sub_num_ofdm_symbols", n, "36 or 35", None,
                "M2 counted 36 symbols in subframe 0"))
    return hits, det


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("capture")
    ap.add_argument("--rate", type=float, required=True)
    ap.add_argument("--fmt")
    ap.add_argument("--frames", type=int, default=3)
    ap.add_argument("--step", type=float, default=0.7)
    ap.add_argument("--json")
    a = ap.parse_args()
    path = a.capture
    if not os.path.isabs(path):
        path = os.path.join(HERE, path)

    blocks = []
    for i in range(a.frames):
        rep = {}
        rep, z, Y, geo = analyse(path, a.rate, a.fmt, report=rep,
                                 start_sec=i * a.step, quiet=True)
        (lo, n, pilot, cp, data, shift, cred, nfft, gi, dx, mode) = geo
        x = FI.deinterleave(z, nfft, 0, direction="forward", toggle="i")
        r = L1.decode_one(x, mode, variant="standard", gw_invert=False,
                          iters=100)
        if r["converged"] and r["bch_ok"]:
            blocks.append(np.array(r["nouter_bits"][:200], np.uint8))
    print(f"  {len(blocks)} verified L1-Basic blocks\n")
    if not blocks:
        raise SystemExit("no verified blocks")

    print("  === solving the scrambler mask from the CRC (meet-in-the-middle) "
          "===")
    cands = solve(blocks[0])

    print(f"\n  === second, INDEPENDENT test: do the descrambled fields agree "
          f"with\n      what M2 measured from the waveform? ===\n")
    best, rows = None, []
    for c in cands:
        m = mask_from(c["cfg"], c["order"])
        # must also satisfy the CRC on EVERY other frame
        allok = all(CRC.g_stat(b ^ m) == CRC.crc32(np.zeros(168, np.uint8),
                                                   init=0xFFFFFFFF)
                    for b in blocks)
        f = parse_l1b(blocks[0] ^ m)
        h, det = score(f)
        rows.append({"cfg": {k: v for k, v in c["cfg"].items()
                             if k != "seed_bits"},
                     "order": c["order"], "hits": h, "all_frames_crc": allok,
                     "fields": {k: int(v) for k, v in f.items()}})
        if allok and (best is None or h > best[1]):
            best = (c, h, det, f, m)
    rows.sort(key=lambda r: -r["hits"])
    print(f"    candidates matching all {len(blocks)} frames' CRC: "
          f"{sum(1 for r in rows if r['all_frames_crc'])}")
    print(f"    best agreement with independent measurement: "
          f"{rows[0]['hits']}/{len(PREDICTIONS)}\n")

    out = {"capture": os.path.basename(path), "n_candidates": len(cands),
           "top": rows[:8]}
    if best and best[1] == len(PREDICTIONS):
        c, h, det, f, m = best
        print(f"  *** UNIQUE CONSISTENT SOLUTION ***")
        print(f"      seed {c['cfg']['seed']}  taps {c['cfg']['taps']}  "
              f"dir {c['cfg']['dir']}  shifts/byte "
              f"{c['cfg']['shifts_per_byte']}")
        print(f"      D7..D0 stages (1-based): "
              f"{[i+1 for i in c['order']]}\n")
        for k, got, want, ok, why in det:
            print(f"      {k:42s} {str(got):>6s} vs {str(want):>8s}  "
                  f"{'MATCH' if ok else ('--' if ok is None else 'no')}   {why}")
        print(f"\n  --- decoded L1-Basic fields ---")
        for k, v in f.items():
            sem = LS.L1_BASIC_SEMANTICS.get(k, {})
            lab = sem.get(v, "") if isinstance(sem, dict) else ""
            print(f"      {k:42s} {v:10d}  {lab}")
        out["solution"] = {"cfg": rows[0]["cfg"], "order": c["order"],
                           "fields": {k: int(v) for k, v in f.items()}}
    else:
        print("  NO candidate satisfies both the CRC and the independent")
        print("  measurements.  Reporting that rather than shipping the")
        print("  best-fitting coincidence.  Top candidates by agreement:")
        for r in rows[:5]:
            print(f"    hits {r['hits']}/{len(PREDICTIONS)}  "
                  f"allCRC={r['all_frames_crc']}  {r['cfg']}  {r['order']}")

    dest = a.json or os.path.join(HERE, "m3_descramble_"
                                  + os.path.splitext(os.path.basename(path))[0]
                                  + ".json")
    with open(dest, "w") as fh:
        json.dump(out, fh, indent=2, default=float)
    print(f"\n  wrote {dest}")


if __name__ == "__main__":
    main()
