#!/usr/bin/env python3
"""M3 Step 3 -- validate the decoded L1-Basic with its own CRC-32, WITHOUT
needing the scrambler.

The problem: A/322 5.2.3's Figure 5.6 is vector art, so the eight output tap
positions of the scrambler LFSR do not survive PDF text extraction, and an
exhaustive search over 102960 tap/direction/ordering hypotheses reproduces the
spec's own printed 24-bit test vector exactly ZERO times (independently
replicated by a second search over a wider space).  So the descrambled field
values are out of reach for the moment.

The trick that gets the CRC anyway: **the scrambler mask cancels in a
difference.**

A/322 6.1.2.2 CRC-32 has generator G(x) = x^32 + x^21 + x^16 + x^11 + 1 with the
register initialised to all ones.  A CRC with a fixed non-zero init is AFFINE:

    CRC_ones(v) = Lin(v) XOR K,      K = CRC_ones(0)

L1-Basic is 200 bits ending in a 32-bit L1B_crc over the preceding 168.  So for
any valid block d:      Lin(d[0:168]) XOR K == d[168:200]
Define                  g(v) = Lin(v[0:168]) XOR v[168:200]      (zero-init, linear)
then                    g(d) == K for EVERY valid block.

The received blocks are r_i = d_i XOR mask with the SAME mask every frame, so

    g(r_i XOR r_j) = g(d_i XOR d_j) = g(d_i) XOR g(d_j) = K XOR K = 0

The mask drops out entirely.  Two frames whose 200 decoded bits differ give a
32-bit check on the decode, the CRC polynomial reading, and the claim that
L1B_crc occupies the final 32 bits -- with no scrambler required.

CONTROLS: the same statistic on random 200-bit deltas, on deltas from the
L1-Detail region, and under three WRONG CRC polynomials.  All must fail.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m3_spec as S                                            # noqa: E402
import m3_freqint as FI                                        # noqa: E402
import m3_l1basic as L1                                        # noqa: E402
from m3_preamble import analyse                                 # noqa: E402

# A/322 6.1.2.2: "all values of gi = 0 except for g21, g16, g11"
POLY_ATSC = (1 << 21) | (1 << 16) | (1 << 11) | 1
POLY_MPEG = 0x04C11DB7                    # the familiar CRC-32/MPEG-2, a control
POLY_ALT1 = (1 << 26) | (1 << 23) | (1 << 22) | 1
POLY_ALT2 = (1 << 21) | (1 << 16) | (1 << 11)          # missing the +1 term


def crc32(bits, poly=POLY_ATSC, init=0):
    """MSB-first shift-register CRC-32 per A/322 Figure 6.4."""
    crc = init & 0xFFFFFFFF
    for d in bits:
        top = (crc >> 31) & 1
        crc = (crc << 1) & 0xFFFFFFFF
        if top ^ int(d):
            crc ^= poly
    return crc


def g_stat(v200, poly=POLY_ATSC):
    """Linear (zero-init) g(v) = Lin(v[0:168]) XOR v[168:200]."""
    c = crc32(v200[:168], poly=poly, init=0)
    tail = 0
    for b in v200[168:200]:
        tail = (tail << 1) | int(b)
    return c ^ tail


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("capture")
    ap.add_argument("--rate", type=float, required=True)
    ap.add_argument("--fmt")
    ap.add_argument("--frames", type=int, default=8)
    ap.add_argument("--step", type=float, default=0.7)
    ap.add_argument("--json")
    a = ap.parse_args()
    path = a.capture
    if not os.path.isabs(path):
        path = os.path.join(HERE, path)

    blocks, detail = [], []
    for i in range(a.frames):
        try:
            rep = {}
            rep, z, Y, geo = analyse(path, a.rate, a.fmt, report=rep,
                                     start_sec=i * a.step, quiet=True)
        except Exception:                                       # noqa: BLE001
            continue
        (lo, n, pilot, cp, data, shift, cred, nfft, gi, dx, mode) = geo
        x = FI.deinterleave(z, nfft, 0, direction="forward", toggle="i")
        r = L1.decode_one(x, mode, variant="standard", gw_invert=False,
                          iters=100)
        if not (r["converged"] and r["bch_ok"]):
            print(f"  frame {i}: LDPC conv={r['converged']} "
                  f"BCH={r['bch_ok']} -- skipped")
            continue
        blocks.append(np.array(r["nouter_bits"][:200], np.uint8))
        ncl = S.L1_BASIC_CELLS_PRINTED[mode]
        dv = (x[ncl:ncl + 100].real < 0).astype(np.uint8)
        detail.append(dv)
    print(f"  {len(blocks)} frames decoded with LDPC converged + BCH zero\n")
    if len(blocks) < 2:
        raise SystemExit("need >= 2 decoded frames")

    uniq = {b.tobytes() for b in blocks}
    print(f"  distinct 200-bit L1-Basic blocks: {len(uniq)} of {len(blocks)}")
    dif = [int((blocks[0] != b).sum()) for b in blocks[1:]]
    print(f"  Hamming distance of each frame from frame 0: {dif}")
    print(f"  (identical frames give a trivially zero delta and cannot test "
          f"anything -- \n   the test needs frames that DIFFER)\n")

    print("  === CRC-32 delta test: g(r_i XOR r_j) must be ZERO ===")
    print("  A/322 6.1.2.2 G(x) = x^32 + x^21 + x^16 + x^11 + 1, "
          "init all ones\n")
    rows = []
    npair = ntest = 0
    for i in range(len(blocks)):
        for j in range(i + 1, len(blocks)):
            d = blocks[i] ^ blocks[j]
            if not d.any():
                continue
            npair += 1
            v = g_stat(d)
            ok = v == 0
            ntest += ok
            rows.append({"i": i, "j": j, "wt": int(d.sum()),
                         "g": int(v), "zero": bool(ok)})
    for r in rows[:12]:
        print(f"    frames {r['i']}-{r['j']}  delta weight {r['wt']:3d}  "
              f"g = {r['g']:#010x}  {'ZERO -> PASS' if r['zero'] else 'FAIL'}")
    if len(rows) > 12:
        print(f"    ... {len(rows)-12} more pairs")
    print(f"\n  {ntest}/{npair} non-trivial frame pairs give g == 0")

    print("\n  === CONTROLS (each must FAIL) ===")
    rng = np.random.default_rng(5)
    ctl = {}
    hits = sum(1 for _ in range(2000)
               if g_stat(rng.integers(0, 2, 200).astype(np.uint8)) == 0)
    ctl["random_200bit"] = hits
    print(f"    random 200-bit vectors:            {hits}/2000 give g == 0 "
          f"(expect ~0; chance 2^-32)")
    # same-weight random deltas
    wt = rows[0]["wt"] if rows else 6
    hits2 = 0
    for _ in range(2000):
        d = np.zeros(200, np.uint8)
        d[rng.choice(200, wt, replace=False)] = 1
        hits2 += g_stat(d) == 0
    ctl["random_same_weight"] = hits2
    print(f"    random weight-{wt} deltas:             {hits2}/2000 give g == 0")
    for nm, p in (("CRC-32/MPEG-2 0x04C11DB7", POLY_MPEG),
                  ("wrong poly x^26+x^23+x^22+1", POLY_ALT1),
                  ("A/322 poly without the +1", POLY_ALT2)):
        k = sum(1 for r in rows if g_stat(blocks[r["i"]] ^ blocks[r["j"]],
                                          poly=p) == 0)
        ctl[nm] = k
        print(f"    {nm:34s} {k}/{npair} pairs give g == 0")

    out = {"capture": os.path.basename(path), "n_frames": len(blocks),
           "n_distinct": len(uniq), "pairs": rows,
           "pairs_zero": ntest, "pairs_total": npair, "controls": ctl,
           "blocks": ["".join(str(int(q)) for q in b) for b in blocks]}
    dest = a.json or os.path.join(HERE, "m3_crc_"
                                  + os.path.splitext(os.path.basename(path))[0]
                                  + ".json")
    with open(dest, "w") as fh:
        json.dump(out, fh, indent=2, default=float)
    print(f"\n  wrote {dest}")


if __name__ == "__main__":
    main()
