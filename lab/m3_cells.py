#!/usr/bin/env python3
"""M3 Step 2 (part 2) -- de-interleave the preamble and isolate L1-Basic.

A/322 7.2.5.2: "L1-Basic cells shall be mapped only to the available cells of
the first Preamble symbol ... L1-Detail cells shall be interleaved and mapped
to the remaining available cells of the first Preamble symbol directly after
the L1-Basic cells".

So in CELL order the first Preamble symbol is:

    [ 0 .. Ncells_L1B-1 ]  L1-Basic, always QPSK for L1-Basic Modes 1-3
    [ Ncells_L1B .. ]      L1-Detail, whose constellation is a DIFFERENT and
                           independently signalled mode (Table 6.17)

Cell order is only reachable through the A/322 7.3 frequency de-interleaver.
That makes this a sharp, falsifiable prediction rather than a fit: a specific
FI-determined subset of 484 carriers -- 10% of the symbol, scattered all over
it -- must be clean QPSK while the complement need not be.  Picking the
"best 484 of 4851" would of course always find a tight subset; the point is
that we do not get to pick.  The subset is dictated by the spec.

CONTROLS run alongside the legal hypothesis:
  * both wire-permutation directions (ambiguity F1) -- only one can be right;
  * a RANDOM 484-cell subset, the same size, which is what "we just found a
    tight cluster because we looked" would score;
  * the LAST 484 cells in cell order, which under the spec are L1-Detail.

Usage:
    python m3_cells.py hit_rf33.cs16 --rate 8e6
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
from m3_preamble import (analyse, QPSK, mer_db, m4_stat,        # noqa: E402
                         qam_points, ascii_scatter)


def describe(z, label, rep=None):
    orders = [("QPSK", QPSK), ("16QAM", qam_points(16)),
              ("64QAM", qam_points(64)), ("256QAM", qam_points(256))]
    m = {n: mer_db(z, p) for n, p in orders}
    m4 = m4_stat(z)
    best = max(m, key=lambda k: m[k])
    print(f"    {label:34s} n={len(z):5d}  "
          + "  ".join(f"{n}:{m[n]:6.2f}" for n in m)
          + f"   |E[z^4]|={abs(m4):.3f}@{np.degrees(np.angle(m4)):+7.1f}deg")
    if rep is not None:
        rep[label] = {"n": int(len(z)), "mer": m,
                      "m4_abs": abs(m4),
                      "m4_deg": float(np.degrees(np.angle(m4)))}
    return m, m4


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("capture")
    ap.add_argument("--rate", type=float, required=True)
    ap.add_argument("--fmt")
    ap.add_argument("--json")
    a = ap.parse_args()

    print("A/322 spec-table verification:")
    r = S.verify(verbose=False)
    print(f"  {sum(1 for _, g, _ in r if g)}/{len(r)} cross-table checks pass")
    if not all(g for _, g, _ in r):
        raise SystemExit("spec tables FAILED -- stopping")

    path = a.capture
    if not os.path.isabs(path):
        path = os.path.join(HERE, path)

    rep = {}
    rep, z, Y, geo = analyse(path, a.rate, a.fmt, report=rep)
    (lo, n, pilot, cp, data, shift, cred, nfft, gi, dx, mode) = geo

    ncells = S.L1_BASIC_CELLS_PRINTED[mode]
    print(f"\n  === A/322 7.3 frequency de-interleave, then isolate the "
          f"first {ncells} cells ===")
    print(f"  (A/322 Table 6.17: L1-Basic Mode {mode} occupies {ncells} cells, "
          f"{S.L1_BASIC_CONSTELLATION[mode]})\n")
    print("    subset                             count   "
          "best-fit MER by order (dB)                 4th-power statistic")

    rng = np.random.default_rng(20260805)
    results = {}
    sub = {}
    for direction in FI.WIRE_DIRECTIONS:
        x = FI.deinterleave(z, nfft, 0, direction=direction, toggle="i")
        m, _ = describe(x[:ncells], f"L1-Basic cells 0..{ncells-1} [{direction}]",
                        results)
        sub[direction] = x
        results[f"_mer_qpsk_{direction}"] = m["QPSK"]

    print()
    best_dir = max(FI.WIRE_DIRECTIONS,
                   key=lambda d: results[f"_mer_qpsk_{d}"])
    x = sub[best_dir]
    describe(x[ncells:], f"L1-Detail cells {ncells}.. [{best_dir}]", results)
    describe(x[-ncells:], f"LAST {ncells} cells (control)", results)
    describe(z[rng.choice(len(z), ncells, replace=False)],
             f"random {ncells} cells (control)", results)
    describe(z, "all cells, carrier order", results)

    rep["freq_interleaver"] = results
    rep["wire_direction"] = best_dir
    rep["l1b_cells"] = ncells

    print(f"\n  --- L1-Basic constellation ({best_dir} wire permutation) ---")
    print(ascii_scatter(x[:ncells], w=57, h=25))

    mq = results[f"L1-Basic cells 0..{ncells-1} [{best_dir}]"]["mer"]["QPSK"]
    mr = results[f"random {ncells} cells (control)"]["mer"]["QPSK"]
    print(f"\n  VERDICT: L1-Basic QPSK MER {mq:.2f} dB vs "
          f"random-subset control {mr:.2f} dB "
          f"({mq - mr:+.2f} dB)")

    out = a.json or os.path.join(
        HERE, "m3_cells_" + os.path.splitext(os.path.basename(path))[0] + ".json")
    with open(out, "w") as fh:
        json.dump(rep, fh, indent=2, default=float)
    print(f"  wrote {out}")


if __name__ == "__main__":
    main()
