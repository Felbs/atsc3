#!/usr/bin/env python3
"""M3 diagnostic -- where along CELL index does the preamble stop being QPSK,
and where along cell index do frames stop agreeing?

Two sharp predictions follow from A/322 if the Step-2 chain is right, and they
are independent of each other and of anything already used to build it.

PREDICTION 1 -- a step at cell 484.
  A/322 7.2.5.2 puts L1-Basic in the FIRST N cells of the first Preamble symbol
  and L1-Detail immediately after.  Table 6.17 fixes N = 484 for L1-Basic
  Mode 3.  L1-Basic Mode 3 is QPSK; L1-Detail is a separately signalled mode
  which for this transmitter is evidently not QPSK.  So a sliding-window QPSK
  MER along cell index must show a CLIFF, and the cliff must sit at 484 -- a
  number we did not measure but looked up.  Nothing in the demodulator knows
  where 484 is; the window is swept blind.

PREDICTION 2 -- a step inside the L1-Basic block, at cell ~184.
  The transmitted L1-Basic block is 968 bits = 484 QPSK cells, and per
  A/322 6.5.2.9 those bits are, in order, the Nouter = 368 systematic bits
  (200 information + 168 BCH parity) followed by the 600 surviving LDPC parity
  bits (12960 parity - 12360 punctured).  368 bits = 184 QPSK cells.
  L1-Basic is near-static frame to frame, so:
      cells   0..183   systematic  -> should agree between frames at ~100%
      cells 184..483   LDPC parity -> flips wholesale if ANY information bit
                                      changes, so ~50% agreement
  A 59% overall agreement is exactly what a mixture like that averages to.  If
  the step is there, at 184, the block boundary predicted by the puncturing
  arithmetic is confirmed from the air -- and the overall 59% stops being a
  failure and becomes a measurement.

Both predictions are falsifiable and neither was used to build the chain.

Usage:
    python m3_profile.py hit_rf33.cs16 --rate 8e6
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
from m3_preamble import analyse, QPSK, mer_db                   # noqa: E402
from m3_replicate import hard_bits                              # noqa: E402


def bar(v, lo, hi, w=44):
    f = (v - lo) / max(hi - lo, 1e-9)
    f = min(max(f, 0.0), 1.0)
    return "#" * int(f * w) + "." * (w - int(f * w))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("capture")
    ap.add_argument("--rate", type=float, required=True)
    ap.add_argument("--fmt")
    ap.add_argument("--frames", type=int, default=8)
    ap.add_argument("--step", type=float, default=0.7)
    ap.add_argument("--win", type=int, default=120)
    ap.add_argument("--json")
    a = ap.parse_args()

    path = a.capture
    if not os.path.isabs(path):
        path = os.path.join(HERE, path)

    xs, bits = [], []
    for i in range(a.frames):
        try:
            rep = {}
            rep, z, Y, geo = analyse(path, a.rate, a.fmt, report=rep,
                                     start_sec=i * a.step, quiet=True)
        except Exception:                                       # noqa: BLE001
            continue
        (lo, n, pilot, cp, data, shift, cred, nfft, gi, dx, mode) = geo
        x = FI.deinterleave(z, nfft, 0, direction="forward", toggle="i")
        xs.append(x)
        bits.append(hard_bits(x))
    if len(xs) < 2:
        raise SystemExit("need >= 2 frames")
    ncells = S.L1_BASIC_CELLS_PRINTED[mode]
    print(f"  {len(xs)} frames demodulated from {os.path.basename(path)}")
    print(f"  A/322 Table 6.17: L1-Basic Mode {mode} = {ncells} cells, "
          f"Nouter {S.L1_BASIC_NOUTER} bits = "
          f"{S.L1_BASIC_NOUTER // 2} cells systematic")

    # ---- PREDICTION 1: sliding QPSK MER along cell index -------------------
    print(f"\n  === PREDICTION 1: QPSK cliff at cell {ncells} ===")
    print(f"  sliding {a.win}-cell window, averaged over {len(xs)} frames\n")
    w = a.win
    cells = min(len(x) for x in xs)
    centres, mers = [], []
    for s in range(0, min(cells, 1600) - w, w // 2):
        vals = [mer_db(x[s:s + w], QPSK) for x in xs]
        centres.append(s + w // 2)
        mers.append(float(np.mean(vals)))
    mlo, mhi = min(mers), max(mers)
    for c, m in zip(centres, mers):
        mark = "  <== A/322 says L1-Basic ends here" if abs(c - ncells) < w // 2 else ""
        print(f"    cells {c-w//2:5d}-{c+w//2:5d}  QPSK MER {m:6.2f} dB  "
              f"{bar(m, mlo, mhi)}{mark}")

    inside = float(np.mean([m for c, m in zip(centres, mers) if c < ncells - w // 2]))
    outside = float(np.mean([m for c, m in zip(centres, mers) if c > ncells + w // 2]))
    print(f"\n    mean QPSK MER before cell {ncells}: {inside:6.2f} dB")
    print(f"    mean QPSK MER after  cell {ncells}: {outside:6.2f} dB")
    print(f"    STEP: {inside - outside:+.2f} dB")

    # ---- PREDICTION 2: frame agreement along cell index --------------------
    print(f"\n  === PREDICTION 2: agreement step at cell "
          f"{S.L1_BASIC_NOUTER // 2} (systematic | LDPC parity) ===")
    B = np.array(bits)                       # frames x (2*ncells_total)
    maj = (B.mean(axis=0) > 0.5).astype(np.uint8)
    agree_bit = np.array([float(np.mean(B[:, j] == maj[j]))
                          for j in range(B.shape[1])])
    # collapse to cells (2 bits per cell)
    agree_cell = agree_bit[:2 * ncells].reshape(-1, 2).mean(axis=1)
    step = 20
    print(f"  agreement with the {len(xs)}-frame majority, {step}-cell bins\n")
    rows = []
    for s in range(0, ncells, step):
        v = float(np.mean(agree_cell[s:s + step]))
        rows.append({"cell": s, "agree": v})
        mark = ("  <== Nouter/2, systematic ends"
                if s <= S.L1_BASIC_NOUTER // 2 < s + step else "")
        print(f"    cells {s:4d}-{min(s+step, ncells)-1:4d}   "
              f"{v*100:6.2f}%  {bar(v, 0.45, 1.0)}{mark}")

    nsys = S.L1_BASIC_NOUTER // 2
    sys_ag = float(np.mean(agree_cell[:nsys]))
    par_ag = float(np.mean(agree_cell[nsys:]))
    print(f"\n    cells   0..{nsys-1:3d} (systematic, 200 info + 168 BCH): "
          f"{sys_ag*100:6.2f}%")
    print(f"    cells {nsys:3d}..{ncells-1:3d} (surviving LDPC parity):      "
          f"{par_ag*100:6.2f}%")
    print(f"    STEP: {(sys_ag - par_ag)*100:+.2f} percentage points")

    out = {"capture": os.path.basename(path), "frames": len(xs),
           "ncells_l1b": ncells, "n_sys_cells": nsys,
           "qpsk_profile": [{"centre": c, "mer": m}
                            for c, m in zip(centres, mers)],
           "qpsk_step_db": inside - outside,
           "agreement_profile": rows,
           "sys_agreement": sys_ag, "parity_agreement": par_ag}
    dest = a.json or os.path.join(
        HERE, "m3_profile_" + os.path.splitext(os.path.basename(path))[0] + ".json")
    with open(dest, "w") as fh:
        json.dump(out, fh, indent=2, default=float)
    print(f"\n  wrote {dest}")


if __name__ == "__main__":
    main()
