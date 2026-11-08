#!/usr/bin/env python3
"""E133 -- the timing plateau of subframe 1: converged Blocks and LDPC iterations against a shift
of the FFT-window start, on banked failing Frames; plus where the cyclic-prefix fold puts the
symbol boundary relative to the live start."""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
import m6_cells as C                                             # noqa: E402
import m9_fast                                                   # noqa: E402
import m13_sf1 as M13                                            # noqa: E402
import m17_plp1_live as M17                                      # noqa: E402
from m2_stage_b import cp_profile, fold_score                    # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", default=r"C:\captures\sf1_fail")
    ap.add_argument("--l1", default=os.path.join(HERE, "m4_l1detail_hit_rf33.json"))
    ap.add_argument("--frames", default="21,15")
    ap.add_argument("--shifts", default="-200,-100,-60,-40,-20,0,20,40,60,100,150,200,300,400,600,800,1000,1200")
    a = ap.parse_args(argv)
    want = [int(x) for x in a.frames.split(",")]
    js = json.load(open(a.l1))
    g0 = C.Geometry.rf33_legacy()
    fd = m9_fast.FrameDecoder(threads=4, backend="cpu", iters=50, cpu_fast=True)
    fd.prewarm()
    dec = M17.SubframeDecoder(fd, g0, js["l1_basic"], js["l1_detail"], 1, comb=True)
    dec.prewarm()
    g = dec.g
    for f in sorted(glob.glob(os.path.join(a.dump, "*.npz"))):
        z = np.load(f)
        if int(z["idx"]) not in want:
            continue
        y = np.asarray(z["y"], np.complex128)
        t0 = int(z["t0"])
        print(f"\nframe {int(z['idx'])}: live t0 {t0}, GI {g.gi}, FFT {g.nfft}")
        # the cyclic-prefix fold over the subframe-1 segment: where is the symbol boundary?
        step = g.nfft + g.gi
        span = g.gi + (g.nsym - 1) * step + g.nfft
        seg = y[t0 + dec.off: t0 + dec.off + span]
        p = cp_profile(seg, g.nfft, g.gi)
        sc, ph, nsy = fold_score(p, step)
        print(f"  CP fold: score {sc:.1f}, peak at {ph} samples into the period (window puts the GI at 0..{g.gi}); "
              f"the fold's best-match region starts {ph} after the live symbol start")
        # the fold around the peak: how wide is the plateau?
        n = (len(p) // step) * step
        fold = p[:n].reshape(-1, step).mean(axis=0)
        thr = fold.max() * 0.9
        above = np.flatnonzero(fold >= thr)
        print(f"  plateau (>= 90 % of peak): {above.min()}..{above.max()} ({above.max() - above.min() + 1} samples)")
        print(f"  {'shift':>6s} {'conv':>5s} {'its':>5s} {'coh':>6s}  bad")
        for sh in [int(x) for x in a.shifts.split(",")]:
            rep = {}
            try:
                dec.decode(y, t0 + sh, rep)
                print(f"  {sh:+6d} {rep['converged']:5d} {rep.get('ldpc_iters_mean', 0):5.1f} "
                      f"{rep.get('data_coherence_mean', 0):6.3f}  {rep.get('bad_fail', [])}")
            except Exception as ex:                                  # noqa: BLE001
                print(f"  {sh:+6d}  err {repr(ex)[:60]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
