#!/usr/bin/env python3
"""E133 -- the live chain's failing PLP 1 Frames, replayed offline through the same SubframeDecoder.

Input: the windows banked by ATSC3_SF1_DUMP (m11_watch worker: the exact samples and t0 the live
decoder used for every Frame that left a PLP 1 Block unconverged).  For each window, decode:

  live      -- the window as the live chain had it (no band limit), the live t0;
  notch     -- the same window through the standalone loader's band limit (firwin 257 taps,
               2.95 MHz: RF34's 8-VSB pilot sits 3.31 MHz above centre, 7 dB above the wanted
               signal -- m2_pilots.load), the live t0;
  re-t0     -- the raw window with the start re-estimated by m10_core.fine_t0;
  notch+t0  -- both.

The standalone decoder (m17 main) does notch + its own t0 and converges 117/117 over the same
seconds; this says which of the two the live chain is missing.

Run: python lab/e133_sf1_replay.py [--dump C:/captures/sf1_fail] [--threads 4]
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import time

import numpy as np
from scipy.signal import firwin, lfilter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import m6_cells as C                                             # noqa: E402
import m9_fast                                                   # noqa: E402
import m13_sf1 as M13                                            # noqa: E402
import m17_plp1_live as M17                                      # noqa: E402
from m10_core import fine_t0                                     # noqa: E402


def notch(y, fs):
    nt = 257
    b = firwin(nt, 2.95e6, fs=fs)
    yy = np.concatenate([y, np.zeros(nt // 2, y.dtype)])
    return lfilter(b, 1.0, yy)[nt // 2:]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", default=r"C:\captures\sf1_fail")
    ap.add_argument("--l1", default=os.path.join(HERE, "m4_l1detail_hit_rf33.json"))
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args(argv)
    files = sorted(glob.glob(os.path.join(a.dump, "*.npz")))
    if a.limit:
        files = files[:a.limit]
    js = json.load(open(a.l1))
    g0 = C.Geometry.rf33_legacy()
    fd = m9_fast.FrameDecoder(threads=a.threads, backend="cpu", iters=50, cpu_fast=True)
    fd.prewarm()
    dec = M17.SubframeDecoder(fd, g0, js["l1_basic"], js["l1_detail"], 1, comb=True)
    dec.prewarm()
    fs = M13.FS_POST
    print(f"{len(files)} banked Frames; decoder: comb on, threads {a.threads}, iters 50")
    print(f"  {'frame':>5s} {'live said':>10s} {'live':>6s} {'notch':>6s} {'re-t0':>6s} {'both':>6s}   t0 live / re-est   bad(live replay)")
    tot = dict(live=0, notch=0, ret0=0, both=0)
    for f in files:
        z = np.load(f)
        y = np.asarray(z["y"], np.complex128)
        t0 = int(z["t0"])
        t0b, _coh = fine_t0(y, g0)
        yn = notch(y, fs)
        res = {}
        for name, yy, tt in (("live", y, t0), ("notch", yn, t0), ("ret0", y, int(t0b)), ("both", yn, int(t0b))):
            rep = {}
            try:
                dec.decode(yy, tt, rep)
                res[name] = (int(rep["converged"]), rep.get("bad_fail", []))
            except Exception as ex:                                 # noqa: BLE001
                res[name] = (-1, repr(ex)[:60])
            tot[name] += max(res[name][0], 0)
        print(f"  {int(z['idx']):5d} {int(z['converged']):4d}/{int(z['n_fec']):<4d} "
              f"{res['live'][0]:6d} {res['notch'][0]:6d} {res['ret0'][0]:6d} {res['both'][0]:6d}   "
              f"{t0} / {int(t0b)}   {res['live'][1]}")
    print(f"\n  totals over {len(files)} Frames: live {tot['live']}  notch {tot['notch']}  re-t0 {tot['ret0']}  both {tot['both']}"
          f"  (all converged would be {117 * len(files)})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
