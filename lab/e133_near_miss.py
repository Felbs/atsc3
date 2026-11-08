#!/usr/bin/env python3
"""E133 -- what finishes a near miss?  On the banked failed batches (ATSC3_LDPC_DUMP_FAIL from
m9_fast.ldpc_batch: the LLRs of every batch that left a Block unconverged, with the verdicts), try
the decoder's levers on the FAILED Blocks only and count what each converges, bytes-checked by the
parity check itself (a converged Block has a zero syndrome by definition):

  * the shipping ladder at 50 sweeps (the reference: should converge nothing new);
  * more sweeps: 100, 200, 400 at each rung;
  * other alphas: 0.95, 0.9, 0.8, 0.7, 0.6 at 50 and 200 sweeps;
  * the two-minimum correction off;
  * the float32 kernel vs the int16 kernel;
  * the flooding kernel (min_sum_decode_batch) at 200.

Run: python lab/e133_near_miss.py [--dump C:/captures/plp1_fail] [--rate 11/15]
"""
from __future__ import annotations

import argparse
import glob
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import m3_ldpc as LD                                             # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", default=r"C:\captures\plp1_fail")
    ap.add_argument("--rate", default="11/15")
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args(argv)
    files = sorted(glob.glob(os.path.join(a.dump, "*.npz")))
    lams, meta = [], []
    for f in files:
        z = np.load(f)
        bad_i = np.flatnonzero(~z["conv"])
        if not len(bad_i):
            continue
        lams.append(np.ascontiguousarray(z["lam"][bad_i], np.float32))
        meta.append((os.path.basename(f), [int(x) for x in z["bad"][bad_i]], [int(x) for x in z["its"][bad_i]]))
    if not lams:
        print("no failed blocks banked")
        return 1
    n = int(np.load(files[0])["ninner"])
    ladder = tuple(float(x) for x in np.load(files[0])["ladder"])
    lam = np.concatenate(lams)
    print(f"{len(files)} failed batches, {len(lam)} failed Blocks of {n} bits, rate {a.rate}, shipping ladder {ladder}")
    for m in meta:
        print(f"  {m[0]}: bad {m[1]} its {m[2]}")
    checks = LD.parity_check(a.rate, n=n)
    ex = ThreadPoolExecutor(max_workers=a.threads)

    def run(name, iters, alphas, corr=True, q16="0", flood=False):
        LD.Q16_DEFAULT = q16
        todo = np.arange(len(lam))
        conv_all = np.zeros(len(lam), bool)
        t = time.perf_counter()
        for al in alphas:
            if not len(todo):
                break
            sub = np.ascontiguousarray(lam[todo])
            if flood:
                r = LD.min_sum_decode_batch(sub, checks, n, iters=iters, alpha=al, dtype=np.float32)
            else:
                r = LD.layered_decode_batch(sub, checks, n, iters=iters, alpha=al, ex=ex, corr=corr,
                                            abort=(0, 0, 0))
            c = np.asarray(r[1], bool)
            conv_all[todo[c]] = True
            todo = todo[~c]
        dt = time.perf_counter() - t
        print(f"  {name:44s} converged {int(conv_all.sum()):3d}/{len(lam)}   {dt:6.1f} s")
        return conv_all

    print("\nlevers on the failed Blocks only (abort rule off so every sweep is spent):")
    run("ladder 50 (reference)", 50, ladder)
    run("ladder 100", 100, ladder)
    run("ladder 200", 200, ladder)
    run("ladder 400", 400, ladder)
    for al in (0.95, 0.9, 0.8, 0.7, 0.6):
        run(f"alpha {al} x 50", 50, (al,))
    for al in (0.95, 0.9, 0.8, 0.7, 0.6):
        run(f"alpha {al} x 200", 200, (al,))
    run("ladder 200, correction OFF", 200, ladder, corr=False)
    run("ladder 200, int16 kernel", 200, ladder, q16="1")
    run("ladder 100, flooding kernel", 100, ladder, flood=True)
    run("fine ladder 1.0,.95,.9,.85,.8,.75,.7 x 100", 100, (1.0, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7))
    return 0


if __name__ == "__main__":
    sys.exit(main())
