#!/usr/bin/env python3
"""size_ldpc_abort.py -- E112's sizing, re-run on any banked LDPC record.

Ground truth = the layered kernel with NO abort (every block to the cap).
Each candidate rule is then judged on the same LLRs: blocks LOST (converged
under truth, not under the rule -- an abort can only lose), lane-sweeps saved,
and the vector-sweeps the eight-at-a-time schedule actually pays.  13.14 did
this once on the Linux box's marginal capture and chose k=8, T=0.36 x nc; this
is the same question asked of a different regime (uniform added noise, where
every failing block is a near miss).

Run: python lab/size_ldpc_abort.py --dump DIR [--rate 9/15] [--threads 8]
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
import m3_ldpc as LD                                                # noqa: E402


def vec_sweeps(its):
    return sum(int(its[lo:lo + 8].max()) for lo in range(0, len(its), 8))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", required=True)
    ap.add_argument("--rate", default="9/15")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--corr", type=int, default=1)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    files = [f for f in sorted(glob.glob(os.path.join(a.dump, "*.npz"))) if np.load(f)["lam"].size]
    if a.limit:
        files = files[:a.limit]
    z0 = np.load(files[0])
    n, iters, alpha = int(z0["ninner"]), int(z0["iters"]), float(z0["alpha"])
    checks = LD.parity_check(a.rate, n=n)
    ex = ThreadPoolExecutor(max_workers=a.threads)
    lams = [np.ascontiguousarray(np.load(f)["lam"], np.float32) for f in files]
    nblk = sum(len(l) for l in lams)
    print(f"{len(files)} batches, {nblk} blocks, rate {a.rate}, cap {iters}, alpha {alpha}, corr {a.corr}")

    def run(abort):
        t = time.perf_counter()
        out = [LD.layered_decode_batch(l, checks, n, iters=iters, alpha=alpha, ex=ex,
                                       corr=a.corr, abort=abort, pack=False) for l in lams]
        return out, time.perf_counter() - t

    truth, t_truth = run((0, 0, 0))
    conv_t = np.concatenate([r[1] for r in truth])
    its_t = np.concatenate([r[2] for r in truth])
    bad_t = np.concatenate([r[3] for r in truth])
    vs_t = sum(vec_sweeps(r[2]) for r in truth)
    print(f"truth (no abort): converged {conv_t.sum()} ({100*conv_t.mean():.1f}%), "
          f"lane-sweeps {its_t.sum()}, vector-sweeps {vs_t}, kernel {t_truth*1000:.0f} ms")
    ci = its_t[conv_t]
    print(f"  converged its p50/p90/p99/max = {np.percentile(ci,[50,90,99]).round(1).tolist()} / {ci.max() if len(ci) else '-'}; "
          f"failed final syndrome p5/p50/p95 = {np.percentile(bad_t[~conv_t],[5,50,95]).round(0).tolist() if (~conv_t).any() else '-'}")
    nc = len(checks) if hasattr(checks, "__len__") else None
    print("\n  rule (k, T x nc, stall)   lost  conv   lane-sweeps  vector-sweeps  kernel ms   (vs truth)")
    rows = []
    for k in (6, 8, 10, 12, 16):
        for T in (0.30, 0.33, 0.36, 0.40):
            for stall in (0, 4, 6, 8):
                out, dt = run((k, T, stall))
                conv = np.concatenate([r[1] for r in out])
                its = np.concatenate([r[2] for r in out])
                lost = int(np.count_nonzero(conv_t & ~conv))
                vs = sum(vec_sweeps(r[2]) for r in out)
                rows.append((k, T, stall, lost, int(conv.sum()), int(its.sum()), vs, dt))
                print(f"  ({k:2d}, {T:.2f}, {stall})            {lost:4d}  {int(conv.sum()):4d}  {int(its.sum()):10d}  {vs:12d}  {dt*1000:8.0f}   "
                      f"({100*(1-vs/vs_t):.0f}% fewer vector-sweeps)")
    zero = [r for r in rows if r[3] == 0]
    if zero:
        best = min(zero, key=lambda r: r[6])
        print(f"\n  best zero-loss rule here: k={best[0]} T={best[1]} stall={best[2]}: "
              f"{100*(1-best[6]/vs_t):.0f}% fewer vector-sweeps than no abort")
    cur = [r for r in rows if r[0] == 8 and r[1] == 0.36 and r[2] == 0][0]
    print(f"  shipping rule (8, 0.36, 0) here: lost {cur[3]}, {100*(1-cur[6]/vs_t):.0f}% fewer vector-sweeps")
    ex.shutdown(wait=True)


if __name__ == "__main__":
    sys.exit(main())
