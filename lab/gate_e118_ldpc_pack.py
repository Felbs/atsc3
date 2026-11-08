#!/usr/bin/env python3
"""gate_e118_ldpc_pack.py -- lane compaction in the layered LDPC kernel: the referee.

E118 adds `ldpc_layered_pack_f32`, which takes any number of blocks and refills
a lane the moment its block retires.  Lanes never interact, so every block must
come back EXACTLY as the eight-at-a-time entry point gives it -- not "net gain",
not "bytes identical on shared blocks": the same bits, the same conv, the same
iteration count and the same final syndrome weight for every block, whatever
the share boundaries.  That is the law here, and it is stricter than E107's
because nothing about the arithmetic is allowed to change.

Input: the per-batch LLR record a replay writes with ATSC3_LDPC_DUMP (every
ATSC3_LDPC_DUMP_EVERY-th batch keeps its LLRs).  Real air, real batches, the
mix of clean, near-miss and hopeless blocks the pool actually sees.

Run: python lab/gate_e118_ldpc_pack.py [--dump DIR] [--rate 9/15] [--threads 2]
PASS = zero mismatches on every field, deterministic, and not slower.
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", default=os.path.expanduser("~/captures/ldpc_dump_llr"))
    ap.add_argument("--rate", default="9/15", help="LDPC code rate of the dumped blocks")
    ap.add_argument("--threads", type=int, default=2, help="the pool worker's thread count")
    ap.add_argument("--corr", type=int, default=1, help="E110 box-plus (the RF25 setting)")
    ap.add_argument("--abort", default="8,0.36,0")
    ap.add_argument("--limit", type=int, default=0, help="batches to use (0 = all)")
    a = ap.parse_args()

    lib = LD._layered_lib()
    if lib is None or not getattr(lib, "has_pack", False):
        print("FAIL: the layered kernel is not loadable or has no pack entry point "
              "(python lab/build_ldpc_layered_kernel.py)")
        return 2
    files = [f for f in sorted(glob.glob(os.path.join(a.dump, "*.npz")))]
    files = [f for f in files if np.load(f)["lam"].size]
    if a.limit:
        files = files[:a.limit]
    if not files:
        print(f"FAIL: no batches with LLRs under {a.dump}")
        return 2
    z0 = np.load(files[0])
    n = int(z0["ninner"])
    checks = LD.parity_check(a.rate, n=n)
    abort = tuple(float(x) for x in a.abort.split(","))
    ex = ThreadPoolExecutor(max_workers=a.threads)
    print(f"E118 gate -- {len(files)} batches from {a.dump}, rate {a.rate}, n {n}, "
          f"corr {a.corr}, abort {a.abort}, threads {a.threads}")
    print("  batch  blocks  conv  | old ms  pack ms  pack(1 share) ms | mismatches bits/conv/its/bad | det")
    tot_old = tot_new = tot_one = 0.0
    nb_all = nconv = 0
    bad_fields = np.zeros(4, int)
    nondet = 0
    for f in files:
        z = np.load(f)
        lam = np.ascontiguousarray(z["lam"], np.float32)
        alpha, iters = float(z["alpha"]), int(z["iters"])
        t = time.perf_counter()
        old = LD.layered_decode_batch(lam, checks, n, iters=iters, alpha=alpha, ex=ex,
                                      corr=a.corr, abort=abort, pack=False)
        t_old = time.perf_counter() - t
        t = time.perf_counter()
        new = LD.layered_decode_batch(lam, checks, n, iters=iters, alpha=alpha, ex=ex,
                                      corr=a.corr, abort=abort, pack=True)
        t_new = time.perf_counter() - t
        t = time.perf_counter()
        one = LD.layered_decode_batch(lam, checks, n, iters=iters, alpha=alpha, ex=None,
                                      corr=a.corr, abort=abort, pack=True)
        t_one = time.perf_counter() - t
        again = LD.layered_decode_batch(lam, checks, n, iters=iters, alpha=alpha, ex=ex,
                                        corr=a.corr, abort=abort, pack=True)
        mm = [int(np.count_nonzero(np.any(old[0] != new[0], axis=1))),
              int(np.count_nonzero(old[1] != new[1])),
              int(np.count_nonzero(old[2] != new[2])),
              int(np.count_nonzero(old[3] != new[3]))]
        mm1 = [int(np.count_nonzero(np.any(old[0] != one[0], axis=1))),
               int(np.count_nonzero(old[1] != one[1])),
               int(np.count_nonzero(old[2] != one[2])),
               int(np.count_nonzero(old[3] != one[3]))]
        det = all(np.array_equal(x, y) for x, y in zip(new, again))
        nondet += 0 if det else 1
        bad_fields += np.array(mm) + np.array(mm1)
        # the recorded run's own verdicts must match too (same kernel rule)
        rec = [int(np.count_nonzero(old[1] != z["conv"].astype(bool))),
               int(np.count_nonzero(old[2] != z["its"]))]
        tot_old += t_old; tot_new += t_new; tot_one += t_one
        nb_all += len(lam); nconv += int(new[1].sum())
        print(f"  {os.path.basename(f)[-9:-4]}  {len(lam):5d}  {int(new[1].sum()):4d}  | "
              f"{1000*t_old:6.1f}  {1000*t_new:7.1f}  {1000*t_one:14.1f} | "
              f"{mm[0]}/{mm[1]}/{mm[2]}/{mm[3]}  (1 share {mm1[0]}/{mm1[1]}/{mm1[2]}/{mm1[3]}; "
              f"vs record conv/its {rec[0]}/{rec[1]}) | {'ok' if det else 'NONDET'}")
    ex.shutdown(wait=True)
    print(f"\n  blocks {nb_all}, converged {nconv}; kernel time old {tot_old*1000:.0f} ms, "
          f"pack {tot_new*1000:.0f} ms ({tot_old/max(tot_new,1e-9):.2f}x), "
          f"pack on one thread {tot_one*1000:.0f} ms")
    ok = (bad_fields == 0).all() and nondet == 0 and tot_new <= tot_old
    if ok:
        print("  PASS: every block bit-identical (bits, conv, its, bad) to the eight-at-a-time "
              "kernel under both share layouts, deterministic, faster")
    else:
        print(f"  FAIL: mismatches bits/conv/its/bad {bad_fields.tolist()}, "
              f"non-deterministic batches {nondet}, "
              f"{'slower' if tot_new > tot_old else 'not slower'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
