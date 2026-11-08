#!/usr/bin/env python3
"""gate_e130_q8.py -- the sixteen-lane int8 LDPC kernel's referee, on banked real-air LLRs.

Against the float32 layered kernel (box-plus on, E112 abort on -- the shipping path), on the
first-rung records of the E121 dumps (RF25 Core, 64800-bit 9/15; a strong bank and the 12 dB
copy), for the int8 kernel ALONE at several quantiser scales and for the HYBRID (int8 first,
float32 on what int8 left):

  1. hybrid: converged count >= float32's at every record (the set can only grow), and the
     Baseband bits of every block both converge are identical;
  2. int8 alone: the loss at the cliff, printed per scale -- this is what the hybrid pays
     the float32 pass to recover;
  3. time: kernel milliseconds per record, each variant, same thread pool.

Run: python lab/gate_e130_q8.py [--dump DIR ...] [--means 8,12,20,32,48] [--limit N]
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


def load(d, limit):
    files = [f for f in sorted(glob.glob(os.path.join(d, "*.npz"))) if np.load(f)["lam"].size]
    alphas = {float(np.load(f)["alpha"]) for f in files}
    first = 0.82 if 0.82 in alphas else 1.0
    files = [f for f in files if float(np.load(f)["alpha"]) == first]
    if limit:
        files = files[:limit]
    lams = [np.ascontiguousarray(np.load(f)["lam"], np.float32) for f in files]
    z0 = np.load(files[0])
    return lams, int(z0["ninner"]), int(z0["iters"]), first


def run(lams, checks, n, iters, alpha, ex, q8, mean, q8_iters=None, q16=None):
    LD.Q8_DEFAULT = q8
    LD.Q8_MEAN = float(mean)
    if q8_iters is not None:
        LD.Q8_ITERS = int(q8_iters)
    LD.Q16_DEFAULT = "0"
    if q16 is not None:
        LD.Q16_DEFAULT = q16
        LD.Q16_MEAN = float(mean)
    t = time.perf_counter()
    out = [LD.layered_decode_batch(l, checks, n, iters=iters, alpha=alpha, ex=ex, corr=True) for l in lams]
    dt = time.perf_counter() - t
    assert all(o is not None for o in out), "kernel refused"
    return out, dt


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", nargs="*", default=[r"C:\captures\ldpc_dump_strong_e121", r"C:\captures\ldpc_dump_n12_e121"])
    ap.add_argument("--rate", default="9/15")
    ap.add_argument("--means", default="8,12,20,32,48")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--q8-iters", default="16,50")
    ap.add_argument("--q16", action="store_true", help="E131: evaluate the int16 kernel instead (means are its scale)")
    a = ap.parse_args(argv)
    lib = LD._layered_lib()
    if lib is None or not getattr(lib, "has_q8", False):
        print("FAIL: no int8 kernel (python lab/build_ldpc_layered_kernel.py)")
        return 1
    ex = ThreadPoolExecutor(max_workers=a.threads)
    ok = True
    for d in a.dump:
        lams, n, iters, alpha = load(d, a.limit)
        checks = LD.parity_check(a.rate, n=n)
        nblk = sum(len(l) for l in lams)
        print(f"\n{d}: {len(lams)} records, {nblk} blocks, rate {a.rate}, alpha {alpha}, "
              f"mean|LLR| {np.mean([np.abs(l).mean() for l in lams]):.2f}")
        ref, t_ref = run(lams, checks, n, iters, alpha, ex, "0", 20)
        conv_ref = sum(int(o[1].sum()) for o in ref)
        print(f"  {'variant':22s} {'conv':>6s} {'vs f32':>7s} {'bytes':>9s} {'ms/rec':>7s}")
        print(f"  {'float32 + box-plus':22s} {conv_ref:6d} {'':>7s} {'':>9s} {1e3 * t_ref / len(lams):7.1f}")
        for mean in [float(x) for x in a.means.split(",")]:
          for qi in [int(x) for x in a.q8_iters.split(",")]:
            for q8 in ("only", "1"):
                if a.q16:
                    out, dt = run(lams, checks, n, iters, alpha, ex, "0", mean, None,
                                  q16=("1" if q8 == "only" else "hybrid"))
                else:
                    out, dt = run(lams, checks, n, iters, alpha, ex, q8, mean, qi)
                conv = sum(int(o[1].sum()) for o in out)
                same = True
                lost = 0
                for o, r in zip(out, ref):
                    both = o[1] & r[1]
                    if both.any() and not np.array_equal(o[0][both], r[0][both]):
                        same = False
                    lost += int((r[1] & ~o[1]).sum())
                name = f"int{'16' if a.q16 else '8'} {'alone' if q8 == 'only' else 'hybrid'} m={mean:g} it={qi}"
                flag = ""
                if q8 == "1" and (lost or not same):
                    ok = False
                    flag = "  <-- FAIL"
                print(f"  {name:22s} {conv:6d} {conv - conv_ref:+7d} {'same' if same else 'DIFFER':>9s} "
                      f"{1e3 * dt / len(lams):7.1f}  lost {lost}{flag}")
    print("\nGATE", "PASS" if ok else "FAIL", "(hybrid: no block lost, bytes identical on shared blocks)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
