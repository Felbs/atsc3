#!/usr/bin/env python3
"""gate_e121_ladder.py -- the LDPC alpha ladder, judged on banked LLR batches.

What the per-block record said (2026-11-07, RF25 banked on the 32-core box,
CORE = QPSK 9/15 64800, E110 correction on): the ladder is (1.0, 0.85, 0.75)
and every failing block climbs it.  On the as-recorded air the first rung
converged 2,449 blocks and the second rung converged 2,556 MORE -- blocks the
first rung had ABORTED AS HOPELESS at sweep 8 -- and the third rung converged
1 for 30% of all sweeps.  At 12 dB of added noise the first rung converged
NOTHING (0 of 4,716) and the second rung 954.  M10 chose 1.0 first for the
64800 rate-6/15 code; this carrier's core is 9/15, and with the correction
on the scale is not the same.

So: each variant below is run on the SAME recorded LLR batches with the
ladder's own semantics (a rung only sees the rows still failing; a rung is
only "taken" where it converges or improves the syndrome), and judged on the
converged set against the shipping ladder: blocks lost, blocks gained,
lane-sweeps, vector-sweeps and kernel time.  A variant may replace the
default only if it loses none and is faster.

Run: python lab/gate_e121_ladder.py --dump DIR [--dump DIR2] [--rate 9/15]
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

# (ladder, near): `near` = E125's gate, a failing Block climbs only if its final
# syndrome is under near x nc; 0 = every failing Block climbs.  The first entry
# is the reference every other row's "lost" is counted against.
VARIANTS = {
    "E121 ship (0.82, 0.85, 1.0)": ((0.82, 0.85, 1.0), 0.0),
    "M10 (1.0, 0.85, 0.75)": ((1.0, 0.85, 0.75), 0.0),
    "E121 + near 0.36": ((0.82, 0.85, 1.0), 0.36),
    "E121 + near 0.10": ((0.82, 0.85, 1.0), 0.10),
    "E121 + near 0.05": ((0.82, 0.85, 1.0), 0.05),
    "E121 + near 0.02": ((0.82, 0.85, 1.0), 0.02),
    "E121 + near 0.01": ((0.82, 0.85, 1.0), 0.01),
    "0.82 only": ((0.82,), 0.0),
}


def vec_sweeps(its):
    return sum(int(its[lo:lo + 8].max()) for lo in range(0, len(its), 8))


def run_ladder(lam, checks, n, iters, corr, abort, ladder, ex, near=0.0):
    """m44_ldm.ldpc_batch's semantics on one batch: -> conv, sweeps, vector-sweeps."""
    nb = len(lam)
    conv = np.zeros(nb, bool)
    bits = np.zeros((nb, n), np.uint8)
    bad = np.full(nb, 1 << 30, np.int32)
    todo = np.arange(nb)
    sw = vs = 0
    for a in ladder:
        if not len(todo):
            break
        sub = lam if len(todo) == nb else lam[todo]
        b, c, its, d = LD.layered_decode_batch(sub, checks, n, iters=iters, alpha=a, ex=ex,
                                               corr=corr, abort=abort, pack=False)
        sw += int(its.sum()); vs += vec_sweeps(its)
        take = c | (d < bad[todo])
        sel = todo[take]
        conv[sel] = c[take]
        bits[sel] = b[take]
        bad[sel] = d[take]
        todo = todo[~c & (d < near * len(checks))] if near > 0 else todo[~c]
    return conv, sw, vs, bits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", action="append", required=True)
    ap.add_argument("--rate", default="9/15")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--corr", type=int, default=1)
    ap.add_argument("--abort", default="8,0.36,0")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    abort = tuple(float(x) for x in a.abort.split(","))
    ex = ThreadPoolExecutor(max_workers=a.threads)
    verdict_ok = True
    for d in a.dump:
        files = [f for f in sorted(glob.glob(os.path.join(d, "*.npz"))) if np.load(f)["lam"].size]
        # only FIRST-rung records carry a batch's full LLRs: alpha 1.0 under M10's
        # ladder (the 13.20 dumps), 0.82 under E121's (the 13.22 dumps)
        alphas = {float(np.load(f)["alpha"]) for f in files}
        first = 0.82 if 0.82 in alphas else 1.0
        files = [f for f in files if float(np.load(f)["alpha"]) == first]
        if a.limit:
            files = files[:a.limit]
        lams = [np.ascontiguousarray(np.load(f)["lam"], np.float32) for f in files]
        z0 = np.load(files[0])
        n, iters = int(z0["ninner"]), int(z0["iters"])
        checks = LD.parity_check(a.rate, n=n)
        nblk = sum(len(l) for l in lams)
        print(f"\n{d}: {len(files)} first-rung batches, {nblk} blocks, rate {a.rate}, corr {a.corr}, abort {a.abort}")
        print(f"  {'variant':26s} {'conv':>5s} {'lost':>5s} {'gain':>5s} {'lane-sw':>8s} {'vec-sw':>7s} {'kernel ms':>9s}")
        ref = None
        for name, (ladder, near) in VARIANTS.items():
            t = time.perf_counter()
            res = [run_ladder(l, checks, n, iters, a.corr, abort, ladder, ex, near=near) for l in lams]
            dt = time.perf_counter() - t
            conv = np.concatenate([r[0] for r in res])
            bits = np.concatenate([r[3] for r in res])
            sw = sum(r[1] for r in res); vs = sum(r[2] for r in res)
            if ref is None:
                ref = (conv, sw, vs, dt, bits)
            lost = int(np.count_nonzero(ref[0] & ~conv)); gain = int(np.count_nonzero(conv & ~ref[0]))
            shared = ref[0] & conv
            bytes_diff = int(np.count_nonzero(np.any(bits[shared] != ref[4][shared], axis=1)))
            print(f"  {name:26s} {int(conv.sum()):5d} {lost:5d} {gain:5d} {sw:8d} {vs:7d} {dt*1000:9.0f}"
                  f"   ({100*(1-vs/max(ref[2],1)):.0f}% fewer vec-sweeps than ship; shared blocks with different bits: {bytes_diff})")
    ex.shutdown(wait=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
