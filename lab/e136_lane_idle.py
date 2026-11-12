#!/usr/bin/env python3
"""e136_lane_idle.py -- how much of the int16 LDPC kernel is lanes with nothing to do?

E118 measured it for the eight-lane float32 kernel (26 % of vector-sweeps) and parked the
compaction at 0.98x.  Since then the default became the sixteen-lane int16 kernel (E131), the
ladder changed (E121) and the near-miss gate (E125) took most later-rung calls away.  This
re-measures the question on the recorded batches (ATSC3_LDPC_DUMP) under the kernel that ships,
and prices the schedules a compacting kernel could run -- from the per-block sweep counts alone,
which no lane layout changes (lanes never interact).

A vector-sweep is one pass of the kernel over all checks with all sixteen lanes, busy or not.

  now        : blocks cut into vectors of 16 in order, each run to its slowest lane (what ships)
  sorted     : the same, blocks first ordered by the weight of their syndrome before any sweep
  pack/2     : a refilling kernel, the batch cut into one share per thread (2), 16 lanes each
  pack/1     : a refilling kernel, the batch as one share
  stream     : a refilling kernel that never drains between batches (the bound for a pool that
               carries lanes across batches)
  floor      : the lanes own sweeps / 16

Run: python lab/e136_lane_idle.py --dump DIR [--rate 9/15] [--lanes 16]
"""
from __future__ import annotations

import argparse
import glob
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import m3_ldpc as LD                                                # noqa: E402


def vectors(its, lanes):
    """In-order vectors of `lanes`: each costs its slowest lane."""
    return sum(int(its[i:i + lanes].max()) for i in range(0, len(its), lanes))


def pack(its, lanes):
    """A refilling kernel on one share: sweeps until the last lane is done.  A block that needs
    0 sweeps retires at load and its lane takes the next one at once."""
    left = [int(x) for x in its if x > 0]
    lane = [0] * lanes
    sweeps = 0
    nxt = 0
    while True:
        for b in range(lanes):
            if lane[b] == 0 and nxt < len(left):
                lane[b] = left[nxt]
                nxt += 1
        if not any(lane):
            return sweeps
        step = min(x for x in lane if x)      # nothing changes until the first lane retires
        sweeps += step
        lane = [x - step if x else 0 for x in lane]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", required=True)
    ap.add_argument("--rate", default="9/15")
    ap.add_argument("--lanes", type=int, default=16)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    files = sorted(glob.glob(os.path.join(a.dump, "*.npz")))
    files = [f for f in files if np.load(f)["lam"].size]
    if a.limit:
        files = files[:a.limit]
    z0 = np.load(files[0])
    n = int(z0["ninner"])
    checks = LD.parity_check(a.rate, n=n)
    cptr, cidx = LD._kernel_tables(checks, n)[:2]
    L = a.lanes
    tot = dict(now=0, sorted=0, pack2=0, pack1=0)
    own = 0
    nblk = nconv = 0
    all_its = []
    t_kernel = 0.0
    hist = np.zeros(4, int)                  # 0, 1-8, 9-49, 50 sweeps
    lone = 0                                 # vector-sweeps run with 1 or 2 busy lanes (now)
    rows = []
    for f in files:
        z = np.load(f)
        lam = np.ascontiguousarray(z["lam"], np.float32)
        t = time.perf_counter()
        r = LD.layered_decode_batch(lam, checks, n, iters=int(z["iters"]), alpha=float(z["alpha"]), ex=None)
        t_kernel += time.perf_counter() - t
        its = np.asarray(r[2], np.int64)
        # the syndrome weight before any sweep: the hard decisions of the LLRs themselves
        hard = (lam < 0).astype(np.uint8)
        bc0 = np.bitwise_xor.reduceat(hard[:, cidx], cptr[:-1], axis=1).sum(axis=1, dtype=np.int64)
        order = np.argsort(bc0, kind="stable")
        half = (len(its) + 1) // 2
        v_now = vectors(its, L)
        v_sorted = vectors(its[order], L)
        v_p2 = pack(its[:half], L) + pack(its[half:], L)
        v_p1 = pack(its, L)
        tot["now"] += v_now; tot["sorted"] += v_sorted; tot["pack2"] += v_p2; tot["pack1"] += v_p1
        own += int(its.sum())
        nblk += len(its); nconv += int(np.count_nonzero(r[1]))
        all_its.extend(its.tolist())
        hist += np.array([np.count_nonzero(its == 0), np.count_nonzero((its > 0) & (its <= 8)),
                          np.count_nonzero((its > 8) & (its < int(z["iters"]))),
                          np.count_nonzero(its >= int(z["iters"]))])
        for i in range(0, len(its), L):
            v = np.sort(its[i:i + L])[::-1]
            lone += int(v[0] - (v[2] if len(v) > 2 else 0))
        cc = np.corrcoef(bc0, its)[0, 1] if its.std() > 0 and bc0.std() > 0 else float("nan")
        rows.append((os.path.basename(f)[-9:-4], len(its), int(np.count_nonzero(r[1])), v_now, v_sorted, v_p2, v_p1,
                     its.sum() / L, cc))
    stream = pack(np.array(all_its), L)
    print(f"E136 -- {len(files)} recorded calls under {a.dump}, rate {a.rate}, n {n}, {L} lanes, "
          f"the shipping kernel path (Q16 {LD.Q16_DEFAULT}, box-plus {int(LD.BOXPLUS_DEFAULT)})")
    print("  call  rows conv |  now sorted pack/2 pack/1 | floor | corr(syndrome0, sweeps)")
    for r in rows:
        print(f"  {r[0]} {r[1]:5d} {r[2]:4d} | {r[3]:4d} {r[4]:6d} {r[5]:6d} {r[6]:6d} | {r[7]:5.1f} | {r[8]:+.2f}")
    fl = own / L
    print(f"\n  blocks {nblk}, converged {nconv}; sweeps per block: 0: {hist[0]}, 1-8: {hist[1]}, "
          f"9-49: {hist[2]}, the cap: {hist[3]}")
    print(f"  kernel time, one thread: {1000 * t_kernel:.0f} ms = {1000 * t_kernel / max(tot['now'], 1):.2f} ms a vector-sweep")
    print(f"  vector-sweeps   now {tot['now']}   (of them {lone} with one or two lanes busy)")
    for k, name in (("sorted", "sorted"), ("pack2", "pack/2"), ("pack1", "pack/1")):
        print(f"                  {name:7s} {tot[k]:6d}  {100 * (tot[k] / tot['now'] - 1):+.1f} %")
    print(f"                  stream  {stream:6d}  {100 * (stream / tot['now'] - 1):+.1f} %")
    print(f"                  floor   {fl:6.0f}  {100 * (fl / tot['now'] - 1):+.1f} %")
    return 0


if __name__ == "__main__":
    sys.exit(main())
