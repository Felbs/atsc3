#!/usr/bin/env python3
"""gate_e136_q16.py -- the E136 int16 kernels against E131 as it shipped, block for block.

ldpc_layered_q16_ref is the E131 kernel, exported unchanged.  ldpc_layered_q16 (the same layout
with a vector load, no memset and a syndrome that stops early where only "zero or not" is asked)
and ldpc_layered_q16c (check-node messages stored as two magnitudes, an index and a sign mask)
claim the same arithmetic.  This gate holds them to it:

  1. banked real-air LLRs (RF25 Core, 64800-bit 9/15; a strong bank and the 12 dB copy), EVERY
     record that kept its LLRs, every rung, in the shipping vectors of sixteen, under several
     abort rules and with the correction on and off: bits, conv, its and bad of every block
     identical to the reference;
  2. the load alone (a call with 0 iterations leaves the quantised posteriors in the work
     buffer): identical int16 for real LLRs and for the values a careless caller could send --
     NaN, infinities, 1e30, exact halves, signed zero;
  3. every lane count 1..16, and a 16200-bit code (16200 is not a multiple of 16: the scalar
     tail of the load);
  4. time, one thread, the three kernels interleaved record by record.

Run: python lab/gate_e136_q16.py [--dump DIR ...] [--limit N] [--threads N] [--no-time]
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

NAMES = ("ref", "dense", "packed")


def kernels(lib):
    return {"ref": lib.ldpc_layered_q16_ref, "dense": lib.ldpc_layered_q16, "packed": lib.ldpc_layered_q16c}


def nwork(name, n, nc, nedges):
    return n * 16 + nc * 48 + nedges if name == "packed" else (n + nedges) * 16


def call(fn, name, llr, tabs, n, iters, alpha, scale, corr, abort, keep_work=False):
    """One kernel call on up to sixteen blocks -> (bits, conv, its, bad[, posteriors])."""
    cptr, cidx, _vptr, _vedge, nedges = tabs
    nc = len(cptr) - 1
    nblk = llr.shape[0]
    bits = np.empty((nblk, n), np.uint8)
    conv = np.zeros(nblk, np.uint8)
    its = np.full(nblk, -7, np.int32)
    bad = np.full(nblk, -7, np.int32)
    work = np.empty(nwork(name, n, nc, nedges), np.int16)
    ak, at, ast = abort
    at = int(round(at * nc)) if 0 < at < 1 else int(at)
    rc = fn(llr, nblk, n, nc, nedges, cptr, cidx, int(iters), np.float32(alpha), np.float32(scale),
            int(corr), int(ak), at, int(ast), bits, conv, its, bad, work)
    assert rc == 0, f"{name} returned {rc}"
    out = (bits, conv, its, bad)
    return out + (work[:n * 16].copy(),) if keep_work else out


def same(a, b):
    return all(np.array_equal(x, y) for x, y in zip(a, b))


def scale_of(llr):
    """The caller's scale, as layered_decode_batch computes it for the whole batch."""
    mean = float(np.mean(np.abs(llr)))
    return np.float32(LD.Q16_MEAN / mean) if mean > 0 else np.float32(1.0)


def records(d, limit):
    out = []
    for f in sorted(glob.glob(os.path.join(d, "*.npz"))):
        z = np.load(f)
        if z["lam"].size:
            out.append((np.ascontiguousarray(z["lam"], np.float32), float(z["alpha"]), int(z["iters"]),
                        int(z["ninner"])))
        if limit and len(out) >= limit:
            break
    return out


def gate_bank(lib, d, rate, limit, ex, settings):
    recs = records(d, limit)
    if not recs:
        print(f"\n{d}: no records with LLRs -- SKIPPED")
        return True, None
    n = recs[0][3]
    tabs = LD._kernel_tables(LD.parity_check(rate, n=n), n)
    K = kernels(lib)
    nblk = sum(len(r[0]) for r in recs)
    print(f"\n{d}: {len(recs)} records, {nblk} blocks, n {n}, rate {rate}")
    ok = True
    for corr, abort in settings:
        jobs = []
        for llr, alpha, iters, _n in recs:
            sc = scale_of(llr)
            for lo in range(0, len(llr), 16):
                jobs.append((llr[lo:lo + 16], alpha, iters, sc))

        def one(job):
            llr, alpha, iters, sc = job
            r = call(K["ref"], "ref", llr, tabs, n, iters, alpha, sc, corr, abort)
            res = {}
            for name in ("dense", "packed"):
                o = call(K[name], name, llr, tabs, n, iters, alpha, sc, corr, abort)
                res[name] = int(sum(not (np.array_equal(o[0][b], r[0][b]) and o[1][b] == r[1][b]
                                         and o[2][b] == r[2][b] and o[3][b] == r[3][b])
                                    for b in range(len(llr))))
            return len(llr), int(r[1].sum()), int(r[2].sum()), res

        got = list(ex.map(one, jobs))
        blocks = sum(g[0] for g in got)
        conv = sum(g[1] for g in got)
        sweeps = sum(g[2] for g in got)
        diff = {name: sum(g[3][name] for g in got) for name in ("dense", "packed")}
        flag = "" if not any(diff.values()) else "   <-- FAIL"
        ok = ok and not any(diff.values())
        print(f"  corr {corr} abort {abort!s:16s} {blocks:5d} blocks, {conv:5d} converged, its sum {sweeps:6d}: "
              f"dense differs on {diff['dense']}, packed on {diff['packed']}{flag}")
    return ok, (recs, tabs, n)


def gate_load(lib, llr_real, tabs, n):
    """The quantised posteriors, element for element, for real LLRs and for hostile ones."""
    K = kernels(lib)
    rng = np.random.default_rng(136)
    cases = [("real LLRs, 16 lanes", llr_real[:16].copy(), scale_of(llr_real))]
    for nb in (1, 2, 3, 7, 8, 9, 15):
        cases.append((f"real LLRs, {nb} lanes", llr_real[:nb].copy(), scale_of(llr_real)))
    odd = llr_real[:16].copy()
    flat = odd.reshape(-1)
    idx = rng.choice(flat.size, 4000, replace=False)
    specials = np.array([np.nan, np.inf, -np.inf, 1e30, -1e30, 0.0, -0.0, 3e9, -3e9, 2.1e9, -2.1e9, 24575.4, 24575.5,
                         24576.5, -24575.5, -24576.5, 32767.4, 32767.6, -32768.4, 65536.2, 1e-30, -1e-30],
                        np.float32)
    flat[idx] = specials[np.arange(len(idx)) % len(specials)]
    cases.append(("NaN, inf, 1e30, range edges, scale 1", odd, np.float32(1.0)))
    halves = (rng.integers(-30000, 30000, size=(16, n)).astype(np.float32) + np.float32(0.5))
    cases.append(("exact halves, scale 1", halves, np.float32(1.0)))
    cases.append(("exact halves, scale 0.5", halves, np.float32(0.5)))
    cases.append(("all zero", np.zeros((5, n), np.float32), np.float32(7.0)))
    ok = True
    print("\nload (posteriors after a 0-iteration call):")
    for label, llr, sc in cases:
        llr = np.ascontiguousarray(llr, np.float32)
        r = call(K["ref"], "ref", llr, tabs, n, 0, 0.82, sc, 1, (0, 0, 0), keep_work=True)
        line = []
        for name in ("dense", "packed"):
            o = call(K[name], name, llr, tabs, n, 0, 0.82, sc, 1, (0, 0, 0), keep_work=True)
            nd = int((o[4] != r[4]).sum())
            good = nd == 0 and same(o[:4], r[:4])
            ok = ok and good
            line.append(f"{name} {'identical' if good else 'DIFFERS on %d values' % nd}")
        print(f"  {label:38s} {', '.join(line)}{'' if 'DIFFERS' not in ' '.join(line) else '   <-- FAIL'}")
    return ok


def gate_lanes(lib, llr_real, tabs, n, alpha, iters):
    K = kernels(lib)
    ok = True
    sc = scale_of(llr_real)
    bad_at = []
    for nb in range(1, 17):
        llr = np.ascontiguousarray(llr_real[:nb])
        r = call(K["ref"], "ref", llr, tabs, n, iters, alpha, sc, 1, (8, 0.36, 0))
        for name in ("dense", "packed"):
            if not same(call(K[name], name, llr, tabs, n, iters, alpha, sc, 1, (8, 0.36, 0)), r):
                ok = False
                bad_at.append((nb, name))
    print(f"\nlane counts 1..16 on one record: {'identical' if ok else 'DIFFER at %s   <-- FAIL' % bad_at}")
    return ok


def gate_short(lib, rate, ex):
    """A 16200-bit code: 16200 = 16 x 1012 + 8, so the last eight values go through the scalar tail."""
    n = 16200
    try:
        tabs = LD._kernel_tables(LD.parity_check(rate, n=n), n)
    except Exception as e:                                       # noqa: BLE001
        print(f"\n16200-bit code at {rate}: not available ({e}) -- SKIPPED")
        return True
    K = kernels(lib)
    rng = np.random.default_rng(1360)
    ok = True
    nb_total = 0
    conv = 0
    for sigma in (0.55, 0.70, 0.80, 0.95):
        for nb in (16, 5):
            llr = np.ascontiguousarray((2.0 / sigma ** 2) * (1.0 + sigma * rng.standard_normal((nb, n))), np.float32)
            sc = scale_of(llr)
            for corr, abort in ((1, (8, 0.36, 0)), (0, (0, 0, 0))):
                r = call(K["ref"], "ref", llr, tabs, n, 50, 0.82, sc, corr, abort, keep_work=True)
                nb_total += nb
                conv += int(r[1].sum())
                for name in ("dense", "packed"):
                    o = call(K[name], name, llr, tabs, n, 50, 0.82, sc, corr, abort, keep_work=True)
                    if not same(o[:4], r[:4]):
                        ok = False
            r0 = call(K["ref"], "ref", llr, tabs, n, 0, 0.82, sc, 1, (0, 0, 0), keep_work=True)
            for name in ("dense", "packed"):
                o0 = call(K[name], name, llr, tabs, n, 0, 0.82, sc, 1, (0, 0, 0), keep_work=True)
                if not np.array_equal(o0[4], r0[4]):
                    ok = False
    print(f"\n16200-bit code at {rate} (synthetic all-zero word in noise, four noise levels): {nb_total} block "
          f"decodes, {conv} converged: {'identical' if ok else 'DIFFER   <-- FAIL'}")
    return ok


def timing(lib, d, recs, tabs, n, rounds):
    K = kernels(lib)
    first = min(r[1] for r in recs)
    recs = [r for r in recs if r[1] == first]
    tot = {name: 0.0 for name in NAMES}
    best = {name: 0.0 for name in NAMES}
    for llr, alpha, iters, _n in recs:
        sc = scale_of(llr)
        per = {name: [] for name in NAMES}
        for _ in range(rounds):
            for name in NAMES:
                t = time.perf_counter()
                for lo in range(0, len(llr), 16):
                    call(K[name], name, llr[lo:lo + 16], tabs, n, iters, alpha, sc, 1, (8, 0.36, 0))
                per[name].append(time.perf_counter() - t)
        for name in NAMES:
            tot[name] += sum(per[name]) / rounds
            best[name] += min(per[name])
    print(f"\ntime, one thread, {os.path.basename(d)}: {len(recs)} first-rung records (alpha {first}), "
          f"{sum(len(r[0]) for r in recs)} blocks, best of {rounds}")
    for name in NAMES:
        print(f"  {name:7s} {1e3 * best[name] / len(recs):7.1f} ms a record   "
              f"{100.0 * (best[name] / best['ref'] - 1.0):+6.1f} % vs ref   (mean {1e3 * tot[name] / len(recs):.1f})")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", nargs="*", default=[r"C:\captures\ldpc_dump_strong_e121", r"C:\captures\ldpc_dump_n12_e121"])
    ap.add_argument("--rate", default="9/15")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--rounds", type=int, default=3)
    ap.add_argument("--no-time", action="store_true")
    a = ap.parse_args(argv)
    lib = LD._layered_lib()
    if lib is None or not getattr(lib, "has_q16c", False):
        print("FAIL: no E136 kernel (python lab/build_ldpc_layered_kernel.py)")
        return 1
    ex = ThreadPoolExecutor(max_workers=a.threads)
    settings = [(1, (8, 0.36, 0)),      # shipping
                (0, (8, 0.36, 0)),      # correction off
                (1, (0, 0, 0)),         # no abort: the syndrome may stop early at every iteration but the last
                (1, (8, 0.36, 3)),      # stall rule on: every count is read
                (1, (3, 0, 0)),         # abort_k alone: nothing ever gives up
                (1, (2, 5000, 0))]      # an absolute threshold, early
    ok = True
    kept = []
    for d in a.dump:
        good, pack = gate_bank(lib, d, a.rate, a.limit, ex, settings)
        ok = ok and good
        if pack:
            kept.append((d, pack))
    if kept:
        recs, tabs, n = kept[0][1]
        ok = gate_load(lib, recs[0][0], tabs, n) and ok
        ok = gate_lanes(lib, recs[0][0], tabs, n, recs[0][1], recs[0][2]) and ok
    ok = gate_short(lib, a.rate, ex) and ok
    if not a.no_time:
        for d, (recs, tabs, n) in kept:
            timing(lib, d, recs, tabs, n, a.rounds)
    print("\nGATE", "PASS" if ok else "FAIL", "(dense and packed: bits, conv, its, bad of every block as the E131 reference)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
