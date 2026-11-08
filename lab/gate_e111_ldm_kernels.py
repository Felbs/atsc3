#!/usr/bin/env python3
"""gate_e111_ldm_kernels.py -- the compiled equaliser and nearest-point kernels
on the LDM path: the referee.

E54's boundary, asserted at the tightest point there is.  The equaliser's
whole output is the cell POOL; the CPE's whole output is the corrected pool.
So, on real LDM Frames (RF25, 189 symbols, SP12_4, the E60 smoothed estimate
live), with the kernels ON against OFF:

  1. cell_pool: the owner map and the coherence report are BIT-identical
     (compared as raw bytes).  The pool is held to a law split by what numpy
     does: BIT-identical for every cell no complex multiply by the E60 gain
     touches (the preamble spare cells; every un-smoothed row -- measured
     0 diffs on real air), and within ULP_BOUND of full scale for the
     smoothed rows.  Why not bit-identical there: numpy's complex128
     multiply reproduces neither the textbook two-rounding form nor an
     FMA-contracted one (measured here -- 2263/6625 cells of one row's
     `H *= cl` differ from the plain-double form -- and before, in
     cellpool_kernel.c's derot_f64 note).  The kernel IS the plain-double
     form, exactly (diag: 0 diffs against a scalar reference).  The bound is
     a few hundred ULP of full scale -- thirteen orders below any real
     error, which is O(1) -- so a wrong index or a wrong gain still fails.
  2. cpe_fast on the SAME input pool: the output is bit-identical.
  3. core_cells end to end (both kernels on vs both off): bit-identical.
  4. Determinism: two kernel-on runs are bit-identical.
  5. Not slower: equalise and cpe each take less time with the kernel.

Run: python lab/gate_e111_ldm_kernels.py [--capture ~/captures/rf25_marginal_45s.cs16]
                                         [--frames 6] [--threads 8]
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m44_ldm as M44                                            # noqa: E402


#: the smoothed-row pool bound, as a fraction of the pool's RMS.  Measured on
#: RF25: max |diff| ~1e-15 on an RMS-order-1 pool, i.e. a few ULP.  1e-14 of
#: full scale is ~45 ULP of a unit cell -- generous for rounding, thirteen
#: orders too tight for a real error.
ULP_BOUND = 1e-14
EPS = 2.220446049250313e-16


def _bits(a):
    return np.ascontiguousarray(a).view(np.uint8)


def same(a, b):
    a, b = _bits(a), _bits(b)
    return a.shape == b.shape and bool(np.array_equal(a, b))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--capture", default=os.path.expanduser("~/captures/rf25_marginal_45s.cs16"))
    ap.add_argument("--rate", type=float, default=6.912e6)
    ap.add_argument("--frames", type=int, default=6)
    ap.add_argument("--start", type=float, default=2.0)
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args(argv)
    if M44._ldm_cellpool_kern() is None or M44._ldm_cpe_kern() is None:
        print("FAIL: a kernel is not built (lab/build_cellpool_kernel.py, lab/build_demap_kernel.py)")
        return 1

    from m2_pilots import load
    from m10_core import fine_t0
    raw = (np.fromfile(a.capture, np.int16, count=int(a.rate)).astype(np.float32)
           .view(np.complex64) / 32768.0)
    plan, info = M44.sniff(raw, a.rate, log=lambda *x: None)
    if plan is None:
        print("FAIL: sniff:", info.get("why")); return 1
    g = plan.g
    ex = ThreadPoolExecutor(a.threads)
    fd = M44.FrameDemod(plan, threads=a.threads, ex=ex)
    fd.prewarm()
    # see which rows took the smoothed estimate: the pool law is split on it
    capt = {}
    _orig_ce = fd._ce_smooth
    def _ce_wrap(Yb, cp):
        r = _orig_ce(Yb, cp); capt["ce"] = r; return r
    fd._ce_smooth = _ce_wrap
    per = plan.frame_samples / a.rate
    print("E111 gate -- compiled equaliser + nearest-point kernels on the LDM path, real air")
    print("=" * 92)
    print(f"  {os.path.basename(a.capture)}: {g.nsym} symbols/frame, {a.frames} frames, {a.threads} threads")
    print(f"\n  {'frame':>5} {'pool':>6} {'ulp':>5} {'sm':>4} {'owner':>6} {'coh':>5} {'cpe':>5} {'e2e':>5} {'det':>5} "
          f"{'eq off':>7} {'eq on':>6} {'cpe off':>8} {'cpe on':>7} {'core off':>9} {'core on':>8}")
    fails = []
    tot = dict(eq0=0.0, eq1=0.0, cpe0=0.0, cpe1=0.0, core0=0.0, core1=0.0)
    for k in range(a.frames):
        y, _fs, _cfo, _h = load(a.capture, a.rate, fmt="cs16", span_sec=per + 0.06,
                                start_sec=a.start + k * per)
        t0, _coh = fine_t0(y, plan.g)       # the Geometry, as the tracker's t0 is

        def run(cell_k, cpe_k):
            fd.use_cellpool_kernel = cell_k
            fd.use_cpe_kernel = cpe_k
            fd.tm.clear()
            rep = {}
            pool, owner = fd.cell_pool(y, t0, rep)
            t_eq = fd.tm["equalise"]
            t = time.perf_counter()
            out = fd.cpe_fast(pool.copy(), owner)
            t_cpe = time.perf_counter() - t
            return pool, owner, rep, out, t_eq, t_cpe

        p0, o0, r0, _c0, eq0, _ = run(False, False)
        ce0 = capt.get("ce")
        p1, o1, r1, _c1, eq1, _ = run(True, False)
        # the CPE leg isolates the search: both run on the OFF pool
        fd.use_cpe_kernel = False
        t = time.perf_counter(); c0 = fd.cpe_fast(p0.copy(), o0); cpe0 = time.perf_counter() - t
        fd.use_cpe_kernel = True
        t = time.perf_counter(); c1 = fd.cpe_fast(p0.copy(), o0); cpe1 = time.perf_counter() - t
        # end to end
        fd.use_cellpool_kernel = False; fd.use_cpe_kernel = False
        t = time.perf_counter(); e0 = fd.core_cells(y, t0, {}); core0 = time.perf_counter() - t
        fd.use_cellpool_kernel = True; fd.use_cpe_kernel = True
        t = time.perf_counter(); e1 = fd.core_cells(y, t0, {}); core1 = time.perf_counter() - t
        e2 = fd.core_cells(y, t0, {})

        # the pool law, split: exact where no complex multiply by the E60 gain
        # is involved; within ULP_BOUND of full scale where one is
        ok_vec = ce0["ok"] if ce0 is not None else np.zeros(g.nsym, bool)
        is_sm = (o0 >= 0) & ok_vec[np.maximum(o0, 0)]
        rms = float(np.sqrt(np.mean(np.abs(p0) ** 2)))
        d = np.abs(p0 - p1)
        ok_exact = same(p0[~is_sm], p1[~is_sm])
        dmax_sm = float(d[is_sm].max()) if is_sm.any() else 0.0
        ok_ulp = dmax_sm <= ULP_BOUND * rms
        ok_pool = ok_exact and ok_ulp
        ok_own = same(o0, o1)
        n_sm = int(ok_vec.sum())
        ok_coh = (r0.get("data_coherence_mean") == r1.get("data_coherence_mean")
                  and r0.get("ce_smoothed_syms") == r1.get("ce_smoothed_syms"))
        ok_cpe = same(c0, c1); ok_e2e = same(e0, e1); ok_det = same(e1, e2)
        for name, ok in (("pool", ok_pool), ("owner", ok_own), ("coherence", ok_coh),
                         ("cpe", ok_cpe), ("core_cells", ok_e2e), ("determinism", ok_det)):
            if not ok:
                fails.append(f"frame {k}: {name} differs with the kernel on")
        if not ok_pool:
            fails.append(f"    pool: exact part identical={ok_exact}; smoothed part max|diff| {dmax_sm:.3e} "
                         f"vs bound {ULP_BOUND * rms:.3e} (rms {rms:.3e}); {int(np.count_nonzero(d))} of {len(p0)} cells differ")
        tot["eq0"] += eq0; tot["eq1"] += eq1; tot["cpe0"] += cpe0; tot["cpe1"] += cpe1
        tot["core0"] += core0; tot["core1"] += core1
        Y = lambda b: "  ok " if b else " DIFF"
        print(f"  {k:>5} {Y(ok_pool)} {dmax_sm / max(rms, 1e-300) / EPS:>5.1f} {n_sm:>4} {Y(ok_own)} {Y(ok_coh)} {Y(ok_cpe)} {Y(ok_e2e)} {Y(ok_det)} "
              f"{eq0*1e3:>7.1f} {eq1*1e3:>6.1f} {cpe0*1e3:>8.1f} {cpe1*1e3:>7.1f} "
              f"{core0*1e3:>9.1f} {core1*1e3:>8.1f}")
    n = a.frames
    print(f"\n  per frame: equalise {tot['eq0']/n*1e3:.1f} -> {tot['eq1']/n*1e3:.1f} ms "
          f"({tot['eq0']/max(tot['eq1'],1e-9):.2f}x)   cpe {tot['cpe0']/n*1e3:.1f} -> {tot['cpe1']/n*1e3:.1f} ms "
          f"({tot['cpe0']/max(tot['cpe1'],1e-9):.2f}x)   core_cells {tot['core0']/n*1e3:.0f} -> {tot['core1']/n*1e3:.0f} ms "
          f"({tot['core0']/max(tot['core1'],1e-9):.2f}x)")
    if tot["eq1"] > tot["eq0"]:
        fails.append("equaliser kernel slower than numpy over the gate")
    if tot["cpe1"] > tot["cpe0"]:
        fails.append("cpe kernel slower than numpy over the gate")
    fd.use_cellpool_kernel = M44._LDM_CELLPOOL_KERNEL
    fd.use_cpe_kernel = M44._LDM_CPE_KERNEL
    if fails:
        print("\n  FAIL:\n    " + "\n    ".join(fails))
        return 1
    print("\n  PASS: pool exact where no gain-multiply is involved and within a few ULP where one is; "
          "owner, coherence, cpe and core_cells bit-identical on every frame; deterministic; faster")
    return 0


if __name__ == "__main__":
    sys.exit(main())
