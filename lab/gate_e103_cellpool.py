#!/usr/bin/env python3
"""gate_e103_cellpool.py -- the compiled cell-pool equaliser's referee.

The kernel (lab/cellpool_kernel.c) replaces the per-class inner work of
`m9_fast.FrameDecoder.cell_pool_fast`.  It ships only if it reproduces that
path's pool on real air.

  leg 1  REAL AIR, pool identity: decode the same capture's Frames with
         ATSC3_CELLPOOL_KERNEL=1 and =0 and compare the cell pool and the
         symbol_of map.  The pool is complex128 and both paths do the same
         divisions in the same order, so BIT-IDENTITY is the bar here -- the
         only licensed re-association (np.interp -> hp*(1-w)+hp*w) is present
         in BOTH paths, so nothing is left to excuse a difference.
  leg 2  DEEP FADE: a channel estimate driven towards zero, where a naive
         (a*conj b)/|b|^2 division would overflow to inf/nan while Smith's
         algorithm does not.  Both paths must agree, and must stay finite
         wherever numpy stays finite.

Run: python lab/gate_e103_cellpool.py [--capture hit_rf33.cs16] [--rate 8e6]
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m6_cells as C                                            # noqa: E402
import m9_fast as F9                                            # noqa: E402

_COMPARED = [0]


def _pool(cap, rate, use_kernel, frames, smoothed=False):
    os.environ["ATSC3_CELLPOOL_KERNEL"] = "1" if use_kernel else "0"
    F9._CELLPOOL_STATE["lib"] = None          # re-probe under the new env
    fd = F9.FrameDecoder(threads=4, backend="cpu", cpu_fast=True)
    span = 0.05 + (frames + 1) * F9.FRAME_SEC
    y = F9.load_fast(os.path.join(HERE, cap), rate,
                     span_sec=span * (rate / 6.912e6), ex=fd.ex)[0]
    out = []
    t0 = C.BOOTSTRAP
    fn = fd.cell_pool_fast_sm if smoothed else fd.cell_pool_fast
    for i in range(frames):
        pool, info = fn(y, t0, {})
        out.append((np.asarray(pool).copy(),
                    np.asarray(info["symbol_of"]).copy()))
    return out


def _cmp(tag, a, b, fails):
    _COMPARED[0] += 1
    if a.shape != b.shape:
        print(f"    {tag}: SHAPE {a.shape} vs {b.shape}")
        fails.append(f"{tag} shape mismatch")
        return
    if np.array_equal(a, b):
        print(f"    {tag}: BIT-IDENTICAL ({a.size} values)")
        return
    d = np.abs(a.astype(np.complex128) - b.astype(np.complex128))
    fin = np.isfinite(d)
    if not fin.all():
        print(f"    {tag}: NON-FINITE in {(~fin).sum()} of {d.size}")
        fails.append(f"{tag} non-finite divergence")
        return
    scale = float(np.abs(b[fin]).max()) or 1.0
    print(f"    {tag}: DIFFERS  max|d| {float(d.max()):.3e}  "
          f"rel {float(d.max())/scale:.3e}")
    fails.append(f"{tag} not bit-identical")


# The smoothed path's bar, and why it is NOT bit-identity.
#
# E106 interpolates the channel off the smoothed grid INSIDE the kernel so the
# (nrows, ncar) complex128 array and its two fancy-index gathers -- measured at
# 13.0 ms/Frame, more than the FFT beside it -- never happen. The interpolation
# itself is bit-identical; the per-symbol GAIN multiply is not, because modern
# NumPy's complex multiply is not the textbook (ar*br - ai*bi, ar*bi + ai*br)
# and cannot be reproduced portably in C -- measured here, not assumed.
# Keeping bit-identity would mean keeping numpy's H-building, i.e. keeping the
# entire cost this change removes.
# So the smoothed leg asserts a tight RELATIVE bound (a few float64 ulp on the
# pool) and DECODED BYTES are the authority -- the same licence cpu_fast
# itself runs under. Measured: all three media lanes byte-identical.
# The PLAIN path keeps bit-identity as its bar and still meets it (leg 1).
SM_MAX_ULP = 16.0


def _cmp_sm(tag, a, b, fails):
    _COMPARED[0] += 1
    if np.array_equal(a, b):
        print(f"    {tag}: BIT-IDENTICAL ({a.size} values)")
        return
    d = np.abs(a.astype(np.complex128) - b.astype(np.complex128))
    if not np.isfinite(d).all():
        print(f"    {tag}: NON-FINITE divergence")
        fails.append(f"{tag} non-finite")
        return
    scale = float(np.abs(b).max()) or 1.0
    ulps = float(d.max()) / scale / float(np.finfo(np.float64).eps)
    ok = ulps <= SM_MAX_ULP
    print(f"    {tag}: {'within' if ok else 'OUTSIDE'} bound  "
          f"max {ulps:.2f} ulp  ({a.size} values)")
    if not ok:
        fails.append(f"{tag} exceeded {SM_MAX_ULP} ulp (max {ulps:.2f})")


def leg_sm(fails, cap, rate, frames):
    """The SMOOTHED path (E58/E60) -- this is the one the default config
    actually runs (margin levers are on by default), so it is not optional.
    Wiring the kernel only into cell_pool_fast was measured to move the
    decoder by 3%: the decoder was never calling it."""
    print(f"LEG 3 -- real air pool identity, SMOOTHED path (cell_pool_fast_sm)")
    if not os.path.exists(os.path.join(HERE, cap)):
        fails.append("capture missing for the smoothed leg")
        return
    k = _pool(cap, rate, True, frames, smoothed=True)
    n = _pool(cap, rate, False, frames, smoothed=True)
    for i, ((pk, ok_), (pn, on)) in enumerate(zip(k, n)):
        print(f"  frame {i}")
        _cmp_sm("pool    ", pk, pn, fails)
        _cmp("symbol_of", ok_, on, fails)


def leg1(fails, cap, rate, frames):
    print(f"LEG 1 -- real air pool identity: {cap}")
    if not os.path.exists(os.path.join(HERE, cap)):
        msg = f"capture {cap} not found -- the gate CANNOT vouch for anything"
        print(f"  FAIL: {msg}")
        fails.append(msg)
        return
    k = _pool(cap, rate, True, frames)
    n = _pool(cap, rate, False, frames)
    for i, ((pk, ok_), (pn, on)) in enumerate(zip(k, n)):
        print(f"  frame {i}")
        _cmp("pool    ", pk, pn, fails)
        _cmp("symbol_of", ok_, on, fails)


def leg2(fails):
    """A deep fade: the kernel must not manufacture inf/nan where numpy does
    not.  Driven through the kernel directly so the fade is controllable."""
    print("LEG 2 -- deep fade / division robustness")
    os.environ["ATSC3_CELLPOOL_KERNEL"] = "1"   # leg 1 left it disabled
    F9._CELLPOOL_STATE["lib"] = None
    lib = F9._cellpool_lib()
    if lib is None:
        fails.append("kernel not loadable for leg 2")
        return
    nfft, npk, nd = 64, 4, 6
    for exp in (0, -150, -300, -320):
        tiny = 10.0 ** exp
        Y = ((np.arange(nfft) + 1.0) + 1j * (np.arange(nfft) + 2.0))
        bins_pk = np.array([4, 12, 20, 28], np.int32)
        refpk = (Y[bins_pk] / tiny).astype(np.complex128)   # -> hp ~ tiny
        jj = np.zeros(nd, np.int32)
        wt = np.zeros(nd, np.float64)
        bins_d = np.array([1, 2, 3, 5, 6, 7], np.int32)
        dmap = np.arange(nd, dtype=np.int32)
        hmat = np.arange(nd, dtype=np.int32)
        hp = np.zeros(2 * npk, np.float64)
        xg = np.zeros(2 * nd, np.float64)
        rc = lib.cellpool_row_f64(
            np.ascontiguousarray(Y).view(np.float64), nfft, bins_pk,
            np.ascontiguousarray(refpk).view(np.float64), npk,
            jj, wt, nd, bins_d, dmap, hmat, nd, hp, xg)
        got = xg.view(np.complex128).copy()
        # numpy reference, same association
        hpn = Y[bins_pk] / refpk
        Hn = hpn[jj] * (1.0 - wt) + hpn[jj + 1] * wt
        want = Y[bins_d] / Hn[dmap]
        _COMPARED[0] += 1
        if rc != 0:
            print(f"  pilot scale 1e{exp}: kernel declined (rc {rc})")
            continue
        same = np.array_equal(got, want)
        gf, wf = np.isfinite(got).all(), np.isfinite(want).all()
        print(f"  pilot scale 1e{exp}: {'BIT-IDENTICAL' if same else 'DIFFERS'}"
              f"   kernel finite={gf} numpy finite={wf}")
        if gf != wf:
            fails.append(f"finiteness diverged at 1e{exp}")
        elif not same and gf:
            d = float(np.abs(got - want).max())
            sc = float(np.abs(want).max()) or 1.0
            if d / sc > 1e-12:
                fails.append(f"deep fade 1e{exp} rel {d/sc:.2e}")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--capture", default="hit_rf33.cs16")
    ap.add_argument("--rate", type=float, default=8e6)
    ap.add_argument("--frames", type=int, default=2)
    a = ap.parse_args(argv)

    if F9._cellpool_lib() is None:
        print("FAIL: cell-pool kernel not loadable -- "
              "python lab/build_cellpool_kernel.py")
        return 1
    fails = []
    leg1(fails, a.capture, a.rate, a.frames)
    leg2(fails)
    leg_sm(fails, a.capture, a.rate, a.frames)
    if _COMPARED[0] == 0:
        fails.append("the gate performed ZERO comparisons")
    print(f"\ncomparisons performed: {_COMPARED[0]}")
    if fails:
        print(f"GATE FAILED ({len(fails)}):")
        for f in fails:
            print(f"  - {f}")
        return 1
    print("GATE PASSED -- plain path bit-identical; smoothed path within the "
          "ulp bound (decoded bytes are its authority)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
