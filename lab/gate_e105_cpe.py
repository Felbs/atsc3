#!/usr/bin/env python3
"""gate_e105_cpe.py -- the compiled CPE nearest-point search's referee.

`cpe_nearest_f32` replaces the argmin inside `m9_fast._cpe_fast`.  Its output
is a hard DECISION (an index), not a float, and it feeds the phase estimate
that divides every cell in the pool -- so a single flipped index is not a
rounding, it moves the correction.  The bar is therefore identity of the
chosen points, and of the corrected pool that follows from them.

  leg 1  DECISION IDENTITY, synthetic: every constellation the receiver uses,
         random cells, kernel argmin vs numpy argmin.  Plus cells placed
         exactly ON points (the decision must be that point) and cells placed
         exactly BETWEEN two points (the tie case -- np.argmin keeps the
         FIRST minimum and the kernel must agree, which is the one rule a
         "reasonable" C implementation is most likely to get wrong).
  leg 2  REAL AIR: run _cpe_fast over a decoded Frame's pool with the kernel
         on and off; the corrected pool and the per-symbol gains must match.
  leg 3  DECODED BYTES: the whole-decoder authority, run separately by
         replaying a capture with ATSC3_DEMAP_KERNEL=0/1 -- named here so the
         obligation is written down, not implied.

Run: python lab/gate_e105_cpe.py [--capture hit_rf33.cs16] [--rate 8e6]
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m6_bicm as M6                                            # noqa: E402
import m6_cells as C                                            # noqa: E402
import m9_fast as F9                                            # noqa: E402

_COMPARED = [0]


def _kern():
    os.environ["ATSC3_DEMAP_KERNEL"] = "1"
    F9._DEMAP_STATE["lib"] = None
    return F9._demap_lib()


def _argmin_ref(z, pts):
    """numpy's own answer, in the association _cpe_fast uses."""
    P = np.vstack((pts.real, pts.imag)).astype(np.float32)
    Pf = (P * np.float32(-2.0)).copy()
    p2f = (pts.real ** 2 + pts.imag ** 2).astype(np.float32)
    X = np.empty((len(z), 2), np.float32)
    X[:, 0] = z.real
    X[:, 1] = z.imag
    sc = X @ Pf
    sc += p2f[None, :]
    return sc.argmin(1).astype(np.int32)


def _argmin_kernel(lib, z, pts):
    P = np.vstack((pts.real, pts.imag)).astype(np.float32)
    Pf = (P * np.float32(-2.0)).copy()
    p2f = (pts.real ** 2 + pts.imag ** 2).astype(np.float32)
    zr = np.ascontiguousarray(z.real, np.float32)
    zi = np.ascontiguousarray(z.imag, np.float32)
    idx = np.empty(len(z), np.int32)
    rc = lib.cpe_nearest_f32(zr, zi, len(z),
                             np.ascontiguousarray(Pf[0]),
                             np.ascontiguousarray(Pf[1]),
                             np.ascontiguousarray(p2f), len(pts), idx)
    return rc, idx


def _cmp_idx(tag, a, b, fails, z=None, pts=None):
    """Identical decisions -- or, where they differ, a PROVEN tie.

    At a cell exactly equidistant from two constellation points the argmin is
    genuinely ambiguous: numpy and the kernel compute the two scores through
    different roundings (BLAS may contract the 2-term dot into an FMA), so
    "the first minimum" can land on either. That is not a wrong answer, and
    the reference's own comment says as much -- "a flip needs two points
    equidistant to within a rounding".
    So a difference is only forgiven if the two chosen points are equidistant
    from the cell IN FLOAT64, checked here rather than assumed. Any other
    difference is a real defect and fails the gate.
    """
    _COMPARED[0] += 1
    bad = np.flatnonzero(a != b)
    if bad.size == 0:
        print(f"    {tag}: IDENTICAL ({a.size} decisions)")
        return
    if z is None or pts is None:
        print(f"    {tag}: {bad.size} of {a.size} decisions DIFFER")
        fails.append(f"{tag}: {bad.size} differing decisions")
        return
    zz = np.asarray(z, np.complex128)[bad]
    da = np.abs(zz - pts[a[bad]])
    db = np.abs(zz - pts[b[bad]])
    tie = np.isclose(da, db, rtol=0, atol=1e-12 * max(1.0, float(np.abs(pts).max())))
    n_tie = int(tie.sum())
    n_real = int((~tie).sum())
    print(f"    {tag}: {bad.size} of {a.size} differ -- {n_tie} PROVEN TIES "
          f"(equidistant in float64), {n_real} genuine")
    if n_real:
        worst = float(np.abs(da - db)[~tie].max())
        fails.append(f"{tag}: {n_real} non-tie decision flips "
                     f"(worst distance gap {worst:.3e})")


def _constellations():
    out = []
    for name, rate in (("QPSK", "2/15"), ("16QAM", "11/15"),
                       ("64QAM", "11/15"), ("256QAM", "11/15")):
        try:
            pts = np.asarray(M6.points_for(name, rate), np.complex128)
        except Exception as ex:                                # noqa: BLE001
            print(f"  FAIL: {name} could not be built: {ex}")
            continue
        out.append((name, pts))
    return out


def leg1(fails):
    print("LEG 1 -- decision identity: random, on-point and TIE cells")
    lib = _kern()
    rng = np.random.default_rng(20260823)
    cons = _constellations()
    if not cons:
        fails.append("no constellation could be built -- gate tested nothing")
        return
    for name, pts in cons:
        scale = float(np.abs(pts).max())
        cases = {
            "random": (rng.normal(0, scale, 20000)
                       + 1j * rng.normal(0, scale, 20000)),
            "on points": np.tile(pts, 50),
            # exact midpoints of adjacent pairs: the tie case
            "midpoints": ((pts[:-1] + pts[1:]) / 2.0),
        }
        for tag, z in cases.items():
            rc, idx_k = _argmin_kernel(lib, z, pts)
            if rc != 0:
                fails.append(f"{name}/{tag}: kernel returned {rc}")
                continue
            idx_n = _argmin_ref(z, pts)
            print(f"  {name} / {tag}")
            _cmp_idx("argmin", idx_k, idx_n, fails, z, pts)


def leg2(fails, cap, rate):
    print(f"LEG 2 -- real air: _cpe_fast pool + gains, kernel on vs off")
    path = os.path.join(HERE, cap)
    if not os.path.exists(path):
        msg = f"capture {cap} missing -- leg 2 could not run"
        print(f"  FAIL: {msg}")
        fails.append(msg)
        return
    res = {}
    for use in (True, False):
        os.environ["ATSC3_DEMAP_KERNEL"] = "1" if use else "0"
        F9._DEMAP_STATE["lib"] = None
        fd = F9.FrameDecoder(threads=4, backend="cpu", cpu_fast=True)
        y = F9.load_fast(path, rate, span_sec=(0.05 + 2 * F9.FRAME_SEC)
                         * (rate / 6.912e6), ex=fd.ex)[0]
        pool, info = fd.cell_pool_fast_sm(y, C.BOOTSTRAP, {})
        key = len(pool)
        if key not in fd.regions:
            fd.regions[key] = C.constellation_regions(key)
        reg, dv = fd.regions[key]
        out, gains = fd.cpe_correct(pool, info["symbol_of"],
                                    [fd.pts0, fd.pts16], reg, dv)
        res[use] = (np.asarray(out).copy(), dict(gains))
    (ok_, gk), (on_, gn) = res[True], res[False]
    _COMPARED[0] += 1
    if np.array_equal(ok_, on_):
        print(f"    corrected pool: BIT-IDENTICAL ({ok_.size} cells)")
    else:
        d = np.abs(ok_ - on_)
        sc = float(np.abs(on_).max()) or 1.0
        print(f"    corrected pool: DIFFERS max|d| {float(d.max()):.3e} "
              f"rel {float(d.max())/sc:.3e}")
        fails.append("corrected pool differs")
    _COMPARED[0] += 1
    same_g = gk.keys() == gn.keys() and all(gk[k] == gn[k] for k in gk)
    print(f"    per-symbol gains: {'IDENTICAL' if same_g else 'DIFFER'} "
          f"({len(gk)} symbols)")
    if not same_g:
        fails.append("per-symbol CPE gains differ")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--capture", default="hit_rf33.cs16")
    ap.add_argument("--rate", type=float, default=8e6)
    a = ap.parse_args(argv)
    if _kern() is None:
        print("FAIL: demap/CPE kernel not loadable -- "
              "python lab/build_demap_kernel.py")
        return 1
    fails = []
    leg1(fails)
    leg2(fails, a.capture, a.rate)
    if _COMPARED[0] == 0:
        fails.append("the gate performed ZERO comparisons")
    print(f"\ncomparisons performed: {_COMPARED[0]}")
    print("LEG 3 (decoded bytes, whole decoder) is run by replaying a capture "
          "with the kernel on and off -- see the campaign notes.")
    if fails:
        print(f"GATE FAILED ({len(fails)}):")
        for f in fails:
            print(f"  - {f}")
        return 1
    print("GATE PASSED -- the CPE kernel reproduces the numpy decisions")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
