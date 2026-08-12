#!/usr/bin/env python3
"""M16 -- the 256QAM demapper as a matrix multiply, exactly.

The M13 follow-up measured the wall for live 1080p and found it was NOT the
LDPC (which fits inside a Frame on the GPU) but the DEMAPPER: PLP 1 is 947,700
cells against 256 constellation points where PLP 0 is 199,800 against 64, so
27x the work, ~619 ms against a 247.11 ms Frame.

WHAT THIS FOUND, INCLUDING THE TWO WRONG GUESSES
------------------------------------------------
Three hypotheses were tried against the same benchmark.  TWO WERE WRONG, and
they are kept here because the wrong ones are the useful part:

  1. "the arithmetic is the cost, so make it a GEMM"      1.2x.  Nearly nothing.
  2. "no -- it is the fancy-index COPIES in the per-bit
     minima, 16 copies of a ~1 GB array"                  SLOWER.  Rewriting
     them as strided views over a reshape made it worse, not better.
  3. "it is simply single-threaded"                       10.8x, and that is
     the whole answer.  Cache-blocking the chunk size from 947,700 down to
     8,000 changed the time by 1%, which rules out a memory-bandwidth story
     outright; 12 threads took it from 15.3 s to 1.70 s.

So this is M9's lesson a third time -- **the parallelism has to be inside the
stage** -- and the GEMM is kept only because it is exact, tidy, and free, not
because it is what made the difference.

THE EXACT TRANSFORM (worth keeping regardless of what it bought)
----------------------------------------------------------------
Candidate pruning is the usual remedy and it is an APPROXIMATION that would
have to be gated against the exhaustive demapper.  An exact transform exists
and should be taken first:

    |y - p|^2  =  |y|^2  -  2*Re(y conj(p))  +  |p|^2
                  ^^^^^
                  common to EVERY point, and max-log LLR is a DIFFERENCE of
                  two such distances -- so this term cancels identically.

What is left is `|p|^2 - 2*Re(y conj(p))`, and the second term is
`[Re y, Im y] @ [Re p, Im p]^T`: a real GEMM instead of the
broadcast-subtract-magnitude over a complex array 256 columns wide.  Exact,
and worth 1.2x on its own -- which is to say, worth having and not worth
believing in.

WHY THIS IS EXACT AND NOT MERELY CLOSE
---------------------------------------
`demap_llr` returns `(min over points with bit=1) - (min over points with
bit=0)`, divided by sigma2.  Adding the same constant `|y|^2` to every distance
shifts both minima equally, so the difference is unchanged -- ALGEBRAICALLY,
not to a tolerance.  What does change is floating-point summation order, so the
gate here is not "close enough": it is (a) the LLR difference against the
exhaustive demapper, and (b) the thing that actually matters, **the decoded
bits, which must be identical**.

One caveat that would have made this silently wrong: `demap_llr` ESTIMATES
sigma2 from `d2.min(1)` when the caller does not supply one, and that estimate
uses the absolute distances -- which this transform does not preserve.  So the
absolute minimum is recovered by adding `|y|^2` back before the estimate.  The
cancellation is exploited only where it is valid.

Usage:
    python m16_demap_fast.py            # gate + benchmark
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m6_bicm as B                                              # noqa: E402
import spec_bitint as BI                                         # noqa: E402


def demap_llr_fast(cells, mod, rate, sigma2=None, dtype=np.float64, ex=None,
                   chunk=20_000):
    """Drop-in for `m6_bicm.demap_llr`, computed as a GEMM.

    `ex` is a ThreadPoolExecutor and it is where nearly all of the speed is --
    see the module docstring.  Without one this runs single-threaded and is
    barely faster than the exhaustive version.
    """
    if ex is not None and len(cells) > chunk:
        nb = BI.MOD_BITS[mod]
        parts = list(ex.map(
            lambda lo: demap_llr_fast(cells[lo:lo + chunk], mod, rate,
                                      sigma2, dtype),
            range(0, len(cells), chunk)))
        return np.concatenate(parts)
    pts = B.points_for(mod, rate)
    nb = BI.MOD_BITS[mod]
    z = np.asarray(cells)
    yy = np.empty((len(z), 2), dtype)
    yy[:, 0], yy[:, 1] = z.real, z.imag
    pp = np.empty((2, len(pts)), dtype)
    pp[0, :], pp[1, :] = pts.real, pts.imag
    # m == |y-p|^2 - |y|^2, and the missing term is common to every column
    m = (np.abs(pts) ** 2).astype(dtype)[None, :] - 2.0 * (yy @ pp)
    if sigma2 is None:
        # the caller's estimator needs ABSOLUTE distances, so put |y|^2 back
        sigma2 = max(float(np.mean(m.min(1) + (np.abs(z) ** 2))), 1e-9)
    out = np.empty((len(z), nb), dtype)
    # THE ACTUAL BOTTLENECK, and it was not the arithmetic.  `m[:, ones[i]]` is
    # fancy indexing, so each of the 16 subset minima COPIES half a
    # (cells x 256) array -- ~15 GB of copying for one Frame, which is why the
    # GEMM alone bought almost nothing.  `bit_masks` labels point k with
    # bit i = (k >> (nb-1-i)) & 1, so each subset is a STRIDE, not a scatter:
    # reshaping to (cells, 2^i, 2, 2^(nb-1-i)) puts bit i on its own axis and
    # the two minima fall out of one reduction over a view, with no copy at all.
    for i in range(nb):
        v = m.reshape(len(z), 1 << i, 2, 1 << (nb - 1 - i))
        mn = v.min(axis=3).min(axis=1)            # (cells, 2): bit=0, bit=1
        out[:, i] = (mn[:, 1] - mn[:, 0]) / sigma2
    return out.ravel().astype(np.float64)


def gate(mod="256QAM", rate="11/15", n=200_000, seed=3, verbose=True):
    """Exhaustive vs GEMM: the LLRs, and then the DECODED BITS."""
    rng = np.random.default_rng(seed)
    pts = B.points_for(mod, rate)
    z = (pts[rng.integers(0, len(pts), n)]
         + 0.06 * (rng.standard_normal(n) + 1j * rng.standard_normal(n)))
    ok = {}
    from concurrent.futures import ThreadPoolExecutor
    exg = ThreadPoolExecutor(max_workers=12)
    for name, s2, xx in (("sigma2 supplied", 0.0036, None),
                         ("sigma2 estimated", None, None),
                         ("sigma2 supplied, 12 threads", 0.0036, exg)):
        a = B.demap_llr(z, mod, rate, sigma2=s2)
        b = demap_llr_fast(z, mod, rate, sigma2=s2, ex=xx)
        same = bool(np.array_equal(a, b))
        d = float(np.abs(a - b).max())
        rel = d / max(float(np.abs(a).max()), 1e-30)
        ok[name] = same or rel < 1e-12
        if verbose:
            print(f"    {'PASS' if ok[name] else 'FAIL'}  {name}: "
                  f"array_equal={same}, max |delta| {d:.3e} "
                  f"(relative {rel:.2e})")
    return ok


def gate_decoded(mod="256QAM", rate="11/15", ninner=64800, snr_db=22.0,
                 seed=4, verbose=True):
    """The gate that matters: same LLR path, same decoded codeword."""
    rng = np.random.default_rng(seed)
    ch = B.PlpChain(ninner, mod, rate, iters=50)
    pay = rng.integers(0, 2, ch.kpayload).astype(np.uint8)
    cells = ch.encode(pay)
    nz = (rng.standard_normal(len(cells))
          + 1j * rng.standard_normal(len(cells))) / np.sqrt(2)
    y = cells + nz * np.sqrt(10 ** (-snr_db / 10.0))
    res = []
    for name, fn in (("exhaustive", B.demap_llr),
                     ("GEMM", demap_llr_fast)):
        q = fn(y, mod, rate)
        lam = np.empty(ch.ninner)
        lam[ch.lam_of_q] = q
        r = ch.decode_llr(lam)
        bb = ch.baseband_packet(r["bits"])
        res.append((name, r["converged"], bb["bch_ok"],
                    np.array_equal(bb["bits"], pay), r["bits"]))
    same = bool(np.array_equal(res[0][4], res[1][4]))
    if verbose:
        for name, cv, bch, exact, _ in res:
            print(f"    {name:11s} converged={cv} bch_ok={bch} "
                  f"payload_exact={exact}")
        print(f"    {'PASS' if same else 'FAIL'}  decoded codewords identical: "
              f"{same}")
    return same and all(r[3] for r in res)


def bench(mod="256QAM", rate="11/15", n=947_700, seed=5, verbose=True):
    rng = np.random.default_rng(seed)
    pts = B.points_for(mod, rate)
    z = (pts[rng.integers(0, len(pts), n)]
         + 0.06 * (rng.standard_normal(n) + 1j * rng.standard_normal(n)))
    out = {}
    from concurrent.futures import ThreadPoolExecutor
    ex = ThreadPoolExecutor(max_workers=12)
    for name, fn in (("exhaustive", lambda: B.demap_llr(z, mod, rate)),
                     ("GEMM f64", lambda: demap_llr_fast(z, mod, rate)),
                     ("GEMM f32", lambda: demap_llr_fast(z, mod, rate,
                                                         dtype=np.float32)),
                     ("f32 x12", lambda: demap_llr_fast(
                         z, mod, rate, dtype=np.float32, ex=ex))):
        t = time.time()
        fn()
        out[name] = time.time() - t
        if verbose:
            print(f"    {name:11s} {out[name]:7.3f} s   "
                  f"{out[name]/0.24711:6.2f}x the 247.11 ms Frame")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(prog="m16 demap fast")
    ap.add_argument("--cells", type=int, default=947_700)
    a = ap.parse_args(argv)
    print("M16 -- the 256QAM demapper as a GEMM")
    print("=" * 72)
    print("\n  1. LLRs vs the exhaustive demapper")
    g1 = gate()
    print("\n  2. the gate that matters -- decoded codewords")
    g2 = gate_decoded()
    print(f"\n  3. one PLP 1 Frame ({a.cells:,} cells x 256 points)")
    b = bench(n=a.cells)
    sp = b["exhaustive"] / b["GEMM f32"]
    print(f"\n    speedup, exhaustive -> GEMM f32: {sp:.1f}x")
    ok = all(g1.values()) and g2
    print("\n" + "=" * 72)
    print(f"  {'ALL GATES PASS' if ok else '*** GATES FAILED ***'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
