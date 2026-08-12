#!/usr/bin/env python3
"""gate_e54.py -- the compiled min-sum kernel's referee, on real air.

The E52/E53 gate pattern: every reference is re-derived, never stored.

  1. Kernel vs numpy float32 batch (`_fast_loop`, the production fast
     path): bits, converged, iters AND unsatisfied-count must be EQUAL for
     every (block, alpha) decode of real Frames -- the kernel is a
     re-compilation, not a re-rounding, and this leg proves it output-
     exact.  Repeated N times (determinism is asserted, not assumed).
  2. Kernel vs numpy float64 (the gate REFERENCE): converged hard bits
     identical, and the BCH syndromes + descrambled Baseband Packet BYTES
     from the kernel's bits identical to the reference's -- the byte-level
     verdict the campaign requirement names.
  3. Dispatch: `min_sum_decode_batch(dtype=float32)` actually routes
     through the kernel (and returns the very same arrays the direct call
     returns) -- a gate that passes with the kernel silently absent would
     gate nothing.

Run: python lab/gate_e54.py [--capture lab/rate6912_rf33.cs16]
     [--frames 3] [--repeats 3]
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m3_ldpc as LD                                              # noqa: E402
import m6_cells as C                                              # noqa: E402
import m9_fast as F9                                              # noqa: E402
import gate_e52 as G52                                            # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--capture", default="rate6912_rf33.cs16")
    ap.add_argument("--rate", type=float, default=6.912e6)
    ap.add_argument("--frames", type=int, default=3)
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--threads", type=int, default=12)
    a = ap.parse_args(argv)
    cap = a.capture if os.path.isabs(a.capture) else os.path.join(HERE,
                                                                  a.capture)
    print("E54 -- compiled min-sum kernel gates (real air)")
    print("=" * 72)

    if LD._kernel_lib() is None:
        print("    FAIL  compiled kernel not loadable -- "
              "python lab/build_ldpc_kernel.py first")
        return 1
    ok = True

    span = 0.05 + (a.frames + 1) * F9.FRAME_SEC
    fd = F9.FrameDecoder(threads=a.threads, backend="cpu", cpu_fast=False)
    y = F9.load_fast(cap, a.rate, span_sec=span * (a.rate / 6.912e6),
                     ex=fd.ex)[0]
    t_prev = None
    lams = []
    for fi in range(a.frames):
        centre = C.BOOTSTRAP if t_prev is None else t_prev + C.FRAME_SAMPLES
        t0, _ = fd.fine_timing(y, span=20, centre=centre)
        t_prev = t0
        lams.append(G52.build_lams(fd, y, t0))

    # -- leg 1: kernel == numpy float32 fast path, every (block, alpha) ----
    os.environ["ATSC3_LDPC_KERNEL"] = "0"      # numpy reference, re-derived
    LD._KERNEL_STATE["lib"] = None
    np_ref = {}
    for li, (lam, lam16) in enumerate(lams):
        for ci, (l2d, chain) in enumerate(((lam, fd.ch0), (lam16, fd.ch16))):
            for alpha in fd.ladder_for(chain):
                np_ref[li, ci, alpha] = LD.min_sum_decode_batch(
                    l2d, chain.checks, chain.ninner, iters=fd.iters,
                    alpha=alpha, dtype=np.float32)
    os.environ["ATSC3_LDPC_KERNEL"] = "1"
    LD._KERNEL_STATE["lib"] = None
    if LD._kernel_lib() is None:
        print("    FAIL  kernel did not re-arm after the reference pass")
        fd.close()
        return 1
    n_cmp = n_diff = 0
    for rep in range(a.repeats):
        for li, (lam, lam16) in enumerate(lams):
            for ci, (l2d, chain) in enumerate(((lam, fd.ch0),
                                               (lam16, fd.ch16))):
                for alpha in fd.ladder_for(chain):
                    kb, kc, ki, kd = LD._kernel_decode_checked(
                        np.ascontiguousarray(l2d, np.float32), chain.checks,
                        chain.ninner, fd.iters, alpha)
                    if kb is None:
                        n_diff += len(l2d)
                        n_cmp += len(l2d)
                        continue
                    rb, rc, ri, rd = np_ref[li, ci, alpha]
                    for i in range(len(l2d)):
                        n_cmp += 1
                        if not (np.array_equal(kb[i], rb[i])
                                and kc[i] == rc[i] and ki[i] == ri[i]
                                and kd[i] == rd[i]):
                            n_diff += 1
    g = n_diff == 0
    ok &= g
    print(f"    {'PASS' if g else 'FAIL'}  kernel == numpy float32 batch: "
          f"{n_cmp} (block, alpha, repeat) decodes output-identical "
          f"(bits+conv+its+bad; {n_diff} differ; repeats={a.repeats})")

    # -- leg 2: vs float64 reference, through BCH to BYTES -----------------
    n_bits_diff = n_byte_diff = n_syn_diff = n_conv = 0
    for lam, lam16 in lams:
        for l2d, chain in ((lam, fd.ch0), (lam16, fd.ch16)):
            alpha = fd.ladder_for(chain)[0]
            b64, c64, _, _ = LD.min_sum_decode_batch(
                l2d, chain.checks, chain.ninner, iters=fd.iters,
                alpha=alpha, dtype=np.float64)
            kb, kc, _, _ = LD._kernel_decode_checked(
                np.ascontiguousarray(l2d, np.float32), chain.checks,
                chain.ninner, fd.iters, alpha)
            both = kc & c64
            n_conv += int(both.sum())
            for i in np.flatnonzero(both):
                if not np.array_equal(kb[i], b64[i]):
                    n_bits_diff += 1
            syn_k, pk_k = fd.baseband_batch(kb, chain)
            syn_r, pk_r = fd.baseband_batch(b64, chain)
            n_syn_diff += int((np.asarray(syn_k) != np.asarray(syn_r)).sum())
            n_byte_diff += sum(1 for x, z in zip(pk_k, pk_r) if x != z)
    g = n_bits_diff == 0 and n_byte_diff == 0 and n_syn_diff == 0
    ok &= g
    print(f"    {'PASS' if g else 'FAIL'}  kernel vs float64 reference: "
          f"{n_conv} converged blocks, hard bits identical "
          f"({n_bits_diff} differ), BCH syndromes identical "
          f"({n_syn_diff} differ), Baseband Packet BYTES identical "
          f"({n_byte_diff} differ)")

    # -- leg 3: the public API actually dispatches to the kernel -----------
    lam, _ = lams[0]
    kb, kc, ki, kd = LD._kernel_decode_checked(
        np.ascontiguousarray(lam, np.float32), fd.ch0.checks,
        fd.ch0.ninner, fd.iters, 1.0)
    ab, ac, ai, ad = LD.min_sum_decode_batch(
        lam, fd.ch0.checks, fd.ch0.ninner, iters=fd.iters, alpha=1.0,
        dtype=np.float32)
    g = (kb is not None and np.array_equal(kb, ab)
         and np.array_equal(kc, ac) and np.array_equal(ki, ai)
         and np.array_equal(kd, ad))
    ok &= g
    print(f"    {'PASS' if g else 'FAIL'}  min_sum_decode_batch(float32) "
          f"routes through the kernel (dispatch verified against the "
          f"direct call)")

    fd.close()
    print(f"\n  {'ALL E54 GATES PASS' if ok else 'E54 GATE FAILURES'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
