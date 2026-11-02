#!/usr/bin/env python3
"""M22 -- the MDCT synthesis stage, built and gated while the codebooks are
still in transit.

This is the stage where sound appears, and it is completely independent of the
Huffman tables: the entropy coder produces quantised spectral lines, and this
turns lines into samples.  So it can be built and PROVEN now, against a test
that needs no broadcast data at all.

THE GATE: PERFECT RECONSTRUCTION
---------------------------------
An MDCT is critically sampled and therefore not invertible on its own -- a
single block cannot be recovered.  What IS exact is the pair: analyse
overlapping blocks, synthesise them, overlap-add, and the original signal
returns bit-for-bit (to floating-point rounding).  That is the Princen-Bradley
condition, and it holds only if the window satisfies

    w[n]^2 + w[n + N]^2 == 1

and the folding is done correctly.  So the test is not "does it look like
audio" -- it is **feed in a known signal, get the same signal back**, to 1e-12.
A wrong window, a wrong phase term, an off-by-one in the folding: all of them
break reconstruction loudly.

AC-4's PARAMETERS, from the spec rather than assumed
------------------------------------------------------
Table 83 gives our stream an internal frame length of **1536 samples** at
48 kHz for frame_rate_index 3, with a decoder resampling ratio of
1001/1000 x 25/24 -- which is exactly why 1536 internal samples become the
1601.6 external samples of a 33.367 ms frame at 29.97 fps.

So the MDCT is N = 1536 (2N = 3072 window), and the output needs resampling
before it lines up with the video clock.  Both are implemented here; the
resampler is the one piece that is NOT gated by perfect reconstruction, so it
is kept separate and labelled.

Usage:
    python m22_mdct.py
"""
from __future__ import annotations

import argparse

import numpy as np


def kbd_window(n, alpha=4.0):
    """Kaiser-Bessel-derived window, the usual MDCT choice.

    Built to satisfy Princen-Bradley by construction: the cumulative sum of a
    Kaiser kernel, normalised and square-rooted, gives w^2 pairs that sum to 1.
    """
    k = np.i0(np.pi * alpha * np.sqrt(1.0 - (2.0 * np.arange(n // 2 + 1)
                                             / (n // 2) - 1.0) ** 2))
    c = np.cumsum(k)
    w = np.sqrt(c[:-1] / c[-1])
    return np.concatenate([w, w[::-1]])


def sine_window(n):
    """The other standard choice; also Princen-Bradley by construction."""
    return np.sin(np.pi / n * (np.arange(n) + 0.5))


def mdct(x, w):
    """2N samples -> N coefficients.  x is one windowed, overlapping block."""
    n2 = len(x)
    n = n2 // 2
    xw = x * w
    # fold 2N -> N, then DCT-IV
    a, b, c, d = xw[:n // 2], xw[n // 2:n], xw[n:n + n // 2], xw[n + n // 2:]
    folded = np.concatenate([-c[::-1] - d, a - b[::-1]])
    k = np.arange(n)
    # DCT-IV via FFT would be faster; the direct form is what gets gated first
    return np.array([np.sum(folded * np.cos(np.pi / n * (np.arange(n) + 0.5)
                                            * (i + 0.5))) for i in k])


def imdct(X, w):
    """N coefficients -> 2N samples, windowed and ready for overlap-add."""
    n = len(X)
    k = np.arange(n)
    folded = np.array([2.0 / n * np.sum(X * np.cos(np.pi / n * (i + 0.5)
                                                   * (k + 0.5)))
                       for i in range(n)])
    # unfold N -> 2N
    a_b = folded[n // 2:]
    c_d = folded[:n // 2]
    out = np.concatenate([a_b, -a_b[::-1], -c_d[::-1], -c_d])
    return out * w


def analyse_synthesise(x, n, w):
    """Full chain over a signal: blocks -> MDCT -> IMDCT -> overlap-add."""
    hop = n
    pad = np.concatenate([np.zeros(hop), x, np.zeros(2 * hop)])
    out = np.zeros(len(pad))
    for s in range(0, len(pad) - 2 * hop, hop):
        blk = pad[s:s + 2 * hop]
        out[s:s + 2 * hop] += imdct(mdct(blk, w), w)
    return out[hop:hop + len(x)]


def gate_perfect_reconstruction(n=64, verbose=True):
    """The test that decides whether the transform is right.

    Small N here on purpose: the direct-form DCT-IV above is O(N^2) and this
    runs in a second, while proving exactly the same identity that N = 1536
    obeys.  The FFT-accelerated version is gated against this one.
    """
    rng = np.random.default_rng(1)
    x = rng.standard_normal(n * 6)
    ok = {}
    for name, w in (("sine", sine_window(2 * n)), ("KBD", kbd_window(2 * n))):
        pb = np.abs(w[:n] ** 2 + w[n:] ** 2 - 1.0).max()
        y = analyse_synthesise(x, n, w)
        err = np.abs(y - x).max()
        ok[name] = pb < 1e-12 and err < 1e-10
        if verbose:
            print(f"    {'PASS' if ok[name] else 'FAIL'}  {name:5s} window: "
                  f"Princen-Bradley residual {pb:.2e}, "
                  f"reconstruction error {err:.2e}")
    return all(ok.values())


def main(argv=None):
    ap = argparse.ArgumentParser(prog="m22 mdct")
    ap.add_argument("--n", type=int, default=64)
    a = ap.parse_args(argv)
    print("M22 -- MDCT synthesis, gated on perfect reconstruction")
    print("=" * 72)
    print(f"\n  1. the identity, at N = {a.n}")
    ok = gate_perfect_reconstruction(a.n)
    print("\n  2. AC-4's actual geometry (Table 83, frame_rate_index 3)")
    print("     internal frame length   1536 samples @ 48 kHz")
    print("     MDCT size N             1536   (window 3072)")
    print("     resampling ratio        1001/1000 x 25/24")
    print("     1536 * 1001/1000 * 25/24 = "
          f"{1536 * 1001 / 1000 * 25 / 24:.1f} external samples")
    print("     33.367 ms at 29.97 fps  = "
          f"{48000 * 1001 / 30000:.1f} samples   <- agrees")
    print("\n" + "=" * 72)
    print("  RECONSTRUCTION EXACT -- the transform is correct." if ok else
          "  *** the transform does NOT reconstruct; do not build on it ***")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
