#!/usr/bin/env python3
"""A/322 7.3 frequency interleaver / de-interleaver.

The FI is a PERMUTATION of the data cells within one OFDM symbol.  It therefore
cannot change the constellation -- only which cell sits on which carrier.  It
matters here for exactly one reason: L1-Basic occupies the FIRST N cells in
CELL order, and cell order is only recoverable by undoing this permutation.

SPEC READING AMBIGUITIES, both tested rather than assumed
---------------------------------------------------------
F1  Wire permutation DIRECTION.  Tables 7.12-7.14 print two rows, "R'i bit
    positions" and "Ri bit positions", without an explicit assignment operator.
    That admits R[row2[j]] = R'[row1[j]]  ("forward") or R[row1[j]] = R'[row2[j]]
    ("reverse").  Both are permutations, so no structural check separates them;
    only the data can.  `WIRE_DIRECTIONS` exposes both.
F2  The MSB toggle argument.  A/322 writes the toggle as "( mod 2)" with the
    variable lost to the PDF's math font; it is either the inner counter i or
    the symbol index l.  Both are exposed as `toggle`.
These are genuine ASSUMPTIONS.  They are resolved in m3_l1basic.py by scoring
the resulting cell order against a criterion that cannot be fitted (whether the
first 484 cells form a clean QPSK cluster), with every other combination acting
as its own control.
"""
from __future__ import annotations

import numpy as np

import m3_spec as S

WIRE_DIRECTIONS = ("forward", "reverse")


def _bits_to_int(bits):
    v = 0
    for j, b in enumerate(bits):
        v |= int(b) << j
    return v


def wire_map(nfft, parity, direction):
    """Return array `m` with R[j] = R'[m[j]] for j in 0..Nr-2."""
    w = S.FI_WIRE[nfft]
    src = w["rprime"]
    dst = w["even" if parity == 0 else "odd"]
    nb = len(src)
    m = np.zeros(nb, int)
    if direction == "forward":
        # R[dst[t]] = R'[src[t]]
        for t in range(nb):
            m[dst[t]] = src[t]
    else:
        # R[src[t]] = R'[dst[t]]
        for t in range(nb):
            m[src[t]] = dst[t]
    return m


def offset_word(nfft, l):
    """G_{floor(l/2)} -- A/322 7.3 symbol offset generator, Nr bits."""
    nr = int(np.log2(S.FI_NMAX[nfft]))
    g = [1] * nr                                  # G_0 = [111...1]
    taps = S.FI_G_TAPS[nfft]
    for _ in range(l // 2):
        fb = 0
        for t in taps:
            fb ^= g[t]
        g = g[1:] + [fb]                          # G_i[nr-2:0] = G_{i-1}[nr-1:1]
    return _bits_to_int(g)


def interleaving_sequence(nfft, l, n_data, direction="forward", toggle="i"):
    """H_l(p) for p = 0 .. n_data-1  (A/322 7.3).

    Returns an int array `H` where, for 8K/16K, A[p] = X[H[p]]: the cell that
    lands on the p-th available carrier is cell index H[p].
    """
    nmax = S.FI_NMAX[nfft]
    nr = int(np.log2(nmax))
    taps = S.FI_R_TAPS[nfft]
    m = wire_map(nfft, l % 2, direction)
    g = offset_word(nfft, l)

    rp = [0] * (nr - 1)                           # R'_0 = R'_1 = 0
    H = np.empty(n_data, int)
    p = 0
    for i in range(nmax):
        if i == 2:
            rp = [0] * (nr - 1)
            rp[0] = 1                             # R'_2 = [0...01]
        elif i > 2:
            fb = 0
            for t in taps:
                fb ^= rp[t]
            rp = rp[1:] + [fb]                    # R'_i[nr-3:0] = R'_{i-1}[nr-2:1]
        r = [rp[m[j]] for j in range(nr - 1)]
        tog = (i if toggle == "i" else l) & 1
        h = ((tog << (nr - 1)) | _bits_to_int(r)) ^ g
        if h < n_data:
            if p < n_data:
                H[p] = h
            p += 1
    if p != n_data:
        raise RuntimeError(f"FI produced {p} valid addresses, expected {n_data}")
    return H


def deinterleave(a, nfft, l, direction="forward", toggle="i"):
    """Given cells in CARRIER order, return them in CELL order."""
    H = interleaving_sequence(nfft, l, len(a), direction, toggle)
    x = np.empty_like(a)
    x[H] = a
    return x


def selftest():
    """The FI must be a bijection on 0..n_data-1 for every variant."""
    ok = True
    for nfft in (8192, 16384):
        for l in (0, 1, 2, 3):
            for d in WIRE_DIRECTIONS:
                for tg in ("i", "l"):
                    n = {8192: 4851, 16384: 9702}[nfft]
                    try:
                        H = interleaving_sequence(nfft, l, n, d, tg)
                        good = np.array_equal(np.sort(H), np.arange(n))
                        note = ""
                    except RuntimeError as e:
                        good, note = False, f"  ({e})"
                    if not good and tg == "i":
                        ok = False
                    print(f"  [{'PASS' if good else 'FAIL'}] FI bijection "
                          f"{nfft//1024}K l={l} {d}/{tg} n={n}{note}")
    print("\n  NOTE: toggle='l' fails structurally -- with the MSB constant it "
          "cannot\n  reach a bijection onto 0..Ndata-1.  Ambiguity F2 is "
          "therefore resolved\n  by construction, not by fitting: the toggle "
          "argument is i.")
    return ok


if __name__ == "__main__":
    print("A/322 7.3 frequency interleaver selftest\n")
    raise SystemExit(0 if selftest() else 1)
