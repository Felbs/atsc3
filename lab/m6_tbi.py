#!/usr/bin/env python3
"""A/322 7.1.5.4 Twisted Block Interleaver -- the SPEC equations, gated on
A/327's printed worked example.

M5 could not read 7.1.5.4: `pdftotext` dropped every variable in the three
diagonal-read equations, so M5 enumerated 24 readings, eliminated what
structure could, and carried the twisting parameter forward as ASSUMPTION T1
to be settled by an LDPC sweep over ST = 0..36.

**THAT PREMISE WAS WRONG, and the sweep it implies is not merely
inconclusive -- it cannot contain the answer.**  There is no free integer:
7.1.5.4 DEFINES the twisting parameter, and the equations are

    R_i     = i mod N_r
    T_i     = R_i mod N_c
    C_i     = (T_i + floor(i / N_r)) mod N_c
    theta_i = N_r * C_i + R_i

with the cell skipped when `theta_i < N_FEC_TI_DUMMY * N_r`.  In M5's own
parameterisation this is ("mod", "div", "colmajor") with ST = 1 -- **the
TRANSPOSE of the reading M5 settled on.**  M5's T-GATE 4 preferred prose
("rightwards along the row") over the equations printed directly beneath that
prose, and picked the reading that advances the column fastest where the spec
advances the ROW fastest.  It matches the spec for NO value of ST, which is
why sweeping ST against the LDPC returned 37 flat failures.

**Record the lesson, it is the same one M4 recorded about the scrambler:
when prose and equations disagree, the equations win -- and "the PDF ate the
equation" is a statement about the extractor, not about the document.**

THE GATE.  A/327 Figure 6.5 prints a worked 4x3 example with one virtual FEC
Block whose expected output is `b g a f d e c h`.  `gold_vector()` reproduces
it 8 of 8.  A test vector with virtual cells exercises the skip rule, the
twist and the linear-array convention at once, and RF33 (N_virtual = 0) does
not exercise the skip rule at all -- so this gate tests strictly more than
the air does.
"""
from __future__ import annotations

import numpy as np

# A/327 Figure 6.5: Nrows 4, Ncols 3, one virtual FEC Block.
GOLD = dict(nrows=4, ncols=3, n_virtual=1, expect="b g a f d e c h".split())


def read_order(nrows, ncols, n_virtual=0):
    """A/322 7.1.5.4.  Memory indices theta_i in TBI OUTPUT order.

    The memory is written column-wise, so linear index theta = Nr*C + R
    holds cell R of FEC Block C, and Blocks 0..N_virtual-1 are virtual.
    """
    i = np.arange(nrows * ncols)
    r = i % nrows
    t = r % ncols
    c = (t + i // nrows) % ncols
    theta = nrows * c + r
    keep = theta >= n_virtual * nrows
    return theta[keep] - n_virtual * nrows


def deinterleave(cells, nrows, ncols, n_virtual=0):
    """Cells in transmitted order -> cells in memory order (FEC Block major)."""
    order = read_order(nrows, ncols, n_virtual)
    if len(order) != len(cells):
        raise ValueError(f"expected {len(order)} cells, got {len(cells)}")
    out = np.empty_like(cells)
    out[order] = cells
    return out


def interleave(mem, nrows, ncols, n_virtual=0):
    return np.asarray(mem)[read_order(nrows, ncols, n_virtual)]


def fec_block(cells, j, nrows, ncols, n_virtual=0):
    """The j-th DATA FEC Block, without materialising the whole memory."""
    r = np.arange(nrows)
    # invert: which i gives C_i = j + n_virtual and R_i = r?
    c = j + n_virtual
    i = nrows * ((c - (r % ncols)) % ncols) + r
    # i is an index into the FULL read sequence; drop the virtual skips
    if n_virtual:
        theta = nrows * ((r % ncols + i // nrows) % ncols) + r
        assert np.array_equal(theta, np.full(nrows, c) * nrows + r)
        skipped = np.searchsorted(np.sort(_virtual_positions(nrows, ncols,
                                                             n_virtual)), i)
        i = i - skipped
    return np.asarray(cells)[i]


def _virtual_positions(nrows, ncols, n_virtual):
    i = np.arange(nrows * ncols)
    r = i % nrows
    c = ((r % ncols) + i // nrows) % ncols
    return i[(nrows * c + r) < n_virtual * nrows]


# ---------------------------------------------------------------------------
# gates
# ---------------------------------------------------------------------------

def gold_vector(verbose=True):
    """A/327 Figure 6.5 -- the printed worked example, reproduced or not."""
    nr, nc, nv = GOLD["nrows"], GOLD["ncols"], GOLD["n_virtual"]
    # memory: virtual Blocks first, then data Blocks a..d, e..h
    mem = [""] * (nr * nv) + list("abcd") + list("efgh")
    got = [mem[nv * nr + t] for t in read_order(nr, nc, nv)]
    ok = got == GOLD["expect"]
    if verbose:
        print(f"    {'PASS' if ok else 'FAIL'}  A/327 Fig 6.5 gold vector "
              f"{nr}x{nc}, {nv} virtual Block: got {' '.join(got)}, "
              f"expected {' '.join(GOLD['expect'])}")
    return ok


def gate_permutation(verbose=True):
    """A bijection for every legal geometry, virtual Blocks included."""
    bad = []
    for nr in (2700, 2025, 4050, 8100, 10800, 16200):
        for nc in (1, 2, 3, 4, 5, 6, 8, 9, 12, 16, 30, 36, 37, 39, 64):
            for nv in (0, 1, nc // 3):
                if nv >= nc:
                    continue
                o = read_order(nr, nc, nv)
                if (len(o) != nr * (nc - nv)
                        or not np.array_equal(np.sort(o),
                                              np.arange(nr * (nc - nv)))):
                    bad.append((nr, nc, nv))
    if verbose:
        print(f"    {'PASS' if not bad else 'FAIL'}  permutation over "
              f"{6*15*3} legal geometries{'' if not bad else f' {bad[:4]}'}")
    return not bad


def gate_roundtrip(verbose=True):
    rng = np.random.default_rng(20260806)
    bad = []
    for nr, nc, nv in ((2700, 37, 0), (2700, 37, 5), (8100, 1, 0),
                       (8100, 39, 3), (2025, 16, 7)):
        x = rng.permutation(nr * (nc - nv))
        if not np.array_equal(deinterleave(interleave(x, nr, nc, nv),
                                           nr, nc, nv), x):
            bad.append((nr, nc, nv))
    if verbose:
        print(f"    {'PASS' if not bad else 'FAIL'}  interleave/de-interleave "
              f"round trip, with and without virtual Blocks")
    return not bad


def gate_fec_block(verbose=True):
    """fec_block(j) must equal slicing the fully de-interleaved memory."""
    rng = np.random.default_rng(7)
    bad = []
    for nr, nc, nv in ((2700, 37, 0), (8100, 39, 3), (2025, 16, 7)):
        x = rng.standard_normal(nr * (nc - nv))
        mem = deinterleave(x, nr, nc, nv)
        for j in range(nc - nv):
            if not np.array_equal(fec_block(x, j, nr, nc, nv),
                                  mem[j * nr:(j + 1) * nr]):
                bad.append((nr, nc, nv, j))
    if verbose:
        print(f"    {'PASS' if not bad else 'FAIL'}  fec_block() agrees with "
              f"the full de-interleave{'' if not bad else f' {bad[:3]}'}")
    return not bad


def selftest(verbose=True):
    if verbose:
        print("  === A/322 7.1.5.4 Twisted Block Interleaver ===")
    return all([gold_vector(verbose), gate_permutation(verbose),
                gate_roundtrip(verbose), gate_fec_block(verbose)])


if __name__ == "__main__":
    raise SystemExit(0 if selftest() else 1)
