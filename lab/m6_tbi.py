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


# ---------------------------------------------------------------------------
# A/322 7.1.5.2 -- the HTI CELL INTERLEAVER
#
# RF33 PLP 0 signals L1D_plp_HTI_cell_interleaver = 0, so this stage was
# bypassed and never written.  field-site RF8 signals 1.  Without it the
# cells inside every FEC Block are in transmitter order, the LDPC is fed a
# permuted codeword, and no amount of signal strength makes it converge --
# which is exactly the "0 of N Blocks at any SNR" this chain showed.
#
# Two generators, per the spec:
#   C_r(j) = [C_0(j) + P(r)] mod N_cells
#   C_0    -- an Nd-bit LFSR word, values >= N_cells discarded
#   P(r)   -- bit-reversed counter over Nd bits, values >= N_cells skipped
# P(r) is gated against the test vector A/322 7.1.5.2 prints in its own text
# (N_cells = 10800, N_d = 14 -> 0, 8192, 4096, 2048, 10240, 6144, 1024, 9216).
# ---------------------------------------------------------------------------

# R_i[Nd-2] feedback taps, on R_{i-1}, per Nd.  A/322 7.1.5.2.
_CI_TAPS = {11: (0, 3), 12: (0, 2), 13: (0, 1, 4, 6),
            14: (0, 1, 4, 5, 9, 11), 15: (0, 1, 2, 12)}
_CI_CACHE = {}


def _ci_nd(ncells):
    nd = int(np.ceil(np.log2(ncells)))
    if nd not in _CI_TAPS:
        raise ValueError(f"no A/322 7.1.5.2 taps for Nd={nd} "
                         f"(N_cells={ncells})")
    return nd


def cell_basic_permutation(ncells):
    """C_0(j): the basic permutation, as an int array of length `ncells`."""
    hit = _CI_CACHE.get(ncells)
    if hit is not None:
        return hit
    nd = _ci_nd(ncells)
    taps = _CI_TAPS[nd]
    out = np.empty(ncells, np.int64)
    q = 0
    prev = np.zeros(nd, np.uint8)
    for i in range(1 << nd):
        r = np.zeros(nd, np.uint8)
        if i == 2:
            r[0] = 1
        elif i > 2:
            r[0:nd - 2] = prev[1:nd - 1]
            b = 0
            for t in taps:
                b ^= int(prev[t])
            r[nd - 2] = b
        r[nd - 1] = i & 1
        v = int(np.dot(r.astype(np.int64), 1 << np.arange(nd, dtype=np.int64)))
        if v < ncells:
            if q < ncells:
                out[q] = v
            q += 1
        prev = r
    if q != ncells:
        raise ValueError(f"C_0 produced {q} values for N_cells={ncells}")
    _CI_CACHE[ncells] = out
    return out


def cell_shift(ncells, nblocks):
    """P(r) for r = 0 .. nblocks-1."""
    nd = _ci_nd(ncells)
    out, k = [], 0
    while len(out) < nblocks:
        v = int(format(k, "0%db" % nd)[::-1], 2)
        if v < ncells:
            out.append(v)
        k += 1
    return np.array(out, np.int64)


def cell_deinterleave(blocks):
    """Undo 7.1.5.2 for a (nblocks, ncells) array of one TI Block's cells.

    The transmitter reads out[j] = in[C_r(j)], so the receiver scatters:
    in[C_r(j)] = out[j].  `r` is the FEC Block index within the TI Block,
    which is why this takes the whole TI Block and not one Block at a time.
    """
    blocks = np.asarray(blocks)
    nblk, ncells = blocks.shape
    c0 = cell_basic_permutation(ncells)
    pr = cell_shift(ncells, nblk)
    out = np.empty_like(blocks)
    for r in range(nblk):
        out[r, (c0 + pr[r]) % ncells] = blocks[r]
    return out


def gate_cell_interleaver(verbose=True):
    """The spec's own printed vector, plus C_0 must be a permutation."""
    want = [0, 8192, 4096, 2048, 10240, 6144, 1024, 9216]
    got = list(cell_shift(10800, 8))
    ok_p = got == want
    if verbose:
        print(f"  P(r) N_cells=10800 Nd=14 -> {got}")
        print(f"    spec 7.1.5.2 printed  {want}   "
              f"{'PASS' if ok_p else '*** FAIL ***'}")
    ok_c = True
    for n in (8100, 10800, 5400, 16200):
        try:
            c0 = cell_basic_permutation(n)
        except ValueError as e:
            if verbose:
                print(f"  C_0 N_cells={n:6d}  skipped ({e})")
            continue
        bij = (len(c0) == n and len(np.unique(c0)) == n
               and c0.min() == 0 and c0.max() == n - 1)
        ok_c &= bij
        if verbose:
            print(f"  C_0 N_cells={n:6d}  bijection over [0,{n})  "
                  f"{'PASS' if bij else '*** FAIL ***'}")
    # round trip: de-interleaving an interleaved block returns it
    rng = np.random.default_rng(0)
    x = rng.standard_normal((5, 8100)) + 1j * rng.standard_normal((5, 8100))
    c0, pr = cell_basic_permutation(8100), cell_shift(8100, 5)
    tx = np.stack([x[r][(c0 + pr[r]) % 8100] for r in range(5)])
    rt = bool(np.allclose(cell_deinterleave(tx), x))
    if verbose:
        print(f"  round trip interleave -> deinterleave        "
              f"{'PASS' if rt else '*** FAIL ***'}")
    return ok_p and ok_c and rt


if __name__ == "__main__" and "--ci" in sys.argv:
    raise SystemExit(0 if gate_cell_interleaver() else 1)
