#!/usr/bin/env python3
"""M6 -- the DATA-PLP BICM chain, end to end, with a synthesized round trip.

This module is the receiver side of A/322 Section 6 for a DATA PLP:

    cells -> constellation demap (LLR)
          -> demultiplexer          6.3.3   (y_i,s = q_{MOD*s+i})
          -> Block de-interleaver   6.2.3   (Type A / Type B)
          -> group-wise de-int      6.2.2   (Annex B tables)
          -> parity de-interleaver  6.2.1   (Type B LDPC codes ONLY)
          -> LDPC decode            6.1.3   (the ORACLE)
          -> BCH syndrome           6.1.2.1
          -> descramble             5.2.3   (the Galois LFSR, solved in M4)
          -> Baseband Packet

WHAT IS GATED HERE AND WHAT IS NOT
----------------------------------
`selftest()` builds the TRANSMITTER from the same equations, runs random data
all the way round, and requires bit-exactness plus LDPC convergence through
AWGN.  That proves THE CODE INVERTS ITSELF.  It does NOT prove the spec
reading -- only the air can, and the air's referee is LDPC convergence with
zero unsatisfied checks.

The one thing a round trip CAN prove about a reading is discrimination: the
selftest also runs the receiver with a DELIBERATELY WRONG twisting parameter
and requires it to fail.  A sweep whose wrong hypotheses also "converge" is
not a sweep.

ASSUMPTIONS INTRODUCED HERE (all falsifiable by the air, all with controls)
    X1  A/322 6.3.3's demultiplexer, whose only definition is Figure 6.11
        (an image), is the trivial serial map y_{i,s} = q_{MOD*s+i}.
        DEFENDED, not guessed: 6.2.3.1 says "the bits from the same group in
        Part 1 shall be mapped onto the bits having the same bit position in
        each modulation symbol", and the Type A Block Interleaver is read
        row-wise with Nc == MOD columns.  A row IS a cell word, and column c
        IS bit level c.  Any other demux would break that sentence.
    X2  group-wise direction.  6.2.2 prints "Y_j = X_{pi(j)}", i.e. output
        group j takes input group pi(j).  The opposite reading is exposed as
        `gw_invert` and run as a control.
    X3  the Type A Block Interleaver equations lose their variables to the
        math font exactly like 7.1.5.4's, but unlike 7.1.5.4 the spec prints
        a WORKED EXAMPLE (256QAM, Ninner 64800: 0, 7920, 15840, ...) and
        `verify_blockint()` reproduces it.  That turns X3 from an assumption
        into a gate.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import m3_ldpc as LD                                            # noqa: E402
import m4_scrambler as SC                                       # noqa: E402
import spec_bicm as BC                                          # noqa: E402
import spec_bitint as BI                                        # noqa: E402
import spec_ldpc as SL                                          # noqa: E402
import spec_nuc as NU                                           # noqa: E402

# A/322 Table 6.2 -- Ninner = 16200, BCH outer
MOUTER_BCH_16200 = 168
KPAYLOAD_BCH_16200 = {r: SL.LDPC_PARAMS[r]["Kldpc"] - MOUTER_BCH_16200
                      for r in BI.RATES}

# A/322 Table C.1.1 -- QPSK.  k = 2*y0 + y1 (y0 = MSB).
QPSK_POINTS = np.array([1 + 1j, -1 + 1j, 1 - 1j, -1 - 1j]) / np.sqrt(2)


# ---------------------------------------------------------------------------
# 6.2 bit interleaver, as an index permutation
# ---------------------------------------------------------------------------

def parity_interleave_map(ninner, rate, explicit=False):
    """A/322 6.2.1.  Returns `lam_of_u[a]` = codeword index carried by u_a.

    **IT IS THE IDENTITY HERE, AND THAT IS A RESULT, NOT AN OVERSIGHT.**

    6.2.1 applies the parity interleaver to Type B codes only:
        u_{K+360t+s} = lam_{K+Qldpc*s+t},  0<=s<360, 0<=t<Qldpc
    while 6.1.3.2 defines the Type B ENCODER as lam_{Nouter+k} = p_k, with no
    permutation at all.  But `m3_ldpc.encode()` -- written in M3 against the
    Type A prose and reused for Type B -- applies exactly that permutation
    inside the encoder.  So m3_ldpc's codeword ordering IS the spec's u, and
    `m3_ldpc.parity_check()` builds H in the SAME ordering.  Applying 6.2.1
    on top of it would apply the permutation TWICE.

    That is verified numerically, not argued: `verify_typeb_parity()` runs
    the 6.1.3.2 encoder from scratch and shows m3_ldpc's codeword equals its
    parity-interleaved image, bit for bit.

    This was one of the two bugs that made the M5 ST sweep return 37 flat
    failures; it is invisible to any round trip, because the transmitter and
    receiver in a round trip share it.  `explicit=True` returns the 6.2.1
    permutation itself, which is kept as the CONTROL that must fail.
    """
    P = SL.LDPC_PARAMS[rate] if ninner == 16200 else SL.LDPC_PARAMS_64800[rate]
    k = P["Kldpc"]
    m = np.arange(ninner)
    if P["type"] == "A" or not explicit:
        return m
    q = P["Qldpc"] or (ninner - k) // 360
    t, s = np.divmod(np.arange(ninner - k), 360)
    m[k + 360 * t + s] = k + q * s + t
    return m


def verify_typeb_parity(rate="11/15", ninner=16200, verbose=True):
    """Gate: m3_ldpc's Type B codeword == 6.2.1(6.1.3.2 codeword)."""
    P = SL.LDPC_PARAMS[rate]
    k, q, mi = P["Kldpc"], P["Qldpc"], ninner - P["Kldpc"]
    rows = SL.LDPC_16200[rate]["rows"]
    rng = np.random.default_rng(2)
    s = rng.integers(0, 2, k).astype(np.uint8)
    # A/322 6.1.3.2 verbatim: no permutation anywhere
    p = np.zeros(mi, np.uint8)
    for kk in range(k):
        if not s[kk]:
            continue
        for x in rows[kk // 360]:
            p[(x + q * (kk % 360)) % mi] ^= 1
    np.bitwise_xor.accumulate(p, out=p)
    spec = np.concatenate([s, p])
    u = spec.copy()
    for t in range(q):
        for ss in range(360):
            u[k + 360 * t + ss] = spec[k + q * ss + t]
    ok = np.array_equal(LD.encode(s, rate, ninner), u)
    if verbose:
        print(f"    {'PASS' if ok else 'FAIL'}  m3_ldpc Type B codeword == "
              f"6.2.1 parity-interleaved 6.1.3.2 codeword ({rate})")
    return ok


def groupwise_map(ninner, mod, rate, invert=False):
    """A/322 6.2.2.  Returns `u_of_v[b]` = index into u carried by v_b.

    Y_j = X_{pi(j)}: output group j holds input group pi(j)  (ASSUMPTION X2).
    """
    perm = BI.GROUPWISE[(ninner, mod, rate)]
    ng = ninner // 360
    assert len(perm) == ng
    p = np.asarray(perm, int)
    if invert:
        p = np.argsort(p)
    off = np.arange(360)
    return (360 * p[:, None] + off[None, :]).ravel()


def blockint_map(ninner, mod, rate):
    """A/322 6.2.3.  Returns `v_of_q[j]` = index into v carried by q_j."""
    nc = BI.MOD_BITS[mod]
    typ = BI.BI_TYPE[(ninner, nc, rate)]
    q_of_v = np.empty(ninner, np.int64)
    if typ == "A":
        nr1, nr2, nc2 = BI.TYPE_A[(ninner, mod)]
        assert nc2 == nc and nr1 + nr2 == ninner // nc
        i = np.arange(nc * nr1)
        c, r = i // nr1, i % nr1
        q_of_v[:nc * nr1] = r * nc + c
        if nr2:
            i2 = np.arange(ninner - nc * nr1)
            c2, r2 = i2 // nr2, nr1 + i2 % nr2
            q_of_v[nc * nr1:] = r2 * nc + c2
    else:
        nq, np1, np2 = BI.TYPE_B[(ninner, mod)]
        assert nq == nc and np1 + np2 == ninner
        # Part 1: each set of NQCB_IG bit groups is written row-wise into a
        # (NQCB_IG x 360) array and read column-wise, 1 bit per row.
        blk = nq * 360
        nblk = np1 // blk
        row, col = np.divmod(np.arange(blk), 360)
        within = col * nq + row
        base = (np.arange(nblk) * blk)[:, None]
        q_of_v[:np1] = (base + within[None, :]).ravel()
        # Part 2: sequential, no interleaving
        q_of_v[np1:] = np.arange(np1, ninner)
    v_of_q = np.empty(ninner, np.int64)
    v_of_q[q_of_v] = np.arange(ninner)
    return v_of_q


def bicm_map(ninner, mod, rate, gw_invert=False):
    """Composite: `lam_of_q[j]` = LDPC codeword bit carried by q_j."""
    lam_of_u = parity_interleave_map(ninner, rate)
    u_of_v = groupwise_map(ninner, mod, rate, gw_invert)
    v_of_q = blockint_map(ninner, mod, rate)
    return lam_of_u[u_of_v[v_of_q]]


# ---------------------------------------------------------------------------
# 6.3 constellation mapping / demapping
# ---------------------------------------------------------------------------

def points_for(mod, rate):
    if mod == "QPSK":
        return QPSK_POINTS
    return NU.get({"16QAM": 16, "64QAM": 64, "256QAM": 256}[mod], rate)


def bit_masks(nbits):
    """masks[i] -> boolean array over the 2**nbits labels, True where y_i = 1.

    y_0 is the MSB (A/322 6.3.4.2).
    """
    k = np.arange(1 << nbits)
    return np.array([((k >> (nbits - 1 - i)) & 1).astype(bool)
                     for i in range(nbits)])


def demap_llr(cells, mod, rate, sigma2=None):
    """Max-log LLRs, q-stream order (ASSUMPTION X1).  llr > 0 means bit 0."""
    pts = points_for(mod, rate)
    nb = BI.MOD_BITS[mod]
    z = np.asarray(cells)
    d2 = np.abs(z[:, None] - pts[None, :]) ** 2
    if sigma2 is None:                       # min-sum is scale invariant, but
        sigma2 = max(float(np.mean(d2.min(1))), 1e-9)   # keep the units sane
    ones = bit_masks(nb)
    out = np.empty((len(z), nb))
    for i in range(nb):
        out[:, i] = (d2[:, ones[i]].min(1) - d2[:, ~ones[i]].min(1)) / sigma2
    return out.ravel()                        # q_{MOD*s+i}


def map_cells(bits, mod, rate):
    """Transmitter side of the same map, for the round trip."""
    pts = points_for(mod, rate)
    nb = BI.MOD_BITS[mod]
    b = np.asarray(bits, np.uint8).reshape(-1, nb)
    k = np.zeros(len(b), int)
    for i in range(nb):
        k |= b[:, i].astype(int) << (nb - 1 - i)
    return pts[k]


# ---------------------------------------------------------------------------
# 6.1.2.1 BCH  (Mouter = 168 for Ninner = 16200)
# ---------------------------------------------------------------------------

def _gen_poly(polys=None):
    g = 1
    for p in (polys or BC.BCH_16200_GEN_POLYS):
        q = sum(1 << e for e in p)
        r, a = 0, g
        while q:
            if q & 1:
                r ^= a
            a <<= 1
            q >>= 1
        g = r
    return g


BCH_GEN = _gen_poly()
BCH_DEG = BCH_GEN.bit_length() - 1

# M10: Ninner = 64800 uses a GF(2^16) BCH -- twelve degree-16 factors, so
# Mouter = 192, not 168.  Both edition witnesses give the same twelve polys
# (spec_bicm's own gates), and the degree sum IS the closed identity: 12x16.
BCH_GEN_BY_N = {16200: BCH_GEN, 64800: _gen_poly(BC.BCH_64800_GEN_POLYS)}
BCH_DEG_BY_N = {n: g.bit_length() - 1 for n, g in BCH_GEN_BY_N.items()}
MOUTER_BCH = {16200: 168, 64800: 192}
assert BCH_DEG_BY_N == MOUTER_BCH, BCH_DEG_BY_N


def _poly_mod(v, ninner=16200):
    g, gd = BCH_GEN_BY_N[ninner], BCH_DEG_BY_N[ninner]
    while v.bit_length() - 1 >= gd and v:
        v ^= g << ((v.bit_length() - 1) - gd)
    return v


def bch_syndrome(bits, ninner=16200):
    v = 0
    for b in np.asarray(bits, np.uint8):
        v = (v << 1) | int(b)
    return _poly_mod(v, ninner)


def bch_encode(info_bits, ninner=16200):
    v = 0
    for b in np.asarray(info_bits, np.uint8):
        v = (v << 1) | int(b)
    deg = BCH_DEG_BY_N[ninner]
    p = _poly_mod(v << deg, ninner)
    par = np.array([(p >> (deg - 1 - i)) & 1 for i in range(deg)], np.uint8)
    return np.concatenate([np.asarray(info_bits, np.uint8), par])


# ---------------------------------------------------------------------------
# one FEC block, receive side
# ---------------------------------------------------------------------------

class PlpChain:
    """Everything a single data PLP needs, built once and reused."""

    #: Normalized-min-sum attenuation ladder, tried in order until the syndrome
    #: is zero.  **This is not a knob that can manufacture a decode**: the
    #: oracle is still 0 unsatisfied checks out of tens of thousands plus a
    #: zero BCH syndrome, and no value of alpha makes a wrong codeword pass.
    #:
    #: M10 found this the expensive way.  alpha = 0.75 is right for the
    #: Ninner = 16200 codes RF33 carries (check degrees 4..11) and is WRONG for
    #: the 64800 rate-6/15 code RF25/RF30 carry, whose check degrees are only
    #: 6..8 -- min-sum's overestimate shrinks with degree, so 0.75
    #: over-attenuates and the decoder parks in a trapping set at 4..16
    #: unsatisfied of 38880 and STAYS THERE at 300 iterations.  That looks
    #: exactly like a link 0.2 dB short.  It is not: at alpha = 1.0 the same
    #: cells converge in 16..43 iterations with zero unsatisfied.
    ALPHA_LADDER = (1.0, 0.85, 0.75)

    def __init__(self, ninner, mod, rate, gw_invert=False, iters=50,
                 alpha=None):
        self.alpha = alpha
        self.ninner, self.mod, self.rate = ninner, mod, rate
        self.mod_bits = BI.MOD_BITS[mod]
        self.cells_per_fec = BI.DEMUX[(ninner, mod)][1]
        self.kldpc = (SL.LDPC_PARAMS[rate]["Kldpc"] if ninner == 16200
                      else SL.LDPC_PARAMS_64800[rate]["Kldpc"])
        self.mouter = MOUTER_BCH[ninner]
        self.kpayload = self.kldpc - self.mouter
        self.lam_of_q = bicm_map(ninner, mod, rate, gw_invert)
        self.checks = LD.parity_check(rate, ninner)
        self.idx, self.mask = LD.pack(self.checks, ninner)
        self.iters = iters

    # -- receive ----------------------------------------------------------
    def llr(self, cells, sigma2=None):
        q = demap_llr(cells, self.mod, self.rate, sigma2)
        lam = np.empty(self.ninner)
        lam[self.lam_of_q] = q
        return lam

    def decode(self, cells, sigma2=None, iters=None, alpha=None):
        lam = self.llr(cells, sigma2)
        return self.decode_llr(lam, iters=iters, alpha=alpha)

    def decode_llr(self, lam, iters=None, alpha=None):
        ladder = ((alpha or self.alpha,) if (alpha or self.alpha)
                  else self.ALPHA_LADDER)
        best = None
        for a in ladder:
            bits, ok, it, bad = LD.min_sum_decode(
                lam, self.checks, self.ninner, iters=iters or self.iters,
                alpha=a)
            r = dict(bits=bits, converged=bool(ok), iters=int(it),
                     unsatisfied=int(bad), alpha=float(a))
            if ok:
                return r
            if best is None or bad < best["unsatisfied"]:
                best = r
        return best

    def baseband_packet(self, cw_bits):
        """Codeword -> BCH syndrome + descrambled Baseband Packet bytes."""
        nouter = cw_bits[:self.kldpc]
        syn = bch_syndrome(nouter, self.ninner)
        payload = SC.descramble(nouter[:self.kpayload])
        return dict(bch_syndrome=int(syn), bch_ok=(syn == 0),
                    bits=payload,
                    bytes=np.packbits(payload).tobytes())

    # -- transmit (round trip only) ---------------------------------------
    def encode(self, payload_bits):
        assert len(payload_bits) == self.kpayload
        scr = np.asarray(payload_bits, np.uint8) ^ SC.sequence(self.kpayload)
        nouter = bch_encode(scr, self.ninner)
        cw = LD.encode(nouter, self.rate, self.ninner)
        q = cw[self.lam_of_q]
        return map_cells(q, self.mod, self.rate)


# ---------------------------------------------------------------------------
# gates
# ---------------------------------------------------------------------------

def verify_blockint(verbose=True):
    """X3: reproduce A/322 6.2.3.1's printed worked example, 256QAM 64800.

    The spec prints the first 16 and last 16 output indices explicitly.  The
    Type A equations lost their variables to the math font, so this example
    is the only thing that can arbitrate the reading -- and it is a gate,
    not an assumption.
    """
    # rate 2/15 is the Type A entry of Table 6.8 for 256QAM/64800 (11/15 is
    # Type B); the block interleaver TYPE is per-rate, the CONFIGURATION
    # (Table 6.10) is not, so the printed example is a Type A example.
    v_of_q = blockint_map(64800, "256QAM", "2/15")
    head = [0, 7920, 15840, 23760, 31680, 39600, 47520, 55440,
            1, 7921, 15841, 23761, 31681, 39601, 47521, 55441]
    tail = [63539, 63719, 63899, 64079, 64259, 64439, 64619, 64799]
    ok1 = list(v_of_q[:16]) == head
    ok2 = list(v_of_q[-8:]) == tail
    mid = [63359, 63360, 63540, 63720, 63900, 64080, 64260, 64440, 64620]
    ok3 = list(v_of_q[8 * 7920 - 1:8 * 7920 + 8]) == mid
    if verbose:
        print(f"    {'PASS' if ok1 else 'FAIL'}  6.2.3.1 example, first 16 "
              f"output indices")
        print(f"    {'PASS' if ok3 else 'FAIL'}  6.2.3.1 example, the Part1 / "
              f"Part2 seam (63359, 63360, 63540, ...)")
        print(f"    {'PASS' if ok2 else 'FAIL'}  6.2.3.1 example, last 8 "
              f"output indices")
    return ok1 and ok2 and ok3


def verify_perm(verbose=True, ninner_want=16200):
    """Every composite BICM map must be a permutation of 0..Ninner-1."""
    bad, n = [], 0
    for (ninner, mod, rate) in sorted(BI.GROUPWISE):
        if ninner != ninner_want:
            continue
        if BI.DEMUX[(ninner, mod)][1] is None:
            continue
        n += 1
        m = bicm_map(ninner, mod, rate)
        if not np.array_equal(np.sort(m), np.arange(ninner)):
            bad.append((ninner, mod, rate))
    if verbose:
        print(f"    {'PASS' if not bad else 'FAIL'}  composite BICM map is a "
              f"permutation for all {n} Ninner={ninner_want} combinations"
              f"{'' if not bad else f' -- {bad[:4]}'}")
    return not bad


def roundtrip(mod, rate, snr_db=None, seed=1, iters=50, verbose=True,
              ninner=16200):
    """Synthesized end-to-end: payload -> cells -> payload, plus AWGN."""
    rng = np.random.default_rng(seed)
    ch = PlpChain(ninner, mod, rate, iters=iters)
    pay = rng.integers(0, 2, ch.kpayload).astype(np.uint8)
    cells = ch.encode(pay)
    if snr_db is not None:
        n = (rng.standard_normal(len(cells))
             + 1j * rng.standard_normal(len(cells))) / np.sqrt(2)
        cells = cells + n * np.sqrt(10 ** (-snr_db / 10))
    r = ch.decode(cells)
    bb = ch.baseband_packet(r["bits"])
    same = bool(np.array_equal(bb["bits"], pay))
    if verbose:
        print(f"    {'PASS' if r['converged'] and bb['bch_ok'] and same else 'FAIL'}"
              f"  {mod:6s} {rate:5s} "
              f"{'noiseless' if snr_db is None else f'{snr_db:.1f} dB'}: "
              f"LDPC {'converged' if r['converged'] else 'FAILED'} in "
              f"{r['iters']:2d} iters ({r['unsatisfied']} unsatisfied), "
              f"BCH {'0' if bb['bch_ok'] else 'NONZERO'}, "
              f"payload {'exact' if same else 'WRONG'}")
    return r["converged"] and bb["bch_ok"] and same


def selftest(verbose=True):
    ok = True
    if verbose:
        print("  === table gates (the modules this chain stands on) ===")
    p, b = BI.verify(verbose=False)
    if verbose:
        print(f"    {'PASS' if not b else 'FAIL'}  spec_bitint: {p} closed "
              f"identities, {b} fail")
    ok &= (b == 0)
    if verbose:
        print("\n  === 6.2.3.1 worked example (ASSUMPTION X3 -> a gate) ===")
    ok &= verify_blockint(verbose)
    if verbose:
        print("\n  === structural ===")
    ok &= verify_perm(verbose)
    ok &= verify_perm(verbose, 64800)
    if verbose:
        print("\n  === synthesized round trip: payload -> cells -> payload ===")
    for mod, rate in (("QPSK", "2/15"), ("64QAM", "11/15"),
                      ("256QAM", "11/15"), ("16QAM", "8/15")):
        ok &= roundtrip(mod, rate, None, verbose=verbose)
    if verbose:
        print("\n  === Ninner = 64800 (M10: the codes RF25/RF30 carry) ===")
    for mod, rate in (("QPSK", "6/15"), ("QPSK", "9/15"),
                      ("64QAM", "6/15"), ("256QAM", "7/15")):
        ok &= roundtrip(mod, rate, None, verbose=verbose, ninner=64800)
    if verbose:
        print("\n  === through AWGN, at roughly the measured air MER ===")
    ok &= roundtrip("QPSK", "2/15", 3.0, verbose=verbose)
    ok &= roundtrip("64QAM", "11/15", 19.0, verbose=verbose)
    if verbose:
        print("\n  === CONTROLS -- these MUST fail ===")
    ch = PlpChain(16200, "64QAM", "11/15")
    rng = np.random.default_rng(9)
    pay = rng.integers(0, 2, ch.kpayload).astype(np.uint8)
    cells = ch.encode(pay)
    for name, c in (("cells rotated by one", np.roll(cells, 1)),
                    ("cells conjugated", np.conj(cells)),
                    ("random cells",
                     cells[rng.permutation(len(cells))])):
        r = ch.decode(c, iters=25)
        good = not r["converged"]
        ok &= good
        if verbose:
            print(f"    {'PASS' if good else 'FAIL'}  {name:24s} -> "
                  f"{'did NOT converge' if not r['converged'] else 'CONVERGED (!!)'}"
                  f"  ({r['unsatisfied']} unsatisfied)")
    gw = PlpChain(16200, "64QAM", "11/15", gw_invert=True)
    r = gw.decode(cells, iters=25)
    good = not r["converged"]
    ok &= good
    if verbose:
        print(f"    {'PASS' if good else 'FAIL'}  {'inverted group-wise (X2)':24s}"
              f" -> {'did NOT converge' if not r['converged'] else 'CONVERGED (!!)'}"
              f"  ({r['unsatisfied']} unsatisfied)")
    dbl = PlpChain(16200, "64QAM", "11/15")
    dbl.lam_of_q = parity_interleave_map(16200, "11/15", explicit=True)[dbl.lam_of_q]
    r = dbl.decode(cells, iters=25)
    good = not r["converged"]
    ok &= good
    if verbose:
        print(f"    {'PASS' if good else 'FAIL'}  6.2.1 applied TWICE (the M5"
              f" bug) -> {'did NOT converge' if not r['converged'] else 'CONVERGED (!!)'}"
              f"  ({r['unsatisfied']} unsatisfied)")
    return ok


if __name__ == "__main__":
    print("M6 -- DATA-PLP BICM chain selftest")
    print("=" * 72)
    good = selftest()
    print(f"\n  {'ALL SELFTESTS PASS' if good else 'SELFTEST FAILURES'}")
    raise SystemExit(0 if good else 1)
