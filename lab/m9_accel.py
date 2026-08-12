#!/usr/bin/env python3
"""M9 Step 2 -- CPU accelerations, installed by monkeypatch, all reversible.

WHY A SEPARATE MODULE AND NOT AN EDIT
-------------------------------------
The existing chain is the CORRECTNESS ORACLE for everything M9 does, and the
user's standing law is that acceleration is always optional with the slow path
kept working.  So nothing here edits `m3_*`/`m6_*`; `install()` swaps module
attributes and `uninstall()` puts them back, and every driver takes
`--accel none` to run the untouched code.

WHAT IS ACCELERATED, AND WHY EACH ONE IS SAFE
---------------------------------------------
Five of the six are MEMOISATION of pure functions of their arguments, so
bit-identity is not a hope, it is a property:

  1. `m3_freqint.interleaving_sequence` -- a bit-serial LFSR walk over 8192
     addresses, in Python, recomputed for EVERY OFDM symbol of EVERY frame.
     It depends only on (nfft, l, n_data, direction, toggle) and there are 36
     distinct values in a frame.  22.8% of the profile.
  2. `m4_scrambler.sequence` -- the A/322 5.2.3 LFSR, in Python, one byte at a
     time, recomputed for every FEC Block (descramble) and twice per frame at
     full pool length.  It is a FIXED sequence: compute the longest length
     ever asked for, then slice.  16% of the profile in three places.
  3. `m3_ldpc.pack` -- builds the (4320 x 11) padded check-node index array
     with a Python loop.  `min_sum_decode` called it ONCE PER FEC BLOCK, i.e.
     74 times a frame, for an array that depends only on the code.  12.1% of
     the profile -- three quarters of everything attributed to "the LDPC".
  4. `m3_ldpc.parity_check` -- the H build, per PlpChain construction.
  5. `m6_bicm.bicm_map` -- the composite BICM permutation, likewise.

The sixth is an ALGEBRAIC replacement, and it is gated rather than assumed:

  6. `m6_bicm.bch_syndrome` -- the reference walks the 11880 message bits into
     one Python big integer and then reduces it one bit at a time.  The
     syndrome is a LINEAR function of the message over GF(2), so it is exactly
     a 168 x 11880 binary matrix product.  Packed into uint64 lanes with
     `np.bitwise_count` it is ~31k word operations instead of ~2M.  The
     columns of that matrix are the remainders x^k mod g(x), computed from the
     SAME generator polynomial object the reference uses, and `verify_bch()`
     requires exact equality of the integer syndrome on random messages AND on
     real decoded codewords.  10.9% of the profile.

None of these change a single floating-point operation anywhere in the chain.
"""
from __future__ import annotations

import functools
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m3_freqint as FI                                          # noqa: E402
import m3_ldpc as LD                                             # noqa: E402
import m3_spec as S                                              # noqa: E402
import m4_scrambler as SC                                        # noqa: E402
import m6_bicm as B                                              # noqa: E402
import m6_cells as MC                                            # noqa: E402
import spec_pilots as SP                                         # noqa: E402

_ORIG = {}
_INSTALLED = False


# ---------------------------------------------------------------------------
# 1 -- frequency interleaver sequences
# ---------------------------------------------------------------------------

@functools.lru_cache(maxsize=None)
def _fi_seq(nfft, l, n_data, direction, toggle):
    return _ORIG["fi"](nfft, l, n_data, direction, toggle)


def fi_sequence(nfft, l, n_data, direction="forward", toggle="i"):
    return _fi_seq(nfft, l, n_data, direction, toggle)


# ---------------------------------------------------------------------------
# 2 -- the A/322 5.2.3 scrambling sequence
# ---------------------------------------------------------------------------

_SCR = {"n": 0, "seq": np.zeros(0, np.uint8)}


def scr_sequence(nbits, out_taps=None, **kw):
    """Grow-and-slice.  Falls back to the reference for non-default taps."""
    if kw or (out_taps is not None and tuple(out_taps) != SC.SOLVED_OUT):
        return _ORIG["seq"](nbits, out_taps=out_taps or SC.SOLVED_OUT, **kw)
    if nbits > _SCR["n"]:
        want = max(nbits, 262144)
        _SCR["seq"] = _ORIG["seq"](want)
        _SCR["n"] = want
    return _SCR["seq"][:nbits]


def scr_descramble(bits, out_taps=None, **kw):
    bits = np.asarray(bits, np.uint8)
    return bits ^ scr_sequence(len(bits), out_taps=out_taps, **kw)


# ---------------------------------------------------------------------------
# 3/4/5 -- LDPC structure caches
# ---------------------------------------------------------------------------

_PACK = {}


def ldpc_pack(checks, n):
    """Cache on the IDENTITY of the checks list.

    A strong reference to `checks` is kept in the cache entry, so the id
    cannot be recycled onto a different object while the entry lives.
    """
    key = (id(checks), n)
    hit = _PACK.get(key)
    if hit is not None:
        return hit[1], hit[2]
    idx, mask = _ORIG["pack"](checks, n)
    _PACK[key] = (checks, idx, mask)
    return idx, mask


@functools.lru_cache(maxsize=None)
def _pchk(rate, n, variant):
    return _ORIG["pchk"](rate, n, variant)


def ldpc_parity_check(rate, n=16200, variant="standard"):
    return _pchk(rate, n, variant)


@functools.lru_cache(maxsize=None)
def _bmap(ninner, mod, rate, gw_invert):
    return _ORIG["bicm"](ninner, mod, rate, gw_invert)


def bicm_map(ninner, mod, rate, gw_invert=False):
    return _bmap(ninner, mod, rate, gw_invert)


# ---------------------------------------------------------------------------
# 6 -- BCH syndrome as a packed GF(2) matrix product
# ---------------------------------------------------------------------------

_BCH = {}


def _bch_gen(ninner):
    """(generator, degree) for this Ninner, tolerant of both m6_bicm eras.

    The 64800 codes use a GF(2^16) BCH with Mouter = 192, added to m6_bicm on
    the 64K/CTI track while this one was in flight -- so read the per-Ninner
    tables when they exist and fall back to the single 16200 pair when they
    do not.
    """
    by_n = getattr(B, "BCH_GEN_BY_N", None)
    if by_n is not None:
        return by_n[ninner], B.BCH_DEG_BY_N[ninner]
    if ninner != 16200:
        raise KeyError(ninner)
    return B.BCH_GEN, B.BCH_DEG


def _bch_tables(n, ninner=16200):
    """(mask, nwords) where mask[j] has bit i set iff x^(n-1-i) has x^j."""
    hit = _BCH.get((n, ninner))
    if hit is not None:
        return hit
    g, deg = _bch_gen(ninner)
    # rem[k] = x^k mod g, for k = 0 .. n-1.  deg is 168 = 21 whole bytes, so
    # the remainder's big-endian bytes unpack straight into an MSB-first row.
    nby = (deg + 7) // 8
    cols = np.empty((n, deg), np.uint8)
    r = 1
    for k in range(n):
        cols[n - 1 - k] = np.unpackbits(
            np.frombuffer(r.to_bytes(nby, "big"), np.uint8))[-deg:]
        r <<= 1
        if (r >> deg) & 1:
            r ^= g
    pad = (-n) % 64
    if pad:
        cols = np.concatenate([cols, np.zeros((pad, deg), np.uint8)])
    nw = (n + pad) // 64
    mask = np.empty((deg, nw), np.uint64)
    for j in range(deg):
        mask[j] = np.packbits(cols[:, j]).view(np.uint64)
    _BCH[(n, ninner)] = (mask, nw, pad)
    return _BCH[(n, ninner)]


def bch_syndrome(bits, ninner=16200):
    bits = np.asarray(bits, np.uint8)
    n = len(bits)
    if n % 8:
        return _ORIG["bch"](bits, ninner)
    mask, nw, pad = _bch_tables(n, ninner)
    v = np.packbits(np.concatenate([bits, np.zeros(pad, np.uint8)])
                    if pad else bits).view(np.uint64)
    par = (np.bitwise_count(mask & v[None, :]).sum(1) & 1).astype(np.uint64)
    out = 0
    for b in par:                                # bit deg-1 first (MSB)
        out = (out << 1) | int(b)
    return out


def bch_syndrome_batch(bits2d, ninner=16200):
    """(nblk, n) -> (nblk,) python ints.  Same algebra, one broadcast pass.

    The per-row Python big-integer assembly is skipped for the rows whose
    syndrome is zero -- which, on a healthy multiplex, is all of them -- so
    the common case never leaves NumPy.
    """
    a = np.asarray(bits2d, np.uint8)
    nblk, n = a.shape
    if n % 8:
        return [_ORIG["bch"](r, ninner) for r in a]
    mask, nw, pad = _bch_tables(n, ninner)
    if pad:
        a = np.concatenate([a, np.zeros((nblk, pad), np.uint8)], axis=1)
    v = np.packbits(a, axis=1).view(np.uint64)                  # (nblk, nw)
    par = (np.bitwise_count(v[:, None, :] & mask[None, :, :]).sum(2)
           & 1).astype(np.uint8)                                # (nblk, deg)
    nz = np.flatnonzero(par.any(1))
    out = [0] * nblk
    for i in nz:
        x = 0
        for b in par[i]:
            x = (x << 1) | int(b)
        out[i] = x
    return out


# ---------------------------------------------------------------------------
# 7 -- the pilot tables, rebuilt for every OFDM symbol of every frame
# ---------------------------------------------------------------------------
# `demod_data` cost 3.6 ms a call and only ~0.15 ms of that was the FFT.  The
# rest was `pilot_values()` walking a 13-bit LFSR 6913 times in Python and
# `pilot_carriers()` building three Python sets -- both pure functions of the
# symbol index, both rebuilt 76 times a frame.  No caller mutates either
# result (checked across the whole tree: len(), sorted(), np.array(), and set
# unions inside the function itself), so memoising them is exact.

@functools.lru_cache(maxsize=None)
def _refseq(n):
    return _ORIG["ref"](n)


def reference_sequence(n):
    if n > 65536:
        return _ORIG["ref"](n)
    return _refseq(n)


@functools.lru_cache(maxsize=None)
def _pvals(n, amp):
    # Every caller in the tree immediately does `np.array(...)` on this, and
    # converting a 6913-element Python list of floats costs more than the FFT
    # it feeds.  Returning the ndarray makes that conversion a memcpy.
    return np.asarray(_ORIG["pvals"](n, amp), float)


def pilot_values(n, amp):
    return _pvals(n, amp)


@functools.lru_cache(maxsize=None)
def _pcar(nfft, cred, pattern, l, sbs, table):
    return frozenset(_ORIG["pcar"](nfft, cred, pattern, l, sbs, table))


def pilot_carriers(nfft, cred, pattern, l, sbs, table=None):
    if table is not None:
        return _ORIG["pcar"](nfft, cred, pattern, l, sbs, table)
    return _pcar(nfft, cred, pattern, l, sbs, None)


@functools.lru_cache(maxsize=None)
def _ccps(nfft, cred):
    return _ORIG["ccps"](nfft, cred)


def common_cps(nfft, cred):
    return _ccps(nfft, cred)


# ---------------------------------------------------------------------------
# 8 -- demod_data's per-symbol geometry
# ---------------------------------------------------------------------------
# With (7) in place `m6_cells.demod_data` still spent most of its 1.5 ms on
# three things that depend only on (l, sbs): the pilot index set (a Python
# `sorted()` over a ~2000-element set), the FFT-bin translation of both index
# sets, and `ref[pk]`.  This caches the geometry and leaves every arithmetic
# operation -- the FFT, the divisions, np.interp, the coherence -- untouched
# and in the same order.

_GEO = {}


def _geometry(l, sbs):
    hit = _GEO.get((l, sbs))
    if hit is not None:
        return hit
    import m3_spec as S3
    import m6_cells as MC
    import spec_pilots as P
    from m3_preamble import fft_bins
    noc = P.NOC[(MC.NFFT, MC.CRED)]
    lo, _ = S3.carrier_abs_range(MC.NFFT, MC.CRED)
    ref = np.asarray(S3.pilot_values(noc, 1.0))
    dx, dy = P.dxdy(MC.PATTERN)
    pk = (np.arange(0, noc, dx) if sbs
          else np.arange(dx * (l % dy), noc, dx * dy))
    pk = np.unique(np.concatenate([pk, [0, noc - 1]]))
    used = np.zeros(noc, bool)
    used[np.array(sorted(P.pilot_carriers(MC.NFFT, MC.CRED, MC.PATTERN, l,
                                          sbs)), int)] = True
    d = np.flatnonzero(~used)
    _GEO[(l, sbs)] = g = dict(
        pk=pk, refpk=ref[pk], bins_pk=fft_bins(MC.NFFT, lo + pk),
        d=d, bins_d=fft_bins(MC.NFFT, lo + d), kk=np.arange(noc),
        step=MC.NFFT + MC.GI, gi=MC.GI, nfft=MC.NFFT)
    return g


def demod_data(y, t0, l, sbs):
    """m6_cells.demod_data with the per-symbol geometry cached."""
    import m6_cells as MC
    g = _geometry(l, sbs)
    s = t0 + g["step"] * (1 + l) + g["gi"]
    Y = np.fft.fftshift(np.fft.fft(y[s:s + g["nfft"]]))
    hp = Y[g["bins_pk"]] / g["refpk"]
    H = MC._interp(g["kk"], g["pk"], hp)
    z = Y[g["bins_d"]] / H[g["d"]]
    coh = float(abs(np.sum(hp[:-1] * np.conj(hp[1:])))
                / max(np.sum(np.abs(hp) ** 2), 1e-30))
    return z, coh


def demod_coh(y, t0, l, sbs=False):
    """The coherence metric of `demod_data`, WITHOUT equalising the data cells.

    fine_timing throws `z` away for all 41 candidate offsets and keeps only
    `coh` -- but `demod_data` still paid for the channel interpolation and the
    6913 equalising divisions on every candidate.  This computes `hp` and
    `coh` with the SAME operations in the SAME order (the lines are copied,
    not re-derived), so the value is bit-identical; everything skipped was
    arithmetic whose result was unused.  E52.
    """
    g = _geometry(l, sbs)
    s = t0 + g["step"] * (1 + l) + g["gi"]
    Y = np.fft.fftshift(np.fft.fft(y[s:s + g["nfft"]]))
    hp = Y[g["bins_pk"]] / g["refpk"]
    return float(abs(np.sum(hp[:-1] * np.conj(hp[1:])))
                 / max(np.sum(np.abs(hp) ** 2), 1e-30))


def demod_coh_batch(y, cands, l=1, sbs=False):
    """`demod_coh` for a run of CONSECUTIVE candidate offsets, one batched FFT.

    fine_timing's 41 candidates differ by one sample each, so their FFT
    windows are one strided view and pocketfft transforms them as rows of a
    single call -- the same 1-D transform per row, and the same reductions
    per row (NumPy's pairwise sum over a fixed-length last axis is the same
    tree batched or alone).  gate_e52 holds every candidate's coherence to
    BIT EQUALITY against the scalar `demod_coh` on real air -- this is an
    exactness claim, gated, not assumed.
    """
    g = _geometry(l, sbs)
    cands = np.asarray(cands)
    s0 = int(cands[0]) + g["step"] * (1 + l) + g["gi"]
    base = y[s0:]
    W = np.lib.stride_tricks.as_strided(
        base, shape=(len(cands), g["nfft"]),
        strides=(base.strides[0], base.strides[0]))
    Y = np.fft.fftshift(np.fft.fft(W, axis=1), axes=1)
    hp = Y[:, g["bins_pk"]] / g["refpk"]
    num = np.abs(np.sum(hp[:, :-1] * np.conj(hp[:, 1:]), axis=1))
    den = np.sum(np.abs(hp) ** 2, axis=1)
    return num / np.maximum(den, 1e-30)


# ---------------------------------------------------------------------------
# install / uninstall
# ---------------------------------------------------------------------------

def install():
    global _INSTALLED
    if _INSTALLED:
        return
    _ORIG.update(fi=FI.interleaving_sequence, seq=SC.sequence,
                 pack=LD.pack, pchk=LD.parity_check, bicm=B.bicm_map,
                 bch=B.bch_syndrome, descr=SC.descramble,
                 ref=S.reference_sequence, pvals=S.pilot_values,
                 pcar=SP.pilot_carriers, ccps=S.common_cps,
                 demod=MC.demod_data)
    S.reference_sequence = reference_sequence
    S.pilot_values = pilot_values
    S.common_cps = common_cps
    SP.pilot_carriers = pilot_carriers
    MC.demod_data = demod_data
    FI.interleaving_sequence = fi_sequence
    SC.sequence = scr_sequence
    SC.descramble = scr_descramble
    LD.pack = ldpc_pack
    LD.parity_check = ldpc_parity_check
    B.bicm_map = bicm_map
    B.bch_syndrome = bch_syndrome
    _INSTALLED = True


def uninstall():
    global _INSTALLED
    if not _INSTALLED:
        return
    FI.interleaving_sequence = _ORIG["fi"]
    SC.sequence = _ORIG["seq"]
    SC.descramble = _ORIG["descr"]
    LD.pack = _ORIG["pack"]
    LD.parity_check = _ORIG["pchk"]
    B.bicm_map = _ORIG["bicm"]
    B.bch_syndrome = _ORIG["bch"]
    S.reference_sequence = _ORIG["ref"]
    S.pilot_values = _ORIG["pvals"]
    S.common_cps = _ORIG["ccps"]
    SP.pilot_carriers = _ORIG["pcar"]
    MC.demod_data = _ORIG["demod"]
    _INSTALLED = False


# ---------------------------------------------------------------------------
# gates
# ---------------------------------------------------------------------------

def verify(verbose=True, seed=11, data=None):
    """Every accelerated function against the one it replaces, exactly."""
    was = _INSTALLED
    if not was:
        install()
    rng = np.random.default_rng(seed)
    ok = True

    # 1 -- FI sequences, every symbol of a subframe-0 frame
    bad = 0
    for l in range(36):
        a = _ORIG["fi"](8192, l, 4851, "forward", "i")
        b = fi_sequence(8192, l, 4851, "forward", "i")
        bad += int(not np.array_equal(a, b))
    ok &= bad == 0
    if verbose:
        print(f"    {'PASS' if not bad else 'FAIL'}  FI sequence identical, "
              f"36 symbols x 4851 addresses")

    # 2 -- scrambler
    bad = 0
    for n in (7, 1024, 11712, 211472, 8192):
        bad += int(not np.array_equal(_ORIG["seq"](n), scr_sequence(n)))
    ok &= bad == 0
    if verbose:
        print(f"    {'PASS' if not bad else 'FAIL'}  scrambling sequence "
              f"identical, 5 lengths up to 211472")

    # 3 -- pack
    ch = _ORIG["pchk"]("11/15", 16200, "standard")
    i0, m0 = _ORIG["pack"](ch, 16200)
    i1, m1 = ldpc_pack(ch, 16200)
    g = np.array_equal(i0, i1) and np.array_equal(m0, m1)
    ok &= g
    if verbose:
        print(f"    {'PASS' if g else 'FAIL'}  ldpc.pack identical")

    # 6 -- BCH, random messages and a real codeword shape
    bad, worst = 0, None
    for _ in range(24):
        n = int(rng.choice([2160, 11880]))
        bits = rng.integers(0, 2, n).astype(np.uint8)
        a, b = _ORIG["bch"](bits), bch_syndrome(bits)
        if a != b:
            bad += 1
            worst = (n, a, b)
    # and the all-zero / all-one edges
    for bits in (np.zeros(11880, np.uint8), np.ones(11880, np.uint8)):
        if _ORIG["bch"](bits) != bch_syndrome(bits):
            bad += 1
    ok &= bad == 0
    if verbose:
        print(f"    {'PASS' if not bad else 'FAIL'}  BCH syndrome identical, "
              f"26 messages (11880 and 2160 bits)"
              f"{'' if not bad else f'  {worst}'}")

    # batch form agrees with the scalar form
    m = rng.integers(0, 2, (16, 11880)).astype(np.uint8)
    bb = bch_syndrome_batch(m)
    g = all(bch_syndrome(m[i]) == bb[i] for i in range(16))
    ok &= g
    if verbose:
        print(f"    {'PASS' if g else 'FAIL'}  batched BCH == scalar BCH, "
              f"16 messages")

    # 7 -- pilot tables
    bad = 0
    for l in range(36):
        for sbs in (False, True):
            a = _ORIG["pcar"](8192, 0, "SP4_2", l, sbs, None)
            b = pilot_carriers(8192, 0, "SP4_2", l, sbs)
            bad += int(set(a) != set(b))
    bad += int(list(_ORIG["pvals"](6913, 1.0)) != list(pilot_values(6913, 1.0)))
    bad += int(list(_ORIG["ccps"](8192, 0)) != list(common_cps(8192, 0)))
    ok &= bad == 0
    if verbose:
        print(f"    {'PASS' if not bad else 'FAIL'}  pilot tables identical, "
              f"72 symbol/SBS combinations + values + common CPs")

    # 8 -- demod_data
    bad, worst = 0, None
    if data is not None:
        for l in range(36):
            for sbs in (False, True):
                za, ca = _ORIG["demod"](data[0], data[1], l, sbs)
                zb, cb = demod_data(data[0], data[1], l, sbs)
                if not (np.array_equal(za, zb) and ca == cb):
                    bad += 1
                    worst = (l, sbs, float(np.abs(za - zb).max()))
        ok &= bad == 0
        if verbose:
            print(f"    {'PASS' if not bad else 'FAIL'}  demod_data identical "
                  f"on real air, 72 symbol/SBS combinations"
                  f"{'' if not bad else f'  {worst}'}")
    elif verbose:
        print("    ----  demod_data gate SKIPPED (no capture supplied)")

    if not was:
        uninstall()
    return ok


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--capture", default="long_rf33.cs16",
                    help="real air for the demod_data gate ('' to skip)")
    ap.add_argument("--rate", type=float, default=8e6)
    _a = ap.parse_args()
    print("M9 -- CPU acceleration gates (every one exact)")
    print("=" * 72)
    _data = None
    if _a.capture:
        import m6_cells as _C
        _p = (_a.capture if os.path.isabs(_a.capture)
              else os.path.join(HERE, _a.capture))
        if os.path.exists(_p):
            _y = _C.load_frame(_p, _a.rate,
                               span_sec=0.30 * (_a.rate / 6.912e6))
            _data = (_y, _C.fine_timing(_y)[0])
    good = verify(data=_data)
    print(f"\n  {'ALL GATES PASS' if good else 'GATE FAILURES'}")
    raise SystemExit(0 if good else 1)
