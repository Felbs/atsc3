#!/usr/bin/env python3
"""M3 Step 3 -- L1-Basic bits from real air: demap, de-interleave, de-puncture,
LDPC decode, BCH check.

No radio.  Offline only.

THE REFEREE
-----------
Everything up to here has been judged by statistics that a careful person could
argue about.  This stage has an arbiter that cannot be argued with: the
rate-3/15 LDPC code either converges to a codeword with all 12960 parity checks
satisfied, or it does not; and the 168-bit BCH syndrome over the resulting 368
bits is either zero or it is not.  Neither can be reached by tuning.  A wrong
table, a wrong permutation direction or a wrong puncturing schedule produces a
loud failure, not a plausible-looking answer.

That is exactly why the remaining spec-reading ambiguities are enumerated here
and RUN, rather than chosen:
    * group-wise parity permutation direction  (Y_j = X_pi(j) vs its inverse)
    * LDPC parity interleave index order       (ASSUMPTION L1 in m3_ldpc.py)
The air arbitrates, with LDPC convergence + BCH syndrome as the referee.  The
control is that the wrong variants must FAIL -- if several variants "work",
the test is worthless and that gets reported.

THE CHAIN, inverted (A/322 6.5.2, L1-Basic Mode 3)
--------------------------------------------------
  484 QPSK cells
    -> 968 LLRs                     Annex C.1.1 mapping, y1 = sign(I), y0 = sign(Q)
    -> bit de-interleave            6.5.2.10, 484 rows x 2 cols, write col, read row
    -> de-puncture / de-shorten     6.5.2.8 / 6.5.2.4, into a 16200-bit LLR frame
       * 2872 shortened info bits are KNOWN ZERO -> LLR +inf, not erasures
       * 12360 punctured parity bits -> LLR 0
    -> LDPC decode                  6.1.3.1 Type A, rate 3/15
    -> 368 Nouter bits              = 200 scrambled L1-Basic bits + 168 BCH parity
    -> BCH syndrome                 6.1.2.1, g(x) = g1(x)...g12(x), degree 168

Usage:
    python m3_l1basic.py hit_rf33.cs16 --rate 8e6
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m3_spec as S                                            # noqa: E402
import m3_freqint as FI                                        # noqa: E402
import m3_ldpc as LD                                           # noqa: E402
import spec_bicm as BI                                         # noqa: E402
from m3_preamble import analyse, QPSK, mer_db                   # noqa: E402

BIG = 60.0                       # LLR magnitude standing for "known bit"


# --------------------------------------------------------------------------
# shortening / puncturing geometry (A/322 6.5.2.4 / 6.5.2.8)
# --------------------------------------------------------------------------

def info_positions(kldpc=3240, nouter=368, pattern=None):
    """Codeword indices carrying the Nouter BCH-encoded bits, ascending.

    A/322 6.5.2.4: Npad whole groups are zero-filled in the order given by the
    shortening pattern, then the FIRST part of the next group in that order is
    additionally padded, then the Nouter bits are mapped sequentially to the
    positions that were NOT padded.
    """
    pattern = pattern or list(BI.SHORTENING_PATTERN["L1-Basic"])
    npad_bits = kldpc - nouter
    npad = npad_bits // 360
    padded = np.zeros(kldpc, bool)
    for j in range(npad):
        g = pattern[j]
        padded[360 * g:360 * (g + 1)] = True
    rem = npad_bits - 360 * npad
    if rem:
        g = pattern[npad]
        padded[360 * g:360 * g + rem] = True
    pos = np.flatnonzero(~padded)
    assert len(pos) == nouter, (len(pos), nouter)
    return pos, padded


def parity_positions(n_keep, kldpc=3240, ninner=16200, groupwise=None,
                     invert=False):
    """Codeword indices of the parity bits that SURVIVE puncturing, in
    transmission order.

    A/322 6.5.2.6 group-wise permutation Y_j = X_{pi_p(j)} over parity groups,
    then 6.5.2.8 punctures the LAST N_punc bits, so the survivors are the FIRST
    n_keep bits of Y.
    """
    gw = groupwise or list(BI.GROUPWISE["L1-Basic"])
    ninfo_g = kldpc // 360
    ngroup = ninner // 360
    # gw is indexed by output group j = ninfo_g .. ngroup-1
    assert len(gw) == ngroup - ninfo_g, (len(gw), ngroup - ninfo_g)
    order = []
    for j in range(ninfo_g, ngroup):
        src = gw[j - ninfo_g]
        if invert:
            # inverse permutation: output group j holds input group gw^-1(j)
            src = gw.index(j) + ninfo_g
        order.extend(range(360 * src, 360 * (src + 1)))
    return np.array(order[:n_keep], dtype=int)


# --------------------------------------------------------------------------
# BCH (A/322 6.1.2.1)
# --------------------------------------------------------------------------

def _gen_poly():
    g = 1
    for p in BI.BCH_16200_GEN_POLYS:
        q = 0
        for e in p:
            q |= 1 << e
        # GF(2) polynomial multiply
        r = 0
        a = g
        b = q
        while b:
            if b & 1:
                r ^= a
            a <<= 1
            b >>= 1
        g = r
    return g


GEN = _gen_poly()
GEN_DEG = GEN.bit_length() - 1


def bch_syndrome(bits):
    """bits[0] is the coefficient of x^(N-1).  Returns the remainder mod g(x)."""
    v = 0
    for b in bits:
        v = (v << 1) | int(b)
    d = len(bits) - 1
    g, gd = GEN, GEN_DEG
    while v.bit_length() - 1 >= gd and v:
        sh = (v.bit_length() - 1) - gd
        v ^= g << sh
    return v


# --------------------------------------------------------------------------
# demap
# --------------------------------------------------------------------------

def cells_to_llr(z, nfec=968, eta=2):
    """A/322 C.1.1 QPSK + 6.5.2.10 block de-interleave -> LLRs, bit 0 positive.

    Annex C.1.1:  00 -> (1+j)/sqrt2, 01 -> (-1+j)/sqrt2,
                  10 -> (1-j)/sqrt2, 11 -> (-1-j)/sqrt2
    so y1 sets the sign of I and y0 sets the sign of Q.
    """
    ncell = nfec // eta
    z = z[:ncell]
    z = z / np.sqrt(np.mean(np.abs(z) ** 2)) * np.sqrt(1.0)
    # noise estimate from the QPSK decision error
    hard = (np.sign(z.real) + 1j * np.sign(z.imag)) / np.sqrt(2)
    n0 = max(np.mean(np.abs(z - hard) ** 2), 1e-6)
    k = 4.0 / n0
    llr_y0 = k * z.imag
    llr_y1 = k * z.real
    out = np.empty(nfec)
    out[:ncell] = llr_y0                     # column 0 of the block interleaver
    out[ncell:] = llr_y1                     # column 1
    return out


def build_frame_llr(bit_llr, nouter=368, kldpc=3240, ninner=16200,
                    groupwise_invert=False):
    n_keep = len(bit_llr) - nouter
    llr = np.zeros(ninner)
    ipos, padded = info_positions(kldpc, nouter)
    llr[:kldpc][padded] = BIG                # shortened bits are known ZERO
    llr[ipos] = bit_llr[:nouter]
    ppos = parity_positions(n_keep, kldpc, ninner, invert=groupwise_invert)
    llr[ppos] = bit_llr[nouter:]
    return llr


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def repetition_bits(mode, nouter=368):
    """A/322 6.5.2.7 Nrepeat for L1-Basic Mode `mode`.

    Nrepeat = 2*floor(C*Nouter) + D, and `spec_bicm.REPETITION` has an entry
    for L1-Basic-1 ONLY -- every other Mode returns 0, which is why one
    multiplex could not reveal this (E82: exactly the shape of the M8 finding,
    one level up).
    """
    r = BI.REPETITION.get("L1-Basic-%d" % mode)
    return 0 if r is None else 2 * int(r["C"] * nouter) + r["D"]


def build_frame_llr_rep(bit_llr, mode, nouter=368, kldpc=3240, ninner=16200,
                        groupwise_invert=False, repeat=True):
    """A/322 6.5.2.9's transmitted word, WITH the 6.5.2.7 parity repetition.

        [ Nouter info+BCH ][ Nrepeat repeated parity ][ Nfec - Nouter tail ]

    and the repeated block is the FIRST Nrepeat bits of the SAME permuted
    parity the tail carries, so their LLRs ADD onto the same codeword
    positions.  This is `m8_l1.decode`'s structure, which M8 derived for
    L1-DETAIL Mode 1 and gated on the air; L1-BASIC Mode 1 is the other half
    of the same clause and was still being read the M4 way.

    THE EVIDENCE THAT THIS IS NOT A GUESS (the same free gate M8 used):
    `S.L1_BASIC_CELLS_PRINTED[1] = 3820` is a number A/322 PRINTS.  Without
    the repetition the arithmetic gives Nfec/eta = 1984 cells; with it,
    (Nfec + Nrepeat)/eta = 3820, exactly.  The receiver was already READING
    3820 cells -- it just mapped 3672 of them onto the wrong codeword
    positions, which is why L1-Basic verified on roughly one RF25 Frame in
    nine while every quality metric read 11-13 dB.

    `repeat=False` reproduces the pre-E82 reading and is kept as the CONTROL
    that must fail (gate_e82 leg 0).  Modes 2..7 have Nrepeat = 0, so for them
    -- including RF33's Mode 3 -- this function is `build_frame_llr` bit for
    bit: `np.add.at` onto a zeroed array at distinct positions IS assignment.
    """
    nrep = repetition_bits(mode, nouter) if repeat else 0
    nlp = ninner - kldpc
    llr = np.zeros(ninner)
    ipos, padded = info_positions(kldpc, nouter)
    llr[:kldpc][padded] = BIG                # shortened bits are known ZERO
    llr[ipos] = bit_llr[:nouter]
    allp = parity_positions(nlp, kldpc, ninner, invert=groupwise_invert)
    o = nouter
    if nrep:
        np.add.at(llr, allp[:nrep], bit_llr[o:o + nrep])
        o += nrep
    tail = bit_llr[o:]
    np.add.at(llr, allp[:len(tail)], tail)
    return llr


def decode_one(z, mode, variant="standard", gw_invert=False, iters=100,
               checks_cache={}, repeat=True):
    ncells = S.L1_BASIC_CELLS_PRINTED[mode]
    nfec, cells, npunc = S.l1_basic_lengths(mode)
    eta = S.L1_BASIC_ETA_MOD[mode]
    bl = cells_to_llr(z[:ncells], nfec=nfec, eta=eta)
    llr = build_frame_llr_rep(bl, mode, groupwise_invert=gw_invert,
                              repeat=repeat)
    key = ("3/15", variant)
    if key not in checks_cache:
        checks_cache[key] = LD.parity_check("3/15", variant=variant)
    checks = checks_cache[key]
    bits, conv, it, bad = LD.min_sum_decode(llr, checks, 16200, iters=iters)
    ipos, _ = info_positions()
    nouter_bits = bits[ipos]
    syn = bch_syndrome(nouter_bits)
    return {"converged": bool(conv), "iters": int(it), "unsatisfied": int(bad),
            "bch_syndrome": int(syn), "bch_ok": syn == 0,
            "nouter_bits": nouter_bits}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("capture")
    ap.add_argument("--rate", type=float, required=True)
    ap.add_argument("--fmt")
    ap.add_argument("--frames", type=int, default=4)
    ap.add_argument("--step", type=float, default=0.7)
    ap.add_argument("--iters", type=int, default=100)
    ap.add_argument("--json")
    a = ap.parse_args()

    path = a.capture
    if not os.path.isabs(path):
        path = os.path.join(HERE, path)
    print(f"A/322 spec tables: "
          f"{sum(1 for _, g, _ in S.verify(verbose=False) if g)}/13 pass")
    print(f"BCH generator degree {GEN_DEG} (A/322 Table 6.3 g1..g12, "
          f"expect 168) -> {'OK' if GEN_DEG == 168 else 'MISMATCH'}")

    zs = []
    for i in range(a.frames):
        try:
            rep = {}
            rep, z, Y, geo = analyse(path, a.rate, a.fmt, report=rep,
                                     start_sec=i * a.step, quiet=True)
        except Exception as e:                                  # noqa: BLE001
            print(f"  frame {i}: skip ({e})")
            continue
        (lo, n, pilot, cp, data, shift, cred, nfft, gi, dx, mode) = geo
        x = FI.deinterleave(z, nfft, 0, direction="forward", toggle="i")
        zs.append((i * a.step, x, mode, rep))
    if not zs:
        raise SystemExit("no frames")
    mode = zs[0][2]
    print(f"\n  L1-Basic Mode {mode}: Ksig 200, Nouter 368, Kldpc 3240, "
          f"Ninner 16200, rate 3/15, QPSK, "
          f"{S.L1_BASIC_CELLS_PRINTED[mode]} cells\n")

    print("  === variant sweep: the air arbitrates the two spec-reading "
          "ambiguities ===")
    print("  (a correct variant converges AND clears BCH; wrong ones are the "
          "controls)\n")
    print("   parity-interleave  groupwise   frame   LDPC conv  iters  unsat"
          "   BCH syndrome")
    results = []
    for variant in LD.PARITY_INTERLEAVE_VARIANTS:
        for gwi in (False, True):
            for (t, x, m, rep) in zs[:2]:
                r = decode_one(x, m, variant=variant, gw_invert=gwi,
                               iters=a.iters)
                results.append({"variant": variant, "gw_invert": gwi,
                                "t": t, **{k: v for k, v in r.items()
                                           if k != "nouter_bits"}})
                print(f"   {variant:16s}  {'inverted' if gwi else 'direct  '}"
                      f"   @{t:4.1f}s   {str(r['converged']):5s}    "
                      f"{r['iters']:4d}  {r['unsatisfied']:5d}   "
                      f"{'ZERO -> PASS' if r['bch_ok'] else hex(r['bch_syndrome'])[:18]}")

    winners = [r for r in results if r["converged"] and r["bch_ok"]]
    print(f"\n  variants that converged AND cleared BCH: {len(winners)} "
          f"of {len(results)}")

    out = {"capture": os.path.basename(path), "sweep": results,
           "n_winners": len(winners)}
    if winners:
        v, gwi = winners[0]["variant"], winners[0]["gw_invert"]
        print(f"  WINNER: parity-interleave '{v}', groupwise "
              f"{'inverted' if gwi else 'direct'}")
        print(f"\n  --- decoding all {len(zs)} frames with the winner ---")
        blocks = []
        for (t, x, m, rep) in zs:
            r = decode_one(x, m, variant=v, gw_invert=gwi, iters=a.iters)
            b = "".join(str(int(q)) for q in r["nouter_bits"])
            blocks.append(b)
            print(f"    @{t:4.1f}s  conv={r['converged']}  "
                  f"BCH={'PASS' if r['bch_ok'] else 'FAIL'}  "
                  f"L1B(200,scrambled)={b[:64]}...")
        same = len(set(blocks)) == 1
        print(f"\n  all {len(blocks)} frames identical: {same}"
              + ("" if same else f"  ({len(set(blocks))} distinct)"))
        out["blocks"] = blocks
        out["identical"] = same
        out["winner"] = {"variant": v, "gw_invert": gwi}
    else:
        print("  NO variant cleared BCH -- reporting the wall rather than a "
              "guess.")

    dest = a.json or os.path.join(
        HERE, "m3_l1basic_" + os.path.splitext(os.path.basename(path))[0] + ".json")
    with open(dest, "w") as fh:
        json.dump(out, fh, indent=2, default=float)
    print(f"\n  wrote {dest}")


if __name__ == "__main__":
    main()
