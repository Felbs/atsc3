#!/usr/bin/env python3
"""A/322 5.2.3 scrambler -- the wall, located precisely.

This is a NEGATIVE result, kept because a precise wall location is worth more
than an optimistic claim.  Everything else in M3 lands; this is the one thing
that does not, and this file says exactly how far it got and why it stopped.

WHAT THE SPEC GIVES CLEANLY
    G(x)          = 1 + X + X^3 + X^6 + X^7 + X^11 + X^12 + X^13 + X^16
    initial load  = 0xF180, reloaded at the start of every block
    Figure 5.6    stage preload X1..X16 = 0 0 0 0 0 0 0 1 1 0 0 0 1 1 1 1
                  (which read X16->X1 is 1111 0001 1000 0000 = 0xF180, so the
                  figure and the stated seed agree)
    procedure     read eight stage outputs D7..D0 as a byte, shift once, repeat
    test vector   first 24 bits = 1100 0000 0110 1101 0011 1111, MSB first

WHAT IS MISSING
    WHICH eight stages are D7..D0.  Figure 5.6 draws them as vector art.

FINDING 1 -- the tap positions ARE recoverable, from glyph columns.
    pdftotext -layout preserves the horizontal positions of the D labels and
    the X labels.  Each D label sits just to the RIGHT of the stage box it taps:

        X stages at columns  0 10 14 22 30 44 48 61 65 76 80 92 101 110 119 124
                             X1 X2 X3 X4 X5 X6 X7 X8 X9 X10 X11 X12 X13 X14 X15 X16
        D labels at columns  3        18  26          57      85  97  106 115
                             D0       D1  D2          D3      D4  D5  D6  D7

    Taking each D from the nearest stage to its LEFT gives

        D0=X1  D1=X3  D2=X4  D3=X7  D4=X11  D5=X12  D6=X13  D7=X14

    and reading those stages out of the Figure 5.6 preload gives
    D7..D0 = 1,1,0,0,0,0,0,0 = 11000000, which is EXACTLY the spec's printed
    first byte.  8 of 8 bits, 1 chance in 256.  Taking each D from the nearest
    stage to its RIGHT instead gives 11101000, which is wrong.  So the tap
    positions are almost certainly right.

FINDING 2 -- but then the printed test vector is unreachable, provably.
    Under a single right-shift, stage X4 of the next state is whatever X3 held,
    and X3 held 0.  So byte 1's D2 bit is FORCED to 0.  The printed byte 1
    (01101101) needs it to be 1.  The same contradiction appears under a single
    left-shift.  No feedback polynomial can rescue this, because feedback only
    writes the END stage -- it cannot touch X4.
    So Figure 5.6 and the printed 24-bit vector are MUTUALLY INCONSISTENT.
    One of them contains an erratum, and the text alone cannot say which.

FINDING 3 -- the air rejects every LFSR hypothesis too.
    Three independent searches, all empty:
      (a) 102960 hypotheses (2 tap sets x 2 directions x C(16,8) subsets x 2
          orderings) against the printed test vector          -> 0 survivors
      (b) 68.7 billion (16 LFSR configurations x 16^8 tap assignments, ALL
          orderings and repeats allowed) against the air's 32-bit CRC, solved
          by meet-in-the-middle rather than enumerated        -> 12 survivors,
          against ~16 expected BY CHANCE, and none of the 12 reproduced the
          independently measured FFT size / GI / pilot pattern / subframe
          count.  Correctly rejected as the coincidences they are.
      (c) constraint propagation from 69 mask bits DERIVED FROM THE AIR --
          47 from L1B_reserved (A/322 3.2.1: "The ATSC default value for
          reserved bits is '1'") plus 22 from fields M2 measured from the
          waveform -- across both the 2024 and 2026 field layouts, reserved
          all-ones and all-zeros, and three readings of the symbol count
                                                              -> 0 survivors

CONCLUSION
    The A/322 L1 scrambler is not a plain 16-stage LFSR with parallel stage
    taps under any convention searched.  Either the inter-byte update is not a
    plain shift, or the outputs are not plain stage taps, or the published
    figure/vector carry an erratum.  Resolving it needs the actual PDF figure
    read by eye, or a known-good A/53 randomizer to compare against -- A/322
    5.2.3 uses the A/53 8-VSB randomizer's polynomial and seed.

    CONSEQUENCE, stated plainly: the 200 L1-Basic bits are recovered and
    PROVEN correct (LDPC converged on 12960 checks, BCH syndrome zero, and 28
    of 28 frame-pair CRC-32 deltas zero), but they remain scrambled, so the
    field VALUES are not yet readable.  Nothing downstream of this is claimed.

    Note that the CRC-32 delta test in m3_crc.py deliberately does not need the
    scrambler at all, which is why Step 3's integrity result stands regardless.
"""
from __future__ import annotations

import json
import os
import sys
from itertools import combinations

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

TEST_VECTOR = "110000000110110100111111"
INIT_FIG = [0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1]     # X1..X16
POLY_EXP = (0, 1, 3, 6, 7, 11, 12, 13, 16)
# Figure 5.6 glyph-column reading, 0-based stage indices
FIGURE_TAPS = {"D0": 0, "D1": 2, "D2": 3, "D3": 6,
               "D4": 10, "D5": 11, "D6": 12, "D7": 13}
TAPSETS = {"spec": tuple(e - 1 for e in POLY_EXP if 1 <= e <= 16),
           "recip": tuple(16 - e for e in POLY_EXP if 1 <= e <= 16)}


def step(s, taps, direction):
    fb = 0
    for t in taps:
        fb ^= s[t]
    return ([fb] + s[:-1]) if direction == "right" else (s[1:] + [fb])


def byte_from(s, order):
    """order[q] = stage supplying bit q of the byte, MSB (D7) first."""
    return "".join(str(s[order[q]]) for q in range(8))


FIG_ORDER = [FIGURE_TAPS[f"D{7-q}"] for q in range(8)]


def finding1():
    print("  FINDING 1 -- tap positions from Figure 5.6 glyph columns")
    print(f"    D0..D7 = X1, X3, X4, X7, X11, X12, X13, X14")
    got = byte_from(INIT_FIG, FIG_ORDER)
    print(f"    byte 0 from the Figure 5.6 preload : {got}")
    print(f"    byte 0 printed in A/322 5.2.3      : {TEST_VECTOR[:8]}")
    print(f"    match: {got == TEST_VECTOR[:8]}  (1 in 256 by chance)")
    right = [FIGURE_TAPS[f"D{7-q}"] + 1 for q in range(8)]
    alt = byte_from(INIT_FIG, [min(15, i + 1) for i in FIG_ORDER])
    print(f"    same read from the stage to the RIGHT instead: {alt} "
          f"-> {'match' if alt == TEST_VECTOR[:8] else 'wrong'}")
    return got == TEST_VECTOR[:8]


def finding2():
    print("\n  FINDING 2 -- byte 1 is then unreachable under any single shift")
    need = TEST_VECTOR[8:16]
    print(f"    printed byte 1: {need}")
    ok_any = False
    for tn, taps in TAPSETS.items():
        for d in ("right", "left"):
            s1 = step(list(INIT_FIG), taps, d)
            got = byte_from(s1, FIG_ORDER)
            print(f"      taps={tn:5s} dir={d:5s} -> {got}  "
                  f"{'MATCH' if got == need else 'no'}")
            ok_any |= got == need
    print("    X4 of the next state is whatever X3 held (=0), so byte 1's D2")
    print("    bit is FORCED to 0; the printed vector needs 1.  Feedback only")
    print("    writes the end stage and cannot touch X4.  Hence Figure 5.6 and")
    print("    the printed vector are mutually inconsistent.")
    return ok_any


def finding3a():
    print("\n  FINDING 3a -- exhaustive subset search vs the printed vector")
    hits = 0
    for taps in TAPSETS.values():
        for d in ("right", "left"):
            for combo in combinations(range(16), 8):
                for order in (list(combo), list(combo)[::-1]):
                    s = list(INIT_FIG)
                    out = ""
                    for _ in range(3):
                        out += byte_from(s, order)
                        s = step(s, taps, d)
                    if out == TEST_VECTOR:
                        hits += 1
    n = 2 * 2 * len(list(combinations(range(16), 8))) * 2
    print(f"    {n} hypotheses searched -> {hits} survivors")
    return hits


def main():
    print(__doc__.split("WHAT THE SPEC GIVES CLEANLY")[0].strip())
    print("\n" + "=" * 74)
    f1 = finding1()
    f2 = finding2()
    f3 = finding3a()
    print("\n" + "=" * 74)
    print(f"  Figure-derived taps reproduce byte 0 : {f1}")
    print(f"  ...and byte 1 under any single shift : {f2}")
    print(f"  Exhaustive subset search survivors   : {f3}")
    print("\n  VERDICT: A/322 5.2.3 scrambler tap geometry UNRESOLVED.")
    print("  The 200 L1-Basic bits are proven correct but stay scrambled.")
    print("  See m3_descramble.py for the air-side searches (3b, 3c).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
