"""
spec_bicm.py -- ATSC A/322 L1-signaling BICM chain constants, extracted from spec text.

PROVENANCE
----------
All tables below were transcribed from pdftotext dumps of the official ATSC A/322
Physical Layer Protocol standard held at lab/spec/ :

    A322_2026_tbl.txt   A/322:2026-04, pdftotext -table   (primary)
    A322_2026.txt       A/322:2026-04, pdftotext -layout  (independent column geometry)
    A322_2024_tbl.txt   A/322:2024-04, pdftotext -table   (independent EDITION)

Every table here was read out of at least two of those three files and required to
agree.  Cross-edition result: 2026-04 and 2024-04 are BYTE-IDENTICAL for every table
in this module (6.3, 6.17, 6.18, 6.19, 6.20, 6.21, 6.22, 6.23, 6.24, C.1.1) and for
the prose of 5.2.3 / 6.5.2.2 / 6.5.2.6 / 6.5.2.10 / 6.5.2.11.  No disagreements found.

KNOWN EXTRACTION HAZARD
-----------------------
The PDF sets all italic math (subscripts, Greek, radicals, ceil/floor brackets) in a
font pdftotext cannot map; those glyphs come out as *nothing at all*.  Parentheses,
digits and roman letters DO survive.  So:
  - table BODIES (all digits) are trustworthy and are permutation-checked below;
  - EQUATIONS lose their subscripts and are reconstructed here, with the evidence
    for each reconstruction spelled out in the relevant docstring;
  - anything that could not be resolved is marked UNRESOLVED and is not guessed.

Page furniture ("ATSC A/322:2026-04", "Physical Layer Protocol", "14 April 2026",
bare page numbers 41..50) is interleaved through the table region in the source text
and was stripped by hand; the permutation self-checks below exist specifically to
catch a swallowed page number.

Run this file to execute all self-checks:  python spec_bicm.py
"""

from fractions import Fraction

__all__ = [
    "SCRAMBLER_POLY_EXPONENTS", "SCRAMBLER_POLY_STR", "SCRAMBLER_INIT",
    "SCRAMBLER_INIT_STAGES", "SCRAMBLER_FIRST_BITS", "SCRAMBLER_PROCEDURE",
    "SCRAMBLER_TAPS_UNRESOLVED",
    "BCH_16200_GEN_POLYS", "BCH_16200_GEN_POLYS_STR", "BCH_64800_GEN_POLYS",
    "MOUTER_16200", "MOUTER_64800",
    "L1_CONFIG", "L1B_KSIG", "L1B_NOUTER", "L1B_KLDPC", "L1B_NINNER",
    "SHORTENING_PATTERN", "SHORTENING_NGROUP",
    "GROUPWISE_L1BASIC", "GROUPWISE_L1D1", "GROUPWISE_L1D2",
    "GROUPWISE_L1D3", "GROUPWISE_L1D4", "GROUPWISE_L1D5",
    "GROUPWISE_L1D6", "GROUPWISE_L1D7",
    "GROUPWISE_FIRST_INDEX", "REPETITION", "PUNCTURING",
    "QPSK_MAP", "QPSK_NORM",
    "l1_basic_frame_params", "groupwise_permute_parity", "demux_qpsk",
]


# ---------------------------------------------------------------------------
# 1.  Section 5.2.3 -- Scrambling of Baseband Packets
#     (Section 6.5.2.2 states the L1 scrambler is THE SAME generator,
#      initialization and operation.)
# ---------------------------------------------------------------------------

# Verbatim from 5.2.3:  "G(x) = 1+X+X3+X6+X7+X11+X12+X13+X16"
SCRAMBLER_POLY_STR = "G(x) = 1 + X + X^3 + X^6 + X^7 + X^11 + X^12 + X^13 + X^16"
SCRAMBLER_POLY_EXPONENTS = (0, 1, 3, 6, 7, 11, 12, 13, 16)

# Verbatim: "The initial sequence (0xF180: 1111 0001 1000 0000) shall be loaded into
#            the shift register at the start of every Baseband Packet."
SCRAMBLER_INIT = 0xF180

# Figure 5.6 draws stages X1..X16 with these preloaded values, left to right.
# 0xF180 fed in LSB-first over X1..X16 reproduces this row exactly, i.e.
#   X[k] = (0xF180 >> (k-1)) & 1
# and the resulting 16-bit word read X16..X1 is 0x018F.
SCRAMBLER_INIT_STAGES = (0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1)

SCRAMBLER_PROCEDURE = """\
1) The initial sequence (0xF180: 1111 0001 1000 0000) shall be loaded into the shift
   register at the start of every Baseband Packet.
2) Eight of the shift register outputs (D7, D6, ..., D0) are used as a randomizing
   byte, which shall then be XOR'd bitwise (MSB to MSB and so on until LSB to LSB)
   with the corresponding byte of Baseband Packet data.
3) The bits in the shift register shall be shifted once.  Go to step 2 above."""

# GOLD -- explicit test vector, verbatim from 5.2.3, identical in the 2024 edition:
#   "The first values of the baseband scrambling sequence are
#    1100 0000 0110 1101 0011 1111 ...
#    (MSB first, or D7, D6, ..., D0, D7, D6, ...)."
SCRAMBLER_FIRST_BITS = "110000000110110100111111"          # 24 bits
SCRAMBLER_FIRST_BYTES = (0xC0, 0x6D, 0x3F)

# --- UNRESOLVED -------------------------------------------------------------
# Figure 5.6 is vector art: the *positions* of the D0..D7 tap points on the
# X1..X16 stage chain are not recoverable from any of the three text dumps, so
# the scrambler is NOT implemented here.  Worse, an exhaustive search proves the
# obvious readings are all wrong -- see self-check `check_scrambler_unresolved`:
#
#   For the published test vector, bit positions D0, D2, D3 and D5 must all carry
#   the same 3-state signature (0,1,1).  Four DISTINCT register stages would have
#   to share that signature.  Starting from the Figure-5.6 state and stepping with
#   G(x) (Fibonacci or Galois, either shift direction, 1 or 8 shifts per byte, and
#   with 0xF180 loaded either end-first), at most TWO stages ever share it.  So the
#   published 24 bits cannot be produced by any fixed-stage byte tap-off of this
#   LFSR under any of those conventions.
#
# Therefore one of the following must hold and must be settled against the actual
# PDF figure (or a known-good A/53 data-randomizer implementation -- this is the
# A/53 8-VSB randomizer polynomial and seed):
#   (a) the D taps are not plain stage outputs, or
#   (b) the update between output bytes is not a plain 1- or 8-step shift, or
#   (c) the printed "first values" line carries an erratum.
# Do NOT ship a guessed scrambler; descramble is a hard gate on L1 decode.
SCRAMBLER_TAPS_UNRESOLVED = True


# ---------------------------------------------------------------------------
# 2.  Section 6.1.2.1 / Table 6.3 -- BCH generator polynomials
#     6.5.2.3: "The systematic BCH code for Ninner = 16200 defined in Section
#     6.1.2.1 shall be used for outer encoding of L1-Basic and L1-Detail."
#     Table 6.18: Mouter = 168 for every L1-Basic and L1-Detail mode.
# ---------------------------------------------------------------------------

MOUTER_16200 = 168          # 12 polynomials x degree 14
MOUTER_64800 = 192          # 12 polynomials x degree 16

# Exponent lists, highest first.  g(x) = g1(x)*g2(x)*...*g12(x), deg g = Mouter.
BCH_16200_GEN_POLYS = (
    (14, 5, 3, 1, 0),                        # g1(x)  = x14+x5+x3+x+1
    (14, 11, 8, 6, 0),                       # g2(x)  = x14+x11+x8+x6+1
    (14, 10, 9, 6, 2, 1, 0),                 # g3(x)  = x14+x10+x9+x6+x2+x+1
    (14, 12, 10, 8, 7, 4, 0),                # g4(x)  = x14+x12+x10+x8+x7+x4+1
    (14, 13, 11, 9, 8, 6, 4, 2, 0),          # g5(x)  = x14+x13+x11+x9+x8+x6+x4+x2+1
    (14, 13, 9, 8, 7, 3, 0),                 # g6(x)  = x14+x13+x9+x8+x7+x3+1
    (14, 13, 11, 10, 7, 6, 5, 2, 0),         # g7(x)  = x14+x13+x11+x10+x7+x6+x5+x2+1
    (14, 11, 10, 9, 8, 5, 0),                # g8(x)  = x14+x11+x10+x9+x8+x5+1
    (14, 10, 9, 3, 2, 1, 0),                 # g9(x)  = x14+x10+x9+x3+x2+x+1
    (14, 12, 11, 9, 6, 3, 0),                # g10(x) = x14+x12+x11+x9+x6+x3+1
    (14, 12, 11, 4, 0),                      # g11(x) = x14+x12+x11+x4+1
    (14, 13, 10, 8, 7, 6, 5, 3, 2, 1, 0),    # g12(x) = x14+x13+x10+x8+x7+x6+x5+x3+x2+x+1
)

BCH_16200_GEN_POLYS_STR = (
    "x14+x5+x3+x+1",
    "x14+x11+x8+x6+1",
    "x14+x10+x9+x6+x2+x+1",
    "x14+x12+x10+x8+x7+x4+1",
    "x14+x13+x11+x9+x8+x6+x4+x2+1",
    "x14+x13+x9+x8+x7+x3+1",
    "x14+x13+x11+x10+x7+x6+x5+x2+1",
    "x14+x11+x10+x9+x8+x5+1",
    "x14+x10+x9+x3+x2+x+1",
    "x14+x12+x11+x9+x6+x3+1",
    "x14+x12+x11+x4+1",
    "x14+x13+x10+x8+x7+x6+x5+x3+x2+x+1",
)

# Retained only as a transcription-quality control (not used for L1).
BCH_64800_GEN_POLYS = (
    (16, 5, 3, 2, 0),
    (16, 8, 6, 5, 4, 1, 0),
    (16, 11, 10, 9, 8, 7, 5, 4, 3, 2, 0),
    (16, 14, 12, 11, 9, 6, 4, 2, 0),
    (16, 12, 11, 10, 9, 8, 5, 3, 2, 1, 0),
    (16, 15, 14, 13, 12, 10, 9, 8, 7, 5, 4, 2, 0),
    (16, 15, 13, 11, 10, 9, 8, 6, 5, 2, 0),
    (16, 14, 13, 12, 9, 8, 6, 5, 2, 1, 0),
    (16, 11, 10, 9, 7, 5, 0),
    (16, 14, 13, 12, 10, 8, 7, 5, 2, 1, 0),
    (16, 13, 12, 11, 9, 5, 3, 2, 0),
    (16, 12, 11, 9, 7, 6, 5, 1, 0),
)

# Encoding rule, 6.1.2.1 (systematic, MSB-first):
#   s(x) = m(x)*x^Mouter - p(x),  p(x) = m(x)*x^Mouter mod g(x)
BCH_ENCODE_RULE = "s(x) = m(x)*x^Mouter + (m(x)*x^Mouter mod g(x)), g = g1*g2*...*g12"


# ---------------------------------------------------------------------------
# 3.  Table 6.17 / 6.18 / 6.19 -- L1 configuration
# ---------------------------------------------------------------------------

L1B_KSIG = 200          # Table 6.17/6.18, fixed for every L1-Basic mode
L1B_NOUTER = 368        # = Ksig + Mouter = 200 + 168
L1B_KLDPC = 3240        # Table 6.19
L1B_NINNER = 16200
L1B_NINFO_GROUP = 9     # = Kldpc/360, Table 6.19
L1B_CODE_RATE = Fraction(3, 15)   # Table 6.17, LDPC Type A

# Table 6.17: (constellation, cells) per L1-Basic mode.  All L1-Basic modes use the
# 16200-length, rate-3/15, Type-A LDPC.
L1_CONFIG = {
    1: ("QPSK",         3820),
    2: ("QPSK",          934),
    3: ("QPSK",          484),
    4: ("NUC_16_8/15",   259),
    5: ("NUC_64_9/15",   163),
    6: ("NUC_256_9/15",  112),
    7: ("NUC_256_13/15",  69),
}

# Table 6.19 -- L1-Detail Kldpc / Ninfo_group
L1D_KLDPC = {1: 3240, 2: 3240, 3: 6480, 4: 6480, 5: 6480, 6: 6480, 7: 6480}
L1D_NINFO_GROUP = {m: L1D_KLDPC[m] // 360 for m in range(1, 8)}

# Table 6.25 -- Kseg for L1-Detail segmentation
L1D_KSEG = {1: 2352, 2: 3072, 3: 6312, 4: 6312, 5: 6312, 6: 6312, 7: 6312}


# ---------------------------------------------------------------------------
# 4.  Table 6.20 -- Shortening Pattern of Information Bit Group to be Padded
#
#     6.5.2.4:  Npad = floor((Kldpc - Nouter)/360); the groups
#     pi_s(0)..pi_s(Npad-1) are fully zero-padded, then group pi_s(Npad) has its
#     first (Kldpc - Nouter - 360*Npad) bits padded.
#
#     CONFIRMED: the L1-Basic row is  4 1 5 2 8 6 0 7 3  (9 entries, Ngroup = 9)
#     -- exactly as read by the caller.  Ngroup is a merged/spanned cell in the
#     PDF (rendered on the "L1-Detail Mode 1" and "L1-Detail Mode 5" rows in
#     -table mode); the -layout dump prints the two values 9 and 18 separately in
#     their own column, and they agree with Ninfo_group in Table 6.19.
# ---------------------------------------------------------------------------

SHORTENING_NGROUP = {
    "L1-Basic":       9,
    "L1-Detail-1":    9,
    "L1-Detail-2":    9,
    "L1-Detail-3":   18,
    "L1-Detail-4":   18,
    "L1-Detail-5":   18,
    "L1-Detail-6":   18,
    "L1-Detail-7":   18,
}

SHORTENING_PATTERN = {
    # L1-Basic, all 7 modes
    "L1-Basic":    (4, 1, 5, 2, 8, 6, 0, 7, 3),
    "L1-Detail-1": (7, 8, 5, 4, 1, 2, 6, 3, 0),
    "L1-Detail-2": (6, 1, 7, 8, 0, 2, 4, 3, 5),
    "L1-Detail-3": (0, 12, 15, 13, 2, 5, 7, 9, 8,
                    6, 16, 10, 14, 1, 17, 11, 4, 3),
    "L1-Detail-4": (0, 15, 5, 16, 17, 1, 6, 13, 11,
                    4, 7, 12, 8, 14, 2, 3, 9, 10),
    "L1-Detail-5": (2, 4, 5, 17, 9, 7, 1, 6, 15,
                    8, 10, 14, 16, 0, 11, 13, 12, 3),
    "L1-Detail-6": (0, 15, 5, 16, 17, 1, 6, 13, 11,
                    4, 7, 12, 8, 14, 2, 3, 9, 10),   # identical to Mode 4
    "L1-Detail-7": (15, 7, 8, 11, 5, 10, 16, 4, 12,
                    3, 0, 6, 9, 1, 14, 17, 2, 13),
}


# ---------------------------------------------------------------------------
# 5.  Section 6.5.2.6 -- Parity Permutation
# ---------------------------------------------------------------------------

PARITY_INTERLEAVER_RULE = """\
6.5.2.6 parity interleaver (u = interleaved LDPC codeword):
    u_i = c_i                                for 0 <= i < Kldpc   (info not interleaved)
    u_{Kldpc + 360*t + s} = c_{Kldpc + 27*s + t}   for 0 <= s < 360, 0 <= t < 27
Used ONLY for L1-Detail Modes 3,4,5,6,7 (Kldpc=6480, Nldpc_parity=9720 = 27*360).
For L1-Basic and L1-Detail Modes 1 and 2 the parity interleaver is NOT used
(u_i = c_i for all i) because it is already folded into the LDPC encoding.
=> for L1-Basic Mode 3 this stage is a no-op."""

# Group-wise parity permutation.  u is split into Ngroup = Ninner/360 groups X_j.
#
#   Y_j = X_j                       for 0 <= j < Kldpc/360      (info groups fixed)
#   Y_j = X_{pi_p(j)}               for Kldpc/360 <= j < Ngroup (parity groups)
#
# DIRECTION EVIDENCE (the subscripts themselves are lost to the math font):
#   The parentheses of "pi_p(j)" DO survive extraction.  In Section 6.2.2, whose
#   rule is unambiguous ("Y_j = X_{pi(j)}"), the line renders as "= ()" -- parens
#   AFTER the equals sign.  The second line of 6.5.2.6 renders identically,
#   "= (), /360   < ", i.e. the parens are on the RIGHT of "=" and the left side
#   carries no parens.  Had the standard written "Y_{pi_p(j)} = X_j" the parens
#   would render to the LEFT of "=".  Both the 2026 -table and -layout dumps and
#   the 2024 dump show parens on the right.
#   Confidence: high, but not certified.  `groupwise_permute_parity(..., invert=)`
#   is provided so a decoder can try the other direction cheaply if L1 refuses to
#   decode with everything else verified.
GROUPWISE_DIRECTION = "Y_j = X_{pi_p(j)}"

# Table 6.21 -- Group-wise Interleaving Pattern for all L1-Basic Modes and
#               L1-Detail Modes 1 and 2.   Ngroup = 45; entries are pi_p(9)..pi_p(44)
#               (36 values, a permutation of {9..44}).
GROUPWISE_FIRST_INDEX = {"L1-Basic": 9, "L1-Detail-1": 9, "L1-Detail-2": 9,
                         "L1-Detail-3": 18, "L1-Detail-4": 18, "L1-Detail-5": 18,
                         "L1-Detail-6": 18, "L1-Detail-7": 18}

GROUPWISE_L1BASIC = (
    20, 23, 25, 32, 38, 41, 18,  9, 10, 11, 31, 24,
    14, 15, 26, 40, 33, 19, 28, 34, 16, 39, 27, 30,
    21, 44, 43, 35, 42, 36, 12, 13, 29, 22, 37, 17,
)

GROUPWISE_L1D1 = (
    16, 22, 27, 30, 37, 44, 20, 23, 25, 32, 38, 41,
     9, 10, 17, 18, 21, 33, 35, 14, 28, 12, 15, 19,
    11, 24, 29, 34, 36, 13, 40, 43, 31, 26, 39, 42,
)

GROUPWISE_L1D2 = (
     9, 31, 23, 10, 11, 25, 43, 29, 36, 16, 27, 34,
    26, 18, 37, 15, 13, 17, 35, 21, 20, 24, 44, 12,
    22, 40, 19, 32, 38, 41, 30, 33, 14, 28, 39, 42,
)

# Table 6.22 -- Group-wise Interleaving Pattern for L1-Detail Modes 3,4,5,6,7.
#               Ngroup = 45; entries are pi_p(18)..pi_p(44)
#               (27 values, a permutation of {18..44}).
GROUPWISE_L1D3 = (
    19, 37, 30, 42, 23, 44, 27, 40, 21, 34, 25, 32, 29, 24,
    26, 35, 39, 20, 18, 43, 31, 36, 38, 22, 33, 28, 41,
)

GROUPWISE_L1D4 = (
    20, 35, 42, 39, 26, 23, 30, 18, 28, 37, 32, 27, 44, 43,
    41, 40, 38, 36, 34, 33, 31, 29, 25, 24, 22, 21, 19,
)

GROUPWISE_L1D5 = (
    19, 37, 33, 26, 40, 43, 22, 29, 24, 35, 44, 31, 27, 20,
    21, 39, 25, 42, 34, 18, 32, 38, 23, 30, 28, 36, 41,
)

GROUPWISE_L1D6 = GROUPWISE_L1D4          # spec rows are identical

GROUPWISE_L1D7 = (
    44, 23, 29, 33, 24, 28, 21, 27, 42, 18, 22, 31, 32, 37,
    43, 30, 25, 35, 20, 34, 39, 36, 19, 41, 40, 26, 38,
)

GROUPWISE = {
    "L1-Basic":    GROUPWISE_L1BASIC,
    "L1-Detail-1": GROUPWISE_L1D1,
    "L1-Detail-2": GROUPWISE_L1D2,
    "L1-Detail-3": GROUPWISE_L1D3,
    "L1-Detail-4": GROUPWISE_L1D4,
    "L1-Detail-5": GROUPWISE_L1D5,
    "L1-Detail-6": GROUPWISE_L1D6,
    "L1-Detail-7": GROUPWISE_L1D7,
}


# ---------------------------------------------------------------------------
# 6.  Table 6.23 (Repetition) and Table 6.24 (Puncturing)
#
#     6.5.2.7 Step 1:  Nrepeat = 2 * floor(C * Nouter) + D
#         (repetition applies ONLY to L1-Basic Mode 1 and L1-Detail Mode 1)
#     6.5.2.8 Step 1:  Npunc_temp = floor(A * (Kldpc - Nouter)) + B
#             Step 2:  Nfec_temp  = Nouter + Nldpc_parity - Npunc_temp
#             Step 3:  Nfec       = floor(Nfec_temp / eta_MOD) * eta_MOD
#             Step 4:  Npunc      = Npunc_temp - (Nfec - Nfec_temp)
#     Transmitted bit count = Nfec + Nrepeat;  cells = (Nfec + Nrepeat)/eta_MOD.
# ---------------------------------------------------------------------------

# Table 6.23: (Nouter, Ksig, Kldpc, C, D, Nldpc_parity, eta_MOD)
REPETITION = {
    "L1-Basic-1":  dict(Nouter=368, Ksig=200, Kldpc=3240,
                        C=Fraction(0), D=+3672, Nldpc_parity=12960, eta=2),
    "L1-Detail-1": dict(Nouter=(368, 2520), Ksig=(200, 2352), Kldpc=3240,
                        C=Fraction(61, 16), D=-508, Nldpc_parity=12960, eta=2),
}

# Table 6.24: (Nouter, Kldpc, A, B, Nldpc_parity, eta_MOD)
PUNCTURING = {
    "L1-Basic-1": dict(Nouter=368, Kldpc=3240, A=Fraction(0), B=9360,  Nldpc_parity=12960, eta=2),
    "L1-Basic-2": dict(Nouter=368, Kldpc=3240, A=Fraction(0), B=11460, Nldpc_parity=12960, eta=2),
    "L1-Basic-3": dict(Nouter=368, Kldpc=3240, A=Fraction(0), B=12360, Nldpc_parity=12960, eta=2),
    "L1-Basic-4": dict(Nouter=368, Kldpc=3240, A=Fraction(0), B=12292, Nldpc_parity=12960, eta=4),
    "L1-Basic-5": dict(Nouter=368, Kldpc=3240, A=Fraction(0), B=12350, Nldpc_parity=12960, eta=6),
    "L1-Basic-6": dict(Nouter=368, Kldpc=3240, A=Fraction(0), B=12432, Nldpc_parity=12960, eta=8),
    "L1-Basic-7": dict(Nouter=368, Kldpc=3240, A=Fraction(0), B=12776, Nldpc_parity=12960, eta=8),

    "L1-Detail-1": dict(Nouter=(368, 2520), Kldpc=3240, A=Fraction(7, 2),    B=0,
                        Nldpc_parity=12960, eta=2),
    "L1-Detail-2": dict(Nouter=(368, 3240), Kldpc=3240, A=Fraction(2),       B=6036,
                        Nldpc_parity=12960, eta=2),
    "L1-Detail-3": dict(Nouter=(368, 6480), Kldpc=6480, A=Fraction(11, 16),  B=4653,
                        Nldpc_parity=9720, eta=2),
    "L1-Detail-4": dict(Nouter=(368, 6480), Kldpc=6480, A=Fraction(29, 32),  B=3200,
                        Nldpc_parity=9720, eta=4),
    "L1-Detail-5": dict(Nouter=(368, 6480), Kldpc=6480, A=Fraction(3, 4),    B=4284,
                        Nldpc_parity=9720, eta=6),
    "L1-Detail-6": dict(Nouter=(368, 6480), Kldpc=6480, A=Fraction(11, 16),  B=4900,
                        Nldpc_parity=9720, eta=8),
    "L1-Detail-7": dict(Nouter=(368, 6480), Kldpc=6480, A=Fraction(49, 256), B=8246,
                        Nldpc_parity=9720, eta=8),
}


def l1_basic_frame_params(mode):
    """Return the full L1-Basic FEC-frame geometry for `mode` (1..7).

    Implements 6.5.2.7 Step 1 and 6.5.2.8 Steps 1-4 verbatim.
    """
    p = PUNCTURING["L1-Basic-%d" % mode]
    Nouter, Kldpc, A, B = p["Nouter"], p["Kldpc"], p["A"], p["B"]
    Nldpc_parity, eta = p["Nldpc_parity"], p["eta"]

    r = REPETITION.get("L1-Basic-%d" % mode)
    if r is None:
        Nrepeat = 0
    else:
        Nrepeat = 2 * int(r["C"] * r["Nouter"]) + r["D"]

    Npunc_temp = int(A * (Kldpc - Nouter)) + B
    Nfec_temp = Nouter + Nldpc_parity - Npunc_temp
    Nfec = (Nfec_temp // eta) * eta
    Npunc = Npunc_temp - (Nfec - Nfec_temp)
    Ntx = Nfec + Nrepeat
    return dict(mode=mode, Ksig=L1B_KSIG, Nouter=Nouter, Kldpc=Kldpc,
                Ninner=L1B_NINNER, Nldpc_parity=Nldpc_parity, eta=eta,
                Nrepeat=Nrepeat, Npunc_temp=Npunc_temp, Nfec_temp=Nfec_temp,
                Nfec=Nfec, Npunc=Npunc, Ntx_bits=Ntx, cells=Ntx // eta)


def groupwise_permute_parity(groups, kind="L1-Basic", invert=False):
    """Apply the 6.5.2.6 group-wise parity permutation to a list of 360-bit groups.

    `groups` is the length-45 list X_0..X_44 (each element anything -- typically a
    360-long sequence).  Info groups (j < Kldpc/360) pass through unchanged.
    Default direction is Y_j = X_{pi_p(j)}; `invert=True` gives Y_{pi_p(j)} = X_j.
    """
    pat = GROUPWISE[kind]
    first = GROUPWISE_FIRST_INDEX[kind]
    out = list(groups)
    for j, pj in enumerate(pat, start=first):
        if invert:
            out[pj] = groups[j]
        else:
            out[j] = groups[pj]
    return out


# ---------------------------------------------------------------------------
# 7.  Section 6.5.2.10 -- Bit Demuxing
# ---------------------------------------------------------------------------

BIT_DEMUX_RULE = """\
6.5.2.10, verbatim in both editions:

  "Following zero padded bit removal, the remaining bits of length NFEC or
   (NFEC + Nrepeat) shall be written serially into the Block Interleaver
   column-wise, where the number of columns shall be the same as the modulation
   order.  In the read operation, the bits for one constellation symbol shall be
   read out sequentially row-wise and fed into the bit demultiplexer block."

  "Depending on modulation order, there are two mapping rules.  In the case of
   QPSK, the reliability of bits in a symbol is equal.  Therefore, a bit group
   read out from the Block Interleaver shall be mapped directly to a QAM symbol
   without any intervening operation.  In the cases of higher order modulations
   a bit group shall be mapped to a QAM symbol with the rule described as
   follows:
       Sdemux_in (i)  = {bi(0), bi(1), ..., bi(MOD-1)}
       Sdemux_out(i)  = {ci(0), ci(1), ..., ci(MOD-1)}
       ci(0) = bi(i % MOD), ci(1) = bi((i+1) % MOD), ...,
       ci(MOD-1) = bi((i + MOD - 1) % MOD)"

So for QPSK the DEMUX ITSELF IS THE IDENTITY -- but the BLOCK INTERLEAVER in
front of it is NOT.  There is no table of per-modulation demux patterns; the
cyclic-shift formula above is the whole of it, and it is skipped for QPSK.

Block interleaver, MOD = eta columns, R = Ntx/eta rows, column-wise write /
row-wise read:  cell s takes bits ( x[s], x[R + s], x[2R + s], ..., x[(eta-1)R + s] )
with x[0] the first bit after zero removal.  For QPSK (eta=2, R = Ntx/2):
    y0 = x[s]           (MSB of the cell word)
    y1 = x[R + s]
"""


def demux_qpsk(bits):
    """Block-interleave + (identity) demux for QPSK.  Returns list of (y0,y1)."""
    n = len(bits)
    assert n % 2 == 0, "QPSK needs an even bit count"
    r = n // 2
    return [(bits[s], bits[r + s]) for s in range(r)]


# ---------------------------------------------------------------------------
# 8.  Section 6.5.2.11 / 6.3.4.1 / Table C.1.1 -- Constellation mapping (QPSK)
#
#     6.5.2.11: "Each demultiplexed LDPC block shall be mapped onto constellation
#     symbols ... using constellations as described in Section 6.3."
#     6.3.4.1:  "A two-bit input (y0,s, y1,s) shall be mapped to data cells
#     following the mapping in Table C.1.1", first column = (y0,s, y1,s),
#     y0,s is the MSB.
#
#     Table C.1.1 extracts as:
#         00 -> (1 + j1)/2        01 -> (-1 + j1)/2
#         10 -> (+1 - j1)/2       11 -> (-1 - j1)/2
#     The radical glyph does not survive pdftotext -- the divisor is sqrt(2), not
#     2.  Both editions and both extraction modes lose it identically, so this is
#     glyph loss, not an edition difference.  Cross-check: every Annex-C NUC
#     position vector is unit-average-power (|w| ~ 1, e.g. 0.7062+j0.7075 for
#     16QAM CR 2/15), and only /sqrt(2) makes QPSK unit average power too.
#     This is Gray labelling: y0 sets the sign of I, y1 sets the sign of Q,
#     bit 0 -> "+", bit 1 -> "-".
# ---------------------------------------------------------------------------

QPSK_NORM = 2 ** -0.5

QPSK_MAP = {
    (0, 0): complex(+1, +1) * QPSK_NORM,
    (0, 1): complex(-1, +1) * QPSK_NORM,
    (1, 0): complex(+1, -1) * QPSK_NORM,
    (1, 1): complex(-1, -1) * QPSK_NORM,
}

QPSK_TABLE_VERBATIM = (
    ("00", "(1 + j1)/sqrt(2)"),
    ("01", "(-1 + j1)/sqrt(2)"),
    ("10", "(+1 - j1)/sqrt(2)"),
    ("11", "(-1 - j1)/sqrt(2)"),
)


def qpsk_demap_llr(z, n0=1.0):
    """Soft demap: returns (llr_y0, llr_y1), positive => bit 0 more likely."""
    return (4.0 * QPSK_NORM * z.real / n0, 4.0 * QPSK_NORM * z.imag / n0)


# ===========================================================================
#                                SELF-CHECKS
# ===========================================================================

def _gf2_mul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        a <<= 1
        b >>= 1
    return r


def _gf2_deg(a):
    return a.bit_length() - 1


def _gf2_mod(a, m):
    dm = _gf2_deg(m)
    while _gf2_deg(a) >= dm and a:
        a ^= m << (_gf2_deg(a) - dm)
    return a


def _gf2_gcd(a, b):
    while b:
        a, b = b, _gf2_mod(a, b)
    return a


def _poly_int(exps):
    v = 0
    for e in exps:
        v |= 1 << e
    return v


def _gf2_powmod(base, e, m):
    r = 1
    base = _gf2_mod(base, m)
    while e:
        if e & 1:
            r = _gf2_mod(_gf2_mul(r, base), m)
        base = _gf2_mod(_gf2_mul(base, base), m)
        e >>= 1
    return r


def _is_irreducible(p):
    n = _gf2_deg(p)
    # x^(2^n) == x mod p
    if _gf2_powmod(2, 1 << n, p) != _gf2_mod(2, p):
        return False
    for q in {2, 7} if n == 14 else {2}:
        if n % q:
            continue
        t = _gf2_powmod(2, 1 << (n // q), p) ^ 2
        if _gf2_deg(_gf2_gcd(p, t)) != 0:
            return False
    return True


def check_permutations():
    print("--- (b) permutation checks ---")
    ok = True
    for name, pat in sorted(SHORTENING_PATTERN.items()):
        n = SHORTENING_NGROUP[name]
        good = (len(pat) == n and sorted(pat) == list(range(n)))
        ok &= good
        print("  Table 6.20  %-13s N=%2d  len=%2d  %s"
              % (name, n, len(pat), "PERMUTATION of range(%d) OK" % n if good else "*** FAIL ***"))
    for name, pat in sorted(GROUPWISE.items()):
        first = GROUPWISE_FIRST_INDEX[name]
        exp = list(range(first, 45))
        good = (len(pat) == len(exp) and sorted(pat) == exp)
        ok &= good
        tbl = "6.21" if first == 9 else "6.22"
        print("  Table %s  %-13s N=%2d  len=%2d  %s"
              % (tbl, name, len(exp), len(pat),
                 "PERMUTATION of {%d..44} OK" % first if good else "*** FAIL ***"))
    assert ok
    return ok


def check_bch():
    print("--- (d) BCH polynomial checks ---")
    ok = True
    assert len(BCH_16200_GEN_POLYS) == 12
    for i, exps in enumerate(BCH_16200_GEN_POLYS, 1):
        d = max(exps)
        good = (d == 14 and min(exps) == 0)
        ok &= good
        print("  g%-2d deg=%2d  %-36s %s" % (i, d, BCH_16200_GEN_POLYS_STR[i - 1],
                                             "" if good else "*** FAIL ***"))
    total = sum(max(e) for e in BCH_16200_GEN_POLYS)
    print("  sum(deg g1..g12) = %d   Mouter (Table 6.18) = %d   %s"
          % (total, MOUTER_16200, "OK" if total == MOUTER_16200 else "*** FAIL ***"))
    ok &= (total == MOUTER_16200)

    ints = [_poly_int(e) for e in BCH_16200_GEN_POLYS]
    irr = [_is_irreducible(p) for p in ints]
    print("  all 12 irreducible over GF(2): %s" % ("YES" if all(irr) else "*** NO ***"))
    ok &= all(irr)
    print("  all 12 distinct: %s" % ("YES" if len(set(ints)) == 12 else "*** NO ***"))
    ok &= (len(set(ints)) == 12)

    # DECISIVE independent check: g1 is the minimal polynomial of a primitive
    # element alpha of GF(2^14).  Build GF(2^14) with g1 as the field polynomial
    # and recompute the minimal polynomials of alpha^1, alpha^3, ..., alpha^23
    # (a 12-bit-error-correcting narrow-sense BCH code needs roots alpha^1..alpha^24).
    # If ANY exponent in the table were mis-transcribed this reconstruction fails.
    m = ints[0]
    order = (1 << 14) - 1
    assert _gf2_powmod(2, order, m) == 1
    # primitivity of g1
    prim = all(_gf2_powmod(2, order // q, m) != 1 for q in (3, 43, 127))  # 16383 = 3*43*127
    print("  g1 is primitive (alpha = x has order 2^14-1 = 16383): %s"
          % ("YES" if prim else "*** NO ***"))
    ok &= prim

    def minpoly(i):
        # cyclotomic coset of i
        cs, k = [], i % order
        while k not in cs:
            cs.append(k)
            k = (k * 2) % order
        # prod (X - alpha^c) as a poly in X with GF(2^14) coefficients
        coeffs = [1]                      # coeffs[k] multiplies X^k, field elements
        for c in cs:
            a = _gf2_powmod(2, c, m)
            new = [0] * (len(coeffs) + 1)
            for k2, cf in enumerate(coeffs):
                new[k2 + 1] ^= cf                       # X * term
                new[k2] ^= _gf2_mod(_gf2_mul(cf, a), m)  # alpha^c * term
            coeffs = new
        assert all(c in (0, 1) for c in coeffs), "minpoly not over GF(2)"
        return _poly_int([k for k, c in enumerate(coeffs) if c])

    want = [minpoly(i) for i in range(1, 24, 2)]
    match = (want == ints)
    print("  minimal polys of alpha^1,alpha^3,...,alpha^23 reconstructed from"
          " GF(2^14)=GF(2)[x]/g1")
    print("    reconstructed set == Table 6.3 (16200) column, IN ORDER: %s"
          % ("YES" if match else "set-equal: %s" % (sorted(want) == sorted(ints))))
    ok &= match
    print("    => designed distance covers alpha^1..alpha^24 => t = 12 errors: %s"
          % ("CONFIRMED" if match else "*** NOT CONFIRMED ***"))

    # 64800 column, degree-only control
    d64 = [max(e) for e in BCH_64800_GEN_POLYS]
    print("  (control) 64800 column: 12 polys, degrees %s, sum=%d (Mouter=%d) %s"
          % (set(d64), sum(d64), MOUTER_64800,
             "OK" if sum(d64) == MOUTER_64800 else "*** FAIL ***"))
    ok &= (sum(d64) == MOUTER_64800)
    assert ok
    return ok


def check_repetition_puncturing():
    print("--- Table 6.23 / 6.24 arithmetic vs Table 6.17 cell counts ---")
    ok = True
    print("  mode  eta  A  B      Nrepeat  Npunc_temp  Nfec  Npunc   cells  Table6.17")
    for mode in range(1, 8):
        p = l1_basic_frame_params(mode)
        want = L1_CONFIG[mode][1]
        good = (p["cells"] == want)
        ok &= good
        print("   %d     %d   0  %-6d %6d   %8d  %5d %6d  %5d  %5d %s"
              % (mode, p["eta"], PUNCTURING["L1-Basic-%d" % mode]["B"],
                 p["Nrepeat"], p["Npunc_temp"], p["Nfec"], p["Npunc"],
                 p["cells"], want, "OK" if good else "*** FAIL ***"))
    print("  All 7 L1-Basic cell counts reproduced from Table 6.24 -> %s"
          % ("PASS" if ok else "FAIL"))
    assert ok
    return ok


def check_scrambler_unresolved():
    """Prove the published 24-bit vector is not reachable with the stated LFSR."""
    print("--- (1) scrambler: test-vector reachability ---")
    import itertools
    poly = [1, 3, 6, 7, 11, 12, 13, 16]
    inits = {
        "LSB-first load (matches Fig 5.6 row)":
            {k: (SCRAMBLER_INIT >> (k - 1)) & 1 for k in range(1, 17)},
        "MSB-first load":
            {k: (SCRAMBLER_INIT >> (16 - k)) & 1 for k in range(1, 17)},
    }
    assert tuple(inits["LSB-first load (matches Fig 5.6 row)"][k]
                 for k in range(1, 17)) == SCRAMBLER_INIT_STAGES

    def step(x, ft, sd):
        x = dict(x)
        fb = 0
        for t in ft:
            fb ^= x[t]
        if sd == "fib-up":
            for k in range(16, 1, -1):
                x[k] = x[k - 1]
            x[1] = fb
        elif sd == "fib-down":
            for k in range(1, 16):
                x[k] = x[k + 1]
            x[16] = fb
        elif sd == "gal-up":
            f = x[16]
            y = {k: x[k - 1] ^ (f if k in ft else 0) for k in range(16, 1, -1)}
            y[1] = f
            x = y
        else:
            f = x[1]
            y = {k: x[k + 1] ^ (f if k in ft else 0) for k in range(1, 16)}
            y[16] = f
            x = y
        return x

    need = [tuple((b >> i) & 1 for b in SCRAMBLER_FIRST_BYTES) for i in range(8)]
    best = 0
    found = []
    for iname, s0 in inits.items():
        for fname, ft in (("G(x)", poly), ("G(x) mirrored", sorted(17 - e for e in poly))):
            for sd in ("fib-up", "fib-down", "gal-up", "gal-down"):
                for nsh in (1, 8):
                    s1 = s0
                    for _ in range(nsh):
                        s1 = step(s1, ft, sd)
                    s2 = s1
                    for _ in range(nsh):
                        s2 = step(s2, ft, sd)
                    sig = {k: (s0[k], s1[k], s2[k]) for k in range(1, 17)}
                    cand = [[k for k in range(1, 17) if sig[k] == need[i]] for i in range(8)]
                    # D0,D2,D3,D5 all need signature (0,1,1): 4 distinct stages
                    best = max(best, len(cand[0]))
                    if any(not c for c in cand):
                        continue
                    # NOTE: no ordering constraint imposed -- any 8 DISTINCT stages
                    # in any order would count as a solution.  There are none.
                    for combo in itertools.product(*cand):
                        if len(set(combo)) == 8:
                            found.append((iname, fname, sd, nsh, combo))
                            break
    print("  first bits (spec, verbatim): %s" % SCRAMBLER_FIRST_BITS)
    print("  D0/D2/D3/D5 all require the same 3-state signature (0,1,1);")
    print("  max number of register stages sharing it over all 32 conventions: %d (need 4)" % best)
    print("  consistent tap assignments (8 distinct stages, ANY order): %d" % len(found))
    print("  => scrambler tap geometry UNRESOLVED (see SCRAMBLER_TAPS_UNRESOLVED)")
    assert len(found) == 0, "unexpected: a consistent assignment exists -> revisit"
    assert SCRAMBLER_TAPS_UNRESOLVED
    return True


def check_constellation():
    print("--- (6) QPSK constellation ---")
    pwr = sum(abs(z) ** 2 for z in QPSK_MAP.values()) / 4.0
    print("  mean |z|^2 with 1/sqrt(2) = %.6f (must be 1.0)" % pwr)
    assert abs(pwr - 1.0) < 1e-12
    for b, s in QPSK_TABLE_VERBATIM:
        y = (int(b[0]), int(b[1]))
        print("   y0y1=%s -> %-20s  = %+.4f%+.4fj" % (b, s, QPSK_MAP[y].real, QPSK_MAP[y].imag))
    # Gray: adjacent labels differ in one bit and are nearest neighbours
    assert QPSK_MAP[(0, 0)].real > 0 and QPSK_MAP[(0, 0)].imag > 0
    assert QPSK_MAP[(1, 1)].real < 0 and QPSK_MAP[(1, 1)].imag < 0
    return True


def check_mode3_chain():
    print("--- L1-Basic Mode 3 (the mode being decoded) ---")
    p = l1_basic_frame_params(3)
    for k in ("Ksig", "Nouter", "Kldpc", "Ninner", "Nldpc_parity", "eta",
              "Nrepeat", "Npunc_temp", "Nfec_temp", "Nfec", "Npunc", "Ntx_bits", "cells"):
        print("   %-13s = %s" % (k, p[k]))
    assert (p["Ksig"], p["Nouter"], p["Kldpc"], p["Ninner"]) == (200, 368, 3240, 16200)
    assert p["Nrepeat"] == 0 and p["cells"] == 484 and p["Nfec"] == 968
    Npad = (L1B_KLDPC - L1B_NOUTER) // 360
    rem = (L1B_KLDPC - L1B_NOUTER) - 360 * Npad
    print("   zero padding: Npad=%d full groups %s, then %d bits of group %d"
          % (Npad, list(SHORTENING_PATTERN["L1-Basic"][:Npad]), rem,
             SHORTENING_PATTERN["L1-Basic"][Npad]))
    assert Npad == 7 and rem == 352
    print("   parity interleaver: NOT USED for L1-Basic (6.5.2.6) -> no-op")
    print("   group-wise parity permutation: %s over groups 9..44" % GROUPWISE_DIRECTION)
    print("   demux: QPSK -> identity; block interleaver 2 cols x %d rows"
          % (p["Ntx_bits"] // 2))
    return True


if __name__ == "__main__":
    print("=" * 78)
    print("spec_bicm.py self-check -- ATSC A/322 L1 BICM tables")
    print("=" * 78)
    check_permutations()
    print()
    check_bch()
    print()
    check_repetition_puncturing()
    print()
    check_scrambler_unresolved()
    print()
    check_constellation()
    print()
    check_mode3_chain()
    print()
    print("ALL SELF-CHECKS PASSED "
          "(scrambler tap geometry remains UNRESOLVED by design -- see module docstring)")
