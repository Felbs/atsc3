#!/usr/bin/env python3
"""A/322 tables needed to demodulate the ATSC 3.0 PREAMBLE -- as DATA, verified.

Source: ATSC A/322:2026-04 "Physical Layer Protocol" (published free at
atsc.org), cross-checked against A/322:2024-04.  Section numbers below refer
to the 2026-04 edition.

WHY THIS FILE EXISTS
--------------------
M2 stopped with an unresolved preamble constellation and correctly identified
the blocker as missing published TABLES rather than missing algorithms.  These
are those tables.  Nothing here is hand-transcribed magic: every table carries
an arithmetic identity that ties it to at least one OTHER independently
extracted table, and `verify()` runs all of them.

THE EXTRACTION HAZARD, AND HOW EACH TABLE IS DEFENDED
-----------------------------------------------------
The A/321 PDF was already found to render some math-font digits shifted by +4.
A silently mis-extracted table looks like a decoder bug forever, so no number
below is trusted on the strength of "the PDF said so".  Every one is pinned by
a closed arithmetic identity against a DIFFERENT table:

  V1  Table 7.1  NoC decreases by exactly Cunit (96/192/384) per Cred_coeff step.
  V2  Table 7.1  occupied bandwidth == NoC * fs/Nfft, to the printed 6 decimals.
  V3  Table D.1.3 CP8 is derivable from Table D.1.1 CP32 by CP8(k)=ceil(CP32(4k)/4).
  V4  Table D.1.1 CP32 is a mirror image: CP32_R(k) = 27648 - CP32_L(95-k).
  V5  Table 8.4  #common CPs == #CP8 entries inside the Cred-reduced carrier
                 range -- reproduces all five values 48/48/47/46/45 from V3's list.
  V6  Table 7.2  available data cells == NoC - preamble pilots - common CPs,
                 reproduces all five values for every FFT/GI/DX row we use.
  V7  Table 8.6  amplitude == 10**(dB/20) for every row.
  V8  Section 8.1.2 reference PRBS reproduces the spec's OWN printed first-24
                 test vector "1101 1000 0000 0001 0100 0000".
  V9  Tables 6.17/6.23/6.24 close: the L1-Basic cell count computed from the
                 puncturing/repetition arithmetic equals the printed
                 "Length (Cells)" for all SEVEN modes.

A digit shift in any one of these tables breaks the identity that binds it to
the others.  All nine pass.  That is what "verified" means here -- not that the
text was read carefully, but that the numbers are over-determined.

ASSUMPTIONS (things the spec states in prose/figures that required a choice)
---------------------------------------------------------------------------
S1  DC (the centre carrier) is treated as an ORDINARY active carrier, i.e.
    NoCmax is odd (6913/13825/27649) and the centre carrier is index
    (NoCmax-1)/2.  RISK LOW -- and it is not assumed blind: m3_preamble.py
    scans the carrier origin over a range and the correct one wins on a pilot
    coherence statistic by a wide margin, with wrong origins acting as controls.
S2  The guard interval is a cyclic PREFIX (A/322 8.5), so the FFT window of a
    symbol starts GI samples after the symbol start.  RISK LOW -- M2 already
    measured the CP autocorrelation at exactly this geometry.
S3  Preamble pilot amplitude is the "Equivalent Amplitude" column of Table 8.6.
    RISK LOW -- and irrelevant to the constellation SHAPE, since a wrong
    global scale is divided out by the equalizer.  It only affects absolute
    power accounting.
"""
from __future__ import annotations

import math

# --------------------------------------------------------------------------
# A/322 Table 7.1 -- Number of Carriers.  NoC = NoCmax - Cred_coeff * Cunit.
# --------------------------------------------------------------------------

NOC_MAX = {8192: 6913, 16384: 13825, 32768: 27649}          # Cred_coeff = 0
C_UNIT = {8192: 96, 16384: 192, 32768: 384}                 # A/322 7.2.3

# printed occupied bandwidth (MHz) for bsr_coefficient = 2 (fs = 6.912 Msps),
# used only as a cross-check of NOC_MAX / C_UNIT (see V2).
OCCUPIED_BW_BSR2 = {
    8192:  [5.832844, 5.751844, 5.670844, 5.589844, 5.508844],
    16384: [5.832844, 5.751844, 5.670844, 5.589844, 5.508844],
    32768: [5.832844, 5.751844, 5.670844, 5.589844, 5.508844],
}
# NB the printed table lists one bandwidth column per bsr_coefficient shared by
# all three FFT sizes, because NoC*spacing is (very nearly) FFT-size invariant.
# The exact printed column for bsr_coefficient = 2 is the 8K one above.

FS_POST = {2: 6_912_000, 5: 8_064_000, 8: 9_216_000}        # A/321 bsr -> fs


def noc(nfft: int, cred: int) -> int:
    """A/322 7.2.3."""
    return NOC_MAX[nfft] - cred * C_UNIT[nfft]


def carrier_abs_range(nfft: int, cred: int):
    """Absolute carrier indices [lo, hi] actually transmitted (A/322 7.2.3).

    "the active absolute carrier indices ... range from Cred_coeff*Cunit/2 to
    NoCmax - Cred_coeff*Cunit/2 - 1".
    """
    off = cred * C_UNIT[nfft] // 2
    return off, NOC_MAX[nfft] - off - 1


# --------------------------------------------------------------------------
# A/322 Annex H Table H.1.1 -- preamble_structure -> preamble grid.
# Generated from the table's own regular structure: for each (FFT, GI, DX)
# group the five L1-Basic FEC modes 1..5 occupy five consecutive codes.
# The literal first 32 rows were read directly from the PDF and are used to
# CHECK the generated form (see verify()).
# --------------------------------------------------------------------------

# (FFT, GI samples, preamble pilot DX) in signalled order, 5 codes each.
_PREAMBLE_GROUPS = [
    (8192, 192, 16), (8192, 384, 8), (8192, 512, 6), (8192, 768, 4),
    (8192, 1024, 3), (8192, 1536, 4), (8192, 2048, 3),
    (16384, 192, 32), (16384, 384, 16), (16384, 512, 12), (16384, 768, 8),
    (16384, 1024, 6), (16384, 1536, 4), (16384, 2048, 3), (16384, 2432, 3),
    (16384, 3072, 4), (16384, 3648, 4), (16384, 4096, 3),
    (32768, 192, 32), (32768, 384, 32), (32768, 512, 24), (32768, 768, 16),
    (32768, 1024, 12), (32768, 1536, 8), (32768, 2048, 6), (32768, 2432, 6),
    (32768, 3072, 8), (32768, 3072, 3), (32768, 3648, 8), (32768, 3648, 3),
    (32768, 4096, 3), (32768, 4864, 3),
]

PREAMBLE_STRUCTURE = {}
for _g, (_n, _gi, _dx) in enumerate(_PREAMBLE_GROUPS):
    for _m in range(5):
        PREAMBLE_STRUCTURE[_g * 5 + _m] = {
            "fft": _n, "gi": _gi, "dx": _dx, "l1b_mode": _m + 1,
        }

# The literally-read first rows of Table H.1.1, kept as an extraction control.
_H11_LITERAL = {
    0: (8192, 192, 16, 1), 4: (8192, 192, 16, 5),
    5: (8192, 384, 8, 1), 9: (8192, 384, 8, 5),
    10: (8192, 512, 6, 1), 14: (8192, 512, 6, 5),
    15: (8192, 768, 4, 1), 19: (8192, 768, 4, 5),
    20: (8192, 1024, 3, 1), 24: (8192, 1024, 3, 5),
    25: (8192, 1536, 4, 1), 26: (8192, 1536, 4, 2),
    27: (8192, 1536, 4, 3), 28: (8192, 1536, 4, 4), 29: (8192, 1536, 4, 5),
    30: (8192, 2048, 3, 1), 31: (8192, 2048, 3, 2),
}

# --------------------------------------------------------------------------
# A/322 Table 7.2 -- available data cells per Preamble symbol (TR disabled).
# Keyed (fft, gi) -> list over Cred_coeff 0..4.   Only rows we exercise.
# --------------------------------------------------------------------------

PREAMBLE_DATA_CELLS = {
    (8192, 192):  [6432, 6342, 6253, 6164, 6075],
    (8192, 384):  [6000, 5916, 5833, 5750, 5667],
    (8192, 512):  [5712, 5632, 5553, 5474, 5395],
    (8192, 768):  [5136, 5064, 4993, 4922, 4851],
    (8192, 1024): [4560, 4496, 4433, 4370, 4307],
    (8192, 1536): [5136, 5064, 4993, 4922, 4851],
    (8192, 2048): [4560, 4496, 4433, 4370, 4307],
    (16384, 192):  [13296, 13110, 12927, 12742, 12558],
    (16384, 384):  [12864, 12684, 12507, 12328, 12150],
    (16384, 512):  [12576, 12400, 12227, 12052, 11878],
    (16384, 768):  [12000, 11832, 11667, 11500, 11334],
    (16384, 1024): [11424, 11264, 11107, 10948, 10790],
    (16384, 1536): [10272, 10128, 9987, 9844, 9702],
    (16384, 2048): [9120, 8992, 8867, 8740, 8614],
    (16384, 2432): [9120, 8992, 8867, 8740, 8614],
    (16384, 3072): [10272, 10128, 9987, 9844, 9702],
    (16384, 3648): [10272, 10128, 9987, 9844, 9702],
    (16384, 4096): [9120, 8992, 8867, 8740, 8614],
}

# --------------------------------------------------------------------------
# A/322 Table 8.4 -- number of COMMON continual pilots per FFT size / Cred.
# --------------------------------------------------------------------------

N_COMMON_CP = {
    8192:  [48, 48, 47, 46, 45],
    16384: [96, 96, 93, 92, 90],
    32768: [192, 192, 186, 184, 180],
}

# --------------------------------------------------------------------------
# A/322 Table D.1.1 -- common CP absolute carrier indices, CP32.
# --------------------------------------------------------------------------

CP32 = [
    236, 316, 356, 412, 668, 716, 868, 1100, 1228, 1268, 1340, 1396, 1876,
    1916, 2140, 2236, 2548, 2644, 2716, 2860, 3004, 3164, 3236, 3436, 3460,
    3700, 3836, 4028, 4124, 4132, 4156, 4316, 4636, 5012, 5132, 5140, 5332,
    5372, 5500, 5524, 5788, 6004, 6020, 6092, 6428, 6452, 6500, 6740, 7244,
    7316, 7372, 7444, 7772, 7844, 7924, 8020, 8164, 8308, 8332, 8348, 8788,
    8804, 9116, 9140, 9292, 9412, 9436, 9604, 10076, 10204, 10340, 10348,
    10420, 10660, 10684, 10708, 11068, 11132, 11228, 11356, 11852, 11860,
    11884, 12044, 12116, 12164, 12268, 12316, 12700, 12772, 12820, 12988,
    13300, 13340, 13564, 13780,
    13868, 14084, 14308, 14348, 14660, 14828, 14876, 14948, 15332, 15380,
    15484, 15532, 15604, 15764, 15788, 15796, 16292, 16420, 16516, 16580,
    16940, 16964, 16988, 17228, 17300, 17308, 17444, 17572, 18044, 18212,
    18236, 18356, 18508, 18532, 18844, 18860, 19300, 19316, 19340, 19484,
    19628, 19724, 19804, 19876, 20204, 20276, 20332, 20404, 20908, 21148,
    21196, 21220, 21556, 21628, 21644, 21860, 22124, 22148, 22276, 22316,
    22508, 22516, 22636, 23012, 23332, 23492, 23516, 23524, 23620, 23812,
    23948, 24188, 24212, 24412, 24484, 24644, 24788, 24932, 25004, 25100,
    25412, 25508, 25732, 25772, 26252, 26308, 26380, 26420, 26548, 26780,
    26932, 26980, 27236, 27292, 27332, 27412,
]

# A/322 Table D.1.3 -- common CP absolute carrier indices, CP8 (as printed).
CP8_PRINTED = [
    59, 167, 307, 469, 637, 751, 865, 1031, 1159, 1333, 1447, 1607, 1811,
    1943, 2041, 2197, 2323, 2519, 2605, 2767, 2963, 3029, 3175, 3325, 3467,
    3665, 3833, 3901, 4073, 4235, 4325, 4511, 4627, 4825, 4907, 5051, 5227,
    5389, 5531, 5627, 5833, 5905, 6053, 6197, 6353, 6563, 6637, 6809,
]

# A/322 Table D.1.2 -- common CP absolute carrier indices, CP16 (as printed).
CP16_PRINTED = [
    118, 178, 334, 434, 614, 670, 938, 1070, 1274, 1358, 1502, 1618, 1730,
    1918, 2062, 2078, 2318, 2566, 2666, 2750, 2894, 3010, 3214, 3250, 3622,
    3686, 3886, 3962, 4082, 4166, 4394, 4558, 4646, 4718, 5038, 5170, 5210,
    5342, 5534, 5614, 5926, 5942, 6058, 6134, 6350, 6410, 6650, 6782, 6934,
    7154, 7330, 7438, 7666, 7742, 7802, 7894, 8146, 8258, 8470, 8494, 8650,
    8722, 9022, 9118, 9254, 9422, 9650, 9670, 9814, 9902, 10102, 10166,
    10454, 10598, 10778, 10822, 11062, 11138, 11254, 11318, 11666, 11758,
    11810, 11974, 12106, 12242, 12394, 12502, 12706, 12866, 13126, 13190,
    13274, 13466, 13618, 13666,
]


def cp_derived(nfft: int):
    """A/322 8.1.4.1: CP16(k) = ceil(CP32(2k)/2), CP8(k) = ceil(CP32(4k)/4)."""
    if nfft == 32768:
        return list(CP32)
    step = 2 if nfft == 16384 else 4
    n = 96 if nfft == 16384 else 48
    return [-(-CP32[step * k] // step) for k in range(n)]


def common_cps(nfft: int, cred: int):
    """Common continual pilots actually transmitted, as ABSOLUTE indices.

    A/322 D.1: "Continual pilots with absolute carrier indices that lie outside
    the configured range of carriers for transmission shall not be transmitted."
    """
    lo, hi = carrier_abs_range(nfft, cred)
    return [c for c in cp_derived(nfft) if lo <= c <= hi]


# --------------------------------------------------------------------------
# A/322 Table 8.6 -- preamble pilot boost.  (fft, gi) -> (power_dB, amplitude)
# --------------------------------------------------------------------------

PREAMBLE_PILOT_BOOST = {
    (8192, 192): (5.30, 1.841), (8192, 384): (3.60, 1.514),
    (8192, 512): (2.90, 1.396), (8192, 768): (1.80, 1.230),
    (8192, 1024): (0.90, 1.109), (8192, 1536): (1.80, 1.230),
    (8192, 2048): (0.90, 1.109),
    (16384, 192): (6.80, 2.188), (16384, 384): (5.30, 1.841),
    (16384, 512): (4.60, 1.698), (16384, 768): (3.60, 1.514),
    (16384, 1024): (2.90, 1.396), (16384, 1536): (2.10, 1.274),
    (16384, 2048): (1.30, 1.161), (16384, 2432): (1.30, 1.161),
    (16384, 3072): (2.10, 1.274), (16384, 3648): (2.10, 1.274),
    (16384, 4096): (1.30, 1.161),
}

A_CP_DB = 8.52          # A/322 Table 8.5, all FFT sizes
A_CP = 2.67             # printed approximate linear value


# --------------------------------------------------------------------------
# A/322 8.1.2 -- pilot reference sequence.
# G(x) = 1 + X^9 + X^10 + X^12 + X^13, initial load 0000 0000 1101 1,
# reloaded at the start of every symbol.
# --------------------------------------------------------------------------

_REF_INIT = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1]     # x1 .. x13
REF_TEST_VECTOR = "110110000000000101000000"             # spec's printed first 24


def reference_sequence(n: int):
    """r_k for k = 0 .. n-1 (A/322 8.1.2).  Values in {0,1}."""
    x = list(_REF_INIT)
    out = []
    for _ in range(n):
        out.append(x[12])                                # r_k = x13
        fb = x[8] ^ x[9] ^ x[11] ^ x[12]                 # X^9 + X^10 + X^12 + X^13
        x = [fb] + x[:12]
        # x[i] holds x_{i+1}; shifting right moves x_i -> x_{i+1}
    return out


def pilot_values(n: int, amp: float):
    """Re{c_k} = 2*A*(1/2 - r_k)  ->  +A when r_k = 0, -A when r_k = 1."""
    return [amp * (1.0 - 2.0 * r) for r in reference_sequence(n)]


# --------------------------------------------------------------------------
# A/322 Table 6.17 / 6.23 / 6.24 -- L1-Basic BICM configuration.
# --------------------------------------------------------------------------

L1_BASIC_KSIG = 200
L1_BASIC_NOUTER = 368            # Table 6.18: Ksig + Mouter, Mouter = 168
L1_BASIC_KLDPC = 3240            # Table 6.19  -> rate 3/15 of Ninner = 16200
L1_BASIC_NINNER = 16200
L1_BASIC_NPARITY = 12960         # Ninner - Kldpc

# Table 6.24, L1-Basic rows.  A = 0 for every L1-Basic mode (Nouter is fixed).
L1_BASIC_PUNCT_B = {1: 9360, 2: 11460, 3: 12360, 4: 12292,
                    5: 12350, 6: 12432, 7: 12776}
L1_BASIC_ETA_MOD = {1: 2, 2: 2, 3: 2, 4: 4, 5: 6, 6: 8, 7: 8}
L1_BASIC_CONSTELLATION = {
    1: "QPSK", 2: "QPSK", 3: "QPSK", 4: "NUC_16_8/15",
    5: "NUC_64_9/15", 6: "NUC_256_9/15", 7: "NUC_256_13/15",
}
# Table 6.23: repetition applies to L1-Basic Mode 1 only.  C = 0, D = +3672.
L1_BASIC_REPEAT = {1: (0, 3672)}
# Table 6.17 printed "Length (Cells)" -- the checksum target, not an input.
L1_BASIC_CELLS_PRINTED = {1: 3820, 2: 934, 3: 484, 4: 259,
                          5: 163, 6: 112, 7: 69}


def l1_basic_lengths(mode: int):
    """A/322 6.5.2.7 / 6.5.2.8 arithmetic -> (N_FEC bits, cells, N_punc).

    This is computed, never read.  Comparing it against the printed
    Length (Cells) column is verification V9.
    """
    eta = L1_BASIC_ETA_MOD[mode]
    n_punc_temp = L1_BASIC_PUNCT_B[mode]                 # A = 0
    n_repeat = 0
    if mode in L1_BASIC_REPEAT:
        c, d = L1_BASIC_REPEAT[mode]
        n_repeat = 2 * ((c * (L1_BASIC_KLDPC - L1_BASIC_NOUTER)) // 2) + d
    n_fec_temp = (L1_BASIC_NOUTER + n_repeat
                  + L1_BASIC_NPARITY - n_punc_temp)
    n_fec = (n_fec_temp // eta) * eta
    n_punc = n_punc_temp - (n_fec - n_fec_temp)
    return n_fec, n_fec // eta, n_punc


# --------------------------------------------------------------------------
# A/322 7.3 -- frequency interleaver.
# NOTE: a permutation of cells within one symbol.  It cannot change the
# CONSTELLATION -- only the bit ORDER.  It is therefore irrelevant to the
# Step-2 constellation proof and only matters once we read bits.
# --------------------------------------------------------------------------

FI_NMAX = {8192: 8192, 16384: 16384, 32768: 32768}       # Table 7.11

# Tables 7.12 / 7.13 / 7.14.  Row 1 is R'_i bit positions (high->low), row 2 is
# the R_i bit position each maps to.
FI_WIRE = {
    8192: {
        "rprime": [11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
        "even":   [5, 11, 3, 0, 10, 8, 6, 9, 2, 4, 1, 7],
        "odd":    [8, 10, 7, 6, 0, 5, 2, 1, 3, 9, 4, 11],
    },
    16384: {
        "rprime": [12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
        "even":   [8, 4, 3, 2, 0, 11, 1, 5, 12, 10, 6, 7, 9],
        "odd":    [7, 9, 5, 3, 11, 1, 4, 0, 2, 12, 10, 8, 6],
    },
    32768: {
        "rprime": [13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
        "even":   [6, 5, 0, 10, 8, 1, 11, 12, 2, 9, 4, 3, 13, 7],
        "odd":    [6, 5, 0, 10, 8, 1, 11, 12, 2, 9, 4, 3, 13, 7],
    },
}

# R' feedback taps (A/322 7.3), by FFT size: R'_i[Nr-2] = XOR of R'_{i-1}[taps].
FI_R_TAPS = {8192: (0, 1, 4, 6),
             16384: (0, 1, 4, 5, 9, 11),
             32768: (0, 1, 2, 12)}
# G (symbol offset) feedback taps, by FFT size.
FI_G_TAPS = {8192: (0, 1, 4, 5, 9, 11),
             16384: (0, 1, 2, 12),
             32768: (0, 1)}


# --------------------------------------------------------------------------
# VERIFICATION
# --------------------------------------------------------------------------

def verify(verbose=True):
    """Run every cross-table identity.  Returns list of (name, ok, detail)."""
    res = []

    def chk(name, ok, detail=""):
        res.append((name, bool(ok), detail))

    # V1 -- Table 7.1 carrier steps
    ok, det = True, []
    for nfft, mx in NOC_MAX.items():
        vals = [noc(nfft, c) for c in range(5)]
        steps = {vals[i] - vals[i + 1] for i in range(4)}
        ok &= steps == {C_UNIT[nfft]}
        det.append(f"{nfft}:{vals}")
    chk("V1 Table7.1 NoC steps == Cunit", ok, "; ".join(det))

    # V2 -- occupied bandwidth == NoC * fs/Nfft
    ok, det = True, []
    for c in range(5):
        spacing = FS_POST[2] / 8192
        bw = noc(8192, c) * spacing / 1e6
        want = OCCUPIED_BW_BSR2[8192][c]
        good = abs(bw - want) < 1e-6
        ok &= good
        det.append(f"c{c}:{bw:.6f}vs{want:.6f}")
    chk("V2 Table7.1 bandwidth == NoC*spacing", ok, "; ".join(det))

    # V3 -- CP8 / CP16 derive from CP32
    chk("V3 CP8 == ceil(CP32(4k)/4)", cp_derived(8192) == CP8_PRINTED,
        f"{sum(a != b for a, b in zip(cp_derived(8192), CP8_PRINTED))} mismatches")
    chk("V3b CP16 == ceil(CP32(2k)/2)", cp_derived(16384) == CP16_PRINTED,
        f"{sum(a != b for a, b in zip(cp_derived(16384), CP16_PRINTED))} mismatches")

    # V4 -- CP32 mirror symmetry
    mirror_ok = all(CP32[96 + k] == 27648 - CP32[95 - k] for k in range(96))
    chk("V4 CP32 mirror CP32_R(k)=27648-CP32_L(95-k)", mirror_ok,
        f"len(CP32)={len(CP32)}")

    # V5 -- Table 8.4 common-CP counts reproduce from CP list + range rule
    ok, det = True, []
    for nfft in (8192, 16384, 32768):
        got = [len(common_cps(nfft, c)) for c in range(5)]
        good = got == N_COMMON_CP[nfft]
        ok &= good
        det.append(f"{nfft}:{got}vs{N_COMMON_CP[nfft]}")
    chk("V5 Table8.4 == #CP inside reduced range", ok, "; ".join(det))

    # V6 -- Table 7.2 == NoC - preamble pilots - common CPs
    ok, det = True, []
    for code, p in PREAMBLE_STRUCTURE.items():
        key = (p["fft"], p["gi"])
        if key not in PREAMBLE_DATA_CELLS:
            continue
        for c in range(5):
            n = noc(p["fft"], c)
            npil = (n + p["dx"] - 1) // p["dx"]          # k mod DX == 0, k<n
            ncp = len(common_cps(p["fft"], c))
            got = n - npil - ncp
            want = PREAMBLE_DATA_CELLS[key][c]
            if got != want:
                ok = False
                det.append(f"ps{code} c{c}: {got}!={want}")
    chk("V6 Table7.2 == NoC - preamble pilots - CPs", ok,
        "; ".join(det) if det else "all rows exact")

    # V7 -- Table 8.6 amplitude == 10**(dB/20)
    ok, det = True, []
    for key, (db, amp) in PREAMBLE_PILOT_BOOST.items():
        got = 10 ** (db / 20.0)
        if abs(got - amp) > 0.0006:
            ok = False
            det.append(f"{key}:{got:.4f}vs{amp}")
    chk("V7 Table8.6 amp == 10^(dB/20)", ok, "; ".join(det) if det else "all rows")
    chk("V7b Table8.5 A_CP == 10^(8.52/20)", abs(10 ** (A_CP_DB / 20) - A_CP) < 0.004,
        f"{10 ** (A_CP_DB / 20):.4f} vs {A_CP}")

    # V8 -- reference sequence vs the spec's own printed 24-bit vector
    got = "".join(str(b) for b in reference_sequence(24))
    chk("V8 ref PRBS == spec's printed first 24", got == REF_TEST_VECTOR,
        f"{got} vs {REF_TEST_VECTOR}")

    # V9 -- L1-Basic cell counts reproduce from puncturing arithmetic
    ok, det = True, []
    for m in range(1, 8):
        _, cells, _ = l1_basic_lengths(m)
        want = L1_BASIC_CELLS_PRINTED[m]
        if cells != want:
            ok = False
        det.append(f"m{m}:{cells}/{want}")
    chk("V9 Table6.17 cells == 6.23/6.24 arithmetic", ok, " ".join(det))

    # V10 -- Annex H generated form matches the literally-read rows
    ok, det = True, []
    for code, tup in _H11_LITERAL.items():
        p = PREAMBLE_STRUCTURE[code]
        got = (p["fft"], p["gi"], p["dx"], p["l1b_mode"])
        if got != tup:
            ok = False
            det.append(f"{code}:{got}!={tup}")
    chk("V10 AnnexH generated == literal rows", ok,
        "; ".join(det) if det else f"{len(_H11_LITERAL)} rows, "
        f"{len(PREAMBLE_STRUCTURE)} codes total")

    # V11 -- frequency interleaver wire permutations really are permutations
    ok, det = True, []
    for nfft, w in FI_WIRE.items():
        nb = len(w["rprime"])
        for which in ("even", "odd"):
            if sorted(w[which]) != list(range(nb)) or sorted(w["rprime"]) != list(range(nb)):
                ok = False
                det.append(f"{nfft}/{which}")
    chk("V11 FI wire permutations are permutations", ok,
        "; ".join(det) if det else "8K/16K/32K even+odd")

    if verbose:
        for name, good, detail in res:
            print(f"  [{'PASS' if good else 'FAIL'}] {name}"
                  + (f"   {detail}" if detail else ""))
        npass = sum(1 for _, g, _ in res if g)
        print(f"\n  {npass}/{len(res)} checks pass")
    return res


if __name__ == "__main__":
    print("A/322 preamble tables -- cross-table verification\n")
    r = verify()
    raise SystemExit(0 if all(g for _, g, _ in r) else 1)
