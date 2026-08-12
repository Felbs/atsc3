#!/usr/bin/env python3
"""M4 Step 1 -- the A/322 5.2.3 scrambler, SOLVED and pinned from the air.

M3 stopped here and called A/322 Figure 5.6 an erratum.  **That was wrong, and
the correction is worth more than the original finding.**

WHAT M3 GOT RIGHT
-----------------
Reading Figure 5.6's glyph columns gives  D0..D7 = X1, X3, X4, X7, X11, X12,
X13, X14, which reproduces the spec's printed FIRST byte 8/8.  M3 then proved
that under a plain FIBONACCI shift, stage X4 of the next state must equal
whatever X3 held (=0), so byte 1's D2 bit is forced to 0 while the printed
vector needs 1 -- and concluded the figure and the vector contradict.

WHAT M3 GOT WRONG
-----------------
That argument assumes the register update is a plain shift with feedback into
one end -- a FIBONACCI LFSR.  It is not.  A/322 5.2.3's register is a **GALOIS**
LFSR: the feedback bit is XORed INTO the tapped stages as the contents move.
Under Galois the next X4 is X3 XOR feedback, so it is NOT forced, and the
contradiction evaporates.  The wall was an assumed convention, not a spec
defect.  (Record this: "the spec is wrong" is the most expensive conclusion
available and must be the last one reached, not the first.)

THE CONVENTION SPACE, ENUMERATED RATHER THAN CHOSEN
---------------------------------------------------
Unknowns: shift direction, Fibonacci vs Galois, printed vs reciprocal tap
exponents, seed orientation.  16 configurations.  Oracle: the spec's own
printed 24-bit sequence 110000000110110100111111.  Exactly ONE configuration
admits any 8-distinct-stage output tap set:

    shift DOWN (X_i <- X_{i+1}),  GALOIS,  taps = reciprocal exponents
    [3,4,5,9,10,13,15,16],  seed 0xF180 loaded with X1 = MSB

...which is precisely the mirror image (j = 17 - i) of the figure's own
coordinates, and `mirror_check()` below shows that in figure coordinates the
solved output stages are D0..D7 = X1, X3, X4, X7, X11, X12, X13, X14 --
**M3's glyph-column reading, exactly.**  The figure's TAPS were right all
along; only the register update was misread.

STEP 1 -- PINNING THE OUTPUT TAPS FROM THE AIR
-----------------------------------------------
24 printed bits give only a 3-bit history per output, so they underdetermine
which stages feed D7..D0: 16875 assignments fit (2160 of them with 8 distinct
stages).  The air settles it.  A/322 6.1.2.2's CRC-32 has a fixed all-ones
init, hence is affine, hence

    g(d) = K  for every valid L1-Basic block d,   d = r XOR mask

is a 32-bit equation the true mask must satisfy and a wrong one satisfies with
probability 2^-32.  Against 16875 candidates that is 4e-6 expected false
survivors.

RESULT: exactly one survivor.  D0..D7 = X16, X14, X13, X10, X6, X5, X4, X3.

Then a SECOND, fully independent test that nothing in the search touched:
descramble and compare the field values against quantities M2 measured off the
waveform months of DSP before any table existed here, and against A/322 3.2.1's
"the ATSC default value for reserved bits is 1" -- 47 reserved bits that the
CRC equation does not constrain.

Usage:
    python m4_scrambler.py                 # uses the banked M3 blocks
    python m4_scrambler.py --capture hit_rf33.cs16 --rate 8e6   # re-decode
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m3_crc as CRC                                            # noqa: E402
import spec_l1syntax as LS                                      # noqa: E402

SEED = 0xF180
POLY_EXP_PRINTED = (1, 3, 6, 7, 11, 12, 13, 16)     # G(x)=1+X+X^3+X^6+X^7+X^11+X^12+X^13+X^16
POLY_EXP_RECIP = tuple(sorted({16 - e for e in POLY_EXP_PRINTED
                               if 0 < 16 - e <= 16} | {16}))
TEST_VECTOR = "110000000110110100111111"
FIG_PRELOAD = [0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1]   # X1..X16

# --- THE SOLVED CONFIGURATION ---------------------------------------------
SOLVED = dict(direction="down", form="galois", taps=POLY_EXP_RECIP,
              seed="X1=MSB")
# D0..D7 -> 1-based stage number.  Pinned by the air (see solve_output_taps()).
SOLVED_OUT = (16, 14, 13, 10, 6, 5, 4, 3)


# --------------------------------------------------------------------------
# the LFSR
# --------------------------------------------------------------------------

def _seed_reg(orient):
    b = [(SEED >> i) & 1 for i in range(16)]                     # b[0] = LSB
    if orient == "X1=LSB":
        return [b[i] for i in range(16)]
    return [b[15 - i] for i in range(16)]


def _step(s, direction, taps, form):
    """s is a list indexed 0..15 for stages X1..X16."""
    n = list(s)
    if form == "fib":
        fb = 0
        for t in taps:
            fb ^= s[t - 1]
        if direction == "up":
            n = [fb] + s[:-1]
        else:
            n = s[1:] + [fb]
    else:                                                       # galois
        if direction == "down":
            out = s[0]
            for i in range(1, 16):
                n[i - 1] = s[i] ^ (out if i in taps else 0)
            n[15] = out
        else:
            out = s[15]
            for i in range(16, 1, -1):
                n[i - 1] = s[i - 2] ^ (out if i in taps else 0)
            n[0] = out
    return n


def states(nbytes, direction=None, taps=None, form=None, seed=None,
           preload=None):
    """State of X1..X16 at the start of each output byte."""
    direction = direction or SOLVED["direction"]
    taps = taps or SOLVED["taps"]
    form = form or SOLVED["form"]
    s = list(preload) if preload is not None else _seed_reg(seed
                                                            or SOLVED["seed"])
    out = []
    for _ in range(nbytes):
        out.append(list(s))
        s = _step(s, direction, taps, form)
    return np.array(out, np.uint8)


def sequence(nbits, out_taps=SOLVED_OUT, **kw):
    """The scrambling sequence, MSB-first (D7, D6, ..., D0, D7, ...)."""
    nb = (nbits + 7) // 8
    st = states(nb, **kw)
    m = np.empty(nb * 8, np.uint8)
    for b in range(nb):
        for q in range(8):                    # q = 0 is the MSB of the byte
            m[8 * b + q] = st[b][out_taps[7 - q] - 1]
    return m[:nbits]


def descramble(bits, out_taps=SOLVED_OUT, **kw):
    bits = np.asarray(bits, np.uint8)
    return bits ^ sequence(len(bits), out_taps=out_taps, **kw)


# --------------------------------------------------------------------------
# gate 1 -- the spec's own printed vector, and the convention enumeration
# --------------------------------------------------------------------------

def enumerate_conventions(verbose=True):
    """16 register conventions vs the printed 24-bit vector."""
    rows = []
    for direction in ("up", "down"):
        for form in ("fib", "galois"):
            for tname, taps in (("printed", POLY_EXP_PRINTED),
                                ("reciprocal", POLY_EXP_RECIP)):
                for seed in ("X1=LSB", "X1=MSB"):
                    st = states(3, direction, taps, form, seed)
                    fits = []
                    for d in range(8):
                        want = tuple(int(TEST_VECTOR[t * 8 + (7 - d)])
                                     for t in range(3))
                        fits.append([s for s in range(16)
                                     if tuple(st[:3, s]) == want])
                    n_any = 1
                    for f in fits:
                        n_any *= len(f)
                    n_dist = sum(1 for c in itertools.product(*fits)
                                 if len(set(c)) == 8) if n_any else 0
                    rows.append(dict(direction=direction, form=form,
                                     taps=tname, seed=seed,
                                     n_assignments=n_any, n_distinct=n_dist))
    if verbose:
        print("  === A/322 5.2.3 register conventions vs the PRINTED 24-bit "
              "vector ===")
        for r in rows:
            ok = r["n_distinct"] > 0
            print(f"    shift {r['direction']:4s}  {r['form']:7s}  taps "
                  f"{r['taps']:10s}  seed {r['seed']:6s}  -> "
                  f"{r['n_distinct']:5d} distinct-stage tap sets"
                  f"{'   <== SOLVES' if ok else ''}")
        n = sum(1 for r in rows if r["n_distinct"])
        print(f"    {n} of {len(rows)} configurations reproduce the vector.")
    return rows


def mirror_check(verbose=True):
    """The solved config, relabelled j = 17 - i, IS Figure 5.6's own reading.

    This is the independent confirmation that matters: the air pinned eight
    output stages with no knowledge of the figure, and under the mirror
    relabelling they land exactly on the glyph-column taps M3 read off the PDF.
    """
    taps_m = tuple(sorted(17 - t for t in SOLVED["taps"]))
    out_m = tuple(17 - s for s in SOLVED_OUT)
    st = states(3, "up", taps_m, "galois", preload=FIG_PRELOAD)
    seq = "".join(str(st[b][out_m[7 - q] - 1])
                  for b in range(3) for q in range(8))
    m3_glyph = (1, 3, 4, 7, 11, 12, 13, 14)
    if verbose:
        print("\n  === mirror check (j = 17 - i): the figure's own "
              "coordinates ===")
        print(f"    seed preload           : Figure 5.6's printed row "
              f"(X1..X16), which reads X16->X1 as 0xF180")
        print(f"    shift UP, GALOIS, taps : {list(taps_m)}")
        print(f"    solved output stages   : D0..D7 = "
              f"{', '.join('X%d' % s for s in out_m)}")
        print(f"    M3 glyph-column read   : D0..D7 = "
              f"{', '.join('X%d' % s for s in m3_glyph)}")
        print(f"    identical              : {out_m == m3_glyph}")
        print(f"    reproduces the printed 24-bit vector: "
              f"{seq == TEST_VECTOR}")
    return out_m == m3_glyph and seq == TEST_VECTOR


# --------------------------------------------------------------------------
# gate 2 -- pin the output taps against the air
# --------------------------------------------------------------------------

def candidate_taps():
    """Every output-stage assignment consistent with the printed 24 bits."""
    st = states(3)
    fits = []
    for d in range(8):
        want = tuple(int(TEST_VECTOR[t * 8 + (7 - d)]) for t in range(3))
        fits.append([s + 1 for s in range(16) if tuple(st[:3, s]) == want])
    return fits


def solve_output_taps(blocks, verbose=True):
    """Filter the candidates with the air's CRC-32.  Repeats ALLOWED."""
    fits = candidate_taps()
    K = CRC.crc32(np.zeros(168, np.uint8), init=0xFFFFFFFF)
    st = states(25)
    surv, n_all, n_dist = [], 0, 0
    for combo in itertools.product(*fits):
        n_all += 1
        n_dist += len(set(combo)) == 8
        m = np.empty(200, np.uint8)
        for b in range(25):
            for q in range(8):
                m[8 * b + q] = st[b][combo[7 - q] - 1]
        if all(CRC.g_stat(np.asarray(bl, np.uint8) ^ m) == K for bl in blocks):
            surv.append(tuple(combo))
    if verbose:
        print("\n  === pinning the eight output stages against the AIR ===")
        print(f"    stages consistent with the printed vector, per output:")
        for d in range(8):
            print(f"      D{d}: {', '.join('X%d' % s for s in fits[d])}")
        print(f"    total assignments tested   : {n_all} "
              f"({n_dist} with 8 distinct stages)")
        print(f"    oracle                     : A/322 6.1.2.2 CRC-32 is "
              f"affine, so g(r XOR mask) = K is 32 hard bits")
        print(f"    expected false survivors   : {n_all / 2**32:.2e}")
        print(f"    SURVIVORS                  : {len(surv)}")
        for c in surv:
            print(f"      D0..D7 = {', '.join('X%d' % s for s in c)}"
                  f"   (distinct: {len(set(c)) == 8})")
    return surv


# --------------------------------------------------------------------------
# gate 3 -- independent measurements, and the reserved-bits test
# --------------------------------------------------------------------------

def parse_l1b(bits):
    out = {}

    def rd(items, i):
        for it in items:
            if it[0] == "field":
                _, name, w = it
                v = 0
                for k in range(w):
                    v = (v << 1) | int(bits[i + k])
                out[name] = v
                i += w
            else:
                i = rd(it[1][out["L1B_frame_length_mode"]], i)
        return i
    n = rd(LS.L1_BASIC_STRUCT, 0)
    assert n == 200, n
    return out


# What M2 measured off the WAVEFORM, before any A/322 table existed in the repo
INDEPENDENT = [
    ("L1B_first_sub_fft_size", 0,
     "8K  -- M2 folded the CP autocorrelation to an 8192-point grid"),
    ("L1B_first_sub_guard_interval", 6,
     "GI6_1536 -- M2 measured CP lag 1536, illegal GIs searched alongside"),
    ("L1B_first_sub_scattered_pilot_pattern", 2,
     "SP4_2 -- M2 measured D_x=4, D_y=2 blind by repeat coherence"),
    ("L1B_first_sub_reduced_carriers", 0,
     "Cred 0 -> NoC 6913 -- M2 counted 6921 active carriers"),
    ("L1B_num_subframes", 1,
     "1 => 2 subframes -- M2 found the grid change at 36 x 9728 samples"),
    ("L1B_first_sub_num_ofdm_symbols", 34,
     "34 => 35 data symbols; +1 preamble = 36 -- M2 counted 36"),
    ("L1B_preamble_num_symbols", 0,
     "0 => 1 preamble symbol -- M2's grid starts data at 1 symbol in"),
    ("L1B_reserved", (1 << 47) - 1,
     "all 47 bits 1 -- A/322 3.2.1 default; the CRC does not constrain these"),
]


def l1_detail_cells_predicted(size_bytes, fec_type_mode):
    """A/322 6.5.2.8 arithmetic: L1-Detail cell count from size + FEC mode.

    Completely independent of L1B_L1_Detail_total_cells, which is a DIFFERENT
    field in the same block.  If both agree the decode is internally closed.
    """
    import spec_bicm as BI
    from fractions import Fraction
    p = BI.PUNCTURING["L1-Detail-%d" % fec_type_mode]
    kldpc = BI.L1D_KLDPC[fec_type_mode]
    kseg = BI.L1D_KSEG[fec_type_mode]
    ksig = size_bytes * 8
    nblocks = 1 if ksig <= kseg else -1          # segmentation, 6.5.2.1
    nouter = ksig + 168                          # BCH parity is 168 bits
    A, B = p["A"], p["B"]
    npunc_t = int(Fraction(A) * (kldpc - nouter)) + B
    nfec_t = nouter + p["Nldpc_parity"] - npunc_t
    eta = p["eta"]
    nfec = (nfec_t // eta) * eta
    nrep = 0
    r = BI.REPETITION.get("L1-Detail-%d" % fec_type_mode)
    if r is not None:
        nrep = 2 * int(Fraction(r["C"]) * nouter) + r["D"]
    return dict(Ksig=ksig, Nblocks=nblocks, Nouter=nouter, Kldpc=kldpc,
                Npunc_temp=npunc_t, Nfec=nfec, Nrepeat=nrep, eta=eta,
                cells=(nfec + nrep) // eta * nblocks)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blocks", default="m3_l1basic_hit_rf33.json")
    ap.add_argument("--json", default="m4_scrambler.json")
    a = ap.parse_args()

    print(__doc__.split("Usage:")[0].split("M4 Step 1 --")[1].strip()[:0] or "")
    print("A/322 5.2.3 SCRAMBLER -- solve and pin\n" + "=" * 70)

    conv = enumerate_conventions()
    mir = mirror_check()

    src = a.blocks if os.path.isabs(a.blocks) else os.path.join(HERE, a.blocks)
    raw = json.load(open(src))["blocks"]
    blocks = [np.array([int(c) for c in b[:200]], np.uint8) for b in raw]
    print(f"\n  {len(blocks)} LDPC-converged, BCH-clean L1-Basic blocks from "
          f"{os.path.basename(src)}")

    surv = solve_output_taps(blocks)
    if len(surv) != 1:
        print(f"\n  AMBIGUOUS: {len(surv)} survivors -- reporting rather than "
              f"picking.")
        return 1
    taps = surv[0]
    assert taps == SOLVED_OUT, (taps, SOLVED_OUT)

    print("\n  === INDEPENDENT test: descrambled fields vs what M2 MEASURED "
          "===")
    f0 = parse_l1b(descramble(blocks[0]))
    hits = 0
    for name, want, why in INDEPENDENT:
        got = f0[name]
        ok = got == want
        hits += ok
        print(f"    {name:38s} {got:>16d} vs {want:>16d}  "
              f"{'MATCH' if ok else 'NO':5s}  {why}")
    print(f"\n    {hits}/{len(INDEPENDENT)} independent predictions matched.")

    print("\n  === INTERNAL closure: two different L1-Basic fields, one "
          "spec identity ===")
    mode = f0["L1B_L1_Detail_fec_type"] + 1
    pred = l1_detail_cells_predicted(f0["L1B_L1_Detail_size_bytes"], mode)
    print(f"    L1B_L1_Detail_size_bytes  = "
          f"{f0['L1B_L1_Detail_size_bytes']} -> Ksig {pred['Ksig']}, "
          f"Nouter {pred['Nouter']}")
    print(f"    L1B_L1_Detail_fec_type    = "
          f"{f0['L1B_L1_Detail_fec_type']} -> Mode {mode}, Kldpc "
          f"{pred['Kldpc']}, eta {pred['eta']}")
    print(f"    A/322 6.5.2.8 arithmetic  -> Nfec {pred['Nfec']}, "
          f"cells {pred['cells']}")
    print(f"    L1B_L1_Detail_total_cells = "
          f"{f0['L1B_L1_Detail_total_cells']}   "
          f"{'MATCH' if pred['cells'] == f0['L1B_L1_Detail_total_cells'] else 'MISMATCH'}")

    print("\n  === per-frame: only the time stamp and its CRC may move ===")
    fs = [parse_l1b(descramble(b)) for b in blocks]
    for i, f in enumerate(fs):
        d = {k: v for k, v in f.items() if v != fs[0][k]}
        print(f"    frame {i}: {d if i else '(reference)'}")
    to = [f["L1B_time_offset"] for f in fs]
    dd = [(to[i + 1] - to[i]) % 65536 for i in range(len(to) - 1)]
    print(f"    L1B_time_offset: {to}")
    print(f"    consecutive deltas mod 2^16: {dd}  "
          f"{'CONSTANT' if len(set(dd)) == 1 else 'varying'}")

    out = {"conventions": conv, "mirror_is_figure56": bool(mir),
           "survivors": [list(s) for s in surv],
           "solved": {"direction": SOLVED["direction"], "form": SOLVED["form"],
                      "taps": list(SOLVED["taps"]), "seed": SOLVED["seed"],
                      "D0_to_D7_stages": list(SOLVED_OUT)},
           "independent_hits": hits, "n_independent": len(INDEPENDENT),
           "l1_detail_cells_predicted": pred["cells"],
           "fields": {k: int(v) for k, v in f0.items()},
           "descrambled_blocks": ["".join(str(int(c)) for c in descramble(b))
                                  for b in blocks]}
    dest = a.json if os.path.isabs(a.json) else os.path.join(HERE, a.json)
    with open(dest, "w") as fh:
        json.dump(out, fh, indent=2, default=float)
    print(f"\n  wrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
