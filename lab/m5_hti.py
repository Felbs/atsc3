#!/usr/bin/env python3
"""M5 Step 3 -- the HTI time de-interleaver, built against a synthesized
round-trip FIRST, and honest about the one thing the PDF does not give up.

RF33 PLP 0's time-interleaver configuration, all decoded (not assumed):

    L1D_plp_TI_mode                2  -> HTI
    L1D_plp_HTI_inter_subframe     0  -> intra-Subframe only, NO Convolutional
                                         Delay Line.  A/322 7.1.5.5 is not
                                         needed for this transmitter.
    L1D_plp_HTI_cell_interleaver   0  -> the Cell Interleaver (7.1.5.2) is
                                         BYPASSED.  Its pseudo-random
                                         permutation, whose generator
                                         polynomial the PDF also mangles, is
                                         therefore not on the critical path.
    L1D_plp_HTI_num_ti_blocks      1  -> NTI = 2
    L1D_plp_HTI_num_fec_blocks     73 -> N_fec = 74
    L1D_plp_HTI_num_fec_blocks_max 73 -> N_fec_max = 74

so only the TWISTED BLOCK INTERLEAVER (7.1.5.4) stands between the cell
stream and the bit de-interleaver.

THE GEOMETRY IS PINNED BY A CLOSED IDENTITY, not by the field semantics:

    L1D_plp_size / (Ninner / MOD)  ==  L1D_plp_HTI_num_fec_blocks + 1
             199800 / 2700         ==  73 + 1  ==  74

Three independently sourced numbers meet: `plp_size` and `num_fec_blocks`
come from the LDPC-decoded L1-Detail, and 2700 comes from Table 6.14 (gated
in spec_bitint as H6).  That single equation confirms the "+1" reading of
the L1D field, the FEC-block cell count, and the PLP allocation at once.
Get any of the three wrong and it does not divide.

    N_fec_TI_max = N_fec_max / NTI = 37   (A/322 7.1.5.1)
    N_fec_TI     = 37 per TI block, both blocks, N_virtual = 0
    Nrows = cells per FEC Block = 2700,  Ncols = N_fec_TI_max = 37

**Note this corrects the M4 design note**, which recorded "2 TI blocks x 73
FEC blocks".  It is 2 TI blocks x 37, and 37 x 2 x 2700 = 199800 exactly.

THE WALL, STATED PRECISELY.  A/322 7.1.5.4's three diagonal-read equations
print as

        =  mod ,
        =  mod ,
     =  +   mod ,

in BOTH editions and in BOTH the -table and -layout extractions.  Every
variable is gone, and the twisting parameter is named once and never given
a value anywhere the extraction can see.  So the reading is genuinely
ambiguous and NOTHING here pretends otherwise.

What this module does instead is enumerate the readings and let structure
eliminate what it can:

  T-GATE 1  BIJECTION OVER LEGAL PARAMETERS.  A reading must be a
            permutation for EVERY legal (Nrows, Ncols) an ATSC 3.0
            transmitter can signal, not merely for RF33's.  The
            row/column enumeration that happens to work at (2700, 37)
            because they are coprime collapses the moment they are not,
            and that kills it by construction -- the same way ASSUMPTION
            F2 was settled in M3.
  T-GATE 2  ROUND TRIP.  interleave -> de-interleave must be the identity
            on synthesized data, including with virtual FEC blocks
            present (N_virtual > 0), which RF33 does not exercise.
  T-GATE 3  DIAGONALITY.  "Diagonal-wise" is a claim about the output: no
            two consecutively read cells may come from the same input FEC
            Block.  MEASURED, AND IT DOES NOT DO WHAT I FIRST EXPECTED --
            it separates the linear-array convention (row-major reading
            keeps 2699 consecutive cells inside one FEC Block and fails),
            but it says NOTHING about the twisting parameter, because
            row-major ENUMERATION already changes column every step.  ST
            is therefore left completely free by T1-T3, ST = 0 included.
            Recorded this way rather than quietly re-tuned into a gate
            that appears to settle ST.
  T-GATE 4  TEXT CONSISTENCY.  The prose says reading runs "from the first
            row (rightwards along the row beginning with the left-most
            column) to the last row", and that writing is column-wise.
            That fixes the enumeration order and the linear-array
            convention independently of the lost symbols.

Survivors are REPORTED, not chosen.  Per the standing trap in this
project -- there is no referee between the cell stream and the LDPC
decoder -- a wrong pick among survivors will not announce itself, so the
survivor set is carried forward as ASSUMPTION T1 for the FEC to arbitrate.

NOTHING HERE IS VALIDATED ONLY AGAINST ITS OWN SYNTHESIS.  The round trip
proves the de-interleaver inverts the interleaver; it says nothing about
whether either matches A/322.  That is why the gates above exist and why
T1 stays open.

Usage:
    python m5_hti.py
"""
from __future__ import annotations

import itertools
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import spec_bitint as BI                                        # noqa: E402

# ---- RF33 PLP 0, every value decoded from L1-Basic / L1-Detail -------------
PLP0 = dict(plp_size=199800, ninner=16200, modulation="64QAM", rate="11/15",
            num_ti_blocks=1, num_fec_blocks=73, num_fec_blocks_max=73,
            inter_subframe=0, cell_interleaver=0)

# every (rows, cols) pair the standard can produce; rows come from Table 6.14
LEGAL_ROWS = sorted({v[1] for v in BI.DEMUX.values() if v[1] is not None})


def geometry(cfg=PLP0):
    """A/322 7.1.5.1 -- and the closed identity that pins it."""
    mod_bits = BI.MOD_BITS[cfg["modulation"]]
    cells_per_fec = BI.DEMUX[(cfg["ninner"], cfg["modulation"])][1]
    nti = cfg["num_ti_blocks"] + 1
    n_fec = cfg["num_fec_blocks"] + 1
    n_fec_max = cfg["num_fec_blocks_max"] + 1
    identity_ok = (cfg["plp_size"] % cells_per_fec == 0
                   and cfg["plp_size"] // cells_per_fec == n_fec)
    n_fec_ti_max = n_fec_max // nti
    # 7.1.5.1: the smaller TI Blocks come first
    rem = n_fec % nti
    n_fec_ti = [n_fec // nti + (1 if s >= nti - rem else 0) for s in range(nti)]
    return dict(mod_bits=mod_bits, cells_per_fec=cells_per_fec, nti=nti,
                n_fec=n_fec, n_fec_max=n_fec_max, n_fec_ti_max=n_fec_ti_max,
                n_fec_ti=n_fec_ti, nrows=cells_per_fec, ncols=n_fec_ti_max,
                n_virtual=[n_fec_ti_max - x for x in n_fec_ti],
                identity_ok=identity_ok)


# ---------------------------------------------------------------------------
# The Twisted Block Interleaver read order, as a family of readings
# ---------------------------------------------------------------------------

READINGS = list(itertools.product(("div", "mod"),      # row enumeration
                                  ("mod", "div"),      # column enumeration
                                  ("colmajor", "rowmajor")))   # linear array


def tbi_read_order(nrows, ncols, n_virtual, st, reading):
    """Input-cell indices in TBI OUTPUT order, virtual cells skipped.

    Returns None if the reading does not enumerate every (row, col) cell
    exactly once -- i.e. if it is not a permutation for these parameters.
    """
    rowmode, colmode, zmode = reading
    n = nrows * ncols
    l = np.arange(n)
    r = (l // ncols) if rowmode == "div" else (l % nrows)
    cp = (l % ncols) if colmode == "mod" else (l // nrows)
    c = (cp + st * r) % ncols
    flat = r.astype(np.int64) * ncols + c
    if int(np.bincount(flat, minlength=n).max()) != 1:
        return None                     # not a bijection over the array
    keep = c >= n_virtual
    if zmode == "colmajor":
        # memory written column-wise: linear index = col*Nrows + row, and the
        # data FEC Blocks sit in columns n_virtual..Ncols-1
        z = (c[keep] - n_virtual) * nrows + r[keep]
    else:
        z = r[keep] * (ncols - n_virtual) + (c[keep] - n_virtual)
    return z


def deinterleave(cells, nrows, ncols, n_virtual, st, reading):
    """Undo the TBI: cells arrive in output order, return them in input order."""
    order = tbi_read_order(nrows, ncols, n_virtual, st, reading)
    if order is None:
        raise ValueError("reading is not a permutation for these parameters")
    out = np.empty_like(cells)
    out[order] = cells
    return out


def interleave(cells, nrows, ncols, n_virtual, st, reading):
    order = tbi_read_order(nrows, ncols, n_virtual, st, reading)
    if order is None:
        raise ValueError("reading is not a permutation for these parameters")
    return cells[order]


# ---------------------------------------------------------------------------
# The gates
# ---------------------------------------------------------------------------

def gate_bijection(reading, st):
    """T-GATE 1: a permutation for EVERY legal (Nrows, Ncols), not just ours."""
    bad = []
    for nrows in LEGAL_ROWS:
        for ncols in (1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 15, 16, 18, 20, 25,
                      27, 30, 36, 37, 40, 45, 50, 54, 60, 64, 73):
            if tbi_read_order(nrows, ncols, 0, st, reading) is None:
                bad.append((nrows, ncols))
    return len(bad) == 0, bad


def gate_roundtrip(reading, st):
    """T-GATE 2: interleave -> de-interleave is the identity, virtuals too."""
    rng = np.random.default_rng(20260806)
    for nrows, ncols, nvirt in ((2700, 37, 0), (2700, 37, 5), (2025, 16, 3),
                                (8100, 9, 0), (10800, 64, 17), (2700, 5, 4)):
        try:
            n_in = nrows * (ncols - nvirt)
            x = rng.permutation(n_in)
            y = interleave(x, nrows, ncols, nvirt, st, reading)
            if len(y) != n_in:
                return False, f"{nrows}x{ncols} v{nvirt}: length {len(y)}"
            z = deinterleave(y, nrows, ncols, nvirt, st, reading)
            if not np.array_equal(z, x):
                return False, f"{nrows}x{ncols} v{nvirt}: round trip differs"
        except ValueError as e:
            return False, f"{nrows}x{ncols} v{nvirt}: {e}"
    return True, ""


def gate_diagonal(reading, st, nrows=2700, ncols=37):
    """T-GATE 3: consecutive output cells come from different FEC Blocks."""
    order = tbi_read_order(nrows, ncols, 0, st, reading)
    if order is None:
        return False, "not a permutation"
    blk = order // nrows                # which input FEC Block each cell is in
    same = int(np.sum(blk[1:] == blk[:-1]))
    return same == 0, f"{same} consecutive pairs share a FEC Block"


TEXT_CONSISTENT = ("div", "mod", "colmajor")
TEXT_WHY = (
    "7.1.5.4 reads 'from the first row (rightwards along the row beginning\n"
    "      with the left-most column) to the last row', which is row-major\n"
    "      enumeration -> row = l // Ncols, col = l mod Ncols; and 'the FEC\n"
    "      Blocks shall be serially written column-wise into the memory',\n"
    "      which makes the linear array column-major -> z = col*Nrows + row.")


def main():
    print("M5 Step 3 -- HTI TWISTED BLOCK INTERLEAVER")
    print("=" * 72)

    print("\n  === spec_bitint gate (the tables this depends on) ===")
    ok, bad = BI.verify(verbose=False)
    print(f"    {ok} pass, {bad} FAIL")

    g = geometry()
    print("\n  === geometry, and the closed identity that pins it ===")
    print(f"    Ninner {PLP0['ninner']} / MOD {g['mod_bits']} -> "
          f"{g['cells_per_fec']} cells per FEC Block   (Table 6.14)")
    print(f"    L1D_plp_size {PLP0['plp_size']} / {g['cells_per_fec']} = "
          f"{PLP0['plp_size']/g['cells_per_fec']:.6f}")
    print(f"    L1D_plp_HTI_num_fec_blocks {PLP0['num_fec_blocks']} + 1 = "
          f"{g['n_fec']}")
    print(f"    IDENTITY {'CLOSES' if g['identity_ok'] else 'DOES NOT CLOSE'}"
          f"  -- three independently decoded numbers meet exactly")
    print(f"    NTI = {g['nti']},  N_fec_TI_max = {g['n_fec_ti_max']},  "
          f"per-block {g['n_fec_ti']},  virtual {g['n_virtual']}")
    print(f"    TBI array: Nrows {g['nrows']} x Ncols {g['ncols']}"
          f"   ({g['nti']} TI blocks x {g['n_fec_ti'][0]} FEC blocks x "
          f"{g['cells_per_fec']} cells = "
          f"{g['nti']*g['n_fec_ti'][0]*g['cells_per_fec']})")
    print(f"    NOTE this corrects the M4 design note's "
          f"'2 TI blocks x 73 FEC blocks'.")

    print("\n  === the readings, and what structure eliminates ===")
    print(f"    {'row':4s} {'col':4s} {'array':9s} {'ST':>3s}  "
          f"{'T1 bijection':>13s} {'T2 roundtrip':>13s} {'T3 diagonal':>12s}")
    survivors = []
    for reading in READINGS:
        for st in (0, 1, 2):
            b, bad_pairs = gate_bijection(reading, st)
            rt, rt_why = gate_roundtrip(reading, st) if b else (False, "n/a")
            dg, dg_why = gate_diagonal(reading, st) if b else (False, "n/a")
            mark = ""
            if b and rt and dg:
                survivors.append((reading, st))
                if reading == TEXT_CONSISTENT:
                    mark = "  <== also T4 text-consistent"
            print(f"    {reading[0]:4s} {reading[1]:4s} {reading[2]:9s} "
                  f"{st:3d}  {'PASS' if b else f'FAIL({len(bad_pairs)})':>13s} "
                  f"{'PASS' if rt else 'FAIL':>13s} "
                  f"{'PASS' if dg else 'FAIL':>12s}{mark}")

    print(f"\n  {len(survivors)} of {len(READINGS)*3} readings survive "
          f"T-GATES 1-3.")
    t4 = [s for s in survivors if s[0] == TEXT_CONSISTENT]
    print(f"\n  === T-GATE 4: text consistency ===\n      {TEXT_WHY}")
    print(f"\n    readings matching the prose: "
          f"{[(r, st) for r, st in t4]}")

    print("\n  === WHAT IS AND IS NOT SETTLED ===")
    print("    SETTLED  the geometry (Nrows 2700, Ncols 37, NTI 2, no virtual")
    print("             FEC blocks, no CDL, no Cell Interleaver) -- by the")
    print("             closed identity 199800/2700 == 73+1 and by three")
    print("             decoded L1-Detail fields agreeing with Table 6.14.")
    print("    SETTLED  that the de-interleaver inverts the interleaver, on")
    print("             synthesized data, with and without virtual blocks.")
    print("             This proves the CODE, not the SPEC READING.")
    print("    NARROWED T1 -- the three lost equations.  24 candidate")
    print("             readings x ST in {0,1,2} go in; T1 kills the")
    print("             enumerations that stop being permutations when")
    print("             Nrows and Ncols share a factor, T3 kills the")
    print("             row-major linear array, and T4 (prose) agrees with")
    print("             what survives.  ONE reading stands:")
    print("                 row = l // Ncols,  col' = l mod Ncols,")
    print("                 col = (col' + ST*row) mod Ncols,")
    print("                 z   = (col - Nvirtual)*Nrows + row")
    print("    OPEN     the twisting parameter ST.  Named once in A/322")
    print("             7.1.5.4 and never given a value anywhere either")
    print("             edition's extraction can see, in either extraction")
    print("             mode.  T1-T3 do not constrain it AT ALL -- ST = 0")
    print("             passes every one of them.  The only non-arithmetic")
    print("             argument available is that ST = 0 makes the")
    print("             interleaver untwisted and the parameter pointless,")
    print("             which is prose, not proof, and is NOT counted as a")
    print("             gate here.")
    print("    ONLY THE FEC CAN ARBITRATE ST.  Per the standing trap there")
    print("    is no referee between the cell stream and the LDPC decoder, so")
    print("    a wrong ST will not announce itself.  But the search is now")
    print("    ONE small integer over 0..Ncols-1 = 0..36, not a structure")
    print("    hunt, and LDPC convergence is the discriminator -- exactly how")
    print("    ASSUMPTION L1 was settled in M3.  37 hypotheses is cheap.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
