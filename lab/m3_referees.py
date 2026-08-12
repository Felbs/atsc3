#!/usr/bin/env python3
"""M3 Step 3 -- separate the two referees, and give each one a real control.

The first sweep in m3_l1basic.py reported "BCH PASS" for every variant,
including variants whose LDPC left 800+ checks unsatisfied.  A result that
cannot fail is not evidence, so this script pins down what each referee is
actually testing, and adds the controls that were missing.

WHY BCH PASSED EVERYWHERE (and why that is good news, not a bug)
----------------------------------------------------------------
The 368 Nouter bits are SYSTEMATIC: they are transmitted directly as QPSK at
~20 dB MER, so their hard decisions are already essentially error-free before
the LDPC does anything.  Changing the LDPC parity permutation cannot move them.
So the BCH syndrome is not testing the LDPC at all -- it is testing everything
UPSTREAM of it:

    BCH zero syndrome  =>  QPSK demapping, the Annex C.1.1 bit-to-IQ
                           convention, the 6.5.2.10 block de-interleave, the
                           6.5.2.4 shortening pattern (WHICH codeword positions
                           carry the 368 bits), the frequency de-interleaver,
                           and the cell ordering are ALL correct.
                           A 168-bit check: 1 in 2^168 by chance.

    LDPC convergence   =>  the Annex A.2.2 parity-check table, the parity
                           interleave order, the 6.5.2.6 group-wise
                           permutation and the 6.5.2.8 puncturing schedule are
                           correct.  12960 simultaneous checks.

Two independent referees covering two disjoint halves of the chain is a better
outcome than one referee covering both -- but only once each is shown to be
capable of failing.  Hence the controls below.

CONTROLS
  C1  wrong shortening pattern  -> BCH must FAIL   (isolates 6.5.2.4)
  C2  swapped QPSK I/Q roles    -> BCH must FAIL   (isolates Annex C.1.1)
  C3  reverse FI wire direction -> BCH must FAIL   (isolates 7.3)
  C4  wrong LDPC variant        -> LDPC must FAIL  (already seen; repeated here)
  C5  random 368 bits           -> BCH must FAIL   (calibrates the referee)
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
import m3_l1basic as L1                                        # noqa: E402
import spec_bicm as BI                                         # noqa: E402
from m3_preamble import analyse                                 # noqa: E402


def nouter_from_cells(x, mode, shortening=None, swap_iq=False):
    """Hard-decision route: no LDPC at all, straight from the QPSK cells."""
    ncells = S.L1_BASIC_CELLS_PRINTED[mode]
    nfec, _, _ = S.l1_basic_lengths(mode)
    z = x[:ncells]
    z = z / np.sqrt(np.mean(np.abs(z) ** 2))
    if swap_iq:
        z = z.imag + 1j * z.real
    bl = L1.cells_to_llr(z, nfec=nfec, eta=S.L1_BASIC_ETA_MOD[mode])
    hard = (bl < 0).astype(np.uint8)
    # place into the 16200 frame and read the Nouter positions
    ipos, _ = L1.info_positions(pattern=shortening)
    return hard[:368], ipos


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("capture")
    ap.add_argument("--rate", type=float, required=True)
    ap.add_argument("--fmt")
    ap.add_argument("--json")
    a = ap.parse_args()
    path = a.capture
    if not os.path.isabs(path):
        path = os.path.join(HERE, path)

    rep = {}
    rep, z, Y, geo = analyse(path, a.rate, a.fmt, report=rep, quiet=True)
    (lo, n, pilot, cp, data, shift, cred, nfft, gi, dx, mode) = geo
    x = FI.deinterleave(z, nfft, 0, direction="forward", toggle="i")
    xrev = FI.deinterleave(z, nfft, 0, direction="reverse", toggle="i")
    print(f"  {os.path.basename(path)}  L1-Basic Mode {mode}  "
          f"pilot SNR {rep['pilot_snr_db']:.1f} dB\n")

    rng = np.random.default_rng(11)
    out = {}

    print("  === REFEREE 1: BCH syndrome over the 368 systematic bits ===")
    print("  (tests demap + bit de-interleave + shortening + FI + cell order;"
          " no LDPC)\n")
    rows = []

    def bch_case(label, bits, expect):
        syn = L1.bch_syndrome(bits)
        ok = syn == 0
        verdict = "PASS" if ok else "FAIL"
        good = (verdict == expect)
        rows.append({"case": label, "syndrome_zero": bool(ok),
                     "expected": expect, "as_expected": good})
        print(f"    {label:44s} {verdict:4s}  expected {expect:4s}   "
              f"{'ok' if good else '<-- UNEXPECTED'}")
        return ok

    b, ipos = nouter_from_cells(x, mode)
    bch_case("as-specified (A/322 chain)", b, "PASS")

    # C1 is NOT a valid BCH control and the first run proved it.  The 368
    # systematic bits are read from their codeword positions in ASCENDING
    # order under any shortening pattern, so the pattern permutes WHICH
    # positions are used but not the ORDER of the 368 bits -- the BCH syndrome
    # is therefore blind to it.  Worse, a pattern whose free groups are
    # disjoint from the true ones reads back all zeros, and the all-zero
    # vector is a legal BCH codeword, so the control "passed" for the wrong
    # reason entirely.  The shortening pattern is tested against the LDPC
    # below, where it does bite.

    # C2 -- swapped I/Q roles in the QPSK map
    b2, _ = nouter_from_cells(x, mode, swap_iq=True)
    bch_case("C2 QPSK I/Q roles swapped", b2, "FAIL")

    # C3 -- reverse FI wire direction
    b3, _ = nouter_from_cells(xrev, mode)
    bch_case("C3 reverse FI wire permutation", b3, "FAIL")

    # C5 -- random
    bch_case("C5 random 368 bits", rng.integers(0, 2, 368).astype(np.uint8),
             "FAIL")
    out["bch"] = rows

    print("\n  === REFEREE 2: LDPC convergence (12960 parity checks) ===")
    print("  (tests Annex A.2.2 + parity interleave + group-wise perm + "
          "puncturing)\n")
    lrows = []
    for variant in LD.PARITY_INTERLEAVE_VARIANTS:
        for gwi in (False, True):
            r = L1.decode_one(x, mode, variant=variant, gw_invert=gwi,
                              iters=100)
            lab = f"{variant}/{'inverted' if gwi else 'direct'}"
            expect = "PASS" if (variant == "standard" and not gwi) else "FAIL"
            got = "PASS" if r["converged"] else "FAIL"
            lrows.append({"variant": lab, "converged": r["converged"],
                          "iters": r["iters"], "unsat": r["unsatisfied"]})
            print(f"    {lab:44s} {got:4s}  iters {r['iters']:3d}  "
                  f"unsatisfied {r['unsatisfied']:4d}"
                  f"   {'ok' if got == expect else '<-- UNEXPECTED'}")
    # C1 -- wrong shortening pattern, tested where it actually matters
    for wrong in ([0, 1, 2, 3, 4, 5, 6, 7, 8],
                  [4, 1, 5, 2, 8, 6, 0, 3, 7]):
        import m3_l1basic as _L1
        _orig = BI.SHORTENING_PATTERN["L1-Basic"]
        try:
            BI.SHORTENING_PATTERN["L1-Basic"] = tuple(wrong)
            r = _L1.decode_one(x, mode, variant="standard", gw_invert=False,
                               iters=100)
        finally:
            BI.SHORTENING_PATTERN["L1-Basic"] = _orig
        got = "PASS" if r["converged"] else "FAIL"
        lrows.append({"variant": f"C1 shortening {wrong}",
                      "converged": r["converged"], "iters": r["iters"],
                      "unsat": r["unsatisfied"]})
        print(f"    C1 wrong shortening {str(wrong[:5]):18s}         "
              f"{got:4s}  iters {r['iters']:3d}  "
              f"unsatisfied {r['unsatisfied']:4d}"
              f"   {'ok' if got == 'FAIL' else '<-- UNEXPECTED'}")

    out["ldpc"] = lrows

    # do the LDPC-corrected bits equal the raw hard decisions?
    r = L1.decode_one(x, mode, variant="standard", gw_invert=False, iters=100)
    agree = int((r["nouter_bits"] == b).sum())
    print(f"\n  LDPC-decoded Nouter vs raw hard decisions: {agree}/368 agree "
          f"({368-agree} bits corrected)")
    print("  -> confirms why BCH passed for every LDPC variant: the systematic"
          "\n     bits arrive clean at ~20 dB MER, so BCH never sees the LDPC.")
    out["ldpc_vs_hard_agree"] = agree

    dest = a.json or os.path.join(HERE, "m3_referees_"
                                  + os.path.splitext(os.path.basename(path))[0]
                                  + ".json")
    with open(dest, "w") as fh:
        json.dump(out, fh, indent=2, default=float)
    print(f"\n  wrote {dest}")


if __name__ == "__main__":
    main()
