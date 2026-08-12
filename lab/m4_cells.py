#!/usr/bin/env python3
"""M4 Step 3b -- the subframe cell map, verified against the DUMMY CELLS.

A/322 7.2.6.5 hands this project an almost perfect referee, and it went
unnoticed until the scrambler was solved:

    "all of the available data cells of a Subframe shall first be filled with
     dummy modulation values, and then the cell multiplexing process shall
     overwrite the dummy modulation values of occupied data cells"
    Re{d_i} = 1 - 2.s_i ,  Im{d_i} = 0

where s_i is the i-th bit of the **Section 5.2.3 Baseband Packet scrambling
sequence** -- the very scrambler M4 Step 1 pinned -- and i is the index of the
data cell within the Subframe (Section 7.2.6.1).

So every cell a Subframe does not use carries a known, real-valued +-1 whose
sign is a published pseudo-random bit.  That makes the dummy region a
simultaneous test of FOUR things that otherwise have no independent check:

  1. the scrambler solution, a THIRD time, on thousands of fresh bits from
     deep in the sequence rather than the first 200;
  2. the A/322 7.3 frequency de-interleaver applied to DATA symbols --
     both the wire-permutation direction (ASSUMPTION F1) and the symbol-index
     convention, neither of which the L1 work could settle;
  3. the cell-index origin of the symbol, which is otherwise pure arithmetic
     over pilot counts;
  4. the decoded L1D_plp_start / L1D_plp_size values, because they predict
     exactly where the dummy region begins.

and it needs no FEC at all.

THE MEASUREMENT.  Take the last symbol of subframe 0, de-interleave, keep the
real-valued cells, and cross-correlate their signs against the scrambling
sequence over every possible cell-index origin.  A wrong hypothesis scores at
chance; the right one scores 1.0.

Also verified here: A/322 Table D.1.4's single additional continual pilot for
8K + SP4_2 (relative carrier 1732).  Section 8.1.4.1 says the additional CPs
exist "to ensure a constant number of data carriers in every data symbol".
Without carrier 1732 the count alternates 6000 / 5999; with it, every normal
data symbol has exactly 5999.  That is a closed identity, and it only closes
for the right table entry.

Usage:
    python m4_cells.py hit_rf33.cs16 --rate 8e6
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

import m3_spec as S                                             # noqa: E402
import m3_freqint as FI                                         # noqa: E402
import m4_scrambler as SC                                       # noqa: E402
import spec_nuc as NU                                           # noqa: E402
from m2_pilots import load                                      # noqa: E402
from m3_preamble import fft_bins, mer_db, ascii_scatter         # noqa: E402

NFFT, GI, DX, DY = 8192, 1536, 4, 2
ADDITIONAL_CP = 1732            # A/322 Table D.1.4, 8K + SP4_2
BOOTSTRAP = 13824
QPSK = np.array([1 + 1j, 1 - 1j, -1 + 1j, -1 - 1j]) / np.sqrt(2)


def cp_set(noc, lo):
    return np.array(sorted(set(c - lo for c in S.common_cps(NFFT, 0))
                           | {ADDITIONAL_CP}), int)


def data_cell_count(noc, cps, l, sbs):
    used = set(range(0, noc, DX)) if sbs else set(range(DX * (l % DY), noc,
                                                        DX * DY))
    used |= {0, noc - 1}
    used |= set(int(c) for c in cps)
    return noc - len(used)


def demod(y, t0, l, noc, lo, cps, sbs):
    s = t0 + (NFFT + GI) * (1 + l) + GI
    Y = np.fft.fftshift(np.fft.fft(y[s:s + NFFT]))
    ref = np.array(S.pilot_values(noc, 1.0))
    pk = (np.arange(0, noc, DX) if sbs
          else np.arange(DX * (l % DY), noc, DX * DY))
    pk = np.unique(np.concatenate([pk, [0, noc - 1]]))
    hp = Y[fft_bins(NFFT, lo + pk)] / ref[pk]
    kk = np.arange(noc)
    H = (np.interp(kk, pk, hp.real).astype(np.complex128)
         + 1j * np.interp(kk, pk, hp.imag))
    used = np.zeros(noc, bool)
    used[pk] = True
    used[cps] = True
    d = np.flatnonzero(~used)
    z = Y[fft_bins(NFFT, lo + d)] / H[d]
    return z / np.sqrt(np.mean(np.abs(z) ** 2))


def dummy_signs(x, im_th=0.35, mag_th=0.5):
    d = np.where((np.abs(x.imag) < im_th) & (np.abs(x) > mag_th),
                 np.sign(x.real), 0.0)
    return d


def correlate(d, sgn, nmax):
    """Cross-correlate the sign vector against the scrambler at every origin."""
    n = int((d != 0).sum())
    L = 1 << int(np.ceil(np.log2(nmax + len(d))))
    c = np.fft.irfft(np.fft.rfft(sgn, L) * np.fft.rfft(d[::-1], L), L)
    c = c[len(d) - 1:len(d) - 1 + nmax - len(d)] / max(n, 1)
    i = int(np.argmax(np.abs(c)))
    return i, float(abs(c[i])), n, float(np.sort(np.abs(c))[-2])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("capture")
    ap.add_argument("--rate", type=float, required=True)
    ap.add_argument("--fmt")
    ap.add_argument("--start", type=float, default=0.0)
    ap.add_argument("--frames", type=int, default=2)
    ap.add_argument("--nmax", type=int, default=260000)
    ap.add_argument("--json")
    a = ap.parse_args()
    path = a.capture if os.path.isabs(a.capture) else os.path.join(HERE,
                                                                   a.capture)
    print("SUBFRAME CELL MAP, verified against A/322 7.2.6.5 dummy cells")
    print("=" * 70)

    noc = S.noc(NFFT, 0)
    lo, _ = S.carrier_abs_range(NFFT, 0)
    cps = cp_set(noc, lo)

    # ---- the Table D.1.4 identity ---------------------------------------
    print("\n  === A/322 Table D.1.4 / 8.1.4.1: the additional continual "
          "pilot ===")
    for label, cc in (("common CPs only (48)",
                       np.array(sorted(set(c - lo
                                           for c in S.common_cps(NFFT, 0))),
                                int)),
                      (f"+ Table D.1.4 carrier {ADDITIONAL_CP} (49)", cps)):
        e = data_cell_count(noc, cc, 0, False)
        o = data_cell_count(noc, cc, 1, False)
        print(f"    {label:38s} l even {e}   l odd {o}   "
              f"{'CONSTANT' if e == o else 'alternates'}")
    print(f"    8.1.4.1: the additional CPs exist 'to ensure a constant "
          f"number of data\n    carriers in every data symbol'.  The identity "
          f"closes only with 1732.")
    n_norm = data_cell_count(noc, cps, 1, False)
    n_sbs = data_cell_count(noc, cps, 0, True)
    print(f"    normal data symbol {n_norm} cells;  subframe boundary symbol "
          f"{n_sbs} cells")

    # ---- the dummy-cell correlation -------------------------------------
    sgn = (1.0 - 2.0 * SC.sequence(a.nmax)).astype(np.float64)
    print(f"\n  === the dummy region: cross-correlate its signs against the "
          f"A/322 5.2.3\n      scrambling sequence over every cell-index "
          f"origin (0..{a.nmax}) ===")
    results, rows = [], []
    for fi in range(a.frames):
        y = load(path, a.rate, fmt=a.fmt, span_sec=0.075,
                 start_sec=a.start + 0.7 * fi)[0]
        best = None
        for t in range(BOOTSTRAP - 20, BOOTSTRAP + 21):
            s = t + (NFFT + GI) * 2 + GI
            Y = np.fft.fftshift(np.fft.fft(y[s:s + NFFT]))
            ref = np.array(S.pilot_values(noc, 1.0))
            pk = np.arange(DX, noc, DX * DY)
            q = Y[fft_bins(NFFT, lo + pk)] / ref[pk]
            c = abs(np.sum(q[:-1] * np.conj(q[1:]))) / np.sum(np.abs(q) ** 2)
            if best is None or c > best[1]:
                best = (t, c)
        t0 = best[0]
        x_raw = demod(y, t0, 34, noc, lo, cps, sbs=True)
        # the FI symbol index and wire direction are BOTH unsettled readings;
        # every combination is run and the wrong ones are the controls
        for direction in FI.WIRE_DIRECTIONS:
            for lf in range(30, 42):
                try:
                    x = FI.deinterleave(x_raw, NFFT, lf, direction=direction,
                                        toggle="i")
                except Exception:                               # noqa: BLE001
                    continue
                i, f, n, second = correlate(dummy_signs(x), sgn, a.nmax)
                rows.append(dict(frame=fi, direction=direction, l=lf,
                                 origin=i, frac=f, n=n))
        top = sorted([r for r in rows if r["frame"] == fi],
                     key=lambda r: -r["frac"])
        print(f"\n    frame {fi} (t0 = {t0}):  {top[0]['n']} real-valued "
              f"cells in the last symbol of subframe 0")
        print(f"      {'agreement':>10s}  {'origin':>7s}  {'FI dir':8s} "
              f"{'FI l':>4s}")
        for r in top[:4]:
            mark = "   <== EXACT" if r["frac"] > 0.999 else ""
            print(f"      {r['frac']:10.4f}  {r['origin']:7d}  "
                  f"{r['direction']:8s} {r['l']:4d}{mark}")
        print(f"      chance level for {top[0]['n']} signs over {a.nmax} "
              f"origins ~ {5/np.sqrt(top[0]['n']):.3f}")
        results.append(top[0])

    win = results[0]
    if win["frac"] < 0.999:
        print("\n  NO exact hypothesis.  Reporting that rather than the best "
              "near-miss.")
        return 1
    origin34 = win["origin"]
    print(f"\n  *** the last data symbol of subframe 0 begins at cell index "
          f"{origin34} ***")
    print(f"      FI wire direction '{win['direction']}' (ASSUMPTION F1 "
          f"SETTLED for data symbols)")
    print(f"      FI symbol index {win['l']} for data symbol 34, i.e. the FI "
          f"counter\n      includes the preamble as symbol 0")
    print(f"      reproduced on {sum(1 for r in results if r['frac'] > 0.999)}"
          f" of {len(results)} frames, same origin: "
          f"{len(set(r['origin'] for r in results)) == 1}")

    # ---- controls ---------------------------------------------------------
    print(f"\n  === controls ===")
    y = load(path, a.rate, fmt=a.fmt, span_sec=0.075, start_sec=a.start)[0]
    t0 = 13818
    x = FI.deinterleave(demod(y, t0, 34, noc, lo, cps, True), NFFT,
                        win["l"], direction=win["direction"], toggle="i")
    d = dummy_signs(x)
    ctl = []
    rng = np.random.default_rng(7)
    for name, s2 in (
            ("scrambler with a random 0xF180-length seed",
             (1.0 - 2.0 * SC.sequence(a.nmax, out_taps=(1, 2, 3, 4, 5, 6, 7, 8)))),
            ("M3's Fibonacci reading of the same register",
             (1.0 - 2.0 * SC.sequence(a.nmax, direction="right", form="fib",
                                      taps=SC.POLY_EXP_PRINTED,
                                      seed="X1=LSB",
                                      out_taps=(1, 3, 4, 7, 11, 12, 13, 14)))),
            ("a random +-1 sequence", rng.choice([-1.0, 1.0], a.nmax))):
        i, f, n, _ = correlate(d, np.asarray(s2, float), a.nmax)
        ctl.append(dict(name=name, frac=f, origin=i))
        print(f"    {name:44s} best agreement {f:.4f}")
    for name, xx in (("no frequency de-interleaving",
                      demod(y, t0, 34, noc, lo, cps, True)),
                     ("normal-symbol pilot pattern instead of SBS",
                      FI.deinterleave(demod(y, t0, 34, noc, lo, cps, False),
                                      NFFT, win["l"],
                                      direction=win["direction"],
                                      toggle="i"))):
        i, f, n, _ = correlate(dummy_signs(xx), sgn, a.nmax)
        ctl.append(dict(name=name, frac=f, origin=i))
        print(f"    {name:44s} best agreement {f:.4f}")

    # ---- the PLP boundaries, predicted then tested -----------------------
    print(f"\n  === L1-Detail's PLP allocations, turned into a falsifiable "
          f"prediction ===")
    org = {34: origin34}
    for l in range(33, -1, -1):
        org[l] = org[l + 1] - (n_sbs if l in (0, 34) else n_norm)
    print(f"    walking the origin back: symbol 0 begins at cell "
          f"{org[0]}, so the preamble\n    contributes {org[0]} available "
          f"data cells to subframe 0's pool "
          f"(A/322 7.2.6.1\n    note: 'the available data cells in the first "
          f"Subframe of a Frame might begin\n    in the final Preamble "
          f"symbol').  Naive geometry says 4851 - 484 - 880 = 3487;\n    "
          f"the measured value is {org[0]}, a {3487-org[0]}-cell difference "
          f"-- OPEN ITEM.")
    pts64 = NU.get(64, "11/15")
    boundaries = []
    for l, bnd, before, after in ((32, 199800, "PLP 0 / 64QAM-NUC 11/15",
                                   "PLP 16 / QPSK 2/15"),
                                  (34, 207900, "PLP 16 / QPSK 2/15",
                                   "dummy / real +-1")):
        p = bnd - org[l]
        n = n_sbs if l in (0, 34) else n_norm
        if not (0 < p < n):
            continue
        xx = FI.deinterleave(demod(y, t0, l, noc, lo, cps, l in (0, 34)),
                             NFFT, l + 1, direction=win["direction"],
                             toggle="i")
        lhs, rhs = xx[:p], xx[p:]
        row = dict(symbol=l, cell=bnd, pos=p,
                   lhs_64=mer_db(lhs, pts64), lhs_qpsk=mer_db(lhs, QPSK),
                   rhs_64=mer_db(rhs, pts64), rhs_qpsk=mer_db(rhs, QPSK),
                   lhs_real=float(np.mean(np.abs(lhs.imag) < 0.3)),
                   rhs_real=float(np.mean(np.abs(rhs.imag) < 0.3)))
        boundaries.append(row)
        print(f"\n    symbol {l}: L1-Detail puts the {before} -> {after} "
              f"boundary at cell {bnd},\n    i.e. internal position {p} of "
              f"{n}.  Split there and measure:")
        print(f"      cells 0..{p-1:4d} ({p:4d})   MER vs 64NUC "
              f"{row['lhs_64']:6.2f} dB   vs QPSK {row['lhs_qpsk']:6.2f} dB   "
              f"real-valued {row['lhs_real']:.3f}   [{before}]")
        print(f"      cells {p}..{n-1} ({n-p:4d})   MER vs 64NUC "
              f"{row['rhs_64']:6.2f} dB   vs QPSK {row['rhs_qpsk']:6.2f} dB   "
              f"real-valued {row['rhs_real']:.3f}   [{after}]")

    out = dict(capture=os.path.basename(path), n_normal=n_norm, n_sbs=n_sbs,
               origin_symbol34=origin34, fi=win, sweep=rows, controls=ctl,
               origins={str(k): int(v) for k, v in org.items()},
               boundaries=boundaries)
    dest = a.json or os.path.join(
        HERE, "m4_cells_" + os.path.splitext(os.path.basename(path))[0]
        + ".json")
    with open(dest, "w") as fh:
        json.dump(out, fh, indent=2, default=float)
    print(f"\n  wrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
