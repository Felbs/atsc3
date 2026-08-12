#!/usr/bin/env python3
"""M5 Step 1 -- the subframe cell map, CLOSED.

M4 left two unexplained numbers:

    subframe 0   the preamble appeared to contribute 3297 available data
                 cells where its spare-cell geometry says 3487   -> 190 off
    subframe 1   the measured cell origin was 946712 where the pilot-count
                 model said 949179                               -> 2467 off

M4 diagnosed both as "too few continual pilots".  That diagnosis was
PARTLY right and mostly wrong, and the correction matters:

  (a) Table D.1.4 really was mis-read -- but it costs only ONE cell per
      symbol in three of every four 16K/SP4_4 symbols (54 of the 2467),
      and NOTHING at all in subframe 0.  See `spec_pilots.py`: the SPx_4
      row groups are three rows tall with a vertically-centred label, so
      the naive reading hands SPx_4's first additional CP to SPx_2.
  (b) The real defect is that the model used the SBS symbol's TOTAL data
      cell count (A/322 Table 7.5/7.6) where cell multiplexing uses the
      ACTIVE count (Annex F).  The difference is L1D_sbs_null_cells.
  (c) And the dummy-cell correlation measures the origin of the SBS
      symbol's FIRST TOTAL DATA CELL, while the available-cell index
      starts N_null_low cells later -- the null cells sit at the two band
      edges (A/322 7.2.6.4, Figure 7.14).

Put together:  190 = 127 + 63  and  2467 = 1609 + 804 + 54.
Every term is now a named quantity, not a residual.

THE GATES, all of which must hold simultaneously:

  A. `spec_pilots.verify()` -- 2420 closed arithmetic identities over every
     legal (FFT, Cred, SP-pattern), with the naive D.1.4 reading kept as a
     control that must FAIL.
  B. Table 7.5/7.6 minus Annex F must reproduce the two L1D_sbs_null_cells
     values the transmitter actually signalled: 127 and 1609.  Those come
     off the air; nothing here can fit them.
  C. The dummy-cell correlation must stay at ~1.0 -- a wrong pilot set
     changes which carriers are data cells and destroys it.
  D. The walked-back cell budget must land on the preamble's spare-cell
     count for subframe 0 (3487, itself = 4851 - 484 - 880 from three
     independently decoded quantities) and on ZERO for subframe 1, which
     inherits nothing from a preamble.

Usage:
    python m5_cellmap.py hit_rf33.cs16 --rate 8e6
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
import spec_pilots as P                                         # noqa: E402
from m2_pilots import load                                      # noqa: E402
from m3_preamble import fft_bins                                # noqa: E402

BOOTSTRAP = 13824
SCRAMBLER_PERIOD = 524280               # 65535 states x 8 output bits

# ---- the RF33 configuration, every field DECODED from L1-Basic/L1-Detail ----
# (see m4_l1detail_hit_rf33.json; nothing below is a guess)
SUBFRAMES = [
    dict(name="subframe 0", nfft=8192, gi=1536, pattern="SP4_2", cred=0,
         nsym=35, spb=1, sbs_null=127, preamble_syms=1,
         # A/322 7.2.5.2: the preamble's cells after L1-Basic and L1-Detail
         # are available data cells of the first Subframe.
         preamble_data_cells=4851, l1b_cells=484, l1d_cells=880),
    dict(name="subframe 1", nfft=16384, gi=1536, pattern="SP4_4", cred=0,
         nsym=75, spb=1, sbs_null=1609, preamble_syms=0,
         preamble_data_cells=0, l1b_cells=0, l1d_cells=0),
]


def subframe_start(sf_index):
    """Sample offset of a Subframe's first symbol, relative to t0."""
    off = 0
    for k in range(sf_index):
        c = SUBFRAMES[k]
        off += (c["nfft"] + c["gi"]) * (c["nsym"] + c["preamble_syms"])
    return off


def pilot_k(cfg, l, sbs):
    """Pilot-bearing relative carrier indices, from spec_pilots."""
    return np.array(sorted(P.pilot_carriers(cfg["nfft"], cfg["cred"],
                                            cfg["pattern"], l, sbs)), int)


def demod(y, t0, cfg, sf_index, l, sbs, amp=1.0):
    """Equalised data cells of symbol l of a Subframe, in CARRIER order."""
    nfft, gi, cred = cfg["nfft"], cfg["gi"], cfg["cred"]
    noc = P.NOC[(nfft, cred)]
    lo, _ = S.carrier_abs_range(nfft, cred)
    s = t0 + subframe_start(sf_index) + (nfft + gi) * (l + cfg["preamble_syms"]) + gi
    Y = np.fft.fftshift(np.fft.fft(y[s:s + nfft]))
    ref = np.array(S.pilot_values(noc, amp))
    dx, dy = P.dxdy(cfg["pattern"])
    pk = (np.arange(0, noc, dx) if sbs
          else np.arange(dx * (l % dy), noc, dx * dy))
    pk = np.unique(np.concatenate([pk, [0, noc - 1]]))
    hp = Y[fft_bins(nfft, lo + pk)] / ref[pk]
    kk = np.arange(noc)
    H = (np.interp(kk, pk, hp.real).astype(np.complex128)
         + 1j * np.interp(kk, pk, hp.imag))
    used = np.zeros(noc, bool)
    used[pilot_k(cfg, l, sbs)] = True
    d = np.flatnonzero(~used)
    z = Y[fft_bins(nfft, lo + d)] / H[d]
    coh = float(abs(np.sum(hp[:-1] * np.conj(hp[1:])))
                / max(np.sum(np.abs(hp) ** 2), 1e-30))
    return z / np.sqrt(np.mean(np.abs(z) ** 2)), coh


def dummy_signs(x, im_th=0.35, mag_th=0.5):
    return np.where((np.abs(x.imag) < im_th) & (np.abs(x) > mag_th),
                    np.sign(x.real), 0.0)


def correlate_period(d, sgn):
    """Circular cross-correlation of the sign vector over one full period.

    `sgn` must be one period plus len(d) wrapped samples.  Returns
    (best_origin_mod_period, agreement, n_signs, runner_up).
    """
    n = int((d != 0).sum())
    nmax = SCRAMBLER_PERIOD
    L = 1 << int(np.ceil(np.log2(len(sgn) + len(d))))
    c = np.fft.irfft(np.fft.rfft(sgn, L) * np.fft.rfft(d[::-1], L), L)
    c = c[len(d) - 1:len(d) - 1 + nmax] / max(n, 1)
    i = int(np.argmax(np.abs(c)))
    return i, float(abs(c[i])), n, float(np.sort(np.abs(c))[-2])


def model(cfg, table=None, use_total=False, null_low="floor"):
    """The available-data-cell model for one Subframe.

    Returns (n_norm, n_sbs_avail, n_null, n_low, cells_before_last_symbol).
    `table`/`use_total`/`null_low` exist so the wrong readings can be run
    as controls.
    """
    nfft, cred, pat, spb = (cfg["nfft"], cfg["cred"], cfg["pattern"],
                            cfg["spb"])
    dx, dy = P.dxdy(pat)
    counts = {P.data_cells(nfft, cred, pat, l, False, table) for l in range(dy)}
    # a control reading can leave this non-constant; take the model's own
    # per-phase counts rather than collapsing them
    tot, act, n_null = P.sbs_cells(nfft, cred, pat, spb)
    n_sbs = tot if use_total else act
    lo_n, hi_n = P.sbs_null_split(n_null)
    if null_low == "ceil":
        lo_n = hi_n
    elif null_low == "zero":
        lo_n = 0
    # symbols 0 .. nsym-2 precede the final (SBS) symbol; symbol 0 is SBS too
    before = n_sbs
    for l in range(1, cfg["nsym"] - 1):
        before += P.data_cells(nfft, cred, pat, l, False, table)
    return counts, n_sbs, n_null, lo_n, before


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("capture")
    ap.add_argument("--rate", type=float, required=True)
    ap.add_argument("--fmt")
    ap.add_argument("--start", type=float, default=0.0)
    ap.add_argument("--frames", type=int, default=2)
    ap.add_argument("--json")
    a = ap.parse_args()
    path = a.capture if os.path.isabs(a.capture) else os.path.join(HERE,
                                                                   a.capture)

    print("M5 Step 1 -- SUBFRAME CELL MAP, closed against A/322 7.2.6.4/7.2.6.5")
    print("=" * 72)

    # ---- GATE A ----------------------------------------------------------
    print("\n  === GATE A: spec_pilots closed identities ===")
    ok, bad = P.verify(verbose=False)
    ok2, bad2 = P.verify(verbose=False, table=P.ADDITIONAL_CP_NAIVE)
    print(f"    Table D.1.4 as re-extracted        {ok:5d} pass  {bad:4d} FAIL")
    print(f"    naive 2-row grouping (CONTROL)     {ok2:5d} pass  {bad2:4d} FAIL"
          f"   <- must fail")
    gateA = (bad == 0 and bad2 > 0)
    print(f"    GATE A {'PASS' if gateA else 'FAIL'}")

    # ---- GATE B ----------------------------------------------------------
    print("\n  === GATE B: Table 7.5/7.6 - Annex F must reproduce the "
          "SIGNALLED\n      L1D_sbs_null_cells (these came off the air) ===")
    gateB = True
    for cfg in SUBFRAMES:
        tot, act, nn = P.sbs_cells(cfg["nfft"], cfg["cred"], cfg["pattern"],
                                   cfg["spb"])
        good = nn == cfg["sbs_null"]
        gateB &= good
        print(f"    {cfg['name']}  {cfg['nfft']//1024}K/{cfg['pattern']}/"
              f"Cred{cfg['cred']}/SPB{cfg['spb']}:  total {tot} - active {act}"
              f" = {nn}   signalled {cfg['sbs_null']}   "
              f"{'MATCH' if good else 'MISMATCH'}")
    print(f"    GATE B {'PASS' if gateB else 'FAIL'}")

    # ---- the measurement --------------------------------------------------
    print("\n  === the dummy-cell correlation, one full scrambler period "
          f"(0..{SCRAMBLER_PERIOD-1}) ===")
    sgn_all = (1.0 - 2.0 * SC.sequence(SCRAMBLER_PERIOD + 20000)
               ).astype(np.float64)

    results, gateC, gateD = [], True, True
    for fi in range(a.frames):
        y = load(path, a.rate, fmt=a.fmt, span_sec=0.26,
                 start_sec=a.start + 0.30 * fi)[0]
        # fine timing, scanned on subframe 0's first normal data symbol
        best = None
        for t in range(BOOTSTRAP - 20, BOOTSTRAP + 21):
            _, coh = demod(y, t, SUBFRAMES[0], 0, 1, False)
            if best is None or coh > best[1]:
                best = (t, coh)
        t0 = best[0]
        print(f"\n    ---- frame {fi}: t0 = {t0} (pilot coherence "
              f"{best[1]:.4f}) ----")

        for si, cfg in enumerate(SUBFRAMES):
            last = cfg["nsym"] - 1
            x_raw, coh = demod(y, t0, cfg, si, last, True)
            x = FI.deinterleave(x_raw, cfg["nfft"], last + cfg["preamble_syms"],
                                direction="forward", toggle="i")
            d = dummy_signs(x)
            org_mod, frac, nsig, second = correlate_period(d, sgn_all)
            counts, n_sbs, n_null, lo_n, before = model(cfg)

            # absolute origin: the model predicts which period we are in
            k = round((before - lo_n - org_mod) / SCRAMBLER_PERIOD)
            org_abs = org_mod + k * SCRAMBLER_PERIOD
            # the correlation locates the FIRST TOTAL data cell; the first
            # AVAILABLE (active) cell is lo_n further in
            first_avail = org_abs + lo_n
            residual = first_avail - before
            preamble_contrib = residual
            expect = (cfg["preamble_data_cells"] - cfg["l1b_cells"]
                      - cfg["l1d_cells"])

            print(f"\n    {cfg['name']}: {cfg['nfft']//1024}K "
                  f"{cfg['pattern']} Cred{cfg['cred']}, {cfg['nsym']} symbols,"
                  f" SBS first+last")
            print(f"      data cells per NORMAL symbol   "
                  f"{sorted(counts)}  (Table 7.3 says "
                  f"{P.avail_data_cells(cfg['nfft'], cfg['cred'], cfg['pattern'])})")
            print(f"      SBS total / ACTIVE / null      "
                  f"{P.SBS_TOTAL[(cfg['nfft'],cfg['cred'],cfg['pattern'])][0]}"
                  f" / {n_sbs} / {n_null}   (null split {lo_n} low, "
                  f"{n_null-lo_n} high)")
            print(f"      dummy signs in last symbol     {nsig}")
            print(f"      correlation agreement          {frac:.4f}   "
                  f"(runner-up {second:.4f}, chance ~{5/np.sqrt(max(nsig,1)):.3f})")
            print(f"      origin mod {SCRAMBLER_PERIOD}          {org_mod}"
                  f"   -> absolute {org_abs}  (period k={k})")
            print(f"      first AVAILABLE cell           {org_abs} + {lo_n} "
                  f"= {first_avail}")
            print(f"      model: cells in symbols 0..{last-1}  {before}")
            print(f"      => cells inherited from the preamble = "
                  f"{preamble_contrib}   EXPECTED {expect}   "
                  f"**{'ZERO DISCREPANCY' if preamble_contrib == expect else f'OFF BY {preamble_contrib-expect}'}**")

            gateC &= frac > 0.95
            gateD &= (preamble_contrib == expect)
            results.append(dict(frame=fi, subframe=cfg["name"], t0=int(t0),
                                agreement=frac, n_signs=nsig,
                                origin_mod=org_mod, origin_abs=org_abs,
                                first_avail=first_avail, model_before=before,
                                preamble_contrib=preamble_contrib,
                                expected=expect,
                                discrepancy=preamble_contrib - expect,
                                n_sbs_active=n_sbs, n_null=n_null,
                                null_low=lo_n, norm_counts=sorted(counts)))

    # ---- controls: every wrong reading, run as a control -------------------
    print("\n  === CONTROLS -- each wrong reading, and what it costs ===")
    y = load(path, a.rate, fmt=a.fmt, span_sec=0.26, start_sec=a.start)[0]
    t0 = results[0]["t0"]
    ctl = []
    for label, kw in (("M4's reading (SBS TOTAL, no null offset, naive D.1.4)",
                       dict(table=P.ADDITIONAL_CP_NAIVE, use_total=True,
                            null_low="zero")),
                      ("SBS TOTAL instead of Annex F ACTIVE",
                       dict(use_total=True)),
                      ("no null-cell offset on the measured origin",
                       dict(null_low="zero")),
                      ("null split CEIL at the low edge",
                       dict(null_low="ceil")),
                      ("naive D.1.4 grouping only",
                       dict(table=P.ADDITIONAL_CP_NAIVE))):
        row = {"reading": label, "discrepancy": {}}
        for si, cfg in enumerate(SUBFRAMES):
            r = next(r for r in results if r["subframe"] == cfg["name"])
            counts, n_sbs, n_null, lo_n, before = model(cfg, **kw)
            expect = (cfg["preamble_data_cells"] - cfg["l1b_cells"]
                      - cfg["l1d_cells"])
            k = round((before - lo_n - r["origin_mod"]) / SCRAMBLER_PERIOD)
            got = r["origin_mod"] + k * SCRAMBLER_PERIOD + lo_n - before
            row["discrepancy"][cfg["name"]] = got - expect
        ctl.append(row)
        print(f"    {label:52s} sf0 {row['discrepancy']['subframe 0']:+7d}   "
              f"sf1 {row['discrepancy']['subframe 1']:+7d}")
    print(f"    {'THE CORRECTED READING':52s} sf0 "
          f"{results[0]['discrepancy']:+7d}   sf1 {results[1]['discrepancy']:+7d}")

    print("\n  === VERDICT ===")
    for g, name in ((gateA, "A  spec_pilots closed identities (+ control)"),
                    (gateB, "B  Annex F reproduces the signalled null counts"),
                    (gateC, "C  dummy-cell correlation held"),
                    (gateD, "D  cell budget closes on both subframes")):
        print(f"    GATE {name:52s} {'PASS' if g else 'FAIL'}")
    allg = gateA and gateB and gateC and gateD
    print(f"\n    {'ALL GATES PASS -- the cell map is closed.' if allg else 'NOT ALL GATES PASS.'}")

    out = dict(capture=os.path.basename(path), gates=dict(A=gateA, B=gateB,
                                                          C=gateC, D=gateD),
               results=results, controls=ctl)
    dest = a.json or os.path.join(
        HERE, "m5_cellmap_" + os.path.splitext(os.path.basename(path))[0]
        + ".json")
    with open(dest, "w") as fh:
        json.dump(out, fh, indent=2, default=float)
    print(f"\n  wrote {dest}")
    return 0 if allg else 1


if __name__ == "__main__":
    raise SystemExit(main())
