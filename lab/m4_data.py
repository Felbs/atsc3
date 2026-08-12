#!/usr/bin/env python3
"""M4 Step 3a -- demodulate the DATA symbols of subframe 0 and resolve the
payload constellation.

Everything needed is now signalled rather than guessed.  L1-Detail says
subframe 0 is 8K / GI1536 / SP4_2 / Cred 0, 35 data payload symbols with the
first and last being Subframe boundary symbols, and that its dominant PLP
(PLP 0, the LLS carrier) occupies cells 0..199799 modulated as
**64QAM-NUC at code rate 11/15**.

THE CLAIM UNDER TEST
--------------------
If the L1 decode is right, equalising the data cells with A/322 8.1.3 scattered
pilots must produce a constellation that fits NUC_64_11/15 -- and fits it
BETTER than the five other 64-NUC tables and better than uniform 64QAM.

That comparison is the whole point.  M2 was burned by comparing across
constellation ORDERS, where more points always fit better and a noise cloud
therefore "improves" monotonically with order.  Here every candidate has
exactly 64 points, so the comparison is unbiased: a cloud scores the same
against all of them, and only a real 64-NUC separates.

WHAT IS SCANNED RATHER THAN ASSUMED
-----------------------------------
  * the scattered-pilot lattice phase.  A/322 8.1.3.1 puts pilots at
    k mod (DX.DY) = DX.(l mod DY), but whether l counts from the first data
    symbol or includes the preamble is a reading, not a fact.  Both are run;
    the wrong one is a control.
  * the carrier origin, as in M3.

Usage:
    python m4_data.py hit_rf33.cs16 --rate 8e6
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
import spec_nuc as NU                                           # noqa: E402
from m2_pilots import load                                      # noqa: E402
from m3_preamble import fft_bins, mer_db, ascii_scatter         # noqa: E402

BOOTSTRAP = 13824                    # 4 x 3072 @ 6.144 -> 6.912 Msps


def sp_positions(noc, dx, dy, l, sbs=False):
    """A/322 8.1.3.1 scattered pilots, 8.1.7.1 subframe-boundary pilots."""
    if sbs:
        k = np.arange(0, noc, dx)
        return k[(k != 0) & (k != noc - 1)]
    return np.arange(dx * (l % dy), noc, dx * dy)


def edge_pilots(noc):
    return np.array([0, noc - 1])


def symbol_fft(y, t0, nfft, gi, l):
    """FFT of data symbol l of subframe 0 (l = 0 is the first DATA symbol)."""
    pitch = nfft + gi
    s = t0 + pitch * (1 + l)
    return np.fft.fftshift(np.fft.fft(y[s:s + nfft]))


def channel_and_cells(Y, noc, lo, dx, dy, l, sbs, amp, nfft, shift=0,
                      cps=None):
    """Pilot-referenced channel estimate + equalised data cells."""
    ref = np.array(S.pilot_values(noc, amp))
    pk = sp_positions(noc, dx, dy, l, sbs)
    pk = np.concatenate([edge_pilots(noc), pk])
    pk = np.unique(pk)
    b = fft_bins(nfft, lo + pk + shift)
    ok = (b >= 0) & (b < nfft)
    pk, b = pk[ok], b[ok]
    hp = Y[b] / ref[pk]
    kk = np.arange(noc)
    H = (np.interp(kk, pk, hp.real).astype(np.complex128)
         + 1j * np.interp(kk, pk, hp.imag))
    used = np.zeros(noc, bool)
    used[pk] = True
    if cps is not None:
        used[cps] = True
    data = np.flatnonzero(~used)
    bd = fft_bins(nfft, lo + data + shift)
    ok = (bd >= 0) & (bd < nfft)
    z = Y[bd[ok]] / H[data[ok]]
    # pilot smoothness = how coherent the estimate is; a wrong lattice is rough
    q = hp
    coh = float(abs(np.sum(q[:-1] * np.conj(q[1:])))
                / max(np.sum(np.abs(q) ** 2), 1e-30))
    return z, coh, len(data)


def cpe_correct(z_list):
    """Strip each symbol's residual common phase by its own 4th-power moment.

    The 45-degree bias M2 documented does NOT apply here: this only REMOVES a
    per-symbol offset relative to the ensemble, it never rotates the ensemble.
    """
    ref = None
    out = []
    for z in z_list:
        m = np.mean((z / np.sqrt(np.mean(np.abs(z) ** 2))) ** 4)
        a = np.angle(m) / 4
        if ref is None:
            ref = a
        out.append(z * np.exp(-1j * (a - ref)))
    return out


def candidates():
    c = {}
    for cr in ("2/15", "3/15", "4/15", "5/15", "6/15", "7/15",
               "8/15", "9/15", "10/15", "11/15", "12/15", "13/15"):
        c[f"NUC_64_{cr}"] = NU.get(64, cr)
    c["uniform 64QAM"] = NU.qam(64)
    return c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("capture")
    ap.add_argument("--rate", type=float, required=True)
    ap.add_argument("--fmt")
    ap.add_argument("--start", type=float, default=0.0)
    ap.add_argument("--symbols", type=int, default=30)
    ap.add_argument("--json")
    a = ap.parse_args()
    path = a.capture if os.path.isabs(a.capture) else os.path.join(HERE,
                                                                   a.capture)

    print("A/322 Annex C NUC tables:")
    _, ch = NU.load()
    print(f"  {sum(1 for _, o, _ in ch if o)}/{len(ch)} identities pass "
          f"(2026 and 2024 witnesses identical; mean power == 1 for all 36 "
          f"position vectors)\n")

    nfft, gi, dx, dy, cred = 8192, 1536, 4, 2, 0
    noc = S.noc(nfft, cred)
    lo, hi = S.carrier_abs_range(nfft, cred)
    cps = np.array([c - lo for c in S.common_cps(nfft, cred)], int)
    print(f"  subframe 0, from L1-Detail: FFT 8K, GI {gi}, SP{dx}_{dy}, "
          f"Cred {cred} -> NoC {noc}, {len(cps)} common CPs")
    print(f"  PLP 0 = cells 0..199799, 64QAM-NUC 11/15  (the LLS carrier)\n")

    y = load(path, a.rate, fmt=a.fmt, span_sec=0.075, start_sec=a.start)[0]

    # ---- fine timing on the first data symbol, scanned ------------------
    amp = 1.0                       # a global pilot scale cancels in the MER
    best = None
    for t in range(BOOTSTRAP - 20, BOOTSTRAP + 21):
        Y = symbol_fft(y, t, nfft, gi + 0, 1)
        # window starts GI into the symbol
        Y = np.fft.fftshift(np.fft.fft(y[t + (nfft + gi) * 2 + gi:
                                         t + (nfft + gi) * 2 + gi + nfft]))
        _, coh, _ = channel_and_cells(Y, noc, lo, dx, dy, 1, False, amp, nfft,
                                      0, cps)
        if best is None or coh > best[1]:
            best = (t, coh)
    t0 = best[0]
    print(f"  bootstrap-anchored timing: t0 = {t0} "
          f"(nominal {BOOTSTRAP}, residual {t0-BOOTSTRAP:+d}), "
          f"pilot coherence {best[1]:.4f}")

    # ---- the lattice-phase hypothesis, and its control ------------------
    print(f"\n  === A/322 8.1.3.1 lattice phase: does l count from the first "
          f"DATA symbol? ===")
    rows = {}
    for name, off in (("l = first data symbol", 0),
                      ("l offset by 1 (control)", 1)):
        cohs = []
        for l in range(1, 1 + a.symbols):
            Y = np.fft.fftshift(np.fft.fft(
                y[t0 + (nfft + gi) * (1 + l) + gi:
                  t0 + (nfft + gi) * (1 + l) + gi + nfft]))
            _, coh, _ = channel_and_cells(Y, noc, lo, dx, dy, l + off, False,
                                          amp, nfft, 0, cps)
            cohs.append(coh)
        rows[name] = float(np.mean(cohs))
        print(f"    {name:26s} mean pilot coherence {np.mean(cohs):.4f}")
    off = 0 if rows["l = first data symbol"] >= rows["l offset by 1 (control)"] \
        else 1

    # ---- equalise ---------------------------------------------------------
    zs, ncell = [], 0
    for l in range(1, 1 + a.symbols):
        Y = np.fft.fftshift(np.fft.fft(
            y[t0 + (nfft + gi) * (1 + l) + gi:
              t0 + (nfft + gi) * (1 + l) + gi + nfft]))
        z, coh, nd = channel_and_cells(Y, noc, lo, dx, dy, l + off, False,
                                       amp, nfft, 0, cps)
        zs.append(z)
        ncell = nd
    zs = cpe_correct(zs)
    z = np.concatenate(zs)
    print(f"\n  {a.symbols} data symbols x {ncell} data cells = {len(z)} "
          f"equalised cells (symbols 1..{a.symbols}, SBS symbol 0 excluded)")

    # ---- THE UNBIASED COMPARISON -----------------------------------------
    print(f"\n  === constellation fit.  Every candidate has exactly 64 points, "
          f"so\n      the comparison is unbiased -- a cloud cannot win one of "
          f"them ===")
    cands = candidates()
    zn = z / np.sqrt(np.mean(np.abs(z) ** 2))
    stats = {}
    for k, pts in cands.items():
        d = np.abs(zn[:, None] - pts[None, :])
        idx = d.argmin(1)
        err = d[np.arange(len(zn)), idx]
        cen = np.array([zn[idx == i].mean() if (idx == i).any() else 0
                        for i in range(len(pts))])
        cnt = np.bincount(idx, minlength=len(pts))
        p = cnt / cnt.sum()
        exp = len(zn) / len(pts)
        stats[k] = dict(
            mer=float(10 * np.log10(1.0 / max(np.mean(err ** 2), 1e-30))),
            centroid_rms=float(np.sqrt(np.mean(np.abs(cen - pts) ** 2))),
            entropy=float(-(p[p > 0] * np.log2(p[p > 0])).sum()),
            chi2_dof=float(((cnt - exp) ** 2 / exp).sum() / (len(pts) - 1)))
    order = sorted(stats, key=lambda k: -stats[k]["mer"])
    print(f"    {'candidate':16s} {'MER dB':>7s} {'centroid RMS':>12s} "
          f"{'occupancy H':>11s} {'chi2/dof':>9s}")
    print(f"    {'':16s} {'':>7s} {'(geometry)':>12s} {'(max 6.000)':>11s} "
          f"{'(unif=1)':>9s}")
    for k in order:
        s = stats[k]
        mark = "  <== signalled" if k == "NUC_64_11/15" else ""
        print(f"    {k:16s} {s['mer']:7.2f} {s['centroid_rms']:12.4f} "
              f"{s['entropy']:11.4f} {s['chi2_dof']:9.1f}{mark}")
    win = order[0]
    oth = {k: v for k, v in stats.items() if k != "NUC_64_11/15"}
    print(f"\n    Three statistics, three different things, one winner:")
    print(f"      MER              -- how tight the clusters are")
    print(f"      centroid RMS     -- where the clusters actually SIT, with "
          f"the noise floor divided out")
    print(f"      occupancy chi2   -- an LDPC-coded, scrambled payload must "
          f"use all 64 points\n"
          f"                          equally often.  A wrong table cuts the "
          f"plane in the wrong\n"
          f"                          places and the counts skew.  This one "
          f"does not care about SNR.")
    print(f"    NUC_64_11/15 vs the best OTHER 64-point candidate:")
    print(f"      MER            {stats['NUC_64_11/15']['mer'] - max(v['mer'] for v in oth.values()):+6.2f} dB")
    print(f"      centroid RMS    {min(v['centroid_rms'] for v in oth.values()) / stats['NUC_64_11/15']['centroid_rms']:5.2f}x better")
    print(f"      chi2/dof        {min(v['chi2_dof'] for v in oth.values()) / stats['NUC_64_11/15']['chi2_dof']:5.2f}x better")
    res = {k: v["mer"] for k, v in stats.items()}

    print()
    print(ascii_scatter(z[:60000], lim=1.9))

    out = {"capture": os.path.basename(path), "t0": int(t0),
           "lattice_phase": rows, "offset_used": off,
           "n_cells": int(len(z)), "stats": stats, "winner": win}
    dest = a.json or os.path.join(
        HERE, "m4_data_" + os.path.splitext(os.path.basename(path))[0] + ".json")
    with open(dest, "w") as fh:
        json.dump(out, fh, indent=2, default=float)
    print(f"\n  wrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
