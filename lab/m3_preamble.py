#!/usr/bin/env python3
"""M3 Step 2 -- demodulate the ATSC 3.0 PREAMBLE using the real A/322 tables.

No radio.  Offline only.  Operates on captures already banked in lab/.

THE CLAIM UNDER TEST
--------------------
M2 ended with an unresolved preamble cloud: best-fit MER rose monotonically
with constellation order (QPSK 4.4 -> 256QAM 15.9 dB), which is the signature
of noise, not of a constellation.  M2 diagnosed the cause as missing published
tables rather than a broken algorithm.  If that diagnosis is right, feeding the
real A/322 tables in must turn the cloud into clean QPSK.  If it does not, the
diagnosis was wrong and we have learned something more valuable than a pass.

WHAT M2 HAD WRONG, per A/322
----------------------------
1.  **NoC.**  A/322 7.2.5.1: "For the first Preamble symbol, the NoC used shall
    be the minimum for the FFT size used for that first Preamble symbol."  The
    minimum for 8K is 6529 (Cred_coeff = 4), NOT the 6913 of the data symbols.
    M2 extracted cells across the full ~6913-carrier span, so it swept in 384
    carriers that the transmitter never modulated and put every relative
    carrier index off by 192.
2.  **Pilot pattern.**  A/322 8.1.6.1: Preamble pilots always use DY = 1 and
    sit at every relative k with k mod DX == 0 -- a quarter of all carriers for
    RF33's DX = 4.  M2 applied the DATA symbol's SP4_2 scattered pattern
    (DY = 2), so three quarters of the true pilots were treated as data and a
    half of the assumed pilot positions were not pilots.
3.  **Pilot values.**  M2 recovered pilot signs blind.  A/322 8.1.2 publishes
    the reference PRBS, so the pilots are now KNOWN, exactly, including sign --
    which also lets M2's blind recovery be scored rather than trusted.

Everything above is a table lookup, and every table is cross-verified in
m3_spec.py (13/13 identities).

DISCIPLINE
----------
The carrier origin is not assumed: it is SCANNED, and wrong origins are the
built-in controls.  The Cred_coeff is also scanned, so "the spec says minimum
NoC" is a hypothesis that has to beat four alternatives on the data rather than
being taken on faith.  A statistic that cannot fail is not evidence.

Usage:
    python m3_preamble.py hit_rf33.cs16 --rate 8e6
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

from atsc3 import bootstrap as bs                              # noqa: E402
from m2_pilots import load                                     # noqa: E402
import m3_spec as S                                            # noqa: E402

BOOTSTRAP_SAMPLES_POST = 13824      # 4 x 3072 @ 6.144 Msps, scaled to 6.912


# --------------------------------------------------------------------------
# cell geometry
# --------------------------------------------------------------------------

def preamble_geometry(nfft, gi, dx, cred):
    """Return (abs0, n_carriers, pilot_rel, cp_rel, data_rel).

    abs0      absolute carrier index of relative carrier 0
    pilot_rel relative indices of preamble pilots (k mod dx == 0)
    cp_rel    relative indices of transmitted common continual pilots
    data_rel  relative indices of available data cells, in ascending order
    """
    lo, hi = S.carrier_abs_range(nfft, cred)
    n = S.noc(nfft, cred)
    assert hi - lo + 1 == n
    pilot = np.arange(0, n, dx)
    cp = np.array([c - lo for c in S.common_cps(nfft, cred)], dtype=int)
    used = np.zeros(n, bool)
    used[pilot] = True
    overlap = int(used[cp].sum()) if len(cp) else 0
    used[cp] = True
    data = np.flatnonzero(~used)
    return lo, n, pilot, cp, data, overlap


def fft_bins(nfft, abs_idx):
    """Absolute carrier index -> index into a fftshift()ed FFT of size nfft.

    ASSUMPTION S1: NoCmax is odd and its centre carrier is DC.
    """
    centre = (S.NOC_MAX[nfft] - 1) // 2
    return (nfft // 2) + (np.asarray(abs_idx) - centre)


# --------------------------------------------------------------------------
# pilot-referenced channel estimate  (pilot VALUES now known from A/322 8.1.2)
# --------------------------------------------------------------------------

def pilot_coherence(Y, nfft, gi, dx, cred, shift=0):
    """Score a carrier-origin hypothesis WITHOUT fitting anything.

    At preamble pilots the transmitted value is a known real +-A.  Dividing the
    received bin by it leaves A*H(k), which is SMOOTH in k because the channel
    is.  If the origin (hence the sign pattern r_k) is wrong, the quotient
    carries random sign flips and is ROUGH.  So

        C = |sum_j Q_j conj(Q_{j+1})| / sum_j |Q_j|^2

    over consecutive pilots is ~1 for the right origin and ~0 otherwise.  It
    fits no parameters, so it cannot be inflated by searching.
    """
    lo, n, pilot, _cp, _d, _o = preamble_geometry(nfft, gi, dx, cred)
    amp = S.PREAMBLE_PILOT_BOOST[(nfft, gi)][1]
    ref = np.array(S.pilot_values(n, amp))
    b = fft_bins(nfft, lo + pilot + shift)
    ok = (b >= 0) & (b < nfft)
    b, pk = b[ok], pilot[ok]
    q = Y[b] / ref[pk]
    num = np.abs(np.sum(q[:-1] * np.conj(q[1:])))
    den = np.sum(np.abs(q) ** 2)
    return float(num / max(den, 1e-30))


def channel_estimate(Y, nfft, gi, dx, cred, shift=0):
    """H over all relative carriers, from the known preamble pilots."""
    lo, n, pilot, cp, data, _o = preamble_geometry(nfft, gi, dx, cred)
    amp = S.PREAMBLE_PILOT_BOOST[(nfft, gi)][1]
    ref = np.array(S.pilot_values(n, amp))
    b = fft_bins(nfft, lo + pilot + shift)
    hp = Y[b] / ref[pilot]
    # linear interpolation on the pilot lattice, real and imaginary separately
    kk = np.arange(n)
    H = (np.interp(kk, pilot, hp.real).astype(np.complex128)
         + 1j * np.interp(kk, pilot, hp.imag))
    return H, hp


# --------------------------------------------------------------------------
# constellation metrics
# --------------------------------------------------------------------------

QPSK = np.array([1 + 1j, 1 - 1j, -1 + 1j, -1 - 1j]) / np.sqrt(2)


def mer_db(z, points):
    """MER against the nearest ideal point, after unit-power normalisation."""
    z = z / np.sqrt(np.mean(np.abs(z) ** 2))
    d = np.abs(z[:, None] - points[None, :])
    err = d.min(axis=1)
    return float(10 * np.log10(1.0 / max(np.mean(err ** 2), 1e-30)))


def m4_stat(z):
    """E[z^4]/E[|z|^4].  For QPSK on the diagonals this is REAL and NEGATIVE
    with magnitude ~1; for a Gaussian cloud it is ~0.  A cloud cannot fake it.
    """
    z = z / np.sqrt(np.mean(np.abs(z) ** 2))
    v = np.mean(z ** 4) / max(np.mean(np.abs(z) ** 4), 1e-30)
    return complex(v)


def qam_points(m):
    r = int(round(np.sqrt(m)))
    lv = np.arange(-(r - 1), r, 2)
    pts = np.array([complex(i, q) for i in lv for q in lv])
    return pts / np.sqrt(np.mean(np.abs(pts) ** 2))


def ascii_scatter(z, w=57, h=25, lim=1.9):
    z = z / np.sqrt(np.mean(np.abs(z) ** 2))
    g = np.zeros((h, w), int)
    xi = ((z.real + lim) / (2 * lim) * (w - 1)).round().astype(int)
    yi = ((lim - z.imag) / (2 * lim) * (h - 1)).round().astype(int)
    m = (xi >= 0) & (xi < w) & (yi >= 0) & (yi < h)
    np.add.at(g, (yi[m], xi[m]), 1)
    ramp = " .:-=+*#%@"
    mx = g.max() or 1
    return "\n".join("".join(ramp[min(len(ramp) - 1, int(v / mx * (len(ramp) - 1) + 0.5))]
                             for v in row) for row in g)


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def analyse(path, rate, fmt=None, shift_span=40, report=None, start_sec=0.0,
            quiet=False):
    rep = report if report is not None else {}
    _p = print if not quiet else (lambda *a, **k: None)
    y, fs_post, cfo, hit = load(path, rate, fmt=fmt, span_sec=0.10,
                                start_sec=start_sec)
    rep["start_sec"] = start_sec
    f = hit["fields"]
    ps = f["preamble_structure"]
    rep["capture"] = os.path.basename(path)
    rep["preamble_structure"] = ps
    rep["fine_cfo_hz"] = cfo
    rep["fs_post"] = fs_post

    if ps not in S.PREAMBLE_STRUCTURE:
        raise SystemExit(f"preamble_structure {ps} not in A/322 Table H.1.1")
    P = S.PREAMBLE_STRUCTURE[ps]
    nfft, gi, dx, mode = P["fft"], P["gi"], P["dx"], P["l1b_mode"]
    rep["A322_TableH11"] = dict(P)
    _p(f"\n  bootstrap says preamble_structure = {ps}")
    _p(f"  A/322 Table H.1.1 -> FFT {nfft//1024}K, GI {gi}, preamble pilot "
          f"DX {dx}, L1-Basic Mode {mode}")
    _p(f"  A/322 Table 6.17  -> {S.L1_BASIC_CONSTELLATION[mode]}, "
          f"{S.L1_BASIC_CELLS_PRINTED[mode]} cells, Ksig 200, rate 3/15")

    # ---- fine timing: the FFT window starts GI samples into the symbol -----
    base = BOOTSTRAP_SAMPLES_POST + gi
    best = None
    for t in range(base - 24, base + 25):
        if t + nfft > len(y):
            continue
        Y = np.fft.fftshift(np.fft.fft(y[t:t + nfft]))
        c = pilot_coherence(Y, nfft, gi, dx, 4, 0)
        if best is None or c > best[1]:
            best = (t, c)
    t0 = best[0]
    rep["fft_window_start"] = t0
    rep["fft_window_nominal"] = base
    _p(f"\n  FFT window start {t0} (nominal {base}, "
          f"residual {t0 - base:+d} samples)")

    Y = np.fft.fftshift(np.fft.fft(y[t0:t0 + nfft]))

    # ---- scan carrier origin AND Cred_coeff.  Wrong values are controls. ---
    _p(f"\n  Scanning carrier origin and Cred_coeff "
          f"(A/322 7.2.5.1 predicts Cred_coeff = 4, i.e. NoC = "
          f"{S.noc(nfft, 4)}):")
    scan = {}
    for cred in range(5):
        row = []
        for sh in range(-shift_span, shift_span + 1):
            row.append((sh, pilot_coherence(Y, nfft, gi, dx, cred, sh)))
        bsh, bc = max(row, key=lambda r: r[1])
        scan[cred] = {"best_shift": bsh, "best_coherence": bc,
                      "noc": S.noc(nfft, cred),
                      "median_coherence": float(np.median([c for _, c in row]))}
        _p(f"    Cred_coeff {cred}  NoC {S.noc(nfft, cred):5d}   "
              f"best coherence {bc:.4f} @ shift {bsh:+d}   "
              f"(median over scan {scan[cred]['median_coherence']:.4f})")
    rep["origin_scan"] = scan

    cred = max(scan, key=lambda c: scan[c]["best_coherence"])
    shift = scan[cred]["best_shift"]
    rep["cred_chosen"] = cred
    rep["shift_chosen"] = shift
    spec_says = 4
    _p(f"\n  WINNER: Cred_coeff {cred}, shift {shift:+d}, "
          f"coherence {scan[cred]['best_coherence']:.4f}")
    _p(f"  A/322 7.2.5.1 predicted Cred_coeff {spec_says} "
          f"({'AGREES' if cred == spec_says else 'DISAGREES'})")

    # ---- equalize ----------------------------------------------------------
    lo, n, pilot, cp, data, overlap = preamble_geometry(nfft, gi, dx, cred)
    want_cells = S.PREAMBLE_DATA_CELLS[(nfft, gi)][cred]
    _p(f"\n  geometry: NoC {n}, preamble pilots {len(pilot)}, "
          f"common CPs {len(cp)} (overlap with pilots {overlap}), "
          f"data cells {len(data)}")
    _p(f"  A/322 Table 7.2 says {want_cells} data cells -> "
          f"{'MATCH' if len(data) == want_cells else 'MISMATCH'}")
    rep["n_data_cells"] = int(len(data))
    rep["n_data_cells_table72"] = want_cells

    H, hp = channel_estimate(Y, nfft, gi, dx, cred, shift)
    b = fft_bins(nfft, lo + data + shift)
    ok = (b >= 0) & (b < nfft)
    z = Y[b[ok]] / H[data[ok]]
    rep["pilot_coherence"] = scan[cred]["best_coherence"]

    # pilot-referenced SNR: scatter of H along the lattice vs its smooth part
    hs = np.convolve(hp, np.ones(9) / 9, mode="same")
    resid = hp[4:-4] - hs[4:-4]
    psnr = 10 * np.log10(np.mean(np.abs(hs[4:-4]) ** 2)
                         / max(np.mean(np.abs(resid) ** 2), 1e-30))
    rep["pilot_snr_db"] = float(psnr)
    _p(f"  pilot-referenced SNR {psnr:.1f} dB")

    # ---- THE VERDICT -------------------------------------------------------
    _p(f"\n  --- constellation of all {len(z)} preamble data cells ---")
    orders = [("QPSK", QPSK), ("16QAM", qam_points(16)),
              ("64QAM", qam_points(64)), ("256QAM", qam_points(256))]
    mers = {}
    for name, pts in orders:
        mers[name] = mer_db(z, pts)
    rep["mer_all_cells"] = mers
    for name, _ in orders:
        _p(f"    best-fit MER vs {name:6s}  {mers[name]:6.2f} dB")

    m4 = m4_stat(z)
    rep["m4_all_cells"] = {"abs": abs(m4), "deg": float(np.degrees(np.angle(m4)))}
    _p(f"    E[z^4]/E[|z|^4] = {abs(m4):.3f} at "
          f"{np.degrees(np.angle(m4)):+.1f} deg "
          f"(QPSK on the diagonals -> magnitude ~1 at 180 deg; "
          f"noise -> magnitude ~0)")

    _p()
    _p(ascii_scatter(z))

    # first N cells = L1-Basic (in CELL order, i.e. after de-interleaving)
    return rep, z, Y, (lo, n, pilot, cp, data, shift, cred, nfft, gi, dx, mode)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("capture")
    ap.add_argument("--rate", type=float, required=True)
    ap.add_argument("--fmt")
    ap.add_argument("--json")
    a = ap.parse_args()

    print("A/322 spec-table verification first:")
    r = S.verify()
    if not all(g for _, g, _ in r):
        raise SystemExit("spec tables FAILED verification -- stopping")

    path = a.capture
    if not os.path.isabs(path):
        path = os.path.join(HERE, path)
    rep, z, Y, geo = analyse(path, a.rate, a.fmt)

    out = a.json or os.path.join(
        HERE, "m3_preamble_" + os.path.splitext(os.path.basename(path))[0] + ".json")
    with open(out, "w") as fh:
        json.dump(rep, fh, indent=2, default=float)
    print(f"\n  wrote {out}")


if __name__ == "__main__":
    main()
