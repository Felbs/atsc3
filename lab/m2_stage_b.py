#!/usr/bin/env python3
"""M2 Stage B -- MEASURE the ATSC 3.0 OFDM grid from the signal itself.

No radio.  Offline only.  This stage deliberately reads NO spec bits: it takes
the samples that follow a bootstrap and measures

    * the FFT size          (cyclic-prefix autocorrelation lag)
    * the guard interval    (the symbol pitch that the CP peaks repeat at)
    * symbol duration, subcarrier spacing, occupied bandwidth, carrier count
    * the scattered-pilot spacing D_x (spectral line in the average |X(k)|^2)

and only THEN compares the result to what Stage A's bootstrap bits claim.
Agreement is mutual validation between two paths that share no code.
Disagreement localizes a bug.  Both numbers are printed either way.

Method notes (why it is built this way)
---------------------------------------
* The grid search folds the CP correlation profile modulo a candidate symbol
  pitch T = N_fft + GI and scores peak/median of the fold.  Peak-PICKING on the
  raw profile was tried first and is far too noisy -- it reported GI = 1521,
  which is on no menu.  Folding averages ~100 symbols and is unambiguous.
* CONTROL LAGS that are not legal ATSC 3.0 FFT sizes are searched alongside the
  legal ones.  If a control wins, the measurement is broken and says so.
* The pilot metric is a spectral LINE test on the average power spectrum, not a
  "best comb phase" ratio.  Best-phase ratios are biased upward in proportion to
  D_x for ANY heavy-tailed spectrum (the first version of this file fell into
  exactly that trap and "found" D_x = 32 in pure payload noise).

Usage:
    python m2_stage_b.py hit_rf33.cs16 --rate 8e6
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from fractions import Fraction

import numpy as np
from scipy.signal import resample_poly

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from atsc3 import bootstrap as bs                              # noqa: E402
from atsc3.capture import guess_format, read_block             # noqa: E402
from m2_stage_a import find_bootstraps                         # noqa: E402

# ASSUMPTION B1: the A/322 guard-interval menu, in samples at the post-bootstrap
# sample rate.  Used to ENUMERATE hypotheses and to LABEL the winner.  Control
# GIs that are on no menu are searched too, so a menu-only answer is not assumed.
GI_MENU = {192: "GI1_192", 384: "GI2_384", 512: "GI3_512", 768: "GI4_768",
           1024: "GI5_1024", 1536: "GI6_1536", 2048: "GI7_2048",
           2432: "GI8_2432", 3072: "GI9_3072", 3648: "GI10_3648",
           4096: "GI11_4096", 4864: "GI12_4864"}
GI_CONTROLS = (256, 640, 896, 1280, 1792, 2560, 3328, 4352)
# ASSUMPTION B2: A/322 allows only these three FFT sizes.
FFT_MENU = (8192, 16384, 32768)
FFT_CONTROLS = (4096, 6144, 12288, 24576)


def resample_to(x, fs_in, fs_out):
    if abs(fs_in - fs_out) < 1.0:
        return np.asarray(x, dtype=np.complex64)
    fr = Fraction(int(round(fs_out)), int(round(fs_in))).limit_denominator(20000)
    return resample_poly(x, fr.numerator, fr.denominator).astype(np.complex64)


def _sliding_sum(x, w):
    c = np.concatenate(([np.zeros((), x.dtype)], np.cumsum(x)))
    return c[w:] - c[:-w]


def cp_profile(x, lag, win):
    """Normalized delay-and-correlate profile at a given lag.

    For the true FFT size the cyclic prefix makes x[n] == x[n+N] over GI
    samples once per OFDM symbol, so this profile spikes once per symbol.
    For any other lag it is noise.  No ATSC 3.0 knowledge is used.
    """
    x = np.asarray(x, dtype=np.complex128)
    n = len(x) - lag
    if n <= win:
        return np.zeros(0)
    prod = x[:n] * np.conj(x[lag:lag + n])
    pw = np.abs(x[:n]) ** 2 + np.abs(x[lag:lag + n]) ** 2
    s = _sliding_sum(prod, win)
    e = _sliding_sum(pw, win)
    return 2.0 * np.abs(s) / np.maximum(np.real(e), 1e-30)


def fold_score(p, period):
    """Fold the profile modulo `period` and score peak/median of the fold."""
    n = (len(p) // period) * period
    if n < 3 * period:
        return 0.0, 0, 0
    f = p[:n].reshape(-1, period).mean(axis=0)
    med = float(np.median(f))
    return float(f.max() / max(med, 1e-30)), int(np.argmax(f)), n // period


def grid_search(y, ffts, gis, min_symbols=4):
    rows = []
    for nfft in ffts:
        for gi in gis:
            p = cp_profile(y, nfft, gi)
            if len(p) == 0:
                continue
            sc, ph, nsym = fold_score(p, nfft + gi)
            if nsym < min_symbols:
                continue
            rows.append({"fft": nfft, "gi": gi, "period": nfft + gi,
                         "score": sc, "phase": ph, "n_folded": nsym,
                         "legal": (nfft in FFT_MENU and gi in GI_MENU)})
    rows.sort(key=lambda r: -r["score"])
    return rows


# Legal A/322 scattered-pilot D_x values, plus CONTROL periods that ATSC 3.0
# never uses.  If a control scores like the legal ones, there is no comb and the
# "detection" is noise -- this is the null the first version of this file lacked.
DX_LEGAL = (3, 4, 6, 8, 12, 16, 24, 32)
DX_CONTROL = (5, 7, 9, 11, 13, 17)
DX_ALL = tuple(sorted(set(DX_LEGAL) | set(DX_CONTROL) |
                      {d * 2 for d in DX_LEGAL} | {d * 4 for d in DX_LEGAL}))


def _detrend(power, k=None):
    """Divide out the smooth channel/filter shape so only the comb survives."""
    if k is None:
        k = max(65, (len(power) // 64) | 1)
    base = np.convolve(power, np.ones(k) / k, mode="same")
    edge = k // 2 + 1
    r = power[edge:-edge] / np.maximum(base[edge:-edge], 1e-30)
    return r / max(r.mean(), 1e-30) - 1.0


def spectral_line(power, periods=DX_ALL):
    """Unbiased comb test: look for a LINE in the DFT of the power spectrum.

    A comb of boosted pilots every D_x carriers puts energy at DFT index
    len/D_x.  Score = |line| / median(|DFT|), which has no D_x-dependent bias.
    The spectrum is de-trended first, because the smooth channel response and
    the band edges otherwise dominate the DFT and bury the line.
    """
    q = _detrend(np.asarray(power, dtype=np.float64))
    Q = np.abs(np.fft.rfft(q * np.hanning(len(q))))
    base = float(np.median(Q[1:]))
    out = []
    for d in periods:
        idx = int(round(len(q) / d))
        if not 1 <= idx < len(Q) - 1:
            continue
        v = float(Q[max(1, idx - 2): idx + 3].max())
        out.append({"dx": d, "line_ratio": v / max(base, 1e-30), "bin": idx,
                    "legal": d in DX_LEGAL})
    out.sort(key=lambda d: -d["line_ratio"])
    return out


def per_symbol_line(syms, lo, hi, periods=None):
    """Comb test averaged over symbols in MAGNITUDE, not power-average.

    Scattered pilots move by D_x each symbol, so averaging the spectra first
    smears the comb to period D_x (the union).  Averaging the |DFT| of each
    symbol's power spectrum keeps the comb line at D_x*D_y no matter where the
    offset sits, and averages the noise down.  This is what distinguishes D_y.
    """
    periods = periods or DX_ALL
    acc, n = None, None
    for s in syms:
        q = _detrend(np.abs(s[lo:hi]) ** 2)
        n = len(q)
        Q = np.abs(np.fft.rfft(q * np.hanning(n)))
        acc = Q if acc is None else acc + Q
    acc /= len(syms)
    base = float(np.median(acc[1:]))
    out = []
    for d in periods:
        idx = int(round(n / d))
        if not 1 <= idx < len(acc) - 1:
            continue
        out.append({"dx": d, "legal": d in DX_LEGAL,
                    "line_ratio": float(acc[max(1, idx - 2):idx + 3].max()
                                        / max(base, 1e-30))})
    out.sort(key=lambda z: -z["line_ratio"])
    return out


DY_LEGAL = (1, 2, 4)
DY_CONTROL = (3, 5)


def pilot_hypothesis(syms, lo, hi, dxs=None, dys=None):
    """Direct (D_x, D_y) hypothesis test -- no DFT, so no leakage ambiguity.

    A/322 scattered pilots: carrier k of symbol l is a pilot iff
        k mod (D_x*D_y) == D_x * (l mod D_y)
    Pilots are power-boosted, so for the TRUE (D_x, D_y) -- and the right
    symbol-index origin -- the mean power on those carriers exceeds the mean
    power on the rest.  Every hypothesis is scored the same way, and ILLEGAL
    D_x / D_y values are scored alongside as the null.
    """
    dxs = dxs or tuple(sorted(set(DX_LEGAL) | set(DX_CONTROL)))
    dys = dys or tuple(sorted(set(DY_LEGAL) | set(DY_CONTROL)))
    S = np.abs(np.array(syms))[:, lo:hi] ** 2
    # divide out the channel/filter shape symbol by symbol
    k = max(65, (S.shape[1] // 64) | 1)
    ker = np.ones(k) / k
    for i in range(S.shape[0]):
        base = np.convolve(S[i], ker, mode="same")
        S[i] = S[i] / np.maximum(base, 1e-30)
    S = S[:, k:-k]
    n = S.shape[1]
    kk = np.arange(n)
    out = []
    for dx in dxs:
        for dy in dys:
            m = dx * dy
            if m > n // 8:
                continue
            kmod = kk % m
            best = None
            for off in range(dy):
                num = den = 0.0
                nn = nd = 0
                for l in range(S.shape[0]):
                    g = (l + off) % dy
                    mask = kmod == (dx * g)
                    num += S[l][mask].sum(); nn += int(mask.sum())
                    den += S[l][~mask].sum(); nd += n - int(mask.sum())
                r = (num / max(nn, 1)) / max(den / max(nd, 1), 1e-30)
                if best is None or r > best:
                    best = r
            out.append({"dx": dx, "dy": dy, "boost": float(best),
                        "legal": dx in DX_LEGAL and dy in DY_LEGAL})
    out.sort(key=lambda z: -z["boost"])
    return out


def occupied(avg, fs, nfft, drop_db=6.0):
    """Active-subcarrier span from the averaged spectrum (fftshifted).

    Must be robust to a single huge out-of-band spike: RF33's neighbour RF34 is
    ATSC 1.0 and its unmodulated 8-VSB pilot lands 3.31 MHz above RF33 centre,
    inside the resampled passband.  Normalizing by max() therefore fails, and
    did (it reported "4 active subcarriers").  Use a smoothed dB spectrum, a
    percentile reference, and the LONGEST CONTIGUOUS run.
    """
    k = max(9, (nfft // 512) | 1)
    sm = np.convolve(avg, np.ones(k) / k, mode="same")
    p = 10.0 * np.log10(sm + 1e-30)
    ref = float(np.percentile(p, 75))
    on = p > ref - drop_db
    # longest contiguous run of True
    d = np.diff(np.r_[0, on.astype(np.int8), 0])
    starts, ends = np.flatnonzero(d == 1), np.flatnonzero(d == -1)
    if not len(starts):
        return 0, 0, 0.0
    j = int(np.argmax(ends - starts))
    lo, hi = int(starts[j]), int(ends[j])
    return lo, hi, (hi - lo) * fs / nfft


def analyse(y, seg_start, seg_len, fs_post, label, rep, ffts, gis):
    seg = y[seg_start: seg_start + seg_len]
    print(f"\n--- {label} ({len(seg)} samples = {len(seg)/fs_post*1e3:.1f} ms) ---")
    rows = grid_search(seg, ffts, gis)
    if not rows:
        print("  segment too short for any candidate")
        return None
    print(f"{'FFT':>7} {'GI':>6} {'pitch':>7} {'legal':>6} {'fold pk/med':>12} "
          f"{'nsym':>5}")
    for r in rows[:10]:
        print(f"{r['fft']:7d} {r['gi']:6d} {r['period']:7d} "
              f"{'yes' if r['legal'] else 'CTRL':>6} {r['score']:12.2f} "
              f"{r['n_folded']:5d}")
    best = rows[0]
    rep.setdefault("segments", {})[label] = {"top": rows[:10]}
    if not best["legal"]:
        print("  *** winner is a CONTROL (illegal) hypothesis -- measurement "
              "is NOT trustworthy here")
        return None
    nfft, gi = best["fft"], best["gi"]
    margin = best["score"] / max(rows[1]["score"], 1e-30)
    print(f"  MEASURED  FFT={nfft} ({GI_MENU[gi]})  pitch={nfft+gi} samples")
    print(f"            symbol duration {(nfft+gi)/fs_post*1e3:.4f} ms   "
          f"subcarrier spacing {fs_post/nfft:.3f} Hz")
    print(f"            margin over runner-up ({rows[1]['fft']}/{rows[1]['gi']}): "
          f"{margin:.2f}x")

    # symbol boundary: fold phase is the START of the CP window
    ph = best["phase"]
    syms, pos = [], seg_start + ph + gi
    while pos + nfft <= len(y) and len(syms) < 96:
        syms.append(np.fft.fftshift(np.fft.fft(y[pos:pos + nfft])))
        pos += nfft + gi
    meas = {"fft": nfft, "gi": gi, "gi_label": GI_MENU[gi], "pitch": nfft + gi,
            "tsym_ms": (nfft + gi) / fs_post * 1e3, "df_hz": fs_post / nfft,
            "margin": margin, "n_symbols_fft": len(syms)}
    if len(syms) >= 4:
        S = np.abs(np.array(syms)) ** 2
        avg = S.mean(axis=0)
        lo, hi, bwhz = occupied(avg, fs_post, nfft)
        print(f"            occupied {bwhz/1e6:.4f} MHz over {hi-lo} "
              f"subcarriers (of {nfft})")
        band = avg[lo:hi]
        lines = spectral_line(band)
        print("  pilot line test on mean |X(k)|^2 (union over symbols): " +
              ", ".join(f"{d['dx']}{'' if d['legal'] else '*'}:"
                        f"{d['line_ratio']:.1f}" for d in lines[:8]))
        one = per_symbol_line(syms, lo, hi)
        print("  pilot line test per-symbol (expect D_x*D_y; * = ILLEGAL "
              "control period): " +
              ", ".join(f"{d['dx']}{'' if d['legal'] else '*'}:"
                        f"{d['line_ratio']:.1f}" for d in one[:8]))
        hyp = pilot_hypothesis(syms[:24], lo, hi)
        print("  (D_x,D_y) boost test [* = ILLEGAL control]: " +
              ", ".join(f"({h['dx']},{h['dy']}){'' if h['legal'] else '*'}"
                        f":{h['boost']:.3f}" for h in hyp[:8]))
        best_legal = next((h for h in hyp if h["legal"]), None)
        best_ctrl = next((h for h in hyp if not h["legal"]), None)
        if best_legal and best_ctrl:
            print(f"    best legal ({best_legal['dx']},{best_legal['dy']}) "
                  f"{best_legal['boost']:.3f}  vs  best control "
                  f"({best_ctrl['dx']},{best_ctrl['dy']}) {best_ctrl['boost']:.3f}"
                  f"  -> {'SEPARATED' if best_legal['boost'] > best_ctrl['boost']*1.15 else 'NOT SEPARATED (inconclusive)'}")
        meas["pilot_hypotheses"] = hyp[:8]
        meas.update({"occupied_hz": bwhz, "active_carriers": int(hi - lo),
                     "union_lines": lines[:5], "single_lines": one[:5]})
    rep["segments"][label]["measured"] = meas
    return meas


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("path")
    ap.add_argument("--rate", type=float, default=8e6)
    ap.add_argument("--start-sec", type=float, default=0.0)
    ap.add_argument("--span-sec", type=float, default=0.30)
    ap.add_argument("--fs-post", type=float, default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    path = a.path if os.path.isabs(a.path) else os.path.join(HERE, a.path)
    fmt = guess_format(path)
    rep = {"path": path, "rate": a.rate}

    print(f"=== M2 STAGE B === {os.path.basename(path)}")
    head = read_block(path, int(a.start_sec * a.rate), int(0.30 * a.rate), fmt)
    y6 = resample_to(head, a.rate, bs.FS)
    hits = find_bootstraps(y6)
    if not hits:
        print("no bootstrap found -- cannot anchor the frame")
        return 1
    h = hits[0]
    claimed = h["fields"]
    fs_post = a.fs_post or claimed["post_bootstrap_sample_rate_hz"]
    cfo = h["fine_cfo_hz"]
    print(f"anchor bootstrap at {h['position']/bs.FS*1e3:.3f} ms into block, "
          f"mratio {h['mean_peak_ratio']:.1f}, CFO {cfo:.1f} Hz")
    print(f"CLAIMED by bootstrap bits: bw={claimed['system_bandwidth_str']}, "
          f"fs_post={fs_post/1e6:g} Msps, preamble_structure="
          f"{claimed['preamble_structure']}")

    f0 = int(a.start_sec * a.rate) + int(round(h["position"] * a.rate / bs.FS))
    x = read_block(path, f0, int(a.span_sec * a.rate), fmt)
    y = resample_to(x, a.rate, fs_post).astype(np.complex128)
    y *= np.exp(-2j * np.pi * cfo * np.arange(len(y)) / fs_post)
    boot_len = int(round(4 * bs.N_SYM * fs_post / bs.FS))
    print(f"block {len(y)} samples @ {fs_post/1e6:g} Msps; bootstrap = first "
          f"{boot_len} ({boot_len/fs_post*1e3:.3f} ms)")
    rep["bootstrap_len_post"] = boot_len
    rep["cfo_hz"] = cfo

    ffts = sorted(set(FFT_MENU) | set(FFT_CONTROLS))
    gis = sorted(set(GI_MENU) | set(GI_CONTROLS))

    # ---- subframe boundary: head-to-head between the two winning grids ----
    # ATSC 3.0 frames are bootstrap + preamble + 1..n SUBFRAMES, and each
    # subframe may carry its OWN FFT size / GI.  Measure the boundary rather
    # than assuming one grid holds for the whole frame.
    def head_to_head(hyps, win_ms=25.0, hop_ms=2.5):
        profs = {h: cp_profile(y, h[0], h[1]) for h in hyps}
        w = int(win_ms * fs_post / 1000.0)
        hop = int(hop_ms * fs_post / 1000.0)
        rows = []
        t = boot_len
        while t + w <= min(len(y), int(0.250 * fs_post)):
            sc = {}
            for h, p in profs.items():
                sc[h] = fold_score(p[t:t + w], h[0] + h[1])[0]
            win_h = max(sc, key=sc.get)
            rows.append({"t_ms": t / fs_post * 1e3,
                         "scores": {f"{a}/{b}": round(v, 2)
                                    for (a, b), v in sc.items()},
                         "winner": f"{win_h[0]}/{win_h[1]}"})
            t += hop
        return rows

    print("\n--- subframe boundary: head-to-head grid scores (25 ms window) ---")
    hyps = [(8192, 1536), (16384, 1536)]
    h2h = head_to_head(hyps)
    for r in h2h:
        bar = "  ".join(f"{k}={v:5.2f}" for k, v in r["scores"].items())
        print(f"  t={r['t_ms']:7.2f} ms   {bar}   -> {r['winner']}")
    rep["head_to_head"] = h2h

    # ---- symbol census: exactly how many symbols does each grid hold? ----
    def census(nfft, gi, t0, t1):
        p = cp_profile(y, nfft, gi)
        T = nfft + gi
        _, ph, _ = fold_score(p[t0:t1], T)
        ph += t0
        while ph - T >= boot_len:
            ph -= T
        vals, k = [], 0
        while ph + k * T + 128 < len(p):
            i = ph + k * T
            vals.append(float(p[max(0, i - 96): i + 96].max()))
            k += 1
        return ph, np.array(vals), float(np.median(p))

    print("\n--- symbol census (CP peak per expected symbol slot) ---")
    rep["census"] = {}
    for nfft, gi, t0, t1 in ((8192, 1536, boot_len, boot_len + int(0.030 * fs_post)),
                             (16384, 1536, int(0.060 * fs_post), int(0.230 * fs_post))):
        ph, v, floor = census(nfft, gi, t0, t1)
        T = nfft + gi
        # a slot counts as a real symbol only if its CP peak stands half-way
        # between the profile's own background and the peak level.  Using a
        # fraction of the median of v alone is far too permissive: it kept
        # counting "symbols" 50 ms past the end of the subframe.
        top = float(np.percentile(v, 90))
        strong = v > (floor + 0.5 * (top - floor))
        runs = np.flatnonzero(np.diff(np.r_[0, strong.astype(np.int8), 0]))
        s, e = runs[::2], runs[1::2]
        j = int(np.argmax(e - s)) if len(s) else 0
        n_sym = int(e[j] - s[j]) if len(s) else 0
        t_start = (ph + s[j] * T) / fs_post * 1e3 if len(s) else float("nan")
        t_end = (ph + e[j] * T) / fs_post * 1e3 if len(s) else float("nan")
        print(f"  FFT {nfft:6d} GI {gi}: {n_sym} contiguous symbols, "
              f"t = {t_start:.3f} .. {t_end:.3f} ms  "
              f"(pitch {T} samples = {T/fs_post*1e3:.4f} ms)")
        rep["census"][f"{nfft}/{gi}"] = {"n_symbols": n_sym, "t_start_ms": t_start,
                                         "t_end_ms": t_end, "phase": int(ph)}

    # ---- frame map: which grid wins in each slice of the frame? ----
    # ATSC 3.0 frames are bootstrap + preamble + 1..n SUBFRAMES, and subframes
    # may each carry their own FFT size / GI.  Measure, do not assume one grid.
    print("\n--- frame map (winning grid per 10 ms slice) ---")
    fmap = []
    win = int(0.010 * fs_post)
    t = boot_len
    while t + win <= min(len(y), int(0.250 * fs_post)):
        rows = grid_search(y[t:t + win], ffts, gis, min_symbols=3)
        if rows:
            r = rows[0]
            fmap.append({"t_ms": t / fs_post * 1e3, "fft": r["fft"],
                         "gi": r["gi"], "score": r["score"],
                         "legal": r["legal"]})
            print(f"  t={t/fs_post*1e3:7.2f} ms  ->  FFT {r['fft']:6d} "
                  f"GI {r['gi']:5d} pitch {r['period']:6d}  "
                  f"score {r['score']:5.2f}  {'' if r['legal'] else '(CTRL)'}")
        t += win
    rep["frame_map"] = fmap

    # payload: a long stretch well inside the frame, ~100 symbols folded
    pay = analyse(y, int(0.060 * fs_post), int(0.180 * fs_post), fs_post,
                  "PAYLOAD 60-240 ms", rep, ffts, gis)
    # preamble: the few symbols immediately after the bootstrap
    pre = analyse(y, boot_len, int(0.020 * fs_post), fs_post,
                  "PREAMBLE+ first 20 ms after bootstrap", rep, ffts, gis)

    print("\n=== CLAIMED (bootstrap bits, Stage A) vs MEASURED (signal, Stage B) ===")
    print(f"  post-bootstrap fs   CLAIMED {fs_post/1e6:g} Msps")
    if pay:
        nc = pay.get("active_carriers", 0)
        print(f"  payload FFT         MEASURED {pay['fft']}  "
              f"GI {pay['gi']} ({pay['gi_label']})  "
              f"Tsym {pay['tsym_ms']:.4f} ms  df {pay['df_hz']:.2f} Hz")
        print(f"  occupied bandwidth  MEASURED {pay.get('occupied_hz',0)/1e6:.4f}"
              f" MHz / {nc} carriers   CLAIMED channel "
              f"{claimed['system_bandwidth_str']}")
    if pre:
        print(f"  preamble region     MEASURED FFT {pre['fft']} GI {pre['gi']} "
              f"({pre['gi_label']})   preamble_structure={claimed['preamble_structure']}")
    rep["claimed"] = claimed
    out = a.out or os.path.join(HERE, "m2_stage_b_" +
                                os.path.splitext(os.path.basename(path))[0] + ".json")
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(rep, fh, indent=1, default=str)
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
