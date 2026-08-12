#!/usr/bin/env python3
"""M2 Stage B (pilots) -- find the A/322 scattered-pilot pattern from air.

No radio.  Offline only.  Companion to m2_stage_b.py.

The idea that makes this work without knowing a single pilot VALUE:

    A scattered pilot on carrier k carries the same reference symbol every
    time it recurs, and it recurs every D_y OFDM symbols.  Over a static
    channel that means X_l(k) == X_{l+D_y}(k), exactly, while a data carrier
    changes every symbol.  So

        M_d(k) = |SUM_l X_l(k) conj(X_{l+d}(k))| / SUM_l |X_l(k)|^2

    is large at pilot carriers when d is a multiple of D_y, and small
    everywhere else.  It is also immune to a constant per-symbol phase
    rotation (residual CFO / common phase error), because that rotation is
    identical for every term in the sum.

    -> the smallest d that lights up IS D_y.
    -> the carrier spacing of the carriers that light up IS D_x.
    -> carriers that light up at EVERY d (including d=1, odd) are the
       CONTINUAL pilots, which are present in every symbol.

Two things this leans on and must therefore be stated:
  ASSUMPTION P1  the channel is static over the analysis window (a few tens of
                 ms).  Violating it lowers M uniformly; it cannot manufacture
                 a lattice.
  ASSUMPTION P2  carrier indices are counted relative to the CENTRE of the
                 measured active band.  A/322 counts them from the centre
                 subcarrier; if our centre is off by c the whole lattice shifts
                 by c, which changes the reported residue but not D_x / D_y.

ILLEGAL control values of D_x and D_y are scored alongside the legal ones, so
"we found a pattern" is never just "we looked for one".

Usage:
    python m2_pilots.py hit_rf33.cs16 --rate 8e6
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
from scipy.signal import firwin, lfilter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

from atsc3 import bootstrap as bs                              # noqa: E402
from atsc3.capture import guess_format, read_block             # noqa: E402
from m2_stage_a import find_bootstraps                         # noqa: E402
from m2_stage_b import resample_to, cp_profile, fold_score, occupied  # noqa: E402

DX_LEGAL = (3, 4, 6, 8, 12, 16, 24, 32)
DX_CONTROL = (5, 7, 9, 10, 11, 13, 14, 20)
DY_LEGAL = (1, 2, 4)
DY_CONTROL = (3, 5, 6)

# ASSUMPTION P3: A/322 active-carrier counts for a 6 MHz channel, "normal"
# carrier mode.  Used only as a CROSS-CHECK against the measured span.
NACT_NOMINAL = {8192: 6913, 16384: 13825, 32768: 27649}


def load(path, rate, fmt=None, span_sec=0.30, start_sec=0.0, notch=True):
    """Read, anchor on a bootstrap, resample to the post-bootstrap rate,
    de-rotate the CFO, and (critically) band-limit."""
    fmt = fmt or guess_format(path)
    head = read_block(path, int(start_sec * rate), int(0.30 * rate), fmt)
    hits = find_bootstraps(resample_to(head, rate, bs.FS))
    if not hits:
        raise RuntimeError("no bootstrap found")
    h = hits[0]
    fs_post = h["fields"]["post_bootstrap_sample_rate_hz"]
    cfo = h["fine_cfo_hz"]
    f0 = int(start_sec * rate) + int(round(h["position"] * rate / bs.FS))
    x = read_block(path, f0, int(span_sec * rate), fmt)
    y = resample_to(x, rate, fs_post).astype(np.complex128)
    y *= np.exp(-2j * np.pi * cfo * np.arange(len(y)) / fs_post)
    if notch:
        # RF33's upper neighbour RF34 is ATSC 1.0 and its unmodulated 8-VSB
        # pilot lands 3.31 MHz above RF33 centre -- INSIDE the 6.912 Msps
        # passband, and 7 dB ABOVE the wanted signal.  Being a pure carrier it
        # correlates at EVERY lag, which simultaneously crushes the CP peak
        # (0.65 instead of 0.999) and lifts the background 7x.  Band-limiting
        # to the ATSC 3.0 occupied width fixes both.  Measured effect on the
        # CP fold ratio: 3.1 -> 40.0.
        nt = 257
        y = lfilter(firwin(nt, 2.95e6, fs=fs_post), 1.0, y)[nt // 2:]
    return y, fs_post, cfo, h


def symbols(y, nfft, gi, t0, t1, n_max=80, skip=0):
    """FFT every OFDM symbol in [t0,t1).  Returns (array, active lo, hi)."""
    p = cp_profile(y, nfft, gi)
    T = nfft + gi
    _, ph, _ = fold_score(p[t0:t1], T)
    ph += t0
    while ph - T >= t0:
        ph -= T
    out, pos = [], ph + gi + skip * T
    while pos + nfft <= len(y) and len(out) < n_max and pos < t1:
        out.append(np.fft.fftshift(np.fft.fft(y[pos:pos + nfft])))
        pos += T
    A = np.array(out)
    avg = (np.abs(A) ** 2).mean(axis=0)
    lo, hi, _ = occupied(avg, 1.0, nfft)
    return A, lo, hi, ph


def m_stat(A, d):
    num = np.abs((A[:-d] * np.conj(A[d:])).sum(axis=0))
    den = (np.abs(A[:-d]) ** 2).sum(axis=0)
    return num / np.maximum(den, 1e-30)


def analyse_pilots(A, lo, hi, nfft, tag, rep):
    ctr = (lo + hi) // 2
    print(f"\n=== {tag}: {len(A)} symbols, active {hi-lo} carriers "
          f"(nominal {NACT_NOMINAL.get(nfft, '?')}), centre bin {ctr} ===")
    k = np.arange(lo, hi) - ctr

    # 1) which lag d lights up?  -> D_y
    Ms = {}
    print("  lag scan (M = repeat coherence):")
    for d in (1, 2, 3, 4, 5, 6, 8, 12):
        if d >= len(A) - 2:
            continue
        M = m_stat(A, d)[lo:hi]
        Ms[d] = M
        print(f"    d={d:2d}  p50 {np.percentile(M,50):.3f}  p90 "
              f"{np.percentile(M,90):.3f}  p99 {np.percentile(M,99):.3f}  "
              f"#(>0.4) {int((M>0.4).sum()):6d}")

    # 2) continual pilots: high at EVERY lag, including odd ones
    odd = [Ms[d] for d in (1, 3, 5) if d in Ms]
    cp_idx = np.array([], dtype=int)
    if odd:
        cont = np.minimum.reduce(odd)
        cp_idx = np.flatnonzero(cont > 0.5)
        print(f"  CONTINUAL pilots (min over odd lags > 0.5): {len(cp_idx)} "
              f"carriers")
        if 0 < len(cp_idx) < 2000:
            rel = k[cp_idx]
            for mod in (2, 3, 4, 6, 8, 12, 16, 24, 32):
                if len(np.unique(rel % mod)) == 1:
                    print(f"    all lie on rel-carrier == {rel[0] % mod} "
                          f"(mod {mod})")
            print(f"    span {rel.min()} .. {rel.max()}, "
                  f"median spacing {int(np.median(np.diff(np.sort(rel))))}")

    # 3) (D_x, D_y) hypothesis test on M_{D_y}, legal AND control
    rows = []
    for dy in sorted(set(DY_LEGAL) | set(DY_CONTROL)):
        if dy not in Ms:
            continue
        M = Ms[dy].copy()
        if len(cp_idx):                    # continual pilots would flatter
            M[cp_idx] = np.nan             # every hypothesis -- exclude them
        good = ~np.isnan(M)
        for dx in sorted(set(DX_LEGAL) | set(DX_CONTROL)):
            best = None
            for r in range(dx):
                sel = good & ((k % dx) == r)
                oth = good & ((k % dx) != r)
                if sel.sum() < 20:
                    continue
                v = float(np.nanmean(M[sel]) / max(np.nanmean(M[oth]), 1e-30))
                if best is None or v > best[0]:
                    best = (v, r)
            if best:
                rows.append({"dx": dx, "dy": dy, "ratio": best[0],
                             "residue": best[1],
                             "legal": dx in DX_LEGAL and dy in DY_LEGAL})
    rows.sort(key=lambda z: -z["ratio"])
    print("  (D_x,D_y) test on M_{D_y}  [* = ILLEGAL control]:")
    for z in rows[:8]:
        print(f"    ({z['dx']:2d},{z['dy']}){'' if z['legal'] else '*'} "
              f"ratio {z['ratio']:6.3f}  residue {z['residue']}")
    bl = next((z for z in rows if z["legal"]), None)
    bc = next((z for z in rows if not z["legal"]), None)
    verdict = "inconclusive"
    if bl and bc:
        sep = bl["ratio"] / max(bc["ratio"], 1e-9)
        verdict = (f"SP{bl['dx']}_{bl['dy']}" if sep > 1.10 else "inconclusive")
        print(f"  best legal SP{bl['dx']}_{bl['dy']} {bl['ratio']:.3f} vs best "
              f"control ({bc['dx']},{bc['dy']}) {bc['ratio']:.3f} -> "
              f"{'PATTERN = SP%d_%d' % (bl['dx'], bl['dy']) if sep > 1.10 else 'NOT SEPARATED'}")
    rep[tag] = {"n_symbols": int(len(A)), "active": int(hi - lo),
                "n_continual": int(len(cp_idx)), "top": rows[:8],
                "verdict": verdict}
    return verdict, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("path")
    ap.add_argument("--rate", type=float, default=8e6)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    path = a.path if os.path.isabs(a.path) else os.path.join(HERE, a.path)
    y, fsp, cfo, h = load(path, a.rate)
    boot = int(round(4 * bs.N_SYM * fsp / bs.FS))
    print(f"=== M2 pilots === {os.path.basename(path)}  fs_post {fsp/1e6:g} "
          f"Msps  CFO {cfo:.1f} Hz  bootstrap {boot} samples")
    rep = {}

    # subframe 0: 8K / GI1536, 2.000 .. 52.667 ms (measured in Stage B)
    A0, lo0, hi0, ph0 = symbols(y, 8192, 1536, boot, boot + int(0.050 * fsp),
                                n_max=36)
    # split it: the PREAMBLE is the first symbols and A/322 gives it its own
    # pilot pattern, so an abrupt change in the pattern LOCATES the preamble.
    analyse_pilots(A0[:6], lo0, hi0, 8192, "SF0 8K symbols 0-5 (preamble?)", rep)
    analyse_pilots(A0[6:], lo0, hi0, 8192, "SF0 8K symbols 6+ (data)", rep)
    analyse_pilots(A0, lo0, hi0, 8192, "SF0 8K all", rep)

    # subframe 1: 16K / GI1536, 52.667 .. 247.111 ms
    A1, lo1, hi1, ph1 = symbols(y, 16384, 1536, int(0.055 * fsp),
                                int(0.235 * fsp), n_max=64)
    analyse_pilots(A1, lo1, hi1, 16384, "SF1 16K", rep)

    out = a.out or os.path.join(HERE, "m2_pilots_" +
                                os.path.splitext(os.path.basename(path))[0] + ".json")
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(rep, fh, indent=1, default=str)
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
