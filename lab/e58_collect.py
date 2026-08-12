#!/usr/bin/env python3
"""E58 -- collect with TIME-SMOOTHED channel estimation (the 2x1D upgrade).

THE MECHANISM
-------------
`m6_cells._demod_data_g` estimates the channel of each OFDM symbol from that
symbol's own pilots alone: SP12_4 on a non-SBS symbol is one pilot every
dx*dy = 48 carriers, each estimate carrying the full channel noise.  At the
1-5 dB SNR this capture lives at, the CHANNEL-ESTIMATE noise is a first-order
part of the demap noise budget.  The channel itself is quasi-static over tens
of symbols (fixed antenna; the breathing is 0.25-1 s ~ 200-800 symbols, one
symbol = 1.26 ms), so pilot estimates may be averaged across symbols:

  * a pilot carrier c = 12j repeats every dy = 4 symbols (and the SBS symbols
    0 / 188 carry ALL multiples of dx = 12) -> a +-W-symbol boxcar average
    holds ~2W/dy + 1 looks at the same channel sample;
  * the smoothed grid is DENSER in frequency too: 12-carrier spacing at every
    symbol instead of 48 (the 2x1D estimator, albacore lever #2's shape).

CHANGE DETECTOR (real fades must still track): per symbol, the residual of
the symbol's OWN raw pilot estimates against the smoothed grid is compared to
the frame median; a symbol whose residual ratio exceeds --detect (default 6x)
is demodulated by the untouched per-symbol path instead.  W = 0 runs the
untouched path for every symbol -- the identity gate.

Everything after equalization (frequency de-interleave, SBS trim, pool
order, CPE, dummy SNR) is the E49 path, unchanged.  Output artifacts have
the same shape as e49_stream collect, so e58_weighted decode runs on them.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m3_spec as S                  # noqa: E402
import m6_cells as C6                # noqa: E402
import m10_core as M10               # noqa: E402
import spec_pilots as P              # noqa: E402
from m3_preamble import fft_bins     # noqa: E402
from m2_pilots import load           # noqa: E402
from e49_dummy_snr import dummy_snr  # noqa: E402


def smoothed_pool(y, t0, g, W=8, detect=6.0, report=None):
    """cell_pool_g with time-smoothed data-symbol channel estimation."""
    rep = report if report is not None else {}
    fi0 = g.np_sym

    parts, owner, cohs = [], [], []
    pre = []
    for s in range(g.np_sym):
        z, coh = C6._demod_preamble_g(y, t0, g, s)
        pre.append(C6.FI.deinterleave(z, g.pre_nfft, s, direction="forward",
                                      toggle="i"))
        cohs.append(coh)
    pre = np.concatenate(pre) if pre else np.zeros(0, complex)
    spare = pre[g.l1b_cells + g.l1d_cells:]
    parts.append(spare)
    owner.append(np.full(len(spare), -1, int))

    noc = P.NOC[(g.nfft, g.cred)]
    lo, _ = S.carrier_abs_range(g.nfft, g.cred)
    ref = np.array(S.pilot_values(noc, 1.0))
    dx, dy = P.dxdy(g.pattern)
    ngrid = (noc - 1) // dx + 1                       # noc-1 is a multiple of dx
    assert (noc - 1) % dx == 0

    # ---- pass 1: FFT every data symbol, raw pilot estimates --------------
    Y = np.empty((g.nsym, g.nfft), np.complex128)
    pk_of, hp_of = [], []
    for l in range(g.nsym):
        w = (t0 + (g.pre_nfft + g.pre_gi) * g.np_sym + (g.nfft + g.gi) * l
             + g.gi)
        Y[l] = np.fft.fftshift(np.fft.fft(y[w:w + g.nfft]))
        sbs = g.is_sbs(l)
        pk = (np.arange(0, noc, dx) if sbs
              else np.arange(dx * (l % dy), noc, dx * dy))
        pk = np.unique(np.concatenate([pk, [0, noc - 1]]))
        hp = Y[l][fft_bins(g.nfft, lo + pk)] / ref[pk]
        pk_of.append(pk)
        hp_of.append(hp)

    # ---- pass 2: the (grid-carrier x symbol) matrix, boxcar in time ------
    Hg = np.zeros((ngrid, g.nsym), np.complex128)
    Mg = np.zeros((ngrid, g.nsym), bool)
    for l in range(g.nsym):
        j = pk_of[l] // dx
        on = pk_of[l] % dx == 0                       # edge noc-1 IS on grid
        Hg[j[on], l] = hp_of[l][on]
        Mg[j[on], l] = True
    # NaN-free boxcar via cumulative sums over the mask
    csH = np.concatenate([np.zeros((ngrid, 1), complex),
                          np.cumsum(np.where(Mg, Hg, 0), axis=1)], axis=1)
    csM = np.concatenate([np.zeros((ngrid, 1)),
                          np.cumsum(Mg, axis=1)], axis=1)
    losym = np.maximum(np.arange(g.nsym) - W, 0)
    hisym = np.minimum(np.arange(g.nsym) + W + 1, g.nsym)
    num = csH[:, hisym] - csH[:, losym]
    cnt = csM[:, hisym] - csM[:, losym]
    ok = cnt > 0
    Hsm = np.where(ok, num / np.maximum(cnt, 1), 0.0)

    # ---- pass 3: equalize each symbol off the smoothed dense grid --------
    dcoh, fell_back, ratios = [], 0, []
    grid_c = np.arange(ngrid) * dx
    res_l = np.empty(g.nsym)
    for l in range(g.nsym):
        sm = Hsm[:, l][ok[:, l]]
        gc = grid_c[ok[:, l]]
        Hf = C6._interp(np.arange(noc), gc, sm)
        res_l[l] = float(np.mean(np.abs(hp_of[l] - Hf[pk_of[l]]) ** 2))
    med = max(float(np.median(res_l)), 1e-30)

    for l in range(g.nsym):
        sbs = g.is_sbs(l)
        ratio = res_l[l] / med
        ratios.append(ratio)
        if W == 0 or ratio > detect:
            z, coh = C6._demod_data_g(y, t0, g, l)    # untouched path
            fell_back += W != 0
        else:
            sm = Hsm[:, l][ok[:, l]]
            gc = grid_c[ok[:, l]]
            Hf = C6._interp(np.arange(noc), gc, sm)
            used = np.zeros(noc, bool)
            used[np.array(sorted(P.pilot_carriers(g.nfft, g.cred, g.pattern,
                                                  l, sbs)), int)] = True
            d = np.flatnonzero(~used)
            z = Y[l][fft_bins(g.nfft, lo + d)] / Hf[d]
            hp = hp_of[l]
            coh = float(abs(np.sum(hp[:-1] * np.conj(hp[1:])))
                        / max(np.sum(np.abs(hp) ** 2), 1e-30))
        dcoh.append(coh)
        x = (C6.FI.deinterleave(z, g.nfft, l + fi0, direction="forward",
                                toggle="i") if g.freq_interleaver
             else np.asarray(z))
        if sbs:
            x = x[g.null_low:len(x) - g.null_high]
        parts.append(x)
        owner.append(np.full(len(x), l, int))

    pool = np.concatenate(parts)
    rep.update(pool=int(len(pool)), preamble_coherence=cohs,
               data_coherence_mean=float(np.mean(dcoh)),
               fell_back=int(fell_back),
               detect_ratio_max=float(np.max(ratios)))
    return pool, dict(symbol_of=np.concatenate(owner), report=rep)


# ---------------------------------------------------------------------------
# whole-capture collect (e49_stream shape)
# ---------------------------------------------------------------------------

_G = {}


def _init(capture, rate, l1, off, W, detect):
    js = json.load(open(l1))
    _G["g"] = M10.geometry_from_json(js, os.path.basename(capture))
    _G.update(capture=capture, rate=rate, off=off, W=W, detect=detect)
    _G["plps"] = {p["id"]: p for p in _G["g"].plps}
    _G["core_total"] = sum(p["size"] for p in _G["g"].plps if p["layer"] == 0)


def _one(n):
    g = _G["g"]
    per = M10.frame_samples(g) / M10.FS_POST
    s = max(0.0, _G["off"] + n * per - 0.01)
    try:
        span = (M10.frame_samples(g) + 40000) / M10.FS_POST
        y, _fs, _cfo, _h = load(_G["capture"], _G["rate"], fmt=None,
                                span_sec=span, start_sec=s)
        t0, coh = M10.fine_t0(y, g)
        rep = {}
        pool, info = smoothed_pool(y, t0, g, W=_G["W"], detect=_G["detect"],
                                   report=rep)
        pool = M10.cpe(pool, info["symbol_of"], g, _G["plps"][0],
                       dummy_start=_G["core_total"])
        snr = dummy_snr(pool, _G["core_total"])
        p0 = _G["plps"][0]
        c0 = pool[p0["start"]:p0["start"] + p0["size"]].astype(np.complex64)
        diag = dict(frame=n, ok=True, t=s, snr_db=snr["snr_db"],
                    data_coh=rep["data_coherence_mean"],
                    fell_back=rep["fell_back"], pool=rep["pool"])
        return n, c0.tobytes(), diag
    except Exception as exc:                                    # noqa: BLE001
        return n, None, dict(frame=n, ok=False, t=s,
                             err=f"{type(exc).__name__}: {exc}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("capture")
    ap.add_argument("--rate", type=float, required=True)
    ap.add_argument("--l1", required=True)
    ap.add_argument("--start-frame", type=int, default=2)
    ap.add_argument("--frames", type=int, required=True)
    ap.add_argument("--window", type=int, default=8)
    ap.add_argument("--detect", type=float, default=6.0)
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--prefix", required=True)
    a = ap.parse_args()
    try:
        import psutil
        psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    except Exception:                                           # noqa: BLE001
        pass

    js = json.load(open(a.l1))
    g = M10.geometry_from_json(js, os.path.basename(a.capture))
    per = M10.frame_samples(g) / M10.FS_POST
    plps = {p["id"]: p for p in g.plps}
    off = a.start_frame * per
    arr0 = np.lib.format.open_memmap(a.prefix + "_plp0.npy", mode="w+",
                                     dtype=np.complex64,
                                     shape=(a.frames, plps[0]["size"]))
    diags = [None] * a.frames
    t0 = time.time()
    done = 0
    with ProcessPoolExecutor(
            max_workers=a.workers, initializer=_init,
            initargs=(os.path.abspath(a.capture), a.rate,
                      os.path.abspath(a.l1), off, a.window, a.detect)) as ex:
        for n, b0, diag in ex.map(_one, range(a.frames), chunksize=4):
            if b0 is not None:
                arr0[n] = np.frombuffer(b0, np.complex64)
            diags[n] = diag
            done += 1
            if done % 25 == 0:
                el = time.time() - t0
                print(f"  {done}/{a.frames} frames ({el:.0f}s, "
                      f"{el/done*1000:.0f} ms/frame)", flush=True)
    arr0.flush()
    ok = [d for d in diags if d.get("ok")]
    v = np.array([d["snr_db"] for d in ok])
    fb = sum(d.get("fell_back", 0) for d in ok)
    print(f"  {len(ok)}/{a.frames} frames  dummy-SNR med {np.median(v):.2f} "
          f"p90 {np.percentile(v,90):.2f} max {v.max():.2f}  "
          f"fallback symbols {fb}")
    json.dump(dict(capture=os.path.basename(a.capture),
                   start_frame=a.start_frame, frames=a.frames,
                   window=a.window, detect=a.detect, diags=diags),
              open(a.prefix + "_diags.json", "w"), indent=1, default=float)
    print(f"  wrote {a.prefix}_plp0.npy / _diags.json")


if __name__ == "__main__":
    main()
