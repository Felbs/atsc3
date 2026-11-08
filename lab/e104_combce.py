#!/usr/bin/env python3
"""E104 -- channel estimation in the delay domain (the pilot comb as an impulse response).

THE IDEA, AND WHERE IT CAME FROM
--------------------------------
Every data symbol's channel response is known only at the pilot carriers, and
the chain fills the carriers between them by LINEAR interpolation (m6_cells
`_interp`, and the same line on the E60 time-smoothed grid).  Linear
interpolation has a cost that does not go away with SNR: every data cell is
equalised by the average of TWO noisy pilots, so it carries their noise on top
of its own.  For the usual pilot boost that is a fixed 0.5 to 1.5 dB.

The pilots are a regular comb (spacing D carriers).  The inverse FFT of that
comb is the channel's impulse response, repeated every N/D samples.  A real
channel has a handful of taps (the direct path and a few echoes) inside the
guard interval; the other delay bins hold only noise.  So: take one period of
the impulse response, zero the bins that are indistinguishable from noise,
transform back -- and the response at EVERY carrier comes from the strong taps
alone, with the pilot noise thrown away in the bins that were zeroed.

The technique is standard OFDM channel estimation (DFT-based / "noise
thresholded CIR"), and we first saw it applied to ATSC 3.0 in OnAir
(github.com/adswill/OnAir, GPL-3.0, core/src/atsc3/atsc3_ofdm.cpp), whose
notes also quote the 0.5-1.5 dB figure.  Nothing here is copied from it; this
is written from the principle above, and it differs in two places on purpose:
the transform length is a multiple of D so the periodic images are exact
copies (no smearing of a non-integer period), and the comb spacing is checked
against the guard interval before the method is trusted at all.

WHEN IT MUST REFUSE
-------------------
The method sees one period, N/D samples.  An echo later than 3/4 of that
period aliases into the next period and is lost.  On RF33's SP4_2 the data
symbols' own comb is spaced 8 (window 1024 < GI 1536): unsafe.  The E60
time-smoothed grid is spaced dx = 4 (window 2048 > GI): safe.  So this estimator
belongs ON the smoothed grid, and `MIN_PERIOD` (the caller sets it to the guard
interval) makes it hand back to linear interpolation anywhere the comb is too
sparse.  It also refuses when more than a quarter of the comb is missing, when
no tap clears the floor, or when "most" bins do (no structure: not a channel).

Everything is gated by the referee the chain already trusts: LDPC convergence
and BCH syndrome zero, block for block, at the SAME injected noise, with linear
interpolation as the control (e104_combce.py main, below).
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

# ---------------------------------------------------------------------------
# the estimator
# ---------------------------------------------------------------------------

MIN_PERIOD = 0          # the caller sets this to the guard interval (samples)
FLOOR_MARGIN = 4.0      # a tap is kept when its power is 6 dB over the floor
DYN_RANGE_DB = 40.0     # ... and never when it is 40 dB under the strongest
LN_4_3 = math.log(4.0 / 3.0)   # lower quartile of an exponential = mean * ln(4/3)

_stats = {"calls": 0, "used": 0, "why": {}}


def _refuse(why):
    _stats["why"][why] = _stats["why"].get(why, 0) + 1
    return None


def h_from_comb(noc, pk, hp, nfft=None):
    """Response at carriers 0..noc-1 from pilots at `pk` (values `hp`), or
    None when linear interpolation should be used instead.  See the module
    docstring for the method and the refusals."""
    _stats["calls"] += 1
    pk = np.asarray(pk, int)
    hp = np.asarray(hp, complex)
    if len(pk) < 16:
        return _refuse("few pilots")
    d = np.diff(pk)
    D = int(np.bincount(d).argmax())              # the comb spacing
    if D < 1:
        return _refuse("no spacing")
    if nfft is None:
        nfft = 1 << int(math.ceil(math.log2(2 * noc)))
    N = D * int(math.ceil(nfft / D))               # period N/D is then exact
    W = N // D
    if MIN_PERIOD and W < MIN_PERIOD:
        return _refuse("comb too sparse for the guard")
    r = int(np.bincount(pk % D).argmax())          # the comb's phase
    on = (pk % D) == r
    comb_k = np.arange(r, noc, D)
    if on.sum() < 3 * len(comb_k) // 4:
        return _refuse("comb mostly missing")
    # fill the comb's holes from their neighbours (the middle carriers, an
    # edge point that is off the comb is simply not used)
    hv = np.interp(comb_k, pk[on], hp[on].real) + 1j * np.interp(comb_k, pk[on], hp[on].imag)
    # one implementation: the batched port in m16_margin (lever 4) does the
    # transform; this reference only regularises the comb and maps the
    # result back onto the carriers when the comb does not start at 0
    import m16_margin as M16
    Hf, used = M16.comb_denoise_grid(hv[:, None], D, nfft, min_period=MIN_PERIOD)
    if not used[0]:
        return _refuse("m16 refused (floor, structure or energy)")
    H = np.empty(noc, complex)
    H[r:] = Hf[0][:noc - r]
    H[:r] = Hf[0][0]
    _stats["used"] += 1
    return H


def selfcheck():
    """Three things that must hold before the estimator touches real air."""
    import spec_pilots as P
    noc, nfft, D = 6913, 8192, 4
    pk = np.arange(0, noc, D)
    pk = np.unique(np.concatenate([pk, [0, noc - 1]]))
    # 1. a flat channel comes back flat, to machine precision
    H = h_from_comb(noc, pk, np.ones(len(pk)), nfft)
    assert H is not None and np.allclose(H, 1.0, atol=1e-10), "flat channel"
    # 2. a two-path channel (echo 700 samples, -8 dB) is recovered at EVERY
    #    carrier from noisy pilots better than linear interpolation is
    rng = np.random.default_rng(7)
    k = np.arange(noc) - (noc - 1) // 2
    a2 = 10 ** (-8 / 20) * np.exp(1j * 0.4)
    Htrue = 1.0 + a2 * np.exp(-2j * np.pi * k * 700 / nfft)
    for snr_db in (10.0, 20.0, 30.0):
        s = 10 ** (-snr_db / 20) / np.sqrt(2)
        hp = Htrue[pk] + s * (rng.standard_normal(len(pk)) + 1j * rng.standard_normal(len(pk)))
        Hc = h_from_comb(noc, pk, hp, nfft)
        Hl = np.interp(np.arange(noc), pk, hp.real) + 1j * np.interp(np.arange(noc), pk, hp.imag)
        ec = np.mean(np.abs(Hc - Htrue) ** 2)
        el = np.mean(np.abs(Hl - Htrue) ** 2)
        print(f"    two-path channel, pilot SNR {snr_db:4.1f} dB: estimate error "
              f"linear {10*np.log10(el):6.1f} dB   comb {10*np.log10(ec):6.1f} dB   "
              f"({10*np.log10(el/ec):+.1f} dB better)")
        assert ec < el, "the comb estimate must beat linear interpolation"
    # 3. pure noise is refused, not 'estimated'
    hn = rng.standard_normal(len(pk)) + 1j * rng.standard_normal(len(pk))
    assert h_from_comb(noc, pk, hn, nfft) is None, "noise must be refused"
    # 4. the guard check refuses a comb that is too sparse
    global MIN_PERIOD
    MIN_PERIOD = 1536
    assert h_from_comb(noc, np.arange(0, noc, 8), np.ones(len(np.arange(0, noc, 8))), nfft) is None
    MIN_PERIOD = 0
    # 5. the batched port in m16_margin (lever 4, what the live chain runs)
    #    equals this per-symbol reference, symbol for symbol
    import m16_margin as M16
    grid_k = np.arange(0, noc, D)
    cols = []
    for sdb in (10.0, 20.0, 30.0):
        s = 10 ** (-sdb / 20) / np.sqrt(2)
        cols.append(Htrue[grid_k] + s * (rng.standard_normal(len(grid_k)) + 1j * rng.standard_normal(len(grid_k))))
    cols.append(rng.standard_normal(len(grid_k)) + 1j * rng.standard_normal(len(grid_k)))   # noise: refused
    Hg = np.array(cols).T
    Hfull, used = M16.comb_denoise_grid(Hg, D, nfft, min_period=1536)
    assert list(used) == [True, True, True, False], f"batched acceptance {list(used)}"
    for l in range(3):
        ref = h_from_comb(noc, grid_k, Hg[:, l], nfft)
        assert ref is not None and np.allclose(Hfull[l], ref, rtol=1e-9, atol=1e-12), f"batched != per-symbol at column {l}"
    print("    selfcheck: all five legs hold")


# ---------------------------------------------------------------------------
# the gate on real air: the same injected noise, two estimators, LDPC/BCH judge
# ---------------------------------------------------------------------------

def _nfft_of(noc):
    import spec_pilots as P
    for (nfft, cred), n in P.NOC.items():
        if n == noc:
            return nfft
    return None


def _interp_comb(kk, pk, hp):
    """Drop-in for m6_cells._interp: the comb estimate, else linear."""
    import m6_cells as C6
    noc = len(kk)
    H = h_from_comb(noc, pk, hp, _nfft_of(noc))
    if H is None:
        return C6._interp_linear(kk, pk, hp)
    return H


def install(on):
    import m6_cells as C6
    if not hasattr(C6, "_interp_linear"):
        C6._interp_linear = C6._interp
    C6._interp = _interp_comb if on else C6._interp_linear


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("capture", nargs="?", default="hit_rf33.cs16")
    ap.add_argument("--rate", type=float, default=8e6)
    ap.add_argument("--fmt", default=None)
    ap.add_argument("--l1", default="m4_l1detail_hit_rf33.json")
    ap.add_argument("--frames", type=int, default=3)
    ap.add_argument("--subframe", type=int, default=1, help="1 = PLP 1 (HTI, 256QAM: the thin margin); 0 = PLP 0 (CTI, partial per frame)")
    ap.add_argument("--plp", type=int, default=None)
    ap.add_argument("--start", type=float, default=0.0)
    ap.add_argument("--snrs", default="none,19,18,17.5,17,16.5,16,15.5,15,14",
                    help="injected noise: signal power over added noise power, dB; 'none' = as recorded")
    ap.add_argument("--W", type=int, default=12, help="E60 time-smoothing window (0 = per-symbol path)")
    ap.add_argument("--iters", type=int, default=50)
    ap.add_argument("--min-period", type=int, default=None,
                    help="smallest delay window (samples) the comb estimate accepts; default = the guard interval")
    ap.add_argument("--threads", type=int, default=12)
    ap.add_argument("--selfcheck", action="store_true")
    ap.add_argument("--json", default=None)
    a = ap.parse_args(argv)

    print("E104 -- delay-domain channel estimate vs linear interpolation")
    print("=" * 72)
    selfcheck()
    if a.selfcheck:
        return 0

    import m6_cells as C6
    import m16_margin as M16
    import m13_sf1 as M13
    from m2_pilots import load, guess_format
    from m10_core import fine_t0, frame_samples
    from concurrent.futures import ThreadPoolExecutor

    path = a.capture if os.path.isabs(a.capture) else os.path.join(HERE, a.capture)
    g0 = C6.Geometry.rf33_legacy()
    if a.subframe == 0:
        # subframe 0 IS the base geometry (its t0 is the Frame's, the Preamble
        # included).  NOTE its PLP 0 rides a CONVOLUTIONAL time interleaver, so
        # one frame alone completes only part of its blocks: the two estimators
        # are still paired, but the counts are not a threshold in themselves.
        g = g0
        plp = dict(id=0, **C6.PLP0)
        bpc = {"QPSK": 2, "16QAM": 4, "64QAM": 6, "256QAM": 8, "1024QAM": 10, "4096QAM": 12}[plp["mod"]]
        plp.setdefault("cells_per_fec", plp["ninner"] // bpc)
        plp.setdefault("inter_subframe", 0)
        off = 0
        fi = None
    else:
        l1p = a.l1 if os.path.isabs(a.l1) else os.path.join(HERE, a.l1)
        js = json.load(open(l1p))
        g = M13.geometry_for_subframe(js["l1_basic"], js["l1_detail"], a.subframe, g0)
        plps = {p["id"]: p for p in g.plps}
        plp = plps[a.plp if a.plp is not None else sorted(plps)[0]]
        off = M13.subframe_offset(g0, a.subframe)
        fi = 0                                         # m13: the interleaver counter resets at the subframe
    fmt = a.fmt or guess_format(path)
    total = frame_samples(g0) + ((g.nfft + g.gi) * g.nsym if a.subframe else 0)
    per = total / M13.FS_POST
    ex = ThreadPoolExecutor(max_workers=a.threads) if a.threads else None

    global MIN_PERIOD
    MIN_PERIOD = int(g.gi) if a.min_period is None else a.min_period
    print(f"  {os.path.basename(path)}: FFT {g.nfft//1024}K GI {g.gi} {g.pattern}, "
          f"subframe {a.subframe} PLP {plp['id']} ({plp['mod']} {plp['rate']}, {plp['size']//plp['cells_per_fec']} blocks/frame), "
          f"smoothing W={a.W}, comb window must be >= {MIN_PERIOD}")

    # the pool comes from the E60 smoothed grid (W=0 is the untouched path)
    def pool_fn(y, t0, gg, shift=0, report=None, fi_offset=None):
        return M16.smoothed_pool_g(y, t0, gg, W=a.W, detect=6.0, report=report)
    C6.cell_pool_g = pool_fn

    snrs = [None if s.strip() == "none" else float(s) for s in a.snrs.split(",")]
    rows = []
    frames = []
    for k in range(a.frames):
        y, _fs, _cfo, _h = load(path, a.rate, fmt=fmt, span_sec=per + 0.06, start_sec=a.start + k * per)
        t0, coh = fine_t0(y, g0)
        frames.append((y, t0))
    psig = float(np.mean([np.mean(np.abs(y) ** 2) for y, _ in frames]))
    print(f"  {len(frames)} frames loaded, signal power {10*np.log10(psig):.2f} dB (as recorded, noise included)")
    print(f"\n  {'inject':>7}  {'estimator':9}  {'blocks ok':>10}  {'BCH ok':>7}  {'median unsat':>12}  {'comb used':>9}  {'s':>5}")
    for snr in snrs:
        for est in ("linear", "comb"):
            install(est == "comb")
            _stats.update(calls=0, used=0, why={})
            conv = bch = n = 0
            unsat = []
            t = time.time()
            for k, (y, t0) in enumerate(frames):
                if snr is None:
                    yn = y
                else:
                    rng = np.random.default_rng(104_000 + k * 97 + int(round(snr * 10)))
                    sn = np.sqrt(psig * 10 ** (-snr / 10) / 2)
                    yn = y + sn * (rng.standard_normal(len(y)) + 1j * rng.standard_normal(len(y)))
                r = M13.decode_subframe(yn, t0 + off, g, plp, fi, iters=a.iters, verbose=False, ex=ex)
                if not r.get("ok"):
                    continue
                conv += r["converged"]; bch += r["bch_ok"]; n += r["n_blocks"]
                if r["converged"] < r["n_blocks"]:
                    unsat.append(r["unsat_median"])
            dt = time.time() - t
            rows.append(dict(inject_db=snr, estimator=est, blocks=n, converged=conv, bch_ok=bch,
                             median_unsat=(int(np.median(unsat)) if unsat else 0),
                             comb_used=_stats["used"], comb_calls=_stats["calls"], refusals=dict(_stats["why"])))
            print(f"  {('as rec' if snr is None else f'{snr:5.1f}'):>7}  {est:9}  {conv:4d}/{n:<5d}  {bch:7d}  "
                  f"{(int(np.median(unsat)) if unsat else 0):12d}  {_stats['used']:4d}/{_stats['calls']:<4d}  {dt:5.1f}")
    install(False)
    if a.json:
        json.dump(dict(capture=os.path.basename(path), W=a.W, frames=a.frames, rows=rows),
                  open(a.json, "w"), indent=1)
    # the verdict: the lowest injection each estimator still decodes >= 95% of its blocks at
    def floor_of(est):
        best = None
        for r in rows:
            if r["estimator"] == est and r["inject_db"] is not None and r["blocks"] and r["converged"] >= 0.95 * r["blocks"]:
                best = r["inject_db"] if best is None else min(best, r["inject_db"])
        return best
    fl, fc = floor_of("linear"), floor_of("comb")
    print(f"\n  95%-of-blocks floor: linear {fl} dB, comb {fc} dB"
          + (f"  ->  the comb estimate buys {fl - fc:+.1f} dB" if fl is not None and fc is not None else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
