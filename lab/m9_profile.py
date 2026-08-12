#!/usr/bin/env python3
"""M9 Step 1 -- WHERE THE TIME GOES.  Instrument the real decode, don't guess.

M7 recorded "the LDPC decoder is the whole cost and it is pure NumPy".  That
is a claim, not a measurement.  This module wraps the actual functions the
m7_route decode path calls -- by monkeypatching the module attributes the
callers look up, so nothing in the chain is edited or reimplemented -- and
reports seconds and percent per stage against ONE reference: the wall-clock
duration of the RF the frames represent.

    1 Frame = 1708032 samples at 6.912 Msps = 0.247111 s of air

so "x-real-time" for a run of N frames is N * 0.247111 / wall.

Nesting is handled by reporting inclusive time for the outer stage and
"of which" lines for the inner ones, so the percentages of the LEAF rows sum
to the total and nothing is counted twice.

Usage:
    python m9_profile.py long_rf33.cs16 --rate 8e6 --frames 8
"""
from __future__ import annotations

import argparse
import collections
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m3_ldpc as LD                                              # noqa: E402
import m4_scrambler as SC                                         # noqa: E402
import m6_bicm as B                                               # noqa: E402
import m6_cells as C                                              # noqa: E402
import m6_payload as P6                                           # noqa: E402
import m6_tbi as T                                                # noqa: E402

FRAME_SEC = C.FRAME_SAMPLES / 6.912e6

T_ACC = collections.Counter()
N_ACC = collections.Counter()


class stage:
    """Context manager + decorator, accumulating into T_ACC."""

    def __init__(self, name):
        self.name = name

    def __enter__(self):
        self.t = time.perf_counter()
        return self

    def __exit__(self, *a):
        T_ACC[self.name] += time.perf_counter() - self.t
        N_ACC[self.name] += 1
        return False


def wrap(mod, attr, name):
    fn = getattr(mod, attr)
    if getattr(fn, "_m9", None):
        return
    def inner(*a, **k):
        t = time.perf_counter()
        try:
            return fn(*a, **k)
        finally:
            T_ACC[name] += time.perf_counter() - t
            N_ACC[name] += 1
    inner._m9 = True
    inner._orig = fn
    setattr(mod, attr, inner)


def install():
    # front end
    wrap(C, "demod_data", "  ofdm.demod_data (FFT+chest+eq)")
    wrap(C, "demod_preamble", "  ofdm.demod_preamble")
    wrap(C, "cpe_correct", "cpe_correct")
    wrap(C, "dummy_check", "dummy_check")
    wrap(C, "constellation_regions", "constellation_regions")
    import m3_freqint as FI
    wrap(FI, "deinterleave", "  freq_deinterleave")
    # BICM
    wrap(T, "fec_block", "time_deinterleave (HTI)")
    wrap(B, "demap_llr", "  demap_llr (max-log)")
    wrap(LD, "min_sum_decode", "LDPC min-sum (inclusive)")
    wrap(LD, "pack", "  ldpc.pack (rebuilt every block!)")
    wrap(LD, "parity_check", "  ldpc.parity_check (H build)")
    wrap(B, "bch_syndrome", "BCH syndrome")
    wrap(SC, "descramble", "descramble")


# ---------------------------------------------------------------------------


def profile_frames(path, rate, fmt, start, n_frames, chain_cache=True):
    """Run the exact m7_route decode path with timers, single process."""
    span = 0.05 + n_frames * FRAME_SEC
    with stage("capture load (disk -> complex)"):
        y = C.load_frame(path, rate, fmt=fmt, start_sec=start,
                         span_sec=span * (rate / 6.912e6))
    t_prev = None
    nfec = conv = bch = 0
    packets = []
    for fi in range(n_frames):
        if (fi + 1) * C.FRAME_SAMPLES + 400000 > len(y):
            break
        centre = C.BOOTSTRAP if t_prev is None else t_prev + C.FRAME_SAMPLES
        with stage("fine_timing (41x demod_data)"):
            t0, coh = C.fine_timing(y, span=20, centre=centre)
        t_prev = t0
        with stage("cell_pool (inclusive)"):
            pool, info = C.cell_pool(y, t0, {})
        reg, dv = C.constellation_regions(len(pool))
        alph = [B.points_for("64QAM", "11/15"), B.QPSK_POINTS]
        pool, _ = C.cpe_correct(pool, info["symbol_of"], alph, reg, dv)
        C.dummy_check(pool)
        with stage("PlpChain build (H + pack)"):
            ch16 = B.PlpChain(16200, "QPSK", "2/15")
            ch0 = B.PlpChain(16200, "64QAM", "11/15")
        r16 = ch16.decode(pool[C.PLP16["start"]:
                               C.PLP16["start"] + C.PLP16["size"]])
        if r16["converged"]:
            ch16.baseband_packet(r16["bits"])
        plp0 = pool[:C.PLP0["size"]]
        nrows, ncols = ch0.cells_per_fec, C.PLP0["n_fec_max"] // C.PLP0["nti"]
        per_ti = len(plp0) // C.PLP0["nti"]
        for ti in range(C.PLP0["nti"]):
            seg = plp0[ti * per_ti:(ti + 1) * per_ti]
            for j in range(ncols):
                blk = T.fec_block(seg, j, nrows, ncols, 0)
                with stage("llr map (demap + bit de-int)"):
                    lam = ch0.llr(blk)
                bits, ok, it, bad = LD.min_sum_decode(
                    lam, ch0.checks, ch0.ninner, iters=ch0.iters)
                nfec += 1
                if ok:
                    conv += 1
                    bb = ch0.baseband_packet(bits)
                    bch += bb["bch_ok"]
                    packets.append(bb["bytes"])
    return dict(frames=fi + 1 if n_frames else 0, nfec=nfec, conv=conv,
                bch=bch, packets=len(packets))


def report(wall, nframes, title="PROFILE"):
    air = nframes * FRAME_SEC
    print(f"\n  === {title} ===")
    print(f"  {nframes} Frames = {air:.3f} s of air; wall {wall:.3f} s "
          f"-> {air/wall:.3f}x real time\n")
    print(f"    {'stage':44s} {'calls':>8s} {'sec':>9s} {'%wall':>7s} "
          f"{'ms/frame':>9s}")
    print("    " + "-" * 80)
    for k in sorted(T_ACC, key=lambda k: -T_ACC[k]):
        print(f"    {k:44s} {N_ACC[k]:8d} {T_ACC[k]:9.3f} "
              f"{100*T_ACC[k]/wall:6.1f}% {1000*T_ACC[k]/max(nframes,1):9.2f}")
    return air / wall


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("capture")
    ap.add_argument("--rate", type=float, required=True)
    ap.add_argument("--fmt")
    ap.add_argument("--start", type=float, default=0.0)
    ap.add_argument("--frames", type=int, default=8)
    ap.add_argument("--accel", default="none", choices=("none", "cpu"))
    a = ap.parse_args()
    path = a.capture if os.path.isabs(a.capture) else os.path.join(HERE,
                                                                   a.capture)
    print("M9 Step 1 -- stage profile of the existing CPU chain")
    print("=" * 72)
    if a.accel == "cpu":
        import m9_accel
        m9_accel.install()
        print("  m9_accel INSTALLED (memoisation + packed BCH)")
    install()
    t0 = time.perf_counter()
    r = profile_frames(path, a.rate, a.fmt, a.start, a.frames)
    wall = time.perf_counter() - t0
    x = report(wall, r["frames"])
    print(f"\n    FEC Blocks {r['conv']}/{r['nfec']} converged, "
          f"{r['bch']} BCH zero")
    print(f"\n    SINGLE PROCESS: {x:.3f}x real time")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
