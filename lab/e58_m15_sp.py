#!/usr/bin/env python3
"""E58 -- the composite threshold (M15 methodology) with the exact-BP rescue.

M15's refine run measured the operating point (QPSK 9/15 64K core +
256QAM-NUC @ 5.0 dB injection, gaussian demapper, own decoder) at 6.80 dB
channel SNR with the normalized-min-sum alpha ladder.  This measures the SAME
threshold with the E58 two-stage decoder: min-sum ladder first, exact
sum-product (float64 log-tanh) only on ladder failures.  Any threshold drop
is decoder-approximation gap recovered -- no demap change, no channel change,
same 3/3 convention, same seeds.

Also reports the ladder-only threshold at the same step so the delta is
engine-paired rather than read against M15's recorded 6.80.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m6_bicm as B                                              # noqa: E402
import m3_ldpc as LD                                             # noqa: E402
from m15_ldm_llr import powers, llr_gaussian                     # noqa: E402
from e58_weighted import sum_product_decode_batch                # noqa: E402


def trial(inj_db, snr_db, seed, mode, iters=100, sp_iters=100, sp_cut=1500):
    rng = np.random.default_rng(seed)
    ch = B.PlpChain(64800, "QPSK", "9/15", iters=iters)
    pay = rng.integers(0, 2, ch.kpayload).astype(np.uint8)
    x_core = ch.encode(pay)
    enh_pts = B.points_for("256QAM", "7/15")
    x_enh = enh_pts[rng.integers(0, len(enh_pts), len(x_core))]
    p_core, p_enh = powers(inj_db)
    n0 = 10 ** (-snr_db / 10.0)
    n = (rng.standard_normal(len(x_core))
         + 1j * rng.standard_normal(len(x_core))) / np.sqrt(2)
    y = (np.sqrt(p_core) * x_core + np.sqrt(p_enh) * x_enh
         + n * np.sqrt(n0))
    q = llr_gaussian(y, ch, p_core, p_enh, n0)
    lam = np.empty(ch.ninner)
    lam[ch.lam_of_q] = q
    # min-sum ladder (float32 batch = the shipping fast engine)
    best_bad = 1 << 30
    for al in ch.ALPHA_LADDER:
        b, c, i, d = LD.min_sum_decode_batch(
            lam[None, :].astype(np.float32), ch.checks, ch.ninner,
            iters=iters, alpha=al, dtype=np.float32)
        if c[0]:
            bb = ch.baseband_packet(b[0])
            return (bool(bb["bch_ok"])
                    and bool(np.array_equal(bb["bits"], pay))), "minsum"
        best_bad = min(best_bad, int(d[0]))
    if mode == "ladder":
        return False, "minsum"
    if best_bad >= sp_cut:
        return False, "minsum"
    bs, cs, its, bad = sum_product_decode_batch(lam[None, :], ch.checks,
                                                ch.ninner, iters=sp_iters)
    if cs[0]:
        bb = ch.baseband_packet(bs[0])
        return (bool(bb["bch_ok"])
                and bool(np.array_equal(bb["bits"], pay))), "sp"
    return False, "sp-fail"


def threshold(mode, lo, hi, step, trials, inj_db):
    snr = lo
    while snr <= hi + 1e-9:
        ok = True
        hows = []
        for k in range(trials):
            good, how = trial(inj_db, snr, 1000 + k, mode)
            hows.append(how)
            if not good:
                ok = False
                break
        print(f"    {mode:8s} SNR {snr:5.2f} dB  "
              f"{'PASS 3/3' if ok else 'fail'}  [{','.join(hows)}]",
              flush=True)
        if ok:
            return snr
        snr += step
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inj", type=float, default=5.0)
    ap.add_argument("--lo", type=float, default=6.0)
    ap.add_argument("--hi", type=float, default=7.6)
    ap.add_argument("--step", type=float, default=0.1)
    ap.add_argument("--trials", type=int, default=3)
    ap.add_argument("--json", default="e58_m15_sp.json")
    a = ap.parse_args()
    try:
        import psutil
        psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    except Exception:                                           # noqa: BLE001
        pass
    out = {}
    for mode in ("ladder", "twostage"):
        t0 = time.time()
        out[mode] = threshold(mode, a.lo, a.hi, a.step, a.trials, a.inj)
        print(f"  -> {mode}: "
              f"{out[mode] if out[mode] is not None else 'NOT REACHED'} dB"
              f"  ({time.time()-t0:.0f}s)")
    if out.get("ladder") is not None and out.get("twostage") is not None:
        out["gain_db"] = out["ladder"] - out["twostage"]
        print(f"  exact-BP rescue moves the composite threshold "
              f"{out['gain_db']:+.2f} dB")
    json.dump(out, open(os.path.join(HERE, a.json), "w"), indent=1)


if __name__ == "__main__":
    main()
