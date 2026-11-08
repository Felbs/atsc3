#!/usr/bin/env python3
"""e124_residual.py -- what the residual after Core cancellation is MADE of.

13.24: across four RF25 captures the E124 Core-layer SNR rose only 2.05 -> 3.40 dB
while convergence went 27% -> 100% and the antenna gained 13 dB.  The residual
after a perfect re-encode cancellation scales with the signal, so it is not
thermal noise -- and 13.25 found it does not move with the estimator's levers
either.  Before anyone touches the receiver for the Enhanced layer, say WHAT
the residual is.  On the same cell bank `e124_ldm_snr.py` reads (every Nth FEC
batch of a replay under ATSC3_FEC_CELL_DUMP), for every BCH-clean Block with
re-encoded Core s, received cells r, LS gain g and residual e = r - g s:

  1. PHASE or AMPLITUDE?  Project the residual on the Core symbol,
     q = e * conj(s) / |s|^2, and compare Var(Im q) against Var(Re q):
     1.0 = isotropic (noise, the Enhanced layer, which is independent of the
     Core); >> 1 = a common-phase error (phase noise, CPE residual); << 1 =
     a gain error.
  2. SLOW or FAST?  Re-fit the gain per SEGMENT of the Block (1, 4, 16, 64,
     256, 1024 segments of the de-interleaved order).  What a piecewise gain
     removes is drift the one-gain fit could not follow.  The de-interleaved
     order is NOT time order (the CTI spreads a Block's cells across ~1024
     rows), so a segment mixes time -- a residual that falls steeply with
     the segment count is nonetheless structured, not white.
  3. GAUSSIAN?  Kurtosis of Re e and Im e (3.0 = Gaussian).  The Enhanced
     layer is 256QAM-NUC, near-Gaussian; heavy tails are something else
     (impulsive interference, a few wrong cells, an equaliser fallback).
  4. HOW MUCH is unexplained?  P_res / P_core per Block (p5 / p50 / p95) and
     the EXCESS over the Enhanced layer's own share 10^(-inj/10) -- the part
     the Enhanced layer does not account for.
  5. Does the excess track the Block's own power (signal-proportional) or
     stay flat (additive)?  Correlation of P_res with P_core across Blocks,
     and the slope of a least-squares line P_res = a + b * P_core.

Run: python lab/e124_residual.py --dump DIR [--inj 10] [--limit N]
"""
from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)
import m6_bicm as B                                                 # noqa: E402

SEGS = (1, 4, 16, 64, 256, 1024)


def db(x):
    return 10 * np.log10(max(float(x), 1e-30))


def kurt(x):
    x = x - x.mean()
    v = float(np.mean(x * x))
    return float(np.mean(x ** 4) / max(v * v, 1e-30))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", required=True)
    ap.add_argument("--inj", type=float, default=10.0)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    files = sorted(glob.glob(os.path.join(a.dump, "*.npz")))
    if a.limit:
        files = files[:a.limit]
    if not files:
        print("no bank under", a.dump)
        return 2
    z0 = np.load(files[0], allow_pickle=True)
    ch = B.PlpChain(int(z0["ninner"]), str(z0["mod"]), str(z0["rate"]))
    kp = int(z0["kpayload"])
    cache = {}
    n_tot = n_clean = 0
    p_core, p_res, ratio = [], [], []
    var_re, var_im = [], []
    k_re, k_im = [], []
    seg_res = {m: [] for m in SEGS}
    per_batch = []
    for f in files:
        z = np.load(f, allow_pickle=True)
        cells = z["cells"]; pkts = z["pkts"]; ok = np.asarray(z["bch_ok"])
        conv = np.asarray(z["conv"], bool) if "conv" in z.files else np.ones(len(cells), bool)
        if ok.ndim == 0:
            ok = conv & np.array([len(pk) > 0 for pk in pkts], bool)
        n_tot += len(cells)
        b_core = b_res = 0.0; nb = 0
        for i in range(len(cells)):
            if not ok[i] or len(pkts[i]) == 0:
                continue
            n_clean += 1
            pay = np.unpackbits(np.frombuffer(pkts[i], np.uint8))[:kp]
            s = cache.get(pkts[i])
            if s is None:
                s = cache[pkts[i]] = np.asarray(ch.encode(pay), np.complex128)
            r = np.asarray(cells[i], np.complex128)
            g = np.vdot(s, r) / np.vdot(s, s)
            e = r - g * s
            pc = float(abs(g) ** 2 * np.mean(np.abs(s) ** 2))
            pr = float(np.mean(np.abs(e) ** 2))
            p_core.append(pc); p_res.append(pr); ratio.append(pr / pc)
            b_core += pc; b_res += pr; nb += 1
            # 1. projection on the Core symbol
            q = e * np.conj(s) / np.maximum(np.abs(s) ** 2, 1e-30)
            var_re.append(float(np.var(q.real))); var_im.append(float(np.var(q.imag)))
            # 3. tails
            k_re.append(kurt(e.real)); k_im.append(kurt(e.imag))
            # 2. piecewise gain
            n = len(s)
            for m in SEGS:
                L = n // m
                ss = s[:L * m].reshape(m, L); rr = r[:L * m].reshape(m, L)
                gm = np.sum(np.conj(ss) * rr, axis=1) / np.sum(np.abs(ss) ** 2, axis=1)
                em = rr - gm[:, None] * ss
                seg_res[m].append(float(np.mean(np.abs(em) ** 2)) / pc)
        if nb:
            per_batch.append((os.path.basename(f), nb, b_res / b_core))
    if not p_core:
        print("no clean Block -- nothing to measure")
        return 1
    pc_ = np.array(p_core); pr_ = np.array(p_res); rt = np.array(ratio)
    enh = 10 ** (-a.inj / 10)
    print(f"{len(files)} batches, {n_tot} Blocks, {n_clean} BCH-clean and re-encoded; "
          f"Core {z0['mod']} {z0['rate']} {int(z0['ninner'])}, injection {a.inj:.1f} dB "
          f"(the Enhanced layer's own share of the residual = {enh:.3f} of P_core)")
    print(f"  P_res/P_core  mean {rt.mean():.4f} ({db(1/rt.mean()):.2f} dB Core SNR)  "
          f"p5/p50/p95 {np.percentile(rt,5):.4f} / {np.percentile(rt,50):.4f} / {np.percentile(rt,95):.4f}")
    print(f"  EXCESS over the Enhanced layer: {rt.mean()-enh:.4f} of P_core = {db(rt.mean()-enh):+.2f} dB "
          f"relative to the Core (an estimation error of that relative power would produce it)")
    vr, vi = np.mean(var_re), np.mean(var_im)
    print(f"\n  1. residual projected on the Core symbol: Var(Im)/Var(Re) = {vi/vr:.3f}  "
          f"(1.0 isotropic; >>1 phase error; <<1 gain error)")
    print(f"     per-Block p5/p50/p95 of the ratio: "
          f"{np.percentile(np.array(var_im)/np.array(var_re),[5,50,95]).round(3).tolist()}")
    print(f"\n  2. piecewise gain -- P_res/P_core after re-fitting the gain per segment of the Block:")
    base = np.mean(seg_res[1])
    for m in SEGS:
        v = np.mean(seg_res[m])
        print(f"     {m:5d} segment(s) ({len(np.asarray(cells[0]))//m:6d} cells each): {v:.4f}  "
              f"({100*(1-v/base):5.1f}% of the one-gain residual removed)")
    print(f"     (a white residual loses only ~m/N of its power to m fitted gains: "
          f"{100*1024/len(np.asarray(cells[0])):.1f}% at 1024)")
    print(f"\n  3. kurtosis of the residual: Re {np.mean(k_re):.2f}  Im {np.mean(k_im):.2f}  "
          f"(Gaussian 3.0; per-Block p95 Re {np.percentile(k_re,95):.2f})")
    cc = np.corrcoef(pc_, pr_)[0, 1]
    bfit, afit = np.polyfit(pc_, pr_, 1)
    print(f"\n  5. P_res against P_core across {len(pc_)} Blocks: corr {cc:+.3f}; "
          f"LS line P_res = {afit:.4f} + {bfit:.4f} x P_core  "
          f"(P_core p5/p95 {np.percentile(pc_,5):.3f} / {np.percentile(pc_,95):.3f}; "
          f"a signal-proportional impairment has intercept ~0 and slope = the ratio above)")
    print(f"\n  per batch (name, clean Blocks, P_res/P_core):")
    for nm, nb, r_ in per_batch:
        print(f"     {nm:28s} {nb:4d}  {r_:.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
