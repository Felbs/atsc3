#!/usr/bin/env python3
"""e124_ldm_snr.py -- E69's cancellation measurement on the live chain's own FEC batches.

Input: a bank written by `ATSC3_FEC_CELL_DUMP=DIR python -m atsc3 watch ...`
(every Nth FEC batch: the de-interleaved Core cells and the decoded Baseband
Packets).  For every BCH-clean Block the Core is re-ENCODED from its payload
(m6_bicm.PlpChain.encode -- the true round trip, E122 PASS), fitted to the
received cells by least squares and subtracted:

    residual = enhanced + noise,       P_enh = P_core / 10^(inj/10)
    SNR_core = P_core / P_residual     SNR_enh = P_enh / (P_residual - P_enh)

plus E69's negative control (a wrong payload must NOT cancel).  This is the
number that says whether the Enhanced layer is reachable from a capture
before anyone builds its decoder into the live chain: E122 measured the
thresholds (RF25's 256QAM 7/15 needs ~12.5 dB on its own layer, i.e. the
injection level + 12.5 of Core SNR; RF30's 64QAM 6/15 at 4 dB needs 13 / 9 dB).

THE INJECTION LEVEL IS AN INDEX.  L1D_plp_ldm_injection_level is 5 bits in
0.5 dB steps (A/331; m3_spec.ldm_injection_db); RF25 signals index 10 =
5.0 dB.  The plan used to print the raw index ("inj10") and 13.20-13.24 ran
this tool with --inj 10, which put the Enhanced layer 18 dB short where it is
8.6 (13.25).  Pass --inj-idx with the signalled index, or --inj with dB.

Run: python lab/e124_ldm_snr.py --dump DIR [--inj 10] [--limit N]
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


def db(x):
    return 10 * np.log10(max(float(x), 1e-30))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", required=True)
    ap.add_argument("--inj", type=float, default=None, help="injection level in dB (default: --inj-idx 10 = 5.0 dB)")
    ap.add_argument("--inj-idx", type=int, default=10,
                    help="the signalled L1D_plp_ldm_injection_level INDEX (RF25: 10 = 5.0 dB; RF30: 8 = 4.0 dB)")
    ap.add_argument("--need", type=float, default=12.5,
                    help="Enhanced layer's own-SNR threshold, dB (E122: RF25 12.5, RF30 9)")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    if a.inj is None:
        import m3_spec as S
        a.inj = S.ldm_injection_db(a.inj_idx)
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
    rows, ctl = [], []
    nblk = nclean = 0
    for f in files:
        z = np.load(f, allow_pickle=True)
        cells = z["cells"]; pkts = z["pkts"]; ok = np.asarray(z["bch_ok"])
        conv = np.asarray(z["conv"], bool) if "conv" in z.files else np.ones(len(cells), bool)
        if ok.ndim == 0:
            # The writer banks bch_ok as ONE flag per batch (measured on the
            # 6-core box's bank, 13.21: every file, parent and workers).  A
            # Block is clean when it converged and carries a packet -- in this
            # chain bch_zero == fec_converged -- so judge per Block on that.
            ok = conv & np.array([len(pk) > 0 for pk in pkts], bool)
        nblk += len(cells)
        for i in range(len(cells)):
            if not ok[i] or len(pkts[i]) == 0:
                continue
            nclean += 1
            pay = np.unpackbits(np.frombuffer(pkts[i], np.uint8))[:kp]
            key = pkts[i]
            s = cache.get(key)
            if s is None:
                s = cache[key] = np.asarray(ch.encode(pay), np.complex128)
            r = np.asarray(cells[i], np.complex128)
            gain = np.vdot(s, r) / np.vdot(s, s)
            res = r - gain * s
            p_core = float(abs(gain) ** 2 * np.mean(np.abs(s) ** 2))
            rows.append((float(np.mean(np.abs(r) ** 2)), p_core, float(np.mean(np.abs(res) ** 2)), abs(gain)))
            if len(ctl) < 40:                      # the control, on a sample
                bad = pay.copy(); bad[::7] ^= 1
                sb = np.asarray(ch.encode(bad), np.complex128)
                gb = np.vdot(sb, r) / np.vdot(sb, sb)
                ctl.append(float(np.mean(np.abs(r - gb * sb) ** 2)) / rows[-1][0])
    print(f"{len(files)} batches, {nblk} Blocks, {nclean} BCH-clean and re-encoded; "
          f"Core {z0['mod']} {z0['rate']} {int(z0['ninner'])}, injection {a.inj:.1f} dB")
    if not rows:
        print("  no clean Block -- nothing to measure")
        return 1
    A = np.array(rows)
    p_tot, p_core, p_res = A[:, 0].mean(), A[:, 1].mean(), A[:, 2].mean()
    p_enh = p_core / 10 ** (a.inj / 10)
    p_noise = p_res - p_enh
    print(f"  fitted |gain| median {np.median(A[:, 3]):.4f}  (1.0 = the Core as transmitted, cells normalised)")
    print(f"  P_total {p_tot:.4f}  P_core {p_core:.4f} ({db(p_core/p_tot):+.2f} dB of total)  "
          f"P_residual {p_res:.4f} ({db(p_res/p_tot):+.2f} dB)  cancellation gain {db(p_tot/p_res):.2f} dB")
    print(f"  negative control: wrong payload leaves {np.mean(ctl):.3f} of the received power "
          f"({'fails to cancel, as required' if np.mean(ctl) > 0.8 else '!! CONTROL DID NOT FAIL'})")
    print(f"  Core-layer SNR  (P_core / P_residual)        {db(p_core/p_res):6.2f} dB")
    if p_noise <= 0:
        print(f"  P_noise <= 0: the residual ({p_res:.4f}) is at or below the predicted Enhanced power "
              f"({p_enh:.4f}) -- noise below the measurement floor; Enhanced SNR unbounded above")
        snr_enh = float("inf")
    else:
        snr_enh = db(p_enh / p_noise)
        print(f"  P_enh (from injection) {p_enh:.4f}   P_noise {p_noise:.4f}   "
              f"ENHANCED-LAYER SNR {snr_enh:6.2f} dB")
    print(f"  Enhanced layer needs {a.need:.1f} dB (E122) -> margin "
          f"{'unbounded' if snr_enh == float('inf') else f'{snr_enh - a.need:+.2f} dB'}  "
          f"({'should decode' if snr_enh >= a.need else 'SHORT'})")
    # per-batch spread, for fades
    per = []
    for f in files:
        pass
    return 0


if __name__ == "__main__":
    sys.exit(main())
