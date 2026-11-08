#!/usr/bin/env python3
"""e127_ce_probe.py -- the channel as the LDM demod's pilots see it (13.25 Part B).

Input: the bank `ATSC3_LDM_CE_DUMP` writes from `FrameDemod._ce_smooth` (every
Nth Frame: the raw scattered-pilot estimates on the (grid x symbol) matrix Hg,
the per-symbol complex gains cl divided out before smoothing, the smoothed grid
Hsm, the per-symbol pilot noise floor sig2_l from second differences, the
fallback mask).  E124 measured the data cells' residual after a perfect Core
cancellation at -3 dB relative to the Core on every capture, isotropic, white
and Gaussian (e124_residual).  Two things the pilots can say that the data
cells cannot:

  1. PILOT SNR.  Per symbol, |Hsm|^2 (the channel power on the grid) against
     sig2_l (the pilot estimate's own noise, channel curvature cancelled) and
     against the raw-minus-smoothed residual.  If the pilots see 15-20 dB
     where the data cells see 3, the impairment is between the pilots and
     the data cells -- the interpolation -- not in the air.
  2. DELAY SPREAD against the grid's reach.  The scattered pilots sample the
     channel every dx carriers, so the estimate can only represent echoes
     inside +-Tu/(2 dx) (one bin of the grid's IFFT = Tu/ngrid).  The IFFT
     of one symbol's grid estimate is the channel impulse response folded
     into that span: the power fraction in the OUTER half of the span, and
     the RMS delay spread, say whether the channel exceeds what the linear
     interpolation between grid points can follow.  An echo beyond the reach
     aliases into every data carrier as a signal-proportional, isotropic,
     white-after-the-CTI error -- E124's residual to the letter.

Run: python lab/e127_ce_probe.py --dump DIR [--limit N]
"""
from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np


def db(x):
    return 10 * np.log10(max(float(x), 1e-30))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", required=True)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--fs", type=float, default=6.912e6, help="sample rate the nfft is quoted at (Tu = nfft/fs)")
    a = ap.parse_args()
    files = sorted(glob.glob(os.path.join(a.dump, "*.npz")))
    if a.limit:
        files = files[:a.limit]
    if not files:
        print("no bank under", a.dump)
        return 2
    rows = []
    for f in files:
        z = np.load(f)
        Hg, Mg, Hsm, cl = z["Hg"].astype(np.complex128), z["Mg"], z["Hsm"].astype(np.complex128), z["cl"]
        res_l, sig2_l, ok = z["res_l"], z["sig2_l"], z["ok"]
        dx, nfft, nsym, noc = int(z["dx"]), int(z["nfft"]), int(z["nsym"]), int(z["noc"])
        ngrid = Hg.shape[0]
        Tu_us = nfft / a.fs * 1e6
        reach_us = Tu_us / dx                       # the grid's unambiguous delay span
        # 1. pilot SNR, three ways, on the symbols the smoother accepted
        Hn = Hg / cl[None, :]                       # per-symbol gain divided out, as the smoother sees it
        fin = np.isfinite(res_l) & ok
        p_ch = np.array([np.mean(np.abs(Hsm[:, l]) ** 2) for l in range(nsym)])
        snr_curv = p_ch[fin] / np.maximum(sig2_l[fin], 1e-30)       # noise from second differences
        snr_resid = p_ch[fin] / np.maximum(res_l[fin], 1e-30)       # raw - smoothed (noise + what the boxcar hides)
        # 2. the impulse response from each symbol's grid (only symbols with every grid point present)
        outer, rms, peak_db = [], [], []
        for l in range(nsym):
            if not (Mg[:, l].all() and ok[l]):
                continue
            h = np.fft.ifft(Hn[:, l])
            p = np.abs(h) ** 2
            p /= p.sum()
            # fold so the strongest tap sits at 0; the span is one grid period
            k0 = int(np.argmax(p))
            p = np.roll(p, -k0)
            tau = np.arange(ngrid); tau[tau > ngrid // 2] -= ngrid        # bins, signed
            t_us = tau * (reach_us / ngrid)
            outer.append(float(p[np.abs(tau) > ngrid // 4].sum()))          # beyond +-reach/4
            rms.append(float(np.sqrt(np.sum(p * t_us ** 2) - np.sum(p * t_us) ** 2)))
            peak_db.append(db(p[0]))
        # 3. the INTERPOLATION the equaliser does: linear between grid points (jg/wtg
        #    in _sm_tables) against the band-limited (DFT) interpolation of the SAME
        #    smoothed grid -- the error, relative to the channel, at the data carriers.
        lin_err = []
        kk = np.arange(noc)
        gc = np.arange(ngrid) * dx
        jg = np.clip(np.searchsorted(gc, kk, side="right") - 1, 0, ngrid - 2)
        wt = (kk - gc[jg]) / dx
        for l in range(nsym):
            if not ok[l]:
                continue
            Hs = Hsm[:, l]
            h = np.fft.ifft(Hs)
            # zero-pad in the delay domain -> the grid's band-limited interpolant at every carrier
            k0 = int(np.argmax(np.abs(h)))
            h = np.roll(h, -k0)
            hz = np.zeros(ngrid * dx, complex)
            half = ngrid // 2
            hz[:ngrid - half] = h[:ngrid - half]
            hz[-half:] = h[-half:]
            Hd = np.fft.fft(hz)                 # Hd[dx*k] == Hs[k] exactly; the rest is the band-limited interpolant
            Hd *= np.exp(-2j * np.pi * k0 * np.arange(ngrid * dx) / (ngrid * dx))   # undo the roll
            Hd = Hd[:noc]
            Hl = Hs[jg] * (1 - wt) + Hs[jg + 1] * wt
            data = (kk % dx) != 0
            lin_err.append(float(np.mean(np.abs(Hl[data] - Hd[data]) ** 2) / np.mean(np.abs(Hd[data]) ** 2)))
        rows.append(dict(f=os.path.basename(f), nsym=nsym, ok=int(ok.sum()), ngrid=ngrid, dx=dx, nfft=nfft,
                         lin=db(np.median(lin_err)) if lin_err else float("nan"),
                         snr_curv=db(np.median(snr_curv)) if len(snr_curv) else float("nan"),
                         snr_resid=db(np.median(snr_resid)) if len(snr_resid) else float("nan"),
                         outer=np.median(outer) if outer else float("nan"),
                         rms=np.median(rms) if rms else float("nan"),
                         peak=np.median(peak_db) if peak_db else float("nan"),
                         reach=reach_us, Tu=Tu_us, nsyms_cir=len(outer)))
    r0 = rows[0]
    print(f"{len(rows)} Frames; nfft {r0['nfft']} (Tu {r0['Tu']:.1f} us), scattered-pilot dx {r0['dx']} -> grid of {r0['ngrid']} "
          f"points, unambiguous delay span Tu/dx = {r0['reach']:.2f} us (bin {r0['reach']/r0['ngrid']*1e3:.1f} ns)")
    print(f"  {'Frame file':22s} {'sym ok':>7s} {'pilot SNR (curv)':>16s} {'pilot SNR (resid)':>17s} "
          f"{'CIR peak':>8s} {'RMS delay us':>12s} {'power beyond +-span/4':>21s} {'lin-interp err':>14s}")
    for r in rows:
        print(f"  {r['f']:22s} {r['ok']:3d}/{r['nsym']:<3d} {r['snr_curv']:16.2f} {r['snr_resid']:17.2f} "
              f"{r['peak']:8.2f} {r['rms']:12.3f} {100*r['outer']:20.2f}% {r['lin']:11.2f} dB")
    li = np.array([r["lin"] for r in rows])
    print(f"  linear interpolation between grid points (dx {r0['dx']}) against the DFT interpolant of the same grid: "
          f"error {np.nanmedian(li):.2f} dB relative to the channel at the data carriers (median over Frames)")
    sc = np.array([r["snr_curv"] for r in rows]); sr = np.array([r["snr_resid"] for r in rows])
    ou = np.array([r["outer"] for r in rows]); rm = np.array([r["rms"] for r in rows])
    print(f"\n  median over Frames: pilot SNR {np.nanmedian(sc):.2f} dB (second differences) / {np.nanmedian(sr):.2f} dB "
          f"(raw - smoothed); CIR RMS delay {np.nanmedian(rm):.3f} us; power beyond +-span/4 {100*np.nanmedian(ou):.2f}%")
    print(f"  E124's data-cell Core SNR on the same captures is 2-3.4 dB.  Pilot SNR far above it = the loss is between "
          f"the pilots and the data cells (interpolation / reach); pilot SNR near it = the air or the front end.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
