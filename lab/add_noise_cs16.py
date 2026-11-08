#!/usr/bin/env python3
"""add_noise_cs16.py -- a banked capture, plus seeded Gaussian noise, as a new cs16.

Why: the cliff gates (E54/E104/E107) add noise to IQ in memory, per Frame.  The
pipeline timing question needs the WHOLE 45 s file in the marginal regime --
the one where most FEC Blocks fail and the LDPC abort does the work -- and
that regime lives on another box's capture.  So make it here: the same air,
a stated number of dB of complex Gaussian noise below the recorded power,
seeded, written once, decodable by `watch --capture` like any other bank.

The noise is relative to the RECORDED power, which already includes the
recording's own noise, so the level is a knob, not an SNR: calibrate it by
the block-convergence fraction it produces, and say so when quoting it.

Run: python lab/add_noise_cs16.py IN.cs16 OUT.cs16 --db 6 [--seed 7]
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("dst")
    ap.add_argument("--db", type=float, required=True,
                    help="noise power this many dB BELOW the recorded power")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--chunk", type=int, default=1 << 22, help="complex samples per chunk")
    a = ap.parse_args()

    nbytes = os.path.getsize(a.src)
    ntot = nbytes // 4
    # pass 1: recorded power (mean |s|^2 in int16 units)
    acc = 0.0
    with open(a.src, "rb") as f:
        while True:
            raw = np.fromfile(f, np.int16, 2 * a.chunk)
            if not raw.size:
                break
            x = raw.astype(np.float64)
            acc += float(np.dot(x, x))
    p_rec = acc / ntot                       # per complex sample (I^2 + Q^2)
    sigma = np.sqrt(p_rec * 10 ** (-a.db / 10) / 2)   # per component
    rng = np.random.default_rng(a.seed)
    peak = 0
    n_clip = 0
    with open(a.src, "rb") as f, open(a.dst, "wb") as g:
        while True:
            raw = np.fromfile(f, np.int16, 2 * a.chunk)
            if not raw.size:
                break
            y = raw.astype(np.float64) + rng.normal(0.0, sigma, raw.size)
            yi = np.rint(y)
            n_clip += int(np.count_nonzero((yi > 32767) | (yi < -32768)))
            yi = np.clip(yi, -32768, 32767).astype(np.int16)
            peak = max(peak, int(np.abs(yi).max()))
            yi.tofile(g)
    side = dict(source=os.path.abspath(a.src), noise_db_below_recorded=a.db,
                seed=a.seed, sigma_per_component_int16=round(float(sigma), 3),
                recorded_power_int16sq=round(p_rec, 1), samples=int(ntot),
                peak_abs=peak, clipped_samples=n_clip)
    # carry the source's sidecar (rate / center) if it has one
    sj = os.path.splitext(a.src)[0] + ".json"
    if os.path.exists(sj):
        try:
            side["source_sidecar"] = json.load(open(sj))
        except Exception:                                      # noqa: BLE001
            pass
    json.dump(side, open(os.path.splitext(a.dst)[0] + ".json", "w"), indent=1)
    print(f"{a.dst}: {ntot} samples, noise {a.db} dB below recorded power "
          f"(sigma {sigma:.1f} per component), peak |s| {peak} "
          f"({20*np.log10(peak/32768):+.1f} dBFS), clipped {n_clip}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
