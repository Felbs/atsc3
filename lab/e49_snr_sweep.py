#!/usr/bin/env python3
"""E49 -- sweep the dummy-cell channel SNR across a whole capture.

One frame every --stride seconds.  Finds the windows (if any) where the
channel clears the composite threshold, and whether the variability is
fading (slow, dB-scale) or glitches (isolated collapses).
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m10_core as M10               # noqa: E402
from e49_dummy_snr import dummy_snr  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("capture")
    ap.add_argument("--rate", type=float, required=True)
    ap.add_argument("--l1", required=True)
    ap.add_argument("--stride", type=float, default=3.0)
    ap.add_argument("--t-max", type=float, default=119.0)
    ap.add_argument("--json")
    a = ap.parse_args()

    js = json.load(open(a.l1))
    g = M10.geometry_from_json(js, os.path.basename(a.capture))
    core_total = sum(p["size"] for p in g.plps if p["layer"] == 0)
    plps = {p["id"]: p for p in g.plps}

    rows = []
    t = 0.0
    while t <= a.t_max:
        try:
            rep = {}
            pool, info, _ = M10.frame_cells(a.capture, a.rate, None, g, t,
                                            plps[0], report=rep)
            pool = M10.cpe(pool, info["symbol_of"], g, plps[0],
                           dummy_start=core_total)
            r = dummy_snr(pool, core_total)
            r["t"] = t
            r["data_coh"] = rep["data_coherence_mean"]
            rows.append(r)
            print(f"  t={t:7.2f}s  SNR {r['snr_db']:6.2f} dB  "
                  f"data-coh {rep['data_coherence_mean']:.3f}")
        except Exception as exc:                       # noqa: BLE001
            print(f"  t={t:7.2f}s  FAILED ({type(exc).__name__}: {exc})")
            rows.append(dict(t=t, snr_db=None))
        t += a.stride

    good = [r for r in rows if r.get("snr_db") is not None]
    if good:
        v = np.array([r["snr_db"] for r in good])
        print(f"\n  {len(good)} frames: median {np.median(v):.2f} dB, "
              f"min {v.min():.2f}, max {v.max():.2f}, "
              f">=7.0 dB: {(v >= 7.0).sum()}, >=5.0 dB: {(v >= 5.0).sum()}")
    if a.json:
        json.dump(rows, open(a.json, "w"), indent=2, default=float)
        print(f"  wrote {a.json}")


if __name__ == "__main__":
    main()
