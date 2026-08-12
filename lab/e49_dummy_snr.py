#!/usr/bin/env python3
"""E49 instrument -- channel SNR off the 51257 dummy cells (A/322 7.2.6.5).

The dummy cells are transmitted as +-1 (real) with the 5.2.3 scrambler's
signs, and they sit OUTSIDE the LDM region, so fitting the known sequence and
measuring the residual gives the CHANNEL SNR directly -- no constellation
hypothesis, no decode needed.  M10 measured 3.37 dB on the old RF25 capture
with this method; reproducing that number here is the calibration gate.

    a     = <Re(d), s> / <s, s>          (real LS gain onto the known signs)
    noise = mean |d - a*s|^2             (both dimensions)
    SNR   = a^2 / noise
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

import m4_scrambler as SC            # noqa: E402
import m10_core as M10               # noqa: E402


def dummy_snr(pool, core_total):
    d = pool[core_total:]
    s = 1.0 - 2.0 * SC.sequence(len(pool))[core_total:]
    a = float(np.dot(d.real, s) / np.dot(s, s))
    noise = float(np.mean(np.abs(d - a * s) ** 2))
    snr = a * a / noise
    return dict(n=len(d), gain=a, noise=noise,
                snr_db=10 * np.log10(snr))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("capture")
    ap.add_argument("--rate", type=float, required=True)
    ap.add_argument("--l1", required=True)
    ap.add_argument("--frames", type=int, default=3)
    ap.add_argument("--json")
    a = ap.parse_args()

    js = json.load(open(a.l1))
    g = M10.geometry_from_json(js, os.path.basename(a.capture))
    core_total = sum(p["size"] for p in g.plps if p["layer"] == 0)
    per = M10.frame_samples(g) / M10.FS_POST
    plps = {p["id"]: p for p in g.plps}

    rows = []
    for n in range(a.frames):
        s = max(0.0, n * per - 0.01)
        rep = {}
        pool, info, _ = M10.frame_cells(a.capture, a.rate, None, g, s,
                                        plps[0], report=rep)
        pool = M10.cpe(pool, info["symbol_of"], g, plps[0],
                       dummy_start=core_total)
        r = dummy_snr(pool, core_total)
        r["frame"] = n
        rows.append(r)
        print(f"  frame {n}  gain {r['gain']:+.4f}  noise {r['noise']:.4f}  "
              f"SNR {r['snr_db']:.2f} dB   ({r['n']} dummy cells)")

    med = float(np.median([r["snr_db"] for r in rows]))
    print(f"  median channel SNR: {med:.2f} dB")
    if a.json:
        json.dump(dict(capture=a.capture, rows=rows, median_snr_db=med),
                  open(a.json, "w"), indent=2, default=float)
        print(f"  wrote {a.json}")


if __name__ == "__main__":
    main()
