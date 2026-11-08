#!/usr/bin/env python3
"""E133 -- probes on the banked failing Frames: which preparation the standalone loader does that
the live window lacks.  For each window: the live result; the start shifted by a few samples; the
carrier offset re-measured on the window's OWN bootstrap (as m2_pilots.load does) and removed before
the decoder sees it; and that plus the band limit."""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
from scipy.signal import firwin, lfilter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
import m6_cells as C                                             # noqa: E402
import m9_fast                                                   # noqa: E402
import m13_sf1 as M13                                            # noqa: E402
import m17_plp1_live as M17                                      # noqa: E402
from atsc3 import bootstrap as bs                                # noqa: E402
from m2_stage_a import find_bootstraps                           # noqa: E402
from m2_stage_b import resample_to                               # noqa: E402


def notch(y, fs):
    nt = 257
    b = firwin(nt, 2.95e6, fs=fs)
    return lfilter(b, 1.0, np.concatenate([y, np.zeros(nt // 2, y.dtype)]))[nt // 2:]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", default=r"C:\captures\sf1_fail")
    ap.add_argument("--l1", default=os.path.join(HERE, "m4_l1detail_hit_rf33.json"))
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--limit", type=int, default=4)
    a = ap.parse_args(argv)
    files = sorted(glob.glob(os.path.join(a.dump, "*.npz")))[: a.limit]
    js = json.load(open(a.l1))
    g0 = C.Geometry.rf33_legacy()
    fd = m9_fast.FrameDecoder(threads=a.threads, backend="cpu", iters=50, cpu_fast=True)
    fd.prewarm()
    dec = M17.SubframeDecoder(fd, g0, js["l1_basic"], js["l1_detail"], 1, comb=True)
    dec.prewarm()
    fs = M13.FS_POST

    def run(y, t0):
        rep = {}
        try:
            dec.decode(y, t0, rep)
            return int(rep["converged"]), rep.get("bad_fail", []), rep.get("cfo_resid_hz", 0.0), rep.get("ldpc_iters_mean", 0.0)
        except Exception as ex:                                      # noqa: BLE001
            return -1, repr(ex)[:50], 0.0, 0.0

    for f in files:
        z = np.load(f)
        y = np.asarray(z["y"], np.complex128)
        t0 = int(z["t0"])
        print(f"\nframe {int(z['idx'])}: live said {int(z['converged'])}/{int(z['n_fec'])}, t0 {t0}, window {len(y)} samples")
        r = run(y, t0)
        print(f"  live window               conv {r[0]:3d}  bad {r[1]}  m17 residual {r[2]:+.2f} Hz  its {r[3]:.1f}")
        for sh in (-40, -20, +20, +40, +100):
            r = run(y, t0 + sh)
            print(f"  t0 {sh:+4d}                   conv {r[0]:3d}  bad {r[1]}")
        # the window's own bootstrap: where it sits and what carrier offset it reads
        head = y[: int(0.05 * fs)]
        hits = find_bootstraps(resample_to(head, fs, bs.FS))
        if hits:
            h = hits[0]
            pos = int(round(h["position"] * fs / bs.FS))
            cfo = float(h["fine_cfo_hz"])
            print(f"  window's bootstrap: position {pos} (live t0 {t0}), fine CFO {cfo:+.2f} Hz")
            yc = y * np.exp(-2j * np.pi * cfo * np.arange(len(y)) / fs)
            r = run(yc, t0)
            print(f"  bootstrap CFO removed     conv {r[0]:3d}  bad {r[1]}  m17 residual after {r[2]:+.2f} Hz  its {r[3]:.1f}")
            r = run(yc, pos)
            print(f"  + bootstrap position      conv {r[0]:3d}  bad {r[1]}")
            r = run(notch(yc, fs), pos)
            print(f"  + band limit              conv {r[0]:3d}  bad {r[1]}  its {r[3]:.1f}")
        else:
            print("  no bootstrap found in the window head")
    return 0


if __name__ == "__main__":
    sys.exit(main())
