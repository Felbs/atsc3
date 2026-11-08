#!/usr/bin/env python3
"""gate_e119_l1_kernel.py -- L1-Basic/Detail LDPC through the compiled kernel: the referee.

`m8_l1.decode` ran the per-block numpy float64 min-sum; E119 offers the
compiled flooding kernel (float32, same alpha) behind ATSC3_L1_KERNEL.  L1
carries its own BCH/CRC, so the law is: on every Frame window where either
path verifies, BOTH verify and the decoded L1 is field-for-field identical;
the kernel path verifies on no fewer windows; and it is faster.  Windows come
off the real bootstrap grid of banked captures -- clean air and the noisy
copies (lab/add_noise_cs16.py) where the cost lives.

Run: python lab/gate_e119_l1_kernel.py --captures A.cs16 B.cs16 ... [--frames 12]
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
import m2_pilots as MP                                             # noqa: E402
import m3_ldpc as LD                                               # noqa: E402
import m44_ldm as M44                                              # noqa: E402
from atsc3 import bootstrap as bs                                  # noqa: E402

RATE = 6.912e6


def load_frames(path, n, start_sec=0.0):
    """n consecutive Frame windows, de-rotated, off the real bootstrap grid
    (gate_e82's loader)."""
    fmt = MP.guess_format(path)
    head = MP.read_block(path, int(start_sec * RATE), int(0.30 * RATE), fmt)
    h = MP.find_bootstraps(MP.resample_to(head, RATE, bs.FS))[0]
    f0 = int(start_sec * RATE) + int(round(h["position"] * RATE / bs.FS))
    ps = h["fields"].get("preamble_structure")
    probe = MP.read_block(path, f0 + 1_600_000, int(0.35 * RATE), fmt)
    h2 = MP.find_bootstraps(MP.resample_to(probe, RATE, bs.FS))[0]
    step = 1_600_000 + int(round(h2["position"] * RATE / bs.FS))
    x = MP.read_block(path, f0, step * n + 500_000, fmt)
    y = x.astype(np.complex128)
    y *= np.exp(-2j * np.pi * h["fine_cfo_hz"] * np.arange(len(y)) / RATE)
    return y, step, ps


def canon(L):
    """The decoded L1 as comparable data (no timings, no floats that are not fields)."""
    if L is None:
        return None
    def fix(v):
        if isinstance(v, dict):
            return {k: fix(x) for k, x in sorted(v.items()) if k not in ("ms", "time", "elapsed")}
        if isinstance(v, (list, tuple)):
            return [fix(x) for x in v]
        if isinstance(v, np.ndarray):
            return v.tolist()
        if isinstance(v, (np.integer,)):
            return int(v)
        if isinstance(v, (np.floating,)):
            return float(v)
        return v
    return fix(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--captures", nargs="+", required=True)
    ap.add_argument("--frames", type=int, default=12)
    ap.add_argument("--start", type=float, default=1.0)
    a = ap.parse_args()
    tot_bad = 0
    n_off = n_on = n_win = 0
    t_off = t_on = 0.0
    for cap in a.captures:
        y, step, ps = load_frames(cap, a.frames, a.start)
        print(f"\n  {os.path.basename(cap)}: {a.frames} Frames, step {step}, ps {ps}")
        for k in range(a.frames):
            w = y[k * step:k * step + 400_000]
            t0, _c = M44.plateau_t0(w, M44.BOOTSTRAP, ps)
            LD.L1_KERNEL = False
            t = time.perf_counter()
            L_off = M44.l1_from_window(w, t0, ps=ps)
            t_off += time.perf_counter() - t
            LD.L1_KERNEL = True
            t = time.perf_counter()
            L_on = M44.l1_from_window(w, t0, ps=ps)
            t_on += time.perf_counter() - t
            LD.L1_KERNEL = False
            ok_off = bool(L_off and L_off.get("ok"))
            ok_on = bool(L_on and L_on.get("ok"))
            n_win += 1; n_off += ok_off; n_on += ok_on
            same = canon(L_off) == canon(L_on)
            bad = (ok_off != ok_on) or (ok_off and not same)
            tot_bad += bad
            print(f"    Frame {k:2d}: numpy {'ok ' if ok_off else '-- '} kernel {'ok ' if ok_on else '-- '} "
                  f"fields {'identical' if same else 'DIFFER'}{'   <-- ' if bad else ''}")
    print(f"\n  windows {n_win}: numpy verified {n_off}, kernel verified {n_on}; "
          f"L1 time numpy {t_off*1000:.0f} ms, kernel {t_on*1000:.0f} ms ({t_off/max(t_on,1e-9):.1f}x)")
    ok = tot_bad == 0 and n_on >= n_off and t_on <= t_off
    print("  PASS: identical L1 wherever either verifies, kernel verifies no fewer, faster"
          if ok else f"  FAIL: {tot_bad} windows differ or verify on one path only"
          + ("" if n_on >= n_off else "; kernel verifies fewer") + ("" if t_on <= t_off else "; slower"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
