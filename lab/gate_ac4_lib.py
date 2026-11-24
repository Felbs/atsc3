#!/usr/bin/env python3
"""GATE: the audio worker's libac4dec path (lab/m45_ac4_lib.py).

On a recorded AC-4 lane (5_X element), with the library found as the worker
finds it ($ATSC3_AC4DEC_LIB or ../ac4-decoder/c/libac4dec.so):

  1. EVERY FRAME OUT ONCE: frames fed == frames emitted after the final
     drain, 1536 samples each, and nothing concealed on a clean lane.
  2. PASSES DO NOT MATTER: the lane fed in pieces of 1, 7 and 60 frames
     gives the SAME samples as one call (np.array_equal) -- the live worker
     feeds whatever a pass read.
  3. STEREO CARRIES THE DIALOGUE: the 2-channel output of a 5_X element,
     regressed on the 6-channel decode's L R C Ls Rs, gives the centre a
     weight of at least 0.5 in both channels (LoRo: 0.707). This is the
     defect of d6b0b9c, lost in the 11-01 snapshot (9d46747).
  4. ON TIME: the 6-channel output's L and C are not shifted against each
     other (the old path's centre led its fronts by 577 samples).
NEGATIVE CONTROL: the L/R pair of the 6-channel decode -- what --channels 2
wrote before -- must FAIL check 3.

    python lab/gate_ac4_lib.py LANE.m4s [--frames N]
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "tools"))

FAILS = []


def check(name, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


def centre_weight(st, six):
    """median over 4 s windows of the C weight in st[:, 0] and st[:, 1]."""
    W = 4 * 46034
    A6 = np.column_stack([six[k] for k in ("L", "R", "C", "Ls", "Rs")])
    w = []
    for s in range(0, min(len(A6), len(st)) - W, W):
        A = A6[s:s + W]
        if np.std(A[:, 2]) < 1.0:
            continue
        for c in range(2):
            co, *_ = np.linalg.lstsq(A, st[s:s + W, c], rcond=None)
            w.append((c, co[2]))
    if not w:
        return None
    return [float(np.median([v for cc, v in w if cc == c])) for c in range(2)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("lane")
    ap.add_argument("--frames", type=int, default=1800)
    a = ap.parse_args()
    try:
        import m45_ac4_lib as M
        M.load_lib()
    except Exception as e:                                     # noqa: BLE001
        print(f"GATE INVALID -- libac4dec not loadable: {e}")
        return 2
    from atsc3_audio import LaneReader
    frames = LaneReader(a.lane).read_new()[:a.frames]
    if len(frames) < 300:
        print(f"GATE INVALID -- {len(frames)} frames in {a.lane}")
        return 2
    n = len(frames)
    print(f"lane {a.lane}: {n} frames, library {M.load_lib()._atsc3_path}")

    d = M.LibAc4Stream("5_X", 6)
    one = d.push(frames, final=True)
    st = d.stats()
    check("every frame out once", all(len(v) == n * 1536 for v in one.values()),
          f"{n} in, {len(one['L']) / 1536:.0f} out")
    check("clean lane: nothing concealed", st["n_concealed"] == 0,
          f"{st['n_concealed']} concealed {d.reasons()}")

    for step in (1, 7, 60):
        d = M.LibAc4Stream("5_X", 6)
        parts = [d.push(frames[i:i + step]) for i in range(0, n, step)]
        parts.append(d.push([], final=True))
        cat = {k: np.concatenate([p[k] for p in parts]) for k in one}
        check(f"pieces of {step} == one call",
              all(np.array_equal(cat[k], one[k]) for k in one))

    d = M.LibAc4Stream("5_X", 2)
    two = d.push(frames, final=True)
    st2 = np.column_stack([two["L"], two["R"]])
    w = centre_weight(st2, one)
    check("stereo carries the centre", w is not None and min(w) >= 0.5,
          f"C weight L {w[0]:.3f} R {w[1]:.3f}" if w else "no centre signal")
    wn = centre_weight(np.column_stack([one["L"], one["R"]]), one)
    check("negative control: L/R alone fails the centre check",
          wn is not None and min(wn) < 0.5,
          f"C weight L {wn[0]:.3f} R {wn[1]:.3f}" if wn else "no centre signal")

    # alignment: the lag of best L-C cross-correlation over the decode
    x, y = one["L"], one["C"]
    m = 2048
    nf = 1 << (len(x) + len(y)).bit_length()
    c = np.fft.irfft(np.fft.rfft(x, nf) * np.conj(np.fft.rfft(y, nf)), nf)
    lags = np.r_[np.arange(0, m + 1), np.arange(-m, 0)]
    lag = int(lags[np.argmax(np.abs(np.r_[c[:m + 1], c[-m:]]))])
    check("L and C on time (|lag| < 100 samples)", abs(lag) < 100, f"lag {lag}")

    print()
    if FAILS:
        print(f"GATE FAIL -- {FAILS}")
        return 1
    print("GATE PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
