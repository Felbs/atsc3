#!/usr/bin/env python3
"""gate_e59.py -- in-position padding + clock-servo plausibility.

Two mechanisms from the 8/09 e31 buzz/desync session:
  1. batch-END padding displaced content after any in-batch hole (E51's
     recorded open item) -- pads must land AT the truncated fragment;
  2. the E50 clock servo acted on one burst-contaminated window and
     slammed its +-1500 ppm clamp, time-stretching VLC 0.15% for the
     rest of the session -- implausible fits must be discarded and the
     first engagement must require two agreeing windows.
Negative controls run the OLD behaviour on the same inputs and fail.

    python lab/gate_e59.py    -> PASS/FAIL, exit 0 iff all pass
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tools"))

from atsc3_audio import pad_insertions                          # noqa: E402
import atsc3_tv as tv                                           # noqa: E402

FAILS = []


def check(name, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" +
          (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


# ---------------------------------------------------------------- gate 1
# assembly simulation: fragment k's frames all carry value k+1; zeros are
# pads. After assembly every fragment's content must sit at its NOMINAL
# 60-frame-grid position.

def assemble(fcounts, ins_fn):
    """Build the batch PCM (1 sample per frame for speed) via `ins_fn`."""
    frames = []
    for k, c in enumerate(fcounts):
        frames += [k + 1] * c
    y = np.array(frames, dtype=float)
    ins = ins_fn(fcounts, 0, 0, len(frames))
    parts, prev = [], 0
    for pos, miss in ins:
        parts.append(y[prev:pos])
        parts.append(np.zeros(miss))
        prev = pos
    parts.append(y[prev:])
    return np.concatenate(parts)


def gate_in_position():
    print("gate 1: pads land at their fragment's position")
    fcounts = [60, 60, 20, 60, 45, 60]          # holes in frags 2 and 4
    out = assemble(fcounts, pad_insertions)
    ok = len(out) == 6 * 60
    for k, c in enumerate(fcounts):
        seg = out[k * 60:(k + 1) * 60]
        ok = ok and (seg[:c] == k + 1).all() and (seg[c:] == 0).all()
    check("every fragment on its nominal 60-grid slot, holes = silence",
          bool(ok), f"len {len(out)}")

    # NEGATIVE CONTROL -- the shipped batch-END padding on the same lane:
    # everything after frag 2's hole is displaced EARLY by the hole size.
    frames = []
    for k, c in enumerate(fcounts):
        frames += [k + 1] * c
    y_old = np.concatenate([np.array(frames, dtype=float),
                            np.zeros(6 * 60 - len(frames))])
    disp3 = 60 * 3 - int(np.argmax(y_old == 4))   # frag 3 nominal - actual
    check("NEGATIVE CONTROL: batch-end padding displaces post-hole content",
          disp3 == 40 and not (y_old == out).all(),
          f"frag 3 early by {disp3} frames (= frag 2's hole)")


# ---------------------------------------------------------------- gate 2
# servo: contaminated window discarded; engagement needs confirmation.

def feed(ct, t0, dur, cushion_fn):
    r = None
    t = t0
    while t < t0 + dur:
        ct.sample(t, cushion_fn(t))
        got = ct.evaluate(t)
        if got is not None:
            r = got
        t += 5.0
    return r, t


def gate_servo():
    print("gate 2: clock servo -- plausibility + confirmation")
    # contaminated window: cushion ramps +0.5 s/s (an append burst)
    ct = tv.ClockTrim(target=60.0)
    r, t = feed(ct, 0.0, 700.0, lambda t: 60.0 + 0.5 * t)
    span = ct.samples[-1][0] - ct.samples[0][0] if ct.samples else 0.0
    check("burst-contaminated window -> NO trim, window discarded",
          r is None and ct.rate == 1.0 and span < 650.0,
          f"rate {ct.rate}, residual window span {span:.0f}s")

    # NEGATIVE CONTROL -- the shipped math on the same window: slope
    # 0.5 s/s -> adj >> clamp -> rate slams 1.0015.
    slope = 0.5
    adj = slope + (60.0 + 0.5 * 700 - 60.0) / 3600.0
    old = max(1.0 - 1500e-6, min(1.0 + 1500e-6, 1.0 + adj))
    check("NEGATIVE CONTROL: old math slams the 1500 ppm clamp",
          old == 1.0015, f"old rate {old}")

    # honest 120 ppm skew: first window -> held pending; second agreeing
    # window -> engages once, near 120 ppm, far from the clamp
    ct = tv.ClockTrim(target=60.0)
    r1, t = feed(ct, 0.0, 700.0, lambda t: 60.0 + 120e-6 * t)
    r2, t = feed(ct, t, 700.0, lambda t: 60.0 + 120e-6 * t)
    engaged = r2 if r2 is not None else r1
    check("real 120 ppm skew engages on the SECOND agreeing window",
          r1 is None and engaged is not None
          and abs((engaged - 1.0) * 1e6 - 120) < 60,
          f"first {r1} second {engaged}")

    ct.hard_reset()
    check("hard_reset returns rate to the new player's 1.0",
          ct.rate == 1.0 and ct.pending is None and not ct.samples)


if __name__ == "__main__":
    gate_in_position()
    gate_servo()
    print(f"\n{'ALL GATES PASS' if not FAILS else 'FAILED: ' + ', '.join(FAILS)}")
    sys.exit(1 if FAILS else 0)
