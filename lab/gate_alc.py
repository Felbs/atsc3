#!/usr/bin/env python3
"""GATE: automatic level control (lab/m46_alc.py) and its use in the audio
worker (tools/atsc3_audio.py --alc on / ATSC3_ALC=1).

SYNTHETIC (always run, no captures needed):
  S1 METER: BS.1770 reference points -- a 997 Hz sine at -23 dBFS on L and R
     reads -23.0 LUFS (EBU Tech 3341 case 1), and in the 5.1 weighting the
     LFE does not count and a surround counts 1.41x.
  S2 PASSES DO NOT MATTER: the same 60 s of 5.1 fed whole and in random
     pieces (1 sample .. 3 s) gives the same output (|diff| < 1e-9).
  S3 NO PUMPING: speech-like bursts at -40 LUFS with 2 s pauses: converged
     to +16 dB within 3 s; across the pauses the gain moves < 0.2 dB.
  S4 CAPS: a -60 LUFS programme gets exactly +max_boost, a -1 LUFS one
     exactly -max_cut.
  S5 LIMITER: a high-crest programme pushed to -8 LUFS never exceeds the
     ceiling; nothing reaches full scale.
  S6 BALANCE: every channel gets the same gain, sample for sample.
  S7 COST: at least 50x real time for 5.1 at 48 kHz on one core.

LANES (optional, recorded AC-4 lanes as banked live dirs):
    python lab/gate_alc.py --lane DIR:ELEMENT:CHANNELS[:REF_PEAK] ...
  runs the worker with ALC off and on (--once), and once more on with
  --max-pass 50 (the multi-pass replay of a live worker), and checks:
  L1 the same number of samples as with ALC off;
  L2 integrated loudness (BS.1770, after the first 10 s) within +-1 LU of
     the target, unless the gain sat on a cap;
  L3 converged in <= 5 s: the gain is within 1 dB of target - (the input's
     integrated loudness over its first 3 s of programme), counted from the
     first programme block;
  L4 no clipping: no full-scale sample, sample peak <= ceiling, true peak
     (4x oversampled) <= 0 dBTP;
  L5 a smooth trace: after the fast start no 100 ms gain step above the
     middle slew (0.4 dB), and none above 1.5 dB at all;
  L6 balance: sample by sample, every channel was multiplied by the same
     gain (on/off ratio spread across channels < 0.05 dB where |off| is
     above 2000 LSB);
  L7 the ALC adds no pass dependence: wherever the ALC-off output of one
     pass and of --max-pass 50 are byte-identical, the ALC-on outputs are
     too. (Where they are not, the worker's REF_PEAK rule moved at a pass
     boundary -- old behaviour, reported, not judged; S2 covers the ALC.)
  A REF_PEAK (e.g. 1.8e6 on a stereo pair) makes a lane quiet on purpose,
  to test the boost side the way broadcaster A's stereo needs it.

    python lab/gate_alc.py                      # synthetic only
    python lab/gate_alc.py --lane ~/lanes/a51:5_X:6 --lane ~/lanes/ap:pair:2
"""
from __future__ import annotations

import argparse
import math
import os
import subprocess
import sys
import tempfile
import time
import wave

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

from m46_alc import (Alc, FS, loudness_integrated,               # noqa: E402
                     true_peak_dbtp)

FAILS = []


def check(name, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}"
          + (f"  ({detail})" if detail else ""), flush=True)
    if not ok:
        FAILS.append(name)


def speechish(n_s, ch, rng, pauses=None, level=-30.0):
    """Band-limited noise with a syllable envelope; optional pause pattern
    (seconds on, seconds off). Scaled to `level` LUFS."""
    n = int(n_s * FS)
    from scipy.signal import butter, sosfilt
    sos = butter(4, [150, 4000], btype="band", fs=FS, output="sos")
    x = sosfilt(sos, rng.standard_normal((n, ch)), axis=0)
    t = np.arange(n) / FS
    env = 0.55 + 0.45 * np.sin(2 * np.pi * 4.0 * t + rng.uniform(0, 6))
    x *= env[:, None]
    if pauses:
        on, off = pauses
        x[(t % (on + off)) >= on] = 0.0
    if ch == 6:
        x[:, 3] *= 0.2
    li = loudness_integrated(x)
    return x * 10 ** ((level - li) / 20)


def run_synthetic():
    rng = np.random.default_rng(1770)
    print("synthetic:")
    # S1
    t = np.arange(20 * FS) / FS
    s = 10 ** (-23 / 20) * np.sin(2 * np.pi * 997 * t)
    l2 = loudness_integrated(np.column_stack([s, s]))
    z = np.zeros_like(s)
    l_lfe = loudness_integrated(np.column_stack([s, s, z, s, z, z]))
    l_sur = loudness_integrated(np.column_stack([z, z, z, z, s, z]))
    l_one = loudness_integrated(np.column_stack([s, z, z, z, z, z]))
    check("S1 meter: 997 Hz at -23 dBFS on L+R = -23.0 LUFS",
          abs(l2 + 23) < 0.1, f"{l2:.2f}")
    check("S1 meter: LFE excluded, surround x1.41",
          abs(l_lfe - l2) < 0.05 and abs(l_sur - l_one - 10 * math.log10(1.41)) < 0.05,
          f"L+R+LFE {l_lfe:.2f}, Ls-L {l_sur - l_one:.2f} dB")
    # S2
    x = speechish(60, 6, rng, pauses=(4, 1), level=-36)
    x[30 * FS:] *= 4.0                      # a step mid-way
    whole = Alc(6).process(x)
    a = Alc(6)
    parts, i = [], 0
    while i < len(x):
        k = int(rng.integers(1, 3 * FS))
        parts.append(a.process(x[i:i + k]))
        i += k
    d = float(np.abs(np.vstack(parts) - whole).max())
    check("S2 passes do not matter", d < 1e-9, f"max |diff| {d:.2e}")
    # S3
    x = speechish(40, 2, rng, pauses=(3, 2), level=-40)
    a = Alc(2)
    g, lv = [], []
    for i in range(0, len(x), 480):           # 10 ms steps
        a.process(x[i:i + 480])
        g.append(a.gain_db)
    g = np.array(g)
    tt = np.arange(len(g)) * 0.01
    conv = next((t_ for t_, v in zip(tt, g) if abs(v - 16.0) <= 1.0), None)
    late = tt >= 10
    t_mod = tt % 5.0
    pause = late & (t_mod >= 3.2) & (t_mod < 5.0)
    move = 0.0
    for c0 in np.arange(10, 40, 5.0):
        m = pause & (tt >= c0) & (tt < c0 + 5)
        if m.any():
            move = max(move, float(np.ptp(g[m])))
    check("S3 converges fast", conv is not None and conv <= 3.0,
          f"within 1 dB of +16 at {conv} s")
    check("S3 no pumping in pauses", move < 0.2,
          f"largest gain move inside a pause {move:.3f} dB; "
          f"final {g[-1]:+.2f} dB")
    # S4
    a = Alc(2)
    a.process(speechish(10, 2, rng, level=-60))
    b = Alc(2)
    b.process(speechish(10, 2, rng, level=-1))
    check("S4 caps", abs(a.gain_db - 20) < 1e-9 and abs(b.gain_db + 20) < 1e-9,
          f"{a.gain_db:+.2f} / {b.gain_db:+.2f} dB")
    # S5
    x = speechish(30, 6, rng, level=-30)
    clicks = rng.integers(0, len(x), 200)
    x[clicks] *= 25.0                       # crest factor way above programme
    a = Alc(6, target=-8.0)
    y = a.process(x)
    pk = float(np.abs(y).max())
    check("S5 limiter holds the ceiling", pk <= 10 ** (-1 / 20) + 1e-12
          and a.n_lim > 0,
          f"peak {20 * math.log10(pk):.3f} dBFS, {a.n_lim} samples limited, "
          f"deepest {a.lim_min:.1f} dB")
    # S6
    nz = np.abs(x) > 1e-6
    r = np.where(nz, y / np.where(nz, x, 1), np.nan)
    spread = np.nanmax(r, axis=1) - np.nanmin(r, axis=1)
    check("S6 one gain for all channels", float(np.nanmax(spread)) < 1e-9,
          f"max spread {float(np.nanmax(spread)):.1e}")
    # S7
    x = speechish(60, 6, rng, level=-30)
    a = Alc(6)
    t0 = time.perf_counter()
    for i in range(0, len(x), 2 * FS):
        a.process(x[i:i + 2 * FS])
    dt = time.perf_counter() - t0
    check("S7 cost", 60 / dt >= 50, f"{60 / dt:.0f}x real time, 5.1")


def read_wav(p):
    w = wave.open(p)
    x = np.frombuffer(w.readframes(w.getnframes()), "<i2")
    return x.reshape(-1, w.getnchannels()).astype(float) / 32768.0


def run_lane(spec, tmp, py):
    parts = spec.split(":")
    d, el, ch = os.path.expanduser(parts[0]), parts[1], int(parts[2])
    ref = parts[3] if len(parts) > 3 else None
    tag = os.path.basename(d.rstrip("/")) + f"_{el}_{ch}" + (f"_r{ref}" if ref else "")
    print(f"lane {tag}:")
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    env.pop("ATSC3_ALC", None)
    tr = os.path.join(tmp, tag + ".trace")
    if os.path.exists(tr):
        os.remove(tr)

    def worker(name, alc, extra=(), trace=None):
        out = os.path.join(tmp, f"{tag}_{name}.wav")
        for suf in (".wav", ".state.json", ".json"):
            p = os.path.splitext(out)[0] + suf
            if os.path.exists(p):
                os.remove(p)
        e = dict(env)
        if trace:
            e["ATSC3_ALC_TRACE"] = trace
        cmd = [py, os.path.join(ROOT, "tools", "atsc3_audio.py"), "--once",
               "--wait", "5", "--live-dir", d, "--element", el,
               "--channels", str(ch), "--out", out, "--alc", alc] + list(extra)
        if ref:
            cmd += ["--ref-peak", ref]
        r = subprocess.run(cmd, env=e, capture_output=True, text=True)
        if r.returncode:
            print(r.stdout[-2000:], r.stderr[-2000:])
            raise SystemExit(f"worker failed on {tag}")
        return read_wav(out), r.stdout
    x_off, _ = worker("off", "off")
    x_on, log_on = worker("on", "on", trace=tr)
    x_mp, _ = worker("mp", "on", ["--max-pass", "50"])
    a = Alc(ch)
    tgt = a.target
    check("L1 same length", len(x_on) == len(x_off) == len(x_mp),
          f"{len(x_on)} samples")
    T = np.loadtxt(tr, ndmin=2)              # t, block, level, gain, lim
    tt, gg = T[:, 0], T[:, 3]
    sk = 10 * FS
    i_in = loudness_integrated(x_off[sk:])
    i_on = loudness_integrated(x_on[sk:])
    i_on_all = loudness_integrated(x_on)
    capped = (tgt - i_in) > a.max_boost or (tgt - i_in) < -a.max_cut
    check("L2 integrated loudness on target",
          capped or abs(i_on - tgt) <= 1.0,
          f"in {i_in:.1f} -> out {i_on:.1f} LUFS after 10 s "
          f"({i_on_all:.1f} whole){' [cap]' if capped else ''}")
    # the input's own level over its first 10 s of programme
    g10 = None
    lev = T[:, 2]
    fin = np.isfinite(lev)
    if fin.any():
        t_first = tt[fin][0]
        seg = x_off[int(max(0, t_first - 0.4) * FS):int((t_first + 2.6) * FS)]
        g10 = min(a.max_boost, max(-a.max_cut, tgt - loudness_integrated(seg)))
        conv = next((t_ for t_, v in zip(tt, gg) if abs(v - g10) <= 1.0), None)
        conv_rel = None if conv is None else conv - t_first + 0.4
    else:
        conv_rel = None
    check("L3 converged within 5 s of programme",
          conv_rel is not None and conv_rel <= 5.0,
          f"{conv_rel if conv_rel is None else round(conv_rel, 1)} s "
          f"to within 1 dB of {g10 if g10 is None else round(g10, 1)} dB")
    pk = float(np.abs(x_on).max())
    full = int((np.abs(x_on) >= 32767 / 32768).sum())
    tp = true_peak_dbtp(x_on)
    check("L4 no clipping", full == 0 and pk <= 10 ** (-1 / 20) + 1 / 32768
          and tp <= 0.0,
          f"peak {20 * math.log10(pk + 1e-12):.2f} dBFS, true peak "
          f"{tp:.2f} dBTP, {full} full-scale; limiter: "
          + (log_on.split("limited")[1].split(")")[0] + ")"
             if "limited" in log_on else "never fired"))
    dg = np.abs(np.diff(gg))
    after = tt[1:] > (tt[fin][0] + 3.5 if fin.any() else 0)
    big = float(dg.max()) if len(dg) else 0.0
    slow = float(dg[after].max()) if after.any() else 0.0
    check("L5 smooth gain trace", big <= 1.5 + 1e-6 and slow <= 0.4 + 1e-6,
          f"largest 100 ms step {big:.2f} dB, after the start {slow:.3f} dB; "
          f"gain {gg.min():+.1f}..{gg.max():+.1f} dB, final {gg[-1]:+.1f}")
    big_ = np.abs(x_off) > 2000 / 32768
    rows = big_.sum(axis=1) >= 2
    if rows.any():
        r = np.where(big_[rows], 20 * np.log10(
            np.abs(x_on[rows]) / np.where(big_[rows], np.abs(x_off[rows]), 1)
            + 1e-12), np.nan)
        spread = float(np.nanmax(np.nanmax(r, axis=1) - np.nanmin(r, axis=1)))
        check("L6 balance kept", spread < 0.05,
              f"largest on/off spread across channels {spread:.4f} dB "
              f"over {int(rows.sum())} samples")
    else:
        check("L6 balance kept", True, "no sample loud in two channels")
    x_off_mp, _ = worker("offmp", "off", ["--max-pass", "50"])
    if np.array_equal(x_off, x_off_mp):
        check("L7 no pass dependence added", np.array_equal(x_on, x_mp),
              "ALC on: one pass == --max-pass 50"
              if np.array_equal(x_on, x_mp) else
              f"differ by up to {np.abs(x_on - x_mp).max() * 32768:.0f} LSB")
    else:
        print("  ----  L7 not judged: with ALC OFF one pass and --max-pass 50 "
              "already differ (REF_PEAK moved at a pass boundary)")
    return dict(tag=tag, i_in=i_in, i_on=i_on, conv=conv_rel,
                gain=float(gg[-1]), tp=tp)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lane", action="append", default=[],
                    help="DIR:ELEMENT:CHANNELS[:REF_PEAK]")
    ap.add_argument("--keep", default=None,
                    help="keep the wavs and traces in this directory")
    a = ap.parse_args()
    run_synthetic()
    rows = []
    if a.lane:
        tmp = a.keep or tempfile.mkdtemp(prefix="gate_alc_")
        os.makedirs(tmp, exist_ok=True)
        for spec in a.lane:
            rows.append(run_lane(spec, tmp, sys.executable))
        print("\n  lane                         in LUFS  out LUFS  converge  "
              "final gain  true peak")
        for r in rows:
            print(f"  {r['tag']:28s} {r['i_in']:7.1f}  {r['i_on']:8.1f}  "
                  f"{(r['conv'] if r['conv'] is not None else float('nan')):6.1f} s"
                  f"  {r['gain']:+8.1f} dB  {r['tp']:6.1f} dBTP")
    print(f"\nGATE {'PASS' if not FAILS else 'FAIL'}"
          + (f" -- {len(FAILS)} failed: {', '.join(FAILS)}" if FAILS else ""))
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
