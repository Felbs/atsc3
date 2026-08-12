#!/usr/bin/env python3
"""gate_e48.py -- audio frame-grid carry: every muxed chunk must contain
EXACTLY as much audio as timeline it spans. The old fixed-288288 slicing
(250.25 mp2 frames) is the NEGATIVE CONTROL: ffmpeg pads the final partial
frame, every chunk carries +18 ms, and the surplus-ppm assertion fails --
proving the gate can fail and that the old code does.

    python lab/gate_e48.py    -> PASS/FAIL, exit 0 iff all pass
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tools"))

import atsc3_tv as tv                                           # noqa: E402

CHUNK = 3 * tv.SPF                # 288288 samples = 6.006 s = 250.25 frames
FAILS = []


def check(name, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" +
          (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


def gate_carry_arithmetic():
    print("gate 1: frame_align is lossless and frame-exact over 12 chunks")
    ramp = np.arange(12 * CHUNK, dtype="<i2") % 20000
    src = np.stack([ramp, ramp], 1)
    carry = {"eng": None, "spa": None}
    outs, leads = [], []
    for k in range(12):
        pcm = src[k * CHUNK:(k + 1) * CHUNK]
        e, s, lead = tv.frame_align(pcm, pcm.copy(), carry)
        outs.append(e)
        leads.append(lead)
    lens = [len(o) for o in outs]
    check("every chunk length is a whole number of mp2 frames",
          all(n % tv.MP2_SPF == 0 for n in lens), f"lens {lens[:6]}...")
    check("250/250/250/251 pattern, 1001 frames per 4 chunks",
          lens[:4] == [288000, 288000, 288000, 289152] and
          sum(lens[:4]) == 4 * CHUNK, f"first 4: {lens[:4]}")
    cat = np.vstack(outs)
    n = len(cat)
    check("concatenated output == input sample-exact (no zeros, no dups)",
          bool((cat[:, 0] == src[:n, 0]).all()) and n == 12 * CHUNK - len(carry["eng"]),
          f"emitted {n}, deferred {len(carry['eng'])}")
    check("itsoffset lead always equals the deferred tail's length",
          all(leads[k + 1] == (k + 1) * CHUNK - sum(lens[:k + 1])
              for k in range(11)), f"leads {leads[:6]}...")


def mp2_frames(wav_path, tmp):
    """Encode a wav to mp2 exactly as the mux does; -> packet count."""
    out = os.path.join(tmp, "probe.ts")
    subprocess.run([tv.ffbin("ffmpeg"), "-y", "-loglevel", "error",
                    "-i", wav_path, "-c:a", "mp2", "-b:a", "224k",
                    "-muxdelay", "0", "-muxpreload", "0",
                    "-f", "mpegts", out], check=True, capture_output=True)
    r = subprocess.run([tv.ffbin("ffprobe"), "-loglevel", "error",
                        "-select_streams", "a:0", "-count_packets",
                        "-show_entries", "stream=nb_read_packets",
                        "-of", "csv=p=0", out],
                       capture_output=True, text=True)
    return int(r.stdout.split()[0].rstrip(","))


def gate_encoder_padding():
    print("gate 2: real ffmpeg -- aligned feed never pads; the old feed does")
    d = tempfile.mkdtemp(prefix="e48_")
    try:
        pcm = (np.random.default_rng(7).integers(-2000, 2000,
               (CHUNK, 2))).astype("<i2")
        old = os.path.join(d, "old.wav")     # 288288 = the shipped slicing
        new = os.path.join(d, "new.wav")     # 288000 = frame-aligned
        tv.write_wav(old, pcm)
        tv.write_wav(new, pcm[:288000])
        n_old = mp2_frames(old, d)
        n_new = mp2_frames(new, d)
        surplus_old = n_old * tv.MP2_SPF - CHUNK
        check("aligned chunk encodes to exactly 250 frames, zero surplus",
              n_new == 250 and n_new * tv.MP2_SPF == 288000,
              f"{n_new} frames")
        # NEGATIVE CONTROL -- the shipped 288288-sample slicing on the
        # same audio: the encoder pads the 251st frame, +864 samples
        # (+18 ms) of surplus per chunk = the measured +2988 ppm drift.
        check("NEGATIVE CONTROL: old slicing yields 251 frames, +864 "
              "samples surplus", n_old == 251 and surplus_old == 864,
              f"{n_old} frames, surplus {surplus_old}")
    finally:
        shutil.rmtree(d, ignore_errors=True)


if __name__ == "__main__":
    gate_carry_arithmetic()
    gate_encoder_padding()
    print(f"\n{'ALL GATES PASS' if not FAILS else 'FAILED: ' + ', '.join(FAILS)}")
    sys.exit(1 if FAILS else 0)
