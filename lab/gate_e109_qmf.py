#!/usr/bin/env python3
"""gate_e109_qmf.py -- the compiled AC-4 QMF bank's referee.

The kernel replaces m33_qmf.analyse and .synthesise, the two Python loops at
the heart of A-SPX high-frequency regeneration.  Their output is audio, so
the bar is stated in audio terms as well as in ulp.

  leg 1  AGREEMENT with the NumPy reference on random signals and on a real
         decoded PCM channel, both directions of the bank.
  leg 2  The bank's OWN round-trip gate, re-run through the kernel: analyse
         then synthesise a known signal and measure the reconstruction SNR
         and group delay, at both constants the reference tests (255, 257).
         This checks the DSP, not merely that two implementations agree.
  leg 3  REAL AUDIO: run apply_hf_pair over real AC-4 frames with the kernel
         on and off and compare the rendered PCM -- the samples a listener
         would actually hear.

Run: python lab/gate_e109_qmf.py [--lane <live_audio_pidN.m4s>]
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m33_qmf as Q                                             # noqa: E402

_COMPARED = [0]


def _use(kernel: bool):
    os.environ["ATSC3_QMF_KERNEL"] = "1" if kernel else "0"
    Q._QMF_STATE["lib"] = None
    return Q._qmf_lib() if kernel else None


def _cmp(tag, a, b, fails, tol=1e-10):
    _COMPARED[0] += 1
    d = float(np.abs(np.asarray(a) - np.asarray(b)).max())
    scale = float(np.abs(np.asarray(b)).max()) or 1.0
    ok = d <= tol * max(1.0, scale)
    print(f"    {tag}: max|d| {d:.3e}  (rel {d/scale:.3e})  "
          f"{'OK' if ok else 'TOO LARGE'}")
    if not ok:
        fails.append(f"{tag} differs by {d:.3e}")


def leg1(fails):
    print("LEG 1 -- agreement with the NumPy reference")
    if _use(True) is None:
        fails.append("kernel not loadable")
        return
    w = Q.qwin()
    M = Q.analysis_matrix()
    N = Q.synthesis_matrix(255)
    rng = np.random.default_rng(20260823)
    cases = {"white noise": rng.standard_normal(64 * 400),
             "silence": np.zeros(64 * 50),
             "impulse": np.eye(1, 64 * 50, 37).ravel(),
             "loud tone": 3.0 * np.sin(np.arange(64 * 400) * 0.31)}
    for tag, x in cases.items():
        _use(True)
        Ak = Q.analyse(x, w, M)
        Sk = Q.synthesise(Ak, w, N)
        _use(False)
        An = Q.analyse(x, w, M)
        Sn = Q.synthesise(An, w, N)
        print(f"  {tag}")
        _cmp("analyse  ", Ak, An, fails)
        _cmp("synthesise", Sk, Sn, fails)


def leg2(fails):
    """The bank's own round-trip gate, driven through the kernel."""
    print("LEG 2 -- round-trip reconstruction (the DSP, not just agreement)")
    if _use(True) is None:
        fails.append("kernel not loadable for leg 2")
        return
    w = Q.qwin()
    for const in (255, 257):
        lag, snr = Q.gate(w, const, verbose=False)
        _COMPARED[0] += 1
        ok = snr >= 40.0
        print(f"  const {const}: SNR {snr:6.1f} dB at group delay {lag}  "
              f"{'OK' if ok else 'TOO LOW'}")
        if not ok:
            fails.append(f"round-trip SNR {snr:.1f} dB at const {const}")


def leg3(fails, lane):
    print(f"LEG 3 -- real audio through apply_hf_pair: {lane}")
    if not lane or not os.path.exists(lane):
        print("  no AC-4 lane given -- SKIPPED (legs 1-2 still bind, but run "
              "this before shipping: it is the only leg a listener hears)")
        return
    import m17_ac4_walk as W
    import m30_filterbank as FB
    import m37_render_hf as R37
    from m42_ac4_stream import Ac4Stream
    fr = W.samples(lane)[:80]
    out = {}
    for kernel in (True, False):
        _use(kernel)
        dec = Ac4Stream(element="5_X")
        wins, lens, grp, cnts = dec.decode_frames(fr, Ac4Stream.CHANS)
        pcm, _ = FB.synthesise_frames(wins, lens, cnts, 1536, states=None,
                                      return_state=True)
        gl = grp["lr"]
        if not any(g is not None for g in gl):
            print("  this lane carries no A-SPX groups -- SKIPPED")
            return
        hf, _n = R37.apply_hf_pair(pcm["L"], pcm["R"], gl, dec.cfg)
        out[kernel] = np.asarray(hf["L"]).copy()
    _cmp("rendered PCM", out[True], out[False], fails, tol=1e-9)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--lane", default="")
    a = ap.parse_args(argv)
    fails = []
    leg1(fails)
    leg2(fails)
    leg3(fails, a.lane)
    if _COMPARED[0] == 0:
        fails.append("the gate performed ZERO comparisons")
    print(f"\ncomparisons performed: {_COMPARED[0]}")
    if fails:
        print(f"GATE FAILED ({len(fails)}):")
        for f in fails:
            print(f"  - {f}")
        return 1
    print("GATE PASSED -- the QMF kernel reproduces the reference bank")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
