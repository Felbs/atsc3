#!/usr/bin/env python3
"""M35 -- the A-SPX HF generator: the patch copy, and nothing after it.

This is deliberately ONE stage, not the whole high-band chain.  A-SPX
reconstructs 12..21 kHz in four steps:

    1. HF GENERATION   copy low QMF subbands up into the A-SPX range   <- here
    2. envelope adjustment   scale them to the transmitted envelopes
    3. noise addition        fill where the envelope says noise, not tone
    4. additional harmonics  insert sinusoids where flagged

Only step 1 is implemented.  Steps 2-4 are what make the result *correct*; the
patched band on its own has the right content at the WRONG LEVEL, typically far
too loud, because the whole point of the envelope stage is that a copy of
7.5-12 kHz is not what 12-21 kHz should sound like.

So this file **does not render audio**.  It gates the patch in the QMF domain,
where the claim is exact and checkable, and stops there.  Rendering a
half-finished high band would produce a file that plays and is wrong, which is
the failure mode this project has guarded against throughout.

WHAT THE PATCH IS (M31, Pseudocode 71)
----------------------------------------
    patch 0: 12 subbands from 20..32 (7500..12000 Hz) -> 12000..16500 Hz
    patch 1: 12 subbands from 20..32 (7500..12000 Hz) -> 16500..21000 Hz

Both source the same block -- the top 4.5 kHz of the MDCT-coded band, ending
exactly at the crossover.

THE GATES
----------
  1. after patching, each destination range holds EXACTLY the energy of its
     source range -- a copy is a copy, so this is equality, not similarity;
  2. the A-SPX range is empty BEFORE the patch (M33's integration gate, which
     is what makes step 1's "before" well defined);
  3. subbands above the A-SPX stop (56 = 21000 Hz) stay empty -- A-SPX does not
     reach there, so a correct patch must not put energy there either.

Usage:
    python m35_hfgen.py
"""
from __future__ import annotations

import argparse
import os
import sys
import wave

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m28_channels as C                                          # noqa: E402
import m31_aspx_bands as A                                        # noqa: E402
import m33_qmf as Q                                               # noqa: E402


def hf_generate(Qm, t):
    """Step 1 only: copy the source blocks into the A-SPX range.

    Q_high[sbx + k] = Q_low[patch_start_sb[p] + k] for each patch p.
    """
    out = Qm.copy()
    for p, (n, st) in enumerate(zip(t["patch_num_sb"], t["patch_start_sb"])):
        dst = t["patches"][p]
        out[dst:dst + n, :] = Qm[st:st + n, :]
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(prog="m35 hfgen")
    ap.add_argument("--audio", default="tv_audio_valid.wav")
    ap.add_argument("--slots", type=int, default=600)
    a = ap.parse_args(argv)
    wav = a.audio if os.path.isabs(a.audio) else os.path.join(HERE, a.audio)

    print("M35 -- A-SPX HF generator (the patch copy only)")
    print("=" * 74)
    if not os.path.exists(wav):
        print(f"  no {os.path.basename(wav)}")
        return 1

    # the A-SPX geometry this stream signals (M31)
    cfg = dict(start_freq=5, stop_freq=1, master_freq_scale=1, noise_sbg=3)
    t = A.sbg_tables(cfg, 0)
    hz = C.QMF_HZ
    print(f"\n  A-SPX range subbands {t['sbx']}..{t['sbx'] + t['num_sb_aspx']}"
          f" = {t['sbx'] * hz:.0f}..{(t['sbx'] + t['num_sb_aspx']) * hz:.0f} Hz")
    for p, (n, st) in enumerate(zip(t["patch_num_sb"], t["patch_start_sb"])):
        d = t["patches"][p]
        print(f"    patch {p}: sb {st}..{st + n} -> sb {d}..{d + n}"
              f"   ({st * hz:.0f}..{(st + n) * hz:.0f} Hz"
              f"  ->  {d * hz:.0f}..{(d + n) * hz:.0f} Hz)")

    w = Q.qwin()
    wv = wave.open(wav)
    n = min(wv.getnframes(), 64 * a.slots)
    d = np.frombuffer(wv.readframes(n), "<i2").reshape(n, wv.getnchannels())
    x = d[:, 0].astype(float) / 32768.0
    x = x[:(len(x) // Q.NSB) * Q.NSB]
    Qm = Q.analyse(x, w, Q.analysis_matrix())
    E0 = (np.abs(Qm) ** 2).sum(axis=1)

    print(f"\n  1. before the patch")
    aspx0 = E0[t["sbx"]:t["sbx"] + t["num_sb_aspx"]].sum()
    g2 = aspx0 / E0.sum() < 1e-5
    print(f"    {'PASS' if g2 else 'FAIL'}  A-SPX range holds "
          f"{100 * aspx0 / E0.sum():.5f} % of the energy (it should be empty)")

    Qh = hf_generate(Qm, t)
    E1 = (np.abs(Qh) ** 2).sum(axis=1)

    print(f"\n  2. after the patch -- a copy is a copy, so this is EQUALITY")
    g1 = True
    for p, (nsb, st) in enumerate(zip(t["patch_num_sb"],
                                      t["patch_start_sb"])):
        dst = t["patches"][p]
        src_e = E0[st:st + nsb].sum()
        dst_e = E1[dst:dst + nsb].sum()
        rel = abs(dst_e - src_e) / max(src_e, 1e-30)
        ok = rel < 1e-12
        g1 &= ok
        print(f"    {'PASS' if ok else 'FAIL'}  patch {p}: source energy "
              f"{src_e:.6e}, destination {dst_e:.6e}, rel diff {rel:.1e}")

    print(f"\n  3. nothing above the A-SPX stop")
    top = t["sbx"] + t["num_sb_aspx"]
    above = E1[top:].sum() / E1.sum()
    g3 = above < 1e-6
    print(f"    {'PASS' if g3 else 'FAIL'}  subbands {top}..63 "
          f"({top * hz:.0f}+ Hz) hold {100 * above:.6f} % -- A-SPX does not "
          f"reach there")

    print(f"\n  the patched high band now carries "
          f"{100 * E1[t['sbx']:top].sum() / E1.sum():.2f} % of the energy, "
          f"at the WRONG level")
    print(f"  until envelope adjustment scales it -- which is why no audio is "
          f"written here.")

    print("\n" + "=" * 74)
    if g1 and g2 and g3:
        print("  HF GENERATOR CORRECT.  The A-SPX range is filled by an exact "
              "copy of the\n  signalled source blocks.  Steps 2-4 (envelope, "
              "noise, harmonics) remain.")
        return 0
    print("  NOT established")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
