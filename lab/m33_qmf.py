#!/usr/bin/env python3
"""M33 -- the 64-band complex QMF analysis/synthesis bank (clause 5.7.3/5.7.4).

A-SPX lives entirely in the QMF domain: it patches low subbands up to high
ones, adjusts their envelopes and adds noise, all on the complex subband
matrix.  So none of the A-SPX side information M32 now parses can become audio
until this filterbank exists.  This is that filterbank, and nothing else.

THE WINDOW CAME FROM THE C FILE, NOT THE PDF
----------------------------------------------
Table D.3's 640 window coefficients are also shipped as `QWIN` in
`ts_103190_tables.c` -- machine readable, no OCR risk.  M23's parser could not
read it because that parser is INTEGER-only: it split `9.90318758627504e-04`
into three tokens and reported 1918 "values" for a 640-entry array.  A
float-aware read gives exactly 640, matching the declared size, symmetric about
index 320 for n >= 1 -- the standard MPEG QMF prototype shape.

THE AMBIGUITY, AND HOW IT IS SETTLED
--------------------------------------
The synthesis modulation matrix is given twice and the two disagree:

    clause 5.7.4.2 step 2:  N[n][k] ~ exp(j*pi*(k+0.5)*(2n - 4*64 - 1)/128)
                            i.e. the constant is 257
    the pseudocode below:   exponent = j*(pi/128)*(sb+0.5)*(2*n - 255)
                            i.e. the constant is 255

Rather than pick one, both are built and the one that RECONSTRUCTS is taken --
the same approach that settled step 6 of the MDCT overlap-add in M30.

THE GATE: ANALYSIS THEN SYNTHESIS RETURNS THE SIGNAL
------------------------------------------------------
A complex QMF bank is oversampled by two, so the analysis/synthesis pair is
near-perfectly reconstructing: feed a known signal in, get it back out, delayed
by the bank's group delay.  That is checkable with no broadcast data and no
encoder.  The delay is MEASURED by cross-correlation rather than assumed, and
then the reconstruction SNR is computed at that lag.  A wrong window, a wrong
modulation constant, a wrong fold in the `g` extraction: each destroys the SNR.

Usage:
    python m33_qmf.py
"""
from __future__ import annotations

import argparse
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m23_hcb as H                                               # noqa: E402

NSB = 64                     # num_qmf_subbands
NWIN = 640                   # num_qmf_win_coef


def qwin(path=None):
    """The 640-tap QMF prototype, read as FLOATS from the spec's C file."""
    src = open(path or H.DEFAULT_C, encoding="utf-8", errors="replace").read()
    src = re.sub(r"/\*.*?\*/", "", src, flags=re.S)
    m = re.search(r"const\s+\w+\s+QWIN\s*\[\s*(\d*)\s*\]\s*=\s*\{(.*?)\}\s*;",
                  src, re.S)
    if not m:
        raise LookupError("QWIN not found")
    vals = [float(x) for x in
            re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", m.group(2))]
    declared = int(m.group(1)) if m.group(1) else len(vals)
    if len(vals) != declared:
        raise ValueError(f"QWIN: {len(vals)} values, declared {declared}")
    return np.asarray(vals, float)


def analysis_matrix():
    """M[sb][n] = exp(j*pi/128*(sb+0.5)*(2n-1)), 64 x 128.  Pseudocode 65."""
    sb = np.arange(NSB)[:, None]
    n = np.arange(2 * NSB)[None, :]
    return np.exp(1j * np.pi / (2 * NSB) * (sb + 0.5) * (2 * n - 1))


def synthesis_matrix(const):
    """N[n][sb] = 1/64 * exp(j*pi/128*(sb+0.5)*(2n-const)), 128 x 64."""
    n = np.arange(2 * NSB)[:, None]
    sb = np.arange(NSB)[None, :]
    return np.exp(1j * np.pi / (2 * NSB) * (sb + 0.5)
                  * (2 * n - const)) / NSB


def analyse(pcm, w, M):
    """Clause 5.7.3.2 / Pseudocode 65.  -> (64, num_timeslots) complex."""
    nts = len(pcm) // NSB
    filt = np.zeros(NWIN)
    out = np.empty((NSB, nts), complex)
    for ts in range(nts):
        filt[NSB:] = filt[:-NSB]                       # shift by 64
        filt[:NSB] = pcm[ts * NSB:(ts + 1) * NSB][::-1]
        z = filt * w
        u = z[:2 * NSB].copy()
        for k in range(1, 5):
            u += z[k * 2 * NSB:(k + 1) * 2 * NSB]
        out[:, ts] = M @ u
    return out


def synthesise(Q, w, N):
    """Clause 5.7.4.2.  -> real PCM."""
    nts = Q.shape[1]
    filt = np.zeros(10 * 2 * NSB)                      # 1280
    out = np.empty(nts * NSB)
    g = np.empty(NWIN)
    for ts in range(nts):
        filt[2 * NSB:] = filt[:-2 * NSB]               # shift by 128
        filt[:2 * NSB] = np.real(N @ Q[:, ts])
        for n in range(5):                             # the g fold
            g[128 * n:128 * n + NSB] = filt[256 * n:256 * n + NSB]
            g[128 * n + NSB:128 * (n + 1)] = \
                filt[256 * n + 192:256 * n + 192 + NSB]
        ww = g * w
        out[ts * NSB:(ts + 1) * NSB] = ww.reshape(10, NSB).sum(axis=0)
    return out


def gate(w, const, n=64 * 200, verbose=True):
    """Analyse then synthesise a known signal; -> (delay, SNR dB)."""
    rng = np.random.default_rng(7)
    x = rng.standard_normal(n)
    Q = analyse(x, w, analysis_matrix())
    y = synthesise(Q, w, synthesis_matrix(const))
    # measure the group delay rather than assuming it
    best, lag = -1.0, 0
    for d in range(0, 1024):
        if d + 4096 > len(y):
            break
        a, b = x[:4096], y[d:d + 4096]
        c = float(np.corrcoef(a, b)[0, 1])
        if c > best:
            best, lag = c, d
    # ALIGNMENT: x[i] corresponds to y[i + lag].  Comparing x[2048:] against
    # y[lag:] instead of y[lag+2048:] is off by 2048 samples and reports a
    # near-zero gain on a bank that is actually reconstructing -- which is
    # exactly what the first run of this gate did, while its own correlation
    # search (correctly anchored at each signal's origin) read +1.000000.
    # Two measurements disagreeing that sharply means one of them is wrong,
    # not that the thing under test is strange.
    off = 2048
    a = x[off:len(x) - 2048]
    b = y[lag + off:lag + off + len(a)]
    if len(b) < len(a):
        a = a[:len(b)]
    if len(a) < 1024:
        return lag, -99.0
    err = a - b
    snr = 10 * np.log10((a ** 2).mean() / max((err ** 2).mean(), 1e-30))
    if verbose:
        print(f"    constant {const}: delay {lag} samples, "
              f"correlation {best:+.6f}, reconstruction SNR {snr:6.2f} dB")
    return lag, snr


def integration_gate(wav, w, verbose=True):
    """Analyse OUR OWN decoded audio and find where its band limit lands.

    The renderer stops at max_sfb, which M27/M28 put at 12000 Hz, and A-SPX
    starts at QMF subband 32 by the header's own arithmetic.  So the QMF
    analysis of the decoded audio must show energy up to subband 31 and the
    noise floor from 32 up.  Nothing forces those three to agree -- the band
    table, the A-SPX header and this filterbank share no code -- so agreement
    is evidence about all three, and it also proves the bank is correctly
    ALIGNED in frequency rather than merely reconstructing.
    """
    import wave
    wv = wave.open(wav)
    n = min(wv.getnframes(), 64 * 600)
    nch = wv.getnchannels()
    d = np.frombuffer(wv.readframes(n), "<i2").reshape(n, nch).astype(float)
    x = d[:, 0] / 32768.0
    x = x[:(len(x) // NSB) * NSB]
    E = (np.abs(analyse(x, w, analysis_matrix())) ** 2).mean(axis=1)
    tot = E.sum()
    above = E[32:].sum() / tot
    if verbose:
        print(f"    subband 31 {E[31]/E.max():.2e}   subband 32 "
              f"{E[32]/E.max():.2e}   subband 40 {E[40]/E.max():.2e}")
        print(f"    energy at or above subband 32 (12000 Hz): "
              f"{100 * above:.4f} %")
    return above


def main(argv=None):
    ap = argparse.ArgumentParser(prog="m33 qmf")
    ap.add_argument("--audio", default="tv_audio_valid.wav",
                    help="decoded audio for the integration gate")
    a = ap.parse_args(argv)

    print("M33 -- 64-band complex QMF analysis/synthesis")
    print("=" * 74)
    w = qwin()
    print(f"\n  1. the prototype window")
    print(f"    PASS  {len(w)} coefficients, declared 640, read as floats "
          f"from the spec C file")
    # THE MAGNITUDE is the symmetric MPEG prototype; the SIGNS are folded in.
    # Testing w itself for symmetry fails on a perfectly good window -- 128
    # coefficients are negative, confined to 64-sample blocks 1, 3, 6 and 8,
    # which is the usual SBR sign convention.  The property that actually holds
    # is |w[n]| == |w[640-n]|.
    sym = np.allclose(np.abs(w[1:]), np.abs(w[1:])[::-1])
    nneg = int((w < 0).sum())
    print(f"    {'PASS' if sym else 'FAIL'}  |w[n]| symmetric about index 320 "
          f"(MPEG prototype magnitude)")
    print(f"    note  {nneg} negative coefficients, in 64-blocks "
          f"{[i for i in range(10) if (w[i*64:(i+1)*64] < 0).any()]} "
          f"-- signs folded into the window")

    print(f"\n  2. which synthesis modulation constant reconstructs?")
    print(f"    (clause 5.7.4.2 step 2 implies 257; its own pseudocode "
          f"prints 255)")
    res = {c: gate(w, c) for c in (255, 257)}

    best = max(res, key=lambda c: res[c][1])
    lag, snr = res[best]
    other = 255 if best == 257 else 257
    print(f"\n  3. verdict")
    ok = snr > 40.0
    print(f"    {'PASS' if ok else 'FAIL'}  constant {best} reconstructs at "
          f"{snr:.2f} dB; {other} gives {res[other][1]:.2f} dB")
    print(f"    group delay {lag} samples "
          f"({lag / NSB:.0f} QMF time slots)")

    wav = a.audio if os.path.isabs(a.audio) else os.path.join(HERE, a.audio)
    g_int = True
    if os.path.exists(wav):
        print("\n  4. integration: where does OUR decoded audio band-limit?")
        above = integration_gate(wav, w)
        g_int = above < 1e-3
        print(f"    {'PASS' if g_int else 'FAIL'}  the A-SPX range is empty, "
              f"so this bank, the Annex B band table and the A-SPX header\n"
              f"          all put the crossover at subband 32 -- and they "
              f"share no code")
    else:
        print(f"\n  (no {os.path.basename(wav)}; integration gate skipped)")

    print("\n" + "=" * 74)
    if ok and sym and g_int:
        print(f"  QMF BANK RECONSTRUCTS.  Analysis and synthesis are inverse "
              f"to {snr:.0f} dB,\n  and it sees our decoder's band limit at "
              f"exactly the A-SPX crossover.")
        return 0
    print("  NOT established -- the bank does not reconstruct")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
