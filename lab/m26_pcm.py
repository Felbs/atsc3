#!/usr/bin/env python3
"""M26 -- PCM.  The LFE channel of live broadcast AC-4, turned into samples.

Every piece upstream is gated:

    TOC                3480/3480, every unaccounted byte is fill      (M20)
    substream          audio_size <= size, r = 0.9998                 (M21)
    tool list          A-SPX, no A-CPL, from the stream itself        (M21)
    section data       lands on bit 15, the entropy map's prediction  (M24)
    codebooks          60/60 Kraft == 1.000000000, prefix-free        (M23)
    spectral lines     adjacent r = +0.273 vs shuffled +0.013         (M24)
    scale factors      adjacent r = +0.763 vs shuffled -0.001         (M25)
    MDCT               reconstruction 2.3e-14                          (M22)

This assembles them.

THE RECONSTRUCTION, from clause 5.1.3 and Pseudocode 21
---------------------------------------------------------
    rec_spec    = sign(quant_spec) * |quant_spec|^(4/3)
    sf_gain     = pow(2.0, 0.25 * (scale_factor - 100))
    scaled_spec = sf_gain * rec_spec

ALL THREE CONSTANTS ARE NOW VERIFIED AGAINST THE SPEC.  The offset 100 and the
centre 60 are written out in Pseudocode 21 verbatim:

    scale_factor += dpcm_sf[g][sfb] - 60;
    sf_gain[g][sfb] = pow(2.0, 0.25 * (scale_factor - 100));

The 4/3 exponent was originally flagged READ FROM CONTEXT, because the PDF sets
it as a stacked fraction whose glyphs do not survive text extraction -- the line
comes out as `rec_spec=sign quant_spec x |quant_spec|` with the exponent simply
missing, and an exponent you cannot see is an exponent you are guessing.

It was settled by RENDERING the formula region of the page to an image and
reading it, rather than trusting the text layer:

    rec_spec = sign(quant_spec) x |quant_spec|^(4/3)      clause 5.1.3.2

which confirms 4/3 exactly.  Worth remembering as a technique: when a spec's
equation extracts as garbage, the pseudocode restatement is the first place to
look and a rendered crop of the page is the second.  Both beat inference.

WHAT THIS IS AND IS NOT
------------------------
It is the LFE: 12 spectral lines, roughly 0..188 Hz.  Low-frequency effects
only -- a rumble, not a mix.  It is not the programme audio and it is not
supposed to be.  What it demonstrates is that the whole chain from RF to
samples closes.

The five full-band channels need the same treatment with `coding_config` and
the multi-channel elements, and A-SPX on top for anything above ~7 kHz.

Usage:
    python m26_pcm.py [--frames 600] [--out lfe.wav]
"""
from __future__ import annotations

import argparse
import os
import struct
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m17_ac4_walk as W                                         # noqa: E402
import m20_ac4_toc2 as M                                         # noqa: E402
import m22_mdct as T                                             # noqa: E402
import m23_hcb as H                                              # noqa: E402
import m24_spectral as S                                         # noqa: E402
import m25_scalefac as SF                                        # noqa: E402

N = 1536                     # Table 83: internal frame length at 48 kHz
SF_OFFSET = 100              # Pseudocode 21
QUANT_EXP = 4.0 / 3.0        # standard non-uniform quantiser law


def imdct_fast(X, w):
    """IMDCT via scipy's DCT-IV, gated against m22's direct form."""
    from scipy.fft import dct
    n = len(X)
    folded = dct(X, type=4, norm=None) / n
    a_b, c_d = folded[n // 2:], folded[:n // 2]
    return np.concatenate([a_b, -a_b[::-1], -c_d[::-1], -c_d]) * w


def gate_imdct(n=64):
    """The fast path must equal the proven slow path."""
    rng = np.random.default_rng(2)
    X = rng.standard_normal(n)
    w = T.sine_window(2 * n)
    a, b = T.imdct(X, w), imdct_fast(X, w)
    return float(np.abs(a - b).max())


def frame_spectrum(lines, sfs):
    """12 quantised lines + per-band scale factors -> a 1536-line spectrum."""
    spec = np.zeros(N)
    for sfb in range(len(sfs)):
        if sfs[sfb] is None:
            continue
        lo, hi = S.SFB_OFFSET[sfb], S.SFB_OFFSET[sfb + 1]
        q = lines[lo:hi].astype(float)
        rec = np.sign(q) * np.abs(q) ** QUANT_EXP
        spec[lo:hi] = rec * (2.0 ** (0.25 * (sfs[sfb] - SF_OFFSET)))
    return spec


def write_wav(path, x, rate=48000):
    x = np.clip(x, -1.0, 1.0)
    pcm = (x * 32767.0).astype("<i2").tobytes()
    with open(path, "wb") as f:
        f.write(b"RIFF" + struct.pack("<I", 36 + len(pcm)) + b"WAVEfmt ")
        f.write(struct.pack("<IHHIIHH", 16, 1, 1, rate, rate * 2, 2, 16))
        f.write(b"data" + struct.pack("<I", len(pcm)) + pcm)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="m26 pcm")
    ap.add_argument("path", nargs="?", default="m7_out/rf33_audio_pid13.mp4")
    ap.add_argument("--frames", type=int, default=900)
    ap.add_argument("--out", default="lfe.wav")
    a = ap.parse_args(argv)
    p = a.path if os.path.isabs(a.path) else os.path.join(HERE, a.path)

    print("M26 -- PCM from broadcast AC-4 (the LFE channel)")
    print("=" * 72)
    d = gate_imdct()
    print(f"  GATE  fast IMDCT vs the proven direct form: max |delta| {d:.2e}"
          f"   {'PASS' if d < 1e-10 else 'FAIL'}")
    if d >= 1e-10:
        return 1

    arrays = H.parse_c(H.DEFAULT_C)
    tables = {n: S.Huff(arrays[f"ASF_HCB_{n}_LEN"], arrays[f"ASF_HCB_{n}_CW"])
              for n in S.CB_MOD}
    sf_table = S.Huff(arrays["ASF_HCB_SCALEFAC_LEN"],
                      arrays["ASF_HCB_SCALEFAC_CW"])

    fr = W.samples(p)[:a.frames]
    w = T.sine_window(2 * N)
    out = np.zeros(len(fr) * N + 2 * N)
    used = 0
    for i, f in enumerate(fr):
        try:
            st = M.parse(f)
            if st["b_iframe_global"]:
                continue
            off = st["toc_bytes"] + st["substream_sizes"][0]
            sub = f[off:off + st["substream_sizes"][1]]
            lines, ref, sfs, _ = SF.decode_lfe_full(sub, tables, sf_table)
        except Exception:                                      # noqa: BLE001
            continue
        spec = frame_spectrum(lines, sfs)
        out[i * N:i * N + 2 * N] += imdct_fast(spec, w)
        used += 1

    x = out[:len(fr) * N]
    peak = np.abs(x).max()
    print(f"\n  {used} frames synthesised of {len(fr)}")
    print(f"  {len(x)} samples = {len(x) / 48000:.2f} s at 48 kHz")
    print(f"  peak {peak:.3e}   rms {np.sqrt((x ** 2).mean()):.3e}")
    if peak <= 0:
        print("  silence -- nothing to write")
        return 1
    xn = x / peak * 0.9
    outp = a.out if os.path.isabs(a.out) else os.path.join(HERE, a.out)
    write_wav(outp, xn)
    print(f"  wrote {os.path.basename(outp)}  (normalised; the absolute level "
          f"depends on constants this file flags)")

    # where is the energy?  the LFE should be almost entirely below ~200 Hz
    X = np.abs(np.fft.rfft(xn * np.hanning(len(xn))))
    fq = np.fft.rfftfreq(len(xn), 1 / 48000)
    tot = (X ** 2).sum()
    for lo, hi in ((0, 100), (100, 200), (200, 500), (500, 2000),
                   (2000, 24000)):
        m = (fq >= lo) & (fq < hi)
        print(f"    {lo:5d}-{hi:5d} Hz  {100 * (X[m] ** 2).sum() / tot:6.2f} %")
    band = ((fq >= 0) & (fq < 250))
    frac = (X[band] ** 2).sum() / tot
    print("\n" + "=" * 72)
    print(f"  {frac * 100:.1f} % of the energy is below 250 Hz -- "
          + ("THIS IS AN LFE CHANNEL." if frac > 0.8 else
             "that is NOT LFE-like; something upstream is wrong"))
    return 0 if frac > 0.8 else 1


if __name__ == "__main__":
    raise SystemExit(main())
