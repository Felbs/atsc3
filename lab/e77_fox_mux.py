#!/usr/bin/env python3
"""E77 -- build fox_with_audio.mp4 and GATE THE AUDIO.

"The decoder ran" is not "the decoder decoded" (E67), and E76 showed the same
trap on the audio side: a run that produced 9.43 s of samples with a peak of
2e12 had decoded nothing.  So this does not just mux -- it measures whether
the waveform behaves like programme audio:

    peak / RMS          sane absolute scale, not 1e12 and not silence
    crest factor        speech and music sit around 3-15; noise ~4 but flat
    non-silent fraction a real programme is mostly not silence
    spectral centroid   broadcast audio is bottom-heavy; white noise is not
    HF/LF ratio         "
    L vs R              identical channels mean a mono source (say so)

The audio is our own AC-4 decode; it is transcoded PCM -> AAC only because
MP4 has no tag for AC-4 (ffmpeg: "Could not find tag for codec ac4").  The
underlying decode is ours, the transcode is for the container.
"""
import argparse
import json
import os
import subprocess
import sys
import wave

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def run(cmd, timeout=600):
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    return p.returncode, p.stdout, p.stderr


def audio_gate(wav):
    w = wave.open(wav)
    n, ch, sr = w.getnframes(), w.getnchannels(), w.getframerate()
    a = np.frombuffer(w.readframes(n), np.int16).reshape(-1, ch).astype(np.float64)
    L = a[:, 0]
    R = a[:, 1] if ch > 1 else L
    rms = float(np.sqrt((L ** 2).mean()))
    peak = float(np.abs(L).max())
    crest = peak / rms if rms > 0 else float("inf")
    nonsilent = float(np.mean(np.abs(L) > peak * 1e-3))
    # spectrum over the middle of the file
    seg = L[len(L) // 4: len(L) // 4 + 1 << 18] if len(L) > (1 << 18) else L
    seg = seg[: (len(seg) // 1024) * 1024].reshape(-1, 1024) * np.hanning(1024)
    S = np.abs(np.fft.rfft(seg, axis=1)).mean(axis=0)
    f = np.fft.rfftfreq(1024, 1.0 / sr)
    cent = float((S * f).sum() / max(S.sum(), 1e-9))
    lo = float(S[f < 1000].sum())
    hi = float(S[f > 6000].sum())
    return dict(seconds=round(n / sr, 2), channels=ch, rate=sr,
                rms=round(rms, 1), peak=round(peak, 1),
                crest=round(crest, 2), nonsilent=round(nonsilent, 4),
                centroid_hz=round(cent, 1),
                hf_lf_ratio=round(hi / max(lo, 1e-9), 4),
                lr_identical=bool(np.array_equal(L, R)),
                rms_L_minus_R=round(float(np.sqrt(((L - R) ** 2).mean())), 2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--video", default=os.path.join(HERE, "e76_out",
                                                    "fox_video.mp4"))
    ap.add_argument("--wav", default=os.path.join(HERE, "e77_fox_eng_full.wav"))
    ap.add_argument("--out", default=os.path.join(HERE, "e77_audio_out",
                                                  "fox_with_audio.mp4"))
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out), exist_ok=True)

    g = audio_gate(a.wav)
    print("AUDIO GATE on %s" % os.path.basename(a.wav))
    for k, v in g.items():
        print("   %-16s %s" % (k, v))
    sane_scale = 1e2 < g["rms"] < 3e4
    ok = (sane_scale and g["nonsilent"] > 0.5 and 2.0 < g["crest"] < 30.0
          and 100 < g["centroid_hz"] < 6000 and g["hf_lf_ratio"] < 0.5)
    print("   VERDICT          %s"
          % ("REAL PROGRAMME AUDIO" if ok else "NOT CREDIBLE AS AUDIO"))

    rc, so, se = run(["ffmpeg", "-y", "-v", "error",
                      "-i", a.video, "-i", a.wav,
                      "-map", "0:v:0", "-map", "1:a:0",
                      "-c:v", "copy", "-c:a", "aac", "-b:a", "160k",
                      "-shortest", a.out])
    print("\nmux rc=%d %s" % (rc, se.strip()[:200]))
    if rc == 0:
        rc2, so2, _ = run(["ffprobe", "-v", "error", "-show_entries",
                           "stream=index,codec_type,codec_name,channels,"
                           "width,height", "-of", "json", a.out])
        print(so2.strip())
        print("  size %d bytes" % os.path.getsize(a.out))
    json.dump(dict(gate=g, mux_rc=rc, out=a.out),
              open(os.path.join(HERE, "e77_fox_mux.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
