#!/usr/bin/env python3
"""E77 -- offline AC-4 decode through the PRODUCTION path, for any lane.

Uses m17_ac4_walk.samples() (the same extractor tools/atsc3_audio.py drives
live) rather than a private slicer, so a regression run here exercises the
code the live stack actually uses.

    python e77_ac4_decode.py <lane.mp4> --element 5_X|pair --channels 2|6 --out x.wav

Reports the honest signal statistics as well as the frame yield: peak, per
channel RMS and non-silent fraction.  A decoder that "runs" and emits silence,
or emits 1e12, has not decoded anything -- this campaign's own law.
"""
import argparse
import os
import sys
import wave

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

import m17_ac4_walk as M17             # noqa: E402
import m30_filterbank as FB            # noqa: E402
from m42_ac4_stream import Ac4Stream   # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("lane")
    ap.add_argument("--element", default="5_X", choices=("5_X", "pair"))
    ap.add_argument("--channels", type=int, default=2, choices=(2, 6))
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-frames", type=int, default=300)
    ap.add_argument("--skip", type=int, default=0)
    ap.add_argument("--lead", type=int, default=4,
                    help="frames of QMF lead-in discarded (A-SPX filter "
                         "memory starts cold; decoding them prints a click)")
    a = ap.parse_args()

    fr = M17.samples(a.lane)
    fr = [f for f in fr if len(f) > 16][a.skip:a.skip + a.max_frames]
    if not fr:
        print("no frames")
        return 1
    chans = Ac4Stream.CHANS if a.channels == 6 else ("L", "R")
    dec = Ac4Stream(element=a.element)
    wins, lens, grp, cnts = dec.decode_frames(fr, chans)
    yield_pc = 100.0 * dec.n_ok / max(dec.n_ok + dec.n_bad, 1)
    print("  frames %d   decoded ok %d  bad %d   YIELD %.1f%%"
          % (len(fr), dec.n_ok, dec.n_bad, yield_pc))
    if dec.n_ok == 0:
        print("  nothing decoded")
        return 1
    pcm, _ = FB.synthesise_frames(wins, lens, cnts, 1536,
                                  states={k: None for k in chans},
                                  return_state=True)
    n = min(len(v) for v in pcm.values())
    lead = a.lead * 2048
    ch = [np.asarray(pcm[k][lead:n], np.float32) for k in chans]
    if ch[0].size == 0:
        print("  no samples after lead-in")
        return 1
    peak = float(max(np.abs(c).max() for c in ch))
    scale = 0.9 / peak if peak > 0 else 1.0
    inter = np.empty(ch[0].size * len(ch), np.int16)
    for i, c in enumerate(ch):
        inter[i::len(ch)] = np.clip(c * scale * 32767, -32768,
                                    32767).astype(np.int16)
    with wave.open(a.out, "wb") as w:
        w.setnchannels(len(ch))
        w.setsampwidth(2)
        w.setframerate(48000)
        w.writeframes(inter.tobytes())
    rms = [float(np.sqrt((c ** 2).mean())) for c in ch]
    nonsilent = 100.0 * float(np.mean(np.abs(ch[0]) > peak * 1e-3))
    print("  %s  %.2f s  peak %.4g" % (a.out, ch[0].size / 48000.0, peak))
    print("  rms " + "  ".join("%s %.5g" % (k, v)
                               for k, v in zip(chans, rms)))
    print("  non-silent %.1f%%" % nonsilent)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
