#!/usr/bin/env python3
"""E76 -- decode Fox 45.100's AC-4 audio with OUR decoder, offline.

tools/atsc3_audio.py is a live worker: it follows a growing lane registry and
carries resume state.  Here the whole track is already on disk as an fMP4, so
this drives the same two stages directly --

    m42_ac4_stream.Ac4Stream.decode_frames   AC-4 frames -> MDCT windows
    m30_filterbank.synthesise_frames         windows      -> PCM

-- after slicing the fMP4's mdat into AC-4 frames using each traf's trun
sample sizes (an AC-4 sample IS one AC-4 frame, TS 103 190-2 Annex E).

Stereo only (channel_pair_element): the ROUTE service's lanes are stereo, per
E47.  The first LEAD frames are dropped because the A-SPX QMF bank has filter
memory and starts cold -- decoding them would print a click, and that click is
exactly the sort of artefact a "the decoder ran" gate would pass.
"""
import argparse
import os
import struct
import sys
import wave

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

from m42_ac4_stream import Ac4Stream   # noqa: E402
import m30_filterbank as FB            # noqa: E402


def boxes(b, s, e):
    o = s
    while o + 8 <= e:
        n = struct.unpack(">I", b[o:o + 4])[0]
        t = b[o + 4:o + 8]
        h = 8
        if n == 1:
            n = struct.unpack(">Q", b[o + 8:o + 16])[0]
            h = 16
        elif n == 0:
            n = e - o
        if n < h or o + n > e:
            return
        yield t, o, o + h, o + n
        o += n


def frames_of(path):
    """Slice every fragment's mdat into AC-4 frames via trun sample sizes."""
    b = open(path, "rb").read()
    out = []
    moof = None
    for t, st, ps, pe in boxes(b, 0, len(b)):
        if t == b"moof":
            moof = (st, ps, pe)
        elif t == b"mdat" and moof is not None:
            sizes, doff = [], None
            for t2, st2, ps2, pe2 in boxes(b, moof[1], moof[2]):
                if t2 != b"traf":
                    continue
                for t3, st3, ps3, pe3 in boxes(b, ps2, pe2):
                    if t3 != b"trun":
                        continue
                    fl = struct.unpack(">I", b[ps3:ps3 + 4])[0] & 0xFFFFFF
                    cnt = struct.unpack(">I", b[ps3 + 4:ps3 + 8])[0]
                    o = ps3 + 8
                    if fl & 0x1:
                        doff = struct.unpack(">i", b[o:o + 4])[0]
                        o += 4
                    if fl & 0x4:
                        o += 4
                    for _ in range(cnt):
                        if fl & 0x100:
                            o += 4
                        if fl & 0x200:
                            sizes.append(struct.unpack(">I", b[o:o + 4])[0])
                            o += 4
                        if fl & 0x400:
                            o += 4
                        if fl & 0x800:
                            o += 4
            pos = (moof[0] + doff) if doff is not None else ps
            for sz in sizes:
                if pos + sz <= pe:
                    out.append(b[pos:pos + sz])
                pos += sz
            moof = None
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mp4", default=os.path.join(HERE, "e76_out",
                                                  "fox_audio1.mp4"))
    ap.add_argument("--out", default=os.path.join(HERE, "e76_out",
                                                  "fox_audio.wav"))
    ap.add_argument("--element", default="pair")
    ap.add_argument("--max-frames", type=int, default=400)
    ap.add_argument("--lead", type=int, default=4)
    a = ap.parse_args()

    fr = frames_of(a.mp4)
    print("AC-4 frames sliced from %s: %d" % (os.path.basename(a.mp4), len(fr)))
    if not fr:
        return 1
    print("  frame sizes: min %d  median %d  max %d"
          % (min(map(len, fr)), int(np.median([len(x) for x in fr])),
             max(map(len, fr))))
    fr = fr[:a.max_frames]
    dec = Ac4Stream(element=a.element)
    chans = ("L", "R")
    wins, lens, grp, cnts = dec.decode_frames(fr, chans)
    print("  decoded: %d ok, %d bad" % (dec.n_ok, dec.n_bad))
    pcm, _st = FB.synthesise_frames(wins, lens, cnts, 1536,
                                    states={k: None for k in chans},
                                    return_state=True)
    n = min(len(v) for v in pcm.values())
    lead = a.lead * 2048
    L = np.asarray(pcm["L"][lead:n], np.float32)
    R = np.asarray(pcm["R"][lead:n], np.float32)
    if L.size == 0:
        print("  no samples")
        return 1
    peak = float(max(np.abs(L).max(), np.abs(R).max()))
    scale = 0.9 / peak if peak > 0 else 1.0
    inter = np.empty(L.size * 2, np.int16)
    inter[0::2] = np.clip(L * scale * 32767, -32768, 32767).astype(np.int16)
    inter[1::2] = np.clip(R * scale * 32767, -32768, 32767).astype(np.int16)
    with wave.open(a.out, "wb") as w:
        w.setnchannels(2)
        w.setsampwidth(2)
        w.setframerate(48000)
        w.writeframes(inter.tobytes())
    dur = L.size / 48000.0
    # honest signal statistics -- silence would prove nothing decoded
    rms_l = float(np.sqrt((L ** 2).mean()))
    rms_r = float(np.sqrt((R ** 2).mean()))
    print("  wrote %s  %.2f s  peak %.4f  rms L %.5f R %.5f"
          % (a.out, dur, peak, rms_l, rms_r))
    print("  non-silent samples: %.1f%%"
          % (100.0 * np.mean(np.abs(L) > peak * 1e-3)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
