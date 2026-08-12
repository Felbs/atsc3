#!/usr/bin/env python3
"""M21 -- which AC-4 tools this broadcast actually uses, read off the air.

The question dolbyTuna has to answer before a line of decoder gets written is
not "how does AC-4 work" -- it is "which of AC-4 does THIS STREAM use".  The
standard is a menu; a broadcast picks a few dishes.  M20's parsed TOC put us at
the substream boundary, and the substream's own header answers it.

WHAT THE STREAM SAYS
--------------------
    channel_mode      5.1        every frame  (M20, and ffprobe agrees: 6 ch)
    n_substreams      2          every frame
    5_X_codec_mode    1 = ASPX   every frame

`5_X_codec_mode` (TS 103 190-1 Table 96) is the first three bits of the audio
payload and it enumerates exactly which of the optional tools are in play:

    0 SIMPLE        no A-SPX, no A-CPL
    1 ASPX          A-SPX only                   <-- THIS STREAM
    2 ASPX_ACPL_1   A-SPX + A-CPL mode 1
    3 ASPX_ACPL_2   A-SPX + A-CPL mode 2
    4 ASPX_ACPL_3   A-SPX + A-CPL mode 3

**A-CPL is not used.**  That is a real scope cut, not a hopeful one: coupling
is one of the two big optional tools, and this stream never switches it on in
116 seconds of audio.

THE GATE THAT MAKES THE SUBSTREAM SPLIT TRUSTWORTHY
----------------------------------------------------
`ac4_substream()` opens with a 15-bit `audio_size`, and everything after the
audio is `metadata()` plus byte alignment.  So `audio_size` must be no larger
than the substream, and should track it closely.  Measured over 3480 frames:

    audio_size <= substream size    3480/3480
    correlation                     +0.9998
    substream - audio_size          4..12 bytes, mean 6.7   (= metadata + fill)

A wrong substream boundary would not produce a 15-bit field that tracks the
size to four decimal places.

Usage:
    python m21_ac4_tools.py [m7_out/rf33_audio_pid13.mp4]
"""
from __future__ import annotations

import argparse
import collections
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m17_ac4_walk as W                                         # noqa: E402
import m20_ac4_toc2 as M                                         # noqa: E402
from m19_ac4_toc import Bits                                     # noqa: E402

CODEC_MODE = {0: "SIMPLE (no A-SPX, no A-CPL)", 1: "ASPX (A-SPX only)",
              2: "ASPX_ACPL_1", 3: "ASPX_ACPL_2", 4: "ASPX_ACPL_3"}


def substreams(frame, st):
    """-> [bytes] each substream's payload, using the TOC's own sizes."""
    out, off = [], st["toc_bytes"]
    for sz in st["substream_sizes"]:
        out.append(frame[off:off + sz])
        off += sz
    return out


def audio_header(sub):
    """ac4_substream()'s opening: -> (audio_size, bit position after it)."""
    b = Bits(sub)
    v = b.u(15)
    if b.u(1):
        v += b.vb(7) << 15
    return v, b


def main(argv=None):
    ap = argparse.ArgumentParser(prog="m21 ac4 tools")
    ap.add_argument("path", nargs="?", default="m7_out/rf33_audio_pid13.mp4")
    a = ap.parse_args(argv)
    p = a.path if os.path.isabs(a.path) else os.path.join(HERE, a.path)
    fr = W.samples(p)
    print("M21 -- which AC-4 tools this stream uses")
    print("=" * 72)

    rows, modes, ch = [], collections.Counter(), collections.Counter()
    fits = 0
    for f in fr:
        st = M.parse(f)
        subs = substreams(f, st)
        ch[tuple(st["channel_modes"])] += 1
        # the LAST substream is the audio one; the small first is presentation
        sz, b = audio_header(subs[-1])
        rows.append((sz, len(subs[-1])))
        fits += sz <= len(subs[-1])
        modes[b.u(3)] += 1

    arr = np.array(rows)
    d = arr[:, 1] - arr[:, 0]
    print(f"  {len(fr)} frames, {len(arr)} audio substreams")
    print(f"\n  GATE  audio_size <= substream size   {fits}/{len(arr)}")
    print(f"    correlation                        "
          f"{np.corrcoef(arr[:, 0], arr[:, 1])[0, 1]:+.4f}")
    print(f"    substream - audio_size             {d.min()}..{d.max()} "
          f"bytes, mean {d.mean():.1f}  (metadata + fill)")
    print(f"\n  channel_mode per frame   "
          f"{ {k: v for k, v in ch.items()} }")
    print(f"  5_X_codec_mode           " +
          str({CODEC_MODE.get(k, k): v for k, v in modes.items()}))

    only = list(modes)
    print("\n" + "=" * 72)
    if len(only) == 1 and only[0] == 1:
        print("  THE TOOL LIST, settled by the stream itself:")
        print("    REQUIRED  spectral front end (ASF/SSF) + quantisation")
        print("    REQUIRED  companding")
        print("    REQUIRED  A-SPX  (spectral extension)")
        print("    NOT USED  A-CPL  (coupling)   <- one whole tool skipped")
        print("    REQUIRED  MDCT synthesis, and it is where sound first"
              " appears")
        return 0
    print(f"  codec modes vary: {dict(modes)} -- scope is wider than one path")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
