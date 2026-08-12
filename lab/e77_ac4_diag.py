#!/usr/bin/env python3
"""E77 -- what is DIFFERENT about Fox's AC-4?

Our decoder produces real programme audio for RF33 every day (5_X English and
the `pair` Spanish simulcast, correlation referees 0.99+).  On Fox 45.100 it
manages 115/300 frames on `pair` and 0/300 on `5_X`.  So this is not a "yield
to nudge upward": something about the bitstream differs, and the TOC is where
a stream declares what it is.

Parses every frame's table of contents (m20, TS 103 190-2) for three streams
and prints the field census side by side:

    Fox 45.100      the stream under test
    RF33 pid13      MMTP programme audio, 5.1 English   (works)
    RF33 WHUT tsi20 ROUTE service audio                 (works)

Any field that is constant within a stream and different between them is a
candidate root cause.
"""
import argparse
import collections
import os
import struct
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

import m20_ac4_toc2 as T2   # noqa: E402


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


def frames_of(path, limit=None):
    """Slice an fMP4's mdat into samples using each traf's trun sizes."""
    b = open(path, "rb").read()
    out, moof = [], None
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
            if limit and len(out) >= limit:
                return out[:limit]
    return out[:limit] if limit else out


def concat_frames(paths, limit):
    out = []
    for p in paths:
        out += frames_of(p)
        if limit and len(out) >= limit:
            break
    return out[:limit] if limit else out


STREAMS = [
    ("FOX 45.100", ["e76_out/fox_audio1.mp4"]),
    ("FOX 45.100 #2", ["e76_out/fox_audio2.mp4"]),
    ("RF33 pid13 (5.1 eng)", ["../data/e30/live_audio_pid13.0004.m4s"]),
    ("RF33 WHUT tsi20", ["m13_out48/obj_239_255_32_1_8321_tsi20_toi225650590.bin",
                         "m13_out48/obj_239_255_32_1_8321_tsi20_toi225650591.bin"]),
    ("RF33 WHUT tsi30", ["m13_out48/obj_239_255_32_1_8321_tsi30_toi225650590.bin",
                         "m13_out48/obj_239_255_32_1_8321_tsi30_toi225650591.bin"]),
]

FIELDS = ("bitstream_version", "fs_index", "frame_rate_index",
          "b_iframe_global", "n_presentations", "total_n_substream_groups",
          "n_substreams", "payload_base", "channel_modes", "b_content_type",
          "n_channels_in_group", "b_channel_coded", "b_ajoc", "b_aspx",
          "b_asf", "sus_ver", "b_hsf_ext")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frames", type=int, default=60)
    a = ap.parse_args()
    for label, rels in STREAMS:
        paths = [os.path.join(HERE, r) for r in rels]
        paths = [p for p in paths if os.path.exists(p)]
        if not paths:
            print("%-24s MISSING" % label)
            continue
        fr = concat_frames(paths, a.frames)
        print("=" * 78)
        print("%-24s  %d frames  sizes min %d med %d max %d"
              % (label, len(fr), min(map(len, fr)),
                 sorted(map(len, fr))[len(fr) // 2], max(map(len, fr))))
        ok, bad = 0, 0
        cens = collections.defaultdict(collections.Counter)
        errs = collections.Counter()
        for f in fr:
            try:
                st = T2.parse(f)
                ok += 1
            except Exception as exc:                            # noqa: BLE001
                bad += 1
                errs["%s: %s" % (type(exc).__name__, str(exc)[:60])] += 1
                continue
            for k, v in st.items():
                if isinstance(v, list):
                    v = tuple(v)
                if isinstance(v, (int, str, tuple)):
                    cens[k][v] += 1
        print("  TOC parsed %d/%d" % (ok, len(fr)))
        for e, n in errs.most_common(3):
            print("    parse error x%-4d %s" % (n, e))
        for k in sorted(cens):
            c = cens[k]
            if k in ("sequence_counter",):
                continue
            vals = c.most_common(4)
            s = "  ".join("%s x%d" % (v, n) for v, n in vals)
            print("    %-26s %s" % (k, s[:110]))


if __name__ == "__main__":
    main()
