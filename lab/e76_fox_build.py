#!/usr/bin/env python3
"""E76 -- assemble Fox 45.100 (WBFFMob) into a playable file, then GATE THE PIXELS.

E67's law, earned the hard way: "the decoder ran" is not "the decoder decoded".
ffprobe parses CENC-encrypted video happily, ffmpeg emits frames and logs zero
errors, and NAL length prefixes tile perfectly -- all on content that renders
as a green field.  So the acceptance test here is flat_frac (the fraction of
pixels in the single most common luma value) plus frame-to-frame correlation,
not "it muxed".

Also reports the sample entries, so an encrypted track would be caught before
anything is claimed.
"""
import json
import os
import re
import struct
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "e76_fox_out")
OUT = os.path.join(HERE, "e76_out")
os.makedirs(OUT, exist_ok=True)
FLOW = "239_255_45_100_5000"
W, H = None, None


def segs(tsi):
    pat = re.compile(r"obj_%s_tsi%d_toi(\d+)\.bin$" % (FLOW, tsi))
    out = []
    for f in sorted(os.listdir(SRC)):
        m = pat.match(f)
        if not m:
            continue
        toi = int(m.group(1))
        if toi == 4294967295:
            continue
        out.append((toi, os.path.join(SRC, f)))
    return [p for _, p in sorted(out)]


def init(tsi):
    p = os.path.join(SRC, "obj_%s_tsi%d_toi4294967295.bin" % (FLOW, tsi))
    return p if os.path.exists(p) else None


def sample_entries(buf):
    """Walk to stsd and return the four-char codes (encv/enca = encrypted)."""
    res = []

    def walk(s, e, depth=0):
        off = s
        while off + 8 <= e:
            size = struct.unpack(">I", buf[off:off + 4])[0]
            typ = buf[off + 4:off + 8].decode("latin1", "replace")
            hdr = 8
            if size == 1:
                size = struct.unpack(">Q", buf[off + 8:off + 16])[0]
                hdr = 16
            elif size == 0:
                size = e - off
            if size < hdr or off + size > e:
                return
            if typ in ("moov", "trak", "mdia", "minf", "stbl"):
                walk(off + hdr, off + size, depth + 1)
            elif typ == "stsd":
                o = off + hdr + 8
                while o + 8 <= off + size:
                    s2 = struct.unpack(">I", buf[o:o + 4])[0]
                    res.append(buf[o + 4:o + 8].decode("latin1", "replace"))
                    if s2 < 8:
                        break
                    o += s2
            off += size
    walk(0, len(buf))
    return res


def build(tsi, name, nseg=24, skip=1):
    """skip=1: the first object of each flow is a partial -- we joined the
    ROUTE session mid-object, so its leading LCT packets are missing and it is
    zero-filled (it starts 00 00 00 00 where every complete segment starts
    with a `sidx` box).  Discarded identically for every track."""
    i = init(tsi)
    if not i:
        return None
    blob = bytearray(open(i, "rb").read())
    ents = sample_entries(bytes(blob))
    ss = segs(tsi)[skip:skip + nseg]
    for p in ss:
        blob += open(p, "rb").read()
    out = os.path.join(OUT, name)
    open(out, "wb").write(bytes(blob))
    return dict(path=out, bytes=len(blob), segments=len(ss),
                sample_entries=ents,
                encrypted=any(e in ("encv", "enca") for e in ents))


def run(cmd, timeout=300):
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    return p.returncode, p.stdout, p.stderr


def picture_gate(mp4, nframes=24):
    rc, so, se = run(["ffprobe", "-v", "error", "-select_streams", "v:0",
                      "-show_entries", "stream=width,height,codec_name",
                      "-of", "json", mp4])
    try:
        st = json.loads(so)["streams"][0]
        w, h = int(st["width"]), int(st["height"])
    except Exception:                                           # noqa: BLE001
        return {"error": "no video stream", "ffprobe": se[:200]}
    p = subprocess.run(["ffmpeg", "-v", "error", "-i", mp4,
                        "-frames:v", str(nframes), "-pix_fmt", "gray",
                        "-f", "rawvideo", "-"], capture_output=True,
                       timeout=300)
    raw = p.stdout
    n = len(raw) // (w * h)
    if n < 2:
        return {"error": "too few frames", "frames": n, "w": w, "h": h}
    a = np.frombuffer(raw[:n * w * h], np.uint8).reshape(n, h, w)
    flats, rs = [], []
    for k in range(n):
        c = np.bincount(a[k].ravel(), minlength=256)
        flats.append(c.max() / a[k].size)
        if k:
            x = a[k - 1].ravel().astype(np.float32)
            y = a[k].ravel().astype(np.float32)
            sx, sy = x.std(), y.std()
            rs.append(float(((x - x.mean()) * (y - y.mean())).mean() / (sx * sy))
                      if sx > 1e-6 and sy > 1e-6 else 0.0)
    return dict(codec=st.get("codec_name"), w=w, h=h, frames=n,
                flat_frac_mean=round(float(np.mean(flats)), 4),
                flat_frac_max=round(float(np.max(flats)), 4),
                temporal_r_mean=round(float(np.mean(rs)), 4),
                luma_mean=round(float(a.mean()), 2),
                luma_std=round(float(a.std()), 2))


def main():
    rep = {}
    v = build(10, "fox_video.mp4")
    a1 = build(20, "fox_audio1.mp4")
    a2 = build(30, "fox_audio2.mp4")
    rep["video"] = v
    rep["audio1"] = a1
    rep["audio2"] = a2
    for k in ("video", "audio1", "audio2"):
        r = rep[k]
        if r:
            print("%-7s %-28s %8d bytes  %2d segs  entries %-14s %s"
                  % (k, os.path.basename(r["path"]), r["bytes"], r["segments"],
                     ",".join(r["sample_entries"]),
                     "ENCRYPTED" if r["encrypted"] else "clear"))
    if v:
        pg = picture_gate(v["path"])
        rep["picture_gate"] = pg
        print()
        print("PICTURE GATE (E67 law -- gate the pixels, not the parse):")
        for k2, val in pg.items():
            print("   %-18s %s" % (k2, val))
        ok = (pg.get("flat_frac_mean", 1) < 0.10
              and pg.get("temporal_r_mean", 0) > 0.80)
        print("   VERDICT           %s"
              % ("REAL VIDEO" if ok else "GARBAGE / NO PICTURE"))
        run(["ffmpeg", "-y", "-v", "error", "-i", v["path"], "-frames:v", "1",
             os.path.join(OUT, "fox_frame.png")])
        run(["ffmpeg", "-y", "-v", "error", "-i", v["path"], "-ss", "2",
             "-frames:v", "1", os.path.join(OUT, "fox_frame2.png")])
    # mux video + first audio into fox.mp4
    if v and a1:
        rc, so, se = run(["ffmpeg", "-y", "-v", "error", "-i", v["path"],
                          "-i", a1["path"], "-c", "copy",
                          os.path.join(OUT, "fox.mp4")])
        rep["mux_rc"] = rc
        rep["mux_err"] = se[:300]
        print("\nmux fox.mp4 rc=%d %s" % (rc, se[:200]))
    json.dump(rep, open(os.path.join(HERE, "e76_fox_build.json"), "w"),
              indent=1, default=str)


if __name__ == "__main__":
    main()
