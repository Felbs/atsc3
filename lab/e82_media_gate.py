#!/usr/bin/env python3
"""E82 -- gate the PIXELS and the SOUND that came out of the live LDM chain.

E67's law, amended by E76: `flat_frac` is a DISCRIMINATOR, not a constant, so
it is only meaningful when known-garbage runs in the same batch.  E67's three
CENC-encrypted RF33 services are that control -- same assembly, same decoder,
content that provably cannot decode without a key -- and E76 read the band off
the data as 0.75 (clear 0.0159-0.5252, garbage 0.9653-0.9913).

Audio is gated the E77 way: our own AC-4 decoder over the lane, reporting the
FRAME YIELD and the peak, because a decoder that "runs" and emits 2e12 has not
decoded anything.

    python e82_media_gate.py e82_out/live200
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import struct
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

FLAT_MAX = 0.75            # E76's amended band
TEMPORAL_MIN = 0.80

CONTROLS = [
    ("e67_out/whut_video.mp4", "RF33 32.1 WHUT clear   (POSITIVE control)"),
    ("e67_out/wttg_video.mp4", "RF33 5.1 WTTG CENC     (NEGATIVE control)"),
    ("e67_out/wrc_video.mp4", "RF33 4.1 WRC CENC      (NEGATIVE control)"),
    ("e67_out/wusa_video.mp4", "RF33 9.1 WUSA CENC     (NEGATIVE control)"),
]


def sample_entries(path):
    """stsd four-char codes -- `encv`/`enca` would mean encrypted."""
    buf = open(path, "rb").read(4 << 20)
    res = []

    def walk(s, e):
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
                walk(off + hdr, off + size)
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


def picture(path, nframes=48):
    p = subprocess.run(["ffprobe", "-v", "error", "-select_streams", "v:0",
                        "-show_entries", "stream=width,height,codec_name",
                        "-of", "json", path], capture_output=True, text=True)
    try:
        st = json.loads(p.stdout)["streams"][0]
        w, h = int(st["width"]), int(st["height"])
    except Exception:                                              # noqa: BLE001
        return dict(error="no video stream")
    q = subprocess.run(["ffmpeg", "-v", "error", "-i", path, "-frames:v",
                        str(nframes), "-pix_fmt", "gray", "-f", "rawvideo",
                        "-"], capture_output=True, timeout=600)
    n = len(q.stdout) // (w * h)
    if n < 2:
        return dict(error="too few frames", frames=n, w=w, h=h)
    a = np.frombuffer(q.stdout[:n * w * h], np.uint8).reshape(n, h, w)
    flats, rs = [], []
    for k in range(n):
        flats.append(np.bincount(a[k].ravel(), minlength=256).max() / a[k].size)
        if k:
            x = a[k - 1].ravel().astype(np.float32)
            y = a[k].ravel().astype(np.float32)
            sx, sy = x.std(), y.std()
            rs.append(float(((x - x.mean()) * (y - y.mean())).mean()
                            / (sx * sy)) if sx > 1e-6 and sy > 1e-6 else 0.0)
    gx = np.abs(np.diff(a[0].astype(np.int16), axis=1)).mean()
    gy = np.abs(np.diff(a[0].astype(np.int16), axis=0)).mean()
    return dict(codec=st.get("codec_name"), w=w, h=h, frames=n,
                flat_frac=round(float(np.mean(flats)), 4),
                temporal_r=round(float(np.mean(rs)), 4),
                gradient=round(float((gx + gy) / 2), 2),
                luma_std=round(float(a.std()), 2))


def audio(path, element="pair", channels=2, max_frames=200):
    """The E77 gate: our own AC-4 decoder, frame yield and honest levels."""
    import m17_ac4_walk as M17
    from m42_ac4_stream import Ac4Stream
    import m30_filterbank as FB
    try:
        fr = [f for f in M17.samples(path) if len(f) > 16][:max_frames]
    except Exception as e:                                         # noqa: BLE001
        return dict(error=f"{type(e).__name__}: {e}")
    if not fr:
        return dict(error="no AC-4 frames in the lane")
    chans = ("L", "R") if channels == 2 else Ac4Stream.CHANS
    dec = Ac4Stream(element=element)
    try:
        wins, lens, grp, cnts = dec.decode_frames(fr, chans)
    except Exception as e:                                         # noqa: BLE001
        return dict(error=f"decode_frames: {type(e).__name__}: {e}",
                    frames=len(fr))
    ok, bad = int(dec.n_ok), int(dec.n_bad)
    if ok == 0:
        return dict(frames=len(fr), frames_ok=0, frames_bad=bad,
                    yield_pc=0.0, entries=sample_entries(path),
                    note="nothing decoded with this element")
    try:
        # e77_ac4_decode's own call, verbatim -- same production stages the
        # live audio worker drives, including the A-SPX lead-in discard
        pcm, _st = FB.synthesise_frames(wins, lens, cnts, 1536,
                                        states={k: None for k in chans},
                                        return_state=True)
        n = min(len(v) for v in pcm.values())
        lead = 4 * 2048
        ch = [np.asarray(pcm[k][lead:n], np.float32) for k in chans]
        peak = float(max(np.abs(c).max() for c in ch))
        rms = [round(float(np.sqrt((c ** 2).mean())), 5) for c in ch]
        nz = round(float(np.mean(np.abs(ch[0]) > peak * 1e-3)), 4)
        secs = round(ch[0].size / 48000.0, 2)
    except Exception as e:                                         # noqa: BLE001
        return dict(error=f"synthesise: {type(e).__name__}: {e}",
                    frames=len(fr), frames_ok=ok, frames_bad=bad)
    return dict(frames=len(fr), frames_ok=ok, frames_bad=bad,
                yield_pc=round(100.0 * ok / max(ok + bad, 1), 1),
                seconds=secs, peak=round(peak, 5), rms=rms, nonsilent=nz,
                entries=sample_entries(path))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("live_dir")
    ap.add_argument("--json", default=os.path.join(HERE, "e82_out",
                                                   "media_gate.json"))
    ap.add_argument("--no-audio", action="store_true")
    a = ap.parse_args()
    meta = json.load(open(os.path.join(a.live_dir, "live.json")))
    rep = dict(lanes={}, controls={})
    print("E82 media gate -- the pixels and the sound off the LDM core layer")
    print("=" * 74)
    print("\n  LANES WRITTEN BY THE LIVE CHAIN")
    for name, ln in sorted(meta["lanes"].items()):
        print(f"    {name:28s} {ln['kind']:14s} pid {ln['pid']:5d}  "
              f"{ln['bytes'] / 1e6:8.3f} MB  {ln['segments']:4d} segments  "
              f"{ln['media_s']:8.3f} s  handler {ln['handler']}")
        rep["lanes"][name] = dict(ln)

    # E93's permanent invariant: a live receiver cannot bank media faster
    # than wall time. foxlive reported 2.6x wall for a month and nothing
    # checked; this one comparison would have caught it the first night,
    # and catches the NEXT encoder whose timescale isn't what we assumed.
    wall = float(meta.get("updated", 0) or 0) - float(meta.get("started", 0) or 0)
    if wall > 0 and meta.get("lanes"):
        worst = max(ln.get("media_s", 0.0) for ln in meta["lanes"].values())
        ok = worst <= wall * 1.02            # 2% grace for rounding
        verdict = "OK" if ok else \
            "VIOLATED - the accounting is lying (E93 read 2.667x here)"
        rep["controls"]["media_le_wall"] = dict(
            wall_s=round(wall, 1), worst_lane_media_s=round(worst, 1), ok=ok)
        print(f"\n  MEDIA<=WALL INVARIANT (E93): wall {wall:.0f} s, "
              f"worst lane {worst:.0f} s -> {verdict}")

    print("\n  PICTURE GATE (E67's law, E76's band: flat_frac < %.2f is a "
          "DISCRIMINATOR)" % FLAT_MAX)
    print("    %-42s %-11s %8s %8s %8s %8s" %
          ("file", "codec/size", "flat", "temporal", "gradient", "lumaStd"))
    rows = []
    for name, ln in sorted(meta["lanes"].items()):
        if ln["handler"] != "vide":
            continue
        pg = picture(ln["path"])
        ents = sample_entries(ln["path"])
        pg["sample_entries"] = ents
        pg["encrypted"] = any(e in ("encv", "enca") for e in ents)
        rep["lanes"][name]["picture"] = pg
        rows.append((name, pg))
    for path, label in CONTROLS:
        p = os.path.join(HERE, path)
        if not os.path.exists(p):
            continue
        pg = picture(p)
        rep["controls"][label] = pg
        rows.append((label, pg))
    for name, pg in rows:
        if "error" in pg:
            print("    %-42s %s" % (name[:42], pg["error"]))
            continue
        print("    %-42s %-11s %8.4f %8.4f %8.2f %8.2f"
              % (name[:42], f"{pg['codec']} {pg['w']}x{pg['h']}",
                 pg["flat_frac"], pg["temporal_r"], pg["gradient"],
                 pg["luma_std"]))

    ours = [(n, pg) for n, pg in rows if n in rep["lanes"] and "error" not in pg]
    garbage = [pg for lbl, pg in rep["controls"].items()
               if "NEGATIVE" in lbl and "error" not in pg]
    ok = bool(ours) and all(pg["flat_frac"] < FLAT_MAX
                            and pg["temporal_r"] > TEMPORAL_MIN
                            and not pg.get("encrypted") for _n, pg in ours)
    sep = None
    if ours and garbage:
        sep = min(pg["flat_frac"] for pg in garbage) - max(
            pg["flat_frac"] for _n, pg in ours)
    print(f"\n    VERDICT  {'REAL PICTURE' if ok else 'GATE FAILED'}"
          + (f"   separation from known garbage: {sep:.4f}"
             if sep is not None else
             "   (no garbage control present -- the band is NOT calibrated "
             "in this run)"))
    rep["picture_pass"] = ok
    rep["separation"] = sep

    if not a.no_audio:
        print("\n  AUDIO GATE (E77: our own AC-4 decoder, frame yield + level)")
        for name, ln in sorted(meta["lanes"].items()):
            if ln["handler"] != "soun":
                continue
            for elem in ("pair", "5_X"):
                r = audio(ln["path"], element=elem)
                rep["lanes"][name].setdefault("audio", {})[elem] = r
                if "error" in r:
                    print(f"    {name:28s} {elem:5s} {r['error'][:64]}")
                elif not r.get("frames_ok"):
                    print(f"    {name:28s} {elem:5s} 0/{r['frames']} frames "
                          f"-- {r.get('note', '')}")
                else:
                    print(f"    {name:28s} {elem:5s} "
                          f"{r['frames_ok']}/{r['frames_ok'] + r['frames_bad']}"
                          f" frames ({r['yield_pc']:.1f}%)  {r['seconds']:.2f}s"
                          f"  peak {r['peak']:.4f}  rms {r['rms']}  "
                          f"non-silent {r['nonsilent']:.3f}  "
                          f"entries {','.join(r['entries'])}")

    json.dump(rep, open(a.json, "w"), indent=1, default=str)
    print(f"\n  wrote {a.json}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
