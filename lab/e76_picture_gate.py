#!/usr/bin/env python3
"""E76 -- the picture gate, RECALIBRATED against known-garbage controls.

E67 set the acceptance band at flat_frac < 0.10 on 1920x1080 broadcast news.
Fox 45.100 is a 960x540 mobile feed and its first frame is a McDonald's advert
on a WHITE background: flat_frac 0.525, which the E67 threshold would have
called garbage.  It is not garbage -- it is a perfect picture, confirmed by
eye.

So the threshold was overfitted to one kind of content.  The fix is not to
move the number by feel, it is to put the KNOWN-GARBAGE cases in the same run
and read the separation.  E67's three CENC-encrypted RF33 services are exactly
that: same extractor, same assembly code, same decoder, content that provably
cannot decode without a key.

Run this and the band reads itself.
"""
import json
import os
import subprocess

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

FILES = [
    ("e76_out/fox_video.mp4", "RF25 45.100 Fox mobile (E76, under test)"),
    ("e67_out/whut_video.mp4", "RF33 32.1 WHUT, clear    (E67 positive ctl)"),
    ("e67_out/wttg_video.mp4", "RF33 5.1 WTTG, CENC      (E67 NEGATIVE ctl)"),
    ("e67_out/wrc_video.mp4", "RF33 4.1 WRC, CENC       (E67 NEGATIVE ctl)"),
    ("e67_out/wusa_video.mp4", "RF33 9.1 WUSA, CENC      (E67 NEGATIVE ctl)"),
]


def stats(path, nframes=24):
    rc = subprocess.run(["ffprobe", "-v", "error", "-select_streams", "v:0",
                         "-show_entries", "stream=width,height",
                         "-of", "json", path], capture_output=True, text=True)
    try:
        st = json.loads(rc.stdout)["streams"][0]
        w, h = int(st["width"]), int(st["height"])
    except Exception:                                           # noqa: BLE001
        return None
    p = subprocess.run(["ffmpeg", "-v", "error", "-i", path, "-frames:v",
                        str(nframes), "-pix_fmt", "gray", "-f", "rawvideo",
                        "-"], capture_output=True, timeout=300)
    n = len(p.stdout) // (w * h)
    if n < 2:
        return dict(w=w, h=h, frames=n, flat=None)
    a = np.frombuffer(p.stdout[:n * w * h], np.uint8).reshape(n, h, w)
    flats = [np.bincount(a[k].ravel(), minlength=256).max() / a[k].size
             for k in range(n)]
    # edge energy: real pictures have structure at many scales; concealment
    # fills large areas with one value and has almost no gradient
    gx = np.abs(np.diff(a[0].astype(np.int16), axis=1)).mean()
    gy = np.abs(np.diff(a[0].astype(np.int16), axis=0)).mean()
    return dict(w=w, h=h, frames=n,
                flat=float(np.mean(flats)),
                grad=float((gx + gy) / 2),
                luma_std=float(a.std()))


def main():
    rows = []
    print("%-46s %-11s %-7s %-8s %-9s %s"
          % ("file", "size", "frames", "flat_frac", "grad", "luma_std"))
    print("-" * 96)
    for rel, label in FILES:
        p = os.path.join(HERE, rel)
        if not os.path.exists(p):
            print("%-46s MISSING" % label)
            continue
        s = stats(p)
        if not s or s["flat"] is None:
            print("%-46s no frames" % label)
            continue
        rows.append((label, s))
        print("%-46s %-11s %-7d %-8.4f %-9.2f %.2f"
              % (label, "%dx%d" % (s["w"], s["h"]), s["frames"], s["flat"],
                 s["grad"], s["luma_std"]))
    clear = [s["flat"] for l, s in rows if "NEGATIVE" not in l]
    garb = [s["flat"] for l, s in rows if "NEGATIVE" in l]
    print()
    if clear and garb:
        print("  clear content   flat_frac %.4f .. %.4f" % (min(clear), max(clear)))
        print("  known garbage   flat_frac %.4f .. %.4f" % (min(garb), max(garb)))
        gap = min(garb) - max(clear)
        mid = (min(garb) + max(clear)) / 2
        print("  SEPARATION      %.4f  -> recalibrated threshold %.2f "
              "(was 0.10, overfitted to 1080p news)" % (gap, mid))
        gclear = [s["grad"] for l, s in rows if "NEGATIVE" not in l]
        ggarb = [s["grad"] for l, s in rows if "NEGATIVE" in l]
        print("  gradient  clear %.2f .. %.2f   garbage %.2f .. %.2f"
              % (min(gclear), max(gclear), min(ggarb), max(ggarb)))
    json.dump([{"label": l, **s} for l, s in rows],
              open(os.path.join(HERE, "e76_picture_gate.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
