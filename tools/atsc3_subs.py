#!/usr/bin/env python3
"""atsc3_subs.py -- turn the live caption lane into a growing .srt.

The cheapest of the three decoupled workers: the subtitle asset is one TTML
document per 2.002 s MPU, a few kB each, so this costs nothing and can run
as often as it likes.

Two things it inherits from m40, both learned the hard way:

  * ROLL-UP COLLAPSE. Live captioning re-emits the whole accumulated line
    every time a character arrives, so ~92 % of cues are prefixes of their
    successor. Collapse on the BOTTOM line -- the one being typed -- because
    a scroll makes the whole block SHORTER and "longest wins" swallows a
    completed line whole.
  * THE TIMES ARE ABSOLUTE wall clock, not media time, despite
    ttp:timeBase="media". Anchor them to the MPU grid.

Usage:
    python tools/atsc3_subs.py                  # follow the default dir
    python tools/atsc3_subs.py --once
"""
from __future__ import annotations

import argparse
import json
import os
import struct
import sys
import time
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

TTML_NS = "{http://www.w3.org/ns/ttml}"
MPU_SECONDS = 120 * 1001 / 60000.0


def lanes(d):
    try:
        return json.load(open(os.path.join(d, "live.json"))).get("lanes", {})
    except (OSError, ValueError):
        return {}


def boxes(b, off=0, end=None):
    end = len(b) if end is None else end
    p = off
    while p + 8 <= end:
        sz = struct.unpack(">I", b[p:p + 4])[0]
        t = b[p + 4:p + 8]
        if sz == 1:
            sz = struct.unpack(">Q", b[p + 8:p + 16])[0]
        elif sz == 0:
            sz = end - p
        if sz < 8 or p + sz > end:
            return
        yield t, p + 8, p + sz
        p += sz


def docs(path):
    """Every TTML document in the growing fMP4 lane, in order."""
    raw = open(path, "rb").read()
    out = []
    for t, s, e in boxes(raw):
        if t == b"mdat":
            body = raw[s:e]
            k = body.rfind(b"</tt>")
            if k > 0:
                out.append(body[:k + 5])
    return out


def text_of(p):
    parts = []

    def walk(el):
        if el.tag == TTML_NS + "br":
            parts.append("\n")
        if el.text:
            parts.append(el.text)
        for c in el:
            walk(c)
            if c.tail:
                parts.append(c.tail)

    if p.text:
        parts.append(p.text)
    for c in p:
        walk(c)
        if c.tail:
            parts.append(c.tail)
    return "".join(parts).strip()


def secs(v):
    v = v.strip()
    return float(v[:-1]) if v.endswith("s") else float(v)


def collapse(cues):
    """Roll-up runs -> one cue per completed bottom line."""
    cues = sorted(cues, key=lambda c: (c[0], c[1]))
    out, cur = [], None
    for b, e, t in cues:
        line = t.splitlines()[-1].strip() if t.strip() else ""
        if not line:
            continue
        if cur is not None and (line.startswith(cur[2]) or line == cur[2]):
            cur[1] = max(cur[1], e)
            cur[2] = line if len(line) > len(cur[2]) else cur[2]
        else:
            if cur:
                out.append(cur)
            cur = [b, e, line]
    if cur:
        out.append(cur)
    merged = []
    for g in out:
        if merged and merged[-1][2] == g[2]:
            merged[-1][1] = max(merged[-1][1], g[1])
        else:
            merged.append(g)
    return merged


def srt_time(t):
    t = max(0.0, t)
    h = int(t // 3600); t -= h * 3600
    m = int(t // 60); t -= m * 60
    s = int(t); ms = int(round((t - s) * 1000))
    if ms == 1000:
        s, ms = s + 1, 0
    return f"{h:02d}:{m:02d}:{s:02d},{ms:03d}"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--live-dir", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--once", action="store_true")
    ap.add_argument("--interval", type=float, default=5.0)
    ap.add_argument("--wait", type=float, default=90.0)
    a = ap.parse_args()
    d = a.live_dir or os.path.join(ROOT, "data", "atsc3_live")

    t0, lane = time.time(), None
    while time.time() - t0 < a.wait:
        cands = [l for l in lanes(d).values() if l.get("kind") == "subs"]
        if cands and os.path.exists(cands[0]["path"]):
            lane = cands[0]
            break
        time.sleep(2.0)
    if not lane:
        print(f"no caption lane in {d} -- run the chain with --assets all",
              file=sys.stderr)
        return 1

    out = a.out or os.path.join(d, "live.srt")
    print(f"caption lane pid {lane['pid']} -> {os.path.basename(lane['path'])}")
    print(f"  writing {out}")
    anchor = None
    try:
        while True:
            cues = []
            try:
                dd = docs(lane["path"])
            except OSError:
                dd = []
            for i, x in enumerate(dd):
                try:
                    root = ET.fromstring(x.decode("utf-8", "replace"))
                except ET.ParseError:
                    continue
                got = []
                for p in root.iter(TTML_NS + "p"):
                    b, e = p.get("begin"), p.get("end")
                    if not b or not e:
                        continue
                    t = text_of(p)
                    if t:
                        got.append((secs(b), secs(e), t))
                if got and anchor is None:
                    # doc i covers media [i*2.002, +2.002); recover the
                    # absolute->media offset from its earliest cue
                    anchor = min(b for b, _, _ in got) - i * MPU_SECONDS
                    # E35: record WHICH MPU SLOT srt time 0 belongs to, so
                    # the muxer can put cues on the lane clock EXACTLY
                    # instead of assuming this worker started at the lane's
                    # fragment 0. The anchor above treats doc index i as
                    # media time i*2.002, so srt t=0 sits i slots before
                    # doc i's OWN slot: anchor_seq = idx[i].seq - i (NOT
                    # first_seq -- a hole before doc i would shift that).
                    try:
                        idxp = (lane.get("idx")
                                or os.path.splitext(lane["path"])[0] + ".idx")
                        seqs = [json.loads(l)["seq"] for l in open(idxp)]
                        if i < len(seqs):
                            json.dump({"anchor_seq": seqs[i] - i,
                                       "srt": os.path.abspath(out)},
                                      open(os.path.join(
                                          d, "live_subs.json"), "w"))
                    except (OSError, ValueError):
                        pass
                cues.extend(got)
            if cues and anchor is not None:
                rel = [(b - anchor, e - anchor, t) for b, e, t in cues]
                merged = [c for c in collapse(rel) if c[1] > c[0] >= 0]
                tmp = out + ".tmp"
                with open(tmp, "w", encoding="utf-8") as fh:
                    for i, (b, e, t) in enumerate(merged, 1):
                        fh.write(f"{i}\n{srt_time(b)} --> {srt_time(e)}\n{t}\n\n")
                os.replace(tmp, out)
                print(f"  {len(dd)} documents, {len(cues)} raw cues -> "
                      f"{len(merged)} readable")
            if a.once:
                break
            time.sleep(a.interval)
    except KeyboardInterrupt:
        print("  stopped")
    return 0


if __name__ == "__main__":
    sys.exit(main())
