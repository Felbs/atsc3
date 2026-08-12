#!/usr/bin/env python3
"""atsc3_meter.py -- quality meters for a live dir: proof without a viewer.

This is how we KNOW picture/sound/captions work without a human watching:

  VIDEO  decode the last N seconds of the video lane into ffmpeg's null
         sink and report frames, errors, and fps (the fleet's standard:
         real decode quality = null-sink fps, not file growth).
  AUDIO  per decoded track (eng / spa wavs): fraction of the last N
         seconds' SLOTS carrying non-silent audio (RMS gate), plus a
         clipping fraction. Walks the slot grid via the audio lane's own
         index -- the same mapping the muxer uses, so it measures what the
         viewer would actually hear.
  CC     cue count over the last N minutes of live.srt, and whether cue
         LANE-times (via the live_subs.json anchor, E35) land inside the
         lane's slot range -- anchor sanity, the exact failure v1 had.

Output: one JSON line on stdout (last line) + a human summary before it.

Usage:
    python tools/atsc3_meter.py --live-dir data/live51
    python tools/atsc3_meter.py --live-dir data/e29 --video-secs 30
"""
from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
import wave

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

MPU_SECONDS = 2.002
SPF = 96096
FPS = 60000 / 1001.0
RMS_SILENCE = 100.0            # int16 units, ~ -50 dBFS
CLIP_LEVEL = 32600


def ffbin(name):
    p = shutil.which(name)
    if p:
        return p
    c = os.path.join(
        os.environ.get("LOCALAPPDATA", "C:"), "Microsoft", "WinGet",
        "Packages", "Gyan.FFmpeg_Microsoft.Winget.Source_8wekyb3d8bbwe",
        "ffmpeg-8.1-full_build", "bin", name + ".exe")
    return c if os.path.exists(c) else name


def read_json(p):
    try:
        return json.load(open(p))
    except (OSError, ValueError):
        return None


def read_idx(p):
    out = []
    try:
        for line in open(p):
            r = json.loads(line)
            out.append((r["seq"], r["off"], r["len"]))
    except (OSError, ValueError):
        pass
    return out


def resolve(live_dir, p):
    if not p:
        return None
    for cand in (p, os.path.join(ROOT, p),
                 os.path.join(live_dir, os.path.basename(p))):
        if os.path.exists(cand):
            return cand
    return None


def lane_of(lanes, kind, pid=None):
    for l in lanes.values():
        if l.get("kind") == kind and (pid is None or l.get("pid") == pid):
            return l
    return None


def meter_video(live_dir, lanes, secs):
    v = lane_of(lanes, "video")
    if not v:
        return {"error": "no video lane"}
    idx = read_idx(os.path.splitext(resolve(live_dir, v["path"]))[0] + ".idx")
    if not idx:
        return {"error": "empty video idx"}
    n = max(1, math.ceil(secs / MPU_SECONDS))
    rows = idx[-n:]
    span_slots = rows[-1][0] - rows[0][0] + 1
    tmp = os.path.join(live_dir, "_tv")
    os.makedirs(tmp, exist_ok=True)
    m4s = os.path.join(tmp, "_meter.m4s")
    src_p = resolve(live_dir, v["path"])
    with open(src_p, "rb") as src, open(m4s, "wb") as dst:
        src.seek(0)
        dst.write(src.read(v["init_bytes"] or 0))
        for _, off, ln in rows:
            src.seek(off)
            dst.write(src.read(ln))
    t0 = time.time()
    r = subprocess.run([ffbin("ffmpeg"), "-hide_banner", "-i", m4s,
                        "-f", "null", "-"], capture_output=True, text=True,
                       errors="replace")
    wall = time.time() - t0
    err = r.stderr or ""
    m = re.findall(r"frame=\s*(\d+)", err)
    frames = int(m[-1]) if m else 0
    err_lines = [l for l in err.splitlines()
                 if re.search(r"error|corrupt|invalid|missing", l, re.I)
                 and "frame=" not in l]
    expected = span_slots * MPU_SECONDS * FPS
    try:
        os.remove(m4s)
    except OSError:
        pass
    return {
        "secs": round(span_slots * MPU_SECONDS, 3),
        "slots": span_slots, "frags_present": len(rows),
        "frames": frames, "frames_expected": round(expected),
        "decode_complete": round(frames / expected, 4) if expected else 0.0,
        "null_sink_fps": round(frames / wall, 1) if wall > 0 else 0.0,
        "decode_errors": len(err_lines),
        "error_sample": err_lines[:3],
    }


def meter_audio_track(live_dir, lanes, sidecar, secs):
    import numpy as np
    am = read_json(os.path.join(live_dir, sidecar))
    if not am:
        return {"present": False, "why": "no sidecar"}
    wav_p = resolve(live_dir, am.get("path"))
    if not wav_p:
        return {"present": False, "why": "wav missing"}
    lane = lane_of(lanes, "audio", am.get("pid"))
    if not lane:
        return {"present": False, "why": f"no lane pid {am.get('pid')}"}
    a_idx = read_idx(os.path.splitext(
        resolve(live_dir, lane["path"]))[0] + ".idx")
    if not a_idx:
        return {"present": False, "why": "empty audio idx"}
    n = max(1, math.ceil(secs / MPU_SECONDS))
    try:
        w = wave.open(wav_p)
        nf, nch = w.getnframes(), w.getnchannels()
    except (OSError, wave.Error) as e:
        return {"present": False, "why": f"wav unreadable ({e})"}
    # ANCHORED MAPPING (E35 addendum): wav sample 0 belongs to the
    # sidecar's first_seq fragment (--start-behind workers), so
    # wav_offset = (frag_index - frag_index(first_seq)) * SPF.
    j0 = next((j for j, (s, _, _) in enumerate(a_idx)
               if s == am.get("first_seq")), None)
    if j0 is None:
        w.close()
        return {"present": False, "why": "first_seq not in lane idx"}
    # the last N seconds of DECODED audio: fragments the wav already covers
    kmax = min(len(a_idx), j0 + nf // SPF)
    rows = list(enumerate(a_idx))[j0:kmax][-n:]
    if not rows:
        w.close()
        return {"present": False, "why": "wav shorter than one slot"}
    nonsilent = clipped_slots = 0
    clip_samples = total_samples = 0
    rmss = []
    for k, (seq, _, _) in rows:
        w.setpos((k - j0) * SPF)
        seg = np.frombuffer(w.readframes(SPF), "<i2").astype("f4")
        rms = float(np.sqrt((seg ** 2).mean())) if len(seg) else 0.0
        rmss.append(rms)
        if rms > RMS_SILENCE:
            nonsilent += 1
        nc = int((np.abs(seg) >= CLIP_LEVEL).sum())
        clip_samples += nc
        total_samples += len(seg)
        if nc:
            clipped_slots += 1
    w.close()
    return {
        "present": True, "pid": am.get("pid"), "channels": nch,
        "wav_secs": round(nf / 48000.0, 1),
        "slots_measured": len(rows),
        "nonsilent_frac": round(nonsilent / len(rows), 4),
        "clip_frac": round(clip_samples / total_samples, 6)
        if total_samples else 0.0,
        "rms_median": round(sorted(rmss)[len(rmss) // 2], 1),
    }


def parse_srt(path):
    out = []
    try:
        txt = open(path, encoding="utf-8", errors="replace").read()
    except OSError:
        return out
    pat = re.compile(r"(\d\d):(\d\d):(\d\d)[,.](\d{3})")

    def sec(g):
        return (int(g[0]) * 3600 + int(g[1]) * 60 + int(g[2])
                + int(g[3]) / 1000.0)

    for blk in txt.split("\n\n"):
        lines = [l for l in blk.strip().split("\n") if l]
        if len(lines) < 3 or "-->" not in lines[1]:
            continue
        m = pat.findall(lines[1])
        if len(m) < 2:
            continue
        out.append((sec(m[0]), sec(m[1]), "\n".join(lines[2:])))
    return out


def meter_cc(live_dir, lanes, mins):
    cues = parse_srt(os.path.join(live_dir, "live.srt"))
    if not cues:
        return {"cues_total": 0, "why": "no srt / no cues"}
    v = lane_of(lanes, "video")
    s = lane_of(lanes, "subs")
    sj = read_json(os.path.join(live_dir, "live_subs.json"))
    anchor = sj.get("anchor_seq") if sj else None
    source = "live_subs.json" if anchor is not None else "subs first_seq"
    if anchor is None and s:
        anchor = s.get("first_seq")
    srt_end = max(e for _, e, _ in cues)
    recent = [c for c in cues if c[1] >= srt_end - mins * 60.0]
    res = {"cues_total": len(cues), "cues_recent": len(recent),
           "window_mins": mins, "srt_span_s": round(srt_end, 1),
           "anchor_source": source, "anchor_seq": anchor}
    if anchor is None or not v:
        res["in_range_frac"] = None
        return res
    lo, hi = v["first_seq"], v["last_seq"] + 1
    inr = sum(1 for b, _, _ in recent
              if lo <= anchor + b / MPU_SECONDS < hi)
    res["lane_slots"] = [lo, hi]
    res["in_range_frac"] = round(inr / len(recent), 4) if recent else None
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--live-dir", required=True)
    ap.add_argument("--video-secs", type=float, default=60.0)
    ap.add_argument("--audio-secs", type=float, default=60.0)
    ap.add_argument("--cc-mins", type=float, default=5.0)
    a = ap.parse_args()
    d = a.live_dir
    lanes = (read_json(os.path.join(d, "live.json")) or {}).get("lanes", {})
    if not lanes:
        print(json.dumps({"dir": d, "error": "no live.json lanes"}))
        return 1

    vid = meter_video(d, lanes, a.video_secs)
    aud = {"eng": meter_audio_track(d, lanes, "live_audio.json",
                                    a.audio_secs),
           "spa": meter_audio_track(d, lanes, "live_audio_spa.json",
                                    a.audio_secs)}
    cc = meter_cc(d, lanes, a.cc_mins)

    print(f"VIDEO  last {vid.get('secs', '?')}s: "
          f"{vid.get('frames', 0)}/{vid.get('frames_expected', '?')} frames "
          f"({100 * vid.get('decode_complete', 0):.1f}%), "
          f"{vid.get('decode_errors', '?')} error lines, "
          f"null-sink {vid.get('null_sink_fps', '?')} fps")
    for k in ("eng", "spa"):
        t = aud[k]
        if t.get("present"):
            print(f"AUDIO  {k}: {t['slots_measured']} slots, "
                  f"non-silent {100 * t['nonsilent_frac']:.1f}%, "
                  f"clip {100 * t['clip_frac']:.3f}%, "
                  f"rms~{t['rms_median']} ({t['channels']}ch, "
                  f"{t['wav_secs']}s wav)")
        else:
            print(f"AUDIO  {k}: absent ({t.get('why')})")
    print(f"CC     {cc.get('cues_recent', 0)} cues in last "
          f"{a.cc_mins:.0f} min (of {cc.get('cues_total', 0)}), "
          f"anchor {cc.get('anchor_source')}, "
          f"in-range {cc.get('in_range_frac')}")
    print(json.dumps({"dir": d, "t": time.time(), "video": vid,
                      "audio": aud, "cc": cc}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
