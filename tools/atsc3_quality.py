#!/usr/bin/env python3
"""atsc3_quality.py -- the measuring apparatus for a live ATSC 3.0 viewing stack.

One sample = one row, every stage graded from what the user actually gets:

  RF/chain   FEC %, real-time factor, SDR overflows, MPUs lost      (chain.log)
  video      codec, WxH, fps, HDR transfer, decoded fps over the last
             10 s of the muxed TS, decode errors, FROZEN? (freezedetect),
             black? (blackdetect)                                     (_tv/live_tv2.ts)
  audio      eng: lane, channels, bad frames, speed, per-channel RMS
             of the last 10 s (the CENTRE channel is the dialogue:
             a 5.1 mix whose centre is silent has no dialogue), silence
             spa: the same for the stereo pair                        (aud_*.log, *.wav)
  captions   cues in the last 60 s, readable/raw ratio, age of the
             newest cue (stale captions = frozen caption path)        (subs.log, live.srt)
  player     VLC state, playhead advancing?, mux stall/dropped/lag    (VLC http, tv.log)

Usage:
  python tools/atsc3_quality.py --live-dir ~/atsc3_live_rf33            # one sample
  python tools/atsc3_quality.py --live-dir ~/atsc3_live_rf33 --watch 30 # every 30 s
  ... --json out.jsonl                                                  # append every sample

Grades: PASS / WARN / FAIL per row and a one-line verdict. Thresholds are at the
top of the file; every number printed can be re-derived by hand from the logs.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import wave
from datetime import datetime

import numpy as np

# ---- thresholds ------------------------------------------------------------
FEC_WARN, FEC_FAIL = 99.9, 99.0          # % converged, latest status line
RT_WARN = 0.98                           # sustained x real time (live caps at 1.0)
FPS_RATIO_WARN = 0.95                    # decoded fps / nominal fps over the window
DIALOG_DBFS_FAIL = -45.0                 # centre-channel RMS below this = no dialogue
SILENCE_DBFS = -60.0                     # a channel below this = silent
CC_STALE_S = 90.0                        # newest cue older than this = caption path stuck
CC_MIN_PER_MIN = 2                       # fewer readable cues/min than this = WARN (ads, music)
WINDOW_S = 10.0                          # analysis window at the tail of each lane


def sh(cmd, timeout=60):
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    return p.returncode, p.stdout, p.stderr


def dbfs(x):
    r = float(np.sqrt(np.mean(np.square(x.astype(np.float64))))) if x.size else 0.0
    return 20 * np.log10(max(r, 1e-9) / 32768.0)


# ---- chain -----------------------------------------------------------------
STATUS_RE = re.compile(r"^\[(\d\d:\d\d:\d\d)\]\s+([\d.]+)s\s+(\d+) Frames\s+([\d.]+)x rt.*?FEC (\d+)/(\d+)", re.M)


def chain(live):
    out = {"status": "n/a"}
    p = os.path.join(live, "chain.log")
    if not os.path.exists(p):
        return out
    txt = open(p, errors="replace").read()
    last = None
    for m in STATUS_RE.finditer(txt):
        last = m
    out["overflows"] = txt.count("sample continuity BREAK")
    out["reacquisitions"] = txt.count("bootstrap acquired") - 1
    bt = re.search(r"bias-T: wrote (\w+), read back (\w+)", txt)
    out["bias_t_readback"] = bt.group(2) if bt else None
    g = re.search(r"gain confirmed by readback: (\S+ \S+)", txt)
    out["gain"] = g.group(1) if g else None
    if last:
        ok, tot = int(last.group(5)), int(last.group(6))
        out.update(t_air=float(last.group(2)), frames=int(last.group(3)),
                   x_rt=float(last.group(4)), fec_ok=ok, fec_total=tot,
                   fec_pct=(100.0 * ok / tot) if tot else None, status="ok")
    return out


# ---- video -----------------------------------------------------------------
def tail_bytes(path, n, tmp):
    sz = os.path.getsize(path)
    with open(path, "rb") as f:
        f.seek(max(0, sz - n))
        data = f.read()
    i = data.find(b"\x47")
    data = data[i:]
    # Start at a CHUNK boundary, not mid picture-group: the viewer's mux writes
    # a PAT (PID 0, payload_unit_start) at the head of every chunk, and every
    # chunk begins on an MPU = a keyframe. A cut anywhere else makes HEVC spew
    # reference errors that are artefacts of the cut, not of the broadcast.
    for k in range(1, len(data) // 188):
        p = data[k * 188: k * 188 + 4]
        if len(p) == 4 and p[0] == 0x47 and (p[1] & 0x40) and ((p[1] & 0x1F) << 8 | p[2]) == 0:
            data = data[k * 188:]
            break
    with open(tmp, "wb") as g:
        g.write(data)
    return len(data)


def lane_tail(live, window, tmp):
    """The last `window` seconds of the chain's own video lane (fMP4: init
    segment + whole fragments, every fragment starting on an MPU keyframe), so
    a decode of it has NO cut artefacts: every error is the broadcast's or the
    chain's. Returns the lane dict from live.json or None."""
    lj = json.load(open(os.path.join(live, "live.json")))
    lane = next((v for v in lj["lanes"].values() if v.get("kind") == "video"), None)
    if not lane or not os.path.exists(lane["path"]):
        return None
    p, ib = lane["path"], int(lane["init_bytes"])
    rate = max(100_000, int(lane["bytes"] / max(lane["media_s"], 1.0)))
    sz = os.path.getsize(p)
    with open(p, "rb") as f:
        init = f.read(ib)
        f.seek(max(ib, sz - int(window * rate * 1.2)))
        tail = f.read()
    i = tail.find(b"moof")
    if i < 4:
        return None
    with open(tmp, "wb") as g:
        g.write(init + tail[i - 4:])
    return lane


def video(live, window):
    out = {"status": "n/a"}
    tmp = os.path.join(tempfile.gettempdir(), "atsc3q_lane.mp4")
    try:
        lane = lane_tail(live, window, tmp)
    except Exception:                                          # noqa: BLE001
        lane = None
    if lane is None:
        return out
    out["lane"] = os.path.basename(lane["path"])
    out["segments"] = lane.get("segments"); out["media_s"] = round(lane.get("media_s", 0), 1)
    rc, so, _ = sh(["ffprobe", "-v", "error", "-select_streams", "v:0", "-show_entries",
                    "stream=codec_name,width,height,r_frame_rate,color_transfer,color_primaries,pix_fmt",
                    "-of", "json", tmp])
    try:
        s = json.loads(so)["streams"][0]
        num, den = s["r_frame_rate"].split("/")
        out.update(codec=s.get("codec_name"), w=s.get("width"), h=s.get("height"),
                   fps_nominal=round(int(num) / int(den), 3), transfer=s.get("color_transfer"),
                   primaries=s.get("color_primaries"), pix_fmt=s.get("pix_fmt"))
        out["hdr"] = s.get("color_transfer") in ("smpte2084", "arib-std-b67")
    except Exception:                                          # noqa: BLE001
        return out
    # decode the window: count frames, errors, freezes, black
    rc, so, se = sh(["ffmpeg", "-v", "error", "-nostdin", "-i", tmp, "-map", "0:v:0",
                     "-vf", "freezedetect=n=-50dB:d=2,blackdetect=d=1.0:pic_th=0.98",
                     "-f", "null", "-"], timeout=120)
    frames_m = re.findall(r"frame=\s*(\d+)", se)
    # -v error hides the frame= progress; count frames with ffprobe instead
    _, so2, _ = sh(["ffprobe", "-v", "error", "-select_streams", "v:0", "-count_frames",
                    "-show_entries", "stream=nb_read_frames,duration", "-of", "json", tmp], timeout=120)
    try:
        st = json.loads(so2)["streams"][0]
        n = int(st.get("nb_read_frames", 0))
    except Exception:                                          # noqa: BLE001
        n = 0
    # fMP4 fragments carry absolute broadcast timestamps, so ffprobe's duration
    # is meaningless here; the window covered is frames / nominal fps, and a
    # fragment that decoded short of its nominal frame count shows up as
    # fps_decoded < fps_nominal.
    frags = tmp and open(tmp, "rb").read().count(b"moof")
    dur = (lane["media_s"] / max(lane["segments"], 1)) * frags if lane.get("segments") else 0.0
    out.update(window_s=round(dur, 2), fragments=frags, frames_decoded=n,
               fps_decoded=round(n / dur, 2) if dur else None)
    out["decode_errors"] = len([l for l in se.splitlines() if l.strip() and "freeze" not in l and "black" not in l])
    out["freezes"] = len(re.findall(r"freeze_start", se))
    out["black"] = len(re.findall(r"black_start", se))
    out["status"] = "ok"
    return out


# ---- audio -----------------------------------------------------------------
WORKER_RE = re.compile(r"(\d+) fr \((\d+) bad\), ([\d.]+) s of audio.*?([\d.]+)x real time")


def audio(live, log, wav, window):
    out = {"status": "n/a"}
    lp, wp = os.path.join(live, log), os.path.join(live, wav)
    if os.path.exists(lp):
        txt = open(lp, errors="replace").read()
        lane = re.search(r"audio lane pid (\d+)", txt)
        out["lane_pid"] = int(lane.group(1)) if lane else None
        last = None
        for m in WORKER_RE.finditer(txt):
            last = m
        if last:
            out.update(frames=int(last.group(1)), bad=int(last.group(2)),
                       audio_s=float(last.group(3)), x_rt=float(last.group(4)))
    if not os.path.exists(wp) or os.path.getsize(wp) < 1000:
        return out
    with open(wp, "rb") as f:
        hdr = f.read(44)
    ch = int.from_bytes(hdr[22:24], "little"); fs = int.from_bytes(hdr[24:28], "little")
    sz = os.path.getsize(wp) - 44
    n = int(window * fs)
    with open(wp, "rb") as f:
        f.seek(44 + max(0, sz - n * ch * 2))
        pcm = np.frombuffer(f.read(n * ch * 2), dtype="<i2")
    pcm = pcm[: (len(pcm) // ch) * ch].reshape(-1, ch)
    names = ["L", "R", "C", "LFE", "Ls", "Rs"] if ch == 6 else ["L", "R"]
    lv = {names[i]: round(float(dbfs(pcm[:, i])), 1) for i in range(ch)}
    out.update(channels=ch, fs=fs, wav_s=round(sz / (ch * 2 * fs), 1), levels_dbfs=lv,
               silent_channels=[k for k, v in lv.items() if v < SILENCE_DBFS])
    if ch == 6:
        out["dialogue_dbfs"] = lv["C"]
    out["status"] = "ok"
    return out


# ---- captions --------------------------------------------------------------
def captions(live):
    out = {"status": "n/a"}
    lp, sp = os.path.join(live, "subs.log"), os.path.join(live, "live.srt")
    if os.path.exists(lp):
        m = None
        for m in re.finditer(r"(\d+) documents, (\d+) raw cues -> (\d+) readable", open(lp, errors="replace").read()):
            pass
        if m:
            out.update(documents=int(m.group(1)), raw=int(m.group(2)), readable=int(m.group(3)))
    if not os.path.exists(sp):
        return out
    txt = open(sp, errors="replace").read()
    cues = re.findall(r"(\d\d):(\d\d):(\d\d),(\d\d\d) --> ", txt)
    ends = [int(h) * 3600 + int(m) * 60 + int(s) + int(ms) / 1000 for h, m, s, ms in cues]
    out["cues_total"] = len(ends)
    out["srt_age_s"] = round(time.time() - os.path.getmtime(sp), 1)
    if ends:
        last = max(ends)
        out["cues_last_60s"] = sum(1 for e in ends if e >= last - 60)
        out["last_cue_text"] = txt.strip().splitlines()[-1][:60]
    out["status"] = "ok"
    return out


# ---- player ----------------------------------------------------------------
def player(live):
    out = {"status": "n/a"}
    tp = os.path.join(live, "tv.log")
    if os.path.exists(tp):
        tel = None
        for tel in re.finditer(r"tel: (\d+)s out .*?stall<=([\d.]+)s\s+dropped (\d+)\s+fails (\d+)\s+lag (\d+)s",
                               open(tp, errors="replace").read()):
            pass
        if tel:
            out.update(out_s=int(tel.group(1)), stall_s=float(tel.group(2)), dropped=int(tel.group(3)),
                       fails=int(tel.group(4)), lag_s=int(tel.group(5)))
    # Is anything being DRAWN? With the monitor powered off (DPMS) the compositor
    # stops presenting VLC's window and almost every picture counts as lost --
    # a measurement of the dark screen, not of the receiver. Record it.
    if shutil.which("hyprctl"):
        _, mon, _ = sh(["hyprctl", "monitors"], timeout=10)
        st = re.findall(r"dpmsStatus: (\d)", mon)
        out["display_on"] = any(s == "1" for s in st) if st else None
    rc, so, _ = sh(["pgrep", "-xa", "vlc"])
    out["vlc_headless"] = "--vout=vdummy" in so
    m = re.search(r"--http-port (\d+) --http-password (\S+)", so)
    if m and shutil.which("curl"):
        def status():
            rc, body, _ = sh(["curl", "-s", "-u", f":{m.group(2)}",
                              f"http://127.0.0.1:{m.group(1)}/requests/status.xml"], timeout=10)
            st = re.search(r"<state>(\w+)", body); t = re.search(r"<time>(\d+)", body)
            return (st.group(1) if st else None), (int(t.group(1)) if t else None)
        s1, t1 = status(); time.sleep(3); s2, t2 = status()
        out.update(vlc_state=s2, vlc_time=t2, playhead_advancing=(t1 is not None and t2 is not None and t2 > t1))
        # VLC's own counters = what reached the screen: pictures displayed vs
        # lost (late/dropped), demux corruption/discontinuities, audio buffers.
        _, body, _ = sh(["curl", "-s", "-u", f":{m.group(2)}",
                         f"http://127.0.0.1:{m.group(1)}/requests/status.xml"], timeout=10)
        st = {k: int(v) for k, v in re.findall(r"<(displayedpictures|lostpictures|demuxcorrupted|demuxdiscontinuity|lostabuffers|playedabuffers)>(\d+)", body)}
        out["vlc_stats"] = st
        if st.get("displayedpictures"):
            out["vlc_lost_pct"] = round(100.0 * st.get("lostpictures", 0) / (st["displayedpictures"] + st.get("lostpictures", 0)), 3)
        hw = sum(1 for l in open(tp, errors="replace") if "hardware acceleration picture allocation failed" in l) if os.path.exists(tp) else 0
        out["vlc_hw_decode_failures"] = hw
        out["status"] = "ok"
    else:
        out["vlc_state"] = "not running"
    return out


# ---- grading ---------------------------------------------------------------
def grade(r):
    g = {}
    c = r["chain"]
    if c.get("status") == "ok":
        g["chain"] = "FAIL" if (c["fec_pct"] or 0) < FEC_FAIL else "WARN" if (c["fec_pct"] or 0) < FEC_WARN or c["x_rt"] < RT_WARN else "PASS"
    v = r["video"]
    if v.get("status") == "ok":
        ratio = (v["fps_decoded"] or 0) / (v["fps_nominal"] or 1)
        g["video"] = "FAIL" if v["freezes"] or v["decode_errors"] > 5 or ratio < 0.8 else "WARN" if ratio < FPS_RATIO_WARN or v["black"] else "PASS"
    a = r["audio_eng"]
    if a.get("status") == "ok":
        bad = a.get("bad", 0); dia = a.get("dialogue_dbfs")
        g["audio_eng"] = "FAIL" if bad > 0 and bad >= 0.05 * max(a.get("frames", 1), 1) or (dia is not None and dia < DIALOG_DBFS_FAIL) \
            else "WARN" if bad or a.get("channels") != 6 or a.get("x_rt", 9) < 1.2 else "PASS"
    s = r["audio_spa"]
    if s.get("status") == "ok":
        g["audio_spa"] = "FAIL" if s.get("bad", 0) >= 0.05 * max(s.get("frames", 1), 1) else "WARN" if s.get("bad") or len(s.get("silent_channels", [])) else "PASS"
    cc = r["captions"]
    if cc.get("status") == "ok":
        g["captions"] = "FAIL" if cc["srt_age_s"] > CC_STALE_S else "WARN" if cc.get("cues_last_60s", 0) < CC_MIN_PER_MIN else "PASS"
    p = r["player"]
    if p.get("status") == "ok":
        if p.get("display_on") is False and not p.get("vlc_headless"):
            g["player"] = "DARK"      # monitor off: the player cannot be graded, and is not the receiver's fault
        else:
            g["player"] = "FAIL" if p["vlc_state"] != "playing" or not p["playhead_advancing"] or p.get("vlc_hw_decode_failures") \
                else "WARN" if p.get("stall_s", 0) > 0 or p.get("dropped", 0) else "PASS"
    worst = "PASS"
    for v in g.values():
        if v == "FAIL": worst = "FAIL"
        elif v == "WARN" and worst == "PASS": worst = "WARN"
    return g, worst


def sample(live, window):
    r = {"time": datetime.now().isoformat(timespec="seconds"),
         "chain": chain(live), "video": video(live, window),
         "audio_eng": audio(live, "aud_eng.log", "live_audio.wav", window),
         "audio_spa": audio(live, "aud_spa.log", "live_audio_spa.wav", window),
         "captions": captions(live), "player": player(live)}
    r["grade"], r["verdict"] = grade(r)
    return r


def show(r):
    g = r["grade"]; c, v, a, s, cc, p = (r[k] for k in ("chain", "video", "audio_eng", "audio_spa", "captions", "player"))
    print(f"=== {r['time']}  VERDICT: {r['verdict']}")
    if c.get("status") == "ok":
        print(f"  chain    [{g.get('chain','-'):4}] FEC {c['fec_pct']:.3f}% ({c['fec_ok']}/{c['fec_total']})  {c['x_rt']:.2f}x rt  "
              f"overflows {c['overflows']}  re-acq {c['reacquisitions']}  bias-T readback {c['bias_t_readback']}  gain {c['gain']}")
    if v.get("status") == "ok":
        print(f"  video    [{g.get('video','-'):4}] {v['codec']} {v['w']}x{v['h']} {v['fps_nominal']} fps  HDR {v['hdr']} ({v['transfer']})  "
              f"decoded {v['fps_decoded']} fps over {v['window_s']} s  errors {v['decode_errors']}  freezes {v['freezes']}  black {v['black']}")
    if a.get("status") == "ok":
        print(f"  eng 5.1  [{g.get('audio_eng','-'):4}] lane pid {a.get('lane_pid')}  {a.get('channels')} ch  bad {a.get('bad')}/{a.get('frames')}  "
              f"{a.get('x_rt')}x rt  levels dBFS {a['levels_dbfs']}  DIALOGUE(C) {a.get('dialogue_dbfs')}")
    if s.get("status") == "ok":
        print(f"  spa      [{g.get('audio_spa','-'):4}] lane pid {s.get('lane_pid')}  bad {s.get('bad')}/{s.get('frames')}  {s.get('x_rt')}x rt  levels {s['levels_dbfs']}")
    if cc.get("status") == "ok":
        print(f"  captions [{g.get('captions','-'):4}] {cc.get('readable')} readable / {cc.get('raw')} raw  last 60 s: {cc.get('cues_last_60s')}  "
              f"srt age {cc['srt_age_s']} s  last: {cc.get('last_cue_text','')!r}")
    if p.get("status") == "ok":
        st = p.get("vlc_stats", {})
        print(f"  player   [{g.get('player','-'):4}] VLC {p['vlc_state']} t={p['vlc_time']} advancing {p['playhead_advancing']}  "
              f"pictures shown {st.get('displayedpictures')} lost {st.get('lostpictures')} ({p.get('vlc_lost_pct')}%)  "
              f"demux corrupt {st.get('demuxcorrupted')} discont {st.get('demuxdiscontinuity')}  audio buffers lost {st.get('lostabuffers')}/{st.get('playedabuffers')}  "
              f"hw-decode failures {p['vlc_hw_decode_failures']}  mux stall {p.get('stall_s')} s dropped {p.get('dropped')} lag {p.get('lag_s')} s")
    else:
        print(f"  player   [----] {p.get('vlc_state')}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--live-dir", default=os.path.expanduser("~/atsc3_live_rf33"))
    ap.add_argument("--watch", type=float, default=0, help="repeat every N seconds")
    ap.add_argument("--json", default="", help="append each sample as one JSON line")
    ap.add_argument("--window", type=float, default=WINDOW_S)
    a = ap.parse_args(argv)
    prev = None
    while True:
        r = sample(os.path.abspath(os.path.expanduser(a.live_dir)), a.window)
        show(r)
        # what the viewer lost SINCE THE LAST SAMPLE: the number that tracks a
        # fix (chunk size, player options), not the lifetime total
        if prev and prev["player"].get("vlc_stats") and r["player"].get("vlc_stats"):
            p0, p1 = prev["player"]["vlc_stats"], r["player"]["vlc_stats"]
            shown = p1.get("displayedpictures", 0) - p0.get("displayedpictures", 0)
            lost = p1.get("lostpictures", 0) - p0.get("lostpictures", 0)
            disc = p1.get("demuxdiscontinuity", 0) - p0.get("demuxdiscontinuity", 0)
            r["player"]["since_last"] = dict(shown=shown, lost=lost, discontinuities=disc,
                                             lost_pct=round(100.0 * lost / max(shown + lost, 1), 3))
            print(f"           since last sample: shown {shown}  lost {lost} ({r['player']['since_last']['lost_pct']}%)  "
                  f"demux discontinuities {disc}")
        prev = r
        if a.json:
            with open(a.json, "a") as f:
                f.write(json.dumps(r) + "\n")
        if not a.watch:
            return 0 if r["verdict"] == "PASS" else 1
        time.sleep(a.watch)


if __name__ == "__main__":
    sys.exit(main())
