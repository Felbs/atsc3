#!/usr/bin/env python3
"""atsc3_play.py -- watch the live file, N seconds behind the write head.

The other half of the decoupling. The chain writes DIR/live.m4s and never
waits for anything; this follows it at a deliberate lag, which is the whole
cushion. If decode stutters for two seconds, a viewer 20 seconds back never
finds out.

WHY YOU CANNOT JUST SEEK INTO IT
    live.m4s is fragmented MP4: `ftyp`+`moov` once, then `moof`+`mdat` per
    2.002 s MPU. A player needs the init boxes to know what the codecs are,
    and a fragment only decodes from its own `moof`. So starting mid-file
    means: send the init from offset 0, then jump to a real fragment
    BOUNDARY at the lag point -- never to an arbitrary byte, which produces
    a player that either refuses the stream or shows garbage.

ON NOT OUTRUNNING THE WALL CLOCK
    STVT learned this as pipe_throttle: dumping history into ffmpeg faster
    than real time makes output PTS jump ahead of the clock and the player
    skips. Here the burst is bounded by construction -- we start `lag`
    seconds back, so only `lag` seconds are ever delivered fast, which is
    exactly the prebuffer we want. After that the writer's own rate limits
    us. `--throttle` is kept for the case where that reasoning is wrong.

Usage:
    python tools/atsc3_play.py                      # default dir, 20 s behind
    python tools/atsc3_play.py --lag 30 --player mpv
    python tools/atsc3_play.py --out clip.mp4       # no player, just capture
"""
from __future__ import annotations

import argparse
import json
import os
import struct
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
MPU_SECONDS = 120 * 1001 / 60000.0          # 2.002


def walk(path, limit=None):
    """Yield (offset, type, size) for top-level boxes."""
    with open(path, "rb") as fh:
        off = 0
        end = limit if limit is not None else os.path.getsize(path)
        while off + 8 <= end:
            fh.seek(off)
            head = fh.read(8)
            if len(head) < 8:
                return
            size = struct.unpack(">I", head[:4])[0]
            typ = head[4:8]
            if size == 1:
                ext = fh.read(8)
                if len(ext) < 8:
                    return
                size = struct.unpack(">Q", ext)[0]
            elif size == 0:
                size = end - off
            if size < 8 or off + size > end:
                return
            yield off, typ, size
            off += size


def layout(path):
    """(init_end, [fragment offsets]) -- fragments start at each `moof`."""
    init_end, frags = 0, []
    for off, typ, size in walk(path):
        if typ in (b"ftyp", b"moov", b"styp"):
            init_end = off + size
        elif typ == b"moof":
            frags.append(off)
    return init_end, frags


def read_lanes(d):
    try:
        return json.load(open(os.path.join(d, "live.json"))).get("lanes", {})
    except (OSError, ValueError):
        return {}


def read_idx(path):
    """[(seq, off, len)] for a lane, from its sidecar index."""
    out = []
    try:
        for line in open(path):
            r = json.loads(line)
            out.append((r["seq"], r["off"], r["len"]))
    except (OSError, ValueError):
        pass
    return out


def align(d, lag_slots):
    """Pick the lanes and put them all on the SAME MPU slot.

    Every asset is carried in 2.002 s MPUs that share one sequence number,
    so alignment is arithmetic once each lane's slots are known. What it is
    NOT is "start each file at its beginning" -- lanes routinely begin on
    different slots (the caption lane loses nothing and so usually starts
    earlier than video), and lining up byte 0 with byte 0 silently offsets
    the whole programme. Nor is it "-shortest", which trims the ends and
    leaves the START misaligned, which is what the first attempt did.
    """
    lns = read_lanes(d)
    vid = next((l for l in lns.values() if l.get("kind") == "video"), None)
    if not vid:
        return None
    v_idx = read_idx(os.path.splitext(vid["path"])[0] + ".idx")
    if not v_idx:
        return None

    aud = ameta = None
    ap = os.path.join(d, "live_audio.json")
    if os.path.exists(ap):
        try:
            ameta = json.load(open(ap))
            if os.path.exists(ameta["path"]):
                aud = ameta
        except (OSError, ValueError):
            pass
    subs = os.path.join(d, "live.srt")
    sub_lane = next((l for l in lns.values() if l.get("kind") == "subs"), None)

    # start `lag_slots` back from the newest video fragment, but never
    # before a slot the audio can also supply
    start = v_idx[max(0, len(v_idx) - lag_slots)][0]
    if aud:
        start = max(start, aud["first_seq"])
    # the video fragment at or after that slot
    v_at = next((i for i, (s, _, _) in enumerate(v_idx) if s >= start),
                len(v_idx) - 1)
    start = v_idx[v_at][0]

    out = dict(video=vid, v_idx=v_idx, v_at=v_at, start=start,
               v_off=v_idx[v_at][1], init_bytes=vid.get("init_bytes", 0))
    if aud:
        out["audio"] = aud["path"]
        out["a_skip"] = max(0.0, (start - aud["first_seq"]) * MPU_SECONDS)
    if os.path.exists(subs) and sub_lane:
        out["subs"] = subs
        # the srt is relative to the SUBS lane's own first slot
        out["s_shift"] = (sub_lane.get("first_seq", start) - start) * MPU_SECONDS
    return out


def player_cmd(kind, title, extra):
    if kind == "none":
        return None
    if kind == "ffplay":
        c = ["ffplay", "-hide_banner", "-loglevel", "error", "-nostats",
             "-autoexit", "-window_title", title, "-i", "pipe:0"]
    elif kind == "mpv":
        c = ["mpv", "--title=" + title, "--cache=yes", "-"]
    else:
        c = kind.split()
    return c + (extra.split() if extra else [])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--live-dir", default=None)
    ap.add_argument("--lag", type=float, default=20.0,
                    help="seconds behind the write head (default 20)")
    ap.add_argument("--player", default="ffplay",
                    choices=("ffplay", "mpv", "none"))
    ap.add_argument("--player-args", default=None)
    ap.add_argument("--out", default=None, help="also write the stream here")
    ap.add_argument("--throttle", type=float, default=0.0,
                    help="cap delivery at N MB/s (0 = off; see module docs)")
    ap.add_argument("--wait", type=float, default=60.0,
                    help="seconds to wait for the chain to produce fragments")
    a = ap.parse_args()
    d = a.live_dir or os.path.join(ROOT, "data", "atsc3_live")
    # The video lane's name comes from live.json, not from convention. The
    # writer names lanes by pid (live_video_pid12.m4s) precisely because this
    # multiplex carries two audio assets; the old fixed "live.m4s" survives
    # only as a fallback for pre-multi-lane recordings. Found 8/07 when the
    # tailing player refused a healthy chain the moment the band reopened.
    path = os.path.join(d, "live.m4s")
    vid = next((l for l in read_lanes(d).values()
                if l.get("kind") == "video"), None)
    if vid and os.path.exists(vid["path"]):
        path = vid["path"]

    # Wait for enough file to exist to start at the requested lag.
    need = max(2, int(a.lag / MPU_SECONDS))
    t0 = time.time()
    while True:
        if os.path.exists(path):
            init_end, frags = layout(path)
            if init_end and len(frags) >= 2:
                if len(frags) >= need or time.time() - t0 > a.wait:
                    break
        if time.time() - t0 > a.wait:
            if not os.path.exists(path):
                print(f"no live file at {path} -- is the chain running?",
                      file=sys.stderr)
                return 1
            break
        time.sleep(1.0)

    init_end, frags = layout(path)
    if not init_end or not frags:
        print("live file has no init or no fragments yet", file=sys.stderr)
        return 1

    start_idx = max(0, len(frags) - need)
    start = frags[start_idx]
    behind = (len(frags) - start_idx) * MPU_SECONDS
    print(f"live file {path}")
    print(f"  {len(frags)} fragments, init {init_end} B")
    print(f"  starting {behind:.1f} s behind the write head "
          f"(fragment {start_idx}/{len(frags)})")

    cmd = player_cmd(a.player, "ATSC 3.0 -- live", a.player_args)
    proc = None
    if cmd:
        proc = subprocess.Popen(cmd, stdin=subprocess.PIPE)
        sink = proc.stdin
    else:
        sink = None
    out_fh = open(a.out, "wb") if a.out else None

    def emit(b):
        if sink:
            sink.write(b); sink.flush()
        if out_fh:
            out_fh.write(b)

    cap = a.throttle * 1e6 if a.throttle else 0.0
    sent, t_start = 0, time.time()
    try:
        with open(path, "rb") as fh:
            fh.seek(0)
            emit(fh.read(init_end))          # the init boxes, always first
            fh.seek(start)
            idle = 0.0
            while True:
                chunk = fh.read(1 << 16)
                if not chunk:
                    if proc and proc.poll() is not None:
                        print("  player exited")
                        break
                    time.sleep(0.2)
                    idle += 0.2
                    if idle > 30:
                        print("  no new data for 30 s -- chain stopped?")
                        break
                    continue
                idle = 0.0
                emit(chunk)
                sent += len(chunk)
                if cap:
                    want = sent / cap
                    slip = want - (time.time() - t_start)
                    if slip > 0:
                        time.sleep(slip)
    except (BrokenPipeError, KeyboardInterrupt):
        print("  stopped")
    finally:
        if sink:
            try:
                sink.close()
            except OSError:
                pass
        if out_fh:
            out_fh.close()
        if proc and proc.poll() is None:
            proc.terminate()
    print(f"  delivered {sent/1e6:.1f} MB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
