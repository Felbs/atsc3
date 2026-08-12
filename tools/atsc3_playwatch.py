#!/usr/bin/env python3
"""atsc3_playwatch.py -- restart the player when the PICTURE stops advancing.

THE GAP THIS FILLS (E93, 8/11).  Fox wedged twice in one hour with every
instrument we own reading green:

    chain FEC        100%, phase locked, 60 min
    video lane       +3.77 MB / 10 s
    viewer TS sink   +11 MB / 20 s
    VLC process      alive
    the screen       FROZEN

`atsc3_tv` detects the player DYING.  `atsc3_tvwatch` detects the VIEWER
wedging -- its test is "lanes advancing AND TS sink stale", which was
correctly silent here because the sink was healthy.  Neither watches the
last hop, and the last hop is the only one the user actually sees.  The
wedge lived downstream of every probe in the stack, so the symptom could
only reach us by the user saying "it froze again".

THE PRIMITIVE -- CORRECTED 8/11 22:45 BY MEASUREMENT.  The first version of
this file watched VLC's <time>, and that was WRONG.  Measured during a real
freeze, 25 s apart:

    time                +25 s     <- the clock runs
    decodedaudio      +2072       <- audio decodes
    decodedvideo          +0      <- THE VIDEO DECODER HAS STOPPED
    displayedpictures     +0
    lostpictures          +0      <- not starved, not dropping: STOPPED

A stalled video decoder does not stall the clock, so a position-based watch
sails straight through the exact failure it was built for.  The honest
primitive is the PICTURE ITSELF: <displayedpictures>.  Frames displayed is
the only number that means "the user can see motion".

VLC's own HTTP interface reports playback position AND frame counters, and
`atsc3_tv.record_player` already writes the port and password to
`_tv/player.pid` -- the docstring there even calls it "warden's
player-progress probe".  The measurement was always available; nothing was
reading it.  A picture that is advancing has a frame count that increases.

THE TEST, and why both halves are required:

    FRAMES STALLED  AND  TS GROWING -> the player is wedged      RESTART
    FRAMES STALLED  AND  TS STALE   -> the feed stopped; the player is
                                       correctly showing the last frame  LEAVE
    FRAMES ADVANCING                -> healthy                    LEAVE

Restarting on stalled frames alone would thrash the player through every
genuine fade, which is the same mistake tvwatch avoids by pairing its two
signals.  A wedge is only a wedge when there is fresh material to play.

Bounded, like every other guard here: at most MAX_PER_HOUR restarts, so a
player that cannot survive its input degrades to a frozen picture and a log
line rather than an infinite respawn storm.

    python tools/atsc3_playwatch.py <live_dir>
"""
from __future__ import annotations

import base64
import json
import os
import re
import sys
import time
import urllib.request

# 20 s, not 60: measured 8/11, the wedge recurs every 7-12 min with a HEALTHY
# chain underneath, so the cost of waiting is a minute of frozen picture per
# occurrence. 20 s still clears any legitimate stutter by a wide margin.
STALL_NEED = 20.0
POLL = 20.0
MAX_PER_HOUR = 10   # the wedge recurs faster than 6/h; a cap below the
                    # real rate silently stops protecting exactly when needed


def log(m):
    print(f"{time.strftime('%H:%M:%S')} {m}", flush=True)


def player_rec(live_dir):
    try:
        return json.load(open(os.path.join(live_dir, "_tv", "player.pid")))
    except (OSError, ValueError):
        return None


def vlc_time(rec):
    """-> (state, displayed_pictures) from VLC, or None.

    Returns the DISPLAYED PICTURE COUNT, not the clock: a stalled video
    decoder leaves the clock advancing (measured), so the clock cannot see
    this wedge at all.

    None means UNREADABLE, which is deliberately not the same as stalled:
    a player that has just been respawned, or is mid-start, must not be
    counted as frozen.  Absence of evidence is not evidence of a wedge.
    """
    if not rec or "port" not in rec:
        return None
    url = f"http://127.0.0.1:{rec['port']}/requests/status.xml"
    auth = base64.b64encode(f":{rec['pw']}".encode()).decode()
    req = urllib.request.Request(url, headers={"Authorization": f"Basic {auth}"})
    try:
        x = urllib.request.urlopen(req, timeout=8).read().decode("utf-8", "replace")
    except Exception:                                        # noqa: BLE001
        return None
    st = re.search(r"<state>([^<]*)</state>", x)
    dp = re.search(r"<displayedpictures>([^<]*)</displayedpictures>", x)
    if not dp:
        return None
    try:
        return (st.group(1) if st else "?", float(dp.group(1)))
    except ValueError:
        return None


def ts_size(live_dir):
    p = os.path.join(live_dir, "_tv", "live_tv2.ts")
    try:
        return os.path.getsize(p)
    except OSError:
        return None


def restart_player(live_dir):
    """Kill VLC; atsc3_tv notices the death and respawns at the live edge."""
    rec = player_rec(live_dir)
    if not rec or "pid" not in rec:
        log("  cannot restart: no player record")
        return False
    os.system(f'taskkill /PID {rec["pid"]} /T /F >nul 2>&1')
    return True


def main():
    live_dir = sys.argv[1] if len(sys.argv) > 1 else "data/foxlive"
    log(f"play-watchdog on {live_dir}: restart the player when DISPLAYED "
        f"FRAMES stop advancing while the TS still grows "
        f"(stall >= {STALL_NEED:.0f}s)")
    t_prev = None
    ts_prev = None
    stalled_since = None
    restarts = []
    while True:
        time.sleep(POLL)
        rec = player_rec(live_dir)
        v = vlc_time(rec)
        ts = ts_size(live_dir)
        if v is None or ts is None:
            t_prev, ts_prev, stalled_since = None, ts, None
            continue
        state, t = v
        ts_growing = ts_prev is not None and ts > ts_prev
        # frames displayed must INCREASE; equal means nothing reached the screen
        pos_stalled = t_prev is not None and t <= t_prev
        now = time.time()
        if pos_stalled and ts_growing and state != "paused":
            stalled_since = stalled_since or now
            held = now - stalled_since
            if held >= STALL_NEED:
                restarts = [r for r in restarts if now - r < 3600]
                if len(restarts) >= MAX_PER_HOUR:
                    log(f"WEDGE at {t:.0f} frames but {len(restarts)} restarts "
                        f"already this hour -- standing down rather than "
                        f"thrashing the player")
                    stalled_since = None
                    continue
                log(f"WEDGE: displayed frames stuck at {t:.0f} for "
                    f"{held:.0f}s while the TS grew {ts - ts_prev} bytes "
                    f"-- restarting player")
                if restart_player(live_dir):
                    restarts.append(now)
                stalled_since = None
                t_prev, ts_prev = None, None
                time.sleep(40)          # let the viewer respawn and settle
                continue
        else:
            stalled_since = None
        t_prev, ts_prev = t, ts


if __name__ == "__main__":
    sys.exit(main())
