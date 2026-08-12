#!/usr/bin/env python3
"""atsc3_playd.py -- keep the live window ALIVE through a churning band.

The night of 8/07 was the stress test nobody scheduled: RF33 opened and
closed on ~10 minute cycles for hours. Under that, the tailing player died
in two different ways and stayed dead in both:

  * `atsc3_play` exits after 30 s of no new data ("chain stopped?"). Right
    call for a dead chain -- wrong ending for a FADE, because nothing brings
    the window back when the band reopens.
  * every supervisor restart ROLLS the lane, and a player tailing the old
    file is holding a cursor into a file that stopped growing (the audio
    worker's E24 lesson, third appearance).

So: a watchdog for the viewing experience, symmetric to the chain's
supervisor. It (re)launches the tailing player ONLY when the lane is
actually advancing -- spawning a player into dead air would open a window
onto a frozen tail every 20 s -- and replaces it when the generation
changes underneath it.

Usage:
    python tools/atsc3_playd.py --live-dir data/e29
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def log(m):
    print(f"{time.strftime('%H:%M:%S')} {m}", flush=True)


def video_lane(live_dir):
    try:
        d = json.load(open(os.path.join(live_dir, "live.json")))
        for ln in d.get("lanes", {}).values():
            if ln.get("kind") == "video":
                return ln
    except (OSError, ValueError):
        pass
    return None


def kill_tree(proc):
    """The player is python + a spawned ffplay; take both."""
    if proc is None or proc.poll() is not None:
        return
    subprocess.run(["taskkill", "/PID", str(proc.pid), "/T", "/F"],
                   capture_output=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--live-dir", required=True)
    ap.add_argument("--lag", type=float, default=30.0)
    ap.add_argument("--check", type=float, default=15.0)
    a = ap.parse_args()

    proc = None
    prev_last = prev_gen = None
    launches = 0
    log(f"watching {a.live_dir}; player relaunches while the lane advances")
    try:
        while True:
            ln = video_lane(a.live_dir)
            if ln:
                last, gen = ln.get("last_seq"), ln.get("generation", 0)
                advancing = (last is not None and last != prev_last)
                alive = proc is not None and proc.poll() is None
                rolled = alive and prev_gen is not None and gen != prev_gen
                if advancing and (not alive or rolled):
                    kill_tree(proc)
                    launches += 1
                    why = "lane rolled" if rolled else \
                        ("player exited" if launches > 1 else "start")
                    log(f"launching player #{launches} ({why}; "
                        f"generation {gen})")
                    proc = subprocess.Popen(
                        [sys.executable, "-u",
                         os.path.join(HERE, "atsc3_play.py"),
                         "--live-dir", a.live_dir,
                         "--lag", str(a.lag), "--player", "ffplay"],
                        cwd=ROOT,
                        stdout=open(os.path.join(a.live_dir, "playd_child.log"),
                                    "ab", buffering=0),
                        stderr=subprocess.STDOUT)
                if last is not None:
                    prev_last = last
                if advancing:
                    prev_gen = gen
            time.sleep(a.check)
    finally:
        kill_tree(proc)


if __name__ == "__main__":
    sys.exit(main())
