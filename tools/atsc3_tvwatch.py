#!/usr/bin/env python3
"""atsc3_tvwatch.py -- restart the viewer when it WEDGES ALIVE.

The viewer's poll loop can hang across a long fade + multi-restart roll
(E42, 8/08): the process stays alive, the band recovers, the lanes grow, but
the muxer's TS sink never resumes. A fresh start always recovers cleanly, so
until the root-cause poll-loop wedge is fixed in daylight, this external
supervisor makes recovery automatic instead of manual -- the STVT pattern:
watch a liveness primitive, restart on a wedge, bounded.

LIVENESS = the pair (lanes advancing) AND (TS sink growing). A TS that is
stale while the video lane's last_seq keeps climbing is the wedge; a stale
TS while the lane is ALSO stale is just a fade -- correct, not a wedge. Only
the first triggers a restart, so a genuine dead band is never thrashed.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = sys.executable
IS_WIN = (os.name == "nt")


def _proc_cmdline(pid):
    try:
        return open(f"/proc/{pid}/cmdline", "rb").read() \
            .replace(b"\0", b" ").decode("utf-8", "replace")
    except OSError:
        return ""


def _proc_comm(pid):
    try:
        return open(f"/proc/{pid}/comm").read().strip()
    except OSError:
        return ""


def video_last_seq(live_dir):
    try:
        d = json.load(open(os.path.join(live_dir, "live.json")))
        v = next(l for l in d["lanes"].values() if l["kind"] == "video")
        return v.get("last_seq")
    except Exception:                                          # noqa: BLE001
        return None


def ts_mtime(live_dir):
    p = os.path.join(live_dir, "_tv", "live_tv2.ts")
    try:
        return os.path.getmtime(p)
    except OSError:
        return None


def viewer_pids():
    """PIDs of python processes running atsc3_tv.py -- never this
    supervisor itself. The 'atsc3_tv.py' literal cannot match our own
    'atsc3_tvwatch.py' cmdline (no '.py' straight after '_tv'), and the
    /proc walk still excludes os.getpid() explicitly (the pgrep
    self-match trap, bitten three times fleet-wide)."""
    if IS_WIN:
        r = subprocess.run(
            ["powershell", "-NoProfile", "-Command",
             "Get-CimInstance Win32_Process -Filter \"Name='python.exe'\" | "
             r"Where-Object { $_.CommandLine -match 'atsc3_tv\.py' } | "
             "Select-Object -ExpandProperty ProcessId"],
            capture_output=True, text=True)
        return [int(x) for x in r.stdout.split() if x.strip().isdigit()]
    me = os.getpid()
    out = []
    for e in os.listdir("/proc"):
        if not e.isdigit() or int(e) == me:
            continue
        cmd = _proc_cmdline(int(e))
        if "atsc3_tv.py" in cmd and "python" in _proc_comm(int(e)).lower():
            out.append(int(e))
    return out


def _kill_pid_tree(pid):
    if IS_WIN:
        subprocess.run(["taskkill", "/PID", str(pid), "/T", "/F"],
                       capture_output=True)
    else:
        import signal
        try:
            os.kill(pid, signal.SIGKILL)
        except (OSError, ProcessLookupError):
            pass


def kill_all_vlc():
    """The restart path kills EVERY vlc on the box (a wedged viewer's
    player is unrecoverable through player.pid alone). Name-EXACT comm
    match, PID census from /proc -- never a substring pgrep."""
    if IS_WIN:
        subprocess.run(["powershell", "-NoProfile", "-Command",
                        "Get-Process vlc -ErrorAction SilentlyContinue | "
                        "Stop-Process -Force"], capture_output=True)
        return
    for e in os.listdir("/proc"):
        if e.isdigit() and _proc_comm(int(e)) == "vlc":
            _kill_pid_tree(int(e))


def restart_viewer(live_dir):
    for pid in viewer_pids():
        _kill_pid_tree(pid)
    kill_all_vlc()
    time.sleep(3)
    log = open(os.path.join(ROOT, live_dir, "tv_watch_child.log"),
               "ab", buffering=0)
    # E45: this said "atsc3_tv[.]py" -- a regex escape that leaked into the
    # PATH, so a fired watchdog killed the viewer and then FAILED to spawn
    # its replacement (FileNotFoundError). Never bitten only because the
    # wedge never triggered while it ran.
    subprocess.Popen(
        [PY, "-u", os.path.join(ROOT, "tools", "atsc3_tv.py"),
         "--live-dir", live_dir], cwd=ROOT, stdout=log,
        stderr=subprocess.STDOUT)


def main():
    live_dir = sys.argv[1] if len(sys.argv) > 1 else "data/e30"
    stale_need = 120.0           # s the TS may lag while lanes advance --
                                 # must EXCEED the viewer's --audio-hold
                                 # (75 s, E45b): a chunk legitimately held
                                 # for its audio is not a wedge
    last_seq_prev = None
    ts_prev = None
    stale_since = None
    print(f"tv-watchdog on {live_dir}: restart viewer if TS wedges while "
          f"lanes advance", flush=True)
    while True:
        seq = video_last_seq(live_dir)
        ts = ts_mtime(live_dir)
        lanes_advancing = (seq is not None and last_seq_prev is not None
                           and seq != last_seq_prev)
        ts_growing = (ts is not None and ts_prev is not None and ts != ts_prev)
        now = time.time()
        if lanes_advancing and not ts_growing:
            stale_since = stale_since or now
            if now - stale_since >= stale_need:
                print(f"{time.strftime('%H:%M:%S')} WEDGE: lanes advancing "
                      f"(seq {seq}) but TS stale {now - stale_since:.0f}s -- "
                      f"restarting viewer", flush=True)
                restart_viewer(live_dir)
                stale_since = None
                time.sleep(30)
        else:
            stale_since = None
        last_seq_prev, ts_prev = seq, ts
        time.sleep(20)


if __name__ == "__main__":
    sys.exit(main())
