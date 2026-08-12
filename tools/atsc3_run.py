#!/usr/bin/env python3
"""atsc3_run.py -- one command to watch live NextGen TV, with recovery.

This is stvt_run.sh's supervisor, ported. Every rule below was paid for once
already on the ATSC 1.0 side; none of it is new invention.

WHAT IT SUPERVISES
    the chain     `atsc3 watch --live-dir DIR` -- restarted if it dies OR
                  WEDGES (alive, burning a core, emitting nothing). On 8/07 a
                  livelocked chain held 98 % of a core for two minutes with
                  frames, FEC and queues all frozen, and nothing was watching.

THE HEALTH SIGNAL HAS TO BE CHEAP AND HONEST
    "the process exists" is not liveness -- the wedge proved that, and it is
    the same lesson as the zombie PID that held the warden. STVT watches
    live.ts GROW; we watch DIR/live.json's byte count, written by LiveWriter
    from inside the emit path. A wedged thread cannot forge it.

THE RULES THAT STOP THE CURE BEING WORSE
    * STRIKES     one bad sample is a relock, not a wedge. Require N
                  consecutive before even calling it a candidate.
    * GRACE       the chain self-heals most droughts internally (re-acquire).
                  stvt_run waits DROUGHT_GRACE_LOOPS before acting; so do we.
                  Restarting during self-heal converts a 15 s blip into a
                  30 s outage plus a fresh acquisition.
    * BOUNDED     MAX_RESTARTS with a COOLDOWN. A genuinely dead SDR must not
                  produce an infinite respawn storm.
    * ONE OF ME   two supervisors see each other's restart gap as a fresh
                  drought and cascade into a restart storm. Observed on the
                  STVT side contaminating a whole stress run. Single-instance
                  lock, always.
    * ADOPT       a healthy chain that is already running is left alone.

Usage:
    python tools/atsc3_run.py --rf 33
    python tools/atsc3_run.py --rf 33 --cpu-isolate 0-7 --secs 3600
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def log(msg):
    line = f"{time.strftime('%H:%M:%S')} {msg}"
    print(line, flush=True)


class SingleInstance:
    """Refuse to start if another supervisor holds the lock."""

    def __init__(self, path):
        self.path = path
        self.fh = None

    def __enter__(self):
        if os.path.exists(self.path):
            try:
                d = json.load(open(self.path))
                if _pid_alive(d.get("pid")):
                    raise SystemExit(
                        f"another supervisor is running (pid {d.get('pid')}); "
                        f"stop it first, or delete {self.path} if it is stale")
            except (OSError, ValueError):
                pass                      # unreadable = stale, take it
        os.makedirs(os.path.dirname(self.path), exist_ok=True)
        json.dump({"pid": os.getpid(), "started": time.time()},
                  open(self.path, "w"))
        return self

    def __exit__(self, *exc):
        try:
            os.remove(self.path)
        except OSError:
            pass


RADIO_YIELD_PAT = ("radio held by", "not seizing", "REFUSING: inside")


FADE_FEC = 20.0            # 1-min FEC mean under this = the air is gone
FADE_MAX_S = 1800.0        # patience bound: probe once per 30 min anyway

# E73: "the air is gone" and "we are overloaded" look IDENTICAL from
# inside -- both are 0% FEC with frames still advancing, and the fade
# branch would have waited out the second one all night. An independent
# tuner on the same antenna tells them apart, so when it says the air is
# fine we stop waiting and change our own front end instead.
HDHR_BIN = r"C:\Program Files\Silicondust\HDHomeRun\hdhomerun_config.exe"
HDHR_ID = os.environ.get("ATSC3_HDHR_ID", "")  # no default; discover if unset
HDHR_TUNER = os.environ.get("ATSC3_HDHR_TUNER", "1")
HDHR_SNQ_OK = 80
GAIN_LADDER = [4, 6, 2, 8]     # tried in order from the configured value
GAIN_MAX_STEPS = 3             # per GAIN_WIN -- never a gain-thrash storm
GAIN_WIN = 3600.0


def hdhr_snq():
    """Independent tuner's quality on the same antenna, or None."""
    if not os.path.exists(HDHR_BIN):
        return None
    try:
        r = subprocess.run([HDHR_BIN, HDHR_ID, "get",
                            f"/tuner{HDHR_TUNER}/status"],
                           capture_output=True, text=True, timeout=10)
    except (OSError, subprocess.SubprocessError):
        return None
    m = re.search(r"lock=(\S+).*?snq=(\d+)", r.stdout or "")
    if not m or m.group(1) == "none":
        return None
    return int(m.group(2))


def swap_rfgain(extra, gain):
    """Return `extra` with --rfgain set to `gain`."""
    if extra and "--rfgain" in extra:
        return re.sub(r"--rfgain\s+\d+", f"--rfgain {gain}", extra)
    return (extra or "") + f" --rfgain {gain}"


def current_rfgain(extra, default=4):
    m = re.search(r"--rfgain\s+(\d+)", extra or "")
    return int(m.group(1)) if m else default


def chain_status(live_dir, since_off=0, n=12):
    """(frames, fec_mean) from the chain's own status lines, or (None, None).

    THE FADE-IS-NOT-A-WEDGE SIGNAL (E63). read_health() watches live.json,
    which only advances when MPUs are actually emitted -- i.e. when FEC
    SUCCEEDS. So every fade reads as a stall and the supervisor restarts a
    perfectly healthy chain on physics: 41 restarts overnight 8/09, each
    costing a ~90 s re-acquisition, which lengthens the very outage it is
    reacting to. The chain's status line carries a liveness signal that
    survives a dead band -- the Frames counter -- plus the FEC percentage
    that says WHY nothing is being emitted. Together they separate the two
    cases the aggregate cannot:
        frames advancing + FEC ~0   -> the air is gone. Wait.
        frames advancing + FEC high -> the chain runs but assets are dead
                                       (the 8/07 flipped-bit case). Restart.
        frames frozen               -> livelocked. Restart.
    """
    p = os.path.join(live_dir, "chain.log")
    try:
        with open(p, "rb") as f:
            f.seek(since_off)
            txt = f.read(400_000).decode("utf-8", "replace")
    except OSError:
        return None, None
    rows = re.findall(r"([\d.]+)s\s+(\d+) Frames.*?\[\s*([\d.]+)% now\]", txt)
    if not rows:
        return None, None
    tail = rows[-n:]
    frames = int(tail[-1][1])
    fec = sum(float(r[2]) for r in tail) / len(tail)
    return frames, fec


def chain_log_size(live_dir):
    try:
        return os.path.getsize(os.path.join(live_dir, "chain.log"))
    except OSError:
        return 0


def radio_yield_reason(live_dir, rc, since_off=0):
    """Why the chain declined the radio, or None if it truly failed.

    Two shapes, both meaning 'someone else owns the antenna right now':
      * rc=3  -- a scheduled reservation refused up front (meteor window);
      * rc=1  -- radio_lock refused mid-open and the RuntimeError reached
                 the log ('radio held by sonde_rx (priority 100); not
                 seizing'). Exit code ALONE cannot tell this from a real
                 crash, so the log is consulted -- but ONLY to downgrade an
                 rc=1, never to excuse any other failure.

    `since_off` is the chain.log size when THIS instance was spawned, so
    only bytes this run actually wrote are read. Without that bound a
    stale refusal from a previous yield would sit in the tail forever and
    silently re-label every later crash as a yield -- a supervisor that
    stops counting real failures is worse than one that miscounts them.
    """
    if rc == 3:
        return "rc=3, scheduled reservation"
    if rc != 1:
        return None
    p = os.path.join(live_dir, "chain.log")
    try:
        with open(p, "rb") as f:
            f.seek(since_off)
            txt = f.read(200_000).decode("utf-8", "replace")
    except OSError:
        return None
    for pat in RADIO_YIELD_PAT:
        if pat in txt:
            line = next((l.strip() for l in reversed(txt.splitlines())
                         if pat in l), pat)
            return f"rc=1, {line[:90]}"
    return None


def child_snapshot(proc):
    """Live children of the chain, for post-mortem reaping (E61)."""
    try:
        import psutil
        return psutil.Process(proc.pid).children(recursive=True)
    except Exception:                                          # noqa: BLE001
        return []


def stop_tree(proc, snapshot=None):
    """Stop the chain AND its whole process tree.

    E61 (8/10, measured): proc.terminate() kills only the parent `-m atsc3
    watch`; its two multiprocessing decode workers survive as orphans.
    One fade night of supervisor restarts leaked 42 zombie decode
    processes (84 threads) that dragged the live chain to 0.41x real time
    -- the leak was the load. Children are snapshotted BEFORE the parent
    dies (afterwards they re-parent and cannot be found by tree walk),
    the tree is killed, and any snapshot survivor is reaped by exact
    psutil identity (pid+create_time, immune to pid reuse). `snapshot` is
    the caller's LAST-KNOWN-ALIVE child list -- required when the chain
    died on its own (a dead pid cannot be walked)."""
    kids = list(snapshot or [])
    try:
        import psutil
        kids = psutil.Process(proc.pid).children(recursive=True)
    except Exception:                                          # noqa: BLE001
        pass
    try:
        if proc.poll() is None:
            if os.name == "nt":
                subprocess.run(["taskkill", "/PID", str(proc.pid),
                                "/T", "/F"], capture_output=True)
            else:
                proc.terminate()
            try:
                proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                proc.kill()
                try:
                    proc.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    pass
    except OSError:
        pass
    for c in kids:
        try:
            if c.is_running():
                c.kill()
                log(f"  reaped orphaned chain child pid {c.pid}")
        except Exception:                                      # noqa: BLE001
            pass


def _pid_alive(pid):
    if not pid:
        return False
    try:
        import psutil
        return psutil.pid_exists(int(pid))
    except Exception:                                          # noqa: BLE001
        return False


def read_health(live_dir):
    """(bytes, updated, media_s, lanes) from the writer's heartbeat, or None.

    `bytes` is the SUM over every lane, and that is exactly how this
    supervisor was fooled on 8/07: a corrupted MPU sequence number killed
    video and both audio assets permanently, the subtitle lane kept writing
    ~3 kB/s, the total kept climbing, and the supervisor reported health for
    twenty minutes over a picture that had stopped.  `lanes` is returned so
    the caller can judge each asset on its OWN growth -- an aggregate is not
    health when any one member can carry it.
    """
    p = os.path.join(live_dir, "live.json")
    try:
        d = json.load(open(p))
        lanes = {k: v.get("bytes", 0) for k, v in (d.get("lanes") or {}).items()}
        return (d.get("bytes", 0), d.get("updated", 0.0),
                d.get("media_s", 0.0), lanes)
    except (OSError, ValueError):
        return None


def find_python(explicit=None):
    """An interpreter that can actually open the radio.

    The chain needs SoapySDR, which on this fleet lives in radioconda and NOT
    in the default python. Spawning `sys.executable` blindly means the
    supervisor faithfully restarts a chain that cannot possibly work, six
    times, and burns the radio window doing it -- which is exactly what
    happened on the first live soak. Check before launching, not after.
    """
    cands = [explicit] if explicit else []
    cands += [sys.executable,
              os.path.expanduser(r"~\radioconda\python.exe"),
              os.path.expanduser("~/radioconda/python.exe"),
              os.path.expanduser("~/radioconda/bin/python")]
    for c in cands:
        if not c or not (os.path.exists(c) or c == sys.executable):
            continue
        try:
            r = subprocess.run([c, "-c", "import SoapySDR"],
                               capture_output=True, timeout=60)
            if r.returncode == 0:
                return c
        except (OSError, subprocess.SubprocessError):
            continue
    return None


def spawn(a):
    cmd = [a.python, "-m", "atsc3", "watch",
           "--rf", str(a.rf), "--live-dir", a.live_dir,
           "--player", "none", "--json", os.path.join(a.live_dir, "run.json")]
    if a.cpu_isolate:
        cmd += ["--cpu-isolate", a.cpu_isolate]
    if a.rate:
        cmd += ["--rate", str(a.rate)]
    if a.ant:
        cmd += ["--ant", a.ant]
    if a.extra:
        cmd += a.extra.split()
    logf = open(os.path.join(a.live_dir, "chain.log"), "ab", buffering=0)
    log(f"  starting chain: {' '.join(cmd[2:])}")
    return subprocess.Popen(cmd, cwd=ROOT, stdout=logf, stderr=subprocess.STDOUT)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rf", type=int, required=True,
                    help="RF channel number (required: a shipped default is wrong for every reader)")
    ap.add_argument("--python", default=None,
                    help="interpreter for the chain (must import SoapySDR); "
                         "auto-detected if omitted")
    ap.add_argument("--live-dir", default=None)
    ap.add_argument("--cpu-isolate", default=None,
                    help="core list for the chain, e.g. 0-7")
    ap.add_argument("--rate", type=float, default=None)
    ap.add_argument("--ant", default=None)
    ap.add_argument("--extra", default=None, help="extra args for `watch`")
    ap.add_argument("--secs", type=float, default=0.0,
                    help="stop after N seconds (0 = until Ctrl-C)")
    ap.add_argument("--check", type=float, default=5.0,
                    help="seconds between health samples")
    ap.add_argument("--strikes", type=int, default=3,
                    help="consecutive stalled samples before it is a candidate")
    ap.add_argument("--grace", type=int, default=4,
                    help="further stalled checks to allow in-chain self-heal")
    ap.add_argument("--max-restarts", type=int, default=30)
    ap.add_argument("--cooldown", type=float, default=5.0)
    ap.add_argument("--startup", type=float, default=45.0,
                    help="grace after a (re)start before health counts -- "
                         "acquisition legitimately emits nothing for a while")
    a = ap.parse_args()
    a.live_dir = a.live_dir or os.path.join(ROOT, "data", "atsc3_live")
    replay = bool(a.extra and "--capture" in a.extra)
    a.python = find_python(a.python) or (sys.executable if replay else None)
    if a.python is None:
        print("no interpreter here can import SoapySDR -- the chain could "
              "never open the radio."
              "\n  Pass --python <path>, or use --extra '--capture FILE' "
              "to replay without a radio.", file=sys.stderr)
        return 2
    log(f"  chain interpreter: {a.python}")
    os.makedirs(a.live_dir, exist_ok=True)

    lock = os.path.join(a.live_dir, "supervisor.lock")
    deadline = time.time() + a.secs if a.secs else None
    restarts = 0

    with SingleInstance(lock):
        log_off = chain_log_size(a.live_dir)
        proc = spawn(a)
        started = time.time()
        last_bytes, strikes, grace_used = -1, 0, 0
        last_lanes = {}
        kids = []                 # last-known-alive children (E61)
        last_frames = None        # chain liveness across a fade (E63)
        fade_since = None
        fade_logged = 0.0
        gs = []                   # gain-step timestamps (E73)
        last_gain = 0.0
        log(f"  supervising -- live dir {a.live_dir}")
        try:
            while True:
                if deadline and time.time() > deadline:
                    log("  time limit reached")
                    break
                time.sleep(a.check)

                died = proc.poll() is not None
                if not died:
                    kids = child_snapshot(proc) or kids
                yield_why = radio_yield_reason(
                    a.live_dir, proc.returncode, log_off) if died else None
                if yield_why:
                    # The chain DECLINED the radio -- a scheduled reservation
                    # (rc=3) or a higher-priority holder (rc=1 + the lock's
                    # refusal on stderr). Yielding is the fleet's single-
                    # tenant law working, never a fault, so it must not spend
                    # the restart budget: measured 8/10, 28 sonde refusals
                    # burned 28 of 200 restarts at a 160 s cadence, and the
                    # meteor window projected 315 spawns against the same
                    # budget. Stand by patiently, spend nothing.
                    log(f"  chain YIELDED the radio ({yield_why}) -- "
                        f"standing by, not counting a restart")
                    time.sleep(max(a.cooldown, 120.0))
                    log_off = chain_log_size(a.live_dir)
                    proc = spawn(a)
                    started = time.time()
                    last_bytes, strikes, grace_used = -1, 0, 0
                    last_lanes = {}
                    kids = []
                    last_frames, fade_since = None, None
                    continue
                if died:
                    log(f"  chain EXITED (rc={proc.returncode})")
                elif time.time() - started < a.startup:
                    continue                      # still acquiring; not a wedge
                else:
                    h = read_health(a.live_dir)
                    if h is None:
                        strikes += 1
                        log(f"  no heartbeat yet ({strikes}/{a.strikes})")
                    else:
                        b, upd, media, lanes = h
                        stale = time.time() - upd
                        # A lane that has EVER grown and then goes quiet while
                        # a sibling keeps writing is a dead asset, not a dead
                        # signal -- if the RF were gone every lane would stop
                        # together.  This is the case the aggregate cannot see.
                        dead = [k for k, v in lanes.items()
                                if v == last_lanes.get(k, -1)]
                        alive = len(lanes) - len(dead)
                        if dead and alive:
                            strikes += 1
                            log(f"  LANE(S) DEAD: {', '.join(sorted(dead))} "
                                f"unchanged while {alive} lane(s) still write "
                                f"({strikes}/{a.strikes})")
                        elif b > last_bytes:
                            if strikes or grace_used:
                                log(f"  recovered -- {media:.1f} s of media, "
                                    f"{b/1e6:.1f} MB written")
                            last_bytes, strikes, grace_used = b, 0, 0
                            last_lanes = lanes
                            continue
                        else:
                            strikes += 1
                            log(f"  STALLED: {b/1e6:.1f} MB unchanged, heartbeat "
                                f"{stale:.0f} s old ({strikes}/{a.strikes})")
                        last_lanes = lanes

                    if strikes < a.strikes:
                        continue

                    # BEFORE calling it a wedge: is the air simply gone?
                    # (E63) A fade stops every lane at once and the chain
                    # keeps processing frames at ~0% FEC. Restarting that
                    # cures nothing and costs a fresh acquisition, so wait
                    # -- but not forever: probe once per FADE_MAX_S in case
                    # our own acquisition, not the sky, is what is stuck.
                    fr, fec = chain_status(a.live_dir, log_off)
                    if fr is not None and fec is not None and fec < FADE_FEC:
                        alive_now = fr != last_frames
                        last_frames = fr
                        if alive_now:
                            fade_since = fade_since or time.time()
                            waited = time.time() - fade_since
                            # E73: before settling in to wait, ask the
                            # independent tuner whether there is anything
                            # to wait FOR. If the air is fine and we alone
                            # cannot decode it, waiting is the wrong move
                            # and the front end is the thing to change.
                            snq = (hdhr_snq() if waited >= 240.0
                                   and time.time() - last_gain > 600.0
                                   else None)
                            if snq is not None and snq >= HDHR_SNQ_OK:
                                gs[:] = [t for t in gs
                                         if time.time() - t < GAIN_WIN]
                                if len(gs) < GAIN_MAX_STEPS:
                                    cur = current_rfgain(a.extra)
                                    nxt = next((g for g in GAIN_LADDER
                                                if g != cur), GAIN_LADDER[0])
                                    order = [g for g in GAIN_LADDER if g != cur]
                                    nxt = order[len(gs) % len(order)]
                                    log(f"  NOT A FADE: independent tuner "
                                        f"reads snq={snq} on the same "
                                        f"antenna while we sit at "
                                        f"{fec:.0f}% -- our front end is "
                                        f"the problem. rfgain_sel "
                                        f"{cur} -> {nxt} "
                                        f"(step {len(gs)+1}/{GAIN_MAX_STEPS})")
                                    a.extra = swap_rfgain(a.extra, nxt)
                                    gs.append(time.time())
                                    last_gain = time.time()
                                    fade_since = None
                                    stop_tree(proc, kids)
                                    kids = []
                                    time.sleep(a.cooldown)
                                    log_off = chain_log_size(a.live_dir)
                                    proc = spawn(a)
                                    started = time.time()
                                    last_bytes, strikes, grace_used = -1, 0, 0
                                    last_lanes = {}
                                    last_frames = None
                                    continue
                                log(f"  gain-step budget spent "
                                    f"({GAIN_MAX_STEPS}/"
                                    f"{GAIN_WIN/60:.0f}min) -- waiting "
                                    f"instead; the air may be fine but I "
                                    f"will not thrash the front end")
                            if waited < FADE_MAX_S:
                                if time.time() - fade_logged > 300:
                                    fade_logged = time.time()
                                    log(f"  BAND DOWN (FEC {fec:.0f}% over "
                                        f"the last minute) but the chain is "
                                        f"processing frames -- a fade, not a "
                                        f"wedge. Waiting "
                                        f"({waited/60:.0f}/"
                                        f"{FADE_MAX_S/60:.0f} min)")
                                strikes = 0
                                grace_used = 0
                                continue
                            log(f"  fade has lasted {waited/60:.0f} min -- "
                                f"probing with one restart in case the "
                                f"acquisition is stuck, not the sky")
                            fade_since = None
                    else:
                        fade_since = None

                    # A candidate wedge. Let the chain try to fix itself
                    # first -- it re-acquires on its own and usually wins.
                    if grace_used < a.grace:
                        grace_used += 1
                        log(f"  candidate wedge -- allowing in-chain "
                            f"self-heal ({grace_used}/{a.grace})")
                        continue
                    log("  WEDGED and did not self-heal -- restarting")

                restarts += 1
                if restarts > a.max_restarts:
                    log(f"  MAX_RESTARTS ({a.max_restarts}) hit -- giving up. "
                        f"Check the SDR.")
                    return 1
                stop_tree(proc, kids)         # E61: the WHOLE tree, always
                kids = []
                time.sleep(a.cooldown)
                log_off = chain_log_size(a.live_dir)
                proc = spawn(a)
                started = time.time()
                last_bytes, strikes, grace_used = -1, 0, 0
                last_lanes = {}
                last_frames, fade_since = None, None
                log(f"  restart {restarts}/{a.max_restarts}")
        except KeyboardInterrupt:
            log("  interrupted")
        finally:
            stop_tree(proc, kids)             # E61: the WHOLE tree, always
            h = read_health(a.live_dir)
            if h:
                log(f"  final: {h[0]/1e6:.1f} MB, {h[2]:.1f} s of media, "
                    f"{restarts} restart(s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
