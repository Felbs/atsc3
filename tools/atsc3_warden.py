#!/usr/bin/env python3
"""atsc3_warden.py -- ONE supervisor owning the whole live-TV stack (E61).

The weekend's failure catalogue (LIVE_TV_LAB E35-E60 + the 8/10 night) proved
that four separately-supervised components fail in the SEAMS: the chain
supervisor watched the chain, tvwatch watched the viewer, the viewer watched
VLC -- and a dead audio worker, a vanished player, a leaked process tree and
a duplicate spawner each fell between two watchers. This is the missing
layer: one loop, one census, every role.

ROLES (all matched by EXACT cmdline census, self-excluded -- the pgrep
self-match trap has bitten this fleet five times):
    run     tools/atsc3_run.py          (chain supervisor; it owns the SDR
                                         and the chain's own restarts)
    eng     tools/atsc3_audio.py --pid 13
    spa     tools/atsc3_audio.py --pid 14
    subs    tools/atsc3_subs.py
    viewer  tools/atsc3_tv.py           (spawns and owns its VLC)
    player  vlc.exe playing OUR sink    (owned by the viewer; the warden
                                         only verifies and unwedges)

Per role: LIVENESS (exactly one process) + PROGRESS (its output moves when
its input moves). Every progress rule is CONDITIONED ON THE UPSTREAM so a
band fade is never called a failure (fades are physics -- RF33 dies nightly
19:56-22:10+): a worker is only wedged if its LANE grew while its output
did not; the viewer only if the VIDEO LANE grew while the TS did not.

Restarts are cooldown-spaced and stampede-gated (the TURBO law): more than
STAMPEDE_N restarts of one role in STAMPEDE_WIN seconds stands that role
down (logged loudly, everything else keeps running) until the window
drains. Every action is logged with its reason to _warden/warden.log.

The warden itself is crash-only: no threads, no daemon magic, state
re-derived from disk every loop; tools/atsc3_warden.bat relaunches it if it
dies, and a single-instance lock keeps relaunches honest. Drop a file at
_warden/pause to make it observe-but-not-act (gates use this).

Usage:
    python tools/atsc3_warden.py --live-dir data/e31 --rf 33
"""
from __future__ import annotations

import argparse
import base64
import json
import os
import re
import subprocess
import sys
import time
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
IS_WIN = (os.name == "nt")

LOOP_S = 15.0
PAUSE_MAX_S = 900.0        # a maintenance pause older than this is stale
STAMPEDE_N = 4                    # restarts ...
STAMPEDE_WIN = 3600.0             # ... per this window = stand down
DUP_CONFIRM = 2                   # consecutive loops before killing extras
# E73 (8/10 20:41, MEASURED live while the chain sat at 0% on air the
# HDHomeRun was decoding at snq=100 on the same antenna):
#     rfgain_sel=2 (the old default)   0.0% FEC, 0.39x
#     rfgain_sel=4                   100.0% FEC, 1.02x
#     rfgain_sel=6                   100.0% FEC, 1.03x
# Two steps less LNA gain is the difference between a dead chain and a
# perfect one. 4 is chosen over 6 because it also measured fine this
# morning (97.3% vs 96.3% at gain 2), so it is the setting with evidence
# in BOTH signal conditions.
CHAIN_EXTRA = ("--assets all --raw-queue 600 --accel cpu "
               "--decode-procs 2 --threads 4 --fe-threads 6 --rfgain 4")
TS_ROLL_BYTES = 30_000_000_000    # maintenance viewer restart (sink reset)

try:
    import psutil
except ImportError:               # the warden cannot run blind
    print("atsc3_warden needs psutil", file=sys.stderr)
    raise


def now_s():
    return time.strftime("%m/%d %H:%M:%S")


class Warden:
    def __init__(self, a):
        self.d = a.live_dir
        self.a = a
        self.wdir = os.path.join(self.d, "_warden")
        os.makedirs(self.wdir, exist_ok=True)
        self.logf = open(os.path.join(self.wdir, "warden.log"), "a",
                         buffering=1, encoding="utf-8")
        self.state_p = os.path.join(self.wdir, "warden_state.json")
        try:
            self.state = json.load(open(self.state_p))
        except (OSError, ValueError):
            self.state = {}
        self.state.setdefault("restarts", {})    # role -> [wall times]
        # in-memory only (crash resets them -- that is fine, they re-derive)
        self.dup_seen = {}          # role -> consecutive dup-loop count
        self.obs = {}               # role -> last observation dict
        self.http_fails = 0
        self.standdown_logged = {}
        self._stale_pause_logged = 0.0

    # ------------------------------------------------------------- infra
    def log(self, msg):
        line = f"{now_s()} {msg}"
        print(line, flush=True)
        self.logf.write(line + "\n")

    def save_state(self):
        try:
            json.dump(self.state, open(self.state_p, "w"))
        except OSError:
            pass

    def paused(self):
        """Is supervision suspended for maintenance?

        THE MARKER EXPIRES (E68). A pause file is a promise by whoever
        wrote it to come back and remove it -- and unattended overnight,
        a maintenance script that dies mid-capture would otherwise leave
        the whole stack unsupervised until someone noticed in the
        morning. Crash-only means the SAFE state must be the one that
        happens when nobody is left to act: after PAUSE_MAX_S the marker
        is ignored, loudly, and supervision resumes on its own."""
        p = os.path.join(self.wdir, "pause")
        try:
            age = time.time() - os.path.getmtime(p)
        except OSError:
            return False
        if age <= PAUSE_MAX_S:
            return True
        if time.time() - self._stale_pause_logged > 300:
            self._stale_pause_logged = time.time()
            self.log(f"STALE PAUSE marker ({age/60:.0f} min old > "
                     f"{PAUSE_MAX_S/60:.0f}) -- IGNORING it and resuming "
                     f"supervision; whoever paused me did not come back")
        return False

    # ------------------------------------------------------------ census
    def census(self):
        """[{pid, name, cmd, t0}] for python + vlc, self-excluded."""
        me = os.getpid()
        out = []
        for p in psutil.process_iter(["pid", "name", "cmdline",
                                      "create_time"]):
            try:
                nm = (p.info["name"] or "").lower()
                if nm not in ("python.exe", "python", "python3", "vlc.exe",
                              "vlc"):
                    continue
                if p.info["pid"] == me:
                    continue
                cmd = " ".join(p.info["cmdline"] or [])
                out.append({"pid": p.info["pid"], "name": nm, "cmd": cmd,
                            "t0": p.info["create_time"], "proc": p})
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                continue
        return out

    def _ours(self, cmd):
        """Does this cmdline's --live-dir point at OUR live dir? A warden
        must never touch another live dir's stack -- nor the gates'
        scratch stacks (measured 8/10: gate_e61's sim viewer matched the
        bare pattern and would have been killed as a 'duplicate')."""
        m = re.search(r"--live-dir\s+(\S+)", cmd)
        if not m:
            return False
        p = m.group(1).strip('"')
        if not os.path.isabs(p):
            p = os.path.join(ROOT, p)
        try:
            return os.path.normcase(os.path.normpath(p)) == \
                os.path.normcase(os.path.normpath(self.d))
        except (OSError, ValueError):
            return False

    def classify(self, census):
        roles = {k: [] for k in
                 ("run", "eng", "spa", "subs", "viewer", "tvwatch", "vlc")}
        sink = os.path.normcase(os.path.normpath(
            os.path.join(self.d, "_tv", "live_tv2.ts")))
        for e in census:
            c = e["cmd"]
            if e["name"].startswith("vlc"):
                # OUR sink only -- by full path, never any other player
                # (gate scratch sinks are also named live_tv2.ts)
                if sink in os.path.normcase(c):
                    roles["vlc"].append(e)
                continue
            if "atsc3_warden.py" in c:
                continue                        # another warden: lock's job
            if "atsc3_tvwatch.py" in c:
                # positional live-dir; deprecated everywhere -- always ours
                roles["tvwatch"].append(e)
                continue
            if not self._ours(c):
                continue
            if "atsc3_run.py" in c:
                roles["run"].append(e)
            elif "atsc3_audio.py" in c:
                roles["spa" if "--pid 14" in c else "eng"].append(e)
            elif "atsc3_subs.py" in c:
                roles["subs"].append(e)
            elif re.search(r"atsc3_tv\.py", c):
                roles["viewer"].append(e)
        return roles

    # ------------------------------------------------------------- kills
    def kill_tree(self, pid, why):
        self.log(f"KILL pid {pid} (tree): {why}")
        if IS_WIN:
            subprocess.run(["taskkill", "/PID", str(pid), "/T", "/F"],
                           capture_output=True)
        else:
            try:
                p = psutil.Process(pid)
                for c in p.children(recursive=True):
                    try:
                        c.kill()
                    except psutil.Error:
                        pass
                p.kill()
            except psutil.Error:
                pass

    def sweep_orphan_forks(self, census):
        """multiprocessing spawn-children whose parent is DEAD are garbage
        by construction (their pipe is gone) -- the 8/10 night leaked 42 of
        them and they dragged the live chain to 0.41x real time."""
        alive = {e["pid"] for e in census}
        alive |= {p.pid for p in psutil.process_iter()}
        for e in census:
            if "multiprocessing" not in e["cmd"] \
                    or "parent_pid=" not in e["cmd"]:
                continue
            m = re.search(r"parent_pid=(\d+)", e["cmd"])
            if m and int(m.group(1)) not in alive:
                self.kill_tree(e["pid"],
                               f"orphaned fork (parent {m.group(1)} dead)")

    # ------------------------------------------------------------ spawns
    def _env(self, blas):
        env = dict(os.environ)
        for k in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                  "OMP_NUM_THREADS"):
            env[k] = str(blas)
        return env

    def _spawn(self, role, cmd, blas, logname):
        lf = open(os.path.join(self.wdir, logname), "ab", buffering=0)
        flags = 0
        if IS_WIN:
            flags = (subprocess.CREATE_NEW_PROCESS_GROUP
                     | subprocess.CREATE_NO_WINDOW)
        p = subprocess.Popen(cmd, cwd=ROOT, stdout=lf,
                             stderr=subprocess.STDOUT, env=self._env(blas),
                             creationflags=flags,
                             start_new_session=not IS_WIN)
        self.log(f"SPAWN {role}: pid {p.pid}: {' '.join(cmd[1:])}")
        return p

    def spawn_role(self, role):
        py = sys.executable
        d = self.d
        if role == "run":
            self._spawn(role, [py, "-u", os.path.join(HERE, "atsc3_run.py"),
                               "--rf", str(self.a.rf), "--secs", "0",
                               "--max-restarts", "200", "--cooldown", "20",
                               "--check", "20", "--startup", "90",
                               "--live-dir", d, "--extra", CHAIN_EXTRA],
                        blas=4, logname="run_w.log")
        elif role == "eng":
            self._spawn(role, [py, "-u", os.path.join(HERE, "atsc3_audio.py"),
                               "--live-dir", d, "--pid", "13",
                               "--interval", "8", "--wait", "900",
                               "--start-behind", "150"],
                        blas=1, logname="eng_w.log")
        elif role == "spa":
            self._spawn(role, [py, "-u", os.path.join(HERE, "atsc3_audio.py"),
                               "--live-dir", d, "--pid", "14",
                               "--element", "pair",
                               "--out", os.path.join(d, "live_audio_spa.wav"),
                               "--interval", "8", "--wait", "900",
                               "--start-behind", "150"],
                        blas=1, logname="spa_w.log")
        elif role == "subs":
            self._spawn(role, [py, "-u", os.path.join(HERE, "atsc3_subs.py"),
                               "--live-dir", d, "--interval", "20",
                               "--wait", "900"],
                        blas=1, logname="subs_w.log")
        elif role == "viewer":
            self._spawn(role, [py, "-u", os.path.join(HERE, "atsc3_tv.py"),
                               "--live-dir", d,
                               "--lag", str(self.a.lag),
                               "--lead", str(self.a.lead)],
                        blas=1, logname="viewer_w.log")

    # --------------------------------------------------------- stampede
    def may_restart(self, role):
        rs = self.state["restarts"].setdefault(role, [])
        now = time.time()
        rs[:] = [t for t in rs if now - t < STAMPEDE_WIN]
        if len(rs) >= STAMPEDE_N:
            last = self.standdown_logged.get(role, 0)
            if now - last > 300:
                self.log(f"STAND-DOWN {role}: {len(rs)} restarts in "
                         f"{STAMPEDE_WIN/60:.0f} min -- not restarting until "
                         f"the window drains (rest of stack unaffected)")
                self.standdown_logged[role] = now
            return False
        return True

    def note_restart(self, role):
        self.state["restarts"].setdefault(role, []).append(time.time())
        self.save_state()

    # ---------------------------------------------------------- helpers
    def lanes(self):
        try:
            return json.load(open(os.path.join(self.d, "live.json"))
                             ).get("lanes", {})
        except (OSError, ValueError):
            return {}

    def lane_bytes(self, pid=None, kind=None):
        for l in self.lanes().values():
            if pid is not None and l.get("pid") == pid:
                return l.get("bytes", 0)
            if kind is not None and l.get("kind") == kind:
                return l.get("bytes", 0)
        return None

    def wav_samples(self, name):
        try:
            return json.load(open(os.path.join(self.d, name))
                             ).get("wav_samples")
        except (OSError, ValueError):
            return None

    def ts_stat(self):
        p = os.path.join(self.d, "_tv", "live_tv2.ts")
        try:
            st = os.stat(p)
            return st.st_mtime, st.st_size
        except OSError:
            return None, None

    def player_http_ok(self):
        """None = no record; True/False = probe result."""
        try:
            rec = json.load(open(os.path.join(self.d, "_tv", "player.pid")))
            port, pw = rec["port"], rec["pw"]
        except (OSError, ValueError, KeyError):
            return None
        rq = urllib.request.Request(
            f"http://127.0.0.1:{port}/requests/status.json",
            headers={"Authorization": "Basic " + base64.b64encode(
                f":{pw}".encode()).decode()})
        try:
            with urllib.request.urlopen(rq, timeout=2.0):
                return True
        except Exception:                                      # noqa: BLE001
            return False

    # ------------------------------------------------- progress verdicts
    def worker_wedged(self, role, lane_pid, state_name):
        """Alive-but-stuck: lane grew >= 1 MB across >= 600 s while
        wav_samples never moved. Conditioned on the lane so a fade can
        never trip it."""
        lb = self.lane_bytes(pid=lane_pid)
        ws = self.wav_samples(state_name)
        o = self.obs.get(role)
        now = time.time()
        if lb is None or ws is None:
            self.obs[role] = None
            return False
        if o is None or ws != o["ws"]:
            self.obs[role] = {"ws": ws, "lb": lb, "t": now}
            return False
        if now - o["t"] >= 600 and lb - o["lb"] >= 1_000_000:
            self.obs[role] = None
            return True
        return False

    def lane_stalled(self, window=180.0):
        """True when the video lane has not advanced for `window` seconds --
        i.e. the chain is not delivering (fade, or the radio is held by a
        higher-priority service). E62: this is the NON-CIRCULAR health
        signal. It comes from the chain, never from the viewer, so it can
        be used to decide whether the viewer owes us a window."""
        v = None
        for l in self.lanes().values():
            if l.get("kind") == "video":
                v = l.get("last_seq")
        now = time.time()
        o = self.obs.setdefault("lane_adv", {"seq": None, "t": now})
        if v is not None and v != o["seq"]:
            o["seq"], o["t"] = v, now
            return False
        return now - o["t"] > window

    def gen_media_s(self):
        """Seconds of media the CURRENT video generation holds. (E72b)

        After a lane roll the viewer must wait for --lag seconds of lane
        before it can even anchor, then bank --lead before a window opens.
        Any rule that demands output from the viewer has to allow for that
        fill window, or it punishes the viewer for arithmetic."""
        for l in self.lanes().values():
            if l.get("kind") == "video":
                f, s = l.get("first_seq"), l.get("last_seq")
                if f is not None and s is not None:
                    return (s - f) * 2.002
        return 0.0

    def viewer_wedged(self):
        """The tvwatch rule, absorbed: video lane advances, TS does not.
        180 s > the 75 s audio-hold and > the governor's 90 s max_hold.

        LANE-AGE GUARD (measured 02:04:27): after a generation roll the
        viewer structurally cannot mux until the NEWBORN lane is older
        than its --lag (120 s), and the governor then withholds up to a
        60 s cushion -- the fill window is ~3 min of stale-TS-while-lane-
        advances that is HEALTH, not a wedge. The first cut killed a
        viewer 90 s before its first burst. Wedge now requires the
        current generation to hold > lag+lead+30 s of media.

        ADVANCING MEANS *STILL* ADVANCING (E72, measured 18:43:18). The
        first cut compared the lane's seq against the ARM-TIME value, so
        once the lane moved past that point and then STOPPED, `v !=
        o["seq"]` stayed true forever and the 180 s timer kept counting as
        though the lane were still running. sonde_rx took the radio at
        ~18:40, the lane correctly froze, the viewer correctly stopped
        emitting -- and the warden killed a perfectly healthy viewer for
        it, truncating a 5 GB sink to zero in the process.

        Same shape as every other bug this campaign found: a stale
        reference standing in for a live measurement. Advancement is now
        judged BETWEEN CONSECUTIVE SAMPLES, so a lane that stops -- fade,
        radio yielded, chain restarting -- re-arms the timer instead of
        running it down."""
        v = first = last = None
        for l in self.lanes().values():
            if l.get("kind") == "video":
                v = l.get("last_seq")
                first, last = l.get("first_seq"), l.get("last_seq")
        mt, _sz = self.ts_stat()
        o = self.obs.get("viewer")
        now = time.time()
        if v is None or mt is None:
            self.obs["viewer"] = None
            return False
        gen_media = (last - first) * 2.002 \
            if (first is not None and last is not None) else 0.0
        if gen_media <= self.a.lag + self.a.lead + 30.0:
            self.obs["viewer"] = None       # post-roll fill window
            return False
        advancing = (o is not None and v != o.get("prev"))
        if o is None or not advancing or mt != o["mt"]:
            # lane STOPPED (fade / radio yielded), or the TS moved: healthy.
            # Re-arm, carrying the current seq as the next comparison point.
            self.obs["viewer"] = {"prev": v, "mt": mt, "t": now}
            return False
        o["prev"] = v                    # still advancing -- keep the timer
        if now - o["t"] >= 180:
            self.obs["viewer"] = None
            return True
        return False

    # -------------------------------------------------------------- run
    def loop(self):
        census = self.census()
        roles = self.classify(census)

        if self.paused():
            self.obs.clear()
            return

        # 0. legacy double-spawner: tvwatch fights the warden for the
        #    viewer role (measured 8/10 00:06: two viewers, each killing
        #    the other's VLC via the shared player.pid record)
        for e in roles["tvwatch"]:
            self.kill_tree(e["pid"], "tvwatch is superseded by the warden")

        # 1. orphaned fork sweep
        self.sweep_orphan_forks(census)

        # 2. one-each (kill extras, keep newest, 2-loop confirmation)
        for role in ("run", "eng", "spa", "subs", "viewer", "vlc"):
            es = roles[role]
            if len(es) > 1:
                self.dup_seen[role] = self.dup_seen.get(role, 0) + 1
                if self.dup_seen[role] >= DUP_CONFIRM:
                    es_sorted = sorted(es, key=lambda e: e["t0"])
                    for e in es_sorted[:-1]:
                        self.kill_tree(e["pid"],
                                       f"duplicate {role} (keeping newest "
                                       f"pid {es_sorted[-1]['pid']})")
                    self.dup_seen[role] = 0
            else:
                self.dup_seen[role] = 0

        lanes = self.lanes()
        have_audio13 = any(l.get("pid") == 13 for l in lanes.values())
        have_audio14 = any(l.get("pid") == 14 for l in lanes.values())
        have_subs = any(l.get("kind") == "subs" for l in lanes.values())

        # 3. liveness (+ spawn preconditions so a dead band does not spin
        #    workers that would only time out waiting for a lane)
        want = {"run": True, "eng": have_audio13, "spa": have_audio14,
                "subs": have_subs, "viewer": True}
        for role, wanted in want.items():
            if roles[role] or not wanted:
                continue
            if self.may_restart(role):
                self.log(f"MISSING {role} -- spawning")
                self.note_restart(role)
                self.spawn_role(role)

        # 4. progress: workers
        for role, lane_pid, st_name in (
                ("eng", 13, "live_audio.state.json"),
                ("spa", 14, "live_audio_spa.state.json")):
            if roles[role] and self.worker_wedged(role, lane_pid, st_name):
                if self.may_restart(role):
                    self.kill_tree(roles[role][0]["pid"],
                                   f"{role} wedged: lane grew, wav did not "
                                   f"(>=600 s)")
                    self.note_restart(role)
                    self.spawn_role(role)

        # 5. progress: viewer (TS vs video lane)
        if roles["viewer"] and self.viewer_wedged():
            if self.may_restart("viewer"):
                self.kill_tree(roles["viewer"][0]["pid"],
                               "viewer wedged: video lane advanced, TS "
                               "stale 180 s")
                for e in roles["vlc"]:
                    self.kill_tree(e["pid"], "viewer's player (wedge cycle)")
                self.note_restart("viewer")
                self.spawn_role("viewer")

        # 6. player checks (viewer owns VLC; the warden unwedges)
        if roles["viewer"] and roles["vlc"]:
            ok = self.player_http_ok()
            self.http_fails = self.http_fails + 1 if ok is False else 0
            if self.http_fails >= 4:        # ~1 min of dead control surface
                self.http_fails = 0
                if self.may_restart("player"):
                    self.note_restart("player")
                    for e in roles["vlc"]:
                        self.kill_tree(e["pid"],
                                       "vlc alive but http dead 4 probes "
                                       "-- wedged player, viewer respawns")
        elif roles["viewer"] and not roles["vlc"]:
            # no player at all: the viewer's own respawn path owns this;
            # the warden acts only when that path has plainly given up --
            # viewer old, TS flowing, and still no window.
            v = roles["viewer"][0]
            mt, _ = self.ts_stat()
            aged = time.time() - v["t0"] > 300
            # E62 (corrected): the original ts_fresh precondition was not
            # merely vacuous, it was CIRCULAR -- it asked the VIEWER'S OWN
            # OUTPUT whether the viewer was healthy, so a broken viewer
            # silenced the check that would have restarted it (measured
            # 8/10: no player 02:00-07:18). The first correction dropped the
            # precondition entirely and thrashed instead: with the radio held
            # by sonde_rx there is no content, the viewer structurally cannot
            # open a window, and the warden killed it every 10 min (07:26,
            # 07:37, 07:47) for failing to do the impossible.
            # The right signal is the CHAIN's output, which is independent of
            # the viewer: demand a player only while the video LANE is
            # advancing (there is something to show). No signal -> no window
            # is not a fault; lane advancing with no window IS.
            o = self.obs.setdefault("no_vlc", {"t": None})
            # E72b: the SAME post-roll fill window viewer_wedged() already
            # allows for. Measured 20:26:10 -- the lane rolled to gen 19 at
            # 20:21, the viewer re-anchored at the new lane's FIRST slot and
            # was still waiting for --lag of content (a brief fade at
            # 20:23-20:26 delayed it further) when this rule killed it for
            # having no window. Its replacement opened VLC in 2 SECONDS,
            # because by then the lane held 300 s -- proof the original was
            # correct and merely early.
            young = self.gen_media_s() <= self.a.lag + self.a.lead + 60.0
            if young:
                o["t"] = None
            if aged and not young and not self.lane_stalled():
                o["t"] = o["t"] or time.time()
                if time.time() - o["t"] > 300 and self.may_restart("viewer"):
                    self.obs["no_vlc"] = {"t": None}
                    fresh = "TS flowing" if (
                        mt is not None and time.time() - mt < 120) else (
                        "band down -- window must wait, not vanish")
                    self.kill_tree(v["pid"],
                                   f"viewer alive but no player for 5 min "
                                   f"({fresh}) -- restarting viewer")
                    self.note_restart("viewer")
                    self.spawn_role("viewer")
            else:
                o["t"] = None
        else:
            self.http_fails = 0

        # 7. maintenance: an unbounded sink eventually fills the disk;
        #    restart the viewer (fresh sink) at TS_ROLL_BYTES. Rare
        #    (~days), logged, and deliberately outside the stampede
        #    counter's semantics -- but still uses it for safety.
        _, sz = self.ts_stat()
        if sz is not None and sz > TS_ROLL_BYTES and roles["viewer"]:
            if self.may_restart("viewer"):
                self.kill_tree(roles["viewer"][0]["pid"],
                               f"sink at {sz/1e9:.1f} GB -- maintenance "
                               f"viewer restart (fresh sink)")
                for e in roles["vlc"]:
                    self.kill_tree(e["pid"], "player follows its sink")
                self.note_restart("viewer")
                self.spawn_role("viewer")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--live-dir", default="data/e31")
    ap.add_argument("--rf", type=int, required=True,
                    help="RF channel number (required: a shipped default is wrong for every reader)")
    ap.add_argument("--lag", type=float, default=120.0)
    ap.add_argument("--lead", type=float, default=60.0)
    ap.add_argument("--once", action="store_true",
                    help="one census/action pass, then exit (gates)")
    a = ap.parse_args()
    if not os.path.isabs(a.live_dir):
        a.live_dir = os.path.join(ROOT, a.live_dir)

    w = Warden(a)
    lock = os.path.join(w.wdir, "warden.lock")
    try:
        d = json.load(open(lock))
        if psutil.pid_exists(d.get("pid", -1)):
            try:
                other = psutil.Process(d["pid"])
                if "atsc3_warden" in " ".join(other.cmdline()):
                    print(f"another warden runs (pid {d['pid']}) -- exiting")
                    return 3
            except psutil.Error:
                pass
    except (OSError, ValueError):
        pass
    json.dump({"pid": os.getpid(), "started": time.time()}, open(lock, "w"))

    w.log(f"warden up: live-dir {a.live_dir} rf {a.rf} "
          f"(loop {LOOP_S:.0f}s, stampede {STAMPEDE_N}/{STAMPEDE_WIN/60:.0f}min)")
    try:
        while True:
            t0 = time.time()
            try:
                w.loop()
            except Exception as e:                             # noqa: BLE001
                # crash-only: one bad loop must not kill the layer
                w.log(f"LOOP ERROR: {type(e).__name__}: {e}")
            if a.once:
                return 0
            time.sleep(max(1.0, LOOP_S - (time.time() - t0)))
    except KeyboardInterrupt:
        w.log("warden stopped (SIGINT)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
