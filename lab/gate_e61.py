#!/usr/bin/env python3
"""gate_e61.py -- the stability campaign's gates (E61).

OFFLINE legs (default; scratch dirs, real archived fragments, no radio):
  1  wav 4GiB wall: header clamps (never raises), close idempotent,
     archive-roll renames the wav+state+sidecar trio, numbering advances
  2  chain tree-kill: stop_tree reaps the child a bare terminate() leaks
     (the leak IS the negative control, measured before the fix each run)
  3  starved-audio bypass ON: dead worker wav, TS still grows near real
     time (the E42 metronome is gone)
  4  starved-audio bypass OFF (--no-audio-starve-bypass): the metronome is
     REPRODUCED -- the gate can fail, and the old code does
  5  E42 fade + multi-roll recovery: viewer idles across a fade with two
     generation rolls, TS resumes within 60 s of recovery (and provably
     did NOT grow during the fade -- the detector can fail)
  6  player-vanish: kill the viewer's (headless) VLC by exact pid; the
     viewer respawns it at the live edge within 75 s and never exits

LIVE legs (--live; requires the deployed warden + stack on data/e31):
  7  fault injection per role -- player, viewer, eng worker, chain child
     -- each restored inside its budget with one-each held throughout
  8  negative control: warden paused -> a killed eng worker STAYS dead
     (proves the restores were the warden's doing), then unpaused ->
     restored

    python lab/gate_e61.py            # offline legs
    python lab/gate_e61.py --live     # fault injection on the real stack
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tools"))

import atsc3_audio as AU                                        # noqa: E402
import atsc3_run as RUN                                         # noqa: E402

try:
    import psutil
except ImportError:
    psutil = None

PY = sys.executable
E31 = os.path.join(ROOT, "data", "e31")
SRC_V = os.path.join(E31, "live_video_pid12.0001.m4s")
SRC_VI = os.path.join(E31, "live_video_pid12.0001.idx")
SRC_A = os.path.join(E31, "live_audio_pid13.0001.m4s")
SRC_AI = os.path.join(E31, "live_audio_pid13.0001.idx")
SCRATCH = os.path.join(ROOT, "lab", "e61_scratch")
RUN_TAG = time.strftime("%m%d_%H%M%S")   # unique per run: NOTHING under
#                                          lab/ is ever deleted (fleet law)

FAILS = []
SKIPPED = []


def band_up(window=12):
    """Is RF33 actually delivering right now? (chain's own FEC line)"""
    p = os.path.join(E31, "chain.log")
    try:
        with open(p, "rb") as f:
            f.seek(max(0, os.path.getsize(p) - 100_000))
            txt = f.read().decode("utf-8", "replace")
    except OSError:
        return False
    vals = [float(v) for v in
            re.findall(r"\[\s*([\d.]+)% now\]", txt)][-window:]
    return bool(vals) and sum(vals) / len(vals) > 60.0


def scratch_dir(name):
    d = os.path.join(SCRATCH, f"{name}_{RUN_TAG}")
    os.makedirs(os.path.join(d, "_tv"), exist_ok=True)
    return d


def check(name, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}"
          + (f"  ({detail})" if detail else ""), flush=True)
    if not ok:
        FAILS.append(name)


def env_pinned():
    env = dict(os.environ)
    for k in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
        env[k] = "1"
    return env


# ---------------------------------------------------------------- leg 1

def leg_wav_wall():
    print("leg 1: the 4GiB wav wall")
    d = os.path.join(SCRATCH, "wall")
    os.makedirs(d, exist_ok=True)
    p = os.path.join(d, "live_audio.wav")
    w = AU.WavAppender(p, ch=2)
    w.n = 2_200_000_000            # past 2^32 bytes of data
    try:
        w._header()
        ok = True
    except OverflowError:
        ok = False
    check("header write past 2^32 bytes never raises", ok)
    w.close()
    w.close()
    check("close is idempotent (the crash re-raised inside close)", True)
    # trio archive
    base = os.path.splitext(p)[0]
    open(base + ".state.json", "w").write("{}")
    open(base + ".json", "w").write("{}")
    tag = AU.archive_roll(p)
    check("archive_roll moves wav+state+sidecar", tag is not None
          and os.path.exists(f"{base}.{tag}.wav")
          and os.path.exists(f"{base}.{tag}.state.json")
          and os.path.exists(f"{base}.{tag}.json")
          and not os.path.exists(p), f"tag {tag}")
    open(p, "wb").write(b"RIFF")
    tag2 = AU.archive_roll(p)
    check("roll numbering advances", tag2 not in (None, tag), f"{tag2}")


# ---------------------------------------------------------------- leg 2

def leg_tree_kill():
    print("leg 2: chain tree-kill vs the orphan leak")
    if psutil is None:
        check("psutil available", False)
        return
    code = ("import subprocess,sys,time;"
            "p=subprocess.Popen([sys.executable,'-c',"
            "'import time;time.sleep(300)']);"
            "print(p.pid,flush=True);time.sleep(300)")

    def spawn_pair():
        par = subprocess.Popen([PY, "-c", code], stdout=subprocess.PIPE,
                               text=True)
        kid = int(par.stdout.readline().strip())
        return par, kid

    # NEGATIVE CONTROL: the pre-E61 kill (terminate parent only)
    par, kid = spawn_pair()
    par.terminate()
    par.wait(timeout=10)
    time.sleep(1.0)
    leaked = psutil.pid_exists(kid)
    check("negative control: bare terminate() LEAKS the child", leaked,
          f"child {kid} survived={leaked}")
    if leaked:
        psutil.Process(kid).kill()
    # THE FIX
    par, kid = spawn_pair()
    RUN.stop_tree(par)
    time.sleep(1.0)
    check("stop_tree reaps parent AND child",
          par.poll() is not None and not psutil.pid_exists(kid))


# ------------------------------------------------------------- lane sim

class LaneSim:
    """Replays a real archived video (+audio) lane into a scratch live
    dir at an accelerated rate, with fades and generation rolls."""

    def __init__(self, d, with_audio=False):
        self.d = d
        os.makedirs(os.path.join(d, "_tv"), exist_ok=True)
        self.vsrc = open(SRC_V, "rb").read()
        self.vidx = [json.loads(l) for l in open(SRC_VI)]
        self.gen = 1
        self.k = 0                 # next fragment index
        self.with_audio = with_audio
        if with_audio:
            self.asrc = open(SRC_A, "rb").read()
            self.aidx = [json.loads(l) for l in open(SRC_AI)]
            self.amap = {r["seq"]: r for r in self.aidx}
        self.vp = os.path.join(d, "live_video_pid12.m4s")
        self.ap = os.path.join(d, "live_audio_pid13.m4s")
        self._fresh_files()

    def _fresh_files(self):
        init = self.vidx[0]["off"]
        open(self.vp, "wb").write(self.vsrc[:init])
        open(os.path.splitext(self.vp)[0] + ".idx", "w").close()
        if self.with_audio:
            ainit = self.aidx[0]["off"]
            open(self.ap, "wb").write(self.asrc[:ainit])
            open(os.path.splitext(self.ap)[0] + ".idx", "w").close()
        self.voff = init
        self.aoff = self.aidx[0]["off"] if self.with_audio else 0
        self.first = None

    def roll(self):
        """Archive current files, bump generation, keep seq continuous."""
        n = self.gen
        for p in ([self.vp, self.ap] if self.with_audio else [self.vp]):
            for ext in ("", ".idx"):
                src = os.path.splitext(p)[0] + ext if ext else p
                if not os.path.exists(src):
                    continue
                dst = f"{os.path.splitext(p)[0]}.{n:04d}{ext or '.m4s'}"
                for _ in range(20):     # viewer may hold it open briefly
                    try:
                        os.rename(src, dst)
                        break
                    except OSError:
                        time.sleep(0.5)
        self.gen += 1
        self._fresh_files()
        self.write_json()

    def feed(self, n=1):
        for _ in range(n):
            if self.k >= len(self.vidx):
                return False
            r = self.vidx[self.k]
            frag = self.vsrc[r["off"]:r["off"] + r["len"]]
            with open(self.vp, "ab") as f:
                rec = {"seq": r["seq"], "off": self.voff, "len": r["len"]}
                f.write(frag)
            with open(os.path.splitext(self.vp)[0] + ".idx", "a") as f:
                f.write(json.dumps(rec) + "\n")
            if self.first is None:
                self.first = r["seq"]
            self.voff += r["len"]
            if self.with_audio and r["seq"] in self.amap:
                ar = self.amap[r["seq"]]
                afrag = self.asrc[ar["off"]:ar["off"] + ar["len"]]
                with open(self.ap, "ab") as f:
                    arec = {"seq": ar["seq"], "off": self.aoff,
                            "len": ar["len"]}
                    f.write(afrag)
                with open(os.path.splitext(self.ap)[0] + ".idx", "a") as f:
                    f.write(json.dumps(arec) + "\n")
                self.aoff += ar["len"]
            self.k += 1
            self.last = r["seq"]
        self.write_json()
        return True

    def write_json(self):
        lanes = {"video_pid12": {
            "bytes": self.voff, "path": self.vp, "kind": "video", "pid": 12,
            "first_seq": self.first or self.vidx[0]["seq"],
            "last_seq": getattr(self, "last", self.vidx[0]["seq"]),
            "init_bytes": self.vidx[0]["off"],
            "idx": os.path.splitext(self.vp)[0] + ".idx",
            "generation": self.gen, "handler": "vide"}}
        if self.with_audio:
            lanes["audio_pid13"] = {
                "bytes": self.aoff, "path": self.ap, "kind": "audio",
                "pid": 13, "first_seq": self.first or self.aidx[0]["seq"],
                "last_seq": getattr(self, "last", self.aidx[0]["seq"]),
                "init_bytes": self.aidx[0]["off"],
                "idx": os.path.splitext(self.ap)[0] + ".idx",
                "generation": self.gen, "handler": "soun"}
        json.dump({"bytes": self.voff, "updated": time.time(),
                   "media_s": 0.0, "lanes": lanes},
                  open(os.path.join(self.d, "live.json"), "w"))


def ts_size(d):
    try:
        return os.path.getsize(os.path.join(d, "_tv", "live_tv2.ts"))
    except OSError:
        return 0


def spawn_viewer(d, extra=(), log_name="viewer.log"):
    lf = open(os.path.join(d, log_name), "ab", buffering=0)
    return subprocess.Popen(
        [PY, "-u", os.path.join(ROOT, "tools", "atsc3_tv.py"),
         "--live-dir", d, "--player", "none", "--lag", "8",
         "--audio-hold", "45"] + list(extra),
        cwd=ROOT, stdout=lf, stderr=subprocess.STDOUT, env=env_pinned())


def make_dead_wav(d, sim):
    """A tiny STATIC eng wav + sidecar: enough anchor for audio_ready to
    answer False (the wav will never cover new slots) -- the dead-worker
    shape."""
    import wave as _w
    p = os.path.join(d, "live_audio.wav")
    w = _w.open(p, "wb")
    w.setnchannels(2); w.setsampwidth(2); w.setframerate(48000)
    w.writeframes(b"\x00\x00" * 2 * 48000)          # 1 s of silence
    w.close()
    json.dump({"first_seq": sim.vidx[0]["seq"], "fs": 48000, "pid": 13,
               "path": p, "wav_base": 0, "generation": sim.gen},
              open(os.path.join(d, "live_audio.json"), "w"))


# ---------------------------------------------------------- legs 3 & 4

def leg_radio_yield():
    print("leg 9: radio yields are not faults (and crashes still are)")
    d = scratch_dir("yield")
    log_p = os.path.join(d, "chain.log")

    def write(txt):
        with open(log_p, "ab") as f:
            f.write(txt.encode())
        return os.path.getsize(log_p)

    off0 = write("[07:52] ATSC 3.0 live -- RF33\n")
    off1 = write("RuntimeError: radio held by sonde_rx (priority 100); "
                 "not seizing\n")
    check("rc=1 + a radio-lock refusal reads as a YIELD",
          RUN.radio_yield_reason(d, 1, off0) is not None,
          str(RUN.radio_yield_reason(d, 1, off0)))
    check("rc=3 reads as a YIELD (scheduled reservation)",
          RUN.radio_yield_reason(d, 3, off1) is not None)
    # NEGATIVE CONTROL 1: a real crash must still be counted
    off2 = write("Traceback (most recent call last):\n"
                 "MemoryError: out of memory in the demapper\n")
    check("negative control: rc=1 with a REAL crash is NOT a yield",
          RUN.radio_yield_reason(d, 1, off1) is None)
    # NEGATIVE CONTROL 2: a STALE refusal must not excuse a later crash --
    # the whole point of the since_off bound
    check("negative control: a stale refusal above the offset is ignored",
          RUN.radio_yield_reason(d, 1, off2) is None)
    check("unbounded read WOULD have been fooled (the bug the bound "
          "prevents)", RUN.radio_yield_reason(d, 1, 0) is not None)


class FakeClock:
    """Drop-in for the warden's `time` module: a clock the gate drives."""

    def __init__(self, t0=1_786_000_000.0):
        self.t = t0

    def time(self):
        return self.t

    def sleep(self, s):
        self.t += s

    def strftime(self, fmt, *a):
        return time.strftime(fmt, *a)


class HarnessWarden:
    """The REAL Warden with its hands tied: census is scripted and the
    two effectful methods record instead of acting. Everything between --
    classify, lane_stalled, the no-player decision -- is production code."""

    def __init__(self, warden_mod, live_dir, clock, lag=120.0, lead=60.0):
        class A:
            pass
        a = A()
        a.live_dir, a.rf, a.lag, a.lead = live_dir, 33, lag, lead
        self.mod = warden_mod
        self.clock = clock
        self.w = warden_mod.Warden(a)
        self.kills, self.spawns = [], []
        self.w.kill_tree = lambda pid, why: self.kills.append((pid, why))
        self.w.spawn_role = self._spawn
        self.w.sweep_orphan_forks = lambda census: None
        self.census = []
        self.w.census = lambda: self.census

    def _spawn(self, role):
        """Record -- and MODEL the consequence: a respawned viewer is a
        NEW process with a new start time. Without this the harness leaves
        a permanently-old viewer in the census and the warden re-fires
        against a restart that never happened (caught by this gate's own
        one-restart assertion)."""
        self.spawns.append(role)
        if role == "viewer":
            for e in self.census:
                if "atsc3_tv.py" in e["cmd"]:
                    e["t0"] = self.clock.t
                    e["pid"] += 1000

    def set_roles(self, live_dir, viewer_t0, with_vlc):
        d = live_dir
        sink = os.path.join(d, "_tv", "live_tv2.ts")
        c = [{"pid": 101, "name": "python.exe", "t0": viewer_t0, "proc": None,
              "cmd": f"python -u tools/atsc3_run.py --live-dir {d}"},
             {"pid": 102, "name": "python.exe", "t0": viewer_t0, "proc": None,
              "cmd": f"python -u tools/atsc3_audio.py --live-dir {d} --pid 13"},
             {"pid": 103, "name": "python.exe", "t0": viewer_t0, "proc": None,
              "cmd": f"python -u tools/atsc3_audio.py --live-dir {d} --pid 14"},
             {"pid": 104, "name": "python.exe", "t0": viewer_t0, "proc": None,
              "cmd": f"python -u tools/atsc3_subs.py --live-dir {d}"},
             {"pid": 105, "name": "python.exe", "t0": viewer_t0, "proc": None,
              "cmd": f"python -u tools/atsc3_tv.py --live-dir {d}"}]
        if with_vlc:
            c.append({"pid": 106, "name": "vlc.exe", "t0": viewer_t0,
                      "proc": None, "cmd": f'vlc.exe "{sink}"'})
        self.census = c


def write_lanes(d, last_seq, gen=1, media_slots=400):
    """live.json with a video lane whose generation is old enough that the
    lane-age guard is satisfied (so the test isolates the player logic)."""
    lanes = {}
    for name, kind, pid in (("video_pid12", "video", 12),
                            ("audio_pid13", "audio", 13),
                            ("audio_pid14", "audio", 14),
                            ("subs_pid15", "subs", 15)):
        lanes[name] = {"bytes": 10_000_000, "kind": kind, "pid": pid,
                       "first_seq": last_seq - media_slots,
                       "last_seq": last_seq, "generation": gen,
                       "path": os.path.join(d, f"live_{name}.m4s"),
                       "idx": os.path.join(d, f"live_{name}.idx")}
    json.dump({"bytes": 10_000_000, "updated": time.time(), "media_s": 800.0,
               "lanes": lanes}, open(os.path.join(d, "live.json"), "w"))


def old_ts_fresh_would_fire(d, now, viewer_t0):
    """The PRE-E62 predicate, verbatim in spirit: restart only when the
    viewer is aged AND the TS is fresh. Kept executable so the gate has a
    real negative control instead of a claim about deleted code."""
    try:
        mt = os.path.getmtime(os.path.join(d, "_tv", "live_tv2.ts"))
    except OSError:
        return False
    return (now - viewer_t0 > 300) and (now - mt < 120)


def leg_e62b():
    print("leg 10: E62b -- the no-player rule keys on the CHAIN, not the TS")
    sys.path.insert(0, os.path.join(ROOT, "tools"))
    import atsc3_warden as W                                    # noqa: E402
    d = scratch_dir("e62b")
    # a TS that exists but is STALE -- exactly the shape a broken viewer
    # leaves behind, and the shape the old predicate could not act on
    ts = os.path.join(d, "_tv", "live_tv2.ts")
    open(ts, "wb").write(b"\x47" * 1000)
    clock = FakeClock()
    real_time = W.time
    W.time = clock
    try:
        h = HarnessWarden(W, d, clock)
        viewer_t0 = clock.t - 3600          # long-established viewer
        os.utime(ts, (clock.t - 4000, clock.t - 4000))   # stale by an hour

        # ---- PHASE A: lane STALLED (radio held by the sonde), no player.
        # Nothing to show: the warden must not touch the viewer.
        write_lanes(d, last_seq=225800900)
        h.set_roles(d, viewer_t0, with_vlc=False)
        for _ in range(40):                 # 40 x 60 s = 40 min of stall
            h.w.loop()
            clock.sleep(60)
        check("lane stalled + no player: viewer NOT restarted (no thrash)",
              not h.kills and not h.spawns,
              f"kills {h.kills}, spawns {h.spawns}")

        # ---- PHASE B: the lane advances and the viewer is MUXING FINE
        # (TS growing) but the window is gone. Nothing else is wrong, so
        # only the no-player rule can act -- the rule under test, isolated.
        h.kills.clear(); h.spawns.clear()
        seq = 225800900
        fired_at = None
        t_start = clock.t
        for _ in range(20):                 # up to 20 min
            seq += 30                       # ~60 s of air per step
            write_lanes(d, last_seq=seq)
            os.utime(ts, (clock.t, clock.t))        # viewer is muxing
            h.w.loop()
            if h.kills and fired_at is None:
                fired_at = clock.t - t_start
                # the respawned viewer does its job: window comes back.
                # (Leaving it windowless would make further restarts
                # CORRECT escalation, which tests nothing about recovery.)
                h.set_roles(d, clock.t, with_vlc=True)
            clock.sleep(60)
        check("lane advancing + TS growing + no player: viewer restarted",
              fired_at is not None,
              f"after {fired_at}s" if fired_at else "never fired")
        check("restore is prompt (<= 10 min of advancing lane)",
              fired_at is not None and fired_at <= 600, f"{fired_at}s")
        check("exactly one restart, and for the no-player reason",
              len(h.kills) == 1 and h.spawns == ["viewer"]
              and "no player" in h.kills[0][1],
              f"kills {[k[1][:40] for k in h.kills]}, spawns {h.spawns}")
        check("once the window is back the warden stops intervening",
              len(h.kills) == 1, f"{len(h.kills)} kills over 20 min")

        # ---- PHASE C: the DISCRIMINATING case -- lane advancing while the
        # TS is stale (a viewer broken outright). This is the 02:00-07:18
        # shape. The warden must act; the pre-E62 predicate cannot, because
        # it asks the broken viewer's own output for permission.
        h.kills.clear(); h.spawns.clear()
        h.set_roles(d, clock.t - 3600, with_vlc=False)
        h.w.obs.pop("no_vlc", None)
        h.w.state["restarts"] = {}          # fresh stampede window
        os.utime(ts, (clock.t - 4000, clock.t - 4000))
        old_fired = False
        acted = False
        for _ in range(12):
            seq += 30
            write_lanes(d, last_seq=seq)
            h.w.loop()
            if h.kills:
                acted = True
            if old_ts_fresh_would_fire(d, clock.t, clock.t - 3600):
                old_fired = True
            clock.sleep(60)
        check("lane advancing + TS stale: the warden ACTS", acted,
              f"kills {[k[1][:45] for k in h.kills]}")
        check("NEGATIVE CONTROL: the pre-E62 ts_fresh rule never fires on "
              "that same trace (the circularity, demonstrated)",
              not old_fired)

        # ---- PHASE D: a healthy stack must be left alone -- the gate has
        # to be able to fail in both directions
        h.kills.clear(); h.spawns.clear()
        h.set_roles(d, clock.t - 3600, with_vlc=True)
        h.w.obs.pop("no_vlc", None)
        h.w.state["restarts"] = {}
        for _ in range(20):
            seq += 30
            write_lanes(d, last_seq=seq)
            os.utime(ts, (clock.t, clock.t))
            h.w.loop()
            clock.sleep(60)
        check("healthy stack (lane advancing, TS growing, player up): "
              "nothing touched", not h.kills and not h.spawns,
              f"kills {h.kills}, spawns {h.spawns}")
    finally:
        W.time = real_time


def leg_judge_clock():
    print("leg 11: the acceptance bar itself can fail")
    sys.path.insert(0, os.path.join(ROOT, "tools"))
    import atsc3_judge as J                                     # noqa: E402
    d = scratch_dir("judge")

    class A:
        pass
    a = A()
    a.live_dir = d
    j = J.Judge(a)

    # ---- QUALIFIED time only accrues while the chain delivers.
    # fec_5min_mean is the band oracle; drive it directly.
    j.fec_5min_mean = lambda: 0.0
    j.qualified = 0.0
    j._last_tick = None
    for _ in range(10):                     # 10 min of DEAD band
        now = time.time()
        dt = 0.0 if j._last_tick is None else 60.0
        j._last_tick = now
        if (j.fec_5min_mean() or 0) > J.FEC_GATE:
            j.qualified += dt
    check("dead band banks NO qualified time", j.qualified == 0.0,
          f"{j.qualified}s")

    j.fec_5min_mean = lambda: 99.0
    for _ in range(10):                     # 10 min of REAL television
        now = time.time()
        dt = 60.0
        j._last_tick = now
        if (j.fec_5min_mean() or 0) > J.FEC_GATE:
            j.qualified += dt
    check("delivering band banks qualified time", j.qualified == 600.0,
          f"{j.qualified}s")

    # ---- a VIOLATION must restart the clock (the bug: it did not)
    banked = j.qualified
    j.violate("test", "synthetic")
    check("a violation RESETS the qualified clock", j.qualified == 0.0,
          f"was {banked}s, now {j.qualified}s")
    check("the violation was counted", j.violations == 1)

    # ---- NEGATIVE CONTROL: the pre-fix violate() reset only soak_start,
    # so the bar it actually reads (qualified) survived every failure.
    class OldJudge(J.Judge):
        def violate(self, what, detail):          # pre-E63 body
            self.violations += 1
            self.soak_start = time.time()
            self.passed_logged = False
    oj = OldJudge(a)
    oj.qualified = 600.0
    oj.violate("test", "synthetic")
    check("NEGATIVE CONTROL: pre-fix violate() left the bar untouched "
          "(6 h 'clean' was reachable after failures)",
          oj.qualified == 600.0, f"{oj.qualified}s survived the violation")

    # ---- dt capping: a suspended judge must not credit the gap
    j2 = J.Judge(a)
    j2.fec_5min_mean = lambda: 99.0
    j2._last_tick = time.time() - 86400        # a day-long suspension
    dt = 0.0 if j2._last_tick is None \
        else min(time.time() - j2._last_tick, 2 * J.SAMPLE_S)
    check("a suspension credits at most one sample window, not the gap",
          dt <= 2 * J.SAMPLE_S, f"{dt:.0f}s credited from an 86400s gap")


def _mk_noise_wav(path, secs, seed, fs=48000):
    subprocess.run(
        [G_FF("ffmpeg"), "-y", "-loglevel", "error", "-f", "lavfi", "-i",
         f"anoisesrc=d={secs}:c=pink:r={fs}:a=0.5:seed={seed}",
         "-ac", "1", "-c:a", "pcm_s16le", path], check=True)


def _mk_ts_from_wav(wav, ts):
    """The real viewer's audio path in miniature: wav -> mp2 224k -> TS."""
    subprocess.run(
        [G_FF("ffmpeg"), "-y", "-loglevel", "error", "-i", wav,
         "-c:a", "mp2", "-b:a", "224k", "-f", "mpegts", ts], check=True)


def G_FF(name):
    sys.path.insert(0, os.path.join(ROOT, "tools"))
    import atsc3_tv as _tv
    return _tv.ffbin(name)


def leg_referee():
    print("leg 12: the A/V referee actually works (it had never once run)")
    sys.path.insert(0, os.path.join(ROOT, "tools"))
    import atsc3_judge as J                                     # noqa: E402
    d = scratch_dir("referee")
    os.makedirs(os.path.join(d, "_warden"), exist_ok=True)
    wav = os.path.join(d, "live_audio.wav")
    ts = os.path.join(d, "_tv", "live_tv2.ts")

    print("  building 30 min of pink noise + its mp2 TS ...", flush=True)
    _mk_noise_wav(wav, 1800, seed=1234)
    _mk_ts_from_wav(wav, ts)
    sz = os.path.getsize(ts)
    check("synthetic TS clears the referee's 40 MB floor", sz > 40e6,
          f"{sz/1e6:.0f} MB")

    class A:
        pass
    a = A()
    a.live_dir = d
    j = J.Judge(a)
    verdict, detail = j.referee()
    check("referee returns a VERDICT on real bytes (not a skip)",
          verdict in ("pass", "fail"), f"{verdict} -- {detail}")
    check("matching audio correlates above the CORR_MIN bar",
          verdict == "pass", detail)

    # ---- NEGATIVE CONTROL: a TS whose audio is DIFFERENT air entirely.
    # If this also 'passed', the invariant would be decorative.
    wav2 = os.path.join(d, "other.wav")
    _mk_noise_wav(wav2, 1800, seed=999)
    _mk_ts_from_wav(wav2, ts)                  # TS no longer matches the wav
    verdict2, detail2 = j.referee()
    check("NEGATIVE CONTROL: mismatched audio FAILS the referee",
          verdict2 == "fail", f"{verdict2} -- {detail2}")

    # ---- silence must SKIP, never fail (a fade is not a desync)
    subprocess.run(
        [G_FF("ffmpeg"), "-y", "-loglevel", "error", "-f", "lavfi", "-i",
         "anullsrc=r=48000:cl=mono:d=1800", "-c:a", "mp2", "-b:a", "224k",
         "-f", "mpegts", ts], check=True)
    verdict3, detail3 = j.referee()
    check("silence SKIPS (a fade is not a desync)", verdict3 == "skip",
          f"{verdict3} -- {detail3}")


def leg_fade_not_wedge():
    print("leg 13: a fade is not a wedge (and a wedge is still a wedge)")
    d = scratch_dir("fade")
    log_p = os.path.join(d, "chain.log")

    def status_line(t, secs, frames, fec):
        return (f"[{t}]   {secs}s   {frames} Frames   0.42x rt  inst  0.44x  "
                f"FEC 0/6216 [ {fec}% now]  q raw/frm/bb 0/6/0\n")

    def write(lines):
        with open(log_p, "ab") as f:
            for l in lines:
                f.write(l.encode())
        return os.path.getsize(log_p)

    # ---- CASE 1: the air is gone, the chain is fine (frames climbing)
    off0 = os.path.getsize(log_p) if os.path.exists(log_p) else 0
    write([status_line("08:2%d:00" % (i % 6), 50 + 5 * i, 100 + 20 * i, "0.0")
           for i in range(12)])
    fr1, fec1 = RUN.chain_status(d, off0)
    check("fade: FEC mean reads ~0", fec1 is not None and fec1 < RUN.FADE_FEC,
          f"fec {fec1}")
    fr1b, _ = RUN.chain_status(d, off0)
    write([status_line("08:26:05", 115, 340, "0.0")])
    fr1c, fec1c = RUN.chain_status(d, off0)
    check("fade: the Frames counter is still ADVANCING (chain alive)",
          fr1c > fr1b, f"{fr1b} -> {fr1c}")
    check("fade verdict: WAIT, do not restart",
          fec1c < RUN.FADE_FEC and fr1c != fr1b)

    # ---- NEGATIVE CONTROL A: frames advancing but FEC HIGH. That is the
    # 8/07 flipped-bit shape (chain runs, assets dead) -- must NOT be
    # excused as a fade.
    off1 = os.path.getsize(log_p)
    write([status_line("08:30:00", 200 + 5 * i, 900 + 20 * i, "98.5")
           for i in range(12)])
    fr2, fec2 = RUN.chain_status(d, off1)
    check("NEGATIVE CONTROL: FEC high + lanes dead is NOT a fade "
          "(the flipped-bit case still restarts)",
          fec2 is not None and not (fec2 < RUN.FADE_FEC), f"fec {fec2}")

    # ---- NEGATIVE CONTROL B: frames FROZEN (livelock) -- the original
    # wedge this supervisor exists for. Same FEC ~0, but no motion.
    off2 = os.path.getsize(log_p)
    write([status_line("08:40:00", 300, 1500, "0.0") for _ in range(12)])
    fr3a, fec3 = RUN.chain_status(d, off2)
    write([status_line("08:40:05", 300, 1500, "0.0")])
    fr3b, _ = RUN.chain_status(d, off2)
    check("NEGATIVE CONTROL: frames FROZEN is a wedge, not a fade "
          "(restart still fires)", fr3a == fr3b, f"{fr3a} == {fr3b}")

    # ---- the since_off bound: an old fade must not colour a new instance
    off3 = os.path.getsize(log_p)
    check("no status lines yet for a fresh instance -> no verdict",
          RUN.chain_status(d, off3) == (None, None))


def leg_rf_referee():
    print("leg 14: the independent RF referee decides air-vs-us")
    sys.path.insert(0, os.path.join(ROOT, "tools"))
    import atsc3_judge as J                                     # noqa: E402
    d = scratch_dir("refrf")

    class A:
        pass
    a = A()
    a.live_dir = d

    def judge_with(ours, ref):
        j = J.Judge(a)
        j.fec_5min_mean = lambda: ours
        j.rf_referee = lambda: ref
        # a unit test must control EVERY input: without this the real radio
        # lock leaks in and the case silently becomes a different case
        # (caught 8/10 -- sonde_rx held the radio and suppressed the very
        # violation this leg exists to prove).
        j.radio_yielded = lambda: None
        j.last_rf = 0.0
        fired = []
        j.violate = lambda w, det: fired.append((w, det))
        for _ in range(J.HDHR_CONFIRM + 1):
            j.last_rf = 0.0                    # force a read every pass
            j.sample()
        return fired

    # OURS DEAD, AIR FINE -> it is us. Must violate.
    f = judge_with(5.0, (97, 66, "atsc3"))
    check("ours 5% + referee snq 97 -> VIOLATION (the fault is ours)",
          any(w == "rf-referee" for w, _ in f), str(f[:1]))

    # NEGATIVE CONTROL 1: both dead -> the air. Must NOT violate.
    f = judge_with(5.0, (0, 20, "none"))
    check("NEGATIVE CONTROL: ours 5% + referee snq 0 -> excused (the air)",
          not any(w == "rf-referee" for w, _ in f), str(f))

    # NEGATIVE CONTROL 2: both fine -> nothing to say.
    f = judge_with(99.0, (97, 66, "atsc3"))
    check("NEGATIVE CONTROL: both healthy -> no violation",
          not any(w == "rf-referee" for w, _ in f), str(f))

    # NEGATIVE CONTROL 3: a MISSING instrument is never a verdict.
    f = judge_with(5.0, None)
    check("NEGATIVE CONTROL: referee unreachable -> skip, never violate",
          not any(w == "rf-referee" for w, _ in f), str(f))

    # confirmation: one disagreement alone must not fire
    j = J.Judge(a)
    j.fec_5min_mean = lambda: 5.0
    j.rf_referee = lambda: (97, 66, "atsc3")
    j.radio_yielded = lambda: None
    fired = []
    j.violate = lambda w, det: fired.append(w)
    j.last_rf = 0.0
    j.sample()
    check("a SINGLE disagreement does not fire (confirmation required)",
          not fired, f"disagree={j.rf_disagree}")

    # E72: NEGATIVE CONTROL -- we are not decoding because we HANDED the
    # antenna to a higher-priority service. Obeying the single-tenant rule
    # must never be scored as a stack fault, however good the air looks.
    def judge_yielded(ours, ref, who):
        jj = J.Judge(a)
        jj.fec_5min_mean = lambda: ours
        jj.rf_referee = lambda: ref
        jj.radio_yielded = lambda: who
        out = []
        jj.violate = lambda w, det: out.append(w)
        for _ in range(J.HDHR_CONFIRM + 2):
            jj.last_rf = 0.0
            jj.sample()
        return out
    f = judge_yielded(5.0, (100, 67, "atsc3"), "sonde_rx (priority 100)")
    check("NEGATIVE CONTROL: radio yielded to a higher priority -> the RF "
          "referee is SUSPENDED, not violated",
          not any(w == "rf-referee" for w, *_ in
                  [(x,) for x in f]) if f else True, str(f))
    # and it must resume the moment the radio comes back
    f = judge_yielded(5.0, (100, 67, "atsc3"), None)
    check("radio ours again -> the referee resumes and DOES violate",
          any(x == "rf-referee" for x in f), str(f))

    # E72b NEGATIVE CONTROL: the 5-min mean still carries a dip we have
    # ALREADY recovered from. An accusation must be present tense.
    def judge_recent(m5, recent):
        jj = J.Judge(a)
        jj.fec_5min_mean = lambda: m5
        jj.fec_recent = lambda w=90: recent
        jj.rf_referee = lambda: (100, 67, "atsc3")
        jj.radio_yielded = lambda: None
        out = []
        jj.violate = lambda w, det: out.append(w)
        for _ in range(J.HDHR_CONFIRM + 2):
            jj.last_rf = 0.0
            jj.sample()
        return out
    check("NEGATIVE CONTROL: 5-min mean 39% but recovered NOW (100%) -> "
          "no violation", not judge_recent(39.0, 100.0))
    check("still failing NOW (both low) -> violation stands",
          any(x == "rf-referee" for x in judge_recent(39.0, 5.0)))


def leg_maintenance_window():
    print("leg 15: the maintenance window expires (it cannot buy silence)")
    sys.path.insert(0, os.path.join(ROOT, "tools"))
    import atsc3_judge as J                                     # noqa: E402
    import atsc3_warden as W                                    # noqa: E402
    d = scratch_dir("mw")
    wdir = os.path.join(d, "_warden")
    os.makedirs(wdir, exist_ok=True)
    pause = os.path.join(wdir, "pause")

    class A:
        pass
    a = A()
    a.live_dir, a.rf, a.lag, a.lead = d, 33, 120.0, 60.0
    j = J.Judge(a)

    # no marker -> judging is ON
    check("no marker: judge is active", j.maintenance() is None)

    # fresh marker -> window honoured, with its reason
    open(pause, "w").write("E68 fade capture")
    why = j.maintenance()
    check("fresh marker: window honoured and carries its reason",
          why == "E68 fade capture", str(why))

    # STALE marker -> ignored. The safety property that matters overnight.
    old = time.time() - (J.PAUSE_MAX_S + 120)
    os.utime(pause, (old, old))
    check("NEGATIVE CONTROL: a STALE marker is IGNORED -- judging resumes "
          "on its own", j.maintenance() is None,
          f"age {(time.time()-old)/60:.0f} min > {J.PAUSE_MAX_S/60:.0f}")

    # the warden must agree, or the two would disagree about reality
    wd = W.Warden(a)
    check("warden also ignores the stale marker (supervision resumes)",
          not wd.paused())
    now = time.time()
    os.utime(pause, (now, now))
    check("warden honours a fresh marker", wd.paused())

    # the qualified clock must not advance across a window
    j.qualified = 100.0
    j._last_tick = time.time() - 60.0
    banked_before = j.qualified
    # simulate the loop's maintenance branch
    j._last_tick = None
    check("qualified clock is frozen during the window (no time banked)",
          j.qualified == banked_before, f"{j.qualified}s")
    os.rename(pause, pause + ".done")       # never delete


def leg_av_skew():
    print("leg 16: the video-vs-audio skew check (the hole E71 fell through)")
    sys.path.insert(0, os.path.join(ROOT, "tools"))
    import atsc3_judge as J                                     # noqa: E402
    d = scratch_dir("skew")
    os.makedirs(os.path.join(d, "_tv"), exist_ok=True)
    os.makedirs(os.path.join(d, "_warden"), exist_ok=True)
    ts = os.path.join(d, "_tv", "live_tv2.ts")

    class A:
        pass
    a = A()
    a.live_dir = d
    j = J.Judge(a)

    def build(voff):
        """A TS whose VIDEO is voff seconds from its own audio."""
        subprocess.run(
            [G_FF("ffmpeg"), "-y", "-loglevel", "error",
             "-itsoffset", str(voff),
             "-f", "lavfi", "-i", "testsrc=size=320x240:rate=30:duration=200",
             "-f", "lavfi", "-i", "sine=frequency=440:duration=200",
             "-c:v", "mpeg2video", "-b:v", "2M", "-c:a", "mp2", "-b:a", "192k",
             "-f", "mpegts", ts], check=True)

    build(0.0)
    sk = j.av_skew()
    check("aligned mux reads ~0 skew", sk is not None and abs(sk) < 1.0,
          f"{sk if sk is None else round(sk, 3)}s")

    # NEGATIVE CONTROL: the exact fault E71 found -- video ahead of audio.
    # If this does not trip, the invariant is decorative.
    build(6.0)
    sk6 = j.av_skew()
    check("NEGATIVE CONTROL: video pushed +6s IS detected",
          sk6 is not None and abs(sk6) > J.SKEW_MAX,
          f"{sk6 if sk6 is None else round(sk6, 3)}s > {J.SKEW_MAX}s")
    check("the sign points the right way (video AHEAD reads positive)",
          sk6 is not None and sk6 > 0, f"{sk6 if sk6 is None else round(sk6,3)}")

    # and the 36 s shape that actually happened
    build(36.061)
    sk36 = j.av_skew()
    check("the real E71 offset (+36.061s) is detected",
          sk36 is not None and abs(sk36) > J.SKEW_MAX,
          f"{sk36 if sk36 is None else round(sk36, 2)}s")

    # confirmation: one bad reading must not fire on its own
    j.skew_bad = 0
    j.last_skew = 0.0
    fired = []
    j.violate = lambda w, det: fired.append(w)
    j.sample()
    check("a SINGLE bad skew does not fire (confirmation required)",
          not fired, f"skew_bad={j.skew_bad}")


def leg_ts_growth_windows():
    print("leg 20: TS-growth must not be judged on time it was not owed")
    sys.path.insert(0, os.path.join(ROOT, "tools"))
    import atsc3_judge as J                                     # noqa: E402
    d = scratch_dir("tsgrow")
    os.makedirs(os.path.join(d, "_tv"), exist_ok=True)
    ts = os.path.join(d, "_tv", "live_tv2.ts")
    open(ts, "wb").write(b"\x47" * 1000)

    class A:
        pass
    a = A()
    a.live_dir = d

    def run(seq, yielded=None):
        """seq: list of (fec, ts_size). -> violations raised."""
        j = J.Judge(a)
        j.rf_referee = lambda: None
        j.radio_yielded = lambda: yielded
        j.referee = lambda: ("skip", "")
        j.av_skew = lambda: 0.0
        j.census_roles = lambda: {k: [(1, time.time())] for k in
                                 ("run", "eng", "spa", "subs", "viewer",
                                  "vlc", "warden")}
        j.lanes = lambda: {"a": {"pid": 13}, "b": {"pid": 14},
                           "c": {"kind": "subs"}}
        out = []
        j.violate = lambda w, det: out.append(w)
        for fec, size in seq:
            j.fec_5min_mean = lambda f=fec: f
            j.ts_stat = lambda s=size: (time.time(), s)
            j.sample()
        return out

    # THE REAL 08:22 TRACE: band down (nothing owed, TS flat), then the
    # radio returns and FEC is instantly 100% while the TS has not yet
    # resumed. The old code violated here.
    trace = [(None, 1000)] * 6 + [(100.0, 1000)] * 3
    check("NEGATIVE CONTROL: band returns after a yield, TS not yet "
          "resumed -> NO violation (the 08:22 false positive)",
          "ts-growth" not in run(trace), str(run(trace)))

    # radio held by someone else, band irrelevant -> never owed
    check("NEGATIVE CONTROL: radio yielded -> NO violation",
          "ts-growth" not in run([(100.0, 1000)] * 8,
                                 yielded="sonde_rx (priority 100)"))

    # A GENUINE STALL: band up and owed throughout, TS flat the whole time
    check("a REAL stall (band up and owed, TS flat) STILL violates",
          "ts-growth" in run([(100.0, 1000)] * 8), "must fire")

    # healthy: band up, TS growing
    grow = [(100.0, 1000 + 500 * i) for i in range(8)]
    check("healthy (band up, TS growing) -> no violation",
          "ts-growth" not in run(grow))


def leg_gain_step():
    print("leg 19: overload is not a fade -- the referee decides")
    d = scratch_dir("gainstep")

    # the gain ladder must actually move, and never to where it already is
    e = "--assets all --accel cpu --rfgain 4"
    check("current_rfgain reads the configured value",
          RUN.current_rfgain(e) == 4, str(RUN.current_rfgain(e)))
    check("swap_rfgain replaces in place (no duplicate flags)",
          RUN.swap_rfgain(e, 6) == "--assets all --accel cpu --rfgain 6",
          RUN.swap_rfgain(e, 6))
    check("swap_rfgain appends when the flag is absent",
          RUN.swap_rfgain("--assets all", 6).endswith("--rfgain 6"))
    cur = 4
    order = [g for g in RUN.GAIN_LADDER if g != cur]
    check("the ladder never re-tries the gain that is already failing",
          cur not in order, str(order))
    check("the ladder is bounded and cycles",
          len(order) == len(RUN.GAIN_LADDER) - 1 and
          [order[i % len(order)] for i in range(RUN.GAIN_MAX_STEPS)]
          == order[:RUN.GAIN_MAX_STEPS], str(order))

    # NEGATIVE CONTROL: a REAL fade -- the referee cannot lock either, so
    # hdhr_snq() returns None and the gain must NOT be touched.
    def would_step(snq):
        return snq is not None and snq >= RUN.HDHR_SNQ_OK
    check("NEGATIVE CONTROL: real fade (referee has no lock) -> NO gain "
          "step, we wait as before", not would_step(None))
    check("NEGATIVE CONTROL: referee also struggling (snq 40) -> NO gain "
          "step", not would_step(40))
    check("referee healthy (snq 100) while we fail -> STEP the gain",
          would_step(100))

    # the stampede bound
    gs = [time.time()] * RUN.GAIN_MAX_STEPS
    gs = [t for t in gs if time.time() - t < RUN.GAIN_WIN]
    check("gain-step budget is enforced (no front-end thrash)",
          not (len(gs) < RUN.GAIN_MAX_STEPS),
          f"{len(gs)}/{RUN.GAIN_MAX_STEPS} used")


def leg_young_generation():
    print("leg 18: a NEWBORN lane owes no window yet")
    sys.path.insert(0, os.path.join(ROOT, "tools"))
    import atsc3_warden as W                                    # noqa: E402
    d = scratch_dir("young")
    ts = os.path.join(d, "_tv", "live_tv2.ts")
    open(ts, "wb").write(b"\x47" * 1000)
    clock = FakeClock()
    real = W.time
    W.time = clock
    try:
        h = HarnessWarden(W, d, clock)
        vt0 = clock.t - 3600                     # long-established viewer
        os.utime(ts, (clock.t, clock.t))
        h.set_roles(d, vt0, with_vlc=False)      # NO player
        seq = 225833804
        # a JUST-ROLLED lane: only ~40 s of media, far below lag+lead
        for _ in range(25):                      # 25 min of wall time
            seq += 15
            write_lanes(d, last_seq=seq, media_slots=20)   # ~40 s of media
            os.utime(ts, (clock.t, clock.t))
            h.w.loop()
            clock.sleep(60)
        check("newborn generation + no player: viewer NOT killed "
              "(the 20:26:10 false kill cannot recur)",
              not h.kills, f"kills {[k[1][:45] for k in h.kills]}")

        # once the generation MATURES past lag+lead, the rule must engage
        h.kills.clear(); h.spawns.clear()
        h.w.obs.pop("no_vlc", None)
        h.w.state["restarts"] = {}
        fired = False
        for _ in range(20):
            seq += 30
            write_lanes(d, last_seq=seq, media_slots=400)  # ~800 s of media
            os.utime(ts, (clock.t, clock.t))
            h.w.loop()
            if h.kills:
                fired = True
                break
            clock.sleep(60)
        check("mature generation + still no player -> the rule DOES engage",
              fired, f"kills {[k[1][:45] for k in h.kills]}")
    finally:
        W.time = real


def leg_lane_stops_midtimer():
    print("leg 17: a lane that STOPS mid-timer must not be called a wedge")
    sys.path.insert(0, os.path.join(ROOT, "tools"))
    import atsc3_warden as W                                    # noqa: E402
    d = scratch_dir("midtimer")
    ts = os.path.join(d, "_tv", "live_tv2.ts")
    open(ts, "wb").write(b"\x47" * 1000)
    clock = FakeClock()
    real = W.time
    W.time = clock
    try:
        h = HarnessWarden(W, d, clock)
        seq = 225830000
        write_lanes(d, last_seq=seq)
        os.utime(ts, (clock.t, clock.t))
        wd = h.w

        # phase 1: lane advancing, TS moving -> healthy, timer never arms
        for _ in range(6):
            seq += 30
            write_lanes(d, last_seq=seq)
            os.utime(ts, (clock.t, clock.t))
            check_silent = wd.viewer_wedged()
            clock.sleep(30)
        check("healthy (lane advancing, TS moving): no wedge",
              not check_silent)

        # phase 2: THE REAL EVENT -- the lane advances a little more, then
        # STOPS DEAD (sonde takes the radio) and the TS correctly freezes.
        for _ in range(3):
            seq += 30
            write_lanes(d, last_seq=seq)          # TS deliberately NOT touched
            wd.viewer_wedged()
            clock.sleep(30)
        frozen_at = seq
        fired = False
        for _ in range(20):                       # 10 minutes of a dead lane
            write_lanes(d, last_seq=frozen_at)    # lane STOPPED
            if wd.viewer_wedged():
                fired = True
                break
            clock.sleep(30)
        check("NEGATIVE CONTROL: lane stops mid-timer -> NO wedge "
              "(the 18:43 false kill cannot recur)", not fired)

        # phase 3: a REAL wedge must still fire -- lane keeps advancing
        # continuously while the TS stays stale
        wd.obs["viewer"] = None
        fired_at = None
        t0 = clock.t
        for _ in range(20):
            seq += 30
            write_lanes(d, last_seq=seq)          # advancing EVERY sample
            if wd.viewer_wedged():
                fired_at = clock.t - t0
                break
            clock.sleep(30)
        check("a REAL wedge (lane advancing every sample, TS stale) still "
              "fires", fired_at is not None,
              f"after {fired_at}s" if fired_at else "never fired")
        check("and it waits the full 180 s first",
              fired_at is not None and fired_at >= 180, f"{fired_at}s")
    finally:
        W.time = real


def leg_starve(bypass):
    n = 3 if bypass else 4
    print(f"leg {n}: dead-worker starvation, bypass "
          f"{'ON' if bypass else 'OFF (negative control)'}")
    d = scratch_dir(f"starve_{'on' if bypass else 'off'}")
    sim = LaneSim(d, with_audio=True)
    make_dead_wav(d, sim)
    sim.feed(10)
    extra = () if bypass else ("--no-audio-starve-bypass",)
    pl = spawn_viewer(d, extra)
    try:
        t0 = time.time()
        # feed 1 fragment/s. The bypass arms only after the wav frontier
        # has been flat for max(hold, 90) s, so the measurement window
        # opens at t=120 -- both branches hold identically before that.
        sz_mark = None
        t_mark = None
        while time.time() - t0 < 210:
            sim.feed(1)
            time.sleep(1.0)
            if time.time() - t0 > 120 and sz_mark is None:
                sz_mark = ts_size(d)
                t_mark = time.time()
        grew = ts_size(d) - (sz_mark or 0)
        dt = time.time() - (t_mark or t0)
        rate = grew / max(dt, 1.0)        # bytes/s of sink growth
        if bypass:
            check("TS grows steadily with a dead wav (no metronome)",
                  rate > 100_000, f"{rate/1e3:.0f} kB/s over {dt:.0f}s")
        else:
            check("negative control: old hold path starves the sink",
                  rate < 60_000, f"{rate/1e3:.0f} kB/s over {dt:.0f}s")
    finally:
        pl.terminate()
        try:
            pl.wait(timeout=10)
        except subprocess.TimeoutExpired:
            pl.kill()


# ---------------------------------------------------------------- leg 5

def leg_e42():
    print("leg 5: E42 -- fade + two generation rolls + recovery")
    d = scratch_dir("e42")
    sim = LaneSim(d, with_audio=False)
    sim.feed(12)
    pl = spawn_viewer(d)
    try:
        for _ in range(25):                       # healthy stretch
            sim.feed(1); time.sleep(1.0)
        sz_pre = ts_size(d)
        check("TS grew before the fade", sz_pre > 1_000_000,
              f"{sz_pre/1e6:.1f} MB")
        time.sleep(30)                            # FADE begins: drain out
        sz_drain = ts_size(d)
        time.sleep(70)                            # deep in the fade
        sz_fade = ts_size(d)
        check("detector: TS did NOT grow across the fade (the gate can "
              "fail)", sz_fade == sz_drain,
              f"+{(sz_fade-sz_drain)/1e3:.0f} kB in 70 s idle")
        sim.roll()                                # multi-restart recovery
        sim.roll()
        t_rec = time.time()
        recovered = False
        while time.time() - t_rec < 90:
            sim.feed(1)
            time.sleep(1.0)
            if ts_size(d) > sz_fade + 500_000:
                recovered = True
                break
        check("TS resumed within 60 s of recovery", recovered
              and time.time() - t_rec <= 62,
              f"{time.time()-t_rec:.0f}s after roll x2")
        logtxt = open(os.path.join(d, "viewer.log"),
                      encoding="utf-8", errors="replace").read()
        check("viewer logged the generation re-anchor",
              "lane generation" in logtxt)
    finally:
        pl.terminate()
        try:
            pl.wait(timeout=10)
        except subprocess.TimeoutExpired:
            pl.kill()


# ---------------------------------------------------------------- leg 6

def leg_player_vanish():
    print("leg 6: player process vanishes -> viewer respawns it")
    d = scratch_dir("vanish")
    sim = LaneSim(d, with_audio=False)
    sim.feed(12)
    lf = open(os.path.join(d, "viewer.log"), "ab", buffering=0)
    pl = subprocess.Popen(
        [PY, "-u", os.path.join(ROOT, "tools", "atsc3_tv.py"),
         "--live-dir", d, "--lag", "8", "--lead", "10",
         "--audio-hold", "0", "--vlc-headless"],
        cwd=ROOT, stdout=lf, stderr=subprocess.STDOUT, env=env_pinned())
    pidf = os.path.join(d, "_tv", "player.pid")
    try:
        vpid = None
        t0 = time.time()
        while time.time() - t0 < 90:
            sim.feed(1); time.sleep(1.0)
            try:
                vpid = json.load(open(pidf)).get("pid")
            except (OSError, ValueError):
                vpid = None
            if vpid and psutil and psutil.pid_exists(vpid):
                break
        check("headless VLC came up", vpid is not None
              and psutil.pid_exists(vpid), f"pid {vpid}")
        if vpid:
            subprocess.run(["taskkill", "/PID", str(vpid), "/F"],
                           capture_output=True)
        t_kill = time.time()
        new_pid = None
        while time.time() - t_kill < 75:
            sim.feed(1); time.sleep(1.0)
            try:
                p2 = json.load(open(pidf)).get("pid")
            except (OSError, ValueError):
                p2 = None
            if p2 and p2 != vpid and psutil.pid_exists(p2):
                new_pid = p2
                break
        check("viewer respawned the player within 75 s",
              new_pid is not None, f"{new_pid} after "
              f"{time.time()-t_kill:.0f}s")
        check("viewer never exited", pl.poll() is None)
    finally:
        pl.terminate()
        try:
            pl.wait(timeout=10)
        except subprocess.TimeoutExpired:
            pl.kill()
        try:
            rec = json.load(open(pidf))
            if psutil.pid_exists(rec.get("pid", -1)):
                subprocess.run(["taskkill", "/PID", str(rec["pid"]), "/F"],
                               capture_output=True)
        except (OSError, ValueError):
            pass


# ------------------------------------------------------------ live legs

def census_role(role):
    out = []
    for p in psutil.process_iter(["pid", "name", "cmdline"]):
        try:
            nm = (p.info["name"] or "").lower()
            c = " ".join(p.info["cmdline"] or [])
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            continue
        if role == "vlc":
            if nm.startswith("vlc") and "live_tv2.ts" in c:
                out.append(p.info["pid"])
            continue
        if not nm.startswith("python"):
            continue
        if role == "run" and "atsc3_run.py" in c:
            out.append(p.info["pid"])
        elif role == "eng" and "atsc3_audio.py" in c and "--pid 14" not in c:
            out.append(p.info["pid"])
        elif role == "spa" and "atsc3_audio.py" in c and "--pid 14" in c:
            out.append(p.info["pid"])
        elif role == "subs" and "atsc3_subs.py" in c:
            out.append(p.info["pid"])
        elif role == "viewer" and re.search(r"atsc3_tv\.py", c):
            out.append(p.info["pid"])
        elif role == "chain" and "-m atsc3" in c and "watch" in c:
            out.append(p.info["pid"])
        elif role == "warden" and "atsc3_warden.py" in c:
            out.append(p.info["pid"])
        elif role == "forks" and "multiprocessing" in c \
                and "parent_pid=" in c:
            m = re.search(r"parent_pid=(\d+)", c)
            if m and not psutil.pid_exists(int(m.group(1))):
                out.append(p.info["pid"])
    return out


def wait_role(role, want_pid_not=None, budget=120):
    t0 = time.time()
    while time.time() - t0 < budget:
        pids = census_role(role)
        if pids and (want_pid_not is None or pids[0] != want_pid_not):
            return time.time() - t0, pids
        time.sleep(5)
    return None, census_role(role)


def kill_pid(pid):
    subprocess.run(["taskkill", "/PID", str(pid), "/T", "/F"],
                   capture_output=True)


def leg_live_player():
    """The one leg that needs LIVE AIR: kill the window, get it back.

    Split out because it is the only outstanding SKIP from the band-down
    run, and re-running the whole injection suite would cost several
    minutes of a scarce delivering window (and roll the lanes)."""
    print("leg 7a (live air): player restore")
    if not band_up():
        check("band is up (a window is owed)", False, "FEC < 60%")
        return
    vp = census_role("vlc")
    if not vp:
        check("player present to inject", False)
        return
    v0 = census_role("viewer")
    kill_pid(vp[0])
    dt, pids = wait_role("vlc", want_pid_not=vp[0], budget=90)
    check("player restored <= 60 s on live air",
          dt is not None and dt <= 62,
          f"{dt if dt is None else round(dt)}s  {vp[0]} -> {pids}")
    check("exactly one player after restore", len(census_role("vlc")) == 1,
          str(census_role("vlc")))
    check("the viewer was NOT restarted (it owns its own player)",
          census_role("viewer") == v0, f"{v0} -> {census_role('viewer')}")


def leg_live():
    print("legs 7-8: LIVE fault injection (deployed stack + warden)")
    ok = True
    for role in ("warden", "run", "eng", "viewer", "chain"):
        pids = census_role(role)
        if len(pids) != 1:
            print(f"  precondition FAILED: {role} x{len(pids)}")
            ok = False
    if not ok:
        check("live preconditions (one of each incl. warden)", False)
        return
    check("live preconditions (one of each incl. warden)", True)

    # 7a player -- only meaningful when there is something to play. With
    # the band dead the viewer CORRECTLY has no window (E62b), so this leg
    # is SKIPPED, never passed: it must be re-run on live air before the
    # gate counts. A leg that quietly passes itself when it cannot run is
    # the vacuity bug this campaign keeps finding.
    vp = census_role("vlc")
    if vp:
        kill_pid(vp[0])
        dt, pids = wait_role("vlc", want_pid_not=vp[0], budget=75)
        check("player restored <= 60 s", dt is not None and dt <= 62,
              f"{dt if dt is None else round(dt)}s -> {pids}")
    elif band_up():
        check("player present to inject (band is up, so a window is owed)",
              False)
    else:
        print("  SKIP  player leg -- band down (FEC < 60%), no window is "
              "owed; MUST be re-run on live air")
        SKIPPED.append("player restore")

    # 7b viewer
    v = census_role("viewer")[0]
    kill_pid(v)
    dt, pids = wait_role("viewer", want_pid_not=v, budget=120)
    check("viewer restored <= 90 s", dt is not None and dt <= 92,
          f"{dt if dt is None else round(dt)}s")

    # 7c eng worker
    e = census_role("eng")
    if e:
        kill_pid(e[0])
        dt, pids = wait_role("eng", want_pid_not=e[0], budget=150)
        check("eng worker restored <= 120 s", dt is not None and dt <= 122,
              f"{dt if dt is None else round(dt)}s")
    else:
        check("eng worker present to inject", False)

    # 7d chain child (atsc3_run's job, warden sweeps the orphans)
    c = census_role("chain")
    if c:
        kill_pid(c[0])
        dt, pids = wait_role("chain", want_pid_not=c[0], budget=120)
        check("chain child restored by atsc3_run", dt is not None,
              f"{dt if dt is None else round(dt)}s")
        time.sleep(45)
        check("no orphaned forks after chain kill",
              len(census_role("forks")) == 0)
    else:
        check("chain child present to inject", False)

    # one-each snapshot after the dust settles
    time.sleep(20)
    # EVERY role, not the four I happened to inject: `spa` and `subs` had
    # no census branch at all, so the invariant was blind to them (a live
    # spa worker read as 0 and nobody noticed).
    dups = {r: census_role(r) for r in
            ("warden", "run", "eng", "spa", "subs", "viewer", "vlc")}
    check("one-each invariant after injections (all roles)",
          all(len(v) <= 1 for v in dups.values()),
          str({k: len(v) for k, v in dups.items()}))
    missing = [r for r in ("warden", "run", "eng", "spa", "subs", "viewer")
               if not dups[r]]
    check("no role went missing across the injections",
          not missing, f"missing: {missing}" if missing else "all present")

    # 8 negative control: paused warden does not restore
    live_dir = os.path.join(ROOT, "data", "e31")
    pause = os.path.join(live_dir, "_warden", "pause")
    open(pause, "w").write("gate_e61 negative control\n")
    try:
        e = census_role("eng")
        if e:
            kill_pid(e[0])
            time.sleep(130)
            check("negative control: eng STAYS dead under paused warden",
                  len(census_role("eng")) == 0)
        else:
            check("eng present for negative control", False)
    finally:
        # never delete under data/ -- rename with a unique suffix
        os.rename(pause, pause + f".off.{RUN_TAG}")
    dt, pids = wait_role("eng", budget=150)
    check("eng restored after unpause <= 120 s",
          dt is not None and dt <= 125,
          f"{dt if dt is None else round(dt)}s")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--live", action="store_true")
    ap.add_argument("--live-player", action="store_true",
                    help="just the live-air player-restore leg")
    ap.add_argument("--only", type=int, default=None)
    a = ap.parse_args()
    os.makedirs(SCRATCH, exist_ok=True)
    if a.live_player:
        leg_live_player()
    elif a.live:
        leg_live()
    else:
        legs = {1: leg_wav_wall, 2: leg_tree_kill,
                3: lambda: leg_starve(True), 4: lambda: leg_starve(False),
                5: leg_e42, 6: leg_player_vanish, 9: leg_radio_yield,
                10: leg_e62b, 11: leg_judge_clock, 12: leg_referee,
                13: leg_fade_not_wedge, 14: leg_rf_referee,
                15: leg_maintenance_window, 16: leg_av_skew,
                17: leg_lane_stops_midtimer, 18: leg_young_generation, 19: leg_gain_step, 20: leg_ts_growth_windows}
        for k, fn in legs.items():
            if a.only and k != a.only:
                continue
            fn()
    if SKIPPED:
        print(f"\nSKIPPED (must be re-run when the band is up): "
              f"{', '.join(SKIPPED)}")
    print(f"\n{'ALL PASS' if not FAILS else 'FAILURES: ' + ', '.join(FAILS)}"
          + (" -- INCOMPLETE, skipped legs outstanding" if SKIPPED else ""))
    return 0 if not FAILS else 1


if __name__ == "__main__":
    sys.exit(main())
