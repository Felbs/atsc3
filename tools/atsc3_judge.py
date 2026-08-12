#!/usr/bin/env python3
"""atsc3_judge.py -- the independent soak referee (E61).

The user is DONE being the failure detector: this process replaces their
eyes. It runs beside (never inside) the warden, samples the stack every
60 s, and holds five invariants. Any violation is logged with evidence and
RESETS THE SOAK CLOCK; the soak passes only after SOAK_HOURS of continuous
wall time with zero violations and zero operator interventions.

    (i)   exactly one of each role (run, eng, spa, subs, viewer), with
          absence tolerated only inside the restore budgets and only for
          roles whose lane exists;
    (ii)  a player: vlc on OUR sink, http control answering -- whenever the
          TS is actually flowing and the viewer is past its banking window;
    (iii) the TS GREW over the last 5 min whenever the chain's FEC 5-min
          mean says the air was decodable (>60%) -- fades are physics and
          are exempt by construction;
    (iv)  no role restarted more than 4x/hour (warden_state.json is the
          ledger);
    (v)   hourly A/V referee: the newest TS audio must correlate > 0.95
          against the worker's wav (FFT normalized cross-correlation on
          8 kHz mono), skipped hon estly when the segment is silence or the
          wav is mid-roll, failed only on two consecutive valid misses.

    python tools/atsc3_judge.py --live-dir data/e31
"""
from __future__ import annotations

import argparse
import base64
import json
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.request
import wave

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

SOAK_HOURS = 6.0
SAMPLE_S = 60.0
FEC_GATE = 60.0                   # 5-min mean 'now' % above this = decodable
BUDGET = {"run": 90.0, "eng": 180.0, "spa": 180.0, "subs": 180.0,
          "viewer": 120.0}        # absence tolerated this long (restore
#                                   budgets + one sample of slack)
PLAYER_BUDGET = 300.0
REFEREE_PERIOD = 3600.0
CORR_MIN = 0.95

# ---- the INDEPENDENT RF referee (E65). An HDHomeRun FLEX 4K sits on the
# same antenna via the splitter and answers over the LAN, so "the band is
# down" stops being an assumption we grant ourselves and becomes a
# measurement someone else took. This closes the campaign's recurring
# vacuity trap from the other side: every band-excused check until now
# trusted OUR chain's own FEC to decide whether to judge OUR chain.
#   both low   -> the air. Excused, as before.
#   ours low, referee HIGH -> it is US. That must violate.
#   referee unreachable -> skip. A missing instrument is never a verdict.
HDHR_BIN = r"C:\Program Files\Silicondust\HDHomeRun\hdhomerun_config.exe"
HDHR_ID = os.environ.get("ATSC3_HDHR_ID", "")  # no default; discover if unset
HDHR_TUNER = os.environ.get("ATSC3_HDHR_TUNER", "1")   # reserved referee


def hdhr_channel_for(rf):
    """The referee must watch THE CARRIER UNDER TEST, not a constant.

    E91, and this is a trap I built myself: E90 gave the referee the power
    to RECLAIM an idle tuner, and at the same moment pinned its channel to
    a hardcoded 587 MHz. Passive-and-wrong is a bad reading; ACTIVE-and-
    wrong retunes the instrument to the wrong carrier and then reports a
    confident snq from it. Run against Fox, the old constant would have
    tuned the referee to RF33 and published a healthy 100 as if it were a
    verdict on Fox -- a FABRICATED AGREEMENT, which is worse than the
    silent zero E90 set out to fix. Making the referee more autonomous
    raised the cost of every assumption baked into it.

    SCOPE LIMIT, equally important: the HDHomeRun is fed from the splitter
    off Old Faithful (Antenna B). Whenever the chain is on a different
    port, the referee is measuring A DIFFERENT ANTENNA. It can then tell
    you THE AIR MOVED; it can never tell you OUR FEED IS HEALTHY. Do not
    quote it as agreement or disagreement across an antenna change."""
    lo = (174 + (rf - 7) * 6) if rf < 14 else (470 + (rf - 14) * 6)
    return "atsc3:%d" % int((lo + 3.0) * 1e6)


HDHR_PERIOD = 300.0        # s between referee reads (cheap, but not free)
HDHR_SNQ_OK = 80           # referee quality that means "the air is fine"
HDHR_CONFIRM = 2           # consecutive disagreements before violating
PAUSE_MAX_S = 900.0        # a maintenance window older than this is stale
SKEW_PERIOD = 900.0        # video-vs-audio skew check cadence (E71)
SKEW_MAX = 2.0             # seconds of video-vs-audio offset tolerated

try:
    import psutil
except ImportError:
    print("atsc3_judge needs psutil", file=sys.stderr)
    raise


def ffbin(name):
    p = shutil.which(name)
    if p:
        return p
    c = os.path.join(
        os.environ.get("LOCALAPPDATA", "C:"), "Microsoft", "WinGet",
        "Packages", "Gyan.FFmpeg_Microsoft.Winget.Source_8wekyb3d8bbwe",
        "ffmpeg-8.1-full_build", "bin", name + ".exe")
    return c if os.path.exists(c) else name


class Judge:
    def __init__(self, a):
        self.d = a.live_dir
        self.a = a
        self.wdir = os.path.join(self.d, "_warden")
        os.makedirs(self.wdir, exist_ok=True)
        self.logf = open(os.path.join(self.wdir, "judge.log"), "a",
                         buffering=1, encoding="utf-8")
        self.hdhr_ch = (os.environ.get("ATSC3_HDHR_CHANNEL")
                        or hdhr_channel_for(getattr(a, "rf", 33)))
        self.soak_start = time.time()
        self.qualified = 0.0       # BAND-DELIVERING seconds banked (E62b)
        self._last_tick = None
        self.absent_since = {}
        self.dup_since = {}
        self.player_bad_since = None
        self.ts_hist = []          # (wall, size)
        self.last_ref = 0.0
        self.ref_misses = 0
        self.last_rf = 0.0            # independent RF referee (E65)
        self.rf_disagree = 0
        self._mw_logged = 0.0         # maintenance-window log throttle
        self._yield_logged = 0.0      # radio-yield log throttle (E72)
        self.last_skew = 0.0          # video-vs-audio skew check (E71)
        self.skew_bad = 0
        self.violations = 0
        self.passed_logged = False

    def maintenance(self):
        """Is the stack under an authorised maintenance window? (E68)

        The judge honours the WARDEN'S OWN pause marker, because the two
        facts are the same fact: if supervision is suspended, the stack is
        not being asked to look after itself, so judging it would score a
        deliberate act as a failure. During the window every invariant is
        skipped AND the qualified clock stops -- the bar must not be
        payable with time nobody was watching.

        The marker EXPIRES (PAUSE_MAX_S). An unattended maintenance script
        that dies must not be able to buy silence for the rest of the
        night; after the window the judge resumes on its own and will
        report whatever it then finds."""
        p = os.path.join(self.wdir, "pause")
        try:
            age = time.time() - os.path.getmtime(p)
        except OSError:
            return None
        if age > PAUSE_MAX_S:
            return None
        try:
            why = open(p, encoding="utf-8", errors="replace").read().strip()
        except OSError:
            why = ""
        return why or "unspecified"

    def log(self, msg):
        line = f"{time.strftime('%m/%d %H:%M:%S')} {msg}"
        print(line, flush=True)
        self.logf.write(line + "\n")

    def violate(self, what, detail):
        """A violation RESTARTS the acceptance clock -- all of it.

        E63 (8/10): when the clock became QUALIFIED time (E62b), this
        still reset only `soak_start`, which no longer feeds the bar --
        `held` reads self.qualified. So violations were counted, logged,
        and then quietly forgiven: the run could reach '6 h clean' having
        failed repeatedly. That is the same shape as the vacuous player
        check, inverted -- a bar that cannot be failed instead of a check
        that cannot fire. Both clocks reset here, and the message quotes
        the QUALIFIED hours actually lost."""
        self.violations += 1
        held = self.qualified / 3600.0
        self.log(f"VIOLATION [{what}]: {detail}  (soak clock RESET; "
                 f"had held {held:.2f} qualified h)")
        self.soak_start = time.time()
        self.qualified = 0.0
        self.passed_logged = False

    # --------------------------------------------------------- sampling
    def _ours(self, cmd):
        m = re.search(r"--live-dir\s+(\S+)", cmd)
        if not m:
            return False
        p = m.group(1).strip('"')
        if not os.path.isabs(p):
            p = os.path.join(ROOT, p)
        return os.path.normcase(os.path.normpath(p)) == \
            os.path.normcase(os.path.normpath(self.d))

    def census_roles(self):
        me = os.getpid()
        roles = {k: [] for k in
                 ("run", "eng", "spa", "subs", "viewer", "vlc", "warden")}
        sink = os.path.normcase(os.path.normpath(
            os.path.join(self.d, "_tv", "live_tv2.ts")))
        for p in psutil.process_iter(["pid", "name", "cmdline",
                                      "create_time"]):
            try:
                nm = (p.info["name"] or "").lower()
                if p.info["pid"] == me or nm not in (
                        "python.exe", "python", "python3", "vlc.exe", "vlc"):
                    continue
                c = " ".join(p.info["cmdline"] or [])
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                continue
            t0 = p.info["create_time"]
            if nm.startswith("vlc"):
                if sink in os.path.normcase(c):
                    roles["vlc"].append((p.info["pid"], t0))
                continue
            if "atsc3_judge.py" in c:
                continue
            if "atsc3_warden.py" not in c and not self._ours(c):
                continue
            if "atsc3_warden.py" in c:
                roles["warden"].append((p.info["pid"], t0))
            elif "atsc3_run.py" in c:
                roles["run"].append((p.info["pid"], t0))
            elif "atsc3_audio.py" in c:
                roles["spa" if "--pid 14" in c else "eng"].append(
                    (p.info["pid"], t0))
            elif "atsc3_subs.py" in c:
                roles["subs"].append((p.info["pid"], t0))
            elif re.search(r"atsc3_tv\.py", c):
                roles["viewer"].append((p.info["pid"], t0))
        return roles

    def lanes(self):
        try:
            return json.load(open(os.path.join(self.d, "live.json"))
                             ).get("lanes", {})
        except (OSError, ValueError):
            return {}

    def fec_5min_mean(self):
        """Mean of the chain log's instantaneous FEC % over ~5 min, or
        None if too few samples (chain down / just started)."""
        p = os.path.join(self.d, "chain.log")
        try:
            with open(p, "rb") as f:
                f.seek(max(0, os.path.getsize(p) - 200_000))
                txt = f.read().decode("utf-8", "replace")
        except OSError:
            return None
        now = time.localtime()
        now_sec = now.tm_hour * 3600 + now.tm_min * 60 + now.tm_sec
        vals = []
        for m in re.finditer(
                r"\[(\d\d):(\d\d):(\d\d)\].*?\[\s*([\d.]+)% now\]", txt):
            t = int(m.group(1)) * 3600 + int(m.group(2)) * 60 + int(m.group(3))
            dt = (now_sec - t) % 86400
            if dt <= 310:
                vals.append(float(m.group(4)))
        return (sum(vals) / len(vals)) if len(vals) >= 20 else None

    def ts_stat(self):
        p = os.path.join(self.d, "_tv", "live_tv2.ts")
        try:
            st = os.stat(p)
            return st.st_mtime, st.st_size
        except OSError:
            return None, None

    def player_http_ok(self):
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

    def radio_yielded(self):
        """Owner of the radio if it is NOT ours, else None. (E72)

        The RF referee asks 'the air is fine, so why aren't WE decoding?'
        -- which is the right question only when we are actually trying.
        The fleet is single-tenant by priority: sonde_rx (100) outranks
        atsc3_watch (60) and our chain YIELDS, correctly, several times a
        night. During a yield we decode nothing while the HDHomeRun --
        which has its own tuner on the splitter -- still reads snq=100,
        and invariant (vi) would blame the stack for obeying the rule it
        is supposed to obey. Same false-positive shape as the warden's
        18:43 kill: a condition measured without asking whether we were
        even in a position to meet it."""
        for p in (os.environ.get("ATSC3_RADIO_LOCK"),
                  r"Z:\the-receiver-repo\radio.lock.json"):
            if not p or not os.path.exists(p):
                continue
            try:
                d = json.load(open(p))
            except (OSError, ValueError):
                continue
            owner, pid = d.get("owner"), d.get("pid")
            if not owner or owner == "atsc3_watch":
                return None
            try:
                if pid and psutil.pid_exists(int(pid)):
                    return f"{owner} (priority {d.get('priority')})"
            except (TypeError, ValueError):
                pass
        return None

    def fec_recent(self, window=90):
        """FEC mean over the last ~90 s -- 'are we failing NOW?'. (E72b)

        The 5-min mean is the right basis for BANKING qualified time (it
        should be conservative), but the wrong basis for ACCUSING the
        stack. Measured 20:26:15: a genuine 95 s dip left the 5-min mean at
        39% for minutes after the chain had fully recovered to 100%, and
        the RF referee -- seeing the HDHomeRun at snq=100 -- began building
        a case that the fault was ours. It was not; it was already over.
        An accusation must rest on the present tense."""
        p = os.path.join(self.d, "chain.log")
        try:
            with open(p, "rb") as f:
                f.seek(max(0, os.path.getsize(p) - 200_000))
                txt = f.read().decode("utf-8", "replace")
        except OSError:
            return None
        now = time.localtime()
        now_s = now.tm_hour * 3600 + now.tm_min * 60 + now.tm_sec
        vals = []
        for m in re.finditer(
                r"\[(\d\d):(\d\d):(\d\d)\].*?\[\s*([\d.]+)% now\]", txt):
            t = (int(m.group(1)) * 3600 + int(m.group(2)) * 60
                 + int(m.group(3)))
            if (now_s - t) % 86400 <= window:
                vals.append(float(m.group(4)))
        return (sum(vals) / len(vals)) if len(vals) >= 6 else None

    def _hdhr(self, *args):
        try:
            r = subprocess.run([HDHR_BIN, HDHR_ID, *args],
                               capture_output=True, text=True, timeout=10)
        except (OSError, subprocess.SubprocessError):
            return None
        return (r.stdout or "").strip()

    def rf_referee(self):
        """(snq, ss, lock) from the reserved HDHomeRun tuner, or None.

        Returns None -- NOT a zero -- whenever the tuner is not locked.
        E88/E89: an unheld tuner reports lock=none/snq=0, which is
        numerically identical to "RF33 is off the air". Twice tonight an
        operator window released tuner1 and it sat reading snq=0 against
        our 100% FEC, manufacturing a referee DISAGREEMENT out of an idle
        instrument. A control that can silently become a constant is worse
        than no control, because it still looks like a measurement.

        SELF-HEAL, added because hand-retuning after every window is a
        step that will eventually be forgotten at 3 a.m.: if the tuner is
        BOTH unlocked AND parked on no channel, it belongs to nobody, so
        claim it and tune it back. If it is parked on some OTHER channel,
        somebody else is using it -- leave it strictly alone and report
        "not measuring". We only ever touch tuner{HDHR_TUNER}; the
        operator's watchers live on the other tuners."""
        if not os.path.exists(HDHR_BIN):
            return None
        out = self._hdhr("get", f"/tuner{HDHR_TUNER}/status")
        if out is None:
            return None
        m = re.search(r"lock=(\S+)\s+ss=(\d+)\s+snq=(\d+)", out)
        if m and m.group(1) not in ("none", "", "0"):
            return int(m.group(3)), int(m.group(2)), m.group(1)

        ch = self._hdhr("get", f"/tuner{HDHR_TUNER}/channel")
        if ch is not None and ch.lower() in ("none", ""):
            self._hdhr("set", f"/tuner{HDHR_TUNER}/channel", self.hdhr_ch)
            time.sleep(5.0)
            out = self._hdhr("get", f"/tuner{HDHR_TUNER}/status") or ""
            m = re.search(r"lock=(\S+)\s+ss=(\d+)\s+snq=(\d+)", out)
            if m and m.group(1) not in ("none", "", "0"):
                self.log(f"  referee: tuner{HDHR_TUNER} was idle -- "
                         f"retuned to {self.hdhr_ch} (RF{self.a.rf}), "
                         f"lock={m.group(1)} "
                         f"snq={m.group(3)}")
                return int(m.group(3)), int(m.group(2)), m.group(1)
        return None

    def av_skew(self, tail_bytes=8_000_000):
        """Median VIDEO-minus-AUDIO PTS skew in the TS tail, or None.

        E71 exposed a hole in invariant (v): the 'A/V referee' correlates
        the TS's AUDIO against the worker's wav, so it proves the audio is
        continuous and correctly anchored -- and is completely BLIND to the
        video being stamped somewhere else entirely. A muxer bug that put
        video 36.061 s ahead of its own audio sailed past it at corr 0.99,
        because both numbers it compares were audio.

        A check that cannot fail for the fault it appears to cover is the
        exact shape this campaign keeps finding. So: read the packet
        timestamps of both streams over the SAME span and compare them.
        Video and audio interleaved in one TS must advance together; a
        persistent offset is a mux bug no audio-only test can see.

        HONEST ABOUT WHAT IT IS: a DETECTOR, not a meter. Comparing medians
        over a fixed byte tail compresses the magnitude -- gated at
        +6.0s it reads 3.05, at +36.061s it reads 18.1 (about half, and
        monotonic, sign correct). So the number is a reliable alarm and a
        direction, and must not be quoted as the true offset."""
        ts_p = os.path.join(self.d, "_tv", "live_tv2.ts")
        try:
            sz = os.path.getsize(ts_p)
            if sz < tail_bytes * 2:
                return None
        except OSError:
            return None
        tmp = os.path.join(self.wdir, "skew_tail.ts")
        try:
            with open(ts_p, "rb") as f:
                f.seek(sz - tail_bytes)
                b = f.read(tail_bytes)
            k = b.find(b"\x47")
            if k < 0:
                return None
            open(tmp, "wb").write(b[k:])
        except OSError:
            return None
        pts = {}
        for sel, key in (("v:0", "v"), ("a:0", "a")):
            r = subprocess.run(
                [ffbin("ffprobe"), "-v", "error", "-select_streams", sel,
                 "-show_entries", "packet=pts_time", "-of", "csv=p=0", tmp],
                capture_output=True, text=True)
            vals = []
            for line in (r.stdout or "").splitlines():
                line = line.strip().rstrip(",")
                try:
                    vals.append(float(line))
                except ValueError:
                    continue
            if len(vals) < 20:
                return None
            pts[key] = vals
        # compare the OVERLAPPING span only: medians of each stream's
        # timestamps over the same bytes should differ by ~0 in a sane mux
        import statistics
        return statistics.median(pts["v"]) - statistics.median(pts["a"])

    # ------------------------------------------------------ A/V referee
    def referee(self):
        """corr(TS tail audio, wav) -> ('pass'|'fail'|'skip', detail)."""
        ts_p = os.path.join(self.d, "_tv", "live_tv2.ts")
        wav_p = os.path.join(self.d, "live_audio.wav")
        try:
            ts_sz = os.path.getsize(ts_p)
            if ts_sz < 40_000_000:
                return "skip", "TS too small"
            if time.time() - os.path.getmtime(wav_p) > 600:
                return "skip", "wav not advancing (fade or mid-roll)"
        except OSError:
            return "skip", "files missing"
        tmp_ts = os.path.join(self.wdir, "ref_tail.ts")
        tmp_a = os.path.join(self.wdir, "ref_ts.wav")
        tmp_w = os.path.join(self.wdir, "ref_wav.wav")
        with open(ts_p, "rb") as f:
            f.seek(ts_sz - 24_000_000)
            b = f.read(24_000_000)
        k = b.find(b"\x47")
        open(tmp_ts, "wb").write(b[k:])
        r = subprocess.run(
            [ffbin("ffmpeg"), "-y", "-loglevel", "error", "-i", tmp_ts,
             "-map", "0:a:0", "-ac", "1", "-ar", "8000", tmp_a],
            capture_output=True)
        if r.returncode:
            return "skip", f"ts demux failed: {r.stderr[-120:]}"
        # last ~8 min of the worker wav, resampled the same way
        try:
            w = wave.open(wav_p)
            nf, ch, fs = w.getnframes(), w.getnchannels(), w.getframerate()
            take = min(nf, 480 * fs)
            w.setpos(nf - take)
            raw = w.readframes(take)
            w.close()
        except Exception as e:                                 # noqa: BLE001
            return "skip", f"wav read failed: {e}"
        open(tmp_w + ".raw", "wb").write(raw)
        r = subprocess.run(
            [ffbin("ffmpeg"), "-y", "-loglevel", "error", "-f", "s16le",
             "-ar", str(fs), "-ac", str(ch), "-i", tmp_w + ".raw",
             "-ac", "1", "-ar", "8000", tmp_w],
            capture_output=True)
        if r.returncode:
            return "skip", "wav resample failed"
        try:
            wa = wave.open(tmp_a)
            probe_all = np.frombuffer(
                wa.readframes(wa.getnframes()), "<i2").astype("f8")
            wa.close()
            wb = wave.open(tmp_w)
            ref = np.frombuffer(
                wb.readframes(wb.getnframes()), "<i2").astype("f8")
            wb.close()
        except Exception as e:                                 # noqa: BLE001
            return "skip", f"probe read failed: {e}"
        if len(probe_all) < 30 * 8000 or len(ref) < 60 * 8000:
            return "skip", "segments too short"
        # 20 s probe from the middle of the TS tail
        mid = len(probe_all) // 2
        probe = probe_all[mid - 10 * 8000: mid + 10 * 8000]
        rms = np.sqrt(np.mean(probe ** 2)) / 32768.0
        if rms < 10 ** (-45 / 20):
            return "skip", f"probe is silence (rms {20*np.log10(max(rms,1e-9)):.0f} dBFS)"
        # normalized cross-correlation via FFT + cumsum energies
        n = int(2 ** np.ceil(np.log2(len(ref) + len(probe))))
        R = np.fft.rfft(ref, n)
        P = np.fft.rfft(probe, n)
        xc = np.fft.irfft(R * np.conj(P), n)[: len(ref) - len(probe) + 1]
        c2 = np.concatenate([[0.0], np.cumsum(ref ** 2)])
        eng = c2[len(probe):] - c2[: len(ref) - len(probe) + 1]
        denom = np.sqrt(eng * np.sum(probe ** 2)) + 1e-12
        ncc = xc / denom
        corr = float(np.max(ncc))
        lag = int(np.argmax(ncc))
        return ("pass" if corr > CORR_MIN else "fail",
                f"corr {corr:.3f} at lag {lag/8000.0:.1f}s "
                f"(probe rms {20*np.log10(rms):.0f} dBFS)")

    # -------------------------------------------------------------- run
    def sample(self):
        now = time.time()
        roles = self.census_roles()
        lanes = self.lanes()
        have = {"eng": any(l.get("pid") == 13 for l in lanes.values()),
                "spa": any(l.get("pid") == 14 for l in lanes.values()),
                "subs": any(l.get("kind") == "subs" for l in lanes.values()),
                "run": True, "viewer": True}

        # (i) one of each
        for role in ("run", "eng", "spa", "subs", "viewer"):
            n = len(roles[role])
            if n == 1 or not have[role]:
                self.absent_since.pop(role, None)
                self.dup_since.pop(role, None)
            elif n == 0:
                t0 = self.absent_since.setdefault(role, now)
                if now - t0 > BUDGET[role]:
                    self.violate("one-each",
                                 f"{role} absent {now - t0:.0f}s "
                                 f"(> budget {BUDGET[role]:.0f}s)")
                    self.absent_since[role] = now
            else:
                t0 = self.dup_since.setdefault(role, now)
                if now - t0 > 90:
                    self.violate("one-each",
                                 f"{n}x {role} for {now - t0:.0f}s "
                                 f"(pids {[p for p, _ in roles[role]]})")
                    self.dup_since[role] = now

        # warden itself must be up (it is the layer under test)
        if not roles["warden"]:
            t0 = self.absent_since.setdefault("warden", now)
            if now - t0 > 120:
                self.violate("warden", f"warden absent {now - t0:.0f}s")
                self.absent_since["warden"] = now
        else:
            self.absent_since.pop("warden", None)

        # (ii) player, conditioned on a flowing TS + banked viewer
        mt, sz = self.ts_stat()
        ts_flowing = mt is not None and now - mt < 90
        viewer_age = (now - min(t for _, t in roles["viewer"])
                      if roles["viewer"] else 0.0)
        # E62: ts_flowing was the wrong precondition -- CIRCULAR (the TS is
        # the viewer's own output, so a broken viewer excused the check that
        # judges it). On 8/10 that logged "4.83/6 h clean" across five hours
        # with no player at all. The correct precondition is the CHAIN's
        # output: demand a player while the air is decodable. A fade excuses
        # the window (nothing to show); a broken viewer never excuses itself.
        band_up = (self.fec_5min_mean() or 0) > 60
        if viewer_age > 420 and band_up:
            bad = (len(roles["vlc"]) == 0) or (self.player_http_ok() is False)
            if bad:
                t0 = self.player_bad_since or now
                self.player_bad_since = t0
                if now - t0 > PLAYER_BUDGET:
                    self.violate("player",
                                 f"no live player {now - t0:.0f}s "
                                 f"(vlc x{len(roles['vlc'])}, "
                                 f"TS {'flowing' if ts_flowing else 'idle'})")
                    self.player_bad_since = None
            else:
                self.player_bad_since = None
        else:
            self.player_bad_since = None
        if len(roles["vlc"]) > 1:
            t0 = self.dup_since.setdefault("vlc", now)
            if now - t0 > 90:
                self.violate("one-each", f"{len(roles['vlc'])} players on "
                                         f"the sink")
                self.dup_since["vlc"] = now
        else:
            self.dup_since.pop("vlc", None)

        # (iii) TS growth vs FEC
        #
        # E77: the two halves of this test used MISALIGNED WINDOWS. "FEC is
        # good" is measured NOW; "the TS grew" is measured over the last
        # 400 s. After the radio comes back from a yield, FEC is instantly
        # 100% while the TS history is still full of samples from the
        # period when the viewer legitimately had nothing to write -- so
        # the invariant fired on a stack that was behaving perfectly
        # (measured 08:22:54: violation raised three minutes after the
        # judge's own log said "ours FEC n/a", with the TS growing 2.5 MB/s
        # by the time I looked).
        #
        # The history may only contain time during which output was
        # actually OWED. Any sample taken while the band was down, or while
        # the radio belonged to someone else, is discarded rather than held
        # against the viewer.
        fec = self.fec_5min_mean()
        owed = (fec is not None and fec > FEC_GATE
                and self.radio_yielded() is None)
        if not owed:
            self.ts_hist.clear()
        elif sz is not None:
            self.ts_hist.append((now, sz))
        self.ts_hist[:] = [(t, s) for t, s in self.ts_hist if now - t < 400]
        if owed and len(self.ts_hist) >= 5:
            grew = self.ts_hist[-1][1] > self.ts_hist[0][1]
            # a fresh sink (viewer restart) resets size; growth from any
            # baseline counts, shrink means a NEW sink -- also fine
            reset = self.ts_hist[-1][1] < self.ts_hist[0][1]
            if not grew and not reset:
                self.violate("ts-growth",
                             f"FEC 5-min mean {fec:.0f}% but TS flat at "
                             f"{sz/1e6:.1f} MB for 5 min")
                self.ts_hist.clear()

        # (iv) restart rate
        try:
            st = json.load(open(os.path.join(self.wdir, "warden_state.json")))
            for role, ts in (st.get("restarts") or {}).items():
                k = len([t for t in ts if now - t < 3600])
                if k > 4:
                    self.violate("restart-rate", f"{role} restarted {k}x/hr")
        except (OSError, ValueError):
            pass

        # (vi) THE INDEPENDENT RF REFEREE (E65): is a bad hour the air, or us?
        if now - self.last_rf >= HDHR_PERIOD:
            self.last_rf = now
            ref = self.rf_referee()
            if ref is None:
                self.rf_disagree = 0          # no instrument, no verdict
            else:
                snq, ss, lock = ref
                ours = self.fec_5min_mean()
                yielded = self.radio_yielded()
                if yielded:
                    # we are not decoding because we handed the antenna to
                    # a higher-priority service -- obeying the rule, not
                    # failing. Never a violation; never banks time either.
                    self.rf_disagree = 0
                    if time.time() - self._yield_logged > 600:
                        self._yield_logged = time.time()
                        self.log(f"  radio yielded to {yielded} -- RF referee "
                                 f"check suspended (we are not trying)")
                elif (ours is not None and ours < FEC_GATE
                      and (self.fec_recent() or 0) < FEC_GATE
                      and snq >= HDHR_SNQ_OK):
                    self.rf_disagree += 1
                    self.log(f"referee DISAGREES ({self.rf_disagree}/"
                             f"{HDHR_CONFIRM}): our FEC {ours:.0f}% but "
                             f"HDHomeRun snq={snq} ss={ss} lock={lock} on the "
                             f"same feed")
                    if self.rf_disagree >= HDHR_CONFIRM:
                        self.violate("rf-referee",
                                     f"our chain at {ours:.0f}% FEC while the "
                                     f"independent tuner reads snq={snq} "
                                     f"lock={lock} on the SAME antenna -- the "
                                     f"fault is OURS, not the air")
                        self.rf_disagree = 0
                else:
                    self.rf_disagree = 0

        # (vii) VIDEO-vs-AUDIO skew (E71): the fault invariant (v) is blind to
        if now - self.last_skew >= SKEW_PERIOD:
            self.last_skew = now
            sk = self.av_skew()
            if sk is None:
                self.skew_bad = 0
            else:
                self.log(f"  a/v skew: video - audio = {sk:+.2f}s")
                if abs(sk) > SKEW_MAX:
                    self.skew_bad += 1
                    if self.skew_bad >= 2:
                        self.violate("av-skew",
                                     f"video is {sk:+.2f}s from its own audio "
                                     f"in the muxed TS (limit {SKEW_MAX}s) -- "
                                     f"an audio-only referee cannot see this")
                        self.skew_bad = 0
                else:
                    self.skew_bad = 0

        # (v) hourly referee
        if now - self.last_ref >= REFEREE_PERIOD:
            self.last_ref = now
            verdict, detail = self.referee()
            self.log(f"referee: {verdict} -- {detail}")
            if verdict == "fail":
                self.ref_misses += 1
                if self.ref_misses >= 2:
                    self.violate("av-referee", detail)
                    self.ref_misses = 0
            elif verdict == "pass":
                self.ref_misses = 0

    def run(self):
        self.log(f"judge up: {SOAK_HOURS:.0f} h soak bar, sampling every "
                 f"{SAMPLE_S:.0f} s, FEC gate {FEC_GATE:.0f}%")
        last_status = 0.0
        while True:
            t0 = time.time()
            mw = self.maintenance()
            if mw is not None:
                if time.time() - self._mw_logged > 120:
                    self._mw_logged = time.time()
                    self.log(f"MAINTENANCE WINDOW ({mw}) -- invariants "
                             f"suspended and the qualified clock is STOPPED")
                self._last_tick = None      # bank nothing across the window
                # clear the per-invariant timers so the window itself is
                # never charged as an outage when checking resumes
                self.absent_since.clear()
                self.dup_since.clear()
                self.player_bad_since = None
                self.ts_hist.clear()
                self.rf_disagree = 0
                time.sleep(max(5.0, SAMPLE_S - (time.time() - t0)))
                continue
            try:
                self.sample()
            except Exception as e:                             # noqa: BLE001
                self.log(f"SAMPLE ERROR: {type(e).__name__}: {e}")
            # E62: QUALIFIED soak time only. Wall-clock hours spent with a
            # dead band prove nothing about stability -- the 8/10 run banked
            # 4.83 "clean" hours while the air was at 0% and nothing was
            # even attempting to play. The clock advances only while the
            # chain is delivering, so "6 h clean" means six hours of REAL
            # television survived, not six hours of silence tolerated.
            now_s = time.time()
            # dt is CAPPED at two sample periods: if this process is ever
            # suspended (or the box sleeps), the wall gap is not time the
            # stack was observed surviving, and crediting it would bank
            # hours nobody watched.
            dt = 0.0 if self._last_tick is None \
                else min(now_s - self._last_tick, 2 * SAMPLE_S)
            self._last_tick = now_s
            if (self.fec_5min_mean() or 0) > FEC_GATE:
                self.qualified += dt
            held = self.qualified / 3600.0
            if time.time() - last_status >= 600:
                last_status = time.time()
                mt, sz = self.ts_stat()
                fec = self.fec_5min_mean()
                ref = self.rf_referee()
                # NOT "snq=0": an unlocked or absent referee is reported as
                # not measuring, never as a reading of zero (E89).
                rtxt = (f"snq={ref[0]} ss={ref[1]} lock={ref[2]}" if ref
                        else "NOT MEASURING (unlocked/unreachable)")
                otxt = "n/a" if fec is None else f"{fec:.0f}%"
                self.log(f"  air: ours FEC {otxt}   referee {rtxt}")
                self.log(f"status: soak {held:.2f}/{SOAK_HOURS:.0f} h clean, "
                         f"{self.violations} violations total, TS "
                         f"{(sz or 0)/1e6:.0f} MB "
                         f"(mtime {time.time()-mt:.0f}s ago)" if mt else
                         f"status: soak {held:.2f} h, no TS yet")
            if held >= SOAK_HOURS and not self.passed_logged:
                self.passed_logged = True
                self.log(f"SOAK PASS: {held:.2f} h continuous, all "
                         f"invariants held, zero interventions")
            time.sleep(max(5.0, SAMPLE_S - (time.time() - t0)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--live-dir", default="data/e31")
    ap.add_argument("--rf", type=int, required=True,
                    help="carrier under test. The referee is tuned to "
                         "THIS, not to a constant -- see "
                         "hdhr_channel_for() for why that matters.")
    a = ap.parse_args()
    if not os.path.isabs(a.live_dir):
        a.live_dir = os.path.join(ROOT, a.live_dir)
    Judge(a).run()
    return 0


if __name__ == "__main__":
    sys.exit(main())
