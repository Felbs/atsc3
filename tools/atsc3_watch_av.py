#!/usr/bin/env python3
"""atsc3_watch_av.py -- one command for live ATSC 3.0 TV with sound (+CC).

This is a thin ORCHESTRATOR. It does not mux anything itself -- an earlier
version did, hand-rolling a raw-PCM-over-TCP bridge into ffmpeg, and it drifted
(no clock servo) and floated the captions in a separate window. The mature
viewer `atsc3_tv.py` already solves all of that: per-chunk MPEG-TS mux with
per-input -itsoffset alignment, the mp2 frame-grid carry (E48/E96), and
captions as a soft DVB-subtitle track the player renders over the UNTOUCHED
broadcast HEVC (E97 -- no re-encode, press 't' to toggle them) -- it just
needs the lanes and the decoded audio/caption sidecars fed to it.

So this wires the standard pipeline and cleans up when the window closes:

    chain            atsc3_run --assets all      -> video + audio + caption lanes
    audio worker     atsc3_audio                 -> live_audio.wav (our AC-4 dec)
    caption worker   atsc3_subs (if --cc)        -> live.srt
    viewer           atsc3_tv --mode v2 ffplay   -> mux (sync!) + HEVC copy + soft CC
                                                   + telemetry (_tv/telemetry.jsonl)

Closing the ffplay window stops everything (--exit-on-player-close).

Radio discipline: the chain (atsc3_run) owns the single-tenant SDR; the other
workers only read files.

Usage:
    python tools/atsc3_watch_av.py --rf 33 --ant "Antenna B"
    python tools/atsc3_watch_av.py --rf 33 --ant "Antenna B" --cc
    python tools/atsc3_watch_av.py --rf 33 --ant "Antenna B" --lang spa
"""
from __future__ import annotations

import argparse
import glob
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# audio-track map for this mux: English 5.1 main vs Spanish SAP stereo pair
LANG = {"eng": (13, "5_X"), "spa": (14, "pair")}


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] av: {m}", flush=True)


def find_python(explicit=None):
    cands = [explicit] if explicit else []
    cands += [os.path.expanduser(r"~\radioconda\python.exe"),
              os.path.expanduser("~/radioconda/python.exe"),
              os.path.expanduser("~/radioconda/bin/python"),
              # Linux distro installs (Arch/Debian python3-soapysdr): the
              # running interpreter and the system python3 are candidates too
              sys.executable,
              shutil.which("python3")]
    for c in cands:
        if c and os.path.isfile(c):
            try:
                subprocess.run([c, "-c", "import SoapySDR"], check=True,
                               capture_output=True, timeout=30)
                return c
            except Exception:                                  # noqa: BLE001
                continue
    return None


def wait_for(path, min_bytes, timeout):
    t0 = time.time()
    while time.time() - t0 < timeout:
        try:
            if os.path.getsize(path) >= min_bytes:
                return True
        except OSError:
            pass
        time.sleep(0.3)
    return False


def teardown(procs):
    """Stop every worker we spawned, children included, on both platforms.

    E99 (8/22, Ubuntu): this used taskkill unconditionally; on Linux the
    command does not exist, the exception was swallowed, and the chain +
    both audio workers + the caption worker lived on after the window
    closed -- holding the single-tenant radio. Windows: taskkill /T.
    Linux: terminate the process and its children (psutil when present,
    else the process itself), then SIGKILL what refuses."""
    for p in reversed(procs):
        try:
            if p.poll() is not None:
                continue
            if sys.platform.startswith("win"):
                subprocess.run(["taskkill", "/PID", str(p.pid), "/T", "/F"],
                               capture_output=True)
                continue
            kids = []
            try:
                import psutil
                kids = psutil.Process(p.pid).children(recursive=True)
            except Exception:                                  # noqa: BLE001
                pass
            for k in kids:
                try:
                    k.terminate()
                except Exception:                              # noqa: BLE001
                    pass
            p.terminate()
            try:
                p.wait(timeout=5)
            except Exception:                                  # noqa: BLE001
                p.kill()
            for k in kids:
                try:
                    if k.is_running():
                        k.kill()
                except Exception:                              # noqa: BLE001
                    pass
        except Exception:                                      # noqa: BLE001
            pass


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rf", type=int, required=True)
    ap.add_argument("--ant", default="Antenna B")
    ap.add_argument("--lang", choices=("eng", "spa"), default="eng",
                    help="which audio programme to decode and open "
                         "on: eng = pid13 5.1 main, spa = pid14 SAP. With "
                         "--both-langs the other one is decoded too and 'a' "
                         "switches between them.")
    ap.add_argument("--cc", action="store_true",
                    help="closed captions as a soft track ('t' toggles)")
    ap.add_argument("--both-langs", dest="both_langs", default=None,
                    action="store_true",
                    help="decode BOTH audio programmes so 'a' switches between "
                         "them. Default: on a box with >=8 logical cores. "
                         "Below that only the --lang programme is decoded, "
                         "because the second worker costs a core the decoder "
                         "needs (E110, measured on a Raspberry Pi 5).")
    ap.add_argument("--one-lang", dest="both_langs", action="store_false",
                    help="decode only the --lang programme")
    ap.add_argument("--stereo", action="store_true",
                    help="decode the main as L/R-only stereo instead of "
                         "full 5.1 (E98: 5.1 runs 1.87x realtime here; use "
                         "this on a box where it cannot keep up)")
    ap.add_argument("--live-dir", default=None)
    ap.add_argument("--tv-extra", default="",
                    help="extra args appended to the atsc3_tv viewer command")
    ap.add_argument("--chain-extra", default="",
                    help="extra args appended to the chain's watch command "
                         "(e.g. \"--decode-procs 2 --threads 2\") -- the "
                         "tuning experiments' hook")
    ap.add_argument("--python", default=None)
    a = ap.parse_args()

    py = find_python(a.python)
    if py is None:
        log("FATAL: no interpreter can import SoapySDR (need radioconda)")
        return 2

    live = a.live_dir or os.path.join(ROOT, "lab", "av_live")
    os.makedirs(live, exist_ok=True)
    # Start clean. A stale lane / wav / resume-state from a PREVIOUS session
    # bleeds through: the audio decoder reads leftover content and you hear a
    # different programme than the picture until fresh data rolls in (the
    # "wrong audio for the first minute" symptom). Remove this session's known
    # live artifacts by glob so every worker tails from byte 0 and shares one
    # origin. Targeted globs only -- never rm -rf a dir that may hold captures.
    for pat in ("live_video_*", "live_audio*", "live_subs_*", "live.json",
                "live.srt", "*.state.json", "*.idx"):
        for f in glob.glob(os.path.join(live, pat)):
            try:
                os.remove(f)
            except OSError:
                pass
    tvdir = os.path.join(live, "_tv")
    if os.path.isdir(tvdir):
        shutil.rmtree(tvdir, ignore_errors=True)

    procs = []

    def spawn(cmd, quiet=True):
        p = subprocess.Popen(
            cmd, cwd=ROOT,
            stdout=subprocess.DEVNULL if quiet else None,
            stderr=subprocess.STDOUT if quiet else None)
        p._av_cmd = list(cmd)          # E111: so the supervisor can respawn
        procs.append(p)
        return p

    try:
        # 1. chain: video + audio + caption lanes; owns the radio
        # E115 (8/23, measured under the FULL stack, 8 min per leg):
        # decoder-alone the E104 default (1 proc x 4 threads) wins, but with
        # the audio worker, mux and player also on the box the contention
        # regime flips: 1x4 0.93x / 2x3 0.95x / 2x2 0.97x sustained. The
        # orchestrator knows it is running the whole stack, so IT carries
        # the small-box override; an explicit --chain-extra still wins.
        # E117 (8/23): at parity (~1.0x sustained) the default 24-deep raw
        # queue sheds whole frames on every transient dip (mux burst, audio
        # pass) -- FEC read a PERFECT 100% because shed frames never reach
        # it, while video+audio lanes lost the SAME slots (a hole every
        # ~40 s, 11% of slots) and the audio worker padded the gaps with
        # silence = the "quieter then louder" report. --raw-queue 600
        # (~74 s cushion) converts dips into latency the >1.0x stretches
        # repay: 16 min measured, ZERO holes past the anchor, pads frozen.
        # The memory only materialises when the queue actually backs up.
        extra = a.chain_extra
        if not extra and (os.cpu_count() or 1) < 8:
            extra = "--decode-procs 2 --threads 2 --raw-queue 600"
            # 2026-09-27: on a small box the E60 margin levers cost ~20 % of the
            # chain on healthy air (see lab/m16_margin.py MarginCfg.auto_off_db);
            # let the chain drop them by itself once the carrier reads 19 dB+.
            os.environ.setdefault("ATSC3_MARGIN_AUTO_OFF", "19")
            # 2026-09-27: software HEVC in mpv is ~half a core on a small box;
            # skip the in-loop deblocking filter, drop late frames at the
            # decoder, and keep lavc to two threads so the chain's workers
            # win the cores. Measured 10-min runs on a 2-core laptop: audio
            # pads 1625 -> 188, re-acquisitions 3 -> 1, no stalls either way.
            os.environ.setdefault(
                "ATSC3_MPV_EXTRA",
                "--vd-lavc-skiploopfilter=all --framedrop=decoder+vo --vd-lavc-threads=2")
        chain = spawn([py, "tools/atsc3_run.py", "--rf", str(a.rf),
                       "--ant", a.ant, "--secs", "0", "--live-dir", live,
                       "--extra", ("--assets all " + extra).strip()])
        log(f"chain up (pid {chain.pid}) -- RF{a.rf} on {a.ant}")

        # 2. audio workers: BOTH languages, each its own process reading
        #    the lane file (E98). The 5.1 main -> live_audio.wav (6 ch,
        #    1.87x realtime on this box; --stereo for L/R only), the SAP
        #    pair -> live_audio_spa.wav. atsc3_tv muxes both; 'a' in the
        #    window switches. A missing SAP lane just yields silence on
        #    track 2 -- it never fails the chunk.
        ncpu = os.cpu_count() or 1
        both = a.both_langs if a.both_langs is not None else ncpu >= 8
        pid_e, el_e = LANG[a.lang]
        other = "spa" if a.lang == "eng" else "eng"
        pid_s, el_s = LANG[other]
        # quiet=False (E116): the worker's pass lines carry n_bad, pads and
        # hf_fail -- the audio-quality telemetry; DEVNULL made them invisible
        # E118 (8/23): --start-behind 60 on EVERY worker spawn. Without it
        # a worker respawned hours into a run (wav roll rc=42, or a crash
        # that lost state.json) anchors at lane fragment 0 and re-decodes
        # the WHOLE lane: gaining only ~0.1x on a 1.0x-growing target it
        # never catches the head, and the new wav hits the 4 GiB wall
        # BEFORE reaching live -- rolling again, forever (observed live:
        # rolls at 14:42 and 15:43, audio ~65 min stale, viewer frozen in
        # a re-anchor loop). At first launch the lane is seconds old, so
        # the flag is a no-op; on a resume (state intact) it is ignored.
        audio = spawn([py, "tools/atsc3_audio.py", "--live-dir", live,
                       "--pid", str(pid_e), "--element", el_e,
                       "--channels", "2" if a.stereo else "6",
                       "--start-behind", "60",
                       "--out", os.path.join(live, "live_audio.wav")],
                      quiet=False)
        log(f"audio worker up (pid {audio.pid}) -- {a.lang} pid {pid_e} "
            f"{'stereo' if a.stereo else '5.1'}")
        if both:
            audio2 = spawn([py, "tools/atsc3_audio.py", "--live-dir", live,
                            "--pid", str(pid_s), "--element", el_s,
                            "--channels", "2",
                            "--start-behind", "60",
                            "--out", os.path.join(live,
                                                  "live_audio_spa.wav")])
            log(f"audio worker up (pid {audio2.pid}) -- {other} pid {pid_s} "
                f"stereo")
        else:
            log(f"second programme ({other}) NOT decoded -- {ncpu} cores; "
                f"--both-langs forces it (E110)")

        # 3. caption worker (optional)
        if a.cc:
            subs = spawn([py, "tools/atsc3_subs.py", "--live-dir", live,
                          "--out", os.path.join(live, "live.srt")])
            log(f"caption worker up (pid {subs.pid})")

        # wait for the chain to actually produce video + audio lanes
        video = os.path.join(live, "live_video_pid12.m4s")
        wav = os.path.join(live, "live_audio.wav")
        log("waiting for video + decoded audio ...")
        if not wait_for(video, 400_000, 120):
            log("FATAL: no video within 120 s -- carrier receivable?")
            return 1
        if not wait_for(wav, 200_000, 120):
            log("FATAL: no decoded audio within 120 s")
            return 1

        # 4. the mature viewer, v2 piped into ffplay (E97): the broadcast
        #    HEVC is COPIED (v1 re-encoded it to burn captions -- a real
        #    fidelity loss), captions ride as a soft DVB track the player
        #    renders, a feed thread keeps a stalled player from wedging the
        #    muxer, and every chunk is logged to _tv/telemetry.jsonl. The
        #    Both languages are decoded (E98): eng is stream a:0 (5.1
        #    AC-3), spa is a:1 (stereo mp2); --ffplay-audio picks the one
        #    the window opens on and 'a' cycles.
        #    It waits on the audio frontier itself (--audio-hold), and
        #    closing the window stops it (--exit-on-player-close).
        # E113: prefer mpv on Linux -- it can use libavcodec hwaccels
        # (ffplay cannot), and on a Pi 5 the stateless V4L2 HEVC hwaccel
        # turns the player from ~30% of a core into ~4%. --hwdec=auto falls
        # back to software wherever no hwaccel serves. Windows keeps ffplay.
        player = ("mpv" if not sys.platform.startswith("win")
                  and shutil.which("mpv") else "ffplay")
        tv = [py, "tools/atsc3_tv.py", "--live-dir", live,
              "--mode", "v2", "--player", player,
              "--subs", "soft" if a.cc else "none",
              "--exit-on-player-close"] + (["--stereo"] if a.stereo else [])             + (a.tv_extra.split() if a.tv_extra else [])
        log(f"starting viewer (atsc3_tv v2/{player}: HEVC copy + soft CC) -- "
            f"close the window to stop")
        tvp = spawn(tv, quiet=False)
        # E111 -- SUPERVISE THE WORKERS (found the hard way: the video played
        # all night and the sound stopped at 05:19). The audio worker rolls
        # its wav before the RIFF 4 GiB wall and exits rc=42, which MEANS
        # "respawn me" -- atsc3_audio has said so since E61 -- and nothing in
        # this orchestrator listened. After ~5.5 h of stereo (sooner in 5.1)
        # the worker exited by design and the mux degraded to silence, which
        # is exactly what it should do when audio is MISSING, and exactly
        # wrong as a permanent state. rc=42 -> respawn, always. Any other
        # death -> respawn too, but stampede-gated (4/hour) so a crash loop
        # cannot eat the box; past the gate it logs loudly and stays down.
        deaths = []
        while tvp.poll() is None:
            time.sleep(2.0)
            for w in list(procs):
                if w is tvp or w.poll() is None:
                    continue
                rc_w = w.returncode
                cmd_w = getattr(w, "_av_cmd", None)
                procs.remove(w)
                name = (os.path.basename(cmd_w[1]) if cmd_w and len(cmd_w) > 1
                        else "worker")
                if cmd_w is None:
                    log(f"{name} exited rc={rc_w} and cannot be respawned")
                    continue
                if rc_w != 42:
                    now = time.time()
                    deaths[:] = [t for t in deaths if now - t < 3600.0]
                    deaths.append(now)
                    if len(deaths) > 4:
                        log(f"** {name} DIED rc={rc_w} -- {len(deaths)} "
                            f"deaths this hour, NOT respawning (its output "
                            f"degrades to silence/absence)")
                        continue
                    log(f"** {name} died rc={rc_w} -- respawning "
                        f"({len(deaths)}/4 this hour)")
                else:
                    log(f"{name} rolled its wav (rc=42) -- respawning, as "
                        f"the contract says")
                spawn(cmd_w)
        rc = tvp.returncode
        log(f"viewer exited (rc={rc}) -- shutting down")
        return 0
    except KeyboardInterrupt:
        return 0
    finally:
        teardown(procs)


if __name__ == "__main__":
    sys.exit(main())
