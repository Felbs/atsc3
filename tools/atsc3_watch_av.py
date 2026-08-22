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
              os.path.expanduser("~/radioconda/bin/python")]
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
    for p in reversed(procs):
        try:
            if p.poll() is None:
                subprocess.run(["taskkill", "/PID", str(p.pid), "/T", "/F"],
                               capture_output=True)
        except Exception:                                      # noqa: BLE001
            pass


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rf", type=int, required=True)
    ap.add_argument("--ant", default="Antenna B")
    ap.add_argument("--lang", choices=("eng", "spa"), default="eng",
                    help="audio track selected when the window opens "
                         "(BOTH are decoded; press 'a' in the window to "
                         "switch): eng = pid13 5.1 main, spa = pid14 SAP")
    ap.add_argument("--cc", action="store_true",
                    help="closed captions as a soft track ('t' toggles)")
    ap.add_argument("--stereo", action="store_true",
                    help="decode the main as L/R-only stereo instead of "
                         "full 5.1 (E98: 5.1 runs 1.87x realtime here; use "
                         "this on a box where it cannot keep up)")
    ap.add_argument("--live-dir", default=None)
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
        procs.append(p)
        return p

    try:
        # 1. chain: video + audio + caption lanes; owns the radio
        chain = spawn([py, "tools/atsc3_run.py", "--rf", str(a.rf),
                       "--ant", a.ant, "--secs", "0", "--live-dir", live,
                       "--extra", "--assets all"])
        log(f"chain up (pid {chain.pid}) -- RF{a.rf} on {a.ant}")

        # 2. audio workers: BOTH languages, each its own process reading
        #    the lane file (E98). The 5.1 main -> live_audio.wav (6 ch,
        #    1.87x realtime on this box; --stereo for L/R only), the SAP
        #    pair -> live_audio_spa.wav. atsc3_tv muxes both; 'a' in the
        #    window switches. A missing SAP lane just yields silence on
        #    track 2 -- it never fails the chunk.
        pid_e, el_e = LANG["eng"]
        pid_s, el_s = LANG["spa"]
        audio = spawn([py, "tools/atsc3_audio.py", "--live-dir", live,
                       "--pid", str(pid_e), "--element", el_e,
                       "--channels", "2" if a.stereo else "6",
                       "--out", os.path.join(live, "live_audio.wav")])
        log(f"audio worker up (pid {audio.pid}) -- eng pid {pid_e} "
            f"{'stereo' if a.stereo else '5.1'}")
        audio2 = spawn([py, "tools/atsc3_audio.py", "--live-dir", live,
                        "--pid", str(pid_s), "--element", el_s,
                        "--channels", "2",
                        "--out", os.path.join(live, "live_audio_spa.wav")])
        log(f"audio worker up (pid {audio2.pid}) -- spa pid {pid_s} stereo")

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
        tv = [py, "tools/atsc3_tv.py", "--live-dir", live,
              "--mode", "v2", "--player", "ffplay",
              "--subs", "soft" if a.cc else "none",
              "--ffplay-audio", a.lang,
              "--exit-on-player-close"] + (["--stereo"] if a.stereo else [])
        log("starting viewer (atsc3_tv v2/ffplay: HEVC copy + soft CC) -- "
            "close the window to stop; 't' toggles captions")
        rc = spawn(tv, quiet=False).wait()
        log(f"viewer exited (rc={rc}) -- shutting down")
        return 0
    except KeyboardInterrupt:
        return 0
    finally:
        teardown(procs)


if __name__ == "__main__":
    sys.exit(main())
