#!/usr/bin/env python3
"""gate_tv2.py -- headless replay gates for the v2 viewer (E35).

For each replay dir (live51 = 26.6 min with real churn holes, catchup =
clean), this:

  1. runs tools/atsc3_tv.py --replay --fast --player none, muxing the whole
     dir into a gate TS (the SAME code path the window uses, minus pacing);
  2. requires zero failed chunks and a clean exit;
  3. ffprobes the TS: EXACTLY hevc + mp2(eng) + mp2(spa) + dvb_subtitle,
     container duration within 1 % of the emitted slot span;
  4. requires the dvbsub stream to actually carry cue packets (>= 2 per
     minute of programme on these dirs -- both have busy captioning).

NEGATIVE CONTROLS (a gate that cannot fail is not a gate): the same checks
are run against (a) a half-truncated copy -- the duration check must FAIL;
(b) a remux with the Spanish track dropped -- the layout check must FAIL.
The gate passes only if every positive passes AND every control fails.

E37 adds the lead-governor cases (the live-window glitch fix): steady
pass-through, starvation withhold+burst, no re-pin after recovery,
max-hold drain, and control D -- the ungoverned pass-through policy on
the same starved schedule must land EVERY post-gap chunk on a starved
player (the measured EOF metronome).

Usage:
    python lab/gate_tv2.py                     # both dirs + controls
    python lab/gate_tv2.py --dirs data/catchup
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
MPU_SECONDS = 2.002


def ffbin(name):
    p = shutil.which(name)
    if p:
        return p
    c = os.path.join(
        os.environ.get("LOCALAPPDATA", "C:"), "Microsoft", "WinGet",
        "Packages", "Gyan.FFmpeg_Microsoft.Winget.Source_8wekyb3d8bbwe",
        "ffmpeg-8.1-full_build", "bin", name + ".exe")
    return c if os.path.exists(c) else name


def probe(ts):
    r = subprocess.run([ffbin("ffprobe"), "-v", "error", "-show_streams",
                        "-show_format", "-of", "json", ts],
                       capture_output=True, text=True, errors="replace")
    try:
        return json.loads(r.stdout)
    except ValueError:
        return {}


def sub_packet_count(ts):
    r = subprocess.run([ffbin("ffprobe"), "-v", "error", "-select_streams",
                        "s", "-show_packets", "-show_entries", "packet=size",
                        "-of", "csv", ts],
                       capture_output=True, text=True, errors="replace")
    return sum(1 for l in (r.stdout or "").splitlines() if l.strip())


def check_layout(ts):
    """-> (ok, detail): exactly hevc + mp2 eng + mp2 spa + dvb_subtitle."""
    st = probe(ts).get("streams", [])
    kinds = [(s.get("codec_type"), s.get("codec_name"),
              (s.get("tags") or {}).get("language")) for s in st]
    want = [("video", "hevc", None),
            ("audio", "mp2", "eng"),
            ("audio", "mp2", "spa"),
            ("subtitle", "dvb_subtitle", None)]
    if len(kinds) != 4:
        return False, f"stream count {len(kinds)} != 4: {kinds}"
    for got, w in zip(kinds, want):
        if got[0] != w[0] or got[1] != w[1]:
            return False, f"stream mismatch: {got} != {w}"
        if w[2] and got[2] != w[2]:
            return False, f"language mismatch: {got} != {w}"
    return True, "hevc + mp2(eng) + mp2(spa) + dvb_subtitle"


def check_duration(ts, expected_s, tol=0.01):
    fmt = probe(ts).get("format", {})
    try:
        dur = float(fmt.get("duration", 0.0))
    except (TypeError, ValueError):
        dur = 0.0
    ok = expected_s > 0 and abs(dur - expected_s) / expected_s <= tol
    return ok, f"duration {dur:.2f}s vs expected {expected_s:.2f}s"


def run_dir(d, chunk):
    print(f"\n=== gate: {d} ===")
    out_ts = os.path.join(d, "_tv", "gate_v2.ts")
    os.makedirs(os.path.dirname(out_ts), exist_ok=True)
    r = subprocess.run(
        [sys.executable, os.path.join(ROOT, "tools", "atsc3_tv.py"),
         "--live-dir", d, "--replay", "--fast", "--player", "none",
         "--mode", "v2", "--chunk", str(chunk), "--out", out_ts],
        capture_output=True, text=True, errors="replace", cwd=ROOT)
    outp = (r.stdout or "") + (r.stderr or "")
    fails = len(re.findall(r"mux failed", outp))
    m = re.search(r"replay complete: ([0-9.]+) s emitted, (\d+) failed",
                  outp)
    emitted = float(m.group(1)) if m else 0.0
    results = []
    results.append(("tool exit 0", r.returncode == 0,
                    f"exit {r.returncode}"))
    results.append(("all chunks muxed", fails == 0 and m is not None,
                    f"{fails} failures, emitted {emitted:.1f}s"))
    # expected span: whole chunks only (the tail partial chunk is by
    # design never emitted)
    lanes = json.load(open(os.path.join(d, "live.json")))["lanes"]
    v = next(l for l in lanes.values() if l["kind"] == "video")
    nslots = v["last_seq"] - v["first_seq"] + 1
    expect = (nslots // chunk) * chunk * MPU_SECONDS
    ok, det = check_layout(out_ts)
    results.append(("stream layout", ok, det))
    ok, det = check_duration(out_ts, expect)
    results.append(("container duration", ok, det))
    nsub = sub_packet_count(out_ts)
    need = max(2, int(2 * expect / 60.0))
    results.append((f"dvbsub packets >= {need}", nsub >= need,
                    f"{nsub} packets"))
    # VIDEO IS ACTUALLY IN THERE. A too-eager corrupt-fragment screen
    # dropped 388/388 fragments during development and every other check
    # above still passed (audio carries the duration). Gate the picture
    # itself: the copied HEVC must deliver nearly every frame the lane's
    # present fragments hold (120 per slot at 59.94).
    r = subprocess.run([ffbin("ffprobe"), "-v", "error", "-select_streams",
                        "v", "-count_packets", "-show_entries",
                        "stream=nb_read_packets", "-of", "csv", out_ts],
                       capture_output=True, text=True, errors="replace")
    m2 = re.search(r"stream,(\d+)", r.stdout or "")
    nvid = int(m2.group(1)) if m2 else 0
    vidx = os.path.splitext(os.path.join(ROOT, v["path"]))[0] + ".idx"
    nfrag = sum(1 for _ in open(vidx))
    vneed = int(0.85 * nfrag * 120)
    results.append((f"video packets >= {vneed}", nvid >= vneed,
                    f"{nvid} packets from {nfrag} fragments"))
    # THE TS ITSELF must carry audible sound on each track that has a
    # decoded wav (mux plumbing end-to-end, not just wav-side RMS). Listen
    # where the WAV says there is sound -- live51's eng wav opens with
    # ~3 min of genuinely silent air (recorded through a fade), and the
    # spa wav may be anchored mid-lane (--start-behind).
    import numpy as np
    sys.path.insert(0, os.path.join(ROOT, "tools"))
    import atsc3_tv as T
    sess = T.Session(d, 0.0)
    for label, stream, track in (("eng", "0:a:0", sess.eng),
                                 ("spa", "0:a:1", sess.spa)):
        if track.path is None or track.first is None:
            continue
        at = None
        for s in range(track.first, track.first + 750, 15):
            pcm = T.audio_slice(sess, d, s, s + 3, track)
            if float(np.sqrt((pcm.astype("f4") ** 2).mean())) > 300.0:
                at = (s - v["first_seq"]) * MPU_SECONDS
                break
        if at is None or at > expect - 10:
            results.append((f"TS {label} track audible", False,
                            "no loud slot found in wav to listen at"))
            continue
        r = subprocess.run([ffbin("ffmpeg"), "-loglevel", "info",
                            "-ss", f"{at:.3f}", "-i", out_ts, "-map", stream,
                            "-t", "30", "-af", "volumedetect",
                            "-f", "null", "-"],
                           capture_output=True, text=True, errors="replace")
        m3 = re.search(r"mean_volume:\s*(-?[0-9.]+) dB", r.stderr or "")
        vol = float(m3.group(1)) if m3 else -91.0
        results.append((f"TS {label} track audible near t={at:.0f}s",
                        vol > -60.0, f"mean_volume {vol} dB"))
    return results, out_ts, expect


def anchored_audio_case():
    """E35 addendum gate, E37-hardened: the mid-lane anchor case is now
    SYNTHESIZED. (It originally leaned on live51's Spanish wav being
    anchored 621 slots in; on 8/08 ~01:36 a worker outside this session
    re-decoded that lane from fragment 0, so the banked dir stopped
    exhibiting its own precondition -- gate-corpus drift. The mapping
    under test is unchanged.) A sidecar anchored 621 audio-idx entries
    into the lane + a 2.0M-frame tone wav: the anchored mapping must
    find the tone 10 entries past the anchor and silence before it; the
    old absolute k*SPF mapping must FAIL (reads past EOF)."""
    import types
    import wave as _w
    import numpy as np
    sys.path.insert(0, os.path.join(ROOT, "tools"))
    import atsc3_tv as T
    d = "data/live51"
    results = []
    sess = T.Session(d, 0.0)
    a_idx = T.read_idx(os.path.join(d, "live_audio_pid14.idx"))
    anchor = a_idx[621][0]
    nf = 2_000_000
    tmpw = os.path.join(d, "_tv", "gate_midlane.wav")
    ww = _w.open(tmpw, "wb")
    ww.setnchannels(2)
    ww.setsampwidth(2)
    ww.setframerate(48000)
    tone = (3000 * np.sin(2 * np.pi * 440 / 48000
                          * np.arange(nf, dtype="f4"))).astype("<i2")
    ww.writeframes(np.stack([tone, tone], 1).tobytes())
    ww.close()
    trk = types.SimpleNamespace(first=anchor, path=tmpw, pid=14, ch=2)
    lane0 = sess.v_idx()[0][0]
    results.append(("synthetic anchor is mid-lane", anchor > lane0,
                    f"first_seq {anchor} vs lane start {lane0}"))
    s_probe = a_idx[631][0]          # 10 entries past the anchor
    pcm = T.audio_slice(sess, d, s_probe, s_probe + 1, trk)
    rms = float(np.sqrt((pcm.astype("f4") ** 2).mean()))
    results.append(("anchored mapping finds audio", rms > 500.0,
                    f"rms {rms:.1f} @ anchor+10 entries"))
    s_pre = a_idx[600][0]            # before the anchor
    pre = T.audio_slice(sess, d, s_pre, s_pre + 1, trk)
    rms_pre = float(np.sqrt((pre.astype("f4") ** 2).mean()))
    results.append(("slots before anchor are silence", rms_pre == 0.0,
                    f"rms {rms_pre:.1f} @ anchor-21 entries"))
    # negative control: absolute mapping on the same slot
    k = 631
    absolute_fails = (k + 1) * T.SPF > nf
    results.append(("control C: absolute k*SPF mapping FAILS mid-lane wav",
                    absolute_fails,
                    f"wants frames up to {(k + 1) * T.SPF}, wav has {nf}"))
    try:
        os.remove(tmpw)
    except OSError:
        pass
    return results


def lead_governor_case():
    """E37 gate: the LeadGovernor (the glitch fix) must
      (1) pass a steady 1x feed through untouched (no withholding),
      (2) after a starvation gap, withhold-and-burst so the modelled
          cushion is rebuilt to `lead` and never re-pins afterwards,
      (3) age-flush a dying feed (no chunk may be held forever).
    NEGATIVE CONTROL: the ungoverned pass-through policy (the shipped
    v2 behaviour the user watched) on the SAME starved schedule arrives
    at every post-gap chunk with zero modelled lead -- the measured
    ~6 s EOF metronome. That control must FAIL the smooth criterion."""
    sys.path.insert(0, os.path.join(ROOT, "tools"))
    import atsc3_tv as T
    C = 3 * MPU_SECONDS              # one chunk of media, 6.006 s
    results = []

    def feed(g, t, n, step=None):
        """n chunks at 1x from wall t; returns (t, outs)."""
        outs = []
        for _ in range(n):
            outs.append(g.push(b"x", C, t))
            if g.play0 is None and g.want_player():
                g.player_started(t)
            t += step if step is not None else C
        return t, outs

    # (1) steady feed: spawn at the 24s bank, 100 chunks, zero rebuilds
    g = T.LeadGovernor(lead=24.0, safety=3.0)
    t, outs = feed(g, 0.0, 104)
    held = sum(1 for o in outs if o == b"")
    lead_end = g.lead_now(t)
    results.append(("governor: steady feed passes through",
                    held == 0 and g.rebuilds == 0 and lead_end > 18.0,
                    f"{held} held, {g.rebuilds} rebuilds, "
                    f"final lead {lead_end:.1f}s"))

    # (2) 40s starvation, then a recovered 1x feed
    g = T.LeadGovernor(lead=24.0, safety=3.0)
    t, _ = feed(g, 0.0, 10)          # healthy: banked + playing
    t += 40.0                        # the fade: nothing arrives
    t, outs = feed(g, t, 4)          # recovery begins
    burst = outs[-1]
    ok_hold = outs[:3] == [b"", b"", b""] and burst == b"xxxx"
    results.append(("governor: gap triggers withhold, 4-chunk burst",
                    ok_hold and g.rebuilds == 1,
                    f"outs {[len(o) for o in outs]}, "
                    f"{g.rebuilds} rebuilds, lead {g.lead_now(t):.1f}s"))
    min_lead = 1e9
    t2, outs2 = t, []
    for _ in range(50):              # feed stays healthy ever after
        outs2.append(g.push(b"x", C, t2))
        min_lead = min(min_lead, g.lead_now(t2))
        t2 += C
    results.append(("governor: cushion never re-pins after recovery",
                    g.rebuilds == 1 and min_lead > 3.0
                    and all(o == b"x" for o in outs2),
                    f"min lead {min_lead:.1f}s, {g.rebuilds} total rebuilds"))

    # (3) age-flush: one chunk withheld, then the feed dies
    g = T.LeadGovernor(lead=24.0, safety=3.0, max_hold=90.0)
    t, _ = feed(g, 0.0, 10)
    t += 40.0
    t, outs = feed(g, t, 1)          # withheld (rebuilding)
    drained = g.tick(t + 91.0)
    results.append(("governor: max_hold drains a dying feed",
                    outs == [b""] and drained == b"x",
                    f"tick returned {len(drained)} bytes after 91s"))

    # NEGATIVE CONTROL: ungoverned pass-through on the same schedule.
    # Model: playhead p advances at 1x, never past appended media;
    # a chunk that arrives when media - p <= 0 lands on a starved player.
    media, p, t_last = 0.0, 0.0, None
    playing = False
    starved_arrivals = 0
    sched = [i * C for i in range(10)]                       # healthy
    resume0 = sched[-1] + C + 40.0
    sched += [resume0 + i * C for i in range(50)]            # post-gap
    for i, tw in enumerate(sched):
        if t_last is not None and playing:
            p = min(media, p + (tw - t_last))
        t_last = tw
        if tw >= resume0 and media - p <= 0.5:   # < half a second of lead
            starved_arrivals += 1
        media += C                   # pass-through appends immediately
        if not playing and media >= 2 * C:
            playing = True           # v2 spawned VLC after 2 chunks
            p = 0.0
    results.append(("control D: ungoverned pass-through metronomes "
                    "(every post-gap chunk lands on a starved player)",
                    starved_arrivals == 50,
                    f"{starved_arrivals}/50 post-gap chunks arrived at "
                    f"zero lead"))
    return results


def negative_controls(good_ts, expect):
    print("\n=== negative controls (must FAIL) ===")
    nd = os.path.dirname(good_ts)
    results = []
    # (a) truncated -> duration must fail
    trunc = os.path.join(nd, "gate_v2_trunc.ts")
    n = os.path.getsize(good_ts) // 2
    with open(good_ts, "rb") as s, open(trunc, "wb") as t:
        t.write(s.read(n))
    ok, det = check_duration(trunc, expect)
    results.append(("control A: truncated duration check FAILS",
                    not ok, det))
    # (b) spanish dropped -> layout must fail
    nospa = os.path.join(nd, "gate_v2_nospa.ts")
    subprocess.run([ffbin("ffmpeg"), "-loglevel", "error", "-y",
                    "-i", good_ts, "-map", "0:0", "-map", "0:1",
                    "-map", "0:3", "-c", "copy", "-f", "mpegts", nospa],
                   capture_output=True)
    ok, det = check_layout(nospa)
    results.append(("control B: no-spa layout check FAILS", not ok, det))
    for f in (trunc, nospa):
        try:
            os.remove(f)
        except OSError:
            pass
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dirs", nargs="*",
                    default=["data/live51", "data/catchup"])
    ap.add_argument("--chunk", type=int, default=3)
    a = ap.parse_args()
    os.chdir(ROOT)
    all_ok = True
    control_ts = None
    control_expect = None
    for d in a.dirs:
        results, out_ts, expect = run_dir(d, a.chunk)
        control_ts, control_expect = out_ts, expect
        for name, ok, det in results:
            print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {det}")
            all_ok &= ok
    if control_ts and os.path.exists(control_ts):
        for name, ok, det in negative_controls(control_ts, control_expect):
            print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {det}")
            all_ok &= ok
    print("\n=== anchored-audio case (E35 addendum) ===")
    for name, ok, det in anchored_audio_case():
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {det}")
        all_ok &= ok
    print("\n=== lead-governor case (E37 glitch fix) ===")
    for name, ok, det in lead_governor_case():
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {det}")
        all_ok &= ok
    print(f"\nGATE {'PASS' if all_ok else 'FAIL'}")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
