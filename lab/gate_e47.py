#!/usr/bin/env python3
"""gate_e47.py -- the 132.1 ROUTE/DASH LIVE path, gated offline on capture.

WHAT IS UNDER TEST
  1. CHAIN SIDE (staged, lab/m11_stream.RouteStream + Transport(route=...)):
     the banked PLP1 baseband stream (m13_plp1_48.bb, 48 frames of real
     RF33 air, 117/117 BCH-clean per frame) is fed through the streaming
     Transport exactly as `atsc3 watch --route 239.255.32.1:8321` would,
     into the real LiveWriter.  REFEREE: the lane files must be
     BYTE-IDENTICAL to tools/atsc3_route.assemble_representation (the
     M14-proven batch path) on the same datagrams -- one transport, two
     schedules, zero drift.
  2. AUDIO WORKER: tools/atsc3_audio.py --kind route_audio --element pair
     on the freshly written lanes -> live_route_audio(.wav/_spa.wav).
  3. VIEWER (tools/atsc3_tv_route.py): --replay over the lanes, headless,
     no player (ONE WINDOW LAW: the user is watching 107.1) -> a TS that
     ffprobe must certify: 1080p HEVC present, two mp2 audio tracks,
     duration within tolerance, video DTS monotonic, audio AUDIBLE.

NEGATIVE CONTROLS (the gate must be able to fail):
  A. The MMTP lane parser on ROUTE input: mmtp_score fails the 0.9 flow
     gate at every ver_ext and MpuStreamer completes ZERO MPUs on the
     8321 datagrams.  The transports are provably not interchangeable.
  B. TODAY'S CHAIN CONFIG (Transport with no `route=`): the same air
     yields ZERO route segments and classifies the flow "other" -- the
     direct, code-level proof that the running chain does NOT emit ROUTE
     lanes and a restart with --route is required.

Never touches the radio, the GPU, the running chain, or anything under
data/.  Output goes to lab/e47_out/run_<stamp>/ (mkdir only -- the 8/08
law: run-dir scripts create, NEVER clear).

    python lab/gate_e47.py
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "tools"))

import m7_objects as O7                                           # noqa: E402
import m7_route as R7                                             # noqa: E402
import m11_stream as ST                                           # noqa: E402
import m11_watch as W                                             # noqa: E402
import atsc3_route as RT                                          # noqa: E402

PY = sys.executable
BB = os.path.join(HERE, "m13_plp1_48.bb")
DG = os.path.join(HERE, "m13_plp1_48.dg")
SLS = os.path.join(HERE, "plp1_out", "obj_239_255_32_1_8321_tsi0_toi4653059.bin")
FLOW = "239.255.32.1:8321"
WANT_ALL = (b"vide", b"soun", b"subt", b"text")

RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok)))
    print(f"    {'PASS' if ok else '*** FAIL ***'}  {name}"
          + (f"  ({detail})" if detail else ""))
    return ok


def ffbin(name):
    import shutil
    p = shutil.which(name)
    if p:
        return p
    c = os.path.join(os.environ.get("LOCALAPPDATA", "C:"), "Microsoft",
                     "WinGet", "Packages",
                     "Gyan.FFmpeg_Microsoft.Winget.Source_8wekyb3d8bbwe",
                     "ffmpeg-8.1-full_build", "bin", name + ".exe")
    return c if os.path.exists(c) else name


def feed_transport(tr, lw=None, chunks=64):
    """The banked bb stream through the STREAMING path, in chunks."""
    import bisect
    segs = []
    for _s0, _fr, stream, bounds in R7.bb_read(BB):
        bl = sorted(set(bounds))
        n = len(stream)
        step = max(1, n // chunks)
        j = 0
        for s in range(0, n, step):
            e = min(n, s + step)
            k = bisect.bisect_left(bl, e)
            out = tr.feed(bytes(stream[s:e]), bl[j:k])
            j = k
            for seg in out:
                segs.append(seg)
                if lw is not None:
                    lw.write(None, seg)
        for seg in tr.flush():
            segs.append(seg)
            if lw is not None:
                lw.write(None, seg)
    return segs


def main():
    print("E47 GATE -- 132.1 ROUTE/DASH live path, offline on banked capture")
    print("=" * 74)
    stamp = time.strftime("%m%d_%H%M%S")
    out_root = os.path.join(HERE, "e47_out")
    live = os.path.join(out_root, f"run_{stamp}")
    os.makedirs(live, exist_ok=True)          # mkdir, NEVER clear (8/08 law)
    print(f"  run dir: {live}")

    # ------------------------------------------------------------------
    print("\n  [1] REFERENCE (M14 batch path, tools/atsc3_route.py verbatim)")
    objs, ndg = RT.route_objects(DG, FLOW)
    sid, reps = RT.parse_sls(open(SLS, "rb").read())
    tsi_of = {r["rep_id"]: r["tsi"] for r in reps}
    ref = {}
    for rep_id in ("TS-4000_1_video", "TS-4000_1_audio_1", "TS-4000_1_audio_2"):
        b, man = RT.assemble_representation(objs, tsi_of[rep_id])
        ref[rep_id] = (b, man)
        print(f"      {rep_id}: {man['n_used']} complete segments, "
              f"TOI {man['used_tois'][0]}..{man['used_tois'][-1]}, "
              f"{len(b)} bytes")
    n_vid_ref = ref["TS-4000_1_video"][1]["n_used"]

    # ------------------------------------------------------------------
    print("\n  [2] CHAIN SIDE, STAGED: Transport(route=...) -> LiveWriter")
    tr = ST.Transport(probe=64, want=WANT_ALL, route=[FLOW])
    lw = W.LiveWriter(live)
    segs = feed_transport(tr, lw)
    W.dump_route_info(tr, live)
    lw.close()
    ip_s, port_s = FLOW.rsplit(":", 1)
    key = (bytes(int(x) for x in ip_s.split(".")), int(port_s))
    check("flow classified ROUTE", tr.kind.get(key, ("?",))[0] == "route",
          f"kind={tr.kind.get(key)}")
    lanes = (json.load(open(os.path.join(live, "live.json")))
             .get("lanes", {}))
    vlane = next((l for l in lanes.values()
                  if l.get("kind") == "route_video"), None)
    alane = {l["pid"]: l for l in lanes.values()
             if l.get("kind") == "route_audio"}
    check("route_video lane exists (pid 110)",
          vlane is not None and vlane["pid"] == 110)
    check("route_audio lanes exist (pid 120 eng, 130 spa)",
          120 in alane and 130 in alane)
    check(f"video lane carries the {n_vid_ref} complete segments",
          vlane and vlane["segments"] == n_vid_ref
          and vlane["first_seq"] == ref["TS-4000_1_video"][1]["used_tois"][0]
          and vlane["last_seq"] == ref["TS-4000_1_video"][1]["used_tois"][-1],
          f"segments={vlane and vlane['segments']} "
          f"seq={vlane and vlane['first_seq']}..{vlane and vlane['last_seq']}")
    # THE REFEREE: byte identity with the batch transport
    for rep_id, pid in (("TS-4000_1_video", 110), ("TS-4000_1_audio_1", 120),
                        ("TS-4000_1_audio_2", 130)):
        lane = next((l for l in lanes.values() if l.get("pid") == pid), None)
        got = open(lane["path"], "rb").read() if lane else b""
        check(f"lane pid {pid} BYTE-IDENTICAL to batch assemble ({rep_id})",
              got == ref[rep_id][0],
              f"{len(got)} vs {len(ref[rep_id][0])} bytes")
    idx = [json.loads(l) for l in open(os.path.splitext(vlane["path"])[0]
                                       + ".idx")]
    seqs = [r["seq"] for r in idx]
    check("video .idx TOIs strictly ascending",
          all(b > a for a, b in zip(seqs, seqs[1:])), f"{seqs}")
    ri = json.load(open(os.path.join(live, "live_route.json")))
    check("live_route.json: serviceId 1, video+2 audio+cc reps",
          ri.get("service_id") == 1 and len(ri.get("reps", [])) == 4,
          f"sid={ri.get('service_id')} reps={len(ri.get('reps', []))}")

    # ------------------------------------------------------------------
    print("\n  [3] NEGATIVE CONTROL A: the MMTP lane parser on ROUTE input")
    pays = [pay for (_s, dst, _sp, dp, pay) in O7.read_dg(DG)
            if f"{dst}:{dp}" == FLOW]
    sc = {v: O7.mmtp_score(pays, v) for v in (0, 2, 4)}
    n_lct = sum(1 for p in pays[:256] if O7.lct_parse(p) is not None)
    mmtp_would_take = any(tot > 0 and ok >= 0.9 * tot and ok > n_lct
                          for ok, tot in sc.values())
    check("flow gate REFUSES ROUTE datagrams as MMTP", not mmtp_would_take,
          f"scores={{{', '.join('%d: %d/%d' % (v, o, t) for v, (o, t) in sc.items())}}}, "
          f"LCT {n_lct}/256")
    ms = ST.MpuStreamer(0)
    for p in pays:
        ms.feed(p)
    # (a few dozen garbage headers PARSE -- the 12-byte MMTP header has no
    # magic to refuse -- but not one MPU ever completes: the reassembly
    # contract, not the header read, is what rejects the wrong transport)
    check("MpuStreamer on ROUTE input completes ZERO MPUs",
          ms.stats["complete"] == 0,
          f"complete={ms.stats['complete']} (header-parsed "
          f"{ms.stats['mmtp']}/{len(pays)}, all discarded)")

    # ------------------------------------------------------------------
    print("\n  [4] NEGATIVE CONTROL B: TODAY'S chain config (no route=)")
    tr2 = ST.Transport(probe=64, want=WANT_ALL)          # the running config
    segs2 = feed_transport(tr2)
    check("default Transport emits ZERO segments from the ROUTE air",
          len(segs2) == 0, f"{len(segs2)} segments")
    check("default Transport classifies the flow 'other' (dropped)",
          tr2.kind.get(key, ("?",))[0] == "other"
          and not getattr(tr2, "routes", {}),
          f"kind={tr2.kind.get(key)}")

    # ------------------------------------------------------------------
    print("\n  [5] AUDIO WORKER on the lanes (--kind route_audio, pair)")
    for pid, out_wav in ((120, "live_route_audio.wav"),
                         (130, "live_route_audio_spa.wav")):
        r = subprocess.run(
            [PY, os.path.join(ROOT, "tools", "atsc3_audio.py"),
             "--live-dir", live, "--kind", "route_audio", "--pid", str(pid),
             "--element", "pair", "--channels", "2", "--once",
             "--out", os.path.join(live, out_wav)],
            capture_output=True, text=True, cwd=ROOT, timeout=600)
        wav = os.path.join(live, out_wav)
        n = (os.path.getsize(wav) - 44) / 4 / 48000 if os.path.exists(wav) \
            else 0.0
        check(f"worker pid {pid} -> {out_wav}",
              r.returncode == 0 and n >= 0.9 * n_vid_ref * 2.002,
              f"rc={r.returncode}, {n:.2f} s of PCM"
              + ("" if r.returncode == 0 else
                 " | " + (r.stdout + r.stderr).strip()[-160:]))

    # ------------------------------------------------------------------
    print("\n  [6] VIEWER: atsc3_tv_route.py --replay, headless, no player")
    ts_path = os.path.join(live, "e47.ts")
    r = subprocess.run(
        [PY, os.path.join(ROOT, "tools", "atsc3_tv_route.py"),
         "--live-dir", live, "--replay", "--fast", "--player", "none",
         "--chunk", "1", "--audio-hold", "0", "--out", ts_path],
        capture_output=True, text=True, cwd=ROOT, timeout=600)
    check("viewer replay exits clean (0 failed chunks)", r.returncode == 0,
          f"rc={r.returncode} | " + r.stdout.strip().splitlines()[-1]
          if r.stdout.strip() else f"rc={r.returncode}")
    pr = subprocess.run(
        [ffbin("ffprobe"), "-v", "error", "-show_entries",
         "stream=codec_type,codec_name,width,height:format=duration",
         "-of", "default=noprint_wrappers=1", ts_path],
        capture_output=True, text=True)
    probe = pr.stdout
    dm = re.search(r"duration=([\d.]+)", probe)
    dur = float(dm.group(1)) if dm else 0.0
    expect = n_vid_ref * 2.002
    check("TS: 1080p HEVC video stream present",
          "codec_name=hevc" in probe and "width=1920" in probe
          and "height=1080" in probe)
    check("TS: two mp2 audio tracks (eng+spa)",
          probe.count("codec_name=mp2") >= 2)
    check(f"TS duration ~ {expect:.3f} s", abs(dur - expect) < 0.6,
          f"{dur:.3f} s")
    pk = subprocess.run(
        [ffbin("ffprobe"), "-v", "error", "-select_streams", "v:0",
         "-show_entries", "packet=dts", "-of", "csv=p=0", ts_path],
        capture_output=True, text=True)
    dts = [int(x) for x in (t.strip().strip(",") for t in pk.stdout.split())
           if x.lstrip("-").isdigit()]
    mono = all(b >= a for a, b in zip(dts, dts[1:]))
    check("video DTS monotonic through the whole TS",
          len(dts) >= 100 and mono, f"{len(dts)} packets")
    ns = subprocess.run(
        [ffbin("ffmpeg"), "-hide_banner", "-i", ts_path, "-f", "null", "-"],
        capture_output=True, text=True)
    fm = None
    for m in re.finditer(r"frame=\s*(\d+)", ns.stderr):
        fm = int(m.group(1))
    exp_frames = round(expect * 60000 / 1001)
    check(f"null-sink decodes ~{exp_frames} HEVC frames",
          fm is not None and abs(fm - exp_frames) <= 2, f"{fm} frames")
    vd = subprocess.run(
        [ffbin("ffmpeg"), "-hide_banner", "-i", ts_path, "-map", "0:a:0",
         "-af", "volumedetect", "-f", "null", "-"],
        capture_output=True, text=True)
    mv = re.search(r"mean_volume:\s*(-?[\d.]+)\s*dB", vd.stderr)
    mv = float(mv.group(1)) if mv else None
    check("eng audio AUDIBLE (mean_volume > -50 dB)",
          mv is not None and mv > -50.0, f"mean {mv} dB")

    # ------------------------------------------------------------------
    n_ok = sum(1 for _n, ok in RESULTS if ok)
    print("\n" + "=" * 74)
    print(f"  E47 GATE: {n_ok}/{len(RESULTS)} checks pass"
          + ("" if n_ok == len(RESULTS) else "  -- FAILING:"))
    for n, ok in RESULTS:
        if not ok:
            print(f"      FAIL  {n}")
    print("=" * 74)
    return 0 if n_ok == len(RESULTS) else 1


if __name__ == "__main__":
    sys.exit(main())
