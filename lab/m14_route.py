#!/usr/bin/env python3
"""M14-ROUTE -- assemble channel 132.1 (ROUTE/DASH) into a played clip.

This is the DRIVER for the milestone.  The reusable transport lives in
`tools/atsc3_route.py` (SLS parse + segment assemble); the AC-4 decoder is
the existing `m42_ac4_stream.Ac4Stream` + `m30_filterbank` + `m38_resample`.
This file only wires them together for 132.1 and GATES the result.

  video  Representation TS-4000_1_video  (tsi 10)  HEVC 1080p  -> fMP4, copy
  audio  Representation TS-4000_1_audio_1 (tsi 20)  AC-4 stereo -> PCM -> AAC

Pipeline:
  1. atsc3_route.route_objects + assemble_representation -> video.mp4, audio.mp4
  2. AC-4 stereo: Ac4Stream(element="pair").decode_frames -> synthesise_frames
     -> m38.resample (1001/960 to the true 48 kHz rate) -> stereo wav
  3. ffmpeg mux: HEVC copy + AAC(from PCM) -> plp1_out/ch132_1.mp4
  4. GATES (each with a primitive that states itself):
       (a) ffprobe: 1920x1080 hevc + an audio stream, duration ~ N*2.002 s
       (b) ffmpeg null-sink: video decodes with ~0 errors at 59.94 fps
       (c) volumedetect: mean_volume > -50 dB  (AUDIBLE, not silence)
  5. NEGATIVE CONTROL: assemble the media segments in SHUFFLED order and show
     the decode's presentation timeline breaks (non-monotonic DTS), proving
     TOI order is load-bearing and the in-order assembly is the right one.
"""
from __future__ import annotations

import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "tools"))

import atsc3_route as RT                                          # noqa: E402
import m17_ac4_walk as W                                          # noqa: E402
import m30_filterbank as FB                                       # noqa: E402
import m37_render_hf as R37                                       # noqa: E402
import m38_resample as RS                                         # noqa: E402
import m29_audio as R29                                           # noqa: E402
from m42_ac4_stream import Ac4Stream                              # noqa: E402

DG = os.path.join(HERE, "m13_plp1_48.dg")
SLS = os.path.join(HERE, "plp1_out",
                   "obj_239_255_32_1_8321_tsi0_toi4653059.bin")
FLOW = "239.255.32.1:8321"
OUTDIR = os.path.join(HERE, "plp1_out")
SEG_SECONDS = 180180 / 90000.0          # 2.002 s per DASH segment (MPD)


def run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


def ffprobe_streams(path):
    r = run(["ffprobe", "-v", "error", "-show_entries",
             "stream=index,codec_type,codec_name,width,height,r_frame_rate:"
             "format=duration", "-of", "default=noprint_wrappers=1", path])
    return r.stdout


def null_sink(path):
    """-> (returncode, n_frames, fps, stderr_tail)."""
    r = run(["ffmpeg", "-hide_banner", "-i", path, "-f", "null", "-"])
    import re
    frames, fps = None, None
    for m in re.finditer(r"frame=\s*(\d+)\s+fps=\s*([\d.]+)", r.stderr):
        frames = int(m.group(1))
    m = re.search(r"([\d.]+)\s*fps,", r.stderr)
    if m:
        fps = float(m.group(1))
    # decode-error lines that ffmpeg emits on broken bitstreams / bad order
    err = [ln for ln in r.stderr.splitlines()
           if any(s in ln for s in ("Invalid", "error", "corrupt",
                                    "non monotonically", "Non-monoton",
                                    "could not", "does not"))]
    return r.returncode, frames, fps, err


def mean_volume(path):
    r = run(["ffmpeg", "-hide_banner", "-i", path, "-map", "0:a",
             "-af", "volumedetect", "-f", "null", "-"])
    import re
    m = re.search(r"mean_volume:\s*(-?[\d.]+)\s*dB", r.stderr)
    mx = re.search(r"max_volume:\s*(-?[\d.]+)\s*dB", r.stderr)
    return (float(m.group(1)) if m else None,
            float(mx.group(1)) if mx else None)


# ---------------------------------------------------------------------------

def assemble(shuffle=False, seed=1):
    """-> (video_path, audio_path, video_manifest, audio_manifest)."""
    objs, ndg = RT.route_objects(DG, FLOW)
    _sid, reps = RT.parse_sls(SLS)
    vtsi = next(r["tsi"] for r in reps if r["rep_id"] == "TS-4000_1_video")
    atsi = next(r["tsi"] for r in reps if r["rep_id"] == "TS-4000_1_audio_1")

    vbytes, vman = RT.assemble_representation(objs, vtsi)
    abytes, aman = RT.assemble_representation(objs, atsi)

    if shuffle:
        # Rebuild the VIDEO stream with its media segments in a deliberately
        # wrong order (init still first). This is the negative control.
        init = objs[(vtsi, RT.INIT_TOI)]
        tois = list(vman["used_tois"])
        rng = np.random.default_rng(seed)
        order = tois.copy()
        while order == tois:                 # guarantee a real permutation
            rng.shuffle(order)
        buf = bytearray(init["data"])
        for toi in order:
            buf += objs[(vtsi, toi)]["data"]
        vpath = os.path.join(OUTDIR, "ch132_1_video_shuffled.mp4")
        open(vpath, "wb").write(bytes(buf))
        vman = dict(vman, used_tois=order, shuffled_from=tois)
        return vpath, None, vman, aman

    vpath = os.path.join(OUTDIR, "ch132_1_video.mp4")
    apath = os.path.join(OUTDIR, "ch132_1_audio.mp4")
    open(vpath, "wb").write(vbytes)
    open(apath, "wb").write(abytes)
    return vpath, apath, vman, aman


def decode_audio(audio_mp4, wav_out):
    """AC-4 stereo fMP4 -> resampled 48 kHz stereo wav. -> (n_frames, seconds)."""
    frames = W.samples(audio_mp4)
    dec = Ac4Stream(element="pair")           # 132.1 audio_1 is a stereo pair
    wins, lens, groups, counts = dec.decode_frames(frames, chans=("L", "R"))
    pcm, _ = FB.synthesise_frames(wins, lens, counts, 1536,
                                  states={"L": None, "R": None},
                                  return_state=True)
    L = np.asarray(pcm["L"], np.float64)
    Rr = np.asarray(pcm["R"], np.float64)
    if dec.cfg is not None:                   # A-SPX high band, if signalled
        out, _n = R37.apply_hf_pair(L, Rr, groups["lr"], dec.cfg)
        L, Rr = np.asarray(out["L"]), np.asarray(out["R"])
    # M38: 1536 internal samples/frame -> 1601.6 output samples, the true
    # presentation rate. Without this the sound is 4.27 % fast against picture.
    stereo = np.column_stack([L, Rr]) / 32768.0
    y = RS.resample(stereo)                   # 1001/960
    peak = np.abs(y).max() or 1.0
    y = np.clip(y / peak * 0.95, -1, 1)
    R29.write_wav(wav_out, [y[:, 0], y[:, 1]])
    return dec.n_ok, dec.n_bad, y.shape[0] / RS.OUT_FS


def mux(video_mp4, wav, out_mp4):
    r = run(["ffmpeg", "-hide_banner", "-y", "-i", video_mp4, "-i", wav,
             "-map", "0:v:0", "-map", "1:a:0", "-c:v", "copy",
             "-c:a", "aac", "-b:a", "192k", "-shortest", out_mp4])
    return r.returncode, r.stderr


def main():
    print("M14-ROUTE -- channel 132.1 ROUTE/DASH assembly")
    print("=" * 74)

    # 1-2. assemble + decode ------------------------------------------------
    vpath, apath, vman, aman = assemble()
    n_v = len(vman["used_tois"])
    n_a = len(aman["used_tois"])
    print(f"  video: {n_v} segments TOI {vman['used_tois']}  init "
          f"{vman['init_state']}")
    print(f"  audio: {n_a} segments TOI {aman['used_tois']}  init "
          f"{aman['init_state']}")
    expect_s = n_v * SEG_SECONDS
    print(f"  expected duration {n_v} x {SEG_SECONDS:.3f} = {expect_s:.3f} s")

    wav = os.path.join(OUTDIR, "ch132_1_audio.wav")
    n_ok, n_bad, a_sec = decode_audio(apath, wav)
    print(f"  AC-4 decode: {n_ok} frames ok, {n_bad} bad -> {a_sec:.3f} s wav")

    out_mp4 = os.path.join(OUTDIR, "ch132_1.mp4")
    rc, mstderr = mux(vpath, wav, out_mp4)
    if rc != 0:
        print("  MUX FAILED\n" + mstderr[-800:])
        return 1
    print(f"  muxed -> {out_mp4}")

    # 4. gates --------------------------------------------------------------
    print("\n  GATES")
    print("  " + "-" * 72)
    probe = ffprobe_streams(out_mp4)
    print("  (a) ffprobe:")
    for ln in probe.strip().splitlines():
        print(f"        {ln}")
    import re
    has_v = "codec_name=hevc" in probe and "width=1920" in probe \
        and "height=1080" in probe
    has_a = "codec_type=audio" in probe
    dm = re.search(r"duration=([\d.]+)", probe)
    dur = float(dm.group(1)) if dm else 0.0
    ga = has_v and has_a and abs(dur - expect_s) / expect_s < 0.02
    print(f"      -> 1080p hevc {has_v}, audio {has_a}, duration {dur:.3f}s "
          f"vs {expect_s:.3f}s ({100*(dur-expect_s)/expect_s:+.2f}%)  "
          f"{'PASS' if ga else 'FAIL'}")

    rc_ns, frames, fps, errs = null_sink(out_mp4)
    exp_frames = round(expect_s * 60000 / 1001)
    gb = rc_ns == 0 and frames is not None and abs(frames - exp_frames) <= 2 \
        and fps is not None and abs(fps - 59.94) < 0.1 and len(errs) == 0
    print(f"  (b) null-sink: rc {rc_ns}, {frames} frames (expect ~{exp_frames})"
          f", {fps} fps, {len(errs)} error lines  {'PASS' if gb else 'FAIL'}")
    for e in errs[:4]:
        print(f"        ! {e.strip()[:90]}")

    mv, mx = mean_volume(out_mp4)
    gc = mv is not None and mv > -50.0
    print(f"  (c) volumedetect: mean_volume {mv} dB, max {mx} dB  "
          f"(> -50 dB = audible)  {'PASS' if gc else 'FAIL'}")

    # 5. negative control ---------------------------------------------------
    print("\n  NEGATIVE CONTROL (segments in WRONG order must break)")
    print("  " + "-" * 72)
    svpath, _a, svman, _am = assemble(shuffle=True)
    print(f"      correct order  {svman['shuffled_from']}")
    print(f"      shuffled order {svman['used_tois']}")
    rc_s, frames_s, fps_s, errs_s = null_sink(svpath)
    broke = len(errs_s) > 0 or rc_s != 0
    print(f"      null-sink on shuffled: rc {rc_s}, {frames_s} frames, "
          f"{len(errs_s)} error lines")
    for e in errs_s[:5]:
        print(f"        ! {e.strip()[:90]}")
    print(f"      -> shuffled stream {'BREAKS' if broke else 'did NOT break'}"
          f"  {'PASS (ordering is load-bearing)' if broke else 'FAIL'}")

    allpass = ga and gb and gc and broke
    print("\n" + "=" * 74)
    print("  M14-ROUTE " + ("COMPLETE -- 132.1 plays picture + sound"
                            if allpass else "INCOMPLETE -- see failing gate"))
    print("=" * 74)
    return 0 if allpass else 1


if __name__ == "__main__":
    raise SystemExit(main())
