#!/usr/bin/env python3
"""atsc3_sync.py -- one slot-aligned file from the live lanes.

`atsc3_play.py` follows the video lane. This is the other job: take video,
decoded audio and captions and put them on ONE timeline, aligned by MPU
slot rather than by file position.

WHY BYTE 0 IS THE WRONG ANCHOR
    Every asset arrives in 2.002 s MPUs carrying a shared sequence number,
    but the lanes do not start on the same one -- the caption lane loses
    nothing and usually begins earlier than video, and either lane can miss
    an MPU and be one fragment short from there on. So neither "start each
    file at its beginning" nor "-shortest" aligns anything: the first
    trims the ends while leaving the START offset, and both SHIFT
    everything after a hole instead of leaving a gap.

    The slot is the anchor. m39 does exactly this for the batch path.

Usage:
    python tools/atsc3_sync.py --out synced.mp4
    python tools/atsc3_sync.py --lag 60 --play
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

from atsc3_play import align, MPU_SECONDS, read_idx, read_lanes  # noqa: E402


def _channels(path):
    try:
        import wave
        return wave.open(path).getnchannels()
    except Exception:                                          # noqa: BLE001
        return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--live-dir", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--lag", type=float, default=0.0,
                    help="seconds back from the write head (0 = whole file)")
    ap.add_argument("--secs", type=float, default=0.0, help="limit duration")
    ap.add_argument("--play", action="store_true")
    a = ap.parse_args()
    d = a.live_dir or os.path.join(ROOT, "data", "atsc3_live")
    out = a.out or os.path.join(d, "synced.mp4")

    lag_slots = int(a.lag / MPU_SECONDS) if a.lag else 10 ** 9
    p = align(d, lag_slots)
    if not p:
        print(f"no aligned lanes in {d} -- run the chain with --assets all "
              f"and the audio/subs workers", file=sys.stderr)
        return 1

    print(f"aligning on MPU slot {p['start']}")
    print(f"  video  {os.path.basename(p['video']['path'])}  "
          f"fragment {p['v_at']}/{len(p['v_idx'])}  byte {p['v_off']}")
    if "audio" in p:
        print(f"  audio  {os.path.basename(p['audio'])}  "
              f"skip {p['a_skip']:.3f} s")
    else:
        print("  audio  MISSING -- run tools/atsc3_audio.py")
    if "subs" in p:
        print(f"  subs   live.srt  shift {p['s_shift']:+.3f} s")
    else:
        print("  subs   MISSING -- run tools/atsc3_subs.py")

    # ---- audio placed PER SLOT ---------------------------------------
    # Aligning only the first slot is not enough. Over two hours the lanes
    # diverged by 37 fragments = 74 s (video 3504, audio 3541), so every
    # unmatched fragment shifts the two apart permanently. A 21-minute
    # sample showed all lanes losing the SAME slot and made that look
    # impossible -- it was a coincidence of the sample.
    #
    # With correct padding each audio fragment contributes exactly 60 AC-4
    # frames = 96096 samples at 48 kHz, so fragment k occupies a known,
    # fixed window of the wav. Walk the VIDEO slots and take the audio for
    # that same slot, silence where the audio lane has none.
    if "audio" in p:
        import wave as _wave
        import numpy as _np
        a_idx = read_idx(os.path.splitext(
            [l for l in read_lanes(d).values()
             if l.get("kind") == "audio"][0]["path"])[0] + ".idx")
        a_slot = {s: k for k, (s, _, _) in enumerate(a_idx)}
        SPF = 96096                          # samples per fragment
        w = _wave.open(p["audio"])
        pcm = _np.frombuffer(w.readframes(w.getnframes()), "<i2")
        pcm = pcm.reshape(-1, w.getnchannels())
        # Walk the full slot RANGE, not the present fragments. The video
        # timeline holds a lost MPU as a 2.002 s freeze (gap held true, E23a),
        # so its PTS spans every slot from first to last whether or not the
        # fragment arrived. Laying audio at present-fragment positions is
        # contiguous -- each video hole slides all later sound one slot
        # early, and E28's 69 video holes summed to sound 136 s AHEAD of
        # picture by the end of the file. This was invisible before the gap
        # fix because the video used to slide identically; two consistent
        # wrongs agreed with each other and disagreed with the air.
        #
        # Walking the range also plays audio through a video freeze whenever
        # the AUDIO lane has that slot -- which is what real television does:
        # the picture stalls, the sound carries on.
        first_s = p["v_idx"][p["v_at"]][0]
        last_s = p["v_idx"][-1][0]
        n_slots = last_s - first_s + 1
        out_a = _np.zeros((n_slots * SPF, pcm.shape[1]), "<i2")
        hit = 0
        for i, s in enumerate(range(first_s, last_s + 1)):
            k = a_slot.get(s)
            if k is None:
                continue                      # silence: audio lost this slot
            seg = pcm[k * SPF:(k + 1) * SPF]
            if len(seg) == SPF:
                out_a[i * SPF:(i + 1) * SPF] = seg
                hit += 1
        print(f"  audio placed per slot: {hit}/{n_slots} slots filled "
              f"({n_slots - hit} silent) across the full slot range")
        tmp_a = os.path.join(d, "_sync_audio.wav")
        ow = _wave.open(tmp_a, "wb")
        ow.setnchannels(pcm.shape[1]); ow.setsampwidth(2)
        ow.setframerate(w.getframerate())
        ow.writeframes(out_a.tobytes()); ow.close()
        p["audio"] = tmp_a
        p["a_skip"] = 0.0                     # already slot-aligned

    # Cut the video lane at a real fragment boundary: init boxes, then from
    # the aligned fragment's byte offset.
    tmp_v = os.path.join(d, "_sync_video.m4s")
    with open(p["video"]["path"], "rb") as src, open(tmp_v, "wb") as dst:
        init_n = p["init_bytes"] or 0
        src.seek(0)
        dst.write(src.read(init_n))
        # `off` in the index is ABSOLUTE -- LiveWriter seeds its byte counter
        # with the init length, so fragment 0 sits at off == init_n. Adding
        # init_n again lands mid-fragment and ffmpeg silently produces a file
        # with no video track at all rather than an error.
        src.seek(p["v_off"])
        while True:
            b = src.read(1 << 20)
            if not b:
                break
            dst.write(b)

    cmd = ["ffmpeg", "-y", "-loglevel", "error", "-i", tmp_v]
    maps = ["-map", "0:v"]
    if "audio" in p:
        cmd += ["-ss", f"{p['a_skip']:.6f}", "-i", p["audio"]]
        maps += ["-map", "1:a"]
    if "subs" in p:
        cmd += ["-itsoffset", f"{p['s_shift']:.6f}", "-i", p["subs"]]
        maps += ["-map", f"{2 if 'audio' in p else 1}:s"]
    # 6-channel wav = 5.1 (L R C LFE Ls Rs is the standard WAV order, chosen
    # by the worker for exactly this reason); give surround a real bitrate
    abr = "384k" if os.path.exists(p.get("audio", "")) and _channels(
        p["audio"]) == 6 else "192k"
    cmd += maps + ["-c:v", "copy", "-c:a", "aac", "-b:a", abr]
    if "subs" in p:
        cmd += ["-c:s", "mov_text", "-metadata:s:s:0", "language=eng"]
    if a.secs:
        cmd += ["-t", str(a.secs)]
    cmd += [out]

    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        print(r.stderr[-1500:], file=sys.stderr)
        return 1
    try:
        os.remove(tmp_v)
    except OSError:
        pass

    pr = subprocess.run(
        ["ffprobe", "-v", "error", "-show_entries",
         "stream=index,codec_type,codec_name,duration", "-of", "csv=p=0", out],
        capture_output=True, text=True)
    print(f"\nwrote {out}")
    for line in pr.stdout.strip().splitlines():
        print("   ", line)

    if a.play:
        subprocess.run(["ffplay", "-hide_banner", "-loglevel", "error",
                        "-autoexit", out])
    return 0


if __name__ == "__main__":
    sys.exit(main())
