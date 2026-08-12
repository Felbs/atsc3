"""M39 -- put the picture and the sound together on the MMTP timeline.

The two streams are carried in MPUs that share a sequence number, and one
MPU is exactly 2.002 s: 120 video frames at 60000/1001 fps, 60 AC-4 frames
at 30000/1001 fps.  That shared number is the only common clock -- the
fragments' own tfdt is 0 in every MPU, because each MPU is an independent
ISOBMFF file that restarts its timeline.

    slots 225641739 .. 225641797            59 slots
    video   57 present   missing 782, 797
    audio   58 present   missing 797
    subs    59 present   missing none

797 is absent from both media assets, so the programme simply ends there.
782 is a real hole in the middle of the picture, at t = 86.1 s.

HOW THE HOLE IS HANDLED
-----------------------
We drop the audio for slot 782 rather than inventing 2 s of filler video.
Both streams then jump at the same instant, which is what actually
happened on air, and the HEVC is passed through untouched -- no re-encode,
no generated frames presented as received ones.  The alternative (freeze
the picture, keep the audio) is what a consumer receiver does for
concealment; it is a nicer viewing experience and a worse record of what
we received, so it lives behind --conceal.
"""
import os, re, sys, wave, subprocess, collections
import numpy as np

LAB = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(LAB, "m7_out")

VIDEO = os.path.join(OUT, "rf33_video.mp4")
AUDIO = os.path.join(LAB, "tv_audio_48k.wav")

FPS_NUM, FPS_DEN = 60000, 1001          # 59.94
V_PER_MPU = 120
A_PER_MPU = 60
MPU_SECONDS = V_PER_MPU * FPS_DEN / FPS_NUM      # 2.002
OUT_FS = 48000
A_SAMPLES_PER_MPU = 96096               # 60 frames * 1601.6


def slots():
    """Which MPU sequence numbers arrived, per asset."""
    got = collections.defaultdict(set)
    for f in os.listdir(OUT):
        m = re.match(r"mpu_.*_(pid\d+)_(\d+)\.seg$", f)
        if m:
            got[m.group(1)].add(int(m.group(2)))
    return got


def plan():
    got = slots()
    v, a = got["pid12"], got["pid13"]
    lo = min(min(v), min(a))
    hi = max(max(v), max(a))
    common_end = max(s for s in range(lo, hi + 1) if s in v and s in a)
    keep = [s for s in range(lo, common_end + 1)]
    v_missing = [s for s in keep if s not in v]
    a_missing = [s for s in keep if s not in a]
    return lo, keep, v_missing, a_missing


def main():
    conceal = "--conceal" in sys.argv
    lo, keep, v_missing, a_missing = plan()
    print(f"timeline slots {keep[0]} .. {keep[-1]}   {len(keep)} slots"
          f"   {len(keep) * MPU_SECONDS:.3f} s")
    print(f"  video holes {v_missing}   audio holes {a_missing}")
    if a_missing:
        print("  audio holes are not handled by this path -- stopping")
        return 1

    w = wave.open(AUDIO)
    n, ch = w.getnframes(), w.getnchannels()
    d = np.frombuffer(w.readframes(n), "<i2").reshape(n, ch)
    w.close()
    print(f"\naudio in  {n} samples ({n / OUT_FS:.3f} s), {ch} ch")

    n_mpu_audio = n // A_SAMPLES_PER_MPU
    print(f"          = {n_mpu_audio} MPUs of {A_SAMPLES_PER_MPU} samples"
          f"   remainder {n % A_SAMPLES_PER_MPU}")

    if conceal:
        print("\n--conceal: keeping all audio, video will be frozen at the hole")
        cut = d
    else:
        drop = sorted(keep.index(s) for s in v_missing)
        print(f"\ndropping audio for slot(s) {v_missing} "
              f"= MPU index {drop} = t {[i * MPU_SECONDS for i in drop]}")
        mask = np.ones(len(d), bool)
        for i in drop:
            s0 = i * A_SAMPLES_PER_MPU
            mask[s0:s0 + A_SAMPLES_PER_MPU] = False
        cut = d[mask]

    print(f"audio out {len(cut)} samples ({len(cut) / OUT_FS:.4f} s)")

    # what the video is
    pr = subprocess.run(
        ["ffprobe", "-v", "error", "-select_streams", "v",
         "-count_frames", "-show_entries",
         "stream=nb_read_frames,r_frame_rate", "-of", "csv=p=0", VIDEO],
        capture_output=True, text=True)
    fields = pr.stdout.strip().split(",")
    nv = int([x for x in fields if x.isdigit()][0])
    v_dur = nv * FPS_DEN / FPS_NUM
    print(f"video in  {nv} frames ({v_dur:.4f} s)")

    # --- GATE: the two must agree to well under one frame ------------
    a_dur = len(cut) / OUT_FS
    delta = a_dur - v_dur
    frame = FPS_DEN / FPS_NUM
    print()
    print(f"GATE 1  audio {a_dur:.4f} s   video {v_dur:.4f} s")
    print(f"        difference {1000 * delta:+.2f} ms = "
          f"{delta / frame:+.3f} video frames"
          f"      {'PASS' if abs(delta) < frame else 'FAIL'}")
    if abs(delta) >= frame:
        print("        refusing to mux streams that do not line up")
        return 1

    tmp = os.path.join(LAB, "_mux_audio.wav")
    o = wave.open(tmp, "wb")
    o.setnchannels(ch); o.setsampwidth(2); o.setframerate(OUT_FS)
    o.writeframes(cut.tobytes()); o.close()

    dst = os.path.join(LAB, "rf33_tv.mp4")
    cmd = ["ffmpeg", "-y", "-loglevel", "error",
           "-i", VIDEO, "-i", tmp,
           "-c:v", "copy",                     # HEVC untouched
           "-c:a", "aac", "-b:a", "192k",
           "-shortest", dst]
    print(f"\nmuxing -> {os.path.basename(dst)}  (video copied, not re-encoded)")
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        print(r.stderr[-2000:])
        return 1
    os.remove(tmp)

    pr = subprocess.run(
        ["ffprobe", "-v", "error", "-show_entries",
         "stream=codec_type,codec_name,duration", "-of", "csv=p=0", dst],
        capture_output=True, text=True)
    print(pr.stdout.strip())
    print(f"\nwrote {dst}   {os.path.getsize(dst) / 1e6:.1f} MB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
