"""M38 -- AC-4 output resampling to the true 48 kHz presentation rate.

THE BUG THIS FIXES
------------------
M29/M37 render 1536 samples per AC-4 frame and label the result 48 kHz.
1536 is the MDCT length -- the INTERNAL frame length.  ETSI TS 103 190-1
Table 83 specifies, for `frame_rate_index = 3` (29.97 fps) at 48 kHz, a
decoder output resampling ratio of

    1001/1000 x 25/24  =  1001/960  =  1.04270833...

so that one frame of 1536 internal samples becomes 1601.6 output samples
occupying 1/29.97 s = 33.367 ms.  Without it the audio plays 4.27 % fast:
0.72 of a semitone sharp, and ~4.8 s of drift across a two-minute programme.
Inaudible on its own; fatal the moment you put it next to picture.

INDEPENDENT CONFIRMATION (does not rely on reading Table 83 correctly)
---------------------------------------------------------------------
Every MMTP MPU on this service carries 120 video frames and 60 AC-4 frames.
Video runs at exactly 60000/1001 fps, so one MPU spans 120 x 1001/60000 =
2.002 s, and therefore one AC-4 frame spans 2.002/60 = 33.367 ms = 1601.6
samples at 48 kHz.  Two derivations, one from the audio spec and one from
the transport layer, agreeing to the digit.

WHAT THIS MEANS FOR THE A-SPX BANDWIDTH CLAIM  (a correction)
--------------------------------------------------------------
Resampling does NOT move physical frequencies -- it only fixes the rate
label -- so every frequency we quoted for this decoder was in the INTERNAL
domain, computed as if the render were 48 kHz.  It was not.  True internal
fs = 1536 / 0.033367 = 46033.97 Hz, i.e. 0.9590 of what we assumed, so
every figure scales by that:

    quoted (internal-normalised)      actual (physical)
    core/A-SPX crossover  12000 Hz    11508 Hz    (line 768 of 1536)
    A-SPX stop            21000 Hz    20140 Hz    (QMF subband 56 of 64)
    QMF band spacing        375 Hz      359.6 Hz

The three-way agreement between the Annex B band table, the A-SPX header
and the QMF measurement still holds -- those subsystems agreed with each
other, and they still do.  What was wrong was the absolute Hz label we
hung on all three, because it came from the sample rate we assumed rather
than the one Table 83 defines.  A shared wrong constant is invisible to
any amount of cross-checking between users of that constant.
"""
import numpy as np
from scipy.signal import resample_poly

UP, DOWN = 1001, 960            # Table 83: 1001/1000 x 25/24
INTERNAL_FRAME = 1536           # MDCT length
OUTPUT_FRAME = 1601.6           # 1536 * 1001/960
FRAME_SECONDS = 1001.0 / 30000.0    # 1/29.97 = 33.367 ms
OUT_FS = 48000


def internal_fs():
    """The rate the renderer's samples are actually at."""
    return INTERNAL_FRAME / FRAME_SECONDS


def resample(x, up=UP, down=DOWN):
    """Rational resample to the presentation rate.

    `x` is (n,) or (n, ch).  resample_poly's polyphase FIR is applied per
    column; the anti-imaging filter is designed for the 1/1001 band so no
    aliasing is introduced (output Nyquist 24000 Hz vs 21000 Hz of content).
    """
    if x.ndim == 1:
        return resample_poly(x, up, down)
    return np.column_stack([resample_poly(x[:, c], up, down)
                            for c in range(x.shape[1])])


def expected_length(n_frames):
    """Exact output sample count for a whole number of AC-4 frames."""
    # 1536 * 1001/960 = 1601.6 -- integral only over groups of 5 frames,
    # so compute on the total rather than per frame.
    total = n_frames * INTERNAL_FRAME
    return total * UP // DOWN


if __name__ == "__main__":
    import wave, sys, hashlib

    src = sys.argv[1] if len(sys.argv) > 1 else "tv_audio_hf_full.wav"
    dst = sys.argv[2] if len(sys.argv) > 2 else "tv_audio_48k.wav"

    w = wave.open(src)
    n, ch, fs = w.getnframes(), w.getnchannels(), w.getframerate()
    d = np.frombuffer(w.readframes(n), "<i2").reshape(n, ch).astype(np.float64)
    w.close()
    print(f"in   {src}   {n} samples x {ch} ch, labelled {fs} Hz")

    n_frames = n // INTERNAL_FRAME
    print(f"     = {n_frames} AC-4 frames of {INTERNAL_FRAME} internal samples")
    print(f"     true internal rate {internal_fs():.2f} Hz")
    print(f"     true duration      {n_frames * FRAME_SECONDS:.3f} s")
    print(f"     as-written         {n / OUT_FS:.3f} s   "
          f"({100 * (OUT_FS / internal_fs() - 1):+.2f} % fast)")
    print()

    y = resample(d / 32768.0)
    want = expected_length(n_frames)
    print(f"out  {y.shape[0]} samples   expected {want}   "
          f"{'PASS' if y.shape[0] == want else 'MISMATCH'}")
    print(f"     duration {y.shape[0] / OUT_FS:.3f} s")

    # --- gates -------------------------------------------------------
    print()
    ok = True

    g1 = y.shape[0] == want
    print(f"GATE 1  exact output length {y.shape[0]} == {want}"
          f"                {'PASS' if g1 else 'FAIL'}")
    ok &= g1

    dur = y.shape[0] / OUT_FS
    g2 = abs(dur - n_frames * FRAME_SECONDS) < 1e-3
    print(f"GATE 2  duration {dur:.4f} s == {n_frames * FRAME_SECONDS:.4f} s "
          f"expected            {'PASS' if g2 else 'FAIL'}")
    ok &= g2

    r_in = np.sqrt((d / 32768.0) ** 2).mean()
    r_out = np.sqrt(y ** 2).mean()
    g3 = abs(r_out / r_in - 1.0) < 0.02
    print(f"GATE 3  rms preserved  in {r_in:.5f}  out {r_out:.5f}  "
          f"ratio {r_out / r_in:.4f}   {'PASS' if g3 else 'FAIL'}")
    ok &= g3

    # band edge should land at 21000 Hz at the OUTPUT rate
    seg = y[:, 0][:1 << 20]
    P = np.abs(np.fft.rfft(seg * np.hanning(len(seg)))) ** 2
    f = np.fft.rfftfreq(len(seg), 1.0 / OUT_FS)
    tot = P.sum()
    b = lambda a, c: 100 * P[(f >= a) & (f < c)].sum() / tot
    print(f"GATE 4  band edge at the output rate:")
    for a, c in ((6000, 12000), (12000, 20000), (20000, 21500),
                 (21500, 24000)):
        print(f"          {a:5d}-{c:5d} Hz   {b(a, c):8.5f} %")
    g4 = b(21500, 24000) < 0.001 and b(12000, 20000) > 0.001
    print(f"        content present below 21 kHz, none above  "
          f"          {'PASS' if g4 else 'FAIL'}")
    ok &= g4

    print()
    print("ALL GATES PASS" if ok else "GATES FAILED")

    pk = np.abs(y).max()
    y16 = np.clip(y / max(pk, 1e-12) * 0.98, -1, 1)
    y16 = (y16 * 32767).astype("<i2")
    out = wave.open(dst, "wb")
    out.setnchannels(y.shape[1]); out.setsampwidth(2); out.setframerate(OUT_FS)
    out.writeframes(y16.tobytes()); out.close()
    print(f"\nwrote {dst}   {y.shape[0] / OUT_FS:.3f} s @ {OUT_FS} Hz")
    print("sha256", hashlib.sha256(open(dst, "rb").read()).hexdigest()[:16])
