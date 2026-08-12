#!/usr/bin/env python3
"""M3 Step 2 (part 3) -- replicate the L1-Basic demodulation across FRAMES.

The point of this test is that it needs NO spec tables beyond Step 2 and no
LDPC decoder, yet it is very hard to pass by accident.

L1-Basic describes the frame's structure: FFT size, guard interval, pilot
pattern, subframe count, L1-Detail size.  Those do not change from one frame to
the next in a running transmitter.  So if the demodulation is real, the 484
QPSK cells -- 968 hard bits -- must come out very nearly IDENTICAL in frames
that are hundreds of milliseconds apart and share no samples.

Two independent noise clouds agree on ~50% of bits.  Two independent
observations of the same signalling agree on ~100%.  There is no way to land in
between by luck, and nothing in the chain is fitted to make it happen: each
frame is demodulated from scratch, with its own bootstrap anchor, its own
timing search and its own channel estimate.

CONTROL: the same comparison is run against cells drawn from the L1-DETAIL
region, which carries interleaved payload signalling and is NOT expected to be
frame-invariant in the same way, and against a bit-shuffled copy.

Usage:
    python m3_replicate.py hit_rf33.cs16 --rate 8e6 --frames 8
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m3_spec as S                                            # noqa: E402
import m3_freqint as FI                                        # noqa: E402
from m3_preamble import analyse, QPSK, mer_db, m4_stat          # noqa: E402


def hard_bits(z):
    """QPSK hard decision -> 2 bits per cell (sign of I, sign of Q)."""
    return np.concatenate([(z.real < 0).astype(np.uint8)[:, None],
                           (z.imag < 0).astype(np.uint8)[:, None]],
                          axis=1).ravel()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("capture")
    ap.add_argument("--rate", type=float, required=True)
    ap.add_argument("--fmt")
    ap.add_argument("--frames", type=int, default=8)
    ap.add_argument("--step", type=float, default=0.7,
                    help="seconds between frames sampled (frame period is "
                         "~0.247 s, so this lands on different frames)")
    ap.add_argument("--json")
    a = ap.parse_args()

    path = a.capture
    if not os.path.isabs(path):
        path = os.path.join(HERE, path)

    print(f"A/322 spec tables: "
          f"{sum(1 for _, g, _ in S.verify(verbose=False) if g)}/13 checks pass")
    print(f"\nDemodulating {a.frames} frames of {os.path.basename(path)}, "
          f"{a.step}s apart -- disjoint samples, independent bootstrap anchors\n")

    frames = []
    for i in range(a.frames):
        t = i * a.step
        try:
            rep = {}
            rep, z, Y, geo = analyse(path, a.rate, a.fmt, report=rep,
                                     start_sec=t, quiet=True)
        except Exception as e:                                  # noqa: BLE001
            print(f"  frame @ {t:5.2f}s   SKIP ({type(e).__name__}: {e})")
            continue
        (lo, n, pilot, cp, data, shift, cred, nfft, gi, dx, mode) = geo
        ncells = S.L1_BASIC_CELLS_PRINTED[mode]
        x = FI.deinterleave(z, nfft, 0, direction="forward", toggle="i")
        l1b = x[:ncells]
        mer = mer_db(l1b, QPSK)
        m4 = m4_stat(l1b)
        frames.append({
            "t": t, "mer": mer, "m4": abs(m4),
            "coh": rep["pilot_coherence"], "psnr": rep["pilot_snr_db"],
            "cred": cred, "shift": shift,
            "bits": hard_bits(l1b),
            "detail_bits": hard_bits(x[ncells:ncells + ncells]),
            "mode": mode, "ncells": ncells,
        })
        print(f"  frame @ {t:5.2f}s   Cred {cred} shift {shift:+d}   "
              f"pilot coh {rep['pilot_coherence']:.4f}   "
              f"pilot SNR {rep['pilot_snr_db']:5.1f} dB   "
              f"L1-Basic QPSK MER {mer:6.2f} dB   |E[z^4]| {abs(m4):.3f}")

    if len(frames) < 2:
        raise SystemExit("need at least 2 good frames")

    ref = frames[0]
    nb = len(ref["bits"])
    print(f"\n  --- bit agreement against frame 0 ({nb} hard bits "
          f"= {ref['ncells']} QPSK cells) ---")
    print("   frame      L1-Basic agreement      L1-Detail region (control)"
          "   shuffled (control)")
    rng = np.random.default_rng(7)
    rows = []
    for f in frames[1:]:
        agree = float(np.mean(f["bits"] == ref["bits"]))
        dctl = float(np.mean(f["detail_bits"] == ref["detail_bits"]))
        sctl = float(np.mean(rng.permutation(f["bits"]) == ref["bits"]))
        rows.append({"t": f["t"], "agree": agree, "detail": dctl,
                     "shuffled": sctl})
        print(f"   @{f['t']:5.2f}s      {agree*100:6.2f}%   "
              f"({int(agree*nb):4d}/{nb})        {dctl*100:6.2f}%"
              f"                    {sctl*100:6.2f}%")

    ag = np.array([r["agree"] for r in rows])
    print(f"\n  L1-Basic bit agreement: mean {ag.mean()*100:.2f}%  "
          f"min {ag.min()*100:.2f}%  max {ag.max()*100:.2f}%")
    print(f"  chance level is 50.00%; a shared-signalling readout is ~100%")

    # majority vote across frames -> the cleanest estimate of the true bits
    B = np.array([f["bits"] for f in frames])
    maj = (B.mean(axis=0) > 0.5).astype(np.uint8)
    per = [float(np.mean(b == maj)) for b in B]
    print(f"\n  agreement with the {len(frames)}-frame majority vote: "
          + " ".join(f"{p*100:.1f}%" for p in per))

    out = {
        "capture": os.path.basename(path),
        "frames": [{k: v for k, v in f.items()
                    if k not in ("bits", "detail_bits")} for f in frames],
        "pairwise": rows,
        "mean_agreement": float(ag.mean()),
        "majority_bits": "".join(str(int(b)) for b in maj),
    }
    dest = a.json or os.path.join(
        HERE, "m3_replicate_" + os.path.splitext(os.path.basename(path))[0] + ".json")
    with open(dest, "w") as fh:
        json.dump(out, fh, indent=2, default=float)
    print(f"\n  wrote {dest}")


if __name__ == "__main__":
    main()
