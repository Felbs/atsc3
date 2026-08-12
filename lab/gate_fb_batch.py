#!/usr/bin/env python3
"""GATE: batched 6-channel filterbank -- bit-identical, and actually faster.

The live 5.1 worker died on speed: 0.62x real time even with six threads,
because the synthesis loop is many small windows and the GIL never lets the
threads overlap (E32). The batch path amortises the Python overhead by
running all six channels through ONE dct call per window.

Three requirements, each falsifiable:
  1. BIT-IDENTICAL to the sequential path on real broadcast frames --
     np.array_equal, not a tolerance (the fleet's gate_lib discipline).
  2. Bit-identical ACROSS chunk boundaries with carried state (the live
     worker decodes in passes; a seam bug would be inaudible-in-gates and
     audible-on-air).
  3. >= 1.3x the sequential rate on the same work, or live 5.1 stays NO-GO.

NEGATIVE CONTROL: perturb one spectral line by one ulp in one channel --
equality must FAIL (proves the comparison has teeth).

Run:  python lab/gate_fb_batch.py
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "tools"))

import m30_filterbank as FB                                      # noqa: E402
from m42_ac4_stream import Ac4Stream                             # noqa: E402
from atsc3_audio import LaneReader                               # noqa: E402


def main():
    lane = os.path.join(os.path.dirname(HERE),
                        "data", "catchup", "live_audio_pid13.m4s")
    r = LaneReader(lane)
    frames = r.read_new()[:3000]                 # 100 s of real air
    dec = Ac4Stream()
    wins, lens, _, cnts = dec.decode_frames(frames)
    chans = list(Ac4Stream.CHANS)
    div = sum(1 for i in range(len(frames))
              if len({tuple(lens[c][sum(cnts[c][:i]):sum(cnts[c][:i + 1])])
                      for c in chans}) > 1)
    print(f"{len(frames)} frames decoded; {div} have divergent framings "
          f"(the grouped path must handle BOTH kinds)")

    # ---- 1. bit-identical, whole run, REAL divergent data -------------
    t0 = time.perf_counter()
    seq = {c: FB.synthesise(wins[c], lens[c], 1536) for c in chans}
    t_seq = time.perf_counter() - t0
    t0 = time.perf_counter()
    bat = FB.synthesise_frames(wins, lens, cnts, 1536)
    t_bat = time.perf_counter() - t0
    ident = all(np.array_equal(seq[c], bat[c]) for c in chans)
    print(f"whole-run identical: {ident}   "
          f"sequential {t_seq:.2f}s vs grouped {t_bat:.2f}s "
          f"= {t_seq / t_bat:.2f}x")

    # ---- 2. seam: chunked with carried state vs one-shot --------------
    hf = len(frames) // 2
    def cut(d, lo, hi):
        w = {}; l = {}; k = {}
        for c in chans:
            a0 = sum(cnts[c][:lo]); a1 = sum(cnts[c][:hi])
            w[c] = wins[c][a0:a1]; l[c] = lens[c][a0:a1]
            k[c] = cnts[c][lo:hi]
        return w, l, k
    w1, l1, k1 = cut(None, 0, hf)
    w2, l2, k2 = cut(None, hf, len(frames))
    p1, s1 = FB.synthesise_frames(w1, l1, k1, 1536, return_state=True)
    p2 = FB.synthesise_frames(w2, l2, k2, 1536, states=s1)
    seam = all(np.array_equal(np.concatenate([p1[c], p2[c]]), bat[c])
               for c in chans)
    print(f"chunked-with-state identical to one-shot: {seam}")

    # ---- 3. negative control ------------------------------------------
    wins2 = {c: [np.array(w, copy=True) for w in wins[c]] for c in chans}
    wins2["C"][7][3] = np.nextafter(wins2["C"][7][3], np.inf)
    bad = FB.synthesise_frames(wins2, lens, cnts, 1536)
    caught = not np.array_equal(bad["C"], bat["C"])
    print(f"negative control (1-ulp perturbation) caught: {caught}")

    speed_ok = (t_seq / t_bat) >= 1.3
    ok = ident and seam and caught
    print()
    if ok:
        verdict = "GO" if speed_ok else "NO-GO (correct but not fast enough)"
        print(f"GATE PASS -- bit-identical incl. seams; live 5.1 verdict: "
              f"{verdict} at {t_seq / t_bat:.2f}x")
        return 0
    print("GATE FAIL")
    return 1


if __name__ == "__main__":
    sys.exit(main())
