#!/usr/bin/env python3
"""e86_snr.py -- median dummy-cell channel SNR of a banked capture.

WHY THIS AND NOT FEC%.  E86 tried to choose a front-end setting for RF25 by
decoding one 45 s capture per setting and comparing "% of FEC Blocks
converged".  That comparison is not sound on this carrier: while the captures
were being taken the HDHomeRun referee on the same splitter read snq 65, 58,
61, 60 and seq 100, 21, 0, 100.  A carrier moving that much between captures
will happily rank a setting first or last on drift alone, and FEC% is the
WORST possible instrument for it -- near threshold it is a cliff, so it turns
a fraction of a dB of drift into a 0%-or-100% verdict and reports the noise
as a decision.

So measure the CHANNEL, not the cliff.  M16.frame_snr_db reads SNR off the
known +-1 dummy cells -- no constellation hypothesis, ~microseconds, and it
degrades smoothly instead of falling off an edge.  It is also the instrument
E76 used (snr_db_median), so numbers from this script are directly comparable
to lab/carrier_gain.json's 18.895 dB for RF25 and to the composite threshold.

Per frame it reports one number, so a run yields a DISTRIBUTION: the median
is the setting's score and the spread is how much the carrier moved while we
looked at it.  If the spread is comparable to the gap between two settings,
the comparison has not earned a winner and this prints that verdict rather
than a ranking.

    python lab/e86_snr.py CAPTURE [--frames 60] [--json OUT]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m11_stream as ST          # noqa: E402
import m16_margin as M16         # noqa: E402
import m44_ldm as M44            # noqa: E402
from e82_ldm_run import _preamble_coh_raw   # noqa: E402


def measure(capture, rate, frames, fe_threads=4):
    from concurrent.futures import ThreadPoolExecutor
    ex = ThreadPoolExecutor(max_workers=fe_threads)
    state = {"plan": None}

    def plan_cb(w, t0_local, ps=None):
        best = None
        for t in range(t0_local - 40, t0_local + 8):
            c = _preamble_coh_raw(w, t, ps)
            if best is None or c > best[1]:
                best = (t, c)
        L = M44.l1_from_window(w, best[0], ps=ps)
        if L is None or not L.get("ok"):
            return None
        plan = M44.LdmPlan.from_l1_result(L, label="snr")
        state["plan"] = plan
        return plan

    fe = ST.FrontEnd(rate, ex=ex, fast=True, plan_cb=plan_cb)
    src = ST.FileSource(capture, rate, realtime=False)
    src.start()
    vals, agree, n = [], [], 0
    fd = None
    try:
        while n < frames:
            x = src.read()
            if x is None:
                break
            if isinstance(x, bytes):
                fe.reacquire()
                continue
            fe.push(x)
            for idx, w, t0, coh in fe.frames():
                if n >= frames:
                    break
                plan = state["plan"]
                if plan is None:
                    continue
                if fd is None:
                    fd = M44.FrameDemod(plan, threads=fe_threads, ex=ex)
                try:
                    pool, owner = fd.cell_pool(w, t0, {})
                    if len(pool) != plan.pool_pred:
                        continue
                    pool = fd.cpe_fast(pool, owner)
                except Exception:
                    continue
                # dummy cells run from core_total to the end of the pool --
                # the same slice dummy_agreement() uses, so the SNR and the
                # agreement referee are reading the identical cells.
                s = M16.frame_snr_db(pool, plan.core_total)
                if s is not None and np.isfinite(s):
                    vals.append(float(s))
                    a = fd.dummy_agreement(pool)
                    if a is not None:
                        agree.append(float(a))
                n += 1
        state["agree"] = agree
    finally:
        try:
            src.stop()
        except Exception:
            pass
        ex.shutdown(wait=False)
    return vals, state.get("agree", []), state["plan"]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("capture")
    ap.add_argument("--rate", type=float, default=6.912e6)
    ap.add_argument("--frames", type=int, default=60)
    ap.add_argument("--json", default=None)
    a = ap.parse_args(argv)

    t0 = time.time()
    vals, agree, plan = measure(a.capture, a.rate, a.frames)
    if not vals:
        print(f"{os.path.basename(a.capture)}: NO frames yielded an SNR "
              f"-- cannot score this setting (that is a result, not a zero)")
        return 1
    v = np.asarray(vals)
    out = dict(capture=os.path.basename(a.capture), n=len(vals),
               median=float(np.median(v)), mean=float(v.mean()),
               p10=float(np.percentile(v, 10)),
               p90=float(np.percentile(v, 90)),
               spread=float(np.percentile(v, 90) - np.percentile(v, 10)),
               lo=float(v.min()), hi=float(v.max()),
               dummy_agreement=(float(np.median(agree)) if agree else None),
               secs=round(time.time() - t0, 1))
    ag = (f"  dummy-agree {out['dummy_agreement']:.3f}"
          if out["dummy_agreement"] is not None else "")
    print(f"{out['capture']:34s} n={out['n']:3d}  median {out['median']:6.2f} dB"
          f"  p10..p90 {out['p10']:6.2f}..{out['p90']:6.2f}"
          f"  (spread {out['spread']:.2f} dB){ag}")
    if a.json:
        json.dump(out, open(a.json, "w"), indent=1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
