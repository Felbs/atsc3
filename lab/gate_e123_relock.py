#!/usr/bin/env python3
"""gate_e123_relock.py -- a lost lock re-locks the SAME way whatever the layout and the timing.

13.20 found that the blind sweep after a lost CTI lock was anchored at
`cti.next_frame` -- wherever the pool happened to be when the asynchronous
FEC verdict arrived -- so the same lossless capture gave six re-locks at six
process counts and two different lane sets from two P=4 runs.  E123 anchors
on stream order instead.  This gate is the referee: the same capture, replayed
at several layouts and twice at one of them, must give IDENTICAL converged
block counts and byte-identical decoded lanes.  It runs the chain as a
subprocess, exactly as a user would.

Run: python lab/gate_e123_relock.py --capture C:/captures/rf25_win_45s_b.cs16 [--layouts 4x2,4x2,8x2,16x2,24x2]
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:16]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--capture", required=True)
    ap.add_argument("--rate", default="6912000")
    ap.add_argument("--layouts", default="4x2,4x2,8x2,16x2,24x2",
                    help="procs x proc-threads per run; repeat one to test run-to-run")
    ap.add_argument("--out", default=os.path.join(ROOT, "e123_out"))
    ap.add_argument("--keep", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1")
    env.setdefault("ATSC3_LDPC_BOXPLUS", "1")
    rows = []
    for i, lay in enumerate(a.layouts.split(",")):
        main_, _, fecpart = lay.lower().partition("+")   # E126: "PxT" or "PxT+FxG"
        p, t = main_.split("x")
        fec_args = []
        if fecpart:
            f_, g_ = fecpart.split("x")
            fec_args = ["--ldm-fec-procs", f_, "--ldm-fec-proc-threads", g_]
        tag = os.path.join(a.out, f"run{i}_{lay}")
        shutil.rmtree(tag, ignore_errors=True)
        cmd = [sys.executable, "-m", "atsc3", "watch", "--capture", a.capture, "--fmt", "cs16",
               "--rate", a.rate, "--accel", "cpu", "--assets", "all", "--plp1", "--player", "none",
               "--ldm-procs", p, "--ldm-proc-threads", t, "--threads", "4", *fec_args,
               "--live-dir", tag, "--json", tag + ".json"]
        t0 = time.time()
        with open(tag + ".log", "w") as lg:
            rc = subprocess.run(cmd, cwd=ROOT, env=env, stdout=lg, stderr=subprocess.STDOUT).returncode
        d = json.load(open(tag + ".json"))
        L = d["ldm"]
        lanes = {os.path.basename(f): (os.path.getsize(f), sha(f)) for f in sorted(glob.glob(os.path.join(tag, "live_*.m4s")))}
        rows.append(dict(layout=lay, rc=rc, x_rt=d["x_real_time"], conv=d["fec_converged"], total=d["fec_total"],
                         state=L["state"], origin=L["origin_frame"], start_row=L["start_row"],
                         lost=L["stats"].get("lock_lost", 0), reacq=L["stats"].get("reacquired_from_window", 0),
                         behind=L["stats"].get("loss_anchor_behind", 0), stalled=L["stats"].get("loss_lag_stalled", 0),
                         lag_wait=L["ms_per_frame"].get("loss_lag_wait", 0.0), lanes=lanes, wall=time.time() - t0))
        r = rows[-1]
        print(f"  {lay:6s} rc={rc} x_rt={r['x_rt']:.3f} conv {r['conv']}/{r['total']} {r['state']} origin {r['origin']} "
              f"start_row {r['start_row']} lost {r['lost']} reacq {r['reacq']} anchor_behind {r['behind']} "
              f"lag_wait {r['lag_wait']:.1f} ms/F  lanes {[(k[5:], v[0]) for k, v in lanes.items()]}  ({r['wall']:.0f} s)",
              flush=True)
    ref = rows[0]
    same_counts = all(r["conv"] == ref["conv"] and r["total"] == ref["total"] for r in rows)
    same_lanes = all(r["lanes"] == ref["lanes"] for r in rows)
    same_lock = all((r["origin"], r["start_row"]) == (ref["origin"], ref["start_row"]) for r in rows)
    no_behind = all(r["behind"] == 0 and r["stalled"] == 0 for r in rows)
    ok = same_counts and same_lanes and same_lock and no_behind and all(r["rc"] == 0 for r in rows)
    print(f"\n  counts identical: {same_counts}; lanes byte-identical: {same_lanes}; "
          f"final lock identical: {same_lock}; anchor never behind / never stalled: {no_behind}")
    print("  PASS: the re-lock is a function of the stream" if ok else "  FAIL")
    if not a.keep and ok:
        for r in rows:
            pass   # keep the artefacts for the record; they are small
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
