#!/usr/bin/env python3
"""size_ldpc_rungs.py -- what the later alpha rungs cost and what they rescue.

Companion to size_ldpc_abort2.py (same ladder re-assembly from the per-call
ATSC3_LDPC_DUMP record: rung k+1's rows are rung k's failing rows, in order).
It answers the two questions E125's near-miss gate stands on, per record:
  * sweeps per rung -- what rungs 2 and 3 cost of all LDPC work;
  * the rung-1 FINAL SYNDROME of every Block a later rung rescued (aborted or
    not), and how many of them lie above the gate (--near x nc): those are the
    Blocks the gate would lose.  13.22 found <= 20 on both 32-core records.
Run: python lab/size_ldpc_rungs.py --dump DIR [--ladder 0.82,0.85,1.0] [--nc 25920] [--near 0.02]
"""
import argparse, glob, os, re
import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", required=True)
    ap.add_argument("--ladder", default="0.82,0.85,1.0")
    ap.add_argument("--nc", type=int, default=25920)
    ap.add_argument("--near", type=float, default=0.02)
    a = ap.parse_args()
    ladder = [float(x) for x in a.ladder.split(",")]
    def key(f):
        m = re.search(r"ldpc_(\d+)_(\d+)\.npz$", f); return (int(m.group(1)), int(m.group(2)))
    files = sorted(glob.glob(os.path.join(a.dump, "*.npz")), key=key)
    recs = []
    for f in files:
        z = np.load(f)
        recs.append(dict(pid=key(f)[0], a=float(z["alpha"]), conv=z["conv"].astype(bool),
                         its=z["its"].astype(int), bad=z["bad"].astype(int)))
    ladders, cur = [], None
    for r in recs:
        if abs(r["a"] - ladder[0]) < 1e-6:
            cur = [r]; ladders.append(cur)
        elif (cur is not None and cur[-1]["pid"] == r["pid"] and len(cur) < len(ladder)
              and abs(r["a"] - ladder[len(cur)]) < 1e-6 and len(r["conv"]) == int((~cur[-1]["conv"]).sum())):
            cur.append(r)
        else:
            cur = None
    sw = np.zeros(len(ladder), int); conv = np.zeros(len(ladder), int); nblk = np.zeros(len(ladder), int)
    bad1_resc, bad1_fail = [], []
    for l in ladders:
        for k, r in enumerate(l):
            sw[k] += int(r["its"].sum()); conv[k] += int(r["conv"].sum()); nblk[k] += len(r["conv"])
        r1 = l[0]; fail1 = np.flatnonzero(~r1["conv"])
        resc = np.zeros(len(fail1), bool)
        if len(l) > 1:
            r2 = l[1]; resc |= r2["conv"]
            if len(l) > 2:
                idx3 = np.flatnonzero(~r2["conv"]); resc[idx3] |= l[2]["conv"]
        bad1_resc += list(r1["bad"][fail1][resc]); bad1_fail += list(r1["bad"][fail1][~resc])
    bad1_resc = np.array(bad1_resc); bad1_fail = np.array(bad1_fail); tot = int(sw.sum())
    print(f"{len(recs)} calls -> {len(ladders)} ladders {ladder}; Blocks entering each rung {nblk.tolist()}; "
          f"converged per rung {conv.tolist()}; nc {a.nc}")
    print("  lane-sweeps per rung: " + ", ".join(f"rung {k+1} {s} ({100*s/max(1,tot):.1f}%)" for k, s in enumerate(sw))
          + f"; total {tot}; rungs 2+ = {100*sw[1:].sum()/max(1,tot):.1f}% of all LDPC sweeps")
    T = int(round(a.near * a.nc))
    if len(bad1_resc):
        q = np.percentile(bad1_resc, [0, 50, 95, 100]).astype(int)
        above = int((bad1_resc >= T).sum())
        print(f"  rung-1 final syndrome of the {len(bad1_resc)} Blocks a later rung RESCUED: "
              f"min/p50/p95/max = {q[0]}/{q[1]}/{q[2]}/{q[3]} ({q[3]/a.nc:.4f} x nc at the top); "
              f"at or above the gate {a.near} x nc = {T}: {above}  <- Blocks the gate would LOSE here")
    else:
        print("  no Block was rescued by a later rung")
    if len(bad1_fail):
        q = np.percentile(bad1_fail, [5, 50, 95]).astype(int)
        print(f"  rung-1 final syndrome of the {len(bad1_fail)} rung-1 failures NEVER rescued: p5/p50/p95 = {q[0]}/{q[1]}/{q[2]}; "
              f"below the gate (would still be retried): {int((bad1_fail < T).sum())}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
