#!/usr/bin/env python3
"""size_ldpc_abort2.py -- sizing E125, the two-level abort, on a banked LDPC record.

With the per-code ladder (E121) a Block that fails rung 1 is re-decoded on
rungs 2 and 3, and a HOPELESS Block pays the E112 abort (sweep 8) on each --
~24 sweeps where E112 alone spent 8; at 27% convergence that is most of the
LDPC (13.21: 658 ms a Frame).  13.20's E121 record showed the two populations
are separable: rung-1 aborts with final syndrome 9,331-10,502 were rescued by
rung 2, truly hopeless Blocks sit flat near 0.39 x nc.

So, from the per-call record (ATSC3_LDPC_DUMP: every call's conv / its / bad),
ladders are re-assembled (rung k+1's rows are rung k's failing rows, in order)
and every rung-1 ABORTED Block is tagged with: its rung-1 final syndrome, and
whether ANY later rung converged it.  For a grid of T2 (fraction of nc):
  dropped  = rung-1 aborts with syndrome >= T2 x nc (would get no later rung)
  lost     = dropped Blocks a later rung DID converge   (must be 0)
  saved    = the later-rung sweeps those dropped Blocks cost
  retried  = rung-1 aborts below T2 (still climb the ladder)
The smallest zero-loss T2 and its saving are the answer.

Run: python lab/size_ldpc_abort2.py --dump DIR [--ladder 0.82,0.85,1.0] [--nc 25920] [--abort-k 8] [--abort-t 0.36]
"""
from __future__ import annotations

import argparse
import glob
import os
import re
import sys

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", required=True)
    ap.add_argument("--ladder", default="0.82,0.85,1.0")
    ap.add_argument("--nc", type=int, default=25920, help="parity checks of the code (64800 x 6/15)")
    ap.add_argument("--abort-k", type=int, default=8)
    ap.add_argument("--abort-t", type=float, default=0.36)
    a = ap.parse_args()
    ladder = [float(x) for x in a.ladder.split(",")]
    files = sorted(glob.glob(os.path.join(a.dump, "*.npz")),
                   key=lambda f: (re.search(r"ldpc_(\d+)_(\d+)", f).group(1),
                                  int(re.search(r"ldpc_(\d+)_(\d+)", f).group(2))))
    recs = []
    for f in files:
        z = np.load(f)
        recs.append(dict(pid=re.search(r"ldpc_(\d+)_", f).group(1), a=float(z["alpha"]),
                         conv=z["conv"].astype(bool), its=z["its"].astype(int), bad=z["bad"].astype(int)))
    # re-assemble ladders: a new one starts at the first rung's alpha; later rungs follow on the same pid
    ladders, cur = [], None
    for r in recs:
        if abs(r["a"] - ladder[0]) < 1e-6:
            cur = [r]; ladders.append(cur)
        elif cur is not None and cur[-1]["pid"] == r["pid"] and len(cur) < len(ladder) \
                and abs(r["a"] - ladder[len(cur)]) < 1e-6 and len(r["conv"]) == int((~cur[-1]["conv"]).sum()):
            cur.append(r)
    T1 = int(round(a.abort_t * a.nc))
    tot_sweeps = sum(int(r["its"].sum()) for r in recs)
    n_blocks = sum(len(l[0]["conv"]) for l in ladders)
    conv_total = sum(int(sum(r["conv"].sum() for r in l)) for l in ladders)
    # per rung-1 failed Block: bad1, aborted?, rescued later?, later sweeps
    bad1, aborted, rescued, later_sw, conv1 = [], [], [], [], 0
    for l in ladders:
        r1 = l[0]
        conv1 += int(r1["conv"].sum())
        fail1 = np.flatnonzero(~r1["conv"])
        resc = np.zeros(len(fail1), bool); lsw = np.zeros(len(fail1), int)
        if len(l) > 1:
            r2 = l[1]
            resc |= r2["conv"]; lsw += r2["its"]
            if len(l) > 2:
                r3 = l[2]
                idx3 = np.flatnonzero(~r2["conv"])          # rung-3 rows map to rung-2's failing rows
                resc[idx3] |= r3["conv"]; lsw[idx3] += r3["its"]
        bad1 += list(r1["bad"][fail1]); aborted += list((r1["its"][fail1] == a.abort_k) & (r1["bad"][fail1] >= T1))
        rescued += list(resc); later_sw += list(lsw)
    bad1 = np.array(bad1); aborted = np.array(aborted, bool); rescued = np.array(rescued, bool); later_sw = np.array(later_sw)
    print(f"{len(recs)} calls -> {len(ladders)} ladders {ladder}, {n_blocks} Blocks, converged {conv_total} "
          f"({100*conv_total/max(1,n_blocks):.1f}%): rung 1 {conv1}, later rungs {conv_total-conv1}; "
          f"lane-sweeps total {tot_sweeps}")
    print(f"rung-1 failures {len(bad1)}: aborted at sweep {a.abort_k} with syndrome >= {T1} (0.36 nc): {aborted.sum()}; "
          f"of those rescued by a later rung: {int((aborted & rescued).sum())}; rescued among non-aborted failures: {int((~aborted & rescued).sum())}")
    if aborted.any():
        ra = bad1[aborted & rescued]
        print(f"  rung-1 syndrome of RESCUED aborts: n={len(ra)} min/p50/max = "
              f"{(ra.min() if len(ra) else '-')}/{(int(np.median(ra)) if len(ra) else '-')}/{(ra.max() if len(ra) else '-')}  "
              f"({ra.max()/a.nc if len(ra) else 0:.4f} x nc at the top)")
        ha = bad1[aborted & ~rescued]
        print(f"  rung-1 syndrome of aborts NEVER rescued: n={len(ha)} p5/p50/p95 = "
              f"{np.percentile(ha,[5,50,95]).astype(int).tolist() if len(ha) else '-'}")
    print(f"\n  {'T2 x nc':>8s} {'T2':>6s} {'dropped once':>12s} {'still retried':>13s} {'lost':>5s} {'sweeps saved':>12s} {'of all LDPC':>11s}")
    best = None
    for frac in np.arange(0.36, 0.421, 0.005):
        T2 = int(round(frac * a.nc))
        drop = aborted & (bad1 >= T2)
        lost = int((drop & rescued).sum())
        saved = int(later_sw[drop].sum())
        print(f"  {T2:8d} {frac:6.3f} {int(drop.sum()):12d} {int((aborted & ~drop).sum()):13d} {lost:5d} {saved:12d} {100*saved/max(1,tot_sweeps):10.1f}%")
        if lost == 0 and best is None:
            best = (frac, T2, int(drop.sum()), saved)
    if best:
        print(f"\n  smallest zero-loss T2 = {best[0]:.3f} x nc ({best[1]}): {best[2]} Blocks dropped after rung 1, "
              f"{best[3]} lane-sweeps saved = {100*best[3]/max(1,tot_sweeps):.1f}% of all LDPC sweeps on this record")
    else:
        print("\n  no zero-loss T2 in the grid")
    return 0


if __name__ == "__main__":
    sys.exit(main())
