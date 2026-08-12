#!/usr/bin/env python3
"""E49 -- correlate per-block FEC verdicts with the dummy-cell SNR timeline.

Each deinterleaved FEC Block's cells come from input positions
[p, p + Nrows*(Nrows-1)] where p is its first cell; with Nrows = 1024 that is
~0.9 Frame of spread, so a block sees roughly the channel of one Frame.
This prints convergence vs per-frame SNR -- the measured composite threshold
on air -- and the contiguous decoded spans (candidate media windows).
"""
from __future__ import annotations

import argparse
import json

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("prefix")
    ap.add_argument("--plp", type=int, default=0)
    a = ap.parse_args()

    blocks = json.load(open(f"{a.prefix}_plp{a.plp}_blocks.json"))
    diags = json.load(open(f"{a.prefix}_diags.json"))
    frames = {d["frame"]: d for d in diags["diags"]}
    per_frame = blocks["per_frame"]
    ncell = blocks["cells_per_fec"]

    rows = []
    for b in blocks["blockmap"]:
        # centre of the block's input span
        centre = b["cell"] + ncell / 2 + 1024 * 1023 / 2
        fr = int(centre // per_frame)
        d = frames.get(fr, {})
        rows.append(dict(block=b["block"], frame=fr,
                         t=d.get("t"), snr=d.get("snr_db"),
                         conv=b["conv"], bch=b["bch"], unsat=b["unsat"]))

    conv = [r for r in rows if r["conv"]]
    print(f"{len(conv)}/{len(rows)} blocks converged "
          f"({100*len(conv)/max(len(rows),1):.1f}%)")

    # convergence vs SNR bins
    print("\n  SNR bin      blocks   converged")
    edges = [-99, -3, 0, 2, 3, 4, 5, 6, 7, 99]
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = [r for r in rows if r["snr"] is not None and lo <= r["snr"] < hi]
        if not sel:
            continue
        c = sum(r["conv"] for r in sel)
        print(f"  [{lo:3d},{hi:3d})   {len(sel):6d}   {c:6d}  "
              f"({100*c/len(sel):5.1f}%)")

    # contiguous converged spans (in blocks, ~35 blocks/frame)
    print("\n  contiguous converged spans (>= 35 blocks ~ 1 Frame):")
    run = []
    spans = []
    for r in rows:
        if r["conv"]:
            run.append(r)
        else:
            if len(run) >= 35:
                spans.append(run)
            run = []
    if len(run) >= 35:
        spans.append(run)
    for s in spans:
        t0 = s[0]["t"] or 0
        t1 = s[-1]["t"] or 0
        print(f"    blocks {s[0]['block']:5d}..{s[-1]['block']:5d}  "
              f"t {t0:7.2f}..{t1:7.2f}s  ({len(s)} blocks, "
              f"{(t1-t0):.2f}s)")
    if not spans:
        print("    none")

    # threshold estimate: median SNR of converged vs failed
    cs = [r["snr"] for r in conv if r["snr"] is not None]
    fs = [r["snr"] for r in rows if not r["conv"] and r["snr"] is not None]
    if cs and fs:
        print(f"\n  median SNR converged {np.median(cs):.2f} dB, "
              f"failed {np.median(fs):.2f} dB")


if __name__ == "__main__":
    main()
