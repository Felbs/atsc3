#!/usr/bin/env python3
"""E49 -- m10_core driver with a frame offset.

WHY THIS EXISTS
---------------
On a fady capture m8_l1's t=0 solve can fail, so its json records the L1 of a
LATER frame (`preamble.start_sec` says which attempt won).  m10_core starts
its frame loop at t=0 regardless, so the json's `L1D_plp_CTI_start_row` then
belongs to a different Subframe than the first cells fed to the CTI -- and a
desynced CTI turns EVERY PLP into ~chance unsatisfied-check counts, which
impersonates a weak link (the M8 lesson, in commutator costume).

The commutator never resets (A/322 9.3.9.1): start_row advances by plp_size
per Subframe.  So the honest fix is to start the frame loop AT the json's own
frame.  `--start-frame k` shifts every load offset by k frames.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m10_core as M10               # noqa: E402
import m7_route as R7                # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("capture")
    ap.add_argument("--rate", type=float, required=True)
    ap.add_argument("--l1", required=True)
    ap.add_argument("--start-frame", type=int, default=0)
    ap.add_argument("--frames", type=int, default=6)
    ap.add_argument("--plp", type=int, default=0)
    ap.add_argument("--iters", type=int, default=50)
    ap.add_argument("--max-blocks", type=int, default=None)
    ap.add_argument("--bb", default=None)
    ap.add_argument("--json")
    a = ap.parse_args()

    js = json.load(open(a.l1))
    g = M10.geometry_from_json(js, os.path.basename(a.capture))
    per = M10.frame_samples(g) / M10.FS_POST
    off = a.start_frame * per

    orig = M10.frame_cells

    def shifted(path, rate, fmt, gg, start_sec, *args, **kw):
        return orig(path, rate, fmt, gg, start_sec + off, *args, **kw)

    M10.frame_cells = shifted
    print(f"E49 -- m10 core chain, frame loop starting at capture frame "
          f"{a.start_frame} (t = {off:.4f} s)")
    r = M10.run(a.capture, a.rate, None, js, a.frames, plp_id=a.plp,
                iters=a.iters, max_blocks=a.max_blocks)
    M10.frame_cells = orig
    if r is None:
        print("  no Frames decoded.")
        return 1

    if a.bb and r["stream"]:
        diags = [dict(d) for d in r["diags"]]
        for k, d in enumerate(diags):
            d["conv"] = r["converged"] if k == 0 else 0
            d["bch"] = r["bch_ok"] if k == 0 else 0
            d["plp16"] = False
            d["coh"] = d.get("data_coh", 0.0)
        with open(a.bb, "wb") as fh:
            R7.bb_write(fh, 0.0, diags, r["stream"], r["bounds"])
        print(f"  wrote {a.bb}  ({len(r['stream'])/1e6:.2f} MB, "
              f"{len(r['bounds'])} ALP anchors)")

    if a.json:
        json.dump(dict(capture=os.path.basename(a.capture), plp=a.plp,
                       start_frame=a.start_frame, frames=r["n_frames"],
                       blocks=r["n_blocks"], converged=r["converged"],
                       bch_ok=r["bch_ok"],
                       stream_bytes=len(r["stream"])),
                  open(a.json, "w"), indent=2, default=float)
        print(f"  wrote {a.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
