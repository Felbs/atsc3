#!/usr/bin/env python3
"""E58 -- CTI-deinterleave a collected PLP0 cell stream (e49_stream's math,
as a standalone step so the smoothed collect plugs into e58_weighted decode).
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

import m10_core as M10               # noqa: E402
import m10_cti as CTI                # noqa: E402
from e49_stream import deinterleave_chunked  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("prefix")
    ap.add_argument("--l1", default="m8_l1_rf25_amped.json")
    ap.add_argument("--plp", type=int, default=0)
    a = ap.parse_args()
    js = json.load(open(a.l1))
    g = M10.geometry_from_json(js, a.prefix)
    plp = {p["id"]: p for p in g.plps}[a.plp]
    nrows = CTI.nrows_of(plp["cti_depth"], bool(plp["cti_extended"]))
    src = np.load(f"{a.prefix}_plp{a.plp}.npy", mmap_mode="r")
    recv = np.ascontiguousarray(src.reshape(-1))
    cells, valid = deinterleave_chunked(recv, nrows, plp["cti_start_row"])
    del recv
    out = f"{a.prefix}_plp{a.plp}_deint.npy"
    np.save(out, cells)
    first_bad = np.flatnonzero(~valid)
    nvalid = int(first_bad[0]) if len(first_bad) else len(valid)
    print(f"  wrote {out}  ({len(cells)} cells, {nvalid} contiguously valid)")


if __name__ == "__main__":
    main()
