#!/usr/bin/env python3
"""M2 Stage A -- decode EVERY bootstrap in a real capture, and arbitrate A1.

No radio.  Offline only.  Reads a banked .cs16, resamples to the fixed
6.144 Msps bootstrap rate in overlapping chunks, finds every bootstrap
occurrence in the record with the 2048-sample coherent acquisition, decodes
the A/321 signaling bits at each one, and reports:

  * per-occurrence bits + decoded fields
  * FIELD STABILITY across all occurrences (the honesty check: a real
    transmitter repeats the same bootstrap every frame; instability means we
    are misreading, not that the station is changing)
  * the observed frame period, cross-checked against min_time_to_next
  * ASSUMPTION A1 arbitration: the same decode under all four PN_VARIANTS.
    Real air, not our own synthesis, picks the winner.

Usage:
    python m2_stage_a.py hit_rf33.cs16 --rate 8e6
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from atsc3 import bootstrap as bs                       # noqa: E402
from atsc3.capture import (guess_format, n_samples, read_block,      # noqa: E402
                           to_bootstrap_rate)

CHUNK_SEC = 1.0
# overlap must hold a whole bootstrap (4*3072) plus the acquisition window
GUARD_BS = 4 * bs.N_SYM + bs.N_FFT + 64


def fine_cfo(y, p, n_symbols=4):
    """Fractional CFO from the lag-2048 part-C/part-A repeat, folded over all
    four bootstrap symbols.  ASSUMPTION-FREE: uses only the A/321 time geometry
    (part C is the last 520 samples of part A), no ZC root and no PN.
    Unambiguous only within +-FS/(2*2048) = +-1500 Hz."""
    acc = 0.0 + 0.0j
    for n in range(n_symbols):
        o = p + n * bs.N_SYM + bs.run_offset(n)
        if o + bs.N_A + bs.N_C > len(y):
            continue
        a = y[o: o + bs.N_C]
        b = y[o + bs.N_A: o + bs.N_A + bs.N_C]
        acc += np.vdot(b, a)          # sum a*conj(b)
    if acc == 0:
        return float("nan")
    return float(-np.angle(acc) * bs.FS / (2.0 * np.pi * bs.N_A))


def find_bootstraps(y, cfo_grid=(0.0,), pn_variant="spec", minor_version=0,
                    q=bs.ZC_ROOT_MAJOR0, top_k=40, min_ratio=8.0,
                    int_search=3):
    """Every bootstrap in a 6.144 Msps block.  Returns list of detect_matched
    dicts with 'position' added.  Acquisition is the symbol-0 template, which
    A/321 5.3.3 fixes at absolute cyclic shift 0 -- a fully known waveform."""
    cands = bs.acquire_matched(y, cfo_grid=cfo_grid, minor_version=minor_version,
                               q=q, pn_variant=pn_variant, top_k=top_k)
    out = []
    for c in cands:
        m = bs.detect_matched(y, c["position"], c.get("cfo_hz") or 0.0,
                              minor_version, q, pn_variant, 4, int_search)
        if not m or "mean_peak_ratio" not in m:
            continue
        if m["mean_peak_ratio"] < min_ratio:
            continue
        m["position"] = c["position"]
        m["acq_ratio"] = c["peak_ratio"]
        m["fine_cfo_hz"] = fine_cfo(y, c["position"])
        out.append(m)
    out.sort(key=lambda d: d["position"])
    return out


def scan(path, fs_in, fmt=None, pn_variant="spec", max_sec=None,
         chunk_sec=CHUNK_SEC, verbose=True):
    fmt = fmt or guess_format(path)
    total = n_samples(path, fmt)
    if max_sec:
        total = min(total, int(fs_in * max_sec))
    step = int(chunk_sec * fs_in)
    guard_in = int(GUARD_BS * fs_in / bs.FS) + 64
    hits = []
    pos = 0
    while pos < total:
        cnt = min(step + guard_in, total - pos)
        if cnt < int(0.05 * fs_in):
            break
        x = read_block(path, pos, cnt, fmt)
        if len(x) < 1000:
            break
        y = to_bootstrap_rate(x, fs_in)
        base = pos * bs.FS / fs_in          # absolute index in the 6.144 domain
        found = find_bootstraps(y, pn_variant=pn_variant)
        for h in found:
            h["abs_pos"] = base + h["position"]
            h["t_sec"] = h["abs_pos"] / bs.FS
            # drop hits that live entirely in the guard of the NEXT chunk
            if h["position"] >= step * bs.FS / fs_in and pos + step < total:
                continue
            hits.append(h)
        if verbose:
            print(f"  chunk t={pos/fs_in:5.2f}s  bootstraps={len(found)}",
                  flush=True)
        pos += step
    hits.sort(key=lambda d: d["abs_pos"])
    return hits


def summarize(hits):
    """Field stability across occurrences."""
    keys = ["ea_wake_up_1", "min_time_to_next", "system_bandwidth",
            "ea_wake_up_2", "bsr_coefficient", "preamble_structure"]
    tally = {k: Counter() for k in keys}
    bitstr = Counter()
    for h in hits:
        f = h.get("fields", {})
        for k in keys:
            if k in f:
                tally[k][f[k]] += 1
        b = h.get("bits", {})
        s = " ".join("".join(str(x) for x in b[n]) for n in sorted(b))
        bitstr[s] += 1
    return tally, bitstr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("path")
    ap.add_argument("--rate", type=float, default=8e6)
    ap.add_argument("--max-sec", type=float, default=None)
    ap.add_argument("--pn-sweep", action="store_true")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    path = a.path if os.path.isabs(a.path) else os.path.join(HERE, a.path)
    print(f"=== M2 STAGE A === {os.path.basename(path)} @ {a.rate/1e6:g} Msps")
    t0 = time.time()
    hits = scan(path, a.rate, max_sec=a.max_sec)
    print(f"\n{len(hits)} bootstrap occurrences in "
          f"{(hits[-1]['t_sec'] - hits[0]['t_sec']) if len(hits) > 1 else 0:.3f} s "
          f"({time.time()-t0:.1f}s wall)")

    print("\n--- per-occurrence ---")
    print(f"{'#':>3} {'t_sec':>8} {'mratio':>8} {'tres':>5} {'cfo_Hz':>8} "
          f"{'shifts':>20}  sym1     sym2     sym3")
    for i, h in enumerate(hits):
        b = h["bits"]
        print(f"{i:3d} {h['t_sec']:8.4f} {h['mean_peak_ratio']:8.1f} "
              f"{h['timing_residual']:5d} {h['fine_cfo_hz']:8.1f} "
              f"{str(h['abs_shifts']):>20}  "
              + "  ".join("".join(str(x) for x in b[n]) for n in sorted(b)))

    if len(hits) > 1:
        d = np.diff([h["abs_pos"] for h in hits]) / bs.FS * 1000.0
        print(f"\nframe period: mean {d.mean():.3f} ms  sd {d.std():.4f} ms  "
              f"min {d.min():.3f}  max {d.max():.3f}  (n={len(d)})")

    tally, bitstr = summarize(hits)
    print("\n--- field stability (value: count) ---")
    for k, c in tally.items():
        flag = "STABLE" if len(c) == 1 else "*** UNSTABLE ***"
        print(f"  {k:22s} {dict(c)}   {flag}")
    print("\n--- raw signaling bitstrings (sym1 sym2 sym3) ---")
    for s, n in bitstr.most_common():
        print(f"  {s}   x{n}")

    if hits:
        f = hits[0]["fields"]
        print("\n--- decoded meaning (modal occurrence) ---")
        for k, v in f.items():
            print(f"  {k:32s} {v}")

    pn_report = {}
    if a.pn_sweep:
        print("\n--- ASSUMPTION A1 arbitration: PN wiring vs REAL AIR ---")
        x = read_block(path, 0, int(1.2 * a.rate), guess_format(path))
        y = to_bootstrap_rate(x, a.rate)
        for v in bs.PN_VARIANTS:
            hv = find_bootstraps(y, pn_variant=v, min_ratio=0.0, top_k=8)
            best = max((h["mean_peak_ratio"] for h in hv), default=0.0)
            n_ok = sum(1 for h in hv if h["mean_peak_ratio"] >= 8.0)
            pn_report[v] = {"best_mean_peak_ratio": best, "n_detected": n_ok,
                            "desc": bs.PN_VARIANTS[v][2]}
            print(f"  {v:15s} best_mean_peak_ratio {best:8.2f}   "
                  f"detected {n_ok}   {bs.PN_VARIANTS[v][2]}")

    out = a.out or os.path.join(HERE, "m2_stage_a_" +
                                os.path.splitext(os.path.basename(path))[0] + ".json")
    with open(out, "w", encoding="utf-8") as fh:
        json.dump({
            "path": path, "rate": a.rate, "n_bootstraps": len(hits),
            "occurrences": [{k: (v if not isinstance(v, dict) else
                                 {str(kk): vv for kk, vv in v.items()})
                             for k, v in h.items()
                             if k not in ("phases",)} for h in hits],
            "field_tally": {k: {str(kk): vv for kk, vv in c.items()}
                            for k, c in tally.items()},
            "pn_arbitration": pn_report,
        }, fh, indent=1, default=str)
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
