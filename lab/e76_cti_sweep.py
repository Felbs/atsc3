#!/usr/bin/env python3
"""E76 -- is g4's 0/400 a dead link, or the interleaver clock?

The two instruments disagree violently on rf25_fox_g4_0810_2254.cs16:

    dummy-cell SNR   median 12.85 dB, 150/150 frames over threshold
    pilot coherence  0.985-0.994 (the E49 control reads 0.887)
    constellation    tighter than the control (p99 2.08 vs 2.81)
    -------------------------------------------------------------
    L1-Basic CRC     fails on ~8 of 9 probes
    PLP 0 decode     0/400 LDPC converged

Every quality-shaped metric says excellent; every correctness-shaped metric
says nothing decodes.  That is not what a weak link looks like -- a weak link
degrades the quality metrics first.  It is what E49 already recorded as a
LAW: "interleaver clock offset impersonates weak SNR."  A CTI commutator one
frame out of phase leaves the constellation pristine and the FEC at zero.

So sweep the frame phase.  e49_stream's --start-frame aligns the collect to
the L1's own Frame; if that alignment is off by k, dropping k frames from the
already-collected array must recover it.  The right k lights up; the wrong k
stays at zero.  The E49 capture at its known-good alignment is the positive
control in the same run -- if it does not converge here, the sweep proves
nothing.
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

import m10_core as M10   # noqa: E402
import m10_cti as CTI    # noqa: E402
import m6_bicm as B      # noqa: E402


def deint(recv, nrows, start_row, chunk=8_000_000):
    L = len(recv)
    out = np.zeros(L, recv.dtype)
    valid = np.zeros(L, bool)
    for lo in range(0, L, chunk):
        hi = min(lo + chunk, L)
        i = np.arange(lo, hi, dtype=np.int64)
        q = i + nrows * ((start_row + i) % nrows)
        ok = q < L
        out[lo:hi][ok] = recv[q[ok]]
        valid[lo:hi] = ok
    return out, valid


def trial(prefix, l1, kdrop, nframes, nblocks, iters):
    js = json.load(open(os.path.join(HERE, l1)))
    g = M10.geometry_from_json(js, prefix)
    plp = {p["id"]: p for p in g.plps}[0]
    nrows = CTI.nrows_of(plp["cti_depth"], bool(plp["cti_extended"]))
    ncell = plp["cells_per_fec"]
    # advancing the collect by k Frames advances the commutator by
    # k * plp_size rows (A/322 9.3.9.1: start_row += plp_size per Subframe)
    sr = (plp["cti_start_row"] + kdrop * plp["size"]) % nrows
    C0 = CTI.solve_C(plp["cti_fec_block_start"], sr, nrows)
    arr = np.load(os.path.join(HERE, f"{prefix}_plp0.npy"), mmap_mode="r")
    seg = np.ascontiguousarray(
        np.asarray(arr[kdrop:kdrop + nframes]).reshape(-1))
    cells, valid = deint(seg, nrows, sr)
    del seg
    nv = int(np.flatnonzero(~valid)[0]) if (~valid).any() else len(valid)
    nb = max(0, (nv - C0) // ncell)
    ch = B.PlpChain(plp["ninner"], plp["mod"], plp["rate"], iters=iters)
    conv = bch = 0
    tried = min(nblocks, nb)
    best = 1 << 30
    for m in range(tried):
        lo = C0 + m * ncell
        r = ch.decode(np.asarray(cells[lo:lo + ncell], np.complex128))
        best = min(best, int(r["unsatisfied"]))
        if r["converged"]:
            conv += 1
            if ch.baseband_packet(r["bits"])["bch_ok"]:
                bch += 1
    return dict(k=kdrop, start_row=sr, C=C0, blocks_avail=nb, tried=tried,
                converged=conv, bch=bch, min_unsatisfied=best)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", default="e76_g4")
    ap.add_argument("--l1", default="e76_l1_g4.json")
    ap.add_argument("--kmax", type=int, default=8)
    ap.add_argument("--frames", type=int, default=16)
    ap.add_argument("--blocks", type=int, default=12)
    ap.add_argument("--iters", type=int, default=60)
    ap.add_argument("--control", action="store_true",
                    help="also run the E49 capture at its known-good phase")
    a = ap.parse_args()

    print("%-10s %-4s %-6s %-7s %-8s %-8s %-8s %s"
          % ("prefix", "k", "srow", "C", "avail", "tried", "conv/BCH",
             "min unsat"))
    print("-" * 74)
    if a.control:
        r = trial("e76_ctl", "m8_l1_rf25_amped.json", 0, a.frames, a.blocks,
                  a.iters)
        print("%-10s %-4d %-6d %-7d %-8d %-8d %-8s %d   <== POSITIVE CONTROL"
              % ("e76_ctl", r["k"], r["start_row"], r["C"], r["blocks_avail"],
                 r["tried"], "%d/%d" % (r["converged"], r["bch"]),
                 r["min_unsatisfied"]))
    out = []
    for k in range(a.kmax + 1):
        r = trial(a.prefix, a.l1, k, a.frames, a.blocks, a.iters)
        out.append(r)
        print("%-10s %-4d %-6d %-7d %-8d %-8d %-8s %d"
              % (a.prefix, r["k"], r["start_row"], r["C"], r["blocks_avail"],
                 r["tried"], "%d/%d" % (r["converged"], r["bch"]),
                 r["min_unsatisfied"]))
    json.dump(out, open(os.path.join(HERE, f"{a.prefix}_cti_sweep.json"), "w"),
              indent=1)


if __name__ == "__main__":
    main()
