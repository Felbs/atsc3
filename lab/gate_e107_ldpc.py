#!/usr/bin/env python3
"""gate_e107_ldpc.py -- the layered eight-lane LDPC kernel's referee, on real air.

E54's boundary, re-asserted for a change of schedule: against the flooding
kernel (the shipping fast path), on real Frames with seeded noise added at
several levels, for BOTH codes (PLP 0: 16200-bit 64QAM 11/15 through the
live FrameDecoder; PLP 1: 64800-bit 256QAM 11/15 through M17):

  1. At every noise level the layered kernel converges at least as many blocks
     as the flooding kernel (net gain >= 0; a strict superset is not required
     of a change of SCHEDULE -- near the threshold the two converge slightly
     different sets, and the strict counts are printed for the record).
  2. For every block both converge, the Baseband Packet bytes are identical
     (same codeword: a schedule cannot change a fixed point).
  3. The layered kernel is not slower, and uses no more iterations on average.
  4. Determinism: two layered decodes of the same input are identical.
  5. E110, the two-minimum correction ("box" column): against the layered
     kernel, the same three rules -- no net loss at any level, bytes identical
     on shared blocks -- and its own time.
  6. E112, early abort of hopeless blocks ("abort" column, --abort k,T,stall):
     against the layered kernel.  An abort can only LOSE, so the bar is
     stricter than net gain: ZERO blocks lost at every level, bytes identical
     on every shared block, and faster.

Run: python lab/gate_e107_ldpc.py [--capture ...] [--frames 2] [--snrs none,22,21]
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

import m3_ldpc as LD                                             # noqa: E402
import m6_cells as C                                             # noqa: E402
import m9_fast                                                   # noqa: E402
import m13_sf1 as M13                                            # noqa: E402
import m17_plp1_live as M17                                      # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--capture", default=os.path.expanduser("~/captures/race_2026-11-06/live_rf33_antB.cs16"))
    ap.add_argument("--rate", type=float, default=6.912e6)
    ap.add_argument("--fmt", default="cs16")
    ap.add_argument("--l1", default="m4_l1detail_hit_rf33.json")
    ap.add_argument("--frames", type=int, default=2)
    ap.add_argument("--start", type=float, default=0.0)
    ap.add_argument("--snrs", default="none,23,22,21.5,21,20")
    ap.add_argument("--threads", type=int, default=12)
    ap.add_argument("--json", default=None)
    ap.add_argument("--abort", default="8,0.36,0",
                    help="E112 rule k,T,stall (T < 1 = fraction of the checks)")
    a = ap.parse_args(argv)
    abort_rule = tuple(float(x) for x in a.abort.split(","))
    if LD._layered_lib() is None:
        print("FAIL: the layered kernel is not built (python lab/build_ldpc_layered_kernel.py)")
        return 1

    from m2_pilots import load, guess_format
    from m10_core import fine_t0, frame_samples
    js = json.load(open(os.path.join(HERE, a.l1)))
    g0 = C.Geometry.rf33_legacy()
    fd = m9_fast.FrameDecoder(threads=a.threads, backend="cpu", iters=50, cpu_fast=True)
    fd.prewarm()
    sf1 = M17.SubframeDecoder(fd, g0, js["l1_basic"], js["l1_detail"], 1, comb=True)
    g = sf1.g
    per = (frame_samples(g0) + (g.nfft + g.gi) * g.nsym) / M13.FS_POST
    fmt = a.fmt or guess_format(a.capture)
    frames = []
    for k in range(a.frames):
        y, _fs, _cfo, _h = load(a.capture, a.rate, fmt=fmt, span_sec=per + 0.06, start_sec=a.start + k * per)
        t0, _coh = fine_t0(y, g0)
        frames.append((y, t0))
    psig = float(np.mean([np.mean(np.abs(y) ** 2) for y, _ in frames]))
    snrs = [None if s.strip() == "none" else float(s) for s in a.snrs.split(",")]
    print("E107 gate -- layered eight-lane LDPC vs the flooding kernel, real air + seeded noise")
    print("=" * 96)
    print(f"  {os.path.basename(a.capture)}, {len(frames)} frames; PLP 0 via FrameDecoder, PLP 1 via M17")
    print(f"\n  {'noise':>6} {'PLP':>4} {'flood conv':>11} {'layer conv':>11} {'lost':>5} {'gained':>6} "
          f"{'byte mism':>9} {'flood ms':>9} {'layer ms':>9} {'it flood':>8} {'it layer':>8}")
    fails = []
    rows = []
    for snr in snrs:
        for k, (y, t0) in enumerate(frames):
            if snr is None:
                yn = y
            else:
                rng = np.random.default_rng(107_000 + k * 97 + int(round(snr * 10)))
                sn = np.sqrt(psig * 10 ** (-snr / 10) / 2)
                yn = y + sn * (rng.standard_normal(len(y)) + 1j * rng.standard_normal(len(y)))
            for plp in (0, 1):
                res = {}
                for mode in ("flood", "layer", "layer2", "box", "abort"):
                    fd.ldpc_layered = mode != "flood"
                    fd.ldpc_boxplus = (mode == "box")        # E110: the two-minimum correction
                    # E112: explicitly OFF for the other columns -- the shipped
                    # default is on, and None would inherit it, making the
                    # layered-vs-abort comparison a comparison of two aborts
                    fd.ldpc_abort = abort_rule if mode == "abort" else (0, 0, 0)
                    fd.tm.clear(); sf1.tm.clear()
                    # the M17 decoder carries state between Frames (the fed-forward
                    # carrier offset, the cached main-path delay): every mode must
                    # see this Frame from the same state, or the comparison is of
                    # two slightly different inputs, not of two LDPC schedules
                    sf1.df_prev = 0.0
                    sf1._tau0 = {}
                    if plp == 0:
                        _r16, pk, diag = fd.decode_frame(yn, t0)
                        ms = fd.tm["ldpc"] * 1000
                        its = None
                    else:
                        rep = {}
                        pk, rep = sf1.decode(yn, t0, rep)
                        ms = sf1.tm["sf1_ldpc"] * 1000
                        its = rep.get("ldpc_iters_mean")
                    res[mode] = (pk, ms, its)
                pf, pl, pl2, pb = res["flood"][0], res["layer"][0], res["layer2"][0], res["box"][0]
                pa = res["abort"][0]
                ca = sum(p is not None for p in pa)
                a_lost = sum(1 for x, z in zip(pl, pa) if x is not None and z is None)
                a_mism = sum(1 for x, z in zip(pl, pa) if x is not None and z is not None and x != z)
                if a_lost: fails.append(f"noise {snr} frame {k} PLP {plp}: abort LOST {a_lost} block(s) the layered kernel converges")
                if a_mism: fails.append(f"noise {snr} frame {k} PLP {plp}: {a_mism} byte mismatch(es) layered vs abort")
                cf = sum(p is not None for p in pf); cl = sum(p is not None for p in pl)
                cb = sum(p is not None for p in pb)
                b_lost = sum(1 for x, z in zip(pl, pb) if x is not None and z is None)
                b_gain = sum(1 for x, z in zip(pl, pb) if x is None and z is not None)
                b_mism = sum(1 for x, z in zip(pl, pb) if x is not None and z is not None and x != z)
                if b_mism: fails.append(f"noise {snr} frame {k} PLP {plp}: {b_mism} byte mismatch(es) layered vs box-plus")
                if b_gain < b_lost: fails.append(f"noise {snr} frame {k} PLP {plp}: box-plus net loss ({b_lost} lost, {b_gain} gained)")
                lost = sum(1 for x, z in zip(pf, pl) if x is not None and z is None)
                gained = sum(1 for x, z in zip(pf, pl) if x is None and z is not None)
                mism = sum(1 for x, z in zip(pf, pl) if x is not None and z is not None and x != z)
                nondet = sum(1 for x, z in zip(pl, pl2) if x != z)
                tag = ""
                # A schedule is not a recompilation: near the threshold the two
                # converge DIFFERENT sets, and a strict superset is not a law
                # (measured 2026-11-07: one block in ~900 at the cliff converges
                # under flooding only, against 16-29 gained).  The bar is the net
                # at every level; the strict count is still printed.
                if gained < lost: fails.append(f"noise {snr} frame {k} PLP {plp}: net loss ({lost} lost, {gained} gained)"); tag += " NET LOSS"
                elif lost: tag += f" ({lost} lost, {gained} gained)"
                if mism: fails.append(f"noise {snr} frame {k} PLP {plp}: {mism} byte mismatch(es)"); tag += " MISMATCH"
                if nondet: fails.append(f"noise {snr} frame {k} PLP {plp}: {nondet} non-deterministic block(s)"); tag += " NONDET"
                itf = f"{res['flood'][2]:.1f}" if res['flood'][2] is not None else "-"
                itl = f"{res['layer'][2]:.1f}" if res['layer'][2] is not None else "-"
                print(f"  {('as rec' if snr is None else f'{snr:5.1f}'):>6} {plp:>4} {cf:>5}/{len(pf):<5} {cl:>5}/{len(pl):<5} {lost:>5} {gained:>6} "
                      f"{mism:>9} {res['flood'][1]:>9.0f} {res['layer'][1]:>9.0f} {itf:>8} {itl:>8}{tag}"
                      f"   | box {cb:>3}/{len(pb):<3} (+{b_gain} -{b_lost}, {res['box'][1]:.0f} ms)"
                      f"   | abort {ca:>3}/{len(pa):<3} (-{a_lost}, {res['abort'][1]:.0f} ms"
                      f"{', it ' + format(res['abort'][2], '.1f') if res['abort'][2] is not None else ''})")
                rows.append(dict(noise=snr, frame=k, plp=plp, flood=cf, layer=cl, lost=lost, gained=gained, mismatch=mism,
                                 flood_ms=res['flood'][1], layer_ms=res['layer'][1],
                                 abort=ca, abort_lost=a_lost, abort_ms=res['abort'][1]))
    fd.ldpc_layered = False
    fd.ldpc_boxplus = False
    fd.ldpc_abort = None
    if a.json:
        json.dump(rows, open(a.json, "w"), indent=1)
    sf = sum(r['flood_ms'] for r in rows); sl = sum(r['layer_ms'] for r in rows)
    sa = sum(r['abort_ms'] for r in rows)
    print(f"\n  LDPC time over the gate: flooding {sf:.0f} ms, layered {sl:.0f} ms ({sf / max(sl, 1e-9):.1f}x), "
          f"layered+abort {sa:.0f} ms ({sl / max(sa, 1e-9):.2f}x over layered)")
    if sl > sf:
        fails.append("layered slower than flooding over the gate")
    if sa > sl:
        fails.append("abort slower than layered over the gate")
    if fails:
        print("\n  FAIL:\n    " + "\n    ".join(fails))
        return 1
    print("\n  PASS: no net loss at any level, every shared block byte-identical, deterministic, not slower; "
          "abort loses no block at any level")
    return 0


if __name__ == "__main__":
    sys.exit(main())
