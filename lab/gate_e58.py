#!/usr/bin/env python3
"""E58 gates -- every lever must earn its numbers, with negative controls.

Legs (each with something that must fail):
  0  structural: the levers live in lab/e58_*.py ONLY -- no shared module of
     the default chain is modified (git diff of tracked files).
  1  engine: the batch float32 kernel ladder agrees with E49's recorded
     float64 scalar verdicts on the first-400 job set (E54's boundary allows
     converged-set deltas ONLY in the direction of extra BCH-clean decodes).
  2  oracle: in every E58 run, converged == BCH-zero, block for block.
     A converged-but-dirty block anywhere fails the campaign.
  3  negative control: blocks attributed wholly below 3.5 dB stay at
     (near-)zero convergence in every run -- no lever manufactures decodes
     from noise.
  4  permuted-sigma control: scrambling the weights within each frame must
     DESTROY the weighted lever's gain (<= plain + slack).  A control that
     also wins would mean the gain never came from the noise model.
  5  W=0 identity: the smoothed collector with the window off reproduces the
     untouched per-symbol path cell for cell (allclose 1e-12) on real air.
  6  exact-BP sanity: decodes a synthetic codeword through noise, and does
     NOT converge on shuffled cells (the syndrome stays an oracle).
"""
from __future__ import annotations

import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

PASS, FAIL = "PASS", "*** FAIL ***"


def leg0():
    r = subprocess.run(["git", "diff", "--name-only", "HEAD"],
                       cwd=os.path.dirname(HERE), capture_output=True,
                       text=True)
    touched = [l for l in r.stdout.splitlines() if l.strip()]
    # The decode chain E58's levers ride on: any diff here would mean a lever
    # leaked into the default path.  tools/* (the viewer stack) and other
    # sessions' lab files are outside E58's blast radius and are only listed.
    chain = [t for t in touched
             if t.startswith(("lab/m", "atsc3/", "spec/"))
             or t in ("selftest.py",)]
    other = [t for t in touched
             if t not in chain
             and not (os.path.basename(t).startswith(("e58_", "gate_e58"))
                      or t.endswith("LIVE_TV_LAB.md"))]
    if other:
        print(f"  [0] (info) non-chain diffs from concurrent work: {other}")
    print(f"  [0] decode-chain modules modified: {chain or 'none'}  "
          f"{PASS if not chain else FAIL}")
    return not chain


def leg1():
    e58 = json.load(open(os.path.join(HERE, "e58_gate_plain400.json")))
    e49 = json.load(open(os.path.join(HERE, "e49_deep_plp0_blocks.json")))
    ref = {e["block"]: e for e in e49["blockmap"]}
    agree = extra = lost = dirty = 0
    for m, r in e58["blocks"].items():
        m = int(m)
        if ref[m]["conv"] == r["conv"]:
            agree += 1
        elif r["conv"]:
            extra += 1
            if not r["bch"]:
                dirty += 1
        else:
            lost += 1
    ok = lost == 0 and dirty == 0 and agree >= 395
    print(f"  [1] engine vs E49 float64: agree {agree}, extra-BCH-clean "
          f"{extra}, lost {lost}, dirty {dirty}  {PASS if ok else FAIL}")
    return ok


def _runs():
    for name in ("e58_run_plain.json", "e58_run_weighted.json",
                 "e58_run_sm_plain.json", "e58_run_sm_weighted.json",
                 "e58_run_sm_weighted_sp.json"):
        p = os.path.join(HERE, name)
        if os.path.exists(p):
            yield name, json.load(open(p))


def leg2():
    ok = True
    for name, r in _runs():
        dirty = sum(1 for b in r["blocks"].values()
                    if b["conv"] and not b["bch"])
        ok &= dirty == 0
        print(f"  [2] {name}: converged {r['converged']}, "
              f"conv-but-BCH-dirty {dirty}  {PASS if dirty == 0 else FAIL}")
    return ok


def leg3():
    from e58_weighted import block_snrs
    ok = True
    for name, r in _runs():
        # comparisons bucket on the BASELINE instrument; but the sub-3.5
        # CONTROL must use each run's OWN instrument -- the smoothed CE
        # legitimately re-measures frames upward (median +1.7 dB), so "below
        # 3.5" means below 3.5 on the collect that produced the cells.
        diags = ("e58_sm_diags.json" if "_sm_" in name
                 else "e49_full_diags.json")
        snr = block_snrs(os.path.join(HERE, diags),
                         os.path.join(HERE, "m8_l1_rf25_amped.json"),
                         os.path.join(HERE, "e49_deep_plp0_blocks.json"))
        low = [m for m in r["blocks"] if snr.get(int(m), -99) < 3.5]
        c = sum(1 for m in low if r["blocks"][m]["conv"])
        good = c <= max(2, len(low) // 200)
        ok &= good
        print(f"  [3] {name}: sub-3.5dB control {c}/{len(low)} converged  "
              f"{PASS if good else FAIL}")
    return ok


def leg4():
    pp = os.path.join(HERE, "e58_gate_permute.json")
    if not os.path.exists(pp):
        print("  [4] permuted-sigma control: NOT RUN  " + FAIL)
        return False
    perm = json.load(open(pp))
    plain = json.load(open(os.path.join(HERE, "e58_run_plain.json")))
    wt = json.load(open(os.path.join(HERE, "e58_run_weighted.json")))
    ids = set(perm["blocks"])
    cp = sum(1 for m in ids if plain["blocks"][m]["conv"])
    cw = sum(1 for m in ids if wt["blocks"][m]["conv"])
    cx = perm["converged"]
    ok = cx <= cp + max(2, cp // 20) and cx < cw
    print(f"  [4] permute control on {len(ids)} jobs: plain {cp}, "
          f"weighted {cw}, permuted {cx}  {PASS if ok else FAIL}")
    return ok


def leg5():
    import m10_core as M10
    from m2_pilots import load
    from e58_collect import smoothed_pool
    js = json.load(open(os.path.join(HERE, "m8_l1_rf25_amped.json")))
    g = M10.geometry_from_json(js, "rf25")
    plps = {p["id"]: p for p in g.plps}
    per = M10.frame_samples(g) / M10.FS_POST
    cap = os.path.join(HERE, "..", "data", "fox25", "rf25_amped_0808.cs16")
    ok = True
    for n in (15, 154):
        s = max(0.0, (2 + n) * per - 0.01)
        span = (M10.frame_samples(g) + 40000) / M10.FS_POST
        y, _, _, _ = load(cap, 6.912e6, fmt=None, span_sec=span, start_sec=s)
        t0, _ = M10.fine_t0(y, g)
        pool0, _, _ = M10.frame_cells(cap, 6.912e6, None, g, s, plps[0])
        poolW0, _ = smoothed_pool(y, t0, g, W=0)
        same = np.allclose(pool0, poolW0, rtol=1e-12, atol=1e-12)
        ok &= same
        print(f"  [5] W=0 identity, frame {n}: {PASS if same else FAIL}")
    return ok


def leg6():
    import m6_bicm as B
    from e58_weighted import sum_product_decode_batch
    rng = np.random.default_rng(3)
    ch = B.PlpChain(64800, "QPSK", "9/15")
    pay = rng.integers(0, 2, ch.kpayload).astype(np.uint8)
    x = ch.encode(pay)
    n = (rng.standard_normal(len(x)) + 1j * rng.standard_normal(len(x))) \
        / np.sqrt(2)
    y = x + n * np.sqrt(10 ** (-3.0 / 10))
    lam = np.empty(ch.ninner)
    lam[ch.lam_of_q] = B.demap_llr(y, "QPSK", "9/15", sigma2=10 ** -0.3)
    b, c, _, _ = sum_product_decode_batch(lam[None, :], ch.checks,
                                          ch.ninner, iters=60)
    good1 = bool(c[0]) and ch.baseband_packet(b[0])["bch_ok"]
    ys = y[rng.permutation(len(y))]
    lam2 = np.empty(ch.ninner)
    lam2[ch.lam_of_q] = B.demap_llr(ys, "QPSK", "9/15", sigma2=10 ** -0.3)
    _, c2, _, bad2 = sum_product_decode_batch(lam2[None, :], ch.checks,
                                              ch.ninner, iters=30)
    good2 = not bool(c2[0])
    print(f"  [6] exact BP decodes AWGN codeword: {PASS if good1 else FAIL};"
          f" shuffled control did {'NOT converge' if good2 else 'CONVERGE (!!)'}"
          f"  {PASS if good2 else FAIL}")
    return good1 and good2


def main():
    print("E58 gates")
    print("=" * 72)
    ok = True
    for leg in (leg0, leg1, leg2, leg3, leg4, leg5, leg6):
        ok &= bool(leg())
    print(f"\n  {'ALL E58 GATES PASS' if ok else 'GATE FAILURES'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
