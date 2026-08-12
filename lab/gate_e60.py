#!/usr/bin/env python3
"""E60 gates -- the E58 levers in the LIVE chain, each with a way to fail.

Legs:
  0  structural: the e58_*/e49_* instruments are UNTOUCHED (they are the
     references m16_margin is held against).
  1  port equality (RF25, real air): m16_margin.smoothed_pool_g == the
     e58_collect instrument cell for cell at W=0 AND W=12; sigma_parts ==
     the banked e58 sigma rows exactly; sum_product_decode_batch(deadline=
     None) == the e58 decoder verdict for verdict on the same LLRs.
  2  live identity, levers OFF: cell_pool_fast_sm in all-fallback mode is
     BIT-identical to cell_pool_fast on real RF33 air, and the --no-margin
     e2e datagram sha equals the banked HEAD sha.
  3  strong-signal identity, levers ON (the port's default): the margin-on
     e2e run decodes the SAME FEC counts, BCH counts and datagram sha as
     HEAD, with the weighted-LLR auto-switch counter RECORDED and (near-)
     zero.  Every full frame of the gate capture measures 17.1-20.0 dB
     dummy-SNR (offline sweep) -- above the 16.0 gate; the chain's 32nd,
     EOF-truncated frame window genuinely reads low and MAY engage (that is
     the switch doing its job on a low-quality frame).  <= 3 engagements
     with an EQUAL sha passes; more, or any byte difference, fails.
  4  fade tracking (negative control included): synthetic step fades on the
     RF33 gate capture, three classes -- amplitude x0.35, phase 90deg, and
     a delay (SHAPE) step, each over 11 of 35 symbols.  The DEFAULT build
     (complex-gain normalisation + relative/absolute change detector) must
     decode within 1 block of the untouched per-symbol path on ALL of them;
     the OLD always-smooth build (no normalisation, no detector -- what E58
     would have shipped raw) must LOSE blocks on the amplitude fade.  If
     that control passes, the gate is vacuous and FAILS.  A clean frame
     must produce (near-)zero fallbacks -- the detector must not fire on
     healthy air.
  5  SP rescue budget (E52 shed law): the exact-BP batch honours its
     deadline on never-converging input (wall <= budget + one iteration)
     and counts what it sheds; with no deadline it still decodes a real
     16200 11/15 codeword through noise and refuses shuffled cells.
  6  the Fox referee: the full RF25 stack re-run THROUGH the ported module
     (e60_fox.py) reproduces E58's 1056 blocks, all BCH-zero, on E49's
     8697-job list.  Skipped with a warning while the background run is
     still in flight.
  7  switch engagement: AWGN added to a strong frame until the dummy-SNR
     sits below the gate -- the auto-switch must ENGAGE (n_wllr_frames = 1)
     and the weighted decode must not lose blocks against the plain decode
     of the same air (the lever is for the cliff; at minimum it must do no
     harm there).
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

PASS, FAIL = "PASS", "*** FAIL ***"
E60D = os.path.join(HERE, "..", "data", "e60")


def leg0():
    r = subprocess.run(["git", "diff", "--name-only", "HEAD"],
                       cwd=os.path.dirname(HERE), capture_output=True,
                       text=True)
    touched = [l for l in r.stdout.splitlines()
               if os.path.basename(l).startswith(("e58_", "e49_",
                                                  "gate_e58"))]
    print(f"  [0] e58/e49 instruments modified: {touched or 'none'}  "
          f"{PASS if not touched else FAIL}")
    return not touched


def leg1():
    import m10_core as M10
    import m16_margin as M16
    from m2_pilots import load
    from e58_collect import smoothed_pool as inst_pool
    from e58_weighted import (sum_product_decode_batch as inst_sp,
                              qpsk_d2min as inst_d2)
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
        for W in (0, 12):
            pa, _ = inst_pool(y, t0, g, W=W)
            pb, _ = M16.smoothed_pool_g(y, t0, g, W=W)
            same = np.array_equal(pa, pb)
            ok &= same
            print(f"  [1] smoothed_pool_g == instrument, frame {n} W={W}: "
                  f"{PASS if same else FAIL}")
    cum = M16.pool_symbol_bounds(g)
    p0 = plps[0]
    src = np.load(os.path.join(HERE, "e58_sm_plp0.npy"), mmap_mode="r")
    tab = np.load(os.path.join(HERE, "e58_sm_sigma.npy"))
    good = True
    for f in (15, 154, 300, 480):
        row = M16.sigma_parts(np.asarray(src[f]), cum, p0["start"],
                              p0["size"])
        ref = tab[f]
        good &= (np.array_equal(np.isnan(row), np.isnan(ref))
                 and np.array_equal(row[~np.isnan(row)],
                                    ref[~np.isnan(ref)]))
    ok &= good
    print(f"  [1] sigma_parts == banked e58 sigma rows (4 frames): "
          f"{PASS if good else FAIL}")
    # d2 helper + exact-BP equivalence on identical inputs
    rng = np.random.default_rng(60)
    z = rng.standard_normal(4096) + 1j * rng.standard_normal(4096)
    good = np.array_equal(M16.qpsk_d2min(z), inst_d2(z))
    ok &= good
    print(f"  [1] qpsk_d2min == instrument: {PASS if good else FAIL}")
    import m6_bicm as B
    ch = B.PlpChain(16200, "64QAM", "11/15")
    lam = rng.standard_normal((2, ch.ninner)) * 2.0
    ra = M16.sum_product_decode_batch(lam, ch.checks, ch.ninner, iters=25)
    rb = inst_sp(lam, ch.checks, ch.ninner, iters=25)
    good = all(np.array_equal(x, y) for x, y in zip(ra, rb))
    ok &= good
    print(f"  [1] sum-product(deadline=None) == instrument (bits/conv/its/"
          f"bad): {PASS if good else FAIL}")
    return ok


def _rf33_setup():
    import m9_fast
    import m16_margin as M16
    y, fs, cfo, h = m9_fast.load_fast(
        os.path.join(HERE, "rate6912_rf33.cs16"), 6912000, span_sec=1.2)
    fd_off = m9_fast.FrameDecoder(threads=6, backend="cpu", cpu_fast=True,
                                  margin=M16.OFF)
    fd_off.prewarm()
    t0, _ = fd_off.fine_timing(y)
    return m9_fast, M16, y, t0, fd_off


def leg2(state):
    m9_fast, M16, y, t0, fd_off = state
    fd_off._cpf = getattr(fd_off, "_cpf", None) or fd_off._cell_pool_tables()
    pa, ia = fd_off.cell_pool_fast(y, t0)
    pb, ib = fd_off.cell_pool_fast_sm(y, t0)     # margin OFF -> all-fallback
    same = (np.array_equal(pa, pb)
            and np.array_equal(ia["symbol_of"], ib["symbol_of"]))
    print(f"  [2] cell_pool_fast_sm(all-fallback) BIT-identical to "
          f"cell_pool_fast on real air: {PASS if same else FAIL}")
    ok = same
    try:
        head = json.load(open(os.path.join(E60D, "head_1.json")))
        nom = json.load(open(os.path.join(E60D, "nomargin_1.json")))
        good = (head["dg"]["sha256"] == nom["dg"]["sha256"]
                and head["fec_converged"] == nom["fec_converged"])
        print(f"  [2] --no-margin e2e sha == HEAD sha "
              f"({head['dg']['sha256'][:16]}...): {PASS if good else FAIL}")
        ok &= good
    except FileNotFoundError:
        print(f"  [2] --no-margin e2e artifacts missing  {FAIL}")
        ok = False
    return ok


def leg3():
    ok = True
    try:
        head = json.load(open(os.path.join(E60D, "head_1.json")))
    except FileNotFoundError:
        print(f"  [3] HEAD baseline missing  {FAIL}")
        return False
    for i in (1, 2, 3):
        p = os.path.join(E60D, f"margin_{i}.json")
        if not os.path.exists(p):
            print(f"  [3] margin_{i}.json missing  {FAIL}")
            ok = False
            continue
        m = json.load(open(p))
        same = (m["dg"]["sha256"] == head["dg"]["sha256"]
                and m["fec_converged"] == head["fec_converged"]
                and m["bch_zero"] == head["bch_zero"])
        # the counter must be RECORDED (missing == -1 fails) and near zero:
        # only the EOF-truncated tail frames may legitimately engage
        nw = m.get("margin_counters", {}).get("n_wllr_frames", -1)
        good = same and 0 <= nw <= 3
        ok &= good
        print(f"  [3] margin-on run {i}: FEC {m['fec_converged']}=="
              f"{head['fec_converged']}, sha "
              f"{'EQUAL' if same else 'DIFFER'}, wllr engaged on {nw} "
              f"frames (EOF-tail)  {PASS if good else FAIL}")
    return ok


def leg4(state):
    import m6_cells as C
    m9_fast, M16, y, t0, fd_off = state
    step = C.NFFT + C.GI
    fade_syms = list(range(12, 23))
    fd_def = m9_fast.FrameDecoder(
        threads=6, backend="cpu", cpu_fast=True,
        margin=dict(ce_w=12, wllr="off", sp=0))
    fd_old = m9_fast.FrameDecoder(
        threads=6, backend="cpu", cpu_fast=True,
        margin=dict(ce_w=12, ce_detect=1e12, ce_norm=0, wllr="off", sp=0))
    fd_old.margin.ce_abs = 1e12
    for fd in (fd_def, fd_old):
        fd.prewarm()

    def inj(alpha=1.0, theta=0.0, shift=0):
        yf = y.copy()
        for l in fade_syms:
            s = t0 + step * (1 + l)
            if shift:
                yf[s:s + step] = y[s - shift:s + step - shift]
            else:
                yf[s:s + step] = y[s:s + step] * (alpha
                                                  * np.exp(1j * theta))
        return yf

    ok = True
    # negative control: the raw always-smooth build must fail
    yf = inj(alpha=0.35)
    _, _, dg_raw = fd_off.decode_frame(yf, t0)
    _, _, dg_old = fd_old.decode_frame(yf, t0)
    ctl = dg_old["converged"] <= dg_raw["converged"] - 3
    ok &= ctl
    print(f"  [4] amp x0.35 fade, OLD always-smooth (no norm/detector): "
          f"{dg_old['converged']}/74 vs raw {dg_raw['converged']}/74 "
          f"(control {'FAILED as required' if ctl else 'PASSED -- vacuous'})"
          f"  {PASS if ctl else FAIL}")
    cases = [("amp x0.35", dict(alpha=0.35)),
             ("phase 90deg", dict(theta=np.pi / 2)),
             ("delay 16 samples", dict(shift=16))]
    for name, kw in cases:
        yf = inj(**kw)
        _, _, dg_raw = fd_off.decode_frame(yf, t0)
        rep = {}
        _, _, dg_def = fd_def.decode_frame(yf, t0, rep)
        good = dg_def["converged"] >= dg_raw["converged"] - 1
        ok &= good
        print(f"  [4] {name}: raw {dg_raw['converged']}/74, DEFAULT "
              f"{dg_def['converged']}/74, fell_back {rep.get('fell_back')}"
              f"/35  {PASS if good else FAIL}")
    rep_clean = {}
    fd_def.decode_frame(y, t0, rep_clean)
    quiet = rep_clean.get("fell_back", 99) <= 2
    ok &= quiet
    print(f"  [4] clean frame: fell_back {rep_clean.get('fell_back')}/35  "
          f"{PASS if quiet else FAIL}")
    fd_def.close()
    fd_old.close()
    return ok


def leg5():
    import m6_bicm as B
    import m16_margin as M16
    ch = B.PlpChain(16200, "64QAM", "11/15")
    rng = np.random.default_rng(6)
    # deadline honoured on never-converging input
    lam = rng.standard_normal((8, ch.ninner)) * 3.0
    M16.sum_product_decode_batch(lam[:1], ch.checks, ch.ninner, iters=1)
    budget = 0.05
    st = {}
    t = time.perf_counter()
    _, conv, _, _ = M16.sum_product_decode_batch(
        lam, ch.checks, ch.ninner, iters=100,
        deadline=time.perf_counter() + budget, stats=st)
    wall = time.perf_counter() - t
    timed = (not conv.any()) and st.get("sp_timeout", 0) == 8
    bounded = wall <= budget + 0.30          # one iteration of slack
    print(f"  [5] deadline: wall {wall*1000:.0f} ms against a {budget*1000:.0f}"
          f" ms budget, shed {st.get('sp_timeout', 0)}/8 counted  "
          f"{PASS if timed and bounded else FAIL}")
    # correctness with no deadline: real codeword through noise + control
    pay = rng.integers(0, 2, ch.kpayload).astype(np.uint8)
    x = ch.encode(pay)
    n = (rng.standard_normal(len(x)) + 1j * rng.standard_normal(len(x)))
    # 16 dB: ~1 dB above this code's measured AWGN threshold (min-sum and
    # exact BP both fail at 14 dB -- the cliff is ~15, not the 11-12 the
    # first draft of this gate assumed)
    y = x + n * np.sqrt(10 ** (-16.0 / 10) / 2)
    lam1 = ch.llr(y)
    b, c, _, _ = M16.sum_product_decode_batch(lam1[None, :], ch.checks,
                                              ch.ninner, iters=60)
    good1 = bool(c[0]) and ch.baseband_packet(b[0])["bch_ok"]
    ys = y[rng.permutation(len(y))]
    lam2 = ch.llr(ys)
    _, c2, _, _ = M16.sum_product_decode_batch(lam2[None, :], ch.checks,
                                               ch.ninner, iters=30)
    good2 = not bool(c2[0])
    print(f"  [5] exact BP decodes a 16200 11/15 codeword: "
          f"{PASS if good1 else FAIL}; shuffled control "
          f"{'did NOT converge' if good2 else 'CONVERGED (!!)'}  "
          f"{PASS if good2 else FAIL}")
    return timed and bounded and good1 and good2


def leg6():
    p = os.path.join(HERE, "e60_run_sm_weighted_sp.json")
    if not os.path.exists(p):
        print(f"  [6] Fox reproduction not yet on disk (background run) -- "
              f"{FAIL}")
        return False
    r = json.load(open(p))
    ref = json.load(open(os.path.join(HERE, "e58_run_sm_weighted_sp.json")))
    dirty = sum(1 for b in r["blocks"].values()
                if b["conv"] and not b["bch"])
    same = (r["converged"] == ref["converged"] == 1056
            and r["bch_ok"] == ref["bch_ok"] and dirty == 0
            and r["jobs"] == ref["jobs"] == 8697)
    print(f"  [6] Fox referee through the port: {r['converged']}/"
          f"{r['jobs']} converged ({r['bch_ok']} BCH zero, {dirty} dirty) "
          f"vs E58's {ref['converged']}  {PASS if same else FAIL}")
    # and the collect/sigma the decode consumed must be the instrument's
    sa = np.load(os.path.join(HERE, "e60_sm_sigma.npy"))
    sb = np.load(os.path.join(HERE, "e58_sm_sigma.npy"))
    g2 = (sa.shape == sb.shape
          and np.array_equal(np.isnan(sa), np.isnan(sb))
          and np.array_equal(sa[~np.isnan(sa)], sb[~np.isnan(sb)]))
    pa = np.load(os.path.join(HERE, "e60_sm_plp0.npy"), mmap_mode="r")
    pb = np.load(os.path.join(HERE, "e58_sm_plp0.npy"), mmap_mode="r")
    g3 = pa.shape == pb.shape and all(
        np.array_equal(np.asarray(pa[f]), np.asarray(pb[f]))
        for f in range(0, pa.shape[0], 7))
    print(f"  [6] ported collect/sigma == banked instrument artifacts: "
          f"sigma {'IDENTICAL' if g2 else 'DIFFER'}, cells(sampled 1/7) "
          f"{'IDENTICAL' if g3 else 'DIFFER'}  {PASS if g2 and g3 else FAIL}")
    return same and g2 and g3


def leg7(state):
    m9_fast, M16, y, t0, fd_off = state
    rng = np.random.default_rng(7)
    fd_w = m9_fast.FrameDecoder(
        threads=6, backend="cpu", cpu_fast=True,
        margin=dict(ce_w=0, wllr="auto", sp=0))
    fd_w.prewarm()
    picked = None
    for nsig in (0.10, 0.14, 0.18, 0.24):
        yn = y + (rng.standard_normal(len(y))
                  + 1j * rng.standard_normal(len(y))) * (
            nsig * np.sqrt(np.mean(np.abs(y[:200000]) ** 2) / 2))
        _, _, dg_p = fd_off.decode_frame(yn, t0)
        _, _, dg_w = fd_w.decode_frame(yn, t0)
        if dg_w.get("snr_db") is not None and dg_w["wllr"]:
            picked = (nsig, dg_p, dg_w)
            if dg_p["converged"] < 74:
                break
    if picked is None:
        print(f"  [7] switch never engaged (noise sweep exhausted)  {FAIL}")
        return False
    nsig, dg_p, dg_w = picked
    ok = dg_w["wllr"] and dg_w["converged"] >= dg_p["converged"] - 1
    print(f"  [7] AWGN-degraded frame (snr {dg_w['snr_db']:.1f} dB < gate): "
          f"switch ENGAGED, weighted {dg_w['converged']}/74 vs plain "
          f"{dg_p['converged']}/74  {PASS if ok else FAIL}")
    fd_w.close()
    return ok


def main():
    print("E60 gates")
    print("=" * 72)
    ok = True
    ok &= bool(leg0())
    ok &= bool(leg1())
    state = _rf33_setup()
    ok &= bool(leg2(state))
    ok &= bool(leg3())
    ok &= bool(leg4(state))
    ok &= bool(leg7(state))
    state[4].close()
    ok &= bool(leg5())
    ok &= bool(leg6())
    print(f"\n  {'ALL E60 GATES PASS' if ok else 'GATE FAILURES'}")
    return 0 if ok else 1


if __name__ == "__main__":
    os.environ.setdefault("M9_NO_TORCH", "1")
    raise SystemExit(main())
