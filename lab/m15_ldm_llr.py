#!/usr/bin/env python3
"""M15 -- alphabet-aware LLRs against the LDM interferer: is the last lever real?

M14's budget measured every candidate for the missing dB on the non-decoding
multiplex and found all but one worth zero.  The survivor:

    the Core layer's demapper treats the Enhanced layer as GAUSSIAN NOISE,
    and it is not -- it is a known finite constellation at a known relative
    power, so the LLRs could marginalise over its alphabet instead.

This measures that, and it is deliberately built so it can come back ZERO.

THE PREDICTION, STATED BEFORE THE RESULT
-----------------------------------------
**I expect a SMALL gain here, possibly none.**  The interferer on this channel
is 256QAM-NUC: 256 points, dense, and a dense constellation's sum is very close
to Gaussian -- which is exactly the condition under which the Gaussian
approximation is already right.  Alphabet-awareness should pay when the
interferer is LOW order and lumpy, not when it is dense.

THE CONTROL THAT MAKES A NULL RESULT MEAN SOMETHING
----------------------------------------------------
A null is only informative if the alphabet-aware demapper is known to work, so
the same comparison is run with a QPSK interferer -- 4 points, maximally lumpy,
where the Gaussian approximation is definitely wrong.  If alphabet-awareness
wins there and not on 256QAM, the mechanism is confirmed and a null on the real
channel means "no gain available", not "the code is broken".  If it wins
NEITHER, the implementation is suspect and the result is not reported as a
finding.

WHAT IS SYNTHESISED
-------------------
The Enhanced layer's own FEC is irrelevant to the Core demapper -- from the
Core layer's point of view its cells are draws from a constellation, and that
is exactly how a real Enhanced layer looks to it.  So the interferer is uniform
random points of its alphabet, which is honest rather than convenient: giving
it real coding could only make it MORE Gaussian-looking, never less.

    y = sqrt(P_core)*x_core + sqrt(P_enh)*x_enh + n
    P_core/P_enh = injection level,  P_core + P_enh = 1,  var(n) = N0

so "channel SNR" below is 1/N0 in dB, directly comparable to the 3.37 dB
measured off the dummy cells on air.

THE TWO DEMAPPERS
-----------------
  gaussian   the shipping one: fold P_enh into the noise and run max-log
             against the Core alphabet alone.
  alphabet   max-log over the JOINT, i.e. minimise over the interferer's
             alphabet rather than averaging it into a variance.

Both feed the same `decode_llr`, so the LDPC, BCH and descrambler are shared
and the only difference in the experiment is the LLR vector.

Usage:
    python m15_ldm_llr.py                 # the real operating point + control
    python m15_ldm_llr.py --control-only
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

import m6_bicm as B                                              # noqa: E402
import spec_bitint as BI                                         # noqa: E402


def powers(inj_db):
    r = 10 ** (-inj_db / 10.0)
    return 1 / (1 + r), r / (1 + r)


def llr_gaussian(y, chain, p_core, p_enh, n0):
    """The shipping demapper: the interferer becomes part of the variance."""
    return B.demap_llr(y / np.sqrt(p_core), chain.mod, chain.rate,
                       sigma2=(p_enh + n0) / p_core)


def llr_alphabet(y, chain, p_core, p_enh, n0, enh_pts, chunk=4096):
    """Max-log over the joint (core symbol, interferer symbol).

    For each transmitted-bit hypothesis the metric is the SMALLEST distance
    over the interferer's alphabet, instead of a distance inflated by its
    average power.  Same max-log convention as `demap_llr`, so the two are
    comparable by construction; `llr > 0` means bit 0.
    """
    pts = B.points_for(chain.mod, chain.rate)
    nb = BI.MOD_BITS[chain.mod]
    ones = B.bit_masks(nb)
    y = np.asarray(y)
    out = np.empty((len(y), nb))
    a, b = np.sqrt(p_core) * pts, np.sqrt(p_enh) * enh_pts
    grid = (a[:, None] + b[None, :]).ravel()          # every joint hypothesis
    lab = np.repeat(np.arange(len(pts)), len(enh_pts))
    for lo in range(0, len(y), chunk):
        z = y[lo:lo + chunk]
        d2 = np.abs(z[:, None] - grid[None, :]) ** 2
        # collapse the interferer: best case over its alphabet, per core point
        best = np.full((len(z), len(pts)), np.inf)
        for k in range(len(pts)):
            best[:, k] = d2[:, lab == k].min(1)
        for i in range(nb):
            out[lo:lo + chunk, i] = (best[:, ones[i]].min(1)
                                     - best[:, ~ones[i]].min(1)) / n0
    return out.ravel()


def llr_exact(y, chain, p_core, p_enh, n0, enh_pts, chunk=2048):
    """TRUE marginalisation over the interferer -- log-sum-exp, not min.

    WHY THIS EXISTS: `llr_alphabet` marginalises with max-log, i.e. it takes the
    BEST interferer point for each hypothesis.  With a dense interferer that is
    a bad approximation and an actively harmful one: every core hypothesis finds
    some interferer symbol that explains the observation, so the LLR magnitudes
    collapse and the LDPC loses the reliability information it runs on.  That is
    an artefact of max-log, NOT a fact about the channel, so the fact has to be
    measured with the marginalisation done properly:

        LLR = logsumexp_{bit=0}(-d^2/N0) - logsumexp_{bit=1}(-d^2/N0)

    which is the exact bit LLR under the true (finite-alphabet) interference
    model, and reduces to `llr_alphabet` when one term dominates each sum.
    """
    from scipy.special import logsumexp
    pts = B.points_for(chain.mod, chain.rate)
    nb = BI.MOD_BITS[chain.mod]
    ones = B.bit_masks(nb)
    y = np.asarray(y)
    out = np.empty((len(y), nb))
    a, b = np.sqrt(p_core) * pts, np.sqrt(p_enh) * enh_pts
    grid = (a[:, None] + b[None, :]).ravel()
    lab = np.repeat(np.arange(len(pts)), len(enh_pts))
    for lo in range(0, len(y), chunk):
        z = y[lo:lo + chunk]
        m = -np.abs(z[:, None] - grid[None, :]) ** 2 / n0
        for i in range(nb):
            sel1 = np.isin(lab, np.flatnonzero(ones[i]))
            out[lo:lo + chunk, i] = (logsumexp(m[:, ~sel1], axis=1)
                                     - logsumexp(m[:, sel1], axis=1))
    return out.ravel()


def trial(core_mod, core_rate, enh_mod, enh_rate, inj_db, snr_db, seed,
          mode, iters=50, ninner=64800):
    """One synthesised codeword through the composite channel. -> bch_ok."""
    rng = np.random.default_rng(seed)
    ch = B.PlpChain(ninner, core_mod, core_rate, iters=iters)
    pay = rng.integers(0, 2, ch.kpayload).astype(np.uint8)
    x_core = ch.encode(pay)
    enh_pts = B.points_for(enh_mod, enh_rate)
    x_enh = enh_pts[rng.integers(0, len(enh_pts), len(x_core))]
    p_core, p_enh = powers(inj_db)
    n0 = 10 ** (-snr_db / 10.0)
    n = (rng.standard_normal(len(x_core))
         + 1j * rng.standard_normal(len(x_core))) / np.sqrt(2)
    y = (np.sqrt(p_core) * x_core + np.sqrt(p_enh) * x_enh
         + n * np.sqrt(n0))
    if mode == "gaussian":
        q = llr_gaussian(y, ch, p_core, p_enh, n0)
    elif mode == "alphabet":
        q = llr_alphabet(y, ch, p_core, p_enh, n0, enh_pts)
    else:
        q = llr_exact(y, ch, p_core, p_enh, n0, enh_pts)
    lam = np.empty(ch.ninner)
    lam[ch.lam_of_q] = q
    r = ch.decode_llr(lam)
    if not r["converged"]:
        return False
    bb = ch.baseband_packet(r["bits"])
    return bool(bb["bch_ok"]) and bool(np.array_equal(bb["bits"], pay))


def threshold(core, enh, inj_db, mode, lo, hi, step=0.25, trials=3,
              verbose=True):
    """Lowest SNR on the grid where all `trials` succeed -- the 3/3 convention
    M9 used for every other threshold in this project, kept so the numbers are
    comparable to the ones already recorded."""
    snr = lo
    while snr <= hi + 1e-9:
        ok = all(trial(core[0], core[1], enh[0], enh[1], inj_db, snr,
                       1000 + k, mode) for k in range(trials))
        if verbose:
            print(f"      {mode:8s} SNR {snr:5.2f} dB  "
                  f"{'PASS 3/3' if ok else 'fail'}", flush=True)
        if ok:
            return snr
        snr += step
    return None


def compare(name, core, enh, inj_db, lo, hi, step, trials, verbose=True):
    print(f"\n  --- {name} ---")
    print(f"    Core {core[0]} {core[1]}   interferer {enh[0]} {enh[1]}   "
          f"injection {inj_db:.1f} dB")
    out = {}
    for mode in ("gaussian", "alphabet", "exact"):
        t0 = time.time()
        out[mode] = threshold(core, enh, inj_db, mode, lo, hi, step, trials,
                              verbose)
        print(f"      -> {mode} threshold "
              f"{('%.2f dB' % out[mode]) if out[mode] is not None else 'NOT REACHED'}"
              f"   ({time.time()-t0:.0f}s)")
    g = out.get("gaussian")
    for mode in ("alphabet", "exact"):
        if g is not None and out.get(mode) is not None:
            out["gain_" + mode] = g - out[mode]
            print(f"    GAIN, {mode:8s} vs gaussian: "
                  f"{out['gain_' + mode]:+.2f} dB")
        else:
            out["gain_" + mode] = None
            print(f"    GAIN, {mode}: not computable -- threshold not reached")
    out["gain_db"] = out.get("gain_exact")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(prog="m15 ldm llr")
    ap.add_argument("--inj", type=float, default=5.0)
    ap.add_argument("--step", type=float, default=0.25)
    ap.add_argument("--trials", type=int, default=3)
    ap.add_argument("--lo", type=float, default=1.0)
    ap.add_argument("--hi", type=float, default=9.0)
    ap.add_argument("--control-only", action="store_true")
    ap.add_argument("--real-only", action="store_true")
    ap.add_argument("--json", default="m15_ldm_llr.json")
    a = ap.parse_args(argv)

    print("M15 -- alphabet-aware LLRs vs the LDM interferer")
    print("=" * 72)
    print("  PREDICTION (before the result): SMALL or ZERO on the real "
          "channel,\n  because a 256-point interferer is already nearly "
          "Gaussian.  The QPSK\n  control exists so that a null means 'no gain "
          "available' rather than\n  'the demapper is broken'.")
    res = {}
    if not a.real_only:
        res["control_qpsk_interferer"] = compare(
            "CONTROL: lumpy interferer (QPSK), where the Gaussian "
            "approximation is definitely wrong",
            ("QPSK", "9/15"), ("QPSK", "9/15"), a.inj, a.lo, a.hi, a.step,
            a.trials)
    if not a.control_only:
        res["real_256qam_interferer"] = compare(
            "THE REAL OPERATING POINT: 256QAM-NUC interferer",
            ("QPSK", "9/15"), ("256QAM", "7/15"), a.inj, a.lo, a.hi, a.step,
            a.trials)

    print("\n" + "=" * 72)
    c = res.get("control_qpsk_interferer", {}).get("gain_db")
    r = res.get("real_256qam_interferer", {}).get("gain_db")
    if c is not None:
        print(f"  control (QPSK interferer)   alphabet-aware gain "
              f"{c:+.2f} dB")
    if r is not None:
        print(f"  real    (256QAM interferer) alphabet-aware gain "
              f"{r:+.2f} dB")
    if c is not None and r is not None:
        if c <= 0.0:
            print("\n  CONTROL DID NOT WIN -- the demapper is not demonstrably "
                  "correct, so\n  the real-channel number is NOT reported as a "
                  "finding.")
        elif r < 0.25:
            print("\n  The mechanism is real (the control shows it) and it is "
                  "NOT AVAILABLE\n  on this channel: the interferer is too "
                  "dense to be worth modelling.\n  The antenna remains the "
                  "only lever.")
        else:
            print("\n  A real gain on the real channel -- worth building into "
                  "the decoder.")
    with open(os.path.join(HERE, a.json), "w") as fh:
        json.dump(res, fh, indent=1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
