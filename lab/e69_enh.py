#!/usr/bin/env python3
"""E69 -- the Enhanced-layer attempt on RF30, with its own threshold measured.

PLP 0 and PLP 1 share cti_depth 3 (Nrows 1024) and cti_start_row 858, and both
span the identical 1,044,402 cells -- which is what LDM means: the layers are
combined BEFORE the time interleaver, so one deinterleave serves both.  Only
the FEC-block phase and size differ (C from cti_fec_block_start 32966,
10800 cells/block for 64QAM vs 32400 for QPSK).

Two modes:

  selftest  Find what OUR OWN 64QAM-NUC 6/15 chain needs, by encoding random
            payloads, adding calibrated AWGN and sweeping SNR.  Without this a
            live failure is unattributable -- "0 blocks converged" could mean
            the link is short OR that the enhanced pipeline is broken.  This is
            the positive control, and E67's lesson is that a gate whose
            positive control has not been demonstrated proves nothing.

  live      Cancel the decoded Core layer out of the received cells and try to
            decode the Enhanced layer underneath.
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


def chains(js, prefix):
    g = M10.geometry_from_json(js, prefix)
    plps = {p["id"]: p for p in g.plps}
    return plps[0], plps[1]


def selftest(enh, a):
    ch = B.PlpChain(enh["ninner"], enh["mod"], enh["rate"], iters=a.iters)
    rng = np.random.default_rng(20260810)
    print("positive control: synthetic %s %s, %d cells/block"
          % (enh["mod"], enh["rate"], ch.cells_per_fec))
    print("  %-8s %-10s %s" % ("SNR dB", "converged", "BCH zero"))
    results = {}
    for snr in a.sweep:
        conv = bch = 0
        for _ in range(a.trials):
            payload = rng.integers(0, 2, ch.kpayload).astype(np.uint8)
            s = np.asarray(ch.encode(payload), np.complex128)
            p = float(np.mean(np.abs(s) ** 2))
            sigma2 = p / (10 ** (snr / 10.0))
            n = (rng.normal(0, np.sqrt(sigma2 / 2), len(s))
                 + 1j * rng.normal(0, np.sqrt(sigma2 / 2), len(s)))
            r = ch.decode(s + n, sigma2=sigma2)
            if r["converged"]:
                conv += 1
                if ch.baseband_packet(r["bits"])["bch_ok"]:
                    bch += 1
        results[snr] = (conv, bch)
        print("  %-8.1f %-10s %s" % (snr, "%d/%d" % (conv, a.trials),
                                     "%d/%d" % (bch, a.trials)))
    ok = [s for s, (c, b) in results.items() if b == a.trials]
    thr = min(ok) if ok else None
    print("  -> our chain's threshold (all %d trials BCH-clean): %s"
          % (a.trials, ("%.1f dB" % thr) if thr is not None else "not reached"))
    return thr, results


def live(core, enh, a):
    ch_c = B.PlpChain(core["ninner"], core["mod"], core["rate"], iters=a.iters)
    ch_e = B.PlpChain(enh["ninner"], enh["mod"], enh["rate"], iters=a.iters)
    nrows = CTI.nrows_of(core["cti_depth"], bool(core["cti_extended"]))
    C_c = CTI.solve_C(core["cti_fec_block_start"], core["cti_start_row"], nrows)
    C_e = CTI.solve_C(enh["cti_fec_block_start"], enh["cti_start_row"], nrows)
    print("  core  C %d, %d cells/block" % (C_c, core["cells_per_fec"]))
    print("  enh   C %d, %d cells/block" % (C_e, enh["cells_per_fec"]))

    cells = np.array(np.load(os.path.join(HERE, f"{a.prefix}_plp0_deint.npy"),
                             mmap_mode="r")[:a.span], np.complex128)
    # ---- cancel the core layer, block by block --------------------------
    ncc = core["cells_per_fec"]
    cancelled = cells.copy()
    done = 0
    m = 0
    while C_c + (m + 1) * ncc <= len(cells):
        lo = C_c + m * ncc
        r = cells[lo:lo + ncc]
        m += 1
        if not np.any(r):
            continue
        d = ch_c.decode(r)
        if not d["converged"]:
            continue
        bb = ch_c.baseband_packet(d["bits"])
        if not bb["bch_ok"]:
            continue
        s = np.asarray(ch_c.encode(bb["bits"]), np.complex128)
        gain = np.vdot(s, r) / np.vdot(s, s)
        cancelled[lo:lo + ncc] = r - gain * s
        done += 1
    print("  cancelled %d core blocks over %d cells" % (done, len(cells)))
    if not done:
        return None

    # normalise so the enhanced constellation sits at unit average power,
    # using the injection level to predict its share of what is left
    p_res = float(np.mean(np.abs(cancelled[C_c:C_c + done * ncc]) ** 2))
    p_enh = p_res / (1 + a.noise_ratio) if a.noise_ratio else p_res * 0.5
    ref = np.asarray(ch_e.encode(np.zeros(ch_e.kpayload, np.uint8)),
                     np.complex128)
    p_ref = float(np.mean(np.abs(ref) ** 2))
    scale = np.sqrt(p_ref / p_enh)
    sigma2 = (p_res - p_enh) * scale ** 2
    print("  P_residual %.5f, assumed P_enh %.5f -> sigma2 %.5f"
          % (p_res, p_enh, sigma2))

    nce = enh["cells_per_fec"]
    conv = bch = tried = 0
    k = 0
    while C_e + (k + 1) * nce <= len(cancelled) and tried < a.blocks:
        lo = C_e + k * nce
        seg = cancelled[lo:lo + nce] * scale
        k += 1
        if not np.any(seg):
            continue
        tried += 1
        d = ch_e.decode(seg, sigma2=sigma2)
        if d["converged"]:
            conv += 1
            if ch_e.baseband_packet(d["bits"])["bch_ok"]:
                bch += 1
    print("  ENHANCED LAYER: %d/%d LDPC converged, %d BCH zero"
          % (conv, tried, bch))
    return dict(tried=tried, converged=conv, bch=bch)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("selftest", "live", "both"))
    ap.add_argument("--prefix", default="e69_smoke")
    ap.add_argument("--l1", default="e69_l1_rf30.json")
    ap.add_argument("--iters", type=int, default=60)
    ap.add_argument("--trials", type=int, default=3)
    ap.add_argument("--blocks", type=int, default=24)
    ap.add_argument("--span", type=int, default=2_000_000)
    ap.add_argument("--noise-ratio", type=float, default=1.074,
                    help="P_noise/P_enh from e69_ldm.py")
    ap.add_argument("--sweep", type=float, nargs="*",
                    default=[4, 6, 7, 8, 9, 10, 12, 15])
    a = ap.parse_args()
    js = json.load(open(os.path.join(HERE, a.l1)))
    core, enh = chains(js, a.prefix)
    out = {}
    if a.mode in ("selftest", "both"):
        thr, res = selftest(enh, a)
        out["threshold_db"] = thr
        out["sweep"] = {str(k): v for k, v in res.items()}
        print()
    if a.mode in ("live", "both"):
        out["live"] = live(core, enh, a)
    json.dump(out, open(os.path.join(HERE, "e69_enh.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
