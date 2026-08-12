#!/usr/bin/env python3
"""E69 -- LDM core cancellation, and the honest Enhanced-layer margin.

RF30's television services are NOT on the Core layer.  The Core layer (PLP 0,
QPSK 6/15) decodes 740/740 BCH-clean and carries only:

    239.255.58.4:5004    RTP/MP2T datacast ("flexstream")
    239.255.14.103:6123  broadSpan RTK GNSS corrections
    224.0.23.60:4937     LLS

while the SLT recovered from that same Core layer places both TV services
elsewhere: 158.1 WIAV at 239.255.58.1:5001 and 158.5 24/7MMT at
239.255.58.5:5006.  Neither flow appears.  So they ride PLP 1 -- the ENHANCED
layer, 64QAM-NUC 6/15, injected 4.0 dB BELOW the core.

m10_core's docstring says why the campaign never built cancellation: "nothing
asked for is there."  Something is asked for now.

WHAT THIS MEASURES
------------------
The Core layer is decoded and re-ENCODED exactly (m6_bicm.PlpChain.encode is a
true round trip: payload -> scramble -> BCH -> LDPC -> bit interleave ->
constellation).  Subtracting it from the received cells leaves

    residual = enhanced + noise

and the injection level closes the system:  P_enh = P_core / 10^(inj/10).
So the Enhanced layer's own SNR is recoverable without ever decoding it:

    P_noise    = P_residual - P_enh
    SNR_enh dB = 10 log10(P_enh / P_noise)

which is the number that decides whether 158.5 is reachable from this capture
or whether the antenna has to move.

THE CONTROL THAT MUST FAIL
--------------------------
Re-encoding a WRONG payload (bits flipped) must leave a residual at full
signal power -- i.e. cancellation must do nothing.  If a wrong codeword
cancelled as well as the right one, the "cancellation" would be fitting noise
and every SNR below would be fiction.
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


def db(x):
    return 10 * np.log10(x) if x > 0 else float("-inf")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("prefix", nargs="?", default="e69_smoke")
    ap.add_argument("--l1", default="e69_l1_rf30.json")
    ap.add_argument("--blocks", type=int, default=48)
    ap.add_argument("--iters", type=int, default=60)
    a = ap.parse_args()

    js = json.load(open(os.path.join(HERE, a.l1)))
    g = M10.geometry_from_json(js, a.prefix)
    plps = {p["id"]: p for p in g.plps}
    core, enh = plps[0], plps[1]
    inj_db = 4.0                      # L1D_plp_ldm_injection_level 8 -> 4.0 dB

    print("RF30 LDM layers, from the decoded L1-Detail")
    print("  core     PLP %d  %-9s %-5s  Ninner %d" %
          (core["id"], core["mod"], core["rate"], core["ninner"]))
    print("  enhanced PLP %d  %-9s %-5s  Ninner %d  injection %.1f dB" %
          (enh["id"], enh["mod"], enh["rate"], enh["ninner"], inj_db))

    ncell = core["cells_per_fec"]
    nrows = CTI.nrows_of(core["cti_depth"], bool(core["cti_extended"]))
    C0 = CTI.solve_C(core["cti_fec_block_start"], core["cti_start_row"], nrows)
    cells = np.load(os.path.join(HERE, f"{a.prefix}_plp0_deint.npy"),
                    mmap_mode="r")

    ch = B.PlpChain(core["ninner"], core["mod"], core["rate"], iters=a.iters)
    rows = []
    ctl = []
    for m in range(a.blocks):
        lo = C0 + m * ncell
        r = np.asarray(cells[lo:lo + ncell], np.complex128)
        if not np.any(r):
            continue
        d = ch.decode(r)
        if not d["converged"]:
            continue
        bb = ch.baseband_packet(d["bits"])
        if not bb["bch_ok"]:
            continue
        s = np.asarray(ch.encode(bb["bits"]), np.complex128)
        if len(s) != len(r):
            print("  !! re-encode length %d != %d" % (len(s), len(r)))
            break
        # complex gain fit, then cancel
        gain = np.vdot(s, r) / np.vdot(s, s)
        res = r - gain * s
        p_core = float(abs(gain) ** 2 * np.mean(np.abs(s) ** 2))
        p_res = float(np.mean(np.abs(res) ** 2))
        p_tot = float(np.mean(np.abs(r) ** 2))
        rows.append((p_tot, p_core, p_res))

        # ---- negative control: cancel with a CORRUPTED payload -----------
        bad = bb["bits"].copy()
        bad[::7] ^= 1
        s_bad = np.asarray(ch.encode(bad), np.complex128)
        gb = np.vdot(s_bad, r) / np.vdot(s_bad, s_bad)
        ctl.append(float(np.mean(np.abs(r - gb * s_bad) ** 2)))

    if not rows:
        print("  no clean blocks -- nothing to measure")
        return 1
    A = np.array(rows)
    p_tot, p_core, p_res = A[:, 0].mean(), A[:, 1].mean(), A[:, 2].mean()
    p_ctl = float(np.mean(ctl))
    p_enh = p_core / (10 ** (inj_db / 10.0))
    p_noise = p_res - p_enh

    print()
    print("  blocks measured: %d" % len(rows))
    print("  P_total    %.5f" % p_tot)
    print("  P_core     %.5f   (%.2f dB of total)" % (p_core, db(p_core / p_tot)))
    print("  P_residual %.5f   (%.2f dB of total)  <- enhanced + noise"
          % (p_res, db(p_res / p_tot)))
    print("  cancellation gain: %.2f dB" % db(p_tot / p_res))
    print("  NEGATIVE CONTROL (corrupted payload re-encoded):")
    print("    P_residual %.5f  -> cancellation gain %.2f dB  %s"
          % (p_ctl, db(p_tot / p_ctl),
             "(FAILS to cancel, as required)" if p_ctl > 0.8 * p_tot
             else "!! CONTROL DID NOT FAIL -- measurement is unsound"))
    print()
    print("  core-layer SNR   (P_core / P_residual)     %6.2f dB" %
          db(p_core / p_res))
    if p_noise <= 0:
        print("  P_noise <= 0: the residual is at or below the predicted "
              "enhanced-layer power -- noise is below the measurement floor")
        snr_enh = float("inf")
    else:
        snr_enh = db(p_enh / p_noise)
        print("  P_enh (predicted from injection) %.5f" % p_enh)
        print("  P_noise (residual - enhanced)    %.5f" % p_noise)
        print("  ENHANCED-LAYER SNR               %6.2f dB" % snr_enh)

    # 64QAM-NUC 6/15 threshold, MEASURED on our own chain by e69_enh.py
    # selftest (0/3 BCH-clean at 8.0 dB, 3/3 at 9.0 dB) -- not a spec number
    # we looked up, so a live failure is attributable to the link and not to
    # an assumed threshold.
    need = 9.0
    print()
    print("  64QAM-NUC 6/15 needs %.1f dB -- MEASURED on our chain "
          "(e69_enh.py selftest), not assumed" % need)
    if snr_enh == float("inf"):
        print("  -> margin cannot be bounded from above; attempt the decode")
    else:
        print("  -> margin %+.2f dB  (%s)"
              % (snr_enh - need,
                 "should decode" if snr_enh >= need else "SHORT"))
    json.dump(dict(p_total=p_tot, p_core=p_core, p_residual=p_res,
                   p_control=p_ctl, p_enh=p_enh, p_noise=p_noise,
                   snr_core_db=db(p_core / p_res),
                   snr_enh_db=(None if snr_enh == float("inf") else snr_enh),
                   blocks=len(rows), inj_db=inj_db, need_db=need),
              open(os.path.join(HERE, "e69_ldm.json"), "w"), indent=1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
