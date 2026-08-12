#!/usr/bin/env python3
"""E76 -- localise WHERE L1 fails on the new RF25 captures.

m8_l1 verifies no L1-Basic in 72 probes across two new captures, while the
same command verifies on E49's rf25_amped_0808.cs16.  The bootstrap detector
meanwhile scores the NEW captures HIGHER (0.93 vs 0.88, peak-ratio 45 vs 44)
with a cleaner shoulder -- but E64's own correction says bootstrap peak-ratio
carries ~20 dB more margin than L1, so "bootstrap is strong" is not "L1 should
decode".  This walks the stages between them and prints, for each capture,
what the preamble acquisition actually got:

    fft_window_start   did it find a Frame at all?
    pilot coherence    is the preamble symbol demodulating?
    l1basic_of         does L1-Basic's own CRC pass?

The old capture is the positive control in every run.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

import m3_freqint as FI       # noqa: E402
import m4_l1detail as D4      # noqa: E402
import m3_spec as S           # noqa: E402
from m3_preamble import analyse  # noqa: E402

CAPS = [
    ("data/fox25/rf25_amped_0808.cs16", "E49 CONTROL (known good)", [0.5]),
    ("data/fox25/rf25_fox_g2_0810_2257.cs16", "NEW gain 2", None),
    ("data/fox25/rf25_fox_g4_0810_2254.cs16", "NEW gain 4", None),
]
RATE = 6.912e6


def probe(path, start_sec):
    rep = {}
    try:
        _r, z0, _Y, geo = analyse(path, RATE, None, report=rep,
                                  start_sec=start_sec, quiet=True)
    except Exception as exc:                                    # noqa: BLE001
        return dict(start=start_sec, err="%s: %s" % (type(exc).__name__, exc))
    (lo, n, pilot, cp, data, shift, cred, nfft, gi, dx, mode) = geo
    out = dict(start=start_sec, nfft=nfft, gi=gi, dx=dx, cred=cred,
               mode=mode, shift=shift,
               t0=rep.get("fft_window_start"),
               cells=int(len(z0)))
    for k in ("pilot_coherence", "coherence", "data_coherence_mean",
              "peak_ratio", "corr_peak"):
        if k in rep:
            out[k] = rep[k]
    x0 = FI.deinterleave(z0, nfft, 0, direction="forward", toggle="i")
    f = D4.l1basic_of(x0, mode)
    out["l1basic"] = "OK" if f else "CRC FAIL"
    if f:
        out["fft_size"] = f.get("L1B_first_sub_fft_size")
        out["num_sub"] = f.get("L1B_num_subframes")
        out["nsym"] = f.get("L1B_first_sub_num_ofdm_symbols")
    return out


def main():
    for rel, label, times in CAPS:
        p = os.path.join(ROOT, rel)
        print("=" * 78)
        print("%-34s %s" % (os.path.basename(rel), label))
        ts = times if times else [0.0, 0.121, 0.243, 0.364, 0.486, 1.0, 5.0,
                                  20.0, 60.0]
        for t in ts:
            r = probe(p, t)
            if "err" in r:
                print("  t=%6.3f  ERROR %s" % (t, r["err"][:70]))
                continue
            coh = r.get("pilot_coherence", r.get("coherence", float("nan")))
            print("  t=%6.3f  nfft=%5d gi=%4d dx=%d cred=%d mode=%s  t0=%-9s "
                  "cells=%5d  coh=%s  L1-Basic %s"
                  % (t, r["nfft"], r["gi"], r["dx"], r["cred"], r["mode"],
                     r["t0"], r["cells"],
                     ("%.4f" % coh) if coh == coh else "n/a", r["l1basic"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
