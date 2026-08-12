#!/usr/bin/env python3
"""E76 -- L1-Basic verification RATE as a link-margin proxy on RF25.

L1-Basic is the most protected payload in the Frame after the bootstrap, and
it either CRC-passes or it does not -- there is no partial credit and no
threshold to argue about.  So the fraction of probed Frames whose L1-Basic
verifies is a blunt but honest margin instrument, and unlike the bootstrap
peak-ratio (which E64 showed carries ~20 dB more margin than L1 and therefore
flatters a dying link) it sits at the layer we actually care about.

Probes each capture at evenly spaced times across its whole duration and
reports the rate, plus the L1-Detail LDPC/BCH outcome where L1-Basic passed.
E49's capture is included in every run as the control.
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

import m3_freqint as FI          # noqa: E402
import m4_l1detail as D4         # noqa: E402
import m8_l1 as M8               # noqa: E402
from m3_preamble import analyse  # noqa: E402

RATE = 6.912e6
CAPS = [
    ("data/fox25/rf25_amped_0808.cs16", "E49 0808 control (no splitter, max gain)"),
    ("data/fox25/rf25_fox_g2_0810_2257.cs16", "0810 gain 2 (splitter)"),
    ("data/fox25/rf25_fox_g4_0810_2254.cs16", "0810 gain 4 (splitter)"),
]


def one(path, t):
    rep = {}
    try:
        _r, z0, _Y, geo = analyse(path, RATE, None, report=rep,
                                  start_sec=t, quiet=True)
    except Exception:                                           # noqa: BLE001
        return None, None, None
    nfft, mode = geo[7], geo[10]
    x0 = FI.deinterleave(z0, nfft, 0, direction="forward", toggle="i")
    f = D4.l1basic_of(x0, mode)
    coh = rep.get("pilot_coherence", rep.get("coherence"))
    if f is None:
        return False, None, coh
    try:
        r = M8.solve(path, RATE, None, t)
        ok = bool(r and r["fec"]["converged"] and r["fec"]["bch_ok"])
    except Exception:                                           # noqa: BLE001
        ok = None
    return True, ok, coh


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--probes", type=int, default=24)
    a = ap.parse_args()
    out = {}
    for rel, label in CAPS:
        p = os.path.join(ROOT, rel)
        if not os.path.exists(p):
            continue
        dur = os.path.getsize(p) / 4 / RATE
        ts = np.linspace(0.5, dur - 1.0, a.probes)
        okb = okd = n = 0
        cohs = []
        for t in ts:
            b, d, coh = one(p, float(t))
            if b is None:
                continue
            n += 1
            okb += bool(b)
            okd += bool(d)
            if coh is not None:
                cohs.append(coh)
        out[os.path.basename(rel)] = dict(
            label=label, probes=n, l1basic_ok=okb, l1detail_ok=okd,
            rate=okb / n if n else None,
            coh_median=float(np.median(cohs)) if cohs else None)
        print("%-34s %-42s  L1-Basic %2d/%2d = %5.1f%%   L1-Detail %2d/%2d   "
              "coh med %s"
              % (os.path.basename(rel), label, okb, n, 100 * okb / max(n, 1),
                 okd, n,
                 ("%.4f" % np.median(cohs)) if cohs else "n/a"))
    json.dump(out, open(os.path.join(HERE, "e76_l1_rate.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
