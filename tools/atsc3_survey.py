#!/usr/bin/env python3
"""atsc3_survey.py -- one pass over the UHF band: WHO transmits ATSC 3.0?

The stress test's first question is "all the channels" -- and the honest
first step is enumerating what exists. For each RF channel this tunes,
captures ~1.6 s, and runs BOTH bootstrap detectors (structural + matched --
the anti-self-consistency pair from the 8/05 campaign; if they disagree the
fault is our spec reading, not the air). The bootstrap is designed for
detection at -6 dB, so this works even in fade conditions where the payload
is undecodable -- which is exactly when the survey is cheapest to run,
because it is not stealing watchable TV time.

Output: one line per channel + survey.json. 8-VSB (ATSC 1.0) carriers show
as no-bootstrap; the STVT side already knows them.

Radio discipline: single-tenant -- do NOT run while the chain owns the SDR.
Honors --force-meteor for explicit user override only.

Usage:
    python tools/atsc3_survey.py --channels 14-36 --force-meteor
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "lab"))
sys.path.insert(0, ROOT)

import SoapySDR                                                   # noqa: E402
from SoapySDR import SOAPY_SDR_RX, SOAPY_SDR_CS16                 # noqa: E402
import m2_pilots as MP                                            # noqa: E402
from atsc3 import bootstrap as bs                                 # noqa: E402

SoapySDR.SoapySDR_setLogLevel(SoapySDR.SOAPY_SDR_FATAL)
FS = 6.912e6


def center_hz(rf):
    lo = (174 + (rf - 7) * 6) if rf < 14 else (470 + (rf - 14) * 6)
    return (lo + 3.0) * 1e6


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def in_meteor_window(now=None):
    now = now or dt.datetime.now()
    mins = now.hour * 60 + now.minute
    return 2 * 60 + 10 <= mins <= 5 * 60 + 40


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--channels", default="14-36")
    ap.add_argument("--secs", type=float, default=1.6)
    ap.add_argument("--ant", default="Antenna B")
    ap.add_argument("--out", default=os.path.join(ROOT, "data",
                                                  "survey.json"))
    ap.add_argument("--force-meteor", action="store_true")
    a = ap.parse_args()
    if in_meteor_window() and not a.force_meteor:
        log("REFUSING: inside the 02:10-05:40 meteor window.")
        return 3
    lo, hi = (int(x) for x in a.channels.split("-"))
    chans = list(range(lo, hi + 1))

    sdr = SoapySDR.Device("driver=sdrplay")
    sdr.setSampleRate(SOAPY_SDR_RX, 0, FS)
    time.sleep(0.2)
    fs = float(sdr.getSampleRate(SOAPY_SDR_RX, 0))
    sdr.setAntenna(SOAPY_SDR_RX, 0, a.ant)
    try:
        sdr.setGainMode(SOAPY_SDR_RX, 0, False)
    except Exception:                                          # noqa: BLE001
        pass
    sdr.setGain(SOAPY_SDR_RX, 0, "IFGR", 32)
    try:
        sdr.writeSetting("rfgain_sel", "2")
    except Exception:                                          # noqa: BLE001
        pass
    st = sdr.setupStream(SOAPY_SDR_RX, SOAPY_SDR_CS16)
    sdr.activateStream(st)
    buf = np.empty(2 * 65536, np.int16)
    results = []
    log(f"surveying RF{lo}-RF{hi} at {fs/1e6:.3f} Msps on {a.ant}")
    try:
        for rf in chans:
            sdr.setFrequency(SOAPY_SDR_RX, 0, center_hz(rf))
            time.sleep(0.25)
            # flush the retune transient
            t0 = time.time()
            while time.time() - t0 < 0.2:
                sdr.readStream(st, [buf], 65536, timeoutUs=300000)
            n_want = int(a.secs * fs)
            iq = np.empty(n_want, np.complex64)
            got = 0
            deadline = time.time() + a.secs * 2 + 3
            while got < n_want and time.time() < deadline:
                r = sdr.readStream(st, [buf], 65536, timeoutUs=500000)
                if r.ret > 0:
                    n = min(r.ret, n_want - got)
                    x = buf[:2 * n].astype(np.float32).view(np.complex64)
                    iq[got:got + n] = x / 32768.0
                    got += n
            y = MP.resample_to(iq[:got], fs, bs.FS)
            hits = MP.find_bootstraps(y)
            in_rms = float(np.sqrt((np.abs(iq[:got]) ** 2).mean()) * 32768)
            if hits:
                h = hits[0]
                pr = h.get("mean_peak_ratio", 0.0)
                results.append(dict(rf=rf, atsc3=True, peak_ratio=pr,
                                    cfo=h.get("fine_cfo_hz"), rms=in_rms,
                                    n_hits=len(hits)))
                log(f"  RF{rf}: ATSC 3.0 BOOTSTRAP  peak {pr:.1f}  "
                    f"({len(hits)} hits)  rms {in_rms:.0f}")
            else:
                results.append(dict(rf=rf, atsc3=False, rms=in_rms))
                log(f"  RF{rf}: no bootstrap  rms {in_rms:.0f}")
    finally:
        try:
            sdr.deactivateStream(st)
            sdr.closeStream(st)
        except Exception:                                      # noqa: BLE001
            pass
        del sdr
    json.dump(dict(when=time.strftime("%Y-%m-%d %H:%M"),
                   ant=a.ant, secs=a.secs, results=results),
              open(a.out, "w"), indent=1)
    n3 = sum(1 for r in results if r.get("atsc3"))
    log(f"survey complete: {n3} ATSC 3.0 carrier(s) of {len(chans)} channels "
        f"-> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
