#!/usr/bin/env python3
"""e87_antenna_sweep.py -- sweep the whole TV band on ONE antenna port.

E87 (8/11).  A Philips antenna went onto ANTENNA PORT A, which the lab's
recorded topology still calls the K-180WLA loop -- so every older port-A
number describes a different antenna and cannot be compared to these.

*** E89 RETRACTION: the Philips was NOT ACTUALLY CONNECTED. ***  The E87
run's port-A column is a measurement of an OPEN COAX PORT and is void as
an antenna result.  It is retained deliberately as a DISCONNECTED
BASELINE (sweep_A.json), which makes the window-4 re-sweep a real A/B:
same tool, same settings, same evening, same port, one variable.

The trap this tool did NOT catch: an open port does not read zero.  It
returned rms 90-340 across the band and a 3/3-round bootstrap detection
on RF33 at peak_ratio 37.4.  The antenna-readback guard below confirms
WHICH PORT the tuner selected; it cannot confirm that anything is
attached to that port.  A readback proves the switch moved, not that the
coax exists.  See CONNECTEDNESS TEST in E89: count distinct carriers,
do not read levels.

WHAT IT MEASURES, per RF channel:
  in_rms        broadband level in the 6 MHz window (a LEVEL, not a quality)
  peak_ratio    A/321 bootstrap correlation -- the ATSC 3.0 detector, and a
                QUALITY measure: it rises with SNR and is meaningless for
                8-VSB, which carries no bootstrap at all
  atsc3         did the gated detector actually find a bootstrap
  cfo_hz        carrier offset, a sanity check that we are on a real signal

WHY BOTH NUMBERS.  E64 measured a case where bias-T left rms UNCHANGED
(1.00x) while peak_ratio moved 64.6 -> 81.0: level and quality are
independent axes and an antenna comparison that reads only one of them will
confidently rank the wrong port.  A dumb antenna pointed at a strong
transmitter can beat a good one pointed away on rms and lose on peak_ratio,
and it is peak_ratio that predicts whether anything decodes.

INTERLEAVING.  Channels are visited in a fixed order and the whole sweep is
repeated `--rounds` times.  A single pass ranks channels partly by when they
were sampled, which on a drifting carrier is how E86 nearly banked a false
gain table.  The median across rounds is the reading; the spread is reported
beside it so a comparison has to earn its winner.

    python lab/e87_antenna_sweep.py --ant "Antenna A" --json out.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import SoapySDR
from SoapySDR import SOAPY_SDR_RX, SOAPY_SDR_CS16

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m2_pilots as MP                                       # noqa: E402
from atsc3 import bootstrap as BS                            # noqa: E402

SoapySDR.SoapySDR_setLogLevel(SoapySDR.SOAPY_SDR_FATAL)
FS = 6.912e6

# US TV after the repack: VHF-hi 7-13, UHF 14-36.  VHF-lo (2-6) is omitted --
# no local ATSC 3.0 lives there and the RSPdx's own noise dominates.
CHANS = list(range(7, 14)) + list(range(14, 37))


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def center_hz(rf):
    lo = (174 + (rf - 7) * 6) if rf < 14 else (470 + (rf - 14) * 6)
    return (lo + 3.0) * 1e6


def measure(sdr, st, buf, rf, secs):
    sdr.setFrequency(SOAPY_SDR_RX, 0, center_hz(rf))
    time.sleep(0.25)
    t0 = time.time()
    while time.time() - t0 < 0.2:                    # flush retune transient
        sdr.readStream(st, [buf], 65536, timeoutUs=300000)
    n_want = int(secs * FS)
    iq = np.empty(n_want, np.complex64)
    got = 0
    deadline = time.time() + secs * 2 + 3
    while got < n_want and time.time() < deadline:
        r = sdr.readStream(st, [buf], 65536, timeoutUs=500000)
        if r.ret > 0:
            n = min(r.ret, n_want - got)
            iq[got:got + n] = (buf[:2 * n].astype(np.float32)
                               .view(np.complex64)) / 32768.0
            got += n
    if got < 1000:
        return dict(rf=rf, rms=0.0, peak_ratio=0.0, atsc3=False, cfo=None)
    x = iq[:got]
    rms = float(np.sqrt((np.abs(x) ** 2).mean()) * 32768)
    peak = float(np.abs(x).max() * 32768)
    try:
        hits = MP.find_bootstraps(MP.resample_to(x, FS, BS.FS))
    except Exception:                                        # noqa: BLE001
        hits = []
    pr = float(hits[0].get("mean_peak_ratio", 0.0)) if hits else 0.0
    cfo = (float(hits[0].get("fine_cfo_hz") or 0.0) if hits else None)
    return dict(rf=rf, rms=rms, peak=peak, peak_ratio=pr,
                atsc3=bool(hits), n_hits=len(hits), cfo=cfo)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--ant", default="Antenna A")
    ap.add_argument("--rfgain", type=int, default=4)
    ap.add_argument("--ifgr", type=int, default=32)
    ap.add_argument("--secs", type=float, default=1.2)
    ap.add_argument("--rounds", type=int, default=3)
    ap.add_argument("--chans", default="")
    ap.add_argument("--json", default=None)
    a = ap.parse_args(argv)

    chans = ([int(c) for c in a.chans.split(",")] if a.chans else CHANS)

    sdr = SoapySDR.Device("driver=sdrplay")
    try:
        avail = list(sdr.listAntennas(SOAPY_SDR_RX, 0))
    except Exception:                                        # noqa: BLE001
        avail = []
    log(f"antennas available: {avail}")
    if avail and a.ant not in avail:
        log(f"REFUSING: '{a.ant}' is not one of {avail} -- a sweep on the "
            f"wrong port is worse than no sweep")
        return 2
    sdr.setSampleRate(SOAPY_SDR_RX, 0, FS)
    sdr.setAntenna(SOAPY_SDR_RX, 0, a.ant)
    try:
        sdr.setGainMode(SOAPY_SDR_RX, 0, False)
    except Exception:                                        # noqa: BLE001
        pass
    sdr.setGain(SOAPY_SDR_RX, 0, "IFGR", a.ifgr)
    sdr.writeSetting("rfgain_sel", str(a.rfgain))
    time.sleep(0.3)
    # readback, never intention (the standing law)
    rb = dict(antenna=sdr.getAntenna(SOAPY_SDR_RX, 0),
              rfgain_sel=str(sdr.readSetting("rfgain_sel")),
              ifgr=float(sdr.getGain(SOAPY_SDR_RX, 0, "IFGR")),
              biasT=str(sdr.readSetting("biasT_ctrl")),
              rate=float(sdr.getSampleRate(SOAPY_SDR_RX, 0)))
    log(f"readback: {rb}")
    if rb["antenna"] != a.ant:
        log(f"REFUSING: antenna readback {rb['antenna']!r} != requested "
            f"{a.ant!r} -- an unverified port is a hope, not a measurement")
        return 2

    st = sdr.setupStream(SOAPY_SDR_RX, SOAPY_SDR_CS16)
    sdr.activateStream(st)
    buf = np.empty(2 * 65536, np.int16)
    rounds = []
    try:
        for r in range(a.rounds):
            log(f"round {r + 1}/{a.rounds} on {rb['antenna']}")
            out = {}
            for rf in chans:
                m = measure(sdr, st, buf, rf, a.secs)
                out[rf] = m
                if m["atsc3"]:
                    log(f"    RF{rf:<3} rms {m['rms']:8.0f}  peak_ratio "
                        f"{m['peak_ratio']:6.1f}  *** ATSC 3.0 ***")
            rounds.append(out)
    finally:
        sdr.deactivateStream(st)
        sdr.closeStream(st)
        del sdr

    # median across rounds; spread is reported so the reader can judge it
    summary = {}
    for rf in chans:
        rms = [rounds[i][rf]["rms"] for i in range(len(rounds))]
        prs = [rounds[i][rf]["peak_ratio"] for i in range(len(rounds))]
        hit = sum(1 for i in range(len(rounds)) if rounds[i][rf]["atsc3"])
        summary[rf] = dict(
            rms_med=float(np.median(rms)),
            rms_spread=float(max(rms) - min(rms)),
            pr_med=float(np.median(prs)),
            pr_spread=float(max(prs) - min(prs)),
            atsc3_rounds=f"{hit}/{len(rounds)}",
            atsc3=hit > len(rounds) // 2)

    print(f"\n  {'RF':>4} {'rms(med)':>10} {'+/-':>7} {'peak_ratio':>11} "
          f"{'+/-':>6}  ATSC3")
    for rf in chans:
        s = summary[rf]
        tag = ("ATSC3 %s" % s["atsc3_rounds"]) if s["atsc3"] else (
            "occupied" if s["rms_med"] > 500 else "-")
        print(f"  {rf:>4} {s['rms_med']:10.0f} {s['rms_spread']:7.0f} "
              f"{s['pr_med']:11.1f} {s['pr_spread']:6.1f}  {tag}")
    found = [rf for rf in chans if summary[rf]["atsc3"]]
    print(f"\n  ATSC 3.0 carriers on {rb['antenna']}: "
          f"{found if found else 'NONE'}")

    if a.json:
        json.dump(dict(readback=rb, rounds=rounds, summary=summary,
                       chans=chans, when=time.strftime("%Y-%m-%d %H:%M:%S")),
                  open(a.json, "w"), indent=1)
        print(f"  wrote {a.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
