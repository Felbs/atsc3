#!/usr/bin/env python3
"""m11_rateprobe.py -- can the SDR give us 6.912 Msps directly?

M9's profile said `scipy.signal.resample_poly` is 70.2 ms of a 243.4 ms Frame
budget -- 29% of the whole pipeline -- and that it exists for exactly one
reason: `m7_capture.py` asks for 8 Msps and the decoder works at 6.912.  If the
radio will simply deliver 6.912 Msps the stage does not get optimised, it
disappears.

This asks, and then it READS BACK, because a `setSampleRate` that is not
read back is a hope (the unverified-control-writes law).  It reports the
readback verbatim -- including the case where the driver quietly rounds to
something else, in which case the answer is "no" and the resampler stays.

It then takes a short capture at whatever rate came back, so the claim
"6.912 Msps decodes" is a decode and not an inference.

Radio discipline: owner atsc3_watch, priority 60, polite wait, release in
finally, heartbeat on a timer, never inside the meteor window.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import SoapySDR
from SoapySDR import SOAPY_SDR_CS16, SOAPY_SDR_RX

HERE = Path(__file__).resolve().parent
_sib = os.environ.get("ATSC3_SIBLING_TOOLS")            # optional sibling toolbox
if _sib:
    sys.path.insert(0, _sib)
try:
    import radio_lock
except Exception:                                              # noqa: BLE001
    radio_lock = None

SoapySDR.SoapySDR_setLogLevel(SoapySDR.SOAPY_SDR_FATAL)

OWNER, PRI = "atsc3_watch", 60
WANT = 6.912e6
ANT = "Antenna B"


def center_hz(rf):
    lo = (174 + (rf - 7) * 6) if rf < 14 else (470 + (rf - 14) * 6)
    return (lo + 3.0) * 1e6


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def in_meteor_window(now=None):
    now = now or dt.datetime.now()
    mins = now.hour * 60 + now.minute
    return 2 * 60 + 10 <= mins <= 5 * 60 + 40


def probe(rf, want, secs, out_path, ant=ANT):
    sdr = SoapySDR.Device("driver=sdrplay")
    rep = {"requested_rate": want}

    try:
        rep["rate_range"] = [f"{r.minimum():.0f}..{r.maximum():.0f}"
                             for r in sdr.getSampleRateRange(SOAPY_SDR_RX, 0)]
    except Exception as e:                                     # noqa: BLE001
        rep["rate_range"] = f"<{e}>"
    try:
        rep["rate_list"] = list(sdr.listSampleRates(SOAPY_SDR_RX, 0))
    except Exception:                                          # noqa: BLE001
        rep["rate_list"] = None
    try:
        rep["bw_list"] = list(sdr.listBandwidths(SOAPY_SDR_RX, 0))
    except Exception:                                          # noqa: BLE001
        rep["bw_list"] = None

    sdr.setSampleRate(SOAPY_SDR_RX, 0, want)
    time.sleep(0.2)
    rep["readback_rate"] = float(sdr.getSampleRate(SOAPY_SDR_RX, 0))
    rep["rate_exact"] = abs(rep["readback_rate"] - want) < 1.0

    sdr.setAntenna(SOAPY_SDR_RX, 0, ant)
    try:
        sdr.setGainMode(SOAPY_SDR_RX, 0, False)
    except Exception:                                          # noqa: BLE001
        pass
    sdr.setGain(SOAPY_SDR_RX, 0, "IFGR", 32)
    try:
        sdr.writeSetting("rfgain_sel", "2")
    except Exception:                                          # noqa: BLE001
        pass
    sdr.setFrequency(SOAPY_SDR_RX, 0, center_hz(rf))
    time.sleep(0.3)
    for k in ("rfgain_sel",):
        try:
            rep[k] = sdr.readSetting(k)
        except Exception:                                      # noqa: BLE001
            rep[k] = "<unreadable>"
    rep["antenna"] = sdr.getAntenna(SOAPY_SDR_RX, 0)
    rep["freq"] = float(sdr.getFrequency(SOAPY_SDR_RX, 0))
    try:
        rep["bandwidth"] = float(sdr.getBandwidth(SOAPY_SDR_RX, 0))
    except Exception:                                          # noqa: BLE001
        rep["bandwidth"] = None
    log(f"    readback: {json.dumps(rep, default=str)}")

    fs = rep["readback_rate"]
    st = sdr.setupStream(SOAPY_SDR_RX, SOAPY_SDR_CS16)
    sdr.activateStream(st)
    n_want = int(secs * fs)
    buf = np.empty(2 * 65536, np.int16)
    got, hb, overflows = 0, [time.time()], 0
    try:
        with open(out_path, "wb") as fh:
            t_warm = time.time()
            while time.time() - t_warm < 0.25:
                sdr.readStream(st, [buf], 65536, timeoutUs=500000)
            t0 = time.time()
            deadline = t0 + secs * 1.5 + 5
            while got < n_want and time.time() < deadline:
                r = sdr.readStream(st, [buf], 65536, timeoutUs=500000)
                now = time.time()
                if radio_lock is not None and now - hb[0] >= 5.0:
                    hb[0] = now                 # TIMER, not per-read
                    radio_lock.heartbeat()
                if r.ret > 0:
                    n = min(r.ret, n_want - got)
                    fh.write(buf[:2 * n].tobytes())
                    got += n
                elif r.ret < 0:
                    overflows += 1
            t1 = time.time()
    finally:
        try:
            sdr.deactivateStream(st)
            sdr.closeStream(st)
        except Exception:                                      # noqa: BLE001
            pass
    wall = t1 - t0
    rep["integrity"] = dict(samples=got, wall_s=round(wall, 3), rate=fs,
                            expected=int(wall * fs),
                            ratio=round(got / (wall * fs), 6),
                            overflow_reads=overflows,
                            got_s=round(got / fs, 3))
    rep["integrity"]["integrity_ok"] = bool(
        0.985 <= rep["integrity"]["ratio"] <= 1.02)
    return rep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rf", type=int, default=33)
    ap.add_argument("--rate", type=float, default=WANT)
    ap.add_argument("--secs", type=float, default=8.0)
    ap.add_argument("--ant", default=ANT)
    ap.add_argument("--out", default=None)
    ap.add_argument("--wait", type=float, default=120.0)
    a = ap.parse_args()

    if in_meteor_window():
        log("REFUSING: inside the 02:10-05:40 meteor window.")
        return 3
    out = Path(a.out) if a.out else HERE / f"rate{int(a.rate/1000)}_rf{a.rf}.cs16"

    if radio_lock is not None:
        if not radio_lock.acquire(OWNER, f"ATSC 3.0 rate probe RF{a.rf}",
                                  PRI, wait_s=a.wait):
            h = radio_lock.status() or {}
            log(f"radio held by {h.get('owner')} (prio {h.get('priority')}); "
                f"NOT seizing.")
            return 2
    try:
        log(f"probing RF{a.rf} @ {center_hz(a.rf)/1e6:.3f} MHz for "
            f"{a.rate/1e6:.4f} Msps")
        rep = probe(a.rf, a.rate, a.secs, out, ant=a.ant)
    finally:
        if radio_lock is not None:
            radio_lock.release(OWNER)

    rep.update(rf=a.rf, center_hz=center_hz(a.rf), out=str(out),
               utc=dt.datetime.now(dt.timezone.utc).isoformat())
    js = Path(str(out.with_suffix("")) + ".json")
    js.write_text(json.dumps(rep, indent=1, default=str), encoding="utf-8")
    print("\n  === RATE READBACK ===")
    print(f"    requested   {a.rate:.1f} Hz")
    print(f"    readback    {rep['readback_rate']:.1f} Hz    "
          f"{'EXACT' if rep['rate_exact'] else '*** DIFFERENT ***'}")
    print(f"    bandwidth   {rep['bandwidth']}")
    i = rep["integrity"]
    print(f"    capture     {i['samples']} samples / {i['wall_s']} s wall "
          f"(ratio {i['ratio']:.4f}, {i['overflow_reads']} overflow reads)  "
          f"{'PASS' if i['integrity_ok'] else '*** VOID ***'}")
    print(f"    -> {out}  and {js.name}")
    return 0 if rep["rate_exact"] and i["integrity_ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
