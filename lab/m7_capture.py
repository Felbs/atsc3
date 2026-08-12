#!/usr/bin/env python3
"""m7_capture.py -- ONE long RF33 capture, for ROUTE object completeness.

Why a longer capture is justified NOW and was not before: M6 closed the whole
chain (74/74 FEC Blocks, ALP -> IPv4/UDP -> ROUTE).  The only thing standing
between us and a playable file is whether a DASH segment finished inside the
recording window, which is purely a function of record LENGTH.  M4/M5's
"we need more margin" was refuted; this is not that request.

Radio discipline (fleet law):
  * radio_lock owner "atsc3", priority 50, polite wait, bounded retries
  * heartbeat on a TIMER inside the read loop, never per-read (8/03 law)
  * release in finally
  * capture-integrity gate: got_samples == wall_seconds * fs, or VOID
  * never between 02:10 and 05:40 local (meteor window)

Settings mirror band_sweep.py exactly, because that is what produced the
capture that decoded: Antenna B (Old Faithful yagi), 8 Msps, IFGR 32,
rfgain_sel 2.
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

OWNER, PRI = "atsc3", 50
RATE = 8e6
ANT = "Antenna B"      # default; --ant overrides (C = discone, the VHF-hi port)
MB = 1 << 20


def center_hz(rf):
    lo = (174 + (rf - 7) * 6) if rf < 14 else (470 + (rf - 14) * 6)
    return (lo + 3.0) * 1e6


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def in_meteor_window(now=None):
    now = now or dt.datetime.now()
    mins = now.hour * 60 + now.minute
    return 2 * 60 + 10 <= mins <= 5 * 60 + 40


def capture(rf, secs, out_path, ant=ANT):
    """Stream to disk in chunks.  Returns the integrity dict."""
    sdr = SoapySDR.Device("driver=sdrplay")
    sdr.setSampleRate(SOAPY_SDR_RX, 0, RATE)
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
    # readback or it did not happen (unverified-control-writes law)
    read_back = {}
    for k in ("rfgain_sel",):
        try:
            read_back[k] = sdr.readSetting(k)
        except Exception:                                      # noqa: BLE001
            read_back[k] = "<unreadable>"
    read_back["antenna"] = sdr.getAntenna(SOAPY_SDR_RX, 0)
    read_back["rate"] = sdr.getSampleRate(SOAPY_SDR_RX, 0)
    read_back["freq"] = None

    sdr.setFrequency(SOAPY_SDR_RX, 0, center_hz(rf))
    time.sleep(0.3)
    read_back["freq"] = sdr.getFrequency(SOAPY_SDR_RX, 0)
    log(f"    readback: {read_back}")

    st = sdr.setupStream(SOAPY_SDR_RX, SOAPY_SDR_CS16)
    sdr.activateStream(st)
    n_want = int(secs * RATE)
    buf = np.empty(2 * 65536, np.int16)
    got, hb, t_wall0 = 0, [time.time()], None
    overflows = 0
    try:
        with open(out_path, "wb") as fh:
            # discard the first ~0.2 s: the stream settles, and a sacrificial
            # first read is the 8/01 meteor lesson.
            t_warm = time.time()
            while time.time() - t_warm < 0.25:
                sdr.readStream(st, [buf], 65536, timeoutUs=500000)
            t_wall0 = time.time()
            deadline = t_wall0 + secs * 1.5 + 5
            while got < n_want and time.time() < deadline:
                r = sdr.readStream(st, [buf], 65536, timeoutUs=500000)
                now = time.time()
                if radio_lock is not None and now - hb[0] >= 5.0:
                    hb[0] = now              # TIMER, not per-read (8/03 law)
                    radio_lock.heartbeat()
                if r.ret > 0:
                    n = min(r.ret, n_want - got)
                    fh.write(buf[:2 * n].tobytes())
                    got += n
                    if got % (RATE * 10) < 65536:
                        log(f"      {got/RATE:6.1f} s / {secs:.0f} s "
                            f"({got*4/MB:.0f} MB)")
                elif r.ret < 0:
                    overflows += 1
            t_wall1 = time.time()
    finally:
        try:
            sdr.deactivateStream(st)
            sdr.closeStream(st)
        except Exception:                                      # noqa: BLE001
            pass
    wall = t_wall1 - t_wall0
    expect = wall * RATE
    integ = dict(samples=got, wall_s=round(wall, 3), rate=RATE,
                 expected=int(expect), ratio=round(got / expect, 6),
                 overflow_reads=overflows,
                 requested_s=secs, got_s=round(got / RATE, 3),
                 readback=read_back)
    # capture-integrity law: samples == wall*fs (within the tail of the last
    # partial read) or the capture is VOID.
    integ["integrity_ok"] = bool(0.985 <= integ["ratio"] <= 1.02)
    return integ


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rf", type=int, default=33)
    ap.add_argument("--secs", type=float, default=120.0)
    ap.add_argument("--out", default=None)
    ap.add_argument("--wait", type=float, default=180.0)
    ap.add_argument("--retries", type=int, default=3)
    ap.add_argument("--ant", default=ANT,
                    help="Antenna A/B/C. B is the UHF yagi and is VHF-deaf; "
                         "C is the discone, which is the reverse.")
    a = ap.parse_args()

    if in_meteor_window():
        log("REFUSING: inside the 02:10-05:40 meteor window.")
        return 3

    out = Path(a.out) if a.out else HERE / f"long_rf{a.rf}.cs16"
    if radio_lock is not None:
        held = False
        for attempt in range(a.retries):
            if radio_lock.acquire(OWNER, f"ATSC 3.0 RF{a.rf} long capture",
                                  PRI, wait_s=a.wait):
                held = True
                break
            h = radio_lock.status() or {}
            log(f"    radio busy: {h.get('owner')} (prio {h.get('priority')})"
                f" -- attempt {attempt+1}/{a.retries}, waiting")
            time.sleep(15)
        if not held:
            log("radio never freed; NOT seizing. Try later.")
            return 2
    try:
        log(f"capturing RF{a.rf} @ {center_hz(a.rf)/1e6:.3f} MHz, "
            f"{a.ant}, {RATE/1e6:.1f} Msps, {a.secs:.0f} s -> {out.name}")
        integ = capture(a.rf, a.secs, out, ant=a.ant)
    finally:
        if radio_lock is not None:
            radio_lock.release(OWNER)

    meta = dict(rf=a.rf, center_hz=center_hz(a.rf), rate=RATE, antenna=a.ant,
                secs=integ["got_s"], fmt="cs16",
                utc=dt.datetime.now(dt.timezone.utc).isoformat(),
                integrity=integ)
    Path(str(out.with_suffix("")) + ".json").write_text(
        json.dumps(meta, indent=1, default=str), encoding="utf-8")
    log(f"    samples {integ['samples']} vs wall*fs {integ['expected']} "
        f"(ratio {integ['ratio']:.4f}), overflow reads "
        f"{integ['overflow_reads']}")
    if not integ["integrity_ok"]:
        log("*** CAPTURE-INTEGRITY GATE FAILED -- capture is VOID ***")
        return 4
    log(f"    CAPTURE-INTEGRITY GATE PASS -- {out} "
        f"({out.stat().st_size/MB:.0f} MB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
