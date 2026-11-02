#!/usr/bin/env python3
"""e64_antenna_probe.py -- is the ANTENNA PATH broken, or is RF33 just off?

E64 (8/10). RF33 has sat at 0.0% FEC for 6.5 h while the bootstrap keeps
arriving at near-normal strength (peak-ratio median 90.3 dead vs 95.7
delivering). That is a QUALITY collapse, not a level collapse -- the
classic signature of a missing low-noise amplifier: the wanted signal
still reaches the antenna, but without the LNA the receiver's own noise
figure dominates and the payload drops below threshold while the
bootstrap (designed for detection near -18 dB) still correlates.

There is a specific recent physical event: the SDR's USB was moved to the
Ubuntu box and back, and the yagi's amp is powered by BIAS-T on Antenna B
ONLY. The atsc3 code NEVER writes biasT_ctrl (verified by grep across the
repo) -- it inherits whatever the device already had. Sibling fleet tools
DO write it (a sibling spectrum tool --biast, Software-TV-Tuner --biast, and
radiotuna's DRM tools which set it explicitly FALSE), so the amp may have
been running on a setting that persisted in the device until the replug
reset it. 'SDR firmware sticks after restarts' cuts both ways.

THE MEASUREMENT, against a real baseline
    data/survey.json (8/08 02:42, Antenna B, IFGR 32, rfgain_sel 2) holds
    in_rms for RF14-36 -- UHF ONLY, no VHF-high, despite the file being
    named "survey". This probe repeats that measurement with
    IDENTICAL settings so the numbers are comparable, and adds 8-VSB
    carriers as the CONTROL: those transmitters are independent of RF33,
    so if they are down too the fault is on our side of the coax.

    Then an explicit bias-T A/B WITH READBACK (an unverified writeSetting
    is a hope, not a change), restoring the as-found state afterwards.

Radio discipline: single-tenant. The caller must stop the chain first;
this refuses to run if another owner holds the radio lock.

    python lab/e64_antenna_probe.py            # ~90 s of radio time
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)
_sib = os.environ.get("ATSC3_SIBLING_TOOLS")            # optional sibling toolbox
if _sib:
    sys.path.insert(0, _sib)

import SoapySDR                                                   # noqa: E402
from SoapySDR import SOAPY_SDR_RX, SOAPY_SDR_CS16                 # noqa: E402
import m2_pilots as MP                                            # noqa: E402
from atsc3 import bootstrap as bs                                 # noqa: E402

SoapySDR.SoapySDR_setLogLevel(SoapySDR.SOAPY_SDR_FATAL)

FS = 6.912e6
SECS = 1.2
# ATSC 3.0 carriers + 8-VSB CONTROLS (independent transmitters; strong in
# the 8/08 HDHomeRun scan: RF31 WETA / RF34 WRC / RF36 WTTG snq 100,
# RF21 MPT snq 89)
# TEN channels, not "the UHF band" -- a hand-picked list. Any null from
# this probe is a null about THESE TEN and nothing more (E88).
CHANS = [15, 21, 25, 29, 30, 31, 33, 34, 35, 36]
VSB_CONTROL = {15, 21, 31, 34, 35, 36}


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def center_hz(rf):
    lo = (174 + (rf - 7) * 6) if rf < 14 else (470 + (rf - 14) * 6)
    return (lo + 3.0) * 1e6


def lock_holder():
    try:
        import radio_lock                                        # noqa: F401
        d = json.load(open(radio_lock.LOCK))
        import psutil
        if psutil.pid_exists(d.get("pid", -1)):
            return d
    except Exception:                                            # noqa: BLE001
        pass
    return None


def measure(sdr, st, buf, rf):
    sdr.setFrequency(SOAPY_SDR_RX, 0, center_hz(rf))
    time.sleep(0.25)
    t0 = time.time()
    while time.time() - t0 < 0.2:                 # flush retune transient
        sdr.readStream(st, [buf], 65536, timeoutUs=300000)
    n_want = int(SECS * FS)
    iq = np.empty(n_want, np.complex64)
    got = 0
    deadline = time.time() + SECS * 2 + 3
    while got < n_want and time.time() < deadline:
        r = sdr.readStream(st, [buf], 65536, timeoutUs=500000)
        if r.ret > 0:
            n = min(r.ret, n_want - got)
            iq[got:got + n] = \
                buf[:2 * n].astype(np.float32).view(np.complex64) / 32768.0
            got += n
    if got < 1000:
        return dict(rf=rf, rms=0.0, peak_ratio=0.0, atsc3=False, n=0)
    in_rms = float(np.sqrt((np.abs(iq[:got]) ** 2).mean()) * 32768)
    hits = MP.find_bootstraps(MP.resample_to(iq[:got], FS, bs.FS))
    pr = hits[0].get("mean_peak_ratio", 0.0) if hits else 0.0
    return dict(rf=rf, rms=in_rms, peak_ratio=pr, atsc3=bool(hits),
                n=len(hits))


def sweep(sdr, st, buf, label):
    log(f"  sweeping {len(CHANS)} channels [{label}] ...")
    out = {}
    for rf in CHANS:
        m = measure(sdr, st, buf, rf)
        out[rf] = m
        tag = "ATSC3" if m["atsc3"] else ("8VSB-ctl" if rf in VSB_CONTROL
                                          else "-")
        log(f"    RF{rf:<3} rms {m['rms']:8.0f}  peak {m['peak_ratio']:6.1f}"
            f"  {tag}")
    return out


def bias_state(sdr):
    """-> the device's own answer, never our intention."""
    try:
        return str(sdr.readSetting("biasT_ctrl"))
    except Exception as e:                                       # noqa: BLE001
        return f"<unreadable: {type(e).__name__}>"


def set_bias(sdr, want):
    """want in ('true','false'). -> (accepted, readback). The literal must
    be 'true'/'false' -- labtuna_doctor lints 0/1/on/off as an ERROR."""
    try:
        sdr.writeSetting("biasT_ctrl", want)
    except Exception as e:                                       # noqa: BLE001
        return False, f"<write failed: {type(e).__name__}: {e}>"
    time.sleep(1.0)                       # LNA power-up settle
    rb = bias_state(sdr)
    return (rb.lower() == want), rb


def main():
    h = lock_holder()
    if h:
        log(f"REFUSING: radio held by {h.get('owner')} "
            f"(priority {h.get('priority')}, pid {h.get('pid')})")
        return 3

    base = {}
    try:
        b = json.load(open(os.path.join(ROOT, "data", "survey.json")))
        base = {r["rf"]: r for r in b["results"]}
        log(f"baseline: data/survey.json  {b['when']}  {b['ant']}")
    except Exception as e:                                       # noqa: BLE001
        log(f"no baseline available ({e})")

    sdr = SoapySDR.Device("driver=sdrplay")
    sdr.setSampleRate(SOAPY_SDR_RX, 0, FS)
    time.sleep(0.2)
    sdr.setAntenna(SOAPY_SDR_RX, 0, "Antenna B")
    try:
        sdr.setGainMode(SOAPY_SDR_RX, 0, False)
    except Exception:                                            # noqa: BLE001
        pass
    sdr.setGain(SOAPY_SDR_RX, 0, "IFGR", 32)          # == survey.py
    try:
        sdr.writeSetting("rfgain_sel", "2")
    except Exception:                                            # noqa: BLE001
        pass
    rb = {k: None for k in ("antenna", "rate", "freq")}
    rb["antenna"] = sdr.getAntenna(SOAPY_SDR_RX, 0)
    rb["rate"] = float(sdr.getSampleRate(SOAPY_SDR_RX, 0))
    log(f"radio readback: antenna={rb['antenna']} rate={rb['rate']/1e6:.3f}M "
        f"IFGR=32 rfgain_sel=2")

    keys = []
    try:
        keys = [k.key for k in sdr.getSettingInfo()]
    except Exception:                                            # noqa: BLE001
        pass
    log(f"device settings exposed: {keys}")
    as_found = bias_state(sdr)
    log(f"bias-T AS FOUND (device readback): {as_found!r}")

    st = sdr.setupStream(SOAPY_SDR_RX, SOAPY_SDR_CS16)
    sdr.activateStream(st)
    buf = np.empty(2 * 65536, np.int16)
    out = {"when": time.strftime("%Y-%m-%d %H:%M"), "ant": rb["antenna"],
           "bias_as_found": as_found, "settings": keys, "phases": {}}
    try:
        out["phases"]["as_found"] = sweep(sdr, st, buf, f"bias={as_found}")

        ok_on, rb_on = set_bias(sdr, "true")
        log(f"bias-T -> true : accepted={ok_on} readback={rb_on!r}")
        out["bias_on_readback"] = rb_on
        if ok_on:
            out["phases"]["bias_on"] = sweep(sdr, st, buf, "bias=true")
        else:
            log("  bias-T could not be enabled -- A/B not possible")

        want_back = as_found if as_found in ("true", "false") else "false"
        ok_r, rb_r = set_bias(sdr, want_back)
        log(f"bias-T restored -> {want_back} : accepted={ok_r} "
            f"readback={rb_r!r}")
        out["bias_restored_readback"] = rb_r
    finally:
        try:
            sdr.deactivateStream(st)
            sdr.closeStream(st)
        except Exception:                                        # noqa: BLE001
            pass
        del sdr

    # ---- the comparison that decides it
    print()
    log("RF   |  8/08 rms  |   now rms  | ratio |  8/08 pk |  now pk | kind")
    a = out["phases"]["as_found"]
    for rf in CHANS:
        b0 = base.get(rf, {})
        r0, r1 = b0.get("rms", 0.0), a[rf]["rms"]
        ratio = (r1 / r0) if r0 else float("nan")
        kind = "8VSB-ctl" if rf in VSB_CONTROL else \
               ("ATSC3" if b0.get("atsc3") or a[rf]["atsc3"] else "-")
        log(f"RF{rf:<3} | {r0:10.0f} | {r1:10.0f} | {ratio:5.2f} | "
            f"{b0.get('peak_ratio', 0.0):8.1f} | {a[rf]['peak_ratio']:7.1f} | "
            f"{kind}")
    ctl = [(a[rf]["rms"], base.get(rf, {}).get("rms", 0.0))
           for rf in CHANS if rf in VSB_CONTROL
           and base.get(rf, {}).get("rms", 0.0) > 0]
    if ctl:
        med = float(np.median([n / o for n, o in ctl]))
        log(f"8-VSB CONTROL median rms ratio now/8-08 = {med:.2f}  "
            f"-> {'ANTENNA PATH degraded' if med < 0.5 else 'antenna path OK'}")
        out["control_ratio"] = med
    if "bias_on" in out["phases"]:
        print()
        log("bias-T A/B (same channels, same gains):")
        for rf in CHANS:
            o, n = a[rf]["rms"], out["phases"]["bias_on"][rf]["rms"]
            log(f"  RF{rf:<3} rms {o:8.0f} -> {n:8.0f}  "
                f"({(n/o if o else float('nan')):.2f}x)   peak "
                f"{a[rf]['peak_ratio']:6.1f} -> "
                f"{out['phases']['bias_on'][rf]['peak_ratio']:6.1f}")

    p = os.path.join(HERE, "e64_antenna_probe.json")
    json.dump(out, open(p, "w"), indent=1)
    log(f"wrote {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
