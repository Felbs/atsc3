#!/usr/bin/env python3
"""band_sweep.py - is there ANY ATSC 3.0 (NextGen TV) in this market?

Tonight's probe answered only RF27/34/36 (all ATSC 1.0). This sweeps
RF7-13 (VHF-high) AND RF14-36 (UHF) for the bootstrap, which is the one
signal every 3.0 station must transmit and which cannot be confused with
8-VSB.

RANGE, STATED WITH THE RESULT (E88): this sweep is RF7-36. It does NOT
cover VHF-low (RF2-6) or RF37+. Any "there is nothing else in this
market" claim sourced from here is bounded by that list and must say so
-- a survey's null is only as wide as its channel list.

The variable below is named UHF but has held RF7-13 since the first
commit. The NAME is the lie, not the list, and that lie cost real
findings: see E88, where RF8's detection sat in this script's own output
file for six days while the market grid was published from a different
instrument.

Two passes, so we don't waste minutes capturing dead air:
  1. POWER SCAN - 0.5 s per channel, note which carry any signal at all.
  2. BOOTSTRAP HUNT - 8 s per occupied channel (>= 6 s guarantees a
     bootstrap even at the worst legal frame spacing of 5300 ms), then run
     the gated detector from ../bootstrap.py.

Radio discipline: takes the fleet lock ONCE for the whole sweep (owner
atsc3_sweep, prio 50, polite wait, heartbeat on a timer inside the read
loop), releases in finally. Captures are deleted after analysis unless a
bootstrap is found - a positive keeps its evidence.
"""
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import SoapySDR
from SoapySDR import SOAPY_SDR_CS16, SOAPY_SDR_RX

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
_sib = os.environ.get("ATSC3_SIBLING_TOOLS")            # optional sibling toolbox
if _sib:
    sys.path.insert(0, _sib)
# Reuse the GATED archive scanner unchanged - same code path that passed
# 32/32 selftest and produced tonight's RF27/34/36 verdicts. Never reach
# into detector internals from here; a second call path is a second bug.
import sweep_archive                               # noqa: E402
try:
    import radio_lock
except Exception:
    radio_lock = None

SoapySDR.SoapySDR_setLogLevel(SoapySDR.SOAPY_SDR_FATAL)

OWNER, PRI = "atsc3_sweep", 50
RATE = 8e6
ANT = "Antenna B"                 # Old Faithful: C is UHF-deaf
UHF = list(range(7, 14)) + list(range(14, 37))   # VHF-hi RF7-13 + UHF RF14-36
#   VHF-hi included 8/05: the user's HDHomeRun confirms ~6 unencrypted 3.0
#   services locally (virtual 107.1/132.1/145.100/154.1/158.1/158.8), and this
#   lab's known-active list includes RF7 and RF9 - do not assume UHF-only.
SCAN_S, HUNT_S = 0.5, 8.0
OUT = HERE / "band_sweep.json"


def center_hz(rf):
    lo = (174 + (rf - 7) * 6) if rf < 14 else (470 + (rf - 14) * 6)
    return (lo + 3.0) * 1e6


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def grab(sdr, st, secs, hb):
    n_want = int(secs * RATE)
    buf = np.empty(2 * 65536, np.int16)
    out = np.empty(2 * n_want, np.int16)
    got = 0
    t0 = time.time()
    while got < n_want and time.time() - t0 < secs * 3 + 5:
        r = sdr.readStream(st, [buf], 65536, timeoutUs=500000)
        if radio_lock is not None and time.time() - hb[0] >= 5.0:
            hb[0] = time.time()          # timer, never per-read (8/03 law)
            radio_lock.heartbeat()
        if r.ret > 0:
            n = min(r.ret, n_want - got)
            out[2 * got:2 * (got + n)] = buf[:2 * n]
            got += n
    return out[:2 * got], got / max(n_want, 1)


def main():
    if radio_lock is not None:
        if not radio_lock.acquire(OWNER, "ATSC 3.0 band sweep", PRI, wait_s=120):
            h = radio_lock.status() or {}
            log(f"radio busy: {h.get('owner')} (prio {h.get('priority')}) - "
                f"not seizing. Run again later.")
            return 2
    try:
        return sweep()
    finally:
        if radio_lock is not None:
            radio_lock.release(OWNER)


def sweep():
    hb = [time.time()]
    sdr = SoapySDR.Device("driver=sdrplay")
    sdr.setSampleRate(SOAPY_SDR_RX, 0, RATE)
    sdr.setAntenna(SOAPY_SDR_RX, 0, ANT)
    try:
        sdr.setGainMode(SOAPY_SDR_RX, 0, False)
    except Exception:
        pass
    sdr.setGain(SOAPY_SDR_RX, 0, "IFGR", 32)
    try:
        sdr.writeSetting("rfgain_sel", "2")
    except Exception:
        pass
    st = sdr.setupStream(SOAPY_SDR_RX, SOAPY_SDR_CS16)
    sdr.activateStream(st)
    results = {}
    try:
        log(f"PASS 1: power scan RF{UHF[0]}-RF{UHF[-1]} ({SCAN_S}s each)")
        occupied = []
        for rf in UHF:
            sdr.setFrequency(SOAPY_SDR_RX, 0, center_hz(rf))
            time.sleep(0.15)
            raw, frac = grab(sdr, st, SCAN_S, hb)
            iq = (raw[0::2].astype(np.float32) +
                  1j * raw[1::2].astype(np.float32)) / 32768.0
            sp = np.abs(np.fft.fftshift(np.fft.fft(iq[:1 << 16] *
                        np.hanning(1 << 16))))
            db = 20 * np.log10(sp + 1e-12)
            peak, med = float(db.max()), float(np.median(db))
            snr = peak - med
            results[rf] = {"scan_snr_db": round(snr, 1),
                           "scan_frac": round(frac, 3)}
            if snr >= 15.0:
                occupied.append(rf)
            log(f"  RF{rf:<3} peak-median {snr:5.1f} dB"
                f"{'   <- occupied' if snr >= 15.0 else ''}")
        log(f"PASS 2: bootstrap hunt on {len(occupied)} occupied channels "
            f"({HUNT_S}s each): {occupied}")
        for rf in occupied:
            sdr.setFrequency(SOAPY_SDR_RX, 0, center_hz(rf))
            time.sleep(0.2)
            raw, frac = grab(sdr, st, HUNT_S, hb)
            tmp = HERE / f"_sweep_rf{rf}.cs16"
            raw.tofile(tmp)
            try:
                r = sweep_archive.scan_file(str(tmp), RATE)
            except Exception as e:
                r = {"verdict": f"SCAN_ERROR: {str(e)[:80]}"}
            results[rf].update({"hunt_frac": round(frac, 3), "bootstrap": r})
            verdict = str(r.get("verdict", r) if isinstance(r, dict) else r)
            log(f"  RF{rf:<3} -> {verdict}")
            if "NO_BOOT" in verdict.upper() or "ERROR" in verdict.upper():
                tmp.unlink(missing_ok=True)   # negatives don't hoard disk
            else:                             # a positive KEEPS its evidence
                tmp.rename(HERE / f"hit_rf{rf}.cs16")
                log(f"     *** ATSC 3.0 CANDIDATE - evidence kept: hit_rf{rf}.cs16")
    finally:
        try:
            sdr.deactivateStream(st)
            sdr.closeStream(st)
        except Exception:
            pass
    OUT.write_text(json.dumps(results, indent=1, default=str), encoding="utf-8")
    hits = [rf for rf, v in results.items()
            if "bootstrap" in v and
            "NO_BOOT" not in str(v["bootstrap"]).upper()]
    log(f"SWEEP DONE -> {OUT.name}. ATSC 3.0 channels found: "
        f"{hits if hits else 'NONE (this market is still all ATSC 1.0)'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
