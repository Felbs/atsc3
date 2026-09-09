#!/usr/bin/env python3
"""atsc3_aim.py -- live antenna-aiming feedback on one carrier.

Aiming needs a number that responds in SECONDS and that you can trust.
The old aim loop ran a full decode probe per position (~75 s at FEC 0,
because LDPC iterates hardest on exactly the garbage signal you are trying
to escape) and then compared a SINGLE peak ratio per position. Both halves
were wrong:

  SLOW    a 75 s answer means you cannot sweep an antenna; you sample four
          positions and give up. This tool runs the bootstrap detector
          only -- no LDPC -- so a probe costs ~3 s and the SDR stays open
          across probes.

  NOISY   at Smith Point on 9/08, back-to-back 4 s probes at IDENTICAL
          settings ranged 29.2 to 48.3 on RF33. A single probe per antenna
          position is indistinguishable from luck, and it already fooled us
          once (a bias-T A/B that "showed" a large gain, and vanished when
          the conditions were interleaved and repeated). So this reports a
          ROLLING MEDIAN over a window of probes, plus the spread, and
          never asks you to act on one sample.

Decode is still the goal, not the metric: when the rolling median crosses
--launch-peak the tool drops the radio, spends one real decode probe, and
execs the full TV stack if FEC clears --launch-fec. Otherwise it goes back
to aiming, having lost ~75 s once rather than on every position.

Reference points measured on this fleet: peak ~80 still gave FEC 0%;
peak ~113 gave FEC 100%. Hence the default launch threshold of 75.

Usage:
    python tools/atsc3_aim.py                      # RF8, Antenna B
    python tools/atsc3_aim.py --rf 33 --ifgr 32
"""
from __future__ import annotations

import argparse
import collections
import os
import statistics
import subprocess
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

# radio_lock lives in the sibling gr-radiotuna toolbox; same discipline as
# atsc3_survey. Degrade to a no-op only if it genuinely cannot be found,
# and SAY SO -- a silent no-op lock is how the arbitration quietly dies.
for _cand in [os.environ.get("ATSC3_SIBLING_TOOLS"),
              os.path.join(os.path.dirname(ROOT), "gr-radiotuna", "tools"),
              os.path.join(os.path.expanduser("~"), "radiotuna", "tools")]:
    if _cand and os.path.isfile(os.path.join(_cand, "radio_lock.py")):
        sys.path.insert(0, _cand)
        break
try:
    import radio_lock                                             # noqa: E402
except Exception:                                                 # noqa: BLE001
    radio_lock = None

SoapySDR.SoapySDR_setLogLevel(SoapySDR.SOAPY_SDR_FATAL)
FS = 6.912e6
OWNER, PRI = "atsc3_aim", 50


def center_hz(rf):
    # US TV allocations are NOT one arithmetic run. Two discontinuities:
    # 72-76 MHz is aeronautical (so ch4 -> ch5 jumps), and high-VHF
    # restarts at 174. FIXED 9/08: the old one-liner treated everything
    # below 14 as high-VHF and returned 144 MHz for ch2 -- off by 87 MHz,
    # which made low-VHF silently unsurveyable rather than merely absent.
    if rf < 5:
        lo = 54 + (rf - 2) * 6          # ch2-4   54-72
    elif rf < 7:
        lo = 76 + (rf - 5) * 6          # ch5-6   76-88
    elif rf < 14:
        lo = 174 + (rf - 7) * 6         # ch7-13  174-216
    else:
        lo = 470 + (rf - 14) * 6        # ch14+   470-608
    return (lo + 3.0) * 1e6


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def bar(v, lo=0.0, hi=120.0, w=32):
    n = max(0, min(w, int(round(w * (v - lo) / (hi - lo)))))
    return "#" * n + "." * (w - n)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rf", type=int, default=8,
                    help="RF channel. Default 8: at Smith Point RF8 "
                         "(VHF-high) beat both UHF carriers roughly 2:1.")
    ap.add_argument("--ant", default="Antenna B")
    ap.add_argument("--rfgain", type=int, default=2)
    ap.add_argument("--ifgr", type=int, default=42,
                    help="Per-channel by necessity: one global IFGR put RF8 "
                         "at -2.2 dBFS and RF25 at -26.5 dBFS. RF8 42, "
                         "RF33 32, RF25 22.")
    ap.add_argument("--secs", type=float, default=2.5, help="per probe")
    ap.add_argument("--window", type=int, default=5,
                    help="probes in the rolling median")
    ap.add_argument("--launch-peak", type=float, default=75.0,
                    help="rolling median that earns one real decode probe")
    ap.add_argument("--launch-fec", type=float, default=90.0)
    ap.add_argument("--live-dir", default=os.path.join(ROOT, "data", "tv_aim"))
    ap.add_argument("--python", default=os.path.join(ROOT, ".venv", "bin",
                                                     "python"),
                    help="interpreter for the TV stack. atsc3_watch_av's "
                         "find_python only probes radioconda paths and dies "
                         "on the Pi without this.")
    ap.add_argument("--no-launch", action="store_true",
                    help="pure meter: never spend a decode probe")
    a = ap.parse_args()

    if radio_lock is not None:
        if not radio_lock.acquire(OWNER, f"RF{a.rf} antenna aiming", PRI,
                                  wait_s=120):
            h = radio_lock.status() or {}
            log(f"radio busy: {h.get('owner')} (prio {h.get('priority')})")
            return 2
    else:
        log("WARNING: radio_lock not found -- opening the SDR WITHOUT "
            "arbitration (set ATSC3_SIBLING_TOOLS).")

    sdr = SoapySDR.Device("driver=sdrplay")
    sdr.setSampleRate(SOAPY_SDR_RX, 0, FS)
    time.sleep(0.2)
    fs = float(sdr.getSampleRate(SOAPY_SDR_RX, 0))
    sdr.setAntenna(SOAPY_SDR_RX, 0, a.ant)
    try:
        sdr.setGainMode(SOAPY_SDR_RX, 0, False)
    except Exception:                                          # noqa: BLE001
        pass
    sdr.setGain(SOAPY_SDR_RX, 0, "IFGR", a.ifgr)
    try:
        sdr.writeSetting("rfgain_sel", str(a.rfgain))
    except Exception:                                          # noqa: BLE001
        pass
    sdr.setFrequency(SOAPY_SDR_RX, 0, center_hz(a.rf))
    st = sdr.setupStream(SOAPY_SDR_RX, SOAPY_SDR_CS16)
    sdr.activateStream(st)
    buf = np.empty(2 * 65536, np.int16)

    owned = True          # do we still hold the SDR + lock?
    hist = collections.deque(maxlen=a.window)
    best = 0.0
    hb = [time.time()]
    log(f"aiming RF{a.rf} ({center_hz(a.rf)/1e6:.3f} MHz) on {a.ant}, "
        f"rfgain_sel={a.rfgain} IFGR={a.ifgr}")
    log(f"move the antenna slowly; judge by the MEDIAN column, not one probe "
        f"(window={a.window})")
    try:
        while True:
            if radio_lock is not None and radio_lock.should_yield():
                log("yielding the radio to a higher-priority waiter.")
                return 0
            t0 = time.time()
            while time.time() - t0 < 0.15:
                sdr.readStream(st, [buf], 65536, timeoutUs=300000)
            n_want = int(a.secs * fs)
            iq = np.empty(n_want, np.complex64)
            got = 0
            deadline = time.time() + a.secs * 2 + 3
            while got < n_want and time.time() < deadline:
                r = sdr.readStream(st, [buf], 65536, timeoutUs=500000)
                if radio_lock is not None and time.time() - hb[0] >= 5.0:
                    hb[0] = time.time()
                    radio_lock.heartbeat()
                if r.ret > 0:
                    n = min(r.ret, n_want - got)
                    x = buf[:2 * n].astype(np.float32).view(np.complex64)
                    iq[got:got + n] = x / 32768.0
                    got += n
            y = MP.resample_to(iq[:got], fs, bs.FS)
            hits = MP.find_bootstraps(y)
            pr = hits[0].get("mean_peak_ratio", 0.0) if hits else 0.0
            rms = float(np.sqrt((np.abs(iq[:got]) ** 2).mean()) * 32768)
            hist.append(pr)
            med = statistics.median(hist)
            spread = (max(hist) - min(hist)) if len(hist) > 1 else 0.0
            flag = ""
            if med > best and len(hist) == hist.maxlen:
                best, flag = med, "  <-- BEST"
            log(f"  peak {pr:6.1f} | median {med:6.1f} (n={len(hist)}, "
                f"spread {spread:5.1f}) rms {rms:6.0f} |{bar(med)}|{flag}")

            if (not a.no_launch and len(hist) == hist.maxlen
                    and med >= a.launch_peak):
                log(f"median {med:.1f} >= {a.launch_peak} -- spending one "
                    f"real decode probe")
                sdr.deactivateStream(st)
                sdr.closeStream(st)
                del sdr
                if radio_lock is not None:
                    radio_lock.release(OWNER)
                owned = False
                cmd = [a.python, "-m", "atsc3", "watch", "--rf", str(a.rf),
                       "--ant", a.ant, "--secs", "20", "--player", "none",
                       "--live-dir", os.path.join(ROOT, "data", "aim_probe"),
                       "--accel", "cpu", "--assets", "none",
                       "--decode-procs", "2", "--threads", "2"]
                out = subprocess.run(cmd, capture_output=True, text=True,
                                     cwd=ROOT).stdout
                fec = 0.0
                for tok in out.split("["):
                    if "% now]" in tok:
                        try:
                            fec = max(fec, float(tok.split("%")[0].strip()))
                        except ValueError:
                            pass
                log(f"decode probe: best interval FEC {fec:.1f}%")
                if fec >= a.launch_fec:
                    log("SIGNAL GOOD -- launching the TV stack")
                    os.execv(a.python, [
                        a.python, os.path.join(HERE, "atsc3_watch_av.py"),
                        "--rf", str(a.rf), "--ant", a.ant, "--cc",
                        "--python", a.python, "--live-dir", a.live_dir])
                log("not decodable yet -- resuming aiming")
                return 1
    except KeyboardInterrupt:
        log(f"stopped. best rolling median this run: {best:.1f}")
    finally:
        # the launch path already dropped both; releasing twice is not safe
        # to assume idempotent, so gate on ownership rather than on luck.
        if owned:
            try:
                sdr.deactivateStream(st)
                sdr.closeStream(st)
            except Exception:                                  # noqa: BLE001
                pass
            if radio_lock is not None:
                radio_lock.release(OWNER)
    return 0


if __name__ == "__main__":
    sys.exit(main())
