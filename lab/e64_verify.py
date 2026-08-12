#!/usr/bin/env python3
"""e64_verify.py -- did moving the USB cable actually help? (E64 verify)

MULTI-PASS BY CONSTRUCTION. The E64 correction exists because a single
sweep caught an interference PEAK and I generalised it. The ingress is
intermittent, so one number is not evidence -- this takes N passes spread
over several minutes in ONE radio session and reports the DISTRIBUTION
(per channel: median/min/max across passes, and how much it moved).

Baseline is data/survey.json (8/08 02:42, Antenna B, IFGR 32,
rfgain_sel 2) -- the same settings, so the ratios mean something.
"""
import argparse, json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, ROOT)
import SoapySDR
from SoapySDR import SOAPY_SDR_RX, SOAPY_SDR_CS16
SoapySDR.SoapySDR_setLogLevel(SoapySDR.SOAPY_SDR_FATAL)
import e64_antenna_probe as P

VACANT = [14, 16, 19, 20, 22, 23, 24, 32]     # noise floor
INGRESS = [16, 19, 20, 22, 23, 24]            # the 485-535 MHz cluster
SIGNAL = [33, 31]                             # target + 8-VSB control

ap = argparse.ArgumentParser()
ap.add_argument("--passes", type=int, default=6)
ap.add_argument("--spacing", type=float, default=45.0)
a = ap.parse_args()

h = P.lock_holder()
if h:
    print(f"REFUSING: radio held by {h.get('owner')} (pid {h.get('pid')})")
    sys.exit(3)
base = {r["rf"]: r for r in json.load(
    open(os.path.join(ROOT, "data", "survey.json")))["results"]}

sdr = SoapySDR.Device("driver=sdrplay")
sdr.setSampleRate(SOAPY_SDR_RX, 0, P.FS); time.sleep(0.2)
sdr.setAntenna(SOAPY_SDR_RX, 0, "Antenna B")
try: sdr.setGainMode(SOAPY_SDR_RX, 0, False)
except Exception: pass
sdr.setGain(SOAPY_SDR_RX, 0, "IFGR", 32)
try: sdr.writeSetting("rfgain_sel", "2")
except Exception: pass
print(f"radio: {sdr.getAntenna(SOAPY_SDR_RX,0)}  bias-T {P.bias_state(sdr)!r}"
      f"  (unchanged from the as-found state)", flush=True)
st = sdr.setupStream(SOAPY_SDR_RX, SOAPY_SDR_CS16); sdr.activateStream(st)
buf = np.empty(2 * 65536, np.int16)
chans = VACANT + SIGNAL
passes = []
try:
    for k in range(a.passes):
        t0 = time.time()
        row = {rf: P.measure(sdr, st, buf, rf) for rf in chans}
        passes.append(row)
        vr = [row[rf]["rms"] / base[rf]["rms"] for rf in VACANT
              if base.get(rf, {}).get("rms", 0) > 0]
        print(f"  pass {k+1}/{a.passes} @{time.strftime('%H:%M:%S')}  "
              f"floor median {np.median(vr):5.2f}x "
              f"({20*np.log10(np.median(vr)):+5.1f} dB)   "
              f"RF33 rms {row[33]['rms']:7.0f} peak {row[33]['peak_ratio']:5.1f}",
              flush=True)
        if k < a.passes - 1:
            time.sleep(max(0.0, a.spacing - (time.time() - t0)))
finally:
    try: sdr.deactivateStream(st); sdr.closeStream(st)
    except Exception: pass
    del sdr

print(f"\n=== NOISE FLOOR vs 8/08 baseline, {a.passes} passes over "
      f"{(a.passes-1)*a.spacing/60:.1f} min ===")
print(f"{'RF':>4} {'8/08':>7} {'median':>8} {'min':>7} {'max':>7} "
      f"{'med dB':>7} {'spread dB':>9}  zone")
allmed = []
for rf in VACANT:
    o = base.get(rf, {}).get("rms", 0.0)
    v = [p[rf]["rms"] for p in passes]
    r = [x / o for x in v] if o else [1]
    med = float(np.median(r)); allmed.append(med)
    spread = 20*np.log10(max(v)/max(min(v), 1e-9))
    print(f"RF{rf:<3} {o:7.0f} {np.median(v):8.0f} {min(v):7.0f} {max(v):7.0f} "
          f"{20*np.log10(med):+7.1f} {spread:9.1f}  "
          f"{'485-535 cluster' if rf in INGRESS else ''}")
m = float(np.median(allmed))
ing = float(np.median([allmed[VACANT.index(rf)] for rf in INGRESS]))
print(f"\n  FLOOR median over all vacant   : {m:.2f}x ({20*np.log10(m):+.1f} dB)")
print(f"  FLOOR median, 485-535 cluster  : {ing:.2f}x ({20*np.log10(ing):+.1f} dB)")
print(f"\n=== RF33 / control across passes ===")
for rf in SIGNAL:
    v = [p[rf]["rms"] for p in passes]; pk = [p[rf]["peak_ratio"] for p in passes]
    b0 = base.get(rf, {})
    print(f"  RF{rf:<3} 8/08 rms {b0.get('rms',0):7.0f} pk {b0.get('peak_ratio',0):5.1f}"
          f"  ->  now rms med {np.median(v):7.0f}  pk med {np.median(pk):5.1f} "
          f"(min {min(pk):5.1f} max {max(pk):5.1f})")
json.dump({"when": time.strftime("%Y-%m-%d %H:%M"), "passes": a.passes,
           "spacing": a.spacing,
           "data": [{str(rf): p[rf] for rf in chans} for p in passes]},
          open(os.path.join(HERE, "e64_verify.json"), "w"), indent=1, default=str)
print(f"\nFLOOR VERDICT: {'BACK NEAR BASELINE' if m < 1.4 else 'STILL ELEVATED'} "
      f"(median {m:.2f}x)")
