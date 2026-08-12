#!/usr/bin/env python3
"""e64_remeasure.py -- the E64 floor measurement, repeated.

WHY (the law: a hit is a POINTER, never the evidence): the 08:49 probe
measured +4.7 dB of noise floor -- but the chain's notch telemetry shows
the interference is INTERMITTENT and PEAKED at 08:20-08:50, exactly the
window the probe ran in (out-of-band/in-band median 0.43 then, 0.075 when
delivering, 0.11-0.16 at 05:40-06:40, 0.13-0.19 now). One measurement
taken at a peak, generalised to a diagnosis, is how instruments lie.
This repeats it against BOTH baselines: 8/08 (healthy) and 08:49 (peak).
"""
import json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, ROOT)
import SoapySDR
from SoapySDR import SOAPY_SDR_RX, SOAPY_SDR_CS16
SoapySDR.SoapySDR_setLogLevel(SoapySDR.SOAPY_SDR_FATAL)
import e64_antenna_probe as P

VACANT = [14, 16, 19, 20, 22, 23, 24, 32]
SIGNAL = [33, 31]                     # target + one 8-VSB control
h = P.lock_holder()
if h:
    print(f"REFUSING: radio held by {h.get('owner')}"); sys.exit(3)
base = {r["rf"]: r for r in json.load(open(os.path.join(ROOT,"data","survey.json")))["results"]}
peak = json.load(open(os.path.join(HERE, "e64_antenna_probe.json")))["phases"]["as_found"]
peak = {int(k): v for k, v in peak.items()}
nf_peak = {d["rf"]: d["rms_now"] for d in json.load(open(os.path.join(HERE,"e64_noisefloor.json")))}

sdr = SoapySDR.Device("driver=sdrplay")
sdr.setSampleRate(SOAPY_SDR_RX, 0, P.FS); time.sleep(0.2)
sdr.setAntenna(SOAPY_SDR_RX, 0, "Antenna B")
try: sdr.setGainMode(SOAPY_SDR_RX, 0, False)
except Exception: pass
sdr.setGain(SOAPY_SDR_RX, 0, "IFGR", 32)
try: sdr.writeSetting("rfgain_sel", "2")
except Exception: pass
print(f"bias-T readback (must still be as-found 'false'): {P.bias_state(sdr)!r}")
st = sdr.setupStream(SOAPY_SDR_RX, SOAPY_SDR_CS16); sdr.activateStream(st)
buf = np.empty(2 * 65536, np.int16)
res = {}
try:
    for rf in VACANT + SIGNAL:
        res[rf] = P.measure(sdr, st, buf, rf)
finally:
    try: sdr.deactivateStream(st); sdr.closeStream(st)
    except Exception: pass
    del sdr

print("\nVACANT (noise floor)   8/08 ->  08:49peak  ->    now      vs8/08   vsPeak")
r08, rpk = [], []
for rf in VACANT:
    o = base.get(rf, {}).get("rms", 0.0); pk = nf_peak.get(rf, 0.0); n = res[rf]["rms"]
    r08.append(n/o if o else 1); rpk.append(n/pk if pk else 1)
    print(f"  RF{rf:<3} {o:9.0f} -> {pk:9.0f} -> {n:9.0f}   "
          f"{20*np.log10(n/o) if o else 0:+6.1f}dB {20*np.log10(n/pk) if pk else 0:+6.1f}dB")
m08, mpk = float(np.median(r08)), float(np.median(rpk))
print(f"\n  MEDIAN vs 8/08 healthy : {m08:.2f}x ({20*np.log10(m08):+.1f} dB)")
print(f"  MEDIAN vs 08:49 peak   : {mpk:.2f}x ({20*np.log10(mpk):+.1f} dB)")
print("\nSIGNAL     8/08 rms/peak ->  08:49 rms/peak ->  now rms/peak")
for rf in SIGNAL:
    b0 = base.get(rf, {}); p0 = peak.get(rf, {})
    print(f"  RF{rf:<3} {b0.get('rms',0):8.0f}/{b0.get('peak_ratio',0):5.1f} -> "
          f"{p0.get('rms',0):8.0f}/{p0.get('peak_ratio',0):5.1f} -> "
          f"{res[rf]['rms']:8.0f}/{res[rf]['peak_ratio']:5.1f}")
json.dump({rf: res[rf] for rf in res}, open(os.path.join(HERE,"e64_remeasure.json"),"w"), indent=1, default=str)
print(f"\nVERDICT: floor is {'STILL ELEVATED' if m08 > 1.4 else 'BACK NEAR BASELINE'} vs 8/08")
