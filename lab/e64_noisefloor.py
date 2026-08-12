#!/usr/bin/env python3
"""e64_noisefloor.py -- did the NOISE FLOOR rise? (E64 follow-up)

The main probe showed: antenna path fine (8-VSB controls at/above the
8/08 baseline), bias-T irrelevant (A/B changed nothing), yet EVERY ATSC
3.0 carrier's bootstrap peak ratio fell ~40% while in-band rms ROSE.
More power, worse quality = interference, not attenuation.

The crux test uses channels with NO transmitter. Their rms IS the noise
floor, measured through the same antenna and the same gains as the 8/08
survey. If the vacant channels are up, broadband noise has been injected
into the receive path -- the classic USB-3 / switching-supply signature,
and the SDR's USB was just replugged."""
import json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, ROOT)
import SoapySDR
from SoapySDR import SOAPY_SDR_RX, SOAPY_SDR_CS16
SoapySDR.SoapySDR_setLogLevel(SoapySDR.SOAPY_SDR_FATAL)
sys.path.insert(0, HERE)
import e64_antenna_probe as P

VACANT = [14, 16, 19, 20, 22, 23, 24, 32]      # no lock in the 8/08 scan
sdr = SoapySDR.Device("driver=sdrplay")
sdr.setSampleRate(SOAPY_SDR_RX, 0, P.FS); time.sleep(0.2)
sdr.setAntenna(SOAPY_SDR_RX, 0, "Antenna B")
try: sdr.setGainMode(SOAPY_SDR_RX, 0, False)
except Exception: pass
sdr.setGain(SOAPY_SDR_RX, 0, "IFGR", 32)
try: sdr.writeSetting("rfgain_sel", "2")
except Exception: pass
st = sdr.setupStream(SOAPY_SDR_RX, SOAPY_SDR_CS16); sdr.activateStream(st)
buf = np.empty(2 * 65536, np.int16)
base = {r["rf"]: r for r in json.load(open(os.path.join(ROOT,"data","survey.json")))["results"]}
rows = []
try:
    for rf in VACANT:
        m = P.measure(sdr, st, buf, rf)
        o = base.get(rf, {}).get("rms", 0.0)
        rows.append((rf, o, m["rms"]))
        print(f"  RF{rf:<3} vacant: 8/08 rms {o:7.0f} -> now {m['rms']:7.0f}"
              f"  ({(m['rms']/o if o else float('nan')):5.2f}x, "
              f"{20*np.log10(m['rms']/o) if o else 0:+5.1f} dB)", flush=True)
finally:
    try: sdr.deactivateStream(st); sdr.closeStream(st)
    except Exception: pass
    del sdr
r = [n/o for _, o, n in rows if o > 0]
med = float(np.median(r))
print(f"\nNOISE FLOOR median ratio now/8-08 = {med:.2f}  "
      f"({20*np.log10(med):+.1f} dB)")
print("VERDICT:", "NOISE FLOOR HAS RISEN -- broadband ingress into the "
      "receive path" if med > 1.4 else "noise floor unchanged")
json.dump([{"rf": a, "rms_0808": b, "rms_now": c} for a, b, c in rows],
          open(os.path.join(HERE, "e64_noisefloor.json"), "w"), indent=1)
