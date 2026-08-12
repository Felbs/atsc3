#!/usr/bin/env python3
"""e64_gain_sweep.py -- is the front end OVERLOADING on the low-UHF ingress?

E64 (8/10). The ingress peaks at 485-535 MHz; RF33 is at 587 MHz. Strong
OUT-OF-BAND energy degrading an IN-BAND link is the signature of front-end
compression / reciprocal mixing, and the cure is LESS front-end gain, not
more filtering. Hypothesis: raising rfgain_sel (less LNA) and/or IFGR
(more IF attenuation) IMPROVES RF33's FEC while the ingress is strong.

DESIGN -- the channel drifts on a ~10 min timescale and the good windows
are ~10 min long, so a slow sequential sweep measures the WEATHER, not the
gain. Instead every setting is captured back-to-back inside ~2 minutes and
decoded OFFLINE afterwards with the real chain (`atsc3 watch --capture`),
so all settings see the same channel. The baseline is captured FIRST and
LAST: if those two disagree, the sweep is void and says so.

Every gain write is READ BACK (the unverified-writeSetting law) and the
as-found state is restored at the end.
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

ap = argparse.ArgumentParser()
ap.add_argument("--secs", type=float, default=12.0)
ap.add_argument("--rf", type=int, default=33)
ap.add_argument("--outdir", default=os.path.join(HERE, "e64_gain"))
ap.add_argument("--ingress", action="store_true",
                help="also capture 3 s at 510 MHz to characterise the ingress")
a = ap.parse_args()

h = P.lock_holder()
if h:
    print(f"REFUSING: radio held by {h.get('owner')} (pid {h.get('pid')})")
    sys.exit(3)
os.makedirs(a.outdir, exist_ok=True)          # mkdir-only, never clears

sdr = SoapySDR.Device("driver=sdrplay")
sdr.setSampleRate(SOAPY_SDR_RX, 0, P.FS); time.sleep(0.2)
sdr.setAntenna(SOAPY_SDR_RX, 0, "Antenna B")
try: sdr.setGainMode(SOAPY_SDR_RX, 0, False)
except Exception: pass

# what does this device actually allow?
info = {}
try:
    for k in sdr.getSettingInfo():
        if k.key == "rfgain_sel":
            info["rfgain_sel"] = {"options": list(k.options),
                                  "range": str(k.range)}
except Exception as e:
    info["err"] = str(e)
try:
    r = sdr.getGainRange(SOAPY_SDR_RX, 0, "IFGR")
    info["IFGR"] = f"{r.minimum()}..{r.maximum()}"
except Exception: pass
print(f"device gain capability: {json.dumps(info)}", flush=True)

AS_FOUND = (2, 32)                 # what the chain runs today
# (rfgain_sel, IFGR). Baseline first AND last = the drift control.
PLAN = [(2, 32), (0, 32), (3, 32), (4, 32), (2, 26), (2, 40), (4, 40), (2, 32)]

def apply(rfsel, ifgr):
    sdr.setGain(SOAPY_SDR_RX, 0, "IFGR", ifgr)
    try: sdr.writeSetting("rfgain_sel", str(rfsel))
    except Exception as e: return None, f"write failed {e}"
    time.sleep(0.35)
    got_if = float(sdr.getGain(SOAPY_SDR_RX, 0, "IFGR"))
    try: got_rf = str(sdr.readSetting("rfgain_sel"))
    except Exception: got_rf = "<unreadable>"
    return (got_rf, got_if), None

st = sdr.setupStream(SOAPY_SDR_RX, SOAPY_SDR_CS16); sdr.activateStream(st)
buf = np.empty(2 * 65536, np.int16)
sdr.setFrequency(SOAPY_SDR_RX, 0, P.center_hz(a.rf))
time.sleep(0.3)
manifest = []
try:
    for i, (rfsel, ifgr) in enumerate(PLAN):
        rb, err = apply(rfsel, ifgr)
        if err:
            print(f"  [{i}] rf{rfsel}/if{ifgr}: {err}"); continue
        t0 = time.time()
        while time.time() - t0 < 0.4:                 # settle
            sdr.readStream(st, [buf], 65536, timeoutUs=300000)
        path = os.path.join(a.outdir, f"s{i}_rf{rfsel}_if{ifgr}.cs16")
        n_want = int(a.secs * P.FS); got = 0
        with open(path, "wb") as f:
            deadline = time.time() + a.secs * 2 + 5
            while got < n_want and time.time() < deadline:
                r = sdr.readStream(st, [buf], 65536, timeoutUs=500000)
                if r.ret > 0:
                    n = min(r.ret, n_want - got)
                    f.write(buf[:2 * n].tobytes()); got += n
        rms = None
        manifest.append(dict(i=i, rfgain_sel=rfsel, ifgr=ifgr, path=path,
                             readback_rfgain=rb[0], readback_ifgr=rb[1],
                             samples=got, secs=got / P.FS))
        print(f"  [{i}] rfgain_sel={rfsel}(rb {rb[0]}) IFGR={ifgr}"
              f"(rb {rb[1]:.0f})  {got} samples = {got/P.FS:.1f}s", flush=True)
    if a.ingress:
        rb, _ = apply(*AS_FOUND)
        sdr.setFrequency(SOAPY_SDR_RX, 0, 510e6)
        time.sleep(0.4)
        t0 = time.time()
        while time.time() - t0 < 0.3:
            sdr.readStream(st, [buf], 65536, timeoutUs=300000)
        path = os.path.join(a.outdir, "ingress_510MHz.cs16")
        n_want = int(3.0 * P.FS); got = 0
        with open(path, "wb") as f:
            deadline = time.time() + 10
            while got < n_want and time.time() < deadline:
                r = sdr.readStream(st, [buf], 65536, timeoutUs=500000)
                if r.ret > 0:
                    n = min(r.ret, n_want - got)
                    f.write(buf[:2 * n].tobytes()); got += n
        print(f"  ingress capture @510 MHz: {got} samples "
              f"({got/P.FS:.1f}s) -> {path}", flush=True)
        manifest.append(dict(i="ingress", freq=510e6, path=path, samples=got))
finally:
    rb, _ = apply(*AS_FOUND)
    print(f"restored as-found rfgain_sel={AS_FOUND[0]} IFGR={AS_FOUND[1]} "
          f"-> readback {rb}", flush=True)
    try: sdr.deactivateStream(st); sdr.closeStream(st)
    except Exception: pass
    del sdr
json.dump(dict(when=time.strftime("%Y-%m-%d %H:%M:%S"), rf=a.rf,
               as_found=AS_FOUND, capability=info, runs=manifest),
          open(os.path.join(a.outdir, "manifest.json"), "w"), indent=1)
print(f"\nwrote {os.path.join(a.outdir,'manifest.json')}")
