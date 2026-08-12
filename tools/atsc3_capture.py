#!/usr/bin/env python3
"""atsc3_capture.py -- write a raw cs16 IQ capture of one RF channel.

Minimal by design: open SDR, tune, stream to .cs16 (interleaved int16 IQ,
the format atsc3.capture already reads). Single-tenant: only run when the
chain is down. mkdir-only, never deletes anything.
"""
import argparse, json, os, sys, time
import numpy as np
import SoapySDR
from SoapySDR import SOAPY_SDR_RX, SOAPY_SDR_CS16

SoapySDR.SoapySDR_setLogLevel(SoapySDR.SOAPY_SDR_FATAL)

def center_hz(rf):
    lo = (174 + (rf - 7) * 6) if rf < 14 else (470 + (rf - 14) * 6)
    return (lo + 3.0) * 1e6

ap = argparse.ArgumentParser()
ap.add_argument("--rf", type=int, required=True)
ap.add_argument("--secs", type=float, default=120.0)
ap.add_argument("--rate", type=float, default=6.912e6)
ap.add_argument("--ant", default="Antenna B")
ap.add_argument("--ifgr", type=int, default=32)
# E67: this tool NEVER set rfgain_sel, so it silently inherited whatever
# the last process left in the device -- and every RF33 capture we own is
# clipped at full scale as a result (hit_rf33.cs16 outright VOID). The LNA
# state is now an explicit, READ-BACK parameter of the capture, because an
# unverified writeSetting is a hope, not a setting.
ap.add_argument("--rfgain", type=int, default=None,
                help="rfgain_sel (LNA state). Higher = MORE attenuation. "
                     "Omit to leave the device as found (the old, unsafe "
                     "behaviour -- the value is still read back and logged).")
ap.add_argument("--allow-clipping", action="store_true",
                help="bank the capture even if it clips (default: refuse)")
# E86: the yagi's LNA on Antenna B is powered by the SDR's bias-T, and the
# atsc3 chain has NEVER written biasT_ctrl -- it inherits whatever the last
# process on the fleet left in the device.  E64 found it 'false' and E86
# measured what that costs: unchanged broadband rms but bootstrap peak
# ratio 64.6 -> 81.0 on RF25 with it on.  An amplifier you did not ask for
# is exactly as unmeasurable as one you did; make it an explicit, read-back
# parameter of the capture like rfgain_sel.
ap.add_argument("--biast", default=None, choices=("on", "off"),
                help="power the Antenna B bias-T (yagi LNA). Omit to leave "
                     "the device as found -- the state is read back and "
                     "recorded either way.")
ap.add_argument("--out", required=True)
a = ap.parse_args()

os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
sdr = SoapySDR.Device("driver=sdrplay")
sdr.setSampleRate(SOAPY_SDR_RX, 0, a.rate)
time.sleep(0.2)
fs = float(sdr.getSampleRate(SOAPY_SDR_RX, 0))
sdr.setAntenna(SOAPY_SDR_RX, 0, a.ant)
try:
    sdr.setGainMode(SOAPY_SDR_RX, 0, False)
except Exception:
    pass
sdr.setGain(SOAPY_SDR_RX, 0, "IFGR", a.ifgr)
if a.rfgain is not None:
    try:
        sdr.writeSetting("rfgain_sel", str(a.rfgain))
    except Exception as e:
        print(f"WARNING: rfgain_sel write failed: {e}")
if a.biast is not None:
    try:
        sdr.writeSetting("biasT_ctrl", "true" if a.biast == "on" else "false")
    except Exception as e:
        print(f"WARNING: biasT_ctrl write failed: {e}")
time.sleep(0.3)
try:
    rb_bias = str(sdr.readSetting("biasT_ctrl"))
except Exception:
    rb_bias = "<unreadable>"
if a.biast is not None and rb_bias != ("true" if a.biast == "on" else "false"):
    print(f"WARNING: bias-T readback {rb_bias!r} != requested {a.biast!r} "
          f"-- the LNA is NOT in the state this capture claims")
# READ BACK what the device actually holds -- never what we asked for
try:
    rb_rf = str(sdr.readSetting("rfgain_sel"))
except Exception:
    rb_rf = "<unreadable>"
rb_if = float(sdr.getGain(SOAPY_SDR_RX, 0, "IFGR"))
rb_ant = sdr.getAntenna(SOAPY_SDR_RX, 0)
print(f"gain readback: rfgain_sel={rb_rf} IFGR={rb_if:.0f} ant={rb_ant} "
      f"rate={fs/1e6:.3f} Msps biasT={rb_bias}")
sdr.setFrequency(SOAPY_SDR_RX, 0, center_hz(a.rf))
st = sdr.setupStream(SOAPY_SDR_RX, SOAPY_SDR_CS16)
sdr.activateStream(st)
buf = np.empty(2 * 65536, np.int16)
n_want = int(a.secs * fs)
got = 0
t0 = time.time()
peak = 0            # max |sample| seen, both I and Q
n_full = 0          # samples at/near full scale (the clipping evidence)
n_over = 0          # SDR overflow events
with open(a.out, "wb") as f:
    # flush retune transient
    while time.time() - t0 < 0.3:
        sdr.readStream(st, [buf], 65536, timeoutUs=300000)
    t_start = time.time()
    while got < n_want:
        r = sdr.readStream(st, [buf], 65536, timeoutUs=500000)
        if r.ret > 0:
            blk = buf[:2 * r.ret]
            f.write(blk.tobytes())
            got += r.ret
            # headroom is measured on the way past -- one pass, no reread
            m = int(np.abs(blk).max())
            if m > peak:
                peak = m
            n_full += int(np.count_nonzero(np.abs(blk) >= 32000))
        elif r.ret == -4:
            n_over += 1                   # overflow: keep going, but count it
        elif r.ret < 0:
            print(f"readStream {r.ret}", flush=True)
wall = time.time() - t_start
sdr.deactivateStream(st); sdr.closeStream(st); del sdr

frac = 100.0 * n_full / max(1, 2 * got)
dbfs = 20.0 * np.log10(peak / 32768.0) if peak else -999.0
print(f"CAPTURED {got} samples ({got/fs:.1f}s) in {wall:.1f}s wall -> {a.out}")
print(f"  headroom: max|s|={peak} ({dbfs:+.1f} dBFS)  "
      f"full-scale samples={frac:.4f}%  overflows={n_over}")
ok = True
# capture-integrity law: samples must match wall clock or the capture is VOID
if abs(got / fs - wall) > 0.05 * wall:
    print(f"  VOID: samples/fs={got/fs:.2f}s != wall={wall:.2f}s (gapped)")
    ok = False
else:
    print(f"  integrity OK: samples/fs={got/fs:.2f}s == wall={wall:.2f}s")
# E67: a CLIPPED capture is void too -- it measures the ADC rail, not the
# air, and every analyser downstream will quietly read the damage as a bad
# link (m8_l1 failed on exactly this and I nearly blamed the transmitter).
if frac > 0.001 or peak >= 32767:
    print(f"  CLIPPED -> VOID. Re-run with more attenuation "
          f"(--rfgain higher, and/or --ifgr higher).")
    ok = False
else:
    print(f"  no clipping")
print(f"  VERDICT: {'USABLE' if ok else 'VOID -- do not bank this'}")

# ---------------------------------------------------------------- sidecar
# E83: A CAPTURE THAT DOES NOT RECORD WHAT IT USED IS HALF A MEASUREMENT.
# carrier_gain.json carries a null-forever entry for exactly this reason:
# a capture was taken, the gain it ran at was never written down, and the
# number could never be attributed afterwards. Every capture now lands with
# a sidecar carrying the CONFIRMED state (readback, not intention), the
# headroom, the integrity arithmetic, and a cheap quality proxy -- so a
# capture found on disk months later can still be ranked against another.
#
# bias-T earns its own field: it does NOT persist across a device close
# (measured 8/11), and the atsc3 chain never writes it, so the live yagi
# LNA has always been unpowered. A capture that does not record the bias
# state cannot be compared with one taken under a different one.
# THE SIDECAR MAY NEVER COST US THE CAPTURE. It runs AFTER the bytes are
# safely on disk and the verdict is decided, and the whole block is
# wrapped: a bug in a bonus feature must not turn a good 5 GB capture into
# a VOID one (fadecap keys on this script's exit code). Worst case the
# sidecar is missing; never that the capture is lost.
side = os.path.splitext(a.out)[0] + ".json"
quality = {}
try:
    try:
        # cheap, no extra radio time: bootstrap correlation off ~1.5 s of
        # the capture we already hold. Not an SNR, but a per-capture number
        # that ranks two captures of the same carrier against each other.
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))), "lab"))
        import m2_pilots as MP
        from atsc3 import bootstrap as _bs
        n_q = min(got, int(1.5 * fs))
        raw = np.fromfile(a.out, dtype=np.int16, count=2 * n_q)
        iq = (raw[0::2].astype(np.float32)
              + 1j * raw[1::2].astype(np.float32)) / 32768.0
        quality["in_rms"] = float(np.sqrt((np.abs(iq) ** 2).mean()) * 32768)
        hits = MP.find_bootstraps(MP.resample_to(iq, fs, _bs.FS))
        quality["bootstrap_hits"] = len(hits)
        quality["peak_ratio"] = (float(hits[0].get("mean_peak_ratio", 0.0))
                                 if hits else 0.0)
        quality["cfo_hz"] = (float(hits[0].get("fine_cfo_hz") or 0.0)
                             if hits else None)
    except Exception as e:
        quality["error"] = f"{type(e).__name__}: {e}"

    json.dump({
        "capture": os.path.basename(a.out),
        "when": time.strftime("%Y-%m-%d %H:%M:%S"),
        "rf": a.rf, "centre_hz": center_hz(a.rf), "rate_hz": fs,
        "requested": {"rfgain_sel": a.rfgain, "ifgr": a.ifgr,
                      "biast": a.biast, "ant": a.ant},
        "readback": {"rfgain_sel": rb_rf, "ifgr": rb_if,
                     "biast": rb_bias, "ant": rb_ant},
        "headroom": {"peak_abs": peak, "peak_dbfs": round(dbfs, 2),
                     "full_scale_pct": round(frac, 6), "overflows": n_over},
        "integrity": {"samples": got, "seconds": round(got / fs, 3),
                      "wall_s": round(wall, 3),
                      "ok": bool(abs(got / fs - wall) <= 0.05 * wall)},
        "quality": quality,
        "verdict": "USABLE" if ok else "VOID",
    }, open(side, "w"), indent=1)
    print(f"  sidecar -> {side}")
    if quality.get("peak_ratio"):
        print(f"  quality: bootstrap peak_ratio {quality['peak_ratio']:.1f}, "
              f"in_rms {quality.get('in_rms', 0):.0f}, "
              f"{quality.get('bootstrap_hits', 0)} hits")
except Exception as _e:
    # a bonus feature must never turn a good capture into a VOID one
    print(f"  ! sidecar step failed ({type(_e).__name__}: {_e}) -- "
          f"the capture itself is unaffected")

sys.exit(0 if (ok or a.allow_clipping) else 1)
