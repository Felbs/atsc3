#!/usr/bin/env python3
"""e64_ingress_id.py -- what IS the 485-535 MHz ingress? (E64)

A signature helps the operator identify the device. Three questions:
  1. spectrum shape -- discrete carriers, a comb, or wideband hash?
  2. if a COMB, what is the line spacing? (a switching supply's comb
     spacing IS its switching frequency; the fleet found 1.5-2 kHz
     pickets on HF this way)
  3. is it pulsed? (envelope periodicity -> a repetition rate)
Offline, no radio.
"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
FS = 6.912e6
p = os.path.join(HERE, "e64_gain", "ingress_510MHz.cs16")
raw = np.fromfile(p, dtype=np.int16)
iq = (raw[0::2].astype(np.float32) + 1j * raw[1::2].astype(np.float32)) / 32768.0
print(f"capture: {len(iq)} samples = {len(iq)/FS:.2f} s @ {FS/1e6:.3f} Msps "
      f"centred 510 MHz")

# ---- 1. spectrum
nfft = 1 << 16
nseg = min(60, len(iq) // nfft)
acc = np.zeros(nfft)
w = np.hanning(nfft)
for k in range(nseg):
    seg = iq[k * nfft:(k + 1) * nfft] * w
    acc += np.abs(np.fft.fftshift(np.fft.fft(seg))) ** 2
psd = 10 * np.log10(acc / nseg + 1e-20)
freqs = np.fft.fftshift(np.fft.fftfreq(nfft, 1 / FS)) + 510e6
med = float(np.median(psd))
peak_idx = int(np.argmax(psd))
print(f"\nPSD: median {med:.1f} dB, peak {psd[peak_idx]:.1f} dB at "
      f"{freqs[peak_idx]/1e6:.4f} MHz  (peak-to-median {psd[peak_idx]-med:.1f} dB)")
# how much of the band is >10 dB above median -> hash vs discrete
frac = float(np.mean(psd > med + 10))
print(f"fraction of bins >10 dB above median: {frac*100:.2f}%  -> "
      f"{'WIDEBAND HASH' if frac > 0.25 else 'DISCRETE lines / comb'}")
# strongest lines
order = np.argsort(psd)[::-1]
lines, seen = [], []
for i in order:
    f = freqs[i]
    if any(abs(f - s) < 20e3 for s in seen):
        continue
    seen.append(f); lines.append((f, psd[i]))
    if len(lines) >= 12:
        break
print("\nstrongest distinct lines (>=20 kHz apart):")
for f, v in lines:
    print(f"   {f/1e6:10.4f} MHz   {v-med:+6.1f} dB above median")
# ---- 2. comb spacing from the line set
fs_sorted = sorted(f for f, _ in lines)
d = np.diff(fs_sorted)
if len(d):
    print(f"\nspacing between adjacent strong lines (kHz): "
          f"{', '.join(f'{x/1e3:.1f}' for x in d)}")
    print(f"  median spacing {np.median(d)/1e3:.1f} kHz")
# ---- 3. pulsed? envelope periodicity
env = np.abs(iq[:int(FS)])                       # 1 s
env -= env.mean()
n = 1 << int(np.ceil(np.log2(len(env))))
ac = np.fft.irfft(np.abs(np.fft.rfft(env, n)) ** 2)[:len(env) // 2]
ac /= (ac[0] + 1e-20)
lo = int(FS / 5000)                              # ignore <lag for >5 kHz
cand = ac[lo:int(FS / 20)]                       # 20 Hz .. 5 kHz
if len(cand):
    k = int(np.argmax(cand)) + lo
    print(f"\nenvelope autocorrelation: strongest period {FS/k:.1f} Hz "
          f"(corr {cand[k-lo]:.3f})  -> "
          f"{'PULSED/periodic' if cand[k-lo] > 0.05 else 'no strong periodicity'}")
np.save(os.path.join(HERE, "e64_ingress_psd.npy"), psd.astype(np.float32))
print(f"\nPSD saved -> lab/e64_ingress_psd.npy")
