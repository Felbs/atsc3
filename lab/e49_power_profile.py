#!/usr/bin/env python3
"""E49 -- time-domain power profile of a cs16 capture, 1 ms resolution.

Distinguishes bursty interference / gain steps / dropouts from stationary
noise, without any demodulation.  Reads the file in chunks; prints a coarse
summary and writes the ms-grid RMS to .npy for later plotting.
"""
import argparse
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("capture")
ap.add_argument("--rate", type=float, required=True)
ap.add_argument("--npy")
a = ap.parse_args()

n_ms = int(a.rate / 1000)          # complex samples per ms
chunk_ms = 2000
rms = []
with open(a.capture, "rb") as fh:
    while True:
        raw = np.fromfile(fh, dtype=np.int16, count=2 * n_ms * chunk_ms)
        if not len(raw):
            break
        raw = raw[: (len(raw) // (2 * n_ms)) * 2 * n_ms]
        if not len(raw):
            break
        x = raw.astype(np.float32).reshape(-1, 2 * n_ms)
        rms.append(np.sqrt(np.mean(x * x, axis=1)))
rms = np.concatenate(rms)
db = 20 * np.log10(np.maximum(rms, 1e-3))
med = np.median(db)
print(f"{len(rms)} ms bins;  median {med:.2f} dB (rel int16 counts)")
print(f"  p1 {np.percentile(db,1):.2f}  p10 {np.percentile(db,10):.2f}  "
      f"p50 {med:.2f}  p90 {np.percentile(db,90):.2f}  "
      f"p99 {np.percentile(db,99):.2f}  max {db.max():.2f}")
hi = np.flatnonzero(db > med + 3)
lo = np.flatnonzero(db < med - 3)
print(f"  bins >3 dB above median: {len(hi)}  "
      f"({len(hi)/len(db)*100:.2f}%)   below: {len(lo)}")
for name, idx in (("HIGH", hi[:20]), ("LOW", lo[:20])):
    if len(idx):
        print(f"  first {name} bins (s): "
              + ", ".join(f"{i/1000:.3f}" for i in idx))
# 1-second aggregate to see slow structure
sec = db[: len(db)//1000*1000].reshape(-1, 1000).mean(axis=1)
print("  per-second mean dB:")
for i in range(0, len(sec), 10):
    row = " ".join(f"{v:6.2f}" for v in sec[i:i+10])
    print(f"    {i:3d}s  {row}")
if a.npy:
    np.save(a.npy, rms)
    print(f"wrote {a.npy}")
