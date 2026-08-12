#!/usr/bin/env python3
"""E76 -- is there an ATSC 3.0 signal in the new RF25 captures at all?

m8_l1 found no verifiable L1-Basic in 40 attempts across 20 s of
rf25_fox_g4_0810_2254.cs16, while the SAME command on E49's own capture
(rf25_amped_0808.cs16) verifies at t=0.5 s.  That is not "intermittent", so
before any margin is reported, go under L1 to the two instruments that cannot
be fooled by a decode-chain assumption:

  1. the SPECTRUM -- occupied bandwidth and where the energy actually sits
  2. the BOOTSTRAP detector (atsc3.bootstrap, gated 32/32 down to -18 dB)

The old capture is the positive control: whatever these instruments say about
the new files, they must say the right thing about the one we know decodes.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from atsc3 import bootstrap as BS  # noqa: E402

FS_CAP = 6.912e6
CAPS = [
    ("data/fox25/rf25_fox_g4_0810_2254.cs16", "NEW gain 4 (180 s)"),
    ("data/fox25/rf25_fox_g2_0810_2257.cs16", "NEW gain 2 (90 s)"),
    ("data/fox25/rf25_amped_0808.cs16", "E49 control (decodes)"),
]


def load(path, sec, start=0.0):
    n = int(sec * FS_CAP)
    off = int(start * FS_CAP) * 4
    with open(path, "rb") as fh:
        fh.seek(off)
        b = fh.read(n * 4)
    a = np.frombuffer(b, dtype="<i2").astype(np.float32)
    return (a[0::2] + 1j * a[1::2]) / 32768.0


def spectrum(x, nfft=4096):
    w = np.hanning(nfft)
    nseg = min(400, len(x) // nfft)
    acc = np.zeros(nfft)
    for k in range(nseg):
        seg = x[k * nfft:(k + 1) * nfft] * w
        acc += np.abs(np.fft.fftshift(np.fft.fft(seg))) ** 2
    acc /= nseg
    return 10 * np.log10(acc + 1e-30)


def main():
    for rel, label in CAPS:
        p = os.path.join(ROOT, rel)
        if not os.path.exists(p):
            print("MISSING", rel)
            continue
        print("=" * 74)
        print("%s   %s" % (os.path.basename(rel), label))
        x = load(p, 2.0, start=1.0)
        pw = float(np.mean(np.abs(x) ** 2))
        print("  mean power %.3e   rms %.5f" % (pw, np.sqrt(pw)))

        S = spectrum(x)
        f = (np.arange(len(S)) - len(S) // 2) * (FS_CAP / len(S)) / 1e6
        # occupied bandwidth: where the PSD is within 10 dB of the in-band median
        peak = np.percentile(S, 99)
        occ = f[S > peak - 10]
        print("  PSD peak %.1f dB, median %.1f dB, peak-median %.1f dB"
              % (peak, np.median(S), peak - np.median(S)))
        if len(occ):
            print("  energy within 10 dB of peak spans %.3f .. %.3f MHz "
                  "(%.3f MHz wide)" % (occ.min(), occ.max(), occ.max() - occ.min()))
        # coarse shape: power in the central 5.8 MHz vs the edges
        cen = np.mean(10 ** (S[np.abs(f) < 2.75] / 10))
        edge = np.mean(10 ** (S[np.abs(f) > 3.1] / 10))
        print("  in-band(+-2.75MHz) / out-of-band(>3.1MHz) = %.1f dB"
              % (10 * np.log10(cen / edge)))

        # ---- bootstrap: resample 6.912 -> 6.144 Msps as the detector wants
        y = load(p, 1.5, start=1.0)
        m = int(len(y) * BS.FS / FS_CAP)
        idx = np.arange(m) * (FS_CAP / BS.FS)
        i0 = idx.astype(np.int64)
        fr = (idx - i0).astype(np.float32)
        i1 = np.minimum(i0 + 1, len(y) - 1)
        yr = y[i0] * (1 - fr) + y[i1] * fr
        try:
            cands = BS.acquire_matched(yr, n_symbols=4, top_k=3)
        except Exception as exc:                                # noqa: BLE001
            print("  bootstrap: ERROR %s: %s" % (type(exc).__name__, exc))
            continue
        if not cands:
            print("  BOOTSTRAP: none found")
        else:
            for c in cands[:3]:
                keys = {k: c[k] for k in list(c)[:7]}
                print("  bootstrap candidate:", keys)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
