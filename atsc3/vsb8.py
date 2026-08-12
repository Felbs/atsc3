"""Minimal ATSC 1.0 8-VSB synthesizer -- used ONLY as a false-alarm test vector.

The fleet has no 8-VSB transmitter (checked: everything ATSC-1.0 on this rig is
receive-side), so this exists purely so the bootstrap detector's false-alarm
rate can be measured reproducibly on any machine against a signal that looks
like the real neighbours of an ATSC 3.0 carrier.  It is NOT a compliant
modulator -- no RS/trellis coding, no interleaver, no real MPEG payload.  It
gets the things a detector could trip on right: symbol rate, segment/field sync
structure, 8-level PAM, RRC pulse shape, pilot, and vestigial sideband.

Parameters from ATSC A/53 Part 2 (cribbed in <archive>/atsc-8vsb\\GRID.md):
  symbol rate 4.5/286 * 684 = 10.762238 MHz
  832 symbols/segment (4 sync + 828 data), 313 segments/field
  8 levels -7..+7, pilot = +1.25 DC offset
  RRC roll-off 0.1152, signal bandwidth 5.38 MHz, pilot 2.69 MHz below center
"""

from __future__ import annotations

from fractions import Fraction

import numpy as np
from scipy.signal import resample_poly, hilbert

SYMBOL_RATE = 4.5e6 / 286.0 * 684.0     # 10_762_237.76 Hz
SEG_SYMBOLS = 832
SEGS_PER_FIELD = 313
PILOT_DC = 1.25
ROLLOFF = 0.1152
BW = 5.38e6


def _rrc(n_taps: int, beta: float, sps: float) -> np.ndarray:
    t = (np.arange(n_taps) - (n_taps - 1) / 2.0) / sps
    out = np.empty_like(t)
    for i, ti in enumerate(t):
        if abs(ti) < 1e-9:
            out[i] = 1.0 - beta + 4.0 * beta / np.pi
        elif abs(abs(4.0 * beta * ti) - 1.0) < 1e-9:
            out[i] = (beta / np.sqrt(2.0)) * (
                (1 + 2 / np.pi) * np.sin(np.pi / (4 * beta))
                + (1 - 2 / np.pi) * np.cos(np.pi / (4 * beta)))
        else:
            num = (np.sin(np.pi * ti * (1 - beta))
                   + 4 * beta * ti * np.cos(np.pi * ti * (1 + beta)))
            den = np.pi * ti * (1 - (4 * beta * ti) ** 2)
            out[i] = num / den
    return out / np.sqrt(np.sum(out ** 2))


def synth_8vsb(n_samples: int, fs_out: float, seed: int = 0,
               snr_db: float | None = None) -> np.ndarray:
    """Return `n_samples` of complex-baseband 8-VSB at `fs_out`, unit RMS."""
    rng = np.random.default_rng(seed)
    sps = 4                                   # symbols oversampled 4x internally
    fs_int = SYMBOL_RATE * sps
    n_sym = int(np.ceil(n_samples * fs_int / fs_out / sps)) + 4 * SEG_SYMBOLS

    n_seg = int(np.ceil(n_sym / SEG_SYMBOLS)) + 1
    syms = rng.choice(np.array([-7, -5, -3, -1, 1, 3, 5, 7], dtype=np.float64),
                      size=n_seg * SEG_SYMBOLS)
    seg = syms.reshape(n_seg, SEG_SYMBOLS)
    seg[:, 0:4] = np.array([5.0, -5.0, -5.0, 5.0])          # segment sync
    # crude field sync: every 313th segment is a deterministic +-5 pattern
    fs_rng = np.random.default_rng(12345)
    fsync = fs_rng.choice(np.array([-5.0, 5.0]), size=SEG_SYMBOLS)
    fsync[0:4] = np.array([5.0, -5.0, -5.0, 5.0])
    seg[::SEGS_PER_FIELD, :] = fsync
    syms = seg.reshape(-1)

    up = np.zeros(len(syms) * sps)
    up[::sps] = syms + PILOT_DC                              # pilot
    x = np.convolve(up, _rrc(129, ROLLOFF, sps), mode="same")

    # real -> vestigial sideband complex envelope, pilot 2.69 MHz below center
    a = hilbert(x)                                           # upper sideband
    t = np.arange(len(a)) / fs_int
    a = a * np.exp(-2j * np.pi * (BW / 2.0) * t)

    # NB: never hand resample_poly the raw sample counts -- 6144000/43048951
    # builds a multi-million-tap filter and hangs.  Rational-approximate first;
    # a sub-ppm symbol-rate error is irrelevant for a false-alarm vector.
    fr = Fraction(fs_out / fs_int).limit_denominator(2000)
    y = resample_poly(a, fr.numerator, fr.denominator)
    y = y[:n_samples]
    if len(y) < n_samples:
        y = np.concatenate([y, np.zeros(n_samples - len(y), dtype=y.dtype)])
    y = y / np.sqrt(np.mean(np.abs(y) ** 2))
    if snr_db is not None:
        y = y + (rng.normal(size=n_samples) + 1j * rng.normal(size=n_samples)) \
            * np.sqrt(10 ** (-snr_db / 10.0) / 2.0)
    return y.astype(np.complex64)
