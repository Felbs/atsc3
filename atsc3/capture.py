"""IQ capture readers and resampling to the fixed 6.144 Msps bootstrap rate.

The fleet convention (see Software-TV-Tuner / a sibling 8-VSB lab lab dirs) is:
  *.cs16  interleaved int16 little-endian  ("ci16_le" in the SigMF sidecars)
  *.cf32  interleaved float32
  *.cu8   interleaved uint8, DC offset 127.5
  *.iq    ambiguous -- caller must say
Sidecars: <name>.json (center_hz / rate / rf) or <name>.sigmf-meta.
"""

from __future__ import annotations

import json
import os
from fractions import Fraction

import numpy as np
from scipy.signal import resample_poly

from .bootstrap import FS

DTYPES = {
    "cs16": (np.int16, 1.0 / 32768.0, 0.0),
    "ci16": (np.int16, 1.0 / 32768.0, 0.0),
    "cf32": (np.float32, 1.0, 0.0),
    "cs8": (np.int8, 1.0 / 128.0, 0.0),
    "cu8": (np.uint8, 1.0 / 128.0, -127.5),
}


def guess_format(path: str) -> str:
    ext = os.path.splitext(path)[1].lower().lstrip(".")
    if ext in DTYPES:
        return ext
    return "cs16"           # fleet default for .iq / unknown


def read_sidecar(path: str) -> dict:
    """Best-effort center frequency / sample rate from a .json or .sigmf-meta."""
    out = {}
    base = os.path.splitext(path)[0]
    for cand, kind in ((base + ".json", "json"), (base + ".sigmf-meta", "sigmf"),
                       (path + ".json", "json")):
        if not os.path.exists(cand):
            continue
        try:
            with open(cand, "r", encoding="utf-8", errors="replace") as fh:
                d = json.load(fh)
        except Exception:
            continue
        if kind == "sigmf":
            g = d.get("global", {})
            if "core:sample_rate" in g:
                out["sample_rate"] = float(g["core:sample_rate"])
            if "core:frequency" in g:
                out["center_hz"] = float(g["core:frequency"])
            for cap in d.get("captures", []) or []:
                if "core:frequency" in cap:
                    out.setdefault("center_hz", float(cap["core:frequency"]))
            if "core:datatype" in g:
                out["datatype"] = g["core:datatype"]
        else:
            for k_src, k_dst in (("center_hz", "center_hz"), ("rate", "sample_rate"),
                                 ("sample_rate", "sample_rate"), ("rf", "rf"),
                                 ("antenna", "antenna"), ("secs", "secs")):
                if k_src in d:
                    out[k_dst] = d[k_src]
        out["sidecar"] = cand
    return out


def n_samples(path: str, fmt: str | None = None) -> int:
    fmt = fmt or guess_format(path)
    dt = np.dtype(DTYPES[fmt][0])
    return os.path.getsize(path) // (dt.itemsize * 2)


def read_block(path: str, start: int, count: int, fmt: str | None = None) -> np.ndarray:
    """Read `count` complex samples starting at complex-sample index `start`."""
    fmt = fmt or guess_format(path)
    dtype, scale, offset = DTYPES[fmt]
    dt = np.dtype(dtype)
    raw = np.fromfile(path, dtype=dt, count=2 * count, offset=start * 2 * dt.itemsize)
    if len(raw) < 2:
        return np.zeros(0, dtype=np.complex64)
    raw = raw[: 2 * (len(raw) // 2)].astype(np.float32)
    if offset:
        raw += offset
    if scale != 1.0:
        raw *= scale
    return (raw[0::2] + 1j * raw[1::2]).astype(np.complex64)


def to_bootstrap_rate(x: np.ndarray, fs_in: float,
                      freq_shift_hz: float = 0.0) -> np.ndarray:
    """Frequency-shift then resample to exactly 6.144 Msps.

    A/321 fixes the bootstrap rate; the detector always works there.  Any
    capture rate is fine as long as it is wide enough to hold the 4.5 MHz
    bootstrap bandwidth -- i.e. fs_in >= ~5 MHz, and >= 6.144 MHz preferred so
    nothing is lost at the band edges.
    """
    if freq_shift_hz:
        t = np.arange(len(x), dtype=np.float64)
        x = x * np.exp(-2j * np.pi * freq_shift_hz * t / fs_in).astype(np.complex64)
    if abs(fs_in - FS) < 1.0:
        return x.astype(np.complex64)
    fr = Fraction(int(round(FS)), int(round(fs_in))).limit_denominator(20000)
    y = resample_poly(x, fr.numerator, fr.denominator)
    return y.astype(np.complex64)
