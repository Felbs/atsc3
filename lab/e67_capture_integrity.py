#!/usr/bin/env python3
"""E67 capture integrity: sample count vs wall x fs, and CLIPPING.

A sibling agent found that tools/atsc3_capture.py never sets rfgain_sel, which
produced clipped captures that made a healthy link look dead (E64/E67 splitter
work).  Before trusting any banked RF33 capture for an encryption verdict,
check it: max|s| near 32767 or a nonzero full-scale fraction means VOID.

Read-only; samples the file rather than loading it whole.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CAPS = [
    ("rate6912_rf33.cs16", 6.912e6),
    ("long_rf33.cs16", 6.912e6),
    ("hit_rf33.cs16", 6.912e6),
]
NCHUNK = 24
CHUNK = 1 << 20  # int16 values per probe chunk


def check(path, fs):
    sz = os.path.getsize(path)
    nsamp = sz // 4  # cs16 = 2 x int16 per complex sample
    rep = {
        "file": os.path.basename(path),
        "bytes": sz,
        "complex_samples": nsamp,
        "assumed_fs": fs,
        "implied_seconds": round(nsamp / fs, 3),
    }
    fs_int = 0
    mx = 0
    tot = 0
    sumsq = 0.0
    with open(path, "rb") as fh:
        for i in range(NCHUNK):
            off = int(i * (sz - CHUNK * 2) / max(NCHUNK - 1, 1)) // 4 * 4
            if off < 0:
                off = 0
            fh.seek(off)
            b = fh.read(CHUNK * 2)
            if not b:
                break
            a = np.frombuffer(b[: len(b) // 2 * 2], dtype="<i2")
            if a.size == 0:
                break
            mx = max(mx, int(np.abs(a.astype(np.int32)).max()))
            fs_int += int(np.count_nonzero(np.abs(a.astype(np.int32)) >= 32767))
            tot += a.size
            sumsq += float(np.square(a.astype(np.float64)).sum())
    rep["max_abs"] = mx
    rep["fullscale_samples"] = fs_int
    rep["fullscale_frac"] = fs_int / tot if tot else None
    rep["rms"] = round((sumsq / tot) ** 0.5, 1) if tot else None
    rep["headroom_dbfs"] = round(20 * np.log10(mx / 32767.0), 2) if mx else None
    rep["probed_int16"] = tot
    # VOID if any sample sits at full scale
    rep["verdict"] = "VOID (clipped)" if fs_int > 0 else "OK (no clipping)"
    return rep


def main():
    out = []
    for name, fs in CAPS:
        p = os.path.join(HERE, name)
        if not os.path.exists(p):
            out.append({"file": name, "verdict": "MISSING"})
            continue
        out.append(check(p, fs))
    with open(os.path.join(HERE, "e67_capture_integrity.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    hdr = "%-22s %-13s %-9s %-9s %-11s %-9s %s" % (
        "capture", "complex", "secs", "max|s|", "fullscale", "peak dBFS", "verdict")
    print(hdr)
    print("-" * len(hdr))
    for r in out:
        if r.get("verdict") == "MISSING":
            print("%-22s %s" % (r["file"], "MISSING"))
            continue
        print("%-22s %-13d %-9s %-9d %-11s %-9s %s" % (
            r["file"], r["complex_samples"], r["implied_seconds"], r["max_abs"],
            r["fullscale_samples"], r["headroom_dbfs"], r["verdict"]))


if __name__ == "__main__":
    main()
