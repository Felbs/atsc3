"""Run the ATSC 3.0 bootstrap detector over ARCHIVED IQ captures.  No radio.

Usage
-----
    python sweep_archive.py --fleet              # the curated fleet manifest
    python sweep_archive.py --fleet --quick      # first 4 s of each file
    python sweep_archive.py FILE [FILE...] --rate 8e6
    python sweep_archive.py FILE --pn-sweep      # try all PN wiring variants

Each capture is read in chunks, resampled to the fixed 6.144 Msps bootstrap
rate, and scanned with both detectors.  The report is per-file:

    peak_ratio   best peak/background of the folded structural metric.  This is
                 the ASSUMPTION-FREE number.  Threshold 6.0; selftest measured
                 max 1.9 on Gaussian noise and 2.4 on 8-VSB over 100 trials.
    matched      best mean correlation peak ratio against the synthesized
                 reference.  Threshold 8.0.
    verdict      NO_BOOTSTRAP / ATSC3_STRUCTURE_ONLY / ATSC3_CONFIRMED

A capture of a real ATSC 1.0 (8-VSB) station is included as a CONTROL group:
those must come back NO_BOOTSTRAP or the detector is not trustworthy.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from atsc3 import bootstrap as bs
from atsc3.capture import (guess_format, n_samples, read_block, read_sidecar,
                           to_bootstrap_rate)

# RF channel -> centre Hz, from
# <archive>/tools\scan_lab\fixtures\manifest.json
RF_HZ = {7: 177e6, 9: 189e6, 15: 479e6, 21: 515e6, 27: 551e6,
         34: 593e6, 35: 599e6, 36: 605e6}

# The curated sweep.  role=target -> the carriers that never decode;
# role=control -> known-good ATSC 1.0, MUST come back negative.
FLEET = [
    # ---- TARGETS: the "strong energy, never decodes" carriers ----
    ("target", r"<archive>/lab\wl_gate2\rf27_capture.cs16", 27, 8e6),
    ("target", r"<archive>/lab\captures\philips_rf36_strong.cs16", 36, 8e6),
    ("target", r"data\offline_test\rf36.cs16", 36, 8e6),
    ("target", r"<archive>/tools\scan_lab\fixtures\rf27.cf32", 27, 8e6),
    ("target", r"<archive>/tools\scan_lab\fixtures\rf36.cf32", 36, 8e6),
    # ---- CONTROLS: known-good 8-VSB that decodes to a TS ----
    ("control", r"<archive>/lab\marginal_iq\rf34_ctrl.cs16", 34, 8e6),
    ("control", r"<archive>/lab\captures\philips_rf15_cliffclean.cs16", 15, 8e6),
    ("control", r"<archive>/lab\captures\night_rabbit_rf34_mer16.1_2227.cs16",
     34, 8e6),
]

CHUNK_SEC = 1.0
GUARD = 4 * bs.N_SYM + bs.N_A + bs.N_C + 64      # overlap between chunks


def characterize(x, fs):
    """What IS this carrier?  Cheap spectral fingerprint of a chunk.

    Distinguishes the three things a UHF TV channel can hold:
      8-VSB   : ~5.38 MHz wide, STRONG narrowband pilot ~2.69 MHz below centre
      ATSC 3.0: ~5.5-5.83 MHz wide, flat, NO pilot
      noise   : full-band, flat, no structure
    """
    n = 1 << 14
    seg = x[: (len(x) // n) * n].reshape(-1, n)
    if len(seg) > 64:
        seg = seg[:64]
    w = np.hanning(n)
    psd = np.mean(np.abs(np.fft.fftshift(np.fft.fft(seg * w, axis=1), axes=1)) ** 2,
                  axis=0)
    f = np.fft.fftshift(np.fft.fftfreq(n, 1.0 / fs))
    p = psd / psd.sum()
    c = np.cumsum(p)
    lo = f[np.searchsorted(c, 0.005)]
    hi = f[np.searchsorted(c, 0.995)]
    band = (f > lo) & (f < hi)
    pb = psd[band]
    flatness = float(np.exp(np.mean(np.log(pb + 1e-30))) / np.mean(pb))
    # pilot hunt: strongest bin vs local median, anywhere in band
    med = np.median(pb)
    i = int(np.argmax(pb))
    fb = f[band]
    # ATSC 1.0 puts a large unmodulated pilot 2.69 MHz below channel centre.
    # ATSC 3.0 has NO pilot at all.  This is a completely independent
    # discriminator from the bootstrap correlation.
    vsb = (np.abs(fb - (-2.69e6)) < 60e3)
    pilot_db = (float(10 * np.log10(pb[vsb].max() / max(med, 1e-30)))
                if vsb.any() else float("nan"))
    return {"occupied_bw_mhz": float((hi - lo) / 1e6),
            "flatness": flatness,
            "peak_db": float(10 * np.log10(pb[i] / max(med, 1e-30))),
            "peak_offset_mhz": float(fb[i] / 1e6),
            "vsb_pilot_db": pilot_db,
            "vsb_pilot": bool(pilot_db > 10.0),
            "in_band_frac": float(pb.sum() / psd.sum())}


def scan_file(path, fs_in, fmt=None, max_sec=None, pn_variant="spec",
              minor_version=0, q=bs.ZC_ROOT_MAJOR0, freq_shift=0.0,
              verbose=False, inject_db=None, deep=False):
    """Chunked scan of one capture.  Returns a summary dict.

    ``inject_db`` adds a synthesized bootstrap into every chunk at that power
    ratio relative to the capture, AFTER resampling -- i.e. through the exact
    same code path the real search uses.  This is the positive control: without
    it, a negative result is an untested instrument reporting nothing.
    """
    fmt = fmt or guess_format(path)
    total = n_samples(path, fmt)
    if max_sec:
        total = min(total, int(fs_in * max_sec))
    chunk = int(fs_in * CHUNK_SEC)
    step = chunk - int(GUARD * fs_in / bs.FS) - 1024

    best = {"peak_ratio": 0.0, "peak_ratio_max": 0.0, "score": 0.0,
            "cfo_hz": None, "position": None, "chunk": None, "t_sec": None,
            "matched_ratio": 0.0, "matched": {}, "verdict": "NO_BOOTSTRAP"}
    ratios = []
    m_ratios = []
    pos = 0
    n_chunks = 0
    fingerprint = None
    while pos + int(GUARD * fs_in / bs.FS) < total:
        raw = read_block(path, pos, min(chunk, total - pos), fmt)
        if len(raw) < GUARD * fs_in / bs.FS:
            break
        x = to_bootstrap_rate(raw, fs_in, freq_shift).astype(np.complex128)
        if fingerprint is None:
            fingerprint = characterize(raw, fs_in)
        del raw
        if inject_db is not None:
            probe, _ = bs.synthesize(fields=dict(
                ea_wake_up_1=0, ea_wake_up_2=1, min_time_to_next=10,
                system_bandwidth=0, bsr_coefficient=2, preamble_structure=170))
            amp = np.sqrt(np.mean(np.abs(x) ** 2) * 10 ** (inject_db / 10.0))
            at = 5000 + n_chunks * 733
            x[at:at + len(probe)] += (probe * amp).astype(x.dtype)
        # Run BOTH detectors on every chunk.  Ranking chunks by the structural
        # score and only matching the winner would miss a bootstrap that is
        # structurally masked by a strong co-channel carrier -- exactly the
        # case that matters here.
        r = bs.detect(x, minor_version=minor_version, q=q, pn_variant=pn_variant,
                      integer_cfo_search=12, deep=deep)
        st = r["structural"]
        mt = r.get("matched") or {}
        n_chunks += 1
        mr = mt.get("mean_peak_ratio", 0.0)
        if st.get("peak_ratio"):
            ratios.append(st["peak_ratio"])
        m_ratios.append(mr)
        if mr > best["matched_ratio"] or best["position"] is None:
            best.update(peak_ratio=st.get("peak_ratio", 0.0),
                        score=st.get("score", 0.0), cfo_hz=st.get("cfo_hz"),
                        position=st.get("position"), chunk=n_chunks - 1,
                        t_sec=pos / fs_in, matched_ratio=mr, matched=mt,
                        verdict=r["verdict"])
        best["peak_ratio_max"] = max(best["peak_ratio_max"], st.get("peak_ratio", 0.0))
        if verbose:
            print(f"      chunk {n_chunks - 1:3d} t={pos / fs_in:7.2f}s  "
                  f"struct={st.get('peak_ratio', 0):6.2f}  matched={mr:7.2f}  "
                  f"{r['verdict']}")
        pos += step

    mt = best["matched"]
    out = {"path": path, "fs_in": fs_in, "fmt": fmt, "chunks": n_chunks,
           "total_sec": total / fs_in,
           "best_peak_ratio": best["peak_ratio_max"], "best_score": best["score"],
           "best_t_sec": best["t_sec"], "cfo_hz": best["cfo_hz"],
           "median_ratio": float(np.median(ratios)) if ratios else None,
           "p95_ratio": float(np.percentile(ratios, 95)) if len(ratios) > 4 else None,
           "structural_detected": best["peak_ratio_max"] >= bs.STRUCT_RATIO_THRESHOLD,
           "fingerprint": fingerprint,
           "matched_ratio": best["matched_ratio"],
           "median_matched_ratio": float(np.median(m_ratios)) if m_ratios else None,
           "matched_detected": bool(mt.get("detected")),
           "timing_residual": mt.get("timing_residual"),
           "matched_cfo_hz": mt.get("cfo_hz"),
           "fields": mt.get("fields", {}),
           "abs_shifts": mt.get("abs_shifts")}

    # Verdict semantics.  NOTE: "matched only" is NOT a bug case on real air.
    # The structural metric is normalized by TOTAL energy, so a strong
    # co-channel carrier suppresses it; the matched filter has 2048x4 coherent
    # gain and cuts through.  Measured on this fleet: a bootstrap injected at
    # -6 dB under a real 8-VSB carrier gives matched ~21 / structural ~3.5.
    out["verdict"] = ("ATSC3_CONFIRMED" if out["structural_detected"] and out["matched_detected"]
                      else "ATSC3_STRUCTURE_ONLY" if out["structural_detected"]
                      else "ATSC3_MATCHED_ONLY" if out["matched_detected"]
                      else "NO_BOOTSTRAP")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="*")
    ap.add_argument("--fleet", action="store_true")
    ap.add_argument("--rate", type=float, default=8e6)
    ap.add_argument("--format", default=None)
    ap.add_argument("--max-sec", type=float, default=None)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--pn-sweep", action="store_true",
                    help="try every PN wiring variant (ASSUMPTION A1 arbitration)")
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--deep", action="store_true",
                    help="also run sliding matched-filter acquisition")
    ap.add_argument("--inject-db", type=float, default=None,
                    help="POSITIVE CONTROL: add a synthetic bootstrap at this "
                         "dB relative to the capture, through the same path")
    ap.add_argument("--json", default=None)
    a = ap.parse_args()

    jobs = []
    if a.fleet:
        for role, path, rf, rate in FLEET:
            jobs.append((role, path, rate, RF_HZ.get(rf), rf))
    for f in a.files:
        sc = read_sidecar(f)
        jobs.append(("target", f, sc.get("sample_rate", a.rate),
                     sc.get("center_hz"), sc.get("rf")))

    max_sec = 4.0 if a.quick else a.max_sec
    print("=" * 100)
    print("ATSC 3.0 BOOTSTRAP -- ARCHIVED CAPTURE SWEEP  (no radio; files only)")
    print(f"structural threshold peak_ratio >= {bs.STRUCT_RATIO_THRESHOLD}, "
          f"matched threshold >= 8.0")
    print("=" * 100)

    results = []
    for role, path, rate, center, rf in jobs:
        if not os.path.exists(path):
            print(f"\n[{role:7s}] {path}\n          MISSING -- skipped")
            continue
        fmt = a.format or guess_format(path)
        mb = os.path.getsize(path) / 1e6
        print(f"\n[{role:7s}] {os.path.basename(path)}   RF{rf} "
              f"{'' if center is None else f'{center/1e6:.0f} MHz'}  "
              f"{mb:.0f} MB  {fmt} @ {rate/1e6:.3f} Msps")
        print(f"          {path}")
        t0 = time.time()
        variants = list(bs.PN_VARIANTS) if a.pn_sweep else ["spec"]
        best_r = None
        for v in variants:
            r = scan_file(path, rate, fmt, max_sec, pn_variant=v,
                          verbose=a.verbose, inject_db=a.inject_db, deep=a.deep)
            r["role"], r["rf"], r["center_hz"], r["pn_variant"] = role, rf, center, v
            if a.pn_sweep:
                print(f"          PN[{v:14s}] structural {r['best_peak_ratio']:6.2f}  "
                      f"matched {r['matched_ratio']:7.2f}")
            if best_r is None or r["matched_ratio"] > best_r["matched_ratio"]:
                best_r = r
        r = best_r
        results.append(r)
        print(f"          scanned {r['total_sec']:.1f} s in {r['chunks']} chunks "
              f"({time.time()-t0:.1f}s)")
        fp = r.get("fingerprint")
        if fp:
            print(f"          spectrum: occupied BW {fp['occupied_bw_mhz']:.2f} MHz  "
                  f"flatness {fp['flatness']:.3f}  "
                  f"peak +{fp['peak_db']:.1f} dB @ {fp['peak_offset_mhz']:+.3f} MHz")
            print(f"          8-VSB pilot @ -2.690 MHz: "
                  f"{'PRESENT' if fp['vsb_pilot'] else 'absent '} "
                  f"({fp['vsb_pilot_db']:+.1f} dB over median)"
                  f"   -> {'ATSC 1.0' if fp['vsb_pilot'] else 'not ATSC 1.0'}")
        print(f"          structural peak_ratio: best {r['best_peak_ratio']:.2f}  "
              f"median {r['median_ratio']:.2f}  "
              f"(normalized corr at peak {r['best_score']:.4f})")
        print(f"          matched ratio:         {r['matched_ratio']:.2f}")
        print(f"          VERDICT: {r['verdict']}")
        if r["verdict"] == "ATSC3_CONFIRMED":
            print(f"          CFO {r['matched_cfo_hz']:+.0f} Hz   "
                  f"abs cyclic shifts {r['abs_shifts']}")
            for k, v in r["fields"].items():
                print(f"             {k:34s} {v}")

    print("\n" + "=" * 100)
    print("SUMMARY")
    print(f"{'role':8s} {'file':46s} {'struct':>8s} {'matched':>9s}  verdict")
    for r in results:
        print(f"{r['role']:8s} {os.path.basename(r['path'])[:46]:46s} "
              f"{r['best_peak_ratio']:8.2f} {r['matched_ratio']:9.2f}  {r['verdict']}")
    ctl = [r for r in results if r["role"] == "control"]
    tgt = [r for r in results if r["role"] == "target"]
    print(f"\ncontrols clean (known 8-VSB must be negative): "
          f"{sum(1 for r in ctl if r['verdict'] == 'NO_BOOTSTRAP')}/{len(ctl)}")
    print(f"targets positive: "
          f"{sum(1 for r in tgt if r['verdict'] != 'NO_BOOTSTRAP')}/{len(tgt)}")
    print("=" * 100)

    if a.json:
        with open(a.json, "w", encoding="utf-8") as fh:
            json.dump(results, fh, indent=2, default=str)
        print(f"wrote {a.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

