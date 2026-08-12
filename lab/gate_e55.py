#!/usr/bin/env python3
"""gate_e55.py -- the FFT window notch's referees.

The E55 notch (WindowNotch in m11_stream.py) replaces the fast path's
time-domain stream FIR with overlap-save FFT filtering of the decode
window only.  It is NOT bit-identical to the FIR (FFT re-association,
~4e-15), so per the E53 precedent it is DECISION-GATED: the FIR remains
the reference (exact mode always; fast mode via ATSC3_NOTCH_FFT=0), and
these legs hold the fast path to the decisions that matter:

  1. RESPONSE + VALUES: the window path's impulse response equals the FIR
     taps and its frequency response equals the FIR's across the band
     (derived FROM the taps -- fft(nb, N) -- so this is a check that the
     plumbing preserves it, not a designed approximation).  On random
     noise chopped at awkward offsets, every window sample equals the
     stream-FIR + trim reference within FFT rounding, including the
     start-of-lock zero-history transient.
  2. DECODE EQUIVALENCE, with a negative control: real gate-capture air
     plus an injected out-of-band interferer (band noise in the notch's
     stopband).  FIR-notched and FFT-notched runs of the m11 streaming
     chain must converge equivalent FEC on the same input; the un-notched
     run must do measurably worse (else the interferer proves nothing).
     MEASURED while building this leg (a finding, recorded): OOB-confined
     noise at the Ubuntu rig's ratios (0.035-0.20) does NOT touch the
     linear decode -- 444/444 with the notch OFF at ratios up to 8.  The
     rectangular-window leakage of the symbol FFT is ~1/(pi^2 d^2) per
     bin, and the nearest data carrier sits 52+ bins from the stopband.
     The un-notched path only starts losing blocks near ratio ~64 (the
     interferer 18 dB ABOVE the signal), so THAT is the gate's setting --
     and it says the rig's 0%-FEC episodes cannot be explained by
     OOB-confined interference alone (in-band contamination or front-end
     nonlinearity must be involved; the OOB ratio is a symptom metric).
  3. BOOKKEEPING (the capture-integrity law): both chains fed identical
     blocks mint the same Frames -- same count, same absolute window
     origins, every window exactly FRAME_WINDOW samples, and the window
     CONTENTS agree within tolerance (any trim/group-delay error would
     shift the FFT windows by >= 1 sample and fail this loudly).

Run: python lab/gate_e55.py [--capture lab/rate6912_rf33.cs16]
                            [--frames 10] [--oob 64] [--keep]
"""
from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
import time

import numpy as np
from scipy.signal import lfilter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m11_stream as ST                                           # noqa: E402

PY = sys.executable


# ---------------------------------------------------------------------------
# leg 0 -- the air-length bound found live: moof_info's trun count
# ---------------------------------------------------------------------------

def leg_moof_bound(verbose=True):
    """The transport wedge of live1d 09:11: a fade-corrupted trun
    sample count (n read off the air) built an n-entry list BEFORE the
    MAX_SAMPLES rejection ran.  The bound must run first; this leg holds
    it to sub-millisecond on a 2^31-1 count and to unchanged parsing on
    a sane moof."""
    import struct
    import m7_objects as O7

    def bx(t, payload):
        return struct.pack(">I", 8 + len(payload)) + t + payload

    huge = bx(b"moof", bx(b"traf", bx(
        b"trun", bytes(1) + b"\x00\x02\x00"
        + struct.pack(">I", 0x7FFFFFFF) + b"\x00" * 64)))
    huge += struct.pack(">I", 16) + b"mdat"
    t = time.perf_counter()
    mi = O7.moof_info(huge)
    dt = time.perf_counter() - t
    g1 = dt < 0.05 and mi["sizes"] == [] and mi["sizes_rejected"] == 0x7FFFFFFF
    sane = bx(b"moof", bx(b"traf", bx(
        b"trun", bytes(1) + b"\x00\x02\x00" + struct.pack(">I", 3)
        + struct.pack(">III", 10, 20, 30))))
    sane += struct.pack(">I", 76) + b"mdat"
    mi2 = O7.moof_info(sane)
    g2 = mi2["sizes"] == [10, 20, 30] and mi2["mdat_declared"] == 76
    ok = g1 and g2
    if verbose:
        print(f"    {'PASS' if g1 else 'FAIL'}  2^31-sample trun rejected "
              f"BEFORE the loop ({dt*1e3:.2f} ms, rejected="
              f"{mi['sizes_rejected']})")
        print(f"    {'PASS' if g2 else 'FAIL'}  sane trun still parses "
              f"(sizes {mi2['sizes']})")
    return ok


# ---------------------------------------------------------------------------
# leg 1 -- response and values vs the FIR reference
# ---------------------------------------------------------------------------

def leg_response(verbose=True):
    ok = True
    nb = ST.firwin(ST.NOTCH_TAPS, ST.NOTCH_CUT, fs=ST.FS_POST)
    wn = ST.WindowNotch(nb)
    m1 = ST.NOTCH_TAPS - 1

    # impulse response: history of zeros + a unit impulse + zeros
    x = np.zeros(4096, np.complex128)
    x[m1] = 1.0
    h = wn.apply(x)[:ST.NOTCH_TAPS]
    g = np.allclose(h, nb, rtol=0, atol=1e-13)
    ok &= g
    if verbose:
        print(f"    {'PASS' if g else 'FAIL'}  impulse response == FIR taps "
              f"(max |delta| {np.max(np.abs(h - nb)):.2e})")

    # frequency response across the band, on the fine grid
    Hf = np.fft.fft(h, 1 << 14)
    Hr = np.fft.fft(nb, 1 << 14)
    g = np.allclose(Hf, Hr, rtol=0, atol=1e-12)
    ok &= g
    if verbose:
        print(f"    {'PASS' if g else 'FAIL'}  frequency response == FIR "
              f"across the band (max |delta| {np.max(np.abs(Hf - Hr)):.2e})")

    # window values vs stream FIR + trim, interior windows AND the
    # start-of-lock zero-history case
    rng = np.random.default_rng(55)
    n = 1_500_000
    raw = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    ref = lfilter(nb, 1.0, raw)      # ref_aligned[t] = ref[t + half]
    half = ST.NOTCH_TAPS // 2
    W = ST.FRAME_WINDOW
    worst = 0.0
    for lo in (0, 37, half, 12345, n - W - half - 1):
        a = lo - half
        zf = max(0, -a)
        buf = raw[a + zf:lo + W + half]
        if zf:
            buf = np.concatenate((np.zeros(zf, buf.dtype), buf))
        w = wn.apply(buf)
        want = ref[lo + half:lo + W + half]
        if len(w) != W:
            ok = False
            print(f"    FAIL  window at lo={lo}: {len(w)} samples != "
                  f"FRAME_WINDOW {W}")
            continue
        worst = max(worst, float(np.max(np.abs(w - want))))
    g = worst < 1e-11
    ok &= g
    if verbose:
        print(f"    {'PASS' if g else 'FAIL'}  window values == stream FIR + "
              f"trim reference, 5 offsets incl. lock start "
              f"(max |delta| {worst:.2e}, tol 1e-11)")
    return ok


# ---------------------------------------------------------------------------
# the contaminated capture
# ---------------------------------------------------------------------------

def oob_ratio(y):
    """Exactly _decide_notch's measurement."""
    n = 1 << 18
    S = np.abs(np.fft.fft(y[:n])) ** 2
    f = np.fft.fftfreq(n, 1.0 / ST.FS_POST)
    return float(S[np.abs(f) > 3.0e6].mean()
                 / max(S[np.abs(f) < 2.7e6].mean(), 1e-30))


def make_contaminated(capture, out_path, seconds, target_ratio, seed=55,
                      verbose=True):
    n = int(seconds * ST.FS_POST)
    raw = np.fromfile(capture, dtype=np.int16, count=2 * n).astype(np.float64)
    x = raw[0::2] + 1j * raw[1::2]
    n = len(x)
    # band noise confined to the notch stopband, both sides -- the shape of
    # the Ubuntu rig's USB hash as _decide_notch sees it
    rng = np.random.default_rng(seed)
    w = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    F = np.fft.fft(w)
    f = np.fft.fftfreq(n, 1.0 / ST.FS_POST)
    F[~((np.abs(f) > 2.96e6) & (np.abs(f) < 3.40e6))] = 0.0
    nz = np.fft.ifft(F)
    # scale so the contaminated OOB/in-band ratio equals target_ratio
    probe = 1 << 18
    Sx = np.abs(np.fft.fft(x[:probe])) ** 2
    Sn = np.abs(np.fft.fft(nz[:probe])) ** 2
    fp = np.fft.fftfreq(probe, 1.0 / ST.FS_POST)
    ob, ib = np.abs(fp) > 3.0e6, np.abs(fp) < 2.7e6
    a2 = (target_ratio * Sx[ib].mean() - Sx[ob].mean()) / Sn[ob].mean()
    a = np.sqrt(max(a2, 0.0))
    y = x + a * nz
    peak = float(np.max(np.abs(np.concatenate((y.real, y.imag)))))
    if peak > 32000:
        y *= 32000.0 / peak          # keep int16 honest; report it
    out = np.empty(2 * n, np.int16)
    out[0::2] = np.round(y.real).astype(np.int16)
    out[1::2] = np.round(y.imag).astype(np.int16)
    out.tofile(out_path)
    r = oob_ratio(y)
    if verbose:
        print(f"    (interferer: stopband noise 2.96-3.40 MHz, amplitude "
              f"{a:.1f}, OOB/in-band {r:.4f} target {target_ratio}, "
              f"peak {peak:.0f}/32767{' CLIPPED' if peak > 32000 else ''})")
    return r


# ---------------------------------------------------------------------------
# leg 2 -- decode equivalence via the real streaming chain
# ---------------------------------------------------------------------------

FEC_RE = re.compile(r"FEC Blocks\s+(\d+)/(\d+)\s+converged")
NOTCH_RE = re.compile(r"notch\s+([0-9.]+)")


def run_watch(capture, frames, env_extra, threads=6):
    env = dict(os.environ)
    env.update(env_extra)
    r = subprocess.run(
        [PY, os.path.join(HERE, "m11_watch.py"), "--capture", capture,
         "--rate", str(ST.FS_POST), "--player", "none", "--accel", "cpu",
         "--threads", str(threads), "--report", "1e9",
         "--frames", str(frames)],
        capture_output=True, text=True, env=env)
    conv = tot = -1
    notch_ms = None
    m = FEC_RE.search(r.stdout)
    if m:
        conv, tot = int(m.group(1)), int(m.group(2))
    for line in r.stdout.splitlines():
        if "front end" in line and "notch" in line:
            m = NOTCH_RE.search(line)
            if m:
                notch_ms = float(m.group(1))
    return conv, tot, notch_ms, r.stdout


def leg_decode(cap_dirty, frames, verbose=True):
    ok = True
    runs = {}
    # E60 note: this leg tests the NOTCH in isolation, so the margin levers
    # are pinned OFF -- with them ON the negative control went vacuous
    # (measured 8/09: the smoothed CE recovered the un-notched interferer
    # damage entirely, OFF 370/370 == notched 370/370).  That recovery is
    # E60's finding; the notch decision machinery is this gate's subject.
    for name, env in (("FIR", {"ATSC3_NOTCH_FFT": "0", "ATSC3_MARGIN": "0"}),
                      ("FFT", {"ATSC3_NOTCH_FFT": "1", "ATSC3_MARGIN": "0"}),
                      ("OFF", {"ATSC3_NOTCH": "0", "ATSC3_MARGIN": "0"})):
        t = time.time()
        conv, tot, notch_ms, out = run_watch(cap_dirty, frames, env)
        runs[name] = (conv, tot, notch_ms)
        if verbose:
            print(f"    {name}: FEC {conv}/{tot} converged, notch stage "
                  f"{notch_ms} ms/Frame  ({time.time()-t:.0f} s)")
        if conv < 0:
            print(f"    FAIL  {name} run produced no FEC summary\n"
                  f"{out[-1500:]}")
            ok = False
    if not ok:
        return ok, runs
    (cf, tf, _), (cx, tx, nx), (co, to, _) = \
        runs["FIR"], runs["FFT"], runs["OFF"]
    # equivalence: within noise of each other on the same input
    tol = max(2, int(0.02 * max(tf, 1)))
    g = tf == tx and abs(cf - cx) <= tol
    ok &= g
    if verbose:
        print(f"    {'PASS' if g else 'FAIL'}  FFT-notched FEC convergence "
              f"equivalent to FIR-notched (|{cf}-{cx}| <= {tol})")
    # the negative control must control: measurably worse, beyond the
    # equivalence tolerance (at --oob 64 measured: OFF 385 vs notched 444)
    g = co <= cf - max(10, tol)
    ok &= g
    if verbose:
        print(f"    {'PASS' if g else 'FAIL'}  un-notched degrades on this "
              f"interferer ({co} vs {cf} converged) -- the notch is doing "
              f"real work")
    return ok, runs


# ---------------------------------------------------------------------------
# leg 3 -- bookkeeping: both chains, identical blocks, same windows
# ---------------------------------------------------------------------------

def leg_bookkeeping(cap_dirty, frames, verbose=True):
    from concurrent.futures import ThreadPoolExecutor
    import m2_pilots as MP
    ok = True
    ex = ThreadPoolExecutor(max_workers=6)
    fes = {}
    wins = {"FIR": [], "FFT": []}
    os.environ["ATSC3_NOTCH_FFT"] = "0"
    fes["FIR"] = ST.FrontEnd(ST.FS_POST, ex=ex, fast=True)
    os.environ["ATSC3_NOTCH_FFT"] = "1"
    fes["FFT"] = ST.FrontEnd(ST.FS_POST, ex=ex, fast=True)
    os.environ.pop("ATSC3_NOTCH_FFT", None)
    assert fes["FFT"].fft_notch and not fes["FIR"].fft_notch
    fmt = MP.guess_format(cap_dirty)
    block = 219_996                  # deliberately awkward block size
    pos = 0
    while min(len(wins["FIR"]), len(wins["FFT"])) < frames:
        x = MP.read_block(cap_dirty, pos, block, fmt)
        if not len(x):
            break
        pos += len(x)
        for k, fe in fes.items():
            fe.push(x)
            for fr in fe.frames():
                if len(wins[k]) < frames:
                    wins[k].append(fr)
    nf = min(len(wins["FIR"]), len(wins["FFT"]))
    g = nf >= frames
    ok &= g
    if verbose:
        print(f"    {'PASS' if g else 'FAIL'}  both chains minted "
              f"{nf} Frames from identical blocks (asked {frames})")
    worst = 0.0
    for i in range(nf):
        ia, wa, ta, _ = wins["FIR"][i]
        ib, wb, tb, _ = wins["FFT"][i]
        # same Frame index; window origin in AIR samples: the FIR chain's
        # tape is trimmed by half a filter (128) relative to the raw tape,
        # so equal air positions mean equal (org-corrected) coordinates.
        same_idx = ia == ib
        same_len = len(wa) == ST.FRAME_WINDOW and len(wb) == ST.FRAME_WINDOW
        d = float(np.max(np.abs(wa - wb))) if same_len else np.inf
        worst = max(worst, d)
        ok &= same_idx and same_len
        if not (same_idx and same_len):
            print(f"    FAIL  frame {i}: idx {ia}/{ib}, "
                  f"len {len(wa)}/{len(wb)}")
    scale = float(np.mean(np.abs(wins["FIR"][0][1]))) if nf else 1.0
    # FFT re-association on int16-scaled air lands ~1e-11 absolute; a
    # one-sample trim error would land at the signal's own magnitude
    g = worst < 1e-6 * max(scale, 1.0)
    ok &= g
    if verbose:
        print(f"    {'PASS' if g else 'FAIL'}  window contents agree, "
              f"{nf} Frames x {ST.FRAME_WINDOW} samples "
              f"(max |delta| {worst:.2e} on mean |x| {scale:.1f}) -- any "
              f"trim/group-delay error would shift a window and fail this")
    n_in = fes["FIR"].n_der
    g = fes["FFT"].n_der == n_in
    ok &= g
    if verbose:
        print(f"    {'PASS' if g else 'FAIL'}  identical sample intake "
              f"({n_in} derotated both chains)")
    ex.shutdown(wait=False)
    return ok


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--capture", default="rate6912_rf33.cs16")
    ap.add_argument("--frames", type=int, default=10)
    ap.add_argument("--oob", type=float, default=64.0,
                    help="target OOB/in-band ratio for the interferer "
                         "(64 = 18 dB over the signal; the un-notched "
                         "linear decode shrugs off anything much less)")
    ap.add_argument("--keep", action="store_true")
    a = ap.parse_args(argv)
    cap = a.capture if os.path.isabs(a.capture) else os.path.join(HERE,
                                                                  a.capture)
    print("E55 -- FFT window notch gates (response / decode / bookkeeping)")
    print("=" * 72)
    ok = True

    print("\n  0. air-length bound: moof_info trun count (the live1d wedge)")
    ok &= leg_moof_bound()

    print("\n  1. response and values vs the FIR reference")
    ok &= leg_response()

    print(f"\n  2. decode equivalence on real air + injected OOB interferer")
    cap_dirty = os.path.join(HERE, "e55_oob.cs16")
    seconds = 0.4 + (a.frames + 2) * ST.FRAME_SEC
    make_contaminated(cap, cap_dirty, seconds, a.oob)
    try:
        g, _ = leg_decode(cap_dirty, a.frames)
        ok &= g

        print("\n  3. bookkeeping: identical blocks, identical windows")
        ok &= leg_bookkeeping(cap_dirty, min(a.frames, 6))
    finally:
        if not a.keep and os.path.exists(cap_dirty):
            os.remove(cap_dirty)

    print(f"\n  {'ALL E55 GATES PASS' if ok else 'E55 GATE FAILURES'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
