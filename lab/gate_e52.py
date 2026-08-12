#!/usr/bin/env python3
"""gate_e52.py -- the CPU real-time campaign's own referees, on real air.

Four legs, each against a re-derived reference (never a stored expectation):

  1. `min_sum_decode_batch` float64 vs `min_sum_decode`, per FEC Block of
     real Frames, at every ladder alpha: bits, converged, iters and
     unsatisfied-count must be EQUAL -- the batching is a re-shaping, not a
     re-rounding, and this leg proves it bit for bit.
  2. `demod_coh` vs `demod_data(...)[1]` over the whole fine-timing candidate
     span: the coherence float must be IDENTICAL (the code is copied lines).
  3. FrameDecoder(cpu_fast=True) vs FrameDecoder(cpu_fast=False) on the same
     Frames: every Baseband Packet byte-identical, PLP16 bits identical,
     convergence identical.  This is the frame-level referee for the float32
     LDPC + score-demap + BLAS-CPE re-associations ON THIS CAPTURE; the m11
     e2e gate then covers the whole capture at datagram level.
  4. float32 batch LDPC informational: convergence pattern and hard-output
     agreement vs float64 (reported; leg 3 is the pass/fail authority).

Run: python lab/gate_e52.py [--capture lab/rate6912_rf33.cs16] [--frames 3]
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m3_ldpc as LD                                              # noqa: E402
import m6_cells as C                                              # noqa: E402
import m9_accel                                                   # noqa: E402
import m9_fast as F9                                              # noqa: E402


def build_lams(fd, y, t0):
    """The lam arrays exactly as decode_frame builds them (same calls)."""
    import m6_tbi as T
    pool, info = fd.cell_pool(y, t0, {})
    key = len(pool)
    if key not in fd.regions:
        fd.regions[key] = C.constellation_regions(key)
    reg, dv = fd.regions[key]
    pool, _ = fd.cpe_correct(pool, info["symbol_of"], [fd.pts0, fd.pts16],
                             reg, dv)
    c16 = pool[C.PLP16["start"]:C.PLP16["start"] + C.PLP16["size"]]
    q16 = fd.demap_batch(c16[None, :], fd.pts16, fd.ones16)
    lam16 = np.empty((1, 16200))
    lam16[:, fd.ch16.lam_of_q] = q16
    plp0 = pool[:C.PLP0["size"]]
    nrows = fd.ch0.cells_per_fec
    ncols = C.PLP0["n_fec_max"] // C.PLP0["nti"]
    per_ti = len(plp0) // C.PLP0["nti"]
    blocks = []
    for ti in range(C.PLP0["nti"]):
        seg = plp0[ti * per_ti:(ti + 1) * per_ti]
        for j in range(ncols):
            blocks.append(T.fec_block(seg, j, nrows, ncols, 0))
    cells = np.array(blocks)
    q = fd.demap_batch(cells, fd.pts0, fd.ones0)
    lam = np.empty((len(cells), 16200))
    lam[:, fd.ch0.lam_of_q] = q
    return lam, lam16


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--capture", default="rate6912_rf33.cs16")
    ap.add_argument("--rate", type=float, default=6.912e6)
    ap.add_argument("--frames", type=int, default=3)
    ap.add_argument("--threads", type=int, default=12)
    a = ap.parse_args(argv)
    cap = a.capture if os.path.isabs(a.capture) else os.path.join(HERE,
                                                                  a.capture)
    print("E52 -- CPU fast-path gates (batch LDPC / demod_coh / fast frame)")
    print("=" * 72)
    ok = True

    span = 0.05 + (a.frames + 1) * F9.FRAME_SEC
    fd = F9.FrameDecoder(threads=a.threads, backend="cpu", cpu_fast=False)
    y = F9.load_fast(cap, a.rate, span_sec=span * (a.rate / 6.912e6),
                     ex=fd.ex)[0]

    # -- leg 2: demod_coh scalar bit-equal; batched same DECISION ----------
    centre = C.BOOTSTRAP
    cands = list(range(centre - 20, centre + 21))
    batch = m9_accel.demod_coh_batch(y, cands, 1, False)
    bad = 0
    refs = []
    for i, t in enumerate(cands):
        ref = m9_accel.demod_data(y, t, 1, False)[1]
        refs.append(ref)
        bad += int(ref != m9_accel.demod_coh(y, t, 1, False))
    g = bad == 0
    ok &= g
    print(f"    {'PASS' if g else 'FAIL'}  demod_coh == demod_data coherence, "
          f"41 candidates, bit-equal floats ({bad} differ)")
    # pocketfft's multi-row path rounds differently in the last ulps, so the
    # batched form is held to the DECISION, not the bits: same argmax, and
    # the drift measured.  (Frame-level bytes are leg 3's job.)
    refs = np.array(refs)
    drift = float(np.max(np.abs(batch - refs) / np.maximum(refs, 1e-30)))
    same_pick = int(np.argmax(batch)) == int(np.argmax(refs))
    ok &= same_pick
    print(f"    {'PASS' if same_pick else 'FAIL'}  batched fine timing picks "
          f"the same t0 (max coherence rel-drift {drift:.2e})")

    # frames + lams off the exact path
    t_prev = None
    lams = []
    t0s = []
    for fi in range(a.frames):
        centre = C.BOOTSTRAP if t_prev is None else t_prev + C.FRAME_SAMPLES
        t0, coh = fd.fine_timing(y, span=20, centre=centre)
        t_prev = t0
        t0s.append(t0)
        lams.append(build_lams(fd, y, t0))

    # -- leg 1: batch float64 == scalar, every ladder alpha ----------------
    n_diff = 0
    n_rows = 0
    for lam, lam16 in lams:
        for l2d, chain in ((lam, fd.ch0), (lam16, fd.ch16)):
            for alpha in fd.ladder_for(chain):
                bb, cb, ib, db = LD.min_sum_decode_batch(
                    l2d, chain.checks, chain.ninner, iters=fd.iters,
                    alpha=alpha, dtype=np.float64)
                for i in range(len(l2d)):
                    b, c, it, d = LD.min_sum_decode(
                        l2d[i], chain.checks, chain.ninner, iters=fd.iters,
                        alpha=alpha)
                    n_rows += 1
                    if not (np.array_equal(b, bb[i]) and c == cb[i]
                            and it == ib[i] and d == db[i]):
                        n_diff += 1
    g = n_diff == 0
    ok &= g
    print(f"    {'PASS' if g else 'FAIL'}  batch float64 LDPC == scalar, "
          f"{n_rows} (block, alpha) decodes bit-identical "
          f"({n_diff} differ)")

    # -- leg 4: float32, informational -------------------------------------
    agree = tot = conv_same = 0
    for lam, lam16 in lams:
        for l2d, chain in ((lam, fd.ch0), (lam16, fd.ch16)):
            alpha = fd.ladder_for(chain)[0]
            b64, c64, _, _ = LD.min_sum_decode_batch(
                l2d, chain.checks, chain.ninner, iters=fd.iters,
                alpha=alpha, dtype=np.float64)
            b32, c32, _, _ = LD.min_sum_decode_batch(
                l2d, chain.checks, chain.ninner, iters=fd.iters,
                alpha=alpha, dtype=np.float32)
            tot += len(l2d)
            conv_same += int(np.array_equal(c64, c32))
            for i in range(len(l2d)):
                if c64[i] and c32[i] and np.array_equal(b64[i], b32[i]):
                    agree += 1
    print(f"    INFO  float32 LDPC: {agree}/{tot} converged blocks "
          f"bit-identical to float64 (frame-level referee is leg 3)")

    # -- leg 3: the fast frame vs the exact frame --------------------------
    # (including the fast chain's own batched fine-timing TRACK: same t0
    # decisions, frame after frame, as the exact scalar track)
    fdf = F9.FrameDecoder(threads=a.threads, backend="cpu", cpu_fast=True)
    t_prev = None
    for fi in range(a.frames):
        centre = C.BOOTSTRAP if t_prev is None else t_prev + C.FRAME_SAMPLES
        t0f, _ = fdf.fine_timing(y, span=20, centre=centre)
        t_prev = t0f
        g = t0f == t0s[fi]
        ok &= g
        if not g:
            print(f"    FAIL  frame {fi}: batched fine timing t0 {t0f} != "
                  f"exact {t0s[fi]}")
    print(f"    PASS  fast fine-timing track: {a.frames} frames, t0 "
          f"decisions identical" if ok else
          "    (fine-timing decisions diverged)")
    # -- E53: the timing-line track.  Its t0s NEED NOT equal the full
    # scan's (the coherence plateau is flat; the full scan itself jitters
    # +-20 samples frame to frame chasing metric noise).  What must hold:
    # every tracked t0 stays within the +-FT_SPAN plateau of the full
    # scan's, and the DECODED BYTES at tracked t0s equal the exact chain's
    # bytes at full-scan t0s.  The plateau is the physics; the bytes are
    # the referee.
    import m11_stream as STm
    ft = STm.FineTrack(fast=True)
    centre = C.BOOTSTRAP
    track_t0s = []
    for fi in range(a.frames):
        t0t, coh = ft.search(y, centre, fd.ex)
        track_t0s.append(t0t)
        centre = t0t + C.FRAME_SAMPLES + ft.drift()
    # Bound: twice the scan's own search span.  The full scan's t0 jitters
    # +-FT_SPAN on the metric plateau (its noise floor), and the tracker
    # may add its own scan-span of dead-reckon error between anchors --
    # beyond that it is walking, not tracking.  MEASURED on this capture:
    # bytes survive ~86 samples of deviation and FEC collapses somewhere
    # past ~90 (the GI/8 bound of the first build let exactly that
    # through).  The BYTES below remain the real referee.
    dev = [abs(tr - fu) for tr, fu in zip(track_t0s, t0s)]
    # 3*FT_SPAN = 60: the anchor scans themselves wander +-FT_SPAN on the
    # plateau and the reckon may add as much again; the measured byte
    # cliff on this capture is past ~90 (86 decoded clean, ~100+ fell to
    # 0/74), so 60 sits a third below the cliff and still trips the
    # FULL_EVERY=8 walk-off (112) that motivated it.
    bound = 3 * STm.FT_SPAN
    g = max(dev) <= bound
    ok &= g
    print(f"    {'PASS' if g else 'FAIL'}  E53 timing line: {a.frames} "
          f"frames, max |track - full scan| = {max(dev)} samples "
          f"(bound 3*FT_SPAN = {bound}), {dict(ft.stats)}")

    n_pkt_diff = n_frames = 0
    for fi, t0 in enumerate(t0s):
        r16a, pka, dga = fd.decode_frame(y, t0)
        # the fast decoder decodes at the TRACKED t0 -- the exact chain at
        # the full-scan t0 -- and the bytes must still agree: that is the
        # whole E53 timing claim in one line.
        r16b, pkb, dgb = fdf.decode_frame(y, track_t0s[fi])
        same_p = len(pka) == len(pkb) and all(x == z for x, z in
                                              zip(pka, pkb))
        same_16 = (bool(r16a["converged"]) == bool(r16b["converged"])
                   and np.array_equal(r16a["bits"], r16b["bits"]))
        same_d = (dga["converged"] == dgb["converged"]
                  and dga["bch_ok"] == dgb["bch_ok"])
        good = same_p and same_16 and same_d
        ok &= good
        n_frames += 1
        n_pkt_diff += sum(1 for x, z in zip(pka, pkb) if x != z)
        print(f"    {'PASS' if good else 'FAIL'}  frame {fi}: fast vs exact "
              f"-- {dgb['converged']}/74 conv, {dgb['bch_ok']}/74 BCH0, "
              f"packets {'BYTE-IDENTICAL' if same_p else 'DIFFER'}, "
              f"PLP16 {'identical' if same_16 else 'DIFFERS'}")
    fdf.close()
    fd.close()
    print(f"    {n_frames} frames, {n_pkt_diff} differing packets")
    print(f"\n  {'ALL E52 GATES PASS' if ok else 'E52 GATE FAILURES'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
