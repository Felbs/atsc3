#!/usr/bin/env python3
"""gate_e82.py -- the LDM/CTI live path, gated, with controls that must fail.

Every leg here answers a question of the form "is the new stage the same as
the stage it replaces" or "does the thing that should be impossible stay
impossible".  Nothing compares against a stored expectation: the RF33
reference is RE-RUN as a subprocess, and the batch de-interleaver, the scalar
BICM chain and the untouched `m6_cells.cell_pool_g` are called live.

    0  L1-Basic 6.5.2.7 repetition: the pre-E82 reading must FAIL on RF25 and
       must be BIT-IDENTICAL on a Mode with no repetition (RF33's Mode 3)
    1  FrameDemod.cell_pool == m6_cells.cell_pool_g, cell for cell
    2  FrameDemod.cpe       == m10_core.cpe, cell for cell
    3  CtiStream (streaming, trimmed, zero-filling) == m10_cti.deinterleave
       (whole-array), cell for cell, for every Block it emits
    4  the batched BICM+LDPC+BCH == m6_bicm.PlpChain scalar, BYTE for BYTE
    5  CONTROLS: wrong commutator phase, wrong Nrows, CTI bypassed -> 0
    6  RF33 through m11_watch: the datagram file SHA-256 is unchanged
    7  phase acquisition with NO hand-set start Frame, from N start offsets
    8  the commutator continuity identity (m10_core P4) on live L1s
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m2_pilots as MP                                             # noqa: E402
import m3_l1basic as L1B                                           # noqa: E402
import m4_l1detail as D4                                           # noqa: E402
import m6_bicm as B                                                # noqa: E402
import m6_cells as C                                               # noqa: E402
import m10_core as M10                                             # noqa: E402
import m10_cti as CTI                                              # noqa: E402
import m3_freqint as FI                                            # noqa: E402
import m44_ldm as M44                                              # noqa: E402
from atsc3 import bootstrap as bs                                  # noqa: E402

PY = sys.executable
FOX = os.path.join(HERE, "..", "data", "fox25",
                   "rf25_fox_g4_0810_2254.cs16")
RF33 = os.path.join(HERE, "rate6912_rf33.cs16")
RATE = 6.912e6

RESULTS = []


def leg(name, ok, detail=""):
    RESULTS.append((name, bool(ok), detail))
    print(f"    {'PASS' if ok else '*** FAIL ***'}  {name}"
          + (f"\n             {detail}" if detail else ""))
    return ok


def load_frames(path, n, start_sec=0.0):
    """n consecutive Frame windows, de-rotated, off the real bootstrap grid."""
    fmt = MP.guess_format(path)
    head = MP.read_block(path, int(start_sec * RATE), int(0.30 * RATE), fmt)
    h = MP.find_bootstraps(MP.resample_to(head, RATE, bs.FS))[0]
    f0 = int(start_sec * RATE) + int(round(h["position"] * RATE / bs.FS))
    ps = h["fields"].get("preamble_structure")
    # measure the Frame grid from the air rather than assuming the nominal
    # length: the sample clocks differ by ~0.8 ppm and the drift is real
    probe = MP.read_block(path, f0 + 1_600_000, int(0.35 * RATE), fmt)
    h2 = MP.find_bootstraps(MP.resample_to(probe, RATE, bs.FS))[0]
    step = 1_600_000 + int(round(h2["position"] * RATE / bs.FS))
    x = MP.read_block(path, f0, step * n + 500_000, fmt)
    y = x.astype(np.complex128)
    y *= np.exp(-2j * np.pi * h["fine_cfo_hz"] * np.arange(len(y)) / RATE)
    return y, step, ps


# ---------------------------------------------------------------------------
def gate_l1basic_repetition():
    print("\n  [0] A/322 6.5.2.7 parity repetition for L1-Basic Mode 1")
    y, step, ps = load_frames(FOX, 6, 0.0)
    ok_true = ok_false = 0
    for k in range(6):
        w = y[k * step:k * step + 400_000]
        t0, _c = M44.plateau_t0(w, M44.BOOTSTRAP, ps)
        for rep, box in ((True, "t"), (False, "f")):
            L = M44.l1_from_window(w, t0, ps=ps) if rep else None
            if not rep:
                # the control has to run the SAME path with only the reading
                # changed, so drive l1basic_of directly
                import m3_spec as S3
                from m3_preamble import channel_estimate, preamble_geometry
                Pt = S3.PREAMBLE_STRUCTURE[ps]
                nfft, gi, dx, mode = Pt["fft"], Pt["gi"], Pt["dx"], Pt["l1b_mode"]
                Y = np.fft.fftshift(np.fft.fft(w[t0 + gi:t0 + gi + nfft]))
                lo, n, pilot, cp, data, _o = preamble_geometry(nfft, gi, dx, 4)
                H, _hp = channel_estimate(Y, nfft, gi, dx, 4, 0)
                from m3_preamble import fft_bins
                bb = fft_bins(nfft, lo + data)
                m = (bb >= 0) & (bb < nfft)
                z0 = Y[bb[m]] / H[data[m]]
                x0 = FI.deinterleave(z0, nfft, 0, direction="forward",
                                     toggle="i")
                f = D4.l1basic_of(x0, mode, repeat=False)
                ok_false += f is not None
            else:
                ok_true += bool(L and L.get("ok"))
    leg("L1-Basic Mode 1 verifies WITH 6.5.2.7 repetition",
        ok_true == 6, f"{ok_true}/6 Frames")
    leg("CONTROL: the pre-E82 reading (no repetition) FAILS",
        ok_false == 0,
        f"{ok_false}/6 Frames verified without it "
        f"(E76 saw ~1 in 9 slip through; a control that sometimes passes is "
        f"still a control that discriminates 6:0 here)")

    # Modes with Nrepeat = 0 must be untouched, bit for bit
    rng = np.random.default_rng(3)
    same = True
    for mode in (2, 3, 4, 5, 6, 7):
        nfec, cells, _np = __import__("m3_spec").l1_basic_lengths(mode)
        bl = rng.standard_normal(nfec)
        a = L1B.build_frame_llr(bl)
        b = L1B.build_frame_llr_rep(bl, mode)
        same &= np.array_equal(a, b)
    leg("Modes 2..7 (Nrepeat = 0, incl. RF33's Mode 3) BIT-IDENTICAL",
        same, "build_frame_llr vs build_frame_llr_rep, np.array_equal")
    nrep = [L1B.repetition_bits(m) for m in range(1, 8)]
    leg("only Mode 1 carries a repetition block", nrep[0] == 3672
        and all(v == 0 for v in nrep[1:]), f"Nrepeat per Mode = {nrep}")


# ---------------------------------------------------------------------------
def _plan_and_demod(y, step, ps, k=0):
    w = y[k * step:k * step + 1_700_000]
    t0, _c = M44.plateau_t0(w, M44.BOOTSTRAP, ps)
    L = M44.l1_from_window(w, t0, ps=ps)
    assert L and L.get("ok"), "L1 did not verify -- gate cannot run"
    plan = M44.LdmPlan.from_l1_result(L, "gate")
    fd = M44.FrameDemod(plan, threads=4)
    return plan, fd, w, t0, L


def gate_demod(nframes=8):
    print("\n  [1][2] the Frame demodulator against the untouched reference")
    y, step, ps = load_frames(FOX, nframes, 0.0)
    plan, fd, w, t0, _L = _plan_and_demod(y, step, ps)
    rep = {}
    # E86c: these two legs assert that m44's RESTRUCTURED demodulator (symbols
    # grouped by pilot class, one batched FFT) is faithful to the reference --
    # a question about the restructuring, not about the margin levers.  CE
    # smoothing deliberately produces DIFFERENT, better cells, so it is forced
    # off here.  Silencing the legs instead would have thrown away the only
    # check that the fast path is honest.  The lever gets its own leg below,
    # and between them the gate is now STRONGER than before: the fast path
    # must match the reference, AND the lever must actually change something.
    _ce_w, M44.CE_W = M44.CE_W, 0
    try:
        poolA, ownerA = fd.cell_pool(w, t0, rep)
    finally:
        M44.CE_W = _ce_w
    poolB, infoB = C.cell_pool_g(w, t0, plan.g)
    d = float(np.abs(poolA - poolB).max())
    leg("cell_pool == m6_cells.cell_pool_g (cell for cell), CE off",
        d < 1e-12 and np.array_equal(ownerA, infoB["symbol_of"])
        and len(poolA) == plan.pool_pred,
        f"max|delta| {d:.3e} over {len(poolA)} cells "
        f"(interp re-association only); symbol_of identical; "
        f"pool {len(poolA)} == predicted {plan.pool_pred}")
    # the lever's own negative control: a smoother that changed nothing would
    # pass every equivalence leg and quietly buy nothing
    rep_sm = {}
    fd_sm = M44.FrameDemod(plan, threads=4)
    poolS, _ownerS = fd_sm.cell_pool(w, t0, rep_sm)
    dsm = float(np.abs(poolS - poolB).max())
    nsm = rep_sm.get("ce_smoothed_syms", 0)
    leg("CE smoothing is NOT a no-op: it moves the cells and reports how many "
        "symbols it touched",
        M44.CE_W > 0 and nsm > 0 and dsm > 1e-9
        and len(poolS) == plan.pool_pred,
        f"CE_W={M44.CE_W}, smoothed {nsm} symbols "
        f"(fallback {rep_sm.get('ce_fallback_syms')}), "
        f"max|delta| vs reference {dsm:.3e}")
    cA = fd.cpe(poolA, ownerA)          # poolA is the CE-off pool, as above
    cB = M10.cpe(poolB, infoB["symbol_of"], plan.g, plan.core,
                 dummy_start=plan.core_total)
    d2 = float(np.abs(cA - cB).max())
    leg("cpe == m10_core.cpe (cell for cell), CE off", d2 < 1e-12,
        f"max|delta| {d2:.3e}")
    ag = fd.dummy_agreement(cA)
    leg("A/322 7.2.6.5 dummy-cell referee agrees with the scrambler",
        ag is not None and ag > 0.99, f"sign agreement {ag}")
    # cpe_fast differs from cpe by ESTIMATOR NOISE, not by algorithm, so the
    # meaningful bound is on the physical quantity it estimates: a per-symbol
    # complex gain.  Decimating 4:1 doubles that estimator's sigma; what has
    # to be shown is that the residual stays far inside the decision margin
    # of the constellation being demapped (QPSK: 45 degrees).  The BYTES are
    # gated separately and are the real referee -- leg 2b.
    cF = fd.cpe_fast(poolA, ownerA)
    ratio = cF / np.where(np.abs(cA) > 1e-12, cA, 1.0)
    dg = float(np.abs(np.abs(ratio) - 1.0).max())
    dp = float(np.abs(np.angle(ratio)).max())
    leg("cpe_fast's per-symbol gain stays far inside the QPSK margin",
        dg < 0.10 and dp < 0.15,
        f"max |gain-1| {dg:.4f}, max phase {dp:.4f} rad "
        f"({np.degrees(dp):.2f} deg vs a 45 deg decision margin) -- "
        f"estimator noise from the 4:1 decimation, not an algorithm change")
    return plan, fd, y, step, ps


# ---------------------------------------------------------------------------
def gate_cti_and_fec(plan, fd, y, step, ps, nframes=4):
    print("\n  [3][4] the streaming CTI and the batched FEC, vs the batch path")
    # per-Frame cells + the L1 of the FIRST Frame (which is the whole point)
    cells = []
    for k in range(nframes):
        w = y[k * step:k * step + plan.frame_window + 4096]
        t0, _c = M44.plateau_t0(w, M44.BOOTSTRAP, ps)
        cells.append(fd.core_cells(w, t0))
    w0 = y[:plan.frame_window + 4096]
    t00, _c = M44.plateau_t0(w0, M44.BOOTSTRAP, ps)
    L0 = M44.l1_from_window(w0, t00, ps=ps)
    d = M44.start_row_of(L0, plan.core["id"])
    sr, fbs = int(d["start_row"]), int(d["fec_block_start"])

    # -- 2b: the fast CPE must not change a single decoded BYTE
    cells_exact = []
    for k in range(nframes):
        w = y[k * step:k * step + plan.frame_window + 4096]
        tk, _c = M44.plateau_t0(w, M44.BOOTSTRAP, ps)
        cells_exact.append(fd.core_cells(w, tk, fast=False))

    # -- streaming
    st = M44.CtiStream(plan)
    st.reset(0, sr, fbs)
    got = []
    for k in range(nframes):
        st.push(k, cells[k])
        got += st.take_blocks()
    # -- batch (m10_cti, whole array)
    recv = np.concatenate(cells)
    outb, valid = CTI.deinterleave(recv, plan.nrows, sr)
    C0 = CTI.solve_C(fbs, sr, plan.nrows)
    same = True
    for bid, seg in got:
        lo = C0 + bid * plan.ncell
        same &= np.array_equal(np.asarray(seg), outb[lo:lo + plan.ncell])
        same &= bool(valid[lo:lo + plan.ncell].all())
    leg("streaming CTI == m10_cti.deinterleave, cell for cell",
        same and len(got) > 0,
        f"{len(got)} FEC Blocks over {nframes} Frames, C = {C0}")

    # -- the FEC: batched vs the scalar PlpChain, byte for byte
    ch = B.PlpChain(plan.core["ninner"], plan.core["mod"], plan.core["rate"],
                    iters=50)
    pipe = M44.LdmPipeline(plan=plan, iters=50, threads=4)
    pipe.fd = fd
    take = got[:6]
    pk, nconv, nbch = pipe.decode_blocks([s for _b, s in take])
    ref = []
    for _b, s in take:
        r = ch.decode(np.asarray(s, np.complex128))
        if not r["converged"]:
            ref.append(None)
            continue
        bb = ch.baseband_packet(r["bits"])
        ref.append(bb["bytes"] if bb["bch_ok"] else None)
    ident = all((a is None) == (b is None) and (a is None or a == b)
                for a, b in zip(pk, ref))
    leg("batched BICM+LDPC+BCH == m6_bicm.PlpChain scalar, BYTE for BYTE",
        ident and nbch == len(take),
        f"{nbch}/{len(take)} Blocks BCH-zero on both paths, "
        f"{sum(1 for a in pk if a is not None)} packets compared")

    st2 = M44.CtiStream(plan)
    st2.reset(0, sr, fbs)
    got2 = []
    for k in range(nframes):
        st2.push(k, cells_exact[k])
        got2 += st2.take_blocks()
    pk2, _nc2, nb2 = pipe.decode_blocks([s for _b, s in got2[:6]])
    leg("2b: cpe_fast decodes the SAME BYTES as the exact CPE",
        pk2 == pk and nb2 == nbch,
        f"{nb2}/{len(got2[:6])} Blocks BCH-zero via the exact CPE, "
        f"Baseband Packets identical: {pk2 == pk}")

    # -- [5] controls
    #
    # A CONTROL DERIVED THROUGH THE SAME IDENTITY AS THE THING IT PERTURBS IS
    # NOT A CONTROL.  The first version of this leg swept `start_row +- 1` and
    # re-derived C from A/322 9.3.9.1 each time -- and every one of them
    # decoded 6/6 BCH-clean.  That is not a leak in the decoder; it is an
    # algebraic near-identity of the CTI map, proved directly:
    #
    #     start_row -> start_row + 1  makes 9.3.9.1 give C -> C - Nrows,
    #     and q(i) = i + Nrows*((start_row+i) mod Nrows) is then UNCHANGED
    #     except where the commutator wraps: 32368 of a Block's 32400 cells
    #     land on the same received cell.  32 wrong cells of 32400 is deep
    #     inside a 64K LDPC's reach at 13 dB, so of course it decodes.
    #
    # So the controls below perturb the phase in ways the signalled pair
    # cannot absorb: the FEC anchor alone, the commutator alone, and the
    # ACTUAL E76 failure -- an L1 that belongs to a different Frame than the
    # cells it is applied to.
    print("\n  [5] CONTROLS -- each of these must decode NOTHING")
    ctl = {}

    def run_ctl(name, sr_x, C_x, nrows_x=None, cell_src=None):
        s2 = M44.CtiStream(plan)
        if nrows_x:
            s2.nrows = nrows_x
        s2.reset(0, sr, fbs)
        s2.start_row = sr_x % s2.nrows
        s2.C0 = C_x
        s2.out_ptr = C_x
        src = cell_src if cell_src is not None else cells
        blocks = []
        for k in range(len(src)):
            s2.push(k, src[k])
            blocks += s2.take_blocks()
        if not blocks:
            ctl[name] = ("produced no complete Block", 0, 0)
            return
        _p, nc, nb = pipe.decode_blocks([s for _b, s in blocks[:6]])
        ctl[name] = ("decoded", nc, nb)

    run_ctl("FEC anchor off by ONE CELL (C+1)", sr, C0 + 1)
    run_ctl("commutator phase off by one row (C held)", sr + 1, C0)
    run_ctl("commutator phase off by 8 rows (C held)", sr + 8, C0)
    run_ctl("wrong Nrows (887)", sr, C0, nrows_x=887)
    # E76's ACTUAL failure: this Frame's L1 applied to the NEXT Frame's cells
    run_ctl("L1 of Frame 0 applied to Frame 1's cells", sr, C0,
            cell_src=cells[1:])
    # CTI bypassed entirely
    raw_blocks = [recv[C0 + m * plan.ncell:C0 + (m + 1) * plan.ncell]
                  for m in range(6)]
    _p, nc0, nb0 = pipe.decode_blocks(raw_blocks)
    ctl["CTI bypassed"] = ("decoded", nc0, nb0)
    allzero = True
    for name, (how, nc, nb) in ctl.items():
        # THE ORACLE IS THE BCH SYNDROME, NOT LDPC CONVERGENCE.  A one-cell
        # anchor error made min-sum settle on *a* codeword once in six, and
        # the outer BCH rejected it -- which is the layering doing its job.
        # The pipeline itself only ever emits a Baseband Packet on a ZERO BCH
        # syndrome, so that is the quantity a control has to drive to zero.
        allzero &= (nb == 0)
        print(f"        {name:44s} {how:26s} {nc} converged, {nb} BCH0")
    leg("every wrong-phase / wrong-interleaver control emits NO PACKET",
        allzero, "a gate whose controls also pass is not a gate; "
                 "LDPC convergence alone is not the oracle, the BCH "
                 "syndrome is (one control converged 1/6 and was rejected)")

    # and the degeneracy itself, stated as a measurement rather than a claim
    i0 = np.arange(C0, C0 + plan.ncell, dtype=np.int64)
    C1 = CTI.solve_C(fbs, sr + 1, plan.nrows)
    i1 = np.arange(C1, C1 + plan.ncell, dtype=np.int64)
    qa = i0 + plan.nrows * ((sr + i0) % plan.nrows)
    qb = i1 + plan.nrows * ((sr + 1 + i1) % plan.nrows)
    ident = int(np.sum(qa == qb))
    leg("MEASURED: (start_row, C) from 9.3.9.1 is degenerate under +1",
        ident == plan.ncell - plan.ncell // plan.nrows - (
            1 if plan.ncell % plan.nrows else 0) or ident > 0.99 * plan.ncell,
        f"start_row+1 with C re-derived shifts C by {C1 - C0} and leaves "
        f"{ident}/{plan.ncell} cells of a FEC Block IDENTICAL "
        f"({100.0 * ident / plan.ncell:.4f}%) -- which is why it is not a "
        f"usable control")
    return got


# ---------------------------------------------------------------------------
def gate_acquire(offsets=(0.0, 17.0, 41.0, 88.0)):
    print("\n  [7][8] phase acquisition with NO hand-set start Frame")
    rows = []
    for off in offsets:
        y, step, ps = load_frames(FOX, 6, off)
        w0 = y[:1_700_000]
        t0, _c = M44.plateau_t0(w0, M44.BOOTSTRAP, ps)
        L = M44.l1_from_window(w0, t0, ps=ps)
        if not (L and L.get("ok")):
            rows.append((off, None, 0, 0))
            continue
        plan = M44.LdmPlan.from_l1_result(L, "acq")
        fd = M44.FrameDemod(plan, threads=4)
        pipe = M44.LdmPipeline(plan=plan, iters=50, threads=4)
        pipe.adopt(plan)
        pipe.fd = fd
        pipe.log = lambda *_a, **_k: None
        for k in range(5):
            w = y[k * step:k * step + plan.frame_window + 4096]
            tk, _cc = M44.plateau_t0(w, M44.BOOTSTRAP, ps)
            pipe.push_frame(k, w, tk, ps=ps)
        rows.append((off, pipe.cti.start_row, pipe.n_blocks, pipe.n_bch))
    ok = all(r[1] is not None and r[2] > 0 and r[3] == r[2] for r in rows)
    leg("acquires and decodes from every start offset, unaided",
        ok, "; ".join(f"t={o:.0f}s start_row={s} {b}/{n} BCH0"
                      for o, s, n, b in rows))

    # P4: the commutator identity, on independently decoded live L1s
    y, step, ps = load_frames(FOX, 10, 0.0)
    srs = []
    plan = None
    for k in range(10):
        w = y[k * step:k * step + 400_000]
        tk, _c = M44.plateau_t0(w, M44.BOOTSTRAP, ps)
        L = M44.l1_from_window(w, tk, ps=ps)
        if L and L.get("ok"):
            if plan is None:
                plan = M44.LdmPlan.from_l1_result(L, "p4")
            srs.append((k, int(M44.start_row_of(L, plan.core["id"])
                               ["start_row"])))
    bad = [(a, s1, b, s2) for (a, s1), (b, s2) in zip(srs, srs[1:])
           if (s1 + (b - a) * plan.plp_size) % plan.nrows != s2]
    leg("P4 commutator continuity across independently decoded L1s",
        len(srs) >= 8 and not bad,
        f"{len(srs)}/10 Frames verified L1, {len(bad)} disagreements")


# ---------------------------------------------------------------------------
def gate_rf33_identity(baseline_sha, frames=24):
    print("\n  [6] RF33 BYTE IDENTITY -- the user's daily television")
    out = os.path.join(HERE, "e82_out", "rf33_gate.dg")
    if os.path.exists(out):
        os.remove(out)
    cmd = [PY, os.path.join(HERE, "m11_watch.py"), "--capture", RF33,
           "--rate", "6912000", "--accel", "cpu", "--frames", str(frames),
           "--player", "none", "--assets", "all", "--threads", "4",
           "--fe-threads", "4", "--decode-procs", "2", "--dump-dg", out]
    t = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    sha = (hashlib.sha256(open(out, "rb").read()).hexdigest()
           if os.path.exists(out) else None)
    fec = [l for l in p.stdout.splitlines() if "FEC Blocks" in l]
    leg("RF33 datagram stream unchanged by everything in E82",
        sha == baseline_sha,
        f"sha256 {sha}\n             baseline  {baseline_sha}\n"
        f"             {fec[0].strip() if fec else p.stdout[-200:]}  "
        f"({time.time() - t:.0f} s)")
    # and the LDM sniff must decline RF33 on the evidence, not on a name
    plan, info = M44.sniff(
        MP.read_block(RF33, 0, int(0.9 * RATE), MP.guess_format(RF33)), RATE)
    declined = plan is None or not plan.uses_cti
    leg("the L1 sniff DECLINES RF33 (no channel number involved)", declined,
        (f"L1 says TI mode {plan.core.get('ti_mode')}, LDM={plan.is_ldm}"
         if plan is not None else f"L1 not read: {info.get('why')}"))


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--baseline-sha", default=None,
                    help="the pre-E82 RF33 datagram SHA-256")
    ap.add_argument("--skip", default="")
    ap.add_argument("--json", default=os.path.join(HERE, "e82_out",
                                                   "gate_e82.json"))
    a = ap.parse_args()
    skip = set(a.skip.split(",")) if a.skip else set()
    print("gate_e82 -- the LDM/CTI live path")
    print("=" * 74)
    if "0" not in skip:
        gate_l1basic_repetition()
    plan = fd = None
    if "1" not in skip:
        plan, fd, y, step, ps = gate_demod()
        if "3" not in skip:
            gate_cti_and_fec(plan, fd, y, step, ps)
    if "7" not in skip:
        gate_acquire()
    if "6" not in skip and a.baseline_sha:
        gate_rf33_identity(a.baseline_sha)
    npass = sum(1 for _n, ok, _d in RESULTS if ok)
    print("\n" + "=" * 74)
    print(f"  {npass}/{len(RESULTS)} legs pass")
    os.makedirs(os.path.dirname(a.json), exist_ok=True)
    json.dump([dict(leg=n, pass_=ok, detail=d) for n, ok, d in RESULTS],
              open(a.json, "w"), indent=1)
    return 0 if npass == len(RESULTS) else 1


if __name__ == "__main__":
    raise SystemExit(main())
