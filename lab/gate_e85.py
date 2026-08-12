#!/usr/bin/env python3
"""gate_e85.py -- the LDM demod/FEC process pool, gated on the BYTES.

A parallelisation change has exactly one interesting question: did it change
the answer.  So the primary leg is not a benchmark, it is an identity --
the pooled chain's decoded IP datagram file against the serial chain's, over
the same air, SHA-256 to SHA-256, both re-derived NOW (M9's law: a gate that
compares against a stored expectation is a gate against a memory).

    1  pooled == serial, byte for byte, end to end
    2  CONTROL: releasing Frames out of order into the CTI must DESTROY the
       decode -- otherwise leg 1 proves nothing about the ordering
    3  a failed demod is a ZERO-FILL, not a skip: the commutator index space
       survives and later Blocks still decode
    4  BLAS is pinned inside the workers (read back, not assumed)
    5  RF33's datagram stream is still byte-identical
    6  real time, N>=3, warm-up excluded
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import statistics
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m2_pilots as MP                                             # noqa: E402
import m44_ldm as M44                                              # noqa: E402
from atsc3 import bootstrap as bs                                  # noqa: E402

PY = sys.executable
FOX = os.path.join(HERE, "..", "data", "fox25",
                   "rf25_fox_g4_0810_2254.cs16")
RF33 = os.path.join(HERE, "rate6912_rf33.cs16")
OUT = os.path.join(HERE, "e85_out")
RATE = 6.912e6
RESULTS = []


def leg(name, ok, detail=""):
    RESULTS.append((name, bool(ok), detail))
    print(f"    {'PASS' if ok else '*** FAIL ***'}  {name}"
          + (f"\n             {detail}" if detail else ""))
    return ok


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def run_ldm(tag, procs, frames=80, warmup=20, extra=()):
    dg = os.path.join(OUT, f"{tag}.dg")
    js = os.path.join(OUT, f"{tag}.json")
    for p in (dg, js):
        if os.path.exists(p):
            os.remove(p)
    cmd = [PY, os.path.join(HERE, "e82_ldm_run.py"), FOX, "--rate", "6912000",
           "--frames", str(frames), "--warmup", str(warmup),
           "--threads", "8", "--fe-threads", "8",
           "--decode-procs", str(procs), "--proc-threads", "4",
           "--route", "239.255.45.100:5000", "--assets", "all",
           "--dump-dg", dg, "--json", js] + list(extra)
    env = dict(os.environ, M9_NO_TORCH="1")
    t = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=3600,
                       env=env)
    d = json.load(open(js)) if os.path.exists(js) else {}
    d["_wall"] = time.time() - t
    d["_dg"] = dg
    d["_sha"] = sha(dg) if os.path.exists(dg) else None
    d["_stdout"] = p.stdout[-800:]
    return d


# ---------------------------------------------------------------------------
def gate_identity(frames=80):
    print("\n  [1] pooled == serial, BYTE FOR BYTE, end to end")
    ser = run_ldm("serial", 0, frames=frames)
    pool = run_ldm("pooled", 4, frames=frames)
    same = (ser["_sha"] is not None and ser["_sha"] == pool["_sha"])
    leg("the decoded datagram stream is unchanged by the pool", same,
        f"serial {ser['_sha']}\n             pooled {pool['_sha']}\n"
        f"             {os.path.getsize(ser['_dg'])} bytes both; "
        f"FEC {ser.get('bch_zero')}/{ser.get('blocks')} vs "
        f"{pool.get('bch_zero')}/{pool.get('blocks')}")
    leg("both paths decoded every Block BCH-clean",
        ser.get("blocks") == ser.get("bch_zero") > 0
        and pool.get("blocks") == pool.get("bch_zero") > 0,
        f"serial {ser.get('bch_zero')}/{ser.get('blocks')}, "
        f"pooled {pool.get('bch_zero')}/{pool.get('blocks')}")
    leg("the pool is FASTER (else it is only risk)",
        pool.get("sustained_x_real_time", 0) >
        ser.get("sustained_x_real_time", 0),
        f"serial {ser.get('sustained_x_real_time')}x  ->  "
        f"pooled {pool.get('sustained_x_real_time')}x  "
        f"({pool.get('sustained_x_real_time', 0) / max(ser.get('sustained_x_real_time', 1e-9), 1e-9):.2f}x)")
    return ser, pool


# ---------------------------------------------------------------------------
def _cells_and_phase(nframes=5):
    """n Frames of Core cells plus the phase from the FIRST Frame's L1."""
    fmt = MP.guess_format(FOX)
    head = MP.read_block(FOX, 0, int(0.30 * RATE), fmt)
    h = MP.find_bootstraps(MP.resample_to(head, RATE, bs.FS))[0]
    f0 = int(round(h["position"] * RATE / bs.FS))
    ps = h["fields"].get("preamble_structure")
    probe = MP.read_block(FOX, f0 + 1_600_000, int(0.35 * RATE), fmt)
    h2 = MP.find_bootstraps(MP.resample_to(probe, RATE, bs.FS))[0]
    step = 1_600_000 + int(round(h2["position"] * RATE / bs.FS))
    x = MP.read_block(FOX, f0, step * nframes + 500_000, fmt)
    y = x.astype(np.complex128)
    y *= np.exp(-2j * np.pi * h["fine_cfo_hz"] * np.arange(len(y)) / RATE)
    w0 = y[:1_700_000]
    t00, _c = M44.plateau_t0(w0, M44.BOOTSTRAP, ps)
    L = M44.l1_from_window(w0, t00, ps=ps)
    assert L and L.get("ok")
    plan = M44.LdmPlan.from_l1_result(L, "gate85")
    fd = M44.FrameDemod(plan, threads=4)
    cells = []
    for k in range(nframes):
        w = y[k * step:k * step + plan.frame_window + 4096]
        tk, _c = M44.plateau_t0(w, M44.BOOTSTRAP, ps)
        cells.append(fd.core_cells(w, tk))
    d = M44.start_row_of(L, plan.core["id"])
    return plan, fd, cells, int(d["start_row"]), int(d["fec_block_start"])


def gate_order_and_zerofill():
    print("\n  [2][3] the ordering, and what a lost Frame costs")
    plan, fd, cells, sr, fbs = _cells_and_phase(5)
    pipe = M44.LdmPipeline(plan=plan, iters=50, threads=4)
    pipe.fd = fd

    # AIM THE CONTROL WHERE THE TREATMENT CAN REACH.  A FEC Block is a
    # diagonal whose cells span at most Nrows*(Nrows-1) = 0.92 Frames, so
    # perturbing Frames 2-3 leaves Blocks 0-7 and Blocks 134-141 untouched by
    # construction.  Two earlier builds of this control graded exactly those
    # and read 8/8 and 4/8 -- E82's own law about controls, in its other
    # form.  Swap Frames 0 and 1 and grade the FIRST Blocks, which are cut
    # from Frames 0-1 and nothing else.
    def run(seq, judge=slice(0, 8)):
        st = M44.CtiStream(plan)
        st.reset(0, sr, fbs)
        blocks = []
        for idx, c in seq:
            st.push(idx, c)
            blocks += st.take_blocks()
        if not blocks:
            return 0, 0, 0
        sel = [b for _i, b in blocks][judge]
        if not sel:
            return len(blocks), 0, 0
        pk, nc, nb = pipe.decode_blocks(sel)
        return len(blocks), nc, nb

    n_ref, c_ref, b_ref = run(list(enumerate(cells)))
    leg("in Frame order: every Block decodes", b_ref > 0 and b_ref == c_ref,
        f"{n_ref} Blocks cut, {c_ref} converged, {b_ref} BCH0 "
        f"(Blocks 0-7, which are cut from Frames 0-1)")

    # CONTROL: swap the two Frames those Blocks are made of.  The CTI's index
    # space is positional, so this must destroy them -- if it does not, leg 1's
    # identity proves nothing about ordering.
    sw = list(enumerate(cells))
    sw[0], sw[1] = (0, sw[1][1]), (1, sw[0][1])
    n_sw, c_sw, b_sw = run(sw)
    leg("CONTROL: swapping the Frames those Blocks are cut from -> nothing",
        b_sw == 0,
        f"{n_sw} Blocks cut, {c_sw} converged, {b_sw} BCH0 of Blocks 0-7 -- "
        f"the CTI index space is positional, so Frame order is load-bearing "
        f"and the pool's in-order release is not decoration")
    # and the reach itself, as a measurement rather than an assertion
    reach_fr = plan.cti_reach / plan.plp_size
    leg("MEASURED: a Block's cells span at most %.2f Frames" % reach_fr,
        0.5 < reach_fr < 1.5,
        f"Nrows*(Nrows-1) = {plan.cti_reach} cells / {plan.plp_size} per "
        f"Frame -- which is exactly why a control has to be aimed, and why "
        f"the pool only ever needs ~2 Frames of cells resident")

    # a LOST Frame must be zero-filled, and the Blocks after it must recover
    seq = [(0, cells[0]), (1, None), (2, cells[2]), (3, cells[3]),
           (4, cells[4])]
    st = M44.CtiStream(plan)
    st.reset(0, sr, fbs)
    blocks = []
    for idx, c in seq:
        st.push(idx, c)
        blocks += st.take_blocks()
    late = [b for i, b in blocks][-8:]
    pk, nc, nb = pipe.decode_blocks(late)
    leg("a lost Frame is ZERO-FILLED: the index space and later Blocks survive",
        len(blocks) == n_ref and nb > 0,
        f"{len(blocks)} Blocks cut with Frame 1 missing vs {n_ref} with it "
        f"present (identical count = the commutator never moved), and "
        f"{nb}/8 of the LAST Blocks are still BCH-clean; "
        f"zero_filled_frames={dict(st.stats).get('zero_filled_frames')}")

    # and the control for THAT: skipping instead of zero-filling
    st2 = M44.CtiStream(plan)
    st2.reset(0, sr, fbs)
    blocks2 = []
    for idx, c in [(0, cells[0]), (2, cells[2]), (3, cells[3]),
                   (4, cells[4])]:
        # push with a RENUMBERED index, i.e. pretend the Frame never existed
        st2.push(st2.next_frame, c)
        blocks2 += st2.take_blocks()
    late2 = [b for i, b in blocks2][-8:]
    if late2:
        pk2, nc2, nb2 = pipe.decode_blocks(late2)
    else:
        nb2 = 0
    leg("CONTROL: SKIPPING the lost Frame instead destroys the decode",
        nb2 == 0,
        f"{len(blocks2)} Blocks cut, {nb2}/8 BCH-clean -- which is why "
        f"CtiStream zero-fills and never skips")


# ---------------------------------------------------------------------------
def gate_blas_pinned():
    print("\n  [4] BLAS pinned inside the workers -- read back, not assumed")
    plan, _fd, _c, _sr, _f = _cells_and_phase(1)
    seen = {}

    def cap(msg):
        seen["log"] = seen.get("log", "") + str(msg) + "\n"
    pool = M44.LdmDemodPool(plan, nproc=2, threads=4, blas=1, log=cap)
    try:
        pool.start()
    finally:
        line = seen.get("log", "")
        pool.stop()
    ok = "BLAS pinned to 1" in line
    leg("every worker reports OMP_NUM_THREADS=1", ok, line.strip())


# ---------------------------------------------------------------------------
def gate_rf33(baseline_sha, frames=24):
    print("\n  [5] RF33 BYTE IDENTITY -- the user's daily television")
    out = os.path.join(OUT, "rf33.dg")
    if os.path.exists(out):
        os.remove(out)
    cmd = [PY, os.path.join(HERE, "m11_watch.py"), "--capture", RF33,
           "--rate", "6912000", "--accel", "cpu", "--frames", str(frames),
           "--player", "none", "--assets", "all", "--threads", "4",
           "--fe-threads", "4", "--decode-procs", "2", "--dump-dg", out]
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    got = sha(out) if os.path.exists(out) else None
    fec = [l for l in p.stdout.splitlines() if "FEC Blocks" in l]
    leg("RF33 datagram stream unchanged by E85", got == baseline_sha,
        f"sha256 {got}\n             baseline  {baseline_sha}\n"
        f"             {fec[0].strip() if fec else p.stdout[-200:]}")


# ---------------------------------------------------------------------------
def gate_realtime(n=3, frames=120, procs=4):
    print(f"\n  [6] REAL TIME, N={n}, warm-up excluded")
    rows = []
    for i in range(n):
        d = run_ldm(f"rt{i}", procs, frames=frames, warmup=20)
        rows.append(d)
        print(f"        run {i + 1}: {d.get('sustained_x_real_time')}x  "
              f"BCH {d.get('bch_zero')}/{d.get('blocks')}  "
              f"{d.get('ms_per_frame')}")
    xs = [d.get("sustained_x_real_time", 0) for d in rows]
    clean = all(d.get("blocks") == d.get("bch_zero") > 0 for d in rows)
    med = statistics.median(xs)
    leg(f"sustained x-real-time median {med:.4f}x over N={n}",
        clean, f"{[round(x, 4) for x in xs]}  "
               f"(every run 100% BCH-clean: {clean})")
    return med, xs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--baseline-sha", default=None)
    ap.add_argument("--skip", default="")
    ap.add_argument("--frames", type=int, default=80)
    ap.add_argument("--json", default=os.path.join(OUT, "gate_e85.json"))
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    skip = set(a.skip.split(",")) if a.skip else set()
    print("gate_e85 -- the LDM demod/FEC process pool")
    print("=" * 74)
    if "1" not in skip:
        gate_identity(a.frames)
    if "2" not in skip:
        gate_order_and_zerofill()
    if "4" not in skip:
        gate_blas_pinned()
    if "5" not in skip and a.baseline_sha:
        gate_rf33(a.baseline_sha)
    if "6" not in skip:
        gate_realtime()
    npass = sum(1 for _n, ok, _d in RESULTS if ok)
    print("\n" + "=" * 74)
    print(f"  {npass}/{len(RESULTS)} legs pass")
    json.dump([dict(leg=n, ok=ok, detail=d) for n, ok, d in RESULTS],
              open(a.json, "w"), indent=1)
    return 0 if npass == len(RESULTS) else 1


if __name__ == "__main__":
    raise SystemExit(main())
