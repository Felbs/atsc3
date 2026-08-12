#!/usr/bin/env python3
"""M10 -- the CORE-layer payload chain for the LDM / CTI / 64K-LDPC multiplexes.

M8 decoded L1 off RF25 and RF30 and then stopped, because those multiplexes are
a different machine from RF33: LDM (two layers on the same cells), CTI instead
of the Hybrid time interleaver, Ninner = 64800 instead of 16200, and a
two-symbol Preamble.  This module is the rest of the chain:

    Frames --> L1 (m8_l1) --> Geometry (m6_cells.Geometry.from_l1)
           --> cell pool, per Frame, driven ENTIRELY by L1
           --> the Core PLP's cell stream, concatenated across Frames
           --> CTI de-interleave (m10_cti)
           --> FEC Blocks anchored on L1D_plp_CTI_fec_block_start
           --> BICM + 64800 LDPC + BCH-192 + descramble (m6_bicm)
           --> Baseband Packets, written in m7_route's `.bb` format so the
               existing ALP / IP / LLS / ROUTE / MMTP chain runs unchanged.

WHY NO LDM CANCELLATION
-----------------------
Both channels signal `L1D_plp_lls_flag = 1` on the **Core** layer, so the
service list lives in the layer that is decodable with the Enhanced layer
treated as interference -- which is exactly what LDM is designed for.  And for
a QPSK Core layer the max-log LLR depends only on the SIGN of I and Q, so the
Core demapper does not even need the injection level: it is scale-invariant.
Cancellation would only be needed to reach the Enhanced layer, and nothing
asked for is there.

GATES (each with something that must fail)
------------------------------------------
P1  **pool == prediction.**  The cell pool built from the pilot tables must
    equal `sum(plp_size) + dummy`.  On RF30 there is exactly ONE PLP and no
    dummy tail, so the identity is `pool == L1D_plp_size` with no slack at
    all: 1044402 from Annex F / D.1.4 pilot arithmetic against 1044402 from
    the LDPC-decoded L1-Detail.
P2  **L1D_sbs_null_cells** computed from the pilot tables == signalled.
P3  **the CTI identity** (A/322 9.3.9.1), gated in m10_cti: only the signalled
    Nrows puts C inside [0, cells_per_FEC_Block).
P4  **the commutator advances by plp_size per Subframe.**  The signalled
    `L1D_plp_CTI_start_row` of Frame n+1 must equal
    `(start_row(n) + plp_size) mod Nrows`.  Two independently decoded L1s.
P5  **LDPC 0 unsatisfied + BCH syndrome zero**, per FEC Block -- the oracle.
P6  **the dummy-cell referee** (A/322 7.2.6.5) where a dummy tail exists.
    RF25 has 51257 dummy cells; RF30 has none, which is itself P1.

CONTROLS
--------
`--control` runs the same decode with (a) the CTI bypassed, (b) the CTI
de-interleaved at the WRONG start_row, and (c) the FEC anchor offset by one
cell.  All three must fail to converge.

Usage:
    python m10_core.py long_rf30.cs16 --rate 8e6 --l1 m8_l1_rf30.json \
        --frames 8 --bb m10_rf30.bb
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m3_spec as S                                              # noqa: E402
import m6_bicm as B                                              # noqa: E402
import m6_cells as C                                             # noqa: E402
import m6_payload as P6                                          # noqa: E402
import m7_route as R7                                            # noqa: E402
import m8_l1 as L8                                               # noqa: E402
import m10_cti as CTI                                            # noqa: E402
import m4_scrambler as SC                                        # noqa: E402
from m2_pilots import load                                       # noqa: E402
from m3_preamble import pilot_coherence                          # noqa: E402

BOOTSTRAP = 13824
FS_POST = 6.912e6


# ---------------------------------------------------------------------------
# geometry + timing
# ---------------------------------------------------------------------------

def geometry_from_json(js, label=""):
    pre = js["preamble"]
    return C.Geometry.from_l1(js["l1_basic"], js["l1_detail"],
                              pre["nfft"], pre["gi"], pre["dx"],
                              pre["l1b_mode"], label=label)


def frame_samples(g):
    return (BOOTSTRAP + (g.pre_nfft + g.pre_gi) * g.np_sym
            + (g.nfft + g.gi) * g.nsym)


def fine_t0(y, g, span=24):
    """FFT-window start of Preamble symbol 0, by pilot coherence.

    Same statistic m3_preamble.analyse uses, but on a buffer we already hold.
    Returns the m6_cells convention: the FIRST sample of the Preamble,
    i.e. the start of its guard interval.
    """
    base = BOOTSTRAP + g.pre_gi
    best = None
    for t in range(base - span, base + span + 1):
        if t + g.pre_nfft > len(y):
            continue
        Y = np.fft.fftshift(np.fft.fft(y[t:t + g.pre_nfft]))
        c = pilot_coherence(Y, g.pre_nfft, g.pre_gi, g.pre_dx, 4, 0)
        if best is None or c > best[1]:
            best = (t, c)
    return best[0] - g.pre_gi, best[1]


# ---------------------------------------------------------------------------
# one Frame -> the Core PLP's cells
# ---------------------------------------------------------------------------

def frame_cells(path, rate, fmt, g, start_sec, plp, fi_offset=None,
                report=None):
    span = (frame_samples(g) + 40000) / FS_POST * (1.0)
    y, _fs, _cfo, _h = load(path, rate, fmt=fmt,
                            span_sec=span, start_sec=start_sec)
    t0, coh = fine_t0(y, g)
    rep = {}
    pool, info = C.cell_pool_g(y, t0, g, report=rep, fi_offset=fi_offset)
    if report is not None:
        report.update(rep, t0=int(t0), preamble_t0_coherence=float(coh))
    return pool, info, rep


def cpe(pool, symbol_of, g, plp, dummy_start=None):
    """Per-symbol residual gain/phase, decision-directed against the CORE
    alphabet (and the known dummy values where they exist)."""
    reg = np.zeros(len(pool), np.int8)
    dv = np.zeros(len(pool), complex)
    if dummy_start is not None and dummy_start < len(pool):
        reg[dummy_start:] = 2
        dv[dummy_start:] = 1.0 - 2.0 * SC.sequence(len(pool))[dummy_start:]
    alph = [B.points_for(plp["mod"], plp["rate"])]
    out, _gains = C.cpe_correct(pool, symbol_of, alph, reg, dv)
    return out


def dummy_agreement(pool, start):
    d = pool[start:]
    if not len(d):
        return None
    sgn = 1.0 - 2.0 * SC.sequence(len(d) + start)[start:start + len(d)]
    real = np.abs(d.imag) < 0.35 * np.sqrt(np.mean(np.abs(d) ** 2))
    if real.sum() < 50:
        return dict(n=int(len(d)), n_real=int(real.sum()), agreement=None)
    return dict(n=int(len(d)), n_real=int(real.sum()),
                real_frac=float(real.mean()),
                agreement=float(np.mean(np.sign(d.real[real]) == sgn[real])))


# ---------------------------------------------------------------------------
# the whole thing
# ---------------------------------------------------------------------------

def run(path, rate, fmt, js, n_frames, start=0.0, plp_id=0, fi_offset=None,
        control=False, iters=50, verbose=True, max_blocks=None):
    g = geometry_from_json(js, os.path.basename(path))
    plps = {p["id"]: p for p in g.plps}
    plp = plps[plp_id]
    nrows = CTI.nrows_of(plp["cti_depth"], bool(plp["cti_extended"]))
    ncell = plp["cells_per_fec"]
    C0 = CTI.solve_C(plp["cti_fec_block_start"], plp["cti_start_row"], nrows)

    spare = g.preamble_spare()
    pool_pred = g.pool_size(spare)
    core_total = sum(p["size"] for p in g.plps if p["layer"] == 0)
    dummy_n = pool_pred - core_total
    per = frame_samples(g) / FS_POST

    if verbose:
        print("  --- geometry, entirely from L1 ---")
        print(f"    FFT {g.nfft//1024}K  GI {g.gi}  Cred {g.cred}  "
              f"{g.pattern}  {g.nsym} symbols  NP {g.np_sym}  "
              f"freq-interleaver {'ON' if g.freq_interleaver else 'BYPASSED'}")
        print(f"    Preamble symbols {g.preamble_cells()[:2]} cells, "
              f"L1-Basic {g.l1b_cells} + L1-Detail {g.l1d_cells} "
              f"-> spare {spare}")
        print(f"    normal {g.n_normal}  SBS active {g.n_sbs_active} "
              f"(null {g.n_null}, signalled {g.sbs_null_signalled})")
        print(f"    [P2] sbs_null computed {g.n_null} == signalled "
              f"{g.sbs_null_signalled}  "
              f"{'PASS' if g.n_null == g.sbs_null_signalled else '*** FAIL ***'}")
        print(f"    [P1] pool {pool_pred} = core PLP cells {core_total} "
              f"+ dummy {dummy_n}   "
              f"{'PASS (exact, no dummy slack)' if dummy_n == 0 else 'PASS'}"
              if dummy_n >= 0 else "*** FAIL: pool < PLP cells ***")
        print(f"    Frame = {frame_samples(g)} samples = {per*1000:.3f} ms")
        print(f"    PLP {plp_id}: {plp['mod']} {plp['rate']} Ninner "
              f"{plp['ninner']}  {ncell} cells/FEC Block")
        print(f"    CTI Nrows {nrows}  start_row {plp['cti_start_row']}  "
              f"fec_block_start {plp['cti_fec_block_start']}  -> C = {C0}")

    # ---- collect the Core PLP's cells across consecutive Frames ----------
    segs, diags = [], []
    for n in range(n_frames):
        s = max(0.0, n * per - 0.01)
        rep = {}
        try:
            pool, info, _ = frame_cells(path, rate, fmt, g, s, plp,
                                        fi_offset=fi_offset, report=rep)
        except Exception as exc:                                # noqa: BLE001
            print(f"    frame {n}: FAILED ({type(exc).__name__}: {exc})")
            break
        if rep["pool"] != pool_pred:
            print(f"    frame {n}: pool {rep['pool']} != predicted "
                  f"{pool_pred}  *** GATE P1 FAIL ***")
            break
        pool = cpe(pool, info["symbol_of"], g, plp,
                   dummy_start=core_total if dummy_n > 0 else None)
        dmy = dummy_agreement(pool, core_total) if dummy_n > 0 else None
        segs.append(pool[plp["start"]:plp["start"] + plp["size"]])
        diags.append(dict(frame=n, t0=rep["t0"],
                          pre_coh=[round(c, 4) for c in
                                   rep["preamble_coherence"]],
                          data_coh=round(rep["data_coherence_mean"], 4),
                          dummy=dmy))
        if verbose:
            d = ("" if dmy is None else
                 f"  dummy {dmy['n_real']}/{dmy['n']} real, agreement "
                 f"{dmy['agreement']:.4f}" if dmy["agreement"] is not None
                 else "  dummy: too few real cells")
            print(f"    frame {n:2d} @{s:6.3f}s  pool {rep['pool']}  "
                  f"pilot-coh pre {rep['preamble_coherence'][0]:.3f} "
                  f"data {rep['data_coherence_mean']:.3f}{d}")
    if not segs:
        return None
    stream_cells = np.concatenate(segs)

    # ---- CTI de-interleave ----------------------------------------------
    sr = plp["cti_start_row"]
    if control == "startrow":
        sr = (sr + 1) % nrows
    if control == "nocti":
        cells, valid = stream_cells, np.ones(len(stream_cells), bool)
    else:
        cells, valid = CTI.deinterleave(stream_cells, nrows, sr)
    anchor = C0 + (1 if control == "anchor" else 0)

    # a FEC Block is usable only if every one of its cells is valid
    nvalid = int(np.flatnonzero(~valid)[0]) if (~valid).any() else len(valid)
    for i in range(len(valid)):                    # first invalid index
        break
    first_bad = np.flatnonzero(~valid)
    nvalid = int(first_bad[0]) if len(first_bad) else len(valid)
    nblk = max(0, (nvalid - anchor) // ncell)
    if max_blocks:
        nblk = min(nblk, max_blocks)
    if verbose:
        print(f"\n  --- CTI: {len(stream_cells)} cells in, {nvalid} "
              f"contiguously valid, {nblk} complete FEC Blocks from C={anchor}")

    ch = B.PlpChain(plp["ninner"], plp["mod"], plp["rate"], iters=iters)
    stream, bounds = bytearray(), []
    conv = bch = 0
    t0 = time.time()
    for m in range(nblk):
        seg = cells[anchor + m * ncell:anchor + (m + 1) * ncell]
        r = ch.decode(seg)
        if not r["converged"]:
            if verbose and m < 12:
                print(f"      block {m:3d}: NOT converged "
                      f"({r['unsatisfied']} unsatisfied of "
                      f"{plp['ninner'] - ch.kldpc})")
            continue
        conv += 1
        bb = ch.baseband_packet(r["bits"])
        bch += bool(bb["bch_ok"])
        if not bb["bch_ok"]:
            continue
        ptr, pay = P6.bb_split(bb["bytes"])
        if ptr is not None:
            bounds.append(len(stream) + ptr)
        stream += pay
    dt = time.time() - t0
    if verbose:
        print(f"  --- FEC: {conv}/{nblk} LDPC converged, {bch}/{nblk} BCH "
              f"zero   ({dt:.1f}s, {dt/max(nblk,1)*1000:.0f} ms/block)")
    return dict(geom=g, plp=plp, nrows=nrows, C=C0, n_frames=len(segs),
                n_blocks=nblk, converged=conv, bch_ok=bch,
                stream=bytes(stream), bounds=bounds, diags=diags,
                pool_pred=pool_pred, dummy_n=dummy_n)


def gate_startrow(path, rate, fmt, js, g, plp, nrows, n=2, verbose=True):
    """P4: the commutator advances by plp_size per Subframe."""
    per = frame_samples(g) / FS_POST
    got = []
    for k in range(n):
        r = L8.solve(path, rate, fmt, max(0.0, k * per - 0.01))
        if r is None or not r.get("ok"):
            got.append(None)
            continue
        d = {e[1]: e[2] for e in r["l1d"]["fields"]}
        got.append(d.get("L1D_plp_CTI_start_row"))
    ok = True
    for k in range(1, n):
        if got[k] is None or got[0] is None:
            ok = False
            continue
        pred = (got[0] + k * plp["size"]) % nrows
        good = pred == got[k]
        ok &= good
        if verbose:
            print(f"    [P4] Frame {k}: start_row signalled {got[k]}, "
                  f"predicted (start_row0 + {k}*plp_size) mod {nrows} = "
                  f"{pred}   {'PASS' if good else '*** FAIL ***'}")
    return ok, got


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("capture")
    ap.add_argument("--rate", type=float, required=True)
    ap.add_argument("--fmt")
    ap.add_argument("--l1", required=True, help="m8_l1_*.json for this channel")
    ap.add_argument("--frames", type=int, default=6)
    ap.add_argument("--plp", type=int, default=0)
    ap.add_argument("--fi-offset", type=int, default=None)
    ap.add_argument("--iters", type=int, default=50)
    ap.add_argument("--max-blocks", type=int, default=None)
    ap.add_argument("--bb", default=None, help="write m7_route .bb cache here")
    ap.add_argument("--controls", action="store_true")
    ap.add_argument("--startrow-gate", action="store_true")
    ap.add_argument("--json")
    a = ap.parse_args()
    path = a.capture if os.path.isabs(a.capture) else os.path.join(HERE,
                                                                   a.capture)
    l1p = a.l1 if os.path.isabs(a.l1) else os.path.join(HERE, a.l1)
    js = json.load(open(l1p))
    print(f"M10 -- Core-layer decode of {os.path.basename(path)}")
    print("=" * 72)

    r = run(path, a.rate, a.fmt, js, a.frames, plp_id=a.plp,
            fi_offset=a.fi_offset, iters=a.iters, max_blocks=a.max_blocks)
    if r is None:
        print("\n  no Frames decoded.")
        return 1

    if a.startrow_gate:
        print("\n  --- [P4] commutator continuity across Frames ---")
        gate_startrow(path, a.rate, a.fmt, js, r["geom"], r["plp"],
                      r["nrows"])

    if a.controls:
        print("\n  --- CONTROLS (these must FAIL) ---")
        for name in ("nocti", "startrow", "anchor"):
            c = run(path, a.rate, a.fmt, js, min(a.frames, 3), plp_id=a.plp,
                    fi_offset=a.fi_offset, control=name, iters=20,
                    verbose=False, max_blocks=8)
            print(f"    {name:10s}: {c['converged']}/{c['n_blocks']} LDPC "
                  f"converged  "
                  f"{'(control failed as required)' if c['converged'] == 0 else '*** CONTROL PASSED -- gate vacuous ***'}")

    if a.bb and r["stream"]:
        bbp = a.bb if os.path.isabs(a.bb) else os.path.join(HERE, a.bb)
        # m7_route's reporting expects the RF33 per-Frame diagnostic shape; the
        # CTI spreads FEC Blocks ACROSS Frames so the counts are per-run, not
        # per-Frame, and they are recorded on the first entry rather than faked
        # per Frame.
        diags = [dict(d) for d in r["diags"]]
        for k, d in enumerate(diags):
            d["conv"] = r["converged"] if k == 0 else 0
            d["bch"] = r["bch_ok"] if k == 0 else 0
            d["plp16"] = False
            d["coh"] = d.get("data_coh", 0.0)
        with open(bbp, "wb") as fh:
            R7.bb_write(fh, 0.0, diags, r["stream"], r["bounds"])
        print(f"\n  wrote {bbp}  ({len(r['stream'])/1e6:.2f} MB Baseband "
              f"stream, {len(r['bounds'])} ALP anchors)")

    dest = a.json or os.path.join(
        HERE, "m10_core_" + os.path.splitext(os.path.basename(path))[0]
        + ".json")
    with open(dest, "w") as fh:
        json.dump(dict(capture=os.path.basename(path), plp=a.plp,
                       nrows=r["nrows"], C=r["C"], frames=r["n_frames"],
                       blocks=r["n_blocks"], converged=r["converged"],
                       bch_ok=r["bch_ok"], pool=r["pool_pred"],
                       dummy=r["dummy_n"], diags=r["diags"],
                       stream_bytes=len(r["stream"])), fh, indent=2,
                  default=float)
    print(f"  wrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
