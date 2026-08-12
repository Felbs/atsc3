#!/usr/bin/env python3
"""E69 -- whole-capture decode for RF30 (LDM, CTI, 64K LDPC).

RF30 is a two-layer LDM multiplex, not a two-subframe one:

    Subframe 0, 174 symbols, 8K FFT, GI 1536
      PLP 0  Core Layer      QPSK      6/15   64K LDPC  CTI  lls_flag=1
      PLP 1  Enhanced Layer  64QAM-NUC 6/15   64K LDPC  CTI  inj 4.0 dB

so both PLPs occupy the SAME 1,044,402 cells, superimposed.  The core layer is
decoded by treating the enhanced layer as noise -- exactly what E49 did on
RF25 -- and it is the layer that carries LLS (the SLT).

This is e49_stream's structure with the hardcoded PLP 16 removed: E49's
instruments are read-only to this campaign, and RF30 has no PLP 16 to be the
continuity referee, so the referee here is the core layer's own BCH syndrome.

  collect   demodulate every Frame once (parallel), store the chosen PLP's
            cell stream as a .npy memmap + per-frame dummy-cell SNR
  decode    CTI-deinterleave and LDPC-decode every FEC Block, writing an
            m7_route .bb Baseband stream of the converged blocks

Usage:
    python e69_stream.py collect ../data/rf30/rf30_158p5_0810.cs16 \
        --rate 6.912e6 --l1 e69_l1_rf30.json --frames 40 --prefix e69_smoke
    python e69_stream.py decode e69_smoke --l1 e69_l1_rf30.json --plp 0 \
        --bb e69_smoke_plp0.bb
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m10_core as M10               # noqa: E402
import m10_cti as CTI                # noqa: E402
import m6_bicm as B                  # noqa: E402
import m6_payload as P6              # noqa: E402
import m7_route as R7                # noqa: E402
from e49_dummy_snr import dummy_snr  # noqa: E402


def be_nice():
    """Stay out of the live soak's way -- it owns the radio and the clock."""
    try:
        import psutil
        psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    except Exception:                                           # noqa: BLE001
        pass


_G = {}


def _collect_init(capture, rate, l1, off, plp_id):
    be_nice()
    js = json.load(open(l1))
    _G["g"] = M10.geometry_from_json(js, os.path.basename(capture))
    _G["capture"], _G["rate"], _G["off"] = capture, rate, off
    _G["plps"] = {p["id"]: p for p in _G["g"].plps}
    _G["core_total"] = sum(p["size"] for p in _G["g"].plps if p["layer"] == 0)
    _G["plp_id"] = plp_id


def _collect_one(n):
    g = _G["g"]
    per = M10.frame_samples(g) / M10.FS_POST
    s = max(0.0, _G["off"] + n * per - 0.01)
    try:
        rep = {}
        anchor = _G["plps"][0]
        pool, info, _ = M10.frame_cells(_G["capture"], _G["rate"], None, g, s,
                                        anchor, report=rep)
        pool = M10.cpe(pool, info["symbol_of"], g, anchor,
                       dummy_start=_G["core_total"])
        snr = dummy_snr(pool, _G["core_total"])
        p = _G["plps"][_G["plp_id"]]
        c = pool[p["start"]:p["start"] + p["size"]].astype(np.complex64)
        diag = dict(frame=n, ok=True, t=s, snr_db=snr["snr_db"],
                    data_coh=rep["data_coherence_mean"],
                    pool=rep["pool"], t0=int(rep["t0"]))
        return n, c.tobytes(), diag
    except Exception as exc:                                    # noqa: BLE001
        return n, None, dict(frame=n, ok=False, t=s,
                             err=f"{type(exc).__name__}: {exc}")


def cmd_collect(a):
    js = json.load(open(a.l1))
    g = M10.geometry_from_json(js, os.path.basename(a.capture))
    per = M10.frame_samples(g) / M10.FS_POST
    plps = {p["id"]: p for p in g.plps}
    n0 = plps[a.plp]["size"]
    off = a.start_frame * per
    nsamp = os.path.getsize(a.capture) // 4
    print(f"  capture {nsamp/a.rate:.2f}s, Frame {per*1000:.2f} ms "
          f"-> {int(nsamp/a.rate/per)} Frames available")
    print(f"  collecting PLP {a.plp}: {n0} cells/frame x {a.frames} frames")

    arr = np.lib.format.open_memmap(f"{a.prefix}_plp{a.plp}.npy", mode="w+",
                                    dtype=np.complex64, shape=(a.frames, n0))
    diags = [None] * a.frames
    t0 = time.time()
    done = 0
    with ProcessPoolExecutor(
            max_workers=a.workers, initializer=_collect_init,
            initargs=(os.path.abspath(a.capture), a.rate,
                      os.path.abspath(a.l1), off, a.plp)) as ex:
        for n, b, diag in ex.map(_collect_one, range(a.frames), chunksize=4):
            if b is not None:
                arr[n] = np.frombuffer(b, np.complex64)
            diags[n] = diag
            done += 1
            if done % 25 == 0:
                el = time.time() - t0
                print(f"  {done}/{a.frames} frames  ({el:.0f}s, "
                      f"{el/done*1000:.0f} ms/frame)", flush=True)
    arr.flush()
    bad = [d for d in diags if not d.get("ok")]
    snrs = [d["snr_db"] for d in diags if d.get("ok")]
    print(f"  collected {a.frames - len(bad)}/{a.frames} frames "
          f"({len(bad)} failed -> zero-filled = LLR erasures)")
    if bad[:3]:
        for d in bad[:3]:
            print(f"    frame {d['frame']}: {d.get('err')}")
    if snrs:
        v = np.array(snrs)
        print(f"  dummy-SNR: median {np.median(v):.2f} dB  "
              f"p10 {np.percentile(v, 10):.2f}  p90 {np.percentile(v, 90):.2f}"
              f"  max {v.max():.2f}  min {v.min():.2f}")
    json.dump(dict(capture=os.path.basename(a.capture),
                   start_frame=a.start_frame, frames=a.frames, plp=a.plp,
                   diags=diags),
              open(f"{a.prefix}_diags.json", "w"), indent=1, default=float)
    print(f"  wrote {a.prefix}_plp{a.plp}.npy / {a.prefix}_diags.json")


def deinterleave_chunked(recv, nrows, start_row, chunk=8_000_000):
    L = len(recv)
    out = np.zeros(L, recv.dtype)
    valid = np.zeros(L, bool)
    for lo in range(0, L, chunk):
        hi = min(lo + chunk, L)
        i = np.arange(lo, hi, dtype=np.int64)
        q = i + nrows * ((start_row + i) % nrows)
        ok = q < L
        out[lo:hi][ok] = recv[q[ok]]
        valid[lo:hi] = ok
    return out, valid


_D = {}


def _decode_init(deint_path, ninner, mod, rate, iters, alpha):
    be_nice()
    _D["cells"] = np.load(deint_path, mmap_mode="r")
    _D["ch"] = B.PlpChain(ninner, mod, rate, iters=iters,
                          alpha=(alpha if alpha > 0 else None))


def _decode_one(job):
    m, lo = job
    ch = _D["ch"]
    seg = np.asarray(_D["cells"][lo:lo + ch.cells_per_fec], np.complex128)
    r = ch.decode(seg)
    if not r["converged"]:
        return m, False, False, int(r["unsatisfied"]), None
    bb = ch.baseband_packet(r["bits"])
    return (m, True, bool(bb["bch_ok"]), 0,
            bb["bytes"].tobytes() if hasattr(bb["bytes"], "tobytes")
            else bytes(bb["bytes"]))


def cmd_decode(a):
    js = json.load(open(a.l1))
    g = M10.geometry_from_json(js, os.path.basename(a.prefix))
    plps = {p["id"]: p for p in g.plps}
    plp = plps[a.plp]
    nrows = CTI.nrows_of(plp["cti_depth"], bool(plp["cti_extended"]))
    ncell = plp["cells_per_fec"]
    C0 = CTI.solve_C(plp["cti_fec_block_start"], plp["cti_start_row"], nrows)
    src = np.load(f"{a.prefix}_plp{a.plp}.npy", mmap_mode="r")
    frames, per_frame = src.shape
    print(f"  PLP {a.plp}: {plp['mod']} {plp['rate']} Ninner {plp['ninner']}"
          f"  CTI Nrows {nrows} start_row {plp['cti_start_row']} -> C {C0}")
    print(f"  {frames} frames x {per_frame} cells, {ncell} cells/FEC Block")

    t0 = time.time()
    recv = np.ascontiguousarray(src.reshape(-1))
    cells, valid = deinterleave_chunked(recv, nrows, plp["cti_start_row"])
    del recv
    first_bad = np.flatnonzero(~valid)
    nvalid = int(first_bad[0]) if len(first_bad) else len(valid)
    nblk = max(0, (nvalid - C0) // ncell)
    print(f"  CTI: {len(cells)} cells, {nvalid} contiguously valid, "
          f"{nblk} FEC Blocks   ({time.time()-t0:.0f}s)")
    deint_path = f"{a.prefix}_plp{a.plp}_deint.npy"
    np.save(deint_path, cells)
    del cells, valid

    jobs = [(m, C0 + m * ncell) for m in range(nblk)]
    if a.max_blocks:
        jobs = jobs[:a.max_blocks]
    res = {}
    t0 = time.time()
    done = conv = 0
    with ProcessPoolExecutor(
            max_workers=a.workers, initializer=_decode_init,
            initargs=(deint_path, plp["ninner"], plp["mod"], plp["rate"],
                      a.iters, a.alpha)) as ex:
        for m, ok, bch_ok, unsat, payload in ex.map(_decode_one, jobs,
                                                    chunksize=4):
            res[m] = (ok, bch_ok, unsat, payload)
            done += 1
            conv += ok
            if done % 100 == 0:
                el = time.time() - t0
                print(f"  {done}/{len(jobs)} blocks, {conv} converged "
                      f"({el:.0f}s)", flush=True)
    el = time.time() - t0
    nconv = sum(1 for r in res.values() if r[0])
    nbch = sum(1 for r in res.values() if r[1])
    print(f"  FEC: {nconv}/{len(jobs)} LDPC converged, {nbch} BCH zero  "
          f"({el:.0f}s, {el/max(len(jobs),1)*1000:.0f} ms/block)")

    stream, bounds, blockmap = bytearray(), [], []
    for m in sorted(res):
        ok, bch_ok, unsat, payload = res[m]
        blockmap.append(dict(block=m, cell=C0 + m * ncell, conv=bool(ok),
                             bch=bool(bch_ok), unsat=int(unsat)))
        if not (ok and bch_ok):
            continue
        ptr, pay = P6.bb_split(payload)
        if ptr is not None:
            bounds.append(len(stream) + ptr)
        stream += pay
    if a.bb:
        diags = [dict(frame=0, conv=nconv, bch=nbch, plp16=False, coh=0.0)]
        with open(a.bb, "wb") as fh:
            R7.bb_write(fh, 0.0, diags, bytes(stream), bounds)
        print(f"  wrote {a.bb}  ({len(stream)/1e6:.2f} MB, "
              f"{len(bounds)} ALP anchors)")
    json.dump(dict(plp=a.plp, nrows=nrows, C=C0, blocks=len(jobs),
                   converged=nconv, bch_ok=nbch, cells_per_fec=ncell,
                   blockmap=blockmap),
              open(f"{a.prefix}_plp{a.plp}_blocks.json", "w"), indent=1)
    print(f"  wrote {a.prefix}_plp{a.plp}_blocks.json")


def main():
    be_nice()
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("collect")
    c.add_argument("capture")
    c.add_argument("--rate", type=float, required=True)
    c.add_argument("--l1", required=True)
    c.add_argument("--start-frame", type=int, default=0)
    c.add_argument("--frames", type=int, required=True)
    c.add_argument("--workers", type=int, default=12)
    c.add_argument("--plp", type=int, default=0)
    c.add_argument("--prefix", required=True)
    d = sub.add_parser("decode")
    d.add_argument("prefix")
    d.add_argument("--l1", required=True)
    d.add_argument("--plp", type=int, default=0)
    d.add_argument("--iters", type=int, default=60)
    d.add_argument("--alpha", type=float, default=0.0)
    d.add_argument("--workers", type=int, default=16)
    d.add_argument("--max-blocks", type=int, default=None)
    d.add_argument("--bb", default=None)
    a = ap.parse_args()
    cmd_collect(a) if a.cmd == "collect" else cmd_decode(a)


if __name__ == "__main__":
    main()
