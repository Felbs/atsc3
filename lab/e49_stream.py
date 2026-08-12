#!/usr/bin/env python3
"""E49 -- whole-capture Core-layer decode for RF25 (LDM core, CTI, 64K LDPC).

Two subcommands:

  collect   demodulate every Frame once (parallel), store the PLP 0 and
            PLP 16 cell streams as .npy memmaps + per-frame diagnostics
            (dummy-cell SNR -- the calibrated instrument -- and coherence).

  decode    CTI-deinterleave a stored stream and LDPC-decode every FEC Block
            (parallel), writing per-block verdicts and an m7_route .bb
            Baseband stream of the converged blocks.

The frame loop starts at the L1 json's own frame (--start-frame, see
e49_core.py for why); PLP 16 (QPSK 2/15, outside the LDM span) decoded across
the whole capture is the continuity referee: its blocks fail only where the
stitching or the channel collapses, and it has ~9 dB more margin than PLP 0.
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

# ---------------------------------------------------------------------------
# collect
# ---------------------------------------------------------------------------

_G = {}


def _collect_init(capture, rate, l1, off):
    js = json.load(open(l1))
    _G["g"] = M10.geometry_from_json(js, os.path.basename(capture))
    _G["capture"], _G["rate"], _G["off"] = capture, rate, off
    _G["plps"] = {p["id"]: p for p in _G["g"].plps}
    _G["core_total"] = sum(p["size"] for p in _G["g"].plps if p["layer"] == 0)


def _collect_one(n):
    g = _G["g"]
    per = M10.frame_samples(g) / M10.FS_POST
    s = max(0.0, _G["off"] + n * per - 0.01)
    try:
        rep = {}
        pool, info, _ = M10.frame_cells(_G["capture"], _G["rate"], None, g, s,
                                        _G["plps"][0], report=rep)
        pool = M10.cpe(pool, info["symbol_of"], g, _G["plps"][0],
                       dummy_start=_G["core_total"])
        snr = dummy_snr(pool, _G["core_total"])
        p0 = _G["plps"][0]
        p16 = _G["plps"][16]
        c0 = pool[p0["start"]:p0["start"] + p0["size"]].astype(np.complex64)
        c16 = pool[p16["start"]:p16["start"] + p16["size"]].astype(np.complex64)
        diag = dict(frame=n, ok=True, t=s, snr_db=snr["snr_db"],
                    data_coh=rep["data_coherence_mean"],
                    pool=rep["pool"], t0=int(rep["t0"]))
        return n, c0.tobytes(), c16.tobytes(), diag
    except Exception as exc:                                    # noqa: BLE001
        return n, None, None, dict(frame=n, ok=False, t=s,
                                   err=f"{type(exc).__name__}: {exc}")


def cmd_collect(a):
    js = json.load(open(a.l1))
    g = M10.geometry_from_json(js, os.path.basename(a.capture))
    per = M10.frame_samples(g) / M10.FS_POST
    plps = {p["id"]: p for p in g.plps}
    n0, n16 = plps[0]["size"], plps[16]["size"]
    off = a.start_frame * per

    arr0 = np.lib.format.open_memmap(a.prefix + "_plp0.npy", mode="w+",
                                     dtype=np.complex64,
                                     shape=(a.frames, n0))
    arr16 = np.lib.format.open_memmap(a.prefix + "_plp16.npy", mode="w+",
                                      dtype=np.complex64,
                                      shape=(a.frames, n16))
    diags = [None] * a.frames
    t0 = time.time()
    done = 0
    with ProcessPoolExecutor(
            max_workers=a.workers, initializer=_collect_init,
            initargs=(os.path.abspath(a.capture), a.rate,
                      os.path.abspath(a.l1), off)) as ex:
        for n, b0, b16, diag in ex.map(_collect_one, range(a.frames),
                                       chunksize=4):
            if b0 is not None:
                arr0[n] = np.frombuffer(b0, np.complex64)
                arr16[n] = np.frombuffer(b16, np.complex64)
            diags[n] = diag
            done += 1
            if done % 25 == 0:
                el = time.time() - t0
                print(f"  {done}/{a.frames} frames  ({el:.0f}s, "
                      f"{el/done*1000:.0f} ms/frame)", flush=True)
    arr0.flush()
    arr16.flush()
    bad = [d for d in diags if not d.get("ok")]
    snrs = [d["snr_db"] for d in diags if d.get("ok")]
    print(f"  collected {a.frames - len(bad)}/{a.frames} frames "
          f"({len(bad)} failed -> zero-filled = LLR erasures)")
    if snrs:
        v = np.array(snrs)
        print(f"  dummy-SNR: median {np.median(v):.2f} dB  "
              f"p90 {np.percentile(v, 90):.2f}  max {v.max():.2f}  "
              f">=6.8dB frames: {(v >= 6.8).sum()}/{len(v)}")
    json.dump(dict(capture=os.path.basename(a.capture),
                   start_frame=a.start_frame, frames=a.frames,
                   diags=diags),
              open(a.prefix + "_diags.json", "w"), indent=1, default=float)
    print(f"  wrote {a.prefix}_plp0.npy / _plp16.npy / _diags.json")


# ---------------------------------------------------------------------------
# decode
# ---------------------------------------------------------------------------

def deinterleave_chunked(recv, nrows, start_row, chunk=8_000_000):
    """CTI deinterleave without building 5 GB of index arrays at once."""
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
    print(f"  {frames} frames x {per_frame} cells")

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
    if a.snr_min is not None:
        # decode only blocks whose source Frame's dummy-cell SNR clears the
        # bar (nothing converges below ~4.5 dB on this channel, measured);
        # every --snr-sample'th other block is still tried, as the control
        # that keeps the filter honest.
        dg = json.load(open(f"{a.prefix}_diags.json"))
        fsnr = {x["frame"]: x.get("snr_db", -99) for x in dg["diags"]}
        spread = nrows * (nrows - 1) / 2
        keep = []
        for m, lo in jobs:
            fr = int((lo + ncell / 2 + spread) // per_frame)
            s = max(fsnr.get(fr, -99), fsnr.get(fr + 1, -99),
                    fsnr.get(fr - 1, -99))
            if s >= a.snr_min or (a.snr_sample and m % a.snr_sample == 0):
                keep.append((m, lo))
        print(f"  SNR filter >= {a.snr_min} dB (+1-in-{a.snr_sample} "
              f"control): {len(keep)}/{len(jobs)} blocks selected")
        jobs = keep
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
            if done % 200 == 0:
                el = time.time() - t0
                print(f"  {done}/{len(jobs)} blocks, {conv} converged "
                      f"({el:.0f}s)", flush=True)
    el = time.time() - t0
    nconv = sum(1 for r in res.values() if r[0])
    nbch = sum(1 for r in res.values() if r[1])
    print(f"  FEC: {nconv}/{len(jobs)} LDPC converged, {nbch} BCH zero  "
          f"({el:.0f}s, {el/max(len(jobs),1)*1000:.0f} ms/block)")

    # ---- assemble the Baseband stream of converged blocks ----------------
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
                   per_frame=per_frame, blockmap=blockmap),
              open(f"{a.prefix}_plp{a.plp}_blocks.json", "w"), indent=1)
    print(f"  wrote {a.prefix}_plp{a.plp}_blocks.json")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("collect")
    c.add_argument("capture")
    c.add_argument("--rate", type=float, required=True)
    c.add_argument("--l1", required=True)
    c.add_argument("--start-frame", type=int, default=0)
    c.add_argument("--frames", type=int, required=True)
    c.add_argument("--workers", type=int, default=16)
    c.add_argument("--prefix", required=True)
    d = sub.add_parser("decode")
    d.add_argument("prefix")
    d.add_argument("--l1", required=True)
    d.add_argument("--plp", type=int, default=0)
    d.add_argument("--iters", type=int, default=60)
    d.add_argument("--alpha", type=float, default=0.0,
                   help="0 = the ALPHA_LADDER; >0 pins a single alpha")
    d.add_argument("--workers", type=int, default=24)
    d.add_argument("--max-blocks", type=int, default=None)
    d.add_argument("--snr-min", type=float, default=None,
                   help="decode only blocks from frames at/above this "
                        "dummy-cell SNR (plus the 1-in-N control)")
    d.add_argument("--snr-sample", type=int, default=20)
    d.add_argument("--bb", default=None)
    a = ap.parse_args()
    if a.cmd == "collect":
        cmd_collect(a)
    else:
        cmd_decode(a)


if __name__ == "__main__":
    main()
