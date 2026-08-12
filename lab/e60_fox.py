#!/usr/bin/env python3
"""E60 -- the Fox-capture reproduction referee, THROUGH the ported module.

E58's stacked result (smoothed CE + per-cell weights + exact-BP rescue on
E49's 8697-job list: 1056 blocks, all BCH-zero) was produced by the
lab/e58_*.py instruments.  E60 wires those levers into m16_margin (shared by
the live m9_fast/m11 chain and the generic m10 path); THIS driver re-runs the
whole offline stack calling ONLY the ported functions, so the 1056-block
number referees the port itself:

    collect  -> m16_margin.smoothed_pool_g   (the lever the live CE uses)
    sigma    -> m16_margin.sigma_parts
    weights  -> m16_margin.build_weights
    ladder   -> m16_margin.ladder_decode     (LD.min_sum_decode_batch f32)
    rescue   -> m16_margin.sum_product_decode_batch (deadline=None)

The e58_* instruments are untouched references; `compare` diffs this run's
artifacts against the banked e58_sm_* files.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m6_bicm as B                  # noqa: E402
import m10_core as M10               # noqa: E402
import m10_cti as CTI                # noqa: E402
import m16_margin as M16             # noqa: E402
from m2_pilots import load           # noqa: E402
from e49_dummy_snr import dummy_snr  # noqa: E402


def below_normal():
    try:
        import psutil
        psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    except Exception:                                           # noqa: BLE001
        pass


# ---------------------------------------------------------------------------
# collect (e58_collect's whole-capture loop, calling the PORTED pool builder)
# ---------------------------------------------------------------------------

_G = {}


def _init(capture, rate, l1, off, W, detect):
    below_normal()
    js = json.load(open(l1))
    _G["g"] = M10.geometry_from_json(js, os.path.basename(capture))
    _G.update(capture=capture, rate=rate, off=off, W=W, detect=detect)
    _G["plps"] = {p["id"]: p for p in _G["g"].plps}
    _G["core_total"] = sum(p["size"] for p in _G["g"].plps if p["layer"] == 0)


def _one(n):
    g = _G["g"]
    per = M10.frame_samples(g) / M10.FS_POST
    s = max(0.0, _G["off"] + n * per - 0.01)
    try:
        span = (M10.frame_samples(g) + 40000) / M10.FS_POST
        y, _fs, _cfo, _h = load(_G["capture"], _G["rate"], fmt=None,
                                span_sec=span, start_sec=s)
        t0, coh = M10.fine_t0(y, g)
        rep = {}
        pool, info = M16.smoothed_pool_g(y, t0, g, W=_G["W"],
                                         detect=_G["detect"], report=rep)
        pool = M10.cpe(pool, info["symbol_of"], g, _G["plps"][0],
                       dummy_start=_G["core_total"])
        snr = dummy_snr(pool, _G["core_total"])
        p0 = _G["plps"][0]
        c0 = pool[p0["start"]:p0["start"] + p0["size"]].astype(np.complex64)
        diag = dict(frame=n, ok=True, t=s, snr_db=snr["snr_db"],
                    data_coh=rep["data_coherence_mean"],
                    fell_back=rep["fell_back"], pool=rep["pool"])
        return n, c0.tobytes(), diag
    except Exception as exc:                                    # noqa: BLE001
        return n, None, dict(frame=n, ok=False, t=s,
                             err=f"{type(exc).__name__}: {exc}")


def cmd_collect(a):
    js = json.load(open(a.l1))
    g = M10.geometry_from_json(js, os.path.basename(a.capture))
    per = M10.frame_samples(g) / M10.FS_POST
    plps = {p["id"]: p for p in g.plps}
    off = a.start_frame * per
    arr0 = np.lib.format.open_memmap(a.prefix + "_plp0.npy", mode="w+",
                                     dtype=np.complex64,
                                     shape=(a.frames, plps[0]["size"]))
    diags = [None] * a.frames
    t0 = time.time()
    done = 0
    with ProcessPoolExecutor(
            max_workers=a.workers, initializer=_init,
            initargs=(os.path.abspath(a.capture), a.rate,
                      os.path.abspath(a.l1), off, a.window, a.detect)) as ex:
        for n, b0, diag in ex.map(_one, range(a.frames), chunksize=4):
            if b0 is not None:
                arr0[n] = np.frombuffer(b0, np.complex64)
            diags[n] = diag
            done += 1
            if done % 50 == 0:
                el = time.time() - t0
                print(f"  {done}/{a.frames} frames ({el:.0f}s, "
                      f"{el/done*1000:.0f} ms/frame)", flush=True)
    arr0.flush()
    ok = [d for d in diags if d.get("ok")]
    v = np.array([d["snr_db"] for d in ok])
    fb = sum(d.get("fell_back", 0) for d in ok)
    print(f"  {len(ok)}/{a.frames} frames  dummy-SNR med {np.median(v):.2f} "
          f"p90 {np.percentile(v,90):.2f} max {v.max():.2f}  "
          f"fallback symbols {fb}")
    json.dump(dict(capture=os.path.basename(a.capture),
                   start_frame=a.start_frame, frames=a.frames,
                   window=a.window, detect=a.detect, diags=diags),
              open(a.prefix + "_diags.json", "w"), indent=1, default=float)
    print(f"  wrote {a.prefix}_plp0.npy / _diags.json")


# ---------------------------------------------------------------------------
# sigma through the port
# ---------------------------------------------------------------------------

def cmd_sigma(a):
    below_normal()
    js = json.load(open(a.l1))
    g = M10.geometry_from_json(js, "rf25")
    cum = M16.pool_symbol_bounds(g)
    p0 = {p["id"]: p for p in g.plps}[0]
    src = np.load(a.plp0, mmap_mode="r")
    frames = src.shape[0]
    table = np.full((frames, len(cum) - 1), np.nan, np.float32)
    for f in range(frames):
        table[f] = M16.sigma_parts(np.asarray(src[f]), cum, p0["start"],
                                   p0["size"])
    np.save(a.out, table)
    good = table[np.isfinite(table)]
    print(f"  wrote {a.out}  shape {table.shape}  "
          f"sigma2 med {np.median(good):.3f}")


# ---------------------------------------------------------------------------
# compare against the banked e58 artifacts (the instruments as references)
# ---------------------------------------------------------------------------

def cmd_compare(a):
    ok = True
    pa = np.load(a.prefix + "_plp0.npy", mmap_mode="r")
    pb = np.load("e58_sm_plp0.npy", mmap_mode="r")
    same = pa.shape == pb.shape
    if same:
        for f in range(0, pa.shape[0], 16):
            if not np.array_equal(np.asarray(pa[f]), np.asarray(pb[f])):
                same = False
                break
        # full byte check in chunks
        if same:
            same = all(np.array_equal(np.asarray(pa[f]), np.asarray(pb[f]))
                       for f in range(pa.shape[0]))
    ok &= same
    print(f"  collect vs e58_sm_plp0.npy: "
          f"{'IDENTICAL' if same else '*** DIFFER ***'}")
    sa = np.load(a.prefix + "_sigma.npy")
    sb = np.load("e58_sm_sigma.npy")
    same = (sa.shape == sb.shape
            and np.array_equal(np.isnan(sa), np.isnan(sb))
            and np.array_equal(sa[~np.isnan(sa)], sb[~np.isnan(sb)]))
    ok &= same
    print(f"  sigma vs e58_sm_sigma.npy:  "
          f"{'IDENTICAL' if same else '*** DIFFER ***'}")
    return 0 if ok else 1


# ---------------------------------------------------------------------------
# decode: E49's job list through the ported ladder + weights + exact BP
# ---------------------------------------------------------------------------

def cmd_decode(a):
    below_normal()
    js = json.load(open(a.l1))
    g = M10.geometry_from_json(js, "rf25")
    plp = {p["id"]: p for p in g.plps}[a.plp]
    nrows = CTI.nrows_of(plp["cti_depth"], bool(plp["cti_extended"]))
    ncell = plp["cells_per_fec"]
    C0 = CTI.solve_C(plp["cti_fec_block_start"], plp["cti_start_row"], nrows)
    cum = M16.pool_symbol_bounds(g)
    per_frame = plp["size"]

    ref = json.load(open(a.jobs_from))
    jobs = [(e["block"], C0 + e["block"] * ncell) for e in ref["blockmap"]]
    if a.max_blocks:
        jobs = jobs[:a.max_blocks]
    print(f"  {len(jobs)} jobs (E49's list), weighted+sp through m16_margin",
          flush=True)

    cells = np.load(a.deint, mmap_mode="r")
    ch = B.PlpChain(plp["ninner"], plp["mod"], plp["rate"])
    table = np.load(a.sigma)
    weights = M16.build_weights(jobs, ncell, nrows, plp["cti_start_row"],
                                per_frame, cum, plp["start"], table)

    def llr_block(mth, lo):
        seg = np.asarray(cells[lo:lo + ncell], np.complex128)
        q = B.demap_llr(seg, ch.mod, ch.rate, sigma2=1.0)
        q = (q.reshape(-1, ch.mod_bits) / weights[mth][:, None]).ravel()
        lam = np.empty(ch.ninner, np.float32)
        lam[ch.lam_of_q] = q
        return lam

    t0 = time.time()
    res = {}
    sp_queue = []

    def run_chunk(chunk):
        lam2d = np.stack([llr_block(mth, lo) for mth, lo in chunk])
        bits, conv, bad, al = M16.ladder_decode(lam2d, ch.checks, ch.ninner,
                                                a.iters)
        out = []
        for k, (mth, lo) in enumerate(chunk):
            bch = False
            if conv[k]:
                bb = ch.baseband_packet(bits[k])
                bch = bool(bb["bch_ok"])
            out.append((mth, bool(conv[k]), bch, int(bad[k]), float(al[k]),
                        lam2d[k] if (not conv[k] and bad[k] < a.sp_cut)
                        else None))
        return out

    chunks = [jobs[i:i + a.chunk] for i in range(0, len(jobs), a.chunk)]
    npr = [0]
    with ThreadPoolExecutor(max_workers=a.threads) as ex:
        for out in ex.map(run_chunk, chunks):
            for mth, conv, bch, bad, al, lam in out:
                res[mth] = dict(conv=conv, bch=bch, unsat=bad, alpha=al,
                                sp=False)
                if lam is not None:
                    sp_queue.append((mth, lam))
            if len(res) - npr[0] >= 800:
                npr[0] = len(res)
                nc = sum(1 for r in res.values() if r["conv"])
                print(f"  {len(res)}/{len(jobs)}  conv {nc}  "
                      f"({time.time()-t0:.0f}s)", flush=True)

    nconv = sum(1 for r in res.values() if r["conv"])
    nbch = sum(1 for r in res.values() if r["bch"])
    print(f"  min-sum ladder: {nconv}/{len(jobs)} converged, {nbch} BCH zero"
          f"  ({time.time()-t0:.0f}s)", flush=True)

    if sp_queue:
        print(f"  exact-BP (ported) on {len(sp_queue)} near-misses ...",
              flush=True)
        t1 = time.time()
        rescued = 0
        for i in range(0, len(sp_queue), a.sp_batch):
            grp = sp_queue[i:i + a.sp_batch]
            lam2d = np.stack([l for _, l in grp]).astype(np.float64)
            bits, conv, its, bad = M16.sum_product_decode_batch(
                lam2d, ch.checks, ch.ninner, iters=a.sp_iters)
            for k, (mth, _) in enumerate(grp):
                if conv[k]:
                    bb = ch.baseband_packet(bits[k])
                    if bb["bch_ok"]:
                        res[mth].update(conv=True, bch=True, unsat=0, sp=True)
                        rescued += 1
            print(f"    {min(i+a.sp_batch, len(sp_queue))}/{len(sp_queue)}"
                  f"  rescued {rescued}  ({time.time()-t1:.0f}s)", flush=True)
        print(f"  exact BP rescued {rescued} blocks the ladder lost")

    nconv = sum(1 for r in res.values() if r["conv"])
    nbch = sum(1 for r in res.values() if r["bch"])
    print(f"  TOTAL: {nconv}/{len(jobs)} converged, {nbch} BCH zero")
    json.dump(dict(mode="weighted", sp=True, permute=False, iters=a.iters,
                   jobs=len(jobs), converged=nconv, bch_ok=nbch,
                   blocks={str(m): res[m] for m in sorted(res)}),
              open(a.json, "w"))
    print(f"  wrote {a.json}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("collect")
    c.add_argument("capture")
    c.add_argument("--rate", type=float, required=True)
    c.add_argument("--l1", default="m8_l1_rf25_amped.json")
    c.add_argument("--start-frame", type=int, default=2)
    c.add_argument("--frames", type=int, required=True)
    c.add_argument("--window", type=int, default=12)
    c.add_argument("--detect", type=float, default=6.0)
    c.add_argument("--workers", type=int, default=10)
    c.add_argument("--prefix", default="e60_sm")
    s = sub.add_parser("sigma")
    s.add_argument("--plp0", default="e60_sm_plp0.npy")
    s.add_argument("--l1", default="m8_l1_rf25_amped.json")
    s.add_argument("--out", default="e60_sm_sigma.npy")
    p = sub.add_parser("compare")
    p.add_argument("--prefix", default="e60_sm")
    d = sub.add_parser("decode")
    d.add_argument("--deint", default="e60_sm_plp0_deint.npy")
    d.add_argument("--l1", default="m8_l1_rf25_amped.json")
    d.add_argument("--jobs-from", default="e49_deep_plp0_blocks.json")
    d.add_argument("--sigma", default="e60_sm_sigma.npy")
    d.add_argument("--plp", type=int, default=0)
    d.add_argument("--sp-cut", type=int, default=1500)
    d.add_argument("--sp-iters", type=int, default=100)
    d.add_argument("--sp-batch", type=int, default=8)
    d.add_argument("--iters", type=int, default=100)
    d.add_argument("--chunk", type=int, default=16)
    d.add_argument("--threads", type=int, default=8)
    d.add_argument("--max-blocks", type=int, default=None)
    d.add_argument("--json", default="e60_run_sm_weighted_sp.json")
    a = ap.parse_args()
    r = {"collect": cmd_collect, "sigma": cmd_sigma, "compare": cmd_compare,
         "decode": cmd_decode}[a.cmd](a)
    return r or 0


if __name__ == "__main__":
    raise SystemExit(main())
