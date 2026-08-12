#!/usr/bin/env python3
"""E58 -- margin recovery on the banked Fox capture: per-cell noise-weighted
LLRs (+ optional exact-BP fallback) against the E49 baseline.

THE MECHANISM (why this is not the albacore CSI null re-run)
------------------------------------------------------------
The CTI (d1024) spreads every 32400-cell FEC Block across ~1.05M received
cells ~ 1 s of air, and this channel BREATHES -11..+7.6 dB on exactly that
timescale (E49).  The shipping demapper divides every cell of a block by ONE
scalar sigma^2, so cells from a -8 dB fade enter the LDPC at the same claimed
confidence as cells from a +7 dB crest.  Per-cell noise weighting restores
the LLR model the decoder assumes.  albacore's CSI weighting measured NEUTRAL
(albacore@0017215) -- but there the variation across a codeword was mild
frequency selectivity; here it is 18 dB of time-selective fading inside one
codeword.  Different cause class, hence the re-try (the failure-library rule).

WHAT THE WEIGHT IS
------------------
sigma^2 per (frame, OFDM symbol), measured from the stored cells themselves:
mean squared distance to the nearest QPSK point over the symbol's PLP0 cells.
That is a direct estimate of the demapper's TOTAL effective noise -- thermal
+ enhanced-layer injection + channel-estimate error -- at the granularity the
fading actually has (symbol = 1.26 ms).  It saturates at deep fades (the
residual cannot exceed the constellation spread), which only COMPRESSES the
weights -- conservative, never inventive.  The LDPC syndrome + BCH stay the
only oracles; a weight cannot manufacture a valid codeword.

ENGINE
------
Batched normalized min-sum (m3_ldpc.min_sum_decode_batch, float32 = the E54
compiled kernel), alpha ladder 1.0/0.85/0.75, iters as E49 (100).  The
baseline leg reruns E49's exact job list with this engine so every comparison
is engine-paired; the float64-vs-float32 delta is reported, not assumed away
(E54's identity boundary: converged rows gate output-identical, non-converged
trajectories may differ).

Optional --sp on failures: exact sum-product (float64, log-tanh) on blocks
whose min-sum best attempt came close (unsatisfied < --sp-cut), the two-stage
rescue.  Cost is offline-fine; the rescue rate is the measurement.

Subcommands:
  sigma   build e58_sigma.npy from e49_full_plp0.npy   (per frame x symbol)
  decode  rerun E49's job list: --mode plain|weighted [--sp] [--permute]
  curve   convergence vs baseline frame-SNR buckets, two runs compared
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m3_ldpc as LD                 # noqa: E402
import m6_bicm as B                  # noqa: E402
import m10_core as M10               # noqa: E402
import m10_cti as CTI                # noqa: E402

SQ2 = np.sqrt(2.0)


def below_normal():
    try:
        import psutil
        psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    except Exception:                                           # noqa: BLE001
        pass


# ---------------------------------------------------------------------------
# geometry helpers
# ---------------------------------------------------------------------------

def pool_symbol_bounds(g):
    """Cumulative pool offsets of [spare, sym0, sym1, ...] -- cell order."""
    lens = [g.preamble_spare()] + [
        (g.n_sbs_active if g.is_sbs(l) else g.n_normal)
        for l in range(g.nsym)]
    return np.cumsum([0] + lens)


def qpsk_d2min(z):
    """min over the 4 QPSK points of |z - p|^2, vectorized (exact)."""
    return (z.real * z.real + z.imag * z.imag + 1.0
            - SQ2 * (np.abs(z.real) + np.abs(z.imag)))


# ---------------------------------------------------------------------------
# sigma: per (frame, symbol-part) effective demap noise
# ---------------------------------------------------------------------------

def cmd_sigma(a):
    js = json.load(open(a.l1))
    g = M10.geometry_from_json(js, "rf25")
    cum = pool_symbol_bounds(g)
    p0 = {p["id"]: p for p in g.plps}[0]
    size = p0["size"]
    # parts fully or partially inside the PLP0 slice [start, start+size)
    src = np.load(a.plp0, mmap_mode="r")
    frames = src.shape[0]
    nparts = len(cum) - 1
    table = np.full((frames, nparts), np.nan, np.float32)
    t0 = time.time()
    for f in range(frames):
        z = np.asarray(src[f])
        d2 = qpsk_d2min(z)
        for k in range(nparts):
            lo = max(int(cum[k]) - p0["start"], 0)
            hi = min(int(cum[k + 1]) - p0["start"], size)
            if hi - lo < 32:
                continue
            table[f, k] = max(float(d2[lo:hi].mean()), 1e-6)
        if (f + 1) % 50 == 0:
            print(f"  {f+1}/{frames} frames ({time.time()-t0:.0f}s)",
                  flush=True)
    np.save(a.out, table)
    good = table[np.isfinite(table)]
    print(f"  wrote {a.out}  shape {table.shape}  "
          f"sigma2 med {np.median(good):.3f}  p10 {np.percentile(good,10):.3f}"
          f"  p90 {np.percentile(good,90):.3f}  max {good.max():.2f}")


# ---------------------------------------------------------------------------
# exact sum-product (float64, log-tanh), batched -- the two-stage rescue
# ---------------------------------------------------------------------------

def sum_product_decode_batch(lam2d, checks, n, iters=100):
    """Exact BP.  Returns (bits2d, conv, its, bad) like the min-sum batch."""
    idx, mask = LD.pack(checks, n)
    nblk = lam2d.shape[0]
    tot = np.empty((nblk, n + 1))
    tot[:, :n] = lam2d
    tot[:, n] = 0.0
    m = np.zeros((nblk,) + idx.shape)
    bits_out = (np.asarray(lam2d) < 0).astype(np.uint8)
    conv_out = np.zeros(nblk, bool)
    its_out = np.full(nblk, iters, np.int32)
    bad_out = np.full(nblk, -1, np.int32)
    act = np.arange(nblk)
    LOGMIN = -700.0
    for it in range(iters + 1):
        gg = tot[:, idx]
        bb = gg < 0
        syn = np.bitwise_xor.reduce(bb & mask, axis=2)
        bad = syn.sum(axis=1, dtype=np.int32)
        done = bad == 0
        if done.any():
            sel = act[done]
            bits_out[sel] = tot[done, :n] < 0
            conv_out[sel] = True
            its_out[sel] = it
            bad_out[sel] = 0
            keep = ~done
            if not keep.any():
                return bits_out, conv_out, its_out, bad_out
            act, tot, m, gg, bad = (act[keep], tot[keep], m[keep],
                                    gg[keep], bad[keep])
        if it == iters:
            bits_out[act] = tot[:, :n] < 0
            bad_out[act] = bad
            break
        v = gg - m
        t = np.tanh(np.clip(v, -38.0, 38.0) * 0.5)
        sg = t < 0
        parity = np.bitwise_xor.reduce(sg & mask, axis=2)
        lt = np.log(np.clip(np.abs(t), 1e-300, 1.0))
        lt[:, ~mask] = 0.0
        lsum = lt.sum(axis=2, keepdims=True)
        le = np.clip(lsum - lt, LOGMIN, -1e-16)      # |prod excl self| < 1
        mag = np.exp(le)
        m_new = 2.0 * np.arctanh(np.clip(mag, 0.0, 1.0 - 1e-16))
        sgn = np.where(parity[:, :, None] ^ sg, -1.0, 1.0)
        m = np.where(mask, sgn * m_new, 0.0)
        # rebuild totals: channel LLR + sum of incoming messages
        na = len(act)
        tot[:, :n] = np.asarray(lam2d)[act]
        tot[:, n] = 0.0
        np.add.at(tot, (np.arange(na)[:, None, None], idx[None, :, :]), m)
        tot[:, n] = 0.0
    return bits_out, conv_out, its_out, bad_out


# ---------------------------------------------------------------------------
# decode
# ---------------------------------------------------------------------------

def build_weights(jobs, ncell, nrows, start_row, per_frame, cum, p0_start,
                  table, permute=False):
    """Per-cell sigma^2 for each job's deinterleaved cell span."""
    if permute:
        rng = np.random.default_rng(58)
        table = table.copy()
        for f in range(table.shape[0]):
            fin = np.isfinite(table[f])
            table[f, fin] = rng.permutation(table[f, fin])
    out = {}
    for (mth, lo) in jobs:
        i = np.arange(lo, lo + ncell, dtype=np.int64)
        q = i + nrows * ((start_row + i) % nrows)
        fr = q // per_frame
        off = q % per_frame
        part = np.searchsorted(cum, off + p0_start, side="right") - 1
        s2 = table[fr, part]
        s2 = np.where(np.isfinite(s2), s2, np.nanmedian(table))
        out[mth] = s2.astype(np.float32)
    return out


def ladder_decode(lam2d, checks, n, iters, alphas=(1.0, 0.85, 0.75)):
    """Batch ladder: try each alpha on the still-failing rows."""
    nblk = lam2d.shape[0]
    bits = np.zeros((nblk, n), np.uint8)
    conv = np.zeros(nblk, bool)
    bad = np.full(nblk, 1 << 30, np.int32)
    alpha_used = np.full(nblk, -1.0)
    rows = np.arange(nblk)
    for al in alphas:
        if not len(rows):
            break
        b, c, i, d = LD.min_sum_decode_batch(lam2d[rows], checks, n,
                                             iters=iters, alpha=al,
                                             dtype=np.float32)
        newly = np.flatnonzero(c)
        if len(newly):
            sel = rows[newly]
            bits[sel] = b[newly]
            conv[sel] = True
            bad[sel] = 0
            alpha_used[sel] = al
        better = d < bad[rows]
        bsel = rows[better]
        bits[bsel] = b[better]
        bad[bsel] = d[better]
        rows = rows[~c]
    return bits, conv, bad, alpha_used


def cmd_decode(a):
    below_normal()
    js = json.load(open(a.l1))
    g = M10.geometry_from_json(js, "rf25")
    plps = {p["id"]: p for p in g.plps}
    plp = plps[a.plp]
    nrows = CTI.nrows_of(plp["cti_depth"], bool(plp["cti_extended"]))
    ncell = plp["cells_per_fec"]
    C0 = CTI.solve_C(plp["cti_fec_block_start"], plp["cti_start_row"], nrows)
    cum = pool_symbol_bounds(g)
    per_frame = plp["size"]

    ref = json.load(open(a.jobs_from))
    jobs = [(e["block"], C0 + e["block"] * ncell) for e in ref["blockmap"]]
    if a.max_blocks:
        jobs = jobs[:a.max_blocks]
    print(f"  {len(jobs)} jobs (E49's list), mode={a.mode}"
          f"{' +sp' if a.sp else ''}{' PERMUTED-CONTROL' if a.permute else ''}")

    cells = np.load(a.deint, mmap_mode="r")
    ch = B.PlpChain(plp["ninner"], plp["mod"], plp["rate"])
    table = np.load(a.sigma) if a.mode == "weighted" else None
    weights = (build_weights(jobs, ncell, nrows, plp["cti_start_row"],
                             per_frame, cum, plp["start"], table,
                             permute=a.permute)
               if a.mode == "weighted" else None)

    def llr_block(mth, lo):
        seg = np.asarray(cells[lo:lo + ncell], np.complex128)
        if weights is None:
            q = B.demap_llr(seg, ch.mod, ch.rate)          # E49's demap
        else:
            q = B.demap_llr(seg, ch.mod, ch.rate, sigma2=1.0)
            q = (q.reshape(-1, ch.mod_bits)
                 / weights[mth][:, None]).ravel()
        lam = np.empty(ch.ninner, np.float32)
        lam[ch.lam_of_q] = q
        return lam

    t0 = time.time()
    res = {}
    lock_print = [0]

    def run_chunk(chunk):
        lam2d = np.stack([llr_block(mth, lo) for mth, lo in chunk])
        bits, conv, bad, al = ladder_decode(lam2d, ch.checks, ch.ninner,
                                            a.iters)
        out = []
        for k, (mth, lo) in enumerate(chunk):
            bch = False
            if conv[k]:
                bb = ch.baseband_packet(bits[k])
                bch = bool(bb["bch_ok"])
            out.append((mth, bool(conv[k]), bch, int(bad[k]),
                        float(al[k]), lam2d[k] if (a.sp and not conv[k]
                                                   and bad[k] < a.sp_cut)
                        else None))
        return out

    chunks = [jobs[i:i + a.chunk] for i in range(0, len(jobs), a.chunk)]
    sp_queue = []
    with ThreadPoolExecutor(max_workers=a.threads) as ex:
        for out in ex.map(run_chunk, chunks):
            for mth, conv, bch, bad, al, lam in out:
                res[mth] = dict(conv=conv, bch=bch, unsat=bad, alpha=al,
                                sp=False)
                if lam is not None:
                    sp_queue.append((mth, lam))
            done = len(res)
            if done - lock_print[0] >= 400:
                lock_print[0] = done
                nc = sum(1 for r in res.values() if r["conv"])
                print(f"  {done}/{len(jobs)}  conv {nc}  "
                      f"({time.time()-t0:.0f}s)", flush=True)

    nconv = sum(1 for r in res.values() if r["conv"])
    nbch = sum(1 for r in res.values() if r["bch"])
    print(f"  min-sum ladder: {nconv}/{len(jobs)} converged, {nbch} BCH zero"
          f"  ({time.time()-t0:.0f}s)")

    if a.sp and sp_queue:
        print(f"  exact-BP fallback on {len(sp_queue)} near-miss failures "
              f"(unsat < {a.sp_cut}) ...")
        t1 = time.time()
        rescued = 0
        for i in range(0, len(sp_queue), a.sp_batch):
            grp = sp_queue[i:i + a.sp_batch]
            lam2d = np.stack([l for _, l in grp]).astype(np.float64)
            bits, conv, its, bad = sum_product_decode_batch(
                lam2d, ch.checks, ch.ninner, iters=a.sp_iters)
            for k, (mth, _) in enumerate(grp):
                if conv[k]:
                    bb = ch.baseband_packet(bits[k])
                    if bb["bch_ok"]:
                        res[mth].update(conv=True, bch=True, unsat=0,
                                        sp=True)
                        rescued += 1
            print(f"    {min(i+a.sp_batch, len(sp_queue))}/{len(sp_queue)}"
                  f"  rescued {rescued}  ({time.time()-t1:.0f}s)",
                  flush=True)
        print(f"  exact BP rescued {rescued} blocks the ladder lost")

    nconv = sum(1 for r in res.values() if r["conv"])
    nbch = sum(1 for r in res.values() if r["bch"])
    print(f"  TOTAL: {nconv}/{len(jobs)} converged, {nbch} BCH zero")
    json.dump(dict(mode=a.mode, sp=bool(a.sp), permute=bool(a.permute),
                   iters=a.iters, jobs=len(jobs), converged=nconv,
                   bch_ok=nbch,
                   blocks={str(m): res[m] for m in sorted(res)}),
              open(a.json, "w"))
    print(f"  wrote {a.json}")


# ---------------------------------------------------------------------------
# curve: convergence vs frame SNR, run A vs run B
# ---------------------------------------------------------------------------

def block_snrs(diags_path, l1_path, jobs_json, plp_id=0):
    """E49's own frame-attribution: max dummy-SNR of frames fr-1..fr+1."""
    js = json.load(open(l1_path))
    g = M10.geometry_from_json(js, "rf25")
    plps = {p["id"]: p for p in g.plps}
    plp = plps[plp_id]
    nrows = CTI.nrows_of(plp["cti_depth"], bool(plp["cti_extended"]))
    ncell = plp["cells_per_fec"]
    C0 = CTI.solve_C(plp["cti_fec_block_start"], plp["cti_start_row"], nrows)
    dg = json.load(open(diags_path))
    fsnr = {x["frame"]: x.get("snr_db", -99) for x in dg["diags"]}
    spread = nrows * (nrows - 1) / 2
    ref = json.load(open(jobs_json))
    out = {}
    for e in ref["blockmap"]:
        lo = C0 + e["block"] * ncell
        fr = int((lo + ncell / 2 + spread) // plp["size"])
        out[e["block"]] = max(fsnr.get(fr, -99), fsnr.get(fr + 1, -99),
                              fsnr.get(fr - 1, -99))
    return out


def cmd_curve(a):
    snr = block_snrs(a.diags, a.l1, a.jobs_from)
    runs = []
    for path in a.runs:
        r = json.load(open(path))
        runs.append((os.path.basename(path), r))
    edges = np.arange(3.5, 8.01, 0.5)
    print(f"  {'bucket':>12s}" + "".join(f"  {name[:24]:>26s}"
                                         for name, _ in runs))
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        line = f"  [{lo:4.1f},{hi:4.1f})"
        row = dict(lo=float(lo), hi=float(hi))
        for name, r in runs:
            ids = [m for m in r["blocks"]
                   if lo <= snr.get(int(m), -99) < hi]
            n = len(ids)
            c = sum(1 for m in ids if r["blocks"][m]["conv"])
            b = sum(1 for m in ids if r["blocks"][m]["bch"])
            assert c == b, f"converged-but-BCH-dirty in {name}!"
            line += f"  {c:5d}/{n:5d} ({100*c/max(n,1):5.1f}%)   "
            row[name] = [c, n]
        rows.append(row)
        print(line)
    # sub-3.5 negative control
    line = "  < 3.5 (ctl)"
    for name, r in runs:
        ids = [m for m in r["blocks"] if snr.get(int(m), -99) < 3.5]
        c = sum(1 for m in ids if r["blocks"][m]["conv"])
        line += f"  {c:5d}/{len(ids):5d}            "
    print(line)
    if a.json:
        json.dump(rows, open(a.json, "w"), indent=1)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("sigma")
    s.add_argument("--plp0", default="e49_full_plp0.npy")
    s.add_argument("--l1", default="m8_l1_rf25_amped.json")
    s.add_argument("--out", default="e58_sigma.npy")
    d = sub.add_parser("decode")
    d.add_argument("--deint", default="e49_full_plp0_deint.npy")
    d.add_argument("--l1", default="m8_l1_rf25_amped.json")
    d.add_argument("--jobs-from", default="e49_deep_plp0_blocks.json")
    d.add_argument("--sigma", default="e58_sigma.npy")
    d.add_argument("--plp", type=int, default=0)
    d.add_argument("--mode", choices=("plain", "weighted"), required=True)
    d.add_argument("--sp", action="store_true")
    d.add_argument("--sp-cut", type=int, default=1500)
    d.add_argument("--sp-iters", type=int, default=100)
    d.add_argument("--sp-batch", type=int, default=8)
    d.add_argument("--permute", action="store_true",
                   help="CONTROL: shuffle the sigma table within each frame")
    d.add_argument("--iters", type=int, default=100)
    d.add_argument("--chunk", type=int, default=16)
    d.add_argument("--threads", type=int, default=8)
    d.add_argument("--max-blocks", type=int, default=None)
    d.add_argument("--json", required=True)
    c = sub.add_parser("curve")
    c.add_argument("runs", nargs="+")
    c.add_argument("--diags", default="e49_full_diags.json")
    c.add_argument("--l1", default="m8_l1_rf25_amped.json")
    c.add_argument("--jobs-from", default="e49_deep_plp0_blocks.json")
    c.add_argument("--json")
    a = ap.parse_args()
    {"sigma": cmd_sigma, "decode": cmd_decode, "curve": cmd_curve}[a.cmd](a)


if __name__ == "__main__":
    main()
