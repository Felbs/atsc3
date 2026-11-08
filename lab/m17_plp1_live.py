#!/usr/bin/env python3
"""M17 -- Subframe 1 / PLP 1 decoded LIVE: the multiplex's other television service.

WHAT THIS IS
------------
Everything the live chain decodes today rides Subframe 0 (PLP 0 + PLP 16, the
MMTP service).  Subframe 1 -- 75 symbols of 16K FFT, PLP 1 at 256QAM 11/15,
117 FEC Blocks of 64800 bits a Frame, the ROUTE/DASH services (WHUT 32.1 on
RF33) -- has only ever been read offline by m13_sf1.py, with the numpy
reference LDPC at ~24 s a Frame.  This module is the live decoder for it:
`SubframeDecoder.decode(y, t0)` takes the SAME Frame window and t0 the PLP 0
decoder takes and returns Baseband Packets in the SAME shape, so m11 can feed
a second ALP walker and the staged E47 ROUTE lanes light up.

It is assembled from parts that already earned their numbers:
  * geometry and the HTI reading     m13_sf1 (geometry_for_subframe, the
                                     3 x 39 twisted read order)
  * the batched CPE, demap, LDPC     m9_fast.FrameDecoder's own methods,
    ladder and BCH                   called on the PLP 0 decoder instance (its
                                     thread pool, its float32 fast path)
  * the channel estimate             E104 (m16_margin lever 4): the symbol's
                                     own pilot comb read in the delay domain.
                                     +1.5 dB on this very signal (PRIOR_ART
                                     §13.4), and this is the path it was
                                     measured on.  `comb=False` is the
                                     linear-interpolation control.
What is NEW here is only the per-symbol pool for a Subframe with no Preamble
(batched FFT over the 75 symbols, cached pilot geometry and FI permutations)
and the glue.

THE GATE
--------
`python m17_plp1_live.py hit_rf33.cs16 --frames 2`: for every Frame, the
Baseband Packet stream this decoder produces must equal m13's offline stream
BYTE FOR BYTE (same converged set, same payload bytes, same pointer bounds),
with the estimator off and with it on.  m13 is the referee because it has
decoded this Subframe 117/117 on clean air and its bytes carried WHUT's
service list (M13).

WHAT IT DOES NOT DO YET
-----------------------
The live front end hands the decoder FRAME_WINDOW = 354,344 samples -- the
Subframe 0 part of the 1,708,032-sample Frame.  Subframe 1 needs the whole
Frame in the window (m11_stream / m11_watch, next step); `decode` checks the
window and says so plainly if it is short.  Time smoothing of the channel
(E60) is not applied here: the lab smoother does not work on this Subframe
(PRIOR_ART §13.3) and the live one is bound to the Subframe 0 geometry.
"""
from __future__ import annotations

import argparse
import collections
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m3_freqint as FI                                          # noqa: E402
import m3_spec as S                                              # noqa: E402
import m6_bicm as B                                              # noqa: E402
import m6_cells as C                                             # noqa: E402
import m6_payload as P6                                          # noqa: E402
import m6_tbi as T                                               # noqa: E402
import m13_sf1 as M13                                            # noqa: E402
import m16_margin as M16                                         # noqa: E402
import spec_bitint as BI                                         # noqa: E402
import spec_pilots as P                                          # noqa: E402
from m3_preamble import fft_bins                                 # noqa: E402


def comb_response(noc, pk, hp, nfft, min_period):
    """E104 for ONE symbol's own pilots: the comb's response at every carrier,
    or None (then the caller interpolates linearly).  The comb may start off
    carrier 0 (a data symbol's pilots start at dx * (l mod dy)); the edge
    pilots m6 adds at 0 and noc-1 are not on the comb and are not used --
    the delay-domain estimate extrapolates to the band edges itself."""
    pk = np.asarray(pk, int)
    d = np.diff(pk)
    if len(d) < 8:
        return None
    D = int(np.bincount(d).argmax())
    r = int(np.bincount(pk % D).argmax())
    on = (pk % D) == r
    comb_k = np.arange(r, noc, D)
    if on.sum() < 3 * len(comb_k) // 4:
        return None
    hv = (np.interp(comb_k, pk[on], hp[on].real)
          + 1j * np.interp(comb_k, pk[on], hp[on].imag))
    Hf, used = M16.comb_denoise_grid(hv[:, None], D, nfft, min_period=min_period)
    if not used[0]:
        return None
    H = np.empty(noc, np.complex128)
    n = min(len(Hf[0]), noc - r)
    H[r:r + n] = Hf[0][:n]
    H[r + n:] = Hf[0][n - 1]
    H[:r] = Hf[0][0]
    return H


class SubframeDecoder:
    """PLP `plp_id` of Subframe `subframe` (> 0), decoded from a whole-Frame window."""

    def __init__(self, fd, g0, l1b, l1d, subframe=1, plp_id=None, comb=True,
                 min_period=None):
        self.fd = fd
        self.g0 = g0
        self.subframe = subframe
        self.g = M13.geometry_for_subframe(l1b, l1d, subframe, g0)
        self.off = M13.subframe_offset(g0, subframe)
        plps = {p["id"]: p for p in self.g.plps}
        self.plp = plps[plp_id if plp_id is not None else sorted(plps)[0]]
        p = self.plp
        self.ch = B.PlpChain(p["ninner"], p["mod"], p["rate"], iters=fd.iters)
        self.pts = B.points_for(p["mod"], p["rate"])
        self.ones = B.bit_masks(BI.MOD_BITS[p["mod"]])
        self.comb = bool(comb)
        g = self.g
        # The symbol's own comb is spaced dx*dy; its delay window is
        # nfft/(dx*dy) samples, under the guard on every pattern that matters
        # (RF33 Subframe 1: 1024 of 1536).  An echo later than 3/4 of that
        # window aliases.  The E104 gate measured +1.5 dB on exactly this
        # window on real air, so the default accepts it; a caller can demand
        # more.
        self.min_period = (2 * g.gi) // 3 if min_period is None else int(min_period)
        noc = P.NOC[(g.nfft, g.cred)]
        lo, _ = S.carrier_abs_range(g.nfft, g.cred)
        self.noc, self.lo = noc, lo
        ref = np.array(S.pilot_values(noc, 1.0))
        dx, dy = P.dxdy(g.pattern)
        self.dy = dy
        self.track_cfo = True        # residual CFO from the Subframe's own pilots (see pool)
        self.fi0 = 0                 # m13: the FI counter resets at the Subframe
        self.syms = []
        for l in range(g.nsym):
            sbs = g.is_sbs(l)
            pk = (np.arange(0, noc, dx) if sbs
                  else np.arange(dx * (l % dy), noc, dx * dy))
            pk = np.unique(np.concatenate([pk, [0, noc - 1]]))
            used = np.zeros(noc, bool)
            used[np.array(sorted(P.pilot_carriers(g.nfft, g.cred, g.pattern, l, sbs)), int)] = True
            d = np.flatnonzero(~used)
            perm = (FI.interleaving_sequence(g.nfft, l + self.fi0, len(d), "forward", "i")
                    if g.freq_interleaver else None)
            nd = len(d)
            self.syms.append(dict(l=l, sbs=sbs, pk=pk, bins_pk=fft_bins(g.nfft, lo + pk),
                                  refpk=ref[pk], d=d, bins_d=fft_bins(g.nfft, lo + d),
                                  perm=perm, nd=nd,
                                  bins_d32=np.ascontiguousarray(fft_bins(g.nfft, lo + d), np.int32),
                                  d32=np.ascontiguousarray(d, np.int32),
                                  hmat=np.ascontiguousarray(perm if perm is not None else np.arange(nd), np.int32),
                                  trim=((g.null_low, nd - g.null_high) if sbs else (0, nd))))
        # E108: the C equaliser row with an identity grid map (every carrier its
        # own grid point; the last leans on its neighbour with weight 1)
        import m9_fast as _F9
        self.kern = _F9._cellpool_lib()
        self.jg_id = np.minimum(np.arange(noc), noc - 2).astype(np.int32)
        self.wt_id = np.zeros(noc, np.float64)
        self.wt_id[-1] = 1.0
        # symbols grouped by pilot pattern, for the batched comb estimate
        self.groups = {}
        for sym in self.syms:
            key = ("sbs",) if sym["sbs"] else ("par", sym["l"] % dy)
            grp = self.groups.get(key)
            if grp is None:
                pk = sym["pk"]
                d = np.diff(pk)
                D = int(np.bincount(d).argmax())
                r = int(np.bincount(pk % D).argmax())
                on = (pk % D) == r
                comb_k = np.arange(r, noc, D)
                # a regular comb with no holes: the grid IS the on-comb pilots
                assert len(comb_k) == int(on.sum()) and np.array_equal(pk[on], comb_k), "irregular pilot comb"
                grp = self.groups[key] = dict(ls=[], bins_pk=sym["bins_pk"], refpk=sym["refpk"], on=on, D=D, r=r)
            grp["ls"].append(sym["l"])
        nrows = p["cells_per_fec"]
        ncols = p["n_fec_max"] // p["nti"]
        per_ti = p["size"] // p["nti"]
        self.nrows, self.ncols = nrows, ncols
        self.hidx = np.array([ti * per_ti + T.fec_block(np.arange(per_ti), j, nrows, ncols, 0)
                              for ti in range(p["nti"]) for j in range(ncols)])
        self.pool_len = g.pool_size(0)
        self.reg, self.dv = C.constellation_regions(self.pool_len, p, None)
        self.tm = collections.Counter()

    def prewarm(self):
        """Build the lazily-made tables now, not under the first live Frame:
        the 64800-bit BCH masks (~130 ms), the LDPC kernels' edge tables, the
        FFT plan.  Tables only -- MEASURED: decoding a silent Frame instead
        cost 3-6 s a worker, four workers warmed for 19 s, and on live air
        the radio's buffer overflowed twice while they did."""
        import m3_ldpc as LD
        import m9_accel as A
        g, p = self.g, self.plp
        A._bch_tables(p["ninner"], p["ninner"])
        LD._kernel_tables(self.ch.checks, self.ch.ninner)
        LD._kernel_lib(); LD._layered_lib()
        try:
            from scipy import fft as _sfft
            _sfft.fft(np.zeros((2, g.nfft), np.complex128), axis=1, workers=1)
        except Exception:                                      # noqa: BLE001
            np.fft.fft(np.zeros(g.nfft, np.complex128))
        self._tau0 = {}
        self.df_prev = 0.0

    # -- the window this Subframe needs ------------------------------------
    def window_needed(self, t0):
        g = self.g
        return t0 + self.off + g.gi + (g.nsym - 1) * (g.nfft + g.gi) + g.nfft

    def pool(self, y, t0, rep=None):
        """Equalised cells of the Subframe in cell order, plus who owns each."""
        g = self.g
        need = self.window_needed(t0)
        if len(y) < need:
            raise ValueError(f"Frame window holds {len(y)} samples; Subframe {self.subframe} "
                             f"needs {need} (the front end must hand over the whole Frame)")
        step = g.nfft + g.gi
        span = g.gi + (g.nsym - 1) * step + g.nfft
        seg = y[t0 + self.off:t0 + self.off + span]

        def ffts(sig):
            base = sig[g.gi:]
            W = np.lib.stride_tricks.as_strided(
                base, shape=(g.nsym, g.nfft), strides=(step * base.strides[0], base.strides[0]))
            try:
                from scipy import fft as _sfft
                Y = _sfft.fft(W, axis=1, workers=max(1, min(self.fd.threads, 8)))
            except Exception:                                  # noqa: BLE001
                Y = np.fft.fft(W, axis=1)
            return np.fft.fftshift(Y, axes=1)

        # FEED-FORWARD: last Frame's residual offset is applied before the
        # first transform (the drift is ~1 Hz/s, a Frame is a quarter second),
        # so the second transform below is needed only when the residual
        # measured afterwards is still over a hertz -- a jump, not steady
        # drift.  Each decode process tracks its own; on the gate's isolated
        # Frames the first one pays the second pass and the next inherit it.
        df_prev = float(getattr(self, "df_prev", 0.0))
        if abs(df_prev) > 0.5:
            n_ = np.arange(span, dtype=np.float64)
            seg = seg * np.exp(-2j * np.pi * df_prev * n_ / M13.FS_POST)
        _pc = time.perf_counter
        _t = _pc()
        Yb = ffts(seg)
        self.tm["sf1_pool_fft"] += _pc() - _t
        _t = _pc()
        # RESIDUAL CARRIER OFFSET, measured on this Subframe's own pilots.
        # The front end fixes the CFO once per lock; the tuner then drifts
        # (measured 2026-11-07 on the race capture: +471 -> +452 Hz over
        # 20 s).  A residual of 19 Hz is 4.5% of a 16K carrier spacing and
        # the inter-carrier interference it makes (~ -22 dB) sinks 256QAM,
        # while the 8K Subframe 0 at half the fraction never notices -- which
        # is exactly how PLP 0 stayed 100% while PLP 1 fell to 0 late in the
        # capture.  Same-parity symbols (dy apart) carry pilots on the same
        # carriers: the common phase step between them, over all pairs, is
        # the residual; rotate the samples by it and transform again.  Only
        # this Subframe's samples are touched; PLP 0's path is unchanged.
        df = 0.0
        if self.track_cfo:
            acc = 0j
            for l in range(self.dy, g.nsym):
                a, b = self.syms[l], self.syms[l - self.dy]
                if a["sbs"] or b["sbs"]:
                    continue
                ha = Yb[l][a["bins_pk"]] / a["refpk"]
                hb = Yb[l - self.dy][b["bins_pk"]] / b["refpk"]
                acc += np.vdot(hb, ha)
            if abs(acc) > 0:
                tsym = step / M13.FS_POST
                df = float(np.angle(acc)) / (2.0 * np.pi * self.dy * tsym)
                if abs(df) > 1.0:          # a jump: rotate the rest out and transform again
                    n_ = np.arange(span, dtype=np.float64)
                    seg = seg * np.exp(-2j * np.pi * df * n_ / M13.FS_POST)
                    Yb = ffts(seg)
                    self.df_prev = df_prev + df
                else:
                    self.df_prev = df_prev + df   # fold the small residual into next Frame's feed-forward
        self.tm["sf1_pool_cfo"] += _pc() - _t
        _t = _pc()
        if rep is not None:
            rep["cfo_resid_hz"] = df
            rep["cfo_applied_hz"] = df_prev
        # THE CHANNEL ESTIMATE FOR EVERY SYMBOL, in four batched transforms
        # instead of seventy-five: symbols with the same pilot pattern (the
        # SBS ones, and each parity l mod dy) share a comb, so their pilot
        # estimates stack as the columns of one grid and m16_margin lever 4
        # denoises them in one call.  Same arithmetic per symbol as the
        # per-symbol reference (gate: M17 == m13 byte for byte).
        Hcomb = {}
        noc, nfft = self.noc, g.nfft
        if self.comb:
            tau_cache = getattr(self, "_tau0", None)
            if tau_cache is None:
                tau_cache = self._tau0 = {}

            def comb_group(item):
                key, grp = item
                ls = grp["ls"]
                hp_all = np.stack([Yb[l][grp["bins_pk"]] / grp["refpk"] for l in ls], axis=1)   # (npk, ng)
                hv = hp_all[grp["on"]]                                                          # the regular comb
                Hf, used_g, tau = M16.comb_denoise_grid(hv, grp["D"], nfft, min_period=self.min_period,
                                                        tau0=tau_cache.get(key), return_tau0=True)
                if tau is not None:
                    tau_cache[key] = int(np.bincount(np.asarray(tau, int) - int(tau.min())).argmax() + tau.min())
                out = []
                r, nfull = grp["r"], Hf.shape[1]
                for gi_, l in enumerate(ls):
                    if not used_g[gi_]:
                        continue
                    H = np.empty(noc, np.complex128)
                    nn = min(nfull, noc - r)
                    H[r:r + nn] = Hf[gi_][:nn]
                    H[r + nn:] = Hf[gi_][nn - 1]
                    H[:r] = Hf[gi_][0]
                    out.append((l, H))
                return out
            items = list(self.groups.items())
            # the groups' transforms release the GIL: run them side by side
            outs = (list(self.fd.ex.map(comb_group, items)) if self.fd.ex is not None
                    else [comb_group(it) for it in items])
            for out in outs:
                for l, H in out:
                    Hcomb[l] = H
        self.tm["sf1_pool_comb"] += _pc() - _t
        _t = _pc()
        noc, nfft = self.noc, g.nfft
        used = np.zeros(g.nsym, bool)

        kern = self.kern

        def one(sym):
            l = sym["l"]
            row = Yb[l]
            hp = row[sym["bins_pk"]] / sym["refpk"]
            H = Hcomb.get(l)
            coh = float(abs(np.sum(hp[:-1] * np.conj(hp[1:])))
                        / max(np.sum(np.abs(hp) ** 2), 1e-30))
            if H is not None and kern is not None:
                # E108: the C equaliser row (E103/E106's kernel) with the comb
                # estimate read through an identity grid map -- the division
                # mirrors numpy's bit for bit (the kernel's cdiv), the scatter
                # is the FI permutation, and ctypes drops the GIL so the
                # symbols really run in parallel.  Python fallback below.
                used[l] = True
                Hc = np.ascontiguousarray(H, np.complex128)
                xg = np.empty(sym["nd"], np.complex128)
                rc = kern.cellpool_row_sm_f64(
                    np.ascontiguousarray(row, np.complex128).view(np.float64), row.shape[0],
                    Hc.view(np.float64), Hc.shape[0], self.jg_id, self.wt_id, noc, 1.0, 0.0,
                    sym["bins_d32"], sym["d32"], sym["hmat"], sym["nd"], xg.view(np.float64))
                if rc == 0:
                    a, b = sym["trim"]
                    return xg[a:b], coh
            if H is None:
                H = C._interp(np.arange(noc), sym["pk"], hp)
            else:
                used[l] = True
            z = row[sym["bins_d"]] / H[sym["d"]]
            if sym["perm"] is not None:
                x = np.empty_like(z)
                x[sym["perm"]] = z
            else:
                x = z
            if sym["sbs"]:
                x = x[g.null_low:len(x) - g.null_high]
            return x, coh

        res = (list(self.fd.ex.map(one, self.syms)) if self.fd.ex is not None
               else [one(s) for s in self.syms])
        pool = np.concatenate([r[0] for r in res])
        owner = np.concatenate([np.full(len(r[0]), s["l"], int) for r, s in zip(res, self.syms)])
        self.tm["sf1_pool_rows"] += _pc() - _t
        if rep is not None:
            rep.update(pool=int(len(pool)), pool_predicted=int(self.pool_len),
                       data_coherence_mean=float(np.mean([r[1] for r in res])),
                       comb_used=int(used.sum()))
        return pool, owner

    def decode(self, y, t0, rep=None):
        """Baseband Packets of the PLP for this Frame: a list, None where the
        Block did not converge -- the shape m11 already consumes for PLP 0."""
        fd, p, tm, pc = self.fd, self.plp, self.tm, time.perf_counter
        rep = rep if rep is not None else {}
        t = pc()
        pool, symbol_of = self.pool(y, t0, rep)
        tm["sf1_pool"] += pc() - t
        if len(pool) != self.pool_len:
            raise ValueError(f"pool {len(pool)} != predicted {self.pool_len}")
        t = pc()
        pool, _ = fd.cpe_correct(pool, symbol_of, [self.pts], self.reg, self.dv)
        tm["sf1_cpe"] += pc() - t
        t = pc()
        cells = pool[p["start"]:p["start"] + p["size"]][self.hidx]
        if p.get("cell_interleaver"):
            nc = len(cells) // p["nti"]
            for ti in range(p["nti"]):
                cells[ti * nc:(ti + 1) * nc] = T.cell_deinterleave(cells[ti * nc:(ti + 1) * nc])
        tm["sf1_hti"] += pc() - t
        t = pc()
        # The live decoder's demapper: the compiled E102 kernel (per cell,
        # twelve subset minima in one sweep, 256QAM in range), threaded in
        # chunks of Blocks, each Block keeping its OWN sigma^2 exactly as m13
        # and the reference chain do.  MEASURED: m16's numpy GEMM demapper is
        # 4x SLOWER than this kernel on this box (2.3 s vs 0.57 s a Frame).
        ldt = np.float32 if fd.cpu_fast else np.float64
        q = fd.demap_batch(cells, self.pts, self.ones)
        tm["sf1_demap"] += pc() - t
        t = pc()
        lam = np.empty((len(cells), p["ninner"]), ldt)
        lam[:, self.ch.lam_of_q] = q.reshape(len(cells), -1)
        tm["sf1_bit_deint"] += pc() - t
        t = pc()
        bits, conv, its, bad = fd.ldpc_batch(lam, self.ch, None)
        tm["sf1_ldpc"] += pc() - t
        rep["ldpc_iters_mean"] = float(np.mean(its))
        rep["ldpc_iters_max"] = int(np.max(its))
        t = pc()
        try:
            syn, by = fd.baseband_batch(bits, self.ch)
        except Exception:                                      # noqa: BLE001
            # the batched BCH may not serve 64800-bit blocks: the chain's own
            # per-block path does, at a cost
            syn, by = [], []
            for i in range(len(bits)):
                bb = self.ch.baseband_packet(bits[i])
                syn.append(0 if bb["bch_ok"] else 1)
                by.append(bb["bytes"])
        tm["sf1_bch"] += pc() - t
        packets, nconv, nbch = [], 0, 0
        for i in range(len(cells)):
            if conv[i]:
                nconv += 1
                nbch += int(syn[i] == 0)
                packets.append(by[i] if syn[i] == 0 else None)
            else:
                packets.append(None)
        rep.update(converged=nconv, bch_ok=nbch, n_fec=len(cells),
                   unsat_median=(int(np.median(bad[~conv])) if (~conv).any() else 0))
        return packets, rep


def stream_of(packets):
    """m11's reading of a packet list: the payload stream and pointer bounds."""
    stream, bounds = bytearray(), []
    for pk in packets:
        if pk is None:
            continue
        ptr, pay = P6.bb_split(pk)
        if ptr is not None:
            bounds.append(len(stream) + ptr)
        stream += pay
    return bytes(stream), bounds


def main(argv=None):
    ap = argparse.ArgumentParser(description="M17 gate: live Subframe 1 decoder vs m13's offline bytes")
    ap.add_argument("capture", nargs="?", default="hit_rf33.cs16")
    ap.add_argument("--rate", type=float, default=8e6)
    ap.add_argument("--fmt", default=None)
    ap.add_argument("--l1", default="m4_l1detail_hit_rf33.json")
    ap.add_argument("--frames", type=int, default=2)
    ap.add_argument("--start", type=float, default=0.0)
    ap.add_argument("--threads", type=int, default=12)
    ap.add_argument("--iters", type=int, default=50)
    ap.add_argument("--no-ref", action="store_true", help="skip m13 (timing only)")
    ap.add_argument("--json", default=None)
    a = ap.parse_args(argv)

    import m9_fast
    from m2_pilots import load, guess_format
    from m10_core import fine_t0, frame_samples
    from concurrent.futures import ThreadPoolExecutor

    path = a.capture if os.path.isabs(a.capture) else os.path.join(HERE, a.capture)
    js = json.load(open(a.l1 if os.path.isabs(a.l1) else os.path.join(HERE, a.l1)))
    g0 = C.Geometry.rf33_legacy()
    fd = m9_fast.FrameDecoder(threads=a.threads, backend="cpu", iters=a.iters, cpu_fast=True)
    fd.prewarm()
    decs = {name: SubframeDecoder(fd, g0, js["l1_basic"], js["l1_detail"], 1, comb=(name == "comb"))
            for name in ("linear", "comb")}
    g = decs["comb"].g
    fmt = a.fmt or guess_format(path)
    total = frame_samples(g0) + (g.nfft + g.gi) * g.nsym
    per = total / M13.FS_POST
    ex = ThreadPoolExecutor(max_workers=a.threads)
    print("M17 -- Subframe 1 live decoder, gated against m13")
    print("=" * 72)
    print(f"  {os.path.basename(path)}: FFT {g.nfft//1024}K GI {g.gi} {g.pattern} {g.nsym} symbols, "
          f"PLP {decs['comb'].plp['id']} {decs['comb'].plp['mod']} {decs['comb'].plp['rate']}, "
          f"{len(decs['comb'].hidx)} blocks/frame; comb window {g.nfft // (P.dxdy(g.pattern)[0]*P.dxdy(g.pattern)[1])} "
          f"samples, accepted from {decs['comb'].min_period}")
    out = []
    allok = True
    for k in range(a.frames):
        y, _fs, _cfo, _h = load(path, a.rate, fmt=fmt, span_sec=per + 0.06, start_sec=a.start + k * per)
        t0, coh = fine_t0(y, g0)
        row = dict(frame=k, t0=int(t0))
        ref = None
        if not a.no_ref:
            t = time.perf_counter()
            r = M13.decode_subframe(y, t0 + decs["comb"].off, g, decs["comb"].plp, 0,
                                    iters=a.iters, verbose=False, ex=ex)
            row["m13"] = dict(converged=r["converged"], bch_ok=r["bch_ok"], secs=round(time.perf_counter() - t, 1))
            ref = (r["stream"], list(r["bounds"]))
        for name, dec in decs.items():
            dec.tm.clear()
            rep = {}
            t = time.perf_counter()
            pk, rep = dec.decode(y, t0, rep)
            secs = time.perf_counter() - t
            st, bd = stream_of(pk)
            same = None if ref is None else (st == ref[0] and bd == ref[1])
            if same is False:
                allok = False
            row[name] = dict(converged=rep["converged"], bch_ok=rep["bch_ok"], comb_used=rep["comb_used"],
                             iters_mean=round(rep.get("ldpc_iters_mean", 0), 1), iters_max=rep.get("ldpc_iters_max", 0),
                             cfo_resid_hz=round(rep.get("cfo_resid_hz", 0), 2),
                             secs=round(secs, 2), bytes=len(st), same_as_m13=same,
                             ms=dict((kk, round(v * 1000, 1)) for kk, v in dec.tm.items()))
        out.append(row)
        m = row.get("m13")
        print(f"\n  Frame {k}: t0 {t0}" + (f"   m13 {m['converged']}/{len(decs['comb'].hidx)} converged, "
                                           f"{m['bch_ok']} BCH ok, {m['secs']} s" if m else ""))
        for name in ("linear", "comb"):
            r = row[name]
            print(f"    {name:7s} {r['converged']:3d} converged, {r['bch_ok']:3d} BCH ok, "
                  f"comb on {r['comb_used']:2d}/{g.nsym} symbols, {r['bytes']} bytes, {r['secs']:.2f} s"
                  + ("" if r['same_as_m13'] is None else
                     ("   == m13 BYTE FOR BYTE" if r['same_as_m13'] else "   *** DIFFERS FROM m13 ***")))
            print("            " + "  ".join(f"{kk[4:]} {v:.0f}" for kk, v in r["ms"].items())
                  + f"  (ms)   LDPC iterations mean {r['iters_mean']} max {r['iters_max']}   residual CFO {r['cfo_resid_hz']:+.2f} Hz")
    print("\n  " + ("PASS: every Frame equals m13 byte for byte, estimator off and on"
                    if allok and not a.no_ref else ("no reference run" if a.no_ref else "FAIL")))
    if a.json:
        json.dump(out, open(a.json, "w"), indent=1)
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
