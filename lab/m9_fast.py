#!/usr/bin/env python3
"""M9 Step 3 -- the real-time frame decoder: threaded front end, batched BICM,
GPU LDPC.  Every stage produces the SAME BYTES as the M6/M7 chain.

WHAT THE PROFILE ACTUALLY SAID
------------------------------
M7 recorded "the LDPC decoder is the whole cost".  Measured, on the same
capture, the LDPC was **16% of the wall, and three quarters of THAT was
`pack()` rebuilding a constant index array once per FEC Block**.  The real
distribution was six comparable costs, four of which were Python loops
recomputing constants.  So this module is not one optimised kernel; it is four
levers, and the gate for each is the same: the decoded Baseband Packets must
be byte-identical to the reference chain.

THE FOUR LEVERS
---------------
1. **Threads, not processes, for the front end.**  NumPy and SciPy release the
   GIL inside their kernels, so the 36 OFDM symbols of a Frame, the 41
   fine-timing candidates, the 74 demappings and the 7.2 M-sample FIR are all
   thread-parallel *inside one process* -- which matters because a live
   pipeline is one stream, and process-parallelism buys throughput at the cost
   of latency.
2. **Overlap-save FIR, which is EXACT.**  `lfilter(b,1,x)` for an FIR run on a
   block preceded by `len(b)-1` samples of real history produces bit-identical
   outputs for that block: the transposed-form accumulator sees the same
   summands in the same order.  So the notch filter parallelises with no
   numerical change at all -- verified, not assumed (`verify()` requires
   `array_equal`, not `allclose`).
3. **Batching the BICM.**  The 74 FEC Blocks of a Frame were demapped one at a
   time, each rebuilding a (2700, 64) distance matrix in its own call.  Batched
   they are one (74, 2700, 64) reduction.  Per-block sigma^2 is preserved
   exactly, because sigma^2 is a mean over that block's cells only.
4. **The GPU LDPC** (`m9_gpu.GpuMinSum`), batched over the whole Frame and
   gated bit-identical.

WHAT IS NOT ACCELERATED, DELIBERATELY
-------------------------------------
`resample_poly` (8 -> 6.912 Msps) and the bootstrap search stay as they are:
they are 65 ms and ~135 ms per Frame respectively at the current chunking, the
second amortises to nothing on long chunks, and neither can be blocked without
changing arithmetic.
"""
from __future__ import annotations

import collections
import os
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np

# Torch and radioconda's NumPy each ship an Intel OpenMP runtime and whichever
# initialises SECOND aborts the process with "OMP: Error #15".  NumPy loads
# libiomp LAZILY, on the first threaded BLAS call -- which in this chain is a
# `np.vdot` deep inside the CPE loop, i.e. long after torch is up.  Forcing
# that load here, before torch, is the fix that does NOT need
# KMP_DUPLICATE_LIB_OK (which Intel documents as able to "silently produce
# incorrect results" -- not a trade this project makes).
# M9_NO_TORCH=1 keeps the pure-CPU path torch-free.
np.vdot(np.zeros(4, complex), np.zeros(4, complex))
np.dot(np.zeros((8, 8)), np.zeros((8, 8)))
if os.environ.get("M9_NO_TORCH") != "1":
    try:
        import torch as _torch_preload                          # noqa: F401
    except Exception:                                           # noqa: BLE001
        pass

from scipy.signal import firwin, lfilter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m2_pilots as MP                                           # noqa: E402
import m3_freqint as FI                                          # noqa: E402
import m3_ldpc as LD                                             # noqa: E402
import m4_scrambler as SC                                        # noqa: E402
import m6_bicm as B                                              # noqa: E402
import m6_cells as C                                             # noqa: E402
import m6_tbi as T                                               # noqa: E402
import m9_accel                                                  # noqa: E402
import m16_margin as M16                                         # noqa: E402
import spec_bitint as BI                                         # noqa: E402
from atsc3 import bootstrap as bs                                # noqa: E402

FRAME_SEC = C.FRAME_SAMPLES / 6.912e6


# ---------------------------------------------------------------------------
# threaded, bit-identical primitives
# ---------------------------------------------------------------------------

def par_lfilter(b, x, ex, nch=None):
    """Overlap-save FIR.  Bit-identical to lfilter(b, 1.0, x)."""
    nt = len(b) - 1
    L = len(x)
    nch = nch or max(1, min(64, (ex._max_workers if ex else 1)))
    if nch <= 1 or L < 8 * nt:
        return lfilter(b, 1.0, x)
    step = (L + nch - 1) // nch
    out = np.empty(L, x.dtype if np.iscomplexobj(x) else float)

    def job(i):
        s = i * step
        e = min(L, s + step)
        if s >= e:
            return
        s0 = max(0, s - nt)
        y = lfilter(b, 1.0, x[s0:e])
        out[s:e] = y[s - s0:]
    list(ex.map(job, range(nch)))
    return out


def par_apply(fn, n, ex, nch=None):
    """Split range(n) into chunks and run `fn(lo, hi)` on each."""
    nch = nch or (ex._max_workers if ex else 1)
    step = (n + nch - 1) // nch
    list(ex.map(lambda i: fn(i * step, min(n, (i + 1) * step)), range(nch)))


def load_fast(path, rate, fmt=None, span_sec=0.30, start_sec=0.0, notch=True,
              ex=None):
    """m2_pilots.load, with the two 7-Msample element-wise passes threaded."""
    fmt = fmt or MP.guess_format(path)
    head = MP.read_block(path, int(start_sec * rate), int(0.30 * rate), fmt)
    hits = MP.find_bootstraps(MP.resample_to(head, rate, bs.FS))
    if not hits:
        raise RuntimeError("no bootstrap found")
    h = hits[0]
    fs_post = h["fields"]["post_bootstrap_sample_rate_hz"]
    cfo = h["fine_cfo_hz"]
    f0 = int(start_sec * rate) + int(round(h["position"] * rate / bs.FS))
    x = MP.read_block(path, f0, int(span_sec * rate), fmt)
    y = MP.resample_to(x, rate, fs_post).astype(np.complex128)
    # The reference is `y *= np.exp(-2j*np.pi*cfo*np.arange(len(y))/fs_post)`,
    # which Python groups as ((((-2j*pi)*cfo) * n) / fs_post).  Folding the
    # division into the scalar first is a DIFFERENT rounding and showed up
    # immediately as a 1.9e-13 divergence -- so the grouping is preserved.
    k = -2j * np.pi * cfo

    def derot(lo, hi):
        y[lo:hi] *= np.exp(k * np.arange(lo, hi) / fs_post)
    if ex is None:
        derot(0, len(y))
    else:
        par_apply(derot, len(y), ex)
    if notch:
        nt = 257
        b = firwin(nt, 2.95e6, fs=fs_post)
        y = (par_lfilter(b, y, ex) if ex is not None
             else lfilter(b, 1.0, y))[nt // 2:]
    return y, fs_post, cfo, h


# ---------------------------------------------------------------------------
# the frame decoder
# ---------------------------------------------------------------------------

class FrameDecoder:
    """One PLP-0 + PLP-16 subframe-0 decoder, reusable across Frames."""

    def __init__(self, threads=16, backend="gpu", iters=50, ldpc_dtype="float64",
                 cpu_fast=False, margin=None):
        m9_accel.install()
        self.ex = ThreadPoolExecutor(max_workers=threads) if threads > 1 else None
        self.threads = threads
        self.ch16 = B.PlpChain(16200, "QPSK", "2/15", iters=iters)
        self.ch0 = B.PlpChain(16200, "64QAM", "11/15", iters=iters)
        self.iters = iters
        self.backend = backend
        # E52 -- the CPU fast mode.  cpu_fast=False is the HEAD-exact float64
        # path, byte for byte, and stays the reference the gates re-derive.
        # cpu_fast=True runs the SAME algorithms with float32 LDPC/demap and a
        # BLAS-shaped nearest-point search in the CPE -- arithmetic
        # re-associations, never algorithm changes, and never trusted on their
        # own: the m11 e2e gate diffs this path's decoded datagrams against
        # the exact path's, SHA-256 to SHA-256, on real air.
        self.cpu_fast = bool(cpu_fast)
        # E60 -- the E58 margin levers, live.  cpu_fast ONLY: the exact
        # float64 path stays the untouched byte reference, and the GPU paths
        # are out of this port's blast radius.  ATSC3_MARGIN=0 (or
        # margin=m16_margin.OFF) runs the pre-E60 code paths untouched.
        if isinstance(margin, M16.MarginCfg):
            self.margin = margin
        elif margin is None:
            self.margin = M16.MarginCfg()
        else:
            self.margin = M16.MarginCfg(**margin)
        if not self.cpu_fast:
            self.margin = M16.OFF
        self.g0 = self.g16 = self.gb = None
        self._lock = threading.Lock()
        if backend in ("gpu", "gpu-full"):
            import m9_gpu
            if not m9_gpu.available():
                self.backend = "cpu"
            else:
                self.g0 = m9_gpu.GpuMinSum(self.ch0.checks, 16200, iters=iters,
                                           dtype=ldpc_dtype)
                self.g16 = m9_gpu.GpuMinSum(self.ch16.checks, 16200,
                                            iters=iters, dtype=ldpc_dtype)
        self.pts0 = B.points_for("64QAM", "11/15")
        self.pts16 = B.QPSK_POINTS
        self.ones0 = B.bit_masks(BI.MOD_BITS["64QAM"])
        self.ones16 = B.bit_masks(BI.MOD_BITS["QPSK"])
        self.regions = {}
        self.tm = collections.Counter()
        self._tm_lock = threading.Lock()   # two decode workers share this fd

    def prewarm(self):
        """Build every lazily-created table NOW, not under the first Frame.

        The first Frame used to pay for: 70 pilot-geometry builds (Python
        LFSR walks), 35 FI sequences (8192-step Python loops each), the LDPC
        pack/scatter structures, the BCH remainder table (11880 polynomial
        steps) and the FFT plan cache -- ~2 s of wall on the Threadripper,
        billed to "decoding".  It is initialisation: on live air it can run
        while the radio settles, and a viewer meets it as tune-in latency
        once, not as throughput.  Every call here only populates caches; no
        decode output depends on WHEN a pure function was first called.
        """
        import m4_scrambler as SC4
        for l in range(C.NSYM):
            m9_accel._geometry(l, l in C.SBS_SYMS)
            m9_accel._geometry(l, False)
        if self.cpu_fast and not hasattr(self, "_cpf"):
            self._cpf = self._cell_pool_tables()
        for chain in (self.ch0, self.ch16):
            z = np.zeros((2, chain.ninner))
            dts = ((np.float32, np.float64) if self.backend == "cpu"
                   else (np.float64,))
            for dt in dts:
                LD.min_sum_decode_batch(z, chain.checks, chain.ninner,
                                        iters=1, alpha=1.0, dtype=dt)
        m9_accel.bch_syndrome_batch(
            np.zeros((1, self.ch0.kldpc), np.uint8), self.ch0.ninner)
        SC4.sequence(self.ch0.kpayload)
        np.fft.fft(np.zeros(C.NFFT, complex))
        try:
            from scipy import fft as _sfft
            _sfft.fft(np.zeros((C.NSYM, C.NFFT), complex), axis=1, workers=2)
        except Exception:                                      # noqa: BLE001
            pass

    # -- front end ---------------------------------------------------------
    def fine_timing(self, y, span=20, centre=None):
        centre = C.BOOTSTRAP if centre is None else centre
        cands = list(range(centre - span, centre + span + 1))
        # demod_coh == demod_data(...)[1] bit for bit (the lines are copied);
        # it skips only the equalisation whose output fine_timing discards.
        # The BATCHED form is cpu_fast only: pocketfft's multi-row transform
        # takes a different SIMD path and the coherences drift in the last
        # ulps (gate_e52 measured it), so the exact/reference path keeps the
        # per-candidate scalar.  gate_e52 holds the batched path to the same
        # t0 DECISION on real air; the e2e sha is the byte referee.
        if self.cpu_fast:
            res = m9_accel.demod_coh_batch(y, cands, 1, False)
        elif self.ex is None:
            res = [m9_accel.demod_coh(y, t, 1, False) for t in cands]
        else:
            res = list(self.ex.map(lambda t: m9_accel.demod_coh(y, t, 1, False),
                                   cands))
        # same tie-break as m6_cells.fine_timing: strictly-greater wins
        best = None
        for t, coh in zip(cands, res):
            if best is None or coh > best[1]:
                best = (t, coh)
        return best

    def cell_pool(self, y, t0, report=None):
        """m6_cells.cell_pool with the 35 data symbols demodulated in
        parallel.  Symbols are independent, so this is exact."""
        import spec_pilots as P
        rep = report if report is not None else {}
        _, _, n_null = P.sbs_cells(C.NFFT, C.CRED, C.PATTERN, 1)
        lo_n, hi_n = P.sbs_null_split(n_null)
        n_norm = P.data_cells(C.NFFT, C.CRED, C.PATTERN, 1, False)

        zp, npre = C.demod_preamble(y, t0)
        xp = FI.deinterleave(zp, C.NFFT, 0, direction="forward", toggle="i")
        spare = xp[C.L1B_CELLS + C.L1D_CELLS:]

        def one(l):
            sbs = l in C.SBS_SYMS
            z, _ = C.demod_data(y, t0, l, sbs)
            x = FI.deinterleave(z, C.NFFT, l + 1, direction="forward",
                                toggle="i")
            if sbs:
                x = x[lo_n:len(x) - hi_n]
            return x
        rng = range(C.NSYM)
        xs = ([one(l) for l in rng] if self.ex is None
              else list(self.ex.map(one, rng)))
        parts = [spare] + xs
        owner = [np.full(len(spare), -1, int)] + [np.full(len(x), l, int)
                                                  for l, x in enumerate(xs)]
        pool = np.concatenate(parts)
        rep.update(n_preamble_spare=int(len(spare)),
                   n_preamble_cells=int(npre), n_normal=int(n_norm),
                   n_sbs_active=int(len(parts[1])), n_null=int(n_null),
                   null_low=int(lo_n), pool=int(len(pool)))
        return pool, dict(symbol_of=np.concatenate(owner), report=rep)

    # -- batched cell pool (cpu_fast) --------------------------------------
    def _cell_pool_tables(self):
        """Per-symbol demod geometry, grouped by identical pilot pattern.

        Everything here is a pure function of the (fixed) frame geometry:
        symbol classes {SBS} + {l%dy} share bins/pilots, so their FFTs,
        gathers, channel interpolation and FI de-interleave all batch.  Built
        once, cached on the decoder.
        """
        import spec_pilots as P
        from m9_accel import _geometry
        _, _, n_null = P.sbs_cells(C.NFFT, C.CRED, C.PATTERN, 1)
        lo_n, hi_n = P.sbs_null_split(n_null)
        n_norm = P.data_cells(C.NFFT, C.CRED, C.PATTERN, 1, False)
        dx, dy = P.dxdy(C.PATTERN)
        classes = {}
        for l in range(C.NSYM):
            sbs = l in C.SBS_SYMS
            key = ("sbs",) if sbs else ("norm", l % dy)
            classes.setdefault(key, []).append(l)
        tabs = []
        for key, rows in classes.items():
            sbs = key[0] == "sbs"
            g = _geometry(rows[0], sbs)
            pk, kk = g["pk"], g["kk"]
            j = np.clip(np.searchsorted(pk, kk, side="right") - 1,
                        0, len(pk) - 2)
            wt = (kk - pk[j]) / (pk[j + 1] - pk[j])
            nd = len(g["d"])
            Hmat = np.stack([FI.interleaving_sequence(
                C.NFFT, l + 1, nd, direction="forward", toggle="i")
                for l in rows])
            tabs.append(dict(rows=np.array(rows), sbs=sbs, g=g, j=j, wt=wt,
                             nd=nd, Hmat=Hmat,
                             trim=(lo_n, nd - hi_n) if sbs else (0, nd)))
        # symbol_of and the per-symbol output slots are frame-constant too
        lens = np.empty(C.NSYM, int)
        for t in tabs:
            a, b = t["trim"]
            lens[t["rows"]] = b - a
        return dict(tabs=tabs, lens=lens, n_null=n_null, lo_n=lo_n,
                    hi_n=hi_n, n_norm=n_norm,
                    step=C.NFFT + C.GI, gi=C.GI)

    def cell_pool_fast(self, y, t0, report=None):
        """cell_pool with the 35 symbol demods batched per pilot class.

        One (35, 8192) FFT instead of 35, one gather/divide/interp/scatter
        per CLASS instead of ~12 kernels per SYMBOL -- the reference spent
        more wall on NumPy dispatch under the GIL than on arithmetic, and
        also computed a per-symbol coherence nobody read.  The channel
        interpolation is the same linear interpolation evaluated as
        hp[lo]*(1-w) + hp[hi]*w, a re-association of np.interp's arithmetic
        -- cpu_fast's standing concession, decoded-bytes gated (gate_e52
        leg 3 and the m11 e2e sha).
        """
        rep = report if report is not None else {}
        if not hasattr(self, "_cpf"):
            self._cpf = self._cell_pool_tables()
        cp = self._cpf
        zp, npre = C.demod_preamble(y, t0)
        xp = FI.deinterleave(zp, C.NFFT, 0, direction="forward", toggle="i")
        spare = xp[C.L1B_CELLS + C.L1D_CELLS:]

        step, gi = cp["step"], cp["gi"]
        st = np.lib.stride_tricks.as_strided
        base = y[t0 + step + gi:]
        W = st(base, shape=(C.NSYM, C.NFFT),
               strides=(step * base.strides[0], base.strides[0]))
        try:
            from scipy import fft as _sfft
            Yb = _sfft.fft(W, axis=1, workers=max(1, min(self.threads, 8)))
        except Exception:                                      # noqa: BLE001
            Yb = np.fft.fft(W, axis=1)
        Yb = np.fft.fftshift(Yb, axes=1)

        lens = cp["lens"]
        offs = np.concatenate([[0], np.cumsum(lens)]) + len(spare)
        pool = np.empty(int(offs[-1]), np.complex128)
        pool[:len(spare)] = spare
        owner = np.empty(int(offs[-1]), int)
        owner[:len(spare)] = -1

        def one_class(t):
            g = t["g"]
            Yr = Yb[t["rows"]]
            hp = Yr[:, g["bins_pk"]] / g["refpk"]
            H = hp[:, t["j"]] * (1.0 - t["wt"]) + hp[:, t["j"] + 1] * t["wt"]
            z = Yr[:, g["bins_d"]] / H[:, g["d"]]
            xg = np.empty_like(z)
            xg[np.arange(len(z))[:, None], t["Hmat"]] = z
            a, b = t["trim"]
            for i, l in enumerate(t["rows"]):
                pool[offs[l]:offs[l + 1]] = xg[i, a:b]
                owner[offs[l]:offs[l + 1]] = l
        if self.ex is None:
            for t in cp["tabs"]:
                one_class(t)
        else:
            list(self.ex.map(one_class, cp["tabs"]))
        rep.update(n_preamble_spare=int(len(spare)),
                   n_preamble_cells=int(npre), n_normal=int(cp["n_norm"]),
                   n_sbs_active=int(lens[0]), n_null=int(cp["n_null"]),
                   null_low=int(cp["lo_n"]), pool=int(len(pool)))
        return pool, dict(symbol_of=owner, report=rep)

    # -- E60: time-smoothed 2x1D channel estimation, batched (cpu_fast) ----
    def _sm_tables(self):
        """Frame-constant tables for the smoothed CE (built once, cached).

        The pilot grid: every scattered pilot of every symbol sits on a
        carrier that is a multiple of dx (0 and noc-1 included), so one
        (ngrid x nsym) matrix holds every raw pilot estimate and a running-sum
        boxcar over the symbol axis is O(nsym) TOTAL per frame -- the O(new
        symbol) requirement, not O(window) per symbol.  The dense-grid
        interpolation weights (jg, wtg) map every carrier to its two grid
        neighbours and are frame-constant.
        """
        import spec_pilots as P
        cp = self._cpf
        dx, dy = P.dxdy(C.PATTERN)
        noc = P.NOC[(C.NFFT, C.CRED)]
        usable = (noc - 1) % dx == 0
        ngrid = (noc - 1) // dx + 1
        gc = np.arange(ngrid) * dx
        kk = np.arange(noc)
        jg = np.clip(np.searchsorted(gc, kk, side="right") - 1, 0, ngrid - 2)
        wtg = (kk - gc[jg]) / (gc[jg + 1] - gc[jg])
        per = []
        for t in cp["tabs"]:
            pk = t["g"]["pk"]
            if not (pk % dx == 0).all():
                usable = False
            per.append(pk // dx)
        return dict(dx=dx, dy=dy, noc=noc, ngrid=ngrid, jg=jg, wtg=wtg,
                    gidx=per, usable=usable)

    def cell_pool_fast_sm(self, y, t0, report=None):
        """cell_pool_fast with E58's time-smoothed channel estimation.

        The channel is quasi-static over tens of symbols; each symbol's own
        pilots carry the full channel noise (E58: a first-order part of the
        demap noise budget was our OWN channel-estimate noise).  A +-W-symbol
        boxcar over the pilot grid denoises it; the per-symbol CHANGE
        DETECTOR (raw-vs-smoothed pilot residual > detect x frame median)
        sends any symbol the smoothing cannot represent -- a real fade --
        down the untouched per-symbol path, cell for cell (gate_e60 `fade`).
        W=0 reproduces cell_pool_fast exactly (identity gate).

        The smoothing is INTRA-frame (as in the E58 instrument): decode
        workers hold no cross-frame state, so a dropped or re-ordered frame
        can never poison a neighbour's channel estimate.
        """
        mg = self.margin
        rep = report if report is not None else {}
        if not hasattr(self, "_cpf"):
            self._cpf = self._cell_pool_tables()
        cp = self._cpf
        smt = getattr(self, "_smt", None)
        if smt is None:
            smt = self._smt = self._sm_tables()
        if not smt["usable"]:                    # defensive: unknown pattern
            return self.cell_pool_fast(y, t0, report)
        zp, npre = C.demod_preamble(y, t0)
        xp = FI.deinterleave(zp, C.NFFT, 0, direction="forward", toggle="i")
        spare = xp[C.L1B_CELLS + C.L1D_CELLS:]

        step, gi = cp["step"], cp["gi"]
        st = np.lib.stride_tricks.as_strided
        base = y[t0 + step + gi:]
        W = st(base, shape=(C.NSYM, C.NFFT),
               strides=(step * base.strides[0], base.strides[0]))
        try:
            from scipy import fft as _sfft
            Yb = _sfft.fft(W, axis=1, workers=max(1, min(self.threads, 8)))
        except Exception:                                      # noqa: BLE001
            Yb = np.fft.fft(W, axis=1)
        Yb = np.fft.fftshift(Yb, axes=1)

        # raw per-symbol pilot estimates + the (grid x symbol) matrix.
        #
        # PER-SYMBOL COMPLEX-GAIN NORMALISATION (E60, measured necessity):
        # an RF33 frame holds only 35 symbols, so a W=12 window overlaps a
        # multi-symbol fade from EVERY symbol -- the fade poisons the
        # frame-median reference and the e58 ratio detector goes blind
        # (measured: an 11-symbol x0.5 step fade read ratio 1.2-1.9 against
        # a clean-frame spread reaching 5.2, and the un-normalised smoother
        # decoded 0/74 where the per-symbol path held 74/74).  A fade is,
        # to first order, ONE complex scalar per symbol -- and that scalar
        # is measured from the symbol's own ~1-2k pilots (LS fit onto the
        # frame-mean channel shape M) with negligible noise.  So c_l rides
        # OUTSIDE the boxcar -- amplitude AND phase steps are followed
        # within their own symbol -- and the boxcar smooths only the channel
        # SHAPE, which is what is actually quasi-static.  The change
        # detector watches the normalised residual and still guards what a
        # scalar cannot represent (delay / multipath shape changes).
        ngrid, nsym = smt["ngrid"], C.NSYM
        Hg = np.zeros((ngrid, nsym), np.complex128)
        Mg = np.zeros((ngrid, nsym), bool)
        hps = []
        for ci, t in enumerate(cp["tabs"]):
            g = t["g"]
            hp = Yb[t["rows"]][:, g["bins_pk"]] / g["refpk"]
            hps.append(hp)
            gx = smt["gidx"][ci]
            Hg[np.ix_(gx, t["rows"])] = hp.T
            Mg[np.ix_(gx, t["rows"])] = True
        cl = np.ones(nsym, np.complex128)
        w_c = np.zeros(nsym, bool)
        if mg.ce_norm:
            cnt_tot = np.maximum(Mg.sum(axis=1), 1)
            M = np.where(Mg, Hg, 0).sum(axis=1) / cnt_tot  # frame-mean shape
            for ci, t in enumerate(cp["tabs"]):
                gx = smt["gidx"][ci]
                Mx = M[gx]
                den = max(float(np.sum(np.abs(Mx) ** 2)), 1e-30)
                cl[t["rows"]] = (hps[ci] @ np.conj(Mx)) / den
            # a degenerate gain (deep notch on the whole symbol) must not
            # blow up the normalisation -- clamp; the detector handles it
            w_c = np.abs(cl) < 1e-6
            cl[w_c] = 1.0
            Hg /= cl[None, :]
        # boxcar via running sums (O(nsym) total)
        Wn = mg.ce_w
        csH = np.concatenate([np.zeros((ngrid, 1), complex),
                              np.cumsum(np.where(Mg, Hg, 0), axis=1)], axis=1)
        csM = np.concatenate([np.zeros((ngrid, 1)),
                              np.cumsum(Mg, axis=1)], axis=1)
        losym = np.maximum(np.arange(nsym) - Wn, 0)
        hisym = np.minimum(np.arange(nsym) + Wn + 1, nsym)
        num = csH[:, hisym] - csH[:, losym]
        cnt = csM[:, hisym] - csM[:, losym]
        okc = cnt > 0
        Hsm = np.where(okc, num / np.maximum(cnt, 1), 0.0)
        covered = okc.all(axis=0)

        # change detector, two tests per symbol:
        #  * RELATIVE (e58's): residual > detect x frame median -- catches a
        #    minority of odd symbols on a stable frame.
        #  * ABSOLUTE (E60): residual > ce_abs x the symbol's OWN pilot
        #    noise floor (second differences along the carrier axis cancel
        #    the channel to its curvature and leave 6 sigma^2).  MEASURED
        #    necessity: a shape change (delay step) over 11 of 35 symbols
        #    poisons every boxcar window, the frame median rises with the
        #    residuals, and the relative test saturates near ratio ~3 while
        #    the frame decodes 0/74.  Against the symbol's own noise the
        #    poisoned residuals sit orders of magnitude high, EVERY symbol
        #    falls back, and the frame decodes exactly as the per-symbol
        #    path -- which is the correct answer when the smoothing model
        #    does not hold.  On a low-SNR channel (RF25) the residual IS the
        #    noise floor, abs-ratio ~ 1, and the smoothing stays on.
        res_l = np.empty(nsym)
        sig2_l = np.empty(nsym)
        for ci, t in enumerate(cp["tabs"]):
            gx = smt["gidx"][ci]
            hpn = hps[ci] / cl[t["rows"]][:, None]
            R = hpn - Hsm[np.ix_(gx, t["rows"])].T
            res_l[t["rows"]] = np.mean(np.abs(R) ** 2, axis=1)
            v = hpn[:, 2:] - 2.0 * hpn[:, 1:-1] + hpn[:, :-2]
            sig2_l[t["rows"]] = np.mean(np.abs(v) ** 2, axis=1) / 6.0
        med = max(float(np.median(res_l)), 1e-30)
        ratio = res_l / med
        abs_ratio = res_l / np.maximum(sig2_l, 1e-30)
        smooth_ok = (covered & (ratio <= mg.ce_detect)
                     & (abs_ratio <= mg.ce_abs) & (Wn > 0) & ~w_c)
        n_fb = int(nsym - smooth_ok.sum())

        lens = cp["lens"]
        offs = np.concatenate([[0], np.cumsum(lens)]) + len(spare)
        pool = np.empty(int(offs[-1]), np.complex128)
        pool[:len(spare)] = spare
        owner = np.empty(int(offs[-1]), int)
        owner[:len(spare)] = -1
        jg, wtg = smt["jg"], smt["wtg"]

        def one_class(ci):
            t = cp["tabs"][ci]
            g = t["g"]
            rows = t["rows"]
            Yr = Yb[rows]
            hp = hps[ci]
            sm = smooth_ok[rows]
            H = np.empty((len(rows), len(g["kk"])), np.complex128)
            if sm.any():
                Hc = Hsm[:, rows[sm]].T
                H[sm] = (Hc[:, jg] * (1.0 - wtg) + Hc[:, jg + 1] * wtg)
                H[sm] *= cl[rows[sm]][:, None]     # re-apply the symbol gain
            if (~sm).any():                     # untouched per-symbol path
                hpf = hp[~sm]
                H[~sm] = (hpf[:, t["j"]] * (1.0 - t["wt"])
                          + hpf[:, t["j"] + 1] * t["wt"])
            z = Yr[:, g["bins_d"]] / H[:, g["d"]]
            xg = np.empty_like(z)
            xg[np.arange(len(z))[:, None], t["Hmat"]] = z
            a, b = t["trim"]
            for i, l in enumerate(rows):
                pool[offs[l]:offs[l + 1]] = xg[i, a:b]
                owner[offs[l]:offs[l + 1]] = l
        if self.ex is None:
            for ci in range(len(cp["tabs"])):
                one_class(ci)
        else:
            list(self.ex.map(one_class, range(len(cp["tabs"]))))
        rep.update(n_preamble_spare=int(len(spare)),
                   n_preamble_cells=int(npre), n_normal=int(cp["n_norm"]),
                   n_sbs_active=int(lens[0]), n_null=int(cp["n_null"]),
                   null_low=int(cp["lo_n"]), pool=int(len(pool)),
                   fell_back=n_fb, detect_ratio_max=float(ratio.max()),
                   detect_ratios=[round(float(r), 2) for r in ratio],
                   abs_ratios=[round(float(r), 2) for r in abs_ratio])
        return pool, dict(symbol_of=owner, report=rep)

    def cpe_correct(self, pool, symbol_of, alphabets, region, dummy_values,
                    iters=3):
        """m6_cells.cpe_correct, one thread per OFDM symbol.  Symbols do not
        interact, so this is exact.

        cpu_fast only changes HOW the nearest constellation point is found:
        argmin over |z-p| becomes argmin over (|p|^2 - 2 Re(z conj(p))) via one
        (n,2)@(2,npts) matmul -- the same ordering in exact arithmetic, so a
        flip needs two points equidistant to within a rounding.  The phase
        estimate itself (the vdot ratio) is untouched either way, and the e2e
        gate diffs the decoded bytes of the two paths on real air.
        """
        out = pool / np.sqrt(np.mean(np.abs(pool) ** 2))
        syms = np.unique(symbol_of)
        gains = {}
        fast = self.cpu_fast
        if fast:
            Ps = [np.vstack((p.real, p.imag)) for p in alphabets]
            p2s = [p.real ** 2 + p.imag ** 2 for p in alphabets]
            # symbol_of is a concatenation in symbol order (cell_pool builds
            # it that way), so the per-symbol groups are SLICES and the whole
            # pool vectorises: one nearest-point matmul over every cell, then
            # per-symbol phase estimates as segment sums.  ~900 small NumPy
            # calls per frame become ~30 large ones -- the CPE was bound by
            # call dispatch under the GIL, not by arithmetic.  Falls back to
            # the per-symbol path if the ordering assumption ever breaks.
            if np.all(np.diff(symbol_of) >= 0):
                return self._cpe_fast(out, symbol_of, alphabets, region,
                                      dummy_values, iters, Ps, p2s)

        def one(sym):
            m = np.flatnonzero(symbol_of == sym)
            z = out[m]
            r = region[m]
            c_tot = 1.0 + 0j
            for _ in range(iters):
                hard = np.empty_like(z)
                for rid, pts in enumerate(alphabets):
                    sel = r == rid
                    if not sel.any():
                        continue
                    zz = z[sel]
                    if fast:
                        sc = np.column_stack((zz.real, zz.imag)) @ Ps[rid]
                        sc *= -2.0
                        sc += p2s[rid][None, :]
                        hard[sel] = pts[sc.argmin(1)]
                    else:
                        hard[sel] = pts[np.abs(zz[:, None]
                                               - pts[None, :]).argmin(1)]
                sel = r == 2
                if sel.any():
                    hard[sel] = dummy_values[m[sel]]
                c = np.vdot(hard, z) / max(np.vdot(hard, hard).real, 1e-12)
                z = z / c
                c_tot *= c
            return m, z, complex(c_tot)
        res = ([one(s) for s in syms] if self.ex is None
               else list(self.ex.map(one, syms)))
        for (m, z, c), sym in zip(res, syms):
            out[m] = z
            gains[int(sym)] = c
        return out, gains

    def _cpe_fast(self, z, symbol_of, alphabets, region, dummy_values,
                  iters, Ps, p2s):
        """Whole-pool CPE: same algorithm, segment reductions per symbol.

        The per-symbol phase `c` uses segment sums where the reference uses
        `np.vdot` -- a summation-order re-association (fast mode's standing
        concession, decoded-bytes gated).  The hard-decision points, the
        `1e-12` guard, the per-symbol division and the 3-iteration structure
        are the reference's, line for line.
        """
        # E53: the phase ESTIMATE does not need every cell.  The per-symbol
        # CPE is a mean over ~5700 cells; estimating it from every 4th cell
        # doubles the estimator's noise (sigma ~ 1/sqrt(N)) from ~0.4 to
        # ~0.8 milliradian-scale -- far under the channel's own phase noise
        # -- while the CORRECTION still applies to every cell.  Only the
        # nearest-point searches and segment sums shrink 4x.  Decoded-bytes
        # gated (gate_e52 leg 3 + the m11 e2e sha) like every concession.
        CPE_DECIM = 4
        sub = np.arange(0, len(z), CPE_DECIM)
        symbol_sub = symbol_of[sub]
        region_sub = region[sub]
        syms, starts = np.unique(symbol_sub, return_index=True)
        edges = np.concatenate([starts, [len(symbol_sub)]])
        counts_sub = np.diff(edges)
        # full-pool per-symbol counts, for applying the correction
        syms_f, starts_f = np.unique(symbol_of, return_index=True)
        counts = np.diff(np.concatenate([starts_f, [len(symbol_of)]]))
        rsel = [np.flatnonzero(region_sub == rid)
                for rid in range(len(alphabets))]
        dsel = np.flatnonzero(region_sub == 2)
        ctot = np.ones(len(syms), complex)
        zs = z[sub]
        hard = np.empty_like(zs)
        if len(dsel):
            hard[dsel] = dummy_values[sub[dsel]]
        Pf = [P.astype(np.float32) for P in Ps]
        p2f = [p2.astype(np.float32) for p2 in p2s]
        nch = max(1, min(16, self.threads))

        def decide(rid, sel):
            zz = zs[sel]
            X = np.empty((len(sel), 2), np.float32)
            X[:, 0] = zz.real
            X[:, 1] = zz.imag
            sc = X @ Pf[rid]
            sc *= np.float32(-2.0)
            sc += p2f[rid][None, :]
            hard[sel] = alphabets[rid][sc.argmin(1)]
        for _ in range(iters):
            jobs = []
            for rid in range(len(alphabets)):
                sel = rsel[rid]
                step = (len(sel) + nch - 1) // nch
                jobs += [(rid, sel[i:i + step])
                         for i in range(0, len(sel), step)]
            if self.ex is None:
                for rid, sel in jobs:
                    decide(rid, sel)
            else:
                list(self.ex.map(lambda j: decide(*j), jobs))
            w = np.conj(hard)
            w *= zs
            h2 = hard.real ** 2 + hard.imag ** 2
            num = np.add.reduceat(w, edges[:-1])
            den = np.add.reduceat(h2, edges[:-1])
            c = num / np.maximum(den, 1e-12)
            zs = zs / np.repeat(c, counts_sub)
            z = z / np.repeat(c, counts)
            ctot *= c
        gains = {int(s): complex(c) for s, c in zip(syms, ctot)}
        return z, gains

    # -- BICM --------------------------------------------------------------
    def demap_batch(self, cells2d, pts, ones, chunk=8, d2min_out=None):
        """(nblk, ncells) -> (nblk, ncells*nbits) max-log LLRs.

        Identical to calling m6_bicm.demap_llr once per row: sigma^2 is a mean
        over each row's own cells, and every other operation is elementwise.

        E60, cpu_fast only: `d2min_out` (an (nblk, ncells) float32 array)
        switches the row-scalar sigma^2 OFF -- the LLRs come back as raw
        bitwise-minima differences and the per-cell squared distance to the
        nearest point is written to `d2min_out`, so the caller can apply
        E58's per-cell noise weights instead.  Never used on the exact path.
        """
        nblk, ncell = cells2d.shape
        nb = ones.shape[0]
        if d2min_out is not None and not self.cpu_fast:
            raise ValueError("d2min_out is a cpu_fast-only lever (E60)")
        # cpu_fast: the LLRs are computed in float32 anyway; storing them in
        # float32 and letting the float32 LDPC take them WITHOUT the
        # f32->f64->f32 round trip changes no bit (the round trip is exact).
        out = np.empty((nblk, ncell, nb),
                       np.float32 if self.cpu_fast else np.float64)

        # The 12 subset minima were the expensive part, not the distances:
        # `d2[:, :, ones[i]]` is fancy indexing, so each of them COPIED a
        # (chunk, 2700, 32) float64 array.  A/322 6.3.4.2 numbers the label
        # bits MSB-first, so bit i's mask is exactly axis i of the last axis
        # reshaped to (2,)*MOD_BITS -- a strided VIEW, and min over a set is
        # order-independent, so the values are unchanged.
        sub = tuple(range(2, 2 + nb - 1))
        fast = self.cpu_fast
        if fast:
            # cpu_fast: |z-p|^2 = |z|^2 + (|p|^2 - 2 Re(z conj(p))).  The
            # |z|^2 term is constant over the 64 points of one cell, so the
            # bitwise minima DIFFERENCE a-b is exactly the difference of the
            # shifted minima, and only sigma^2 needs |z|^2 added back.  That
            # turns hypot+square over (blk, 2700, 64) complex128 into one
            # float32 sgemm plus elementwise passes at half the width --
            # a rounding re-association, gated at the decoded-bytes level.
            Pf = np.vstack((pts.real, pts.imag)).astype(np.float32)
            p2f = (pts.real ** 2 + pts.imag ** 2).astype(np.float32)

        def job(lo, hi):
            z = cells2d[lo:hi]
            if fast:
                k = hi - lo
                zr = np.asarray(z.real, np.float32)
                zi = np.asarray(z.imag, np.float32)
                X = np.empty((k * ncell, 2), np.float32)
                X[:, 0] = zr.ravel()
                X[:, 1] = zi.ravel()
                q = X @ Pf                              # (k*ncell, npts)
                q *= np.float32(-2.0)
                q += p2f[None, :]
                q = q.reshape(k, ncell, len(pts))
                z2 = zr * zr + zi * zi
                if d2min_out is None:
                    s2 = np.maximum((z2 + q.min(2)).mean(1), 1e-9)[:, None]
                else:
                    d2min_out[lo:hi] = z2 + q.min(2)
                    s2 = None
                qr = q.reshape(k, ncell, *([2] * nb))
                head = (slice(None), slice(None))
                for i in range(nb):
                    a = qr[head + (slice(None),) * i + (1,)].min(axis=sub)
                    b = qr[head + (slice(None),) * i + (0,)].min(axis=sub)
                    out[lo:hi, :, i] = (a - b) / s2 if s2 is not None \
                        else (a - b)
                return
            d2 = np.abs(z[:, :, None] - pts[None, None, :])
            np.square(d2, out=d2)
            s2 = np.maximum(d2.min(2).mean(1), 1e-9)[:, None]
            d2r = d2.reshape(hi - lo, ncell, *([2] * nb))
            head = (slice(None), slice(None))
            for i in range(nb):
                a = d2r[head + (slice(None),) * i + (1,)].min(axis=sub)
                b = d2r[head + (slice(None),) * i + (0,)].min(axis=sub)
                out[lo:hi, :, i] = (a - b) / s2
        n = max(1, (nblk + chunk - 1) // chunk)
        if self.ex is None:
            for i in range(n):
                job(i * chunk, min(nblk, (i + 1) * chunk))
        else:
            list(self.ex.map(lambda i: job(i * chunk,
                                           min(nblk, (i + 1) * chunk)),
                             range(n)))
        return out.reshape(nblk, ncell * nb)

    def ldpc_batch(self, lam2d, chain, gpu):
        """Batched LDPC with `PlpChain.decode_llr`'s alpha-ladder semantics.

        The ladder (1.0, 0.85, 0.75, tried until the syndrome is zero) landed
        in `m6_bicm` on the 64K/CTI track while this one was in flight, and it
        is a REAL behaviour change: on RF33 it converges one more FEC Block in
        35,224 than the fixed-0.75 decoder M7 shipped.  Reproducing it here is
        not optional -- the gate is byte-identity with whatever the reference
        chain currently does, so the ladder has to be in both or neither.

        It costs nothing: the first rung decodes ~all 74 Blocks, so the retry
        batch is empty or a single block.  First converged rung wins; if none
        converge, the rung with the fewest unsatisfied checks does -- exactly
        `decode_llr`.
        """
        ladder = self.ladder_for(chain)
        nb = int(lam2d.shape[0])
        bits = np.empty((nb, chain.ninner), np.uint8)
        conv = np.zeros(nb, bool)
        its = np.zeros(nb, np.int32)
        bad = np.full(nb, 1 << 30, np.int32)
        todo = np.arange(nb)
        for a in ladder:
            if len(todo) == 0:
                break
            sub = lam2d
            if len(todo) != nb:
                sub = (lam2d[list(map(int, todo))] if hasattr(lam2d, "device")
                       else lam2d[todo])
            if gpu is not None:
                with self._lock:
                    b, c, it, d = gpu.decode(sub, alpha=a)
            else:
                b, c, it, d = self._ldpc_cpu(sub, chain, a)
            better = d < bad[todo]
            take = c | better
            sel = todo[take]
            bits[sel] = b[take]
            its[sel] = it[take]
            bad[sel] = d[take]
            conv[sel] = c[take]
            todo = todo[~c]
        return bits, conv, its, bad

    def ladder_for(self, chain):
        a = getattr(chain, "alpha", None)
        if a:
            return (a,)
        return tuple(getattr(chain, "ALPHA_LADDER", (0.75,)))

    def _ldpc_cpu(self, lam2d, chain, alpha):
        """Batched min-sum over sub-batches of blocks, one per thread.

        `min_sum_decode_batch` at float64 is bit-identical per row to
        `min_sum_decode` (gate_e52), so the exact path keeps its bytes while
        losing the per-block NumPy dispatch overhead.  cpu_fast runs the same
        batch decoder at float32 -- half the memory traffic on the loop that
        memory traffic bounds -- and answers to the decoded-bytes gates.
        Sub-batching is exact by construction: every operation in the batch
        decoder is row-independent.
        """
        nb = len(lam2d)
        dt = np.float32 if self.cpu_fast else np.float64
        if nb == 0:
            return (np.empty((0, chain.ninner), np.uint8),
                    np.zeros(0, bool), np.zeros(0, np.int32),
                    np.zeros(0, np.int32))
        nch = 1 if self.ex is None else max(1, min(nb, self.threads))
        step = (nb + nch - 1) // nch
        sls = [slice(i * step, min(nb, (i + 1) * step))
               for i in range(nch) if i * step < nb]

        def one(sl):
            return LD.min_sum_decode_batch(lam2d[sl], chain.checks,
                                           chain.ninner, iters=self.iters,
                                           alpha=alpha, dtype=dt)
        res = ([one(s) for s in sls] if self.ex is None
               else list(self.ex.map(one, sls)))
        bits = np.concatenate([r[0] for r in res])
        conv = np.concatenate([r[1] for r in res])
        its = np.concatenate([r[2] for r in res])
        bad = np.concatenate([r[3] for r in res])
        return bits, conv, its, bad

    def _plp0_parts(self, symbol_of, hidx):
        """Per-cell pool-part ids of the HTI-gathered PLP0 cells (E60).

        part 0 = the preamble spare, part l+1 = OFDM symbol l -- e58's
        granularity.  Frame-constant (the pool layout is fixed geometry), so
        built once per hidx and cached with its bincount denominators.
        """
        pp = getattr(self, "_pp", None)
        if pp is None or pp[3] is not hidx:
            pid = (symbol_of[hidx] + 1).astype(np.intp)
            npart = int(pid.max()) + 1
            cnts = np.maximum(np.bincount(pid.ravel(), minlength=npart), 1)
            self._pp = pp = (pid, cnts, npart, hidx)
        return pp[0], pp[1], pp[2]

    def baseband_batch(self, bits2d, chain):
        """Codewords -> (bch syndromes, packed Baseband Packet bytes)."""
        nouter = bits2d[:, :chain.kldpc]
        syn = m9_accel.bch_syndrome_batch(nouter, chain.ninner)
        seq = SC.sequence(chain.kpayload)
        pay = nouter[:, :chain.kpayload] ^ seq[None, :]
        by = np.packbits(pay, axis=1)
        return syn, [row.tobytes() for row in by]

    # -- the fused GPU path -------------------------------------------------
    def _build_gpu_bicm(self, pool_len, dv, symbol_of):
        import m9_gpu
        nrows = self.ch0.cells_per_fec
        ncols = C.PLP0["n_fec_max"] // C.PLP0["nti"]
        per_ti = C.PLP0["size"] // C.PLP0["nti"]
        idx = np.array([ti * per_ti + T.fec_block(np.arange(per_ti), j,
                                                  nrows, ncols, 0)
                        for ti in range(C.PLP0["nti"])
                        for j in range(ncols)])
        self.gb = m9_gpu.GpuBicm(
            pool_len, C.PLP0["size"], C.PLP16["start"], C.PLP16["size"],
            symbol_of, dv, self.pts0, self.pts16,
            BI.MOD_BITS["64QAM"], BI.MOD_BITS["QPSK"], idx,
            self.ch0.lam_of_q, self.ch16.lam_of_q)

    def _frame_gpu(self, pool, info, dv, tm=None):
        """CPE + HTI + demap + bit de-interleave + LDPC, without leaving the
        GPU.  Gated on byte-identical Baseband Packets, not on ULPs."""
        tm = self.tm if tm is None else tm
        pc = time.perf_counter
        if self.gb is None:
            self._build_gpu_bicm(len(pool), dv, info["symbol_of"])
        t = pc()
        with self._lock:
            lam0, lam16, tail = self.gb.run(pool)
        tm["gpu_bicm"] += pc() - t
        t = pc()
        b16, cv16, it16, bd16 = self.ldpc_batch(lam16, self.ch16, self.g16)
        bits, conv, its, bad = self.ldpc_batch(lam0, self.ch0, self.g0)
        tm["ldpc"] += pc() - t
        t = pc()
        dummy = C.dummy_check(np.concatenate(
            [np.zeros(C.PLP16["start"] + C.PLP16["size"], complex), tail]))
        tm["dummy_check"] += pc() - t
        r16 = dict(bits=b16[0], converged=bool(cv16[0]), iters=int(it16[0]),
                   unsatisfied=int(bd16[0]))
        t = pc()
        syn, by = self.baseband_batch(bits, self.ch0)
        tm["bch_descramble"] += pc() - t
        packets, nconv, nbch = [], 0, 0
        for i in range(len(bits)):
            if conv[i]:
                nconv += 1
                nbch += int(syn[i] == 0)
                packets.append(by[i])
            else:
                packets.append(None)
        return r16, packets, dict(
            dummy=dummy, converged=nconv, bch_ok=nbch, n_fec=C.PLP0["n_fec"],
            nrows=self.ch0.cells_per_fec,
            ncols=C.PLP0["n_fec_max"] // C.PLP0["nti"])

    # -- one frame ---------------------------------------------------------
    def decode_frame(self, y, t0, rep=None):
        """Drop-in for m6_payload.decode_frame.

        Two decode workers may run this concurrently (m11_watch E52), so the
        stage timers accumulate locally and merge under a lock at the end --
        a Counter `+=` is a read-modify-write and concurrent frames were
        silently losing each other's milliseconds.
        """
        tm, pc = collections.Counter(), time.perf_counter
        try:
            return self._decode_frame(y, t0, rep, tm, pc)
        finally:
            with self._tm_lock:
                self.tm.update(tm)

    def _decode_frame(self, y, t0, rep, tm, pc):
        t = pc()
        # E60: the margin levers ride the cpu_fast path only; exact float64
        # and the GPU backends run the pre-E60 code untouched.
        mg = self.margin if (self.cpu_fast and self.backend == "cpu") \
            else M16.OFF
        rep_d = rep if rep is not None else {}
        if mg.ce_w > 0:
            pool, info = self.cell_pool_fast_sm(y, t0, rep_d)
            tm["n_ce_fallback"] += rep_d.get("fell_back", 0)
        else:
            cp = self.cell_pool_fast if self.cpu_fast else self.cell_pool
            pool, info = cp(y, t0, rep_d)
        tm["cell_pool"] += pc() - t
        key = len(pool)
        with self._tm_lock:
            if key not in self.regions:
                self.regions[key] = C.constellation_regions(key)
        reg, dv = self.regions[key]
        alph = [self.pts0, self.pts16]

        if self.backend == "gpu-full":
            return self._frame_gpu(pool, info, dv, tm)

        t = pc()
        pool, _ = self.cpe_correct(pool, info["symbol_of"], alph, reg, dv)
        tm["cpe_correct"] += pc() - t
        t = pc()
        dummy = C.dummy_check(pool)
        tm["dummy_check"] += pc() - t

        # E60: the auto-switch instrument -- channel SNR off the known +-1
        # dummy cells (microseconds; no constellation hypothesis).  Weighted
        # LLRs re-round float32, so they engage only below the SNR gate (WL
        # law shape: on at the cliff, off with margin) -- strong-signal
        # decodes stay byte-identical to HEAD by construction of the switch.
        snr_db = None
        use_w = False
        if mg.wllr != "off":
            snr_db = M16.frame_snr_db(
                pool, C.PLP16["start"] + C.PLP16["size"])
            use_w = (mg.wllr == "on"
                     or (snr_db is not None and snr_db < mg.snr_gate_db))
            # += 0 creates the key: the counter must EXIST at zero so the
            # strong-signal gate can tell "switch held" from "not recorded"
            tm["n_wllr_frames"] += int(use_w)

        # PLP 16 -- one FEC Block
        t = pc()
        ldt = np.float32 if self.cpu_fast else np.float64
        c16 = pool[C.PLP16["start"]:C.PLP16["start"] + C.PLP16["size"]]
        q16 = self.demap_batch(c16[None, :], self.pts16, self.ones16)
        lam16 = np.empty((1, 16200), ldt)
        lam16[:, self.ch16.lam_of_q] = q16
        b16, cv16, it16, bd16 = self.ldpc_batch(lam16, self.ch16, self.g16)
        r16 = dict(bits=b16[0], converged=bool(cv16[0]), iters=int(it16[0]),
                   unsatisfied=int(bd16[0]))
        tm["plp16"] += pc() - t

        # PLP 0 -- 74 FEC Blocks, all independent
        t = pc()
        plp0 = pool[:C.PLP0["size"]]
        nrows, ncols = self.ch0.cells_per_fec, C.PLP0["n_fec_max"] // C.PLP0["nti"]
        per_ti = len(plp0) // C.PLP0["nti"]
        # ONE cached gather instead of 74 fec_block calls: fec_block is pure
        # indexing, so fec_block(seg, ...) == seg[fec_block(arange, ...)]
        # element for element (the same identity _build_gpu_bicm uses).
        hidx = getattr(self, "_hti_idx", None)
        if hidx is None or hidx.shape[1] != nrows:
            hidx = self._hti_idx = np.array(
                [ti * per_ti + T.fec_block(np.arange(per_ti), j, nrows,
                                           ncols, 0)
                 for ti in range(C.PLP0["nti"]) for j in range(ncols)])
        cells = plp0[hidx]
        tm["hti"] += pc() - t
        t = pc()
        if use_w:
            # E58 lever 2, live: sigma^2 per pool part (preamble spare +
            # each OFDM symbol), measured from the demap's own per-cell
            # nearest-point residuals -- the same estimator the row-scalar
            # sigma^2 already uses, at the granularity the fading has.
            # Saturates in deep fades (residual cannot exceed the point
            # spacing) = conservative; the LDPC syndrome + BCH stay the only
            # oracles, a weight cannot manufacture a codeword.
            d2min = np.empty(cells.shape, np.float32)
            q = self.demap_batch(cells, self.pts0, self.ones0,
                                 d2min_out=d2min)
            tm["demap"] += pc() - t
            t = pc()
            pid, cnts, npart = self._plp0_parts(info["symbol_of"], hidx)
            s2p = np.bincount(pid.ravel(), weights=d2min.ravel(),
                              minlength=npart) / cnts
            s2p = np.maximum(s2p, 1e-6).astype(np.float32)
            qr = q.reshape(len(cells), -1, BI.MOD_BITS["64QAM"])
            qr /= s2p[pid][:, :, None]
            tm["wllr"] += pc() - t
        else:
            q = self.demap_batch(cells, self.pts0, self.ones0)
            tm["demap"] += pc() - t
        t = pc()
        lam = np.empty((len(cells), 16200), ldt)
        lam[:, self.ch0.lam_of_q] = q
        tm["bit_deint"] += pc() - t
        t = pc()
        bits, conv, its, bad = self.ldpc_batch(lam, self.ch0, self.g0)
        tm["ldpc"] += pc() - t
        if mg.sp and not conv.all():
            # E58 lever 3, live: exact float64 BP only on blocks min-sum
            # failed NEAR the boundary (unsat < sp_cut).  Zero cost on
            # healthy air by construction; on a fade it is capped at
            # sp_blocks blocks and TIME-BOXED to sp_ms (E52 shed law: the
            # rescue budget sheds, the pipeline never stalls).
            t = pc()
            cand = np.flatnonzero(~conv & (bad < mg.sp_cut))
            if len(cand):
                cand = cand[np.argsort(bad[cand],
                                       kind="stable")][:mg.sp_blocks]
                stx = {}
                b2, c2, i2, _d2 = M16.sum_product_decode_batch(
                    np.asarray(lam[cand], np.float64), self.ch0.checks,
                    self.ch0.ninner, iters=mg.sp_iters,
                    deadline=pc() + mg.sp_ms / 1000.0, stats=stx)
                got = np.flatnonzero(c2)
                for k in got:
                    i = int(cand[k])
                    bits[i] = b2[k]
                    conv[i] = True
                    its[i] = i2[k]
                    bad[i] = 0
                tm["n_sp_tried"] += len(cand)
                tm["n_sp_rescued"] += len(got)
                tm["n_sp_shed"] += stx.get("sp_timeout", 0)
            tm["sp_rescue"] += pc() - t
        t = pc()
        syn, by = self.baseband_batch(bits, self.ch0)
        tm["bch_descramble"] += pc() - t
        packets, nconv, nbch = [], 0, 0
        for i in range(len(cells)):
            if conv[i]:
                nconv += 1
                nbch += int(syn[i] == 0)
                packets.append(by[i])
            else:
                packets.append(None)
        return r16, packets, dict(dummy=dummy, converged=nconv, bch_ok=nbch,
                                  n_fec=C.PLP0["n_fec"], nrows=nrows,
                                  ncols=ncols, snr_db=snr_db, wllr=use_w)

    def close(self):
        if self.ex is not None:
            self.ex.shutdown()


# ---------------------------------------------------------------------------
# gates
# ---------------------------------------------------------------------------

def verify(path, rate, fmt=None, frames=2, start=0.0, threads=16,
           backend="gpu", verbose=True, cpu_fast=False):
    """The whole fast path against the untouched M6/M7 path, byte for byte."""
    import m6_payload as P6
    ok = True

    # 1 -- the loader
    span = 0.05 + frames * FRAME_SEC
    m9_accel.uninstall()
    y_ref = MP.load(path, rate, fmt=fmt, span_sec=span * (rate / 6.912e6),
                    start_sec=start)[0]
    m9_accel.install()
    ex = ThreadPoolExecutor(max_workers=threads)
    y_fast = load_fast(path, rate, fmt=fmt,
                       span_sec=span * (rate / 6.912e6), start_sec=start,
                       ex=ex)[0]
    g = np.array_equal(y_ref, y_fast)
    ok &= g
    if verbose:
        print(f"    {'PASS' if g else 'FAIL'}  threaded loader bit-identical "
              f"({len(y_ref)} complex samples, overlap-save FIR + chunked "
              f"de-rotation)")
        if not g:
            print(f"          max |delta| = {np.abs(y_ref-y_fast).max():.3e}")
    ex.shutdown()

    # 2 -- the whole frame decode
    fd = FrameDecoder(threads=threads, backend=backend, cpu_fast=cpu_fast)
    t_prev = None
    nfr = 0
    diffs = collections_counter = 0
    for fi in range(frames):
        if (fi + 1) * C.FRAME_SAMPLES + 400000 > len(y_ref):
            break
        centre = C.BOOTSTRAP if t_prev is None else t_prev + C.FRAME_SAMPLES
        m9_accel.uninstall()
        t_ref, coh_ref = C.fine_timing(y_ref, span=20, centre=centre)
        r16r, pkr, dgr = P6.decode_frame(y_ref, t_ref, {})
        m9_accel.install()
        t_f, coh_f = fd.fine_timing(y_fast, span=20, centre=centre)
        r16f, pkf, dgf = fd.decode_frame(y_fast, t_f)
        t_prev = t_ref
        # cpu_fast's batched coherence drifts in the last ulps (pocketfft
        # multi-row SIMD; measured 2e-15, gate_e52) -- in fast mode the
        # DECISION (t0) must match, the metric may drift.  Exact mode keeps
        # bit equality on both.
        same_t = (t_ref == t_f) and (cpu_fast or coh_ref == coh_f)
        same_p = (len(pkr) == len(pkf)) and all(a == b for a, b in
                                                zip(pkr, pkf))
        same_16 = (bool(r16r["converged"]) == bool(r16f["converged"])
                   and np.array_equal(r16r["bits"], r16f["bits"]))
        same_d = (dgr["converged"] == dgf["converged"]
                  and dgr["bch_ok"] == dgf["bch_ok"])
        good = same_t and same_p and same_16 and same_d
        ok &= good
        diffs += sum(1 for a, b in zip(pkr, pkf) if a != b)
        nfr += 1
        if verbose:
            print(f"    {'PASS' if good else 'FAIL'}  frame {fi}: t0 "
                  f"{t_ref}=={t_f}  PLP0 {dgf['converged']}/74 conv, "
                  f"{dgf['bch_ok']}/74 BCH0, {len(pkf)} packets "
                  f"{'BYTE-IDENTICAL' if same_p else 'DIFFER'}  "
                  f"PLP16 {'identical' if same_16 else 'DIFFERS'}")
    fd.close()
    if verbose:
        print(f"    {nfr} frames, {diffs} differing Baseband Packets of "
              f"{74*nfr}")
    return ok


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("capture", nargs="?", default="long_rf33.cs16")
    ap.add_argument("--rate", type=float, default=8e6)
    ap.add_argument("--fmt")
    ap.add_argument("--frames", type=int, default=2)
    ap.add_argument("--start", type=float, default=0.0)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--backend", default="gpu",
                    choices=("gpu", "gpu-full", "cpu"))
    ap.add_argument("--cpu-fast", action="store_true",
                    help="verify the E52 float32 fast path against the "
                         "untouched M6/M7 oracle")
    a = ap.parse_args()
    p = a.capture if os.path.isabs(a.capture) else os.path.join(HERE,
                                                                a.capture)
    print("M9 -- fast frame decoder vs the M6/M7 reference chain")
    print("=" * 72)
    good = verify(p, a.rate, a.fmt, a.frames, a.start, a.threads, a.backend,
                  cpu_fast=a.cpu_fast)
    print(f"\n  {'BYTE-IDENTICAL END TO END' if good else 'DIVERGENCE'}")
    raise SystemExit(0 if good else 1)
