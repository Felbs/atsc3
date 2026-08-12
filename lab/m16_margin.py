#!/usr/bin/env python3
"""M16 -- the E58 margin-recovery levers as a SHARED module (E60 port).

E58 proved three levers offline on the banked RF25 Fox capture (+1.61 dB
stacked effective threshold, 176 -> 1056 blocks, all BCH-zero, gate_e58
7 legs).  They lived in lab/e58_*.py as opt-in instruments.  This module is
the port target: the SAME arithmetic, importable by the LIVE chain (m9_fast /
m11) and by the generic L1-driven path (m10 / the e49 RF25 chain), so both
inherit one implementation.

The e58_* instruments are NOT modified and remain the references; gate_e60
holds this module to numerical equality against them (smoothed collect
cell-for-cell, sigma table byte-for-byte, exact-BP verdict-for-verdict).

The three levers
----------------
1. `smoothed_pool_g`   -- time-smoothed 2x1D channel estimation with the
                          per-symbol change detector (e58_collect's math,
                          geometry-driven; W=0 is the identity).
2. `sigma_parts` /     -- per (frame x symbol-part) effective demap noise and
   `build_weights`        the per-cell LLR weights through the CTI
                          (e58_weighted's math).
3. `sum_product_decode_batch`
                       -- exact float64 log-tanh BP for the near-miss rescue,
                          e58_weighted's decoder plus an optional DEADLINE:
                          live is real time, so the rescue budget obeys the
                          E52 shed law (time-boxed, shed beyond, counted).
                          deadline=None reproduces e58 exactly.

`MarginCfg` is the live chain's one switch panel: every lever has a flag, an
env override, and a documented default.  ATSC3_MARGIN=0 turns the whole port
off and the chain runs the pre-E60 code paths untouched.
"""
from __future__ import annotations

import os
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
import sys
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m3_ldpc as LD                 # noqa: E402

SQ2 = np.sqrt(2.0)


# ---------------------------------------------------------------------------
# configuration -- flags, env overrides, defaults
# ---------------------------------------------------------------------------

class MarginCfg:
    """The E60 switch panel.  Defaults are the DEFAULT-ON-WHERE-SAFE choices:

      ce_w = 12       CE smoothing on for all channels (E58's window; the
                      change detector keeps real fades on the per-symbol
                      path, gate_e60 leg `fade`).  0 disables.
      wllr = "auto"   per-cell noise-weighted LLRs engage only when the
                      frame's dummy-cell SNR is BELOW snr_gate_db -- the WL
                      law shape: on at the cliff, off with margin -- because
                      re-scaling LLRs re-rounds float32 and strong-signal
                      byte identity is the gate.  "on"/"off" force it.
      snr_gate_db     the switch point, MEASURED both ways on the RF33 gate
                      capture: healthy air reads 17.1-20.0 dB dummy-SNR
                      (31 frames), and the AWGN-degraded cliff sits at
                      ~14.5-15.4 dB (74/74 at 15.4, 41/74 at 14.5, 0/74 at
                      13.4).  16.0 dB arms the lever ABOVE the cliff (WL
                      law: on at the cliff) while every healthy frame stays
                      >1.1 dB clear (strong-signal byte identity, gated).
      sp = True       exact-BP rescue of near-miss failures (unsat < sp_cut),
                      zero cost on healthy air by construction (only fires on
                      blocks min-sum already failed), time-boxed to sp_ms per
                      frame and capped at sp_blocks blocks per frame.
    """

    def __init__(self, ce_w=None, ce_detect=None, wllr=None, snr_gate_db=None,
                 sp=None, sp_ms=None, sp_cut=None, sp_iters=None,
                 sp_blocks=None, ce_norm=None):
        e = os.environ.get

        def num(v, key, dflt, cast):
            if v is not None:
                return cast(v)
            s = e(key)
            if s in (None, ""):
                return dflt
            try:
                return cast(s)
            except ValueError:
                return dflt
        master = e("ATSC3_MARGIN", "1") != "0"
        self.ce_w = num(ce_w, "ATSC3_CE_W", 12, int) if master else 0
        self.ce_detect = num(ce_detect, "ATSC3_CE_DETECT", 6.0, float)
        # per-symbol complex-gain normalisation (the RF33 short-frame
        # adaptation, gate_e60 leg 4).  ce_norm=0 is the e58-style absolute
        # grid -- kept ONLY as the gate's negative control.
        self.ce_norm = bool(int(num(ce_norm, "ATSC3_CE_NORM", 1, int)))
        # absolute detector: residual vs the symbol's own noise floor.
        # MEASURED (gate capture, complex-norm): clean frames sit at 1-4;
        # a whole-frame-poisoning shape change sits orders of magnitude up.
        self.ce_abs = num(None, "ATSC3_CE_ABS", 10.0, float)
        w = (wllr if wllr is not None else e("ATSC3_WLLR", "auto")).lower()
        self.wllr = w if w in ("auto", "on", "off") else "auto"
        if not master:
            self.wllr = "off"
        self.snr_gate_db = num(snr_gate_db, "ATSC3_SNR_GATE", 16.0, float)
        self.sp = (bool(int(num(sp, "ATSC3_SP", 1, int))) if master else False)
        self.sp_ms = num(sp_ms, "ATSC3_SP_MS", 120.0, float)
        self.sp_cut = num(sp_cut, "ATSC3_SP_CUT", 1500, int)
        self.sp_iters = num(sp_iters, "ATSC3_SP_ITERS", 100, int)
        self.sp_blocks = num(sp_blocks, "ATSC3_SP_BLOCKS", 8, int)

    def active(self):
        return self.ce_w > 0 or self.wllr != "off" or self.sp

    def asdict(self):
        return dict(ce_w=self.ce_w, ce_detect=self.ce_detect,
                    ce_norm=self.ce_norm, ce_abs=self.ce_abs, wllr=self.wllr,
                    snr_gate_db=self.snr_gate_db, sp=self.sp,
                    sp_ms=self.sp_ms, sp_cut=self.sp_cut,
                    sp_iters=self.sp_iters, sp_blocks=self.sp_blocks)


OFF = MarginCfg(ce_w=0, wllr="off", sp=0)


# ---------------------------------------------------------------------------
# lever 1 -- time-smoothed 2x1D channel estimation (generic geometry path)
# ---------------------------------------------------------------------------
# The arithmetic is e58_collect.smoothed_pool's, line for line: gate_e60
# holds this to cell-for-cell equality against the instrument on the Fox
# capture, and W=0 to equality with the untouched per-symbol path.

def smoothed_pool_g(y, t0, g, W=12, detect=6.0, report=None):
    """cell_pool_g with time-smoothed data-symbol channel estimation."""
    import m3_spec as S
    import m6_cells as C6
    import spec_pilots as P
    from m3_preamble import fft_bins
    rep = report if report is not None else {}
    fi0 = g.np_sym

    parts, owner, cohs = [], [], []
    pre = []
    for s in range(g.np_sym):
        z, coh = C6._demod_preamble_g(y, t0, g, s)
        pre.append(C6.FI.deinterleave(z, g.pre_nfft, s, direction="forward",
                                      toggle="i"))
        cohs.append(coh)
    pre = np.concatenate(pre) if pre else np.zeros(0, complex)
    spare = pre[g.l1b_cells + g.l1d_cells:]
    parts.append(spare)
    owner.append(np.full(len(spare), -1, int))

    noc = P.NOC[(g.nfft, g.cred)]
    lo, _ = S.carrier_abs_range(g.nfft, g.cred)
    ref = np.array(S.pilot_values(noc, 1.0))
    dx, dy = P.dxdy(g.pattern)
    ngrid = (noc - 1) // dx + 1                       # noc-1 is a multiple of dx
    assert (noc - 1) % dx == 0

    # ---- pass 1: FFT every data symbol, raw pilot estimates --------------
    Y = np.empty((g.nsym, g.nfft), np.complex128)
    pk_of, hp_of = [], []
    for l in range(g.nsym):
        w = (t0 + (g.pre_nfft + g.pre_gi) * g.np_sym + (g.nfft + g.gi) * l
             + g.gi)
        Y[l] = np.fft.fftshift(np.fft.fft(y[w:w + g.nfft]))
        sbs = g.is_sbs(l)
        pk = (np.arange(0, noc, dx) if sbs
              else np.arange(dx * (l % dy), noc, dx * dy))
        pk = np.unique(np.concatenate([pk, [0, noc - 1]]))
        hp = Y[l][fft_bins(g.nfft, lo + pk)] / ref[pk]
        pk_of.append(pk)
        hp_of.append(hp)

    # ---- pass 2: the (grid-carrier x symbol) matrix, boxcar in time ------
    Hg = np.zeros((ngrid, g.nsym), np.complex128)
    Mg = np.zeros((ngrid, g.nsym), bool)
    for l in range(g.nsym):
        j = pk_of[l] // dx
        on = pk_of[l] % dx == 0                       # edge noc-1 IS on grid
        Hg[j[on], l] = hp_of[l][on]
        Mg[j[on], l] = True
    # NaN-free boxcar via cumulative sums over the mask
    csH = np.concatenate([np.zeros((ngrid, 1), complex),
                          np.cumsum(np.where(Mg, Hg, 0), axis=1)], axis=1)
    csM = np.concatenate([np.zeros((ngrid, 1)),
                          np.cumsum(Mg, axis=1)], axis=1)
    losym = np.maximum(np.arange(g.nsym) - W, 0)
    hisym = np.minimum(np.arange(g.nsym) + W + 1, g.nsym)
    num = csH[:, hisym] - csH[:, losym]
    cnt = csM[:, hisym] - csM[:, losym]
    ok = cnt > 0
    Hsm = np.where(ok, num / np.maximum(cnt, 1), 0.0)

    # ---- pass 3: equalize each symbol off the smoothed dense grid --------
    dcoh, fell_back, ratios = [], 0, []
    grid_c = np.arange(ngrid) * dx
    res_l = np.empty(g.nsym)
    for l in range(g.nsym):
        sm = Hsm[:, l][ok[:, l]]
        gc = grid_c[ok[:, l]]
        Hf = C6._interp(np.arange(noc), gc, sm)
        res_l[l] = float(np.mean(np.abs(hp_of[l] - Hf[pk_of[l]]) ** 2))
    med = max(float(np.median(res_l)), 1e-30)

    for l in range(g.nsym):
        sbs = g.is_sbs(l)
        ratio = res_l[l] / med
        ratios.append(ratio)
        if W == 0 or ratio > detect:
            z, coh = C6._demod_data_g(y, t0, g, l)    # untouched path
            fell_back += W != 0
        else:
            sm = Hsm[:, l][ok[:, l]]
            gc = grid_c[ok[:, l]]
            Hf = C6._interp(np.arange(noc), gc, sm)
            used = np.zeros(noc, bool)
            used[np.array(sorted(P.pilot_carriers(g.nfft, g.cred, g.pattern,
                                                  l, sbs)), int)] = True
            d = np.flatnonzero(~used)
            z = Y[l][fft_bins(g.nfft, lo + d)] / Hf[d]
            hp = hp_of[l]
            coh = float(abs(np.sum(hp[:-1] * np.conj(hp[1:])))
                        / max(np.sum(np.abs(hp) ** 2), 1e-30))
        dcoh.append(coh)
        x = (C6.FI.deinterleave(z, g.nfft, l + fi0, direction="forward",
                                toggle="i") if g.freq_interleaver
             else np.asarray(z))
        if sbs:
            x = x[g.null_low:len(x) - g.null_high]
        parts.append(x)
        owner.append(np.full(len(x), l, int))

    pool = np.concatenate(parts)
    rep.update(pool=int(len(pool)), preamble_coherence=cohs,
               data_coherence_mean=float(np.mean(dcoh)),
               fell_back=int(fell_back),
               detect_ratio_max=float(np.max(ratios)))
    return pool, dict(symbol_of=np.concatenate(owner), report=rep)


# ---------------------------------------------------------------------------
# lever 2 -- per-cell noise-weighted LLRs (e58_weighted's math)
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


def sigma_parts(z, cum, p0_start, p0_size, min_cells=32):
    """One frame's per-symbol-part sigma^2 table row (e58 cmd_sigma's inner
    loop): mean squared distance of the PLP slice's cells to the nearest
    QPSK point, per pool part.  NaN where a part holds < min_cells cells."""
    nparts = len(cum) - 1
    row = np.full(nparts, np.nan, np.float32)
    d2 = qpsk_d2min(z)
    for k in range(nparts):
        lo = max(int(cum[k]) - p0_start, 0)
        hi = min(int(cum[k + 1]) - p0_start, p0_size)
        if hi - lo < min_cells:
            continue
        row[k] = max(float(d2[lo:hi].mean()), 1e-6)
    return row


def build_weights(jobs, ncell, nrows, start_row, per_frame, cum, p0_start,
                  table, permute=False):
    """Per-cell sigma^2 for each job's deinterleaved cell span (CTI path)."""
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


def frame_snr_db(pool, dummy_start):
    """Channel SNR off the known +-1 dummy cells (e49_dummy_snr's math).

    Free of any constellation hypothesis, ~3.6k cells on RF33 -- microseconds.
    This is the auto-switch's instrument.
    """
    import m4_scrambler as SC
    d = pool[dummy_start:]
    if len(d) < 64:
        return None
    s = 1.0 - 2.0 * SC.sequence(len(pool))[dummy_start:]
    a = float(np.dot(d.real, s) / np.dot(s, s))
    noise = float(np.mean(np.abs(d - a * s) ** 2))
    if noise <= 0 or a == 0.0:
        return None
    return float(10.0 * np.log10(a * a / noise))


# ---------------------------------------------------------------------------
# lever 3 -- exact sum-product on near-misses, deadline-capable
# ---------------------------------------------------------------------------

def sum_product_decode_batch(lam2d, checks, n, iters=100, deadline=None,
                             stats=None):
    """Exact BP (float64, log-tanh).  Returns (bits2d, conv, its, bad).

    deadline: perf_counter() value after which no NEW iteration is started
    (E52 shed law: live is real time, a fade must not stall the pipeline).
    Rows already converged keep their verdicts; the rest return their
    current hard decisions with conv=False and stats["sp_timeout"] += 1.
    deadline=None is exactly e58_weighted.sum_product_decode_batch.
    """
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
        if it == iters or (deadline is not None
                           and time.perf_counter() > deadline):
            if it < iters and stats is not None:
                stats["sp_timeout"] = stats.get("sp_timeout", 0) + len(act)
            bits_out[act] = tot[:, :n] < 0
            bad_out[act] = bad
            its_out[act] = it
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


def ladder_decode(lam2d, checks, n, iters, alphas=(1.0, 0.85, 0.75),
                  dtype=np.float32):
    """Batch alpha ladder (e58_weighted's, = m9_fast.ldpc_batch semantics)."""
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
                                             dtype=dtype)
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
