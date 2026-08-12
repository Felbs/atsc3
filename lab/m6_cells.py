#!/usr/bin/env python3
"""M6 -- build subframe 0's AVAILABLE-DATA-CELL pool for RF33, in cell order.

This is rung 1 of the payload chain and it is the one M5 closed.  Nothing
here is new physics; it is the M5 cell map turned from an arithmetic proof
into an actual array of cells.

The pool, exactly as `m5_cellmap.py` gated it (all four gates PASS, both
discrepancies ZERO, three frames):

    index      0 ..   3486   3487 cells inherited from the PREAMBLE symbol
                            (its data cells 1364..4850: 4851 total, minus
                             484 L1-Basic and 880 L1-Detail)
            3487 ..   8495   data symbol 0, an SBS symbol: 5136 TOTAL data
                             cells minus 63 null at the low band edge and
                             64 at the high edge = 5009 ACTIVE
            8496 .. 206462   data symbols 1..33, 5999 each
          206463 .. 211471   data symbol 34, SBS again, 5009 ACTIVE

    PLP  0 = pool[     0 : 199800]   64QAM-NUC 11/15, 74 FEC Blocks
    PLP 16 = pool[199800 : 207900]   QPSK 2/15, 1 FEC Block
    dummy  = pool[207900 : 211472]   3572 cells of known +-1  (A/322 7.2.6.5)

The dummy tail is not decoration: it is a per-run CHECK that the pool is
aligned, and it is checked here rather than assumed.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m3_freqint as FI                                          # noqa: E402
import m3_spec as S                                              # noqa: E402
import m4_scrambler as SC                                        # noqa: E402
import spec_pilots as P                                          # noqa: E402
from m2_pilots import load                                       # noqa: E402
from m3_preamble import fft_bins                                 # noqa: E402

BOOTSTRAP = 13824
NFFT, GI, PATTERN, CRED = 8192, 1536, "SP4_2", 0
NSYM = 35                       # data symbols in subframe 0
SBS_SYMS = (0, NSYM - 1)
PREAMBLE_CRED = 4               # A/322 7.2.5.1: minimum NoC for the FFT size
PREAMBLE_DX = 4
L1B_CELLS, L1D_CELLS = 484, 880

PLP0 = dict(start=0, size=199800, mod="64QAM", rate="11/15", ninner=16200,
            nti=2, n_fec=74, n_fec_max=74)
PLP16 = dict(start=199800, size=8100, mod="QPSK", rate="2/15", ninner=16200,
             nti=1, n_fec=1, n_fec_max=1)


# ---------------------------------------------------------------------------
# demodulation
# ---------------------------------------------------------------------------

def _interp(kk, pk, hp):
    return (np.interp(kk, pk, hp.real).astype(np.complex128)
            + 1j * np.interp(kk, pk, hp.imag))


def demod_data(y, t0, l, sbs):
    """Equalised data cells of DATA symbol l of subframe 0, in CARRIER order."""
    noc = P.NOC[(NFFT, CRED)]
    lo, _ = S.carrier_abs_range(NFFT, CRED)
    s = t0 + (NFFT + GI) * (1 + l) + GI
    Y = np.fft.fftshift(np.fft.fft(y[s:s + NFFT]))
    ref = np.array(S.pilot_values(noc, 1.0))
    dx, dy = P.dxdy(PATTERN)
    pk = (np.arange(0, noc, dx) if sbs
          else np.arange(dx * (l % dy), noc, dx * dy))
    pk = np.unique(np.concatenate([pk, [0, noc - 1]]))
    hp = Y[fft_bins(NFFT, lo + pk)] / ref[pk]
    H = _interp(np.arange(noc), pk, hp)
    used = np.zeros(noc, bool)
    used[np.array(sorted(P.pilot_carriers(NFFT, CRED, PATTERN, l, sbs)), int)] = True
    d = np.flatnonzero(~used)
    z = Y[fft_bins(NFFT, lo + d)] / H[d]
    coh = float(abs(np.sum(hp[:-1] * np.conj(hp[1:])))
                / max(np.sum(np.abs(hp) ** 2), 1e-30))
    return z, coh


def demod_preamble(y, t0):
    """Equalised data cells of the PREAMBLE symbol, in CARRIER order."""
    noc = S.noc(NFFT, PREAMBLE_CRED)
    lo, _ = S.carrier_abs_range(NFFT, PREAMBLE_CRED)
    Y = np.fft.fftshift(np.fft.fft(y[t0 + GI:t0 + GI + NFFT]))
    amp = S.PREAMBLE_PILOT_BOOST[(NFFT, GI)][1]
    ref = np.array(S.pilot_values(noc, amp))
    pk = np.arange(0, noc, PREAMBLE_DX)
    hp = Y[fft_bins(NFFT, lo + pk)] / ref[pk]
    H = _interp(np.arange(noc), pk, hp)
    used = np.zeros(noc, bool)
    used[pk] = True
    used[[c - lo for c in S.common_cps(NFFT, PREAMBLE_CRED)]] = True
    d = np.flatnonzero(~used)
    z = Y[fft_bins(NFFT, lo + d)] / H[d]
    return z, len(d)


# bootstrap + subframe 0 (36 x 9728) + subframe 1 (75 x 17920), at 6.912 Msps
FRAME_SAMPLES = 13824 + (8192 + 1536) * 36 + (16384 + 1536) * 75


def fine_timing(y, span=20, centre=None):
    """Scan the bootstrap-anchored window on a normal data symbol."""
    best = None
    centre = BOOTSTRAP if centre is None else centre
    for t in range(centre - span, centre + span + 1):
        _, coh = demod_data(y, t, 1, False)
        if best is None or coh > best[1]:
            best = (t, coh)
    return best


# ---------------------------------------------------------------------------
# the pool
# ---------------------------------------------------------------------------

def cell_pool(y, t0, report=None):
    """Subframe 0's available data cells, in cell order.  Returns (pool, info).

    `info["symbol_of"]` maps every pool index to its OFDM symbol (-1 for the
    preamble), which is what the per-symbol phase correction needs.
    """
    rep = report if report is not None else {}
    _, _, n_null = P.sbs_cells(NFFT, CRED, PATTERN, 1)
    lo_n, hi_n = P.sbs_null_split(n_null)
    n_norm = P.data_cells(NFFT, CRED, PATTERN, 1, False)

    zp, npre = demod_preamble(y, t0)
    xp = FI.deinterleave(zp, NFFT, 0, direction="forward", toggle="i")
    spare = xp[L1B_CELLS + L1D_CELLS:]
    parts, owner = [spare], [np.full(len(spare), -1, int)]

    for l in range(NSYM):
        sbs = l in SBS_SYMS
        z, _ = demod_data(y, t0, l, sbs)
        x = FI.deinterleave(z, NFFT, l + 1, direction="forward", toggle="i")
        if sbs:
            x = x[lo_n:len(x) - hi_n]
        parts.append(x)
        owner.append(np.full(len(x), l, int))

    pool = np.concatenate(parts)
    rep.update(n_preamble_spare=int(len(spare)), n_preamble_cells=int(npre),
               n_normal=int(n_norm), n_sbs_active=int(len(parts[1])),
               n_null=int(n_null), null_low=int(lo_n), pool=int(len(pool)))
    return pool, dict(symbol_of=np.concatenate(owner), report=rep)


def constellation_regions(pool_len, plp0=PLP0, plp16=PLP16):
    """Per-cell reference: which alphabet each pool cell is known to carry.

    Returns (region_id, dummy_values).  region 0 = PLP 0's 64QAM-NUC,
    1 = PLP 16's QPSK, 2 = A/322 7.2.6.5 dummy cells, whose values are
    KNOWN exactly (+-1 from the scrambler), not merely constrained.
    """
    reg = np.full(pool_len, 2, np.int8)
    reg[plp0["start"]:plp0["start"] + plp0["size"]] = 0
    reg[plp16["start"]:plp16["start"] + plp16["size"]] = 1
    dstart = plp16["start"] + plp16["size"]
    dv = np.zeros(pool_len, complex)
    if pool_len > dstart:
        dv[dstart:] = (1.0 - 2.0 * SC.sequence(pool_len)[dstart:pool_len])
    return reg, dv


def cpe_correct(pool, symbol_of, alphabets, region, dummy_values, iters=3):
    """Per-symbol residual complex gain, decision-directed.

    TWO REAL PARAMETERS PER OFDM SYMBOL, fitted from thousands of cells
    against the alphabet each cell is ALREADY known to carry (from the
    LDPC-decoded PLP allocation).  This is ordinary common-phase-error
    removal, not a search: it cannot rescue a wrong constellation
    hypothesis, and the LDPC remains the only oracle.  `--no-cpe` runs
    without it so its contribution is measured rather than assumed.
    """
    out = pool / np.sqrt(np.mean(np.abs(pool) ** 2))
    gains = {}
    for sym in np.unique(symbol_of):
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
                hard[sel] = pts[np.abs(zz[:, None] - pts[None, :]).argmin(1)]
            sel = r == 2
            if sel.any():
                hard[sel] = dummy_values[m[sel]]
            c = np.vdot(hard, z) / max(np.vdot(hard, hard).real, 1e-12)
            z = z / c
            c_tot *= c
        out[m] = z
        gains[int(sym)] = complex(c_tot)
    return out, gains


def dummy_check(pool, start=207900):
    """A/322 7.2.6.5: the tail cells must be +-1 with the scrambler's signs."""
    d = pool[start:]
    sgn = 1.0 - 2.0 * SC.sequence(len(d) + start)[start:start + len(d)]
    real = np.abs(d.imag) < 0.35
    agree = float(np.mean(np.sign(d.real[real]) == sgn[real]))
    return dict(n=int(len(d)), n_real=int(real.sum()),
                real_frac=float(real.mean()), sign_agreement=agree)


def load_frame(path, rate, fmt=None, start_sec=0.0, span_sec=0.26):
    return load(path, rate, fmt=fmt, span_sec=span_sec, start_sec=start_sec)[0]


# ===========================================================================
# M10 -- the same cell map, DRIVEN FROM L1 INSTEAD OF FROM THE CONSTANTS ABOVE
# ===========================================================================
#
# Everything above this line is RF33's geometry as module constants: one FFT
# size, one guard interval, one pilot pattern, two Subframes, a one-symbol
# Preamble, the frequency interleaver always on, and a hard-coded PLP split.
# That was fine for one multiplex and it is exactly what made the next one a
# rebuild (M8's finding).
#
# `Geometry` carries the same information but reads it out of the decoded
# L1-Basic / L1-Detail, so a new multiplex costs a JSON file rather than a
# module.  The legacy path is kept and `selftest_geometry()` requires the
# L1-driven path to reproduce it CELL FOR CELL on RF33 -- a refactor that
# changes a number is a bug, and this is how it would be caught.
#
# What is parameterised: FFT size, guard interval, Cred/NoC, pilot pattern,
# SBS first/last and the null-cell split, symbol count, the number of Preamble
# symbols and their two different NoCs (7.2.5.1), the L1-Basic/L1-Detail cell
# cost, the frequency interleaver ENABLE flag (RF30 bypasses it), and the PLP
# start/size/modulation/rate/Ninner/time-interleaver table.

import spec_l1syntax as _LS                                       # noqa: E402

_MOD_NAME = {0: "QPSK", 1: "16QAM", 2: "64QAM", 3: "256QAM",
             4: "1024QAM", 5: "4096QAM"}
_RATE_NAME = {i: "%d/15" % (i + 2) for i in range(12)}
_FFT_SIZE = {0: 8192, 1: 16384, 2: 32768}
_MOD_BITS = {"QPSK": 2, "16QAM": 4, "64QAM": 6, "256QAM": 8,
             "1024QAM": 10, "4096QAM": 12}


class Geometry:
    """One Subframe's cell map, derived from L1 (or from the RF33 constants)."""

    def __init__(self, *, nfft, gi, cred, pattern, nsym, sbs_first, sbs_last,
                 np_sym, pre_nfft, pre_gi, pre_dx, pre_cred, l1b_cells,
                 l1d_cells, freq_interleaver, plps, sbs_null_signalled=None,
                 label=""):
        self.nfft, self.gi, self.cred, self.pattern = nfft, gi, cred, pattern
        self.nsym = nsym
        self.sbs_first, self.sbs_last = bool(sbs_first), bool(sbs_last)
        self.np_sym = np_sym
        self.pre_nfft, self.pre_gi = pre_nfft, pre_gi
        self.pre_dx, self.pre_cred = pre_dx, pre_cred
        self.l1b_cells, self.l1d_cells = l1b_cells, l1d_cells
        self.freq_interleaver = bool(freq_interleaver)
        self.plps = plps
        self.sbs_null_signalled = sbs_null_signalled
        self.label = label
        # --- derived, from spec_pilots (M5's gated tables) -----------------
        self.n_normal = P.data_cells(nfft, cred, pattern, 1, False)
        tot, act, nnull = P.sbs_cells(nfft, cred, pattern, 1)
        self.n_sbs_total, self.n_sbs_active, self.n_null = tot, act, nnull
        self.null_low, self.null_high = P.sbs_null_split(nnull)

    # -- symbol classification ------------------------------------------
    def is_sbs(self, l):
        return (self.sbs_first and l == 0) or (self.sbs_last and
                                               l == self.nsym - 1)

    def n_sbs(self):
        return int(self.sbs_first) + int(self.sbs_last)

    # -- the closed identity that gates the whole thing ------------------
    def pool_size(self, preamble_spare):
        """Available data cells in the Subframe, INCLUDING the Preamble spare.

        This is the number that must equal sum(plp_size) + dummy cells, and on
        a single-PLP Subframe it must equal L1D_plp_size EXACTLY.
        """
        nsbs = self.n_sbs()
        return (preamble_spare
                + (self.nsym - nsbs) * self.n_normal
                + nsbs * self.n_sbs_active)

    def preamble_cells(self):
        """(cells in Preamble symbol 0, cells in each later Preamble symbol)."""
        n0 = S.noc(self.pre_nfft, 4)                      # 7.2.5.1 minimum NoC
        first = _preamble_data_count(self.pre_nfft, self.pre_gi, self.pre_dx, 4)
        rest = _preamble_data_count(self.pre_nfft, self.pre_gi, self.pre_dx,
                                    self.pre_cred)
        return first, rest, n0

    def preamble_spare(self):
        first, rest, _ = self.preamble_cells()
        total = first + rest * (self.np_sym - 1)
        return total - self.l1b_cells - self.l1d_cells

    # -- constructors -----------------------------------------------------
    @classmethod
    def from_l1(cls, l1b, l1d_fields, pre_nfft, pre_gi, pre_dx, l1b_mode,
                label=""):
        """Build from a decoded L1-Basic dict + L1-Detail field list.

        `l1d_fields` is the [(path, name, value)] / [{"name","value"}] shape
        that m8_l1.py emits, so this consumes the JSON directly.
        """
        f = {}
        plps, cur, sub = [], None, {}
        for e in l1d_fields:
            n = e["name"] if isinstance(e, dict) else e[1]
            v = e["value"] if isinstance(e, dict) else e[2]
            if n == "L1D_plp_id":
                if cur:
                    plps.append(cur)
                cur = {"id": v}
            elif cur is not None:
                cur[n] = v
            else:
                sub[n] = v
            f[n] = v
        if cur:
            plps.append(cur)
        core = {}
        out = []
        for p in plps:
            mod = _MOD_NAME[p["L1D_plp_mod"]]
            if p.get("L1D_plp_layer", 0) == 0 and "L1D_plp_CTI_depth" in p:
                core = p
            out.append(dict(
                id=p["id"], layer=p.get("L1D_plp_layer", 0),
                lls=bool(p.get("L1D_plp_lls_flag", 0)),
                start=p["L1D_plp_start"], size=p["L1D_plp_size"],
                mod=mod, rate=_RATE_NAME[p["L1D_plp_cod"]],
                ninner=64800 if p["L1D_plp_fec_type"] else 16200,
                cells_per_fec=(64800 if p["L1D_plp_fec_type"] else 16200)
                // _MOD_BITS[mod],
                ti_mode=p.get("L1D_plp_TI_mode"),
                cti_depth=p.get("L1D_plp_CTI_depth",
                                core.get("L1D_plp_CTI_depth")),
                cti_start_row=p.get("L1D_plp_CTI_start_row",
                                    core.get("L1D_plp_CTI_start_row")),
                cti_fec_block_start=p.get("L1D_plp_CTI_fec_block_start"),
                cti_extended=p.get("L1D_plp_TI_extended_interleaving",
                                   core.get("L1D_plp_TI_extended_interleaving",
                                            0)),
                ldm_injection=p.get("L1D_plp_ldm_injection_level")))
        return cls(
            nfft=_FFT_SIZE[l1b["L1B_first_sub_fft_size"]],
            gi=_LS.GUARD_INTERVAL_SAMPLES[l1b["L1B_first_sub_guard_interval"]],
            cred=l1b["L1B_first_sub_reduced_carriers"],
            pattern=P.PATTERNS[l1b["L1B_first_sub_scattered_pilot_pattern"]],
            nsym=l1b["L1B_first_sub_num_ofdm_symbols"] + 1,
            sbs_first=l1b["L1B_first_sub_sbs_first"],
            sbs_last=l1b["L1B_first_sub_sbs_last"],
            np_sym=l1b["L1B_preamble_num_symbols"] + 1,
            pre_nfft=pre_nfft, pre_gi=pre_gi, pre_dx=pre_dx,
            pre_cred=l1b["L1B_preamble_reduced_carriers"],
            l1b_cells=S.L1_BASIC_CELLS_PRINTED[l1b_mode],
            l1d_cells=l1b["L1B_L1_Detail_total_cells"],
            freq_interleaver=sub.get("L1D_frequency_interleaver", 1),
            sbs_null_signalled=sub.get("L1D_sbs_null_cells"),
            plps=out, label=label)

    @classmethod
    def rf33_legacy(cls):
        """The module constants above, expressed as a Geometry."""
        return cls(nfft=NFFT, gi=GI, cred=CRED, pattern=PATTERN, nsym=NSYM,
                   sbs_first=True, sbs_last=True, np_sym=1,
                   pre_nfft=NFFT, pre_gi=GI, pre_dx=PREAMBLE_DX,
                   pre_cred=PREAMBLE_CRED, l1b_cells=L1B_CELLS,
                   l1d_cells=L1D_CELLS, freq_interleaver=True,
                   plps=[dict(id=0, layer=0, lls=True, **{
                             k: v for k, v in PLP0.items()
                             if k in ("start", "size", "mod", "rate",
                                      "ninner")}),
                         dict(id=16, layer=0, lls=False, **{
                             k: v for k, v in PLP16.items()
                             if k in ("start", "size", "mod", "rate",
                                      "ninner")})],
                   label="RF33 (module constants)")


def _preamble_data_count(nfft, gi, dx, cred):
    _lo, _n, _pilot, _cp, data, _amp = _pre_geo(nfft, gi, dx, cred)
    return len(data)


def _pre_geo(nfft, gi, dx, cred):
    from m3_preamble import preamble_geometry
    return preamble_geometry(nfft, gi, dx, cred)


# ---------------------------------------------------------------------------
# generic demodulation, driven by a Geometry
# ---------------------------------------------------------------------------

def _demod_preamble_g(y, t0, g, s, shift=0):
    """Equalised data cells of Preamble symbol `s`.  A/322 7.2.5.1: symbol 0
    runs at the MINIMUM NoC for its FFT size, the rest at the signalled Cred."""
    cred = 4 if s == 0 else g.pre_cred
    lo, n, pilot, cp, data, _ = _pre_geo(g.pre_nfft, g.pre_gi, g.pre_dx, cred)
    w = t0 + (g.pre_nfft + g.pre_gi) * s + g.pre_gi
    Y = np.fft.fftshift(np.fft.fft(y[w:w + g.pre_nfft]))
    amp = S.PREAMBLE_PILOT_BOOST[(g.pre_nfft, g.pre_gi)][1]
    ref = np.array(S.pilot_values(n, amp))
    hp = Y[fft_bins(g.pre_nfft, lo + pilot + shift)] / ref[pilot]
    H = _interp(np.arange(n), pilot, hp)
    z = Y[fft_bins(g.pre_nfft, lo + data + shift)] / H[data]
    coh = float(abs(np.sum(hp[:-1] * np.conj(hp[1:])))
                / max(np.sum(np.abs(hp) ** 2), 1e-30))
    return z, coh


def _demod_data_g(y, t0, g, l, shift=0):
    """Equalised data cells of DATA symbol `l`, in CARRIER order."""
    sbs = g.is_sbs(l)
    noc = P.NOC[(g.nfft, g.cred)]
    lo, _ = S.carrier_abs_range(g.nfft, g.cred)
    w = (t0 + (g.pre_nfft + g.pre_gi) * g.np_sym
         + (g.nfft + g.gi) * l + g.gi)
    Y = np.fft.fftshift(np.fft.fft(y[w:w + g.nfft]))
    ref = np.array(S.pilot_values(noc, 1.0))
    dx, dy = P.dxdy(g.pattern)
    pk = (np.arange(0, noc, dx) if sbs
          else np.arange(dx * (l % dy), noc, dx * dy))
    pk = np.unique(np.concatenate([pk, [0, noc - 1]]))
    hp = Y[fft_bins(g.nfft, lo + pk + shift)] / ref[pk]
    H = _interp(np.arange(noc), pk, hp)
    used = np.zeros(noc, bool)
    used[np.array(sorted(P.pilot_carriers(g.nfft, g.cred, g.pattern, l, sbs)),
                  int)] = True
    d = np.flatnonzero(~used)
    z = Y[fft_bins(g.nfft, lo + d + shift)] / H[d]
    coh = float(abs(np.sum(hp[:-1] * np.conj(hp[1:])))
                / max(np.sum(np.abs(hp) ** 2), 1e-30))
    return z, coh


def cell_pool_g(y, t0, g, shift=0, report=None, fi_offset=None):
    """Subframe cell pool in cell order, for ANY geometry.

    `fi_offset` is the A/322 7.3 frequency-interleaver symbol counter's origin
    for the first DATA symbol.  M4's ASSUMPTION D2 settled it as "the Preamble
    counts as symbol 0", i.e. offset == NP; it is exposed so it can be swept
    against the dummy-cell referee rather than assumed on a new multiplex.
    """
    rep = report if report is not None else {}
    fi0 = g.np_sym if fi_offset is None else fi_offset

    parts, owner, cohs = [], [], []
    pre = []
    for s in range(g.np_sym):
        z, coh = _demod_preamble_g(y, t0, g, s, shift)
        pre.append(FI.deinterleave(z, g.pre_nfft, s, direction="forward",
                                   toggle="i"))
        cohs.append(coh)
    # np_sym == 0 is legitimate, not an edge case: a Subframe after the first
    # carries no Preamble, so it contributes no spare and no L1 cells.
    pre = np.concatenate(pre) if pre else np.zeros(0, complex)
    spare = pre[g.l1b_cells + g.l1d_cells:]
    parts.append(spare)
    owner.append(np.full(len(spare), -1, int))

    dcoh = []
    for l in range(g.nsym):
        z, coh = _demod_data_g(y, t0, g, l, shift)
        dcoh.append(coh)
        x = (FI.deinterleave(z, g.nfft, l + fi0, direction="forward",
                             toggle="i") if g.freq_interleaver
             else np.asarray(z))
        if g.is_sbs(l):
            x = x[g.null_low:len(x) - g.null_high]
        parts.append(x)
        owner.append(np.full(len(x), l, int))

    pool = np.concatenate(parts)
    rep.update(n_preamble_spare=int(len(spare)),
               n_preamble_total=int(len(pre)),
               n_normal=int(g.n_normal), n_sbs_active=int(g.n_sbs_active),
               n_null=int(g.n_null), null_low=int(g.null_low),
               pool=int(len(pool)),
               pool_predicted=int(g.pool_size(g.preamble_spare())),
               preamble_coherence=cohs,
               data_coherence_mean=float(np.mean(dcoh)),
               freq_interleaver=g.freq_interleaver, fi_offset=int(fi0))
    return pool, dict(symbol_of=np.concatenate(owner), report=rep)


def selftest_geometry(verbose=True):
    """The refactor gate: closed arithmetic, plus the RF33 path reproduced.

    A refactor that changes a number is a bug.  This checks the arithmetic
    identities that do not need a capture, and `--rf33` in m10_core.py checks
    the cell-for-cell identity against `cell_pool()` on real air.
    """
    ok = True
    g = Geometry.rf33_legacy()
    checks = [
        ("preamble spare == 4851 - 484 - 880", g.preamble_spare() == 3487),
        ("normal data cells == 5999", g.n_normal == 5999),
        ("SBS total/active/null == 5136/5009/127",
         (g.n_sbs_total, g.n_sbs_active, g.n_null) == (5136, 5009, 127)),
        ("null split == 63/64", (g.null_low, g.null_high) == (63, 64)),
        ("pool == 211472", g.pool_size(g.preamble_spare()) == 211472),
        ("PLP0+PLP16 == 207900 leaves 3572 dummy",
         g.pool_size(g.preamble_spare()) - 199800 - 8100 == 3572),
    ]
    for name, c in checks:
        ok &= c
        if verbose:
            print("    %s  RF33 via Geometry: %s" % ("PASS" if c else "FAIL",
                                                     name))
    return ok
