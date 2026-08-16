#!/usr/bin/env python3
"""M44 (E82) -- the LDM / CTI core layer as a LIVE stage, not an offline tool.

WHAT WAS MISSING
----------------
E49/E69/E76 decoded the Layered-Division-Multiplexing CORE layer of RF25 (Fox
45.100) and RF30, perfectly, from banked captures -- `m10_core.py`,
`e49_core.py`, `e69_stream.py`.  Every one of those tools is a BATCH tool:
load the whole capture, demodulate every Frame, concatenate the cell stream,
de-interleave it once, decode every FEC Block.  The live chain
(`m11_stream.FrontEnd` -> `m9_fast.FrameDecoder`) is the opposite shape: one
Frame in, one Frame's Baseband Packets out, forever, and it is hard-wired to
RF33's single-layer geometry and its INTRA-frame Hybrid time interleaver.

The gap is not the demodulator.  It is the interleaver:

    HTI (RF33)   one Frame's FEC Blocks live entirely inside that Frame.
                 Frames are independent -> trivially streamable, trivially
                 parallel, and a lost Frame costs exactly that Frame.
    CTI (RF25)   a convolutional interleaver whose delay lines are up to
                 Nrows*(Nrows-1) = 1,047,552 cells deep -- 0.92 of a Frame --
                 and whose commutator NEVER RESETS.  A FEC Block is a diagonal
                 through ~two Frames, and the whole stream is one continuous
                 index space anchored on a PHASE that is signalled once per
                 Frame in L1 and can be lost.

So this module is three things: a geometry-driven Frame demodulator (any
FFT/GI/pilot pattern/PLP layout, from L1), a STREAMING CTI de-interleaver, and
-- the part the campaign has actually been burned by -- a PHASE TRACKER.

WHY THE PHASE TRACKER IS THE WHOLE JOB
--------------------------------------
E76's near-miss, quoted from the lab log because it is the requirement:

    e49_stream --start-frame   FEC (10 blocks)
        245 / 246 / 247        0/10
        248                   10/10   <- BCH zero on all ten

One Frame -- 242 ms -- between "RF25 is dead" and "RF25 decodes perfectly",
with every quality metric (pilot coherence 0.985-0.994, dummy-cell SNR
12.85 dB, constellation p99 tighter than the control) reading EXCELLENT on
both sides of the line.  E49 wrote the law: *interleaver clock offset
impersonates weak SNR*.

A human found 248 by sweeping.  A tuner cannot sweep, so the phase has to be
ACQUIRED, and the acquisition has to be re-runnable at any moment because a
live front end re-acquires its bootstrap whenever the samples break.

Three sources, in order, and the FEC syndrome is the arbiter of all three:

  1. **L1, on the Frame we are actually collecting.**  This is the structural
     fix for E76: the batch tools took `L1D_plp_CTI_start_row` from whatever
     Frame m8_l1 happened to verify on (t = 60.0 s) and applied it to a cell
     stream that started at t = 0.  Here the L1 is decoded FROM THE SAME
     WINDOW whose cells open the stream, so the two cannot disagree.  There is
     no `--start-frame` and there is nothing to hand-set.
  2. **Dead reckoning.**  A/322 9.3.9.1: the commutator advances by plp_size
     rows per Subframe, so start_row(n) = (start_row(n0) + (n-n0)*plp_size)
     mod Nrows.  The front end counts Frames, so once ONE (frame, start_row)
     pair is known every later Frame's phase is known WITHOUT decoding L1
     again -- which matters, because L1 verifies on only ~1 Frame in 9 on this
     carrier even at 13 dB (E76's own open question).
  3. **A bounded sweep of the same axis, arbitrated by the syndrome.**  If (1)
     and (2) are both unavailable (a bootstrap re-acquisition throws the Frame
     count away), sweep k Frames of commutator advance around the dead-
     reckoned estimate and let LDPC convergence pick.  This is E76's sweep,
     automated -- and it is the FALLBACK, not the mechanism.

A phase is never *believed*: a candidate is PROVISIONAL until `LOCK_BLOCKS`
FEC Blocks have been decoded at it and at least `LOCK_FRAC` of them converged
with a zero BCH syndrome.  A wrong phase produces ~chance unsatisfied-check
counts and converges nothing, which is what makes the syndrome a competent
referee here and not merely a quality metric.

WHY THE CORE LAYER NEEDS NO CANCELLATION (m10_core's argument, unchanged)
------------------------------------------------------------------------
Both RF25 and RF30 signal `L1D_plp_lls_flag = 1` on the CORE layer, so the
service list rides the layer that is decodable with the Enhanced layer treated
as interference -- which is what LDM is for.  For a QPSK core the max-log LLR
depends only on the SIGNS of I and Q, so the demapper is scale-invariant and
does not even need the injection level.

NOTHING HERE IS CHANNEL-SPECIFIC.  There is no "if rf == 25".  The layer, the
modulation, the code rate, Ninner, the CTI depth and the commutator phase all
come from the decoded L1-Detail, so RF30 (also LDM, also CTI, E69) runs the
same code with different numbers.
"""
from __future__ import annotations

import collections
import os
import sys
import threading
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m3_freqint as FI                                           # noqa: E402
import m3_ldpc as LD                                              # noqa: E402
import m3_spec as S                                               # noqa: E402
import m4_scrambler as SC                                         # noqa: E402
import m6_bicm as B                                               # noqa: E402
import m6_cells as C                                              # noqa: E402
import m10_cti as CTI                                             # noqa: E402
import m16_margin as M16                                          # noqa: E402
import spec_bitint as BI                                          # noqa: E402
import spec_pilots as P                                           # noqa: E402
from m3_preamble import fft_bins                                  # noqa: E402

FS_POST = 6.912e6
BOOTSTRAP = 13824

# --- E86: the exact-BP rescue, ported from E58/E60's m9_fast work ----------
# Env-tunable so the frozen-capture A/B can flip exactly one variable, and
# so a live run can shed the cost under a deadline (E52's shed law: a fade
# must never stall the pipeline).  SP_BUDGET_S=None means "no deadline",
# which is right offline and wrong live.
# DEFAULT OFF on this path, and that is a MEASURED choice, not a doubt.
# On rf25_live_1210 (4.00 dB) the rescue alone lifted 4.2% -> 13.0%, but with
# CE smoothing on it added EXACTLY NOTHING (13.1% either way) while costing
# ~8% of the run.  The two levers were not additive because they were curing
# the same illness from opposite ends: the near-misses SP was rescuing exist
# because the channel estimate was noisy, and smoothing the estimate makes
# min-sum close them by itself.  Re-test before assuming this holds at other
# operating points -- a recorded negative is a reference, not an authority.
SP_RESCUE = os.environ.get("ATSC3_LDM_SP", "0") not in ("0", "off", "no")
SP_NEAR_MAX = int(os.environ.get("ATSC3_LDM_SP_NEAR", "1500"))
SP_ITERS = int(os.environ.get("ATSC3_LDM_SP_ITERS", "100"))
SP_MAX_BLOCKS = int(os.environ.get("ATSC3_LDM_SP_BLOCKS", "8"))
_sp_budget = os.environ.get("ATSC3_LDM_SP_BUDGET", "0.12")
SP_BUDGET_S = float(_sp_budget) if _sp_budget not in ("", "none") else None

# --- E86c: CE smoothing, ported from E58/E60's m9_fast work ----------------
# The biggest of the three levers (+1.05 dB in E58) and the one that was NOT
# a function call: m9_fast fuses it into its own channel estimator, so this
# is a reimplementation against m44's grouped filterbank, not a wrapper.
# Same knobs, same measured defaults, so the two paths can be compared.
CE_W = int(os.environ.get("ATSC3_LDM_CE_W", "12"))          # 0 disables
CE_DETECT = float(os.environ.get("ATSC3_LDM_CE_DETECT", "6.0"))
CE_ABS = float(os.environ.get("ATSC3_LDM_CE_ABS", "10.0"))
CE_NORM = os.environ.get("ATSC3_LDM_CE_NORM", "1") not in ("0", "off", "no")


# ===========================================================================
# 1 -- THE PLAN: everything this module does, derived from a decoded L1
# ===========================================================================

class LdmPlan:
    """The whole configuration, read off L1-Detail.  No channel constants.

    `is_ldm` is true when L1 signals a PLP on layer 1; `uses_cti` when the
    chosen PLP signals TI mode 1 (the Convolutional Time Interleaver).  The
    live chain asks those two questions and nothing else -- it never asks
    which RF channel it is tuned to.
    """

    def __init__(self, g, core, label=""):
        self.g = g
        self.core = core
        self.label = label
        # A multiplex that does NOT use the CTI has no Nrows to read, and
        # asking for one raises.  This class must be CONSTRUCTIBLE for any
        # multiplex, because `uses_cti` is the question the live chain asks
        # it -- a plan that throws on RF33 turns "identify the multiplex"
        # into "crash on the user's television", which is exactly what the
        # first build did.
        self.uses_cti = core.get("ti_mode") == 1 and \
            core.get("cti_depth") is not None
        self.nrows = (CTI.nrows_of(core["cti_depth"],
                                   bool(core["cti_extended"]))
                      if self.uses_cti else None)
        self.ncell = core["cells_per_fec"]
        self.plp_size = core["size"]
        self.plp_start = core["start"]
        self.frame_samples = int(BOOTSTRAP
                                 + (g.pre_nfft + g.pre_gi) * g.np_sym
                                 + (g.nfft + g.gi) * g.nsym)
        self.frame_sec = self.frame_samples / FS_POST
        # A decode window must cover [t0, t0 + frame_samples) plus the fine
        # timing search on each side, plus slack -- the same shape as
        # m11_stream.FRAME_WINDOW, computed instead of written down.
        self.ft_span = 20
        self.frame_window = self.frame_samples + 2 * self.ft_span + 4096
        self.pool_pred = int(g.pool_size(g.preamble_spare()))
        self.core_total = int(sum(p["size"] for p in g.plps
                                  if p["layer"] == 0))
        self.dummy_n = self.pool_pred - self.core_total
        self.is_ldm = any(p["layer"] == 1 for p in g.plps)
        self.mod_bits = BI.MOD_BITS[core["mod"]]
        # the CTI's own reach: an output cell may be fed by an input cell up
        # to Nrows*(Nrows-1) cells LATER in the received stream
        self.cti_reach = (self.nrows * (self.nrows - 1)
                          if self.nrows else None)

    # -- constructors ------------------------------------------------------
    @classmethod
    def _build(cls, l1b, fields, nfft, gi, dx, mode, label, plp_id):
        g = C.Geometry.from_l1(l1b, fields, nfft, gi, dx, mode, label=label)
        self = cls(g, cls.pick_core(g, plp_id), label=label)
        # E85: the plan has to cross a process boundary to reach a demod
        # worker, and a Geometry is not picklable in any way I would trust to
        # stay true.  Keep the DECODED L1 that produced it instead -- plain
        # ints and strings -- and rebuild the plan on the far side from the
        # same constructor.  A worker that reconstructs the geometry from the
        # same signalling cannot silently hold a different one.
        self.spec = dict(l1b={k: int(v) for k, v in l1b.items()},
                         fields=[[p, n, int(v)] for p, n, v in
                                 ([f if not isinstance(f, dict)
                                   else (f["path"], f["name"], f["value"])
                                   for f in fields])],
                         nfft=int(nfft), gi=int(gi), dx=int(dx),
                         mode=int(mode), label=label, plp_id=plp_id)
        return self

    @classmethod
    def from_spec(cls, spec):
        """Rebuild a plan from `to_spec()`, in another process."""
        return cls._build(spec["l1b"], [tuple(f) for f in spec["fields"]],
                          spec["nfft"], spec["gi"], spec["dx"], spec["mode"],
                          spec["label"], spec["plp_id"])

    def to_spec(self):
        return self.spec

    @classmethod
    def from_l1_json(cls, js, label="", plp_id=None):
        return cls._build(js["l1_basic"], js["l1_detail"],
                          js["preamble"]["nfft"], js["preamble"]["gi"],
                          js["preamble"]["dx"], js["preamble"]["l1b_mode"],
                          label, plp_id)

    @classmethod
    def from_l1_result(cls, r, label="", plp_id=None):
        """From an in-memory m8_l1-shaped result (see `l1_from_window`)."""
        return cls._build(r["l1b"], r["l1d"]["fields"], r["ctx"]["nfft"],
                          r["ctx"]["gi"], r["ctx"]["dx"],
                          r["ctx"]["l1b_mode"], label, plp_id)

    @staticmethod
    def pick_core(g, plp_id=None):
        """The layer the service list rides, chosen by what L1 SIGNALS.

        A/322: `L1D_plp_layer` 0 is the Core Layer.  `L1D_plp_lls_flag` marks
        the PLP carrying Low Level Signalling (the SLT), which is the one a
        receiver must decode to find any service at all.  Ties break on size,
        so the data PLP wins over a small signalling PLP such as RF25's
        PLP 16.  An explicit `plp_id` overrides, for instruments.
        """
        by_id = {p["id"]: p for p in g.plps}
        if plp_id is not None:
            return by_id[plp_id]
        cand = [p for p in g.plps if p["layer"] == 0 and p["lls"]]
        if not cand:
            cand = [p for p in g.plps if p["layer"] == 0]
        if not cand:
            raise ValueError("L1 signals no Core-Layer PLP")
        return max(cand, key=lambda p: p["size"])

    def describe(self):
        g = self.g
        return (
            f"{self.label or 'multiplex'}: FFT {g.nfft // 1024}K GI {g.gi} "
            f"Cred {g.cred} {g.pattern} {g.nsym} symbols NP {g.np_sym} -> "
            f"Frame {self.frame_samples} samples "
            f"({self.frame_sec * 1000:.3f} ms)\n"
            f"    layers: "
            + ", ".join(f"PLP{p['id']}(L{p['layer']} {p['mod']} {p['rate']}"
                        + (f" inj{p['ldm_injection']}"
                           if p["ldm_injection"] is not None else "") + ")"
                        for p in g.plps)
            + f"\n    CORE = PLP{self.core['id']}: {self.core['mod']} "
              f"{self.core['rate']} Ninner {self.core['ninner']}, "
              f"{self.plp_size} cells/Frame, {self.ncell} cells/FEC Block "
              f"({self.plp_size / self.ncell:.2f} Blocks/Frame)\n"
            f"    time interleaver: "
            + (f"CTI Nrows {self.nrows} (reach {self.cti_reach} cells = "
               f"{self.cti_reach / self.plp_size:.2f} Frames)"
               if self.uses_cti else f"TI mode {self.core.get('ti_mode')}")
            + f"\n    pool {self.pool_pred} = core {self.core_total} + dummy "
              f"{self.dummy_n}   LDM={self.is_ldm}")


# ===========================================================================
# 2 -- L1 FROM A DECODE WINDOW (not from a file)
# ===========================================================================

def preamble_coherence(y, t0, plan_or_geom):
    """Pilot coherence of Preamble symbol 0 at FFT window t0 + pre_gi.

    This is the fine-timing metric for a general multiplex.  m9_accel's
    `demod_coh` is RF33's (module constants, DATA symbol 1); the Preamble is
    the right symbol here because it is the ONE symbol whose pilot geometry is
    fixed by A/322 7.2.5.1 rather than signalled -- so it works before L1 has
    ever been decoded, which is exactly when acquisition needs it.
    """
    g = getattr(plan_or_geom, "g", plan_or_geom)
    nfft, gi, dx = g.pre_nfft, g.pre_gi, g.pre_dx
    w = t0 + gi
    if w < 0 or w + nfft > len(y):
        return -1.0
    from m3_preamble import pilot_coherence
    Y = np.fft.fftshift(np.fft.fft(y[w:w + nfft]))
    return float(pilot_coherence(Y, nfft, gi, dx, 4, 0))


class PreambleCoh:
    """Batched preamble pilot coherence over a window of candidate t0s.

    The scalar form costs one 8K FFT per candidate; 41 of them is 116 ms of a
    242 ms Frame, which would make fine timing alone the wall.  The candidates
    are 41 unit-strided windows of the same buffer, so one strided (41, NFFT)
    transform serves them all -- m9_accel.demod_coh_batch's trick, for the
    Preamble instead of a data symbol.  Tables are frame-constant and cached.
    """

    def __init__(self, geom):
        g = getattr(geom, "g", geom)
        self.nfft, self.gi, self.dx = g.pre_nfft, g.pre_gi, g.pre_dx
        from m3_preamble import preamble_geometry
        lo, n, pilot, _cp, _d, _o = preamble_geometry(self.nfft, self.gi,
                                                      self.dx, 4)
        amp = S.PREAMBLE_PILOT_BOOST[(self.nfft, self.gi)][1]
        ref = np.asarray(S.pilot_values(n, amp))
        b = fft_bins(self.nfft, lo + pilot)
        ok = (b >= 0) & (b < self.nfft)
        # fftshift is a pure reindex; fold it into the bin list instead
        self.bins = ((b[ok] + self.nfft // 2) % self.nfft).astype(np.intp)
        self.ref = ref[pilot[ok]]

    def scan(self, y, centre, span, workers=1):
        """-> (t0, coherence) with the project's strictly-greater tie-break."""
        cands = np.arange(centre - span, centre + span + 1)
        lo = int(cands[0]) + self.gi
        hi = int(cands[-1]) + self.gi + self.nfft
        if lo < 0 or hi > len(y):
            best = None
            for t in cands:
                c = preamble_coherence(y, int(t), self)
                if best is None or c > best[1]:
                    best = (int(t), c)
            return best
        seg = np.ascontiguousarray(y[lo:hi])
        st = seg.strides[0]
        W = np.lib.stride_tricks.as_strided(
            seg, shape=(len(cands), self.nfft), strides=(st, st))
        # MATERIALISE THE OVERLAP BEFORE THE TRANSFORM.  These 41 rows are
        # unit-strided views of ONE buffer, so every row overlaps the next by
        # nfft-1 samples.  A multi-threaded pocketfft call over that view
        # segfaulted the chain, reproducibly, two Frames after the CTI locked
        # -- a shared read that each worker also wants to stage.  The copy is
        # 5.4 MB and ~1 ms; m9_accel's RF33 twin sidesteps the same hazard by
        # using single-threaded np.fft.  (Law, again: an unbounded/aliased
        # buffer handed to a thread pool is not an optimisation.)
        Wc = np.ascontiguousarray(W)
        try:
            from scipy import fft as _sfft
            Y = _sfft.fft(Wc, axis=1, workers=max(1, workers))
        except Exception:                                          # noqa: BLE001
            Y = np.fft.fft(Wc, axis=1)
        q = Y[:, self.bins] / self.ref
        num = np.abs(np.sum(q[:, :-1] * np.conj(q[:, 1:]), axis=1))
        den = np.maximum(np.sum(np.abs(q) ** 2, axis=1), 1e-30)
        coh = num / den
        best = None
        for t, c in zip(cands, coh):
            if best is None or c > best[1]:
                best = (int(t), float(c))
        return best


class DataCoh:
    """Fine timing off a DATA symbol's scattered pilots, batched.

    WHY NOT THE PREAMBLE.  The Preamble metric is what identification needs --
    its geometry is fixed by A/322 7.2.5.1 and therefore exists before L1 --
    but it is a poor RULER.  Its pilots sit every DX = 6 carriers, so a timing
    error tau rotates adjacent pilots by only 2*pi*tau*6/N, and on RF25 (GI
    512) the coherence plateau measured 200 SAMPLES WIDE and flat to four
    decimals across it.  Inside the tracker's +-20 window that is not a peak
    at all: argmax picks noise, t0 random-walks, and when it walks past the
    right-hand edge of the plateau the FFT window reaches into the next symbol
    and the FEC collapses.  Measured, on this capture: 100% BCH-clean for 50
    Frames, then 87% by Frame 75, with every other metric still perfect.
    (The same costume E49 and E76 both recorded -- excellent quality
    indicators, correctness gone.)

    A DATA symbol's scattered pilots for SP12_4 are DX*DY = 48 carriers apart
    -- eight times the phase slope per sample of timing error -- so the same
    statistic on the same air has a real peak.  This is the metric m9_accel
    uses for RF33; here it is built from the Geometry instead of from module
    constants.
    """

    def __init__(self, plan, l=None):
        g = plan.g
        self.g = g
        # the first NON-SBS data symbol: an SBS symbol has a denser pilot
        # pattern and a different cell count, and one class is all we need
        self.l = l if l is not None else next(
            i for i in range(g.nsym) if not g.is_sbs(i))
        sg = _sym_geometry(g, self.l, False)
        self.bins = ((sg["bins_pk"] + g.nfft // 2) % g.nfft).astype(np.intp)
        self.ref = sg["refpk"]
        self.nfft = g.nfft
        self.off = ((g.pre_nfft + g.pre_gi) * g.np_sym
                    + (g.nfft + g.gi) * self.l + g.gi)

    def scan(self, y, centre, span, workers=1):
        cands = np.arange(centre - span, centre + span + 1)
        lo = int(cands[0]) + self.off
        hi = int(cands[-1]) + self.off + self.nfft
        if lo < 0 or hi > len(y):
            return None
        seg = np.ascontiguousarray(y[lo:hi])
        st = seg.strides[0]
        W = np.lib.stride_tricks.as_strided(
            seg, shape=(len(cands), self.nfft), strides=(st, st))
        Wc = np.ascontiguousarray(W)        # see PreambleCoh.scan
        try:
            from scipy import fft as _sfft
            Y = _sfft.fft(Wc, axis=1, workers=max(1, workers))
        except Exception:                                          # noqa: BLE001
            Y = np.fft.fft(Wc, axis=1)
        q = Y[:, self.bins] / self.ref
        num = np.abs(np.sum(q[:, :-1] * np.conj(q[:, 1:]), axis=1))
        den = np.maximum(np.sum(np.abs(q) ** 2, axis=1), 1e-30)
        coh = num / den
        best = None
        for t, c in zip(cands, coh):
            if best is None or c > best[1]:
                best = (int(t), float(c))
        return best


def fine_timing_g(y, plan, centre, span=None, cache={}):
    """(t0, coherence) -- the general-geometry twin of m6_cells.fine_timing."""
    span = plan.ft_span if span is None else span
    g = getattr(plan, "g", plan)
    key = (g.pre_nfft, g.pre_gi, g.pre_dx)
    pc = cache.get(key)
    if pc is None:
        pc = cache[key] = PreambleCoh(g)
    return pc.scan(y, centre, span)


def l1_from_window(y, t0, ps=None, quiet=True):
    """Decode L1-Basic + L1-Detail from ONE decode window, in memory.

    `m8_l1.solve()` takes a PATH and re-reads and re-bootstraps the file; a
    live chain already holds the samples and already knows where the Frame
    starts.  Everything below is m8_l1's own code -- `demod_preamble_symbol`,
    `l1basic_of`, `geometry`, `decode`, `parse_l1d` -- fed from a window.

    Returns None if L1-Basic does not verify, or a dict with `ok` False if
    L1-Detail does not; the caller must treat both as "no phase from L1 this
    Frame" and fall back, never as an error.
    """
    import m3_l1basic as L1B                                       # noqa: F401
    import m4_l1detail as D4
    import m8_l1 as L8
    import m3_crc as CRC
    from m3_preamble import channel_estimate, preamble_geometry

    # A/322 Table H.1.1 is keyed by the bootstrap's preamble_structure, which
    # the front end has already consumed.  The Preamble's own pilot geometry
    # is fixed at Cred 4 / the minimum NoC for its FFT size (7.2.5.1), and the
    # bootstrap fixes FFT and GI, so scan the small menu the bootstrap allows
    # and let pilot coherence choose -- the same statistic m3_preamble.analyse
    # uses, on a window instead of a file.
    from m3_preamble import pilot_coherence
    menu = (S.PREAMBLE_STRUCTURE if ps is None
            else {ps: S.PREAMBLE_STRUCTURE[ps]})
    best = None
    for _ps, Pt in menu.items():
        nfft, gi, dx, mode = Pt["fft"], Pt["gi"], Pt["dx"], Pt["l1b_mode"]
        w = t0 + gi
        if w < 0 or w + nfft > len(y):
            continue
        try:
            Y = np.fft.fftshift(np.fft.fft(y[w:w + nfft]))
            c = pilot_coherence(Y, nfft, gi, dx, 4, 0)
        except KeyError:
            # a Table H.1.1 row this build has no pilot-boost entry for; it
            # cannot be the structure we are locked to, because the bootstrap
            # that got us here already named one we can demodulate
            continue
        if best is None or c > best[0]:
            best = (c, nfft, gi, dx, mode, Y)
    if best is None:
        return None
    coh, nfft, gi, dx, mode, Y = best
    cred, shift = 4, 0
    lo, n, pilot, cp, data, _ = preamble_geometry(nfft, gi, dx, cred)
    H, hp = channel_estimate(Y, nfft, gi, dx, cred, shift)
    bins = fft_bins(nfft, lo + data + shift)
    ok = (bins >= 0) & (bins < nfft)
    z0 = Y[bins[ok]] / H[data[ok]]
    x0 = FI.deinterleave(z0, nfft, 0, direction="forward", toggle="i")
    f = D4.l1basic_of(x0, mode)
    if f is None:
        return None
    npre = f["L1B_preamble_num_symbols"] + 1
    nb = S.L1_BASIC_CELLS_PRINTED[mode]
    tot = f["L1B_L1_Detail_total_cells"]
    parts, cohs = [x0[nb:]], []
    credp = f["L1B_preamble_reduced_carriers"]
    for s in range(1, npre):
        zs, cs, _nd = L8.demod_preamble_symbol(
            y, t0 + gi + (nfft + gi) * s, nfft, gi, dx, credp, shift)
        parts.append(FI.deinterleave(zs, nfft, s, direction="forward",
                                     toggle="i"))
        cohs.append(cs)
    stream = np.concatenate(parts)[:tot]
    if len(stream) < tot:
        return None
    g = L8.geometry(f["L1B_L1_Detail_size_bytes"], f["L1B_L1_Detail_fec_type"])
    if g["cells"] != tot:
        # the free gate M8 earned: a computed cell count that disagrees with
        # the signalled one means the geometry reading is wrong, not the air
        return dict(ok=False, why="cells %d != signalled %d"
                    % (g["cells"], tot))
    xd = L8.l1d_deinterleave(stream, tot, npre)
    r = L8.decode(xd, g)
    ctx = dict(nfft=nfft, gi=gi, dx=dx, l1b_mode=mode, NP=npre,
               preamble_coherence=float(coh), extra_coherence=cohs)
    if not (r["converged"] and r["bch_ok"]):
        return dict(ok=False, why="L1-Detail FEC", fec=r, ctx=ctx, l1b=f)
    d = SC.descramble(np.array(r["nouter"][:g["Ksig"]], np.uint8))
    crc_ok = bool(CRC.crc32(d[:g["Ksig"] - 32], init=0xFFFFFFFF)
                  == int("".join(str(int(c)) for c in d[g["Ksig"] - 32:]), 2))
    if not crc_ok:
        return dict(ok=False, why="L1D_crc", fec=r, ctx=ctx, l1b=f)
    parsed = None
    for mm in (1, 0):
        try:
            p = D4.parse_l1d(d, f, mm)
        except Exception:                                          # noqa: BLE001
            continue
        if p["reserved_all_ones"]:
            parsed = p
            break
    if parsed is None:
        parsed = D4.parse_l1d(d, f, 1)
    return dict(ok=True, l1b=f, ctx=ctx, fec=r, crc_ok=True,
                l1d=dict(fields=parsed["fields"]))


def coh_raw(y, t0, ps=None):
    """Preamble pilot coherence with NO geometry assumed beyond the bootstrap.

    A/322 Table H.1.1 keys FFT / GI / pilot DX / L1-Basic Mode off the
    bootstrap's own `preamble_structure`, and 7.2.5.1 fixes the first Preamble
    symbol at the minimum NoC -- so this metric exists BEFORE L1 does, which
    is exactly when identification needs it.
    """
    from m3_preamble import pilot_coherence
    menu = (S.PREAMBLE_STRUCTURE if ps is None
            else {ps: S.PREAMBLE_STRUCTURE[ps]})
    best = -1.0
    for _ps, Pt in menu.items():
        nfft, gi, dx = Pt["fft"], Pt["gi"], Pt["dx"]
        w = t0 + gi
        if w < 0 or w + nfft > len(y):
            continue
        try:
            Y = np.fft.fftshift(np.fft.fft(y[w:w + nfft]))
            best = max(best, pilot_coherence(Y, nfft, gi, dx, 4, 0))
        except KeyError:
            continue
    return best


def plateau_t0(y, t0_nominal, ps=None, lo=-40, hi=8):
    """Frame start by the Preamble's coherence plateau.

    The plateau is as wide as the guard interval minus the channel's delay
    spread (measured 200 samples on RF25's GI of 512) and being EARLY inside
    it is harmless while being LATE is inter-symbol interference, so the scan
    is deliberately asymmetric around nominal.
    """
    best = None
    for t in range(t0_nominal + lo, t0_nominal + hi):
        c = coh_raw(y, t, ps)
        if best is None or c > best[1]:
            best = (t, c)
    return best


def sniff(raw, rate, log=print, tries=8):
    """Identify a multiplex from raw capture-rate samples, before any decode.

    Returns `(plan, info)`.  `plan is None` means "this is not a multiplex
    this module owns" -- either L1 did not verify, or L1 verified and said the
    PLPs use the Hybrid time interleaver, which is m9_fast's machine.  The
    caller must be able to fall through to the existing chain unchanged, so
    this function reads samples and returns; it changes nothing.
    """
    import m2_pilots as MP
    from atsc3 import bootstrap as bs
    x = np.asarray(raw)
    head_n = min(len(x), int(0.30 * rate))
    try:
        hits = MP.find_bootstraps(MP.resample_to(x[:head_n], rate, bs.FS))
    except Exception as e:                                         # noqa: BLE001
        return None, dict(why=f"bootstrap search failed: {e}")
    if not hits:
        return None, dict(why="no bootstrap in the sniff window")
    h = hits[0]
    fs_post = h["fields"]["post_bootstrap_sample_rate_hz"]
    ps = h["fields"].get("preamble_structure")
    f0 = int(round(h["position"] * rate / bs.FS))
    y = x[f0:]
    if abs(rate - fs_post) > 1.0:
        y = MP.resample_to(y, rate, fs_post)
    y = np.asarray(y, np.complex128)
    y = y * np.exp(-2j * np.pi * h["fine_cfo_hz"]
                   * np.arange(len(y)) / fs_post)
    info = dict(preamble_structure=ps, cfo_hz=float(h["fine_cfo_hz"]),
                peak_ratio=float(h.get("mean_peak_ratio", 0.0)))
    off = 0
    for k in range(tries):
        t0, coh = plateau_t0(y[off:], BOOTSTRAP, ps)
        L = l1_from_window(y[off:], t0, ps=ps)
        if L is not None and L.get("ok"):
            plan = LdmPlan.from_l1_result(L, label="sniffed")
            info.update(l1_frames_tried=k + 1, coherence=round(coh, 4))
            return plan, info
        # step to the next bootstrap; the Frame length is not known yet, so
        # it is FOUND rather than assumed
        rest = y[off + BOOTSTRAP:]
        if len(rest) < int(0.35 * fs_post):
            break
        try:
            nh = MP.find_bootstraps(MP.resample_to(
                rest[:int(3.0 * fs_post)], fs_post, bs.FS))
        except Exception:                                          # noqa: BLE001
            break
        if not nh:
            break
        off += BOOTSTRAP + int(round(nh[0]["position"] * fs_post / bs.FS))
    info["why"] = f"L1 did not verify in {tries} Frames"
    return None, info


def start_row_of(l1res, plp_id):
    for pth, name, v in l1res["l1d"]["fields"]:
        pass
    cur, out = None, {}
    for pth, name, v in l1res["l1d"]["fields"]:
        if name == "L1D_plp_id":
            cur = v
        elif cur == plp_id and name == "L1D_plp_CTI_start_row":
            out["start_row"] = v
        elif cur == plp_id and name == "L1D_plp_CTI_fec_block_start":
            out["fec_block_start"] = v
    return out


# ===========================================================================
# 3 -- THE FRAME DEMODULATOR, geometry-driven and batched
# ===========================================================================

def _sym_geometry(g, l, sbs):
    """Per-DATA-symbol demod tables for an arbitrary Geometry.

    This is `m9_accel._geometry` with the RF33 module constants replaced by a
    Geometry object -- same arithmetic, same spec_pilots tables, so a symbol's
    cells come out in the same order the batch chain (`m6_cells._demod_data_g`)
    produces them.  `gate_e82.py` checks that cell for cell.
    """
    noc = P.NOC[(g.nfft, g.cred)]
    lo, _ = S.carrier_abs_range(g.nfft, g.cred)
    ref = np.asarray(S.pilot_values(noc, 1.0))
    dx, dy = P.dxdy(g.pattern)
    pk = (np.arange(0, noc, dx) if sbs
          else np.arange(dx * (l % dy), noc, dx * dy))
    pk = np.unique(np.concatenate([pk, [0, noc - 1]]))
    used = np.zeros(noc, bool)
    used[np.array(sorted(P.pilot_carriers(g.nfft, g.cred, g.pattern, l, sbs)),
                  int)] = True
    d = np.flatnonzero(~used)
    return dict(pk=pk, refpk=ref[pk], bins_pk=fft_bins(g.nfft, lo + pk),
                d=d, bins_d=fft_bins(g.nfft, lo + d), kk=np.arange(noc))


class FrameDemod:
    """One Frame's decode window -> the CORE PLP's cells, in cell order.

    The structure is `m9_fast.cell_pool_fast`'s -- symbols grouped by pilot
    class, one batched FFT, one gather/interp/scatter per class -- because
    that is what makes it real-time-shaped: RF25 carries 189 data symbols per
    Frame against RF33's 35, so a per-symbol NumPy dispatch loop pays 189
    times the overhead for the same arithmetic.  Everything is a pure function
    of the geometry, so the tables are built once and cached.
    """

    def __init__(self, plan, threads=4, ex=None):
        self.plan = plan
        self.threads = threads
        self.ex = ex
        self._tab = None
        self._dv = None
        self.tm = collections.Counter()
        self._lock = threading.Lock()

    # -- tables ------------------------------------------------------------
    def tables(self):
        if self._tab is not None:
            return self._tab
        with self._lock:
            if self._tab is not None:
                return self._tab
            g = self.plan.g
            dx, dy = P.dxdy(g.pattern)
            classes = {}
            for l in range(g.nsym):
                sbs = g.is_sbs(l)
                key = ("sbs",) if sbs else ("norm", l % dy)
                classes.setdefault(key, []).append(l)
            fi0 = g.np_sym
            tabs = []
            for key, rows in classes.items():
                sbs = key[0] == "sbs"
                sg = _sym_geometry(g, rows[0], sbs)
                pk, kk = sg["pk"], sg["kk"]
                j = np.clip(np.searchsorted(pk, kk, side="right") - 1,
                            0, len(pk) - 2)
                wt = (kk - pk[j]) / (pk[j + 1] - pk[j])
                nd = len(sg["d"])
                if g.freq_interleaver:
                    Hmat = np.stack([FI.interleaving_sequence(
                        g.nfft, l + fi0, nd, direction="forward", toggle="i")
                        for l in rows]).astype(np.int32)
                else:
                    Hmat = np.tile(np.arange(nd, dtype=np.int32),
                                   (len(rows), 1))
                trim = ((g.null_low, nd - g.null_high) if sbs else (0, nd))
                tabs.append(dict(rows=np.array(rows), sbs=sbs, g=sg, j=j,
                                 wt=wt, nd=nd, Hmat=Hmat, trim=trim))
            lens = np.empty(g.nsym, int)
            for t in tabs:
                a, b = t["trim"]
                lens[t["rows"]] = b - a
            self._tab = dict(tabs=tabs, lens=lens,
                             step=g.nfft + g.gi, gi=g.gi,
                             data_off=(g.pre_nfft + g.pre_gi) * g.np_sym)
            return self._tab

    def prewarm(self):
        """Build every table before the clock starts.

        The frequency-interleaver sequences are 8192-step Python LFSR walks,
        one per OFDM symbol: 189 of them is ~1.5 M interpreted steps and
        several seconds.  It is initialisation, not decoding -- on live air it
        runs while the radio settles, exactly as m9_fast.prewarm does.
        """
        self.tables()
        self.scrambler()
        self.dummy_values()
        pl = self.plan
        ch = self.chain()
        z = np.zeros((2, ch.ninner))
        LD.min_sum_decode_batch(z, ch.checks, ch.ninner, iters=1, alpha=1.0,
                                dtype=np.float32)
        SC.sequence(max(ch.kpayload, pl.pool_pred))
        np.fft.fft(np.zeros(pl.g.nfft, complex))

    def chain(self, iters=50):
        ch = getattr(self, "_chain", None)
        if ch is None:
            c = self.plan.core
            ch = self._chain = B.PlpChain(c["ninner"], c["mod"], c["rate"],
                                          iters=iters)
        return ch

    def scrambler(self):
        """A/322 5.2.3's +-1 dummy-cell sequence over the whole pool, CACHED.

        `m4_scrambler.sequence` is a two-level Python loop -- 1.2 M interpreted
        steps for an RF25 pool, ~0.9 s.  Calling it per Frame (the first build
        did, from `dummy_agreement`) was two thirds of the Frame decode and
        invisible in the stage timers, because it sat inside the stage that
        looked cheap.  It is a pure function of the pool length: build once.
        """
        s = getattr(self, "_seq", None)
        if s is None:
            s = self._seq = 1.0 - 2.0 * SC.sequence(self.plan.pool_pred)
        return s

    def dummy_values(self):
        if self._dv is None:
            pl = self.plan
            dv = np.zeros(pl.pool_pred, np.complex128)
            if pl.dummy_n > 0:
                dv[pl.core_total:] = self.scrambler()[pl.core_total:]
            self._dv = dv
        return self._dv

    # -- the pool ----------------------------------------------------------
    def _sm_tables(self):
        """Frame-constant tables for the smoothed CE (built once, cached).

        E86c, ported from m9_fast._sm_tables.  Every scattered pilot of every
        symbol sits on a carrier that is a multiple of dx, so one
        (ngrid x nsym) matrix holds every raw pilot estimate and a running-sum
        boxcar over the SYMBOL axis costs O(nsym) per Frame, not O(window)
        per symbol -- which is the only reason this is affordable on RF25,
        where a Frame is 189 data symbols against RF33's 35.

        DIFFERENCE FROM m9_fast, and it matters here: RF25's Frames carry
        SBS (subframe-boundary) symbols whose pilot pattern is NOT on the dx
        grid.  m9_fast gives up on the whole Frame in that case (`usable`
        False).  Refusing to smooth 189 good symbols because a couple of
        boundary symbols have a different pattern would throw the lever away
        on exactly the carrier it was ported for, so grid-incompatible
        CLASSES are instead excluded from the smoother and stay on the
        per-symbol path.  They still decode; they just do not contribute to,
        or draw from, the boxcar.
        """
        g = self.plan.g
        dx, dy = P.dxdy(g.pattern)
        noc = P.NOC[(g.nfft, g.cred)]
        ngrid = (noc - 1) // dx + 1
        gc = np.arange(ngrid) * dx
        kk = np.arange(noc)
        jg = np.clip(np.searchsorted(gc, kk, side="right") - 1, 0, ngrid - 2)
        wtg = (kk - gc[jg]) / (gc[jg + 1] - gc[jg])
        per, ok = [], []
        for t in self.tables()["tabs"]:
            pk = t["g"]["pk"]
            fits = bool((pk % dx == 0).all()) and len(pk) > 2
            ok.append(fits)
            per.append((pk // dx) if fits else None)
        usable = ((noc - 1) % dx == 0) and any(ok)
        return dict(dx=dx, dy=dy, noc=noc, ngrid=ngrid, jg=jg, wtg=wtg,
                    gidx=per, cls_ok=ok, usable=usable)

    def _ce_smooth(self, Yb, cp):
        """-> (Hsm, cl, smooth_ok) or None when the lever is off/unusable.

        The channel SHAPE is quasi-static; the per-symbol complex GAIN is
        not.  So a scalar c_l is fitted per symbol (LS onto the frame-mean
        shape) and divided out BEFORE the boxcar, then re-applied after --
        amplitude and phase steps are followed within their own symbol while
        the boxcar averages only what is genuinely slow.  Two change
        detectors then decide, per symbol, whether the smoothing model holds
        at all; where it does not, that symbol falls back to the per-symbol
        estimate, which is the correct answer rather than a failure.
        """
        if CE_W <= 0:
            return None
        smt = getattr(self, "_smt", None)
        if smt is None:
            smt = self._smt = self._sm_tables()
        if not smt["usable"]:
            return None
        g = self.plan.g
        nsym, ngrid = g.nsym, smt["ngrid"]
        Hg = np.zeros((ngrid, nsym), np.complex128)
        Mg = np.zeros((ngrid, nsym), bool)
        hps, use = [], []
        for ci, tb in enumerate(cp["tabs"]):
            if not smt["cls_ok"][ci]:
                hps.append(None)
                continue
            rows = tb["rows"]
            hp = Yb[rows][:, tb["g"]["bins_pk"]] / tb["g"]["refpk"]
            hps.append(hp)
            use.append(ci)
            gx = smt["gidx"][ci]
            Hg[np.ix_(gx, rows)] = hp.T
            Mg[np.ix_(gx, rows)] = True
        if not use:
            return None

        cl = np.ones(nsym, np.complex128)
        w_c = np.zeros(nsym, bool)
        if CE_NORM:
            cnt_tot = np.maximum(Mg.sum(axis=1), 1)
            M = np.where(Mg, Hg, 0).sum(axis=1) / cnt_tot   # frame-mean shape
            for ci in use:
                rows = cp["tabs"][ci]["rows"]
                Mx = M[smt["gidx"][ci]]
                den = max(float(np.sum(np.abs(Mx) ** 2)), 1e-30)
                cl[rows] = (hps[ci] @ np.conj(Mx)) / den
            # a degenerate gain (a notch across the whole symbol) must not
            # blow up the normalisation; the detector handles that symbol
            w_c = np.abs(cl) < 1e-6
            cl[w_c] = 1.0
            Hg /= cl[None, :]

        Wn = CE_W
        csH = np.concatenate([np.zeros((ngrid, 1), complex),
                              np.cumsum(np.where(Mg, Hg, 0), axis=1)], axis=1)
        csM = np.concatenate([np.zeros((ngrid, 1)),
                              np.cumsum(Mg, axis=1)], axis=1)
        lo = np.maximum(np.arange(nsym) - Wn, 0)
        hi = np.minimum(np.arange(nsym) + Wn + 1, nsym)
        num = csH[:, hi] - csH[:, lo]
        cnt = csM[:, hi] - csM[:, lo]
        okc = cnt > 0
        Hsm = np.where(okc, num / np.maximum(cnt, 1), 0.0)
        covered = okc.all(axis=0)

        # RELATIVE test (E58): residual vs the frame median -- catches a
        # minority of odd symbols on an otherwise stable Frame.
        # ABSOLUTE test (E60): residual vs the symbol's OWN pilot noise
        # floor, where second differences along the carrier axis cancel the
        # channel to its curvature and leave 6 sigma^2.  The relative test
        # alone goes blind when a shape change poisons every window at once,
        # because then the median rises with the residuals.
        res_l = np.full(nsym, np.inf)
        sig2_l = np.ones(nsym)
        for ci in use:
            rows = cp["tabs"][ci]["rows"]
            gx = smt["gidx"][ci]
            hpn = hps[ci] / cl[rows][:, None]
            R = hpn - Hsm[np.ix_(gx, rows)].T
            res_l[rows] = np.mean(np.abs(R) ** 2, axis=1)
            v = hpn[:, 2:] - 2.0 * hpn[:, 1:-1] + hpn[:, :-2]
            sig2_l[rows] = np.mean(np.abs(v) ** 2, axis=1) / 6.0
        fin = np.isfinite(res_l)
        med = max(float(np.median(res_l[fin])) if fin.any() else 0.0, 1e-30)
        ratio = res_l / med
        abs_ratio = res_l / np.maximum(sig2_l, 1e-30)
        smooth_ok = (covered & (ratio <= CE_DETECT)
                     & (abs_ratio <= CE_ABS) & ~w_c & fin)
        # classes that are not on the grid never smooth
        for ci, tb in enumerate(cp["tabs"]):
            if not smt["cls_ok"][ci]:
                smooth_ok[tb["rows"]] = False
        return dict(Hsm=Hsm, cl=cl, ok=smooth_ok, smt=smt,
                    n_fallback=int(nsym - smooth_ok.sum()))

    def cell_pool(self, y, t0, report=None):
        pl, g = self.plan, self.plan.g
        cp = self.tables()
        rep = report if report is not None else {}
        pc = time.perf_counter
        t = pc()

        pre, cohs = [], []
        for s in range(g.np_sym):
            z, coh = C._demod_preamble_g(y, t0, g, s, 0)
            pre.append(FI.deinterleave(z, g.pre_nfft, s, direction="forward",
                                       toggle="i"))
            cohs.append(coh)
        pre = np.concatenate(pre) if pre else np.zeros(0, complex)
        spare = pre[g.l1b_cells + g.l1d_cells:]

        step, gi = cp["step"], cp["gi"]
        base = y[t0 + cp["data_off"] + gi:]
        need = (g.nsym - 1) * step + g.nfft
        if len(base) < need:
            raise RuntimeError(f"window short: {len(base)} < {need}")
        W = np.lib.stride_tricks.as_strided(
            base, shape=(g.nsym, g.nfft),
            strides=(step * base.strides[0], base.strides[0]))
        try:
            from scipy import fft as _sfft
            Yb = _sfft.fft(W, axis=1, workers=max(1, min(self.threads, 8)))
        except Exception:                                          # noqa: BLE001
            Yb = np.fft.fft(W, axis=1)
        Yb = np.fft.fftshift(Yb, axes=1)
        self.tm["fft"] += pc() - t
        t = pc()

        lens = cp["lens"]
        offs = np.concatenate([[0], np.cumsum(lens)]) + len(spare)
        pool = np.empty(int(offs[-1]), np.complex128)
        pool[:len(spare)] = spare
        owner = np.empty(int(offs[-1]), np.int32)
        owner[:len(spare)] = -1
        dcoh = np.empty(g.nsym)

        t_ce = pc()
        ce = self._ce_smooth(Yb, cp)
        self.tm["ce_smooth"] += pc() - t_ce

        def one_class(tb):
            sg = tb["g"]
            rows = tb["rows"]
            Yr = Yb[rows]
            hp = Yr[:, sg["bins_pk"]] / sg["refpk"]
            if ce is None:
                # untouched path -- byte-identical to pre-E86c by construction
                H = (hp[:, tb["j"]] * (1.0 - tb["wt"])
                     + hp[:, tb["j"] + 1] * tb["wt"])
            else:
                sm = ce["ok"][rows]
                jg, wtg = ce["smt"]["jg"], ce["smt"]["wtg"]
                H = np.empty((len(rows), len(sg["kk"])), np.complex128)
                if sm.any():
                    Hc = ce["Hsm"][:, rows[sm]].T
                    H[sm] = Hc[:, jg] * (1.0 - wtg) + Hc[:, jg + 1] * wtg
                    H[sm] *= ce["cl"][rows[sm]][:, None]   # re-apply the gain
                if (~sm).any():
                    hpf = hp[~sm]
                    H[~sm] = (hpf[:, tb["j"]] * (1.0 - tb["wt"])
                              + hpf[:, tb["j"] + 1] * tb["wt"])
            z = Yr[:, sg["bins_d"]] / H[:, sg["d"]]
            xg = np.empty_like(z)
            xg[np.arange(len(z))[:, None], tb["Hmat"]] = z
            num = np.abs(np.sum(hp[:, :-1] * np.conj(hp[:, 1:]), axis=1))
            den = np.maximum(np.sum(np.abs(hp) ** 2, axis=1), 1e-30)
            dcoh[rows] = num / den
            a, b = tb["trim"]
            for i, l in enumerate(rows):
                pool[offs[l]:offs[l + 1]] = xg[i, a:b]
                owner[offs[l]:offs[l + 1]] = l
        if self.ex is None:
            for tb in cp["tabs"]:
                one_class(tb)
        else:
            list(self.ex.map(one_class, cp["tabs"]))
        self.tm["equalise"] += pc() - t
        rep.update(pool=int(len(pool)), pool_predicted=pl.pool_pred,
                   preamble_coherence=[round(float(c), 4) for c in cohs],
                   data_coherence_mean=float(np.mean(dcoh)),
                   n_preamble_spare=int(len(spare)))
        # E86c: a lever that silently no-ops is worse than one that is off,
        # because the cost is still paid and the number still moves for some
        # other reason.  Report how many symbols actually took the smoothed
        # estimate, so "CE on" can never mean "CE ran and fell back 189/189".
        if ce is None:
            rep["ce"] = "off"
        else:
            rep["ce_smoothed_syms"] = int(g.nsym - ce["n_fallback"])
            rep["ce_fallback_syms"] = int(ce["n_fallback"])
        return pool, owner

    # -- CPE ---------------------------------------------------------------
    #: E53's lever, ported.  The per-symbol phase is a mean over ~6000 cells;
    #: estimating it from every 4th multiplies the estimator's noise by 2
    #: (sigma ~ 1/sqrt(N)) -- far under the channel's own phase noise -- while
    #: the CORRECTION still reaches every cell.  On RF25 the CPE is over a
    #: 1.2 M-cell pool and was 397 ms of an 807 ms Frame, so this is not a
    #: micro-optimisation; it is the difference between a tool and a tuner.
    #: Gated at the DECODED BYTES (gate_e82 leg 2b), never on ULPs.
    CPE_DECIM = 4

    def cpe_fast(self, pool, owner):
        """`cpe` with the phase estimated on a decimated subset and the total
        correction applied once.

        Two re-associations against `cpe`, both of the standing kind:
        the estimator sees every CPE_DECIM-th cell, and the three per-iteration
        divisions become one division by their product.  The hard-decision
        rule, the dummy-cell handling and the 1e-12 guard are `cpe`'s, line for
        line.
        """
        pl = self.plan
        pts = B.points_for(pl.core["mod"], pl.core["rate"])
        dv = self.dummy_values()
        n = len(pool)
        out = pool / np.sqrt(np.mean(np.abs(pool) ** 2))
        d0 = min(pl.core_total, n)
        pr = np.asarray(pts.real, np.float32)
        pi = np.asarray(pts.imag, np.float32)
        p2f = (pts.real ** 2 + pts.imag ** 2).astype(np.float32)

        idx = np.arange(0, n, self.CPE_DECIM)
        zs = out[idx].copy()
        os_ = owner[idx]
        ncore = int(np.searchsorted(idx, d0))     # idx is sorted
        hard = np.empty(len(idx), np.complex128)
        if ncore < len(idx):
            hard[ncore:] = dv[idx[ncore:]]
        syms, starts = np.unique(os_, return_index=True)
        edges = np.concatenate([starts, [len(idx)]])
        counts_s = np.diff(edges)
        ctot = np.ones(len(syms), np.complex128)

        nch = max(1, min(8, self.threads))
        step = max(1, (ncore + nch - 1) // nch)
        chunks = [(i, min(ncore, i + step)) for i in range(0, ncore, step)]

        def decide(lo, hi):
            z = zs[lo:hi]
            sc = _point_metric(np.asarray(z.real, np.float32),
                               np.asarray(z.imag, np.float32), pr, pi, p2f)
            hard[lo:hi] = pts[sc.argmin(1)]

        for _ in range(3):
            if self.ex is None or len(chunks) < 2:
                decide(0, ncore)
            else:
                list(self.ex.map(lambda c: decide(*c), chunks))
            w = np.conj(hard) * zs
            h2 = hard.real ** 2 + hard.imag ** 2
            num = np.add.reduceat(w, edges[:-1])
            den = np.add.reduceat(h2, edges[:-1])
            c = num / np.maximum(den.real, 1e-12)
            zs /= np.repeat(c, counts_s)
            ctot *= c
        # one pass over the full pool instead of three
        syms_f, starts_f = np.unique(owner, return_index=True)
        counts_f = np.diff(np.concatenate([starts_f, [n]]))
        out /= np.repeat(ctot, counts_f)
        return out

    def cpe(self, pool, owner):
        """Per-symbol residual complex gain, decision-directed.

        `m6_cells.cpe_correct`'s algorithm exactly -- 3 iterations, hard
        decision against the CORE alphabet (and against the KNOWN +-1 dummy
        values where they exist), c = <hard, z> / <hard, hard>, the same 1e-12
        guard -- with the per-symbol Python loop replaced by segment
        reductions, which is legal because `owner` is non-decreasing by
        construction of the pool.  m9_fast made the same trade for RF33 and
        the same caveat applies: the summation ORDER differs, so this is a
        rounding re-association, gated on the DECODED BYTES (gate_e82 leg 2
        diffs this path's FEC Blocks against m10_core's).
        """
        pl = self.plan
        pts = B.points_for(pl.core["mod"], pl.core["rate"])
        dv = self.dummy_values()
        n = len(pool)
        out = pool / np.sqrt(np.mean(np.abs(pool) ** 2))
        d0 = min(pl.core_total, n)
        pr = np.asarray(pts.real, np.float32)
        pi = np.asarray(pts.imag, np.float32)
        p2f = (pts.real ** 2 + pts.imag ** 2).astype(np.float32)
        syms, starts = np.unique(owner, return_index=True)
        edges = np.concatenate([starts, [n]])
        counts = np.diff(edges)
        hard = np.empty(n, np.complex128)
        if n > d0:
            hard[d0:] = dv[d0:n]
        for _ in range(3):
            zz = out[:d0]
            sc = _point_metric(np.asarray(zz.real, np.float32),
                               np.asarray(zz.imag, np.float32), pr, pi, p2f)
            hard[:d0] = pts[sc.argmin(1)]
            w = np.conj(hard) * out
            h2 = hard.real ** 2 + hard.imag ** 2
            num = np.add.reduceat(w, edges[:-1])
            den = np.add.reduceat(h2, edges[:-1])
            c = num / np.maximum(den.real, 1e-12)
            out = out / np.repeat(c, counts)
        return out

    def dummy_agreement(self, pool):
        """A/322 7.2.6.5 referee, and the reason E76 could tell a phase error
        from a dead link: it fits the KNOWN scrambler signs, so it is
        independent of the time interleaver entirely."""
        pl = self.plan
        if pl.dummy_n <= 0:
            return None
        d = pool[pl.core_total:]
        sgn = self.scrambler()[pl.core_total:len(pool)]
        rms = np.sqrt(np.mean(np.abs(d) ** 2))
        real = np.abs(d.imag) < 0.35 * rms
        if real.sum() < 50:
            return None
        return float(np.mean(np.sign(d.real[real]) == sgn[real]))

    def core_cells(self, y, t0, report=None, fast=True, dummy=True):
        """window -> the CORE PLP's cells for this Frame (complex64)."""
        pc = time.perf_counter
        rep = report if report is not None else {}
        pool, owner = self.cell_pool(y, t0, rep)
        if len(pool) != self.plan.pool_pred:
            raise RuntimeError(f"pool {len(pool)} != predicted "
                               f"{self.plan.pool_pred}")
        t = pc()
        pool = self.cpe_fast(pool, owner) if fast else self.cpe(pool, owner)
        self.tm["cpe"] += pc() - t
        # The dummy-cell referee is the instrument that told E76 a clean
        # capture from a dead one, so it stays ON in the live path -- but it
        # is a diagnostic, not a decode step, and it is timed as one.
        t = pc()
        rep["dummy_agreement"] = self.dummy_agreement(pool) if dummy else None
        self.tm["dummy_ref"] += pc() - t
        s, k = self.plan.plp_start, self.plan.plp_size
        return pool[s:s + k].astype(np.complex64)


# ===========================================================================
# 4 -- THE STREAMING CTI DE-INTERLEAVER
# ===========================================================================

class CtiStream:
    """Received cells in (one Frame at a time) -> complete FEC Blocks out.

    The index map is A/322 9.3.9.1's own, exactly as m10_cti derives it:

        q(i) = i + Nrows * ((start_row + i) mod Nrows)

    with `i` counted from the FIRST cell of the Frame whose L1 signalled
    `start_row`.  So an output cell at i needs a received cell as far ahead as
    i + Nrows*(Nrows-1) -- for RF25 that is 1,047,552 cells, 0.92 of a Frame.
    The buffer therefore holds a rolling window a little over two Frames deep
    and is trimmed from behind; nothing is ever re-read.

    A MISSING Frame is ZERO-FILLED, not skipped.  That is the property that
    makes this survivable live: the commutator's phase is a function of the
    CELL INDEX, so keeping the index space intact turns a shed Frame into
    ~35 blocks of erasure instead of a permanent loss of lock.  Skipping it
    would be exactly E76's bug, self-inflicted once per dropped Frame.
    """

    def __init__(self, plan, cap_frames=4):
        self.plan = plan
        self.nrows = plan.nrows
        self.ncell = plan.ncell
        self.cap = max(cap_frames, 3) * plan.plp_size
        self.buf = np.zeros(self.cap, np.complex64)
        self.base = 0          # absolute cell index of buf[0]
        self.end = 0           # absolute cell index one past the last written
        self.origin_frame = None
        self.next_frame = None
        self.start_row = None
        self.C0 = None
        self.out_ptr = 0       # next output cell index to emit
        self.block_id = 0
        self.stats = collections.Counter()

    def reset(self, origin_frame, start_row, fec_block_start):
        self.buf[:] = 0
        self.base = self.end = 0
        self.origin_frame = int(origin_frame)
        self.next_frame = int(origin_frame)
        self.start_row = int(start_row) % self.nrows
        self.C0 = CTI.solve_C(int(fec_block_start), self.start_row, self.nrows)
        if not (0 <= self.C0 < self.ncell):
            # m10_core gate P3 -- the signalled pair must place C inside one
            # FEC Block.  Out of range means the L1 reading is wrong, and a
            # phase built on it would be a guess wearing a spec citation.
            raise ValueError(f"C={self.C0} outside [0,{self.ncell}) -- "
                             f"start_row {self.start_row} / fec_block_start "
                             f"{fec_block_start} disagree")
        self.out_ptr = self.C0
        self.block_id = 0

    # -- input -------------------------------------------------------------
    def push(self, frame_idx, cells):
        """Append one Frame's cells.  `cells` None = a Frame we never got."""
        if self.origin_frame is None:
            return
        n = self.plan.plp_size
        while self.next_frame < frame_idx:      # zero-fill the gap
            self._append(np.zeros(n, np.complex64))
            self.next_frame += 1
            self.stats["zero_filled_frames"] += 1
        if frame_idx < self.next_frame:
            self.stats["late_frame_dropped"] += 1
            return
        if cells is None:
            self.stats["zero_filled_frames"] += 1
        self._append(np.zeros(n, np.complex64) if cells is None
                     else np.asarray(cells, np.complex64))
        self.next_frame += 1

    def _append(self, a):
        n = len(a)
        if self.end - self.base + n > self.cap:
            keep_from = max(self.base, self.out_ptr)
            k = keep_from - self.base
            if k > 0:
                m = self.end - keep_from
                self.buf[:m] = self.buf[k:k + m]
                self.base = keep_from
            if self.end - self.base + n > self.cap:
                # the consumer has fallen more than the buffer behind; the
                # honest answer is to say so, not to silently reuse cells
                drop = self.end - self.base + n - self.cap
                self.buf[:self.cap - drop] = self.buf[drop:]
                self.base += drop
                self.out_ptr = max(self.out_ptr, self.base)
                self.stats["buffer_overrun"] += 1
        o = self.end - self.base
        self.buf[o:o + n] = a
        self.end += n

    # -- output ------------------------------------------------------------
    def available_blocks(self):
        """How many complete FEC Blocks the buffered cells can produce."""
        if self.origin_frame is None:
            return 0
        reach = self.plan.cti_reach
        limit = self.end - reach            # highest output index we can fill
        return max(0, (limit - self.out_ptr) // self.ncell)

    def take_blocks(self, nmax=None):
        """-> list of (block_id, cells) for every complete Block available."""
        nb = self.available_blocks()
        if nmax is not None:
            nb = min(nb, nmax)
        if nb <= 0:
            return []
        i = np.arange(self.out_ptr, self.out_ptr + nb * self.ncell,
                      dtype=np.int64)
        q = i + self.nrows * ((self.start_row + i) % self.nrows)
        cells = self.buf[q - self.base]
        out = [(self.block_id + m,
                cells[m * self.ncell:(m + 1) * self.ncell])
               for m in range(nb)]
        self.out_ptr += nb * self.ncell
        self.block_id += nb
        return out

    def rewind_to(self, start_row, fec_block_start, keep_frames=2):
        """Re-anchor the phase WITHOUT throwing the buffered cells away.

        A candidate that fails must be replaced by the next candidate over the
        SAME air, or the sweep costs a Frame of latency per try and can never
        finish inside a viewer's patience.  The cells are already here; only
        the index map changes.
        """
        self.start_row = int(start_row) % self.nrows
        self.C0 = CTI.solve_C(int(fec_block_start), self.start_row, self.nrows)
        if not (0 <= self.C0 < self.ncell):
            return False
        lo = self.base + self.C0
        # start as far back as the buffer allows, on a Block boundary
        k = 0
        while lo + (k + 1) * self.ncell + self.plan.cti_reach <= self.end:
            k += 1
        self.out_ptr = self.base + self.C0
        self.block_id = 0
        return True


# ===========================================================================
# 5 -- THE PHASE TRACKER (the part E76 says must not need a human)
# ===========================================================================

class PhaseTracker:
    """Acquire, verify and re-acquire the CTI commutator phase, live.

    STATES
        cold        nothing known; needs an L1 (or a remembered anchor)
        probing     a candidate phase is set and is being judged by the FEC
        locked      the candidate met LOCK_FRAC over LOCK_BLOCKS

    The ONLY evidence that promotes a candidate is LDPC convergence with a
    zero BCH syndrome, because that is the only statistic a wrong phase cannot
    fake.  Pilot coherence, dummy-cell SNR and constellation tightness are all
    excellent at a wrong phase -- that is precisely the E76 trap, and it is
    why none of them appear in this decision.
    """

    LOCK_BLOCKS = 4            # blocks judged before a candidate is promoted
    LOCK_FRAC = 0.5            # and the fraction of them that must be clean
    LOSS_WINDOW = 24           # rolling window used to detect a lost lock
    LOSS_FRAC = 0.10           # under this, the lock is not a lock any more
    SWEEP_K = 6                # +-Frames of commutator advance swept blind

    def __init__(self, plan, log=print):
        self.plan = plan
        self.log = log
        self.state = "cold"
        self.anchor = None      # (frame_idx, start_row) known good
        self.fec_block_start = plan.core.get("cti_fec_block_start")
        self.cand = None        # (start_row, k) currently on probation
        self.queue = []         # remaining candidates for this acquisition
        self.probe = collections.Counter()
        self.hist = collections.deque(maxlen=self.LOSS_WINDOW)
        self.stats = collections.Counter()
        self.acq_id = 0

    # -- knowledge ---------------------------------------------------------
    def learn_from_l1(self, frame_idx, l1res):
        """The primary source: L1 decoded FROM THE FRAME WE ARE COLLECTING."""
        d = start_row_of(l1res, self.plan.core["id"])
        if "start_row" not in d:
            return False
        self.anchor = (int(frame_idx), int(d["start_row"]))
        if "fec_block_start" in d:
            self.fec_block_start = int(d["fec_block_start"])
        self.stats["l1_anchor"] += 1
        return True

    def predict(self, frame_idx):
        """A/322 9.3.9.1: start_row advances by plp_size rows per Subframe."""
        if self.anchor is None:
            return None
        f0, sr0 = self.anchor
        return int((sr0 + (frame_idx - f0) * self.plan.plp_size)
                   % self.plan.nrows)

    def candidates(self, frame_idx):
        """The phases to try at `frame_idx`, best first.

        With an anchor this is a ONE-element list and the sweep never runs --
        which is the point.  The sweep exists for the case where the anchor
        itself was lost (a bootstrap re-acquisition renumbers the Frames), and
        it sweeps exactly E76's axis: whole Frames of commutator advance.
        """
        sr = self.predict(frame_idx)
        if sr is None:
            return []
        n, step = self.plan.nrows, self.plan.plp_size
        out = [(sr, 0)]
        if self.stats["blind"]:
            for k in range(1, self.SWEEP_K + 1):
                out.append((int((sr + k * step) % n), k))
                out.append((int((sr - k * step) % n), -k))
        return out

    # -- the acquisition state machine ------------------------------------
    def begin(self, frame_idx, blind=False):
        """Start (or restart) an acquisition anchored at `frame_idx`."""
        if blind:
            self.stats["blind"] += 1
        self.acq_id += 1
        self.queue = self.candidates(frame_idx)
        self.probe.clear()
        self.hist.clear()
        if not self.queue:
            self.state = "cold"
            self.cand = None
            return None
        self.cand = self.queue.pop(0)
        self.state = "probing"
        return self.cand

    def next_candidate(self):
        self.probe.clear()
        if not self.queue:
            self.state = "cold"
            self.cand = None
            return None
        self.cand = self.queue.pop(0)
        self.stats["candidate_rejected"] += 1
        return self.cand

    def observe(self, n_blocks, n_clean):
        """Feed the syndrome verdict for a batch of Blocks.

        Returns one of: None (carry on), "promote", "reject", "lost".
        """
        if self.state == "probing":
            self.probe["n"] += n_blocks
            self.probe["ok"] += n_clean
            if self.probe["n"] < self.LOCK_BLOCKS:
                return None
            if self.probe["ok"] >= self.LOCK_FRAC * self.probe["n"]:
                self.state = "locked"
                self.hist.clear()
                self.stats["locked"] += 1
                return "promote"
            return "reject"
        if self.state == "locked":
            self.hist.append((n_blocks, n_clean))
            if len(self.hist) < self.hist.maxlen:
                return None
            tot = sum(a for a, _ in self.hist)
            ok = sum(b for _, b in self.hist)
            if tot and ok < self.LOSS_FRAC * tot:
                self.state = "cold"
                self.stats["lock_lost"] += 1
                return "lost"
        return None


# ===========================================================================
# 6 -- THE PIPELINE: Frame windows in, Baseband stream out
# ===========================================================================

def _demod_worker(shm_win, shm_cell, shm_blk, nslots, win_len, cell_len,
                  blk_len, in_q, out_q, spec, cfg):
    """One demod worker PROCESS (E85).

    WHY PROCESSES.  E52 measured it on the RF33 path and the answer has not
    changed: Python threads cannot scale this decoder because it is
    gather/scatter-heavy and NumPy holds the GIL for advanced indexing --
    1 -> 6 threads bought 36% while processes bought 2.2x.  E82 then measured
    that the LDM Frame demodulator is 230 ms of a 242 ms budget, and that it
    is a pure function of (window, t0): per-Frame independent, therefore
    divisible.  Nothing else in the LDM pipeline is -- the CTI is one
    continuous index space and stays in the parent.

    WHAT THIS WORKER IS NOT ALLOWED TO DECIDE.  It receives a Frame window and
    returns that Frame's Core-PLP cells.  It holds no interleaver state, no
    phase and no Frame counter, so it cannot get them wrong; the parent
    re-imposes FIFO order before anything reaches the CTI, which is what makes
    a pool byte-identical to the serial path (gated: gate_e85 leg 1).
    """
    import os
    import traceback
    from concurrent.futures import ThreadPoolExecutor
    from multiprocessing import shared_memory
    import numpy as np
    try:
        os.environ.setdefault("M9_NO_TORCH", "1")
        import m44_ldm as M44
        shm_w = shared_memory.SharedMemory(name=shm_win)
        shm_c = shared_memory.SharedMemory(name=shm_cell)
        shm_b = shared_memory.SharedMemory(name=shm_blk)
        win = np.ndarray((nslots, win_len), np.complex128, buffer=shm_w.buf)
        cel = np.ndarray((nslots, cell_len), np.complex64, buffer=shm_c.buf)
        blk = np.ndarray((nslots, blk_len), np.complex64, buffer=shm_b.buf)
        plan = M44.LdmPlan.from_spec(spec)
        nth = int(cfg["threads"])
        ex = ThreadPoolExecutor(max_workers=nth) if nth > 1 else None
        fd = M44.FrameDemod(plan, threads=nth, ex=ex)
        fd.prewarm()
        # A FEC worker is an LdmPipeline with no stream state at all: no CTI,
        # no phase, no Frame counter.  It is handed complete FEC Blocks and
        # returns Baseband Packets, which is a pure function -- so which
        # worker ran a Block cannot be visible in the output.
        fec = M44.LdmPipeline(plan=plan, iters=cfg["iters"], threads=nth,
                              ex=ex, accel="cpu", log=lambda *_a: None,
                              fec_threads=nth)
        fec.fd = fd
        ncell = plan.ncell
        out_q.put(("ready", (os.getpid(),
                             os.environ.get("OMP_NUM_THREADS"))))
        while True:
            msg = in_q.get()
            if msg is None:
                break
            kind = msg[0]
            t = time.perf_counter()
            if kind == "demod":
                _k, wslot, cslot, idx, n, t0 = msg
                rep = {}
                ok = True
                try:
                    cells = fd.core_cells(win[wslot, :n], t0, rep)
                    cel[cslot, :len(cells)] = cells
                except Exception as e:                         # noqa: BLE001
                    ok = False
                    rep["err"] = f"{type(e).__name__}: {e}"
                out_q.put(("cells", (idx, wslot, cslot, ok,
                                     dict(coh=rep.get("data_coherence_mean"),
                                          dummy=rep.get("dummy_agreement"),
                                          err=rep.get("err")),
                                     time.perf_counter() - t)))
            elif kind == "fec":
                _k, bslot, first_id, nb = msg
                segs = [blk[bslot, m * ncell:(m + 1) * ncell]
                        for m in range(nb)]
                try:
                    pk, nconv, nbch = fec.decode_blocks(segs)
                    err = None
                except Exception as e:                         # noqa: BLE001
                    pk, nconv, nbch = [None] * nb, 0, 0
                    err = f"{type(e).__name__}: {e}"
                out_q.put(("fec", (first_id, bslot, nb, pk, nconv, nbch, err,
                                   time.perf_counter() - t)))
        tm = dict(fd.tm)
        tm.update({("w_" + k): v for k, v in fec.tm.items()})
        out_q.put(("done", tm))
        if ex is not None:
            ex.shutdown(wait=False)
        shm_w.close()
        shm_c.close()
        shm_b.close()
    except Exception as e:                                     # noqa: BLE001
        try:
            out_q.put(("err", f"{type(e).__name__}: {e}\n"
                              f"{traceback.format_exc()}"))
        except Exception:                                      # noqa: BLE001
            pass


class LdmDemodPool:
    """Frame windows in by shared memory, Core-PLP cells out, STRICTLY ORDERED.

    Two shared arenas, because the two directions are different sizes and
    different dtypes: a Frame WINDOW is 1,684,520 complex128 (27 MB) and a
    Frame's Core cells are 1,133,282 complex64 (9 MB).  Both are memcpy at
    ~10 GB/s, so the round trip is ~4 ms of a 242 ms budget -- the reason
    windows travel by shm and not by queue.

    ORDER IS NOT OPTIONAL HERE, and that is the difference from the RF33 pool.
    There, out-of-order completion only costs the Baseband stream's ordering,
    which the parent restores.  Here a Frame IS a position in the CTI's index
    space: hand the interleaver Frame 7 before Frame 6 and the commutator is
    wrong from that point on, permanently.  So results park in a dict and are
    released only from the head of the dispatch deque -- and a Frame that
    FAILED to demodulate is still released, as a None, so the CtiStream can
    zero-fill it and keep the index space intact (E82's rule: a shed Frame is
    a zero-fill, never a skip).
    """

    #: FEC Blocks per dispatched batch.  A Frame yields ~35, and the batch is
    #: what the LDPC's own batching wants anyway; larger batches are split.
    MAX_BLK = 48

    def __init__(self, plan, nproc=4, threads=4, blas=1, iters=50, log=print):
        self.plan = plan
        self.nproc = int(nproc)
        self.threads = int(threads)
        self.blas = int(blas)
        self.iters = int(iters)
        self.log = log
        self.win_len = plan.frame_window + 4096
        self.cell_len = plan.plp_size
        self.blk_len = self.MAX_BLK * plan.ncell
        self.nslots = self.nproc + 2
        self.procs = []
        self.free_w = list(range(self.nslots))
        self.free_c = list(range(self.nslots))
        self.free_b = list(range(self.nslots))
        self.order = collections.deque()
        self.parked = {}
        self.forder = collections.deque()
        self.fparked = {}
        self.stats = collections.Counter()
        self.tm = collections.Counter()
        self.started = False

    # -- lifecycle ---------------------------------------------------------
    def start(self):
        import multiprocessing as mp
        from multiprocessing import shared_memory
        # PIN BLAS IN THE CHILDREN, FROM THE PARENT, BEFORE THEY EXIST.
        #
        # `spawn` hands the child the parent's environment, and the child
        # imports NumPy during that spawn -- so the only moment the setting
        # can take effect is now.  Setting it inside the worker body is too
        # late (NumPy is already up) and setting it in a launcher script only
        # fixes the copy of the process that remembered to.
        #
        # It is not a tuning knob, it is a correctness one.  E82 deadlocked
        # this exact chain twice by letting a multi-threaded BLAS call run
        # inside its own thread pool, and the Ubuntu port measured a 12-thread
        # OpenBLAS pool making an audio worker 100x SLOWER on a loop of tiny
        # matmuls -- "a BLAS pool under a Python-loop of tiny ops is a
        # THROTTLE dressed as parallelism".  The fan-out here is processes and
        # our own thread pool; BLAS gets one thread and no opinions.
        for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                  "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
            os.environ[v] = str(self.blas)
        self.shm_w = shared_memory.SharedMemory(
            create=True, size=self.nslots * self.win_len * 16)
        self.shm_c = shared_memory.SharedMemory(
            create=True, size=self.nslots * self.cell_len * 8)
        self.shm_b = shared_memory.SharedMemory(
            create=True, size=self.nslots * self.blk_len * 8)
        self.win = np.ndarray((self.nslots, self.win_len), np.complex128,
                              buffer=self.shm_w.buf)
        self.cel = np.ndarray((self.nslots, self.cell_len), np.complex64,
                              buffer=self.shm_c.buf)
        self.blk = np.ndarray((self.nslots, self.blk_len), np.complex64,
                              buffer=self.shm_b.buf)
        ctx = mp.get_context("spawn")
        self.in_q = ctx.Queue()
        self.out_q = ctx.Queue()
        cfg = dict(threads=self.threads, iters=self.iters)
        spec = self.plan.to_spec()
        self.procs = [ctx.Process(target=_demod_worker,
                                  args=(self.shm_w.name, self.shm_c.name,
                                        self.shm_b.name, self.nslots,
                                        self.win_len, self.cell_len,
                                        self.blk_len, self.in_q, self.out_q,
                                        spec, cfg), daemon=True)
                      for _ in range(self.nproc)]
        t0 = time.time()
        for p in self.procs:
            p.start()
        ready = 0
        omp = None
        while ready < self.nproc:
            kind, payload = self.out_q.get(timeout=300)
            if kind == "err":
                raise RuntimeError(f"demod worker failed to start: {payload}")
            if kind == "ready":
                ready += 1
                omp = payload[1]
        self.started = True
        self.log(f"  LDM demod pool: {self.nproc} processes x {self.threads} "
                 f"threads, BLAS pinned to {omp}, {self.nslots} shared "
                 f"window slots ({self.nslots * self.win_len * 16 >> 20} MB) "
                 f"+ cell slots ({self.nslots * self.cell_len * 8 >> 20} MB), "
                 f"warmed in {time.time() - t0:.1f} s")

    def stop(self):
        if not self.procs:
            return
        for _ in self.procs:
            try:
                self.in_q.put(None)
            except Exception:                                  # noqa: BLE001
                break
        deadline = time.time() + 8
        got = 0
        while got < len(self.procs) and time.time() < deadline:
            try:
                kind, payload = self.out_q.get(timeout=0.5)
            except Exception:                                  # noqa: BLE001
                continue
            if kind == "done":
                self.tm.update(payload)
                got += 1
        for p in self.procs:
            p.join(timeout=5)
            if p.is_alive():
                p.terminate()
        self.procs = []
        for s in (self.shm_w, self.shm_c, self.shm_b):
            try:
                s.close()
                s.unlink()
            except Exception:                                  # noqa: BLE001
                pass

    # -- dispatch / collect ------------------------------------------------
    def in_flight(self):
        return len(self.order)

    def can_submit(self):
        return bool(self.free_w) and bool(self.free_c)

    def submit(self, idx, window, t0):
        """-> True if dispatched.  The caller MUST retry a False."""
        if not self.can_submit():
            return False
        w = self.free_w.pop()
        c = self.free_c.pop()
        n = len(window)
        if n > self.win_len:
            raise ValueError(f"window {n} > slot {self.win_len}")
        self.win[w, :n] = window
        self.in_q.put(("demod", w, c, int(idx), int(n), int(t0)))
        self.order.append(int(idx))
        self.stats["submitted"] += 1
        return True

    def _pump(self, block, timeout=0.5):
        try:
            kind, payload = (self.out_q.get(timeout=timeout) if block
                             else self.out_q.get_nowait())
        except Exception:                                      # noqa: BLE001
            return False
        if kind == "err":
            raise RuntimeError(f"demod worker: {payload}")
        if kind == "done":
            self.tm.update(payload)
            return False
        if kind == "cells":
            idx, wslot, cslot, ok, diag, dt_ = payload
            self.free_w.append(wslot)
            self.parked[idx] = (cslot, ok, diag, dt_)
            self.stats["returned"] += 1
        elif kind == "fec":
            first_id, bslot, nb, pk, nconv, nbch, err, dt_ = payload
            self.free_b.append(bslot)
            self.fparked[first_id] = (nb, pk, nconv, nbch, err, dt_)
            self.stats["fec_returned"] += 1
        return True

    # -- FEC Blocks --------------------------------------------------------
    def can_submit_fec(self):
        return bool(self.free_b)

    def submit_fec(self, first_id, segs):
        """Dispatch a contiguous run of FEC Blocks.  -> True if dispatched.

        Blocks are independent given their cells, so this is the same kind of
        job as a Frame demod: a pure function with no stream state.  What is
        NOT independent is their ORDER in the Baseband stream, so batches are
        released from a deque head exactly as Frames are.
        """
        if not self.free_b or not segs:
            return False
        nb = len(segs)
        if nb > self.MAX_BLK:
            raise ValueError(f"batch {nb} > MAX_BLK {self.MAX_BLK}")
        b = self.free_b.pop()
        ncell = self.plan.ncell
        for m, seg in enumerate(segs):
            self.blk[b, m * ncell:(m + 1) * ncell] = seg
        self.in_q.put(("fec", b, int(first_id), int(nb)))
        self.forder.append(int(first_id))
        self.stats["fec_submitted"] += 1
        return True

    def fec_in_flight(self):
        return len(self.forder)

    def collect_fec(self, on_batch, block=False, drain=False):
        """Release completed FEC batches IN ORDER: on_batch(first_id, packets,
        nconv, nbch, err)."""
        n = 0
        if block and self.forder:
            self._pump(True)
        while self._pump(False):
            pass
        if drain:
            deadline = time.time() + 30
            while self.forder and self.forder[0] not in self.fparked \
                    and time.time() < deadline:
                self._pump(True)
                while self._pump(False):
                    pass
        while self.forder and self.forder[0] in self.fparked:
            fid = self.forder.popleft()
            nb, pk, nconv, nbch, err, dt_ = self.fparked.pop(fid)
            self.tm["worker_fec"] += dt_
            on_batch(fid, pk, nconv, nbch, err)
            n += 1
        return n

    def discard_fec(self):
        """Throw away every FEC batch in flight.  Used when the phase is lost:
        those Blocks were cut from an index space we no longer believe."""
        deadline = time.time() + 20
        while self.forder and time.time() < deadline:
            self._pump(True)
            while self._pump(False):
                pass
            while self.forder and self.forder[0] in self.fparked:
                fid = self.forder.popleft()
                self.fparked.pop(fid)
                self.stats["fec_discarded"] += 1

    def collect(self, on_frame, block=False, drain=False):
        """Release completed Frames IN ORDER, calling `on_frame(idx, cells,
        diag)` while the shared slot is still owned -- so the CtiStream copies
        straight out of shared memory and the slot is freed immediately after,
        with no intermediate 9 MB copy.

        `cells` is None for a Frame the worker could not demodulate; the
        caller must zero-fill it, never skip it.
        """
        n = 0
        if block and self.order:
            self._pump(True)
        while self._pump(False):
            pass
        if drain:
            deadline = time.time() + 30
            while self.order and self.order[0] not in self.parked \
                    and time.time() < deadline:
                self._pump(True)
                while self._pump(False):
                    pass
        while self.order and self.order[0] in self.parked:
            idx = self.order.popleft()
            cslot, ok, diag, dt_ = self.parked.pop(idx)
            self.tm["worker_demod"] += dt_
            try:
                on_frame(idx, self.cel[cslot, :self.cell_len] if ok else None,
                         diag)
            finally:
                self.free_c.append(cslot)
            n += 1
        return n


class LdmPipeline:
    """`m9_fast.FrameDecoder`'s role for an LDM/CTI multiplex.

    The interface m11_watch needs is deliberately the same shape as the RF33
    one -- push a Frame, get Baseband bytes -- but the INSIDE cannot be, and
    the reason is the interleaver.  RF33's Frames are independent, so its
    decoder is a pure function of (window, t0).  Here a FEC Block is a
    diagonal across ~two Frames of a single continuous index space, so this
    object is a STREAM with state: a phase, a cell buffer, and a lock.

    Output is `(bytes, bounds)` per call -- exactly what `m11_watch._decode_done`
    already builds and hands to the transport, so everything downstream of
    this class (ALP walker, IP reassembly, ROUTE/MMTP, the lanes, the viewer)
    is the RF33 chain, untouched.
    """

    #: how many Frames to keep trying L1 on before admitting we cannot anchor
    L1_TRIES = 40
    #: Frames between free P4 continuity checks (m10_core's cross-Frame gate,
    #: live).  An L1 decode is ~30-600 ms depending on whether it converges,
    #: so this is a real cost: at 100 it is 6 ms/Frame amortised and still
    #: samples the commutator identity every 24 seconds.  0 disables it.
    L1_CHECK_EVERY = 100

    def __init__(self, plan=None, iters=50, threads=4, ex=None, accel="cpu",
                 log=print, alpha=None, procs=0, proc_threads=4,
                 fec_threads=None):
        self.plan = plan
        self.iters = iters
        self.threads = threads
        self.ex = ex
        self.accel = accel
        self.log = log
        self.alpha = alpha
        # E85: `procs` > 0 moves the Frame demodulator into worker PROCESSES.
        # Everything else -- the CTI, the phase, the FEC -- stays here, and
        # the decoded bytes are identical either way (gate_e85 leg 1).
        self.procs = int(procs)
        self.proc_threads = int(proc_threads)
        self.fec_threads = int(fec_threads or threads)
        self.pool = None
        self.fd = None
        self.cti = None
        self.ph = None
        self.gpu = None
        self.origin = None
        self.l1_tries = 0
        self.bb_abs = 0        # absolute byte offset into the Baseband stream
        self.n_frames = 0
        self.n_blocks = 0
        self.n_conv = 0
        self.n_bch = 0
        self.last_l1 = None
        self.tm = collections.Counter()
        self.stats = collections.Counter()
        self.diag = collections.deque(maxlen=64)

    # -- setup -------------------------------------------------------------
    def adopt(self, plan):
        self.plan = plan
        self.fd = FrameDemod(plan, threads=self.threads, ex=self.ex)
        self.cti = CtiStream(plan)
        self.ph = PhaseTracker(plan, log=self.log)
        self.log("  LDM plan adopted --\n    " + plan.describe())
        # On an LDM carrier the LDPC runs inside the demod POOL, and pool
        # workers cannot reach the GPU -- they use the compiled C kernel or
        # the numpy fallback. Measured 8/15, same carrier, same antenna, six
        # minutes apart: fallback 0.04x, kernel 1.01x. So here the kernel
        # is not "optional, ~7x on one stage"; it is the difference between
        # television and a stopped clock, and the one-line stderr warning
        # from m3_ldpc is not proportionate to that. Say it where it bites.
        if not LD._kernel_lib():
            self.log("  ** LDM carrier WITHOUT the compiled LDPC kernel: expect "
                     "~0.04x real time (measured), i.e. no watchable picture. "
                     "Build it once per box: python lab/build_ldpc_kernel.py")

    def prewarm(self):
        if self.fd is not None:
            self.fd.prewarm()
            if self.accel in ("gpu", "gpu-full"):
                self._make_gpu()
        if self.procs > 0 and self.pool is None and self.plan is not None:
            self.pool = LdmDemodPool(self.plan, nproc=self.procs,
                                     threads=self.proc_threads,
                                     iters=self.iters, log=self.log)
            self.pool.start()

    def close(self):
        if self.pool is not None:
            self.pool.stop()
            self.pool = None

    def _make_gpu(self):
        try:
            import m9_gpu
            if not m9_gpu.available():
                self.log("  GPU LDPC unavailable -- CPU min-sum (USER LAW: "
                         "the GPU is a lever, never a requirement)")
                return
            ch = self.fd.chain(self.iters)
            self.gpu = m9_gpu.GpuMinSum(ch.checks, ch.ninner, iters=self.iters,
                                        dtype="float32")
            self.log(f"  GPU LDPC engaged for Ninner {ch.ninner}")
        except Exception as e:                                     # noqa: BLE001
            self.log(f"  GPU LDPC unavailable ({type(e).__name__}: {e}) -- CPU")

    def reset(self):
        """A sample-continuity break: the Frame count and the phase both die."""
        if self.cti is not None:
            self.cti.origin_frame = None
        if self.ph is not None:
            self.ph.state = "cold"
            self.ph.anchor = None
        self.origin = None
        self.l1_tries = 0
        self.stats["pipeline_reset"] += 1

    # -- one Frame ---------------------------------------------------------
    def push_frame(self, frame_idx, y, t0, ps=None):
        """-> list of (baseband bytes, ALP bound offsets) ready for transport."""
        pc = time.perf_counter
        # 1 -- L1, if we still need a phase.  Decoded from THIS window, so it
        #      describes THIS Frame; that identity is the whole E76 fix.
        # THE CONDITION IS "no live cell stream", not "no anchor".  A rejected
        # candidate leaves the anchor in place on purpose (it is the only
        # thing a sweep can be centred on), so keying this branch off the
        # anchor would leave the pipeline with a dead CtiStream and no path
        # that ever re-decodes L1 -- a wedge, not a fallback.
        need_phase = (self.plan is None or self.ph is None or self.cti is None
                      or self.cti.origin_frame is None)
        if need_phase:
            t = pc()
            L = l1_from_window(y, t0, ps=ps)
            self.tm["l1"] += pc() - t
            self.l1_tries += 1
            if L is not None and L.get("ok"):
                if self.plan is None:
                    self.adopt(LdmPlan.from_l1_result(L, label="live"))
                self.last_l1 = L
                if self.ph.learn_from_l1(frame_idx, L):
                    self.stats["l1_verified"] += 1
                    self._begin(frame_idx)
            elif (self.l1_tries >= self.L1_TRIES and self.plan is not None
                  and self.ph.anchor is not None):
                # L1 will not verify but the Frame count is still continuous,
                # so the dead-reckoned start_row is still meaningful and only
                # its Frame OFFSET is in doubt.  Sweep it.  (Since the 6.5.2.7
                # fix L1 verifies on every Frame of both banked captures, so
                # this path is REACHABLE BUT UNTESTED ON AIR -- said plainly.)
                self.stats["l1_gave_up"] += 1
                self.l1_tries = 0
                self._begin(frame_idx, blind=True)
        elif self.L1_CHECK_EVERY and frame_idx % self.L1_CHECK_EVERY == 0:
            # FREE CONTINUITY GATE (m10_core P4, live): whenever L1 verifies
            # again, its signalled start_row must equal the dead-reckoned one.
            # A disagreement means the Frame counter and the commutator have
            # drifted apart -- which is exactly the failure this class exists
            # to prevent, so it is checked continuously rather than assumed.
            t = pc()
            L = l1_from_window(y, t0, ps=ps)
            self.tm["l1"] += pc() - t
            if L is not None and L.get("ok"):
                d = start_row_of(L, self.plan.core["id"])
                pred = self.ph.predict(frame_idx)
                if "start_row" in d:
                    if int(d["start_row"]) == pred:
                        self.stats["p4_agree"] += 1
                        self.ph.anchor = (int(frame_idx), int(d["start_row"]))
                    else:
                        self.stats["p4_disagree"] += 1
                        self.log(f"  *** L1 says start_row "
                                 f"{d['start_row']} at Frame {frame_idx}, "
                                 f"dead reckoning says {pred} -- "
                                 f"re-anchoring ***")
                        self.ph.anchor = (int(frame_idx), int(d["start_row"]))
                        self.ph.fec_block_start = int(d["fec_block_start"])
                        self._begin(frame_idx)

        if self.plan is None or self.cti.origin_frame is None:
            return []

        if self.pool is not None:
            return self._push_pooled(frame_idx, y, t0)

        # 2 -- demodulate this Frame's core cells
        t = pc()
        rep = {}
        try:
            cells = self.fd.core_cells(y, t0, rep)
        except Exception as e:                                     # noqa: BLE001
            self.stats["demod_failed"] += 1
            self.log(f"  Frame {frame_idx}: demod failed "
                     f"({type(e).__name__}: {e}) -- zero-filled")
            cells = None
        self.tm["demod"] += pc() - t
        self._frame_done(frame_idx, dict(coh=rep.get("data_coherence_mean"),
                                         dummy=rep.get("dummy_agreement")))

        # 3 -- the CTI, and the FEC
        self.cti.push(frame_idx, cells)
        return self.drain()

    # -- E85: the same Frame, demodulated in a worker process --------------
    def _frame_done(self, frame_idx, diag):
        self.n_frames += 1
        self.diag.append(dict(frame=frame_idx, coh=diag.get("coh"),
                              dummy=diag.get("dummy")))

    def _push_pooled(self, frame_idx, y, t0):
        """Dispatch this Frame, then release whatever finished, IN ORDER.

        The release callback is where the CTI is fed, and it runs from the
        head of the dispatch deque only -- so the interleaver sees Frames in
        exactly the order the serial path fed them, which is what makes the
        two paths byte-identical.  A Frame the worker could not demodulate
        arrives as None and is ZERO-FILLED by CtiStream, never skipped: the
        commutator's phase is a function of the cell INDEX, and skipping is
        how you lose it.
        """
        pc = time.perf_counter
        out = []

        def on_frame(idx, cells, diag):
            if cells is None:
                self.stats["demod_failed"] += 1
                if diag.get("err"):
                    self.log(f"  Frame {idx}: demod failed in worker "
                             f"({diag['err']}) -- zero-filled")
            self._frame_done(idx, diag)
            self.cti.push(idx, cells)
            out.extend(self.drain())

        # back-pressure: a full pool is the rate limit, and waiting here is
        # correct -- the front end's own queue absorbs it and sheds a FRAME if
        # it must, which costs 242 ms of picture and never a sample.
        t = pc()
        deadline = time.time() + 30
        while not self.pool.can_submit():
            self.pool.collect(on_frame, block=True)
            if time.time() > deadline:
                self.stats["pool_stalled"] += 1
                break
        self.tm["pool_wait"] += pc() - t
        if self.pool.can_submit():
            self.pool.submit(frame_idx, y, t0)
        else:
            # never observed; say so rather than silently losing the Frame
            self.log(f"  *** demod pool did not free a slot in 30 s -- "
                     f"Frame {frame_idx} dropped ***")
        self.pool.collect(on_frame)
        return out

    def flush(self):
        """Release every Frame still inside the pool.  End of stream only."""
        if self.pool is None:
            return []
        out = []

        def on_frame(idx, cells, diag):
            self._frame_done(idx, diag)
            self.cti.push(idx, cells)
            out.extend(self.drain())
        while self.pool.in_flight():
            if self.pool.collect(on_frame, block=True, drain=True) == 0:
                break

        def on_batch(first_id, pkts, nconv, nbch, err):
            nb = len(pkts)
            self.n_blocks += nb
            self.n_conv += nconv
            self.n_bch += nbch
            self.ph.observe(nb, nbch)
            seg = self._packets_to_stream(pkts)
            if seg is not None:
                out.append(seg)
        while self.pool.fec_in_flight():
            if self.pool.collect_fec(on_batch, block=True, drain=True) == 0:
                break
        return out

    def _begin(self, frame_idx, blind=False):
        """Adopt the next phase candidate, skipping any that fail the P3 gate.

        E86.  CtiStream enforces "C must land inside one FEC Block" in two
        places and, until now, with two different consequences: rewind_to()
        returned False and the sweep carried on, while reset() raised -- so a
        candidate rejected DURING a sweep was routine and the identical
        candidate rejected while STARTING one killed the receiver.  The blind
        re-begin path reaches reset() with a dead-reckoned start_row and a
        fec_block_start from an older anchor, and those two disagreeing is a
        normal event on a marginal signal, not a corrupt program state.

        The consequence was backwards: the receiver died precisely when the
        air was bad enough to need re-acquisition.  A candidate that fails the
        range gate is now what it always was in the sweep -- REJECTED, counted,
        and replaced by the next one.  The gate keeps its authority (no phase
        is ever built on an out-of-range C); it just no longer takes the whole
        run down with it, and it no longer truncates the evidence of the fade
        that produced it.
        """
        cand = self.ph.begin(frame_idx, blind=blind)
        if cand is None:
            return
        rejected = 0
        while cand is not None:
            sr, k = cand
            try:
                self.cti.reset(frame_idx, sr, self.ph.fec_block_start)
            except ValueError as e:
                rejected += 1
                self.stats["phase_rejected_range"] += 1
                self.log(f"  CTI phase candidate start_row {sr} (k={k:+d}) "
                         f"fails the range gate before use -- {e}")
                cand = self.ph.next_candidate()
                continue
            self.origin = frame_idx
            self.log(f"  CTI phase candidate: Frame {frame_idx}, start_row "
                     f"{sr} (k={k:+d}), C={self.cti.C0} -- PROVISIONAL until "
                     f"{self.ph.LOCK_BLOCKS} FEC Blocks judge it")
            return
        # every candidate the tracker could offer was out of range: this is
        # the same dead end the FEC-rejection path reaches, so end in the same
        # state rather than inventing a second one.
        self.log(f"  *** all {rejected} CTI phase candidates fail the range "
                 f"gate -- cold, waiting for a fresh L1 ***")
        self.cti.origin_frame = None
        self.l1_tries = 0

    def drain(self):
        """CTI -> FEC -> Baseband, with the FEC either here or in the pool.

        THE SPLIT IS DELIBERATE AND IT IS ABOUT AUTHORITY, not speed.  While
        the phase is PROVISIONAL the FEC verdict is the thing that promotes or
        rejects a candidate, and the tracker has to be able to rewind the
        CtiStream over the same air before any more Blocks are cut.  A batch
        in flight during that decision would be judging a phase that no longer
        exists.  So probing is SYNCHRONOUS, here, exactly as it has always
        been -- a handful of Blocks, once per acquisition -- and only a LOCKED
        phase goes asynchronous.  E82's rule about coordinate systems applies
        to time as well as to offsets.
        """
        if self.pool is not None and self.ph.state == "locked":
            return self._drain_pooled()
        out = []
        blocks = self.cti.take_blocks()
        if not blocks:
            return out
        pc = time.perf_counter
        t = pc()
        pkts, nconv, nbch = self.decode_blocks([b for _, b in blocks])
        self.tm["fec"] += pc() - t
        self.n_blocks += len(blocks)
        self.n_conv += nconv
        self.n_bch += nbch
        verdict = self.ph.observe(len(blocks), nbch)
        if verdict == "reject":
            self.stats["phase_rejected"] += 1
            nxt = self.ph.next_candidate()
            if nxt is None:
                self.log("  *** every CTI phase candidate rejected by the FEC "
                         "-- cold, waiting for a fresh L1 ***")
                self.cti.origin_frame = None
                # the anchor is NOT cleared: a fresh L1 overrides it anyway,
                # and until one arrives it is the only thing a sweep can be
                # centred on
                self.l1_tries = 0
            else:
                sr, k = nxt
                self.log(f"  CTI phase rejected ({nbch}/{len(blocks)} Blocks "
                         f"clean) -- next candidate start_row {sr} (k={k:+d})")
                # E86: rewind_to reports an out-of-range C by returning False,
                # and ignoring that left the stream anchored on a C the P3 gate
                # had just refused -- a rejected phase still cutting Blocks.
                # A candidate that cannot be anchored is not a candidate.
                if not self.cti.rewind_to(sr, self.ph.fec_block_start):
                    self.stats["phase_rejected_range"] += 1
                    self.log("  *** that candidate is out of range too -- "
                             "cold, waiting for a fresh L1 ***")
                    self.cti.origin_frame = None
                    self.l1_tries = 0
            return []
        if verdict == "promote":
            self.log(f"  *** CTI PHASE LOCKED: start_row {self.cti.start_row}, "
                     f"C={self.cti.C0}, {nbch}/{len(blocks)} Blocks BCH-clean "
                     f"***")
        if verdict == "lost":
            # A lost lock with the Frame count INTACT is the case the blind
            # sweep exists for, and it is reachable: the front end is still
            # tracking, so the dead-reckoned start_row is still meaningful and
            # only its Frame OFFSET is in doubt.  Sweep +-SWEEP_K Frames of
            # commutator advance around it -- E76's axis, automated -- and let
            # the FEC choose.  The anchor is deliberately KEPT; throwing it
            # away would leave nothing to sweep around and make the next
            # verified L1 the only way back.
            self.log("  *** CTI lock LOST (FEC collapsed over the last "
                     f"{self.ph.LOSS_WINDOW} batches) -- sweeping "
                     f"+-{self.ph.SWEEP_K} Frames of commutator advance ***")
            self.stats["lock_lost"] += 1
            self._begin(self.cti.next_frame, blind=True)
            return []
        # 4 -- Baseband Packets -> the transport's own format.
        #
        # THE ANCHORS ARE ABSOLUTE.  `StreamAlpWalker.feed(chunk, bounds_abs)`
        # indexes into the WHOLE Baseband stream, not into this batch -- the
        # first build emitted per-batch offsets and the walker, quite
        # correctly, then refused to parse anything at all: 1716 FEC Blocks
        # BCH-clean and zero segments out.  A correct decode with a wrong
        # coordinate system produces silence, not an error.
        import m6_payload as P6
        stream, bounds = bytearray(), []
        for p in pkts:
            if p is None:
                continue
            ptr, pay = P6.bb_split(p)
            if ptr is not None:
                bounds.append(self.bb_abs + len(stream) + ptr)
            stream += pay
        if stream:
            self.bb_abs += len(stream)
            out.append((bytes(stream), bounds))
        return out

    # -- E85: the locked-phase FEC, in the pool ----------------------------
    def _packets_to_stream(self, pkts):
        import m6_payload as P6
        stream, bounds = bytearray(), []
        for p in pkts:
            if p is None:
                continue
            ptr, pay = P6.bb_split(p)
            if ptr is not None:
                bounds.append(self.bb_abs + len(stream) + ptr)
            stream += pay
        if not stream:
            return None
        self.bb_abs += len(stream)
        return (bytes(stream), bounds)

    def _drain_pooled(self):
        pc = time.perf_counter
        out = []
        lost = [False]

        def on_batch(first_id, pkts, nconv, nbch, err):
            if err:
                self.stats["worker_fec_error"] += 1
                self.log(f"  FEC batch at Block {first_id} failed in worker: "
                         f"{err}")
            nb = len(pkts)
            self.n_blocks += nb
            self.n_conv += nconv
            self.n_bch += nbch
            v = self.ph.observe(nb, nbch)
            if v == "lost":
                lost[0] = True
                return
            seg = self._packets_to_stream(pkts)
            if seg is not None:
                out.append(seg)

        # dispatch everything the CTI can cut, in Block order
        t = pc()
        blocks = self.cti.take_blocks()
        self.tm["cti_gather"] += pc() - t
        t = pc()
        i = 0
        while i < len(blocks):
            chunk = blocks[i:i + self.pool.MAX_BLK]
            deadline = time.time() + 30
            while not self.pool.can_submit_fec():
                self.pool.collect_fec(on_batch, block=True)
                if lost[0] or time.time() > deadline:
                    break
            if lost[0]:
                break
            if not self.pool.submit_fec(chunk[0][0], [c for _b, c in chunk]):
                self.stats["fec_submit_failed"] += 1
                break
            i += len(chunk)
        self.tm["fec_dispatch"] += pc() - t
        t = pc()
        self.pool.collect_fec(on_batch)
        self.tm["fec_collect"] += pc() - t
        if lost[0]:
            # Blocks still in flight were cut from an index space we have just
            # stopped believing in.  Throw them away rather than emit them:
            # a Baseband stream spliced across a phase change is worse than a
            # gap, because the ALP walker cannot tell it happened.
            self.pool.discard_fec()
            self.log("  *** CTI lock LOST (FEC collapsed over the last "
                     f"{self.ph.LOSS_WINDOW} batches) -- sweeping "
                     f"+-{self.ph.SWEEP_K} Frames of commutator advance ***")
            self.stats["lock_lost"] += 1
            self._begin(self.cti.next_frame, blind=True)
            return out
        return out

    # -- FEC ---------------------------------------------------------------
    def decode_blocks(self, cell_list):
        """A batch of FEC Blocks -> Baseband Packet bytes (None where lost)."""
        ch = self.fd.chain(self.iters)
        n = len(cell_list)
        if not n:
            return [], 0, 0
        pc = time.perf_counter
        t = pc()
        cells = np.asarray(cell_list, np.complex128)
        self.tm["fec_widen"] += pc() - t
        t = pc()
        q = demap_batch_qam(cells, ch, self.ex, self.fec_threads)
        self.tm["fec_demap"] += pc() - t
        t = pc()
        lam = np.empty((n, ch.ninner), np.float32)
        lam[:, ch.lam_of_q] = q
        self.tm["fec_bitdeint"] += pc() - t
        t = pc()
        bits, conv, bad = ldpc_batch(lam, ch, self.gpu, self.iters, self.ex,
                                     self.fec_threads, alpha=self.alpha)
        self.tm["fec_ldpc"] += pc() - t
        t = pc()
        import m9_accel
        nouter = bits[:, :ch.kldpc]
        syn = m9_accel.bch_syndrome_batch(nouter, ch.ninner)
        self.tm["fec_bch"] += pc() - t
        t = pc()
        seq = SC.sequence(ch.kpayload)
        pay = nouter[:, :ch.kpayload] ^ seq[None, :]
        by = np.packbits(pay, axis=1)
        self.tm["fec_descramble"] += pc() - t
        out, nconv, nbch = [], 0, 0
        for i in range(n):
            if not conv[i]:
                out.append(None)
                continue
            nconv += 1
            if syn[i] == 0:
                nbch += 1
                out.append(by[i].tobytes())
            else:
                out.append(None)
        return out, nconv, nbch

    def summary(self):
        s = dict(frames=self.n_frames, blocks=self.n_blocks,
                 converged=self.n_conv, bch_zero=self.n_bch,
                 stats=dict(self.stats), phase=dict(self.ph.stats)
                 if self.ph else {},
                 cti=dict(self.cti.stats) if self.cti else {},
                 state=self.ph.state if self.ph else "cold",
                 start_row=self.cti.start_row if self.cti else None,
                 C=self.cti.C0 if self.cti else None,
                 origin_frame=self.origin)
        nf = max(self.n_frames, 1)
        s["ms_per_frame"] = {k: round(1000 * v / nf, 2)
                             for k, v in self.tm.items()}
        s["demod_ms_per_frame"] = {k: round(1000 * v / nf, 2)
                                   for k, v in (self.fd.tm.items()
                                                if self.fd else [])}
        if self.pool is not None:
            s["pool"] = dict(procs=self.pool.nproc,
                             threads=self.pool.threads,
                             blas=self.pool.blas,
                             stats=dict(self.pool.stats))
            # worker_demod is WALL INSIDE THE WORKERS summed over all of
            # them, so it is the cost the pool is hiding, not the cost the
            # Frame pays.  Reported as both, because confusing them is how a
            # pool gets credited with work it merely moved.
            s["worker_demod_ms_per_frame"] = round(
                1000 * self.pool.tm.get("worker_demod", 0.0) / nf, 2)
            s["worker_fec_ms_per_frame"] = round(
                1000 * self.pool.tm.get("worker_fec", 0.0) / nf, 2)
            s["pool_demod_ms_per_frame"] = {
                k: round(1000 * v / nf, 2) for k, v in self.pool.tm.items()
                if not k.startswith("worker_")}
        return s


# ---------------------------------------------------------------------------
# batched BICM + LDPC for the 64K codes (the same shape m9_fast uses at 16K)
# ---------------------------------------------------------------------------

def _point_metric(zr, zi, pr, pi, p2):
    """|p|^2 - 2 Re(z conj(p)) for every point, WITHOUT calling BLAS.

    m9_fast writes this as one float32 sgemm, and on RF33 that is right.  Here
    it deadlocked the chain: the LDM path drives its own thread pool, and a
    multi-threaded BLAS `X @ P` called from inside eight of those workers
    wedged every one of them (stack dumps: eight threads parked in the same
    `job`, no CPU, no progress -- twice, once as a hang and once as an access
    violation).  Nested parallelism is the hazard, and pinning BLAS through
    the environment only fixes the copy of the process that remembered to.

    A Core-Layer alphabet is QPSK -- FOUR points.  A (n,2) x (2,4) product is
    not work a gemm should ever have been asked to do; broadcasting is both
    faster here and structurally incapable of nesting a thread pool inside a
    thread pool.  The arithmetic is identical term for term.
    """
    q = zr[:, None] * pr[None, :]
    q += zi[:, None] * pi[None, :]
    q *= np.float32(-2.0)
    q += p2[None, :]
    return q


def demap_batch_qam(cells2d, ch, ex=None, threads=4, chunk=4):
    """(nblk, ncell) -> (nblk, ncell*nbits) max-log LLRs, float32.

    Identical to calling `m6_bicm.demap_llr` per row: sigma^2 is a mean over
    that row's own cells and everything else is elementwise.  The distance is
    evaluated as |p|^2 - 2 Re(z conj(p)) with the constant |z|^2 folded back
    only into sigma^2 -- m9_fast's re-association, which is exact for the
    bitwise-minima DIFFERENCE.
    """
    pts = B.points_for(ch.mod, ch.rate)
    nb = ch.mod_bits
    nblk, ncell = cells2d.shape
    out = np.empty((nblk, ncell, nb), np.float32)
    pr = np.asarray(pts.real, np.float32)
    pi = np.asarray(pts.imag, np.float32)
    p2f = (pts.real ** 2 + pts.imag ** 2).astype(np.float32)
    sub = tuple(range(2, 2 + nb - 1))

    def job(lo, hi):
        z = cells2d[lo:hi]
        k = hi - lo
        zr = np.asarray(z.real, np.float32)
        zi = np.asarray(z.imag, np.float32)
        q = _point_metric(zr.ravel(), zi.ravel(), pr, pi, p2f)
        q = q.reshape(k, ncell, len(pts))
        z2 = zr * zr + zi * zi
        s2 = np.maximum((z2 + q.min(2)).mean(1), 1e-9)[:, None]
        qr = q.reshape(k, ncell, *([2] * nb))
        head = (slice(None), slice(None))
        for i in range(nb):
            a = qr[head + (slice(None),) * i + (1,)].min(axis=sub)
            b = qr[head + (slice(None),) * i + (0,)].min(axis=sub)
            out[lo:hi, :, i] = (a - b) / s2
    nj = max(1, (nblk + chunk - 1) // chunk)
    if ex is None:
        for i in range(nj):
            job(i * chunk, min(nblk, (i + 1) * chunk))
    else:
        list(ex.map(lambda i: job(i * chunk, min(nblk, (i + 1) * chunk)),
                    range(nj)))
    return out.reshape(nblk, ncell * nb)


def ldpc_batch(lam2d, ch, gpu=None, iters=50, ex=None, threads=4, alpha=None):
    """`PlpChain.decode_llr`'s alpha-ladder semantics, batched.

    The ladder matters MORE here than on RF33: m6_bicm records that alpha 0.75
    is right for the 16200 codes and WRONG for the 64800 rate-6/15 code, where
    it parks in a trapping set at 4-16 unsatisfied of 38880 and stays there --
    "exactly like a link 0.2 dB short".  First converged rung wins; a rung is
    only tried on the rows that are still failing, so on healthy air the
    ladder costs nothing.
    """
    ladder = (alpha,) if alpha else tuple(ch.ALPHA_LADDER)
    nb = int(lam2d.shape[0])
    bits = np.empty((nb, ch.ninner), np.uint8)
    conv = np.zeros(nb, bool)
    bad = np.full(nb, 1 << 30, np.int32)
    todo = np.arange(nb)
    for a in ladder:
        if not len(todo):
            break
        sub = lam2d if len(todo) == nb else lam2d[todo]
        if gpu is not None:
            b, c, _it, d = gpu.decode(sub, alpha=a)
        else:
            b, c, _it, d = _ldpc_cpu(sub, ch, a, iters, ex, threads)
        take = c | (d < bad[todo])
        sel = todo[take]
        bits[sel] = b[take]
        bad[sel] = d[take]
        conv[sel] = c[take]
        todo = todo[~c]
    if len(todo) and SP_RESCUE:
        # E86: the E58/E60 margin work landed in m9_fast and NEVER reached
        # this path, so the LDM carrier -- the WEAKER of the two -- has been
        # decoding against the old threshold while RF33 enjoys the improved
        # one.  This is the cleanest of the three levers to port: exact
        # sum-product on the rows min-sum could not close, which is where a
        # trapping set (not a lack of information) is what stopped it.
        #
        # Only NEAR MISSES are offered.  A block 4000 checks short is not
        # short of arithmetic, it is short of signal, and spending float64
        # BP on it buys nothing and costs a live deadline.
        near = todo[bad[todo] <= SP_NEAR_MAX]
        if len(near) > SP_MAX_BLOCKS:
            # E60's shape: take the CLOSEST misses, not the first ones, so a
            # capped batch spends its budget where it is most likely to pay.
            near = near[np.argsort(bad[near], kind="stable")[:SP_MAX_BLOCKS]]
        if len(near):
            dl = (None if SP_BUDGET_S is None
                  else time.perf_counter() + SP_BUDGET_S)
            b, c, _it, d = M16.sum_product_decode_batch(
                lam2d[near], ch.checks, ch.ninner, iters=SP_ITERS,
                deadline=dl)
            got = np.flatnonzero(c)
            if len(got):
                sel = near[got]
                bits[sel] = b[got]
                conv[sel] = True
                bad[sel] = 0
    return bits, conv, bad


def _ldpc_cpu(lam2d, ch, alpha, iters, ex=None, threads=4):
    nb = len(lam2d)
    if nb == 0:
        return (np.empty((0, ch.ninner), np.uint8), np.zeros(0, bool),
                np.zeros(0, np.int32), np.zeros(0, np.int32))
    nch = 1 if ex is None else max(1, min(nb, threads))
    step = (nb + nch - 1) // nch
    sls = [slice(i * step, min(nb, (i + 1) * step))
           for i in range(nch) if i * step < nb]

    def one(sl):
        return LD.min_sum_decode_batch(lam2d[sl], ch.checks, ch.ninner,
                                       iters=iters, alpha=alpha,
                                       dtype=np.float32)
    res = ([one(s) for s in sls] if ex is None else list(ex.map(one, sls)))
    return (np.concatenate([r[0] for r in res]),
            np.concatenate([r[1] for r in res]),
            np.concatenate([r[2] for r in res]),
            np.concatenate([r[3] for r in res]))
