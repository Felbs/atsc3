#!/usr/bin/env python3
"""M13 -- Subframe 1 / PLP 1 on RF33: the multiplex's OTHER television service.

Everything decoded so far on RF33 lives in Subframe 0.  The Frame has two:

    bootstrap  13824
    Preamble    1 x (8192 + 1536)      <- L1-Basic + L1-Detail live here
    Subframe 0 35 x (8192 + 1536)      <- PLP 0 (64QAM 11/15) + PLP 16
    Subframe 1 75 x (16384 + 1536)     <- PLP 1, and nothing has ever read it

and Subframe 1 is the larger half of the Frame by a factor of four.  M6 and M7
both named it as the next service and both stopped at the same recorded
blocker: "A/327 Section 6 specifies a frequency de-interleaver RESET rule for
the second Subframe which this project has not implemented; missing it would
break PLP 1 specifically."

THE PARAMETERS ARE NOT GUESSED -- THEY ARE READ
-----------------------------------------------
Straight out of the decoded L1-Detail, scoped to `i=1/`:

    L1D_fft_size 1 -> 16K      L1D_guard_interval 6 -> 1536
    L1D_num_ofdm_symbols 74 -> 75 symbols     pattern 3      Cred 0
    sbs_first 1  sbs_last 1    L1D_sbs_null_cells 1609
    L1D_frequency_interleaver 1  (ON)
    PLP 1: start 0, size 947700, fec_type 1 -> 64800, mod 3 -> 256QAM,
           cod 9 -> 11/15, TI_mode 2 -> HTI,
           HTI_num_ti_blocks 2, HTI_num_fec_blocks 116, cell_interleaver 0

and the closed arithmetic agrees before a single cell is demodulated:
947700 cells / (64800/8 bits) = **117 FEC Blocks exactly**, 117 = 3 x 39.

THE +1 CONVENTION IS TAKEN FROM A WORKING DECODER, NOT FROM THE SPEC TEXT
-------------------------------------------------------------------------
`L1D_plp_HTI_num_fec_blocks` could plausibly be per-TI-Block or per
interleaving frame, and the two readings differ by 3x.  PLP 0 settles it
without an argument: its L1 says num_ti_blocks 1 / num_fec_blocks 73, and the
decoder that has produced watchable video for two days uses nti=2, n_fec=74,
ncols = n_fec_max // nti = 37.  So the field is (value + 1) and it counts the
WHOLE interleaving frame.  Applied to PLP 1: nti=3, n_fec=117, ncols=39 --
which is exactly the "3 TI Blocks x 39 columns" M6 measured off the waveform.

THE RESET RULE IS A HYPOTHESIS WITH AN ORACLE, NOT A BLOCKER
-------------------------------------------------------------
`cell_pool_g` already exposes `fi_offset`, the frequency interleaver's symbol
counter origin.  So "does the counter reset at the Subframe boundary?" is a
one-line sweep with a decisive referee: **LDPC convergence and BCH syndrome
zero over 117 FEC Blocks**.  A wrong interleaver is not subtly worse, it is
0/117.  Note also that Subframe 1 runs a DIFFERENT FFT size, so it necessarily
uses a different interleaving sequence -- which is a reason to suspect a reset,
and not a reason to assume one.

Usage:
    python m13_sf1.py long_rf33.cs16 --rate 8e6 --l1 m4_l1detail_hit_rf33.json
    python m13_sf1.py ... --sweep-fi          # the reset-rule experiment
    python m13_sf1.py ... --frames 8 --bb m13_plp1.bb
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

import m2_pilots as MP                                           # noqa: E402
import m3_spec as S                                              # noqa: E402
import m4_scrambler as SC                                        # noqa: E402
import m6_bicm as B                                              # noqa: E402
import m6_cells as C                                             # noqa: E402
import m6_payload as P6                                          # noqa: E402
import m6_tbi as TBI                                             # noqa: E402
import m7_route as R7                                            # noqa: E402
import m16_demap_fast as DF                                      # noqa: E402
import spec_l1syntax as LS                                       # noqa: E402
import spec_pilots as P                                          # noqa: E402
from m2_pilots import load                                       # noqa: E402
from m10_core import fine_t0, frame_samples                      # noqa: E402

FS_POST = 6.912e6
_FFT = {0: 8192, 1: 16384, 2: 32768}
_MOD = {0: "QPSK", 1: "16QAM", 2: "64QAM", 3: "256QAM",
        4: "1024QAM", 5: "4096QAM"}
_BITS = {"QPSK": 2, "16QAM": 4, "64QAM": 6, "256QAM": 8,
         "1024QAM": 10, "4096QAM": 12}
_RATE = {i: "%d/15" % (i + 2) for i in range(12)}


# ---------------------------------------------------------------------------
# L1 -> a Subframe's geometry
# ---------------------------------------------------------------------------

def scope(l1d, i, j=None):
    """The L1-Detail fields belonging to Subframe i (and PLP j within it).

    `from_l1` flattens the whole field list by NAME, which silently lets
    Subframe 1 overwrite Subframe 0 and puts every PLP of every Subframe into
    one list.  That is harmless while only Subframe 0 is decoded and wrong the
    moment it is not, so this reads the loop path instead.
    """
    pre = f"i={i}/" + (f"j={j}/" if j is not None else "")
    out = {}
    for e in l1d:
        p = e["path"] if isinstance(e, dict) else e[0]
        n = e["name"] if isinstance(e, dict) else e[1]
        v = e["value"] if isinstance(e, dict) else e[2]
        if p == pre:
            out[n] = v
    return out


def plp_record(f):
    mod = _MOD[f["L1D_plp_mod"]]
    ninner = 64800 if f["L1D_plp_fec_type"] else 16200
    nti = f.get("L1D_plp_HTI_num_ti_blocks", 0) + 1
    n_fec = f.get("L1D_plp_HTI_num_fec_blocks", 0) + 1
    n_fec_max = f.get("L1D_plp_HTI_num_fec_blocks_max", n_fec - 1) + 1
    return dict(
        id=f["L1D_plp_id"], layer=f.get("L1D_plp_layer", 0),
        lls=bool(f.get("L1D_plp_lls_flag", 0)),
        start=f["L1D_plp_start"], size=f["L1D_plp_size"],
        mod=mod, rate=_RATE[f["L1D_plp_cod"]], ninner=ninner,
        cells_per_fec=ninner // _BITS[mod],
        ti_mode=f.get("L1D_plp_TI_mode"), nti=nti, n_fec=n_fec,
        n_fec_max=n_fec_max,
        cell_interleaver=f.get("L1D_plp_HTI_cell_interleaver", 0),
        inter_subframe=f.get("L1D_plp_HTI_inter_subframe", 0))


def geometry_for_subframe(l1b, l1d, i, g0):
    """Geometry for Subframe i>0.  `g0` supplies the Preamble parameters.

    np_sym is ZERO here and that is the substantive part: the Preamble belongs
    to the Frame, not to this Subframe, so this Subframe contributes no
    Preamble spare and carries no L1 cells.  Its t0 is its own first sample.
    """
    f = scope(l1d, i)
    nplp = f["L1D_num_plp"] + 1
    plps = [plp_record(scope(l1d, i, j)) for j in range(nplp)]
    return C.Geometry(
        nfft=_FFT[f["L1D_fft_size"]],
        gi=LS.GUARD_INTERVAL_SAMPLES[f["L1D_guard_interval"]],
        cred=f["L1D_reduced_carriers"],
        pattern=P.PATTERNS[f["L1D_scattered_pilot_pattern"]],
        nsym=f["L1D_num_ofdm_symbols"] + 1,
        sbs_first=f["L1D_sbs_first"], sbs_last=f["L1D_sbs_last"],
        np_sym=0,
        pre_nfft=g0.pre_nfft, pre_gi=g0.pre_gi, pre_dx=g0.pre_dx,
        pre_cred=g0.pre_cred,
        l1b_cells=0, l1d_cells=0,
        freq_interleaver=f.get("L1D_frequency_interleaver", 1),
        sbs_null_signalled=f.get("L1D_sbs_null_cells"),
        plps=plps, label=f"subframe {i}")


def subframe_offset(g0, i):
    """Samples from the Frame's t0 (the Preamble's first sample) to Subframe i.

    Subframe 0 starts after the Preamble; Subframe 1 after Subframe 0.  The
    same expression `_demod_data_g` uses, evaluated at l == nsym.
    """
    off = (g0.pre_nfft + g0.pre_gi) * g0.np_sym
    if i > 0:
        off += (g0.nfft + g0.gi) * g0.nsym
    return off


# ---------------------------------------------------------------------------
# the decode
# ---------------------------------------------------------------------------

def hti_blocks(cells, plp):
    """HTI de-interleave -> the FEC Blocks, in order.

    Identical in shape to what m9_fast does for PLP 0, driven by this PLP's own
    L1: `ncols = n_fec_max // nti` columns of `cells_per_fec` rows per TI
    Block, read out by the twisted read order.
    """
    nrows = plp["cells_per_fec"]
    ncols = plp["n_fec_max"] // plp["nti"]
    per_ti = len(cells) // plp["nti"]
    out = []
    for ti in range(plp["nti"]):
        seg = cells[ti * per_ti:(ti + 1) * per_ti]
        for j in range(ncols):
            out.append(TBI.fec_block(seg, j, nrows, ncols, 0))
    return np.array(out), nrows, ncols


def decode_subframe(y, t0_sf, g, plp, fi_offset, iters=50, verbose=True,
                    max_blocks=None, ex=None):
    rep = {}
    pool, info = C.cell_pool_g(y, t0_sf, g, report=rep, fi_offset=fi_offset)
    pool_pred = g.pool_size(0)
    dummy_n = pool_pred - plp["size"]
    if verbose:
        print(f"    pool {len(pool)} (predicted {pool_pred}), "
              f"PLP {plp['id']} size {plp['size']}, dummy {dummy_n}, "
              f"data coherence {rep['data_coherence_mean']:.4f}")
    if len(pool) < plp["start"] + plp["size"]:
        return dict(ok=False, why="pool shorter than the PLP", pool=len(pool),
                    pool_pred=pool_pred, dummy=dummy_n)

    # per-symbol residual gain/phase against this PLP's own alphabet
    alph = [B.points_for(plp["mod"], plp["rate"])]
    reg = np.zeros(len(pool), np.int8)
    dv = np.zeros(len(pool), complex)
    if dummy_n > 0:
        start = plp["start"] + plp["size"]
        reg[start:] = 2
        dv[start:] = 1.0 - 2.0 * SC.sequence(len(pool))[start:]
    pool, _g = C.cpe_correct(pool, info["symbol_of"], alph, reg, dv)

    cells = pool[plp["start"]:plp["start"] + plp["size"]]
    blocks, nrows, ncols = hti_blocks(cells, plp)
    n = len(blocks) if max_blocks is None else min(len(blocks), max_blocks)
    ch = B.PlpChain(plp["ninner"], plp["mod"], plp["rate"], iters=iters)
    stream, bounds = bytearray(), []
    conv = bch = 0
    unsat = []
    t = time.time()
    # M16's demapper, threaded ACROSS FEC Blocks rather than within one.  Each
    # Block keeps its OWN sigma2 -- `demap_llr` estimates it per call, and
    # demapping the Frame as one array would silently make that estimate global,
    # which is a numerical change dressed up as an optimisation.
    lams = None
    if ex is not None:
        def one(k):
            q = DF.demap_llr_fast(blocks[k], plp["mod"], plp["rate"],
                                  dtype=np.float32)
            lam = np.empty(ch.ninner)
            lam[ch.lam_of_q] = q
            return lam
        lams = list(ex.map(one, range(n)))
    for m in range(n):
        r = (ch.decode_llr(lams[m]) if lams is not None
             else ch.decode(blocks[m]))
        if not r["converged"]:
            unsat.append(int(r["unsatisfied"]))
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
    dt = time.time() - t
    if verbose:
        print(f"    HTI {plp['nti']} TI Blocks x {ncols} cols x {nrows} rows"
              f"  ->  {len(blocks)} FEC Blocks")
        print(f"    FEC {conv}/{n} LDPC converged, {bch}/{n} BCH zero"
              f"   ({dt:.1f}s)"
              + (f"   median unsatisfied {int(np.median(unsat))}"
                 if unsat else ""))
    return dict(ok=True, converged=conv, bch_ok=bch, n_blocks=n,
                stream=bytes(stream), bounds=bounds, pool=len(pool),
                pool_pred=pool_pred, dummy=dummy_n,
                coherence=float(rep["data_coherence_mean"]),
                unsat_median=(int(np.median(unsat)) if unsat else 0))


def main(argv=None):
    ap = argparse.ArgumentParser(prog="m13 subframe 1")
    ap.add_argument("capture")
    ap.add_argument("--rate", type=float, default=8e6)
    ap.add_argument("--fmt", default=None)
    ap.add_argument("--l1", default="m4_l1detail_hit_rf33.json")
    ap.add_argument("--subframe", type=int, default=1)
    ap.add_argument("--plp", type=int, default=None)
    ap.add_argument("--frames", type=int, default=1)
    ap.add_argument("--start", type=float, default=0.0)
    ap.add_argument("--fi-offset", type=int, default=None)
    ap.add_argument("--sweep-fi", action="store_true",
                    help="the reset-rule experiment: try every plausible "
                         "symbol-counter origin and let LDPC/BCH referee")
    ap.add_argument("--iters", type=int, default=50)
    ap.add_argument("--threads", type=int, default=12,
                    help="demap threads; 0 uses the original per-Block path")
    ap.add_argument("--max-blocks", type=int, default=None)
    ap.add_argument("--bb", default=None)
    ap.add_argument("--json", default=None)
    a = ap.parse_args(argv)
    path = a.capture if os.path.isabs(a.capture) else os.path.join(HERE,
                                                                   a.capture)
    l1p = a.l1 if os.path.isabs(a.l1) else os.path.join(HERE, a.l1)
    js = json.load(open(l1p))
    l1b, l1d = js["l1_basic"], js["l1_detail"]
    fmt = a.fmt or MP.guess_format(path)

    g0 = C.Geometry.rf33_legacy()
    g = geometry_for_subframe(l1b, l1d, a.subframe, g0)
    plps = {p["id"]: p for p in g.plps}
    plp = plps[a.plp if a.plp is not None else sorted(plps)[0]]

    from concurrent.futures import ThreadPoolExecutor
    EX = ThreadPoolExecutor(max_workers=a.threads) if a.threads else None
    print(f"M13 -- Subframe {a.subframe} of {os.path.basename(path)}")
    print("=" * 72)
    print(f"  geometry from L1: FFT {g.nfft//1024}K  GI {g.gi}  Cred {g.cred}"
          f"  {g.pattern}  {g.nsym} symbols  SBS {int(g.sbs_first)}/"
          f"{int(g.sbs_last)}  FI {'ON' if g.freq_interleaver else 'OFF'}")
    print(f"    null cells computed {g.n_null} vs signalled "
          f"{g.sbs_null_signalled}  "
          f"{'PASS' if g.n_null == g.sbs_null_signalled else '*** FAIL ***'}")
    print(f"  PLP {plp['id']}: {plp['mod']} {plp['rate']}  Ninner "
          f"{plp['ninner']}  {plp['cells_per_fec']} cells/Block  "
          f"size {plp['size']} = {plp['size']//plp['cells_per_fec']} Blocks")
    print(f"    HTI: {plp['nti']} TI Blocks, n_fec {plp['n_fec']} "
          f"(max {plp['n_fec_max']}), cell interleaver "
          f"{plp['cell_interleaver']}, inter-subframe {plp['inter_subframe']}")

    # frame_samples(g0) is Subframe 0's Frame -- bootstrap + Preamble + 35
    # symbols -- and Subframe 1 is FOUR TIMES that on its own.  Loading the
    # short one reads the FFT window off the end of the buffer, which is how
    # this was caught.
    total = frame_samples(g0) + (g.nfft + g.gi) * g.nsym
    per = total / FS_POST
    span = per + 0.06
    off = subframe_offset(g0, a.subframe)
    print(f"  Frame {total} samples ({per*1000:.3f} ms); Subframe "
          f"{a.subframe} starts +{off}, runs "
          f"{(g.nfft + g.gi) * g.nsym} samples")

    results = []
    for k in range(a.frames):
        y, _fs, _cfo, _h = load(path, a.rate, fmt=fmt, span_sec=span,
                                start_sec=a.start + k * per)
        t0, coh = fine_t0(y, g0)
        print(f"\n  Frame {k}: t0 {t0}, preamble coherence {coh:.4f}")
        cands = ([a.fi_offset] if a.fi_offset is not None else
                 ([0, g0.np_sym, g0.np_sym + g0.nsym, g0.nsym,
                   g0.np_sym + g0.nsym + 1] if a.sweep_fi else [0]))
        for fi in cands:
            if a.sweep_fi:
                print(f"    -- fi_offset {fi} "
                      f"({'RESET at the Subframe' if fi == 0 else 'continued'})")
            r = decode_subframe(y, t0 + off, g, plp, fi, iters=a.iters,
                                max_blocks=a.max_blocks, ex=EX)
            r.update(frame=k, fi_offset=fi)
            results.append(r)
            if a.sweep_fi and r.get("ok") and r["bch_ok"]:
                print(f"    *** fi_offset {fi} DECODES -- "
                      f"{r['bch_ok']}/{r['n_blocks']} BCH zero ***")
                break

    good = [r for r in results if r.get("ok") and r.get("bch_ok")]
    print("\n" + "=" * 72)
    if good:
        tot = sum(r["bch_ok"] for r in good)
        blk = sum(r["n_blocks"] for r in good)
        print(f"  PLP {plp['id']} DECODES: {tot}/{blk} FEC Blocks BCH zero "
              f"across {len(good)} Frame(s), fi_offset "
              f"{sorted({r['fi_offset'] for r in good})}")
    else:
        print(f"  PLP {plp['id']} did NOT decode -- "
              f"0 BCH-clean Blocks in {len(results)} attempt(s)")
        for r in results[:6]:
            print(f"    fi_offset {r.get('fi_offset')}: "
                  f"{r.get('why', '')} converged {r.get('converged')}, "
                  f"median unsatisfied {r.get('unsat_median')}")

    if a.bb and good:
        # ONE .bb record for the whole run, with the bounds re-based, because
        # the ALP walk is a single stream: writing one record per Frame would
        # invite a walker restart at every Frame boundary, which is exactly the
        # bug M7 paid for once already.
        stream, bounds = bytearray(), []
        for r in good:
            bounds += [b + len(stream) for b in r["bounds"]]
            stream += r["stream"]
        out = a.bb if os.path.isabs(a.bb) else os.path.join(HERE, a.bb)
        with open(out, "wb") as fh:
            # the per-Frame records m7_route's reporting expects
            fr = [dict(frame=r["frame"], conv=r["converged"],
                       bch=r["bch_ok"], n=r["n_blocks"], coh=r["coherence"])
                  for r in good]
            R7.bb_write(fh, a.start, fr, bytes(stream), bounds)
        print(f"  wrote {a.bb}: {len(stream)} bytes, {len(bounds)} anchors")
    if a.json:
        with open(a.json if os.path.isabs(a.json)
                  else os.path.join(HERE, a.json), "w") as fh:
            json.dump([{k: v for k, v in r.items() if k != "stream"}
                       for r in results], fh, indent=1, default=str)
    return 0 if good else 1


if __name__ == "__main__":
    sys.exit(main())
