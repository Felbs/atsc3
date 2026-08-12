#!/usr/bin/env python3
"""M8 Step 1 -- L1 off ANY multiplex: multi-symbol preambles and Mode-1 repetition.

M4/M6/M7 decoded one multiplex end to end.  Two of its choices turned out to be
that multiplex's, not the standard's, and both are load-bearing:

  1. **NP = 1.**  That multiplex signals `L1B_preamble_num_symbols = 0`, so its
     whole Preamble is one OFDM symbol and A/322 7.2.5.2's L1-Detail Block
     Interleaver is Lc = 1 column, i.e. the identity.  m4_l1detail.py refuses
     anything else.  Two of the other multiplexes here use **NP = 2**, because
     their L1-Basic is Mode 1, which costs 3820 of the first Preamble symbol's
     4851 cells and leaves nowhere near enough room for L1-Detail.

  2. **L1-Detail Mode 3.**  A/322 6.5.2.7 applies parity REPETITION to
     *L1-Basic Mode 1 and L1-Detail Mode 1 only*.  Mode 3 therefore never
     exercised it, and `m4_l1detail.geometry()` has no Nrepeat term at all.

WHY THE SECOND ONE WAS NOT A GUESS
----------------------------------
It was signalled.  `L1B_L1_Detail_total_cells` is a field the transmitter
publishes, and the receiver can compute the same number independently:

    Nrepeat = 2*floor(C*Nouter) + D                       (6.5.2.7 Step 1)
    cells   = (Nfec + Nrepeat) / eta                      (6.5.2.8 + 6.5.2.10)

Without the repetition term the two numbers disagree -- 1944 computed against
3611 signalled on RF30, 2232 against 4387 on RF25 -- and that disagreement was
sitting in plain sight in M4's own output, printed as "MISMATCH", for two
different channels, before any of this was understood.  With it, both match
EXACTLY.  **A signalled field that the receiver can also derive is a free gate,
and this one had been failing silently.**

The LDPC agreed second, not first: 98 unsatisfied checks of 12960 -- close
enough to look like a link problem and not close enough to decode.  Reading
"nearly converged" as noise rather than as a missing term would have sent this
back to the antenna.

THE GATES (each with controls that must fail)
---------------------------------------------
1. Computed `cells` == signalled `L1B_L1_Detail_total_cells`, on every channel.
2. LDPC convergence + zero BCH syndrome + L1D_crc.
3. The Lc = NP Block Interleaver's column height Lr is SWEPT, not assumed:
   `--sweep-lr` scores Lr over a +-14 window.  Only floor(total/Lc) works, and
   every neighbour fails, so the "first Lr*Lc cells interleaved, remainder
   passed through" reading of 7.2.5.2 is settled by the air.
4. The QPSK cliff.  L1-Detail ends part-way through the LAST Preamble symbol
   and the cells after it are data-PLP cells of a much higher order.  Where
   that MER cliff falls is predicted by the mapping and measured independently
   of it (`--cliff`).

Usage:
    python m8_l1.py long_rf30.cs16 --rate 8e6 --json m8_l1_rf30.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m3_spec as S                                             # noqa: E402
import m3_freqint as FI                                         # noqa: E402
import m3_ldpc as LD                                            # noqa: E402
import m3_l1basic as L1                                         # noqa: E402
import m3_crc as CRC                                            # noqa: E402
import m4_l1detail as D4                                        # noqa: E402
import m4_scrambler as SC                                       # noqa: E402
import spec_bicm as BI                                          # noqa: E402
import spec_l1syntax as LS                                      # noqa: E402
from m2_pilots import load                                      # noqa: E402
from m3_preamble import (analyse, fft_bins, preamble_geometry,   # noqa: E402
                         mer_db, QPSK)

SPAN_SEC = 0.10


# ---------------------------------------------------------------------------
# A/322 7.2.5.1 -- Preamble symbols after the first
# ---------------------------------------------------------------------------

def demod_preamble_symbol(y, tw, nfft, gi, dx, cred, shift=0):
    """Equalised data cells of ONE Preamble symbol, in CARRIER order.

    Identical to the first symbol's treatment except for NoC: 7.2.5.1 fixes the
    FIRST Preamble symbol at the minimum NoC for its FFT size and lets L1-Basic
    signal the rest via L1B_preamble_reduced_carriers.  The pilots are the same
    (every DX-th carrier, DY = 1) and their VALUES come from A/322 8.1.2, so a
    wrong NoC collapses the pilot coherence rather than passing quietly.
    """
    lo, n, pilot, cp, data, _ = preamble_geometry(nfft, gi, dx, cred)
    Y = np.fft.fftshift(np.fft.fft(y[tw:tw + nfft]))
    amp = S.PREAMBLE_PILOT_BOOST[(nfft, gi)][1]
    ref = np.array(S.pilot_values(n, amp))
    hp = Y[fft_bins(nfft, lo + pilot + shift)] / ref[pilot]
    H = (np.interp(np.arange(n), pilot, hp.real).astype(np.complex128)
         + 1j * np.interp(np.arange(n), pilot, hp.imag))
    z = Y[fft_bins(nfft, lo + data + shift)] / H[data]
    coh = float(abs(np.sum(hp[:-1] * np.conj(hp[1:])))
                / max(np.sum(np.abs(hp) ** 2), 1e-30))
    return z, coh, len(data)


def l1d_deinterleave(yv, total, lc, lr=None):
    """Inverse of the A/322 7.2.5.2 L1-Detail Block Interleaver.

    Forward: the first Lr*Lc cells are written row-wise into an Lr x Lc array
    and read out column by column; cells past Lr*Lc pass through unchanged.
    Lr = floor(total/Lc) -- swept and gated, see --sweep-lr.
    """
    if lc == 1:
        return np.array(yv)
    lr = total // lc if lr is None else lr
    x = np.array(yv, dtype=complex)
    for i in range(lc):
        for j in range(lr):
            m, nn = lc * j + i, lr * i + j
            if m < total and nn < total:
                x[m] = yv[nn]
    return x


# ---------------------------------------------------------------------------
# A/322 6.5.2.7 + 6.5.2.8 -- L1-Detail FEC geometry WITH repetition
# ---------------------------------------------------------------------------

def geometry(size_bytes, fec_type):
    """Every term straight from A/322's own equations; nothing fitted.

    The only difference from m4_l1detail.geometry() is Nrepeat, which is zero
    for every mode except 1 -- which is why one multiplex could not reveal it.
    """
    mode = fec_type + 1
    kind = "L1-Detail-%d" % mode
    p = BI.PUNCTURING[kind]
    kldpc, eta, nlp = p["Kldpc"], p["eta"], p["Nldpc_parity"]
    kseg = BI.L1D_KSEG[mode]
    ksig = size_bytes * 8
    if ksig > kseg:
        raise SystemExit(f"Ksig {ksig} > Kseg {kseg}: L1-Detail is segmented "
                         f"into multiple FEC Frames -- not implemented")
    nouter = ksig + 168
    r = BI.REPETITION.get(kind)
    nrep = 0 if r is None else 2 * int(r["C"] * nouter) + r["D"]
    npunc_t = int(p["A"] * (kldpc - nouter)) + p["B"]
    nfec_t = nouter + nlp - npunc_t
    nfec = (nfec_t // eta) * eta
    npunc = npunc_t - (nfec - nfec_t)
    ntx = nfec + nrep
    return dict(mode=mode, kind=kind, Ksig=ksig, Nouter=nouter, Kldpc=kldpc,
                Ninner=16200, Nldpc_parity=nlp, eta=eta, Nrepeat=nrep,
                Nfec=nfec, Npunc=npunc, n_parity_kept=nlp - npunc,
                Ntx=ntx, cells=ntx // eta,
                rate={1: "3/15", 2: "3/15"}.get(mode, "6/15"))


def decode(cells, g, iters=150, no_repeat=False, _c={}):
    """A/322 6.5.2.9's transmitted word is

        [ Nouter info+BCH ][ Nrepeat repeated parity ][ Nldpc_parity - Npunc ]

    and the repeated block is the FIRST Nrepeat bits of the *same* permuted
    parity the tail carries -- so their LLRs ADD onto the same codeword
    positions.  `no_repeat=True` reproduces the M4 reading as a control.
    """
    bl = L1.cells_to_llr(cells[:g["cells"]], nfec=g["Ntx"], eta=g["eta"])
    llr = np.zeros(g["Ninner"])
    ipos, padded = L1.info_positions(
        g["Kldpc"], g["Nouter"],
        pattern=list(BI.SHORTENING_PATTERN[g["kind"]]))
    llr[:g["Kldpc"]][padded] = L1.BIG
    llr[ipos] = bl[:g["Nouter"]]
    allp = L1.parity_positions(g["Nldpc_parity"], g["Kldpc"], g["Ninner"],
                               groupwise=list(BI.GROUPWISE[g["kind"]]))
    o = g["Nouter"]
    if g["Nrepeat"] and not no_repeat:
        np.add.at(llr, allp[:g["Nrepeat"]], bl[o:o + g["Nrepeat"]])
        o += g["Nrepeat"]
    tail = bl[o:o + g["n_parity_kept"]]
    np.add.at(llr, allp[:len(tail)], tail)
    key = (g["rate"], g["Ninner"])
    if key not in _c:
        _c[key] = LD.parity_check(g["rate"], n=g["Ninner"])
    bits, conv, it, bad = LD.min_sum_decode(llr, _c[key], g["Ninner"],
                                            iters=iters)
    nouter = bits[ipos]
    syn = L1.bch_syndrome(nouter)
    return dict(converged=bool(conv), iters=int(it), unsatisfied=int(bad),
                bch_ok=syn == 0, bch_syndrome=int(syn), nouter=nouter)


# ---------------------------------------------------------------------------
# the whole Preamble, however many symbols it has
# ---------------------------------------------------------------------------

def preamble_cells(path, rate, fmt, start_sec, quiet=True):
    """Return (l1b_fields, l1d_stream, ctx) for one Frame.

    `l1d_stream` is the L1-Detail cell region in TRANSMISSION order: the cells
    of the first Preamble symbol after L1-Basic, then every cell of Preamble
    symbols 1..NP-1, truncated to L1B_L1_Detail_total_cells.
    """
    rep = {}
    _r, z0, _Y, geo = analyse(path, rate, fmt, report=rep,
                              start_sec=start_sec, quiet=quiet)
    (lo, n, pilot, cp, data, shift, cred, nfft, gi, dx, mode) = geo
    x0 = FI.deinterleave(z0, nfft, 0, direction="forward", toggle="i")
    f = D4.l1basic_of(x0, mode)
    if f is None:
        return None, None, None
    npre = f["L1B_preamble_num_symbols"] + 1
    nb = S.L1_BASIC_CELLS_PRINTED[mode]
    tot = f["L1B_L1_Detail_total_cells"]
    parts, cohs, sym_cells = [x0[nb:]], [], [len(x0)]
    ctx = dict(nfft=nfft, gi=gi, dx=dx, cred0=cred, shift=shift,
               l1b_mode=mode, NP=npre, l1b_cells=nb, total_cells=tot,
               t0=rep["fft_window_start"], start_sec=start_sec)
    if npre > 1:
        y, _fs, _cfo, _hit = load(path, rate, fmt=fmt, span_sec=SPAN_SEC,
                                  start_sec=start_sec)
        credp = f["L1B_preamble_reduced_carriers"]
        for s in range(1, npre):
            zs, coh, nd = demod_preamble_symbol(
                y, ctx["t0"] + (nfft + gi) * s, nfft, gi, dx, credp, shift)
            parts.append(FI.deinterleave(zs, nfft, s, direction="forward",
                                         toggle="i"))
            cohs.append(coh)
            sym_cells.append(nd)
    ctx["extra_symbol_coherence"] = cohs
    ctx["symbol_cells"] = sym_cells
    stream = np.concatenate(parts)
    ctx["available"] = int(len(stream))
    ctx["last_symbol_spare"] = int(len(stream) - tot)
    return f, stream[:tot], ctx


def solve(path, rate, fmt, start_sec, quiet=True):
    f, yv, ctx = preamble_cells(path, rate, fmt, start_sec, quiet)
    if f is None:
        return None
    g = geometry(f["L1B_L1_Detail_size_bytes"], f["L1B_L1_Detail_fec_type"])
    xd = l1d_deinterleave(yv, ctx["total_cells"], ctx["NP"])
    r = decode(xd, g)
    out = dict(l1b=f, geom=g, ctx=ctx, fec=r, l1d_stream=yv, l1d_cells=xd)
    if not (r["converged"] and r["bch_ok"]):
        out["ok"] = False
        return out
    d = SC.descramble(np.array(r["nouter"][:g["Ksig"]], np.uint8))
    out["crc_ok"] = bool(
        CRC.crc32(d[:g["Ksig"] - 32], init=0xFFFFFFFF)
        == int("".join(str(int(c)) for c in d[g["Ksig"] - 32:]), 2))
    out["bits"] = d
    best = None
    for mm in (1, 0):
        try:
            p = D4.parse_l1d(d, f, mm)
        except Exception:                                       # noqa: BLE001
            continue
        if p["reserved_all_ones"] and best is None:
            best = (mm, p)
    if best is None:
        best = (1, D4.parse_l1d(d, f, 1))
    out["mimo_mixed"], out["l1d"] = best
    out["ok"] = out["crc_ok"]
    return out


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("capture")
    ap.add_argument("--rate", type=float, required=True)
    ap.add_argument("--fmt")
    ap.add_argument("--frames", type=int, default=2)
    ap.add_argument("--step", type=float, default=0.5)
    ap.add_argument("--sweep-lr", action="store_true",
                    help="score the Block Interleaver column height Lr")
    ap.add_argument("--cliff", action="store_true",
                    help="MER profile of the last Preamble symbol")
    ap.add_argument("--json")
    a = ap.parse_args()
    path = a.capture if os.path.isabs(a.capture) else os.path.join(HERE,
                                                                   a.capture)
    print(f"M8 -- L1 from {os.path.basename(path)}\n" + "=" * 70)

    got = []
    for i in range(a.frames):
        try:
            r = solve(path, a.rate, a.fmt, i * a.step)
        except Exception as exc:                                # noqa: BLE001
            print(f"  frame {i}: skip ({exc})")
            continue
        if r is None:
            print(f"  frame {i}: L1-Basic did not verify")
            continue
        c = r["ctx"]
        g = r["geom"]
        match = g["cells"] == c["total_cells"]
        print(f"  frame {i} @{i*a.step:.1f}s: NP={c['NP']} "
              f"L1-Basic Mode {c['l1b_mode']} ({c['l1b_cells']} cells)  "
              f"L1-Detail Mode {g['mode']} Ksig {g['Ksig']}")
        print(f"           Nrepeat {g['Nrepeat']}  cells computed "
              f"{g['cells']} vs signalled {c['total_cells']}  "
              f"{'MATCH' if match else '*** MISMATCH ***'}")
        print(f"           LDPC {'CONVERGED' if r['fec']['converged'] else 'no'}"
              f" ({r['fec']['unsatisfied']} unsatisfied)  BCH "
              f"{'ZERO' if r['fec']['bch_ok'] else 'FAIL'}  L1D_crc "
              f"{'PASS' if r.get('crc_ok') else 'FAIL'}")
        if c["NP"] > 1:
            print(f"           extra Preamble symbols: cells {c['symbol_cells']}"
                  f" pilot-coherence "
                  f"{['%.4f' % v for v in c['extra_symbol_coherence']]}"
                  f"  spare in last {c['last_symbol_spare']}")
        if r["ok"]:
            got.append(r)
    if not got:
        print("\n  L1 did NOT verify on any Frame.  Reporting the wall.")
        return 1
    r = got[0]
    g, c = r["geom"], r["ctx"]

    # ---- control: the M4 reading MUST fail -------------------------------
    # The faithful M4 control is not "the same LLRs minus the repeat block":
    # it is M4's whole geometry, which believes the FEC Frame is Nfec/eta cells
    # long and reads only that many.  Anything softer is not the reading that
    # was actually wrong.
    if g["Nrepeat"]:
        g4 = D4.geometry(r["l1b"]["L1B_L1_Detail_size_bytes"],
                         r["l1b"]["L1B_L1_Detail_fec_type"])
        ctl = D4.decode(r["l1d_cells"][:g4["cells"]], g4, iters=100)
        print(f"\n  === control: M4's geometry, with no 6.5.2.7 repetition ===")
        print(f"    it reads {g4['cells']} cells where the transmitter "
              f"signals {c['total_cells']}")
        print(f"    conv={ctl['converged']} unsatisfied={ctl['unsatisfied']} "
              f"BCH={'ZERO' if ctl['bch_ok'] else 'fail'}   "
              + ("(control FAILED as required)" if not ctl["bch_ok"]
                 else "*** control PASSED -- this gate is vacuous ***"))

    # ---- gate: sweep the interleaver column height -----------------------
    # Scored under M4's geometry ON PURPOSE.  With repetition the transmitted
    # word is 7222 bits for 336 information bits -- rate 0.047 -- and the LDPC
    # shrugs off a permutation of half the cells, so the correct reading makes
    # this gate VACUOUS (measured: 10 of 13 wrong Lr values still decode).
    # A gate that cannot fail is not evidence, so it is run where it can.
    if a.sweep_lr and c["NP"] > 1:
        g4 = D4.geometry(r["l1b"]["L1B_L1_Detail_size_bytes"],
                         r["l1b"]["L1B_L1_Detail_fec_type"])
        lr0 = c["total_cells"] // c["NP"]
        print(f"\n  === gate: Block Interleaver column height Lr "
              f"(A/322 7.2.5.2 says floor(total/Lc) = {lr0}) ===")
        print(f"    scored under M4's geometry, where the LDPC still "
              f"discriminates.")
        rows = []
        for lr in range(lr0 - 6, lr0 + 7):
            xd = l1d_deinterleave(r["l1d_stream"], c["total_cells"], c["NP"],
                                  lr=lr)
            q = D4.decode(xd[:g4["cells"]], g4, iters=60)
            rows.append((lr, q["unsatisfied"]))
            print(f"    Lr {lr:6d}  unsatisfied {q['unsatisfied']:6d}")
        bl_ = min(rows, key=lambda t: t[1])
        others = [u for lr, u in rows if lr != bl_[0]]
        print(f"    -> best Lr {bl_[0]} at {bl_[1]} unsatisfied; every other "
              f"value {min(others)}..{max(others)}  "
              + ("(SPEC VALUE WINS)" if bl_[0] == lr0 else "*** NOT the spec value ***"))

    # ---- referee: the QPSK cliff in the last Preamble symbol -------------
    if a.cliff and c["NP"] > 1:
        _f, _s, _c2 = preamble_cells(path, a.rate, a.fmt, c["start_sec"])
        y, _fs, _cfo, _h = load(path, a.rate, fmt=a.fmt, span_sec=SPAN_SEC,
                                start_sec=c["start_sec"])
        s = c["NP"] - 1
        zs, _coh, _nd = demod_preamble_symbol(
            y, c["t0"] + (c["nfft"] + c["gi"]) * s, c["nfft"], c["gi"],
            c["dx"], r["l1b"]["L1B_preamble_reduced_carriers"], c["shift"])
        xs = FI.deinterleave(zs, c["nfft"], s, direction="forward", toggle="i")
        pred = c["total_cells"] - (c["available"] - len(xs))
        print(f"\n  === referee: the QPSK cliff in Preamble symbol {s} ===")
        print(f"    the mapping predicts L1-Detail ends at cell {pred}; "
              f"after it are data-PLP cells of a much higher order.")
        W = 200
        for k in range(0, len(xs) - W, W):
            m = mer_db(xs[k:k + W], QPSK)
            mk = "  <-- predicted end" if k <= pred < k + W else ""
            print(f"      {k:5d}-{k+W:5d}  {m:5.2f} dB  "
                  f"{'#' * max(0, int((m - 1) * 3))}{mk}")

    # ---- the configuration ----------------------------------------------
    print(f"\n  --- L1-DETAIL, decoded (L1B_first_sub_mimo_mixed="
          f"{r['mimo_mixed']}) ---")
    for pathname, name, v in r["l1d"]["fields"]:
        sem = LS.L1_DETAIL_SEMANTICS.get(name, {})
        lab = sem.get(v, "") if isinstance(sem, dict) else ""
        print(f"    {pathname:12s}{name:36s} {v:12d}   {lab}")
    print(f"    {'':12s}{'L1D_reserved':36s} {r['l1d']['reserved_len']:12d} "
          f"bits, all ones: {r['l1d']['reserved_all_ones']}")

    res = dict(capture=os.path.basename(path), rate=a.rate,
               l1_basic={k: int(v) for k, v in r["l1b"].items()},
               geometry={k: (v if isinstance(v, (int, float)) else str(v))
                         for k, v in g.items()},
               preamble={k: v for k, v in c.items() if k != "t0"},
               cells_match=bool(g["cells"] == c["total_cells"]),
               fec={k: v for k, v in r["fec"].items() if k != "nouter"},
               crc_ok=r["crc_ok"], mimo_mixed=r["mimo_mixed"],
               l1_detail=[{"path": p_, "name": n_, "value": int(v_)}
                          for p_, n_, v_ in r["l1d"]["fields"]],
               frames_ok=len(got))
    dest = a.json or os.path.join(
        HERE, "m8_l1_" + os.path.splitext(os.path.basename(path))[0] + ".json")
    with open(dest, "w") as fh:
        json.dump(res, fh, indent=2, default=float)
    print(f"\n  wrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
