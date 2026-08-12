#!/usr/bin/env python3
"""M4 Step 2 -- L1-Detail off real air: the PLP configuration.

With the A/322 5.2.3 scrambler pinned (m4_scrambler.py), L1-Basic's field
VALUES are readable, and they hand over every parameter L1-Detail needs:

    L1B_L1_Detail_size_bytes  = 64   -> Ksig 512, Nouter 680
    L1B_L1_Detail_fec_type    = 2    -> Mode 3: Kldpc 6480, Ninner 16200,
                                        rate 6/15 (Type B), QPSK, eta 2
    L1B_L1_Detail_total_cells = 880
    L1B_preamble_num_symbols  = 0    -> NP = 1 preamble symbol, so the
                                        7.2.5.2 L1-Detail block interleaver
                                        is Lc = 1 column, i.e. IDENTITY
    L1B_L1_Detail_additional_parity_mode = 0 -> no additional parity cells

so L1-Detail is exactly cells 484 .. 1363 of the frequency-de-interleaved
first preamble symbol, and nothing about its geometry is guessed.

THE REFEREES (three, each with controls that must fail)
-------------------------------------------------------
1. LDPC rate 6/15 must converge -- 9720 simultaneous parity checks.
2. The 168-bit BCH syndrome over the 680 Nouter bits must be zero.
3. A/322 9.3's L1D_crc must check over the descrambled 480 + 32 bits.  This
   one is decisive: unlike L1-Basic's, it is evaluated with a KNOWN scrambler,
   so it tests the scrambler solution again on completely different data.

THE WAITING CROSS-CHECK
-----------------------
M3 measured the QPSK MER cliff at cell ~1320 and predicted
L1B_L1_Detail_total_cells ~= 1320 - 484 = 836.  The decoded field says 880
(end of L1-Detail at cell 1364).  This tool re-measures the cliff finely so
that prediction can be scored honestly rather than quietly dropped.

Usage:
    python m4_l1detail.py hit_rf33.cs16 --rate 8e6
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
import spec_bicm as BI                                          # noqa: E402
import spec_l1syntax as LS                                      # noqa: E402
import m4_scrambler as SC                                       # noqa: E402
from m3_preamble import analyse, QPSK, mer_db, m4_stat          # noqa: E402

L1B_CELLS = 484                      # A/322 Table 6.17, L1-Basic Mode 3


# --------------------------------------------------------------------------
# L1-Detail FEC geometry, straight from the decoded L1-Basic fields
# --------------------------------------------------------------------------

def geometry(size_bytes, fec_type):
    """A/322 6.5.2.1/6.5.2.4/6.5.2.8 for L1-Detail.  Nothing here is fitted."""
    mode = fec_type + 1
    kind = "L1-Detail-%d" % mode
    p = BI.PUNCTURING[kind]
    kldpc = BI.L1D_KLDPC[mode]
    kseg = BI.L1D_KSEG[mode]
    ksig = size_bytes * 8
    if ksig > kseg:
        raise SystemExit(f"Ksig {ksig} > Kseg {kseg}: L1-Detail is segmented "
                         f"into multiple FEC blocks -- not implemented")
    nouter = ksig + 168
    npunc_t = int(p["A"] * (kldpc - nouter)) + p["B"]
    nfec_t = nouter + p["Nldpc_parity"] - npunc_t
    eta = p["eta"]
    nfec = (nfec_t // eta) * eta
    return dict(mode=mode, kind=kind, Ksig=ksig, Nouter=nouter, Kldpc=kldpc,
                Ninner=16200, Nldpc_parity=p["Nldpc_parity"], eta=eta,
                Nfec=nfec, Npunc=npunc_t - (nfec - nfec_t),
                n_parity_kept=nfec - nouter, cells=nfec // eta,
                rate={1: "3/15", 2: "3/15"}.get(mode, "6/15"))


def build_llr(cells, g, shorten_pattern=None, gw=None, gw_invert=False):
    bl = L1.cells_to_llr(cells[:g["cells"]], nfec=g["Nfec"], eta=g["eta"])
    llr = np.zeros(g["Ninner"])
    ipos, padded = L1.info_positions(
        g["Kldpc"], g["Nouter"],
        pattern=list(shorten_pattern or BI.SHORTENING_PATTERN[g["kind"]]))
    llr[:g["Kldpc"]][padded] = L1.BIG
    llr[ipos] = bl[:g["Nouter"]]
    ppos = L1.parity_positions(g["n_parity_kept"], g["Kldpc"], g["Ninner"],
                               groupwise=list(gw or BI.GROUPWISE[g["kind"]]),
                               invert=gw_invert)
    llr[ppos] = bl[g["Nouter"]:]
    return llr, ipos


def decode(cells, g, iters=100, gw_invert=False, shorten_pattern=None,
           gw=None, _cache={}):
    llr, ipos = build_llr(cells, g, shorten_pattern, gw, gw_invert)
    key = (g["rate"], g["Ninner"])
    if key not in _cache:
        _cache[key] = LD.parity_check(g["rate"], n=g["Ninner"])
    checks = _cache[key]
    bits, conv, it, bad = LD.min_sum_decode(llr, checks, g["Ninner"],
                                            iters=iters)
    nouter = bits[ipos]
    syn = L1.bch_syndrome(nouter)
    return dict(converged=bool(conv), iters=int(it), unsatisfied=int(bad),
                bch_ok=syn == 0, bch_syndrome=int(syn), nouter=nouter)


# --------------------------------------------------------------------------
# A/322 Table 9.9 syntax -- CORRECTED
#
# The extracted `spec_l1syntax.L1_DETAIL_FIELDS` drops four branches, because
# the PDF prints `}else {` on its own line and the extractor treated the
# following block as unconditional or absent.  Running the tree as extracted
# derails immediately after the first PLP (PLP 2 came back with
# scrambler_type = 3 "Reserved" and fec_type = 8 "Reserved", i.e. nonsense).
# The missing branches, restored verbatim from A/322:2026-04 Table 9.9:
#
#   1. else if (L1D_plp_TI_mode = 01) { L1D_plp_CTI_fec_block_start  22 }
#   2. else if (L1D_plp_TI_mode = 10) { the six HTI fields }
#   3. else                           { L1D_plp_ldm_injection_level   5 }
#      (the else of `if (L1D_plp_layer = 0)` -- Enhanced-layer LDM PLPs)
#   4. the whole second `for i` loop is the "L1D_version >= 2 loop"
#      (A/322 9.3 footnote 2) and is ABSENT when L1D_version < 2.
#
# ASSUMPTION D1: inside the HTI branch the PDF's brace placement leaves it
# ambiguous whether L1D_plp_HTI_cell_interleaver sits inside the
# inter_subframe != 0 arm or after the if/else.  Read here as AFTER, i.e.
# always present in HTI mode.  RISK LOW -- it is 1 bit and the all-ones
# L1D_reserved run arbitrates: a wrong reading shifts every later field.
# --------------------------------------------------------------------------

_HTI = [
    ("field", "L1D_plp_HTI_inter_subframe", 1),
    ("field", "L1D_plp_HTI_num_ti_blocks", 4),
    ("field", "L1D_plp_HTI_num_fec_blocks_max", 12),
    ("if", "L1D_plp_HTI_inter_subframe == 0",
     [("field", "L1D_plp_HTI_num_fec_blocks", 12)]),
    ("if", "L1D_plp_HTI_inter_subframe != 0",
     [("for", "k", "0 .. L1D_plp_HTI_num_ti_blocks",
       [("field", "L1D_plp_HTI_num_fec_blocks", 12)])]),
    ("field", "L1D_plp_HTI_cell_interleaver", 1),
]

_PLP = [
    ("field", "L1D_plp_id", 6),
    ("field", "L1D_plp_lls_flag", 1),
    ("field", "L1D_plp_layer", 2),
    ("field", "L1D_plp_start", 24),
    ("field", "L1D_plp_size", 24),
    ("field", "L1D_plp_scrambler_type", 2),
    ("field", "L1D_plp_fec_type", 4),
    ("if", "L1D_plp_fec_type in (0,1,2,3,4,5)",
     [("field", "L1D_plp_mod", 4), ("field", "L1D_plp_cod", 4)]),
    ("field", "L1D_plp_TI_mode", 2),
    ("if", "L1D_plp_TI_mode == 0", [("field", "L1D_plp_fec_block_start", 15)]),
    ("if", "L1D_plp_TI_mode == 1",
     [("field", "L1D_plp_CTI_fec_block_start", 22)]),          # RESTORED
    ("if", "L1D_num_rf > 0", [
        ("field", "L1D_plp_num_channel_bonded", 3),
        ("if", "L1D_plp_num_channel_bonded > 0", [
            ("field", "L1D_plp_channel_bonding_format", 2),
            ("for", "k", "0 .. L1D_plp_num_channel_bonded",
             [("field", "L1D_plp_bonded_rf_id", 3)])])]),
    ("if", "(i == 0 and L1B_first_sub_mimo == 1) or (i > 0 and L1D_mimo == 1)",
     [("field", "L1D_plp_mimo_stream_combining", 1),
      ("field", "L1D_plp_mimo_IQ_interleaving", 1),
      ("field", "L1D_plp_mimo_PH", 1)]),
    ("if", "L1D_plp_layer == 0", [
        ("field", "L1D_plp_type", 1),
        ("if", "L1D_plp_type == 1",
         [("field", "L1D_plp_num_subslices", 14),
          ("field", "L1D_plp_subslice_interval", 24)]),
        ("if", "L1D_plp_TI_mode in (1,2) and L1D_plp_mod == 0",
         [("field", "L1D_plp_TI_extended_interleaving", 1)]),
        ("if", "L1D_plp_TI_mode == 1",
         [("field", "L1D_plp_CTI_depth", 3),
          ("field", "L1D_plp_CTI_start_row", 11)]),
        ("if", "L1D_plp_TI_mode == 2", _HTI)]),                # RESTORED
    ("if", "L1D_plp_layer != 0",
     [("field", "L1D_plp_ldm_injection_level", 5)]),           # RESTORED
]

L1D_SYNTAX = [
    ("field", "L1D_version", 4),
    ("field", "L1D_num_rf", 3),
    ("for", "rf", "1 .. L1D_num_rf",
     [("field", "L1D_bonded_bsid", 16), ("field", "L1D_rf_reserved", 3)]),
    ("if", "L1B_time_info_flag != 0", [
        ("field", "L1D_time_sec", 32),
        ("field", "L1D_time_msec", 10),
        ("if", "L1B_time_info_flag != 1", [
            ("field", "L1D_time_usec", 10),
            ("if", "L1B_time_info_flag != 2",
             [("field", "L1D_time_nsec", 10)])])]),
    ("for", "i", "0 .. L1B_num_subframes", [
        ("if", "i > 0", [
            ("field", "L1D_mimo", 1), ("field", "L1D_miso", 2),
            ("field", "L1D_fft_size", 2),
            ("field", "L1D_reduced_carriers", 3),
            ("field", "L1D_guard_interval", 4),
            ("field", "L1D_num_ofdm_symbols", 11),
            ("field", "L1D_scattered_pilot_pattern", 5),
            ("field", "L1D_scattered_pilot_boost", 3),
            ("field", "L1D_sbs_first", 1), ("field", "L1D_sbs_last", 1)]),
        ("if", "L1B_num_subframes > 0",
         [("field", "L1D_subframe_multiplex", 1)]),
        ("field", "L1D_frequency_interleaver", 1),
        ("if", "(i == 0 and (L1B_first_sub_sbs_first or "
                "L1B_first_sub_sbs_last)) or (i > 0 and (L1D_sbs_first or "
                "L1D_sbs_last))", [("field", "L1D_sbs_null_cells", 13)]),
        ("field", "L1D_num_plp", 6),
        ("for", "j", "0 .. L1D_num_plp", _PLP)]),
    ("field", "L1D_bsid", 16),
    ("if", "L1D_version >= 2", [                               # RESTORED gate
        ("for", "i", "0 .. L1B_num_subframes", [
            ("if", "i > 0", [("field", "L1D_mimo_mixed", 1)]),
            ("if", "(i == 0 and L1B_first_sub_mimo_mixed == 1) or "
                   "(i > 0 and L1D_mimo_mixed == 1)",
             [("for", "j", "0 .. L1D_num_plp", [
                 ("field", "L1D_plp_mimo", 1),
                 ("if", "L1D_plp_mimo == 1",
                  [("field", "L1D_plp_mimo_stream_combining", 1),
                   ("field", "L1D_plp_mimo_IQ_interleaving", 1),
                   ("field", "L1D_plp_mimo_PH", 1)])])])])]),
]


class Parser:
    """Evaluate the extracted L1-Detail syntax tree against a bit array."""

    def __init__(self, bits, l1b):
        self.b = bits
        self.i = 0
        self.env = dict(l1b)
        self.out = []                    # (path, name, value)
        self.loop = {}

    def _read(self, name, w):
        if self.i + w > len(self.b):
            raise IndexError(f"{name}: ran off the end at bit {self.i}")
        v = 0
        for k in range(w):
            v = (v << 1) | int(self.b[self.i + k])
        self.i += w
        return v

    def _cond(self, expr):
        e = (expr.replace("0b00", "0").replace("0b01", "1")
                 .replace("0b10", "2").replace("0b11", "3")
                 .replace("0b0000", "0").replace("!=", "!=").replace("=", "="))
        e = expr
        for a, b in (("0b0000", "0"), ("0b00", "0"), ("0b01", "1"),
                     ("0b10", "2"), ("0b11", "3")):
            e = e.replace(a, b)
        env = dict(self.env)
        env.update(self.loop)
        try:
            return bool(eval(e, {"__builtins__": {}}, env))    # noqa: S307
        except Exception as exc:                                # noqa: BLE001
            raise SystemExit(f"cannot evaluate condition {expr!r}: {exc}")

    def run(self, items, path=""):
        for it in items:
            if it[0] == "field":
                _, name, w = it
                if w == "as needed":
                    continue
                v = self._read(name, w)
                self.out.append((path, name, v))
                self.env[name] = v
                self.loop[name] = v
            elif it[0] == "if":
                if self._cond(it[1]):
                    self.run(it[2], path)
            elif it[0] == "for":
                _, var, rng, body = it
                lo_s, hi_s = [x.strip() for x in rng.split("..")]
                env = dict(self.env)
                env.update(self.loop)
                lo = int(eval(lo_s, {"__builtins__": {}}, env))   # noqa: S307
                hi = int(eval(hi_s, {"__builtins__": {}}, env))   # noqa: S307
                for n in range(lo, hi + 1):
                    self.loop[var] = n
                    self.env[var] = n
                    self.run(body, f"{path}{var}={n}/")
            else:
                raise SystemExit(f"unknown syntax node {it[0]}")


def parse_l1d(bits, l1b, mimo_mixed):
    """Parse, then use the trailing reserved bits as the arbiter.

    A/322 3.2.1: reserved bits default to 1.  L1D_reserved is 'as needed'
    padding between the last real field and the 32-bit L1D_crc, so a correct
    parse leaves a run of ONES there and an incorrect one does not.  That makes
    the reserved run a genuine test rather than filler.
    """
    env = dict(l1b)
    env["L1B_first_sub_mimo_mixed"] = mimo_mixed
    p = Parser(bits, env)
    p.run(L1D_SYNTAX)
    used = p.i
    crc_bits = bits[-32:]
    resv = bits[used:len(bits) - 32]
    return dict(fields=p.out, bits_used=used, reserved=resv,
                reserved_all_ones=bool(len(resv) and np.all(resv == 1)),
                reserved_len=int(len(resv)),
                crc=int("".join(str(int(c)) for c in crc_bits), 2))


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def l1basic_of(x, mode, repeat=True):
    # E82: `repeat` carries A/322 6.5.2.7's parity repetition, which applies
    # to L1-Basic Mode 1 (Nrepeat = 3672 of 7640 transmitted bits) and to no
    # other Mode.  repeat=False is the pre-E82 reading, kept as a control.
    r = L1.decode_one(x, mode, variant="standard", gw_invert=False, iters=100,
                      repeat=repeat)
    if not (r["converged"] and r["bch_ok"]):
        return None
    d = SC.descramble(np.array(r["nouter_bits"][:200], np.uint8))
    K = CRC.crc32(np.zeros(168, np.uint8), init=0xFFFFFFFF)
    if CRC.g_stat(d) != K:
        return None
    return SC.parse_l1b(d)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("capture")
    ap.add_argument("--rate", type=float, required=True)
    ap.add_argument("--fmt")
    ap.add_argument("--frames", type=int, default=3)
    ap.add_argument("--step", type=float, default=0.7)
    ap.add_argument("--json")
    a = ap.parse_args()
    path = a.capture if os.path.isabs(a.capture) else os.path.join(HERE,
                                                                   a.capture)
    print("L1-DETAIL from real air\n" + "=" * 70)

    frames = []
    for i in range(a.frames):
        try:
            rep, z, Y, geo = analyse(path, a.rate, a.fmt, report={},
                                     start_sec=i * a.step, quiet=True)
        except Exception as exc:                                # noqa: BLE001
            print(f"  frame {i}: skip ({exc})")
            continue
        mode = geo[-1]
        nfft = geo[-4]
        x = FI.deinterleave(z, nfft, 0, direction="forward", toggle="i")
        f = l1basic_of(x, mode)
        if f is None:
            print(f"  frame {i}: L1-Basic did not verify -- skipped")
            continue
        frames.append((i * a.step, x, f))
    if not frames:
        raise SystemExit("no verified frames")
    t0, x0, l1b = frames[0]
    print(f"  {len(frames)} frames with a CRC-clean, descrambled L1-Basic\n")

    g = geometry(l1b["L1B_L1_Detail_size_bytes"], l1b["L1B_L1_Detail_fec_type"])
    print(f"  L1-Detail geometry, entirely from L1-Basic field values:")
    for k in ("mode", "Ksig", "Nouter", "Kldpc", "Ninner", "rate", "eta",
              "Nfec", "Npunc", "n_parity_kept", "cells"):
        print(f"    {k:16s} {g[k]}")
    tot = l1b["L1B_L1_Detail_total_cells"]
    print(f"    L1B_L1_Detail_total_cells (signalled) {tot}  "
          f"{'MATCH' if tot == g['cells'] else 'MISMATCH'}")
    npre = l1b["L1B_preamble_num_symbols"] + 1
    print(f"    NP = {npre} preamble symbol(s) -> 7.2.5.2 interleaver has "
          f"Lc = {npre} column(s)"
          + ("  (identity)" if npre == 1 else "  -- NOT IMPLEMENTED"))
    if npre != 1:
        raise SystemExit("multi-symbol preamble de-interleave not implemented")
    lo, hi = L1B_CELLS, L1B_CELLS + tot

    # ---- the cliff cross-check M3 left waiting ---------------------------
    print(f"\n  === the cross-check M3 left waiting ===")
    print(f"    M3 measured the QPSK MER cliff at cell ~1320 and predicted")
    print(f"    L1B_L1_Detail_total_cells ~= 1320 - 484 = 836.")
    print(f"    DECODED VALUE: {tot}  -> L1-Detail ends at cell {hi}.")
    print(f"    prediction error {836 - tot:+d} cells "
          f"({100*(836-tot)/tot:+.1f}%)\n")
    win = 60
    print(f"    fine re-measure of the cliff (MER of {win}-cell windows):")
    prof = []
    for c in range(1150, min(1550, len(x0)) - win, 20):
        m = mer_db(x0[c:c + win], QPSK)
        prof.append((c, m))
    for c, m in prof:
        bar = "#" * max(0, int((m - 2) * 2))
        mark = "  <-- signalled end of L1-Detail" if c <= hi < c + 20 else ""
        print(f"      cells {c:5d}-{c+win:5d}  {m:5.2f} dB  {bar}{mark}")

    # ---- decode ----------------------------------------------------------
    print(f"\n  === decoding L1-Detail: cells {lo}..{hi-1} ===")
    rows = []
    for (t, x, f) in frames:
        r = decode(x[lo:hi], g)
        rows.append(dict(t=t, converged=r["converged"], iters=r["iters"],
                         unsatisfied=r["unsatisfied"], bch_ok=r["bch_ok"]))
        print(f"    @{t:4.1f}s  LDPC {'CONVERGED' if r['converged'] else 'no'}"
              f" in {r['iters']:3d} it, {r['unsatisfied']:5d} unsatisfied of "
              f"{g['Ninner'] - g['Kldpc']}   BCH "
              f"{'ZERO -> PASS' if r['bch_ok'] else hex(r['bch_syndrome'])[:14]}")
        if t == t0:
            first = r

    # ---- controls: each MUST fail ---------------------------------------
    print(f"\n  === controls (each must FAIL) ===")
    ctl = []
    for name, cells, kw in (
            ("cells shifted +40", x0[lo + 40:hi + 40], {}),
            ("cells shifted -40", x0[lo - 40:hi - 40], {}),
            ("L1-Basic's own cells", x0[0:tot], {}),
            ("tail of the symbol", x0[-tot:], {}),
            ("wrong shortening pattern (L1-Detail-4)", x0[lo:hi],
             dict(shorten_pattern=BI.SHORTENING_PATTERN["L1-Detail-4"])),
            ("wrong group-wise pattern (L1-Detail-4)", x0[lo:hi],
             dict(gw=BI.GROUPWISE["L1-Detail-4"])),
            ("group-wise permutation inverted", x0[lo:hi],
             dict(gw_invert=True)),
    ):
        r = decode(cells, g, iters=50, **kw)
        ctl.append(dict(name=name, converged=r["converged"],
                        unsatisfied=r["unsatisfied"], bch_ok=r["bch_ok"]))
        print(f"    {name:42s} conv={str(r['converged']):5s} "
              f"unsat={r['unsatisfied']:5d}  BCH="
              f"{'PASS' if r['bch_ok'] else 'fail'}")

    if not (first["converged"] and first["bch_ok"]):
        print("\n  L1-Detail did NOT verify.  Reporting the wall.")
        return 1

    # ---- descramble + CRC ------------------------------------------------
    print(f"\n  === descramble (A/322 6.5.2.2 = the 5.2.3 scrambler, "
          f"restarted) + L1D_crc ===")
    out_frames = []
    for (t, x, f), row in zip(frames, rows):
        r = decode(x[lo:hi], g)
        if not (r["converged"] and r["bch_ok"]):
            continue
        d = SC.descramble(np.array(r["nouter"][:g["Ksig"]], np.uint8))
        ok = CRC.crc32(d[:g["Ksig"] - 32], init=0xFFFFFFFF) == int(
            "".join(str(int(c)) for c in d[g["Ksig"] - 32:]), 2)
        print(f"    @{t:4.1f}s  L1D_crc over {g['Ksig']-32} bits -> "
              f"{'PASS' if ok else 'FAIL'}")
        out_frames.append((t, d, f, ok))
    if not out_frames or not out_frames[0][3]:
        print("  L1D_crc failed -- stopping rather than reporting fields.")
        return 1

    # ---- parse -----------------------------------------------------------
    t, d, f, _ = out_frames[0]
    print(f"\n  === A/322 9.3 parse.  L1B_first_sub_mimo_mixed is the one "
          f"reading the\n      L1-Basic block cannot settle (it abuts the "
          f"all-ones reserved run),\n      so BOTH are run and the trailing "
          f"reserved bits arbitrate ===")
    best = None
    for mm in (1, 0):
        try:
            p = parse_l1d(d, f, mm)
        except Exception as exc:                                # noqa: BLE001
            print(f"    mimo_mixed={mm}: parse failed ({exc})")
            continue
        print(f"    mimo_mixed={mm}: {p['bits_used']:4d} bits used, "
              f"{p['reserved_len']:3d} reserved bits, all ones: "
              f"{p['reserved_all_ones']}")
        if p["reserved_all_ones"] and best is None:
            best = (mm, p)
    if best is None:
        print("    NEITHER reading leaves an all-ones reserved run -- "
              "reporting the ambiguity.")
        mm, p = 1, parse_l1d(d, f, 1)
    else:
        mm, p = best
        print(f"    -> L1B_first_sub_mimo_mixed = {mm}")

    # ---- external clock cross-check -------------------------------------
    fv = {n: v for _, n, v in p["fields"]}
    if "L1D_time_sec" in fv:
        import datetime
        TAI_UTC = 37                     # seconds, constant since 2017-01-01
        t_tai = fv["L1D_time_sec"]
        t_utc = t_tai - TAI_UTC
        mt = os.stat(path).st_mtime
        print(f"\n  === EXTERNAL cross-check: the transmitter's clock vs "
              f"the capture file's ===")
        print(f"    A/322 9.3: L1D_time_sec is the TAI seconds at the first "
              f"sample of the bootstrap.")
        print(f"    decoded  {t_tai} TAI  -> {t_utc} UTC = "
              f"{datetime.datetime.fromtimestamp(t_utc, datetime.timezone.utc)}")
        print(f"    capture file mtime      {mt:.1f} = "
              f"{datetime.datetime.fromtimestamp(mt, datetime.timezone.utc)}"
              f"   (8 s record, so it began ~{mt-8:.1f})")
        print(f"    decoded frame time - capture start = {t_utc - (mt-8):+.1f} s"
              f"   (must land inside the 8 s record)")
        print(f"    sub-second: .{fv.get('L1D_time_msec',0):03d}"
              f"{fv.get('L1D_time_usec',0):03d}"
              f"{fv.get('L1D_time_nsec',0):03d}")

    print(f"\n  --- L1-DETAIL, decoded ---")
    for pathname, name, v in p["fields"]:
        sem = LS.L1_DETAIL_SEMANTICS.get(name, {})
        lab = sem.get(v, "") if isinstance(sem, dict) else ""
        print(f"    {pathname:12s}{name:34s} {v:12d}   {lab}")
    print(f"    {'':12s}{'L1D_reserved':34s} {p['reserved_len']:12d} bits, "
          f"all ones: {p['reserved_all_ones']}")
    print(f"    {'':12s}{'L1D_crc':34s} {p['crc']:12d}   VERIFIED")

    res = {"capture": os.path.basename(path), "geometry": {
        k: (str(v) if not isinstance(v, (int, float)) else v)
        for k, v in g.items()},
        "l1_basic": {k: int(v) for k, v in l1b.items()},
        "frames": rows, "controls": ctl,
        "cliff_profile": prof,
        "mimo_mixed": mm,
        "l1_detail": [{"path": a_, "name": b_, "value": int(c_)}
                      for a_, b_, c_ in p["fields"]],
        "reserved_len": p["reserved_len"],
        "reserved_all_ones": p["reserved_all_ones"],
        "l1d_bits": "".join(str(int(c)) for c in d)}
    dest = a.json or os.path.join(
        HERE, "m4_l1detail_" + os.path.splitext(os.path.basename(path))[0]
        + ".json")
    with open(dest, "w") as fh:
        json.dump(res, fh, indent=2, default=float)
    print(f"\n  wrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
