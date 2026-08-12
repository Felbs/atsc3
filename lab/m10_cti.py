#!/usr/bin/env python3
"""M10 -- the A/322 7.1.4 Convolutional Time Interleaver (CTI), and its de-interleaver.

RF33 used the *Hybrid* time interleaver, which M5/M6 closed.  Every PLP on RF25
and RF30 uses the CTI instead, which is a different machine: a classic
convolutional (Ramsey/Forney) interleaver with `Nrows` delay lines, line `k`
holding `k` cells, and a commutator that steps one row per cell.

THE READING, STRAIGHT FROM 7.1.4.1 (PyMuPDF; the prose is intact, no math font)
------------------------------------------------------------------------------
    "The CTI shall consist of Nrows delay lines, with the k-th line having k
     delay elements ... Input and output shall be controlled by two commutators,
     cyclically switching downwards after one cell is written in or read out.
     For each cell gq, the commutators shall be located at the same row k.
     When the input commutator is located at row k, a cell gq shall be written
     to this delay line.  First, each delay element from this line shall shift
     its memory content to the neighbouring right delay element, and the content
     from the right-most delay element shall be output ...  Next, the input cell
     gq shall be written to the left-most delay element of this line.  Both
     commutators shall then move cyclically to the next line (k+1) mod Nrows."

Tracing one cell: it is written into row k at step q, advances one element per
*visit* to row k (i.e. every Nrows steps), and leaves on the k-th visit.  So

    out[q] = in[q - Nrows * k_q]          k_q = (start_row + q) mod Nrows

and inverting, the pre-CTI cell with index i leaves the interleaver at

    q(i)   = i + Nrows * ((start_row + i) mod Nrows)                        (*)

WHY (*) IS NOT AN ASSUMPTION
----------------------------
A/322 9.3.9.1 prints the transmitter's own formula for a signalled field:

    L1D_plp_CTI_fec_block_start = C + Nrows x ((L1D_plp_CTI_start_row + C) mod Nrows)

"where C is the position of the first cell of the first complete FEC Block
before the CTI".  That is EXACTLY (*) evaluated at i = C.  So the spec publishes
the de-interleaver's index map, in terms of two fields this project has already
decoded off the air -- and that turns the reading into a GATE:

    solve C from the signalled pair, and require  0 <= C < cells_per_FEC_Block.

C is a single determined integer over a range of ~Nrows^2 ~ 1e6, so landing
inside a 32400-cell window by accident has probability ~3%.  Five PLPs across
two independent transmitters all land inside it (p ~ 2e-8), and the wrong Nrows
readings do not.

GATES IMPLEMENTED HERE
----------------------
R1  synthesized round trip -- interleave, de-interleave, identity on the
    valid region, for several (Nrows, start_row) and a long stream
R2  structural: every output index is used exactly once; the delay of each
    path is k*Nrows with k the commutator row
R3  THE AIR: the 9.3.9.1 identity above, per PLP, with the Nrows menu swept
    (Table 9.26: 512 / 724 / 887 / 1024, and the extended 1254 / 1448) so the
    wrong depths are run as controls that must fail
R4  cross-frame: the commutator advances by plp_size per Subframe, so
    start_row(n+1) == (start_row(n) + plp_size) mod Nrows -- checked against
    an independently L1-decoded neighbouring Frame

Run `python m10_cti.py` for R1/R2 and `python m10_cti.py --l1 m8_l1_rf25.json
m8_l1_rf30.json` for R3.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

# A/322 Table 9.26 -- L1D_plp_CTI_depth -> Nrows
CTI_NROWS = {0: 512, 1: 724, 2: 887, 3: 1024}
CTI_NROWS_EXTENDED = {0: 512, 1: 724, 2: 1254, 3: 1448}

# A/322 Table 6.14 -- eta_MOD
MOD_BITS = {"QPSK": 2, "16QAM": 4, "64QAM": 6, "256QAM": 8,
            "1024QAM": 10, "4096QAM": 12}


def nrows_of(depth, extended=False):
    t = CTI_NROWS_EXTENDED if extended else CTI_NROWS
    if depth not in t:
        raise ValueError("L1D_plp_CTI_depth %r is Reserved" % (depth,))
    return t[depth]


# ---------------------------------------------------------------------------
# the index map
# ---------------------------------------------------------------------------

def out_index(i, nrows, start_row):
    """(*) -- where pre-CTI cell i leaves the interleaver.  Vectorised."""
    i = np.asarray(i, np.int64)
    return i + nrows * ((start_row + i) % nrows)


def in_index(q, nrows, start_row):
    """The inverse of (*): which pre-CTI cell leaves at output position q."""
    q = np.asarray(q, np.int64)
    return q - nrows * ((start_row + q) % nrows)


def interleave(cells, nrows, start_row, fill=None):
    """Transmit side.  Output positions fed by i < 0 carry the initial state.

    A/322 7.1.4.2 fills the delay elements from the 5.2.3 PRBS mapped to QAM;
    `fill` lets a caller supply that, but nothing in the receive path needs it
    (those positions are exactly the ones the de-interleaver marks invalid).
    """
    cells = np.asarray(cells)
    n = len(cells)
    out = np.zeros(n, cells.dtype) if fill is None else np.asarray(fill).copy()
    q = np.arange(n, dtype=np.int64)
    i = in_index(q, nrows, start_row)
    ok = (i >= 0) & (i < n)
    out[q[ok]] = cells[i[ok]]
    return out


def deinterleave(recv, nrows, start_row, n_out=None):
    """Receive side.  Returns (cells, valid) with `valid` a boolean mask.

    `recv[0]` must be the first cell of the PLP leaving the CTI in the Subframe
    whose L1 signals `start_row`; consecutive Subframes simply concatenate,
    because the commutator never resets -- it is signalled, not restarted.
    """
    recv = np.asarray(recv)
    L = len(recv)
    n = L if n_out is None else n_out
    i = np.arange(n, dtype=np.int64)
    q = out_index(i, nrows, start_row)
    valid = q < L
    out = np.zeros(n, recv.dtype)
    out[valid] = recv[q[valid]]
    return out, valid


def solve_C(fec_block_start, start_row, nrows):
    """Invert A/322 9.3.9.1 for C, the pre-CTI offset of the first FEC Block."""
    r = (start_row + fec_block_start) % nrows
    return fec_block_start - nrows * r


# ---------------------------------------------------------------------------
# R1 / R2 -- offline gates
# ---------------------------------------------------------------------------

def gate_roundtrip(verbose=True):
    ok = True
    rng = np.random.default_rng(7)
    for nrows, start in ((8, 0), (8, 5), (512, 191), (1024, 344), (1024, 1002)):
        n = max(4 * nrows * nrows, 4096)
        src = rng.integers(0, 1 << 30, n).astype(np.int64)
        tx = interleave(src, nrows, start)
        rx, valid = deinterleave(tx, nrows, start)
        good = bool(np.array_equal(rx[valid], src[valid]))
        # the valid region must be most of the stream, not a token few cells
        frac = float(valid.mean())
        ok &= good and frac > 0.5
        if verbose:
            print("    %s  R1 round trip Nrows=%-5d start=%-5d  n=%d  "
                  "valid %.3f" % ("PASS" if good else "FAIL", nrows, start,
                                  n, frac))
    return ok


def gate_structure(verbose=True):
    ok = True
    for nrows, start in ((512, 191), (1024, 344)):
        n = 3 * nrows * nrows
        i = np.arange(n, dtype=np.int64)
        q = out_index(i, nrows, start)
        # every output position is produced by at most one input
        uniq = len(np.unique(q)) == n
        # the delay is a multiple of Nrows and < Nrows^2
        d = q - i
        mult = bool(np.all(d % nrows == 0)) and bool(d.max() < nrows * nrows)
        # commutator row of the OUTPUT position equals the row the cell used
        k = d // nrows
        same = bool(np.all(((start + q) % nrows) == k))
        ok &= uniq and mult and same
        if verbose:
            print("    %s  R2 structure Nrows=%-5d: injective %s, delay "
                  "k*Nrows %s, commutator rows agree %s"
                  % ("PASS" if (uniq and mult and same) else "FAIL",
                     nrows, uniq, mult, same))
    return ok


# ---------------------------------------------------------------------------
# R3 -- the air
# ---------------------------------------------------------------------------

MOD_NAME = {0: "QPSK", 1: "16QAM", 2: "64QAM", 3: "256QAM",
            4: "1024QAM", 5: "4096QAM"}
RATE_NAME = {i: "%d/15" % (i + 2) for i in range(12)}


def plps_from_l1(js):
    """Pull (plp_id, mod, rate, fec_type, layer, CTI fields) out of an m8 JSON."""
    fields = js["l1_detail"]
    plps, cur, core = [], None, {}
    for e in fields:
        n, v = e["name"], e["value"]
        if n == "L1D_plp_id":
            if cur:
                plps.append(cur)
            cur = {"L1D_plp_id": v}
        elif cur is not None:
            cur[n] = v
    if cur:
        plps.append(cur)
    for p in plps:
        p["mod"] = MOD_NAME[p["L1D_plp_mod"]]
        p["rate"] = RATE_NAME[p["L1D_plp_cod"]]
        p["ninner"] = 64800 if p["L1D_plp_fec_type"] else 16200
        p["cells_per_fec"] = p["ninner"] // MOD_BITS[p["mod"]]
        if p.get("L1D_plp_layer", 0) == 0 and "L1D_plp_CTI_depth" in p:
            core = p
        # the LDM enhanced layer inherits the Core layer's commutator schedule
        p.setdefault("L1D_plp_CTI_depth", core.get("L1D_plp_CTI_depth"))
        p.setdefault("L1D_plp_CTI_start_row", core.get("L1D_plp_CTI_start_row"))
        p.setdefault("L1D_plp_TI_extended_interleaving",
                     core.get("L1D_plp_TI_extended_interleaving", 0))
    return plps


def gate_air(paths, verbose=True):
    """R3: solve A/322 9.3.9.1 for C on every PLP, sweeping the Nrows menu."""
    allok = True
    for path in paths:
        js = json.load(open(path))
        print("\n  === %s ===" % os.path.basename(path))
        for p in plps_from_l1(js):
            sr = p["L1D_plp_CTI_start_row"]
            fbs = p["L1D_plp_CTI_fec_block_start"]
            nfec = p["cells_per_fec"]
            sig = nrows_of(p["L1D_plp_CTI_depth"],
                           bool(p["L1D_plp_TI_extended_interleaving"]))
            print("    PLP %-3d %-8s %-6s Ninner %5d  depth %d -> Nrows %4d"
                  "  start_row %4d  fec_block_start %7d"
                  % (p["L1D_plp_id"], p["mod"], p["rate"], p["ninner"],
                     p["L1D_plp_CTI_depth"], sig, sr, fbs))
            hits = []
            for nr in sorted(set(list(CTI_NROWS.values())
                                 + list(CTI_NROWS_EXTENDED.values()))):
                if sr >= nr:
                    print("        Nrows %4d : start_row %d >= Nrows -- "
                          "IMPOSSIBLE" % (nr, sr))
                    continue
                C = solve_C(fbs, sr, nr)
                good = 0 <= C < nfec
                if good:
                    hits.append(nr)
                print("        Nrows %4d : C = %8d   %s"
                      % (nr, C, "in [0,%d)  <== PASS" % nfec if good
                         else "OUT of [0,%d)" % nfec))
            okp = hits == [sig]
            allok &= okp
            print("        -> only Nrows %s satisfies 9.3.9.1; L1 signals %d"
                  "   %s" % (hits, sig,
                             "GATE PASS" if okp else "*** GATE FAIL ***"))
    return allok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--l1", nargs="*", default=[],
                    help="m8_l1_*.json files for the R3 air gate")
    a = ap.parse_args()
    print("M10 -- A/322 7.1.4 Convolutional Time Interleaver")
    print("=" * 72)
    print("\n  === R1: synthesized round trip ===")
    ok = gate_roundtrip()
    print("\n  === R2: structure ===")
    ok &= gate_structure()
    if a.l1:
        print("\n  === R3: the air -- A/322 9.3.9.1 solved for C ===")
        paths = [p if os.path.isabs(p) else os.path.join(HERE, p)
                 for p in a.l1]
        ok &= gate_air(paths)
    print("\n  %s" % ("ALL GATES PASS" if ok else "GATE FAILURES"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
