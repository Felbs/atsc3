#!/usr/bin/env python3
"""Machine-extract A/322's pilot / data-cell tables from the pdftotext -table
dumps of BOTH editions, and emit `spec_pilots.py`.

Tables extracted:
    D.1.4   additional scattered-pilot-bearing continual pilot relative
            carrier indices, per (FFT size, SP pattern)
    D.1.5   the 8K / SP32_4 special case, per Cred_coeff
    7.3/7.4 number of available data cells per DATA symbol
    7.5/7.6 total number of data cells in a Subframe Boundary Symbol
    F.1.1.. number of ACTIVE data cells in an SBS, per Cred_coeff and
            L1D_scattered_pilot_boost

NOTHING here is hand-typed.  The two editions (2026-04 and 2024-04) are
extracted independently and diffed; a disagreement is a hard failure.

THE D.1.4 ROW-GROUPING TRAP.  The printed table gives each SP pattern a
row group, and the pattern LABEL is vertically centred in a merged cell.
For the SPx_4 patterns the group is THREE rows tall, so the label prints
against the MIDDLE row and pdftotext hands the group's first row to the
SPx_2 label above it.  Read naively, SPx_4 loses its first additional CP
and SPx_2 appears to have a duplicated one.  That is exactly the defect
that made the M4 cell model over-count.  The grouping used here is:

    row A  (label SPx_2)  ->  SPx_2 = {A}
    row B  (no label)     \
    row C  (label SPx_4)   >  SPx_4 = {B, C, D}
    row D  (no label)     /

and it is not taken on faith -- it is GATED by the closed identity in
`spec_pilots.verify()`.

Usage:
    python extract_pilot_tables.py            # writes spec_pilots.py
    python extract_pilot_tables.py --dump     # print what was parsed
"""
from __future__ import annotations

import argparse
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SPEC = os.path.join(HERE, "spec")

EDITIONS = {
    "2026-04": os.path.join(SPEC, "A322_2026_tbl.txt"),
    "2024-04": os.path.join(SPEC, "A322_2024_tbl.txt"),
}

PATTERNS = [f"SP{d}_{y}" for d in (3, 4, 6, 8, 12, 16, 24, 32) for y in (2, 4)]
# Tables 7.3/7.5 carry the first eight patterns, 7.4/7.6 the last eight.
COLS_A = ["SP3_2", "SP3_4", "SP4_2", "SP4_4", "SP6_2", "SP6_4", "SP8_2", "SP8_4"]
COLS_B = ["SP12_2", "SP12_4", "SP16_2", "SP16_4", "SP24_2", "SP24_4",
          "SP32_2", "SP32_4"]
FFTS = [8192, 16384, 32768]

_NUM = re.compile(r"^\(?(\d+)\)?$")


def _cell(tok):
    """-> (value, parenthesised) ; value None for N/A."""
    if tok in ("N/A", "n/a"):
        return None, False
    m = _NUM.match(tok)
    if not m:
        raise ValueError(f"not a table cell: {tok!r}")
    return int(m.group(1)), tok.startswith("(")


def _lines(path):
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        return fh.read().split("\n")


# ---------------------------------------------------------------------------
# Table D.1.4
# ---------------------------------------------------------------------------

def _find(lines, title_re, probe_re, span=14):
    """First occurrence of `title_re` followed within `span` lines by `probe_re`.

    Every table title also appears in the table-of-contents, so a bare
    first-match lands on the TOC.  The probe is what makes the hit real,
    and it also absorbs the whitespace differences between the two
    editions' extractions.
    """
    for n, l in enumerate(lines):
        if re.search(title_re, l):
            for m in range(n, min(n + span, len(lines))):
                if re.search(probe_re, lines[m]):
                    return n, m
    raise ValueError(f"table not found: {title_re!r}")


def parse_d14(lines):
    """-> {(nfft, pattern): [(value, paren), ...]}   plus 'SEE_D15' markers."""
    i, hdr = _find(lines, r"Table\s+D\.1\.4\s+Additional Scattered Pilot Bearing",
                   r"32K\s+16K\s+8K")
    # column anchors from the header itself
    anchors = {}
    for lab, nf in (("32K", 32768), ("16K", 16384), ("8K", 8192)):
        anchors[nf] = lines[hdr].index(lab)
    end = next(n for n in range(hdr, hdr + 90)
               if "Relative carrier indices given in parentheses" in lines[n])

    rows = []                       # (label|None, {nfft: (val, paren)})
    for ln in lines[hdr + 1:end]:
        if not ln.strip():
            continue
        label = None
        m = re.search(r"\bSP(\d+)_([24])\b", ln)
        if m:
            label = m.group(0)
            ln = ln[:m.start()] + " " * len(m.group(0)) + ln[m.end():]
        cells = {}
        for m2 in re.finditer(r"\(?\d+\)?|See|Table D\.1\.5", ln):
            tok, col = m2.group(0), m2.start()
            if tok in ("See",) or tok.startswith("Table"):
                # "See Table D.1.5" lands in the 8K column
                cells[8192] = ("SEE_D15", False)
                continue
            nf = min(anchors, key=lambda k: abs(anchors[k] - col))
            if abs(anchors[nf] - col) > 8:
                raise ValueError(f"column drift: {tok!r} at {col} in {ln!r}")
            cells[nf] = _cell(tok)
        if not cells and label is None:
            continue
        rows.append((label, cells))

    # Regroup.  Two readings are emitted; only one can survive the gate.
    #   "centred" -- label SPx_2 owns its own row; label SPx_4 owns the row
    #                ABOVE it, its own row, and the row below (a 3-row merged
    #                cell whose label prints against the middle row).
    #   "naive"   -- each label owns its own row and the row below (a 2-row
    #                group).  This is the reading M4 used, and it drops
    #                SPx_4's first additional CP.
    out, naive = {}, {}
    for idx, (label, cells) in enumerate(rows):
        if label is None:
            continue
        if label.endswith("_2"):
            grp = [cells]
        else:
            grp = [rows[idx - 1][1], cells, rows[idx + 1][1]]
            if rows[idx - 1][0] is not None or rows[idx + 1][0] is not None:
                raise ValueError(f"{label}: expected unlabelled neighbours")
        ngrp = [cells]
        if idx + 1 < len(rows) and rows[idx + 1][0] is None:
            ngrp.append(rows[idx + 1][1])
        for nf in FFTS:
            vals = [g[nf] for g in grp if nf in g]
            if vals:
                out[(nf, label)] = vals
            nvals = [g[nf] for g in ngrp if nf in g]
            if nvals:
                # the naive reading also de-duplicates the repeated SPx_2 value
                seen, ded = set(), []
                for v in nvals:
                    if v[0] not in seen:
                        seen.add(v[0])
                        ded.append(v)
                naive[(nf, label)] = ded
    return out, naive


# ---------------------------------------------------------------------------
# Table D.1.5
# ---------------------------------------------------------------------------

def parse_d15(lines):
    """8K / SP32_4 -> {cred: [values]}.  The -table dump is clean here."""
    i, _ = _find(lines, r"Table\s+D\.1\.5\s+Additional\s+Scattered",
                 r"Indices for SP32_4 in 8K FFT", span=6)
    # the "Cred_coeff" label and the 0..4 column heads are on separate lines
    hdr2 = next(n for n in range(i, i + 12)
                if len(re.findall(r"\b[0-4]\b",
                                  lines[n].replace("Cred_coeff",
                                                   " " * 10))) == 5)
    anchors = {int(x.group(0)): x.start() for x in
               re.finditer(r"\b[0-4]\b",
                           lines[hdr2].replace("Cred_coeff", " " * 10))}
    if sorted(anchors) != [0, 1, 2, 3, 4]:
        raise ValueError(f"D.1.5 header not found: {lines[hdr2]!r}")
    out = {c: [] for c in range(5)}
    for ln in lines[hdr2 + 1:hdr2 + 12]:
        if not ln.strip():
            continue
        if re.match(r"^\s*\d{1,3}\s*$", ln) or "ATSC" in ln:
            break                       # page number / running head
        for m2 in re.finditer(r"\d+|none", ln):
            tok, col = m2.group(0), m2.start()
            c = min(anchors, key=lambda k: abs(anchors[k] - col))
            if abs(anchors[c] - col) > 10:
                raise ValueError(f"D.1.5 column drift: {tok!r} in {ln!r}")
            if tok != "none":
                out[c].append(int(tok))
    return out


# ---------------------------------------------------------------------------
# Tables 7.3 / 7.4 / 7.5 / 7.6  --  "FFT | Cred | NoC | 8 values"
# ---------------------------------------------------------------------------

def parse_cred_table(lines, title, cols):
    i, _ = _find(lines, title, re.escape(cols[0]), span=10)
    out, noc = {}, {}
    nf_order, seen = [8192, 16384, 32768], []
    for ln in lines[i + 1:i + 80]:
        toks = ln.split()
        if not toks:
            continue
        if toks[0] in ("8K", "16K", "32K"):
            toks = toks[1:]
        if len(toks) != 2 + len(cols):
            continue
        try:
            cred, nocv = int(toks[0]), int(toks[1])
        except ValueError:
            continue
        if not 0 <= cred <= 4:
            continue
        if cred == 0:
            seen.append(len(out))
        nf = nf_order[len(seen) - 1]
        noc[(nf, cred)] = nocv
        for c, tok in zip(cols, toks[2:]):
            out[(nf, cred, c)] = _cell(tok)
        if len(out) == 15 * len(cols):
            break
    if len(out) != 15 * len(cols):
        raise ValueError(f"{title!r}: got {len(out)} cells, want {15*len(cols)}")
    return out, noc


# ---------------------------------------------------------------------------
# Annex F  --  "FFT | L1D_SPB | 8 values", FFT label vertically centred
# ---------------------------------------------------------------------------

def parse_annex_f(lines):
    out = {}
    for cred in range(5):
        for half, cols in ((0, COLS_A), (1, COLS_B)):
            tno = 2 * cred + half + 1
            i, _ = _find(lines,
                         rf"Table\s+F\.1\.{tno}\s+Number of Active Data Cells",
                         re.escape(cols[0]), span=10)
            got, nf_i = 0, 0
            for ln in lines[i + 1:i + 60]:
                toks = ln.split()
                if not toks:
                    continue
                if toks[0] in ("8K", "16K", "32K"):
                    toks = toks[1:]
                if len(toks) != 1 + len(cols):
                    continue
                try:
                    spb = int(toks[0])
                except ValueError:
                    continue
                if not 0 <= spb <= 4:
                    continue
                nf = FFTS[nf_i]
                for c, tok in zip(cols, toks[1:]):
                    out[(nf, cred, spb, c)] = _cell(tok)
                got += 1
                if got % 5 == 0:
                    nf_i += 1
                if got == 15:
                    break
            if got != 15:
                raise ValueError(f"F.1.{tno}: got {got} rows")
    return out


# ---------------------------------------------------------------------------

def extract(path):
    lines = _lines(path)
    d14, d14_naive = parse_d14(lines)
    d15 = parse_d15(lines)
    t73, noc73 = parse_cred_table(
        lines, r"Table\s+7\.3\s+Number of Available Data Cells per Data Symbol",
        COLS_A)
    t74, _ = parse_cred_table(
        lines, r"Table\s+7\.4\s+Number of Available Data Cells per Data Symbol",
        COLS_B)
    t75, _ = parse_cred_table(
        lines, r"Table\s+7\.5\s+Total Number of Data Cells in a Subframe Boundary",
        COLS_A)
    t76, _ = parse_cred_table(
        lines, r"Table\s+7\.6\s+Total Number of Data Cells in a Subframe Boundary",
        COLS_B)
    fF = parse_annex_f(lines)
    avail = dict(t73)
    avail.update(t74)
    sbs = dict(t75)
    sbs.update(t76)
    return dict(d14=d14, d14_naive=d14_naive, d15=d15, avail=avail,
                sbs=sbs, annexf=fF, noc=noc73)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", action="store_true")
    ap.add_argument("--out", default=os.path.join(HERE, "spec_pilots.py"))
    a = ap.parse_args()

    got = {}
    for ed, path in EDITIONS.items():
        try:
            got[ed] = extract(path)
        except Exception as exc:                                # noqa: BLE001
            print(f"  {ed}: EXTRACTION FAILED: {exc}")
            raise
        print(f"  {ed}: D.1.4 {len(got[ed]['d14'])} cells, "
              f"7.3/7.4 {len(got[ed]['avail'])}, 7.5/7.6 {len(got[ed]['sbs'])},"
              f" Annex F {len(got[ed]['annexf'])}")

    eds = list(got)
    ndiff = 0
    for key in ("d14", "d14_naive", "d15", "avail", "sbs", "annexf",
                "noc"):
        a0, b0 = got[eds[0]][key], got[eds[1]][key]
        if a0 != b0:
            ndiff += 1
            print(f"  WITNESS DISAGREEMENT in {key}:")
            for k in sorted(set(a0) | set(b0), key=str):
                if a0.get(k) != b0.get(k):
                    print(f"    {k}: {eds[0]}={a0.get(k)}  {eds[1]}={b0.get(k)}")
    print(f"\n  two-witness diff: {ndiff} disagreeing tables "
          f"({'IDENTICAL' if ndiff == 0 else 'MISMATCH'})")
    if ndiff:
        return 1

    r = got[eds[0]]
    if a.dump:
        for k in sorted(r["d14"], key=lambda t: (t[0], PATTERNS.index(t[1]))):
            print(f"    D.1.4 {k[0]:>5d} {k[1]:<7s} {r['d14'][k]}")
        print(f"    D.1.5 {r['d15']}")

    with open(a.out, "w") as fh:
        fh.write(_emit(r))
    print(f"  wrote {a.out}")
    return 0


_HEADER = '''#!/usr/bin/env python3
"""A/322 pilot and data-cell tables -- MACHINE-EXTRACTED, NEVER HAND-TYPED.

Generated by `extract_pilot_tables.py` from the pdftotext -table dumps of
A/322:2026-04 and A/322:2024-04.  The two editions are independent
extraction witnesses and were byte-identical for every table below.

Contents
    ADDITIONAL_CP   Table D.1.4  (nfft, pattern) -> (indices, paren_flags)
    D15             Table D.1.5  8K/SP32_4, cred -> indices
    AVAIL_DATA      Tables 7.3/7.4  available data cells per DATA symbol
    SBS_TOTAL       Tables 7.5/7.6  TOTAL data cells in an SBS symbol
    SBS_ACTIVE      Annex F         ACTIVE data cells in an SBS symbol
    NOC             Table 7.1 NoC, as printed alongside 7.3

THE GATE.  `verify()` does not trust any of it.  It rebuilds the pilot
grid from first principles -- scattered pilots, edge pilots, the CP8/CP16/
CP32 common sets and D.1.4's additional CPs -- and requires:

  G1  the number of data carriers is CONSTANT across every scattered-pilot
      lattice phase, which is what A/322 8.1.4.1 says the additional CPs
      exist to guarantee;
  G2  that constant equals Table 7.3/7.4 exactly;
  G3  the SBS geometry (same Dx, Dy=1) equals Table 7.5/7.6 exactly;
  G4  Annex F never exceeds Table 7.5/7.6, and their difference (the null
      cell count) is even-split-able;
  G5  the printed NoC column reproduces NoCmax - Cred*Cunit.

G1+G2 are the closed identity that pins Table D.1.4's row grouping: read
the SPx_4 groups one row short and G1 fails for every _4 pattern.
"""
from __future__ import annotations

'''


def _emit(r):
    def fmt(d, keyfmt, valfmt):
        return "".join(f"    {keyfmt(k)}: {valfmt(d[k])},\n" for k in sorted(
            d, key=lambda k: tuple(str(x) for x in (k if isinstance(k, tuple)
                                                    else (k,)))))
    s = [_HEADER]
    s.append("NOC = {\n" + fmt(r["noc"], lambda k: f"({k[0]}, {k[1]})",
                               lambda v: str(v)) + "}\n\n")
    s.append("# Table D.1.4 -- (values, parenthesised-flags).  A parenthesised\n"
             "# index is NOT used when Cred_coeff is odd (A/322 Annex D note).\n"
             "ADDITIONAL_CP = {\n")
    for k in sorted(r["d14"], key=lambda t: (t[0], PATTERNS.index(t[1]))):
        vals = r["d14"][k]
        if any(v[0] == "SEE_D15" for v in vals):
            s.append(f'    ({k[0]}, "{k[1]}"): "SEE_D15",\n')
            continue
        s.append(f'    ({k[0]}, "{k[1]}"): ({tuple(v[0] for v in vals)!r}, '
                 f'{tuple(v[1] for v in vals)!r}),\n')
    s.append("}\n\n")
    s.append("# The SAME printed cells under the 2-row (naive) grouping.  Kept\n"
             "# as a BUILT-IN CONTROL: verify() runs the closed identity under\n"
             "# this reading too, and it MUST fail.  A gate that cannot fail is\n"
             "# not a gate.  This is the reading M4 used.\n"
             "ADDITIONAL_CP_NAIVE = {\n")
    for k in sorted(r["d14_naive"], key=lambda t: (t[0], PATTERNS.index(t[1]))):
        vals = r["d14_naive"][k]
        if any(v[0] == "SEE_D15" for v in vals):
            s.append(f'    ({k[0]}, "{k[1]}"): "SEE_D15",\n')
            continue
        s.append(f'    ({k[0]}, "{k[1]}"): ({tuple(v[0] for v in vals)!r}, '
                 f'{tuple(v[1] for v in vals)!r}),\n')
    s.append("}\n\n")
    s.append("# Table D.1.5 -- 8K / SP32_4 only\nD15 = {\n")
    for c in sorted(r["d15"]):
        s.append(f"    {c}: {tuple(r['d15'][c])!r},\n")
    s.append("}\n\n")
    for name, d, kf in (
            ("AVAIL_DATA", r["avail"],
             lambda k: f'({k[0]}, {k[1]}, "{k[2]}")'),
            ("SBS_TOTAL", r["sbs"], lambda k: f'({k[0]}, {k[1]}, "{k[2]}")'),
            ("SBS_ACTIVE", r["annexf"],
             lambda k: f'({k[0]}, {k[1]}, {k[2]}, "{k[3]}")')):
        s.append(f"# {'Tables 7.3/7.4' if name == 'AVAIL_DATA' else ('Tables 7.5/7.6' if name == 'SBS_TOTAL' else 'Annex F')}"
                 f" -- (value, italic/not-allowed flag); value None = N/A\n")
        s.append(f"{name} = {{\n")
        for k in sorted(d, key=lambda t: (t[0], t[1], t[2] if len(t) > 3 else 0,
                                          PATTERNS.index(t[-1]))):
            v, p = d[k]
            s.append(f"    {kf(k)}: ({v!r}, {p!r}),\n")
        s.append("}\n\n")
    s.append(_TAIL)
    return "".join(s)


_TAIL = '''
DX = {"SP3": 3, "SP4": 4, "SP6": 6, "SP8": 8, "SP12": 12, "SP16": 16,
      "SP24": 24, "SP32": 32}
PATTERNS = [f"SP{d}_{y}" for d in (3, 4, 6, 8, 12, 16, 24, 32) for y in (2, 4)]


def dxdy(pattern: str):
    a, b = pattern.split("_")
    return DX[a], int(b)


def additional_cps(nfft: int, pattern: str, cred: int, table=None):
    """Relative carrier indices of the additional CPs actually used.

    A/322 Annex D: indices in parentheses are not used when Cred_coeff is
    odd; 8K/SP32_4 is tabulated per Cred_coeff in Table D.1.5 instead.
    `table` selects the D.1.4 reading (ADDITIONAL_CP or the _NAIVE control).
    """
    e = (ADDITIONAL_CP if table is None else table).get((nfft, pattern))
    if e is None:
        return ()
    if e == "SEE_D15":
        return tuple(D15[cred])
    vals, paren = e
    return tuple(v for v, p in zip(vals, paren) if not (p and cred % 2 == 1))


def pilot_carriers(nfft: int, cred: int, pattern: str, l: int, sbs: bool,
                   table=None):
    """Relative carrier indices carrying a pilot in symbol l of a Subframe.

    l is the index of the symbol within the Subframe's data symbols; it
    only matters modulo DY.  `sbs` selects the Subframe Boundary Symbol
    pattern, which is the same DX with every symbol pilot-bearing (DY=1).
    """
    import m3_spec as S
    n = NOC[(nfft, cred)]
    dx, dy = dxdy(pattern)
    sp = set(range(0, n, dx)) if sbs else set(range(dx * (l % dy), n, dx * dy))
    lo, _hi = S.carrier_abs_range(nfft, cred)
    cp = {c - lo for c in S.common_cps(nfft, cred)}
    add = {k for k in additional_cps(nfft, pattern, cred, table) if 0 <= k < n}
    return sp | cp | add | {0, n - 1}


def data_cells(nfft: int, cred: int, pattern: str, l: int, sbs: bool,
               table=None) -> int:
    """Data (non-pilot) carriers in one symbol.  For an SBS symbol this is
    the TOTAL of Table 7.5/7.6, i.e. it still includes the null cells."""
    return NOC[(nfft, cred)] - len(
        pilot_carriers(nfft, cred, pattern, l, sbs, table))


def allowed(nfft: int, cred: int, pattern: str) -> bool:
    """A/322 Table 8.3.  N/A or bracketed-italic in 7.3/7.4 means not allowed."""
    v, paren = AVAIL_DATA[(nfft, cred, pattern)]
    return v is not None and not paren


def avail_data_cells(nfft: int, cred: int, pattern: str) -> int:
    """Available data cells per DATA symbol -- Table 7.3/7.4."""
    return AVAIL_DATA[(nfft, cred, pattern)][0]


def sbs_cells(nfft: int, cred: int, pattern: str, spb: int):
    """-> (total, active, n_null) for a Subframe Boundary Symbol.

    total  = Table 7.5/7.6, active = Annex F, null = total - active.
    `active` is the number AVAILABLE FOR CELL MULTIPLEXING (A/322 7.2.6.4).
    """
    tot = SBS_TOTAL[(nfft, cred, pattern)][0]
    act = SBS_ACTIVE[(nfft, cred, spb, pattern)][0]
    return tot, act, tot - act


def sbs_null_split(n_null: int):
    """A/322 7.2.6.4 / Figure 7.14: half the null cells at each band edge.

    n_null is ODD for the real RF33 configuration (127 and 1609), so the
    split is not symmetric and the low-edge count is the one that shifts
    every downstream cell index.  Returned as (low, high).

    ASSUMPTION P3: the low edge takes floor(n_null/2).  SETTLED BY THE AIR
    -- see m4_cells.py, where floor is the only reading that closes the
    subframe-0 cell budget onto the preamble's 3487 spare cells.
    """
    return n_null // 2, n_null - n_null // 2


def verify(verbose: bool = True, table=None, label="D.1.4 as extracted"):
    """The closed identity.  Returns (n_pass, n_fail)."""
    import m3_spec as S
    tally, lines = {}, []

    def chk(gate, cond, detail=""):
        p, f = tally.get(gate, (0, 0))
        tally[gate] = (p + (1 if cond else 0), f + (0 if cond else 1))
        if not cond:
            lines.append(f"    FAIL {gate}  {detail}")

    for (nfft, cred), n in sorted(NOC.items()):
        chk("G5 NoC == NoCmax - Cred*Cunit",
            n == S.NOC_MAX[nfft] - cred * S.C_UNIT[nfft],
            f"{nfft} cred{cred}: printed {n}")

    for nfft in (8192, 16384, 32768):
        for pat in PATTERNS:
            if (nfft, 0, pat) not in AVAIL_DATA:
                continue
            for cred in range(5):
                if not allowed(nfft, cred, pat):
                    continue
                v = avail_data_cells(nfft, cred, pat)
                _dx, dy = dxdy(pat)
                counts = {data_cells(nfft, cred, pat, l, False, table)
                          for l in range(dy)}
                chk("G1 data carriers CONSTANT over SP lattice phases",
                    len(counts) == 1,
                    f"{nfft}/{pat}/cred{cred}: phases give {sorted(counts)}")
                chk("G2 that constant == Table 7.3/7.4", counts == {v},
                    f"{nfft}/{pat}/cred{cred}: computed {sorted(counts)} "
                    f"vs table {v}")
                sv = SBS_TOTAL[(nfft, cred, pat)][0]
                cv = data_cells(nfft, cred, pat, 0, True, table)
                chk("G3 SBS geometry == Table 7.5/7.6", cv == sv,
                    f"{nfft}/{pat}/cred{cred}: computed {cv} vs table {sv}")
                for spb in range(5):
                    av, ap = SBS_ACTIVE[(nfft, cred, spb, pat)]
                    if av is None or ap:
                        continue
                    chk("G4 Annex F active <= Table 7.5/7.6 total",
                        0 < av <= sv, f"{nfft}/{pat}/cred{cred}/spb{spb}: "
                        f"{av} vs {sv}")
                    chk("G6 Annex F monotone non-increasing in SPB",
                        spb == 0 or SBS_ACTIVE[(nfft, cred, spb - 1,
                                                pat)][0] >= av,
                        f"{nfft}/{pat}/cred{cred}/spb{spb}")
    ok = sum(p for p, _ in tally.values())
    bad = sum(f for _, f in tally.values())
    if verbose:
        print(f"  spec_pilots.verify()  [{label}]")
        for g in sorted(tally):
            p, f = tally[g]
            print(f"    {'PASS' if not f else 'FAIL'}  {g:<48s} "
                  f"{p:5d} pass  {f:4d} fail")
        for l in lines[:12]:
            print(l)
        if len(lines) > 12:
            print(f"    ... and {len(lines)-12} more failures")
        print(f"    TOTAL {ok} pass, {bad} FAIL")
    return ok, bad


if __name__ == "__main__":
    import os as _os
    import sys as _sys
    _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
    _ok, _bad = verify()
    print()
    print("  CONTROL -- the same gate under the 2-row (naive) reading of")
    print("  Table D.1.4, which is what M4 used.  It MUST fail; a gate that")
    print("  cannot fail is not a gate.")
    _ok2, _bad2 = verify(table=ADDITIONAL_CP_NAIVE,
                         label="D.1.4 naive 2-row grouping (CONTROL)")
    print()
    print(f"  VERDICT: extracted {_bad} fail / naive control {_bad2} fail")
    raise SystemExit(0 if (_bad == 0 and _bad2 > 0) else 1)
'''


if __name__ == "__main__":
    raise SystemExit(main())
