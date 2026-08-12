#!/usr/bin/env python3
"""Machine-extract A/322's DATA-PLP bit-interleaver tables, both editions.

`spec_bicm.py` carries only the L1 versions of these (Tables 6.20-6.24).
The data-PLP chain (A/322 6.2) is a different set entirely:

    6.2.1  parity interleaver      -- Type B codes only, Qldpc from Table 6.7
    6.2.2  group-wise interleaver  -- permutations in ANNEX B (B.1 for
                                      Ninner=64800/Ngroup=180, B.2 for
                                      Ninner=16200/Ngroup=45)
    6.2.3  block interleaver       -- Type A (Table 6.10) or Type B
                                      (Table 6.11), selected by Table 6.8
                                      (64800) / Table 6.9 (16200)
    6.3.3  demultiplexing          -- Table 6.14

A CORRECTION TO THE M4 DESIGN NOTE.  It listed "group-wise interleaver
(Tables 6.10/6.11), block interleaver (Tables 6.12/6.13)".  That mapping
is wrong on both counts.  6.10/6.11 are the BLOCK interleaver
configurations; the group-wise permutations are not in section 6 at all,
they are the Annex B table series.  And 6.12/6.13 are the MANDATORY
MODULATION AND CODING COMBINATIONS tables -- a conformance checklist of
check-marks, not decoder data, and the check-mark glyphs do not survive
pdftotext at all.  They are extracted here only to record that they carry
nothing a receiver needs.

Usage:
    python extract_bitint_tables.py            # writes spec_bitint.py
    python extract_bitint_tables.py --dump
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

RATES = [f"{n}/15" for n in range(2, 14)]
MOD_BITS = {"QPSK": 2, "16QAM": 4, "64QAM": 6, "256QAM": 8,
            "1024QAM": 10, "4096QAM": 12}
_SKIP = re.compile(r"ATSC|Physical Layer|Annex [A-Z]|Order\s+of\s+Group|"
                   r"Ninner|PERMUTATION|SEQUENCES", re.I)


def _lines(path):
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        return fh.read().split("\n")


# ---------------------------------------------------------------------------
# Annex B -- the group-wise interleaving permutations
# ---------------------------------------------------------------------------

def parse_annex_b(lines):
    """-> {(ninner, modulation, rate): tuple(perm)}  and the identity headers."""
    # two printed variants: "(Ninner = 64800)" and "(Code length = 16200 bits)"
    title = re.compile(r"Table\s+B\.([12])\.(\d+)\s+([0-9A-Za-z]+)\s*"
                       r"\(\s*(?:Ninner|Code\s+length)\s*=\s*(\d+)")
    hits = [(n, m) for n, l in enumerate(lines) for m in [title.search(l)] if m]
    # every table title also appears in the table-of-contents; a real one is
    # followed within a few lines by an actual data row
    hits = [(n, m) for n, m in hits
            if any(len(re.findall(r"\b\d+\b", lines[j])) >= 15
                   for j in range(n + 1, min(n + 12, len(lines))))]
    if not hits:
        raise ValueError("Annex B tables not found")
    out, headers, labelmap = {}, {}, {}
    for k, (n, m) in enumerate(hits):
        mod, ninner = m.group(3), int(m.group(4))
        ngroup = ninner // 360
        end = hits[k + 1][0] if k + 1 < len(hits) else len(lines)
        # the next Annex begins with its own "Annex X:" banner
        for j in range(n + 1, end):
            if re.match(r"\s*Annex [C-Z]:", lines[j]):
                end = j
                break
        toks, labels = [], []
        for ln in lines[n + 1:end]:
            if _SKIP.search(ln):
                continue
            lab = re.findall(r"\b(\d+)\s*/\s*15\b", ln)
            ln = re.sub(r"\b\d+\s*/\s*15\b", " ", ln)      # code-rate labels
            ln = re.sub(r"\b(Code|Rate|MOD|CR)\b", " ", ln)
            nums = re.findall(r"\b\d+\b", ln)
            if len(nums) < 5:            # page numbers / stray captions
                continue
            # record which BLOCK each printed rate label fell in.  The label
            # is vertically centred in its block, so its position is taken
            # from the running token count, which pins the block ORDER
            # (ASSUMPTION B1) instead of assuming it.
            for x in lab:
                labels.append((f"{x}/15", len(toks) // ngroup))
            toks.extend(int(x) for x in nums)
        if len(toks) != ngroup * (len(RATES) + 1):
            raise ValueError(f"B.{m.group(1)}.{m.group(2)} {mod}: got "
                             f"{len(toks)} values, want "
                             f"{ngroup*(len(RATES)+1)}")
        blocks = [toks[i * ngroup:(i + 1) * ngroup]
                  for i in range(len(RATES) + 1)]
        headers[(ninner, mod)] = tuple(blocks[0])
        labelmap[(ninner, mod)] = tuple(labels)
        for r, b in zip(RATES, blocks[1:]):
            out[(ninner, mod, r)] = tuple(b)
    return out, headers, labelmap


# ---------------------------------------------------------------------------
# Tables 6.8 / 6.9  --  block interleaver TYPE
# ---------------------------------------------------------------------------

def parse_bi_type(lines, table, ninner, nmod):
    i = next(n for n, l in enumerate(lines)
             if re.search(rf"Table\s+{table}\s+Block\s+Interleaver\s+Type", l)
             and any("2/15" in lines[m] for m in range(n, n + 8)))
    out, got = {}, 0
    for ln in lines[i + 1:i + 60]:
        toks = ln.split()
        if len(toks) != 1 + len(RATES):
            continue
        try:
            mod = int(toks[0])
        except ValueError:
            continue
        if mod not in (2, 4, 6, 8, 10, 12) or set(toks[1:]) - {"A", "B"}:
            continue
        for r, t in zip(RATES, toks[1:]):
            out[(ninner, mod, r)] = t
        got += 1
        if got == nmod:
            break
    if got != nmod:
        raise ValueError(f"Table {table}: got {got} rows, want {nmod}")
    return out


# ---------------------------------------------------------------------------
# Table 6.10 (Type A) / 6.11 (Type B) / 6.14 (demux)
# ---------------------------------------------------------------------------

def _rowsafe(tok):
    return None if tok in ("N/A", "n/a") else int(tok)


def parse_610(lines):
    i = next(n for n, l in enumerate(lines)
             if re.search(r"Table\s+6\.10\s+Type\s+A\s+Block\s+Interleaver", l)
             and any("QPSK" in lines[m] for m in range(n, n + 8)))
    out = {}
    for ln in lines[i + 1:i + 40]:
        toks = ln.split()
        if len(toks) != 6 or toks[0] not in MOD_BITS:
            continue
        mod = toks[0]
        out[(64800, mod)] = (_rowsafe(toks[1]), _rowsafe(toks[3]),
                             int(toks[5]))
        out[(16200, mod)] = (_rowsafe(toks[2]), _rowsafe(toks[4]),
                             int(toks[5]))
        if len(out) == 12:
            break
    if len(out) != 12:
        raise ValueError(f"Table 6.10: got {len(out)} entries")
    return out


def parse_611(lines):
    i = next(n for n, l in enumerate(lines)
             if re.search(r"Table\s+6\.11\s+Parameters\s+for\s+Type\s+B", l)
             and any("QPSK" in lines[m] for m in range(n, n + 10)))
    out = {}
    for ln in lines[i + 1:i + 40]:
        toks = ln.split()
        if len(toks) != 6 or toks[0] not in MOD_BITS:
            continue
        mod, nq = toks[0], int(toks[1])
        out[(64800, mod)] = (nq, _rowsafe(toks[2]), _rowsafe(toks[4]))
        out[(16200, mod)] = (nq, _rowsafe(toks[3]), _rowsafe(toks[5]))
        if len(out) == 12:
            break
    if len(out) != 12:
        raise ValueError(f"Table 6.11: got {len(out)} entries")
    return out


def parse_614(lines):
    i = next(n for n, l in enumerate(lines)
             if re.search(r"Table\s+6\.14\s+Parameters\s+for\s+Bit-Mapping", l)
             and any("QPSK" in lines[m] for m in range(n, n + 10)))
    out = {}
    for ln in lines[i + 1:i + 40]:
        toks = ln.split()
        if len(toks) != 4 or toks[0] not in MOD_BITS:
            continue
        out[(64800, toks[0])] = (int(toks[1]), _rowsafe(toks[2]))
        out[(16200, toks[0])] = (int(toks[1]), _rowsafe(toks[3]))
        if len(out) == 12:
            break
    if len(out) != 12:
        raise ValueError(f"Table 6.14: got {len(out)} entries")
    return out


def parse_612_613(lines):
    """Tables 6.12/6.13 -- mandatory combinations.  The check-mark glyphs do
    not survive pdftotext, so this records WHAT SURVIVED, honestly."""
    out = {}
    for table, ninner in (("6.12", 64800), ("6.13", 16200)):
        i = next(n for n, l in enumerate(lines)
                 if re.search(rf"Table\s+{table}\s+Mandatory\s+Modulation", l)
                 and any("QPSK" in lines[m] for m in range(n, n + 10)))
        rows = []
        for ln in lines[i + 1:i + 30]:
            toks = ln.split()
            if not toks or toks[0] not in MOD_BITS:
                continue
            rows.append((toks[0], len(toks) - 1))
        out[ninner] = rows
    return out


def extract(path):
    lines = _lines(path)
    b, hdr, lab = parse_annex_b(lines)
    t = dict(parse_bi_type(lines, r"6\.8", 64800, 6))
    t.update(parse_bi_type(lines, r"6\.9", 16200, 4))
    return dict(annexb=b, annexb_hdr=hdr, annexb_lab=lab, bitype=t, t610=parse_610(lines),
                t611=parse_611(lines), t614=parse_614(lines),
                t612=parse_612_613(lines))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", action="store_true")
    ap.add_argument("--out", default=os.path.join(HERE, "spec_bitint.py"))
    a = ap.parse_args()

    got = {}
    for ed, path in EDITIONS.items():
        got[ed] = extract(path)
        print(f"  {ed}: Annex B {len(got[ed]['annexb'])} permutations, "
              f"6.8/6.9 {len(got[ed]['bitype'])}, 6.10 {len(got[ed]['t610'])}, "
              f"6.11 {len(got[ed]['t611'])}, 6.14 {len(got[ed]['t614'])}")

    eds = list(got)
    ndiff = 0
    for key in ("annexb", "annexb_hdr", "annexb_lab", "bitype", "t610", "t611", "t614"):
        a0, b0 = got[eds[0]][key], got[eds[1]][key]
        if a0 != b0:
            ndiff += 1
            print(f"  WITNESS DISAGREEMENT in {key}:")
            for k in sorted(set(a0) | set(b0), key=str):
                if a0.get(k) != b0.get(k):
                    print(f"    {k}: differs")
    print(f"\n  two-witness diff: {ndiff} disagreeing tables "
          f"({'IDENTICAL' if ndiff == 0 else 'MISMATCH'})")
    if ndiff:
        return 1

    r = got[eds[0]]
    if a.dump:
        for k in sorted(r["t610"]):
            print(f"    6.10 {k} Nr1/Nr2/Nc = {r['t610'][k]}")
        for k in sorted(r["t611"]):
            print(f"    6.11 {k} NQCB_IG/Npart1/Npart2 = {r['t611'][k]}")
        for k in sorted(r["annexb"])[:3]:
            print(f"    B    {k} -> {r['annexb'][k][:12]}...")
        print(f"    6.12/6.13 (check-mark tables): {r['t612']}")

    with open(a.out, "w") as fh:
        fh.write(_emit(r))
    print(f"  wrote {a.out}")
    return 0


_HEADER = '''#!/usr/bin/env python3
"""A/322 DATA-PLP bit-interleaver tables -- MACHINE-EXTRACTED, NEVER TYPED.

Generated by `extract_bitint_tables.py` from the pdftotext -table dumps of
A/322:2026-04 and A/322:2024-04.  The two editions are independent
extraction witnesses and were identical for every table below.

This is the DATA-PLP chain of A/322 6.2.  `spec_bicm.py` carries the L1
chain (Tables 6.20-6.24); the two share no tables.

    GROUPWISE   Annex B  group-wise interleaving permutation, keyed
                (Ninner, modulation, code_rate) -> tuple of Ngroup indices
    BI_TYPE     Tables 6.8/6.9  "A" or "B" block interleaver
    TYPE_A      Table 6.10  (Nr1, Nr2, Nc)
    TYPE_B      Table 6.11  (NQCB_IG, Npart1, Npart2)
    DEMUX       Table 6.14  (MOD, output data cells)

THE GATES (`verify()`), all closed arithmetic, none of them "the PDF says
so":

  H1  every Annex B sequence is a PERMUTATION of 0..Ngroup-1 -- exact,
      and a single mis-extracted digit breaks it
  H2  the identity header row above each Annex B block reads 0..Ngroup-1
      in order (it is printed as data and extracted the same way, so it
      is a free per-table checksum of the column ordering)
  H3  Table 6.10:  Nr1 + Nr2 == Ninner / Nc, Nc == MOD, and
                   Nr1 == floor((Ninner/Nc)/360) * 360
  H4  Table 6.11:  Npart1 + Npart2 == Ninner, NQCB_IG == MOD, and
                   Npart1 == floor(Ninner/(360*NQCB_IG)) * 360 * NQCB_IG
  H5  CROSS-TABLE: Nr1 * Nc == Npart1 for every (Ninner, modulation).
      Tables 6.10 and 6.11 are printed pages apart and describe different
      interleavers; that they pin each other is the strongest check here.
  H6  Table 6.14 output cells == Ninner / MOD
  H7  Tables 6.8/6.9 contain only "A"/"B" over all 12 code rates

  H8  the printed code-rate labels land in blocks 1..12 IN ORDER.  The
      labels have to be stripped before parsing (they tokenise as "2" and
      "15" and would corrupt the sequence), so their positions are
      recorded on the way out and checked afterwards.  This is what pins
      the block ORDER; without it H1 would pass just as happily with two
      rate blocks swapped.

ASSUMPTIONS
  B1  (was: Annex B block order.)  SETTLED BY H8 -- the labels are read
      back positionally in all 10 tables rather than assumed.  What
      remains unpinned is only whether a label printed against the middle
      row of its block could belong to a NEIGHBOURING block; that would
      require the label to sit >=4 rows from its own block, which the
      printed layout does not do.  Falsifiable off the air regardless:
      a swapped block cannot make PLP 0's LDPC converge.
  B2  Tables 6.12/6.13 carry no extractable content (check-mark glyphs).
      They are a conformance checklist, not decoder data, so nothing
      downstream depends on them.  RECORDED, NOT WORKED AROUND.
"""
from __future__ import annotations

'''


def _emit(r):
    s = [_HEADER]
    s.append("RATES = " + repr(RATES) + "\n")
    s.append("MOD_BITS = " + repr(MOD_BITS) + "\n\n")
    s.append("# Table 6.10 -- Type A block interleaver: (Nr1, Nr2, Nc)\n"
             "TYPE_A = {\n")
    for k in sorted(r["t610"]):
        s.append(f'    ({k[0]}, "{k[1]}"): {r["t610"][k]!r},\n')
    s.append("}\n\n")
    s.append("# Table 6.11 -- Type B block interleaver: "
             "(NQCB_IG, Npart1, Npart2)\nTYPE_B = {\n")
    for k in sorted(r["t611"]):
        s.append(f'    ({k[0]}, "{k[1]}"): {r["t611"][k]!r},\n')
    s.append("}\n\n")
    s.append("# Table 6.14 -- (MOD, output data cells per FEC frame)\n"
             "DEMUX = {\n")
    for k in sorted(r["t614"]):
        s.append(f'    ({k[0]}, "{k[1]}"): {r["t614"][k]!r},\n')
    s.append("}\n\n")
    s.append("# Tables 6.8 / 6.9 -- block interleaver type\nBI_TYPE = {\n")
    for k in sorted(r["bitype"]):
        s.append(f'    ({k[0]}, {k[1]}, "{k[2]}"): "{r["bitype"][k]}",\n')
    s.append("}\n\n")
    s.append("# Tables 6.12/6.13 -- mandatory modulation/coding combinations.\n"
             "# The check-mark glyphs do NOT survive pdftotext; what remains is\n"
             "# only the row labels.  Recorded so nobody re-extracts them\n"
             "# expecting decoder data -- they are a conformance checklist and\n"
             "# carry nothing a receiver needs.\n"
             f"MANDATORY_COMBINATIONS_UNREADABLE = {r['t612']!r}\n\n")
    s.append("# Annex B -- group-wise interleaving permutation Pi(j)\n"
             "GROUPWISE = {\n")
    for k in sorted(r["annexb"], key=lambda t: (t[0], MOD_BITS[t[1]],
                                                RATES.index(t[2]))):
        s.append(f'    ({k[0]}, "{k[1]}", "{k[2]}"): {r["annexb"][k]!r},\n')
    s.append("}\n\n")
    s.append("# the identity header row printed above each Annex B block\n"
             "GROUPWISE_HDR = {\n")
    for k in sorted(r["annexb_hdr"], key=lambda t: (t[0], MOD_BITS[t[1]])):
        s.append(f'    ({k[0]}, "{k[1]}"): {r["annexb_hdr"][k]!r},\n')
    s.append("}\n\n")
    s.append("# where each printed code-rate label fell, as (label, block\n"
             "# index).  Block 0 is the identity header, so label i must sit\n"
             "# in block i+1.  This is what pins the block ORDER (B1).\n"
             "GROUPWISE_LABELS = {\n")
    for k in sorted(r["annexb_lab"], key=lambda t: (t[0], MOD_BITS[t[1]])):
        s.append(f'    ({k[0]}, "{k[1]}"): {r["annexb_lab"][k]!r},\n')
    s.append("}\n\n")
    s.append(_TAIL)
    return "".join(s)


_TAIL = '''
def verify(verbose: bool = True):
    """The closed identities.  Returns (n_pass, n_fail)."""
    tally, msgs = {}, []

    def chk(gate, cond, detail=""):
        p, f = tally.get(gate, (0, 0))
        tally[gate] = (p + (1 if cond else 0), f + (0 if cond else 1))
        if not cond:
            msgs.append(f"    FAIL {gate}  {detail}")

    for (ninner, mod, rate), perm in GROUPWISE.items():
        ng = ninner // 360
        chk("H1 Annex B sequence is a permutation of 0..Ngroup-1",
            sorted(perm) == list(range(ng)), f"{ninner}/{mod}/{rate}")
    for (ninner, mod), hdr in GROUPWISE_HDR.items():
        ng = ninner // 360
        chk("H2 Annex B identity header reads 0..Ngroup-1",
            list(hdr) == list(range(ng)), f"{ninner}/{mod}")

    for (ninner, mod), (nr1, nr2, nc) in TYPE_A.items():
        if nr1 is None:
            continue
        chk("H3 Table 6.10 Nr1+Nr2 == Ninner/Nc and Nc == MOD",
            nc == MOD_BITS[mod] and nr1 + nr2 == ninner // nc
            and nr1 == (ninner // nc) // 360 * 360,
            f"{ninner}/{mod}: {nr1}+{nr2} vs {ninner//nc}, Nc {nc}")
    for (ninner, mod), (nq, np1, np2) in TYPE_B.items():
        if np1 is None:
            continue
        chk("H4 Table 6.11 Npart1+Npart2 == Ninner and NQCB_IG == MOD",
            nq == MOD_BITS[mod] and np1 + np2 == ninner
            and np1 == ninner // (360 * nq) * (360 * nq),
            f"{ninner}/{mod}: {np1}+{np2} vs {ninner}, NQCB_IG {nq}")
    for k, (nr1, _nr2, nc) in TYPE_A.items():
        if nr1 is None:
            continue
        chk("H5 CROSS-TABLE Nr1*Nc == Npart1 (6.10 pins 6.11)",
            nr1 * nc == TYPE_B[k][1], f"{k}: {nr1*nc} vs {TYPE_B[k][1]}")
    for (ninner, mod), (m, cells) in DEMUX.items():
        if cells is None:
            continue
        chk("H6 Table 6.14 cells == Ninner/MOD",
            m == MOD_BITS[mod] and cells == ninner // m, f"{ninner}/{mod}")
    for k, v in BI_TYPE.items():
        chk("H7 Tables 6.8/6.9 entries are A or B", v in ("A", "B"), str(k))
    for k, labs in GROUPWISE_LABELS.items():
        chk("H8 printed rate labels land in blocks 1..12, in order",
            [x[0] for x in labs] == RATES
            and [x[1] for x in labs] == list(range(1, len(RATES) + 1)),
            f"{k}: {labs}")

    ok = sum(p for p, _ in tally.values())
    bad = sum(f for _, f in tally.values())
    if verbose:
        print("  spec_bitint.verify()")
        for g in sorted(tally):
            p, f = tally[g]
            print(f"    {'PASS' if not f else 'FAIL'}  {g:<52s} "
                  f"{p:5d} pass  {f:4d} fail")
        for m in msgs[:20]:
            print(m)
        print(f"    TOTAL {ok} pass, {bad} FAIL")
    return ok, bad


def plan(ninner: int, mod: str, rate: str):
    """The bit-interleaver configuration for one (Ninner, modulation, rate)."""
    t = BI_TYPE[(ninner, MOD_BITS[mod], rate)]
    return dict(ninner=ninner, modulation=mod, rate=rate, bi_type=t,
                mod_bits=MOD_BITS[mod],
                groupwise=GROUPWISE[(ninner, mod, rate)],
                type_a=TYPE_A[(ninner, mod)], type_b=TYPE_B[(ninner, mod)],
                out_cells=DEMUX[(ninner, mod)][1])


if __name__ == "__main__":
    _ok, _bad = verify()
    print()
    p = plan(16200, "64QAM", "11/15")
    print("  RF33 PLP 0 -- Ninner 16200, 64QAM, 11/15 (from L1-Detail):")
    print(f"    block interleaver type   {p['bi_type']}")
    print(f"    Table 6.10 (Nr1,Nr2,Nc)  {p['type_a']}")
    print(f"    Table 6.11 (NQ,Np1,Np2)  {p['type_b']}")
    print(f"    output data cells        {p['out_cells']}")
    print(f"    group-wise perm (45)     {p['groupwise'][:10]} ...")
    raise SystemExit(0 if _bad == 0 else 1)
'''


if __name__ == "__main__":
    raise SystemExit(main())
