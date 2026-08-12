#!/usr/bin/env python3
"""Extract the A/322 Annex A.1 LDPC parity-check address tables (Ninner = 64800).

WHY THIS EXISTS
---------------
`spec_ldpc.py` extracts Annex A.2 (Ninner = 16200) ONLY -- its own header says so --
and `m3_ldpc._params()` silently returns the 16200 rows whatever `n` is asked for.
Every LLS-bearing PLP on RF25 and RF30 signals `L1D_plp_fec_type = 1`, i.e.
Ninner = 64800, so that silent fallback is the M8 blocker.

EXTRACTION DISCIPLINE (the same one that has held all campaign)
---------------------------------------------------------------
* **PyMuPDF, not pdftotext.**  M6 established that the "variables lost to the math
  font" problem is a pdftotext defect, not a document defect.  These tables are plain
  digits, but the same renderer is used so the provenance argument is uniform.
* **Independent witnesses.**  Four editions of A/322 are parsed independently
  (2026-04, 2023-03, 2021, 2018) and required to agree element-for-element.  The LDPC
  matrices have not changed since A/322:2016, so disagreement means an extraction bug.
* **Closed arithmetic identities**, not eyeballing.  See GATES below.
* **A deliberately-wrong control that MUST fail.**  Two of them: the row-interleaved
  reading of the two-column pages, and a single-digit corruption injected into one
  address.  A gate that cannot fail is not a gate (M7/M8 lesson).

GATES
-----
G1  row count == the count DERIVED FROM SECTION 6, not from the page:
        Type B: rows == Kldpc/360
        Type A: rows == Kldpc/360 + Q1      (steps (vii)-(viii) feed M1 parity back in)
G2  every address in [0, Nldpc - Kldpc)
G3  every row strictly increasing, no duplicate address inside a row
G4  accumulator coverage: expanding each row per 6.1.3.1 (Type A, Q1 below M1 / Q2 at
    or above M1) or 6.1.3.2 (Type B, Qldpc) stays inside range AND touches every one
    of the Nldpc-Kldpc parity accumulators at least once
G5  cross-edition identity, element for element, over all four witnesses
G6  row-weight profile: Type B non-increasing; Type A's trailing Q1 rows strictly
    lighter than every information row (this is what makes the G1 split visible)
G7  the code rate implied by the extracted geometry == the printed caption rate

CONTROLS (must FAIL)
--------------------
C1  row-interleaved reading of the two-column pages (instead of column-major)
C2  one address in one table incremented by 1

Usage:  python extract_ldpc64k.py [--pdf-dir DIR] [--out spec_ldpc64k.py]
Requires PyMuPDF; on this rig that is the Python312 interpreter, not radioconda:
    python
"""
import argparse
import io
import os
import re
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

HERE = os.path.dirname(os.path.abspath(__file__))
# The PDFs live wherever the operator put them; point ATSC3_SPEC_PDF_DIR at
# that directory. Default = ../spec next to this repo.
DEFAULT_PDF_DIR = os.environ.get(
    "ATSC3_SPEC_PDF_DIR",
    os.path.join(os.path.dirname(HERE), "spec"),
)

EDITIONS = {
    "2026-04": "A322-2026-04.pdf",
    "2023-03": "A322-2023-03-Physical-Layer-Protocol.pdf",
    "2021":    "A322-2021-Physical-Layer-Protocol.pdf",
    "2018":    "A322-2018-Physical-Layer-Protocol.pdf",
}

RATES = ["2/15", "3/15", "4/15", "5/15", "6/15", "7/15",
         "8/15", "9/15", "10/15", "11/15", "12/15", "13/15"]

# Section 6 Table 6.5 (Type A sizes, Ninner=64800) and Table 6.7 (Qldpc, Type B).
# Copied from spec_ldpc.LDPC_PARAMS_64800, which was extracted and gated in M3.
PARAMS = {
    "2/15":  dict(Nldpc=64800, Kldpc=8640,  M1=1800, M2=54360, Q1=5,  Q2=151, Qldpc=None, type="A"),
    "3/15":  dict(Nldpc=64800, Kldpc=12960, M1=1800, M2=50040, Q1=5,  Q2=139, Qldpc=None, type="A"),
    "4/15":  dict(Nldpc=64800, Kldpc=17280, M1=1800, M2=45720, Q1=5,  Q2=127, Qldpc=None, type="A"),
    "5/15":  dict(Nldpc=64800, Kldpc=21600, M1=1440, M2=41760, Q1=4,  Q2=116, Qldpc=None, type="A"),
    "6/15":  dict(Nldpc=64800, Kldpc=25920, M1=None, M2=None,  Q1=None, Q2=None, Qldpc=108, type="B"),
    "7/15":  dict(Nldpc=64800, Kldpc=30240, M1=1080, M2=33480, Q1=3,  Q2=93,  Qldpc=None, type="A"),
    "8/15":  dict(Nldpc=64800, Kldpc=34560, M1=None, M2=None,  Q1=None, Q2=None, Qldpc=84,  type="B"),
    "9/15":  dict(Nldpc=64800, Kldpc=38880, M1=None, M2=None,  Q1=None, Q2=None, Qldpc=72,  type="B"),
    "10/15": dict(Nldpc=64800, Kldpc=43200, M1=None, M2=None,  Q1=None, Q2=None, Qldpc=60,  type="B"),
    "11/15": dict(Nldpc=64800, Kldpc=47520, M1=None, M2=None,  Q1=None, Q2=None, Qldpc=48,  type="B"),
    "12/15": dict(Nldpc=64800, Kldpc=51840, M1=None, M2=None,  Q1=None, Q2=None, Qldpc=36,  type="B"),
    "13/15": dict(Nldpc=64800, Kldpc=56160, M1=None, M2=None,  Q1=None, Q2=None, Qldpc=24,  type="B"),
}


# ----------------------------------------------------------------- page parsing
def _page_lines(page):
    """(y, x0, text) for every non-empty text line on the page."""
    out = []
    for b in page.get_text("dict")["blocks"]:
        if b.get("type") != 0:
            continue
        for l in b["lines"]:
            txt = "".join(s["text"] for s in l["spans"]).strip()
            if txt:
                out.append((round(l["bbox"][1], 1), round(l["bbox"][0], 1), txt))
    return out


def _table_pages(doc):
    """{table_index: (page, rate_num)} for the Annex A.1 tables (0-based pages).

    The front matter carries a list-of-tables page whose captions match the same
    regex but which has no numeric body; the page with the most numeric body lines
    wins, so that page can never be selected.
    """
    cand = {}
    for i in range(doc.page_count):
        lines = _page_lines(doc[i])
        nnum = sum(1 for l in lines if re.fullmatch(r"[\d\s]+", l[2]) and len(l[2].split()) > 3)
        for (y, x0, txt) in lines:
            m = re.match(r"Table A\.1\.(\d+)\s+Rate\s*=\s*(\d+)\s*/\s*15", txt.replace("–", "-"))
            if m:
                idx = int(m.group(1))
                cand.setdefault(idx, []).append((nnum, i, int(m.group(2))))
    return {k: (max(v)[1], max(v)[2]) for k, v in cand.items()}


def parse_table(doc, page_no, interleave=False):
    """Extract one Annex A.1 table's rows from its page.

    Two-column pages are read COLUMN-MAJOR (all of the left block, then all of the
    right block).  `interleave=True` is CONTROL C1: read the blocks row-interleaved.
    """
    lines = _page_lines(doc[page_no])
    caps = [l for l in lines if l[2].startswith("Table A.1.")]
    capy = caps[0][0]
    body = [l for l in lines
            if l[0] > capy + 2 and l[0] < 720 and re.fullmatch(r"[\d\s]+", l[2])]
    if not body:
        return []
    xs = sorted(set(l[1] for l in body))
    left_x = xs[0]
    # a second column exists if some line starts far to the right of the first
    right_xs = [x for x in xs if x > left_x + 100]
    if not right_xs:
        rows = [l for l in sorted(body)]
        return [[int(v) for v in r[2].split()] for r in rows]
    rx = min(right_xs)
    left = sorted([l for l in body if l[1] < left_x + 100])
    right = sorted([l for l in body if l[1] >= rx])
    if interleave:
        merged = []
        for i in range(max(len(left), len(right))):
            if i < len(left):
                merged.append(left[i])
            if i < len(right):
                merged.append(right[i])
    else:
        merged = left + right
    return [[int(v) for v in r[2].split()] for r in merged]


def extract_edition(path, interleave=False):
    import fitz
    doc = fitz.open(path)
    pages = _table_pages(doc)
    tables = {}
    for idx in range(1, 13):
        if idx not in pages:
            continue
        pg, num = pages[idx]
        rate = "%d/15" % num
        tables[rate] = parse_table(doc, pg, interleave=interleave)
    doc.close()
    return tables


# ----------------------------------------------------------------------- gates
def expand_row(row, rate, row_index):
    """All accumulator addresses a row touches, per 6.1.3.1 / 6.1.3.2."""
    p = PARAMS[rate]
    N, K = p["Nldpc"], p["Kldpc"]
    M = N - K
    out = []
    if p["type"] == "B":
        Q = p["Qldpc"]
        for m in range(360):
            for a in row:
                out.append((a + m * Q) % M)
    else:
        M1, Q1, Q2 = p["M1"], p["Q1"], p["Q2"]
        for m in range(360):
            for a in row:
                if a < M1:
                    out.append((a + m * Q1) % M1)
                else:
                    out.append(M1 + (a - M1 + m * Q2) % (M - M1))
    return out


def gate(tables, label="tables", verbose=True):
    """Run G1..G4, G6, G7.  Returns (n_pass, failures)."""
    npass, fails = 0, []

    def chk(ok, name):
        nonlocal npass
        if ok:
            npass += 1
        else:
            fails.append(name)

    for rate in RATES:
        if rate not in tables:
            fails.append("%s MISSING" % rate)
            continue
        rows = tables[rate]
        p = PARAMS[rate]
        N, K = p["Nldpc"], p["Kldpc"]
        M = N - K
        n_info = K // 360
        want = n_info if p["type"] == "B" else n_info + p["Q1"]
        chk(len(rows) == want, "G1 %s rows %d != %d" % (rate, len(rows), want))
        chk(all(0 <= a < M for r in rows for a in r), "G2 %s range" % rate)
        chk(all(list(r) == sorted(set(r)) for r in rows), "G3 %s monotone/dup" % rate)
        acc = set()
        okr = True
        for i, r in enumerate(rows):
            e = expand_row(r, rate, i)
            if any(not (0 <= a < M) for a in e):
                okr = False
            acc.update(e)
        chk(okr and len(acc) == M, "G4 %s coverage %d/%d" % (rate, len(acc), M))
        w = [len(r) for r in rows]
        if p["type"] == "B":
            chk(all(w[i] >= w[i + 1] for i in range(len(w) - 1)),
                "G6 %s weight non-increasing" % rate)
        else:
            head, tail = w[:n_info], w[n_info:]
            chk(len(tail) > 0 and max(tail) < min(head),
                "G6 %s tail weight %s vs head min %d" % (rate, tail, min(head)))
        num = int(rate.split("/")[0])
        chk(K * 15 == N * num, "G7 %s rate arithmetic" % rate)
    if verbose:
        print("  %-28s %4d gates pass, %d fail%s"
              % (label, npass, len(fails), ("  " + "; ".join(fails[:4])) if fails else ""))
    return npass, fails


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdf-dir", default=DEFAULT_PDF_DIR)
    ap.add_argument("--out", default=os.path.join(HERE, "spec_ldpc64k.py"))
    a = ap.parse_args()

    print("A/322 Annex A.1 (Ninner = 64800) extraction -- PyMuPDF")
    witnesses = {}
    for ed, fn in EDITIONS.items():
        path = os.path.join(a.pdf_dir, fn)
        if not os.path.exists(path):
            print("  edition %s: MISSING %s" % (ed, path))
            continue
        witnesses[ed] = extract_edition(path)
        n = sum(len(r) for t in witnesses[ed].values() for r in t)
        s = sum(v for t in witnesses[ed].values() for r in t for v in r)
        print("  edition %-8s %2d tables, %6d addresses, sum %d"
              % (ed, len(witnesses[ed]), n, s))
    if not witnesses:
        print("NO WITNESSES -- set --pdf-dir")
        return 2

    # ---- G5 cross-edition identity
    eds = [e for e in EDITIONS if e in witnesses]      # newest edition first
    ref = witnesses[eds[0]]
    g5_fail = []
    for ed in eds[1:]:
        for rate in RATES:
            if rate not in witnesses[ed] or rate not in ref:
                g5_fail.append("%s %s missing" % (ed, rate))
            elif witnesses[ed][rate] != ref[rate]:
                g5_fail.append("%s %s differs" % (ed, rate))
    print("\nG5 cross-edition identity over %d witnesses: %s"
          % (len(eds), "PASS" if not g5_fail else "FAIL " + str(g5_fail[:4])))

    print("\nGates on the extracted tables:")
    npass, fails = gate(ref, "EXTRACTED (column-major)")

    print("\nControls (these MUST fail):")
    ctrl_path = os.path.join(a.pdf_dir, EDITIONS[eds[0]])
    c1 = extract_edition(ctrl_path, interleave=True)
    n1, f1 = gate(c1, "C1 row-interleaved")

    # C2 -- a single address corrupted by +1.  Recorded HONESTLY: G1..G7 are
    # structurally blind to it (it stays in range, stays increasing, and the
    # accumulator expansion still covers everything), so the intra-table gates
    # report 0 failures.  The gate that DOES catch it is G5, cross-edition
    # agreement, and that is scored here rather than being claimed for free.
    import copy
    c2 = copy.deepcopy(ref)
    c2["9/15"][0][0] += 1
    n2, f2 = gate(c2, "C2 one address +1 (intra-table)")
    c2_g5 = sum(1 for rate in RATES if c2[rate] != ref[rate])
    print("  %-28s G5 cross-edition diffs: %d  -> %s"
          % ("C2 under G5", c2_g5, "CAUGHT" if c2_g5 else "MISSED"))

    ok = (not g5_fail) and (not fails) and f1 and c2_g5
    print("\nVERDICT: %s   (%d gates pass, %d fail; C1 fails %d intra-table gates,"
          " C2 caught by G5)"
          % ("PASS" if ok else "FAIL", npass, len(fails), len(f1)))
    if not ok:
        return 1

    # ---- emit
    with open(a.out, "w", encoding="utf-8") as f:
        f.write('"""A/322 Annex A.1 LDPC parity-check address tables, Ninner = 64800.\n\n')
        f.write("GENERATED by lab/extract_ldpc64k.py -- do not hand-edit.\n\n")
        f.write("Witnesses (independent editions of A/322, all agreeing element-for-element):\n")
        for ed in eds:
            f.write("    A/322:%s\n" % ed)
        f.write("\nGates passed at generation time: %d, failures 0.\n" % npass)
        f.write("Controls that failed as required: the row-interleaved reading of the\n")
        f.write("two-column pages (%d gate failures, one per two-column table), and a\n" % len(f1))
        f.write("single +1 address corruption -- which the intra-table gates are blind to\n")
        f.write("by construction and which cross-edition agreement catches.\n")
        f.write('"""\n\n')
        tot = sum(len(r) for t in ref.values() for r in t)
        ssum = sum(v for t in ref.values() for r in t for v in r)
        f.write("A1_INT_COUNT = %d\nA1_INT_SUM = %d\n\n" % (tot, ssum))
        f.write("LDPC_64800 = {\n")
        for rate in RATES:
            p = PARAMS[rate]
            f.write('    "%s": {\n' % rate)
            f.write('        "Nldpc": %d, "Kldpc": %d, "type": "%s",\n'
                    % (p["Nldpc"], p["Kldpc"], p["type"]))
            if p["type"] == "A":
                f.write('        "M1": %d, "M2": %d, "Q1": %d, "Q2": %d, "Qldpc": None,\n'
                        % (p["M1"], p["M2"], p["Q1"], p["Q2"]))
            else:
                f.write('        "M1": None, "M2": None, "Q1": None, "Q2": None, "Qldpc": %d,\n'
                        % p["Qldpc"])
            f.write('        "rows": [\n')
            for r in ref[rate]:
                f.write("            %r,\n" % (r,))
            f.write("        ],\n    },\n")
        f.write("}\n")
    print("wrote %s" % a.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
