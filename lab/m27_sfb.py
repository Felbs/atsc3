#!/usr/bin/env python3
"""M27 -- Table B.4: the scale factor band offsets, and what max_sfb means.

M24 hard-coded `SFB_OFFSET = [0, 4, 8, 12]` because the LFE only needs four
entries.  The full-band channels report `max_sfb = 43`, so the whole column has
to be read.  This reads it, gates it, and states the number.

THIS FILE WAS WRONG ONCE.  THE POSTMORTEM IS THE POINT.
--------------------------------------------------------
The first version of this table was contaminated from sfb 40 upward and claimed
the transform has 49 scale factor bands.  It passed six gates.  It was wrong,
and it silently broke the five full-band channels in M28 while leaving the LFE
working perfectly -- which is exactly why it survived so long.

**Failure 1: the table boundary is not the page boundary.**  Table B.4's last
row (sfb 55) is the FIRST row on page 259, and Table B.5 begins immediately
below it on that same page, its own sfb column restarting at 0.  Reading "pages
258 and 259" therefore overwrote B.4's rows 40..55 with B.5's rows of the same
index.  I had already been bitten by contamination once in this extraction and
"fixed" it by narrowing the page range -- which is not a fix, because the
structure I was reading does not respect pages.  The scan now stops when the
sfb column stops increasing: the table's own structure ending, not a page
number.  (This is the same lesson as anchoring a map to the structure's origin
rather than the container's, arrived at from the other direction.)

**Failure 2: the corroborating evidence was manufactured by the bug.**  The bad
table's widths ended `32x14, 128x6`, and `max_sfb = 43` landed exactly on that
32->128 boundary.  I wrote that up as independent structural confirmation: the
encoder "stops exactly where the fine resolution stops."  It was a coincidence
produced by the contamination.  The real widths end `52x2, 64x7` and 43 lands
nowhere special.  A satisfying story that explains your result is not evidence
for it -- and the more satisfying it is, the harder it is to notice that you
never tested it.

**Failure 3: every gate was internal.**  Monotone, multiples of 4, sums to
1536, no gaps -- the contaminated table passed all of them, because it was
built from a DIFFERENT VALID TABLE and so inherited every structural property a
band layout has.  Internal consistency cannot detect a coherently wrong table.

**What actually caught it: Table B.1.**  The spec states num_sfb per transform
length in its own table, and `asf_section_data`'s NOTE points at it explicitly.
num_sfb(1536 @ 48 kHz) = 55.  My table said 49 -- which is num_sfb(1536 @ 96
kHz), a real number from a real table, just not ours.  I had derived "49" from
where my own column crossed 1536, then used it to validate that same column:
circular, and the circle closed cleanly enough that nothing complained.  Table
B.1 is the independent source, it was three pages earlier, and I had not read
it.  The first gate in gate_table() is now that check.

WHAT THE TABLE ACTUALLY IS
---------------------------
Header: 44,1/48 kHz with transform 2048/1920/1536; or 96 kHz with
4096/3840/3072; or 192 kHz with 8192/7680/6144.  Three value columns, one per
family.  And over sfb 0..55 **all three columns are identical** -- they diverge
only above 55, where the 1536 family no longer applies.  So the anxious
question of "which column is 1536?" that drove the first attempt did not even
need answering in our range; what needed answering was where the table STOPS.

Printed two table rows per line (sfb n left, sfb n+56 right), with thousands
set as `1 024` -- two word tokens.  So values are binned by x-coordinate and
concatenated within a bin; token counting reads a different column depending on
how many values on that line happen to be four digits.

THE ANSWER
-----------
    num_sfb = 55 bands, sfb_offset[55] = 1536      (Table B.1 and B.4 agree)
    max_sfb = 43  ->  sfb_offset[43] = 768 lines

    CODED BANDWIDTH  0 .. 12000 Hz  (50 % of the spectrum)
    A-SPX            12000 .. 24000 Hz

12000 Hz was also what the contaminated table gave, because sfb 43 is one of
the rows where the wrong table coincidentally agreed with the right one.  The
right answer for the wrong reason is still worth flagging: it is the number
that stood while everything around it was false.

`max_sfb` is an EXCLUSIVE end -- bands 0..42 are coded.  That is not a
convenience: M24 depends on it, decoding `SFB_OFFSET[min(max_sfb, 3)]` = 12
lines for the LFE's max_sfb 3, which is what produced sound.

WHY THE LFE NEVER NOTICED
--------------------------
The LFE uses bands 0..2, and the tables agree up to sfb 10.  So M24/M25/M26
were correct throughout and reproduce their numbers unchanged after this fix
(r = +0.2730 / -0.3039, 100.0 % below 250 Hz).  A bug that spares the component
you are testing is invisible until you build the next thing on top of it.

Usage:
    python m27_sfb.py [--pdf ../spec/ts_10319001v010301p.pdf]
"""
from __future__ import annotations

import argparse
import collections
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_PDF = os.path.join(os.path.dirname(HERE), "spec",
                           "ts_10319001v010301p.pdf")

# Table B.4, sfb 0..55, for transform length 1536 at 48 kHz.
# Extracted by coordinate binning and gated by gate_table() below.
#
# CORRECTED.  The first version of this list was contaminated from sfb 40 up
# and claimed 49 bands; see the postmortem in the module docstring.  The
# discriminator is Table B.1: num_sfb(1536 @ 48 kHz) = 55, NOT 49.
SFB_OFFSET_1536 = [
    0, 4, 8, 12, 16, 20, 24, 28, 32, 36,
    40, 44, 52, 60, 68, 76, 84, 92, 100, 108,
    116, 124, 136, 148, 160, 172, 188, 204, 220, 240,
    260, 284, 308, 336, 364, 396, 432, 468, 508, 552,
    600, 652, 704, 768, 832, 896, 960, 1024, 1088, 1152,
    1216, 1280, 1344, 1408, 1472, 1536,
]

# Table B.1: number of scale factor bands at 44,1 / 48 kHz.
NUM_SFB_B1 = {2048: 63, 1920: 61, 1536: 55, 1024: 49, 960: 49, 768: 43,
              512: 36, 480: 36, 384: 33, 256: 20, 240: 20, 192: 18,
              128: 14, 120: 14, 96: 12}

N_1536 = 1536
NUM_SFB_1536 = NUM_SFB_B1[1536]        # 55 -- from Table B.1, independently
NYQUIST = 24000.0

# x ranges of the printed columns.  The table is two rows wide; L is the left
# half (sfb n), R the right (sfb n+56).  c3 is the 1536 family.
BINS = [("sfbL", 50, 70), ("c1L", 70, 120), ("c2L", 120, 190),
        ("c3L", 190, 260),
        ("sfbR", 295, 315), ("c1R", 320, 360), ("c2R", 395, 430),
        ("c3R", 465, 500)]
PAGES = (258, 259)                      # PDF indices of Table B.4


def extract(pdf=DEFAULT_PDF, pages=PAGES):
    """-> {column: {sfb: offset}} by binning word x-coordinates.

    THE TABLE BOUNDARY IS NOT A PAGE BOUNDARY.  Table B.4's last row (sfb 55)
    is the first row on page 259, and Table B.5 starts immediately below it on
    the SAME page -- with its own sfb column restarting at 0.  Reading "pages
    258 and 259" therefore silently overwrites B.4's rows 40..55 with B.5's,
    which is exactly the bug that produced the first, wrong table.

    So rows are consumed in document order and the scan STOPS as soon as the
    left-hand sfb stops increasing.  That is the table's own structure ending,
    which is the only reliable signal available.
    """
    import fitz
    d = fitz.open(pdf)
    cols = {1: {}, 2: {}, 3: {}}
    last_sfb = -1
    for pg in pages:
        words = [(x0, round(y0, 1), t)
                 for x0, y0, x1, y1, t, *_ in d[pg].get_text("words")]
        by_y = collections.defaultdict(list)
        for x0, y, t in words:
            by_y[y].append((x0, t))
        for y in sorted(by_y):
            items = sorted(by_y[y])
            # a table row is numbers and the separating dash, nothing else
            if len(items) < 6 or not all(re.fullmatch(r"-|\d+", t)
                                         for _, t in items):
                continue
            cell = {}
            for name, lo, hi in BINS:
                # concatenate within the bin: "1" + "024" -> 1024
                tok = [t for x, t in items if lo <= x < hi and t != "-"]
                if tok:
                    cell[name] = int("".join(tok))
            sfb_l = cell.get("sfbL")
            if sfb_l is None or sfb_l <= last_sfb:
                return cols                     # a new table has started
            last_sfb = sfb_l
            for side in "LR":
                sfb = cell.get("sfb" + side)
                if sfb is None:
                    continue
                for i in (1, 2, 3):
                    v = cell.get(f"c{i}{side}")
                    if v is not None:
                        cols[i][sfb] = v
    return cols


def split_tables(pdf=DEFAULT_PDF, pages=(258, 259, 260, 261, 262)):
    """Every scale-factor-band table in Annex B, split at its own boundaries.

    Tables B.4..B.7 run back-to-back and start mid-page, so the split is done
    where the left-hand sfb column RESETS -- the structure's own boundary, not
    a page number.  -> [ {column: {sfb: offset}} ] in document order.
    """
    import fitz
    d = fitz.open(pdf)
    tables, cur, last = [], {1: {}, 2: {}, 3: {}}, -1
    for pg in pages:
        words = [(x0, round(y0, 1), t)
                 for x0, y0, x1, y1, t, *_ in d[pg].get_text("words")]
        by_y = collections.defaultdict(list)
        for x0, y, t in words:
            by_y[y].append((x0, t))
        for y in sorted(by_y):
            items = sorted(by_y[y])
            if len(items) < 6 or not all(re.fullmatch(r"-|\d+", t)
                                         for _, t in items):
                continue
            cell = {}
            for name, lo, hi in BINS:
                tok = [t for x, t in items if lo <= x < hi and t != "-"]
                if tok:
                    cell[name] = int("".join(tok))
            sfb_l = cell.get("sfbL")
            if sfb_l is None:
                continue
            if sfb_l <= last:                      # a new table begins
                tables.append(cur)
                cur = {1: {}, 2: {}, 3: {}}
            last = sfb_l
            for side in "LR":
                sfb = cell.get("sfb" + side)
                if sfb is None:
                    continue
                for i in (1, 2, 3):
                    v = cell.get(f"c{i}{side}")
                    if v is not None:
                        cur[i][sfb] = v
    tables.append(cur)
    return tables


# Table B.7 does NOT share the two-half layout of B.4..B.6.  It serves six
# transform lengths, so it prints SIX value columns side by side with no
# separator and no right-hand half -- 37 rows, sfb 0..36.  Reading it with the
# two-half bins silently mixes column 4 into the "right-hand sfb" slot, which
# is why 96 could not be found and why 192 was right only by coincidence.
B7_PAGE = 261
B7_BINS = [(50, 70), (70, 120), (120, 190), (190, 260),
           (290, 320), (380, 420), (455, 500)]


def extract_b7(pdf=DEFAULT_PDF, page=B7_PAGE):
    """-> [ {sfb: offset} ] for Table B.7's six columns, in printed order."""
    import fitz
    d = fitz.open(pdf)
    cols = [{} for _ in range(6)]
    by_y = collections.defaultdict(list)
    for x0, y0, x1, y1, t, *_ in d[page].get_text("words"):
        by_y[round(y0, 1)].append((x0, t))
    for y in sorted(by_y):
        items = sorted(by_y[y])
        if len(items) < 4 or not all(re.fullmatch(r"\d+", t)
                                     for _, t in items):
            continue
        cell = []
        for lo, hi in B7_BINS:
            tok = [t for x, t in items if lo <= x < hi]
            cell.append(int("".join(tok)) if tok else None)
        if cell[0] is None:
            continue
        for i in range(6):
            if cell[i + 1] is not None:
                cols[i][cell[0]] = cell[i + 1]
    return cols


# WHICH ANNEX B TABLE AND COLUMN SERVES EACH TRANSFORM LENGTH, from the table
# headers, at 44,1/48 kHz.  (table index into split_tables(), column index)
#
# THIS REPLACED A SEARCH, AND THE SEARCH WAS WRONG.  The first version picked
# any column satisfying `col[num_sfb] == length`, reasoning that only the right
# column could pass.  It is not unique: the band tables are NESTED, so Table
# B.4's 1536 column also has col[43] = 768 -- the same value at the same index
# that identifies B.5's 768 column.  B.4 is searched first, so every 768-line
# partial block silently got the 1536 band layout (sfb 12 -> 52 instead of 56,
# sfb 20 -> 116 instead of 132).  Sections still tiled, because section data
# does not use offsets; only the spectral decode was wrong, and only in short
# frames.
#
# Same shape of bug as the original Table B.4 contamination: a selection rule
# that looked decisive, was under-constrained, and produced a coherently wrong
# table.  The rule is kept below, demoted from selector to GATE.
TABLE_FOR = {
    2048: (0, 0), 1920: (0, 1), 1536: (0, 2),      # Table B.4
    1024: (1, 0), 960: (1, 1), 768: (1, 2),        # Table B.5
    512: (2, 0), 480: (2, 1), 384: (2, 2),         # Table B.6
    256: (3, 0), 240: (3, 1), 192: (3, 2),         # Table B.7, six columns
    128: (3, 3), 120: (3, 4), 96: (3, 5),
}


BANK = os.path.join(HERE, "spec_bank", "ac4_sfb_offsets.json")
_BANK_CACHE = None


def _banked_offsets(length):
    """Band offsets from the committed numeric artifact, or None if absent.

    Gated on load, not trusted on faith: a layout must be strictly increasing
    and end exactly at its transform length, which is what a transposed digit
    would break.
    """
    global _BANK_CACHE
    if _BANK_CACHE is None:
        if not os.path.exists(BANK):
            return None
        import json
        with open(BANK, encoding="utf-8") as fh:
            _BANK_CACHE = json.load(fh)
    off = _BANK_CACHE.get(str(length))
    if off is None:
        return None
    if off[-1] != length or any(off[i] >= off[i + 1] for i in range(len(off) - 1)):
        raise LookupError(f"banked band layout for {length} failed its gate")
    return list(off)


def offsets_for(length, tables=None):
    """Scale factor band offsets for a transform length at 48 kHz.

    The table and column come from TABLE_FOR (the printed headers).  Table
    B.1's independent num_sfb is then used to CHECK the result rather than to
    find it.
    """
    n_sfb = NUM_SFB_B1[length]
    where = TABLE_FOR.get(length)
    if where is None:
        raise LookupError(f"no Annex B table registered for length {length}")
    # A fresh clone has no ETSI document -- it is copyrighted and not ours to
    # redistribute. The band OFFSETS are integers, so they ship as a banked
    # artifact and the audio decodes out of the box. The document wins when it
    # is present, because parsing it is the primary evidence.
    if tables is None and not os.path.exists(DEFAULT_PDF):
        off = _banked_offsets(length)
        if off is not None:
            return off
    ti, ci = where
    if ti == 3:                                   # Table B.7, six columns
        col = extract_b7()[ci]
    else:
        col = (tables if tables is not None else split_tables())[ti][ci + 1]
    out = [col.get(k) for k in range(n_sfb + 1)]
    if None in out:
        raise LookupError(f"table B.{ti + 4} column {ci} is short for "
                          f"{length} (needs sfb 0..{n_sfb})")
    # the old rule, demoted from SELECTOR to GATE -- see the note above
    if out[-1] != length or not all(out[k] < out[k + 1] for k in range(n_sfb)):
        raise LookupError(f"table B.{ti + 4} column {ci} does not end at "
                          f"{length} on sfb {n_sfb}")
    return out


def gate_table(off, n=N_1536, verbose=True):
    """Properties a correct band layout has.  -> True if all hold.

    The FIRST check is the one that matters and the one the earlier version of
    this file did not have: the band count must equal Table B.1's num_sfb.
    Every other check here is internal consistency, and the contaminated table
    passed all of them -- monotone, multiples of 4, summing to 1536 -- while
    being wrong from sfb 11 up.  Internal consistency cannot catch a table that
    is coherently wrong; only an independent source can.
    """
    w = [off[k + 1] - off[k] for k in range(len(off) - 1)]
    checks = [
        (f"band count == Table B.1 num_sfb ({NUM_SFB_1536})  <-- INDEPENDENT",
         len(off) - 1 == NUM_SFB_1536),
        ("all bands present, none silently dropped",
         len(off) == NUM_SFB_1536 + 1),
        ("strictly increasing", all(off[k] < off[k + 1]
                                    for k in range(len(off) - 1))),
        (f"last band ends at exactly {n}", off[-1] == n),
        ("widths non-decreasing", all(w[i] <= w[i + 1]
                                      for i in range(len(w) - 1))),
        ("every width a multiple of 4", all(v % 4 == 0 for v in w)),
        (f"widths sum to {n}", sum(w) == n),
    ]
    for label, ok in checks:
        if verbose:
            print(f"    {'PASS' if ok else 'FAIL'}  {label}")
    return all(ok for _, ok in checks)


def bandwidth(max_sfb, off=None, n=N_1536):
    """Coded bandwidth in Hz for a given max_sfb."""
    off = off or SFB_OFFSET_1536
    return off[max_sfb] / n * NYQUIST


def main(argv=None):
    ap = argparse.ArgumentParser(prog="m27 sfb")
    ap.add_argument("--pdf", default=DEFAULT_PDF)
    ap.add_argument("--max-sfb", type=int, default=43)
    a = ap.parse_args(argv)

    print("M27 -- Table B.4 and the coded bandwidth")
    print("=" * 72)

    print("\n  1. the banked table")
    ok = gate_table(SFB_OFFSET_1536)

    print("\n  2. re-extract from the PDF and compare")
    if os.path.exists(a.pdf):
        try:
            cols = extract(a.pdf)
            live = [cols[3][k] for k in range(NUM_SFB_1536 + 1)]
            same = live == SFB_OFFSET_1536
            print(f"    {'PASS' if same else 'FAIL'}  extraction reproduces "
                  f"the banked table exactly")
            # the column choice, decided rather than assumed
            # the three columns are IDENTICAL over sfb 0..55 and diverge
            # only above it, where the 1536 family no longer applies
            same = all(cols[1][k] == cols[2][k] == cols[3][k]
                       for k in range(NUM_SFB_1536 + 1))
            print(f"    {'PASS' if same else 'FAIL'}  all three columns agree "
                  f"over sfb 0..{NUM_SFB_1536} (they diverge only above it)")
            for i in (1, 2, 3):
                hit = next((k for k in sorted(cols[i])
                            if cols[i][k] == N_1536), None)
                print(f"      column {i}: reaches {N_1536} at sfb {hit}"
                      + (f"   <-- Table B.1 says {NUM_SFB_1536}"
                         if hit == NUM_SFB_1536 else ""))
            ok = ok and same
        except ImportError:
            print("    (PyMuPDF not available -- banked table only)")
    else:
        print(f"    (no PDF at {a.pdf} -- banked table only)")

    print("\n  3. band structure")
    w = [SFB_OFFSET_1536[k + 1] - SFB_OFFSET_1536[k] for k in range(49)]
    runs = []
    for v in w:
        if runs and runs[-1][0] == v:
            runs[-1][1] += 1
        else:
            runs.append([v, 1])
    print("    widths  " + "  ".join(f"{v}x{c}" for v, c in runs))
    # `fine` is the index of the LAST fine band; `max_sfb` is an exclusive end,
    # so the comparison below is against fine + 1.  (I got this backwards on the
    # first run and printed "unexpected" for a result that was exactly right.)
    fine = next(i for i in range(len(w) - 1, -1, -1) if w[i] <= 32)
    print(f"    bands 0..{fine} are the fine-resolution region; the last "
          f"{len(w) - fine - 1} bands are 128 lines wide")

    print("\n  4. what our stream's max_sfb means")
    m = a.max_sfb
    bw = bandwidth(m)
    print(f"    max_sfb = {m}   ->  sfb_offset[{m}] = {SFB_OFFSET_1536[m]} "
          f"of {N_1536} lines")
    print(f"    coded    0 .. {bw:.0f} Hz      "
          f"({100 * SFB_OFFSET_1536[m] / N_1536:.0f} % of the spectrum)")
    print(f"    A-SPX    {bw:.0f} .. {NYQUIST:.0f} Hz  "
          f"({N_1536 - SFB_OFFSET_1536[m]} lines, not entropy coded)")
    edge = m == fine + 1
    print(f"    {'PASS' if edge else 'note'}  max_sfb is an exclusive end, so "
          f"bands 0..{m - 1} are coded -- exactly the fine region 0..{fine}"
          if edge else
          f"    note  max_sfb {m} does not land on the fine/coarse boundary "
          f"({fine + 1})")

    print("\n" + "=" * 72)
    print(f"  CODED BANDWIDTH = {bw:.0f} Hz." if ok else
          "  *** table did not gate; the bandwidth above is not trustworthy ***")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
