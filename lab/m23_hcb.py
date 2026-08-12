#!/usr/bin/env python3
"""M23 -- the AC-4 Huffman codebooks, parsed from the spec's own C file and
gated individually on the Kraft equality.

Annex A does not print the tables; it ships them as `ts_103190_tables.c` inside
`ts_10319001v010301p0.zip`.  That is a better source than a PDF in every way --
machine readable, explicit `_LEN` / `_CW` pairs, no OCR ambiguity -- but "better
source" is not "verified", so every table gets a test it can fail.

THE GATE: KRAFT EQUALITY
-------------------------
A binary prefix code with codeword lengths L satisfies

    sum(2^-L) == 1

exactly, when the code is COMPLETE -- which every Huffman code is.  Less than 1
means the table is missing codewords; more than 1 means it is over-subscribed
and cannot be a prefix code at all.  Either way the sum is a single number that
catches a truncated read, an off-by-one in the array length, or a mis-parsed
entry, and it does so per table rather than in aggregate.

A SECOND, INDEPENDENT CHECK
----------------------------
The codewords themselves must actually form a prefix code: no codeword may be a
prefix of another.  Kraft is necessary but not sufficient (a set of lengths can
satisfy it while the assigned codewords collide), so both are run.  And the
declared `codebook_length` from Annex A -- 289 for HCB_11, 81 for HCB_1, and so
on -- must match the array size, which ties the C file back to the PDF.

Usage:
    python m23_hcb.py [--c spec/ts_103190_tables.c]
"""
from __future__ import annotations

import argparse
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_C = os.path.join(os.path.dirname(HERE), "spec", "ts_103190_tables.c")

# codebook_length values printed in Annex A, for the cross-check
ANNEX_A = {
    "ASF_HCB_SCALEFAC": 121, "ASF_HCB_1": 81, "ASF_HCB_2": 81,
    "ASF_HCB_3": 81, "ASF_HCB_4": 81, "ASF_HCB_5": 81, "ASF_HCB_6": 81,
    "ASF_HCB_7": 64, "ASF_HCB_8": 64, "ASF_HCB_9": 169, "ASF_HCB_10": 169,
    "ASF_HCB_11": 289, "ASF_HCB_SNF": 22,
}

ARRAY = re.compile(
    r"const\s+\w+\s+([A-Za-z0-9_]+)\s*\[\s*(\d*)\s*\]\s*=\s*\{(.*?)\}\s*;",
    re.S)


BANK = os.path.join(HERE, "spec_bank", "ac4_huffman.json")


def parse_c(path):
    """-> {name: [ints]} for every const array in the file.

    A fresh clone has no ETSI source file -- the document is copyrighted and
    not ours to redistribute. The codebooks themselves are integers, so they
    ship as a banked artifact and the AC-4 decoder works out of the box.
    Without this the receiver decoded video and no sound.
    """
    if not os.path.exists(path) and os.path.exists(BANK):
        import json
        with open(BANK, encoding="utf-8") as fh:
            arrays = json.load(fh)
        # Gated on load, not trusted: a complete prefix code has Kraft sum
        # 1.0, which a transposed digit breaks.
        for name, lengths in arrays.items():
            if name.endswith("_LEN") and abs(kraft(lengths) - 1.0) > 1e-9:
                raise SystemExit(f"banked codebook {name} failed its Kraft gate")
        return arrays
    src = open(path, encoding="utf-8", errors="replace").read()
    src = re.sub(r"/\*.*?\*/", "", src, flags=re.S)
    out = {}
    for m in ARRAY.finditer(src):
        name, _decl, body = m.group(1), m.group(2), m.group(3)
        # base 0 chokes on the file's zero-padded decimals ("-004" reads as
        # octal and raises), so hex and decimal are separated explicitly
        vals = []
        for v in re.findall(r"-?(?:0[xX][0-9A-Fa-f]+|\d+)", body):
            vals.append(int(v, 16) if v.lower().lstrip("-").startswith("0x")
                        else int(v, 10))
        out[name] = vals
    return out


def kraft(lengths):
    """sum(2^-L).  Exactly 1.0 for a complete prefix code."""
    return sum(2.0 ** -L for L in lengths if L > 0)


def prefix_free(lengths, words):
    """No codeword may be a prefix of another.  -> (ok, first collision)."""
    seen = {}
    for i, (L, w) in enumerate(zip(lengths, words)):
        if L <= 0:
            continue
        bits = format(w & ((1 << L) - 1), f"0{L}b")
        for k in range(1, len(bits) + 1):
            p = bits[:k]
            if p in seen and (k == len(bits) or seen[p][1] == k):
                return False, (i, bits, seen[p][0])
        seen[bits] = (i, len(bits))
    # a real prefix test: walk every codeword against the set of prefixes
    words_set = set(seen)
    for bits in list(words_set):
        for k in range(1, len(bits)):
            if bits[:k] in words_set:
                return False, (seen[bits][0], bits, seen[bits[:k]][0])
    return True, None


def main(argv=None):
    ap = argparse.ArgumentParser(prog="m23 hcb")
    ap.add_argument("--c", default=DEFAULT_C)
    ap.add_argument("--show", default=None, help="dump one codebook")
    a = ap.parse_args(argv)
    if not os.path.exists(a.c):
        print(f"  no tables file at {a.c}")
        return 2
    arrays = parse_c(a.c)
    print("M23 -- AC-4 Huffman codebooks, Kraft-gated")
    print("=" * 78)
    print(f"  {len(arrays)} const arrays in {os.path.basename(a.c)}")

    pairs = sorted({n[:-4] for n in arrays if n.endswith("_LEN")}
                   & {n[:-3] for n in arrays if n.endswith("_CW")})
    print(f"  {len(pairs)} codebooks (a _LEN and a _CW each)\n")
    bad = 0
    print(f"  {'codebook':32s} {'n':>5s} {'Kraft':>12s}  prefix  annex")
    for nm in pairs:
        L, W = arrays[nm + "_LEN"], arrays[nm + "_CW"]
        k = kraft(L)
        pf, _ = prefix_free(L, W)
        ann = ANNEX_A.get(nm)
        ann_ok = "-" if ann is None else ("ok" if ann == len(L) else "MISMATCH")
        ok = abs(k - 1.0) < 1e-9 and pf and len(L) == len(W) and ann_ok != \
            "MISMATCH"
        bad += not ok
        print(f"  {'PASS' if ok else 'FAIL'} {nm:27s} {len(L):5d} "
              f"{k:12.9f}  {'yes' if pf else 'NO ':6s} {ann_ok}")
    print("\n" + "=" * 78)
    if a.show and a.show + "_LEN" in arrays:
        L, W = arrays[a.show + "_LEN"], arrays[a.show + "_CW"]
        print(f"  {a.show}: first 12 (len, codeword)")
        for i in range(min(12, len(L))):
            print(f"    [{i:3d}] len {L[i]:2d}  cw {W[i]:#07x} = "
                  f"{format(W[i] & ((1 << L[i]) - 1), f'0{L[i]}b')}")
    print(f"  {len(pairs) - bad}/{len(pairs)} codebooks pass Kraft AND the "
          f"prefix test AND the Annex A length")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
