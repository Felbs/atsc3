#!/usr/bin/env python3
"""Bank the parsed physical-layer numeric tables as committed artifacts.

WHY THIS EXISTS
---------------
The decoder needs A/322's non-uniform constellation (NUC) tables to build its
demapper. Those numbers live inside the published standard, and until now a
fresh clone could not decode until the reader downloaded the standard by hand
and ran an extraction pass. That was the single biggest obstacle between
"cloned it" and "watching television".

WHAT IS SHIPPED, AND WHAT IS NOT
--------------------------------
This banks the **numbers only** -- the constellation points, as a compact .npz.
A parity-check matrix and a constellation are mathematical facts; facts are not
copyrightable. The standards **documents** are copyrighted by their publishers
and are NOT redistributed by this project: no PDF, no extracted prose, no table
text or layout, ever. `spec/` and `lab/spec/` stay gitignored.

REGENERATING
------------
Only needed if you are changing what is banked; the artifact is committed.

    python lab/bank_tables.py            # writes lab/spec_bank/nuc_a322.npz
    python lab/bank_tables.py --verify   # re-parse the spec and compare

Running this requires the extracted A/322 text in lab/spec/ (i.e. a developer
box, not a fresh clone).
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BANK_DIR = os.path.join(HERE, "spec_bank")
NUC_BANK = os.path.join(BANK_DIR, "nuc_a322.npz")
SFB_BANK = os.path.join(BANK_DIR, "ac4_sfb_offsets.json")
QWIN_BANK = os.path.join(HERE, "spec_bank", "ac4_qwin.json")
NOISE_BANK = os.path.join(HERE, "spec_bank", "ac4_aspx_noise.json")
HCB_BANK = os.path.join(BANK_DIR, "ac4_huffman.json")


def _key(table: str, cr: str) -> str:
    """npz keys must be filename-safe: 'C.1.2' + '6/15' -> 'C_1_2__6_15'."""
    return f"{table.replace('.', '_')}__{cr.replace('/', '_')}"


def build(verify: bool = False) -> int:
    sys.path.insert(0, HERE)
    import spec_nuc

    if not os.path.isdir(spec_nuc.SPEC):
        print(f"error: {spec_nuc.SPEC} not found -- this needs the extracted "
              f"A/322 text, which only exists on a development box.")
        return 2

    payload: dict[str, np.ndarray] = {}
    for table, m, cols in spec_nuc.TABLES:
        path = os.path.join(spec_nuc.SPEC, "A322_2026_tbl.txt")
        vecs = spec_nuc._parse(path, table, m, len(cols))
        for cr, w in zip(cols, vecs):
            payload[_key(table, cr)] = np.asarray(w, dtype=np.complex128)

    # A banked table that is wrong is worse than no banked table, so gate the
    # payload on the same physical identity the parser gates on: a normalised
    # constellation has unit mean power.
    bad = []
    for table, m, cols in spec_nuc.TABLES:
        tol = {16: 1.5e-3, 64: 2e-4, 256: 2e-4}[m]
        for cr in cols:
            x = spec_nuc.constellation(payload[_key(table, cr)])
            p = float(np.mean(np.abs(x) ** 2))
            if abs(p - 1.0) > tol:
                bad.append(f"NUC_{m}_{cr}: mean power {p:.6f}")
    if bad:
        print("REFUSING TO BANK -- failed the unit-power gate:")
        for b in bad:
            print("   ", b)
        return 1

    if verify:
        if not os.path.exists(NUC_BANK):
            print(f"error: no bank at {NUC_BANK}")
            return 2
        have = np.load(NUC_BANK)
        diffs = [k for k in payload
                 if k not in have or not np.allclose(payload[k], have[k])]
        if diffs:
            print(f"MISMATCH in {len(diffs)} table(s): {diffs[:4]}")
            return 1
        print(f"verify: {len(payload)} banked tables match a fresh parse of the spec")
        return 0

    os.makedirs(BANK_DIR, exist_ok=True)
    np.savez_compressed(NUC_BANK, **payload)
    size = os.path.getsize(NUC_BANK)
    print(f"banked {len(payload)} NUC tables -> {NUC_BANK} ({size:,} bytes)")
    print("gate: all constellations unit mean power")
    return build_sfb(verify=verify)


def build_sfb(verify: bool = False) -> int:
    """Bank the AC-4 scale-factor-band offsets (ETSI TS 103 190 Annex B).

    Same principle as the constellations: the OFFSETS are integers -- facts --
    and they ship. The document they were read out of does not.
    """
    import json
    import m27_sfb

    if not os.path.exists(m27_sfb.DEFAULT_PDF):
        print(f"note: {os.path.basename(m27_sfb.DEFAULT_PDF)} not present -- "
              f"skipping the AC-4 band bank (needs a development box)")
        return 0

    tables = m27_sfb.split_tables()
    out: dict[str, list[int]] = {}
    for length in sorted(m27_sfb.TABLE_FOR):
        out[str(length)] = [int(v) for v in m27_sfb.offsets_for(length, tables=tables)]

    # Gate: a band layout is monotonic and ends exactly at the transform length.
    bad = [k for k, off in out.items()
           if off[-1] != int(k) or any(off[i] >= off[i + 1] for i in range(len(off) - 1))]
    if bad:
        print(f"REFUSING TO BANK -- band layout gate failed for {bad}")
        return 1

    if verify:
        if not os.path.exists(SFB_BANK):
            print(f"error: no bank at {SFB_BANK}")
            return 2
        have = json.load(open(SFB_BANK, encoding="utf-8"))
        if have != out:
            print("MISMATCH: banked AC-4 band offsets differ from a fresh parse")
            return 1
        print(f"verify: {len(out)} banked AC-4 band tables match a fresh parse")
        return 0

    with open(SFB_BANK, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
        fh.write("\n")
    print(f"banked {len(out)} AC-4 band tables -> {SFB_BANK} "
          f"({os.path.getsize(SFB_BANK):,} bytes)")
    print("gate: every layout monotonic and ends at its transform length")
    rc = build_hcb(verify=verify)
    return rc or build_floats(verify=verify)


def build_floats(verify: bool = False) -> int:
    """Bank the two AC-4 FLOAT tables the A-SPX path needs (E108).

    The codebooks and band layouts were banked on 8/12 so a fresh clone could
    decode; these two were missed, and they are what high-frequency
    regeneration needs. Without them a clone decoded core-band-only audio and
    said nothing about it. Each is gated on a property the spec itself states,
    so a bad transcription cannot be banked quietly.
    """
    import json
    import numpy as np
    import m23_hcb as H
    if not os.path.exists(H.DEFAULT_C):
        print(f"note: {os.path.basename(H.DEFAULT_C)} not present -- skipping "
              f"the AC-4 float tables (the committed bank stands)")
        return 0
    import m33_qmf as Q
    import m36_envadj as E
    qw = Q.qwin(H.DEFAULT_C)
    nt = E.noise_table(H.DEFAULT_C)
    if qw.size != 640:
        print(f"REFUSING TO BANK: QWIN has {qw.size} taps, expected 640")
        return 1
    mag = np.abs(nt)
    if nt.size != 512 or not np.allclose(mag, 1.0, atol=1e-6):
        print(f"REFUSING TO BANK: ASPX_NOISE {nt.size} entries, magnitude "
              f"{mag.min():.6f}..{mag.max():.6f} (the spec says unit)")
        return 1
    pairs = ((QWIN_BANK, {"name": "QWIN", "n": int(qw.size),
                          "coef": [float(x) for x in qw]}),
             (NOISE_BANK, {"name": "ASPX_NOISE", "n": int(nt.size),
                           "re": [float(x) for x in nt.real],
                           "im": [float(x) for x in nt.imag]}))
    for path, payload in pairs:
        if verify:
            if not os.path.exists(path):
                print(f"error: no bank at {path}")
                return 1
            have = json.load(open(path, encoding="utf-8"))
            for k in payload:
                if have.get(k) != payload[k]:
                    print(f"verify FAILED: {os.path.basename(path)} field {k}")
                    return 1
            print(f"verify: {os.path.basename(path)} matches a fresh parse")
            continue
        payload["note"] = ("Numeric facts only -- the ETSI document itself is "
                           "not ours to redistribute.")
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            json.dump(payload, fh, indent=1)
            fh.write("\n")
        print(f"banked {payload['name']} -> {path} "
              f"({os.path.getsize(path):,} bytes)")
    print("gate: QWIN 640 taps; ASPX_NOISE 512 entries of unit magnitude")
    return 0


def build_hcb(verify: bool = False) -> int:
    """Bank the AC-4 Huffman codebooks (ETSI TS 103 190).

    Without these the AC-4 decoder cannot run at all, so a clone got video and
    no sound. Same rule as the rest: the NUMBERS ship, the document does not.
    """
    import json
    import m23_hcb as H

    if not os.path.exists(H.DEFAULT_C):
        print(f"note: {os.path.basename(H.DEFAULT_C)} not present -- skipping "
              f"the AC-4 codebook bank (needs a development box)")
        return 0

    arrays = H.parse_c(H.DEFAULT_C)

    # Gate: a complete prefix code has Kraft sum exactly 1.0 and no codeword
    # that prefixes another. A transposed digit breaks both.
    bad = []
    for name in sorted(arrays):
        if not name.endswith("_LEN"):
            continue
        cw = arrays.get(name[:-4] + "_CW")
        if cw is None:
            continue
        lengths = arrays[name]
        k = H.kraft(lengths)
        if abs(k - 1.0) > 1e-9:
            bad.append(f"{name}: Kraft sum {k!r} != 1.0")
            continue
        ok, where = H.prefix_free(lengths, cw)
        if not ok:
            bad.append(f"{name}: not prefix-free at {where}")
    if bad:
        print("REFUSING TO BANK -- the codebooks failed their own gate:")
        for b in bad[:6]:
            print("   ", b)
        return 1

    if verify:
        if not os.path.exists(HCB_BANK):
            print(f"error: no bank at {HCB_BANK}")
            return 2
        have = json.load(open(HCB_BANK, encoding="utf-8"))
        if have != arrays:
            print("MISMATCH: banked AC-4 codebooks differ from a fresh parse")
            return 1
        print(f"verify: {len(arrays)} banked AC-4 codebook arrays match a fresh parse")
        return 0

    with open(HCB_BANK, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(arrays, fh, sort_keys=True, separators=(",", ":"))
        fh.write("\n")
    print(f"banked {len(arrays)} AC-4 codebook arrays -> {HCB_BANK} "
          f"({os.path.getsize(HCB_BANK):,} bytes)")
    print("gate: every codebook Kraft-complete and prefix-free")
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--verify", action="store_true",
                    help="re-parse the spec and compare against the committed bank")
    raise SystemExit(build(verify=ap.parse_args().verify))
