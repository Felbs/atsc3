#!/usr/bin/env python3
"""M18 -- find AC-4's substream size field by SEARCH, not by recall.

M17 parsed the TOC as far as the container could referee it: `fs_index` and
`frame_rate_index` read out of every frame and matched what `dac4` declares
independently.  Past that point the next fields are the presentation and
substream tables, and my knowledge of their exact syntax is recall, not
reference.  Recall is where this project has been wrong before, and the
correction has always come from the air rather than from thinking harder.

So this does not parse forward from memory.  It SEARCHES for a field with a
property no coincidence should have:

    somewhere in the header there is a value that PREDICTS THE FRAME SIZE.

Every AC-4 frame carries its substreams' sizes, and the container tells us each
frame's true length from the `trun`.  So for a candidate (bit offset, width),
read the value on thousands of frames and ask whether it tracks the real length
with a CONSTANT offset.  A field that satisfies `value + k == length` on 3480
frames of varying size is the size field; nothing else will do that by accident,
because the lengths vary by hundreds of bytes.

This is the same move as M13's frequency-de-interleaver sweep: replace a
question I cannot answer from memory with an experiment whose oracle is the
data itself.

WHAT A NEGATIVE RESULT WOULD MEAN
----------------------------------
If no (offset, width) tracks the length, that is informative too: it means the
sizes are variable-length coded (AC-4 uses `variable_bits()` in places), and
the next step is a different search shape rather than a decoder.  Reported
honestly either way.

Usage:
    python m18_ac4_substream.py [--frames 400] [--max-bits 260]
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m17_ac4_walk as W                                         # noqa: E402


def bit_matrix(frames, nbits):
    """(nframes, nbits) uint8 of the first `nbits` bits of each frame."""
    n = len(frames)
    out = np.zeros((n, nbits), np.uint8)
    need = (nbits + 7) // 8
    for i, f in enumerate(frames):
        b = np.frombuffer(f[:need].ljust(need, b"\0"), np.uint8)
        out[i] = np.unpackbits(b)[:nbits]
    return out


def values_at(bits, off, width):
    """Big-endian unsigned value of bits[:, off:off+width], vectorised."""
    seg = bits[:, off:off + width].astype(np.int64)
    w = (1 << np.arange(width - 1, -1, -1)).astype(np.int64)
    return seg @ w


def search(frames, lengths, max_bits=260, widths=range(8, 21), verbose=True):
    """-> [(offset, width, k, how many frames satisfy value + k == length)]"""
    bits = bit_matrix(frames, max_bits + max(widths) + 1)
    n = len(frames)
    hits = []
    for width in widths:
        for off in range(0, max_bits):
            v = values_at(bits, off, width)
            d = lengths - v
            # a size field has a CONSTANT relationship to the frame length
            k, cnt = np.unique(d, return_counts=True)
            best = int(cnt.max())
            if best >= int(0.98 * n):
                hits.append((off, width, int(k[cnt.argmax()]), best))
    hits.sort(key=lambda h: (-h[3], h[1]))
    if verbose:
        if not hits:
            print("    no (offset, width) tracks the frame length")
        for off, width, k, cnt in hits[:8]:
            print(f"    bit {off:3d}  width {width:2d}  value + {k:5d} == "
                  f"length on {cnt}/{n} frames "
                  f"({'ALL' if cnt == n else f'{cnt/n*100:.1f}%'})")
    return hits


def main(argv=None):
    ap = argparse.ArgumentParser(prog="m18 ac4 substream")
    ap.add_argument("path", nargs="?",
                    default="m7_out/rf33_audio_pid13.mp4")
    ap.add_argument("--frames", type=int, default=400,
                    help="frames to search over before verifying on all")
    ap.add_argument("--max-bits", type=int, default=260)
    a = ap.parse_args(argv)
    p = a.path if os.path.isabs(a.path) else os.path.join(HERE, a.path)

    print("M18 -- find the AC-4 size field by search")
    print("=" * 72)
    fr = W.samples(p)
    L = np.array([len(f) for f in fr], np.int64)
    print(f"  {len(fr)} frames, length min {L.min()} max {L.max()} "
          f"mean {L.mean():.0f}, {len(np.unique(L))} distinct lengths")
    if len(np.unique(L)) < 20:
        print("  WARNING: too few distinct lengths for this search to "
              "discriminate")

    sub = fr[:a.frames]
    print(f"\n  searching the first {a.max_bits} bit positions, widths 8..20, "
          f"over {len(sub)} frames")
    hits = search(sub, L[:len(sub)], a.max_bits)

    if not hits:
        print("\n  NEGATIVE -- no fixed-position field predicts the length.")
        print("  That is a real answer: the sizes are almost certainly")
        print("  variable-length coded, so the next search is over VLC shapes,")
        print("  not over more bit offsets.  No decoder work should start on a")
        print("  guessed layout.")
        return 1

    off, width, k, _ = hits[0]
    v = values_at(bit_matrix(fr, a.max_bits + 21), off, width)
    exact = int(np.sum(v + k == L))
    print(f"\n  VERIFY the best candidate on ALL {len(fr)} frames")
    print(f"    bit {off}, width {width}, value + {k} == length on "
          f"{exact}/{len(fr)} ({exact/len(fr)*100:.2f}%)")
    ok = exact >= int(0.999 * len(fr))
    print("\n" + "=" * 72)
    if ok:
        print(f"  PASS -- the frame carries its own size at bit {off}, width "
              f"{width}.")
        print("  That is the anchor the substream table hangs off, found from")
        print("  the data rather than from memory.")
    else:
        print("  INCONCLUSIVE on the full set -- the candidate did not hold.")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
