#!/usr/bin/env python3
"""GATE: does a LOST MPU leave a hole, or does everything after it slide early?

A damaged MPU that is partly usable already holds the clock true -- the
`repair_unusable` path advances `base_time` across it. An MPU that never
arrived at all did not, so the next fragment was placed 2.002 s early and
every fragment after it inherited the shift.

That is cumulative and per-lane, which is why E20's lanes ended 37 fragments
apart over two hours. It is not drift; it is a sum of un-accounted holes, and
the lanes that lose the most slide the furthest.

The test is the only one that matters: where does fragment N sit on the
timeline, with and without a loss before it?  Those must be the SAME number --
a lost MPU should cost a freeze, not a shift.

`m7_play.retime` is stubbed to a constant 180180 ticks (E19 measured that as
the broadcaster's exact per-MPU total) so the bookkeeping under test is
isolated from box parsing.

Run:  python lab/gate_gap_clock.py
"""
from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import m7_play as PL             # noqa: E402
import m11_stream as S           # noqa: E402

DUR = 180180                     # ticks per MPU, from the broadcaster's trun
PID = 12


SEEN = []


def stub_retime(body, frag, base_time):
    # Record the base_time `retime` was ACTUALLY handed. Reading `st` before
    # calling `_segment` measures the wrong instant -- the gap advance happens
    # inside, ahead of the retime -- and made this gate report a shift on the
    # first fragment after the hole and none on any later one, which is not a
    # shape the bug can even produce. Measure where the fragment LANDED.
    SEEN.append(base_time)
    return b"seg", DUR


def run(drop=()):
    """Feed MPUs 0..19, skipping `drop`. Return {seq: base_time it landed at}."""
    PL.retime = stub_retime
    del SEEN[:]
    tr = S.Transport.__new__(S.Transport)          # no radio, no threads
    tr.assets = {}
    tr.stats = __import__("collections").Counter()
    tr.want = (b"vide",)
    tr.force_pid = None

    placed = {}
    for i in range(20):
        if i in drop:
            continue
        m = dict(pid=PID, seq=1000 + i, body=b"", n_samples=120,
                 meta=b"vide-init" if not tr.assets else None)
        n = len(SEEN)
        seg = S.Transport._segment(tr, m)
        if seg is not None and len(SEEN) > n:
            placed[i] = SEEN[n]          # where retime actually put it
    return placed, tr.stats


def main():
    # moov_handler is only consulted on the first MPU of an asset; make it
    # answer 'video' without needing a real moov box.
    S.moov_handler = lambda meta: b"vide"

    clean, _ = run()
    lossy, stats = run(drop=(5,))

    print(f"  clean run  : {len(clean)} fragments placed")
    print(f"  lost MPU 5 : {len(lossy)} fragments placed, "
          f"gap {stats['gap']}, gap_mpus {stats['gap_mpus']}, "
          f"gap_ticks {stats['gap_ticks']}")
    print()
    print("  seq   clean base_time   after-loss base_time   shift")
    bad = 0
    for i in sorted(lossy):
        if i not in clean:
            continue
        d = lossy[i] - clean[i]
        if i in (4, 6, 10, 19):
            print(f"  {i:3d}   {clean[i]:14d}   {lossy[i]:20d}   {d:+7d}")
        if d != 0:
            bad += 1

    print()
    if bad == 0:
        print("GATE PASS -- the lost MPU left a hole; nothing after it moved")
        return 0
    print(f"GATE FAIL -- {bad} fragments slid earlier by {DUR} ticks "
          f"({DUR / 90000:.3f} s) per lost MPU. Loss becomes DESYNC, and it "
          f"accumulates for the rest of the run.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
