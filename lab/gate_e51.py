#!/usr/bin/env python3
"""gate_e51.py -- fresh-start anchoring on a lane WITH truncated history.

A worker started mid-stream with --start-behind must anchor its wav
exactly at the first non-skipped fragment, pad only the holes it actually
consumes, and stamp first_seq from the fragment index. The NEGATIVE
CONTROL runs the shipped arithmetic (trim = skip0*60, pad baseline =
whole-lane shortfall, anchor slot = cursor//60) on the same lane and must
mis-anchor -- proving the gate can fail and the old code does.

    python lab/gate_e51.py    -> PASS/FAIL, exit 0 iff all pass
"""
from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tools"))

from atsc3_audio import fresh_anchor                            # noqa: E402

FAILS = []


def check(name, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" +
          (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


# a lane shaped like the 8/08 SDR-overflow afternoon: healthy head and
# tail, a diseased stretch of truncated fragments in the SKIPPED region,
# and a few truncated fragments in the CONSUMED region
fcounts = [60] * 40 + [30] * 40 + [60] * 8 + [50, 50, 50, 50] + [60] * 8
#          0..39      40..79      80..87     88..91 (consumed)  92..99
skip0 = 88                       # --start-behind lands at fragment 88
seqs = list(range(7000, 7100))   # idx: seq per fragment position


def frag_of_frame(frame):
    """Which fragment holds `frame`, walking the true counts."""
    tot = 0
    for j, c in enumerate(fcounts):
        if frame < tot + c:
            return j
        tot += c
    return len(fcounts)


def gate_new():
    print("gate 1: E51 anchoring on the truncated-history lane")
    trim, k, pad_base = fresh_anchor(fcounts, skip0)
    check("trim equals the frames the skipped fragments declare",
          trim == sum(fcounts[:88]) == 40 * 60 + 40 * 30 + 8 * 60,
          f"trim {trim}")
    check("wav sample 0 lands exactly on the anchor fragment",
          frag_of_frame(trim) == k == 88, f"frag {frag_of_frame(trim)}")
    check("first_seq stamped from fragment index matches the content",
          seqs[k] == 7088)
    check("pad baseline covers ONLY the skipped region's shortfall",
          pad_base == 40 * 30, f"{pad_base}")
    # post-anchor: consumed fragments 88..91 are 10 frames short each --
    # with the new baseline their 40 frames DO get padded
    total_short = sum(60 - c for c in fcounts if c < 60)
    check("consumed-region holes remain eligible for padding",
          total_short - pad_base == 40, f"{total_short - pad_base}")


def gate_old_negative():
    print("gate 2: NEGATIVE CONTROL -- the shipped arithmetic mis-anchors")
    trim_old = skip0 * 60                       # frames assumed 60/frag
    landed = frag_of_frame(trim_old)
    err_frags = landed - skip0
    check("old trim overshoots into later fragments (gate CAN fail)",
          err_frags > 0 and trim_old - sum(fcounts[:88]) == 1200,
          f"anchor claimed frag 88, content starts frag {landed} "
          f"(+{err_frags} frags = {err_frags * 2.002:.0f}s stale shift)")
    # old anchor slot: seqs[cursor//60] claims fragment 88's slot while the
    # trimmed content starts at `landed` (here: past the WHOLE lane -- the
    # old trim discards more frames than the lane holds)
    claimed = seqs[min(trim_old // 60, len(seqs) - 1)]
    actual = seqs[landed] if landed < len(seqs) else None
    check("old first_seq claim disagrees with actual content",
          actual != claimed, f"claimed {claimed}, actual {actual}")
    # old pad baseline = whole-lane shortfall at first pass: swallows the
    # consumed region's 40 real hole-frames -> they are never padded
    pad_base_old = sum(60 - c for c in fcounts if c < 60)
    check("old baseline un-pads the consumed region's real holes",
          pad_base_old - (40 * 30) == 40 and pad_base_old > 40 * 30,
          f"old baseline {pad_base_old} vs true skipped {40 * 30}")


if __name__ == "__main__":
    gate_new()
    gate_old_negative()
    print(f"\n{'ALL GATES PASS' if not FAILS else 'FAILED: ' + ', '.join(FAILS)}")
    sys.exit(1 if FAILS else 0)
