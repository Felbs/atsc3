#!/usr/bin/env python3
"""e68_diurnal.py -- what fraction of the day is RF33 actually watchable?

THE QUESTION: did the splitter's +10.7 dB flatten the night floor? If RF33
now decodes THROUGH the night, "stable video indefinitely" stops being
antenna-blocked and the campaign changes shape.

The evidence is already on disk: data/e31/chain.log runs from 8/09 13:36
through today, so it spans BOTH ERAS -- pre-splitter and post-splitter --
recorded by the same receiver, same antenna, same settings. That makes
this a controlled before/after rather than a comparison against memory.

Timestamps in chain.log carry no date, so days are reconstructed from
midnight wraps (HH going 23 -> 00). Offline; touches nothing.

    python lab/e68_diurnal.py [--split "08-10 10:25"]
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DECODABLE = 60.0          # % FEC that counts as "watchable"


def parse_chain(path):
    """-> [(abs_minutes, fec_pct)] with days reconstructed from wraps."""
    out = []
    day = 0
    prev = None
    pat = re.compile(r"\[(\d\d):(\d\d):(\d\d)\][^\[]*\[\s*([\d.]+)% now\]")
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            m = pat.search(line)
            if not m:
                continue
            hh, mm, ss = int(m.group(1)), int(m.group(2)), int(m.group(3))
            t = hh * 60 + mm
            if prev is not None and t + 120 < prev:      # midnight wrap
                day += 1
            prev = t
            out.append((day * 1440 + t, float(m.group(4))))
    return out


def frac_decodable(rows):
    if not rows:
        return None, 0
    # one value per MINUTE (the log emits ~12/min); a minute counts as
    # decodable if its mean clears the bar
    per = {}
    for t, v in rows:
        per.setdefault(t, []).append(v)
    mins = {t: sum(v) / len(v) for t, v in per.items()}
    good = sum(1 for v in mins.values() if v >= DECODABLE)
    return 100.0 * good / len(mins), len(mins)


def by_hour(rows):
    per = {}
    for t, v in rows:
        per.setdefault(t, []).append(v)
    mins = {t: sum(v) / len(v) for t, v in per.items()}
    hh = {}
    for t, v in mins.items():
        hh.setdefault((t // 60) % 24, []).append(v)
    return {h: (100.0 * sum(1 for x in v if x >= DECODABLE) / len(v), len(v))
            for h, v in sorted(hh.items())}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chain", default=os.path.join(ROOT, "data", "e31",
                                                    "chain.log"))
    ap.add_argument("--split-min", type=int, default=None,
                    help="absolute minute index where the splitter went in; "
                         "default = auto (the last long 0%% run before a "
                         "sustained recovery)")
    a = ap.parse_args()
    rows = parse_chain(a.chain)
    if not rows:
        print("no FEC samples found")
        return 1
    t0, t1 = rows[0][0], rows[-1][0]
    print(f"chain.log: {len(rows)} FEC samples spanning "
          f"{(t1-t0)/60:.1f} h of log time")

    # the splitter went in at 10:25 on the LAST day of the log
    last_day = t1 // 1440
    split = a.split_min if a.split_min is not None else last_day * 1440 + 625
    pre = [r for r in rows if r[0] < split]
    post = [r for r in rows if r[0] >= split]
    fp, np_ = frac_decodable(pre)
    fq, nq = frac_decodable(post)
    print(f"\n=== PRE-SPLITTER  ({np_} minutes observed) ===")
    print(f"  decodable (FEC >= {DECODABLE:.0f}%): {fp:.1f}% of observed time")
    print(f"\n=== POST-SPLITTER ({nq} minutes observed) ===")
    if nq:
        print(f"  decodable (FEC >= {DECODABLE:.0f}%): {fq:.1f}% of observed time")
    print(f"\n=== BY HOUR OF DAY (% of observed minutes decodable) ===")
    hp, hq = by_hour(pre), by_hour(post)
    print(f"  {'hh':>3} {'PRE':>12} {'POST':>12}")
    for h in range(24):
        a_, b_ = hp.get(h), hq.get(h)
        sa = f"{a_[0]:5.1f}% n={a_[1]:<4d}" if a_ else "        --  "
        sb = f"{b_[0]:5.1f}% n={b_[1]:<4d}" if b_ else "        --  "
        print(f"  {h:02d}: {sa} {sb}")
    if nq:
        print(f"\n  NOTE: post-splitter coverage is only {nq} minutes so far; "
              f"the night hours are the ones that decide this, and they have "
              f"not happened yet.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
