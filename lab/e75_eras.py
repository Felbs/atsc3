#!/usr/bin/env python3
"""e75_eras.py -- three eras, one receiver, one antenna, one parser. (E75)

THE QUESTION E73 FORCED. We have always called RF33's bad hours "band
fades". E73 proved that at least some of them were FRONT-END COMPRESSION
-- our own gain setting, not the sky: on the same air, rfgain_sel=2 gave
0.0% FEC and rfgain_sel=4 gave 100.0%. So the historical "71.8% of the
day decodable" may be partly self-inflicted, and the honest comparison is
three-way:

    (a) PRE-SPLITTER     direct feed,     rfgain_sel=2
    (b) SPLITTER ONLY    splitter in,     rfgain_sel=2
    (c) SPLITTER + GAIN  splitter in,     rfgain_sel=4   <- tonight

THE CONTROL THAT MAKES IT MEAN ANYTHING: the HDHomeRun shares the same
antenna through the splitter and has its own front end. If OUR decodable
fraction rises across eras while ITS snq stays flat, the improvement was
ours, not the sky's. That separation is only possible because both
receivers see one feed, and it is the whole point of the exercise.

Sources: data/e31/chain.log (ours, all eras -- days reconstructed from
midnight wraps since the log carries no dates) and
data/e31/_warden/diurnal_*.jsonl (both receivers, tonight only).

Offline. Touches nothing.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DECODABLE = 60.0

# Era boundaries. chain.log carries NO DATES, so days are reconstructed by
# counting midnight wraps: the log opens 8/09 13:36 = day 0, so 8/10 = day
# 1 and 8/11 = day 2.
#
# E75 CORRECTION: these were originally anchored to "the last day in the
# log", which silently BROKE the moment the log crossed into 8/11 -- the
# boundaries moved to a day that had not happened, every sample fell into
# era (a), and the hourly table began averaging pre-splitter 5% with
# post-fix 100% into a meaningless 50%. An era boundary is a FIXED point
# in history and must never be defined relative to a moving end.
EVENT_DAY = 1               # 8/10, the day both changes were made
SPLITTER_HHMM = (10, 25)    # splitter inserted
GAINFIX_HHMM = (20, 47)     # rfgain_sel 2 -> 4 adopted

# OPERATOR WINDOWS -- minutes when the SDR was DELIBERATELY OFF RF33 (RF25
# tests on RF25). During these the chain reads 0% while the HDHomeRun
# referee still reads snq~100 on RF33, which is EXACTLY the E73 "the air
# was fine and we failed" signature. Counting them would libel the very
# number that showed the night floor was ours -- a self-inflicted absence
# scored as a decode failure. They are excluded by construction, and
# listed here rather than hidden in a filter so the exclusion is auditable.
#   (day, start_hhmm, end_hhmm)
EXCLUDED_WINDOWS = [
    (2, (12, 4), (12, 37)),    # 8/11 RF25 window 1  -- RF33 outage only
    (2, (15, 40), (15, 46)),   # 8/11 RF25 window 2  -- RF33 outage only
]

# ANTENNA-CHANGE WINDOWS are excluded for a STRONGER reason and kept in a
# separate list so the distinction survives. The others are minutes when
# the SDR was merely pointed elsewhere; during these the PHYSICAL ANTENNA
# was different, so the minutes are not comparable to the era baseline
# even once the RF33 outage is accounted for. Every era number in this
# file describes Antenna B (Old Faithful) through the splitter; nothing
# measured on another port belongs in the same column.
#   (day, start_hhmm, end_hhmm_or_None)  -- None = still open
ANTENNA_WINDOWS = [
    (2, (16, 40), (17, 10)),   # 8/11 window 3: Philips antenna on PORT A.
    #   Boundaries are the OPERATOR's (16:40 / 17:10), not my stand-down
    #   clock (16:44) -- the antenna was already off B before I stopped
    #   contending, so my timestamp would have leaked 4 Philips minutes
    #   into era c. Close came as an explicit message, as agreed; the
    #   chain reappearing at 17:10 was NOT the signal (that inference
    #   cost a correct bring-up at 15:46, E84).
    #   *** VERDICT VOID (E89). *** The Philips was NOT PLUGGED IN. The
    #   "1 carrier vs 4, 1/22 the level, blind to RF25" numbers are a
    #   measurement of an OPEN COAX PORT, not of an antenna. They must
    #   never be quoted as an antenna comparison. The minutes stay
    #   excluded -- whatever port A was, it was not the yagi -- but the
    #   reason is now "unknown/disconnected input", not "worse antenna".
    (2, (17, 54), (18, 14)),   # 8/11 window 4: port A re-seated, RE-SWEEP.
    #   A genuine A/B against window 3, which is now a MEASURED
    #   DISCONNECTED BASELINE for this port: same tool, settings, evening
    #   and port, one variable (a coax actually attached).
    #   RESULT: the coax was the whole story. Connected port A saw RF25
    #   AND RF30 -- a second carrier, which is exactly the tell the
    #   connectedness test predicted -- and a LIVE SERVICE DECODED at
    #   94.5% FEC for 200 s, 1280x720, our own AC-4 doing the sound.
    #   THE SERVICE WAS WNUV 54.1 (The CW), NOT Fox -- see E94. RF25 is
    #   WBFF's transmitter; WBFF is the HOST LICENSEE, not the service
    #   we locked onto. The measurements stand; only the name changes.
    #   These minutes stay excluded: the chain was not on RF33/B.
    (2, (18, 35), None),       # 8/11 window 5: THE USER IS WATCHING TV.
    #   Live stack moved to RF25/WNUV 54.1 on Antenna A (Philips), later
    #   re-pointed at WBFFMob (the actual Fox service). Full viewer;
    #   viewer. OPEN-ENDED by design -- minutes drop out as they accrue
    #   and the close is an explicit message, never inferred (E84).
    #   Not an outage and not a fade: for this whole span there is no
    #   RF33 measurement to make, because nothing is pointed at RF33.
    #   Era c's decodable fraction simply does not cover these minutes.
    #   18:40-~20:20 inside this window the FLEET RADIO was held by the
    #   scheduled sonde hunt at priority 100, so there is no RF33 (or
    #   RF25) measurement for that span either -- and for a third distinct
    #   reason: not a fade, not our antenna, but a higher-priority owner.
    #   "No data" here means the radio was legitimately somebody else's.
]


def parse_chain(path):
    """-> [(abs_minute, fec)] with days reconstructed from midnight wraps."""
    out, day, prev = [], 0, None
    pat = re.compile(r"\[(\d\d):(\d\d):(\d\d)\][^\[]*\[\s*([\d.]+)% now\]")
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            m = pat.search(line)
            if not m:
                continue
            hh, mm = int(m.group(1)), int(m.group(2))
            t = hh * 60 + mm
            if prev is not None and t + 120 < prev:
                day += 1
            prev = t
            out.append((day * 1440 + t, float(m.group(4))))
    return out


def minutes(rows):
    """collapse ~12 samples/min into one value per minute"""
    per = {}
    for t, v in rows:
        per.setdefault(t, []).append(v)
    return {t: sum(v) / len(v) for t, v in per.items()}


def frac(mins):
    if not mins:
        return None, 0
    good = sum(1 for v in mins.values() if v >= DECODABLE)
    return 100.0 * good / len(mins), len(mins)


def referee_by_era(live_dir, split_abs, gain_abs, day_of_last):
    """HDHomeRun snq per era, from the diurnal log (tonight only)."""
    out = {}
    for p in sorted(glob.glob(os.path.join(live_dir, "_warden",
                                           "diurnal_*.jsonl"))):
        for line in open(p, encoding="utf-8", errors="replace"):
            try:
                r = json.loads(line)
            except ValueError:
                continue
            h = r.get("hdhr") or {}
            snq = h.get("snq")
            if snq is None:
                continue
            # E88: A RELEASED TUNER IS NOT A ZERO-SIGNAL READING. When
            # nothing holds tuner1 it reports lock=none / snq=0, which is
            # numerically identical to "RF33 is off the air" and would drag
            # the control's floor down as if the sky had done it. Window 3
            # released the tuner and era c's referee min duly read 0 while
            # our own FEC sat at 100% -- an apparent referee DISAGREEMENT
            # manufactured entirely by an idle instrument. Only count
            # samples where the tuner was actually locked to something.
            lock = str(h.get("lock") or "none").lower()
            if lock in ("none", "", "0"):
                continue
            hhmm = r.get("hhmm", "")
            m = re.search(r"(\d\d):(\d\d)$", hhmm)
            if not m:
                continue
            t = int(m.group(1)) * 60 + int(m.group(2))
            day = 1 if hhmm.startswith("08-10") else 2
            era = ("c_gainfix" if (day > 1 or t >= gain_abs % 1440)
                   else "b_splitter")
            out.setdefault(era, []).append(snq)
            out.setdefault("_ours_" + era, []).append(r.get("fec_5min") or 0.0)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--live-dir", default=os.path.join(ROOT, "data", "e31"))
    ap.add_argument("--event-day", type=int, default=EVENT_DAY,
                    help="reconstructed day index of 8/10")
    a = ap.parse_args()
    rows = parse_chain(os.path.join(a.live_dir, "chain.log"))
    if not rows:
        print("no samples")
        return 1
    mins = minutes(rows)
    # drop the operator windows BEFORE any era arithmetic
    excl, ant_excl = set(), set()
    for day, (h0, m0), (h1, m1) in EXCLUDED_WINDOWS:
        excl.update(range(day * 1440 + h0 * 60 + m0,
                          day * 1440 + h1 * 60 + m1 + 1))
    end_of_log = max(mins) if mins else 0
    for day, (h0, m0), end in ANTENNA_WINDOWS:
        start = day * 1440 + h0 * 60 + m0
        stop = (day * 1440 + end[0] * 60 + end[1]) if end else end_of_log
        ant_excl.update(range(start, stop + 1))
    n_before = len(mins)
    mins = {t: v for t, v in mins.items() if t not in excl and t not in ant_excl}
    d_op = n_before - len(mins)
    if d_op:
        print(f"excluded {d_op} minutes: operator windows (SDR off RF33) and "
              f"ANTENNA-CHANGE windows (different physical antenna --\n"
              f"  not comparable to the era baseline, which is Antenna B "
              f"through the splitter throughout)\n")
    split_abs = (a.event_day * 1440 + SPLITTER_HHMM[0] * 60
                 + SPLITTER_HHMM[1])
    gain_abs = (a.event_day * 1440 + GAINFIX_HHMM[0] * 60 + GAINFIX_HHMM[1])
    ndays = max(mins) // 1440 + 1
    print(f"log spans {ndays} reconstructed day(s); era boundaries pinned to "
          f"day {a.event_day} (8/10)\n")

    eras = {
        "a  PRE-SPLITTER   (direct feed, rfgain 2)":
            {t: v for t, v in mins.items() if t < split_abs},
        "b  SPLITTER ONLY  (splitter,    rfgain 2)":
            {t: v for t, v in mins.items() if split_abs <= t < gain_abs},
        "c  SPLITTER+GAIN  (splitter,    rfgain 4)":
            {t: v for t, v in mins.items() if t >= gain_abs},
    }
    print(f"chain.log: {len(rows)} samples, {len(mins)} distinct minutes\n")
    print(f"{'era':<44} {'decodable':>10} {'minutes':>8}")
    for name, m in eras.items():
        f, n = frac(m)
        print(f"{name:<44} {('%.1f%%' % f) if f is not None else '   n/a':>10}"
              f" {n:>8}")

    print(f"\n=== BY HOUR (% of observed minutes decodable) ===")
    print(f"  {'hh':>3} {'a pre':>12} {'b splitter':>12} {'c +gain':>12}")
    for h in range(24):
        cells = []
        for name, m in eras.items():
            vals = [v for t, v in m.items() if (t // 60) % 24 == h]
            if vals:
                g = 100.0 * sum(1 for v in vals if v >= DECODABLE) / len(vals)
                cells.append(f"{g:5.1f}% n={len(vals):<3d}")
            else:
                cells.append("       --   ")
        print(f"  {h:02d}: {cells[0]} {cells[1]} {cells[2]}")

    ref = referee_by_era(a.live_dir, split_abs, gain_abs, a.event_day)
    print(f"\n=== THE CONTROL: what did the REFEREE see? ===")
    print(f"  (its own front end, same antenna via the splitter)")
    for era in ("b_splitter", "c_gainfix"):
        s = ref.get(era) or []
        o = ref.get("_ours_" + era) or []
        if not s:
            continue
        s_sorted = sorted(s)
        print(f"  {era:<12} referee snq median {s_sorted[len(s)//2]:5.1f} "
              f"(min {min(s)}, max {max(s)}, n={len(s)})   "
              f"ours FEC median "
              f"{sorted(o)[len(o)//2] if o else float('nan'):5.1f}%")
    print("\n  READ IT THIS WAY: if the referee's snq is FLAT across eras "
          "while\n  our decodable fraction rises, the improvement was OURS, "
          "not the sky's.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
