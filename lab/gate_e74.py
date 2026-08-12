#!/usr/bin/env python3
"""gate_e74.py -- an MPU's self-declared duration must never price a gap.

MEASURED 8/10 on e31 (video lane .0019, seq 225821274): a damaged MPU
arrived carrying 93 of its 120 frames, and its media `trun` sample
durations summed to 720720 ticks = EXACTLY 4.0000 MPU instead of 180180.
LiveWriter._segment then did BOTH of the things that turn one bad field
into a timeline error:

    st["base_time"] += dur                      # 4 slots, once
    st["last_dur"]   = dur                      # and it becomes the ruler
    ...
    st["base_time"] += missing * st["last_dur"] # 5 missing x 4 = 20 slots

so 6 sequence numbers of air advanced the lane clock by 24 slots --
+18.0000 slots = +36.036 s injected, permanently, into every later
fragment. The viewer then stamped video 36 s off its own audio (E71) and
VLC discarded 82.7 % of pictures. The audio lane over the identical
sequence range shows +0.025 s, which is what a healthy lane looks like.

This is the ATSC 3.0 unprotected-header law again: the trun is a field
read from the air, and an unbounded field read from the air is a bug
waiting for a bad day. The fix bounds it -- a per-MPU duration is
CLAMPED to the service's nominal cadence, and gaps are priced with that
nominal, never with whatever the last MPU happened to claim.

    python lab/gate_e74.py    -> PASS/FAIL, exit 0 iff all pass

Nothing here is imported by the live stack; this file is the reference
implementation of the PROPOSED fix plus its gate. The patch to
lab/m11_stream.py is NOT applied (data/e31 is frozen).
"""
from __future__ import annotations

import sys

MPU_TICKS = 180180          # one MPU on this service, measured constant
TOL = 0.25                  # a whole MPU may honestly vary by this much


# --------------------------------------------------------------- shipped
def advance_shipped(events, nominal=MPU_TICKS):
    """The CURRENT LiveWriter accounting. events: (seq, dur_ticks) for the
    MPUs that were actually written; gaps are inferred from seq."""
    base, next_seq, last_dur = 0, None, 0
    stamps = []
    for seq, dur in events:
        if next_seq is not None and seq > next_seq:
            missing = seq - next_seq
            if last_dur:
                base += missing * last_dur       # <-- priced by ONE sample
        stamps.append((seq, base))
        base += dur
        last_dur = dur
        next_seq = seq + 1
    return stamps


# -------------------------------------------------------------- proposed
def advance_bounded(events, nominal=MPU_TICKS, tol=TOL):
    """PROPOSED: bound the field, and price gaps with the SERVICE's cadence.

    Two independent changes, either of which alone would have prevented
    the live fault, kept together because they answer different questions:
      * `dur` is clamped to the nominal cadence when it deviates beyond
        tol -- a fragment whose own sample table is damaged must not move
        the clock further than a whole MPU can;
      * gaps are priced with `nominal`, never with `last_dur`, so a bad
        measurement can never be MULTIPLIED by the size of a hole.
    Returns (stamps, n_clamped) so a caller can log the anomaly rather
    than absorb it silently."""
    base, next_seq = 0, None
    stamps, clamped = [], 0
    for seq, dur in events:
        if next_seq is not None and seq > next_seq:
            base += (seq - next_seq) * nominal
        use = dur
        if nominal and not (nominal * (1 - tol) <= dur <= nominal * (1 + tol)):
            use = nominal
            clamped += 1
        stamps.append((seq, base))
        base += use
        next_seq = seq + 1
    return stamps, clamped


FAILS = []


def check(name, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"  ({detail})" if detail
                                                     else ""))
    if not ok:
        FAILS.append(name)


def live_case():
    """The measured event: healthy run, one 4x MPU, then a 5-MPU hole."""
    ev = [(225821270 + k, MPU_TICKS) for k in range(5)]     # ..270-274
    ev[-1] = (225821274, 4 * MPU_TICKS)                     # the damaged one
    ev.append((225821280, MPU_TICKS))                       # after the hole
    return ev


def gate_live_event():
    print("gate 1: the measured 8/10 event")
    ev = live_case()
    a = dict(advance_shipped(ev))
    got = a[225821280] - a[225821274]
    check("NEGATIVE CONTROL: shipped accounting injects +18.0000 slots",
          got == 24 * MPU_TICKS,
          f"advanced {got} ticks = {got/MPU_TICKS:.4f} slots across "
          f"6 sequence numbers (excess {(got-6*MPU_TICKS)/MPU_TICKS:+.4f})")

    stamps, clamped = advance_bounded(ev)
    b = dict(stamps)
    got2 = b[225821280] - b[225821274]
    check("bounded accounting advances exactly 6 slots",
          got2 == 6 * MPU_TICKS,
          f"advanced {got2} ticks = {got2/MPU_TICKS:.4f} slots")
    check("the damaged MPU is reported, not silently absorbed", clamped == 1,
          f"clamped {clamped}")


def gate_healthy_unchanged():
    print("gate 2: healthy air is untouched (no behaviour change)")
    ev = [(1000 + k, MPU_TICKS) for k in range(50)]
    s_old = dict(advance_shipped(ev))
    s_new, clamped = advance_bounded(ev)
    check("stamps identical to the shipped path on a clean lane",
          s_old == dict(s_new) and clamped == 0)

    # a clean lane WITH an honest hole: both must hold the clock true
    ev = [(2000 + k, MPU_TICKS) for k in range(10)]
    ev += [(2020 + k, MPU_TICKS) for k in range(10)]        # 10-MPU hole
    s_old = dict(advance_shipped(ev))
    s_new, clamped = advance_bounded(ev)
    check("a real hole is still priced at one MPU each (E20's lesson kept)",
          s_old == dict(s_new) and clamped == 0
          and s_old[2020] - s_old[2009] == 11 * MPU_TICKS,
          f"{(s_old[2020]-s_old[2009])/MPU_TICKS:.1f} slots across a 10-hole")

    # the benign per-MPU jitter measured in BOTH lanes (+-258 ticks) must
    # pass through unclamped -- it is real and it is common-mode
    ev = [(3000, MPU_TICKS), (3001, MPU_TICKS + 258),
          (3002, MPU_TICKS - 257), (3003, MPU_TICKS)]
    _, clamped = advance_bounded(ev)
    check("+-258-tick jitter is NOT clamped (it is honest air)", clamped == 0)


def gate_amplification():
    print("gate 3: a bad duration must never be multiplied by a hole")
    worst = 0
    for hole in (1, 5, 20, 100):
        ev = [(10, MPU_TICKS), (11, 4 * MPU_TICKS), (11 + 1 + hole, MPU_TICKS)]
        a = dict(advance_shipped(ev))
        err_old = (a[11 + 1 + hole] - a[11]) - (1 + hole) * MPU_TICKS
        s, _ = advance_bounded(ev)
        b = dict(s)
        err_new = (b[11 + 1 + hole] - b[11]) - (1 + hole) * MPU_TICKS
        worst = max(worst, err_old / MPU_TICKS)
        if err_new != 0:
            check(f"bounded stays exact across a {hole}-MPU hole", False,
                  f"err {err_new}")
    check("bounded is exact for every hole size tried", True,
          "holes 1/5/20/100")
    check("NEGATIVE CONTROL: shipped error grows WITH the hole",
          worst >= 300, f"worst-case shipped error {worst:.0f} slots "
          f"({worst*2.002:.0f} s) on a 100-MPU hole")


if __name__ == "__main__":
    gate_live_event()
    gate_healthy_unchanged()
    gate_amplification()
    print(f"\n{'ALL GATES PASS' if not FAILS else 'FAILED: ' + ', '.join(FAILS)}")
    sys.exit(1 if FAILS else 0)
