#!/usr/bin/env python3
"""gate_e86.py -- a rejected CTI phase candidate must not kill the receiver.

THE BUG.  CtiStream enforces A/322's "C must land inside one FEC Block" (the
m10_core P3 gate) in two places, and until E86 the two disagreed about what
failing it MEANS:

    rewind_to()  ->  returns False, the sweep moves to the next candidate
    reset()      ->  raises ValueError, which nothing caught

_begin() calls reset().  The blind re-begin path reaches it with a
dead-reckoned start_row and a fec_block_start left over from an older L1
anchor, and on a marginal carrier those two disagree routinely.  So the SAME
event was housekeeping during a sweep and fatal while starting one, and the
receiver died exactly when the air got bad enough to need re-acquisition --
taking the evidence of the fade with it (E82's "a crash truncates the
evidence" again, one layer down).

Observed on live RF25 (Fox), 2026-08-11:
    ValueError: C=531124 outside [0,32400) -- start_row 720 /
                fec_block_start 928436 disagree

WHAT THIS GATE PINS.  The fix must make rejection ROUTINE without making the
gate TOOTHLESS, so both directions are checked, and every check is built to
be able to fail:

  G1  an in-range candidate is still ADOPTED           (fix didn't eat everything)
  G2  an out-of-range candidate does not raise         (the crash is gone)
  G3  ... and is not adopted either -- the stream is left cold rather than
      anchored on a C the gate just refused            (no silent corruption)
  G4  a later in-range candidate is reached after earlier ones are rejected
                                                       (the sweep continues)
  G5  reset() ITSELF still raises on a bad pair        (gate keeps authority)
  G6  rewind_to() still reports False on a bad pair, and its caller now acts
      on that instead of ignoring it

NEGATIVE CONTROL: G1/G4 fail if the fix over-rejects; G3/G5 fail if it
under-rejects.  A change that merely swallowed the exception passes G2 and
fails G3.
"""
from __future__ import annotations

import collections
import os
import sys
import types

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m10_cti as CTI            # noqa: E402
import m44_ldm as M44            # noqa: E402

NROWS, NCELL, PLP = 1024, 32400, 1133282


def _plan():
    return types.SimpleNamespace(nrows=NROWS, ncell=NCELL, plp_size=PLP,
                                 cti_reach=NROWS * (NROWS - 1))


def _find_fbs(start_row, want_in_range):
    """A fec_block_start whose solved C is inside / outside [0, NCELL)."""
    for fbs in range(0, 4_000_000, 7):
        C = CTI.solve_C(fbs, start_row, NROWS)
        if (0 <= C < NCELL) == want_in_range:
            return fbs, C
    raise AssertionError("no such fec_block_start")


class _Ph:
    """Minimal PhaseTracker stand-in: hands out candidates in order."""
    LOCK_BLOCKS = 4

    def __init__(self, cands, fec_block_start):
        self._c = list(cands)
        self.fec_block_start = fec_block_start
        self.state = "cold"

    def begin(self, frame_idx, blind=False):
        return self._c.pop(0) if self._c else None

    def next_candidate(self):
        return self._c.pop(0) if self._c else None


class _Pipe:
    """Just enough of LdmPipeline for the real _begin to run against."""

    def __init__(self, ph, cti):
        self.ph, self.cti = ph, cti
        self.stats = collections.Counter()
        self.origin = None
        self.l1_tries = 3
        self.lines = []

    def log(self, msg):
        self.lines.append(msg)

    _begin = M44.LdmPipeline._begin


def _fresh(cands, fbs):
    cti = M44.CtiStream(_plan())
    return _Pipe(_Ph(cands, fbs), cti)


def main():
    ok = True

    def check(name, cond, detail=""):
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}"
              + (f"   {detail}" if detail else ""))

    print("gate_e86 -- a rejected CTI phase candidate must not kill the run")

    sr_ok = 344
    fbs_ok, C_ok = _find_fbs(sr_ok, True)
    sr_bad = 720
    fbs_bad, C_bad = _find_fbs(sr_bad, False)
    print(f"  in-range  pair: start_row {sr_ok} fbs {fbs_ok} -> C {C_ok}")
    print(f"  out-range pair: start_row {sr_bad} fbs {fbs_bad} -> C {C_bad}")

    # G1 -- a good candidate is still adopted
    p = _fresh([(sr_ok, 0)], fbs_ok)
    p._begin(5)
    check("G1 in-range candidate ADOPTED", p.origin == 5
          and p.cti.origin_frame == 5 and p.cti.C0 == C_ok,
          f"origin={p.origin} C0={p.cti.C0}")

    # G2/G3 -- a bad candidate is rejected, not fatal, and not adopted
    p = _fresh([(sr_bad, 0)], fbs_bad)
    raised = None
    try:
        p._begin(9)
    except Exception as e:                                   # noqa: BLE001
        raised = e
    check("G2 out-of-range candidate does NOT raise", raised is None,
          f"raised={raised!r}")
    check("G3 ... and is NOT adopted (left cold)",
          p.origin is None and p.cti.origin_frame is None
          and p.stats["phase_rejected_range"] == 1,
          f"origin={p.origin} cti.origin={p.cti.origin_frame} "
          f"rejected={p.stats['phase_rejected_range']}")

    # G4 -- the sweep keeps going and reaches a good candidate behind bad ones
    p = _fresh([(sr_bad, 0), (sr_bad, +1), (sr_ok, +2)], fbs_bad)
    p.ph.fec_block_start = fbs_bad
    # make only the third candidate solvable, by swapping the anchor when it
    # is reached: emulate the real tracker, whose fec_block_start is fixed --
    # so instead pick a start_row that IS in range under fbs_bad
    sr_ok2 = None
    for s in range(NROWS):
        if 0 <= CTI.solve_C(fbs_bad, s, NROWS) < NCELL:
            sr_ok2 = s
            break
    p = _fresh([(sr_bad, 0), (sr_bad, +1), (sr_ok2, +2)], fbs_bad)
    p._begin(11)
    check("G4 sweep continues past bad candidates to a good one",
          sr_ok2 is not None and p.origin == 11
          and p.stats["phase_rejected_range"] == 2,
          f"good start_row={sr_ok2} origin={p.origin} "
          f"rejected={p.stats['phase_rejected_range']}")

    # G5 -- the gate itself still has teeth
    cti = M44.CtiStream(_plan())
    try:
        cti.reset(0, sr_bad, fbs_bad)
        raised5 = None
    except ValueError as e:
        raised5 = e
    check("G5 reset() STILL raises on a bad pair (gate keeps authority)",
          raised5 is not None, f"raised={type(raised5).__name__ if raised5 else None}")

    # G6 -- rewind_to reports the same condition as a return value
    cti = M44.CtiStream(_plan())
    cti.reset(0, sr_ok, fbs_ok)
    r_bad = cti.rewind_to(sr_bad, fbs_bad)
    cti2 = M44.CtiStream(_plan())
    cti2.reset(0, sr_ok, fbs_ok)
    r_ok = cti2.rewind_to(sr_ok, fbs_ok)
    check("G6 rewind_to False on bad pair / True on good",
          r_bad is False and r_ok is True, f"bad={r_bad} good={r_ok}")

    print(f"\ngate_e86: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
