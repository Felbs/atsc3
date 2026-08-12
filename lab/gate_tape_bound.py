#!/usr/bin/env python3
"""GATE: can the front-end tape grow without bound when the cursor is lost?

On 8/07 the chain "wedged" twice, ~24 minutes into each run. It was not a
deadlock. `atsc3 watch` RSS went 1.2 GB -> 34.9 GB in about three minutes
while the Frame counter sat frozen at 457 and x-real-time decayed 0.74 ->
0.44. The supervisor correctly called it a wedge and restarted it; nobody
knew why.

The mechanism: `FrontEnd.frames()` trims the tape ONLY inside its loop. When
`self.centre` is unreachable the loop breaks on the first test and trims
nothing, while `_stream` keeps appending complex128 at ~110 MB/s -- which
matches the observed growth rate.

`self.centre` becomes unreachable because `reacquire()` runs on the DECODE
thread and replaces both `self.tape` and `self.centre` while `frames()` runs
on the FRONT thread. The identical race on `self.acq` was already found and
fixed with a generation snapshot (see `_acquire`); the same reset's effect on
the tape cursor was not.

This gate covers the BOUND rather than the race: force an unreachable cursor
and require that the tape stops growing. A race is timing-dependent and a
gate for it would be flaky; a bound is not, and it is the property that
actually matters -- a buffer fed by a radio must never be unbounded, whatever
the reason it stopped draining.

Run:  python lab/gate_tape_bound.py
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import m11_stream as S            # noqa: E402


def main():
    fe = S.FrontEnd.__new__(S.FrontEnd)
    fe.ready = True
    fe.gen = 0
    fe.stats = __import__("collections").Counter()
    fe.tape = S.SampleTape()
    fe.tm = __import__("collections").Counter()
    fe.ex = None

    # The lost cursor: `centre` measured against a tape that no longer exists.
    fe.centre = 10 ** 9

    reacquired = []
    fe.reacquire = lambda: reacquired.append(fe.tape.end - fe.tape.org)

    # Feed 20 s of air, a block at a time, the way `_stream` would.
    block = np.zeros(1 << 16, np.complex128)
    fed = 0
    peak = 0
    for _ in range(int(20.0 * S.FS_POST) // len(block)):
        fe.tape.append(block)
        fed += len(block)
        fe.frames()
        peak = max(peak, fe.tape.end - fe.tape.org)
        if reacquired:
            break

    held_s = peak / S.FS_POST
    fed_s = fed / S.FS_POST
    print(f"  fed {fed_s:5.1f} s of air with the cursor unreachable")
    print(f"  peak tape held      {held_s:5.1f} s "
          f"({peak * 16 / 1e6:.0f} MB as complex128)")
    print(f"  tape_overrun fired  {fe.stats['tape_overrun']}  "
          f"reacquire called {len(reacquired)}")

    ok = bool(reacquired) and held_s <= 9.0
    print()
    if ok:
        print(f"GATE PASS -- tape capped at {held_s:.1f} s and the front end "
              f"re-acquired instead of growing")
        return 0
    print(f"GATE FAIL -- tape reached {held_s:.1f} s "
          f"({peak * 16 / 1e9:.1f} GB) and nothing stopped it. At 110 MB/s "
          f"this reaches 34 GB in about five minutes.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
