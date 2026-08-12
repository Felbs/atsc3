#!/usr/bin/env python3
"""atsc3_lockstat.py -- measure how long RF33 holds a decodable lock.

A lock = payload actually decoding. We call it locked when interval FEC is
>= 80% and lost when it drops < 20% or the chain re-acquires. Reports every
lock's duration + the running distribution. This is the honest answer to
"how long can we stay on 107.1" -- a measurement, not an assertion.
"""
import json, os, re, sys, time

D = sys.argv[1] if len(sys.argv) > 1 else "data/e30"
CH = os.path.join("", D, "chain.log")
locks = []            # completed lock durations (s)
lock_start = None
last_wall = None
print(f"lockstat on {D}: waiting for RF33 payload to decode...", flush=True)
while True:
    try:
        lines = open(CH, encoding="utf-8", errors="replace").read().splitlines()
    except OSError:
        time.sleep(5); continue
    # most recent interval FEC %
    fec = None
    for ln in reversed(lines[-40:]):
        m = re.search(r"\[ *([0-9.]+)% now\]", ln)
        if m:
            fec = float(m.group(1)); break
    now = time.time()
    if fec is not None:
        if fec >= 80.0 and lock_start is None:
            lock_start = now
            print(f"{time.strftime('%H:%M:%S')} LOCKED (FEC {fec:.0f}%) -- 107.1 decoding", flush=True)
        elif fec < 20.0 and lock_start is not None:
            dur = now - lock_start
            locks.append(dur)
            lock_start = None
            mean = sum(locks)/len(locks)
            print(f"{time.strftime('%H:%M:%S')} LOST after {dur:.0f}s "
                  f"| locks={len(locks)} mean={mean:.0f}s longest={max(locks):.0f}s", flush=True)
    time.sleep(10)
