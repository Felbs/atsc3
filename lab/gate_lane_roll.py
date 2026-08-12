#!/usr/bin/env python3
"""GATE: does a worker survive the chain restarting underneath it?

On 8/07 the supervisor did its job for the first time -- it detected a
genuinely wedged chain (heartbeat stale 159 s, all lanes frozen) and restarted
it. The chain came back at 100 % FEC. Everything downstream did not:

  * `LiveWriter` opened each lane "wb", discarding 433.8 MB of recorded TV;
  * `LaneReader.read_new` compared `size <= self.pos` and returned [] on every
    pass thereafter, so the audio worker sat at 1277.3 s reporting `behind 0`
    -- a cursor past EOF reads as "caught up" -- and never emitted again.

The subtitle worker survived the same event untouched, because it re-reads
documents instead of holding a byte cursor. The component that survived was
the one not holding fragile state.

This gate covers the reader half: shrink the file underneath a live reader and
require that it notices and recovers. The writer half (rolling rather than
truncating) is checked by `gate_lane_roll_writer` below.

Run:  python lab/gate_lane_roll.py
"""
from __future__ import annotations

import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "tools"))

from atsc3_audio import LaneReader          # noqa: E402


def real_frags(n=14):
    """n real fragments (+ the init) sliced out of a lane recorded off air.

    Hand-built boxes were not rich enough for `trun_scan` to accept -- the
    first version of this gate reported 0 fragments and correctly refused to
    grade the code on a harness that fed it nothing. Real broadcaster bytes
    remove the whole question.
    """
    root = os.path.dirname(HERE)
    import glob
    import json as _json
    for d in sorted(glob.glob(os.path.join(root, "data", "*")), reverse=True):
        lane = os.path.join(d, "live_audio_pid13.m4s")
        idx = os.path.join(d, "live_audio_pid13.idx")
        meta = os.path.join(d, "live.json")
        if not (os.path.exists(lane) and os.path.exists(idx)):
            continue
        try:
            rows = [_json.loads(l) for l in open(idx)]
            init_n = _json.load(open(meta))["lanes"]["audio_pid13"]["init_bytes"]
        except Exception:                                      # noqa: BLE001
            continue
        if len(rows) < n:
            continue
        raw = open(lane, "rb").read()
        init = raw[:init_n]
        frags = [raw[r["off"]:r["off"] + r["len"]] for r in rows[:n]]
        if all(frags):
            return init, frags, d
    return None, None, None


def main():
    init, frags, src = real_frags()
    if init is None:
        print("GATE INVALID -- no recorded lane to slice real fragments from")
        return 2
    print(f"  fragments from {os.path.basename(src)}")
    d = tempfile.mkdtemp()
    path = os.path.join(d, "live_audio_pid13.m4s")

    # ---- generation 1: the reader consumes some fragments ----------------
    with open(path, "wb") as f:
        f.write(init + b"".join(frags[:10]))
    r = LaneReader(path)
    r.read_new()
    pos1, nfrag1 = r.pos, r.nfrag
    print(f"  generation 1 : cursor {pos1}, {nfrag1} fragments seen")
    if nfrag1 == 0:
        print("GATE INVALID -- the harness never fed a readable fragment")
        return 2

    # ---- the chain restarts: the lane is recreated, SHORTER --------------
    with open(path, "wb") as f:
        f.write(init + b"".join(frags[:3]))
    size2 = os.path.getsize(path)
    print(f"  restart      : file is now {size2} bytes, cursor is {pos1}")

    r.read_new()
    print(f"  after roll   : cursor {r.pos}, {r.nfrag} fragments, "
          f"rolls={r.rolls}")

    # ---- and it must keep following the NEW file -------------------------
    with open(path, "ab") as f:
        f.write(b"".join(frags[3:7]))
    r.read_new()
    print(f"  still following: {r.nfrag} fragments after 4 more appended")

    ok = r.rolls == 1 and r.nfrag == 7 and r.pos == os.path.getsize(path)
    print()
    if ok:
        print("GATE PASS -- reader noticed the roll and followed the new lane")
        return 0
    print(f"GATE FAIL -- reader is stranded: rolls={r.rolls}, "
          f"nfrag={r.nfrag} (want 7), cursor={r.pos} vs file "
          f"{os.path.getsize(path)}. A restart silently ends the audio.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
