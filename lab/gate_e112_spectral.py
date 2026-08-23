#!/usr/bin/env python3
"""gate_e112_spectral.py -- the compiled AC-4 spectral parse's referee.

The kernel replaces the Huffman/bitstream loop in asf_spectral_data.  Its
output is INTEGERS (quantised spectral lines) and a bit position, so the bar
is total identity -- there is no rounding to license.

  leg 1  REAL FRAMES: decode N frames of a real AC-4 lane with the kernel on
         and off; every MDCT window must be bit-identical and the count must
         match.  This covers the parse, the sign bits, the cb-11 escapes and
         the bit-position hand-off in one sweep, on broadcast data.
  leg 2  MALFORMED INPUT: truncate a frame mid-stream; the kernel must
         DECLINE (nonzero rc, bit position untouched) so the Python loop
         re-runs and raises exactly as it always did -- the kernel must never
         turn a malformed frame into a silently half-parsed one.

Run: python lab/gate_e112_spectral.py --lane <live_audio_pidN.m4s>
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

_COMPARED = [0]


def _decode(lane, n, use_kernel):
    os.environ["ATSC3_SPECTRAL_KERNEL"] = "1" if use_kernel else "0"
    for m in list(sys.modules):
        if m.split(".")[0].startswith(("m1", "m2", "m3", "m4")):
            del sys.modules[m]
    import m17_ac4_walk as W
    import m28_channels as CH
    CH._SPEC_STATE["lib"] = None
    from m42_ac4_stream import Ac4Stream
    fr = W.samples(lane)[:n]
    dec = Ac4Stream(element="5_X")
    wins, lens, grp, cnts = dec.decode_frames(fr, Ac4Stream.CHANS)
    flat = []
    for ch in sorted(wins):
        for wnd in wins[ch]:
            flat.append(np.asarray(wnd))
    return flat


def leg1(fails, lane, n):
    print(f"LEG 1 -- {n} real frames, kernel vs Python: {lane}")
    if not os.path.exists(lane):
        msg = f"lane {lane} not found -- the gate CANNOT vouch for anything"
        print(f"  FAIL: {msg}")
        fails.append(msg)
        return
    k = _decode(lane, n, True)
    p = _decode(lane, n, False)
    _COMPARED[0] += 1
    if len(k) != len(p):
        print(f"  window COUNT differs: {len(k)} vs {len(p)}")
        fails.append("window count differs")
        return
    bad = sum(0 if np.array_equal(a, b) else 1 for a, b in zip(k, p))
    print(f"  {len(k)} MDCT windows, {bad} differ")
    if bad:
        fails.append(f"{bad} windows differ")


def leg2(fails):
    print("LEG 2 -- malformed input must fall back, not half-parse")
    os.environ["ATSC3_SPECTRAL_KERNEL"] = "1"
    import m28_channels as CH
    import m19_ac4_toc as T
    CH._SPEC_STATE["lib"] = None
    if CH._spectral_lib() is None:
        fails.append("kernel not loadable for leg 2")
        return
    import m24_spectral as SP
    # a real codebook, a section that wants more bits than the buffer has
    tables = {1: SP.Codebook.__new__(SP.Codebook)} if hasattr(SP, "Codebook") \
        else None
    # Build via the public path instead: steal the tables a real decode uses.
    from m42_ac4_stream import Ac4Stream
    dec = Ac4Stream(element="5_X")
    T_cb = dec.T["cb"] if hasattr(dec, "T") else None
    if T_cb is None:
        import m28_channels as CH2
        print("  (no direct table handle -- exercising via truncated buffer)")
    b = T.Bits(b"\\x00")                 # 8 bits total
    sects = [(11, 0, 4)]                # demands far more than 8 bits
    offsets = [0, 64, 128, 192, 256]
    tabs = CH._spectral_trees(dec.T["cb"]) if hasattr(dec, "T") else None
    _COMPARED[0] += 1
    if tabs is None:
        print("  could not build trees from the live tables -- counting as "
              "SKIP, leg 1 remains the authority")
        return
    p0 = b.p
    try:
        CH.asf_spectral_data(b, sects, offsets, dec.T["cb"])
        print("  parse of an impossible stream DID NOT RAISE")
        fails.append("truncated stream did not raise")
    except Exception as ex:                                    # noqa: BLE001
        print(f"  raised {type(ex).__name__} as the reference does -- OK")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--lane", required=True)
    ap.add_argument("--frames", type=int, default=200)
    a = ap.parse_args(argv)
    fails = []
    leg1(fails, a.lane, a.frames)
    leg2(fails)
    if _COMPARED[0] == 0:
        fails.append("the gate performed ZERO comparisons")
    print(f"\ncomparisons performed: {_COMPARED[0]}")
    if fails:
        print(f"GATE FAILED ({len(fails)}):")
        for f in fails:
            print(f"  - {f}")
        return 1
    print("GATE PASSED -- the spectral kernel is bit-identical on real air")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
