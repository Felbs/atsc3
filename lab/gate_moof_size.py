#!/usr/bin/env python3
"""GATE: does a corrupted sample size in a moof allocate the machine away?

The wedge hunt on 8/07 cost two wrong diagnoses before this. What the
evidence finally forced:

  * `atsc3 watch` RSS 968 MB -> 32 GB in ~100 s. That is ~310 MB/s, and the
    radio only delivers 55 MB/s, so NOTHING that accumulates a stream can
    produce it. It has to be a single oversized allocation.
  * Identical on `--accel gpu-full` and `--accel cpu`, so it is not the GPU.
  * `tape_overrun 0`, `frames_raced 0`, so it is not the front-end tape.
  * It only ever happens during re-acquisition storms -- i.e. when the signal
    is bad and corrupted headers get through.

`MmtpFlow.build` does `bytearray(want)` TWICE for every entry of the trun's
sample-size table, and that table is read straight out of a `moof` carried by
MMTP, which has no integrity protection at all. One flipped bit asks for up
to 4 GB. This is E21's exposure exactly -- an unvalidated field from an
unprotected header -- with a length in place of a sequence number.

This gate is honest about its own limits: it proves the MECHANISM is lethal
and that the bound stops it. It does NOT by itself prove this is what wedged
the chain. That needs a live run where `sizes_rejected` fires while RSS stays
flat.

Run:  python lab/gate_moof_size.py
"""
from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import m7_objects as O7          # noqa: E402


def trun(sizes):
    """A moof carrying `sizes` as its media trun sample-size table."""
    n = len(sizes)
    flags = 0x200                                   # sample-size present
    body = b"\x00" + flags.to_bytes(3, "big") + n.to_bytes(4, "big")
    body += b"".join(s.to_bytes(4, "big") for s in sizes)
    tr = (len(body) + 8).to_bytes(4, "big") + b"trun" + body
    tfhd_body = b"\x00\x00\x00\x08" + (1).to_bytes(4, "big") + \
                (0).to_bytes(4, "big")
    tf = (len(tfhd_body) + 8).to_bytes(4, "big") + b"tfhd" + tfhd_body
    traf = (len(tf) + len(tr) + 8).to_bytes(4, "big") + b"traf" + tf + tr
    moof = (len(traf) + 8).to_bytes(4, "big") + b"moof" + traf
    return moof


def main():
    good = [40000] * 120                     # a real 2.002 s video MPU
    mi = O7.moof_info(trun(good))
    n_good = len(mi["sizes"])
    print(f"  healthy moof   : {n_good} samples accepted, "
          f"rejected={mi['sizes_rejected']}")
    if n_good != 120:
        print("GATE INVALID -- the harness cannot even express a good moof")
        return 2

    # ONE flipped bit in the top byte of a single sample size.
    poisoned = list(good)
    poisoned[7] = good[7] | 0x4000_0000
    mi = O7.moof_info(trun(poisoned))
    would = sum(mi["sizes"]) * 2 if mi["sizes"] else 0
    print(f"  one bit flipped: {len(mi['sizes'])} samples accepted, "
          f"rejected={mi['sizes_rejected']}")
    print(f"                   build() would have allocated "
          f"{would / 1e9:.2f} GB")

    # A corrupted sample COUNT: the loop length itself comes off the air.
    many = O7.moof_info(trun([1000] * 9000))
    print(f"  9000 samples   : {len(many['sizes'])} accepted, "
          f"rejected={many['sizes_rejected']}")

    ok = (n_good == 120 and not mi["sizes"] and mi["sizes_rejected"]
          and not many["sizes"] and many["sizes_rejected"])
    print()
    if ok:
        print("GATE PASS -- healthy tables pass, impossible ones are refused "
              "before anything is allocated")
        return 0
    print("GATE FAIL -- an unbounded sample size still reaches the allocator")
    return 1


if __name__ == "__main__":
    sys.exit(main())
