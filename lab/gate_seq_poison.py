#!/usr/bin/env python3
"""GATE: does ONE corrupted mpu_sequence_number kill an asset permanently?

E20 (8/07) lost video and both audio lanes during a fade while the subtitle
lane -- same chain, same FEC, same second -- kept writing.  The RF recovered
to 100% interval FEC and the lanes never came back.  Two latches in the
receive path are monotone functions of a field that arrives with NO integrity
protection at all:

    m7_objects.mpu_payload:  seq = int.from_bytes(b[4:8], "big")   # unchecked
    m11_stream.MpuStreamer:  self.hi[pid] = max(self.hi[pid], seq)
    m11_stream.Transport:    st["next_seq"] = m["seq"] + 1

MMTP headers carry no CRC.  A single bit flip in the top bit of that field
adds 2**31.  This gate asks whether that is survivable, by driving the real
MpuStreamer with synthetic MPUs and patching ONLY the box parsing/building --
the latch logic under test is the shipped code, untouched.

Run:  python lab/gate_seq_poison.py
"""
from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import m7_objects as O7          # noqa: E402
import m11_stream as S           # noqa: E402

PID = 12
WANT = 10                        # bytes of media per synthetic MPU


# ---- fakes: only the byte-level parsing, never the latch logic ------------
def fake_mmtp_parse(p, ver_ext):
    return dict(ver=0, typ=0, pid=p["pid"], psn=p["psn"], body=p)


def fake_mpu_payload(b):
    du = b"\x00" * (O7.DU_HDR + WANT) if b["ft"] == 2 else b"\x00" * 8
    return dict(ft=b["ft"], timed=1, fi=0, agg=0, fc=0, seq=b["seq"], dus=[du])


def fake_moof_info(mfmd):
    return dict(sizes=[WANT], mdat_declared=WANT, hint_size=0, hint_n=0)


def fake_build(self, pid, seq, e, partial=False):
    if e.get("sizes") is None or e["got"] < e["want"]:
        return None                      # incomplete -> genuinely lost
    return dict(seq=seq, truncate_to=None, n_samples=1, pid=pid)


def install():
    O7.mmtp_parse = fake_mmtp_parse
    O7.mpu_payload = fake_mpu_payload
    O7.moof_info = fake_moof_info
    S.MpuStreamer._build = fake_build


def mpu(st, seq, psn):
    """Feed one complete MPU (moof fragment, then media fragment)."""
    out = []
    out += st.feed(dict(pid=PID, psn=psn, seq=seq, ft=1))
    out += st.feed(dict(pid=PID, psn=psn + 1, seq=seq, ft=2))
    return out


def main():
    install()
    base = 225_698_000

    # ---- 1. baseline: a clean run emits every MPU ------------------------
    st = S.MpuStreamer(ver_ext=0)
    psn = 0
    emitted = 0
    for i in range(10):
        emitted += len(mpu(st, base + i, psn))
        psn += 2
    print(f"  clean       : {emitted}/10 MPUs emitted")
    if emitted != 10:
        print("GATE INVALID -- the harness cannot even pass a clean stream")
        return 2

    # ---- 2. one bit flip in the top bit of mpu_sequence_number -----------
    st = S.MpuStreamer(ver_ext=0)
    psn = 0
    before = 0
    for i in range(10):
        before += len(mpu(st, base + i, psn))
        psn += 2

    poison = (base + 10) | 0x8000_0000          # ONE bit
    st.feed(dict(pid=PID, psn=psn, seq=poison, ft=1))
    psn += 1
    hi_after = st.hi[PID]

    after = 0
    for i in range(10, 40):                      # 30 more good MPUs = 60 s
        after += len(mpu(st, base + i, psn))
        psn += 2

    print(f"  poisoned    : {before}/10 before, {after}/30 after one bit flip")
    print(f"  hi[{PID}]      : {hi_after}  (real seq was {base + 10})")
    print(f"  late_fragment {st.stats['late_fragment']}  "
          f"lost {st.stats['lost']}  lost_no_moof {st.stats['lost_no_moof']}")

    # ---- 3. the other half: a REAL counter reset must still be FOLLOWED ---
    # Refusing wild jumps is only half a fix.  A guard that also refuses a
    # genuine broadcaster restart would trade a rare permanent death for a
    # certain one, and a gate that tested only the poison case would call that
    # a pass -- the same shape as the AC-4 balance bug, where L was exemplary
    # and R was 20x too loud and three gates saw only L.
    st = S.MpuStreamer(ver_ext=0)
    psn = 0
    for i in range(10):
        mpu(st, base + i, psn)
        psn += 2
    new = 77_000                              # broadcaster restarts its count
    resync = 0
    for i in range(20):
        resync += len(mpu(st, new + i, psn))
        psn += 2
    print(f"  real reset  : {resync}/20 MPUs emitted after a genuine restart "
          f"(resyncs {st.stats['seq_resync']})")

    ok_poison = after == 30
    ok_reset = resync >= 20 - S.MpuStreamer.SEQ_RESYNC
    print()
    if not ok_poison:
        print(f"GATE FAIL -- asset is PERMANENTLY DEAD: {30 - after} of 30 "
              f"subsequent MPUs never emitted, from ONE flipped bit")
    if not ok_reset:
        print(f"GATE FAIL -- the guard also blocks a REAL counter reset: only "
              f"{resync}/20 recovered. Rejecting poison is not enough.")
    if ok_poison and ok_reset:
        print("GATE PASS -- survives a corrupted sequence number AND still "
              "follows a genuine counter reset")
    return 0 if (ok_poison and ok_reset) else 1


if __name__ == "__main__":
    sys.exit(main())
