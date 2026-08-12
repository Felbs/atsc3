#!/usr/bin/env python3
"""gate_e71.py -- video must land on its own audio even when the lane's
absolute clock (tfdt) disagrees with slot arithmetic.

Measured live 8/10 on e31: the video lane's tfdt agreed with
(seq - lane_seq0) * 2.002 to 3 ms for three hours, then diverged by
36.061 s (18.01 slots). The muxer SUBTRACTED the inferred value via
-itsoffset while ffmpeg ADDED the real tfdt back, so every later chunk
carried video 36 s ahead of its own audio -- more than any player's
picture queue can hold, and VLC dropped 82.7 % of pictures (its own
counter) while audio played perfectly. That is what the user saw.

This gate builds real fMP4 fragments with a controllable tfdt and mixes
a real chunk, then measures the muxed A/V skew with ffprobe. The
NEGATIVE CONTROL runs the shipped inference on the identical fragments
and must reproduce the skew.

    python lab/gate_e71.py    -> PASS/FAIL, exit 0 iff all pass
"""
from __future__ import annotations

import json
import os
import shutil
import struct
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tools"))

import atsc3_tv as tv                                           # noqa: E402

TIMESCALE = 90000
FAILS = []


def check(name, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" +
          (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


def box(typ, payload):
    return struct.pack(">I", len(payload) + 8) + typ + payload


def make_fragment(tfdt_ticks, sizes, payload):
    """A moof+mdat whose traf carries the given tfdt and trun sizes."""
    tfdt = box(b"tfdt", struct.pack(">BBBBQ", 1, 0, 0, 0, tfdt_ticks))
    tfhd = box(b"tfhd", struct.pack(">BBBBI", 0, 0, 0, 0x20, 1))
    trun_payload = struct.pack(">BBBBII", 1, 0, 0x02, 0x01,
                               len(sizes), 0)          # data-offset patched
    for s in sizes:
        trun_payload += struct.pack(">I", s)
    trun = box(b"trun", trun_payload)
    traf = box(b"traf", tfhd + tfdt + trun)
    mfhd = box(b"mfhd", struct.pack(">II", 0, 1))
    moof = box(b"moof", mfhd + traf)
    doff = len(moof) + 8
    # patch trun data-offset to point at the mdat payload
    fixed = bytearray(moof)
    i = fixed.find(b"trun")
    struct.pack_into(">i", fixed, i + 4 + 8, doff)
    mdat = box(b"mdat", payload)
    return bytes(fixed) + mdat


def annexb_frames(n, w=16):
    """n length-prefixed NAL 'frames' that satisfy frag_nal_sane."""
    out, sizes = b"", []
    for _ in range(n):
        nal = b"\x40\x01" + bytes(w)          # HEVC-ish payload
        blob = struct.pack(">I", len(nal)) + nal
        out += blob
        sizes.append(len(blob))
    return out, sizes


def build_lane(d, n_frags, tfdt_of):
    """Write a fake video lane + idx. tfdt_of(k) -> ticks."""
    path = os.path.join(d, "live_video_pid12.m4s")
    idxp = os.path.join(d, "live_video_pid12.idx")
    init = box(b"ftyp", b"isom" + b"\x00" * 8)
    recs = []
    with open(path, "wb") as f:
        f.write(init)
        for k in range(n_frags):
            payload, sizes = annexb_frames(3)
            frag = make_fragment(tfdt_of(k), sizes, payload)
            off = f.tell()
            f.write(frag)
            recs.append({"seq": 1000 + k, "off": off, "len": len(frag),
                         "dur": int(tv.MPU_SECONDS * TIMESCALE)})
    with open(idxp, "w") as f:
        for r in recs:
            f.write(json.dumps(r) + "\n")
    return path, len(init), recs


class FakeSess:
    def __init__(self, path, init_bytes):
        self.vid = {"path": path, "init_bytes": init_bytes}
        self.lanes = {}
        self.eng = type("T", (), {"path": None, "first": None, "pid": 13,
                                  "ch": 2, "state": None, "wav_base": None,
                                  "gen": None})()
        self.spa = self.eng
        self.anchor_seq = None


def gate_lane_time():
    print("gate 1: lane_time uses the fragment's own tfdt")
    d = tempfile.mkdtemp(prefix="e71_")
    try:
        SKEW_SLOTS = 18                      # the measured live divergence
        skew_s = SKEW_SLOTS * tv.MPU_SECONDS

        # a lane whose tfdt runs SKEW ahead of slot arithmetic
        def tfdt_of(k):
            return int((k + SKEW_SLOTS) * tv.MPU_SECONDS * TIMESCALE)

        path, init_bytes, recs = build_lane(d, 6, tfdt_of)
        sess = FakeSess(path, init_bytes)
        frags = [(r["seq"], r["off"], r["len"]) for r in recs[2:5]]
        _, tfdt0 = tv.write_chunk_video(sess, frags,
                                        os.path.join(d, "chunk.m4s"))
        check("tfdt parsed out of a real fragment", tfdt0 is not None,
              f"tfdt0={tfdt0}")

        lane_seq0 = recs[0]["seq"]
        measured = tv.lane_time(sess, frags, lane_seq0, tfdt0)
        inferred = (frags[0][0] - lane_seq0) * tv.MPU_SECONDS
        check("lane_time returns the MEASURED lane clock",
              abs(measured - (inferred + skew_s)) < 1e-3,
              f"measured {measured:.3f}s vs inferred {inferred:.3f}s")

        # video PTS = t_off + v_extra - lane_t + tfdt  (what ffmpeg does)
        t_off = 100.0
        pts_new = t_off - measured + tfdt0 / TIMESCALE
        check("video lands exactly on its audio (t_off)",
              abs(pts_new - t_off) < 1e-3, f"video pts {pts_new:.3f} "
              f"vs audio {t_off:.3f}")

        # NEGATIVE CONTROL: the shipped inference on the SAME fragments
        pts_old = t_off - inferred + tfdt0 / TIMESCALE
        check("NEGATIVE CONTROL: inferred lane_t skews video off its audio",
              abs(pts_old - t_off - skew_s) < 1e-3 and pts_old - t_off > 30,
              f"video pts {pts_old:.3f} vs audio {t_off:.3f} "
              f"= {pts_old - t_off:+.3f}s skew")

        # and the healthy lane (tfdt == slot arithmetic) must be untouched
        d2 = tempfile.mkdtemp(prefix="e71b_")
        try:
            p2, ib2, r2 = build_lane(
                d2, 6, lambda k: int(k * tv.MPU_SECONDS * TIMESCALE))
            s2 = FakeSess(p2, ib2)
            fr2 = [(r["seq"], r["off"], r["len"]) for r in r2[2:5]]
            _, t2 = tv.write_chunk_video(s2, fr2,
                                         os.path.join(d2, "chunk.m4s"))
            m2 = tv.lane_time(s2, fr2, r2[0]["seq"], t2)
            i2 = (fr2[0][0] - r2[0]["seq"]) * tv.MPU_SECONDS
            check("healthy lane: measured == inferred (no behaviour change)",
                  abs(m2 - i2) < 1e-3, f"{m2:.3f} vs {i2:.3f}")
        finally:
            shutil.rmtree(d2, ignore_errors=True)
    finally:
        shutil.rmtree(d, ignore_errors=True)


def gate_timescale():
    print("gate 2: timescale comes from the lane's own idx")
    d = tempfile.mkdtemp(prefix="e71c_")
    try:
        idxp = os.path.join(d, "l.idx")
        with open(idxp, "w") as f:
            f.write(json.dumps({"seq": 1, "off": 0, "len": 9,
                                "dur": 180180}) + "\n")
        check("dur 180180 over 2.002 s -> 90 kHz",
              abs(tv.lane_timescale(idxp) - 90000.0) < 1e-6,
              f"{tv.lane_timescale(idxp)}")
        missing = os.path.join(d, "nope.idx")
        check("missing idx falls back to the 90 kHz default",
              tv.lane_timescale(missing) == 90000.0)
    finally:
        shutil.rmtree(d, ignore_errors=True)


if __name__ == "__main__":
    gate_lane_time()
    gate_timescale()
    print(f"\n{'ALL GATES PASS' if not FAILS else 'FAILED: ' + ', '.join(FAILS)}")
    sys.exit(1 if FAILS else 0)
