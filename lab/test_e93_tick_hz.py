"""E93 verification: the lane timeline must divide by the TRACK'S unit.

The bug: LiveWriter divided every trun dur by an assumed 90 kHz. foxlive's
ATEME encoder muxes at 240000 ticks/s, so every lane over-counted media_s
by exactly 240000/90000 = 8/3 = 2.667x - and all four lanes "agreed",
because one shared accounting path is one shared bug (E93, LIVE_TV_LAB.md
~8864). The fix stamps tick_hz at the producers and divides by it in the
accountant.

This test is entirely offline (banked captures only, zero radio):

  1. mdhd_timescale reads the truth from the air's own bytes:
     foxlive video init -> 240000, atsc3_live video init -> 90000.
     Both are read through the HARDENED trak selection (foxlive's init
     carries a second `hint` trak; "first trak" is a latent bug).
  2. The REAL LiveWriter.write(), fed synthetic segments:
     240 kHz seg with tick_hz stamped   -> 2.002 s   (the fix)
     90 kHz seg with tick_hz stamped    -> 2.002 s   (negative control A:
                                           the lanes that were always right
                                           must STAY right)
  3. NEGATIVE CONTROL B - the discriminator: a 240 kHz seg WITHOUT a
     tick_hz stamp must reproduce the BUG (5.3387 s) through the fallback,
     and must log the loud E93 warning. If this "fails" by reading 2.002,
     the fallback grew a guess and the test proves nothing.
  4. The media<=wall invariant against the banked foxlive live.json: the
     HISTORIC (pre-fix) numbers must VIOLATE it (2.6x). That the invariant
     fires on known-bad data is what makes it worth shipping.

Run:  python lab/test_e93_tick_hz.py        (exit 0 all green)
"""
import json
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data")
sys.path.insert(0, HERE)

import m7_play as PL                                          # noqa: E402
import m11_watch as MW                                        # noqa: E402

RESULTS = []
LOGS = []


class _StubST:
    @staticmethod
    def log(msg):
        LOGS.append(str(msg))


def rec(name, ok, detail=""):
    RESULTS.append((name, ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f" - {detail}" if detail else ""))


def read_init(path, cap=4 << 20):
    with open(path, "rb") as f:
        return f.read(cap)


def seg(pid, dur, tick_hz, handler=b"vide"):
    d = dict(seq=1, bytes=b"x" * 64, dur=dur, pid=pid, handler=handler,
             init=b"", samples=None, kept=None)
    if tick_hz is not None:
        d["tick_hz"] = tick_hz
    return d


def main():
    print("E93 tick_hz verification (banked data + real LiveWriter, no radio)")
    MW.ST = _StubST()                     # capture lane logs, no side effects

    # ---- 1. the truth from the air's own bytes -------------------------
    cases = [
        ("foxlive video init", os.path.join(DATA, "foxlive",
                                            "live_video_pid5017.m4s"), 240000),
        ("atsc3_live video init", os.path.join(DATA, "atsc3_live",
                                               "live_video_pid12.m4s"), 90000),
    ]
    for name, path, want in cases:
        if not os.path.exists(path):
            rec(f"mdhd_timescale({name})", False, f"missing {path}")
            continue
        got = PL.mdhd_timescale(read_init(path), b"vide")
        rec(f"mdhd_timescale({name}) == {want}", got == want, f"got {got}")

    # ---- 2 & 3. the real accountant ------------------------------------
    tmp = tempfile.mkdtemp(prefix="e93_")
    try:
        w = MW.LiveWriter(tmp)
        w.write(b"", seg(101, 480480, 240000))          # ATEME, stamped
        w.write(b"", seg(102, 180180, 90000))           # TS-4000, stamped
        w.write(b"", seg(103, 480480, None))            # unstamped = the bug
        m101 = w.lanes[101]["media"]
        m102 = w.lanes[102]["media"]
        m103 = w.lanes[103]["media"]
        rec("240 kHz stamped seg -> 2.002 s", abs(m101 - 2.002) < 1e-6,
            f"got {m101:.4f}")
        rec("90 kHz stamped seg -> 2.002 s (stays right)",
            abs(m102 - 2.002) < 1e-6, f"got {m102:.4f}")
        buggy = 480480 / 90000.0                        # 5.3387
        rec("UNSTAMPED 240 kHz seg reproduces the BUG (discriminator)",
            abs(m103 - buggy) < 1e-4, f"got {m103:.4f} want {buggy:.4f}")
        rec("unstamped lane logged the loud E93 warning",
            any("E93" in l and "tick_hz" in l for l in LOGS),
            f"logs={LOGS[-2:] if LOGS else []}")
        # the idx rows must carry the unit so no consumer infers it again
        idxf = os.path.splitext(w.lanes[101]["path"])[0] + ".idx"
        row = json.loads(open(idxf).read().splitlines()[0])
        rec("idx row carries tick_hz", row.get("tick_hz") == 240000,
            f"row={row}")
    finally:
        for ln in getattr(w, "lanes", {}).values():
            try:
                ln["fh"].close()
                ln["idx"].close()
            except Exception:
                pass
        shutil.rmtree(tmp, ignore_errors=True)

    # ---- 4. the invariant must FIRE on the banked pre-fix numbers ------
    fl = os.path.join(DATA, "foxlive", "live.json")
    if os.path.exists(fl):
        meta = json.load(open(fl))
        wall = float(meta["updated"]) - float(meta["started"])
        worst = max(ln.get("media_s", 0.0) for ln in meta["lanes"].values())
        rec("media<=wall VIOLATED on banked foxlive (invariant fires)",
            worst > wall * 1.02, f"wall {wall:.0f}s worst {worst:.0f}s "
            f"ratio {worst / wall:.2f}")
    else:
        rec("banked foxlive live.json present", False, fl)

    print("=" * 62)
    bad = [n for n, ok in RESULTS if not ok]
    print(f"{len(RESULTS) - len(bad)}/{len(RESULTS)} passed")
    for n in bad:
        print(f"  FAILED: {n}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
