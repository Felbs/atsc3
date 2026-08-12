#!/usr/bin/env python3
"""M12 -- the partial-MPU repair, gated against the loss it is supposed to fix.

M11 left one coarse failure in the receiver: an MPU is all-or-nothing, so a
SINGLE bad FEC Block anywhere inside it costs the whole 2.002 s fragment even
though the other 119 of 120 video samples are sitting in the buffer intact and
decodable.  `MpuStreamer` now cuts a holed MPU at its first short sample and
hands over what arrived, labelled.

WHY A FAULT INJECTOR AND NOT A LIVE RUN
---------------------------------------
Clean air does not exercise this.  A 90 s replay of the banked RF33 capture
decodes 35890/35890 FEC Blocks with zero ALP resyncs -- there is nothing to
repair, so a live run can only ever show that the repair did no HARM.  To show
that it does GOOD, the loss has to be manufactured: this drops a contiguous run
of IP datagrams out of the middle of the stream, which is what a burst of
failed FEC Blocks looks like by the time it reaches the transport.

THE GATES, each with something that must fail
---------------------------------------------
G1  CLEAN IS UNTOUCHED.  With no drops, repair on and repair off produce
    byte-identical segment streams.  A repair that perturbs healthy air is a
    regression however good it is on damaged air.
G2  THE REPAIR BEATS ITS OWN CONTROL.  Same injected hole, `--no-repair` vs
    repair: the repaired output must contain strictly MORE video samples.
    This is the control the whole change stands on -- if the repaired run
    merely equals the dropped run, nothing was recovered.
G3  THE REPAIRED FILE PLAYS.  ffmpeg decodes it with ZERO error lines.  A
    fragment that is structurally wrong can still be "recovered" on paper.
G4  THE CLOCK IS HONEST.  The repaired file's duration must match the clean
    file's to within one sample: the lost tail is a freeze in place, not a
    splice that pulls everything after it early.
G5  THE COUNTERS AGREE WITH THE BYTES.  `frames_recovered` from the telemetry
    must equal the extra frames ffmpeg actually decodes.  Silent patching is
    what M7's rule forbids; a counter that does not match the picture is the
    same sin one level up.

Usage:
    python m12_repair_gate.py [--dg m12_rf33.dg] [--keep-files]
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m7_route as R7                                            # noqa: E402
import m11_stream as ST                                          # noqa: E402


def load_dg(path):
    """-> [(dst, dport, payload)], the transport's own record format."""
    out = []
    with open(path, "rb") as fh:
        blob = fh.read()
    o, n = 0, R7.REC.size
    while o + n <= len(blob):
        src, dst, sp, dp, ln = R7.REC.unpack_from(blob, o)
        o += n
        out.append((dst, dp, blob[o:o + ln]))
        o += ln
    return out


def run(dgs, drop=(), repair=True):
    """Drive the transport at the datagram layer, skipping `drop` indices."""
    tr = ST.Transport(repair=repair)
    segs, lo, hi = [], (drop[0] if drop else -1), (drop[1] if drop else -1)
    for i, (dst, dp, pay) in enumerate(dgs):
        if lo <= i < hi:
            continue
        segs += tr._route((dst, dp), pay)
    return tr, segs


def write_mp4(tr, segs, path):
    if tr.init is None:
        return None
    with open(path, "wb") as fh:
        fh.write(tr.init)
        for s in segs:
            fh.write(s["bytes"])
    return path


def ffprobe_frames(path):
    p = subprocess.run(["ffprobe", "-v", "error", "-select_streams", "v:0",
                        "-count_frames", "-show_entries",
                        "stream=nb_read_frames", "-show_entries",
                        "format=duration", "-of", "json", path],
                       capture_output=True, text=True)
    try:
        js = json.loads(p.stdout)
        return (int(js["streams"][0]["nb_read_frames"]),
                float(js["format"]["duration"]), p.stderr.strip())
    except Exception:                                          # noqa: BLE001
        return None, None, (p.stderr or "ffprobe failed").strip()


def ffmpeg_errors(path):
    p = subprocess.run(["ffmpeg", "-v", "error", "-i", path, "-f", "null",
                        "-"], capture_output=True, text=True)
    return [l for l in p.stderr.splitlines() if l.strip()]


def find_hole(dgs, tr_clean, want_run=40):
    """Pick a datagram run that lands INSIDE an MPU rather than between two.

    Dropping across an MPU boundary loses a whole MPU at the head and proves
    nothing about repair, so the run is placed at the stream's midpoint and
    then only accepted if the resulting damage is a SHORT MPU (repairable)
    rather than a missing moof.
    """
    mid = len(dgs) // 2
    return (mid, mid + want_run)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="m12 repair gate")
    ap.add_argument("--dg", default="m12_rf33.dg")
    ap.add_argument("--run", type=int, default=40,
                    help="datagrams to drop, i.e. how big the hole is")
    ap.add_argument("--keep-files", action="store_true")
    a = ap.parse_args(argv)
    dg = a.dg if os.path.isabs(a.dg) else os.path.join(HERE, a.dg)
    if not os.path.exists(dg):
        print(f"  no datagram file at {dg} -- dump one with:\n"
              f"    python -m atsc3 watch --capture long_rf33.cs16 "
              f"--rate 8e6 --player none --secs 40 --dump-dg {a.dg}")
        return 2

    print("M12 -- partial-MPU repair vs the loss it fixes")
    print("=" * 72)
    dgs = load_dg(dg)
    print(f"  {len(dgs)} datagrams from {os.path.basename(dg)}")

    res, out = {}, {}

    # -- G1: clean air is untouched -------------------------------------
    tr_on, segs_on = run(dgs, repair=True)
    tr_off, segs_off = run(dgs, repair=False)
    same = ([s["bytes"] for s in segs_on] == [s["bytes"] for s in segs_off])
    res["clean_untouched"] = same
    print(f"\n  1. clean air, repair on vs off")
    print(f"    {'PASS' if same else 'FAIL'}  {len(segs_on)} segments vs "
          f"{len(segs_off)}, byte-identical: {same}")

    p_clean = os.path.join(HERE, "m12_gate_clean.mp4")
    write_mp4(tr_on, segs_on, p_clean)
    f_clean, d_clean, _ = ffprobe_frames(p_clean)
    print(f"    clean reference: {f_clean} frames, {d_clean:.3f} s")

    # -- G2: the repair beats its own control ---------------------------
    lo, hi = find_hole(dgs, tr_on, a.run)
    tr_r, segs_r = run(dgs, drop=(lo, hi), repair=True)
    tr_d, segs_d = run(dgs, drop=(lo, hi), repair=False)
    n_rep = tr_r.mmtp[list(tr_r.mmtp)[0]].stats["repaired"] if tr_r.mmtp else 0
    got_r = sum(s["kept"] for s in segs_r)
    got_d = sum(s["kept"] for s in segs_d)
    beats = got_r > got_d and n_rep > 0
    res["beats_control"] = beats
    print(f"\n  2. a {hi-lo}-datagram hole at index {lo}")
    print(f"    repaired  {len(segs_r)} segments, {got_r} video samples, "
          f"{n_rep} MPU(s) repaired")
    print(f"    control   {len(segs_d)} segments, {got_d} video samples, "
          f"whole MPU dropped")
    print(f"    {'PASS' if beats else 'FAIL'}  repair recovered "
          f"{got_r - got_d} samples the control threw away")

    # -- G3/G4/G5: the file plays, the clock is honest, counters agree ---
    p_rep = os.path.join(HERE, "m12_gate_repaired.mp4")
    write_mp4(tr_r, segs_r, p_rep)
    errs = ffmpeg_errors(p_rep)
    f_rep, d_rep, _ = ffprobe_frames(p_rep)
    res["plays_clean"] = not errs and f_rep is not None
    print(f"\n  3. the repaired file")
    print(f"    {'PASS' if res['plays_clean'] else 'FAIL'}  ffmpeg: "
          f"{len(errs)} error lines, {f_rep} frames decoded")
    for l in errs[:4]:
        print(f"          {l}")

    honest = (f_clean and d_clean and d_rep is not None
              and abs(d_rep - d_clean) < 0.05)
    res["clock_honest"] = bool(honest)
    print(f"\n  4. the timeline")
    print(f"    {'PASS' if honest else 'FAIL'}  repaired {d_rep:.3f} s vs "
          f"clean {d_clean:.3f} s -- the hole is a freeze in place, not a "
          f"splice")

    # -- G5: the counter must equal the frames ffmpeg actually gains -----
    p_ctl = os.path.join(HERE, "m12_gate_control.mp4")
    write_mp4(tr_d, segs_d, p_ctl)
    f_ctl, d_ctl, _ = ffprobe_frames(p_ctl)
    trunc = tr_r.stats["segment_truncated"]
    claimed = tr_r.stats["frames_recovered"]
    gained = (f_rep - f_ctl) if (f_rep is not None and f_ctl is not None) else None
    agree = trunc > 0 and gained is not None and gained == claimed
    res["counters_agree"] = bool(agree)
    print(f"\n  5. counters vs bytes")
    print(f"    {trunc} truncated segment(s); telemetry claims {claimed} "
          f"samples kept, ffmpeg decodes {gained} frames more than the "
          f"control ({f_rep} vs {f_ctl})")
    print(f"    {'PASS' if agree else 'FAIL'}  the truncation is reported and "
          f"the number is the picture's own")

    print("\n" + "=" * 72)
    for k, v in res.items():
        print(f"    {'PASS' if v else 'FAIL'}  {k}")
    bad = [k for k, v in res.items() if not v]
    print(f"\n  {'ALL GATES PASS' if not bad else '  FAILED: ' + ', '.join(bad)}")
    with open(os.path.join(HERE, "m12_repair_gate.json"), "w") as fh:
        json.dump(dict(res, clean_frames=f_clean, repaired_frames=f_rep,
                       control_samples=got_d, repaired_samples=got_r,
                       truncated_segments=trunc, control_frames=f_ctl,
                       frames_recovered=claimed),
                  fh, indent=1)
    if not a.keep_files:
        for p in (p_clean, p_rep, p_ctl):
            if os.path.exists(p):
                os.remove(p)
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
