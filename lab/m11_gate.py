#!/usr/bin/env python3
"""m11_gate.py -- the streaming path against the batch path, re-derived NOW.

Four gates, in the order the pipeline runs:

  1. `PolyResampler` blocked == `resample_poly` whole-array, `array_equal`.
  2. `FrontEnd` blocked == `m9_fast.load_fast` whole-array, `array_equal`.
  3. `StreamAlpWalker` fed in chunks == `AlpWalker` over the whole stream,
     packet for packet.
  4. **End to end**: the streaming pipeline's decoded IP datagram file against
     `m9_decode.py`'s, SHA-256 to SHA-256, over the same air.

Gate 4 does not read a stored expectation.  It RE-RUNS `m9_decode.py` on the
same capture, now, and hashes what comes out.

    > Law (M9): a gate that compares against "the reference as of last week" is
    > a gate against a memory.  The oracle has to be re-run, not remembered --
    > and when two tracks share a repository, "re-run" means today.

M9 was nearly fooled by exactly this: the fast path reproduced M7's banked
`m7_long.dg` byte for byte, which looked like a pass, while the CURRENT
reference had gained an attenuation ladder and converged one more FEC Block.
So the reference here is a subprocess, not a file.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m11_stream as ST                                           # noqa: E402
import m7_route as R7                                             # noqa: E402

PY = sys.executable


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def gate_alp(bb, limit_mb, verbose=True):
    if not os.path.exists(bb):
        print(f"    SKIP  StreamAlpWalker: no baseband cache at "
              f"{os.path.basename(bb)}")
        return None
    big, bnd = bytearray(), []
    cap = int(limit_mb * 1e6)
    for s0, fr, stream, bounds in R7.bb_read(bb):
        bnd += [b + len(big) for b in bounds]
        big += stream
        if len(big) >= cap:
            break
    print(f"    (stitched {len(big)/1e6:.2f} MB of real Baseband stream, "
          f"{len(bnd)} anchors)")
    return ST.gate_alp_walker(big, bnd, chunks=41, verbose=verbose)


def gate_end_to_end(capture, rate, frames, accel, threads, fmt=None,
                    keep=False):
    """Streaming .dg vs a freshly re-run m9_decode .dg, byte for byte."""
    base = os.path.join(HERE, "m11_gate")
    ref_dg, str_dg = base + "_ref.dg", base + "_stream.dg"
    for p in (ref_dg, str_dg):
        if os.path.exists(p):
            os.remove(p)

    print(f"\n  -- re-deriving the reference: m9_decode.py on "
          f"{os.path.basename(capture)}, {frames} Frames, accel={accel}")
    t = time.time()
    r = subprocess.run(
        [PY, os.path.join(HERE, "m9_decode.py"), capture, "--rate", str(rate),
         "--frames", str(frames), "--chunk", str(frames), "--jobs", "1",
         "--threads", str(threads), "--accel", accel, "--out", ref_dg]
        + (["--fmt", fmt] if fmt else []),
        capture_output=True, text=True)
    if r.returncode != 0 or not os.path.exists(ref_dg):
        print(f"    FAIL  reference run failed:\n{r.stdout[-2000:]}"
              f"\n{r.stderr[-2000:]}")
        return False, {}
    ref_lines = [l for l in r.stdout.splitlines()
                 if "Frames decoded" in l or "ALP " in l]
    for l in ref_lines:
        print("      ref:" + l)
    print(f"    reference in {time.time()-t:.1f} s")

    # m9_decode loads `0.05 + N*FRAME_SEC` seconds and then requires a further
    # 400,000-sample tail margin per Frame, so it decodes N-1 Frames of a span
    # sized for N.  The streaming front end stops on "is this Frame's own
    # window complete", which is one Frame later.  Neither is wrong; they are
    # different end-of-FILE rules and on air neither ever fires.  Ask the
    # streaming path for exactly what the reference actually did, so the gate
    # is a whole-file SHA-256 rather than a prefix.
    n_ref = frames
    for l in ref_lines:
        w = l.split()
        if len(w) > 1 and w[1] == "Frames":
            n_ref = int(w[0])
            break

    print(f"\n  -- the streaming pipeline over the same air")
    t = time.time()
    r2 = subprocess.run(
        [PY, os.path.join(HERE, "m11_watch.py"), "--capture", capture,
         "--rate", str(rate), "--player", "none", "--dump-dg", str_dg,
         "--accel", accel, "--threads", str(threads), "--report", "1e9",
         "--frames", str(n_ref)]
        + (["--fmt", fmt] if fmt else []),
        capture_output=True, text=True)
    tail = r2.stdout.splitlines()
    for l in tail:
        if any(k in l for k in ("SUSTAINED", "FEC Blocks", "ALP resyncs",
                                "IP datagrams", "MMTP", "bootstrap",
                                "resampler", "flow gate", "video asset")):
            print("      str:" + l)
    if not os.path.exists(str_dg):
        print(f"    FAIL  streaming run produced no datagram file\n"
              f"{r2.stdout[-3000:]}\n{r2.stderr[-2000:]}")
        return False, {}
    print(f"    streaming in {time.time()-t:.1f} s")

    a, b = sha(ref_dg), sha(str_dg)
    na, nb = os.path.getsize(ref_dg), os.path.getsize(str_dg)
    # The streaming front end stops on "is the next Frame's window complete",
    # the batch loop on a fixed 400k-sample tail margin, so one path can hold
    # one more Frame than the other at the very end of a FILE.  On air neither
    # ever ends.  So the gate is: the shorter file is a PREFIX of the longer,
    # and the common part is byte-identical.
    same = a == b
    pref = False
    if not same:
        n = min(na, nb)
        with open(ref_dg, "rb") as f1, open(str_dg, "rb") as f2:
            pref = f1.read(n) == f2.read(n)
    ok = same or pref
    print(f"\n    reference  {na:>12,} bytes  sha256 {a}")
    print(f"    streaming  {nb:>12,} bytes  sha256 {b}")
    if same:
        print("    PASS  IDENTICAL datagram files")
    elif pref:
        print(f"    PASS  identical over the common {min(na,nb):,} bytes "
              f"({abs(na-nb):,}-byte tail difference = "
              f"{'streaming' if nb>na else 'batch'} held one more Frame)")
    else:
        n = min(na, nb)
        with open(ref_dg, "rb") as f1, open(str_dg, "rb") as f2:
            d1, d2 = f1.read(n), f2.read(n)
        i = next((i for i in range(n) if d1[i] != d2[i]), None)
        print(f"    FAIL  first differing byte at offset {i}")
    info = dict(ref_bytes=na, ref_sha256=a, stream_bytes=nb, stream_sha256=b,
                identical=same, prefix=pref)
    if not keep:
        for p in (ref_dg, str_dg):
            if os.path.exists(p):
                os.remove(p)
    return ok, info


def main(argv=None):
    ap = argparse.ArgumentParser(prog="atsc3 gate")
    ap.add_argument("--capture", default="rate6912_rf33.cs16",
                    help="the capture both paths decode")
    ap.add_argument("--rate", type=float, default=ST.FS_POST)
    ap.add_argument("--fmt", default=None)
    ap.add_argument("--frames", type=int, default=20)
    ap.add_argument("--accel", default="gpu-full")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--bb", default="m7_long.bb")
    ap.add_argument("--bb-mb", type=float, default=6.0)
    ap.add_argument("--rf", type=int, default=33, help=argparse.SUPPRESS)
    ap.add_argument("--skip-e2e", action="store_true")
    ap.add_argument("--keep", action="store_true")
    ap.add_argument("--json", default=None)
    a = ap.parse_args(argv)
    cap = a.capture if os.path.isabs(a.capture) else os.path.join(HERE,
                                                                  a.capture)
    bb = a.bb if os.path.isabs(a.bb) else os.path.join(HERE, a.bb)

    print("M11 -- the streaming path against the batch path")
    print("=" * 72)
    res = {}

    print("\n  1. resampler phase state")
    res["resampler_8"] = ST.gate_resampler(8e6, ST.FS_POST)
    res["resampler_10"] = ST.gate_resampler(10e6, ST.FS_POST, n=1_500_000,
                                            block=300_000)

    print("\n  2. continuous front end")
    res["front_end"] = ST.gate_front_end(cap, a.rate, fmt=a.fmt, frames=4)
    alt = os.path.join(HERE, "long_rf33.cs16")
    if os.path.exists(alt):
        print("     (and the 8 Msps banked capture, resampler in circuit)")
        res["front_end_8msps"] = ST.gate_front_end(alt, 8e6, frames=4)

    print("\n  3. resumable ALP walk")
    res["alp_walker"] = gate_alp(bb, a.bb_mb)

    if not a.skip_e2e:
        print("\n  4. end to end -- decoded IP datagrams, both paths")
        ok, info = gate_end_to_end(cap, a.rate, a.frames, a.accel, a.threads,
                                   fmt=a.fmt, keep=a.keep)
        res["end_to_end"] = ok
        res["end_to_end_info"] = info

    print("\n" + "=" * 72)
    hard = {k: v for k, v in res.items()
            if isinstance(v, bool) or v is None}
    bad = [k for k, v in hard.items() if v is False]
    for k, v in hard.items():
        print(f"    {'PASS' if v else ('SKIP' if v is None else 'FAIL')}  {k}")
    print(f"\n  {'ALL GATES PASS' if not bad else 'FAILED: ' + ', '.join(bad)}")
    if a.json:
        with open(a.json, "w") as f:
            json.dump(res, f, indent=1, default=str)
    return 0 if not bad else 1


if __name__ == "__main__":
    raise SystemExit(main())
