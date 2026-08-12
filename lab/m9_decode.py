#!/usr/bin/env python3
"""M9 Step 4 -- the end-to-end driver, measured in units of REAL TIME.

Same job as `m7_route.py` (Frames -> IP datagrams with the boundaries kept,
written to a `.dg`) and the same ALP/IPv4 code, so the output is directly
diffable against M7's.  What changes is the decode: `m9_fast.FrameDecoder`
instead of `m6_payload.decode_frame`, with `--accel` selecting

    none   the untouched M6/M7 chain          (the oracle, and the fallback)
    cpu    m9_accel + threads, no CUDA
    gpu    m9_accel + threads + batched CUDA LDPC

One Frame is 1708032 samples at 6.912 Msps = 0.247111 s of air, so a run of N
Frames in W seconds of wall clock is N*0.247111/W times real time.  That is
the ONLY number this file reports as a headline; kernel speedups are printed
as a breakdown underneath it.

Usage:
    python m9_decode.py long_rf33.cs16 --rate 8e6 --frames 476 --accel gpu
"""
from __future__ import annotations

import argparse
import collections
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m6_cells as C                                             # noqa: E402
import m6_payload as P6                                          # noqa: E402
import m7_route as R7                                            # noqa: E402

FRAME_SEC = C.FRAME_SAMPLES / 6.912e6

_JOB = {}


def _init(path, rate, fmt, accel, threads, iters, dtype, cpu_fast=False):
    import m9_fast
    _JOB.update(path=path, rate=rate, fmt=fmt, accel=accel, threads=threads)
    if accel == "none":
        _JOB["fd"] = None
        return
    _JOB["fd"] = m9_fast.FrameDecoder(
        threads=threads, backend=("cpu" if accel == "cpu" else accel),
        iters=iters, ldpc_dtype=dtype, cpu_fast=cpu_fast)


def decode_chunk(start_sec, n_frames, timers=None):
    """-> (stream, bounds, per-frame diagnostics).  Same contract as m7."""
    import m9_fast
    tm = timers if timers is not None else collections.Counter()
    fd = _JOB["fd"]
    path, rate, fmt = _JOB["path"], _JOB["rate"], _JOB["fmt"]
    span = 0.05 + n_frames * FRAME_SEC
    t = time.perf_counter()
    if fd is None:
        y = C.load_frame(path, rate, fmt=fmt, start_sec=start_sec,
                         span_sec=span * (rate / 6.912e6))
    else:
        y = m9_fast.load_fast(path, rate, fmt=fmt, start_sec=start_sec,
                              span_sec=span * (rate / 6.912e6), ex=fd.ex)[0]
    tm["load"] += time.perf_counter() - t

    stream, bounds, frames = bytearray(), [], []
    t_prev = None
    for fi in range(n_frames):
        if (fi + 1) * C.FRAME_SAMPLES + 400000 > len(y):
            break
        centre = C.BOOTSTRAP if t_prev is None else t_prev + C.FRAME_SAMPLES
        t = time.perf_counter()
        if fd is None:
            t0, coh = C.fine_timing(y, span=20, centre=centre)
        else:
            t0, coh = fd.fine_timing(y, span=20, centre=centre)
        tm["fine_timing"] += time.perf_counter() - t
        t_prev = t0
        t = time.perf_counter()
        if fd is None:
            r16, pkts, diag = P6.decode_frame(y, t0, {})
        else:
            r16, pkts, diag = fd.decode_frame(y, t0)
        tm["decode_frame"] += time.perf_counter() - t
        good = [p for p in pkts if p is not None]
        for p in good:
            ptr, pay = P6.bb_split(p)
            if ptr is not None:
                bounds.append(len(stream) + ptr)
            stream += pay
        frames.append(dict(t0=int(t0), coh=round(float(coh), 4),
                           conv=diag["converged"], bch=diag["bch_ok"],
                           plp16=bool(r16["converged"])))
    return stream, bounds, frames, tm


def _job(args):
    s0, n = args
    try:
        stream, bounds, fr, tm = decode_chunk(s0, n)
    except Exception as e:                                     # noqa: BLE001
        import traceback
        return s0, None, None, None, f"{type(e).__name__}: {e}\n" \
                                     f"{traceback.format_exc()}"
    return s0, bytes(stream), bounds, fr, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("capture")
    ap.add_argument("--rate", type=float, required=True)
    ap.add_argument("--fmt")
    ap.add_argument("--frames", type=int, default=40)
    ap.add_argument("--chunk", type=int, default=16)
    ap.add_argument("--start", type=float, default=0.0)
    ap.add_argument("--out", default="m9.dg")
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--accel", default="gpu",
                    choices=("none", "cpu", "gpu", "gpu-full"))
    ap.add_argument("--iters", type=int, default=50)
    ap.add_argument("--ldpc-dtype", default="float64")
    ap.add_argument("--cpu-fast", action="store_true",
                    help="E52 float32 CPU fast path (default: exact float64, "
                         "which is what the m11 gate re-derives as reference)")
    ap.add_argument("--no-write", action="store_true",
                    help="throughput measurement only; skip the .dg")
    a = ap.parse_args()
    path = a.capture if os.path.isabs(a.capture) else os.path.join(HERE,
                                                                   a.capture)
    outp = a.out if os.path.isabs(a.out) else os.path.join(HERE, a.out)
    print(f"M9 -- ATSC 3.0 decode, accel={a.accel}, {a.jobs} process(es) x "
          f"{a.threads} threads")
    print("=" * 72)

    plan, k = [], 0
    while k < a.frames:
        n = min(a.chunk, a.frames - k)
        plan.append((a.start + k * FRAME_SEC, n))
        k += n

    t_start = time.perf_counter()
    tm = collections.Counter()
    if a.jobs > 1:
        import concurrent.futures as cf
        ex = cf.ProcessPoolExecutor(
            max_workers=a.jobs, initializer=_init,
            initargs=(path, a.rate, a.fmt, a.accel, a.threads, a.iters,
                      a.ldpc_dtype, a.cpu_fast))
        results = list(ex.map(_job, plan))
        ex.shutdown()
        t_init = 0.0
    else:
        t0 = time.perf_counter()
        _init(path, a.rate, a.fmt, a.accel, a.threads, a.iters, a.ldpc_dtype,
              a.cpu_fast)
        t_init = time.perf_counter() - t0
        results = []
        for s0, n in plan:
            stream, bounds, fr, _ = decode_chunk(s0, n, tm)
            results.append((s0, bytes(stream), bounds, fr, None))
    t_decode = time.perf_counter() - t_start

    # stitch, walk ALP, write datagrams -- byte-for-byte the M7 code
    t0 = time.perf_counter()
    big, bnd, frames = bytearray(), [], []
    for s0, stream, bounds, fr, err in results:
        if err:
            print(f"  chunk at {s0:.3f}s FAILED: {err}")
            continue
        if not fr:
            continue
        bnd += [b + len(big) for b in bounds]
        big += stream
        frames += fr
    alps, st = R7.AlpWalker(bytearray(big), bnd).run()
    ipr = R7.IpReasm()
    flows = collections.Counter()
    fh = None if a.no_write else open(outp, "wb")
    nip = 0
    for t, p in alps:
        if t != 0:
            continue
        for src, dst, sp, dp, pay in ipr.feed(p):
            if fh:
                fh.write(R7.REC.pack(src, dst, sp, dp, len(pay)) + pay)
            flows[(".".join(map(str, dst)), dp)] += 1
            nip += 1
    if fh:
        fh.close()
    t_transport = time.perf_counter() - t0
    wall = time.perf_counter() - t_start

    nfr = len(frames)
    air = nfr * FRAME_SEC
    conv = sum(f["conv"] for f in frames)
    bch = sum(f["bch"] for f in frames)
    print(f"\n  {nfr} Frames decoded, {conv}/{74*nfr} FEC Blocks converged, "
          f"{bch}/{74*nfr} BCH zero")
    print(f"  ALP {len(alps)} packets, {st['resync']} resyncs -> {nip} IP "
          f"datagrams over {len(flows)} flows")
    print(f"\n  === TIMING ===")
    print(f"    air decoded        {air:9.3f} s")
    print(f"    decode wall        {t_decode:9.3f} s   "
          f"({air/t_decode:6.3f}x real time)")
    print(f"    transport wall     {t_transport:9.3f} s")
    print(f"    TOTAL wall         {wall:9.3f} s   "
          f"**{air/wall:6.3f}x REAL TIME**")
    _fd = _JOB.get("fd")
    if _fd is not None and getattr(_fd, "tm", None):
        for k2, v in _fd.tm.items():
            tm["  ." + k2] += v
    if tm:
        print(f"\n    within decode (single process):")
        for k2 in sorted(tm, key=lambda k2: -tm[k2]):
            print(f"      {k2:16s} {tm[k2]:8.3f} s  "
                  f"{1000*tm[k2]/max(nfr,1):8.2f} ms/frame  "
                  f"{100*tm[k2]/t_decode:5.1f}%")
        print(f"      (decoder build/CUDA init {t_init:.2f} s, "
              f"not in the per-frame rows)")
    js = os.path.splitext(outp)[0] + "_summary.json"
    with open(js, "w") as f:
        json.dump(dict(capture=os.path.basename(path), accel=a.accel,
                       jobs=a.jobs, threads=a.threads, frames=nfr,
                       conv=conv, bch=bch, alp=len(alps), ip=nip,
                       resync=int(st["resync"]), air_s=air, wall_s=wall,
                       x_real_time=air / wall,
                       flows={f"{k2[0]}:{k2[1]}": v
                              for k2, v in flows.items()}), f, indent=1,
                  default=float)
    print(f"\n  wrote {js}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
