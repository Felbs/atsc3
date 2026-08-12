#!/usr/bin/env python3
"""E82 -- drive the LIVE chain's front end over an LDM/CTI multiplex.

This is `m11_watch` minus the radio and minus the player: the real
`m11_stream.FrontEnd` (bootstrap lock, continuous de-rotation, notch decision,
tracked fine timing) feeding the real `m44_ldm.LdmPipeline` (geometry from L1,
batched Frame demod, streaming CTI, phase tracker, 64K LDPC).

There is NO `--start-frame`, no `--l1 file.json`, and no channel number
anywhere in the decode path.  The only inputs are a capture and a sample rate;
everything else -- FFT size, guard interval, symbol count, PLP layout, which
layer is the Core, the CTI depth and the commutator phase -- is read off the
air.  That is the requirement E82 was set: a phase search that only works
offline, with a human choosing the start Frame, is not done.

    python e82_ldm_run.py ../data/fox25/rf25_fox_g4_0810_2254.cs16 \
        --rate 6912000 --frames 200 --dump-dg e82_out/fox.dg
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m11_stream as ST                                            # noqa: E402
import m44_ldm as M44                                              # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("capture")
    ap.add_argument("--rate", type=float, default=6.912e6)
    ap.add_argument("--fmt", default=None)
    ap.add_argument("--start", type=float, default=0.0)
    ap.add_argument("--frames", type=int, default=40)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--fe-threads", type=int, default=4)
    ap.add_argument("--iters", type=int, default=50)
    ap.add_argument("--accel", default="cpu", choices=("cpu", "gpu"))
    ap.add_argument("--decode-procs", type=int, default=0,
                    help="E85: demodulate Frames in N worker PROCESSES.  The "
                         "Frame demod is per-Frame independent and was 230 ms "
                         "of a 242 ms budget; the CTI and the phase are not "
                         "and stay in the parent.  0 = the serial path.")
    ap.add_argument("--proc-threads", type=int, default=4)
    ap.add_argument("--fec-threads", type=int, default=None)
    ap.add_argument("--block", type=int, default=854_000)
    ap.add_argument("--dump-dg", default=None)
    ap.add_argument("--live-dir", default=None)
    ap.add_argument("--route", default=None)
    ap.add_argument("--assets", default="all")
    ap.add_argument("--json", default=None)
    ap.add_argument("--warmup", type=int, default=0,
                    help="Frames excluded from the SUSTAINED figure.  Table "
                         "building, the L1 probe and the CTI's own fill are "
                         "one-off tune-in costs; a 'sustained' rate that "
                         "amortises them over a short run measures the run "
                         "length, not the chain.")
    ap.add_argument("--force-phase", type=int, default=None,
                    help="CONTROL: override the acquired start_row with this "
                         "value (a wrong phase must produce nothing)")
    ap.add_argument("--no-cti", action="store_true",
                    help="CONTROL: bypass the CTI de-interleaver entirely")
    a = ap.parse_args(argv)

    from concurrent.futures import ThreadPoolExecutor
    ex = ThreadPoolExecutor(max_workers=a.fe_threads)
    pipe = M44.LdmPipeline(iters=a.iters, threads=a.threads, ex=ex,
                           accel=a.accel, log=ST.log,
                           procs=a.decode_procs,
                           proc_threads=a.proc_threads,
                           fec_threads=a.fec_threads)
    if a.no_cti:
        M44.CTI.out_index = lambda i, n, s: np.asarray(i, np.int64)
        ST.log("  CONTROL: CTI BYPASSED")

    def plan_cb(w, t0_local, ps=None):
        # the Frame start is not yet known to sub-sample accuracy, so scan the
        # Preamble's own coherence plateau the way m3_preamble.analyse does
        best = None
        for t in range(t0_local - 40, t0_local + 8):
            c = _preamble_coh_raw(w, t, ps)
            if best is None or c > best[1]:
                best = (t, c)
        L = M44.l1_from_window(w, best[0], ps=ps)
        if L is None or not L.get("ok"):
            return None
        plan = M44.LdmPlan.from_l1_result(L, label="live")
        if not plan.uses_cti:
            ST.log("  L1 says this multiplex does NOT use the CTI -- "
                   "the m9_fast path owns it, not this one")
        pipe.adopt(plan)
        pipe.prewarm()
        return plan

    fe = ST.FrontEnd(a.rate, ex=ex, fast=True, plan_cb=plan_cb)
    src = ST.FileSource(a.capture, a.rate, fmt=a.fmt, block=a.block,
                        realtime=False, start=a.start)
    dg_fh = open(a.dump_dg, "wb") if a.dump_dg else None
    want = {"video": (b"vide",), "av": (b"vide", b"soun"),
            "all": (b"vide", b"soun", b"subt", b"text")}[a.assets]
    tr = ST.Transport(probe=64, want=want, dg_fh=dg_fh,
                      route=([a.route] if a.route else None))
    live = None
    if a.live_dir:
        import m11_watch as W11
        live = W11.LiveWriter(a.live_dir)

    src.start()
    t_start = time.time()
    t_warm = None
    n_warm = 0
    nseg = 0
    nframe = 0
    t_fe = 0.0
    try:
        while nframe < a.frames:
            x = src.read()
            if x is None:
                break
            if isinstance(x, bytes):
                fe.reacquire()
                pipe.reset()
                continue
            t = time.perf_counter()
            fe.push(x)
            frames = fe.frames()
            t_fe += time.perf_counter() - t
            for idx, w, t0, coh in frames:
                if nframe >= a.frames:
                    break
                if a.force_phase is not None and pipe.cti is not None \
                        and pipe.cti.origin_frame is not None:
                    pipe.cti.start_row = a.force_phase % pipe.plan.nrows
                for stream, bounds in pipe.push_frame(idx, w, t0):
                    for seg in tr.feed(stream, bounds):
                        nseg += 1
                        if live is not None:
                            live.write(tr.init, seg)
                nframe += 1
                if a.warmup and nframe == a.warmup:
                    t_warm = time.time()
                    n_warm = nframe
                    for d in (pipe.tm, pipe.fd.tm):
                        d.clear()
                    t_fe = 0.0
                if nframe % 25 == 0:
                    wall = time.time() - t_start
                    air = nframe * (pipe.plan.frame_sec if pipe.plan else 0)
                    ST.log(f"  {nframe:4d} Frames  {air / max(wall, 1e-9):5.3f}x rt"
                           f"  FEC {pipe.n_conv}/{pipe.n_blocks}"
                           f"  BCH {pipe.n_bch}  phase {pipe.ph.state}"
                           f"  seg {nseg}")
    except KeyboardInterrupt:
        ST.log("  interrupted")
    finally:
        src.close()
        # E85: drain the pool before the transport, or the tail Frames are
        # decoded and then thrown away -- the same shed the m11 e2e gate
        # caught on the RF33 path as a short datagram file.
        try:
            for stream, bounds in pipe.flush():
                for seg in tr.feed(stream, bounds):
                    nseg += 1
                    if live is not None:
                        live.write(tr.init, seg)
        except Exception as e:                                 # noqa: BLE001
            ST.log(f"  pool flush failed: {type(e).__name__}: {e}")
        for seg in tr.flush():
            nseg += 1
            if live is not None:
                live.write(tr.init, seg)
        if dg_fh:
            dg_fh.close()
        if live is not None:
            live.close()
    wall = time.time() - t_start
    s = pipe.summary()
    s["wall_s"] = round(wall, 3)
    s["air_s"] = round(nframe * (pipe.plan.frame_sec if pipe.plan else 0), 3)
    s["x_real_time"] = round(s["air_s"] / wall, 4) if wall else 0
    nmeas = nframe - n_warm
    if t_warm is not None and nmeas > 0:
        wall_s = time.time() - t_warm
        s["sustained_x_real_time"] = round(
            nmeas * pipe.plan.frame_sec / wall_s, 4)
        s["sustained_frames"] = nmeas
        s["sustained_wall_s"] = round(wall_s, 3)
        s["warmup_frames"] = n_warm
    denom = max(nmeas if t_warm is not None else nframe, 1)
    s["front_end_ms_per_frame"] = round(1000 * t_fe / denom, 2)
    # the stage timers were cleared at the warm-up mark, so they must be
    # divided by the MEASURED Frames, not by every Frame the run saw
    s["ms_per_frame"] = {k: round(1000 * v / denom, 2)
                         for k, v in pipe.tm.items()}
    s["demod_ms_per_frame"] = {k: round(1000 * v / denom, 2)
                               for k, v in pipe.fd.tm.items()}
    s["front_end"] = dict(fe.stats)
    s["ft"] = dict(fe.ft.stats)
    s["transport"] = dict(tr.stats)
    s["segments"] = nseg
    s["frame_budget_ms"] = round(1000 * (pipe.plan.frame_sec
                                         if pipe.plan else 0), 2)
    print("\n  === E82 LDM RUN ===")
    print(f"    Frames             {s['frames']}  ({s['air_s']:.3f} s of air "
          f"over {wall:.3f} s of wall)")
    print(f"    whole run          {s['x_real_time']:.4f}x real time "
          f"(includes tune-in)")
    if "sustained_x_real_time" in s:
        print(f"    SUSTAINED          {s['sustained_x_real_time']:.4f}x REAL "
              f"TIME  ({s['sustained_frames']} Frames after a "
              f"{s['warmup_frames']}-Frame warm-up, "
              f"{s['sustained_wall_s']:.1f} s wall)")
    print(f"    CTI phase          {s['state']}  start_row {s['start_row']} "
          f"C {s['C']}  origin Frame {s['origin_frame']}")
    print(f"    FEC Blocks         {s['converged']}/{s['blocks']} LDPC "
          f"converged, {s['bch_zero']} BCH ZERO "
          f"({100.0 * s['bch_zero'] / max(s['blocks'], 1):.1f}%)")
    print(f"    segments           {nseg}")
    print(f"    where the time goes, ms/Frame (budget "
          f"{s['frame_budget_ms']} ms):")
    print(f"      front end  {s['front_end_ms_per_frame']:7.1f}")
    for k, v in sorted(s["ms_per_frame"].items(), key=lambda kv: -kv[1]):
        print(f"      {k:10s} {v:7.1f}")
    for k, v in sorted(s["demod_ms_per_frame"].items(), key=lambda kv: -kv[1]):
        print(f"        (demod) {k:8s} {v:7.1f}")
    if "pool" in s:
        print(f"    demod pool         {s['pool']['procs']} procs x "
              f"{s['pool']['threads']} threads, BLAS {s['pool']['blas']}  "
              f"{s['pool']['stats']}")
        print(f"      worker fec   {s.get('worker_fec_ms_per_frame',0):7.1f} ms/Frame of WORKER wall")
        print(f"      worker demod {s['worker_demod_ms_per_frame']:7.1f} ms/Frame "
              f"of WORKER wall (summed over all of them -- the cost the pool "
              f"HIDES, not the cost the Frame pays)")
    print(f"    phase stats        {s['phase']}")
    print(f"    pipeline stats     {s['stats']}")
    print(f"    cti stats          {s['cti']}")
    if a.dump_dg and os.path.exists(a.dump_dg):
        h = hashlib.sha256(open(a.dump_dg, "rb").read()).hexdigest()
        s["dg"] = dict(bytes=os.path.getsize(a.dump_dg), sha256=h)
        print(f"    datagram file      {os.path.getsize(a.dump_dg)} bytes  "
              f"sha256 {h[:32]}...")
    if a.json:
        json.dump(s, open(a.json, "w"), indent=1, default=float)
        print(f"  wrote {a.json}")
    pipe.close()
    return 0


_PLANLESS = None


def _preamble_coh_raw(y, t0, ps=None, _c={}):
    """Preamble pilot coherence before any geometry is known.

    A/322 Table H.1.1 keys FFT/GI/DX off the bootstrap's own
    `preamble_structure`, and 7.2.5.1 fixes the first Preamble symbol at the
    minimum NoC -- so this metric exists BEFORE L1, which is exactly when
    acquisition needs it.  The menu is swept and the best row wins.
    """
    import m3_spec as S3
    from m3_preamble import pilot_coherence
    best = -1.0
    menu = (S3.PREAMBLE_STRUCTURE if ps is None
            else {ps: S3.PREAMBLE_STRUCTURE[ps]})
    for _ps, Pt in menu.items():
        nfft, gi, dx = Pt["fft"], Pt["gi"], Pt["dx"]
        w = t0 + gi
        if w < 0 or w + nfft > len(y):
            continue
        try:
            Y = np.fft.fftshift(np.fft.fft(y[w:w + nfft]))
            best = max(best, pilot_coherence(Y, nfft, gi, dx, 4, 0))
        except KeyError:
            continue
    return best


if __name__ == "__main__":
    raise SystemExit(main())
