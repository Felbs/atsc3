#!/usr/bin/env python3
"""race3.py -- the three-way race: OnAir vs ours-before vs ours-with-their-ideas.

Lanes, all on the SAME recording, the SAME machine, one after another with
nothing else running:

  A  OnAir            build/atsc3rx of ~/src/OnAir, 12 threads (their own decoder)
  B  ours, before     this repository at the commit before 2026-11-06 (a worktree),
                      with its original scalar kernels: PLP 0 only, linear channel
                      estimate, flooding LDPC
  C  ours, with ideas today's tree: --plp1 (M17 live Subframe 1 with the E104
                      estimate and residual-CFO tracking), the E107 layered LDPC,
                      the E106 vectorised demapper
  D  ours, ideas off  today's tree with the borrowed ideas switched off where a
                      switch exists: --no-plp1-comb, ATSC3_LDPC_LAYERED=0 (PLP 1
                      live itself is idea 2 and has no "off": D isolates ideas 1
                      and 3 on top of it)

Measures: wall time for the 20 s of air (lossless replay), PLP 0 and PLP 1 blocks,
WHUT 32.1 video frames and gaps, real-time headroom (a --realtime --shed replay),
and for C vs D the PLP 1 cliff with seeded noise (blocks converged per level).

Run: python lab/race3.py --head /path/to/worktree [--skip-sweep] [--json out.json]
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PY = sys.executable
CAP = os.path.expanduser("~/captures/race_2026-11-06/live_rf33_antB.cs16")
ONAIR = os.path.expanduser("~/src/OnAir/build/atsc3rx")
ROUTE = "239.255.32.1:8321"


def run(cmd, env=None, cwd=None, timeout=1800):
    e = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    if env:
        e.update(env)
    t = time.perf_counter()
    r = subprocess.run(cmd, env=e, cwd=cwd, capture_output=True, text=True, timeout=timeout)
    return r, time.perf_counter() - t


def video_stats(path):
    if not os.path.exists(path):
        return dict(frames=0, gaps=None, seconds=0.0)
    r = subprocess.run(["ffprobe", "-v", "error", "-select_streams", "v", "-show_entries",
                        "packet=pts_time", "-of", "csv=p=0", path], capture_output=True, text=True)
    pts = sorted(float(x.rstrip(",")) for x in r.stdout.split() if x.rstrip(","))
    gaps = sum(1 for a, b in zip(pts, pts[1:]) if b - a > 0.05)
    return dict(frames=len(pts), gaps=gaps, seconds=round(pts[-1] - pts[0], 2) if len(pts) > 1 else 0.0)


def lane_ours(root, name, extra_args, env, outdir, realtime=False):
    live = os.path.join(outdir, name + ("_rt" if realtime else ""))
    shutil.rmtree(live, ignore_errors=True)
    os.makedirs(live)
    js = live + ".json"
    cmd = [PY, "-m", "atsc3", "watch", "--capture", CAP, "--fmt", "cs16", "--rate", "6912000",
           "--accel", "cpu", "--assets", "all", "--player", "none", "--live-dir", live, "--json", js] + extra_args
    if realtime:
        cmd += ["--realtime", "--shed", "--decode-procs", "3"]
    r, secs = run(cmd, env=env, cwd=root)
    out = dict(lane=name, realtime=realtime, wall_s=round(secs, 1), rc=r.returncode)
    try:
        j = json.load(open(js))
        out.update(frames=j["frames"], plp0=f"{j['fec_converged']}/{j['fec_total']}",
                   plp0_ok=j["fec_converged"], plp0_n=j["fec_total"], x_rt=j["x_real_time"],
                   plp1=(f"{j['plp1']['fec_converged']}/{j['plp1']['fec_total']}" if j.get("plp1") else "not decoded"),
                   plp1_ok=(j["plp1"]["fec_converged"] if j.get("plp1") else None),
                   plp1_n=(j["plp1"]["fec_total"] if j.get("plp1") else None))
    except Exception as ex:                                     # noqa: BLE001
        out["error"] = f"{type(ex).__name__}: {ex}; stderr tail: {r.stderr[-400:]}"
    out["whut"] = video_stats(os.path.join(live, "live_route_video_pid110.m4s"))
    out["wjla"] = video_stats(os.path.join(live, "live_video_pid12.m4s"))
    return out


def lane_onair(outdir):
    ts = os.path.join(outdir, "onair.ts")
    r, secs = run([ONAIR, CAP, "--rate", "6.912e6", "--format", "cs16", "--out", ts], env={"ATSC3_THREADS": "12"})
    o = r.stdout
    m0 = re.search(r"PLP 0 over the run: FEC blocks (\d+), decoded (\d+)", o)
    m1 = re.search(r"PLP 1 over the run: FEC blocks (\d+), decoded (\d+)", o)
    out = dict(lane="A OnAir", wall_s=round(secs, 1), rc=r.returncode,
               plp0=(f"{m0.group(2)}/{m0.group(1)}" if m0 else "?"), plp1=(f"{m1.group(2)}/{m1.group(1)}" if m1 else "?"),
               plp0_ok=int(m0.group(2)) if m0 else None, plp0_n=int(m0.group(1)) if m0 else None,
               plp1_ok=int(m1.group(2)) if m1 else None, plp1_n=int(m1.group(1)) if m1 else None)
    mt = re.search(r"read ([\d.]+) s of signal in ([\d.]+) s \(([\d.]+)x real time\)", o)
    out["x_rt"] = float(mt.group(3)) if mt else None
    out["whut"] = video_stats(ts)
    out["wjla"] = dict(frames=0, gaps=None, seconds=0.0, note="no MMTP in OnAir")
    return out


def cliff_c_vs_d(frames_n=2, levels=(23.0, 22.0, 21.5, 21.0, 20.5)):
    """PLP 1 blocks converged per noise level, for the four combinations of the
    two switchable ideas (E104 estimate, E107 layered LDPC) on the M17 decoder."""
    sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
    import numpy as np
    import m6_cells as C
    import m9_fast
    import m13_sf1 as M13
    import m17_plp1_live as M17
    from m2_pilots import load
    from m10_core import fine_t0, frame_samples
    js = json.load(open(os.path.join(HERE, "m4_l1detail_hit_rf33.json")))
    g0 = C.Geometry.rf33_legacy()
    fd = m9_fast.FrameDecoder(threads=12, backend="cpu", iters=50, cpu_fast=True)
    fd.prewarm()
    decs = {c: M17.SubframeDecoder(fd, g0, js["l1_basic"], js["l1_detail"], 1, comb=c) for c in (True, False)}
    g = decs[True].g
    per = (frame_samples(g0) + (g.nfft + g.gi) * g.nsym) / M13.FS_POST
    frames = []
    for k in range(frames_n):
        y, _fs, _cfo, _h = load(CAP, 6.912e6, fmt="cs16", span_sec=per + 0.06, start_sec=10.0 + k * per)
        t0, _ = fine_t0(y, g0)
        frames.append((y, t0))
    psig = float(np.mean([np.mean(np.abs(y) ** 2) for y, _ in frames]))
    rows = []
    for snr in levels:
        row = dict(noise_db=snr)
        for comb in (True, False):
            for layered in (True, False):
                fd.ldpc_layered = layered
                conv = n = 0
                for k, (y, t0) in enumerate(frames):
                    rng = np.random.default_rng(3_000 + k * 97 + int(round(snr * 10)))
                    sn = np.sqrt(psig * 10 ** (-snr / 10) / 2)
                    yn = y + sn * (rng.standard_normal(len(y)) + 1j * rng.standard_normal(len(y)))
                    rep = {}
                    _pk, rep = decs[comb].decode(yn, t0, rep)
                    conv += rep["converged"]; n += rep["n_fec"]
                row[f"comb{int(comb)}_layer{int(layered)}"] = f"{conv}/{n}"
        rows.append(row)
        print("   ", row, flush=True)
    fd.ldpc_layered = True
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--head", required=True, help="worktree of the commit before the ideas (lane B)")
    ap.add_argument("--out", default=os.path.join(HERE, "e104_out", "race3"))
    ap.add_argument("--skip-sweep", action="store_true")
    ap.add_argument("--skip-realtime", action="store_true")
    ap.add_argument("--json", default=None)
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    res = {}
    print("RACE 3 -- OnAir vs ours-before vs ours-with-their-ideas, same air, same PC")
    print("=" * 80, flush=True)
    print("  lane A: OnAir", flush=True)
    res["A"] = lane_onair(a.out); print("   ", res["A"], flush=True)
    print("  lane B: ours before (worktree)", flush=True)
    res["B"] = lane_ours(a.head, "B_before", [], {"ATSC3_LDPC_LAYERED": "0"}, a.out); print("   ", res["B"], flush=True)
    print("  lane C: ours with their ideas", flush=True)
    res["C"] = lane_ours(ROOT, "C_ideas", ["--plp1", "--route", ROUTE], {"ATSC3_LDPC_LAYERED": "1"}, a.out); print("   ", res["C"], flush=True)
    print("  lane D: ours, ideas 1 and 3 off", flush=True)
    res["D"] = lane_ours(ROOT, "D_ideas_off", ["--plp1", "--no-plp1-comb", "--route", ROUTE], {"ATSC3_LDPC_LAYERED": "0"}, a.out); print("   ", res["D"], flush=True)
    if not a.skip_realtime:
        print("  real-time replays (3 decode processes, shedding allowed)", flush=True)
        res["B_rt"] = lane_ours(a.head, "B_before", [], {"ATSC3_LDPC_LAYERED": "0"}, a.out, realtime=True); print("   ", res["B_rt"], flush=True)
        res["C_rt"] = lane_ours(ROOT, "C_ideas", ["--plp1", "--route", ROUTE], {"ATSC3_LDPC_LAYERED": "1"}, a.out, realtime=True); print("   ", res["C_rt"], flush=True)
    if not a.skip_sweep:
        print("  PLP 1 cliff, C vs D (ideas 1 and 3 on/off), seeded noise", flush=True)
        res["cliff"] = cliff_c_vs_d()
    out = a.json or os.path.join(a.out, "race3.json")
    json.dump(res, open(out, "w"), indent=1)
    print(f"\n  wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
