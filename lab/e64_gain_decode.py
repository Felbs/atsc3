#!/usr/bin/env python3
"""e64_gain_decode.py -- decode every gain-sweep capture and rank by FEC.

Offline: no radio. Each capture is replayed through the REAL chain
(`atsc3 watch --capture`) so the metric is FEC block convergence, not rms.
Runs at reduced parallelism and BelowNormal priority so the live chain
keeps its delivering window.
"""
import json, os, re, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, "e64_gain")
PY_RC = os.environ.get(
    "ATSC3_RADIOCONDA",
    os.path.join(os.path.expanduser("~"), "radioconda", "python.exe"))
man = json.load(open(os.path.join(OUT, "manifest.json")))
runs = [r for r in man["runs"] if r.get("i") != "ingress"]

def decode(path, tag):
    scratch = os.path.join(OUT, f"live_{tag}")
    os.makedirs(scratch, exist_ok=True)          # mkdir-only
    cmd = [PY_RC, "-m", "atsc3", "watch", "--rf", str(man["rf"]),
           "--capture", path, "--player", "none", "--live-dir", scratch,
           "--accel", "cpu", "--decode-procs", "1", "--threads", "2",
           "--assets", "none"]
    env = dict(os.environ)
    for k in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
        env[k] = "2"
    flags = subprocess.BELOW_NORMAL_PRIORITY_CLASS if os.name == "nt" else 0
    try:
        r = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True,
                           env=env, creationflags=flags,
                           timeout=900)  # pipe-ok: TimeoutExpired salvaged below
    except subprocess.TimeoutExpired as e:  # keep the evidence (7/31 balloon law)
        r = subprocess.CompletedProcess(cmd, -1, e.stdout or "", e.stderr or "")
    txt = (r.stdout or "") + (r.returncode and (r.stderr or "") or "")
    # cumulative FEC counter from the status lines: FEC good/total
    tot = re.findall(r"FEC (\d+)/(\d+)", txt)
    now = [float(x) for x in re.findall(r"\[\s*([\d.]+)% now\]", txt)]
    good, total = (int(tot[-1][0]), int(tot[-1][1])) if tot else (0, 0)
    frames = re.findall(r"(\d+) Frames", txt)
    return dict(good=good, total=total,
                pct=(100.0 * good / total if total else 0.0),
                now_mean=(sum(now) / len(now) if now else 0.0),
                now_max=(max(now) if now else 0.0),
                n_samples=len(now), frames=int(frames[-1]) if frames else 0)

print(f"{'set':>3} {'rfgain':>6} {'IFGR':>5} {'FEC good/total':>18} "
      f"{'cum%':>7} {'inst-mean%':>11} {'frames':>7}")
res = []
for r in runs:
    tag = f"s{r['i']}_rf{r['rfgain_sel']}_if{r['ifgr']}"
    t0 = time.time()
    d = decode(r["path"], tag)
    d.update(i=r["i"], rfgain_sel=r["rfgain_sel"], ifgr=r["ifgr"])
    res.append(d)
    print(f"{r['i']:>3} {r['rfgain_sel']:>6} {r['ifgr']:>5} "
          f"{d['good']:>8}/{d['total']:<9} {d['pct']:>6.1f}% "
          f"{d['now_mean']:>10.1f}% {d['frames']:>7}   "
          f"({time.time()-t0:.0f}s)", flush=True)
json.dump(res, open(os.path.join(OUT, "decode_results.json"), "w"), indent=1)

base = [d for d in res if d["rfgain_sel"] == 2 and d["ifgr"] == 32]
print("\n=== DRIFT CONTROL (baseline captured FIRST and LAST) ===")
for d in base:
    print(f"  set {d['i']}: cum {d['pct']:.1f}%  inst-mean {d['now_mean']:.1f}%")
if len(base) >= 2:
    spread = abs(base[0]["now_mean"] - base[-1]["now_mean"])
    print(f"  spread {spread:.1f} points -> "
          f"{'USABLE (channel held still)' if spread < 15 else 'VOID: the channel drifted during the sweep'}")
    ref = (base[0]["now_mean"] + base[-1]["now_mean"]) / 2.0
    print(f"\n=== RANKED vs baseline mean {ref:.1f}% ===")
    for d in sorted(res, key=lambda x: -x["now_mean"]):
        tag = "  <-- as-found" if (d["rfgain_sel"], d["ifgr"]) == (2, 32) else ""
        print(f"  rfgain_sel={d['rfgain_sel']} IFGR={d['ifgr']:2d}  "
              f"inst-mean {d['now_mean']:6.1f}%  ({d['now_mean']-ref:+6.1f}){tag}")
