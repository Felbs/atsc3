#!/usr/bin/env python3
"""e73_gain_live.py -- LIVE gain A/B while the chain is failing. (E73)

The RF referee proved the fault is ours: HDHomeRun snq=100 lock=atsc3 on
the same antenna while our chain sits at 0.0% FEC and burns 6 cores to do
it. Signal strength rose (ss 65-68 -> 73) and the chain's notch went from
0.004 to 0.082 -- the front end looks compressed again, which is exactly
the condition the splitter relieved this morning.

Runs the REAL chain briefly at several front-end gains and reports the
FEC each one achieves. Scratch live-dir so the production lanes are never
touched. Nothing is adopted here -- this only measures.
"""
import json, os, re, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PY_RC = os.environ.get("ATSC3_RADIOCONDA",
                       os.path.join(os.path.expanduser("~"), "radioconda",
                                    "python.exe"))
SECS = float(sys.argv[1]) if len(sys.argv) > 1 else 75.0
GAINS = [int(x) for x in (sys.argv[2].split(",") if len(sys.argv) > 2
                          else ["2", "4", "6"])]
scratch = os.path.join(ROOT, "data", "e73_gain")
os.makedirs(scratch, exist_ok=True)          # mkdir-only

print(f"live gain A/B: {SECS:.0f}s per setting, gains {GAINS}", flush=True)
res = []
for g in GAINS:
    d = os.path.join(scratch, f"rf{g}")
    os.makedirs(d, exist_ok=True)
    cmd = [PY_RC, "-m", "atsc3", "watch", "--rf", "33", "--secs", str(SECS),
           "--player", "none", "--live-dir", d, "--accel", "cpu",
           "--decode-procs", "2", "--threads", "4", "--fe-threads", "6",
           "--assets", "none", "--rfgain", str(g)]
    t0 = time.time()
    r = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True,
                       timeout=SECS * 3 + 120)
    txt = (r.stdout or "") + (r.stderr or "")
    now = [float(x) for x in re.findall(r"\[\s*([\d.]+)% now\]", txt)]
    rates = [float(x) for x in re.findall(r"inst\s+([\d.]+)x", txt)]
    tail = now[-12:] if len(now) >= 12 else now
    row = dict(rfgain=g, n=len(now),
               mean_all=(sum(now) / len(now)) if now else 0.0,
               mean_tail=(sum(tail) / len(tail)) if tail else 0.0,
               best=max(now) if now else 0.0,
               inst=(sum(rates[-8:]) / len(rates[-8:])) if rates else 0.0)
    res.append(row)
    print(f"  rfgain_sel={g}: FEC tail-mean {row['mean_tail']:6.1f}%  "
          f"best {row['best']:6.1f}%  all-mean {row['mean_all']:6.1f}%  "
          f"inst {row['inst']:.2f}x  ({len(now)} samples, "
          f"{time.time()-t0:.0f}s)", flush=True)
json.dump(res, open(os.path.join(HERE, "e73_gain_live.json"), "w"), indent=1)
best = max(res, key=lambda r: r["mean_tail"]) if res else None
print()
if best:
    print(f"BEST: rfgain_sel={best['rfgain']} at {best['mean_tail']:.1f}% "
          f"(as-found rfgain_sel=2 got "
          f"{[r['mean_tail'] for r in res if r['rfgain']==2][0]:.1f}%)"
          if any(r['rfgain'] == 2 for r in res) else "")
