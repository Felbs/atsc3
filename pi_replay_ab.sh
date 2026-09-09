#!/bin/sh
# pi_replay_ab.sh <label>  -- decode the SAME banked capture with the current
# code and report the stage profile. Deterministic input: differences between
# runs are the CODE, not the air.
cd ~/atsc3 || exit 1
LAB=${1:?label}
CAP=/home/felbs/atsc3/captures/rf33_bench.cs16
PY=.venv/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
D=data/replay_$LAB; mkdir -p "$D"
$PY -m atsc3 watch --capture "$CAP" --live-dir "$D" --player none \
    --accel cpu --decode-procs 1 --threads 4 --fe-threads 4 \
    --json "$D/run.json" > "$D/out.log" 2>&1
$PY - "$D" "$LAB" <<'PYEOF'
import os, re, sys
d, lab = sys.argv[1], sys.argv[2]
txt = ""
for f in ("out.log", "chain.log"):
    p = os.path.join(d, f)
    if os.path.exists(p):
        txt += open(p, errors="replace").read()
fe   = [l for l in txt.splitlines() if "front end " in l and "acquisi" not in l]
dec  = [l for l in txt.splitlines() if re.search(r"\bdecode\s+[0-9]", l)]
rt   = [float(m.group(1)) for l in txt.splitlines()
        for m in [re.search(r"([0-9.]+)x rt", l)] if m]
fec  = [(int(m.group(1)), int(m.group(2))) for l in txt.splitlines()
        for m in [re.search(r"FEC (?:Blocks\s+)?(\d+)/(\d+)", l)] if m]
tot  = fec[-1] if fec else (0, 0)
print(f"REPLAY {lab}: peak {max(rt) if rt else 0:.3f}x  "
      f"FEC {tot[0]}/{tot[1]}")
if fe:  print("   fe :", fe[-1].strip()[:150])
if dec: print("   dec:", dec[-1].strip()[:170])
PYEOF
