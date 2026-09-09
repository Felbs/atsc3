#!/bin/sh
# pi_exp.sh "<label>|<extra args>" ...   -- one live RF33 decode per config,
# same carrier/antenna/duration, printing sustained rate + stage profile so
# configs are comparable. One radio user at a time (serialized by the loop).
cd ~/atsc3 || exit 1
SECS=${SECS:-75}
PY=.venv/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
OUT=data/pi_exp_$(date +%H%M); mkdir -p "$OUT"
for spec in "$@"; do
    LABEL=${spec%%|*}; EXTRA=${spec#*|}
    D="$OUT/$LABEL"; mkdir -p "$D"
    echo "######## $LABEL : $EXTRA"
    $PY tools/atsc3_run.py --rf 33 --ant "Antenna B" --secs "$SECS" \
        --live-dir "$D" --extra "--assets all --accel cpu $EXTRA" \
        > "$D/run.out" 2>&1
    $PY - "$D" "$LABEL" <<'PYEOF'
import os, re, sys
d, label = sys.argv[1], sys.argv[2]
p = os.path.join(d, "chain.log")
log = open(p, errors="replace").read().splitlines() if os.path.exists(p) else []
rt = [float(m.group(1)) for l in log for m in [re.search(r"([0-9.]+)x rt", l)] if m]
inst = [float(m.group(1)) for l in log for m in [re.search(r"inst\s+([0-9.]+)x", l)] if m]
tail = rt[-6:] or [0]; itail = inst[-6:] or [0]
fecs = [(int(m.group(1)), int(m.group(2))) for l in log
        for m in [re.search(r"FEC (\d+)/(\d+)", l)] if m]
fp = 100.0*fecs[-1][0]/max(1,fecs[-1][1]) if fecs else float("nan")
reacq = sum(1 for l in log if "re-acquisition" in l or "re-acquiring" in l)
prof = [l for l in log if re.search(r"\bdecode\s+[0-9]", l)]
fe   = [l for l in log if "front end " in l and "acquisi" not in l]
pool = [l for l in log if "decode pool" in l or "cpu mode" in l]
print(f"RESULT {label}: sustained {sum(tail)/len(tail):.3f}x  inst {sum(itail)/len(itail):.3f}x  "
      f"FEC {fp:.1f}%  reacq {reacq}")
for l in pool[:2]: print("   cfg:", l.strip()[:120])
if fe:   print("   fe :", fe[-1].strip()[:150])
if prof: print("   dec:", prof[-1].strip()[:170])
PYEOF
done
echo "######## EXPERIMENTS DONE: $OUT"
