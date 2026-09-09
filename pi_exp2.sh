#!/bin/sh
# pi_exp2.sh "<label>|<extra args>|<env assignments>" ...
# Same as pi_exp.sh but each config carries its own env (for the BLAS-pin test).
cd ~/atsc3 || exit 1
SECS=${SECS:-75}
PY=.venv/bin/python
OUT=data/pi_exp2_$(date +%H%M); mkdir -p "$OUT"
for spec in "$@"; do
    LABEL=$(echo "$spec" | cut -d'|' -f1)
    EXTRA=$(echo "$spec" | cut -d'|' -f2)
    ENVS=$(echo "$spec" | cut -d'|' -f3)
    D="$OUT/$LABEL"; mkdir -p "$D"
    echo "######## $LABEL : [$EXTRA] env[$ENVS]"
    env $ENVS $PY tools/atsc3_run.py --rf 33 --ant "Antenna B" --secs "$SECS" \
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
print(f"RESULT {label}: sustained {sum(tail)/len(tail):.3f}x  inst {sum(itail)/len(itail):.3f}x  "
      f"FEC {fp:.1f}%  reacq {reacq}")
if fe:   print("   fe :", fe[-1].strip()[:140])
if prof: print("   dec:", prof[-1].strip()[:165])
PYEOF
done
echo "######## BATCH DONE: $OUT"
