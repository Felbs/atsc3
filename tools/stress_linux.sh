#!/bin/sh
# stress_ubuntu.sh -- mirror of the 8/22 Windows stress test, for the Ubuntu rig.
#   usage: lab/stress_ubuntu.sh "<antenna port>" [secs-per-carrier]
# Runs: doctor --radio, a bootstrap survey over the known carriers, then a
# headless decode per carrier through the SAME levers the launcher uses
# (BLAS pinned, accel autodetected -> cpu), and prints one line per carrier:
# sustained rt factor, FEC %, MPUs. mkdir-only: never clears anything (E38).
cd "$(dirname "$0")/.." || exit 1
ANT=${1:?antenna port, e.g. "Antenna B"}
SECS=${2:-150}
PY=${PY:-.venv/bin/python}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
OUT=data/stress_$(date +%m%d)
mkdir -p "$OUT"
echo "== governor: $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor)"
echo "== doctor (with radio)"
$PY tools/atsc3_doctor.py --radio 2>&1 | grep -E "FAIL|radio|SDR|SUMMARY"
echo "== survey on $ANT (bootstrap peak per carrier)"
# atsc3_survey takes a RANGE (lo-hi), not a list -- one call per carrier
for C in 8 25 29 30 33; do
    $PY tools/atsc3_survey.py --channels "$C-$C" --secs 6 --ant "$ANT" \
        --out "$OUT/survey_rf$C.json" 2>&1 | grep -vE "^\s*$" | tail -2
done
for RF in ${CARRIERS:-33 25 30 8 29}; do
    D="$OUT/rf${RF}_live"; mkdir -p "$D"
    echo "== RF$RF decode ${SECS}s -> $D"
    $PY tools/atsc3_run.py --rf "$RF" --ant "$ANT" --secs "$SECS" \
        --live-dir "$D" --extra "--assets all" > "$D/run.out" 2>&1
    tail -2 "$D/chain.log" 2>/dev/null | cut -c1-150
    $PY - "$D" <<'EOF'
import json, os, sys, re
d = sys.argv[1]
log = open(os.path.join(d, "chain.log"), errors="replace").read().splitlines() if os.path.exists(os.path.join(d, "chain.log")) else []
rt = [float(m.group(1)) for l in log for m in [re.search(r"([0-9.]+)x rt", l)] if m]
fec = [l for l in log if "FEC" in l]
last = fec[-1] if fec else "(no FEC line)"
m = re.search(r"FEC (?:Blocks\s+)?(\d+)/(\d+)", last)
hr = [l for l in log if "slowest stage" in l]
if hr:
    last = last + "  " + hr[-1].strip()
pct = (100.0 * int(m.group(1)) / max(1, int(m.group(2)))) if m else float("nan")
vid = [f for f in os.listdir(d) if "video" in f and f.endswith(".m4s")]
sz = sum(os.path.getsize(os.path.join(d, f)) for f in vid)
print(f"   RESULT RF{os.path.basename(d)[2:4]}: sustained {sum(rt[-12:])/max(1,len(rt[-12:])):.2f}x  FEC {pct:.2f}%  video lanes {len(vid)} ({sz/1e6:.0f} MB)")
EOF
done
echo "== done: $OUT"
