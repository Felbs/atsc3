#!/bin/bash
# atsc3_stress.sh -- decode every ATSC 3.0 carrier the survey found, both PLPs,
# and on the marginal ones run the borrowed ideas OFF as well (E104 estimate,
# E107 layered LDPC) so the comparison is on the same air minutes apart.
#
#   tools/atsc3_stress.sh SURVEY.json OUTDIR [SECS] [ANT]
#
# For each RF with a bootstrap: watch --rf N --plp1 --assets all --player none
# for SECS, with the per-Frame PLP 1 log; then, if the carrier is marginal
# (bootstrap peak ratio under 40, or the ideas-on run lost blocks, shed Frames
# or re-acquired), the same again with ATSC3_LDPC_LAYERED=0 --no-plp1-comb.
# Single-tenant radio: stop the TV first (~/watch_tv_both.sh stop).
set -u
SURVEY=$1; OUT=$2; SECS=${3:-45}; ANT=${4:-Antenna B}
REPO=${REPO:-$(cd "$(dirname "$0")/.." && pwd)}   # the checkout this script lives in
PY=$HOME/Projects/.venv/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ATSC3_BIAST=1 ATSC3_METEOR_WINDOW=off ATSC3_PLP1_LOG=1
mkdir -p "$OUT"
cd "$REPO" || exit 1

RFS=$($PY -c '
import json, sys
j = json.load(open(sys.argv[1]))
hits = [(int(r["rf"]), float(r.get("peak_ratio", 0))) for r in j["results"] if r.get("atsc3")]
print(" ".join(f"{rf}:{pr:.0f}" for rf, pr in sorted(hits)))
' "$SURVEY")
echo "carriers with a bootstrap: $RFS"

clean_on_run() {   # exit 0 when the ideas-on run of RF $1 was strong and lost nothing
  $PY -c '
import json, sys
try:
    j = json.load(open(sys.argv[1]))
except Exception:
    sys.exit(1)
p = j.get("plp1") or {}
clean = (j["fec_converged"] == j["fec_total"] and p.get("fec_converged") == p.get("fec_total")
         and not j["pipeline"].get("frame_dropped", 0) and not j["front_end"].get("reacquire", 0))
sys.exit(0 if clean else 1)
' "$OUT/rf$1_on.json"
}

summarise() { $PY tools/atsc3_stress_summary.py "$OUT" 2>/dev/null | awk -v n="$(basename "$1" .json)" '$0==n{p=1;print;next} /^rf/{p=0} p'; }

for item in $RFS; do
  rf=${item%%:*}; pr=${item##*:}
  for mode in on off; do
    if [ $mode = off ] && [ "$pr" -ge 40 ] && clean_on_run $rf; then
      echo "== RF$rf: strong and clean with the ideas on; skipping the ideas-off run"; continue
    fi
    d=$OUT/rf${rf}_$mode; rm -rf "$d"; mkdir -p "$d"
    if [ $mode = off ]; then export ATSC3_LDPC_LAYERED=0; EXTRA="--no-plp1-comb"; else export ATSC3_LDPC_LAYERED=1; EXTRA=""; fi
    echo "== RF$rf ideas $mode (bootstrap peak $pr)"
    timeout $((SECS + 120)) $PY -m atsc3 watch --rf $rf --ant "$ANT" --rate 6912000 --accel cpu --assets all \
      --plp1 $EXTRA --player none --secs $SECS --live-dir "$d" --json "$OUT/rf${rf}_$mode.json" > "$OUT/rf${rf}_$mode.log" 2>&1
    summarise "$OUT/rf${rf}_$mode.json" "$OUT/rf${rf}_$mode.log"
    for f in "$d"/live*video*.m4s; do [ -s "$f" ] || continue; printf '   %s: ' "$(basename $f)"; ffprobe -v error -select_streams v:0 -show_entries packet=pts_time -of csv=p=0 "$f" 2>/dev/null | sort -n | awk 'NR>1 && $1-p>0.05{g++} {p=$1} NR==1{f=$1} END{printf "%d frames, %.1f s, %d gaps\n", NR, p-f, g}'; done
  done
done
echo "stress done: $OUT"
