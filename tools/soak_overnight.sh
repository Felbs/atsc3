#!/bin/bash
# Overnight soak of the whole live stack (13.27): start now, close the window at STOP_AT (HH:MM),
# sample the stack's working set every 5 minutes, leave the telemetry and logs under LIVE_DIR.
# Scheduled by Task Scheduler ("atsc3 overnight soak") between the wxTuna sonde hunts
# (18:40-20:20 and 06:40-08:20 hold the RSPdx): 20:30 -> 06:30.
#   bash tools/soak_overnight.sh [STOP_AT] [ANT] [LIVE_DIR]
set -u
unset PYTHONPATH   # 13.27: a PothosSDR py3.9 PYTHONPATH breaks radioconda's SoapySDR import
STOP_AT=${1:-06:30}
ANT=${2:-Antenna A}
D=${3:-/c/captures/soak_night}
cd "$(dirname "$0")/.." || exit 1
export OPENBLAS_NUM_THREADS=1 ATSC3_LDPC_BOXPLUS=1 ATSC3_METEOR_WINDOW=off
export ATSC3_WORKER_SPREAD=8   # 13.29: one demod worker per core complex; refereed on replay, bytes identical, +33 %
PY=${PY:-$HOME/radioconda/python.exe}
LOG=$D.log; MEM=${D}_mem.txt
rm -rf "$D"; mkdir -p "$(dirname "$D")"
echo "soak start $(date '+%F %T') ant '$ANT' stop at $STOP_AT -> $D" > "$MEM"
# RF and CHAIN_EXTRA may come from the environment (13.30: the same script stress-tests RF33)
RF=${RF:-25}
CHAIN_EXTRA=${CHAIN_EXTRA:---accel cpu --plp1 --ldm-procs 8 --ldm-proc-threads 2 --threads 4 --raw-queue 600}
"$PY" tools/atsc3_watch_av.py --rf "$RF" --ant "$ANT" --cc --live-dir "$D" \
  --chain-extra "$CHAIN_EXTRA" > "$LOG" 2>&1 &
ORCH=$!
stop_epoch=$(date -d "today $STOP_AT" +%s); now=$(date +%s)
[ "$stop_epoch" -le "$now" ] && stop_epoch=$(date -d "tomorrow $STOP_AT" +%s)
while [ "$(date +%s)" -lt "$stop_epoch" ]; do
  if ! kill -0 $ORCH 2>/dev/null; then echo "orchestrator exited early at $(date +%T)" >> "$MEM"; break; fi
  sleep 300
  tot=$(powershell -NoProfile -Command "(Get-CimInstance Win32_Process | Where-Object { \$_.CommandLine -match 'atsc3_|atsc3 watch|ffplay|multiprocessing' -and \$_.CommandLine -notmatch 'Get-CimInstance' } | ForEach-Object { (Get-Process -Id \$_.ProcessId -ErrorAction SilentlyContinue).WorkingSet64 } | Measure-Object -Sum).Sum / 1MB")
  echo "$(date +%T) stack ${tot%.*} MB" >> "$MEM"
done
# close the window the way a person does: ffplay's pid from the viewer's JSON sidecar
if [ -f "$D/_tv/player.pid" ]; then
  pp=$("$PY" -c "import json,sys; print(json.load(open(sys.argv[1]))['pid'])" "$D/_tv/player.pid" 2>/dev/null)
  echo "closing ffplay pid $pp at $(date +%T)" >> "$MEM"
  [ -n "$pp" ] && powershell -NoProfile -Command "Stop-Process -Id $pp -Force -Confirm:\$false -ErrorAction SilentlyContinue"
fi
for i in $(seq 1 30); do kill -0 $ORCH 2>/dev/null || break; sleep 2; done
if kill -0 $ORCH 2>/dev/null; then echo "orchestrator still up 60 s after the close -- killing" >> "$MEM"; kill $ORCH; fi
sleep 3
left=$(powershell -NoProfile -Command "(Get-CimInstance Win32_Process | Where-Object { \$_.CommandLine -match 'atsc3_|atsc3 watch|ffplay' -and \$_.CommandLine -notmatch 'Get-CimInstance' }).Count")
echo "soak end $(date '+%F %T') leftovers $left" >> "$MEM"
