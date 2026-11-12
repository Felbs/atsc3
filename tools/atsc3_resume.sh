#!/bin/sh
# Resume the full live ATSC 3.0 stack on RF33/data/e30. Called when the
# balloon receiver releases the radio. mkdir-only (E38 law); one-each workers.
cd /z/src/atsc3
D=data/e30
mkdir -p $D
rm -f $D/supervisor.lock
PY="${PY:-$LOCALAPPDATA/Programs/Python/Python312/python.exe}"
# 1) the chain hunter (radioconda python inside atsc3_run's find_python)
$PY -u tools/atsc3_run.py --rf 33 --secs 21600 --max-restarts 200 --cooldown 20 \
    --check 20 --startup 60 --live-dir $D --extra "--assets all --raw-queue 600" \
    > $D/sup_day.log 2>&1 &
sleep 50
# 2) workers: one each, near-head
$PY -u tools/atsc3_audio.py --live-dir $D --interval 15 --wait 600 --channels 2 --start-behind 120 > $D/aud_day.log 2>&1 &
$PY -u tools/atsc3_audio.py --live-dir $D --pid 14 --element pair --channels 2 --start-behind 120 --interval 15 --wait 600 --out $D/live_audio_spa.wav > $D/aud_spa_day.log 2>&1 &
$PY -u tools/atsc3_subs.py --live-dir $D --interval 20 --wait 600 > $D/sub_day.log 2>&1 &
# 3) viewer + wedge supervisor
$PY -u tools/atsc3_tv.py --live-dir $D > $D/tv_day.log 2>&1 &
$PY -u tools/atsc3_tvwatch.py $D > $D/tvwatch_day.log 2>&1 &
echo "RESUMED"
