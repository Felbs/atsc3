#!/bin/sh
# Test-bench stack (8/09): new CPU engine + DEEP cushion for stable video.
# lag 120 = the mux reads 2 min behind the decode head, so a band fade
# shorter than ~2 min never reaches the screen; lead 60 banks the player
# cushion on top. Total ~3 min behind air -- the price of smooth.
cd /z/src/atsc3
D=data/e31
mkdir -p $D  # NEVER rm (E38 law)
# E57 LAW (measured BOTH platforms): unpinned BLAS pools spin-wait the box --
# 100% CPU that pins to 4% with the SAME 1.02x decode. Always pin.
export OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 OMP_NUM_THREADS=4
PY="$(cygpath "$LOCALAPPDATA")/Programs/Python/Python312/python.exe"
python -u tools/atsc3_run.py --rf 33 --secs 43200 --max-restarts 200 --cooldown 20 --check 20 --startup 90 \
  --live-dir $D --extra "--assets all --raw-queue 600 --accel cpu --decode-procs 4 --threads 6 --fe-threads 12" \
  > $D/sup.log 2>&1 &
sleep 60
$PY -u tools/atsc3_audio.py --live-dir $D --interval 15 --wait 900 --channels 2 --start-behind 60 > $D/aud_eng.log 2>&1 &
$PY -u tools/atsc3_audio.py --live-dir $D --pid 14 --element pair --channels 2 --start-behind 60 --interval 15 --wait 900 --out $D/live_audio_spa.wav > $D/aud_spa.log 2>&1 &
$PY -u tools/atsc3_subs.py --live-dir $D --interval 20 --wait 900 > $D/subs.log 2>&1 &
sleep 120
$PY -u tools/atsc3_tv.py --live-dir $D --lag 120 --lead 60 > $D/tv.log 2>&1 &
$PY -u tools/atsc3_tvwatch.py $D > $D/tvwatch.log 2>&1 &
echo BENCH-STACK-UP
