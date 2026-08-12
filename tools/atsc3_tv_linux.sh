#!/bin/sh
# atsc3_tv_linux.sh -- the full live VIEWING stack on a Linux box (E56).
#
# Starts what atsc3_resume.sh starts on Windows, MINUS the chain: the
# chain owns the radio and is launched/supervised separately (it is
# already running when you want a picture). One-each discipline; every
# process nohup-detached with its log inside the live dir; mkdir-only
# (E38 law: run scripts create, NEVER clear).
#
#   tools/atsc3_tv_linux.sh [live-dir]      # default data/live1
#
# The viewer resolves the GUI environment (WAYLAND_DISPLAY / DISPLAY /
# XAUTHORITY, per-boot auth suffix) itself at VLC spawn time -- nothing
# display-related needs to be set here.
cd "$(dirname "$0")/.." || exit 1
D=${1:-data/live1}
PY=${PY:-.venv/bin/python}
mkdir -p "$D"
# MEASURED (E56, 8/09, 1600X): the QMF synthesis loop issues thousands of
# TINY matmuls; a default 12-thread OpenBLAS pool thrashes on them and the
# audio worker ran 0.09x real time (34 OS threads, 480% CPU, going
# nowhere). Single-threaded BLAS: 9.56x. The python-level channel
# threading is unaffected -- it is the BLAS pool that must not fan out.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
# one-each: refuse to double-start any stage (bracket patterns so this
# script's own argv can never self-match -- the pgrep trap)
for pat in 'atsc3_audio[.]py' 'atsc3_subs[.]py' 'atsc3_tv[.]py' \
           'atsc3_tvwatch[.]py'; do
    if pgrep -f "$pat" >/dev/null 2>&1; then
        echo "already running: $pat -- stop it first (one-each law)" >&2
        exit 1
    fi
done
# workers: eng (pid13 5.1 element, stereo out), spa (pid14 pair), subs
nohup "$PY" -u tools/atsc3_audio.py --live-dir "$D" --interval 15 \
    --wait 600 --channels 2 --start-behind 120 \
    > "$D/aud_eng.log" 2>&1 &
nohup "$PY" -u tools/atsc3_audio.py --live-dir "$D" --pid 14 \
    --element pair --channels 2 --start-behind 120 --interval 15 \
    --wait 600 --out "$D/live_audio_spa.wav" \
    > "$D/aud_spa.log" 2>&1 &
nohup "$PY" -u tools/atsc3_subs.py --live-dir "$D" --interval 20 \
    --wait 600 > "$D/subs.log" 2>&1 &
# viewer (spawns VLC itself once the cushion is banked) + wedge supervisor
nohup "$PY" -u tools/atsc3_tv.py --live-dir "$D" > "$D/tv.log" 2>&1 &
nohup "$PY" -u tools/atsc3_tvwatch.py "$D" > "$D/tvwatch.log" 2>&1 &
echo "VIEW STACK UP on $D: eng+spa audio workers, subs worker, viewer, tvwatch"
echo "logs: $D/{aud_eng,aud_spa,subs,tv,tvwatch}.log"
