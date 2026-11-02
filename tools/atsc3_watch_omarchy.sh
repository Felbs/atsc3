#!/bin/bash
# tools/atsc3_watch_omarchy.sh [RF] -- ATSC 3.0 (NextGen TV) on this Omarchy box: picture + 5.1 English (dialogue
# included) + Spanish + live captions, in VLC. Verified 2026-11-01 on RF33 / Antenna B / bias-T.
#
#   tools/atsc3_watch_omarchy.sh          RF33 on Antenna B (the yagi + LNA; bias-T is switched on and read back)
#   tools/atsc3_watch_omarchy.sh stop     stop everything gently (chain last, so the SDRplay service never wedges)
#
# What it starts (one of each, like tools/atsc3_tv_linux.sh on the old rig):
#   chain   python -m atsc3 watch ... --assets all   (owns the radio; writes the video/audio/subs lanes)
#   eng     atsc3_audio.py --pid 13 --element 5_X --channels 6
#                                                     (full 5.1 -- the centre channel carries the dialogue;
#                                                      2-channel mode drops it. 3x real time on this box.
#                                                      --pid 13 explicitly: without it the worker follows
#                                                      whichever lane locked first, once pid 14 = silence)
#   spa     atsc3_audio.py --pid 14 --element pair    (Spanish stereo simulcast)
#   subs    atsc3_subs.py                             (captions -> live.srt)
#   viewer  atsc3_tv.py --mode v2 --player vlc        (VLC: Audio menu = eng/spa, captions default ON)
#   tvwatch atsc3_tvwatch.py                          (restarts the viewer if it wedges)
set -u
RF=${1:-33}
ANT=${ATSC3_ANT:-"Antenna B"}
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
PY=${PY:-$HOME/Projects/.venv/bin/python}   # the shared venv (Arch blocks system pip); override with PY=
D=$HOME/atsc3_live_rf$RF
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1   # the AC-4 QMF loop thrashes a threaded BLAS (E56)
# bracket patterns everywhere: this script's own argv can never match itself (the pgrep trap)
HELPERS='atsc3_tvwat[c]h\.py|atsc3_t[v]\.py|atsc3_aud[i]o\.py|atsc3_su[b]s\.py'
CHAIN='atsc3 w[a]tch'

wait_gone() {  # wait_gone PATTERN SECONDS -> 0 if nothing matches
  for _ in $(seq 1 "$2"); do pgrep -f "$1" >/dev/null || return 0; sleep 1; done
  pgrep -f "$1" >/dev/null && return 1 || return 0
}

stop_all() {
  # helpers first: INT, then TERM, then KILL (python workers finish a pass on INT; don't wait forever)
  pkill -INT  -f "$HELPERS" 2>/dev/null; pkill -TERM -x vlc 2>/dev/null; pkill -TERM -x ffplay 2>/dev/null
  wait_gone "$HELPERS" 10 || { pkill -TERM -f "$HELPERS" 2>/dev/null; wait_gone "$HELPERS" 10; } \
                          || { pkill -KILL -f "$HELPERS" 2>/dev/null; sleep 1; }
  pkill -KILL -x vlc 2>/dev/null; pkill -KILL -x ffplay 2>/dev/null
  # the chain last, and gently: it is streaming from the RSPdx (INT = its own clean shutdown, TERM as backup)
  pkill -INT -f "$CHAIN" 2>/dev/null
  wait_gone "$CHAIN" 30 || { pkill -TERM -f "$CHAIN" 2>/dev/null; wait_gone "$CHAIN" 15; }
  left=$(pgrep -f "$HELPERS|$CHAIN" | wc -l)
  echo "stopped (left running: $left, vlc: $(pgrep -xc vlc))"
}

if [ "$RF" = "stop" ]; then stop_all; exit 0; fi

cd "$REPO" || exit 1
if pgrep -f "$HELPERS|$CHAIN" >/dev/null || pgrep -x vlc >/dev/null; then
  echo "something is already running -- stopping it first"; stop_all
fi
mkdir -p "$D"
# fresh audio/caption state: the viewer reads the wav channel count once, at start
rm -f "$D/live_audio.wav" "$D/live_audio.json" "$D/live_audio.state.json" \
      "$D/live_audio_spa.wav" "$D/live_subs.json" "$D/live.srt"

echo "chain: RF$RF on $ANT, bias-T on"
ATSC3_BIAST=1 systemd-run --user --scope -q -p MemoryMax=10G \
  "$PY" -m atsc3 watch --rf "$RF" --ant "$ANT" --rate 6912000 --accel cpu --assets all \
  --player none --secs 0 --live-dir "$D" --json "$D/report.json" > "$D/chain.log" 2>&1 &
CHAIN_PID=$!
# wait for the audio lane -- but give up if the chain refused (e.g. the meteor window) or died
until grep -qE 'live lane: audio pid 13|Traceback|NOT CONFIRMED|REFUSING' "$D/chain.log" 2>/dev/null; do
  kill -0 "$CHAIN_PID" 2>/dev/null || break
  sleep 2
done
grep -E 'bias-T|gain confirmed|acquired|NOT CONFIRMED|Traceback|REFUSING' "$D/chain.log" | head -4
if ! grep -q 'live lane: audio pid 13' "$D/chain.log"; then
  echo "chain did not start (see $D/chain.log) -- nothing else launched"; exit 1
fi

nohup "$PY" -u tools/atsc3_audio.py --live-dir "$D" --pid 13 --element 5_X --channels 6 \
  --interval 15 --wait 600 --start-behind 30 > "$D/aud_eng.log" 2>&1 &
nohup "$PY" -u tools/atsc3_audio.py --live-dir "$D" --pid 14 --element pair --channels 2 \
  --start-behind 30 --interval 15 --wait 600 --out "$D/live_audio_spa.wav" > "$D/aud_spa.log" 2>&1 &
nohup "$PY" -u tools/atsc3_subs.py --live-dir "$D" --interval 20 --wait 600 > "$D/subs.log" 2>&1 &
# the viewer decides stereo-vs-5.1 from the eng wav ONCE at start: wait until both wavs exist
echo "waiting for the audio decoders to bank the first seconds..."
until [ -s "$D/live_audio.wav" ] && [ -s "$D/live_audio_spa.wav" ]; do sleep 3; done
# ATSC3_TV_EXTRA="--chunk 6" / ATSC3_VLC_EXTRA="--no-drop-late-frames": A/B hooks for tools/atsc3_quality.py
# Unattended runs: with the monitor powered off (Hyprland DPMS) the compositor stops drawing VLC,
# ~96 % of pictures count as "lost" and the viewer pauses it -- a measurement of the dark screen,
# not the receiver. ~/.atsc3_headless_when_dark (a flag file) makes such runs use a headless VLC.
if [ -f "$HOME/.atsc3_headless_when_dark" ] && hyprctl monitors 2>/dev/null | grep -q 'dpmsStatus: 0' \
   && ! hyprctl monitors 2>/dev/null | grep -q 'dpmsStatus: 1'; then
  ATSC3_TV_EXTRA="${ATSC3_TV_EXTRA:-} --vlc-headless"; echo "display is off -> headless VLC"
fi
nohup "$PY" -u tools/atsc3_tv.py --live-dir "$D" --mode v2 --player vlc ${ATSC3_TV_EXTRA:-} > "$D/tv.log" 2>&1 &
nohup "$PY" -u tools/atsc3_tvwatch.py "$D" > "$D/tvwatch.log" 2>&1 &
echo "VIEW STACK UP on $D -- VLC opens in ~1-2 minutes (it banks a cushion first)."
echo "logs: $D/{chain,aud_eng,aud_spa,subs,tv,tvwatch}.log     stop: tools/atsc3_watch_omarchy.sh stop"
