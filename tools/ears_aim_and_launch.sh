#!/bin/sh
# ears_aim_and_launch.sh -- aim the antenna on RF8 with live feedback, and
# open the TV by itself when a position finally decodes.
#
# 9/08 rewrite. The RF33-on-Antenna-A version this replaces had two flaws
# that made it nearly useless for actual aiming:
#   * it ran a FULL decode probe per position (~75 s at FEC 0, because LDPC
#     iterates hardest on the worst signal), so you could try four positions
#     in five minutes;
#   * it compared a SINGLE peak ratio per position, and peak ratio here is
#     noisy enough (29.2-48.3 back-to-back at identical settings) that one
#     sample is indistinguishable from luck.
# tools/atsc3_aim.py fixes both: bootstrap-only probes (~3 s, SDR held open)
# and a ROLLING MEDIAN with the spread shown. It spends a real decode probe
# only once the median says it is worth one.
#
# Target is RF8 (VHF-high, 183 MHz): at Smith Point it beat both UHF
# carriers roughly 2:1 (peak 55 vs 27/24). Note it is VHF -- the element
# length and orientation that peak it are NOT the ones you used for UHF.
# Override with: ears_aim_and_launch.sh --rf 33 --ifgr 32
cd "$(dirname "$0")/.." || exit 1
exec .venv/bin/python tools/atsc3_aim.py \
    --rf 8 --ant "Antenna B" --rfgain 2 --ifgr 42 \
    --live-dir data/tv_rf8 "$@"
