@echo off
REM Live ATSC 3.0 TV: English audio + closed captions (overlay). Close window to stop.
REM  - press  t  in the video window to toggle captions on/off,  a  to switch English 5.1 / Spanish
REM  - 5.1 too heavy on this PC? add  --stereo  to the command below
REM  - per-chunk health log: Z:\src\atsc3\lab\av_live\_tv\telemetry.jsonl
REM  - don't want captions? delete  --cc  below
REM  - Spanish audio? use the "Watch ATSC3 TV (Spanish)" shortcut
cd /d Z:\src\atsc3
"C:\Users\emane\radioconda\python.exe" tools\atsc3_watch_av.py --rf 33 --ant "Antenna B" --lang eng --cc
pause
