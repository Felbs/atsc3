@echo off
REM Live ATSC 3.0 TV with the SPANISH (SAP) audio track. Close window to stop.
cd /d Z:\src\atsc3
"%USERPROFILE%\radioconda\python.exe" tools\atsc3_watch_av.py --rf 33 --ant "Antenna B" --lang spa
pause
