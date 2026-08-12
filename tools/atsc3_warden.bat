@echo off
rem atsc3_warden.bat -- relaunch loop for the stack warden (E61).
rem The warden is crash-only by design; THIS is its supervisor. Single-
rem instance safety lives in the warden's own lock (a second copy exits
rem rc=3 immediately), so running this wrapper twice is harmless.
rem Launch detached:  Start-Process -WindowStyle Hidden tools\atsc3_warden.bat
cd /d "%~dp0.."
:loop
python -u tools\atsc3_warden.py --live-dir data\e31 --rf 33 >> data\e31\_warden\warden_boot.log 2>&1
rem rc 3 = another warden holds the lock; do not spin against it
if %errorlevel%==3 (
  timeout /t 300 /nobreak >nul
) else (
  timeout /t 20 /nobreak >nul
)
goto loop
