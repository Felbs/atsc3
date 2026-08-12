#!/usr/bin/env python3
"""Relaunch atsc3_tv while no stop-file exists. The viewing stack's PID 1."""
import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def main():
    d = sys.argv[1] if len(sys.argv) > 1 else "data/e29"
    stop = os.path.join(ROOT, d, "tv_stop")
    log = open(os.path.join(ROOT, d, "tv.log"), "ab", buffering=0)
    while not os.path.exists(stop):
        p = subprocess.Popen(
            [sys.executable, "-u",
             os.path.join(ROOT, "tools", "atsc3_tv.py"),
             "--live-dir", d, "--lag", "60", "--subs", "burn"],
            cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
        while p.poll() is None:
            if os.path.exists(stop):
                p.terminate()
                return 0
            time.sleep(3)
        time.sleep(5)
    return 0

if __name__ == "__main__":
    sys.exit(main())
