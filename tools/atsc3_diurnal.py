#!/usr/bin/env python3
"""atsc3_diurnal.py -- one record every 5 minutes, all night. (E68)

THE QUESTION THIS ANSWERS: did the splitter's +10.7 dB flatten RF33's
night floor? 8/09-8/10 measured 71.8% of 24 h decodable, with a nightly
floor ~19:56-22:10 and a pre-dawn floor. If RF33 now decodes THROUGH the
night, "stable video indefinitely" stops being antenna-blocked and the
whole campaign changes shape.

Deliberately the dumbest process in the stack: it NEVER touches the
radio, never kills anything, and only ever APPENDS. It cannot disturb the
soak it is measuring, which is the point -- the 8/10 lesson is that
instruments which act on the thing they measure end up measuring
themselves.

Each record carries BOTH receivers on the same antenna:
    ours    the chain's FEC 5-min mean, its instantaneous x-real-time,
            the video lane's generation and last_seq
    theirs  the HDHomeRun's ss/snq/seq on 107.1 (read-only, reserved tuner)
    viewer  whether the TS actually grew since the previous record
so a bad hour can be attributed to the air or to us, after the fact,
without anyone having been awake for it.

    python tools/atsc3_diurnal.py --live-dir data/e31
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

PERIOD = 300.0
HDHR_BIN = r"C:\Program Files\Silicondust\HDHomeRun\hdhomerun_config.exe"
HDHR_ID = os.environ.get("ATSC3_HDHR_ID", "")  # no default; discover if unset
HDHR_TUNER = os.environ.get("ATSC3_HDHR_TUNER", "1")


def fec_5min_mean(live_dir):
    """Mean of the chain's instantaneous FEC % over ~5 min, or None."""
    p = os.path.join(live_dir, "chain.log")
    try:
        with open(p, "rb") as f:
            f.seek(max(0, os.path.getsize(p) - 300_000))
            txt = f.read().decode("utf-8", "replace")
    except OSError:
        return None, None
    now = time.localtime()
    now_s = now.tm_hour * 3600 + now.tm_min * 60 + now.tm_sec
    vals, rates = [], []
    for m in re.finditer(
            r"\[(\d\d):(\d\d):(\d\d)\].*?inst\s+([\d.]+)x.*?\[\s*([\d.]+)% now\]",
            txt):
        t = int(m.group(1)) * 3600 + int(m.group(2)) * 60 + int(m.group(3))
        if (now_s - t) % 86400 <= 310:
            rates.append(float(m.group(4)))
            vals.append(float(m.group(5)))
    if not vals:
        return None, None
    return (sum(vals) / len(vals), sum(rates) / len(rates))


def hdhr():
    if not os.path.exists(HDHR_BIN):
        return None
    try:
        r = subprocess.run([HDHR_BIN, HDHR_ID, "get",
                            f"/tuner{HDHR_TUNER}/status"],
                           capture_output=True, text=True, timeout=10)
    except (OSError, subprocess.SubprocessError):
        return None
    m = re.search(r"ch=(\S+)\s+lock=(\S+)\s+ss=(\d+)\s+snq=(\d+)\s+seq=(\d+)",
                  r.stdout or "")
    if not m:
        return None
    return dict(ch=m.group(1), lock=m.group(2), ss=int(m.group(3)),
                snq=int(m.group(4)), seq=int(m.group(5)))


def lanes(live_dir):
    try:
        d = json.load(open(os.path.join(live_dir, "live.json")))
    except (OSError, ValueError):
        return None
    v = next((l for l in (d.get("lanes") or {}).values()
              if l.get("kind") == "video"), None)
    return dict(updated_age=round(time.time() - d.get("updated", 0), 1),
                gen=(v or {}).get("generation"),
                last_seq=(v or {}).get("last_seq"),
                bytes=d.get("bytes"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--live-dir", default="data/e31")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    if not os.path.isabs(a.live_dir):
        a.live_dir = os.path.join(ROOT, a.live_dir)
    out = a.out or os.path.join(a.live_dir, "_warden",
                                f"diurnal_{time.strftime('%m%d')}.jsonl")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    ts_p = os.path.join(a.live_dir, "_tv", "live_tv2.ts")
    prev_ts = None
    print(f"diurnal logger -> {out} (every {PERIOD/60:.0f} min, radio never "
          f"touched)", flush=True)
    while True:
        t0 = time.time()
        fec, rate = fec_5min_mean(a.live_dir)
        try:
            sz = os.path.getsize(ts_p)
        except OSError:
            sz = None
        rec = dict(t=time.time(), hhmm=time.strftime("%m-%d %H:%M"),
                   fec_5min=fec, inst_xrt=rate, hdhr=hdhr(),
                   lanes=lanes(a.live_dir), ts_bytes=sz,
                   ts_grew=(None if (sz is None or prev_ts is None)
                            else sz > prev_ts))
        prev_ts = sz if sz is not None else prev_ts
        try:
            with open(out, "a", encoding="utf-8") as f:
                f.write(json.dumps(rec) + "\n")
        except OSError as e:
            print(f"append failed: {e}", flush=True)
        h = rec["hdhr"] or {}
        print(f"{rec['hhmm']}  ours {('%.0f%%' % fec) if fec is not None else 'n/a':>5}"
              f"  {('%.2fx' % rate) if rate else '  -  '}"
              f"  referee snq={h.get('snq', '-')} ss={h.get('ss', '-')}"
              f"  ts_grew={rec['ts_grew']}", flush=True)
        time.sleep(max(5.0, PERIOD - (time.time() - t0)))


if __name__ == "__main__":
    sys.exit(main())
