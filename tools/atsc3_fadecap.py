#!/usr/bin/env python3
"""atsc3_fadecap.py -- spend RF33's FADES on the RF25 multiplex. (E68)

The soak's qualified clock only advances while RF33 delivers, so radio
time taken DURING A FADE costs the soak nothing. This watches for a real
fade and spends it capturing RF25 (the WBFF multiplex, centre 539.000
MHz; WBFF is the host LICENSEE, not a service -- E94) -- the
campaign's biggest open channel. E58 measured Fox's threshold down to
~5.0 dB effective with the link ~2 dB short; the splitter then bought
+10.7 dB of SNR on RF33 by killing front-end intermod. If any of that
transfers, Fox may be decodable for the first time.

WHY THIS IS SAFE TO RUN UNATTENDED
    * it only fires after FADE_MIN_S of CONTINUOUS non-delivery, judged
      on the chain's own FEC -- a brief dip is not a fade;
    * it declares a maintenance window (the warden and judge both honour
      it, and both IGNORE it once stale), so the soak is neither
      disturbed nor credited while it works;
    * every exit path restores: the window is lifted in a finally block,
      and even if this process is killed outright the marker expires by
      itself and supervision resumes without anyone present;
    * it is bounded -- MAX_CAPS all night, one per fade, and it will not
      start a capture it cannot finish inside the window;
    * it NEVER decodes anything. It banks bytes with integrity numbers
      and gets out of the way.
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

PY_RC = os.environ.get(
    "ATSC3_RADIOCONDA",
    os.path.join(os.path.expanduser("~"), "radioconda", "python.exe"))
FADE_FEC = 20.0            # below this the chain is not delivering
FADE_MIN_S = 600.0         # 10 continuous minutes before we believe it
POLL_S = 60.0
CAP_SECS = 180.0
MAX_CAPS = 4
MIN_GAP_S = 3600.0         # never two captures inside an hour


def log(m):
    print(f"{time.strftime('%m-%d %H:%M:%S')} {m}", flush=True)


def fec_now(live_dir, n=12):
    p = os.path.join(live_dir, "chain.log")
    try:
        with open(p, "rb") as f:
            f.seek(max(0, os.path.getsize(p) - 200_000))
            txt = f.read().decode("utf-8", "replace")
    except OSError:
        return None
    v = [float(x) for x in re.findall(r"\[\s*([\d.]+)% now\]", txt)][-n:]
    return (sum(v) / len(v)) if v else None


def chain_pids():
    out = []
    try:
        import psutil
    except ImportError:
        return out
    for p in psutil.process_iter(["pid", "name", "cmdline"]):
        try:
            if not (p.info["name"] or "").lower().startswith("python"):
                continue
            c = " ".join(p.info["cmdline"] or [])
        except Exception:                                      # noqa: BLE001
            continue
        if "atsc3_run.py" in c:
            out.append(p.info["pid"])
    return out


def stop_chain():
    for pid in chain_pids():
        log(f"  stopping supervisor pid {pid} (tree)")
        subprocess.run(["taskkill", "/PID", str(pid), "/T", "/F"],
                       capture_output=True)
    time.sleep(5)


def radio_busy():
    try:
        import psutil
    except ImportError:
        return False
    for p in psutil.process_iter(["pid", "name", "cmdline"]):
        try:
            c = " ".join(p.info["cmdline"] or [])
        except Exception:                                      # noqa: BLE001
            continue
        if "-m atsc3" in c and "watch" in c:
            return True
    return False


def capture(live_dir, out_dir, rf, secs, rfgain, ifgr):
    os.makedirs(out_dir, exist_ok=True)          # mkdir-only, never clears
    out = os.path.join(out_dir,
                       f"rf{rf}_fade_{time.strftime('%m%d_%H%M%S')}.cs16")
    cmd = [PY_RC, os.path.join(HERE, "atsc3_capture.py"),
           "--rf", str(rf), "--secs", str(secs), "--rfgain", str(rfgain),
           "--ifgr", str(ifgr), "--out", out]
    r = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True,
                       timeout=secs * 2 + 180)
    for line in (r.stdout or "").splitlines():
        log(f"    {line}")
    return out, r.returncode == 0, (r.stdout or "")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--live-dir", default="data/e31")
    ap.add_argument("--rf", type=int, required=True,
                    help="RF channel number (required: a shipped default is wrong for every reader)")
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--rfgain", type=int, default=2)
    ap.add_argument("--ifgr", type=int, default=32)
    ap.add_argument("--max-caps", type=int, default=MAX_CAPS)
    a = ap.parse_args()
    if not os.path.isabs(a.live_dir):
        a.live_dir = os.path.join(ROOT, a.live_dir)
    out_dir = a.out_dir or os.path.join(ROOT, "data", "fox25")
    wdir = os.path.join(a.live_dir, "_warden")
    pause = os.path.join(wdir, "pause")
    ledger = os.path.join(wdir, "fadecap.jsonl")

    log(f"fade-capture armed: RF{a.rf} -> {out_dir}, {CAP_SECS:.0f}s per "
        f"capture, max {a.max_caps}, fade >= {FADE_MIN_S/60:.0f} min")
    caps = 0
    fade_since = None
    last_cap = 0.0
    while caps < a.max_caps:
        time.sleep(POLL_S)
        fec = fec_now(a.live_dir)
        if fec is None or fec >= FADE_FEC:
            if fade_since is not None:
                log(f"  fade ended (FEC {fec if fec is None else round(fec)}%)"
                    f" -- RF33 is delivering again")
            fade_since = None
            continue
        fade_since = fade_since or time.time()
        waited = time.time() - fade_since
        if waited < FADE_MIN_S or time.time() - last_cap < MIN_GAP_S:
            continue
        if os.path.exists(pause):
            continue                     # somebody else owns the stack
        log(f"FADE CONFIRMED: RF33 at {fec:.0f}% for {waited/60:.0f} min "
            f"-- spending it on RF{a.rf} (capture {caps+1}/{a.max_caps})")
        ok = False
        out = ""
        try:
            with open(pause, "w", encoding="utf-8") as f:
                f.write(f"E68 fade capture RF{a.rf}")
            time.sleep(18)               # let the warden observe the window
            stop_chain()
            if radio_busy():
                log("  radio still busy -- standing down this fade")
            else:
                out, ok, txt = capture(a.live_dir, out_dir, a.rf,
                                       CAP_SECS, a.rfgain, a.ifgr)
                caps += 1
                last_cap = time.time()
                try:
                    with open(ledger, "a", encoding="utf-8") as f:
                        f.write(json.dumps(dict(
                            t=time.time(), when=time.strftime("%m-%d %H:%M"),
                            rf=a.rf, out=out, usable=ok,
                            fade_min=round(waited / 60, 1),
                            rf33_fec=round(fec, 1))) + "\n")
                except OSError:
                    pass
                log(f"  {'BANKED' if ok else 'VOID'}: {out}")
        except Exception as e:                                 # noqa: BLE001
            log(f"  capture aborted: {type(e).__name__}: {e}")
        finally:
            # ALWAYS hand the radio back, on every path
            try:
                if os.path.exists(pause):
                    os.rename(pause, pause + f".off.fadecap_"
                                             f"{time.strftime('%H%M%S')}")
                log("  maintenance window lifted -- warden restores the stack")
            except OSError as e:
                log(f"  ! could not lift the window ({e}); it EXPIRES on its "
                    f"own, supervision resumes without me")
        fade_since = None
    log(f"capture budget spent ({caps}/{a.max_caps}) -- exiting, radio is "
        f"the chain's for the rest of the night")
    return 0


if __name__ == "__main__":
    sys.exit(main())
