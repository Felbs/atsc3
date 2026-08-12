#!/usr/bin/env python3
"""Log RSS, disk and lane growth every 30 s during a soak.

Duration is the axis nothing has been tested along. A 21-minute run cannot
show a leak, a disk filling, or a clock drifting -- all three are silent
until they are fatal.
"""
import json, os, sys, time
import psutil

d = sys.argv[1] if len(sys.argv) > 1 else r"data\atsc3_live"
# Inside the run's own directory, not its parent. Writing to the parent put
# every run's samples into one shared data/monitor.jsonl, so E13's appended
# alongside E5's and the two could only be told apart by timestamp.
out = os.path.join(d, "monitor.jsonl")
t0 = time.time()
while True:
    procs = {}
    for p in psutil.process_iter(["pid", "cmdline", "memory_info", "cpu_percent"]):
        cl = " ".join(p.info["cmdline"] or [])
        for tag in ("atsc3 watch", "atsc3_audio", "atsc3_subs", "atsc3_run"):
            if tag in cl and "process_iter" not in cl:
                try:
                    procs[tag] = dict(rss_mb=p.info["memory_info"].rss / 1e6,
                                      cpu=p.cpu_percent(), pid=p.info["pid"])
                except Exception:
                    pass
    lanes = {}
    try:
        for k, v in json.load(open(os.path.join(d, "live.json")))["lanes"].items():
            lanes[k] = v["bytes"]
    except Exception:
        pass
    du = 0
    for root, _, fs in os.walk(d):
        for f in fs:
            try:
                du += os.path.getsize(os.path.join(root, f))
            except OSError:
                pass
    rec = dict(t=time.time() - t0, procs=procs, lanes=lanes, dir_mb=du / 1e6,
               free_gb=psutil.disk_usage(d).free / 1e9,
               ram_pct=psutil.virtual_memory().percent)
    with open(out, "a") as fh:
        fh.write(json.dumps(rec) + "\n")
    time.sleep(30)
