#!/usr/bin/env python3
"""atsc3_vs -- head-to-head scoreboard: our software tuner vs the FLEX 4K.

The operator's framing: *"make sure our software tuner goes tic for tac with
the HDHomeRun silicon tuner."*  This is the instrument that answers it with
numbers, per channel and per capability.

WHY THE COMPARISON IS FAIR NOW
------------------------------
Since 8/10 both receivers share ONE antenna through a passive splitter (E66),
so every paired reading is same-air and same-instant.  There is no cable to
swap and no propagation excuse available to either side.

THE SIX DIMENSIONS
    D1 coverage      every carrier/service each side can actually present
    D2 concurrent    same-instant delivery, from the paired diurnal log
    D3 sensitivity   *** NOT COMPARABLE -- and this file says so ***
    D4 acquisition   cold tune-to-first-payload, both sides measured
    D5 stability     who holds through the night, from the same paired log
    D6 capability    what each can DO, each row carrying its evidence

THE ANTI-FLATTERY RULES, BUILT IN
---------------------------------
This campaign's culture is that a favourable number you cannot defend is worse
than an honest loss.  Four rules are therefore enforced in code, not in
intentions:

1.  **A missing sample of OUR side counts as our downtime, never as a skip.**
    The diurnal logger samples on a clock; if our chain was restarting, its
    FEC is null.  Dropping those rows would silently delete exactly the
    moments we lost.  `_ours_up()` maps null to DOWN.

2.  **A verdict that depends on a threshold must be shown across thresholds.**
    D2 and D5 rest on "delivering" cut-offs for two different meters (our FEC
    %, their snq %).  The scoreboard sweeps them and reports whether the
    verdict is stable; an unstable verdict is downgraded to NOT-COMPARABLE.

3.  **A dimension the method cannot support is reported as NOT COMPARABLE,
    with the reason** -- never quietly omitted and never softened into a tie.
    D3 is the case: we can attenuate ourselves in software via `--rfgain`,
    the FLEX has no such control, so no equal-footing threshold sweep exists.

4.  **The reserved tuner is never retuned.**  tuner1 is the live soak's RF
    referee.  Every live probe here picks an ATSC-3.0-capable tuner that is
    NOT the reserved one and NOT already in use, and restores it to `none`
    on every exit path including exceptions.

MEASURED HARDWARE FACT worth stating up front: on this FLEX 4K only **tuner0
and tuner1** accept the `atsc3` modulation.  tuner2 and tuner3 answer
"ERROR: invalid channel" to `atsc3:...` and HTTP 503 "806 Tune Failed" to an
ATSC 3.0 vchannel, while locking 8VSB on the same feed in under a second.
So the device has 4 tuners but only 2 ATSC 3.0 front ends -- which also means
the reserved-tuner rule leaves exactly one usable for probing.

    python tools/atsc3_vs.py --probe      # take live FLEX readings, bank them
    python tools/atsc3_vs.py              # scoreboard from banked + our logs
    python tools/atsc3_vs.py --json
"""
from __future__ import annotations

import argparse
import contextlib
import glob
import json
import os
import re
import statistics as st
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
LAB = os.path.join(ROOT, "lab")
for p in (HERE, LAB):
    if p not in sys.path:
        sys.path.insert(0, p)

import atsc3_inspect as I   # noqa: E402

BUILD = "E81-2026-08-11"

OUT = os.path.join(LAB, "e81_out")
PROBE_JSON = os.path.join(OUT, "flex_probe.json")

# --- the FLEX, and the rule that protects the soak -------------------------
HDHR_BIN = r"C:\Program Files\Silicondust\HDHomeRun\hdhomerun_config.exe"
HDHR_ID = os.environ.get("ATSC3_HDHR_ID", "")  # no default; discover if unset
HDHR_HTTP = os.environ.get("ATSC3_HDHR_HTTP", "")     # resolved by discover()
RESERVED_TUNER = os.environ.get("ATSC3_HDHR_TUNER", "1")   # the soak's referee
N_TUNERS = 4

# "delivering" cut-offs.  HDHR_SNQ_OK matches tools/atsc3_run.py so the
# scoreboard and the supervisor cannot disagree about what a good link is.
OURS_FEC_OK = 95.0
THEIRS_SNQ_OK = 80
FEC_SWEEP = (80.0, 90.0, 95.0, 99.0)
SNQ_SWEEP = (60, 70, 80, 90)

LIVE_DIR = os.path.join(ROOT, "data", "e31")

WIN, LOSE, TIE, NC = "WIN", "LOSE", "TIE", "NOT-COMPARABLE"


def _rel(p):
    try:
        return os.path.relpath(p, ROOT).replace("\\", "/")
    except Exception:
        return str(p)


# ========================================================== FLEX plumbing

def _hd(*args, timeout=15, with_id=True):
    """Run hdhomerun_config.  `discover` is the one subcommand that must NOT
    be given a device id -- passing one makes it fail silently, which cost a
    probe run that reported 'FLEX not discovered' while the device was up."""
    if not os.path.exists(HDHR_BIN):
        return None
    cmd = [HDHR_BIN] + ([HDHR_ID] if with_id else []) + list(args)
    try:
        r = subprocess.run(cmd, capture_output=True, text=True,
                           timeout=timeout)
        return (r.stdout or "") + (r.stderr or "")
    except (OSError, subprocess.SubprocessError):
        return None


def discover():
    """The FLEX's IPv4 address.  It is DHCP-mobile (.48 and .27 both seen),
    so it is discovered every run rather than pinned in a constant."""
    global HDHR_HTTP
    if HDHR_HTTP:
        return HDHR_HTTP
    out = _hd("discover", with_id=False) or ""
    m = re.search(r"found at (\d+\.\d+\.\d+\.\d+)", out)
    HDHR_HTTP = m.group(1) if m else ""
    return HDHR_HTTP


def tuner_status(n):
    out = _hd("get", "/tuner%s/status" % n) or ""
    m = re.search(r"ch=(\S+)\s+lock=(\S+)\s+ss=(\d+)\s+snq=(\d+)\s+seq=(\d+)"
                  r"(?:\s+bps=(\d+))?", out)
    if not m:
        return None
    return {"ch": m.group(1), "lock": m.group(2), "ss": int(m.group(3)),
            "snq": int(m.group(4)), "seq": int(m.group(5)),
            "bps": int(m.group(6) or 0)}


def atsc3_capable(n):
    """Does tuner n accept the atsc3 modulation at all?

    Probe-and-restore: set a known ATSC 3.0 frequency and read the error.
    tuner2/tuner3 on this box answer 'ERROR: invalid channel'.
    """
    out = _hd("set", "/tuner%s/channel" % n, "atsc3:587000000")
    ok = out is not None and "invalid" not in (out or "").lower()
    _hd("set", "/tuner%s/channel" % n, "none")
    return ok


def pick_tuner():
    """An ATSC-3.0-capable tuner that is neither reserved nor busy.

    Returns (n, why) or (None, why-not).  Refusing is a valid outcome: a
    measurement we cannot take safely is not one we invent.
    """
    for n in range(N_TUNERS):
        if str(n) == str(RESERVED_TUNER):
            continue
        s = tuner_status(n)
        if s is None:
            continue
        if s["ch"] != "none":
            continue                       # somebody else is using it
        if not atsc3_capable(n):
            continue
        return n, "tuner%d: idle, ATSC 3.0 capable, not the reserved tuner" % n
    return None, ("no ATSC-3.0-capable tuner is free (reserved=tuner%s is "
                  "the soak's referee and is never retuned)" % RESERVED_TUNER)


@contextlib.contextmanager
def borrowed(n):
    """Borrow a tuner and ALWAYS hand it back, exceptions included."""
    try:
        yield n
    finally:
        _hd("set", "/tuner%s/channel" % n, "none")


def probe_carrier(n, freq_hz):
    """Lock, plpinfo and streaminfo for one carrier, from the FLEX itself."""
    _hd("set", "/tuner%s/channel" % n, "auto:%d" % freq_hz)
    best = None
    for _ in range(14):
        s = tuner_status(n)
        if s and s["lock"] != "none":
            best = s if (best is None or s["seq"] > best["seq"]) else best
            if s["seq"] >= 100:
                break
        time.sleep(0.6)
    plp_raw = _hd("get", "/tuner%s/plpinfo" % n) or ""
    si_raw = _hd("get", "/tuner%s/streaminfo" % n) or ""
    plps = []
    for m in re.finditer(r"^\s*(\d+):\s*(.*)$", plp_raw, re.M):
        body = m.group(2)
        if "=" not in body:
            continue
        d = dict(re.findall(r"(\w+)=(\S+)", body))
        d["plp"] = m.group(1)
        plps.append(d)
    bsid = None
    mb = re.search(r"bsid=(\S+)", plp_raw)
    if mb:
        bsid = int(mb.group(1), 16) if mb.group(1).startswith("0x") \
            else int(mb.group(1))
    svcs = []
    for m in re.finditer(r"^\s*(\d+):\s*([\d.]+)\s+(\S+)(\s+\(no data\))?\s*$",
                         si_raw, re.M):
        svcs.append({"program": m.group(1), "vchannel": m.group(2),
                     "name": m.group(3), "no_data": bool(m.group(4))})
    return {"freq_hz": freq_hz, "status": best, "bsid": bsid,
            "plps": plps, "services": svcs,
            "plpinfo_raw": plp_raw.strip(), "streaminfo_raw": si_raw.strip()}


def probe_acquisition(n, vchannel, runs=5):
    """Cold tune-to-first-payload-byte, the FLEX side.

    Directly comparable to our chain's start -> first MPU: both are "from the
    instant the receiver is told to tune, how long until it hands over usable
    payload".  The tuner is forced back to `none` before every run so each is
    genuinely cold.
    """
    ip = discover()
    if not ip:
        return {"error": "device not discovered"}
    url = "http://%s:5004/tuner%s/v%s" % (ip, n, vchannel)
    times, codes = [], []
    for _ in range(runs):
        _hd("set", "/tuner%s/channel" % n, "none")
        time.sleep(2.0)
        r = subprocess.run(
            ["curl", "-s", "-o", os.devnull, "-m", "25",
             "--max-filesize", "200000",
             "-w", "%{http_code} %{time_starttransfer}", url],
            capture_output=True, text=True, timeout=40)
        parts = (r.stdout or "").split()
        if len(parts) == 2:
            codes.append(parts[0])
            if parts[0] == "200":
                times.append(float(parts[1]))
    return {"url": url, "runs": runs, "codes": codes, "ttfb_s": times,
            "median_s": round(st.median(times), 3) if times else None,
            "min_s": round(min(times), 3) if times else None,
            "max_s": round(max(times), 3) if times else None,
            "metric": "time to first payload byte over HTTP, cold tune"}


def lineup(ip):
    try:
        import urllib.request
        with urllib.request.urlopen("http://%s/lineup.json" % ip,
                                    timeout=10) as r:
            return json.load(r)
    except Exception:
        return None


# Carriers we have SLTs for, so the two sides can be lined up by bsid.
CARRIER_FREQ = {33: 587000000, 30: 569000000, 25: 539000000}


def run_probe(runs=5):
    os.makedirs(OUT, exist_ok=True)
    ip = discover()
    rep = {"build": BUILD, "when": time.time(), "device_ip": ip,
           "reserved_tuner": RESERVED_TUNER}
    if not ip:
        rep["error"] = "FLEX not discovered"
        json.dump(rep, open(PROBE_JSON, "w"), indent=2)
        return rep
    rep["lineup"] = lineup(ip)
    rep["tuner_atsc3_capable"] = {}
    for n in range(N_TUNERS):
        if str(n) == str(RESERVED_TUNER):
            rep["tuner_atsc3_capable"][str(n)] = "reserved -- not probed"
            continue
        s = tuner_status(n)
        rep["tuner_atsc3_capable"][str(n)] = (
            atsc3_capable(n) if s and s["ch"] == "none" else "busy")
    n, why = pick_tuner()
    rep["tuner_used"] = n
    rep["tuner_choice"] = why
    rep["reserved_tuner_before"] = tuner_status(RESERVED_TUNER)
    if n is None:
        json.dump(rep, open(PROBE_JSON, "w"), indent=2)
        return rep
    with borrowed(n):
        rep["carriers"] = {}
        for rf, hz in sorted(CARRIER_FREQ.items()):
            rep["carriers"][str(rf)] = probe_carrier(n, hz)
        rep["acquisition"] = probe_acquisition(n, "107.1", runs=runs)
    rep["tuner_after"] = tuner_status(n)
    rep["reserved_tuner_after"] = tuner_status(RESERVED_TUNER)
    json.dump(rep, open(PROBE_JSON, "w"), indent=2)
    return rep


# ===================================================== our side telemetry

def gain_fix_line(lines):
    """Index of the first rfgain 4 READBACK -- the era boundary, in lines."""
    for i, l in enumerate(lines):
        if '"rfgain_sel": "4"' in l:
            return i
    return 0


def our_acquisition(chain_log=None, window=120, from_line=0):
    """Cold start -> bootstrap, and cold start -> first MPU, from chain.log.

    Only starts that reach lock inside `window` are timed; the ones that do
    not are counted and reported as a failure rate, because a tuner that
    sometimes never acquires cannot be scored on its good runs alone.
    """
    chain_log = chain_log or os.path.join(LIVE_DIR, "chain.log")
    if not os.path.isfile(chain_log):
        return None
    lines = open(chain_log, encoding="utf-8", errors="replace").read().splitlines()

    def ts(l):
        m = re.match(r"\[(\d\d):(\d\d):(\d\d)\]", l)
        return None if not m else (int(m.group(1)) * 3600 +
                                   int(m.group(2)) * 60 + int(m.group(3)))
    starts = [i for i, l in enumerate(lines)
              if "ATSC 3.0 live --" in l and i >= from_line]
    boot, mpu, nolock = [], [], 0
    for s in starts:
        t0 = ts(lines[s])
        if t0 is None:
            continue
        tb = tm = None
        for j in range(s + 1, len(lines)):
            t = ts(lines[j])
            if t is None:
                continue
            d = (t - t0) % 86400
            if d > window or "ATSC 3.0 live --" in lines[j]:
                break
            if tb is None and "bootstrap acquired" in lines[j]:
                tb = d
            if tm is None:
                m = re.search(r"MPU (\d+) seg", lines[j])
                if m and int(m.group(1)) >= 1:
                    tm = d
            if tb is not None and tm is not None:
                break
        if tb is None:
            nolock += 1
        else:
            boot.append(tb)
            if tm is not None:
                mpu.append(tm)

    def stats(v):
        if not v:
            return None
        sv = sorted(v)
        return {"n": len(v), "median_s": st.median(v), "min_s": min(v),
                "max_s": max(v), "p90_s": sv[min(len(sv) - 1,
                                                 9 * len(sv) // 10)]}
    return {"source": _rel(chain_log), "cold_starts": len(starts),
            "window_s": window,
            "no_lock_within_window": nolock,
            "no_lock_pct": round(100 * nolock / len(starts), 1) if starts else None,
            "start_to_bootstrap": stats(boot),
            "start_to_first_mpu": stats(mpu)}


def paired_samples(live_dir=None):
    """The same-air, same-instant record: our FEC beside their snq.

    Written every 5 minutes by tools/atsc3_diurnal.py while both receivers
    watch the same antenna through the splitter.  This is the only dataset in
    the campaign where "who was delivering at 03:14" is answerable for both
    sides without an assumption.
    """
    live_dir = live_dir or LIVE_DIR
    rows, srcs = [], []
    for f in sorted(glob.glob(os.path.join(live_dir, "_warden",
                                           "diurnal_*.jsonl"))):
        srcs.append(_rel(f))
        for line in open(f, encoding="utf-8", errors="replace"):
            line = line.strip()
            if not line:
                continue
            try:
                rows.append(json.loads(line))
            except ValueError:
                continue
    rows.sort(key=lambda r: r.get("t", 0))
    return rows, srcs


def _ours_up(r, fec_ok):
    """ANTI-FLATTERY RULE 1: a null FEC is our downtime, not a missing row.

    The logger samples on a wall clock. If our chain was restarting or wedged
    it has no FEC to report -- exactly the moments we lost. Treating null as
    "no data, skip" would delete our failures from our own scorecard."""
    v = r.get("fec_5min")
    return (v is not None) and (v >= fec_ok)


def _theirs_up(r, snq_ok):
    h = r.get("hdhr")
    if not h:
        return None                      # instrument silent: excluded, both ways
    return h.get("lock") not in (None, "none") and (h.get("snq") or 0) >= snq_ok


def contingency(rows, fec_ok, snq_ok):
    both = ours_only = theirs_only = neither = skipped = 0
    for r in rows:
        t = _theirs_up(r, snq_ok)
        if t is None:
            skipped += 1
            continue
        o = _ours_up(r, fec_ok)
        if o and t:
            both += 1
        elif o:
            ours_only += 1
        elif t:
            theirs_only += 1
        else:
            neither += 1
    n = both + ours_only + theirs_only + neither
    return {"n": n, "both_up": both, "ours_only": ours_only,
            "theirs_only": theirs_only, "neither": neither,
            "referee_silent_excluded": skipped,
            "ours_up_pct": round(100 * (both + ours_only) / n, 2) if n else None,
            "theirs_up_pct": round(100 * (both + theirs_only) / n, 2) if n else None}


# ============================================================== verdicts

def verdict(ours, theirs, higher_is_better=True, tol=0.0):
    if ours is None or theirs is None:
        return NC
    if abs(ours - theirs) <= tol:
        return TIE
    better = ours > theirs if higher_is_better else ours < theirs
    return WIN if better else LOSE


def dim(name, question, v, ours, theirs, method, caveats=None, sources=None,
        detail=None):
    return {"dimension": name, "question": question, "verdict": v,
            "ours": ours, "theirs": theirs, "method": method,
            "caveats": caveats or [], "sources": sources or [],
            "detail": detail or {}}


# ------------------------------------------------------------ D1 coverage

# What "we can present it" means, per service, in descending strength.
OURS_RANK = {"DECODED-LIVE": 4, "DECODED-OFFLINE": 3, "SIGNALLING-ONLY": 1,
             "LOCKED-DRM": 0, "BROADBAND-ONLY": 0, "DATA": 1}
THEIRS_RANK = {"PLAYS": 4, "NO-DATA": 1, "LOCKED-DRM": 0, "ABSENT": 0}


def ours_state(svc):
    p = svc["playable"]
    if svc["protection"]["state"] == "PROTECTED":
        return "LOCKED-DRM"
    # A service the live chain is decoding RIGHT NOW outranks every inference
    # about banked objects.  Without this, WJLA -- which our chain has been
    # decoding continuously for hours -- scored SIGNALLING-ONLY, because the
    # MMTP lane files are not named like ROUTE objects and _has_media missed
    # them.  Understating our own side is as much an error as overstating it.
    if svc.get("live"):
        return "DECODED-LIVE"
    if svc.get("broadbandAccessRequired") == "true":
        return "BROADBAND-ONLY"
    if svc.get("serviceCategory") == "4":
        return "DATA"
    if p["playable"]:
        # DECODED-LIVE requires POSITIVE evidence of a running lane, which is
        # the `live` branch above.  The inspector's `realtime` flag only says
        # "this PLP is not known to be too slow", and for a carrier absent
        # from the PLP table that defaults to true -- which labelled Fox
        # 45.100 (decoded offline from a banked capture, E76) as live.
        # Absence of a known blocker is not evidence of real-time operation.
        return "DECODED-OFFLINE"
    return "SIGNALLING-ONLY"


def flex_vchannel(v):
    """Our virtual channel -> the FLEX's numbering for the same service.

    MEASURED, not assumed: the FLEX renumbers every ATSC 3.0 service by
    +100 on the major -- 32.1->132.1, 7.1->107.1, 45.1->145.1, 58.1->158.1
    on all three carriers we hold SLTs for.  It is not a coincidence and it
    is not our channel; matching without it silently pairs the wrong rows.
    """
    try:
        maj, minr = str(v).split(".")
        return "%d.%s" % (int(maj) + 100, minr)
    except Exception:
        return None


def theirs_state(svc, bsid, rf, probe):
    """What the FLEX can present, matched CARRIER-SCOPED and BY CHANNEL.

    Two real bugs lived in the earlier name-based version, and both moved the
    scoreboard in the FLEX's favour:

    1.  Appearing in lineup.json is NOT evidence of playability.  WHUT 32.1 is
        in its lineup with no DRM flag, while its own streaminfo marks it
        "(no data)" -- it cannot present ROUTE/DASH as a programme at all
        (E67).  Reading lineup first credited the FLEX with a service its own
        instrument says it cannot show.

    2.  Names collide ACROSS STANDARDS.  "WBFF" is both the ATSC 3.0 service
        on RF25 (145.1, "(no data)") and an ATSC 1.0 MPEG2 service in the
        lineup (45.1 WBFF45).  Matching by name scored the 3.0 service as
        PLAYS on the strength of an entirely different 1.0 broadcast.

    So: find the FLEX's reading of THIS carrier (by bsid, falling back to RF),
    then match the service by its renumbered virtual channel.  lineup.json is
    consulted only for the DRM flag, which is what it is actually good for.
    """
    car = None
    for k, c in (probe.get("carriers") or {}).items():
        if bsid is not None and c.get("bsid") is not None \
                and str(c["bsid"]) == str(bsid):
            car = c
            break
        if rf is not None and str(k) == str(rf):
            car = c
    if car is None:
        return "ABSENT"
    want = flex_vchannel(svc.get("vchannel"))
    hit = None
    for s in car.get("services") or []:
        if want and s.get("vchannel") == want:
            hit = s
            break
    if hit is None:
        # data/broadband services the FLEX never enumerates
        return "ABSENT"
    if hit.get("no_data"):
        return "NO-DATA"
    for e in (probe.get("lineup") or []):
        if e.get("GuideNumber") == want:
            return "LOCKED-DRM" if e.get("DRM") else "PLAYS"
    return "PLAYS"


def gain_fix_epoch(chain_log=None, live_dir=None):
    """When the current front-end configuration began (E73's rfgain 4).

    Taken from the chain's own READBACK line, not from the flag we asked for
    -- `unverified_control_writes_lie`.  Everything before this timestamp
    describes a receiver that no longer exists, so the scoreboard reports the
    era split rather than letting an obsolete configuration stand in for the
    current one (or, worse, quietly dropping it).
    """
    chain_log = chain_log or os.path.join(live_dir or LIVE_DIR, "chain.log")
    if not os.path.isfile(chain_log):
        return None
    day = time.localtime(os.path.getmtime(chain_log))
    first = None
    for line in open(chain_log, encoding="utf-8", errors="replace"):
        if '"rfgain_sel": "4"' not in line:
            continue
        m = re.match(r"\[(\d\d):(\d\d):(\d\d)\]", line)
        if not m:
            continue
        first = (int(m.group(1)), int(m.group(2)), int(m.group(3)))
        break
    if not first:
        return None
    # the log is one day's worth; anchor on the file's own date
    t = time.mktime((day.tm_year, day.tm_mon, day.tm_mday,
                     first[0], first[1], first[2], 0, 0, -1))
    if t > time.time():
        t -= 86400
    return t


def d1_coverage(model, probe):
    rows = []
    for bsid, c in model["carriers"].items():
        rf = (c.get("rf") or {}).get("rf")
        for svc in c["services"]:
            o = ours_state(svc)
            t = theirs_state(svc, bsid, rf, probe)
            # PRESENTABILITY, not rank arithmetic.  Ranking LOCKED-DRM (0)
            # below NO-DATA (1) scored every Widevine service as a LOSS for
            # us -- as if the FLEX printing "(no data)" beside a channel
            # neither receiver may decrypt were an advantage. Neither side can
            # present it; that is a TIE, and our guide+captions edge on those
            # services is scored in D6 rather than smuggled in here.
            po = OURS_RANK.get(o, 0) >= 3
            pt = THEIRS_RANK.get(t, 0) >= 3
            v = TIE if po == pt else (WIN if po else LOSE)
            rows.append({"rf": rf, "bsid": bsid, "name": svc.get("name"),
                         "vchannel": svc.get("vchannel"),
                         "ours": o, "theirs": t, "verdict": v})
    rows.sort(key=lambda r: (-(r["rf"] or 0), r["vchannel"] or "zz"))
    w = sum(1 for r in rows if r["verdict"] == WIN)
    l = sum(1 for r in rows if r["verdict"] == LOSE)
    tie = sum(1 for r in rows if r["verdict"] == TIE)
    return dim(
        "coverage",
        "which services can each side actually present?",
        WIN if w > l else (LOSE if l > w else TIE),
        "%d services presentable" % sum(1 for r in rows
                                        if OURS_RANK.get(r["ours"], 0) >= 3),
        "%d services presentable" % sum(1 for r in rows
                                        if THEIRS_RANK.get(r["theirs"], 0) >= 3),
        "our state from the four-layer inspector (E70); theirs from its own "
        "lineup.json DRM flag and per-carrier streaminfo '(no data)' marker",
        caveats=[
            "'(no data)' is a TRANSPORT flag, not a DRM flag (E67): the FLEX "
            "cannot present ROUTE/DASH services as programmes at all, which "
            "is why we win the ROUTE rows rather than out-receiving it",
            "our DECODED-OFFLINE rows are real decodes that are NOT real "
            "time (PLP1 is ~140x short, E47) -- counted, but never as live",
        ],
        sources=["lab/*/lls_*SLT*.xml", "FLEX lineup.json + streaminfo"],
        detail={"rows": rows, "wins": w, "losses": l, "ties": tie})


# ------------------------------------------- D2 concurrent, D5 stability

def era_note(gain_fix):
    if not gain_fix:
        return "no era split available"
    return ("ERA SPLIT at %s, the first rfgain 4 READBACK (E73): everything "
            "before it describes a front end we no longer run"
            % time.strftime("%m-%d %H:%M", time.localtime(gain_fix)))


def d2_concurrent(rows, srcs, gain_fix=None):
    if not rows:
        return dim("concurrent delivery", "who was delivering, same instant?",
                   NC, None, None, "no paired samples found",
                   caveats=["tools/atsc3_diurnal.py has not run"], sources=srcs)
    base = contingency(rows, OURS_FEC_OK, THEIRS_SNQ_OK)
    # ANTI-FLATTERY RULE 2: is the verdict an artifact of the cut-offs?
    sweep = []
    verdicts = set()
    for f in FEC_SWEEP:
        for s in SNQ_SWEEP:
            c = contingency(rows, f, s)
            v = verdict(c["ours_up_pct"], c["theirs_up_pct"], tol=0.5)
            sweep.append({"fec_ok": f, "snq_ok": s,
                          "ours_up_pct": c["ours_up_pct"],
                          "theirs_up_pct": c["theirs_up_pct"], "verdict": v})
            verdicts.add(v)
    stable = len(verdicts) == 1
    v = verdict(base["ours_up_pct"], base["theirs_up_pct"], tol=0.5)
    if not stable:
        v = NC
    # The current configuration, reported beside the whole window.  The
    # HEADLINE verdict stays on the full window -- the conservative, less
    # flattering choice -- but a reader who wants to know what the machine
    # does TODAY must not have to take that on trust.
    cur = None
    if gain_fix:
        post = [r for r in rows if r.get("t", 0) >= gain_fix]
        if post:
            c = contingency(post, OURS_FEC_OK, THEIRS_SNQ_OK)
            cur = {"since": gain_fix, "n": c["n"],
                   "ours_up_pct": c["ours_up_pct"],
                   "theirs_up_pct": c["theirs_up_pct"],
                   "verdict": verdict(c["ours_up_pct"], c["theirs_up_pct"],
                                      tol=0.5)}
    caveats = [
        "FEC% and snq% are DIFFERENT METERS on different scales; this "
        "compares 'was each side above its own good-link threshold', not "
        "the numbers themselves",
        "verdict %s across the %dx%d threshold sweep" %
        ("STABLE" if stable else "UNSTABLE -> downgraded to NOT-COMPARABLE",
         len(FEC_SWEEP), len(SNQ_SWEEP)),
        "a null FEC sample is counted as OUR downtime, never skipped",
        era_note(gain_fix),
    ]
    if cur:
        caveats.append(
            "current configuration only (n=%d): ours %.2f%% vs theirs "
            "%.2f%% -> %s. The headline verdict above deliberately uses the "
            "FULL window, which is the less flattering number."
            % (cur["n"], cur["ours_up_pct"], cur["theirs_up_pct"],
               cur["verdict"]))
    return dim(
        "concurrent delivery",
        "over the same minutes on the same antenna, who was up?",
        v,
        "%.2f%% of samples delivering" % (base["ours_up_pct"] or 0),
        "%.2f%% of samples delivering" % (base["theirs_up_pct"] or 0),
        "paired 5-minute samples; ours = chain FEC >= %.0f%%, theirs = "
        "lock + snq >= %d" % (OURS_FEC_OK, THEIRS_SNQ_OK),
        caveats=caveats,
        sources=srcs, detail={"base": base, "sweep": sweep,
                              "threshold_stable": stable,
                              "current_config": cur})


def _runs_down(rows, fec_ok, snq_ok, side):
    """Longest consecutive stretch (minutes) with a side not delivering."""
    worst = cur = 0
    prev_t = None
    for r in rows:
        up = _ours_up(r, fec_ok) if side == "ours" else _theirs_up(r, snq_ok)
        if up is None:
            continue
        t = r.get("t", 0)
        gap = ((t - prev_t) / 60.0) if prev_t else 5.0
        prev_t = t
        if up:
            cur = 0
        else:
            cur += gap
            worst = max(worst, cur)
    return round(worst, 1)


def d5_stability(rows, srcs, gain_fix=None):
    if not rows:
        return dim("stability", "who holds lock through the fades?", NC,
                   None, None, "no paired samples", sources=srcs)
    span_h = (rows[-1]["t"] - rows[0]["t"]) / 3600.0 if len(rows) > 1 else 0
    ours_worst = _runs_down(rows, OURS_FEC_OK, THEIRS_SNQ_OK, "ours")
    theirs_worst = _runs_down(rows, OURS_FEC_OK, THEIRS_SNQ_OK, "theirs")
    v = verdict(ours_worst, theirs_worst, higher_is_better=False, tol=5.0)
    cur = None
    if gain_fix:
        post = [r for r in rows if r.get("t", 0) >= gain_fix]
        if len(post) > 1:
            cur = {"since": gain_fix, "samples": len(post),
                   "span_hours": round((post[-1]["t"] - post[0]["t"]) / 3600, 2),
                   "ours_longest_outage_min":
                       _runs_down(post, OURS_FEC_OK, THEIRS_SNQ_OK, "ours"),
                   "theirs_longest_outage_min":
                       _runs_down(post, OURS_FEC_OK, THEIRS_SNQ_OK, "theirs")}
            cur["verdict"] = verdict(cur["ours_longest_outage_min"],
                                     cur["theirs_longest_outage_min"],
                                     higher_is_better=False, tol=5.0)
    caveats = [
        "the window is %.1f h and includes the pre-gain-fix era; our earlier "
        "hours are front-end overload (E73/E75), not the sky" % span_h,
        "one night is not a stability claim; E75 says so explicitly and the "
        "PRE-DAWN floor has not been tested yet",
        era_note(gain_fix),
    ]
    if cur:
        caveats.append(
            "current configuration only (%.1f h, n=%d): our longest outage "
            "%.0f min vs theirs %.0f min -> %s"
            % (cur["span_hours"], cur["samples"],
               cur["ours_longest_outage_min"],
               cur["theirs_longest_outage_min"], cur["verdict"]))
    return dim(
        "stability",
        "who holds through the fades, over the logged window?",
        v,
        "longest outage %.0f min" % ours_worst,
        "longest outage %.0f min" % theirs_worst,
        "longest consecutive not-delivering run over %.1f h of paired "
        "5-minute samples" % span_h,
        caveats=caveats,
        sources=srcs,
        detail={"span_hours": round(span_h, 2), "samples": len(rows),
                "ours_longest_outage_min": ours_worst,
                "theirs_longest_outage_min": theirs_worst,
                "current_config": cur})


# ------------------------------------------------------- D3 sensitivity

def d3_sensitivity():
    """The dimension we REFUSE to score, and exactly why.

    A sensitivity comparison needs both receivers walked down the same
    attenuation ramp until each fails.  We can attenuate ourselves in software
    (`--rfgain`, and E73 did precisely that live).  The FLEX exposes no gain
    or attenuation control at all -- `/sys/features` lists channelmaps and
    modulations and nothing else.  Without a common ramp there is no
    equal-footing threshold, and the honest report is that the method cannot
    support the claim.

    What is NOT a substitute, and why:
      * comparing our SNR to its snq -- different meters, undocumented scale;
      * inserting a physical attenuator -- that touches the shared antenna
        feed, i.e. the live soak's signal path, which is frozen;
      * using the fade record as a threshold -- fades move both receivers at
        once, so it measures WHO FAILED FIRST on the ramp the sky provided,
        not the level at which each fails. That is a real and useful number,
        and it is reported under D2/D5 as delivery, not relabelled here.
    """
    return dim(
        "sensitivity / threshold",
        "at what signal level does each receiver fail?",
        NC, "not measurable on equal footing", "no gain/attenuation control",
        "refused: no common attenuation ramp exists",
        caveats=[
            "we CAN attenuate ourselves in software (--rfgain); the FLEX "
            "exposes no such control, so a matched sweep is impossible",
            "a physical attenuator would touch the shared feed the frozen "
            "soak is using, so it is not available either",
            "the closest defensible proxy -- who failed first during real "
            "fades -- is reported as D2/D5 delivery and deliberately NOT "
            "relabelled as a threshold",
        ],
        sources=["/sys/features (no gain control listed)"])


# ------------------------------------------------------- D4 acquisition

def d4_acquisition(ours, probe, ours_cur=None, gain_fix=None):
    acq = (probe or {}).get("acquisition") or {}
    theirs = acq.get("median_s")
    o = (ours or {}).get("start_to_first_mpu") or {}
    ours_med = o.get("median_s")
    v = verdict(ours_med, theirs, higher_is_better=False, tol=0.3)
    boot = (ours or {}).get("start_to_bootstrap") or {}
    return dim(
        "acquisition",
        "cold tune to first usable payload -- how long?",
        v,
        "%.1f s median (n=%s)" % (ours_med, o.get("n")) if ours_med else None,
        "%.2f s median (n=%s)" % (theirs, acq.get("runs")) if theirs else None,
        "ours: chain start -> first MPU segment, from chain.log; theirs: "
        "HTTP tune -> first payload byte, cold (tuner forced to none between "
        "runs)",
        caveats=[
            "our figure INCLUDES process start and decode-pool warm (~5.2 s "
            "of the ~%.0f s); bootstrap alone is %.1f s median" %
            (ours_med or 0, boot.get("median_s") or 0),
            "%s%% of our cold starts never locked within %s s -- those are "
            "excluded from the median and reported separately; a receiver "
            "that sometimes never acquires cannot be scored on its good runs "
            "alone" % ((ours or {}).get("no_lock_pct"),
                       (ours or {}).get("window_s")),
            "their n is small (%s cold runs in one sitting)" % acq.get("runs"),
            era_note(gain_fix),
        ] + ([
            "current configuration only: %s%% of cold starts fail to lock "
            "within %ss (n=%s starts), median start->first MPU %s s -- the "
            "no-lock rate above is dominated by the front-end overload era"
            % ((ours_cur or {}).get("no_lock_pct"),
               (ours_cur or {}).get("window_s"),
               (ours_cur or {}).get("cold_starts"),
               ((ours_cur or {}).get("start_to_first_mpu") or {}).get("median_s"))
        ] if ours_cur else []),
        sources=[(ours or {}).get("source"), _rel(PROBE_JSON)],
        detail={"ours": ours, "ours_current_config": ours_cur, "theirs": acq})


# -------------------------------------------------------- D6 capability

def d6_capability(model, probe):
    """What each side can DO.  Each row names the evidence for our column."""
    rows = [
        ("HEVC video decode", "yes", "yes", TIE,
         "both decode 1080p hvc1; theirs in silicon, ours in software"),
        ("AC-4 audio decode", "yes -- our own decoder, 5.1 discrete verified",
         "pass-through only", WIN,
         "E-AC4 campaign: the FLEX hands AC-4 to a TV to decode; we decode it "
         "ourselves (lab AC-4 decoder, 5.1 discrete-verified)"),
        ("second audio language", "yes (eng+spa lanes run concurrently)",
         "stream carries it; device does not select", TIE,
         "both deliver the spa track; selection is a player concern"),
        ("captions on CLEAR services", "yes (stpp -> live.srt)", "yes", TIE,
         "data/e31/live.srt"),
        ("captions on ENCRYPTED services", "yes", "no", WIN,
         "E70: stpp is unencrypted on Widevine services; the FLEX presents "
         "nothing at all for a ROUTE service"),
        ("guide on ENCRYPTED services", "yes (now/next from the ESG)", "no",
         WIN, "E70 grid, E40 SGDU decoder"),
        ("service guide (ESG) decode", "yes -- service guide (ESG) decoded",
         "no ESG exposed over the API", WIN, "lab/m41_sgdu.py"),
        ("ROUTE/DASH services", "yes (WHUT 32.1, WBFFMob 45.100)",
         "no -- '(no data)'", WIN, "E67 / E76"),
        ("MMTP services", "yes (WJLA 107.1 live)", "yes", TIE, "data/e31"),
        ("per-PLP detail", "yes, but attributed from decode artifacts",
         "yes -- plpinfo gives mod/cod/layer directly", LOSE,
         "the FLEX reports 'plp1 mod=qam256 cod=7/15 layer=enhanced' from the "
         "L1 itself; we infer the PLP from which baseband file the flow came "
         "out of"),
        ("real-time PLP1 / enhanced layer", "no (~140x short, E47)",
         "yes", LOSE, "E47 measured the working PLP1 decoder at 34.6 s per "
         "0.247 s frame"),
        ("DRM playback", "no -- and never", "no", TIE,
         "both locked out; neither attempts circumvention"),
        ("front-end gain control", "yes (--rfgain, live A/B in E73)",
         "no", WIN, "the control that fixed our night floor does not exist "
         "on the FLEX -- which is also why D3 is unscoreable"),
        ("cold acquisition speed", "~13 s", "~1.8 s", LOSE, "D4"),
        ("tuners available", "1 SDR", "4, of which 2 do ATSC 3.0", LOSE,
         "measured: tuner2/3 reject the atsc3 modulation"),
    ]
    out = [{"capability": a, "ours": b, "theirs": c, "verdict": d,
            "evidence": e} for a, b, c, d, e in rows]
    w = sum(1 for r in out if r["verdict"] == WIN)
    l = sum(1 for r in out if r["verdict"] == LOSE)
    return dim("capability parity", "what can each side DO?",
               WIN if w > l else (LOSE if l > w else TIE),
               "%d capability wins" % w, "%d capability wins" % l,
               "declared matrix; every row in our column names its evidence",
               caveats=["a declared matrix is an argument, not a measurement "
                        "-- each row is only as good as the evidence named"],
               detail={"rows": out, "wins": w, "losses": l})


# =================================================== independent checks

# Claims this campaign made from ITS OWN decodes, which the FLEX's plpinfo can
# check with a different receiver's L1 parser.  Not scored -- a corroboration
# is not a contest -- but it is the cheapest audit available to us.
CORROBORATE = [
    (33, "540", "bsid", "540",
     "our RF33 SLT (lab/m7_out/lls_1_SLT.xml)"),
    (33, "1", "plp1 subframe", "sfi=1",
     "E67: the four ROUTE services ride PLP1, SUBFRAME 1"),
    (33, "1", "plp1 modulation", "qam256",
     "E67: PLP1 is 256QAM-NUC where PLP0 is 16QAM"),
    (30, "1", "plp1 layer", "enhanced",
     "E69: both RF30 TV services ride the LDM ENHANCED layer"),
    (25, "1", "plp1 layer", "enhanced",
     "E76/E69: Fox RF25's enhanced layer"),
    (25, "1408", "bsid", "1408",
     "our RF25 SLT (lab/e75_fox_out/lls_1_SLT.xml), bsid 0x580"),
]


def corroboration(model, probe):
    out = []
    for rf, key, what, expect, ours_src in CORROBORATE:
        car = (probe.get("carriers") or {}).get(str(rf)) or {}
        got = None
        if what == "bsid":
            got = str(car.get("bsid")) if car.get("bsid") is not None else None
        else:
            for p in car.get("plps") or []:
                if p.get("plp") == key:
                    got = p.get("sfi") and ("sfi=" + p["sfi"]) \
                        if what.endswith("subframe") else \
                        (p.get("mod") if "modulation" in what else p.get("layer"))
        out.append({"rf": rf, "check": what, "ours_claim": expect,
                    "flex_reports": got,
                    "agrees": (got == expect) if got is not None else None,
                    "our_source": ours_src})
    return out


# ============================================================ scoreboard

def scoreboard(sources=None, probe=None, live_dir=None):
    probe = probe if probe is not None else load_probe()
    model = I.build_lineup(sources or (I.DEFAULT_SOURCES + ["lab/e75_fox_out"]))
    rows, srcs = paired_samples(live_dir)
    chain = os.path.join(live_dir or LIVE_DIR, "chain.log")
    gain_fix = gain_fix_epoch(chain)
    ours_acq = our_acquisition(chain)
    ours_acq_cur = None
    if os.path.isfile(chain):
        lines = open(chain, encoding="utf-8", errors="replace").read().splitlines()
        gl = gain_fix_line(lines)
        if gl:
            ours_acq_cur = our_acquisition(chain, from_line=gl)
    dims = [d1_coverage(model, probe),
            d2_concurrent(rows, srcs, gain_fix),
            d3_sensitivity(),
            d4_acquisition(ours_acq, probe, ours_acq_cur, gain_fix),
            d5_stability(rows, srcs, gain_fix),
            d6_capability(model, probe)]
    tally = {WIN: 0, LOSE: 0, TIE: 0, NC: 0}
    for d in dims:
        tally[d["verdict"]] = tally.get(d["verdict"], 0) + 1
    return {"build": BUILD, "when": time.time(),
            "probe_when": (probe or {}).get("when"),
            "gain_fix": gain_fix,
            "dimensions": dims, "tally": tally,
            "corroboration": corroboration(model, probe),
            "headline": headline(tally)}


def headline(t):
    if t[LOSE] > t[WIN]:
        return "BEHIND on more dimensions than ahead"
    if t[WIN] > t[LOSE]:
        return "AHEAD on more dimensions than behind"
    return "LEVEL"


def load_probe():
    try:
        return json.load(open(PROBE_JSON))
    except Exception:
        return {}


# ---------------------------------------------------------------- report

MARK = {WIN: "WIN ", LOSE: "LOSE", TIE: "TIE ", NC: "n/c "}


def report(sb, out=sys.stdout, verbose=False):
    w = out.write
    w("OUR SOFTWARE TUNER  vs  HDHomeRun FLEX 4K      build %s\n" % sb["build"])
    if sb.get("probe_when"):
        w("FLEX probed %s   (one antenna, passive splitter -- same air, "
          "same instant)\n" % time.strftime("%m-%d %H:%M",
                                            time.localtime(sb["probe_when"])))
    w("=" * 78 + "\n")
    for d in sb["dimensions"]:
        w("\n[%s] %-22s %s\n" % (MARK[d["verdict"]], d["dimension"],
                                 d["question"]))
        w("       ours   : %s\n" % d["ours"])
        w("       theirs : %s\n" % d["theirs"])
        w("       method : %s\n" % d["method"])
        for c in d["caveats"]:
            w("       caveat : %s\n" % c)
    w("\n" + "=" * 78 + "\n")
    t = sb["tally"]
    w("TALLY  win %d   lose %d   tie %d   not-comparable %d      %s\n"
      % (t[WIN], t[LOSE], t[TIE], t[NC], sb["headline"]))
    cov = next(d for d in sb["dimensions"] if d["dimension"] == "coverage")
    w("\nPER-SERVICE COVERAGE\n")
    w("  %-4s %-8s %-7s %-16s %-12s %s\n"
      % ("RF", "service", "vchan", "ours", "theirs", ""))
    for r in cov["detail"]["rows"]:
        w("  %-4s %-8s %-7s %-16s %-12s %s\n"
          % (r["rf"] or "-", r["name"], r["vchannel"] or "-",
             r["ours"], r["theirs"], MARK[r["verdict"]]))
    cor = sb.get("corroboration") or []
    if cor:
        w("\nINDEPENDENT CORROBORATION"
          " (the FLEX's own L1 parser checking our claims)\n")
        for c in cor:
            mark = "agrees" if c["agrees"] else (
                "n/a" if c["agrees"] is None else "DISAGREES")
            w("  RF%-3s %-16s ours %-10s flex %-10s %s\n"
              % (c["rf"], c["check"], c["ours_claim"], c["flex_reports"], mark))
    cap = next(d for d in sb["dimensions"]
               if d["dimension"] == "capability parity")
    w("\nCAPABILITY\n")
    for r in cap["detail"]["rows"]:
        w("  %s %-34s ours: %-42s theirs: %s\n"
          % (MARK[r["verdict"]], r["capability"], r["ours"][:42], r["theirs"]))
    return sb


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--probe", action="store_true",
                    help="take live FLEX readings and bank them")
    ap.add_argument("--runs", type=int, default=5)
    ap.add_argument("--source", action="append", default=None)
    ap.add_argument("--live-dir", default=None)
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--verbose", action="store_true")
    a = ap.parse_args(argv)
    I.be_polite()
    if a.probe:
        rep = run_probe(runs=a.runs)
        print("probe -> %s  (tuner used: %s; reserved tuner%s untouched: %s)"
              % (_rel(PROBE_JSON), rep.get("tuner_used"), RESERVED_TUNER,
                 (rep.get("reserved_tuner_after") or {}).get("ch")))
        return 0
    sb = scoreboard(a.source, live_dir=a.live_dir)
    if a.json:
        json.dump(sb, sys.stdout, indent=2, default=str)
        sys.stdout.write("\n")
    else:
        report(sb, verbose=a.verbose)
    return 0


if __name__ == "__main__":
    sys.exit(main())
