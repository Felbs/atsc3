#!/usr/bin/env python3
"""carrier_gain.py -- the right front-end gain is a property OF THE CARRIER.

WHY THIS EXISTS (measured 2026-08-10, two facts one setting cannot serve):
  * RF33 REQUIRES rfgain_sel 4. At 2 the chain read 0.0% FEC while an
    HDHomeRun on the SAME feed decoded at snq=100 -- total, not marginal
    (E73).
  * RF25/Fox PREFERS rfgain_sel 2 by +6.04 dB (18.89 vs 12.85 dB median
    SNR, 150 frames each, captures three minutes apart differing in
    nothing else). It does not overload at 2 because it is a far weaker
    carrier: -12.6 dBFS peak, zero full-scale samples (E75/E76).
The tuner is single-tenant, so whichever carrier we are on dictates the
setting. A global constant is therefore not merely inconvenient, it is
WRONG for one of these two channels whichever value it takes.

THREE RULES THIS MODULE ENFORCES
  1. LOOK IT UP, do not remember it. `lookup(rf)` always answers, and says
     whether the answer was measured or is the documented default.
  2. WRITE THEN READ BACK. `apply_gain` re-reads rfgain_sel and IFGR from
     the device and raises if either disagrees -- an unverified
     writeSetting is a hope, not a setting (fleet law, bought the hard way
     more than once).
  3. PROVENANCE OR IT DID NOT HAPPEN. `calibrate` returns an entry that
     carries when/how/what-number, ready to be written into the table;
     it never returns a bare integer.

This module is import-safe with no SDR present: SoapySDR is touched only
inside apply_gain/calibrate_live, so the table, the chooser and the gate
all run offline.

    python lab/carrier_gain.py                 # show the table
    python lab/carrier_gain.py --rf 33         # what would we use, and why
"""
from __future__ import annotations

import json
import os
import statistics
import time

HERE = os.path.dirname(os.path.abspath(__file__))
TABLE_PATH = os.path.join(HERE, "carrier_gain.json")


try:                                   # the real value when Soapy is present
    from SoapySDR import SOAPY_SDR_RX as _RX
except Exception:                                              # noqa: BLE001
    _RX = 1        # SoapySDR's RX direction. Duplicated ONLY so the readback
    #                path stays gateable on a box without SoapySDR installed;
    #                the import above wins whenever the real library exists,
    #                so this constant can never silently disagree with it in
    #                production (the shared-constant law: a copy nobody can
    #                cross-check is the dangerous kind -- this one is checked
    #                by the import itself).


class GainReadbackError(RuntimeError):
    """The device did not confirm what we asked for."""


# ------------------------------------------------------------------ table

def load(path=None):
    with open(path or TABLE_PATH, encoding="utf-8") as f:
        return json.load(f)


def lookup(rf, table=None):
    """-> dict(rfgain_sel, ifgr, measured, source, entry).

    ALWAYS answers. `measured` False means the value is the documented
    default and the carrier has never been measured -- callers should say
    so out loud rather than let an unmeasured carrier look calibrated.
    An entry present but carrying rfgain_sel null (we decoded it once but
    never recorded the gain) also falls back, and is reported distinctly
    so it can be told from a carrier we have never seen at all."""
    t = table or load()
    ent = (t.get("carriers") or {}).get(str(int(rf)))
    dflt = t["default"]
    if ent and ent.get("rfgain_sel") is not None:
        return dict(rfgain_sel=int(ent["rfgain_sel"]),
                    ifgr=int(ent.get("ifgr", dflt["ifgr"])),
                    measured=bool(ent.get("measured")),
                    source="table", entry=ent)
    return dict(rfgain_sel=int(dflt["rfgain_sel"]),
                ifgr=int(dflt.get("ifgr", 32)),
                measured=False,
                source="default-known-carrier" if ent else "default",
                entry=ent or dflt)


def describe(rf, table=None):
    """One human line for the log, including the provenance."""
    g = lookup(rf, table)
    ent = g["entry"] or {}
    if g["measured"]:
        m = ent.get("metric", {}) or {}
        how = f"{ent.get('evidence', '?')} {ent.get('method', '')}".strip()
        num = ""
        if m.get("kind") == "snr_db_median" and m.get("delta_db") is not None:
            num = (f", chosen over {_other(ent)} by "
                   f"{m['delta_db']:+.2f} dB")
        elif m.get("kind") == "fec_pct_mean":
            num = (f", FEC {m.get('at_2')}% at 2 vs "
                   f"{m.get('at_4')}% at 4")
        return (f"RF{rf}: rfgain_sel={g['rfgain_sel']} IFGR={g['ifgr']} "
                f"(measured {ent.get('measured_utc', '?')}, {how}{num})")
    why = "never measured"
    if g["source"] == "default-known-carrier":
        why = "decoded once but the gain was not recorded"
    return (f"RF{rf}: rfgain_sel={g['rfgain_sel']} IFGR={g['ifgr']} "
            f"** DEFAULT, {why} -- run calibration **")


def _other(ent):
    m = ent.get("metric", {}) or {}
    for k in m:
        if k.startswith("at_") and not k.endswith(str(ent.get("rfgain_sel"))):
            return f"rfgain_sel {k[3:]}"
    return "the alternative"


# ------------------------------------------------------- apply (with SDR)

def apply_gain(sdr, rf, table=None, log=print, settle=0.3,
               rfgain=None, ifgr=None):
    """Set IFGR + rfgain_sel for this carrier and CONFIRM BY READBACK.

    `rfgain`/`ifgr` override the table when given (an operator asking for
    a specific value on the command line must still get it -- the table
    decides only what we do when nobody said). The override is announced
    as an override, so a hand-set gain is never mistaken in the log for a
    calibrated one.

    Raises GainReadbackError if the device does not report what we asked
    for. The caller may catch it and stand down; what it must not do is
    proceed believing a write that never landed."""
    g = lookup(rf, table)
    if rfgain is not None or ifgr is not None:
        g = dict(g, rfgain_sel=int(rfgain if rfgain is not None
                                   else g["rfgain_sel"]),
                 ifgr=int(ifgr if ifgr is not None else g["ifgr"]),
                 measured=False, source="explicit")
        log(f"  gain OVERRIDE (command line): RF{rf} "
            f"rfgain_sel={g['rfgain_sel']} IFGR={g['ifgr']} "
            f"-- table says: {describe(rf, table)}")
    else:
        log(f"  gain table: {describe(rf, table)}")
    try:
        sdr.setGainMode(_RX, 0, False)
    except Exception:                                          # noqa: BLE001
        pass
    sdr.setGain(_RX, 0, "IFGR", g["ifgr"])
    sdr.writeSetting("rfgain_sel", str(g["rfgain_sel"]))
    time.sleep(settle)

    rg_raw = sdr.readSetting("rfgain_sel")
    try:
        rg = int(str(rg_raw).strip())
    except (TypeError, ValueError):
        raise GainReadbackError(
            f"RF{rf}: rfgain_sel read back as {rg_raw!r}, not a number")
    ifgr_raw = sdr.getGain(_RX, 0, "IFGR")
    ifgr = int(round(float(ifgr_raw)))
    bad = []
    if rg != g["rfgain_sel"]:
        bad.append(f"rfgain_sel asked {g['rfgain_sel']} got {rg}")
    if ifgr != g["ifgr"]:
        bad.append(f"IFGR asked {g['ifgr']} got {ifgr}")
    if bad:
        raise GainReadbackError(f"RF{rf}: " + "; ".join(bad))
    log(f"  gain confirmed by readback: rfgain_sel={rg} IFGR={ifgr}")
    return dict(rfgain_sel=rg, ifgr=ifgr, measured=g["measured"],
                source=g["source"])


# ------------------------------------------------------------- calibrate

def choose(observations):
    """Pick the best rfgain_sel from candidate observations.

    observations: [{"rfgain_sel": int, "fec_pct": float|None,
                    "snr_db": [floats]|None, "decoded": bool|None}, ...]

    ORDER OF EVIDENCE, strongest first -- because these are not
    interchangeable. A carrier that does not decode at all is not "a bit
    worse", it is off the air (RF33 at 2 read 0.0% FEC while the referee
    saw snq=100). So decode evidence outranks any margin number, and
    margin only breaks ties among settings that work.
      1. FEC percentage, when measured live;
      2. explicit decode success;
      3. median SNR.
    Ties go to the LOWER rfgain_sel index only when the metric is equal,
    since nothing else distinguishes them. -> (best, ranked_table)"""
    rows = []
    for o in observations:
        snr = o.get("snr_db")
        med = statistics.median(snr) if snr else None
        rows.append(dict(rfgain_sel=int(o["rfgain_sel"]),
                         fec_pct=o.get("fec_pct"),
                         decoded=o.get("decoded"),
                         snr_med=med))
    have_fec = [r for r in rows if r["fec_pct"] is not None]
    if have_fec:
        key = lambda r: (r["fec_pct"] if r["fec_pct"] is not None else -1,
                         -r["rfgain_sel"])
        ranked = sorted(have_fec, key=key, reverse=True)
        return ranked[0], ranked
    have_dec = [r for r in rows if r["decoded"] is not None]
    if have_dec and any(r["decoded"] for r in have_dec):
        ok = [r for r in have_dec if r["decoded"]]
        if any(r["snr_med"] is not None for r in ok):
            ranked = sorted(ok, key=lambda r: (r["snr_med"] or -1e9,
                                               -r["rfgain_sel"]), reverse=True)
            return ranked[0], ranked
        ranked = sorted(ok, key=lambda r: -r["rfgain_sel"], reverse=True)
        return ranked[0], ranked
    have_snr = [r for r in rows if r["snr_med"] is not None]
    if not have_snr:
        return None, rows
    ranked = sorted(have_snr, key=lambda r: (r["snr_med"], -r["rfgain_sel"]),
                    reverse=True)
    return ranked[0], ranked


def entry_from(rf, observations, method, evidence, when=None, ifgr=32,
               captures=None, notes=None):
    """Build a table entry WITH provenance from candidate observations."""
    best, ranked = choose(observations)
    if best is None:
        raise ValueError(f"RF{rf}: no usable observation to choose from")
    metric = {}
    fec = {r["rfgain_sel"]: r["fec_pct"] for r in ranked
           if r["fec_pct"] is not None}
    snr = {r["rfgain_sel"]: r["snr_med"] for r in ranked
           if r["snr_med"] is not None}
    if fec:
        metric = dict(kind="fec_pct_mean",
                      **{f"at_{k}": v for k, v in sorted(fec.items())})
    elif snr:
        metric = dict(kind="snr_db_median",
                      **{f"at_{k}": round(v, 3) for k, v in sorted(snr.items())})
        if len(snr) >= 2:
            vals = sorted(snr.values(), reverse=True)
            metric["delta_db"] = round(vals[0] - vals[1], 3)
    return {
        "rfgain_sel": best["rfgain_sel"],
        "ifgr": ifgr,
        "measured": True,
        "measured_utc": when or time.strftime("%Y-%m-%dT%H:%M"),
        "method": method,
        "evidence": evidence,
        "captures": captures or {},
        "metric": metric,
        "notes": notes or "",
    }


def calibrate_from_diags(rf, candidates, method, evidence, **kw):
    """Populate an entry from BANKED capture diagnostics.

    candidates: {rfgain_sel: path-to-*_diags.json}. Uses the median of the
    per-frame snr_db over frames that locked -- the same statistic the
    E76 A/B was decided on, so the table is reproducible from the bank."""
    obs = []
    for g, path in candidates.items():
        with open(path, encoding="utf-8") as f:
            d = json.load(f)
        snr = [x["snr_db"] for x in d.get("diags", [])
               if x.get("ok") and "snr_db" in x]
        obs.append(dict(rfgain_sel=int(g), snr_db=snr,
                        decoded=bool(snr) or None))
    return entry_from(rf, obs, method, evidence, **kw)


def save(table, path=None):
    with open(path or TABLE_PATH, "w", encoding="utf-8") as f:
        json.dump(table, f, indent=2)
        f.write("\n")


def put(rf, entry, path=None):
    """Write one entry into the table, keeping the rest untouched."""
    t = load(path)
    t.setdefault("carriers", {})[str(int(rf))] = entry
    save(t, path)
    return t


# ------------------------------------------------------------------ main

def main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rf", type=int, default=None)
    a = ap.parse_args()
    t = load()
    if a.rf is not None:
        print(describe(a.rf, t))
        return 0
    print(f"carrier gain table v{t.get('version')}  "
          f"default rfgain_sel={t['default']['rfgain_sel']} "
          f"IFGR={t['default']['ifgr']}")
    for rf in sorted(t.get("carriers", {}), key=int):
        print("  " + describe(int(rf), t))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
