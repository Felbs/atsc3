#!/usr/bin/env python3
"""E81 gate -- the head-to-head scoreboard, and the controls that can fail it.

House law: a gate that cannot fail proves nothing.  The specific risk with a
scoreboard is not that it crashes -- it is that it quietly flatters us.  So
most of these legs are negative controls first and assertions second: feed the
harness a case we KNOW we lose and require it to say LOSE.

    leg 1  reports a LOSS when we lose      (the coordinator's own ask)
    leg 2  reports a WIN when we win, and a TIE when neither can present
    leg 3  null FEC is OUR downtime, never a skipped row
    leg 4  a threshold-dependent verdict is swept, and instability downgrades
    leg 5  D3 stays NOT-COMPARABLE and cannot be talked into a number
    leg 6  the reserved tuner is never retuned, even when the probe throws
    leg 7  carrier-scoped matching resists the cross-standard name collision
    leg 8  the FLEX's own L1 parser corroborates our claims (unscored)
    leg 9  courtesy: the soak was not disturbed

No radio.  The only device traffic is read-only status polls plus, in leg 6,
a deliberately failing borrow against a FAKE tuner id so nothing real moves.
"""
from __future__ import annotations

import copy
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
TOOLS = os.path.join(ROOT, "tools")
for p in (TOOLS, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)

import atsc3_inspect as I     # noqa: E402
import atsc3_vs as V          # noqa: E402

OUT = os.path.join(HERE, "e81_out")
CHAIN = os.path.join(ROOT, "data", "e31", "chain.log")

legs = []


def leg(name, ok, detail, control=None):
    legs.append({"leg": len(legs) + 1, "name": name,
                 "verdict": "PASS" if ok else "FAIL",
                 "detail": detail, "negative_control": control})
    print("LEG %-2d %-28s %s" % (len(legs), name, "PASS" if ok else "FAIL"))
    for line in detail.splitlines():
        print("        " + line)
    if control:
        for line in control.splitlines():
            print("        control: " + line)
    return ok


def main():
    I.be_polite()
    os.makedirs(OUT, exist_ok=True)
    inst_before = I.soak_instantaneous(CHAIN)
    if inst_before is not None and inst_before < 0.95:
        print("STAND DOWN: live chain inst %.2fx < 0.95x" % inst_before)
        return 2

    probe = V.load_probe()
    sb = V.scoreboard()
    allok = True
    dims = {d["dimension"]: d for d in sb["dimensions"]}

    # ---- leg 1: it says LOSE when we lose --------------------------------
    # Acquisition is the case we KNOW we lose: ~13 s against ~1.8 s.
    acq = dims["acquisition"]
    real_lose = acq["verdict"] == V.LOSE
    # and the machinery is not hardcoded to say LOSE there: hand it a world in
    # which we acquire faster and require the verdict to flip.
    fast = V.d4_acquisition(
        {"source": "synthetic", "cold_starts": 10, "window_s": 120,
         "no_lock_within_window": 0, "no_lock_pct": 0.0,
         "start_to_bootstrap": {"n": 10, "median_s": 0.2, "min_s": 0.2,
                                "max_s": 0.2, "p90_s": 0.2},
         "start_to_first_mpu": {"n": 10, "median_s": 0.4, "min_s": 0.4,
                                "max_s": 0.4, "p90_s": 0.4}},
        probe)
    wnuv = [r for r in dims["coverage"]["detail"]["rows"]
            if r["name"] == "WNUV"]
    ok1 = real_lose and fast["verdict"] == V.WIN and \
        bool(wnuv) and wnuv[0]["verdict"] == V.LOSE
    allok &= leg(
        "reports a LOSS when we lose", ok1,
        "acquisition: ours %s vs theirs %s -> %s\n"
        "coverage: WNUV 54.1 (the FLEX plays an RF25 MMTP service we have "
        "not decoded) -> %s"
        % (acq["ours"], acq["theirs"], acq["verdict"],
           wnuv[0]["verdict"] if wnuv else "MISSING"),
        "the same function fed a synthetic 0.4 s acquisition returns %s, so "
        "the LOSE above is a reading and not a constant" % fast["verdict"])

    # ---- leg 2: WIN and TIE are equally reachable -------------------------
    cov = dims["coverage"]["detail"]["rows"]
    whut = [r for r in cov if r["name"] == "WHUT"]
    locked = [r for r in cov if r["ours"] == "LOCKED-DRM"]
    ok2 = (bool(whut) and whut[0]["verdict"] == V.WIN
           and bool(locked) and all(r["verdict"] == V.TIE for r in locked))
    allok &= leg(
        "WIN and TIE are both reachable", ok2,
        "WHUT 32.1: ours %s / theirs %s -> %s\n"
        "%d DRM-locked services, all -> TIE (neither side may decrypt them)"
        % (whut[0]["ours"], whut[0]["theirs"], whut[0]["verdict"],
           len(locked)),
        "an earlier revision ranked LOCKED-DRM BELOW '(no data)' and scored "
        "every Widevine row as a LOSS for us; presentability comparison is "
        "what makes them ties")

    # ---- leg 3: a null FEC sample is our downtime -------------------------
    up_null = V._ours_up({"fec_5min": None}, 95.0)
    up_low = V._ours_up({"fec_5min": 10.0}, 95.0)
    up_good = V._ours_up({"fec_5min": 100.0}, 95.0)
    rows, _ = V.paired_samples()
    nulls = sum(1 for r in rows if r.get("fec_5min") is None)
    c_all = V.contingency(rows, 95.0, 80)
    c_drop = V.contingency([r for r in rows if r.get("fec_5min") is not None],
                           95.0, 80)
    ok3 = (up_null is False and up_low is False and up_good is True
           and nulls > 0 and c_drop["ours_up_pct"] > c_all["ours_up_pct"])
    allok &= leg(
        "null FEC counts as OUR downtime", ok3,
        "%d of %d paired samples have a null FEC (our chain had nothing to "
        "report)\nkept: ours %.2f%% up" % (nulls, len(rows),
                                           c_all["ours_up_pct"]),
        "dropping those rows instead would report ours %.2f%% up -- "
        "+%.2f points of pure self-flattery, which is exactly the bug this "
        "leg exists to prevent"
        % (c_drop["ours_up_pct"],
           c_drop["ours_up_pct"] - c_all["ours_up_pct"]))

    # ---- leg 4: threshold sweep, and instability downgrades ---------------
    d2 = dims["concurrent delivery"]
    swept = d2["detail"]["sweep"]
    stable = d2["detail"]["threshold_stable"]
    # negative control: a synthetic dataset engineered to flip with the
    # threshold must come back NOT-COMPARABLE rather than picking a winner.
    synth = []
    for i in range(40):
        synth.append({"t": 1786400000 + i * 300, "fec_5min": 92.0,
                      "hdhr": {"lock": "atsc3", "snq": 85}})
    d2s = V.d2_concurrent(synth, ["synthetic"])
    ok4 = (len(swept) == len(V.FEC_SWEEP) * len(V.SNQ_SWEEP)
           and stable and d2s["verdict"] == V.NOT_COMPARABLE
           if hasattr(V, "NOT_COMPARABLE") else
           (len(swept) == len(V.FEC_SWEEP) * len(V.SNQ_SWEEP) and stable
            and d2s["verdict"] == V.NC))
    allok &= leg(
        "threshold sweep and downgrade", ok4,
        "%d threshold combinations swept; the real verdict is %s across all "
        "of them" % (len(swept), "STABLE" if stable else "UNSTABLE"),
        "a synthetic set at FEC 92 / snq 85 -- which flips as the cut-offs "
        "cross it -- returns %s instead of a winner" % d2s["verdict"])

    # ---- leg 5: D3 refuses to be scored -----------------------------------
    d3 = V.d3_sensitivity()
    ok5 = (d3["verdict"] == V.NC and len(d3["caveats"]) >= 3
           and dims["sensitivity / threshold"]["verdict"] == V.NC)
    allok &= leg(
        "sensitivity stays NOT-COMPARABLE", ok5,
        "verdict %s with %d stated reasons" % (d3["verdict"],
                                               len(d3["caveats"])),
        "there is no code path that assigns D3 a WIN/LOSE/TIE: the function "
        "returns the constant, so no data can talk it into a number")

    # ---- leg 6: the reserved tuner survives an exception ------------------
    calls = []
    real_hd = V._hd

    def spy(*a, **k):
        calls.append(a)
        return ""
    V._hd = spy
    fake = "99"
    try:
        with V.borrowed(fake):
            raise RuntimeError("probe blew up mid-measurement")
    except RuntimeError:
        pass
    finally:
        V._hd = real_hd
    restored = any(a[:2] == ("set", "/tuner%s/channel" % fake) and
                   a[2] == "none" for a in calls if len(a) >= 3)
    touched_reserved = any(V.RESERVED_TUNER in str(a) and "set" in str(a)
                           for a in calls)
    live_t1 = V.tuner_status(V.RESERVED_TUNER)
    ok6 = (restored and not touched_reserved and live_t1 is not None
           and live_t1.get("lock") not in (None, "none"))
    allok &= leg(
        "reserved tuner never retuned", ok6,
        "borrow() restored tuner%s to none even though the body raised; "
        "tuner%s (the soak's referee) is %s lock=%s snq=%s"
        % (fake, V.RESERVED_TUNER, (live_t1 or {}).get("ch"),
           (live_t1 or {}).get("lock"), (live_t1 or {}).get("snq")),
        "the exception path is the control -- an ordinary return would hide a "
        "missing restore; no `set` was issued against the reserved tuner: %s"
        % (not touched_reserved))

    # ---- leg 7: the cross-standard name collision stays fixed -------------
    # "WBFF" names BOTH the ATSC 3.0 service on RF25 and an unrelated ATSC 1.0
    # MPEG2 broadcast in lineup.json. Name matching scored the 3.0 service as
    # PLAYS on the strength of the 1.0 one.
    svc_wbff = {"vchannel": "45.1", "name": "WBFF"}
    got = V.theirs_state(svc_wbff, "1408", 25, probe)
    lineup_names = [e.get("GuideName") for e in (probe.get("lineup") or [])
                    if (e.get("GuideName") or "").upper().startswith("WBFF")]
    whut_state = V.theirs_state({"vchannel": "32.1", "name": "WHUT"},
                                "540", 33, probe)
    ok7 = got == "NO-DATA" and whut_state == "NO-DATA" and \
        V.flex_vchannel("45.1") == "145.1"
    allok &= leg(
        "carrier-scoped matching", ok7,
        "WBFF 45.1 (ATSC 3.0, RF25) -> %s;  WHUT 32.1 -> %s;  "
        "our 45.1 maps to the FLEX's %s"
        % (got, whut_state, V.flex_vchannel("45.1")),
        "lineup.json also contains %s (ATSC 1.0 MPEG2); name matching "
        "returned PLAYS from those rows, and lineup presence alone credited "
        "the FLEX with WHUT against its own streaminfo" % lineup_names)

    # ---- leg 8: the independent parser agrees (unscored) ------------------
    cor = sb.get("corroboration") or []
    checked = [c for c in cor if c["agrees"] is not None]
    agree = [c for c in checked if c["agrees"]]
    ok8 = bool(checked) and len(agree) == len(checked)
    allok &= leg(
        "independent L1 corroboration", ok8,
        "%d/%d of our own claims confirmed by the FLEX's L1 parser:\n%s"
        % (len(agree), len(checked),
           "\n".join("  RF%s %s: ours %s, flex %s"
                     % (c["rf"], c["check"], c["ours_claim"],
                        c["flex_reports"]) for c in checked)),
        "a DISAGREES row would fail this leg; the checks are on values we "
        "published BEFORE reading plpinfo (E67 sfi=1/qam256, E69 and E76 "
        "enhanced layer), so they were falsifiable when written")

    # ---- leg 9: courtesy ---------------------------------------------------
    inst_after = I.soak_instantaneous(CHAIN)
    ok9 = (inst_after is None or inst_after >= 0.95)
    allok &= leg("soak undisturbed", ok9,
                 "live chain inst %sx before, %sx after" %
                 (inst_before, inst_after),
                 "fails below 0.95x")

    print()
    print("E81 GATE: %s  (%d legs)" % ("ALL PASS" if allok else "FAILURES",
                                       len(legs)))
    json.dump({"build": V.BUILD, "when": time.time(),
               "all_pass": bool(allok), "tally": sb["tally"], "legs": legs},
              open(os.path.join(HERE, "e81_gate.json"), "w"), indent=2)
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
