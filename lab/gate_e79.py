#!/usr/bin/env python3
"""gate_e79.py -- the per-carrier gain table, gated against BANKED captures.

The table's claims are re-derived here from the evidence on disk rather
than trusted: the chooser is run over the same 8/10 diagnostics that
decided E73 and E76, and must land on the values the table ships.

NEGATIVE CONTROLS (house rule -- a gate that cannot fail proves nothing):
  1. ONE GLOBAL CONSTANT CANNOT SERVE BOTH CARRIERS. Whichever single
     rfgain_sel is chosen, it must be wrong for RF33 or for RF25 on the
     banked numbers. This is the whole justification for the table, so it
     is asserted, not assumed.
  2. A SILENT WRITE MUST BE CAUGHT. A fake radio that accepts the write
     and reports a different value back must make apply_gain raise; the
     compliant fake must not.
  3. An unmeasured carrier must NOT come back looking measured.

    python lab/gate_e79.py     -> PASS/FAIL, exit 0 iff all pass
"""
from __future__ import annotations

import json
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import carrier_gain as CG                                      # noqa: E402

FAILS = []


def check(name, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"  ({detail})" if detail
                                                     else ""))
    if not ok:
        FAILS.append(name)


def _diag_snr(path):
    with open(path, encoding="utf-8") as f:
        d = json.load(f)
    return [x["snr_db"] for x in d.get("diags", [])
            if x.get("ok") and "snr_db" in x]


def gate_rf25_from_bank():
    print("gate 1: RF25 -- re-derived from the banked A/B captures")
    g2 = os.path.join(HERE, "e75_g2_diags.json")
    g4 = os.path.join(HERE, "e75_g4_diags.json")
    if not (os.path.exists(g2) and os.path.exists(g4)):
        check("banked RF25 captures present", False, "missing diags")
        return
    s2, s4 = _diag_snr(g2), _diag_snr(g4)
    m2, m4 = statistics.median(s2), statistics.median(s4)
    check("both gains locked every frame", len(s2) == len(s4) == 150,
          f"{len(s2)} / {len(s4)} frames")
    check("gain 2 measures ~6.04 dB better than gain 4",
          abs((m2 - m4) - 6.049) < 0.05,
          f"{m2:.3f} - {m4:.3f} = {m2-m4:+.3f} dB")
    ent = CG.calibrate_from_diags(
        25, {2: g2, 4: g4},
        method="banked-capture A/B, 150 frames each",
        evidence="E75/E76", when="2026-08-10T22:57")
    check("the chooser picks rfgain_sel 2 from the evidence alone",
          ent["rfgain_sel"] == 2, f"chose {ent['rfgain_sel']}")
    check("the derived entry matches what the table ships",
          CG.lookup(25)["rfgain_sel"] == ent["rfgain_sel"] == 2)
    check("the entry carries provenance, not a bare integer",
          bool(ent.get("measured_utc")) and bool(ent.get("evidence"))
          and ent["metric"]["kind"] == "snr_db_median"
          and abs(ent["metric"]["delta_db"] - 6.049) < 0.05,
          f"delta {ent['metric'].get('delta_db')} dB")


def gate_rf33_from_bank():
    print("gate 2: RF33 -- re-derived from the live A/B FEC record")
    p = os.path.join(HERE, "e73_gain_live.json")
    if not os.path.exists(p):
        check("banked RF33 A/B present", False, "missing e73_gain_live.json")
        return
    with open(p, encoding="utf-8") as f:
        rows = json.load(f)
    obs = [dict(rfgain_sel=r["rfgain"], fec_pct=r["mean_all"]) for r in rows]
    best, ranked = CG.choose(obs)
    at2 = next(r["fec_pct"] for r in ranked if r["rfgain_sel"] == 2)
    check("rfgain_sel 2 is a TOTAL loss on RF33, not a degradation",
          at2 == 0.0, f"FEC {at2}% at rfgain_sel 2")
    check("the chooser picks a gain that decodes", best["fec_pct"] == 100.0,
          f"chose {best['rfgain_sel']} at {best['fec_pct']}%")
    check("table ships 4 (the lower of the two that work)",
          CG.lookup(33)["rfgain_sel"] == 4
          and {r["rfgain_sel"] for r in ranked if r["fec_pct"] == 100.0}
          == {4, 6})


def gate_no_global_constant():
    print("gate 3: NEGATIVE CONTROL -- no single global gain serves both")
    fec33 = {r["rfgain"]: r["mean_all"] for r in
             json.load(open(os.path.join(HERE, "e73_gain_live.json"),
                            encoding="utf-8"))}
    m2 = statistics.median(_diag_snr(os.path.join(HERE, "e75_g2_diags.json")))
    m4 = statistics.median(_diag_snr(os.path.join(HERE, "e75_g4_diags.json")))
    snr25 = {2: m2, 4: m4}
    verdicts = []
    for g in (2, 4):
        rf33_ok = fec33.get(g, 0.0) == 100.0
        cost25 = max(snr25.values()) - snr25[g]
        verdicts.append((g, rf33_ok, cost25))
    bad = [v for v in verdicts if not v[1] or v[2] > 1.0]
    check("every single global value is wrong for one of the two carriers",
          len(bad) == len(verdicts),
          "; ".join(f"g{g}: RF33 {'ok' if ok else 'DEAD'}, "
                    f"RF25 -{c:.2f} dB" for g, ok, c in verdicts))
    # and the table gets both right at once, which is the point
    check("the TABLE serves both carriers simultaneously",
          CG.lookup(33)["rfgain_sel"] == 4 and CG.lookup(25)["rfgain_sel"] == 2)


class FakeSdr:
    """Minimal Soapy stand-in. `liar` makes it accept a write and report
    something else -- the failure the readback law exists to catch."""

    def __init__(self, liar=False, ifgr_liar=False):
        self.liar, self.ifgr_liar = liar, ifgr_liar
        self.rfgain, self.ifgr = None, None

    def setGainMode(self, *a):
        pass

    def setGain(self, d, ch, name, val):
        self.ifgr = val

    def writeSetting(self, k, v):
        if k == "rfgain_sel":
            self.rfgain = int(v)

    def readSetting(self, k):
        if k != "rfgain_sel":
            return ""
        return str((self.rfgain or 0) + 1) if self.liar else str(self.rfgain)

    def getGain(self, d, ch, name):
        return float((self.ifgr or 0) + 5 if self.ifgr_liar else self.ifgr)


def gate_readback():
    print("gate 4: write-then-readback confirmation")
    quiet = lambda *a, **k: None
    ok = CG.apply_gain(FakeSdr(), 33, log=quiet)
    check("compliant radio: applies the table's value and confirms it",
          ok["rfgain_sel"] == 4 and ok["ifgr"] == 32, str(ok))

    # NEGATIVE CONTROL: the silent-write failure must be caught
    raised = False
    try:
        CG.apply_gain(FakeSdr(liar=True), 33, log=quiet)
    except CG.GainReadbackError:
        raised = True
    check("NEGATIVE CONTROL: rfgain_sel readback mismatch RAISES", raised)

    raised = False
    try:
        CG.apply_gain(FakeSdr(ifgr_liar=True), 33, log=quiet)
    except CG.GainReadbackError:
        raised = True
    check("NEGATIVE CONTROL: IFGR readback mismatch RAISES", raised)


    # an explicit override must WIN, be announced, and still be verified
    fs = FakeSdr()
    ok = CG.apply_gain(fs, 25, log=quiet, rfgain=4)
    check("explicit override beats the table and is still read back",
          ok["rfgain_sel"] == 4 and ok["source"] == "explicit"
          and CG.lookup(25)["rfgain_sel"] == 2, str(ok))
    raised = False
    try:
        CG.apply_gain(FakeSdr(liar=True), 25, log=quiet, rfgain=4)
    except CG.GainReadbackError:
        raised = True
    check("NEGATIVE CONTROL: an overridden write is verified too", raised)


def gate_unmeasured():
    print("gate 5: an unmeasured carrier must not look measured")
    g = CG.lookup(99)
    check("unknown carrier falls back to the documented default",
          g["rfgain_sel"] == 4 and g["measured"] is False
          and g["source"] == "default", str({k: g[k] for k in
                                             ("rfgain_sel", "measured",
                                              "source")}))
    check("its description says so out loud", "DEFAULT" in CG.describe(99)
          and "run calibration" in CG.describe(99))
    g30 = CG.lookup(30)
    check("RF30 decoded but gain unrecorded -> default, flagged distinctly",
          g30["measured"] is False
          and g30["source"] == "default-known-carrier"
          and "not recorded" in CG.describe(30), CG.describe(30))
    check("a measured carrier reports its provenance",
          "measured" in CG.describe(33) and "E73" in CG.describe(33))


if __name__ == "__main__":
    gate_rf25_from_bank()
    gate_rf33_from_bank()
    gate_no_global_constant()
    gate_readback()
    gate_unmeasured()
    print(f"\n{'ALL GATES PASS' if not FAILS else 'FAILED: ' + ', '.join(FAILS)}")
    sys.exit(1 if FAILS else 0)
