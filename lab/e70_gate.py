#!/usr/bin/env python3
"""E70 gate -- the locked-channel inspector and the EPG grid, with controls.

House law: a gate that cannot fail proves nothing.  Every leg below carries at
least one NEGATIVE control -- an input that SHOULD make the check fail -- and
the leg only passes if the positive passes AND the negative is caught.

    leg 1  classification    132.1 WHUT clear, 5.1/4.1/9.1 locked, from fields
    leg 2  field-driven      controls: mutate the SIGNALLING, not the name
    leg 3  never playable    control: a regressed decide_playable is caught
    leg 4  the UI cannot lie control: a play affordance on a locked row is caught
    leg 5  PICTURE GATE      E67 reused: flat_frac discriminates, "decoded" does not
    leg 6  captions          locked services yield cues from a CLEAR stpp track
    leg 7  no fabricated now control: outside the guide span, nothing is claimed
    leg 8  courtesy          the live soak was not disturbed

No radio.  Banked captures and read-only reads of the running chain's files.
"""
from __future__ import annotations

import copy
import json
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
TOOLS = os.path.join(ROOT, "tools")
for p in (TOOLS, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)

import atsc3_inspect as I          # noqa: E402
import atsc3_grid as G             # noqa: E402

OUT = os.path.join(HERE, "e70_out")
CHAIN = os.path.join(ROOT, "data", "e31", "chain.log")

# The four services E67 settled, by virtual channel.  105.1/104.1/109.1 in the
# task brief are the HDHomeRun's numbering; the SLT calls them 5.1/4.1/9.1.
CLEAR_TRUTH = {"32.1"}
LOCKED_TRUTH = {"5.1", "4.1", "9.1"}

legs = []


def leg(name, ok, detail, control=None):
    legs.append({"leg": len(legs) + 1, "name": name,
                 "verdict": "PASS" if ok else "FAIL",
                 "detail": detail, "negative_control": control})
    print("LEG %-2d %-22s %s" % (len(legs), name, "PASS" if ok else "FAIL"))
    for line in detail.splitlines():
        print("        " + line)
    if control:
        for line in control.splitlines():
            print("        control: " + line)
    return ok


def svc_by_chan(model, ch):
    for c in model["carriers"].values():
        for s in c["services"]:
            if s.get("vchannel") == ch:
                return s
    return None


def soak():
    v = I.soak_instantaneous(CHAIN)
    return v


# ------------------------------------------------------------------- run

def main():
    I.be_polite()
    os.makedirs(OUT, exist_ok=True)
    inst_before = soak()
    if inst_before is not None and inst_before < 0.95:
        print("STAND DOWN: live chain inst %.2fx < 0.95x" % inst_before)
        return 2

    model = I.build_lineup()
    allok = True

    # ---- leg 1: classification -------------------------------------------
    got = {}
    for ch in sorted(CLEAR_TRUTH | LOCKED_TRUTH):
        s = svc_by_chan(model, ch)
        got[ch] = None if s is None else s["protection"]["state"]
    ok = all(got[ch] == "CLEAR" for ch in CLEAR_TRUTH) and \
        all(got[ch] == "PROTECTED" for ch in LOCKED_TRUTH)
    # and the deepest layer -- the bytes -- must agree with the signalling.
    # E67 measured 98.3% of the video payload protected on the locked three;
    # if this drifts, either the reader broke or the broadcaster changed.
    senc = {}
    for ch in sorted(CLEAR_TRUTH | LOCKED_TRUTH):
        s = svc_by_chan(model, ch).get("_senc") or {}
        senc[ch] = None if s.get("absent", True) else s.get("protected_frac")
    ok = ok and senc["32.1"] is None and \
        all((senc[ch] or 0) > 0.90 for ch in LOCKED_TRUTH)
    fields = {}
    for ch in sorted(CLEAR_TRUTH | LOCKED_TRUTH):
        s = svc_by_chan(model, ch)
        fields[ch] = sorted({"%s/%s" % (e["layer"], e["field"])
                             for e in s["protection"]["evidence"]})[:4]
    allok &= leg(
        "classification",
        ok,
        "verdicts: " + ", ".join("%s=%s(%d layers)" %
                                 (ch, got[ch],
                                  svc_by_chan(model, ch)["protection"]["layers_read"])
                                 for ch in sorted(got)) +
        "\nWHUT read from: " + ", ".join(fields["32.1"]) +
        "\nWTTG read from: " + ", ".join(fields["5.1"]) +
        "\nsenc protected fraction (layer 4, the bytes): " +
        ", ".join("%s=%s" % (ch, "absent" if senc[ch] is None
                             else "%.1f%%" % (100 * senc[ch]))
                  for ch in sorted(senc)),
        "leg 2 mutates those fields and requires the verdicts to move")

    # ---- leg 2: the verdict follows the FIELDS, not the callsign ----------
    # A hardcoded "WTTG is locked" table would pass leg 1 and fail here.
    wttg = copy.deepcopy(svc_by_chan(model, "5.1"))
    whut = copy.deepcopy(svc_by_chan(model, "32.1"))

    # control A -- rename the locked service to the clear one's callsign.
    ren = copy.deepcopy(wttg)
    ren["name"] = "WHUT"
    ren["vchannel"] = "32.1"
    a_state = I.classify_protection(ren)["state"]

    # control B -- strip every protection field from the locked service and
    # remove its media evidence: it must fall back to CLEAR.
    stripped = copy.deepcopy(wttg)
    stripped["_slt_protected_attr"] = None
    stripped["_slt_drmSystemID"] = None
    if stripped.get("_sls"):
        stripped["_sls"]["content_protection"] = []
    stripped["_inits"] = {}
    stripped["_senc"] = None
    b_state = I.classify_protection(stripped)["state"]

    # control C -- plant @protected on the CLEAR service: it must go locked.
    planted = copy.deepcopy(whut)
    planted["_slt_protected_attr"] = "true"
    planted["_slt_drmSystemID"] = "urn:uuid:" + I.WIDEVINE_UUID
    c_state = I.classify_protection(planted)["state"]

    # control D -- an SLT we never read at all is UNKNOWN, not CLEAR.
    blind = copy.deepcopy(whut)
    blind["slt_source"] = None
    blind["_sls"] = None
    blind["_inits"] = {}
    blind["_senc"] = None
    d_state = I.classify_protection(blind)["state"]

    ok2 = (a_state == "PROTECTED" and b_state == "CLEAR" and
           c_state == "PROTECTED" and d_state == "UNKNOWN")
    allok &= leg(
        "field-driven verdict", ok2,
        "the verdict tracks the signalling, not the callsign",
        "A rename WTTG->WHUT           -> %s (must stay PROTECTED)\n"
        "B strip protection fields     -> %s (must become CLEAR)\n"
        "C plant @protected on WHUT    -> %s (must become PROTECTED)\n"
        "D no signalling read at all   -> %s (must be UNKNOWN, never CLEAR)"
        % (a_state, b_state, c_state, d_state))

    # ---- leg 3: a protected service is NEVER playable ---------------------
    bad = [s["vchannel"] for c in model["carriers"].values()
           for s in c["services"]
           if s["protection"]["state"] == "PROTECTED" and
           s["playable"]["playable"]]
    # negative control: the regression this leg exists to catch -- a
    # decide_playable that consults media availability BEFORE protection.
    def regressed(svc):
        if svc.get("_has_media"):
            return {"playable": True, "reason": "CLEAR", "detail": "regressed"}
        return I.decide_playable(svc)
    caught = regressed(wttg)["playable"] is True     # the mutant DOES lie
    ok3 = (not bad) and caught
    allok &= leg(
        "protected is never playable", ok3,
        "protected services claiming playable: %s" % (bad or "none"),
        "a decide_playable() that checks media before protection returns "
        "playable=%s for WTTG -- this leg would fail against it, so the leg "
        "can fail" % regressed(wttg)["playable"])

    # ---- leg 4: the rendered UI cannot present a locked service as watchable
    g = G.grid_json(at="guide")
    static = G.render_static(g, os.path.join(OUT, "grid_guide.html"))
    G.render_static(G.grid_json(), os.path.join(OUT, "grid_now.html"))
    page = open(static, encoding="utf-8").read()
    locked_rows = [r for r in g["rows"] if r["state"] == "PROTECTED"]
    badge_ok = all(("ENCRYPTED" in r["badge"] and not r["playable"])
                   for r in locked_rows)
    # every locked row must actually appear in the HTML, badge and all
    rendered = all(re.search(
        r'<span class="vch">%s</span>.{0,400}ENCRYPTED' % re.escape(r["vchannel"]),
        page, re.S) for r in locked_rows)
    # ...and no play AFFORDANCE.  Scan the interactive surface -- links,
    # buttons and event handlers -- not the prose.  The first version of this
    # check grepped the whole page for /decrypt/ and flagged the page's own
    # sentence "Nothing here attempts or assists decryption", which is the
    # opposite of an offence.  An affordance is something you can click.
    offenders = _affordances(page)
    # and a locked row must carry no clickable element at all
    row_clickable = [r["vchannel"] for r in locked_rows
                     if _row_html(page, r["vchannel"]) and
                     re.search(r"<a\b|<button\b|onclick=",
                               _row_html(page, r["vchannel"]), re.I)]
    # negative control: the same scan over a page that DOES offer one
    planted = page.replace("</tbody>",
                           '<tr><td><a href="#" onclick="play()">Watch now'
                           '</a></td></tr></tbody>')
    control_catches = bool(_affordances(planted))
    ok4 = (badge_ok and rendered and not offenders and not row_clickable
           and control_catches)
    allok &= leg(
        "the UI cannot lie", ok4,
        "%d locked rows rendered with the ENCRYPTED badge; play affordances "
        "in the page: %s; clickable elements inside a locked row: %s" %
        (len(locked_rows), offenders or "none", row_clickable or "none"),
        "the same scan over a page with a planted 'Watch now' link finds it: "
        "%s" % control_catches)

    # ---- leg 5: THE PICTURE GATE (E67, reused not rewritten) --------------
    ok5, pic_detail, pic_control = picture_leg()
    allok &= leg("picture gate (E67)", ok5, pic_detail, pic_control)

    # ---- leg 6: captions off a locked service, from a CLEAR track ---------
    cap = {}
    for ch in sorted(LOCKED_TRUTH | CLEAR_TRUTH):
        s = svc_by_chan(model, ch)
        text_entry = None
        for tsi, info in (s.get("_inits") or {}).items():
            if info.get("kind") == "text":
                text_entry = info.get("sample_entry")
        video_entry = None
        for tsi, info in (s.get("_inits") or {}).items():
            if info.get("kind") == "video":
                video_entry = info.get("sample_entry")
        cap[ch] = {"cues": len(s.get("captions") or []),
                   "text_entry": text_entry, "video_entry": video_entry,
                   "status": s.get("captions_status")}
    # the product claim: a LOCKED service yields readable captions...
    have = [ch for ch in LOCKED_TRUTH if cap[ch]["cues"] > 0]
    # ...from a track that is genuinely NOT encrypted, on a service whose
    # video track genuinely IS.  Both halves, or the claim is hollow.
    clean_text = all(cap[ch]["text_entry"] == "stpp" for ch in LOCKED_TRUTH)
    locked_video = all(cap[ch]["video_entry"] == "encv" for ch in LOCKED_TRUTH)
    ok6 = bool(have) and clean_text and locked_video
    allok &= leg(
        "captions of a locked service", ok6,
        "cues recovered on %s; per-service tracks: %s" %
        (", ".join(sorted(have)) or "none",
         ", ".join("%s text=%s video=%s cues=%d" %
                   (ch, cap[ch]["text_entry"], cap[ch]["video_entry"],
                    cap[ch]["cues"]) for ch in sorted(cap))),
        "the same services' VIDEO sample entries are encv (encrypted), so the "
        "caption text cannot be coming from an accidentally-clear service: "
        "locked_video=%s, text_all_stpp=%s" % (locked_video, clean_text))

    # ---- leg 7: nothing is invented outside the guide's span --------------
    span = model["esg"]["span"]
    far = span[1] + 30 * 86400
    g_far = G.grid_json(at=str(far))
    filled_far = sum(1 for r in g_far["rows"] for c in r["cells"]
                     if not c.get("empty"))
    g_in = G.grid_json(at="guide")
    filled_in = sum(1 for r in g_in["rows"] for c in r["cells"]
                    if not c.get("empty"))
    ok7 = filled_far == 0 and filled_in > 0 and g_far["status"] == "stale"
    allok &= leg(
        "no fabricated now/next", ok7,
        "inside the guide window %d cells are filled; 30 days past its end "
        "%d are, and the status is '%s'" %
        (filled_in, filled_far, g_far["status"]),
        "the positive half is the control: a grid that always says something "
        "would show filled cells in BOTH cases")

    # ---- leg 8: an INDEPENDENT receiver agrees ---------------------------
    ok8a, ref_detail, ref_control = referee_leg(model)
    allok &= leg("independent receiver", ok8a, ref_detail, ref_control)

    # ---- leg 9: courtesy ---------------------------------------------------
    inst_after = soak()
    ok8 = (inst_after is None or inst_after >= 0.95)
    allok &= leg(
        "soak undisturbed", ok8,
        "live chain inst %sx before, %sx after this gate" %
        (inst_before, inst_after),
        "the leg fails below 0.95x, which is the threshold the task set")

    print()
    print("E70 GATE: %s  (%d legs)" % ("ALL PASS" if allok else "FAILURES",
                                       len(legs)))
    rep = {"build": I.BUILD, "grid_build": G.BUILD, "when": time.time(),
           "all_pass": bool(allok), "legs": legs}
    with open(os.path.join(HERE, "e70_gate.json"), "w") as fh:
        json.dump(rep, fh, indent=2)
    return 0 if allok else 1


HDHR = os.environ.get("ATSC3_HDHR_HOST", "")  # tuner IP via env; referee leg skips when unset


def referee_leg(model):
    """A commercial receiver, on the same antenna, answering independently.

    Our four layers are all readings of documents WE parsed.  The HDHomeRun
    FLEX 4K parsed the same air with somebody else's code and publishes its
    own verdict as a `DRM` flag in /lineup.json.  If it disagrees with us,
    one of us is wrong and we want to know.

    Strictly read-only HTTP: /lineup.json is the stored scan result, so no
    tuner is engaged.  tuner1 is the live soak's referee and is NOT touched --
    the leg records its state to prove that.

    E65's rule applies to this leg too: an unreachable instrument SKIPS.  A
    missing instrument is never a verdict, in either direction.
    """
    import urllib.request
    try:
        with urllib.request.urlopen("http://%s/lineup.json" % HDHR,
                                    timeout=8) as r:
            lineup = json.load(r)
        with urllib.request.urlopen("http://%s/status.json" % HDHR,
                                    timeout=8) as r:
            status = json.load(r)
    except Exception as exc:
        return True, ("SKIP -- the FLEX 4K did not answer (%r). A missing "
                      "instrument is never a verdict." % exc), \
            "nothing was concluded from its silence"
    with open(os.path.join(OUT, "e70_hdhr_lineup.json"), "w") as fh:
        json.dump({"lineup": lineup, "status": status}, fh, indent=2)

    theirs = {}
    for e in lineup:
        theirs[(e.get("GuideName") or "").upper()] = bool(e.get("DRM"))
    ours, rows = {}, []
    for c in model["carriers"].values():
        for s in c["services"]:
            nm = (s.get("name") or "").upper()
            if nm not in theirs:
                continue
            mine = s["protection"]["state"] == "PROTECTED"
            ours[nm] = mine
            rows.append("%-6s ours=%-9s theirs DRM=%d %s" %
                        (nm, s["protection"]["state"], theirs[nm],
                         "agree" if mine == theirs[nm] else "DISAGREE"))
    agree = all(ours[n] == theirs[n] for n in ours)
    # and the leg must be able to notice a disagreement
    mutant = dict(theirs)
    for k in mutant:
        mutant[k] = not mutant[k]
    would_catch = not all(ours[n] == mutant[n] for n in ours)
    t1 = [t for t in status if t.get("Resource") == "tuner1"]
    ok = agree and bool(ours) and would_catch
    detail = ("%d services matched by callsign\n" % len(ours)) + \
        "\n".join(sorted(rows)) + \
        "\ntuner1 (the soak's referee) left alone: %s" % (t1[0] if t1 else "?")
    control = ("inverting the FLEX's DRM flags makes this leg fail (%s), so "
               "agreement is a result and not a tautology; and WHUT is the "
               "case that matters -- the FLEX scan calls it '(no data)' yet "
               "flags it DRM=0, which is E67's finding reproduced by somebody "
               "else's receiver" % would_catch)
    return ok, detail, control


PLAY_WORDS = re.compile(
    r"(?i)(watch|\bplay\b|\btune\b|decrypt|unlock|widevine|"
    r"\b(vlc|mpv|rtp|udp|rtsp)://)")


def _affordances(page):
    """Every clickable thing on the page that offers playback.

    Looks only at anchors, buttons and inline handlers -- the things a user can
    act on -- so the page is free to EXPLAIN encryption in prose without the
    gate mistaking the explanation for an offer.
    """
    hits = []
    for m in re.finditer(r"<a\b[^>]*>(.*?)</a>|<button\b[^>]*>(.*?)</button>",
                         page, re.S | re.I):
        frag = m.group(0)
        if PLAY_WORDS.search(frag):
            hits.append(re.sub(r"\s+", " ", frag)[:80])
    for m in re.finditer(r"(?:href|onclick|onmousedown)\s*=\s*"
                         r"[\"']([^\"']*)[\"']", page, re.I):
        if PLAY_WORDS.search(m.group(1)):
            hits.append(m.group(0)[:80])
    return hits


def _row_html(page, vchannel):
    m = re.search(r"<tr[^>]*>(?:(?!</tr>).)*?<span class=\"vch\">%s</span>"
                  r"(?:(?!</tr>).)*?</tr>" % re.escape(vchannel), page, re.S)
    return m.group(0) if m else ""


# ------------------------------------------------------------ picture leg

def picture_leg():
    """E67's picture gate, re-run and re-argued.

    Two things must hold:
      POSITIVE  flat_frac still separates the clear service's picture from the
                locked services' concealment, by the margin E67 measured.
      NEGATIVE  the checks that DID NOT discriminate still do not -- ffprobe
                parses all four, and mdat entropy is ~8 bits/byte on all four.
                If a future contributor swaps flat_frac for one of those, this
                leg must fail.  And this tool renders no thumbnail at all.
    """
    import e67_picture_gate as PG
    clips = {"whut": ("32.1", "clear"), "wttg": ("5.1", "locked"),
             "wrc": ("4.1", "locked"), "wusa": ("9.1", "locked")}
    e67out = os.path.join(HERE, "e67_out")
    stats, probe = {}, {}
    for name in clips:
        mp4 = os.path.join(e67out, "%s_video.mp4" % name)
        if not os.path.isfile(mp4):
            return False, "missing %s" % mp4, None
        stats[name] = PG.picture_stats(mp4)
        probe[name] = _ffprobe(mp4)
    clear = stats["whut"]["flat_frac_mean"]
    locked = [stats[n]["flat_frac_mean"] for n in ("wttg", "wrc", "wusa")]
    positive = clear < 0.10 and all(v > 0.80 for v in locked)

    # the non-discriminators, still not discriminating
    all_parse = all(p.get("codec") == "hevc" for p in probe.values())
    all_same_size = len({p.get("size") for p in probe.values()}) == 1
    temporal = [stats[n]["temporal_r_mean"] for n in ("wttg", "wrc", "wusa")]
    temporal_useless = max(temporal) > 0.80        # garbage correlates fine
    negative = all_parse and all_same_size and temporal_useless

    # and the inspector refuses to claim a picture at all
    claim = I.picture_claim({})
    refuses = claim["claims_picture"] is False
    g = G.grid_json(at="guide")
    no_thumbs = g["picture_policy"]["renders_thumbnails"] is False
    page = open(os.path.join(OUT, "grid_guide.html"), encoding="utf-8").read()
    no_img = "<img" not in page.lower() and "data:image" not in page.lower()

    ok = positive and negative and refuses and no_thumbs and no_img
    detail = ("flat_frac  WHUT %.4f (clear)  vs  WTTG %.4f / WRC %.4f / "
              "WUSA %.4f (locked) -- %.0fx separation\n"
              "grid renders thumbnails: %s;  <img> in page: %s;  "
              "picture_claim(): %s"
              % (clear, stats["wttg"]["flat_frac_mean"],
                 stats["wrc"]["flat_frac_mean"], stats["wusa"]["flat_frac_mean"],
                 min(locked) / clear if clear else 0,
                 g["picture_policy"]["renders_thumbnails"],
                 "<img" in page.lower(), claim["claims_picture"]))
    control = ("the metrics that LIED in E67 still lie, so swapping them back "
               "in would fail this leg:\n"
               "  ffprobe reports hevc on all four: %s\n"
               "  identical frame size on all four: %s\n"
               "  temporal_r on locked clips: %s (max %.3f > 0.80, so it "
               "cannot discriminate)"
               % (all_parse, all_same_size,
                  ["%.3f" % t for t in temporal], max(temporal)))
    return ok, detail, control


def _ffprobe(mp4):
    try:
        p = subprocess.run(
            ["ffprobe", "-v", "error", "-select_streams", "v:0",
             "-show_entries", "stream=codec_name,width,height",
             "-of", "json", mp4], capture_output=True, timeout=120)
        j = json.loads(p.stdout or b"{}")
        s = (j.get("streams") or [{}])[0]
        return {"codec": s.get("codec_name"),
                "size": "%sx%s" % (s.get("width"), s.get("height"))}
    except Exception as exc:
        return {"error": repr(exc)}


if __name__ == "__main__":
    sys.exit(main())
