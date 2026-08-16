#!/usr/bin/env python3
"""atsc3_inspect -- the LOCKED-CHANNEL INSPECTOR  (E70, task #65)

Point it at a live directory or a banked ROUTE/MMTP object directory and it
reports, per service:

    name + virtual channel        transport (ROUTE/MMTP) + PLP
    PROTECTION STATE, and the FIELD IT WAS READ FROM -- never inferred
    resolution / codec / bitrate / audio languages
    current + next programme, from the broadcaster's own service guide
    a tail of the live captions

It works IDENTICALLY on clear and on protected services.  The only difference
is the `playable` verdict: a protected service is LOCKED and its media is
never presented as watchable.  Everything the broadcaster sends in the clear
-- the guide, the manifest, the track list, and the (unencrypted) subtitles --
is surfaced for both.

DESIGN RULES THIS FILE ENFORCES
-------------------------------
1.  **Protection is READ, never inferred.**  Every protection verdict carries
    a list of `Evidence` records naming the layer, the field, the value and
    the file it came from.  Four independent layers are read:

        L1  SLT       Service/@protected, Service/@drmSystemID
        L2  MPD       ContentProtection/@schemeIdUri, cenc:default_KID
        L3  init seg  stsd sample entry (encv/enca), sinf/schm, tenc
        L4  bytes     senc subsample map -> protected byte fraction

    Absence is evidence too, but only when the carrying document was actually
    parsed: A/331 defaults @protected to false, so "attribute absent in an SLT
    we read" is a real negative.  "We never found an SLT" is UNKNOWN.

2.  **`playable` can never be true for a protected service.**  There is one
    function, `decide_playable()`, it is the only place the field is set, and
    it short-circuits on protection before it looks at anything else.
    `lab/e70_gate.py` holds a negative control against exactly this.

3.  **"The decoder ran" is not "the decoder decoded"** (E67's law).  This tool
    never claims a picture.  If a caller wants a thumbnail it must go through
    `lab/e67_picture_gate.py` -- flat_frac, not frame count.  See
    `picture_claim()`, which refuses by construction.

4.  No circumvention, ever.  `default_KID` is an IDENTIFIER and is printed as
    provenance for the protection verdict; no key material is involved and
    none is sought.

USAGE
    python tools/atsc3_inspect.py                       # default banked sources
    python tools/atsc3_inspect.py --source data/e31     # a live directory
    python tools/atsc3_inspect.py --service WTTG        # one service, verbose
    python tools/atsc3_inspect.py --json                # the whole model
"""
from __future__ import annotations

import argparse
import datetime as dt
import glob
import gzip
import io
import json
import os
import re
import struct
import sys
import xml.etree.ElementTree as ET
from collections import OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
LAB = os.path.join(ROOT, "lab")
if LAB not in sys.path:
    sys.path.insert(0, LAB)

BUILD = "E70-2026-08-10"

# NTP epoch: the service guide's PresentationWindow times are NTP seconds.
NTP_EPOCH = 2208988800

# ------------------------------------------------------------------ courtesy


def be_polite():
    """Drop to BelowNormal.  A live 6-hour soak owns this machine."""
    try:
        import ctypes
        # BELOW_NORMAL_PRIORITY_CLASS
        ctypes.windll.kernel32.SetPriorityClass(
            ctypes.windll.kernel32.GetCurrentProcess(), 0x00004000)
        return "BelowNormal"
    except Exception:
        try:
            os.nice(10)
            return "nice+10"
        except Exception:
            return "unchanged"


def soak_instantaneous(chain_log):
    """Last `inst` figure from a chain log, or None.  Used to stand down."""
    try:
        with open(chain_log, "rb") as fh:
            fh.seek(0, io.SEEK_END)
            fh.seek(max(0, fh.tell() - 8192))
            tail = fh.read().decode("utf-8", "replace")
    except Exception:
        return None
    hits = re.findall(r"inst\s+([0-9.]+)x", tail)
    return float(hits[-1]) if hits else None


# ------------------------------------------------------------------ evidence

def ev(layer, field, value, source, note=""):
    """One provenance record.  `source` is always a real path we opened."""
    return OrderedDict(layer=layer, field=field, value=value,
                       source=_rel(source), note=note)


def _rel(p):
    if not p:
        return p
    try:
        return os.path.relpath(p, ROOT).replace("\\", "/")
    except Exception:
        return str(p)


# -------------------------------------------------------------- XML helpers

def _ln(tag):
    """local-name: SLTs come both with and without a default namespace."""
    return tag.rsplit("}", 1)[-1]


def _find(el, name):
    for c in el.iter():
        if _ln(c.tag) == name:
            return c
    return None


def _iter(el, name):
    for c in el.iter():
        if _ln(c.tag) == name:
            yield c


# --------------------------------------------------------------- ISO-BMFF

def boxes(buf, start, end):
    off = start
    while off + 8 <= end:
        if off + 8 > len(buf):
            return
        size = struct.unpack(">I", buf[off:off + 4])[0]
        typ = buf[off + 4:off + 8].decode("latin1", "replace")
        hdr = 8
        if size == 1:
            if off + 16 > end:
                return
            size = struct.unpack(">Q", buf[off + 8:off + 16])[0]
            hdr = 16
        elif size == 0:
            size = end - off
        if size < hdr or off + size > end:
            return
        yield typ, off + hdr, off + size
        off += size


SAMPLE_ENTRY_VISUAL = {"encv", "hvc1", "hev1", "avc1"}
SAMPLE_ENTRY_AUDIO = {"enca", "mp4a", "ac-4", "ac-3", "ec-3"}
SAMPLE_ENTRY_TEXT = {"stpp", "wvtt", "tx3g"}
# What a track IS, taken from its sample entry rather than its hdlr: an MMT
# MPU init carries a second `hint` trak, and reading the last hdlr in the file
# labelled the video track "hint".  The sample entry cannot lie about this.
KIND_BY_ENTRY = ({e: "video" for e in SAMPLE_ENTRY_VISUAL} |
                 {e: "audio" for e in SAMPLE_ENTRY_AUDIO} |
                 {e: "text" for e in SAMPLE_ENTRY_TEXT})
TEXT_CONTENT_TYPES = {"text", "subtitles", "subtitle", "captions", "caption"}


def parse_init(path_or_bytes):
    """Read a fragmented-MP4 init segment: sample entry, sinf/schm, tenc.

    This is protection layer 3 and it is INDEPENDENT of the manifest: it is
    written by the packager, the MPD is written by the encoder toolchain.
    Two tools agreeing on the same default_KID is a real cross-check (E67).
    """
    if isinstance(path_or_bytes, (bytes, bytearray)):
        buf = bytes(path_or_bytes)
        src = "<bytes>"
    else:
        src = path_or_bytes
        with open(path_or_bytes, "rb") as fh:
            buf = fh.read(1 << 20)
    traks = []
    cur = {}

    def blank():
        return {"sample_entry": None, "handler": None, "width": None,
                "height": None, "schm": None, "frma": None, "tenc": None}

    def rec(s, e, depth=0):
        nonlocal cur
        for typ, ps, pe in boxes(buf, s, e):
            if depth > 10:
                return
            if typ == "trak":
                cur = blank()
                traks.append(cur)
                rec(ps, pe, depth + 1)
            elif typ in ("moov", "mdia", "minf", "stbl", "sinf", "schi",
                         "mvex", "edts", "udta", "dinf"):
                rec(ps, pe, depth + 1)
            elif typ == "hdlr":
                if cur:
                    cur["handler"] = buf[ps + 8:ps + 12].decode(
                        "latin1", "replace")
            elif typ == "stsd":
                rec(ps + 8, pe, depth + 1)      # skip version/flags + count
            elif typ in KIND_BY_ENTRY:
                if not cur:
                    cur = blank()
                    traks.append(cur)
                if cur["sample_entry"] is None:
                    cur["sample_entry"] = typ
                if typ in SAMPLE_ENTRY_VISUAL and pe - ps >= 78:
                    w, h = struct.unpack(">HH", buf[ps + 24:ps + 28])
                    if w and h:
                        cur["width"], cur["height"] = w, h
                    rec(ps + 78, pe, depth + 1)
                elif typ in SAMPLE_ENTRY_AUDIO and pe - ps >= 28:
                    rec(ps + 28, pe, depth + 1)
            elif typ == "schm" and cur:
                cur["schm"] = buf[ps + 4:ps + 8].decode("latin1", "replace")
            elif typ == "frma" and cur:
                cur["frma"] = buf[ps:ps + 4].decode("latin1", "replace")
            elif typ == "tenc" and cur:
                cur["tenc"] = {
                    "default_isProtected": buf[ps + 6],
                    "default_Per_Sample_IV_Size": buf[ps + 7],
                    "default_KID": buf[ps + 8:ps + 24].hex(),
                }
    rec(0, len(buf))
    real = [t for t in traks if t["sample_entry"]]
    out = dict(real[0]) if real else blank()
    out["source"] = _rel(src)
    out["found"] = bool(real)
    out["kind"] = KIND_BY_ENTRY.get(out.get("sample_entry") or "")
    out["traks"] = len(traks)
    return out


def parse_senc(seg, iv_size=16):
    """Protection layer 4: the bytes.  Delegates to E67's own implementation
    so the inspector and the gate cannot drift apart."""
    import e67_picture_gate
    return e67_picture_gate.parse_senc(seg, iv_size=iv_size)


# ------------------------------------------------------------------- TTML

def ttml_cues(seg_bytes):
    """Pull (begin, end, text) out of an stpp media segment's mdat.

    stpp rides UNENCRYPTED even on Widevine-protected services (E67 layer 3:
    tsi40 is plain `stpp`, never `encv`).  Captions of a locked channel are
    legitimately displayable, and this is where they come from.
    """
    i = seg_bytes.find(b"mdat")
    if i < 0:
        return []
    xml = seg_bytes[i + 4:]
    j = xml.find(b"<")
    if j < 0:
        return []
    txt = xml[j:].decode("utf-8", "replace")
    k = txt.rfind("</tt>")
    if k > 0:
        txt = txt[:k + 5]
    try:
        root = ET.fromstring(txt)
    except ET.ParseError:
        return []
    cues = []
    for p in _iter(root, "p"):
        body = re.sub(r"\s+", " ", _p_text(p)).strip()
        if not body:
            continue
        cues.append({"begin": p.get("begin"), "end": p.get("end"),
                     "text": body})
    # Roll-up captions re-send the whole two-line render on every character,
    # so a segment holds ~24 progressively longer copies of the same thing.
    # The LAST one in a segment is the finished render; the rest are typing.
    return cues[-1:] if cues else []


def _p_text(p):
    """Flatten a TTML <p>, turning <br/> into a space.

    itertext() concatenates the two caption rows with nothing between them,
    which welds the last word of row 1 to the first of row 2 ("KEEP ANEYE").
    """
    parts = []

    def walk(el, top=False):
        if _ln(el.tag) == "br":
            parts.append(" ")
        if el.text:
            parts.append(el.text)
        for c in el:
            walk(c)
        if not top and el.tail:
            parts.append(el.tail)
    walk(p, top=True)
    return "".join(parts)


def _ttml_seconds(v):
    if not v:
        return None
    m = re.match(r"^([0-9.]+)s$", v)
    if m:
        return float(m.group(1))
    m = re.match(r"^(\d+):(\d+):(\d+(?:\.\d+)?)$", v)
    if m:
        return int(m.group(1)) * 3600 + int(m.group(2)) * 60 + float(m.group(3))
    return None


def dedupe_cues(cues):
    """Roll-up captions repeat the growing line; keep only real advances."""
    out = []
    for c in cues:
        t = c["text"]
        if out and (t.startswith(out[-1]["text"]) or out[-1]["text"].endswith(t)):
            if len(t) >= len(out[-1]["text"]):
                out[-1] = dict(c)
            continue
        out.append(dict(c))
    return out


def srt_tail(path, n=8):
    """Last n cues of a live .srt written by the running chain (read-only)."""
    try:
        with open(path, "rb") as fh:
            fh.seek(0, io.SEEK_END)
            fh.seek(max(0, fh.tell() - 16384))
            blob = fh.read().decode("utf-8", "replace")
    except Exception:
        return []
    # Split on blank lines FIRST.  A lazy `(.+?)(?:\n\n|\Z)` with re.S ran
    # straight through the separator on CRLF text and glued the next cue's
    # index and timecode onto the caption line.
    blob = blob.replace("\r\n", "\n").replace("\r", "\n")
    cues = []
    for block in re.split(r"\n\s*\n", blob):
        m = re.search(r"(\d{2}:\d{2}:\d{2},\d{3})\s*-->\s*"
                      r"(\d{2}:\d{2}:\d{2},\d{3})\s*\n(.+)", block, re.S)
        if not m:
            continue
        cues.append({"begin": m.group(1), "end": m.group(2),
                     "text": re.sub(r"\s+", " ", m.group(3)).strip()})
    return dedupe_cues(cues)[-n:]


# ------------------------------------------------------------------- SLT

def parse_slt(path):
    """A carrier and its services, straight out of the SLT.

    The protection fields live HERE and nowhere else in the SLT: @protected
    (A/331 defaults it to false when absent) and @drmSystemID.
    """
    try:
        root = ET.parse(path).getroot()
    except Exception:
        return None
    if _ln(root.tag) != "SLT":
        return None
    carrier = {"bsid": root.get("bsid"), "source": _rel(path), "services": []}
    for s in _iter(root, "Service"):
        sig = None
        for c in s:
            if _ln(c.tag) == "BroadcastSvcSignaling":
                sig = c
        maj, mnr = s.get("majorChannelNo"), s.get("minorChannelNo")
        prot_attr = s.get("protected")
        svc = {
            "serviceId": s.get("serviceId"),
            "name": s.get("shortServiceName"),
            "vchannel": ("%s.%s" % (maj, mnr)) if maj and mnr else None,
            "serviceCategory": s.get("serviceCategory"),
            "globalServiceID": s.get("globalServiceID"),
            "configuration": s.get("configuration"),
            "broadbandAccessRequired": s.get("broadbandAccessRequired"),
            "slsProtocol": sig.get("slsProtocol") if sig is not None else None,
            "ip": sig.get("slsDestinationIpAddress") if sig is not None else None,
            "port": sig.get("slsDestinationUdpPort") if sig is not None else None,
            "slt_source": _rel(path),
            "_slt_protected_attr": prot_attr,
            "_slt_drmSystemID": s.get("drmSystemID"),
        }
        carrier["services"].append(svc)
    return carrier


SLS_PROTOCOL = {"1": "ROUTE/DASH", "2": "MMTP/MPU"}
SERVICE_CATEGORY = {"1": "linear A/V", "2": "linear audio",
                    "3": "app-based", "4": "ESG / data"}


# ------------------------------------------------------------------- SLS

def split_multipart(blob):
    txt = blob.decode("utf-8", "replace")
    m = re.search(r"boundary=\"([^\"]+)\"", txt) or re.search(r"--(\S+)", txt)
    if not m:
        return [("", txt)]
    bound = m.group(1).rstrip("-")
    parts = []
    for chunk in txt.split("--" + bound):
        chunk = chunk.strip("\r\n-")
        if not chunk.strip():
            continue
        if "\r\n\r\n" in chunk:
            hdr, body = chunk.split("\r\n\r\n", 1)
        elif "\n\n" in chunk:
            hdr, body = chunk.split("\n\n", 1)
        else:
            hdr, body = "", chunk
        parts.append((hdr.strip(), body))
    return parts


def parse_sls(path):
    """The ROUTE SLS multipart: MPD (tracks + ContentProtection) and S-TSID
    (which tsi carries which representation)."""
    with open(path, "rb") as fh:
        blob = fh.read()
    out = {"source": _rel(path), "mpd": None, "tracks": [],
           "content_protection": [], "tsi_map": {}, "mpd_kids": []}
    for hdr, body in split_multipart(blob):
        loc = ""
        m = re.search(r"Content-Location:\s*(\S+)", hdr, re.I)
        if m:
            loc = m.group(1)
        if "MPD" in body[:400] or loc.startswith("mpd"):
            out["mpd"] = loc or "mpd"
            _parse_mpd(body, out)
        elif "S-TSID" in body[:400] or loc.startswith("stsid"):
            _parse_stsid(body, out)
    return out


def _parse_mpd(body, out):
    try:
        root = ET.fromstring(body)
    except ET.ParseError:
        return
    for a in _iter(root, "AdaptationSet"):
        reps = [r for r in a if _ln(r.tag) == "Representation"]
        bw = None
        rep_id = None
        for r in reps:
            if r.get("bandwidth"):
                bw = int(r.get("bandwidth"))
            rep_id = r.get("id") or rep_id
        ctype = a.get("contentType") or (a.get("mimeType") or "").split("/")[0]
        t = {
            "contentType": ctype,
            "codecs": a.get("codecs"),
            "lang": a.get("lang"),
            "width": int(a.get("width")) if a.get("width") else None,
            "height": int(a.get("height")) if a.get("height") else None,
            "frameRate": a.get("frameRate"),
            "bandwidth_bps": bw,
            "repId": rep_id,
            "roles": [r.get("value") for r in a
                      if _ln(r.tag) == "Role" and r.get("value")],
        }
        for acc in a:
            if _ln(acc.tag) == "Accessibility" and acc.get("value"):
                t.setdefault("accessibility", []).append(acc.get("value"))
        out["tracks"].append(t)
    for cp in _iter(root, "ContentProtection"):
        rec = {"schemeIdUri": cp.get("schemeIdUri"),
               "value": cp.get("value")}
        for k, v in cp.attrib.items():
            if _ln(k) == "default_KID":
                rec["default_KID"] = v
                out["mpd_kids"].append(v)
        out["content_protection"].append(rec)


def _parse_stsid(body, out):
    try:
        root = ET.fromstring(body)
    except ET.ParseError:
        return
    for ls in _iter(root, "LS"):
        tsi = ls.get("tsi")
        mi = None
        for c in ls.iter():
            if _ln(c.tag) == "MediaInfo":
                mi = c
        if tsi:
            out["tsi_map"][tsi] = {
                "repId": mi.get("repId") if mi is not None else None,
                "contentType": mi.get("contentType") if mi is not None else None,
            }


# --------------------------------------------------- ROUTE object directory

OBJ_RE = re.compile(
    r"^obj_(\d+)_(\d+)_(\d+)_(\d+)_(\d+)_tsi(\d+)_toi(\d+)(_FDT)?\.(bin|xml)$")
INIT_TOI = 4294967295


def scan_route_dir(d):
    """-> {"ip:port": {"tsi": {"init": path, "media": [paths], "sls": path}}}"""
    flows = {}
    try:
        names = sorted(os.listdir(d))
    except OSError:
        return flows
    for n in names:
        m = OBJ_RE.match(n)
        if not m:
            continue
        a, b, c, e, port, tsi, toi, fdt, _ext = m.groups()
        key = "%s.%s.%s.%s:%s" % (a, b, c, e, port)
        p = os.path.join(d, n)
        f = flows.setdefault(key, {"dir": d, "tsi": {}})
        t = f["tsi"].setdefault(tsi, {"init": None, "media": [], "fdt": []})
        if fdt:
            t["fdt"].append(p)
        elif int(toi) == INIT_TOI:
            t["init"] = p
        elif tsi == "0":
            t.setdefault("sls", []).append(p)
        else:
            t["media"].append(p)
    return flows


# -------------------------------------------------------------------- ESG

def load_esg(dirs):
    """Decode the broadcaster's service guide out of banked SGDU objects.

    Reuses lab/m41_sgdu.py's container reverse-engineering verbatim -- that
    decoder was gated when it was built and this file does not re-derive it.
    """
    import m41_sgdu as m41
    OMA = m41.OMA
    SA = m41.SA
    services, contents, schedules = {}, {}, []
    units = 0
    srcs = []
    for d in dirs:
        if not os.path.isdir(d):
            continue
        for f in sorted(glob.glob(os.path.join(d, "*"))):
            if f.endswith((".json", ".srt", ".mp4", ".png", ".ts", ".wav",
                           ".log", ".md")):
                continue
            if not os.path.isfile(f):
                continue
            try:
                raw = open(f, "rb").read()
            except OSError:
                continue
            if len(raw) > 4 << 20:
                continue
            try:
                parts = m41.parse_sgdu(m41.unwrap(raw))
            except Exception:
                continue
            if not parts:
                continue
            units += 1
            if _rel(d) not in srcs:
                srcs.append(_rel(d))
            for _fid, x in parts:
                try:
                    root = ET.fromstring(x.decode("utf-8", "replace"))
                except ET.ParseError:
                    continue
                tag = _ln(root.tag)
                if tag == "Service":
                    nm = root.find(OMA + "Name")
                    maj = root.find(".//" + SA + "MajorChannelNum")
                    mnr = root.find(".//" + SA + "MinorChannelNum")
                    services[root.get("id")] = {
                        "name": nm.get("text") if nm is not None else "?",
                        "globalServiceID": root.get("globalServiceID"),
                        "vchannel": ("%s.%s" % (maj.text, mnr.text)
                                     if maj is not None and mnr is not None
                                     else None),
                    }
                elif tag == "Content":
                    nm = root.find(OMA + "Name")
                    de = root.find(OMA + "Description")
                    rat = [r.text.strip()
                           for r in root.iter(SA + "RatingDescription")
                           if r.text and r.text.strip()]
                    contents[root.get("id")] = {
                        "name": nm.get("text") if nm is not None else "?",
                        "desc": de.get("text") if de is not None else "",
                        "rating": "; ".join(rat),
                    }
                elif tag == "Schedule":
                    svc = root.find(OMA + "ServiceReference")
                    for cr in root.iter(OMA + "ContentReference"):
                        pw = cr.find(OMA + "PresentationWindow")
                        if pw is None:
                            continue
                        schedules.append({
                            "svc": svc.get("idRef") if svc is not None else "?",
                            "cid": cr.get("idRef"),
                            "start": int(pw.get("startTime", 0)) - NTP_EPOCH,
                            "end": int(pw.get("endTime", 0)) - NTP_EPOCH,
                        })
    rows = []
    for s in schedules:
        c = contents.get(s["cid"])
        if not c:
            continue
        svc = services.get(s["svc"], {})
        rows.append({
            "gsid": svc.get("globalServiceID"),
            "svc_name": svc.get("name"),
            "vchannel": svc.get("vchannel"),
            "start": s["start"], "end": s["end"],
            "name": c["name"], "desc": c["desc"], "rating": c["rating"],
        })
    rows.sort(key=lambda r: (r["gsid"] or "", r["start"]))
    span = None
    if rows:
        span = (min(r["start"] for r in rows), max(r["end"] for r in rows))
    return {"rows": rows, "units": units, "sources": srcs,
            "services": len({r["gsid"] for r in rows}), "span": span}


def now_next(esg, gsid, name, at):
    """Current and next programme for one service, at unix time `at`.

    Returns (now, next, status).  `status` is 'live' when `at` lies inside the
    guide's span, 'stale' when the banked guide does not reach the present --
    a guide from 8/07 cannot tell you what is on now, and pretending it can is
    the kind of thing this campaign logs as a bug.
    """
    rows = [r for r in esg["rows"]
            if (gsid and r["gsid"] == gsid) or
               (not gsid and name and r["svc_name"] == name)]
    if not rows:
        return None, None, "no-guide"
    rows.sort(key=lambda r: r["start"])
    lo, hi = rows[0]["start"], rows[-1]["end"]
    status = "live" if lo <= at <= hi else "stale"
    if status == "stale":
        return None, None, "stale"
    cur = nxt = None
    for r in rows:
        if r["start"] <= at < r["end"]:
            cur = r
        elif r["start"] > at and nxt is None:
            nxt = r
    return cur, nxt, status


# ------------------------------------------------------------ PLP evidence
#
# The PLP a service rides is NOT carried in any signalling document we can
# read (it lives in L1-Detail, which our object extractors do not preserve).
# So it is attributed the only honest way: a flow is on PLP n because it was
# recovered from a decode of PLP n, and the decode artifact is named here.
# Anything not in this table reports "unknown" -- never a guess.

PLP_ATTRIBUTION = {
    # flow                  plp   artifact the flow was recovered from
    "239.255.32.1:8321":   ("1", "lab/m13_plp1_48.bb -> lab/m13_out48/"),
    "239.255.5.1:8051":    ("1", "lab/m13_plp1_48.bb -> lab/m13_out48/"),
    "239.255.4.1:8041":    ("1", "lab/m13_plp1_48.bb -> lab/m13_out48/"),
    "239.255.9.1:8091":    ("1", "lab/m13_plp1_48.bb -> lab/m13_out48/"),
    "239.255.7.1:8071":    ("0", "lab/m7_long.bb (PLP0) -> lab/m7_out/"),
    "239.255.0.255:8000":  ("0", "lab/m7_long.bb (PLP0) -> lab/m7_out/"),
}
PLP_BASIS = ("flow recovered from the named PLP's baseband; the PLP is "
             "recorded in the decode artifact, not in a signalling field")


# ------------------------------------------------------- protection verdict

WIDEVINE_UUID = "edef8ba9-79d6-4ace-a3c8-27dcd51d21ed"
DRM_NAMES = {
    WIDEVINE_UUID: "Widevine",
    "9a04f079-9840-4286-ab92-e65be0885f95": "PlayReady",
    "94ce86fb-07ff-4f43-adb8-93d2fa968ca2": "FairPlay",
}


def drm_name(uri):
    if not uri:
        return None
    u = uri.lower().replace("urn:uuid:", "").strip()
    return DRM_NAMES.get(u, "unrecognised DRM system")


def classify_protection(svc):
    """The whole point of this tool.  Read four layers, cite every one.

    Returns {"state": CLEAR|PROTECTED|UNKNOWN, "drm": name|None,
             "evidence": [...], "layers_read": n}
    """
    evs = []
    positive = False
    layers_read = 0
    drm = None

    # ---- L1: the SLT ----------------------------------------------------
    if svc.get("slt_source"):
        layers_read += 1
        pa = svc.get("_slt_protected_attr")
        did = svc.get("_slt_drmSystemID")
        if pa is None:
            evs.append(ev("SLT", "Service/@protected", "(absent)",
                          svc["slt_source"],
                          "A/331 defaults @protected to false when absent"))
        else:
            evs.append(ev("SLT", "Service/@protected", pa, svc["slt_source"]))
            if str(pa).lower() == "true":
                positive = True
        if did:
            drm = drm or drm_name(did)
            evs.append(ev("SLT", "Service/@drmSystemID", did,
                          svc["slt_source"], drm or ""))
            positive = True
        else:
            evs.append(ev("SLT", "Service/@drmSystemID", "(absent)",
                          svc["slt_source"]))

    # ---- L2: the MPD inside the SLS -------------------------------------
    sls = svc.get("_sls")
    if sls:
        layers_read += 1
        cps = sls.get("content_protection") or []
        if not cps:
            evs.append(ev("MPD", "ContentProtection", "(no element)",
                          sls["source"],
                          "an unprotected DASH manifest carries none"))
        for cp in cps:
            positive = True
            evs.append(ev("MPD", "ContentProtection/@schemeIdUri",
                          cp.get("schemeIdUri"), sls["source"],
                          cp.get("value") or ""))
            if cp.get("default_KID"):
                drm = drm or drm_name(cp.get("schemeIdUri"))
                evs.append(ev("MPD", "ContentProtection/cenc:default_KID",
                              cp["default_KID"], sls["source"],
                              "key IDENTIFIER, not key material"))

    # ---- L3: the init segments ------------------------------------------
    inits = svc.get("_inits") or {}
    if inits:
        layers_read += 1
    for tsi, info in sorted(inits.items()):
        se = info.get("sample_entry")
        if se is None:
            continue
        # The note must describe THIS entry, not the general rule -- it read
        # "encv/enca = an encrypted sample entry" beside tsi40's plain `stpp`,
        # which says the opposite of the truth about the caption track.
        enc = se in ("encv", "enca")
        evs.append(ev("init", "tsi%s stsd sample entry" % tsi, se,
                      info["source"],
                      "an ENCRYPTED sample entry" if enc
                      else "a plain sample entry -- this track is not encrypted"))
        if enc:
            positive = True
        if info.get("schm"):
            evs.append(ev("init", "tsi%s sinf/schm@scheme_type" % tsi,
                          info["schm"], info["source"]))
            if info["schm"] in ("cenc", "cbcs", "cens", "cbc1"):
                positive = True
        if info.get("tenc"):
            t = info["tenc"]
            evs.append(ev("init", "tsi%s tenc.default_isProtected" % tsi,
                          t["default_isProtected"], info["source"]))
            evs.append(ev("init", "tsi%s tenc.default_KID" % tsi,
                          t["default_KID"], info["source"],
                          "key IDENTIFIER, not key material"))
            if t["default_isProtected"]:
                positive = True

    # ---- L4: the bytes ---------------------------------------------------
    senc = svc.get("_senc")
    if senc is not None:
        layers_read += 1
        if senc.get("absent"):
            evs.append(ev("bytes", "moof/traf/senc", "(absent)",
                          senc["source"],
                          "no per-sample encryption boxes in the media"))
        else:
            evs.append(ev("bytes", "senc protected byte fraction",
                          "%.1f%%" % (100 * (senc.get("protected_frac") or 0)),
                          senc["source"],
                          "%d/%d IVs unique/non-zero" %
                          (senc.get("ivs_unique", 0),
                           senc.get("ivs_nonzero", 0))))
            if (senc.get("protected_frac") or 0) > 0:
                positive = True

    if layers_read == 0:
        state = "UNKNOWN"
    elif positive:
        state = "PROTECTED"
    else:
        state = "CLEAR"
    return {"state": state, "drm": drm, "evidence": evs,
            "layers_read": layers_read}


# ------------------------------------------------------------- playability

def decide_playable(svc):
    """THE ONLY place `playable` is decided.  Protection short-circuits.

    A protected service is LOCKED.  Not "maybe", not "if we had a key", not
    "playable=True but greyed out in the UI" -- False, with a reason string
    the UI is required to show.  lab/e70_gate.py asserts that no code path
    can produce playable=True while protection.state == 'PROTECTED'.
    """
    prot = svc["protection"]["state"]
    if prot == "PROTECTED":
        return {
            "playable": False,
            "reason": "LOCKED",
            "detail": ("encrypted by the broadcaster (%s); guide, manifest "
                       "and captions are in the clear and shown, the media "
                       "is not decodable and is never presented as watchable"
                       % (svc["protection"]["drm"] or "DRM")),
        }
    if prot == "UNKNOWN":
        return {"playable": False, "reason": "UNKNOWN",
                "detail": "no protection signalling was readable for this "
                          "service; refusing to claim it is watchable"}
    if svc.get("broadbandAccessRequired") == "true":
        return {"playable": False, "reason": "BROADBAND",
                "detail": "signalled broadbandAccessRequired=true: the media "
                          "does not come off the air at all"}
    if svc.get("serviceCategory") == "4":
        return {"playable": False, "reason": "DATA",
                "detail": "service category 4 (ESG / data): no A/V component"}
    if not svc.get("_has_media"):
        return {"playable": False, "reason": "NO-MEDIA",
                "detail": "clear, but no media objects are banked for it here"}
    out = {"playable": True, "reason": "CLEAR",
           "detail": "signalled in the clear at every layer read, and media "
                     "objects are present"}
    # Playable is not the same as playable-in-real-time.  Say which.
    plp = (svc.get("plp") or {}).get("plp")
    if plp not in (None, "0", "16"):
        out["realtime"] = False
        out["realtime_note"] = (
            "PLP%s: our live chain decodes PLP0+PLP16 only, and E47 measured "
            "the working PLP%s decoder ~140x short of real time -- clear and "
            "decodable, but OFFLINE ONLY today" % (plp, plp))
    else:
        out["realtime"] = True
    return out


def picture_claim(_svc):
    """E67's law, encoded: this module never claims a picture.

    Frame count, exit code, error count, container parse and NAL tiling all
    pass on CENC-encrypted video.  The only thing that discriminated was
    flat_frac.  So no caller gets a thumbnail or a "has picture" flag from
    here -- if you want one, run lab/e67_picture_gate.py and carry ITS
    verdict.  Returns the refusal, deliberately.
    """
    return {"claims_picture": False,
            "gate": "lab/e67_picture_gate.py",
            "why": "'the decoder ran' is not 'the decoder decoded' (E67)"}


# -------------------------------------------------------------- live dirs

def scan_live_dir(d):
    """A running chain's output directory: what it is tuned to, its lanes,
    and its captions.  Read-only; nothing here writes or signals."""
    out = {"dir": _rel(d), "rf": None, "transport": None, "lanes": {},
           "captions": [], "chain": {}}
    log = os.path.join(d, "chain.log")
    if os.path.isfile(log):
        try:
            first = open(log, "r", errors="replace").readline()
        except OSError:
            first = ""
        m = re.search(r"\bRF(\d+)\b", first)
        if m:
            out["rf"] = int(m.group(1))
            out["rf_source"] = _rel(log) + ":1"
        blob = ""
        try:
            with open(log, "rb") as fh:
                blob = fh.read(60000).decode("utf-8", "replace")
        except OSError:
            pass
        if re.search(r"->\s*MMTP", blob):
            out["transport"] = "MMTP/MPU"
            out["transport_source"] = _rel(log) + " flow gate"
        inst = soak_instantaneous(log)
        if inst is not None:
            out["chain"]["inst_x"] = inst
        mfec = re.findall(r"FEC \d+/\d+ \[([0-9.]+)% now\]",
                          _tail_text(log, 4096))
        if mfec:
            out["chain"]["fec_now_pct"] = float(mfec[-1])
    lj = os.path.join(d, "live.json")
    if os.path.isfile(lj):
        try:
            j = json.load(open(lj))
        except Exception:
            j = {}
        for name, lane in (j.get("lanes") or {}).items():
            p = lane.get("path")
            rec = {"kind": lane.get("kind"), "pid": lane.get("pid"),
                   "media_s": lane.get("media_s"),
                   "bytes": lane.get("bytes"), "path": _rel(p or "")}
            if p and os.path.isfile(p) and lane.get("init_bytes"):
                try:
                    with open(p, "rb") as fh:
                        head = fh.read(int(lane["init_bytes"]))
                    rec["init"] = parse_init(head)
                    rec["init"]["source"] = _rel(p) + " [first %d bytes = init]" \
                        % int(lane["init_bytes"])
                except OSError:
                    pass
            out["lanes"][name] = rec
    srt = os.path.join(d, "live.srt")
    if os.path.isfile(srt):
        out["captions"] = srt_tail(srt, 8)
        out["captions_source"] = _rel(srt)
    return out


def _tail_text(path, n):
    try:
        with open(path, "rb") as fh:
            fh.seek(0, io.SEEK_END)
            fh.seek(max(0, fh.tell() - n))
            return fh.read().decode("utf-8", "replace")
    except OSError:
        return ""


# ------------------------------------------------------------ the lineup

DEFAULT_SOURCES = [
    "lab/m7_out",              # RF33 PLP0: SLT, WJLA MMTP, the ESG
    "lab/stress_0807/live_out2",   # RF33 PLP0 again: a fresher ESG
    "lab/m13_out48",           # RF33 PLP1: WHUT + the three locked services
    "lab/plp1_out",            # RF33 PLP1: 132.1 media
    "lab/m10_rf30_out",        # RF30: a second carrier's SLT
    "data/e31",                # the live chain (read-only)
]


def build_lineup(sources=None, at=None):
    sources = sources or DEFAULT_SOURCES
    at = at or dt.datetime.now(dt.timezone.utc).timestamp()
    dirs = [s if os.path.isabs(s) else os.path.join(ROOT, s) for s in sources]
    dirs = [d for d in dirs if os.path.isdir(d)]

    carriers = OrderedDict()
    flows = {}
    live_dirs = []
    non_slt = []
    for d in dirs:
        # SLTs
        for f in sorted(glob.glob(os.path.join(d, "lls_*.xml"))):
            c = parse_slt(f)
            if not c:
                # An LLS table that is not an SLT. 8/15, RF30: the ONLY table
                # on that carrier was table_id 255 (reserved/private, a
                # datacast), and this tool printed a header and nothing else
                # -- indistinguishable from "no signalling decoded". Record
                # what WAS there so the report can say so.
                try:
                    tag = _ln(ET.parse(f).getroot().tag)
                except Exception:                                # noqa: BLE001
                    tag = "(unparseable)"
                non_slt.append((_rel(f), tag))
                continue
            key = c["bsid"]
            if key not in carriers:
                c["rf"] = carrier_rf(f)
                c["center_mhz"] = rf_center_mhz(
                    (c["rf"] or {}).get("rf"))
                carriers[key] = c
        # ROUTE objects
        for k, v in scan_route_dir(d).items():
            cur = flows.setdefault(k, {"tsi": {}, "dirs": []})
            cur["dirs"].append(_rel(d))
            for tsi, t in v["tsi"].items():
                ct = cur["tsi"].setdefault(
                    tsi, {"init": None, "media": [], "sls": [], "fdt": []})
                ct["init"] = ct["init"] or t.get("init")
                ct["media"].extend(t.get("media") or [])
                ct["sls"].extend(t.get("sls") or [])
                ct["fdt"].extend(t.get("fdt") or [])
        if os.path.isfile(os.path.join(d, "chain.log")):
            live_dirs.append(d)

    esg = load_esg(dirs)
    lives = [scan_live_dir(d) for d in live_dirs]

    # attach evidence to every service
    for bsid, c in carriers.items():
        for svc in c["services"]:
            key = "%s:%s" % (svc.get("ip"), svc.get("port"))
            svc["flow"] = key
            f = flows.get(key)
            svc["_has_media"] = False
            svc["_inits"] = {}
            svc["_senc"] = None
            svc["_sls"] = None
            svc["captions"] = []
            svc["captions_source"] = None

            if f:
                sls_paths = sorted(f["tsi"].get("0", {}).get("sls", []))
                if sls_paths:
                    try:
                        svc["_sls"] = parse_sls(sls_paths[-1])
                    except Exception:
                        svc["_sls"] = None
                for tsi, t in sorted(f["tsi"].items()):
                    if tsi == "0":
                        continue
                    if t.get("init"):
                        try:
                            svc["_inits"][tsi] = parse_init(t["init"])
                        except Exception:
                            pass
                    if t.get("media"):
                        svc["_has_media"] = True
                # layer 4 + captions: pick the video tsi and the text tsi
                tsi_map = (svc["_sls"] or {}).get("tsi_map", {})
                for tsi, t in sorted(f["tsi"].items()):
                    ctype = (tsi_map.get(tsi) or {}).get("contentType")
                    if ctype in TEXT_CONTENT_TYPES:
                        ctype = "text"
                    if ctype is None and tsi in svc["_inits"]:
                        ctype = svc["_inits"][tsi].get("kind")
                    media = sorted(t.get("media") or [],
                                   key=lambda p: int(re.search(r"toi(\d+)",
                                                               p).group(1)))
                    if ctype == "video" and media and svc["_senc"] is None:
                        # E67's discipline: the FIRST media object of a flow is
                        # a partial -- we joined the flow mid-object -- so it
                        # has no moof/senc and reading it reports "(absent)" on
                        # a service that is in fact 98.3% encrypted.  Skip it,
                        # and keep trying until a whole object answers.
                        ivs = 16
                        ti = svc["_inits"].get(tsi, {}).get("tenc")
                        if ti:
                            ivs = ti["default_Per_Sample_IV_Size"] or 16
                        cands = media[1:] or media
                        last = None
                        for p in cands[:4]:
                            s = parse_senc(open(p, "rb").read(), iv_size=ivs)
                            last = p
                            if s is not None:
                                s = dict(s)
                                s["absent"] = False
                                s["source"] = _rel(p)
                                svc["_senc"] = s
                                break
                        if svc["_senc"] is None and last:
                            svc["_senc"] = {
                                "absent": True, "source": _rel(last),
                                "tried": len(cands[:4]),
                                "note": "first (partial) object skipped"}
                    if ctype == "text" and media:
                        cues = []
                        for p in media[:6]:
                            cues.extend(ttml_cues(open(p, "rb").read()))
                        cues = dedupe_cues(cues)
                        svc["captions_source"] = _rel(os.path.dirname(media[0]))
                        svc["captions_kind"] = "stpp/TTML (unencrypted)"
                        if cues:
                            svc["captions"] = cues[-8:]
                            svc["captions_status"] = "cues"
                        else:
                            # not a parse failure: the broadcaster sent a
                            # well-formed TTML document with an empty <div>.
                            # Say so rather than showing a blank cell.
                            svc["captions_status"] = "empty"

            svc["protection"] = classify_protection(svc)
            plp = PLP_ATTRIBUTION.get(key)
            svc["plp"] = {"plp": plp[0] if plp else None,
                          "source": plp[1] if plp else None,
                          "basis": PLP_BASIS if plp else "not attributed"}
            svc["transport"] = SLS_PROTOCOL.get(svc.get("slsProtocol"),
                                                svc.get("slsProtocol"))
            svc["transport_source"] = (
                _rel(svc["slt_source"]) + " BroadcastSvcSignaling/@slsProtocol")
            svc["category"] = SERVICE_CATEGORY.get(svc.get("serviceCategory"),
                                                   svc.get("serviceCategory"))
            svc["media"] = _media_summary(svc)
            svc["playable"] = decide_playable(svc)
            svc["picture"] = picture_claim(svc)
            cur, nxt, st = now_next(esg, svc.get("globalServiceID"),
                                    svc.get("name"), at)
            svc["now"] = cur
            svc["next"] = nxt
            svc["guide_status"] = st

    # bind live directories to the service they are actually carrying
    for lv in lives:
        lv["service"] = _bind_live(lv, carriers)
        if lv.get("service") and lv.get("captions"):
            for c in carriers.values():
                for svc in c["services"]:
                    if svc["serviceId"] == lv["service"]["serviceId"] and \
                            c["bsid"] == lv["service"]["bsid"]:
                        svc["captions"] = lv["captions"]
                        svc["captions_source"] = lv.get("captions_source")
                        svc["captions_kind"] = "live chain .srt (from stpp)"
                        svc["live"] = {
                            "dir": lv["dir"], "inst_x": lv["chain"].get("inst_x"),
                            "fec_now_pct": lv["chain"].get("fec_now_pct"),
                            "lanes": {k: {"kind": v["kind"], "pid": v["pid"],
                                          "media_s": v["media_s"]}
                                      for k, v in lv["lanes"].items()},
                        }
                        # the live init is a fourth independent protection read
                        for k, lane in lv["lanes"].items():
                            if lane.get("init"):
                                svc["_inits"]["live_pid%s" % lane["pid"]] = \
                                    lane["init"]
                        svc["protection"] = classify_protection(svc)
                        svc["media"] = _media_summary(svc)
                        svc["playable"] = decide_playable(svc)

    return {"build": BUILD, "at": at,
            "sources": [_rel(d) for d in dirs],
            "carriers": carriers, "esg": esg, "live": lives,
            "non_slt_lls": non_slt}


# RF channel -> centre frequency, for display only (A/321 US UHF/VHF plan).
def rf_center_mhz(rf):
    if rf is None:
        return None
    if 2 <= rf <= 4:
        return 57 + 6 * (rf - 2)
    if 5 <= rf <= 6:
        return 79 + 6 * (rf - 5)
    if 7 <= rf <= 13:
        return 177 + 6 * (rf - 7)
    if 14 <= rf <= 36:
        return 473 + 6 * (rf - 14)
    return None


def carrier_rf(slt_path):
    """Which RF channel is this carrier on?  Derived, with the chain shown.

    No signalling document we bank names the RF channel -- the SLT does not
    carry it.  But the SLT was decoded out of a named capture, and the capture
    filename records the channel.  So: SLT -> the decode summary sitting beside
    its object directory -> "capture": "long_rf33.cs16" -> RF33.  Returned with
    the summary path so the reader can check the link themselves.
    """
    d = os.path.dirname(slt_path)
    base = os.path.basename(d)
    stem = base[:-4] if base.endswith("_out") else base
    stem = re.sub(r"_(long|test|48)$", "", stem)
    for cand in sorted(glob.glob(os.path.join(os.path.dirname(d),
                                              stem + "*summary*.json"))):
        try:
            j = json.load(open(cand))
        except Exception:
            continue
        cap = j.get("capture") or ""
        m = re.search(r"rf(\d+)", cap, re.I)
        if m:
            return {"rf": int(m.group(1)), "capture": cap,
                    "source": _rel(cand),
                    "basis": "the SLT was decoded from this capture; the RF "
                             "channel is recorded in the capture name"}
    return None


def _bind_live(lv, carriers):
    """Which service is this live directory carrying?

    Stated as a derivation, not asserted, in three steps that are each
    checkable: (1) the chain log's own first line names the RF channel; (2)
    each banked carrier is tied to an RF channel through the capture its SLT
    was decoded from; (3) if exactly one broadcast service on THAT carrier
    uses the transport the flow gate locked, the binding is forced.  More than
    one candidate and we refuse rather than pick.
    """
    if lv.get("rf") is None or lv.get("transport") is None:
        return None
    cands = []
    for bsid, c in carriers.items():
        rf = c.get("rf")
        if not rf or rf.get("rf") != lv["rf"]:
            continue
        for svc in c["services"]:
            if SLS_PROTOCOL.get(svc.get("slsProtocol")) != lv["transport"]:
                continue
            if svc.get("configuration") == "Broadband":
                continue
            cands.append((bsid, c, svc))
    if len(cands) != 1:
        return None
    bsid, c, svc = cands[0]
    return {"bsid": bsid, "serviceId": svc["serviceId"], "name": svc["name"],
            "vchannel": svc["vchannel"],
            "derivation": ("%s says RF%d; %s ties bsid %s to RF%d; the flow "
                           "gate says %s and exactly one broadcast service on "
                           "that carrier uses it" %
                           (lv.get("rf_source", "chain.log"), lv["rf"],
                            c["rf"]["source"], bsid, c["rf"]["rf"],
                            lv["transport"]))}


def _media_summary(svc):
    sls = svc.get("_sls")
    out = {"video": None, "audio": [], "text": [], "source": None}
    if sls and sls.get("tracks"):
        out["source"] = sls["source"] + " (MPD)"
        for t in sls["tracks"]:
            if t["contentType"] == "video" and out["video"] is None:
                out["video"] = {
                    "resolution": ("%dx%d" % (t["width"], t["height"]))
                    if t["width"] and t["height"] else None,
                    "codec": t["codecs"], "frameRate": t["frameRate"],
                    "bitrate_bps": t["bandwidth_bps"]}
            elif t["contentType"] == "audio":
                out["audio"].append({"lang": t["lang"], "codec": t["codecs"],
                                     "bitrate_bps": t["bandwidth_bps"],
                                     "roles": t.get("roles") or []})
            elif t["contentType"] in ("text", "application"):
                out["text"].append({"lang": t["lang"], "codec": t["codecs"],
                                    "roles": t.get("roles") or []})
        return out
    # no MPD (MMTP): read the init segments we do have
    for tsi, info in sorted((svc.get("_inits") or {}).items()):
        out["source"] = out["source"] or (info["source"] + " (init segment)")
        kind = info.get("kind")
        se = info.get("sample_entry")
        if kind == "video" and out["video"] is None:
            out["video"] = {
                "resolution": ("%dx%d" % (info["width"], info["height"]))
                if info.get("width") else None,
                "codec": info.get("frma") or se, "frameRate": None,
                "bitrate_bps": None}
        elif kind == "audio":
            out["audio"].append({"lang": None,
                                 "codec": info.get("frma") or se,
                                 "bitrate_bps": None, "roles": []})
        elif kind == "text":
            out["text"].append({"lang": None, "codec": se, "roles": []})
    return out


# ---------------------------------------------------------------- reporting

def _fmt_bps(b):
    if not b:
        return "-"
    return "%.2f Mbps" % (b / 1e6) if b >= 1e6 else "%d kbps" % (b / 1e3)


def _fmt_time(ts):
    if not ts:
        return "-"
    return dt.datetime.fromtimestamp(ts, dt.timezone.utc).strftime(
        "%m-%d %H:%MZ")


def report(model, only=None, verbose=False, out=sys.stdout):
    w = out.write
    at = model["at"]
    w("ATSC 3.0 LOCKED-CHANNEL INSPECTOR   build %s\n" % model["build"])
    w("reference time %s   sources: %s\n" %
      (_fmt_time(at), ", ".join(model["sources"])))
    e = model["esg"]
    if e["rows"]:
        w("service guide: %d programmes, %d services, %s .. %s  (%s)\n" %
          (len(e["rows"]), e["services"], _fmt_time(e["span"][0]),
           _fmt_time(e["span"][1]), ", ".join(e["sources"])))
    w("\n")
    if not model["carriers"]:
        ns = model.get("non_slt_lls") or []
        if ns:
            w("NO SLT (service list) in any source -- no television services "
              "are signalled here.\nLLS tables that WERE present:\n")
            for path, tag in ns:
                w("    %s   root element <%s>\n" % (path, tag))
            w("A carrier that signals only reserved/private LLS tables is a "
              "datacast (a data service riding\nATSC 3.0), not a television "
              "multiplex: nothing here is watchable, and nothing is "
              "encrypted either.\n")
        else:
            w("NO SLT and no other LLS table in any source -- nothing was "
              "decoded to inspect.\n")

    for bsid, c in model["carriers"].items():
        w("=" * 78 + "\n")
        rf = c.get("rf")
        w("CARRIER bsid %s   %s   (%s)\n" %
          (bsid,
           ("RF%d / %s MHz   [%s]" % (rf["rf"], c.get("center_mhz"),
                                      rf["source"])) if rf else "RF unknown",
           c["source"]))
        w("=" * 78 + "\n")
        for svc in c["services"]:
            if only and only.upper() not in (
                    (svc.get("name") or "").upper(),
                    (svc.get("vchannel") or ""),
                    str(svc.get("serviceId"))):
                continue
            _report_service(svc, w, verbose)
    return model


def _report_service(svc, w, verbose):
    p = svc["protection"]
    play = svc["playable"]
    badge = {"PROTECTED": "[ ENCRYPTED -- guide and captions only ]",
             "CLEAR": "[ clear ]",
             "UNKNOWN": "[ protection UNKNOWN ]"}[p["state"]]
    w("\n%-8s %-7s  sid %-5s %s\n" %
      (svc.get("name") or "?", svc.get("vchannel") or "-",
       svc.get("serviceId"), badge))
    w("  %-14s %s\n" % ("category", svc.get("category")))
    w("  %-14s %s   (%s)\n" % ("transport", svc.get("transport"),
                               svc.get("transport_source")))
    plp = svc["plp"]
    w("  %-14s %s   (%s)\n" % ("PLP", plp["plp"] or "unknown",
                               plp["source"] or plp["basis"]))
    w("  %-14s %s -- %s\n" % ("PLAYABLE",
                              "YES" if play["playable"] else "NO  (%s)" %
                              play["reason"], play["detail"]))
    if play.get("realtime") is False:
        w("  %-14s %s\n" % ("", play["realtime_note"]))
    w("  %-14s %s%s   [%d layer(s) read]\n" %
      ("protection", p["state"],
       ("  " + p["drm"]) if p["drm"] else "", p["layers_read"]))
    for e in p["evidence"] if verbose else _key_evidence(p["evidence"]):
        w("      %-6s %-38s = %-42s  %s\n" %
          (e["layer"], e["field"], str(e["value"])[:42], e["source"]))
    m = svc.get("media") or {}
    if m.get("video"):
        v = m["video"]
        w("  %-14s %s %s %s %s\n" % ("video", v.get("resolution") or "?",
                                     v.get("codec") or "?",
                                     ("@%s" % v["frameRate"]) if v.get("frameRate") else "",
                                     _fmt_bps(v.get("bitrate_bps"))))
    for a in m.get("audio") or []:
        w("  %-14s %-5s %-16s %-10s %s\n" %
          ("audio", a.get("lang") or "-", a.get("codec") or "?",
           _fmt_bps(a.get("bitrate_bps")), ",".join(a.get("roles") or [])))
    for t in m.get("text") or []:
        w("  %-14s %-5s %s %s\n" % ("captions", t.get("lang") or "-",
                                    t.get("codec") or "?",
                                    ",".join(t.get("roles") or [])))
    if m.get("source"):
        w("  %-14s %s\n" % ("", m["source"]))
    gs = svc.get("guide_status")
    if svc.get("now"):
        n = svc["now"]
        w("  %-14s %s  %s .. %s  %s\n" % ("NOW", n["name"],
                                          _fmt_time(n["start"]),
                                          _fmt_time(n["end"]), n["rating"]))
        if n.get("desc"):
            w("  %-14s %s\n" % ("", n["desc"][:70]))
    if svc.get("next"):
        n = svc["next"]
        w("  %-14s %s  %s\n" % ("NEXT", n["name"], _fmt_time(n["start"])))
    if gs == "stale":
        w("  %-14s banked guide does not reach the reference time "
          "(no now/next claimed)\n" % "guide")
    elif gs == "no-guide":
        w("  %-14s no guide fragments for this service\n" % "guide")
    if svc.get("live"):
        lv = svc["live"]
        w("  %-14s %s  inst %sx  FEC %s%%  lanes %s\n" %
          ("LIVE", lv["dir"], lv.get("inst_x"), lv.get("fec_now_pct"),
           ", ".join("%s(pid%s)" % (v["kind"], v["pid"])
                     for v in lv["lanes"].values())))
    caps = svc.get("captions") or []
    if caps:
        w("  %-14s %s   (%s)\n" % ("CAPTIONS", svc.get("captions_kind", ""),
                                   svc.get("captions_source")))
        for c in caps[-5:]:
            w("      %-14s %s\n" % (c.get("begin") or "", c["text"][:64]))
    elif svc.get("captions_status") == "empty":
        w("  %-14s stpp track present and unencrypted, but the banked "
          "segments carry an empty <div> -- no cues in this window (%s)\n" %
          ("CAPTIONS", svc.get("captions_source")))


KEY_FIELDS = ("Service/@protected", "Service/@drmSystemID", "ContentProtection",
              "senc protected byte fraction")


def _key_evidence(evs):
    keep, seen = [], set()
    for e in evs:
        f = e["field"]
        if not (any(f.startswith(k) or k in f for k in KEY_FIELDS)
                or "stsd sample entry" in f or "schm" in f):
            continue
        sig = (e["layer"], f, str(e["value"]), e["source"])
        if sig in seen:
            continue
        seen.add(sig)
        keep.append(e)
    return keep[:9]


def strip_private(o):
    """Drop the _-prefixed working fields before serialising."""
    if isinstance(o, dict):
        return {k: strip_private(v) for k, v in o.items()
                if not k.startswith("_")}
    if isinstance(o, list):
        return [strip_private(v) for v in o]
    return o


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--source", action="append", default=None,
                    help="a live dir or banked object dir (repeatable)")
    ap.add_argument("--service", default=None,
                    help="restrict to one service (name, vchannel or sid)")
    ap.add_argument("--at", type=float, default=None,
                    help="reference unix time for now/next (default: now)")
    ap.add_argument("--json", action="store_true", help="emit the model")
    ap.add_argument("--verbose", action="store_true",
                    help="every evidence record, not just the key ones")
    a = ap.parse_args(argv)
    be_polite()
    model = build_lineup(a.source, at=a.at)
    if a.json:
        json.dump(strip_private(model), sys.stdout, indent=2, default=str)
        sys.stdout.write("\n")
    else:
        report(model, only=a.service, verbose=a.verbose)
    return 0


if __name__ == "__main__":
    sys.exit(main())
