#!/usr/bin/env python3
"""atsc3_route.py -- the ROUTE/DASH transport layer for an ATSC 3.0 service.

The MMTP path (tools/atsc3_tv.py) rebuilds a service whose media rides MMTP.
A ROUTE service is different in EXACTLY two places and identical everywhere
else, and this module is those two places:

  1. SLS PARSE.  The Service Layer Signalling arrives as ONE object on the
     service's LCT tsi=0 -- a MIME `multipart/related` bundle carrying the
     MPD (DASH) and the S-TSID (ROUTE).  The MPD says which Representations
     exist and their @codecs / SegmentTemplate; the S-TSID binds each
     Representation to an LCT `tsi` and names its init-segment and media
     `$TOI$` file template.  Parse both, JOIN them by Representation id, and
     you know: video-rep -> which tsi, audio-rep -> which tsi.

  2. SEGMENT ASSEMBLE.  A DASH Representation is a byte stream:
        init.mp4  +  media-$N$.mp4  +  media-$(N+1)$.mp4  +  ...
     in segment-number order.  Over ROUTE the init segment is the object
     with TOI = 0xFFFFFFFF and each media segment is the object with
     TOI = the DASH `$Number$`.  So "assemble the Representation" is: take
     the RouteFlow objects for this tsi, put the init first, then the media
     objects sorted by TOI, and concatenate.  That is a playable fMP4.

Everything upstream -- the physical layer, the ALP walk, the LCT reassembly
into (tsi, toi) objects -- is UNCHANGED and reused: this imports
`m7_objects.read_dg` + `RouteFlow` verbatim.  Nothing here re-derives the
transport it sits on.

LIVE NOTE.  In the live chain this slots beside the MMTP path exactly where
`atsc3_tv.py` chooses a transport per flow (m7_objects already gates
ROUTE-vs-MMTP per flow with the LCT-parse / MMTP-length-agreement referee).
For a ROUTE service the live loop would, each time a new SLS object or a new
media TOI completes: feed datagrams to a persistent RouteFlow, and whenever a
media object for a rep's tsi reaches `state == "complete"`, append its bytes
to that rep's lane file (init once, then media in TOI order) -- the same
"drop-a-frame-not-a-sample, append-complete-segments-only" lane model the
MMTP live player already uses.  This file's `assemble_representation` is that
append step; only the driving loop differs.

Usage:
    # parse the SLS bundle and print the service structure
    python tools/atsc3_route.py sls lab/plp1_out/obj_..._8321_tsi0_toi4653059.bin

    # assemble a Representation's fMP4 from a banked datagram cache
    python tools/atsc3_route.py assemble --dg lab/m13_plp1_48.dg \\
        --sls lab/plp1_out/obj_..._8321_tsi0_toi4653059.bin \\
        --flow 239.255.32.1:8321 --rep TS-4000_1_video --out video.mp4
"""
from __future__ import annotations

import argparse
import os
import re
import sys
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
LAB = os.path.join(ROOT, "lab")
sys.path.insert(0, LAB)
sys.path.insert(0, ROOT)

import m7_objects as OBJ                                          # noqa: E402

INIT_TOI = 0xFFFFFFFF            # A/331: the ROUTE init segment's reserved TOI


# ---------------------------------------------------------------------------
# 1. SLS parse  (multipart -> MPD + S-TSID -> a per-Representation table)
# ---------------------------------------------------------------------------

def split_multipart(data):
    """A `multipart/related` bundle -> {Content-Location: body_bytes}.

    The bundle is a tiny MIME document: a preamble naming the boundary, then
    parts separated by `--boundary`, each with its own headers and body.  We
    key the parts by their `Content-Location` because the envelope refers to
    them by that name (metadataURI), not by order.
    """
    if isinstance(data, str):
        data = open(data, "rb").read()
    m = re.search(rb'boundary="?([^";\r\n]+)"?', data)
    if not m:
        raise ValueError("no MIME boundary in the SLS bundle")
    boundary = b"--" + m.group(1)
    parts = data.split(boundary)
    out = {}
    for part in parts:
        part = part.strip(b"\r\n")
        if not part or part == b"--":
            continue
        if b"\r\n\r\n" not in part:
            continue
        head, body = part.split(b"\r\n\r\n", 1)
        loc = None
        for line in head.split(b"\r\n"):
            if line.lower().startswith(b"content-location:"):
                loc = line.split(b":", 1)[1].strip().decode("ascii", "replace")
        if loc:
            out[loc] = body.rstrip(b"\r\n")
    return out


def _localname(tag):
    return tag.rsplit("}", 1)[-1]


def parse_mpd(mpd_bytes):
    """DASH MPD -> [Representation dicts] with SegmentTemplate resolved.

    A Representation inherits its AdaptationSet's @codecs / @mimeType /
    @contentType and (usually) its SegmentTemplate; we resolve that
    inheritance here so each returned dict stands on its own.
    """
    root = ET.fromstring(mpd_bytes)
    reps = []
    for period in root.iter():
        if _localname(period.tag) != "Period":
            continue
        for aset in period:
            if _localname(aset.tag) != "AdaptationSet":
                continue
            a = aset.attrib
            aset_tmpl = None
            for child in aset:
                if _localname(child.tag) == "SegmentTemplate":
                    aset_tmpl = child.attrib
            for rep in aset:
                if _localname(rep.tag) != "Representation":
                    continue
                r = rep.attrib
                tmpl = dict(aset_tmpl or {})
                for child in rep:
                    if _localname(child.tag) == "SegmentTemplate":
                        tmpl.update(child.attrib)
                reps.append(dict(
                    rep_id=r.get("id"),
                    content_type=a.get("contentType"),
                    codecs=r.get("codecs", a.get("codecs")),
                    mime_type=r.get("mimeType", a.get("mimeType")),
                    lang=a.get("lang"),
                    bandwidth=r.get("bandwidth"),
                    width=a.get("width"), height=a.get("height"),
                    frame_rate=a.get("frameRate"),
                    audio_sampling_rate=r.get("audioSamplingRate"),
                    init=tmpl.get("initialization"),
                    media=tmpl.get("media"),
                    start_number=tmpl.get("startNumber"),
                    duration=tmpl.get("duration"),
                    timescale=tmpl.get("timescale")))
    return reps


def parse_stsid(stsid_bytes):
    """ROUTE S-TSID -> [LS dicts]: tsi <-> repId + init/media file templates."""
    root = ET.fromstring(stsid_bytes)
    out = []
    for rs in root.iter():
        if _localname(rs.tag) != "RS":
            continue
        rs_ip = rs.attrib.get("dIpAddr")
        rs_port = rs.attrib.get("dPort")
        for ls in rs:
            if _localname(ls.tag) != "LS":
                continue
            tsi = ls.attrib.get("tsi")
            rep_id = content_type = init_loc = file_template = None
            for node in ls.iter():
                ln = _localname(node.tag)
                if ln == "FDT-Instance":
                    for k, v in node.attrib.items():
                        if _localname(k) == "fileTemplate":
                            file_template = v
                elif ln == "File":
                    init_loc = node.attrib.get("Content-Location", init_loc)
                elif ln == "MediaInfo":
                    rep_id = node.attrib.get("repId", rep_id)
                    content_type = node.attrib.get("contentType", content_type)
            out.append(dict(tsi=int(tsi) if tsi is not None else None,
                            rep_id=rep_id, content_type=content_type,
                            init_location=init_loc,
                            file_template=file_template,
                            ip=rs_ip, port=rs_port))
    return out


def parse_sls(sls_bytes):
    """The whole SLS bundle -> (service_id, [joined-Representation dicts]).

    Each returned dict is one Representation with its DASH view (codecs,
    SegmentTemplate) and its ROUTE view (tsi, init/media file templates)
    JOINED on the Representation id.
    """
    parts = split_multipart(sls_bytes)
    mpd = next((v for k, v in parts.items() if k.endswith(".xml")
                and b"<MPD" in v), None)
    stsid = next((v for k, v in parts.items()
                  if b"<S-TSID" in v), None)
    usbd = next((v for k, v in parts.items()
                 if b"BundleDescriptionROUTE" in v), None)
    if mpd is None or stsid is None:
        raise ValueError("SLS bundle missing MPD or S-TSID "
                         f"(parts: {sorted(parts)})")
    reps = parse_mpd(mpd)
    ls_list = parse_stsid(stsid)
    by_rep = {ls["rep_id"]: ls for ls in ls_list}

    service_id = None
    if usbd is not None:
        m = re.search(rb'serviceId="(\d+)"', usbd)
        if m:
            service_id = int(m.group(1))

    joined = []
    for r in reps:
        ls = by_rep.get(r["rep_id"], {})
        j = dict(r)
        j["tsi"] = ls.get("tsi")
        j["init_location"] = ls.get("init_location")
        j["file_template"] = ls.get("file_template")
        j["ip"] = ls.get("ip")
        j["port"] = ls.get("port")
        joined.append(j)
    return service_id, joined


def print_service(service_id, reps, ip=None, port=None):
    print("=" * 74)
    print(f"  ROUTE SERVICE  serviceId={service_id}"
          + (f"   flow {ip}:{port}" if ip else ""))
    print("=" * 74)
    for r in reps:
        seg_s = None
        if r.get("duration") and r.get("timescale"):
            seg_s = float(r["duration"]) / float(r["timescale"])
        print(f"  Representation {r['rep_id']!r}")
        print(f"      contentType {r['content_type']}   codecs {r['codecs']}"
              f"   mime {r['mime_type']}"
              + (f"   lang {r['lang']}" if r.get("lang") else ""))
        dims = ""
        if r.get("width"):
            dims = f"   {r['width']}x{r['height']} @ {r.get('frame_rate')}"
        elif r.get("audio_sampling_rate"):
            dims = f"   {r['audio_sampling_rate']} Hz"
        print(f"      LCT tsi = {r['tsi']}{dims}"
              + (f"   segment {seg_s:.3f}s" if seg_s else ""))
        print(f"      init  {r['init_location']}")
        print(f"      media {r['file_template']}"
              f"   (init TOI = 0x{INIT_TOI:08X})")
    print("=" * 74)


# ---------------------------------------------------------------------------
# 2. Segment assemble  (RouteFlow objects -> a Representation's fMP4)
# ---------------------------------------------------------------------------

def route_objects(dg_path, flow):
    """Reassemble a flow's LCT objects.  -> {(tsi,toi): assembled-object}.

    Verbatim reuse of the M7 transport layer: read_dg boundary-preserving
    datagrams, feed the flow's packets to RouteFlow, assemble.
    """
    dst, port = flow.rsplit(":", 1)
    port = int(port)
    recs = OBJ.read_dg(dg_path)
    fl = OBJ.RouteFlow(flow)
    n = 0
    for src, d, sp, dp, pay in recs:
        if d == dst and dp == port:
            fl.feed(pay)
            n += 1
    objs = {(o["tsi"], o["toi"]): o for o in fl.assemble()}
    return objs, n


def assemble_representation(objs, tsi, require_complete=True):
    """init + media segments in TOI order -> (fmp4_bytes, manifest).

    `objs` is the {(tsi,toi): object} map from route_objects.  The init
    segment is the object at TOI 0xFFFFFFFF; media segments are every other
    TOI for this tsi, laid down in ascending TOI order -- which IS DASH
    segment-number order because the ROUTE fileTemplate maps $Number$ -> TOI.

    A media object is used only if COMPLETE (its reassembled length equals the
    transport-declared length with no gaps); a segment with a hole at the
    front has zero bytes where its `styp` should be and would poison the
    concatenation.  Partial/leading/trailing segments are reported, not used.
    """
    init = objs.get((tsi, INIT_TOI))
    media_keys = sorted(k for k in objs
                        if k[0] == tsi and k[1] != INIT_TOI)
    used, skipped = [], []
    buf = bytearray()
    if init is not None and init["state"] == "complete":
        buf += init["data"]
    manifest = dict(tsi=tsi, init_bytes=(len(init["data"]) if init else 0),
                    init_state=(init["state"] if init else "MISSING"),
                    segments=[])
    for (t, toi) in media_keys:
        o = objs[(t, toi)]
        entry = dict(toi=toi, got=o["got"], declared=o["declared"],
                     state=o["state"], first_box=_first_box(o["data"]))
        if require_complete and o["state"] != "complete":
            skipped.append(entry)
            continue
        buf += o["data"]
        used.append(entry)
        manifest["segments"].append(entry)
    manifest["n_used"] = len(used)
    manifest["n_skipped"] = len(skipped)
    manifest["skipped"] = skipped
    manifest["used_tois"] = [e["toi"] for e in used]
    return bytes(buf), manifest


def _first_box(data):
    if len(data) >= 8 and data[4:8].isalpha():
        return data[4:8].decode("ascii", "replace")
    return ""


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(prog="atsc3_route")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_sls = sub.add_parser("sls", help="parse an SLS bundle and print it")
    p_sls.add_argument("bundle")

    p_asm = sub.add_parser("assemble", help="assemble a Representation fMP4")
    p_asm.add_argument("--dg", required=True)
    p_asm.add_argument("--sls", required=True)
    p_asm.add_argument("--flow", required=True, help="dst:port")
    p_asm.add_argument("--rep", required=True, help="Representation id")
    p_asm.add_argument("--out", required=True)
    p_asm.add_argument("--allow-partial", action="store_true")

    a = ap.parse_args(argv)

    if a.cmd == "sls":
        sid, reps = parse_sls(a.bundle)
        ip = reps[0].get("ip") if reps else None
        port = reps[0].get("port") if reps else None
        print_service(sid, reps, ip, port)
        return 0

    if a.cmd == "assemble":
        sid, reps = parse_sls(a.sls)
        rep = next((r for r in reps if r["rep_id"] == a.rep), None)
        if rep is None:
            print(f"  no Representation {a.rep!r}; have "
                  f"{[r['rep_id'] for r in reps]}")
            return 1
        tsi = rep["tsi"]
        print(f"  Representation {a.rep!r}  ->  LCT tsi {tsi}  "
              f"({rep['content_type']}, {rep['codecs']})")
        objs, n = route_objects(a.dg, a.flow)
        print(f"  flow {a.flow}: {n} datagrams -> {len(objs)} objects")
        fmp4, man = assemble_representation(
            objs, tsi, require_complete=not a.allow_partial)
        out = a.out if os.path.isabs(a.out) else os.path.join(os.getcwd(),
                                                              a.out)
        with open(out, "wb") as fh:
            fh.write(fmp4)
        print(f"  init {man['init_state']} ({man['init_bytes']} B); "
              f"used {man['n_used']} media segments "
              f"(TOI {man['used_tois']}), skipped {man['n_skipped']}")
        for s in man["skipped"]:
            print(f"      skipped TOI {s['toi']}: {s['state']} "
                  f"({s['got']}/{s['declared']}, box {s['first_box']!r})")
        print(f"  wrote {out}  ({len(fmp4)} bytes)")
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
