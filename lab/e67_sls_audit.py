#!/usr/bin/env python3
"""E67: audit RF33 SLS + media objects banked by M13 (lab/plp1_out).

Question: are 105.1 WTTG / 104.1 WRC / 109.1 WUSA decodable by us, or
genuinely DRM-protected?  Read the signalling that is actually on the air
(SLT @protected/@drmSystemID, MPD ContentProtection, ISO-BMFF sinf/pssh)
and the media bytes themselves.  No radio; banked captures only.
"""
import json
import os
import re
import struct
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "plp1_out")

# svc-ip -> (serviceId, name, major.minor) from the RF33 SLT (bsid 540)
SERVICES = {
    "239_255_32_1_8321": ("1", "WHUT", "32.1"),
    "239_255_5_1_8051": ("3", "WTTG", "5.1"),
    "239_255_4_1_8041": ("4", "WRC", "4.1"),
    "239_255_9_1_8091": ("5", "WUSA", "9.1"),
}

# ---------------------------------------------------------------- SLS parsing


def split_multipart(blob):
    """ROUTE SLS is a multipart/related MIME body: return [(headers, body)]."""
    txt = blob.decode("utf-8", "replace")
    m = re.search(r"--(\S+)", txt)
    if not m:
        return [("", txt)]
    bound = m.group(1).rstrip("-")
    parts = []
    for chunk in txt.split("--" + bound):
        chunk = chunk.strip("\r\n-")
        if not chunk or chunk.strip() == "":
            continue
        if "\r\n\r\n" in chunk:
            hdr, body = chunk.split("\r\n\r\n", 1)
        elif "\n\n" in chunk:
            hdr, body = chunk.split("\n\n", 1)
        else:
            hdr, body = "", chunk
        parts.append((hdr.strip(), body))
    return parts


def sls_report(path):
    blob = open(path, "rb").read()
    rep = {"bytes": len(blob), "parts": []}
    for hdr, body in split_multipart(blob):
        loc = ""
        mm = re.search(r"Content-Location:\s*(\S+)", hdr, re.I)
        if mm:
            loc = mm.group(1)
        root = ""
        rm = re.search(r"<([A-Za-z0-9_:]+)[\s>]", body)
        if rm:
            root = rm.group(1)
        p = {"location": loc, "root": root, "bytes": len(body)}
        # DRM / encryption signalling inside the SLS
        cp = re.findall(
            r"<[^>]*ContentProtection[^>]*>", body, re.I)
        if cp:
            p["ContentProtection"] = cp
        for attr in ("schemeIdUri", "cenc:default_KID", "default_KID"):
            v = re.findall(attr + r'="([^"]*)"', body)
            if v:
                p[attr] = sorted(set(v))
        for tag in ("pssh", "ContentProtection", "cenc", "Widevine", "PlayReady"):
            if re.search(tag, body, re.I):
                p.setdefault("mentions", []).append(tag)
        # codecs advertised
        codecs = sorted(set(re.findall(r'codecs="([^"]*)"', body)))
        if codecs:
            p["codecs"] = codecs
        mimes = sorted(set(re.findall(r'mimeType="([^"]*)"', body)))
        if mimes:
            p["mimeTypes"] = mimes
        rep["parts"].append(p)
    return rep


# ------------------------------------------------------------- ISO-BMFF walk

CONTAINERS = {b"moov", b"trak", b"mdia", b"minf", b"stbl", b"moof", b"traf",
              b"mvex", b"edts", b"stsd", b"sinf", b"schi", b"encv", b"enca",
              b"mp4a", b"hvc1", b"hev1", b"ac-4", b"udta", b"dinf"}
# sample-entry boxes carry 8 bytes of reserved + fields before children
SAMPLE_ENTRY = {b"encv", b"enca", b"hvc1", b"hev1", b"mp4a", b"ac-4"}


def walk(buf, start, end, depth, acc, maxdepth=8):
    off = start
    while off + 8 <= end:
        size = struct.unpack(">I", buf[off:off + 4])[0]
        typ = buf[off + 4:off + 8]
        hdr = 8
        if size == 1:
            if off + 16 > end:
                break
            size = struct.unpack(">Q", buf[off + 8:off + 16])[0]
            hdr = 16
        elif size == 0:
            size = end - off
        if size < hdr or off + size > end:
            break
        acc.append((depth, typ.decode("latin1"), size))
        if typ in CONTAINERS and depth < maxdepth:
            cstart = off + hdr
            if typ in SAMPLE_ENTRY:
                # visual: 78 bytes of fields; audio: 28. skip conservatively
                cstart += 78 if typ in (b"encv", b"hvc1", b"hev1") else 28
            walk(buf, cstart, off + size, depth + 1, acc, maxdepth)
        off += size
    return acc


def entropy(b):
    if not b:
        return 0.0
    c = Counter(b)
    n = len(b)
    import math
    return -sum((v / n) * math.log2(v / n) for v in c.values())


def media_report(path):
    buf = open(path, "rb").read()
    acc = walk(buf, 0, len(buf), 0, [])
    types = Counter(t for _, t, _ in acc)
    rep = {
        "bytes": len(buf),
        "boxes_parsed": len(acc),
        "box_types": dict(types),
        "first_boxes": [t for _, t, _ in acc[:6]],
        "encrypted_sample_entries": [t for t in ("encv", "enca") if types.get(t)],
        "protection_boxes": {t: types[t] for t in ("sinf", "schm", "schi",
                                                  "tenc", "senc", "pssh",
                                                  "saiz", "saio", "frma")
                             if types.get(t)},
    }
    # scheme_type from schm if present
    schm = re.search(rb"schm.{4}(....)", buf, re.S)
    if types.get("schm"):
        i = buf.find(b"schm")
        rep["schm_scheme_type"] = buf[i + 8:i + 12].decode("latin1")
    if types.get("frma"):
        i = buf.find(b"frma")
        rep["frma_original_format"] = buf[i + 8:i + 12].decode("latin1")
    if types.get("pssh"):
        i = buf.find(b"pssh")
        rep["pssh_systemid"] = buf[i + 12:i + 28].hex()
    # mdat entropy (encrypted payload is ~8.0 bits/byte and has no start codes)
    i = buf.find(b"mdat")
    if i > 0:
        mdat = buf[i + 4:i + 4 + 65536]
        rep["mdat_entropy_bits"] = round(entropy(mdat), 3)
        rep["mdat_startcode_count"] = mdat.count(b"\x00\x00\x01")
    return rep


def main():
    result = {}
    for key, (sid, name, chan) in SERVICES.items():
        entry = {"serviceId": sid, "name": name, "channel": chan}
        for f in sorted(os.listdir(OUT)):
            if not f.startswith("obj_" + key):
                continue
            p = os.path.join(OUT, f)
            if f.endswith("FDT.xml"):
                continue
            if "_tsi0_" in f:
                entry["sls"] = {"file": f, **sls_report(p)}
            elif "_tsi10_" in f:
                entry["media"] = {"file": f, **media_report(p)}
        result[name] = entry
    print(json.dumps(result, indent=2))
    with open(os.path.join(HERE, "e67_sls_audit.json"), "w") as fh:
        json.dump(result, fh, indent=2)


if __name__ == "__main__":
    main()
