#!/usr/bin/env python3
"""E69 -- recover RF30's SLT from the LLS_table_id 254 multi-table container.

m7_objects reported table 254 as "[payload did not gunzip to XML]" because it
tried to gunzip the container itself.  The container is not gzip; it HOLDS
gzip members.  Its body starts:

    02 | 01 a6 02 15 | 1f 8b 08 08 .. 00 03 53 4c 54 00 ...
    ^^   ^^ ^^ ^^^^^   ^^^^^                  ^^^^^^^^
    |    |  |  length  gzip magic             gzip's stored filename: "SLT"
    |    |  version
    |    LLS_table_id 1 = SLT
    table count

so it is {count, then count x (table_id, version, length_be16, payload)}.  The
gzip members even carry their original filenames ("SLT", "UserDefined"), which
is an independent confirmation of the field decode -- the id byte says 1 and
the payload says "SLT" without being asked.
"""
import gzip
import io
import json
import os
import struct
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

import m7_objects as M7O  # noqa: E402

NAMES = {1: "SLT", 2: "RRT", 3: "SystemTime", 4: "AEAT", 5: "OnscreenMessage",
         6: "CertificationData", 254: "multi-table", 255: "UserDefined"}


def gz_name(b):
    """The stored filename in a gzip member, if FNAME is set."""
    if len(b) < 10 or b[0] != 0x1F or b[1] != 0x8B:
        return None
    if not (b[3] & 0x08):
        return None
    end = b.find(b"\x00", 10)
    return b[10:end].decode("latin1") if end > 0 else None


def parse_multi(body):
    """Return [(table_id, version, raw, xml_or_None, gzip_name)]."""
    out = []
    n = body[0]
    o = 1
    for _ in range(n):
        if o + 4 > len(body):
            break
        tid, ver = body[o], body[o + 1]
        ln = struct.unpack(">H", body[o + 2:o + 4])[0]
        o += 4
        raw = body[o:o + ln]
        o += ln
        name = gz_name(raw)
        xml = None
        try:
            xml = gzip.decompress(raw)
        except Exception:                                       # noqa: BLE001
            try:
                xml = gzip.GzipFile(fileobj=io.BytesIO(raw)).read()
            except Exception:                                   # noqa: BLE001
                xml = None
        out.append((tid, ver, raw, xml, name))
    return out, o


def main():
    dgp = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "e69_smoke.dg")
    outdir = os.path.join(HERE, "e69_out")
    os.makedirs(outdir, exist_ok=True)
    seen = {}
    for src, dst, sp, dp, p in M7O.read_dg(dgp):
        if dst != "224.0.23.60" or len(p) < 5:
            continue
        seen.setdefault((p[0], p[3]), p)
    wrote = []
    for (tid, ver), p in sorted(seen.items()):
        body = p[4:]
        print("LLS_table_id %d (%s) version %d, %d bytes"
              % (tid, NAMES.get(tid, "?"), ver, len(p)))
        if tid != 254:
            try:
                xml = gzip.decompress(body)
            except Exception:                                   # noqa: BLE001
                xml = None
            if xml:
                fn = os.path.join(outdir, "e69_lls%d_v%d.xml" % (tid, ver))
                open(fn, "wb").write(xml)
                wrote.append(fn)
                print("   -> %s (%d bytes)" % (os.path.basename(fn), len(xml)))
            continue
        tabs, consumed = parse_multi(body)
        print("   multi-table: %d entries, %d/%d bytes consumed"
              % (len(tabs), consumed, len(body)))
        for t_id, t_ver, raw, xml, name in tabs:
            print("     id %3d (%-12s) ver %3d  %5d bytes  gzip-name %r  -> %s"
                  % (t_id, NAMES.get(t_id, "?"), t_ver, len(raw), name,
                     "%d bytes XML" % len(xml) if xml else "FAILED"))
            if xml:
                fn = os.path.join(outdir, "e69_lls%d_%s_v%d.xml"
                                  % (t_id, NAMES.get(t_id, "x"), t_ver))
                open(fn, "wb").write(xml)
                wrote.append(fn)
    print()
    for fn in wrote:
        print("=" * 70)
        print(os.path.basename(fn))
        print("=" * 70)
        sys.stdout.buffer.write(open(fn, "rb").read())
        print()


if __name__ == "__main__":
    main()
