#!/usr/bin/env python3
"""E69 probe -- identify the transport on RF30's unparsed flow, and recover
the LLS table that m7_objects could not gunzip.

m7_objects reported flow 239.255.58.4:5004 (1364 of 1446 datagrams -- the bulk
of the multiplex) as "LCT-parses 0/200, MMTP length-field agreement 0/0 ==>
ROUTE/LCT" and then rejected all 1364.  A gate that classifies a flow and then
fails to parse a single packet of it has not classified anything.  Look at the
bytes.
"""
import binascii
import collections
import gzip
import json
import os
import struct
import sys
import zlib

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

import m7_objects as M7O  # noqa: E402


def hexdump(b, n=64):
    return " ".join("%02x" % c for c in b[:n])


def try_inflate(b):
    for name, fn in (("raw", lambda x: x),
                     ("gzip", gzip.decompress),
                     ("deflate-raw", lambda x: zlib.decompress(x, -15)),
                     ("zlib", zlib.decompress)):
        try:
            out = fn(b)
            if b"<" in out[:2048]:
                return name, out
        except Exception:                                       # noqa: BLE001
            continue
    # try skipping a few header bytes before the gzip magic
    i = b.find(b"\x1f\x8b")
    if i >= 0:
        try:
            return "gzip@%d" % i, gzip.decompress(b[i:])
        except Exception:                                       # noqa: BLE001
            pass
    return None, None


def mmtp_parse(p):
    """A/331 MMTP packet header (Table 9.2 / ISO 23008-1).

    byte0: V(2) C(1) FEC(2) r(1) X(1) R(1)   -- version in the top 2 bits
    byte1: RES(2) type(6)
    then packet_id(16), timestamp(32), packet_sequence_number(32)
    """
    if len(p) < 20:
        return None
    b0, b1 = p[0], p[1]
    ver = (b0 >> 6) & 3
    ptype = b1 & 0x3F
    pid, ts, seq = struct.unpack(">HII", p[2:12])
    return {"ver": ver, "type": ptype, "packet_id": pid,
            "timestamp": ts, "seq": seq}


def main():
    dg = os.path.join(HERE, "e69_smoke.dg")
    target = ("239.255.58.4", 5004)
    pkts = []
    lls = []
    for rec in M7O.read_dg(dg):
        src, dst, sport, dport, payload = (
            rec[0], rec[1], rec[2], rec[3], rec[4]) if len(rec) >= 5 else (
            None, None, None, None, None)
        if dst == target[0] and dport == target[1]:
            pkts.append(payload)
        elif dst == "224.0.23.60":
            lls.append(payload)
    print("flow %s:%d -- %d datagrams captured" % (target + (len(pkts),)))
    if pkts:
        lens = collections.Counter(len(p) for p in pkts)
        print("  payload lengths (top 6):", lens.most_common(6))
        for i, p in enumerate(pkts[:4]):
            print("  [%d] len=%4d  %s" % (i, len(p), hexdump(p, 32)))
        # MMTP hypothesis
        h = [mmtp_parse(p) for p in pkts]
        h = [x for x in h if x]
        vers = collections.Counter(x["ver"] for x in h)
        types = collections.Counter(x["type"] for x in h)
        pids = collections.Counter(x["packet_id"] for x in h)
        print("  MMTP hypothesis: version", dict(vers),
              " type", dict(types))
        print("  packet_id histogram:", pids.most_common(8))
        seqs = collections.defaultdict(list)
        for p, x in zip(pkts, h):
            seqs[x["packet_id"]].append(x["seq"])
        for pid, s in list(seqs.items())[:6]:
            d = [b - a for a, b in zip(s, s[1:])]
            inc1 = sum(1 for v in d if v == 1)
            print("    pid %5d: %5d pkts, seq %d..%d, %d/%d consecutive "
                  "deltas == 1" % (pid, len(s), s[0], s[-1], inc1, len(d)))

    print()
    print("LLS datagrams: %d" % len(lls))
    seen = {}
    for p in lls:
        if len(p) < 4:
            continue
        tid, gid, gcount, ver = p[0], p[1], p[2], p[3]
        key = (tid, ver)
        seen.setdefault(key, []).append(p)
    for (tid, ver), ps in sorted(seen.items()):
        body = ps[0][4:]
        name, out = try_inflate(body)
        print("  table_id %3d version %3d  x%d  body %d bytes  -> %s" %
              (tid, ver, len(ps), len(body), name or "NOT XML"))
        if out:
            fn = os.path.join(HERE, "e69_out", "e69_lls_%d_v%d.xml" % (tid, ver))
            open(fn, "wb").write(out)
            print("      wrote", os.path.basename(fn), len(out), "bytes")
            print("     ", out[:300].decode("utf-8", "replace").replace("\n", " "))
        else:
            print("      head:", hexdump(body, 48))


if __name__ == "__main__":
    main()
