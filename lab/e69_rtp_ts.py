#!/usr/bin/env python3
"""E69 -- RF30's main flow is RTP/MP2T, not MMTP and not ROUTE.

m7_objects' transport gate called 239.255.58.4:5004 "ROUTE/LCT" and then
rejected 1364 of 1364 datagrams.  The bytes say something the gate has no
branch for:

    80 21 58 ed 12 af 8a b1 7f bf b1 c5 | 47 1f ff 10 ...
    ^^ ^^                                 ^^
    |  RTP payload type 33 = MP2T         MPEG-2 TS sync byte
    RTP version 2

1328 bytes = 12-byte RTP header + 7 x 188-byte TS packets, exactly.  So RF30's
core layer carries a plain MPEG-2 Transport Stream over RTP -- the datacast
shape, not ATSC 3.0's own MMTP or ROUTE.  (My first probe called this "MMTP
version 2 type 33" because RTP's V/PT fields land in the same bits as MMTP's
version/type.  Two different protocols agreeing on a byte layout is exactly
how a transport gate gets fooled; the discriminator is the 0x47 sync byte at
offset 12 and the 7x188 arithmetic.)

Writes a .ts and reports the PAT/PMT/PID census.
"""
import collections
import os
import struct
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

import m7_objects as M7O  # noqa: E402

TS = 188


def rtp_split(p):
    """Return the RTP payload if this is a well-formed RTP/MP2T packet."""
    if len(p) < 12:
        return None
    if (p[0] >> 6) != 2:                 # version must be 2
        return None
    if (p[1] & 0x7F) != 33:              # payload type 33 = MP2T
        return None
    cc = p[0] & 0x0F
    ext = (p[0] >> 4) & 1
    off = 12 + 4 * cc
    if ext:
        if len(p) < off + 4:
            return None
        off += 4 + 4 * struct.unpack(">H", p[off + 2:off + 4])[0]
    pay = p[off:]
    if len(pay) % TS or not pay or pay[0] != 0x47:
        return None
    return struct.unpack(">H", p[2:4])[0], pay      # (seq, payload)


def main():
    dgp = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "e69_smoke.dg")
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "e69_rf30.ts")
    flow = sys.argv[3] if len(sys.argv) > 3 else "239.255.58.4"
    port = int(sys.argv[4]) if len(sys.argv) > 4 else 5004

    got, bad, seqs = [], 0, []
    for src, dst, sp, dp, payload in M7O.read_dg(dgp):
        if dst != flow or dp != port:
            continue
        r = rtp_split(payload)
        if r is None:
            bad += 1
            continue
        seqs.append(r[0])
        got.append(r[1])
    print("RTP/MP2T on %s:%d" % (flow, port))
    print("  %d datagrams parsed as RTP/MP2T, %d rejected" % (len(got), bad))
    if not got:
        return 1
    # sequence continuity (16-bit wrapping)
    gaps = 0
    for a, b in zip(seqs, seqs[1:]):
        if ((b - a) & 0xFFFF) != 1:
            gaps += 1
    print("  RTP seq %d..%d, %d discontinuities in %d steps"
          % (seqs[0], seqs[-1], gaps, len(seqs) - 1))

    blob = b"".join(got)
    npkt = len(blob) // TS
    print("  %d TS packets (%.2f MB)" % (npkt, len(blob) / 1e6))

    # PID census + sync integrity
    pids = collections.Counter()
    sync_ok = 0
    scrambled = collections.Counter()
    for i in range(npkt):
        p = blob[i * TS:(i + 1) * TS]
        if p[0] != 0x47:
            continue
        sync_ok += 1
        pid = ((p[1] & 0x1F) << 8) | p[2]
        pids[pid] += 1
        scrambled[(p[3] >> 6) & 3] += 1
    print("  sync byte present on %d/%d packets" % (sync_ok, npkt))
    print("  transport_scrambling_control: %s"
          % dict(scrambled))
    print("  PID census (top 12):")
    for pid, n in pids.most_common(12):
        tag = ""
        if pid == 0:
            tag = "  <== PAT"
        elif pid == 0x1FFF:
            tag = "  <== null padding"
        elif pid == 0x1FFB:
            tag = "  <== ATSC PSIP base"
        print("    pid 0x%04X (%5d): %6d packets%s" % (pid, pid, n, tag))

    with open(out, "wb") as fh:
        fh.write(blob)
    print("  wrote %s" % out)

    # ---- PAT ----------------------------------------------------------
    for i in range(npkt):
        p = blob[i * TS:(i + 1) * TS]
        pid = ((p[1] & 0x1F) << 8) | p[2]
        if pid != 0 or not (p[1] & 0x40):
            continue
        o = 4 + (1 + p[4] if (p[3] & 0x20) else 0)
        o += 1 + p[o]                      # pointer_field
        if p[o] != 0x00:
            continue
        slen = ((p[o + 1] & 0x0F) << 8) | p[o + 2]
        body = p[o + 8:o + 3 + slen - 4]
        print("  PAT: %d program entries" % (len(body) // 4))
        for k in range(0, len(body) - 3, 4):
            prog = struct.unpack(">H", body[k:k + 2])[0]
            pmt = struct.unpack(">H", body[k + 2:k + 4])[0] & 0x1FFF
            print("    program %d -> PMT pid 0x%04X" % (prog, pmt))
        break
    else:
        print("  PAT: none found")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
