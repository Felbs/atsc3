#!/usr/bin/env python3
"""M6 -- RF33 PLP 0 and PLP 16, from cells to IP datagrams.

This is the rung M4 named and M5 stopped short of: the payload chain, end to
end, for a real ATSC 3.0 transmitter.

    cell pool (M5, gated)  ->  HTI Twisted Block de-interleaver (7.1.5.4)
      ->  bit de-interleaver (6.2, Annex B)  ->  64QAM-NUC demap (Annex C)
      ->  LDPC 16200 11/15  ->  BCH  ->  descramble  ->  Baseband Packets
      ->  ALP depacketization (A/330)  ->  IPv4/UDP  ->  LLS / ROUTE

THE REFEREE PROBLEM, AND HOW IT WAS SOLVED
------------------------------------------
M4 and M5 both recorded the same trap: between the cell stream and the LDPC
there is no intermediate check, so three rungs must all be right at once and
a wrong one cannot announce itself.

**PLP 16 breaks that deadlock and it was sitting in the L1-Detail all along.**
Its configuration is QPSK 2/15, `plp_size` 8100 = exactly ONE FEC Block,
`num_ti_blocks` 0 -> NTI 1, `num_fec_blocks` 0 -> 1 Block, so the Twisted
Block Interleaver is a 8100x1 array and **no twisting is possible.**  PLP 16
therefore tests the cell de-multiplex, the bit de-interleaver, the demapper,
the LDPC, the BCH and the descrambler with the time interleaver removed as a
variable -- and QPSK 2/15 decodes ~20 dB below the available MER, so a
failure there cannot be blamed on the link.

It decodes.  That converts the M5 wall from "three rungs, no referee" into a
single-variable question about 7.1.5.4, which is exactly what it turned out
to be.

Usage:
    python m6_payload.py hit_rf33.cs16 --rate 8e6 --frames 8
"""
from __future__ import annotations

import argparse
import collections
import gzip
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m6_bicm as B                                              # noqa: E402
import m6_cells as C                                             # noqa: E402
import m6_tbi as T                                               # noqa: E402

LLS_ADDR, LLS_PORT = "224.0.23.60", 4937
LLS_TABLES = {1: "SLT", 2: "RRT", 3: "SystemTime", 4: "AEAT",
              5: "OnscreenMessageNotification", 6: "CertificationData",
              7: "SignedMultiTable"}


# ---------------------------------------------------------------------------
# A/322 5.2.2 Baseband Packet header
# ---------------------------------------------------------------------------

def bb_split(b):
    """A/322 5.2.2 Figure 5.4 -> (pointer, payload).

    Base Field   MODE(1) + Pointer LSB(7)
    if MODE: byte 1 = Pointer MSB(6) + OFI(2), pointer is 13 bits
       OFI 00  no Optional/Extension Field           header 2
       OFI 01  short:  EXT_TYPE(3) + EXT_LEN(5)      header 3 + EXT_LEN
       OFI 10  long:   EXT_LEN is 13 bits            header 4 + EXT_LEN
       OFI 11  mixed:  NUM_EXT(3) + EXT_LEN(13)      header 4 + EXT_LEN

    Getting this wrong shifts the payload by a byte or thirty and shreds the
    ALP layer while leaving the FEC perfectly happy -- which is how it
    presented before it was fixed.  Pointer == all ones means no ALP packet
    starts in this Baseband Packet.
    """
    ptr, off = b[0] & 0x7F, 1
    if b[0] >> 7:
        ptr |= (b[1] >> 2) << 7
        ofi = b[1] & 0x03
        if ofi == 0:
            off = 2
        elif ofi == 1:
            off = 3 + (b[2] & 0x1F)
        else:
            off = 4 + (((b[3] << 5) | (b[2] & 0x1F)))
        none = ptr == 0x1FFF
    else:
        none = ptr == 0x7F
    return (None if none else ptr), b[off:]


# ---------------------------------------------------------------------------
# A/330 ALP
# ---------------------------------------------------------------------------

def alp_parse(stream, boundaries):
    """Walk the concatenated Baseband Packet payloads.

    `boundaries` are the (offset + pointer) positions at which A/322 5.2.2
    says an ALP packet STARTS.  They are used to RESYNC: a capture is a
    sequence of Frames and any gap between them desynchronises the walk, so
    on a malformed header we jump to the next signalled boundary instead of
    giving up.  Resyncs are counted and reported -- a silent resync would
    hide exactly the kind of chain error this project keeps finding.
    """
    bounds = sorted(set(boundaries))
    out, pos, resync = [], bounds[0] if bounds else 0, 0
    while pos + 2 <= len(stream):
        h = (stream[pos] << 8) | stream[pos + 1]
        ptype, pc, hm, ln = h >> 13, (h >> 12) & 1, (h >> 11) & 1, h & 0x7FF
        hdr = 2
        bad = pc == 1
        if not bad and hm:                      # additional header carries
            ln |= stream[pos + 2] << 11         # the length MSBs
            hdr = 3
        if not bad and (ln == 0 or pos + hdr + ln > len(stream)):
            bad = True
        if bad:
            nxt = [b for b in bounds if b > pos]
            if not nxt:
                break
            pos, resync = nxt[0], resync + 1
            continue
        out.append((ptype, bytes(stream[pos + hdr:pos + hdr + ln])))
        pos += hdr + ln
    return out, resync


def ip_udp(p):
    if len(p) < 28 or (p[0] >> 4) != 4 or p[9] != 17:
        return None
    ihl = (p[0] & 0xF) * 4
    src = ".".join(str(x) for x in p[12:16])
    dst = ".".join(str(x) for x in p[16:20])
    u = p[ihl:]
    if len(u) < 8:
        return None
    return dict(src=src, dst=dst, sport=(u[0] << 8) | u[1],
                dport=(u[2] << 8) | u[3], payload=u[8:])


# ---------------------------------------------------------------------------
# one frame
# ---------------------------------------------------------------------------

def decode_frame(y, t0, rep):
    """-> (plp16 result, list of 74 PLP-0 Baseband Packets, diagnostics)."""
    pool, info = C.cell_pool(y, t0, rep)
    reg, dv = C.constellation_regions(len(pool))
    alph = [B.points_for("64QAM", "11/15"), B.QPSK_POINTS]
    pool, _ = C.cpe_correct(pool, info["symbol_of"], alph, reg, dv)
    dummy = C.dummy_check(pool)

    ch16 = B.PlpChain(16200, "QPSK", "2/15")
    r16 = ch16.decode(pool[C.PLP16["start"]:C.PLP16["start"] + C.PLP16["size"]])
    r16["bb"] = ch16.baseband_packet(r16["bits"]) if r16["converged"] else None

    ch0 = B.PlpChain(16200, "64QAM", "11/15")
    plp0 = pool[:C.PLP0["size"]]
    nrows, ncols = ch0.cells_per_fec, C.PLP0["n_fec_max"] // C.PLP0["nti"]
    per_ti = len(plp0) // C.PLP0["nti"]
    packets, conv, bch = [], 0, 0
    for ti in range(C.PLP0["nti"]):
        seg = plp0[ti * per_ti:(ti + 1) * per_ti]
        for j in range(ncols):
            r = ch0.decode(T.fec_block(seg, j, nrows, ncols, 0))
            if r["converged"]:
                conv += 1
                bb = ch0.baseband_packet(r["bits"])
                bch += bb["bch_ok"]
                packets.append(bb["bytes"])
            else:
                packets.append(None)
    return r16, packets, dict(dummy=dummy, converged=conv, bch_ok=bch,
                              n_fec=C.PLP0["n_fec"], nrows=nrows, ncols=ncols)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("capture")
    ap.add_argument("--rate", type=float, required=True)
    ap.add_argument("--fmt")
    ap.add_argument("--start", type=float, default=0.0)
    ap.add_argument("--frames", type=int, default=4)
    ap.add_argument("--step", type=float, default=0.25)
    ap.add_argument("--outdir")
    ap.add_argument("--json")
    a = ap.parse_args()
    path = a.capture if os.path.isabs(a.capture) else os.path.join(HERE,
                                                                   a.capture)
    print("M6 -- RF33 PLP 0 / PLP 16 payload chain")
    print("=" * 72)

    print("\n  === GATE: the modules this stands on ===")
    ok_t = T.selftest(verbose=True)
    print()
    ok_b = B.selftest(verbose=False)
    print(f"    {'PASS' if ok_b else 'FAIL'}  m6_bicm selftest (round trips, "
          f"6.2.3.1 worked example, Type B parity gate, 5 controls)")

    # ONE load, then step by the EXACT Frame period so the Baseband Packet
    # stream is genuinely contiguous.  Stepping by a rounded wall-clock
    # interval drops ~1% of every Frame and shreds the ALP layer.
    span = 0.05 + a.frames * C.FRAME_SAMPLES / 6.912e6
    y_all = C.load_frame(path, a.rate, fmt=a.fmt, start_sec=a.start,
                         span_sec=span * (a.rate / 6.912e6))
    print(f"\n  Frame period {C.FRAME_SAMPLES} samples "
          f"({C.FRAME_SAMPLES/6.912e6*1000:.2f} ms); loaded "
          f"{len(y_all)} samples for {a.frames} Frames")
    stream, all_pkts, frames, bounds = bytearray(), [], [], []
    base = None
    for fi in range(a.frames):
        if (fi + 1) * C.FRAME_SAMPLES + 400000 > len(y_all):
            print(f"    (capture exhausted after {fi} Frames)")
            break
        y = y_all
        t0, coh = C.fine_timing(y, centre=(base or C.BOOTSTRAP)
                                + fi * C.FRAME_SAMPLES)
        if base is None:
            base = t0
        rep = {}
        r16, pkts, diag = decode_frame(y, t0, rep)
        good = [p for p in pkts if p is not None]
        all_pkts += good
        for p in good:
            ptr, pay = bb_split(p)
            if ptr is not None:
                bounds.append(len(stream) + ptr)
            stream += pay
        frames.append(dict(frame=fi, t0=int(t0), coherence=coh,
                           plp16_converged=r16["converged"],
                           plp16_bch=bool(r16["bb"] and r16["bb"]["bch_ok"]),
                           **{k: v for k, v in diag.items() if k != "dummy"},
                           dummy_agreement=diag["dummy"]["sign_agreement"]))
        print(f"\n  ---- frame {fi} (t0 {t0}, pilot coherence {coh:.4f}) ----")
        print(f"    dummy cells (A/322 7.2.6.5): {diag['dummy']['n_real']} of "
              f"{diag['dummy']['n']} real-valued, sign agreement "
              f"{diag['dummy']['sign_agreement']:.4f}")
        bchs = "n/a" if not r16["bb"] else str(r16["bb"]["bch_syndrome"])
        print(f"    PLP 16  QPSK 2/15, 1 FEC Block, NO twisting possible:  "
              f"LDPC {'CONVERGED' if r16['converged'] else 'FAILED'}, "
              f"BCH syndrome {bchs}")
        print(f"    PLP 0   64QAM-NUC 11/15, TBI {diag['nrows']}x"
              f"{diag['ncols']} x NTI 2:  "
              f"**{diag['converged']}/{diag['n_fec']} FEC Blocks converged, "
              f"{diag['bch_ok']}/{diag['n_fec']} BCH syndromes zero**")

    print(f"\n  === ALP depacketization (A/330) ===")
    alps, resync = alp_parse(stream, bounds)
    kinds = collections.Counter(t for t, _ in alps)
    print(f"    {len(all_pkts)} Baseband Packets -> {len(stream)} payload "
          f"bytes -> {len(alps)} ALP packets, {resync} resyncs")
    print(f"    packet_type histogram: "
          f"{ {('IPv4' if k == 0 else k): v for k, v in kinds.items()} }")

    flows, lls, route = collections.Counter(), [], collections.defaultdict(list)
    for t, p in alps:
        if t != 0:
            continue
        d = ip_udp(p)
        if not d:
            continue
        flows[(d["src"], d["dst"], d["dport"])] += 1
        if d["dst"] == LLS_ADDR and d["dport"] == LLS_PORT:
            lls.append(d["payload"])
        else:
            route[(d["dst"], d["dport"])].append(d["payload"])

    print(f"\n  === IPv4/UDP flows ===")
    for (s, dst, dp), n in flows.most_common(20):
        tag = "   <== LLS (well-known)" if (dst, dp) == (LLS_ADDR,
                                                         LLS_PORT) else ""
        print(f"    {s:15s} -> {dst:15s} :{dp:<6d} {n:5d} datagrams{tag}")

    slt_xml = None
    print(f"\n  === LLS / SLT ===")
    if not lls:
        print(f"    NO datagram to {LLS_ADDR}:{LLS_PORT} in "
              f"{a.frames} frames.  The service list is NOT available from "
              f"this slice; LLS repeats on its own cycle, so try more frames.")
    for d in lls[:4]:
        tid, gid, gc, ver = d[0], d[1], d[2], d[3]
        name = LLS_TABLES.get(tid, f"reserved({tid})")
        print(f"    LLS_table_id {tid} ({name}), group {gid}, "
              f"group_count {gc + 1}, version {ver}, {len(d)} bytes")
        try:
            xml = gzip.decompress(d[4:]).decode("utf-8", "replace")
        except Exception as e:                                  # noqa: BLE001
            print(f"      gunzip failed: {e}")
            continue
        if tid == 1:
            slt_xml = xml
        print("      " + xml[:2000].replace("\n", "\n      "))

    if a.outdir:
        os.makedirs(a.outdir, exist_ok=True)
        for (dst, dp), pl in route.items():
            with open(os.path.join(a.outdir, f"route_{dst}_{dp}.bin"),
                      "wb") as fh:
                for x in pl:
                    fh.write(x)
        if slt_xml:
            with open(os.path.join(a.outdir, "slt.xml"), "w",
                      encoding="utf-8") as fh:
                fh.write(slt_xml)
        print(f"\n  wrote ROUTE payloads to {a.outdir}")

    out = dict(capture=os.path.basename(path), frames=frames,
               gates=dict(tbi=ok_t, bicm=ok_b),
               n_alp=len(alps), flows={f"{k[0]}->{k[1]}:{k[2]}": v
                                       for k, v in flows.items()},
               lls=len(lls), slt=bool(slt_xml))
    dest = a.json or os.path.join(
        HERE, "m6_payload_" + os.path.splitext(os.path.basename(path))[0]
        + ".json")
    with open(dest, "w") as fh:
        json.dump(out, fh, indent=2, default=float)
    print(f"\n  wrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
