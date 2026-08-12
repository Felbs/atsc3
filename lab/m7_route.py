#!/usr/bin/env python3
"""m7_route.py -- RF33 payload chain, part 1: Frames -> IP DATAGRAMS on disk.

M6 got as far as concatenating every UDP payload of a flow into one .bin.
That threw away the one thing the transport layer needs: **where each
datagram ended.**  LCT and MMTP both carry an explicit header but NOT an
explicit payload length -- the payload length is the UDP length.  Concatenate
the payloads and the object reassembler has nothing to walk.

So this stage re-runs the M6 chain and writes a datagram-boundary-preserving
cache (`.dg`), which `m7_objects.py` then walks.  Three things are new:

  1. **ALP payload_configuration == 1** (A/330 5.2) -- segmentation AND
     concatenation.  M6 treated both as malformed and resynced past them.
     Small packets (the LLS/SLT is ~300 bytes) are exactly what a broadcaster
     concatenates, so this is the leading candidate for "no SLT in 20 frames".
  2. **IPv4 fragment reassembly.**  M6's ip_udp() read bytes 0..7 of every
     fragment as a UDP header, so a fragmented datagram produced one garbled
     record and lost the rest.
  3. **Chunked decode**, so a 120 s capture can be walked without holding
     8 GB of complex128.

Every new reading is GATED, not assumed: the concatenation reading is scored
by "do all the components parse as IPv4 with self-consistent lengths", and
the score is printed.  A reading that cannot beat its control is reported as
failed, not quietly used.

Usage:
    python m7_route.py long_rf33.cs16 --rate 8e6 --frames 120 --out m7.dg
"""
from __future__ import annotations

import argparse
import bisect
import collections
import json
import os
import struct
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m6_cells as C                                             # noqa: E402
import m6_payload as P6                                          # noqa: E402

REC = struct.Struct("<4s4sHHI")
BBHDR = struct.Struct("<dIII")          # start_sec, n_frames, n_stream, n_bnd


def bb_write(fh, s0, fr, stream, bounds):
    fh.write(BBHDR.pack(s0, len(fr), len(stream), len(bounds)))
    fh.write(bytes(stream))
    fh.write(struct.pack(f"<{len(bounds)}I", *bounds))
    fh.write(json.dumps(fr).encode() + b"\n")


def bb_read(path):
    d = open(path, "rb").read()
    o = 0
    while o + BBHDR.size <= len(d):
        s0, nf, ns, nb = BBHDR.unpack_from(d, o)
        o += BBHDR.size
        stream = d[o:o + ns]
        o += ns
        bounds = list(struct.unpack_from(f"<{nb}I", d, o))
        o += 4 * nb
        e = d.index(b"\n", o)
        fr = json.loads(d[o:e])
        o = e + 1
        yield s0, fr, stream, bounds


# ---------------------------------------------------------------------------
# A/330 5.2 -- the full ALP header, including payload_configuration == 1
# ---------------------------------------------------------------------------

class AlpWalker:
    """Walk the concatenated Baseband Packet payload stream.

    Returns (packet_type, payload) tuples in order, plus a statistics dict.
    Segmented ALP packets are reassembled by segment_sequence_number and only
    emitted when the last_segment_indicator arrives.

    TWO REFEREES, because a length field that is read wrong is silent:

    A. **The anchor contradiction.**  Every Baseband Packet header carries a
       pointer to the first ALP packet that STARTS inside it.  So an anchor
       strictly inside a candidate ALP packet's span is a proof that the
       candidate is wrong -- an ALP packet cannot begin inside another one.
       This is free, it comes from a different layer, and it catches exactly
       the failure that a plausible-looking length cannot announce.

    B. **The IPv4 gate.**  packet_type 0 payloads must begin 0x4_ and carry a
       total_length that fits.  A cross-layer check, again from outside.

    Without (A), one mis-read `header_mode == 1` length swallowed a quarter of
    a megabyte of stream in a single "valid" packet and threw away 62% of the
    media flow while reporting one resync.
    """

    MAX_ALP = 65600            # an IP datagram cannot exceed 65535 + headers

    def __init__(self, stream, boundaries):
        self.s = stream
        self.bounds = sorted(set(boundaries))
        self.bset = set(self.bounds)
        self.stats = collections.Counter()
        self.segbuf = {}

    def _anchor_clean(self, start, end):
        """No anchor may fall strictly inside [start, end)."""
        i = bisect.bisect_right(self.bounds, start)
        return not (i < len(self.bounds) and self.bounds[i] < end)

    def _ok_payload(self, ptype, pay):
        if ptype != 0:
            return True
        if len(pay) < 20 or (pay[0] >> 4) != 4:
            return False
        tot = (pay[2] << 8) | pay[3]
        return 20 <= tot <= len(pay)

    # -- header readings ----------------------------------------------------

    def _single(self, pos):
        """payload_configuration == 0.  M6's reading, kept verbatim."""
        b = self.s
        h = (b[pos] << 8) | b[pos + 1]
        ptype = h >> 13
        hm = (h >> 11) & 1
        ln = h & 0x7FF
        hdr = 2
        if hm:
            if pos + 2 >= len(b):
                return None
            ln |= b[pos + 2] << 11
            hdr = 3
            self.stats["hm1"] += 1
        if ptype == 4:
            # A/330 link-layer signalling: `length` counts the SIGNALLING
            # PAYLOAD only.  The 5-byte signalling sub-header (signalling_type
            # 8, signalling_type_extension 16, signalling_version 8,
            # signalling_format/id 8) sits between the base header and it.
            # Reading it as part of the payload leaves the walk 5 bytes short
            # exactly once per Link Mapping Table -- and the LMT repeats about
            # once a second, so at 8 Frames it costs nothing and at 476 Frames
            # it costs a resync (and an MFU) every single second.  Proven at
            # the byte level: the next ALP header lands at +5, not +0.
            hdr += 5
            self.stats["signalling"] += 1
        if ln == 0 or ln > self.MAX_ALP or pos + hdr + ln > len(b):
            self.stats["rej_len"] += 1
            return None
        if not self._anchor_clean(pos, pos + hdr + ln):
            self.stats["rej_anchor"] += 1
            return None
        pay = bytes(b[pos + hdr:pos + hdr + ln])
        if not self._ok_payload(ptype, pay):
            self.stats["rej_ip"] += 1
            return None
        return [(ptype, pay)], hdr + ln

    def _concat(self, pos):
        """payload_configuration == 1, segmentation_concatenation == 1.

        A/330 5.2:  length_MSB(4) count(3) HEF(1), then (count-1) 12-bit
        component_length fields, byte-aligned with 4 bits of stuffing when
        (count-1) is odd.  The LAST component's length is the remainder.
        """
        b = self.s
        if pos + 3 > len(b):
            return None
        ptype = b[pos] >> 5
        ln = ((b[pos] & 0x07) << 8) | b[pos + 1]
        msb, count, _hef = b[pos + 2] >> 4, (b[pos + 2] >> 1) & 0x07, \
            b[pos + 2] & 1
        total = (msb << 11) | ln
        if count < 2:
            return None
        nfield = count - 1
        nbytes = (nfield * 12 + 7) // 8
        hdr = 3 + nbytes
        if total == 0 or total > self.MAX_ALP or pos + hdr + total > len(b):
            return None
        if not self._anchor_clean(pos, pos + hdr + total):
            self.stats["rej_anchor"] += 1
            return None
        bits = int.from_bytes(bytes(b[pos + 3:pos + 3 + nbytes]), "big")
        pad = nbytes * 8 - nfield * 12
        bits >>= pad
        lens = [(bits >> (12 * (nfield - 1 - i))) & 0xFFF
                for i in range(nfield)]
        if sum(lens) >= total:
            return None
        lens.append(total - sum(lens))
        out, off = [], pos + hdr
        for L in lens:
            pay = bytes(b[off:off + L])
            if not self._ok_payload(ptype, pay):
                self.stats["rej_concat_ip"] += 1
                return None
            out.append((ptype, pay))
            off += L
        self.stats["concat_pkts"] += 1
        self.stats["concat_components"] += count
        return out, hdr + total

    def _segment(self, pos):
        """payload_configuration == 1, segmentation_concatenation == 0.

        A/330 5.2: segment_sequence_number(5) last_segment_indicator(1)
        SIF(1) HEF(1), then optional sub_stream_identifier(8).
        """
        b = self.s
        if pos + 3 > len(b):
            return None
        ptype = b[pos] >> 5
        ln = ((b[pos] & 0x07) << 8) | b[pos + 1]
        ssn = b[pos + 2] >> 3
        last = (b[pos + 2] >> 2) & 1
        sif = (b[pos + 2] >> 1) & 1
        hdr = 3 + (1 if sif else 0)
        if ln == 0 or ln > self.MAX_ALP or pos + hdr + ln > len(b):
            return None
        if not self._anchor_clean(pos, pos + hdr + ln):
            self.stats["rej_anchor"] += 1
            return None
        key = b[pos + 3] if sif else 0
        buf = self.segbuf.setdefault(key, {})
        buf[ssn] = bytes(b[pos + hdr:pos + hdr + ln])
        self.stats["seg_pkts"] += 1
        out = []
        if last:
            keys = sorted(buf)
            if keys == list(range(len(keys))):
                out = [(ptype, b"".join(buf[k] for k in keys))]
                self.stats["seg_reassembled"] += 1
            else:
                self.stats["seg_incomplete"] += 1
            self.segbuf[key] = {}
        return out, hdr + ln

    # -- the walk -----------------------------------------------------------

    def run(self):
        b, out = self.s, []
        pos = self.bounds[0] if self.bounds else 0
        while pos + 2 <= len(b):
            pc = (b[pos] >> 4) & 1
            if pc == 0:
                r = self._single(pos)
                self.stats["pc0"] += 1
            else:
                sc = (b[pos] >> 3) & 1
                r = self._concat(pos) if sc else self._segment(pos)
                self.stats["pc1"] += 1
            if r is None:
                nxt = [x for x in self.bounds if x > pos]
                if not nxt:
                    break
                pos = nxt[0]
                self.stats["resync"] += 1
                continue
            pkts, consumed = r
            out += pkts
            pos += consumed
        return out, self.stats


# ---------------------------------------------------------------------------
# IPv4 -- with fragment reassembly
# ---------------------------------------------------------------------------

class IpReasm:
    def __init__(self):
        self.frag = collections.defaultdict(dict)
        self.stats = collections.Counter()

    def feed(self, p):
        """-> list of complete (src, dst, sport, dport, udp_payload)."""
        if len(p) < 20 or (p[0] >> 4) != 4:
            self.stats["not_ip"] += 1
            return []
        ihl = (p[0] & 0xF) * 4
        tot = (p[2] << 8) | p[3]
        ident = (p[4] << 8) | p[5]
        flags = (p[6] >> 5) & 0x7
        foff = (((p[6] & 0x1F) << 8) | p[7]) * 8
        proto = p[9]
        src, dst = bytes(p[12:16]), bytes(p[16:20])
        if proto != 17:
            self.stats["not_udp"] += 1
            return []
        body = p[ihl:tot] if 0 < tot <= len(p) else p[ihl:]
        mf = (flags >> 0) & 1
        if not mf and foff == 0:
            return self._udp(src, dst, body)
        # fragment
        self.stats["frag"] += 1
        key = (src, dst, ident, proto)
        d = self.frag[key]
        d[foff] = bytes(body)
        if not mf:
            d["_end"] = foff + len(body)
        end = d.get("_end")
        if end is None:
            return []
        off, chunks = 0, []
        while off < end:
            if off not in d:
                return []
            chunks.append(d[off])
            off += len(d[off])
        del self.frag[key]
        self.stats["frag_reassembled"] += 1
        return self._udp(src, dst, b"".join(chunks))

    def _udp(self, src, dst, u):
        if len(u) < 8:
            return []
        ulen = (u[4] << 8) | u[5]
        pay = u[8:ulen] if 8 <= ulen <= len(u) else u[8:]
        self.stats["udp"] += 1
        return [(src, dst, (u[0] << 8) | u[1], (u[2] << 8) | u[3], pay)]


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------

def decode_chunk(path, rate, fmt, start_sec, n_frames, log):
    """Decode up to n_frames Frames starting at start_sec.  -> (stream, bounds,
    per-frame diagnostics)."""
    span = 0.05 + n_frames * C.FRAME_SAMPLES / 6.912e6
    y = C.load_frame(path, rate, fmt=fmt, start_sec=start_sec,
                     span_sec=span * (rate / 6.912e6))
    stream, bounds, frames = bytearray(), [], []
    t_prev = None
    for fi in range(n_frames):
        if (fi + 1) * C.FRAME_SAMPLES + 400000 > len(y):
            break
        centre = C.BOOTSTRAP if t_prev is None else t_prev + C.FRAME_SAMPLES
        t0, coh = C.fine_timing(y, span=20, centre=centre)
        t_prev = t0
        r16, pkts, diag = P6.decode_frame(y, t0, {})
        good = [p for p in pkts if p is not None]
        for p in good:
            ptr, pay = P6.bb_split(p)
            if ptr is not None:
                bounds.append(len(stream) + ptr)
            stream += pay
        frames.append(dict(t0=int(t0), coh=round(float(coh), 4),
                           conv=diag["converged"], bch=diag["bch_ok"],
                           plp16=bool(r16["converged"])))
    return stream, bounds, frames


_JOB = {}


def _job_init(path, rate, fmt):
    _JOB.update(path=path, rate=rate, fmt=fmt)


def _job(args):
    s0, n = args
    try:
        stream, bounds, fr = decode_chunk(_JOB["path"], _JOB["rate"],
                                          _JOB["fmt"], s0, n, print)
    except Exception as e:                                     # noqa: BLE001
        return s0, None, None, None, f"{type(e).__name__}: {e}"
    return s0, bytes(stream), bounds, fr, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("capture")
    ap.add_argument("--rate", type=float, required=True)
    ap.add_argument("--fmt")
    ap.add_argument("--frames", type=int, default=40)
    ap.add_argument("--chunk", type=int, default=16)
    ap.add_argument("--start", type=float, default=0.0)
    ap.add_argument("--out", default="m7.dg")
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--no-stitch", dest="stitch", action="store_false",
                    help="walk each chunk's Baseband stream separately "
                         "(loses one ALP packet per chunk boundary)")
    ap.add_argument("--bb", default=None,
                    help="baseband-stream cache; written if absent, READ if "
                         "present (skips the expensive FEC decode)")
    a = ap.parse_args()
    path = a.capture if os.path.isabs(a.capture) else os.path.join(HERE,
                                                                   a.capture)
    outp = a.out if os.path.isabs(a.out) else os.path.join(HERE, a.out)
    print("m7 -- RF33 Frames -> IP datagrams (boundaries preserved)")
    print("=" * 72)

    frame_s = C.FRAME_SAMPLES / 6.912e6
    done, chunks, allstats = 0, [], collections.Counter()
    ipr = IpReasm()
    flows = collections.Counter()
    fh = open(outp, "wb")
    t_start = time.time()

    plan = []
    k = 0
    while k < a.frames:
        n = min(a.chunk, a.frames - k)
        plan.append((a.start + k * frame_s, n))
        k += n

    bbp = None if not a.bb else (a.bb if os.path.isabs(a.bb)
                                 else os.path.join(HERE, a.bb))
    if bbp and os.path.exists(bbp):
        print(f"  reading baseband cache {bbp} "
              f"({os.path.getsize(bbp)/1e6:.1f} MB) -- no FEC decode")
        results = [(s0, stream, bounds, fr, None)
                   for s0, fr, stream, bounds in bb_read(bbp)]
    else:
        if a.jobs > 1:
            import concurrent.futures as cf
            for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                      "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
                os.environ.setdefault(v, "2")
            ex = cf.ProcessPoolExecutor(max_workers=a.jobs,
                                        initializer=_job_init,
                                        initargs=(path, a.rate, a.fmt))
            results = list(ex.map(_job, plan))
            ex.shutdown()
        else:
            _job_init(path, a.rate, a.fmt)
            results = [_job(p) for p in plan]
        if bbp:
            with open(bbp, "wb") as bf:
                for s0, stream, bounds, fr, err in results:
                    if err or not fr:
                        continue
                    bb_write(bf, s0, fr, stream, bounds)
            print(f"  wrote baseband cache {bbp} "
                  f"({os.path.getsize(bbp)/1e6:.1f} MB)")

    if a.stitch:
        # Chunk k covers Frames [k*chunk, (k+1)*chunk) and chunk k+1 starts on
        # the next Frame, so the Baseband Packet streams are CONTIGUOUS and
        # must be walked as one.  Walking them separately throws away one ALP
        # packet per chunk boundary -- invisible at 8 Frames, and at 476
        # Frames it was one 1400-byte hole in essentially EVERY MPU.
        big, bnd, frames = bytearray(), [], []
        for s0, stream, bounds, fr, err in results:
            if err or not fr:
                continue
            bnd += [b + len(big) for b in bounds]
            big += stream
            frames += fr
        results = [(a.start, bytes(big), bnd, frames, None)]
        print(f"  stitched {len(frames)} Frames into one "
              f"{len(big)/1e6:.2f} MB Baseband stream")

    for s0, stream, bounds, fr, err in results:
        if err:
            print(f"  chunk at {s0:.3f}s FAILED: {err}")
            continue
        if not fr:
            print(f"  capture exhausted at {s0:.3f}s")
            continue
        alps, st = AlpWalker(bytearray(stream), bounds).run()
        allstats += st
        nip = 0
        for t, p in alps:
            if t != 0:
                allstats[f"ptype{t}"] += 1
                continue
            for src, dst, sp, dp, pay in ipr.feed(p):
                fh.write(REC.pack(src, dst, sp, dp, len(pay)) + pay)
                flows[(".".join(map(str, dst)), dp)] += 1
                nip += 1
        conv = sum(f["conv"] for f in fr)
        bch = sum(f["bch"] for f in fr)
        cohs = [f["coh"] for f in fr]
        print(f"  chunk {len(chunks):3d}  t={s0:7.2f}s  {len(fr):3d} frames  "
              f"LDPC {conv}/{74*len(fr)}  BCH {bch}/{74*len(fr)}  "
              f"coh {min(cohs):.3f}-{max(cohs):.3f}  "
              f"ALP {len(alps)}  IP {nip}  resync {st['resync']}")
        chunks.append(dict(t=s0, frames=len(fr), conv=conv, bch=bch,
                           alp=len(alps), ip=nip, resync=st["resync"],
                           coh_min=min(cohs)))
        done += len(fr)
    fh.close()

    print(f"\n  {done} Frames in {time.time()-t_start:.0f}s -> {outp} "
          f"({os.path.getsize(outp)/1e6:.2f} MB)")
    print(f"\n  === ALP statistics (A/330 5.2) ===")
    for k in sorted(allstats):
        print(f"    {k:22s} {allstats[k]}")
    print(f"  IP: {dict(ipr.stats)}")
    print(f"\n  === flows ===")
    for (dst, dp), n in flows.most_common(30):
        tag = ""
        if (dst, dp) == (P6.LLS_ADDR, P6.LLS_PORT):
            tag = "   <== LLS (well-known)"
        print(f"    {dst:16s} :{dp:<6d} {n:6d} datagrams{tag}")

    js = os.path.splitext(outp)[0] + "_summary.json"
    with open(js, "w") as f:
        json.dump(dict(capture=os.path.basename(path), frames=done,
                       chunks=chunks, alp=dict(allstats),
                       ip=dict(ipr.stats),
                       flows={f"{k[0]}:{k[1]}": v for k, v in flows.items()}),
                  f, indent=1, default=float)
    print(f"\n  wrote {js}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
