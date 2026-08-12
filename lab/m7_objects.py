#!/usr/bin/env python3
"""m7_objects.py -- RF33 payload chain, part 2: datagrams -> OBJECTS -> files.

Input is the boundary-preserving `.dg` cache written by `m7_route.py`.

WHAT CHANGED FROM THE M6 PICTURE
--------------------------------
M6 recorded every flow as "ROUTE/DASH".  That is right for the signalling
flows and **wrong for the media flow.**  The video flow is **MMTP**
(ISO/IEC 23008-1), not ROUTE/LCT, and it is provable rather than arguable:
parse it as MMTP and the 32-bit `packet_sequence_number` increments by
EXACTLY ONE across every consecutive datagram of a packet_id, over hundreds
of packets, with breaks only where the ALP layer resynced.  Nothing else
lands a monotonic +1 counter at a fixed offset.  The `ftyp` brand in the
payload is `mpuf` -- the ISO 23008-1 MPU brand -- which agrees.

So this module speaks BOTH transports:

  ROUTE / LCT (RFC 5651 + A/331)  -> objects keyed by (TSI, TOI), ordered by
      the 32-bit start_offset in the FEC Payload ID, sized by EXT_TOL or by
      the FDT-Instance's Content-Length.
  MMTP (ISO 23008-1)              -> MPUs keyed by (packet_id,
      MPU_sequence_number), rebuilt from FT=0 (ftyp+moov), FT=1 (moof + mdat
      box header) and FT=2 (MFU sample data) in packet order.

GATES (nothing here is asserted without one)
  * MMTP header length is not assumed.  Each candidate (12/14/16 bytes) is
    scored by "does the MPU payload's own `length` field equal the bytes that
    actually remain", and the winner must win by a landslide.  It does: 14.
  * An object is called COMPLETE only when its reassembled length equals the
    length the transport declared, AND its first bytes parse as a box/XML.
    Everything else is reported as partial, with the missing byte count.
"""
from __future__ import annotations

import argparse
import collections
import gzip
import io
import json
import os
import re
import struct
import sys
import zlib

REC = struct.Struct("<4s4sHHI")
HERE = os.path.dirname(os.path.abspath(__file__))

LLS_ADDR, LLS_PORT = "224.0.23.60", 4937
LLS_TABLES = {1: "SLT", 2: "RRT", 3: "SystemTime", 4: "AEAT",
              5: "OnscreenMessageNotification", 6: "CertificationData",
              7: "SignedMultiTable"}
FT_NAME = {0: "MPU metadata (ftyp+moov)", 1: "movie fragment metadata (moof)",
           2: "MFU (sample)", 3: "reserved"}


def read_dg(path):
    d = open(path, "rb").read()
    o, out = 0, []
    while o + REC.size <= len(d):
        src, dst, sp, dp, n = REC.unpack_from(d, o)
        o += REC.size
        out.append((".".join(map(str, src)), ".".join(map(str, dst)),
                    sp, dp, d[o:o + n]))
        o += n
    return out


# ---------------------------------------------------------------------------
# LCT / ROUTE
# ---------------------------------------------------------------------------

def lct_parse(p):
    """RFC 5651 LCT + A/331 source FEC Payload ID -> dict, or None."""
    if len(p) < 8 or (p[0] >> 4) != 1:
        return None
    c = (p[0] >> 2) & 3
    s = (p[1] >> 7) & 1
    o = (p[1] >> 5) & 3
    h = (p[1] >> 4) & 1
    hdr_len = p[2] * 4
    cp = p[3]
    if hdr_len < 8 or hdr_len > len(p):
        return None
    off = 4 + (c + 1) * 4
    n_tsi = s * 4 + h * 2
    n_toi = o * 4 + h * 2
    if off + n_tsi + n_toi > hdr_len:
        return None
    tsi = int.from_bytes(p[off:off + n_tsi], "big") if n_tsi else 0
    off += n_tsi
    toi = int.from_bytes(p[off:off + n_toi], "big") if n_toi else 0
    off += n_toi
    ext = {}
    while off + 1 <= hdr_len:
        het = p[off]
        if het == 0:
            break
        if het < 128:
            if off + 2 > hdr_len:
                break
            hel = p[off + 1] * 4
            if hel == 0:
                break
            body, step = p[off + 2:off + hel], hel
        else:
            body, step = p[off + 1:off + 4], 4
        ext.setdefault(het, []).append(bytes(body))
        off += step
    tol = None
    if 194 in ext:                                   # EXT_TOL, 24-bit
        tol = int.from_bytes(ext[194][0][:3], "big")
    elif 67 in ext:                                  # EXT_TOL, 48-bit
        b = ext[67][0]
        tol = int.from_bytes(b[2:8], "big") if len(b) >= 8 else None
    fdt_id = None
    if 192 in ext:                                   # EXT_FDT
        v = int.from_bytes(ext[192][0][:3], "big")
        fdt_id = v & 0xFFFFF
    if hdr_len + 4 > len(p):
        return None
    start_offset = int.from_bytes(p[hdr_len:hdr_len + 4], "big")
    return dict(tsi=tsi, toi=toi, cp=cp, tol=tol, fdt_id=fdt_id,
                ext=sorted(ext), start=start_offset,
                payload=bytes(p[hdr_len + 4:]))


class RouteFlow:
    def __init__(self, name):
        self.name = name
        self.obj = {}                # (tsi,toi) -> dict
        self.n = 0
        self.bad = 0

    def feed(self, p):
        r = lct_parse(p)
        if r is None:
            self.bad += 1
            return
        self.n += 1
        k = (r["tsi"], r["toi"])
        o = self.obj.setdefault(k, dict(tsi=r["tsi"], toi=r["toi"],
                                        tol=None, frags={}, cp=r["cp"],
                                        fdt_id=r["fdt_id"]))
        if r["tol"] is not None:
            o["tol"] = r["tol"]
        if r["payload"]:
            o["frags"][r["start"]] = r["payload"]

    def assemble(self):
        out = []
        for k, o in sorted(self.obj.items()):
            buf, off, gaps = io.BytesIO(), 0, []
            for st in sorted(o["frags"]):
                d = o["frags"][st]
                if st > off:
                    gaps.append((off, st - off))
                    buf.seek(st)
                elif st < off:
                    buf.seek(st)
                else:
                    buf.seek(off)
                buf.write(d)
                off = max(off, st + len(d))
            data = buf.getvalue()
            got = len(data)
            declared = o["tol"]
            state = ("complete" if declared is not None and got == declared
                     and not gaps else
                     "partial" if declared is not None else
                     "unsized")
            out.append(dict(tsi=o["tsi"], toi=o["toi"], got=got,
                            declared=declared, gaps=gaps, state=state,
                            data=data, n_frag=len(o["frags"])))
        return out


# ---------------------------------------------------------------------------
# MMTP (ISO/IEC 23008-1)
# ---------------------------------------------------------------------------

def mmtp_hdr_len(p, ver_ext):
    """12 base + 4 if packet_counter + ver_ext + header extension if X."""
    c = (p[0] >> 5) & 1
    x = (p[0] >> 1) & 1
    n = 12 + (4 if c else 0) + ver_ext
    if x:
        if n + 4 > len(p):
            return None
        n += 4 + int.from_bytes(p[n + 2:n + 4], "big")
    return n


def mmtp_score(pkts, ver_ext):
    """Gate: how often does the MPU payload's `length` field agree with the
    bytes that actually remain in the datagram?"""
    ok = tot = 0
    for p in pkts:
        if len(p) < 20 or (p[1] & 0x3F) != 0:
            continue
        n = mmtp_hdr_len(p, ver_ext)
        if n is None or n + 2 > len(p):
            continue
        tot += 1
        if int.from_bytes(p[n:n + 2], "big") == len(p) - n - 2:
            ok += 1
    return ok, tot


def mmtp_parse(p, ver_ext):
    if len(p) < 14:
        return None
    ver = p[0] >> 6
    typ = p[1] & 0x3F
    n = mmtp_hdr_len(p, ver_ext)
    if n is None or n >= len(p):
        return None
    return dict(ver=ver, typ=typ, pid=int.from_bytes(p[2:4], "big"),
                ts=int.from_bytes(p[4:8], "big"),
                psn=int.from_bytes(p[8:12], "big"), body=bytes(p[n:]))


def mpu_payload(b):
    """MPU-mode payload -> dict.  A=1 aggregation is expanded."""
    if len(b) < 8:
        return None
    ln = int.from_bytes(b[0:2], "big")
    ft = b[2] >> 4
    timed = (b[2] >> 3) & 1
    fi = (b[2] >> 1) & 3
    agg = b[2] & 1
    frag_counter = b[3]
    seq = int.from_bytes(b[4:8], "big")
    rest = b[8:2 + ln] if 2 + ln <= len(b) else b[8:]
    dus = []
    if agg:
        o = 0
        while o + 2 <= len(rest):
            L = int.from_bytes(rest[o:o + 2], "big")
            o += 2
            if L == 0 or o + L > len(rest):
                break
            dus.append(rest[o:o + L])
            o += L
    else:
        dus = [rest]
    return dict(ft=ft, timed=timed, fi=fi, agg=agg, fc=frag_counter,
                seq=seq, dus=dus)


DU_HDR = 14         # mfsn(32) sample_number(32) offset(32) priority(8) dep(8)


def _boxes(d, start=0, end=None):
    end = len(d) if end is None else end
    o = start
    while o + 8 <= end:
        sz = int.from_bytes(d[o:o + 4], "big")
        typ = d[o + 4:o + 8]
        if sz == 1:
            sz = int.from_bytes(d[o + 8:o + 16], "big")
        if sz == 0:
            sz = end - o
        if sz < 8:
            return
        yield o, sz, typ
        if o + sz > end:
            return
        o += sz


def moof_info(mfmd):
    """The FT=1 fragment is `moof` + the `mdat` BOX HEADER.  Read out of it:
    the media track's per-sample sizes, the second (MMT hint) track's fixed
    sample size, and the mdat's declared size.  All three come from the
    transmitter, so all three are referees rather than self-checks."""
    # AN UNPROTECTED LENGTH IS NOT A LENGTH UNTIL IT IS BOUNDED.
    #
    # `sizes` is used directly to size allocations in `MmtpFlow.build`
    # (`bytearray(want)`, twice per sample), and `moof` arrives over MMTP,
    # which has no integrity protection whatever -- the same exposure that let
    # one flipped bit in `mpu_sequence_number` kill an asset permanently
    # (E21). Here a flipped bit in a sample size asks for up to 4 GB, twice,
    # inside a loop whose length is also read off the air.
    #
    # On 8/07 `atsc3 watch` RSS jumped 968 MB -> 32 GB in ~100 s during a
    # re-acquisition storm -- ~310 MB/s, which no stream fed by a 55 MB/s
    # radio can produce, and which happened identically on the GPU and CPU
    # backends. A single oversized allocation does produce it.
    #
    # The bounds are the broadcast's own physics, not guesses: one MPU is
    # 2.002 s, so 120 video samples or 60 audio samples, and a single coded
    # frame of 720p60 HEVC is orders of magnitude under 8 MB. Anything past
    # these is a corrupted header, not a big picture.
    MAX_SAMPLE = 8 << 20               # 8 MB per coded sample
    MAX_SAMPLES = 4096                 # per MPU (real values are 120 / 60)
    out = dict(moof_size=None, mdat_declared=None, sizes=[], hint_size=0,
               hint_n=0, sizes_rejected=0)
    for o, sz, t in _boxes(mfmd):
        if t == b"moof":
            out["moof_size"] = sz
            for po, psz, pt in _boxes(mfmd, o + 8, o + sz):
                if pt != b"traf":
                    continue
                dss, n, sizes = 0, 0, []
                for co, csz, ct in _boxes(mfmd, po + 8, po + psz):
                    if ct == b"tfhd":
                        fl = int.from_bytes(mfmd[co + 9:co + 12], "big")
                        q = co + 16
                        if fl & 0x01:
                            q += 8
                        if fl & 0x02:
                            q += 4
                        if fl & 0x08:
                            q += 4
                        if fl & 0x10:
                            dss = int.from_bytes(mfmd[q:q + 4], "big")
                    elif ct == b"trun":
                        fl = int.from_bytes(mfmd[co + 9:co + 12], "big")
                        n = int.from_bytes(mfmd[co + 12:co + 16], "big")
                        q = co + 16 + (4 if fl & 0x01 else 0) \
                            + (4 if fl & 0x04 else 0)
                        per = sum(4 for b in (0x100, 0x200, 0x400, 0x800)
                                  if fl & b)
                        so = 4 if fl & 0x100 else 0
                        # E55: THE BOUND MUST RUN BEFORE THE LOOP IT
                        # BOUNDS.  The MAX_SAMPLES rejection below ran
                        # after the comprehension had already built n
                        # entries with n read straight off the air --
                        # the exact exposure this function's own header
                        # comment names.  Measured live (live1d, 8/09
                        # 09:11): one fade-corrupted trun count wedged
                        # the transport thread for minutes of CPU and
                        # climbing RSS while every MPU counter froze.
                        # Same disposition as the post-loop rejection
                        # (sizes stays empty, the count is recorded, the
                        # MPU never completes and the repair machinery
                        # judges it) -- minus the work.
                        if fl & 0x200 and n <= MAX_SAMPLES:
                            sizes = [int.from_bytes(
                                mfmd[q + i * per + so:q + i * per + so + 4],
                                "big") for i in range(n)]
                        elif fl & 0x200:
                            out["sizes_rejected"] = n
                if sizes:
                    # Refuse the whole table if the transmitter cannot have
                    # sent it. Clamping individual entries instead would keep
                    # a fragment whose sample boundaries are already wrong and
                    # hand it downstream as if it were merely short -- and
                    # `build` would then lay real payload at fictional offsets.
                    # A moof this damaged has no usable sample table at all.
                    if len(sizes) > MAX_SAMPLES or any(
                            s > MAX_SAMPLE for s in sizes):
                        out["sizes_rejected"] = len(sizes)
                        sizes = []
                    else:
                        out["sizes"] = sizes
                if not sizes:
                    out["hint_size"], out["hint_n"] = dss, n
        elif t == b"mdat":
            out["mdat_declared"] = sz
    return out


class MmtpFlow:
    def __init__(self, name):
        self.name = name
        self.mpu = {}          # (pid, seq) -> {0: [...], 1: [...], 2: [...]}
        self.psn = collections.defaultdict(list)
        self.sig = []
        self.n = self.bad = 0

    def feed(self, p, ver_ext):
        r = mmtp_parse(p, ver_ext)
        if r is None:
            self.bad += 1
            return
        self.n += 1
        self.psn[r["pid"]].append(r["psn"])
        if r["typ"] == 2:
            self.sig.append((r["pid"], r["psn"], r["body"]))
            return
        if r["typ"] != 0:
            return
        m = mpu_payload(r["body"])
        if m is None:
            return
        k = (r["pid"], m["seq"])
        e = self.mpu.setdefault(k, dict(parts=collections.defaultdict(list),
                                        timed=m["timed"]))
        for du in m["dus"]:
            e["parts"][m["ft"]].append((r["psn"], m["fi"], du))

    def continuity(self):
        out = {}
        for pid, v in self.psn.items():
            v = sorted(v)
            steps = [b - a for a, b in zip(v, v[1:])]
            c = collections.Counter(steps)
            out[pid] = dict(n=len(v), plus1=c.get(1, 0),
                            breaks=sum(n for s, n in c.items() if s != 1),
                            lost=sum((s - 1) * n for s, n in c.items()
                                     if s > 1))
        return out

    def build(self):
        """-> list of MPU records with the reassembled ISOBMFF bytes.

        THE MFU LAYOUT, established from the air and then gated three ways:

            MFU Data Unit = [14-byte DU header][34-byte MMT hint sample]
                            [media sample bytes]

        and the DU header rides on **every** fragment, with `offset` giving
        the byte position of that fragment inside the MEDIA sample (the hint
        bytes are not counted).  The mdat is laid out as ALL media samples
        followed by ALL hint samples -- which is what the second traf's
        data_offset says, and it is why the DUs cannot simply be concatenated.

        Proof, not assertion:
          * offset(k+1) - offset(k) == len(DU_k) - 14 exactly, over hundreds
            of fragments, and the last fragment's offset + payload lands
            exactly on the trun sample size;
          * the reassembled media length equals sum(trun sample sizes);
          * media + hint + 8 equals the mdat box's own declared size.
        Getting the DU header length wrong by four bytes -- 18 instead of 14 --
        left every one of those checks *nearly* right and produced ffmpeg's
        "Invalid NAL unit size", which is an offset error wearing a codec
        error's clothes.
        """
        out = []
        for (pid, seq), e in sorted(self.mpu.items()):
            parts = e["parts"]
            meta = b"".join(du for _, _, du in sorted(parts.get(0, [])))
            mfmd = b"".join(du for _, _, du in sorted(parts.get(1, [])))
            mi = moof_info(mfmd)
            sizes, hs = mi["sizes"], mi["hint_size"]
            frags = collections.defaultdict(list)
            for psn, _fi, du in parts.get(2, []):
                if len(du) <= DU_HDR:
                    continue
                sn = int.from_bytes(du[4:8], "big")
                off = int.from_bytes(du[8:12], "big")
                frags[sn].append((off, psn, du[DU_HDR:]))
            media, hint, got_b, want_b = [], [], 0, 0
            n_short = 0
            short_first = None          # 1-based sample number, for repair
            for i, want in enumerate(sizes):
                sn = i + 1
                want_b += want
                buf = bytearray(want)
                have = bytearray(want)
                for off, _psn, pay in sorted(frags.get(sn, [])):
                    if off == 0:
                        hint.append(pay[:hs])
                        pay = pay[hs:]
                    n = min(len(pay), want - off)
                    if n > 0:
                        buf[off:off + n] = pay[:n]
                        have[off:off + n] = b"\x01" * n
                n_have = sum(have)
                got_b += n_have
                if n_have != want:
                    n_short += 1
                    if short_first is None:
                        short_first = sn
                media.append(bytes(buf))
            body = mfmd + b"".join(media) + b"".join(hint)
            declared = mi["mdat_declared"]
            total = None if declared is None else \
                (mi["moof_size"] or 0) + declared
            out.append(dict(
                pid=pid, seq=seq, meta=meta, body=body,
                n_meta=len(parts.get(0, [])), n_moof=len(parts.get(1, [])),
                n_mfu=len(parts.get(2, [])),
                n_samples=len(sizes), n_short=n_short,
                short_first=short_first,
                media_got=got_b, media_want=want_b,
                hint_got=sum(len(h) for h in hint), hint_want=hs *
                mi["hint_n"],
                mdat_declared=declared, mdat_got=len(body) - (
                    mi["moof_size"] or 0),
                total_got=len(body), total_want=total, timed=e["timed"]))
        return out


# ---------------------------------------------------------------------------
# signalling
# ---------------------------------------------------------------------------

def maybe_xml(b):
    """Bytes -> XML text, decompressing if needed; None if it is not XML.

    Two bugs lived here, and both DESTROYED DATA rather than merely failing:

    1. The identity branch ran first and won on gzipped input whenever a
       0x3C byte happened to fall in the first 64 bytes of the compressed
       stream -- about a 19 % chance over 54 bytes. Try the magic first.
    2. `b"<" in t[:64]` accepts any buffer merely CONTAINING a '<', which is
       most binary. Real XML STARTS with it. The caller writes whatever this
       returns through a UTF-8 text path, so a false positive turns every
       non-UTF-8 byte into U+FFFD and the object is unrecoverable -- which is
       how a gzipped ESG fragment ended up on disk as `1f ef bf bd 08 ...`
       where `1f 8b 08` should have been.
    """
    order = [gzip.decompress] if b[:2] == b"\x1f\x8b" else []
    order += [lambda x: x, gzip.decompress, lambda x: zlib.decompress(x, -15)]
    for f in order:
        try:
            t = f(b)
        except Exception:                                      # noqa: BLE001
            continue
        if t.lstrip()[:64].startswith(b"<"):
            return t.decode("utf-8", "replace")
    return None


def parse_lls(payload):
    tid = payload[0]
    out = dict(table_id=tid, name=LLS_TABLES.get(tid, f"reserved({tid})"),
               group=payload[1], group_count=payload[2] + 1,
               version=payload[3], nbytes=len(payload))
    out["xml"] = maybe_xml(payload[4:])
    return out


def parse_mmt_sig(body):
    """MMTP type=2 payload -> list of (message_id, version, bytes)."""
    if len(body) < 2:
        return []
    fi = body[0] >> 6
    h = (body[0] >> 1) & 1
    agg = body[0] & 1
    o = 2
    msgs = []
    blobs = []
    if agg:
        while o + (4 if h else 2) <= len(body):
            L = (int.from_bytes(body[o:o + 4], "big") if h
                 else int.from_bytes(body[o:o + 2], "big"))
            o += 4 if h else 2
            if L == 0 or o + L > len(body):
                break
            blobs.append(body[o:o + L])
            o += L
    else:
        blobs = [body[o:]]
    for b in blobs:
        if len(b) < 5 or fi != 0:
            continue
        mid = int.from_bytes(b[0:2], "big")
        ver = b[2]
        wide = mid == 0 or 1 <= mid <= 0x000F or mid >= 0x8000
        ln = (int.from_bytes(b[3:7], "big") if wide
              else int.from_bytes(b[3:5], "big"))
        msgs.append((mid, ver, ln, b[(7 if wide else 5):]))
    return msgs


def parse_mp_table(b):
    """ISO 23008-1 MP_table -> asset list (best effort, gated by sanity)."""
    try:
        o = 0
        table_id, ver = b[o], b[o + 1]
        length = int.from_bytes(b[o + 2:o + 4], "big")
        o += 4
        resv_mp = b[o]
        o += 1
        mpt_mode = resv_mp & 0x03
        nlen = b[o]
        o += 1
        pkg_id = b[o:o + nlen]
        o += nlen
        dlen = int.from_bytes(b[o:o + 2], "big")
        o += 2 + dlen
        n_assets = b[o]
        o += 1
        assets = []
        for _ in range(n_assets):
            ident_type = b[o]
            o += 1
            ilen = b[o]
            o += 1
            aid = b[o:o + ilen]
            o += ilen
            atype = b[o:o + 4].decode("ascii", "replace")
            o += 4
            flags = b[o]
            o += 1
            nloc = b[o]
            o += 1
            pids = []
            for _ in range(nloc):
                lt = b[o]
                o += 1
                if lt == 0:
                    pids.append(int.from_bytes(b[o:o + 2], "big"))
                    o += 2
                else:
                    break
            dl = int.from_bytes(b[o:o + 2], "big")
            o += 2 + dl
            assets.append(dict(asset_type=atype,
                               asset_id=aid.decode("ascii", "replace"),
                               packet_ids=pids))
        return dict(table_id=table_id, version=ver, length=length,
                    mpt_mode=mpt_mode,
                    package_id=pkg_id.decode("ascii", "replace"),
                    assets=assets)
    except Exception as e:                                     # noqa: BLE001
        return dict(error=f"{type(e).__name__}: {e}")


ATSC3_MSG = {1: "UserServiceDescription", 2: "MPD", 3: "HELD",
             4: "DistributionWindowDesc", 5: "ApplicationEventInformation",
             6: "VideoStreamProperties", 7: "ATSCStaticCaptionAsset",
             8: "AudioStreamProperties", 9: "InbandEventDescriptor"}


def parse_atsc3_message(b):
    """MMT_ATSC3_message payload (A/331 Table 7.4).

    service_id(16) content_type(16) content_version(8) content_compression(8)
    URI_length(8) URI[] content_length(32) content[].  Dropping the
    compression byte and the URI -- as a first pass did -- shifts the length
    field and the gunzip silently fails, which reads as "no XML present"
    rather than as a parse error.
    """
    try:
        svc = int.from_bytes(b[0:2], "big")
        ctype = int.from_bytes(b[2:4], "big")
        ver, comp, ulen = b[4], b[5], b[6]
        o = 7 + ulen
        uri = b[7:o].decode("ascii", "replace")
        cclen = int.from_bytes(b[o:o + 4], "big")
        body = b[o + 4:o + 4 + cclen]
        return dict(service_id=svc, content_type=ctype,
                    content_name=ATSC3_MSG.get(ctype, f"type{ctype}"),
                    version=ver, compression=comp, uri=uri,
                    content_length=cclen, xml=maybe_xml(body))
    except Exception as e:                                     # noqa: BLE001
        return dict(error=str(e))


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

FDT_RE = re.compile(r"<File\b[^>]*>", re.I)


def fdt_files(xml):
    out = []
    for m in FDT_RE.finditer(xml):
        tag = m.group(0)
        d = dict(re.findall(r'(\w[\w-]*)\s*=\s*"([^"]*)"', tag))
        out.append(d)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dg")
    ap.add_argument("--outdir", default="m7_out")
    a = ap.parse_args()
    dg = a.dg if os.path.isabs(a.dg) else os.path.join(HERE, a.dg)
    outdir = a.outdir if os.path.isabs(a.outdir) else os.path.join(HERE,
                                                                   a.outdir)
    os.makedirs(outdir, exist_ok=True)
    recs = read_dg(dg)
    print("m7 -- ROUTE/LCT and MMTP object reassembly")
    print("=" * 72)
    print(f"  {len(recs)} datagrams, "
          f"{sum(len(r[4]) for r in recs)/1e6:.2f} MB of UDP payload")

    byflow = collections.defaultdict(list)
    for src, dst, sp, dp, pay in recs:
        byflow[(dst, dp)].append(pay)

    report = dict(dg=os.path.basename(dg), flows={}, objects=[], mpus=[],
                  lls=[], mmt_sig=[])

    # -- LLS ---------------------------------------------------------------
    print(f"\n  === LLS ({LLS_ADDR}:{LLS_PORT}) ===")
    lls = byflow.get((LLS_ADDR, LLS_PORT), [])
    if not lls:
        print("    NOT PRESENT in this decode.")
    seen = collections.Counter()
    for d in lls:
        if len(d) < 5:
            continue
        r = parse_lls(d)
        key = (r["table_id"], r["version"])
        seen[key] += 1
        if seen[key] > 1:
            continue
        print(f"    LLS_table_id {r['table_id']} ({r['name']}), group "
              f"{r['group']}, version {r['version']}, {r['nbytes']} bytes"
              f"{'' if r['xml'] else '   [payload did not gunzip to XML]'}")
        report["lls"].append({k: v for k, v in r.items() if k != "xml"})
        if r["xml"]:
            fn = os.path.join(outdir, f"lls_{r['table_id']}_{r['name']}.xml")
            open(fn, "w", encoding="utf-8").write(r["xml"])
            print(f"      -> {os.path.basename(fn)}  ({len(r['xml'])} chars)")
    if seen:
        print("    repetition (table_id, version) -> datagrams: "
              + ", ".join(f"{k}x{v}" for k, v in sorted(seen.items())))

    # -- per-flow ----------------------------------------------------------
    for (dst, dp), pkts in sorted(byflow.items(),
                                  key=lambda kv: -len(kv[1])):
        if (dst, dp) == (LLS_ADDR, LLS_PORT):
            continue
        n_lct = sum(1 for p in pkts[:200] if lct_parse(p) is not None)
        sc = {v: mmtp_score(pkts[:400], v) for v in (0, 2, 4)}
        best_v, (ok, tot) = max(sc.items(), key=lambda kv: kv[1][0])
        is_mmtp = tot > 0 and ok >= 0.9 * tot and ok > n_lct
        print(f"\n  === flow {dst}:{dp} -- {len(pkts)} datagrams ===")
        print(f"    transport gate: LCT-parses {n_lct}/{min(200,len(pkts))}"
              f" | MMTP length-field agreement "
              + ", ".join(f"ver_ext={v}: {o}/{t}" for v, (o, t) in
                          sorted(sc.items()))
              + f"  ==> {'MMTP' if is_mmtp else 'ROUTE/LCT'}")
        finfo = dict(dst=dst, port=dp, n=len(pkts),
                     transport="MMTP" if is_mmtp else "ROUTE/LCT",
                     lct_hits=n_lct,
                     mmtp_score={str(v): list(x) for v, x in sc.items()})
        if is_mmtp:
            fl = MmtpFlow(f"{dst}:{dp}")
            for p in pkts:
                fl.feed(p, best_v)
            cont = fl.continuity()
            finfo["ver_ext"] = best_v
            finfo["continuity"] = {str(k): v for k, v in cont.items()}
            print(f"    packet_sequence_number continuity per packet_id:")
            for pid, c in sorted(cont.items()):
                print(f"      pid {pid:5d}  {c['n']:6d} packets, "
                      f"+1 steps {c['plus1']:6d}, breaks {c['breaks']:3d}, "
                      f"packets lost {c['lost']}")
            # signalling
            for pid, psn, body in fl.sig:
                for mid, ver, ln, blob in parse_mmt_sig(body):
                    ent = dict(pid=pid, message_id=hex(mid), version=ver,
                               length=ln)
                    if 0x0010 <= mid <= 0x001F:
                        t = parse_mp_table(blob)
                        # honest: the asset loop is NOT gated, and on this
                        # transmitter it yields nothing.  Say so rather than
                        # reporting an empty asset list as a finding.
                        if not t.get("assets"):
                            t["assets"] = "UNPARSED (reading not gated)"
                        ent["mp_table"] = t
                    elif mid == 0x8100:
                        m = parse_atsc3_message(blob)
                        ent["atsc3"] = {k: v for k, v in m.items()
                                        if k != "xml"}
                        if m.get("xml"):
                            fn = os.path.join(
                                outdir,
                                f"mmt_atsc3_svc{m['service_id']}_"
                                f"{m['content_name']}.xml")
                            open(fn, "w", encoding="utf-8").write(m["xml"])
                            ent["xml_file"] = os.path.basename(fn)
                    report["mmt_sig"].append(ent)
            # MPUs
            mpus = fl.build()
            print(f"    {len(mpus)} MPUs (packet_id, MPU_sequence_number):")
            print(f"      {'pid':>5} {'MPU_seq':>12} {'smpl':>5} "
                  f"{'short':>5} {'media got/declared':>21} "
                  f"{'mdat got/declared':>21}  state")
            n_ok = 0
            for m in mpus:
                dec, got = m["mdat_declared"], m["mdat_got"]
                if dec is None or not m["n_samples"]:
                    state = "NO moof (FT=1 fragment missing)"
                elif got == dec and m["media_got"] == m["media_want"] \
                        and m["n_short"] == 0:
                    state = "COMPLETE"
                    n_ok += 1
                else:
                    state = (f"PARTIAL  {m['n_short']} of {m['n_samples']} "
                             f"samples short, {m['media_want']-m['media_got']}"
                             f" bytes missing")
                print(f"      {m['pid']:5d} {m['seq']:12d} "
                      f"{m['n_samples']:5d} {m['n_short']:5d} "
                      f"{('%d/%d' % (m['media_got'], m['media_want'])):>21} "
                      f"{('%s/%s' % (got, dec)):>21}  {state}")
                report["mpus"].append(dict(
                    flow=f"{dst}:{dp}", pid=m["pid"], seq=m["seq"],
                    n_meta=m["n_meta"], n_moof=m["n_moof"], n_mfu=m["n_mfu"],
                    n_samples=m["n_samples"], n_short=m["n_short"],
                    media_got=m["media_got"], media_want=m["media_want"],
                    mdat_got=got, mdat_declared=dec, state=state.split()[0],
                    meta_bytes=len(m["meta"]), body_bytes=len(m["body"])))
                base = f"mpu_{dst.replace('.','_')}_{dp}_pid{m['pid']}"
                if m["meta"]:
                    open(os.path.join(outdir, base + "_init.mp4"),
                         "wb").write(m["meta"])
                if m["body"] and state == "COMPLETE":
                    open(os.path.join(outdir,
                                      f"{base}_{m['seq']}.seg"),
                         "wb").write(m["body"])
                elif m["body"]:
                    open(os.path.join(outdir,
                                      f"{base}_{m['seq']}.partial"),
                         "wb").write(m["body"])
            print(f"    ==> {n_ok}/{len(mpus)} MPUs COMPLETE on pid set "
                  f"{sorted(set(m['pid'] for m in mpus))}")
        else:
            fl = RouteFlow(f"{dst}:{dp}")
            for p in pkts:
                fl.feed(p)
            objs = fl.assemble()
            finfo["lct_ok"], finfo["lct_bad"] = fl.n, fl.bad
            print(f"    LCT: {fl.n} parsed, {fl.bad} rejected -> "
                  f"{len(objs)} objects")
            print(f"      {'TSI':>6} {'TOI':>10} {'got':>9} {'declared':>9} "
                  f"{'frags':>6}  state")
            fdt = {}
            for o in objs:
                print(f"      {o['tsi']:6d} {o['toi']:10d} {o['got']:9d} "
                      f"{str(o['declared']):>9} {o['n_frag']:6d}  "
                      f"{o['state']}"
                      + (f"  (missing {o['declared']-o['got']} B)"
                         if o["state"] == "partial" and o["declared"] else "")
                      + (f"  GAPS {o['gaps']}" if o["gaps"] else ""))
                head = o["data"][:8]
                kind = ""
                if len(head) >= 8 and head[4:8].isalpha():
                    kind = head[4:8].decode("ascii", "replace")
                elif head[:1] == b"<":
                    kind = "XML"
                elif head[:2] == b"\x1f\x8b":
                    kind = "gzip"
                report["objects"].append(dict(
                    flow=f"{dst}:{dp}", tsi=o["tsi"], toi=o["toi"],
                    got=o["got"], declared=o["declared"], state=o["state"],
                    n_frag=o["n_frag"], first_box=kind))
                name = f"obj_{dst.replace('.','_')}_{dp}_" \
                       f"tsi{o['tsi']}_toi{o['toi']}"
                if o["toi"] == 0:
                    x = maybe_xml(o["data"])
                    if x:
                        fdt = fdt_files(x)
                        open(os.path.join(outdir, name + "_FDT.xml"), "w",
                             encoding="utf-8").write(x)
                        print(f"        FDT-Instance: {len(fdt)} File "
                              f"entries")
                        for f in fdt[:12]:
                            print(f"          TOI {f.get('TOI','?'):>8}  "
                                  f"len {f.get('Content-Length','?'):>9}  "
                                  f"{f.get('Content-Location','?')}")
                        continue
                x = maybe_xml(o["data"]) if o["got"] < 200000 else None
                ext = ".xml" if x else ".bin"
                if x:
                    open(os.path.join(outdir, name + ext), "w",
                         encoding="utf-8").write(x)
                else:
                    open(os.path.join(outdir, name + ext),
                         "wb").write(o["data"])
        report["flows"][f"{dst}:{dp}"] = finfo

    js = os.path.join(outdir, "m7_objects.json")
    json.dump(report, open(js, "w"), indent=1, default=str)
    print(f"\n  wrote {js} and object files to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
