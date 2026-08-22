#!/usr/bin/env python3
"""E67 decode gate: point OUR normal ROUTE/DASH assembly at all four RF33
video services from the SAME banked capture (lab/m13_out48) and see which
ones actually yield decodable HEVC.

WHUT (32.1) is the positive control -- it is signalled clear and we have
already decoded it (lab/plp1_out/ch132_1.mp4).  WTTG/WRC/WUSA are the
services under test.  Nothing here attempts to decrypt anything: we run the
ordinary pipeline and record what it produces.

Three independent readings per service:
  1. init segment structure   -- is the sample entry 'hvc1' (clear) or 'encv'
                                 (ISO/IEC 23001-7 Common Encryption)?
  2. NAL tiling of the mdat   -- do the 4-byte length prefixes tile the
                                 samples exactly, with legal HEVC NAL types?
  3. ffmpeg/ffprobe           -- does a real decoder produce frames?
"""
import json
import os
import re
import struct
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "m13_out48")
OUT = os.path.join(HERE, "e67_out")
os.makedirs(OUT, exist_ok=True)

SERVICES = [
    # (key, serviceId, name, channel, first media toi)
    ("239_255_32_1_8321", 1, "WHUT", "32.1", 225650589),
    ("239_255_5_1_8051", 3, "WTTG", "5.1", 225642448),
    ("239_255_4_1_8041", 4, "WRC", "4.1", 234326165),
    ("239_255_9_1_8091", 5, "WUSA", "9.1", 225642081),
]
# The FIRST object of every service is a partial (the capture joins the flow
# mid-object, so its leading LCT packets are missing and m7_route zero-fills).
# That is true symmetrically for all four services, so skip it; the LAST one is
# truncated by the end of capture, so stop before it.  Segments first+1..first+5.
SKIP_FIRST = 1
NSEG = 5
TSI = 10  # video

# ---------------------------------------------------------------- box walking


def boxes(buf, start, end, depth=0, path=""):
    """Yield (path, type, payload_start, payload_end)."""
    off = start
    while off + 8 <= end:
        size = struct.unpack(">I", buf[off:off + 4])[0]
        typ = buf[off + 4:off + 8].decode("latin1")
        hdr = 8
        if size == 1:
            size = struct.unpack(">Q", buf[off + 8:off + 16])[0]
            hdr = 16
        elif size == 0:
            size = end - off
        if size < hdr or off + size > end:
            return
        yield (path + "/" + typ, typ, off + hdr, off + size)
        off += size


def find_stsd_entries(buf):
    """Return the four-char sample entry codes inside stsd, plus any
    sinf/frma/schm details."""
    res = {"sample_entries": [], "frma": None, "schm": None, "tenc_kid": None}

    def rec(s, e, path):
        for p, typ, ps, pe in boxes(buf, s, e, path=path):
            if typ in ("moov", "trak", "mdia", "minf", "stbl"):
                rec(ps, pe, p)
            elif typ == "stsd":
                # 4 bytes version/flags + 4 bytes entry_count
                rec2(ps + 8, pe, p)

    def rec2(s, e, path):
        for p, typ, ps, pe in boxes(buf, s, e, path=path):
            res["sample_entries"].append(typ)
            # visual sample entry: 78 bytes of fixed fields before children
            body = ps + 78 if typ in ("encv", "hvc1", "hev1") else ps + 28
            for p2, t2, s2, e2 in boxes(buf, body, pe, path=p):
                if t2 == "sinf":
                    for p3, t3, s3, e3 in boxes(buf, s2, e2, path=p2):
                        if t3 == "frma":
                            res["frma"] = buf[s3:s3 + 4].decode("latin1")
                        elif t3 == "schm":
                            res["schm"] = buf[s3 + 4:s3 + 8].decode("latin1")
                        elif t3 == "schi":
                            for p4, t4, s4, e4 in boxes(buf, s3, e3, path=p3):
                                if t4 == "tenc":
                                    res["tenc_kid"] = buf[s4 + 8:s4 + 24].hex()

    rec(0, len(buf), "")
    return res


def trun_samples(buf):
    """Return (data_offset, [sample_sizes]) from the first traf/trun."""
    out = []

    def rec(s, e):
        for p, typ, ps, pe in boxes(buf, s, e):
            if typ in ("moof", "traf"):
                rec(ps, pe)
            elif typ == "tfhd":
                flags = struct.unpack(">I", buf[ps - 8 + 8:ps - 8 + 12])[0] & 0xFFFFFF
                out.append(("tfhd", flags))
            elif typ == "trun":
                vf = struct.unpack(">I", buf[ps:ps + 4])[0]
                flags = vf & 0xFFFFFF
                cnt = struct.unpack(">I", buf[ps + 4:ps + 8])[0]
                o = ps + 8
                doff = None
                if flags & 0x000001:
                    doff = struct.unpack(">i", buf[o:o + 4])[0]
                    o += 4
                if flags & 0x000004:
                    o += 4
                sizes = []
                for _ in range(cnt):
                    if flags & 0x000100:
                        o += 4
                    if flags & 0x000200:
                        sizes.append(struct.unpack(">I", buf[o:o + 4])[0])
                        o += 4
                    if flags & 0x000400:
                        o += 4
                    if flags & 0x000800:
                        o += 4
                out.append(("trun", doff, sizes))
    rec(0, len(buf))
    return out


HEVC_NAL_NAMES = {0: "TRAIL_N", 1: "TRAIL_R", 19: "IDR_W_RADL", 20: "IDR_N_LP",
                  21: "CRA_NUT", 32: "VPS", 33: "SPS", 34: "PPS", 35: "AUD",
                  39: "PREFIX_SEI", 40: "SUFFIX_SEI"}


def nal_tiling(seg):
    """Walk the media data as length-prefixed NAL units, sample by sample.

    trun@data_offset is measured from the first byte of the enclosing moof
    box, per ISO/IEC 14496-12.
    """
    tr = [t for t in trun_samples(seg) if t[0] == "trun"]
    mdat = None
    moof_start = None
    for p, typ, ps, pe in boxes(seg, 0, len(seg)):
        if typ == "mdat":
            mdat = (ps, pe)
        elif typ == "moof":
            moof_start = ps - 8
    if not tr or mdat is None or moof_start is None:
        return {"error": "no trun/moof/mdat"}
    _, doff, sizes = tr[0]
    ok_samples = 0
    bad_samples = 0
    nal_types = {}
    pos = moof_start + doff if doff is not None else mdat[0]
    for sz in sizes:
        s, e = pos, pos + sz
        pos = e
        if e > mdat[1]:
            bad_samples += 1
            continue
        o = s
        good = True
        local = []
        while o + 4 <= e:
            ln = struct.unpack(">I", seg[o:o + 4])[0]
            if ln == 0 or o + 4 + ln > e:
                good = False
                break
            nt = (seg[o + 4] >> 1) & 0x3F
            fbit = (seg[o + 4] >> 7) & 1
            if fbit != 0 or nt > 63:
                good = False
                break
            local.append(nt)
            o += 4 + ln
        if good and o == e and local:
            ok_samples += 1
            for nt in local:
                nal_types[nt] = nal_types.get(nt, 0) + 1
        else:
            bad_samples += 1
    return {
        "samples": len(sizes),
        "nal_tiling_ok": ok_samples,
        "nal_tiling_bad": bad_samples,
        "nal_types": {("%d:%s" % (k, HEVC_NAL_NAMES.get(k, "?"))): v
                      for k, v in sorted(nal_types.items())},
    }


def run(cmd):
    try:
        p = subprocess.run(cmd, capture_output=True, text=True,
                           timeout=180)  # pipe-ok: TimeoutExpired salvaged below
    except subprocess.TimeoutExpired as e:  # keep the evidence (7/31 balloon law)
        return -1, e.stdout or "", (e.stderr or "") + "\n[timeout 180s]"
    return p.returncode, p.stdout, p.stderr


def main():
    report = {}
    for key, sid, name, chan, first in SERVICES:
        rec = {"serviceId": sid, "name": name, "channel": chan}
        init_p = os.path.join(SRC, "obj_%s_tsi%d_toi4294967295.bin" % (key, TSI))
        init = open(init_p, "rb").read()
        rec["init"] = find_stsd_entries(init)
        rec["init_bytes"] = len(init)

        # assemble init + N media segments, exactly as our player does
        mp4 = os.path.join(OUT, "%s_video.mp4" % name.lower())
        blob = bytearray(init)
        segs = []
        for i in range(SKIP_FIRST, SKIP_FIRST + NSEG):
            p = os.path.join(SRC, "obj_%s_tsi%d_toi%d.bin" % (key, TSI, first + i))
            if os.path.exists(p):
                b = open(p, "rb").read()
                blob += b
                segs.append(b)
        open(mp4, "wb").write(bytes(blob))
        rec["assembled_bytes"] = len(blob)
        rec["segments"] = len(segs)

        # reading 2: NAL tiling on the first media segment
        rec["nal"] = nal_tiling(segs[0]) if segs else {"error": "none"}

        # reading 3: real decoder
        rc, so, se = run(["ffprobe", "-v", "error", "-print_format", "json",
                          "-show_streams", "-show_format", mp4])
        try:
            pr = json.loads(so) if so.strip() else {}
        except Exception:
            pr = {}
        rec["ffprobe_rc"] = rc
        rec["ffprobe_streams"] = [
            {k: v for k, v in s.items()
             if k in ("codec_name", "codec_tag_string", "width", "height",
                      "nb_frames", "profile")}
            for s in pr.get("streams", [])]
        rec["ffprobe_err"] = se.strip()[:300]

        png = os.path.join(OUT, "%s_frame.png" % name.lower())
        rc2, so2, se2 = run(["ffmpeg", "-y", "-v", "error", "-i", mp4,
                             "-frames:v", "1", png])
        rec["ffmpeg_frame_rc"] = rc2
        rec["ffmpeg_frame_err"] = se2.strip()[:400]
        rec["frame_png_bytes"] = os.path.getsize(png) if os.path.exists(png) else 0

        # decoded-frame count via null sink (the real quality metric)
        rc3, so3, se3 = run(["ffmpeg", "-y", "-nostdin", "-i", mp4,
                             "-f", "null", "-"])
        rec["ffmpeg_null_rc"] = rc3
        nframes = 0
        for ln in se3.splitlines():
            m = re.search(r"frame=\s*(\d+)", ln)
            if m:
                nframes = max(nframes, int(m.group(1)))
        rec["frames_decoded"] = nframes
        errs = [l for l in se3.splitlines()
                if "Error" in l or "error" in l or "Invalid" in l
                or "corrupt" in l or "Could not" in l]
        rec["decoder_errors"] = len(errs)
        rec["decoder_error_sample"] = errs[:4]

        report[name] = rec

    with open(os.path.join(HERE, "e67_decode_gate.json"), "w") as fh:
        json.dump(report, fh, indent=2)

    # verdict table
    hdr = ("%-6s %-6s %-6s %-6s %-11s %-16s %-8s %-8s %-7s %s" %
           ("svc", "ch", "entry", "schm", "NALtile", "ffprobe", "frames",
            "errs", "png", "verdict"))
    print(hdr)
    print("-" * len(hdr))
    for name, r in report.items():
        se_ = ",".join(r["init"]["sample_entries"])
        nal = r["nal"]
        tile = "%d/%d" % (nal.get("nal_tiling_ok", 0), nal.get("samples", 0))
        st = r["ffprobe_streams"]
        pv = ("%s %sx%s" % (st[0].get("codec_name"), st[0].get("width"),
                            st[0].get("height"))) if st else "-"
        good = (r["frames_decoded"] > 0 and r["decoder_errors"] == 0
                and r["frame_png_bytes"] > 10000
                and nal.get("samples", 0) > 0
                and nal.get("nal_tiling_ok", 0) == nal.get("samples", -1))
        print("%-6s %-6s %-6s %-6s %-11s %-16s %-8s %-8s %-7s %s" %
              (name, r["channel"], se_, str(r["init"]["schm"]), tile, pv,
               r["frames_decoded"], r["decoder_errors"], r["frame_png_bytes"],
               "DECODED" if good else "NO DECODE"))


if __name__ == "__main__":
    main()
