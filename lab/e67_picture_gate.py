#!/usr/bin/env python3
"""E67 picture gate + senc quantification.

The naive decode gate said all four RF33 video services "decoded": ffprobe
parsed them, ffmpeg emitted 600 frames and logged zero errors.  Looking at the
frames killed that verdict -- WHUT is a photograph of a man in a study, WUSA is
a green field with a band of noise.  That is the campaign's own law again: a
parse closing is never sufficient.  So gate the PICTURE, not the parse.

Two measurements per service:

  A. senc / subsample accounting.  ISO/IEC 23001-7 Common Encryption with
     subsample ranges leaves the NAL length prefixes and slice headers CLEAR
     and encrypts the residual.  That is exactly why the NAL-tiling check
     passed 120/120 on encrypted services -- it is not a discriminator.  The
     senc box tells us directly how many payload bytes are protected.

  B. Picture content.  Decode frames to raw luma and measure
       flat_frac  -- fraction of pixels in the single most common value
                     (decoder concealment paints large uniform areas)
       temporal_r -- Pearson r between consecutive frames.  Real video is
                     strongly correlated frame to frame; garbage is not.
"""
import json
import os
import struct
import subprocess

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "e67_out")
SRC = os.path.join(HERE, "m13_out48")

SERVICES = [
    ("239_255_32_1_8321", 1, "WHUT", "32.1", 225650589),
    ("239_255_5_1_8051", 3, "WTTG", "5.1", 225642448),
    ("239_255_4_1_8041", 4, "WRC", "4.1", 234326165),
    ("239_255_9_1_8091", 5, "WUSA", "9.1", 225642081),
]
W, H = 1920, 1080
NFRAMES = 24


def boxes(buf, start, end):
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
        yield (typ, off + hdr, off + size)
        off += size


def parse_tenc(init):
    """Pull default_Per_Sample_IV_Size / isProtected / KID out of tenc.

    tenc payload: [0:4] version+flags, [4] reserved, [5] crypt/skip byte block,
    [6] default_isProtected, [7] default_Per_Sample_IV_Size, [8:24] default_KID.
    """
    i = init.find(b"tenc")
    if i < 0:
        return None
    p = i + 4
    return {
        "isProtected": init[p + 6],
        "iv_size": init[p + 7],
        "default_KID": init[p + 8:p + 24].hex(),
    }


def parse_senc(seg, iv_size=16):
    """Return per-sample IVs and subsample (clear, protected) byte counts."""
    loc = None

    def rec(s, e):
        nonlocal loc
        for typ, ps, pe in boxes(seg, s, e):
            if typ in ("moof", "traf"):
                rec(ps, pe)
            elif typ == "senc":
                loc = (ps, pe)
    rec(0, len(seg))
    if loc is None:
        return None
    ps, pe = loc
    vf = struct.unpack(">I", seg[ps:ps + 4])[0]
    flags = vf & 0xFFFFFF
    cnt = struct.unpack(">I", seg[ps + 4:ps + 8])[0]
    o = ps + 8
    ivs, clear_tot, prot_tot, nsub = [], 0, 0, 0
    for _ in range(cnt):
        if o + iv_size > pe:
            break
        ivs.append(seg[o:o + iv_size])
        o += iv_size
        if flags & 0x000002:  # subsample encryption present
            if o + 2 > pe:
                break
            n = struct.unpack(">H", seg[o:o + 2])[0]
            o += 2
            for _ in range(n):
                if o + 6 > pe:
                    break
                c = struct.unpack(">H", seg[o:o + 2])[0]
                p = struct.unpack(">I", seg[o + 2:o + 6])[0]
                clear_tot += c
                prot_tot += p
                nsub += 1
                o += 6
    nonzero = sum(1 for iv in ivs if any(iv))
    uniq = len(set(ivs))
    return {
        "senc_flags": flags,
        "subsample_encryption": bool(flags & 0x2),
        "sample_count": cnt,
        "ivs_parsed": len(ivs),
        "ivs_nonzero": nonzero,
        "ivs_unique": uniq,
        "subsamples": nsub,
        "clear_bytes": clear_tot,
        "protected_bytes": prot_tot,
        "protected_frac": (round(prot_tot / (clear_tot + prot_tot), 4)
                           if (clear_tot + prot_tot) else None),
    }


def picture_stats(mp4):
    cmd = ["ffmpeg", "-v", "error", "-i", mp4, "-frames:v", str(NFRAMES),
           "-pix_fmt", "gray", "-f", "rawvideo", "-"]
    p = subprocess.run(cmd, capture_output=True, timeout=300)
    raw = p.stdout
    n = len(raw) // (W * H)
    if n < 2:
        return {"frames": n, "error": "too few frames"}
    a = np.frombuffer(raw[:n * W * H], dtype=np.uint8).reshape(n, H, W)
    flats, rs = [], []
    for i in range(n):
        f = a[i]
        counts = np.bincount(f.ravel(), minlength=256)
        flats.append(counts.max() / f.size)
        if i:
            x = a[i - 1].ravel().astype(np.float32)
            y = f.ravel().astype(np.float32)
            sx, sy = x.std(), y.std()
            rs.append(float(((x - x.mean()) * (y - y.mean())).mean() / (sx * sy))
                      if sx > 1e-6 and sy > 1e-6 else 0.0)
    return {
        "frames": n,
        "flat_frac_mean": round(float(np.mean(flats)), 4),
        "flat_frac_max": round(float(np.max(flats)), 4),
        "temporal_r_mean": round(float(np.mean(rs)), 4),
        "temporal_r_min": round(float(np.min(rs)), 4),
        "luma_mean": round(float(a.mean()), 2),
        "luma_std": round(float(a.std()), 2),
    }


def main():
    rep = {}
    for key, sid, name, chan, first in SERVICES:
        seg = open(os.path.join(
            SRC, "obj_%s_tsi10_toi%d.bin" % (key, first + 1)), "rb").read()
        init = open(os.path.join(
            SRC, "obj_%s_tsi10_toi4294967295.bin" % key), "rb").read()
        r = {"serviceId": sid, "channel": chan}
        r["tenc"] = parse_tenc(init)
        ivs = r["tenc"]["iv_size"] if r["tenc"] else 16
        r["senc"] = parse_senc(seg, iv_size=ivs)
        r["picture"] = picture_stats(os.path.join(OUT, "%s_video.mp4" % name.lower()))
        rep[name] = r
    with open(os.path.join(HERE, "e67_picture_gate.json"), "w") as fh:
        json.dump(rep, fh, indent=2)

    hdr = ("%-6s %-6s %-9s %-9s %-14s %-11s %-11s %s" %
           ("svc", "ch", "senc", "prot%", "IVs uniq/nz", "flat_frac",
            "temporal_r", "PICTURE"))
    print(hdr)
    print("-" * len(hdr))
    for name, r in rep.items():
        s = r["senc"]
        p = r["picture"]
        if s:
            senc = "present"
            pf = "%.1f%%" % (100 * s["protected_frac"])
            iv = "%d/%d" % (s["ivs_unique"], s["ivs_nonzero"])
        else:
            senc, pf, iv = "absent", "0.0%", "-"
        # a real picture: little concealment, high frame-to-frame correlation
        ok = p.get("flat_frac_mean", 1) < 0.10 and p.get("temporal_r_mean", 0) > 0.80
        print("%-6s %-6s %-9s %-9s %-14s %-11s %-11s %s" %
              (name, r["channel"], senc, pf, iv,
               p.get("flat_frac_mean"), p.get("temporal_r_mean"),
               "REAL VIDEO" if ok else "GARBAGE / NO PICTURE"))


if __name__ == "__main__":
    main()
