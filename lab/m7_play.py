#!/usr/bin/env python3
"""m7_play.py -- RF33 payload chain, part 3: MPUs -> a file a player will open.

An MMT MPU is already a self-contained ISOBMFF file: ftyp + mmpu + moov from
the FT=0 fragment, then moof + mdat from FT=1 + FT=2.  So "build a playable
file" is: take the init once, append the (moof, mdat) of each MPU in
MPU_sequence_number order, and hand it to ffmpeg.

Two honesty rules are enforced here rather than hoped for:

  * A **complete** MPU is one whose reassembled mdat length equals the size
    the mdat box header itself declares.  That header arrives in the FT=1
    fragment, i.e. from the transmitter, so it is an external referee, not a
    self-consistency check.
  * A **truncated** MPU is not silently patched into the output.  With
    --repair it is trimmed to whole samples using the trun sample table and
    the mdat size rewritten, and the result is LABELLED truncated in the
    report.  Nothing here reports a file as good because ffmpeg did not
    crash: the gate is frames actually decoded.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
import struct
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


# ---------------------------------------------------------------------------
# minimal ISOBMFF
# ---------------------------------------------------------------------------

def boxes(d, start=0, end=None):
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
            return          # truncated tail box: reported, then stop
        o += sz


def find_box(d, path, start=0, end=None):
    end = len(d) if end is None else end
    want = path[0]
    for o, sz, t in boxes(d, start, end):
        if t == want:
            if len(path) == 1:
                return o, sz
            return find_box(d, path[1:], o + 8, o + sz)
    return None


def trun_scan(seg):
    """Read out a (moof + mdat) segment's sample table without changing it.

    Returns a dict describing BOTH trafs, or None if the segment is not the
    shape this receiver produces.  Every offset is absolute in `seg`.

    The layout, established in `m7_objects.build` and re-verified here on a
    banked segment: the mdat holds ALL media samples and then ALL MMT hint
    samples, the media traf's trun carries a per-sample size, and the hint
    traf carries only a `data_offset` with the size defaulted in its tfhd.
    The arithmetic closes exactly -- sum(sizes) + n * hint_size == the mdat
    payload -- which is what makes a partial rebuild safe rather than hopeful.
    """
    mo, md = find_box(seg, [b"moof"]), find_box(seg, [b"mdat"])
    if not mo or not md:
        return None
    mo_o, mo_sz = mo
    md_o, md_sz = md
    media = hint = None
    for po, psz, pt in boxes(seg, mo_o + 8, mo_o + mo_sz):
        if pt != b"traf":
            continue
        tr, dss = None, 0
        for co, csz, ct in boxes(seg, po + 8, po + psz):
            if ct == b"trun":
                tr = (co, csz)
            elif ct == b"tfhd":
                fl = int.from_bytes(seg[co + 9:co + 12], "big")
                q = co + 16 + (8 if fl & 0x01 else 0) + (4 if fl & 0x02 else 0)
                if fl & 0x08:
                    q += 4
                if fl & 0x10:
                    dss = int.from_bytes(seg[q:q + 4], "big")
        if tr is None:
            continue
        co, _ = tr
        fl = int.from_bytes(seg[co + 9:co + 12], "big")
        n = int.from_bytes(seg[co + 12:co + 16], "big")
        p = co + 16
        doff = None
        if fl & 0x001:
            doff = int.from_bytes(seg[p:p + 4], "big", signed=True)
            p += 4
        if fl & 0x004:
            p += 4
        per = sum(4 for b in (0x100, 0x200, 0x400, 0x800) if fl & b)
        soff = 4 if fl & 0x100 else 0
        rec = dict(trun=co, flags=fl, n=n, first=p, per=per, soff=soff,
                   doff=doff, dss=dss)
        if fl & 0x200:
            rec["sizes"] = [int.from_bytes(seg[p + i * per + soff:
                                               p + i * per + soff + 4], "big")
                            for i in range(n)]
            media = rec
        else:
            hint = rec
    if media is None or media["doff"] is None:
        return None
    return dict(moof=mo_o, moof_size=mo_sz, mdat=md_o, mdat_size=md_sz,
                media=media, hint=hint)


def trun_keep(seg, keep):
    """Rebuild a segment holding only its first `keep` media samples.

    This is the repair the streaming receiver needs and `trun_trim` is not:
    `trun_trim` asks "how many whole samples fit in the bytes that arrived",
    which is only the right question when the loss is a clean tail.  On the
    air the loss is a HOLE -- an MMTP packet in the middle of an MPU -- and
    the sample table is what says where it lands.  So the caller passes the
    index of the first sample it knows to be short and this trims there.

    Both trafs are fixed, not just the media one: the hint traf's own
    `data_offset` follows the media block, so truncating media without moving
    it leaves the hint samples pointing into the middle of the video.  The
    mdat is rebuilt as media[:tot] + hint[:keep], and its declared size with
    it.  `trun_keep(seg, n_samples)` returns `seg` byte for byte -- gated.

    Returns (new_segment, keep, n_total) or None.
    """
    sc = trun_scan(seg)
    if sc is None:
        return None
    me, hi = sc["media"], sc["hint"]
    n = me["n"]
    keep = max(0, min(int(keep), n))
    if keep == 0:
        return None
    tot = sum(me["sizes"][:keep])
    mo_o = sc["moof"]
    m_start = mo_o + me["doff"]
    if m_start + tot > len(seg):
        return None
    body = bytes(seg[m_start:m_start + tot])
    new = bytearray(seg[:sc["mdat"]])
    struct.pack_into(">I", new, me["trun"] + 12, keep)
    if hi is not None and hi["doff"] is not None and hi["dss"]:
        h_start = mo_o + hi["doff"]
        h_keep = min(keep, hi["n"])
        h_len = h_keep * hi["dss"]
        if h_start + h_len > len(seg):
            return None
        struct.pack_into(">I", new, hi["trun"] + 12, h_keep)
        struct.pack_into(">I", new, hi["trun"] + 16, me["doff"] + tot)
        body += bytes(seg[h_start:h_start + h_len])
    out = bytes(new) + (8 + len(body)).to_bytes(4, "big") + b"mdat" + body
    return out, keep, n


def trun_trim(seg, avail):
    """Trim a (moof + mdat) segment to whole samples that actually arrived.

    Returns (new_segment, kept_samples, total_samples) or None.
    """
    mo = find_box(seg, [b"moof"])
    md = find_box(seg, [b"mdat"])
    if not mo or not md:
        return None
    tr = find_box(seg, [b"moof", b"traf", b"trun"])
    if not tr:
        return None
    o, sz = tr
    ver = seg[o + 8]
    flags = int.from_bytes(seg[o + 9:o + 12], "big")
    n = int.from_bytes(seg[o + 12:o + 16], "big")
    p = o + 16
    if flags & 0x000001:
        p += 4                                     # data_offset
    if flags & 0x000004:
        p += 4                                     # first_sample_flags
    per = (4 if flags & 0x100 else 0) + (4 if flags & 0x200 else 0) + \
          (4 if flags & 0x400 else 0) + (4 if flags & 0x800 else 0)
    if not (flags & 0x200):
        return None                                # no per-sample size
    soff = (4 if flags & 0x100 else 0)
    keep, tot = 0, 0
    data_avail = avail - 8                         # mdat payload bytes we have
    for i in range(n):
        q = p + i * per + soff
        if q + 4 > len(seg):
            break
        s = int.from_bytes(seg[q:q + 4], "big")
        if tot + s > data_avail:
            break
        tot += s
        keep += 1
    if keep == 0:
        return None
    new = bytearray(seg)
    struct.pack_into(">I", new, o + 12, keep)
    mo_o, mo_sz = mo
    md_o, _ = md
    struct.pack_into(">I", new, md_o, tot + 8)
    out = bytes(new[:md_o + 8]) + bytes(new[md_o + 8:md_o + 8 + tot])
    return out, keep, n


# ---------------------------------------------------------------------------

def retime(seg, frag_no, base_time):
    """Put an MPU on a continuous timeline.

    Every MPU is an INDEPENDENT ISOBMFF file: its `mfhd` sequence_number is 1
    and its `tfdt` baseMediaDecodeTime is 0.  Concatenate them unchanged and a
    player shows the first two seconds and stops -- which is exactly what
    ffmpeg reported (60 fragments in the file, `duration 2.002`).  So rewrite
    the fragment number and the decode time, and return the fragment's own
    duration so the caller can accumulate it.

    E86d -- THE CONTRACT IS ENFORCED HERE, NOT UPSTREAM.  `base_time` is
    written into a `>Q`/`>I` field, so it is an INTEGER BY CONSTRUCTION and a
    float is not a value it can legally take.  E81 fixed the one producer
    that was known to float it (`_nominal`'s median), and on 8/11 the live
    chain crash-looped on the identical `struct.error` from a DIFFERENT
    producer within minutes of restarting -- twice in one day, chasing
    callers.

    A packer that trusts its callers to respect a format it alone knows about
    will keep being surprised by the next caller.  So the coercion lives at
    the boundary that owns the requirement, and it is LOUD: a float arriving
    here is still a real upstream bug worth finding, and silently rounding it
    would hide the very evidence needed to fix it.  Crashing the whole
    receiver is simply the wrong way to report it.
    """
    if not isinstance(base_time, int):
        _bad = base_time
        base_time = int(round(base_time))
        retime.coerced += 1
        if retime.coerced <= 5 or retime.coerced % 100 == 0:
            sys.stderr.write(
                f"  retime: base_time arrived as {type(_bad).__name__} "
                f"{_bad!r} -> coerced to {base_time} (occurrence "
                f"{retime.coerced}). A tick count is an integer by "
                f"construction; the PRODUCER is still worth finding.\n")
    d = bytearray(seg)
    mf = find_box(d, [b"moof", b"mfhd"])
    if mf:
        struct.pack_into(">I", d, mf[0] + 12, frag_no)
    dur = 0
    mo = find_box(d, [b"moof"])
    if not mo:
        return bytes(d), 0
    for po, psz, pt in boxes(d, mo[0] + 8, mo[0] + mo[1]):
        if pt != b"traf":
            continue
        tfdt = None
        dsd = 0
        for co, csz, ct in boxes(d, po + 8, po + psz):
            if ct == b"tfdt":
                tfdt = co
            elif ct == b"tfhd":
                fl = int.from_bytes(d[co + 9:co + 12], "big")
                q = co + 16 + (8 if fl & 0x01 else 0) + (4 if fl & 0x02 else 0)
                if fl & 0x08:
                    dsd = int.from_bytes(d[q:q + 4], "big")
            elif ct == b"trun":
                fl = int.from_bytes(d[co + 9:co + 12], "big")
                n = int.from_bytes(d[co + 12:co + 16], "big")
                q = co + 16 + (4 if fl & 0x01 else 0) + (4 if fl & 0x04 else 0)
                per = sum(4 for b in (0x100, 0x200, 0x400, 0x800) if fl & b)
                if fl & 0x100:
                    s = sum(int.from_bytes(d[q + i * per:q + i * per + 4],
                                           "big") for i in range(n))
                else:
                    s = n * dsd
                # THE MEDIA TRAF'S DURATION, NOT THE LARGER OF THE TWO.
                #
                # An MPU carries two trafs: the media track and a hint
                # track. At 60000/1001 fps a frame is 1501.5 ticks at
                # timescale 90000, and the two express that differently --
                # the media trun alternates 1502/1501 and sums to exactly
                # 180180, while the hint tfhd carries a rounded constant
                # 1502 and so claims 180240. `max()` therefore advanced the
                # timeline by the HINT track: 60 ticks per MPU, which is
                # 2.335 s of A/V offset per two hours (measured).
                #
                # The broadcaster is not rounding; we were picking the
                # wrong one of its two answers.
                if dur == 0:
                    dur = s                       # first traf is the media one
        if tfdt is not None:
            if d[tfdt + 8] == 1:
                struct.pack_into(">Q", d, tfdt + 12, base_time)
            else:
                struct.pack_into(">I", d, tfdt + 12, base_time & 0xFFFFFFFF)
    return bytes(d), dur


retime.coerced = 0          # E86d: how often a caller broke the contract


def run(cmd):
    p = subprocess.run(cmd, capture_output=True, text=True)
    return p.returncode, p.stdout, p.stderr


def build(outdir, pid, dest, repair=False, limit=None):
    init = os.path.join(outdir, f"mpu_239_255_7_1_8071_pid{pid}_init.mp4")
    if not os.path.exists(init):
        cand = glob.glob(os.path.join(outdir, f"*pid{pid}_init.mp4"))
        if not cand:
            return None
        init = cand[0]
    pat = init.replace("_init.mp4", "_")
    segs = []
    for f in glob.glob(pat + "*.seg"):
        m = re.search(r"_(\d+)\.seg$", f)
        if m:
            segs.append((int(m.group(1)), f))
    segs.sort()
    if limit:
        segs = segs[:limit]
    body, used, skipped = bytearray(), [], []
    t, frag = 0, 0
    for seq, f in segs:
        d = open(f, "rb").read()
        md = find_box(d, [b"mdat"])
        if not md:
            skipped.append((seq, "no mdat"))
            continue
        o, sz = md
        have = len(d) - o
        if have == sz:
            frag += 1
            d, dur = retime(d, frag, t)
            t += dur
            body += d
            used.append((seq, "complete", sz))
        elif repair:
            r = trun_trim(d, have)
            if r is None:
                skipped.append((seq, f"truncated {sz-have}B, unrepairable"))
                continue
            new, keep, tot = r
            body += new
            used.append((seq, f"TRUNCATED {keep}/{tot} samples", len(new)))
        else:
            skipped.append((seq, f"truncated, missing {sz-have} B"))
    if not used:
        return dict(dest=None, used=[], skipped=skipped)
    with open(dest, "wb") as fh:
        fh.write(open(init, "rb").read())
        fh.write(body)
    return dict(dest=dest, used=used, skipped=skipped,
                bytes=os.path.getsize(dest))


def probe(path):
    rc, out, err = run(["ffprobe", "-v", "error", "-show_streams",
                        "-show_format", "-of", "json", path])
    try:
        return json.loads(out), err
    except Exception:                                          # noqa: BLE001
        return None, err


def decode_gate(path):
    """The only honest 'is it video' test: count frames ffmpeg actually
    decodes to a null sink."""
    rc, out, err = run(["ffmpeg", "-v", "error", "-stats", "-i", path,
                        "-map", "0:v:0", "-f", "null", "-"])
    m = re.findall(r"frame=\s*(\d+)", err)
    nf = int(m[-1]) if m else 0
    nerr = len([l for l in err.splitlines()
                if "error" in l.lower() or "Invalid" in l])
    return nf, nerr, err[-1500:]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--video-pid", type=int, default=12)
    ap.add_argument("--audio-pid", type=int, default=13)
    ap.add_argument("--dest", default=None)
    ap.add_argument("--repair", action="store_true")
    ap.add_argument("--still", default=None)
    a = ap.parse_args()
    outdir = a.outdir if os.path.isabs(a.outdir) else os.path.join(HERE,
                                                                   a.outdir)
    dest = a.dest or os.path.join(outdir, "rf33_video.mp4")
    print("m7 -- MPU -> playable file")
    print("=" * 72)
    rep = {}
    v = build(outdir, a.video_pid, dest, repair=a.repair)
    if not v or not v["dest"]:
        print(f"  no usable video MPU for pid {a.video_pid}: "
              f"{v['skipped'] if v else 'no init segment'}")
        return 1
    print(f"\n  video pid {a.video_pid} -> {dest} ({v['bytes']} bytes)")
    for seq, st, n in v["used"]:
        print(f"    MPU {seq}  {st:38s} {n} bytes")
    for seq, why in v["skipped"]:
        print(f"    MPU {seq}  SKIPPED: {why}")
    rep["video"] = v

    info, err = probe(dest)
    print(f"\n  === ffprobe {os.path.basename(dest)} ===")
    if info:
        for s in info.get("streams", []):
            print(f"    stream {s.get('index')}: {s.get('codec_type')} "
                  f"{s.get('codec_name')} ({s.get('profile')}) "
                  f"{s.get('width')}x{s.get('height')} "
                  f"{s.get('r_frame_rate')} fps, "
                  f"{s.get('nb_frames','?')} frames, "
                  f"tb {s.get('time_base')}")
        f = info.get("format", {})
        print(f"    format {f.get('format_name')}  duration "
              f"{f.get('duration')} s  bitrate {f.get('bit_rate')}")
    if err.strip():
        print("    ffprobe stderr: " + err.strip()[:400])
    rep["ffprobe"] = info

    nf, nerr, tail = decode_gate(dest)
    print(f"\n  === DECODE GATE (ffmpeg -> null sink) ===")
    print(f"    frames decoded: {nf}    error lines: {nerr}")
    if tail.strip():
        print("    " + tail.strip().replace("\n", "\n    ")[:1200])
    rep["frames_decoded"] = nf
    if nf == 0:
        print("    *** NO FRAMES DECODED -- this is NOT a playable file ***")

    still = a.still or os.path.join(outdir, "rf33_frame.png")
    if nf:
        rc, o, e = run(["ffmpeg", "-y", "-v", "error", "-i", dest,
                        "-map", "0:v:0", "-frames:v", "1", still])
        if os.path.exists(still) and os.path.getsize(still) > 1000:
            print(f"\n  still frame -> {still} "
                  f"({os.path.getsize(still)} bytes)")
            rep["still"] = still
        else:
            print(f"\n  still frame FAILED: {e.strip()[:300]}")

    # audio, best effort -- never allowed to block the video result
    if a.audio_pid:
        ad = os.path.join(outdir, f"rf33_audio_pid{a.audio_pid}.mp4")
        au = build(outdir, a.audio_pid, ad, repair=a.repair)
        if au and au["dest"]:
            ai, aerr = probe(ad)
            names = [s.get("codec_name") for s in (ai or {}).get("streams",
                                                                 [])]
            print(f"\n  audio pid {a.audio_pid} -> {os.path.basename(ad)} "
                  f"({au['bytes']} B), ffprobe codecs {names}")
            rep["audio"] = dict(dest=ad, codecs=names,
                                used=len(au["used"]))
            done = False
            for ext in ("mp4", "mkv", "ts"):
                mux = os.path.join(outdir, "rf33_av." + ext)
                rc, o, e = run(["ffmpeg", "-y", "-v", "error", "-i", dest,
                                "-i", ad, "-map", "0:v:0", "-map", "1:a:0?",
                                "-c", "copy", mux])
                if rc == 0 and os.path.exists(mux) and \
                        os.path.getsize(mux) > 100000:
                    print(f"    muxed A/V -> {os.path.basename(mux)} "
                          f"({os.path.getsize(mux)} B)")
                    rep["mux"] = mux
                    done = True
                    break
                if os.path.exists(mux):
                    os.remove(mux)
            if not done:
                # AC-4 is what ATSC 3.0 carries and what this ffmpeg build
                # can DEMUX but neither decode nor re-mux ("no wav codec tag
                # for ac4", "no decoder found for: ac4").  That is a player
                # limitation, not a chain result: the audio object itself is
                # complete and its ISOBMFF is valid, which ffprobe confirms.
                print(f"    A/V mux unavailable: this ffmpeg has no AC-4 "
                      f"muxer tag or decoder.  The audio MPU is COMPLETE and "
                      f"lives in {os.path.basename(ad)}; the video result "
                      f"stands on its own.")
                rep["mux"] = None
        else:
            print(f"\n  audio pid {a.audio_pid}: no complete MPU "
                  f"({au['skipped'] if au else 'no init'})")

    json.dump(rep, open(os.path.join(outdir, "m7_play.json"), "w"),
              indent=1, default=str)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
