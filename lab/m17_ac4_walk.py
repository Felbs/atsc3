#!/usr/bin/env python3
"""M17 -- can we READ the AC-4 bitstream?  The gate before dolbyTuna starts.

Writing an AC-4 decoder is months of work.  Before spending any of it, one
question has to be answered with evidence: **can we walk the bitstream at all,
or are we guessing at bit offsets?**

THE TEST HAD TO BE ABLE TO FAIL
--------------------------------
"Each sample is one AC-4 frame and its size is in the `trun`" is true BY
CONSTRUCTION -- checking it proves nothing about our understanding of AC-4.  So
this gates on something internal to the bitstream instead, and picks the field
that is hardest to fake:

    ac4_toc() {
        bitstream_version   2 bits    (escape to variable_bits if == 3)
        sequence_counter   10 bits
        ...

**`sequence_counter` must increment by exactly 1, modulo 1024, from one frame
to the next, across thousands of frames.**  If our bit offsets are wrong by
even one bit, that field is garbage and the ramp collapses immediately.  A
clean ramp over 6,960 frames cannot happen by accident: there are 1024 possible
values and a wrong offset would have to reproduce a +1 sequence 6,959 times.

And `bitstream_version` read from every frame's TOC must equal the value the
`dac4` box declares (2), which is an independent cross-check between the
container's config and the elementary stream's own header.

CONTROLS, so a pass means something
------------------------------------
The same walk is re-run with the bit cursor deliberately shifted by 1, 2 and 3
bits.  Those MUST fail.  A gate whose wrong hypotheses also pass is not a gate
-- the same rule this project applied to the frequency de-interleaver sweep.

Usage:
    python m17_ac4_walk.py [m7_out/rf33_audio_pid13.mp4]
"""
from __future__ import annotations

import argparse
import collections
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m7_play as PL                                             # noqa: E402


def samples(path):
    """-> [bytes] one raw AC-4 frame per media sample, in order."""
    d = open(path, "rb").read()
    out = []
    o = 0
    while o + 8 <= len(d):
        sz = int.from_bytes(d[o:o + 4], "big")
        typ = d[o + 4:o + 8]
        if sz < 8:
            break
        if typ == b"moof":
            frag = d[o:]
            sc = PL.trun_scan(frag)
            if sc is not None:
                me = sc["media"]
                p = me["doff"]
                for s in me["sizes"]:
                    out.append(frag[p:p + s])
                    p += s
        o += sz
    return out


class BR:
    def __init__(self, d, skip=0):
        self.d, self.p = d, skip

    def u(self, n):
        v = 0
        for _ in range(n):
            if (self.p >> 3) >= len(self.d):
                return None
            v = (v << 1) | ((self.d[self.p >> 3] >> (7 - (self.p & 7))) & 1)
            self.p += 1
        return v


def variable_bits(r, n_bits):
    """TS 103 190-1 Table 3, verbatim -- and NOT what I had guessed.

    I had written the obvious thing, `value = (value << n) | chunk`.  The spec
    is not that: on each continuation it shifts AND adds an offset, so the
    ranges of successive lengths do not overlap:

        value = 0
        do {
            value += read(n_bits)
            if (b_read_more) { value <<= n_bits ; value += (1 << n_bits) }
        } while (b_read_more)

    It never fired on this stream -- `bitstream_version` is 2, so the escape at
    3 is untaken -- which is exactly why it is worth fixing now rather than
    when it silently mis-parses something later.
    """
    value = 0
    while True:
        chunk = r.u(n_bits)
        more = r.u(1)
        if chunk is None or more is None:
            return None
        value += chunk
        if not more:
            return value
        value <<= n_bits
        value += (1 << n_bits)


def toc_head(frame, skip=0):
    """-> (bitstream_version, sequence_counter) from ac4_toc, or None."""
    r = BR(frame, skip)
    bv = r.u(2)
    if bv is None:
        return None
    if bv == 3:
        v = variable_bits(r, 2)
        if v is None:
            return None
        bv += v
    sc = r.u(10)
    if sc is None:
        return None
    return bv, sc


def toc_deep(frame, skip=0):
    """Further into ac4_toc, as far as the fields that can SELF-CHECK.

    TS 103 190-1 4.2.3, continuing past sequence_counter:

        b_wait_frames        1
        if b_wait_frames:  wait_frames 3 ; if wait_frames > 0: br_code 2
        fs_index             1
        frame_rate_index     4
        b_iframe_global      1
        b_single_presentation 1

    The value of this stopping point is that `fs_index` and `frame_rate_index`
    are ALSO carried in the `dac4` box, derived from a completely different
    parse of a completely different structure.  If reading them here at these
    offsets reproduces the container's 48 kHz and 29.97 Hz -- on every frame --
    then the intervening field widths are right, because a wrong width would
    shift them into garbage.  Five bits of agreement, 3480 times over.
    """
    r = BR(frame, skip)
    bv = r.u(2)
    if bv is None:
        return None
    if bv == 3:
        v = variable_bits(r, 2)
        if v is None:
            return None
        bv += v
    sc = r.u(10)
    wait = r.u(1)
    if wait is None:
        return None
    if wait:
        wf = r.u(3)
        if wf is None:
            return None
        if wf > 0 and r.u(2) is None:   # 2 bits, "reserved" in V1.3.1
            return None
    fs = r.u(1)
    fri = r.u(4)
    ifr = r.u(1)
    single = r.u(1)
    if single is None:
        return None
    return dict(bv=bv, seq=sc, wait=wait, fs=fs, fri=fri,
                iframe=ifr, single_pres=single)


def walk_deep(frames, expect_fs, expect_fri, verbose=True):
    """-> did every frame reproduce the container's fs and frame rate?"""
    fs_c, fri_c, ifr_c, sp_c, wf_c = (collections.Counter() for _ in range(5))
    bad = 0
    for f in frames:
        t = toc_deep(f)
        if t is None:
            bad += 1
            continue
        fs_c[t["fs"]] += 1
        fri_c[t["fri"]] += 1
        ifr_c[t["iframe"]] += 1
        sp_c[t["single_pres"]] += 1
        wf_c[t["wait"]] += 1
    FS = {0: 44100, 1: 48000}
    FR = {0: 23.976, 1: 24, 2: 25, 3: 29.97, 4: 30, 5: 47.95, 6: 48, 7: 50,
          8: 59.94, 9: 60, 10: 100, 11: 119.88, 12: 120, 13: 23.44}
    ok = (len(fs_c) == 1 and next(iter(fs_c)) == expect_fs
          and len(fri_c) == 1 and next(iter(fri_c)) == expect_fri and not bad)
    if verbose:
        f1, r1 = list(fs_c)[0] if fs_c else None, list(fri_c)[0] if fri_c else None
        print(f"    b_wait_frames      {dict(wf_c)}")
        print(f"    fs_index           {dict(fs_c)}  -> {FS.get(f1)} Hz")
        print(f"    frame_rate_index   {dict(fri_c)}  -> {FR.get(r1)} Hz")
        print(f"    b_iframe_global    {dict(ifr_c)}")
        print(f"    b_single_presentation {dict(sp_c)}")
        if bad:
            print(f"    {bad} frames ran out of bits before the header ended")
    return ok, dict(fs=dict(fs_c), fri=dict(fri_c), iframe=dict(ifr_c),
                    single=dict(sp_c), wait=dict(wf_c), bad=bad)


def walk(frames, skip=0, expect_bv=None, verbose=True):
    bvs = collections.Counter()
    ramp = broken = 0
    prev = None
    first_break = None
    for i, f in enumerate(frames):
        t = toc_head(f, skip)
        if t is None:
            broken += 1
            continue
        bv, sc = t
        bvs[bv] += 1
        if prev is not None:
            if sc == (prev + 1) % 1024:
                ramp += 1
            else:
                broken += 1
                if first_break is None:
                    first_break = (i, prev, sc)
        prev = sc
    n = max(len(frames) - 1, 1)
    frac = ramp / n
    ok = frac > 0.999 and (expect_bv is None or
                           (len(bvs) == 1 and next(iter(bvs)) == expect_bv))
    if verbose:
        top = bvs.most_common(3)
        print(f"    skip {skip} bit(s): sequence_counter +1 on "
              f"{ramp}/{n} transitions ({frac*100:6.2f}%), "
              f"bitstream_version {dict(top)}")
        if first_break and frac < 0.999:
            i, a, b = first_break
            print(f"        first break at frame {i}: {a} -> {b}")
    return ok, frac, dict(bvs)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="m17 ac4 walk")
    ap.add_argument("path", nargs="?",
                    default="m7_out/rf33_audio_pid13.mp4")
    ap.add_argument("--bv", type=int, default=2,
                    help="bitstream_version the dac4 box declares")
    a = ap.parse_args(argv)
    p = a.path if os.path.isabs(a.path) else os.path.join(HERE, a.path)

    print("M17 -- can we read the AC-4 bitstream?")
    print("=" * 72)
    fr = samples(p)
    tot = sum(len(f) for f in fr)
    print(f"  {len(fr)} media samples, {tot:,} bytes "
          f"({tot*8/max(len(fr),1)/1000:.1f} kbit per frame average)")
    if not fr:
        print("  no samples -- nothing to walk")
        return 2

    print("\n  the walk, at the correct bit offset")
    ok, frac, bvs = walk(fr, 0, expect_bv=a.bv)
    print(f"    dac4 declared bitstream_version {a.bv}; "
          f"TOC says {list(bvs)} -> "
          f"{'AGREE' if list(bvs) == [a.bv] else '*** DISAGREE ***'}")

    print("\n  CONTROLS -- these must FAIL")
    ctl = []
    for s in (1, 2, 3):
        c_ok, c_frac, _ = walk(fr, s)
        ctl.append(c_ok)

    print("\n  DEEPER -- fields the container ALSO carries, so they self-check")
    deep_ok, deep = walk_deep(fr, expect_fs=1, expect_fri=3)
    print(f"    dac4 says 48000 Hz / 29.97 Hz; the TOC of all {len(fr)} frames "
          f"says {'THE SAME' if deep_ok else '*** SOMETHING ELSE ***'}")

    print("\n" + "=" * 72)
    passed = ok and not any(ctl)
    if passed:
        print(f"  PASS -- the sequence_counter ramps cleanly over "
              f"{len(fr)-1} transitions and every shifted control collapses.")
        print("  We can read the bitstream.  dolbyTuna is worth starting.")
    elif ok and any(ctl):
        print("  INCONCLUSIVE -- a shifted control also passed, so the ramp is "
              "not evidence.")
    else:
        print(f"  FAIL -- only {frac*100:.2f}% of transitions ramp.  We are "
              "not parsing the TOC correctly;")
        print("  no decoder should be started on this basis.")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
