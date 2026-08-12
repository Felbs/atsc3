"""M40 -- closed captions: IMSC1/TTML out of the stpp asset, onto our timeline.

The pid15 asset carries one TTML document per MPU, one sample each, 59 of
them -- the only asset with no losses at all.  Two things make it more than
an XML dump:

1. THE TIMES ARE ABSOLUTE.  `begin="1786023159.304s"` is wall clock, not
   media time, even though ttp:timeBase says "media".  They have to be
   anchored to the MPU grid, and the anchor is recovered (not assumed) by
   fitting each document's earliest cue against its own slot's 2.002 s
   window.  The residual of that fit is the gate: if the anchor were wrong
   the cues would not land inside their own MPUs.

2. THEY ARE ROLL-UP CAPTIONS.  Live captioning re-emits the entire
   accumulated line every time a character or two arrives:

       NO.\\nSHE'S  ->  NO.\\nSHE'S W  ->  NO.\\nSHE'S WON  -> ...

   so the 59 documents hold thousands of cues that are almost all prefixes
   of their successor.  Converted naively you get unreadable strobing.  We
   collapse each growth run to the text it settled on, and start a new cue
   when the text stops being an extension of what came before (a scroll).

The dropped video MPU is applied here too, so the captions stay aligned
with the muxed file rather than with the original transmission.
"""
import os, re, sys, struct, subprocess
import xml.etree.ElementTree as ET

LAB = os.path.dirname(os.path.abspath(__file__))
# Directory of MPU segments.  Defaults to the 8/06 bank; pass another to run
# against a live capture rebuilt by m7_objects.
OUT = os.path.join(LAB, sys.argv[1] if len(sys.argv) > 1 else "m7_out")

TTML_NS = "{http://www.w3.org/ns/ttml}"
MPU_SECONDS = 120 * 1001 / 60000.0        # 2.002
# BASE_SLOT and the video hole were hardcoded for the 8/06 bank.  Both are
# properties of whichever capture we are handed, so derive them.
BASE_SLOT = None
DROPPED_VIDEO_SLOT = None


def survey():
    """Base slot, and which slots the VIDEO asset is missing."""
    got = {}
    for f in os.listdir(OUT):
        m = re.match(r"mpu_.*_(pid\d+)_(\d+)\.seg$", f)
        if m:
            got.setdefault(m.group(1), set()).add(int(m.group(2)))
    if not got:
        raise SystemExit(f"no MPU segments in {OUT}")
    # Anchor to the VIDEO, not to the earliest slot of any asset. Captions
    # have to line up with the file we actually muxed, and the muxer starts
    # at the first video MPU -- the subtitle asset routinely starts a slot or
    # two earlier (it loses nothing, video does), and anchoring to that put
    # every cue 2 slots = 4.004 s early.
    vid = got.get("pid12", set())
    base = min(vid) if vid else min(min(v) for v in got.values())
    hi = max(vid) if vid else max(max(v) for v in got.values())
    holes = [s for s in range(base, hi + 1) if s not in vid]
    return base, holes


def boxes(b, off=0, end=None):
    if end is None:
        end = len(b)
    p = off
    while p + 8 <= end:
        sz = struct.unpack(">I", b[p:p + 4])[0]
        t = b[p + 4:p + 8]
        if sz == 0:
            sz = end - p
        elif sz == 1:
            sz = struct.unpack(">Q", b[p + 8:p + 16])[0]
        if sz < 8 or p + sz > end:
            break
        yield t, p + 8, p + sz
        p += sz


def mdat(seg):
    for t, s, e in boxes(seg):
        if t == b"mdat":
            return seg[s:e]
    return b""


def text_of(p):
    """Flatten a <p>, turning <br/> into newline."""
    parts = []

    def walk(el):
        if el.tag == TTML_NS + "br":
            parts.append("\n")
        if el.text:
            parts.append(el.text)
        for c in el:
            walk(c)
            if c.tail:
                parts.append(c.tail)

    if p.text:
        parts.append(p.text)
    for c in p:
        walk(c)
        if c.tail:
            parts.append(c.tail)
    return "".join(parts).strip()


def secs(v):
    v = v.strip()
    if v.endswith("s"):
        return float(v[:-1])
    if ":" in v:                       # hh:mm:ss.ms
        bits = [float(x) for x in v.split(":")]
        t = 0.0
        for b in bits:
            t = t * 60 + b
        return t
    return float(v)


def load():
    """Every cue from every subtitle MPU, with the slot it came from."""
    files = sorted(
        (int(re.search(r"_(\d+)\.seg$", f).group(1)), f)
        for f in os.listdir(OUT) if re.match(r"mpu_.*_pid15_\d+\.seg$", f))
    cues, per_mpu = [], {}
    for seq, f in files:
        raw = mdat(open(os.path.join(OUT, f), "rb").read())
        end = raw.rfind(b"</tt>")
        if end < 0:
            continue
        try:
            root = ET.fromstring(raw[:end + 5].decode("utf-8", "replace"))
        except ET.ParseError as ex:
            print(f"  seq {seq}: parse error {ex}")
            continue
        got = []
        for p in root.iter(TTML_NS + "p"):
            b, e = p.get("begin"), p.get("end")
            if not b or not e:
                continue
            txt = text_of(p)
            if txt:
                got.append((secs(b), secs(e), txt))
        per_mpu[seq] = got
        cues.extend((seq, b, e, t) for b, e, t in got)
    return cues, per_mpu


def anchor_of(per_mpu):
    """Recover absolute->media offset by fitting each MPU's own window."""
    est = []
    for seq, got in per_mpu.items():
        if not got:
            continue
        i = seq - BASE_SLOT
        est.append(min(b for b, _, _ in got) - i * MPU_SECONDS)
    est.sort()
    return est[len(est) // 2], est


def collapse(cues):
    """Roll-up runs -> one cue per completed LINE.

    Tracking the whole caption block is wrong: a roll-up scrolls, so the
    block gets SHORTER when the top line falls off, and a run keyed on
    "longest text wins" swallows the new line and throws it away (this is
    how AGO and ANGRY went missing).  The unit that actually grows
    monotonically is the bottom line -- the one being typed.  A run ends
    when that line stops being an extension of itself, which is exactly
    when the caption operator finished it.
    """
    cues = sorted(cues, key=lambda c: (c[1], c[2]))
    out, cur = [], None
    for seq, b, e, t in cues:
        line = t.splitlines()[-1].strip() if t.strip() else ""
        if not line:
            continue
        if cur is not None and (line.startswith(cur[2]) or line == cur[2]):
            cur[1] = max(cur[1], e)
            cur[2] = line if len(line) > len(cur[2]) else cur[2]
        else:
            if cur:
                out.append(cur)
            cur = [b, e, line]
    if cur:
        out.append(cur)

    merged = []
    for g in out:
        if merged and merged[-1][2] == g[2]:
            merged[-1][1] = max(merged[-1][1], g[1])
        else:
            merged.append(g)
    return merged


def srt_time(t):
    if t < 0:
        t = 0.0
    h = int(t // 3600); t -= h * 3600
    m = int(t // 60); t -= m * 60
    s = int(t)
    ms = int(round((t - s) * 1000))
    if ms == 1000:
        s, ms = s + 1, 0
    return f"{h:02d}:{m:02d}:{s:02d},{ms:03d}"


def main():
    global BASE_SLOT, DROPPED_VIDEO_SLOT
    BASE_SLOT, holes = survey()
    DROPPED_VIDEO_SLOT = holes[0] if holes else None
    print(f"capture: base slot {BASE_SLOT}, video holes {holes}")
    cues, per_mpu = load()
    print(f"{len(per_mpu)} subtitle MPUs parsed, {len(cues)} raw cues")

    anchor, est = anchor_of(per_mpu)
    spread = max(est) - min(est)
    print(f"\nanchor (absolute -> media) {anchor:.3f} s")
    print(f"  per-MPU estimates span {spread:.3f} s over {len(est)} MPUs")

    # GATE 1: with that anchor, does every cue land inside its own MPU?
    bad = 0
    for seq, b, e, t in cues:
        i = seq - BASE_SLOT
        m = b - anchor
        if not (i * MPU_SECONDS - 0.1 <= m < (i + 1) * MPU_SECONDS + 0.1):
            bad += 1
    print(f"\nGATE 1  cues landing inside their own MPU window: "
          f"{len(cues) - bad}/{len(cues)}"
          f"   {'PASS' if bad == 0 else 'FAIL (' + str(bad) + ' outside)'}")

    med = collapse([(s, b - anchor, e - anchor, t) for s, b, e, t in cues])
    print(f"\ncollapsed {len(cues)} raw cues -> {len(med)} readable cues"
          f"   ({100 * (1 - len(med) / max(len(cues), 1)):.1f} % were roll-up)")

    # GATE 2: collapsing must not lose any COMPLETED word.
    #
    # Comparing against every token in the raw cues is the wrong test -- most
    # of them are half-typed fragments (ABOU, ALWAY, ANGR) whose whole
    # purpose is to disappear.  A token is a real word only once typing has
    # moved past it, i.e. it appears somewhere with something after it on the
    # same line.  A fragment is always the last token of its cue and never
    # anything else, so this distinguishes them without reference to how the
    # collapse works -- the gate stays independent of the thing it checks.
    # Only the BOTTOM line is being typed; the lines above it are history
    # that scrolled up. On a mid-stream capture (i.e. live) the top lines
    # contain words typed BEFORE we tuned in, which collapse never had the
    # chance to keep -- demanding them fails the gate for something that is
    # unrecoverable by construction rather than for a defect. Restrict to
    # words that were observed being typed inside this capture.
    complete = set()
    for _, _, _, t in cues:
        lines = t.splitlines()
        if lines:
            complete |= set(lines[-1].split()[:-1])
    kept = set()
    for _, _, t in med:
        kept |= set(t.replace("\n", " ").split())
    lost = complete - kept
    print(f"GATE 2  completed words preserved: "
          f"{len(complete - lost)}/{len(complete)}"
          f"   {'PASS' if not lost else 'FAIL'}")
    if lost:
        print(f"        lost: {sorted(lost)[:20]}")

    # apply the dropped video MPU so captions match the muxed file
    cut_at = ((DROPPED_VIDEO_SLOT - BASE_SLOT) * MPU_SECONDS
              if DROPPED_VIDEO_SLOT is not None else float("inf"))
    n_slots = max(s for s in per_mpu) - BASE_SLOT + 1
    dur = (n_slots - (1 if DROPPED_VIDEO_SLOT is not None else 0)) * MPU_SECONDS
    out, past_end, in_hole, before_start = [], 0, 0, 0
    for b, e, t in med:
        if b >= cut_at + MPU_SECONDS:
            b, e = b - MPU_SECONDS, e - MPU_SECONDS
        elif e > cut_at:
            in_hole += 1                # cue sat in the slot we dropped
            continue
        if b >= dur:
            past_end += 1               # subtitle outruns the media
            continue
        if b < 0:
            # The subtitle asset starts before the first video MPU, so these
            # cues belong to air we did not keep. Clamping them to 0 stacked
            # four cues on the first frame; drop them instead.
            before_start += 1
            continue
        out.append((b, min(e, dur), t))
    print(f"\ndropped-MPU adjustment at t={cut_at:.3f} s: {len(med)} -> "
          f"{len(out)} cues  ({in_hole} inside the hole, "
          f"{past_end} past the end, {before_start} before it starts)")

    # GATE 3: monotone and inside the programme
    mono = all(out[i][0] <= out[i + 1][0] for i in range(len(out) - 1))
    inside = all(0 <= b < dur and b < e <= dur for b, e, _ in out)
    print(f"GATE 3  monotone {mono}, all within 0..{dur:.1f} s {inside}"
          f"   {'PASS' if mono and inside else 'FAIL'}")

    dst = os.path.join(LAB, "rf33_tv.srt")
    with open(dst, "w", encoding="utf-8") as fh:
        for i, (b, e, t) in enumerate(out, 1):
            if e <= b:
                e = b + 0.5
            fh.write(f"{i}\n{srt_time(b)} --> {srt_time(e)}\n{t}\n\n")
    print(f"\nwrote {os.path.basename(dst)}   {len(out)} cues")

    print("\nfirst lines of dialogue:")
    for b, e, t in out[:6]:
        print(f"  {srt_time(b)}  {t.splitlines()[-1][:60]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
