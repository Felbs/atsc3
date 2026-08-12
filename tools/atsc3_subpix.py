#!/usr/bin/env python3
"""atsc3_subpix.py -- text cues -> PGS (.sup) bitmap subtitles, for dvbsub.

WHY THIS EXISTS (E35). The v2 viewer carries captions as a real, toggleable
subtitle TRACK in the MPEG-TS instead of burning them into the picture. The
TS-native subtitle codec is DVB subtitles -- a BITMAP format -- and ffmpeg
refuses text->bitmap: "Subtitle encoding currently only possible from text
to text or bitmap to bitmap" (measured, 8/08). So we render the text
ourselves (PIL) into the simplest bitmap subtitle container ffmpeg can read
back, PGS (.sup), and let ffmpeg do the bitmap->bitmap transcode to dvbsub.

Two timing quirks, both measured, both load-bearing:

  * The pgssub decoder leaves each cue's end OPEN (PGS ends a cue with a
    separate clear-set), so the dvbsub encoder saw end_display_time =
    UINT32_MAX ms and stamped packets ~4295 s late (the literal 2^32).
    `-fix_sub_duration` on the .sup INPUT clamps every event at the next
    one -- with the clear-sets we emit, that reproduces the srt exactly.
  * ffmpeg rebases each input by its own start time, so a .sup whose first
    event is a cue at 7.3 s would slide every cue 7.3 s early. Every chunk
    therefore begins with a fully TRANSPARENT dummy cue at t=0 -- which
    also guarantees the dvbsub stream EXISTS in cue-less chunks, keeping
    the PID layout of every chunk identical (the concat contract).

Palette: 0 transparent, 1 white text, 2 dark box. Canvas is 720x576 (the
dvbsub encoder's SD world); VLC/ffmpeg scale it onto the 1280x720 video.
"""
from __future__ import annotations

import struct

CANVAS_W, CANVAS_H = 720, 576
FONT_CANDIDATES = (
    # Windows
    r"C:\Windows\Fonts\arialbd.ttf",
    r"C:\Windows\Fonts\arial.ttf",
    r"C:\Windows\Fonts\segoeui.ttf",
    # Linux (E56 port): DejaVu ships with Ubuntu, Liberation as backup
    "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
    "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",
    "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
)
FONT_SIZE = 26
WRAP_COLS = 46


def _seg(pts, typ, payload):
    t = int(round(pts * 90000))
    return b"PG" + struct.pack(">IIBH", t, t, typ, len(payload)) + payload


def _rle_line(row):
    """One bitmap row of palette indices -> PGS run-length code."""
    out = bytearray()
    i, n = 0, len(row)
    while i < n:
        c = row[i]
        j = i
        while j < n and row[j] == c:
            j += 1
        run = j - i
        if c == 0:
            while run > 0:
                r = min(run, 16383)
                if r < 64:
                    out += bytes((0, r))
                else:
                    out += bytes((0, 0x40 | (r >> 8), r & 0xFF))
                run -= r
        else:
            while run > 0:
                r = min(run, 16383)
                if r <= 2:
                    out += bytes((c,)) * r
                    run = 0
                elif r < 64:
                    out += bytes((0, 0x80 | r, c))
                    run -= r
                else:
                    out += bytes((0, 0xC0 | (r >> 8), r & 0xFF, c))
                    run -= r
        i = j
    out += bytes((0, 0))
    return bytes(out)


class SupRenderer:
    """Renders cue text to PGS display sets, with a small bitmap cache
    (roll-up captions repeat the same line across many chunks)."""

    def __init__(self):
        from PIL import ImageFont
        self.font = None
        for cand in FONT_CANDIDATES:
            try:
                self.font = ImageFont.truetype(cand, FONT_SIZE)
                break
            except OSError:
                continue
        if self.font is None:
            raise OSError("no usable font for subtitle rendering")
        self._cache = {}

    def _wrap(self, text):
        lines = []
        for raw in text.splitlines():
            raw = raw.strip()
            while len(raw) > WRAP_COLS:
                k = raw.rfind(" ", 0, WRAP_COLS)
                k = k if k > 0 else WRAP_COLS
                lines.append(raw[:k])
                raw = raw[k:].strip()
            if raw:
                lines.append(raw)
        return lines or [" "]

    def _bitmap(self, text):
        """text -> (w, h, palette-index bytes). Cached."""
        hit = self._cache.get(text)
        if hit is not None:
            return hit
        from PIL import Image, ImageDraw
        lines = self._wrap(text)
        pad = 8
        img = Image.new("L", (4, 4))
        dr = ImageDraw.Draw(img)
        boxes = [dr.textbbox((0, 0), l, font=self.font) for l in lines]
        lh = max(b[3] - b[1] for b in boxes) + 6
        w = min(CANVAS_W - 8, max(b[2] - b[0] for b in boxes) + 2 * pad + 4)
        h = lh * len(lines) + 2 * pad
        img = Image.new("P", (w, h), 2)
        dr = ImageDraw.Draw(img)
        for k, l in enumerate(lines):
            bb = dr.textbbox((0, 0), l, font=self.font)
            x = max(2, (w - (bb[2] - bb[0])) // 2 - bb[0])
            y = pad + k * lh - bb[1] + 3
            dr.text((x, y), l, font=self.font, fill=1)
        out = (w, h, img.tobytes())
        if len(self._cache) > 128:
            self._cache.clear()
        self._cache[text] = out
        return out

    def _display_set(self, pts, comp, w, h, bitmap, x, y):
        pcs = (struct.pack(">HHBHBBBB", CANVAS_W, CANVAS_H, 0x10,
                           comp & 0xFFFF, 0x80, 0, 0, 1)
               + struct.pack(">HBBHH", 0, 0, 0, x, y))
        wds = struct.pack(">BBHHHH", 1, 0, x, y, w, h)
        pal = (bytes((0, 0))
               + bytes((0, 16, 128, 128, 0))          # 0: transparent
               + bytes((1, 235, 128, 128, 255))       # 1: white text
               + bytes((2, 16, 128, 128, 200)))       # 2: dark box
        rle = b"".join(_rle_line(bitmap[r * w:(r + 1) * w]) for r in range(h))
        ods = (struct.pack(">HBB", 0, 0, 0xC0)
               + struct.pack(">I", len(rle) + 4)[1:]
               + struct.pack(">HH", w, h) + rle)
        return (_seg(pts, 0x16, pcs) + _seg(pts, 0x17, wds)
                + _seg(pts, 0x14, pal) + _seg(pts, 0x15, ods)
                + _seg(pts, 0x80, b""))

    def _clear_set(self, pts, comp, wds_echo):
        pcs = struct.pack(">HHBHBBBB", CANVAS_W, CANVAS_H, 0x10,
                          comp & 0xFFFF, 0x00, 0, 0, 0)
        return _seg(pts, 0x16, pcs) + _seg(pts, 0x17, wds_echo) \
            + _seg(pts, 0x80, b"")

    def sup_bytes(self, cues, base=0.0):
        """cues: [(begin_s, end_s, text)], absolute, >= base, sorted.

        Always begins with the transparent dummy at t=base (see module
        doc); the caller muxes with -copyts so `base` is the chunk's
        session offset and every timestamp is carried literally.
        Overlapping cues are serialised -- each show clamps the previous
        one -- which matches how the collapsed roll-up srt behaves.
        """
        out = bytearray()
        comp = 0
        # dummy: 8x8 transparent object at origin, shown base .. base+0.2
        w = h = 8
        bm = bytes(w * h)
        out += self._display_set(base, comp, w, h, bm, 0, 0)
        wds = struct.pack(">BBHHHH", 1, 0, 0, 0, w, h)
        out += self._clear_set(base + 0.2, comp + 1, wds)
        comp += 2
        for b, e, text in cues:
            if e <= b:
                continue
            b = max(b, base + 0.25)  # never collide with the dummy's window
            w, h, bm = self._bitmap(text)
            x = (CANVAS_W - w) // 2
            y = CANVAS_H - h - 40
            out += self._display_set(b, comp, w, h, bm, x, y)
            wds = struct.pack(">BBHHHH", 1, 0, x, y, w, h)
            out += self._clear_set(e, comp + 1, wds)
            comp += 2
        return bytes(out)
