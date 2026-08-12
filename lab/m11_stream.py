#!/usr/bin/env python3
"""m11_stream.py -- M11: the four pieces the M9 design note named, built.

M9 stopped at a design plus a measured budget rather than half-built plumbing:

    STAGE 1  continuous front end   135.7 ms/Frame vs a 247.1 ms budget
    STAGE 2  frame decoder          107.7 ms/Frame per worker
    STAGE 3  transport                1.06 ms/Frame
    STAGE 4  MPU -> player          2.002 s of INHERENT MPU latency

and it named the four genuinely new pieces.  All four are here:

  1. `PolyResampler` -- `resample_poly` with phase state carried across block
     boundaries, so a block seam is not a discontinuity.  Gated bit-identical
     against `scipy.signal.resample_poly` on the whole array.
  2. `FrontEnd` -- ONE bootstrap search, then t0 tracks by +FRAME_SAMPLES with
     the existing +-20 fine-timing search.  De-rotation carries a continuous
     sample index; the notch FIR carries `len(b)-1` samples of real history,
     which is the same overlap-save identity M9 proved for `par_lfilter`.
  3. `StreamAlpWalker` -- a resumable A/330 walk.  M7 already paid for the
     lesson here: walking chunks separately lost one ALP packet per boundary
     and punched a 1400-byte hole in essentially every MPU, and before that a
     mis-read length swallowed 62% of the payload because no anchor could
     contradict it.  So this walker does NOT guess: it refuses to parse at all
     until it holds `MAX_ALP + 8` bytes of lookahead **and** the boundary list
     that covers them, which makes every reading it takes provably the same
     reading the whole-stream walker takes.  Nothing is parsed speculatively
     and then un-parsed -- `_segment` mutates state, so a speculative parse
     would be a bug, not an inefficiency.
  4. `MpuStreamer` -- emit an MPU the instant it is complete.  The completeness
     test and the bytes emitted are `m7_objects.MmtpFlow.build()` itself,
     called on one MPU, so "the streaming reassembly agrees with the batch
     reassembly" is true by construction rather than by comparison.

THE RESAMPLER IS NORMALLY NOT USED, AND THAT IS THE POINT
---------------------------------------------------------
M9's profile said `resample_poly` was 70.2 ms of a 243.4 ms Frame -- 29% of
the pipeline -- and that it existed only because the capture was at 8 Msps.
The radio was asked for 6.912 Msps and the readback says 6912000.0 Hz exactly,
so on live air the stage is not optimised, it is absent.  `PolyResampler` is
kept and gated anyway because banked 8 Msps captures still have to replay
through the same front end, and because "we removed a stage" deserves a
fallback that is known to work rather than a deleted code path.
"""
from __future__ import annotations

import bisect
import collections
import statistics
import math
import os
import queue
import sys
import threading
import time
from fractions import Fraction

import numpy as np
from scipy.signal import firwin, lfilter, upfirdn
from scipy import fft as SFFT

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m2_pilots as MP                                            # noqa: E402
import m6_cells as C                                              # noqa: E402
import m9_fast as F9                                              # noqa: E402
import m7_objects as O7                                           # noqa: E402
import m7_play as PL                                              # noqa: E402
import m7_route as R7                                             # noqa: E402
from atsc3 import bootstrap as bs                                 # noqa: E402
from atsc3.capture import DTYPES                                  # noqa: E402

FS_POST = 6.912e6
FRAME_SAMPLES = C.FRAME_SAMPLES
FRAME_SEC = FRAME_SAMPLES / FS_POST
NOTCH_TAPS = 257
NOTCH_CUT = 2.95e6

# demod_data of the last data symbol (l = NSYM-1 = 34) reads
#   y[t0 + (NFFT+GI)*(1+l) + GI : ... + NFFT] = y[t0+342016 : t0+350208]
# so a Frame decode needs [t0, t0+350208).  fine_timing scans t0 over
# +-FT_SPAN around the tracked centre, and everything it touches is inside
# that.  A little slack on top, and this is the whole per-Frame window.
FT_SPAN = 20
FRAME_WINDOW = 350208 + 2 * FT_SPAN + 4096


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


# ---------------------------------------------------------------------------
# 1 -- the resampler, with phase state
# ---------------------------------------------------------------------------

class PolyResampler:
    """`scipy.signal.resample_poly`, fed a block at a time.

    scipy computes `y_full[m] = sum_j h[m*down - j*up] * x[j]` on a filter it
    zero-pads so the output samples land in the centre, then keeps
    `y_full[n_pre_remove:]`.  Blocking that exactly needs two things:

      * `len(h)` samples of real input history, so no output sample is ever
        missing a tap -- the same overlap-save argument as the notch FIR, and
        exact for the same reason (each output is one dot product and it sees
        the same summands in the same order);
      * a block boundary at which the polyphase PHASE is zero.  Output m sits
        in branch `(m*down) % up`, so a local call starting at input index
        `n0 - H` reproduces branch and order iff `(n0 - H) * up` is divisible
        by `down`.  With `gcd(up, down) == 1` that is `n0 % down == 0` and
        `H % down == 0`, which is why blocks are constrained to multiples of
        `down` rather than to a convenient power of two.

    A resampler that just re-called `resample_poly` per block would restart
    the filter from zeros every time: every seam becomes a step, 5000 times a
    second.  This is that bug, not fixed but made unable to occur.
    """

    def __init__(self, fs_in, fs_out, window=("kaiser", 5.0)):
        fr = Fraction(int(round(fs_out)),
                      int(round(fs_in))).limit_denominator(20000)
        up, down = fr.numerator, fr.denominator
        g = math.gcd(up, down)
        self.up, self.down = up // g, down // g
        self.bypass = (self.up == 1 and self.down == 1)
        if self.bypass:
            self.hist = None
            return
        max_rate = max(self.up, self.down)
        half_len = 10 * max_rate
        # THE TAP DTYPE IS PART OF THE ANSWER, not an implementation detail.
        # `resample_poly` does `h = asarray(h, dtype=x.dtype)` and only THEN
        # `h *= up`, so on a complex64 stream it convolves with complex64 taps
        # and accumulates in single precision.  Designing in float64 and
        # scaling before the cast -- the obvious way to write this -- is more
        # accurate and therefore WRONG here: it left 1.58 M of 1.73 M output
        # samples differing in the last float32 digits, which `array_equal`
        # correctly refused to call bit-identical.  Match scipy's order.
        h = np.asarray(firwin(2 * half_len + 1, 1.0 / max_rate, window=window),
                       dtype=np.complex64) * np.complex64(self.up)
        n_pre_pad = self.down - half_len % self.down
        self.n_pre_remove = (half_len + n_pre_pad) // self.down
        self.h = np.concatenate((np.zeros(n_pre_pad, h.dtype), h))
        # history: a multiple of `down`, long enough that every tap is real
        self.H = self.down * max(1, -(-len(self.h) // (self.up * self.down)))
        self.hist = np.zeros(self.H, np.complex64)
        self.off = self.H * self.up // self.down
        self.dropped = 0

    def block_multiple(self):
        return 1 if self.bypass else self.down

    def feed(self, x):
        """x: complex64, len(x) % down == 0.  -> complex64, same dtype path as
        `m2_stage_b.resample_to` (which casts the complex128 upfirdn result
        back to complex64)."""
        if self.bypass:
            return np.asarray(x, np.complex64)
        if len(x) % self.down:
            raise ValueError("block length must be a multiple of down")
        buf = np.concatenate((self.hist, x))
        y = upfirdn(self.h, buf, self.up, self.down)
        n_out = len(x) * self.up // self.down
        out = y[self.off:self.off + n_out].astype(np.complex64)
        self.hist = buf[-self.H:]
        if self.dropped < self.n_pre_remove:
            k = min(self.n_pre_remove - self.dropped, len(out))
            self.dropped += k
            out = out[k:]
        return out


def gate_resampler(fs_in=8e6, fs_out=FS_POST, n=2_000_000, block=250_000,
                   seed=7, verbose=True):
    """Blocked == whole-array, bit for bit.  `array_equal`, not `allclose`."""
    from m2_stage_b import resample_to
    rng = np.random.default_rng(seed)
    x = (rng.standard_normal(n) + 1j * rng.standard_normal(n)).astype(
        np.complex64)
    ref = resample_to(x, fs_in, fs_out)
    r = PolyResampler(fs_in, fs_out)
    outs = [r.feed(x[i:i + block]) for i in range(0, n, block)]
    got = np.concatenate(outs) if outs else np.zeros(0, np.complex64)
    k = len(got)
    ok = k > 0 and np.array_equal(ref[:k], got)
    if verbose:
        print(f"    {'PASS' if ok else 'FAIL'}  PolyResampler {fs_in/1e6:g} ->"
              f" {fs_out/1e6:g} Msps: {k} of {len(ref)} output samples "
              f"bit-identical to resample_poly (up={r.up} down={r.down}, "
              f"{block}-sample blocks)")
        if not ok and k:
            d = np.abs(ref[:k] - got)
            print(f"          first mismatch at {int(np.argmax(d>0))}, "
                  f"max |delta| {d.max():.3e}")
    return ok


# ---------------------------------------------------------------------------
# a tape: append at the end, drop from the front, index ABSOLUTELY
# ---------------------------------------------------------------------------

class SampleTape:
    """Complex sample stream with absolute indices and a trimmable head."""

    def __init__(self):
        self.blocks = []
        self.org = 0
        self.end = 0

    def append(self, a):
        if len(a):
            self.blocks.append(a)
            self.end += len(a)

    def take(self, lo, n):
        if lo < self.org or lo + n > self.end:
            raise IndexError(f"take({lo},{n}) outside [{self.org},{self.end})")
        out = np.empty(n, self.blocks[0].dtype)
        p = self.org
        w = 0
        for b in self.blocks:
            s, e = p, p + len(b)
            p = e
            if e <= lo:
                continue
            if s >= lo + n:
                break
            a = max(lo, s) - s
            z = min(lo + n, e) - s
            out[w:w + z - a] = b[a:z]
            w += z - a
        return out

    def trim(self, lo):
        while self.blocks and self.org + len(self.blocks[0]) <= lo:
            self.org += len(self.blocks.pop(0))


class ByteTape:
    """bytearray with absolute indices, so `AlpWalker`'s readings -- which are
    written against one whole-stream bytearray -- run unmodified over a window
    that is being trimmed from the front."""

    def __init__(self):
        self.b = bytearray()
        self.org = 0

    def __len__(self):
        return self.org + len(self.b)

    def __getitem__(self, i):
        if isinstance(i, slice):
            return self.b[(i.start or 0) - self.org:
                          (len(self) if i.stop is None else i.stop) - self.org]
        return self.b[i - self.org]

    def extend(self, d):
        self.b += d

    def trim(self, upto):
        k = upto - self.org
        if k > 0:
            del self.b[:k]
            self.org += k


# ---------------------------------------------------------------------------
# 3 -- the resumable ALP walk
# ---------------------------------------------------------------------------

class StreamAlpWalker(R7.AlpWalker):
    """`m7_route.AlpWalker`, resumable, and provably taking the same readings.

    The header readings (`_single`, `_concat`, `_segment`) are the parent's,
    unmodified.  What is new is WHEN they are allowed to run.

    A whole-stream walker decides at `pos` with the whole stream and the whole
    anchor list in hand.  A streaming walker that decided with less would
    sometimes take a different reading -- and worse, `_segment` writes into
    `self.segbuf`, so a speculative parse that is later abandoned corrupts
    state rather than merely wasting time.

    So: no speculation.  A reading is attempted only once the tape holds
    `MAX_ALP + 8` bytes past `pos` AND the last known Baseband anchor is at
    least that far past `pos` too.  Then every length check has the same
    answer it has on the whole stream (nothing legal can fail to fit), and
    every `_anchor_clean` sees every anchor that could contradict the
    candidate -- because anchors arrive in increasing order, so "all anchors
    below the last one" is "all anchors".  Cost: about 65 kB of lookahead,
    which is two thirds of a Frame, and it disappears inside the 2.002 s an
    MPU takes to exist at all.
    """

    LOOKAHEAD = R7.AlpWalker.MAX_ALP + 8

    def __init__(self):
        super().__init__(ByteTape(), [])
        self.bounds = []
        self.pos = None
        self.emitted = 0

    def feed(self, chunk, bounds_abs):
        """chunk: bytes appended at the end.  bounds_abs: absolute anchor
        positions inside it, increasing.  -> list of (packet_type, payload)."""
        self.s.extend(chunk)
        for b in bounds_abs:
            if not self.bounds or b > self.bounds[-1]:
                self.bounds.append(b)
        return self.pump()

    def pump(self, final=False):
        out = []
        if not self.bounds:
            return out
        if self.pos is None:
            self.pos = self.bounds[0]
        b = self.s
        limit = len(b) if final else self.bounds[-1]
        while self.pos + 2 <= len(b):
            if not final and limit - self.pos < self.LOOKAHEAD:
                break
            pos = self.pos
            pc = (b[pos] >> 4) & 1
            if pc == 0:
                r = self._single(pos)
                self.stats["pc0"] += 1
            else:
                sc = (b[pos] >> 3) & 1
                r = self._concat(pos) if sc else self._segment(pos)
                self.stats["pc1"] += 1
            if r is None:
                i = bisect.bisect_right(self.bounds, pos)
                if i >= len(self.bounds):
                    break
                self.pos = self.bounds[i]
                self.stats["resync"] += 1
                continue
            pkts, consumed = r
            out += pkts
            self.pos += consumed
        self.emitted += len(out)
        # release everything behind the cursor
        self.s.trim(self.pos)
        i = bisect.bisect_left(self.bounds, self.pos)
        if i > 64:
            del self.bounds[:i]
        return out


def gate_alp_walker(stream, bounds, chunks=37, verbose=True):
    """The resumable walk must equal the whole-stream walk, packet for packet.

    The reference is re-derived here from the same bytes, not read from a
    stored expectation -- M9's own lesson about gating against a memory.
    """
    ref, rst = R7.AlpWalker(bytearray(stream), list(bounds)).run()
    w = StreamAlpWalker()
    got = []
    n = len(stream)
    step = max(1, n // chunks)
    bl = sorted(set(bounds))
    j = 0
    for s in range(0, n, step):
        e = min(n, s + step)
        k = bisect.bisect_left(bl, e)
        got += w.feed(bytes(stream[s:e]), bl[j:k])
        j = k
    got += w.pump(final=True)
    ok = len(ref) == len(got) and all(a == b for a, b in zip(ref, got))
    if verbose:
        print(f"    {'PASS' if ok else 'FAIL'}  StreamAlpWalker over "
              f"{chunks} chunks: {len(got)} ALP packets vs {len(ref)} from the "
              f"whole-stream walk, {w.stats['resync']} vs {rst['resync']} "
              f"resyncs")
        if not ok:
            for i, (a, b) in enumerate(zip(ref, got)):
                if a != b:
                    print(f"          first difference at packet {i}")
                    break
    return ok


# ---------------------------------------------------------------------------
# 4 -- MPUs, emitted the moment they complete
# ---------------------------------------------------------------------------

def moov_handler(meta):
    """The `hdlr` handler_type of an MPU metadata fragment: b'vide', b'soun'.

    This is how the video asset is identified.  Reading it off the moov the
    transmitter sent is the only way that does not hard-code a packet_id --
    and a hard-coded packet_id is exactly the sort of thing that works on one
    multiplex and silently picks the audio on the next one.
    """
    h = PL.find_box(meta, [b"moov", b"trak", b"mdia", b"hdlr"])
    if not h:
        return None
    o = h[0]
    return bytes(meta[o + 16:o + 20])


class MpuStreamer:
    """MMTP datagrams in, complete MPUs out, at the instant they complete.

    Reassembly is not reimplemented: the bytes come from
    `m7_objects.MmtpFlow.build()` run over exactly one MPU.  That keeps this
    honest by construction -- there is no second reassembler to drift.
    """

    # How far ahead of what we have seen an mpu_sequence_number is allowed to
    # jump before we treat it as a corrupted field rather than a real one.
    # 1024 MPUs is 34 minutes of air -- far longer than any outage we could
    # sit through and still be receiving, and far short of the 2**31 that one
    # flipped bit produces.
    SEQ_GUARD = 1024
    SEQ_RESYNC = 3           # packets that must AGREE before we believe a jump

    def __init__(self, ver_ext, keep=3, repair=True):
        self.ver_ext = ver_ext
        self.keep = keep
        self.repair = repair
        self.open = {}
        self.done = set()
        self.hi = {}
        self.stats = collections.Counter()
        self.psn_last = {}
        self.resync = {}

    def _seq_ok(self, pid, seq):
        """Is this mpu_sequence_number believable?

        THE MMTP HEADER HAS NO CRC.  `mpu_payload` reads the field straight out
        of the datagram, so a single bit error inside a fade yields a sequence
        number up to 2**31 away from the truth -- and both consumers of it are
        MONOTONE: `hi[pid]` is a max(), and `Transport.next_seq` only moves
        forward.  A latched wild value makes `_reap` cut above every real MPU,
        so each one is reaped the instant it opens and every later fragment of
        it is dropped as `late_fragment`.  The asset stops emitting FOREVER
        while the chain reports perfect FEC.

        E20 lost video and both audio lanes this way in one fade on 8/07; the
        subtitle lane simply never drew a bad packet.  `gate_seq_poison.py`
        reproduces it: 10/10 MPUs before the flip, 0/30 after.

        A real broadcaster counter reset looks identical in ONE packet, so the
        jump is not refused outright -- it has to be corroborated by
        SEQ_RESYNC packets agreeing on the new base before we move.  One bad
        bit cannot do that; a genuine reset does it in three MPUs.
        """
        hi = self.hi.get(pid)
        if hi is None:
            return True
        d = (seq - hi) & 0xFFFFFFFF          # 32-bit field: read it signed
        if d > 0x7FFFFFFF:
            d -= 1 << 32
        if -self.SEQ_GUARD <= d <= self.SEQ_GUARD:
            self.resync.pop(pid, None)
            return True
        base, n = self.resync.get(pid, (None, 0))
        if base is None or abs(seq - base) > self.SEQ_GUARD:
            self.resync[pid] = (seq, 1)
            self.stats["seq_wild"] += 1
            return False
        n += 1
        if n < self.SEQ_RESYNC:
            self.resync[pid] = (seq, n)
            self.stats["seq_wild"] += 1
            return False
        self.resync.pop(pid, None)
        self.stats["seq_resync"] += 1
        log(f"  pid {pid}: sequence resync {hi} -> {seq} "
            f"({self.SEQ_RESYNC} packets agreed)")
        self.hi[pid] = seq - 1
        return True

    def feed(self, p):
        r = O7.mmtp_parse(p, self.ver_ext)
        if r is None:
            self.stats["bad"] += 1
            return []
        self.stats["mmtp"] += 1
        pid, psn = r["pid"], r["psn"]
        prev = self.psn_last.get(pid)
        if prev is not None and psn != (prev + 1) & 0xFFFFFFFF:
            self.stats["psn_break"] += 1
            # The mask already forces a non-negative result, so the old
            # `max(0, ...)` guard was dead code and a REORDERED packet (psn
            # behind prev) was counted as ~2^32 losses -- one late datagram
            # reported 4,294,971,143 lost in a run that saw 5,411 packets.
            # Read the gap as signed: forward is loss, backward is reorder.
            d = (psn - prev - 1) & 0xFFFFFFFF
            if d > 0x7FFFFFFF:
                self.stats["psn_reorder"] += 1
            else:
                self.stats["psn_lost"] += d
        self.psn_last[pid] = psn
        if r["typ"] == 2:
            self.stats["signalling"] += 1
            return []
        if r["typ"] != 0:
            return []
        m = O7.mpu_payload(r["body"])
        if m is None:
            return []
        if not self._seq_ok(pid, m["seq"]):
            return []
        k = (pid, m["seq"])
        if k in self.done:
            self.stats["late_fragment"] += 1
            return []
        e = self.open.get(k)
        if e is None:
            e = self.open[k] = dict(
                parts=collections.defaultdict(list), timed=m["timed"],
                sizes=None, want=0, got=0, born=time.time())
        for du in m["dus"]:
            e["parts"][m["ft"]].append((psn, m["fi"], du))
            if m["ft"] == 2 and len(du) > O7.DU_HDR:
                e["got"] += len(du) - O7.DU_HDR
        if m["ft"] == 1 and e["sizes"] is None:
            mfmd = b"".join(d for _, _, d in sorted(e["parts"][1]))
            mi = O7.moof_info(mfmd)
            if mi["sizes"] and mi["mdat_declared"]:
                e["sizes"] = mi["sizes"]
                e["want"] = sum(mi["sizes"]) + mi["hint_size"] * mi["hint_n"]
        self.hi[pid] = max(self.hi.get(pid, -1), m["seq"])
        out = []
        if e["sizes"] is not None and e["got"] >= e["want"]:
            rec = self._build(pid, m["seq"], e)
            if rec is not None:
                # An MPU's packets all precede the next MPU's, so every still
                # open MPU older than this one has stopped growing and can be
                # judged NOW.  Repairs are emitted here, ahead of the completed
                # MPU, because a fragment handed over after its own successor
                # is out of order and the Transport drops it -- which is how a
                # repair turns into a silent no-op.
                out += self._repair_older(pid, m["seq"])
                out.append(rec)
                self.done.add(k)
                del self.open[k]
                self.stats["complete"] += 1
        out += self._reap(pid)
        return out

    def _build(self, pid, seq, e, partial=False):
        fl = O7.MmtpFlow("stream")
        fl.mpu[(pid, seq)] = dict(parts=e["parts"], timed=e["timed"])
        m = fl.build()[0]
        dec, got = m["mdat_declared"], m["mdat_got"]
        if dec is not None and m["n_samples"] and got == dec \
                and m["media_got"] == m["media_want"] and not m["n_short"]:
            m["truncate_to"] = None
            return m
        if not partial or dec is None or not m["n_samples"]:
            return None
        # `build` lays every sample out at its declared size and zero-fills the
        # holes, so the box offsets are still the transmitted ones and a cut at
        # the FIRST short sample leaves only bytes that actually arrived.
        keep = (m["short_first"] or 1) - 1
        if keep <= 0:
            return None
        m["truncate_to"] = keep
        return m

    def _repair_older(self, pid, seq):
        """Emit what arrived of the MPUs this one has overtaken.

        Whole-MPU loss is the receiver's coarsest failure: one bad FEC Block
        costs 2.002 s of video even when 119 of 120 frames are sitting in the
        buffer intact.  Everything up to the first short sample is playable and
        is handed over, labelled -- nothing is patched silently (M7's rule).
        """
        out = []
        for k in sorted(k for k in self.open
                        if k[0] == pid and k[1] < seq):
            e = self.open.pop(k)
            self.done.add(k)
            rec = self._build(pid, k[1], e, partial=self.repair)
            if rec is None:
                self.stats["lost"] += 1
                if e["sizes"] is None:
                    self.stats["lost_no_moof"] += 1
                else:
                    self.stats["lost_bytes"] += max(0, e["want"] - e["got"])
                continue
            self.stats["repaired"] += 1
            self.stats["repaired_samples"] += rec["truncate_to"]
            self.stats["repaired_dropped"] += (rec["n_samples"]
                                               - rec["truncate_to"])
            out.append(rec)
        return out

    def _reap(self, pid):
        """An MPU that newer MPUs have overtaken is never going to complete.

        Dropped and COUNTED, because "the picture froze for two seconds" has to
        show up as a number somewhere or the report is decoration.
        """
        cut = self.hi[pid] - self.keep
        # Reaching here means `keep` successors have arrived and none of them
        # COMPLETED, so nothing newer has been emitted and these repairs are
        # still in order.  Whatever cannot be repaired is dropped and COUNTED,
        # because "the picture froze for two seconds" has to show up as a
        # number somewhere or the report is decoration.
        out = self._repair_older(pid, cut)
        if len(self.done) > 4096:
            self.done = {k for k in self.done if k[1] >= cut - 64}
        return out


class RouteStream:
    """One ROUTE/LCT flow -> complete DASH media segments, the moment they
    complete.  (E47 -- STAGED for 132.1; inert unless Transport is built with
    `route=` naming this flow.)

    The batch shape of this transport is `m7_objects.RouteFlow` (feed all,
    assemble at the end) and `tools/atsc3_route.assemble_representation`
    (init + complete media objects in TOI order).  This is the same
    arithmetic run incrementally, and the E47 gate holds the two to BYTE
    IDENTITY on the banked 132.1 datagrams -- there is no second transport
    to drift, only a second schedule.

    Discipline mirrors MpuStreamer one-for-one where the problems are the
    same shape:
      * COMPLETE OR NOTHING.  A DASH media segment with a hole has zero
        bytes where its `styp` belongs (measured, M14) and poisons the
        concatenation, so partial objects are dropped whole and counted --
        drop-a-SEGMENT-not-a-sample, the lane freezes for 2.002 s instead
        of feeding the player garbage.
      * EMIT IN TOI ORDER, and never backwards.  The ROUTE fileTemplate
        maps DASH $Number$ -> LCT TOI, so ascending TOI IS presentation
        order (M14's negative control: shuffled TOIs = 360 non-monotonic
        DTS errors).
      * BOUND EVERY FIELD READ FROM AIR (the 8/07 header law).  `tol`
        caps object memory; an object whose declared length exceeds
        OBJ_MAX is refused outright, not allocated.
      * The SLS bundle (tsi 0) is parsed with `tools/atsc3_route.parse_sls`
        -- the M14 parser, imported, not re-derived -- to bind each tsi to
        its Representation (handler, kind, lang, duration).
    """

    INIT_TOI = 0xFFFFFFFF
    OBJ_MAX = 16 << 20            # no 2.002 s segment is 16 MB; a "tol" that
    #                               says otherwise is a corrupt field, refused
    PENDING_MAX = 16              # complete segments held while the SLS/init
    #                               are still awaited (carousel period ~2 s)
    PID_BASE = 100                # lane pid = PID_BASE + tsi: distinct from
    #                               the MMTP pids (12..15) by two decades

    _HANDLER = {"video": (b"vide", "route_video"),
                "audio": (b"soun", "route_audio"),
                "text": (b"subt", "route_subs")}

    def __init__(self, flow):
        self.flow = flow                    # "dst_ip:port", for the record
        self.obj = {}                       # (tsi, toi) -> open object
        self.done = set()                   # (tsi, toi) already emitted/lost
        self.init = {}                      # tsi -> init-segment bytes
        self.reps = None                    # tsi -> rep record, from the SLS
        self.service_id = None
        self.next_toi = {}                  # tsi -> next TOI to emit
        self.pending = {}                   # tsi -> {toi: bytes} pre-SLS
        self.stats = collections.Counter()
        self.info_ver = 0                   # bumped when the SLS (re)parses

    # -- LCT in -------------------------------------------------------------
    def feed(self, pay):
        """One UDP payload -> list of lane-segment dicts (usually empty)."""
        r = O7.lct_parse(pay)
        if r is None:
            self.stats["lct_bad"] += 1
            return []
        self.stats["lct"] += 1
        tsi, toi = r["tsi"], r["toi"]
        if tsi == 0 and toi == 0:           # FDT carousel; TOI==$Number$
            self.stats["fdt"] += 1          # makes it redundant here
            return []
        k = (tsi, toi)
        if k in self.done:
            self.stats["late_packet"] += 1
            return []
        o = self.obj.get(k)
        if o is None:
            o = self.obj[k] = dict(tol=None, frags={}, got=0, hi=0)
        if r["tol"] is not None:
            if r["tol"] > self.OBJ_MAX:     # a field read from air, bounded
                self.stats["tol_wild"] += 1
                self.done.add(k)
                del self.obj[k]
                return []
            o["tol"] = r["tol"]
        p = r["payload"]
        if p and r["start"] not in o["frags"]:
            o["frags"][r["start"]] = p
            o["got"] += len(p)
            o["hi"] = max(o["hi"], r["start"] + len(p))
            if o["got"] > self.OBJ_MAX:     # frags without a believable tol
                self.stats["obj_overrun"] += 1
                self.done.add(k)
                del self.obj[k]
                return []
        if o["tol"] is None or o["got"] < o["tol"] or o["hi"] != o["tol"]:
            return []
        data = self._assemble(o)
        self.done.add(k)
        del self.obj[k]
        if data is None:                    # byte count closed over a hole
            self.stats["obj_gapped"] += 1   # (overlapping frags): not usable
            return []
        return self._object_done(tsi, toi, data)

    def _assemble(self, o):
        """Contiguous concatenation, or None.  Same test RouteFlow.assemble
        applies: length equal AND no gaps -- `complete`, not `plausible`."""
        buf = bytearray()
        off = 0
        for st in sorted(o["frags"]):
            d = o["frags"][st]
            if st > off:
                return None
            if st < off:                    # overlap: keep the earlier bytes
                d = d[off - st:]
            buf += d
            off += len(d)
        return bytes(buf) if off == o["tol"] else None

    # -- objects out --------------------------------------------------------
    def _object_done(self, tsi, toi, data):
        if toi == self.INIT_TOI:
            if self.init.get(tsi) != data:
                self.init[tsi] = data
                self.stats["init_seg"] += 1
            return self._drain(tsi)
        if tsi == 0:
            return self._sls(data)
        self.pending.setdefault(tsi, {})[toi] = data
        self.stats["media_obj"] += 1
        return self._drain(tsi)

    def _sls(self, data):
        """The SLS multipart bundle: MPD + S-TSID, joined by the M14 parser."""
        try:
            sys.path.insert(0, os.path.join(os.path.dirname(HERE), "tools"))
            import atsc3_route as RT
            sid, reps = RT.parse_sls(data)
        except Exception as e:                                 # noqa: BLE001
            self.stats["sls_bad"] += 1
            log(f"  ROUTE {self.flow}: SLS bundle failed to parse "
                f"({type(e).__name__}: {e})")
            return []
        table = {}
        for r in reps:
            if r.get("tsi") is None:
                continue
            hd, kind = self._HANDLER.get(r.get("content_type"),
                                         (b"text", "route_data"))
            dur = 180180
            if r.get("duration") and r.get("timescale"):
                dur = int(round(float(r["duration"]) / float(r["timescale"])
                                * 90000.0))
            table[r["tsi"]] = dict(rep_id=r.get("rep_id"), handler=hd,
                                   kind=kind, dur=dur, lang=r.get("lang"),
                                   codecs=r.get("codecs"),
                                   content_type=r.get("content_type"))
        if table and table != (self.reps or {}):
            self.reps = table
            self.service_id = sid
            self.info_ver += 1
            self.stats["sls"] += 1
            log(f"  ROUTE {self.flow}: SLS parsed -- serviceId {sid}, "
                f"{len(table)} Representations "
                f"({', '.join(str(t) + ':' + v['kind'] for t, v in sorted(table.items()))})")
        out = []
        for tsi in sorted(self.pending):
            out += self._drain(tsi)
        return out

    def _drain(self, tsi):
        """Emit this tsi's complete segments in ascending TOI order.

        Requires the SLS row AND the init segment; until both exist the
        segments wait, bounded by PENDING_MAX (a carousel period is ~2 s,
        so a healthy flow never accumulates more than a couple)."""
        pend = self.pending.get(tsi)
        if not pend:
            return []
        rep = (self.reps or {}).get(tsi)
        if rep is None or tsi not in self.init:
            while len(pend) > self.PENDING_MAX:
                drop = min(pend)
                del pend[drop]
                self.stats["dropped_waiting"] += 1
            return []
        out = []
        nt = self.next_toi.get(tsi)
        for toi in sorted(pend):
            data = pend[toi]
            if nt is not None and toi < nt:
                self.stats["out_of_order"] += 1
                continue
            if nt is not None and toi > nt:
                # LEAVE A HOLE THE SIZE OF WHAT IS MISSING (the E20 law,
                # ROUTE edition): the lane's .idx records the TOI, so the
                # viewer freezes across the gap instead of shifting.
                self.stats["gap"] += 1
                self.stats["gap_segs"] += toi - nt
            # everything older and still open has been overtaken: it can
            # never complete (the carousel does not revisit media TOIs)
            for kk in [kk for kk in self.obj
                       if kk[0] == tsi and kk[1] < toi]:
                del self.obj[kk]
                self.done.add(kk)
                self.stats["lost"] += 1
            out.append(dict(seq=toi, bytes=data, dur=rep["dur"],
                            pid=self.PID_BASE + tsi, handler=rep["handler"],
                            init=self.init[tsi], transport="route",
                            lane_kind=rep["kind"], samples=None, kept=None))
            self.stats["segment"] += 1
            nt = toi + 1
        self.next_toi[tsi] = nt
        pend.clear()
        if len(self.done) > 4096:
            keep = min(self.next_toi.values(), default=0) - 64
            self.done = {kk for kk in self.done
                         if kk[1] >= keep or kk[1] == self.INIT_TOI}
        return out

    def route_info(self):
        """The viewer-facing service description (live_route.json)."""
        if self.reps is None:
            return None
        reps = []
        for t, v in sorted(self.reps.items()):
            r = {k: (x.decode("ascii") if isinstance(x, bytes) else x)
                 for k, x in v.items()}
            r.update(tsi=t, pid=self.PID_BASE + t)
            reps.append(r)
        return dict(ver=self.info_ver, flow=self.flow,
                    service_id=self.service_id, reps=reps)


# ---------------------------------------------------------------------------
# transport: ALP -> IP -> flow classification -> MPU -> segments
# ---------------------------------------------------------------------------

class Transport:
    """Baseband stream in, retimed ISOBMFF segments out.

    Flow classification is the same gate `m7_objects.main` applies -- the MPU
    payload's own `length` field against the bytes that actually remain -- run
    on the first `probe` datagrams of each flow instead of on the whole file.
    Nothing is hard-coded: not the address, not the port, not the packet_id.
    """

    # ISOBMFF handler types we know how to carry.  `vide` alone reproduces
    # the old video-only behaviour exactly.
    HANDLERS = (b"vide", b"soun", b"subt", b"text")

    def __init__(self, probe=64, want=(b"vide",), pid=None, dg_fh=None,
                 repair=True, route=None):
        self.repair = repair
        self.walker = StreamAlpWalker()
        self.ip = R7.IpReasm()
        self.probe = probe
        self.want = (want,) if isinstance(want, bytes) else tuple(want)
        self.force_pid = pid
        self.dg_fh = dg_fh
        self.pending = collections.defaultdict(list)
        self.kind = {}
        self.mmtp = {}
        # ROUTE (E47, staged): `route` is an explicit allow-list of
        # "dst_ip:port" flows to carry as ROUTE/DASH lanes.  Explicit
        # because three of the four LCT flows on this multiplex are
        # Widevine-locked services (law: untouched); only a flow the
        # operator NAMES is ever reassembled.  None = today's behaviour,
        # byte for byte: LCT flows classify "other" and are dropped.
        self.route_want = set(route) if route else None
        self.routes = {}          # flow key -> RouteStream
        self.stats = collections.Counter()
        self.assets = {}          # pid -> {init, handler, frag, base_time, ...}

    @property
    def init(self):
        """The video asset's init segment.

        Kept as a property so the single-asset player path reads the same as
        it always did, rather than every caller having to learn about the
        asset table.
        """
        for st in self.assets.values():
            if st["handler"] == b"vide":
                return st["init"]
        return next(iter(self.assets.values()))["init"] if self.assets else None

    @property
    def init_pid(self):
        for pid, st in self.assets.items():
            if st["handler"] == b"vide":
                return pid
        return next(iter(self.assets), None)

    def feed(self, chunk, bounds):
        segs = []
        for t, p in self.walker.feed(chunk, bounds):
            if t != 0:
                self.stats[f"ptype{t}"] += 1
                continue
            for src, dst, sp, dp, pay in self.ip.feed(p):
                self.stats["datagram"] += 1
                if self.dg_fh is not None:
                    self.dg_fh.write(R7.REC.pack(src, dst, sp, dp, len(pay))
                                     + pay)
                segs += self._route((dst, dp), pay)
        return segs

    def _route(self, key, pay):
        k = self.kind.get(key)
        if k is None:
            buf = self.pending[key]
            buf.append(pay)
            if len(buf) < self.probe:
                return []
            n_lct = sum(1 for p in buf if O7.lct_parse(p) is not None)
            sc = {v: O7.mmtp_score(buf, v) for v in (0, 2, 4)}
            best_v, (ok, tot) = max(sc.items(), key=lambda kv: kv[1][0])
            is_mmtp = tot > 0 and ok >= 0.9 * tot and ok > n_lct
            # ROUTE only for a flow the operator NAMED (see __init__), and
            # only if it actually LCT-parses -- the same referee the m13
            # inventory used (the 239.255.32.100:9000 flow parses 0/200
            # and stays "other" even if someone names it).
            # the IP layer keys flows by RAW 4-byte address; the operator
            # names them dotted -- normalise before comparing
            dst = key[0]
            if isinstance(dst, (bytes, bytearray)):
                dst = ".".join(map(str, dst))
            key_s = f"{dst}:{key[1]}"
            is_route = (not is_mmtp and self.route_want is not None
                        and key_s in self.route_want
                        and n_lct >= max(1, len(buf) // 2))
            k = self.kind[key] = (("mmtp", best_v) if is_mmtp else
                                  ("route", 0) if is_route else ("other", 0))
            self.stats["flow_" + k[0]] += 1
            log(f"  flow gate: {len(buf)} datagrams -> "
                f"{'MMTP (ver_ext=%d, %d/%d length agreement)' % (best_v, ok, tot) if is_mmtp else ('ROUTE/LCT (%d/%d parse) -- carrying as DASH lanes' % (n_lct, len(buf))) if is_route else 'not MMTP (LCT hits %d)' % n_lct}")
            if is_mmtp:
                self.mmtp[key] = MpuStreamer(best_v, repair=self.repair)
            elif is_route:
                self.routes[key] = RouteStream(key_s)
            out = []
            for p in buf:
                out += self._deliver(key, p)
            del self.pending[key]
            return out
        return self._deliver(key, pay)

    def _deliver(self, key, pay):
        k0 = self.kind[key][0]
        if k0 == "route":
            # RouteStream segments are already complete, self-timed DASH
            # fragments: no retime, no _segment -- straight to the lanes.
            return self.routes[key].feed(pay)
        if k0 != "mmtp":
            return []
        out = []
        for m in self.mmtp[key].feed(pay):
            s = self._segment(m)
            if s is not None:
                out.append(s)
        return out

    # E74: an MPU's own trun is an UNPROTECTED FIELD READ FROM THE AIR, and
    # this receiver bounds every one of those. Measured 8/10: a damaged MPU
    # (93 of its 120 frames) declared a media trun summing to exactly
    # 4.0000 MPU. Trusted, it moved the lane clock 4 slots for itself and
    # then priced the 5-MPU hole behind it at 4 slots each -- +18.0000 slots
    # injected, permanently, into every later fragment.
    DUR_TOL = 0.25        # a whole MPU may honestly vary by this much
    DUR_HIST = 16         # rolling window the nominal is taken from

    def _nominal(self, st):
        """The asset's own MPU cadence: the MEDIAN of recent durations.

        Median, not last-seen and not mean, so a single wild declaration
        cannot become the ruler -- which is exactly how E74 happened. Needs
        3 samples before it will speak; before that the caller falls back to
        the previous behaviour, and in practice a lane delivers hundreds of
        clean MPUs before its first hole (measured: 8 banked lanes, zero
        discontinuities)."""
        h = st.get("dur_hist") or []
        if len(h) < 3:
            return None
        # INTEGER TICKS, ALWAYS. statistics.median() returns a FLOAT on an
        # even-length window even when every sample is identical -- 16
        # copies of 180180 give 180180.0 -- and DUR_HIST is 16. The clamp
        # path happened to coerce (`dur = int(nom)`), but the GAP-PRICING
        # path did not, so `st["base_time"] += missing * nom` turned the
        # lane clock into a float and the next retime died on
        #     struct.error: required argument is not an integer
        # inside struct.pack_into(">Q", ...). It needed a gap AFTER the
        # 16-sample history filled, which is why a patch that gated clean
        # crash-looped the live chain ~2 h later (E81, 8/11: rc=1 every
        # ~2 min, 23 restarts, a lane roll each time).
        # A tick count is an integer by construction; round rather than
        # truncate so the cadence cannot drift low.
        return int(round(statistics.median(h)))

    def _bound(self, st, dur, seq, pid):
        """Clamp one MPU's declared duration to the asset's cadence.

        Returns the duration to use and RECORDS the observation. A clamp is
        logged, never silent: a damaged sample table is a fact about the air
        worth seeing, and absorbing it quietly is how a 36-second timeline
        error hid behind healthy-looking FEC for seven hours."""
        nom = self._nominal(st)
        if nom and not (nom * (1 - self.DUR_TOL) <= dur
                        <= nom * (1 + self.DUR_TOL)):
            self.stats["dur_clamped"] += 1
            self.stats["dur_clamped_ticks"] += dur - nom
            log(f"  pid {pid}: MPU {seq} declares {dur} ticks "
                f"({dur / nom:.2f}x the {nom:.0f}-tick cadence) -- clamping. "
                f"A damaged sample table must not move the lane clock.")
            dur = int(nom)
        h = st.setdefault("dur_hist", [])
        h.append(dur)
        del h[:-self.DUR_HIST]
        return dur

    def _segment(self, m):
        """One complete MPU -> one retimed (moof + mdat), or None if it is not
        an asset we are carrying.

        See `_bound`: an MPU's self-declared duration is a field read from
        the air, and this receiver bounds every one of those.

        MULTI-ASSET. This used to lock exactly ONE pid and skip the rest, so
        live was video-only: the AC-4 audio and the TTML captions were being
        demodulated, reassembled, and then thrown away at this line. Each
        asset now keeps its OWN fragment counter and media clock, because
        `retime` rewrites both per fragment and sharing them across assets
        would interleave two timelines into one and desynchronise everything.
        """
        pid = m["pid"]
        st = self.assets.get(pid)
        if st is None:
            if not m["meta"]:
                return None
            h = moov_handler(m["meta"])
            if self.force_pid is not None:
                if pid != self.force_pid:
                    return None
            elif h not in self.want:
                self.stats["skip_" + (h or b"none").decode("ascii",
                                                           "replace")] += 1
                return None
            st = self.assets[pid] = dict(init=m["meta"], handler=h, frag=0,
                                         base_time=0, next_seq=None,
                                         last_dur=0, dur_hist=[])
            log(f"  asset locked: pid {pid}, handler "
                f"{(h or b'?').decode('ascii', 'replace')}, "
                f"MPU metadata {len(m['meta'])} bytes")
        if st["next_seq"] is not None and m["seq"] < st["next_seq"]:
            self.stats["out_of_order"] += 1
            return None
        # `next_seq` is the second monotone latch on the unprotected sequence
        # number (MpuStreamer.hi is the first).  The streamer's SEQ_GUARD
        # should mean nothing wild ever reaches here; say so out loud rather
        # than trusting it silently, because when this latch does move too far
        # the asset goes quiet permanently and NOTHING logs it.
        if st["next_seq"] is not None:
            ahead = m["seq"] - st["next_seq"]
            if ahead > MpuStreamer.SEQ_GUARD:
                self.stats["seq_wild_transport"] += 1
                log(f"  pid {pid}: refusing a {ahead}-MPU forward jump "
                    f"(next_seq {st['next_seq']}, saw {m['seq']})")
                return None
        if st["next_seq"] is not None and m["seq"] > st["next_seq"]:
            missing = m["seq"] - st["next_seq"]
            self.stats["gap"] += 1
            self.stats["gap_mpus"] += missing
            # LEAVE A HOLE THE SIZE OF WHAT IS MISSING.
            #
            # A damaged-but-partly-usable MPU already holds the clock true
            # (`repair_unusable` below advances base_time across it). An MPU
            # that never arrived at all did not, so the next fragment was
            # placed 2.002 s early and EVERYTHING after it inherited the
            # shift -- per lost MPU, per lane, permanently. That is why E20's
            # lanes ended 37 fragments apart over two hours: not drift, an
            # accumulated sum of un-accounted holes.
            #
            # The missing MPUs' own durations are gone with them, so use the
            # asset's NOMINAL cadence. The broadcaster is constant-rate
            # (E19: 180180 ticks every time), which makes that exact in
            # practice and honest when it is not: the picture freezes for the
            # length of what was lost instead of sliding early by it.
            #
            # NOT `last_dur` (E74, 8/10). One damaged MPU declared a media
            # trun summing to exactly 4.0000 MPU; it became `last_dur` and
            # then PRICED THE HOLE THAT FOLLOWED IT -- 5 missing x 4 = 20
            # slots, +1 for itself, so 6 sequence numbers of air advanced the
            # lane clock by 24. The +18.0000 slots (36.036 s) rode every
            # later fragment, the viewer stamped video that far off its own
            # audio, and VLC discarded 82.7 % of pictures. A single bad
            # measurement must never become the ruler: the error is
            # multiplied by the size of every hole it prices, so a bigger
            # hole would have injected minutes.
            nom = self._nominal(st) or st.get("last_dur")
            if nom:
                st["base_time"] += missing * nom
                self.stats["gap_ticks"] += missing * nom
        keep = m.get("truncate_to")
        if keep:
            # (durations from both retime paths below pass through _bound)
            # Retime FIRST, then cut: `retime` reads the untouched sample table,
            # so `dur` is the MPU's FULL 2.002 s.  The timeline therefore keeps
            # advancing by the air's own clock and the lost tail shows as a
            # short freeze in place, not as everything after it sliding early.
            seg, dur = PL.retime(m["body"], st["frag"] + 1, st["base_time"])
            dur = self._bound(st, dur, m["seq"], pid)
            r = PL.trun_keep(seg, keep)
            if r is None:
                self.stats["repair_unusable"] += 1
                st["next_seq"] = m["seq"] + 1
                st["base_time"] += dur     # hold the clock true across the hole
                st["last_dur"] = dur
                return None
            seg = r[0]
            self.stats["segment_truncated"] += 1
            self.stats["frames_recovered"] += r[1]
        else:
            seg, dur = PL.retime(m["body"], st["frag"] + 1, st["base_time"])
            dur = self._bound(st, dur, m["seq"], pid)
        st["next_seq"] = m["seq"] + 1
        st["frag"] += 1
        st["base_time"] += dur
        st["last_dur"] = dur           # what a whole MPU is worth on this asset
        self.stats["segment"] += 1
        return dict(seq=m["seq"], bytes=seg, dur=dur, samples=m["n_samples"],
                    kept=keep or m["n_samples"],
                    pid=pid, handler=st["handler"], init=st["init"])

    def flush(self):
        segs = []
        for t, p in self.walker.pump(final=True):
            if t != 0:
                continue
            for src, dst, sp, dp, pay in self.ip.feed(p):
                self.stats["datagram"] += 1
                if self.dg_fh is not None:
                    self.dg_fh.write(R7.REC.pack(src, dst, sp, dp, len(pay))
                                     + pay)
                segs += self._route((dst, dp), pay)
        return segs


# ---------------------------------------------------------------------------
# 1/2 -- the continuous front end
# ---------------------------------------------------------------------------

def fine_timing_par(y, ex, span=FT_SPAN, centre=None, fast=False, plan=None):
    """`m9_fast.FrameDecoder.fine_timing`, verbatim, including the tie-break.

    Strictly-greater wins, which is `m6_cells.fine_timing`'s rule; a `>=` here
    would pick a different t0 on a tie and the whole decode would be a
    different decode.

    E82: `plan` routes an LDM/CTI multiplex to the general-geometry metric
    (`m44_ldm.PreambleCoh` -- the PREAMBLE's pilot coherence, whose geometry
    A/322 7.2.5.1 fixes rather than signals, so it works before L1 exists).
    `plan=None` is RF33's path, byte for byte: m9_accel's data-symbol
    coherence over the module constants, untouched.
    """
    centre = C.BOOTSTRAP if centre is None else centre
    if plan is not None:
        import m44_ldm as M44
        nw = (getattr(ex, "_max_workers", 1) if ex else 1)
        dc = getattr(plan, "_dcoh", None)
        if dc is None:
            dc = plan._dcoh = M44.DataCoh(plan)
        r = dc.scan(y, centre, span, workers=nw)
        if r is not None:
            return r
        # too near the edge of the window for a data symbol -- fall back to
        # the Preamble metric, which needs far fewer samples
        pc = getattr(plan, "_coh", None)
        if pc is None:
            pc = plan._coh = M44.PreambleCoh(plan.g)
        return pc.scan(y, centre, span, workers=nw)
    cands = list(range(centre - span, centre + span + 1))
    # E52: demod_coh is demod_data's coherence with the lines copied -- the
    # same value bit for bit, minus the equalisation nobody here reads.  The
    # batched form (one FFT for all 41 candidates) is FAST MODE ONLY: its
    # coherences drift in the last ulps, so it is gated on the t0 decision
    # (gate_e52) and on the e2e datagram sha, not assumed.
    import m9_accel as A9
    if fast:
        res = A9.demod_coh_batch(y, cands, 1, False)
    elif ex is None:
        res = [A9.demod_coh(y, t, 1, False) for t in cands]
    else:
        res = list(ex.map(lambda t: A9.demod_coh(y, t, 1, False), cands))
    best = None
    for t, coh in zip(cands, res):
        if best is None or coh > best[1]:
            best = (t, coh)
    return best


class WindowNotch:
    """E55: the notch as overlap-save FFT filtering of the DECODE WINDOW.

    Two ideas, and the second is the bigger one:

      1. Overlap-save instead of the time-domain FIR: FFT a segment,
         multiply by `fft(taps, N)` -- the FIR's OWN frequency response,
         so the mask matches the reference by construction, not by a
         designed approximation -- and IFFT.  Each valid output is the
         same 257-summand convolution the FIR computes, re-associated by
         the FFT; the values land ~4e-15 off the FIR's, so this path is
         DECISION-GATED (fast mode only; the FIR remains the reference)
         exactly like the E52/E53 fast-mode concessions.  Referee:
         gate_e55 -- response match, decode equivalence, bookkeeping.
      2. Filter the WINDOW, not the stream.  A Frame is 1,708,032 samples
         but the decoder reads only FRAME_WINDOW = 354,344 of them (the
         first subframe; PLP0+PLP16): the continuous notch spends 79.5%
         of its work on samples nobody decodes.  Windows of consecutive
         Frames never overlap (they advance by FRAME_SAMPLES), so each
         window is filtered independently from real tape history --
         which also makes this path STATELESS across re-acquisitions:
         there is no primed/hist state for a cross-thread reacquire()
         to reset out from under us (the class of race that produced
         the live1b TypeError in the de-rotation ramp).

    Alignment (the trim/group-delay accounting): the continuous FIR path
    filters the stream and then eats NOTCH_TAPS//2 = 128 samples once per
    lock, so its tape holds zero-phase-aligned samples.  Here the tape
    holds RAW samples and the alignment is done per window: output t is
    sum_j taps[j] * raw[t + 128 - j], i.e. the SAME zero-phase sample the
    FIR path's tape holds at index t.  Sample count: exactly FRAME_WINDOW
    out per window; the tape append count equals the input count (trim=0
    on this path) -- samples in == samples out, the capture-integrity law,
    checked by gate_e55 leg 3.

    Measured (TR 2990WX): the 1,708,032-sample continuous FIR is 133 ms
    single-threaded (446 ms on the 1600X, a 247 ms budget); this path is
    ~12 ms single-threaded, ~7 ms on the front-end pool, per Frame.
    """

    NFFT = 8192                       # V = NFFT-256 = 7936; ~45 seg/window

    def __init__(self, taps, nfft=NFFT):
        self.m = len(taps)            # 257
        self.n = int(nfft)
        self.v = self.n - self.m + 1  # valid outputs per segment
        self.H = np.fft.fft(taps, self.n)   # the FIR's response, exactly

    def apply(self, buf, ex=None):
        """buf: complex128 with `m-1` samples of REAL history in front.
        -> len(buf) - (m-1) filtered samples (complex128), each the full
        257-tap convolution -- same summands as the FIR, FFT-grouped."""
        m1 = self.m - 1
        n = len(buf) - m1
        if n <= 0:
            return np.zeros(0, np.complex128)
        nseg = -(-n // self.v)
        need = nseg * self.v + m1
        if len(buf) < need:
            # zero TAIL pad only: affects outputs beyond n, all discarded
            buf = np.concatenate((buf, np.zeros(need - len(buf), buf.dtype)))
        elif not buf.flags.c_contiguous:
            buf = np.ascontiguousarray(buf)
        st = buf.strides[0]
        W = np.lib.stride_tricks.as_strided(
            buf, shape=(nseg, self.n), strides=(self.v * st, st))
        out = np.empty(n, np.complex128)

        def work(r0, r1):
            Y = SFFT.fft(W[r0:r1], axis=1, workers=1)
            Y *= self.H
            z = SFFT.ifft(Y, axis=1, workers=1)
            a, b = r0 * self.v, min(r1 * self.v, n)
            out[a:b] = z[:, m1:].reshape(-1)[:b - a]

        if ex is None or nseg < 8:
            work(0, nseg)
        else:
            nw = max(1, getattr(ex, "_max_workers", 4))
            step = max(1, -(-nseg // nw))
            futs = [ex.submit(work, i, min(i + step, nseg))
                    for i in range(0, nseg, step)]
            for f in futs:
                f.result()
        return out


class FineTrack:
    """E53: fine timing as a TIMING LINE, not a 41-point search every Frame.

    What the full scan actually measures on real air (gate capture): t0
    advances by FRAME_SAMPLES + d with d jittering +-20 samples FRAME TO
    FRAME while its mean sits near +7 (4 ppm of sample clock).  The jitter
    is metric noise on a flat-topped coherence plateau -- the first +-1
    tracker DIVERGED because it chased it, and the decoder itself does not
    care: every t0 on the plateau decodes 74/74 (that is what a cyclic
    prefix is FOR).  The physical quantity is the clock-offset LINE, so
    estimate that instead:

      * dead-reckon t0 along the line (FRAME_SAMPLES + d per Frame, with
        the FRACTIONAL part of d carried, never dropped -- the E48 lesson);
      * every FULL_EVERY-th Frame, run the full +-FT_SPAN scan and update
        d from the measured advance since the last scan (EMA);
      * every dead-reckoned Frame still pays ONE FFT: the coherence at the
        predicted t0.  If it falls under COH_GUARD of the running average,
        the plateau moved -- full scan NOW.

    Cost: ~1.2 FFTs/Frame amortised vs 41.  Fast mode only; the exact path
    keeps HEAD's full scan.  The referee is not the trajectory (two points
    on a plateau are equally right) but the DECODED BYTES: gate_e52 decodes
    at tracked t0s against the exact chain at full-scan t0s, and the m11
    e2e gate holds the whole datagram file.
    """

    # Anchors every 4 Frames: the dead-reckon's deviation from the scan
    # grows with anchor spacing, and the byte-plateau measured on the gate
    # capture ends somewhere past ~90 samples (38/74 -> 0/74 FEC at ~100+;
    # E53's first build used 8 and walked off it).  At 4 the deviation
    # stays in the tens.  The coherence guard below is a WEAK guard --
    # measured: it never tripped while FEC fell to 0/74 -- so the decode
    # loop can also demand a full rescan directly (request_full), which
    # makes the LDPC syndrome the tracker's real backstop.
    FULL_EVERY = 4
    COH_GUARD = 0.6

    def __init__(self, fast, ema=0.5, plan=None):
        self.fast = fast
        self.plan = plan          # E82: None = RF33's constants, unchanged
        self.frame_samples = (C.FRAME_SAMPLES if plan is None
                              else plan.frame_samples)
        self.ema = ema
        self.d = None                 # per-Frame advance beyond FRAME_SAMPLES
        self.frac = 0.0               # fractional-sample carry of d
        self.coh_ema = None
        self.anchor_t0 = None         # last full-scan t0
        self.anchor_n = 0             # Frames since that scan
        self.force_full = False
        self.stats = collections.Counter()

    def request_full(self):
        """Downstream (FEC) demands a full rescan on the next Frame -- the
        syndrome is the tracker's real backstop; the coherence guard was
        measured blind while FEC fell to 0/74."""
        self.force_full = True

    def reset(self):
        self.d = None
        self.frac = 0.0
        self.coh_ema = None
        self.anchor_t0 = None
        self.anchor_n = 0
        self.force_full = False

    def _full(self, y, centre, ex, abs_off):
        t0, coh = fine_timing_par(y, ex, FT_SPAN, centre, fast=self.fast,
                                  plan=self.plan)
        # The anchor arithmetic runs on ABSOLUTE sample indices.  The first
        # build fed it window-relative t0s: the rate came out as
        # (20 - 9)/1 - FRAME_SAMPLES = -1,708,021 samples/Frame and the
        # cursor advanced +31 per "Frame" -- frames() then minted ~100k
        # windows per push and the front thread looked hung inside np.roll.
        # A track without a coordinate system is a random walk with
        # confidence.
        t0_abs = abs_off + t0
        # E55: snapshot -- reset() runs on the DECODE thread (via
        # reacquire()) and rebinds anchor_t0 to None; same race class as
        # the de-rotation ramp that killed live1b.
        anchor = self.anchor_t0
        if anchor is not None and self.anchor_n > 0:
            rate = (t0_abs - anchor) / self.anchor_n - self.frame_samples
            if abs(rate) <= FT_SPAN:      # a sane clock offset; else re-anchor
                self.d = rate if self.d is None \
                    else self.ema * self.d + (1.0 - self.ema) * rate
            else:
                self.stats["rate_wild"] += 1
                self.d = None
        self.anchor_t0 = t0_abs
        self.anchor_n = 0
        return t0, coh

    def search(self, y, centre, ex, abs_off=0):
        """-> (t0, coh) in the CALLER's coordinates; `abs_off + t0` must be
        the absolute sample index (abs_off=0 when y is the whole stream)."""
        self.anchor_n += 1
        if not (self.fast and self.d is not None):
            t0, coh = self._full(y, centre, ex, abs_off)
            self.stats["full_scan"] += 1
        elif self.force_full:
            self.force_full = False
            t0, coh = self._full(y, centre, ex, abs_off)
            self.stats["fec_rescan"] += 1
        elif self.anchor_n >= self.FULL_EVERY:
            t0, coh = self._full(y, centre, ex, abs_off)
            self.stats["full_rescan"] += 1
        else:
            t0 = centre
            coh = float(fine_timing_par(y, ex, 0, centre, fast=True,
                                        plan=self.plan)[1])
            if self.coh_ema is not None and coh < self.COH_GUARD * self.coh_ema:
                t0, coh = self._full(y, centre, ex, abs_off)
                self.stats["guard_rescan"] += 1
            else:
                self.stats["tracked"] += 1
        self.coh_ema = coh if self.coh_ema is None \
            else 0.9 * self.coh_ema + 0.1 * coh
        return t0, coh

    def drift(self):
        """Integer samples to add to FRAME_SAMPLES for the next prediction,
        with the fractional remainder CARRIED, not dropped.  0 on the exact
        path -- its cursor arithmetic stays HEAD's."""
        d = self.d       # E55: reset() (decode thread) may rebind d = None
        if not (self.fast and d is not None):
            return 0
        x = d + self.frac
        step = int(np.floor(x + 0.5))
        self.frac = x - step
        return step


class FrontEnd:
    """Raw capture-rate samples in, per-Frame decode windows out.

    ONE bootstrap search.  M9's design said it plainly: re-detecting per batch
    is 6 ms/Frame of pure waste and, worse, it re-anchors -- two adjacent
    batches can disagree about where the Frame starts.  Here the search runs
    once, fixes `f0` and the CFO, and after that t0 is a running quantity:
    `centre = t_prev + FRAME_SAMPLES`, refined by the same +-20 fine-timing
    scan the batch chain uses, so the tracking IS the reference's tracking.
    """

    def __init__(self, fs_in, ex=None, notch=True, fast=False, plan=None,
                 plan_cb=None):
        self.fs_in = fs_in
        self.fast = fast   # E52: batched fine timing (decision-gated)
        # E82: `plan` (an m44_ldm.LdmPlan) makes the Frame geometry a
        # PARAMETER instead of RF33's module constants.  plan=None keeps
        # every number below exactly what it has always been -- the RF33
        # byte-identity gate (gate_e82 leg 1) is what holds that claim.
        self.plan = plan
        # E82: on an unknown multiplex the Frame LENGTH is itself signalled,
        # so it cannot be a constructor argument -- `plan_cb` is called once,
        # on the first window after the bootstrap lock, and returns the plan.
        # A radio cannot be rewound, so this runs FORWARD off the tape rather
        # than re-reading; if L1 does not verify the cursor steps to the next
        # bootstrap (found in the tape, not assumed) and tries again.
        self.plan_cb = plan_cb
        self.plan_tries = 0
        self.frame_samples = (FRAME_SAMPLES if plan is None
                              else plan.frame_samples)
        self.frame_window = (FRAME_WINDOW if plan is None
                             else plan.frame_window)
        self.ft = FineTrack(fast, plan=plan)   # E53: timing-line tracked search
        self.notch_decided = not fast   # E53: measured notch decision
        self._ramp = None               # E53: cached de-rotation ramp
        self.ex = ex
        self.notch = notch
        self.rs = PolyResampler(fs_in, FS_POST)
        self.acq = []
        self.acq_n = 0
        self.gen = 0          # bumped by reacquire(); see _acquire()
        self.acq_need = int(0.30 * fs_in)
        self.ready = False
        self.pend = np.zeros(0, np.complex64)
        self.tape = SampleTape()
        self.nb = firwin(NOTCH_TAPS, NOTCH_CUT, fs=FS_POST)
        # E55: fast mode filters the decode WINDOW by overlap-save FFT
        # (WindowNotch) instead of running the FIR over the whole stream;
        # ATSC3_NOTCH_FFT=0 restores the E53 stream-FIR in fast mode (the
        # A/B lever gate_e55 uses).  The exact path is untouched.
        self.fft_notch = bool(fast) and \
            os.environ.get("ATSC3_NOTCH_FFT", "1") != "0"
        self.osn = WindowNotch(self.nb) if self.fft_notch else None
        self.hist = np.zeros(NOTCH_TAPS - 1, np.complex128)
        self.notch_primed = False
        self.n_der = 0
        # the FFT path's tape holds RAW samples (alignment is per window),
        # so it never eats the FIR group delay off the stream
        self.trim = 0 if self.fft_notch else NOTCH_TAPS // 2
        self.cfo = None
        self.f0 = None
        self.ps = None
        self.centre = C.BOOTSTRAP
        self.frame = 0
        self.stats = collections.Counter()
        self.tm = collections.Counter()

    # -- acquisition -------------------------------------------------------
    def _acquire(self):
        # `self.acq` is rebound to [] by reacquire(), which runs on the DECODE
        # thread (m11_watch bad_run path) while we are on the FRONT thread.
        # find_bootstraps() below takes milliseconds, so reading self.acq once
        # before it and again after it is a race with a wide window: the
        # second read used to raise "need at least one array to concatenate"
        # every time a re-acquire storm coincided with an acquisition, which
        # is precisely when the signal is bad.  Snapshot once, and carry a
        # generation so a reset that lands mid-search is honoured instead of
        # being overwritten with a result derived from discarded samples.
        acq, gen = self.acq, self.gen
        if not acq:
            return False
        buf = np.concatenate(acq)
        hits = MP.find_bootstraps(
            MP.resample_to(buf[:self.acq_need], self.fs_in, bs.FS))
        if not hits:
            return False
        if self.gen != gen:
            # reacquire() fired while we were searching; these samples are
            # the ones it just decided to throw away.
            self.stats["acquire_raced"] += 1
            return False
        h = hits[0]
        fs_post = h["fields"]["post_bootstrap_sample_rate_hz"]
        if abs(fs_post - FS_POST) > 1.0:
            raise RuntimeError(f"post-bootstrap rate {fs_post} != {FS_POST}")
        self.cfo = h["fine_cfo_hz"]
        # E82: the bootstrap NAMES the preamble structure (A/322 Table H.1.1
        # -> FFT / GI / pilot DX / L1-Basic Mode).  Carrying it forward means
        # the geometry probe never has to guess which row it is looking at --
        # guessing produced a silent wrong-FFT read on the first build.
        self.ps = h["fields"].get("preamble_structure")
        rel = int(round(h["position"] * self.fs_in / bs.FS))
        self.f0 = rel
        self.k = -2j * np.pi * self.cfo
        rest = buf[rel:]
        self.acq = []
        self.acq_n = 0
        self.ready = True
        self.stats["acquire"] += 1
        log(f"  bootstrap acquired: position {h['position']} @ "
            f"{bs.FS/1e6:g} Msps -> input sample {rel}, CFO "
            f"{self.cfo:+.2f} Hz, peak ratio {h.get('mean_peak_ratio', 0):.1f}")
        self._stream(rest)
        return True

    def reacquire(self):
        """Sample continuity was broken (an overflow, or the decode fell off).
        Everything downstream of `f0` is meaningless, so drop it and search
        again rather than pretend the tracking survived."""
        self.gen += 1         # invalidate any acquisition in flight
        self.ready = False
        self.acq = []
        self.acq_n = 0
        self.pend = np.zeros(0, np.complex64)
        self.tape = SampleTape()
        self.rs = PolyResampler(self.fs_in, FS_POST)
        self.hist = np.zeros(NOTCH_TAPS - 1, np.complex128)
        self.notch_primed = False
        self.n_der = 0
        self.notch = self.notch or self.fast   # E53: back on until re-measured
        self._ramp = None               # E53: CFO changes with the new lock
        self.notch_decided = not self.fast
        self.trim = 0 if self.fft_notch else NOTCH_TAPS // 2
        self.centre = C.BOOTSTRAP
        self.ft.reset()             # E53: the track died with the lock
        self.plan_tries = 0         # E82: a new lock gets a fresh probe budget
        self.stats["reacquire"] += 1

    # -- steady state ------------------------------------------------------
    def push(self, x):
        """x: complex64 at the capture rate."""
        if not self.ready:
            self.acq.append(x)
            self.acq_n += len(x)
            if self.acq_n >= self.acq_need:
                if not self._acquire():
                    self.acq = []
                    self.acq_n = 0
                    self.stats["acquire_miss"] += 1
            return
        self._stream(x)

    # E53: the notch's job is adjacent-channel/interferer rejection, and on
    # a clean channel it is 51 ms/Frame of work that changes no decoded
    # byte (measured: 74/74 FEC and BYTE-IDENTICAL datagrams with it off,
    # out-of-band power 24 dB down).  So in fast mode, MEASURE instead of
    # assume: one PSD ratio on the first post-acquisition block decides the
    # whole lock, re-decided at every re-acquisition.  The threshold sits
    # 7 dB above the clean-channel reading and well below an interferer
    # that would actually threaten 64QAM.  Exact mode always notches.
    NOTCH_OOB_THRESH = 0.02          # mean PSD >3.0 MHz / mean PSD <2.7 MHz

    # -- E82: runtime geometry acquisition ---------------------------------
    PROBE_WINDOW = 96_000        # bootstrap + a 2-symbol Preamble + slack
    PROBE_MAX = 24               # Frames tried before giving up on this lock

    def adopt_plan(self, plan):
        self.plan = plan
        self.frame_samples = plan.frame_samples
        self.frame_window = plan.frame_window
        self.ft.plan = plan
        self.ft.frame_samples = plan.frame_samples
        self.ft.reset()

    def _probe_plan(self):
        """-> True once the Frame geometry is known."""
        lo = self.centre - FT_SPAN
        if lo < self.tape.org:
            return False
        if self.tape.end < lo + self.PROBE_WINDOW:
            return False
        w = self.tape.take(lo, self.PROBE_WINDOW)
        plan = None
        try:
            # index 0 of `w` is absolute `lo`, so the predicted Frame start
            # sits at FT_SPAN inside it -- the same convention `frames()` uses
            plan = self.plan_cb(w, FT_SPAN, self.ps)
        except Exception as e:                                 # noqa: BLE001
            log(f"  geometry probe raised {type(e).__name__}: {e}")
        self.plan_tries += 1
        if plan is not None:
            self.adopt_plan(plan)
            log(f"  Frame geometry acquired from L1 after "
                f"{self.plan_tries} probe(s): {plan.frame_samples} samples "
                f"({plan.frame_sec * 1000:.3f} ms), window "
                f"{plan.frame_window}")
            return True
        # No L1 this Frame.  We do not know the Frame length yet, so the
        # cursor cannot be advanced by arithmetic -- find the NEXT bootstrap
        # in the tape and move there.  Costs one Frame per try and needs no
        # assumption about a multiplex we have not identified.
        nxt = None
        avail = self.tape.end - (lo + 1000)
        if avail < int(0.35 * FS_POST):
            self.plan_tries -= 1          # not a try: the tape is just short
            # A tape that grows with no Frame produced is the 34.9 GB failure
            # E53 recorded; the probe path must obey the same bound as the
            # decode path rather than being an exemption from it.
            if self.tape.end - self.tape.org > self.TAPE_MAX:
                self.stats["tape_overrun"] += 1
                log("  *** front-end tape full while probing the geometry -- "
                    "re-acquiring ***")
                self.reacquire()
            return False
        seek = self.tape.take(lo + 1000, min(avail, 3_000_000))
        try:
            hits = MP.find_bootstraps(MP.resample_to(seek, FS_POST, bs.FS))
            if hits:
                nxt = lo + 1000 + int(round(hits[0]["position"] * FS_POST
                                            / bs.FS))
        except Exception:                                      # noqa: BLE001
            nxt = None
        log(f"  geometry probe {self.plan_tries}: L1 did not verify on this "
            f"Frame -- stepping to the next bootstrap")
        if nxt is None or nxt <= lo:
            self.stats["plan_probe_stalled"] += 1
            if self.tape.end - self.tape.org > self.TAPE_MAX:
                self.stats["tape_overrun"] += 1
                self.reacquire()
            return False
        self.centre = nxt + C.BOOTSTRAP
        self.tape.trim(nxt - 8)
        self.stats["plan_probe_retry"] += 1
        if self.plan_tries >= self.PROBE_MAX:
            log(f"  *** L1 did not verify in {self.plan_tries} Frames -- "
                f"re-acquiring the bootstrap ***")
            self.reacquire()
        return False

    def _decide_notch(self, y):
        n = 1 << 18
        if len(y) < n:
            return                    # wait for a block big enough
        S = np.abs(np.fft.fft(y[:n])) ** 2
        f = np.fft.fftfreq(n, 1.0 / FS_POST)
        ratio = float(S[np.abs(f) > 3.0e6].mean()
                      / max(S[np.abs(f) < 2.7e6].mean(), 1e-30))
        self.notch_decided = True
        if os.environ.get("ATSC3_NOTCH") == "0":
            # E55: operator kill-switch. On the 1600X the FIR costs 446 ms
            # of a 247 ms budget -- when the notch is the difference between
            # "decodes with interference" and "cannot keep up at all", the
            # only way to know which side of that trade the channel is on
            # is to run without it and read the FEC.
            self.notch = False
            self.trim = 0
            log(f"  notch FORCED OFF (ATSC3_NOTCH=0): out-of-band {ratio:.4f}"
                f" of in-band power would have armed it")
        elif ratio < self.NOTCH_OOB_THRESH:
            self.notch = False
            self.trim = 0             # no FIR, no group delay to eat
            log(f"  notch BYPASSED for this lock: out-of-band {ratio:.4f} "
                f"of in-band power (< {self.NOTCH_OOB_THRESH}); "
                f"re-measured at every re-acquisition")
        else:
            log(f"  notch ACTIVE: out-of-band {ratio:.4f} of in-band power "
                f"(>= {self.NOTCH_OOB_THRESH})"
                + (", FFT window path" if self.fft_notch else ", stream FIR"))

    def _stream(self, x):
        mult = self.rs.block_multiple()
        self.pend = np.concatenate((self.pend, x)) if len(self.pend) else x
        n = (len(self.pend) // mult) * mult
        if n == 0:
            return
        blk, self.pend = self.pend[:n], self.pend[n:]
        pc = time.perf_counter
        t = pc()
        y = self.rs.feed(blk).astype(np.complex128)
        self.tm["resample"] += pc() - t
        if not len(y):
            return
        if self.fast and not self.notch_decided:
            self._decide_notch(y)
        base = self.n_der
        self.n_der += len(y)

        # Both of these are M9 primitives, used for M9's reason: NumPy and
        # SciPy release the GIL inside their kernels, and a live pipeline is
        # ONE stream, so the parallelism has to be inside the stage.  Run
        # single-threaded, the de-rotation and the notch alone are ~0.7x real
        # time and the decoder starves -- measured, on the first build of this
        # file, before this line existed.
        #
        # The de-rotation reference groups the exponent as
        # ((((-2j*pi)*cfo) * n) / fs_post); folding the division into the
        # scalar is a different rounding (M9 measured 1.9e-13).  `n` is the
        # sample index from f0, which is what makes the phase continuous
        # across blocks -- the whole point of carrying `base`.
        #
        # E53 fast mode: the CFO is FIXED for the whole lock, so the ramp
        # exp(k*i/FS) is one cached table and the per-chunk phasor is one
        # scalar exp: exp(k*(lo+base+i)) becomes exp(k*(lo+base)) * R[i].
        # Two multiplies per sample instead of a complex exp (~10x the
        # cost) per sample -- a 1-ulp re-association, decoded-bytes gated
        # like every other fast-mode concession.  Exact mode unchanged.
        if self.fast:
            need = len(y)
            # E55 CRASH FIX (live1b, "'NoneType' object is not
            # subscriptable" after 51 re-acquisitions): `reacquire()` runs
            # on the DECODE thread (m11_watch bad_run path) and rebinds
            # self._ramp = None -- the cached table dies with the lock's
            # CFO -- while THIS thread sits between the rebuild guard and
            # derot() reading self._ramp on the executor.  On a weak
            # channel the bad-run trigger fires every few seconds, and
            # eventually one lands inside that window.  Snapshot the ramp
            # into a LOCAL: the block in flight keeps the old lock's ramp
            # (its samples belong to the old lock anyway; the generation
            # guards downstream discard them), and a rebind is atomic.
            ramp = self._ramp
            if ramp is None or len(ramp) < need:
                ramp = np.exp(self.k * np.arange(need) / FS_POST)
                self._ramp = ramp

            def derot(lo, hi):
                seg = y[lo:hi]
                seg *= ramp[:hi - lo]
                seg *= np.exp(self.k * (lo + base) / FS_POST)
        else:
            def derot(lo, hi):
                y[lo:hi] *= np.exp(self.k * np.arange(lo + base, hi + base)
                                   / FS_POST)
        t = pc()
        if self.ex is None:
            derot(0, len(y))
        else:
            F9.par_apply(derot, len(y), self.ex)
        self.tm["derotate"] += pc() - t

        t = pc()
        if self.notch and not self.fft_notch:
            # Overlap-save is exact once the history is REAL: an FIR run on a
            # block preceded by len(b)-1 genuine samples sees the same summands
            # in the same order.  The one place that argument does not hold is
            # the very first block, where the history is 256 fabricated zeros.
            # `lfilter` on the whole array starts from an IMPLICIT zero state
            # and so sums FEWER terms for the first 256 outputs; padding makes
            # it sum all 257, and the two groupings round differently.  That --
            # and only that -- was the 85-sample, 1.11e-16 divergence the M11
            # front-end gate reported: every mismatch fell inside the filter's
            # startup transient.  So prime the first block the way the
            # reference does instead of padding it.
            if not self.notch_primed:
                # the history is the filter's INPUT tail, never its output
                tail = np.concatenate((self.hist, y))[-(NOTCH_TAPS - 1):]
                y = (F9.par_lfilter(self.nb, y, self.ex)
                     if self.ex is not None else lfilter(self.nb, 1.0, y))
                self.hist = tail
                self.notch_primed = True
            else:
                buf = np.concatenate((self.hist, y))
                y = (F9.par_lfilter(self.nb, buf, self.ex)
                     if self.ex is not None
                     else lfilter(self.nb, 1.0, buf))[NOTCH_TAPS - 1:]
                self.hist = buf[-(NOTCH_TAPS - 1):]
        self.tm["notch"] += pc() - t
        if self.trim:
            k = min(self.trim, len(y))
            self.trim -= k
            y = y[k:]
        self.tape.append(y)

    # A tape this long without producing a Frame is not a slow decoder, it is
    # a lost cursor: 4 s is over 16 Frames of slack, far more than any
    # scheduling hiccup needs. The tape holds complex128, so 4 s is 442 MB --
    # against the 34.9 GB actually observed. (First written as 8 s on the
    # arithmetic of an 8-byte sample; the gate printed 885 MB and corrected
    # it. The tape is 16 bytes per sample, not 8.)
    TAPE_MAX = int(4.0 * FS_POST)   # E82: >= 3 Frames on every geometry
                                    # here (RF25's Frame is 242 ms)

    def frames(self):
        """-> list of (index, window, local_t0, coh).  The window is a copy, so
        the tape behind it can be released immediately."""
        out = []
        # SNAPSHOT THE GENERATION, exactly as `_acquire` does for `self.acq`.
        #
        # `reacquire()` runs on the DECODE thread and replaces `self.tape` and
        # `self.centre` while THIS runs on the front thread. The acq race was
        # found and fixed; this one, on the same two threads and the same
        # reset, was not. If the reset lands mid-loop, `centre` keeps a value
        # measured against the OLD tape while the new tape starts empty, the
        # `tape.end < lo + FRAME_WINDOW` test below is then false forever,
        # nothing is ever trimmed, and `_stream` keeps appending complex128 at
        # ~110 MB/s. Measured 8/07: RSS 1.2 GB -> 34.9 GB in three minutes with
        # the Frame counter frozen at 457 -- which the supervisor correctly
        # called a wedge without anyone knowing why.
        gen = self.gen
        if self.plan_cb is not None and self.plan is None:
            if not self._probe_plan():
                return out
        # E55: the FFT window notch needs NOTCH_TAPS//2 = 128 samples of
        # raw margin on EACH side of the window (zero-phase alignment +
        # real convolution history).  Read the armed/decided flags once --
        # reacquire() on the decode thread may flip them mid-loop, and a
        # half-updated pair must not change this loop's arithmetic.
        half = NOTCH_TAPS // 2
        fftn = self.fast and self.fft_notch and self.notch
        while self.ready:
            if self.gen != gen:
                self.stats["frames_raced"] += 1
                return out
            lo = self.centre - FT_SPAN
            if lo < self.tape.org:
                # cannot happen unless trimming got ahead of the cursor
                raise RuntimeError("front end lost the Frame cursor")
            if self.tape.end < lo + self.frame_window + (half if fftn else 0):
                # BOUND THE TAPE regardless of WHY the cursor is unreachable.
                # The generation guard above fixes the race we identified; this
                # makes any other cause of a lost cursor cost a re-acquisition
                # instead of the machine. An unbounded buffer fed by a radio is
                # never acceptable, whatever the reason it stopped draining.
                if self.tape.end - self.tape.org > self.TAPE_MAX:
                    self.stats["tape_overrun"] += 1
                    held = (self.tape.end - self.tape.org) / FS_POST
                    log(f"  *** front-end tape {held:.1f} s with no Frame "
                        f"(cursor at {self.centre}, tape ends "
                        f"{self.tape.end}) -- re-acquiring rather than "
                        f"growing ***")
                    self.reacquire()
                break
            pc = time.perf_counter
            if fftn:
                # zero-phase window filtering: output t is
                # sum_j nb[j] * raw[t + half - j] -- the exact sample the
                # stream-FIR path's tape (filter + trim half) holds at t.
                # Take [lo-half, lo+FRAME_WINDOW+half): 2*half = 256 extra
                # raw samples, the first 256 of which are the overlap-save
                # history.  At the very start of a lock the margin can
                # reach behind the tape origin; zeros there reproduce the
                # reference's own unprimed-lfilter startup transient.
                t = pc()
                a = lo - half
                zf = max(0, self.tape.org - a)
                raw = self.tape.take(a + zf,
                                     self.frame_window + 2 * half - zf)
                if zf:
                    raw = np.concatenate(
                        (np.zeros(zf, raw.dtype), raw))
                self.tm["window"] += pc() - t
                t = pc()
                w = self.osn.apply(raw, self.ex)
                self.tm["notch"] += pc() - t
            else:
                t = pc()
                w = self.tape.take(lo, self.frame_window)
                self.tm["window"] += pc() - t
            t = pc()
            t0, coh = self.ft.search(w, FT_SPAN, self.ex, abs_off=lo)
            self.tm["fine_timing"] += pc() - t
            out.append((self.frame, w, t0, coh))
            if len(out) > 256:
                # a cursor minting >63 s of Frames from one push is broken,
                # whatever broke it -- re-acquire, never spin (E53's own
                # coordinate bug produced exactly this; bound the loop).
                self.stats["frame_storm"] += 1
                log("  *** front end minted 256 Frames in one push -- "
                    "re-acquiring ***")
                self.reacquire()
                break
            self.frame += 1
            # In fast mode the tracker's drift estimate joins the cursor
            # prediction (t0 is ABSOLUTE either way; only the guess of where
            # to look next improves).  drift() is 0 on the exact path.
            self.centre = (lo + t0) + self.frame_samples + self.ft.drift()
            # fast mode retains 2*half extra samples UNCONDITIONALLY (not
            # just when fftn): under-retaining after a notch decision flip
            # would make the margin take() raise; 256 samples is 590 us.
            self.tape.trim(self.centre - FT_SPAN - 8
                           - (2 * half if self.fast else 0))
        return out


def gate_front_end(path, rate, fmt=None, frames=6, block=250_000,
                   threads=16, verbose=True):
    """The streaming front end against `m9_fast.load_fast`, bit for bit.

    Same air, same bootstrap, same CFO, same notch -- and the streaming path
    runs it in 36 ms blocks while the reference runs it in one array.  If the
    resampler phase or the FIR history is wrong by one sample this fails, and
    it fails on `array_equal`, not on a tolerance.
    """
    from concurrent.futures import ThreadPoolExecutor
    import m9_fast
    fmt = fmt or MP.guess_format(path)
    span = 0.05 + frames * FRAME_SEC
    ex = ThreadPoolExecutor(max_workers=threads)
    y_ref = m9_fast.load_fast(path, rate, fmt=fmt,
                              span_sec=span * (rate / FS_POST), ex=ex)[0]
    # gate the THREADED front end, because that is the one that runs live
    fe = FrontEnd(rate, ex=ex)
    mult = fe.rs.block_multiple()
    block = (block // mult) * mult
    pos, need = 0, int(span * rate * 1.05)
    while pos < need and fe.tape.end < len(y_ref):
        x = MP.read_block(path, pos, block, fmt)
        if not len(x):
            break
        pos += len(x)
        fe.push(x)
    n = min(fe.tape.end, len(y_ref))
    got = fe.tape.take(fe.tape.org, n - fe.tape.org)
    ok = n > 0 and np.array_equal(y_ref[fe.tape.org:n], got)
    ex.shutdown()
    if verbose:
        print(f"    {'PASS' if ok else 'FAIL'}  FrontEnd vs load_fast: "
              f"{n} complex samples bit-identical "
              f"({block}-sample blocks, one bootstrap search, "
              f"{'resampled' if not fe.rs.bypass else 'no resampler'})")
        if not ok and n:
            d = np.abs(y_ref[fe.tape.org:n] - got)
            print(f"          first mismatch at {int(np.argmax(d > 0))}, "
                  f"max |delta| {d.max():.3e}")
    return ok


# ---------------------------------------------------------------------------
# sources
# ---------------------------------------------------------------------------

class FileSource:
    """A banked capture, optionally throttled to real time so a replay is a
    fair rehearsal for live air rather than a benchmark in disguise."""

    def __init__(self, path, rate, fmt=None, block=250_000, realtime=False,
                 start=0.0):
        self.path = path
        self.rate = rate
        self.fmt = fmt or MP.guess_format(path)
        self.block = block
        self.realtime = realtime
        self.pos = int(start * rate)
        self.t0 = None
        self.n = 0
        self.stats = collections.Counter()

    def start(self):
        self.t0 = time.time()

    def read(self):
        x = MP.read_block(self.path, self.pos, self.block, self.fmt)
        if not len(x):
            return None
        self.pos += len(x)
        self.n += len(x)
        if self.realtime:
            d = self.t0 + self.n / self.rate - time.time()
            if d > 0:
                time.sleep(d)
        return x

    def close(self):
        pass


class SdrSource:
    """The radio, read continuously.

    Fleet discipline, and every clause of it is load-bearing:
      * radio_lock owner `atsc3_watch`, priority 60 -- a human is watching, so
        above lab background and below a satellite pass or a window hold;
      * polite wait and bounded retries, never a seizure;
      * **heartbeat on a TIMER inside the read loop, never per read.**  The
        8/03 law was paid for in a decoder that reported ~97 wpm because
        per-read file I/O was chopping the stream out from under it;
      * `should_yield()` checked on the same timer, so a higher-priority
        holder gets the radio without anyone having to kill anything;
      * release in `finally`;
      * refuses to run inside the 02:10-05:40 meteor window.
    """

    def __init__(self, rf, rate=FS_POST, ant="Antenna B", ifgr=None,
                 rfgain=None, owner="atsc3_watch", pri=60, wait=120.0,
                 block=250_000, radio_lock=None):
        self.rf = rf
        self.want_rate = rate
        self.ant = ant
        self.ifgr = ifgr          # None = consult the gain table
        self.rfgain = rfgain      # None = consult the gain table
        self.owner = owner
        self.pri = pri
        self.wait = wait
        self.block = block
        self.rl = radio_lock
        self.sdr = None
        self.st = None
        self.held = False
        self.rate = None
        self.readback = {}
        self.gain_source = None
        self.stats = collections.Counter()
        self._hb = 0.0
        self.yield_reason = None
        scale = DTYPES["cs16"][1]
        self.scale = np.float32(scale)
        self._buf = np.empty(2 * 65536, np.int16)
        self._acc = []
        self._acc_n = 0

    @staticmethod
    def center_hz(rf):
        lo = (174 + (rf - 7) * 6) if rf < 14 else (470 + (rf - 14) * 6)
        return (lo + 3.0) * 1e6

    def open(self):
        import SoapySDR
        from SoapySDR import SOAPY_SDR_CS16, SOAPY_SDR_RX
        SoapySDR.SoapySDR_setLogLevel(SoapySDR.SOAPY_SDR_FATAL)
        if self.rl is not None:
            if not self.rl.acquire(self.owner, f"ATSC 3.0 live watch RF{self.rf}",
                                   self.pri, wait_s=self.wait):
                h = self.rl.status() or {}
                raise RuntimeError(
                    f"radio held by {h.get('owner')} "
                    f"(priority {h.get('priority')}); not seizing")
            self.held = True
        self.sdr = SoapySDR.Device("driver=sdrplay")
        self.sdr.setSampleRate(SOAPY_SDR_RX, 0, self.want_rate)
        time.sleep(0.2)
        # a write without a readback is a hope, not a setting
        self.rate = float(self.sdr.getSampleRate(SOAPY_SDR_RX, 0))
        self.sdr.setAntenna(SOAPY_SDR_RX, 0, self.ant)
        try:
            self.sdr.setGainMode(SOAPY_SDR_RX, 0, False)
        except Exception:                                      # noqa: BLE001
            pass
        # E79: the right gain is a property OF THE CARRIER, looked up, not
        # a constant carried in an argument default. RF33 needs rfgain_sel 4
        # (at 2 it reads 0.0% FEC while an HDHomeRun on the same feed
        # decodes at snq=100) and RF25 wants 2 (+6.05 dB); the tuner is
        # single-tenant so one global value is WRONG for one of them
        # whichever it is. An explicit --rfgain still wins, and either way
        # the write is confirmed by readback before we believe it.
        import carrier_gain as CG
        self.sdr.setFrequency(SOAPY_SDR_RX, 0, self.center_hz(self.rf))
        try:
            g = CG.apply_gain(self.sdr, self.rf, log=log,
                              rfgain=self.rfgain, ifgr=self.ifgr)
            rg, self.ifgr = g["rfgain_sel"], g["ifgr"]
            self.gain_source = g["source"]
        except CG.GainReadbackError as e:
            # NOT a warning to scroll past: the front end is not where we
            # think it is, and on this carrier that can mean 0.0% FEC with
            # every other indicator looking healthy.
            log(f"  *** GAIN NOT CONFIRMED: {e} ***")
            raise
        time.sleep(0.3)
        self.readback = dict(
            rate=self.rate, rate_exact=abs(self.rate - self.want_rate) < 1.0,
            freq=float(self.sdr.getFrequency(SOAPY_SDR_RX, 0)),
            antenna=self.sdr.getAntenna(SOAPY_SDR_RX, 0), rfgain_sel=rg)
        try:
            self.readback["bandwidth"] = float(
                self.sdr.getBandwidth(SOAPY_SDR_RX, 0))
        except Exception:                                      # noqa: BLE001
            pass
        self.st = self.sdr.setupStream(SOAPY_SDR_RX, SOAPY_SDR_CS16)
        self.sdr.activateStream(self.st)
        # a sacrificial first read is the 8/01 meteor lesson; half a second
        # rather than a quarter because a startup overflow costs a bootstrap
        # re-acquisition, which is 0.3 s of air thrown away
        t = time.time()
        while time.time() - t < 0.5:
            self.sdr.readStream(self.st, [self._buf], 65536, timeoutUs=500000)
        self._hb = time.time()
        return self.readback

    def start(self):
        pass

    def read(self):
        """One block of `self.block` complex64 samples, or None on a break.

        Nothing here touches the filesystem except through the 5-second
        heartbeat timer.
        """
        while self._acc_n < self.block:
            r = self.sdr.readStream(self.st, [self._buf], 65536,
                                    timeoutUs=500000)
            now = time.time()
            if now - self._hb >= 5.0:
                self._hb = now
                if self.rl is not None:
                    self.rl.heartbeat()
                    y = self.rl.should_yield()
                    if y:
                        self.yield_reason = y
                        return None
            if r.ret > 0:
                iq = self._buf[:2 * r.ret].astype(np.float32) * self.scale
                self._acc.append(iq[0::2] + 1j * iq[1::2])
                self._acc_n += r.ret
                self.stats["samples"] += r.ret
            elif r.ret < 0:
                self.stats["overflow_reads"] += 1
                self.stats["breaks"] += 1
                self._acc = []
                self._acc_n = 0
                return b"BREAK"
        a = np.concatenate(self._acc).astype(np.complex64)
        out, rest = a[:self.block], a[self.block:]
        self._acc = [rest] if len(rest) else []
        self._acc_n = len(rest)
        return out

    def close(self):
        try:
            if self.st is not None:
                self.sdr.deactivateStream(self.st)
                self.sdr.closeStream(self.st)
        except Exception:                                      # noqa: BLE001
            pass
        if self.held and self.rl is not None:
            self.rl.release(self.owner)
            self.held = False


# ---------------------------------------------------------------------------
# the player
# ---------------------------------------------------------------------------

class Player:
    """A player on stdin, and an HONEST occupancy meter in front of it.

    A pipe cannot tell us the player is starving, so we do not ask it.  Each
    segment is a known number of seconds of media; once playback has started,
    the player has consumed one second of media per second of wall clock.  So

        occupancy = seconds pushed - seconds of wall since playback started

    is what the player still holds, and an occupancy that reaches zero is an
    underrun whether or not anything on screen looked wrong.  Prebuffering
    raises the starting occupancy; it does not change the measurement.
    """

    def __init__(self, cmd=None, prebuffer=5.0, out_path=None, timescale=90000,
                 quiet=False):
        self.cmd = cmd
        self.prebuffer = prebuffer
        self.out_path = out_path
        self.timescale = timescale
        self.quiet = quiet
        self.proc = None
        self.fh = None
        self.q = queue.Queue(maxsize=256)
        self.thread = None
        self.stop = threading.Event()
        self.started = False
        self.t_play = None
        self.media_s = 0.0
        self.pushed_s = 0.0
        self.hold = []
        self.hold_s = 0.0
        self.stats = collections.Counter()
        self.min_occ = None
        self.bytes = 0
        self.dead = False

    def open(self, init):
        import subprocess
        if self.out_path:
            self.fh = open(self.out_path, "wb")
            self.fh.write(init)
        if self.cmd:
            self.proc = subprocess.Popen(
                self.cmd, stdin=subprocess.PIPE,
                stdout=subprocess.DEVNULL if self.quiet else None,
                stderr=subprocess.DEVNULL if self.quiet else None)
            self.proc.stdin.write(init)
            self.proc.stdin.flush()
        self.bytes += len(init)
        self.thread = threading.Thread(target=self._pump, daemon=True)
        self.thread.start()

    def _pump(self):
        while not self.stop.is_set():
            try:
                d = self.q.get(timeout=0.5)
            except queue.Empty:
                continue
            if d is None:
                break
            try:
                if self.proc is not None and self.proc.poll() is None:
                    self.proc.stdin.write(d)
                    self.proc.stdin.flush()
                elif self.proc is not None:
                    self.dead = True
                if self.fh is not None:
                    self.fh.write(d)
            except Exception:                                  # noqa: BLE001
                self.dead = True
            self.bytes += len(d)

    def occupancy(self, now=None):
        if not self.started:
            return self.hold_s
        return self.pushed_s - ((now or time.time()) - self.t_play)

    def push(self, seg):
        dur = seg["dur"] / self.timescale
        self.media_s += dur
        if not self.started:
            self.hold.append(seg)
            self.hold_s += dur
            if self.hold_s < self.prebuffer:
                return
            self.started = True
            self.t_play = time.time()
            for s in self.hold:
                self.q.put(s["bytes"])
                self.pushed_s += s["dur"] / self.timescale
            self.hold = []
            self.min_occ = self.pushed_s
            return
        occ = self.occupancy()
        if self.min_occ is None or occ < self.min_occ:
            self.min_occ = occ
        if occ < 0:
            # The player ran dry for exactly -occ seconds and then resumed, so
            # the playback clock is that much behind the wall clock from here
            # on.  Shifting `t_play` by the deficit records the stall instead
            # of either forgetting it (every later segment then looks early)
            # or carrying it forward (every later segment then looks late).
            self.stats["underrun"] += 1
            self.stats["underrun_s"] += -occ
            self.t_play += -occ
        self.q.put(seg["bytes"])
        self.pushed_s += dur
        self.stats["segment"] += 1

    def close(self):
        self.stop.set()
        try:
            self.q.put_nowait(None)
        except Exception:                                      # noqa: BLE001
            pass
        if self.thread is not None:
            self.thread.join(timeout=3)
        if self.fh is not None:
            self.fh.close()
        if self.proc is not None:
            try:
                self.proc.stdin.close()
            except Exception:                                  # noqa: BLE001
                pass
            try:
                self.proc.wait(timeout=3)
            except Exception:                                  # noqa: BLE001
                self.proc.kill()
