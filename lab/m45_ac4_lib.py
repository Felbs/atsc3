"""M45 -- the receiver's AC-4 through libac4dec (the ac4-decoder C library).

WHY. The receiver's own decoder (m17-m43, driven by m42.Ac4Stream) is the
code `ac4-decoder` was copied from at c4b296e. Since the copy, ac4-decoder
fixed, in its Python AND its C port, each with a test that fails without it:

  * full SAP (sap_mode 3) -- 22-34 % of 5.1 frames on air were decoded as if
    the pair were L/R (the receiver passed them through);
  * A-SPX previous_stop_pos under aspx_balance (a parse that went off the
    rails on balanced frames -- the "A-SPX bursts" core_for_bursts papers
    over);
  * the A-CPL upmix of ASPX_ACPL_3 (RF25's 5.1: the receiver writes C, Ls,
    Rs, LFE silent);
  * the IMDCT scale of short blocks, A-SPX envelope dequantisation (1.5 dB
    steps), the QMF time grid (TS_OFFSET), A-SPX on the centre, companding,
    the noise index, limiter groups, transient envelopes, spectral noise
    fill;
  * every channel through the same 577-sample QMF delay (the receiver's
    centre and LFE lead its fronts by 12 ms);
  * damaged-frame concealment and the true-peak output guard;
  * loudness / dialnorm, DRC, downmixes with the stream's coefficients.

Back-porting that piecemeal would re-write ac4-decoder inside the receiver
-- two copies of a decoder diverge (AC4_SPLIT_PLAN.md says so). Instead the
worker calls the library, and keeps m42 as the fallback when the library is
missing, refuses to open, or is switched off (ATSC3_AC4_BACKEND=python).

THE LIBRARY. Found at $ATSC3_AC4DEC_LIB, else next to this repository at
../ac4-decoder/c/libac4dec.so (build it there with `make -C c libac4dec.so`).
Nothing is vendored. The C decoder streams (no decode passes), at 18-25x real
time for 5.1 on one core.

WHAT THIS CLASS GIVES THE WORKER. The same shape of output the worker built
from m42 + m30 + m37: per channel, float64 at the INTERNAL rate (1536
samples a frame at 29.97 fps), decoder scale (the worker's REF_PEAK rule
applies unchanged), whole frames only. Two delays are taken out so the
worker's frame bookkeeping (padding, resume, the sidecar's first_seq) stays
true:

  * the lookahead (ac4dec_delay, >= 8 frames): a frame comes out that many
    calls after it went in. The class counts what it has emitted, so the
    worker knows which LANE frame its output starts at (`out_frame`).
  * the 577-sample QMF delay of every channel: the first 577 samples after
    an open are dropped, so output sample 0 is the first frame's sample 0
    (the old path left its A-SPX pairs 577 late against C / LFE / video).

Hence the output lags the input by lookahead + 1 frames (the remainder of
the 577 is held for the next call) -- 9 frames = 0.3 s at the minimum
lookahead. `flush()` drains everything at the end of a recording.
"""
from __future__ import annotations

import ctypes
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

FRAME = 1536            # internal samples per frame at 29.97 fps
QMF_DELAY = 577         # every channel's delay through the QMF bank
MIN_LOOKAHEAD = 8       # AC4DEC_MIN_LOOKAHEAD

# ac4dec_element / ac4dec_layout (c/include/ac4dec.h)
EL_AUTO, EL_5X, EL_PAIR = 0, 1, 2
LAYOUT_NATIVE, LAYOUT_2_0 = 0, 1
STATUS_OK = 0
MAX_CH = 16


class _Config(ctypes.Structure):
    _fields_ = [("element", ctypes.c_int), ("channels", ctypes.c_int),
                ("internal_rate", ctypes.c_int), ("lookahead", ctypes.c_int),
                ("no_guard", ctypes.c_int), ("presentation", ctypes.c_int),
                ("presentation_id", ctypes.c_long),
                ("language", ctypes.c_char * 16),
                ("loudness", ctypes.c_int), ("target_level", ctypes.c_double),
                ("drc", ctypes.c_int), ("ch_mode", ctypes.c_int),
                ("de_level", ctypes.c_double), ("substream", ctypes.c_int),
                ("no_mix", ctypes.c_int), ("assoc_gain", ctypes.c_double),
                ("layout", ctypes.c_int), ("present", ctypes.c_int),
                ("ajoc_core", ctypes.c_int), ("stereo_dmx", ctypes.c_int)]


class _Frame(ctypes.Structure):
    _fields_ = [("status", ctypes.c_int), ("reason", ctypes.c_char_p),
                ("index", ctypes.c_longlong), ("n_channels", ctypes.c_int),
                ("n_samples", ctypes.c_int),
                ("pcm", ctypes.POINTER(ctypes.c_double) * MAX_CH),
                ("sequence_counter", ctypes.c_int), ("iframe", ctypes.c_int),
                ("seq_jump", ctypes.c_int), ("sample_rate", ctypes.c_int)]


class _Stats(ctypes.Structure):
    _fields_ = [("n_frames", ctypes.c_long), ("n_out", ctypes.c_long),
                ("n_ok", ctypes.c_long), ("n_concealed", ctypes.c_long),
                ("n_bad", ctypes.c_long), ("n_guard", ctypes.c_long),
                ("n_aspx_level", ctypes.c_long), ("n_late", ctypes.c_long),
                ("n_burst", ctypes.c_long), ("n_seq_jumps", ctypes.c_long),
                ("element", ctypes.c_int), ("acpl3", ctypes.c_int),
                ("frame_len", ctypes.c_int), ("out_rate", ctypes.c_int),
                ("aspx", ctypes.c_int)]


# The structs are handed over inside zeroed buffers this much larger than
# ours: a newer library whose config grew reads zeros (= its defaults) past
# our fields instead of whatever follows in memory, and one whose frame or
# stats grew writes into slack instead of over the heap.
_SLACK = 1024


def _boxed(struct_cls):
    buf = ctypes.create_string_buffer(ctypes.sizeof(struct_cls) + _SLACK)
    return buf, struct_cls.from_buffer(buf)


def lib_candidates():
    env = os.environ.get("ATSC3_AC4DEC_LIB", "").strip()
    out = [env] if env else []
    out.append(os.path.join(os.path.dirname(ROOT), "ac4-decoder", "c",
                            "libac4dec.so"))
    return out


_LIB = None


def load_lib():
    """-> the ctypes library, or raise OSError saying why not."""
    global _LIB
    if _LIB is not None:
        return _LIB
    why = []
    for path in lib_candidates():
        if not os.path.exists(path):
            why.append(f"{path}: not found")
            continue
        try:
            lib = ctypes.CDLL(path)
            lib.ac4dec_open.restype = ctypes.c_void_p
            lib.ac4dec_open.argtypes = [ctypes.c_void_p]
            lib.ac4dec_close.restype = None
            lib.ac4dec_close.argtypes = [ctypes.c_void_p]
            lib.ac4dec_delay.restype = ctypes.c_int
            lib.ac4dec_delay.argtypes = [ctypes.c_void_p]
            lib.ac4dec_channels.restype = ctypes.c_int
            lib.ac4dec_channels.argtypes = [ctypes.c_void_p]
            lib.ac4dec_decode_frame.restype = ctypes.c_int
            lib.ac4dec_decode_frame.argtypes = [
                ctypes.c_void_p, ctypes.c_char_p, ctypes.c_size_t,
                ctypes.c_void_p]
            lib.ac4dec_flush.restype = ctypes.c_int
            lib.ac4dec_flush.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
            lib.ac4dec_get_stats.restype = None
            lib.ac4dec_get_stats.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
            lib.ac4dec_reason.restype = ctypes.c_char_p
            lib.ac4dec_reason.argtypes = [
                ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ctypes.c_long)]
            lib.ac4dec_ref_peak.restype = ctypes.c_double
            lib.ac4dec_ref_peak.argtypes = [ctypes.c_int]
            # this binding's contract: the 5_X / pair reference scale the
            # worker's REF_PEAK rule is written for
            if (abs(lib.ac4dec_ref_peak(EL_5X) - 9e5) > 1
                    or abs(lib.ac4dec_ref_peak(EL_PAIR) - 2e5) > 1):
                raise OSError("unexpected decoder scale (ac4dec_ref_peak)")
            lib._atsc3_path = path
            _LIB = lib
            return lib
        except (OSError, AttributeError) as e:
            why.append(f"{path}: {e}")
    raise OSError("; ".join(why) or "no candidate path")


class LibAc4Stream:
    """A live AC-4 decoder on libac4dec; see the module docstring.

    element: "5_X" or "pair" (the worker's --element); channels: 6 (5.1,
    L R C LFE Ls Rs) or 2. Two channels of a 5_X element are the standard's
    LoRo downmix with the stream's own coefficients (centre and surrounds
    folded in, TS 103 190-1 6.2.17.4) -- never L/R alone, which drops the
    dialogue; an ASPX_ACPL_3 stream gives its transmitted stereo core.
    """

    CHANS6 = ("L", "R", "C", "lfe", "Ls", "Rs")

    def __init__(self, element="5_X", channels=6, lookahead=None):
        self.lib = load_lib()
        self.element = element
        self.n_out_ch = 2 if (element == "pair" or channels == 2) else 6
        self.chans = (("L", "R") if self.n_out_ch == 2 else self.CHANS6)
        self.lookahead = int(lookahead or os.environ.get(
            "ATSC3_AC4_LOOKAHEAD", MIN_LOOKAHEAD))
        self._cfg_buf, cfg = _boxed(_Config)
        cfg.element = EL_PAIR if element == "pair" else EL_5X
        if element == "pair":
            cfg.channels = 0
        elif channels == 2:
            cfg.channels = 0
            cfg.layout = LAYOUT_2_0          # LoRo (stereo_dmx 0)
        else:
            cfg.channels = 6
        cfg.internal_rate = 1                # the worker resamples itself
        cfg.lookahead = max(MIN_LOOKAHEAD, self.lookahead)
        cfg.loudness = 0                     # decoder scale: REF_PEAK rule
        self._fr_buf, self._fr = _boxed(_Frame)
        self._st_buf, self._st = _boxed(_Stats)
        self.d = None
        self.cfg = None                      # (m42's A-SPX config: not ours)
        self.base = 0                        # lane frame of the first frame fed
        self.out_frame = 0                   # lane frame of the next output
        self._open()

    # -- lifecycle -------------------------------------------------------
    def _open(self):
        d = self.lib.ac4dec_open(ctypes.addressof(self._cfg_buf))
        if not d:
            raise OSError("ac4dec_open refused the configuration")
        self.d = d
        self._skip = QMF_DELAY               # samples still to drop
        self._fifo = [np.zeros(0) for _ in range(self.n_out_ch)]
        self._chunks = [[] for _ in range(self.n_out_ch)]
        self._n_in = 0
        self._n_stream_out = 0               # frames emitted since open
        self._ended = False

    def fresh(self):
        """True while nothing has been fed since the (re)open."""
        return self._n_in == 0

    def align(self, lane_frame):
        """Before the first frame of a stream: say which lane frame it is."""
        if self._n_in == 0:
            self.base = self.out_frame = int(lane_frame)

    def reset(self, base_frame=0):
        """A new stream (lane roll, resume): close, reopen, and say which
        lane frame the next frame fed is."""
        self._keep_stats()
        self.close()
        self._open()
        self.base = self.out_frame = int(base_frame)

    def close(self):
        if self.d:
            self.lib.ac4dec_close(self.d)
            self.d = None

    def __del__(self):
        try:
            self.close()
        except Exception:                                      # noqa: BLE001
            pass

    # -- counters (summed over reopenings) -------------------------------
    _prev = None

    def _keep_stats(self):
        s = self.stats()
        self._prev = s

    def stats(self):
        out = dict(n_frames=0, n_out=0, n_ok=0, n_concealed=0, n_bad=0,
                   n_guard=0, n_aspx_level=0, n_late=0, n_burst=0,
                   n_seq_jumps=0, acpl3=0)
        if self.d:
            self.lib.ac4dec_get_stats(self.d, ctypes.addressof(self._st_buf))
            for k in out:
                out[k] = int(getattr(self._st, k))
        if self._prev:
            for k in out:
                out[k] = (max(out[k], self._prev[k]) if k == "acpl3"
                          else out[k] + self._prev[k])
        return out

    def reasons(self):
        out = {}
        if not self.d:
            return out
        n = ctypes.c_long(0)
        i = 0
        while True:
            r = self.lib.ac4dec_reason(self.d, i, ctypes.byref(n))
            if r is None:
                break
            out[r.decode("ascii", "replace")] = int(n.value)
            i += 1
        return out

    @property
    def n_bad(self):
        return self.stats()["n_concealed"]

    @property
    def acpl3(self):
        return bool(self.stats()["acpl3"])

    # -- decoding --------------------------------------------------------
    def _take(self):
        fr = self._fr
        n = int(fr.n_samples)
        nc = int(fr.n_channels)
        for c in range(self.n_out_ch):
            x = np.zeros(FRAME)
            m = min(n, FRAME)               # the worker's ruler is 1536
            if c < nc and m > 0:
                x[:m] = np.ctypeslib.as_array(fr.pcm[c], shape=(n,))[:m]
            self._chunks[c].append(x)

    def _emit(self, final=False):
        if self._chunks[0]:
            self._fifo = [np.concatenate([f] + ch)
                          for f, ch in zip(self._fifo, self._chunks)]
            self._chunks = [[] for _ in range(self.n_out_ch)]
        if self._skip:
            k = min(self._skip, len(self._fifo[0]))
            self._fifo = [f[k:] for f in self._fifo]
            self._skip -= k
        avail = len(self._fifo[0])
        if final and avail % FRAME:
            pad = FRAME - avail % FRAME
            self._fifo = [np.concatenate([f, np.zeros(pad)])
                          for f in self._fifo]
            avail += pad
        # never emit more frames than were fed (the tail padding above
        # stands in for the 577 samples the QMF never gives back)
        nf = min(avail // FRAME, self._n_in - self._n_stream_out)
        ns = nf * FRAME
        out = {k: self._fifo[c][:ns] for c, k in enumerate(self.chans)}
        self._fifo = [f[ns:] for f in self._fifo]
        self._n_stream_out += nf
        self.out_frame += nf
        return out

    def push(self, frames, final=False):
        """Feed AC-4 frames (bytes). -> {channel: float64 array} of whole
        frames, starting at lane frame `out_frame` as it was BEFORE the
        call. final=True drains the decoder (end of a recording); the next
        push then starts a new stream."""
        if self._ended and frames:
            self.reset(self.out_frame)
        for f in frames:
            f = bytes(f)
            got = self.lib.ac4dec_decode_frame(self.d, f, len(f),
                                               ctypes.addressof(self._fr_buf))
            self._n_in += 1
            if got:
                self._take()
        if final:
            while self.lib.ac4dec_flush(self.d,
                                        ctypes.addressof(self._fr_buf)):
                self._take()
        out = self._emit(final=final)
        if final:
            self._ended = True
        return out
