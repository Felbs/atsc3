#!/usr/bin/env python3
"""atsc3_audio.py -- decode the live AC-4 lane, behind the write head.

The third decoupled stage. The chain writes live_audio_pidN.m4s and never
waits; this tails it, decodes with our AC-4 decoder, and appends 48 kHz PCM
to live_audio.wav. It runs at ~1.28x real time -- enough to keep up, not
enough to be in anyone's critical path, which is exactly why it is a
separate process reading a file rather than a stage inside the chain.

DECODING A GROWING STREAM IN CHUNKS
-----------------------------------
`m37.decode_all` is batch: it collects every MDCT window and synthesises
once. Live needs to emit as it goes, and you cannot simply resume mid-frame
because two things carry state across frames:

  * the MDCT overlap -- each output block needs the previous window;
  * A-SPX -- envelopes are delta-coded against the previous frame.

Both reset at an I-FRAME (`b_iframe_global`), which the encoder emits
regularly. So each pass restarts at the last i-frame at or before the
cursor, decodes forward, and DISCARDS the samples before the cursor. That
is correct by construction rather than by convergence, and cheap: the
re-decoded lead is only back to the most recent i-frame.

Usage:
    python tools/atsc3_audio.py                     # follow the default dir
    python tools/atsc3_audio.py --pid 13 --once     # one pass, then exit
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import wave

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "lab"))
sys.path.insert(0, ROOT)

OUT_FS = 48000
FRAME_SECONDS = 1001.0 / 30000.0
from concurrent.futures import ThreadPoolExecutor              # noqa: E402
_EX = ThreadPoolExecutor(max_workers=6)


def lanes(live_dir):
    try:
        d = json.load(open(os.path.join(live_dir, "live.json")))
        return d.get("lanes", {})
    except (OSError, ValueError):
        return {}


def frag_seq(idx_path, frag_k):
    """MPU sequence number of fragment number `frag_k` (idx position).

    E51: the old form of this helper took a FRAME index and divided by 60
    -- which mislocates the anchor whenever any earlier fragment is
    truncated (a 2026-08-08 SDR-overflow afternoon put 286 sub-60
    fragments in one lane; //60 was off by ~100 fragments = minutes).
    Fragment position is the only unit the idx actually speaks."""
    try:
        seqs = [json.loads(l)["seq"] for l in open(idx_path)]
    except (OSError, ValueError):
        return None
    return seqs[frag_k] if 0 <= frag_k < len(seqs) else \
        (seqs[0] if seqs else None)


def pad_insertions(fcounts, anchor_frag, cursor, n_new, per_mpu=60,
                   y_first=None):
    """Where must silence be inserted INSIDE this batch? -> [(pos, miss)].

    E59: batch-END padding (E51's recorded open item) time-compresses
    everything after a hole WITHIN the batch -- content pulled up to the
    batch's whole pad total early, then a silence block, then a snap back
    to alignment. On a dipping band that is an audible artifact per pad
    event (the "audio cuts out ... then resumes" report) and a local
    displacement the slot->offset arithmetic cannot see. The cure is to
    put each truncated fragment's missing frames AT that fragment's
    position: pos is in FRAMES from the start of the emitted batch
    (y_first, default cursor), one entry per truncated fragment whose
    last frame lands in this batch, pre-anchor fragments excluded."""
    if y_first is None:
        y_first = cursor
    csum = [0]
    for c in fcounts:
        csum.append(csum[-1] + c)
    out = []
    import bisect
    k = max(0, bisect.bisect_right(csum, cursor) - 1)
    end_batch = cursor + n_new
    while k < len(fcounts) and csum[k] < end_batch:
        fr_end = csum[k + 1]
        miss = per_mpu - fcounts[k]
        if miss > 0 and k >= anchor_frag and cursor < fr_end <= end_batch \
                and fr_end - y_first >= 0:
            out.append((fr_end - y_first, miss))
        k += 1
    return out


def fresh_anchor(fcounts, skip0, per_mpu=60):
    """FRESH start-behind anchoring -> (trim_frames, anchor_frag, pad_base).

    E51: the shipped arithmetic trimmed skip0*60 FRAMES for skip0 skipped
    FRAGMENTS and baselined pads to the whole-lane shortfall. Both are only
    correct when every fragment declares exactly 60 frames. On the 8/08
    SDR-overflow lane (286 truncated fragments) the trim overshot by 6497
    frames = 217 s -- the wav's sample 0 landed minutes past the sidecar's
    first_seq, the viewer read stale sport under a live commercial, and the
    whole-lane baseline silently un-padded the consumed fragments' real
    holes. Frames are trimmed by what the skipped fragments DECLARE; pads
    are baselined to the SKIPPED region's shortfall only."""
    k = min(skip0, len(fcounts))
    trim = sum(fcounts[:k])
    pad_base = sum(per_mpu - c for c in fcounts[:k] if c < per_mpu)
    return trim, k, pad_base


def count_frags(idx_path):
    try:
        return sum(1 for _ in open(idx_path))
    except OSError:
        return 0


def lane_shortfall(path, per_mpu=60):
    """Frames missing from fragments that ARE present, from ONE read.

    Comparing a frame count against a fragment count read separately is a
    race on a growing file: the two reads are at different instants and the
    skew is the same order as the shortfall being measured. Reading them in
    one order gave 158, the other order gave 80, and the true value cannot
    decrease -- so both were wrong.

    No second read is needed. Every fragment DECLARES its own sample count
    in its `trun`, and a truncated one declares fewer (measured: 809
    fragments declaring 60, one declaring 22). The shortfall is therefore
    sum(per_mpu - declared), computed from the same bytes the frames came
    from, and cannot be skewed by anything the writer does afterwards.
    """
    try:
        d = open(path, "rb").read()
    except OSError:
        return 0, 0
    short = nfrag = 0
    o = 0
    while o + 8 <= len(d):
        sz = int.from_bytes(d[o:o + 4], "big")
        typ = d[o + 4:o + 8]
        if sz < 8 or o + sz > len(d):
            break
        if typ == b"moof":
            # first traf's trun is the media one; the second is the hint traf
            p, end, got = o + 8, o + sz, None
            while p + 8 <= end and got is None:
                s2 = int.from_bytes(d[p:p + 4], "big")
                if s2 < 8:
                    break
                if d[p + 4:p + 8] == b"traf":
                    q, e2 = p + 8, p + s2
                    while q + 8 <= e2:
                        s3 = int.from_bytes(d[q:q + 4], "big")
                        if s3 < 8:
                            break
                        if d[q + 4:q + 8] == b"trun":
                            got = int.from_bytes(d[q + 12:q + 16], "big")
                            break
                        q += s3
                p += s2
            if got is not None:
                nfrag += 1
                short += max(0, per_mpu - got)
        o += sz
    return nfrag, short


HF_DELAY = 577      # samples the A-SPX path delays a channel against its core
#                     (measured 11/10: 602 at 48 kHz, correlation 0.997)
HF_TAIL = 2         # frames after a burst frame that also play core-band


def core_for_bursts(hf_a, hf_b, core_a, core_b, lead, limit, owed):
    """One pair's HF-rendered pass, with the core put back wherever A-SPX
    has run away. -> (a, b, frames replaced, frames still owed).

    Measured 11/10 on air: the Ls/Rs pair, digitally silent before and
    after, came out at 5.6e7 for three frames -- 60 times the loudest
    programme sample ever seen. With A-SPX off the same frames peak at 366;
    their envelope indices read 54 where the frames either side read 6, the
    signal gain sits on its 1.58e5 cap. That is a decode error in the
    envelope chain (balance 1 -> 0 at that frame; NOT yet found), and no
    programme. A frame whose peak is above `limit` is therefore replaced by
    its core, delayed as the A-SPX path delays it, together with the
    HF_TAIL frames after it (the filter bank is causal, its memory reaches
    forward only; with one frame the next still peaked at 2.2e4 on a
    channel that idles at 100). `owed` is how many of those frames the
    previous pass left for this one, so the result does not depend on
    where a pass ends. Both channels go together: the pair is rendered
    jointly.
    """
    n = min(len(hf_a), len(hf_b)) // 1536
    if not n:
        return hf_a, hf_b, 0, owed
    pk = np.maximum(np.abs(hf_a[:n * 1536]).reshape(n, 1536).max(axis=1),
                    np.abs(hf_b[:n * 1536]).reshape(n, 1536).max(axis=1))
    hot = pk > limit
    bad = np.zeros(n, dtype=bool)
    left = int(owed)
    for f in range(n):
        if hot[f]:
            bad[f], left = True, HF_TAIL
        elif left:
            bad[f], left = True, left - 1
    if not bad.any():
        return hf_a, hf_b, 0, 0
    out = []
    for hf, core in ((hf_a, core_a), (hf_b, core_b)):
        d = np.concatenate([np.zeros(HF_DELAY), core])[lead:lead + len(hf)]
        hf = np.array(hf, dtype=float)
        for f in np.flatnonzero(bad):
            hf[f * 1536:(f + 1) * 1536] = d[f * 1536:(f + 1) * 1536]
        out.append(hf)
    return out[0], out[1], int(bad.sum()), left


def _size(path):
    try:
        return os.path.getsize(path)
    except OSError:
        return -1


class LaneReader:
    """Incremental reader for an append-only fMP4 lane.

    `m17.samples()` and `lane_shortfall()` each read and parse the ENTIRE
    lane on every call. That is fine for a file and ruinous for a stream:
    the lane grows ~1.4 MB/min, so by 65 minutes each pass was parsing
    ~90 MB twice, and per-pass throughput decayed 1.06x -> 0.79x while the
    decode work per pass stayed identical. The worker fell 12 minutes
    behind and the deficit grew linearly.

    This keeps a byte cursor and parses only what is new. A moof is only
    consumed once its following mdat is COMPLETE, so a fragment caught
    half-written is simply left for the next pass rather than mis-parsed.
    The shortfall accumulates as fragments are seen, so it needs no second
    read (E8) and no rescan (E17).
    """

    def __init__(self, path, per_mpu=60):
        self.path = path
        self.per_mpu = per_mpu
        self.pos = 0
        self.nfrag = 0
        self.short = 0
        self.fcounts = []       # frames DECLARED by each fragment (E51):
        #                         the only honest ruler when fragments are
        #                         truncated -- every frame<->fragment
        #                         conversion must walk this, never //60
        self.rolls = 0
        self.rolled = False
        self.no_progress = 0

    def force_roll(self, why):
        print(f"  lane reader reset ({why})")
        self.pos = 0
        self.nfrag = 0
        self.short = 0
        self.fcounts = []
        self.no_progress = 0
        self.rolls += 1
        self.rolled = True

    def read_new(self):
        """-> list of new AC-4 frames, in order."""
        try:
            size = os.path.getsize(self.path)
            # A SHORTER FILE IS NOT AN EMPTY ONE.
            #
            # `size <= self.pos` folds two completely different situations
            # into one answer: "nothing new yet" and "the lane restarted
            # underneath me". On 8/07 the supervisor restarted a wedged chain,
            # the lane was recreated, and this returned [] on every pass for
            # the rest of the run -- the worker sat at 1277.3 s of audio
            # reporting `behind 0`, because a cursor past EOF reads as caught
            # up. The subtitle worker survived the same event only because it
            # re-reads documents instead of holding a byte cursor.
            #
            # Third time today that clamping a signed quantity destroyed the
            # evidence of backwards motion (psn_lost's dead max(0,...), the
            # E21 monotone latch, this). Read the sign.
            if size < self.pos:
                print(f"  lane ROLLED (size {size} < cursor {self.pos}) -- "
                      f"restarting from the new file")
                self.pos = 0
                self.nfrag = 0
                self.short = 0
                self.fcounts = []
                self.tail = b""
                self.rolls = getattr(self, "rolls", 0) + 1
                self.rolled = True
            elif size == self.pos:
                return []
            with open(self.path, "rb") as f:
                f.seek(self.pos)
                d = f.read()
        except OSError:
            return []
        import m7_play as PL                                   # noqa: E402
        out, o, consumed = [], 0, 0
        while o + 8 <= len(d):
            sz = int.from_bytes(d[o:o + 4], "big")
            if sz < 8 or o + sz > len(d):
                break
            if d[o + 4:o + 8] == b"moof":
                m2 = o + sz
                if m2 + 8 > len(d):
                    break                       # mdat not here yet
                sz2 = int.from_bytes(d[m2:m2 + 4], "big")
                if sz2 < 8 or m2 + sz2 > len(d):
                    break                       # mdat still being written
                frag = d[o:m2 + sz2]
                sc = PL.trun_scan(frag)
                if sc is not None:
                    me = sc["media"]
                    p = me["doff"]
                    for s in me["sizes"]:
                        out.append(frag[p:p + s])
                        p += s
                    self.nfrag += 1
                    self.short += max(0, self.per_mpu - len(me["sizes"]))
                    self.fcounts.append(len(me["sizes"]))
                o = consumed = m2 + sz2
                continue
            o += sz
            consumed = o
        self.pos += consumed
        # A POISONED CURSOR CANNOT HEAL ITSELF. The shrink test above is only
        # a PROXY for roll detection: when the chain rolls the lane and the
        # NEW file grows past the old cursor between polls, size never dips
        # below pos, the cursor points mid-mdat in a DIFFERENT file, and the
        # box walk reads audio samples as a size field -- consuming nothing,
        # forever, on a file a fresh reader handles fine. That is precisely
        # how the 23:14 worker froze at 7 frames for 70 minutes while the
        # lane grew to 1229 fragments. Consuming nothing from a large unread
        # region N times running is not patience, it is the poisoned state.
        if out or consumed:
            self.no_progress = 0
        elif len(d) > (1 << 20):
            self.no_progress += 1
            if self.no_progress >= 8:
                self.force_roll(f"no progress at pos {self.pos} with "
                                f"{len(d)} unread bytes -- cursor poisoned")
        return out


def pick_audio(lns, pid=None, kind="audio"):
    # `kind` is matched EXACTLY: "audio" is the MMTP programme audio,
    # "route_audio" the ROUTE/DASH service's (E47) -- exact equality is
    # what keeps two services' lanes from ever being confused.
    cands = [l for l in lns.values() if l.get("kind") == kind]
    if pid is not None:
        cands = [l for l in cands if l.get("pid") == pid]
    # Deterministic: lowest pid is the primary programme audio. The
    # multiplex carries two `soun` assets and lock order is not stable.
    return sorted(cands, key=lambda l: l.get("pid", 1 << 30))[0] if cands else None


# THE 4 GiB WALL (E61, measured 8/10 00:17): RIFF sizes are 32-bit. At
# 48 kHz stereo the wav crosses 2^32 bytes after ~6.2 h and the header
# write raised OverflowError -- the worker DIED, close() re-raised the
# same overflow, and every relaunch that resumed the same wav died again
# within seconds (aud_eng3.err). Two defences, both required:
#   * the header CLAMPS at 0xFFFFFFFF so no code path can ever raise;
#   * the worker ROLLS the wav (archive + exit for respawn) at
#     WAV_ROLL_BYTES, well before the clamp could matter.
WAV_ROLL_BYTES = int(os.environ.get("ATSC3_WAV_ROLL_BYTES", 3_800_000_000))
ROLL_EXIT_RC = 42                 # exit code meaning "rolled, respawn me"


def archive_roll(out):
    """Rename-archive the wav + its sidecars to the next free .rollNNNN
    slot (NEVER delete under data/ -- the fleet law). The sidecar json is
    archived too: a stale sidecar pointing at a FRESH wav would map old
    slots onto new samples; with it gone the viewer plays silence until
    the respawned worker writes a true one. -> the suffix used, or None.

    The viewer opens these files briefly every chunk, so a rename can hit
    a sharing violation on Windows -- retry over ~10 s before giving up
    (a failed archive falls back to fresh-truncate, which the resume
    logic already treated as authoritative)."""
    base = os.path.splitext(out)[0]

    def pairs(tag):
        return [(out, f"{base}.{tag}.wav"),
                (base + ".state.json", f"{base}.{tag}.state.json"),
                (base + ".json", f"{base}.{tag}.json")]
    k = 1
    while k < 10000 and any(os.path.exists(d) for _, d in pairs(f"roll{k:04d}")):
        k += 1
    tag = f"roll{k:04d}"
    ok = True
    for src, dst in pairs(tag):
        if not os.path.exists(src):
            continue
        for _ in range(20):
            try:
                os.rename(src, dst)
                break
            except OSError:
                time.sleep(0.5)
        else:
            print(f"  ! could not archive {src} (still open elsewhere)")
            ok = False
    return tag if ok else None


class WavAppender:
    """A WAV whose header is rewritten as it grows, so it stays playable."""

    def __init__(self, path, ch=2, fs=OUT_FS, resume=0):
        self.path = path
        self.ch, self.fs = ch, fs
        self.n = resume
        if resume:
            self.fh = open(path, "r+b")
            self.fh.seek(44 + resume * ch * 2)
        else:
            self.fh = open(path, "wb")
        self._header()

    def bytes_on_disk(self):
        return 44 + self.n * self.ch * 2

    def _header(self):
        d = self.n * self.ch * 2
        # CLAMP, never raise: past 2^32 the header is wrong either way,
        # but an exception here killed the worker (E61) while a clamped
        # header merely under-reports length to naive readers. The roll
        # threshold keeps honest files from ever reaching the clamp.
        c1 = min(36 + d, 0xFFFFFFFF)
        c2 = min(d, 0xFFFFFFFF)
        h = (b"RIFF" + c1.to_bytes(4, "little") + b"WAVEfmt " +
             (16).to_bytes(4, "little") + (1).to_bytes(2, "little") +
             self.ch.to_bytes(2, "little") + self.fs.to_bytes(4, "little") +
             (self.fs * self.ch * 2).to_bytes(4, "little") +
             (self.ch * 2).to_bytes(2, "little") + (16).to_bytes(2, "little") +
             b"data" + c2.to_bytes(4, "little"))
        pos = self.fh.tell()
        self.fh.seek(0); self.fh.write(h)
        self.fh.seek(max(pos, len(h)))

    def append(self, pcm):
        pcm = np.clip(pcm, -1.0, 1.0)
        self.fh.write((pcm * 32767).astype("<i2").tobytes())
        self.n += len(pcm)
        self.fh.flush()
        self._header()
        self.fh.flush()

    def close(self):
        if self.fh.closed:
            return                # roll path already closed us (E61)
        self._header()
        self.fh.close()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--live-dir", default=None)
    ap.add_argument("--pid", type=int, default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--once", action="store_true")
    ap.add_argument("--interval", type=float, default=8.0,
                    help="seconds between catch-up passes")
    ap.add_argument("--follow", action="store_true",
                    help="start a pass as soon as the lane has grown instead "
                         "of every --interval seconds (which stays as the "
                         "longest wait). For a low-latency viewer: the "
                         "decoded sound is then a pass behind the lane, not "
                         "up to an interval plus a pass")
    ap.add_argument("--wait", type=float, default=90.0)
    ap.add_argument("--element", choices=("5_X", "pair"), default="5_X",
                    help="AC-4 element type: 5_X for the 5.1 main (pid13), "
                         "pair for the stereo simulcast (pid14, Spanish)")
    ap.add_argument("--kind", default="audio",
                    help="lane kind to follow: audio (MMTP programme, "
                         "default) or route_audio (the ROUTE/DASH service's "
                         "AC-4 lanes, E47 -- those are stereo "
                         "channel_pair_elements, so use --element pair)")
    ap.add_argument("--start-behind", type=float, default=0.0,
                    help="on a FRESH start, skip to this many seconds behind "
                         "the lane head instead of decoding from the lane's "
                         "beginning. A live viewer needs sound NEAR NOW; "
                         "decoding an hour of backlog at 1.17x means hours "
                         "of silence at the lag point (measured 8/08 00:26). "
                         "DVR renders want 0; live wants ~120.")
    ap.add_argument("--channels", type=int, choices=(2, 6), default=2,
                    help="2 = stereo pair, keeps up live (1.27x); 6 = full "
                         "5.1, measured 0.62x even threaded -- recordings "
                         "only until the filterbank is faster")
    ap.add_argument("--max-pass", type=int, default=1800,
                    help="most AC-4 frames decoded in one pass (1800 = 60 s; "
                         "bounds memory after a stall). With --once, a small "
                         "value replays a lane the way a live worker with a "
                         "short --interval sees it -- the multi-pass gate")
    ap.add_argument("--ref-peak", type=float, default=0.0,
                    help="decoder-scale level that is written at 0.9 of full "
                         "scale; 0 = the default for the element (see "
                         "REF_PEAK in the code)")
    ap.add_argument("--standby", action="store_true",
                    help="start, load the decoder, then wait for the real "
                         "argument list as one JSON line on stdin (an empty "
                         "line or end of input = not needed, exit 0). The "
                         "orchestrator starts its workers this way while the "
                         "chain is still tuning")
    a = ap.parse_args()
    if a.standby:
        # A worker started when its lane appears spends 1.4 s getting ready
        # (the interpreter, numpy, and scipy.signal alone is 0.94 s) while
        # the first two seconds of sound sit in the lane -- measured 11/10,
        # 2.0 s from the lane to a usable wav, on every channel change. In
        # standby that is paid while the chain tunes, which takes 12 s
        # anyway. What a worker decodes cannot be known until the chain has
        # read the multiplex, so the arguments come afterwards.
        import m17_ac4_walk                                     # noqa: F401
        import m37_render_hf                                    # noqa: F401
        import m38_resample                                     # noqa: F401
        import m30_filterbank                                   # noqa: F401
        import m42_ac4_stream                                   # noqa: F401
        line = sys.stdin.readline()
        if not line.strip():
            return 0
        a = ap.parse_args([str(x) for x in json.loads(line)])
    d = a.live_dir or os.path.join(ROOT, "data", "atsc3_live")

    t0 = time.time()
    lane = None
    while time.time() - t0 < a.wait:
        lane = pick_audio(lanes(d), a.pid, a.kind)
        if lane and os.path.exists(lane["path"]):
            break
        time.sleep(0.2)
    if not lane:
        print(f"no audio lane in {d} -- is the chain running with "
              f"--assets av|all ?", file=sys.stderr)
        return 1

    src = lane["path"]
    out = a.out or os.path.join(d, "live_audio.wav")
    print(f"audio lane pid {lane['pid']} -> {os.path.basename(src)}")
    print(f"  decoding to {out}")

    import m17_ac4_walk as W                                    # noqa: E402
    import m37_render_hf as R37                                 # noqa: E402
    import m38_resample as RS                                   # noqa: E402
    import m30_filterbank as FB                                 # noqa: E402
    from m42_ac4_stream import Ac4Stream                        # noqa: E402

    idx_path = os.path.splitext(src)[0] + ".idx"
    # RESUME. A worker that always starts at cursor 0 re-decodes the whole
    # lane: fine at 27 minutes, and after four hours it never catches up --
    # so the chain survives a restart and the audio does not. The state
    # needed is small and is written every pass.
    state_p = os.path.splitext(out)[0] + ".state.json"
    # E61: a wav already at/over the roll threshold must NEVER be resumed
    # (the header write dies) nor truncated (that deletes recorded air).
    # Archive it and start fresh -- exactly what the runtime roll does.
    try:
        if os.path.exists(out) and os.path.getsize(out) >= WAV_ROLL_BYTES:
            sz = os.path.getsize(out)
            tag = archive_roll(out)
            print(f"  wav at the 32-bit wall ({sz} bytes >= {WAV_ROLL_BYTES})"
                  f" -- archived as .{tag}, starting fresh")
    except OSError:
        pass
    resume_n = 0
    ref_resume = 0.0
    try:
        stv = json.load(open(state_p))
        _ch_ok = False
        if os.path.exists(out):
            try:
                import wave as _w
                _ch_ok = _w.open(out).getnchannels() == a.channels
            except Exception:                                  # noqa: BLE001
                _ch_ok = False
        # same path AND same generation: after a chain restart the lane has
        # the same name and different content, and the old cursor is not in it
        if (_ch_ok and stv.get("lane") == os.path.abspath(src)
                and stv.get("generation") in (None, lane.get("generation", 0))
                and os.path.getsize(out) >= 44 + stv["wav_samples"] * 2 * a.channels):
            resume_n = stv["wav_samples"]
            cursor0, pad0 = stv["cursor"], stv["pad_done"]
            # Rewind the LANE cursor a few frames and discard that much of
            # the first pass's output. Resuming cold would restart the
            # filterbank against zeros -- the very discontinuity the
            # resumable-synthesise gate exists to prevent -- so pay a few
            # frames to let the overlap and the A-SPX state settle.
            cursor0 = max(0, stv["cursor"] - 8)
            resume_skip = stv["cursor"] - cursor0
            # first_seq belongs to the wav's SAMPLE 0, not to wherever we
            # resumed. Recomputing it from the resumed cursor put it 3 slots
            # (6 s) late, which would have offset the whole muxed programme.
            resume_first = stv.get("first_seq")
            ref_resume = float(stv.get("ref_peak") or 0.0)
            print(f"  resuming: {resume_n/OUT_FS:.1f} s already written, "
                  f"lane frame {stv['cursor']} (re-running {resume_skip} "
                  f"frames to settle state)")
        else:
            cursor0 = pad0 = resume_skip = 0; resume_first = None
    except (OSError, ValueError, KeyError):
        cursor0 = pad0 = resume_skip = 0; resume_first = None
    skip0 = 0                       # fragments to skip on a FRESH start
    if not cursor0 and a.start_behind > 0:
        nfrag_now = count_frags(idx_path)
        skip0 = max(0, nfrag_now - int(a.start_behind / 2.002) - 1)
        if skip0:
            resume_skip = 0
            print(f"  fresh start {a.start_behind:.0f}s behind head: "
                  f"skipping {skip0} fragments")
            # E51: the frame count of those fragments is NOT skip0*60 --
            # that arithmetic mis-trimmed by 6497 frames (217 s) on the
            # 8/08 SDR-overflow afternoon, anchoring the whole wav minutes
            # off and playing stale sport under a live commercial. The
            # true trim is summed from what each fragment DECLARES, which
            # is only known once the reader has walked them (first pass).
    reader = LaneReader(src)
    pending = []                    # frames read but not yet decoded
    # SIX CHANNELS. The stream decoder always decoded the full 5.1 element;
    # until 8/07 this worker synthesized only the first pair and the rest
    # went to silence. Order L R C LFE Ls Rs = the standard WAV 6-channel
    # layout, so the file IS a 5.1 file without any channel map.
    chans = Ac4Stream.CHANS if a.channels == 6 else ("L", "R")
    acpl3_said = [False]
    wav = WavAppender(out, ch=a.channels, resume=resume_n)
    dec = Ac4Stream(element=a.element)
    LEAD = 4                        # frames of QMF lead-in, discarded
    MAX_PASS = max(8, a.max_pass)   # 60 s of audio; caps peak RSS
    # ONE LEVEL FOR THE WHOLE STREAM.
    #
    # Until 11/10 every pass was scaled so that ITS OWN loudest sample sat at
    # 0.9. A pass with a loud moment in it was therefore turned down and a
    # quiet one turned up: the level stepped at every pass boundary. Measured
    # against a one-pass decode of the same lane: 5-10 % steps with 8 s
    # passes, several dB with 2 s passes -- the sound "going quieter then
    # louder" that E116 put down to failed HF passes (hf_fail was 0).
    #
    # The decoder's own scale is fixed, so the gain is now a constant:
    # 0.9 / ref_peak. ref_peak starts at REF_PEAK and only ever goes UP, when
    # a pass is louder than anything before it, so nothing clips. It is
    # logged when it moves and kept in the state file, so a resumed worker
    # does not start louder than it stopped.
    #
    # REF_PEAK is set a little ABOVE the loudest sample seen in ~9 minutes of
    # this market's air (11/10: 5.1 element 8.23e5, stereo pair 1.69e5), so
    # that an ordinary programme never moves it. The first values (7.5e5 and
    # 1.8e5) were passed on air within 20 s: the gain dropped 0.8 dB at a
    # moment that depended on when the worker started, which is the very
    # run-to-run difference this is here to remove. The cost against the
    # per-pass scaling is a quieter track: about 3 dB on the 5.1, 2 dB on
    # the stereo. The decoder's true full scale is not established.
    #
    # OPEN: the stereo element sits 12 dB under the 5.1 on the same scale.
    # Whether that is the broadcast or a scale difference in the pair path is
    # not established; two defaults keep both tracks where they were.
    #
    # The reference follows the PROGRAMME, never a decode error. On air,
    # 11/10, three frames of the Ls/Rs pair came out 60 times louder than
    # anything before (see core_for_bursts): the reference followed them,
    # and the whole track played 36 dB down from that moment to the end of
    # the run. So it now moves only as far as BURST, twice the starting
    # value (6 dB: more than that above the loudest sample of 9 minutes is
    # not a programme); A-SPX frames above BURST play core-band; and
    # anything still above it is clipped where it stands and reported.
    REF_PEAK = {"pair": 2.0e5}.get(a.element, 9.0e5)
    REF_PEAK = a.ref_peak or REF_PEAK
    BURST = 2.0 * REF_PEAK
    ref_peak = min(max(REF_PEAK, ref_resume), BURST)
    hf_hot = {"lr": 0, "sr": 0}
    n_burst = 0
    # one resampler for the life of the worker: no seam at pass boundaries
    rs_stream = RS.StreamResampler()
    fb_state = {k: None for k in chans}
    # per-PAIR lead buffers: the two A-SPX pairs carry independent envelope
    # state and cannot share a lead-in
    lead_pcm = {k: np.zeros(0) for k in ("L", "R", "Ls", "Rs")}
    lead_groups = {"lr": [], "sr": []}
    # A-SPX state per pair (noise-table index + the delta-coded signal and
    # noise envelopes): as it stood at the START of each of the last LEAD
    # frames, then after the last frame. apply_hf_pair used to start every
    # call from nothing, so each pass (a) restarted the noise table at 0 --
    # the noise a frame got depended on where the pass happened to begin --
    # and (b) decoded TIME-coded envelopes against zeros until the next
    # FREQ-coded one, up to 9 frames (0.3 s) into the pass. Measured on a
    # 5.1 lane, 2 s passes against one pass: the 6-24 kHz band of L/R/Ls/Rs
    # was only -7..-12 dB alike. None = the start of a lane.
    hf_state = {"lr": [None], "sr": [None]}
    self_hf_fail = 0
    consumed_resume = not (cursor0 or skip0)
    pad_done = pad0
    n_padded = 0
    first_slot = resume_first
    anchor_frag = 0                 # idx position of the wav's sample-0 frag
    cursor = cursor0                # frames already consumed
    gen_seen = None                 # lane generation last observed
    total_wall = total_media = 0.0
    try:
        while True:
            t_pass = time.time()
            size_read = _size(src)      # before the read, so nothing that
            #                             lands during the pass is missed
            try:
                # IDENTITY-BASED roll detection, ahead of the byte-level
                # proxy: the writer publishes `generation` in live.json, and
                # a generation change IS a roll whether or not the new file
                # happens to be shorter than the cursor at poll time.
                ln_now = pick_audio(lanes(d), a.pid, a.kind)
                gen_now = ln_now.get("generation", 0) if ln_now else None
                if gen_seen is not None and gen_now is not None                         and gen_now != gen_seen:
                    reader.force_roll(f"generation {gen_seen} -> {gen_now}")
                if gen_now is not None:
                    gen_seen = gen_now
                # O(new), not O(total): only the bytes appended since the
                # last pass are read and parsed. Fragments still being
                # written are left for next time rather than mis-parsed.
                pending.extend(reader.read_new())
                if skip0 and reader.nfrag and not consumed_resume:
                    # FRESH start-behind (E51): trim by the frames the
                    # skipped fragments actually DECLARE (fcounts), never
                    # skip0*60 -- truncated fragments make those differ by
                    # minutes and every surplus frame trimmed mis-anchors
                    # the entire wav. Baseline pad_done to the SKIPPED
                    # region's shortfall only: the old reader.short
                    # baseline also swallowed the shortfall of fragments
                    # we DO consume from the first pass, silently
                    # un-padding them (the mirror image of the overpad).
                    trim, k, pad_done = fresh_anchor(
                        reader.fcounts, skip0, reader.per_mpu)
                    pending = pending[trim:]
                    cursor = trim
                    anchor_frag = k
                    consumed_resume = True
                    print(f"  anchored at fragment {k}: trimmed {trim} "
                          f"frames, baselined {pad_done} pad frames")
                elif cursor0 and reader.nfrag and not consumed_resume:
                    # a resumed worker must skip what the wav already holds
                    pending = pending[cursor0:]
                    consumed_resume = True
                frames = pending
                # Snapshot the fragment count AT THE SAME INSTANT as the
                # frames. Reading it after the decode instead -- ~19 s later
                # -- counts fragments the writer appended in the meantime,
                # which the frame snapshot cannot contain, so the "shortfall"
                # includes a whole pass of air that simply has not been read
                # yet. Measured: 1118 frames padded where the true shortfall
                # was 38. That inserts silence and drifts audio LATE, which
                # is worse than the drift it was added to cure.
                nfrag_now, true_short = reader.nfrag, reader.short
            except Exception as e:                             # noqa: BLE001
                print(f"  cannot read lane yet ({type(e).__name__})")
                frames, nfrag_now, true_short = [], reader.nfrag, reader.short
            if reader.rolled:
                # The chain restarted and rolled the lane. Everything the
                # cursor referred to belongs to the previous file, so keep the
                # wav (it is the record of what we decoded) and restart the
                # LANE side cleanly. Without this the worker holds a cursor
                # into a file that no longer contains it and never emits
                # another sample -- silently, while reporting `behind 0`.
                reader.rolled = False
                pending = list(frames)
                cursor = 0
                consumed_resume = True
                fb_state = {k: None for k in chans}
                lead_pcm = {k: np.zeros(0) for k in ("L", "R", "Ls", "Rs")}
                lead_groups = {"lr": [], "sr": []}
                hf_state = {"lr": [None], "sr": [None]}
                # `first_slot`, not `first_seq` -- the state FILE calls it
                # first_seq but the loop variable is first_slot, and resetting
                # the wrong name would silently leave the old lane's anchor in
                # place and mis-align every later slot.
                first_slot = None
                anchor_frag = 0     # new lane consumed from its fragment 0
                # pad_done tracks reader.short, which force_roll just reset
                # to 0 -- a stale pad_done blocks ALL shortfall padding until
                # short re-climbs past it, so every post-roll lane hole would
                # compress the wav timeline (audio drifting early). E45.
                pad_done = 0
                print(f"  lane roll #{reader.rolls}: audio cursor reset; "
                      f"wav keeps {wav.samples if hasattr(wav, 'samples') else '?'} "
                      f"samples of earlier air")
            behind = 0
            if len(pending) > 2:
                # TRULY INCREMENTAL: decode only the frames that are new.
                #
                # Resuming from the last i-frame was still O(n) per pass --
                # i-frames are far enough apart that `start` stayed near the
                # beginning and each pass re-decoded most of the programme.
                # Measured 0.49x real time over 320 s of air, so the worker
                # fell permanently behind after ~140 s.
                #
                # Three things carry state across a frame boundary and all
                # three are now held, so nothing needs re-deriving:
                #   fstate   A-SPX envelopes are delta-coded
                #   fb_state the MDCT overlap (m30.synthesise, gated
                #            bit-identical to the batch path)
                #   cfg      the A-SPX configuration, i-frames only
                # BOUND THE BATCH. Decoding every backlogged frame in one
                # pass allocates in proportion to the backlog: a cold restart
                # with 150k frames of lane peaked at 9.8 GB of spectra, and
                # the same happens after any stall. Cap the pass and let the
                # loop come round again -- progress is identical, peak memory
                # is bounded by MAX_PASS instead of by how far behind we are.
                new = pending[:MAX_PASS]
                behind = len(pending) - len(new)
                wins6, lens6, grp6, cnts6 = dec.decode_frames(new, chans)
                # SIX filterbanks in parallel. Sequential synthesis measured
                # 0.77x real time -- and a worker under 1.0x falls behind
                # FOREVER (E13). NumPy releases the GIL inside its kernels
                # (the same fact the whole m11 pipeline is built on), so the
                # six independent channels thread cleanly.
                # GROUPED BATCH synthesis (E39): channels sharing a
                # frame's framing run through ONE dct call. The threaded
                # version measured 0.62x -- the GIL serialised the many
                # small windows; batching amortises them instead.
                pcm, fb_state = FB.synthesise_frames(
                    wins6, lens6, cnts6, 1536,
                    states=fb_state, return_state=True)
                # 13.25: an ASPX_ACPL_3 programme (RF25's 5.1) is a two-
                # channel core plus A-CPL; the decoder gives L/R only.  The
                # wav keeps the channel count it was opened with, so the
                # four channels that do not exist as coded channels are
                # written silent (A-CPL upmix is not built).
                if len(pcm) < len(chans):
                    if not acpl3_said[0]:
                        acpl3_said[0] = True
                        print(f"  stream is ASPX_ACPL_3 (2-ch core + A-CPL): "
                              f"rendering L/R; {[k for k in chans if k not in pcm]} "
                              f"written silent", flush=True)
                    n_ref = len(next(iter(pcm.values())))
                    for k in chans:
                        if k not in pcm:
                            pcm[k] = np.zeros(n_ref)
                # Six independent framings must still agree on TOTAL samples
                # -- every frame is one frame long on every channel. If they
                # ever disagree, crop to the shortest and SAY so; a silent
                # crop here would be a slow A/V drift downstream.
                ns = {len(v) for v in pcm.values()}
                if len(ns) != 1:
                    n_min = min(ns)
                    print(f"  ! channel lengths disagree {sorted(ns)} -- "
                          f"cropping to {n_min}")
                    pcm = {k: v[:n_min] for k, v in pcm.items()}
                # A-SPX per PAIR: L/R rides group stream "lr", Ls/Rs rides
                # "sr", each with its own lead-in -- the QMF bank has filter
                # memory and the envelopes are delta-coded per pair. 4 frames
                # of lead is 133 ms against a 577-sample group delay. The
                # centre stays core-band until a mono HF path exists and is
                # gated (its 1-ch block IS parsed, keeping state warm).
                if dec.cfg is not None:
                    for (ka, kb), gkey in ((("L", "R"), "lr"),
                                           (("Ls", "Rs"), "sr")):
                        if ka not in chans:
                            continue
                        gl = grp6[gkey]
                        if any(g is not None for g in gl):
                            lead = min(len(lead_pcm[ka]), LEAD * 1536)
                            pa = np.concatenate([lead_pcm[ka][-lead:],
                                                 pcm[ka]])
                            pb = np.concatenate([lead_pcm[kb][-lead:],
                                                 pcm[kb]])
                            gg = lead_groups[gkey][-LEAD:] + gl
                            n_lead = len(gg) - len(gl)
                            hs, slog = hf_state[gkey], []
                            try:
                                hf, _n = R37.apply_hf_pair(
                                    pa, pb, gg, dec.cfg,
                                    state0=hs[max(0, len(hs) - 1 - n_lead)],
                                    state_log=slog)
                                pcm[ka], pcm[kb], nb_, hf_hot[gkey] = \
                                    core_for_bursts(
                                        hf["L"][lead:], hf["R"][lead:],
                                        pa, pb, lead, BURST, hf_hot[gkey])
                                if nb_:
                                    n_burst += nb_
                                    print(f"  ! A-SPX burst on {ka}/{kb}: "
                                          f"{nb_} frame(s) above "
                                          f"{BURST:.3g} play CORE-BAND "
                                          f"({n_burst} so far)", flush=True)
                                hf_state[gkey] = slog[-(LEAD + 1):]
                            except Exception as _hfe:          # noqa: BLE001
                                hf_state[gkey] = (
                                    hs + [hs[-1]] * len(gl))[-(LEAD + 1):]
                                self_hf_fail += 1
                                if self_hf_fail <= 5 or                                         self_hf_fail % 50 == 0:
                                    print(f"  ! HF pass failed "
                                          f"({type(_hfe).__name__}: "
                                          f"{str(_hfe)[:90]}) -- this pass "
                                          f"plays CORE-BAND only "
                                          f"(#{self_hf_fail})", flush=True)
                            lead_pcm[ka] = pa[-LEAD * 1536:]
                            lead_pcm[kb] = pb[-LEAD * 1536:]
                        else:
                            lead_pcm[ka] = pcm[ka][-LEAD * 1536:]
                            lead_pcm[kb] = pcm[kb][-LEAD * 1536:]
                            hf_state[gkey] = (
                                hf_state[gkey] + [hf_state[gkey][-1]]
                                * len(gl))[-(LEAD + 1):]
                        lead_groups[gkey] = (lead_groups[gkey] + gl)[-LEAD:]
                y = np.column_stack([pcm[k] for k in chans])
                rs0 = resume_skip
                if resume_skip:
                    y = y[resume_skip * 1536:]
                    resume_skip = 0
                # GAP, NEVER SHIFT -- AND IN PLACE, NEVER AT THE END.
                #
                # A wav has no timeline. The video container expresses a lost
                # MPU as a gap and holds the clock true either side of it; a
                # flat PCM file can only be SHORTER, so every frame the lane
                # is missing slides all subsequent audio earlier. E51 padded
                # the total at the BATCH END, which kept the aggregate true
                # but displaced everything after a hole within the batch --
                # audible as garble/cut-out/snap per pad event (E59, user-
                # reported on a dipping band). Each truncated fragment's
                # missing frames now land AT that fragment's position.
                ins = pad_insertions(reader.fcounts, anchor_frag, cursor,
                                     len(new), reader.per_mpu,
                                     y_first=cursor + rs0)
                if ins:
                    parts, prev = [], 0
                    for pos, miss in ins:
                        pos = min(pos, len(y) // 1536)
                        parts.append(y[prev * 1536:pos * 1536])
                        parts.append(np.zeros((miss * 1536, len(chans))))
                        n_padded += miss
                        pad_done += miss
                        prev = pos
                    parts.append(y[prev * 1536:])
                    y = np.vstack(parts)
                n_wav0 = wav.n          # wav offset where THIS batch begins
                if len(y):
                    pk = float(np.abs(y).max())
                    if pk > BURST:
                        n_over = int((np.abs(y) > BURST).any(axis=1).sum())
                        print(f"  ! level: {n_over} samples up to {pk:.4g}, "
                              f"above anything a programme reaches "
                              f"({BURST:.3g}) -- clipped, the reference "
                              f"stays at {ref_peak:.4g}", flush=True)
                        y = np.clip(y, -ref_peak, ref_peak)
                        pk = float(np.abs(y).max())
                    if pk > ref_peak:
                        print(f"  level: peak {pk:.4g} is above the "
                              f"reference {ref_peak:.4g} -- gain lowered "
                              f"{20 * np.log10(pk / ref_peak):.1f} dB from "
                              f"here on", flush=True)
                        ref_peak = pk
                    y = rs_stream.push(y * (0.9 / ref_peak))
                    wav.append(y)
                    total_media += len(y) / OUT_FS
                # Record WHICH MPU SLOT sample 0 of this wav belongs to.
                # Without it the wav is just "some audio" and the player can
                # only guess the offset -- which is how the first mux ended
                # up trimming with -shortest instead of aligning.
                # E45: also record WHERE IN THE WAV that slot lives. After a
                # lane roll the wav keeps every earlier generation's air, so
                # first_seq's fragment sits wav_base samples in, NOT at 0 --
                # the viewer that assumed sample 0 played 17-minute-old audio
                # against current video. wav_base makes the anchor roll-proof.
                if first_slot is None:
                    # E51: anchor by FRAGMENT position, never cursor//60 --
                    # truncated fragments make the two disagree by minutes.
                    first_slot = frag_seq(idx_path, anchor_frag)
                    # 13.26: how many channels the STREAM codes.  An
                    # ASPX_ACPL_3 programme is a stereo core (the other
                    # four are written silent above, A-CPL not built), and
                    # the viewer should mux it as the stereo it is rather
                    # than a 5.1 whose C/Ls/Rs/LFE are zeros.
                    json.dump(dict(first_seq=first_slot, fs=OUT_FS,
                                   pid=lane["pid"], path=out,
                                   wav_base=n_wav0, generation=gen_seen,
                                   channels_coded=(2 if acpl3_said[0]
                                                   else len(chans))),
                              open(os.path.splitext(out)[0] + ".json", "w"))
                cursor += len(new)
                del pending[:len(new)]
            try:
                json.dump(dict(lane=os.path.abspath(src), cursor=cursor,
                               pad_done=pad_done, wav_samples=wav.n,
                               first_seq=first_slot, generation=gen_seen,
                               ref_peak=ref_peak),
                          open(state_p, "w"))
            except OSError:
                pass
            dt_ = time.time() - t_pass
            total_wall += dt_
            # E116: hf_fail was counted and never SHOWN -- and each failed
            # apply_hf_pair mutes the whole pass's top octave, which the
            # listener hears as the sound "going quieter then louder".
            # A counter nobody reads is indistinguishable from zero.
            print(f"  {cursor} fr ({dec.n_bad} bad), "
                  f"{total_media:.1f} s of audio, "
                  f"pad {n_padded}, behind {behind}, "
                  f"hf_fail {self_hf_fail}, burst {n_burst}, "
                  f"pass {dt_:.1f} s"
                  + (f", {total_media/total_wall:.2f}x real time"
                     if total_wall > 0 else ""))
            # E61: ROLL before the 32-bit wall. Close first (we hold the
            # handle), archive wav+state+sidecar together, then exit with
            # the roll code -- the warden respawns us and the fresh start
            # re-anchors near the head with a true new sidecar. The viewer
            # sees the sidecar vanish and plays silence (never stale
            # audio, never a starve) until the new anchor lands.
            if wav.bytes_on_disk() >= WAV_ROLL_BYTES:
                print(f"  wav reached roll threshold "
                      f"({wav.bytes_on_disk()} >= {WAV_ROLL_BYTES}) -- "
                      f"archiving and exiting rc={ROLL_EXIT_RC} for respawn")
                wav.close()
                tag = archive_roll(out)
                print(f"  archived as .{tag}")
                return ROLL_EXIT_RC
            if a.once and not behind:
                break
            if a.follow and not a.once:
                # A pass on a clock is a pass at a random moment of the 2 s
                # an MPU takes to arrive: the sound for a slot was ready
                # anywhere from one pass (0.6 s) to an interval plus a pass
                # (2.6 s at --interval 2) after the slot landed, and which,
                # depended on when the worker happened to start. Measured
                # 11/10 with a muxer at --lag 0: one run held every chunk a
                # whole slot, the next put the first chunk out and the second
                # 4.3 s later. The lane is appended a whole fragment at a
                # time, so its size is the event to wait for; a fragment
                # caught half-written changes the size again when it is
                # complete.
                t_end = t_pass + max(0.5, a.interval)
                while (not behind and time.time() < t_end
                       and _size(src) == size_read):
                    time.sleep(0.05)
            elif not a.once:
                time.sleep(max(0.5, a.interval - dt_))
    except KeyboardInterrupt:
        print("  stopped")
    finally:
        wav.close()
    print(f"  wrote {out}  {total_media:.1f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
