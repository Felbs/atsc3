#!/usr/bin/env python3
"""atsc3_tv.py -- ONE window: live picture WITH sound, from the lanes.

v2 (E35): USER-TOGGLEABLE TRACKS. The TS now carries
    video    HEVC, -c:v copy (no re-encode -- captions are no longer burned)
    audio 0  English  (live_audio.wav)     mp2 stereo, language=eng
    audio 1  Spanish  (live_audio_spa.wav) mp2 stereo, language=spa
             (silence for any slot the Spanish wav cannot cover -- a
             missing second language must never fail the chunk)
    subs     DVB subtitles (bitmap, via the PGS bridge in atsc3_subpix --
             ffmpeg cannot do text->bitmap, measured E35)
and the player is VLC, whose UI gives the user audio-track and subtitle
toggle buttons. VLC on this box CANNOT read a piped TS (stdin '-'/fd://0:
"cannot peek ... no demux modules matched", measured E35), so v2 appends
chunks to a growing .ts FILE and VLC tails that -- measured: VLC drains to
EOF, waits, and resumes when the file grows. --mode v1 keeps the old
ffplay+burn path unchanged.

WHY CHUNKS, WHY TS
    A growing MP4 cannot be tail-played reliably (moov/moof placement), and
    ffmpeg cannot mux "a bit more of the same file" incrementally. But
    MPEG-TS is a byte-concatenable stream by DESIGN -- it is what the air
    itself uses. So each K-slot span becomes its own tiny mux, offset onto
    the shared clock, and appended to the sink. Every chunk carries the
    same four streams in the same order, so the PID layout never changes
    mid-stream -- that uniformity is the concat contract.

THE CLOCK IS THE SLOT GRID (m39's rule, again)
    video chunk t = (slot - slot0) * 2.002 exactly; the audio slice for the
    SAME slots comes from the wav via the audio lane's own index, ANCHORED
    at the wav sidecar's first_seq (the workers may --start-behind, so wav
    sample 0 is first_seq's fragment, not the lane's fragment 0):
        wav_offset(slot) = (frag_index(slot) - frag_index(first_seq)) * SPF
    -- walking the slot RANGE, never the present fragments, so a lost MPU
    is a freeze and silence, not a shift (the 8/07 sync lesson).

    AUDIO MAPPING VALIDITY: the mapping is still generation-bound. After a
    future lane roll it is broken until the worker restarts and writes a
    fresh sidecar -- the sidecars carry no generation field yet. Known,
    not solved here.

CAPTION TIME -> LANE TIME (the v1 bug, closed)
    live.srt's clock is the SUBS WORKER's own start, not the lane clock --
    v1 missed every cue because of it. The worker now records the MPU slot
    of srt t=0 in live_subs.json {anchor_seq}; the exact conversion is
        cue_lane_t = cue_t + (anchor_seq - video_lane_seq0) * 2.002
    Fallback for dirs recorded before the sidecar existed: the subs lane's
    first_seq (right whenever the worker began at the lane's first doc).

ROLL-AWARE by construction: when the video lane's generation changes, the
session re-anchors -- the output clock keeps counting up from where it was,
so the player never sees time go backwards.

ONE WINDOW LAW: before spawning a player this tool kills the player IT
(or a previous run of it) spawned -- tracked by PID in _tv/player.pid and
verified by image name -- never any unrelated player on the box.

Usage:
    python tools/atsc3_tv.py --live-dir data/e29             # live v2, VLC
    python tools/atsc3_tv.py --live-dir data/live51 --replay # from start
    python tools/atsc3_tv.py --live-dir data/live51 --replay \
        --player none --out out.ts --fast                    # headless mux
    python tools/atsc3_tv.py --live-dir data/e29 --mode v1   # old ffplay
"""
from __future__ import annotations

import argparse
import base64
import json
import os
import secrets
import shutil
import socket
import subprocess
import sys
import time
import urllib.request
import wave

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

MPU_SECONDS = 2.002
SPF = 96096                       # audio samples per fragment at 48 kHz
IS_WIN = (os.name == "nt")
VLC_DEFAULT = r"C:\Program Files\VideoLAN\VLC\vlc.exe"   # Windows fallback


def player_image(name):
    """Platform image name for a player binary ('vlc' -> 'vlc.exe' on
    Windows, 'vlc' on Linux). One-window-law records compare against it."""
    return name + ".exe" if IS_WIN else name


def display_env():
    """Environment for spawning a WINDOWED player (E56, Linux port).

    The viewer usually runs nohup-detached from ssh, where DISPLAY/
    WAYLAND_DISPLAY are absent and a spawned VLC dies with no display.
    Resolve the session's GUI environment at SPAWN time -- the Xwayland
    auth-file suffix CHANGES per boot, so it can never be baked into a
    launcher script. Only MISSING variables are filled in; a desktop
    launch passes through untouched. Windows: env is inherited as-is."""
    env = dict(os.environ)
    if IS_WIN:
        return env
    uid = os.getuid()
    run = env.setdefault("XDG_RUNTIME_DIR", f"/run/user/{uid}")
    if "WAYLAND_DISPLAY" not in env \
            and os.path.exists(os.path.join(run, "wayland-0")):
        env["WAYLAND_DISPLAY"] = "wayland-0"
    env.setdefault("DISPLAY", ":0")
    if "XAUTHORITY" not in env:
        import glob
        auth = sorted(glob.glob(os.path.join(run, ".mutter-Xwaylandauth.*")),
                      key=os.path.getmtime, reverse=True)
        if auth:
            env["XAUTHORITY"] = auth[0]
    return env


def log(m):
    print(f"{time.strftime('%H:%M:%S')} {m}", flush=True)


def ffbin(name):
    p = shutil.which(name)
    if p:
        return p
    c = os.path.join(
        os.environ.get("LOCALAPPDATA", "C:"), "Microsoft", "WinGet",
        "Packages", "Gyan.FFmpeg_Microsoft.Winget.Source_8wekyb3d8bbwe",
        "ffmpeg-8.1-full_build", "bin", name + ".exe")
    return c if os.path.exists(c) else name


def read_json(p):
    try:
        return json.load(open(p))
    except (OSError, ValueError):
        return None


def read_idx(p):
    out = []
    try:
        for line in open(p):
            r = json.loads(line)
            out.append((r["seq"], r["off"], r["len"]))
    except (OSError, ValueError):
        pass
    return out


def resolve(live_dir, p):
    """Sidecar paths may be ROOT-relative; make them findable."""
    if not p:
        return None
    for cand in (p, os.path.join(ROOT, p),
                 os.path.join(live_dir, os.path.basename(p))):
        if os.path.exists(cand):
            return cand
    return None


class AudioTrack:
    """One decoded-wav track: sidecar json + wav + the lane that made it."""

    def __init__(self, live_dir, sidecar):
        self.first = self.path = self.pid = None
        self.wav_base = None          # sample offset of first_seq (E45)
        self.state = None             # worker state.json, the fallback anchor
        self.ch = 2
        am = read_json(os.path.join(live_dir, sidecar))
        self.gen = None
        if am:
            self.first = am.get("first_seq")
            self.pid = am.get("pid")
            self.wav_base = am.get("wav_base")
            self.gen = am.get("generation")
            self.path = resolve(live_dir, am.get("path"))
        if self.path:
            st = read_json(os.path.splitext(self.path)[0] + ".state.json")
            if st and "cursor" in st and "wav_samples" in st:
                self.state = st
            try:
                self.ch = wave.open(self.path).getnchannels()
            except Exception:                                  # noqa: BLE001
                self.path = None


# E95: which lane family to play.  One ATSC 3.0 carrier can deliver several
# SERVICES over DIFFERENT transports at once -- RF25 carries WNUV 54.1 over
# MMTP (kind "video") and WBFF's 145.100 over ROUTE/DASH (kind "route_video")
# simultaneously.  Picking the first lane whose kind == "video" silently
# chooses the MMTP service, which is exactly how six hours of work got
# labelled Fox while decoding The CW.  Default unchanged; opt in explicitly.
VIDEO_KIND = "video"
KIND_SUFFIX = ""          # "" for MMTP lanes, "route_" for ROUTE lanes


class Session:
    """One continuous timeline over one lane generation."""

    def __init__(self, live_dir, t_base):
        lanes = (read_json(os.path.join(live_dir, "live.json")) or {}
                 ).get("lanes", {})
        self.lanes = lanes
        self.vid = next((l for l in lanes.values()
                         if l.get("kind") == VIDEO_KIND), None)
        self.gen = self.vid.get("generation", 0) if self.vid else None
        self.t_base = t_base          # output clock where this session began
        self.slot0 = None             # first slot actually emitted
        self.subs_seq0 = next((l.get("first_seq") for l in lanes.values()
                               if l.get("kind") == KIND_SUFFIX + "subs"), None)
        # caption anchor: exact if the subs worker recorded it (E35)
        sj = read_json(os.path.join(live_dir, "live_subs.json"))
        self.anchor_seq = sj.get("anchor_seq") if sj else None
        self.anchor_exact = self.anchor_seq is not None
        if self.anchor_seq is None:
            self.anchor_seq = self.subs_seq0
        self.eng = AudioTrack(live_dir, "live_audio.json")
        self.spa = AudioTrack(live_dir, "live_audio_spa.json")

    def ok(self):
        return self.vid is not None

    def v_idx(self):
        return read_idx(os.path.splitext(self.vid["path"])[0] + ".idx")


def _downmix(seg):
    """(n, 6) int16 [L R C LFE Ls Rs] -> (n, 2) int16 ITU-style."""
    import numpy as np
    f = seg.astype("f4")
    l = f[:, 0] + 0.707 * f[:, 2] + 0.707 * f[:, 4]
    r = f[:, 1] + 0.707 * f[:, 2] + 0.707 * f[:, 5]
    out = np.stack([l, r], 1) * 0.7
    return np.clip(out, -32768, 32767).astype("<i2")


def audio_slice(sess, live_dir, s_lo, s_hi, track):
    """Stereo int16 PCM for slots [s_lo, s_hi): slot-aligned, silence where
    the wav (or the whole track) is absent -- a missing second language
    must never fail the chunk.

    ANCHORED MAPPING (E35 addendum): the workers now support
    --start-behind, so wav sample 0 belongs to the sidecar's first_seq
    fragment, NOT to the lane's fragment 0. The map is
        wav_offset(slot s) = (frag_index(s) - frag_index(first_seq)) * SPF
    with frag_index() = the fragment's position in the audio lane's .idx.
    The old absolute k*SPF mapping reads past EOF on a mid-lane wav and
    produces pure silence (that is why the running v1 window lost sound).
    Slots before first_seq, lane holes, and slots past the decode head are
    silence. Still generation-bound: after a lane roll the mapping breaks
    until the worker restarts (sidecars carry no generation field yet)."""
    import numpy as np
    n_slots = s_hi - s_lo
    out = np.zeros((n_slots * SPF, 2), "<i2")
    if track.path is None or track.first is None:
        return out
    lane = next((l for l in sess.lanes.values()
                 if l.get("kind") == "audio" and l.get("pid") == track.pid),
                None)
    if lane is None:
        return out
    a_idx = read_idx(os.path.splitext(lane["path"])[0] + ".idx")
    a_slot = {s: j for j, (s, _, _) in enumerate(a_idx)}
    # WHERE IN THE WAV does fragment j live? wav_anchor answers, best
    # source first. E45 (8/08): the old "first_seq = wav sample 0" answer
    # is WRONG the moment the lane rolls -- the worker KEEPS the wav
    # across rolls and rewrites first_seq to the NEW lane's first slot,
    # so sample 0 is hours-old air. Measured live: 1212 s of wav, 183 s
    # of it from the current generation -- the viewer was playing audio
    # 17 MINUTES stale against current video ("can't even tell if it's
    # the same channel").
    #   1. sidecar wav_base: exact, written by a roll-aware worker;
    #   2. worker state.json: {cursor, wav_samples} pins frame `cursor`
    #      of the CURRENT lane to the wav's end, so fragment j (frames
    #      j*60..) sits at wav_samples - (cursor - j*60)*SPF/60;
    #   3. legacy sample-0 mapping (no roll has happened yet).
    anchor = wav_anchor(track, a_slot, lane.get("generation"))
    if anchor is None:
        return out
    base, jref, jmin = anchor
    try:
        w = wave.open(track.path)
        nf = w.getnframes()
        nch = w.getnchannels()
        for i, s in enumerate(range(s_lo, s_hi)):
            j = a_slot.get(s)
            if j is None or j < jmin:
                continue            # lane hole, or before the wav's anchor
            off = base + (j - jref) * SPF
            if off < 0 or off + SPF > nf:
                continue            # wav not decoded that far yet
            w.setpos(off)
            seg = np.frombuffer(w.readframes(SPF), "<i2").reshape(-1, nch)
            if len(seg) != SPF:
                continue
            if nch == 2:
                out[i * SPF:(i + 1) * SPF] = seg
            elif nch == 6:
                out[i * SPF:(i + 1) * SPF] = _downmix(seg)
        w.close()
    except Exception as e:                                     # noqa: BLE001
        # NEVER a silent swallow: an exception here writes SILENCE into
        # the programme, and a run where every chunk hits it produces a
        # perfectly-formed, perfectly-mute TS (the shape a blanket
        # `except: pass` hid until the TS-audibility gate caught it).
        log(f"  ! audio_slice(pid {track.pid}) fell back to silence: "
            f"{type(e).__name__}: {e}")
    return out


def wav_anchor(track, a_slot, lane_gen=None):
    """(base, jref, jmin) mapping fragment j -> wav sample offset
    base + (j - jref) * SPF, or None if the track carries no anchor.
    Shared by audio_slice (reading) and audio_ready (holding). E45.

    An anchor stamped with a generation is only valid for THAT lane
    generation: in the window between a chain roll and the worker's own
    roll detection (up to one worker pass) the state pair still describes
    the OLD lane, and trusting it would replay the stale-audio bug for
    those chunks. No generation stamp (pre-E45 worker) = accept as-is."""
    def gen_ok(g):
        return g is None or lane_gen is None or g == lane_gen
    j0 = a_slot.get(track.first)
    if track.wav_base is not None and j0 is not None and gen_ok(track.gen):
        return track.wav_base, j0, j0          # exact (roll-aware worker)
    if track.state is not None and gen_ok(track.state.get("generation")):
        stw, stc = track.state["wav_samples"], track.state["cursor"]
        return (stw - (stc * SPF) // 60, 0,
                a_slot.get(track.state.get("first_seq"), 0))
    if track.wav_base is None and track.state is None and j0 is not None:
        return 0, j0, j0                       # legacy sample-0 mapping
    return None


def audio_ready(sess, track, s_hi):
    """Does the decoded wav already cover the chunk ending at s_hi?
    None = unanswerable (no track/anchor -- do not hold on it).

    E45b (8/08, measured): with the anchor CORRECT the wav frontier
    trails the air by 20-60 s (worker interval + pass time + write
    cadence), while the mux lag is 25 s -- so a chunk muxed the moment
    it cleared the video lag baked SILENCE for slots the worker simply
    had not reached (every chunk of the 09:23 session: eng RMS -inf).
    The old stale-anchor path never showed this only because it read
    minutes-old samples that always existed. The muxer therefore HOLDS
    an emit-ready chunk until its audio exists, bounded by --audio-hold
    so a dead worker degrades to silence, never to frozen video."""
    if track.path is None:
        return None
    lane = next((l for l in sess.lanes.values()
                 if l.get("kind") == "audio" and l.get("pid") == track.pid),
                None)
    if lane is None:
        return None
    a_idx = read_idx(os.path.splitext(lane["path"])[0] + ".idx")
    a_slot = {s: j for j, (s, _, _) in enumerate(a_idx)}
    anchor = wav_anchor(track, a_slot, lane.get("generation"))
    if anchor is None:
        return None
    base, jref, _ = anchor
    js = [j for s, j in ((s, a_slot.get(s)) for s in (s_hi - 1, s_hi - 2))
          if j is not None]
    if not js:
        return None                  # lane hole at the boundary: nothing to wait for
    try:
        nf = wave.open(track.path).getnframes()
    except Exception:                                          # noqa: BLE001
        return None
    return base + (max(js) - jref) * SPF + SPF <= nf


def parse_srt(path):
    """[(begin_s, end_s, text)] from an srt file."""
    import re
    out = []
    try:
        txt = open(path, encoding="utf-8", errors="replace").read()
    except OSError:
        return out
    pat = re.compile(r"(\d\d):(\d\d):(\d\d)[,.](\d{3})")

    def sec(g):
        return (int(g[0]) * 3600 + int(g[1]) * 60 + int(g[2])
                + int(g[3]) / 1000.0)

    for blk in txt.split("\n\n"):
        lines = [l for l in blk.strip().split("\n") if l]
        if len(lines) < 3 or "-->" not in lines[1]:
            continue
        m = pat.findall(lines[1])
        if len(m) < 2:
            continue
        out.append((sec(m[0]), sec(m[1]), "\n".join(lines[2:])))
    return out


def srt_stamp(t):
    t = max(0.0, t)
    h, r = divmod(t, 3600)
    m, r = divmod(r, 60)
    return f"{int(h):02d}:{int(m):02d}:{int(r):02d},{int((r % 1) * 1000):03d}"


def chunk_srt(cues, lo_t, hi_t, path):
    """Cues intersecting [lo_t, hi_t) rebased to the chunk. -> count."""
    n = 0
    with open(path, "w", encoding="utf-8") as f:
        for b, e, txt in cues:
            if e <= lo_t or b >= hi_t:
                continue
            n += 1
            f.write(f"{n}\n{srt_stamp(b - lo_t)} --> "
                    f"{srt_stamp(min(e, hi_t) - lo_t)}\n{txt}\n\n")
    return n


def chunk_cues(cues_srt, sess, lane_seq0, s_lo, s_hi):
    """srt cues -> [(b, e, text)] in CHUNK time via the lane clock.

    cue_lane_t = cue_t + (anchor_seq - lane_seq0) * 2.002  (E35, exact
    when live_subs.json exists; subs-lane first_seq otherwise)."""
    if sess.anchor_seq is None or lane_seq0 is None:
        return []
    shift = (sess.anchor_seq - lane_seq0) * MPU_SECONDS
    lo_t = (s_lo - lane_seq0) * MPU_SECONDS
    span = (s_hi - s_lo) * MPU_SECONDS
    out = []
    for b, e, txt in cues_srt:
        bl, el = b + shift - lo_t, e + shift - lo_t
        if el <= 0 or bl >= span:
            continue
        out.append((max(0.0, bl), min(span, el), txt))
    return out


class LeadGovernor:
    """Keeps the tailing player's file-lead from ever pinning at zero.

    E35 believed VLC on a growing file "drains to EOF, waits, resumes".
    E37 MEASURED the truth: when its file read-ahead (~7 s) touches EOF
    -- one churn hole, one slow-chain stretch, one lane roll -- VLC
    3.0.23 REWINDS to zero and replays the session file with state still
    "playing", then crosses every PTS hole in real time. And a viewer
    that merely idles at EOF can never rebuild lead against a real-time
    feed: every later chunk lands on a starved player (the 6-s stutter
    metronome). Baseline run: 2 rewinds, 55% of the watch frozen.

    The governor (a) delays the window until `lead` seconds of media are
    banked in the sink, and (b) models the playhead's upper bound: it
    advances at 1x wall time and never past appended media, so
    `media - p_est` is a LOWER bound on the player's true lead. When that
    bound falls to `safety`, appends are WITHHELD in RAM until a full
    `lead` is banked again, then flushed as one burst -- the player,
    frozen at EOF anyway, resumes with a rebuilt cushion instead of
    metronoming at zero. `max_hold` age-flushes a partial buffer so a
    dying feed still drains to the screen.

    Pure logic, no I/O -- gated with synthetic schedules plus an
    ungoverned negative control in lab/gate_tv2.py."""

    def __init__(self, lead=24.0, safety=10.0, max_hold=90.0):
        # safety must beat VLC's file READ-AHEAD, not just its playhead:
        # measured (E37 run B) the EOF touch -- and the rewind -- fired
        # while the playhead still had ~7 s of lead. 10 s clears it with
        # the 1 s supervision period to spare.
        self.lead, self.safety, self.max_hold = lead, safety, max_hold
        self.media = 0.0          # media seconds appended to the sink
        self.play0 = None         # wall time the player opened
        self.p_est = 0.0          # playhead upper bound, media seconds
        self.t_last = None
        self.paused = False       # playhead known not to advance
        self.buf = []             # withheld [(ts_bytes, media_s, wall)]
        self.buf_media = 0.0
        self.rebuilding = False
        self.rebuilds = 0
        self.events = []          # human-readable state transitions

    def reanchor(self):
        """A lane roll happened: the session re-anchored and the output
        clock jumped across the gap. The governor's pause/buffer/cushion
        state all describe the OLD generation and would deadlock the new
        one -- E41 (8/08): a roll that landed while the governor was paused
        left it withholding the fresh generation's chunks forever, waiting
        for a cushion signal from a playhead stuck on content that stopped
        growing. A roll IS new data by definition, never something to wait
        out, so drop the withheld buffer and re-arm from empty. The VLC
        window keeps playing what it has; the next chunks flow straight
        through until the fresh cushion is banked."""
        self.buf = []
        self.buf_media = 0.0
        self.paused = False
        self.rebuilding = False
        # keep media/play0/p_est: the output clock is continuous, so the
        # player's timeline does not restart -- only the withholding does.
        self.events.append("lane roll -- governor re-armed, buffer cleared")

    def _advance(self, now):
        if self.play0 is not None and self.t_last is not None \
                and not self.paused:
            self.p_est = min(self.media, self.p_est + (now - self.t_last))
        self.t_last = now

    def lead_now(self, now):
        self._advance(now)
        return self.media - self.p_est

    def want_player(self):
        return self.play0 is None and self.media >= self.lead

    def player_started(self, now, p0=0.0):
        self.play0, self.t_last, self.p_est = now, now, p0
        self.paused = False

    def observe_playhead(self, t_meas, now, paused=False):
        """Feed a REAL playhead reading (VLC http). Replaces the model
        when sane. -> True iff the reading is the EOF-rewind signature
        (E37 measured: at EOF of the growing file VLC does NOT wait --
        it rewinds to zero and replays the session file)."""
        self._advance(now)
        self.paused = paused
        rewound = (self.p_est - t_meas > 5.0 and t_meas < 3.0)
        if not rewound and t_meas <= self.p_est + 2.0:
            # never trust a reading AHEAD of the model. E45: nor a reading
            # carrying the REWIND SIGNATURE -- E37's "small is conservative
            # for the lead" was exactly backwards there: clamping p_est to a
            # rewound/booting VLC's zeros made `media - p_est` read as the
            # WHOLE session file (the logged cushion 234s == media), which
            # disabled the pause guard and armed the respawn storm.
            self.p_est = min(float(t_meas), self.media)
        return rewound

    def push(self, ts, media_s, now):
        """One muxed chunk in -> bytes to append to the sink now
        (b'' while withholding, a multi-chunk burst on flush)."""
        self._advance(now)
        if self.play0 is None:
            self.media += media_s
            return ts
        if not self.rebuilding and self.media - self.p_est <= self.safety:
            self.rebuilding = True
            self.rebuilds += 1
            self.events.append(
                f"lead {self.media - self.p_est:.1f}s <= safety "
                f"{self.safety:.0f}s -- withholding appends to rebuild "
                f"a {self.lead:.0f}s cushion")
        if not self.rebuilding:
            self.media += media_s
            return ts
        self.buf.append((ts, media_s, now))
        self.buf_media += media_s
        if self.buf_media >= self.lead \
                or now - self.buf[0][2] > self.max_hold:
            return self._flush(now)
        return b""

    def tick(self, now):
        """Idle-loop call: age-flush so a stalled feed still drains."""
        self._advance(now)
        if self.buf and now - self.buf[0][2] > self.max_hold:
            return self._flush(now)
        return b""

    def _flush(self, now):
        out = b"".join(t for t, _, _ in self.buf)
        self.events.append(
            f"cushion rebuilt: flushing {len(self.buf)} chunks "
            f"({self.buf_media:.1f}s) in one burst")
        self.media += self.buf_media
        self.buf, self.buf_media, self.rebuilding = [], 0.0, False
        return out

    def pop_events(self):
        ev, self.events = self.events, []
        return ev


class RespawnGuard:
    """Decides when a playhead reading justifies killing/respawning VLC.

    E45 (8/08): three different lies reached the old watchdog as
    "playhead 0" and it treated every one as an EOF rewind:
      1. an OPENING VLC reports time=0/state=stopped for ~19 s (measured
         on the 103 MB session file) -- a boot, not a rewind; with the
         15 s cooldown shorter than the boot, the watchdog respawned VLC
         before it ever finished opening, forever (the 08:52-08:56 storm,
         one visible black flash per respawn);
      2. our own respawn's --start-time seek settling;
      3. a transient demux clock re-base (stale subtitle PTS, also E45).
    So a respawn now requires: a real state (main() only feeds playing/
    paused readings), spawn-grace expiry, CONFIRM consecutive rewind
    signatures, and a cooldown a boot cannot beat. The genuinely WEDGED
    player -- stalled at EOF of the growing file with a fat cushion,
    playhead motionless while state stays "playing" (measured 08:48-08:52,
    playhead pinned ~13 s for 3.7 min) -- earns a respawn too, after
    FROZEN_S of no motion, long enough to clear the ~20 s real-time
    crossings of fade-gap video holes (measured on the same file)."""

    GRACE = 25.0                  # s after any spawn: readings are boot noise
    COOLDOWN = 30.0               # s between respawns: longer than a boot
    CONFIRM = 3                   # consecutive rewind readings required
    FROZEN_S = 45.0               # playing-but-motionless tolerance
    STOPPED_CONFIRM = 4           # consecutive stopped readings after grace

    def __init__(self, lead, safety):
        self.lead, self.safety = lead, safety
        self.spawn_t = None
        self.last_respawn = 0.0
        self.rew_n = 0
        self.stop_n = 0
        self.frozen_t = None
        self.frozen_since = None

    def note_spawn(self, now):
        self.spawn_t = now
        self.last_respawn = now
        self.rew_n = 0
        self.stop_n = 0
        self.frozen_t = self.frozen_since = None

    def stopped(self, now):
        """A 'stopped' status from a live VLC process -> respawn reason or
        None. E46 (8/08, observed live 14:04): E45's fix filtered stopped-
        state readings out entirely (they are boot noise for ~19 s), which
        made a genuinely ENDED player invisible -- VLC touched EOF, ended
        playback, and sat stopped for minutes while the stream stayed
        healthy; no playing/paused reading would ever arrive to trigger
        the other paths. Boot noise is still excluded by GRACE; a run of
        stopped readings past it means playback ended and only a respawn
        brings it back. No cushion condition: a stopped player shows
        nothing, so the live-edge respawn (which the pause guard then
        manages) is strictly better than dead air."""
        if self.spawn_t is not None and now - self.spawn_t < self.GRACE:
            return None
        self.stop_n += 1
        if (self.stop_n >= self.STOPPED_CONFIRM
                and now - self.last_respawn >= self.COOLDOWN):
            return f"player ended (stopped x{self.stop_n} readings)"
        return None

    def reading(self, now, t_meas, state, rewound, cushion):
        """One playing/paused reading in -> reason-string to respawn, or
        None. The caller owns the actual kill/relaunch (and must call
        note_spawn afterwards)."""
        if self.spawn_t is not None and now - self.spawn_t < self.GRACE:
            return None
        self.stop_n = 0           # a real playhead reading: not stopped
        can = (now - self.last_respawn >= self.COOLDOWN
               and cushion > self.safety)
        if rewound:
            self.rew_n += 1
            if self.rew_n >= self.CONFIRM and can:
                return f"EOF rewind confirmed ({self.rew_n} readings)"
            return None
        self.rew_n = 0
        if state == "playing":
            if t_meas != self.frozen_t:
                self.frozen_t, self.frozen_since = t_meas, now
            elif (self.frozen_since is not None
                  and now - self.frozen_since >= self.FROZEN_S
                  and cushion > self.lead and can):
                return (f"playhead frozen {now - self.frozen_since:.0f}s "
                        f"at {t_meas}s with {cushion:.0f}s banked")
        else:
            self.frozen_t = self.frozen_since = None
        return None


class ClockTrim:
    """E50: the three-clock problem, corrected by feedback instead of seeks.

    The TS is stamped on the BROADCASTER's clock (media time advances at
    the transmitter's rate); VLC consumes at ITS clock (in practice the
    audio card's); the SDR crystal sits in between. Any ppm disagreement
    makes the cushion C(t) = media_appended - playhead drift LINEARLY, and
    until now the pause/rebuild machinery "fixed" that with a lurch every
    few minutes, forever -- a watchdog treating a deterministic rate error
    as a fault. The deterministic correction: estimate dC/dt by least
    squares over a rolling window of (wall, cushion) samples taken only
    during clean playback, and trim VLC's playback rate by that ppm.
    cushion' = media_rate - playhead_rate, so playhead_rate must move BY
    +dC/dt to zero the slope; a small level term recenters the cushion on
    its target over ~30 min. Discontinuities (pause, respawn, roll, hold)
    poison the fit, so the window resets on every one -- the controller
    only ever acts on evidence from an unbroken stretch of playback."""

    WINDOW = 1500.0       # s of samples behind the fit. VLC reports the
    #                       playhead in WHOLE SECONDS (+-0.5 s quantization);
    #                       least-squares slope noise scales as 1/T^1.5, so
    #                       resolving tens of ppm honestly needs tens of
    #                       minutes -- which is fine: crystal ppm is constant
    #                       over hours, the pauses cover the first minutes,
    #                       and the servo owns the long run.
    MIN_SPAN = 600.0      # first (coarse) correction after 10 min
    DEADBAND = 60e-6      # |adjust| under 60 ppm is measurement noise
    MAX_TRIM = 1500e-6    # never trim more than 1500 ppm total
    PLAUSIBLE = 400e-6    # E59: no pair of crystals disagrees by more --
    #                       a fit demanding it measured APPEND BURSTS
    #                       (deep-lag catch-up feeds media faster than 1x),
    #                       not clock skew. Measured 8/09: the first eval
    #                       under --lag 120 demanded >=1500 ppm and slammed
    #                       the clamp; VLC then time-stretched 0.15% for
    #                       the rest of the session. Implausible windows
    #                       are DISCARDED, never acted on.
    PERIOD = 120.0        # re-evaluate cadence
    LEVEL_GAIN = 1.0 / 3600.0   # cushion-error (s) -> ppm-equivalent pull

    def __init__(self, target):
        self.target = target        # desired cushion, s (the lead)
        self.samples = []           # (wall, cushion)
        self.rate = 1.0
        self.last_eval = 0.0
        self.pending = None         # first-engagement confirmation (E59)

    def reset(self, why=""):
        self.samples.clear()

    def hard_reset(self):
        """E59: a RESPAWNED VLC starts at rate 1.0 -- carrying the old
        trim forward would make the next correction a blind jump from a
        rate the player is not actually running."""
        self.samples.clear()
        self.rate = 1.0
        self.pending = None

    def sample(self, now, cushion):
        self.samples.append((now, cushion))
        t0 = now - self.WINDOW
        while self.samples and self.samples[0][0] < t0:
            self.samples.pop(0)

    def evaluate(self, now):
        """-> new rate if a trim should be applied, else None."""
        if now - self.last_eval < self.PERIOD:
            return None
        self.last_eval = now
        if len(self.samples) < 8:
            return None
        span = self.samples[-1][0] - self.samples[0][0]
        if span < self.MIN_SPAN:
            return None
        # least-squares slope of cushion vs wall time
        n = len(self.samples)
        mt = sum(t for t, _ in self.samples) / n
        mc = sum(c for _, c in self.samples) / n
        num = sum((t - mt) * (c - mc) for t, c in self.samples)
        den = sum((t - mt) ** 2 for t, _ in self.samples)
        if den <= 0:
            return None
        slope = num / den                     # s of cushion per s of wall
        level_err = self.samples[-1][1] - self.target
        adj = slope + self.LEVEL_GAIN * level_err
        if abs(adj) < self.DEADBAND:
            self.pending = None               # clean window: drop suspicion
            return None
        if abs(adj) > self.PLAUSIBLE:
            # contaminated window (append burst / catch-up ramp): discard,
            # touch nothing -- the next unbroken window decides. E59.
            self.samples.clear()
            return None
        if self.rate == 1.0:
            # FIRST engagement gates on two consecutive agreeing windows:
            # one window is one measurement, and the 8/09 clamp was one
            # measurement acted on. Same sign and within 3x = agreement.
            if self.pending is None or self.pending * adj <= 0 \
                    or not (1 / 3 <= adj / self.pending <= 3):
                self.pending = adj
                self.samples.clear()
                return None
            self.pending = None
        new = self.rate + adj
        new = max(1.0 - self.MAX_TRIM, min(1.0 + self.MAX_TRIM, new))
        if abs(new - self.rate) < 20e-6:
            return None
        self.rate = new
        # The window now holds data measured under the OLD rate; a fit that
        # spans a rate change measures neither. Start clean -- each window
        # then measures exactly the residual of the rate it played under,
        # so corrections converge instead of integrating into a limit
        # cycle (the E50 gate's no-hunting case caught precisely that).
        self.samples.clear()
        return new


def drop_stale_subs(ts, t_off):
    """Strip dvbsub PES packets stamped BEFORE this chunk's clock offset.

    MEASURED (E45): ffmpeg's -fix_sub_duration EOF-flush still emits the
    occasional held event with a stale, chunk-LOCAL, unoffset timestamp
    (0.151 s inside the t_off 6.006 chunk; 26.250 inside a ~36 s chunk) --
    a BACKWARD PTS in the appended stream, which a playing VLC can read as
    a timestamp discontinuity and re-base its clock to ~0 (the playhead-0-
    while-playing signature). Every legitimate v2 sub packet carries
    pts >= t_off by construction (chunk_cues clips to [0, span] then adds
    t_off; the dummy sits at exactly t_off), so pts < t_off - 0.25 s
    identifies the stale flush exactly. Only stream_id 0xBD (private_
    stream_1, dvbsub) packets are ever touched; video/audio/PSI pass
    through byte-identical. A packet whose PES header we cannot parse is
    KEPT -- never drop blind."""
    if not ts or len(ts) % 188:
        return ts
    lim = (t_off - 0.25) * 90000.0
    out = bytearray()
    dropping = set()
    for o in range(0, len(ts), 188):
        p = ts[o:o + 188]
        if p[0] != 0x47:
            out += p
            continue
        pid = ((p[1] & 0x1F) << 8) | p[2]
        if p[1] & 0x40:                        # payload_unit_start
            dropping.discard(pid)
            af = (p[3] >> 4) & 0x3
            q = 5 + p[4] if af == 0x3 else 4   # skip adaptation field
            if (af & 0x1 and q + 14 <= 188
                    and p[q:q + 3] == b"\x00\x00\x01" and p[q + 3] == 0xBD
                    and (p[q + 7] & 0x80)):
                b = p[q + 9:q + 14]
                pts = (((b[0] >> 1) & 7) << 30 | b[1] << 22
                       | (b[2] >> 1) << 15 | b[3] << 7 | b[4] >> 1)
                if pts < lim:
                    dropping.add(pid)
        if pid not in dropping:
            out += p
    return bytes(out)


_VENC = None


def video_encoder():
    """NVENC if the box has it, x264 otherwise -- probed ONCE. GPU stays
    optional (fleet law): the fallback is always real."""
    global _VENC
    if _VENC is None:
        r = subprocess.run([ffbin("ffmpeg"), "-loglevel", "error",
                            "-f", "lavfi", "-i", "color=black:s=64x64:d=0.1",
                            "-c:v", "h264_nvenc", "-f", "null", "-"],
                           capture_output=True)
        _VENC = ["-c:v", "h264_nvenc", "-preset", "p4", "-b:v", "6M"] \
            if r.returncode == 0 else \
            ["-c:v", "libx264", "-preset", "ultrafast", "-crf", "23"]
        log(f"  video encoder for burn-in: {_VENC[1]}")
    return _VENC


_SUP = None


def sup_renderer():
    global _SUP
    if _SUP is None:
        sys.path.insert(0, HERE)
        import atsc3_subpix
        _SUP = atsc3_subpix.SupRenderer()
    return _SUP


_TRUN_SCAN = None


def _trun_scan():
    """lab/m7_play.trun_scan, the fleet's fMP4 fragment referee -- the
    mdat holds media samples AND a hint-track tail, so only the first
    traf's trun says where the video bytes actually are (the same fact
    atsc3_audio leans on)."""
    global _TRUN_SCAN
    if _TRUN_SCAN is None:
        try:
            sys.path.insert(0, os.path.join(ROOT, "lab"))
            import m7_play
            _TRUN_SCAN = m7_play.trun_scan
        except Exception:                                      # noqa: BLE001
            _TRUN_SCAN = False
            log("  ! m7_play unavailable -- corrupt-fragment screen OFF")
    return _TRUN_SCAN


def frag_nal_sane(frag):
    """True iff every media SAMPLE (per the trun) tiles exactly into
    4-byte length-prefixed NAL units.

    WHY (E35): catchup slot 225703306 carries a corrupted HEVC NAL
    ("Invalid NAL unit size -950566461") banked off real churn. Burn mode
    survived it -- the DECODER tolerates damage -- but v2's -c:v copy path
    runs the mp4->annexB converter instead, which aborts the ENTIRE chunk
    mux with AVERROR_INVALIDDATA. So corruption is caught at chunk build
    and the fragment is dropped: a corrupt MPU becomes a LOST MPU (a
    freeze, gap held true) instead of a dead chunk.

    NOT a naive mdat walk: the mdat carries the hint traf's bytes after
    the media samples, so tiling the whole box flags EVERY fragment
    (measured: 388/388) -- only the trun's [doff, sizes] region is video.
    If the referee is unavailable the screen is OFF, never drop-blind."""
    scan = _trun_scan()
    if not scan:
        return True
    try:
        sc = scan(frag)
    except Exception:                                          # noqa: BLE001
        return False
    if sc is None or "media" not in sc:
        return False
    me = sc["media"]
    p = me["doff"]
    for s in me["sizes"]:
        q, e = p, p + s
        if e > len(frag):
            return False
        while q < e:
            if q + 4 > e:
                return False
            ln = int.from_bytes(frag[q:q + 4], "big")
            if ln < 1 or q + 4 + ln > e:
                return False
            q += 4 + ln
        p = e
    return True


def frag_tfdt(b):
    """Absolute media time (ticks) from a fragment's moof/traf/tfdt, or None.

    E71: this number is the fragment's OWN statement of when it belongs on
    the lane clock. The muxer used to INFER it as (seq - lane_seq0)*2.002
    and subtract that via -itsoffset, while ffmpeg added the real tfdt back
    -- so the moment the two disagree the video is stamped off by exactly
    their difference. Measured 8/10: they agreed to 3 ms for three hours,
    then diverged by 36.061 s (18.01 slots) after a lane roll, and every
    later chunk carried video 36 s ahead of its own audio."""
    o = 0
    while o + 8 <= len(b):
        sz = int.from_bytes(b[o:o + 4], "big")
        if sz < 8:
            return None
        if b[o + 4:o + 8] == b"moof":
            p = o + 8
            while p + 8 <= o + sz:
                s2 = int.from_bytes(b[p:p + 4], "big")
                if s2 < 8:
                    return None
                if b[p + 4:p + 8] == b"traf":
                    q = p + 8
                    while q + 8 <= p + s2:
                        s3 = int.from_bytes(b[q:q + 4], "big")
                        if s3 < 8:
                            return None
                        if b[q + 4:q + 8] == b"tfdt":
                            ver = b[q + 8]
                            if ver and q + 20 <= len(b):
                                return int.from_bytes(b[q + 12:q + 20], "big")
                            if not ver and q + 16 <= len(b):
                                return int.from_bytes(b[q + 12:q + 16], "big")
                            return None
                        q += s3
                p += s2
            return None
        o += sz
    return None


def lane_time(sess, frags, lane_seq0, tfdt0):
    """Seconds of lane clock already inside this chunk's first fragment.

    The mux subtracts this via -itsoffset because ffmpeg will add the
    fragment's own tfdt back; getting it wrong offsets the video against
    its own audio by exactly the error (E71: 36.061 s, live). So MEASURE
    it from the fragment when the fragment says, and only fall back to
    the sequence-number inference when it does not."""
    if not frags:
        return 0.0
    inferred = ((frags[0][0] - lane_seq0) * MPU_SECONDS
                if lane_seq0 is not None else 0.0)
    if tfdt0 is None:
        return inferred
    ts = getattr(sess, "timescale", None)
    if ts is None:
        ts = lane_timescale(os.path.splitext(sess.vid["path"])[0] + ".idx")
        sess.timescale = ts
    measured = tfdt0 / ts
    if abs(measured - inferred) > MPU_SECONDS:
        log(f"  ! lane clock disagrees with slot arithmetic: tfdt "
            f"{measured:.3f}s vs inferred {inferred:.3f}s "
            f"({measured - inferred:+.3f}s) -- using the fragment's own "
            f"tfdt so video stays on its audio")
    return measured


def lane_timescale(idx_path, default=90000.0):
    """Ticks per second on this lane's clock, from the idx's own dur field
    (dur is one MPU = MPU_SECONDS), so the tfdt is read in the units the
    writer actually used rather than an assumed 90 kHz."""
    try:
        for line in open(idx_path):
            dur = json.loads(line).get("dur")
            if dur:
                return dur / MPU_SECONDS
    except (OSError, ValueError):
        pass
    return default


def write_chunk_video(sess, frags, path):
    """-> (frags actually written, first kept fragment's tfdt ticks).

    NAL-corrupt fragments are dropped; the tfdt belongs to the first
    fragment that SURVIVED, which is the one whose time the chunk starts on."""
    kept = []
    tfdt0 = None
    with open(sess.vid["path"], "rb") as src, open(path, "wb") as dst:
        src.seek(0)
        dst.write(src.read(sess.vid["init_bytes"] or 0))
        for seq, off, ln in frags:
            src.seek(off)
            b = src.read(ln)
            if not frag_nal_sane(b):
                log(f"  dropping corrupt video fragment slot {seq} "
                    f"({ln} bytes) -- gap, never a dead chunk")
                continue
            if not kept:
                tfdt0 = frag_tfdt(b)
            dst.write(b)
            kept.append((seq, off, ln))
    return kept, tfdt0


def write_wav(path, pcm):
    ww = wave.open(path, "wb")
    ww.setnchannels(pcm.shape[1])
    ww.setsampwidth(2)
    ww.setframerate(48000)
    ww.writeframes(pcm.tobytes())
    ww.close()


MP2_SPF = 1152                    # samples per mp2 frame at 48 kHz


def frame_align(pcm_e, pcm_s, carry):
    """Align this chunk's audio to the mp2 FRAME GRID via a carry.

    E48 (8/08, measured): a chunk span is 288288 samples = 250.25 mp2
    frames. Feeding that per-chunk makes ffmpeg's encoder pad the final
    partial frame with 864 zero samples, so EVERY chunk carried 251
    frames = 6.024 s of audio in a 6.006 s slot: +2988 ppm of surplus,
    18 ms of double audio at every boundary (154/154 chunks measured).
    A player draining samples continuously runs ~180 ms/min away from
    the picture and glitches at each correction -- the "getting out of
    sync ... starting to skip" report. Invisible pre-E45 only because
    the audio content was minutes stale.

    The cure is exact arithmetic, not tolerance: prepend the samples
    deferred from the previous chunk, emit the largest multiple of 1152,
    defer the remainder. Chunks land 250/250/250/251 (1001 frames per
    24.024 s -- the 1001 shows up again), the encoder never pads, and
    boundaries butt-join sample-exact. -> (pcm_e, pcm_s, lead_samples):
    the audio itsoffset moves EARLIER by lead_samples/48000."""
    import numpy as np
    lead = len(carry["eng"]) if carry["eng"] is not None else 0
    if lead:
        pcm_e = np.vstack([carry["eng"], pcm_e])
        pcm_s = np.vstack([carry["spa"], pcm_s])
    keep = (len(pcm_e) // MP2_SPF) * MP2_SPF
    carry["eng"], carry["spa"] = pcm_e[keep:].copy(), pcm_s[keep:].copy()
    return pcm_e[:keep], pcm_s[:keep], lead


def mux_chunk_v2(sess, live_dir, frags, s_lo, s_hi, t_off, tmp_dir,
                 lane_seq0=None, cues_srt=None, carry=None):
    """frags: [(seq, off, ln)] with s_lo <= seq < s_hi -> TS bytes.

    THE FRAGMENTS ALREADY CARRY TIME. `retime` writes each fragment's
    absolute lane clock into its tfdt (gap-true), so a chunk's internal
    video time is (seq - lane_seq0) * 2.002 BEFORE any offset; the video
    itsoffset SUBTRACTS that lane time (the v1 double-count lesson).
    Audio and subs are chunk-local and get the plain session offset
    (frame-grid aligned via `carry`, E48)."""
    tmp_dir = os.path.abspath(tmp_dir)
    vm = os.path.join(tmp_dir, "chunk.m4s")
    frags, tfdt0 = write_chunk_video(sess, frags, vm)  # corrupt frags dropped
    ae = os.path.join(tmp_dir, "chunk_eng.wav")
    asp = os.path.join(tmp_dir, "chunk_spa.wav")
    pcm_e = audio_slice(sess, live_dir, s_lo, s_hi, sess.eng)
    pcm_s = audio_slice(sess, live_dir, s_lo, s_hi, sess.spa)
    a_lead = 0
    if carry is not None:
        pcm_e, pcm_s, a_lead = frame_align(pcm_e, pcm_s, carry)
    a_off = t_off - a_lead / 48000.0
    write_wav(ae, pcm_e)
    write_wav(asp, pcm_s)
    sup = os.path.join(tmp_dir, "chunk.sup")
    # ABSOLUTE sup times, NO -itsoffset on this input, -copyts on the mux.
    # Two measured ffmpeg subtitle-pipeline traps force this shape:
    #   * -itsoffset + -fix_sub_duration flushed every event still held at
    #     input EOF with ONE stale, unoffset timestamp (chunk 2's last
    #     seven events all stamped 5.412533 = 17.4245 minus the 12.012
    #     offset) -- later-chunk cues were silently LOST;
    #   * without -copyts, ffmpeg rebases the sup input by its own first
    #     event, so the per-chunk dummy always lands at pts 0 -- and a TS
    #     sprinkled with pts-0 packets breaks timestamp-bisection seeks
    #     (input -ss landed at the silent head; the audibility gate FAILED
    #     on a TS whose whole-file loudness was -25 dB).
    # With -copyts every stream carries its literal timestamps: video via
    # the tfdt-compensating itsoffset, audio via plain t_off, subs baked.
    with open(sup, "wb") as f:
        f.write(sup_renderer().sup_bytes(
            [(b + t_off, e + t_off, txt) for b, e, txt in
             chunk_cues(cues_srt or [], sess, lane_seq0, s_lo, s_hi)],
            base=t_off))
    # The VIDEO inside the chunk starts at (frags[0].seq - s_lo) into the
    # span (a leading lost MPU means audio starts first).
    lane_t = lane_time(sess, frags, lane_seq0, tfdt0)
    v_extra = (frags[0][0] - s_lo) * MPU_SECONDS if frags else 0.0
    r = subprocess.run(
        [ffbin("ffmpeg"), "-loglevel", "error",
         "-itsoffset", f"{t_off + v_extra - lane_t:.6f}", "-i", vm,
         "-itsoffset", f"{a_off:.6f}", "-i", ae,
         "-itsoffset", f"{a_off:.6f}", "-i", asp,
         # -fix_sub_duration closes the open-ended PGS cues (E35: without
         # it the dvbsub packets land 2^32 ms late -- measured); no
         # -itsoffset here, the sup carries absolute times (see above)
         "-fix_sub_duration", "-i", sup,
         "-map", "0:v", "-map", "1:a", "-map", "2:a", "-map", "3:s",
         "-c:v", "copy",
         "-c:a", "mp2", "-b:a", "224k",
         "-c:s", "dvbsub",
         "-metadata:s:a:0", "language=eng",
         "-metadata:s:a:1", "language=spa",
         "-metadata:s:s:0", "language=eng",
         "-copyts", "-muxdelay", "0", "-muxpreload", "0",
         "-f", "mpegts", "pipe:1"],
        capture_output=True, cwd=tmp_dir)
    if r.returncode:
        log(f"  mux failed: {r.stderr[-200:].decode('ascii', 'replace')}")
        return b""
    if not r.stdout:
        # E61: rc=0 with EMPTY output was the one silent way a chunk
        # could vanish -- cursor advanced, nothing logged, TS frozen.
        log(f"  mux produced EMPTY output rc=0 (frags {len(frags)}): "
            f"{r.stderr[-200:].decode('ascii', 'replace')}")
        return b""
    out = drop_stale_subs(r.stdout, t_off)
    if not out:
        log(f"  drop_stale_subs emptied the chunk ({len(r.stdout)} bytes in)")
    return out


def mux_chunk_v1(sess, live_dir, frags, s_lo, s_hi, t_off, tmp_dir,
                 lane_seq0=None, cues=None):
    """The validated v1 path (burned captions, single audio), unchanged.

    Kept as --mode v1: the fallback must stay real (fleet law)."""
    tmp_dir = os.path.abspath(tmp_dir)
    vm = os.path.join(tmp_dir, "chunk.m4s")
    frags, tfdt0 = write_chunk_video(sess, frags, vm)  # corrupt frags dropped
    aw = os.path.join(tmp_dir, "chunk.wav")
    pcm = audio_slice(sess, live_dir, s_lo, s_hi, sess.eng)
    write_wav(aw, pcm)
    lane_t = lane_time(sess, frags, lane_seq0, tfdt0)   # E71: measured
    v_extra = (frags[0][0] - s_lo) * MPU_SECONDS if frags else 0.0
    vargs = ["-c:v", "copy"]
    vf = []
    if cues is not None:
        # CAPTIONS ARE BURNED, and therefore EVERY chunk re-encodes -- a TS
        # that alternates copied HEVC with encoded H.264 changes codec
        # mid-stream and players give up. Uniformity is the contract.
        # (v1 keeps its historical subs-worker-starts-at-lane-0 assumption.)
        sub_shift = ((lane_seq0 - sess.subs_seq0) * MPU_SECONDS
                     if sess.subs_seq0 is not None else 0.0)
        lo_t = (s_lo - lane_seq0) * MPU_SECONDS + sub_shift
        hi_t = (s_hi - lane_seq0) * MPU_SECONDS + sub_shift
        n = chunk_srt(cues, lo_t, hi_t, os.path.join(tmp_dir, "chunk.srt"))
        vargs = list(video_encoder())
        if n:
            # ffmpeg runs with cwd=tmp_dir so the subtitles filter can name
            # a bare "chunk.srt" -- Windows drive colons break its parser
            vf = ["-vf", "subtitles=chunk.srt"]
    r = subprocess.run(
        [ffbin("ffmpeg"), "-loglevel", "error",
         "-itsoffset", f"{t_off + v_extra - lane_t:.6f}", "-i", vm,
         "-itsoffset", f"{t_off:.6f}", "-i", aw,
         "-map", "0:v", "-map", "1:a"] + vf + vargs +
        ["-c:a", "mp2", "-b:a", "224k",
         "-muxdelay", "0", "-muxpreload", "0",
         "-f", "mpegts", "pipe:1"],
        capture_output=True, cwd=tmp_dir)
    if r.returncode:
        log(f"  mux failed: {r.stderr[-200:].decode('ascii', 'replace')}")
        return b""
    return r.stdout


# ---------------------------------------------------------------- players

def _image_name(pid):
    """Image name for a PID, or None. PID-exact -- never a name census
    (the pgrep self-match trap, bitten three times fleet-wide).
    Windows: tasklist. Linux: /proc/<pid>/comm."""
    if not IS_WIN:
        try:
            return open(f"/proc/{pid}/comm").read().strip().lower()
        except OSError:
            return None
    try:
        r = subprocess.run(["tasklist", "/FI", f"PID eq {pid}", "/FO", "CSV",
                            "/NH"], capture_output=True, text=True)
        line = (r.stdout or "").strip().splitlines()
        if line and line[0].startswith('"'):
            return line[0].split('","')[0].strip('"').lower()
    except OSError:
        pass
    return None


def kill_tree(pid):
    """Kill a player process (and any children) by PID, both platforms.
    Windows: taskkill /T. Linux: psutil children-first, then the process
    group (players are spawned with start_new_session=True, so the pgid
    is theirs alone -- never ours)."""
    if IS_WIN:
        subprocess.run(["taskkill", "/PID", str(pid), "/T", "/F"],
                       capture_output=True)
        return
    import signal
    try:
        import psutil
        p = psutil.Process(pid)
        for c in p.children(recursive=True):
            try:
                c.kill()
            except psutil.Error:
                pass
        p.kill()
        return
    except Exception:                                          # noqa: BLE001
        pass
    for sender in (os.killpg, os.kill):
        try:
            sender(pid, signal.SIGKILL)
            return
        except (OSError, ProcessLookupError):
            continue


def kill_prior_player(tmp_dir):
    """ONE WINDOW LAW: kill only the player a run of THIS tool spawned,
    identified by recorded PID and verified by image name."""
    pf = os.path.join(tmp_dir, "player.pid")
    rec = read_json(pf)
    if not rec:
        return
    pid, want = rec.get("pid"), (rec.get("image") or "").lower()
    if pid and want in ("vlc.exe", "ffplay.exe", "vlc", "ffplay") \
            and _image_name(pid) == want:
        log(f"  one-window law: killing prior {want} pid {pid}")
        kill_tree(pid)
    try:
        os.remove(pf)
    except OSError:
        pass


def record_player(tmp_dir, proc, image, http=None):
    rec = {"pid": proc.pid, "image": image}
    if http is not None:
        rec["port"], rec["pw"] = http     # warden's player-progress probe
    try:
        json.dump(rec, open(os.path.join(tmp_dir, "player.pid"), "w"))
    except OSError:
        pass


class VlcRemote:
    """Minimal localhost http control for the spawned VLC (E37).

    Why it exists: MEASURED on this box (VLC 3.0.23, file input), a
    tailing VLC that touches EOF of the growing TS does not wait -- it
    REWINDS to zero and replays the session file, then crosses every
    PTS hole in real time. The remote lets the muxer (a) read the true
    playhead, (b) pause the player BEFORE it can touch EOF, (c) resume
    it once the lead cushion is rebuilt. All failures are soft: if the
    interface is unreachable the tool falls back to the modelled lead."""

    def __init__(self, port, pw):
        self.port, self.pw = port, pw
        self.fails = 0

    def _req(self, query=""):
        rq = urllib.request.Request(
            f"http://127.0.0.1:{self.port}/requests/status.json{query}",
            headers={"Authorization": "Basic " + base64.b64encode(
                f":{self.pw}".encode()).decode()})
        with urllib.request.urlopen(rq, timeout=1.0) as f:
            return json.load(f)

    def status(self):
        """-> (playhead_s or None, state or None)."""
        try:
            st = self._req()
            self.fails = 0
            return st.get("time"), st.get("state")
        except Exception:                                      # noqa: BLE001
            self.fails += 1
            return None, None

    def set_pause(self, want):
        try:
            self._req("?command=pl_forcepause" if want
                      else "?command=pl_forceresume")
            return True
        except Exception:                                      # noqa: BLE001
            return False

    def set_rate(self, rate):
        """E50: sub-percent playback-rate trim (VLC time-stretches, so a
        few hundred ppm is inaudible). Used by ClockTrim only."""
        try:
            self._req(f"?command=rate&val={rate:.6f}")
            return True
        except Exception:                                      # noqa: BLE001
            return False


def free_port():
    s = socket.socket()
    s.bind(("127.0.0.1", 0))
    p = s.getsockname()[1]
    s.close()
    return p


def spawn_vlc(out_ts, title, headless=False, http=None, start_time=None):
    """http=(port, pw) arms the localhost control interface; headless
    swaps the window for dummy interface/outputs (probes and gates)."""
    vlc = shutil.which("vlc") or \
        (VLC_DEFAULT if IS_WIN and os.path.exists(VLC_DEFAULT) else None)
    if not vlc:
        return None
    # file input, NOT stdin: VLC-on-Windows cannot demux a pipe (E35).
    # NOTE (E37): VLC does NOT tail a growing file across EOF -- at EOF it
    # rewinds to zero (measured). The caller must keep it off EOF via the
    # lead governor + VlcRemote pause.
    cmd = [vlc, "--file-caching=3000", "--meta-title", title]
    if not IS_WIN:
        # fresh Linux installs pop a first-run privacy/network modal over
        # the video; nobody is at the keyboard of a spawned player (E56).
        # NOT --no-qt-updates-notif: distro builds compile without update
        # checking and the option does not EXIST there -- VLC refuses to
        # start on an unknown option (measured, first spawn on the box).
        cmd += ["--no-qt-privacy-ask"]
    if headless:
        # --dummy-quiet only exists on the Windows build (it suppresses
        # the console window); the Linux build refuses unknown options
        cmd += ["-I", "dummy"] + (["--dummy-quiet"] if IS_WIN else []) \
            + ["--vout=vdummy", "--aout=adummy"]
    if http is not None:
        cmd += ["--extraintf", "http", "--http-host", "127.0.0.1",
                "--http-port", str(http[0]), "--http-password", http[1]]
    if start_time:
        cmd += [f"--start-time={start_time:.1f}"]
    cmd.append(out_ts)
    return subprocess.Popen(cmd, env=display_env(),
                            start_new_session=not IS_WIN)


def spawn_ffplay(title):
    return subprocess.Popen(
        [ffbin("ffplay"), "-hide_banner", "-loglevel", "error", "-nostats",
         "-window_title", title, "-i", "pipe:0"],
        stdin=subprocess.PIPE, env=display_env(),
        start_new_session=not IS_WIN)


# ------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--live-dir", required=True)
    ap.add_argument("--transport", choices=("mmtp", "route"), default="mmtp",
                    help="which SERVICE family to play when the carrier "
                         "delivers several. mmtp = lanes named video/audio/"
                         "subs; route = route_video/route_audio/route_subs. "
                         "One carrier can run both at once (RF25 = WNUV over "
                         "MMTP + WBFFMob over ROUTE), and the old behaviour "
                         "silently picked MMTP.")
    ap.add_argument("--lag", type=float, default=25.0,
                    help="seconds behind the write head (must exceed the "
                         "audio worker's pass interval)")
    ap.add_argument("--lead", type=float, default=24.0,
                    help="media cushion (s) the sink keeps ahead of the "
                         "tailing player: the window opens only once this "
                         "much is banked, and after any starvation the "
                         "muxer re-banks it before appending again "
                         "(E37: kills the post-churn EOF metronome)")
    ap.add_argument("--chunk", type=int, default=3, help="slots per chunk")
    ap.add_argument("--replay", action="store_true",
                    help="static dir: stream from the start at real time")
    ap.add_argument("--mode", choices=("v2", "v1"), default="v2",
                    help="v2 = soft subs + dual audio + VLC; "
                         "v1 = burned subs + ffplay (the validated fallback)")
    ap.add_argument("--player", choices=("auto", "vlc", "ffplay", "none"),
                    default="auto")
    ap.add_argument("--vlc-headless", action="store_true",
                    help="spawn VLC with dummy interface/outputs (probes "
                         "and gates; the control interface still runs)")
    ap.add_argument("--vlc-http", default=None, metavar="PORT:PW",
                    help="pin the player control interface (default: "
                         "random free port + random password)")
    ap.add_argument("--subs", choices=("soft", "burn", "none"),
                    default=None,
                    help="default: soft in v2, burn in v1")
    ap.add_argument("--audio-hold", type=float, default=75.0,
                    help="max seconds to hold an emit-ready chunk waiting "
                         "for its eng audio to be decoded (E45b: the wav "
                         "frontier trails the air by 20-60 s; muxing before "
                         "it arrives bakes permanent silence). 0 disables.")
    ap.add_argument("--out", default=None,
                    help="v2 TS sink path (default <live-dir>/_tv/live_tv2.ts)")
    ap.add_argument("--fast", action="store_true",
                    help="no replay pacing (headless gates only)")
    ap.add_argument("--max-seconds", type=float, default=None,
                    help="stop (and close the spawned player) after this "
                         "much wall time")
    ap.add_argument("--no-audio-starve-bypass", action="store_true",
                    help="pre-E61 behaviour: per-chunk audio-hold even "
                         "when the worker wav is dead (the E42 metronome; "
                         "kept real as the gate's negative control)")
    ap.add_argument("--exit-on-player-close", action="store_true",
                    help="pre-E61 behaviour: a vanished VLC stops the "
                         "viewer. Default now RESPAWNS the player at the "
                         "live edge (rate-gated) -- a stack meant to run "
                         "indefinitely never treats a dead player as "
                         "operator intent.")
    a = ap.parse_args()
    # E95: select the SERVICE FAMILY before anything reads a lane.
    global VIDEO_KIND, KIND_SUFFIX
    if getattr(a, "transport", "mmtp") == "route":
        KIND_SUFFIX, VIDEO_KIND = "route_", "route_video"
    d = a.live_dir
    if a.subs is None:
        a.subs = "soft" if a.mode == "v2" else "burn"
    if a.player == "auto":
        a.player = "vlc" if a.mode == "v2" else "ffplay"

    tmp = os.path.join(d, "_tv")
    os.makedirs(tmp, exist_ok=True)
    kill_prior_player(tmp)

    pl = None
    sink = None
    out_ts = None
    if a.mode == "v1":
        if a.player == "none":
            out_ts = a.out or os.path.join(tmp, "live_tv1.ts")
            sink = open(out_ts, "wb")
        else:
            pl = spawn_ffplay("LIVE 107.1 + SOUND (atsc3_tv v1)")
            record_player(tmp, pl, player_image("ffplay"))
            sink = pl.stdin
            log("ffplay window opened; muxing begins")
    else:
        out_ts = a.out or os.path.join(tmp, "live_tv2.ts")
        sink = open(out_ts, "wb")
        log(f"v2 sink: {out_ts}" + ("" if a.player == "none"
                                    else " (VLC follows first chunk)"))

    mux = mux_chunk_v2 if a.mode == "v2" else mux_chunk_v1

    # the governor exists exactly when a player tails the growing sink
    gov = LeadGovernor(lead=a.lead) \
        if (a.mode == "v2" and a.player == "vlc") else None

    vlc_http = None               # (port, password), allocated once
    if a.vlc_http:
        port_s, _, pw = a.vlc_http.partition(":")
        vlc_http = (int(port_s), pw)
    remote = None                 # VlcRemote for the spawned player
    paused_by_us = False
    guard = RespawnGuard(lead=a.lead, safety=10.0)   # E45 respawn gate
    trim = ClockTrim(target=a.lead)                  # E50 clock-rate servo
    player_deaths = []            # E61 stampede gate for vanished players
    player_gone = None            # wall time the player was found dead

    def launch_vlc(start_time=None):
        nonlocal vlc_http, remote
        if vlc_http is None:
            vlc_http = (free_port(), secrets.token_hex(8))
        p = spawn_vlc(out_ts, "LIVE 107.1 (atsc3_tv v2)",
                      headless=a.vlc_headless, http=vlc_http,
                      start_time=start_time)
        if p is not None:
            remote = VlcRemote(*vlc_http)
            record_player(tmp, p, player_image("vlc"), http=vlc_http)
            log(f"  VLC control interface on 127.0.0.1:{vlc_http[0]}")
        return p

    t_clock = 0.0                 # the shared output clock, only ever grows
    sess = Session(d, t_clock)
    if a.mode == "v2":
        log(f"  eng wav: {'yes' if sess.eng.path else 'NO (silence)'}   "
            f"spa wav: {'yes' if sess.spa.path else 'NO (silence)'}   "
            f"cc anchor: "
            f"{'exact' if sess.anchor_exact else 'first_seq fallback'}")
    cursor = None                 # next slot to emit
    audio_wait0 = None            # wall time we began holding for audio (E45b)
    a_frontier = {"sig": None, "t": time.time()}     # eng wav progress (E61)
    stall_probe = {"t": None, "head": None}          # E61 freeze forensics
    a_carry = {"eng": None, "spa": None}   # mp2 frame-grid carry (E48)
    t_wall0 = time.time()
    emitted_s = 0.0
    idle_rounds = 0
    n_fail = 0
    try:
        while True:
            if a.max_seconds and time.time() - t_wall0 >= a.max_seconds:
                log("max-seconds reached -- stopping")
                return 0
            if pl is not None and pl.poll() is not None:
                # E61 (8/10): a VANISHED player process was the uncovered
                # third state -- E45 covered the wedged player, E46 the
                # stopped one, but a vlc.exe that is simply GONE (crash,
                # OOM, external kill) ended the whole viewer here, and
                # during a fade nothing downstream ever restarted it:
                # dead air until the operator noticed. The player is an
                # OUTPUT; losing it must never stop the muxer. Respawn at
                # the live edge, cooldown-spaced, with a stampede gate
                # (the TURBO law: a player dying 8x/hr is a real fault --
                # stop respawning, keep muxing, warden owns recovery).
                rc = pl.returncode
                if a.exit_on_player_close or a.mode != "v2" \
                        or a.player != "vlc" or gov is None:
                    log(f"player closed (rc={rc}) -- stopping")
                    return 0
                now = time.time()
                # count the corpse ONCE (the first cut of this code left
                # pl pointing at the dead process and re-counted it every
                # loop second -- 8 'deaths' in 8 s, instant stand-down;
                # measured by gate leg 6, 01:00:41)
                pl = None
                player_gone = now
                player_deaths.append(now)
                player_deaths[:] = [t for t in player_deaths
                                    if now - t < 3600.0]
                log(f"  player process GONE (rc={rc}) -- "
                    f"death {len(player_deaths)} this hour")
                if len(player_deaths) >= 8:
                    log(f"  player died {len(player_deaths)}x in an hour "
                        f"-- STANDING DOWN respawn; muxing continues "
                        f"headless (warden owns recovery)")
                    player_gone = None
                    remote = None
                    gov = None            # exactly the --player none path
                    a.player = "none"
            if pl is None and player_gone is not None and gov is not None \
                    and a.player == "vlc":
                now = time.time()
                if now - guard.last_respawn >= guard.COOLDOWN:
                    log("  respawning player at the live edge")
                    st0 = max(0.0, gov.media - gov.lead)
                    pl = launch_vlc(start_time=st0)
                    if pl is not None:
                        player_gone = None
                        gov.player_started(now, p0=st0)
                        guard.note_spawn(now)
                        paused_by_us = False
                        trim.hard_reset()
                    else:
                        log("  VLC not found -- continuing headless")
                        player_gone = None
                        remote = None
                        gov = None
                        a.player = "none"
            # ---- player supervision (E37): the growing-file VLC REWINDS
            # to zero if it ever touches EOF (measured), so keep it off
            # EOF: observe the true playhead, pause at the safety floor,
            # resume once the governor rebuilt the cushion, and respawn
            # at the live edge if a rewind slipped through anyway.
            if pl is not None and gov is not None and remote is not None:
                now = time.time()
                t_meas, state = remote.status()
                # E45: a reading is only a PLAYHEAD when VLC is actually
                # playing (or paused by us). An OPENING VLC reports time=0
                # with state stopped for ~19 s (measured on the session
                # file); the old code fed those zeros to the model, the
                # cushion inflated to the whole file, and the watchdog
                # respawned VLC every 15 s while it was still booting.
                # RespawnGuard owns grace/confirmation/cooldown/wedge.
                why = None
                rew = False
                cushion = gov.lead_now(now)
                if t_meas is not None and state in ("playing", "paused"):
                    rew = gov.observe_playhead(
                        t_meas, now, paused=(state == "paused"))
                    cushion = gov.lead_now(now)
                    why = guard.reading(now, t_meas, state, rew, cushion)
                elif state == "stopped" and not paused_by_us:
                    # E46: a live vlc.exe reporting stopped past boot-grace
                    # has ENDED playback -- no playhead reading will ever
                    # arrive; only the guard's stopped path can rescue it.
                    rew = False
                    why = guard.stopped(now)
                if why and not paused_by_us:
                    log(f"  {why} -- respawning at the live edge "
                        f"(playhead {t_meas}s, cushion {cushion:.0f}s)")
                    kill_tree(pl.pid)
                    st0 = max(0.0, gov.media - gov.lead)
                    pl = launch_vlc(start_time=st0)
                    if pl is not None:
                        gov.player_started(now, p0=st0)
                        guard.note_spawn(now)
                    paused_by_us = False
                    trim.hard_reset()   # E59: the new VLC runs at rate 1.0
                elif rew and not paused_by_us and cushion <= gov.safety:
                    # rewind-at-EOF but nothing to seek to = a fade. Pause
                    # and let the cushion rebuild; do NOT thrash VLC.
                    if remote.set_pause(True):
                        paused_by_us = True
                        gov.paused = True
                        trim.reset("fade pause")
                elif state == "playing" and not paused_by_us and not rew:
                    # E50: only unbroken clean playback feeds the clock
                    # servo -- every discontinuity path resets it instead.
                    trim.sample(now, cushion)
                    nr = trim.evaluate(now)
                    if nr is not None and remote.set_rate(nr):
                        log(f"  clock-trim: cushion slope corrected, "
                            f"rate {nr:.6f} ({(nr - 1.0) * 1e6:+.0f} ppm)")
                if pl is not None and not paused_by_us \
                        and gov.lead_now(now) <= gov.safety:
                    if remote.set_pause(True):
                        paused_by_us = True
                        gov.paused = True
                        trim.reset("safety pause")
                        log(f"  lead {gov.lead_now(now):.1f}s <= safety "
                            f"{gov.safety:.0f}s -- player paused while "
                            f"the cushion rebuilds")
                elif pl is not None and paused_by_us and not gov.rebuilding \
                        and gov.lead_now(now) > gov.safety:
                    if remote.set_pause(False):
                        paused_by_us = False
                        gov.paused = False
                        log(f"  cushion {gov.lead_now(now):.1f}s -- "
                            f"player resumed")
            ns = Session(d, t_clock)
            if not ns.ok():
                time.sleep(2.0)
                continue
            if sess.gen is None or ns.gen != sess.gen:
                if cursor is not None:
                    log(f"lane generation {sess.gen} -> {ns.gen}; "
                        f"re-anchoring at t={t_clock:.1f}s")
                    # E41: the governor's pause/buffer state belongs to the
                    # OLD generation and would withhold the new one forever.
                    # Re-arm it in the SAME breath as the session re-anchor --
                    # the two self-heals must fire together, not each alone.
                    if gov is not None:
                        gov.reanchor()
                        trim.reset("lane roll")
                        paused_by_us = False
                else_note = "" if cursor is not None else " (no prior cursor)"
                if cursor is None and sess.gen is not None:
                    log(f"  session re-adopted: gen {sess.gen} -> {ns.gen}"
                        f"{else_note}")
                sess = ns
                cursor = None
                audio_wait0 = None      # hold state belongs to the old lane
                a_carry = {"eng": None, "spa": None}   # E48: so does this
            else:
                # refresh growing state, keep session anchors
                sess.vid = ns.vid or sess.vid
                sess.lanes = ns.lanes
                sess.eng, sess.spa = ns.eng, ns.spa
                if ns.anchor_exact and not sess.anchor_exact:
                    sess.anchor_seq = ns.anchor_seq
                    sess.anchor_exact = True
            idx = sess.v_idx()
            if not idx:
                time.sleep(2.0)
                continue
            # E61 -- THE E42 ROOT CAUSE, closed. The audio-hold (E45b) is
            # only honest while the WORKER IS ALIVE: with the worker dead
            # the old per-chunk hold turned the muxer into a metronome --
            # one 6 s chunk per 75 s hold expiry, forever. Downstream that
            # trickle could never rebank the governor's cushion, the sink
            # went stale >120 s while the lanes advanced, tvwatch declared
            # the E42 "poll wedge" and cycled viewer+VLC every ~7.5 min
            # (the 8/10 00:xx "vanishing VLC" was exactly this loop). A
            # hold is a bet that the wav frontier will advance; when the
            # frontier has been FLAT longer than the hold budget the bet
            # is void -- emit at full rate with silence and let holds
            # re-arm the moment the wav moves again. Act on ABSENCE (the
            # E46 law, applied to the worker).
            now_f = time.time()
            f_sig = None
            if sess.eng.path is not None:
                try:
                    f_sig = (sess.eng.path, os.path.getsize(sess.eng.path))
                except OSError:
                    f_sig = None
            if f_sig != a_frontier["sig"]:
                a_frontier = {"sig": f_sig, "t": now_f}
            eng_starved = (not a.no_audio_starve_bypass) and \
                (now_f - a_frontier["t"] > max(a.audio_hold, 90.0))
            head = idx[-1][0]
            if cursor is None:
                if a.replay:
                    cursor = idx[0][0]
                else:
                    lag_slots = max(2, int(a.lag / MPU_SECONDS))
                    cursor = idx[max(0, len(idx) - lag_slots)][0]
                sess.slot0 = cursor
                # fires once per session/roll in health; a STREAM of these
                # is the re-anchor-every-pass freeze (E61 hunt)
                log(f"  anchored at slot {cursor} (head {head}, "
                    f"gen {sess.gen})")
            # emit while a full chunk sits behind the lag point
            limit = head + 1 if a.replay else \
                head - int(a.lag / MPU_SECONDS) + 1
            made = 0
            while cursor + a.chunk <= limit:
                if a.max_seconds and time.time() - t_wall0 >= a.max_seconds:
                    break
                s_lo, s_hi = cursor, cursor + a.chunk
                # E45b: hold the chunk until its ENG audio has actually been
                # decoded -- bounded, so a dead worker means silence resumes
                # after --audio-hold, never frozen video. Spanish is not
                # waited on (a missing second language never holds a chunk).
                if not a.replay and a.audio_hold > 0 and not eng_starved \
                        and audio_ready(sess, sess.eng, s_hi) is False:
                    if audio_wait0 is None:
                        audio_wait0 = time.time()
                    if time.time() - audio_wait0 < a.audio_hold:
                        break        # let the audio frontier catch up
                    log(f"  audio-hold expired ({a.audio_hold:.0f}s) at "
                        f"slot {s_lo} -- emitting with silence")
                audio_wait0 = None
                frags = [f for f in idx if s_lo <= f[0] < s_hi]
                cues = (parse_srt(os.path.join(d, "live.srt"))
                        if a.subs in ("burn", "soft") else
                        ([] if a.mode == "v2" else None))
                ts = mux(sess, d, frags, s_lo, s_hi, t_clock, tmp,
                         lane_seq0=idx[0][0], cues_srt=cues,
                         carry=a_carry) \
                    if a.mode == "v2" else \
                    mux(sess, d, frags, s_lo, s_hi, t_clock, tmp,
                        lane_seq0=idx[0][0],
                        cues=cues if a.subs == "burn" else None)
                if ts:
                    out_b = gov.push(ts, a.chunk * MPU_SECONDS,
                                     time.time()) if gov is not None else ts
                    for ev in (gov.pop_events() if gov is not None else ()):
                        log(f"  {ev}")
                    try:
                        if out_b:
                            sink.write(out_b)
                            sink.flush()
                    except OSError:
                        log("sink closed -- stopping")
                        return 0
                else:
                    n_fail += 1
                t_clock += a.chunk * MPU_SECONDS
                emitted_s += a.chunk * MPU_SECONDS
                cursor = s_hi
                made += 1
                if (a.mode == "v2" and pl is None and a.player == "vlc"
                        and gov is not None and gov.want_player()):
                    pl = launch_vlc()
                    if pl is None:
                        log("  VLC not found -- continuing headless "
                            "(TS keeps growing)")
                        a.player = "none"
                        gov = None      # nothing tails the sink now
                    else:
                        gov.player_started(time.time())
                        guard.note_spawn(time.time())
                        log(f"  VLC window opened on the growing TS "
                            f"({gov.media:.1f}s cushion banked)")
                # replay pacing: never run ahead of a real-time watcher
                # (with a governed player, run ahead by the lead cushion)
                if a.replay and not a.fast:
                    cap = 6.0 + (gov.lead if gov is not None else 0.0)
                    ahead = emitted_s - (time.time() - t_wall0)
                    if ahead > cap:
                        time.sleep(ahead - cap)
            if os.environ.get("ATSC3_TV_DEBUG"):
                log(f"  dbg cursor={cursor} head={head} limit={limit} "
                    f"made={made} starved={eng_starved} "
                    f"wait0={audio_wait0 is not None} nfail={n_fail} "
                    f"gen={sess.gen} t_clock={t_clock:.1f}")
            # E61 STALL PROBE: silent in health; if no chunk lands for 60 s
            # WHILE the head advances, print the emit decision variables --
            # the exact state every freeze theory needed and no log had.
            if made:
                stall_probe = {"t": None, "head": head}
            elif stall_probe["t"] is None:
                stall_probe = {"t": time.time(), "head": head}
            elif time.time() - stall_probe["t"] > 60.0 \
                    and head != stall_probe["head"]:
                ar = audio_ready(sess, sess.eng, cursor + a.chunk) \
                    if cursor is not None else "n/a"
                log(f"  STALL-PROBE cursor={cursor} head={head} "
                    f"limit={limit} starved={eng_starved} "
                    f"wait0={audio_wait0} ready={ar} nfail={n_fail} "
                    f"gen={sess.gen} idx_n={len(idx)}")
                stall_probe = {"t": time.time(), "head": head}
            if made:
                idle_rounds = 0
            else:
                idle_rounds += 1
                if a.replay and idle_rounds >= 3:
                    # static dir fully drained (tail shorter than one chunk
                    # is never emitted -- the concat contract needs whole
                    # chunks)
                    if gov is not None and gov.buf:
                        try:
                            sink.write(gov._flush(time.time()))
                            sink.flush()
                        except OSError:
                            pass
                    log(f"replay complete: {emitted_s:.1f} s emitted, "
                        f"{n_fail} failed chunks")
                    if pl is not None:
                        # a player we paused must never be left paused
                        if paused_by_us and remote is not None:
                            remote.set_pause(False)
                        pl.wait()
                    return 0 if n_fail == 0 else 1
            # a stalled feed must still drain any withheld chunks
            if gov is not None:
                out_b = gov.tick(time.time())
                for ev in gov.pop_events():
                    log(f"  {ev}")
                if out_b:
                    try:
                        sink.write(out_b)
                        sink.flush()
                    except OSError:
                        log("sink closed -- stopping")
                        return 0
            time.sleep(0.2 if (a.replay and a.fast) else 1.0)
    finally:
        try:
            if a.mode == "v1" and pl is not None:
                pl.stdin.close()
            elif sink is not None:
                sink.close()
        except OSError:
            pass
        # --max-seconds means a supervised run (gates): close our window
        if a.max_seconds and pl is not None and pl.poll() is None:
            log("  closing spawned player (one-window law cleanup)")
            kill_tree(pl.pid)


if __name__ == "__main__":
    sys.exit(main())
