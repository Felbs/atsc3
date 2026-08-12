#!/usr/bin/env python3
"""atsc3_tv_route.py -- ONE window: channel 132.1 (ROUTE/DASH), with sound.

The ROUTE sibling of atsc3_tv.py v2, for the SECOND service on RF33:
serviceId 1, "132.1", 1080p60 HEVC + AC-4 stereo (eng+spa), carried as
ROUTE/DASH on PLP1 -- where 107.1 is MMTP on PLP0.  The watchdog machinery
is IMPORTED from atsc3_tv, never forked: LeadGovernor (E37/E41),
RespawnGuard incl. the E46 stopped() path, VlcRemote, the one-window-law
kill/record helpers, and the E45 anchored audio mapping
(AudioTrack/wav_anchor/audio_slice/audio_ready) all run here verbatim.

LANE CONTRACT (written by the staged chain, `atsc3 watch --route IP:PORT`)
    live.json lanes with kind "route_video" / "route_audio" / "route_subs",
    pid = 100 + LCT tsi (110 video, 120 eng, 130 spa, 140 cc), one file per
    lane (live_route_video_pid110.m4s + .idx), idx `seq` = the DASH segment
    number == LCT TOI, one complete 2.002 s DASH segment per idx row --
    COMPLETE or absent, the chain never writes a holed ROUTE segment.
    live_route.json carries the SLS view: serviceId + per-Representation
    tsi/pid/kind/codecs/lang.  Kinds are matched EXACTLY everywhere, so
    these lanes are invisible to the MMTP viewer and vice versa.

WHAT IS DIFFERENT FROM THE MMTP VIEWER, AND WHY
  * NO retime happened chain-side: ROUTE segments are self-timed DASH
    fragments whose tfdt is the broadcaster's absolute clock (90 kHz UTC,
    ~4.5e8 s).  Each chunk REBASES tfdt by the session anchor (the lane's
    first segment) before muxing, so the mux sees small times and the v2
    itsoffset arithmetic carries over unchanged -- and the rebase is a
    constant subtraction, so a lane hole stays a HOLE (gap-true, the E20
    law).
  * NO per-sample NAL screen (v2's frag_nal_sane).  That screen exists for
    partially-REPAIRED MPUs, a class the ROUTE chain never emits: a lane
    segment is a complete LCT object off BCH-clean FEC blocks or it was
    dropped whole.  (m7_play.trun_scan also mis-indexes styp-prefixed DASH
    segments -- doff is moof-relative -- so reusing it here would silently
    drop every fragment.)
  * The audio lanes are AC-4 stereo channel_pair_elements (M14): the
    worker is `atsc3_audio.py --kind route_audio --element pair`, wavs
    live_route_audio.wav / live_route_audio_spa.wav.
  * Captions (tsi 40, stpp) are NOT muxed yet -- noted in E47, the stpp ->
    dvbsub bridge is the one v2 stage with no ROUTE twin.

ONE WINDOW LAW: this tool tracks its player in <live-dir>/_tv_route/ and
kills only that.  It does NOT touch the 107.1 viewer's player -- switching
services is: stop the MMTP viewer (and tvwatch, which would respawn it and
kills every vlc.exe when it fires), then start this.  Exact commands in
LIVE_TV_LAB.md E47.

Usage:
    python tools/atsc3_tv_route.py --live-dir data/e30            # live
    python tools/atsc3_tv_route.py --live-dir lab/e47_out/live \
        --replay --fast --player none --out gate.ts               # gate
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import atsc3_tv as TV                                             # noqa: E402

MPU_SECONDS = TV.MPU_SECONDS
SPF = TV.SPF
log = TV.log


# ---------------------------------------------------------------------------
# ISOBMFF helpers (read-only walks + the one constant-field patch)
# ---------------------------------------------------------------------------

def _boxes(d, start=0, end=None):
    end = len(d) if end is None else end
    o = start
    while o + 8 <= end:
        sz = int.from_bytes(d[o:o + 4], "big")
        if sz < 8 or o + sz > end:
            return
        yield o, sz, bytes(d[o + 4:o + 8])
        o += sz


def mdhd_timescale(init):
    """moov/trak/mdia/mdhd -> timescale, or None."""
    for o, sz, t in _boxes(init):
        if t != b"moov":
            continue
        for o2, sz2, t2 in _boxes(init, o + 8, o + sz):
            if t2 != b"trak":
                continue
            for o3, sz3, t3 in _boxes(init, o2 + 8, o2 + sz2):
                if t3 != b"mdia":
                    continue
                for o4, sz4, t4 in _boxes(init, o3 + 8, o3 + sz3):
                    if t4 == b"mdhd":
                        ver = init[o4 + 8]
                        p = o4 + 12 + (16 if ver else 8)
                        return int.from_bytes(init[p:p + 4], "big")
    return None


def first_tfdt(seg):
    """baseMediaDecodeTime of the segment's first traf, or None."""
    for o, sz, t in _boxes(seg):
        if t != b"moof":
            continue
        for o2, sz2, t2 in _boxes(seg, o + 8, o + sz):
            if t2 != b"traf":
                continue
            for o3, sz3, t3 in _boxes(seg, o2 + 8, o2 + sz2):
                if t3 == b"tfdt":
                    ver = seg[o3 + 8]
                    if ver:
                        return int.from_bytes(seg[o3 + 12:o3 + 20], "big")
                    return int.from_bytes(seg[o3 + 12:o3 + 16], "big")
    return None


def rebase_tfdt(seg, delta):
    """Subtract `delta` ticks from every tfdt in this segment.

    A fixed-width field patched in place -- no repack, every other byte
    identical.  The subtraction is the SAME constant for every segment of a
    session, so inter-segment spacing (including the spacing across a lost
    segment) is exactly what the broadcaster transmitted."""
    b = bytearray(seg)
    for o, sz, t in _boxes(b):
        if t != b"moof":
            continue
        for o2, sz2, t2 in _boxes(b, o + 8, o + sz):
            if t2 != b"traf":
                continue
            for o3, sz3, t3 in _boxes(b, o2 + 8, o2 + sz2):
                if t3 == b"tfdt":
                    ver = b[o3 + 8]
                    if ver:
                        v = int.from_bytes(b[o3 + 12:o3 + 20], "big")
                        b[o3 + 12:o3 + 20] = max(0, v - delta).to_bytes(8, "big")
                    else:
                        v = int.from_bytes(b[o3 + 12:o3 + 16], "big")
                        b[o3 + 12:o3 + 16] = max(0, v - delta).to_bytes(4, "big")
    return bytes(b)


# ---------------------------------------------------------------------------
# session
# ---------------------------------------------------------------------------

class RouteSession:
    """One continuous timeline over one route_video lane generation.

    Presents the ROUTE lanes under the plain kinds ("video"/"audio") so the
    IMPORTED atsc3_tv audio machinery (audio_slice, audio_ready, wav_anchor
    -- the whole E45 anchored mapping) runs on them unchanged.  MMTP lanes
    are EXCLUDED here for the same reason they exclude us: exact-kind
    matching is what keeps two services from ever crossing."""

    def __init__(self, live_dir, t_base):
        lj = TV.read_json(os.path.join(live_dir, "live.json")) or {}
        self.lanes = {k: {**l, "kind": l["kind"][len("route_"):]}
                      for k, l in (lj.get("lanes") or {}).items()
                      if str(l.get("kind", "")).startswith("route_")}
        self.vid = next((l for l in self.lanes.values()
                         if l.get("kind") == "video"), None)
        self.gen = self.vid.get("generation", 0) if self.vid else None
        self.t_base = t_base
        self.slot0 = None
        self.info = TV.read_json(os.path.join(live_dir, "live_route.json"))
        self.eng = TV.AudioTrack(live_dir, "live_route_audio.json")
        self.spa = TV.AudioTrack(live_dir, "live_route_audio_spa.json")
        self.tfdt0 = None             # session video-clock anchor (ticks)
        self.timescale = None

    def ok(self):
        return self.vid is not None

    def v_idx(self):
        return TV.read_idx(os.path.splitext(self.vid["path"])[0] + ".idx")

    def anchor(self, idx):
        """tfdt of the lane's FIRST segment + the track timescale, once."""
        if self.tfdt0 is not None or not idx:
            return self.tfdt0 is not None
        try:
            with open(self.vid["path"], "rb") as f:
                init = f.read(self.vid.get("init_bytes") or 0)
                seq0, off0, ln0 = idx[0]
                f.seek(off0)
                seg = f.read(ln0)
        except OSError:
            return False
        t0 = first_tfdt(seg)
        if t0 is None:
            return False
        self.tfdt0 = t0
        self.timescale = mdhd_timescale(init) or 90000
        log(f"  video clock anchored: tfdt0 {t0} @ {self.timescale} Hz "
            f"(segment {seq0})")
        return True


def write_chunk_video(sess, frags, path):
    """init + the chunk's segments, tfdt rebased to the session anchor."""
    with open(sess.vid["path"], "rb") as src, open(path, "wb") as dst:
        src.seek(0)
        dst.write(src.read(sess.vid.get("init_bytes") or 0))
        for seq, off, ln in frags:
            src.seek(off)
            dst.write(rebase_tfdt(src.read(ln), sess.tfdt0))
    return frags


def mux_chunk_route(sess, live_dir, frags, s_lo, s_hi, t_off, tmp_dir,
                    lane_seq0=None):
    """frags: [(seq, off, ln)] with s_lo <= seq < s_hi -> TS bytes.

    Same clock discipline as v2: the chunk's internal video time is read
    off the REBASED tfdt (exact, not assumed), the itsoffset subtracts it
    and adds the session clock; audio is chunk-local wav at plain t_off."""
    tmp_dir = os.path.abspath(tmp_dir)
    vm = os.path.join(tmp_dir, "chunk.m4s")
    write_chunk_video(sess, frags, vm)
    ae = os.path.join(tmp_dir, "chunk_eng.wav")
    asp = os.path.join(tmp_dir, "chunk_spa.wav")
    TV.write_wav(ae, TV.audio_slice(sess, live_dir, s_lo, s_hi, sess.eng))
    TV.write_wav(asp, TV.audio_slice(sess, live_dir, s_lo, s_hi, sess.spa))
    internal_t = 0.0
    if frags:
        try:
            with open(sess.vid["path"], "rb") as f:
                f.seek(frags[0][1])
                t = first_tfdt(f.read(frags[0][2]))
            if t is not None:
                internal_t = (t - sess.tfdt0) / float(sess.timescale or 90000)
        except OSError:
            pass
    v_extra = (frags[0][0] - s_lo) * MPU_SECONDS if frags else 0.0
    r = subprocess.run(
        [TV.ffbin("ffmpeg"), "-loglevel", "error",
         "-itsoffset", f"{t_off + v_extra - internal_t:.6f}", "-i", vm,
         "-itsoffset", f"{t_off:.6f}", "-i", ae,
         "-itsoffset", f"{t_off:.6f}", "-i", asp,
         "-map", "0:v", "-map", "1:a", "-map", "2:a",
         "-c:v", "copy",
         "-c:a", "mp2", "-b:a", "224k",
         "-metadata:s:a:0", "language=eng",
         "-metadata:s:a:1", "language=spa",
         # -copyts, as in mux_chunk_v2: without it ffmpeg re-bases every
         # chunk's clock to zero (measured here: DTS sawtoothed 0..2.002 s
         # per chunk in the appended TS), and the concat contract needs
         # every chunk to carry its LITERAL session timestamps.
         "-copyts", "-muxdelay", "0", "-muxpreload", "0",
         "-f", "mpegts", "pipe:1"],
        capture_output=True, cwd=tmp_dir)
    if r.returncode:
        log(f"  mux failed: {r.stderr[-200:].decode('ascii', 'replace')}")
        return b""
    return r.stdout


# ------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--live-dir", required=True)
    ap.add_argument("--lag", type=float, default=25.0,
                    help="seconds behind the write head")
    ap.add_argument("--lead", type=float, default=24.0,
                    help="media cushion (s) banked ahead of the tailing "
                         "player (see atsc3_tv.LeadGovernor)")
    ap.add_argument("--chunk", type=int, default=3, help="slots per chunk")
    ap.add_argument("--replay", action="store_true",
                    help="static dir: stream from the start at real time")
    ap.add_argument("--player", choices=("auto", "vlc", "none"),
                    default="auto")
    ap.add_argument("--vlc-headless", action="store_true")
    ap.add_argument("--vlc-http", default=None, metavar="PORT:PW")
    ap.add_argument("--audio-hold", type=float, default=75.0,
                    help="max seconds to hold an emit-ready chunk for its "
                         "eng audio (E45b); 0 disables")
    ap.add_argument("--out", default=None,
                    help="TS sink (default <live-dir>/_tv_route/live_tv_route.ts)")
    ap.add_argument("--fast", action="store_true",
                    help="no replay pacing (headless gates only)")
    ap.add_argument("--max-seconds", type=float, default=None)
    a = ap.parse_args()
    d = a.live_dir
    if a.player == "auto":
        a.player = "vlc"

    tmp = os.path.join(d, "_tv_route")
    os.makedirs(tmp, exist_ok=True)
    TV.kill_prior_player(tmp)      # only ever a player THIS tool spawned

    out_ts = a.out or os.path.join(tmp, "live_tv_route.ts")
    sink = open(out_ts, "wb")
    log(f"route sink: {out_ts}" + ("" if a.player == "none"
                                   else " (VLC follows first chunk)"))

    gov = TV.LeadGovernor(lead=a.lead) if a.player == "vlc" else None
    vlc_http = None
    if a.vlc_http:
        port_s, _, pw = a.vlc_http.partition(":")
        vlc_http = (int(port_s), pw)
    remote = None
    pl = None
    paused_by_us = False
    guard = TV.RespawnGuard(lead=a.lead, safety=10.0)

    def launch_vlc(start_time=None):
        nonlocal vlc_http, remote
        if vlc_http is None:
            import secrets
            vlc_http = (TV.free_port(), secrets.token_hex(8))
        p = TV.spawn_vlc(out_ts, "LIVE 132.1 (atsc3_tv_route)",
                         headless=a.vlc_headless, http=vlc_http,
                         start_time=start_time)
        if p is not None:
            remote = TV.VlcRemote(*vlc_http)
            TV.record_player(tmp, p, "vlc.exe")
            log(f"  VLC control interface on 127.0.0.1:{vlc_http[0]}")
        return p

    t_clock = 0.0
    sess = RouteSession(d, t_clock)
    log(f"  eng wav: {'yes' if sess.eng.path else 'NO (silence)'}   "
        f"spa wav: {'yes' if sess.spa.path else 'NO (silence)'}   "
        f"service: "
        f"{(sess.info or {}).get('service_id', '?')} "
        f"({len((sess.info or {}).get('reps', []))} reps)")
    cursor = None
    audio_wait0 = None
    t_wall0 = time.time()
    emitted_s = 0.0
    idle_rounds = 0
    n_fail = 0
    try:
        while True:
            if a.max_seconds and time.time() - t_wall0 >= a.max_seconds:
                log("max-seconds reached -- stopping")
                return 0
            # ---- player supervision: IMPORTED E37/E45/E46 machinery ----
            if pl is not None and gov is not None and remote is not None:
                now = time.time()
                t_meas, state = remote.status()
                why = None
                rew = False
                cushion = gov.lead_now(now)
                if t_meas is not None and state in ("playing", "paused"):
                    rew = gov.observe_playhead(
                        t_meas, now, paused=(state == "paused"))
                    cushion = gov.lead_now(now)
                    why = guard.reading(now, t_meas, state, rew, cushion)
                elif state == "stopped" and not paused_by_us:
                    rew = False
                    why = guard.stopped(now)          # E46 path, kept
                if why and not paused_by_us:
                    log(f"  {why} -- respawning at the live edge "
                        f"(playhead {t_meas}s, cushion {cushion:.0f}s)")
                    subprocess.run(["taskkill", "/PID", str(pl.pid),
                                    "/T", "/F"], capture_output=True)
                    st0 = max(0.0, gov.media - gov.lead)
                    pl = launch_vlc(start_time=st0)
                    if pl is not None:
                        gov.player_started(now, p0=st0)
                        guard.note_spawn(now)
                    paused_by_us = False
                elif rew and not paused_by_us and cushion <= gov.safety:
                    if remote.set_pause(True):
                        paused_by_us = True
                        gov.paused = True
                if pl is not None and not paused_by_us \
                        and gov.lead_now(now) <= gov.safety:
                    if remote.set_pause(True):
                        paused_by_us = True
                        gov.paused = True
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
            ns = RouteSession(d, t_clock)
            if not ns.ok():
                time.sleep(2.0)
                continue
            if sess.gen is None or ns.gen != sess.gen:
                if cursor is not None:
                    log(f"lane generation {sess.gen} -> {ns.gen}; "
                        f"re-anchoring at t={t_clock:.1f}s")
                    if gov is not None:
                        gov.reanchor()               # E41: same breath
                        paused_by_us = False
                sess = ns
                cursor = None
                audio_wait0 = None
            else:
                sess.vid = ns.vid or sess.vid
                sess.lanes = ns.lanes
                sess.eng, sess.spa = ns.eng, ns.spa
                sess.info = ns.info or sess.info
            idx = sess.v_idx()
            if not idx:
                time.sleep(2.0)
                continue
            if not sess.anchor(idx):
                time.sleep(2.0)
                continue
            head = idx[-1][0]
            if cursor is None:
                if a.replay:
                    cursor = idx[0][0]
                else:
                    lag_slots = max(2, int(a.lag / MPU_SECONDS))
                    cursor = idx[max(0, len(idx) - lag_slots)][0]
                sess.slot0 = cursor
            limit = head + 1 if a.replay else \
                head - int(a.lag / MPU_SECONDS) + 1
            made = 0
            while cursor + a.chunk <= limit:
                if a.max_seconds and time.time() - t_wall0 >= a.max_seconds:
                    break
                s_lo, s_hi = cursor, cursor + a.chunk
                if not a.replay and a.audio_hold > 0 \
                        and TV.audio_ready(sess, sess.eng, s_hi) is False:
                    if audio_wait0 is None:
                        audio_wait0 = time.time()
                    if time.time() - audio_wait0 < a.audio_hold:
                        break
                    log(f"  audio-hold expired ({a.audio_hold:.0f}s) at "
                        f"slot {s_lo} -- emitting with silence")
                audio_wait0 = None
                frags = [f for f in idx if s_lo <= f[0] < s_hi]
                ts = mux_chunk_route(sess, d, frags, s_lo, s_hi, t_clock,
                                     tmp, lane_seq0=idx[0][0])
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
                if pl is None and a.player == "vlc" \
                        and gov is not None and gov.want_player():
                    pl = launch_vlc()
                    if pl is None:
                        log("  VLC not found -- continuing headless "
                            "(TS keeps growing)")
                        a.player = "none"
                        gov = None
                    else:
                        gov.player_started(time.time())
                        guard.note_spawn(time.time())
                        log(f"  VLC window opened on the growing TS "
                            f"({gov.media:.1f}s cushion banked)")
                if a.replay and not a.fast:
                    cap = 6.0 + (gov.lead if gov is not None else 0.0)
                    ahead = emitted_s - (time.time() - t_wall0)
                    if ahead > cap:
                        time.sleep(ahead - cap)
            if made:
                idle_rounds = 0
            else:
                idle_rounds += 1
                if a.replay and idle_rounds >= 3:
                    if gov is not None and gov.buf:
                        try:
                            sink.write(gov._flush(time.time()))
                            sink.flush()
                        except OSError:
                            pass
                    log(f"replay complete: {emitted_s:.1f} s emitted, "
                        f"{n_fail} failed chunks")
                    # A replay SHORTER than the lead cushion never trips
                    # want_player() -- the 10 s gate dir would end with a
                    # perfect TS and no window. The dir is static now, so
                    # a plain VLC on the finished file is the right player
                    # and needs no governor at all.
                    if pl is None and a.player == "vlc" and emitted_s > 0:
                        sink.flush()
                        pl = launch_vlc()
                        if pl is not None:
                            log("  VLC window opened on the finished replay")
                    if pl is not None:
                        if paused_by_us and remote is not None:
                            remote.set_pause(False)
                        pl.wait()
                    return 0 if n_fail == 0 else 1
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
            sink.close()
        except OSError:
            pass
        if a.max_seconds and pl is not None and pl.poll() is None:
            log("  closing spawned player (one-window law cleanup)")
            subprocess.run(["taskkill", "/PID", str(pl.pid), "/T", "/F"],
                           capture_output=True)


if __name__ == "__main__":
    sys.exit(main())
