#!/usr/bin/env python3
"""gate_e45.py -- gates for the E45 fixes (stale audio anchor after lane
roll; VLC respawn storm on boot-phase zeros; stale backward-PTS subtitle
packets), each with a NEGATIVE CONTROL that runs the OLD behaviour on the
SAME input and must FAIL -- proving every assertion here can fail at all.

    python lab/gate_e45.py            -> PASS/FAIL per gate, exit 0 iff all

Pure offline: synthetic lanes/wav/TS, no SDR, no live dir touched.
"""
from __future__ import annotations

import json
import os
import shutil
import sys
import tempfile
import types

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tools"))

import atsc3_tv as tv                                           # noqa: E402

SPF = tv.SPF
FAILS = []


def check(name, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" +
          (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


# ---------------------------------------------------------------- gate 1
# audio anchor: a post-roll wav (old-generation air banked at the front)
# must map the CURRENT lane's slots to the CURRENT generation's samples.

def synth_roll_dir(d):
    """A live dir shaped like the measured 09:03 state: the wav holds
    5 fragments of OLD-generation air (constant -3000), then 10 fragments
    of current-lane audio (fragment k = constant (k+1)*100), while the
    sidecar says first_seq = the CURRENT lane's first slot -- exactly what
    the running worker writes after a roll."""
    seqs = list(range(1000, 1010))          # current lane: 10 fragments
    lane = os.path.join(d, "live_audio_pid13.m4s")
    open(lane, "wb").write(b"")             # audio_slice never reads it
    with open(os.path.join(d, "live_audio_pid13.idx"), "w") as f:
        for k, s in enumerate(seqs):
            f.write(json.dumps({"seq": s, "off": k, "len": 1}) + "\n")
    old = np.full((5 * SPF, 2), -3000, "<i2")
    cur = np.concatenate([np.full((SPF, 2), (k + 1) * 100, "<i2")
                          for k in range(10)])
    wavp = os.path.join(d, "live_audio.wav")
    tv.write_wav(wavp, np.vstack([old, cur]))
    json.dump({"first_seq": 1000, "fs": 48000, "pid": 13, "path": wavp},
              open(os.path.join(d, "live_audio.json"), "w"))
    # worker state as measured live: cursor = 600 frames (10 fragments)
    # of the current lane, wav_samples = the whole wav
    json.dump({"lane": lane, "cursor": 600, "pad_done": 0,
               "wav_samples": 15 * SPF, "first_seq": 1000},
              open(os.path.join(d, "live_audio.state.json"), "w"))
    sess = types.SimpleNamespace(lanes={
        "audio_pid13": {"kind": "audio", "pid": 13, "path": lane}})
    return sess


def gate_audio_anchor():
    print("gate 1: audio wav anchor survives a lane roll")
    d = tempfile.mkdtemp(prefix="e45a_")
    try:
        sess = synth_roll_dir(d)
        want = [300, 400, 500]              # fragments j=2,3,4

        # NEW code path (state-based re-anchor, no wav_base in sidecar --
        # i.e. the worker that is running RIGHT NOW, unrestarted)
        trk = tv.AudioTrack(d, "live_audio.json")
        pcm = tv.audio_slice(sess, d, 1002, 1005, trk)
        got = [int(pcm[i * SPF, 0]) for i in range(3)]
        check("state-based re-anchor returns current-generation audio",
              got == want, f"got {got} want {want}")

        # NEW code path with an E45 worker sidecar (wav_base explicit)
        json.dump({"first_seq": 1000, "fs": 48000, "pid": 13,
                   "path": os.path.join(d, "live_audio.wav"),
                   "wav_base": 5 * SPF, "generation": 9},
                  open(os.path.join(d, "live_audio.json"), "w"))
        trk = tv.AudioTrack(d, "live_audio.json")
        pcm = tv.audio_slice(sess, d, 1002, 1005, trk)
        got = [int(pcm[i * SPF, 0]) for i in range(3)]
        check("sidecar wav_base anchor returns current-generation audio",
              got == want, f"got {got} want {want}")

        # NEGATIVE CONTROL -- the OLD mapping (first_seq = wav sample 0)
        # on the SAME input: no state file, no wav_base = exactly the code
        # that shipped. It must return the OLD generation's -3000 junk,
        # proving this gate distinguishes right from wrong.
        json.dump({"first_seq": 1000, "fs": 48000, "pid": 13,
                   "path": os.path.join(d, "live_audio.wav")},
                  open(os.path.join(d, "live_audio.json"), "w"))
        os.remove(os.path.join(d, "live_audio.state.json"))
        trk = tv.AudioTrack(d, "live_audio.json")
        pcm = tv.audio_slice(sess, d, 1002, 1005, trk)
        got = [int(pcm[i * SPF, 0]) for i in range(3)]
        check("NEGATIVE CONTROL: legacy sample-0 mapping plays stale air",
              got == [-3000] * 3 and got != want, f"got {got}")
    finally:
        shutil.rmtree(d, ignore_errors=True)


# ---------------------------------------------------------------- gate 2
# governor model: a rewind-signature reading must not poison p_est.

def gate_observe():
    print("gate 2: rewound reading does not poison the lead model")
    g = tv.LeadGovernor(lead=24.0)
    g.media = 234.0
    g.player_started(1000.0, p0=210.0)
    rew = g.observe_playhead(0, 1001.0)
    lead = g.lead_now(1001.0)
    check("rewind flagged, p_est held, lead stays honest",
          rew and g.p_est > 200.0 and lead < 40.0,
          f"p_est {g.p_est:.0f} lead {lead:.0f}")

    # NEGATIVE CONTROL -- the OLD update rule (E37: accept any reading
    # <= p_est+2, "small is conservative") on the same reading: p_est
    # collapses to 0 and the lead inflates to the whole file -- the
    # logged cushion-234s==media signature that disabled every guard.
    p_est_old = 210.0
    t_meas = 0.0
    if t_meas <= p_est_old + 2.0:
        p_est_old = min(t_meas, 234.0)
    lead_old = 234.0 - p_est_old
    check("NEGATIVE CONTROL: old clamp inflates lead to whole file",
          p_est_old == 0.0 and lead_old > 200.0,
          f"old p_est {p_est_old:.0f} old lead {lead_old:.0f}")


# ---------------------------------------------------------------- gate 3
# respawn guard vs the measured boot trace: ~19 s of time=0/state=stopped
# after every (re)spawn, then honest playback.

def boot_trace():
    """(t, t_meas, state) readings shaped like the measured 08:52-08:56
    storm: a booting VLC starves its own http thread (readings FAIL for
    stretches, letting the old model's p_est climb off the clamped zero),
    and the readings that do land say time=0 with a non-playing state
    for ~19 s, then honest playback. state None = the http read failed."""
    tr = []
    for t in range(1, 20):                       # booting, http flaky
        tr.append((float(t), 0, None if t % 8 else "stopped"))
    tr += [(19.0 + k, 210 + k, "playing") for k in range(60)]
    return tr


def gate_respawn_guard():
    print("gate 3: no respawn on boot-phase zeros; real events still fire")
    g = tv.RespawnGuard(lead=24.0, safety=10.0)
    g.note_spawn(0.0)
    fired = 0
    for t, t_meas, state in boot_trace():
        if state not in ("playing", "paused"):    # main()'s state gate
            continue
        rew = (210.0 - t_meas > 5.0 and t_meas < 3.0)
        if g.reading(t, t_meas, state, rew, cushion=234.0):
            fired += 1
    check("boot zeros + healthy playback -> zero respawns", fired == 0,
          f"fired {fired}")

    # NEGATIVE CONTROL -- the OLD watchdog (no state gate, no grace, no
    # confirmation, cooldown 15 s) driven exactly like the shipped code:
    # every reading that ARRIVES is believed, whatever its state; failed
    # readings advance the model at 1x. On the same trace it respawns
    # during boot -- the measured storm, cushion == media at each firing.
    last_respawn, fired_old = 0.0, 0
    p_est, media = 210.0, 234.0
    for t, t_meas, state in boot_trace():
        p_est = min(media, p_est + 1.0)           # _advance, 1 s cadence
        if state is None:
            continue                              # http failed: no reading
        rew = (p_est - t_meas > 5.0 and t_meas < 3.0)
        if t_meas <= p_est + 2.0:
            p_est = min(float(t_meas), media)     # the old clamp
        cushion = media - p_est
        if rew and t - last_respawn >= 15.0 and cushion > 10.0:
            fired_old += 1
            last_respawn = t
            p_est = 210.0                         # respawned at live edge
    check("NEGATIVE CONTROL: old watchdog respawns during boot",
          fired_old >= 1, f"fired {fired_old}")

    # a GENUINE persistent rewind (VLC replaying from 0, state playing)
    # must still respawn once confirmed and out of cooldown
    g = tv.RespawnGuard(lead=24.0, safety=10.0)
    g.note_spawn(0.0)
    fired = sum(1 for t in range(40, 50)
                if g.reading(float(t), 1, "playing", True, cushion=100.0))
    check("confirmed real rewind still respawns", fired >= 1,
          f"fired {fired}")

    # the frozen-playhead wedge (stalled at EOF with a fat cushion)
    g = tv.RespawnGuard(lead=24.0, safety=10.0)
    g.note_spawn(0.0)
    fired = sum(1 for t in range(30, 130)
                if g.reading(float(t), 50, "playing", False, cushion=60.0))
    check("50s-frozen playhead with fat cushion respawns", fired >= 1,
          f"fired {fired}")


# ---------------------------------------------------------------- gate 4
# stale-subtitle filter: backward sub PTS dropped, everything else kept.

def mk_pkt(pid, stream_id=None, pts=None, pusi=False):
    p = bytearray(188)
    p[0] = 0x47
    p[1] = (0x40 if pusi else 0) | (pid >> 8)
    p[2] = pid & 0xFF
    p[3] = 0x10                                   # payload only
    if pusi and stream_id is not None:
        q = 4
        p[q:q + 3] = b"\x00\x00\x01"
        p[q + 3] = stream_id
        p[q + 6] = 0x80
        p[q + 7] = 0x80                           # PTS present
        p[q + 8] = 5
        t = int(pts * 90000)
        p[q + 9] = 0x21 | ((t >> 29) & 0x0E)
        p[q + 10] = (t >> 22) & 0xFF
        p[q + 11] = 0x01 | ((t >> 14) & 0xFE)
        p[q + 12] = (t >> 7) & 0xFF
        p[q + 13] = 0x01 | ((t << 1) & 0xFE)
    return bytes(p)


def sub_pts_list(ts):
    """All stream-id 0xBD PTS in packet order."""
    out = []
    for o in range(0, len(ts), 188):
        p = ts[o:o + 188]
        if p[0] != 0x47 or not (p[1] & 0x40):
            continue
        af = (p[3] >> 4) & 0x3
        q = 5 + p[4] if af == 0x3 else 4
        if (q + 14 <= 188 and p[q:q + 3] == b"\x00\x00\x01"
                and p[q + 3] == 0xBD and (p[q + 7] & 0x80)):
            b = p[q + 9:q + 14]
            out.append((((b[0] >> 1) & 7) << 30 | b[1] << 22
                        | (b[2] >> 1) << 15 | b[3] << 7 | b[4] >> 1)
                       / 90000.0)
    return out


def gate_sub_filter():
    print("gate 4: stale backward-PTS subtitle packets are dropped")
    t_off = 6.006
    chunk = b"".join([
        mk_pkt(256, 0xE0, t_off + 0.1, pusi=True),      # video: keep
        mk_pkt(257, 0xC0, t_off + 0.05, pusi=True),     # audio: keep
        mk_pkt(259, 0xBD, t_off, pusi=True),            # sub dummy: keep
        mk_pkt(259, 0xBD, 0.151, pusi=True),            # STALE sub: drop
        mk_pkt(259),                                    # its tail: drop
        mk_pkt(259, 0xBD, t_off + 2.0, pusi=True),      # later sub: keep
    ])
    got = tv.drop_stale_subs(chunk, t_off)
    ok_len = len(got) == 4 * 188
    subs = sub_pts_list(got)
    mono = all(b2 >= a2 for a2, b2 in zip(subs, subs[1:]))
    check("filter drops exactly the stale PES (+tail), keeps the rest",
          ok_len and mono and abs(min(subs) - t_off) < 1e-3,
          f"kept {len(got)//188}/6 pkts, sub pts {subs}")
    vid_kept = got[:188] == chunk[:188]
    check("video/audio packets byte-identical", vid_kept and
          got[188:376] == chunk[188:376])

    # NEGATIVE CONTROL -- the unfiltered chunk (the old code appended
    # r.stdout as-is) fails the same monotonicity assertion.
    subs_old = sub_pts_list(chunk)
    mono_old = all(b2 >= a2 for a2, b2 in zip(subs_old, subs_old[1:]))
    check("NEGATIVE CONTROL: unfiltered chunk has backward sub PTS",
          not mono_old, f"sub pts {subs_old}")


# ---------------------------------------------------------------- gate 5
# E45b: a chunk must WAIT for its audio, not bake silence -- and an anchor
# from the wrong lane generation must be rejected, not trusted.

def gate_audio_hold():
    print("gate 5: audio_ready holds a chunk until the wav covers it")
    d = tempfile.mkdtemp(prefix="e45b_")
    try:
        sess = synth_roll_dir(d)          # wav covers fragments 0..9
        trk = tv.AudioTrack(d, "live_audio.json")
        ready_late = tv.audio_ready(sess, trk, 1012 + 1)   # frag 12: not yet
        ready_ok = tv.audio_ready(sess, trk, 1008 + 1)     # frag 8: decoded
        # extend the idx so slot 1012 exists (air arrived, audio not decoded)
        with open(os.path.join(d, "live_audio_pid13.idx"), "a") as f:
            for k, s in enumerate(range(1010, 1014)):
                f.write(json.dumps({"seq": s, "off": 99, "len": 1}) + "\n")
        ready_late = tv.audio_ready(sess, trk, 1013)
        check("undecoded slot -> not ready; decoded slot -> ready",
              ready_late is False and ready_ok is True,
              f"late {ready_late} ok {ready_ok}")

        # NEGATIVE CONTROL -- the old no-hold path on the same input:
        # audio_slice for the undecoded slot returns pure silence, which
        # is exactly what every chunk of the 09:23 live session baked
        # (measured eng RMS -inf across the whole TS).
        pcm = tv.audio_slice(sess, d, 1012, 1013, trk)
        check("NEGATIVE CONTROL: muxing without the hold bakes silence",
              int(np.abs(pcm).max()) == 0)

        # generation mismatch: an anchor stamped for lane gen 8 must be
        # rejected when the lane says gen 9 (the roll window)
        st = json.load(open(os.path.join(d, "live_audio.state.json")))
        st["generation"] = 8
        json.dump(st, open(os.path.join(d, "live_audio.state.json"), "w"))
        sess.lanes["audio_pid13"]["generation"] = 9
        trk = tv.AudioTrack(d, "live_audio.json")
        a_idx = tv.read_idx(os.path.join(d, "live_audio_pid13.idx"))
        a_slot = {s: j for j, (s, _, _) in enumerate(a_idx)}
        anc = tv.wav_anchor(trk, a_slot, 9)
        check("stale-generation anchor rejected (roll window -> silence, "
              "never stale audio)", anc is None, f"anchor {anc}")
    finally:
        shutil.rmtree(d, ignore_errors=True)


if __name__ == "__main__":
    gate_audio_anchor()
    gate_observe()
    gate_respawn_guard()
    gate_sub_filter()
    gate_audio_hold()
    print(f"\n{'ALL GATES PASS' if not FAILS else 'FAILED: ' + ', '.join(FAILS)}")
    sys.exit(1 if FAILS else 0)
