# TV2_READY -- atsc3_tv v2 (toggleable tracks), validated 8/08 ~00:53
# + E37 glitch fix (lead governor + player supervision), 8/08 ~01:45

## STATE: ALL GATES PASS (incl. E37 governor cases). Ready to deploy;
NOT deployed. The production viewer (tvd watchdog on data/e29) was not
touched. Deployment is the orchestrator's call.

## THE E37 GLITCH FIX (why the user saw stutter/freezing)

MEASURED on a real-time replay of live51's churn: VLC on the growing TS
does NOT wait at EOF -- it REWINDS to zero and replays the session file
(state never leaves "playing"), then crosses every PTS hole in real
time. One churn stall was enough to throw the viewer minutes back and
freeze 55% of the watch. v2 now:

  * banks a `--lead` cushion (default 24 s) in the sink before the
    window opens;
  * models + measures the true playhead (VLC http control interface on
    localhost, random port/password) and PAUSES the player before its
    read-ahead can touch EOF (safety floor 10 s -- the EOF touch fires
    ~7 s of lead early, measured);
  * withholds appends during starvation and flushes them as one burst
    once a full cushion is rebuilt, then resumes the player -- churn
    becomes ONE honest pause instead of a 6-s stutter metronome;
  * if a rewind slips through anyway, detects the signature (playhead
    drop to ~0) and respawns the window at the live edge (<1 s). The
    pause is BEST-EFFORT (VLC's read-ahead after its startup probe can
    exceed the floor); the respawn watchdog is the guarantee.

A/B on the same scripted churn (hole + 0.4x slow-chain + lane roll),
E37: baseline 2 permanent rewind-throwbacks, 9 stutter events, 147 s
frozen of 267 s watched, metronome forever after churn; fixed: churn =
one honest pause then a 24 s burst + resume, recovery phase 0 stalls,
one auto-corrected 0.3 s respawn, viewer at the live edge (playhead 234
vs 179) through a generation roll. Viewer latency at the screen is now
~lag+lead (~50 s behind air).

## DEPLOY COMMAND (live, e29)

    python tools/atsc3_tv.py --live-dir data/e29

That is the whole thing: v2 mode, soft dvbsub captions, dual audio
(eng/spa), VLC on a growing TS are all defaults. The tool enforces the
one-window law itself (kills only the player recorded in
`<live-dir>/_tv/player.pid`, PID + image-name verified). Useful knobs:

    --lag 25            seconds behind the write head (default)
    --lead 24           media cushion banked ahead of the player (E37);
                        bigger = smoother through churn, more latency
    --mode v1           the old validated ffplay + burned-captions path
    --player none --out X.ts --fast    headless mux (what the gates run)
    --max-seconds N     supervised run: stops and closes its own window
    --vlc-headless      dummy-output VLC (probes/gates)
    --vlc-http P:PW     pin the player control interface (probes)

Stop the running v1 viewer first (or let the tool's pidfile kill do it if
that viewer was tool-spawned); tvd should be pointed at v2 only after the
orchestrator decides.

DEPLOY-TIME NOTE (8/08 ~01:50): the user's ORIGINAL glitchy window --
VLC pid 4608, spawned by the 00:56 v2 run on e29's now-frozen TS -- is
still open and ORPHANED (its player.pid is gone), so a fresh tool run
will NOT auto-kill it. Close it at deploy:
    taskkill /PID 4608 /T /F     (verify it is still the vlc.exe from
                                  00:56:27 before killing)
Also: e29 currently has no lanes (chain restarting into the pre-meteor
dead band); a deployed viewer will correctly idle until lanes return.

## WHAT SHIPPED

* `tools/atsc3_tv.py` v2 -- video `-c:v copy` (no re-encode), TWO mp2
  stereo tracks tagged eng/spa (spa emits silence wherever its wav has no
  coverage -- never fails a chunk), soft DVB subtitles, VLC player. v1
  path kept intact under `--mode v1`.
* `tools/atsc3_subpix.py` -- srt text -> PGS `.sup` bitmaps (PIL), the
  bridge that makes `-c:s dvbsub` possible (ffmpeg cannot text->bitmap).
* `tools/atsc3_subs.py` -- now records `live_subs.json {anchor_seq}` (the
  MPU slot of srt t=0) the moment it sets its anchor; atsc3_tv converts
  cue time to lane time exactly: `cue_lane_t = cue_t + (anchor_seq -
  video_lane_seq0) * 2.002`. Fallback for old dirs: subs-lane first_seq.
* `tools/atsc3_meter.py` -- headless proof: video null-sink fps + errors,
  per-track non-silent/clipping slot fractions (first_seq-ANCHORED wav
  mapping), CC cue count + anchor sanity. One JSON line + summary.
* `lab/gate_tv2.py` -- the full gate suite below, self-negative-controlled.

## GATE RESULTS (8/08 ~01:45 re-run, all PASS)

E37 additions: five lead-governor cases (steady pass-through, gap ->
withhold+burst, no re-pin after recovery, max-hold drain, and the
ungoverned-pass-through negative control that reproduces the metronome
50/50). NOTE a gate-corpus drift found during the re-run: at ~01:36 a
worker outside this session began re-decoding live51's Spanish lane
from fragment 0, so the banked spa wav is no longer mid-lane-anchored;
the anchored-audio case now SYNTHESIZES its mid-lane fixture instead of
leaning on that file (mapping under test unchanged).

## ORIGINAL GATE RESULTS (8/08 ~00:53)

data/live51 (26.6 min, real churn) and data/catchup (15.2 min):
  tool exit 0; 0 failed chunks; layout exactly hevc + mp2(eng) +
  mp2(spa) + dvb_subtitle; duration 1591.61s/912.96s vs expected
  1591.59s/912.91s (<<1%); dvbsub 2340/822 packets; video packets
  82577/42134 (>=85% of fragments x 120); TS AUDIBLE per track where the
  wav says there is sound (eng -25.2 dB @240s, spa -16.6 dB @1243s,
  catchup eng -27.3 dB).
Negative controls (all correctly FAIL the checks): truncated TS fails
duration; spa-stripped TS fails layout; absolute k*SPF audio mapping
fails the mid-lane wav (wants frame 58.1M of a 11.0M-frame wav).
Anchored-audio case: live51 spa wav anchored 621 slots into the lane
lands rms 4487 at first_seq+10 and exact silence before the anchor.

DEPLOY GATE: 79 s VLC window on live51 replay (opened, played, closed by
the tool). Frames extracted from the OUTPUT TS: cue at 7.95s and 15.9s
render the exact srt text on schedule; 30.6s (srt gap) correctly blank.

Meter over live51: video last 108s = 2658 frames, 0 error lines,
null-sink 2166 fps (41% of slot-span -- the recording's tail is genuine
fade churn); eng 100% non-silent; spa 100% non-silent (229s pair-element
wav); CC 132 cues/5 min, anchor exact, in-range 0.97.

## THE THREE FFMPEG TRAPS (read before touching the mux)

1. Text->bitmap is impossible; PGS `.sup` + `-fix_sub_duration` is the
   bridge. Without fix_sub_duration, dvbsub packets land 2^32 ms late.
2. `-itsoffset` + `-fix_sub_duration` on the sup input LOSES every event
   still queued at input EOF (one stale unoffset timestamp). Sup times
   are therefore ABSOLUTE (session clock baked in), no itsoffset.
3. Without `-copyts`, ffmpeg rebases the sup input so the per-chunk dummy
   lands at pts 0 -- pts-0 packets sprinkled through the TS break
   timestamp-bisection seeks (input -ss reads the file head; a TS whose
   whole-file loudness was -25 dB probed as -91 dB). The mux runs -copyts
   and every stream carries literal session time.

## CAVEATS / KNOWN LIMITS

* AUDIO MAPPING IS GENERATION-BOUND: wav sample 0 = the sidecar's
  first_seq fragment of the CURRENT lane generation. After a lane roll
  the mapping is broken until the worker restarts (sidecars carry no
  generation field yet). Same limitation as v1, now documented in code.
* A cue spanning a chunk boundary blinks ~250 ms at the seam (the
  transparent dummy owns the first 0.2 s of each chunk). Cosmetic;
  roll-up cues are mostly short.
* ffmpeg's dvbsub decoder prints a handful of "Error decoding subtitles"
  lines when decoding the concatenated TS (chunk-boundary page resets);
  every extracted frame renders correctly, VLC uses its own decoder.
* VLC stdin is DEAD on Windows ("cannot peek", both `-` and fd://0) --
  the growing-file player is not a fallback, it is the design.
  CORRECTED by E37: the E35 claim that VLC "drains to EOF, waits,
  resumes" is FALSE for a real starvation -- at EOF it REWINDS to zero
  and replays (measured, VLC 3.0.23). The lead governor + pause
  supervision exists precisely so the player never touches EOF.
* live51's eng wav genuinely opens with ~3 min of silent air (recorded
  through a fade); meters/gates listen where the wav has sound.
* The spa lane needs the pair-element worker
  (`--pid 14 --element pair`); a 5_X worker on pid14 produces a
  perfectly-formed all-silence wav (0/42502 frames closed).
