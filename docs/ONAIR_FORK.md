# The OnAir fork: our receiver's ideas, offered to their shell

Started 2026-11-08.  Goal, in Felbs' words: a version of OnAir with our code, working, shown to the
OnAir team to see if they want to adopt it.  This file is the plan and the running record; the
measurements go to PRIOR_ART.md as they are made.

## The rule that shapes everything

OnAir is GPL-3.0-or-later.  This repository is Apache-2.0.  Apache code may be carried INTO a
GPL work; GPL code may never come back.  So:

* the fork is a GPL work and lives in its own checkout (`../onair-upstream/OnAir`, later Felbs'
  GitHub fork); anything written there is GPL and stays there;
* nothing from OnAir is copied here -- ideas only, credited in NOTICE (that rule predates this file);
* what we carry into the fork is written fresh in their language (C++20) from OUR design and the
  standards, not transliterated from their files.  Our Python is the reference the fork's C++ is
  gated against (same capture, same Blocks, same bytes).

## What they have, measured on 2026-11-08 (upstream 25034214, 2026-11-04)

| | OnAir | ours |
|---|---|---|
| shape | C++20 library (`core/`) + ImGui app (`app/`), Rust protocol layers (required to build), one binary, 200+ modes | Python + C kernels, one command (`tools/atsc3_watch_av.py`), ATSC 3.0 only |
| ATSC 3.0 chain | 8.7k lines: bootstrap with CFO, Preamble, L1, every Subframe, CTI, LDPC 16200/64800 (CPU, GPU on macOS/Windows), ALP, IP, ROUTE -> their MPEG-TS remux -> FFmpeg player | bootstrap -> L1 -> every Subframe -> CTI/HTI -> layered LDPC -> ALP -> IP -> ROUTE and MMTP -> fMP4 lanes -> our mux/viewer |
| MMTP | **not supported** (`atsc3_receiver.cpp`: ROUTE only; UI says "not supported yet") | live television since M11 (RF33 WJLA, RF25 Fox) |
| LDM | the L1 field `ldmInjectionLevel` is parsed; no cancellation, no Enhanced layer | Core layer live; canceller measured (E124), Enhanced on synthesized signal (E122) |
| AC-4 | through libavcodec only when that FFmpeg has it | our own decoder, 5.1 and the ASPX_ACPL_3 stereo |
| real air | tested on their generator + impairments (`test_atsc3_realworld.cpp`); README asks for recordings | every number in PRIOR_ART.md is off air |
| radios | HackRF native, SDRplay API, SoapySDR; the engine watches unplug/reconnect | RSPdx through SoapySDR, single-tenant lock, the supervisor |
| gain | AUTOGAIN.md: a guarded servo on MER/SNR | gr-rxtune: decoder-truth dial, one-shot search (the servo ideas are the next job, see memory) |
| head to head (13.2, RF33, same PC) | 18.0 s for 20 s of air, PLP 1 9046/9360 with two bursts, WHUT 32.1 plays without audio | 10.8 s, 9360/9360, WJLA 7.1 plays with both audio tracks and captions |

Their CLI `tools/atsc3rx` takes `--format cs8|cu8|cf32` (not cs16); our banks are converted with
numpy into `C:\captures\onair\*.cf32`.

## What they would adopt, in order

1. **Real-air recordings and results.**  Their README asks for them.  Which recordings can be
   shared is a consent question (a capture identifies the market); decided by Felbs before anything
   leaves this box.  The 13.2 race numbers go with them.
2. **MMTP.**  The biggest gap for a US user: here, WJLA 7.1 and RF25's Core-layer service
   (WNUV 54.1, by its SLT) are MMTP.  Their receiver
   already ends in an HLS/MPD pipeline for ROUTE (`Atsc3Receiver::mpd()`), so MMTP's MPUs can feed
   the same player as fMP4.  Our m11 design (MMTP packet -> MPU reassembly -> init + media segments,
   the signalling in MP table / MPI) is the reference; the C++ is new.
3. **AC-4** -- their player is mute on most US stations without it.  Our decoder is 10k lines of
   Python with a C QMF kernel; a port is a project of its own, after MMTP.
4. **LDM.**  The Core layer already decodes in their chain as a plain PLP (the Enhanced layer is
   noise at -5 dB to it) -- to be measured on the RF25 bank below.  The canceller and the Enhanced
   path (E82-E126, docs/ARCHITECTURE.md s3b) are what we would add.
5. **The supervisor's lessons**: fade is not a wedge (E63), the LDM status line (13.27), the
   single-tenant radio, the kill-on-close tree.  Their engine has its own watchdog; these are
   review notes, not code.

## Steps

- [x] clone upstream, read ARCHITECTURE.md, the receiver API, the tests, the build notes
- [x] install the Windows toolchain they document (MSYS2 UCRT64: gcc, cmake, ninja, pkgconf,
      ffmpeg, glfw, hackrf, soapysdr, libiconv, rust) -- `C:\captures\msys2_install.log`
- [x] build `atsc3rx` and the app on this box (12 minutes, 1391 targets, 0 errors;
      `C:\captures\onair\build.log`)
- [x] run `atsc3rx` on our banks (below)
- [ ] Felbs forks the repository on GitHub; the work moves to a branch there
- [x] why their receiver sees 0 service list tables on RF25: every LLS table there is a
      SignedMultiTable (0xFE) and `parseLls()` knew tables 1-5 -- fixed on branch `felbs/real-air`
      (699015b6, unit test, 25/25 ATSC 3.0 tests pass): 12 services found.  Still no transport
      stream, because the ROUTE services' signalling addresses are not in the Core-layer
      datagrams: they ride the Enhanced layer; the Core carries one MMTP service (WNUV 54.1)
- [ ] MMTP in C++, gated against our lanes byte for byte on the same capture
- [ ] the write-up for the OnAir team: what was measured, what was added, what is offered

## Results

**2026-11-08, upstream 25034214 built here, same two 20 s banks, one after the other
(PRIOR_ART.md s13.28).**  Their `atsc3rx` at 12 threads; ours at 8 processes.

| bank | | OnAir | ours |
|---|---|---|---|
| RF33, 20 s, Antenna B (`rf33_win_20s`) | wall | **7.1 s (2.8x)** | 13.4 s (1.47x) |
| | PLP 0 / PLP 1 / PLP 16 | 5920 / **9360** / 80, 15,360 baseband packets, 0 bad, 0 frames failed | 5920 / 9333 of 9360 / -- |
| | services | 10 found; 9 ROUTE to a 7.07 MB TS; WJLA 7.1 "MMTP (not supported)" | WJLA 7.1 MMTP lanes, two audio languages, captions |
| RF25, 20 s, Antenna A (`rf25_s3_antA`) | wall | 23.7 s (0.8x; "frame cut and resampling" 15.5 s of it) | **9.1 s (2.19x)**, resampler bypassed |
| | Core layer | 3,130 baseband packets, 0 bad, 0 frames failed -- the plain-PLP path carries the Core out | 2,835 / 2,835 FEC Blocks, CTI phase locked on Frame 0 |
| | services | **0 service list tables, 0 bytes of TS, exit 2** | 37 MPUs, lanes written |

What changed since the 13.2 race (2026-11-06): their PLP 1 bursts are gone (9046 -> 9360 of
9360) and they are now 1.9x faster than us on RF33.  What has not changed: MMTP, so WJLA does
not play; no AC-4.  New: on RF25 their physical layer decodes the LDM Core layer as a plain PLP
(as predicted -- the Enhanced layer is noise at -5 dB to it) but the signalling never resolves,
so nothing plays.  Ours is 2.6x faster on that carrier because it does not resample 6.912 Msps.
Ours left 27 PLP 1 Blocks on the table on RF33 that they did not -- an item for our list.

**First fork commit, 2026-11-08 14:13** (`felbs/real-air` 699015b6): SignedMultiTable in
`parseLls()`.  Their receiver on the RF25 bank before / after: 0 / 20 service list tables,
0 / 12 services; transport stream still 0 bytes -- the ROUTE services live on the Enhanced
layer (their SLS addresses 239.255.N.1:8000 are absent from the 10,040 Core datagrams), and the
Core's one service is MMTP.  On this carrier OnAir plays nothing until MMTP exists.

**The referee at the cliff, 2026-11-11 (PRIOR_ART.md s13.34).**  The same `atsc3rx` build on the
45 s RF25 bank where the Core layer sits at its threshold, and on its noisy copies:

| bank | OnAir, baseband packets good / bad | of which Core, about | ours, Core Blocks converged / judged |
|---|---|---|---|
| `rf25_win_45s_b` | 802 / 6342 | 95 | **4887 / 5986** |
| `rf25_win_n20` | 708 / 6436 | 1 | **4570 / 5986** |
| `rf25_win_n15` | 707 / 6437 | 0 | **3701 / 5920** |
| `rf25_win_n12` | 707 / 6437 | 0 | **1551 / 5250** |

On the clean bank above their plain-PLP path carries the whole Core out.  Here it carries about
1.5 % of it and ours about three quarters: the box-plus decoder and the ladder (s13.20) are what
we have to offer on this carrier, besides LDM itself and MMTP.  Their run is 0.3x real time on
these files (resampling), ours 1.6-2.4x.  One thing to hand over with the write-up: their
IQ-file probe finds a WAV header's offset and the playback start does not use it, so the header
is played as samples.

## The product route, 2026-11-09: our receivers behind their window, no C++ first

Felbs, 2026-11-09: "our version of OnAir, ATSC 1 and 3, so we can show their devs how to do live
TV properly ... it will hopefully work with their UI when we're finished."  Two routes, in order.

**Route 1 (done the same day): OnAir as the player, through its IPTV mode.**  OnAir's Internet TV
mode (docs/IPTV.md) plays M3U playlists and MPEG-TS over HTTP through the same player as its
antenna modes, with an XMLTV guide; fragmented MP4 and DASH it does not take yet.  tvTuna
(Z:/src/tvTuna, GPL, a separate repository) is the glue: one command per job over the Software TV
Tuner's 8-VSB chain, this receiver, gr-rxtune's search and servo, and the site radio warden; its
`serve` command is a network tuner (playlist, `/ch/<rf>.ts` at the live edge, `/guide.xml` from the
decoded PSIP).  ATSC 1.0 streams the chain's own transport stream; ATSC 3.0 streams this
receiver's muxer output (`atsc3_tv.py --player none --out`: HEVC copied, AC-3 5.1 and MP2 with
their languages, captions as a soft track).  Measured 2026-11-09 with ffmpeg as the client: ATSC
1.0 flowing 2 s after the first request, ATSC 3.0 after 38 s, 60 s read at 59.94 fps with FEC
100 %.  OnAir's app on this box, pointed at the playlist, is the demo: their shell, our receivers.

**Route 2 (next): a native engine mode in their core.**  What the IPTV route cannot show is our
numbers in their panels: the gain servo's state, decoder truth, frames and interval FEC, the
fade-versus-wedge call, the warden.  The shape, from ARCHITECTURE.md: a receiver in `core/src/<mode>/`
is fed by the Engine and emits a transport stream to `ts.cpp` / the Player; the app draws what the
library reports.  The plan is an **external-engine receiver**: a mode whose "receiver" spawns and
supervises an external process (tvTuna's engine over a local control socket: tune, gains, status
JSON) and hands its MPEG-TS to the Player, with the status fields mapped onto their quality and
service panels.  New C++ in the fork, GPL, written from this design; no OnAir code comes here.
Steps: (a) a `tvtuna engine` command with a line-protocol control channel (JSON per line: tune /
status / stop) and the TS on a second pipe or a local TCP port; (b) the fork's `external_rx`
receiver + mode-table entry + a quality adapter; (c) the panels: held setting, MER or SNR with
cliff, pictures or FEC, warden holder; (d) gate: the same air through both routes gives the
same stream bytes.

Open on the way: OnAir's IPTV player and CMAF (if they add fMP4 HLS, the ATSC 3.0 lanes can be
served as they are, no remux); their player's handling of our TS when the chain re-acquires
(discontinuity flags); the guide for ATSC 3.0 (SLT names, ESG) which neither side parses yet.

**Route 2 built the same afternoon (fork branch `felbs/real-air`, 329b7c63):** an "External engine"
mode (family 241) in their core and app -- `core/src/extengine` spawns the engine with pipes, reads
its JSON status lines, pulls the transport stream from its TCP socket through their `iptv::TsAligner`
into their Player; `app/extengine/extengine_ui.cpp` is the screen (the Internet TV screen's shape):
servo state, held setting, MER against the cliff, pictures per second or frames and FEC, packet and
continuity counts, who holds the radio, then their player and controls.  Two tests in their tree
(status lines recorded from tvTuna; numbers / row / screen / feed branch), both pass with
`mode_numbers` and `iptv`.  Seen in their window on real air (screenshot kept off the repositories):
"External engine up MER 18.6 dB", 1080p30 MPEG-2 + AC-3 5.1, 2478 pictures, 0 late, 0 errors, the
engine holding the site warden and exiting with the app.  Engine-side lesson: a receiver spawned
from another program must put its own interpreter's folders first on PATH, or it loads the host's
SoapySDR and finds no radio.  The fork's `docs/EXTERNAL-ENGINE.md` carries the protocol.
What remains for the pitch: the write-up itself and Felbs' pushes (tvTuna has no remote yet).

**Later the same day: the picture was a postage stamp, and three faults behind it.**  Felbs could
not see the picture in the first screenshot.  He was right: the mode's tab gave it about 95 x 55
pixels of a 1500 x 900 window, under the settings and the measurements.  Fork 14fcbe35 moves the
measurements to the row under the tabs (where the TV modes draw constellations) and folds the
engine's settings away while a channel plays: about 325 x 183 pixels, the full height of the tab.
It also plays the multiplex's lowest programme number instead of whichever PMT arrived first (one
run had shown the 1080 main, the next a 480 subchannel).  Timing the stop on air then found:

- stop took exactly 90 s: the engine joined its worker while holding the lock the worker needed to
  finish (tvTuna 5802261).  OnAir kills an engine after 90 s, so a chain could have been orphaned
  holding the radio.  With the servo's rest now honouring `until` (gr-rxtune 917d44e): 15.5 s when
  the quit lands in the search, 0.3 s in a hold, radio released each time.
- a proven hold equal to the default start setting read as "no hold", so a known channel was
  searched from scratch in front of the viewer, dead cells included (same commit).
- still open: the confirmation search is seven cells and about 180 s, visible as "searching" while
  the picture plays; the stream's first seconds carry a continuity burst that varied 61 to 2737
  errors across five starts, then 10 to 20 errors every 5 to 13 s (the standing STVT question);
  OnAir's player counted 44 decode errors in one 90 s run.  The on-screen continuity count is
  cumulative from the chain's start, so it reads high beside a clean picture.
