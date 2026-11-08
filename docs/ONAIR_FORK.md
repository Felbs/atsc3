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
