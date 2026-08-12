# atsc3 — an open-source ATSC 3.0 (NextGen TV) receiver

A complete software receiver for ATSC 3.0, built from the published standards:
**bootstrap → L1 → LDPC/BCH → ALP → IP → MMTP/ROUTE → a picture on the screen.**

It runs **in real time on an ordinary CPU**. No GPU, no FPGA, no demodulator
chip. The reference live run is a 6-core desktop processor from 2016.

There was no open-source ATSC 3.0 receiver that got to a picture, and there was
no open-source AC-4 decoder at all — ffmpeg demuxes AC-4 and cannot decode it,
VLC cannot parse the sample entry. Both were written here, from the specs.

**Project site:** [felbs.software](https://felbs.software) · **Contact:** [E@felbs.software](mailto:E@felbs.software) <!-- scrub-allow: HOST_USER public project domain and contact address, not the box login -->

---

## What it does

* **Live television.** Tune a NextGen TV channel, demodulate it, decode the
  physical layer, reassemble MMTP media, decode the audio, and play it — with
  a viewer that survives fades, respawns, and six-hour nights.
* **Two audio languages**, decoded by our own AC-4 decoder, switchable in the
  player.
* **Closed captions**, muxed as a real subtitle track and toggleable.
* **Both transports.** MMTP *and* ROUTE/DASH. One of the ROUTE services on the
  local multiplex is one a commercial NextGen tuner on the same antenna
  **cannot present at all** — it reports `(no data)` and offers no programme.
  We have its picture.
* **LDM.** A layered-division-multiplex core layer, decoded from a second
  transmitter, with sound.
* **The service guide**, decoded out of the broadcast ESG container.
* **DRM-locked services are enumerated, labeled, and left alone.** See
  [Scope limit](#scope-limit).

## What it does not do

Stated first, and in full, because a receiver that oversells itself wastes the
reader's evening.

* **Not every service is live.** One clear ROUTE service decodes **offline
  only** — its physical-layer subframe is not yet real-time, and making it so
  is a campaign, not a flag.
* **The LDM path is not live yet.** It decodes perfectly and runs at **0.48×**
  real time: about 2× short in one process. That subframe carries **5.9× the
  cells per frame** as the live one — 1,216,939 against 207,900 — in the same
  242 ms. Work in progress.
* **Cold tune-in is about 9× slower than commercial silicon**: 13.0 s against
  1.41 s. See the [scoreboard](#the-head-to-head-scoreboard).
* **Live 5.1 is not available.** 5.1 is fully decoded and verified as discrete
  broadcast surround — but only on **recordings**. Six-channel synthesis runs
  0.62–0.91× on live air, so the live setting is stereo (1.59×). The bottleneck
  is a Python-loop-bound filterbank, not the maths.
* **The full 5.1 upmix of the second broadcaster's audio is not implemented.**
  What we render there is the coded stereo core — which *is* the broadcast
  downmix, not something we invented. See [the AC-4 decoder](#the-ac-4-decoder).
* **The guide goes stale.** The ESG arrives on a ROUTE flow; the live chain
  decodes the MMTP service. There is no live path to a fresh guide, and the UI
  says so with a banner rather than pretending.
* **Some services are simply out of reach.** On one nearby multiplex both TV
  services ride the LDM *enhanced* layer and measure **9.3 dB short**. The
  levers we have are worth about 1 dB against that.
* **One tuner.** The commercial box has four (two ATSC 3.0). We have one SDR,
  so we cannot record one service while watching another.
* **We do not report per-PLP modulation from our own L1** the way a commercial
  tuner does; we attribute a flow to a PLP by remembering which baseband file
  it came out of.
* **We ship the numeric tables, not the standards.** The physical-layer
  constellations and the AC-4 band layouts are shipped as small numeric
  artifacts (`lab/spec_bank/`, ~23 KB), so a fresh clone decodes picture *and*
  sound with no extra downloads. The standards **documents** are copyrighted by
  their publishers and are not redistributed here. If you do have them, the
  decoder parses them and prefers that path, because two independent editions
  cross-checking each other is stronger evidence than one file we wrote. See
  [INSTALL.md](INSTALL.md).

---

## Measured performance

All figures are medians of N ≥ 3 unless noted. "Real-time factor" is decode
throughput over wall clock; 1.0× means it keeps up with the air.

| machine | engine | sustained | notes |
|---|---|---|---|
| Threadripper 2990WX (32C/64T) | CPU only | **2.66×** | 4 decode processes, ~4% of the box |
| Threadripper 2990WX | GPU (optional) | 1.48× | the GPU path is *slower* here |
| Ryzen 5 1600X (6C/12T, 2016) | CPU only | **1.02× live** | full TV including the viewer |
| GTX 1060 (6 GB, 2016) | GPU | 0.896× | on the 1600X box |

**The GPU is an optional accelerator and never a requirement.** On a large CPU
it loses to the CPU path outright.

How the CPU path got there, since the honest starting line matters: **0.43×** →
1.41× (decode in **processes**, not threads — the GIL caps decode threads at
~1.2× on any core count) → 1.66× (front-end fast paths: fine timing 72.2 → 2.4
ms/frame, notch 72.0 → 0.0, derotate 50.0 → 3.3) → **2.66×** with an optional
compiled C min-sum LDPC kernel (LDPC stage ~131 → ~20 ms/frame). After that the
wall is no longer LDPC — it is the demapper's sgemm.

### The live run

```
263 s wall · 1048 frames · 77,547 of 77,552 FEC blocks converged · 0 underruns
254 s of 720p60 HEVC that re-decodes with zero error lines
```

The decoded IP datagrams are **byte-identical (same SHA-256)** to the batch
chain re-run on the same air.

### Byte-identical across machines

32 frames, decoded on Windows and Linux, AMD desktop and AMD workstation, CPU
and GPU: **2368/2368 FEC converged, 2368 BCH zero, 0 resyncs, 2768 IP
datagrams, 13 MPUs** — the same numbers on every machine and mode combination.

### The stability soak

```
6.02 h continuous · 0 invariant violations · 0 interventions · 0 maintenance windows
A/V correlation referee   6 runs, 6 pass, 0.980–0.994
video/audio skew         25 checks, worst +1.73 s against a 2.0 s limit
independent RF referee   35 readings, 0 disagreements
warden actions            2, both designed log-file rolls
```

Two caveats the run itself recorded:

* The clock is **qualified time**, not wall time — it accrues only while the air
  is actually delivering. Six qualified hours can take longer than six hours.
* The invariants were **exercised, not merely armed**. That distinction is the
  reason the soak exists; a monitor that never fired is not evidence.

And the honest sequel: the same run **did** violate at 7.72 h — a false positive,
the fifth of its family. The chain was healthy; the *judge* compared "FEC is
good now" against "the file grew over the last 400 s" and the misaligned windows
convicted a stack that was behaving perfectly. Fixed by only counting time the
output was actually owed.

---

## The head-to-head scoreboard

Against a commercial NextGen TV tuner, **on the same antenna through a passive
splitter**, so every paired reading is same-air and same-instant.

| dimension | verdict | measured |
|---|---|---|
| coverage | **TIE** | 1 win, 1 loss, 23 ties |
| concurrent delivery | **LOSE** | 83.77% vs 98.71% of paired samples |
| sensitivity | **NOT COMPARABLE** | no common attenuation ramp exists |
| acquisition | **LOSE** | 13.0 s vs 1.41 s cold |
| stability | **LOSE** | 105 min vs 5 min longest outage |
| capability | **WIN** | 6 capability wins to 4 |

**Win 1 · lose 3 · tie 1 · not comparable 1. We are behind on more dimensions
than we are ahead.**

The unflattering numbers, unpacked rather than buried:

* **Acquisition.** Of our 13.0 s, ~5.2 s is decode-pool warm-up before any RF
  work happens at all; the bootstrap lands at 7.0 s (a very tight 7–8 s p10–p90)
  and the first media unit at 13.0 s. Pre-warming the pool should recover ~5 s
  and a warm re-tune should be 2–3 s. **It will not reach 1.4 s** — their
  demodulator is a dedicated ASIC that is already running when the request
  arrives. *Caveat: n=2 cold starts on the current configuration, and their n=5
  from one sitting. That is a good sign, not yet a measurement.*
* **Concurrent delivery.** The headline quotes the whole 12.7 h window. Since
  the gain fix that closed a self-inflicted fault it reads **97.78% against
  100.00%**. Deleting the minutes we spent restarting would move the headline
  from 83.77% to 95.56% — **+11.79 points for free**, which is exactly why the
  rule is that a missing sample of our side is our downtime, never a skipped row.
* **Stability.** Same story: 105 min is the full-window number and it is the
  worse one. On the current configuration our longest outage is 5 min against
  their 0 min, which is a tie within tolerance.
* **Sensitivity was refused, not lost.** It is the number everyone wants. The
  commercial tuner exposes no gain or attenuation control, and a physical
  attenuator would touch the shared feed. There is no dataset that can talk the
  comparison into a verdict, so the code has **no path that returns one**.
  Comparing our SNR to their proprietary quality metric, or relabelling a fade
  record as a threshold, were both considered and both rejected.

What we win on: presenting a service their silicon cannot present at all; a
locked-channel inspector that surfaces guide, manifest, track list and live
captions for services **neither of us can decrypt and neither of us tries**;
AC-4 decoded in software rather than passed through to a TV; a decoded ESG,
which their API does not expose; and a gain control — the one that fixed our own
night floor — which they do not have.

**The comparison also audits us.** Their L1 readings agree with ours 6 for 6
across two transmitters. And building the scoreboard found five bugs in the
scoreboard itself: three had been flattering them, two had been understating us.

---

## The AC-4 decoder

Written from ETSI TS 103 190, because no open-source AC-4 decoder existed. It
decodes real broadcast audio from more than one encoder vendor:

* `5_X` element — 5.1, all six channels, verified **discrete** broadcast
  surround rather than an upmix (channel correlations: front L–R +0.887,
  surround Ls–Rs +0.999, L–C −0.023).
* `channel_pair_element` — stereo.
* **A-SPX** spectral extension, complete: 64-band complex QMF (78.49 dB
  reconstruction), envelope decode, HF generation, envelope adjustment.
* `ASPX_ACPL_3` core mode, and the `payload_base` field — the two fixes that
  took a second broadcaster's audio from 0.0% and 3.5% of frames to **100.0%
  and 100.0%**.
* The 1001/960 output resampler. Without it the audio runs **4.27% fast** — and
  three independent subsystems agreed on the wrong constant, because a *shared*
  constant is invisible to cross-checking.

**The honest gap.** For the second broadcaster's 5.1 programme we render a
correct **stereo** decode — the stereo core *is* the coded downmix, not a
fold-down we computed. The **5.1 upmix** needs `acpl_data_2ch()` (framing plus
five entropy-coded parameter sets) and the QMF-domain A-CPL synthesis stage.
Scoped, not hand-waved: one bitstream parser and one synthesis stage.

**The law this decoder taught, which generalises past AC-4:** a decoder that has
only ever met one broadcaster has only ever been tested against one set of
encoder choices. Every `if this flag is set` branch we have not implemented is a
station we cannot hear yet.

> The AC-4 decoder is intended to ship as a **standalone package**. See
> [AC4_SPLIT_PLAN.md](AC4_SPLIT_PLAN.md).

---

## Scope limit

**Encrypted services are detected, enumerated and labeled as locked — never
attacked, never circumvented.** This is a hard rule, and it is enforced in the
code, not just promised in a README: the channel grid renders a locked service
with an `ENCRYPTED — guide and captions only` badge and **zero play
affordances**, and a gate fails the build if a "Watch now" control is ever
planted on one.

Where a locked service leaves its guide and its captions in the clear, we read
those, because they are in the clear. The media stays locked.

We also do not work around access controls on the standards documents. The
specs must be downloaded by hand.

---

## How this was built, and what it cost

Roughly two weeks, spec-first, with an experiment notebook. Some of what is in
that notebook is more useful than the code:

* **A bug can wear the costume of a hardware limit.** A carrier written off as
  "needs a better antenna" decoded perfectly once a **one-frame interleaver
  phase offset** was corrected. 242 milliseconds of bookkeeping stood between
  "this transmitter is dead" and "this transmitter is perfect", and every
  quality metric read *excellent* the whole time while every correctness metric
  read *dead*.
* **Header fields in this standard have no integrity protection.** One flipped
  bit killed a media lane permanently and allocated 42 GB. Every field read from
  the air is now bounds-checked — that is a class to audit, not two bugs to fix.
* **Gate the symmetry, not the part.** A stereo balance bug passed three gates
  with the left channel exemplary and the right 20× too loud, because all three
  gates measured the left.
* **A shared constant is invisible to cross-checking.** Three subsystems agreed.
  All three were wrong by 4.27%.
* **Check the instrument before you believe the negative.** Re-running a control
  on old data moved a result by 1.94 dB — the *instrument* had drifted, and
  without the control this project would have published a "+11.5 dB
  improvement" that was partly its own measurement moving underneath it.
* **A crash truncates the evidence.** A "monotone decay" turned out to be the
  acquisition transient; the steady state only existed after the crash was fixed.
* **A missing spec clause is not a weak signal.** L1 verification sat at 0/30
  because a repetition term — 48% of the codeword — was missing. The arithmetic
  had never matched the value the standard prints, and nothing was checking.

Five entries in that notebook are formal **retractions** of earlier entries. A
recorded failure is a reference, not an authority.

---

## Documentation

| file | what |
|---|---|
| [INSTALL.md](INSTALL.md) | hardware, dependencies, first tune — target: watching in under 15 minutes |
| [USAGE.md](USAGE.md) | how to use every feature — watch, record, audio languages, captions, the guide, MMTP/ROUTE, LDM, replay, the gate |
| [How it works](https://felbs.software/atsc3) | the science, from first principles — the bootstrap, OFDM, LDPC, transport and codecs, an educational walk through the radio arts | <!-- scrub-allow: HOST_USER public project domain, not the box login -->
| `docs/ARCHITECTURE.md` | the signal path and the supervision layer |
| `DESIGN_NOTE.md` | the full build ladder, every result and every named bug |
| `PRIOR_ART.md` | what else exists and how this differs |
| [AC4_SPLIT_PLAN.md](AC4_SPLIT_PLAN.md) | the standalone AC-4 decoder |
| `tools/atsc3_doctor.py` | diagnose an install before blaming the sky |

## Licence

See [LICENSE](LICENSE) and the reasoning in
[LICENSE_RATIONALE.md](LICENSE_RATIONALE.md). Short version: Apache-2.0, chosen
for its patent grant, with an explicit notice that **implementing a broadcast
standard may require patent licences this project does not and cannot grant.**

---

*Not affiliated with, endorsed by, or sponsored by ATSC, ETSI, Dolby
Laboratories, or any broadcaster or receiver manufacturer. All trademarks belong
to their owners and are used only to say what this software is compatible with.*
