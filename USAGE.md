# Using the receiver

Everything the receiver does is reachable from one command with different
options. If a command misbehaves, run the doctor first — it knows the common
failures by name:

```
python tools/atsc3_doctor.py
```

> **`--rf N` is required when tuning a radio.** There is no default channel; a
> shipped default would be wrong for everyone but the person who shipped it.
> Find your local NextGen TV (ATSC 3.0) channel with a channel scan or a
> coverage map first. Replaying a file with `--capture` needs no `--rf`.

---

## 1. Watch live television

```
python -m atsc3 watch --rf N
```

Tune channel *N*, detect the bootstrap, decode the physical layer, reassemble
the media, decode the audio, and open a player — running until you press
`Ctrl-C`. On a 6-core-or-better CPU this keeps up with the air in real time; no
GPU is needed.

Common options:

| option | what it does |
|---|---|
| `--secs N` | stop automatically after *N* seconds |
| `--record OUT.mp4` | save what was played to a file |
| `--player none` | decode headless — no window (pairs well with `--record`) |
| `--player-args "…"` | pass extra arguments through to the player |
| `--rfgain N` / `--ifgr N` | set SDR RF gain / IF gain reduction (per-site; tune for your antenna) |
| `--ant NAME` | choose the SDR antenna port |
| `--json OUT.json` | write the full run telemetry as JSON |
| `--report N` | progress-line interval in seconds (default 5) |
| `--live-dir DIR` | where the media lanes are written (video/audio/subtitle tracks) |

Example — record 5 minutes headless, keep the telemetry:

```
python -m atsc3 watch --rf N --secs 300 --player none \
    --record tonight.mp4 --json tonight.json
```

---

## 2. Choose the audio language

NextGen TV multiplexes can carry more than one audio program — commonly a
primary language and a second (e.g. Spanish). The receiver decodes both with its
own AC-4 decoder and offers them as **selectable audio tracks in the player** —
switch tracks the same way you would for any multi-audio file. Both are muxed
into `--record` output, so a recording keeps every language.

---

## 3. Turn captions on and off

Closed captions (IMSC1 / TTML) are decoded and muxed as a **real, toggleable
subtitle track**, time-aligned to the picture — not burned into it. Toggle them
with your player's subtitle control; in a recording they remain a separate track
you can enable or disable later.

---

## 4. Read the on-air program guide

The broadcast Electronic Service Guide (ESG) is decoded straight out of its
over-the-air container — service names, and the schedule of what is on and what
is coming. Note that the guide is delivered on a separate flow from the live
media, so a guide captured during a live session can be **stale**; the UI says
so with a banner rather than pretending it is fresh.

---

## 5. Inspect a locked (encrypted) service

Some services are DRM-encrypted. The receiver **never attempts to decrypt or
circumvent them** — a hard rule enforced in the code. What it *does* do, where a
locked service leaves them in the clear, is surface the service's **guide,
manifest, track list, and live captions**, with the media itself clearly
labeled `ENCRYPTED — guide and captions only` and zero "play" controls. It is an
inspector, not a key.

---

## 6. Pick the transport: MMTP or ROUTE

ATSC 3.0 delivers media over one of two transports, and this receiver decodes
both:

```
python -m atsc3 watch --rf N            # MMTP (default when present)
python -m atsc3 watch --rf N --route    # ROUTE/DASH service
```

Some ROUTE services on a multiplex are ones a commercial NextGen tuner will not
present at all. Use `--route` to select the ROUTE path when a multiplex carries
one.

---

## 7. Decode the LDM (enhanced) layer

Layered-Division Multiplexing stacks two signals on one channel — a robust core
layer and a higher-capacity enhanced layer. The enhanced layer is heavier to
decode:

```
python -m atsc3 watch --rf N --ldm
```

This path is **not real-time yet** on the reference hardware — it decodes
correctly but below 1× — so treat it as a recording/analysis mode for now.

---

## 8. Replay a captured signal (no radio required)

Any raw IQ capture can be decoded offline — useful for testing, for sharing a
signal, or for working without the SDR connected:

```
python -m atsc3 watch --capture signal.iq --rate 6.144e6 --player none --record out.mp4
```

`--rate` is the capture's sample rate. The replay path is byte-for-byte
identical to a live decode of the same air.

---

## 9. Prove it: the gate

The gate re-decodes the same signal through an independent offline path and
checks the streaming receiver against it — the receiver's own honesty test:

The gate runs on a **capture**, not on live air — both paths have to decode the
same samples for the comparison to mean anything. Record one first:

```
python tools/atsc3_capture.py --rf N --secs 12 --out signal.cs16
python -m atsc3 gate --capture signal.cs16 --rate 6912000
```

`atsc3_capture.py` refuses to bank a clipped or short capture (`VERDICT: VOID`)
— if it says so, add attenuation with `--rfgain`/`--ifgr` and record again.

The leg that matters is **`end_to_end`**: it passes only when the streaming
receiver's IP datagrams are **byte-identical** (matching SHA-256) to the
reference batch decode of the same samples. That is what separates a picture
that is *correct* from one that merely *looks* correct.

> **Known issue — `resampler_8` / `resampler_10` may FAIL.** These legs assert
> that the blocked resampler is *bit*-identical to `scipy.signal.resample_poly`
> over a synthetic signal, and some SciPy builds round the last bit differently
> (observed: ~1.2e-06 max delta on SciPy 1.15.2, about 10 float32 ULP). It does
> not affect decoding — `end_to_end` still reports identical datagrams. The
> check is deliberately left strict rather than loosened to hide the
> difference.

---

## Advanced / performance tuning

These affect throughput and resource use, not what gets decoded. Defaults are
sensible; reach for these only when you are chasing real-time on a specific box.

| option | what it does |
|---|---|
| `--decode-procs N` | number of decode **processes** (decoding is process-parallel; the GIL caps threads) |
| `--fe-threads N` | front-end worker threads |
| `--accel cpu\|gpu` | select the decode engine — **GPU is optional and often slower on big CPUs** |
| `--iters N` | LDPC belief-propagation iteration cap |
| `--no-repair` | disable optional error-repair stages (diagnostic) |
| `--realtime` | pace decode to wall-clock instead of running flat out |
| `--prebuffer N` | seconds of media to buffer before the player starts |
| `--cpu-isolate` / `--exact-cpu` | pin work to specific cores |

Run `python -m atsc3 watch --help` for the complete list.

---

## When something goes wrong

1. `python tools/atsc3_doctor.py` — it names the failure (missing spec tables,
   no SDR driver, player missing, and so on).
2. No bootstrap detected? It is almost always the **antenna** — NextGen TV is
   usually on UHF and rarely the strongest signal in the band. Try a real UHF
   antenna with gain before suspecting the software.
3. Decodes but cannot watch? You have the decoder but no player/ffmpeg — see
   [INSTALL.md](INSTALL.md) step 2.
4. Live stutters but a `--capture` replay is clean? That is a real-time /
   throughput limit on your CPU, not a decode error — see the tuning table above.

For the science of *why* each of these stages exists, see the
[how-it-works explainer](https://felbs.software/atsc3). <!-- scrub-allow: HOST_USER public project domain, not the box login -->
