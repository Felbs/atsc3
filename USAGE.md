# Using the receiver

Everything the receiver does is reachable from one command with different
options. If a command misbehaves, run the doctor first — it knows the common
failures by name:

```
python tools/atsc3_doctor.py
```

> **`--rf N` is required.** There is no default channel; a shipped default would
> be wrong for everyone but the person who shipped it. Find your local NextGen TV
> (ATSC 3.0) channel with a channel scan or a coverage map first.

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
| `--report OUT.json` *(alias `--json`)* | write the full run telemetry as JSON |

Example — record 5 minutes headless, keep the telemetry:

```
python -m atsc3 watch --rf N --secs 300 --player none \
    --record tonight.mp4 --report tonight.json
```

---

## 2. Choose the audio language

NextGen TV multiplexes can carry more than one audio program — commonly a
primary language and a second (e.g. Spanish). The receiver decodes both with its
own AC-4 decoder and offers them as **selectable audio tracks in the player** —
switch tracks the same way you would for any multi-audio file. Both are muxed
into `--record` output, so a recording keeps every language.

**Which AC-4 decoder.** The audio worker uses `libac4dec`, the C library of the
separate `ac4-decoder` project, when it can load it (full SAP, the A-CPL 5.1
upmix of ASPX_ACPL_3 programmes, the A-SPX fixes, damaged-frame concealment),
and falls back to the receiver's own decoder otherwise. It looks for
`$ATSC3_AC4DEC_LIB`, then `../ac4-decoder/c/libac4dec.so` next to this
repository (build it with `make -C ../ac4-decoder/c libac4dec.so`). The first
lines of the audio log say which one is running. `ATSC3_AC4_BACKEND=python`
forces the old decoder. Stereo from a 5.1 programme is a proper downmix with
the centre (dialogue) and surrounds folded in.

**Level control (ALC).** Services differ a lot in level (one broadcaster's
stereo decodes ~17 dB under its own loudness metadata), so changing service
used to jump. The audio worker can level what it decodes: it measures the
loudness of the decoded sound as ITU-R BS.1770 / EBU R128 do (K-weighted,
gated: silence and speech pauses do not count), steers one gain for all
channels toward a target, and limits peaks so a boost never clips. It
settles in 1-3 s after a service starts, then rides slowly (a step in the
programme's level is followed at up to 2 dB/s down, 1 dB/s up). The limiter
has no lookahead, so nothing is delayed. `lab/m46_alc.py` has the details;
`python lab/gate_alc.py` checks it.

The live launchers (`tools/atsc3_watch_av.py`, `tools/atsc3_watch_omarchy.sh`,
`tools/atsc3_tv_linux.sh`) turn it **on**; a worker started by hand leaves it
off unless asked (its output is then byte-identical to before), so gates and
recordings rendered offline keep the decoder's own level.

| setting | what it does |
|---|---|
| `ATSC3_ALC=1` / `0` (worker: `--alc on\|off`) | level control on / off |
| `ATSC3_ALC_TARGET=-24` (worker: `--alc-target`) | target loudness, LUFS (ATSC A/85 says -24) |
| `ATSC3_ALC_MAX_BOOST=20`, `ATSC3_ALC_MAX_CUT=20` | the most it turns up / down, dB |
| `ATSC3_ALC_CEILING=-1` | limiter ceiling, dBFS |
| `ATSC3_ALC_TAU=15` | how many seconds of programme the slow level estimate averages |
| `ATSC3_ALC_TRACE=FILE` | one line per 100 ms: time, block loudness, level, gain, limiter |

The audio log's progress lines end with the gain in use, e.g.
`alc +17.9 dB (in -41.9 LUFS)`.

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

```
python -m atsc3 gate --rf N
```

A pass means the live streaming path produced **byte-identical** output (matching
SHA-256) to the reference batch decode. This is the check that separates a
picture that is *correct* from one that merely *looks* correct.

---

## Advanced / performance tuning

These affect throughput and resource use, not what gets decoded. Defaults are
sensible; reach for these only when you are chasing real-time on a specific box.

**Small boxes (< 8 CPUs): the chain drops its margin levers by itself on strong carriers.**
The E60 margin levers (CE smoothing, weighted LLRs, exact-BP rescue) buy decode
margin at the cliff but cost ~12-20 % of the chain on healthy air (measured N=3 on
a 2-core laptop: 1.28x -> 1.44x offline, FEC identical). `ATSC3_MARGIN_AUTO_OFF=<dB>`
makes the chain skip the lever machinery once a frame's dummy-cell SNR reads at or
above that value, and bring it back 2 dB below it (hysteresis); the chain log says
`margin levers OFF/ON` when it switches. The couch launcher (`atsc3_watch_av.py`)
sets 19 dB on boxes with fewer than 8 CPUs; `0` (the default elsewhere) keeps every
existing result byte-identical. On that laptop, 10-minute live runs went from a
drifting lag (22 -> 60 s, 5368 audio pads) to a steady one (22-30 s, ~2800 pads).

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
5. The supervisor (`tools/atsc3_run.py`) restarts a chain that dies. When the
   chain died with a Python traceback it first writes `crash_<time>.txt` in the
   live directory (the traceback and the 40 log lines before it) and raises a
   desktop notification; `ATSC3_CRASH_STOP=1` makes it stop instead of
   restarting, so the bug can be looked at.

For the science of *why* each of these stages exists, see the
[how-it-works explainer](https://felbs.software/atsc3). <!-- scrub-allow: HOST_USER public project domain, not the box login -->
