# INSTALL

**Goal: a stranger clones this and is watching television in under 15 minutes.**

That is the bar. Where a step cannot meet it, this document says so instead of
pretending.

| step | time | blocker if it fails |
|---|---|---|
| 1. Python + packages | 2 min | nothing decodes |
| 2. ffmpeg + a player | 2 min | decodes, cannot watch |
| 3. *(optional)* the standards documents | 0 min — **tables now ship** | nothing; only needed to re-derive tables |
| 4. SDR driver + SoapySDR | 5 min | replay works, live does not |
| 5. optional C LDPC kernel | 2 min | ~7× slower LDPC stage |
| 6. first tune | 1 min | — |

**You do not need a GPU, and you do not need `torch`.** Nothing in that table
mentions one because the receiver decodes on an ordinary CPU; `torch` only
enables an OPTIONAL accelerator. On a 6-core desktop the tuned CPU path
measured **1.02x** real time against **0.896x** for that same box's GTX 1060 —
the GPU lost. If the chain is running slowly, the cause is far more likely to
be the antenna, the CPU governor, or the unbuilt LDPC kernel in step 5 than a
missing GPU. `python tools/atsc3_doctor.py` checks all three.

**Step 3 is no longer a blocker.** The numeric tables ship with the repo, so a
fresh clone decodes picture and sound immediately. It is kept below for anyone
who wants to re-derive the tables from the published standards themselves.

At any point:

```
python tools/atsc3_doctor.py
```

The doctor knows every failure below by name and tells you which one you have.

---

## 0. What you need

**Hardware**

* A general-coverage SDR that can deliver **≥ 6.144 Msps** of complex baseband
  over the UHF TV band. The bootstrap symbol rate is fixed by the standard, so
  that floor is not negotiable. This project's reference runs use **6.912 Msps**
  (native) and **8 Msps** (which resamples to 6.144 by an exact 96/125 ratio).
  
* A **UHF antenna** with real gain — a yagi, not a set-top loop. NextGen TV is
  usually on UHF and usually not the strongest signal in the band.
* **A CPU with 6 or more cores** if you want live TV. Fewer will decode
  recordings fine, just not in real time.
* A **short, direct USB 3.0** connection. Not a hub, not a long cable, and route
  it away from the antenna coax — a cable-tie bundle of USB and coax is an
  antenna, and USB 3.0 is a notorious UHF noise radiator.

**Software**

* Python 3.10+ (3.12 and 3.14 are both exercised)
* numpy, scipy — required
* psutil — strongly recommended (process supervision)
* pymupdf — required for **audio** only
* SoapySDR + a driver for your SDR — required for **live** only
* ffmpeg
* VLC (mpv and ffplay work as fallbacks, with less of the UI)
* a C compiler — optional, for the fast LDPC kernel

---

## 1. Python

```
python -m venv .venv
. .venv/bin/activate            # Windows: .venv\Scripts\activate
pip install numpy scipy psutil pymupdf
python tools/atsc3_doctor.py
```

> **If SoapySDR lives in a different interpreter than everything else** — which
> is normal, because SDR stacks often ship their own Python — that is fine and
> supported. The live supervisor probes interpreters, finds the one that can
> `import SoapySDR`, and refuses to start if none can. Just know that the
> interpreter you run the doctor in may not be the one that opens the radio.

---

## 2. ffmpeg and a player

```
# Debian/Ubuntu
sudo apt-get install -y ffmpeg vlc
# macOS
brew install ffmpeg vlc
# Windows
winget install Gyan.FFmpeg ; winget install VideoLAN.VLC
```

**A trap worth 90 seconds of your life: VLC's option set differs per *build*,
not per version, and an unknown option is fatal.** Same 3.0.23 on two boxes, and
a flag that is a harmless no-op on one kills the other about a second after the
window opens — which looks exactly like a decode failure and is not.
`--no-qt-updates-notif` does not exist in distro builds (no update checker is
compiled in) and `--dummy-quiet` is Windows-only.

The doctor probes every flag we pass against *your* build:

```
python tools/atsc3_doctor.py       # look for "VLC flags on THIS build"
```

**ffmpeg has no AC-4 decoder**, in any build. That is expected — it is why this
project has its own — and it is why proof clips are transcoded to AAC (MP4 has
no tag to write AC-4 back out with).

---

## 3. The standards documents — optional

**A fresh clone cannot decode until you do this.** The symptom is unmistakable:

```
KeyError: '2026'
```

...raised at *import*, before any RF work happens at all.

### Why

The physical-layer modules parse constellation, LDPC, pilot and L1-syntax tables
**out of the published standards at import time**. Those tables are runtime
inputs, not documentation — the decoder cannot construct its constellations
without them. The AC-4 decoder does the same for its scalefactor-band tables,
reading them straight out of the ETSI PDF with pymupdf.

They are not in this repository because **they are free to download and not ours
to redistribute.** ATSC and ETSI both publish these documents at no charge; that
is what makes this project possible at all. Neither grants us the right to
re-host them.

ETSI returns HTTP 403 to automated fetching. **Working around an access control
is not on the table** — download by hand.

### What to fetch

| document | where it goes | needed for |
|---|---|---|
| ATSC A/322 (Physical Layer Protocol), 2024 **and** 2026 revisions | `lab/spec/` | video, everything |
| ATSC A/331 (Signaling, Delivery, Synchronization) | `lab/spec/` | services, guide |
| ETSI TS 103 190-1 (AC-4 Part 1) | `spec/` | audio only |

Both A/322 revisions are wanted, not one: they are cross-checked against each
other, and 78 parse identities have to agree before the tables are trusted.

### Extracting them

The parsers read text extracted two ways — `pdftotext -table` and
`pdftotext -layout` — into these filenames:

```
lab/spec/A322_2026_tbl.txt      pdftotext -table  A322-2026.pdf  A322_2026_tbl.txt
lab/spec/A322_2026.txt          pdftotext -layout A322-2026.pdf  A322_2026.txt
lab/spec/A322_2024_tbl.txt      (same, 2024 revision)
lab/spec/A322_2024.txt
lab/spec/A331_tbl.txt           pdftotext -table  A331.pdf       A331_tbl.txt
spec/ts_10319001v010301p.pdf    the ETSI PDF itself, unextracted
```

`pdftotext` ships with poppler (`apt install poppler-utils`, `brew install
poppler`, or the Windows poppler build).

Verify:

```
python tools/atsc3_doctor.py            # "spec tables (lab/spec)  5/5 present"
python lab/spec_ldpc.py --gate          # 78/78 parse identities
```

> **Open question for this release.** Committing the *parsed numeric tables* —
> arrays of integers, which are facts rather than ETSI prose — would delete this
> entire step and make a clone self-contained. It is a real licensing question
> and it has not been answered yet. Until it is, this step stands.

---

## 4. The SDR

Install your SDR's driver and SoapySDR support, then:

```
python tools/atsc3_doctor.py --radio
```

`--radio` only **enumerates**. It never tunes, because this repository is often
running an unattended session and the SDR is single-tenant.

Failures, in the order they actually happen:

**0 devices, but the vendor's own app works.** The vendor API directory is off
your `PATH`. Confirm the OS enumerates the device first, then fix the path.

**`enumerate()` raises.** Usually a wedged USB device.
* Linux: `kill -9` any stale chain, then restart the vendor's service
  (`systemctl restart <vendor>`), then relaunch. Killing a chain while it holds
  the device can take the driver session down with it and leave the *next* one
  wedged inside SoapySDR teardown.
* Windows: a ghost-enumerated device is **reboot-only**. Replugging does not
  clear it. Do not spend an hour on it.

**It worked yesterday and now nothing decodes.** Before you touch code: check
that nothing else is holding the radio (`python tools/atsc3_doctor.py` lists
other processes — it *reports*, it never kills, because an unfiltered process
census matches the tool doing the censusing).

### Gain — the setting that is worth 100% of your FEC

This is the single highest-value paragraph in this document.

A front end driven into overload produces **0.0% FEC** while a commercial tuner
on the *same antenna* decodes the same air perfectly. It does not look like
overload from the inside; it looks like a weak signal, because frames keep
advancing and nothing errors. Measured live, same minute, same air:

```
rfgain_sel = 2      0.0% FEC     0.39x     (dead)
rfgain_sel = 4    100.0% FEC     1.02x     (perfect)
rfgain_sel = 6    100.0% FEC     1.03x
```

**Two steps less LNA gain was the difference between a dead chain and a perfect
one.** If your FEC is at or near zero while the chain works hard, **try less
gain before you try anything else.**

Three things this project learned the expensive way:

1. **Gain is per carrier, not global.** Two transmitters on this rig want
   settings measured **6.04 dB apart**. One global value cannot serve both.
2. **Read the gain back.** An unverified control write is a hope, not a setting.
3. **On an overloaded front end, throwing signal away is therapeutic.** Adding a
   passive splitter here *improved* SNR by **+10.7 dB net** — signal −4.3 dB,
   noise −15.0 dB — because a strong nearby signal had been driving the front
   end into intermodulation. An LNA never helped this site, for the same reason.
   The starvation model says the opposite and is wrong in this case.

### Captures: check them before you trust them

If you bank IQ to a file, check it before you draw any conclusion from it:

```
python tools/atsc3_doctor.py --capture yourfile.cs16 --rate 6.912e6
```

* **Any samples at full scale ⇒ VOID.** A clipped capture makes every downstream
  analyser lie. This project once retracted an entire L1 verdict — "the link is
  dead" — that was really a measurement of its own clipping.
* **RMS in the hundreds ⇒ a disconnected input**, not a weak signal. Rename it
  `VOID_...` rather than deleting it; a banked negative is worth keeping.
* **Sample count must equal wall clock × sample rate**, or the file is not what
  it says it is.

---

## 5. The optional C LDPC kernel

Worth roughly **7× on the LDPC stage**. Everything works without it.

```
python lab/build_ldpc_kernel.py
python lab/gate_e54.py
```

Needs any C compiler — gcc, clang, `cc`, or MSVC. On Windows it finds MSYS2's
gcc or MSVC's build tools itself.

* **Build it per box. Never copy the binary between machines.** It is pure C ABI
  via ctypes with no `Python.h`, so one build serves every CPython *on that box*
  — and none on any other.
* The compiler flags are **IEEE-pinned on purpose**. Do not "optimize" them:
  `-ffast-math` would let the compiler reassociate float32 operations and
  un-gate the kernel against the numpy reference. gcc and MSVC builds are
  measured output-identical.
* MSYS2 trap: `C:\msys64\mingw64\bin` must **lead** your `PATH`, or `cc1` dies
  with no error message at all.
* `ATSC3_LDPC_KERNEL=0` forces the numpy fallback.

---

## 6. Threads — the trap that looks like slow hardware

**A BLAS thread pool underneath a Python loop of tiny operations is a throttle
wearing the costume of parallelism.**

Measured: an audio worker sat **eleven minutes at 480% CPU** — 34 OS threads for
one Python thread — inside the QMF synthesis loop, running at ~0.09× real time
on a machine that does 1.59×. The loop issues thousands of *tiny* matrix
multiplies per second and the BLAS fans each one across the pool; it is pure
handoff thrash. With `OPENBLAS_NUM_THREADS=1` the same pass took 12.6 s: **9.56×**.

```
export OPENBLAS_NUM_THREADS=1      # audio workers
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
```

The shipped launchers export their own pins (workers 1, chain 4). You need this
only if you run a worker by hand.

The doctor measures it on your box, both ways, and shows you the ratio:

```
[  ok  ] BLAS tiny-op behaviour   20k tiny sgemms: unpinned 0.04s vs pinned 0.04s (1.02x)
```

**A corollary that pins cannot fix:** environment variables only fix the process
that remembered to read them. A multi-threaded BLAS call from inside a worker
pool nests a pool inside a pool, and that does not merely run slow — it has
deadlocked eight workers, twice, once as a hang and once as an access violation.
On a hot path, prefer broadcasting to a matrix multiply that a pool could grab.

---

## 7. First tune

Replay first, if you have a banked capture — it proves the whole chain with no
RF variables:

```
python -m atsc3 watch --capture yourfile.cs16 --player none --secs 60
```

Then live:

```
python -m atsc3 watch --rf <channel>
```

Find your local NextGen TV channel first — a channel that carries it will have a
bootstrap, and the detector answers that question definitively without decoding
anything:

```
python selftest.py                 # 32 checks, ~8 s, no radio needed
python sweep_archive.py FILE       # is there a bootstrap in this capture?
```

Useful options: `--secs N`, `--record OUT.mp4`, `--player none`,
`--json OUT.json`.

---

## When it does not work

Run the doctor first. Then, in this order:

1. **Is the gain wrong?** (§4) 0% FEC with the chain working hard is the
   signature. Try less gain.
2. **Is something else holding the radio?** Single-tenant. The doctor lists
   other processes.
3. **Are the spec tables there?** `KeyError: '2026'` is §3, every time.
4. **Is the capture clipped?** (§4) A clipped file lies to every analyser
   downstream of it.
5. **Is a thread pool eating you?** (§6) 480% CPU on one Python thread.
6. **Did VLC die a second after opening?** (§2) An unsupported flag, not a
   decode failure.
7. **Wedged USB?** Reboot-only on Windows; service restart on Linux.
8. **Only now, suspect the sky.** And when you do, get an independent opinion —
   a second receiver on the same antenna is the instrument that caught the gain
   fault above. Without one you cannot distinguish "the air is bad" from "we are
   broken", and this project guessed wrong about that for an entire evening.

**A capture that is correctly acquired and cleanly demodulated but decodes
nothing is not evidence about your antenna.** A one-frame interleaver phase
offset produced exactly that picture here — every quality metric excellent,
every correctness metric dead — and cost a transmitter being written off as
unreachable for two weeks.
