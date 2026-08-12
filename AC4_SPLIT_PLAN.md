# The `ac4-decoder` split — plan

The AC-4 decoder was always intended to ship separately. This is the proposal
for how, and what it may honestly claim.

**Nothing here has been done.** No files were moved. No repository was created.

---

## Why split at all

Three independent reasons, and the third is the strongest:

1. **The audiences do not overlap.** People who want an AC-4 decoder mostly do
   not own an SDR. Making them clone a software television receiver — with its
   physical-layer tables, its SDR dependency and its 15-minute install — to get
   a codec is a tax on the wrong people.
2. **The dependencies do not overlap.** The decoder needs numpy, scipy and a PDF
   reader. It does not need SoapySDR, an SDR, an antenna, ffmpeg, VLC, or the
   ATSC standards documents.
3. **There is no open-source AC-4 decoder.** ffmpeg demuxes AC-4 and cannot
   decode it; VLC cannot even parse the sample entry. That gap is worth filling
   in a form somebody can actually depend on, and "vendored inside a hobby TV
   receiver" is not that form.

---

## The seam is already clean — measured, not assumed

The AC-4 cluster is 25 modules with **exactly one import that crosses out of
it**, and that import uses **exactly one function**:

```
lab/m17_ac4_walk.py:50   import m7_play as PL
lab/m17_ac4_walk.py:65       sc = PL.trun_scan(frag)
```

`trun_scan` is an fMP4 fragment-run scanner — an ISO-BMFF utility, not a
receiver concept. Everything else in the cluster imports only other cluster
members, numpy, scipy, and the standard library.

So the split does not need a refactor. It needs a boundary drawn where one
already is.

---

## What moves

**25 modules**, currently `lab/m*.py`, renamed to what they are:

| now | becomes | role |
|---|---|---|
| `m17_ac4_walk` | `ac4/frames.py` | frame walker, `samples()` |
| `m18_ac4_substream` | `ac4/substream.py` | substream sizing |
| `m19_ac4_toc`, `m20_ac4_toc2` | `ac4/toc.py` | table of contents |
| `m21_ac4_tools` | `ac4/tools.py` | tool flags |
| `m22_mdct` | `ac4/mdct.py` | MDCT |
| `m23_hcb` | `ac4/huffman.py` | Huffman codebooks |
| `m24_spectral` | `ac4/spectral.py` | spectral data |
| `m25_scalefac` | `ac4/scalefactors.py` | scalefactors |
| `m26_pcm` | `ac4/pcm.py` | PCM reconstruction |
| `m27_sfb` | `ac4/tables/sfb.py` | scalefactor-band tables (reads the spec) |
| `m28_channels` | `ac4/elements.py` | framing, channel elements, `_acpl3_core` |
| `m29_audio` | `ac4/audio.py` | dequant, M/S, noise fill |
| `m30_filterbank` | `ac4/filterbank.py` | KBD windows, block switching |
| `m31_aspx_bands` | `ac4/aspx/bands.py` | A-SPX band layout |
| `m32_aspx_parse` | `ac4/aspx/parse.py` | A-SPX side info |
| `m33_qmf` | `ac4/qmf.py` | 64-band complex QMF |
| `m34_aspx_env` | `ac4/aspx/envelope.py` | envelope decode |
| `m35_hfgen` | `ac4/aspx/hfgen.py` | HF generation |
| `m36_envadj` | `ac4/aspx/adjust.py` | envelope adjustment |
| `m37_render_hf` | `ac4/aspx/render.py` | joint L/R render |
| `m38_resample` | `ac4/resample.py` | the 1001/960 output resampler |
| `m42_ac4_stream` | `ac4/stream.py` | `Ac4Stream`, the driver |
| `m43_ac4_pair` | `ac4/elements_pair.py` | `channel_pair_element` |

Plus:

* `tools/atsc3_audio.py` → `ac4/cli.py` (`ac4-decode`), minus its
  receiver-specific lane bookkeeping.
* `trun_scan` and the box-walking helpers from `lab/m7_play.py` → `ac4/containers/mp4.py`.
  **Copied, not moved** — the receiver keeps its own copy, which it needs for
  video and captions anyway. This is a ~60-line duplication that buys a clean
  one-way dependency; that trade is worth it.
* `lab/M19_AC4_FIELDS.md` → the new repo's `docs/BITSTREAM.md`. It is a field-by-field
  map of the bitstream and it is the most valuable document in the split.
* The AC-4 gates and their negative controls (see below).

## What stays

* Everything physical-layer: bootstrap, L1, LDPC/BCH, interleavers, demappers.
* ALP, IP, MMTP, ROUTE, MPU assembly.
* The lane bookkeeping, the viewer, the supervision layer, the warden, the judge.
* `lab/m39_mux.py` — muxing decoded audio into a transport stream is the
  receiver's job, not the codec's.
* The whole experiment notebook. The AC-4 entries get *cited* from the new repo,
  not copied — and they need the scrub pass in the project scrub audit first.

## What is deliberately NOT in the split

**Anything that turns the codec into a receiver.** The temptation will be to let
`ac4-decoder` grow a "just tune a channel" convenience. It should not. The
dependency runs one way, and the moment it does not, both repos have to be
released together forever.

---

## How the two repos depend on each other

```
    ac4-decoder                 atsc3 (receiver)
    ───────────                 ────────────────
    numpy, scipy                numpy, scipy, SoapySDR, ffmpeg, VLC
    ETSI TS 103 190-1           ATSC A/322, A/331  (+ ac4-decoder)
         ▲
         └──────────────────────────── depends on ────┘
```

**One direction, always.** `atsc3` depends on `ac4-decoder`; `ac4-decoder`
depends on nothing of ours. A circular dependency between a codec and a receiver
would mean neither can be released without the other, which defeats the split.

The receiver's audio worker becomes:

```python
from ac4 import Ac4Stream          # was: import m42_ac4_stream
```

Pinned by version, so a codec change cannot silently alter the receiver's
byte-identical regression baselines. Those baselines — three RF33 md5s, verified
unchanged across both of the E77 fixes — should be **run in both repos**: in
`ac4-decoder` as its own regression, and in `atsc3` as an integration gate. The
same three numbers, checked from two directions.

**Development ergonomics:** during the transition, `pip install -e ../ac4-decoder`
in the receiver's environment. Do not vendor a copy; two copies of a decoder
diverge, and the fleet already has a law about a fix that reached only the
published copy.

---

## The API the new repo should expose

Small, and shaped by what the receiver actually needs:

```python
from ac4 import Ac4Stream, decode_frame, __version__

# 1. raw frames in, PCM out -- the core contract, no container involved
pcm, rate = decode_frame(frame_bytes, element="5_X")

# 2. a stream driver, which is what a real consumer wants
s = Ac4Stream(presentation=0, element="5_X")
for frame in frames:
    pcm = s.push(frame)          # None until a frame completes

# 3. optional container convenience
from ac4.containers.mp4 import ac4_frames
for frame in ac4_frames("segment.m4s"):
    ...
```

**The core contract takes bytes, not filenames.** The receiver has its frames in
memory already; making it write them to a file to decode them would be absurd.

---

## What its README may honestly claim

Held to the same standard as the receiver's: measured, with the caveats attached.

### It may claim

* **It decodes real broadcast AC-4** — not synthetic test vectors, not a
  round-trip against our own encoder. Off-air audio from **two different
  broadcasters using different encoder configurations**.
* **`5_X` (5.1) elements**, all six channels, verified as **discrete broadcast
  surround** rather than an upmix. The evidence is the channel correlation
  structure: front L–R +0.887, surrounds Ls–Rs +0.999, L–C −0.023.
* **`channel_pair_element`** (stereo).
* **`ASPX_ACPL_3` core mode** (`5_X_codec_mode` 4), and the `payload_base` field.
  Both were needed for the second broadcaster and neither had ever come up on the
  first: 0.0% → **100.0%** of frames on one programme, 3.5% → **100.0%** on the
  other.
* **A-SPX spectral extension, complete**: 64-band complex QMF with **78.49 dB**
  reconstruction and 577-sample group delay, envelope decode, HF generation,
  envelope adjustment.
* **The 1001/960 output resampler.** Without it the audio runs **4.27% fast**.
* **Regression-gated**: three real-air baselines byte-identical across changes,
  200/200 frames before and after.
* **Gated with negative controls that fail for the right reasons** — the passing
  case and two deliberately broken ones, where a broken decode fails on
  non-silence and crest factor while white noise fails on spectral centroid and
  HF/LF ratio. Two controls failing the same way would mean the gate only
  measures one thing.
* **Written from ETSI TS 103 190**, clause references given throughout, so any
  claim can be checked against the standard.

### It must state as a gap

* **No full 5.1 upmix for A-CPL-coded content.** What it renders is a correct
  **stereo** decode — the coded downmix, which is what the broadcaster sent, not
  a fold-down we computed. The upmix needs `acpl_data_2ch()` (framing plus five
  entropy-coded parameter sets over N parameter bands) and the QMF-domain A-CPL
  synthesis stage. **Scoped, not hand-waved:** one bitstream parser, modelled on
  the existing A-SPX Huffman machinery, and one synthesis stage.
* **Companding is signalled and not applied.** Additional harmonics are missing
  in under 1.5% of channel-frames. Limiter tie-breaking is partial. The SSF
  speech frontend is never exercised by the streams we have.
* **Two broadcasters is not "AC-4".** State the law that this decoder taught,
  because it is the honest frame for the whole project: *a decoder that has only
  ever met one broadcaster has only ever been tested against one set of encoder
  choices.* Every unimplemented `if this flag is set` branch is a station it
  cannot hear yet. The corollary is a roadmap, and publishing it as one invites
  exactly the bug reports the decoder needs.
* **It is Python and numpy.** It decodes faster than real time for stereo (1.59×
  measured) and **slower than real time for six channels** (0.62–0.91×), because
  the filterbank is a Python loop over many small windows. Anyone wanting live
  5.1 needs to vectorise it.
* **It reads the ETSI PDF at runtime** for scalefactor-band tables. Same
  redistribution question as the receiver — see `LICENSE_RATIONALE.md`.

### It must not claim

* Not "a complete AC-4 decoder". It is a decoder for the AC-4 that two US
  broadcasters actually transmit.
* No conformance or certification claim of any kind.
* Nothing about Dolby beyond naming the format descriptively.

### The bugs worth publishing in its README

These are the reason to trust it, not a reason to doubt it:

* **The balance bug.** `aspx_balance=1` means channel B carries a *pan ratio*,
  not a level. Three gates passed while the right channel was **20× too loud**,
  because all three measured the left. **Gate the symmetry, not the part.**
* **The shared constant.** Three subsystems agreed the frame length was 1536.
  1536 is the MDCT length. Audio ran 4.27% fast and cross-checking could not see
  it, because **a shared constant is invisible to cross-checking.**
* **The contaminated table.** The blocking bug on the first 5.1 decode was a
  corrupted Annex B table that passed **six internal gates**. What caught it was
  the printed table in the standard — an *independent* check, not a better
  internal one.
* **One byte.** `payload_base` is a single byte, and missing it shifts the entire
  element by eight bits: every field after it reads noise. The parse still
  "closed" — which is why *"the parse closed"* is never sufficient. Short frames
  once raised closure from 64.6% to 76.1% while every correctness gate collapsed.

---

## Naming, and the trademark care it needs

**`ac4-decoder`.** Descriptive, nominative, says what it decodes.

**"Dolby" must not appear in the name of the repository, the package, the module,
or the CLI.** A component in this fleet was already renamed for exactly this
reason. Naming a format to say what you are compatible with is nominative fair
use; putting a company's mark in your product name is not.

Open questions for the operator:

* **Is `ac4-decoder` itself too close?** "AC-4" is a Dolby mark. Referring to it
  is normal and necessary — ffmpeg has an `ac4` demuxer. Naming a *package*
  after it is a step further than referring to it. A more cautious alternative
  exists (a neutral project name, with "decodes AC-4 audio" as the description),
  at the cost of discoverability, which for a codec is most of the value.
  **Recommendation: keep `ac4-decoder`, add a disclaimer line in the package
  description and the README** — see the `NOTICE` draft in
  `LICENSE_RATIONALE.md`. But it is a judgement call and it is the operator's.
* **PyPI name.** Check availability and squatting before announcing anything.
* **Same licence as the receiver** (Apache-2.0), for the same patent reasons —
  more so, if anything: audio codecs are the most patent-dense software there is.

---

## Sequencing — and what must not happen

Nothing below is done. Ordered so nothing is public until the audit is.

1. **Scrub first.** The AC-4 modules and `M19_AC4_FIELDS.md` go through
   `tools/scrub_audit.py` before they move anywhere. They are lower risk than the
   receiver — a codec has no antenna — but `M19_AC4_FIELDS.md` carries absolute
   paths and the field docs cite specific broadcasts.
2. Create the new repository **private**.
3. Move + rename, carve `trun_scan` into `containers/mp4.py`, keep the git
   history if the operator wants attribution continuity — **but note that
   history carries the audit findings with it**, so the receiver's history
   rewrite (the project scrub audit) must happen *before* any history is copied,
   or the split repo inherits the leak.
4. Stand up the gates in the new repo. **The gates move with the code.** A
   decoder published without its negative controls is an assertion.
5. Point the receiver at it, pinned. Re-run the receiver's full gate set and
   confirm the three md5 baselines are unchanged. **If they move, the split
   broke something** — that is precisely what those baselines are for.
6. Only then consider making either repository public, and only through
   the publish checklist.

**The failure mode to avoid:** publishing `ac4-decoder` first because it is
smaller and cleaner, while the receiver still holds the unscrubbed history. The
two repos would share commits, an author identity, and a writing style. Splitting
does not split the exposure.
