# AC-4 DECODER — STATUS as of 8/07

The sections below this one are a chronological research log, appended as the
work happened. This is the summary: what the decoder does, what proves it, and
what it does not do.

## WHAT IT DOES

**The MDCT core is complete and produces audio.** Every frame of the RF33
capture decodes, all six channels of the 5.1 element, and renders to
continuous, real-time-exact PCM.

    3480 / 3480 frames   0 overruns, 0 errors
    111.4 s continuous stereo, real-time exact, band-limited at 12 kHz
    long frames, short frames (54 % of transform windows are partial
    blocks), and i-frames all decode

**The A-SPX side information is completely parsed**, though it does not yet
make sound.

    3480 / 3480 frames close EXACTLY to the encoder's declared audio_size

## WHAT PROVES IT

Gates that a plausible wrong answer could not also pass:

    audio_data closes to audio_size          every bit accounted for, no slack
    per-channel adjacent-frame correlation   lfe +0.30  L +0.43  R +0.75
                                             Ls +0.73  Rs +0.50  C +0.58
                                             (shuffled controls all ~0)
    pair coupling L/R                        -0.43, control ~0 (mid/side, so
                                             the sign is expected)
    rendered audio band limit                100.00 % below 12 kHz
    L/R in the render                        +0.9829, separation 0.105
    MDCT block switching                     S @ S.T == I to 4.4e-16 across
                                             every switching sequence
    QMF analysis/synthesis                   78.49 dB reconstruction
    A-SPX envelopes                          adjacent +0.7605, shuffled +0.0244
    HF patch copy                            destination == source, rel 0.0e+00

**Three subsystems that share no code agree the crossover is 12000 Hz**: the
Annex B band table (max_sfb 43 -> line 768 of 1536), the A-SPX header
(start_freq 5 -> QMF subband 32), and the QMF bank measuring our own decoded
audio's band edge. Nothing forces that agreement.

## WHAT IT DOES NOT DO

**A-SPX makes no sound.** 12–21 kHz is parsed, its geometry derived and its
envelopes decoded and validated, and the HF patch copy is implemented and
gated — but the four-step high-band chain stops after step 1:

    1. HF generation        DONE and gated (m35_hfgen)
    2. envelope adjustment  NOT STARTED — Pseudocode 90, 91, 94, 95
    3. noise addition       NOT STARTED — Pseudocode 83, 94
    4. additional harmonics NOT STARTED — Pseudocode 92, 93
       plus the limiter     NOT STARTED — Pseudocode 72 (sbg_lim), 96

The patched band currently sits at the WRONG LEVEL, which is why m35 writes no
audio: a copy of 7.5–12 kHz is not what 12–21 kHz should sound like, and
rendering it would produce a file that plays and is wrong.

Also not implemented, and counted rather than hidden: **A-CPL** (not signalled
in this stream), **companding** (signalled, not applied), **sap_mode 3 full SAP**
(91 frames, passed through un-decorrelated), **spectral noise fill** (implemented
but `b_snf_data_exists = 0` on every frame — a structural zero here), and the
**SSF speech frontend** (never used by this stream).

One value is inferred rather than read: the A-SPX **F0 codebook offset**. Annex A
gives no `cb_off` for the A-SPX books, so the delta books are centred from their
code lengths and the F0 books used raw. It shifts absolute high-band level, not
shape, and does not affect any gate above.

## 8/07 — A-SPX, PART 7: ENVELOPE ADJUSTMENT AND NOISE (m36_envadj.py)

    GATE 1  gain reproduces the transmitted envelope
            relative error median 8.84e-09, max 5.16e-07          PASS
    GATE 2  limited+boosted / unlimited gain worst case 1.417572
            MAX_BOOST_FACT (Pseudocode 100)          1.584893     PASS
    GATE 3  envelope error signal only        0.2028
            envelope error signal + noise      0.1108             PASS

Gate 1 is a closed loop with no free parameters: estimate the patched band
(Pseudocode 90), dequantise the transmitted envelope (82), map it onto QMF
subbands (91), compute the gain (95), apply it, then RE-MEASURE. Target and
achieved agree to float64 precision.

**THE BLOCKER WAS MY OWN PIPELINE, NOT THE SPEC.** The first run read 0.9998.
m29 peak-normalises its render (0.9 / peak), which destroys exactly the
absolute level the A-SPX envelopes are defined against; Pseudocode 95's
`EPSILON = 1.0` then dominates the denominator, and with est_sig ~ 1e-4 the
gain saturates and the band comes out ~1e6 too quiet. Re-rendering the core at
its NATIVE scale took the error to 0.2028, and isolating the limiter took the
gain path to zero.

**TWO GATES I WROTE WRONG, BOTH FAILING CORRECT CODE.**

- A 1e-9 threshold is tighter than float64 through a sqrt and a square can
  deliver. 8.84e-09 IS exact.
- I asserted "the limiter only ever reduces the gain". It does not: the boost
  is computed PER LIMITER SUBBAND GROUP and applied to every subband in it, so
  an individual gain can exceed its unlimited value while the group's energy is
  restored. What is actually bounded is the boost, at MAX_BOOST_FACT, and the
  measured worst case sat exactly on it.

That is the fourth and fifth time in this project a gate encoding my
expectation has failed correct code (after M24's ">20 % zeros", M29's stereo
sign, and M33's window symmetry).

**NOISE ADDITION** (Pseudocode 83/94/97, Table D.2). The noise envelopes are
reconstructed with the same delta scheme as the signal ones, then dequantised
as `2 ** (NOISE_FLOOR_OFFSET - qscf)` -- note the SIGN, a larger transmitted
index means LESS noise. `ASPX_NOISE[512][2]` from the spec C file is verified:
**mean |z|^2 = 1.000000 exactly**, every entry unit magnitude, phase spread
1.787 rad against uniform's 1.814.

Adding it halves the residual envelope error, 0.2028 -> 0.1108, and it also
relaxed the boost (worst case 1.5849 -> 1.4176) because the noise now supplies
energy the boost previously had to invent.

**Remaining:** additional harmonics (<1.5 % of channel-frames), the full
tie-breaking in the limiter subband table (Pseudocode 72), and wiring the
adjusted high band back through QMF synthesis into the render. Audio is still
written band-limited at 12 kHz -- no high band has been rendered to a file.

## 8/07 — A-SPX, PART 8: THE HIGH BAND IS RENDERED (m37_render_hf.py)

The chain is wired end to end and writes audio containing 12..21 kHz:
decode the core at NATIVE scale -> QMF analyse -> HF generate -> envelope
adjust + noise -> QMF synthesise -> write.

    GATE 1  12..21 kHz   core 0.0000 %   with A-SPX 0.1408 %      PASS
    GATE 2  21..24 kHz   core 0.0000 %   with A-SPX 0.0000 %      PASS
    GATE 3  correlation with the core, delay-aligned  +0.9996     PASS

    delivered spectrum (38.4 s):
       6000-12000 Hz  2.2070 %     core, MDCT
      12000-16500 Hz  0.1416 %     A-SPX patch 0
      16500-21000 Hz  0.0017 %     A-SPX patch 1
      21000-24000 Hz  0.0008 %     nothing, as A-SPX stops at 21 kHz

Gate 2 matters as much as gate 1: A-SPX stops at QMF subband 56, so a correct
high band is BOUNDED at 21 kHz. Energy above it would mean the patch or the
synthesis is wrong, and "the file sounds fuller" would not have caught that.

Gate 3 guards the other failure mode -- it is easy to produce a fuller-sounding
file by breaking the core. Correlation with the core render is +0.9996 once the
QMF's 577-sample delay is applied.

The rolloff is continuous across 12 kHz with no step at the crossover, and
patch 1 is far quieter than patch 0, which matches the -2.5 steps/subband-group
envelope slope M34 measured independently.

**THE ALIGNMENT MISTAKE, A THIRD TIME.** Gate 3 first read a correlation of
-0.0178 -- essentially zero -- on a render that was fine, because I compared
`core[:n]` against `hf[:n]` while the QMF round trip delays the output by 577
samples. That is the same error as M33's SNR gate and M29's run-interior
measurement. The delay is a MEASURED property of the filterbank, not an
unknown, so it is now applied rather than searched for.

**Outstanding.** The A-SPX control data for frame i is applied to that frame's
nominal QMF slots without compensating for the analysis delay, so the high band
can be smeared by up to ~9 time slots against the core. That is a TEMPORAL
error, not a spectral one; the gates above are spectral and insensitive to it,
and M36's envelope-accuracy gate uses the same slots for measurement and
application so it is unaffected. Also outstanding: additional harmonics
(<1.5 % of channel-frames) and the full limiter-table tie-breaking.

## 8/07 — THE BALANCE BUG: L PERFECT, R BROKEN, THREE GATES GREEN

The first full high-band render passed gates 1, 2 and 3 and was WRONG. Its
right channel was 20x too loud and uncorrelated with its own core.

    L output vs L core, delay-aligned   +0.9999    rms 3.23e4  (core 3.23e4)
    R output vs R core, delay-aligned   +0.0519    rms 6.28e5  (core 3.26e4)
    output L/R correlation              +0.0504    (core L/R +0.9713)

**Cause.** When `aspx_balance = 1` the pair is coded as a SUM channel and a
BALANCE channel, and channel B's transmitted scale factors are a PAN RATIO, not
a level. Pseudocode 82 -- the dequantisation I was using -- says in its own
preamble that it covers "a 1-channel element and a 2-channel element with
aspx_balance = 0". Pseudocode 84 is the joint version:

    PAN_OFFSET = 12
    nom     = 2 ** (qscf_a/a + 1) * num_qmf_subbands
    denom_a = 1 + 2 ** (PAN_OFFSET - qscf_b/a)
    denom_b = 1 + 2 ** (qscf_b/a - PAN_OFFSET)

Running Pseudocode 82 on channel B treats a pan ratio as an absolute level.
**89.7 % of frames in this stream use balance = 1**, so the right channel's
high band was wrong nearly everywhere.

**Why the gates missed it.** Gates 1-3 all measure the LEFT channel: does the
high band appear, does it stop at 21 kHz, is the core undisturbed. Channel A is
the sum channel and decodes correctly under Pseudocode 82, so every
single-channel measurement looked perfect. The bug was invisible by
construction.

The fix required restructuring the renderer: the pair cannot be processed one
channel at a time, because Pseudocode 84 needs BOTH channels' indices at once.

    after the fix   L/R correlation +0.9710, rms 4166 vs 4194

**GATE 4 now checks the SYMMETRY, not a channel.** L/R correlation and the R/L
rms ratio. That is the general lesson and it is worth stating plainly: when a
system has two parts that should mirror each other, gate the RELATIONSHIP
between them. Gating each part separately can pass while the pair is broken --
here one part was not merely correct, it was exemplary, which is what made the
green gates so convincing.

This is the sixth gate in this project that failed to constrain what it looked
like it constrained (after M24's ">20 % zeros", M29's stereo sign, M33's window
symmetry, and M36's two thresholds) -- but the first where the gates PASSED on
broken output rather than failing on correct output. That is the more dangerous
direction.

## 8/07 — THE SECOND RENDERER BUG: PER-CHANNEL FRAMING

After the balance fix, the full render reported "A-SPX applied to 3366 frames"
of 3480 -- 114 frames silently fell through to core-only. A per-frame data
diagnostic said all 3480 frames had usable A-SPX data, which localised it: the
DATA was fine, the RENDERER was wrong.

Cause: with `aspx_balance = 0` (about 10 % of frames) the two channels have
INDEPENDENT framings, so channel 1 can have a different envelope-border list
from channel 0. The renderer used channel 0's borders for both, and the shape
check then rejected the frame. Each channel now uses its own borders; with
`balance = 1` they legitimately share one, because the framing is transmitted
once.

    before   3366 / 3480 frames received A-SPX
    after     600 /  600 in the verification run, all four gates green

Worth noting how this was found: the count "3366 of 3480" was printed by the
renderer itself and did not fail any gate. A gate on COVERAGE -- every frame
that has data should get the treatment -- would have caught it directly, and
that is the same class as the earlier lesson about coverage vs correctness,
pointed the other way: there, coverage rose while correctness fell; here,
correctness was fine while coverage quietly dropped.

## VERIFICATION (8/07)

Two properties checked that the gates above do not cover.

**The stateful parse self-synchronises.** Both the A-SPX framing borders and
the envelope reconstruction carry state between frames (`previous_stop_pos`,
`qscf_sig_sbg_prev`), which is the first thing in this decoder that cannot be
evaluated on a single frame. Parsing from frame 0 and parsing from frame 900
give **bit-identical end positions for all 2280 frames from index 1200 on** --
zero differences. So the cross-frame state re-anchors itself and the decoder
can be started mid-stream, which is exactly the spliced-bitstream case clause
5.7.6.3.1 warns about.

**The render is deterministic.** Re-rendering the whole file reproduces
`tv_audio_valid.wav` byte-for-byte (sha256 8fa45c10d52deefd, 21381164 bytes),
so the delivered audio is exactly what the current code produces -- no drift
from the KBD-window fix, the short-frame fix or the i-frame fix, and no hidden
nondeterminism from the noise-fill RNG (which is seeded per frame from
sequence_counter).

## FILES

    m22 MDCT          m23 codebooks     m24 spectral      m25 scale factors
    m26 LFE PCM       m27 Annex B       m28 channels      m29 audio render
    m30 filterbank    m31 A-SPX bands   m32 A-SPX parse   m33 QMF bank
    m34 A-SPX env     m35 HF generator

All 14 pass their gates as of this writing. Nothing has been pushed —
the atsc3 repo still needs the history rewrite for the three commits carrying
multicast addresses, and your go-ahead.

---

# M19 — AC-4 TOC field widths, mined from the spec

Working notes for the parser, written as **field names and bit widths** — the
facts needed to implement — rather than as a copy of the standard's text. The
spec itself lives in `spec/` and is gitignored: free to download from ETSI, not
ours to redistribute. Section numbers are given so anything here can be checked
against the source.

Source: ETSI TS 103 190-1 V1.3.1, clause 4.2. Our stream is
`bitstream_version = 2`, `fs_index = 1` (48 kHz), `frame_rate_index = 3`
(29.97 Hz), `b_single_presentation = 1`.

## raw_ac4_frame (4.2.1)

    ac4_toc()
    fill_area, byte_align
    for s in range(n_substreams): ac4_substream_data()   # stays byte aligned
    fill_area, byte_align

**This is the oracle for the whole parse**: the substream sizes from the TOC
must account for the frame, and the container gives us each frame's true
length. Get the TOC wrong and the sizes will not add up.

## variable_bits(n) (4.2.2)

    value = 0
    loop:
        value += read(n)
        more = read(1)
        if not more: return value
        value <<= n
        value += (1 << n)

**Not** a plain shift-or accumulate — the offset on each continuation is what
keeps successive ranges disjoint. I had this wrong from memory.

## ac4_toc (4.2.3.1)

| field | bits | note |
|---|---|---|
| bitstream_version | 2 | escape: `== 3` then `+= variable_bits(2)` |
| sequence_counter | 10 | +1 mod 1024 per frame — the walk gate |
| b_wait_frames | 1 | |
| wait_frames | 3 | only if b_wait_frames |
| reserved | 2 | only if wait_frames > 0 (**not** br_code) |
| fs_index | 1 | |
| frame_rate_index | 4 | |
| b_iframe_global | 1 | our bit 23 |
| b_single_presentation | 1 | 1 ⇒ n_presentations = 1 |
| b_more_presentations | 1 | only if not single; then `variable_bits(2) + 2` |
| b_payload_base | 1 | |
| payload_base_minus1 | 5 | if b_payload_base; `payload_base = v + 1` |
| (payload_base escape) | | if payload_base == 0x20: `+= variable_bits(3)` |
| ac4_presentation_info() | | × n_presentations |
| substream_index_table() | | |
| byte_align | 0–7 | |

## ac4_presentation_info (4.2.3.2)

    b_single_substream                        1
    if not single: presentation_config        3   (== 7 -> += variable_bits(2))
    presentation_version()
    if not single and presentation_config == 6:
        b_add_emdf_substreams = 1
    else:
        mdcompat                              3
        b_belongs_to_presentation_id          1
        if set: presentation_id = variable_bits(2)
        frame_rate_multiply_info(frame_rate_index)
        emdf_info()
        if b_single_substream: ac4_substream_info()
        else: b_hsf_ext (1), then per presentation_config a sequence of
              ac4_substream_info() / ac4_hsf_ext_substream_info()

## presentation_version (4.2.3.3)

Unary: count 1-bits until a 0. Value = number of 1s.

## frame_rate_multiply_info (4.2.3.4)

    frame_rate_index in {2,3,4}:  b_multiplier (1); if set: multiplier_bit (1)
    frame_rate_index in {0,1,7,8,9}: b_multiplier (1)
    otherwise: nothing

**Ours is 3**, so: 1 bit, plus 1 more if that bit is set.

## emdf_info (4.2.3.5)

    emdf_version        2    (== 3 -> += variable_bits(2))
    key_id              3    (== 7 -> += variable_bits(3))
    b_emdf_payloads_substream_info  1
        if set: emdf_payloads_substream_info()
    emdf_protection()

## emdf_payloads_substream_info (4.2.3.10)

    substream_index     2    (== 3 -> += variable_bits(2))

## emdf_protection (4.3, Table 80)

    protection_length_primary     2
    protection_length_secondary   2
    protection_bits_primary       8 / 32 / 128
    protection_bits_secondary     0 / 8 / 32 / 128

**UNKNOWN #1**: the mapping from the 2-bit length codes to those bit counts is
not yet extracted. Assumption to test: `{0:0, 1:8, 2:32, 3:128}`.

## ac4_substream_info (4.2.3.6)

    channel_mode        1/2/4/7   variable-length code
                                  (== 0b1111111 -> += variable_bits(2))
    if fs_index == 1:  b_sf_multiplier (1); if set: sf_multiplier (1)
    b_bitrate_info      1 ; if set: bitrate_indicator  3/5
    if channel_mode in {0b1111010, 0b1111011, 0b1111100, 0b1111101}:
        add_ch_base     1
    b_content_type      1 ; if set: content_type()
    for i in range(frame_rate_factor): b_iframe (1)
    substream_index     2   (== 3 -> += variable_bits(2))

**UNKNOWN #2**: the `channel_mode` code table (which prefixes are 1, 2, 4 or 7
bits, and what they mean). **UNKNOWN #3**: how `frame_rate_factor` is derived —
it controls how many `b_iframe` bits follow, and our entropy map found FOUR
bits carrying the I-frame flag, so this is very likely where three of them come
from.

## substream_index_table (4.2.3.11)

    if n_substreams == 0:  (2 bits)   n_substreams = variable_bits(2) + 4
    if n_substreams == 1:  b_size_present (1)   else b_size_present = 1
    if b_size_present:
        for s in range(n_substreams):
            b_more_bits        1
            substream_size[s]  10
            if b_more_bits: substream_size[s] += (variable_bits(2) << 10)

This is exactly why M18's search for a fixed-offset size field found nothing:
the size is 10 bits plus a variable-length extension, and it sits after a
variable-length presentation_info.

## The three unknowns, and the gate

Resolve UNKNOWN #1–#3 and the parse closes. The test is not opinion: **the
substream sizes must sum to the frame's byte count**, on all 3480 frames, and
the four `b_iframe` bits must land at 23, 95, 108, 147 where the entropy map
already found them. A parse that misses either is wrong.

---

# CORRECTION — part 1's presentation syntax does not apply to our stream

The three unknowns resolved cleanly:

* **protection_length_primary** (Table 175): `00` Reserved, `01` 8, `10` 32,
  `11` 128. **secondary** (Table 176): `00` 0, `01` 8, `10` 32, `11` 128.
  My assumption was right for secondary and wrong for primary — `00` is
  Reserved there, not zero.
* **channel_mode** (Table 88) is a prefix code: `0` Mono, `10` Stereo,
  `1100`/`1101`/`1110` = 3.0/5.0/5.1, `1111000`…`1111101` = the six 7.x
  layouts, `1111110` Reserved, `1111111` escapes to `variable_bits(2)`.
* **frame_rate_factor** (Table 87): for `frame_rate_index` 2/3/4 — ours is 3 —
  it is 1 when `b_multiplier` is 0, else `multiplier_bit` picks 2 or 4. It is
  the loop count for `b_iframe` in each `ac4_substream_info()`.

**Then the parser was built and GATE 1 failed on frame 0**, and the failure
diagnosed itself. The trace decoded `presentation_config = 6`, and part 1's
Table 85 says **`presentation_config ≥ 6` is Reserved** — so the stream was not
exotic, the parse was misaligned.

**The cause: part 1's `ac4_presentation_info()` is not our syntax.** TS 103
190-2 clause 6.2.1 is explicit — the presentation information is carried in
`ac4_presentation_info()` **if bitstream_version ≤ 1** and in
`ac4_presentation_v1_info()` **otherwise**. Our `bitstream_version` is 2.

Part 2 also carries its own `ac4_toc()` (6.2.1.1) with a
`total_n_substream_groups` loop over `ac4_substream_group_info()` — a structure
that does not exist in part 1 at all.

**This is worth recording as more than a bug.** The original plan named
`ac4_presentation_v1_info` from memory. Then part 1 arrived, part 1 has a
function with almost the same name, and I used it — trading a correct recollection
for an authoritative-looking wrong one. Having the document is not the same as
reading the right clause of it. The gate is what caught it, which is exactly
why the gate was written before the parser.

**Next: implement part 2's `ac4_toc` + `ac4_presentation_v1_info` +
`ac4_substream_group_info`.** The two gates are unchanged — sizes must sum to
the frame, and the `b_iframe` bits must land at 23, 95, 108, 147.

---

# LFE-first, and the re-anchored grid (user's idea, 8/06 night)

## Is the signal even unencrypted?  Verified, three ways

Worth asking before spending a night on it.

1. **The container carries no encryption at all.** No `enca`, `sinf`, `schm`,
   `schi`, `tenc`, `senc`, `saiz`, `saio` or `pssh`. The sample entry is plain
   **`ac-4`**. Under MPEG Common Encryption an encrypted audio track *must* be
   `enca` wrapping a `sinf` — it is not optional.
2. **A field INSIDE the audio payload is constant.** `5_X_codec_mode` sits
   within `audio_data_chan`, past any encryption boundary, and reads 1 on all
   3480 frames. Ciphertext would be ~uniform over 8 values; P(constant) = 8^-3479
   ≈ 10^-3142.
3. **The video of the same service decodes** — 114 s, zero ffmpeg errors.

Noted for honesty: payload byte entropy is ~8 bits, but that is true of
compressed *and* encrypted data, so it is not diagnostic either way.

## The numbers the LFE path needed

* **Table 83**: `frame_rate_index 3` (29.97 fps) → **internal frame length
  1536 samples** at 48 kHz, with decoder resampling ratio 1001/1000 × 25/24.
  1536 is the MDCT length, and the resampler is why 1536 internal samples
  become 1601.6 external ones per 33.367 ms frame.
* **Table 105**: transform length 1536 → `n_msfb_bits` 6, `n_side_bits` 5,
  **`n_msfbl_bits` 3**.
* Read with 3 bits, the LFE's `max_sfb` is **3 on 3183 of 3248** non-I-frames.
  (I-frames excluded: `aspx_config()` sits between `5_X_codec_mode` and
  `mono_data`.)

## THE GRID, re-anchored — and why the first one was half-blind

The original entropy map indexed by **absolute position from the frame start**.
That works only while fields sit at fixed offsets, which is true of the TOC and
false of everything after it. Re-anchoring to **each structure's own origin**
is what makes the payload legible:

    entropy vs bit offset from audio_data_chan start, 3248 non-I-frames

    +  0   0.00 0.00 0.00 0.00 0.08 0.11 0.08 0.04 0.05 0.16 0.00 0.00 0.00 0.14 0.08
    + 15   0.87 0.96 0.97 0.98 1.00 1.00 1.00 ...  and 1.00 from there on

**The structured header is exactly bits 0..14. The entropy coder begins at bit
15.** A one-bit boundary, located statistically with no parsing.

Accounted so far: `5_X_codec_mode` (3) + LFE `max_sfb` (3) = 6 bits, leaving
**9 structured bits** before the coder starts — which is where
`asf_section_data()` begins (`sect_cb` is 4 bits, then section-length bits).

That is the AC-4 analogue of the OFDM grid: not a picture of the waveform, but
a map that says *where the addressable structure is*, so a parse can be checked
against something other than its own confidence.

---

# The 9 structured bits, and where the Huffman tables actually live

## The grid's boundary is confirmed by the syntax, to the bit

`asf_section_data()` reads `sect_cb` (4 bits) then `sect_len_incr`
(`n_sect_bits`, which is 5 for a long frame). So one section costs 9 bits:

    5_X_codec_mode   3
    LFE max_sfb      3
    sect_cb          4
    sect_len_incr    5
    -------------------
                    15   <- and the entropy map, built with no spec at all,
                            put the entropy coder's start at bit 15.

Measured over 3248 non-I-frames:

    sect_cb           {11: 3160, 10: 50, 2: 10, 8: 7, 6: 7, 1: 6, ...}
    sect_len (1+incr) {3: 3183, 2: 33, 1: 32}
    sect_len covers max_sfb   3248/3248

Two independent methods — a statistic and a syntax — agreeing on the same bit,
and a semantic check (the section must span `max_sfb`) passing on every frame.

## The codebook this stream actually uses

`sect_cb = 11` on 97 % of frames, so **ASF_HCB_11 is the first table dolbyTuna
needs**, and its shape is fully specified in Annex A even though the table is
not:

| codebook | length | shape | values |
|---|---|---|---|
| ASF_HCB_1..4 | 81 = 3^4 | quads (CB_DIM 4) | -1..+1 (cb_off 1) or 0..2 |
| ASF_HCB_5,6 | 81 = 9^2 | pairs | -4..+4 (cb_off 4) |
| ASF_HCB_7,8 | 64 = 8^2 | pairs | 0..7, unsigned + sign bits |
| ASF_HCB_9,10 | 169 = 13^2 | pairs | 0..12, unsigned |
| **ASF_HCB_11** | **289 = 17^2** | **pairs** | **0..16, unsigned** |
| ASF_HCB_SCALEFAC | 121 | scalefactors | |
| ASF_HCB_SNF | 22 | | |
| ASPX_HCB_* | 13..141 | the A-SPX envelope/noise books | |

The codebook_length values are all perfect powers, which is itself a check on
the extraction: 289 = 17^2 and 81 = 3^4 are not coincidences, they are the
alphabet size raised to CB_DIM.

## THE TABLES ARE NOT IN THE PDF — and that is good news

Annex A, verbatim: *"All Huffman codebook tables are available in the file
`ts_103190_tables.c` contained in archive `ts_10319001v010301p0.zip` which
accompanies the present document."*

A C source file beats PDF table extraction in every way — machine readable,
complete, no OCR ambiguity, and each table arrives as an explicit
`_LEN` / `_CW` pair. The archive sits next to the PDF on the ETSI server:

    .../10319001/01.03.01_60/ts_10319001v010301p0.zip

**It needs one more manual download** (ETSI 403s automated fetches), into
`spec\`. Once it is there the Kraft check becomes a per-table
pass/fail: a valid prefix code satisfies `sum(2^-len) == 1` exactly, so a
mis-read table is caught immediately rather than producing confident noise.

---

# Past the LFE: what the five full-band channels actually are

Parsing continues past `mono_data(1)` through `asf_snf_data`,
`companding_control(5)` and `coding_config`.  All three come out CONSTANT over
3160 frames, which is what a consistent broadcast encoder produces:

    b_snf_data_exists   {0: 3160}              no spectral noise fill
    companding          sync 0, all 5 OFF      COMPANDING IS NOT USED
    coding_config       {0: 3160}              = 2ch + 2ch + mono

**Two more tools fall out of scope.** M21 listed companding as required
because `5_X_codec_mode == ASPX` selects `companding_control(5)` — but the
control itself signals every channel OFF on every frame, so the companding
*transform* is never applied.  And `b_snf_data_exists` is 0 throughout, so
spectral noise fill is never applied either.  Both were read from the stream,
not assumed.

`asf_snf_data` also costs exactly 1 bit here for a reason worth noting: its
loop body only fires for bands that are EMPTY (`sfb_cb == 0` or
`max_quant_idx == 0`), and all three LFE bands are populated.  So the flag is
read and the loop does nothing.

## The remaining structure, from `coding_config = 0`

    2ch_mode                1 bit
    two_channel_data()      L / R
    two_channel_data()      Ls / Rs
    mono_data(0)            C

and `two_channel_data()` is:

    b_enable_mdct_stereo_proc   1
    if set:   sf_info(ASF,0,0) ; chparam_info()
    else:     sf_info(ASF,0,0) ; sf_info(ASF,0,0)
    sf_data(ASF)
    sf_data(ASF)

So the full-band channels reuse `sf_data` — the same section/spectral/scalefac/
snf chain already decoded and gated for the LFE.  What is NEW is only:

  * `sf_info(ASF, ...)` = `asf_transform_info()` + `asf_psy_info()`, i.e. the
    long/short window decision and `max_sfb` at `n_msfb_bits = 6` (Table 105
    for transform length 1536) rather than the LFE's 3-bit field;
  * `chparam_info()` when MDCT stereo processing is on;
  * `mono_data(0)`, which reads a 1-bit `spec_frontend` and can select SSF
    instead of ASF — the LFE path forces ASF and skips that choice;
  * A-SPX afterwards, for content above the coded band.

**The band count changes everything about size.**  The LFE used `max_sfb = 3`
= 12 lines.  A full-band channel at `n_msfb_bits = 6` can carry up to 63 bands,
which at this transform length is most of 1536 lines — so these channels are
where the bitrate actually goes, and where a bit error would be far more
audible than in a 12-line rumble.

---

# sf_info for the full-band channels

`n_grp_bits` (4.3.6.2.4) is **0 when `b_long_frame` is true and
`frame_len_base >= 1536`** — which is us.  So a long-frame `sf_info(ASF,0,0)`
is just:

    b_long_frame     1 bit   (asf_transform_info)
    max_sfb[0]       6 bits  (asf_psy_info, n_msfb_bits = 6 per Table 105)
    -------------------------
                     7 bits, and nothing else

Read on the first full-band channel, over the frames that get that far:

    2ch_mode                   {0: 2301}
    b_enable_mdct_stereo_proc  {1: 2301}     -> sf_info once + chparam_info()
    b_long_frame               {1: 2301}
    max_sfb                    {43: 2235, 42: 29, 41: 26, 37..40: 11}

**`max_sfb = 43` is the headline and it validates.**  Table B.4 gives transform
length 1536 exactly **49 scale factor bands** (sfb 49 -> offset 1536, the whole
transform), so 43 is inside the legal range and near the top — which is what a
full-band channel should look like, with A-SPX carrying whatever sits above the
coded edge.  Compare the LFE's `max_sfb = 3`.

**859 frames report `b_long_frame = 0`** — genuine short/transient frames, 27 %
of the total.  Those need the other branch: two `transf_length` fields, up to
two `max_sfb` values, and up to 15 grouping bits from Table 108.  Not
implemented yet, and counted rather than hidden.

**Honest gap — CLOSED, then REOPENED, then closed properly (M27/M28).**

The first "close" was wrong and is worth recording in full, because it passed
its gates and broke something else silently.

**What I got right.** Table B.4 is printed two table rows per printed line
(sfb n left, sfb n+56 right) and the spec sets thousands as `1 024`, i.e. two
word tokens — so an index-based parse reads a different column depending on how
many values on that line happen to be four digits. Binning by x-coordinate
fixes that, and it was the right call.

**What I got wrong.** I restricted the scan to "pages 258–259" to avoid
contamination from neighbouring tables. But Table B.4's last row (sfb 55) is
the *first* row on page 259, and **Table B.5 starts immediately below it on the
same page**, its sfb column restarting at 0. So B.5's rows overwrote B.4's from
sfb 40 up. The structure I was reading does not respect page boundaries, and
narrowing the page range is not a fix for that — the scan now stops when the
sfb column stops increasing.

**Why it survived.** Three reasons, all instructive:

- *Every gate was internal.* Monotone, all widths multiples of 4, sums to 1536,
  no gaps. The contaminated table passed all six — because it was assembled
  from another **valid** band table and inherited every structural property a
  band layout has. Internal consistency cannot detect a coherently wrong table.
- *The corroboration was manufactured by the bug.* The bad widths ended
  `32x14, 128x6` and `max_sfb = 43` landed exactly on the 32→128 boundary. I
  wrote that up as independent confirmation — "the encoder stops where the fine
  resolution stops." Pure coincidence. Real widths end `52x2, 64x7`; 43 lands
  nowhere special. A satisfying story that explains the result is not evidence
  for it.
- *The LFE never noticed.* The tables agree up to sfb 10 and the LFE uses bands
  0–2, so M24/M25/M26 stayed correct and kept passing. A bug that spares the
  component under test is invisible until something is built on top of it.

**What caught it: Table B.1**, three pages earlier, which `asf_section_data`'s
own NOTE points at. It states num_sfb per transform length directly:

    num_sfb(1536 @ 48 kHz) = 55        <- ours
    num_sfb(1536 @ 96 kHz) = 49        <- what I had

My 49 was a real number from a real table, just not our sampling rate — and I
had derived it from where my own column crossed 1536, then used it to validate
that same column. Circular, and the circle closed cleanly enough that nothing
complained. That check is now the first gate.

**The corrected table** (all three columns of B.4 are *identical* over sfb
0..55 and diverge only above it, so the "which column is 1536?" question that
drove the first attempt never needed answering in our range):

    0 4 8 12 16 20 24 28 32 36 40 44 52 60 68 76 84 92 100 108 116 124 136 148
    160 172 188 204 220 240 260 284 308 336 364 396 432 468 508 552 600 652
    704 768 832 896 960 1024 1088 1152 1216 1280 1344 1408 1472 1536

    sfb_offset[43] = 768 of 1536 lines

    CODED BANDWIDTH  0 .. 12000 Hz   (50 % of the spectrum)
    A-SPX            12000 .. 24000 Hz

12000 Hz is unchanged from the wrong table — sfb 43 is one of the rows where
the two coincidentally agree. The right answer for the wrong reason, and worth
flagging as such.

`max_sfb` is an EXCLUSIVE end: bands 0..42 are coded. M24 already depends on
this, decoding `SFB_OFFSET[min(max_sfb, 3)]` = 12 lines for the LFE.

**The fix's real payoff — all six channels (M28).** With the corrected offsets
the entire 5.1 channel element closes:

    2098 / 2098 long frames    all six channels, 0 overruns, 0 errors
    short frames               1150 (b_long_frame = 0), counted not hidden

    time structure       adjacent r / shuffled control
      LFE   +0.3157 / +0.0075      L   +0.6224 / -0.0200
      R     +0.8177 / -0.0132      Ls  +0.7535 / +0.0038
      Rs    +0.5275 / +0.0181      C   +0.7636 / +0.0147

    pair coupling        L/R  -0.5914  (control -0.0006)
                         Ls/Rs +0.1071 (control +0.0071)

L/R is strongly **negative** and that is the expected sign: with MDCT stereo
processing on (`sap_mode = 2`, M/S on all bands) the two `sf_data` blocks are
not left and right but sum and difference, and an encoder holding a bitrate
spends on one at the other's expense. My first gate demanded r > +0.3 and
failed a correct decode reading −0.44 — the same shape of error as M24's
">20 % zeros" gate. What proves the second block decoded correctly is a strong
relationship where the shuffled control has none, in *either* direction.

The per-channel gate thresholds were also wrong at first: a fixed
`|r_shuffled| < 0.06` demands less noise than a 249-frame sample can deliver
(its standard error is 0.063). The control is a measurement, not a constant,
so the test now compares against `3 × max(|control|, 1/sqrt(n))`.

**Field widths confirmed from the spec, not guessed:** Table 105 gives
n_msfb_bits = 6, n_side_bits = 5, n_msfbl_bits = 3 for transform length 1536 —
so the LFE's max_sfb is a 3-bit field and a full-band channel's is 6 bits.
Clause 4.3.6.2.4: n_grp_bits = 0 when b_long_frame is true and
frame_len_base >= 1536. Part 2's Table 72 lists `mono_data`, `two_channel_data`,
`sf_data`, `asf_section_data`, `chparam_info` and `companding_control` as
specified in part 1 and unchanged — checked, because reading part 1 syntax for
a bitstream_version-2 stream is a mistake already made once here.

**Element walk, as decoded (5.1, ASPX, coding_config 0):**

    5_X_codec_mode (3) = 1 ASPX
    mono_data(1)              LFE       max_sfb (3), then sf_data
    companding_control(5)
    coding_config (2) = 0
    2ch_mode (1)
    two_channel_data()        L / R     one sf_info + chparam_info, two sf_data
    two_channel_data()        Ls / Rs
    mono_data(0)              C         spec_frontend (1), sf_info, sf_data
    aspx_data_2ch x2, aspx_data_1ch     (not decoded; ~1730 bits remain)

`sf_data` is section → spectral → scalefac → **snf**, in that order. The LFE
decoder in M24 skipped `asf_snf_data`, which was harmless there because nothing
followed it in that substream and fatal here, because every later channel reads
from the bit position it leaves behind.


## M29 — PROGRAMME AUDIO (8/06)

The five full-band channels, dequantised, un-M/S'd, inverse-transformed and
overlap-added. This is the actual programme, not the LFE rumble.

    2098 / 3480 frames rendered (60.3 %)     the rest are short frames
    850 frames of alias-cancelled audio      27.2 s, written to tv_audio_valid.wav

    GATE 1   below 12000 Hz   99.96 %
             above 12000 Hz    0.01 %     PASS
    GATE 2   L/R correlation  +0.9918, mean|L-R|/mean|L| = 0.082   PASS

    spectrum   20-300 Hz 19.9 %   300-1000 Hz 48.1 %   1-3 kHz 21.1 %
               3-6 kHz 9.6 %      6-12 kHz 1.2 %       above 12 kHz 0.01 %

**Gate 1 is the one that cannot be fudged.** Only lines 0..767 were decoded, so
the output must be band-limited at exactly 12 kHz — a physical consequence of
what was decoded, checkable on the output alone. It fails loudly if the IMDCT
size, the line mapping or the overlap-add is wrong. It reads 0.01 %.

**An MDCT block alone is not a signal.** Each block carries time-domain
aliasing that cancels only when its neighbour is overlap-added. Where a short
frame was skipped, the neighbours' aliasing never cancels — those samples
reconstruct nothing. Measuring a run *including its edges* made the stereo gate
read r = −0.21 on audio that is actually +0.9918. Only run interiors are valid,
and only those are measured and written.

**Mid/side is not optional.** Table 113: `sap_mode = 2` = M/S in all bands, used
in 98 % of frames. L = M+S, R = M−S. Skipping it does not give slightly wrong
stereo, it gives the mid signal in one ear and the difference signal in the
other — which is why the two blocks' energies are anti-correlated (−0.59).

**Still missing, and counted rather than hidden:**

- **short frames, 40 %** — `b_long_frame = 0` needs two `transf_length` fields,
  up to two `max_sfb`, and up to 15 grouping bits (Table 108). This is now the
  single biggest gap: it is what makes the render choppy, and the longest
  unbroken run is only 15 frames.
- **A-SPX, 12–24 kHz** — the entire top half of the spectrum is parametric and
  not decoded. The audio is genuinely band-limited, not merely dull.
- **sap_mode 3 (full SAP), 31 frames** — passed through un-decorrelated.
- **companding** — signalled, not applied.
- the 4/3 quantiser exponent is still flagged READ-FROM-CONTEXT (M26).


## SHORT FRAMES — IMPLEMENTED, NOT PROVEN, OFF BY DEFAULT (`--short`)

The short-block branch is written and partially works. It is disabled by
default because turning it on makes the decoder *look* better and *be* worse,
which is worth recording precisely.

**What was built.** Table 99 (transf_length index → 96/192/384/768 at 1536@48),
Table 108 (n_grp_bits), Pseudocode 3 (num_windows, num_window_groups,
window_to_group), Pseudocode 4 (num_win_in_group and the packed
`sect_sfb_offset[g][sfb] = group_offset + sfb_offset[sfb] * num_win_in_group[g]`),
Pseudocode 5 (`get_max_sfb(g)`), and Annex B offset tables for all five
transform lengths — each extracted and gated against Table B.1's independent
band count.

Table 108 is transcribed verbatim, but it also *derives*: with block counts
n_i = 768 / length_i, n_grp_bits is `n0 + n1 - 1` when the halves share a length
and `(n0-1) + (n1-1)` when they differ. That reproduces all 16 entries, and
Pseudocode 3's `num_windows = n_grp_bits + 1` (plus one more for different
framing) is the same statement.

**Two real bugs found and fixed along the way:**

- *Stage nesting.* Each of the four `sf_data` stages loops over groups
  internally, so it is all sections, then all spectral, then all scalefac, then
  all snf — NOT one group end-to-end at a time. Measured: sections-first closes
  70.3 % of short frames, interleaved closes 1.7 %. A long frame has one group
  and parses identically either way, so this bug would have hidden indefinitely.
- *`n_msfb_bits` is indexed by transform length, not by frame.* Table 105 gives
  6 bits at 384 and above, 5 at 192, 4 at 96. Using 6 everywhere parses long
  frames and the 384/768 short frames and desynchronises exactly the rest —
  which is what the first run showed, with failures confined to the framings
  involving 96 or 192.

**Why it is still off.** With `--short`:

    frames closing   2471 / 3248  (76.1 %)   vs  2098 / 3248 (64.6 %) without
    L  adjacent r    +0.1220                 vs  +0.6224 without
    R  adjacent r    +0.1475                 vs  +0.8177 without
    Ls adjacent r    +0.0333                 vs  +0.7535 without
    L/R coupling     +0.0796                 vs  -0.5914 without

**More frames closed and the audio got worse.** The extra frames are closing
without being right, and they dilute every statistic. This is the sharpest
demonstration yet of something already written down here twice: a complete
prefix code decodes *anything*, so "the element closed" is necessary and never
sufficient. Had the closure rate been the only gate, this would have shipped as
an improvement.

Still unresolved in that branch: whether `get_transf_length(g)` in the
`n_sect_bits` test means the length or the index (both give 5 for our dominant
(3,3) case, so the stream cannot currently distinguish them), and the short
-block window shapes — AC-4 needs transition windows between long and short
blocks, and the overlap-add for partial blocks is not yet read from clause 5.
Rendering short frames without that would produce audio that plays and is
wrong, so M29 does not attempt it.


## 8/07 — THE FLAGGED CONSTANTS, CLOSED

**The 4/3 quantiser exponent is confirmed verbatim.** M26 carried it as READ
FROM CONTEXT because the PDF sets it as a stacked fraction whose glyphs do not
survive text extraction — clause 5.1.3.2 comes out as
`rec_spec=sign quant_spec x |quant_spec|`, exponent simply missing, and an
exponent you cannot see is an exponent you are guessing.

Settled by RENDERING the formula region of page 131 to a 500-dpi image and
reading it instead of trusting the text layer:

    rec_spec = sign(quant_spec) x |quant_spec|^(4/3)

Technique worth keeping: when a spec equation extracts as garbage, the
pseudocode restatement is the first place to look and a rendered crop of the
page is the second. Both beat inference.

Pseudocode 21 confirms the other two verbatim, so all three M26 constants are
now spec-sourced rather than assumed:

    scale_factor += dpcm_sf[g][sfb] - 60;                    -> SF_CENTRE = 60
    sf_gain[g][sfb] = pow(2.0, 0.25 * (scale_factor - 100)); -> SF_OFFSET = 100

**SNF_CENTRE was wrong: 17, not 11.** Pseudocode 23 states
`delta = dpcm_snf[g][sfb] - 17` directly. The original 11 was the midpoint of a
22-entry alphabet — right for a symmetric DPCM table, and `ASF_HCB_SNF` is not
symmetric: its shortest codewords sit at indices 13..16. The code-length profile
raised the suspicion; the spec settled it. Values only, not bit counts, so no
sync was ever affected.

**Spectral noise fill: implemented, and a STRUCTURAL ZERO on this stream.**
Clause 5.1.4 / Pseudocode 22-23 are implemented in `m29_audio.noise_fill()`:
`previous_rms` is log2 of the band mean-square (1.44269504 = 1/ln 2), the
transmitted delta adds in that log domain, and the inserted amplitude is
`2 ** (0.5 * noise_rms)`. Measured on the broadcast:

    b_snf_data_exists = 0 on every frame, every channel (L, R, C: 249/249)
    bands carrying an snf delta: 0

So this encoder never turns noise fill on, and enabling the tool changes the
output by nothing — the rendered gates are identical to four decimal places
with and without it. Kept anyway (correct for other streams), defaulted on,
disableable with `--no-noise-fill`.

This also retrospectively explains a flat experiment from the night before: the
snf structural variants (no snf / no flag / different conditions) all scored the
same ~2-5 % because they differed only in a code path that never executes on
this stream. A variant sweep over a dead branch measures nothing, and the sweep
gave no hint of that — the giveaway was available all along in the flag itself.

One honest deviation is recorded in the code: the spec's `GetRandomNoiseValue()`
is a specific table-driven generator seeded from `sequence_counter`
(Pseudocode 24) whose table lives in the speech front end, which is not
implemented. The spec states the output is normal with unit variance and zero
mean, so a normal RNG seeded from the same counter is used. Levels exact,
statistics as specified, realisation different — perceptually equivalent, NOT
bit-exact, and must not be called conformant.


## 8/07 — THE FILTERBANK: KBD, BLOCK SWITCHING, AND A WRONG WINDOW SHIPPED

Clause 5.5 (`IMDCT equations and block switching`) is now implemented and gated
in `m30_filterbank.py`.

**Table 185, alpha by transform length (44,1 / 48 kHz):**

    2048 1920 1536 -> 3     1024 960 768 -> 4     512 480 384 -> 4.5
     256  240  192 -> 5      128 120  96 -> 6

**The KBD kernel**, rendered from page 172 because the text layer mangles it
into `I pa1.0 -2/(N -1)`:

    W(N,n,a)      = I0(pi*a*sqrt(1.0 - (2n/N - 1)^2)) / I0(pi*a)   0 <= n < N
    KBD_LEFT(N,n) = sqrt( sum(W, p=0..n) / sum(W, p=0..N) )

Denominator is **N**, not N-1, and the normalising sum runs to p = N inclusive
— which is exactly what `m22.kbd_window` already computed, so that primitive
needed no change. (Second use of the render-the-page technique tonight, and the
second time it settled something inference would have got wrong or left open.)

**The transition window** (step 5): `NW = min(N, Nprev)`, `Nskip = (N-NW)/2`,
and `w` is zero for `Nskip`, KBD_LEFT(NW) for `NW`, then flat 1. The transition
region is always the SHORTER of the two blocks, centred, and it degenerates to
a plain KBD when `N == Nprev`.

**THE AMBIGUITY, AND HOW IT WAS SETTLED WITHOUT DERIVING ANYTHING.** Step 6
windows the previous block's stored half as
`for (n = 0; n < Nprev; n++) overlap[nskip_prev + n] *= w[n];` where `w` is the
CURRENT block's left window of length N. When `Nprev > N` that indexes past the
end of `w`, and with an ascending `w` it cannot satisfy TDAC, which needs the
previous block's right window to be the time-REVERSE of the next block's left.

So both readings were implemented and tested against a property that admits no
argument: **a lapped orthogonal transform satisfies `S @ S.T == I`.** The whole
synthesis chain is a linear map; it was built column by column by pushing unit
impulses through the real synthesis code, and the identity checked.

    uniform  [16]*6                    literal 5.62e-01   reversed 6.66e-16
    switch   [16,16,8,8,16,16]         literal 7.56e-01   reversed 4.44e-16
    switch   [16,8,8,4,4,4,4,16,16]    literal 8.27e-01   reversed 4.44e-16
    switch   [16,4,4,8,16,16]          literal 7.77e-01   reversed 4.44e-16

The `reversed` reading reconstructs EXACTLY across every switching sequence,
including a 4x jump and a nested 8->4->8. No encoder, no broadcast data, and no
hand-derivation trusted — one matrix identity settles it.

**The first run failed by a constant, and the constant was the answer.** Both
readings initially failed the uniform case at 8.75e-01 = 1 - 2/16, i.e.
`S @ S.T` came out as `(2/N) * I`: the construction was right and merely not
orthoNORMAL, because m22's IMDCT carries a 2/N. That is not cosmetic once
blocks switch — each block length would take a different gain, so short blocks
would sit at the wrong level against long ones and every transition would step
in loudness. Scaling by `sqrt(N/2)` makes every block unit-gain.

**A WRONG WINDOW WAS ALREADY IN THE SHIPPED AUDIO.** M29 rendered long frames
with a SINE window. AC-4 specifies KBD. The sine window also satisfies
Princen-Bradley, so it reconstructs *something* and **every gate in M29 kept
passing** — band-limiting at 12 kHz, L/R correlation +0.9918, the spectrum
shape. Measured against the corrected render:

    rms difference / signal   11.56 %
    window-change SNR         18.74 dB
    correlation               +0.9934

So the wrong window was injecting about -19 dB of uncancelled aliasing —
audible, roughly three bits of distortion, and invisible to every test that was
being run. Another instance of the standing lesson: gates that a *plausible*
wrong answer also passes are not gates. The band-limit test constrains WHICH
lines were decoded; nothing was testing the window against the encoder's.
`m29_audio.py` now uses KBD alpha=3 and the audio is regenerated.

**Still blocking short-frame rendering:** the filterbank is now proven for
block switching, so the remaining blocker is not the windows — it is the
short-frame PARSE, which still fails its correctness gates (see the section
above). Rendering an unproven parse through a proven filterbank would produce
audio that plays and is wrong.


## 8/07 — SHORT FRAMES FIXED: 3248 / 3248, AND CONTINUOUS AUDIO

The short-frame branch is on by default and every non-i-frame in the file now
decodes all six channels.

    GATE 1  3248 / 3248 decodable frames (100.0 %), 0 overruns, 0 errors
    GATE 2  lfe +0.3079  L +0.4316  R +0.7525
            Ls +0.7309   Rs +0.4995  C +0.5732     (controls all ~0)
    GATE 3  L/R -0.4352   Ls/Rs +0.1501            (controls ~0)

**BUG 1 — the wrong Annex B table, again.** `offsets_for()` picked a column by
searching for `col[num_sfb] == length`, on the reasoning that only the right
column could pass. It is not unique: the band tables are **nested**, so Table
B.4's 1536 column also has `col[43] = 768` — the same value at the same index
that identifies Table B.5's 768 column. B.4 is searched first, so every
768-sample partial block silently got the 1536 band layout (sfb 12 -> 52
instead of 56, sfb 20 -> 116 instead of 132).

This is the SAME SHAPE as the original Table B.4 contamination: a selection
rule that looked decisive, was under-constrained, and produced a coherently
wrong table. The fix maps each length to its table and column from the printed
headers (`TABLE_FOR`), and demotes the old rule from selector to **gate**.

Two independent signals had been pointing at it and I had read both as noise:

- *the bandwidth was inconsistent.* A-SPX has one crossover frequency, so the
  MDCT must code to the same Hz regardless of block size. Long frames gave
  12000 Hz; 768-sample blocks gave 8875 Hz. After the fix the 768 blocks read
  **exactly 12000 Hz**. (384/192/96 read 13000/13500/14000 — the nearest band
  edge at or above the crossover, which is what coarser resolution forces.)
- *the parse under-consumed.* Remainder after six channels, against a
  long-frame baseline of median 263 bits:

        before   pair1 SHORT            median  568, p95 2764, 6.5 % implausible
                 pair1+pair2 SHORT      median 1111, p95 2797, 17.4 %
        after    pair1 SHORT            median  271, p95  472, 0.3 %
                 pair1+pair2 SHORT      median  314, p95  673, 0.0 %

**BUG 2 — `chparam_info` loops per window group.** Both the `ms_used` loop
(sap_mode 1) and `sap_data` (sap_mode 3) run `for g` in the spec, `max_sfb_g`
is fetched per group, and `delta_code_time` exists **only when
num_window_groups != 1**. A long frame has exactly one group, so a
single-group implementation parses long frames perfectly and under-reads every
short frame using M/S-per-band or full SAP.

The failure fingerprint named it outright: of the 40 remaining errors, **33 were
sap_mode 3 and 5 were sap_mode 1** — while only 2 were sap_mode 2, the one mode
that reads no bits at all. Fixing it took 98.8 % to 100 %, and the 13 frames
that had been reported as "SSF (speech frontend)" vanished with it: they were
desyncs, and this stream never uses SSF.

**AUDIO IS NOW CONTINUOUS.** With block switching decoded and the M30
filterbank, there are no gaps and no splices:

    103.9 s continuous stereo
    5177 transform windows on L, 2813 of them PARTIAL blocks (54 %)
    GATE 1  100.00 % below 12 kHz, 0.00 % above
    GATE 2  L/R +0.9829, mean|L-R|/mean|L| = 0.105

Each channel carries its OWN framing and gets its own window sequence — the LFE
is always long, and the pairs and centre are signalled independently, so two
channels can have the same window COUNT and different lengths within a frame
((768,192,192,192,192) vs (192,192,192,192,768)). Forcing them onto one
sequence is what broke the first attempt.

**REMAINING GAP, stated plainly:** 232 of 3480 frames (6.7 %) are i-frames,
which carry `aspx_config()` after `5_X_codec_mode`, and they are still skipped.
That does not leave silence — it TIME-COMPRESSES the output, so the render is
103.9 s of 111.4 s of real time. Parsing `aspx_config` (Table 50) is the next
item, and it is needed for A-SPX anyway.


## 8/07 — I-FRAMES: 3480 / 3480, AND THE RENDER IS REAL-TIME EXACT

The last structural gap is closed. `aspx_config()` (Table 50) is **fifteen
fixed-width bits with no conditionals**:

    aspx_quant_mode_env 1   aspx_start_freq 3   aspx_stop_freq 2
    aspx_master_freq_scale 1   aspx_interpolation 1   aspx_preflat 1
    aspx_limiter 1   aspx_noise_sbg 2   aspx_num_env_bits_fixfix 1
    aspx_freq_res_mode 2

It sits immediately after `5_X_codec_mode` in an i-frame, and not parsing it is
why 232 of 3480 frames (6.7 %) were skipped. Reading those 15 bits:

    i-frames        232 / 232 decode, 0 failures
    bits remaining  median 320, p5 282, p95 434  (baseline was 263)

    WHOLE FILE      3480 / 3480 frames, 0 overruns, 0 errors
                    3479 adjacent pairs -- the stream is now fully contiguous
    time structure  lfe +0.3025  L +0.4333  R +0.7483
                    Ls +0.7311   Rs +0.4952  C +0.5764   (controls ~0)
    pair coupling   L/R -0.4282   Ls/Rs +0.1442          (controls ~0)

    AUDIO           111.4 s, continuous, REAL-TIME EXACT
                    100.00 % below 12 kHz, L/R +0.9829

Skipping i-frames had not been leaving silence, it had been TIME-COMPRESSING
the output — 103.9 s of 111.4 s of real content. That is now gone: the render
matches the capture second for second.

**A number recorded but NOT converted.** Every i-frame carries
`aspx_start_freq = 5`, `aspx_stop_freq = 1`, `aspx_master_freq_scale = 1`
(the high-resolution template). It is tempting to turn start_freq into a
crossover frequency and check it against the 12 kHz the MDCT implies -- clause
4.3.10.1.2 says an aspx_start_freq of 1 points to subband 20 and the index
moves in steps of 2 subbands, which would put 5 at subband 28. But the index
runs into the scale factor SUBBAND GROUP table (sbg_template_highres, clause
5.7.6.3.1.1), whose entries are not claimed to be uniformly spaced, and that
table has not been read. So the Hz value is left unstated rather than inferred
from an arithmetic that may not hold. It is A-SPX work, and A-SPX is next.


## 8/07 — THE A-SPX CROSSOVER, AND WHY THE NUMBER WAS WORTH NOT GUESSING

Last tick I recorded `aspx_start_freq = 5` and explicitly refused to convert it
to Hz, because clause 4.3.10.1.2's prose ("an index into the template tables
starting from the first QMF subband and moving upwards with steps of 2") gives
an arithmetic I could not verify. That refusal was correct: **the prose reading
is wrong.**

Clause 5.7.6.3.1.1 gives the static templates:

    sbg_template_highres = 18 19 20 21 22 23 24 26 28 30 32 34 36 38 40 42 44
                           47 50 53 56 59 62
    sbg_template_lowres  = 10 11 12 13 14 15 16 17 18 19 20 22 24 26 28 30 32
                           35 38 42 46

and Pseudocode 67 is unambiguous about what the "step of 2" indexes:

    sbg_master[sbg] = sbg_template_highres[2 * aspx_start_freq + sbg]

The step of 2 is in the **template index**, not in QMF subbands — and the
templates are NOT uniformly spaced above index 6. Reading it as "subband
18 + 2*start_freq" gives subband 28 = 10500 Hz for start_freq 5. The pseudocode
gives `template[10]` = subband **32 = 12000 Hz**. Both are plausible numbers,
one is right, and prose alone could not tell them apart.

The derivation reproduces all three of the spec's own worked examples:
start_freq 1 -> subband 20, stop_freq 2 -> subband 50, and lowres start 2 ->
subband 14 = 5250 Hz (Table 189).

**GATE 4 — a cross-check between subsystems that share no bits.**

    A-SPX  start_freq 5 -> QMF subband 32 = 12000 Hz   (i-frame header +
                                                        subband templates)
    A-SPX  stop_freq  1 -> QMF subband 56 = 21000 Hz
    MDCT   max_sfb 43   -> line 768/1536  = 12000 Hz   (Annex B band table)

A 64-band QMF at 48 kHz puts each subband at 24000/64 = 375 Hz, and 32 x 375 is
exactly 12000. **A-SPX begins precisely where the MDCT stops.** Nothing in the
decoder forces those two numbers to agree — one comes from Table 50 and the
subband group templates, the other from `max_sfb` and Table B.4 — so the
agreement is evidence about both.

It is also an independent confirmation of the short-block table fix: with the
wrong (nested B.4) table, 768-sample partial blocks coded to 8875 Hz, which
would NOT have met A-SPX's 12000 Hz start. With the correct B.5 column they
code to exactly 12000.

**And it bounds what is missing.** A-SPX spans 12000..21000 Hz. Above 21000 Hz
this stream reproduces nothing at all, so a complete decoder would be
band-limited at 21 kHz, not 24 kHz. The current renderer stops at 12 kHz; the
gap left to close is those 9 kHz, not 12.


## 8/07 — A-SPX, PART 1: THE FREQUENCY SKELETON (m31_aspx_bands.py)

A-SPX is the remaining half of the spectrum. It is a large tool, so it is being
taken in gated pieces rather than in one go. This piece is the frequency
skeleton: the subband group tables that `aspx_framing`, `aspx_hfgen_iwc_*` and
`aspx_ec_data` are all indexed by. Without these counts the A-SPX payload
cannot even be SKIPPED, let alone decoded.

**Everything derives from four header fields and one payload field.**
`aspx_xover_subband_offset` is the first element of `aspx_data_2ch` in an
i-frame, and since the six channel elements now close exactly, its bit position
is known rather than searched: read 3 bits at `r["bits"]`. Measured **0 in all
232 i-frames**.

    aspx_config: start_freq 5, stop_freq 1, master_freq_scale 1, noise_sbg 3
    aspx_xover_subband_offset: 0

    master  (10 groups)  32 34 36 38 40 42 44 47 50 53 56   12000..21000 Hz
    sig hi  (10 groups)  32 34 36 38 40 42 44 47 50 53 56   (xover 0, so equal)
    sig lo  ( 5 groups)  32 36 40 44 50 56
    noise   ( 2 groups)  32 40 56
    sbx 32 = 12000 Hz, num_sb_aspx 24 subbands = 9000 Hz wide

**num_aspx_timeslots = 12**, from Table 191 and Pseudocode 75a: a 1536-sample
frame is 1536/64 = 24 QMF timeslots, and frame_length 1536 has num_ts_in_ats 2,
so 24/2 = 12. That matters for parsing, not just bookkeeping: `aspx_framing`'s
Note 1 shrinks several relative-border fields from 2 bits to 1 when
num_aspx_timeslots <= 8. Ours is 12, so they stay at **2 bits** — a width that
would otherwise have had to be guessed, and guessed wrong 50 % of the time.

**Gates** (structural properties the spec states independently, since the
tables are pure arithmetic on stream values):

    PASS  num_sbg_noise = 2 <= 5            clause 5.7.6.3.1.3 states it as a shall
    PASS  every table strictly increasing
    PASS  signal range nested in master, noise nested in lowres
    PASS  master starts AND ends on an even QMF subband   (5.7.6.3.1.2 says so,
          and the arithmetic does not force it, so it tests the transcription)
    PASS  sbg_sig_lowres is a genuine subsequence of sbg_sig_highres
    PASS  num_sbg_sig_lowres follows the floor rule
    PASS  A-SPX starts at the MDCT coded edge (12000 Hz)

**WHAT REMAINS FOR A-SPX, stated so the scope is not understated.** This piece
is the skeleton only. The payload parse still needs `aspx_int_class` (a 1..3
bit VLC, Table 53), `aspx_delta_dir` (Table 54), and `aspx_huff_data`
(Table 58) with the 24 ASPX_HCB_* codebooks — those are already parsed and
Kraft-gated by M23, which helps. Beyond parsing, actually SYNTHESISING the high
band needs a 64-band complex QMF analysis/synthesis bank, HF patching
(Pseudocode 71, partially read), envelope adjustment, noise addition and the
additional-harmonics generator. That is a substantial project in its own right
and should not be described as nearly done: what exists today is the skeleton
plus a verified crossover, and the audio remains band-limited at 12 kHz.


## 8/07 — A-SPX, PART 2: THE PAYLOAD PARSES, AND audio_data CLOSES EXACTLY

`m32_aspx_parse.py` walks the A-SPX bits after the six channel elements:
`aspx_data_2ch()` (L/R), `aspx_data_2ch()` (Ls/Rs), `aspx_data_1ch()` (C).

    3021 / 3021 FIXFIX frames close EXACTLY to audio_size
    459 frames use a variable interval class: FIXVAR 225, VARFIX 224, VARVAR 10
        -- initially counted and excluded; now implemented, see below

**THE GATE WAS WRONG BEFORE THE PARSE WAS.** The first version demanded the
parse reach the end of the substream and failed by a stubborn ~42 bits. Table
16 explains it:

    ac4_substream() {
        audio_size = audio_size_value;  15
        if (b_more_bits) audio_size += variable_bits(7) << 15;
        byte_align;
        audio_data(channel_mode, b_iframe);
        fill_bits;  byte_align;  metadata(b_iframe);  byte_align;
    }

What follows `audio_data` is `fill_bits` and a `metadata()` payload — not slack.
So "the substream ends" was never the right target. `audio_size` is: the
encoder writes down how long audio_data is, so the parse must land within 0..7
bits below `start + audio_size*8`, the 0..7 being byte alignment. Measured
median -4, min -7, max 0, uniformly spread — exactly what alignment produces.

That is the strongest gate in the project so far. It is not "a piece decoded";
it is *every bit of the audio payload accounted for against a length the
encoder declared*, with no slack for a mis-sized field to hide in.

**Field notes that decided bit counts, all read rather than inferred:**

- `aspx_int_class` (Table 125) is a prefix code `0 / 10 / 110 / 111` =
  FIXFIX / FIXVAR / **VARFIX** / **VARVAR** — and note the last two are the
  REVERSE of the order Table 53 lists its switch cases in. Easy to transcribe
  backwards from the syntax table alone.
- `num_aspx_timeslots = 12 > 8`, so Table 53's relative-border fields stay at
  2 bits (Note 1). A coin-flip if guessed.
- `aspx_freq_res_mode = 2`, so `aspx_freq_res` is NOT in the bitstream —
  Pseudocode 77 derives it. For FIXFIX `aspx_tsg_ptr` is 0, so the first clause
  is vacuous and the test reduces to envelope duration >
  num_aspx_timeslots/6 + 3.25 = 5.25. With 12 timeslots, 1 envelope spans 12
  and 2 span 6, so both are HIGH resolution (10 subband groups, not 5). That
  choice alone changes the codeword count per envelope by a factor of two.
- `get_aspx_hcb` (Pseudocode 79) expands to
  `ASPX_HCB_ENV_<LEVEL|BALANCE>_<15|30>_<F0|DF|DT>` for signal and
  `ASPX_HCB_NOISE_<LEVEL|BALANCE>_<F0|DF|DT>` for noise, with quant_mode 0/1 to
  15/30. The 18 A-SPX codebooks were already Kraft-gated in M23.

**Still open for A-SPX:** the three variable interval classes (14 % of frames)
need the interval border construction of Pseudocode 76; and everything above is
still only PARSING — the decoded envelopes are not yet turned into audio, which
needs the 64-band complex QMF bank, HF patching, envelope adjustment, noise
addition and the harmonics generator. The rendered audio remains band-limited
at 12 kHz.


## 8/07 — A-SPX PARSE COMPLETE: 3480 / 3480

The three variable interval classes are implemented, via Pseudocode 76's border
construction feeding Pseudocode 77's duration test.

    3480 / 3480 frames close EXACTLY to audio_size
    (end of parse) - (start + audio_size*8):  median -4, min -7, max 0

**Every bit of audio_data in the entire capture is now accounted for** -- the
MDCT core (six channels, long and short frames, i-frames) plus the complete
A-SPX side information -- measured against the length the encoder itself wrote
into the substream header.

**Why the borders were needed at all.** It would be easy to assume the variable
classes only change WHERE envelopes sit, not how many bits they cost. They
change both: Pseudocode 77 picks each envelope's frequency resolution from its
DURATION, and that selects between num_sbg_sig_highres (10 groups) and
num_sbg_sig_lowres (5) for that envelope's Huffman run. Getting a border wrong
halves or doubles a codeword count, and the audio_size gate catches it
immediately.

**VARFIX and VARVAR are STATEFUL.** On a non-i-frame the leading border is
`previous_stop_pos - num_aspx_timeslots`, where previous_stop_pos is the
trailing border of the SAME channel group in the PREVIOUS frame. So the A-SPX
parse is sequential across the stream, not per-frame independent -- the first
structure in this decoder that cannot be evaluated on a single frame in
isolation. It is threaded per (channel group, channel).

**Bit widths that had to be read rather than assumed**, each of which would
have been a coin flip:

    aspx_num_rel_left/right   2 bits (num_aspx_timeslots 12 > 8), not 1
    aspx_rel_bord_*           the READ value maps as 2*tmp + 2, not raw
    aspx_tsg_ptr              ceil(log2(num_env+2)) bits, and the value is
                              tmp - 1, not tmp
    aspx_var_bord_*           2 bits, used RAW in the border arithmetic

**What A-SPX still does not do: make sound.** Everything above is parsing. The
decoded envelopes are not yet applied. Synthesis needs the 64-band complex QMF
analysis/synthesis bank, HF patching (Pseudocode 71), envelope adjustment,
noise addition and the additional-harmonics generator. The rendered audio
remains band-limited at 12 kHz, and A-SPX would extend it to 21 kHz.


## 8/07 — A-SPX, PART 3: THE 64-BAND COMPLEX QMF BANK (m33_qmf.py)

A-SPX operates entirely in the QMF domain, so none of the side information M32
parses can become audio until this filterbank exists. Clause 5.7.3 / 5.7.4,
implemented and gated.

    analysis + synthesis reconstruction   78.49 dB
    measured group delay                  577 samples (~9 QMF time slots)

**The window came from the C file, not the PDF.** Table D.3's 640 coefficients
are also shipped as `QWIN` in `ts_103190_tables.c`. M23's parser could not read
it because that parser is INTEGER-only -- it split `9.90318758627504e-04` into
three tokens and reported 1918 "values" for a 640-entry array. A float-aware
read gives exactly 640, matching the declared size.

**The spec contradicts itself on the synthesis modulation, and the pseudocode
wins.** Clause 5.7.4.2 step 2 gives the matrix as

    N[n][k] ~ exp(j*pi*(k+0.5)*(2n - 4*num_qmf_subbands - 1) / 128)   -> 257

while the pseudocode immediately below it prints

    exponent = j*(pi/128)*(sb+0.5)*(2*n - 255)                        -> 255

Both were built and the reconstruction decided: **255 gives 78.49 dB, 257 gives
43.41 dB.** Same resolution method as M30's step-6 ambiguity, and the same
lesson as the 4/3 exponent and the aspx_start_freq arithmetic -- when a spec
states something twice and they disagree, the pseudocode is the one to trust.

**TWO OF MY OWN MEASUREMENTS WERE WRONG BEFORE THE CODE WAS.**

- *The SNR was computed on misaligned segments.* x[i] corresponds to
  y[i + lag], but the gate compared `x[2048:]` against `y[lag:]` instead of
  `y[lag+2048:]` -- off by 2048 samples. It reported -2.98 dB and a best-fit
  gain of -0.004 on a bank that was reconstructing at 78 dB. The tell was that
  the SAME run's correlation search, correctly anchored at each signal's
  origin, read +1.000000: two measurements disagreeing that sharply means one
  of them is wrong, not that the thing under test is strange.
- *The window symmetry gate tested a property the window does not have.* 128
  coefficients are negative, confined to 64-sample blocks 1, 3, 6 and 8 -- the
  usual SBR convention of folding signs into the window. What actually holds is
  |w[n]| == |w[640-n]|, and that passes. Third time this project has failed a
  correct thing with a wrong expectation (after M24's ">20 % zeros" and M29's
  stereo-sign gate).

**Where A-SPX stands now.** Frequency skeleton (M31), full payload parse
closing to audio_size (M32), and a reconstructing QMF bank (M33). What remains
before the high band makes sound: HF patching (Pseudocode 71), envelope
adjustment, noise addition and the additional-harmonics generator. The rendered
audio is still band-limited at 12 kHz.


## 8/07 — A-SPX, PART 4: THE HF PATCH TABLE (Pseudocode 71)

A-SPX does not transmit the high band. It COPIES low subbands upward and then
reshapes them, and the patch table decides which source block feeds which part
of the A-SPX range. Implemented in `m31_aspx_bands.patches()`.

    HF patches (2): borders [32, 44, 56]
      patch 0: 12 subbands from 20..32 (7500..12000 Hz) -> 12000..16500 Hz
      patch 1: 12 subbands from 20..32 (7500..12000 Hz) -> 16500..21000 Hz

**The result corroborates itself.** Both patches source from subbands 20..32 --
which is 7500..12000 Hz, the top 4.5 kHz of the MDCT-coded band, ending EXACTLY
at the crossover. Nothing in the implementation imposes that the source block
end at sbx; it falls out of the spec's own arithmetic (the `odd` parity term,
the `sba - source_band_low + msb - odd` search bound, and the
`sbg_master[sbg] - sb < 3` escape). Getting any of those wrong moves the source
block off the crossover, so landing on it exactly is evidence the transcription
is right.

The two patches tile 12000..21000 Hz exactly, which is the A-SPX range M31
already derived independently from the subband group templates.

**Three new gates**, all passing:

    num_sbg_patches = 2 <= 5          clause 5.7.6.3.1.4 states it as a shall
    patches tile the A-SPX range exactly
    every patch sources from BELOW the crossover

**The A-SPX structural work is now complete**: frequency skeleton (M31),
payload parse closing to audio_size (M32), reconstructing QMF bank (M33), and
the patch table. What remains is signal processing rather than structure --
envelope adjustment, noise addition and the additional-harmonics generator --
plus applying the decoded envelopes M32 already reads. The rendered audio is
still band-limited at 12 kHz.


## 8/07 — A-SPX, PART 5: THE ENVELOPES ARE DECODED (m34_aspx_env.py)

M32 proved the A-SPX bit COUNT is right. That is not the same as having the
VALUES, and the envelopes are what actually shape the high band. Pseudocode 80
reconstructs them from the Huffman deltas:

    delta = 2 if (ch == 1 and aspx_balance == 1) else 1
    FREQ:  qscf[sbg][atsg] = delta * sum(aspx_data_sig[atsg][0..sbg])
    TIME:  qscf[sbg][atsg] = qscf_prev[sbg] + delta * aspx_data_sig[atsg][sbg]

with the TIME branch remapping sbg through sbg_idx_low2high / high2low whenever
the frequency resolution CHANGES between consecutive envelopes, so a resolution
switch does not corrupt the running value.

**Stateful across frames**, like the framing borders: `qscf_sig_sbg_prev` is the
last envelope of the PREVIOUS A-SPX interval, so a TIME-coded first envelope
depends on the frame before it. Second structure here that cannot be evaluated
on a single frame in isolation.

    3480 frames, envelopes reconstructed
    mean level  min -6.7  max 20.3  mean 3.8   (3.0 dB quantiser step)

    GATE 1  adjacent frames r = +0.7605, shuffled +0.0244 (se 0.017)   PASS
    GATE 2  mean slope -2.501 steps/subband-group                      PASS

Gate 1 is the same test that proved the MDCT scale factors: envelopes are
LEVELS, loudness is continuous in time, and random bits through a complete
prefix code produce plausible values but not time structure. Gate 2 is
independent of the arithmetic -- real programme material loses energy with
frequency, and -2.5 steps x 3.0 dB is -7.5 dB per subband group across the
A-SPX range.

**One inference flagged in the code.** Annex A lists only `codebook_length` for
the 18 A-SPX books -- no `cb_off`, unlike the ASF spectrum books. The delta
mapping is therefore taken from the code lengths: every DF and DT book has ODD
length with its single shortest codeword at index (n-1)/2, which is what a
zero-centred delta alphabet looks like. The F0 books are not centred, consistent
with an absolute level index, so they are used raw. The F0 mapping is inferred
from structure, not read from a clause; it offsets every envelope by a
constant, so it affects absolute high-band LEVEL, not shape and not the time
structure the gates test.

(Table A.15 also turned up free confirmation of two earlier inferences:
UNSIGNED_CB and CB_DIM for the ASF spectrum books match exactly what M24
derived from `cb_off == 0` back at the start.)

**What is left before the high band makes sound:** applying these envelopes as
gains to the patched QMF subbands, plus noise addition and the additional-
harmonics generator. Every input to that now exists and is gated -- patch table
(M31), QMF bank (M33), envelopes (M34).


## 8/07 — INTEGRATION: THREE SUBSYSTEMS AGREE ON THE CROSSOVER

The QMF bank is now pointed at our own decoded audio (M33 gate 4):

    subband 31  1.41e-04      subband 32  4.23e-06      subband 40  1.68e-08
    energy at or above subband 32 (12000 Hz): 0.0003 %

The decoded audio falls off a cliff at subband 31 and hits the numerical noise
floor from 33 up. So three subsystems that **share no code** all place the
crossover at QMF subband 32 = 12000 Hz:

    the Annex B band table    max_sfb 43 -> line 768 of 1536   (M27/M28)
    the A-SPX header          start_freq 5 -> subband 32       (M31)
    this QMF filterbank       measured band edge at subband 32 (M33)

Nothing in the code forces that agreement, so it is evidence about all three at
once -- and specifically it shows the QMF bank is correctly ALIGNED in
frequency, which reconstruction alone does not test: a bank with a mis-indexed
modulation would still reconstruct perfectly while putting the band edge in the
wrong place.

**Dequantisation read** (Pseudocode 82), ready for the gain stage:

    a = 2 if aspx_qmode_env == 0 (1.5 dB) else 1 (3.0 dB)
    scf_sig_sbg[sbg][atsg] = num_qmf_subbands * pow(2, qscf_sig_sbg / a)

so scf is an ENERGY scale factor -- each qscf step doubles power at 3 dB, which
is the consistency check on the exponent.

**HONEST SCOPE.** The remaining HF synthesis is envelope estimation
(Pseudocode 90), scale factor mapping (91), gain computation, the limiter,
noise addition and the additional-harmonics generator. That is a multi-stage
signal-processing chain and it has NOT been started. A partial implementation
would produce audio with energy above 12 kHz that plays and is wrong, which is
the failure mode this project has explicitly guarded against all night, so it
is being left until it can be done and gated whole. The rendered audio remains
band-limited at 12 kHz.


## 8/07 — A-SPX, PART 6: THE HF GENERATOR (m35_hfgen.py)

A-SPX reconstructs 12..21 kHz in four steps. **Only the first is implemented:**

    1. HF GENERATION        copy low QMF subbands into the A-SPX range   DONE
    2. envelope adjustment  scale them to the transmitted envelopes
    3. noise addition       fill where the envelope says noise, not tone
    4. additional harmonics insert sinusoids where flagged

    patch 0: sb 20..32 -> sb 32..44   (7500..12000 Hz -> 12000..16500 Hz)
    patch 1: sb 20..32 -> sb 44..56   (7500..12000 Hz -> 16500..21000 Hz)

    GATE 1  A-SPX range empty before the patch      0.00025 %      PASS
    GATE 2  destination energy == source energy     rel diff 0.0   PASS
    GATE 3  nothing above the A-SPX stop (sb 56)    0.000007 %     PASS

Gate 2 is EQUALITY, not similarity, because a copy is a copy -- the destination
energy matched the source to 0.0e+00 relative difference on both patches.

**THIS FILE DELIBERATELY WRITES NO AUDIO.** The patched band has the right
content at the WRONG LEVEL -- typically far too loud, because the entire point
of step 2 is that a copy of 7.5-12 kHz is not what 12-21 kHz should sound like.
Rendering it would produce a file that plays, sounds "fuller", and is wrong.
That is precisely the failure mode this project has guarded against since the
short-frame branch, so the stage is gated in the QMF domain -- where the claim
is exact and checkable -- and stops there.

For scale: the patched high band carries 0.34 % of the total energy as copied.
What the correct level is, only the envelopes decide.

**Remaining, with clause references so it is directly implementable:**
Pseudocode 90 (estimate the actual envelope of Q_high), 91 (map scf_sig_sbg and
scf_noise_sbg onto QMF subbands), 92 (sinusoid markers -- read), then the gain
computation, the limiter, and noise addition. Pseudocode 82 for dequantisation
is already read: `scf = num_qmf_subbands * pow(2, qscf/a)`, a = 2 at 1.5 dB and
1 at 3.0 dB.
