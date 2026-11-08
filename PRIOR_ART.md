# PRIOR ART — ATSC 3.0 receivers in open-source SDR, and what would finish ours

Research date: **2026-08-06**. No radio was opened. No code in this repo was changed.
Everything below is desk research plus direct inspection of published source code.

Convention used throughout:
**[FOUND]** = I read it directly (source file, spec PDF, or primary page).
**[CLAIM]** = a project's or vendor's own self-description, not independently verified.
**[INFER]** = my reasoning, not a citation.
**[NOT FOUND]** = searched, did not find; recorded rather than guessed.

---

## 0. The four findings that matter most

If you read nothing else:

1. **The ST sweep is unnecessary and, as currently parameterised, would have failed for
   all 37 values.** A/322 §7.1.5.4 does not contain a free integer `ST`. The "twisting
   parameter" is `T_i`, and the spec *defines* it as `T_i = R_i mod N_c`. The reading
   we settled on has the row/column enumeration transposed relative to the spec, so no
   value of ST can rescue it. Full equations, a worked test vector, and a numerical
   proof of the mismatch are in §1. **This is the single highest-value item in this
   document.**

2. **The "lost to the math font" problem is an extractor artifact, not a property of the
   PDF.** `pdftotext` drops every math-font glyph in A/322:2026-04. **PyMuPDF (`fitz`)
   reads the same file, same page, with every variable intact.** Demonstrated in §2 on
   our own current edition. This potentially reopens *every* assumption in DESIGN_NOTE.md
   that was closed with "unreadable in both editions and both extraction modes".

3. **The encryption picture is better than DESIGN_NOTE.md assumes.** Measured today:
   **34.4% of ATSC 3.0 subchannels are encrypted**, the figure has plateaued since
   December 2025, it is concentrated in NBC (76%) and CBS (68%), and **PBS encrypts
   nothing (0 of 37)**. Two thirds of services are legitimately decodable, and the
   signaling, lineup, PLP structure and emergency alerts are in the clear for *all* of
   them (§7).

4. **A complete, ATSC-V&V-verified ATSC 3.0 *transmitter* exists in GNU Radio and ships
   every table we have been extracting by hand** — LDPC, Annex B group-wise interleaver,
   NUC constellations, frequency-interleaver bit permutations, pilot tables, BCH,
   scrambler, bootstrap, cell interleaver, TBI, CDL. It is `drmpeg/gr-atsc3`. It also
   ships a 281 MB test transport stream that is **live and downloadable today**, which
   means we can synthesise a known-good ATSC 3.0 waveform end-to-end and finally have
   the referee we lack between the cell stream and the LDPC (§6).

---

## 1. THE ST QUESTION — ANSWERED, from the primary source

### 1.1 The equations, verbatim

**[FOUND]** ATSC A/322 §7.1.5.4 "Twisted Block Interleaver". Text below extracted with
PyMuPDF from four separate editions (2016, 2018, 2021, 2023-03) and from our own
**2026-04** edition; the wording and equations are **identical in all five**.

> In the Twisted Block Interleaver, the number of rows `N_r` shall be equal to the number
> of cells in a FEC Block while the number of columns `N_c` shall be set to
> `N_FEC_TI_MAX`. […] the number of virtual FEC Blocks in a TI Block is denoted by
> `N_FEC_TI_DUMMY(n,s) = N_FEC_TI_MAX − N_FEC_TI(n,s)`. Note that any virtual FEC Blocks
> that are included in a TI Block **shall be ahead of** data FEC Blocks in the same TI
> Block […]
>
> The FEC Blocks shall be **serially written column-wise** into the Twisted Block
> Interleaver memory […] Then, cells shall be read out diagonal-wise […] During the
> reading process, virtual cells belonging to virtual FEC Blocks shall be skipped.
>
> In a block interleaving array, the diagonal-wise reading can be performed by
> calculating the position for data and virtual cells with a coordinate `(R_i, C_i)`
> (for `i = 0, ⋯, N_r·N_c − 1`):
>
> ```
> R_i = i mod N_r
> T_i = R_i mod N_c
> C_i = ( T_i + floor( i / N_r ) ) mod N_c
> ```
>
> where `R_i` and `C_i` indicate the row and column indexes, respectively, and **`T_i` is
> a twisting parameter**. Assuming cells are read out sequentially from a linear memory
> array, the cell position can be calculated as `θ_i = N_r·C_i + R_i`. Note that virtual
> cells shall be skipped during the reading process if the condition of
> `θ_i ≥ N_FEC_TI_DUMMY(n,s) · N_r` is **not** satisfied.

Sources (all public, all direct PDF downloads):
- A/322:2026-04 — https://www.atsc.org/wp-content/uploads/2026/04/A322-2026-04-Physical-Layer-Protocol.pdf (§7.1.5.4, PDF page 81)
- A/322:2023-03 — https://www.atsc.org/wp-content/uploads/2023/09/A322-2023-03-Physical-Layer-Protocol.pdf (PDF page 77)
- A/322:2018 — https://www.atsc.org/wp-content/uploads/2016/10/A322-2018-Physical-Layer-Protocol.pdf (PDF page 77)
- A/322:2021 — http://www.atsc.org/wp-content/uploads/2021/04/A322-2021-Physical-Layer-Protocol.pdf (PDF page 77)
- A/322:2016 — https://www.atsc.org/wp-content/uploads/2016/10/A322-2016-Physical-Layer-Protocol.pdf (PDF page 77)

### 1.2 What this means for ASSUMPTION T1

**There is no free parameter.** DESIGN_NOTE.md M5 Step 3 records the twisting parameter
as "named once and never given a value anywhere the extraction can see", and narrows the
work to "one small integer over 0..36". That framing is an artifact of the broken
extraction: the sentence that names `T_i` as "a twisting parameter" survived `pdftotext`,
but the equation *defining* it (`T_i = R_i mod N_c`) did not. `T_i` is not a constant to
be found — it is a per-cell quantity fully determined by `R_i` and `N_c`.

In the vocabulary of DESIGN_NOTE.md, the effective twist is **ST = 1** (each successive
row is rotated by one additional column), but the more accurate statement is that the
twist is defined, not parameterised.

### 1.3 The reading in DESIGN_NOTE.md is transposed — and no ST rescues it

DESIGN_NOTE.md's surviving candidate is:

```
row = l // Ncols
col = (l mod Ncols + ST*row) mod Ncols
z   = (col - Nvirtual)*Nrows + row
```

The spec's enumeration is the **transpose** of this: `R_i = i mod N_r` (row advances
fastest, because `N_r` is the *cells-per-FEC-block* dimension = 2700), and the column
advance comes from `floor(i / N_r)`. Ours advances the column fastest.

**[FOUND — verified numerically, not asserted.]** For the A/327 worked-example geometry
`N_r = 4, N_c = 3`:

```
A/322 7.1.5.4 (spec equations):  [0, 5, 10, 3, 4, 9, 2, 7, 8, 1, 6, 11]

our reading, ST=0:               [0, 4,  8, 1, 5, 9, 2, 6, 10, 3, 7, 11]   != spec
our reading, ST=1:               [0, 4,  8, 5, 9, 1, 10, 2, 6, 3, 7, 11]   != spec
our reading, ST=2:               [0, 4,  8, 9, 1, 5, 6, 10, 2, 3, 7, 11]   != spec
```

Every one of ours is a valid permutation — which is exactly why gates T1 and T2 could
not eliminate it, and why T3 (recorded honestly in DESIGN_NOTE.md as "did not do what I
expected") had nothing to say. **A permutation test cannot distinguish a transpose.**

> **Therefore: the planned "sweep ST 0..36 against LDPC convergence" would have returned
> 37 failures, and the standing link-budget caveat would have made that result
> uninterpretable — "chain wrong" vs "SNR short" — exactly the ambiguity DESIGN_NOTE.md
> warned about. Do not run that sweep. Replace the TBI with the spec equations above.**

### 1.3b The correct reading is ALREADY in `m5_hti.py` — T4 selected the wrong one

**[FOUND — computed against our own code.]** `lab/m5_hti.py` enumerates
`READINGS = product(("div","mod"), ("mod","div"), ("colmajor","rowmajor"))`. The spec's
reading is inside that space. Searching all 8 readings × all ST for the one that
reproduces §1.1's permutation exactly:

| geometry | reading(s) matching the spec |
|---|---|
| `Nr=4, Nc=3, Nvirt=1` | `('mod','mod','colmajor') ST=0` **and** `('mod','div','colmajor') ST=1` |
| `Nr=4, Nc=3, Nvirt=0` | same two |
| `Nr=5, Nc=7, Nvirt=3` | **`('mod','div','colmajor') ST=1`** only |
| `Nr=2700, Nc=37, Nvirt=0` | **`('mod','div','colmajor') ST=1`** only |
| `Nr=2700, Nc=37, Nvirt=2` | **`('mod','div','colmajor') ST=1`** only |

Unique answer:

```python
reading = ("mod", "div", "colmajor")     # NOT ("div", "mod", "colmajor")
st      = 1                              # fixed by the spec, not swept
```

(The `('mod','mod',…) ST=0` twin appears only at 4×3, where `l % nrows` and `l % ncols`
degenerate because 4 and 3 are both small — it dies immediately at 5×7 and at our real
geometry. Worth noting as a reminder that a single small worked example is not enough to
pin a reading; that is why §1.5's vector should be used *together* with a second
geometry.)

`m5_hti.py` currently hardcodes `TEXT_CONSISTENT = ("div", "mod", "colmajor")` — the
first element differs. **The fix is one line, plus deleting the sweep.**

### 1.4 Where our T4 gate went wrong

DESIGN_NOTE.md gate **T4** reads: "the prose fixes enumeration order and the linear-array
convention." The prose is the trap. A/322 §7.1.5.4 says cells are read

> "diagonal-wise **from the first row (rightwards along the row beginning with the
> left-most column)** to the last row"

which reads as row-major traversal — and that is the reading we adopted. But the
**equations directly beneath it** advance the row fastest, not the column. The prose and
the equations disagree, and **the equations are what the reference implementation and the
worked example follow** (§1.5, §1.6). Worth adding to the bug ledger: *when prose and
equations in the same subsection disagree, the equations win, and a gate built on prose is
not a gate.*

### 1.5 An independent published test vector — and our equations reproduce it exactly

**[FOUND]** ATSC A/327:2025-06 "Guidelines for the Physical Layer Protocol", §6.3.2.2.2
"Twisted Block De-Interleaver", Figure 6.5(a), gives a worked example with
`N_r × N_c = 4 × 3`, `N_FEC_TI_MAX = 3`, `N_FEC_TI = 2` (i.e. **one virtual FEC block**,
exercising the skip logic). The published output cell order after TBI is:

```
b0  g0  a0  f0  d0  e0  c0  h0
```

Implementing §1.1's equations verbatim — including the
`θ_i ≥ N_FEC_TI_DUMMY · N_r` skip — and labelling the two data FEC blocks `a..d` and
`e..h` reproduces:

```
b0  g0  a0  f0  d0  e0  c0  h0        <- 8/8 exact
```

**This is a gold test vector for the TBI, from ATSC's own documentation, that requires no
radio, no FEC and no capture.** Build the de-interleaver against it before touching the
air. Source: https://www.atsc.org/wp-content/uploads/2025/07/A327-2025-06-Physical-Layer-RP.pdf
(PDF pages 89–92).

### 1.5b Sanity-checked at our real geometry

**[FOUND — computed]** The corrected equations at RF33 PLP 0's actual geometry
`N_r = 2700, N_c = N_FEC_TI_MAX = 37`:

```
N_FEC_TI = 37 (N_virtual = 0):  emits 99,900 cells (= 2700 x 37), bijection = True
N_FEC_TI = 35 (N_virtual = 2):  emits 94,500 cells (= 2700 x 35), bijection = True
consecutive output cells drawn from the same FEC Block: 0 in both cases
```

Note the last line: the corrected reading satisfies gate **T3** — and so did ours. That
is a second confirmation that T3 carried no discriminating information here, exactly as
DESIGN_NOTE.md recorded when it refused to claim T3 had done more than it did.

### 1.6 Confirmed a third time, by an implementation that passed ATSC's V&V suite

**[FOUND]** `drmpeg/gr-atsc3`, `lib/framemapper_cc_impl.cc` (and identically in
`subframemapper_cc_impl.cc`, `tdmframemapper_cc_impl.cc`, `fdmframemapper_cc_impl.cc`,
`ldmframemapper_cc_impl.cc`):

```c
for (int n = 0; n < fec_cells * Nfec_ti_max; n++) {
  Ri = n % fec_cells;
  Ti = Ri % Nfec_ti_max;
  Ci = (Ti + (n / fec_cells)) % Nfec_ti_max;
  if ((fec_cells * Ci) + Ri >= (Nfec_ti_max - HtimeNfec[x]) * fec_cells) {
    Htime[q++] = (fec_cells * Ci) + Ri;
  }
}
```

Variable-for-variable identical to §1.1, with `fec_cells = N_r` and
`Nfec_ti_max = N_c`. The project states its baseline "has been verified against the ATSC
3.0 validation and verification suite" **[CLAIM — the verification itself is not public,
but the claim is made by the author in the repo README]**.

Note in particular the **virtual-block convention**: gr-atsc3 skips source addresses below
`(N_c − N_FEC_TI)·N_r` and then subtracts that offset when reading the compacted input
stream (`in[H[n] - virtual_offset]`). That matches DESIGN_NOTE.md's `(col - Nvirtual)`
handling, so that part of our reading was right.

### 1.7 A fourth, independent, plain-language confirmation

**[FOUND]** Sungho Jeon et al., *"Physical Layer Time Interleaving for the ATSC 3.0
System"*, IEEE Transactions on Broadcasting (2016) — describes the TBI as one that
"**rotates each row to the left** and reads it out **column by column**", which is
precisely a per-row twist of one column with a column-major read.
https://ieeexplore.ieee.org/document/7378942/ (open PDF mirror:
https://scispace.com/pdf/physical-layer-time-interleaving-for-the-atsc-3-0-system-1i0adahoh7.pdf)

### 1.8 Bonus: the cell interleaver, and a printed test vector for it

RF33 PLP 0 has `cell_interleaver = 0` so this is bypassed there — but PLP 1 may not be,
and A/322 §7.1.5.2 turns out to carry its **own printed test vector** that extracts
cleanly with PyMuPDF:

> "under the condition of `N_cells = 10800` and `N_d = 14`, the shift value `P(r)` to be
> added to the basic permutation (for `r = 0,1,2,3, etc.`) would be
> **0, 8192, 4096, 2048, 10240, 6144, 1024, 9216**, etc."

That is a bit-reversal-based shift sequence, and gr-atsc3's `HtimePr` loop generates
exactly it. Free gate, no air time. **[FOUND]** (A/322:2016 PDF page 76; same text in
current editions.)

---

## 2. THE EXTRACTION UNLOCK — `pdftotext` is the bug, not the PDF

**[FOUND — demonstrated on our own current edition.]**

`pdftotext` on **A/322:2026-04** (the edition in `lab/`), §7.1.5.4:

```
 =  mod ,
 =  mod ,
 =  +   mod ,
where  and  indicate the row and column indexes, respectively, and  is a twisting parameter.
```

That is character-for-character the failure recorded in DESIGN_NOTE.md ("`= mod ,` /
`= mod ,` / `= + mod ,` — every variable lost to the math font, in BOTH editions and BOTH
extraction modes"). `pdftotext -layout` fails identically.

**PyMuPDF (`fitz`) on the same file, same page:**

```
𝑅𝑅𝑖𝑖= 𝑖𝑖 mod 𝑁𝑁𝑟𝑟,
𝑇𝑇𝑖𝑖= 𝑅𝑅𝑖𝑖 mod 𝑁𝑁𝑐𝑐,
𝐶𝐶𝑖𝑖= ቀ𝑇𝑇𝑖𝑖+ ቔ𝑖𝑖/𝑁𝑁𝑟𝑟ቕ ቁmod 𝑁𝑁𝑐𝑐,
```

Every variable present. (The glyphs are doubled — a known artifact of the embedded
subset font — and stacked fractions arrive as PUA glyphs `ቔ ቕ ቀ ቁ`, so a small
normaliser is needed: collapse doubled characters, map the bracket PUA codepoints to
`floor(...)`. Both are mechanical.)

**Caveats, stated honestly:**
- **[FOUND]** the doubling and PUA artifacts are real; extraction is not clean text, it
  is *recoverable* text. A normaliser must be written and gated.
- **[INFER]** the 2021 edition shows a different subscript corruption (`N_i`/`N_l` in
  place of `N_r`/`N_c`), so cross-edition diffing is still worth doing — but now the
  diff has content on both sides.

**Consequences worth chasing immediately:**
- **Assumption T1 dissolves** (§1).
- **Assumption B2** ("Tables 6.12/6.13 carry no extractable content — check-mark
  glyphs") should be **re-tested** with PyMuPDF before staying closed.
- Any other place DESIGN_NOTE.md says an equation was unreadable is now a re-run, not a
  reverse-engineering job. Given the standing law that a recorded failure is a reference
  and not an authority, this is a *cause-class* reopen: the cause was the extractor, and
  the re-test trigger is "re-extract with PyMuPDF".

Also worth noting: **eight editions of A/322 are freely downloadable** (2016, 2017, 2017a,
2018, 2020, 2021, 2023-03, 2026-04), so the "two witnesses" gate can become an
eight-witness gate at zero cost.

**[FOUND — retrieved and diffed]** **A/322:2026-06 exists and is newer than our 2026-04
copy**: https://www.atsc.org/wp-content/uploads/2026/06/A322-2026-06-Physical-Layer-Protocol.pdf
(7,936,979 bytes). §7.1.5.4 is **byte-identical in content** to 2026-04 (same page 81,
same equations), so the TBI conclusion is unaffected — but the repo should pick up the
current edition, and the rest of the document has not been diffed.

---

## 3. Prior art census

### 3.1 Findings table

| Project | What it actually does | Lang | License | Maintained | Covers our gaps? |
|---|---|---|---|---|---|
| [drmpeg/gr-atsc3](https://github.com/drmpeg/gr-atsc3) | **Full ATSC 3.0 transmitter** PHY (A/321+A/322): bootstrap, LDPC, BCH, bit interleaver, NUC modulator, cell interleaver, **TBI**, CDL, pilots, freq interleaver, framing, ALP BB header, PAPR. TX only. | C++ | GPL-3.0 | **Yes** — last commit 2026-06-18; 51★ | **HTI/TBI: YES (definitive).** ALP: header only. ROUTE/DASH: no. |
| [nextgenbroadcast / libatsc3](https://github.com/kansonkong/libatsc3) (also [junhuac/libatsc3](https://github.com/junhuac/libatsc3)) | A/331 and above: LLS/SLT parsing, ALP, ROUTE/ALC object recovery, MMT, NRT; Android/AppleTV sample apps. **Starts at IP/ALP, no PHY.** | C (+Java/ObjC) | MIT | **Marginal** — main repos last pushed 2023-07 / 2019 | **ALP: YES. ROUTE/DASH: YES.** HTI: no. |
| [lukefay/ATSC3_Rx](https://github.com/lukefay/ATSC3_Rx) | Receiver **above IP**: finds SLT at 224.0.23.60:4937, runs ROUTE to collect DASH segments, serves to a browser player; handles AST alignment. Derived from Stockhammer's ATSC_ROUTE. | Python + PHP | unspecified | Active-ish — 2026-06 | **ROUTE/DASH/SLT: YES.** PHY: no. |
| [haudiobe/ATSC_ROUTE](https://github.com/haudiobe/ATSC_ROUTE) | Reference ROUTE receiver + Apache/PHP playback setup (Thomas Stockhammer). | C | unspecified | **No** — 2020-05 | ROUTE: yes (dated). |
| [awanga/usrp-atsc3](https://github.com/awanga/usrp-atsc3) | *Self-described* "full ATSC 3.0 physical-layer receiver". **In reality a scaffold** — see §3.2. | C++ | none | Commits 2026; 0★ | **No.** |
| [zrcheng1991/atsc-3.0](https://github.com/zrcheng1991/atsc-3.0) | ATSC 3.0 **LDPC only**: parity-check matrix generation (Type A/B), encoder, sum-product decoder. | MATLAB | none | No — 2025-08 | Useful LDPC cross-check only. |
| [philburr/atsc3](https://github.com/philburr/atsc3) | "ATSC 3.0 Encoder" | — | GPL-3.0 | **No** — 2020-05, 6★ | Unassessed; superseded by gr-atsc3. |
| [Silicondust/atsc3-sqa](https://github.com/Silicondust/atsc3-sqa) | Generates **TS-over-ALP test data as PCAP**. Hosted sample PCAPs are now 404, but the generator source is present. | C | MPL-2.0 | No — 2020-04 | **ALP test vectors: YES (by generation).** |
| [LiamRPower/wireshark-atsc3](https://github.com/LiamRPower/wireshark-atsc3) | Wireshark build with **STLTP / LLS / ROUTE / MMT dissectors**. | C | Wireshark (GPL-2.0) | No — 2023-03 | Debugging tool for M5/M6 — high practical value. |
| [junhuac/atsc-3.0-mmt-pcaps](https://github.com/junhuac/atsc-3.0-mmt-pcaps) | **Real in-market ATSC 3.0 IP captures** (ROUTE-DASH and MMT), CES 2019 + Hunt Valley. Git-LFS. | data | none | No — 2019 | **IP-layer test data: YES, and live (§6.3).** |
| [nekohkr/danttoUHD](https://github.com/nekohkr/danttoUHD) | ATSC 3.0 → MPEG-2 TS converter for a **Korean LG UHD STB container**. Includes decryption functionality. | C++ | MIT | Yes — 2025-09 | **Out of scope — do not use** (§7.4). |
| [johnb-7/hdhr-ac4](https://github.com/johnb-7/hdhr-ac4) | Emulates an HDHomeRun tuner, proxies a real HDHR5-4K, **transcodes AC-4 → AC-3**. | Python | Apache-2.0 | No — 2024-05, 98★ | Playback-layer lesson (§8). |
| [ferrellsl/VideoPlayer](https://github.com/ferrellsl/VideoPlayer) | Qt/SDL2 player using an **experimental AC-4 ffmpeg branch**. | C | none | No — 2023 | Playback-layer lesson (§8). |
| [yanglikkk/ATSC3.0Parser](https://github.com/yanglikkk/ATSC3.0Parser) | ATSC 3.0 stream parser (Chinese). | Python | none | No — 2018 | Unassessed, likely signaling-only. |
| Dionisio & Akamine, *"ATSC 3.0 implementation in GNU Radio Companion"*, IEEE 2017 | Academic ATSC 3.0 **modulator** in GRC, validated against a DekTec DTA-2131 + Atsc3Xpress receiver. TX only. | GRC | paper | n/a | No. |
| GNU Radio in-tree `gr-atsc` / `gr-atsc3` (core) | **ATSC 1.0 only.** Core GNU Radio cannot transmit or receive ATSC 3.0. | C++ | GPL-3.0 | Yes (1.0) | No. |
| SDRangel / SDR++ | **[NOT FOUND]** — no ATSC 3.0 demodulator plugin located in either. SDRangel's DATV demod covers DVB-S/S2 and analog ATV. | — | — | — | No. |

### 3.2 `awanga/usrp-atsc3` — assessed and dismissed, with evidence

This is the only repository on GitHub whose description claims *"a full ATSC 3.0
physical-layer receiver"*. It does not deliver one, and it is worth recording why so the
claim does not get recycled.

**[FOUND — read from the repository's own source]:**
- `lib/fec/ldpc_decoder.cc` contains, in the function that is supposed to load the A/322
  parity-check matrices:
  ```c
  // TODO: Load actual ATSC 3.0 parity matrices from JSON files
  // For now, generate a simplified test matrix
  generate_test_matrix();
  ```
  and inside that: *"This is NOT the actual ATSC 3.0 matrix, but allows testing the
  algorithm"*. There is no A/322 LDPC in this project.
- `lib/ofdm/time_deinterleaver.cc` implements **only** the convolutional (CTI) delay-line
  mode. There is no HTI, no TBI, no twist. Its header cites "ATSC A/322 Section 8.2 (Time
  Interleaving)" — time interleaving is §7.1.5; §8 is the waveform/pilot section. The
  LDPC file likewise cites "Section 9 (LDPC Coding)"; §9 is L1 signaling.
- `test/captures/*.sigmf-data` are **135-byte Git-LFS pointers**, not captures.
- The repo contains a `CLAUDE.md` and a `TASKS.md`, and the entire body of work lands in
  a two-week burst in May 2026. 0 stars, no issues, no external contributors.

**[INFER]** This is an AI-scaffolded project with correct-looking architecture and no
verified decode. It is not evidence that the problem has been solved, and nothing in it
is worth reusing. Recording it here so that a future search that turns up its README does
not mistake the description for a result.

### 3.3 The shape of the census

**[FOUND]** Ron Economos (`drmpeg`) is the most prolific author of SDR digital-television
implementations in the GNU Radio ecosystem — gr-dvbs2 (163★), gr-paint, gr-dvbt2,
gr-dvbs, gr-dvbc, and notably **gr-dvbs2rx, "DVB-S2 and DVB-T2 *receiver* blocks"**. So
this author demonstrably builds receivers when he chooses to. For ATSC 3.0 he built a
**transmitter and stopped**. **[INFER]** That asymmetry is the strongest single signal
about the difficulty and the payoff: the receive side is where the work is, and the
encryption picture (§7) blunts the reward.

---

## 4. What in gr-atsc3 is directly reusable — and the licensing caveat

**[FOUND]** by direct inspection of the cloned repository (29,763 lines across `lib/`):

| We extracted by hand | gr-atsc3 has it as | Reuse value |
|---|---|---|
| A/322 Annex A LDPC parity tables | `lib/ldpc_bb_impl.h` — `ldpc_tab_2_15N[29][21]` … all rates, N and S | Independent second reading of every table |
| Annex B group-wise permutations (our 120) | `lib/interleaver_bb_impl.h` — `group_tab_<rate>_<mod>[180]` per rate/constellation | Direct diff against `spec_bitint.py` |
| Annex C NUC constellations (our 36) | `lib/modulator_bc_impl.h` — `mod_table_*QAM`, up to `4096QAM[12][32]` | Direct diff against `spec_nuc.py` |
| Frequency de-interleaver | `lib/freqinterleaver_cc_impl.h` — `bitperm8keven/odd`, `16keven/odd`, `bitperm32k` | Confirms our FI reading |
| Pilot/cell-map tables | `lib/pilotgenerator_cc_impl.cc` (2156 lines), `scattered_power_table[16][5]` | Cross-check on the 2420-identity gate |
| Baseband scrambler | `lib/bbscrambler_bb_impl.cc` | Confirms the Galois-LFSR finding |
| BCH | `lib/bch_bb_impl.cc` | — |
| Bootstrap (incl. A/321) | `lib/bootstrap_cc_impl.cc` (1539 lines) | Confirms M1 |
| **HTI: cell interleaver + TBI + CDL** | `lib/framemapper_cc_impl.cc` ~line 1810–1945 | **The M5 blocker** |
| ALP base-band header | `lib/alpbbheader_bb_impl.cc` | Partial M4→M5 |

**Why this matters epistemically.** DESIGN_NOTE.md is explicit that two editions of the
same PDF agreeing "only proves the extraction is faithful, not that the READING is
right". gr-atsc3 is a **genuinely independent reading** — a different person, from the
same public spec, whose output was checked against ATSC's own validation and verification
suite. Diffing our extracted tables against it is a materially stronger gate than any
second-witness PDF, and it costs an afternoon.

**Licensing caveat — [INFER], not legal advice.** gr-atsc3 is **GPL-3.0**. Numeric tables
transcribed from a published standard are a weak candidate for copyright in themselves,
but *copying them out of GPL source* is the kind of thing that creates an obligation and
an argument. The clean move, which is also the epistemically better one:

> **Use gr-atsc3 as a comparison oracle, not as a source. Keep extracting from the spec
> PDF (now with PyMuPDF), then diff against gr-atsc3 and record the diff count as a gate.
> Never paste.**

This preserves the clean-room character of `` and turns a licensing hazard
into an additional referee.

---

## 5. Anything that covers ALP → IP → ROUTE/DASH?

Yes, and this is the second-largest saving after §1.

- **`libatsc3`** **[FOUND, repo inspected; functionality is CLAIM from README]** — C
  library, MIT, covering exactly the layers M5/M6 needs: LLS/SLT parsing, ALP, ROUTE/ALC
  object recovery (including EFDT), MMT, NRT delivery, plus the Low-Level Signaling
  multicast conventions. Originated from the NGBP / ngbp.org effort associated with
  Sinclair-adjacent ATSC 3.0 open-source work (from the NGBP / ngbp.org effort).
  **Maintenance is the weak point**: the most active mirror (`kansonkong/libatsc3`) last
  pushed 2023-07; `junhuac/libatsc3` 2019. MIT licence means we may vendor or port freely.
- **`lukefay/ATSC3_Rx`** **[FOUND]** — the closest thing to our M5+M6 written in Python,
  and it makes exactly the assumptions our chain will produce: it listens for the SLT on
  `224.0.23.60:4937`, collects ROUTE objects into DASH segments, fixes up the DASH
  `availabilityStartTime` against segment arrival, and hands to a browser MSE player. It
  also documents the practical codec-support checks
  (`MediaSource.isTypeSupported('audio/mp4; codecs=ac-4.02.01.00')`) that turn out to
  matter a lot (§8).
- **`wireshark-atsc3`** **[FOUND]** — once we emit IP packets, this gives instant
  visibility into LLS/ROUTE/MMT without writing a parser first. Highest
  debugging-value-per-hour item in the transport milestones.

**[INFER]** Taken together, M5 and M6 are much less "unbuilt work" than DESIGN_NOTE.md's
ladder implies. The honest re-estimate is that once ALP packets exist, there is a
well-trodden path with reference code at every step, and the remaining risk is
operational (capture length, timing) rather than cryptographic or reverse-engineering.

---

## 6. Test vectors, reference streams, and the referee we lack

### 6.1 Does ATSC publish bit-exact A/322 test vectors publicly?

**[NOT FOUND] for bit-exact PHY vectors.** The ATSC 3.0 **Validation & Verification
(V&V)** effort is real and documented — see Verification and validation of the physical
layer ATSC 3.0 standard, IEEE BMSB 2016
(https://ieeexplore.ieee.org/document/7521999/) which describes software and hardware
tests across plug-fests in Shanghai (Oct 2015) and Baltimore (Mar 2016) — but the V&V
suite itself appears to be an ATSC member/participant artifact and I could not locate a
public download. Conformance test suites exist commercially (Eurofins/StreamWise, DekTec)
**[CLAIM — vendor descriptions]**.

**What IS public and useful:**

- **A/327 worked examples** — the TBI/TBDI example in §1.5 is a genuine, citable,
  bit-level test vector. **[FOUND]** A/327 also contains Annex C "Physical Layer
  Parameters" reference configurations (single-PLP, 2-PLP TDM, 2-PLP LDM, 3-PLP TDM,
  4-PLP TDM, 2-subframe) — i.e. published, self-consistent parameter sets to test a
  signalling parser against.
- **A/322's own printed vectors** — the cell-interleaver `P(r)` sequence (§1.8), the
  A/321 Annex B Gray table already used in M1, and the scrambler byte already used in M2.
  Now that PyMuPDF works, sweep the spec for every other printed "for example …" value.
- **A/325 PHY Lab Performance Test Plan** — **[FOUND]**, public, but its C/N tables are
  *blank templates* for a lab to fill in, not reference data.
  https://www.atsc.org/wp-content/uploads/2023/04/A325-2023-04-Lab-Performance-Test-Plan.pdf

### 6.2 The referee we actually lack — and how to build it this week

gr-atsc3 ships **V&V test flowgraphs** `examples/vv031.grc`, `vv320.grc`, `vv410.grc`,
`vv503.grc`, `vv504.grc` (named after V&V test cases), and its README points at a test
transport stream.

**[FOUND — verified live, 2026-08-06]:**
```
https://www.w6rz.net/advatsc3.ts    HTTP 200    281,912,392 bytes
```

`vv031.grc` is a fully specified configuration: 8K FFT, GI_5_1024, PILOT_SP3_4, SPB_4,
256QAM, rate 9/15, FECFRAME_NORMAL, 2 preamble syms + 72 payload syms, CRED_0, MISO off,
PAPR off, `timode: TI_MODE_OFF`, L1 FEC mode 1 — matched to ~15 dB S/N at 21.34 Mbps.

**This closes the gap DESIGN_NOTE.md identifies as the whole problem with M2/M4:** *"a
self-built encoder tests your encoder against your decoder, not against the spec."*
gr-atsc3 is **not** our encoder. Running it produces a waveform from a
V&V-verified independent implementation, with known input bits. That gives us, for the
first time, a stage-by-stage referee between the cell stream and the LDPC:

1. Build gr-atsc3, feed `advatsc3.ts`, dump complex samples to a file (no radio, no USRP
   — replace the `uhd_usrp_sink` with a file sink).
2. Run our chain on that file. Every intermediate is now checkable against a known truth,
   at arbitrary SNR, with the geometry we choose.
3. **Crucially: configure a run with `TI_MODE_HYBRID` and a non-trivial
   `tiblocks`/`tifecblocks`** so the TBI is actually exercised — `vv031` has
   `TI_MODE_OFF`, so it will *not* test the piece we care about most. This is a
   one-parameter change in the flowgraph.

**[INFER]** This is probably worth more than any other single item in this document apart
from §1, because it converts "a total failure is ambiguous between chain-wrong and
SNR-short" — the exact trap flagged as blocking the M5 gate — into a decisive test.

### 6.3 Real off-air IP-layer data, verified downloadable

**[FOUND — verified live, 2026-08-06]** `junhuac/atsc-3.0-mmt-pcaps` stores real
in-market captures via Git-LFS, and the LFS objects are still served. Probing one:

```
2019-01-09-LAS-587mhz-ROUTE-DASH/00.43.46.pcap
  oid  sha256:b4b5eb28…  size 781,448,521 bytes   -> download URL issued OK
```

Multiple 15-minute ROUTE-DASH captures from Las Vegas (CES 2019, 587 MHz) plus MMT
captures from Hunt Valley. **[INFER]** These are 2019-vintage and predate widescale
encryption, which makes them *more* useful as clear-content development data for the
ROUTE/DASH/SLT stack than anything we could capture off air today.

Also: `Silicondust/atsc3-sqa` (MPL-2.0) generates TS-over-ALP PCAPs locally. Its hosted
sample files are **404 as of today [FOUND]**, but the generator compiles.

---

## 7. The encryption reality

**Headline: DESIGN_NOTE.md's M6 risk statement — "Many US NextGen TV broadcasts encrypt
the media" — is true but pessimistically framed. About one third of subchannels are
encrypted, the figure has plateaued, it is concentrated in NBC and CBS network feeds, and
PBS encrypts nothing at all. Roughly two thirds of real services are legitimately
decodable, and everything below the media sample is in the clear permanently, by design.**

Nothing in this section describes or links to any circumvention method, and none was
researched. This is strictly about what is legitimately decodable.

### 7.1 There are two gates, and they are different shapes

**[FOUND]** A3SA (the ATSC 3.0 Security Authority) runs two separate programs, and the
distinction matters a great deal to us. From A3SA's own broadcaster document:

> "Signal and application signing are **required** in ATSC 3. Content encryption is
> **optional**, but is commonly required by media companies in their distribution
> agreements."

— https://a3sa.com/wp-content/uploads/2023/12/A-Short-Introduction-to-ATSC-3-Security-Systems-for-Broadcasters-2022.03.24.pdf

| | **Encryption** | **Signing ("High Noon")** |
|---|---|---|
| Mechanism | CENC / Widevine over the DASH media samples | X.509/CMS signatures over LLS and SLS |
| Reach | ~34% of subchannels | 100% of stations |
| Carried | media segments only | **in the clear**, verifiable by anyone |
| Effect on us | media samples undecodable | **none** |

**"High Noon"** is A3SA's name for the date after which *certified* receivers must refuse
to display (or must warn on) any station lacking a paid, annually renewed signing
certificate — encrypted or not. Originally 2025-06-30, delayed in March 2025. Reported
cost $998/station/year from a single practical source. Weigel Broadcasting (MeTV) raised
this with the FCC in GN 16-142 (letter 2025-08-27), and the FCC's own FNPRM ¶35 n.134
cites that *"A3SA will not even discuss issuing a license to broadcasters that do not
sign a non-disclosure agreement."*
Sources: https://blog.lon.tv/2025/09/02/atsc-3-update-high-noon-a-secret-broadcaster-plan-to-take-over-the-public-airwaves/ ·
https://www.newscaststudio.com/2025/09/04/weigel-broadcasting-challenges-atsc-3-0-security-authoritys-control-over-nextgen-tv/

**[INFER — and this is the useful conclusion]** High Noon is a *contractual obligation on
certified receivers*. It changes nothing in the RF. Signing data stays in the clear and
is verifiable by anyone (`libatsc3` already parses it). An uncertified receiver is simply
not bound by a rule that says "refuse to display". **High Noon does not raise the
encryption ceiling for us.**

### 7.2 How much is actually encrypted — measured, not quoted

**[FOUND — measured by this research, reproducible]** Published reporting on this stops
in February 2024, so the numbers were counted directly from the live
[RabbitEars ATSC 3.0 station list](https://www.rabbitears.info/market.php?request=atsc3)
plus 20 Wayback snapshots.

**184 of 535 over-the-air ATSC 3.0 subchannels encrypted = 34.4%** (2026-08-06).

| Snapshot | Encrypted / OTA subch. | % |
|---|---|---|
| 2022-12-26 | 0 / 381 | **0%** |
| 2023-02-13 | 17 / 387 | 4.4% |
| 2023-07-07 | 61 / 414 | 14.7% |
| 2023-12-02 | 99 / 451 | 22.0% |
| 2024-06-10 | 123 / 492 | 25.0% |
| 2025-08-01 | 171 / 530 | 32.3% |
| 2025-12-15 | 176 / 530 | 33.2% |
| **2026-08-06** | **184 / 535** | **34.4%** |

The series independently reproduces both published datapoints — "roughly 16 percent"
(Jul 2023) and "nearly 24 percent" (Feb 2024) from Jared Newman at TechHive (both
original URLs now 410; use the Wayback copies) — which is a good sign the counting method
is sound.

**It is network-driven, not group-driven:**

| Network | Encrypted / Total | Rate |
|---|---|---|
| NBC | 54 / 71 | **76%** |
| Telemundo | 7 / 10 | 70% |
| CBS | 44 / 65 | **68%** |
| ABC | 26 / 69 | 38% |
| FOX | 23 / 72 | 32% |
| CW | 9 / 49 | 18% |
| **PBS** | **0 / 37** | **0%** |

**And growth has plateaued** — near-linear 4%→33% between Feb 2023 and Dec 2025, then only
33.2%→34.4% in the eight months since. On-air encryption started in a tight window
between **2022-12-27 and 2023-02-13** (zero streams in the December snapshot, 17 in the
February one).

**[NOT FOUND]** No published source gives the current percentage; the table above is
original measurement and should be labelled as such if it is ever quoted.

### 7.3 What a receiver legitimately sees when a service is protected

**[FOUND]** Encryption is applied at the **DASH media-sample layer** (CENC common
encryption, Widevine). Everything below it is untouched:

- **Bootstrap, L1-Basic, L1-Detail, PLP structure** — clear (they must be; the receiver
  cannot demodulate otherwise).
- **ALP, IP, UDP** — clear.
- **LLS / SLT / SLS / S-TSID / MPD** — clear. **The channel lineup is always readable.**
- **ROUTE/LCT packet flow, object boundaries, segment arrival timing** — clear and
  observable.
- **Signing data** — clear by construction.
- **AEA emergency alerts** — an LLS table, always in the clear.

So the honest deliverable DESIGN_NOTE.md already specifies — *"present, identified,
locked"* — is exactly right, and it is achievable for 100% of services. For ~66% of them
we can go further and get video.

**[FOUND — a real-world caveat worth recording]** Weigel's 2026 bench tests found that
certified receivers (ZapperBox, HDHomeRun) on *encrypted* signals with the internet
disconnected showed a **"Secure Content License Server" error instead of the EAS monthly
test**, in four markets (WNDU South Bend, WFRV Green Bay, WMAQ Chicago, WSOC Charlotte),
while ATSC 1.0 on the same stations carried the alert.
https://www.tvtechnology.com/regulatory-legal/digital-alert-systems-details-atsc-3-0-eas-capabilities-to-fcc
**[INFER]** The alert reached the receiver; the receiver could not render the service
carrying it. A receiver of ours that reads AEA straight out of LLS would not have that
failure mode — which is a small but genuine argument for the thing we are building.

### 7.4 Is there a path for an open-source receiver?

**[NOT FOUND]** — no public path. A3SA licensing requires an NDA before discussion (per
the FCC's own footnote above), and DRM client certification is fundamentally incompatible
with distributing source. SiliconDust has filed repeatedly with the FCC on the
consequences, including a **2026-05-27** ex parte warning that a rule barring encryption
of only the *primary* channel would be a loophole, since subchannels could still be
locked. The FCC proceeding (GN 16-142) remains **open with no decision** as of July 2026.

**[INFER]** So the position DESIGN_NOTE.md already takes is the only available one, and it
is the correct one: enumerate encrypted services, label them locked, never attack them.

### 7.5 An unexpected opportunity: BPS is entirely below the DRM layer

**[FOUND]** The **Broadcast Positioning System** rides in the **bootstrap and preamble**
(`L1B_time_info_flag = 11b`, plus L1-Basic/L1-Detail fields) with a robust data PLP
carrying tower location. Bootstrap decodes at ~−12 dB SNR, preamble ~−9 dB; the location
PLP is the weak link at ~−5.72 dB.
https://www.nab.org/bps/Broadcast_Positioning_System_Using_ATSC30.pdf

On air since 2024 (KWGN Denver, ch. 34, 593 MHz). NIST measured tens-of-nanoseconds
peak-to-peak deviation over 50 days at up to 106 km:
https://www.nist.gov/publications/time-transfer-performance-broadcast-positioning-system-bps

**[INFER]** This sits entirely below CENC — no DRM touches it, and **any receiver that can
demodulate L1 can use it**. Our M1+M3 already decode bootstrap and L1. Given this lab's
existing GPS/timing work, ATSC 3.0 time transfer is a genuinely reachable secondary
result that is *immune to the encryption problem entirely*. Position fixes need ≥3
geographically separated towers; we have four ATSC 3.0 transmitters in range (RF25, 29,
30, 33). Worth a note in the ladder even if it is not pursued now.

### 7.6 Source-quality warning

**[FOUND]** A cluster of sites — antennaland.com, smarttvs.org, brainvoyage.blog,
factually.co — dominate search results for ATSC 3.0 DRM and produce plausible-sounding
but uncorroborated claims, including a **fabricated** "FCC voted in May 2026 for an
18-month mandatory migration with an Oct 1, 2027 deadline" that contradicts every primary
and legal-trade source. **That claim is false.** Everything in §7 is anchored to FCC
documents, A3SA's own publications, TV Tech, Lon.TV, NewscastStudio, or direct
measurement. Recorded here because this topic is unusually polluted and a future search
will hit these sites first.

---

## 8. After first convergence — what a naive chain omits

The single most important discovery in this section: **A/327 §6 is titled "Guidelines for
Receiver Implementation".** ATSC wrote a receiver cookbook and we have been building
M3–M5 without it. It covers bootstrap correlation topologies, channel estimation, PAPR
removal, frequency/time de-interleaving (including the two-memory ping-pong for 8K/16K
and the single-memory scheme for 32K), inverse CDL / TBI / cell de-interleavers, LDM
cancellation, **exact and Max-Log LLR formulas for 1D-NUC and 2D-NUC (§6.5.1)**, bit
de-interleaving, inner/outer decoding, and L1 decoding.

### 8.1 The gap is quantified in the standard: ~3 dB, and it is not optional

**[FOUND]** A/327 Annex B ("ATSC 3.0 Receiver C/N Model") budgets the distance between
the ideal BICM threshold and what a real receiver needs:

```
C/N_total = -10 log10( 10^-((C/N_BICM + A + Δpilot + IM)/10) - 10^(BN/10) )
```

- `A = 0.5 dB` — BER 1e-6 → quasi-error-free (1e-11)
- `Δpilot` ≈ 0.2–0.8 dB — boosted pilots steal power from data carriers
- `IM = L_ce + 0.5 dB` — real channel-estimation loss, plus imperfect LDPC, sync, and
  fixed-point losses
- `BN` — backstop noise, −33 dBc (QPSK…256QAM) / −38 dBc (1024/4096QAM)

**`L_ce` alone ranges from 2.43 dB with no pilot boosting down to ~0.31 dB with maximum
boost on a sparse pattern** (Table B.5.1). A/327's own worked example lands at
**ideal 22.2 dB → receiver spec 25.7 dB**.

**[FOUND]** And it is confirmed by measurement, not just modelled — A/327 Annex A gives
Simulation / Lab / Field C/N for a real commercial receiver:

| ModCod (64800, AWGN) | Simulation | Lab | Field |
|---|---|---|---|
| QPSK 2/15 | −6.1 | −4.0 | −3.9 |
| 64QAM 11/15 | 14.4 | 15.0 | 15.2 |
| 256QAM 11/15 | 18.9 | 19.7 | 19.9 |
| 4096QAM 13/15 | 33.1 | 33.4 | 34.0 |

Note the shape, which is counter-intuitive and matters to us: **the implementation gap is
~2 dB at the robust end and shrinks to ~0.3–0.8 dB at high SNR.** Sync loops and channel
estimation degrade fastest exactly where you would expect an easy first decode.

### 8.2 The things that will bite this chain specifically

Ordered by how likely they are to hit us, given what DESIGN_NOTE.md already records.

**(a) LLRs that ignore channel-estimate error.** **[FOUND]** "Estimation of channel MSE
for ATSC 3.0 receiver and its applications", ICT Express 8(2), 2022
(https://www.sciencedirect.com/science/article/pii/S2405959521001132). Our LLR divides by
σ². If we use the *channel* noise variance and ignore that Ĥ is itself noisy, the LLRs
are **overconfident**, and an LDPC decoder fed confident-but-wrong LLRs does not degrade
gracefully — it converges to the wrong codeword. The fix is one line:
`σ²_eff = σ²_noise + |x|²·MSE_ce`. The paper reports ~0.2 dB gain for 2-D estimation but
describes the correction as **essential to decode at all in the 1-D case** — and our
current estimator is 1-D (DESIGN_NOTE.md records that 2-D interpolation across the Dy
pair "gained nothing", which is itself consistent with a missing MSE term).

**(b) Min-sum flavour is worth up to 1.9 dB.** **[FOUND]** A/327 Table 6.3, 50
iterations, long codes, dB from Shannon:

| Rate | SPA | MSA | **OMSA** | NMSA |
|---|---|---|---|---|
| 6/15 | 0.73 | 1.87 | **1.05** | 2.59 (floor) |
| 11/15 | 0.46 | 1.05 | **0.56** | 0.82 |
| 13/15 | 0.41 | 0.85 | **0.48** | 0.52 (floor) |

A/327 explicitly warns that **normalized** min-sum has "inevitable error floors due to the
channel mismatch effect" and recommends **offset** min-sum: "the OMSA provides stable
performance without any error floor in all FER regions." **Table 6.1 gives the optimized
per-rate offset and scaling values** (64800: offsets 0.44–0.60; 16200: 0.31–0.63) — use
them rather than guessing 0.5. Short codes cost a further 0.15–0.64 dB vs long.

**(c) The frequency de-interleaver reset.** **[FOUND]** A/327 §6.3.1: the FI is reset at
the first preamble symbol of the frame **and at the first symbol of the 2nd and
subsequent subframes.** RF33 has two subframes. Miss this and subframe 1 is garbage while
subframe 0 looks perfect — a failure mode that would be very easy to misread as "PLP 1
needs more SNR".

**(d) SBS is a full-rate pilot symbol.** **[FOUND]** A/327 §5.4.1.5: SBS uses the same Dx
but **Dy is effectively 1**, so it carries the union of all scattered-pilot phases. It
must be special-cased in the estimator, not treated as a data symbol. (Our cell-map work
already handles SBS *cell counting* correctly; this is about the *channel estimate*.)

**(e) Channel-estimation holes at every boundary.** **[FOUND]** A/327 §5.3.2.2: going
sparse→dense across a subframe boundary, "the previous subframe's frequency-interpolated
carriers cannot be used… the result may be even worse than beginning channel estimation
from the subframe boundary itself." Expect a burst of bad LLRs at each subframe boundary
and at frame start. Signature: decodes fine, then a periodic glitch at the frame rate.

**(f) Time-interleaver warm-up.** **[INFER, mechanism-grounded]** A CTI with 1024 rows
holds 1024·1023/2 cells of state; the first ~200 ms after acquisition produces garbage no
matter how clean the signal. HTI has an analogous fill. **Do not read the first TI block
as a decoder bug.** (For our intra-subframe HTI case the fill is one TI block, which is
much shorter — but it is not zero.)

**(g) PAPR tone reservation.** **[FOUND]** A/327 §6.2.2: if `L1B_papr_reduction` signals
Tone Reservation, the receiver must **discard** the reserved carriers; if it signals ACE,
no action is needed. Treating TR carriers as data misaligns the entire cell-to-carrier
map — a total failure, not a degradation. Worth checking what RF33 signals before
debugging anything else.

**(h) SFO is what breaks first over minutes.** **[INFER — flagged as such by the
research; A/327 §6.1 covers acquisition only and gives no steady-state tracking
algorithm.]** Continual pilots exist for this: A/327 §6.2.1 says they are "usually used
for system frequency synchronization, rather than equalization" — so they should **not**
go into the equalizer interpolator. The standard approach is to fit a line to
continual-pilot phase across carrier index per symbol: **intercept = common phase error
(residual CFO), slope = sample clock offset**. At 16K/32K the outer carriers accumulate
phase thousands of times faster than the centre, so the signature is a constellation that
looks clean near DC and smears at the band edges, worsening across the frame — and
because the bootstrap re-acquires every ~250 ms, a naive chain sawtooths rather than
failing cleanly. **Recommendation: log CPE and SFO-slope per symbol from the very first
payload decode.** Two scalars, and they are the ATSC 3.0 analogue of the existing MER
dial.

### 8.3 The transport layer is bigger than the PHY, and AL-FEC will not save a weak link

**[FOUND]** A/331:2025-06 §8.1.1.6 states plainly:

> "The application of AL-FEC is of marginal benefit for linear streaming Service content
> that has object physical layer durations in the range of the soft-decision FEC Frame in
> the physical layer."

and gives the arithmetic that sets the real target:

> "Users may tolerate up to 1% lost time for streaming video… **a per-FEC Frame error
> rate of ~4e-5 is required to achieve 1% lost time.** For this FEC Frame loss rate, the
> probability of a successful recording of an hour-long show with no lost seconds is
> ~2e-16 without AL-FEC."

**This is the number that replaces "it decoded once".** Per-FEC-frame error rate ≲ 4e-5,
which is directly analogous to the existing STVT delivery-rate/quality dial. RaptorQ
repair flows (A/331 Annex A.4) exist mainly for NRT files, DVR recording and ESG objects
— and **may not be present in a given broadcast at all**.

Other transport realities worth knowing before M6 starts:
- **[FOUND]** A/331 §8.1.1.4–5: three buffers in series (ROUTE transport buffer → ROUTE
  output buffer → DASH elementary stream buffer), and the system model is defined such
  that none may overflow and the decoder may not stall. Notably, **"a notable aspect of
  the ROUTE/DASH System Model is the lack of a physical layer buffer"** — the PHY delivers
  discrete Data Delivery Events, not a byte stream.
- **[FOUND]** §8.1.1.1: the timing anchor is `EXT_ROUTE_PRESENTATION_TIME` on the first
  LCT packet of an MDE block containing a RAP, relative to `SCT` from the `EXT_TIME`
  extension. `(EXT_ROUTE_PRESENTATION_TIME − SCT)` is exactly how long to hold that block
  for stall-free playback. This is the piece `lukefay/ATSC3_Rx` implements as its AST
  alignment (§5).
- **[FOUND]** Media segments are designed to be decoded **partially** — MDE blocks at
  video-frame granularity. A "wait for the whole segment" receiver adds ~2 s of latency
  and is *more* fragile, not less. Exception: **IMSC1 captions require the complete
  segment** before decoding (§8.1.1.7).
- **[INFER]** Loss propagation is brutal by design: a failed BCH discards a whole baseband
  frame → the ALP packets it carried are lost → A/330 ALP segmentation means one lost
  segment kills the whole IP packet → UDP checksum drops the datagram → an LCT packet is
  missing → the ROUTE object is incomplete. **One failed 64800-bit FEC frame ≈ several
  contiguous kilobytes of IP payload gone.**

### 8.4 And then the codecs

**[FOUND]** Even with perfect bits, playback is not free. `johnb-7/hdhr-ac4` (98★) exists
solely because **ATSC 3.0 audio is Dolby AC-4**, and its author's note is the whole story:
*"The signal and video quality were much better, but the AC4 audio every program was using
was supported by nothing I owned."* `ferrellsl/VideoPlayer` and `plazareff/VideoPlayer`
both depend on an **experimental AC-4 branch of ffmpeg**. Video is HEVC (fine).

**[INFER]** Budget for this: our chain can be bit-perfect and still produce silent video.
Check ffmpeg/mpv AC-4 support on the fleet *before* M6, and note that
`lukefay/ATSC3_Rx` documents the browser-side probe
(`MediaSource.isTypeSupported('audio/mp4; codecs=ac-4.02.01.00')`) as the practical test.

### 8.5 Architectural note that validates our milestone split

**[FOUND]** Sony's ATSC 3.0 demodulator LSIs (CXD2878ER, CXD6801GL) present a host
interface of **"TS / ALP"** — the PHY chip hands off ALP packets and the SoC runs the
entire A/331 stack
(https://www.atsc.org/wp-content/uploads/2021/01/f-36-26-13345252_W2dnWCqD_sony_ATSC3.0_receiverLSI_rev1.0.pdf).
**[INFER]** That is exactly the M4/M5 seam in the ladder, and it is the right place to
hand off to `libatsc3`. Our architecture matches commercial silicon.

---

## 9. Link budget: the numbers our capture should be judged against

**[FOUND]** A/327:2025-06 Annex B.6 publishes **expected AWGN threshold C/N** for every
FFT size / GI / modulation / code rate. These are the numbers DESIGN_NOTE.md's
"18.8 dB MER is about the AWGN threshold for 64QAM 11/15" was estimating.

For our actual RF33 configuration:

| PLP | Config | Table | **Threshold C/N** | Our measurement | Margin |
|---|---|---|---|---|---|
| PLP 0 | 8K FFT, GI 1536 (Dx=4), 64QAM, **16200**, 11/15 | B.6.3 | **17.5 dB** | 18.78 dB MER | **+1.3 dB** |
| PLP 1 | 16K FFT, GI 1536 (Dx=4), 256QAM, **64800**, 11/15 | B.6.9 | **22.0 dB** | not yet measured | — |

Two consequences:

- **PLP 0 has positive margin, not zero margin.** DESIGN_NOTE.md's framing ("AT the
  cliff, not comfortably above it") is right in spirit but the actual figure is 1.3 dB in
  hand. **[INFER]** With the TBI corrected per §1, a correct chain has a real chance of
  converging on the existing capture — so getting more antenna margin, while still
  worthwhile, may not be the hard prerequisite the note assumes.
- **PLP 1 needs 22.0 dB** — 4.5 dB harder than PLP 0. Attack PLP 0 first, which is
  already the plan.

**[INFER — caveat]** Threshold C/N is a channel measurement; our 18.78 dB is a
post-equaliser MER. Implementation loss means these are not the same quantity and MER
typically flatters the link slightly. Treat +1.3 dB as "plausibly enough", not "proven
enough".

### 9.1 Two different tables — do not mix them up

A/327 publishes **two** sets of numbers and they differ by ~3 dB:

| | what it is | 64QAM 11/15 (16200) |
|---|---|---|
| **Tables 4.3–4.6** | **ideal BICM** threshold, perfect CSI, 50 iterations, BER 1e-6 | 14.52 dB (AWGN) |
| **Annex B.6** | **modelled receiver** expected C/N, per FFT/GI, incl. `A + Δpilot + IM + BN` | **17.5 dB** |

The ~3 dB between them is exactly the Annex B budget in §8.1. **Annex B.6 is the one to
compare a real measurement against** — which is what §9's table uses. (The Table 4.3–4.6
figures quoted here are second-hand and the source PDF's column layout is awkward to
extract; treat them as **[CLAIM]** pending our own PyMuPDF extraction, which is now
trivial.)

### 9.2 Why L1 decodes cleanly while the payload may not — and why that is normal

**[FOUND]** A/327 Table 4.10 gives L1 signaling thresholds at FER 1e-4. For FEC mode 1
they are around **−9 dB** (L1-Basic −9.2 AWGN, L1-Detail −9.0), and the bootstrap
threshold is **≈ −9.5 dB**. Our own M1 result — detection to −18 dB with processing gain
— sits comfortably below even that.

So the correct ordering is: **bootstrap locks ~9 dB below L1, and L1 locks ~24 dB below a
64QAM 11/15 payload.** A chain that reads L1-Detail with a clean CRC and fails to
converge on the payload is **behaving exactly as designed**, and that observation alone
does not indicate a bug. It also means the L1 success in M3/M4 — impressive as it is —
carries almost no information about whether the payload chain is correct. Only the
payload can test the payload.

A/327 §4.2.10 adds that **L1 has no time diversity** unless L1-Detail Additional Parity
is enabled, which is one of the ten features gr-atsc3 does not implement.

Source: https://www.atsc.org/wp-content/uploads/2025/07/A327-2025-06-Physical-Layer-RP.pdf
(Tables B.6.3 and B.6.9, PDF pages 129 and 135).

---

## 10. Honest novelty assessment

**The "world first" claim is not supportable. A completed, documented, end-to-end
hobbyist SDR→video ATSC 3.0 chain would still be genuinely rare and, as far as this
search can establish, without a public open-source peer.**

What the evidence supports, precisely:

1. **ATSC 3.0 receivers exist in quantity — as silicon and products.** Sony ships
   demodulator parts (CXD2878ER, CXD6801GL); SiliconDust HDHomeRun, ZapperBox, and
   NextGen TV sets are consumer products; Sencore and DekTec sell professional receivers.
   Nobody should describe decoding ATSC 3.0 as unsolved. **[FOUND/CLAIM — vendor pages.]**

2. **In open source, the PHY receive path is genuinely unoccupied.**
   - Core GNU Radio: ATSC 1.0 only. **[FOUND]**
   - `gr-atsc3`: transmitter, explicitly. **[FOUND]**
   - The academic GRC implementation (Dionisio & Akamine 2017): modulator, validated
     against a commercial receiver. **[FOUND]**
   - `libatsc3` / `ATSC3_Rx` / `ATSC_ROUTE`: everything **above** IP. **[FOUND]**
   - The one repo claiming a full PHY receiver is a scaffold with a synthetic LDPC
     matrix and no HTI (§3.2). **[FOUND]**
   - A GitHub **code search for `Nfec_ti_max` returns 10 hits, all of them files inside
     `drmpeg/gr-atsc3`** — i.e. exactly one public codebase on GitHub implements the
     ATSC 3.0 hybrid time interleaver at all, and it is the transmitter. **[FOUND]**

3. **So the accurate claim is narrow and defensible:**
   > *No public open-source ATSC 3.0 receiver taking RF to decoded video appears to
   > exist. A working one would be, as far as I can find, the first published
   > open-source implementation of the ATSC 3.0 receive PHY — and the first
   > independent-of-vendor demonstration that the standard can be received with a
   > general-purpose SDR and open code.*

   That is a real result. It is **not** "world first ATSC 3.0 reception", and the
   difference matters: the first is checkable and will survive contact with someone who
   works in broadcast; the second will not.

4. **[INFER] Why the space is empty is not mainly difficulty.** DVB-T2 — comparably
   complex, same LDPC/BCH family, same OFDM-with-pilots structure — *does* have open
   receivers (`gr-dvbs2rx` covers DVB-T2). The differences specific to ATSC 3.0 are
   (a) it is US-only in practice, so the hobbyist population is much smaller than DVB's,
   and (b) **the payoff is partly capped by encryption** (§7) — about a third of
   subchannels, concentrated in exactly the NBC/CBS network feeds a hobbyist would most
   want. That is a real disincentive, and the *perception* of it is probably larger than
   the measured 34%, since the widely-read coverage of NextGen TV DRM does not report a
   percentage at all.

   **Correction to my own first framing:** I initially wrote that "a large fraction of the
   video you would decode is DRM-protected". The measurement in §7.2 does not support
   "large fraction" at the subchannel level — it is ~34%, it has plateaued, and PBS is at
   zero. Two thirds of services are decodable. The disincentive is real but smaller than
   the reputation.

5. **What is unambiguously novel regardless of the above:** the *method* record. A
   milestone-gated, control-driven reverse-engineering of A/322 from public PDFs, with
   falsifiable gates at every rung and a written ledger of the wrong turns, is not
   something the census turned up anywhere. **[INFER]** DESIGN_NOTE.md is arguably the
   more publishable artifact than the decoder.

---

## 11. Prioritised — what would actually save us time

Ranked by value. Items 1–3 are the ones that change what happens next.

### 1. Replace the TBI with A/322's own equations. Do not run the ST sweep. *(saves: the sweep, plus an unbounded debugging tail)*
Implement §1.1 exactly:
```python
R = i % Nr                     # Nr = cells per FEC block (2700)
T = R % Nc                     # Nc = N_FEC_TI_MAX (37)
C = (T + i // Nr) % Nc
theta = Nr * C + R
if theta >= N_FEC_TI_DUMMY * Nr:   # else skip: virtual cell
    emit(theta - N_FEC_TI_DUMMY * Nr)
```
Gate it on the A/327 Figure 6.5(a) vector (§1.5) **before** any air data:
`4×3, one virtual block → b g a f d e c h`. Then invert for the de-interleaver.
Close ASSUMPTION T1 as *settled by the spec, reading corrected*, and record in the bug
ledger that a gate built on prose contradicted by the adjacent equations is not a gate.

### 2. Switch spec extraction to PyMuPDF and re-open every "unreadable" item. *(saves: an unknown but large amount of future reverse-engineering)*
`pip install pymupdf`; write the normaliser (collapse doubled glyphs; map the PUA
bracket codepoints to `floor`/paren); re-extract A/322 §7.1.5.x and then sweep the whole
document for equations previously lost. Re-test assumption **B2** (Tables 6.12/6.13)
specifically. Add a gate that the PyMuPDF extraction of every already-gated table
reproduces the `pdftotext` result, so the new extractor is itself validated against
known-good output before it is trusted on unknown output.

### 3. Build gr-atsc3 and make it our referee. *(saves: the "chain wrong vs SNR short" ambiguity, permanently)*
```
git clone https://github.com/drmpeg/gr-atsc3
curl -O https://www.w6rz.net/advatsc3.ts        # 281 MB, verified live
```
Open `examples/vv031.grc`, swap `uhd_usrp_sink` → file sink, and in the
`atsc3_framemapper_cc` block set:

| param | vv031 default | set to | why |
|---|---|---|---|
| `timode` | `TI_MODE_OFF` | **`TI_MODE_HYBRID`** | vv031 does **not** exercise the TBI at all |
| `tiblocks` | 2 | 2 | matches RF33 PLP 0 (`NTI = 2`) |
| `tifecblocksmax` | 14 | 37 | matches RF33 PLP 0 `N_fec_TI_max` |
| `tifecblocks` | 14 | **e.g. 35** | `< max` ⇒ **non-zero virtual FEC blocks** |

That last row matters more than it looks: RF33 PLP 0 has `N_virtual = 0`, so the
virtual-cell skip path in §1.1 has only ever been tested by our own synthesis (gate T2).
Setting `tifecblocks < tifecblocksmax` exercises it against an independent implementation.
Now every stage of our chain has a known-truth input at a chosen SNR. Offline, no radio,
no fleet lock, no saboteur risk.

### 4. Diff our extracted tables against gr-atsc3's — as an oracle, never as a source. *(saves: latent silent-wrong-table bugs)*
`spec_ldpc.py` ↔ `lib/ldpc_bb_impl.h`; `spec_bitint.py` ↔ `lib/interleaver_bb_impl.h`;
`spec_nuc.py` ↔ `lib/modulator_bc_impl.h`; frequency de-interleaver ↔
`bitperm{8,16,32}k*`; pilot tables ↔ `pilotgenerator_cc_impl.cc`. Record the diff count
as a gate. **Keep the clean-room boundary: read to compare, never paste** (GPL-3.0, §4).

### 5. Read A/327 §6 — it is a receiver cookbook and we have been ignoring it. *(saves: rediscovering receiver-side geometry, repeatedly)*
https://www.atsc.org/wp-content/uploads/2025/07/A327-2025-06-Physical-Layer-RP.pdf —
185 pages, public, extracts cleanly with PyMuPDF. **§6 is literally titled "Guidelines
for Receiver Implementation".** It contains: bootstrap correlation topologies and the
peak-to-average validation gate (§6.1); channel estimation and the continual-pilot role
(§6.2.1); PAPR/TR carrier removal (§6.2.2); the frequency de-interleaver **reset rule**
(§6.3.1); time de-interleaving including the single-memory TBDI, inverse CDL and cell
de-interleaver (§6.3.2 — the section that answers §1 of this document); LDM cancellation
(§6.4); **exact and Max-Log LLR formulas for 1D-NUC and 2D-NUC (§6.5.1)**; bit
de-interleaving (§6.5.2); LDPC decoder algorithm comparison and per-rate offset/scaling
tables (§6.5.3). Annex A has lab/field measured C/N; Annex B has the receiver C/N model
and the B.6 threshold tables (§9); Annex C has reference PLP configurations; Annex E has
channel models. **We have been building M3–M6 without the document ATSC wrote to explain
M3–M6.**

### 6. Collect the cheap dB before blaming the antenna. *(saves: an antenna campaign that may not be needed)*
Four items, all small code changes, all cited in §8.2:
- **Fold channel-estimate MSE into the LLR denominator** (`σ²_eff = σ²_noise + |x|²·MSE_ce`).
  Cited as *essential* for 1-D estimation, which is what we have.
- **Use offset min-sum with A/327 Table 6.1's per-rate offsets**, layered scheduling, 50
  iterations, early termination on syndrome + BCH. Worth up to ~1.9 dB vs plain min-sum,
  and A/327 explicitly warns off *normalized* min-sum (error floors).
- **Check `L1B_papr_reduction`** — if Tone Reservation is signalled, those carriers must
  be discarded or the entire cell map shifts.
- **Reset the frequency de-interleaver at subframe 1** — RF33 has two subframes, and
  missing this breaks PLP 1 only, which looks exactly like "PLP 1 needs more SNR".

Given §9's +1.3 dB of existing margin on PLP 0, these may be the difference between a
converging chain and another antenna campaign.

### 7. Bank the free test vectors now, before they are needed. *(cheap, and they expire)*
- The A/327 TBI example (§1.5) — into the selftest today.
- A/322 §7.1.5.2 cell-interleaver `P(r)` sequence (§1.8) — free gate on the cell
  interleaver, needed for PLP 1.
- `junhuac/atsc-3.0-mmt-pcaps` — real pre-encryption ROUTE/DASH captures, LFS objects
  verified live today (§6.3). **Mirror them locally now**; a 2019 repo with a dead
  maintainer is one LFS-quota change away from gone. Note the standing law that
  `.gitignore` hides data from backups — keep these where the backup counts them.

### 8. Do not rebuild M5/M6 from scratch. *(saves: most of two milestones)*
`libatsc3` (MIT — vendorable) for ALP/LLS/ROUTE/MMT; `lukefay/ATSC3_Rx` for the
SLT→ROUTE→DASH→player pattern in Python, including the DASH `availabilityStartTime`
alignment that is easy to get wrong; `wireshark-atsc3` for visibility the moment IP
packets exist. Read all three before writing M5.

### 9. Do not touch `danttoUHD`. *(risk avoidance)*
It is scoped to a Korean LG STB container and **includes decryption functionality**.
It is outside this project's stated scope and outside the standing position that
encrypted services are labelled locked, never attacked. Recorded here only so it is
not rediscovered later and mistaken for a useful ATSC 3.0 remuxer.

### 10. Correct the novelty language in README.md before anyone reads it. *(credibility)*
Replace any "world first" framing with §10's narrow claim. The narrow claim is stronger
precisely because it survives scrutiny from someone who works in broadcast.

---

## 12. Sources

**Standards (all public, all direct downloads, all extract with PyMuPDF)**
- **A/322:2026-06 (current)** — https://www.atsc.org/wp-content/uploads/2026/06/A322-2026-06-Physical-Layer-Protocol.pdf
- A/322:2026-04 Physical Layer Protocol — https://www.atsc.org/wp-content/uploads/2026/04/A322-2026-04-Physical-Layer-Protocol.pdf
- A/322:2023-03 — https://www.atsc.org/wp-content/uploads/2023/09/A322-2023-03-Physical-Layer-Protocol.pdf
- A/322:2021 — http://www.atsc.org/wp-content/uploads/2021/04/A322-2021-Physical-Layer-Protocol.pdf
- A/322:2018 — https://www.atsc.org/wp-content/uploads/2016/10/A322-2018-Physical-Layer-Protocol.pdf
- A/322:2016 — https://www.atsc.org/wp-content/uploads/2016/10/A322-2016-Physical-Layer-Protocol.pdf
- **A/327:2025-06 Guidelines for the Physical Layer Protocol** — https://www.atsc.org/wp-content/uploads/2025/07/A327-2025-06-Physical-Layer-RP.pdf
- A/327:2023-03 — https://www.atsc.org/wp-content/uploads/2023/04/A327-2023-03-Physical-Layer-RP.pdf
- A/325:2023-04 PHY Lab Performance Test Plan — https://www.atsc.org/wp-content/uploads/2023/04/A325-2023-04-Lab-Performance-Test-Plan.pdf
- A/326:2025-07 Field Test Plan — https://www.atsc.org/wp-content/uploads/2025/08/A326-2025-07-Field-Test-Plan.pdf
- **A/331:2025-06 Signaling, Delivery, Synchronization, and Error Protection** — https://www.atsc.org/wp-content/uploads/2025/06/A331-2025-06-Signaling-Delivery-Sync-FEC.pdf
- A/330:2025-07 Link-Layer Protocol (ALP) — https://www.atsc.org/wp-content/uploads/2025/08/A330-2025-07-Link-Layer-Protocol.pdf
- ATSC 3.0 standards index — https://www.atsc.org/atsc-documents/type/3-0-standards/

**Papers**
- S. Jeon et al., "Physical Layer Time Interleaving for the ATSC 3.0 System", IEEE Trans. Broadcasting — https://ieeexplore.ieee.org/document/7378942/ · PDF: https://scispace.com/pdf/physical-layer-time-interleaving-for-the-atsc-3-0-system-1i0adahoh7.pdf
- "Verification and validation of the physical layer ATSC 3.0 standard", IEEE BMSB 2016 — https://ieeexplore.ieee.org/document/7521999/
- V. M. Dionisio, C. Akamine, "ATSC 3.0 implementation in GNU Radio Companion", IEEE 2017 — https://ieeexplore.ieee.org/document/7986220/
- "Estimation of channel MSE for ATSC 3.0 receiver and its applications", ICT Express 8(2) 2022 — https://www.sciencedirect.com/science/article/pii/S2405959521001132
- Ahn et al., "Evaluation of ATSC 3.0 and 3GPP Rel-17 5G Broadcasting Systems for Mobile Handheld Applications", IEEE Trans. Broadcasting — https://mys.mapyourshow.com/mys_shared/nab23/handouts/Evaluation_of_ATSC_3.0_and_3GPP_Rel-17_5G_Broadcasting_Systems_for_Mobile_Handheld_Applications.pdf
- Fay et al., "An Overview of the ATSC 3.0 Physical Layer Specification", IEEE Trans. Broadcasting 62(1), 2016 — https://www.researchgate.net/publication/290472101_An_Overview_of_the_ATSC_30_Physical_Layer_Specification
- "Performance Analysis of Practical QC-LDPC Codes: From DVB-S2 to ATSC 3.0", IEEE Trans. Comm. — https://ieeexplore.ieee.org/document/8554116/
- "Efficient Decoding of LDM Core Layer at Fixed Receivers in ATSC 3.0", IEEE Trans. Broadcasting — https://ieeexplore.ieee.org/document/8000374/
- Sony ATSC 3.0 receiver LSI overview (host interface is "TS / ALP") — https://www.atsc.org/wp-content/uploads/2021/01/f-36-26-13345252_W2dnWCqD_sony_ATSC3.0_receiverLSI_rev1.0.pdf

**Encryption / policy**
- A3SA, "A Short Introduction to ATSC 3 Security Systems for Broadcasters" — https://a3sa.com/wp-content/uploads/2023/12/A-Short-Introduction-to-ATSC-3-Security-Systems-for-Broadcasters-2022.03.24.pdf
- RabbitEars ATSC 3.0 station list (the encryption dataset in §7.2) — https://www.rabbitears.info/market.php?request=atsc3
- Lon.TV on "High Noon" — https://blog.lon.tv/2025/09/02/atsc-3-update-high-noon-a-secret-broadcaster-plan-to-take-over-the-public-airwaves/
- NewscastStudio on the Weigel FCC filing — https://www.newscaststudio.com/2025/09/04/weigel-broadcasting-challenges-atsc-3-0-security-authoritys-control-over-nextgen-tv/
- TV Tech, Digital Alert Systems / EAS on encrypted ATSC 3.0 — https://www.tvtechnology.com/regulatory-legal/digital-alert-systems-details-atsc-3-0-eas-capabilities-to-fcc
- NAB, "Broadcast Positioning System Using ATSC 3.0" — https://www.nab.org/bps/Broadcast_Positioning_System_Using_ATSC30.pdf
- NIST, BPS time-transfer performance — https://www.nist.gov/publications/time-transfer-performance-broadcast-positioning-system-bps

**Code and data**
- https://github.com/drmpeg/gr-atsc3 · test stream https://www.w6rz.net/advatsc3.ts
- https://github.com/drmpeg/dtv-utils (includes `atsc3rate.c` bit-rate calculator)
- https://github.com/kansonkong/libatsc3 · https://github.com/junhuac/libatsc3 · https://www.ngbp.org/
- https://github.com/lukefay/ATSC3_Rx · https://github.com/haudiobe/ATSC_ROUTE
- https://github.com/LiamRPower/wireshark-atsc3
- https://github.com/junhuac/atsc-3.0-mmt-pcaps
- https://github.com/Silicondust/atsc3-sqa
- https://github.com/zrcheng1991/atsc-3.0 · https://github.com/philburr/atsc3
- https://github.com/awanga/usrp-atsc3 (assessed, §3.2)
- https://github.com/johnb-7/hdhr-ac4 · https://github.com/ferrellsl/VideoPlayer
- GNU Radio ATSC wiki (1.0 only) — https://wiki.gnuradio.org/index.php/ATSC

---

## 13. Addendum 2026-11-06 — OnAir, a second independent receiver

Research date: **2026-11-06**. This section follows the §4 rule: **comparison oracle,
not a source. Never paste.** OnAir is GPL-3.0 and this project is Apache-2.0, so the
only things that cross over are ideas, written here from the principle and gated by our
own referee. Credit is recorded in `NOTICE` ("Techniques learned from other projects").

### 13.0 Provenance — what was ours before we had seen it

**[FOUND — this repository's and Software-TV-Tuner's git history]**

| Date | Milestone | Evidence |
|---|---|---|
| 2026-04-29 | The TV tuning journey begins: ATSC 1.0, a gr-atsc fork with five rewritten blocks | Software-TV-Tuner, initial commit |
| 2026-05-01 | First broadcast decoded by our own chain ("HDHR-quality SDR decode", RF34) | Software-TV-Tuner, README milestone commit |
| 2026-05-02/03 | First stable, unattended decoding: continuous playback, watchdogs, seamless recovery | Software-TV-Tuner commits of those days |
| 2026-08-05 | ATSC 3.0 begins, alongside ATSC 1.0: bootstrap detector, signalling and L1-Basic decoded from real RF33 air the same day (M1–M3) | this repository, first commits |
| 2026-08-22/23 | Live ATSC 3.0 television with 5.1 audio, on a Raspberry Pi 5 | this repository |
| **2026-11-06** | **First idea taken from OnAir** (the delay-domain channel estimate, §13.3–13.4) | this section |

OnAir is not mentioned anywhere in either history before 2026-11-06. Everything above
that line was built independently, from the standards; everything from that line on
says so where it applies.

### 13.1 Census row

| Project | What it actually does | Lang | License | Maintained | Covers our gaps? |
|---|---|---|---|---|---|
| [adswill/OnAir](https://github.com/adswill/OnAir) | Multi-mode SDR receiver with an **ATSC 3.0 receive chain built from the spec**: bootstrap, L1, every subframe (8K/16K/32K), CTI and HTI, LDPC 16200 and 64800 (a batched 8-bit layered min-sum), ALP, LLS/SLT, ROUTE/DASH to a TS player. No MMTP, no AC-4. Its README says its ATSC 3.0 was tested only on simulated signals. | C++ | GPL-3.0 | Yes — 2026-11 | **Subframe 1 / PLP 1 live: YES.** Delay-domain channel estimate: YES. MMTP: no. |

### 13.2 Head to head, same 20 s of RF33 air, same PC (2026-11-06)

**[FOUND — measured]** Recording `~/captures/race_2026-11-06/live_rf33_antB.cs16`
(Antenna B, bias-T, rfgain_sel 4, IFGR 32, 6.912 Msps, no overflows).

| | ours (`python -m atsc3 watch --capture`) | OnAir (`atsc3rx`, 12 threads) |
|---|---|---|
| wall time for 20 s of air | 10.8 s | 18.0 s |
| PLP 0 FEC blocks | 5994/5994 | 5920/5920 |
| PLP 1 FEC blocks | not decoded live (m13 offline: 819/819 on the frames checked) | 9046/9360 — all 314 failures in two bursts (frames 20-21, 28-30) while power, flutter and PLP 0 were steady: a receiver fault, not the air |
| service played | WJLA 7.1 (MMTP), 2 audio tracks, captions | WHUT 32.1 (ROUTE), 1080p60, no audio (no AC-4) |

### 13.3 What is worth learning from it, in order

1. **Channel estimation in the delay domain** (their `atsc3_ofdm.cpp`, a standard
   DFT-based estimate). Ours interpolates linearly between pilots. Re-implemented as
   `lab/e104_combce.py`; see §13.4 for the gate.
2. **Every subframe decoded live.** Our live path is subframe 0 only; m13 proves PLP 1
   offline. The port is wiring (m13 + m16 demapper + the staged E47 ROUTE lanes). Found
   on the way: `m16_margin.smoothed_pool_g` (E60, lab copy) does NOT work on subframe 1
   (0/234 blocks on clean air at W=12); the live copy in `m9_fast` derotates per symbol
   and the lab copy does not. The port must use the live copy.
3. **Batched 8-bit layered LDPC**, 16 blocks per vector. A speed lever that only pays on
   weak signals where blocks run many iterations; our kernel retires clean blocks on the
   syndrome check with no work, so on clean air there is nothing to gain.

Not worth learning: their PLP 1 tracking (the bursts above), and they have no MMTP or AC-4.

### 13.4 E104 gate — delay-domain estimate vs linear interpolation

Protocol (E58's, adapted to real air): seeded complex Gaussian noise added to the IQ of
`hit_rf33.cs16` at stated levels below the recorded power, subframe 1 / PLP 1 (256QAM
11/15, 117 blocks per frame, 2 frames), decoded twice on the per-symbol path with the SAME
noise — linear interpolation (control) and the comb estimate — judged by LDPC convergence
and BCH syndrome zero. Self-check legs: flat channel exact; two-path channel 8.5-11 dB
more accurate than linear at 10/20/30 dB pilot SNR; pure noise refused; a comb too
sparse for the guard refused.

| added noise, dB below signal | linear | comb |
|---|---|---|
| none / 24 | 234/234 | 234/234 |
| 23.5 | 233/234 | 234/234 |
| 23.0 | 198/234 | 234/234 |
| 22.5 | 161/234 | 234/234 |
| 22.0 | 31/234 | **234/234** |
| 21.5 | 0/234 | **198/234** |
| 21.0 | 0/234 | **108/234** |
| 20 and below | 0/234 | 0/234 |

Result: linear interpolation holds 95% of the blocks down to 23.5 dB, the comb estimate
holds 100% down to 22.0 dB: **+1.5 dB of margin** on this signal (the 50%-of-blocks
point moves by about the same). Raw rows: `lab/e104_out/`. Caveat: the
"added noise" axis is relative to the recorded power, which already includes the
recording's own noise, so only the DIFFERENCE between the columns is a measurement.
### 13.5 Into the live chain — what broke first, and the fix (2026-11-06, evening)

Wired as `MarginCfg.ce_comb` (env `ATSC3_CE_COMB`, OFF by default): on the E60
time-smoothed grid the accepted symbols get a full-resolution delay-domain row, which the
C equaliser kernel reads through an identity grid map (every carrier its own grid point),
so `cellpool_kernel.c` is untouched. `ATSC3_CE_COMB_LOG=1` reports how many symbols the
estimator accepted; `ATSC3_SNR_LOG=1` prints the dummy-cell SNR of every frame.

**First result: it LOST.** 18 dB of added noise, PLP 0, 32 frames: 494/2368 blocks with
the estimate on against 1259/2368 off, every symbol accepted (35/35). The estimator was
right and the integration was right; the recipe was wrong for this regime.

**Diagnosis (synthetic grid, spacing 4, measured):** the standard recipe's cap of "never
keep a tap more than 40 dB under the strongest" leaves a bias floor near −31 dB, because
the band-edge bridge leaks a little energy into sidelobe skirts around every real tap and
the cap throws the skirts away. On the per-symbol path the pilot noise is far above that
floor and the gate showed +1.5 dB; on the smoothed grid the noise is already BELOW it on
a mild channel (linear interpolation reaches −41 dB at a 40 dB grid), so the estimate
replaced a small error with a larger one. Two changes remove the floor: (1) rotate the
grid so the main path sits at zero delay before bridging (the edges then match), and
(2) cap at 60 dB and keep 8 bins of skirt either side of every kept tap, the refusals
still judged on the taps that cleared the floor on their own. After that the estimate is
better than or equal to linear interpolation in every synthetic case tried (3–7 dB on a
mild channel at a 40 dB grid, 17–23 dB with a strong echo) — and the live chain at 18 dB
turned round to 1471/2368 against 1228/2368.

**A referee lesson.** Block counts at the cliff are NOT a stable referee for the live
chain: a re-acquisition ("8 frames under 60/74") throws away seconds, so identical input
gave 1471 and 805 blocks on two runs. The live gate therefore reads the per-frame
dummy-cell SNR at a level where both settings decode fully (the estimate's effect on the
equalised cells, free of acquisition luck), plus repeated runs at the cliff for the
trend. Results:

**Live chain, dummy-cell SNR per frame (32 frames each, same frames, PLP 0 on the
E60 smoothed grid, `ATSC3_SNR_LOG=1`):**

| recording | estimate off | estimate on | gain |
|---|---|---|---|
| 20 dB added noise | 16.13 dB (min 14.58) | 16.23 dB (min 14.66) | **+0.10 dB** |
| 19 dB added noise | 15.55 dB (min 14.20) | 15.69 dB (min 14.30) | **+0.14 dB** |
| as recorded | 15.75 dB (min 13.79) | 15.90 dB (min 13.94) | **+0.15 dB** |

Cell-pool stage: 20.4 ms/frame off, 55.6 ms/frame on (of a 247 ms budget).

**Cliff repeats** (blocks of 2368, PLP 0): 17 dB — off 1813/1686/1686, on 2131/1663/1663;
18 dB — off 1259/803/1469, on 805/805/297 (re-acquisitions 1–3 in every run). Identical
run to run where nothing re-acquired; otherwise acquisition luck, as said above.

**Verdict.** On subframe 0 the E60 time smoothing has already taken most of the pilot
noise, and the delay-domain estimate adds only **+0.1 to +0.15 dB** on top of it,
consistently, for about +35 ms a frame. `ce_comb` therefore STAYS OFF by default
there. Its real value in this chain is the **per-symbol path**, where the gate measured
**+1.5 dB** (§13.4): the fallback symbols the change detector sends that way, and all of
**subframe 1 / PLP 1**, which has no working time smoothing yet (§13.3 item 2). The PLP 1
live port should take this estimator as its channel estimate from the start.

Not yet done: RF25 (SP12_4, where a single symbol's comb is too sparse and the estimate
must sit on the smoothed grid), and a live A/B on air.

### 13.6 M17 — Subframe 1 / PLP 1 decoded live (2026-11-07, 00:25–00:58)

What OnAir does that we did not: decode every Subframe live. Done as `lab/m17_plp1_live.py`
(`SubframeDecoder`: m13's geometry and HTI reading, the live decoder's threaded CPE / demap /
LDPC-ladder / BCH batches, the E104 estimate as its channel estimate) plus `--plp1` in
`m11_watch` (whole-Frame window and 26 MB shared slots, the second decoder in both decode
paths, PLP 1's own Baseband stream tagged for the transport) and one ALP walker per PLP in
`m11_stream.Transport`. The staged E47 ROUTE lanes then light up unchanged.

**Gates.**
1. Byte identity: `python lab/m17_plp1_live.py hit_rf33.cs16 --frames 2` — the Baseband
   stream equals m13's offline stream byte for byte, estimator off and on (1.1 s a Frame
   against m13's 13 s).
2. Offline over the whole race capture (`--frames 81 --no-ref`): 117/117 on every Frame,
   both estimators.
3. Live (`watch --capture live_rf33_antB.cs16 --plp1 --route 239.255.32.1:8321 --assets all`):
   see the table.

**What broke, and the finding.** The first live run lost PLP 1 from Frame ~50 to the end:
linear interpolation fell to 0/117 for the last 27 Frames, the E104 estimate held 30–116.
PLP 0 never dropped a block. Offline the same Frames were 117/117. Cause: the front end
fixes the carrier offset once per lock (`m11_stream`, "the CFO is FIXED for the whole
lock"), and the tuner drifted +471 → +452 Hz over the 20 s. 19 Hz is 4.5% of a 16K carrier
spacing — inter-carrier interference near −22 dB, fatal for 256QAM — and 2.2% of an 8K
spacing, which 64QAM never notices. Fix (M17 only; PLP 0's path and its byte gates are
untouched): the residual is measured on the Subframe's own pilots as the common phase step
between same-parity symbols, the Subframe's samples are rotated by it before the FFT, and
`cfo_resid_hz` is reported per Frame. The same drift costs PLP 0 margin silently
(−28 dB ICI floor); tracking it in the front end is a separate, gated change.

| live on the race capture, 80 Frames | PLP 1 blocks | PLP 0 blocks | WHUT video lane |
|---|---|---|---|
| linear estimate, CFO fixed at lock | 5920/9009 (0/117 from Frame 55 on) | 100% | — |
| E104 estimate, CFO fixed at lock | 8230/8892 | 100% | 12 s |
| **E104 estimate + residual CFO tracked** | **9360/9360** | **5920/5920** | **18 s, 1080 frames, no gap** |
| OnAir on the same air (§13.2) | 9046/9360, two bursts | 5920/5920 | 4 s hole |

Speed, this box, lossless replay: ~1.1 s a Frame for Subframe 1 (demap 570 ms, LDPC
200–400 ms, pool+CPE ~200 ms) against the 247 ms Frame. Live air needs a compiled demapper
(the `cellpool_kernel.c` pattern) and a third decode process; the m16 numpy demapper is
slower still here (140 ms a Block). Instruments: `ATSC3_PLP1_LOG=1` (per-Frame PLP 1
blocks), `ATSC3_WIN_LOG=1` (window sizes handed to the workers).

### 13.7 The speed pass — E106 (demapper) and E107 (layered LDPC), 2026-11-07 01:00–02:00

Subframe 1 cost ~1.1 s a Frame against the 247 ms budget; two stages were the bulk of it.

**E106 — the demapper, vectorised.** `demap_kernel.c` (E102) was scalar C, built without
`-march`, with a branchy inner loop the compiler could not vectorise. Added an AVX2 path in the
same file (eight constellation points per vector, branch-free subset minima, chosen at run
time by `__builtin_cpu_supports`, the scalar path kept as the fallback). Per lane it is the
scalar path's arithmetic in the scalar path's association, and a minimum is the same number
whatever order it is taken in, so the output is bit-identical: `gate_e102_demap.py` passes
unchanged and M17 still equals m13 byte for byte. **570 ms → 66 ms a Frame** on PLP 1; the
same kernel serves PLP 0. (Our m16 numpy GEMM demapper was 4x slower than E102 on this box —
140 ms a Block — because numpy's strided minima do not behave here as they did where it was
measured; it stays as the exact reference it always was.)

**E107 — the layered, eight-lane LDPC kernel** (OnAir's idea; §13.3 item 3). A new kernel,
`ldpc_layered_kernel.c`, beside the flooding one: the same normalised min-sum rule, alpha and
first-minimum tie rule, the syndrome checked before the first update, but each check's
messages enter the totals at once (layered), and eight codewords ride the eight lanes of one
AVX2 vector. Python dispatch `m3_ldpc.layered_decode_batch`, hooked into `m9_fast._ldpc_cpu`
for the cpu_fast path only (`ATSC3_LDPC_LAYERED=0` restores flooding; the exact float64 path
never routes through it). Gate `gate_e107_ldpc.py`, E54's boundary on real air with seeded
noise, both codes (PLP 0 16200-bit through the live decoder, PLP 1 64800-bit through M17),
six noise levels, two Frames:

| | result |
|---|---|
| blocks converged under flooding but not layered | **0** |
| blocks converged under layered but not flooding | **28** (all at 21 dB: 73→96 and 109→114 of 117) |
| Baseband byte mismatches on shared blocks | 0 |
| non-deterministic blocks | 0 |
| mean iterations, PLP 1 (clean / 23 / 22 / 21.5 dB) | 3.5 / 11 / 17 / 24 → 2.2 / 6 / 10 / 12 |
| LDPC time over the gate | 22.5 s → 12.1 s (1.9x); clean-air PLP 1 134–381 → 74–79 ms |

So it is faster AND converges more at the cliff; on by default for cpu_fast.

*Amended 08:40, after the gate was made to reset the M17 decoder's per-Frame state between
its three decodes (the fed-forward carrier offset and the cached delay had made the "same
input" slightly different each time, which showed up as one "non-deterministic" block):*
with identical inputs the kernel is deterministic, and the strict count at the cliff is
**17 gained, 1 lost** on one Frame at 21 dB (92 vs 76 of 117), 5 gained on the other, 15 and
29 gained at 20.5 dB where flooding has 0. A change of schedule converges a slightly
different set near the threshold; a strict superset is not a law, and the gate's rule is now
net gain at every level (plus byte identity on shared blocks and determinism). The claim
"never loses a block" above is withdrawn in favour of "+16 to +29 a Frame at the cliff, one
block in ~900 the other way".

**Two more trims in M17 (07:00–07:15).** The residual carrier offset is fed FORWARD: last
Frame's value rotates the samples before the first transform, and the second transform runs
only when the residual measured afterwards is still over a hertz (a jump), so in steady drift
the pool pays one FFT pass, not two (the gate's second Frame: −5.83 Hz inherited, −0.44 Hz
left). The comb estimate runs in four batched transforms (symbols grouped by pilot parity)
instead of seventy-five, and `comb_denoise_grid` rotates by one vector per distinct
main-path delay instead of one per symbol (the per-symbol exponentials were ~100 ms a Frame
at 16K). M17 still equals m13 byte for byte. Standalone PLP 1 Frame now ≈ 490 ms (pool
166–188 with the comb / 95–122 without, CPE 90, demap 66, LDPC 80, BCH 66).

**Where the Frame stands after both:** PLP 1 ≈ 520–590 ms (pool 190 with the comb estimate and
the CFO second transform, CPE 90, demap 66, LDPC 75–150, BCH 68, rest ~35) plus PLP 0 ≈ 130 ms,
per decode process. Three processes put the effective Frame time at the budget on this 12-core
box; the real-time replay numbers are below. The next fat, if needed: the BCH syndrome (68 ms
for 117 x 64800 bits) and the pool's second FFT pass.

### 13.8 The three-way race (2026-11-07, 07:45) — `lab/race3.py`, `lab/e104_out/race3/`

Same 20 s of RF33 air (`~/captures/race_2026-11-06/live_rf33_antB.cs16`), same 12-core PC,
one lane after another, nothing else running.

| lane | what | wall for 20 s of air (lossless) | PLP 0 blocks | PLP 1 blocks | WHUT 32.1 video | WJLA 7.1 video |
|---|---|---|---|---|---|---|
| A | OnAir `atsc3rx`, 12 threads | 18.1 s | 5920/5920 | 9046/9360 | 720 frames, 16 s, **1 gap** | none (no MMTP) |
| B | ours, before 2026-11-06 (worktree at HEAD, original scalar kernels) | **9.8 s** | 5994/5994 | not decoded | none | 1080 frames, 18 s, 0 gaps |
| C | ours with their ideas (E104 + M17 + E107, plus E106) | 37.1 s | 5920/5920 | **9360/9360** | **1080 frames, 18 s, 0 gaps** | 1080 frames, 18 s, 0 gaps |
| D | C with the switchable ideas off (`--no-plp1-comb`, flooding LDPC) | 40.7 s | 5920/5920 | 9360/9360 | 1080 frames, 0 gaps | 1080 frames, 0 gaps |

(B decodes 81 Frames and C/D 80: C/D's whole-Frame window costs the last partial Frame of
the file. B's 5994 vs 5920 is that extra Frame.)

**Real-time headroom** (`--realtime --shed`): B, 3 processes: 1.19x, 77 of 81 Frames, WJLA
945 frames with one gap. C, 3 processes x 16 threads (the harness default, which thrashes):
0.77x, 42 Frames, no ROUTE lane written. C, 4 processes x 3 threads (`C_rt_4x3`): 0.91x, 47 of
81 Frames, WHUT 480 frames (8 s, no gap) and WJLA 480 frames (18 s, one gap). Under shedding
a ROUTE lane keeps only the DASH segments whose every Frame arrived (a segment is a whole
object), so PLP 1 on live air wants the full budget; PLP 0's MMTP lane degrades more
gracefully.

**PLP 1 cliff, seeded noise, two Frames** (blocks of 234; idea 1 = E104 estimate, idea 3 =
E107 layered LDPC; M17 live decoder throughout):

| noise below signal | both on | E104 only | E107 only | both off |
|---|---|---|---|---|
| 23 dB | 234 | 234 | 234 | 234 |
| 22 dB | 234 | 234 | 234 | 234 |
| 21.5 dB | **234** | 234 | 145 | 78 |
| 21 dB | **225** | 219 | 0 | 0 |
| 20.5 dB | **28** | 1 | 0 | 0 |

Read: the E104 estimate is worth about a dB at the cliff on this signal (0 → 219–225 blocks at
21 dB); the layered LDPC adds a few tenths on top (219 → 225 at 21 dB, 1 → 28 at 20.5) and
halves the decode time. Idea 2 (PLP 1 live) has no "off": without it the WHUT column is empty,
as lane B shows.

**Verdict.** Against OnAir on the same air: more blocks (9360 vs 9046), the whole programme
without a hole (1080 vs 720 frames), with sound (AC-4) and the MMTP station beside it — but
2x their wall time for the two PLPs on this box, and not yet real time with both PLPs
(0.90x; PLP 0 alone 1.19x). Against ourselves: the ideas cost wall time (two PLPs instead of
one) and bought a new station, a dB of margin and a faster LDPC.

### 13.9 Real time with both PLPs — E108 and the replay clock (2026-11-07, 07:50–08:30)

Asked first: is our code slower now with their ideas? Like for like (today's tree, `--plp1`
off, the job lane B did): **8.7 s lossless (4.4x) against B's 9.8 s (3.7x); real time 1.22x,
79 Frames, against B's 1.19x, 77.** About 12% faster on the same work. The race's 37 s was
the second PLP, four times the bits.

Then the budget. Single-thread CPU per Subframe 1 Frame after E106/E107 was ~0.6 s, so the
box had the capacity twice over; what did not scale was the Python inside each decode
process. Four changes, each gated (M17 == m13 byte for byte after every one):

* **E108a — the equaliser rows in C.** Each comb symbol's divide-and-scatter goes through
  `cellpool_row_sm_f64` with an identity grid map (the E104 trick in `m9_fast`); the kernel's
  division mirrors numpy's bit for bit. Linear-fallback symbols keep the Python path.
* **E108b — the BCH syndrome in C** (`bch_kernel.c`, `build_bch_kernel.py`, hooked in
  `m9_accel.bch_syndrome_batch`, `ATSC3_BCH_KERNEL=0` restores numpy): the same popcounts
  over the same masks, integer arithmetic, equal by construction and asserted. Raw call 11 ms
  against numpy's 59 on 117 x 64800 bits — once compiled for the POPCNT instruction; the
  generic builtin was SLOWER than numpy (244 ms), a lesson for every kernel here.
* **E108c — the comb estimate's first transform cached.** Its first pass only finds the main
  path's delay, which is static; `comb_denoise_grid(tau0=...)` takes last Frame's and checks
  it (a stale hint, strongest tap away from zero, is measured afresh). The four pilot groups
  run on the thread pool. Comb 85 → 33 ms.
* **E108d — `SubframeDecoder.prewarm()`** in both decode paths: one silent Frame builds the
  64800-bit BCH masks, the LDPC edge tables and the FFT plans before "ready" (3 s a worker).

Subframe 1 standalone: ~1.1 s a Frame at first light → **~360 ms** (pool 120, CPE 71, LDPC 80,
BCH 34, demap 27, rest ~30).

**The replay clock.** The first real-time runs decoded 48 of 81 Frames at a reported 1.0x
with 32 "dropped"; prewarming made it 23 of 81. Not throughput: when the chain sniffs the
multiplex it starts the replay source's real-time clock INSIDE the sniff, before the decode
pool warms, and the reader never re-anchored it — so the warm-up's worth of air arrived as a
back-dated burst that the queues shed. `FileSource.rebase()` re-anchors the pacing when the
reader takes over (pacing only; the sample stream is unchanged). Lane B's 77 of 81 in §13.8
was the same artefact.

| real-time replay, 4 processes x 3 threads, shedding allowed | Frames | dropped | PLP 0 | PLP 1 | WHUT | WJLA |
|---|---|---|---|---|---|---|
| before the clock fix, no prewarm | 48/81 | 32 | 100% | 100% | — | — |
| **final** | **80/81** | **0** | 5920/5920 | 9360/9360 | 1080 frames, no gap | 1080 frames, no gap |

Instantaneous rate 1.11–1.34x on this 12-core box with nothing else running. Both PLPs of RF33
are decodable live. Not yet tried: live air with the radio (the next session's first run).

**Final gates over this state (08:35):** E102 demap identity PASS; E104 self-check PASS (5 legs);
M17 == m13 byte for byte PASS; E107 PASS under the amended rule (§13.7); and a lossless replay
of the race capture whose WHUT video, WJLA video and WHUT audio lanes are **byte-identical** to
this morning's 00:58 run — every speed change since then left the chain's output untouched.

### 13.10 Live air — E109, and the first television from both PLPs (2026-11-07, 07:55–08:50)

RF33 off the RSPdx, Antenna B with the bias-T, the gain table's rfgain_sel 4 / IFGR 32,
`--plp1 --route 239.255.32.1:8321 --assets all`, 60 s runs, this 12-core PC.

**First run (4 x 3): it worked, and it starved itself.** Lock at +464 Hz, PLP 0 100% of what it
decoded, a clean 1080p frame of WHUT's Saturday-morning cartoon — but 138 Frames of ~243, 107
dropped, two SDR overflows, two re-acquisitions, 0.53x. Cause: `SdrSource.open()` set up AND
activated the stream inside `build()`, and the four workers then warmed for 19 s (each decoding
a silent Frame in `SubframeDecoder.prewarm`) while the radio's 1.3 s ring buffer overflowed.
The PLP 0 chain had the same wart at a smaller scale (2 s warm-up; the "startup overflow costs
a re-acquisition" comment in the source).

**E109.** (1) `prewarm()` builds tables only (BCH masks, LDPC edge tables, FFT plan): 0.2 s a
worker instead of 3–6. (2) `SdrSource.open()` configures but does not stream; `start()`
streams and `pause()` stops. The sniff starts and pauses the radio; the reader starts it again
once the pool is warm and pushes ONE deliberate continuity break, so the front end re-acquires
once on purpose instead of twice by overflow. (`FileSource` keeps its clock `rebase()`, §13.9.)

| live RF33, 60 s | Frames (of ~243) | dropped | PLP 0 | PLP 1 | WHUT lane | WJLA lane | overall / instantaneous |
|---|---|---|---|---|---|---|---|
| 4 x 3, before E109 | 138 | 107 | 10212/10212 | 15178/16146 | 2 s | 73 s, 10 gaps | 0.53x |
| 4 x 3, after | 227 | 16 | 16798/16798 | 26549/26559 | 40 s, no gap | 52 s, 8 gaps | 0.88x / 0.98–0.99x |
| **5 x 2** | **239** | **4** | **17686/17686** | **27945/27963** | **56 s, 2 gaps** | **58 s, 3 gaps** | 0.93x / 0.90–0.93x |
| 6 x 2 | 178 | 65 | 13162/13172 | 20779/20826 | 38 s | 51 s, 16 gaps | 0.68x |
| 4 x 2 | 137 | 106 | 10133/10138 | 15612/16029 | 58 s (720 frames) | 56 s | 0.53x |

On live air the main process also runs the radio reader and the front end, so the decode pool
gets ten threads, not twelve: **5 processes x 2 threads** is this box's live layout, 98.4% of
Frames, both stations' pictures nearly continuous. The replay's 4 x 3 oversubscribes by one.
A few Frames a minute still shed at the margin; the next CPU to recover is the CPE's per-symbol
Python (71 ms wall, ~118 CPU-ms a Frame) and the comb's remaining transforms.

The live run also makes the point the race could not: on air, OnAir's PLP 1 bursts (§13.2) and
our fixed-at-lock CFO loss (§13.6) are the same class of fault, and the chain now tracks the
drift it measured here (+464 Hz at lock, +447 Hz a second later).

### 13.11 E110 — the two-minimum min-sum correction, and the marginal-carrier win (2026-11-07, 09:00-09:30)

The third thing OnAir does that we did not: the two-minimum (box-plus) correction of the
min-sum check rule (their minsum_corr.h; the textbook Chen-Fossorier result). Every edge of
a check that does not hold the smallest input gets, instead of min1, the box-plus of the two
smallest, min1 (+) min2 = min1 + softplus(-(min1+min2)) - softplus(-(min2-min1)) -- closer to
the belief-propagation message, needing one number the check already has. Written from the
formula, not their code (NOTICE). In `ldpc_layered_kernel.c` as an AVX2 path (ABI 2, the
`corr` argument; softplus(-t) as the max of four tangents, which vectorises and is a lower
bound within ~0.02); switch `ATSC3_LDPC_BOXPLUS` / `fd.ldpc_boxplus`, OFF by default.

**Gate (`gate_e107_ldpc.py`, the `box` column):** against the layered kernel, byte-identical
on shared blocks, net gain at every noise level, deterministic. At 21 dB of added noise:
layered 92/117 and 114/117 -> box 110/117 and 117/117. At 20.5 dB: layered 15 and 29 -> box
32 and 58 (roughly doubles). On clean air it is free (no extra iterations); it only costs
time on frames already failing (7 ms clean, up to ~2.5 s a frame grinding at the deepest
cliff).

**On the air it matters most where it is weakest.** Banked 45 s of RF25 (the market's
weakest ATSC 3.0 carrier, bootstrap peak 77, -34 dBFS, LDM+CTI, PLP0 QPSK 9/15), replayed
losslessly three ways on identical samples:

| setting | CTI commutator | core blocks | video |
|---|---|---|---|
| flooding min-sum | never locks (0/35 clean on every candidate) | 86/2013 | none |
| layered | never locks (0/35 clean) | 115/2013 | none |
| **layered + E110 correction** | **locks** | **1426/5296** | **a lane written** |

The CTI phase search confirms a candidate start_row only on 4 clean FEC Blocks. At this level
min-sum decodes 0/35 on every candidate and never locks; the correction decodes just enough
more that a candidate passes, the interleaver locks, and the decode runs away. The only
difference between the layered and box runs is the one flag, so the lock is the correction's.
This is the marginal-signal payoff OnAir's notes pointed at -- larger than the ~0.2 dB it is
worth on strong carriers, because it crosses the threshold the lock logic needs.

**Why still OFF by default.** On a live marginal signal the per-frame cost (up to ~2.5 s
grinding) can blow the real-time budget and shed frames -- the opposite of help -- where the
layered kernel is always faster. The clean fix is the E52 shed treatment the SP rescue
already has (time-box the extra iterations, shed beyond, count), after which it can go on by
default. Today it earns the offline/marginal-recovery win above; the deadline-box is the
follow-up. Also still staged: the three ATSC 3.0 carriers here (RF25/RF30 LDM+CTI, RF33 not)
and a quiet-box redo of the live stress timing (the 09:00 stress run's rate is contaminated
by concurrent builds; its per-frame convergence counts stand).

### 13.12 LDM speed — the profile that sets the plan (2026-11-07, ~10:00)

RF25 (marginal LDM, 189 symbols/frame) decodes with the E110 correction but at ~0.3x real
time. Profiling the OFDM demod directly (`core_cells` on the banked capture, lab/e104_out
has the script) per frame: fft 82 ms, equalise 52 ms, cpe 50 ms, ce_smooth 10 ms = 193 ms
(budget 242). Thread scaling 1/2/4/6/10 = 275/225/198/196/195 ms: it PLATEAUS at ~4 threads.
The equalise and cpe are numpy loops that hold the GIL; scipy's FFT releases it and scales a
little. So:

* **Threads do not scale one channel's demod (GIL).** The 32-core's lever is more decode
  PROCESSES (the demod is frame-independent; m44's LdmDemodPool already spawns worker
  processes), not threads. The CTI + FEC stay serial.
* **C kernels for equalise and cpe are a double win** -- faster AND they release the GIL, so
  threads rescale. The `cellpool_kernel.c` equaliser and `demap_kernel.c` cpe_nearest both
  exist on the m9_fast path and are not wired into the m44 demod. This is the next build.
* **FFT (82 ms) is the GPU's job** (batched 8K cuFFT); it already releases the GIL and
  scales. Parked until the GPU is sorted (Pascal GTX 1060 + Python 3.14 torch has no working
  wheel) or tested on the bigger box.

Honest bottom line for this 6-core box: the C kernels buy a worthwhile ~1.4x but not
real-time; real-time for an LDM channel needs more cores (the 32-core, via processes) or the
GPU for the FFT. The decoding itself is solved (E110); this is purely throughput.

### 13.13 E111 — the compiled equaliser and nearest-point search on the LDM path (2026-11-07)

The E103/E106 equaliser (`cellpool_kernel.c`) and the E105 nearest-point search
(`demap_kernel.c`), both already bit-exact on the m9_fast path, wired into the LDM
demod (`m44_ldm.FrameDemod`): `cellpool_row_sm_f64` / `cellpool_row_f64` per row in
`cell_pool`, `cpe_nearest_f32` in `cpe_fast`.  Off: `ATSC3_LDM_CELLPOOL_KERNEL=0`,
`ATSC3_LDM_CPE_KERNEL=0`.  Referee: `lab/gate_e111_ldm_kernels.py` on RF25 (189
symbols, every row on the E60 smoothed estimate).

**Measured, per frame (8 threads):** equalise 52.9 -> 29.6 ms (1.79x), cpe 53.4 -> 48.3 ms
(1.11x), core_cells 192 -> 163 ms (1.18x).  Owner map, coherence, cpe output,
core_cells (complex64) and determinism bit-identical on every frame.  The cpe win is
small because the search was a minority of `cpe_fast`: the numpy reductions over the
1.2 M-cell pool (two full passes, three reduceat/repeat rounds) are the rest.

**Finding -- the pool is NOT bit-identical on a smoothed row, and it is numpy, not the
kernel.**  ~6 ULP of full scale in a third of the cells.  Isolated: interpolation
bit-exact, division bit-exact (0/6026 with the same H); the `H *= cl` complex multiply
is the whole difference (2263/6625 of one row).  A scalar plain-double reference agrees
with the kernel in 0 of 23,432 cells and disagrees with numpy -- numpy 2.5's complex128
multiply reproduces neither the textbook two-rounding form nor an FMA one, which the
`derot_f64` note in the same file had already measured.  So the gate's law is split, and
is derot's: bit-identical wherever no gain multiply is involved (preamble cells, every
un-smoothed row -- measured 0 diffs), within 1e-14 of full scale on smoothed rows (a few
hundred ULP; a real error is O(1), thirteen orders away), bit-identity on everything
downstream.  Recorded in `cellpool_row_sm_f64`'s comment.

**Thread scaling, kernels on (1/2/4/6/10 threads): 177/163/164/163/164 ms** against
275/225/198/196/195 off.  The "C drops the GIL so threads rescale" claim of 13.12 is
retracted in its strong form: the equaliser does drop it, but it is 30 ms of the frame;
the plateau -- now reached at 2 threads instead of 4 -- is set by the FFT (82 ms, on
scipy's own workers) and the numpy CPE passes (48 ms, GIL-held), ~130 ms threads cannot
touch.  At ONE thread the kernels are worth 1.55x, which is the number that matters for
a many-process layout.  Next walls, in order: the FFT (float32 fast path with a
decoded-bytes gate, or the GPU) and fusing the CPE's numpy passes into one C pass.

Checkpoint committed (e40abef) and pushed before the replay referee ran -- the rule from
here on: save before any new C runs.

**Replay referee (lossless RF25, 2 procs x 4 threads): PASS.**  Every decoded lane
byte-identical.  `live.json`, the one file that differs, differs in exactly 9 of its 70
fields -- seven output paths (the two run directories) and the two wall-clock stamps --
and in no decoded or statistical field.  Core 1426/5296 in both runs, 184 frames, 0
dropped; the on-run log carries no fallback message, so the kernels ran inside the
spawned demod workers (E86c: the lever engaged, and says so in `rep["cellpool_kernel"]`).
Pipeline: 0.309 -> 0.322x real time, +4% only.  Amdahl on a saturated six-core: the demod
is one share of the CPU, and `fec_dispatch` (389 -> 405 ms) names the other -- on this
marginal carrier 73% of blocks never converge and each burns all 50 LDPC iterations in
the FEC pool the main thread then waits on.  That hidden cost is the next lever (E112:
early abort of hopeless blocks, under the E107 net-gain rule), ahead of the FFT.

### 13.14 E112 — early abort of hopeless LDPC blocks (2026-11-07)

**Why.** 13.13's replay named the cost E111 could not reach: `fec_dispatch` 389 -> 405 ms
of a 1258 ms LDM frame is the main thread waiting on the FEC pool.  Sized with a per-block
record (`ATSC3_LDPC_DUMP`: every call's iterations-to-converge and final syndrome weight,
every 8th call's LLR block kept) over the RF25 replay: **94% of all LDPC work -- 661,050 of
699,750 iterations -- is burned by blocks that never converge** (9.7% converge).  Blocks
recorded 14,647; the CTI phase search re-submits, so more than the 5,296 counted.

**The shape of a hopeless block, measured.**  Re-decoding 1,829 dumped 64800-bit blocks
(alpha 0.85, nc = 25,920 checks) at fifteen iteration caps with the layered kernel:

* A block that never converges is **flat from sweep 1** -- syndrome weight p50 10,075 at
  cap 1 and 10,088 at cap 50; p5 above 9,350 from sweep 4.  It is not slowly failing; it
  was never going to decode.
* A block that converges **descends**: p50 7047 (sweep 1) -> 3047 (4) -> 343 (8) -> 11
  (12) -> 0.  Half are done by sweep 25, 88% by 40, the last at 50 (so the cap itself is
  not the lever; the cap is right).
* The two populations are separated in LEVEL: at sweep 4 the worst eventually-converging
  block is at 8,208; the 5th-percentile hopeless block at 9,352.  A threshold between them
  is a zero-loss rule.

**Rules, evaluated on the trajectories (ground truth = the kernel's own 50-sweep result):**
threshold "at sweep k, give up if syndrome >= T": k=4, T=8563 aborts 1,721 of 1,739 hopeless
blocks, saves 89.6% of LDPC work, loses 0; every zero-loss row from k=4 to 20 saves 58-90%.
Stall "give up if no improvement since sweep k-m": at most 39% saved, and k=6,m=2 LOSES two
blocks -- a hopeless block jitters about its flat level, so "did not improve" fires on only
half of them.  Threshold is the rule.  Shipping choice **k=8, T=0.36 x nc** (9,331 here):
14% above the worst converging block at sweep 8, ~80% saved; margin bought with savings,
because 90 converging blocks on one carrier is a thin sample and the referee below runs on
two other codes.

**Mechanism.**  `ldpc_layered_kernel.c` ABI 3: `abort_k, abort_t, abort_stall` after `corr`
(0 disables each; all zero is ABI 2, bit-identical).  A lane that fires keeps its hard
decisions, conv=0, its=the sweep it stopped at, bad=its syndrome -- a failed block's
meaning.  `ATSC3_LDPC_ABORT="k,T,stall"` (T<1 = fraction of the checks), `layered_decode_batch(abort=)`,
`fd.ldpc_abort`.  Committed (0b303e6) and pushed BEFORE the first build of the new C: the
rule from here on.  Referee: `gate_e107_ldpc.py --abort` -- an abort can only lose, so the
law is stricter than E107's net gain: **zero blocks lost at every noise level**, bytes
identical, faster; then the RF25 replay must return exactly the 1,426 blocks with
byte-identical lanes.

**Verdict.**  Gate (RF33, 16200-bit 64QAM 11/15 and 64800-bit 256QAM 11/15; as received and
23/22/21.5/21/20 dB seeded): **abort loses no block at any level** -- at the 21.0 dB cliff,
where the layered kernel converges 92/117 and 114/117, abort converges 92/117 and 114/117;
bytes identical on every shared block; LDPC time over the gate flooding 22,822 ms, layered
12,385, layered+abort **6,621 ms (1.87x over layered, 3.4x over flooding)**; at 20 dB, where
nothing converges, 2,627 -> 1,211 ms.  RF25 replay, abort on: **core 1426/5296, identical**,
184 frames, 0 dropped, 6 of 6 lanes byte-identical to the E111 run, no fallback; decode
1233 -> 805 ms/frame (`fec_dispatch` 406 -> 219; `pool_wait` 669 -> 465, because the freed
FEC pool also unblocks the demod workers that share it); **x rt 0.322 -> 0.452, +40%, and
+46% over the day's start (0.309)**.  The largest single pipeline gain of the day, for a
fifteen-line change chosen by measurement.  Default ON from this commit
(`ATSC3_LDPC_ABORT=0,0,0` turns it off); the gate passes (0,0,0) explicitly for its
non-abort columns so they stay the pure kernel.  The frame is now 805 ms: pool_wait 465 is
the demod again (FFT 82 ms of a 163 ms core_cells), and that is the next wall.

### 13.15 E113 — the "FFT" that was a Python loop (2026-11-07)

**The mislabel.**  13.12 profiled `core_cells` as fft 82 / equalise 52 / cpe 50 / ce_smooth
10 ms and called the FFT the next wall, a GPU job.  Sizing it before building anything: a
(189 x 8192) complex128 FFT is **14.8 ms**, flat from 1 to 8 workers (this batch does not
parallelise); the strided-to-contiguous copy 3.3 ms; `fftshift` a 25 MB copy, 3.5 ms; the
whole fft+shift stage **18.3 ms**.  complex64 would buy 2.7 ms and break bit-exactness --
the float32 FFT path is dead.  So where were the other 64 ms?  `cell_pool`'s "fft" timer
starts at the top of the function and also covers the PREAMBLE loop: `_demod_preamble_g`
(2.8 ms per symbol) and **`FI.deinterleave` -- 30 ms per symbol, two symbols, 60 ms per
Frame.**  The timer is now split so this cannot hide again.

**The cause.**  `m3_freqint.interleaving_sequence` -- the A/322 7.3 frequency-interleaver
address generator -- is a pure-Python LFSR loop over nmax (8192) iterations, and it was
regenerated on every call.  The sequence is FRAME-CONSTANT: a function of (nfft, l,
n_data, direction, toggle) only.  The data symbols' sequences were already built once in
`tables()` (`Hmat`); the preamble path alone called it fresh each Frame.  60 ms of a 163 ms
`core_cells`, in a loop that computes the same 5,395 integers every 242 ms.

**The fix, E113.**  `interleaving_sequence` is an `lru_cache` over the generator, and every
caller receives its OWN COPY (45 KB, 0.026 ms measured against 30 ms) -- so no caller can
observe another's array and the cached call is bit-identical to the uncached one by
construction (checked bitwise against the raw generator).  One place, every caller: m9_fast,
M17, m13 and m44 all go through it.  Committed and pushed before the referees ran.

**The CPE, sized while there (48 ms in situ):** normalise 9.7, decimate+unique 5.5, the C
search 10.8, conj-mult+reduceat 13.8, repeat-divide 4.0, `unique(owner)` 8.3, final
repeat-divide 7.6.  The two `np.unique(owner)` calls (~14 ms) are frame-constant too -- a
bit-exact follow-up; the two full-pool passes (~17 ms) fuse into one C pass but change
rounding, so they need the decoded-bytes gate.  Neither is urgent beside E113.

**Verdict.**  gate_e111 PASS, pool bit-identical (the preamble cells are in it);
`core_cells` 134 -> 103 ms with the kernels on, against 192 -> 163 before E113 -- the
predicted 60 ms off both columns.  gate_e107 PASS, bytes identical on both PLPs (RF33's
path goes through the same cache).  RF25 replay: core 1426/5296 identical, 6 of 6 lanes
byte-identical, x rt 0.452 -> 0.462: **+2% only.**

**The finding that matters: `pool_wait` did not move -- 465 -> 464 ms -- while the demod it
waits on got 37% faster.**  The pool's wall time per Frame is no longer the demod's
arithmetic.  It is the plumbing: 2 x 4 threads on six cores (the demod plateaus at 2
threads, 13.13, so the other two are contention), demod and FEC sharing the same worker
processes, the 25 MB window and the 19 MB pool copied through shared memory every Frame,
and the serial main thread (`fec_dispatch` 253 ms) that must drain a slot before the pool
can take the next window.  Per-Frame math is now cheap enough -- demod 192 -> 103 ms,
LDPC 1.9x -- that further kernels will not move this box's real-time rate.  The lever is
the layout and the structure.  Next measurement, zero code: a procs x threads sweep of the
replay.  The CPE's frame-constant `unique(owner)` (~14 ms, bit-exact) is demoted for the
same reason.

### 13.16 The layout sweep, the frames that were "lost", and what actually scales (2026-11-07)

**The sweep.**  13.15 ended on "the lever is the layout".  `--ldm-procs x --ldm-proc-threads`
on the RF25 replay (the LDM pool's own flags; `--decode-procs`/`--threads` are the RF33
path's and do not reach it -- this box defaults to 2 x 4): 2x4 0.461, 2x2 0.398, 3x2 0.477,
4x1 0.486, 4x2 0.492, 5x1 **0.652** -- and not one of those is a throughput number, because
the block counts moved with the layout on a LOSSLESS replay (1426/5296 at 2x4 and 2x2,
1550/5433 at 3x2, 1329/5050 at 4x1, 708/3308 at 5x1) and so did the decoded lanes (2x2 and
3x2 byte-identical to 2x4; 4x1, 4x2 and 5x1 not).  5x1 was "fast" by doing 63% of the work.

**The wrong fix, measured and reverted.**  `late_frame_dropped` was 7 of 184 Frames even at
2x4 and 11 at 5x1, and CtiStream.push zero-fills a gap then drops the real Frame when it
arrives -- so E114 fed the CTI through a reorder buffer.  The referee was flat: identical
counts, `gap_filled` 0 on every layout.  Traced (`ATSC3_CTI_LOG`): the feed is strictly in
Frame order -- the buffer never had anything to do.  Reverted; a lever that silently no-ops
is worse than one that is off (E86c).

**What the "late" Frames are.**  Every one is `idx = origin - 1` immediately after
"CTI phase candidate: Frame N": reset() sets the origin to the Frame whose L1 was just read
while Frame N-1 is still in flight, and it lands behind the new origin.  Frame 5 -> 4,
10 -> 9, 146 -> 145, 151 -> 150, 156 -> 155, 161 -> 160, 166 -> 165: seven adoptions, seven
"late" Frames.  They are PRE-ORIGIN Frames the new stream does not want.  Benign; the
counter is misnamed, not the data lost.

**What actually scales with the layout: re-acquisition cadence.**  On this marginal signal a
candidate is adopted, judged over 35 Blocks, rejected (0/35 clean), the CTI goes cold --
and the next fresh L1 is read FIVE Frames later (146, 151, 156, 161, 166), because the
pool already holds nslots = nproc + 2 = 4 submitted Frames when it goes cold and the L1
is read from the next window to enter.  At 5 procs that is 7 in flight, ~8 Frames per
rejected candidate instead of ~5, and a different set of Frames gets its L1 read, so the
search walks a different path and locks at a different time: hundreds of Blocks either way
on a 45 s capture, and 3x2's 1550 was the cycle landing on a better candidate.  The
consequence is general and matters for the bigger box, whose default is 6 processes:
**on a marginal carrier, more processes re-acquire SLOWER.**  The old E85 scaling note
(4 procs 1.49x, 6 procs 2.01x) was measured on a locked signal and is not contradicted,
but it does not cover this.  Decoupling the L1 re-read from the pipeline depth -- reading it
from the most recent completed window the moment the CTI goes cold -- is the next build,
and until it lands the layout sweep is only meaningful after lock.

### 13.17 E115 — re-acquire from the held window; the layout sweep, honestly (2026-11-07)

**Mechanism.**  `push_frame` keeps a copy of each window until the CTI has consumed its
Frame (at most the in-flight Frames, nslots x ~27 MB); at either FEC-verdict cold site
(`every candidate rejected`, `out of range too`) `_reacquire_now` reads L1 from the oldest
unconsumed window and anchors there.  Not hooked at the range-gate cold inside `_begin`
(same window, same candidates: it would only recurse).  Counter `ldm.stats.reacquired_from_window`.

**Referee (RF25 lossless replay, 2x4 / 3x2 / 5x1).**  The lever engaged (14 / 14 / 13
re-acquisitions).  **Decoded lanes byte-identical across all three layouts** (6 of 6 each;
before E115, 5x1 matched 2 of 4), and 2x4 and 3x2 agree block for block (1571/5763).
Blocks 1426 -> 1571 on the default layout: the in-flight Frames are kept.  Lane bytes
identical in size to before -- on this capture the kept Frames added accounting, not picture.

**The cost, attributed.**  x rt 0.461 -> 0.420 on 2x4.  Not the copies, not the L1 reads
(`l1` 5.4 -> 8.1 ms/frame); it is `fec` 43 -> 77 ms/frame, because **candidate attempts
went 9 -> 16**: a rejected candidate used to be followed by ~5 Frames of forced waiting, and
re-acquiring at once runs the adopt-judge-reject cycle nearly twice as fast on a carrier
that mostly fails to lock, ~35 Blocks of FEC per attempt.  So: on a carrier that locks and
stays locked, E115 costs nothing (no cold events); where a rejection is followed by a
successful re-acquire it saves ~5 Frames (1.2 s) of blackout per event, the live-TV win;
where nothing locks it spends ~9% trying harder.  Kept ON, because the output no longer
depends on the process layout -- the correctness point, and what makes the bigger box's
six-process default behave like this box's two.  Follow-up if live testing shows the
treadmill hurting: throttle held-window re-acquisition to one per cold event.

**The layout sweep, now meaningful.**  2x4 0.420, 3x2 0.434, 5x1 0.433 -- within 3%.  On
this six-core box the layout barely matters: 5x1 does cut `pool_wait` (458 -> 365 ms) but
the main thread's serial work (`fec` 77 + `fec_dispatch` 244 + CTI gather) absorbs it --
the cap the E85 note found "past 6 processes", seen from here.  No default change for this
box class.  The next structural lever for throughput is that serial main thread; the
32-core test will show where it caps.

**Day's ledger, RF25 on this box:** 0.309 -> 0.420 x real time with layout-independent
output (0.461 before paying for E115's determinism); demod 192 -> 103 ms; LDPC 1.9x;
every decoded byte identical throughout E111-E113; E114 reverted on evidence.

### 13.18 Live, RF25, the whole day's stack (2026-11-07 12:03)

Antenna B, 90 s, E111-E113 + E112 abort + E115 on, the E110 correction on (ATSC3_LDPC_BOXPLUS=1),
default 2x4 layout, against the morning's live baseline on the same carrier (pre-E111: 79
Frames, 1770/4400 blocks, 0.22x, no video lane):

**125 Frames (+58%), 2255/6922 blocks (+27%), 0.323x real time (+47%); phase-locked from
10.3 s and locked through the end (17 status lines); E115 re-acquired 3 times across 4
candidate adoptions -- no treadmill on air, so no throttle on this evidence.**  Still not
watchable: at 0.32x the front end shed 246 Frames and no complete video object was
assembled.  The live frame is 1182 ms -- pool_wait 661 + fec_dispatch 384 = 88% waiting,
worse than the 833 ms replay because the front end, tracker and reader share the same six
cores with the pool.  On this box the LDM carrier now decodes and holds lock on air; what it
lacks is cores, and the serial main thread will be the cap once it has them.  That is the
32-core test: `--ldm-procs` 6 / 8 / 12 / 16 with `--ldm-proc-threads 2`, block count and
lanes unchanged by the layout (E115), and the rate's ceiling named.

### 13.19 The GPU path, refereed (2026-11-07, torch 2.14.1+cu126 on the GTX 1060)

**torch.**  Python 3.14 + Pascal has exactly one working wheel: 2.14.1+cu126 (2.15 is built
cu132 and dropped Pascal; the cu126 index stops at 2.14.1).  Real kernel on the 1060: 3,540
GFLOP/s.  The wheel is kept at `~/.cache/torch_wheel/` -- the download drops constantly --
and the install needs `TMPDIR` on /home, not the RAM-backed /tmp.

**Referee: RF33 lossless replay (`live_rf33_antB.cs16`, 20 s), `--accel cpu` vs the GPU modes.**
cpu: 80 Frames, 5920/5920 + 9360/9360, 1.03x.  `gpu-full`: 46 Frames before the 300 s
timeout, 0.04x.  `gpu` (LDPC only): 36 Frames, 0.03x; with `--decode-procs 5`: 35, 0.03x.
**Every GPU lane is an exact byte-prefix of the CPU lane** (video 3,128,687 bytes identical
on gpu-full; 4 of 4 lanes on both gpu legs) and every block it touched converged: the GPU
LDPC is bit-correct.  Engaged, by its own timers (`gpu_bicm`).  **On its own stage it is
fast: ldpc 80 ms vs the CPU's 225 (2.8x).**

**Why 0.03x anyway -- a mis-wiring, not the GPU.**  Beside the GPU LDPC, every CPU stage ran
5-10x slower than in cpu mode: cell_pool 470 vs 44, demap 494 vs 101, cpe_correct 490 vs
182, the front end 763 vs 39 (notch 466 vs 0).  Those are exactly the stages the E102-E108
kernels serve, and `m11_watch` gates them off: `cpu_fast = (accel == "cpu" and not
exact_cpu)` at lines 397, 963, 1020 and 1526 -- so a GPU accel runs the HEAD-exact float64
path with no kernels at all -- and `use_procs` requires `accel == "cpu"` (974), so the
explicit five processes were ignored.  The GPU modes also print no mode line (E86c).  The
only tracebacks, one per timed-out run, are shutdown races at the kill (`cannot schedule new
futures after shutdown`, `write to closed file`): benign.

**What it would take, and what it would pay.**  E116: let `cpu_fast` depend only on
`--exact-cpu`, let the GPU backend replace only the LDPC (and CPE), allow processes in gpu
mode, and gate it as E107 was gated -- the GPU is a flooding min-sum against the CPU's
layered kernel, so the law is net gain >= 0 at every seeded level and bytes identical on
shared blocks.  Payoff on this box: RF33 ldpc 225 -> 80 of a ~1320 ms Frame, ~11%, which is
headroom for both PLPs live.  **Not an LDM lever here:** the LDM LDPC is 60 ms with the
E112 abort, which the GPU min-sum does not have, so on a marginal carrier it would run all
50 sweeps on every hopeless block and lose.  The GPU is a modest RF33 / multi-channel lever
on this box, nothing more; the LDM wall remains cores and the serial main thread.

### 13.20 The 32-core box (2026-11-07, Windows 11, Threadripper 32c/64t, MSYS2 gcc 15.2, Python 3.12)

The question 13.18 set: with cores, does one LDM channel reach real time, and what caps
it.  Setup as INSTALL says; all seven kernels built and smoke-tested with the MSYS2 gcc the
build scripts find on their own.  The Linux box's marginal capture never arrived, so RF25
was banked HERE on the RSPdx (`tools/atsc3_capture.py`, Antenna B, rfgain_sel 2, IFGR 32,
bias-T on, 0 overflows, no clipping, `C:\captures\rf25_win_45s_b.cs16`), and it is a
STRONG RF25: 80% of blocks converge as recorded.  The marginal regime was then made from
it with `lab/add_noise_cs16.py` (seeded complex Gaussian, stated dB below the recorded
power): 20 dB = 74% of blocks, 15 dB = 59%, **12 dB = 21% -- the Linux capture's 27%**,
and 9 dB and below never lock.  PothosSDR's SoapySDR binding is Python-3.9-only, so 3.9
was installed for the capture tool; radioconda's 3.12 (which has SoapySDR) drives the live
chain (`tools/atsc3_watch_av.py`), where the one-viewer lock's `fcntl` import killed the
viewer on Windows -- fixed with `msvcrt` (E117 commit).  WJLA 7.1 with 5.1 sound and
captions ran 7 minutes from this box at 8 x 2, 0 drops, 0 fails.

**The sweep (lossless replay, `--ldm-proc-threads 2 --threads 4`, as-recorded RF25):**

| P | x rt | fec | lanes vs P=4 |
|---|---|---|---|
| 4 | 0.807 | 4943/6161 | reference |
| 6 | 0.978 | 4943/6161 | 4/4 byte-identical |
| 8 | 1.042 | 4908/6126 | 4/4 byte-identical |
| 12 | 1.121 | 4698/5951 | video a shorter exact prefix |
| 16 | 1.210 | 3824/5040 | all four differ from 0x15F93C |
| 24 | 1.205 | 4970/6300 | all four differ |
| 32 | 1.209 | 4138/5426 | differ |

Real time from 8 processes; a plateau at 1.21x from 16.  At P=24 `pool_wait` is 0.0 and
`fec_dispatch` 9, and what is left is the serial main thread -- `cti_gather` 93 + front end
79 + `fec` 7 + `l1` 3 of the 242.5 ms budget -- which CONFIRMS the old m11_watch note
("past 6 the parent's own serial work is the limit"), with one addition: both of those grow
with P (cti_gather 36 -> 93, front end 28 -> 79 ms), so part of the cap is contention from
48 worker threads, not fixed work (24 x 1 thread: cti_gather 60, front end 48, the same
1.19x because the pool then waits).  The non-LDM path on a fresh 20 s RF33 bank: decode-procs
4/8/12/16 = 0.94/1.15/1.24/1.31x, PLP 0 5920/5920, lanes byte-identical across D; gate_e107
PASS on that bank.

**E115 does not hold here, and the reason is a race.**  Every run loses the CTI lock on the
same air ("FEC collapsed over the last 24 batches"); the recovery is the blind sweep, and
`_begin(self.cti.next_frame, blind=True)` anchors it at the CTI's position WHEN THE
ASYNCHRONOUS FEC VERDICT ARRIVES.  More processes in flight move that position (k=0
start_row 298/298/--/344/556/722 for P=4/6/8/12/16/24), a different anchor gives a different
dead-reckoned candidate set, and the search locks at a different start_row and Frame.  A P=4
run with only the LDPC dump switched on (file writes inside the FEC) produced different
lanes from two clean P=4 runs: same layout, different timing, different bytes.  Lanes are
byte-identical up to the lock loss in every run, so demod and FEC are layout-independent;
E115 fixed the cold path and the lost-lock path still depends on pipeline depth AND on
scheduling.  The fix has to anchor on stream order (the batch that tripped the loss window),
not on `cti.next_frame`; it is not in this commit because every lock-path gate is soft until
it lands, and the throughput work below was refereed on recorded LLRs instead.

**E117 -- the pool's own timeline.**  `LdmPipeline.close()` dropped the pool before
m11_watch called `summary()`, so every pooled JSON ever written had
`demod_ms_per_frame: {}`: the workers' timers were collected and thrown away.  Kept now,
plus queue-wait, out-queue latency and ordering-wait around each job (perf_counter stamps
cross the process boundary).  Lanes byte-identical at P=4 and P=8 to the uninstrumented
runs.  What it says, per Frame, summed over workers, P=4:

| regime | demod in worker | FEC batch (LDPC) | demod job queued behind FEC | ordering wait | copies |
|---|---|---|---|---|---|
| as recorded (80%) | 212 | 877 (706) | 856 | 151 | 4 + 7 + 13 |
| 15 dB (59%) | 243 | 974 (827) | 963 | 207 | ~same |
| 12 dB (21%) | 243 | 895 (837) | 879 | 404 | ~same |

The copies are nothing.  The pool is four-fifths FEC in EVERY regime on this capture, the
demod (132 ms standalone, gate_e111 PASS here) runs at 212-243 inside a worker because it
shares the process with the FEC, and in the marginal regime the main thread adds its own:
the synchronous probing FEC (137 ms a Frame -- 35 Blocks decoded to judge a candidate that
LOCK_BLOCKS = 4 would judge) and the L1 re-reads (77 ms a Frame) of the E115 treadmill, 30
re-acquisitions in 45 s.  Where nothing locks at all (9 dB and below) the treadmill alone
holds the rate at 0.36-0.49x: 49-81 candidates, `fec` 212-277 + `l1` 169-573 ms a Frame.

**E118 -- lane compaction, measured and parked.**  The per-block record (ATSC3_LDPC_DUMP,
378 batches, 11,633 blocks) showed each 8-lane vector running to its SLOWEST lane: 37,323
vector-sweeps where the lanes' own counts sum to 27,674 vectors' worth, the slow lane a
near miss running all 50 sweeps while seven retired lanes ride along.  Built as a second
entry point `ldpc_layered_pack_f32` (ABI 4; any nblk, lanes refilled as blocks retire, the
sweep-0 check at load so a clean block still retires with its = 0).  `gate_e118_ldpc_pack.py`:
0 mismatches on bits, conv, its and bad over 1,162 banked blocks under both share layouts,
deterministic -- and 0.98x, because with ~17 blocks per thread-share the near-miss lanes end
up ALONE at the tail of each share (26% of the pack's sweeps run with one or two lanes
active) and a refill costs about a sweep.  The recorded schedule says the gain needs the
tail amortised ACROSS batches (one 8-lane share per batch, two batches per worker at once:
24% off; 70-block merged batches: 20%), a pool change.  OFF by default (`ATSC3_LDPC_PACK=1`),
kept for that change; gate_e107 PASS with identical verdicts under it.

**E119 -- the L1 decoders on the kernel.**  Profiled: an L1 read from a held window was 316
ms, 73 of it L1-Basic's numpy min-sum of which 42 ms was `pack()` rebuilding the
code-constant check tables on EVERY call (E113's lesson, again).  L1-Basic and L1-Detail now
run their one block through the compiled flooding kernel (float32, alpha 0.75, the same
rule) and `pack()` is cached.  `gate_e119_l1_kernel.py`, 36 windows off the bootstrap grid
of the clean, 15 dB and 12 dB banks: L1 fields identical on every window where either path
verifies, the same 31 verified, **316 -> 20 ms an L1 read (8.8x)**.  On by default.  Replay
at 12 dB: 0.473 -> 0.561x (+19%), 958/4541 both; clean air unchanged and byte-identical.

**E120 -- the probe size, a knob.**  `ATSC3_LDM_PROBE_BLOCKS=k` judges a provisional phase
on the first k Blocks instead of every Block available.  k=8 at 12 dB: 0.561 -> 0.639x
(+14%), but 932 vs 958 blocks, and -4% at 15 dB, -0.7% clean -- the verdict's statistics
change.  Off by default; Felbs' call.

**E121 -- the alpha ladder is per code, and this code had it backwards.**  The record,
grouped into ladders (alpha 1.0 -> 0.85 -> 0.75, each rung only on the rows still failing):
as recorded, rung 1 converged 2,449 blocks and rung 2 converged 2,556 MORE -- blocks rung 1
had ABORTED AS HOPELESS at sweep 8 (their rung-1 syndrome 9,331-10,502, the E112 threshold
being 9,331, sized at alpha 0.85) -- and rung 3 converged 1 for 30% of all sweeps; at 12 dB
rung 1 converged 0 of 4,716 and rung 2 954.  M10 chose 1.0 first on the 64800 rate-6/15
code; RF25's core is 9/15 and the E110 correction changes the scale.  `gate_e121_ladder.py`
runs each variant on the recorded first-rung LLRs with the ladder's own semantics (1,136
blocks, two noise levels), judged on lost / gained / sweeps / bits on shared blocks:

| ladder | as recorded: conv lost gain, vec-sweeps | 12 dB: conv lost gain, vec-sweeps |
|---|---|---|
| ship (1.0, 0.85, 0.75) | 623 0 0, 4653 | 142 0 0, 4816 |
| (0.85, 1.0) | 623 0 0, 2884 (-38%) | 142 0 0, 2798 (-42%) |
| (0.82, 1.0) | 624 2 3, 3069 | 208 0 66, 2811 |
| **(0.82, 0.85, 1.0)** | **626 0 3, 3982 (-14%)** | **208 0 66, 4150 (-14%)** |
| 1.0 only | 239 384 0 | 0 142 0 |

Bits identical on every shared block for every variant.  Shipped: **(0.82, 0.85, 1.0) for
the 64800 rate-9/15 code with the correction on** -- no block lost at either level, +46% of
blocks at 12 dB, 14% fewer sweeps; every other code keeps M10's ladder; `ATSC3_LDPC_LADDER`
overrides.  (0.85, 1.0) is the pure-speed alternative if margin is not wanted.  Replay:
P=4, E119 already on: 12 dB 0.561 -> 0.598x and **958 -> 1,389 blocks (+45%)**, re-acquisitions 30 -> 19; 15 dB 0.642 -> 0.728x, 3,292 -> 3,447; as recorded 0.817 -> 0.913x, 4,943 -> 5,023, LDPC 702 -> 577 ms a batch.  gate_e107 PASS after it (the RF33 codes keep their ladder).  Day's ledger at 12 dB on this box, P=4: 0.473 -> 0.598x with +45% of blocks (0.704x with the probe knob); as recorded 0.814 -> 0.913x.

**Abort rule re-sized here** (`lab/size_ldpc_abort.py`, truth = no abort, per-file alpha):
on the as-recorded first rung the shipping (8, 0.36, 0) loses 0 and saves 75% of vector
sweeps; (6, 0.33, 0) saves 80% at 0 lost; nothing to change.

**Where this leaves the 6-core box.**  Its regime is the 12 dB row.  Levers in this commit
that reach it: L1 on the kernel (+19% here), the ladder (+46% blocks at the cliff, -14%
LDPC); knobs it can try: probe 8 (+14%, costs blocks), pack (nothing until the pool packs
across batches).  What this box could not test: its contention (the worker demod at 212-243
vs 132 standalone is the FEC sharing the process), which argues for demod-only and FEC-only
workers, and the lost-lock anchor.  The next run there is `--ldm-procs 2 --ldm-proc-threads
4` on the marginal capture with this commit, and the JSON's `pool_ms_per_frame` says the
rest.

**Live RF25 from this box (16:35, Antenna B, 8 x 2, 90 s, `--player none`).**  The first
LDM television off this machine: 358 Frames, 8,139/11,360 blocks, 0.93x sustained (live is
capped at 1.0x by the air; the deficit is the two acquisitions), 13 Frames shed, 0 SDR
overflows, all four lanes written (video 15.3 MB, 64 s), locked twice and lost twice, E115
re-acquired 11 times.  The pool kept up: `pool_wait` 5 and `fec_dispatch` 24 ms a Frame
against the 6-core box's 661 + 384 in 13.18.  **The marginal regime with E119 + E121**
(12 dB file): P=4 0.598x, 8 0.726, 16 0.788, 24 0.812, converged blocks identical (1,389)
on all four.  The cap there is the serial main thread again, 304 ms at P=24 -- front end
134 and `cti_gather` 93 inflated by contention, the probing `fec` 75 from 18
re-acquisitions -- so the marginal LDM channel is 0.81x on 32 cores and the next lever is
structural: the probe (E120 showed +14%), the treadmill, and the worker split.

**E122 -- the LDM round trip DESIGN_NOTE asked for, built and PASSED.**  `lab/gate_e122_ldm_roundtrip.py`:
Core and Enhanced FEC Blocks synthesized through the receiver's own transmit path (`PlpChain.encode`),
superposed per A/322 7.2 at the signalled injection (unit-power layers, g = 10^(-inj/20)), the Enhanced
stream at its own block offset, seeded noise at stated Core SNRs, and leg B through the 1024-row CTI
(`m10_cti.interleave`, de-interleaved by the receiver's function, cancelled in that domain -- E69's
plan).  Receiver = the live chain's demapper (right sigma^2), the layered kernel with the per-code
ladder and the E112 abort, BCH.  Canceller = E69's (re-encode, least-squares gain, subtract), run as
ORACLE (transmitted payload) and REAL (the decoded payload).  Both configurations PASS: noiseless, every
Core and Enhanced Block decodes after REAL cancellation with payloads equal, |gain| 0.9995-1.0005,
residual power g^2/(1+g^2) as predicted; a wrong payload leaves 0.84-0.94 of the received power (the
control fails to cancel, as required); through the CTI the de-interleaved cells are exact and the
Enhanced layer decodes 35/35 and 47/47.  Thresholds, synthetic AWGN: **RF30's Enhanced layer (64QAM-NUC
6/15, 4 dB injection) decodes down to 13 dB of Core SNR = 9 dB on its own layer -- E69's measured
number -- and is 0/35 at 12**; **RF25's (256QAM-NUC 7/15, 10 dB injection) needs ~22.5 dB of Core SNR
(47/47 at 23, 14/47 at 22, 0 at 21) = ~12.5 dB on its own layer.**  So the canceller inverts itself; the
August RF30 failure was the link (E69 measured ~6 dB on the Enhanced layer against 9 needed), not a rung,
and "3 dB short" is now a result.  RF25's Enhanced layer is out of reach of the yagi by ~15 dB (its Core
converges 80% at a QPSK 9/15 threshold of a few dB); the Philips, where the gain table calls RF25 strong,
is the antenna to try.  Still untested on air: the Enhanced PLP's own C anchor and the NUC alphabet at a
live SNR, which the round trip cannot see (M6's law).

**E123 -- the re-lock is a function of the stream.**  Three changes, each found by
`lab/gate_e123_relock.py` failing on the previous one: (1) the blind sweep after a LOST lock
anchors at the Frame of the batch that tripped the loss window + LOSS_LAG + 1 (default 12
Frames), with the parent settling every FEC verdict older than LOSS_LAG Frames before cutting
more Blocks so the anchor is always ahead of the CTI (`loss_lag_wait` 0 at 4x2, 77-129 ms a
Frame at 16-24x2, `loss_anchor_behind` 0 everywhere); (2) a COLD floor -- no fresh L1 from a
Frame below the point the stream went cold, because a deep pool still delivers those Frames
and a shallow one has consumed them (the 15 dB file re-locked at Frame 35 on 24x2 and 38 on
4x2 for exactly that reason); (3) ONE cold L1 read per Frame through a single routine for held
and arriving Frames, from the floor up (E115 read every held window at once, then the arriving
Frames again).  Gate, the same capture replayed as a user would: clean RF25 at 4x2 twice, 8x2,
16x2, 24x2; the 15 dB and 12 dB copies at 4x2 / 16x2 / 24x2 -- **converged counts identical,
lanes byte-identical, final lock identical on all three captures.**  Verdicts landing after
the loss tripped are no longer counted (a 35-Block fec_total wobble).  Throughput unchanged
(clean 24x2 1.22x).  The lock-path numbers of 13.16-13.20 can be re-measured with a referee
now; the E117 sweep's "lanes differ from P=12" row is this commit's before.

### 13.21 The 6-core box, with 13.20's levers (2026-11-07 evening, Omarchy, Ryzen 1600X 6c/12t, gcc 16, Python 3.14.8)

The handoff's question: what does this box do now, and where does its time go.  Same box
and install as 13.12-13.19 (the "previous OS install" clause did not apply): the existing
clone pulled to 9ad6d74, the existing venv (numpy 2.5.4, scipy 1.18.1), both captures on
disk.  All seven kernels REBUILT from the pulled source -- `ldpc_layered_kernel.c` had
changed (E118, ABI 4) and the morning's ABI-3 binary would have been refused and silently
replaced by flooding -- every smoke test PASS, the layered library loading at ABI 4.  The
radio is on the other box; replay only.  `OPENBLAS_NUM_THREADS=1`, `ATSC3_LDPC_BOXPLUS=1`.

**A. The number of record, `--ldm-procs 2 --ldm-proc-threads 4 --threads 4`, lossless.**

| | 13.17 (E111-E115) | now (+ E117-E124) | |
|---|---|---|---|
| x real time | 0.420 | **0.451** | +7% |
| converged / total | 1571 / 5763 | **1916 / 5719** | **+22% blocks** |
| frames, dropped | 184, 0 | 184, 0 | |
| lock lost, re-acquired | -- | 1, 8 | |

`ldm.stats`: l1_verified 10, lock_lost 1, reacquired_from_window 8, phase_rejected 8,
phase_rejected_range 20, pre_floor_skipped 3 (E123's cold floor), p4_agree 1;
`loss_anchor_behind` and `loss_lag_stalled` are written only when non-zero -- both 0;
`fec_after_loss` likewise (0 here; 3 at 5x1 in E).  The one loss: 19:31:44 "CTI lock LOST (FEC
collapsed over the last 24 batches) -- sweeping +-6 Frames from Frame 150", re-locked by the
blind sweep at origin Frame 164, start_row 198, 19 candidates rejected on the way.  The lanes
are, as 13.20 said they would be, not the old record's bytes: E121 changes which blocks
converge and E123 the re-lock.

**Where the time goes, ms per Frame (budget 242.5).**  Main thread: front end 30.9; decode
853.4 = pool_wait 450.0 + fec_dispatch 269.2 + fec 48.2 (fec_ldpc 40.0) + cti_gather 27.5 +
win_hold_copy 8.9 + l1 **1.3** (E119: 8.1 in 13.17, 17.9 live in 13.18) + the rest under 4;
HEADROOM 0.28x.  **The pool (E117, summed over the two workers, per Frame):** worker_qwait
1393.1, worker_demod 202.1, out_q_latency 45.7, park_wait 24.3, submit_copy 6.5;
fec_worker_qwait 1093.4, **worker_fec 780.4** (w_fec_ldpc **657.5**, w_fec_demap 81.3,
w_fec_descramble 30.5, w_fec_bch 7.9, w_fec_bitdeint 6.7, w_fec_widen 6.4), fec_park_wait
14.3, fec_out_q_latency 0.8, fec_submit_copy 2.4.  Demod stages inside the worker: cpe 92.5,
fft 50.9, equalise 48.6, ce_smooth 22.5, dummy_ref 1.2 (= the 202).

What that says, in order.  (1) **The pool is 79% FEC** (780 of 982) -- 13.20's "four-fifths
FEC in every regime", reproduced on the box it was predicted for.  (2) **The demod runs at
202 ms inside a worker against 103 standalone** (gate_e111 here, 13.15): the FEC sharing the
process, exactly as E117 measured on the 32-core (212-243 vs 132).  The copies are nothing
(6.5 + 2.4).  (3) **The LDPC is 657 ms a Frame summed over workers at 27% convergence,** and
that is the price of E121's blocks: with the ladder, a hopeless block is aborted at sweep 8 on
EACH of three rungs -- ~24 sweeps where E112 alone spent 8 -- while 13.20's "-14% sweeps" was
measured on the blocks a rung keeps.  Here 73% of blocks are hopeless.  The lever this names
is a two-level abort: a block flat at ~0.39 x nc is hopeless for every rung and should be
dropped once; the near misses just above the E112 threshold (9,331-10,502, the ones rung 2
rescued) should get the next rung.  The E124 bank written in C below is the data to size it.
(4) The main thread still spends 84% waiting; the L1 reads that were a treadmill in 13.17
are 1.3 ms.

**B. `gate_e123_relock.py --layouts 2x4,2x4,4x1,5x1`: PASS.**

| layout | x rt | converged / total | lock | lost | re-acq | anchor behind | lag wait |
|---|---|---|---|---|---|---|---|
| 2x4 | 0.452 | 1916 / 5719 | origin 164, start_row 198 | 1 | 8 | 0 | 0.0 |
| 2x4 | 0.452 | 1916 / 5719 | origin 164, start_row 198 | 1 | 8 | 0 | 0.0 |
| 4x1 | 0.455 | 1916 / 5719 | origin 164, start_row 198 | 1 | 8 | 0 | 0.0 |
| 5x1 | 0.469 | 1916 / 5719 | origin 164, start_row 198 | 1 | 9 | 0 | 0.0 |

Counts identical, lanes byte-identical, final lock identical, anchor never behind and never
stalled: the re-lock is a function of the stream on this box too, and 13.16's layout sweep --
unmeasurable then because the lock path depended on timing -- now has a referee behind it:
the layout still barely matters here (5x1 +4%).

**C. The marginal carrier's real SNR, E69's way (`ATSC3_FEC_CELL_DUMP`, every 10th batch; `e124_ldm_snr.py --inj 10 --need 12.5`).**
The dump replay reproduced A exactly (0.449x, 1916/5719; 18 batches banked, 141 MB, parent
and workers).  The analysis first CRASHED -- and the crash is a finding about the bank, on
both boxes: the writer stored `bch_ok=np.asarray(syn == 0)`, but `bch_syndrome_batch` returns
a Python LIST of ints, so `syn == 0` is the scalar `False` and every bank ever written carries
a 0-d `False` for `bch_ok` -- no per-Block BCH flag at all; `ok[i]` fails on the first file.
Fixed in both places: the writer now banks `[s == 0 for s in syn]`; the reader, meeting a
batch-scalar flag, judges a Block clean when it converged and carries a packet (exact in this
chain: `bch_zero == fec_converged`, 1916 = 1916).  On this bank, then: 562 Blocks, 219
BCH-clean and re-encoded through `PlpChain.encode`, fitted |gain| median 0.930, P_core -2.10 dB
of the total, P_residual -4.16 dB (cancellation gain 4.16 dB), negative control leaves 0.888 of
the received power (fails to cancel, as required).  **Core-layer SNR 2.05 dB** -- a QPSK 9/15
Core at its cliff, which is what 27-33% convergence looks like.  **Enhanced-layer SNR -7.19 dB
against the 12.5 dB E122 measured it needs: margin -19.7 dB, SHORT.**  13.20 estimated "out of
reach of the yagi by ~15 dB"; measured, it is 20.  RF25's Enhanced layer is not a decoder
problem on this antenna; it is a link problem, and the Philips -- where the gain table calls
RF25 strong -- is the test, on the box that has the radio.

**D. `gate_e107_ldpc.py` on the RF33 bank: PASS.**  Flooding 22,890 ms, layered 12,428 (1.8x),
layered+abort 6,600 (1.88x over layered); no net loss at any level, every shared Block
byte-identical, deterministic; the abort loses no Block at any level.  With the ABI-4 kernel.

**E. Layouts.**

| layout | x rt | converged / total | lanes vs 2x4 | worker_fec (sum) | w_fec_ldpc (sum) | park_wait (sum) |
|---|---|---|---|---|---|---|
| 2x4 | 0.451 | 1916 / 5719 | reference | 780 | 658 | 24 |
| 3x2 | 0.471 | 1916 / 5719 | 4/4 byte-identical | 1171 | 992 | 147 |
| 5x1 | 0.466 | 1916 / 5719 | 4/4 byte-identical | 1959 | 1679 | 502 |

Identical Blocks and bytes at every layout (E123 holds; 13.16's sweep is finally a
measurement), and the rate does not move: 0.45-0.47x.  The summed worker timers grow with the
worker count -- LDPC 658 -> 992 -> 1679, park 24 -> 502 -- while the output and the rate stay
put: that is contention on six cores (five single-thread workers plus the parent) inflating
each worker's LDPC, 13.20's finding at P=24 in miniature.  Main thread at 5x1: pool_wait 384,
fec_dispatch 228, cti_gather 33.

**Where this box stands.**  0.451x and 1916 Blocks on the marginal capture, every byte
layout-independent, the Core at 2 dB and the Enhanced layer 20 dB away.  The pool is 79% FEC;
the LDPC's 650-plus ms a Frame at 27% convergence is the ladder paying three aborts per
hopeless Block, and the two-level abort named in A is the next lever that reaches this box.
The cell bank for sizing it is on disk (`~/captures/cells_rf25`, 141 MB), now with a real
per-Block BCH flag for the next one written.

### 13.22 The 32-core box, second session: the Enhanced layer closed as a link, E125 sized and built (2026-11-07 night)

Pulled to 096d030 (13.21's two Python fixes; no kernel change, the ABI-4 layered library
loads with the pack entry point).  `OPENBLAS_NUM_THREADS=1`, `ATSC3_LDPC_BOXPLUS=1`.

**Part 1 -- the Enhanced layer, E69's way, on the strong RF25 (`rf25_win_45s_b.cs16`, 80% of
blocks as recorded).**  `ATSC3_FEC_CELL_DUMP` every 10th batch at 8 x 2 (the replay itself
1.143x, 4848/5986 -- E123's bytes), 22 batches, **736 Blocks, 679 BCH-clean and re-encoded**
(the per-Block flag is real now), fitted |gain| median 0.934, P_core -2.03 dB of the total,
P_residual -4.28 dB (cancellation gain 4.28 dB), negative control leaves 0.887 of the received
power (fails to cancel, as required).  **Core-layer SNR 2.25 dB.  Enhanced-layer SNR -6.95 dB
against 12.5: margin -19.45 dB, SHORT.**  Against 13.21's marginal capture (2.05 dB, -7.19,
-19.7): the two captures that look so different to the FEC -- 80% against 27% of blocks --
sit 0.2 dB apart, because QPSK 9/15's cliff is that steep, and both are 19-20 dB from the
Enhanced layer.  That is the effective SNR the decoder sees (the residual carries the
equaliser's own error too), and it is the number that matters.

The other antenna could not be measured: the RSPdx has been streaming at half rate since an
SDRplay API service restart earlier in the evening (45 s of samples in 90 s of wall, 0 overflows,
the capture tool's integrity law marks it VOID), and by the time of the port probe it answered
"no sdrplay device found" on two of the three ports and delivered 0 bytes on the third.  Device
and service restarts from software did not recover it; it needs a USB replug.  **So the one
number the question wanted is half-answered: on the yagi, both RF25 captures put the Enhanced
layer ~19.5 dB out of reach, and that is a link problem, not a decoder problem.  The Philips
measurement is one 45 s bank and one replay away once the radio is replugged; until it is
done, nothing gets built for the Enhanced layer.**

**Part 2 -- E125, sized first.**  13.21 asked for a two-level abort: drop the hopeless Blocks
(flat at ~0.39 x nc) after rung 1 instead of re-aborting them on rungs 2 and 3.
`lab/size_ldpc_abort2.py` re-assembles the ladders from the per-call record (rung k+1's rows
are rung k's failing rows) and tags every rung-1 abort with its syndrome and whether any later
rung converged it.  Two records, both the 12 dB file (the 6-core regime):

* Under M10's ladder (1.0 first; this morning's record, 4,716 Blocks): the populations do NOT
  separate.  Rung 1 converged nothing; the 958 Blocks rung 2 rescued had rung-1 syndromes
  9,502-10,404, inside the never-rescued range (p5 9,958, p50 10,211).  Smallest zero-loss T2
  = 0.405 x nc, saving 1%.  Alpha 1.0 is a bad judge of hopelessness on this code.
* Under E121's ladder (0.82 first; a fresh record, EVERY=8, 5,285 Blocks, 1,389 converged):
  **rung 1 converges 1,386 of the 1,389; rungs 2 and 3 cost 42% of all LDPC sweeps (102,503 +
  31,253 of 321,556) and converge 3 Blocks -- each of which rung 1 had left at a syndrome of
  2.**  The 533 rung-1 aborts (>= 0.36 nc) are never rescued, but dropping them saves only
  2.7%: the cost is not the hopeless Blocks, it is the 3,366 NEAR MISSES below the abort
  threshold (final syndrome p5/p50/p95 = 259 / 8,766 / 9,451 after 50 sweeps) that climb to
  rungs 2 and 3 and never converge there either.

So E125 is not the kernel's second abort threshold; it is a gate on the LADDER, in Python,
no ABI change: a Block climbs to the next rung only if the rung it failed left it close --
final syndrome under LADDER_NEAR x nc (`ATSC3_LDPC_LADDER_NEAR`).  At 0.02 x nc (518) the
record retries 220 of the 3,366 near misses and all 3 rescues (syndrome 2); the hopeless
abort is a subset.  **Gated under both laws.**  `gate_e121_ladder.py` (now with the `near` variants, the E121 ladder
as the reference, bits on shared Blocks compared) on three LLR records -- the 12 dB E121-ladder
record (605 Blocks), the as-recorded E121-ladder record (948) and 13.20's as-recorded record (773):
near 0.02 loses **0** Blocks on all three, bits identical on every shared Block, **35% / 22% / 23%
fewer vector-sweeps**; near 0.10 and 0.05 are the same; the hopeless-only rule (near 0.36) saves
3-4%; M10's ladder loses 64 / 15 / 3 against E121's.  `gate_e107_ldpc.py` PASS.  **Replay at
8 x 2: the 12 dB file 0.679 -> 0.950x (+40%) with the converged count 1,389/5,250 and the lanes
byte-identical; as recorded 1.143 -> 1.322x (+16%), 4,848/5,986, lanes identical**; LDPC summed over
the workers 1,545 -> 958 and 969 -> 781 ms a Frame.  Default ON at 0.02 (`ATSC3_LDPC_LADDER_NEAR=0`
restores E121 as shipped).  For the 6-core box this is the lever 13.21 asked for, in a form its
data agrees with: at 27% convergence its 658 ms of LDPC a Frame should drop by about the same
40%, before any pool restructuring.

**Part 3 -- the worker split, sized, not built.**  `LdmDemodPool` has one job queue and one
worker type; each worker's `nth` threads serve both its demod (a ThreadPoolExecutor the
FrameDemod shares) and its FEC batch.  Both boxes measured the price: the demod runs at
202-243 ms inside a worker against 103-132 standalone, and a demod job waits 850-1,400 ms
(summed) behind FEC batches.  The arithmetic for the 6-core box on the marginal capture, per
242.5 ms Frame: demod 103 ms single-process at 2 threads; FEC 780 ms a batch at 4 threads
(LDPC 658), ~470 with E125 at 0.02 (the later rungs gone, 42% of sweeps); the serial main
thread 110 (front end 31, cti_gather 28, probing fec 48, l1 1).  With 1 demod process x 2
threads and 2 FEC processes x 2 threads, the pipeline's slowest stage would be the FEC at
~235-300 ms a Frame (two batches in flight), i.e. **~0.8-1.0x where it is 0.45x today**, on
3-4 of the box's 6 cores, with the demod back to its standalone cost because nothing shares its
process.  Without E125 the same split is FEC-bound at ~390 ms = 0.6x.  The build is a second
job kind on a second queue (FEC-only workers take only "fec", demod-only take only "demod"),
`--ldm-fec-procs`, and the E117 timers already measure every piece of it; it is the next
session's work, after the E125 default is settled on both boxes.

### 13.23 The 6-core box with E125 (2026-11-07 night, Omarchy, Ryzen 1600X 6c/12t)

Pulled to 61f751b on the 13.21 clone and venv; no kernel source changed since 13.21 (the
diff is PRIOR_ART, gate_e121_ladder, m44_ldm, size_ldpc_abort2), so nothing was rebuilt --
verified once: the layered library loads at ABI 4 and all seven `.so` are present.  Radio on
the other box; lossless replay of `rf25_marginal_45s.cs16`.  `OPENBLAS_NUM_THREADS=1`,
`ATSC3_LDPC_BOXPLUS=1`.  The question: what E125 does to the marginal capture on six cores,
where 13.21 measured 0.451x / 1916 Blocks with the LDPC at 658 ms a Frame summed over the
workers, 79% of the pool.

**A and B -- E125 on and off, 2x4, same session, same box.**

| | E125 off (`LADDER_NEAR=0`) | **E125 on (default)** | |
|---|---|---|---|
| x real time | 0.451 | **0.601** | **+33%** |
| converged / total | 1916 / 5719 | **1916 / 5719** | identical |
| lanes | reference | **4 of 4 byte-identical** | |
| LDPC, summed over workers | 656 | **430** | -34% |
| worker_fec / worker_demod (summed) | 780 / 203 | 552 / 196 | pool 79% -> 74% FEC |
| main thread: pool_wait + fec_dispatch | 449 + 266 | 337 + 161 | waiting 716 -> 498 |
| HEADROOM (slowest stage) | 0.28x | 0.40x | |

B reproduces 13.21 to the digit (0.451x, 1916/5719, LDPC 656 vs 658; `ldm.stats` identical:
l1_verified 10, lock_lost 1, reacquired_from_window 8, phase_rejected 8 + 20 range,
pre_floor_skipped 3) -- E123 holds across sessions -- and A changes no byte of the output.
The rest of A's pool line: worker_qwait 1043, fec_worker_qwait 851, out_q_latency 29,
park_wait 20, submit_copy 6; demod stages in the worker cpe 89, equalise 49, fft 48,
ce_smooth 23 (= the 196); main thread fec 33 (fec_ldpc 25), cti_gather 28, l1 1.4.  What
E125 removed is exactly what 13.22 said it would: the later rungs' sweeps on near misses that
never converge.  The demod is untouched (196 vs 203) and still twice its standalone cost
inside a worker that shares the process with FEC -- the worker split's argument, intact.

**C. `gate_e123_relock.py --layouts 2x4,2x4,4x1,5x1`, E125 on: PASS.**  0.601 / 0.596 /
0.615 / 0.641x; 1916/5719 at every layout, lanes byte-identical, final lock origin 164 start_row
198 everywhere, anchor never behind, lag-wait 0.

**D. The ladder's own law on this box's Blocks.**  Record with E125 OFF (every 8th call keeps
its LLRs): 505 calls, 446 MB, re-assembled into 176 ladders (0.82, 0.85, 1.0) over 5,824
Blocks.  `size_ldpc_abort2.py` and the companion written here, `lab/size_ldpc_rungs.py`
(the same re-assembly; it prints the per-rung sweeps and the rung-1 syndrome of every rescued
Block, which the first tool reports only for rescued aborts):

* Blocks entering each rung 5,824 / 3,990 / 3,908; **converged per rung 1,834 / 82 / 0** --
  rung 3 converges nothing on this capture.  **Lane-sweeps per rung 227,115 / 132,249 /
  31,441 = 58.1% / 33.8% / 8.0%: rungs 2 and 3 are 41.9% of all LDPC sweeps**, 13.22's 42%
  to the point.
* **The 82 Blocks a later rung rescued had rung-1 final syndromes min / p50 / p95 / max =
  1 / 2 / 4 / 49** (0.0019 x nc at the top): **none at or above 518 = 0.02 x nc.**  E125 loses
  nothing here, with a tenfold margin under its gate; 13.22's records topped out at 20.  Of
  the 3,908 failures never rescued (p5/p50/p95 = 64 / 8,605 / 9,427), 369 lie below the gate
  and are still retried -- so E125 retries 451 of 3,990 rung-1 failures and skips 3,539.
* The hopeless abort 13.21 proposed: 148 rung-1 aborts (>= 0.36 nc; p5/p50/p95 9,335 /
  9,401 / 9,581), 0 ever rescued, smallest zero-loss T2 = 0.36 x nc saving 2,408 lane-sweeps
  = **0.6%**.  Dead on this box too, as 13.22 found; the cost was never the hopeless Blocks.
* `gate_e121_ladder.py --threads 4` on 25 first-rung batches, 774 Blocks (rate 9/15, corr 1,
  abort 8,0.36,0):

| variant | conv | lost | gain | vec-sweeps | vs ship | shared bits differ |
|---|---|---|---|---|---|---|
| E121 ship (0.82, 0.85, 1.0) | 314 | 0 | 0 | 8,999 | -- | 0 |
| M10 (1.0, 0.85, 0.75) | 283 | **31** | 0 | 9,433 | -5% | 0 |
| E121 + near 0.36 | 314 | 0 | 0 | 8,573 | 5% fewer | 0 |
| E121 + near 0.10 / 0.05 | 314 | 0 | 0 | 6,047 | 33% fewer | 0 |
| **E121 + near 0.02 (shipped)** | **314** | **0** | 0 | **5,947** | **34% fewer** | **0** |
| E121 + near 0.01 | 314 | 0 | 0 | 5,889 | 35% fewer | 0 |
| 0.82 only | 301 | 13 | 0 | 5,142 | 43% fewer | 0 |

Near 0.02 loses 0 against the E121 ladder with bits identical on every shared Block and 34%
fewer vector-sweeps (13.22: 35 / 22 / 23% on its three records); M10's ladder would lose 31
Blocks of 314 here, 0.82-only 13.

**E. Layouts, E125 on.**  3x2 **0.643x**, 5x1 **0.641x**, both 1916/5719 with lanes byte-identical
to 2x4.  With the LDPC smaller the layout now matters a little (+7%; pool_wait 337 -> 287 -> 250),
where 13.21 saw +4%; 3x2 is this box's layout for now.

**Where this box stands.**  0.60-0.64x on the marginal capture, +33-43% over 13.21 for a
Python gate on the ladder that changes no byte and is now gated on this box's own record.
The pool is 74% FEC at 552 ms summed; the demod is 196 in-worker against 103 standalone;
the main thread waits 498 of 603.  Every number 13.22's Part 3 arithmetic used for the worker
split is now measured here with E125 in place (FEC ~550 summed at 2 x 4 threads, demod 103
standalone, serial main thread ~95), and that split -- demod-only and FEC-only processes on
their own queues, `--ldm-fec-procs` -- is the next lever on this box, with 13.22's estimate of
0.8-1.0x to be measured, not assumed.  The LDPC record for sizing it is on disk
(`~/captures/ldpc_marginal_e121`, 446 MB).

### 13.24 The 6-core box, night: the radio was here; both antennas measured; the worker split (2026-11-07)

**Where the radio is.**  The 13.22 handoff said "the radio is on the other box".  `lsusb` and
SoapySDR on this box: RSPdx 230204CB48, plugged in HERE.  The half-rate / "no device" fault
13.22 fought was that box's USB or API state, not this radio's: both banks below streamed 45.0 s
of samples in 45.0 s of wall, 0 overflows, no clipping (the capture tool's integrity law).  And
`tools/atsc3_capture.py`'s own help says `--biast` powers "the Antenna B bias-T (yagi LNA)",
so Antenna B is the yagi with its LNA; the morning's note had the labels the other way round.
Neither label decided anything below -- the measurement did.

**Part 1 -- RF25 on both antennas, E69's way.**  Banked tonight, 45 s each, IFGR 32, rfgain_sel 2:
Antenna B with the bias-T on, max|s| 4,004 = -18.3 dBFS; Antenna A without, max|s| 17,348 =
**-5.5 dBFS** -- some 13 dB more raw RF25 at the ADC on A, even against B's LNA.  Lossless replay
at this box's 3x2 (E125 on), `ATSC3_FEC_CELL_DUMP` every 10th batch, `e124_ldm_snr --inj 10
--need 12.5` (the carrier signals `PLP1(L1 256QAM 7/15 inj10)`, so 10 is the right injection;
it also carries PLP16, QPSK 2/15, core layer):

| capture | convergence | x rt at 3x2 | Core-layer SNR | Enhanced SNR | margin vs 12.5 | control |
|---|---|---|---|---|---|---|
| 13.21 marginal (afternoon, B) | 27% (1916/5719) | 0.643 | 2.05 dB | -7.19 | -19.7 | 0.888 |
| 13.22 strong (32-core's bank) | 80% | -- | 2.25 dB | -6.95 | -19.45 | 0.887 |
| **tonight, Antenna B** | **99.5% (6371/6403)** | **1.136** | **2.72 dB** | -6.38 | -18.88 | 0.890 |
| **tonight, Antenna A** | **100% (6402/6402)** | **1.391** | **3.40 dB** | -5.52 | -18.02 | 0.875 |

Two results.  First, **the LDM carrier runs in real time on this six-core box when the air is
good**: 1.136x and 1.391x at 3x2 with every lane written and every block converged -- the
afternoon's 0.64x was the air, not the decoder, and the marginal capture stays the stress case.
Second, the Enhanced layer is 18-20 dB short on every antenna and every capture -- but the
table's spread is the finding.  The Core-layer SNR as E124 measures it (P_core over the
residual after a perfect re-encode cancellation, i.e. the SNR the Core's demapper experiences,
Enhanced interference included) rises only **2.05 -> 3.40 dB across 27% -> 100% convergence**,
and **13 dB more signal on Antenna A buys 0.7 dB of it**.  The residual is therefore not
thermal noise.  It is a **signal-proportional impairment** -- channel-estimate error, phase
noise, residual multipath, in some mix -- that caps the effective SNR near 3.4 dB however much
antenna is added; the numbers are coherent (a QPSK 9/15 Core including a -10 dB layer as
interference does decode at 2.7-3.4 dB, and fails at 2.05), and the control fails to cancel every
time, so the method is honest about what it sees.  13.22 called the Enhanced layer "a link
problem, not a decoder problem"; the four-capture series says the opposite half: **more antenna
cannot reach it, because the floor scales with the signal; a better estimator might.**  E122's
round trip cancelled perfectly with the true channel; on air the channel estimate's error is in
the residual.  The next step for the Enhanced layer is therefore receiver-side and testable: the
residual against channel-estimate quality (E104's comb estimate on and off, E60's smoothing), and
the round trip re-run with a modelled estimation error, to see how much of the 18 dB is ours.
Until then nothing gets built for the Enhanced layer -- but the reason has changed.

**Part 2 -- E126, the worker split, built and refereed.**  `LdmDemodPool` gains FEC-only
processes on their own queue (`--ldm-fec-procs N`, `--ldm-fec-proc-threads`); the demod workers
then take only demod jobs.  With N = 0 the FEC queue IS the demod queue and every worker has
role "both" -- today's pool, byte for byte.  `gate_e123_relock` accepts "PxT+FxG".  Committed
(efc1a07) before the first run.

*Identity (the marginal capture, `--layouts 3x2,1x2+2x2,2x2+2x2,1x2+3x1`): PASS* -- 1916/5719 at
every layout, lanes byte-identical, final lock origin 164 / start_row 198 everywhere, anchor
never behind.  `gate_e107` PASS (the RF33 path never touches the pool).  Strong-bank lanes
identical too.  The split changes no byte.

*Throughput, marginal capture (E117 timers summed over workers, ms per Frame):*

| layout | x rt | worker_demod | worker_fec (ldpc) | demod qwait | fec qwait | park_wait | HEADROOM |
|---|---|---|---|---|---|---|---|
| **3x2** (single kind) | **0.634** | 220 | 809 (630) | 1027 | 766 | 97 | 0.42x |
| 1x2+2x2 | 0.600 | 210 | 708 (539) | 449 | 1023 | 285 | 0.35x |
| 2x2+2x2 | 0.565 | 281 | 747 (582) | 272 | 1463 | 686 | 0.33x |
| 1x2+3x1 | 0.555 | 189 | 1149 (875) | 458 | 1120 | 594 | 0.32x |

*Strong Antenna-A bank:* 3x2 **1.398x** (worker_demod 176, worker_fec 319, fec qwait 449) against
1x2+2x2 1.311x (176, 317, fec qwait 60, demod qwait 563 -> 661): with one demod process the demod
becomes the cap (176 ms a Frame in one process is 1.37x at most).

**The split loses on this box -- 5% marginal, 6% strong -- and the timers say why, which is the
result.**  13.22 Part 3's premise was that the demod runs at ~200 ms inside a worker because it
shares the process with FEC, and would return to its ~103 ms standalone cost on its own.  It did
not: **210 ms in a process with no FEC in it** (220 with; 176-189 elsewhere), against 103
standalone.  The inflation is box-level contention -- six or eight worker threads plus the
parent and the front end on six cores -- not intra-process sharing; the hypothesis is
falsified.  And on the marginal capture the split takes threads away from the stage that is
actually short of them: FEC is the pool's bulk (809 of 1,029 summed), the FEC queue wait rises
766 -> 1023, and the single demod worker's output parks (97 -> 285) waiting for a parent that is
waiting for FEC.  **Decision: `--ldm-fec-procs 0` stays the default and 3x2 stays this box's
layout.**  The split is kept, off and bit-identical, for the 32-core box, where cores are not
scarce and the FEC dispatch wait -- not contention -- was the pool's cost; that is its test, not
this box's.

**Where this box stands tonight.**  A strong LDM carrier decodes in real time here -- 1.14x and
1.40x at 3x2 on tonight's two banks, every block, every lane -- so the live-TV goal is met for
RF25 when the air is good, on six cores, with every lever gated and every byte accounted for.
The marginal capture is at 0.634x (0.309x this morning), its remaining cost the FEC of a Core at
its cliff, and the next lever there is FEC work, not layout.  The Enhanced layer is capped by a
signal-proportional impairment that a better estimator, not a better antenna, would have to
move (Part 1).  Banks on disk for all of it: `rf25_antA_45s.cs16`, `rf25_antB_45s.cs16`,
`cells_rf25_antA/B`, `ldpc_marginal_e121`.

**Live, this box, 21:31 (Antenna A, no bias-T, rfgain_sel 2 / IFGR 32, 3x2, 90 s, `--player none`).**
One sample-continuity break in the first 0.3 s (the reader racing the pool's 8 s warm-up; it
re-acquired the bootstrap unaided), lock at 5.3 s, one E115 re-acquisition at 10 s, re-locked
with 35/35 Blocks BCH-clean -- and then **1.00x real time sustained to the end: 371 Frames,
11,653/12,492 Blocks (93%, the deficit being the ten seconds before lock), 0 Frames shed, 0 SDR
overflow reads, all four lanes written (video 29.5 MB, audio 1.05 + 0.59 MB, captions 126 KB)**,
instantaneous rate 1.01-1.02x, locked on 17 of 18 status lines.  The main thread's decode is
37 ms a Frame; the pool carries the rest; live is capped at 1.0x by the air.  **The six-core
box does live LDM television on a strong night**, at the edge of real time with nothing shed;
13.18's live run on the same carrier this morning was 0.32x with 246 Frames shed, on marginal
air and before E111-E125.  The 32-core box does the same with headroom (13.20: WJLA seven
minutes, 0 drops).  So the radio question is about margin, not capability: for watching, the
box with headroom; for this box, whose remaining work -- the marginal capture's FEC and the
Enhanced layer's estimator -- is replay on banks already here, the radio is not needed.

### 13.25 The 32-core box, third session: the radio back here, E126 refereed where it was aimed, the Enhanced layer's "impairment" named (2026-11-08)

Pulled to 477e791 (13.23, 13.24, E126); no kernel source changed since the 13.20 build, so nothing
was rebuilt -- verified once: the layered library loads at ABI 4 with the pack entry point and all
seven DLLs load.  `OPENBLAS_NUM_THREADS=1`, `ATSC3_LDPC_BOXPLUS=1`, every run.

**Step 0 -- the radio.**  The RSPdx enumerated (`SDRplay Dev0 RSPdx 230204CB48`), opened, took its
settings and read them back -- and streamed NOTHING: every `readStream` returned -1 (timeout) for
the whole 110 s of the first attempt, 0 samples.  So 13.22's fault was never the unit (13.24 proved
that) and not the USB port either: one `Restart-Service SDRplayAPIService` (it works without
elevation; stopping the service's process directly is denied) and the next probe streamed 5.00 s in
4.99 s.  What it took: one service restart, no replug, no port change.  Then the two banks, 45 s
each at IFGR 32 / rfgain_sel 2, both **45.01 s of samples in 44.99 s of wall, 0 overflows, no
clipping, USABLE**: `rf25_s3_antB.cs16` (bias-T on, max|s| 4,093 = -18.1 dBFS) and
`rf25_s3_antA.cs16` (bias-T off, 4,626 = -17.0 dBFS).  **Here the two ports are within 1 dB of
each other at the ADC**, where the 6-core box saw 13 dB between them (13.24); so the antennas are
not wired the same way at the two desks, or not the same antennas -- kept both, decided nothing by
the label.  Replayed at 8 x 2 (E125 on): Antenna A **2.001x, 6,436/6,437**; Antenna B **2.142x,
6,438/6,438**.  Both banks decode every Block; B is the stronger by E124 (below) and is the bank
the Part B table was built on.

**Part A -- E126 on the box it was aimed at.**  `gate_e123_relock.py --layouts
8x2,4x2+4x2,6x2+6x2,2x2+8x2`, on the 12 dB file and on the as-recorded strong file (13.22's bank):

*Identity: PASS on both.*  12 dB: 1,389/5,250 at every layout, lanes byte-identical, final lock
cold / origin 80 / start_row 676 everywhere, anchor never behind, 0 stalls.  As recorded:
4,848/5,986, locked / origin 74 / start_row 344, lanes identical (video 10,102,507 bytes at all
four).  The split changes no byte here either.  (One stat does move: `reacquired_from_window` is
19 at 4x2+4x2 against 20 at the other three on the 12 dB file -- a held-window re-acquisition
that happened one Frame later as a live read; the counts and bytes are the referee and they agree.)

*Throughput, E117 timers summed over workers, ms per Frame; main thread's own waits at the right:*

| file | layout | x rt | worker_demod | worker_fec (ldpc) | demod qwait | fec qwait | park / fec park | main: pool_wait | fec_dispatch |
|---|---|---|---|---|---|---|---|---|---|
| 12 dB | **8x2** | **0.935** | 310 | 1,142 (996) | 507 | 282 | 365 / 169 | 143 | 57 |
| 12 dB | 4x2+4x2 | 0.822 | 327 | 814 (674) | 187 | 1,135 | 708 / 42 | 175 | 121 |
| 12 dB | 6x2+6x2 | 0.933 | 432 | 1,004 (847) | 91 | 572 | 584 / 104 | 0 | 3 |
| 12 dB | 2x2+8x2 | 0.946 | 294 | 1,214 (1,051) | 741 | 251 | 298 / 111 | 94 | 3 |
| strong | **8x2** | **1.314** | 297 | 946 (758) | 437 | 184 | 203 / 275 | 128 | 17 |
| strong | 4x2+4x2 | 1.160 | 294 | 722 (537) | 196 | 1,009 | 476 / 73 | 149 | 85 |
| strong | 6x2+6x2 | 1.320 | 432 | 914 (722) | 209 | 449 | 540 / 150 | 71 | 3 |
| strong | 2x2+8x2 | 1.342 | 280 | 1,003 (831) | 1,024 | 89 | 252 / 217 | 145 | 3 |

Run-to-run on the strong file, same session: 8x2 again **1.332x**, 2x2+8x2 again **1.351x**, 4,848/5,986 and the lanes identical both times
(13.22 measured 8x2 at 1.322 and 0.950 on these files).  The two questions, answered by the
timers:

* **Does `worker_demod` fall toward its standalone ~132 ms when no FEC shares the process?  No.**
  280-327 ms with 2-4 demod-only processes, **432 with 6** -- it rises with the number of
  processes in the box, exactly as on six cores (13.24: 210 with no FEC in the process).  On 32
  cores and 64 threads that is not a shortage of cores; the Frame demod (an 8K FFT over 189
  symbols, the pilot grid, the equaliser over 1.2 M cells) is memory-bound, and twelve processes
  doing it and the LDPC at once share one memory system.  The inflation is box-level contention
  on both boxes; the split's premise is falsified here too.
* **Does the demod queue wait collapse?  Only by moving.**  With 4 demod + 4 FEC processes the
  demod jobs wait 187-196 ms instead of 437-507 -- and the FEC jobs then wait 1,009-1,135 instead
  of 184-282, the parent's `fec_dispatch` rises 57 -> 121 and the output parks 365 -> 708.  With
  2 demod + 8 FEC the FEC waits vanish (89-251) and the two demod workers queue 741-1,024 ms.
  The main thread's `fec_dispatch` does collapse to 3 ms in every split with >= 6 FEC processes
  -- 13.20's "pool cost" -- and the throughput does not follow, because the parent then waits in
  `pool_wait` instead (94-145 ms): the dispatch wait was a symptom of the pool's total rate, not
  its cause.

**Decision: `--ldm-fec-procs` stays 0 and the ncpu ladder does not set it.**  The best split
(2x2+8x2) is +1.2% / +2.1% over 8x2 on the two files, inside the run-to-run spread (13.22 to
today: 0.950 -> 0.935, 1.322 -> 1.314 at 8x2); 4x2+4x2 loses 12%.  The split is kept, off and
bit-identical, as 13.24 left it; its test is done on both boxes and it does not pay on either.
What would: fewer, fatter demod processes (the inflation is per process) and the LDPC itself,
which is 75-85% of every worker-second on every layout above.

**Part B -- the Enhanced layer's "impairment", measured, and what it turned out to be.**

*The brief's table first.*  13.24 asked whether the residual after a perfect Core cancellation
moves with the channel estimator.  Two of the four named levers do not exist on this path:
`ATSC3_CE_COMB` and `ATSC3_MARGIN` are read by `m16_margin` / `m9_fast` (the RF33 chain); the LDM
`FrameDemod` has its own smoother (`_ce_smooth`, E86c) with its own knobs, `ATSC3_LDM_CE_W` (the
boxcar half-width in symbols, 12; 0 = per-symbol estimates, no smoothing), `ATSC3_LDM_CE_NORM`
(the per-symbol complex gain divided out before the boxcar), `_CE_DETECT`, `_CE_ABS`.  So the
table runs the brief's two (expecting no-ops, and proving it by bytes) and the path's own.
Antenna B's fresh bank, 8 x 2, `ATSC3_FEC_CELL_DUMP` every 10th batch, E124 on 807-842 BCH-clean
Blocks per run, everything else identical:

| setting | x rt | converged | lanes vs default | fitted gain | P_res / P_core | Core-layer SNR | negative control |
|---|---|---|---|---|---|---|---|
| **default** (CE_W 12, NORM on) | 2.142 | 6438/6438 | -- | 0.9628 | 0.4446 | **3.52 dB** | 0.875 |
| `ATSC3_LDM_CE_W=0` (no smoothing) | 2.065 | 6436/6438 | differ (2 Blocks) | 0.9603 | 0.4769 | 3.22 | 0.878 |
| `ATSC3_LDM_CE_W=4` | 2.093 | 6436/6438 | identical | 0.9620 | 0.4538 | 3.43 | 0.876 |
| `ATSC3_LDM_CE_W=24` | 2.168 | 6438/6438 | identical | 0.9630 | 0.4436 | 3.53 | 0.875 |
| `ATSC3_LDM_CE_NORM=0` | 2.146 | 6438/6438 | identical | 0.9627 | 0.4456 | 3.51 | 0.875 |
| `ATSC3_CE_COMB=1` | 2.143 | 6438/6438 | identical | 0.9627 | 0.4453 | 3.51 | 0.875 |
| `ATSC3_CE_COMB=0` | 2.142 | 6438/6438 | identical | 0.9627 | 0.4453 | 3.51 | 0.875 |
| `ATSC3_MARGIN=0` | 2.146 | 6438/6438 | identical | 0.9627 | 0.4458 | 3.51 | 0.875 |

**Flat: every lever leaves the Core-layer SNR within 3.51-3.53 dB; the only ones that move it,
move it DOWN (less smoothing, -0.1 and -0.3 dB).**  The brief's two are byte-identical no-ops on
the LDM path, as expected.  The estimator is not where the 18 dB were.

*So the structure of the residual instead* -- `lab/e124_residual.py`, on the same banks: the
residual projected on the Core symbol has Var(Im)/Var(Re) = **1.002** (1.014 on 13.22's bank;
a phase error would be >> 1, a gain error << 1); re-fitting the gain per segment of the Block
removes **3.2% at 1,024 segments = exactly the white expectation** (0.0-0.8% below that); kurtosis
**3.02 / 3.02** (Gaussian 3.0); P_res against P_core across Blocks has no positive slope.  The
residual is isotropic, white, Gaussian and signal-proportional.  That is not an estimator's
signature, nor phase noise, nor multipath.  It is the signature of an independent Gaussian-like
layer riding on the Core -- **which is what the Enhanced layer is** (256QAM-NUC, 8,100 cells a
Block, spread by the CTI).

*And the pilots' own view* -- `ATSC3_LDM_CE_DUMP` (every 20th Frame's raw scattered-pilot grid,
per-symbol gains, smoothed grid, pilot noise floor) and `lab/e127_ce_probe.py`: the carrier is
**8K, SP12_4** (grid of 553 points at dx 12, 140 of them per data symbol, unambiguous delay span
Tu/dx = 98.8 us); the impulse response from each symbol's grid puts **87-93% of the power within
1 us of the FFT window start, 5-11% at 1-5 us, < 1% before the window, 2% beyond +-25 us**;
the equaliser's linear interpolation between grid points against the DFT interpolant of the same
grid is **-29.2 / -29.6 dB** relative to the channel (A / B) -- the interpolation is not the loss
either; and the **pilot-estimate SNR is 14.0-14.7 dB (A) and 15.5-16.5 dB (B)**, second
differences and raw-minus-smoothed agreeing within 0.9 dB (the scattered-pilot boost is in that
number; the data cells see it less the boost).  A channel this benign, pilots this clean, and a
Core "SNR" of 3.5 dB on the data cells.

**The answer: the injection level.**  `L1D_plp_ldm_injection_level` is a 5-bit INDEX in 0.5 dB
steps (A/331): RF25 signals **index 10 = 5.0 dB**, and the repo always knew it -- `m14_budget.py`
("injection level 5.0 dB SIGNALLED in L1-Detail"), `e58_m15_sp.py`, `e69_ldm.py` ("index 8 ->
4.0 dB" for RF30).  The LDM plan's `describe()` printed the raw index, `inj10`, and E122 (the
round trip's RF25 config), E124 (`--inj 10`) and every Enhanced-layer margin in 13.20-13.24 took
it as 10 dB -- an Enhanced layer at a tenth of its real power.  Re-read with the signalled 5.0 dB,
nothing else changed:

| capture | P_res / P_core | Core-layer SNR (unchanged) | Enhanced share 10^(-0.5) | **noise = P_res - P_enh** | **Core vs noise** | **Enhanced own SNR** | margin vs 12.5 |
|---|---|---|---|---|---|---|---|
| 13.22 strong (80%) | 0.596 | 2.25 dB | 0.316 | 0.280 | 5.5 dB | **0.5 dB** | -12.0 |
| tonight A (100%) | 0.490 | 3.10 | 0.316 | 0.173 | 7.6 dB | **2.6 dB** | -9.9 |
| tonight B (100%) | 0.4446 | 3.52 | 0.316 | 0.129 | **8.9 dB** | **3.9 dB** | **-8.6** |

**The "signal-proportional impairment that no antenna can move" (13.24) is the Enhanced layer at
its signalled level: 0.316 of the Core's power, 71% of the residual on the best bank.**  The
Core-layer SNR is capped at the injection level, 5.0 dB, reachable only as N -> 0 -- `m14_budget`
said exactly this in words ("a hard ceiling on Core-layer SINR") -- and the four-capture series
2.05 -> 3.40 -> 3.52 dB is that ceiling being approached, not an impairment at all.  The part that
IS noise, 0.13-0.28 of P_core, does fall with antenna (13.22's bank 0.280, tonight's 0.129: 3.4 dB
of the 13 dB the ADC saw -- the rest went into the ceiling); the residual's isotropic white
Gaussian structure is the 256QAM layer plus that noise; and the pilots' 14-16 dB is consistent
with a data-cell noise SNR of 8-9 dB once the SP12_4 boost is taken off.  Everything measured
tonight agrees once the index is read as an index.

*The round trip, re-run at 5.0 dB (`gate_e122_ldm_roundtrip.py --config rf25`, PASS, CONFIGS
corrected):* the Enhanced layer decodes at Core SNR 18 dB and not at 16 (oracle and real alike),
so it needs **~17 dB of Core SNR = ~12 dB on its own layer** (the 12.5 of 13.20 was measured at
the wrong injection and happens to be the same code threshold), and the Core decodes down to 6 dB
and below.  **Tonight's best bank is 8.6 dB short of the Enhanced layer, not 18-19.7; and the
shortfall is NOISE, which antenna, LNA and site do move** -- 13.22's "link problem" was the right
call for the wrong reason, and 13.24's "a better estimator might" was the wrong call for the
right reason.

*The modelled estimate error, as the brief asked (`--ce-err`; the noiseless synthetic stream
multiplied cell by cell by (1 + e), e complex Gaussian at a stated power relative to the channel,
E124's arithmetic on the result):*

| estimate error | P_res / P_core | E124 "Core SNR" | implied Enhanced SNR | Enhanced decodes |
|---|---|---|---|---|
| none | 0.3162 | 5.00 dB | unbounded | 47/47 |
| -30 dB | 0.3114 | 5.07 | unbounded | 47/47 |
| -20 dB | 0.3230 | 4.91 | 16.7 dB | 47/47 |
| -15 dB | 0.3507 | 4.55 | 9.6 | **0/47** |
| -12 dB | 0.3928 | 4.06 | 6.2 | 0/47 |
| **-10 dB** | **0.4410** | **3.56** | 4.0 | 0/47 |
| -8 dB | 0.5177 | 2.86 | 2.0 | 0/47 |
| -6 dB | 0.6377 | 1.95 | -0.1 | 0/47 |

An estimate error of -10 dB relative to the channel would reproduce tonight's P_res/P_core
(0.4410 against the measured 0.4446) by itself -- that is the budget: the Enhanced layer tolerates
an estimator no worse than about -17 dB, and everything in the air's residual that is not the
Enhanced layer is -8.9 dB relative to the Core.  The estimator's share of that is bounded by the
table above: the levers move the residual by at most 0.3 dB, the interpolation error is -29 dB,
and the pilots see 14-16 dB.  So the estimator is a minor term in the 8.6 dB; the noise is the
major one.  **Nothing is built for the Enhanced layer, as the brief required -- but the number to
beat is now 8.6 dB of noise on the best bank, and the lever is the link again.**  Fixed:
`m3_spec.ldm_injection_db()`, the plan prints `inj 5.0 dB (idx 10)`, gate_e122's RF25 config,
`e124_ldm_snr --inj-idx` (default the signalled index).  13.20-13.24's Enhanced-layer margins
should be read with 9.5 dB added (the arithmetic of the index) and their noise terms recomputed as
above; their Core-layer SNRs stand.

**Part C -- live LDM television, 90 s, 8 x 2, Antenna A, no bias-T (03:29, `--accel cpu`).**
The chain refused to start at first: `atsc3_run` / the watch command sit inside a
"02:10-05:40 meteor window" guard (`ATSC3_METEOR_WINDOW=off` disables it on a box with no meteor
job) -- so both live runs below carry that variable.  With the viewer (`tools/atsc3_watch_av.py`):
the chain locked and ran at **0.99-1.00x real time, FEC 10,778/11,653 at 85 s, all four lanes
written (video 27.5 MB)** -- and the viewer never opened, because the orchestrator waits for
`live_video_pid12.m4s` (RF33's video PID, hard-coded at `atsc3_watch_av.py:269`) and this
multiplex's video is PID 5017; it gave up at 120 s with "no video" while the video lane sat
beside it.  Its chain also ran at `accel=gpu-full` (the `atsc3_run` default), i.e. the GPU LDPC,
not the kernel path the replays use.  Both are viewer plumbing, not decoder work, and neither
was fixed on air.  Then once with `--player none` for the record:

| | 13.20 live here (16:35, B, 8x2) | 13.24 live, 6-core (A, 3x2) | **tonight (A, 8x2)** |
|---|---|---|---|
| sustained | 0.93x | 1.00x | **0.991x (372 Frames, 90.2 s of air in 91.0 s of wall)** |
| Blocks | 8,139/11,360 | 11,653/12,492 | **11,686/12,527 (93%)** |
| Frames shed / SDR overflow reads | 13 / 0 | 0 / 0 | **0 / 0** (raw/frm/bb queues 0/0/0 on every status line) |
| lock | locked twice, lost twice, E115 x11 | lost once at start, E115 x1 | **one break at 0.3 s (the reader racing the pool's 27 s warm-up), lock lost once at 6 s, E115 re-acquired once (held Frame 39), locked for the remaining 84 s** |
| lanes | 4 (video 15.3 MB) | 4 (29.5 MB) | **4 (video 28.5 MB, audio 1.05 + 0.57 MB, captions 148 KB)** |
| main thread per Frame | -- | 37 ms | **62 ms (front end 19, decode 44: cti_gather 31); pool_wait 0.0, fec_dispatch 1.5** |
| HEADROOM | -- | ~1.0x | **5.58x slowest stage, 3.90x summed** |

The 7% Block deficit is the first ten seconds (warm-up, the break, the re-lock), the same as 13.24;
after that every status line is `phase locked` with no backlog.  So yes: **this box now holds
lock and 1.0x with margin** -- 5.6x of it by the slowest stage -- where 13.20 had 0.93x with two
lock losses and 13 shed Frames; E123 and E125 are the difference, and the air tonight is the same
Antenna A that banked at -17 dBFS.

**Where this box stands.**  The LDM carrier replays at 2.0-2.1x and runs live at 1.0x with 5.6x of
headroom; the worker split is refereed on both boxes and stays off; the Enhanced layer is 8.6 dB
short on the best bank, of NOISE, after the index that put it 18 dB short is read as an index --
the estimator moves it by 0.3 dB at most and the pilots see 14-16 dB.  Next for the Enhanced
layer: the link (an outdoor antenna, or the LNA on the port that has the signal), then the
round trip's 17 dB of Core SNR is the number to bank.  Next for the decoder: the viewer's PID
plumbing for RF25 (two hard-coded RF33 PIDs), and fewer, fatter demod processes, since the
demod's cost per process rises with their count on both boxes.

**Addendum, 04:16 -- live LDM television WITH the viewer, on this box.**  "Do it" meant the viewer
plumbing; three pieces, each refereed on banked lanes before the air saw it:

* `atsc3_watch_av.py` now discovers the lanes from the chain's `live.json` (video = kind "video",
  the MMTP service the viewer plays, E95; audio PIDs sorted, lowest = the main programme's 5_X
  element, next = the simulcast pair) instead of RF33's hard-coded `live_video_pid12.m4s` and
  PIDs 13/14, which is why it gave up beside a 27 MB video lane above.  On air: `lanes: video
  live_video_pid5017.m4s; audio pid 5018 (5_X), pid 5019 (pair)` 50 s after the chain came up.
* The main-programme audio worker then died with an IndexError in `m30.synthesise_frames`: RF25's
  5.1 is ASPX_ACPL_3, a two-channel core plus A-CPL (E77 knew), and `Ac4Stream.decode_frames`
  narrowed its channel set to L/R on seeing that while the six-channel output lists it had already
  built kept growing for L/R only, so the grouped filterbank indexed C/Ls/Rs/lfe past their frame
  counts.  It now drops the channels the stream cannot give, and `atsc3_audio` writes those four
  silent in the 6-channel wav it opened (the A-CPL upmix is not built).  Referee: RF33's banked
  pid 13 (540 frames, 0 bad) decodes to a byte-identical wav before and after; RF25's banked pid
  5018 (2,329 frames) decodes 0 bad where it crashed, L/R rms 3,031 / 2,841, the rest 0.  (The
  TOC of that lane: 2 substreams, channel_mode 5.1, parsed 2,329/2,329; the fill gate passes on
  217 of them against 540/540 on RF33's, so the ACPL-3 substream carries bytes the size sum does
  not account for -- noted, not chased; the frames decode.)
* The caption worker died five times in ten seconds and was parked: RF25's IMSC1 cues are TTML
  clock-time with the hours counted from an epoch (`498368:29:28.234`) and `secs()` knew only
  `12.5s`.  It parses clock-time and the offset suffixes now; main()'s anchor needs one clock
  across the cues, which either form gives.  Referee: RF33's banked lanes -> the same 7 cues,
  srt byte-identical; RF25's banked lanes -> 40 documents, 330 raw cues, 15 readable.

**Then 4 minutes of live RF25 with the viewer (Antenna A, `--accel cpu`, 8 x 2, `--cc`):** the
window opened at 53 s (the pool's warm-up plus the 22 s mux lag), `eng wav: yes 5.1 -> AC-3, spa
wav: yes, cc anchor: exact`; the mux telemetry logged **28 chunks, 0 dropped, 0 fails, stall 0.0 s,
player up on every chunk, lag held at 22.0 s**; `live.srt` 316 lines and the soft-subtitle track
built; no worker died until the run was stopped (the rc=143 lines are the stop).  **Live LDM
television with picture, sound and captions, off this box, through the viewer, for the first time**
-- 13.20's live RF25 here had no viewer at all.  Commits 99180e0 (lanes + ACPL-3) and e082a75
(captions).  Still RF33-shaped and not touched: the `LANG` map's language names (on RF25 "eng" is
simply the lower PID), and `atsc3_run`'s `accel=gpu-full` default (the runs above passed
`--accel cpu` through `--chain-extra`).

### 13.26 After the picture: the fix list, worked through (2026-11-08)

The viewer ran on air from 07:15 to 07:32 on Antenna A at 8 x 2 while it was watched ("it's
playing"): 895 s of mux out, 0 dropped chunks, 0 fails, 0 stalls, lag 22-24 s; the window was
closed and 13 processes became 0 in 2 s.  Then the list from 13.25's "where it stands", in the
order it was given, each item refereed on a banked capture before anything went back on air.
Code: 8a411c6, c31e379, and the commit that carries this section.

**1. Teardown.**  The window-close path was always clean (above).  13.25's two orphaned audio
workers belonged to the OTHER path: the orchestrator killed from outside (a `timeout` wrapper,
Task Manager, a closed terminal) -- `TerminateProcess` runs no `finally`, and in the seconds
before it the supervisor had respawned two workers whose lanes had stopped growing.  Two
changes in `atsc3_watch_av.py`: every child is assigned to a Windows Job Object with
`KILL_ON_JOB_CLOSE` (the kernel kills the tree when the orchestrator's handle closes, however
it closes; Linux gets a SIGTERM handler routed into the same `finally`), and the supervisor
loop stops when the chain -- the radio's owner -- is gone, respawning nothing during a
shutdown.  Refereed through the whole stack on the banked RF25 capture (`atsc3_run
--extra '--capture ...'` under the orchestrator; `atsc3_run` now treats a replay chain exiting
rc=0 as "replay finished", not a wedge to restart 30 times): lanes found, viewer up, chain at
end of file -> "the radio owner is gone; shutting down, respawning nothing", 0 processes
left.  Then `Stop-Process` on the orchestrator itself 38 s into a second replay with 7 children
alive (chain, `atsc3 watch`, two audio workers, captions, `atsc3_tv`, ffplay): **0 left after
3 s**.

**2. Language, from the bytes we already had.**  The orchestrator's `LANG` map was RF33's PIDs
(13 eng, 14 spa) and on RF25 "eng" meant "the lower pid".  The MMTP service names no language
in anything the chain parsed -- but every lane's init segment carries it in the track's
`mdhd`: **5018 eng, 5019 spa, 5020 eng** (video `und`).  `m7_play.mdhd_lang` reads it beside
`mdhd_timescale` (same trak selection); `m11_watch` stores it on the lane and `live.json` says
it; a ROUTE segment carries its Representation's `@lang` instead; the orchestrator picks the
`--lang` track by declared language and falls back to pid order only when the lanes are silent
about it.  RF33 bank replayed: **13 eng, 14 spa, 15 eng** -- the hard-coded map, now derived
(5994/5994, lane bytes unchanged).  RF25 bank through the orchestrator: `audio pid 5018 (5_X,
eng), pid 5019 (pair, spa)`.

**3. The stereo that is.**  RF25's "5.1" is ASPX_ACPL_3: a two-channel core plus A-CPL, and
the A-CPL upmix is not built (TS 103 190-1's parametric coupling; `m28` decodes the core and
nothing parses the coupling parameters -- a decoder feature, days, not an afternoon).  13.25
wrote C/Ls/Rs/LFE silent into a 6-channel wav and the viewer encoded a 5.1 AC-3 of four zeros.
Now the worker's sidecar says how many channels the stream CODES (`channels_coded`), and the
viewer muxes a 2-coded 6-channel wav as stereo mp2 -- L/R as decoded, not the ITU fold (which,
with nothing to fold, would only have dropped it 3 dB).  Worker on the banked RF33 lane:
`channels_coded 6`, 540 frames 0 bad (path unchanged); on the banked RF25 lane:
`channels_coded 2`, 1283 frames 0 bad.  Orchestrator replay: `eng wav: yes stereo mp2
(ASPX_ACPL_3 core)`, `eng_ch 2` on every chunk, 0 dropped, 0 fails.

**4. E128 -- fewer, fatter demod processes** (the lever 13.25 left open), strong bank
(`rf25_win_45s_b`), `gate_e123_relock`, same session, E117 timers in ms per Frame summed over
workers:

| layout | x rt | worker_demod | worker_fec (ldpc) | demod qwait | fec qwait | main: pool_wait |
|---|---|---|---|---|---|---|
| **8x2** | **1.681 / 1.699** | 283 / 287 | 691 (526) / 684 (527) | 353 / 335 | 117 / 118 | 97 / 90 |
| 6x4 | 1.754 | 255 | 468 (350) | 326 | 190 | 81 |
| 4x4 | 1.529 | 204 | 379 (262) | 479 | 370 | 122 |
| 3x8 | 1.470 | 187 | 270 (178) | 469 | 370 | 136 |
| 2x8 | 1.211 | 159 | 217 (131) | 518 | 452 | 173 |

Identity PASS at every layout (4850/5986, lanes byte-identical, lock identical).  The premise
holds -- `worker_demod` does fall as processes get fewer, 283 -> 159 ms -- and it does not
pay: the FEC jobs then queue behind too few processes (fec qwait 117 -> 452) and the Frame rate
follows the LDPC, not the demod.  6x4 is the one layout above 8x2, by 3-4%, inside the
run-to-run spread.  **Decision: 8x2 stays the recipe; nothing changes.**  The LDPC is the
limit on every layout, as 13.25 said.

**5. The night, not the code.**  8x2 on this file read 1.322 (13.22) and 1.314 (13.25); today
it reads **1.681, 1.699, 1.687, 1.703** across radioconda (numpy 2.2.3) and Python312 (numpy
2.4.4) -- so it is not the interpreter and not today's commits (no kernel changed; the DLLs
are dated 11-07).  Per kind: `worker_fec` 946 -> 682-691, its LDPC 758 -> 526-531, `worker_demod`
297 -> 283-293.  The LDPC was slower by a third on those two nights and nothing in this tree
explains it; what the box was doing at the time (the 06:40 sonde hunt was found holding the
radio this morning; what ran overnight is not recorded) is the open variable.  Recorded here so
that no run-to-run comparison across sessions is read as a code effect again; within a
session the spread is 1-2%.  One more thing the two interpreters show: Python312 converges
**4848**/5986 and radioconda **4850** on the same file, each identical to itself at every
layout -- a two-Block float difference between numpy builds, not a defect, but the number to
expect when a gate is re-run under the other interpreter.

**6. The 4090.**  `_default_accel()` picks `gpu-full` on this box (torch 2.6, CUDA present),
which is why the recipe carries `--accel cpu`.  Measured on the strong bank with `--plp1` and
the default pool (6 x 4): default/`gpu`/`gpu-full` all run the LDM path (locked, 185 Frames)
at **1.568 / 1.583 / 1.592x**, timers identical -- the LDM pool ignores the engine, only the
front end changes -- against `--accel cpu` at the same pool **1.705x**.  The GPU front end is
7% slower here and converges **4864**/5986 against the CPU's 4850: its float path differs
and, on this file, in the favourable direction by 14 Blocks.  The recipe keeps `--accel cpu`;
a GPU that carried the LDPC would be a different project.

**7. The Enhanced layer** needs an antenna with 8.6 dB more gain, not code (13.25).  Not
measured today: the yagi test needs hands on the mast.

**8. Soak.**  The whole stack on Antenna A at 8 x 2 with `--cc`, 08:01:41 to 09:37:57 -- **95 minutes**
(the driver's 60-minute close failed: `_tv/player.pid` holds JSON, `{"pid": N, ...}`, not a number;
the window was closed by hand at 09:37 and the stack went to 0 processes).  Chain: **23,638 Frames at
1.00x, 0 shed**, FEC **749,522 / 815,616 = 91.9%** converged over the run -- 98-99% per five-minute
span from 08:17 to 08:47, then down to 80-85% by 09:30 as the morning carrier faded; **25 CTI lock
losses**, every one in a fade (six between 08:12:11 and 08:13:35, four around 08:50, the rest after
09:07), every one re-acquired.  Viewer: **950 chunks, 5,700 s of mux out, 0 dropped, 0 fails**; eng
muxed stereo throughout; one player respawn at 09:16:48 (`player stalled 20s (q 4) -- respawning
ffplay at the live edge`, the viewer's own catch-up rule).  Audio workers: 132,979 frames, 18 bad,
11,201 pad frames (the holes the fades left, filled with silence by design), 0 respawns.  Captions:
6,600-line srt.  Memory: the stack's working set 4.4 GB at 08:06, 6.4 GB from 08:16 on and **flat to
the end** (6.3-6.7 GB over the last 75 minutes) -- no growth.

One thing the soak taught about the telemetry: from 08:16 to 08:50 every chunk's `player write`
read 4-12 s and the stall column 4-13 s with 0 dropped.  That is not the player misbehaving and
not the codec (the A/B switch `--acpl-as-51` was written for the question and not needed): the
08:12 fade left holes, the mux held, then caught up in a burst (17 chunks in the 08:13 minute, lag
64 s) that filled ffplay's input pipe, and from then on each write waited for playback -- a FULL
pipe, steady, until the next fade drained it.  The picture was ~35 s further behind air than the
mux-side `lag` column said for those 34 minutes.  The 09:16 respawn is the policy that ends such a
state when it exceeds 20 s with four chunks queued.  A "stall" with 0 dropped after a fade is the
viewer working; a stall with no fade in the chain log would be something else.

**Where it stands.**  Branch `onair-techniques-plp1-live` is 43 commits past `origin/master`;
the merge is a decision, not a fix.  Repo-root run outputs (`out*/`, `gate_*.log`, `e12x_*/`)
are ignored since 8a411c6.  Left open: the A-CPL upmix (a decoder feature), the LDPC (the limit
on every layout), and the Enhanced layer's antenna.

### 13.27 Two of the three tests: the front end swept, Antenna B soaked -- and the supervisor caught restarting a healthy chain (2026-11-08)

Felbs is away from the antennas, so the yagi stays where it points and the third test (aim it)
waits.  The other two ran on the 32-core box while he was out: E129, a front-end sweep asking
what the SDR's own knobs can still give the Enhanced layer with the antennas as they are, and a
60-minute whole-stack soak on Antenna B, the morning's Antenna A soak (13.26) as the control.
Code: the commit that carries this section (tools/atsc3_run.py, tools/atsc3_watch_av.py,
tools/atsc3_capture.py); the architecture charts were brought to 13.26 in 4aa6648 (docs
s3b draws the LDM path).

**E129 -- the front end, 13 banks of 45 s.**  Antenna B (yagi + bias-T LNA) across IFGR 20-44
and rfgain_sel 0-4 with the bias-T on and off, Antenna A as the control at rfgain_sel 0, 2, 4;
each bank replayed at 8 x 2 with the E124 cell dump (every 10th batch) and the Core cancellation
measured on it.  10:44-11:19, the carrier as it was then.

| bank | antenna | IFGR | rfgain_sel | bias-T | Blocks (replay) | fit &#124;gain&#124; | Core-vs-noise | Enhanced own SNR | margin vs 12.5 |
|---|---|---|---|---|---|---|---|---|---|
| B_i32_g2_on | B | 32 | 2 | on | 3912/6370 | 0.933 | 5.36 dB | 0.36 dB | **-12.14 dB** |
| B_i32_g2_off | B | 32 | 2 | off | 5086/6437 | 0.939 | 5.94 dB | 0.94 dB | -11.56 dB |
| B_i20_g2_on | B | 20 | 2 | on | 2597/5600 | 0.942 | 6.13 dB | 1.13 dB | -11.37 dB |
| B_i26_g2_on | B | 26 | 2 | on | 645/3595 | 0.941 | 6.14 dB | 1.15 dB | -11.35 dB |
| B_i38_g2_on | B | 38 | 2 | on | 4511/5885 | 0.947 | 6.57 dB | 1.57 dB | -10.93 dB (best) |
| B_i44_g2_on | B | 44 | 2 | on | 212/3314 | 0.927 | 4.87 dB | -0.13 dB | -12.63 dB |
| B_i32_g0_on | B | 32 | 0 | on | **0/13542** | -- | no cells | -- | -- |
| B_i32_g1_on | B | 32 | 1 | on | 5181/6403 | 0.935 | 5.69 dB | 0.69 dB | -11.81 dB |
| B_i32_g3_on | B | 32 | 3 | on | 329/3414 | 0.932 | 4.64 dB | -0.36 dB | -12.86 dB |
| B_i32_g4_on | B | 32 | 4 | on | 422/3651 | 0.925 | 4.72 dB | -0.28 dB | -12.78 dB |
| A_i32_g2_off | A | 32 | 2 | off | 6197/6437 | 0.944 | 6.20 dB | 1.20 dB | -11.30 dB |
| A_i32_g0_off | A | 32 | 0 | off | **0/13542** | -- | no cells | -- | -- |
| A_i32_g4_off | A | 32 | 4 | off | **6377/6402** | 0.937 | 5.75 dB | 0.76 dB | -11.74 dB |

Three things the table says.  (a) The margin is **flat**: -10.9 to -12.9 dB across every IFGR,
every LNA state, bias-T on or off, both antennas -- a 2 dB spread on an 11-13 dB shortfall.  The
front end is not the limit; the aperture is, as 13.25 said and E129 now measures.  (b) The
carrier was weaker today than on the 13.25 bank: Core-vs-noise 4.6-6.6 dB against 8.9 dB then,
Enhanced own SNR around 0-1.6 dB against 3.9.  The Core still decodes at 99.6 % on Antenna A at
rfgain_sel 4 and 80 % on B at rfgain_sel 1; on B the LNA state matters for the Core (rfgain_sel
3-4 = more attenuation: 10-12 % of Blocks) and not at all for the Enhanced margin.  (c) **rfgain_sel
0 on either port decodes nothing** (0 of 13,542 Blocks) and the capture tool passed both banks
as "no clipping": peak 32,753 of 32,767, -0.0 dBFS, yet under 0.001 % of samples near the rail.
That is not clipping at the ADC; it is the front end overloading before it (E65 saw the same
shape).  `tools/atsc3_capture.py` now prints a WARNING when a bank's peak sits within ~1 dB of
the rail (>= 29,000) whatever the near-rail fraction says; the verdict line is unchanged.

**The Antenna B soak, 11:27:35-12:32:54, Antenna A's morning soak as the control.**  Same
script (`tools/soak_overnight.sh`, 8 x 2, raw queue 600, captions on), closed by the script at
12:32:48 through `_tv/player.pid`: 13 processes to **0 leftovers** in 6 s.  The first attempt
died at birth -- the E129 driver had exported the PothosSDR Python 3.9 `PYTHONPATH` for the
capture tool, and under it radioconda could not import SoapySDR (`FATAL: no interpreter can
import SoapySDR`); the soak runs with `env -u PYTHONPATH`.  The rerun:

| | Antenna A, 08:00-09:38 (13.26) | Antenna B, 11:27-12:32 |
|---|---|---|
| chain restarts | **0** | **4** (11:43:44, 11:48:16, 11:55:54, 12:29:46) |
| Frames | 23,638 at 1.00x | 15,154 in five runs (3827 / 931 / 1658 / 8148 / 590) at 0.97-1.00x |
| FEC | 749,522/815,616 = 91.9 % | 380,101/496,572 = **76.5 %** |
| CTI lock losses | 25, all re-acquired | 53 |
| viewer | 950 chunks, 5,700 s out, **0 dropped**, 0 fails | 552 chunks, 3,255 s out, **21 dropped**, 0 fails |
| worst player write / stall | 16.6 s / 19.2 s (full pipe after a fade, 13.26) | 12.4 s / 12.1 s; lag peaked at **202 s** |
| audio (eng) | 132,979 frames, 18 bad, 11,201 pads | 36,159 frames, 17 bad, 16,462 pads, 6.0x |
| captions | 6,600-line srt | 136-line srt; the caption worker died once (rc=1, 11:41:35), respawned |
| memory (stack working set) | 6.3-6.7 GB, flat | 5.4-7.0 GB, flat |
| leftovers at close | 0 | 0 |

The weaker carrier is the sky's; everything else in the right-hand column follows from the four
restarts, and those were ours.  Each one cost a black screen of 59, 70, 75 and 78 s (restart to
the viewer's `lane generation N -> N+1`) -- 4.7 of the 65 minutes -- and the 21 dropped chunks
are one event: at 11:54:56 the viewer found itself 154 s behind the head after the second and
third restarts stacked, re-anchored (E114) and dropped the backlog.  Before each restart the
chain was **processing Frames at 0.97-0.99x and probing the CTI phase** -- the E115 tracker doing
what it does on a fade -- and the supervisor killed it.

**Why: the fade detector never saw an LDM status line.**  `atsc3_run.py` learnt on 8/09 (E63, 41
restarts in one night) that a fade is not a wedge: when the lanes stop growing it reads the
chain's own status lines and, if the Frames advance while the last minute's FEC is under 20 %,
waits instead of restarting.  Its parser looked for the non-LDM line's `[xx.x% now]` field.
The LDM line (`lab/m11_watch.py`) prints `FEC conv/total BCH0 [phase locked]` and no `% now`
-- so on RF25, since E82, `chain_status()` returned `(None, None)`, the fade branch was skipped,
and every plateau of the lanes longer than 3 strikes + 4 grace at 5 s = **35 s** was "WEDGED and
did not self-heal -- restarting".  Antenna A's 25 lock losses all re-acquired inside 35 s;
Antenna B's did not, four times.  And we could not read the supervisor's verdicts afterwards,
because the orchestrator spawned it (and the caption worker) with stdout on DEVNULL -- the four
"why"s and the caption worker's traceback went nowhere.

Two fixes, both refereed off air.  (1) `chain_status()` parses the LDM counters too: the
interval FEC is the difference of `conv/total` over the window, and the window is the quiet
stretch being judged, (strikes + grace) x check seconds, not a fixed 60 s.  Refereed on the B
soak's own chain.log, cut at each of the four restart instants, asking what the supervisor would
have read:

| window | restart 1 | restart 2 | restart 3 | restart 4 | verdict |
|---|---|---|---|---|---|
| before 13.27 | (None, None) | (None, None) | (None, None) | (None, None) | wedge path, all four |
| 60 s (12 lines) | 29.5 % | 20.6 % | 17.0 % | 26.2 % | one "fade -> wait" of four |
| **35 s (7 lines)** | 21.0 % | 6.7 % | 7.6 % | 14.3 % | **three of four wait** |

The locked minute before a 35 s plateau dilutes a 60 s average above the 20 % line; the 35 s
window reads the plateau itself.  Restart 1 stays on the wedge path at 21.0 % -- and its lanes
WERE growing in the last 10 s before the kill (MPU 1166 -> 1168), so a different branch of the
supervisor fired there, most likely `LANE(S) DEAD` on a lane that was quiet while a sibling
wrote.  That one is not resolved by this section; it will be read, because of fix (2): quiet
children now append to `<live-dir>/<worker>.log` (`atsc3_run.log`, `atsc3_subs.log`) instead of
DEVNULL.  The non-LDM `% now` path is unchanged (synthesized log: 93.0 %); a log with no status
lines still returns `(None, None)`.  Gate: the whole stack on a replay of A_i32_g4_off through
the orchestrator -- `atsc3_run.log` carries the supervisor's four lines ending `replay finished
(chain rc=0) -- done` and `0 restart(s)`, `atsc3_subs.log` exists, the viewer and workers came
down with the chain, rc 0, 0 leftovers.

**Where it stands.**  The Enhanced layer is antenna, now with a 13-bank table under the word.
The supervisor no longer restarts an LDM chain for re-acquiring; tonight's overnight soak
(Antenna A, 20:30-06:30, the same script) runs the fixed supervisor and keeps its log, and the
first thing to read tomorrow is `soak_night/atsc3_run.log`.  Left open: restart 1's branch, the
caption worker's rc=1 (both answered by the logs from now on), the A-CPL upmix, the LDPC, and
the yagi, which needs a person.

### 13.28 The OnAir fork begins: their November receiver built here and raced on both carriers (2026-11-08)

Felbs' order of work from here: the OnAir fork first, then OnAir's gain-servo ideas into
gr-rxtune (hold while locked, re-run on triggers, revert a step that hurt -- and NOT their
MER/SNR objective), then ATSC 1.0 with that policy.  The plan and the licence rule are in
`docs/ONAIR_FORK.md`; this section is the first measurement.

**Survey** (upstream 25034214, 2026-11-04; 3,305 commits, one author on `core/src/atsc3`):
C++20 receiver library + ImGui app, Rust protocol layers required to build, every receiver ends
in an MPEG-TS for their FFmpeg player.  ATSC 3.0 is 8.7k lines: bootstrap with CFO, L1, every
Subframe, CTI, LDPC 16200/64800 (GPU on Windows and macOS), ALP, IP, ROUTE.  **No MMTP** ("not
supported yet"), **no LDM** beyond parsing the injection field, AC-4 only through libavcodec,
tests on their own generator with impairments, and a README asking for real recordings.  Built on
this box the way their README says (MSYS2 UCRT64 + Rust): 12 minutes, 1,391 targets, no errors.

**The race, same two 20 s banks, one run after the other** (their `atsc3rx` 12 threads, ours 8
processes; `C:\captures\onair\race.txt`):

| bank | | OnAir 2026-11-04 | ours |
|---|---|---|---|
| RF33 (`rf33_win_20s`, Antenna B) | wall | **7.1 s, 2.8x** | 13.4 s, 1.47x |
| | Blocks PLP 0 / 1 / 16 | 5920 / **9360/9360** / 80 -- 15,360 baseband packets, 0 bad | 5920 / **9333/9360** / -- |
| | services | 10 found, 9 ROUTE into a 7.07 MB TS; WJLA 7.1 "MMTP (not supported)" | WJLA 7.1 playing lanes, eng + spa, captions |
| RF25 (`rf25_s3_antA`, Antenna A, LDM) | wall | 23.7 s, 0.8x ("frame cut and resampling" 15.5 s) | **9.1 s, 2.19x** |
| | Core layer | 3,130 baseband packets, 0 bad, 0 frames failed | 2,835/2,835 Blocks, phase locked on Frame 0 |
| | services | **0 service list tables, 0 bytes, exit 2** | 37 MPUs, lanes written |

Three facts.  (1) Since the 13.2 race two days ago their PLP 1 bursts are gone (9046 -> 9360)
and they are 1.9x faster than us on RF33 -- our speed claim from 13.2 no longer holds against
their current code, and we left 27 PLP 1 Blocks that they did not: an item for our list.
(2) On RF25 their plain-PLP path carries the whole Core layer out (the Enhanced layer is -5 dB
noise to it, as 13.25 predicted) and then resolves no service list, so nothing plays; the fastest
path to RF25 on OnAir is signalling and MMTP, not physics.  (3) Their time on RF25 goes to
resampling 6.912 Msps into their frame cutter; ours bypasses the resampler at the native rate
and is 2.6x faster there.  Where the fork earns its keep is unchanged: MMTP (WJLA and the RF25
Core-layer service are both MMTP here), AC-4, LDM's Enhanced layer, and recordings from real air.

**The first fork change, an hour later.**  Their 0 service list tables on RF25 were not the
physics: the 20 s hold 40 LLS datagrams (224.0.23.60:4937), the same count as RF33, but on RF25
every one is a **SignedMultiTable** (LLS_table_id 0xFE, 20 of them) or UserDefined (0xFF, 20), and
their `parseLls()` accepted tables 1-5 only.  Their generator never makes a 0xFE, so no test could
have shown it; our E69 read the same container on RF30 in August (count, then count x {id,
version, 16-bit length, gzip}, then a 16-bit signature length and a CMS signature).  Written
fresh in C++ for the fork (branch `felbs/real-air` of the upstream clone, commit 699015b6, GPL):
`parseLls()` walks the container and returns the SLT inside, `makeLlsMulti()` builds one for the
tests, the unit test covers the SLT inside, a container without an SLT, truncation (and caught
my first version, which stopped at the SLT and misread the signature length).  Their 25 ATSC 3.0
tests pass.  Rerun on the bank: **20 service list tables, 12 services** (45.1 WBFF, 54.1 WNUV,
2.1 WMAR, 11.1 WBAL, 22.1 WMPT, ...) where there were none -- and still 0 bytes of transport
stream, for a reason the datagrams state: the SLT puts every ROUTE service's signalling at
239.255.N.1:8000 (WBFF 239.255.45.1, WMAR 239.255.2.1, WBAL 239.255.11.1), and **none of those
addresses occurs among the 10,040 Core-layer datagrams**.  The ROUTE services ride the Enhanced
layer.  The Core carries one MMTP service -- 239.255.54.1:8000, 7,373 datagrams -- plus a few
small ROUTE flows (239.255.45.100:5000 TSI 10, 239.255.54.1:5100), and the SLT names that MMTP
flow: **service 5007, 54.1 WNUV**.  So the picture this receiver has shown from RF25 since
8/08 is WNUV's (the CW affiliate) on WBFF's lighthouse, not WBFF's own; where these notes say
"Fox" for the RF25 service, read "the service on WBFF's RF25 Core layer, WNUV 54.1".  For OnAir
the consequence is exact: on this carrier it can play nothing until it has either MMTP (the
Core) or the Enhanced layer (the ROUTE services), and MMTP is the one that does not need 12 dB
more antenna.

### 13.29 The speed items: where OnAir's 1.9x came from, and two fixes that take it back (2026-11-08)

Felbs: "get started with 1, 2 and 3 in that order" -- a C kernel for the subframe-1 per-cell
stages, an int8 sixteen-lane LDPC, then more processes.  The first measurement withdrew item 1
and reordered the rest; this section is the measurements, in the order they were made.

**Standalone, the PLP 1 decode is already fast.**  `lab/m17_plp1_live.py` alone on the RF33 bank,
two threads: 0.40-0.43 s per Frame -- pool 160-195 ms, CPE 69, demap 52, LDPC 62-65 (117 Blocks
of 64800, 2.2 iterations mean), BCH 13.  The same stages inside the live pool at 8 x 2 read pool
280, CPE 134, demap 66, **LDPC 450-480**.  So the live chain's cost is not Python around C; the C
itself runs 2-7x slower when eight copies run.  (Note for the record: the box is a Threadripper
**3970X**, 32 cores, uniform memory, 128 GB -- not the 2990WX earlier notes assumed.)

**Independent processes, same curve.**  N standalone decoders at once, no shared queue, no shared
memory, median per-Frame stage times:

| decoders running | pool | CPE | demap | LDPC |
|---|---|---|---|---|
| 1 | 160 | 69 | 52 | 62 |
| 4 | 197 | 90 | 57 | 124 |
| 8 | 250 | 124 | 64 | 226 |
| 16 | 402 | 226 | 69 | 356 |

Four decoders pinned one per core complex: LDPC 74 ms.  The same four packed into one complex
(16 logical CPUs): 134.  Into four cores: 204.  Pinning the live pool to each 16-CPU group in turn
gave the same 16-18 s on every group, so no die is "far" from memory: it is cache capacity, and
the float32 eight-lane kernel's working set is 1.4 MB per 64800-bit Block (`tot` n x 8 x 4 B +
`msg` 279k edges x 8 x 4 B, shared by eight) -- 8 workers x 2 groups do not fit.  Total CPU for
the 20 s bank (psutil over every descendant): 20 s cut 260 s, 10 s cut 201 s -> **1.5 s of CPU per
Frame for the decode** and ~140 s of start-up (eight workers each building 110 frequency-
interleaver sequences in pure Python, 9.5 s apiece -- a separate item, tune-in latency).  OnAir's
`atsc3rx`: 71 s of CPU for the same Frames.

**Fix 1 -- placement (committed 130ba75, 8a73e71).**  `ATSC3_WORKER_SPREAD=G` pins worker i to
logical-CPU group i mod G, in both pools (`m11_watch._decode_proc_main`, `m44_ldm._demod_worker`).
The scheduler was not spreading them.

| bank, 8 x 2 | unpinned | spread 4 | spread 8 |
|---|---|---|---|
| RF33 20 s | 13.4 s (1.47x) | 11.9 s | **10.5-10.6 s (1.87x)** |
| RF25 45 s (`rf25_win_45s_b`) | 26.7 s (1.68x) | 23.2 s | **20.1 s (2.23x)**, twice |

Blocks and the video lane's sha identical in every layout.  Placement only; `tools/soak_overnight.sh`
exports G=8 for tonight.

**Fix 2 -- narrower lanes (E130, E131).**  Item 2 as given was int8 sixteen lanes, OnAir's layout
(their `atsc3_ldpc.cpp`; idea credited in NOTICE).  Built (`ldpc_layered_q8`, ABI 5; posteriors
+-96, messages +-31, no two-minimum correction) and refereed on the E121 first-rung records with
`lab/gate_e130_q8.py`: at its best scale int8 alone converges **689 of the float32 kernel's 772**
on the strong records and 129 of 230 at 12 dB -- the 1-unit step and the missing correction cost
11-44 % at the cliff.  A hybrid (int8 first, float32 on the rest) loses nothing and keeps every
byte, but near the cliff the float32 pass dominates: RF33 chain 10.6 -> 9.0 s, RF25 12 dB bank
29.9 -> 31.9 s.  Kept in the file, off.

So **int16** (`ldpc_layered_q16`, ABI 6): sixteen lanes in one AVX2 vector, the batch's mean |LLR|
at 1024 units (a step of 0.1 % of the mean), posteriors +-24576 and messages +-8191 so their
difference never overflows, saturating arithmetic, alpha in Q15, and the E110 correction kept --
computed per check on the two minima in float32 (16 values, not per edge) and rounded back.  Half
the float32 footprint, 0.69 MB per Block.  The referee:

| records (RF25 Core 9/15, alpha 0.82) | float32 + box-plus | int16 at 1024 | lost | bytes on shared | ms per record |
|---|---|---|---|---|---|
| strong, 948 Blocks | 772 | **775** | 0 | identical | 101 -> **70** |
| 12 dB, 605 Blocks | 230 | **265** | 0 | identical | 159 -> **102** |

More Blocks, none lost, faster even single-process (less memory per sweep).  At 4096 units the
posterior clamp costs 178 Blocks; 256 is a hair worse than 1024; 1024 it is.  Whole chain, 8 x 2
spread over 8, every converged Block BCH-clean:

| bank | float32 kernel | int16 kernel |
|---|---|---|
| RF33 20 s | 10.6 s (1.87x), PLP 0 5920 | **9.7 s (2.04x)**, 5920, video bytes identical |
| RF25 45 s strong | 20.1 s (2.23x), 4850/5986 | **17.9 s (2.51x), 4887/5986** |
| RF25 45 s at 12 dB | 29.9 s (1.50x), 1390/5250 | **27.5 s (1.63x), 1551/5250** |

`ATSC3_LDPC_Q16=1` is the default from 8dc5355; `0` restores float32.  The LDM pool's Block count
is the same at 8, 12 and 16 workers (4887), i.e. deterministic; the RF33 PLP 1 count is not (next
paragraph).

**Item 3 -- more processes.**  With both fixes: RF33 8 x 2 9.7 s, 12 x 2 10.6, 16 x 2 11.2;
RF25 8 x 2 17.9 s, 12 x 2 19.4, 16 x 2 20.5.  Eight remains the number on this box; the footprint
still bounds it.  The next lever is the demod itself: the samples arrive complex64 and are widened
to complex128 at `m9_fast` line 141, and every pool stage (FFT, CFO, estimate combine, the f64
equaliser kernel) runs in 16 bytes per cell -- the pool stage is now the largest at 8 workers
(RF33 sf1_pool 265 ms vs LDPC 207).  A complex64 demod halves that footprint; it changes the
LLRs in the last bits, so it needs the converged-set / shared-bytes gate, not the sha.  Not started.

**Where we stand against OnAir on this box.**  Their `atsc3rx` on the RF33 bank: 7.1-10.8 s run
to run (GPU on or off makes no difference; the spread is theirs).  Ours: 9.7 s, stable to 0.2 s.
The gap of the morning (13.4 vs 7.1) is 9.7 vs 7.1-10.8.

**Not a speed item: the PLP 1 count moves.**  9317-9336 of 9360 across today's RF33 runs with the
same input (OnAir: 9360).  Two pieces of per-worker memory in `m17_plp1_live.py` make a Frame's
result depend on which Frames that worker saw before: the residual-CFO feed-forward `df_prev`
("each decode process tracks its own") and the comb estimate's `_tau0` seed cache.  Both break the
rule the LDM pool is built on (E85: a worker decides nothing it carries between Frames).  The fix
is a self-sufficient Frame -- estimate from its own pilots, rotate, transform again; seed the
delay search from a fixed default -- about one FFT per Frame, with the E123 gate (same bank, two
layouts, two runs, identical counts and bytes).  Only then can the 27 missing Blocks be read.

### 13.30 PLP 1 made deterministic, then whole: 9360 of 9360 (2026-11-08)

Felbs: "do the PLP 1 fix".  13.29 had found the count moving run to run (9317-9336 of 9360 on
one bank) and named the suspects: two pieces of per-worker memory in `lab/m17_plp1_live.py`.
Code: 6377d98 (E132-E134), the m44 probe knob in the commit that carries this section.

**E132 -- a worker decides nothing it carries between Frames.**  The residual-CFO feed-forward
(`df_prev`: last Frame's measured offset applied before this Frame's first transform, per decode
process) and the comb denoiser's per-group delay seed (`_tau0`) are gone.  Every Frame measures
its own residual on the first transform and, when it exceeds 0.5 Hz (`CFO_REDO_HZ`), rotates and
transforms again itself; the denoiser seeds its delay search from the Frame's own transform.  Four
runs over three layouts (8 x 2 twice, 4 x 4, 16 x 1): **one per-Frame sequence**, and the count
rose to 9344 -- the carried correction had been hurting as well as wandering.

**E133 -- the 16 that remained, read instead of guessed.**  The per-Frame log line now says why
(residual and applied CFO, LDPC iterations, pilot coherence, the syndrome weight of each failed
Block): all 16 were near misses, 2-10 unsatisfied checks at the 50-sweep cap, in Frames 15-30,
with the carrier offset corrected and the pilots coherent -- and the LDPC's mean iterations there
3.6-5.3 against 1.9-2.2 elsewhere, so the cells were worse, not the decoder.  The standalone
decoder (`m17_plp1_live.py` main, its own loader) over the same seconds: 117/117 on all sixteen
Frames at 1.9 iterations.  So the live chain's *input* to the Subframe differed.  The failing
Frames' exact windows and start samples were banked from the live workers (`ATSC3_SF1_DUMP`) and
replayed offline (`lab/e133_sf1_replay.py`, `e133_sf1_probe.py`, `e133_sf1_scan.py`):

| on the banked live windows | converged of 1287 (11 Frames) |
|---|---|
| as the live chain had them | 1271 (the live result, reproduced exactly) |
| through the standalone loader's 2.95 MHz band limit | 1275 |
| FFT windows moved **later** by +20 / +40 / +100 / +150 / +200 samples (two Frames) | 117/117 at every shift from +40; iterations 3.6-4.0 -> 1.8-2.0 at +150 |
| moved earlier by -20 / -40 / -100 / -200 | 115-116 / 114-116 / 110-112 / 96-102 |
| +300 and beyond | 0 (out of the guard) |

The cyclic-prefix fold of each Frame's own Subframe-1 samples (x[n] conj x[n+N] over the guard,
folded modulo the symbol) put its peak -- the true guard start -- **177 and 195 samples after**
the start the decoder assumed, with the region above 90 % of the peak only ~310 samples wide of a
1536-sample guard: a channel whose delay spread uses most of the guard, so the window must sit
near the true guard end, and the assumed one sat ~180 samples early at the edge of the usable
zone.  Not the LDPC, not the CFO, not the band limit: the FFT window.

**E134 -- the Subframe places its own windows.**  Per Frame, the fold (computed as folded sums:
one product pass over the segment, then a 17,920-sample circular sliding sum, ~10 ms) gives the
peak, and the FFT windows start at peak + GI - 32 (`TIMING_MARGIN`; `ATSC3_PLP1_TIMING=0` turns
it off).  A pure function of the Frame.  RF33 bank, 8 x 2 twice and 4 x 4: **9360/9360**, one
per-Frame sequence, iterations 2.06 mean in the once-failing Frames, PLP 0 bytes identical, wall
9.8 s (2.01x) -- the fold's first version cost 140 ms a Frame and 10.6 s; the folded-sum version
costs the 0.1 s.  Across the 80 Frames the shift ran -65 to +183 samples (median +45): the
assumed boundary wanders Frame to Frame, which the per-Frame fold absorbs.  OnAir's 9360 is
matched; it was never the decoder's cliff.

**E135 -- the same question asked of the LDM path, answered no.**  `ATSC3_LDM_TSHIFT` (a probe
knob on the data windows in `m44_ldm.cell_pool`) swept on the strong RF25 bank: 0 is the optimum
(4887 Blocks; -30 4852, +30 4670, +60 3721, +200 537) -- the LDM demod's pilot-coherence fine
timing already sits on its 512-sample guard's peak.  RF25's 81 % is the air (13.25: Core-layer
SNR), not timing.

**Where ATSC 3.0 stands tonight.**  RF33: both PLPs 100 % on the bank (5920 + 9360), 2.0x on
eight workers, deterministic from any layout, live television with both audio languages and
captions.  RF25: the Core layer live at 2.5x on the bank with every converged Block BCH-clean,
the Enhanced layer still 11-13 dB of antenna away.  Open and known: the demod's complex128 (speed),
each worker's ~17 s of pure-Python start-up (tune-in latency), the A-CPL upmix, the yagi.  Next,
by Felbs' order: gr-rxtune takes OnAir's three servo ideas.

### 13.31 Everything from today, live on both carriers before the sonde hunt (2026-11-08, 17:14-18:35)

Felbs: "why don't we use all of the changes we made in the soak right now?  I kinda wanna stress
test the new changes.  Make sure they work on every channel.  And doesn't break any of our
telemetry."  So, the whole stack with every change of the day -- worker placement (13.29), the
int16 LDPC (E131), the deterministic and self-timing PLP 1 (E132-E134), the supervisor's fade
detector on the LDM line and its own log (13.27) -- on RF33 and then RF25, each the full orchestrator
with viewer, both audio workers and captions, closed by the script, the radio handed back at 18:35
for the 18:40 sonde hunt.  `tools/soak_overnight.sh` took `RF` and `CHAIN_EXTRA` from the
environment for it (d6753ff).

**A false start, mine.**  The first RF33 run was pointed at Antenna A: bootstrap locked at peak
ratio 66-86, FEC 0 of 35,668, re-acquiring every eight Frames, the orchestrator gave up at 150 s
with no video lane.  RF33 has only ever been received on Antenna B with its LNA (`watch_tv.bat`);
a strong bootstrap with zero payload is what a wrong antenna or an overloaded front end looks like
(E129), and no code was involved.  Killed by PID, restarted on Antenna B.

| | RF33, Antenna B, 17:14-17:54 | RF25, Antenna A, 17:54-18:35 |
|---|---|---|
| chain | 9,642 Frames at **1.00x**, 8 x 2 spread over 8 | 9,784 Frames at **1.00x**, LDM pool 8 x 2 spread over 8 |
| PLP 0 / Core | 713,487 / 713,508 Blocks (99.997 %) | 319,534 / 333,930 (95.7 %) |
| PLP 1 (per-Frame log) | 1,123,902 / 1,129,986 (99.46 %); 116 of 9,658 Frames short, in two fades (coherence 0.992-0.995, dozens of Blocks down together); timing shift -2 to -31 samples tonight | -- |
| lock / acquisition | one bootstrap re-acquisition at start (an SDR overflow), none after | 11 CTI lock losses, every one 17:55-18:00 (a six-minute fade), all re-acquired |
| supervisor (`atsc3_run.log`) | 0 restarts, 0 strikes | **0 restarts**; `LANE(S) DEAD` struck 5 times, each cleared at the next tick -- see below |
| viewer | 392 chunks, 2,348 s out, **0 dropped, 0 fails**, no write over 1 s, lag 22-28 s | 387 chunks, 2,318 s, 4 dropped, 0 fails; one re-anchor at 18:00:41 (126 s behind after the fade) and one player respawn at 18:04:39 (the 20 s / q 4 rule), both the viewer's own recovery |
| audio | 60,349 frames, 0 bad, 191 pads, 0 respawns; 5.1 as AC-3 | 56,037 frames, 0 bad, 1,023 pads, 0 respawns; stereo mp2 (ASPX_ACPL_3) |
| captions | 3,008-line srt; `atsc3_subs.log` 350 lines, no error | 3,024-line srt; no error |
| memory | 8.6-8.8 GB, flat | 4.4-4.5 GB, flat |
| close | 0 leftovers in 8 s | 0 leftovers in 6 s |

**What the telemetry showed that the morning could not.**  With the supervisor's log kept, its
decisions during the RF25 fade are on record: at 18:00:36 it reached "candidate wedge -- allowing
in-chain self-heal (1/4)", at 18:00:41 "LANE(S) DEAD: audio_pid5018, audio_pid5019 unchanged while
2 lane(s) still write (4/3) ... (2/4)", at 18:00:46 "recovered".  So the branch behind the
morning's unexplained restart 1 (13.27) is named: in a degraded carrier the audio lanes pause
while video and captions still write, and the rule written for the 8/07 flipped-bit case (FEC at
100 %, one asset dead) read it as a dead asset.  It recovered inside the grace here; in the Antenna
B soak it did not.  Fixed in 8746d09: the branch asks the chain's interval FEC first and, under
90 % (`LANE_DEAD_FEC`), logs the quiet lane as the air and counts no strike.  Not exercised live
yet (the running supervisor predates it); tonight's soak runs it.

**Where it stands.**  Both carriers, the whole stack, every change of the day, 1.00x each, nothing
broken in the telemetry -- the telemetry got richer (the per-Frame PLP 1 line, the supervisor's and
the caption worker's logs).  The overnight soak (Antenna A, 20:30-06:30) runs the same code plus
the lane rule.
