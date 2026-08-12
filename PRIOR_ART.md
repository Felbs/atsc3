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
