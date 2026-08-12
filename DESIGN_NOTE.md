# ATSC 3.0 / NextGen TV — design note (task #37)

Status: **M1 built and gated.** M2–M7 are planned below with honest difficulty
and risk estimates. Nothing here has touched a radio.

## Why this is a separate repo (``)

Not folded into `Software-TV-Tuner`. ATSC 3.0 shares *nothing* with 8-VSB below
the transport layer — different modulation (OFDM vs single-carrier VSB),
different FEC (LDPC+BCH vs RS+trellis), different framing, different transport
(ROUTE/DASH over IP vs MPEG-TS). Dropping a multi-thousand-line PHY into STVT
would pollute four branches (`main`, `main-universal`, `main-linux`,
`pi-port-stvt`), three Windows checkouts, and the Pi/Ubuntu rigs, for code that
cannot be tested on most of them yet. M7 is the merge milestone; until then this
stays a clean-room repo with its own gate.

---

## The ladder

Difficulty is 1–5. "Risk" is the thing most likely to eat a week.

### M1 — Bootstrap detect + signaling decode ✅ DONE
**Difficulty 2.** ATSC A/321 is 24 pages, fully public, and completely fixed:
6.144 Msps, 2048-FFT, 4 symbols, Zadoff-Chu root 137, published PN seeds. It can
be synthesized bit-exactly, so the detector is gateable with zero air time.

**Biggest risk (realized, and mitigated):** the self-consistent-but-wrong trap —
a wrong ZC root or PN wiring gives a detector that passes its own selftest and
fails on real air, with no way to tell which. Mitigated three ways: (1) a second,
**assumption-free structural detector** that keys only on the A/B/C time geometry
and needs no ZC/PN at all; (2) `PN_VARIANTS` + `--pn-sweep` so a real capture can
arbitrate the one genuinely ambiguous reading (ASSUMPTION A1); (3) the A/321
Annex B Gray-code table used as a **gold test vector** from the spec itself.

### M2 — L1-Basic / L1-Detail (the preamble)
**Difficulty 4.** A/322 §6. This is where the difficulty steps up hard. The
preamble that follows the bootstrap is itself an OFDM symbol set whose FFT size,
guard interval and pilot pattern are signaled *by* `preamble_structure` (the 8
bits M1 already decodes — see A/322 Table 9.5). L1-Basic is 200 bits, LDPC-coded
at rate 3/15 with QPSK, and carries the parameters needed to find L1-Detail,
which is variable-length (up to ~4500 bits) and describes every PLP in the frame.

**Biggest risk:** L1-Basic needs a *working LDPC decoder and the exact 
A/322 parity-check matrices* before you can read a single field — there is no
partial credit and no intermediate signal to debug against. Unlike M1 you cannot
synthesize your way to confidence, because a self-built encoder tests your
encoder against your decoder, not against the spec. Mitigation: implement the
A/322 LDPC encoder from the published tables, generate a frame, and *also* verify
against a real off-air capture as early as possible — L1-Basic has enough fixed
and range-limited fields (FFT size, GI, pilot pattern, version) that a wrong
decode is usually obvious.

### M3 — Grid sync, channel estimation, equalization, MER telemetry
**Difficulty 3.** Conventional OFDM receiver work: fine timing from the guard
interval, CFO/SFO tracking, scattered-pilot interpolation, per-subcarrier
equalization. The rig already has the MER-dial apparatus from STVT
(`mer_dial_universal_algorithm`) and the equalizer research platform, and the
concepts transfer directly.

**Biggest risk:** ATSC 3.0's pilot patterns are a large parameterized family
(D_x ∈ {3,4,6,8,12,16,24,32}, D_y ∈ {2,4}) selected by L1. Getting scattered-pilot
indexing wrong produces an equalizer that *almost* works — plausible constellations,
garbage output — which is expensive to diagnose. Build the pilot-position
generator as a separately unit-tested function with the A/322 tables as vectors.

### M4 — LDPC → BCH → ALP → IP
**Difficulty 4.** The data-path FEC: LDPC (16200/64800 bits, 12 code rates), BCH
outer, bit interleaving, constellation demapping (QPSK…4096NUQAM, plus non-uniform
constellations), then ALP (A/330) framing to recover IP packets.

**Biggest risk:** raw throughput. A 6 MHz ATSC 3.0 channel is ~25 Mbit/s and
64800-bit LDPC at 50 iterations is genuinely expensive; a naive Python decoder
will run thousands of times slower than real time. This milestone is where the
project either gets a C extension / numba / GPU path or becomes offline-only.
Per the standing GPU law, any GPU path must stay optional with a CPU fallback.
Recommend: prove correctness offline on a short capture first, optimize second.

### M5 — LLS / SLT service list + ESG (the TV guide)
**Difficulty 2.** Once IP packets exist this is mostly parsing. Low-Level
Signaling arrives on a well-known multicast address (224.0.23.60:4937); the SLT
is XML and yields the service list — the ATSC 3.0 equivalent of a channel scan.
ESG (A/344) is more XML.

**Biggest risk:** low. The main hazard is assuming the SLT is always present in a
short capture — LLS repeats on its own schedule, so capture length matters more
than cleverness.

### M6 — ROUTE/DASH → fMP4 → player
**Difficulty 3.** ROUTE (A/331) carries DASH segments as LCT/ALC objects;
reassemble to fMP4 and hand to ffmpeg/mpv. The fleet already knows how to drive
ffmpeg and mpv from STVT.

**Biggest risk:** **A3SA encryption.** Many US NextGen TV broadcasts encrypt the
media with A3SA DRM. This is explicitly **out of scope**: encrypted services will
be detected, enumerated and **labeled as locked** — never attacked, circumvented,
or decrypted. The honest deliverable for an encrypted service is "present,
identified, locked", not video. Expect a meaningful fraction of real services to
land there, and design the UI to say so plainly rather than looking broken.

### M7 — Fold into STVT as an HD/NextGen mode
**Difficulty 3.** Scanner learns to label ATSC 3.0 carriers (M1 already gives a
definitive yes/no, replacing the `tv_tuner.py` flat-OFDM heuristic), the DVR
learns a second decode path, the panel gets an ATSC 3.0 lane.

**Biggest risk:** the single-tuner rule and the launcher law. STVT's chain
assumes one SDR and one decode path; adding a second PHY must not let the daemon
stack tuners. Also portability — the community's #1 barrier is already
portability, and a heavyweight new dependency would make that worse.

---

## What M1 actually proves, and what it does not

- **Proves:** if a bootstrap is present at ≥ −15 dB SNR (isolated) or within
  ~12 dB below a co-channel carrier, this tool finds it, and reads the system
  version, EAS wake-up bits, frame cadence, bandwidth, post-bootstrap sample rate
  and preamble structure — 100% of the time with every signaling bit correct.
  Zero false alarms in 100 Gaussian-noise trials, 100 synthetic 8-VSB trials, and
  3 real off-air 8-VSB captures. The assumption-free structural detector alone
  reaches −12 dB.
- **Does not prove:** anything about the *current* state of the air. The archived
  captures swept in `SWEEP_RESULTS.txt` are from May–July 2026 and carry 8-VSB
  pilots; they answer "was RF27/RF36 ATSC 3.0 back then" (no), not "is it now".
  That needs the fresh capture specified in `README.md`.

## Bugs this milestone surfaced (kept as a record — each was a real defect)

1. **Circular cross-correlation sign.** The peak of `ifft(fft(rx)·conj(fft(ref)))`
   lands at −M, not M. Caught because G4 compares decoded fields against the
   *input* fields rather than checking self-consistency.
2. **Brittle `shift0 == 0` gate.** Discarded perfectly good decodes over ±1
   sample of timing jitter. A/321 encodes the bits *differentially* precisely so
   timing error cancels; the correct test is a tolerance, and `shifts[0]` is
   actually a free readout of residual timing error.
3. **Self-contaminated noise floor.** The background estimate excluded only
   ±½ symbol, so the bootstrap's own partial alignments raised its own floor.
   This compressed `peak_ratio` at *high* SNR and made the detector
   non-monotonic — worse at 0 dB than at −3 dB. Widening the exclusion to the
   full 4-symbol footprint improved sensitivity *and* specificity at once.
4. **Absolute-threshold detection statistic.** Thresholding the normalized
   correlation caps detection near 0 dB no matter how much processing gain is
   applied, because that statistic saturates at SNR/(1+SNR). The decision has to
   be made on peak-to-background ratio.
5. **`resample_poly(a, 6144000, 43048951)`** built a multi-million-tap filter and
   hung. Always rational-approximate a resample ratio first.
6. **Deep candidates appended after truncation.** `--deep` silently did nothing
   because its candidates went behind 12 structural ones and were cut. Caught
   only because deep and non-deep results were *byte-identical* — a suspicious
   equality is evidence.

## Assumptions ledger

See the module docstring in `atsc3/bootstrap.py` for A1–A5 with risk ratings.
The one that matters is **A1 (PN LFSR wiring)** — the only reading in A/321 that
a figure leaves genuinely ambiguous. `--pn-sweep` exists precisely so a real
off-air capture settles it instead of us guessing. A2 (BCA part-B time offset)
was initially uncertain and is now settled three independent ways; the residual
freedom is a constant phase on part B, which no detector here ever looks at.

---

# M2 results (2026-08-05) — real air, RF33 + three cross-checks

Everything below is offline work on captures already banked in `lab/`.
No radio was opened. Tools: `lab/m2_stage_a.py`, `lab/m2_stage_b.py`,
`lab/m2_pilots.py`, `lab/m2_stage_c.py`.

## Stage A — bootstrap signaling, and the A1 verdict

`m2_stage_a.py hit_rf33.cs16 --rate 8e6 --pn-sweep` found **32 bootstraps in
8 s**, decoded all 32, and every field was identical in all 32:

    sym1 00001100   sym2 00000010   sym3 00011011
    ea_wake_up            0,0  (no emergency alert)
    min_time_to_next      3    -> 200 ms  (a LOWER bound; see below)
    system_bandwidth      0    -> 6 MHz
    bsr_coefficient       2    -> post-bootstrap fs = 6.912 Msps
    preamble_structure    27
    abs cyclic shifts     [0, 68, 96, 244]   timing residual 0 every time

Frame period measured bootstrap-to-bootstrap: **247.111 ms, sd 0.0001 ms**
over 31 intervals. Consistent with min_time_to_next = 200 ms as a floor.
Fine CFO 471-477 Hz (0.80 ppm at 587 MHz), stable across the record.

Cross-channel (3 s each) - four independent transmitters:

| RF | frame period | preamble_structure | bw | bsr | EAS |
|----|--------------|--------------------|----|-----|-----|
| 25 | 242.518 ms   | 10                 | 6 MHz | 2 | 0 |
| 29 | 249.704 ms   | 27                 | 6 MHz | 2 | 0 |
| 30 | 249.703 ms   | 25                 | 6 MHz | 2 | 0 |
| 33 | 247.111 ms   | 27                 | 6 MHz | 2 | 0 |

Three distinct `preamble_structure` values and three distinct frame periods,
each rock-stable within its own channel. A decoder reading noise would give
either identical or unstable values; it gives neither.

### ASSUMPTION A1 (PN LFSR wiring) - SETTLED by real air

    spec            best mean peak ratio  114.54   detected
    spec_revseed                            3.38   not detected
    recip                                 114.54   detected
    recip_revseed                           3.38   not detected

**Verdict: the spec reading (out = r_0, feedback into r_15, seed MSB = r_15)
is correct.** 114.5 vs 3.4 is not a close call, and it is arbitrated by a real
transmitter, not by our own synthesis.

**But `--pn-sweep` was weaker than it looked, and that is worth recording.**
`spec` and `recip` returned *byte-identical* numbers because the A/321
polynomial x^16+x^15+x^14+x+1 has a PALINDROMIC tap set {0,1,14,15}, so
`_TAPS_RECIP` equals `_TAPS_SPEC`. Worse, a shift-left LFSR with tap set T is
provably identical to a shift-right LFSR with tap set mirror(T) on the
reversed register - and with T palindromic that collapses variant 3 onto
variant 2 and variant 4 onto variant 1. So the four "plausible wirings" span
only **two** distinct sequences. That is not a bug in the answer (the space
really is binary here) but the mitigation advertised 4 hypotheses and tested 2.
A suspicious equality is evidence - same lesson as bug #6 above.

## Stage B — the grid, MEASURED (no spec bits used)

The measurement folds the cyclic-prefix autocorrelation modulo a candidate
symbol pitch, with ILLEGAL FFT sizes and ILLEGAL guard intervals searched
alongside the legal ones as a null.

**First attempt was wrong and the null caught it**: peak-PICKING on the raw CP
profile reported GI = 1521, which is on no A/322 menu. Folding ~100 symbols
instead of picking peaks gave an unambiguous answer.

**A bigger trap, found by looking at absolute numbers rather than ratios:**
RF33's upper neighbour RF34 is ATSC 1.0, and its unmodulated 8-VSB pilot lands
3.31 MHz above RF33 centre - *inside* the 6.912 Msps passband, and measured at
**7 dB ABOVE** the wanted ATSC 3.0 signal. Being a pure carrier it correlates
at EVERY lag, so it simultaneously crushed the CP peak and lifted the
correlation floor 7x. Band-limiting to 2.95 MHz fixed both:

    CP fold peak/median   3.1  ->  40.0
    CP peak (normalized)  0.65 ->  0.999

Measured frame structure of RF33 (all from the signal):

    bootstrap        4 x 3072 @ 6.144 Msps = 13824 @ 6.912 Msps = 2.0000 ms
    subframe 0       FFT  8192, GI 1536, pitch  9728, 36 symbols = 50.6667 ms
    subframe 1       FFT 16384, GI 1536, pitch 17920, 75 symbols = 194.4444 ms
    -----------------------------------------------------------------------
    total                        1708032 samples @ 6.912 Msps = 247.1111 ms

**That total equals the Stage A bootstrap-to-bootstrap period (247.111 ms,
sd 0.0001 ms) exactly.** Two measurements that share no code close the frame
budget to the sample. This is the strongest cross-check in the milestone.

    subcarrier spacing   843.750 Hz (8K)      421.875 Hz (16K)
    symbol duration      1.4074 ms            2.5926 ms
    occupied bandwidth   5.8396 MHz           5.8392 MHz
    active subcarriers   6921 measured        13841 measured
                         (A/322 nominal 6913) (A/322 nominal 13825)

### CLAIMED vs MEASURED

| quantity | claimed (bootstrap bits) | measured (signal only) | verdict |
|---|---|---|---|
| system bandwidth | 6 MHz | 5.839 MHz occupied | agree |
| post-bootstrap fs | 6.912 Msps | CP lags land on exact powers of two, and the carrier count lands on A/322's 6913/13825, only at this rate | agree |
| preamble_structure | 27 | 8K / GI 1536 grid immediately after the bootstrap | consistent |

Note the frame carries **two subframes with different FFT sizes** - an 8K
subframe followed by a 16K one. Nothing in the bootstrap says that; it was
found by a head-to-head grid score in sliding windows, and the boundary lands
at 52.667 ms = exactly 36 x 9728 samples after the bootstrap.

### Pilots

Determined without knowing a single pilot VALUE, using
`M_d(k) = |SUM_l X_l(k) conj(X_{l+d}(k))| / SUM_l |X_l(k)|^2`. Pilots repeat
the same reference every D_y symbols, data does not; and M_d is immune to a
constant per-symbol phase rotation. The smallest d that lights up IS D_y.

    subframe 0 (8K):  d = 2,4,6,8,12 light up, d = 1,3,5 do not  -> D_y = 2
                      carriers that light up are spaced 4        -> D_x = 4
                      53 continual pilots
    subframe 1 (16K): d = 4,8,12 light up, d = 1,2,3,5,6 do not -> D_y = 4
                      best (D_x,D_y) = (4,4) at 3.12, best illegal control
                      (20,4) at 2.18                            -> SP4_4
                      101 continual pilots, ALL at rel-carrier == 0 (mod 2)

Two earlier pilot metrics were **discarded as biased** and the record is kept:
a "best comb phase" power ratio is biased upward in proportion to D_x for ANY
heavy-tailed spectrum and happily "found" D_x = 32 in pure data; and a
straight amplitude-boost test could not separate legal from illegal D_x at all.
Only the repeat-coherence statistic separated. The D_y control values are also
known-weak: any multiple of the true D_y trivially also lights up, so D_y is
read off the LAG SCAN (smallest lit lag), not off the control comparison.

## Stage C — climbed, gated, and stopped short. Honestly.

Reached and gated:

* cells extracted on the measured grid;
* continual pilots found blind, then used to strip each symbol's common phase
  error AND residual timing slope (without this the two scattered-pilot
  residue groups sit at different phases and the channel is not smooth in k);
* scattered-pilot channel estimate with **no A/322 pilot tables**: the pilot
  reference is +-1 and the channel is smooth, so the sign pattern separates
  from the channel by phase continuity along the lattice. Pilot coherence
  0.95, pilot-referenced SNR 31 dB on RF33;
* `selftest_equalizer()` - synthetic OFDM, known QPSK, SP4_2 pilots with a
  random +-1 reference, 3-tap multipath (47 dB |H| range), 30 dB SNR - comes
  back through the SAME code at **23.9 dB QPSK MER**;
* `min_sum_decode()` - normalized min-sum BP - reaches zero bit errors on a
  synthesized rate-1/5 QC-LDPC. **This validates decoder MACHINERY ONLY.** It
  is our encoder tested against our decoder and says nothing about A/322.

**That gate earned its keep immediately.** It failed at 2.7 dB on the first
run and exposed a real bug that looked completely benign on air: the 4th-power
blind phase estimator has a **45 degree bias**. For QPSK and for every square
QAM, `sum(z**4)` is real and NEGATIVE when the constellation is already
correctly oriented, so `angle(mean(z**4))/4` returns 45 degrees and
"correcting" by it rotates a good constellation onto the axes where it matches
nothing. On real air this produced a plausible-looking noisy cloud and cost
~18 dB of MER with no symptom whatsoever. Two other candidate causes were
eliminated first, each by direct measurement rather than assumption: the blind
sign recovery is fine (it is correct up to ~20 piecewise flips, harmless
because every square constellation is 180-degree symmetric), and the lattice
interpolation is fine (linear interpolation reconstructs the true channel to
-43 dB).

Not reached: **the equalized RF33 cells do not resolve into a constellation.**
Best-fit MER rises monotonically with order - QPSK 4.4 dB, 16QAM 9.8 dB,
64QAM 13.6 dB, 256QAM 15.9 dB - which is the signature of an unresolved cloud,
not of 256QAM. In particular the first preamble symbol is not QPSK, so the
cell mapping assumed for it is wrong. Shortening the channel-estimation window
to 2-16 symbols made it worse, so time variation is not the cause.

**The wall, precisely.** The decoder ALGORITHM is no longer the blocker. What
is missing is published TABLES this repo does not have:

  (a) the A/322 cell-mapping order for the preamble - which cells of which
      preamble symbol carry L1-Basic - and the frequency interleaver
      permutation. This alone explains the unresolved preamble constellation
      and must be fixed BEFORE any LDPC work is meaningful;
  (b) the L1-Basic LDPC parity-check address tables (16200-bit mother code,
      rate 3/15) plus its shortening and puncturing schedule;
  (c) the L1 scrambler and bit interleaver.

None of (a)-(c) is derivable from the signal.

## Cross-check against the HDHomeRun ground truth

The user's HDHomeRun reports six unencrypted ATSC 3.0 services locally.
The virtual channel list is the user's location and is kept in gitignored
`lab/ground_truth.txt`, not here. M2 cannot name a
service - service names live in the SLT, which is M5, four milestones of FEC
and IP reassembly away. What M2 *can* say and does: four RF channels
(25/29/30/33) each carry a genuine, stable, 6 MHz ATSC 3.0 physical layer,
which is consistent with 6 services spread over 4 RF carriers. That is a
consistency check, not a confirmation, and it is recorded as such.

## Assumptions added in M2

    B1  the A/322 guard-interval menu {192,384,512,768,1024,1536,2048,2432,
        3072,3648,4096,4864} in post-bootstrap samples.  RISK LOW -- used to
        enumerate and to LABEL, never to constrain; illegal GIs are searched
        alongside and the winner beat them by 2.6x.
    B2  FFT size is one of {8192,16384,32768}.  RISK LOW -- same treatment;
        illegal lags 4096/6144/12288/24576 searched as controls and lost.
    B3  A/322 active-carrier counts for 6 MHz normal mode: 6913 / 13825 /
        27649.  RISK LOW -- used ONLY as a cross-check, and the measurement
        landed within 16 carriers of both without being told them.
    P1  the channel is static over the pilot analysis window (tens of ms).
        RISK LOW -- violating it lowers the coherence metric uniformly; it
        cannot manufacture a lattice.
    P2  carrier indices are counted from the CENTRE of the measured active
        band.  RISK LOW -- a wrong origin shifts the reported residue, not
        D_x or D_y.
    C1  the scattered-pilot lattice origin is found by maximizing measured
        pilot coherence over all D_x*D_y shifts.  RISK LOW -- but note the
        first version searched only D_x shifts, which cannot reach a swapped
        symbol parity, and silently returned coherence 0.269 instead of 0.946.

## Next step for M3/M4 - sharpest first

The single highest-value move is **not** more DSP. It is to obtain the A/322
tables (the standard is public at atsc.org) and load them as DATA into this
repo, in this order:

  1. **preamble cell mapping + frequency interleaver** (A/322 sec. 8.2, 9.2).
     This is what is blocking the constellation right now, and it is testable
     the moment it lands: a correct mapping turns the preamble cloud into
     clean QPSK, and `selftest_equalizer()` already proves the chain will show
     it if it is there.
  2. the scattered/continual pilot REFERENCE sequence - our blind sign
     recovery can then be checked against it, which independently validates
     both. Our measured SP4_2 / SP4_4 and our 53/101 continual-pilot positions
     are testable predictions against those tables today.
  3. only then the LDPC address tables. The BP decoder is already written and
     gated.

Second-highest: run `m2_stage_b.py` over RF25/29/30 as well. Their
`preamble_structure` values differ (10, 27, 25), so measuring their grids
gives four independent (preamble_structure -> grid) data points and lets the
A/322 table entry be *checked* rather than trusted when it arrives.

---

# M3 results (2026-08-05) — the A/322 tables, and L1-Basic off real air

Offline work on captures already banked in `lab/`. **No radio was opened.**
Tools: `lab/m3_spec.py`, `m3_preamble.py`, `m3_freqint.py`, `m3_cells.py`,
`m3_replicate.py`, `m3_profile.py`, `m3_ldpc.py`, `m3_l1basic.py`,
`m3_referees.py`, `m3_crc.py`, `m3_descramble.py`, `m3_scrambler.py`, plus the
extracted tables `spec_ldpc.py`, `spec_bicm.py`, `spec_l1syntax.py`.

M2 ended by naming its own blocker: *published tables, not algorithms.* That
diagnosis is now confirmed. The tables were fetched from atsc.org
(A/322:2026-04, cross-checked against A/322:2024-04) and the preamble cloud
resolved on the first run with them.

## Step 1 — the tables, and how each is defended

The A/321 math-font digit-shift hazard means a mis-extracted table looks like a
decoder bug forever. So **no number is trusted because the PDF said so.** Each
is pinned by a closed arithmetic identity against a *different* table:

| identity | what it pins |
|---|---|
| NoC steps by exactly Cunit (96/192/384) | Table 7.1 |
| occupied BW == NoC x fs/Nfft to 6 decimals | Table 7.1 |
| CP8 == ceil(CP32(4k)/4), CP16 == ceil(CP32(2k)/2) | Tables D.1.1/2/3 |
| CP32 mirror-symmetric about 27648 | Table D.1.1 |
| common-CP count == #CP inside the reduced range (48/48/47/46/45) | Table 8.4 |
| data cells == NoC - preamble pilots - CPs, every 8K/16K row | Table 7.2 |
| amplitude == 10^(dB/20), every row | Table 8.6 |
| the reference PRBS reproduces the spec's OWN printed 24-bit vector | Sec 8.1.2 |
| all 7 printed cell counts reproduce from the puncturing arithmetic | Tables 6.17/6.23/6.24 |
| GF(2^14) minimal polys of alpha^1..alpha^23 reproduce g1..g12 in order | Table 6.3 BCH |
| Tables 6.20/6.21/6.22 are exact permutations of their ranges | bit interleavers |
| 1819 Annex A.2 integers identical across 4 witnesses, sum 6754166 | LDPC matrices |

`m3_spec.verify()` runs 13 of these on every invocation: **13/13 pass.**

**Two premises the extraction corrected, rather than accepting:**

1. Table A.2.1 is rate **2/15**, not 3/15 — the L1-Basic code is **A.2.2**.
2. The Type A row count is **Kldpc/360 + Q1**, not Kldpc/360. Rate 3/15 has 12
   rows, not 9. The extra Q1 rows exist because encoding steps (vii)-(viii)
   feed the M1 dual-diagonal parity back in as further information. Confirmed
   independently: in *every* Type A code the tail rows contain only addresses
   >= M1, exactly as that structure requires.

## Step 2 — PASS. The cloud is QPSK.

What M2 had wrong was three table lookups:

* **NoC.** A/322 7.2.5.1 — the *first* preamble symbol uses the **minimum** NoC
  for its FFT size, **6529** for 8K, not the data symbols' 6913.
* **Pilot pattern.** A/322 8.1.6.1 — preamble pilots are DY=1 at every
  `k mod DX == 0`, a quarter of all carriers. M2 used the data SP4_2 pattern.
* **Pilot values.** A/322 8.1.2 publishes the reference PRBS, so pilots are now
  known exactly instead of sign-recovered blind.

Carrier origin and Cred_coeff were **scanned, not assumed**, so wrong values are
built-in controls:

    Cred 0  NoC 6913   pilot coherence 0.0446
    Cred 1  NoC 6817                   0.0515
    Cred 2  NoC 6721                   0.0515
    Cred 3  NoC 6625                   0.0507
    Cred 4  NoC 6529                   0.9948   <- spec-predicted, shift +0

Geometry then closes exactly: 6529 - 1633 pilots - 45 CPs = **4851** = Table 7.2.

After the A/322 7.3 frequency de-interleaver, the first 484 cells in CELL order
(Table 6.17: L1-Basic Mode 3, QPSK):

| subset | QPSK MER | 4th-power statistic |
|---|---|---|
| **L1-Basic cells 0..483** | **20.65 dB** | **0.966 @ -178.6 deg** |
| same, reverse wire permutation | 6.04 dB | 0.214 |
| L1-Detail cells 484.. | 5.49 dB | 0.064 |
| last 484 cells (control) | 4.39 dB | 0.113 |
| random 484 cells (control) | 6.10 dB | 0.167 |

**+14.5 dB over the random-subset control on a subset we do not get to pick.**
Replicated over 10 frames spanning 6 s: 17.0-26.3 dB, 4th-power 0.92-0.99.

Annex H Table H.1.1 validated four ways. `preamble_structure` 27 -> 8K/GI1536/
DX4, which M2 had **measured blind**. RF25's `preamble_structure` 10 -> 8K/
GI512/**DX6** is a genuine out-of-sample hit: never measured before, and its
geometry closes at 5395 data cells = Table 7.2 exactly.

The constellation resolves only on RF33. RF25/29/30 are 3-second cross-check
grabs at pilot SNR 12.0/6.9/0.4 dB vs RF33's 23.5 — **capture-limited, not
decoder-limited.** Their *tables* still validate.

## Step 3 — L1-Basic bits recovered and proven; field values blocked

Three referees, each given controls that must fail:

**Referee 1 — BCH syndrome zero** over the 368 systematic bits (1 in 2^168).
Covers QPSK demap, Annex C.1.1 bit-to-IQ, 6.5.2.10 de-interleave, the frequency
de-interleaver and cell order. Controls: I/Q swapped FAIL, reverse FI direction
FAIL, random bits FAIL.

**Referee 2 — LDPC convergence**, 12960 simultaneous checks:

    standard / direct        CONVERGED, 4 iterations
    standard / inverted      729 checks unsatisfied
    swapped  / direct        864 unsatisfied
    swapped  / inverted      959 unsatisfied
    wrong shortening pattern 374 / 636 unsatisfied

This settles ASSUMPTION L1 (the parity interleave index order, whose subscripts
the PDF math font destroys) **from the air rather than by choice.**

**Referee 3 — CRC-32, and it needs no scrambler.** A/322 6.1.2.2's CRC has a
fixed all-ones init, so it is affine, so the unknown scrambler mask **cancels in
a difference**: `g(r_i XOR r_j) = 0` for any two valid blocks. Over 8 frames
(all distinct, Hamming 17-20 apart): **28 of 28 pairs give zero.** Controls:
0/2000 random vectors, 0/2000 random same-weight deltas, and 0/28 under
CRC-32/MPEG-2, under x^26+x^23+x^22+1, and under the A/322 polynomial with the
+1 dropped.

**A suspicious equality, run down rather than accepted.** The first sweep
reported "BCH PASS" for all eight variants including ones leaving 800+ checks
unsatisfied. Cause: the 368 Nouter bits are *systematic* and arrive at ~20 dB
MER, so LDPC-decoded vs raw hard decisions agree **368/368, zero bits
corrected** — BCH never sees the LDPC. The two referees cover disjoint halves
of the chain. A control of mine was also wrong and is kept: a shortening pattern
disjoint from the true one reads back all zeros, and the all-zero vector is a
legal BCH codeword, so it "passed" for the wrong reason entirely.

**Field layout confirmed from the air.** Across 8 frames the only bits that vary
are 11-14 and **168-199 — exactly the 32-bit L1B_crc field**, precisely where
the extracted Section 9.2 layout puts it, with nothing between 15 and 167
varying. A raw (unscrambled) parse gives illegal values (pilot pattern 28 =
Reserved, GI 3072 illegal at 8K), confirming the block really is scrambled.

## The wall, precisely: the A/322 5.2.3 scrambler

The polynomial, the 0xF180 seed, the procedure and a printed 24-bit test vector
all extract cleanly. **Which eight register stages are D7..D0 does not** —
Figure 5.6 is vector art.

* **The taps ARE recoverable from glyph columns.** Each D label sits just right
  of the stage box it taps, giving `D0..D7 = X1, X3, X4, X7, X11, X12, X13, X14`,
  which reproduces the spec's printed **first byte exactly, 8/8** (1 in 256).
  Reading from the stage to the right instead gives the wrong byte.
* **But then the printed vector is provably unreachable.** Under a single shift,
  X4 of the next state is whatever X3 held (=0), so byte 1's D2 bit is *forced*
  to 0 while the printed vector needs 1. Feedback writes only the end stage and
  cannot touch X4. **Figure 5.6 and the printed 24-bit vector are mutually
  inconsistent** — one carries an erratum and the text cannot say which.
* **Three independent searches, all empty.** (a) 102 960 hypotheses vs the
  printed vector -> 0. (b) 68.7 billion (16 configs x 16^8 tap assignments, all
  orderings and repeats) vs the air's 32-bit CRC, solved by meet-in-the-middle
  -> 12 survivors against ~16 expected *by chance*, none reproducing the
  independently measured FFT size / GI / pilot pattern / subframe count;
  correctly rejected. (c) constraint propagation from **69 mask bits derived
  from the air** — 47 from L1B_reserved (A/322 3.2.1: "The ATSC default value
  for reserved bits is 1") plus 22 from M2-measured fields — across both the
  2024 and 2026 field layouts, reserved all-ones and all-zeros, three
  symbol-count readings -> **0**.

**Consequence, stated plainly:** the 200 L1-Basic bits are recovered and proven
correct, but they stay scrambled, so the field *values* are not readable and
**Step 4 (L1-Detail / the PLP configuration) was not reached.** No guessed
scrambler is shipped.

## Assumptions added in M3

    S1  DC is an ordinary active carrier; NoCmax odd, centre = (NoCmax-1)/2.
        RISK LOW -- not assumed blind, the origin scan picks it by 20x margin.
    S2  the guard interval is a cyclic PREFIX.  RISK LOW -- M2 measured the CP
        autocorrelation at exactly this geometry.
    S3  preamble pilot amplitude = Table 8.6's "Equivalent Amplitude".
        RISK LOW -- a wrong global scale divides out in the equalizer.
    F1  FI wire-permutation direction (Tables 7.12-7.14 print two bit-position
        rows with no assignment operator).  SETTLED BY DATA: 20.65 vs 6.04 dB.
    F2  the FI MSB toggle argument.  SETTLED BY CONSTRUCTION: toggle=l cannot
        produce a bijection onto 0..Ndata-1 (1510 of 4851 addresses).
    L1  the LDPC parity interleave index order (6.1.3.1 steps vi/viii, whose
        subscripts the PDF destroys).  SETTLED BY THE AIR: the alternative
        leaves 729-959 of 12960 checks unsatisfied.

## Next step for M4 — sharpest first

1. **Resolve the scrambler.** Everything else is built and gated behind it.
   Cheapest routes: read Figure 5.6 by eye from the PDF; or compare against a
   known-good ATSC A/53 randomizer, since A/322 5.2.3 reuses A/53's polynomial
   and seed; or file the Figure-5.6-vs-test-vector inconsistency with ATSC via
   https://www.atsc.org/feedback/ — it looks like a genuine erratum.
2. Then L1-Detail follows immediately: same machinery, size signalled by
   L1-Basic, and the QPSK cliff already measured at cell ~1320 predicts
   `L1B_L1_Detail_total_cells` ~= 836, which the decoded field must match.
   That cross-check is set up and waiting.
3. Only then the payload PLPs and the Annex C NUC point tables.

---

# M4 results (2026-08-06) — the scrambler, L1-Detail, and the payload path

Offline work on captures already banked in `lab/`. **No radio was opened.**
The one bounded live capture this session was permitted was deliberately NOT
taken; see "On not taking the radio" at the end.

Tools: `lab/m4_scrambler.py`, `m4_l1detail.py`, `m4_data.py`, `m4_cells.py`,
plus the extracted table module `spec_nuc.py`.

## Step 1 — the A/322 5.2.3 scrambler. M3's erratum call was WRONG.

M3 concluded that Figure 5.6 contradicts the spec's own printed 24-bit test
vector and called it an erratum. That conclusion was wrong, and the correction
is more valuable than the original finding.

M3's proof was: under a single shift, stage X4's next value is whatever X3 held
(= 0), so byte 1's D2 bit is forced to 0 while the printed vector needs 1;
feedback writes only the end stage and cannot touch X4. **That argument silently
assumes a FIBONACCI register.** A/322 5.2.3's register is **GALOIS** — the
feedback bit is XORed *into* the tapped stages as the contents move — so X4's
next value is X3 XOR feedback, it is not forced, and the contradiction
evaporates. The wall was an assumed convention, not a spec defect.

**Record the general lesson: "the spec is wrong" is the most expensive
conclusion available and must be the last one reached, not the first.**

Enumerating the 16 register conventions against the printed vector leaves
exactly one:

    shift DOWN (X_i <- X_{i+1}), GALOIS form,
    feedback taps = the RECIPROCAL exponents [3,4,5,9,10,13,15,16],
    seed 0xF180 loaded with X1 = MSB

The 24 printed bits give only a 3-bit history per output, so they underdetermine
which stages feed D7..D0: 16875 assignments fit, 2160 of them with 8 distinct
stages. **The air pins them.** A/322 6.1.2.2's CRC-32 has a fixed all-ones init,
so it is affine, so `g(r XOR mask) = K` is 32 hard bits, and 16875 candidates
against 32 bits gives 3.9e-06 expected false survivors.

    SURVIVORS: 1.    D0..D7 = X16, X14, X13, X10, X6, X5, X4, X3

Three confirmations the search could not have fitted:

1. **The mirror check.** Relabel j = 17 - i and the solved output stages become
   `D0..D7 = X1, X3, X4, X7, X11, X12, X13, X14` — byte for byte, M3's
   glyph-column reading of Figure 5.6, recovered from the air with no knowledge
   of the figure. **The figure's taps were right all along.**
2. **Eight independent measurements, 8/8.** The descrambled fields against
   quantities M2 measured off the waveform before any A/322 table existed in
   this repo: FFT 8K, GI 1536, SP4_2, Cred 0, 2 subframes, 36 symbols
   (34 -> 35 data + 1 preamble), 1 preamble symbol — plus `L1B_reserved` = all
   47 bits set, the A/322 3.2.1 default, on bits the CRC equation does not
   constrain at all.
3. **Internal closure.** `L1B_L1_Detail_size_bytes` = 64 and `fec_type` = Mode 3
   run through the Table 6.24 puncturing arithmetic predict exactly 880 cells;
   the separate field `L1B_L1_Detail_total_cells` reads 880.

Per frame, only `L1B_time_offset` and its CRC move, and the offset advances by
exactly 12288 every frame.

## Step 2 — L1-Detail, and the PLP configuration of RF33

Referees: LDPC rate 6/15 **converged in 5 iterations with 0 of 9720 checks
unsatisfied**; BCH syndrome zero over the 680 Nouter bits; A/322 9.3's
`L1D_crc` passes over the descrambled 480 + 32 bits — and unlike L1-Basic's,
that CRC is evaluated *with* the scrambler, so it re-tests Step 1 on completely
different data. `L1D_reserved` is 7 trailing bits, all ones.

All seven controls fail to converge: cells shifted +-40, L1-Basic's own cells,
the tail of the symbol, wrong shortening pattern, wrong group-wise pattern,
group-wise inverted.

**Four extraction bugs found and fixed** in `spec_l1syntax`'s L1-Detail tree,
all caused by the PDF printing `}else {` on its own line. Run as extracted the
parse derails right after the first PLP (PLP 2 came back with
`scrambler_type` = "Reserved", `fec_type` = "Reserved" — nonsense that was easy
to *see* precisely because the semantics tables label reserved values).
Restored: the `TI_mode=01` -> `L1D_plp_CTI_fec_block_start` branch; the
`TI_mode=10` -> six HTI fields branch; the `else` of `if (plp_layer = 0)` ->
`L1D_plp_ldm_injection_level`; and the gate on the second `for i` loop, which
A/322 9.3 footnote 2 calls the "L1D_version >= 2 loop" and which is **absent
here because this transmitter signals L1D_version = 1**. That last point also
dissolves the one ambiguity L1-Basic could not settle on its own:
`L1B_first_sub_mimo_mixed` abuts the all-ones reserved run, but at
L1D_version 1 it gates nothing and both readings parse identically.

### The configuration — BSID 540, 3 PLPs over 2 subframes

    subframe 0   8K / GI1536 / SP4_2 / 35 data symbols + 1 preamble
                 SBS first and last, 127 SBS null cells, FI enabled
      PLP  0   Core layer, LLS CARRIER, cells 0..199799,
               BCH + 16K LDPC, 64QAM-NUC 11/15,
               HTI intra-subframe, 2 TI blocks x 73 FEC blocks
      PLP 16   Core layer, cells 199800..207899,
               BCH + 16K LDPC, QPSK 2/15, HTI, 1 TI block, 0 FEC blocks
    subframe 1   16K / GI1536 / SP4_4 / 75 symbols
                 SBS first and last, 1609 SBS null cells, FI enabled
      PLP  1   Core layer, cells 0..947699,
               BCH + 64K LDPC, 256QAM-NUC 11/15,
               HTI intra-subframe, 3 TI blocks x 116 FEC blocks

**Five more independent measurement matches, from a block M2 never saw.** M2
measured subframe 1 blind off the waveform: FFT 16384, GI 1536, D_x = 4,
D_y = 4, 75 symbols, ~13841 carriers. L1-Detail signals FFT 16K, GI6_1536,
SP4_4, `num_ofdm_symbols` 74 (= 75), Cred 0 (NoC 13825). Five for five.

**External clock check.** `L1D_time_sec` is TAI seconds at the first bootstrap
sample. Decoded 1785979209 TAI, minus the 37 s TAI-UTC offset, is
2026-08-06 01:19:32 UTC. The capture file's own mtime puts the start of the 8 s
record at 01:19:31.7. **The decoded transmitter clock lands +0.3 s into the
record** — a 32-bit field carried through bootstrap, OFDM, frequency
de-interleave, LDPC, BCH, descramble and syntax parse, agreeing with the
filesystem clock.

### The ~836-cell prediction, scored honestly — it MISSED

M3 read the QPSK MER cliff at cell ~1320 and predicted
`L1B_L1_Detail_total_cells` ~= 836. The decoded field is **880**, so L1-Detail
ends at cell 1364: the prediction was **44 cells (5.0%) LOW**. A fine
re-measure shows why — MER holds ~20 dB through the 1290-1350 window and only
collapses in 1310-1370, so the cliff genuinely sits near 1364 and M3's coarse
read landed early. The prediction was in the right neighbourhood but it is not
a hit. The FEC is the real arbiter: LDPC will not converge with the boundary
off by even 40 cells, and it converges at exactly 880.

## Step 3a — the payload constellation

A/322 Annex C NUC tables extracted (`lab/spec_nuc.py`), never hand-typed.
78 identities: both witnesses (2026-04 and 2024-04) give identical tables for
all six definition tables, and the mean power of every one of the 36
reconstructed constellations is 1.0 — to 7e-4 for the 16-NUCs, 3e-5 for the
64-NUCs, 1.5e-5 for the 256-NUCs. **That ordering is itself evidence the parse
is clean:** it is exactly the 4-decimal quantisation error from averaging 4, 16
and 64 distinct magnitudes.

Subframe 0 data symbols demodulated with A/322 8.1.3 scattered pilots. The
lattice-phase reading is scanned, not assumed: `l` counting from the first data
symbol gives mean pilot coherence 0.9922, the off-by-one control gives 0.0324.

**The comparison is deliberately unbiased.** M2 was burned comparing across
constellation *orders*, where more points always fit better and a noise cloud
"improves" monotonically. Every candidate here has exactly 64 points.

                       MER dB   centroid RMS   occupancy H   chi2/dof
      NUC_64_11/15      18.78         0.0229        5.9991        3.6  <- signalled
      NUC_64_12/15      18.09         0.0394        5.9962       14.5
      NUC_64_10/15      17.29         0.0523        5.9838       64.2
      uniform 64QAM     15.22         0.0887        5.9148      275.4
      NUC_64_2/15        8.19         0.2643        5.4763     1986.9

Three statistics measuring three different things, one winner, and it is the
one L1-Detail signalled. The sharpest is the occupancy chi-square, which does
not depend on SNR at all: an LDPC-coded, scrambled payload must use all 64
points equally often, and a wrong table cuts the plane in the wrong places.

**Link-budget note: 18.8 dB MER is about the AWGN threshold for 64QAM 11/15.
PLP 0 on this capture is AT the cliff, not comfortably above it.** 2-D pilot
interpolation across the DY = 2 symbol pair was tried and gained nothing
(18.76 vs 18.78 dB), so this is link-limited, not estimator-limited.

## Step 3b — the cell map, and A/322's gift of a referee

A/322 7.2.6.5 has been sitting in the spec the whole time and only becomes
usable once the scrambler is solved: every data cell a Subframe does not use is
filled with `Re{d_i} = 1 - 2.s_i`, `Im{d_i} = 0`, where `s_i` is the i-th bit of
the **Section 5.2.3 scrambling sequence** and `i` is the cell's index within the
Subframe. **Cells nobody is using carry known values, and no FEC is required to
read them.**

The dummy cells were found by looking, not by assuming: sweeping `|Im(z)|`
across all 35 data symbols, symbols 0-33 look like ordinary modulation and
symbol 34 is 71% real-valued at `|Re|` = 1.0009. Their signs were then
cross-correlated against the scrambler over every possible cell-index origin
and every open reading of the frequency de-interleaver:

    agreement   origin   FI direction   FI symbol index
       1.0000   206400   forward        35        <== EXACT, 3572 of 3572
       0.0885    45160   reverse        30
       0.0862   232790   reverse        34
    chance level for 3572 signs swept over 260000 origins ~ 0.084

Reproduced on two frames, identical origin. Controls all at chance: a scrambler
with different output taps 0.1652, M3's Fibonacci reading 0.0739, a random
sequence 0.0756, no frequency de-interleaving 0.0773, normal-symbol pilot
pattern instead of SBS 0.0665.

**One correlation settled four things at once:**

1. **the scrambler, a third time** — now on 3572 fresh bits from around index
   206400 of the sequence rather than the first 200, with zero errors;
2. **A/322 7.3 frequency de-interleaving of DATA symbols.** ASSUMPTION F1 (wire
   permutation direction) is settled for data symbols as well as the preamble:
   *forward*. And the FI symbol counter **includes the preamble as symbol 0**,
   so data symbol 34 is FI symbol 35 — a reading nothing in the L1 work could
   have arbitrated;
3. **the cell-index origin**, measured rather than accumulated;
4. **the decoded `L1D_plp_start` / `L1D_plp_size` values**, via the predictions
   below.

**Table D.1.4, and why it matters.** A/322 8.1.4.1 says the additional continual
pilots exist "to ensure a constant number of data carriers in every data
symbol". With the 48 common CPs alone the count alternates 6000 / 5999. Add
Table D.1.4's single extra CP for 8K + SP4_2 — relative carrier **1732** — and
every normal data symbol has exactly **5999**. A closed identity that closes
only for the right table entry. SBS symbols have 5136.

**The PLP allocations became a falsifiable prediction and passed.** With the
origin measured, `plp_start`/`plp_size` say exactly which cell each modulation
change falls on:

    symbol 32, internal position 5398:  MER vs QPSK 4.60 dB before, 16.17 after
    symbol 34, internal position 1500:  real-valued fraction 0.042 before,
                                        0.999 after

Symbol 33 in isolation is pure QPSK at 20.94 dB — PLP 16 at QPSK 2/15, exactly
as signalled — while its neighbours 31 and 32 are 64QAM-NUC at 20.92 / 17.51.
And symbol 34 contains exactly **127** near-zero cells, which is
`L1D_sbs_null_cells` to the cell.

The preamble's spare cells were also verified directly: cells 0..483 and
484..1363 are QPSK (L1-Basic, L1-Detail) and cells 1364..4850 are 64QAM-NUC
11/15 at ~20.3 dB — PLP 0 data in the preamble symbol, exactly as A/322 7.2.5.2
requires.

### OPEN ITEM — the 190-cell preamble handover discrepancy

Walking the measured origin (206400) back through 33 normal symbols at 5999 and
one SBS symbol at 5136 puts symbol 0 at cell 3297, i.e. the preamble contributes
3297 available data cells to subframe 0's pool. But the preamble demonstrably
carries PLP 0 data from cell 1364 to 4850, which is **3487** cells.

**190 cells are unaccounted for**, and this is recorded rather than fudged.
Equivalently, the first SBS symbol has 4946 available cells where the pilot
arithmetic says 5136. Note 190 = 127 + 63 and the SBS null count is 127, which
is suggestive but not a derivation.

This does not block the payload chain — PLP 0's cell stream can be anchored at
its measured END (cell 199799 = symbol 32 internal position 5397) and counted
back 199800 cells — but it means something in the preamble-to-subframe cell
handover is still misunderstood, and **it should be run down before the payload
chain is trusted end to end**, because a 190-cell offset would poison the time
de-interleaver silently.

## Addendum — subframe 1, and what it says about the open item

The same dummy-cell test was run on subframe 1 (16K / GI1536 / SP4_4 / 75
symbols, PLP 1 at 256QAM-NUC 11/15). Three results:

**1. Subframe 1 demodulates cleanly.** Symbols 1, 2 and 66-73 give
21.7-22.9 dB MER against NUC_256_11/15 — the constellation L1-Detail signalled,
and slightly better than subframe 0's. (256QAM 11/15 also needs roughly this
much SNR, so PLP 1 is at its cliff too.)

**2. The dummy mechanism confirms at 16K.** The last symbol (74, an SBS symbol)
is **98.4% real-valued** and contains **1611 near-zero cells** against a
signalled `L1D_sbs_null_cells` of **1609** — two cells apart. Correlating the
8507 dummy signs against the scrambler:

    agreement   origin    FI direction   FI symbol index
       0.9951   422432    forward        74
       0.0575   379712    forward        77
    chance ~ 0.054

A real lock (18x chance), with ~21 sign mismatches out of 8507 attributable to
low-SNR cells. The A/322 7.3 frequency de-interleaver reading — direction
*forward*, and the symbol counter — therefore holds at 16K as well as 8K.

**3. And it sharpens the 190-cell open item into a diagnosis.** The scrambling
sequence has period 65535 states = 524280 bits, so the correlation fixes the
cell origin only modulo 524280. The plausible absolute candidate is 946712.
The pilot-count model predicts 949179. **They disagree by 2467 cells**, or about
34 per symbol.

Subframe 1 inherits no cells from the preamble, so this cannot be a
preamble-handover problem. Both discrepancies are therefore the same defect:
**the available-data-cell model is incomplete — too few continual pilots.**
Consistent with that, A/322's "constant number of data carriers in every data
symbol" identity CLOSES at 8K/SP4_2 with Table D.1.4's single extra carrier
1732, but does NOT close at 16K/SP4_4 with the two carriers (5768, 11452) as
extracted: the counts still alternate 12862 / 12861 / 12862 / 12862. The data
implies roughly 35 additional continual pilots per 16K symbol, not 2.

**So the next concrete task is not the time de-interleaver. It is to re-extract
Table D.1.4 (and check whether it has more rows than the two-line-per-pattern
layout suggests), using the constant-data-carrier identity as the gate — it is
a closed arithmetic test that either passes or does not, for every FFT size and
pilot pattern.** Until it passes, every cell index downstream is uncertain by
tens of cells per symbol, which would poison the time de-interleaver silently.
That is a far better place for the next session to start than three rungs of
unverifiable interleaver code.

## Where M4 stops, precisely

Reached and gated: L1-Basic values, L1-Detail / full PLP configuration, the NUC
tables, data-symbol OFDM demodulation and equalization, the payload
constellation, the frequency de-interleaver for data symbols, the subframe cell
map and the PLP boundaries.

**Not reached: payload BITS.** The remaining rungs for PLP 0, in order, with
what each still needs:

  1. **Cell de-multiplex to a PLP cell stream.** Mostly done — blocked only by
     the 190-cell open item above.
  2. **Time de-interleaver (HTI).** `L1D_plp_TI_mode` = 10, `inter_subframe` =
     0, `num_ti_blocks` = 1 (so NTI = 2), `num_fec_blocks` = 73,
     `cell_interleaver` = 0. Needs A/322 7.1.5's Twisted Block Interleaver
     geometry implemented. **This is the largest single piece of unbuilt work.**
  3. **Bit de-interleaver for data PLPs.** A/322 6.2: parity interleaver,
     group-wise interleaver (Tables 6.10/6.11), block interleaver
     (Tables 6.12/6.13). **None of these tables is extracted yet** —
     `spec_bicm.py` carries only the L1 versions.
  4. **LDPC 16200 rate 11/15** (the code itself IS already in `spec_ldpc.py`
     and the BP decoder is built and gated) + BCH outer + descramble.
  5. ALP depacketization -> IP -> LLS/SLT -> ROUTE/DASH -> fMP4.

There is no intermediate referee between rung 1 and rung 4: LDPC convergence is
the first thing that can say "yes". Three rungs must therefore be built
correctly before anything can be checked, which is exactly the situation M2 and
M3 warned about, and it is the main reason M4 stops here rather than shipping a
half-verified chain.

**And a link-budget caveat that will matter when it is built:** at 18.8 dB MER
PLP 0 sits at the 64QAM 11/15 threshold, so even a perfect implementation may
not converge on this capture. A failure at rung 4 will be ambiguous between
"chain wrong" and "SNR short" unless a higher-SNR record exists to separate
them. **Getting more margin on RF33 is therefore a prerequisite for a clean
M4 gate, not an optimization.**

## On not taking the radio

One bounded live capture was permitted this session (RF33, Antenna B, 8 Msps,
up to 120 s, before the 02:10 meteor-window handover). It was not taken, and the
reasoning is recorded because the option has now expired for the night:

* the existing `hit_rf33.cs16` was **already** taken on Antenna B at 8 Msps, so
  a fresh capture buys **length, not SNR** — and SNR is the binding constraint
  identified above;
* a longer record is not the current blocker. The blocker is three unbuilt
  chain rungs and two unextracted table sets;
* 8 s = 32 frames is ample for everything M4 can currently check, and LLS
  repeats fast enough that the SLT is very likely already in this record.

**The right capture to take next is not simply longer — it is higher SNR.**
Antenna B with the bias-T LNA on RF33, or a re-aim, measured against the
existing record by the `m4_data.py` MER number, which is now a calibrated
instrument. 120 s at +3 dB would unblock both the FEC gate and M5/M6.

## Assumptions added in M4

    SC1 A/322 5.2.3's register is GALOIS, not Fibonacci.  SETTLED BY THE
        SPEC'S OWN VECTOR: 1 of 16 conventions reproduces the printed 24 bits.
    SC2 the eight output stages D7..D0.  SETTLED BY THE AIR: 1 survivor of
        16875 against a 32-bit CRC, then confirmed by 8 independent
        measurements, by the Figure 5.6 mirror, and by 3572 dummy cells.
    D1  inside the HTI branch of Table 9.9 the PDF's brace placement leaves it
        ambiguous whether L1D_plp_HTI_cell_interleaver sits inside the
        inter_subframe != 0 arm or after the if/else.  Read as AFTER.
        RISK LOW -- 1 bit, and the all-ones L1D_reserved run arbitrates: a
        wrong reading shifts every later field and the run stops being ones.
    D2  the FI symbol counter includes the preamble as symbol 0.  SETTLED BY
        THE AIR: the dummy-cell correlation is 1.0000 at FI symbol 35 for data
        symbol 34 and at chance for every other index.
    N1  A/322 Annex C 2D-NUC quadrant rule x = [w, -conj(w), conj(w), -w] with
        index = the decimal value of (y0..y_MOD-1), y0 the MSB.  RISK LOW --
        it is printed explicitly in 6.3.4.2 and the mean-power identity closes
        for all 36 position vectors.

---

# M5 results (2026-08-06) — the cell map CLOSED, and the next two rungs' tables

Offline work on captures already banked in `lab/`. **No radio was opened.**
Tools: `lab/extract_pilot_tables.py` -> `spec_pilots.py`, `lab/m5_cellmap.py`,
`lab/extract_bitint_tables.py` -> `spec_bitint.py`, `lab/m5_hti.py`.

## Step 1 — Table D.1.4, and the 190 / 2467 discrepancies go to ZERO

M4 ended by naming its own blocker: *re-extract Table D.1.4, gated on the
constant-data-carrier identity.* That is done, and the diagnosis M4 attached
to it was **partly right and mostly wrong.** All three causes are now named,
and the arithmetic closes exactly:

    190  = 127 + 63  + 0
    2467 = 1609 + 804 + 54

**(a) Table D.1.4 really was mis-read — but it costs 54 of the 2467 and
nothing at all in subframe 0.** The printed table gives each SP pattern a row
group with the label vertically centred in a merged cell. The `SPx_4` groups
are THREE rows tall, so the label prints against the middle row and
`pdftotext` hands the group's first row to the `SPx_2` label above it. Read
that way, `SPx_4` silently loses its first additional continual pilot — and
that lost carrier is exactly the one congruent to `Dx mod (Dx*4)`, i.e. the
one lattice phase that needed it. **16K/SP4_4 is {3460, 5768, 11452}, not
{5768, 11452}.** M4's "roughly 35 additional CPs per 16K symbol" estimate was
off by an order of magnitude; the true pilot deficit is 1 cell in 3 symbols
out of every 4.

**(b) The real defect** was using the SBS symbol's TOTAL data-cell count
(Table 7.5/7.6) where A/322 7.2.6.4 says cell multiplexing uses the ACTIVE
count (Annex F). Their difference *is* `L1D_sbs_null_cells`.

**(c) And the dummy-cell correlation** measures the origin of the SBS
symbol's first TOTAL data cell, while available-cell indexing starts
`N_null_low` cells later — the null cells sit at the two band edges
(7.2.6.4, Figure 7.14). With `N_null` odd the split is 63/64 and 804/805.

### The gate: 2420 closed identities, with a control that must fail

`extract_pilot_tables.py` parses D.1.4, D.1.5, 7.3/7.4, 7.5/7.6 and all ten
Annex F tables from BOTH editions and diffs them: **1717 cells per edition, 0
disagreements.** Two witnesses agreeing only proves the extraction is
faithful, not that the READING is right, so the pilot grid is then rebuilt
from first principles and required to satisfy:

| gate | what it pins | n |
|---|---|---|
| G1 | data carriers CONSTANT over every SP lattice phase | 185 |
| G2 | that constant == Table 7.3/7.4 | 185 |
| G3 | SBS geometry == Table 7.5/7.6 | 185 |
| G4 | Annex F active <= Table 7.5/7.6 total | 925 |
| G5 | NoC == NoCmax - Cred*Cunit | 15 |
| G6 | Annex F monotone non-increasing in SPB | 925 |

**2420 pass, 0 fail.** The naive 2-row grouping is kept *inside the module*
as a built-in control and fails **143** of them — every single `SPx_4`
combination. A gate that cannot fail is not a gate.

### Off the air (`m5_cellmap.py`, reproduced on 3 frames)

* **GATE B.** Table 7.5/7.6 minus Annex F reproduces BOTH signalled
  `L1D_sbs_null_cells` values — `5136-5009 = 127` and `10272-8663 = 1609`.
  Those came out of the LDPC-decoded L1-Detail; nothing in the pilot
  arithmetic can fit them.
* **GATE C.** The dummy-cell correlation HOLDS: **1.0000** on subframe 0
  (3572 signs, runner-up 0.163) and **0.9953-0.9969** on subframe 1 (8505
  signs, runner-up 0.147).
* **GATE D.** The walked-back cell budget lands on **3487** for subframe 0 —
  which is `4851 - 484 - 880` from three independently decoded quantities —
  and on **ZERO** for subframe 1, which inherits nothing from a preamble.
  **Both discrepancies are 0, on every frame.**

Every wrong reading is run as a control, and each reproduces its own defect:

    M4's reading (TOTAL, no null offset, naive D.1.4)  sf0  -190  sf1 -2467
    SBS TOTAL instead of Annex F ACTIVE               sf0  -127  sf1 -1609
    no null-cell offset on the measured origin        sf0   -63  sf1  -804
    null split CEIL at the low edge                   sf0    +1  sf1    +1
    naive D.1.4 grouping only                         sf0    +0  sf1   -54
    THE CORRECTED READING                             sf0    +0  sf1    +0

The first control regenerating M4's -190 / -2467 *to the cell* is the point:
this is a decomposition, not a re-tune.

## Step 2 — the data-PLP bit interleaver tables, and a corrected map

**M4's rung-3 table list was wrong**, and it would have sent this session
hunting in the wrong place. It read "group-wise interleaver (Tables
6.10/6.11), block interleaver (Tables 6.12/6.13)". In fact:

* 6.10/6.11 are the **BLOCK** interleaver configurations (Type A / Type B);
* the group-wise permutations are **not in Section 6 at all** — they are the
  **Annex B** table series, B.1 for `Ninner` = 64800 (Ngroup 180) and B.2 for
  16200 (Ngroup 45);
* **6.12/6.13 are the mandatory modulation/coding combination tables** — a
  conformance checklist of check-marks that carry nothing a receiver needs,
  and whose glyphs do not survive `pdftotext` at all. Recorded as unreadable
  rather than worked around.

`spec_bitint.py` now holds 120 Annex B permutations, Tables 6.8/6.9 (120),
6.10 (12), 6.11 (12) and 6.14 (12). Two witnesses identical. **300 gates
pass, 0 fail** (H1 bijectivity x120, H2 identity headers, H3/H4 the two
block-interleaver arithmetics, **H5 the cross-table `Nr1*Nc == Npart1`**, H6
demux counts, H7 type letters, **H8 the printed rate labels landing in blocks
1..12 in order**).

H5 is the sharpest: 6.10 and 6.11 are printed two pages apart and describe
two *different* interleavers, yet every `Nr1*Nc` reproduces the other table's
`Npart1`. H8 exists because H1 alone is not enough — a permutation stays a
permutation if two rate blocks are swapped, so the stripped code-rate labels'
positions are recorded on the way out and checked, turning what would have
been an assumption into a gate.

## Step 3 — the HTI Twisted Block Interleaver: geometry closed, ST open

RF33 PLP 0 is the easy HTI case: `inter_subframe = 0` (no Convolutional Delay
Line) and `cell_interleaver = 0` (the Cell Interleaver is bypassed), so only
the TBI stands between the cell stream and the bit de-interleaver.

**A closed identity pins the geometry:**

    L1D_plp_size / (Ninner/MOD) == L1D_plp_HTI_num_fec_blocks + 1
             199800 / 2700      == 73 + 1 == 74

Three independently sourced numbers meet — `plp_size` and `num_fec_blocks`
from the LDPC-decoded L1-Detail, 2700 from Table 6.14. So `NTI = 2`,
`N_fec_TI_max = 37`, **37 FEC blocks per TI block** and `N_virtual = 0`;
`37 x 2 x 2700 = 199800` exactly. **This corrects the M4 note's "2 TI blocks
x 73 FEC blocks".**

**The wall, precisely.** A/322 7.1.5.4's three diagonal-read equations print
as `= mod ,` / `= mod ,` / `= + mod ,` — every variable lost to the math
font, in BOTH editions and BOTH extraction modes — and the twisting parameter
is named once and never given a value anywhere the extraction can see.

Rather than guess, `m5_hti.py` enumerates 24 readings x ST in {0,1,2} and lets
structure eliminate:

* **T1** a reading must be a permutation for EVERY legal (Nrows, Ncols), not
  just RF33's. The enumerations that only work because 2700 and 37 are
  coprime die the moment they share a factor — 202-225 failures each.
* **T2** synthesized round trip, including with virtual FEC blocks present
  (which RF33 does not exercise).
* **T3** no two consecutive output cells from the same input FEC Block.
  **This did not do what I expected and is recorded as measured**: it kills
  the row-major linear array, but says nothing about ST, because row-major
  *enumeration* already changes column every step.
* **T4** the prose fixes enumeration order and the linear-array convention.

One reading stands: `row = l // Ncols`, `col = (l mod Ncols + ST*row) mod
Ncols`, `z = (col - Nvirtual)*Nrows + row`.

**OPEN — ASSUMPTION T1: the twisting parameter ST.** T1-T3 do not constrain
it at all; `ST = 0` passes every one of them. The only argument against
`ST = 0` is that it makes the interleaver untwisted and the parameter
pointless — prose, not proof, and **not counted as a gate here**. The search
is now one small integer over 0..36, with LDPC convergence as the
discriminator, exactly how ASSUMPTION L1 was settled in M3.

## Assumptions added in M5

    P3  with an odd L1D_sbs_null_cells the low band edge takes floor(n/2).
        SETTLED BY THE AIR -- the ceil reading misses by exactly +1 on both
        subframes.
    B1  Annex B code-rate block order.  SETTLED BY H8 -- the labels are read
        back positionally in all 10 tables rather than assumed.
    B2  Tables 6.12/6.13 carry no extractable content (check-mark glyphs).
        Conformance checklist, not decoder data; nothing depends on them.
    T1  the A/322 7.1.5.4 twisting parameter ST.  OPEN.  Narrowed from a
        structure hunt to a 37-value integer sweep; only the FEC can settle
        it, and a wrong value will not announce itself.

## Where M5 stops, and the sharpest next step

Rung 1 (cell de-multiplex) is now **closed and gated** — that was the whole
point of this session, because a 190- or 2467-cell offset would have poisoned
the time de-interleaver silently and looked like an FEC failure.

Remaining for PLP 0: the TBI ST sweep, then the bit de-interleaver (tables
now in hand and gated), then LDPC 16200 rate 11/15 (already in
`spec_ldpc.py`, decoder already gated) + BCH + descramble.

**Sharpest next step:** wire rungs 1-4 end to end and sweep ST 0..36 against
LDPC convergence. That single sweep either produces the first payload bits or
proves the chain wrong — and note the standing link-budget caveat: at 18.8 dB
MER PLP 0 sits at the 64QAM 11/15 threshold, so a total failure across all 37
ST values would remain ambiguous between "chain wrong" and "SNR short".
**Getting more margin on RF33 is still a prerequisite for a clean gate, not
an optimization.**

---

# M6 results (2026-08-06) — the payload chain CLOSED: PLP 0 decoded off real air

Offline work on captures already banked in `lab/`. **No radio was opened.**
Tools: `lab/m6_tbi.py`, `m6_cells.py`, `m6_bicm.py`, `m6_payload.py`.

## The headline

    PLP 16   LDPC converged, BCH syndrome 0                      every frame
    PLP  0   **74/74 FEC Blocks converged, 0 unsatisfied checks**
             **74/74 BCH syndromes zero**                        every frame
    -> descrambled Baseband Packets -> ALP -> IPv4/UDP -> ROUTE/DASH

Reproduced on 20 consecutive Frames of `hit_rf33.cs16`. The first ATSC 3.0
payload bits, and then IP datagrams, off real air in this project.

## M5's plan was wrong, and the way it was wrong is the lesson

M5 ended with a single named unknown — the HTI twisting parameter — and a
plan: sweep `ST = 0..36` against LDPC convergence. **That sweep was run. It
returned 37 flat failures, and it could not have contained the answer.**

A/322 7.1.5.4 does not leave the twisting parameter free. It *defines* it:

    R_i     = i mod N_r
    T_i     = R_i mod N_c
    C_i     = (T_i + floor(i / N_r)) mod N_c
    theta_i = N_r * C_i + R_i          (skip while theta_i < N_virtual * N_r)

M5's surviving reading is the **transpose** of this: it advances the column
fastest where the spec advances the row fastest, and it matches the spec for
*no* value of ST. M5's T-GATE 4 chose it because the prose says the read runs
"rightwards along the row" — and preferred that prose to the equations
printed directly beneath it.

**Record it alongside M4's scrambler lesson, because it is the same shape:**

* M4: *"the spec is wrong" is the most expensive conclusion available and
  must be the last one reached.*
* M6: **when prose and equations disagree, the equations win** — and *"the
  PDF ate the equation"* is a statement about the extractor, not about the
  document. The variables M5 called unreadable in both editions and both
  extraction modes are read intact by PyMuPDF.

A gate has been added that M5 could not have written: A/327 Figure 6.5 prints
a **worked 4x3 example with one virtual FEC Block**, expected output
`b g a f d e c h`. `m6_tbi.py` reproduces it 8 of 8, and it exercises the
virtual-cell skip rule — which RF33, with `N_virtual = 0`, never exercises.
**A test vector that tests more than the air does.**

## The second bug: the parity interleaver, applied twice

Fixing 7.1.5.4 alone still fails. A/322 6.1.3.2's Type B encoder has **no
permutation at all** (`lambda_{Nouter+k} = p_k`), but `m3_ldpc.encode()` —
written in M3 against the Type A prose and reused for Type B — applies
6.2.1's parity interleave *inside* the encoder, and `parity_check()` builds
H to match. So m3_ldpc's codeword ordering already **is** the spec's `u`, and
applying 6.2.1 again in the receiver double-counts it.

`verify_typeb_parity()` proves this numerically against a from-scratch
6.1.3.2 encoder rather than arguing it.

**This class of bug is invisible to a round trip, because the transmitter and
the receiver share it.** M6's synthesized round trip passed at every stage
while the air failed at every stage — which is precisely the failure mode M2
and M3 warned about and is worth stating as a law:

> **A round trip proves the code inverts itself. It cannot prove a spec
> reading. Only an external referee can — printed test vectors, an
> independent implementation, or the air.**

## The referee that broke the "three rungs, no umpire" deadlock

M4 and M5 both recorded the same trap: nothing sits between the cell stream
and the LDPC, so three rungs must be right at once and a wrong one is silent.

**PLP 16 is that umpire, and it was in the L1-Detail the whole time.** Its
decoded configuration is QPSK 2/15, `plp_size` 8100 = **exactly one FEC
Block**, `num_ti_blocks` 0 -> NTI 1, so the Twisted Block Interleaver is an
8100x1 array and **no twisting is possible**. PLP 16 therefore tests the cell
de-multiplex, the bit de-interleaver, the demapper, the LDPC, the BCH and the
descrambler with the time interleaver removed as a variable — and QPSK 2/15
decodes ~20 dB below the available MER, so a failure cannot be blamed on the
link.

It decoded first, on the first attempt. That converted a three-unknown
problem into a one-unknown problem, and is the single most useful structural
move in this milestone.

## The link-budget caveat is REFUTED, measured two ways

M4 and M5 both recorded that "getting more margin on RF33 is a prerequisite
for a clean FEC gate, not an optimization", on the grounds that 18.8 dB MER
sits at the 64QAM 11/15 threshold. **That was wrong.**

1. **Our own decoder's AWGN threshold, measured**: 64QAM-NUC 11/15 converges
   4/4 down to **15.0 dB** and 3/4 at 14.5 dB. A/327 Annex B.6 puts the
   published threshold at 17.5 dB. The measured MER on this capture is
   **19.7 dB** — 2 to 5 dB of margin depending on which figure you use, not
   zero.
2. **The capture then decoded 74/74.** A hypothesis that says the data cannot
   decode is refuted by the data decoding.

Independently, the user's HDHomeRun decodes RF33 cleanly on the *same
antenna*, which exonerated the link before either measurement above.
**Recorded as a reversed finding: the fault was ours, twice, and "we need a
better capture" was a comfortable answer that delayed the real one.**

## What the payload actually contains

`m6_payload.py` walks A/322 5.2.2 Baseband Packet headers -> A/330 ALP ->
IPv4/UDP. Over 8 Frames: 592 Baseband Packets -> 612 ALP packets (605 IPv4),
2 resyncs.

    <src> -> <mcast A>  :<pA>    fMP4 (mdia / trex / tref boxes)
    <src> -> <lls>      :<pL>    ROUTE/ESG: S-TSID, USBD,
                                 OMA BCAST SGDU, service1..9
    <src> -> <mcast B>  :<pB>
    <src> -> <mcast C>  :<pC>
    <src> -> <mcast D>  :<pD>
    <src> -> <mcast E>  :<pE>

(Addresses are written generically on purpose. An ATSC 3.0 broadcaster's
multicast addresses and ports are DERIVED FROM ITS VIRTUAL CHANNEL NUMBERS, so
a list of them is a channel lineup in disguise and identifies a metropolitan
market. Nothing in this project hard-codes them either -- the receiver
discovers flows by classifying them, which is the whole point of the transport
gate below.)

The ESG session carries a complete `<S-TSID>` and `<BundleDescriptionROUTE>`
(serviceId 65024) plus SGDU fragments naming `service1` through `service9`,
so **RF33 carries at least nine services**.

**A third bug was found here and it is worth recording**: A/322 5.2.2's
Baseband Packet header needs the full Optional/Extension Field decode
(OFI 00/01/10/11, 13-bit pointer). A two-byte guess left the FEC perfectly
happy — 74/74 still converged — and silently shredded **74% of the payload**
at the ALP layer. *The FEC gate says nothing about the layers above it.*

## Assumptions added in M6

    X1  A/322 6.3.3's demultiplexer (defined only by Figure 6.11, an image)
        is the serial map y_{i,s} = q_{MOD*s+i}.  SETTLED BY THE AIR --
        PLP 16 and PLP 0 both decode with it, and all 720 permutations of
        the six sub-streams were swept: only the identity works.
    X2  group-wise direction Y_j = X_{pi(j)}.  SETTLED BY THE AIR; the
        inverted reading is kept as a control and fails.
    X3  Type A Block Interleaver equations.  GATED, not assumed --
        6.2.3.1's printed worked example (256QAM/64800) is reproduced.
    T1  (was: the twisting parameter ST.)  **CLOSED, and the question was
        malformed.**  7.1.5.4 defines it; gated on A/327 Fig 6.5.
    L1  the LDPC parity interleave index order.  SETTLED BY THE AIR for
        Type A (PLP 16, rate 2/15) and for Type B (PLP 0, rate 11/15).

## Where M6 stops, and the sharpest next step

Reached and gated: the full BICM receive chain for data PLPs, the HTI, ALP
depacketization, IPv4/UDP, and the ROUTE session structure.

**Not reached: the SLT.** No datagram to the well-known LLS address
224.0.23.60:4937 appeared in the Frames decoded so far, so the virtual
channel numbers for RF33 are still unknown and cannot yet be checked against
the six unencrypted services the user's HDHomeRun reports (list in gitignored
`lab/ground_truth.txt`).

**Sharpest next step, in order:**

1. **Find the LLS.** Decode more Frames (the capture holds ~32) and, if it
   still does not appear, handle ALP segmentation/concatenation
   (`payload_configuration = 1`) and IP fragmentation, either of which would
   drop an SLT while leaving 99% of the media flow intact.
2. **Reassemble ROUTE/LCT objects by TOI** — the ESG SGDU fragments already
   in hand name nine services and are gzip XML; they may carry the service
   names without the SLT.
3. **fMP4 -> playable file** from the video flow once objects reassemble.
4. **Subframe 1 / PLP 1** (256QAM-NUC 11/15, 64800 LDPC, 3 TI Blocks x 39
   columns) — and note A/327 Section 6 specifies a **frequency de-interleaver
   RESET rule for the second Subframe** which this project has not
   implemented; missing it would break PLP 1 specifically.

---

# M7 results (2026-08-06) — WATCHABLE VIDEO. 114 seconds of 720p HEVC off real air.

Tools: `lab/m7_capture.py` (the one radio session), `lab/m7_route.py`
(Frames -> IP datagrams), `lab/m7_objects.py` (datagrams -> objects),
`lab/m7_play.py` (objects -> a file a player opens).

## The headline

    ffprobe   hevc (Main 10) 1280x720, 60000/1001 fps, 114.152 s
    ffmpeg    6840 frames decoded to a null sink, ZERO error lines
    still     rf33_frame.png -- a clean, full-colour, artefact-free picture

Plus, from the same capture:

    LLS / SLT     FOUND.  bsid 540, 9 services + the ESG service.
    ROUTE objects 38 of 38 COMPLETE (every FDT, S-TSID, MPD and SGDU)
    MMT MPUs      232 of 239 COMPLETE; the 7 that are not are the first and
                  last MPU of each flow, i.e. the two ends of the recording
    audio         AC-4, 6 channels, 48 kHz, 116.1 s, MPUs COMPLETE, muxed
                  with the video into an MPEG-TS

**"First decoded frame" and "watchable video" are different achievements.
This is the second one.** It is not, however, a *player*: nothing here runs
in real time, and the audio codec is one this ffmpeg build can demux but not
decode.

## The single biggest correction: the media flow is MMTP, not ROUTE

M6 recorded all six flows as "ROUTE/DASH". Five of them are. The sixth --
the one carrying the video -- is **MMTP (ISO/IEC 23008-1)**, and the SLT then
confirmed it in writing: that service's `BroadcastSvcSignaling` carries
`slsProtocol="2"` (MMTP) while every other service carries `slsProtocol="1"`
(ROUTE). A/331 allows both and this multiplex uses both **at the same time.**

It was provable before the SLT arrived, which matters more than the answer:
parse the flow as MMTP and the 32-bit `packet_sequence_number` steps by
EXACTLY ONE across 30,504 of 30,506 consecutive datagrams. Nothing else puts
a monotonic +1 counter at a fixed offset. The `ftyp` brand agrees -- `mpuf`
is the ISO 23008-1 MPU brand, not a DASH brand.

**Lesson to keep: "which transport is this?" is a question with a cheap,
external, statistical answer, and M6 answered it by assumption instead.**

## Three bugs, and all three were silent

### 1. The `header_mode == 1` swallow — cost 62% of the payload

M6's ALP walk accepted any length that fit in the buffer. A mis-read
`header_mode` bit produced a 19-bit length of a quarter of a megabyte, the
walk consumed it as one "valid" packet, and 62% of the media flow vanished
while the run reported *one resync*. The FEC was 74/74 the whole time.

The fix is not a better length reading. It is **a referee from another
layer**, and there were two available for free:

* **The anchor contradiction.** Every Baseband Packet header points at the
  first ALP packet that STARTS inside it. An anchor strictly inside a
  candidate ALP packet is a *proof* the candidate is wrong -- an ALP packet
  cannot begin inside another one. One anchor per ~1464 bytes bounds the
  damage of any mis-read to a single Baseband Packet.
* **The IPv4 gate.** `packet_type` 0 payloads must begin `0x4_` and carry a
  `total_length` that fits.

With those two, 268 IP datagrams per 8 Frames became 743, and the LLS
appeared in the *same* 8 Frames that had shown none.

### 2. The link-layer signalling sub-header — cost one MFU per second

A/330 packet_type 4's `length` counts the **signalling payload only**; the
5-byte signalling sub-header (`signalling_type`, `..._type_extension`,
`..._version`, `..._format`) sits between the base header and it. Reading it
as payload leaves the walk 5 bytes short, exactly once per Link Mapping
Table. The LMT repeats about once a second.

**At 8 Frames this costs nothing. At 476 Frames it punched a 1400-byte hole
in essentially every single MPU** -- 239 of 239 came back PARTIAL, each
missing exactly one packet. Proven at the byte level: the next ALP header
lands at +5, not +0, and the bytes at +5 are `45 00 02 f6` (an IPv4 header
whose total_length is the datagram we were about to lose).

Fixed: 118 resyncs over 476 Frames became **1**.

### 3. The MFU Data Unit layout — cost a "codec error" that was an offset error

    MFU Data Unit = [14-byte DU header][34-byte MMT hint sample][media sample]

and the DU header rides on **every** fragment, with `offset` giving the byte
position of that fragment inside the MEDIA sample (hint bytes not counted).
The mdat is laid out as ALL media samples followed by ALL hint samples --
which is what the second `traf`'s `data_offset` says -- so the DUs cannot
simply be concatenated in arrival order.

Getting the DU header length wrong by four bytes (18 instead of 14) left
every check *nearly* right and produced

    [hevc] Invalid NAL unit size (1000870736 > 590)

**An offset error wearing a codec error's clothes.** Worth naming as a class:
when a demuxer reports impossible sizes, suspect the byte offset before the
bitstream.

The reading is now gated three ways, all external:

    offset(k+1) - offset(k) == len(DU_k) - 14, over hundreds of fragments
    last fragment's offset + payload lands exactly on the trun sample size
    reassembled media == sum(trun sample sizes)          533387 == 533387
    media + hint + 8  == the mdat box's own declared size 537475 == 537475

### 4 (not a bug, but it looks like one): every MPU starts at time zero

An MPU is an independent ISOBMFF file: `mfhd` sequence_number 1, `tfdt`
baseMediaDecodeTime 0. Concatenate 57 of them unchanged and ffprobe reports
`duration 2.002` -- the player shows two seconds and stops. `m7_play.retime()`
renumbers the fragment and accumulates the decode time from the trun's own
sample durations. 2.002 s -> **114.152 s**.

## The LLS was found, and BOTH of M6's named causes were wrong

M6 named two candidates for the missing SLT, in order: ALP
`payload_configuration = 1` (segmentation/concatenation), and IP
fragmentation. **Neither was the cause.** Over 476 Frames the whole capture
contains **one** `payload_configuration = 1` packet and **zero** IP
fragments. The SLT was being destroyed by bug #1, one flow at a time, and
"decode more Frames" would not have fixed it either -- the same 8 Frames that
had shown no LLS produced three LLS datagrams the moment the walk was gated.

**Record it as a reversed finding.** Two plausible spec-shaped hypotheses, a
sensible ordering, and the actual cause was a missing gate in our own walker.
It is the M4/M6 lesson again: the expensive explanation was reached for
before the cheap one had been ruled out.

The SLT (`lab/m7_out/lls_1_SLT.xml`, gitignored) gives bsid 540 and ten
entries. Against the six unencrypted services the user's HDHomeRun reports
(the list is location-identifying and lives in gitignored
`lab/ground_truth.txt`; the HDHomeRun adds 100 to the major channel for 3.0
services), **RF33 carries exactly two of them**:

    serviceId 2                       slsProtocol 2 (MMTP)   NOT protected
                  the flow decoded above
    serviceId 1                       slsProtocol 1 (ROUTE)  NOT protected
    serviceId 3,4,5                   protected="true", DRM system ID present
    serviceId 6,7                     broadbandAccessRequired="true"
    serviceId 8,9 (M6's "promising")  serviceCategory 3 -- app-based, not
                                      linear video; 8 is broadband-required
    serviceId 65024                   serviceCategory 4 -- the ESG

So the guess that serviceIds 8 and 9 were the watchable unencrypted pair was
wrong in an instructive way: *absence of encryption markers in an MPD is not
presence of a linear video service.* The service that turned out to be both
unencrypted and watchable is the one whose signalling is not DASH at all.

The remaining four HDHomeRun entries are not on RF33 and must live on other
RF channels.

## What the longer capture bought, precisely

    120.000 s captured, capture-integrity gate PASS
      960,000,000 samples vs wall*fs 959,896,940  (ratio 1.000107)
      0 overflow reads, control readback verified
    476 Frames decoded, 35,223 of 35,224 FEC Blocks converged and BCH-zero
      (99.997%; the one failure is a single Block in one Frame)

It bought completeness and nothing else, which is what was predicted:

    8 Frames  (2 s):  0 complete video MPUs -- both are cut by the window
    476 Frames (118 s): 57 of 60 complete; the 3 that are not are the first
                        and the last, i.e. the two ends of the recording

An MPU is 2 s of video, so **any capture shorter than ~3 s cannot contain a
complete one, and every extra second past that is one more.** The SLT repeats
roughly every second and the ESG carousel is on a ~5 s cycle, so the
signalling never needed 120 s -- it needed the walker fixed. The video did
need the length.

## Where M7 stops, and the sharpest next step

Reached and gated: ALP with its signalling sub-header, IPv4 (+ fragment
reassembly, unused here), LLS/SLT, ROUTE/LCT object reassembly with FDT
completeness, MMTP/MPU reassembly with three independent size referees, and
an ISOBMFF writer that puts the MPUs on one timeline. 38/38 ROUTE objects and
232/239 MPUs complete.

**Sharpest next step, in order:**

1. **PLP 1 / subframe 1.** The other unencrypted service on this multiplex
   (serviceId 1, ROUTE) has its SLS at an address that never appeared in any
   Frame, so its media is almost certainly on the PLP this project does not
   decode. A/327 Section 6's **frequency de-interleaver RESET rule for the
   second Subframe** is still unimplemented and would break PLP 1
   specifically.
2. **Real time.** 476 Frames took 222 s of wall clock on 14 cores -- about
   1.9x slower than real time for one PLP. The LDPC decoder is the whole
   cost and it is pure NumPy.
3. **AC-4.** The audio objects are complete and valid; this ffmpeg build has
   no AC-4 decoder. A Dolby-capable player, or an AC-4 decoder, would close
   the last gap between "video file" and "television".
4. **Don't trust a gate that cannot fail.** `mdat_got == mdat_declared`
   became vacuous the moment the reassembler started zero-filling missing
   samples to their declared size. `media_got == media_want` is the gate that
   still bites. Audit the others for the same defect.

---

# M8 results (2026-08-06) — the other multiplexes: L1 decoded, and the pipeline's reach measured honestly

Tools: `lab/m7_capture.py` (three radio sessions, reused unchanged),
`lab/m8_l1.py` (the new L1 path).

## The headline

    RF25   L1-Basic + L1-Detail DECODED.  LDPC 0 unsatisfied, BCH zero,
           L1D_crc PASS.  bsid 1408.
    RF30   L1-Basic + L1-Detail DECODED.  LDPC 0 unsatisfied, BCH zero,
           L1D_crc PASS.
    RF29   NO SIGNAL.  No bootstrap anywhere in 120 s; 25 dB below RF30.
    ----
    And the SLT on RF25/RF30 is NOT reached, because those multiplexes are
    not the same machine as RF33.

**The task was framed as "an application of a working pipeline, not a
rebuild". It is not, and finding that out is this milestone's result.**
RF33 is non-LDM, HTI-interleaved, 16K LDPC, single Preamble symbol. RF25 and
RF30 are **LDM** (two layers superimposed on the same cells), **CTI**
(convolutional time interleaving, not the twisted block interleaver M5/M6
closed), **64K LDPC**, and **two Preamble symbols**. Three of those four are
rungs this project has never built.

## Two things RF33 could not have taught, both load-bearing

### 1. NP = 2 — the Preamble is not always one symbol

RF33 signals `L1B_preamble_num_symbols = 0`, so A/322 7.2.5.2's L1-Detail
Block Interleaver has `Lc = NP = 1` column and is the identity;
`m4_l1detail.py` refuses anything else outright. RF25 and RF30 use **L1-Basic
Mode 1**, which costs **3820** of the first Preamble symbol's 4851 cells and
leaves nowhere near enough room for L1-Detail, so both use NP = 2.

7.2.5.1 then matters: the FIRST Preamble symbol runs at the minimum NoC for
its FFT size, and the remaining ones at the NoC signalled by
`L1B_preamble_reduced_carriers` — a different symbol width inside one
Preamble. On RF30 that is 4851 cells, then 5136.

### 2. Parity REPETITION — and it had been failing in plain sight

A/322 6.5.2.7 applies parity repetition to **L1-Basic Mode 1 and L1-Detail
Mode 1 only**. RF33's L1-Detail is Mode 3, so `m4_l1detail.geometry()` has no
`Nrepeat` term at all. Both new multiplexes are Mode 1.

    Nrepeat = 2*floor(C*Nouter) + D        C = 61/16, D = -508   (Table 6.23)
    cells   = (Nfec + Nrepeat) / eta

**This did not have to be guessed. The transmitter signals the answer.**
`L1B_L1_Detail_total_cells` is a field the receiver can independently derive,
and without the repetition term the two numbers disagree:

    RF30   computed 1944   signalled 3611      RF25   computed 2232   signalled 4387
    with 6.5.2.7:   3611 == 3611                       4387 == 4387

M4 already printed that disagreement, as the word "MISMATCH", for two
different channels, and it was read as "this channel is odd" rather than as
"our geometry is missing a term".

> **A signalled field the receiver can also derive is a free gate. This one
> had been failing silently, in output we had already read.**

The LDPC agreed only second, and it nearly lied on the way: under M4's
geometry the decoder settles at **98 unsatisfied checks of 12960** — close
enough to look like a marginal link, not close enough to decode. Reading
"nearly converged" as noise would have sent this back to the antenna. It is
the M4/M6 lesson in a third costume: *the expensive explanation (weak signal)
was available before the cheap one (a missing term) had been ruled out.*

## The gates, and one of them is honestly vacuous

1. **cells computed == cells signalled**, exactly, on two independent
   channels with different Ksig. Above.
2. **LDPC 0 unsatisfied + BCH syndrome zero + L1D_crc**, on every Frame tried.
3. **The M4 control fails.** M4's own geometry on the same cells: RF30
   127 unsatisfied / BCH fail, RF25 4308 / BCH fail.
4. **The interleaver column height is swept, not assumed.** `--sweep-lr`
   scores `Lr` over a plus/minus 6 window. On RF30: **144 unsatisfied at the
   spec's floor(total/Lc) = 1805, and 3637-4289 at every other value.** That
   settles 7.2.5.2's "the first Lr x Lc cells are interleaved, the remainder
   passes through" reading against 12 controls.
   *On RF25 the same sweep is only marginal (3771 vs 3897-4783) — reported
   because it is what happened, not because it helps.*
5. **The QPSK cliff.** L1-Detail ends part-way through the last Preamble
   symbol and data-PLP cells of much higher order follow. On RF30 the
   measured MER cliff — 5.5-6.4 dB before, 3.8-4.5 dB after — falls in the
   200-cell window containing the predicted cell 2580. The prediction comes
   from the mapping; the measurement does not.
6. **The transmitter's own clock.** `L1D_time_sec` decodes to 14:38:29 UTC on
   RF25 and 14:42:33 UTC on RF30; those captures began at 14:38:27 and
   14:42:32. Two transmitters, two seconds, an external referee for the parse
   alignment that owes nothing to our code.

**Gate 4 is vacuous under the CORRECT reading and is deliberately scored
under the wrong one.** With repetition the transmitted word is 7222 bits
carrying 336 information bits — rate 0.047 — and the LDPC shrugs off a
permutation of half the cells: 10 of 13 wrong `Lr` values still decode
cleanly. A gate that cannot fail is not evidence, so it is run where it can
still fail. Same defect class as M7's `mdat_got == mdat_declared`.

## What the other multiplexes actually are

    RF25   1 Subframe, 189 OFDM symbols, 8K FFT, GI 512, Cred 3
           PLP 0  Core Layer      QPSK       9/15   64K LDPC   CTI   LLS
           PLP 1  Enhanced Layer  256QAM-NUC 7/15   64K LDPC   CTI   LLS
                  LDM injection level 5.0 dB, same cells as PLP 0
           PLP 16 Core Layer      QPSK       2/15   16K LDPC   CTI
           bsid 1408, frequency interleaver ENABLED

    RF30   1 Subframe, 174 OFDM symbols, 8K FFT, GI 1536, Cred 0
           PLP 0  Core Layer      QPSK       6/15   64K LDPC   CTI   LLS
           PLP 1  Enhanced Layer  64QAM-NUC  6/15   64K LDPC   CTI
                  LDM injection level 4.0 dB, same cells as PLP 0
           frequency interleaver BYPASSED for this Subframe
           L1D_version 0 -- see the caveat below

Both carry `L1D_plp_lls_flag = 1` on the **core** layer, which is the good
news: the SLT does not need the LDM enhanced layer cancelled. The core layer
is demapped with the enhanced layer treated as interference, which is exactly
what LDM is designed for, and QPSK 6/15 has metres of margin for it.

### A caveat recorded rather than buried

RF30 reports **`L1D_version = 0`** — the original A/322 L1-Detail structure.
`spec_l1syntax.py` implements version 1 (2024) with an added block for
version 2 (2026), so RF30's field values *downstream of the version field*
are **PROVISIONAL**. Two facts bound the doubt: the bits themselves are
certain (LDPC, BCH and L1D_crc all pass), and `L1D_time_sec` lands inside the
capture window, so the parse is right at least that far. But its trailing
`L1D_bsid` decodes to 65019 against RF25's plausible 1408 and RF33's 540
(which the SLT independently confirmed), and its reserved run is **0 bits**
where a correct parse leaves a run of ones — RF25 leaves 5, RF33 leaves 7.
**Two of the three parse-health indicators are bad, so RF30's bsid is not
reported as a result.**

## Where M8 stops, precisely

Reached and gated: multi-symbol Preambles, A/322 7.2.5.2's L1-Detail Block
Interleaver, 6.5.2.7 repetition, and therefore the complete PLP configuration
of two multiplexes that previously produced nothing at all.

**Not reached: the SLT, and therefore not a single virtual channel number.**
The wall is four named pieces, in dependency order:

1. **64K LDPC.** `L1D_plp_fec_type = 1` on every LLS-bearing PLP here.
   `spec_ldpc.py` extracts the parity-check address tables for
   `Ninner = 16200` only — its own header says so — and `m3_ldpc._params()`
   reads `LDPC_16200[rate]["rows"]` whatever `n` is asked for. A/322 Annex
   A.1 (13 tables) is in the extracted spec text and has not been mined.
2. **CTI — Convolutional Time Interleaving (A/322 7.1.4).** Every PLP on both
   channels uses it. M5/M6 closed the *Hybrid* time interleaver; this is a
   different machine. `L1D_plp_CTI_start_row` and
   `L1D_plp_CTI_fec_block_start` are decoded and hand over the alignment, so
   the geometry is not a search.
3. **The single-Subframe cell pool.** `m6_cells.py` hard-codes RF33's grid —
   FFT, GI, pilot pattern, symbol count, two Subframes and the PLP split are
   all module constants. It needs to be driven from L1 instead. Both new
   channels are simpler (one Subframe), and RF30 bypasses the frequency
   interleaver entirely.
4. **LDM core-layer demapping.** The smallest of the four: demap the core
   QPSK with the enhanced layer as interference. No cancellation needed for
   the SLT.

**RF29 is not on that list, and not a decoder problem.** It was present on
2026-08-05 at ~01:18 UTC (12 bootstraps, peak ratio 22-28 — weak even then)
and is absent on 2026-08-06 at 14:42 UTC: no bootstrap at any of 13 offsets
across the whole 120 s, and 25 dB less power than RF30 on the same antenna in
the same session. The correct next action for RF29 is another capture at
night, not another line of code.

## Capture-integrity, all three

    RF25   960,000,000 samples vs wall*fs 959,800,416   ratio 1.0002   0 overflows
    RF29   960,000,000                    959,942,623         1.0001   0 overflows
    RF30   960,000,000                    959,818,550         1.0002   0 overflows

Antenna B, 8 Msps, 120 s each, radio_lock owner "atsc3" priority 50, released
in `finally`, heartbeat on a 5 s timer inside the read loop. Outside the
02:10-05:40 meteor window. All three PASS.

## Addendum — the VHF-hi candidates, and the one that is a LINK problem

The 2026-08-05 band sweep flagged seven ATSC 3.0 candidates, and M8 above only
addresses the four UHF ones. The other three are VHF-hi (RF7, RF8, RF13) and
were swept on **Antenna B, which is the UHF yagi and is VHF-deaf.** Antenna C
is the roof wideband vertical antenna, which is the reverse. So the sweep's verdict on those
three was never a statement about the air.

Re-probed on Antenna C, 20 s each, integrity gate PASS on all three
(ratios 1.0014, 1.0014, 1.0015; 0 overflows):

    RF7    no bootstrap on either antenna            -> not an ATSC 3.0 signal
    RF13   no bootstrap on either antenna            -> not an ATSC 3.0 signal
    RF8    bootstrap on BOTH, same fields, and the
           wideband vertical antenna is worth 3.7x in peak ratio       -> REAL, and under-received

RF8's bootstrap decodes to `preamble_structure = 15` — 8K FFT, GI 768,
preamble DX 4, L1-Basic Mode 1 — identically at four offsets and on two
different antennas, which is not something noise does. It is a genuine ATSC
3.0 multiplex.

**It does not decode, and the reason is measured rather than guessed:**

    RF8  (Antenna C)   bootstrap peak ratio 50   pilot coherence 0.72
                       L1-Basic region MER 0.7-1.6 dB   LDPC 46 unsatisfied
    RF30 (Antenna B)   bootstrap peak ratio 90   pilot coherence 0.85
                       L1-Basic region MER 5.0-5.2 dB   LDPC 0, BCH zero

**RF8 is about 4 dB short at the same rung RF30 clears.** That is a link
budget, and no amount of code closes it. The wideband vertical antenna already bought most of
the way there (peak ratio 13 -> 50 versus the yagi); the remaining gap wants a
VHF-hi directional antenna, or a night pass, and not another spec reading.

Recorded because the failure mode matters: **RF7 and RF13 stop at "no
signal", RF8 stops at "not enough signal", and RF25/RF30 stop at "not enough
receiver".** Only the last of those three is ours to fix at a keyboard.

## Sharpest next step, in order

1. **Mine A/322 Annex A.1** (the `Ninner = 64800` LDPC tables) the way
   `extract_bitint_tables.py` mined D.1.4, and gate the result the way M5 did:
   a structural bijection check plus a printed worked example, before the air.
2. **CTI**, gated on `L1D_plp_CTI_start_row` — a field that is decoded, so a
   wrong reading contradicts a known number rather than merely failing.
3. **Drive `m6_cells.py` from L1** instead of from module constants. This is
   the piece that stops the next multiplex from being a rebuild too.
4. Then RF25's PLP 0 -> ALP -> the SLT, and the M7 chain is unchanged
   downstream of that.
5. **Re-capture RF29 at night**, when it was last seen.

---

# M9 results (2026-08-06) — REAL TIME. 0.52x became 2.99x, and M7's diagnosis of why was wrong.

Offline work on captures already banked in `lab/`. **No radio was opened.**
Tools: `lab/m9_profile.py` (the measurement), `lab/m9_accel.py` (CPU),
`lab/m9_gpu.py` (CUDA), `lab/m9_fast.py` (the frame decoder),
`lab/m9_decode.py` (the driver). This is the throughput track and it ran in
parallel with M10's 64K-LDPC/CTI track; **it edits no shared module** —
everything is installed by reversible monkeypatch, and `--accel none` runs the
untouched chain.

## The headline

    476 Frames = 117.625 s of air, on `long_rf33.cs16`

    BEFORE  reference chain, 14 processes    227.206 s wall   0.518x real time
    AFTER   1 process, 16 threads, CUDA       23.401 s / 96 Frames  1.014x
    AFTER   8 proc x 6 thr, NO CUDA AT ALL    23.570 s / 96 Frames  1.006x
    AFTER   8 processes x 6 threads, CUDA     39.375 s wall   **2.987x**

    output: 46,700,035-byte .dg, sha256 a1b3507c181ebf1674dbe7d584779e9b...
            IDENTICAL from both paths, 35,224/35,224 FEC Blocks, 0 resyncs,
            41,426 IP datagrams, same 7 flows

**5.77x end to end, and the decoded bytes are the same bytes.** The gate is
not a checksum of a summary: it is the SHA-256 of the whole 46.7 MB datagram
file that `m7_objects.py` consumes. With the attenuation ladder pinned to the
fixed 0.75 M7 shipped, the fast path also reproduces **M7's own `m7_long.dg`**
byte for byte (46,698,583 bytes, sha256 7013bf58...) — so it matches the
reference chain both as it is today and as it was.

Single-stream matters separately from aggregate: **one process now runs
1.014x**, which is what a live pipeline needs, because process-parallelism
buys throughput at the cost of latency.

## M7 said "the LDPC decoder is the whole cost". It was 17%.

M7's note recorded: *"476 Frames took 222 s of wall clock on 14 cores — about
1.9x slower than real time for one PLP. The LDPC decoder is the whole cost and
it is pure NumPy."* The first sentence is a measurement and it reproduces
(227.2 s here). The second is an impression, and it is wrong by a factor of
six. `m9_profile.py` wraps the functions the decode path actually calls and
reports seconds against one reference — the wall-clock duration of the RF the
Frames represent:

    reference chain, 1 process, 12 Frames             ms/Frame    % wall
    ------------------------------------------------------------------------
    freq de-interleave  (m3_freqint)                    1318.5     23.9%
    LDPC min-sum, inclusive                              914.9     16.6%
      ...of which m3_ldpc.pack(), rebuilt per Block      699.1     12.6%
    descramble          (m4_scrambler, per FEC Block)    748.7     13.5%
    BCH syndrome        (Python big-integer reduction)   630.5     11.4%
    capture load        (resample + de-rotate + notch)   556.2     10.1%
    CPE correction                                       334.4      6.1%
    demod_data          (FFT + channel est + equalise)   275.3      5.0%
    demap (max-log LLR)                                  252.4      4.6%
    dummy_check + constellation_regions                  373.0      6.7%
    fine_timing (41 x demod_data)                        144.8      2.6%
    PlpChain build (H + pack)                             95.4      1.7%
    HTI time de-interleave                                 6.5      0.1%
    ------------------------------------------------------------------------
    TOTAL                                               5460       0.045x RT

**Three quarters of everything attributed to "the LDPC" was `pack()`** — a
Python loop that builds the constant (4320 x 16) padded check-node index array
and that `min_sum_decode` called ONCE PER FEC BLOCK, 74 times a Frame, for an
array that depends only on the code. The belief-propagation arithmetic itself,
the part a GPU can help with, was **3.9% of the wall.**

And the top item was not the FEC at all. It was the A/322 7.3 frequency
de-interleaver: a bit-serial LFSR walk over 8192 addresses, in Python,
recomputed for every OFDM symbol of every Frame, when there are 36 distinct
sequences in a Frame and they never change.

> **Law, and it is the M6/M7 lesson pointed at ourselves: a profile is a
> measurement and "where the time goes" is not something a milestone note gets
> to assert. M7's sentence was written after watching the LDPC run slowly, not
> after timing anything.**

## Six of the seven levers were memoisation, and memoisation cannot change a bit

`m9_accel.install()` swaps eight module attributes and `uninstall()` puts them
back. Seven are pure-function caches, so bit-identity is a property rather
than a hope:

    1  m3_freqint.interleaving_sequence   36 distinct values per Frame
    2  m4_scrambler.sequence              one fixed sequence; grow and slice
    3  m3_ldpc.pack                       constant per code
    4  m3_ldpc.parity_check               constant per code
    5  m6_bicm.bicm_map                   constant per (Ninner, mod, rate)
    6  m3_spec.pilot_values / reference_sequence / common_cps,
       spec_pilots.pilot_carriers         constant per symbol index
    7  m6_cells.demod_data's GEOMETRY     pilot index sets, FFT-bin maps and
                                          ref[pk], cached per (l, sbs); every
                                          arithmetic operation left untouched

The eighth is algebraic and therefore gated rather than assumed:

    8  m6_bicm.bch_syndrome  The reference walks 11880 message bits into one
       Python big integer and reduces it a bit at a time.  The syndrome is a
       LINEAR map over GF(2), so it is exactly a 168 x 11880 binary matrix
       product; packed into uint64 lanes with np.bitwise_count that is ~31k
       word operations instead of ~2M.  The matrix columns are the remainders
       x^k mod g(x) built from the SAME generator polynomial object the
       reference uses.

    m9_accel.verify(), all exact, all required to be array_equal:
      PASS  FI sequence identical, 36 symbols x 4851 addresses
      PASS  scrambling sequence identical, 5 lengths up to 211472
      PASS  ldpc.pack identical
      PASS  BCH syndrome identical, 26 messages (11880 and 2160 bits)
      PASS  batched BCH == scalar BCH, 16 messages
      PASS  pilot tables identical, 72 symbol/SBS combinations
      PASS  demod_data identical ON REAL AIR, 72 symbol/SBS combinations

These alone took the single process from 0.045x to 0.141x — **3.1x for code
that computes exactly the same numbers.**

## Threads, not processes, and an FIR that parallelises EXACTLY

NumPy and SciPy release the GIL inside their kernels, so a single process can
use many cores without the latency cost of chunking. Three of the remaining
costs are embarrassingly parallel within one Frame: the 36 OFDM symbols, the
41 fine-timing candidates, and the 74 demappings.

The interesting one is the 257-tap notch FIR over 7.2 M complex samples
(1129 ms serial). `lfilter(b, 1, x)` run on a block preceded by `len(b)-1`
samples of real history produces **bit-identical** outputs for that block —
the transposed-form accumulator sees the same summands in the same order — so
overlap-save with a thread pool is exact, not approximate:

    lfilter serial       1128.7 ms
    lfilter 16 threads    115.6 ms      np.array_equal -> True
    np.exp serial         341.1 ms
    np.exp 16 threads      46.5 ms      np.array_equal -> True

**A floating-point trap, caught by the gate in the first run.** The reference
de-rotation is `y *= np.exp(-2j*np.pi*cfo*np.arange(len(y))/fs_post)`, which
Python groups as `((((-2j*pi)*cfo) * n) / fs_post)`. Folding the division into
the scalar — the obvious thing to do when you are about to call this a million
times — is a *different rounding*, and the loader gate reported
`max |delta| = 1.9e-13` immediately. Worth naming as a class:

> **Algebraically equal is not numerically equal. Re-associating a scalar
> chain is an optimisation with a numerical cost, and the only reason it was
> not shipped silently is that the gate compared with `array_equal` and not
> `allclose`.**

## The GPU LDPC: batch the Frame, not the codeword

One Frame of PLP 0 carries **74 FEC Blocks that are completely independent**,
so the natural CUDA unit is a Frame: a (74, 69120) message array, five gathers
and one accumulation per iteration. `m9_gpu.GpuMinSum` is normalized min-sum
in torch, and it is **bit-identical to the NumPy decoder**, which took three
specific design decisions:

1. **No scatter-add.** Every operation in min-sum is exactly-rounded
   elementwise IEEE-754 or an exact integer reduction *except one*: the
   variable-node sum. `np.add.at` visits edges in raveled (check-major,
   slot-minor) order and accumulates strictly left to right; a CUDA
   scatter-add with atomics would not reproduce that, because floating-point
   addition is not associative. So the decoder builds a variable-node
   incidence table whose columns are each variable's edges **sorted by that
   same raveled edge index**, pads with an edge whose message is permanently
   +0.0, and adds the 12 columns in a Python loop — left to right. Same
   summands, same order, same double.
2. **The first minimum, computed deterministically.** NumPy's `argmin` returns
   the first minimal position; `torch.min`'s index tie-break is not contractual
   on CUDA. The position is computed as `argmax((mag == m1) * (dmax - arange))`,
   whose maximum is unique by construction.
3. **The second minimum with multiplicity.** NumPy takes
   `np.partition(mag,1)[:,1]`; here it is the minimum after masking out only
   the first minimal position, which is the same value even when the two
   smallest are equal.

    m9_gpu.verify(), 74 blocks per case, SAME LLRs into both decoders:
      PASS  64QAM 11/15 @ 15.5 dB (near threshold, iters 8..50, 73/74 conv)
            0 bit differences of 1,198,800   CPU 1728.5 ms  GPU 264.6 ms   6.5x
      PASS  64QAM 11/15 @ 19.0 dB (air-like MER, iters 2..5)
            0 bit differences                CPU  408.5 ms  GPU  14.6 ms  28.0x
      PASS  QPSK 2/15 @ 1.0 dB  (the PLP 16 shape, iters 4..5)
            0 bit differences                CPU  756.8 ms  GPU  12.1 ms  62.8x

Every field is compared, not just the bits: the converged flag, the iteration
count and the number of unsatisfied checks. The near-threshold case matters
most — it is the one that runs 50 iterations, and the one block that fails on
the CPU fails on the GPU, identically.

**Reduced precision, measured rather than assumed.** The expectation going in
was that fp16/bf16 would roughly double throughput:

    float64   14.6 ms   BIT-IDENTICAL
    float32   14.7 ms   BIT-IDENTICAL on all three cases (still not a proof)
    bfloat16  13.0 ms   DIVERGED -- one block converged at a different
                        ITERATION COUNT.  Rejected.

**The expectation was wrong in both directions.** fp32 is not faster, because
at 74 blocks the kernel is latency- and bandwidth-bound and nowhere near the
4090's fp64 throughput limit; and bf16, which *is* marginally faster, changes
the answer. The gate caught the bf16 divergence as an iteration-count
difference with zero bit differences — i.e. it would have been invisible to a
bits-only comparison. Comparing every field is what made it visible.

## The one thing that is NOT bit-identical, and the gate that replaces it

`m9_gpu.GpuBicm` moves CPE, the HTI gather, the max-log demapper and the bit
de-interleaver onto the GPU so the LLRs never cross the bus. That took the
per-Frame decode from 545 ms to 159 ms — but it **cannot** be bit-identical,
and the reason is worth recording precisely:

> Both the CPE hard-decision and the demapper are built on `abs(complex)`.
> NumPy's `npy_cabs` and CUDA's complex magnitude are different
> implementations of `hypot`, and they disagree in the last unit in the last
> place on about 40% of inputs. Measured: `torch.abs` disagrees with `np.abs`
> **on the CPU too**, so this is not a GPU artefact and there is no flag that
> fixes it.

So the strict gate is unavailable here and the honest response is to change
the gate, not the claim. The gate that IS available is the one that matters:
**the decoded Baseband Packets must be byte-identical to the reference chain**,
and over the full capture that is 35,224 FEC Blocks and a 46.7 MB datagram
file with a matching SHA-256. The bit-identical path (`--accel gpu`, CPU front
end, GPU LDPC only) stays selectable and is the *default* for that reason;
`--accel gpu-full` is the fast one and says what it costs.

    --accel none      the untouched M6/M7 chain   (oracle AND fallback)
    --accel cpu       memoisation + threads, no CUDA
    --accel gpu       + batched CUDA LDPC, bit-identical everywhere
    --accel gpu-full  + CUDA CPE/demap, gated on byte-identical output

Per the standing GPU-optional law, `--accel none` and `--accel cpu` need no
CUDA at all, and `M9_NO_TORCH=1` keeps torch out of the process entirely.
**That is not a token fallback: measured with `M9_NO_TORCH=1` set, so torch is
never imported, `--accel cpu` on 8 processes x 6 threads runs 1.006x real
time**, byte-identical, on the CPU alone. The GPU is what turns real time into
3x real time; it is not what makes real time possible.

## Where the time goes now

    1 process, 16 threads, 96 Frames        ms/Frame   % wall
    ---------------------------------------------------------
    capture load                              135.7    55.8%
      ...resample_poly 8 -> 6.912 Msps         70.2    28.9%
      ...notch FIR (16 threads)                24.7    10.2%
      ...read_block (disk)                     20.6     8.5%
      ...de-rotate (16 threads)                14.9     6.1%
      ...bootstrap search (once per chunk)      6.0     2.5%
    cell_pool (36 x demod + freq de-int)       26.1    10.7%
    GPU BICM (CPE + HTI + demap + de-int)      19.5     8.0%
    LDPC (74 Blocks, batched, CUDA)            12.5     5.1%
    fine_timing (41 candidates, threaded)      11.7     4.8%
    BCH + descramble (batched)                  9.2     3.8%
    ---------------------------------------------------------
    TOTAL                                     243.4     1.014x real time
    (one Frame of air = 247.1 ms)

**The new bottleneck is `scipy.signal.resample_poly`** — the 8 -> 6.912 Msps
rational resampler — at 70 ms/Frame, 29% of the whole pipeline. It is a
single-threaded C loop, and unlike the FIR it cannot be blocked without
replicating scipy's internal polyphase phase alignment and padding, which
would change the arithmetic. Three ways out, in increasing order of cost:
capture at 6.912 Msps and skip it entirely (the SDR supports it); reimplement
`upfirdn` blocking with the phase constraint that a block boundary be a
multiple of `down` and gate it byte-identical; or move it to CUDA and accept
the same hypot-class gate the fused front end accepts.

## A concurrent branch changed the oracle mid-flight

M10 added an attenuation LADDER to `PlpChain.decode_llr` — alpha tried at
1.0, 0.85, 0.75 until the syndrome is zero — while this milestone was in
flight. It is a real behaviour change: on RF33 it converges **one more FEC
Block in 35,224** than the fixed-0.75 decoder, and that one Block is one more
ALP packet and one more IP datagram.

The fast path was written against the pre-ladder reference and reproduced
`m7_long.dg` exactly, which looked like a pass. Running the CURRENT reference
end to end is what exposed it: 35,224 vs 35,223. `m9_fast.ldpc_batch` now
reproduces `decode_llr`'s semantics exactly — first converged rung wins,
otherwise the rung with the fewest unsatisfied checks — and it costs nothing,
because the first rung decodes all 74 Blocks and the retry batch is empty.

> **Law: a gate that compares against "the reference as of last week" is a
> gate against a memory. The oracle has to be re-run, not remembered — and
> when two tracks share a repository, "re-run" means today.**

## Streaming: the design, and a latency budget in measured numbers

Not built, deliberately — a clear design plus a measured budget is worth more
than broken plumbing. What the numbers say:

**Throughput is no longer the problem. Structure is.** The decode is currently
`capture -> decode all -> stitch -> walk ALP -> file`, and three of those
stages assume the whole recording is in hand.

    STAGE 1  continuous front end, ONE state-carrying process
             read -> resample -> de-rotate -> notch -> Frame queue
             measured 135.7 ms/Frame against a 247.1 ms budget = 1.82x
             NEW WORK: the resampler needs phase state across blocks (the FIR
             already blocks exactly via overlap-save, and the de-rotation
             needs only a continuous sample index).  The bootstrap search must
             run ONCE, not per chunk; after that t0 tracks by
             +FRAME_SAMPLES with the existing +-20 sample search.

    STAGE 2  frame decoder, N worker processes off the queue
             measured 107.7 ms/Frame (decode_frame + fine_timing)
             against 247.1 ms = 2.29x on ONE worker.
             Frames are independent once t0 is known, so N = 1 suffices and
             N = 2 is the headroom.

    STAGE 3  transport, one process
             measured 0.504 s for 476 Frames = 1.06 ms/Frame = 233x.
             NEW WORK: AlpWalker must become resumable across block
             boundaries.  M7 already proved why -- walking chunks separately
             lost one ALP packet per boundary and punched a 1400-byte hole in
             essentially every MPU.  The fix is a carry buffer, not a rewrite.

    STAGE 4  MMTP -> MPU -> player
             NEW WORK: emit an MPU the moment it completes instead of at end
             of file, and make m7_play.retime() incremental (it currently
             accumulates decode time over the whole recording in one pass).

    LATENCY BUDGET (steady state, N=2 workers, 1-Frame granularity)
      front end                                 0.136 s
      decode                                    0.108 s
      transport + object reassembly             0.002 s
      **MPU completeness -- INHERENT**          2.002 s
      player buffer (typical)                   0.5   s
      -----------------------------------------------
      total                                     ~2.75 s

**The dominant term is not ours and cannot be optimised away: an MPU is one
2.002-second ISOBMFF fragment and it is not decodable until its last fragment
arrives.** M7 measured this from the other side — "any capture shorter than
~3 s cannot contain a complete MPU". So ~3 s is the floor for this multiplex
whatever the decoder does, and the decoder's own contribution to latency is
about a quarter of a second. That is the number that says the streaming
rewrite is worth doing: **it is plumbing, not physics.**

The AC-4 audio gap M7 recorded is unchanged and is orthogonal.

## Assumptions added in M9

    None.  Every change is either a memoisation of a pure function, an exact
    re-blocking of an LTI filter, or a re-implementation gated to
    bit-identity against the function it replaces.  The single exception is
    named above (CUDA hypot != NumPy hypot in the fused front end) and is
    gated on byte-identical decoded output instead, over 35,224 FEC Blocks.

## Where M9 stops, and the sharpest next step

Reached and gated: 2.987x real time aggregate, 1.014x single-stream, decoded
output byte-identical to the reference chain in both of its configurations,
with the untouched chain still selectable as `--accel none`.

**Sharpest next step, in order:**

1. **Capture at 6.912 Msps.** The single biggest remaining cost, 29% of the
   pipeline, exists only because the capture is at 8 Msps and has to be
   rationally resampled. `m7_capture.py` sets the rate; if the SDR will run
   at 6.912 Msps the stage disappears rather than being optimised. Verify the
   readback (per the unverified-writeSetting law) and re-gate the whole chain
   against a resampled 8 Msps capture of the same air.
2. **The streaming rewrite**, in the four stages above. Stage 1's resampler
   state and stage 3's resumable walker are the only genuinely new pieces.
3. **Apply M9 to M10's paths.** The 64800 LDPC has 4.5x the edges of the
   16200 code and RF25/RF30 need 16..43 iterations where RF33 needs 2..5, so
   the LDPC's share of THAT pipeline is far larger than 3.9% and the GPU
   decoder should be worth much more there. `GpuMinSum` takes any
   `(checks, n)` and needs no change; only the gate has to be re-run.
4. **`m3_ldpc.min_sum_decode` is still the fallback and is still slow.** If
   the CPU path is ever the one running live, `pack()` is now cached but the
   per-block Python loop over 50 iterations is not batched. Batching the CPU
   decoder over 74 blocks with the same NumPy it already uses is a
   half-day's work and would raise the no-CUDA floor.

---

# M10 results (2026-08-06) — the SLT off a second multiplex, and a decoder constant that had been impersonating a weak signal

Offline work on captures already banked in `lab/`. **No radio was opened.**
Tools: `lab/extract_ldpc64k.py` -> `spec_ldpc64k.py`, `lab/m10_cti.py`,
`lab/m10_core.py`, plus additive changes to `m3_ldpc.py`, `m6_bicm.py` and
`m6_cells.py`.

(M9 is a separate, concurrent milestone — GPU acceleration — in `lab/m9_*.py`.)

## The headline

    A/322 Annex A.1 (Ninner = 64800)  12 tables extracted, 4 witnesses agree,
                                      72 gates pass, 2 controls fail
    A/322 7.1.4 CTI                   implemented and gated on the
                                      transmitter's OWN published formula
    m6_cells.py                       now driven from L1, not constants
    ----
    RF30 CORE LAYER DECODED.  218/224 FEC Blocks, every one BCH-zero.
    ALP -> IP -> **LLS -> SLT**.  Three services, categorised, and mapped
    against the user's list.
    ----
    RF25 does NOT decode, and the reason is measured rather than argued:
    its PLP 16 -- which sits OUTSIDE the LDM region -- decodes 15/15, so the
    whole chain is proven on that channel; PLP 0 is about 4 dB short.

## The bug worth the whole milestone: `alpha = 0.75`

The first RF30 run left the LDPC at **4 to 16 unsatisfied checks of 38880**,
on every FEC Block, and *stayed there at 300 iterations*. That is the exact
picture M8 warned about — "a near-converging LDPC is NOT evidence of weak
signal" — and it is even more convincing than M8's case, because 9 of 38880
is 0.02%.

It was not the link and it was not a spec reading. It was the **normalized
min-sum attenuation constant**, which M2 set to 0.75 and which had been right
for every code this project had ever run. The Ninner = 16200 codes RF33 carries
have check-node degrees 4..11; the 64800 rate-6/15 code RF25 and RF30 carry has
check degrees **6..8**. Min-sum's overestimate of the check message shrinks with
degree, so 0.75 over-attenuates a low-degree code and parks the decoder in a
trapping set it cannot leave.

    alpha 0.60   363 / 514 / 503 unsatisfied
    alpha 0.75     9 /  13 /  14        <- the M2 default, three blocks
    alpha 0.85     0 /   2 /   6
    alpha 1.00     0 /   0 /   0        converged in 16 / 30 / 43 iterations
    sum-product    0 /   2 /   6        (and BCH zero even at 2 and 6)

`PlpChain.ALPHA_LADDER = (1.0, 0.85, 0.75)` now tries them in order. **This
cannot manufacture a decode**: the oracle is still zero unsatisfied checks out
of tens of thousands plus a zero BCH syndrome, and no value of alpha makes a
wrong codeword pass.

> **Lesson, and it is a new shape.** M4's was "the spec is wrong is the most
> expensive conclusion available". M6's was "when prose and equations
> disagree, the equations win". This one is: **a tuning constant that was
> correct for every case you have ever run is still an untested assumption on
> the next one — and it fails by looking exactly like 0.2 dB of missing
> signal.**

## Job 1 — the Ninner = 64800 LDPC tables

`spec_ldpc.py` extracts Annex A.2 (16200) only — its own header says so — and
`m3_ldpc._params()` returned the 16200 rows *whatever `n` was asked for*, so a
64800 request got silently wrong data. Every LLS-bearing PLP on both new
channels signals `fec_type = 1`.

`extract_ldpc64k.py` mines Annex A.1 with **PyMuPDF** (M6's finding: the
"variables lost to the math font" problem is a pdftotext defect) from **four
independent editions** of A/322 — 2026-04, 2023-03, 2021 and 2018.

    all four editions:  12 tables, 7865 addresses, arithmetic sum 113641173
    G5 cross-edition identity, element for element:  PASS
    G1..G7 on the extracted tables:                  72 pass, 0 fail

The gates are closed arithmetic, not eyeballing:

| gate | what it pins |
|---|---|
| G1 | row count == Kldpc/360 (Type B) or Kldpc/360 + Q1 (Type A), derived from Section 6 and not from the page |
| G2 | every address inside [0, Nldpc - Kldpc) |
| G3 | every row strictly increasing, no duplicate inside a row |
| G4 | the 6.1.3.1 / 6.1.3.2 accumulator expansion covers every one of the N-K parity accumulators |
| G6 | Type B row weights non-increasing; Type A's trailing Q1 rows strictly lighter than every information row |
| G7 | Kldpc*15 == Nldpc*rate_numerator |

G1 landing on all twelve tables first time is the sharpest single result: the
page line-counts are 29/41/53/64/72/87/96/108/120/132/144/156 and the Section 6
arithmetic predicts exactly those, independently.

**Controls.** Six of the twelve tables are printed as two side-by-side column
blocks. Read row-interleaved instead of column-major they fail G6 — and they
fail it **exactly six times, one per two-column table**. The second control is
a single address incremented by one: **the intra-table gates are blind to it by
construction** (it stays in range, stays increasing, still covers every
accumulator), and that is recorded rather than papered over. What catches it is
G5, cross-edition agreement, and the script scores it there.

`m3_ldpc.py` now raises on an unservable 64800 request instead of lying, and
its encoder/parity-check cross-validation passes at 64800 for rates 2/15, 6/15,
7/15 and 9/15 (Type A and Type B both).

`m6_bicm.py` gained the **GF(2^16) BCH**: Ninner = 64800 uses twelve degree-16
generator polynomials, so **Mouter = 192**, not 168. The degree sum is the
identity that pins it. The full synthesized round trip — payload to cells and
back, with LDPC, BCH and descrambler — now closes at 64800 for QPSK 6/15,
QPSK 9/15, 64QAM 6/15 and 256QAM 7/15.

## Job 2 — the Convolutional Time Interleaver

A/322 7.1.4's prose is intact under PyMuPDF, and tracing one cell through it
gives an index map rather than a simulation:

    out[q] = in[q - Nrows * k_q]        k_q = (start_row + q) mod Nrows
    q(i)   = i + Nrows * ((start_row + i) mod Nrows)

**That second line did not have to be trusted, because the spec prints it.**
A/322 9.3.9.1 defines a signalled field as

    L1D_plp_CTI_fec_block_start = C + Nrows x ((L1D_plp_CTI_start_row + C) mod Nrows)

which is the same equation at `i = C`, in terms of two fields this project had
already decoded off the air. So the reading became a gate: solve for C and
require `0 <= C < cells_per_FEC_Block`.

    channel  PLP  Ninner  depth->Nrows  start_row  fec_block_start        C
    A         0   64800   3 -> 1024        344        497784       22648  of 32400
    A         1   64800   3 -> 1024        344        670000        6448  of  8100
    A        16   16200   0 ->  512        192         98304           0  of  8100
    B         0   64800   3 -> 1024       1002        773878       20214  of 32400
    B         1   64800   3 -> 1024       1002        189638        9414  of 10800

Every other entry on Table 9.26's menu (512 / 724 / 887 / 1024, plus the
extended 1254 / 1448) is run as a control and every one lands outside the
window. C is a single determined integer over a range of ~Nrows^2, so a chance
hit is about 3%; five PLPs on two independent transmitters all landing inside
is about 2e-8.

Plus a synthesized round trip (identity over the valid region for five
(Nrows, start_row) pairs), a structural check that the map is injective with
delay exactly `k*Nrows`, and — off the air — **the commutator advances by
`plp_size` per Subframe**, so `start_row` of Frame 1 must be
`(start_row(0) + plp_size) mod Nrows`. Signalled 924, predicted 924, from two
independently LDPC-decoded L1-Details.

## Job 3 — `m6_cells.py` is now driven from L1

`m6_cells.Geometry` carries FFT size, guard interval, Cred/NoC, pilot pattern,
SBS first/last and the null-cell split, symbol count, the number of Preamble
symbols **and their two different NoCs (7.2.5.1)**, the L1-Basic/L1-Detail cell
cost, the **frequency-interleaver ENABLE flag**, and the PLP table — all read
out of the decoded L1. `Geometry.rf33_legacy()` expresses the old module
constants in the same object and `selftest_geometry()` requires the L1-driven
path to reproduce all six of RF33's numbers exactly, so a refactor that changes
a number is caught rather than shipped.

It paid for itself immediately, with **two closed identities that nothing in
the pilot arithmetic could have fitted**:

* **RF30 has exactly ONE PLP and the pool has NO dummy tail.** Walking the
  pilot tables gives `2556 + 172*5999 + 2*5009 = 1044402`; the LDPC-decoded
  `L1D_plp_size` is **1044402**. Exact, with no slack anywhere.
* **`L1D_sbs_null_cells` computed == signalled** on both channels: 127 and
  **747**. The 747 is a value this project had never seen; it falls out of the
  Annex F / Table 7.5 arithmetic for 8K / Cred 3 / SP12_4.

And on RF25 the leftover **51257 dummy cells** became the referee that settled
the one genuinely open reading — the A/322 7.3 frequency-interleaver symbol
counter with a two-symbol Preamble:

    fi_offset 0   dummy sign agreement 0.5060
    fi_offset 1                        0.5060
    fi_offset 2                        0.9913   <- offset == NP
    fi_offset 3                        0.5085
    fi_offset 4                        0.5104

So M4's ASSUMPTION D2 generalises: the FI symbol counter includes **every**
Preamble symbol, not just the first. One correlation settles the FI reading,
the cell-pool origin, the PLP boundary, the SBS null split and the scrambler
at once — the same move M4 Step 3b made, now on a channel with a different
Cred, a different pilot pattern and a two-symbol Preamble.

## What the two multiplexes gave up

**RF30 — the SLT, and the M8 caveat resolved in both directions.**
218 of 224 FEC Blocks converged with zero unsatisfied checks and a zero BCH
syndrome, over 8 Frames (2 s). The ALP / IPv4 / LLS chain from M7 then ran
**unchanged**, and the LLS arrived on the well-known address in the first
2 seconds of decoded payload.

The LLS is *not* a bare SLT. It is a **multi-table container** (LLS_table_id
254) carrying `[count][id][version][length][gzip]` records plus a CMS
signature, and a second table (id 255) carrying a vendor-defined XML service
description. Inside the container are LLS_table_id 1 (**the SLT**) and
LLS_table_id 3 (SystemTime, `currentUtcOffset` 37 — which agrees with the
TAI-UTC offset M4 used to check `L1D_time_sec`).

The SLT lists **three** services: one linear (`serviceCategory` 1) ROUTE
service with **no `protected` attribute**, the ESG (`serviceCategory` 4), and
one linear MMTP service. Neither linear service is protected, so nothing here
is LOCKED. Virtual channel numbers and service names are location-identifying
and live in gitignored `lab/ground_truth.txt`.

**M8's `L1D_version = 0` caveat is now settled in both directions.** M8
recorded that RF30's trailing fields were provisional because `L1D_bsid`
decoded to an implausible 65019 and the reserved run was 0 bits. That was right
about the *trailing* fields — the SLT carries a different, plausible bsid — and
the doubt is now bounded from the other side too: **every field up to and
including the CTI parameters is validated by the payload decoding**, because
`plp_size`, `sbs_null_cells`, `CTI_depth`, `CTI_start_row` and
`CTI_fec_block_start` all had to be exactly right for 218 FEC Blocks to
converge.

**But the Core layer carries no video, and that is a finding, not a failure.**
The three flows on RF30's Core layer are the LLS, an RTP/MPEG-TS stream that is
**entirely PID 0x1FFF null packets**, and a third-party RTK GNSS correction
service. Neither linear service's SLS address appears anywhere on the Core
layer, so their media is on the **LDM Enhanced layer** (PLP 1, 64QAM-NUC 6/15) —
the one layer this milestone deliberately did not build, on the grounds that
the LLS did not need it. That reasoning was correct and the conclusion is
narrower than hoped: **the service list needed only the Core layer; the video
needs the Enhanced layer.**

**RF25 — and the umpire says it is the link.**
RF25's PLP 0 leaves ~11800 of 25920 checks unsatisfied, which is 46% — not
"nearly", but *nothing*. Two independent measurements say why, and neither is
an inference drawn from a failure:

1. **PLP 16 decodes 15 of 15, BCH zero.** It sits at cells 1133282..1165682,
   which is **outside** the LDM Enhanced layer's span, so it sees the channel
   with no injected interference — and it exercises the L1-driven cell pool,
   the frequency de-interleaver, the CTI (Nrows 512, start_row 192, C = 0), the
   BICM chain, the LDPC, the BCH and the descrambler. **The chain is proven on
   this channel.** This is M6's PLP-16-as-umpire move again, and it converted
   "RF25 fails" into "RF25's PLP 0 fails".
2. **The SNR is measured against a KNOWN transmitted sequence.** The 51257
   dummy cells are exactly `+-1` with the 5.2.3 scrambler's signs, so fitting
   and subtracting them gives the noise directly: **3.37 dB**. With the
   signalled LDM injection of 5.0 dB (Table 9.24) the Core layer's SINR is
   `0.760 / (0.240 + 0.457)` = **0.35 dB**.

And the threshold it is measured against is **this decoder's own**, not a
published figure — the same discipline M6 used to refute its link-budget
caveat. Three trials each, AWGN, Ninner = 64800:

    QPSK   6/15   3/3 at  0.5 dB      (RF30's Core layer)
    QPSK   9/15   3/3 at  2.5 dB      (RF25's Core layer)
    64QAM  6/15   3/3 at  9.0 dB      (RF30's Enhanced layer)

So RF25's Core layer is **about 2 dB short** against our own decoder, and about
4 dB short against A/327's published ~4.6 dB for the mode. Either way it is the
link. RF30's Core layer, at ~12 dB channel SNR, has SINR
`0.715 / (0.285 + 0.07)` = 3.0 dB against a 0.5 dB threshold — a 2.5 dB
margin, which is what 218 of 224 Blocks looks like.

Two further checks, because "the link is short" is exactly the comfortable
answer M6 had to reverse: 2-D pilot interpolation across the SP12_4 lattice — a
real suspicion, since that pattern leaves only 48-carrier spacing inside a
single symbol — buys **0.08 dB** (3.37 -> 3.45), so the estimator is not the
limit; and the SBS symbols' **null cells are a second, assumption-free noise
meter** (they are transmitted as zero), giving ~0.6-1.2 dB at the band edges on
RF25 against 11.4-12.1 dB on RF30.

## Two measurements this milestone leaves behind as instruments

Both are exact, and neither depends on a constellation hypothesis — which
matters, because M2 was burned comparing fits across constellation orders, and
a 256-point set fits a noise cloud better than a 4-point set by construction.
(That trap reappeared here: an LDM "composite constellation" fit reported
12.5 dB on RF30 and then made the decoder strictly worse when used as a
demapper. It is recorded as a control that failed.)

* **Dummy cells (A/322 7.2.6.5).** Known `+-1` values wherever a Subframe's
  pool is not fully used. RF25 has 51257 of them per Frame.
* **SBS null cells (A/322 7.2.6.4).** Transmitted as zero, so their received
  power IS the noise. Present on every channel that signals
  `sbs_first`/`sbs_last`. Slightly pessimistic at the band edges, because they
  sit where the receive filter rolls off.

## Where M10 stops, precisely

Reached and gated: the Ninner = 64800 LDPC tables and the GF(2^16) BCH, the
Convolutional Time Interleaver, an L1-driven cell map, the Core-layer payload
chain for an LDM multiplex, and **the SLT off a second multiplex** — which is
what the milestone was for.

Not reached: **video off either new multiplex.** Precisely:

1. **RF30's linear services are on the LDM Enhanced layer.** The canceller was
   built and run inside this session — see the addendum below — and it does
   not reach. After perfect cancellation the Enhanced layer would see
   `0.285 / 0.07` = 6.1 dB against this decoder's measured 9.0 dB threshold,
   so about 3 dB is missing; but that arithmetic is NOT what the addendum
   concludes, because the Enhanced path has no referee and therefore cannot
   be told apart from an implementation error.
2. **RF25's PLP 0 is 2-4 dB short**, depending on whose threshold you use. No
   line of code closes that.

## Sharpest next step, in order

1. **Get more UHF signal, and measure it with the instruments above rather
   than by whether something decodes.** One re-aim, or a bias-T LNA on
   Antenna B, then re-run `m10_core.py` — the dummy-cell SNR on RF25 and the
   SBS-null SNR on RF30 are calibrated numbers now, so the antenna experiment
   has a scoreboard that does not need a successful decode to read.
   +4 dB on RF25 unlocks its SLT; +2-3 dB on RF30 unlocks its video.
2. **Build the LDM canceller** — the smallest remaining rung, and testable
   *now* on RF30 even while it fails, because a correct canceller must reduce
   the residual power to exactly `a_enh^2 + N` and a wrong one cannot.
3. **Audit the other decoder constants the way `alpha` should have been
   audited.** `iters`, the `sigma2` estimator and the `|Im| < 0.35` dummy-cell
   threshold are all M2/M4-era defaults that have only ever been exercised on
   one multiplex.
4. **Re-capture RF29 at night**, still outstanding from M8.

## Assumptions added in M10

    A1  Annex A.1's two-column tables are read COLUMN-MAJOR.  GATED, not
        assumed -- the row-interleaved reading fails G6 on exactly the six
        two-column tables and on no others.
    K1  Ninner = 64800 uses Mouter = 192 (a GF(2^16) BCH).  GATED by the
        degree sum of the twelve extracted generator polynomials and by the
        synthesized round trip.
    C1  A/322 7.1.4's cell trace, out[q] = in[q - Nrows*k_q].  SETTLED BY THE
        SPEC'S OWN PUBLISHED FORMULA (9.3.9.1) on five PLPs across two
        transmitters, with the whole Table 9.26 menu run as controls.
    D2' the A/322 7.3 frequency-interleaver symbol counter includes EVERY
        Preamble symbol (offset == NP), generalising M4's D2 from NP = 1.
        SETTLED BY THE AIR: dummy-cell agreement 0.9913 at offset 2 and
        0.506-0.510 at every other offset.
    M1  the min-sum attenuation is not a constant of the receiver.  It is now
        a ladder terminated by the syndrome, and the syndrome is still the
        only oracle.

## Addendum — the LDM canceller was built, and it stops UNGATED, not "3 dB short"

The task said not to build LDM cancellation unless the Enhanced layer turned
out to carry a channel the user wants. It does, so it was built.

The mechanism is small, because the two layers share a time interleaver: both
PLPs signal the same `CTI_depth` and the same `CTI_start_row`, so **one
de-interleave permutation serves both**, and cancellation can be done entirely
in the de-interleaved domain:

    residual[i] = cells[i] - a * reencode(decode(cells[i]))
    Enhanced FEC Blocks then sit at i = C_enh + m*10800,
    where C_enh comes from the ENHANCED PLP's own L1D_plp_CTI_fec_block_start
    solved through the same 9.3.9.1 identity (9414, and it passes that gate).

The Core reconstruction is real: six Core FEC Blocks decode, re-encode, and fit
the received cells at |a| = 0.903, leaving a residual of power 0.75-0.81.

**The Enhanced layer does not decode: 18300-18600 of 38880 unsatisfied, which
is 47% — chance, not "nearly".** And the residual will not resolve into the
signalled constellation either. All four candidates below have exactly 64
points, so this is not M2's constellation-order bias:

    candidate                MER dB    occupancy H   chi2/dof
    NUC_64_6/15 <- signalled   9.06        5.8612       101.6
    NUC_64_11/15               9.35        5.8892        83.6
    uniform 64QAM              8.70        5.7959       161.6
    NUC_64_2/15                5.42        5.2303       703.2
    Gaussian noise (floor)        --            --        57.1

Compare M4 Step 3a, where the signalled constellation won all three statistics
and the runner-up was 4x worse on chi-square. Here **the signalled table loses
to a table the transmitter is not using**, and the whole family sits close to
the pure-noise floor. The only thing this test rejects is 2/15, which is
grossly different in shape.

**So the honest statement is: the Enhanced layer is UNGATED.** Nothing has
verified the Enhanced path — not the constellation, not the anchor, not the
canceller's residual — and 47% unsatisfied is exactly what *both* "buried in
noise" and "a wrong rung" produce. The link arithmetic (6.1 dB available vs a
measured 9.0 dB threshold) is consistent with the first, but the project's own
history says an arithmetic that is merely consistent is not a diagnosis: M8's
98-of-12960 was consistent with a weak link too, and it was a missing term.

**What would tell them apart, and it is the sharpest next step for this rung:**
a synthesized LDM round trip — generate Core QPSK 6/15 and Enhanced
64QAM-NUC 6/15 at the signalled injection level, push both through the same
CTI, cancel, and require the Enhanced layer to decode noiselessly. That is a
round trip and therefore proves only that the code inverts itself (M6's law),
but it is precisely the discrimination needed here, and it costs no signal.
Until it exists, "RF30's Enhanced layer is 3 dB short" is a hypothesis, not a
result, and it is recorded as one.

---

# M11 results (2026-08-06) — **LIVE TV.** The pipeline is continuous, and the radio was asked for the rate that deletes a stage.

The radio WAS opened: `radio_lock` owner `atsc3_watch`, priority 60, polite
wait, heartbeat on a timer inside the read loop, released in `finally`,
outside the meteor window. Tools: `lab/m11_stream.py` (the four pieces),
`lab/m11_watch.py` (the driver and the telemetry), `lab/m11_gate.py` (the
gates), `lab/m11_rateprobe.py` (the rate readback), `atsc3/__main__.py` (the
front door).

## The headline

    python -m atsc3 watch --rf 33

    LIVE, 2026-08-06 13:11:39 -> 13:16:02
      wall                    263.217 s
      Frames                  1048        air 258.972 s
      FEC Blocks              77,547 / 77,552 converged, all BCH zero
      ALP resyncs             5
      IP datagrams            91,216
      MPUs                    512 complete, 8 lost (1 of them video)
      to the player           125 segments = 254.339 s of 720p60 HEVC
      **UNDERRUNS             0**        minimum buffer 1.757 s
      queues (raw/frame/bb)   0 / 0 / 0 for the whole run

    the recording, re-decoded by ffmpeg to a null sink:
      15,240 frames, 254.338667 s, hevc Main 10 1280x720 60000/1001,
      **zero error lines**

**Four minutes and fourteen seconds of continuous watchable video off the air,
with nothing buffered ahead and nothing dropped.**

## x-real-time on a radio is capped at 1.000x, and saying so matters

The run reports 0.9839x sustained. That is not a shortfall: **the air arrives
at 1x and that is the entire supply.** A live pipeline cannot decode faster
than the transmitter transmits, so "sustained x-real-time" — the number M9
optimised — stops being the right instrument the moment the source is an
antenna. The 4.2 s deficit is acquisition (one bootstrap search, plus one
startup overflow that cost a re-acquisition), not backlog: every queue in the
pipeline sat at zero for the whole run.

The honest headroom number is the stage cost against the Frame budget, and it
is measured live:

    ms per Frame, live, budget 247.11 ms
    -------------------------------------------------------------
    front end        69.5   notch 34.2  derotate 15.1
                            fine_timing 14.6  resample 4.7  window 1.0
    decode           80.0   gpu_bicm 29.2  ldpc 23.0  cell_pool 15.6
                            bch_descramble 10.7  dummy_check 1.5
    -------------------------------------------------------------
    HEADROOM         slowest stage 3.09x, both stages summed 1.65x

Compare M9's single-stream table, which totalled 243.4 ms against the same
247.1 ms budget — 1.014x, i.e. no headroom at all. Two things bought the
difference: the resampler is gone, and the front end and the decoder are now
separate stages that overlap instead of one serial pass.

## The free 29% was free. The radio does 6.912 Msps.

M9 named this as the sharpest next step: `resample_poly` was 70.2 ms of a
243.4 ms Frame, 29% of the pipeline, and it existed for exactly one reason —
`m7_capture.py` asks for 8 Msps and the decoder works at 6.912.

Asked, and **read back**, because a `setSampleRate` that is not read back is a
hope:

    requested   6912000.0 Hz
    readback    6912000.0 Hz          EXACT
    bandwidth   6000000.0 Hz
    capture     55,296,000 samples / 7.985 s wall, ratio 1.0019,
                0 overflow reads      CAPTURE-INTEGRITY PASS

    that capture decoded: 27 Frames, 1998/1998 FEC Blocks, 1998/1998 BCH
                          zero, 0 ALP resyncs, 7 flows

Worth recording: `listSampleRates` does **not** offer 6912000 — it lists a
tidy ladder from 62.5 kHz to 10 MHz. `getSampleRateRange` reports a continuous
2 MHz–10.66 MHz span, and the readback of the off-ladder request is exact. A
receiver that had trusted the enumeration would have concluded the rate was
unavailable and kept a stage it did not need.

So the resampler is not optimised, it is **absent**. It is still built and
still gated, because banked 8 Msps captures replay through the same front end
and because a deleted fallback is not a fallback.

## The four pieces, and the gate on each

M9 said the streaming rewrite was plumbing rather than physics and named four
genuinely new pieces. All four are built, and each is gated against the batch
path **re-run today**, not against a stored expectation.

### 1. Resampler phase state

`scipy.signal.resample_poly` computes `y[m] = sum_j h[m*down - j*up] * x[j]`
and keeps `y[n_pre_remove:]`. Blocking it exactly needs `len(h)` samples of
real input history — the overlap-save argument M9 already proved for the FIR —
**and a block boundary at which the polyphase phase is zero.** Output `m` sits
in branch `(m*down) % up`, so a local call starting at input index `n0 - H`
reproduces both the branch and the summation order iff `(n0 - H) * up` is
divisible by `down`; with `gcd(up, down) == 1` that is `n0 % down == 0` and
`H % down == 0`. Blocks are therefore constrained to multiples of `down`
rather than to a convenient power of two, and that constraint is the whole
correctness argument.

    PASS  8 -> 6.912 Msps:  1,727,989 of 1,728,000 samples bit-identical
    PASS  10 -> 6.912 Msps: 1,036,789 of 1,036,800 samples bit-identical

(The 11-sample shortfall is `n_pre_remove` exactly: the streaming output is a
prefix of the batch output, which is what a stream should be.)

### 2. One-shot bootstrap, then t0 tracking

The search runs **once**, fixes `f0` and the CFO, and after that
`centre = t_prev + FRAME_SAMPLES` refined by the same ±20 fine-timing scan the
batch chain uses. The de-rotation carries a continuous sample index from `f0`
— which is what makes the phase continuous across a block seam — and the
notch carries `len(b)-1` samples of real history.

    PASS  FrontEnd vs load_fast, 6.912 Msps:  6,226,315 complex samples
          bit-identical, 854,000-sample blocks, one bootstrap search
    PASS  FrontEnd vs load_fast, 8 Msps:      6,326,393 complex samples
          bit-identical, resampler in circuit

### 3. The resumable ALP walk, and the trap that was avoided rather than fixed

M7 paid twice on this layer: walking chunks separately lost one ALP packet per
boundary and punched a 1400-byte hole in essentially every MPU, and before
that a mis-read `header_mode == 1` length swallowed 62% of the media flow
inside one "valid" packet.

A streaming walker invites a third failure that neither of those was: it can
be forced to *guess*. If it parses at `pos` with only part of the stream in
hand, three of `AlpWalker`'s own checks change answer — the length fit, the
anchor contradiction, and the resync target. And `_segment` **writes into
`self.segbuf`**, so a speculative parse that is later abandoned is not wasted
work, it is corrupted state.

So the walker does not speculate at all. It refuses to attempt a reading until
it holds `MAX_ALP + 8` bytes past the cursor **and** the Baseband anchor list
covering them. Under that precondition every length check has the answer it
has on the whole stream (nothing legal can fail to fit), and every
`_anchor_clean` sees every anchor that could contradict the candidate, because
anchors arrive in increasing order — so "all anchors below the last one" is
"all anchors". The cost is about 65 kB of lookahead, two thirds of a Frame,
and it vanishes inside the 2.002 s an MPU takes to exist at all.

    PASS  6.30 MB of real Baseband stream, 4361 anchors, fed in 41 chunks:
          5589 ALP packets vs 5589 from the whole-stream walk,
          0 resyncs vs 0

> **Law: a streaming parser's correctness argument is about what it refuses to
> do before it has enough, not about how cleverly it recovers afterwards. If
> the parse mutates state, "try it and back out" is not available.**

### 4. Emit-MPU-on-complete

The completeness test is cheap (`moof` known, and the accumulated MFU payload
has reached the declared total) but the **bytes are not reimplemented**:
`m7_objects.MmtpFlow.build()` is called on a single MPU. There is no second
reassembler, so there is nothing to drift. Retiming is incremental —
`m7_play.retime()` per fragment with a running `base_time` — instead of one
pass over the whole recording.

The video asset is identified by reading the `moov` `hdlr` handler_type off
the air and requiring `vide`. **No address, port or packet_id is hard-coded
anywhere in the streaming path**, which is the difference between a receiver
and a script that works on one multiplex.

### The end-to-end gate

    reference (m9_decode.py, re-run now)  2,217,999 bytes
                                          sha256 3bd78cec0ed24053...
    streaming pipeline, same air          2,217,999 bytes
                                          sha256 3bd78cec0ed24053...
    IDENTICAL

The reference is a **subprocess**, not a file. M9 recorded why: it reproduced
M7's banked `m7_long.dg` byte for byte, which looked like a pass, while the
current reference had gained an attenuation ladder and converged one more FEC
Block. `m11_gate.py` also reads the frame count back out of the reference run
and asks the streaming path for exactly that many, because the two paths have
different end-of-FILE rules — `m9_decode` requires a 400,000-sample tail
margin, the streaming front end asks only whether this Frame's own window is
complete — and a tail difference that is a file artefact should not be allowed
to look like a decode difference.

## Two things that were wrong first, and the measurements that said so

### The front end has to be threaded INSIDE the stage

The first working build ran at **0.558x real time** with every downstream
queue empty: the decoder was starving. The notch and the de-rotation were
plain `lfilter` and a plain `np.exp`, single-threaded, 82.6 ms and 56.0 ms per
Frame against a 247.1 ms budget for the whole pipeline. Restoring M9's
`par_lfilter` and `par_apply` — the same overlap-save FIR and the same chunked
de-rotation, still bit-identical — took the same file to 1.41x.

This is M9's "threads, not processes" argument arriving at its real
destination. It was a throughput argument there; here it is structural. **A
live pipeline is one stream, so the only parallelism available is the
parallelism inside a stage.**

Block size mattered more than expected and was swept rather than guessed:

    block (samples)   x real time   front end ms/Frame
    250,000            1.155        178.1
    500,000            1.283        154.8
    854,016 (half)     1.409        137.0
    1,708,032 (Frame)  1.356        125.5

Half a Frame wins: smaller blocks pay thread dispatch on every stage, larger
ones buy nothing and cost latency. 854,000 is the default (a multiple of 125,
so the 8 Msps resampler path stays aligned).

### A 2.0 s player buffer is exactly one MPU, and that is not a buffer

The first live run held real time and still logged **4 underruns**. The cause
was not throughput — it was that media arrives in indivisible 2.002 s units,
so with a 2.0 s prebuffer the occupancy sawtooths between 0 and 2 s and sits
on zero at the bottom of every tooth. One lost MPU, or a few tens of
milliseconds of jitter, and it goes negative.

Raising the prebuffer to 5.0 s — three MPUs — took the *same air, the same
signal quality* from 4 underruns to **0**, with the minimum buffer at 1.757 s.

> **Law: a jitter buffer smaller than one delivery unit is not a small buffer,
> it is no buffer. Size it in units of what actually arrives, not in seconds
> that feel comfortable.**

The occupancy meter is worth stating because it is the honest kind. A pipe
cannot tell us the player is starving, so it is not asked: each segment is a
known number of seconds of media, and once playback has begun the player has
consumed one second per second of wall clock, so

    occupancy = seconds pushed - seconds of wall since playback started

is what the player still holds. When it goes negative the player stalled for
exactly that long, and `t_play` is shifted by the deficit — recording the
stall rather than forgetting it (every later segment then looks early) or
carrying it forward (every later segment then looks late). An earlier version
reset the clock instead and double-counted: after a stall the occupancy sits
exactly on zero, so the *next* on-time segment re-triggered. That is why the
test is `< 0` and not `<= 0`.

## Where the losses actually come from

Over 263 s: 5 FEC Blocks of 77,552 failed to converge, and they cost 5 ALP
resyncs, 459 lost MMTP packets and 8 incomplete MPUs — **one of them video.**
The amplification is worth naming: one failed FEC Block is one missing
Baseband Packet, which is a hole in the ALP stream, which the walker correctly
refuses to parse across, which loses ~1400 bytes of an MPU, which makes a
2.002-second video fragment incomplete. **A 0.006% Block failure rate becomes
a 2-second gap.**

The 5 s buffer absorbed the one video gap without an underrun, which is what a
buffer is for. But the ratio says where the next work is: the decoder is not
the limit and neither is the plumbing — the limit is that MPU completeness is
all-or-nothing.

## Radio discipline, as run

    owner        atsc3_watch, priority 60 (a human is watching: above lab
                 background, below a satellite pass or a window hold)
    acquire      polite wait, never a seizure; refuses if outranked
    heartbeat    on a 5-second TIMER inside the read loop, never per read
                 (the 8/03 law -- per-read file I/O chopped a stream and
                 produced a fake 97 wpm)
    yield        should_yield() checked on the same timer; the run ends
                 cleanly rather than being killed
    release      in `finally`
    window       refuses inside 02:10-05:40

Two SDR overflow reads occurred across 4m20s (one at stream start, one during
shutdown). An overflow is a hole in the sample stream and t0 tracking cannot
survive one, so it forces a bootstrap re-acquisition rather than being
counted and ignored — and the same is true of a raw block dropped because a
queue was full, which is why that path raises a continuity break instead of
incrementing a counter. The startup one is why the sacrificial warm-up read is
now 0.5 s rather than 0.25 s.

## Assumptions added in M11

    None.  Every stage is either bit-identical to the function it replaces
    (resampler, de-rotation, notch, fine timing, frame decode) or is the
    existing function called on a smaller unit (MPU reassembly).  The one new
    behavioural choice -- the walker's lookahead precondition -- is argued to
    produce the SAME readings rather than acceptable ones, and gated on
    5589/5589 packets against the whole-stream walk.

## Where M11 stops, and the sharpest next step

Reached and gated: a continuous receiver that tunes, decodes and plays live
video for as long as it is left running, at 3.09x headroom on its slowest
stage, with output byte-identical to the batch chain re-run today.

**Sharpest next step, in order:**

1. **Partial-MPU repair.** One failed FEC Block currently costs a whole
   2.002 s fragment. `m7_play.trun_trim` already trims a fragment to whole
   samples; driving it from the *first short sample* instead of from byte
   availability would turn a 2-second freeze into a fraction of a second of
   lost motion. It must be labelled truncated in the telemetry, per M7's rule
   that nothing is patched silently.
2. **Audio.** The AC-4 gap M7 recorded is unchanged and is orthogonal: the
   audio MPUs reassemble COMPLETE and their ISOBMFF is valid, and this
   ffmpeg has neither an AC-4 decoder nor a muxer tag. The receiver is not
   the blocker.
3. **Channel change without a restart.** Re-tuning currently means a new
   process. The front end already re-acquires cleanly on a continuity break,
   which is most of the mechanism.
4. **The other multiplexes.** M10 drove `m6_cells` from L1; the streaming
   front end still uses RF33's module constants for `FRAME_SAMPLES` and the
   Frame window. Wiring `Geometry` through is what makes `watch --rf 25` mean
   anything.

---

# M12 — the two-second freeze, and a gate that was right to fail

M11 named partial-MPU repair as the sharpest next step. This is it, plus a
finding that arrived unasked: **two of M11's six gates do not reproduce.**

## The gate that was right to fail

Re-running `m11_gate.py` on the committed M11 tree — before touching anything —
fails `resampler_8` and `resampler_10`. Not marginally: **1,580,070 of
1,727,989 output samples differ**, at ~1 ulp of float32. `array_equal` was
doing its job; the reported PASS does not reproduce and the cause is real.

`PolyResampler` designed its taps the obvious way — `firwin(...) * up` in
float64, cast later. `scipy.signal.resample_poly` does the opposite:

```python
h = firwin(2 * half_len + 1, f_c, window=window)
h = xp.asarray(h, dtype=x.dtype)     # complex64 for a complex64 stream
h *= up                              # ...and only THEN scale
```

So scipy convolves with **complex64** taps and accumulates in single precision.
Ours accumulated in double. **Our version was more accurate, and therefore
wrong** — the gate asks whether the streaming path reproduces the batch path,
not whether it improves on it. Matching scipy's dtype *and its order* takes the
mismatch to **0 of 1,727,989**, both rates.

This is the same lesson M9 recorded when folding `/fs_post` into a scalar moved
results by 1.9e-13, one level lower down: **precision is part of the
specification when bit-identity is the contract.**

### …and the 85 samples, closed the same afternoon

`front_end` and `front_end_8msps` also failed: **85 samples of 9,726,315,
max |Δ| 1.11e-16.** Bisected rather than guessed, and two hypotheses died on
the way — it is not the threading (single-threaded reproduces it exactly) and
not the CFO (both acquisition paths return bit-identical `fine_cfo_hz`).
De-rotation alone is **0 mismatches in 10,368,000**. The notch is the whole of
it, and the mismatches all fall inside the **first 257 output samples**:

`lfilter` on the whole array starts from an *implicit* zero state, so its first
256 outputs sum FEWER terms. The streaming path prepended 256 fabricated zeros,
making it sum all 257 — the same summands with the exact zeros included, but a
different pairwise grouping, and therefore a different last bit. Reproduced
synthetically at 185 of 1,000,000, every one inside the transient; priming the
first block from the implicit state instead of padding it gives **0**.

So overlap-save is exact once the history is REAL — the M11 comment was right
about steady state and silent about the one block where the history is
invented. All six gates now pass, and that is the first time the claim has
been true rather than reported.

## The repair

One failed FEC Block used to cost a whole 2.002 s fragment: `MpuStreamer`
dropped any MPU whose media had a hole. But `m7_objects.build` already lays
every sample out at its declared size and zero-fills what is missing — so the
box offsets are the transmitted ones, and everything before the first short
sample is playable video sitting in the buffer.

`build` now reports `short_first`; the streamer cuts there and hands over what
arrived. Three things had to be right:

1. **Cut at the first short SAMPLE, not at byte availability.** `trun_trim`
   asks how many whole samples fit in the bytes received, which is the right
   question only for a clean tail. On air the loss is a *hole*. `trun_keep`
   is the sample-indexed cut.
2. **Fix BOTH trafs.** The mdat is [all media][all MMT hint], and the hint
   traf's `data_offset` follows the media block — truncate media without
   moving it and the hint samples point into the middle of the video.
   `trun_keep(seg, n_samples)` returns the segment byte for byte, which is
   what proves the rewrite is a no-op when nothing is cut.
3. **Emit the repair BEFORE the MPU that overtook it.** MMTP delivers an MPU's
   packets before the next one's, so a still-open older MPU has stopped
   growing the moment a newer one completes. Judging it at reap time instead
   would hand the Transport a fragment that arrives after its own successor —
   dropped as out-of-order, and the repair becomes a silent no-op.

And the timeline stays the air's: `retime` runs on the untouched sample table,
so `base_time` advances by the full 2.002 s. The missing tail is a **freeze in
place**, not a splice that pulls everything after it early.

## Gated against the loss it fixes — `m12_repair_gate.py`

Clean air cannot test this: a 90 s replay of the banked RF33 capture decodes
35,890/35,890 FEC Blocks with zero resyncs. So the loss is manufactured — a
contiguous run of IP datagrams dropped out of the middle, which is what a burst
of failed Blocks looks like by the time it reaches the transport.

| gate | result |
|---|---|
| clean air, repair on == repair off, byte-identical | PASS, 32 segments |
| repaired output > its own `--no-repair` control | PASS |
| ffmpeg decodes the repaired file | PASS, **0 error lines** |
| repaired duration == clean duration | PASS, 64.085 s vs 64.085 s |
| `frames_recovered` == frames ffmpeg actually gains | PASS, 17 == 17 |

Swept over 9 hole positions, recovery is **17 to 113 of 120 frames, mean
61.4** — i.e. a uniformly-placed hole costs on average *half* a fragment
instead of all of it. A 2.002 s freeze becomes ~0.98 s of lost motion.

`--no-repair` is kept as a CLI flag because a control that cannot be re-run is
not a control.

## Where M12 stops

Unchanged from M11: AC-4 audio (a tooling gap, not a receiver gap), channel
change without a restart, and the other multiplexes — the streaming front end
still uses RF33's module constants where M10 drives `m6_cells` from L1.
Added: the 85-sample front-end divergence, and the fact that an unrepairable
MPU still advances the clock in `_segment` but a reaped one never reaches it,
so gap accounting is honest per-segment and approximate per-outage.

## Addendum — the gate-vacuousness audit M7 asked for, done

M7's fourth next-step was "don't trust a gate that cannot fail — `mdat_got ==
mdat_declared` became vacuous the moment the reassembler started zero-filling
missing samples to their declared size. Audit the others for the same defect."
Audited; **no vacuous gate found**, and the reasoning is worth recording so it
is not re-done:

* **MMTP completeness** never rests on the zero-filled quantity. `mdat_got ==
  mdat_declared` survives only as one conjunct of `media_got == media_want and
  n_short == 0`, and those two are computed from the per-sample `have` bitmap,
  which zero-filling cannot flatter.
* **ROUTE object completeness** (`MmtpFlow`-adjacent `assemble`) writes into a
  `BytesIO` that zero-fills on `seek` past the end — so `got == declared` alone
  would be exactly the same defect. It is conjoined with `not gaps`, and `gaps`
  is accumulated explicitly at each discontinuity. Bites.
* **M10's P1** looks one-sided where a dummy tail exists: it only tests
  `pool - core_total >= 0`. It is not, because `dummy_agreement` (P6)
  independently referees that tail against the scrambler sequence. The identity
  is two-sided; it is just split across two gates, as P1/P6's own text says.
  On RF30, where there is no dummy tail, P1 alone is exact.
* **M10 already prints `*** CONTROL PASSED -- gate vacuous ***`** if any of its
  three deliberate controls succeeds, which is the discipline itself in code.

What the audit did NOT cover, and should not be assumed clean: the M11/M12
gates were audited the hard way instead — four of them were observed FAILING
today, which is a stronger proof of falsifiability than reading them.

---

# M13 — Subframe 1 / PLP 1: the blocker was a hypothesis nobody had tested

Every service decoded before this lived in Subframe 0. The Frame has two, and
the second is **four times the size of the first**:

```
bootstrap  13824
Preamble    1 x ( 8192 + 1536)     L1-Basic + L1-Detail
Subframe 0 35 x ( 8192 + 1536)     PLP 0 (64QAM 11/15) + PLP 16
Subframe 1 75 x (16384 + 1536)     PLP 1  -- never read until now
```

M6 and M7 both named Subframe 1 as the next service and both stopped at the
same recorded sentence: *"A/327 Section 6 specifies a frequency de-interleaver
RESET rule for the second Subframe which this project has not implemented;
missing it would break PLP 1 specifically."* That sentence was carried forward
through two milestones as a known blocker. **It was never tested.**

## The parameters were read, not guessed

`Geometry.from_l1` flattens the L1-Detail field list **by name**, which lets
Subframe 1 silently overwrite Subframe 0 and puts every PLP of every Subframe
into one undifferentiated list. Harmless while only Subframe 0 is decoded;
wrong the moment it is not. `m13_sf1.scope()` reads the loop PATH instead
(`i=1/`, `i=1/j=0/`), and out falls: 16K FFT, GI 1536, 75 symbols, SP4_4,
Cred 0, SBS both ends, frequency interleaver ON, one PLP at 256QAM 11/15 with
Ninner 64800.

Two identities closed before a single cell was demodulated, and both are the
kind that cannot be fudged:

* 947700 cells / (64800/8) = **117 FEC Blocks exactly**, and 117 = 3 x 39
* **null cells computed 1609 == signalled 1609**

and one more on first contact: **pool 956179 == predicted 956179**, dummy 8479.

## The `+1` convention came from a working decoder, not from the spec text

`L1D_plp_HTI_num_fec_blocks` could plausibly count per TI Block or per
interleaving frame, and the readings differ by 3x. PLP 0 settles it without an
argument: its L1 says `num_ti_blocks 1 / num_fec_blocks 73`, and the decoder
that has been producing watchable video uses `nti=2, n_fec=74,
ncols = n_fec_max // nti = 37`. So the fields are (value + 1) and count the
whole interleaving frame. Applied to PLP 1 that gives 3 TI Blocks x 39 columns
— which is independently what M6 measured off the waveform. **This is an
inference from a working system, not from the standard, and is labelled as
such.**

## The reset rule, settled by experiment

`cell_pool_g` already exposed `fi_offset`, the frequency interleaver's symbol
counter origin, so the question was never a blocker — it was a one-line sweep
with a decisive referee: LDPC convergence and BCH syndrome over 117 Blocks.

| symbol-counter origin | result |
|---|---|
| **0 — RESET at the Subframe** | **117/117 converged, 117/117 BCH zero** |
| 1 (Preamble counts as 0, as in Subframe 0) | 0/8, median 7950 unsatisfied |
| 35 | 0/8, median 7949 |
| 36 (counter continues across the Frame) | 0/3, median 8018 |
| 37 | 0/3, median 7951 |

The counter **resets at the Subframe boundary**. A wrong interleaver is not
subtly worse, it is total — which is what makes the controls worth more than
the result. Over 4 Frames: **468/468 FEC Blocks BCH zero**, then 1731 ALP
packets with **0 resyncs** and 1731 UDP datagrams.

## What is actually on it

Five flows, all **ROUTE/LCT** — a different transport from Subframe 0's MMTP,
and the reason M7 could not find them: it was looking on the wrong Subframe,
not at the wrong addresses. The service list table confirms the split exactly:
one service in this multiplex signals `slsProtocol = 2` (MMTP) and it is the
one already playing; **every other service signals `slsProtocol = 1` (ROUTE)
and lives on PLP 1.** Four additional television services, each 1920x1080p60
HEVC (`hvc1.2.6.L123`) with two AC-4 audio languages and captions — higher
resolution than the MMTP service.

**Encryption is per service and must be checked per service.** The first flow
to yield a complete SLS carries `ContentProtection` with `cenc:default_KID`
and `pssh` — Widevine. That one is LOCKED and is not touched, per the standing
rule. Whether the others are clear is a question for their own signalling, not
an assumption in either direction.

## Where M13 stops

The receiver work for Subframe 1 is done and gated, and it is done for ALL
PLPs on it, not one. What remains is not physical layer: per-service SLS
capture (an SLS appears about once a second, so it needs seconds of air, not
Frames), then ROUTE/DASH object assembly for any service that proves clear.

## Addendum — the per-service answer, and a second watchable channel

Sixteen Frames (**1872/1872 FEC Blocks BCH zero**) gave every ROUTE service on
PLP 1 time to repeat its SLS — an SLS appears roughly once a second, so this
was a question of seconds of air, not of decoding. All four captured:

| service | ContentProtection | KID | pssh | verdict |
|---|---|---|---|---|
| A | 0 | 0 | 0 | **CLEAR** |
| B | 6 | 6 | 6 | Widevine — LOCKED, untouched |
| C | 6 | 6 | 6 | Widevine — LOCKED, untouched |
| D | 6 | 6 | 6 | Widevine — LOCKED, untouched |

Three of the four are encrypted and are not pursued. **The fourth is in the
clear**, and it is the one an independent commercial receiver also reports as
unencrypted — two instruments agreeing, which is worth more than either alone.

Init segment + one complete media segment, concatenated: **1920x1080 HEVC,
zero ffmpeg error lines.** Stitching a second, truncated segment through
M12's `trun_keep` brought it to **145 frames** — so the partial-MPU repair
built for MMTP generalises unchanged to ROUTE/DASH, which was not designed for
and is worth recording: both transports carry the same ISOBMFF fragment shape,
so the repair lives at the right layer.

This multiplex therefore yields **two watchable services on two different
transports**: the MMTP one at 720p60 that `watch` already plays, and this
ROUTE/DASH one at 1080p60. The remaining work for the second is plumbing, not
physics: object assembly into a continuous stream, and the same AC-4 audio gap
as everywhere else.

---

# M14 — what actually buys the missing dB (analysis, `lab/m14_budget.py`)

M10 recorded that the non-decoding multiplex's Core layer is "about 2 dB short
against our own decoder, and about 4 dB against A/327's published figure."
Both numbers are right and both are **SINR shortfalls** — and quoting either as
"how much antenna do I need" is wrong, because of a structural fact the phrase
hides.

## The ceiling

That multiplex is LDM. Two layers share the same cells and the Core layer is
decoded FIRST, with the Enhanced layer treated as interference; there is no
ordering in which the Enhanced layer can be cancelled beforehand, because
cancellation runs the other way. So

    SINR_core = P_core / (P_enh + N)

and **P_enh does not shrink when the antenna improves.** As `N -> 0` the SINR
approaches `P_core / P_enh`, which is exactly the signalled injection level:

    injection 5.0 dB -> P_core 0.7597, P_enh 0.2403
    channel SNR 3.37 dB (measured off the dummy cells) -> N 0.4603
    Core SINR now      +0.35 dB      (reproduces M10's figure exactly)
    CEILING as N -> 0  +5.00 dB      the injection level, and unimprovable

## The exchange rate, and why it decides the shopping list

Only the `N` term responds to the antenna, and `N` is barely two thirds of the
denominator — so **every dB of antenna buys less than a dB of SINR**, and the
closer the target sits to the ceiling the worse the exchange rate gets:

| threshold | SINR short | antenna actually needed | headroom under ceiling |
|---|---|---|---|
| our own decoder, 2.5 dB | 2.15 dB | **3.91 dB** | 2.5 dB |
| A/327 published, 4.6 dB | 4.25 dB | **12.98 dB** | 0.4 dB |

So **which threshold is correct changes the answer, not just the margin.**
Against our own decoder the channel is reachable with a real but ordinary
antenna improvement. Against the published figure it is effectively
unreachable, because only 0.4 dB of headroom exists beneath the ceiling.

## The menu

| lever | worth |
|---|---|
| antenna / preamp | **THE lever** — the only one that moves `N` at all; ~3.9 dB needed |
| matched filter | **0 dB** — for OFDM the FFT already IS the matched filter per subcarrier |
| more FEC / iterations | **0 dB, measured** — already 50 iterations, already 16–43 needed here; on the waterfall, not short of effort |
| 2-D pilot interpolation | **0.08 dB, measured** — the estimator is not the limit |
| widely-linear equaliser | **predicted NEGATIVE** — WL wins at a multipath cliff and loses in near-AWGN below ~17 dB; this is 3.37 dB of nearly pure AWGN |
| LDM cancellation | **0 dB for the Core layer, by ordering** — it helps the layer we do not need |
| alphabet-aware LLRs | **UNKNOWN — the only software candidate not yet measured to be worth zero** |

The last one is the only open question and it deserves an experiment rather
than an estimate. Standard max-log demapping treats the Enhanced layer as
Gaussian. It is not: it is a known finite constellation at a known relative
power, so marginalising over its alphabet gives correctly-scaled LLRs into the
LDPC. `m6_bicm.roundtrip` already synthesises a full chain through AWGN, so the
experiment is: synthesise Core + Enhanced at the measured injection, find each
demapper's threshold by the same 3/3 convention, report the difference.
Radio-free, offline, and falsifiable — and it should be run **before** any
hardware is bought, because it is free and everything else on the menu has
already been measured to be worth nothing.

---

# M15 — the last lever, measured: alphabet-aware LLRs buy nothing here

M14 left exactly one candidate unmeasured: the Core demapper treats the LDM
Enhanced layer as Gaussian noise, and it is not — it is a known finite
constellation at a known relative power. `lab/m15_ldm_llr.py` measures it, with
a control designed so a null result would still mean something.

**Prediction, recorded before the run:** small or zero, because a 256-point
interferer is already nearly Gaussian and that is precisely the condition under
which the existing approximation is right.

## Three demappers, same LDPC

`gaussian` folds the interferer's power into the noise variance (what ships).
`alphabet` marginalises with **max-log** — the minimum distance over the
interferer's alphabet, matching the convention used everywhere else in this
project. `exact` marginalises properly, `logsumexp(-d²/N0)` over the alphabet.
All three feed the same `decode_llr`, so the only difference in the experiment
is the LLR vector. Thresholds by the same 3/3-trials convention as every other
threshold recorded here.

| | gaussian | alphabet (max-log) | exact (log-sum-exp) |
|---|---|---|---|
| **CONTROL** — QPSK interferer (lumpy) | 7.00 dB | 6.50 dB → **+0.50** | 6.50 dB → **+0.50** |
| **REAL** — 256QAM-NUC interferer | 6.80 dB | 8.00 dB → **−1.2** | 6.80 dB → **+0.00** |

## Three findings, in order of usefulness

**1. The answer is no, and it is a tight no.** On the real channel the exact
finite-alphabet demapper and the Gaussian approximation threshold at the *same*
0.2 dB grid point. The gain is below 0.2 dB. The prediction was right for the
stated reason: 256QAM-NUC is Gaussian enough that modelling it exactly wins
nothing. **M14's menu now has no unmeasured entries, and the antenna is the
only lever.**

**2. Max-log marginalisation over a dense interferer is actively HARMFUL — not
neutral, −1.2 dB.** Taking the *best* interferer symbol for each hypothesis
means every core hypothesis finds some symbol that explains the observation, so
the LLR magnitudes collapse and the LDPC loses the reliability information it
runs on. This is worth keeping as a general lesson: max-log is a fine
approximation for demapping a constellation and a bad one for marginalising a
nuisance parameter. It also nearly produced a wrong headline — the first run
reported "alphabet-aware is 1 dB worse", which would have been an
implementation choice reported as a fact about the channel. The control is what
caught it: the mechanism plainly worked on QPSK, so a *negative* result on
256QAM had to be an artefact rather than physics.

**3. A better number for the antenna, measured instead of inferred.** M14
computed the requirement indirectly — AWGN threshold 2.5 dB plus an
interference-as-variance model — and got 3.91 dB. M15 measures the composite
channel end to end: the threshold is **6.80 dB of channel SNR**, against the
**3.37 dB** measured on air, so the deficit is **3.43 dB**. The two agree to
within half a dB, and the gap is in the expected direction and for the expected
reason: a finite-alphabet interferer is slightly *less* harmful than Gaussian
noise of the same power. **3.43 dB is the number to shop against.**

## Addendum — what live 1080p would cost, measured (M13 follow-up)

M13 called the remaining work for the clear ROUTE service "plumbing, not
physics". True, but the plumbing has a speed wall, and it is worth measuring
before anyone assumes `watch` can simply be pointed at PLP 1: **the 64800 path
is pure NumPy, which is why a PLP 1 Frame currently takes ~33 s of wall clock
for 247 ms of air.**

**The LDPC is not the problem — it already fits.** `m9_gpu.GpuMinSum` is
generic over (checks, Ninner); only `FrameDecoder` hardcodes 16200. Built
against PLP 1's own chain (256QAM 11/15, Ninner 64800) and run on a full
Frame's 117 Blocks:

| operating point | CPU, 117 Blocks | GPU, 117 Blocks | vs the 247.11 ms Frame |
|---|---|---|---|
| 6 LDPC iterations | 11.0 s | **0.22 s** | 0.89x — fits |
| 3 LDPC iterations | 7.1 s | **0.06 s** | 0.25x — 4x headroom |

(A non-converging worst case at 50 iterations is 1.19 s, 4.8x over budget, so
the margin depends on how hard the air is working the decoder.)

**The demapper becomes the new bottleneck, and by a lot.** PLP 1 is 947,700
cells against a 256-point constellation where PLP 0 is 199,800 against 64
points — **27x the work**, measured. Scaling M11's live figure for PLP 0's GPU
BICM (22.8 ms) gives **~619 ms per Frame for PLP 1, 2.5x over budget on its
own.**

So the honest statement is: **the LDPC fits, the demapper does not, and the
total is unmeasured until both are built.** The standard remedy is candidate
pruning — find the nearest constellation point and evaluate only its
neighbourhood instead of all 256 — but the textbook per-axis decomposition does
NOT apply here, because a NUC constellation is not separable into independent I
and Q. Pruning would have to be built and gated against the exhaustive
demapper's bit-identical output, which is a real task and not a tweak.

Nothing here is a reason PLP 1 cannot be watched live; it is the measurement
that says which stage to work on, and it is not the one that looked obvious.

## Addendum — the demapper, and two wrong guesses about why it was slow

The bottleneck for live 1080p was measured as the demapper (~619 ms against a
247.11 ms Frame). `lab/m16_demap_fast.py` attacks it, and the useful part of
the result is which attacks failed.

| hypothesis | outcome |
|---|---|
| the arithmetic is the cost — make it a GEMM | **1.2x**, nearly nothing |
| the per-bit minima's fancy-index COPIES are the cost (16 copies of a ~1 GB array) | **SLOWER** — strided views over a reshape made it worse |
| it is simply single-threaded | **10.8x**, and that is the whole answer |

Cache-blocking settles it: dropping the chunk from 947,700 cells to 8,000
changed the time by 1%, which rules out a memory-bandwidth story outright. Then
12 threads took 15.3 s to **1.70 s**.

    exhaustive, 1 thread   18.4 s   74.4x Frame
    GEMM f32,   1 thread   15.3 s   62.1x
    GEMM f32,  12 threads   1.70 s   6.9x

This is M9's lesson for the third time — **the parallelism has to be inside the
stage** — and it is now the third stage where the obvious structural
explanation was wrong and threading was the answer.

The GEMM is kept because it is exact and free, not because it mattered: the
identity `|y-p|² = |y|² - 2Re(y·conj(p)) + |p|²` drops a term that is common to
every constellation point, and max-log LLR is a *difference* of two such
distances, so it cancels **algebraically**. One trap avoided: `demap_llr`
estimates sigma2 from absolute minima when the caller supplies none, and that
estimate does *not* survive the transform, so `|y|²` is added back before it.

Gated three ways: LLRs against the exhaustive demapper (relative delta 7e-16),
the same for the threaded path, and — the one that matters — **the decoded
codewords are identical**.

Still 6.9x over budget on CPU, so this does not make PLP 1 live on its own. The
GPU BICM path remains the candidate, and its true scaling is UNMEASURED: the
~619 ms figure is a linear extrapolation from PLP 0's 22.8 ms, which may well
be pessimistic if that measurement was launch-latency dominated rather than
throughput dominated.

## Addendum — the GPU demapper measured, and my extrapolation was 12.5x wrong

The previous addendum put PLP 1's GPU demap at ~619 ms, 2.5x over the Frame
budget, and flagged it as a LINEAR extrapolation from PLP 0's 22.8 ms that
"may well be pessimistic". It was pessimistic by **12.5x**. Measured, same
kernel `m9_gpu.GpuBicm.demap` uses, 5 repetitions after a warm-up:

| shape | precision | median | vs 247.11 ms Frame |
|---|---|---|---|
| PLP 0, 74 x 2700 x 64 | float64 | 4.27 ms | 0.017x |
| **PLP 1, 117 x 8100 x 256** | float64 | **49.4 ms** | **0.20x** |
| PLP 1, 117 x 8100 x 256 | float32 | 24.8 ms | 0.10x |

PLP 1 is 27x the *work* of PLP 0 but only **11.6x the time** — the GPU absorbs
the extra parallelism, which is exactly what a linear extrapolation from a
small kernel cannot capture. (The first run of this benchmark reported PLP 0 at
0.0 ms, which was a warm-up artefact; repetitions after an explicit warm-up and
`cuda.synchronize` give 4.27 ms with a 0.01 ms spread.)

**So the demapper is not the wall after all.** Assembling what is now measured,
for decoding PLP 1 alone:

    front end        69.5 ms   MEASURED live (M11) -- covers the whole Frame,
                               so it does not grow for Subframe 1
    cell pool       ~67 ms     ESTIMATED, not measured: 75 symbols at 16K
                               against PLP 0's 35 at 8K, scaled from 15.6 ms
    GPU demap        49.4 ms   MEASURED here
    GPU LDPC      60-220 ms    MEASURED (3 to 6 iterations)
    BCH etc         ~10 ms     MEASURED for PLP 0
    ------------------------------------------------------------------
    total        ~255-415 ms   against a 247.11 ms budget

which is **parity at the good end and 1.7x over at the bad end** — a tuning
problem, not a rebuild. The float32 demap would take another 25 ms off, but it
is a numerical change and would need the usual bit-identity gate before it
counts.

The honest status of each number matters more than the total: everything above
is measured except the cell pool, and that one is the largest remaining
unknown. It is also the obvious next measurement.

## Addendum — the AC-4 gap, bounded: our audio is correct, the decoder does not exist here

"No audio" has been carried since M7 as a tooling gap. It is now measured on
both sides rather than assumed.

**Our audio output is structurally correct and complete.** Parsing the audio
track's own boxes:

    mdhd version 1, timescale 90000 Hz
    116 fragments, 6960 AC-4 frames, 10,450,440 ticks = 116.12 s
    1501.5 ticks per frame = 16.683 ms
    a 59.94 Hz video frame is 16.683 ms  ->  AUDIO IS FRAME-LOCKED TO VIDEO

6960 frames at exactly one video-frame period each, against 114.15 s of video
on the same capture. `ffprobe` independently identifies the stream as `ac4`,
tag `ac-4`, 48 kHz, 6 channels. So the receiver is not dropping, mis-framing or
truncating anything — what comes out is a valid AC-4 elementary stream that is
frame-locked to the picture.

**No decoder is available on this machine, and that is now tested rather than
assumed:**

* `ffmpeg 8.1-full` — demuxes and identifies AC-4, has **no AC-4 decoder**
  (`-decoders | grep ac4` is empty). It knows what the stream is and cannot
  turn it into samples.
* `VLC 3.0.23` — worse: `unknown box type dac4`, no packetizer match, and a
  zero-byte output file. It cannot even parse the sample entry.

So the gap is real, external, and precisely located: **between a valid AC-4
elementary stream and a decoder that does not ship in any free tool installed
here.** Writing one is not a next step — AC-4 (ETSI TS 103 190) is a codec on
the scale of the physical layer this whole project just built, and it would be
its own campaign. The realistic paths are a player that gains AC-4 support, or
hardware that already has a Dolby decoder.

The useful part of this is what it rules OUT: no amount of work on the receiver
will produce sound, because the receiver is already producing exactly what it
should.

## Addendum — dolbyTuna: scoping an AC-4 decoder from what is actually on the air

Before writing a codec, read the config. The `dac4` box in our own captured
audio (36 bytes) parses as `ac4_dsi_v1`:

    ac4_dsi_version    1
    bitstream_version  2          (TS 103 190-2, the modern bitstream)
    fs_index           1  ->  48000 Hz
    frame_rate_index   3
    n_presentations    1

Two of those are cross-checked against measurements taken independently of this
parse: 48 kHz agrees with `ffprobe`, and `n_presentations = 1` is consistent
with a single audio track. **`frame_rate_index = 3` does NOT agree** with the
16.683 ms (59.94 Hz) frame period measured from the `mdhd` timing — my table
maps index 3 to 29.97 Hz, exactly half. The likely reconciliation is AC-4's
sub-frame multiplier (`b_sf_multiplier`), which lets the audio run at 2x or 4x
the signalled base rate, and 29.97 x 2 = 59.94 exactly. **That is a plausible
reading and not a verified one** — it is precisely the kind of thing the spec
has to settle, and it is recorded as open rather than smoothed over.

**What makes this project possible at all:** ETSI publishes TS 103 190-1 and
-2 free of charge. The bitstream is documented. That is not true of every
codec, and it is the single biggest enabler here.

**What makes it a real project:** there is no open-source AC-4 decoder.
FFmpeg 8.1 demuxes and identifies AC-4 and stops; VLC 3.0.23 cannot parse the
sample entry. So this would be built from the spec, not ported.

**What the config buys us in scope:** `n_presentations = 1` removes the entire
presentation/alternate-audio layer. One bitstream version, one sample rate.
The remaining core is an MDCT front end with companding, plus whatever spectral
extension and coupling tools the substream actually signals — comparable in
size to AAC with SBR and PS. Months of evenings, not a weekend, and not
comparable to a tweak.

**The cheap decisive next step, before any decoder work:** walk the raw AC-4
frames and parse each frame's TOC, then check the frame sizes chain exactly
against the `trun` sample sizes over all 6,960 frames. That is a hard,
falsifiable gate on whether we can read the bitstream at all — and if it fails,
no decoder is worth starting. It needs no new DSP and no radio.

## Addendum — dolbyTuna: what the AC-4 header gives up without the spec

Two experiments, one negative, and a structural finding neither was looking for.

**The frame length is not carried at a fixed offset.** Every bit position
0..259 at widths 8..20, searched for a value with a constant relationship to
the frame's true length (which the container supplies). 3480 frames, 384
distinct lengths spanning 617..1424 bytes — the search had ample power and
found nothing. **The substream sizes are variable-length coded.** Recorded so
the search is not repeated.

**The per-bit entropy map identifies fields by their statistics.** Over 3480
frames, every field parsed from the TOC lands where the entropy says it should:
`bitstream_version` H=0.000, `sequence_counter` H=0.999, `br_code` H=1.000
(so `wait_frames > 0`), `fs_index` and `frame_rate_index` H=0.000. The sharpest
is `b_iframe_global` at bit 23: **measured H = 0.3534, and the entropy
predicted from the independently counted 232 I-frames in 3480 is 0.3534.**

**And then the finding: FOUR bits carry the I-frame flag.**

    bits identical to b_iframe_global on all 3480 frames: 23, 95, 108, 147
    gaps between them:                                    72, 13, 39

Not "correlated" — *identical*, on every one of 3480 frames. Together with the
constant-bit runs (24 bits at 24..47, 38 bits at 109..146, 23 bits at 165..187,
last constant bit at 187), this says the header has a **FIXED-LAYOUT PREFIX
through at least bit 187**, and that something structural repeats four times
inside it. Four substreams, or four substream groups each carrying its own
random-access flag, is the obvious reading — and it is a hypothesis for the
spec to confirm, not a conclusion.

That is not a contradiction of the negative result above: the header prefix is
fixed, and it is the substream SIZES, later, that are variably coded.

**Why this stops here.** The remaining syntax is recall rather than reference,
and recall is what this project has been wrong about repeatedly. ETSI publishes
TS 103 190-1 and -2 free of charge but returns HTTP 403 to automated fetching,
and working around an access control is not on the table. The spec has to be
fetched by hand. When it is, this map inverts its role: it stops substituting
for the spec and becomes the **test of whether the spec has been read
correctly** — a parse that does not put four identical I-frame bits at 23, 95,
108 and 147, and H=0.3534 at bit 23, is a misreading.
