# ATSC 3.0 → television : architecture

Built 8/05–8/07. Two charts: the whole journey from photons to a watchable
file, then the AC-4 audio decoder on its own, because that is where most of
the new work went.

Solid arrows carry data. Dashed arrows carry control, timing or gates.
Boxes marked **NEW** landed 8/07.

> **This file is the 8/07 "first television" snapshot**, kept because it is the
> clearest single picture of the offline path that produced `rf33_tv.mp4`.
> For the system as it stands now — live chain, supervision, LDM, clocks — see
> [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md); for the standard itself, see
> [docs/HOW-ATSC3-WORKS.md](docs/HOW-ATSC3-WORKS.md). Claims below that time
> has overtaken are corrected in place and marked.

---

## 1. RF to television

```mermaid
flowchart TB
%% ============================================================================
%%  ATSC 3.0 (NextGen TV) -- from antenna to picture, sound and captions.
%%
%%  Read top-to-bottom as the signal's journey.  The left spine is the
%%  receiver proper (bootstrap up to MPUs); the four columns below it are the
%%  four assets the multiplex actually carries -- we worked with two of them
%%  for a day before noticing the other two existed.
%%
%%  The thing that makes it TELEVISION rather than two files is the MPU
%%  sequence number: it is the ONLY common clock, because each MPU is an
%%  independent ISOBMFF file whose own tfdt restarts at zero.
%% ============================================================================

    subgraph RX["Receiver -- A/321 up to MMTP  (8/05-8/06)"]
        direction TB
        ANT>"Antenna<br/>UHF yagi"]
        BOOT["<b>bootstrap</b><br/><i>atsc3/bootstrap.py</i><br/>6.144 MHz fixed &middot; ZC root 137<br/>gated 32/32 to &minus;18 dB"]
        L1["<b>L1 signalling</b><br/><i>m8_l1.py</i><br/>L1-Basic &middot; L1-Detail &middot; CRC"]
        FEC["<b>LDPC + BCH</b><br/><i>m9, m10, m15, m16</i><br/>77,547 of 77,552 blocks converged"]
        IP["<b>ALP &rarr; IP &rarr; MMTP</b><br/><i>m11_stream.py</i><br/>byte-identical to the batch chain"]
        MPU["<b>MPU reassembly</b><br/><i>m7_objects.py</i><br/>FT=0 init &middot; FT=1 moof &middot; FT=2 mfu"]
        ANT --> BOOT --> L1 --> FEC --> IP --> MPU
    end

    subgraph ASSETS["What the multiplex actually carries -- four assets, not three"]
        direction LR
        V["<b>pid12 &middot; hvc1</b><br/>video 1280&times;720 @ 60000/1001<br/>57 MPUs &times; 120 frames<br/><i>2 MPUs lost off air</i>"]
        A1["<b>pid13 &middot; ac-4</b><br/>audio 5.1<br/>58 MPUs &times; 60 frames"]
        A2["<b>pid14 &middot; ac-4</b><br/>audio #2 = SPANISH<br/><i>channel_pair &mdash; decoded since E35</i>"]
        S["<b>pid15 &middot; stpp</b><br/>IMSC1 / TTML captions<br/>59 MPUs &middot; <i>zero losses</i>"]
    end

    MPU --> V & A1 & A2 & S

    AC4["<b>AC-4 decoder</b><br/><i>m19 &hellip; m37</i> &mdash; see chart 2<br/>all six 5.1 channels + A-SPX high band"]
    RS["<b>output resample</b> <b>NEW</b><br/><i>m38_resample.py</i><br/>1001/960 per Table 83<br/><i>without it: 4.27 % fast</i>"]
    CAP["<b>caption builder</b> <b>NEW</b><br/><i>m40_captions.py</i><br/>986 raw cues &rarr; 69 readable<br/>92.7 % were roll-up re-typing"]
    MUX["<b>mux on the MPU grid</b> <b>NEW</b><br/><i>m39_mux.py</i><br/>audio and video agree to <b>0.00 ms</b><br/>HEVC copied, never re-encoded"]
    OUT(["<b>rf33_tv.mp4</b><br/>114.114 s &middot; 720p60 &middot; AC-4 &middot; captions"])

    A1 --> AC4 --> RS --> MUX
    S --> CAP --> MUX
    V --> MUX --> OUT

    CLK{{"<b>mpu_sequence_number</b><br/>the only common clock<br/>1 MPU = 120 video frames<br/>= 60 AC-4 frames = <b>2.002 s</b><br/><i>every tfdt is 0</i>"}}
    CLK -.-> MUX
    CLK -.-> CAP
    CLK -.-> RS

    HOLE["slot 225641782 lost<br/><i>drop the audio, do not invent video</i>"]
    HOLE -.-> MUX

    ATLAS["<b>radio-grid-atlas/atsc-3.0</b> <b>NEW</b><br/><i>measure.py</i> &middot; PUBLIC<br/>frame period 247.111 ms<br/>0.000 ms spread / 11 gaps"]
    BOOT -.-> ATLAS

    classDef new fill:#0b5,stroke:#063,color:#fff
    classDef gap fill:#a33,stroke:#611,color:#fff
    classDef clock fill:#048,stroke:#024,color:#fff
    class RS,CAP,MUX,ATLAS new
    class HOLE gap
    class CLK clock
```

---

## 2. The AC-4 decoder

```mermaid
flowchart TB
%% ============================================================================
%%  AC-4 audio, built from ETSI TS 103 190-1/-2 with no reference decoder.
%%
%%  Left column is the CORE coder (MDCT up to the crossover).  Right column
%%  is A-SPX, which rebuilds everything above it from a handful of envelope
%%  values.  They meet in the QMF domain at the bottom.
%%
%%  Two bugs are called out because both PASSED gates that looked green:
%%  the balance bug (three gates, all measuring the left channel) and the
%%  missing resampler (three subsystems agreeing on a shared wrong constant).
%% ============================================================================

    BITS["<b>ac4_substream</b><br/><i>m18, m19, m42</i><br/>audio_size &middot; b_more_bits<br/>closes 3480/3480 exactly<br/><b>+ payload_base</b> (E77) &mdash; RF33 never<br/>sets it, Fox always does"]

    subgraph CORE["Core coder -- up to the crossover"]
        direction TB
        FRAME["<b>framing</b><br/><i>m28_channels.py</i><br/>4 interval classes &middot; window groups<br/><i>short frames: per-group loops</i><br/><b>codec_mode 4 = ASPX_ACPL_3</b> (E77):<br/>a TWO-channel core + parametric coupling"]
        SFB["<b>band tables</b><br/><i>m27_sfb.py</i><br/>Annex B &middot; num_sfb(1536)=55<br/><i>tables NEST -- map explicitly</i>"]
        SPEC["<b>spectral data</b><br/><i>m29_audio.py</i><br/>dequant &middot; ungroup &middot; M/S &middot; noise fill"]
        FB["<b>filterbank</b><br/><i>m30_filterbank.py</i><br/>KBD windows &middot; block switching<br/>S&middot;S<sup>T</sup> = I to 4.4e&minus;16"]
        FRAME --> SPEC
        SFB -.-> SPEC
        SPEC --> FB
    end

    subgraph ASPX["A-SPX -- everything above the crossover"]
        direction TB
        APARSE["<b>side info</b><br/><i>m32_aspx_parse.py</i><br/>closes to audio_size, 3480/3480"]
        ABAND["<b>band layout</b><br/><i>m31_aspx_bands.py</i><br/>Pseudocode 67-71 &middot; 2 patches"]
        AENV["<b>envelope</b><br/><i>m34_aspx_env.py</i><br/>adjacent r = +0.76 vs +0.02 shuffled"]
        AHF["<b>HF patch</b><br/><i>m35_hfgen.py</i><br/>copy only &middot; writes no audio by design"]
        AADJ["<b>envelope adjust + noise</b><br/><i>m36_envadj.py</i><br/>Pseudocode 82/83/<b>84</b>/90-101"]
        APARSE --> AENV --> AADJ
        ABAND -.-> AHF --> AADJ
    end

    QMF["<b>64-band complex QMF</b><br/><i>m33_qmf.py</i><br/>78.49 dB reconstruction<br/>group delay 577 samples"]
    REND["<b>joint render</b><br/><i>m37_render_hf.py</i><br/>L and R together &middot; TS_OFFSET = 4"]
    PCM(["<b>PCM</b><br/>0 &hellip; 20.1 kHz physical"])

    BITS --> FRAME
    BITS --> APARSE
    FB --> QMF --> AADJ --> REND --> PCM

    BAL["<b>balance bug</b><br/>aspx_balance=1 &rArr; channel B is a<br/>PAN RATIO, not a level<br/><i>R was 20&times; too loud</i><br/>gates 1-3 all measured L"]
    RSMP["<b>missing resampler</b><br/>1536 is the MDCT length,<br/>not the frame length<br/><i>audio ran 4.27 % fast</i>"]
    BAL -.-> AADJ
    RSMP -.-> PCM

    classDef bug fill:#a33,stroke:#611,color:#fff
    classDef out fill:#0b5,stroke:#063,color:#fff
    class BAL,RSMP bug
    class PCM out
```

---

## What the charts do not show

Counted, not hidden. Rows struck through were true on 8/07 and have since been
closed; they are kept rather than deleted so the record stays legible.

| gap | size |
|---|---|
| additional harmonics | &lt;1.5 % of channel-frames &mdash; `SineTable` is not in the spec document we have |
| limiter tie-breaking | Pseudocode 72, partial |
| A-CPL | ~~not signalled on this service~~ &mdash; **superseded (E77)**: Fox signals `5_X_codec_mode 4` = ASPX_ACPL_3. We now parse the two-channel core and render it as correct stereo; `acpl_data_2ch` and a QMF-domain A-CPL synthesis stage are still missing, so the 5.1 upmix of *that* stream is not reachable yet |
| companding | signalled, not applied |
| `sap_mode 3` full SAP | 91 frames |
| SSF speech frontend | never used by this stream |
| ~~**pid14, the second audio**~~ | **CLOSED (E35)**: it is the Spanish programme, a `channel_pair` element, decoded and offered as a selectable audio track in the live viewer |
| ~~**32.1**~~ | **CLOSED (M13)**: subframe 1 / PLP 1 is demodulated and gated — 468/468 FEC blocks BCH-zero over 4 frames, later 1872/1872 over 16. It carries five ROUTE/LCT flows and four additional 1080p60 services. Live playback of PLP 1 is still blocked on real-time cost, not on the physical layer |

Two further limits worth naming, both measured after this document was
written:

| limit | number |
|---|---|
| full 5.1 AC-4 in real time | discrete-verified on all six channels, but **0.62x** even threaded — recordings only. Live runs the stereo pair at 1.27x |
| PLP 1 live | the 64800 path is the cost; the demapper, not the LDPC, becomes the bottleneck |
