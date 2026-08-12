"""ATSC 3.0 bootstrap synthesizer + detector.

Task #37 / milestone M1.  The bootstrap (ATSC A/321) is the fixed, universal,
deliberately-discoverable preamble that begins every ATSC 3.0 physical layer
frame.  It is transmitted identically no matter how the rest of the frame is
configured, at a FIXED 6.144 Msample/s and 4.5 MHz bandwidth regardless of the
channel bandwidth.  That makes it the one part of ATSC 3.0 you can synthesize
bit-exactly and validate against, with no radio.

Spec source
-----------
ATSC A/321:2024-04 "System Discovery and Signaling", Sections 4-6 + Annex A/B.
https://www.atsc.org/wp-content/uploads/2024/04/A321-2024-04-System-Discovery-and-Signaling.pdf

Everything in the CONSTANTS block below is quoted from that document.  Items
that required interpretation are collected in ASSUMPTIONS at the bottom of this
docstring and are individually tagged ``# ASSUMPTION`` at their use site.

Two independent detectors
-------------------------
This module deliberately ships TWO detectors, because a detector built only
from a synthesized reference can be self-consistently wrong (right code, wrong
spec => passes its own selftest, fails on real air, and you cannot tell which).

  1. ``detect_structural`` -- ASSUMPTION-FREE.  Keys only on the A/B/C time
     geometry: every bootstrap symbol contains a 520-sample run that repeats
     exactly 2048 samples later, and those runs recur every 3072 samples.  It
     does not know or care about the Zadoff-Chu root, the PN seed, or the
     signaling bits.  If this fires on real air, the carrier IS an ATSC 3.0
     bootstrap (or something structurally identical).

  2. ``detect_matched`` -- uses the synthesized reference.  Gives cyclic
     shifts, hence the signaling bits, hence the actual decoded fields.
     Depends on the ZC root and PN sequence being right.

Reading them together is the diagnostic:

    structural  matched   meaning
    ----------  -------   -------
    NO          NO        not ATSC 3.0 (or too weak / wrong center freq)
    YES         YES       ATSC 3.0, and our spec implementation is correct
    YES         NO        ATSC 3.0 bootstrap IS present but our ZC/PN
                          implementation is WRONG -- go fix ASSUMPTIONS
    NO          YES       should not happen; treat as a bug

ASSUMPTIONS (things A/321 states in prose or figures that required a choice)
---------------------------------------------------------------------------
A1  PN LFSR wiring.  A/321 Fig 5.2 draws registers r_{l-1}..r_0 with taps
    g_{l-1}..g_0 and says the register is "clocked one position to the right",
    with p(0) the output before any clocking.  We read that as: output = r_0,
    feedback into r_{l-1} = XOR of g_i*r_i, i.e. p(k+16) = p(k+15) ^ p(k+14) ^
    p(k+1) ^ p(k).  r_init = {r_15..r_0} so the MSB of the table value is r_15.
    RISK: MEDIUM.  The mirror convention (reciprocal polynomial) is equally
    drawable.  Mitigated by ``PN_VARIANTS`` + ``--sweep-pn`` which tries all
    plausible wirings against a real capture; and by detector 1 which does not
    use the PN at all.

A2  BCA part-B time offset.  The A/321 PDF renders the BCA part-B exponent in
    a symbol font whose digits are offset by +4 per digit in text extraction.
    The raw glyphs read "964"; Figure 5.6's caption reads "exp(-j2*pi*f_delta*
    (t - 520*Ts))".  520 is confirmed three ways: (a) the +4 digit shift maps
    964 -> 520, verified independently on the generator polynomial where
    "5:","59","58" map to 16,15,14; (b) the figure caption; (c) it is the only
    value that makes the closed form agree with A/321's own alternative
    definition (part B = one-subcarrier-down frequency shift of part A).
    RISK: LOW.

A3  IFFT scaling is 1/sqrt(N_ZC-1) per A/321 5.2.4.  Irrelevant to detection
    (all correlations here are normalized) but implemented for fidelity.

A4  Bootstrap symbol 0 is NOT phase inverted.  A/321 inverts only symbol
    n = N_S - 1.  With the mandated N_S >= 4 this is symbol 3.  RISK: NONE.

A5  Signaling bit packing within a symbol is MSB-first in field order as
    tabled (Tables 6.2/6.4/6.5).  RISK: LOW -- stated explicitly.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np

# --------------------------------------------------------------------------
# CONSTANTS -- all from ATSC A/321:2024-04
# --------------------------------------------------------------------------

FS = 6_144_000              # A/321 5.1: fixed bootstrap sampling rate, Hz
BW = 4_500_000              # A/321 5.1: fixed bootstrap bandwidth, Hz
N_FFT = 2048                # A/321 5.1
DF = FS / N_FFT             # 3000 Hz subcarrier spacing
N_ZC = 1499                 # A/321 5.2.1: largest prime fitting 4.5 MHz @ 3 kHz
N_H = (N_ZC - 1) // 2       # 749
N_A, N_B, N_C = 2048, 504, 520      # A/321 5.4
N_SYM = N_A + N_B + N_C             # 3072 samples = 500 us

ZC_ROOT_MAJOR0 = 137        # A/321 6.1: q for bootstrap_major_version = 0
                            # 0..136 and 138..1498 are Reserved for future majors

# A/321 5.2.2: g = {g_16,...,g_0}; p(x) = x^16 + x^15 + x^14 + x + 1
PN_ORDER = 16
PN_G = (1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1)   # g_16 .. g_0

# A/321 Table 6.1: r_init per bootstrap_minor_version (major version 0)
MINOR_SEEDS = {
    0: 0x019D, 1: 0x00ED, 2: 0x01E8, 3: 0x00E8,
    4: 0x00FB, 5: 0x0021, 6: 0x0054, 7: 0x00EC,
}

# A/321 Table 6.3 -- min_time_to_next, index -> milliseconds (31 = N/A)
def min_time_to_next_ms(x: int):
    """A/321 6.1.1.1 sliding-scale mapping.  Returns ms, or None for 31."""
    if not 0 <= x < 32:
        raise ValueError("min_time_to_next must be 0..31")
    if x == 31:
        return None                       # spec: "shall not be indicated"
    if x < 8:
        return 50 * x + 50
    if x < 16:
        return 100 * (x - 8) + 500
    if x < 24:
        return 200 * (x - 16) + 1300
    return 400 * (x - 24) + 2900


SYSTEM_BANDWIDTH = {0: "6 MHz", 1: "7 MHz", 2: "8 MHz", 3: ">8 MHz (reserved)"}


def bsr_sample_rate_hz(n: int) -> float:
    """A/321 6.1: post-bootstrap sample rate = (N + 16) * 0.384 MHz."""
    return (n + 16) * 384_000.0


# --------------------------------------------------------------------------
# PN sequence (A/321 5.2.2)
# --------------------------------------------------------------------------

# Variant wirings.  'spec' is our reading of Figure 5.2 (see ASSUMPTION A1).
# The others are the plausible mirror readings, kept so that a real capture can
# arbitrate rather than us guessing.  Each entry is
#   (taps, seed_msb_is_r_high, description)
# taps = tuple of register indices XORed to form the feedback bit.
_TAPS_SPEC = tuple(i for i in range(PN_ORDER) if PN_G[PN_ORDER - i])   # g_i, i<16
_TAPS_RECIP = tuple(PN_ORDER - 1 - i for i in _TAPS_SPEC)

PN_VARIANTS = {
    "spec":        (_TAPS_SPEC,  True,  "Fig 5.2 as read: out=r0, fb->r15, seed MSB=r15"),
    "spec_revseed": (_TAPS_SPEC, False, "same taps, seed bit order reversed"),
    "recip":       (_TAPS_RECIP, True,  "reciprocal polynomial wiring"),
    "recip_revseed": (_TAPS_RECIP, False, "reciprocal taps, seed reversed"),
}


@lru_cache(maxsize=64)
def _pn_cached(seed: int, length: int, variant: str) -> np.ndarray:
    taps, msb_high, _ = PN_VARIANTS[variant]
    # r[i] holds register r_i.  r_init = {r_{l-1},...,r_0}: MSB of seed -> r_15.
    if msb_high:
        r = [(seed >> i) & 1 for i in range(PN_ORDER)]          # r_i = bit i
    else:
        r = [(seed >> (PN_ORDER - 1 - i)) & 1 for i in range(PN_ORDER)]
    out = np.empty(length, dtype=np.uint8)
    for k in range(length):
        out[k] = r[0]                    # ASSUMPTION A1: output tapped at r_0
        fb = 0
        for t in taps:
            fb ^= r[t]
        r = r[1:] + [fb]                 # clock "one position to the right"
    out.setflags(write=False)
    return out


def pn_sequence(seed: int, length: int, variant: str = "spec") -> np.ndarray:
    """Generate p(k), k = 0..length-1, as uint8 in {0,1}.

    A/321 5.2.2: p(0) is the generator output after loading the seed and
    BEFORE any clocking; each subsequent clock emits the next bit.
    """
    return _pn_cached(int(seed), int(length), variant)


def pn_pm1(seed: int, length: int, variant: str = "spec") -> np.ndarray:
    """c(k) = 1 - 2*p(k), values in {+1,-1}.  A/321 5.2.3."""
    return 1.0 - 2.0 * pn_sequence(seed, length, variant).astype(np.float64)


# --------------------------------------------------------------------------
# Zadoff-Chu sequence (A/321 5.2.1)
# --------------------------------------------------------------------------

def zc_sequence(q: int = ZC_ROOT_MAJOR0, n_zc: int = N_ZC) -> np.ndarray:
    """z_q(k) = exp(-j*pi*q*k*(k+1)/N_ZC), k = 0..N_ZC-1."""
    if not 1 <= q <= n_zc - 1:
        raise ValueError("ZC root q must be in 1..N_ZC-1")
    k = np.arange(n_zc, dtype=np.float64)
    # keep the phase argument small before the exponential to avoid float blowup
    ph = np.mod(q * k * (k + 1.0), 2.0 * n_zc) / n_zc
    return np.exp(-1j * np.pi * ph)


# --------------------------------------------------------------------------
# Gray coded signaling bits <-> relative cyclic shift (A/321 5.3.2, Annex A)
# --------------------------------------------------------------------------

MAX_SIG_BITS = 11           # log2(N_FFT)


def bits_to_relative_shift(bits) -> int:
    """A/321 5.3.2 Gray map.  ``bits`` is b_0..b_{Nb-1}, MSB (most robust) first."""
    bits = [int(b) & 1 for b in bits]
    nb = len(bits)
    if not 1 <= nb <= MAX_SIG_BITS:
        raise ValueError("1..11 signaling bits")
    pivot = MAX_SIG_BITS - 1 - nb        # index i where m_i is forced to 1
    d = 0
    for i in range(MAX_SIG_BITS):
        if i > pivot:
            m = 0
            for k in range(0, MAX_SIG_BITS - 1 - i + 1):
                m ^= bits[k]
        elif i == pivot:
            m = 1
        else:
            m = 0
        d |= m << i
    return d


def relative_shift_to_bits(d: int, n_bits: int) -> list:
    """A/321 Annex A Gray de-map.  Returns b_0..b_{n_bits-1}, MSB first."""
    m = [(d >> i) & 1 for i in range(MAX_SIG_BITS)]
    out = []
    for i in range(n_bits):
        out.append(m[MAX_SIG_BITS - 1] if i == 0
                   else m[MAX_SIG_BITS - i] ^ m[MAX_SIG_BITS - 1 - i])
    return out


# --------------------------------------------------------------------------
# Symbol synthesis
# --------------------------------------------------------------------------

def freq_domain_symbol(n: int, pn_bits: np.ndarray, q: int = ZC_ROOT_MAJOR0,
                       invert: bool = False) -> np.ndarray:
    """S'_n(k) mapped onto N_FFT bins (numpy FFT bin order).  A/321 5.2.3."""
    z = zc_sequence(q)
    c = 1.0 - 2.0 * pn_bits.astype(np.float64)
    s = np.zeros(N_FFT, dtype=np.complex128)
    k_neg = np.arange(-N_H, 0)
    k_pos = np.arange(1, N_H + 1)
    s[k_neg % N_FFT] = z[k_neg + N_H] * c[(n + 1) * N_H + k_neg]
    s[k_pos] = z[k_pos + N_H] * c[(n + 1) * N_H - k_pos]
    s[0] = 0.0                            # DC subcarrier null (A/321 5.2.3)
    if invert:                            # A/321 5.2.3 termination signaling
        s = -s
    return s


def time_domain_A(s: np.ndarray, shift: int) -> np.ndarray:
    """IFFT then apply the absolute cyclic shift: A_n(t) = A~_n((t+M_n) mod N)."""
    a_tilde = np.fft.ifft(s) * N_FFT / np.sqrt(N_ZC - 1)     # ASSUMPTION A3
    return np.roll(a_tilde, -int(shift) % N_FFT)


def assemble_cab(a: np.ndarray) -> np.ndarray:
    """A/321 5.4.1 -- symbol 0 (sync symbol)."""
    t = np.arange(2568, 3072)
    out = np.empty(N_SYM, dtype=np.complex128)
    out[0:520] = a[1528:2048]                       # C = last 520 of A
    out[520:2568] = a[0:2048]                       # A
    # B = last 504 of A, +f_delta shift and e^{-j*pi} phase, both folded into
    # exp(j*2*pi*f_delta*t) with t the absolute sample index in the symbol.
    out[2568:3072] = a[t - 1024] * np.exp(1j * 2 * np.pi * t / N_FFT)
    return out


def assemble_bca(a: np.ndarray) -> np.ndarray:
    """A/321 5.4.2 -- symbols 1..N_S-1."""
    t = np.arange(0, 504)
    out = np.empty(N_SYM, dtype=np.complex128)
    # B = first 504 of C, -f_delta shift.  Offset 520*Ts -- see ASSUMPTION A2.
    out[0:504] = a[t + 1528] * np.exp(-1j * 2 * np.pi * (t - 520) / N_FFT)
    out[504:1024] = a[1528:2048]                    # C = last 520 of A
    out[1024:3072] = a[0:2048]                      # A
    return out


# Field layouts for bootstrap_major_version=0 / minor_version=0 (A/321 6.1.1.1)
SIG_BITS_PER_SYMBOL = 8


def pack_signaling(symbol: int, fields: dict) -> list:
    """Return the 8 signaling bits (MSB first) for symbol 1, 2 or 3."""
    def u(v, n):
        v = int(v)
        if not 0 <= v < (1 << n):
            raise ValueError(f"field out of range for {n} bits: {v}")
        return [(v >> (n - 1 - i)) & 1 for i in range(n)]

    if symbol == 1:
        return (u(fields.get("ea_wake_up_1", 0), 1)
                + u(fields.get("min_time_to_next", 10), 5)
                + u(fields.get("system_bandwidth", 0), 2))
    if symbol == 2:
        return (u(fields.get("ea_wake_up_2", 0), 1)
                + u(fields.get("bsr_coefficient", 2), 7))
    if symbol == 3:
        return u(fields.get("preamble_structure", 0), 8)
    raise ValueError("symbol must be 1..3")


def unpack_signaling(bits_by_symbol: dict) -> dict:
    """Inverse of pack_signaling, plus human-readable interpretation."""
    def val(bits):
        v = 0
        for b in bits:
            v = (v << 1) | (b & 1)
        return v

    out = {}
    b1 = bits_by_symbol.get(1)
    if b1 is not None:
        out["ea_wake_up_1"] = b1[0]
        out["min_time_to_next"] = val(b1[1:6])
        out["min_time_to_next_ms"] = min_time_to_next_ms(out["min_time_to_next"])
        out["system_bandwidth"] = val(b1[6:8])
        out["system_bandwidth_str"] = SYSTEM_BANDWIDTH[out["system_bandwidth"]]
    b2 = bits_by_symbol.get(2)
    if b2 is not None:
        out["ea_wake_up_2"] = b2[0]
        out["bsr_coefficient"] = val(b2[1:8])
        out["post_bootstrap_sample_rate_hz"] = bsr_sample_rate_hz(out["bsr_coefficient"])
    b3 = bits_by_symbol.get(3)
    if b3 is not None:
        out["preamble_structure"] = val(b3)
    if "ea_wake_up_1" in out and "ea_wake_up_2" in out:
        out["ea_wake_up"] = (out["ea_wake_up_1"] << 1) | out["ea_wake_up_2"]
    return out


def synthesize(minor_version: int = 0, n_symbols: int = 4,
               q: int = ZC_ROOT_MAJOR0, pn_variant: str = "spec",
               fields: dict | None = None, amplitude: float = 1.0):
    """Build a complete, spec-correct bootstrap.

    Returns (iq, info) where iq is complex128 at 6.144 Msps of length
    n_symbols*3072, and info carries the cyclic shifts and the fields used.
    """
    fields = dict(fields or {})
    seed = MINOR_SEEDS[minor_version]
    pn = pn_sequence(seed, (n_symbols + 1) * N_H, pn_variant)

    shifts, bits_by_symbol, parts = [], {}, []
    m_prev = 0
    for n in range(n_symbols):
        if n == 0:
            m = 0                                   # A/321 5.3.3
        else:
            bits = pack_signaling(n, fields) if n <= 3 else [0] * SIG_BITS_PER_SYMBOL
            bits_by_symbol[n] = bits
            d = bits_to_relative_shift(bits)
            m = (m_prev + d) % N_FFT
        shifts.append(m)
        m_prev = m
        invert = (n == n_symbols - 1)               # ASSUMPTION A4
        s = freq_domain_symbol(n, pn, q, invert)
        a = time_domain_A(s, m)
        parts.append(assemble_cab(a) if n == 0 else assemble_bca(a))

    iq = np.concatenate(parts)
    iq *= amplitude / np.sqrt(np.mean(np.abs(iq) ** 2))
    info = {
        "minor_version": minor_version, "major_version": 0, "zc_root": q,
        "pn_seed": seed, "pn_variant": pn_variant, "n_symbols": n_symbols,
        "abs_shifts": shifts, "bits": bits_by_symbol,
        "fields": unpack_signaling(bits_by_symbol),
    }
    return iq, info


def reference_A(n: int, minor_version: int = 0, q: int = ZC_ROOT_MAJOR0,
                pn_variant: str = "spec", n_symbols: int = 4) -> np.ndarray:
    """Un-shifted part A (A~_n) for symbol n -- the matched-filter template."""
    pn = pn_sequence(MINOR_SEEDS[minor_version], (n_symbols + 1) * N_H, pn_variant)
    s = freq_domain_symbol(n, pn, q, invert=(n == n_symbols - 1))
    return np.fft.ifft(s) * N_FFT / np.sqrt(N_ZC - 1)


# --------------------------------------------------------------------------
# Detector 1 -- structural, ASSUMPTION-FREE
# --------------------------------------------------------------------------

def _sliding_sum(x: np.ndarray, w: int) -> np.ndarray:
    """Sum of each length-w window; result length len(x)-w+1."""
    c = np.concatenate(([0.0 + 0.0j] if np.iscomplexobj(x) else [0.0],
                        np.cumsum(x)))
    return c[w:] - c[:-w]


def structural_metric(x: np.ndarray, win: int = N_C, lag: int = N_A):
    """Delay-and-correlate at lag 2048 over a 520-sample window.

    Every bootstrap symbol -- CAB or BCA -- contains a 520-sample run that is
    repeated verbatim 2048 samples later (part C is the last 520 samples of
    part A).  This keys on nothing but that geometry: no ZC root, no PN seed,
    no signaling bits.

    Returns (corr, energy).  ``corr`` is complex; its magnitude is the match
    strength and its angle carries the carrier frequency offset.  Callers fold
    corr COHERENTLY across the 4 symbols -- the correlation phase depends only
    on CFO and the fixed lag, not on the symbol index, so coherent folding is
    valid at any CFO and buys the full 4x processing gain.
    """
    n = len(x) - lag
    if n < win:
        return np.zeros(0, dtype=np.complex128), np.zeros(0)
    prod = x[:n] * np.conj(x[lag:lag + n])
    pw = np.abs(x[:n]) ** 2 + np.abs(x[lag:lag + n]) ** 2
    return _sliding_sum(prod, win), _sliding_sum(pw, win)


def run_offset(n: int) -> int:
    """Where inside symbol n the 520-sample lag-2048 repeat starts.

    CAB (symbol 0):  part C occupies [0,520) and is repeated at [2048,2568).
    BCA (symbols>0): part C occupies [504,1024) and is repeated at [2552,3072).
    Getting this right is what lets the 4 symbols be coherently folded.
    """
    return 0 if n == 0 else 504


# Default detection threshold on the peak-to-background RATIO of the folded
# structural metric.  Calibrated in selftest.py G6 against Gaussian noise and
# 8-VSB; see the printed max-ratio numbers there.  An absolute threshold on the
# normalized correlation cannot work -- that statistic saturates at SNR/(1+SNR)
# no matter how much processing gain you apply, which caps detection near 0 dB.
#
# Calibration (selftest G6, 100 trials x ~20k candidate positions each):
#   Gaussian noise    max peak_ratio 5.13
#   synthetic 8-VSB   max peak_ratio 2.94
#   real off-air 8-VSB (20 s captures, 123M positions)  max 3.29
# 6.0 would leave only a 17% margin over the noise maximum, and the maximum
# grows with the number of positions searched -- a 20 s capture searches 60x
# more than one selftest trial.  8.0 is the conservative choice, which is the
# right call for this detector's ROLE: it is the assumption-free corroborator,
# so it should speak only when certain.  The matched filter is both more
# sensitive AND has the healthier margin (noise max 3.41 vs threshold 8.0), and
# is the detector of record for sensitivity.
STRUCT_RATIO_THRESHOLD = 8.0


def detect_structural(x: np.ndarray, n_symbols: int = 4,
                      threshold: float = STRUCT_RATIO_THRESHOLD,
                      top_k: int = 12):
    """Coherently fold the lag-2048 correlation over the 3072-sample pitch.

    ``position`` is the index of the first sample of bootstrap symbol 0.
    ``score`` is the normalized correlation (0..1) at the peak; ``peak_ratio``
    is peak/background and is what the decision is made on.
    """
    corr, energy = structural_metric(x)
    span = (n_symbols - 1) * N_SYM + run_offset(n_symbols - 1)
    if len(corr) <= span:
        return {"detected": False, "score": 0.0, "reason": "capture too short",
                "position": None, "cfo_hz": None, "noise_floor": None,
                "peak_ratio": 0.0, "threshold": threshold}
    usable = len(corr) - span
    fcorr = np.zeros(usable, dtype=np.complex128)
    fener = np.zeros(usable)
    for j in range(n_symbols):
        o = j * N_SYM + run_offset(j)
        fcorr += corr[o: o + usable]           # COHERENT across symbols
        fener += energy[o: o + usable]
    folded = 2.0 * np.abs(fcorr) / np.maximum(fener, 1e-30)

    pos = int(np.argmax(folded))
    score = float(folded[pos])
    # Background estimate.  The exclusion zone must cover the WHOLE bootstrap
    # footprint (n_symbols * 3072), not just half a symbol: partial alignments
    # of the same bootstrap raise the folded metric for thousands of samples
    # either side of the true peak.  Excluding only +-N_SYM/2 lets the signal
    # contaminate its own noise floor, which compresses peak_ratio at HIGH SNR
    # and made the detector non-monotonic (worse at 0 dB than at -3 dB).
    guard = n_symbols * N_SYM
    mask = np.ones(len(folded), dtype=bool)
    mask[max(0, pos - guard): min(len(folded), pos + guard)] = False
    if mask.sum() >= 1000:
        floor = float(np.median(folded[mask]))
    else:
        # short block: not enough clean samples, fall back to a low percentile
        floor = float(np.percentile(folded, 10))
    ratio = score / max(floor, 1e-12)
    # CFO from the phase of the lag-2048 correlation
    cfo = float(-np.angle(fcorr[pos]) * FS / (2.0 * np.pi * N_A))

    # Top-K candidates.  Even when the peak does not clear threshold, the true
    # symbol-0 position is usually still ranked highly -- and the matched
    # filter (2048x4 coherent gain) can confirm it from a mere hint.  Handing
    # the matched stage only argmax throws that away and makes the whole
    # detector acquisition-limited.
    cands = []
    work = folded.copy()
    for _ in range(max(1, top_k)):
        i = int(np.argmax(work))
        if work[i] <= 0:
            break
        cands.append({"position": i, "score": float(folded[i]),
                      "peak_ratio": float(folded[i] / max(floor, 1e-12)),
                      "cfo_hz": float(-np.angle(fcorr[i]) * FS / (2.0 * np.pi * N_A))})
        work[max(0, i - 512): i + 512] = 0.0

    return {
        "detected": ratio >= threshold,
        "score": score,
        "noise_floor": floor,
        "peak_ratio": ratio,
        "position": pos,          # index of first sample of bootstrap symbol 0
        "cfo_hz": cfo,            # unambiguous only within +-1.5 kHz
        "threshold": threshold,
        "candidates": cands,
    }


def acquire_matched(x, cfo_grid=(-2000.0, -1000.0, 0.0, 1000.0, 2000.0),
                    minor_version: int = 0, q: int = ZC_ROOT_MAJOR0,
                    pn_variant: str = "spec", n_symbols: int = 4, top_k: int = 8):
    """DEEP acquisition: slide the symbol-0 template over the whole record.

    A/321 5.3.3 fixes bootstrap symbol 0 at absolute cyclic shift 0, so its
    part A is a fully known 2048-sample waveform.  Correlating it against every
    sample position gives 2048-sample coherent processing gain for ACQUISITION
    (not merely confirmation), which is what lets the detector find a bootstrap
    that the structural stage cannot see through co-channel interference.

    Returns a list of candidate dicts (position = start of symbol 0), best first.
    """
    from scipy.signal import fftconvolve

    x = np.asarray(x, dtype=np.complex64)
    ref = reference_A(0, minor_version, q, pn_variant, n_symbols).astype(np.complex64)
    ref = ref / np.sqrt(np.sum(np.abs(ref) ** 2))
    mf = np.conj(ref[::-1])
    # running energy over the 2048-sample window, for normalization
    e = _sliding_sum(np.abs(x.astype(np.complex128)) ** 2, N_FFT)

    best = None
    for f in cfo_grid:
        y = x if f == 0 else (x * np.exp(-2j * np.pi * f * np.arange(len(x))
                                         / FS).astype(np.complex64))
        c = fftconvolve(y, mf, mode="valid")          # len = len(x)-2047
        n = min(len(c), len(e))
        m = np.abs(c[:n]) / np.sqrt(np.maximum(e[:n], 1e-30))
        if best is None or m.max() > best[0].max():
            best = (m, f)
    m, f = best
    floor = float(np.median(m))
    cands = []
    work = m.copy()
    for _ in range(top_k):
        i = int(np.argmax(work))
        if work[i] <= 0:
            break
        cands.append({"position": int(i - 520),      # part A starts 520 into CAB
                      "score": float(m[i]),
                      "peak_ratio": float(m[i] / max(floor, 1e-12)),
                      "cfo_hz": f})
        work[max(0, i - 4096): i + 4096] = 0.0
    return [c for c in cands if c["position"] >= 0]


def detect(x: np.ndarray, fs_in: float | None = None, n_symbols: int = 4,
           minor_version: int = 0, q: int = ZC_ROOT_MAJOR0,
           pn_variant: str = "spec",
           struct_threshold: float = STRUCT_RATIO_THRESHOLD,
           matched_threshold: float = 8.0, integer_cfo_search: int = 10,
           deep: bool = False, max_candidates: int = 12) -> dict:
    """Run both detectors on a block of IQ and merge the verdicts.

    ``fs_in`` may be omitted if ``x`` is already at 6.144 Msps.  Returns a dict
    with ``verdict`` in {NO_BOOTSTRAP, ATSC3_CONFIRMED, ATSC3_STRUCTURE_ONLY,
    MATCHED_ONLY_SUSPECT}.
    """
    if fs_in is not None and abs(fs_in - FS) > 1.0:
        from .capture import to_bootstrap_rate
        x = to_bootstrap_rate(x, fs_in)
    x = np.asarray(x, dtype=np.complex128)

    st = detect_structural(x, n_symbols, struct_threshold)
    out = {"structural": st, "matched": None}
    if st.get("position") is None:
        out["verdict"] = "NO_BOOTSTRAP"
        return out

    cands = list(st.get("candidates") or [{"position": st["position"],
                                           "cfo_hz": st["cfo_hz"]}])
    if deep:
        # Deep candidates come from a 2048-sample coherent correlation and are
        # far better hypotheses than structural peaks -- they must go FIRST, or
        # max_candidates truncates them away and deep mode silently does nothing.
        cands = acquire_matched(x, minor_version=minor_version, q=q,
                                pn_variant=pn_variant, n_symbols=n_symbols) + cands
        max_candidates = max(max_candidates, 20)
    mt = None
    for c in cands[:max_candidates]:
        m = detect_matched(x, c["position"], c.get("cfo_hz") or 0.0,
                           minor_version, q, pn_variant, n_symbols,
                           integer_cfo_search)
        if not m or "mean_peak_ratio" not in m:
            continue
        m["detected"] = (m["mean_peak_ratio"] >= matched_threshold
                         and m.get("sym0_shift_ok", False))
        if mt is None or m["mean_peak_ratio"] > mt["mean_peak_ratio"]:
            mt = m
        if m["detected"]:
            break
    out["matched"] = mt
    if mt is None:
        out["verdict"] = "NO_BOOTSTRAP"
        return out

    s_ok, m_ok = st["detected"], mt["detected"]
    # "matched only" is a legitimate positive, not a bug: the structural metric
    # is normalized by total energy so a strong co-channel carrier suppresses
    # it, while the matched filter's coherent gain cuts through.
    out["verdict"] = ("ATSC3_CONFIRMED" if (s_ok and m_ok) else
                      "ATSC3_STRUCTURE_ONLY" if s_ok else
                      "ATSC3_MATCHED_ONLY" if m_ok else "NO_BOOTSTRAP")
    if s_ok and m_ok:
        out["fields"] = mt.get("fields", {})
        out["cfo_hz"] = mt.get("cfo_hz")
    return out


# --------------------------------------------------------------------------
# Detector 2 -- matched filter against the synthesized reference
# --------------------------------------------------------------------------

def _cyclic_xcorr(rx: np.ndarray, ref: np.ndarray) -> np.ndarray:
    """Circular cross-correlation of two length-N_FFT blocks, via FFT."""
    return np.fft.ifft(np.fft.fft(rx) * np.conj(np.fft.fft(ref)))


def _symbol_A_offset(n: int) -> int:
    """Offset of part A within a symbol: CAB -> 520, BCA -> 1024."""
    return 520 if n == 0 else 1024


def detect_matched(x: np.ndarray, sym0_start: int, cfo_hz: float = 0.0,
                   minor_version: int = 0, q: int = ZC_ROOT_MAJOR0,
                   pn_variant: str = "spec", n_symbols: int = 4,
                   integer_cfo_search: int = 10, sym0_shift_tol: int = 64):
    """Correlate each symbol's part A against the reference.

    ``sym0_start`` is the index of the first sample of bootstrap symbol 0
    (start of part C).  ``cfo_hz`` is the fractional estimate from the
    structural detector; integer multiples of the 3 kHz subcarrier spacing are
    resolved by searching +-``integer_cfo_search`` bins.

    Returns a dict with per-symbol peak ratios, absolute cyclic shifts, the
    decoded signaling bits and the interpreted fields.
    """
    need = sym0_start + n_symbols * N_SYM
    if sym0_start < 0 or need > len(x):
        return {"detected": False, "reason": "symbol window outside capture"}

    refs = [reference_A(n, minor_version, q, pn_variant, n_symbols)
            for n in range(n_symbols)]

    best = None
    for k in range(-integer_cfo_search, integer_cfo_search + 1):
        f = cfo_hz + k * DF
        t = np.arange(need - sym0_start)
        seg = x[sym0_start:need] * np.exp(-1j * 2 * np.pi * f * t / FS)
        peaks, shifts, ratios, phases = [], [], [], []
        for n in range(n_symbols):
            off = n * N_SYM + _symbol_A_offset(n)
            blk = seg[off:off + N_FFT]
            r = _cyclic_xcorr(blk, refs[n])
            mag = np.abs(r)
            idx = int(np.argmax(mag))
            pk = float(mag[idx])
            # noise floor of the correlation: RMS excluding a guard around peak
            g = np.ones(N_FFT, dtype=bool)
            g[(np.arange(-8, 9) + idx) % N_FFT] = False
            fl = float(np.sqrt(np.mean(mag[g] ** 2)))
            peaks.append(pk)
            ratios.append(pk / max(fl, 1e-30))
            # ifft(fft(rx)*conj(fft(ref)))[m] = sum_t rx[t]*conj(ref[t-m]).
            # A_n(t) = A~_n((t+M) mod N), so the peak lands at m = -M mod N.
            shifts.append((-idx) % N_FFT)
            phases.append(float(np.angle(r[idx])))
        mean_ratio = float(np.mean(ratios))
        if best is None or mean_ratio > best["mean_peak_ratio"]:
            best = {"cfo_hz": f, "int_bin": k, "peak_ratios": ratios,
                    "abs_shifts": shifts, "phases": phases,
                    "mean_peak_ratio": mean_ratio}

    shifts = best["abs_shifts"]
    bits_by_symbol = {}
    for n in range(1, n_symbols):
        d = (shifts[n] - shifts[n - 1]) % N_FFT
        bits_by_symbol[n] = relative_shift_to_bits(d, SIG_BITS_PER_SYMBOL)
    best["rel_shifts"] = {n: (shifts[n] - shifts[n - 1]) % N_FFT
                          for n in range(1, n_symbols)}
    best["bits"] = bits_by_symbol
    best["fields"] = unpack_signaling(bits_by_symbol) if n_symbols >= 4 else {}
    # A/321 5.3.3: symbol 0 always has absolute cyclic shift 0.  A timing error
    # of d samples moves EVERY absolute shift by the same d, so shifts[0] is a
    # direct readout of the residual timing error -- and the signaling bits,
    # which come from DIFFERENCES of absolute shifts, are immune to it.  That
    # is precisely why A/321 encodes them differentially.  So test shifts[0]
    # against a tolerance and report it, never against exact zero.
    resid = shifts[0] if shifts[0] < N_FFT // 2 else shifts[0] - N_FFT
    best["timing_residual"] = resid
    best["sym0_shift_ok"] = abs(resid) <= sym0_shift_tol
    # the final symbol is phase inverted; its correlation phase should be ~pi
    # away from the others once the reference inversion is included, which it
    # is, so all phases should agree.
    best["detected"] = best["mean_peak_ratio"] >= 8.0 and best["sym0_shift_ok"]
    return best
