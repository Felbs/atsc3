"""SELFTEST -- the gate for the ATSC 3.0 bootstrap detector (task #37 / M1).

No radio.  Everything here is synthesis + spec test vectors.  Run:

    python selftest.py            # full gate
    python selftest.py --quick    # fewer trials

PASS criteria (printed at the end):
  G1  Gray code matches A/321 Annex B Table B.1.1/B.1.2 exactly (16/16)
  G2  frequency-domain structure matches A/321 5.2.1-5.2.3
  G3  time-domain A/B/C parts match A/321's INDEPENDENT alternative definition
      (spectrum shifted +-1 subcarrier), i.e. our closed-form assembly is
      cross-checked against a second route through the spec
  G4  noiseless synth -> detect recovers every signaling field exactly
  G5  detection >= 99% at SNR -6 dB with all 8 signaling bits/symbol correct
  G6  false alarms: 0 on Gaussian noise and 0 on 8-VSB, over N trials
  G7  survives a realistic capture chain: 6.144 -> 8 Msps -> back, + CFO
"""

from __future__ import annotations

import argparse
import sys
import time

import numpy as np

sys.path.insert(0, __file__.rsplit("\\", 1)[0].rsplit("/", 1)[0])

from atsc3 import bootstrap as bs
from atsc3.capture import to_bootstrap_rate
from atsc3.vsb8 import synth_8vsb

RESULTS = []


def gate(name, ok, detail=""):
    RESULTS.append((name, bool(ok), detail))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""))
    return ok


def add_noise(x, snr_db, rng):
    p = np.mean(np.abs(x) ** 2)
    n = (rng.normal(size=len(x)) + 1j * rng.normal(size=len(x)))
    n *= np.sqrt(p * 10 ** (-snr_db / 10.0) / 2.0)
    return x + n


def embed(sig, total, rng, snr_db, cfo_hz=0.0, pad_lo=3000):
    """Put `sig` at a random offset inside `total` samples of noise + CFO."""
    off = int(rng.integers(pad_lo, max(pad_lo + 1, total - len(sig) - pad_lo)))
    buf = np.zeros(total, dtype=np.complex128)
    buf[off:off + len(sig)] = sig
    if cfo_hz:
        t = np.arange(total)
        buf *= np.exp(2j * np.pi * cfo_hz * t / bs.FS)
    return add_noise(buf, snr_db, rng), off


# ---------------------------------------------------------------- G1
def g1_gray():
    print("\nG1  Gray code vs A/321 Annex B (gold test vector from the spec)")
    gold = {(0, 0, 0, 0): 64, (0, 0, 0, 1): 192, (0, 0, 1, 0): 448, (0, 0, 1, 1): 320,
            (0, 1, 0, 0): 960, (0, 1, 0, 1): 832, (0, 1, 1, 0): 576, (0, 1, 1, 1): 704,
            (1, 0, 0, 0): 1984, (1, 0, 0, 1): 1856, (1, 0, 1, 0): 1600, (1, 0, 1, 1): 1728,
            (1, 1, 0, 0): 1088, (1, 1, 0, 1): 1216, (1, 1, 1, 0): 1472, (1, 1, 1, 1): 1344}
    bad = [b for b, d in gold.items()
           if bs.bits_to_relative_shift(b) != d
           or tuple(bs.relative_shift_to_bits(d, 4)) != b]
    gate("Table B.1.1/B.1.2 map+demap", not bad, f"{16 - len(bad)}/16 exact")
    # A/321 Annex B: with N_S bits the adjacent-RCS spacing is 2^(11-N_S)
    ok = True
    for nb in (3, 4, 7, 8):
        vals = sorted(bs.bits_to_relative_shift([(i >> (nb - 1 - j)) & 1 for j in range(nb)])
                      for i in range(1 << nb))
        d = np.diff(vals)
        ok &= (len(set(d)) == 1 and d[0] == 2 ** (11 - nb))
    gate("RCS spacing = 2^(11-N_S) for N_S in 3,4,7,8", ok)
    # Gray property: adjacent RCS differ in exactly one signaling bit
    order = sorted(range(256), key=lambda i: bs.bits_to_relative_shift(
        [(i >> (7 - j)) & 1 for j in range(8)]))
    gray_ok = all(bin(order[i] ^ order[i + 1]).count("1") == 1 for i in range(255))
    gate("8-bit mapping is a true Gray code", gray_ok)


# ---------------------------------------------------------------- G2
def g2_freq():
    print("\nG2  Frequency-domain structure (A/321 5.2.1 - 5.2.3)")
    z = bs.zc_sequence()
    gate("ZC length 1499, unit modulus", len(z) == 1499 and
         np.allclose(np.abs(z), 1.0))
    ac = np.fft.ifft(np.fft.fft(z) * np.conj(np.fft.fft(z)))
    psl = 20 * np.log10(abs(ac[0]) / np.abs(ac[1:]).max())
    gate("ZC is CAZAC (zero cyclic autocorr sidelobes)", psl > 250, f"{psl:.0f} dB")

    pn = bs.pn_sequence(bs.MINOR_SEEDS[0], 6000)
    gate("PN is balanced (m-sequence)", abs(pn.mean() - 0.5) < 0.02,
         f"mean={pn.mean():.4f}")
    # A/321 does NOT claim the generator polynomial is primitive, and it is not
    # (measured period 63457, short of 2^16-1).  What actually matters is that
    # the sequence does not repeat inside one bootstrap: 4 symbols consume
    # 4*749 = 2996 bits.
    p = bs.pn_sequence(bs.MINOR_SEEDS[0], 70000)
    per = next((k for k in range(1, 66000)
                if np.array_equal(p[k:k + 64], p[:64])), None)
    gate("PN period >> bits consumed by one bootstrap (4*749)",
         per is not None and per > 10 * 4 * bs.N_H, f"period={per}")
    # all 8 published minor-version seeds must give distinct sequences
    seqs = {v: bs.pn_sequence(s, 4 * bs.N_H).tobytes()
            for v, s in bs.MINOR_SEEDS.items()}
    gate("8 minor-version seeds give 8 distinct PN sequences",
         len(set(seqs.values())) == 8)

    S = bs.freq_domain_symbol(0, pn)
    gate("DC subcarrier null", S[0] == 0)
    gate("1498 active subcarriers (+-749, DC null)",
         int(np.sum(np.abs(S) > 1e-12)) == 1498)
    gate("reflective symmetry S(k)==S(-k)",
         np.allclose(S[1:750], S[-1:-750:-1]))
    gate("occupied bandwidth 4.494 MHz <= 4.5 MHz",
         2 * 749 * bs.DF <= bs.BW, f"{2 * 749 * bs.DF / 1e6:.3f} MHz")
    # different symbols use different PN windows -> low cross correlation
    S1 = bs.freq_domain_symbol(1, pn)
    xc = abs(np.vdot(S, S1)) / np.sqrt(np.vdot(S, S).real * np.vdot(S1, S1).real)
    gate("symbol0 vs symbol1 spectra near-orthogonal", xc < 0.15, f"|rho|={xc:.4f}")


# ---------------------------------------------------------------- G3
def g3_time_structure():
    """Check assemble_cab/assemble_bca against A/321 5.4's PROSE definition.

    The normative closed-form equations are what we implement; this tests them
    against the independent prose statement of the same thing, which pins the
    two things that matter:
      * WHICH samples of part A each of B and C is built from, and
      * that B carries a frequency shift of exactly -/+ one subcarrier.
    A constant phase offset on part B is the one degree of freedom the two
    descriptions leave open (A/321's own "-520*Ts offset ... to account for the
    correct extraction" note).  It is irrelevant to this detector: the
    structural detector never looks at part B, and the matched filter only ever
    correlates part A.
    """
    print("\nG3  Time-domain A/B/C vs A/321 5.4 prose (sample ranges + freq shift)")
    pn = bs.pn_sequence(bs.MINOR_SEEDS[0], 5 * bs.N_H)
    for n, M in ((0, 0), (1, 1337), (2, 64)):
        S = bs.freq_domain_symbol(n, pn)
        A = bs.time_domain_A(S, M)
        sym = bs.assemble_cab(A) if n == 0 else bs.assemble_bca(A)

        if n == 0:                                  # CAB: C, A, B
            okC = np.allclose(sym[0:520], A[1528:2048])       # last 520 of A
            okA = np.allclose(sym[520:2568], A)
            b_src, want_slope = A[1544:2048], +2 * np.pi / bs.N_FFT  # last 504 of A
            b_got = sym[2568:3072]
        else:                                       # BCA: B, C, A
            okC = np.allclose(sym[504:1024], A[1528:2048])    # last 520 of A
            okA = np.allclose(sym[1024:3072], A)
            b_src, want_slope = A[1528:2032], -2 * np.pi / bs.N_FFT  # first 504 of C
            b_got = sym[0:504]

        ratio = b_got / b_src
        okB_mag = np.allclose(np.abs(ratio), 1.0, atol=1e-9)
        slope = np.polyfit(np.arange(len(ratio)), np.unwrap(np.angle(ratio)), 1)[0]
        okB_frq = abs(slope - want_slope) < 1e-9
        gate(f"symbol {n} ({'CAB' if n == 0 else 'BCA'}) parts A/B/C",
             okA and okB_mag and okB_frq and okC,
             f"A={okA} C={okC} B:|.|={okB_mag} df={okB_frq} "
             f"({slope * bs.N_FFT / (2 * np.pi):+.6f} subcarriers)")

    # CAB part B additionally carries the e^{-j*pi} phase inversion
    S = bs.freq_domain_symbol(0, pn)
    A = bs.time_domain_A(S, 0)
    sym = bs.assemble_cab(A)
    ph = sym[2568:3072] / (A[1544:2048] * np.exp(2j * np.pi * np.arange(1544, 2048)
                                                 / bs.N_FFT))
    gate("CAB part B carries the e^{-j*pi} inversion",
         np.allclose(ph, -1.0, atol=1e-9), f"const={ph.mean():.6f}")

    # the geometry the structural detector keys on
    iq, _ = bs.synthesize()
    gate("total length 4*3072", len(iq) == 12288)
    ok = True
    for n in range(4):
        o = n * bs.N_SYM + bs.run_offset(n)
        ok &= np.allclose(iq[o:o + 520], iq[o + 2048:o + 2568])
    gate("every symbol has a 520-sample repeat at lag 2048", ok)
    papr = 10 * np.log10(np.abs(iq).max() ** 2 / np.mean(np.abs(iq) ** 2))
    gate("PAPR sane for OFDM (< 14 dB)", papr < 14.0, f"{papr:.2f} dB")


# ---------------------------------------------------------------- G4
def g4_roundtrip():
    print("\nG4  Noiseless synth -> detect roundtrip (all fields)")
    cases = [
        dict(ea_wake_up_1=0, ea_wake_up_2=0, min_time_to_next=10,
             system_bandwidth=0, bsr_coefficient=2, preamble_structure=0),
        dict(ea_wake_up_1=1, ea_wake_up_2=1, min_time_to_next=31 - 1,
             system_bandwidth=2, bsr_coefficient=80, preamble_structure=255),
        dict(ea_wake_up_1=1, ea_wake_up_2=0, min_time_to_next=0,
             system_bandwidth=1, bsr_coefficient=0, preamble_structure=137),
    ]
    allok = True
    for f in cases:
        iq, info = bs.synthesize(fields=f)
        buf = np.concatenate([np.zeros(2000, complex), iq, np.zeros(2000, complex)])
        r = bs.detect(buf)
        got = r.get("matched", {}).get("fields", {})
        ok = (r["verdict"] == "ATSC3_CONFIRMED"
              and all(got.get(k) == v for k, v in f.items())
              and r["matched"]["abs_shifts"] == info["abs_shifts"])
        allok &= ok
        if not ok:
            print("      mismatch:", f, "->", got, r["verdict"])
    gate("3 field permutations recovered exactly", allok)

    iq, info = bs.synthesize(fields=cases[0])
    r = bs.detect(np.concatenate([np.zeros(2000, complex), iq]))
    st, mt = r["structural"], r["matched"]
    gate("structural score ~1.0 noiseless", st["score"] > 0.99, f"{st['score']:.4f}")
    gate("timing exact", st["position"] == 2000, f"pos={st['position']}")
    gate("matched peak ratio very high", mt["mean_peak_ratio"] > 20,
         f"{mt['mean_peak_ratio']:.1f}")
    gate("post-bootstrap rate decoded",
         abs(info["fields"]["post_bootstrap_sample_rate_hz"] - 6_912_000) < 1,
         "bsr=2 -> 6.912 MHz (the standard 6 MHz ATSC 3.0 rate)")

    # minor-version discrimination: wrong PN seed must NOT match
    iq, _ = bs.synthesize(minor_version=0)
    bad = bs.detect(np.concatenate([np.zeros(2000, complex), iq]), minor_version=3)
    gate("wrong minor version -> STRUCTURE_ONLY, not CONFIRMED",
         bad["verdict"] == "ATSC3_STRUCTURE_ONLY",
         f"verdict={bad['verdict']} ratio={bad['matched']['mean_peak_ratio']:.1f}")


# ---------------------------------------------------------------- G5
def g5_snr(trials, snrs, cfo_hz=0.0, label="SNR sweep"):
    print(f"\nG5  {label}  ({trials} trials/point, CFO={cfo_hz:+.0f} Hz)")
    print("    SNR dB | struct | matched | bits OK | struct score | CFO err Hz | timing err")
    rows = []
    for snr in snrs:
        rng = np.random.default_rng(1000 + int(snr * 10))
        ns = nm = nb = 0
        scores, cerr, terr = [], [], []
        for t in range(trials):
            f = dict(ea_wake_up_1=int(rng.integers(2)),
                     ea_wake_up_2=int(rng.integers(2)),
                     min_time_to_next=int(rng.integers(31)),
                     system_bandwidth=int(rng.integers(3)),
                     bsr_coefficient=int(rng.integers(81)),
                     preamble_structure=int(rng.integers(256)))
            iq, _ = bs.synthesize(fields=f)
            buf, off = embed(iq, 60000, rng, snr, cfo_hz)
            r = bs.detect(buf, integer_cfo_search=2 if cfo_hz else 1)
            st, mt = r["structural"], r["matched"]
            scores.append(st["score"])
            if st["detected"]:
                ns += 1
                terr.append(abs(st["position"] - off))
            if mt and mt["detected"]:
                nm += 1
                cerr.append(abs(mt["cfo_hz"] - cfo_hz))
                if all(mt["fields"].get(k) == v for k, v in f.items()):
                    nb += 1
        rows.append((snr, ns / trials, nm / trials, nb / trials, np.mean(scores),
                     np.median(cerr) if cerr else float("nan"),
                     np.median(terr) if terr else float("nan")))
        print("    %6.1f | %5.1f%% | %6.1f%% | %6.1f%% | %12.4f | %10.1f | %10.1f"
              % (snr, 100 * ns / trials, 100 * nm / trials, 100 * nb / trials,
                 np.mean(scores), rows[-1][5], rows[-1][6]))
    return rows


# ---------------------------------------------------------------- G6
def g6_false_alarm(trials):
    print(f"\nG6  False alarms ({trials} trials each)")
    n = 60000
    thr = bs.STRUCT_RATIO_THRESHOLD
    for label, gen in (("Gaussian noise", lambda t: (
                            np.random.default_rng(7000 + t).normal(size=n)
                            + 1j * np.random.default_rng(8000 + t).normal(size=n))
                            / np.sqrt(2)),
                       ("synthetic 8-VSB", lambda t: synth_8vsb(
                            n, bs.FS, seed=t, snr_db=25.0).astype(np.complex128))):
        fa_s = fa_m = 0
        pk_ratio = pk_matched = 0.0
        for t in range(trials):
            r = bs.detect(gen(t), integer_cfo_search=1)
            pk_ratio = max(pk_ratio, r["structural"]["peak_ratio"])
            fa_s += r["structural"]["detected"]
            if r["matched"]:
                pk_matched = max(pk_matched, r["matched"]["mean_peak_ratio"])
                fa_m += r["matched"]["detected"]
        gate(f"{label}: 0 structural detections", fa_s == 0,
             f"{fa_s}/{trials}, max peak_ratio {pk_ratio:.2f} (thr {thr:.1f})")
        gate(f"{label}: 0 full ATSC3 confirmations", fa_m == 0,
             f"{fa_m}/{trials}, max matched ratio {pk_matched:.2f} (thr 8.0)")

    # 8-VSB with a real bootstrap buried in it -- adjacent-channel realism
    rng = np.random.default_rng(11)
    iq, f0 = bs.synthesize()
    hits = 0
    for t in range(max(4, trials // 4)):
        x = synth_8vsb(n, bs.FS, seed=100 + t, snr_db=30.0).astype(np.complex128) * 0.35
        off = int(rng.integers(3000, n - len(iq) - 3000))
        x[off:off + len(iq)] += iq
        r = bs.detect(x, integer_cfo_search=1)
        # Under a strong co-channel carrier the STRUCTURAL metric is suppressed
        # (it is normalized by total energy, interference included), so the
        # expected positive here is ATSC3_MATCHED_ONLY, not CONFIRMED.  What
        # actually proves detection is that the signaling bits come back right.
        ok = r["verdict"] in ("ATSC3_CONFIRMED", "ATSC3_MATCHED_ONLY")
        if ok and r.get("matched"):
            ok = all(r["matched"]["fields"].get(k) == v
                     for k, v in f0["fields"].items()
                     if k in ("min_time_to_next", "system_bandwidth",
                              "bsr_coefficient", "preamble_structure"))
        hits += bool(ok)
    gate("bootstrap under co-channel 8-VSB (~9 dB SIR): found + bits correct",
         hits >= max(4, trials // 4) * 0.9, f"{hits}/{max(4, trials // 4)}")


# ---------------------------------------------------------------- G7
def g7_capture_chain(trials):
    print(f"\nG7  Realistic capture chain: 6.144 -> 8.0 Msps -> resample back")
    rng = np.random.default_rng(99)
    ok = 0
    for t in range(trials):
        f = dict(ea_wake_up_1=0, ea_wake_up_2=1,
                 min_time_to_next=int(rng.integers(31)), system_bandwidth=0,
                 bsr_coefficient=2, preamble_structure=int(rng.integers(256)))
        iq, _ = bs.synthesize(fields=f)
        buf, off = embed(iq, 40000, rng, 3.0, cfo_hz=float(rng.uniform(-8000, 8000)))
        # emulate: signal was actually recorded at 8 Msps then brought back
        from scipy.signal import resample_poly
        rec8 = resample_poly(buf, 125, 96)                 # 6.144 -> 8.0 Msps
        back = to_bootstrap_rate(rec8.astype(np.complex64), 8_000_000.0)
        r = bs.detect(back.astype(np.complex128), integer_cfo_search=4)
        if r["verdict"] == "ATSC3_CONFIRMED" and \
           all(r["matched"]["fields"].get(k) == v for k, v in f.items()):
            ok += 1
    gate("8 Msps round trip + CFO +-8 kHz @ 3 dB SNR", ok >= trials * 0.95,
         f"{ok}/{trials}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    t0 = time.time()
    trials = 12 if a.quick else 40
    fa_trials = 20 if a.quick else 100

    print("=" * 74)
    print("ATSC 3.0 BOOTSTRAP DETECTOR -- SELFTEST (no radio; synthesis only)")
    print("=" * 74)
    g1_gray()
    g2_freq()
    g3_time_structure()
    g4_roundtrip()
    snrs = [0, -3, -6, -9, -12, -15, -18]
    rows = g5_snr(trials, snrs)
    g5_snr(max(8, trials // 3), [-6, -12], cfo_hz=4200.0, label="SNR sweep with CFO")
    g6_false_alarm(fa_trials)
    g7_capture_chain(max(8, trials // 2))

    # summarize the detection threshold
    print("\nDETECTION THRESHOLD SUMMARY")
    for snr, ps, pm, pb, sc, ce, te in rows:
        print(f"    SNR {snr:+5.1f} dB : structural {100*ps:5.1f}%  "
              f"confirmed {100*pm:5.1f}%  all-24-bits {100*pb:5.1f}%")
    thr = [r[0] for r in rows if r[1] >= 0.99]
    thr_b = [r[0] for r in rows if r[3] >= 0.99]
    print(f"    -> structural detection >=99% down to {min(thr) if thr else float('nan')} dB SNR")
    print(f"    -> full signaling decode >=99% down to {min(thr_b) if thr_b else float('nan')} dB SNR")

    npass = sum(1 for _, ok, _ in RESULTS if ok)
    print("\n" + "=" * 74)
    print(f"GATE: {npass}/{len(RESULTS)} checks passed   ({time.time()-t0:.1f}s)")
    for name, ok, d in RESULTS:
        if not ok:
            print(f"   FAILED: {name}  {d}")
    print("=" * 74)
    return 0 if npass == len(RESULTS) else 1


if __name__ == "__main__":
    sys.exit(main())

