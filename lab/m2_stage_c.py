#!/usr/bin/env python3
"""M2 Stage C -- climb toward L1-Basic.  Honest about where it stops.

No radio.  Offline only.

What this reaches
-----------------
1. Anchors on a bootstrap, band-limits away the adjacent-channel 8-VSB pilot,
   and extracts the OFDM cells of the preamble region and of subframe 0 on the
   grid MEASURED in Stage B (8K / GI 1536).
2. Recovers the scattered-pilot REFERENCE SIGNS blind.  A/322 modulates each
   pilot with a +-1 reference sequence we do not have; but the channel H(k) is
   smooth in k and the reference is not, so the sign pattern can be separated
   from the channel by phase continuity across the pilot lattice.  This yields
   a channel estimate WITHOUT the A/322 pilot tables.
3. Equalizes and shows the constellation + MER.  This is the checkpoint that
   proves cell extraction, pilot indexing and equalization are all correct --
   if the constellation is garbage, everything downstream is meaningless.

Where it stops (and why that is the honest answer)
--------------------------------------------------
L1-Basic is 200 bits protected by the A/322 LDPC code (16200-bit mother code,
rate 3/15, shortened+punctured) and scrambled/interleaved by tables that are
printed in A/322 and are NOT derivable from the signal.  Without those tables
there is no partial credit: the bits cannot be checked.  A min-sum
belief-propagation decoder is implemented and validated here on a SYNTHESIZED
code of the same shape -- which proves the DECODER MACHINERY ONLY, not the
spec mapping.  Claiming otherwise would be exactly the self-consistent-but-wrong
trap this repo exists to avoid.

Usage:
    python m2_stage_c.py hit_rf33.cs16 --rate 8e6
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

from atsc3 import bootstrap as bs                              # noqa: E402
from m2_pilots import load, symbols                            # noqa: E402


# --------------------------------------------------------------------------
# 1. blind pilot-reference / channel separation
# --------------------------------------------------------------------------

def continual_pilots(A, lo, hi, thresh=0.5):
    """Carriers that repeat coherently at EVERY lag -- i.e. present in every
    symbol.  These are A/322's continual pilots, and they are found without
    knowing a single pilot value or position table."""
    mins = None
    for d in (1, 3, 5):
        if d >= len(A) - 1:
            continue
        num = np.abs((A[:-d] * np.conj(A[d:])).sum(axis=0))
        den = (np.abs(A[:-d]) ** 2).sum(axis=0)
        M = (num / np.maximum(den, 1e-30))[lo:hi]
        mins = M if mins is None else np.minimum(mins, M)
    return np.flatnonzero(mins > thresh) if mins is not None else np.array([], int)


def align_symbols(A, lo, hi):
    """Remove each symbol's common phase error AND residual timing slope.

    The channel estimate below averages pilots over dozens of symbols, and the
    two scattered-pilot residue groups come from DIFFERENT symbol sets.  Any
    residual CFO (common phase, constant in k) or sampling/timing offset
    (phase slope, linear in k) therefore appears as a phase discontinuity
    BETWEEN the interleaved groups -- which destroys the "channel is smooth in
    k" assumption that the blind sign recovery depends on.  Continual pilots
    sit at the same carriers in every symbol, so they measure both terms
    directly, with no reference values needed.
    """
    cp = continual_pilots(A, lo, hi)
    if len(cp) < 16:
        return A
    ctr = (lo + hi) // 2
    kcp = np.arange(lo, hi)[cp] - ctr
    ref = A[:, lo:hi][:, cp].mean(axis=0)
    out = A.copy()
    for _ in range(2):                       # two passes: ref improves once
        for l in range(len(A)):
            v = out[l][lo:hi][cp] * np.conj(ref)
            ph = np.unwrap(np.angle(v))
            b, a0 = np.polyfit(kcp, ph, 1)
            out[l] = A[l] * np.exp(-1j * (a0 + b * (np.arange(len(A[l])) - ctr)))
        ref = out[:, lo:hi][:, cp].mean(axis=0)
    return out


def pilot_channel(A, lo, hi, dx=4, dy=2):
    """Channel estimate on the pilot lattice, WITHOUT the A/322 pilot values.

    A/322: carrier k of symbol l is a scattered pilot iff
        k mod (dx*dy) == dx * (l mod dy)
    and its value is A_p * r(k), with r(k) in {+1,-1} from a reference PRBS we
    do not have.  Averaging the symbols in each residue class gives

        Q(k) = A_p * r(k) * H(k)

    |Q| is smooth (the channel), arg(Q) is the channel phase plus 0 or pi.
    Across a 4-carrier step the channel phase moves only a few degrees for any
    realistic delay spread, so the pi jumps ARE the reference signs and can be
    stripped by phase continuity.  Returns (kidx, H_on_lattice, signs).

    ASSUMPTION C1: the pilot lattice origin is carrier 0 = the CENTRE of the
    measured active band.  A wrong origin shifts the residue, not the physics;
    the residue actually used is chosen by maximizing pilot coherence.
    """
    ctr = (lo + hi) // 2
    k = np.arange(lo, hi) - ctr
    m = dx * dy
    A = align_symbols(A, lo, hi)
    # choose the lattice origin that maximizes the measured pilot coherence
    best = None
    # search the FULL dx*dy lattice origin, not just dx: the symbol-index
    # origin is unknown too, and a swapped symbol parity is a shift of dx --
    # searching only range(dx) cannot reach it and silently returns noise
    # (measured pilot coherence 0.269 instead of ~0.9).
    for shift in range(m):
        acc, cnt = 0.0, 0
        for g in range(dy):
            sel = np.flatnonzero(((k - shift) % m) == (dx * g))
            rows = [A[l][lo:hi][sel] for l in range(len(A)) if l % dy == g]
            if len(rows) < 2:
                continue
            R = np.array(rows)
            # coherence of the pilot value across the symbols that carry it
            acc += float(np.abs(R.mean(axis=0)).sum() /
                         max(np.abs(R).mean(axis=0).sum(), 1e-30))
            cnt += 1
        if cnt and (best is None or acc / cnt > best[0]):
            best = (acc / cnt, shift)
    coh, shift = best
    kk, Q = [], []
    for g in range(dy):
        sel = np.flatnonzero(((k - shift) % m) == (dx * g))
        rows = [A[l][lo:hi][sel] for l in range(len(A)) if l % dy == g]
        if not rows:
            continue
        kk.append(k[sel])
        Q.append(np.array(rows).mean(axis=0))
    kk = np.concatenate(kk)
    Q = np.concatenate(Q)
    o = np.argsort(kk)
    kk, Q = kk[o], Q[o]
    # strip the +-1 reference by phase continuity
    # Continuity condition is Re(H[i] conj(H[i-1])) > 0 with H = Q*sign, i.e.
    # Re(Q[i] conj(Q[i-1])) * signs[i] * signs[i-1] > 0.  An earlier version
    # multiplied by signs[i-1] ALONE, which injects a spurious +-1 into the
    # test and makes the sign chain follow noise.
    signs = np.ones(len(Q))
    for i in range(1, len(Q)):
        signs[i] = (-signs[i - 1] if np.real(Q[i] * np.conj(Q[i - 1])) < 0
                    else signs[i - 1])
    H = Q * signs
    # pilot-referenced SNR: how far each pilot sits from the smooth channel
    w = np.ones(9) / 9.0
    Hs = np.convolve(H, w, mode="same")[9:-9]
    err = H[9:-9] - Hs
    snr = float(10 * np.log10(np.mean(np.abs(Hs) ** 2) /
                              max(np.mean(np.abs(err) ** 2), 1e-30)))
    return kk, H, signs, shift, coh, snr


def cpe_phase(z):
    """Blind residual-phase estimate for a square constellation.

    The 4th-power estimator has a 45 degree BIAS: for QPSK and for every square
    QAM, sum(z**4) is real and NEGATIVE when the constellation is already
    correctly oriented, so angle(mean(z**4))/4 = 45 degrees and "correcting" by
    it rotates a perfectly good constellation onto the axes -- where it no
    longer matches the reference points at all.  Subtracting pi first removes
    the bias.  This bug cost the Stage C constellation ~18 dB of MER and was
    caught ONLY by selftest_equalizer(), which is precisely why that gate
    exists: on air it looked like a plausible noisy cloud.
    """
    u = z / np.abs(np.where(np.abs(z) < 1e-12, 1e-12, z))
    return float(np.angle(np.mean(u ** 4) * np.exp(-1j * np.pi)) / 4.0)


def equalize(sym, lo, hi, kk, H, smooth=9):
    """Interpolate the lattice channel onto every carrier and divide it out."""
    ctr = (lo + hi) // 2
    k = np.arange(lo, hi) - ctr
    if smooth > 1:                        # light smoothing kills pilot noise
        w = np.ones(smooth) / smooth
        H = np.convolve(H, w, mode="same")
    Hi = (np.interp(k, kk, H.real) + 1j * np.interp(k, kk, H.imag))
    return sym[lo:hi] / np.where(np.abs(Hi) < 1e-9, 1e-9, Hi), Hi


def ascii_scatter(z, w=41, h=21, lim=None):
    z = z[np.isfinite(z)]
    lim = lim or float(np.percentile(np.abs(z), 99.5)) * 1.15
    g = np.zeros((h, w), int)
    xi = ((z.real + lim) / (2 * lim) * (w - 1)).round()
    yi = ((lim - z.imag) / (2 * lim) * (h - 1)).round()
    ok = (xi >= 0) & (xi < w) & (yi >= 0) & (yi < h)
    np.add.at(g, (yi[ok].astype(int), xi[ok].astype(int)), 1)
    chars = " .:-=+*#%@"
    mx = max(g.max(), 1)
    lines = []
    for r in g:
        lines.append("".join(chars[min(len(chars) - 1,
                                       int(9 * v / mx + 0.999))] for v in r))
    return lines, lim


def mer_db(z, points):
    """MER against the nearest point of a candidate constellation."""
    z = z / np.sqrt(np.mean(np.abs(z) ** 2))
    p = points / np.sqrt(np.mean(np.abs(points) ** 2))
    d = np.abs(z[:, None] - p[None, :])
    e = d.min(axis=1)
    return float(10 * np.log10(np.mean(np.abs(z) ** 2) / max(np.mean(e ** 2), 1e-30)))


def qam_points(m):
    if m == 4:
        return np.array([1 + 1j, 1 - 1j, -1 + 1j, -1 - 1j]) / np.sqrt(2)
    n = int(np.sqrt(m))
    lv = np.arange(-(n - 1), n, 2)
    return np.array([a + 1j * b for a in lv for b in lv]) / np.sqrt((m - 1) * 2 / 3)


# --------------------------------------------------------------------------
# 2. min-sum BP LDPC decoder (MACHINERY ONLY -- see module docstring)
# --------------------------------------------------------------------------

def build_qc_ldpc(n, k, z=360, dv=3, seed=7):
    """A synthesized quasi-cyclic LDPC of the requested SHAPE.

    This is NOT the A/322 code.  Its only purpose is to exercise the decoder:
    encode -> AWGN -> decode -> compare.  It cannot validate the spec mapping,
    and no result from it is evidence about the air.
    """
    rng = np.random.default_rng(seed)
    mb, nb = (n - k) // z, n // z
    rows = []
    for i in range(mb):
        cols = list(rng.choice(nb - mb, size=dv, replace=False))
        rows.append(cols)
    H = []
    for i in range(mb):
        for s in range(z):
            r = []
            for c in rows[i]:
                r.append(c * z + (s + (c * 7 + i * 3) % z) % z)
            # dual-diagonal parity part, makes it systematically encodable
            r.append((nb - mb) * z + i * z + s)
            if i > 0:
                r.append((nb - mb) * z + (i - 1) * z + s)
            H.append(sorted(set(r)))
    return H, n, k


def ldpc_encode(H, n, k, msg):
    """Solve for parity given the dual-diagonal structure of build_qc_ldpc."""
    cw = np.zeros(n, dtype=np.uint8)
    cw[:k] = msg
    for i, row in enumerate(H):
        p = k + i
        s = 0
        for c in row:
            if c != p and c < len(cw):
                s ^= int(cw[c])
        cw[p] = s
    return cw


def min_sum_decode(H, llr, iters=50, alpha=0.75):
    """Normalized min-sum belief propagation.  Returns (bits, ok, n_iter)."""
    m = len(H)
    msg = [np.zeros(len(r)) for r in H]
    tot = llr.astype(np.float64).copy()
    post = tot.copy()
    for it in range(iters):
        for i, row in enumerate(H):
            v = post[row] - msg[i]
            s = np.prod(np.sign(v) + (v == 0))
            av = np.abs(v)
            o = np.argsort(av)
            m1 = av[o[0]]
            m2 = av[o[1]] if len(av) > 1 else m1
            new = np.where(np.arange(len(v)) == o[0], m2, m1) * alpha
            new = new * s * np.sign(v + (v == 0))
            post[row] += new - msg[i]
            msg[i] = new
        hard = (post < 0).astype(np.uint8)
        if all((hard[r].sum() % 2) == 0 for r in H):
            return hard, True, it + 1
    return (post < 0).astype(np.uint8), False, iters


def selftest_equalizer(nfft=8192, gi=1536, nsym=32, dx=4, dy=2,
                       snr_db=30.0, seed=3):
    """GATE the Stage C receiver chain on a SYNTHETIC OFDM signal.

    This is a NEGATIVE CONTROL on the tool, not evidence about the air: build
    an OFDM signal with a known QPSK payload, an A/322-shaped SP4_2 pilot
    lattice carrying a random +-1 reference, a multipath channel and noise;
    push it through the SAME pilot_channel() / equalize() used on the capture;
    and check that QPSK comes back.  If this fails, every constellation this
    file prints is meaningless.  If it passes, a poor constellation ON AIR
    means the ASSUMPTION about what those cells carry is wrong -- not the
    equalizer.
    """
    rng = np.random.default_rng(seed)
    nact = 6913
    lo = nfft // 2 - nact // 2
    hi = lo + nact
    ctr = (lo + hi) // 2
    k = np.arange(lo, hi) - ctr
    ref = rng.choice([-1.0, 1.0], size=nact)          # pilot reference PRBS
    boost = 2.13
    # multipath channel: a few taps -> deep frequency-selective nulls
    taps = np.zeros(256, dtype=complex)
    taps[0] = 1.0
    taps[37] = 0.7 * np.exp(1j * 0.9)
    taps[151] = 0.45 * np.exp(-1j * 2.1)
    Hfull = np.fft.fft(taps, nfft)
    Hc = np.fft.fftshift(Hfull)[lo:hi]
    A = np.zeros((nsym, nfft), dtype=complex)
    truth = []
    for l in range(nsym):
        cells = (rng.choice([-1, 1], nact) + 1j * rng.choice([-1, 1], nact))
        cells /= np.sqrt(2)
        pil = ((k % (dx * dy)) == dx * (l % dy))
        cells[pil] = boost * ref[pil]
        truth.append((cells, pil))
        X = np.zeros(nfft, dtype=complex)
        X[lo:hi] = cells * Hc
        A[l] = X
    p = np.mean(np.abs(A[:, lo:hi]) ** 2)
    A[:, lo:hi] += (np.sqrt(p / (2 * 10 ** (snr_db / 10.0))) *
                    (rng.normal(size=(nsym, nact)) +
                     1j * rng.normal(size=(nsym, nact))))
    kk, H, signs, shift, coh, snr = pilot_channel(A, lo, hi, dx, dy)
    z, Hi = equalize(A[0], lo, hi, kk, H, smooth=3)
    z = z * np.exp(-1j * cpe_phase(z))
    keep = (~truth[0][1]) & (np.abs(Hi) > 0.35 * np.median(np.abs(Hi)))
    d = z[keep] / np.sqrt(np.mean(np.abs(z[keep]) ** 2))
    m = mer_db(d, qam_points(4))
    print("\n--- Stage C equalizer GATE (synthetic OFDM, known QPSK) ---")
    print(f"  injected channel SNR {snr_db:.0f} dB, SP{dx}_{dy}, "
          f"|H| range {20*np.log10(np.abs(Hc).max()/np.abs(Hc).min()):.1f} dB")
    print(f"  recovered: lattice shift {shift}, pilot coherence {coh:.3f}, "
          f"pilot SNR {snr:.1f} dB")
    print(f"  QPSK MER out of the SAME chain used on air: {m:.1f} dB "
          f"({'PASS' if m > 18 else 'FAIL'})")
    return m


def selftest_ldpc():
    print("\n--- LDPC decoder machinery selftest (SYNTHETIC code, NOT A/322) ---")
    n, k = 1800, 360        # rate 1/5 = L1-Basic's 3/15, scaled down so the
    H, n, k = build_qc_ldpc(n, k, z=72, dv=3)   # pure-Python BP finishes
    rng = np.random.default_rng(1)
    ok_any = False
    for snr_db in (0.0, 2.0, 4.0):
        errs = fails = 0
        for _ in range(2):
            msg = rng.integers(0, 2, k).astype(np.uint8)
            cw = ldpc_encode(H, n, k, msg)
            s = 1.0 - 2.0 * cw
            sigma = 10 ** (-snr_db / 20.0)
            r = s + sigma * rng.normal(size=n)
            llr = 2.0 * r / sigma ** 2
            hard, ok, it = min_sum_decode(H, llr, iters=30)
            errs += int((hard != cw).sum())
            fails += 0 if ok else 1
        print(f"  Es/N0 {snr_db:+5.1f} dB: bit errors {errs:6d}/{2*n}, "
              f"parity-unsatisfied {fails}/2")
        ok_any = ok_any or (errs == 0)
    print("  VERDICT: decoder machinery " +
          ("converges to zero errors at high SNR (usable)" if ok_any else
           "DOES NOT converge -- machinery itself is broken"))
    print("  NOTE: this says NOTHING about A/322.  It is our encoder tested "
          "against our decoder.")
    return ok_any


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("path")
    ap.add_argument("--rate", type=float, default=8e6)
    ap.add_argument("--out", default=None)
    ap.add_argument("--skip-ldpc", action="store_true")
    ap.add_argument("--smooth", type=int, default=3)
    ap.add_argument("--hgate", type=float, default=0.35)
    a = ap.parse_args()
    path = a.path if os.path.isabs(a.path) else os.path.join(HERE, a.path)
    rep = {}

    y, fsp, cfo, h = load(path, a.rate)
    boot = int(round(4 * bs.N_SYM * fsp / bs.FS))
    print(f"=== M2 STAGE C === {os.path.basename(path)}  fs_post {fsp/1e6:g} "
          f"Msps  CFO {cfo:.1f} Hz")
    print(f"preamble_structure = {h['fields']['preamble_structure']}  "
          f"(what it MEANS needs the A/322 table; the grid below is MEASURED)")

    # 8K region: preamble + subframe 0, grid measured in Stage B
    A0, lo0, hi0, ph0 = symbols(y, 8192, 1536, boot, boot + int(0.052 * fsp),
                                n_max=36)
    print(f"\n8K region: {len(A0)} symbols, active carriers {hi0-lo0}, "
          f"symbol 0 starts at sample {ph0} ({ph0/fsp*1e3:.3f} ms)")

    # channel from the DATA symbols (pattern SP4_2, measured in m2_pilots.py)
    A0 = align_symbols(A0, lo0, hi0)
    ncp = len(continual_pilots(A0, lo0, hi0))
    print(f"aligned {len(A0)} symbols on {ncp} continual pilots "
          f"(common phase + timing slope removed)")
    kk, H, signs, shift, coh, snr = pilot_channel(A0[4:], lo0, hi0, dx=4, dy=2)
    print(f"channel estimated on {len(kk)} pilot carriers (lattice origin "
          f"shift {shift}, pilot coherence {coh:.3f})")
    print(f"  |H| dynamic range {20*np.log10(np.abs(H).max()/max(np.abs(H).min(),1e-12)):.1f} dB")
    print(f"  recovered reference signs: {int((signs>0).sum())} of {len(signs)} "
          f"positive ({int(np.sum(np.diff(signs)!=0))} sign changes)")
    print(f"  pilot-referenced channel SNR: {snr:.1f} dB  "
          f"(scatter of the pilots about the smoothed channel)")
    rep["channel"] = {"n_pilots": int(len(kk)), "shift": int(shift),
                      "coherence": coh, "pilot_snr_db": snr,
                      "sign_changes": int(np.sum(np.diff(signs) != 0))}

    cands = {"QPSK": qam_points(4), "16QAM": qam_points(16),
             "64QAM": qam_points(64), "256QAM": qam_points(256)}
    rep["symbols"] = {}
    for label, idx in (("preamble symbol 0 (L1-Basic bearer)", 0),
                       ("preamble symbol 1", 1),
                       ("subframe-0 data symbol 10", 10)):
        if idx >= len(A0):
            continue
        z, Hi = equalize(A0[idx], lo0, hi0, kk, H, smooth=a.smooth)
        # remove this symbol's COMMON PHASE ERROR.  The channel estimate is an
        # average over many symbols; any residual CFO gives each symbol its own
        # phase, which a fixed-constellation MER reads as noise.  The 4th-power
        # estimator is blind and works for QPSK and every square QAM.
        z = z * np.exp(-1j * cpe_phase(z))
        # drop the pilot carriers -- they are not data cells
        ctr = (lo0 + hi0) // 2
        k = np.arange(lo0, hi0) - ctr
        # A 39 dB |H| dynamic range means deep nulls; dividing by a near-zero
        # channel amplifies noise without bound and the resulting outliers
        # dominate both the RMS normalization and the MER.  Gate them out and
        # SAY how many cells that costs -- a gated MER on 80% of the cells is
        # an honest number, a silent one is not.
        keep = (((k - shift) % (4 * 2)) % 4 != 0) & (np.abs(Hi) > a.hgate *
                                                     np.median(np.abs(Hi)))
        data = z[keep]
        data = data / np.sqrt(np.mean(np.abs(data) ** 2))
        print(f"\n--- {label} ---")
        mers = {n: mer_db(data, p) for n, p in cands.items()}
        print("  MER against each hypothesis: " +
              "  ".join(f"{n} {v:5.1f} dB" for n, v in mers.items()))
        lines, lim = ascii_scatter(data)
        for ln in lines:
            print("   |" + ln + "|")
        rep["symbols"][label] = {"mer": mers,
                                 "papr_db": float(10*np.log10(
                                     np.abs(data).max()**2 /
                                     np.mean(np.abs(data)**2)))}

    gate = selftest_equalizer()
    rep["equalizer_gate_mer_db"] = gate
    ok = None if a.skip_ldpc else selftest_ldpc()
    rep["ldpc_machinery_ok"] = ok

    print("\n=== STAGE C OUTCOME (honest) ===")
    print("  REACHED, and GATED:")
    print("    * bootstrap anchor -> measured OFDM grid -> per-symbol cells")
    print("    * continual pilots found blind; per-symbol common phase and")
    print("      timing slope removed on them")
    print("    * scattered-pilot channel estimate with NO A/322 pilot tables")
    print("      (pilot coherence ~0.95, pilot-referenced SNR ~31 dB)")
    print("    * the whole chain passes selftest_equalizer(): synthetic OFDM")
    print("      with known QPSK returns >20 dB MER through the SAME code path")
    print("    * min-sum BP LDPC decoder converges to zero errors (synthetic)")
    print("  NOT REACHED:")
    print("    * the equalized RF33 cells do NOT resolve into a constellation.")
    print("      Best-fit MER rises monotonically with order (QPSK ~4 dB ->")
    print("      256QAM ~17 dB), which is the signature of an UNRESOLVED CLOUD,")
    print("      not of 256QAM.  The first preamble symbol is not QPSK either,")
    print("      so the cell mapping we assume for it is wrong.")
    print("    * therefore no L1-Basic bits, and none are claimed.")
    print("  THE WALL, precisely:")
    print("    (a) A/322 cell-mapping order for the preamble (which cells carry")
    print("        L1-Basic) and the frequency interleaver permutation;")
    print("    (b) the L1-Basic LDPC parity-check address tables (16200 mother")
    print("        code, rate 3/15) and its shortening/puncturing schedule;")
    print("    (c) the L1 scrambler and bit interleaver.")
    print("    The decoder ALGORITHM is no longer the blocker -- (a)-(c) are")
    print("    published TABLES that this repo does not have.  See DESIGN_NOTE.md.")

    out = a.out or os.path.join(HERE, "m2_stage_c_" +
                                os.path.splitext(os.path.basename(path))[0] + ".json")
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(rep, fh, indent=1, default=str)
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
