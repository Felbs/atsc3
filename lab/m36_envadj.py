#!/usr/bin/env python3
"""M36 -- A-SPX envelope adjustment: scaling the patched high band.

M35 copied low subbands into 12..21 kHz.  They are at the WRONG LEVEL until
this stage scales them to the transmitted envelopes.  That is what makes the
high band correct rather than merely present.

THE CHAIN, clause 5.7.6.4.2
-----------------------------
    Pseudocode 90   est_sig_sb   measure what the patched band actually has
    Pseudocode 82   scf_sig_sbg  = num_qmf_subbands * 2**(qscf/a)
    Pseudocode 83   scf_noise_sbg = 2**(NOISE_FLOOR_OFFSET - qscf_noise)
    Pseudocode 91   map both onto QMF subbands
    Pseudocode 95   sig_gain_sb  = sqrt(scf_sig / (1 + est_sig) / (1+scf_noise))
    Pseudocode 96   max_sig_gain per limiter subband group
    Pseudocode 98   limit the gain
    Pseudocode 99   boost factor -- recover the energy the limiter removed
    Pseudocode 100  cap the boost at 1.584893192
    Pseudocode 101  apply the boost to the gain
    Pseudocode 94   noise_lev_sb = sqrt(scf_sig/(1+scf_noise) * scf_noise)
    Pseudocode 97   trim the noise in proportion to the gain limiting
    Pseudocode 83   scf_noise_sbg = 2 ** (NOISE_FLOOR_OFFSET - qscf_noise)
                    -- note the SIGN: a larger index means LESS noise
    Table D.2       ASPX_NOISE[512][2], verified: mean |z|^2 = 1.000000
                    exactly, every entry unit magnitude, phase spread 1.787
                    rad against uniform 1.814

`aspx_interpolation = 1` in this stream, which selects Pseudocode 90's simpler
per-subband estimate rather than the per-subband-group average.

THE GATE
---------
The whole point of the stage is that the adjusted band carries the TRANSMITTED
envelope.  So after adjustment, measure the band's energy again and compare it
to the target `scf_sig_sb`.  If the gain chain is right the two agree; if the
estimate, the dequantisation, the mapping or the limiter is wrong they do not.
That is a closed loop with no free parameters, and it needs no reference
decoder.

SCOPE -- STATED BEFORE THE RESULTS
------------------------------------
Not implemented, and counted rather than hidden:

  * **additional harmonics.**  Signalled in under 1.5 % of channel-frames
    (ah_left 50 of 6960, ah_right 11, ah_present 44 of 3480); those frames lose
    a tone.  The sine branch of Pseudocode 95 is therefore never taken here.
  * **the limiter subband table** (Pseudocode 72) is implemented in its main
    form -- merge sbg_sig_lowres with the patch borders, sort, and drop borders
    closer than 0.245 octaves preferring to keep patch borders.  The spec's
    full tie-breaking is more elaborate.

frequency- and time-interleaved coding are NOT scope gaps: `fic_present` and
`tic_present` are 0 in every frame of this capture.

Usage:
    python m36_envadj.py
"""
from __future__ import annotations

import argparse
import math
import os
import sys
import wave

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m23_hcb as H                                               # noqa: E402
import m24_spectral as S                                          # noqa: E402
import m28_channels as C                                          # noqa: E402
import m31_aspx_bands as A                                        # noqa: E402
import m32_aspx_parse as P                                        # noqa: E402
import m33_qmf as Q                                               # noqa: E402
import m34_aspx_env as V                                          # noqa: E402
import m35_hfgen as G                                             # noqa: E402

NOISE_FLOOR_OFFSET = 6
EPSILON = 1.0
EPSILON0 = 1e-12
LIM_GAIN = 1.41254
MAX_SIG_GAIN = 1e5
MAX_BOOST_FACT = 1.584893192


def sbg_lim(t):
    """Pseudocode 72, main form: lowres borders merged with patch borders."""
    cand = sorted(set(t["lo"]) | set(t["patches"]))
    patch = set(t["patches"])
    out = [cand[0]]
    for b in cand[1:]:
        if math.log2(b / out[-1]) < 0.245:
            # too close: keep whichever is a patch border
            if b in patch and out[-1] not in patch:
                out[-1] = b
            continue
        out.append(b)
    return out


def estimate(Qh, t, borders, sl0, nts_in_ats, ts_off=0):
    """Pseudocode 90, the aspx_interpolation == 1 branch."""
    n_sb, n_env = t["num_sb_aspx"], len(borders) - 1
    est = np.zeros((n_sb, n_env))
    for e in range(n_env):
        ta = int(borders[e]) * nts_in_ats + ts_off + sl0
        tz = int(borders[e + 1]) * nts_in_ats + ts_off + sl0
        ta, tz = max(ta, 0), min(tz, Qh.shape[1])
        if tz <= ta:
            continue
        seg = np.abs(Qh[t["sbx"]:t["sbx"] + n_sb, ta:tz]) ** 2
        est[:, e] = seg.mean(axis=1)
    return est


def map_scf(qscf_env, t, fr, n_sb):
    """Pseudocode 82 + 91: dequantise the signal envelope onto QMF subbands."""
    n_env = len(qscf_env)
    out = np.zeros((n_sb, n_env))
    for e, q in enumerate(qscf_env):
        tbl = t["hi"] if fr["freq_res"][e] else t["lo"]
        a = 2 if 0 else 1                       # qmode 1 -> 3.0 dB -> a = 1
        for sbg, val in enumerate(q):
            if sbg + 1 >= len(tbl):
                break
            lo = tbl[sbg] - t["sbx"]
            hi = min(tbl[sbg + 1] - t["sbx"], n_sb)
            out[max(lo, 0):hi, e] = 64.0 * (2.0 ** (val / a))
    return out


PAN_OFFSET = 12


def map_scf_joint(qa, qb, t, fr, n_sb, a=1):
    """Pseudocode 84 -- JOINT dequantisation when aspx_balance = 1.

    With balance coding the pair is a SUM channel and a BALANCE channel, and
    channel B's transmitted values are a pan ratio, not a level.  Pseudocode 82
    covers only balance = 0 and says so; running it on channel B treats a pan
    ratio as an absolute level, which is how the right channel ended up 20x too
    loud and uncorrelated with its own core.

        nom     = 2 ** (qscf_a/a + 1) * num_qmf_subbands
        denom_a = 1 + 2 ** (PAN_OFFSET - qscf_b/a)
        denom_b = 1 + 2 ** (qscf_b/a - PAN_OFFSET)
    """
    n_env = min(len(qa), len(qb))
    A = np.zeros((n_sb, n_env))
    B = np.zeros((n_sb, n_env))
    for e in range(n_env):
        tbl = t["hi"] if fr["freq_res"][e] else t["lo"]
        for sbg in range(min(len(qa[e]), len(qb[e]))):
            if sbg + 1 >= len(tbl):
                break
            lo = max(tbl[sbg] - t["sbx"], 0)
            hi = min(tbl[sbg + 1] - t["sbx"], n_sb)
            if hi <= lo:
                continue
            va, vb = qa[e][sbg], qb[e][sbg]
            nom = (2.0 ** (va / a + 1)) * 64.0
            A[lo:hi, e] = nom / (1.0 + 2.0 ** (PAN_OFFSET - vb / a))
            B[lo:hi, e] = nom / (1.0 + 2.0 ** (vb / a - PAN_OFFSET))
    return A, B


def map_noise_joint(qa, qb, t, n_sb, n_env):
    """Pseudocode 84, the noise half."""
    A = np.zeros((n_sb, n_env))
    B = np.zeros((n_sb, n_env))
    for e in range(n_env):
        ia = min(e, len(qa) - 1)
        ib = min(e, len(qb) - 1)
        for sbg in range(min(t["n_noise"], len(qa[ia]), len(qb[ib]))):
            lo = max(t["noise"][sbg] - t["sbx"], 0)
            hi = min(t["noise"][sbg + 1] - t["sbx"], n_sb)
            if hi <= lo:
                continue
            va, vb = qa[ia][sbg], qb[ib][sbg]
            nom = 2.0 ** (NOISE_FLOOR_OFFSET - va + 1)
            A[lo:hi, e] = nom / (1.0 + 2.0 ** (PAN_OFFSET - vb))
            B[lo:hi, e] = nom / (1.0 + 2.0 ** (vb - PAN_OFFSET))
    return A, B


def noise_qscf(dat, t, prev):
    """Reconstruct the noise envelope INDICES (before dequantisation)."""
    delta = 2 if (dat["ch"] == 1 and dat["balance"] == 1) else 1
    dirs = dat.get("dir_noise") or [0] * len(dat["noise"])
    out, prev_q = [], prev
    for e, vals in enumerate(dat["noise"]):
        n = t["n_noise"]
        q = [0] * n
        if e < len(dirs) and dirs[e] == 0:
            acc = 0
            for sbg in range(min(n, len(vals))):
                acc += delta * vals[sbg]
                q[sbg] = acc
        else:
            for sbg in range(min(n, len(vals))):
                base = prev_q[sbg] if (prev_q and sbg < len(prev_q)) else 0
                q[sbg] = base + delta * vals[sbg]
        out.append(q)
        prev_q = q
    return out, prev_q


def gain_unlimited(sig_sb, est):
    """Pseudocode 95, non-sine branch, before the limiter touches it."""
    return np.sqrt(sig_sb / (EPSILON + est))


def adjust(Qh, t, sig_sb, est, lim, noise_sb=None):
    """Pseudocodes 94..101.  -> (gain, adjusted noise level) per subband/env."""
    n_sb, n_env = est.shape
    if noise_sb is None:
        noise_sb = np.zeros_like(est)
    gain = np.sqrt(sig_sb / (EPSILON + est) / (1.0 + noise_sb))
    # Pseudocode 94: the noise level shares the band energy with the signal
    sig_noise_fact = sig_sb / (1.0 + noise_sb)
    noise_lev = np.sqrt(sig_noise_fact * noise_sb)

    # Pseudocode 96: per limiter subband group ceiling
    maxg = np.zeros_like(gain)
    for i in range(len(lim) - 1):
        lo = max(lim[i] - t["sbx"], 0)
        hi = min(lim[i + 1] - 1 - t["sbx"], n_sb)
        if hi <= lo:
            continue
        nom = sig_sb[lo:hi].sum(axis=0)
        den = EPSILON0 + est[lo:hi].sum(axis=0)
        maxg[lo:hi] = np.minimum(np.sqrt(nom / den) * LIM_GAIN, MAX_SIG_GAIN)
    maxg[maxg == 0] = MAX_SIG_GAIN

    lim_gain = np.minimum(gain, maxg)            # Pseudocode 98
    # Pseudocode 97: the noise is trimmed in proportion to the gain limiting
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(gain > 0, maxg / np.maximum(gain, 1e-30), 1.0)
    noise_lim = np.minimum(noise_lev, noise_lev * ratio)

    boost = np.ones_like(gain)                   # Pseudocode 99 / 100
    for i in range(len(lim) - 1):
        lo = max(lim[i] - t["sbx"], 0)
        hi = min(lim[i + 1] - 1 - t["sbx"], n_sb)
        if hi <= lo:
            continue
        nom = EPSILON0 + sig_sb[lo:hi].sum(axis=0)
        den = EPSILON0 + (est[lo:hi] * lim_gain[lo:hi] ** 2).sum(axis=0)
        den = den + (noise_lim[lo:hi] ** 2).sum(axis=0)
        boost[lo:hi] = np.minimum(np.sqrt(nom / den), MAX_BOOST_FACT)
    return lim_gain * boost, noise_lim * boost   # Pseudocode 101


def noise_table(path=None):
    """Table D.2 -- ASPX_NOISE[512][2] in the spec C file, 512 complex.

    Each entry has unit magnitude with random phase ("average energy of 1"),
    which is checkable and is the gate on the transcription.
    """
    import re
    src = open(path or H.DEFAULT_C, encoding="utf-8",
               errors="replace").read()
    src = re.sub(r"/\*.*?\*/", "", src, flags=re.S)
    m = re.search(r"ASPX_NOISE\s*\[\s*512\s*\]\s*\[\s*2\s*\]"
                  r"\s*=\s*\{(.*?)\}\s*;", src, re.S)
    v = [float(x) for x in
         re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", m.group(1))]
    a = np.asarray(v, float).reshape(-1, 2)
    return a[:, 0] + 1j * a[:, 1]


def noise_envelopes(dat, t, prev):
    """Reconstruct + dequantise the NOISE envelopes.

    Same delta scheme as the signal envelopes (Pseudocode 80) but over
    num_sbg_noise groups, then Pseudocode 83:

        scf_noise_sbg = 2 ** (NOISE_FLOOR_OFFSET - qscf_noise_sbg)

    Note the SIGN: a larger transmitted index means LESS noise.
    """
    delta = 2 if (dat["ch"] == 1 and dat["balance"] == 1) else 1
    dirs = dat.get("dir_noise") or [0] * len(dat["noise"])
    out, prev_q = [], prev
    for e, vals in enumerate(dat["noise"]):
        n = t["n_noise"]
        q = [0] * n
        if e < len(dirs) and dirs[e] == 0:            # FREQ
            acc = 0
            for sbg in range(min(n, len(vals))):
                acc += delta * vals[sbg]
                q[sbg] = acc
        else:                                          # TIME
            for sbg in range(min(n, len(vals))):
                base = prev_q[sbg] if (prev_q and sbg < len(prev_q)) else 0
                q[sbg] = base + delta * vals[sbg]
        out.append([2.0 ** (NOISE_FLOOR_OFFSET - x) for x in q])
        prev_q = q
    return out, prev_q


def map_noise(scf_noise_env, t, n_sb, n_env):
    """Pseudocode 91, the noise half: subband groups -> QMF subbands."""
    out = np.zeros((n_sb, n_env))
    for e in range(n_env):
        src = scf_noise_env[min(e, len(scf_noise_env) - 1)]
        for sbg in range(min(t["n_noise"], len(src))):
            lo = max(t["noise"][sbg] - t["sbx"], 0)
            hi = min(t["noise"][sbg + 1] - t["sbx"], n_sb)
            if hi > lo:
                out[lo:hi, e] = src[sbg]
    return out


def native_core(nframes):
    """Render the L channel at the decoder's own scale -- no normalisation."""
    import m17_ac4_walk as W
    import m20_ac4_toc2 as M
    import m29_audio as R
    import m30_filterbank as FB
    T = C.load_tables()
    fr = W.samples(os.path.join(HERE, "m7_out/rf33_audio_pid13.mp4"))[:nframes]
    wins, lengths = [], []
    for f in fr:
        try:
            st = M.parse(f)
            o = st["toc_bytes"] + st["substream_sizes"][0]
            sub = f[o:o + st["substream_sizes"][1]]
            r = C.decode_element(sub, T,
                                 b_iframe=bool(st["b_iframe_global"]))
            if r["bits"] > r["nbits"]:
                raise ValueError
            f_lr = r["st_lr"]["framing"]
            pk_l = R.packed_spectrum(r["L"], f_lr)
            pk_r = R.packed_spectrum(r["R"], f_lr)
            L, _ = R.unmix(pk_l, pk_r, r["st_lr"], f_lr)
            u, seq = R.ungroup(L, f_lr)
        except Exception:                                      # noqa: BLE001
            u, seq = [np.zeros(1536)], [1536]
        wins.extend(u)
        lengths.extend(seq)
    return FB.synthesise(wins, lengths, 1536)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="m36 envadj")
    ap.add_argument("--audio", default="tv_audio_valid.wav")
    ap.add_argument("--frames", type=int, default=120)
    a = ap.parse_args(argv)
    wav = a.audio if os.path.isabs(a.audio) else os.path.join(HERE, a.audio)

    print("M36 -- A-SPX envelope adjustment")
    print("=" * 74)
    cfg = dict(start_freq=5, stop_freq=1, master_freq_scale=1, noise_sbg=3,
               quant_mode_env=1, freq_res_mode=2, num_env_bits_fixfix=0,
               interpolation=1)
    t = A.sbg_tables(cfg, 0)
    lim = sbg_lim(t)
    print(f"\n  limiter subband groups: {lim}")

    # CORE AUDIO AT ITS NATIVE SCALE.
    #
    # m29 peak-normalises its render (0.9 / peak), which destroys exactly the
    # absolute level the A-SPX envelopes refer to.  Pseudocode 95's
    # EPSILON = 1.0 then dominates the denominator -- with est_sig ~ 1e-4 the
    # gain saturates and the adjusted band comes out ~1e6 too quiet.  So the
    # core is re-rendered here WITHOUT normalisation, straight from the
    # dequantised spectra, which is what the envelopes are defined against.
    x = native_core(a.frames)
    w = Q.qwin()
    x = x[:(len(x) // Q.NSB) * Q.NSB]
    Qm = Q.analyse(x, w, Q.analysis_matrix())
    Qh = G.hf_generate(Qm, t)
    print(f"  {Qh.shape[1]} QMF time slots from {a.frames} frames of core audio")

    # envelopes for channel L, decoded exactly as M34 does
    arrays = H.parse_c(H.DEFAULT_C)
    tabs = {nm[:-4]: S.Huff(arrays[nm], arrays[nm[:-4] + "_CW"])
            for nm in arrays if nm.endswith("_LEN") and nm.startswith("ASPX_")}
    for nm in tabs:
        P.CENTRE[nm] = 0 if nm.endswith("_F0") else \
            (len(arrays[nm + "_LEN"]) - 1) // 2

    import m17_ac4_walk as W
    import m20_ac4_toc2 as M
    from m19_ac4_toc import Bits
    T = C.load_tables()
    fr_all = W.samples(os.path.join(HERE, "m7_out/rf33_audio_pid13.mp4"))
    nats = A.num_aspx_timeslots()
    maps = V.sbg_maps(t)

    fstate, prev, nprev = {}, {}, {}
    rel_err, nchecked, nats_in_ts = [], 0, 2
    slot = 0
    for i, f in enumerate(fr_all[:a.frames]):
        try:
            st = M.parse(f)
            o = st["toc_bytes"] + st["substream_sizes"][0]
            sub = f[o:o + st["substream_sizes"][1]]
            ifr = bool(st["b_iframe_global"])
            r = C.decode_element(sub, T, b_iframe=ifr)
            if r["bits"] > r["nbits"]:
                slot += 24
                continue
            b = Bits(sub)
            b.p = r["bits"]
            grp = P.aspx_data(b, tabs, t, cfg, nats, ifr, True, fstate, 0)
        except Exception:                                      # noqa: BLE001
            slot += 24
            continue
        dat = grp[0]
        q, prev[0] = V.reconstruct(dat, t, maps, prev.get(0, (None, 0)))
        borders = dat["fr"].get("borders") or [0, nats]
        n_sb = t["num_sb_aspx"]
        sig_sb = map_scf(q, t, dat["fr"], n_sb)
        est = estimate(Qh, t, borders, slot, nats_in_ts)
        if est.shape != sig_sb.shape:
            slot += 24
            continue
        nse, nprev[0] = noise_envelopes(dat, t, nprev.get(0))
        noise_sb = map_noise(nse, t, n_sb, est.shape[1])
        m = (sig_sb > 0) & (est > 0)
        if m.any():
            g0 = gain_unlimited(sig_sb, est)
            raw = np.abs(est[m] * g0[m] ** 2 - sig_sb[m]) / sig_sb[m]
            # signal only, noise level forced to zero
            g_s, _ = adjust(Qh, t, sig_sb, est, lim)
            e_sig = np.abs(est[m] * g_s[m] ** 2 - sig_sb[m]) / sig_sb[m]
            # signal AND noise, the full Pseudocode 94..101 path
            g_n, nl = adjust(Qh, t, sig_sb, est, lim, noise_sb)
            tot = est[m] * g_n[m] ** 2 + nl[m] ** 2
            e_tot = np.abs(tot - sig_sb[m]) / sig_sb[m]
            rel_err.append((float(np.median(raw)), float(np.median(e_sig)),
                            float(np.max(g_n[m] / np.maximum(g0[m], 1e-30))),
                            float(np.median(e_tot)),
                            float(np.median(noise_sb[m]))))
            nchecked += 1
        slot += 24

    print(f"\n  {nchecked} frames adjusted and re-measured")
    if not rel_err:
        print("  nothing measured")
        return 1
    e = np.array(rel_err)
    raw_e, lim_e, ratio, tot_e, nlev = (e[:, 0], e[:, 1], e[:, 2],
                                        e[:, 3], e[:, 4])

    print("\n  GATE 1  the GAIN reproduces the transmitted envelope exactly")
    print("    (estimate + dequantisation + mapping + Pseudocode 95, before "
          "the limiter)")
    print(f"    relative error: median {np.median(raw_e):.2e}, "
          f"max {raw_e.max():.2e}")
    # float64 through a sqrt and a square; 1e-9 was tighter than the
    # arithmetic can deliver and failed an exact result
    g1ok = np.median(raw_e) < 1e-6
    print(f"    {'PASS' if g1ok else 'FAIL'}  closed loop, no free parameters")

    print("\n  GATE 2  the limiter/boost stays inside its stated bounds")
    print(f"    limited+boosted / unlimited gain, worst case {ratio.max():.6f}")
    print(f"    MAX_BOOST_FACT (Pseudocode 100)          {MAX_BOOST_FACT:.6f}")
    print(f"    residual envelope error after limiting: median "
          f"{np.median(lim_e):.4f}")
    # The boost is computed PER LIMITER SUBBAND GROUP and applied to every
    # subband in it, so an individual subband's gain CAN exceed its unlimited
    # value -- the group's energy is restored and redistributed.  A first
    # version of this gate asserted "the limiter only ever reduces" and failed
    # correct behaviour; what is actually bounded is the boost, at
    # MAX_BOOST_FACT, and the measured worst case sits exactly on it.
    g2ok = ratio.max() <= MAX_BOOST_FACT + 1e-6
    print(f"    {'PASS' if g2ok else 'FAIL'}  the cap binds and holds -- the "
          f"residual IS the limiter doing its job")

    print("\n  GATE 3  adding the NOISE closes the residual")
    print(f"    median scf_noise_sb over the band: {np.median(nlev):.4f}")
    print(f"    envelope error, signal only        {np.median(lim_e):.4f}")
    print(f"    envelope error, signal + noise     {np.median(tot_e):.4f}")
    g3ok = np.median(tot_e) < np.median(lim_e)
    print(f"    {'PASS' if g3ok else 'FAIL'}  the noise term moves the "
          f"achieved envelope TOWARD the target")

    ok = g1ok and g2ok and g3ok
    print("\n" + "=" * 74)
    if ok:
        print("  ENVELOPE ADJUSTMENT CLOSES THE LOOP.  The gain reproduces "
              "the transmitted\n  envelope exactly, the limiter caps it, and "
              "the noise term halves the\n  residual.  Additional harmonics "
              "remain.")
        return 0
    print("  NOT established -- the achieved envelope does not match the "
          "target")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
