#!/usr/bin/env python3
"""M37 -- render WITH the A-SPX high band: 0..21 kHz instead of 0..12 kHz.

Everything up to here produced either core audio band-limited at 12 kHz (m29)
or high-band stages gated in the QMF domain without writing a file (m35, m36).
This wires them together and writes audio that actually contains 12..21 kHz.

THE CHAIN
----------
    1. decode the MDCT core at NATIVE scale          m28 + m29 machinery
    2. QMF analyse                                    m33
    3. HF generate -- patch low subbands upward       m35
    4. envelope adjust + noise, per frame             m36
    5. QMF synthesise                                 m33
    6. write

Step 1 must be at native scale: m29's own renderer peak-normalises, which
destroys the absolute level the A-SPX envelopes are defined against and makes
the gain saturate (see M36).

THE GATES
----------
  1. energy APPEARS in 12..21 kHz, where the core render has none;
  2. energy does NOT appear above 21 kHz -- A-SPX stops at QMF subband 56, so
     a correct high band is bounded there and a wrong one is not;
  3. below 12 kHz the output still matches the core render closely -- the
     high-band work must not disturb the part that was already verified.

Gate 3 matters: it is easy to produce a "fuller" file by breaking the core.

ALIGNMENT -- MEASURED, NOT ASSUMED
-----------------------------------
The control data for frame i has to be applied to the QMF slots that actually
carry frame i's audio, and that is two delays deep.  Both were measured with
impulses rather than reasoned about:

    MDCT render delay   0 frames.  A spectrum placed in frame 6 alone peaks in
                        output frame 6 (spreading into 7 through the overlap),
                        so the render is frame-aligned.
    QMF analysis delay  4.0 slots.  An impulse at input slot 10 has its
                        subband-energy centroid at slot 14.0.

So TS_OFFSET = 4.  The first version of this file used 0 and smeared the high
band ~4 slots (~5.3 ms) early against the core.  Two statistical attempts to
find the offset by correlating the transmitted envelope against the core's
energy both failed first -- one gave a broad monotone curve with no peak, the
other (differencing to sharpen it) gave noise, because the transmitted envelope
is piecewise CONSTANT and differencing annihilates it.  A deterministic impulse
measurement settled in one step what two statistical ones could not.

Usage:
    python m37_render_hf.py [--frames 600] [--out tv_audio_hf.wav]
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m17_ac4_walk as W                                          # noqa: E402
import m20_ac4_toc2 as M                                          # noqa: E402
import m23_hcb as H                                               # noqa: E402
import m24_spectral as S                                          # noqa: E402
import m28_channels as C                                          # noqa: E402
import m29_audio as R                                             # noqa: E402
import m30_filterbank as FB                                       # noqa: E402
import m31_aspx_bands as A                                        # noqa: E402
import m32_aspx_parse as P                                        # noqa: E402
import m33_qmf as Q                                               # noqa: E402
import m34_aspx_env as V                                          # noqa: E402
import m35_hfgen as G                                             # noqa: E402
import m36_envadj as E                                            # noqa: E402

FS = 48000

# QMF ANALYSIS DELAY, MEASURED not assumed: an impulse at input slot 10 has its
# subband-energy centroid at slot 14.0, so frame i's audio occupies QMF slots
# [i*24 + TS_OFFSET, ...).  The MDCT render itself was measured to add ZERO
# frames of delay (frame 6's spectrum peaks in output frame 6), so this is the
# whole correction.  Applying the control at offset 0 -- the previous default --
# smears the high band ~4 slots (~5.3 ms) early against the core.
TS_OFFSET = 4


def decode_all(nframes, start=0):
    """-> native-scale L and R PCM, plus the per-frame A-SPX group data."""
    T = C.load_tables()
    arrays = H.parse_c(H.DEFAULT_C)
    tabs = {n[:-4]: S.Huff(arrays[n], arrays[n[:-4] + "_CW"])
            for n in arrays if n.endswith("_LEN") and n.startswith("ASPX_")}
    for nm in tabs:
        P.CENTRE[nm] = 0 if nm.endswith("_F0") else \
            (len(arrays[nm + "_LEN"]) - 1) // 2

    # The capture directory is a parameter, not a constant: the same chain
    # has to run on a fresh capture, not only on the 8/06 bank.
    fr = W.samples(os.path.join(HERE, os.environ.get(
        "AC4_SRC", "m7_out/rf33_audio_pid13.mp4")))
    if nframes:
        fr = fr[:nframes]
    # A live decoder resumes from an I-FRAME, where both the MDCT overlap
    # and the A-SPX delta chain reset. Slicing here keeps every downstream
    # index frame-relative, so callers do not have to think about it.
    if start:
        fr = fr[start:]
    nats = A.num_aspx_timeslots()
    cfg, fstate = None, {}
    wins = {"L": [], "R": []}
    lengths = []
    groups = []
    for f in fr:
        got = None
        try:
            st = M.parse(f)
            o = st["toc_bytes"] + st["substream_sizes"][0]
            sub = f[o:o + st["substream_sizes"][1]]
            ifr = bool(st["b_iframe_global"])
            r = C.decode_element(sub, T, b_iframe=ifr)
            if r["bits"] > r["nbits"]:
                raise ValueError
            if r.get("aspx"):
                cfg = r["aspx"]
            fl = r["st_lr"]["framing"]
            pl = R.packed_spectrum(r["L"], fl)
            pr = R.packed_spectrum(r["R"], fl)
            L, Rr = R.unmix(pl, pr, r["st_lr"], fl)
            uL, seq = R.ungroup(L, fl)
            uR, _ = R.ungroup(Rr, fl)
            from m19_ac4_toc import Bits
            b = Bits(sub)
            b.p = r["bits"]
            got = P.aspx_data(b, tabs, t_for(cfg), cfg, nats, ifr, True,
                              fstate, 0) if cfg else None
        except Exception:                                      # noqa: BLE001
            uL = uR = [np.zeros(1536)]
            seq = [1536]
        wins["L"].extend(uL)
        wins["R"].extend(uR if len(uR) == len(seq) else
                         [np.zeros(x) for x in seq])
        lengths.extend(seq)
        groups.append(got)
    pcm = {k: FB.synthesise(wins[k], lengths, 1536) for k in ("L", "R")}
    return pcm, groups, cfg


_T_CACHE = {}


def t_for(cfg):
    key = (cfg["start_freq"], cfg["stop_freq"], cfg["master_freq_scale"],
           cfg["noise_sbg"])
    if key not in _T_CACHE:
        _T_CACHE[key] = A.sbg_tables(cfg, 0)
    return _T_CACHE[key]


def apply_hf_pair(pcmL, pcmR, groups, cfg):
    """Steps 2-5 for the L/R pair, JOINTLY.

    The pair cannot be processed one channel at a time: when aspx_balance = 1
    -- 89.7 % of frames here -- channel B carries a PAN RATIO and its scale
    factors only become levels through Pseudocode 84, which needs both
    channels' indices at once.  Treating B's values as levels (Pseudocode 82,
    which explicitly covers only balance = 0) made the right channel 20x too
    loud and uncorrelated with its own core, while the left channel looked
    perfect -- so a single-channel gate could never have caught it.
    """
    t = t_for(cfg)
    lim = E.sbg_lim(t)
    nt = E.noise_table()
    w = Q.qwin()
    maps = V.sbg_maps(t)
    n_sb, sbx = t["num_sb_aspx"], t["sbx"]

    out = {}
    Qh = {}
    for k, x in (("L", pcmL), ("R", pcmR)):
        xx = x[:(len(x) // Q.NSB) * Q.NSB]
        Qh[k] = G.hf_generate(Q.analyse(xx, w, Q.analysis_matrix()), t)

    prev = {0: (None, 0), 1: (None, 0)}
    nprev = {0: None, 1: None}
    idx_noise = 0
    slot = 0
    applied = 0
    for got in groups:
        if got is None or len(got) < 2:
            slot += 24
            continue
        try:
            qa, prev[0] = V.reconstruct(got[0], t, maps, prev[0])
            qb, prev[1] = V.reconstruct(got[1], t, maps, prev[1])
            na, nprev[0] = E.noise_qscf(got[0], t, nprev[0])
            nb, nprev[1] = E.noise_qscf(got[1], t, nprev[1])
            fr0 = got[0]["fr"]
            borders = fr0.get("borders") or [0, 12]
            if got[0].get("balance") == 1:
                sigA, sigB = E.map_scf_joint(qa, qb, t, fr0, n_sb)
                nsA, nsB = E.map_noise_joint(na, nb, t, n_sb, sigA.shape[1])
            else:
                sigA = E.map_scf(qa, t, fr0, n_sb)
                sigB = E.map_scf(qb, t, got[1]["fr"], n_sb)
                dqa = [[2.0 ** (E.NOISE_FLOOR_OFFSET - v) for v in e]
                       for e in na]
                dqb = [[2.0 ** (E.NOISE_FLOOR_OFFSET - v) for v in e]
                       for e in nb]
                nsA = E.map_noise(dqa, t, n_sb, sigA.shape[1])
                nsB = E.map_noise(dqb, t, n_sb, sigB.shape[1])
        except Exception:                                      # noqa: BLE001
            slot += 24
            continue
        # WITH balance = 0 THE TWO CHANNELS HAVE INDEPENDENT FRAMINGS, so each
        # gets its OWN borders.  Using channel 0's for both silently skipped
        # 114 of 3480 frames (they failed a shape check and fell through to
        # core-only), which the per-frame data diagnostic could not see because
        # the data was fine -- the renderer was wrong.
        b0 = borders
        b1 = (borders if got[0].get("balance") == 1
              else (got[1]["fr"].get("borders") or borders))
        ok = True
        for k, sig, ns, bd in (("L", sigA, nsA, b0), ("R", sigB, nsB, b1)):
            est = E.estimate(Qh[k], t, bd, slot, 2, ts_off=TS_OFFSET)
            if est.shape != sig.shape:
                ok = False
                break
            g, nl = E.adjust(Qh[k], t, sig, est, lim, ns[:, :est.shape[1]])
            for e in range(est.shape[1]):
                ta = slot + int(bd[e]) * 2 + TS_OFFSET
                tz = slot + int(bd[e + 1]) * 2 + TS_OFFSET
                ta, tz = max(ta, 0), min(tz, Qh[k].shape[1])
                if tz <= ta:
                    continue
                Qh[k][sbx:sbx + n_sb, ta:tz] *= g[:, e:e + 1]
                base = idx_noise if k == "L" else idx_noise
                for ts in range(ta, tz):
                    kk = (base + n_sb * (ts - ta) + np.arange(n_sb) + 1) % 512
                    Qh[k][sbx:sbx + n_sb, ts] += nl[:, e] * nt[kk]
        if ok:
            applied += 1
            span = sum(max(0, min(slot + int(borders[e + 1]) * 2 + TS_OFFSET,
                                  Qh["L"].shape[1])
                           - max(slot + int(borders[e]) * 2 + TS_OFFSET, 0))
                       for e in range(len(borders) - 1))
            idx_noise = (idx_noise + n_sb * span) % 512
        slot += 24

    for k in ("L", "R"):
        out[k] = Q.synthesise(Qh[k], w, Q.synthesis_matrix(255))
    return out, applied


def band(x, lo, hi):
    X = np.abs(np.fft.rfft(x * np.hanning(len(x))))
    f = np.fft.rfftfreq(len(x), 1.0 / FS)
    m = (f >= lo) & (f < hi)
    return float((X[m] ** 2).sum()), float((X ** 2).sum())


def main(argv=None):
    ap = argparse.ArgumentParser(prog="m37 render hf")
    ap.add_argument("--src", default=None,
                    help="AC-4 asset (default m7_out/rf33_audio_pid13.mp4)")
    ap.add_argument("--frames", type=int, default=600)
    ap.add_argument("--out", default="tv_audio_hf.wav")
    a = ap.parse_args(argv)
    if a.src:
        os.environ["AC4_SRC"] = a.src

    print("M37 -- render with the A-SPX high band")
    print("=" * 74)
    pcm, groups, cfg = decode_all(a.frames)
    print(f"  decoded {a.frames} frames, core at native scale")
    if cfg is None:
        print("  no aspx_config")
        return 1

    out, n = apply_hf_pair(pcm["L"], pcm["R"], groups, cfg)
    print(f"  A-SPX applied jointly to {n} frames, {len(out['L'])} samples")

    core = pcm["L"][:len(out["L"])]
    hf = out["L"]
    lo_c, tot_c = band(core, 20, 12000)
    hi_c, _ = band(core, 12000, 21000)
    top_c, _ = band(core, 21000, 24000)
    lo_h, tot_h = band(hf, 20, 12000)
    hi_h, _ = band(hf, 12000, 21000)
    top_h, _ = band(hf, 21000, 24000)

    print(f"\n  GATE 1  the high band APPEARS")
    print(f"    12..21 kHz   core {100 * hi_c / tot_c:8.4f} %   "
          f"with A-SPX {100 * hi_h / tot_h:8.4f} %")
    g1 = (hi_h / tot_h) > 100 * (hi_c / tot_c) and (hi_h / tot_h) > 1e-4
    print(f"    {'PASS' if g1 else 'FAIL'}")

    print(f"\n  GATE 2  and STOPS at the A-SPX stop frequency")
    print(f"    21..24 kHz   core {100 * top_c / tot_c:8.4f} %   "
          f"with A-SPX {100 * top_h / tot_h:8.4f} %")
    g2 = (top_h / tot_h) < 0.01 * (hi_h / tot_h)
    print(f"    {'PASS' if g2 else 'FAIL'}  bounded at 21 kHz, as A-SPX is")

    print(f"\n  GATE 3  the core band is NOT disturbed")
    # THE QMF ROUND TRIP DELAYS BY 577 SAMPLES (measured in M33).  Comparing
    # core[:n] against hf[:n] ignores that and reads a correlation of ~0 on a
    # render that is actually fine.  Same alignment mistake as M33's SNR gate;
    # the delay is a property of the filterbank, not an unknown, so it is
    # applied rather than searched for.
    DELAY = 577
    n = min(len(core) - DELAY, len(hf) - DELAY)
    r = float(np.corrcoef(core[:n], hf[DELAY:DELAY + n])[0, 1])
    print(f"    correlation with the core, aligned by the QMF's {DELAY}-sample "
          f"delay: {r:+.4f}")
    print(f"    below 12 kHz  core {100 * lo_c / tot_c:.2f} %  "
          f"with A-SPX {100 * lo_h / tot_h:.2f} %")
    g3 = r > 0.95
    print(f"    {'PASS' if g3 else 'FAIL'}  the high band was added, not "
          f"traded for the core")

    print("\n  GATE 4  the two channels are still a stereo PAIR")
    nn = min(len(out["L"]), len(out["R"]))
    r_lr = float(np.corrcoef(out["L"][:nn], out["R"][:nn])[0, 1])
    rms_l = float(np.sqrt((out["L"] ** 2).mean()))
    rms_r = float(np.sqrt((out["R"] ** 2).mean()))
    bal = rms_r / max(rms_l, 1e-30)
    print(f"    L/R correlation {r_lr:+.4f}   rms ratio R/L {bal:.3f}")
    # This gate exists because a single-channel check CANNOT catch the
    # balance-coding bug: with Pseudocode 84 missing, L was perfect (+0.9999
    # against its core) while R was 20x too loud and uncorrelated.  Gates 1-3
    # all passed on that render.
    g4 = r_lr > 0.5 and 0.5 < bal < 2.0
    print(f"    {'PASS' if g4 else 'FAIL'}  correlated and level-matched "
          f"(a broken balance decode shows up here and nowhere else)")

    peak = max(np.abs(out[k]).max() for k in out)
    g = 0.9 / peak if peak > 0 else 1.0
    outp = a.out if os.path.isabs(a.out) else os.path.join(HERE, a.out)
    R.write_wav(outp, [np.clip(out["L"] * g, -1, 1),
                       np.clip(out["R"] * g, -1, 1)])
    print(f"\n  wrote {os.path.basename(outp)}  "
          f"{len(out['L']) / FS:.1f} s stereo, 0..21 kHz")

    print("\n" + "=" * 74)
    if g1 and g2 and g3 and g4:
        print("  HIGH BAND RENDERED.  The file now carries 12..21 kHz from "
              "A-SPX,\n  bounded where A-SPX stops, with the core intact.")
        return 0
    print("  NOT established -- see the failing gate")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
