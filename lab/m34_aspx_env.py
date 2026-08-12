#!/usr/bin/env python3
"""M34 -- decoding the A-SPX signal envelopes, not just counting their bits.

M32 walks the A-SPX payload and proves it consumes exactly `audio_size`.  But
proving the bit COUNT is right is not the same as having the VALUES: the
envelopes are what actually shape the high band, and until they are
reconstructed nothing above 12 kHz can be synthesised.

PSEUDOCODE 80 -- RECONSTRUCTION
--------------------------------
    delta = 2 if (ch == 1 and aspx_balance == 1) else 1

    FREQ (aspx_sig_delta_dir == 0):
        qscf[sbg][atsg] = delta * sum(aspx_data_sig[atsg][0..sbg])
    TIME:
        qscf[sbg][atsg] = qscf_prev[sbg] + delta * aspx_data_sig[atsg][sbg]

where the TIME branch remaps sbg through sbg_idx_low2high / high2low when the
frequency resolution CHANGES between consecutive envelopes -- so a resolution
switch does not corrupt the running value.

**This is stateful across frames.**  `qscf_sig_sbg_prev` is the last envelope of
the PREVIOUS A-SPX interval, so a TIME-coded first envelope depends on the
frame before it.  Second structure in this decoder that cannot be evaluated on
one frame in isolation (the framing borders were the first).

THE GATE
---------
Envelopes are LEVELS, quantised in 1.5 or 3.0 dB steps, so the same test that
proved the MDCT scale factors applies: loudness is continuous in time, so the
mean envelope of frame N must track frame N+1, with a shuffled control at zero.
Random bits pushed through a complete prefix code produce plausible-looking
values; they do not produce time structure.

Second gate: envelope level should FALL with frequency across the A-SPX range.
That is a property of real programme material, not of the arithmetic.

Usage:
    python m34_aspx_env.py
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
import m31_aspx_bands as A                                        # noqa: E402
import m32_aspx_parse as P                                        # noqa: E402
from m19_ac4_toc import Bits                                      # noqa: E402


def sbg_maps(t):
    """Pseudocode 80's index maps between the high- and low-res tables."""
    hi, lo = t["hi"], t["lo"]
    high2low = [0] * t["n_hi"]
    low2high = [0] * (t["n_lo"] + 1)
    sbg_low = 0
    for sbg in range(t["n_hi"]):
        if sbg_low + 1 < len(lo) and lo[sbg_low + 1] == hi[sbg]:
            sbg_low += 1
            low2high[sbg_low] = sbg
        high2low[sbg] = sbg_low
    return high2low, low2high


def reconstruct(dat, t, maps, prev):
    """Pseudocode 80.  -> qscf[atsg][sbg], and the new `prev` state."""
    high2low, low2high = maps
    fr = dat["fr"]
    delta = 2 if (dat["ch"] == 1 and dat["balance"] == 1) else 1
    out = []
    prev_q, prev_res = prev
    for atsg, vals in enumerate(dat["sig"]):
        res = fr["freq_res"][atsg]
        n_sbg = t["n_hi"] if res else t["n_lo"]
        q = [0] * n_sbg
        if dat["dir_sig"][atsg] == 0:                     # FREQ
            acc = 0
            for sbg in range(min(n_sbg, len(vals))):
                acc += delta * vals[sbg]
                q[sbg] = acc
        else:                                             # TIME
            for sbg in range(min(n_sbg, len(vals))):
                if prev_q is None:
                    base = 0
                elif res == prev_res:
                    base = prev_q[sbg] if sbg < len(prev_q) else 0
                elif res == 0 and prev_res == 1:
                    i = low2high[sbg] if sbg < len(low2high) else 0
                    base = prev_q[i] if i < len(prev_q) else 0
                else:
                    i = high2low[sbg] if sbg < len(high2low) else 0
                    base = prev_q[i] if i < len(prev_q) else 0
                q[sbg] = base + delta * vals[sbg]
        out.append(q)
        prev_q, prev_res = q, res
    return out, (prev_q, prev_res)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="m34 aspx env")
    ap.add_argument("path", nargs="?",
                    default="m7_out/rf33_audio_pid13.mp4")
    a = ap.parse_args(argv)
    p = a.path if os.path.isabs(a.path) else os.path.join(HERE, a.path)

    print("M34 -- A-SPX signal envelopes, decoded")
    print("=" * 74)
    T = C.load_tables()
    arrays = H.parse_c(H.DEFAULT_C)
    tabs = {n[:-4]: S.Huff(arrays[n], arrays[n[:-4] + "_CW"])
            for n in arrays if n.endswith("_LEN") and n.startswith("ASPX_")}
    for nm in tabs:
        P.CENTRE[nm] = 0 if nm.endswith("_F0") else \
            (len(arrays[nm + "_LEN"]) - 1) // 2

    fr = W.samples(p)
    nats = A.num_aspx_timeslots()
    cfg = None
    fstate = {}
    prev = {}
    ix, mean_lvl, slopes = [], [], []
    ok = 0
    for i, f in enumerate(fr):
        try:
            st = M.parse(f)
            o = st["toc_bytes"] + st["substream_sizes"][0]
            sub = f[o:o + st["substream_sizes"][1]]
            ifr = bool(st["b_iframe_global"])
            r = C.decode_element(sub, T, b_iframe=ifr)
            if r["bits"] > r["nbits"]:
                continue
        except Exception:                                      # noqa: BLE001
            continue
        if r.get("aspx"):
            cfg = r["aspx"]
        if cfg is None:
            continue
        t = A.sbg_tables(cfg, 0)
        maps = sbg_maps(t)
        b = Bits(sub)
        b.p = r["bits"]
        try:
            groups = []
            groups += P.aspx_data(b, tabs, t, cfg, nats, ifr, True, fstate, 0)
            groups += P.aspx_data(b, tabs, t, cfg, nats, ifr, True, fstate, 1)
            groups += P.aspx_data(b, tabs, t, cfg, nats, ifr, False, fstate, 2)
        except Exception:                                      # noqa: BLE001
            continue
        lv = []
        for gi, dat in enumerate(groups):
            key = (gi,)
            q, prev[key] = reconstruct(dat, t, maps, prev.get(key, (None, 0)))
            for env in q:
                if env:
                    lv.append(float(np.mean(env)))
                    if len(env) > 2:
                        x = np.arange(len(env))
                        slopes.append(float(np.polyfit(x, env, 1)[0]))
        if not lv:
            continue
        ok += 1
        ix.append(i)
        mean_lvl.append(float(np.mean(lv)))

    print(f"\n  {ok} frames with reconstructed envelopes")
    if ok < 200:
        print("  too few to gate")
        return 1
    E = np.array(mean_lvl)
    I = np.array(ix)
    print(f"  mean envelope level: min {E.min():.1f}  max {E.max():.1f}  "
          f"mean {E.mean():.1f}  (quantiser step "
          f"{'3.0' if cfg['quant_mode_env'] else '1.5'} dB)")

    adj = np.array([j for j in range(len(I) - 1) if I[j + 1] == I[j] + 1])
    rng = np.random.default_rng(0)
    q = rng.permutation(len(E))
    r_adj = float(np.corrcoef(E[adj], E[adj + 1])[0, 1])
    r_shuf = float(np.corrcoef(E[q[:-1]], E[q[1:]])[0, 1])
    se = 1.0 / np.sqrt(max(len(adj), 2))
    print(f"\n  GATE 1  envelopes are levels, so they are continuous in time")
    print(f"    adjacent frames   r = {r_adj:+.4f}   ({len(adj)} pairs)")
    print(f"    shuffled control  r = {r_shuf:+.4f}   (se {se:.3f})")
    g1 = r_adj > 0.15 and abs(r_adj) > 3 * max(abs(r_shuf), se)
    print(f"    {'PASS' if g1 else 'FAIL'}")

    sl = np.array(slopes)
    print(f"\n  GATE 2  high-band level falls with frequency in real material")
    print(f"    mean slope across the A-SPX range {sl.mean():+.3f} "
          f"steps/subband-group  ({len(sl)} envelopes)")
    g2 = sl.mean() < 0
    print(f"    {'PASS' if g2 else 'FAIL'}  negative slope")

    print("\n" + "=" * 74)
    if g1 and g2:
        print("  A-SPX ENVELOPES DECODED.  They carry time structure and fall "
              "with\n  frequency -- what the high band will be shaped by.")
        return 0
    print("  NOT established")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
