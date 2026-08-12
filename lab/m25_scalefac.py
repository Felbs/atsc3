#!/usr/bin/env python3
"""M25 -- the LFE's scale factors: the exponents that make the lines real.

M24 decoded the quantised spectral lines -- integers.  They are not amplitudes
yet: each scale factor band carries an exponent, and the reconstruction is

    coefficient = sign(q) * |q|^(4/3) * 2^((scale_factor - offset)/4)

so without the scale factors the spectrum has the right shape and the wrong
loudness, per band.  This reads them.

THE SYNTAX (4.2.8.5)
---------------------
    reference_scale_factor         8 bits
    first_scf_found = 0
    for each band with sfb_cb != 0 AND max_quant_idx > 0:
        if first_scf_found:  dpcm_sf = huff_decode(ASF_HCB_SCALEFAC)
        else:                first_scf_found = 1     # this band USES the reference

So the first populated band spends no bits at all -- it takes the 8-bit
reference -- and every later band is a Huffman-coded DELTA from the running
value.  `ASF_HCB_SCALEFAC` has 121 entries, which is 2*60+1, so the delta is
`index - 60`.  That centring is not a guess: a DPCM alphabet of odd size is
symmetric about its middle, and 121 admits exactly one centre.

WHAT MAKES A DECODE BELIEVABLE HERE
-------------------------------------
The same trap as M24 applies -- ASF_HCB_SCALEFAC is a complete prefix code, so
any bits decode to some delta, and the delta histogram will look like the
codebook whether or not the bits are ours.  Static shape proves nothing.

Two things random bits cannot fake:

  * **Scale factors are loudness, and loudness is continuous in time.**  The
    reference scale factor of frame N should track frame N+1, with a shuffled
    control at zero.  This is M24's test applied to a different quantity, and
    it is stronger here because the reference is read as a plain 8-bit field
    rather than through a codebook -- so a bit-offset error moves it wildly.
  * **Scale factors and quantised magnitudes are coupled.**  A band whose lines
    are large does not usually also carry a large exponent; the encoder trades
    them off to hit a bitrate.  A non-zero rank correlation between the two is
    evidence they came from the same encoder decision, not from noise.

Usage:
    python m25_scalefac.py
"""
from __future__ import annotations

import argparse
import collections
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m17_ac4_walk as W                                         # noqa: E402
import m20_ac4_toc2 as M                                         # noqa: E402
import m23_hcb as H                                              # noqa: E402
import m24_spectral as S                                         # noqa: E402
from m19_ac4_toc import Bits                                     # noqa: E402

SF_CENTRE = 60          # ASF_HCB_SCALEFAC is 121 = 2*60+1 entries


def decode_lfe_full(sub, tables, sf_table):
    """-> (lines, scale_factors, bit position after asf_scalefac_data)."""
    b = Bits(sub)
    b.u(15)
    if b.u(1):
        b.vb(7)
    b.u(3)                                   # 5_X_codec_mode
    max_sfb = b.u(3)
    sect_cb = b.u(4)
    b.u(5)                                   # sect_len_incr
    if sect_cb != 11:
        raise ValueError("not codebook 11")
    hb = tables[sect_cb]
    mod, off = S.CB_MOD[sect_cb], S.CB_OFF[sect_cb]
    end = S.SFB_OFFSET[min(max_sfb, 3)]

    lines = []
    while len(lines) < end:
        idx = hb.decode(b)
        vals = [idx // mod - off, idx % mod - off]
        vals = [(-v if (v and b.u(1)) else v) for v in vals]
        for j, v in enumerate(vals):
            if abs(v) == 16:
                n = 0
                while b.u(1):
                    n += 1
                mag = (1 << (n + 4)) + b.u(n + 4)
                vals[j] = -mag if v < 0 else mag
        lines.extend(vals)
    lines = np.array(lines[:end])

    # --- asf_scalefac_data -------------------------------------------------
    ref = b.u(8)
    sfs, cur, first = [], ref, False
    for sfb in range(min(max_sfb, 3)):
        band = lines[S.SFB_OFFSET[sfb]:S.SFB_OFFSET[sfb + 1]]
        if sect_cb != 0 and np.abs(band).max() > 0:
            if first:
                cur += sf_table.decode(b) - SF_CENTRE
            else:
                first = True
            sfs.append(cur)
        else:
            sfs.append(None)
    return lines, ref, sfs, b.p


def main(argv=None):
    ap = argparse.ArgumentParser(prog="m25 scalefac")
    ap.add_argument("path", nargs="?", default="m7_out/rf33_audio_pid13.mp4")
    a = ap.parse_args(argv)
    p = a.path if os.path.isabs(a.path) else os.path.join(HERE, a.path)
    arrays = H.parse_c(H.DEFAULT_C)
    tables = {n: S.Huff(arrays[f"ASF_HCB_{n}_LEN"], arrays[f"ASF_HCB_{n}_CW"])
              for n in S.CB_MOD}
    sf_table = S.Huff(arrays["ASF_HCB_SCALEFAC_LEN"],
                      arrays["ASF_HCB_SCALEFAC_CW"])

    print("M25 -- the LFE's scale factors")
    print("=" * 72)
    fr = W.samples(p)
    refs, ix, mags, allsf = [], [], [], []
    ok = err = 0
    first_err = None
    for i, f in enumerate(fr):
        try:
            st = M.parse(f)
            if st["b_iframe_global"]:
                continue
            off = st["toc_bytes"] + st["substream_sizes"][0]
            sub = f[off:off + st["substream_sizes"][1]]
            lines, ref, sfs, _ = decode_lfe_full(sub, tables, sf_table)
            ok += 1
            refs.append(ref)
            ix.append(i)
            mags.append(float(np.abs(lines).mean()))
            allsf.extend([s for s in sfs if s is not None])
        except Exception as e:                                 # noqa: BLE001
            err += 1
            if first_err is None:
                first_err = (i, f"{type(e).__name__}: {e}")

    print(f"  decoded {ok} frames, {err} skipped/failed")
    if first_err:
        print(f"    first: frame {first_err[0]}: {first_err[1]}")
    if ok < 100:
        return 1
    R = np.array(refs, float)
    Mg = np.array(mags)
    I = np.array(ix)
    print(f"\n  reference_scale_factor  min {R.min():.0f} max {R.max():.0f} "
          f"mean {R.mean():.1f}  ({len(set(refs))} distinct)")
    sf = np.array(allsf, float)
    print(f"  all band scale factors  min {sf.min():.0f} max {sf.max():.0f} "
          f"mean {sf.mean():.1f}")

    adj = np.array([j for j in range(len(I) - 1) if I[j + 1] == I[j] + 1])
    r_adj = float(np.corrcoef(R[adj], R[adj + 1])[0, 1])
    rng = np.random.default_rng(0)
    q = rng.permutation(len(R))
    r_shuf = float(np.corrcoef(R[q[:-1]], R[q[1:]])[0, 1])
    print(f"\n  GATE 1  reference_scale_factor is loudness, so it must be "
          f"continuous in time")
    print(f"    adjacent frames   r = {r_adj:+.4f}   ({len(adj)} pairs)")
    print(f"    shuffled control  r = {r_shuf:+.4f}")

    rank = lambda v: np.argsort(np.argsort(v))
    r_couple = float(np.corrcoef(rank(R), rank(Mg))[0, 1])
    print(f"\n  GATE 2  scale factor and quantised magnitude are traded off by "
          f"the encoder")
    print(f"    rank correlation  r = {r_couple:+.4f}")

    print("\n" + "=" * 72)
    good = r_adj > 0.3 and abs(r_shuf) < 0.05 and abs(r_couple) > 0.1
    print("  SCALE FACTORS ARE REAL -- continuous in time and coupled to the "
          "lines." if good else
          "  NOT established")
    return 0 if good else 1


if __name__ == "__main__":
    raise SystemExit(main())
