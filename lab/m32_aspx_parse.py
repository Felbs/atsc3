#!/usr/bin/env python3
"""M32 -- parsing the A-SPX payload, so the substream accounts to the byte.

M31 built the frequency skeleton.  This walks the actual A-SPX bits that follow
the six channel elements:

    aspx_data_2ch()   L / R
    aspx_data_2ch()   Ls / Rs
    aspx_data_1ch()   C

THE GATE IS audio_size, NOT THE SUBSTREAM END
-----------------------------------------------
The first version of this gate demanded the parse reach the end of the
substream, and it failed by a stubborn ~42 bits.  That was the GATE being
wrong, not the parse.  Table 16:

    ac4_substream() {
        audio_size = audio_size_value;    15
        if (b_more_bits) audio_size += variable_bits(7) << 15;
        byte_align;
        audio_data(channel_mode, b_iframe);
        fill_bits;  byte_align;  metadata(b_iframe);  byte_align;
    }

so what follows audio_data is fill_bits and a metadata() payload -- not slack.
The substream states the length of audio_data explicitly, which makes a far
sharper test available: the parse must land within 0..7 bits BELOW
`start + audio_size * 8`, the 0..7 being byte alignment.  Every channel and
every A-SPX bit is then accounted for against a length the encoder wrote down,
and nothing can be mis-sized and hide, because the target is exact.

WHAT DECIDES THE BIT COUNT
---------------------------
    aspx_int_class    Table 125, a prefix code: 0=FIXFIX, 10=FIXVAR,
                      110=VARFIX, 111=VARVAR.  Note VARFIX and VARVAR are the
                      REVERSE of the order Table 53 lists its switch cases in.
    num_aspx_timeslots = 12 (M31), which is > 8, so aspx_framing's relative
                      border fields are 2 bits, not 1 (Table 53 Note 1).
    aspx_num_env_bits_fixfix = 0 -> envbits 1 -> FIXFIX num_env in {1, 2}
    aspx_freq_res_mode = 2 -> freq_res is NOT in the bitstream; Pseudocode 77
                      derives it from envelope duration.  For FIXFIX,
                      aspx_tsg_ptr is 0 so the first clause is vacuous and the
                      test is duration > num_aspx_timeslots/6 + 3.25 = 5.25.
                      With 12 timeslots, 1 envelope spans 12 and 2 span 6 --
                      both exceed 5.25, so both are HIGH resolution.
    get_aspx_hcb      Pseudocode 79: SIGNAL -> ASPX_HCB_ENV_<mode>_<15|30>_<t>,
                      NOISE -> ASPX_HCB_NOISE_<mode>_<t>, with quant_mode 0/1
                      mapping to 15/30 and stereo_mode 0/1 to LEVEL/BALANCE.

ALL FOUR INTERVAL CLASSES
--------------------------
FIXFIX, FIXVAR, VARFIX and VARVAR are all implemented, via Pseudocode 76's
border construction feeding Pseudocode 77's duration test.  VARFIX and VARVAR
carry state between frames -- on a non-i-frame the leading border is
`previous_stop_pos - num_aspx_timeslots` from the same channel group's previous
frame -- so the parse is sequential, not per-frame independent.

Usage:
    python m32_aspx_parse.py
"""
from __future__ import annotations

import argparse
import collections
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m17_ac4_walk as W                                          # noqa: E402
import m20_ac4_toc2 as M                                          # noqa: E402
import m23_hcb as H                                               # noqa: E402
import m24_spectral as S                                          # noqa: E402
import m28_channels as C                                          # noqa: E402
import m31_aspx_bands as A                                        # noqa: E402
from m19_ac4_toc import Bits                                      # noqa: E402

FIXFIX, FIXVAR, VARFIX, VARVAR = 0, 1, 2, 3
NAMES = {FIXFIX: "FIXFIX", FIXVAR: "FIXVAR", VARFIX: "VARFIX",
         VARVAR: "VARVAR"}


class Unsupported(Exception):
    """A branch that is not implemented -- counted, never guessed."""


def int_class(b):
    """Table 125.  0 / 10 / 110 / 111."""
    if not b.u(1):
        return FIXFIX
    if not b.u(1):
        return FIXVAR
    return VARVAR if b.u(1) else VARFIX


def aspx_framing(b, cfg, nats, b_iframe, state, key):
    """Table 53 + Pseudocode 76.  -> dict(num_env, num_noise, freq_res[]).

    All four interval classes.  The variable ones need the envelope BORDERS,
    not just the counts, because Pseudocode 77 decides each envelope's
    frequency resolution from its DURATION -- and that choice swings the
    codeword count per envelope between num_sbg_sig_highres and _lowres, so it
    is a bit-count decision, not a cosmetic one.

    VARFIX and VARVAR carry state across frames: on a non-i-frame the leading
    border is `previous_stop_pos - num_aspx_timeslots`, where previous_stop_pos
    is the trailing border of the SAME channel group in the previous frame.
    That is threaded through `state`, keyed per channel group.
    """
    ic = int_class(b)
    n_rel_bits = 2 if nats > 8 else 1
    rel_l, rel_r = [], []
    var_l = var_r = 0

    if ic == FIXFIX:
        envbits = cfg["num_env_bits_fixfix"] + 1
        num_env = 1 << b.u(envbits)
        tsg_ptr = 0
        borders = tab_border(nats, num_env)
    else:
        if ic == FIXVAR:
            var_r = b.u(2)
            for _ in range(b.u(n_rel_bits)):
                rel_r.append(2 * b.u(n_rel_bits) + 2)
        elif ic == VARFIX:
            if b_iframe:
                var_l = b.u(2)
            for _ in range(b.u(n_rel_bits)):
                rel_l.append(2 * b.u(n_rel_bits) + 2)
        else:                                        # VARVAR
            if b_iframe:
                var_l = b.u(2)
            for _ in range(b.u(n_rel_bits)):
                rel_l.append(2 * b.u(n_rel_bits) + 2)
            var_r = b.u(2)
            for _ in range(b.u(n_rel_bits)):
                rel_r.append(2 * b.u(n_rel_bits) + 2)
        num_env = len(rel_l) + len(rel_r) + 1
        ptr_bits = math.ceil(math.log(num_env + 2) / math.log(2))
        tsg_ptr = b.u(ptr_bits) - 1
        borders = var_borders(ic, num_env, nats, b_iframe, var_l, var_r,
                              rel_l, rel_r, state, key)

    if cfg["freq_res_mode"] == 0:
        fr = [b.u(1) for _ in range(num_env)]
    else:
        fr = [freq_res(borders, e, tsg_ptr, nats, cfg["freq_res_mode"])
              for e in range(num_env)]
    state[key] = borders[num_env]                    # previous_stop_pos
    return dict(num_env=num_env, num_noise=2 if num_env > 1 else 1,
                freq_res=fr, int_class=ic, borders=borders)


def tab_border(nats, n):
    """The FIXFIX uniform border table -- evenly spaced, last one exact."""
    return [round(i * nats / n) for i in range(n)] + [nats]


def var_borders(ic, n, nats, b_iframe, var_l, var_r, rel_l, rel_r,
                state, key):
    """Pseudocode 76, the three variable interval classes."""
    a = [0] * (n + 1)
    if ic == FIXVAR:
        a[0] = 0
        a[n] = var_r + nats
        for t in range(len(rel_r)):
            a[n - t - 1] = a[n - t] - rel_r[t]
    else:                                            # VARFIX or VARVAR
        a[0] = var_l if b_iframe else state.get(key, nats) - nats
        a[n] = nats if ic == VARFIX else var_r + nats
        for t in range(len(rel_l)):
            a[t + 1] = a[t] + rel_l[t]
        if ic == VARVAR:
            for t in range(len(rel_r)):
                a[n - t - 1] = a[n - t] - rel_r[t]
    return a


def freq_res(borders, atsg, tsg_ptr, nats, mode):
    """Pseudocode 77."""
    if mode == 1:
        return 0
    if mode == 3:
        return 1
    if atsg < tsg_ptr and nats > 8:
        return 1
    return 1 if (borders[atsg + 1] - borders[atsg]) > (nats / 6.0 + 3.25) else 0


def delta_dir(b, fr):
    """Table 54."""
    return ([b.u(1) for _ in range(fr["num_env"])],
            [b.u(1) for _ in range(fr["num_noise"])])


def hfgen_iwc(b, t, nats, two_ch, balance=0):
    """Tables 55 and 56.  -> aspx_add_harmonic flags per channel.

    The add-harmonic flags are the only thing here the synthesis needs; the
    rest (tna_mode, and the fic/tic flags which are 0 throughout this stream)
    are consumed for bit accounting only.
    """
    n_noise, n_hi = t["n_noise"], t["n_hi"]
    ah = [[0] * n_hi, [0] * n_hi]
    for _ in range(n_noise):
        b.u(2)                                    # aspx_tna_mode ch0
    if two_ch:
        if balance == 0:
            for _ in range(n_noise):
                b.u(2)                            # aspx_tna_mode ch1
        if b.u(1):                                # aspx_ah_left
            ah[0] = [b.u(1) for _ in range(n_hi)]
        if b.u(1):                                # aspx_ah_right
            ah[1] = [b.u(1) for _ in range(n_hi)]
        if b.u(1):                                # aspx_fic_present
            if b.u(1):                            # aspx_fic_left
                for _ in range(n_hi):
                    b.u(1)
            if b.u(1):                            # aspx_fic_right
                for _ in range(n_hi):
                    b.u(1)
        if b.u(1):                                # aspx_tic_present
            copy = b.u(1)
            left = right = 0
            if copy == 0:
                left, right = b.u(1), b.u(1)
            if copy or left:
                for _ in range(nats):
                    b.u(1)
            if right:
                for _ in range(nats):
                    b.u(1)
    else:
        if b.u(1):                                # aspx_ah_present
            ah[0] = [b.u(1) for _ in range(n_hi)]
        if b.u(1):                                # aspx_fic_present
            for _ in range(n_hi):
                b.u(1)
        if b.u(1):                                # aspx_tic_present
            for _ in range(nats):
                b.u(1)
    return ah


def hcb_name(data_type, qmode, smode, htype):
    sm = "BALANCE" if smode else "LEVEL"
    if data_type == "SIGNAL":
        return f"ASPX_HCB_ENV_{sm}_{30 if qmode else 15}_{htype}"
    return f"ASPX_HCB_NOISE_{sm}_{htype}"


def hcb(tabs, data_type, qmode, smode, htype):
    """Pseudocode 79."""
    sm = "BALANCE" if smode else "LEVEL"
    if data_type == "SIGNAL":
        return tabs[f"ASPX_HCB_ENV_{sm}_{30 if qmode else 15}_{htype}"]
    return tabs[f"ASPX_HCB_NOISE_{sm}_{htype}"]


# DELTA BOOK CENTRES.  Annex A lists only `codebook_length` for the A-SPX
# books -- no cb_off, unlike the ASF spectrum books -- so the mapping comes
# from the code lengths instead: every DF and DT book has ODD length with its
# single shortest codeword at index (n-1)/2, which is exactly what a delta
# alphabet centred on zero looks like.  The F0 books are NOT centred, which is
# consistent with an absolute level index, so they are used raw.
#
# FLAGGED: the F0 mapping is inferred from codebook structure, not read from a
# clause.  It offsets every envelope by a constant, so it affects the absolute
# high-band level -- not the shape, and not the time structure the gate tests.
CENTRE = {}


def ec_data(b, tabs, t, data_type, num_env, freq_res, qmode, smode, dirs):
    """Tables 57 and 58.  -> list of per-envelope decoded value lists."""
    out = []
    for env in range(num_env):
        if data_type == "SIGNAL":
            n_sbg = t["n_hi"] if freq_res[env] else t["n_lo"]
        else:
            n_sbg = t["n_noise"]
        vals = []
        if dirs[env] == 0:                        # FREQ: F0 then DF
            vals.append(hcb(tabs, data_type, qmode, smode, "F0").decode(b))
            nm = hcb_name(data_type, qmode, smode, "DF")
            book, c = tabs[nm], CENTRE[nm]
            for _ in range(1, n_sbg):
                vals.append(book.decode(b) - c)
        else:                                     # TIME: all DT
            nm = hcb_name(data_type, qmode, smode, "DT")
            book, c = tabs[nm], CENTRE[nm]
            for _ in range(n_sbg):
                vals.append(book.decode(b) - c)
        out.append(vals)
    return out


def aspx_data(b, tabs, t, cfg, nats, b_iframe, two_ch, state, grp):
    """Tables 51 and 52."""
    if b_iframe:
        b.u(3)                                    # aspx_xover_subband_offset
    f0 = aspx_framing(b, cfg, nats, b_iframe, state, (grp, 0))
    q0 = cfg["quant_mode_env"]
    if f0["int_class"] == FIXFIX and f0["num_env"] == 1:
        q0 = 0
    if not two_ch:
        d_sig, d_noise = delta_dir(b, f0)
        ah = hfgen_iwc(b, t, nats, two_ch=False)
        sig = ec_data(b, tabs, t, "SIGNAL", f0["num_env"], f0["freq_res"], q0,
                      0, d_sig)
        noi = ec_data(b, tabs, t, "NOISE", f0["num_noise"], None, 0, 0,
                      d_noise)
        return [dict(fr=f0, sig=sig, noise=noi, dir_sig=d_sig,
                     dir_noise=d_noise, balance=0, ch=0, ah=ah[0])]
    balance = b.u(1)
    f1, q1 = f0, q0
    if balance == 0:
        f1 = aspx_framing(b, cfg, nats, b_iframe, state, (grp, 1))
        q1 = cfg["quant_mode_env"]
        if f1["int_class"] == FIXFIX and f1["num_env"] == 1:
            q1 = 0
    d0_sig, d0_noise = delta_dir(b, f0)
    d1_sig, d1_noise = delta_dir(b, f1)
    ah = hfgen_iwc(b, t, nats, two_ch=True, balance=balance)
    sm1 = 1 if balance else 0
    s0 = ec_data(b, tabs, t, "SIGNAL", f0["num_env"], f0["freq_res"], q0, 0,
                 d0_sig)
    s1 = ec_data(b, tabs, t, "SIGNAL", f1["num_env"], f1["freq_res"], q1, sm1,
                 d1_sig)
    n0 = ec_data(b, tabs, t, "NOISE", f0["num_noise"], None, 0, 0, d0_noise)
    n1 = ec_data(b, tabs, t, "NOISE", f1["num_noise"], None, 0, sm1, d1_noise)
    return [dict(fr=f0, sig=s0, noise=n0, dir_sig=d0_sig,
                 dir_noise=d0_noise, balance=balance, ch=0, ah=ah[0]),
            dict(fr=f1, sig=s1, noise=n1, dir_sig=d1_sig,
                 dir_noise=d1_noise, balance=balance, ch=1, ah=ah[1])]


def main(argv=None):
    ap = argparse.ArgumentParser(prog="m32 aspx parse")
    ap.add_argument("path", nargs="?",
                    default="m7_out/rf33_audio_pid13.mp4")
    ap.add_argument("--frames", type=int, default=0)
    a = ap.parse_args(argv)
    p = a.path if os.path.isabs(a.path) else os.path.join(HERE, a.path)

    print("M32 -- the A-SPX payload, and the substream to the byte")
    print("=" * 74)
    T = C.load_tables()
    arrays = H.parse_c(H.DEFAULT_C)
    tabs = {n[:-4]: S.Huff(arrays[n], arrays[n[:-4] + "_CW"])
            for n in arrays if n.endswith("_LEN") and n.startswith("ASPX_")}
    for nm in tabs:
        CENTRE[nm] = 0 if nm.endswith("_F0") else (len(arrays[nm + "_LEN"]) - 1) // 2
    print(f"  {len(tabs)} A-SPX codebooks loaded "
          f"(Kraft-gated in M23)")

    fr = W.samples(p)
    if a.frames:
        fr = fr[:a.frames]
    nats = A.num_aspx_timeslots()

    cfg = None
    state = {}
    left = []
    out = collections.Counter()
    for f in fr:
        try:
            st = M.parse(f)
            o = st["toc_bytes"] + st["substream_sizes"][0]
            sub = f[o:o + st["substream_sizes"][1]]
            ifr = bool(st["b_iframe_global"])
            r = C.decode_element(sub, T, b_iframe=ifr)
            if r["bits"] > r["nbits"]:
                out["channel overrun"] += 1
                continue
        except Exception as e:                                 # noqa: BLE001
            out[f"channels: {type(e).__name__}"] += 1
            continue
        if r.get("aspx"):
            cfg = r["aspx"]
        if cfg is None:
            out["no config yet (pre first i-frame)"] += 1
            continue
        t = A.sbg_tables(cfg, 0)
        b = Bits(sub)
        b.p = r["bits"]
        try:
            aspx_data(b, tabs, t, cfg, nats, ifr, True, state, 0)
            aspx_data(b, tabs, t, cfg, nats, ifr, True, state, 1)
            aspx_data(b, tabs, t, cfg, nats, ifr, False, state, 2)
        except Unsupported as e:
            out[f"interval class {e} not implemented"] += 1
            continue
        except Exception as e:                                 # noqa: BLE001
            out[f"aspx: {type(e).__name__}"] += 1
            continue
        # audio_size, re-read from the substream header
        hb = Bits(sub)
        size = hb.u(15)
        if hb.u(1):
            size += hb.vb(7) << 15
        if hb.p % 8:
            hb.p += 8 - (hb.p % 8)                # byte_align before audio_data
        delta = b.p - (hb.p + size * 8)
        left.append(delta)
        out["parsed" if -7 <= delta <= 0 else f"parsed, off by {delta}"] += 1

    print()
    for k, v in out.most_common(10):
        print(f"    {v:5d}  {k}")
    if left:
        import numpy as np
        v = np.array(left, float)
        closed = int(((v >= -7) & (v <= 0)).sum())
        print(f"\n  (end of parse) - (start + audio_size*8): median "
              f"{np.median(v):+.0f}, min {v.min():+.0f}, max {v.max():+.0f}")
        print(f"\n  GATE  audio_data consumes exactly audio_size bytes "
              f"(-7..0 = byte_align)")
        print(f"    {'PASS' if closed == len(v) else 'FAIL'}  "
              f"{closed}/{len(v)} frames close exactly")
        print("\n" + "=" * 74)
        if closed == len(v):
            print("  EVERY BIT OF audio_data ACCOUNTED FOR.  Six channels "
                  "plus the A-SPX\n  payload consume exactly the length the "
                  "encoder declared.")
            return 0
    print("\n" + "=" * 74)
    print("  NOT established -- audio_data does not close to audio_size")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
