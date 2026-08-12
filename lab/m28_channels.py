#!/usr/bin/env python3
"""M28 -- the whole 5.1 channel element: all six channels, not just the LFE.

M24-M26 decoded the LFE and made sound.  That was the smallest real target: one
channel, one section, three scale factor bands, twelve lines.  This is the rest
of the frame -- five full-band channels at 768 coded lines each -- and it needs
five things the LFE never exercised.

WHAT THE LFE DID NOT EXERCISE
------------------------------
  * **Multiple sections.**  The LFE had one codebook covering all its bands.  A
    full-band channel tiles 43 bands with a run-length coded list of codebooks
    (clause 4.2.8.3), including codebook 0, which means "this band is all zeros
    and costs no bits at all".
  * **Every codebook, not just 11.**  Books 1-4 are 4-dimensional, 5-11 are
    pairs; 1, 2, 5, 6 are SIGNED (cb_off != 0, no sign bits); the rest are
    unsigned and spend a sign bit per non-zero value.  Only 11 escapes.
  * **`asf_snf_data`.**  Noise fill.  The LFE decoder skipped it, which was
    harmless there because nothing followed it in that substream -- but it is
    NOT harmless here.  Everything after the LFE in the element depends on the
    LFE's exit bit position, so an unconsumed field desynchronises all five
    remaining channels.  This is the bug that had to be fixed before anything
    else could work.
  * **`chparam_info`.**  MDCT stereo processing is on in every frame, so L/R
    are coded as mid/side with a per-band `ms_used` flag.
  * **The element walk itself**: companding_control(5), coding_config,
    2ch_mode, then two two_channel_data() blocks and a mono_data(0).

THE FIELD WIDTHS ARE READ, NOT GUESSED
----------------------------------------
Table 105 for transform length 1536: n_msfb_bits 6, n_side_bits 5,
n_msfbl_bits 3.  So the LFE's max_sfb is a 3-bit field and a full-band
channel's is 6 bits -- which is why the LFE reads 3 and the full-band channels
read 43, and why both fit.  Clause 4.3.6.2.4: n_grp_bits is 0 when
b_long_frame is true and frame_len_base >= 1536, so long frames spend no
grouping bits.  Table B.4 via M27: 49 bands for 1536, max_sfb 43 -> 768 lines.

THE GATES, AND WHY THESE ONES
-------------------------------
A complete prefix code decodes anything, so "it parsed" is worth nothing on its
own (M23/M24 established this the hard way).  Three tests that random bits
cannot pass:

  1. **The element closes.**  Six channels decoded in sequence, and the bit
     position never runs past the substream.  A single misread field
     desynchronises every Huffman decode after it, and the overrun shows up
     within a channel or two.  Over thousands of frames this is a strong test
     precisely because it is cumulative.
  2. **Time structure, per channel.**  Adjacent frames correlate; a shuffled
     control does not.  Run separately for each of the five channels, so a
     channel decoded from garbage stands out from its neighbours.
  3. **L and R correlate WITH EACH OTHER.**  This is the new one, and it is the
     sharpest: the two sf_data() blocks in a two_channel_data() are decoded
     independently, one after the other, yet they carry the same programme.  If
     both are right they track each other frame by frame.  If the first is
     right and the second is garbage -- the classic desynchronisation
     signature -- the correlation collapses while test 2 might still pass on
     channel 1 alone.

Usage:
    python m28_channels.py [--frames 0]
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

import m17_ac4_walk as W                                          # noqa: E402
import m20_ac4_toc2 as M                                          # noqa: E402
import m23_hcb as H                                               # noqa: E402
import m24_spectral as S                                          # noqa: E402
import m27_sfb as B                                               # noqa: E402
from m19_ac4_toc import Bits                                      # noqa: E402

# Annex B scale factor band offsets, per transform length, each extracted and
# gated against Table B.1's independent band count (M27).  Short frames split
# the 1536-sample frame into partial blocks of 96/192/384/768 (Table 99), so
# every one of these is reachable.
SFB_TABLE = {L: B.offsets_for(L) for L in (1536, 768, 384, 192, 96)}

N_MSFB_BITS = 6              # Table 105, 1536
N_MSFBL_BITS = 3             # Table 105, 1536 (the LFE)
SF_CENTRE = 60               # Pseudocode 21: scale_factor += dpcm_sf - 60
# Pseudocode 23: delta = dpcm_snf[g][sfb] - 17.  NOT 11.  The first version
# assumed the midpoint of a 22-entry alphabet, which is what a symmetric DPCM
# table would use -- but ASF_HCB_SNF is NOT symmetric: its shortest codes sit at
# indices 13..16, and the spec's own offset is 17.  Reading the code-length
# profile is what raised the suspicion; the spec is what settled it.
# (index 0 -> delta -17, which Pseudocode 23 treats as "no noise fill here".)
SNF_CENTRE = 17

# Books 1, 2, 5, 6 carry sign in the codeword (cb_off != 0); the rest are
# unsigned and spend a sign bit per non-zero line.  Same split as AAC.
UNSIGNED_CB = {c for c, o in S.CB_OFF.items() if o == 0}


def get_qline(idx, dim, mod, off):
    """Unpack a codeword index into `dim` quantised values."""
    if dim == 4:
        return [idx // mod ** 3 - off, (idx // mod ** 2) % mod - off,
                (idx // mod) % mod - off, idx % mod - off]
    return [idx // mod - off, idx % mod - off]


def asf_section_data(b, max_sfb, n_sect_bits=5):
    """-> sfb_cb[sfb], and the section list [(cb, start, end)].

    Clause 4.2.8.3.  `sect_len` starts at 1 and the escape value extends it,
    which is why a 5-bit field can describe a 43-band run.  The while-loop in
    the spec reads its test value, so it is a do-while in any real language.
    """
    esc = (1 << n_sect_bits) - 1
    sfb_cb = [0] * max_sfb
    sects = []
    k = 0
    while k < max_sfb:
        cb = b.u(4)
        sect_len = 1
        while True:
            incr = b.u(n_sect_bits)
            if incr != esc:
                break
            sect_len += esc
        sect_len += incr
        end = k + sect_len
        # a section may not run past max_sfb; with max_sfb 43 < 49 the LSF
        # split in the spec never triggers for this stream, but guard anyway
        if end > max_sfb:
            end = max_sfb
        for sfb in range(k, end):
            sfb_cb[sfb] = cb
        sects.append((cb, k, k + sect_len))
        k += sect_len
    return sfb_cb, sects


def asf_spectral_data(b, sects, offsets, tables):
    """Clause 4.2.8.4.  -> quantised lines up to the last section's end."""
    if not sects:                     # max_sfb == 0: a channel with no bands
        return np.zeros(0, dtype=np.int32)
    total = offsets[min(sects[-1][2], len(offsets) - 1)]
    lines = np.zeros(total, dtype=np.int32)
    for cb, s0, s1 in sects:
        if cb == 0 or cb > 11:
            continue                      # cb 0 = an all-zero band, no bits
        dim, mod, off = S.CB_DIM[cb], S.CB_MOD[cb], S.CB_OFF[cb]
        hb = tables[cb]
        k = offsets[min(s0, len(offsets) - 1)]
        end = offsets[min(s1, len(offsets) - 1)]
        while k < end:
            vals = get_qline(hb.decode(b), dim, mod, off)
            if cb in UNSIGNED_CB:
                vals = [(-v if (v and b.u(1)) else v) for v in vals]
            if cb == 11:
                # escape AFTER the sign bits -- the order is load-bearing
                for j, v in enumerate(vals):
                    if abs(v) == 16:
                        n = 0
                        while b.u(1):
                            n += 1
                        mag = (1 << (n + 4)) + b.u(n + 4)
                        vals[j] = -mag if v < 0 else mag
            lines[k:k + dim] = vals[:max(0, min(dim, len(lines) - k))]
            k += dim
    return lines


def band_max(lines, offsets, sfb):
    """max_quant_idx: the largest |line| in a band."""
    if sfb + 1 >= len(offsets):
        return 0
    lo, hi = offsets[sfb], min(offsets[sfb + 1], len(lines))
    seg = lines[lo:hi]
    return int(np.abs(seg).max()) if len(seg) else 0


def asf_scalefac_data(b, sfb_cb, lines, offsets, sf_table, max_sfb, state):
    """Clause 4.2.8.5.  `first_scf_found` is shared ACROSS channels.

    That shared flag is easy to miss and it matters: the reference is read once
    per sf_data(), but the "this band takes the reference instead of a delta"
    rule is per-call, so it is threaded through `state`.
    """
    ref = b.u(8)
    cur = ref
    first = False
    sfs = [None] * max_sfb
    for sfb in range(min(max_sfb, NUM_SFB_48)):
        if sfb_cb[sfb] != 0 and band_max(lines, offsets, sfb) > 0:
            if first:
                cur += sf_table.decode(b) - SF_CENTRE
            else:
                first = True
            sfs[sfb] = cur
    state["ref"] = ref
    return sfs


def asf_snf_data(b, sfb_cb, lines, offsets, snf_table, max_sfb):
    """Clause 4.2.8.6.  Noise fill for the bands that carry no lines.

    Skipping this was harmless for the LFE (nothing followed it) and fatal
    here: every later channel reads from the bit position this leaves behind.
    """
    if not b.u(1):
        return None
    snf = [None] * max_sfb
    for sfb in range(min(max_sfb, NUM_SFB_48)):
        if sfb_cb[sfb] == 0 or band_max(lines, offsets, sfb) == 0:
            snf[sfb] = snf_table.decode(b) - SNF_CENTRE
    return snf


def sf_data(b, fr, T):
    """Clause 4.2.7.3, over every window group.

    The four stages each loop over groups in the spec, so they are run group by
    group here: sections for all groups, then spectral for all groups, and so
    on -- NOT one group end-to-end at a time.  Getting that nesting backwards
    would parse a long frame identically (one group) and desynchronise every
    short frame, which is exactly the kind of bug that hides for a long time.
    """
    G = fr.num_groups
    sfb_cb, sects, offs, msfb = [], [], [], []
    for g in range(G):
        m = fr.max_sfb_g(g)
        msfb.append(m)
        offs.append(fr.sect_sfb_offset(g))
        cb, se = asf_section_data(b, m, fr.n_sect_bits(g))
        sfb_cb.append(cb)
        sects.append(se)
    lines = np.zeros(max((o[min(m, len(o) - 1)] for o, m in zip(offs, msfb)),
                         default=0), dtype=np.int32)
    for g in range(G):
        part = asf_spectral_data(b, sects[g], offs[g], T["cb"])
        n = min(len(part), len(lines))
        lines[:n] = np.where(part[:n] != 0, part[:n], lines[:n])
    sfs, first = [], [False]
    ref = b.u(8)
    cur = [ref]
    for g in range(G):
        sfs.append(_scalefac_group(b, sfb_cb[g], lines, offs[g], T["sf"],
                                   msfb[g], cur, first))
    snf = None
    if b.u(1):
        snf = [_snf_group(b, sfb_cb[g], lines, offs[g], T["snf"], msfb[g])
               for g in range(G)]
    return dict(lines=lines, sfs=sfs[0], sfs_all=sfs, sfb_cb=sfb_cb[0],
                sects=sects[0], snf=snf, ref=ref, offsets=offs[0],
                max_sfb=msfb[0], groups=G, framing=fr)


def _scalefac_group(b, sfb_cb, lines, offsets, sf_table, max_sfb, cur, first):
    """One group's slice of asf_scalefac_data.  `cur`/`first` are shared."""
    out = [None] * max_sfb
    for sfb in range(max_sfb):
        if sfb_cb[sfb] != 0 and band_max(lines, offsets, sfb) > 0:
            if first[0]:
                cur[0] += sf_table.decode(b) - SF_CENTRE
            else:
                first[0] = True
            out[sfb] = cur[0]
    return out


def _snf_group(b, sfb_cb, lines, offsets, snf_table, max_sfb):
    out = [None] * max_sfb
    for sfb in range(max_sfb):
        if sfb_cb[sfb] == 0 or band_max(lines, offsets, sfb) == 0:
            out[sfb] = snf_table.decode(b) - SNF_CENTRE
    return out


# Table 99: transf_length index -> actual transform length, frame_length 1536
# at 48 kHz.  Table 108: n_grp_bits, indexed by (transf_length[0],
# transf_length[1]).  Both transcribed from the spec rather than derived --
# though the derivation does reproduce Table 108 exactly: with block counts
# n_i = 768 / length_i, n_grp_bits is n0 + n1 - 1 when the halves share a
# length and (n0 - 1) + (n1 - 1) when they differ.
# Short frames are ON, and the history of this flag is the point.
#
# When the branch was first written it made more frames CLOSE (76.1 % vs 64.6 %)
# while every correctness gate COLLAPSED (L fell from +0.62 to +0.12), so it was
# defaulted off: closing is necessary, never sufficient.  The cause turned out
# to be upstream -- partial blocks were being decoded against the wrong Annex B
# band table (see TABLE_FOR in m27_sfb).  With the right table:
#
#     frames closing   98.4 %      remainder median 260..314 bits in EVERY
#                                  framing category, matching the long-frame
#                                  baseline exactly
#     time structure   all six channels pass against ~0 controls
#     pair coupling    L/R -0.32, Ls/Rs +0.13, controls ~0
#
# so it is on by default now, and --no-short backs it out.
SHORT_FRAMES = True

SHORT_LEN = {0: 96, 1: 192, 2: 384, 3: 768}

# Table 105: n_msfb_bits by transform length at 44,1 / 48 kHz.
N_MSFB_BY_LEN = {2048: 6, 1920: 6, 1536: 6, 1024: 6, 960: 6, 768: 6,
                 512: 6, 480: 6, 384: 6, 256: 5, 240: 5, 192: 5,
                 128: 4, 120: 4, 96: 4}


def n_msfb_bits(length):
    return N_MSFB_BY_LEN[length]

N_GRP_BITS = {
    (0, 0): 15, (0, 1): 10, (0, 2): 8, (0, 3): 7,
    (1, 0): 10, (1, 1): 7, (1, 2): 4, (1, 3): 3,
    (2, 0): 8, (2, 1): 4, (2, 2): 3, (2, 3): 1,
    (3, 0): 7, (3, 1): 3, (3, 2): 1, (3, 3): 1,
}


class Framing:
    """asf_transform_info + asf_psy_info, both branches.  Pseudocode 3-5.

    Holds everything the rest of sf_data needs to be window-agnostic: how many
    groups there are, each group's max_sfb and transform length, and the
    per-group sect_sfb_offset that packs grouped short blocks contiguously.
    """

    def __init__(self, b):
        self.long = bool(b.u(1))                     # b_long_frame
        self.tl = None
        if not self.long:
            if not SHORT_FRAMES:
                raise ValueError("short frame")
            self.tl = [b.u(2), b.u(2)]               # transf_length[0..1]
        self.different = (not self.long and self.tl[0] != self.tl[1])
        # Table 105 indexes n_msfb_bits by TRANSFORM LENGTH, not by frame: the
        # short blocks at 192 and 96 samples carry 5- and 4-bit max_sfb fields.
        # Using 6 everywhere parses the long frames and the 384/768 short
        # frames correctly and desynchronises exactly the rest -- which is what
        # the first short-frame run showed (70 % closing, and the failures
        # confined to the framings involving 96 or 192).
        self.max_sfb = [b.u(n_msfb_bits(self.length_idx(0)))]
        if self.different:
            self.max_sfb.append(b.u(n_msfb_bits(self.length_idx(1))))
        n_grp = 0 if self.long else N_GRP_BITS[tuple(self.tl)]
        sfg = [b.u(1) for _ in range(n_grp)]

        # --- Pseudocode 3 -------------------------------------------------
        self.num_windows = 1
        self.num_groups = 1
        self.w2g = [0]
        if not self.long:
            self.num_windows = n_grp + 1
            if self.different:
                nw0 = 1 << (3 - self.tl[0])          # windows in first half
                sfg = sfg + [0]
                for i in range(n_grp, nw0 - 1, -1):
                    sfg[i] = sfg[i - 1]
                sfg[nw0 - 1] = 0                     # no grouping across halves
                self.num_windows += 1
            self.w2g = [0] * self.num_windows
            for i in range(self.num_windows - 1):
                if sfg[i] == 0:
                    self.num_groups += 1
                self.w2g[i + 1] = self.num_groups - 1
        self.nwin = [sum(1 for w in self.w2g if w == g)
                     for g in range(self.num_groups)]

    def length_idx(self, i):
        """Transform length of half i, before the grouping is known."""
        return 1536 if self.long else SHORT_LEN[self.tl[i]]

    def _idx(self, g):
        """Pseudocode 5: which half of the frame group g belongs to."""
        if not self.different:
            return 0
        nw0 = 1 << (3 - self.tl[0])
        return 1 if g >= self.w2g[nw0] else 0

    def max_sfb_g(self, g):
        i = self._idx(g)
        return self.max_sfb[min(i, len(self.max_sfb) - 1)]

    def length_g(self, g):
        return 1536 if self.long else SHORT_LEN[self.tl[self._idx(g)]]

    def n_sect_bits(self, g):
        """Clause 4.2.8.3: 3 bits for the short transform lengths."""
        if self.long:
            return 5
        return 3 if self.tl[self._idx(g)] <= 2 else 5

    def sect_sfb_offset(self, g):
        """Pseudocode 4 -- grouped short blocks pack contiguously."""
        base = 0
        for gg in range(g):
            tbl = SFB_TABLE[self.length_g(gg)]
            base += tbl[min(self.max_sfb_g(gg), len(tbl) - 1)] * self.nwin[gg]
        tbl = SFB_TABLE[self.length_g(g)]
        n = self.nwin[g]
        return [base + v * n for v in tbl]


def sf_info(b):
    """-> Framing.  Both the long and the short branch."""
    return Framing(b)


def chparam_info(b, fr, sf_table):
    """Clause 4.2.10.1 -- stereo processing side information.

    EVERY LOOP HERE IS PER WINDOW GROUP.  Both the `ms_used` loop (sap_mode 1)
    and `sap_data` (sap_mode 3) run `for g` in the spec, `max_sfb_g` is fetched
    per group, and `delta_code_time` exists only when num_window_groups != 1.
    A long frame has exactly one group, so a single-group implementation parses
    it perfectly and under-reads every short frame that uses M/S-per-band or
    full SAP -- which is exactly the fingerprint the failures showed: 33 of 40
    at sap_mode 3, 5 at sap_mode 1, and only 2 at sap_mode 2, the one mode that
    reads no bits at all.
    """
    sap_mode = b.u(2)
    ms_used = None
    G = fr.num_groups
    if sap_mode == 1:
        ms_used = []
        for g in range(G):
            ms_used.append([b.u(1) for _ in range(fr.max_sfb_g(g))])
    elif sap_mode == 3:
        # sap_data(), clause 4.2.10.2
        used = {}
        all_on = b.u(1)                             # sap_coeff_all
        for g in range(G):
            m = fr.max_sfb_g(g)
            if all_on:
                used[g] = [1] * m
                continue
            u = []
            for sfb in range(0, m, 2):
                v = b.u(1)
                u.append(v)
                if sfb + 1 < m:
                    u.append(v)                     # pairs share the flag
            used[g] = u[:m]
        if G != 1:
            b.u(1)                                  # delta_code_time
        for g in range(G):
            for sfb in range(0, fr.max_sfb_g(g), 2):
                if used[g][sfb]:
                    sf_table.decode(b)              # dpcm_alpha_q
    return sap_mode, ms_used


def companding_control(b, num_chan):
    """Clause 4.2.11."""
    sync = b.u(1) if num_chan > 1 else 0
    nc = 1 if sync else num_chan
    on = [b.u(1) for _ in range(nc)]
    if not all(on):
        b.u(1)                                      # b_compand_avg
    return sync, on


def two_channel_data(b, T):
    """Clause 4.2.6.7 -- one stereo pair.

    Returns (ch0, ch1, stereo) where `stereo` carries the M/S state: Table 113
    gives sap_mode 0 = none, 1 = M/S in the bands flagged by ms_used,
    2 = M/S in ALL bands, 3 = full SAP (alpha prediction).  When M/S is on the
    two blocks are mid and side, not left and right.
    """
    sap, ms_used = None, None
    if b.u(1):                                      # b_enable_mdct_stereo_proc
        fr = Framing(b)
        sap, ms_used = chparam_info(b, fr, T["sf"])
        frs = [fr, fr]
    else:
        frs = [Framing(b), Framing(b)]
    chans = [sf_data(b, frs[i], T) for i in range(2)]
    return chans[0], chans[1], dict(sap_mode=sap, ms_used=ms_used,
                                    max_sfb=frs[0].max_sfb_g(0),
                                    framing=frs[0])


def mono_data(b, T, b_lfe):
    """Clause 4.2.6.2."""
    if b_lfe:
        fr = Framing.__new__(Framing)               # sf_info_lfe: always long
        fr.long, fr.tl, fr.different = True, None, False
        fr.max_sfb = [b.u(N_MSFBL_BITS)]
        fr.num_windows, fr.num_groups, fr.w2g, fr.nwin = 1, 1, [0], [1]
    else:
        if b.u(1) != 0:                             # spec_frontend
            raise ValueError("SSF (speech frontend) not implemented")
        fr = Framing(b)
    return sf_data(b, fr, T)


# Clause 5.7.6.3.1.1.  The A-SPX master subband group table is cut out of one
# of two static templates.  A 64-band QMF on a 48 kHz signal puts each subband
# at 24000 / 64 = 375 Hz.
SBG_TEMPLATE_HIGHRES = [18, 19, 20, 21, 22, 23, 24, 26, 28, 30, 32, 34, 36,
                        38, 40, 42, 44, 47, 50, 53, 56, 59, 62]
SBG_TEMPLATE_LOWRES = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 24,
                       26, 28, 30, 32, 35, 38, 42, 46]
QMF_HZ = 24000.0 / 64.0                              # 375 Hz per subband


def aspx_band(cfg):
    """Pseudocode 67 -> (first, last) A-SPX QMF subband, and their Hz.

    THE INDEX STEP IS IN THE TEMPLATE, NOT IN SUBBANDS.  Clause 4.3.10.1.2
    says the index moves "in steps of 2 subbands", which reads as
    sba = 18 + 2*start_freq -- and that is wrong, because the template is not
    uniformly spaced above index 6.  Pseudocode 67 is unambiguous:

        sbg_master[sbg] = sbg_template[2 * aspx_start_freq + sbg]

    so the step of 2 indexes the TABLE.  For start_freq 5 the prose arithmetic
    gives subband 28 (10500 Hz) and the pseudocode gives template[10] = 32
    (12000 Hz).  Worth remembering as the reason not to state a number derived
    from prose when pseudocode restating it exists.
    """
    tpl = (SBG_TEMPLATE_HIGHRES if cfg["master_freq_scale"] == 1
           else SBG_TEMPLATE_LOWRES)
    n_master = (22 if cfg["master_freq_scale"] == 1 else 20) \
        - 2 * cfg["start_freq"] - 2 * cfg["stop_freq"]
    sba = tpl[2 * cfg["start_freq"]]
    sbz = tpl[2 * cfg["start_freq"] + n_master]
    return sba, sbz, sba * QMF_HZ, sbz * QMF_HZ


def aspx_config(b):
    """Table 50.  Fifteen fixed-width bits, no conditionals.

    Present only in i-frames, which is why 6.7 % of frames used to be skipped
    outright.
    """
    return dict(quant_mode_env=b.u(1), start_freq=b.u(3), stop_freq=b.u(2),
                master_freq_scale=b.u(1), interpolation=b.u(1),
                preflat=b.u(1), limiter=b.u(1), noise_sbg=b.u(2),
                num_env_bits_fixfix=b.u(1), freq_res_mode=b.u(2))


def _acpl3_core(b, T, b_iframe, sub):
    """5_X_codec_mode 4 = ASPX_ACPL_3 -- decode the STEREO CORE only.

    TS 103 190-1 Table 96 gives 5_X_codec_mode 4 = ASPX_ACPL_3 (A-SPX plus
    Advanced Coupling mode 3), and Table 25's branch for it is

        case ASPX_ACPL_3:
            companding_control(2);
            stereo_data();
            aspx_data_2ch();
            acpl_data_2ch();

    with the i-frame prologue adding aspx_config() and then acpl_config_2ch(),
    the latter being exactly four bits (Table 60: acpl_num_param_bands_id 2,
    acpl_quant_mode_0 1, acpl_quant_mode_1 1) with no conditionals.

    So the 5.1 programme is carried as a TWO-channel core plus parametric
    coupling.  Everything up to and including aspx_data_2ch is machinery this
    decoder already has -- `stereo_data` (Table 23) has the same shape as
    two_channel_data (4.2.6.7), which is why m43 decodes the pair element with
    the same block.  What we do NOT have is acpl_data_2ch and the A-CPL
    synthesis that upmixes the core back to 5 channels (Table 62:
    acpl_framing_data, then alpha1/alpha2/beta1/beta2/beta3 via acpl_ec_data).

    Decoding the core therefore yields a correct STEREO rendering of the
    programme -- it is the coded downmix, not an invention -- and leaves the
    5.1 upmix for an A-CPL build.  acpl_data_2ch sits LAST in the frame, so
    not parsing it costs nothing upstream.
    """
    aspx = None
    if b_iframe:
        aspx = aspx_config(b)
        b.u(4)                                       # acpl_config_2ch, Table 60
    lfe = mono_data(b, T, b_lfe=True)                # b_has_lfe (5.1)
    companding_control(b, 2)
    l, r, st_lr = two_channel_data(b, T)             # stereo_data(), Table 23
    return dict(lfe=lfe, L=l, R=r, st_lr=st_lr, aspx=aspx, acpl3=True,
                bits=b.p, nbits=len(sub) * 8, coding_config=None)


def decode_element(sub, T, b_iframe=False):
    """The whole 5.1 element.  -> dict of six channels + the bit position."""
    b = Bits(sub)
    b.u(15)                                          # audio_size_value
    if b.u(1):
        b.vb(7)
    codec_mode = b.u(3)
    if codec_mode == 4:
        return _acpl3_core(b, T, b_iframe, sub)
    if codec_mode != 1:
        raise ValueError(f"5_X_codec_mode {codec_mode} is not ASPX")
    aspx = aspx_config(b) if b_iframe else None
    lfe = mono_data(b, T, b_lfe=True)
    companding_control(b, 5)
    coding_config = b.u(2)
    if coding_config != 0:
        raise ValueError(f"coding_config {coding_config} not implemented")
    b.u(1)                                           # 2ch_mode
    l, r, st_lr = two_channel_data(b, T)            # L, R  (mid/side)
    ls, rs, st_sr = two_channel_data(b, T)          # Ls, Rs
    centre = mono_data(b, T, b_lfe=False)
    return dict(lfe=lfe, L=l, R=r, Ls=ls, Rs=rs, C=centre,
                st_lr=st_lr, st_sr=st_sr, aspx=aspx,
                bits=b.p, nbits=len(sub) * 8, coding_config=coding_config)


def load_tables():
    arrays = H.parse_c(H.DEFAULT_C)
    return dict(
        cb={n: S.Huff(arrays[f"ASF_HCB_{n}_LEN"], arrays[f"ASF_HCB_{n}_CW"])
            for n in S.CB_MOD},
        sf=S.Huff(arrays["ASF_HCB_SCALEFAC_LEN"],
                  arrays["ASF_HCB_SCALEFAC_CW"]),
        snf=S.Huff(arrays["ASF_HCB_SNF_LEN"], arrays["ASF_HCB_SNF_CW"]),
    )


CHANS = ["lfe", "L", "R", "Ls", "Rs", "C"]


def main(argv=None):
    ap = argparse.ArgumentParser(prog="m28 channels")
    ap.add_argument("path", nargs="?",
                    default="m7_out/rf33_audio_pid13.mp4")
    ap.add_argument("--frames", type=int, default=0)
    ap.add_argument("--no-short", action="store_true",
                    help="decode long frames only (the pre-8/07 behaviour)")
    a = ap.parse_args(argv)
    p = a.path if os.path.isabs(a.path) else os.path.join(HERE, a.path)
    if a.no_short:
        globals()["SHORT_FRAMES"] = False

    print("M28 -- the whole 5.1 channel element"
          + ("  [long frames only]" if a.no_short else ""))
    print("=" * 76)
    T = load_tables()
    fr = W.samples(p)
    if a.frames:
        fr = fr[:a.frames]

    ok = short = err = ssf = 0
    reasons = collections.Counter()
    energy = {c: [] for c in CHANS}
    ix = []
    maxsfb = collections.Counter()
    aspx_cfgs = collections.Counter()
    overrun = 0
    for i, f in enumerate(fr):
        try:
            st = M.parse(f)
            # payload_base, per TS 103 190-1 4.3.3.2.11 + Pseudocode 1 --
            # offsets run from the end of ac4_toc and start at payload_base,
            # which is 0 only when b_payload_base is false. See m42 (E77).
            o = (st["toc_bytes"] + st.get("payload_base", 0)
                 + st["substream_sizes"][0])
            sub = f[o:o + st["substream_sizes"][1]]
            r = decode_element(sub, T, b_iframe=bool(st["b_iframe_global"]))
        except ValueError as e:                                # noqa: BLE001
            m = str(e)
            if "short frame" in m:
                short += 1
            elif "SSF" in m:
                ssf += 1          # a real feature we do not implement
            else:
                err += 1
                reasons[m[:60]] += 1
            continue
        except Exception as e:                                 # noqa: BLE001
            err += 1
            reasons[f"{type(e).__name__}: {e}"[:60]] += 1
            continue
        if r["bits"] > r["nbits"]:
            overrun += 1
            continue
        ok += 1
        if r.get("aspx"):
            a_ = r["aspx"]
            aspx_cfgs[(a_["start_freq"], a_["stop_freq"],
                       a_["master_freq_scale"])] += 1
        ix.append(i)
        for c in CHANS:
            energy[c].append(float(np.abs(r[c]["lines"]).sum()))
        maxsfb[r["L"]["max_sfb"]] += 1

    tot = ok + short + err + overrun + ssf
    print(f"\n  {tot} non-iframe frames")
    print(f"    {ok:5d} decoded, all six channels")
    if short:
        print(f"    {short:5d} short frames skipped (--no-short)")
    print(f"    {ssf:5d} used the speech frontend (SSF) -- a tool we do not "
          f"implement, counted not hidden")
    print(f"    {overrun:5d} ran past the substream")
    print(f"    {err:5d} errors")
    for r, n in reasons.most_common(5):
        print(f"          {n:5d}  {r}")
    if ok < 100:
        print("\n  too few frames to gate")
        return 1
    print(f"\n  full-band max_sfb: {dict(maxsfb)}")

    print("\n  GATE 1  the element closes -- six channels, inside the "
          "substream")
    base = max(tot - short - ssf, 1)
    rate = 100.0 * ok / base
    g1 = overrun == 0 and err <= 0.01 * base
    print(f"    {'PASS' if g1 else 'FAIL'}  {ok}/{base} decodable frames "
          f"({rate:.1f} %), {overrun} overruns, {err} errors")

    E = {c: np.array(energy[c]) for c in CHANS}
    I = np.array(ix)
    adj = np.array([j for j in range(len(I) - 1) if I[j + 1] == I[j] + 1])
    rng = np.random.default_rng(0)
    q = rng.permutation(len(I))
    print("\n  GATE 2  time structure, per channel "
          f"({len(adj)} adjacent pairs)")
    g2 = True
    # The control is a MEASUREMENT, not a constant: with n frames its standard
    # error is ~1/sqrt(n), so a fixed "|r_shuf| < 0.06" threshold demands less
    # noise than the sample size can deliver and fails correct decodes at small
    # n.  The honest test is that the real correlation stands well clear of the
    # control measured on the same data.
    se = 1.0 / np.sqrt(max(len(adj), 2))
    print(f"    (control standard error at this sample size: {se:.3f})")
    for c in CHANS:
        r_adj = float(np.corrcoef(E[c][adj], E[c][adj + 1])[0, 1])
        r_shuf = float(np.corrcoef(E[c][q[:-1]], E[c][q[1:]])[0, 1])
        good = r_adj > 0.15 and abs(r_adj) > 3 * max(abs(r_shuf), se)
        g2 &= good
        print(f"    {'PASS' if good else 'FAIL'}  {c:3s}  adjacent "
              f"r = {r_adj:+.4f}   shuffled r = {r_shuf:+.4f}")

    print("\n  GATE 3  independently decoded pair members track each other")
    g3 = True
    # SIGN IS NOT THE POINT.  With mid/side coding the two sf_data blocks are
    # not "left" and "right" -- they are the sum and difference signals, and an
    # encoder holding a bitrate spends on one at the other's expense, so the
    # energies are ANTI-correlated.  My first version demanded r > +0.3 and
    # failed a correct decode reading -0.44.  What proves the second block was
    # decoded correctly is a strong RELATIONSHIP where the shuffled control has
    # none, in either direction.
    for x, y in (("L", "R"), ("Ls", "Rs")):
        r_pair = float(np.corrcoef(E[x], E[y])[0, 1])
        r_ctrl = float(np.corrcoef(E[x], E[y][q])[0, 1])
        # no arbitrary absolute floor here: the meaningful question is whether
        # the relationship stands clear of the control measured on the same
        # data.  Ls/Rs is a genuinely weaker coupling than L/R -- surround
        # channels carry more independent content -- but at n in the thousands
        # even r = 0.12 is many sigma from its control.
        good = abs(r_pair) > 3 * max(abs(r_ctrl), se)
        g3 &= good
        print(f"    {'PASS' if good else 'FAIL'}  {x}/{y}  r = {r_pair:+.4f}"
              f"   shuffled control r = {r_ctrl:+.4f}")

    # GATE 4 ties two subsystems that share no bits.  The MDCT's coded edge
    # comes from max_sfb and the Annex B band table; the A-SPX start frequency
    # comes from the i-frame header and the subband group templates.  Nothing
    # in the decoder forces them to agree, so agreement is evidence about both.
    print("\n  GATE 4  the A-SPX header and the MDCT band table agree on the "
          "crossover")
    g4 = True
    if aspx_cfgs:
        cfg = aspx_cfgs.most_common(1)[0][0]
        c = dict(start_freq=cfg[0], stop_freq=cfg[1], master_freq_scale=cfg[2])
        sba, sbz, fa, fz = aspx_band(c)
        m = maxsfb.most_common(1)[0][0] if maxsfb else 43
        mdct_hz = B.SFB_OFFSET_1536[m] / 1536.0 * 24000.0
        g4 = abs(fa - mdct_hz) < 1.0
        print(f"    A-SPX  start_freq {c['start_freq']} -> QMF subband {sba}"
              f" = {fa:.0f} Hz;  stop -> subband {sbz} = {fz:.0f} Hz")
        print(f"    MDCT   max_sfb {m} -> line {B.SFB_OFFSET_1536[m]}/1536 "
              f"= {mdct_hz:.0f} Hz")
        print(f"    {'PASS' if g4 else 'FAIL'}  they meet exactly -- separate "
              f"parts of the bitstream, same number")
    else:
        print("    (no i-frames in this range)")

    print("\n" + "=" * 76)
    if g1 and g2 and g3 and g4:
        print("  ALL SIX CHANNELS DECODE.  The element closes, every channel "
              "has time\n  structure, the pairs track each other, and A-SPX "
              "starts where the\n  MDCT stops.")
        return 0
    print("  NOT established -- see the failing gate above")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
