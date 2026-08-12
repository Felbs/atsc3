"""M42 -- a resumable AC-4 decoder front half, for live use.

`m37.decode_all` reads a whole file and synthesises once. A live decoder has
to emit as it goes, and the only way to do that without re-deriving the past
is to hold the state that actually crosses a frame boundary:

    T, tabs   the Annex B tables and Huffman books -- constant, but building
              them per pass was pure waste
    cfg       the A-SPX configuration; it arrives on I-FRAMES ONLY, so a
              decoder that forgets it goes silent above the crossover until
              the next one
    fstate    A-SPX envelopes are DELTA CODED against the previous frame

This class holds all three and turns frames into spectra. The MDCT overlap
is the caller's to carry (m30.synthesise takes and returns it), because the
caller is the one who knows where a block boundary is.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import m20_ac4_toc2 as M                                         # noqa: E402
import m28_channels as C                                         # noqa: E402
import m29_audio as R                                            # noqa: E402
import m32_aspx_parse as P                                       # noqa: E402
import m23_hcb as H                                              # noqa: E402
import m24_spectral as S                                         # noqa: E402
import m31_aspx_bands as A                                       # noqa: E402
import m37_render_hf as R37                                      # noqa: E402


class Ac4Stream:
    def __init__(self, element="5_X"):
        # "5_X" = the 5.1 programme element (English main). "pair" = the
        # stereo channel_pair_element (Table 22) -- the Spanish simulcast.
        self.element = element
        self.T = C.load_tables()
        arrays = H.parse_c(H.DEFAULT_C)
        self.tabs = {n[:-4]: S.Huff(arrays[n], arrays[n[:-4] + "_CW"])
                     for n in arrays
                     if n.endswith("_LEN") and n.startswith("ASPX_")}
        for nm in self.tabs:
            P.CENTRE[nm] = 0 if nm.endswith("_F0") else \
                (len(arrays[nm + "_LEN"]) - 1) // 2
        self.nats = A.num_aspx_timeslots()
        self.cfg = None
        self.fstate = {}
        self.n_ok = 0
        self.n_bad = 0

    # 5.1 channel order everywhere downstream: matches m29's proven batch
    # render AND the standard WAV 6-channel layout (FL FR FC LFE BL BR), so a
    # plain 6-channel wav is a correct 5.1 file with no channel map needed.
    CHANS = ("L", "R", "C", "lfe", "Ls", "Rs")

    def decode_frames(self, frames, chans=None):
        """-> (wins per channel, lengths per channel, groups {"lr","sr"}).

        ALL SIX CHANNELS. `decode_element` always decoded the whole 5.1
        element -- this method used to keep only the first stereo pair and
        throw away C, LFE, Ls, Rs after their bits were already spent
        (found 8/07 via the spec index: TS 103 190-1 Table 25).

        EACH CHANNEL CARRIES ITS OWN FRAMING (m29's lesson, kept): the two
        pairs, the centre and the LFE are signalled independently and can
        differ in window count and lengths, so nothing here assumes one
        window sequence. Every frame still contributes exactly one frame of
        samples per channel, which is what keeps the six lanes equal-length
        and the timeline frame-accurate through bad frames.

        A-SPX: Table 25 puts THREE blocks in sequence -- aspx_data_2ch (L/R),
        aspx_data_2ch (Ls/Rs), aspx_data_1ch (C) -- and the parser's state is
        keyed by group index, so the pairs get groups 0 and 1. The centre's
        1-ch block is parsed to keep its delta-state warm, but not yet
        rendered: only apply_hf_pair exists, and a rushed mono HF port would
        put its mistakes on DIALOG. Core-band centre until that is built and
        gated.
        """
        # `chans` limits which channels are PACKED and SYNTHESISED -- the
        # element decode itself always walks all six (the bits are one
        # stream). ("L","R") restores the old stereo cost: measured 1.27x
        # real time against 0.62x for six channels even threaded, because
        # the filterbank is many small windows and the GIL never lets the
        # threads overlap. Live keeps up in stereo; recordings render 5.1.
        chans = tuple(chans) if chans else self.CHANS
        if self.element == "pair":
            chans = ("L", "R")           # a pair element HAS only a pair
        wins = {k: [] for k in chans}
        lengths = {k: [] for k in chans}
        counts = {k: [] for k in chans}   # windows per FRAME, per channel --
        # the grouped filterbank needs frame boundaries because channels'
        # window lists do not align positionally when framings diverge
        groups = {"lr": [], "sr": []}
        for f in frames:
            try:
                st = M.parse(f)
                # PAYLOAD BASE (TS 103 190-1 4.3.3.2.10/.11 + Pseudocode 1):
                # substream offsets are measured from the end of the
                # byte-aligned ac4_toc and START at payload_base --
                #     substream_n_offset = payload_base
                #     for (s = 0; s < n; s++) offset += substream_size[s]
                # -- and payload_base is 0 ONLY when b_payload_base is false
                # ("the first ac4_substream_data shall immediately follow the
                # ac4_toc element").  Omitting it was invisible for two years
                # because every RF33 lane signals b_payload_base = false; Fox
                # RF25 signals payload_base = 1 on every frame, so the element
                # was read one byte early and decoded as noise (E77).
                o = (st["toc_bytes"] + st.get("payload_base", 0)
                     + st["substream_sizes"][0])
                sub = f[o:o + st["substream_sizes"][1]]
                ifr = bool(st["b_iframe_global"])
                if self.element == "pair":
                    import m43_ac4_pair as PAIR
                    r = PAIR.decode_pair(sub, self.T, b_iframe=ifr)
                else:
                    r = C.decode_element(sub, self.T, b_iframe=ifr)
                if r["bits"] > r["nbits"]:
                    raise ValueError("element overran its substream")
                if r.get("aspx"):
                    self.cfg = r["aspx"]
                # noise fill seeded by the stream's own sequence counter,
                # exactly as the batch render does
                rng = np.random.default_rng(st["sequence_counter"])
                fm = {"L": r["st_lr"]["framing"], "R": r["st_lr"]["framing"]}
                # ASPX_ACPL_3 (E77) carries the 5.1 programme as a TWO-channel
                # core plus A-CPL; Ls/Rs/C do not exist as coded channels, so
                # this frame has only a pair to give no matter what `chans`
                # asked for. Narrow the request instead of KeyError-ing.
                if r.get("acpl3"):
                    chans = ("L", "R")
                elif self.element != "pair":
                    fm.update({"Ls": r["st_sr"]["framing"],
                               "Rs": r["st_sr"]["framing"],
                               "C": r["C"]["framing"],
                               "lfe": r["lfe"]["framing"]})
                pk = {k: R.packed_spectrum(r[k], fm[k], rng)
                      for k in chans}
                if "L" in pk:
                    pk["L"], pk["R"] = R.unmix(pk["L"], pk["R"],
                                               r["st_lr"], fm["L"])
                if "Ls" in pk:
                    pk["Ls"], pk["Rs"] = R.unmix(pk["Ls"], pk["Rs"],
                                                 r["st_sr"], fm["Ls"])
                u, seq = {}, {}
                for k in chans:
                    u[k], seq[k] = R.ungroup(pk[k], fm[k])
                got_lr = got_sr = None
                if self.cfg:
                    from m20_ac4_toc2 import Bits
                    b = Bits(sub)
                    b.p = r["bits"]
                    t = R37.t_for(self.cfg)
                    try:
                        got_lr = P.aspx_data(b, self.tabs, t, self.cfg,
                                             self.nats, ifr, True,
                                             self.fstate, 0)
                        if "Ls" in pk:
                            got_sr = P.aspx_data(b, self.tabs, t, self.cfg,
                                                 self.nats, ifr, True,
                                                 self.fstate, 1)
                            P.aspx_data(b, self.tabs, t, self.cfg,
                                        self.nats, ifr, False,
                                        self.fstate, 2)
                    except Exception:                          # noqa: BLE001
                        # a later block mis-parsing must not cost an earlier
                        # one: whatever was assigned before the throw stands
                        pass
                self.n_ok += 1
            except Exception:                                  # noqa: BLE001
                # A bad frame must not desynchronise the timeline: emit
                # silence of the right length so every downstream index
                # stays frame-accurate.
                u = {k: [np.zeros(1536)] for k in chans}
                seq = {k: [1536] for k in chans}
                got_lr = got_sr = None
                self.n_bad += 1
            for k in chans:
                wins[k].extend(u[k])
                lengths[k].extend(seq[k])
                counts[k].append(len(seq[k]))
            groups["lr"].append(got_lr)
            groups["sr"].append(got_sr)
        return wins, lengths, groups, counts
