#!/usr/bin/env python3
"""M43 -- channel_pair_element: the STEREO AC-4 element (Table 22).

Built 8/08 for the Spanish track. RF33's pid14 carries language 'spa' in its
own moov, and its audio substream is ~365 bytes against English's ~796 --
a stereo simulcast. `m28.decode_element` parses only the 5_X (5.1) element
and raised on every pid14 frame: 118,361 frames, 118,361 bad, decoded at
10x real time into perfect silence. A decoder that fails fast on the wrong
element type LOOKS like a healthy fast decoder in the throughput column --
the bad-frame counter was the only honest number on that line.

Syntax (TS 103 190-1 Table 22, channel_pair_element):
    stereo_codec_mode : 2 bits   (0 SIMPLE, 1 ASPX, 2/3 ASPX_ACPL)
    if (b_iframe && mode uses A-SPX)  aspx_config()
    SIMPLE: stereo_data()
    ASPX:   companding_control(2); stereo_data(); aspx_data_2ch()

`stereo_data` (Table 23) has the same shape as the 5_X pair block m28
already decodes (4.2.6.7 two_channel_data): one b_enable_mdct_stereo_proc
bit, then shared-or-separate framing + two sf_data. Reused, not rewritten.
The trailing aspx_data_2ch is the caller's to parse (m42 does it at
r["bits"], exactly as for 5_X).
"""
from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from m20_ac4_toc2 import Bits                                     # noqa: E402
from m28_channels import (aspx_config, companding_control,        # noqa: E402
                          two_channel_data)


def decode_pair(sub, T, b_iframe=False):
    """One channel_pair_element substream -> dict like decode_element's,
    with L/R only. Raises on ACPL modes (not seen on this air; counted by
    the caller as bad frames rather than mis-decoded)."""
    b = Bits(sub)
    b.u(15)                                          # audio_size_value
    if b.u(1):
        b.u(7)                                       # audio_size extension
    mode = b.u(2)                                    # stereo_codec_mode
    if mode >= 2:
        raise ValueError(f"stereo_codec_mode {mode} (ACPL) not implemented")
    aspx = None
    if b_iframe and mode == 1:
        aspx = aspx_config(b)
    if mode == 1:
        companding_control(b, 2)
    l, r, st = two_channel_data(b, T)
    return dict(L=l, R=r, st_lr=st, aspx=aspx,
                bits=b.p, nbits=len(sub) * 8)
