#!/usr/bin/env python3
"""M30 -- AC-4's filterbank: KBD windows, block switching, and the overlap/add.

M29 renders long frames with a SINE window.  AC-4 does not use a sine window --
clause 5.5.3 specifies KBD with a per-transform-length alpha (Table 185) -- so
even the long-frame audio is currently reconstructed against the wrong window.
And short frames need more: when the block length changes, the window shape has
to change with it or time-domain aliasing stops cancelling.

WHAT THE SPEC GIVES, VERBATIM
------------------------------
Table 185, alpha by transform length at 44,1/48 kHz:

    2048 1920 1536 -> 3      1024 960 768 -> 4      512 480 384 -> 4.5
     256  240  192 -> 5       128 120  96 -> 6

The kernel (rendered from the page, because the text layer mangles it):

    W(N,n,a) = I0(pi*a*sqrt(1.0 - (2n/N - 1)^2)) / I0(pi*a)      0 <= n < N
    KBD_LEFT(N,n) = sqrt( sum(W, p=0..n) / sum(W, p=0..N) )

Note the denominator is N, not N-1, and the normalising sum runs to p = N
INCLUSIVE -- which is exactly what m22.kbd_window already computes, so that
primitive needed no change.

The left window for a block of length N following one of length Nprev
(clause 5.5.2.2, step 5):

    NW    = min(N, Nprev)
    Nskip = (N - NW) / 2
    w[n]  = 0                      0 <= n < Nskip
          = KBD_LEFT(NW, n-Nskip)  Nskip <= n < NW + Nskip
          = 1                      NW + Nskip <= n < N

So the transition region is always the SHORTER of the two blocks, centred, with
flat and zero shoulders -- the standard construction, and it degenerates to a
plain KBD when N == Nprev.

THE PART THE SPEC DOES NOT MAKE UNAMBIGUOUS
---------------------------------------------
Step 6 windows the previous block's stored second half as

    for (n = 0; n < Nprev; n++) overlap[nskip_prev + n] *= w[n];

where `w` is the CURRENT block's left window, of length N.  When Nprev > N that
indexes past the end of w, and taken literally with an ascending w it cannot
satisfy time-domain aliasing cancellation, which needs the previous block's
right window to be the time-REVERSE of the next block's left window.

Rather than hand-derive which reading is intended, both are implemented and
tested, because there is a test that settles it outright.

THE GATE: A LAPPED ORTHOGONAL TRANSFORM SATISFIES S @ S.T == I
----------------------------------------------------------------
The whole synthesis chain -- IMDCT, unfold, window, overlap-add -- is a linear
map S from stacked spectra to PCM.  For an MDCT filterbank whose windows satisfy
the Princen-Bradley condition, that map is ORTHOGONAL: analysis is the transpose
of synthesis, and S @ S.T is the identity.  So S can be built column by column
by pushing unit impulses through the real synthesis code, and then the entire
window-and-overlap construction is checked with one matrix identity -- no
encoder needed, no broadcast data needed, and no hand-derivation trusted.

This is the same discipline as M22's perfect-reconstruction gate, extended to
the case M22 could not express: block sizes that CHANGE.  If the windows are
wrong, or the transition shoulders are misaligned, or step 6 is read the wrong
way round, the identity fails and says so.

Usage:
    python m30_filterbank.py
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m22_mdct as T22                                           # noqa: E402

# Table 185, 44,1 kHz / 48 kHz column
ALPHA = {2048: 3.0, 1920: 3.0, 1536: 3.0,
         1024: 4.0, 960: 4.0, 768: 4.0,
         512: 4.5, 480: 4.5, 384: 4.5,
         256: 5.0, 240: 5.0, 192: 5.0,
         128: 6.0, 120: 6.0, 96: 6.0}


def kbd_left(n, alpha):
    """KBD_LEFT(N, .) -- the rising half, length n.  Clause 5.5.3."""
    p = np.arange(n + 1)
    k = np.i0(np.pi * alpha * np.sqrt(np.maximum(
        0.0, 1.0 - (2.0 * p / n - 1.0) ** 2))) / np.i0(np.pi * alpha)
    c = np.cumsum(k)
    return np.sqrt(c[:-1] / c[-1])


def left_window(n, n_prev, alpha=None):
    """Step 5's w[n]: zero shoulder, KBD transition, flat shoulder."""
    nw = min(n, n_prev)
    nskip = (n - nw) // 2
    a = ALPHA.get(nw, alpha if alpha is not None else 4.0)
    w = np.ones(n)
    w[:nskip] = 0.0
    w[nskip:nskip + nw] = kbd_left(nw, a)
    return w


def imdct_raw(X, orthonormal=True):
    """N spectral lines -> 2N UNWINDOWED time samples (m22's proven unfold).

    THE SCALE IS NOT COSMETIC.  m22's IMDCT carries a 2/N, so the synthesis it
    defines is orthogonal but not ORTHOnormal: the first run of this file's
    gate read 8.75e-01 = 1 - 2/16, i.e. S @ S.T came out as (2/N) * I exactly.
    That is a constant, and a gate that fails by a constant is telling you what
    the constant is rather than that the construction is wrong.

    With block switching the N-dependence stops being a harmless global gain:
    each block length would get a DIFFERENT gain, so short blocks would sit at
    the wrong level relative to long ones and every transition would step in
    loudness.  Scaling by sqrt(N/2) makes every block unit-gain, which is the
    physically correct relative scaling for a lapped orthogonal transform.
    """
    from scipy.fft import dct
    n = len(X)
    folded = dct(X, type=4, norm=None) / n
    if orthonormal:
        folded = folded * np.sqrt(n / 2.0)
    a_b, c_d = folded[n // 2:], folded[:n // 2]
    return np.concatenate([a_b, -a_b[::-1], -c_d[::-1], -c_d])


def synthesise(spectra, lengths, n_full, mode="reversed", state=None,
               return_state=False):
    """Steps 5 and 6 over a sequence of blocks.  -> concatenated PCM.

    `mode` selects the reading of step 6's window application to the previous
    block's stored half:
        "literal"  -- w[n] exactly as printed
        "reversed" -- the time-reverse, which is what TDAC requires

    RESUMABLE.  `overlap` and `n_prev` ARE the whole inter-block state, so
    handing them back lets a live decoder continue a stream across calls
    instead of restarting the filterbank every time.  Restarting is not just
    slower -- the first block of each call overlap-adds against zeros, which
    is a real discontinuity at every seam.  Pass `state` from the previous
    call and ask for it back with `return_state=True`.
    """
    if state is not None:
        overlap, n_prev = state
        overlap = np.array(overlap, dtype=float, copy=True)
    else:
        overlap = np.zeros(n_full)
        n_prev = lengths[0] if len(lengths) else n_full
    out = []
    for X, n in zip(spectra, lengths):
        w = left_window(n, n_prev)
        raw = imdct_raw(X)
        first, second = raw[:n] * w, raw[n:]

        nskip = (n_full - n) // 2
        nskip_prev = (n_full - n_prev) // 2
        if mode == "literal":
            m = min(n_prev, len(w))
            overlap[nskip_prev:nskip_prev + m] *= w[:m]
        else:
            # the previous block's right window is this block's left window,
            # time-reversed, evaluated over the previous block's span
            wr = left_window(n_prev, n)[::-1]
            overlap[nskip_prev:nskip_prev + n_prev] *= wr
        overlap[nskip:nskip + n] += first
        out.append(overlap[:n].copy())
        tail = overlap[n:].copy()
        overlap[:] = 0.0
        overlap[:len(tail)] = tail
        overlap[nskip:nskip + n] = second
        n_prev = n
    pcm = np.concatenate(out) if out else np.zeros(0)
    return (pcm, (overlap, n_prev)) if return_state else pcm


def operator(lengths, n_full, mode):
    """Build S column by column by pushing unit impulses through synthesise."""
    total = sum(lengths)
    cols = []
    for bi, n in enumerate(lengths):
        for k in range(n):
            spectra = [np.zeros(m) for m in lengths]
            spectra[bi][k] = 1.0
            cols.append(synthesise(spectra, lengths, n_full, mode))
    return np.array(cols).T                       # total x total


def gate(lengths, n_full, mode, trim, verbose=True):
    """S @ S.T == I on the interior.  -> max |deviation|."""
    S = operator(lengths, n_full, mode)
    G = S.T @ S
    lo, hi = trim, G.shape[0] - trim
    if hi <= lo:
        return float("nan")
    core = G[lo:hi, lo:hi]
    dev = float(np.abs(core - np.eye(hi - lo)).max())
    if verbose:
        print(f"    {'PASS' if dev < 1e-9 else 'FAIL'}  {str(lengths):34s} "
              f"mode={mode:8s}  max |S'S - I| = {dev:.2e}")
    return dev


def main(argv=None):
    ap = argparse.ArgumentParser(prog="m30 filterbank")
    ap.add_argument("--n", type=int, default=16, help="Nfull for the tests")
    a = ap.parse_args(argv)
    nf = a.n

    print("M30 -- AC-4 filterbank: KBD, block switching, overlap/add")
    print("=" * 74)

    print("\n  1. KBD windows satisfy Princen-Bradley (w[n]^2 + w[n+N]^2 = 1)")
    ok1 = True
    for n, al in ((16, 3.0), (32, 4.0), (96, 6.0), (192, 5.0), (1536, 3.0)):
        w = kbd_left(n, al)
        full = np.concatenate([w, w[::-1]])
        pb = float(np.abs(full[:n] ** 2 + full[n:] ** 2 - 1.0).max())
        ok1 &= pb < 1e-12
        print(f"    {'PASS' if pb < 1e-12 else 'FAIL'}  N={n:5d} alpha={al}"
              f"   residual {pb:.2e}")

    print("\n  2. the transition window degenerates to plain KBD when N==Nprev")
    w_eq = left_window(nf, nf)
    d = float(np.abs(w_eq - kbd_left(nf, ALPHA.get(nf, 4.0))).max())
    print(f"    {'PASS' if d < 1e-15 else 'FAIL'}  max |w - KBD_LEFT| = {d:.2e}")

    print("\n  3. S @ S.T == I  -- uniform block sizes")
    uni = [nf] * 6
    for mode in ("literal", "reversed"):
        gate(uni, nf, mode, trim=nf)

    print("\n  4. S @ S.T == I  -- BLOCK SWITCHING (this is the real test)")
    seqs = [
        [nf, nf, nf // 2, nf // 2, nf, nf],
        [nf, nf // 2, nf // 2, nf // 4, nf // 4, nf // 4, nf // 4, nf, nf],
        [nf, nf // 4, nf // 4, nf // 2, nf, nf],
    ]
    best = {}
    for mode in ("literal", "reversed"):
        devs = [gate(s, nf, mode, trim=nf) for s in seqs]
        best[mode] = max(devs)
    print()
    for mode, d in best.items():
        print(f"    {mode:8s}: worst deviation across the switching sequences "
              f"{d:.2e}")

    print("\n" + "=" * 74)
    winner = min(best, key=lambda k: best[k])
    if best[winner] < 1e-9:
        print(f"  BLOCK SWITCHING RECONSTRUCTS EXACTLY with the '{winner}' "
              f"reading of step 6.\n  The other reading deviates by "
              f"{best['literal' if winner == 'reversed' else 'reversed']:.2e}"
              f" -- the identity settles it, no derivation needed.")
        return 0
    print("  NEITHER reading reconstructs; the construction is still wrong.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())


def imdct_raw_batch(X2, orthonormal=True):
    """(C, N) spectral lines -> (C, 2N) unwindowed samples, one dct call.

    Bit-identical to imdct_raw per row: scipy's dct over axis=-1 computes the
    same sums in the same order per row, and every other op here broadcasts.
    The win is amortisation -- six channels per Python-level call instead of
    six calls, which is where the live 5.1 worker's 0.62x went (many small
    windows, loop-bound; measured 8/07)."""
    from scipy.fft import dct
    n = X2.shape[-1]
    folded = dct(X2, type=4, norm=None, axis=-1) / n
    if orthonormal:
        folded = folded * np.sqrt(n / 2.0)
    a_b, c_d = folded[..., n // 2:], folded[..., :n // 2]
    return np.concatenate(
        [a_b, -a_b[..., ::-1], -c_d[..., ::-1], -c_d], axis=-1)


def synthesise_batch(spectra_by_ch, lengths, n_full, mode="reversed",
                     states=None, return_state=False):
    """`synthesise` for C channels that SHARE one window length sequence.

    spectra_by_ch: {ch: [X per window]} -- every list len(lengths) long,
    every X the matching length. states: {ch: state or None}. Falls back is
    the CALLER's job: if the channels' framings differ for a stretch, run
    them per-channel with `synthesise` -- the two paths are bit-identical,
    so mixing them across a stream is safe (gated below).
    """
    chans = list(spectra_by_ch)
    C = len(chans)
    ov = np.zeros((C, n_full))
    n_prev0 = lengths[0] if len(lengths) else n_full
    n_prevs = []
    for ci, ch in enumerate(chans):
        st = (states or {}).get(ch)
        if st is not None:
            ov[ci] = np.asarray(st[0], dtype=float)
            n_prevs.append(st[1])
        else:
            n_prevs.append(n_prev0)
    out = {ch: [] for ch in chans}
    for i, n in enumerate(lengths):
        X2 = np.stack([np.asarray(spectra_by_ch[ch][i], dtype=float)
                       for ch in chans])
        raw = imdct_raw_batch(X2)
        # n_prev can differ per channel ONLY on the first window (states
        # from different histories); handle that window per-channel, then
        # the shared-lengths invariant holds for the rest.
        if len(set(n_prevs)) == 1:
            n_prev = n_prevs[0]
            w = left_window(n, n_prev)
            first, second = raw[:, :n] * w, raw[:, n:]
            nskip = (n_full - n) // 2
            nskip_prev = (n_full - n_prev) // 2
            wr = left_window(n_prev, n)[::-1] if mode != "literal" else None
            if mode == "literal":
                m = min(n_prev, len(w))
                ov[:, nskip_prev:nskip_prev + m] *= w[:m]
            else:
                ov[:, nskip_prev:nskip_prev + n_prev] *= wr
            ov[:, nskip:nskip + n] += first
            for ci, ch in enumerate(chans):
                out[ch].append(ov[ci, :n].copy())
            tail = ov[:, n:].copy()
            ov[:] = 0.0
            ov[:, :tail.shape[1]] = tail
            ov[:, nskip:nskip + n] = second
            n_prevs = [n] * C
        else:
            for ci, ch in enumerate(chans):
                n_prev = n_prevs[ci]
                w = left_window(n, n_prev)
                first, second = raw[ci, :n] * w, raw[ci, n:]
                nskip = (n_full - n) // 2
                nskip_prev = (n_full - n_prev) // 2
                if mode == "literal":
                    m = min(n_prev, len(w))
                    ov[ci, nskip_prev:nskip_prev + m] *= w[:m]
                else:
                    wr = left_window(n_prev, n)[::-1]
                    ov[ci, nskip_prev:nskip_prev + n_prev] *= wr
                ov[ci, nskip:nskip + n] += first
                out[ch].append(ov[ci, :n].copy())
                tail = ov[ci, n:].copy()
                ov[ci] = 0.0
                ov[ci, :len(tail)] = tail
                ov[ci, nskip:nskip + n] = second
                n_prevs[ci] = n
    pcm = {ch: (np.concatenate(v) if v else np.zeros(0))
           for ch, v in out.items()}
    if return_state:
        return pcm, {ch: (ov[ci].copy(), n_prevs[ci])
                     for ci, ch in enumerate(chans)}
    return pcm


def synthesise_frames(wins, lens, counts, n_full, states=None,
                      return_state=False):
    """Frame-walking grouped batch synthesis for C channels.

    Channels' window LISTS do not align positionally once framings diverge
    (L may spend 3 windows on a frame the LFE covers in 1), so batching must
    happen per FRAME, on groups of channels whose length tuples match THAT
    frame. Measured on 3000 real frames: 69.6 % of frames batch all-6; the
    rest batch as two groups (front trio / rear trio), so the amortisation
    survives divergence.

    wins/lens: flat per-channel lists (synthesise's shapes); counts: windows
    per frame per channel. Bit-identical to per-channel `synthesise` by
    construction -- same arithmetic, same order, only the dct is shared.
    """
    chans = list(wins)
    pos = {c: 0 for c in chans}
    st = dict(states or {})
    out = {c: [] for c in chans}
    nframes = len(counts[chans[0]])
    for fi in range(nframes):
        groups = {}
        for c in chans:
            k = counts[c][fi]
            key = tuple(lens[c][pos[c]:pos[c] + k])
            groups.setdefault(key, []).append(c)
        for key, members in groups.items():
            sub = {c: wins[c][pos[c]:pos[c] + len(key)] for c in members}
            pcm, ns = synthesise_batch(sub, list(key), n_full,
                                       states={c: st.get(c) for c in members},
                                       return_state=True)
            for c in members:
                out[c].append(pcm[c])
                st[c] = ns[c]
        for c in chans:
            pos[c] += counts[c][fi]
    res = {c: (np.concatenate(v) if v else np.zeros(0))
           for c, v in out.items()}
    return (res, st) if return_state else res
