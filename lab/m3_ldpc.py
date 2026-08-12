#!/usr/bin/env python3
"""A/322 6.1.3.1 Type A / 6.1.3.2 Type B LDPC -- encoder, parity-check matrix,
and a normalized min-sum BP decoder.

WHY BOTH AN ENCODER AND AN H MATRIX
-----------------------------------
M2 already had a working min-sum decoder, gated on a SYNTHESIZED code.  That
proved decoder machinery only: our encoder against our decoder, which says
nothing about A/322.  Here the encoder (procedural, straight from the prose of
6.1.3.1) and the parity-check matrix (analytical, built from the same Annex A
table by a different route) are two INDEPENDENT code paths.  Requiring
H . encode(s) == 0 for random s catches construction bugs in either.

It still does not prove the TABLE is right -- only the air can do that, and the
BCH syndrome plus the 32-bit CRC in the decoded L1-Basic are the referees.

THE TWO-STAGE TYPE A STRUCTURE (why the table has Kldpc/360 + Q1 rows)
----------------------------------------------------------------------
Steps (i)-(vi) accumulate the Kldpc INFORMATION bits and produce the M1
dual-diagonal parity bits.  Steps (vii)-(viii) then feed those M1 parity bits
back in as if they were more information, using Q1 EXTRA table rows, to produce
the M2 identity-matrix parity bits.  This is why the printed rate-3/15 table has
12 rows and not 9, and it is confirmed structurally: in every Type A code the
tail rows contain only addresses >= M1, which is exactly what is required if
they may not disturb the already-finished dual-diagonal part.

ASSUMPTION L1 -- the parity interleave index order of steps (vi) and (viii).
    The PDF renders these two equations with the subscripts stripped by its
    math font, so the surviving text is unreadable as printed ("+360 + = 1+").
    Implemented as the only assignment that is a valid bijection:
        lambda_{K + 360*t + s}      = p_{Q1*s + t}          0<=s<360, 0<=t<Q1
        lambda_{K + M1 + 360*t + s} = p_{M1 + Q2*s + t}     0<=s<360, 0<=t<Q2
    The literal reading (360*s + t on the left) overflows M1 and is impossible,
    so this is forced rather than chosen -- but it is still an ASSUMPTION and
    `PARITY_INTERLEAVE_VARIANTS` exposes the alternative so the air can arbitrate.
    RISK LOW-MEDIUM: a wrong choice permutes parity and the decoder simply will
    not converge, which is a loud failure rather than a silent one.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import spec_ldpc as SL                                         # noqa: E402

try:                                                            # M10: Ninner = 64800
    import spec_ldpc64k as SL64                                 # noqa: E402
except ImportError:                                             # pragma: no cover
    SL64 = None

PARITY_INTERLEAVE_VARIANTS = ("standard", "swapped")


# --------------------------------------------------------------------------
# accumulator addressing (A/322 6.1.3.1 step ii)
# --------------------------------------------------------------------------

def _addr(x, mm, M1, M2, Q1, Q2):
    if x < M1:
        return (x + mm * Q1) % M1
    return M1 + ((x - M1 + mm * Q2) % M2)


def _pi(a, M1, M2, Q1, Q2, K, variant="standard"):
    """p-index -> lambda-index (ASSUMPTION L1)."""
    if a < M1:
        if variant == "standard":
            s, t = divmod(a, Q1)
            return K + 360 * t + s
        t, s = divmod(a, 360)
        return K + Q1 * s + t
    b = a - M1
    if variant == "standard":
        s, t = divmod(b, Q2)
        return K + M1 + 360 * t + s
    t, s = divmod(b, 360)
    return K + M1 + Q2 * s + t


def _params(rate, n=16200):
    """Parameters + Annex A address rows for (rate, Ninner).

    BUG FIXED HERE (M10): this used to read ``SL.LDPC_16200[rate]["rows"]``
    unconditionally, so a caller asking for Ninner = 64800 got the 16200-bit
    rows back with no complaint at all.  Every LLS-bearing PLP on RF25 and RF30
    signals fec_type = 1 (64800), so that silent fallback was the M8 blocker.
    Now a 64800 request that cannot be served raises instead of lying.
    """
    if n == 16200:
        P = SL.LDPC_PARAMS[rate]
        T = SL.LDPC_16200[rate]["rows"]
    elif n == 64800:
        if SL64 is None:
            raise RuntimeError(
                "Ninner=64800 requested but spec_ldpc64k.py is not present; "
                "run lab/extract_ldpc64k.py (needs PyMuPDF) to generate it.")
        P = SL64.LDPC_64800[rate]
        T = P["rows"]
    else:
        raise ValueError("Ninner must be 16200 or 64800, got %r" % (n,))
    K = P["Kldpc"]
    if P["type"] == "A":
        M1, M2, Q1, Q2 = P["M1"], P["M2"], P["Q1"], P["Q2"]
    else:
        M1, M2 = 0, n - K
        Q1, Q2 = 0, P["Qldpc"] or (n - K) // 360
    return P, T, K, M1, M2, Q1, Q2


# --------------------------------------------------------------------------
# encoder
# --------------------------------------------------------------------------

def encode(s, rate, n=16200, variant="standard"):
    """A/322 6.1.3.1 (Type A) / 6.1.3.2 (Type B).  s is a length-Kldpc {0,1} array."""
    P, T, K, M1, M2, Q1, Q2 = _params(rate, n)
    s = np.asarray(s, np.uint8)
    assert len(s) == K, f"expected {K} info bits, got {len(s)}"
    lam = np.zeros(n, np.uint8)
    lam[:K] = s
    p = np.zeros(M1 + M2, np.uint8)

    ninfo = K // 360
    for m in range(K):
        if not s[m]:
            continue
        row = T[m // 360]
        mm = m % 360
        for x in row:
            p[_addr(x, mm, M1, M2, Q1, Q2)] ^= 1

    if P["type"] == "A":
        # step (v): running XOR over the dual-diagonal part
        np.bitwise_xor.accumulate(p[:M1], out=p[:M1])
        # step (vi): interleave -> lambda_K .. lambda_{K+M1-1}
        for a in range(M1):
            lam[_pi(a, M1, M2, Q1, Q2, K, variant)] = p[a]
        # steps (vii): feed those M1 bits back through the extra Q1 rows
        for i in range(M1):
            b = lam[K + i]
            if not b:
                continue
            row = T[ninfo + i // 360]
            mm = i % 360
            for x in row:
                p[_addr(x, mm, M1, M2, Q1, Q2)] ^= 1
        # step (viii)
        for a in range(M1, M1 + M2):
            lam[_pi(a, M1, M2, Q1, Q2, K, variant)] = p[a]
    else:
        np.bitwise_xor.accumulate(p, out=p)
        for a in range(M1 + M2):
            lam[_pi(a, M1, M2, Q1, Q2, K, variant)] = p[a]
    return lam


# --------------------------------------------------------------------------
# analytical parity-check matrix (a DIFFERENT route from the same table)
# --------------------------------------------------------------------------

def parity_check(rate, n=16200, variant="standard"):
    """Return `rows`: list of lists of column indices, one per check equation."""
    P, T, K, M1, M2, Q1, Q2 = _params(rate, n)
    ninfo = K // 360
    checks = [[] for _ in range(M1 + M2)]

    for m in range(K):
        row = T[m // 360]
        mm = m % 360
        for x in row:
            checks[_addr(x, mm, M1, M2, Q1, Q2)].append(m)

    if P["type"] == "A":
        for i in range(M1):
            # step (vii) feeds back the CODEWORD bits lambda_{K+i} in codeword
            # order -- not the p-indexed accumulators.  Using _pi(i) here is
            # the natural-looking mistake and it makes H . encode(s) != 0.
            row = T[ninfo + i // 360]
            mm = i % 360
            for x in row:
                checks[_addr(x, mm, M1, M2, Q1, Q2)].append(K + i)
        for a in range(M1):
            checks[a].append(_pi(a, M1, M2, Q1, Q2, K, variant))
            if a >= 1:
                checks[a].append(_pi(a - 1, M1, M2, Q1, Q2, K, variant))
        for a in range(M1, M1 + M2):
            checks[a].append(_pi(a, M1, M2, Q1, Q2, K, variant))
    else:
        for a in range(M1 + M2):
            checks[a].append(_pi(a, M1, M2, Q1, Q2, K, variant))
            if a >= 1:
                checks[a].append(_pi(a - 1, M1, M2, Q1, Q2, K, variant))

    # a column appearing twice in a check cancels over GF(2)
    out = []
    for c in checks:
        c.sort()
        d = []
        i = 0
        while i < len(c):
            if i + 1 < len(c) and c[i] == c[i + 1]:
                i += 2
            else:
                d.append(c[i])
                i += 1
        out.append(d)
    return out


def flatten(checks, n):
    """CSR-ish flat arrays for the BP decoder, plus the variable-node index."""
    cn = np.array([len(c) for c in checks], np.int32)
    cptr = np.concatenate([[0], np.cumsum(cn)]).astype(np.int32)
    cidx = np.fromiter((v for c in checks for v in c), np.int32,
                       int(cn.sum()))
    return cptr, cidx


# --------------------------------------------------------------------------
# normalized min-sum BP
# --------------------------------------------------------------------------

def pack(checks, n):
    """Rectangular padded view of the check node structure (degrees vary 4..11).

    Returns (idx, mask) of shape (n_checks, max_degree).  Padded slots point at
    a dummy column n and are masked out, which lets the whole BP iteration run
    as array operations instead of a 12960-iteration Python loop.
    """
    nc = len(checks)
    dmax = max(len(c) for c in checks)
    idx = np.full((nc, dmax), n, np.int32)
    mask = np.zeros((nc, dmax), bool)
    for i, c in enumerate(checks):
        idx[i, :len(c)] = c
        mask[i, :len(c)] = True
    return idx, mask


def min_sum_decode(llr, checks, n, iters=50, alpha=0.75):
    """Normalized min-sum belief propagation.  llr > 0 means bit 0.

    Returns (bits, converged, n_iter, n_unsatisfied).
    """
    idx, mask = pack(checks, n)
    llr = np.asarray(llr, np.float64)
    tot = np.empty(n + 1)                        # last slot = dummy
    tot[:n] = llr
    tot[n] = 0.0
    m = np.zeros(idx.shape)
    BIG = 1e9
    bits = (llr < 0).astype(np.uint8)
    bad = -1
    for it in range(iters):
        v = tot[idx] - m
        v = np.where(mask, v, BIG)               # padding must not win the min
        sg = v < 0
        parity = np.bitwise_xor.reduce(sg & mask, axis=1)
        mag = np.abs(v)
        part = np.partition(mag, 1, axis=1)
        m1 = part[:, 0:1]
        m2 = part[:, 1:2]
        out = np.where(mag == m1, m2, np.broadcast_to(m1, mag.shape))
        # ties: only the FIRST occurrence of the minimum may take m2
        first = np.argmin(mag, axis=1)
        out2 = np.broadcast_to(m1, mag.shape).copy()
        out2[np.arange(len(first)), first] = m2[:, 0]
        out = out2
        sgn = np.where(parity[:, None] ^ sg, -1.0, 1.0)
        m = np.where(mask, alpha * sgn * out, 0.0)
        tot[:n] = llr
        tot[n] = 0.0
        np.add.at(tot, idx.ravel(), m.ravel())
        tot[n] = 0.0
        bits = (tot[:n] < 0).astype(np.uint8)
        bb = np.concatenate([bits, [0]]).astype(np.uint8)
        syn = np.bitwise_xor.reduce(np.where(mask, bb[idx], 0), axis=1)
        bad = int(syn.sum())
        if bad == 0:
            return bits, True, it + 1, 0
    return bits, False, iters, bad


# --------------------------------------------------------------------------
# E52 -- the same decoder, over ALL of a Frame's FEC Blocks at once
# --------------------------------------------------------------------------
# `min_sum_decode` runs ~15 NumPy kernels per iteration over one block's
# (4320 x 16) edge array -- small enough that dispatch overhead and cold
# caches dominate, and one Python thread per block leaves every core waiting
# on the same tiny arrays.  The blocks of a Frame are INDEPENDENT codewords,
# so the identical arithmetic runs with a leading (n_blocks, ...) axis: same
# elementwise operations, same per-row reductions, same scatter-add order
# within each row.  Rows leave the batch at the iteration their syndrome
# first reads zero, exactly where the scalar decoder returns.
#
# dtype=float64 is BIT-IDENTICAL to the scalar decoder per row (gated by
# lab/gate_e52.py on real air, not assumed).  dtype=float32 halves the memory
# traffic -- which is what actually bounds this loop -- and is a REAL
# arithmetic change: it is only ever trusted through the e2e gates, which
# compare DECODED BYTES against the float64 reference on the same capture.
# The syndrome stays the oracle either way: a converged row is a valid
# codeword by construction, in either precision.

_BATCH_SCATTER = {}


def _batch_scatter(idx, n):
    """(perm, starts) for an order-preserving segment-sum of the messages.

    perm sorts the flattened edge slots by destination variable, KEEPING flat
    order inside each variable, so summing each segment left-to-right sees the
    same addends in the same order `np.add.at` does.  (float32 path only; the
    float64 path keeps np.add.at itself.)
    """
    key = (id(idx), n)
    hit = _BATCH_SCATTER.get(key)
    if hit is not None:
        return hit[1], hit[2]
    flat = idx.ravel()
    perm = np.argsort(flat, kind="stable").astype(np.int64)
    counts = np.bincount(flat, minlength=n + 1)
    starts = np.concatenate([[0], np.cumsum(counts)[:-1]]).astype(np.int64)
    # every real variable has >= 1 edge (gated by selftest); only the dummy
    # column n can be empty, and reduceat cannot take an out-of-range start.
    # Its segment sum lands in s[:, n], which no reader ever touches.
    starts = np.minimum(starts, max(0, flat.size - 1))
    _BATCH_SCATTER[key] = (idx, perm, starts, counts)
    return perm, starts


def _fast_loop(llr, tot, m, idx, mask, notmask, perm, starts, act,
               bits_out, conv_out, its_out, bad_out, n, iters, al, BIG, dt):
    """The float32 min-sum loop with the syndrome FUSED into the gather.

    The syndrome of the current hard decisions is exactly the sign pattern
    of `tot[:, idx]` -- the same gather the next iteration's `v` needs.  So
    check FIRST on the shared gather, then subtract the messages in place:
    one (na, nc, dmax) fancy gather per iteration instead of two, and a
    block whose CHANNEL hard decisions already satisfy every check (PLP16
    on clean air: 0 raw errors, measured) retires before any BP work at
    all.  E53.  Semantics: a row retires after k message updates exactly as
    in the scalar decoder; only the `its` bookkeeping of the 0-update case
    differs (0 here, 1 there), which nothing gates or ships.
    """
    nc, dmax = idx.shape
    lam_act = llr
    keep = None
    for it in range(iters + 1):
        g = tot[:, idx]                            # syndrome AND v
        syn = np.bitwise_xor.reduce(g < 0, axis=2)  # pad: BIG > 0 -> False
        bad = syn.sum(axis=1, dtype=np.int32)
        done = bad == 0
        if done.any():
            sel = act[done]
            bits_out[sel] = tot[done, :n] < 0
            conv_out[sel] = True
            its_out[sel] = it
            bad_out[sel] = 0
            keep = ~done
            if not keep.any():
                return bits_out, conv_out, its_out, bad_out
            act = act[keep]
            tot = tot[keep]
            m = m[keep]
            lam_act = lam_act[keep]
            g = g[keep]
        else:
            keep = None
        na = len(act)
        if it == iters:
            bits_out[act] = tot[:, :n] < 0
            bad_out[act] = bad[keep] if keep is not None else bad
            break
        v = g
        v -= m
        sg = v < 0
        parity = np.bitwise_xor.reduce(sg, axis=2)
        mag = np.abs(v, out=v)
        part = np.partition(mag, 1, axis=2)
        m1 = part[:, :, 0:1]
        first = np.argmin(mag, axis=2)
        out2 = np.broadcast_to(m1, mag.shape).copy()
        out2[np.arange(na)[:, None], np.arange(nc)[None, :], first] = \
            part[:, :, 1]
        out2 *= al
        np.negative(out2, out=out2, where=parity[:, :, None] ^ sg)
        m = out2
        m[:, notmask] = 0.0
        s = np.add.reduceat(m.reshape(na, -1)[:, perm], starts, axis=1)
        tot[:, :n] = lam_act + s[:, :n]
        tot[:, n] = BIG
    return bits_out, conv_out, its_out, bad_out


# --------------------------------------------------------------------------
# E54 -- optional compiled kernel (lab/ldpc_kernel.c via ctypes)
# --------------------------------------------------------------------------
# E53 measured the numpy batch floor (~130 ms/frame @6T) to be DISPATCH-
# bound, not math- or bandwidth-bound.  The compiled kernel is that lever:
# the same float32 fused-syndrome min-sum as `_fast_loop`, as one C loop
# nest.  Gated output-IDENTICAL (bits, conv, its, bad) to the numpy float32
# path on real air (lab/gate_e54.py), and byte-identical through BCH.
#
# Fleet laws: OPTIONAL -- if lab/ldpc_kernel.{dll,so} is absent (build it
# with lab/build_ldpc_kernel.py, never copy it between boxes) every call
# falls back to the numpy batch below, one log line, and the float64 numpy
# path remains the gate reference either way.  ATSC3_LDPC_KERNEL=0
# disables the kernel outright.  The kernel is single-threaded per call
# and ctypes drops the GIL, so m9_fast's existing thread-pool sub-batching
# parallelizes it with no OpenMP dependency.

import threading as _threading

_KERNEL_STATE = {"lib": None, "logged": False}   # lib: None=untried, False=absent
_KERNEL_LOCK = _threading.Lock()
_KERNEL_TABLES = {}
_KERNEL_MAX_DEG = 128                            # must match ldpc_kernel.c


def _kernel_lib():
    with _KERNEL_LOCK:
        lib = _KERNEL_STATE["lib"]
        if lib is None:
            lib = False
            name = "ldpc_kernel.dll" if os.name == "nt" else "ldpc_kernel.so"
            path = os.path.join(HERE, name)
            if os.environ.get("ATSC3_LDPC_KERNEL", "1") != "0":
                try:
                    import ctypes
                    cand = ctypes.CDLL(path)
                    cand.ldpc_kernel_abi.restype = ctypes.c_int32
                    if int(cand.ldpc_kernel_abi()) != 1:
                        raise OSError("kernel ABI mismatch")
                    i32 = np.ctypeslib.ndpointer(np.int32, flags="C")
                    f32 = np.ctypeslib.ndpointer(np.float32, flags="C")
                    u8 = np.ctypeslib.ndpointer(np.uint8, flags="C")
                    cand.ldpc_minsum_f32.restype = ctypes.c_int32
                    cand.ldpc_minsum_f32.argtypes = [
                        f32, ctypes.c_int32, ctypes.c_int32, ctypes.c_int32,
                        ctypes.c_int32, i32, i32, i32, i32, ctypes.c_int32,
                        ctypes.c_float, u8, u8, i32, i32, f32]
                    lib = cand
                except OSError as ex:
                    if not _KERNEL_STATE["logged"]:
                        _KERNEL_STATE["logged"] = True
                        print(f"[m3_ldpc] compiled LDPC kernel unavailable "
                              f"({name}: {ex}); numpy batch fallback in use "
                              f"-- python lab/build_ldpc_kernel.py",
                              file=sys.stderr)
            _KERNEL_STATE["lib"] = lib
        return lib or None


def _kernel_tables(checks, n):
    """CSR check structure + variable-major edge ids, cached per table.

    vedge lists each variable's edge ids in check-major (flat) order --
    the same addend order as the numpy path's stable-argsort perm, so the
    kernel's per-variable segment sums see what np.add.reduceat sees.
    """
    key = (id(checks), n)
    hit = _KERNEL_TABLES.get(key)
    if hit is not None:
        return hit[1]
    cn = np.array([len(c) for c in checks], np.int32)
    if int(cn.max(initial=0)) > _KERNEL_MAX_DEG:
        tabs = None                       # kernel cannot serve this code
    else:
        cptr = np.concatenate([[0], np.cumsum(cn)]).astype(np.int32)
        cidx = np.fromiter((v for c in checks for v in c), np.int32,
                           int(cn.sum()))
        vedge = np.argsort(cidx, kind="stable").astype(np.int32)
        counts = np.bincount(cidx, minlength=n)
        vptr = np.concatenate([[0], np.cumsum(counts)]).astype(np.int32)
        tabs = (cptr, cidx, vptr, vedge, int(cidx.size))
    _KERNEL_TABLES[key] = (checks, tabs)  # keep `checks` alive: id() is key
    return tabs


def _kernel_decode_checked(lam2d, checks, n, iters, alpha):
    """Kernel decode, or (None,)*4 when the kernel cannot serve the call."""
    lib = _kernel_lib()
    if lib is None:
        return None, None, None, None
    tabs = _kernel_tables(checks, n)
    if tabs is None:
        return None, None, None, None
    cptr, cidx, vptr, vedge, nedges = tabs
    llr = np.ascontiguousarray(lam2d, np.float32)
    nblk = int(llr.shape[0])
    bits = np.empty((nblk, n), np.uint8)
    conv = np.zeros(nblk, np.uint8)
    its = np.empty(nblk, np.int32)
    bad = np.empty(nblk, np.int32)
    work = np.empty(n + nedges, np.float32)
    rc = lib.ldpc_minsum_f32(llr, nblk, n, len(cptr) - 1, nedges,
                             cptr, cidx, vptr, vedge, int(iters),
                             np.float32(alpha), bits, conv, its, bad, work)
    if rc != 0:                                        # pragma: no cover
        return None, None, None, None
    return bits, conv.view(bool), its, bad


def min_sum_decode_batch(lam2d, checks, n, iters=50, alpha=0.75,
                         dtype=np.float64):
    """Batched `min_sum_decode`: (nblk, n) LLRs -> (bits2d, conv, its, bad).

    Per-row semantics identical to the scalar decoder: a row's bits are taken
    at the first iteration whose syndrome is zero; a row that never converges
    returns its final bits and its last unsatisfied-check count.
    """
    dt = np.dtype(dtype)
    exact = dt == np.float64
    lam2d = np.asarray(lam2d)
    nblk = int(lam2d.shape[0])
    if not exact and dt == np.float32 and nblk:
        # E54: the compiled kernel serves the fast path; float64 stays numpy
        # (the gate reference is never routed through the kernel).
        b, c, i, d = _kernel_decode_checked(lam2d, checks, n, iters, alpha)
        if b is not None:
            return b, c, i, d
    idx, mask = pack(checks, n)
    llr = np.ascontiguousarray(lam2d, dtype=dt)
    nc, dmax = idx.shape
    BIG = dt.type(1e9)
    neg1, pos1 = dt.type(-1.0), dt.type(1.0)
    al = float(alpha) if exact else dt.type(alpha)
    notmask = ~mask
    flat = idx.ravel()
    if not exact:
        perm, starts = _batch_scatter(idx, n)

    bits_out = np.zeros((nblk, n), np.uint8)
    conv_out = np.zeros(nblk, bool)
    its_out = np.full(nblk, iters, np.int32)
    bad_out = np.full(nblk, -1, np.int32)
    if nblk == 0:
        return bits_out, conv_out, its_out, bad_out

    act = np.arange(nblk)                       # active row -> original row
    tot = np.empty((nblk, n + 1), dt)
    tot[:, :n] = llr
    # The fast path never writes the padding slots of `v` at all: the dummy
    # variable n HOLDS BIG and the padded message slots hold 0, so the gather
    # itself produces BIG there -- one full-size masked write per iteration
    # gone, exactly (BIG - 0 is exact in either dtype).
    tot[:, n] = 0.0 if exact else BIG
    m = np.zeros((nblk, nc, dmax), dt)
    lam_act = llr
    bits_out[:] = (llr < 0)
    if not exact:
        return _fast_loop(llr, tot, m, idx, mask, notmask, perm, starts,
                          act, bits_out, conv_out, its_out, bad_out,
                          n, iters, al, BIG, dt)
    for it in range(iters):
        na = len(act)
        v = tot[:, idx] - m
        if exact:
            v[:, notmask] = BIG                 # padding must not win the min
            sg = v < 0
            parity = np.bitwise_xor.reduce(sg & mask, axis=2)
            mag = np.abs(v)
        else:
            sg = v < 0                          # padding: BIG -> False
            parity = np.bitwise_xor.reduce(sg, axis=2)
            mag = np.abs(v, out=v)              # v is dead; reuse its memory
        part = np.partition(mag, 1, axis=2)
        m1 = part[:, :, 0:1]
        # ties: only the FIRST occurrence of the minimum may take m2
        first = np.argmin(mag, axis=2)
        out2 = np.broadcast_to(m1, mag.shape).copy()
        out2[np.arange(na)[:, None], np.arange(nc)[None, :], first] = \
            part[:, :, 1]
        if exact:
            sgn = np.where(parity[:, :, None] ^ sg, neg1, pos1)
            m = np.where(mask, al * sgn * out2, dt.type(0.0))
            tot[:, :n] = lam_act
            tot[:, n] = 0.0
            mf = m.reshape(na, -1)
            for i in range(na):
                np.add.at(tot[i], flat, mf[i])
            tot[:, n] = 0.0
        else:
            # same values, fewer temporaries: |m| in place, negate in place
            # where the extrinsic sign says so, zero the padding slots only.
            out2 *= al
            np.negative(out2, out=out2, where=parity[:, :, None] ^ sg)
            m = out2
            m[:, notmask] = 0.0
            # segment-sum in np.add.at's per-variable order; only the
            # (llr + sum) association differs, which is what float32 already
            # concedes -- the decoded-bytes gate is the referee here.
            s = np.add.reduceat(m.reshape(na, -1)[:, perm], starts, axis=1)
            tot[:, :n] = lam_act + s[:, :n]
            tot[:, n] = BIG
        bits = tot[:, :n] < 0
        if exact:
            bb = np.concatenate([bits, np.zeros((na, 1), bool)], axis=1)
            syn = np.bitwise_xor.reduce(bb[:, idx] & mask, axis=2)
        else:
            bb = np.zeros((na, n + 1), bool)
            bb[:, :n] = bits
            syn = np.bitwise_xor.reduce(bb[:, idx], axis=2)
        bad = syn.sum(axis=1, dtype=np.int32)
        done = bad == 0
        if done.any():
            sel = act[done]
            bits_out[sel] = bits[done]
            conv_out[sel] = True
            its_out[sel] = it + 1
            bad_out[sel] = 0
            keep = ~done
            if not keep.any():
                return bits_out, conv_out, its_out, bad_out
            act = act[keep]
            tot = tot[keep]
            m = m[keep]
            lam_act = lam_act[keep]
        else:
            keep = None
        if it == iters - 1:
            live = slice(None) if keep is None else keep
            bits_out[act] = bits[live] if keep is not None else bits
            bad_out[act] = bad[live] if keep is not None else bad
    return bits_out, conv_out, its_out, bad_out


# --------------------------------------------------------------------------
# selftests
# --------------------------------------------------------------------------

def selftest(rates=("3/15",), variant="standard", seed=3, n=16200):
    rng = np.random.default_rng(seed)
    ok = True
    for rate in rates:
        P, T, K, M1, M2, Q1, Q2 = _params(rate, n)
        checks = parity_check(rate, n=n, variant=variant)
        # G1 -- structural: the right number of checks, right column coverage
        g1 = len(checks) == n - K
        cols = np.zeros(n, np.int64)
        for c in checks:
            for v in c:
                cols[v] += 1
        g2 = bool((cols > 0).all())
        # G2 -- encoder vs analytical H, three random messages
        g3 = True
        for _ in range(3):
            s = rng.integers(0, 2, K).astype(np.uint8)
            lam = encode(s, rate, n=n, variant=variant)
            for c in checks:
                if np.bitwise_xor.reduce(lam[c]):
                    g3 = False
                    break
            if not g3:
                break
        # G3 -- systematic
        s = rng.integers(0, 2, K).astype(np.uint8)
        g4 = bool((encode(s, rate, n=n, variant=variant)[:K] == s).all())
        for nm, g in (("check count == N-K", g1),
                      ("every column covered", g2),
                      ("H . encode(s) == 0", g3),
                      ("code is systematic", g4)):
            print(f"  [{'PASS' if g else 'FAIL'}] n={n} {rate} {variant}: {nm}")
            ok &= g
        wts = np.array([len(c) for c in checks])
        print(f"        checks={len(checks)} col-degree "
              f"min/max={cols.min()}/{cols.max()} "
              f"row-degree min/max={wts.min()}/{wts.max()}")
    return ok


if __name__ == "__main__":
    print("A/322 LDPC encoder / parity-check cross-validation\n")
    good = True
    for v in PARITY_INTERLEAVE_VARIANTS:
        good &= selftest(("3/15",), variant=v)
        print()
    if SL64 is not None:
        print("Ninner = 64800 (Annex A.1) -- the codes RF25/RF30 need\n")
        # 9/15 and 6/15 are the two LLS-bearing core-layer rates; 2/15 and 7/15
        # exercise the Type A path at 64800 as well.
        good &= selftest(("2/15", "6/15", "7/15", "9/15"), n=64800)
        print()
    raise SystemExit(0 if good else 1)
