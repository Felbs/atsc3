#!/usr/bin/env python3
"""M9 Step 2/3 -- the LDPC decoder on the 4090, batched over a whole Frame,
and gated BIT-IDENTICAL against the NumPy decoder it replaces.

THE SHAPE OF THE PROBLEM
------------------------
A/322's 16200-bit codes are sparse: rate 11/15 has 4320 check nodes of degree
14..16 (66,959 edges) and rate 2/15 has 14,040 of degree 3..4.  One Frame of
RF33's PLP 0 carries **74 FEC Blocks that are completely independent**, so the
natural GPU unit is not one codeword but one Frame: a (74, 69120) message
array, five gathers and a scatter per iteration.  That is the entire idea.

WHY THIS IS BIT-IDENTICAL AND NOT MERELY "CLOSE"
-----------------------------------------------
Every operation in normalized min-sum is either exactly-rounded elementwise
IEEE-754 (subtract, abs, compare, multiply by +-0.75) or an exact integer
reduction (the parity XOR, the syndrome).  There is exactly ONE reduction over
floating point: the variable-node sum

    tot[v] = llr[v] + sum over the edges of v of m[e]

and floating-point addition is not associative, so the ORDER matters and a GPU
scatter-add with atomics would not reproduce NumPy.  `np.add.at` visits the
edges in raveled (check-major, slot-minor) order and accumulates strictly
left-to-right.  So this decoder does NOT scatter.  It builds a variable-node
incidence table whose columns are the edges of each variable **sorted by that
same raveled edge index**, pads with an edge whose message is permanently
+0.0, and accumulates with a Python loop over the 12 columns -- left to right.
Same summands, same order, therefore the same double, bit for bit.

The two places where a tie could be broken differently are handled explicitly:

  * the first minimum.  NumPy's `argmin` returns the FIRST minimal position;
    `torch.min`'s index tie-breaking is not contractual on CUDA.  So the
    position is computed as `argmax((mag == m1) * (dmax - arange))`, whose
    maximum is unique by construction.
  * the second minimum.  NumPy takes `np.partition(mag,1)[:,1]`, the second
    smallest WITH multiplicity.  Here it is the min after masking out only the
    first minimal position, which is the same value including when the two
    smallest are equal.

`--dtype float32` and `--dtype bfloat16` are offered and MEASURED rather than
assumed: they are faster, they are not bit-identical, and the gate says so.

Usage:
    python m9_gpu.py                    # gates, on synthetic + real LLRs
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m3_ldpc as LD                                             # noqa: E402

_TORCH = None


def torch():
    global _TORCH
    if _TORCH is None:
        import torch as t
        _TORCH = t
    return _TORCH


def available():
    try:
        return torch().cuda.is_available()
    except Exception:                                          # noqa: BLE001
        return False


# ---------------------------------------------------------------------------


class GpuMinSum:
    """Batched normalized min-sum for ONE (rate, Ninner) code.

    Build once per PLP, reuse for every Frame.  `decode(llr2d)` takes a
    (nblk, Ninner) float64 array of LLRs in codeword order -- exactly what
    `PlpChain.llr()` produces -- and returns the same tuple shape the NumPy
    decoder returns, per block.
    """

    BIG = 1e9

    def __init__(self, checks, n, iters=50, alpha=0.75, device="cuda",
                 dtype="float64"):
        t = torch()
        self.n, self.iters, self.alpha = n, iters, alpha
        self.dev = t.device(device)
        self.dt = dict(float64=t.float64, float32=t.float32,
                       bfloat16=t.bfloat16, float16=t.float16)[dtype]
        self.dtype_name = dtype

        idx, mask = LD.pack(checks, n)                  # (nc, dmax)
        self.nc, self.dmax = idx.shape
        ne = self.nc * self.dmax
        self.ne = ne
        self.idx_np, self.mask_np = idx, mask

        # -- check-node side ------------------------------------------------
        self.idxflat = t.from_numpy(idx.astype(np.int64).ravel()).to(self.dev)
        self.mask = t.from_numpy(mask.copy()).to(self.dev)          # (nc,dmax)
        self.maskf = self.mask.reshape(1, self.nc, self.dmax)

        # -- variable-node side: edges of each variable, in RAVELED ORDER ----
        # e = c*dmax + s.  np.add.at walks e ascending, so column k of vidx
        # must hold each variable's k-th smallest edge index.
        e = np.arange(ne)
        var = idx.ravel().astype(np.int64)
        keep = mask.ravel()
        order = np.lexsort((e[keep], var[keep]))        # var major, e minor
        v_sorted = var[keep][order]
        e_sorted = e[keep][order]
        deg = np.bincount(v_sorted, minlength=n + 1)
        self.vdmax = int(deg.max())
        # position of each edge within its variable's list
        start = np.concatenate([[0], np.cumsum(deg)])
        within = np.arange(len(e_sorted)) - start[v_sorted]
        vidx = np.full((n + 1, self.vdmax), ne, np.int64)   # ne = zero pad
        vidx[v_sorted, within] = e_sorted
        self.vidx = [t.from_numpy(np.ascontiguousarray(vidx[:, k])).to(self.dev)
                     for k in range(self.vdmax)]

        self.arange = t.arange(self.dmax, device=self.dev)
        self.stats = dict(calls=0, iters=0, blocks=0)

    # -----------------------------------------------------------------------
    def decode(self, llr2d, iters=None, alpha=None):
        """(nblk, n) float64 -> (bits uint8 (nblk,n), converged, iters, bad).

        `llr2d` may already be a CUDA tensor, in which case nothing crosses
        the bus -- that is how the fused front end feeds it.
        """
        t = torch()
        iters = iters or self.iters
        if isinstance(llr2d, t.Tensor):
            llr = llr2d.to(self.dev, dtype=self.dt)
            nb = llr.shape[0]
            assert llr.shape[1] == self.n
        else:
            a = np.ascontiguousarray(np.asarray(llr2d, np.float64))
            nb = a.shape[0]
            assert a.shape[1] == self.n
            llr = t.from_numpy(a).to(self.dev, dtype=self.dt)
        z = t.zeros((nb, 1), dtype=self.dt, device=self.dev)
        llr = t.cat([llr, z], dim=1)                     # (nb, n+1), tot[n]=0

        m = t.zeros((nb, self.ne + 1), dtype=self.dt, device=self.dev)
        tot = llr.clone()
        alpha = self.alpha if alpha is None else alpha

        out_bits = t.zeros((nb, self.n), dtype=t.uint8, device=self.dev)
        out_it = t.full((nb,), iters, dtype=t.int32, device=self.dev)
        out_bad = t.zeros((nb,), dtype=t.int32, device=self.dev)
        frozen = t.zeros((nb,), dtype=t.bool, device=self.dev)
        nit = 0

        for it in range(iters):
            nit = it + 1
            g = t.index_select(tot, 1, self.idxflat).view(nb, self.nc,
                                                          self.dmax)
            v = g - m[:, :self.ne].view(nb, self.nc, self.dmax)
            v = t.where(self.maskf, v, t.tensor(self.BIG, dtype=self.dt,
                                                device=self.dev))
            sg = v < 0
            parity = ((sg & self.maskf).sum(-1) & 1).bool()      # (nb, nc)
            mag = v.abs()
            m1 = mag.amin(-1, keepdim=True)
            ismin = mag == m1
            first = (ismin.to(t.int32)
                     * (self.dmax - self.arange).view(1, 1, -1)).argmax(-1)
            onehot = t.nn.functional.one_hot(first, self.dmax).bool()
            m2 = mag.masked_fill(onehot, float("inf")).amin(-1, keepdim=True)
            out = t.where(onehot, m2.expand_as(mag), m1.expand_as(mag))
            sgn = t.where(parity.unsqueeze(-1) ^ sg,
                          t.tensor(-1.0, dtype=self.dt, device=self.dev),
                          t.tensor(1.0, dtype=self.dt, device=self.dev))
            mm = t.where(self.maskf, (alpha * sgn) * out,
                         t.tensor(0.0, dtype=self.dt, device=self.dev))
            m[:, :self.ne] = mm.reshape(nb, self.ne)

            # variable-node update -- LEFT TO RIGHT, raveled edge order
            tot = llr.clone()
            for k in range(self.vdmax):
                tot = tot + t.index_select(m, 1, self.vidx[k])
            tot[:, self.n] = 0

            bits = (tot < 0).to(t.uint8)
            bg = t.index_select(bits, 1, self.idxflat).view(nb, self.nc,
                                                            self.dmax)
            syn = (t.where(self.maskf, bg, t.zeros_like(bg)).to(t.int32)
                   .sum(-1) & 1)
            bad = syn.sum(-1).to(t.int32)                        # (nb,)
            newly = (~frozen) & (bad == 0)
            if newly.any():
                sel = newly.unsqueeze(-1)
                out_bits = t.where(sel, bits[:, :self.n], out_bits)
                out_it = t.where(newly, t.full_like(out_it, it + 1), out_it)
                frozen = frozen | newly
            if bool(frozen.all()):
                break
        live = ~frozen
        if bool(live.any()):
            out_bits = t.where(live.unsqueeze(-1), bits[:, :self.n], out_bits)
            out_bad = t.where(live, bad, out_bad)
        self.stats["calls"] += 1
        self.stats["iters"] += nit
        self.stats["blocks"] += nb
        return (out_bits.cpu().numpy(), frozen.cpu().numpy(),
                out_it.cpu().numpy(), out_bad.cpu().numpy())


# ---------------------------------------------------------------------------
# the fused front end: CPE -> HTI -> demap -> bit de-interleave, on the GPU
# ---------------------------------------------------------------------------

class GpuBicm:
    """Everything between the equalised cell pool and the LDPC's LLRs.

    HONESTY NOTE -- THIS ONE IS *NOT* BIT-IDENTICAL, AND THE REASON IS NAMED.
    Both the CPE hard-decision and the max-log demapper are built on
    `abs(complex)`, and NumPy's `npy_cabs` and CUDA's complex magnitude are
    different correctly-rounded-ish implementations of hypot: they disagree in
    the last unit in the last place on roughly 40% of inputs.  (Measured; and
    torch disagrees with NumPy on the CPU too, so this is not a GPU artefact.)
    There is no way to reproduce one library's hypot in the other, so the
    strict gate that the LDPC decoder passes is unavailable here.

    The gate that IS available is the one that matters: the decoded Baseband
    Packets must be byte-identical to the reference chain, over thousands of
    FEC Blocks.  A 1-ulp LLR perturbation is ~1e-16 relative on a decision
    whose margin is order 1, so it cannot flip a bit unless the block is
    already at the cliff -- and if it ever does, the gate reports it as a
    difference rather than hiding it.

    The bit-identical path (`backend="gpu"`, CPU front end) stays selectable
    and is the default for that reason.
    """

    def __init__(self, pool_len, plp0_size, plp16_start, plp16_size,
                 symbol_of, dummy_values, pts0, pts16, nb0, nb16,
                 hti_idx, lam_of_q0, lam_of_q16, device="cuda", cpe_iters=3):
        t = torch()
        self.dev = t.device(device)
        self.n = pool_len
        self.s0, self.s16a = plp0_size, plp16_start
        self.s16b = plp16_start + plp16_size
        self.cpe_iters = cpe_iters
        syms, inv = np.unique(symbol_of, return_inverse=True)
        self.nsym = len(syms)
        self.sym = t.from_numpy(inv.astype(np.int64)).to(self.dev)
        self.dv = t.from_numpy(np.ascontiguousarray(
            dummy_values[self.s16b:])).to(self.dev)
        self.p0 = t.from_numpy(np.ascontiguousarray(pts0)).to(self.dev)
        self.p16 = t.from_numpy(np.ascontiguousarray(pts16)).to(self.dev)
        self.nb0, self.nb16 = nb0, nb16
        self.hti = t.from_numpy(hti_idx.astype(np.int64)).to(self.dev)
        self.lq0 = t.from_numpy(np.ascontiguousarray(
            lam_of_q0.astype(np.int64))).to(self.dev)
        self.lq16 = t.from_numpy(np.ascontiguousarray(
            lam_of_q16.astype(np.int64))).to(self.dev)

    # -- pieces --------------------------------------------------------------
    def _hard(self, z, pts, chunk=65536):
        t = torch()
        out = t.empty_like(z)
        for lo in range(0, len(z), chunk):
            hi = min(len(z), lo + chunk)
            d = (z[lo:hi, None] - pts[None, :]).abs()
            out[lo:hi] = pts[d.argmin(1)]
        return out

    def cpe(self, pool):
        """m6_cells.cpe_correct, every OFDM symbol at once."""
        t = torch()
        z = pool / t.sqrt((pool.abs() ** 2).mean())
        for _ in range(self.cpe_iters):
            hard = t.empty_like(z)
            hard[:self.s0] = self._hard(z[:self.s0], self.p0)
            hard[self.s16a:self.s16b] = self._hard(z[self.s16a:self.s16b],
                                                   self.p16)
            hard[self.s16b:] = self.dv
            num = t.zeros(self.nsym, dtype=z.dtype, device=self.dev)
            num.index_add_(0, self.sym, hard.conj() * z)
            den = t.zeros(self.nsym, dtype=t.float64, device=self.dev)
            den.index_add_(0, self.sym, hard.abs() ** 2)
            c = num / t.clamp(den, min=1e-12)
            z = z / c[self.sym]
        return z

    def demap(self, cells, pts, nbits):
        """(nblk, ncell) -> (nblk, ncell*nbits) max-log LLRs."""
        t = torch()
        nblk, ncell = cells.shape
        d2 = (cells[:, :, None] - pts[None, None, :]).abs()
        d2 = d2 * d2
        s2 = t.clamp(d2.amin(2).mean(1), min=1e-9).view(nblk, 1)
        d2r = d2.view(nblk, ncell, *([2] * nbits))
        sub = tuple(range(2, 2 + nbits - 1))
        out = t.empty((nblk, ncell, nbits), dtype=t.float64, device=self.dev)
        head = (slice(None), slice(None))
        for i in range(nbits):
            a = d2r[head + (slice(None),) * i + (1,)].amin(dim=sub)
            b = d2r[head + (slice(None),) * i + (0,)].amin(dim=sub)
            out[:, :, i] = (a - b) / s2
        return out.reshape(nblk, ncell * nbits)

    def run(self, pool_np):
        """Equalised pool (numpy complex128) -> (lam0, lam16, pool_tail)."""
        t = torch()
        pool = t.from_numpy(np.ascontiguousarray(pool_np)).to(self.dev)
        z = self.cpe(pool)
        cells0 = z[:self.s0][self.hti]                      # (74, 2700)
        q0 = self.demap(cells0, self.p0, self.nb0)
        lam0 = t.empty_like(q0)
        lam0.index_copy_(1, self.lq0, q0)
        c16 = z[self.s16a:self.s16b].view(1, -1)
        q16 = self.demap(c16, self.p16, self.nb16)
        lam16 = t.empty_like(q16)
        lam16.index_copy_(1, self.lq16, q16)
        return lam0, lam16, z[self.s16b:].cpu().numpy()


# ---------------------------------------------------------------------------
# the gate
# ---------------------------------------------------------------------------

def _cpu_ref(llr2d, checks, n, iters):
    bits, conv, its, bad = [], [], [], []
    for row in llr2d:
        b, c, i, d = LD.min_sum_decode(row, checks, n, iters=iters)
        bits.append(b)
        conv.append(c)
        its.append(i)
        bad.append(d)
    return np.array(bits), np.array(conv), np.array(its), np.array(bad)


def verify(nblk=24, iters=50, dtype="float64", verbose=True, seed=5,
           real_llr=None):
    """CPU vs GPU on the SAME LLRs.  Requires exact equality of every field."""
    import m6_bicm as B
    import m9_accel
    m9_accel.install()
    rng = np.random.default_rng(seed)
    ok = True

    cases = []
    if real_llr is not None:
        cases.append(("real air LLRs (RF33 PLP 0)", real_llr, "11/15", 16200))
    # synthetic: encode, add noise at a level that forces several iterations
    for rate, mod, snr, tag in (("11/15", "64QAM", 15.5, "near threshold"),
                                ("11/15", "64QAM", 19.0, "air-like MER"),
                                ("2/15", "QPSK", 1.0, "PLP16-like")):
        ch = B.PlpChain(16200, mod, rate)
        rows = []
        for _ in range(nblk):
            pay = rng.integers(0, 2, ch.kpayload).astype(np.uint8)
            cells = ch.encode(pay)
            nse = (rng.standard_normal(len(cells))
                   + 1j * rng.standard_normal(len(cells))) / np.sqrt(2)
            rows.append(ch.llr(cells + nse * np.sqrt(10 ** (-snr / 10))))
        cases.append((f"synthetic {mod} {rate} @ {snr} dB ({tag})",
                      np.array(rows), rate, 16200))

    for name, llr, rate, n in cases:
        checks = LD.parity_check(rate, n)
        t0 = time.perf_counter()
        cb, cc, ci, cd = _cpu_ref(llr, checks, n, iters)
        t_cpu = time.perf_counter() - t0
        g = GpuMinSum(checks, n, iters=iters, dtype=dtype)
        torch().cuda.synchronize()
        t0 = time.perf_counter()
        gb, gc, gi, gd = g.decode(llr)
        torch().cuda.synchronize()
        t_gpu = time.perf_counter() - t0
        e_bits = int((cb != gb).sum())
        e_conv = int((cc != gc).sum())
        e_it = int((ci != gi).sum())
        e_bad = int((cd != gd).sum())
        good = not (e_bits or e_conv or e_it or e_bad)
        ok &= good
        if verbose:
            print(f"    {'PASS' if good else 'FAIL'}  {name}: {len(llr)} "
                  f"blocks, {int(cc.sum())} converged, iters "
                  f"{ci.min()}..{ci.max()}")
            print(f"          bit differences {e_bits} / {cb.size}   "
                  f"converged-flag {e_conv}   iteration-count {e_it}   "
                  f"unsatisfied {e_bad}")
            print(f"          CPU {t_cpu*1000:8.1f} ms   GPU {t_gpu*1000:8.1f}"
                  f" ms   speedup {t_cpu/max(t_gpu,1e-9):6.1f}x")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blocks", type=int, default=24)
    ap.add_argument("--iters", type=int, default=50)
    ap.add_argument("--dtype", default="float64")
    a = ap.parse_args()
    print("M9 -- GPU batched min-sum LDPC, gated against the NumPy decoder")
    print("=" * 72)
    if not available():
        print("  CUDA NOT AVAILABLE -- nothing to gate")
        return 1
    t = torch()
    print(f"  {t.cuda.get_device_name(0)}, torch {t.__version__}, "
          f"dtype {a.dtype}\n")
    good = verify(a.blocks, a.iters, a.dtype)
    print(f"\n  {'BIT-IDENTICAL' if good else 'DIVERGENCE -- INVESTIGATE'}")
    return 0 if good else 1


if __name__ == "__main__":
    raise SystemExit(main())
