#!/usr/bin/env python3
"""gate_e122_ldm_roundtrip.py -- the synthesized LDM round trip DESIGN_NOTE asked for.

The question it answers (DESIGN_NOTE, "the LDM canceller was built, and it
stops UNGATED"): when the Enhanced layer fails to decode after the Core is
cancelled, is it buried or is a rung wrong?  47% unsatisfied looks the same
either way.  A round trip proves only that the code inverts itself (M6's law),
but that is exactly the discrimination needed here, and it costs no signal.

What is synthesized, per configuration (RF30's layers and RF25's):
  * Core PLP FEC Blocks through the receiver's own transmit path
    (m6_bicm.PlpChain.encode: payload -> scramble -> BCH -> LDPC -> bit
    interleave -> constellation), and Enhanced PLP Blocks the same way, each
    constellation normalised to unit power;
  * the LDM superposition of A/322 7.2:  x = (core + g * enh) / sqrt(1 + g^2),
    g = 10^(-injection_dB / 20), the Enhanced stream sitting at a block offset
    of its own (its FEC Blocks do not align with the Core's, as on air);
  * seeded complex Gaussian noise at stated Core SNRs;
  * leg B: the whole stream through the shared CTI (m10_cti.interleave at
    RF25's 1024 rows), de-interleaved by the receiver's function, cancelled in
    the de-interleaved domain -- E69's plan.

The receiver side is the live chain's: demap_llr with the right sigma^2, the
layered LDPC kernel with the per-code alpha ladder and the E112 abort, BCH.
Cancellation is E69's: re-encode the DECODED Core payload, fit the gain by
least squares, subtract; the Enhanced layer is demapped from the residual.
Two cancellers are run: ORACLE (re-encode the transmitted payload) and REAL
(re-encode what the Core decoder produced; a Core Block that did not decode is
left in place, so the Enhanced Blocks under it see it as interference).

PASS requires, for each configuration:
  1. noiseless: every Core and every Enhanced Block decodes after REAL
     cancellation, and the Enhanced payloads equal the transmitted bits;
  2. the negative control: re-encoding a WRONG payload leaves a residual of
     at least 0.8 of the received power (cancellation must do nothing);
  3. leg B: the de-interleaved cells equal the pre-CTI stream exactly, and
     the Enhanced layer decodes through it noiselessly;
  4. the noise sweep is PRINTED (Core and Enhanced convergence per SNR, oracle
     and real) with the 50% crossings -- the thresholds the air must clear.
  5. (13.25, --ce-err) a MODELLED channel-estimate error, printed not judged:
     the noiseless stream is multiplied, cell by cell, by (1 + e) with e
     complex Gaussian at a stated power relative to the channel (-30 ... -3
     dB) -- an equaliser whose estimate is wrong by that much -- and the
     table reports what E124 would then measure (P_res/P_core, "Core-layer
     SNR", the implied Enhanced SNR) and whether the Enhanced layer still
     decodes.  Read against the air's E124 numbers it says how much of the
     measured residual an estimation error of each size would account for.

Run: python lab/gate_e122_ldm_roundtrip.py [--config rf30|rf25|both] [--blocks 12] [--ce-err -20,-15,-10]
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)
import m3_ldpc as LD                                                # noqa: E402
import m6_bicm as B                                                 # noqa: E402
import m10_cti as CTI                                               # noqa: E402

NINNER = 64800
CONFIGS = {
    # from the decoded L1 of each carrier (DESIGN_NOTE / 13.20)
    "rf30": dict(core=("QPSK", "6/15"), enh=("64QAM", "6/15"), inj_db=4.0),
    # 13.25: L1D_plp_ldm_injection_level is a 5-bit INDEX (A/331: 0.5 dB steps
    # to 10 dB); the plan printed the raw index 10 and 13.20-13.24 read it as
    # 10 dB.  The signalled level is 5.0 dB (m3_spec.ldm_injection_db, m14).
    "rf25": dict(core=("QPSK", "9/15"), enh=("256QAM", "7/15"), inj_db=5.0),
}
SNRS = (30, 24, 20, 18, 16, 14, 13, 12, 11, 10, 9, 8, 7, 6)


def unit_power(mod, rate):
    pts = B.points_for(mod, rate)
    return float(np.mean(np.abs(pts) ** 2))


def decode_blocks(chain, blocks, sigma2, ex):
    """Receiver: cells -> (converged, bch_ok, payload bits) per block, the live
    chain's rule: layered kernel, per-code alpha ladder, E112 abort."""
    n = chain.ninner
    lam = np.stack([chain.llr(b, sigma2) for b in blocks]).astype(np.float32)
    nb = len(blocks)
    conv = np.zeros(nb, bool)
    bits = np.zeros((nb, n), np.uint8)
    bad = np.full(nb, 1 << 30, np.int32)
    todo = np.arange(nb)
    for a in chain.ALPHA_LADDER:
        if not len(todo):
            break
        sub = lam if len(todo) == nb else lam[todo]
        r = LD.layered_decode_batch(sub, chain.checks, n, iters=chain.iters, alpha=a, ex=ex,
                                    corr=LD.BOXPLUS_DEFAULT)
        if r is None:
            raise SystemExit("layered kernel unavailable: python lab/build_ldpc_layered_kernel.py")
        b, c, _its, d = r
        take = c | (d < bad[todo])
        sel = todo[take]
        conv[sel] = c[take]; bits[sel] = b[take]; bad[sel] = d[take]
        todo = todo[~c]
    ok = np.zeros(nb, bool)
    payload = [None] * nb
    for i in range(nb):
        if conv[i]:
            bb = chain.baseband_packet(bits[i])
            ok[i] = bool(bb["bch_ok"])
            payload[i] = bb["bits"]
    return conv, ok, payload


def cancel(rx_scaled, core_blocks_cells, ncell_c, which, cache, chain_c, payload_tx, payload_rx, ok_rx):
    """E69's canceller on the scaled stream: for each Core Block, re-encode the
    payload (`which` = 'oracle': transmitted; 'real': decoded, if BCH-clean),
    fit the gain, subtract.  Returns the residual, the mean |gain| and how many
    Core Blocks were cancelled."""
    res = rx_scaled.copy()
    gains = []
    n_c = 0
    for k in range(len(core_blocks_cells)):
        if which == "oracle":
            p = payload_tx[k]
        else:
            if not ok_rx[k]:
                continue
            p = payload_rx[k]
        key = p.tobytes()
        s = cache.get(key)
        if s is None:
            s = cache[key] = np.asarray(chain_c.encode(p), np.complex128)
        lo = k * ncell_c
        r = res[lo:lo + ncell_c]
        gain = np.vdot(s, r) / np.vdot(s, s)
        res[lo:lo + ncell_c] = r - gain * s
        gains.append(abs(gain)); n_c += 1
    return res, (float(np.mean(gains)) if gains else 0.0), n_c


def run_config(name, cfg, nblk, seed, ex, snrs=SNRS, verbose=True, ce_err=()):
    rng = np.random.default_rng(seed)
    (mc, rc), (me, re_), inj = cfg["core"], cfg["enh"], cfg["inj_db"]
    chc = B.PlpChain(NINNER, mc, rc, iters=50)
    che = B.PlpChain(NINNER, me, re_, iters=50)
    pc, pe = unit_power(mc, rc), unit_power(me, re_)
    g = 10 ** (-inj / 20.0)
    norm = np.sqrt(1 + g * g)
    ncell_c, ncell_e = chc.cells_per_fec, che.cells_per_fec
    L = nblk * ncell_c
    # -- transmit: Core Blocks back to back; Enhanced Blocks at their own offset
    pay_c = [rng.integers(0, 2, chc.kpayload, dtype=np.uint8) for _ in range(nblk)]
    core = np.concatenate([np.asarray(chc.encode(p), np.complex128) for p in pay_c]) / np.sqrt(pc)
    off = 1234 % ncell_e
    me_n = (L - off) // ncell_e
    pay_e = [rng.integers(0, 2, che.kpayload, dtype=np.uint8) for _ in range(me_n)]
    enh = np.zeros(L, np.complex128)
    enh[off:off + me_n * ncell_e] = np.concatenate(
        [np.asarray(che.encode(p), np.complex128) for p in pay_e]) / np.sqrt(pe)
    x = (core + g * enh) / norm                       # A/322 7.2: unit-power superposition
    cache_tx = {}
    print(f"\n== {name}: Core {mc} {rc} ({ncell_c} cells/Block, {nblk} Blocks), "
          f"Enhanced {me} {re_} ({ncell_e} cells/Block, {me_n} Blocks at offset {off}), "
          f"injection {inj:.1f} dB (g = {g:.4f}), ladders core {chc.ALPHA_LADDER} enh {che.ALPHA_LADDER}, "
          f"box-plus {int(LD.BOXPLUS_DEFAULT)}")

    def receive(rx, label, snr_db):
        """Core decode, both cancellations, Enhanced decode.  -> dict of counts."""
        t = time.perf_counter()
        sn2 = 0.0 if snr_db is None else 10 ** (-snr_db / 10.0)   # noise var per unit-power Core
        y = rx * norm                                              # back to Core units (+ g*enh + n')
        # Core demap: everything that is not Core is interference, in constellation units
        sig_c = (g * g + sn2 * norm * norm) * pc
        cb = [y[k * ncell_c:(k + 1) * ncell_c] * np.sqrt(pc) for k in range(nblk)]
        conv_c, ok_c, pay_rx = decode_blocks(chc, cb, max(sig_c, 1e-9), ex)
        out = dict(core_conv=int(conv_c.sum()), core_ok=int(ok_c.sum()))
        for which in ("oracle", "real"):
            res, gmean, ncanc = cancel(y * np.sqrt(pc), cb, ncell_c, which, cache_tx, chc, pay_c, pay_rx, ok_c)
            # Enhanced units: residual / (g * sqrt(pc)) * sqrt(pe); noise var there = sn2*norm^2/g^2 * pe
            e_units = res / (g * np.sqrt(pc)) * np.sqrt(pe)
            sig_e = (sn2 * norm * norm / (g * g)) * pe
            eb = [e_units[off + m * ncell_e: off + (m + 1) * ncell_e] for m in range(me_n)]
            conv_e, ok_e, pay_e_rx = decode_blocks(che, eb, max(sig_e, 1e-9), ex)
            same = sum(1 for m in range(me_n) if ok_e[m] and np.array_equal(pay_e_rx[m], pay_e[m]))
            out[which] = dict(enh_conv=int(conv_e.sum()), enh_ok=int(ok_e.sum()), enh_payload_match=same,
                              gain=gmean, cancelled=ncanc,
                              p_res=float(np.mean(np.abs(res) ** 2)), p_rx=float(np.mean(np.abs(y * np.sqrt(pc)) ** 2)))
        out["ms"] = 1000 * (time.perf_counter() - t)
        return out

    verdicts = []
    # ---- 1. noiseless, leg A
    r0 = receive(x, "noiseless", None)
    ok1 = (r0["core_ok"] == nblk and r0["real"]["enh_payload_match"] == me_n)
    verdicts.append(("noiseless: every Core and Enhanced Block decodes after REAL cancellation, payloads equal",
                     ok1, f"core {r0['core_ok']}/{nblk}, enhanced {r0['real']['enh_payload_match']}/{me_n} "
                          f"(oracle {r0['oracle']['enh_payload_match']}/{me_n}), |gain| {r0['real']['gain']:.4f}, "
                          f"residual/received power {r0['real']['p_res']/r0['real']['p_rx']:.4f} "
                          f"(expected g^2/(1+g^2) = {g*g/(1+g*g):.4f})"))
    # ---- 2. the negative control: a WRONG payload must not cancel
    y = x * norm * np.sqrt(pc)
    bad_pay = [p.copy() for p in pay_c]
    for p in bad_pay:
        p[::7] ^= 1
    res_bad, gbad, _ = cancel(y, [None] * nblk, ncell_c, "oracle", {}, chc, bad_pay, None, None)
    ratio = float(np.mean(np.abs(res_bad) ** 2) / np.mean(np.abs(y) ** 2))
    ok2 = ratio >= 0.8
    verdicts.append(("negative control: a wrong payload leaves >= 0.8 of the received power",
                     ok2, f"residual/received {ratio:.3f}, fitted |gain| {gbad:.4f}"))
    # ---- 3. leg B: through the shared CTI at RF25's depth, de-interleaved, cancelled there
    nrows = 1024
    reach = nrows * (nrows - 1)                        # the delay-line fill, a lead-in of filler
    lead = int(np.ceil(reach / ncell_c)) * ncell_c
    # a cell entering at i leaves at i + nrows*((start_row+i) % nrows): the test
    # Blocks need filler BEFORE them (the delay lines' initial state) and AFTER
    # them (their own cells leave up to nrows^2 later), so both ends get a lead
    filler = lambda n: (rng.standard_normal(n) + 1j * rng.standard_normal(n)) / np.sqrt(2)
    stream = np.concatenate([filler(lead), x, filler(lead)])
    start_row = 4
    tx = CTI.interleave(stream, nrows, start_row)
    deint, valid = CTI.deinterleave(tx, nrows, start_row, n_out=len(stream))
    seg = slice(lead, lead + L)
    exact = bool(np.all(valid[seg]) and np.array_equal(deint[seg], stream[seg]))
    rB = receive(np.asarray(deint[seg], np.complex128), "cti", None) if exact else None
    ok3 = exact and rB is not None and rB["real"]["enh_payload_match"] == me_n
    verdicts.append(("leg B (CTI 1024 rows): de-interleaved cells exact, Enhanced decodes noiselessly",
                     ok3, f"exact {exact}" + (f", enhanced {rB['real']['enh_payload_match']}/{me_n}" if rB else "")))
    # ---- 4. the noise sweep
    print(f"  {'Core SNR dB':>11s} {'core ok':>8s} {'enh ok oracle':>13s} {'enh ok real':>11s} {'enh match':>9s} {'|gain|':>7s}  ms")
    cross_c = cross_e = None
    for snr in snrs:
        noise = (rng.standard_normal(L) + 1j * rng.standard_normal(L)) * np.sqrt(10 ** (-snr / 10.0) / 2)
        r = receive(x + noise / norm, f"{snr} dB", snr)   # noise var relative to unit-power Core: 10^(-snr/10)
        print(f"  {snr:11d} {r['core_ok']:5d}/{nblk:<3d} {r['oracle']['enh_ok']:8d}/{me_n:<3d} "
              f"{r['real']['enh_ok']:7d}/{me_n:<3d} {r['real']['enh_payload_match']:9d} {r['real']['gain']:7.4f}  {r['ms']:.0f}")
        if cross_c is None and r["core_ok"] < 0.5 * nblk:
            cross_c = snr
        if cross_e is None and r["oracle"]["enh_ok"] < 0.5 * me_n:
            cross_e = snr
    print(f"  50% crossings: Core below {cross_c} dB; Enhanced (oracle cancellation) below {cross_e} dB "
          f"-> Enhanced needs about {cross_e} dB of Core SNR, i.e. {cross_e - inj if cross_e else '-'} dB on its own layer")
    # ---- 5. a modelled channel-estimate error on the noiseless stream (13.25)
    if ce_err:
        print(f"  modelled estimate error (noiseless air, per-cell complex Gaussian (1 + e) on the stream):")
        print(f"  {'err dB':>7s} {'P_res/P_core':>12s} {'Core SNR':>9s} {'Enh SNR':>8s} {'core ok':>8s} "
              f"{'enh ok real':>11s} {'enh match':>9s} {'|gain|':>7s}")
        p_enh_share = g * g                                 # the Enhanced layer's own share of P_core
        for lvl in ce_err:
            sig = np.sqrt(10 ** (lvl / 10.0) / 2)
            err = (rng.standard_normal(L) + 1j * rng.standard_normal(L)) * sig
            r = receive(x * (1 + err), f"ce {lvl} dB", None)
            # E124's arithmetic on the REAL cancellation: P_core = |gain|^2 x Core power (pc), residual as measured
            p_core_m = (r["real"]["gain"] ** 2) * pc
            rat = r["real"]["p_res"] / max(p_core_m, 1e-30)
            p_noise = rat - p_enh_share
            enh_snr = (10 * np.log10(p_enh_share / p_noise) if p_noise > 0 else float("inf"))
            print(f"  {lvl:7.1f} {rat:12.4f} {10*np.log10(1/rat):9.2f} {enh_snr:8.2f} {r['core_ok']:5d}/{nblk:<3d} "
                  f"{r['real']['enh_ok']:7d}/{me_n:<3d} {r['real']['enh_payload_match']:9d} {r['real']['gain']:7.4f}")
        print(f"  (no error: P_res/P_core = g^2 = {p_enh_share:.4f}, Core SNR {10*np.log10(1/p_enh_share):.2f} dB = the injection)")
    ok = all(v[1] for v in verdicts)
    for text, good, detail in verdicts:
        print(f"  [{'PASS' if good else 'FAIL'}] {text}\n         {detail}")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="both", choices=("rf30", "rf25", "both"))
    ap.add_argument("--blocks", type=int, default=12, help="Core Blocks per trial")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--snrs", default=None, help="comma list of Core SNRs in dB (default the sweep)")
    ap.add_argument("--ce-err", default=None,
                    help="13.25: comma list of modelled channel-estimate error levels, dB relative to the channel")
    a = ap.parse_args()
    snrs = SNRS if a.snrs is None else tuple(int(x) for x in a.snrs.split(","))
    ce_err = () if a.ce_err is None else tuple(float(x) for x in a.ce_err.split(","))
    ex = ThreadPoolExecutor(max_workers=a.threads)
    names = list(CONFIGS) if a.config == "both" else [a.config]
    ok = True
    for nm in names:
        ok &= run_config(nm, CONFIGS[nm], a.blocks, a.seed, ex, snrs=snrs, ce_err=ce_err)
    ex.shutdown(wait=True)
    print("\n  OVERALL:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
