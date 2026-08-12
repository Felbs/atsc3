#!/usr/bin/env python3
"""M14 -- the dB budget for the multiplex that does NOT decode, computed.

The question this answers is the practical one: "it is about 4 dB short, so
what buys 4 dB?"  The answer turns out to depend on a structural fact that the
phrase "4 dB short" hides, so the budget is computed here rather than asserted
in prose, from measurements that already exist:

    injection level  5.0 dB   SIGNALLED in L1-Detail (Table 9.24)
    channel SNR      3.37 dB  MEASURED against a KNOWN transmitted sequence --
                              the 51257 dummy cells are +-1 with the 5.2.3
                              scrambler's signs, so fitting and subtracting
                              them gives the noise directly.  This is not a
                              constellation fit and does not depend on one.
    threshold        2.5 dB   OUR OWN decoder, QPSK 9/15 Ninner 64800, 3/3
                              trials.  A/327 publishes ~4.6 dB for the mode;
                              both are carried through, because which one is
                              right changes the ANSWER, not just the margin.

THE STRUCTURAL FACT
-------------------
This multiplex is LDM: two layers share the same cells, and the Core layer is
decoded FIRST, with the Enhanced layer treated as interference.  There is no
ordering in which the Enhanced layer can be cancelled before the Core layer is
read -- cancellation runs the other way.

So the Core layer's SINR is

    SINR = P_core / (P_enh + N)

and P_enh does not shrink when the antenna improves.  **The injection level is
therefore a hard ceiling on Core-layer SINR, reachable only as N -> 0.**  That
single fact does two things: it makes the required antenna improvement almost
DOUBLE the apparent shortfall, and it means the whole question of whether the
channel is reachable at all is decided by which threshold is correct.

WHAT THIS SCRIPT IS NOT
-----------------------
It is arithmetic over measurements, not a new measurement.  The one number it
cannot settle -- how much of the loss is recoverable by demapping the Core
layer against the Enhanced layer's FINITE ALPHABET instead of against a
Gaussian -- is named as an experiment, with the harness that would run it
(`m6_bicm.roundtrip` already synthesizes a full chain through AWGN), and is
deliberately NOT given a number here.

Usage:  python m14_budget.py [--inj 5.0] [--snr 3.37]
"""
from __future__ import annotations

import argparse
import math


def budget(inj_db, snr_db, thresholds):
    """-> the Core layer's SINR now, its ceiling, and what each threshold costs.

    Powers are normalised so P_core + P_enh == 1, which is what makes N
    directly comparable to the measured channel SNR.
    """
    r = 10 ** (-inj_db / 10.0)
    p_core, p_enh = 1 / (1 + r), r / (1 + r)
    n = 10 ** (-snr_db / 10.0)
    sinr = p_core / (p_enh + n)
    out = dict(p_core=p_core, p_enh=p_enh, n=n,
               sinr_db=10 * math.log10(sinr),
               ceiling_db=10 * math.log10(p_core / p_enh), rows=[])
    for name, thr in thresholds:
        t = 10 ** (thr / 10.0)
        need_n = p_core / t - p_enh
        row = dict(name=name, thr=thr,
                   sinr_short=thr - out["sinr_db"],
                   headroom=inj_db - thr)
        row["reachable"] = need_n > 0
        row["snr_needed_db"] = (10 * math.log10(n / need_n) if need_n > 0
                                else float("inf"))
        out["rows"].append(row)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(prog="m14 budget")
    ap.add_argument("--inj", type=float, default=5.0,
                    help="signalled LDM injection level, dB")
    ap.add_argument("--snr", type=float, default=3.37,
                    help="channel SNR measured off the dummy cells, dB")
    a = ap.parse_args(argv)
    b = budget(a.inj, a.snr,
               [("our own decoder", 2.5), ("A/327 published", 4.6)])

    print("M14 -- what actually buys the missing dB")
    print("=" * 72)
    print(f"  injection {a.inj:.1f} dB -> P_core {b['p_core']:.4f}, "
          f"P_enh {b['p_enh']:.4f};  channel SNR {a.snr:.2f} dB -> "
          f"N {b['n']:.4f}")
    print(f"  Core-layer SINR now        {b['sinr_db']:+.2f} dB")
    print(f"  CEILING as N -> 0          {b['ceiling_db']:+.2f} dB"
          f"   <- exactly the injection level, and unimprovable by ANY antenna")
    print()
    for r in b["rows"]:
        if not r["reachable"]:
            print(f"  vs {r['name']:16s} ({r['thr']:.1f} dB): UNREACHABLE at "
                  f"any antenna -- threshold is above the ceiling")
            continue
        print(f"  vs {r['name']:16s} ({r['thr']:.1f} dB): SINR short "
              f"{r['sinr_short']:.2f} dB, but needs "
              f"{r['snr_needed_db']:.2f} dB MORE CHANNEL SNR")
        print(f"       {'':18s}  headroom under the ceiling: "
              f"{r['headroom']:.1f} dB")
    print("""
  WHY THE TWO NUMBERS DIFFER, AND WHY IT MATTERS
  ----------------------------------------------
  A 2.15 dB SINR shortfall needs 3.91 dB of antenna because only the N term
  responds to the antenna and N is barely two thirds of the denominator.  Every
  dB of antenna buys LESS than a dB of SINR, and the closer the target sits to
  the ceiling the worse the exchange rate gets.  Against A/327's threshold the
  exchange rate is catastrophic: 0.4 dB of headroom needs ~13 dB of antenna.

  THE MENU, with what each lever is actually worth
  ------------------------------------------------
    antenna / preamp        THE lever -- the only one that moves N at all.
                            ~3.9 dB needed against our own threshold.
    matched filter          0 dB.  For OFDM the FFT already IS the matched
                            filter per subcarrier; that gain is collected.
    more FEC / iterations   0 dB, MEASURED (the marginal-FEC study).  This
                            decoder already runs 50 iterations and this
                            multiplex already needs 16-43 of them, so it is
                            sitting on the waterfall, not short of effort.
    2-D pilot interpolation 0.08 dB, MEASURED.  The estimator is not the limit.
    widely-linear equaliser PREDICTED NEGATIVE here, and worth stating as a
                            prediction: WL wins at a multipath cliff edge and
                            loses in near-AWGN below ~17 dB, and this channel
                            is 3.37 dB of nearly pure AWGN, not a cliff.
    LDM cancellation        0 dB for the Core layer, by ordering -- it helps
                            the Enhanced layer, which is the layer we do not
                            need.
    alphabet-aware LLRs     UNKNOWN, and the only software lever with a real
                            case.  Standard max-log treats the Enhanced layer
                            as Gaussian; it is not, it is a known finite
                            constellation at a known relative power.
                            Marginalising over it gives correctly-scaled LLRs
                            into the LDPC.  NOT given a number here -- it is an
                            experiment: synthesise Core + Enhanced at the
                            measured injection through `m6_bicm.roundtrip`,
                            find each demapper's threshold, and report the
                            difference.  Radio-free, offline, falsifiable.

  THE HONEST SUMMARY
  ------------------
  The "4 dB" is real and the antenna is the lever -- but the reason is the LDM
  interference floor, not a plain link deficit, and that changes the shopping
  decision: the improvement has to land on N, and it has to be bigger than the
  shortfall suggests.  Before buying anything, the alphabet-aware-LLR
  experiment is worth running, because it is free and it is the only candidate
  that has not been measured to be worth zero.""")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
