#!/usr/bin/env python3
"""Gate: the array form of the A/322 7.3 frequency-interleaver sequence
against the loop it replaced.

`m3_freqint._interleaving_sequence_impl` is what every decoder path calls;
`_interleaving_sequence_ref` is the loop it was, one step of the spec per
line.  The two must return the same integers, in the same order, with the
same dtype, and fail the same way -- for every FFT size, both wire directions
and both toggle readings, not only the combination the receiver settled on.

    python lab/gate_fi_sequence.py            # ~1 min on 8 processes

Symbols 0..255 are walked for the combination in use (forward, toggle i);
the other three combinations get symbols 0..7, which covers both parities and
four offset words.  Exit status 0 = identical everywhere.
"""
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np                                               # noqa: E402

N_DATA = {8192: (4851, 6913), 16384: (9702, 13825), 32768: (19404, 27649)}


def cases():
    for nfft in sorted(N_DATA):
        for n in N_DATA[nfft]:
            for d in ("forward", "reverse"):
                for tg in ("i", "l"):
                    top = 256 if (d, tg) == ("forward", "i") else 8
                    for l in range(top):
                        yield nfft, l, n, d, tg


def one(case):
    import m3_freqint as FI

    def run(f):
        # a reading of the spec that is not a permutation (toggle "l") is
        # refused with a RuntimeError; the refusal is part of the contract
        t = time.perf_counter()
        try:
            r = f(*case)
        except RuntimeError as e:
            r = str(e)
        return r, time.perf_counter() - t

    a, t_ref = run(FI._interleaving_sequence_ref)
    b, t_new = run(FI._interleaving_sequence_impl)
    if isinstance(a, str) or isinstance(b, str):
        both = isinstance(a, str) and isinstance(b, str)
        return case, both and a == b, None, t_ref, t_new
    same = (a.dtype == b.dtype and a.shape == b.shape
            and bool(np.array_equal(a, b)))
    perm = bool(np.array_equal(np.sort(b), np.arange(case[2])))
    return case, same, perm, t_ref, t_new


def errors_match():
    """Asking for more cells than the FFT has addresses must raise the same
    error, with the same count in it, from both."""
    import m3_freqint as FI
    out = []
    for nfft in sorted(N_DATA):
        msgs = []
        for f in (FI._interleaving_sequence_ref, FI._interleaving_sequence_impl):
            try:
                f(nfft, 3, nfft + 5, "forward", "i")
                msgs.append("no error")
            except RuntimeError as e:
                msgs.append(str(e))
        out.append((nfft, msgs[0] == msgs[1] and msgs[0] != "no error", msgs[0]))
    return out


def main():
    todo = list(cases())
    t0 = time.time()
    bad = notperm = refused = 0
    t_ref = t_new = 0.0
    with ProcessPoolExecutor(max_workers=min(8, os.cpu_count() or 1)) as ex:
        for case, same, perm, tr, tn in ex.map(one, todo, chunksize=16):
            t_ref += tr
            t_new += tn
            if not same:
                bad += 1
                print("  DIFFERENT", case)
            if perm is None:
                refused += 1
            elif not perm:
                notperm += 1
                print("  NOT A PERMUTATION", case)
    print(f"  {len(todo)} cases compared in {time.time() - t0:.0f} s: "
          f"{bad} different, {notperm} returned but not a permutation, "
          f"{refused} refused by both with the same error")
    print(f"  loop {1e3 * t_ref / len(todo):.1f} ms a sequence, array form "
          f"{1e3 * t_new / len(todo):.2f} ms (the first call per FFT size "
          f"walks the LFSR once)")
    errs = errors_match()
    for nfft, ok, msg in errs:
        print(f"  {nfft}: over-long request -> {'same error' if ok else 'ERRORS DIFFER'}: {msg}")
    ok = bad == 0 and notperm == 0 and all(e[1] for e in errs)
    print("  GATE", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
