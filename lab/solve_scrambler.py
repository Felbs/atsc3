#!/usr/bin/env python3
"""solve_scrambler.py - settle A/322 5.2.3 from the PRINTED TEXT, not the figure.

The M3 run hit a wall reading Figure 5.6 (vector art; the text extraction keeps
only mangled glyph columns, and column spacing in an extracted figure is NOT
faithful geometry). But 5.2.3's prose carries everything unambiguously:

  seed        0xF180  ("loaded into the shift register")
  polynomial  G(x) = 1 + X + X^3 + X^6 + X^7 + X^11 + X^12 + X^13 + X^16
  cadence     "The bits in the shift register shall be shifted once."
  output      D7,D6,...,D0 as a byte, MSB first
  vector      1100 0000  0110 1101  0011 1111  ...

So: do NOT guess the D taps from pixel columns. SOLVE for them against the
printed vector. For each stage we compute its value at t=0,1,2; a D position
is satisfied by any stage whose (v0,v1,v2) equals the bits that position must
emit. If exactly one stage fits each position, the figure is decoded - by
arithmetic, from the spec's own test vector.
"""
SEED_HEX = 0xF180
POLY = [1, 3, 6, 7, 11, 12, 13, 16]          # printed G(x) exponents
VEC = "110000000110110100111111"             # printed, MSB-first D7..D0


def load_seed():
    """Figure 5.6's glyph row reads 0000000110001111 = 0xF180 MSB-first,
    REVERSED -> X1 is the LSB end, X16 the MSB end. Verified below."""
    bits = [(SEED_HEX >> i) & 1 for i in range(16)]     # bits[0] = LSB = X1
    return {i + 1: bits[i] for i in range(16)}


def step(R, direction, fb_taps):
    fb = 0
    for t in fb_taps:
        fb ^= R[t]
    N = {}
    if direction == "up":            # X_{i+1} <- X_i, X1 <- feedback
        for i in range(16, 1, -1):
            N[i] = R[i - 1]
        N[1] = fb
    else:                            # X_i <- X_{i+1}, X16 <- feedback
        for i in range(1, 16):
            N[i] = R[i + 1]
        N[16] = fb
    return N


def stage_histories(direction, fb_taps, n=3):
    R = load_seed()
    hist = {i: [] for i in range(1, 17)}
    for _ in range(n):
        for i in range(1, 17):
            hist[i].append(R[i])
        R = step(R, direction, fb_taps)
    return hist


def solve(direction, fb_taps):
    hist = stage_histories(direction, fb_taps)
    need = {}                                   # D index -> required (v0,v1,v2)
    for t in range(3):
        byte = VEC[t * 8:(t + 1) * 8]           # MSB first = D7..D0
        for j, ch in enumerate(byte):
            d = 7 - j
            need.setdefault(d, []).append(int(ch))
    out = {}
    for d in range(8):
        want = tuple(need[d])
        fits = [s for s in range(1, 17) if tuple(hist[s]) == want]
        out[d] = fits
    return out


def verify(direction, fb_taps, taps, nbytes=3):
    """Full forward check: generate the byte stream and compare to VEC."""
    R = load_seed()
    bits = ""
    for _ in range(nbytes):
        for d in range(7, -1, -1):
            bits += str(R[taps[d]])
        R = step(R, direction, fb_taps)
    return bits


if __name__ == "__main__":
    seed = load_seed()
    glyph = "0000000110001111"                  # what Figure 5.6's row shows
    got = "".join(str(seed[i]) for i in range(1, 17))
    print(f"seed check: figure row {glyph} vs 0xF180-reversed {got} "
          f"-> {'MATCH' if got == glyph else 'MISMATCH'}\n")

    for direction in ("up", "down"):
        res = solve(direction, POLY)
        uniq = all(len(v) >= 1 for v in res.values())
        print(f"--- shift '{direction}', feedback = printed G(x) {POLY}")
        for d in range(7, -1, -1):
            names = ["X%d" % s for s in res[d]]
            print(f"   D{d}: {names if names else 'NO STAGE FITS'}")
        if uniq:
            # try every combination of the fitting stages, verify forward
            import itertools
            cands = [res[d] for d in range(8)]
            ok = []
            for combo in itertools.product(*cands):
                if len(set(combo)) != 8:
                    continue
                taps = {d: combo[d] for d in range(8)}
                if verify(direction, POLY, taps) == VEC:
                    ok.append(taps)
            print(f"   -> {len(ok)} tap set(s) reproduce all 24 printed bits")
            for t in ok[:4]:
                print("      D7..D0 =",
                      [f"X{t[d]}" for d in range(7, -1, -1)])
        print()
