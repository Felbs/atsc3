/* bch_kernel.c -- E108 (2026-11-07): the batched BCH syndrome in C.
 *
 * m9_accel.bch_syndrome_batch forms, for every block and every parity bit j,
 * the parity of (block words AND mask[j]) -- the same algebra as the LFSR
 * remainder, expressed as a popcount over precomputed masks.  numpy does it
 * as bitwise_count over an (nblk, deg, nw) broadcast: 117 x 192 x 1013
 * uint64 for PLP 1, ~22M popcounts and a 180 MB temporary, 66 ms a Frame on
 * one thread and twice that under load.  This loop does the same popcounts
 * from registers with no temporary.  Integer arithmetic: the parities are
 * EQUAL, not merely close (build script smoke + gate_e108 assert it).
 *
 * words : nblk*nw uint64 (the packed, zero-padded codeword bits, MSB first)
 * mask  : deg*nw uint64 (m9_accel._bch_tables' masks)
 * par   : nblk*deg uint8 out, each 0 or 1
 */
#include <stdint.h>
#include <stddef.h>
#ifdef _WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif
EXPORT int32_t bch_kernel_abi(void) { return 1; }

static inline int pop64(uint64_t x) {
#if defined(__GNUC__)
    return __builtin_popcountll(x);
#else
    x = x - ((x >> 1) & 0x5555555555555555ULL);
    x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
    x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
    return (int)((x * 0x0101010101010101ULL) >> 56);
#endif
}

#if defined(__x86_64__) || defined(__i386__)
/* The builtin only becomes the POPCNT instruction when the function is
 * compiled for it; without this the generic bit-trick fallback ran SLOWER
 * than numpy (measured: 244 ms vs 81 ms on 117 x 64800 bits). */
__attribute__((target("popcnt")))
static int32_t parity_popcnt(const uint64_t *words, int32_t nblk, int32_t nw,
                             const uint64_t *mask, int32_t deg, uint8_t *par)
{
    for (int32_t b = 0; b < nblk; ++b) {
        const uint64_t *w = words + (size_t)b * nw;
        uint8_t *o = par + (size_t)b * deg;
        for (int32_t j = 0; j < deg; ++j) {
            const uint64_t *m = mask + (size_t)j * nw;
            int acc = 0;
            for (int32_t k = 0; k < nw; ++k) acc += __builtin_popcountll(w[k] & m[k]);
            o[j] = (uint8_t)(acc & 1);
        }
    }
    return 0;
}
#endif

EXPORT int32_t bch_parity_u64(const uint64_t *words, int32_t nblk, int32_t nw,
                              const uint64_t *mask, int32_t deg, uint8_t *par)
{
    if (nblk < 0 || nw <= 0 || deg <= 0) return 1;
#if defined(__x86_64__) || defined(__i386__)
    if (__builtin_cpu_supports("popcnt"))
        return parity_popcnt(words, nblk, nw, mask, deg, par);
#endif
    for (int32_t b = 0; b < nblk; ++b) {
        const uint64_t *w = words + (size_t)b * nw;
        uint8_t *o = par + (size_t)b * deg;
        for (int32_t j = 0; j < deg; ++j) {
            const uint64_t *m = mask + (size_t)j * nw;
            int acc = 0;
            for (int32_t k = 0; k < nw; ++k) acc += pop64(w[k] & m[k]);
            o[j] = (uint8_t)(acc & 1);
        }
    }
    return 0;
}
