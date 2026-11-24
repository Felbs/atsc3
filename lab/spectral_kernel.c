/* spectral_kernel.c -- E112: the AC-4 spectral Huffman parse, in C.
 *
 * WHAT THIS IS
 * ------------
 * m28_channels.asf_spectral_data is the single hottest region of the AC-4
 * decoder: profiled on a Raspberry Pi 5 it is 58% of decode_frames, which is
 * itself 72% of the audio worker.  The Python loop reads the bitstream ONE
 * BIT AT A TIME through Bits.u(1) (590k calls per hundred frames) and walks
 * a dict of (length, codeword) per Huffman symbol.  This kernel does the
 * same walk over the same bits with the codebooks flattened into arrays.
 *
 * SEMANTICS (mirrored from the Python, order and all; the gate is the
 * referee and the decoded WAV is the authority):
 *   - bits are MSB-first over a byte buffer, exactly Bits.u;
 *   - per section (cb, s0, s1): cb 0 and cb > 11 are skipped entirely;
 *     k    = offsets[min(s0, noffs-1)]
 *     end  = offsets[min(s1, noffs-1)]
 *     while k < end: decode one codeword, unpack dim values,
 *       - unsigned codebooks (offset 0) then read ONE SIGN BIT PER NONZERO
 *         value, in value order -- a zero line consumes no sign bit;
 *       - codebook 11's escape runs AFTER the sign bits ("the order is
 *         load-bearing", says the reference): |v| == 16 -> count 1-bits as
 *         n, then mag = (1 << (n+4)) + u(n+4), signed by the original v;
 *       - lines[k : k+dim] takes the values, CLIPPED at nlines exactly as
 *         the numpy slice assignment clips;
 *       k += dim  (even when clipped -- the Python does the same).
 *   - any read past the end of the buffer, or a codeword that fails to
 *     resolve within its book's max length, aborts with a nonzero rc and
 *     WITHOUT committing the bit position: the caller re-runs the Python
 *     path from the saved position, so a malformed frame degrades to the
 *     old behaviour instead of a half-parsed one.
 *
 * The codebooks arrive as WALK TREES built once in Python from the same
 * (length, codeword) -> symbol maps the reference decodes with: node i is a
 * pair of int32 (left = bit 0, right = bit 1); a negative entry is a leaf
 * carrying -(symbol+1); zero is "no code here" (the reference would raise).
 * A tree reproduces ANY prefix code exactly -- the real books are not
 * canonical, so no canonical-decode shortcut is safe.
 *
 * Build: python lab/build_spectral_kernel.py
 */

#include <stddef.h>
#include <stdint.h>

#ifdef _WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif

EXPORT int32_t spectral_kernel_abi(void) { return 2; }   /* 2: + scale-factor / SNF group decoders (2026-09-27) */

typedef struct {
    const uint8_t *d;
    int64_t nbits;
    int64_t p;
} bits_t;

static inline int32_t bit1(bits_t *b, int32_t *err)
{
    if (b->p >= b->nbits) { *err = 1; return 0; }
    const int32_t v = (b->d[b->p >> 3] >> (7 - (b->p & 7))) & 1;
    b->p += 1;
    return v;
}

static inline int64_t un(bits_t *b, int32_t n, int32_t *err)
{
    if (b->p + n > b->nbits) { *err = 1; return 0; }
    int64_t v = 0;
    for (int32_t i = 0; i < n; ++i) {
        v = (v << 1) | ((b->d[b->p >> 3] >> (7 - (b->p & 7))) & 1);
        b->p += 1;
    }
    return v;
}

/* buf/nbytes : the bitstream; bitpos: in/out, committed only on rc == 0
 * sects      : nsect triples (cb, s0, s1)
 * offsets    : noffs int32
 * cb_dim/mod/off : indexed by cb, 12 entries each; off == 0 marks unsigned
 * tree       : flattened trees, pairs of int32; troot[cb] = node index of
 *              cb's root (-1 = book absent)
 * lines      : nlines int32, PRE-ZEROED by the caller
 */
EXPORT int32_t asf_spectral_i32(const uint8_t *buf, int64_t nbytes,
                                int64_t *bitpos,
                                const int32_t *sects, int32_t nsect,
                                const int32_t *offsets, int32_t noffs,
                                const int32_t *cb_dim, const int32_t *cb_mod,
                                const int32_t *cb_off,
                                const int32_t *tree, const int32_t *troot,
                                int32_t *lines, int32_t nlines)
{
    if (noffs <= 0 || nlines < 0 || *bitpos < 0) return 2;
    bits_t b = { buf, nbytes * 8, *bitpos };
    int32_t err = 0;

    for (int32_t si = 0; si < nsect; ++si) {
        const int32_t cb = sects[3 * si];
        if (cb == 0 || cb > 11)
            continue;
        const int32_t s0 = sects[3 * si + 1], s1 = sects[3 * si + 2];
        const int32_t dim = cb_dim[cb], mod = cb_mod[cb], off = cb_off[cb];
        const int32_t root = troot[cb];
        if (root < 0 || dim <= 0 || dim > 4) return 3;
        const int32_t uns = (off == 0);
        int32_t k = offsets[s0 < noffs ? s0 : noffs - 1];
        const int32_t end = offsets[s1 < noffs ? s1 : noffs - 1];

        while (k < end) {
            /* one codeword: walk the tree bit by bit */
            int32_t node = root;
            int32_t idx = -1;
            for (int32_t depth = 0; depth < 64; ++depth) {
                const int32_t bit = bit1(&b, &err);
                if (err) return 10;
                const int32_t nxt = tree[2 * node + bit];
                if (nxt == 0) return 11;        /* no codeword matched */
                if (nxt < 0) { idx = -nxt - 1; break; }
                node = nxt;
            }
            if (idx < 0) return 12;

            int32_t vals[4];
            if (dim == 4) {
                const int32_t m2 = mod * mod, m3 = m2 * mod;
                vals[0] = idx / m3 - off;
                vals[1] = (idx / m2) % mod - off;
                vals[2] = (idx / mod) % mod - off;
                vals[3] = idx % mod - off;
            } else {
                vals[0] = idx / mod - off;
                vals[1] = idx % mod - off;
            }
            if (uns) {
                for (int32_t j = 0; j < dim; ++j) {
                    if (vals[j]) {
                        const int32_t s = bit1(&b, &err);
                        if (err) return 13;
                        if (s) vals[j] = -vals[j];
                    }
                }
            }
            if (cb == 11) {
                for (int32_t j = 0; j < dim; ++j) {
                    const int32_t v = vals[j];
                    if (v == 16 || v == -16) {
                        int32_t n = 0;
                        while (bit1(&b, &err)) {
                            if (err) return 14;
                            n += 1;
                            if (n > 40) return 15;
                        }
                        if (err) return 14;
                        const int64_t mag = ((int64_t)1 << (n + 4))
                                            + un(&b, n + 4, &err);
                        if (err) return 16;
                        vals[j] = (int32_t)(v < 0 ? -mag : mag);
                    }
                }
            }
            for (int32_t j = 0; j < dim && k + j < nlines; ++j)
                lines[k + j] = vals[j];
            k += dim;
        }
    }
    *bitpos = b.p;
    return 0;
}

/* ---------------------------------------------------------------------------
 * 2026-09-27: the scale-factor path.  Profiled on a 2-core laptop, the two
 * Python loops below (m28_channels._scalefac_group / _snf_group) were ~40 %
 * of the whole AC-4 decode: one Huffman symbol per band through the Python
 * bit reader, plus a NumPy max() over a handful of lines per band, ~1.3 M
 * calls each for four minutes of audio.  Same tree walk as the spectral
 * kernel above; the arithmetic is Pseudocode 21 unchanged.
 *   out[sfb] = INT32_MIN marks "no scale factor for this band" (Python None).
 *   cur / first are in/out and committed only on rc == 0, like bitpos.
 */
static inline int32_t band_max_i32(const int32_t *lines, int32_t nlines,
                                   const int32_t *offsets, int32_t noffs,
                                   int32_t sfb)
{
    if (sfb + 1 >= noffs) return 0;
    int32_t lo = offsets[sfb], hi = offsets[sfb + 1];
    if (hi > nlines) hi = nlines;
    int32_t m = 0;
    for (int32_t i = lo; i < hi; ++i) {
        const int32_t a = lines[i] < 0 ? -lines[i] : lines[i];
        if (a > m) m = a;
    }
    return m;
}

static inline int32_t walk_one(bits_t *b, const int32_t *tree, int32_t root,
                               int32_t *rc)
{
    int32_t node = root, err = 0;
    for (int32_t depth = 0; depth < 64; ++depth) {
        const int32_t bit = bit1(b, &err);
        if (err) { *rc = 10; return -1; }
        const int32_t nxt = tree[2 * node + bit];
        if (nxt == 0) { *rc = 11; return -1; }
        if (nxt < 0) return -nxt - 1;
        node = nxt;
    }
    *rc = 12;
    return -1;
}

EXPORT int32_t asf_sf_group_i32(const uint8_t *buf, int64_t nbytes,
                                int64_t *bitpos,
                                const int32_t *sfb_cb, int32_t max_sfb,
                                const int32_t *lines, int32_t nlines,
                                const int32_t *offsets, int32_t noffs,
                                const int32_t *tree, int32_t root,
                                int32_t centre, int32_t *cur, int32_t *first,
                                int32_t *out)
{
    if (root < 0 || max_sfb < 0 || *bitpos < 0) return 2;
    bits_t b = { buf, nbytes * 8, *bitpos };
    int32_t c = *cur, f = *first, rc = 0;
    for (int32_t sfb = 0; sfb < max_sfb; ++sfb) {
        out[sfb] = INT32_MIN;
        if (sfb_cb[sfb] != 0 &&
            band_max_i32(lines, nlines, offsets, noffs, sfb) > 0) {
            if (f) {
                const int32_t sym = walk_one(&b, tree, root, &rc);
                if (sym < 0) return rc;
                c += sym - centre;
            } else {
                f = 1;
            }
            out[sfb] = c;
        }
    }
    *bitpos = b.p; *cur = c; *first = f;
    return 0;
}

EXPORT int32_t asf_snf_group_i32(const uint8_t *buf, int64_t nbytes,
                                 int64_t *bitpos,
                                 const int32_t *sfb_cb, int32_t max_sfb,
                                 const int32_t *lines, int32_t nlines,
                                 const int32_t *offsets, int32_t noffs,
                                 const int32_t *tree, int32_t root,
                                 int32_t centre, int32_t *out)
{
    if (root < 0 || max_sfb < 0 || *bitpos < 0) return 2;
    bits_t b = { buf, nbytes * 8, *bitpos };
    int32_t rc = 0;
    for (int32_t sfb = 0; sfb < max_sfb; ++sfb) {
        out[sfb] = INT32_MIN;
        if (sfb_cb[sfb] == 0 ||
            band_max_i32(lines, nlines, offsets, noffs, sfb) == 0) {
            const int32_t sym = walk_one(&b, tree, root, &rc);
            if (sym < 0) return rc;
            out[sfb] = sym - centre;
        }
    }
    *bitpos = b.p;
    return 0;
}

