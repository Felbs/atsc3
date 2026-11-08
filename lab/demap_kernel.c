/* demap_kernel.c -- E102: compiled max-log BICM demapper, float32, per cell.
 *
 * WHAT THIS IS
 * ------------
 * The inner loop of m9_fast.FrameDemod.demap_batch's cpu_fast path, written
 * as plain C.  The numpy path is not dispatch-bound like E53's LDPC was --
 * it is MEMORY-bound: it materialises q = (k*ncell, npts) float32 (for 64QAM
 * at chunk 8: 8*2700*64*4 = 5.5 MB, and 256QAM is 4x that), then walks the
 * whole array again for `+= p2f`, once for `q.min(2)`, and twice more per
 * label bit for the two subset minima.  Nine streaming passes over an array
 * that cannot live in cache.
 *
 * This kernel never materialises q at all.  Each cell's npts distances are
 * computed into registers/L1, the twelve subset minima (nb bits x {0,1}) and
 * the global minimum are accumulated in the same sweep, and only the nb
 * differences plus one d2 scalar are written out.  DRAM traffic per cell
 * falls from ~npts*4 bytes written + re-read several times to nb*4 + 4 bytes
 * written once.  Measured on a Raspberry Pi 5 (4x Cortex-A76) -- see E102 in
 * the campaign notes.
 *
 * SEMANTICS (mirrored from the numpy fast path; the gate is the referee):
 *   - q_j = zr*pfx[j] + zi*pfy[j] + p2f[j], evaluated in exactly that
 *     association.  pfx/pfy arrive ALREADY scaled by -2 (E101 folds the
 *     constant into the gemm matrix), so this is |p_j|^2 - 2 Re(z conj p_j),
 *     the |z|^2-shifted squared distance.
 *   - d2 = (zr*zr + zi*zi) + min_j q_j   -- the true squared distance to the
 *     nearest point, written out so the CALLER computes sigma^2 with numpy's
 *     own pairwise mean.  The statistics deliberately stay in numpy: a C
 *     sequential sum would not match np.mean's pairwise association, and
 *     sigma^2 divides every LLR.
 *   - out[i] = min_{j: bit i of j == 1} q_j  -  min_{j: bit i of j == 0} q_j.
 *     A/322 6.3.4.2 numbers label bits MSB-first, so bit i of point j is
 *     (j >> (nb-1-i)) & 1 -- the same mapping the numpy path gets from
 *     reshaping the last axis to (2,)*nb and taking axis i.
 *     Minima are order-independent (no NaNs in this data), so accumulating
 *     them in one sweep gives the same values as numpy's separate reductions.
 *   - The caller applies `/ sigma^2` (or not, for the E60 per-cell-weight
 *     lever) exactly as before, on the values this kernel produced.
 *
 * DELIBERATELY NOT HERE: threads (demap_batch's ThreadPoolExecutor already
 * sub-batches by block and ctypes drops the GIL), Python.h (pure C ABI, so
 * one build serves any CPython on that box -- but never copy the binary
 * between boxes), and -ffast-math (the build pins -ffp-contract=off so no
 * FMA contraction can change a float32 rounding).
 *
 * Build: python lab/build_demap_kernel.py
 */

#include <stddef.h>
#include <stdint.h>

#ifdef _WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif

#define DK_MAX_PTS  256         /* 256QAM is the largest A/322 constellation */
#define DK_MAX_BITS 8           /* log2(256) */

/* ABI 2 adds the CPE nearest-point search; a stale ABI-1 build is
 * refused so the loader falls back to numpy rather than missing a
 * symbol at argtypes time. */
EXPORT int32_t demap_kernel_abi(void) { return 2; }

/* Demap k*ncell cells.
 *
 * zr, zi : k*ncell float32, C-contiguous -- cell real/imag parts
 * k      : blocks in this chunk
 * ncell  : cells per block
 * npts   : constellation size (<= DK_MAX_PTS)
 * nb     : label bits (<= DK_MAX_BITS), npts == 1 << nb
 * pfx    : npts float32 -- constellation real parts, ALREADY scaled by -2
 * pfy    : npts float32 -- constellation imag parts, ALREADY scaled by -2
 * p2f    : npts float32 -- |p_j|^2
 * out    : k*ncell*nb float32 -- receives the raw bitwise-minima differences
 * d2     : k*ncell   float32 -- receives |z|^2 + min_j q_j
 *
 * Returns 0 on success, non-zero if the shape is outside the kernel's range
 * (the caller then falls back to the numpy path).
 */
static int32_t demap_scalar(const float *zr, const float *zi,
                            int32_t k, int32_t ncell,
                            int32_t npts, int32_t nb,
                            const float *pfx, const float *pfy,
                            const float *p2f,
                            float *out, float *d2)
{

    /* Which subset each point belongs to, per label bit.  Hoisted out of the
     * cell loop: it depends only on the constellation, and recomputing the
     * shift per point per cell is the one avoidable integer cost in here. */
    uint8_t bit_of[DK_MAX_PTS][DK_MAX_BITS];
    for (int32_t j = 0; j < npts; ++j)
        for (int32_t i = 0; i < nb; ++i)
            bit_of[j][i] = (uint8_t)((j >> (nb - 1 - i)) & 1);

    const size_t ncells = (size_t)k * (size_t)ncell;
    for (size_t c = 0; c < ncells; ++c) {
        const float x = zr[c], y = zi[c];

        float amin[DK_MAX_BITS], bmin[DK_MAX_BITS];
        for (int32_t i = 0; i < nb; ++i) {
            amin[i] = 3.402823466e+38f;      /* FLT_MAX, no <float.h> needed */
            bmin[i] = 3.402823466e+38f;
        }
        float qmin = 3.402823466e+38f;

        for (int32_t j = 0; j < npts; ++j) {
            /* exactly ((x*pfx) + (y*pfy)) + p2f -- see SEMANTICS */
            const float qj = x * pfx[j] + y * pfy[j] + p2f[j];
            if (qj < qmin) qmin = qj;
            const uint8_t *bo = bit_of[j];
            for (int32_t i = 0; i < nb; ++i) {
                if (bo[i]) { if (qj < amin[i]) amin[i] = qj; }
                else       { if (qj < bmin[i]) bmin[i] = qj; }
            }
        }

        d2[c] = x * x + y * y + qmin;
        float *o = out + c * (size_t)nb;
        for (int32_t i = 0; i < nb; ++i)
            o[i] = amin[i] - bmin[i];
    }
    return 0;
}

/* E106 (2026-11-07): the same sweep, eight constellation points per AVX2
 * vector, branch-free.  Per lane the arithmetic is the scalar path's, in the
 * scalar path's association (two multiplies, two adds, never contracted:
 * the build pins -ffp-contract=off and these are explicit mul/add
 * intrinsics), and a minimum is the same number whatever order it is taken
 * in -- so the output is BIT-IDENTICAL to demap_scalar, which
 * gate_e102_demap.py re-asserts against the numpy path.  Label bit i of
 * point j is (j >> (nb-1-i)) & 1: for the top bits that is a property of
 * the whole 8-point vector (j0 >> s), for the bottom three it is a fixed
 * lane pattern inside every vector.  Chosen at run time by
 * __builtin_cpu_supports, so one build of this file runs anywhere. */
#if defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
__attribute__((target("avx2")))
static int32_t demap_avx2(const float *zr, const float *zi,
                          int32_t k, int32_t ncell,
                          int32_t npts, int32_t nb,
                          const float *pfx, const float *pfy,
                          const float *p2f,
                          float *out, float *d2)
{
    const __m256 inf = _mm256_set1_ps(3.402823466e+38f);
    const int32_t nv = npts / 8;
    const size_t ncells = (size_t)k * (size_t)ncell;
    for (size_t c = 0; c < ncells; ++c) {
        const float x = zr[c], y = zi[c];
        const __m256 xv = _mm256_set1_ps(x), yv = _mm256_set1_ps(y);
        __m256 amin[DK_MAX_BITS], bmin[DK_MAX_BITS], qmin = inf;
        for (int32_t i = 0; i < nb; ++i) { amin[i] = inf; bmin[i] = inf; }
        /* The three LOW label bits of point j depend only on j mod 8 -- the
         * lane -- so their subset minima are lane subsets of ONE lane-wise
         * minimum over all vectors (qmin), taken once per cell.  Only the
         * nb-3 HIGH bits, which depend on the vector index, cost a minimum
         * per vector.  (MEASURED: 270 -> ~120 ms a Frame single-threaded.) */
        const int32_t ntop = nb - 3;
        for (int32_t v = 0; v < nv; ++v) {
            const int32_t j0 = v * 8;
            const __m256 q = _mm256_add_ps(
                _mm256_add_ps(_mm256_mul_ps(xv, _mm256_loadu_ps(pfx + j0)),
                              _mm256_mul_ps(yv, _mm256_loadu_ps(pfy + j0))),
                _mm256_loadu_ps(p2f + j0));
            qmin = _mm256_min_ps(qmin, q);
            for (int32_t i = 0; i < ntop; ++i) {
                if ((v >> (ntop - 1 - i)) & 1) amin[i] = _mm256_min_ps(amin[i], q);
                else                           bmin[i] = _mm256_min_ps(bmin[i], q);
            }
        }
        float t[8];
        _mm256_storeu_ps(t, qmin);
        float qm = t[0];
        for (int l = 1; l < 8; ++l) if (t[l] < qm) qm = t[l];
        d2[c] = x * x + y * y + qm;
        float *o = out + c * (size_t)nb;
        for (int32_t i = 0; i < ntop; ++i) {
            float a, b, u[8];
            _mm256_storeu_ps(u, amin[i]);
            a = u[0]; for (int l = 1; l < 8; ++l) if (u[l] < a) a = u[l];
            _mm256_storeu_ps(u, bmin[i]);
            b = u[0]; for (int l = 1; l < 8; ++l) if (u[l] < b) b = u[l];
            o[i] = a - b;
        }
        for (int32_t i = ntop; i < nb; ++i) {
            const int32_t sh = nb - 1 - i;          /* 2, 1, 0 */
            float a = 3.402823466e+38f, b = 3.402823466e+38f;
            for (int l = 0; l < 8; ++l) {
                if ((l >> sh) & 1) { if (t[l] < a) a = t[l]; }
                else               { if (t[l] < b) b = t[l]; }
            }
            o[i] = a - b;
        }
    }
    return 0;
}
#endif

EXPORT int32_t demap_llr_f32(const float *zr, const float *zi,
                             int32_t k, int32_t ncell,
                             int32_t npts, int32_t nb,
                             const float *pfx, const float *pfy,
                             const float *p2f,
                             float *out, float *d2)
{
    if (npts <= 0 || npts > DK_MAX_PTS || nb <= 0 || nb > DK_MAX_BITS)
        return 1;
    if ((int32_t)1 << nb != npts)
        return 2;
    if (k <= 0 || ncell <= 0)
        return 0;
#if defined(__x86_64__) || defined(__i386__)
    if (npts >= 8 && (npts % 8) == 0 && nb >= 3 && __builtin_cpu_supports("avx2"))
        return demap_avx2(zr, zi, k, ncell, npts, nb, pfx, pfy, p2f, out, d2);
#endif
    return demap_scalar(zr, zi, k, ncell, npts, nb, pfx, pfy, p2f, out, d2);
}

/* E105: nearest constellation point per cell, for the common-phase-error
 * correction (m9_fast._cpe_fast).
 *
 * Same disease the demapper had: the numpy path builds an (n, npts) float32
 * score array, scales it, adds the point energies and argmins along it --
 * three streaming passes over ~13 MB, three times per Frame.  Here each
 * cell's scores live in registers and only the winning INDEX is written.
 * The gather (alphabet[idx]) stays in numpy: it is one indexed copy, and
 * keeping it there keeps this kernel to one job.
 *
 * pfx/pfy arrive ALREADY scaled by -2, as in the demapper, so the score is
 * |p|^2 - 2 Re(z conj p) -- a monotone function of |z-p|^2, which is why its
 * argmin is the nearest point.
 *
 * TIE RULE: strict `<` keeps the FIRST minimum, exactly what np.argmin
 * returns.  Ties are not hypothetical here -- the reference's own comment
 * notes that "a flip needs two points equidistant to within a rounding" --
 * so the rule has to MATCH, not merely be reasonable.
 */
#if defined(__x86_64__) || defined(__i386__)
/* E106: eight points per vector.  Lane l tracks the first strict minimum
 * among the points j = l mod 8 in ascending j (the vectors ascend); across
 * lanes, the smallest INDEX among the lanes holding the global minimum is
 * the first minimum in j order -- np.argmin's rule, exactly. */
__attribute__((target("avx2")))
static int32_t cpe_nearest_avx2(const float *zr, const float *zi, int32_t n,
                                const float *pfx, const float *pfy,
                                const float *p2f, int32_t npts, int32_t *idx)
{
    const int32_t nv = npts / 8;
    const __m256i lane = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    for (int32_t c = 0; c < n; ++c) {
        const __m256 xv = _mm256_set1_ps(zr[c]), yv = _mm256_set1_ps(zi[c]);
        __m256 best = _mm256_set1_ps(3.402823466e+38f);
        __m256i arg = _mm256_setzero_si256();
        for (int32_t v = 0; v < nv; ++v) {
            const int32_t j0 = v * 8;
            const __m256 sc = _mm256_add_ps(
                _mm256_add_ps(_mm256_mul_ps(xv, _mm256_loadu_ps(pfx + j0)),
                              _mm256_mul_ps(yv, _mm256_loadu_ps(pfy + j0))),
                _mm256_loadu_ps(p2f + j0));
            const __m256 lt = _mm256_cmp_ps(sc, best, _CMP_LT_OQ);
            best = _mm256_blendv_ps(best, sc, lt);
            arg = _mm256_blendv_epi8(arg, _mm256_add_epi32(_mm256_set1_epi32(j0), lane),
                                     _mm256_castps_si256(lt));
        }
        float b[8]; int32_t a[8];
        _mm256_storeu_ps(b, best);
        _mm256_storeu_si256((__m256i *)a, arg);
        float gmin = b[0];
        for (int l = 1; l < 8; ++l) if (b[l] < gmin) gmin = b[l];
        int32_t win = 0x7fffffff;
        for (int l = 0; l < 8; ++l) if (b[l] == gmin && a[l] < win) win = a[l];
        idx[c] = win;
    }
    return 0;
}
#endif

EXPORT int32_t cpe_nearest_f32(const float *zr, const float *zi, int32_t n,
                               const float *pfx, const float *pfy,
                               const float *p2f, int32_t npts,
                               int32_t *idx)
{
    if (npts <= 0 || npts > DK_MAX_PTS || n < 0) return 1;
#if defined(__x86_64__) || defined(__i386__)
    if (npts >= 8 && (npts % 8) == 0 && __builtin_cpu_supports("avx2"))
        return cpe_nearest_avx2(zr, zi, n, pfx, pfy, p2f, npts, idx);
#endif
    for (int32_t c = 0; c < n; ++c) {
        const float x = zr[c], y = zi[c];
        float best = x * pfx[0] + y * pfy[0] + p2f[0];
        int32_t arg = 0;
        for (int32_t j = 1; j < npts; ++j) {
            const float sc = x * pfx[j] + y * pfy[j] + p2f[j];
            if (sc < best) { best = sc; arg = j; }
        }
        idx[c] = arg;
    }
    return 0;
}
