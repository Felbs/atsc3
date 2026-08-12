/* ldpc_kernel.c -- E54: compiled normalized min-sum, float32, per-block.
 *
 * WHAT THIS IS
 * ------------
 * The same algorithm as m3_ldpc._fast_loop (the float32 batched min-sum
 * with the syndrome FUSED before the message update), written as plain C
 * so the ~15 NumPy kernel dispatches per iteration become one loop nest.
 * E53 measured the numpy floor (~130 ms/frame @6T) to be DISPATCH-bound,
 * not math- or bandwidth-bound; this kernel is the lever it named.
 *
 * SEMANTICS (mirrored line for line from _fast_loop; the gate is the
 * referee, not this comment):
 *   - iteration it = 0..iters: FIRST check the syndrome of the current
 *     hard decisions (sign pattern of tot).  A block whose channel hard
 *     decisions already satisfy every check retires with its = 0 and NO
 *     BP work (PLP16 on clean air, measured).
 *   - message update: v = tot[var] - m_prev per edge; per check the
 *     normalized min-sum rule with the FIRST-minimum tie rule (only the
 *     first edge attaining min1 takes min2 -- numpy argmin semantics).
 *   - tot rebuild: per-variable segment sum of messages in check-major
 *     edge order, then tot = llr + s -- the np.add.reduceat association
 *     of the float32 path, not np.add.at's.
 *   - a row that never converges returns its final bits and its last
 *     unsatisfied-check count, its = iters, conv = 0.
 *
 * DELIBERATELY NOT HERE: threads (the caller's ThreadPoolExecutor
 * sub-batches blocks and ctypes releases the GIL -- the existing decode
 * topology parallelizes this kernel with no OpenMP dependency), python.h
 * (pure C ABI via ctypes: the same binary serves any CPython version --
 * though the fleet law stands: never copy the binary, rebuild per box),
 * and -ffast-math anywhere (float32 ops in source order; the build pins
 * -ffp-contract=off so no FMA contraction changes a rounding).
 *
 * Build: python lab/build_ldpc_kernel.py   (MSVC or gcc on Windows,
 *        cc/gcc on Linux; see PORT.md).
 */

#include <stddef.h>   /* size_t -- MSVC pulls it in transitively, gcc 15 does not */
#include <stdint.h>
#include <math.h>

#ifdef _WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif

#define MAX_ROW_DEG 128         /* loader verifies every check fits; 128
                                 * covers every A/322 table incl. 64800-bit
                                 * 7/15 (row degree 86); 16200 11/15 is 16 */

EXPORT int32_t ldpc_kernel_abi(void) { return 1; }

/* Decode nblk independent codewords, serially (caller parallelizes).
 *
 * llr    : nblk*n float32, C-contiguous
 * cptr   : nc+1   edge offsets per check
 * cidx   : nedges variable index of each edge (check-major order)
 * vptr   : n+1    edge-id offsets per variable
 * vedge  : nedges edge ids grouped by variable, check-major inside
 * work   : n + nedges float32 scratch (one per concurrent call)
 * bits   : nblk*n uint8 out    conv/its/bad : nblk out
 * returns 0, or -1 if a check degree exceeds MAX_ROW_DEG (nothing written)
 */
EXPORT int32_t ldpc_minsum_f32(
    const float *llr, int32_t nblk, int32_t n, int32_t nc, int32_t nedges,
    const int32_t *cptr, const int32_t *cidx,
    const int32_t *vptr, const int32_t *vedge,
    int32_t iters, float alpha,
    uint8_t *bits, uint8_t *conv, int32_t *its, int32_t *bad,
    float *work)
{
    float *tot = work;            /* n     */
    float *m   = work + n;        /* nedges */
    float v[MAX_ROW_DEG];
    int32_t c;

    for (c = 0; c < nc; c++)
        if (cptr[c + 1] - cptr[c] > MAX_ROW_DEG)
            return -1;

    for (int32_t b = 0; b < nblk; b++) {
        const float *lam = llr + (size_t)b * n;
        uint8_t *bout = bits + (size_t)b * n;
        int32_t i, e, it;

        for (i = 0; i < n; i++) tot[i] = lam[i];
        for (e = 0; e < nedges; e++) m[e] = 0.0f;
        conv[b] = 0;
        its[b] = iters;
        bad[b] = -1;

        for (it = 0; it <= iters; it++) {
            /* fused syndrome: sign pattern of tot at the check gathers */
            int32_t badct = 0;
            for (c = 0; c < nc; c++) {
                int par = 0;
                for (e = cptr[c]; e < cptr[c + 1]; e++)
                    par ^= (tot[cidx[e]] < 0.0f);
                badct += par;
            }
            if (badct == 0) {
                for (i = 0; i < n; i++) bout[i] = tot[i] < 0.0f;
                conv[b] = 1;
                its[b] = it;
                bad[b] = 0;
                break;
            }
            if (it == iters) {
                for (i = 0; i < n; i++) bout[i] = tot[i] < 0.0f;
                bad[b] = badct;           /* conv stays 0, its stays iters */
                break;
            }

            /* normalized min-sum update, first-minimum tie rule */
            for (c = 0; c < nc; c++) {
                int32_t e0 = cptr[c], e1 = cptr[c + 1];
                float min1 = INFINITY, min2 = INFINITY;
                int32_t first = e0;
                int par = 0;
                for (e = e0; e < e1; e++) {
                    float vv = tot[cidx[e]] - m[e];
                    v[e - e0] = vv;
                    par ^= (vv < 0.0f);
                    float mag = fabsf(vv);
                    if (mag < min1) { min2 = min1; min1 = mag; first = e; }
                    else if (mag < min2) { min2 = mag; }
                }
                for (e = e0; e < e1; e++) {
                    float val = ((e == first) ? min2 : min1) * alpha;
                    if (par ^ (v[e - e0] < 0.0f)) val = -val;
                    m[e] = val;
                }
            }

            /* tot = llr + per-variable segment sum (reduceat association) */
            for (i = 0; i < n; i++) {
                float s = 0.0f;
                for (int32_t k = vptr[i]; k < vptr[i + 1]; k++)
                    s += m[vedge[k]];
                tot[i] = lam[i] + s;
            }
        }
    }
    return 0;
}
