/* ldpc_layered_kernel.c -- E107 (2026-11-07): layered normalised min-sum,
 * eight codewords per AVX2 vector.
 *
 * WHAT THIS IS
 * ------------
 * The same check rule as ldpc_kernel.c (normalised min-sum, alpha, the
 * first-minimum tie rule, the syndrome checked BEFORE the first update so a
 * clean block retires with its = 0) on two changes of SCHEDULE, neither of
 * which is an approximation of the rule:
 *
 *   1. LAYERED.  The flooding kernel updates every check from the previous
 *      iteration's totals and only then rebuilds every total.  Here each
 *      check's new messages go into the totals AT ONCE, so the next check in
 *      the sweep already sees them.  Same fixed points; about half the
 *      iterations to reach them (MEASURED on this air, see the gate).  The
 *      order of the sweep is the table's own check order.
 *   2. BATCHED.  Eight independent codewords ride the eight lanes of one
 *      vector, so every gather, subtract, compare and select serves eight
 *      blocks.  Lanes never interact; an unused lane holds the all-zero word
 *      and is satisfied from the start.
 *
 * The idea of running the layered schedule on several blocks per vector we
 * learned from OnAir (adswill/OnAir, GPL-3.0, atsc3_ldpc.cpp, which does it
 * on sixteen 8-bit lanes).  This is written from the rule, not from their
 * code; it keeps our float32 LLRs, our alpha ladder and our tie rule, so
 * that a block both kernels decode yields the SAME codeword, and the gate
 * (gate_e107_ldpc.py) holds it to E54's boundary: the converged set may only
 * grow, and the Baseband bytes of every converged block are identical.
 *
 * DELIBERATELY NOT HERE: threads (the caller hands each group of eight to a
 * thread; ctypes drops the GIL), Python.h, -ffast-math (the build pins
 * -ffp-contract=off; every product and sum here is an explicit intrinsic in
 * the scalar kernel's association).  x86-64 only: on any other machine the
 * loader keeps the flooding kernel.
 *
 * Build: python lab/build_ldpc_layered_kernel.py
 */
#include <stddef.h>
#include <stdint.h>
#include <string.h>
#ifdef _WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif
#define LANES 8
#define MAX_ROW_DEG 128

EXPORT int32_t ldpc_layered_abi(void) { return 6; }   /* 4: + pack_f32 (E118); 5: + q8 (E130); 6: + q16 (E131) */

#if defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>

/* E110 -- the two-minimum correction of the check rule.
 *
 * Min-sum sends every edge of a check the smallest OTHER input magnitude; the
 * exact (belief-propagation) message is the box-plus of all the other inputs,
 * which min-sum approximates from above, and `alpha` scales it down to
 * compensate on average.  A better per-edge estimate costs one more number
 * the check already has: the edges that do not hold the smallest input see
 * the box-plus of the two smallest,
 *     min1 (+) min2 = min1 + ln(1 + e^-(min1+min2)) - ln(1 + e^-(min2-min1)),
 * which is closer to what belief propagation would send (Chen & Fossorier's
 * two-minimum correction; the textbook result).  It needs the LLRs in their
 * true scale -- ours are max-log LLRs over the demapper's own sigma^2, which
 * is that scale -- and we learned to try it here from OnAir (minsum_corr.h),
 * whose notes put it at 0.1-0.3 dB on the DVB-S2 codes.  Written from the
 * formula, not from their code.  ln(1 + e^-t) for t >= 0 is convex and sits
 * at or above every one of its tangents, so the largest of a few tangents is
 * a lower bound within ~0.02 of it and vectorises as max of affine functions:
 * tangents at t = 0.25, 1.0, 2.0 and 3.5, coefficients from the derivative
 * -1/(1 + e^t).  Off unless `corr` is set; gate_e107_ldpc.py --box is the referee. */
__attribute__((target("avx2")))
static inline __m256 softplus_neg(__m256 t)
{
    /* ln(1+e^-t) ~ max(a_i - b_i t) over tangents at t0 = 0.25, 1, 2, 3.5 */
    const __m256 a0 = _mm256_set1_ps(0.67415f), b0 = _mm256_set1_ps(0.43782f);
    const __m256 a1 = _mm256_set1_ps(0.58220f), b1 = _mm256_set1_ps(0.26894f);
    const __m256 a2 = _mm256_set1_ps(0.36534f), b2 = _mm256_set1_ps(0.11920f);
    const __m256 a3 = _mm256_set1_ps(0.13262f), b3 = _mm256_set1_ps(0.02931f);
    __m256 v = _mm256_max_ps(_mm256_sub_ps(a0, _mm256_mul_ps(b0, t)),
                             _mm256_sub_ps(a1, _mm256_mul_ps(b1, t)));
    v = _mm256_max_ps(v, _mm256_sub_ps(a2, _mm256_mul_ps(b2, t)));
    v = _mm256_max_ps(v, _mm256_sub_ps(a3, _mm256_mul_ps(b3, t)));
    return _mm256_max_ps(v, _mm256_setzero_ps());
}

__attribute__((target("avx2")))
static int32_t layered_avx2(const float *llr, int32_t nblk, int32_t n,
                            int32_t nc, int32_t nedges,
                            const int32_t *cptr, const int32_t *cidx,
                            int32_t iters, float alpha, int32_t corr,
                            int32_t abort_k, int32_t abort_t, int32_t abort_stall,
                            uint8_t *bits, uint8_t *conv, int32_t *its,
                            int32_t *bad, float *work)
{
    float *tot = work;                      /* n      x LANES */
    float *msg = work + (size_t)n * LANES;  /* nedges x LANES */
    const __m256 signmask = _mm256_set1_ps(-0.0f);
    const __m256 inf = _mm256_set1_ps(3.402823466e+38f);
    const __m256 alphav = _mm256_set1_ps(alpha);
    const __m256i one = _mm256_set1_epi32(1);
    __m256 v[MAX_ROW_DEG];
    int32_t c, e, i, b;
    for (c = 0; c < nc; c++)
        if (cptr[c + 1] - cptr[c] > MAX_ROW_DEG) return -1;
    if (nblk <= 0 || nblk > LANES) return -2;
    /* lanes: the blocks, then all-zero words */
    for (i = 0; i < n; i++)
        for (b = 0; b < LANES; b++)
            tot[(size_t)i * LANES + b] = b < nblk ? llr[(size_t)b * n + i] : 0.0f;
    memset(msg, 0, (size_t)nedges * LANES * sizeof(float));
    uint8_t done[LANES];
    for (b = 0; b < LANES; b++) { done[b] = b >= nblk; }
    int32_t last_bc[LANES], stall[LANES];
    for (b = 0; b < LANES; b++) { last_bc[b] = 0x7fffffff; stall[b] = 0; }
    for (b = 0; b < nblk; b++) { conv[b] = 0; its[b] = iters; bad[b] = -1; }
    int32_t ndone = LANES - nblk;
    for (int32_t it = 0; it <= iters; it++) {
        /* syndrome of the current hard decisions, per lane */
        __m256i badcnt = _mm256_setzero_si256();
        for (c = 0; c < nc; c++) {
            __m256i par = _mm256_setzero_si256();
            for (e = cptr[c]; e < cptr[c + 1]; e++) {
                const __m256 t = _mm256_loadu_ps(tot + (size_t)cidx[e] * LANES);
                par = _mm256_xor_si256(par, _mm256_srli_epi32(_mm256_castps_si256(t), 31));
            }
            badcnt = _mm256_add_epi32(badcnt, par);
        }
        int32_t bc[LANES];
        _mm256_storeu_si256((__m256i *)bc, badcnt);
        for (b = 0; b < nblk; b++) {
            if (done[b]) continue;
            /* E112: a block that is not improving, or is still far from the
             * code after enough sweeps, is given up: it keeps its hard
             * decisions, conv = 0, its = the sweep it stopped at, bad = its
             * syndrome weight (a failed block's meaning, as before).  With
             * abort_k = 0 nothing here fires and the kernel is ABI 2's. */
            int32_t give_up = 0;
            if (abort_k > 0 && it >= abort_k && bc[b] != 0) {
                if (abort_t > 0 && bc[b] >= abort_t) give_up = 1;
                if (abort_stall > 0 && stall[b] >= abort_stall) give_up = 1;
            }
            if (bc[b] >= last_bc[b]) stall[b]++; else stall[b] = 0;
            last_bc[b] = bc[b];
            if (bc[b] == 0 || it == iters || give_up) {
                uint8_t *bout = bits + (size_t)b * n;
                for (i = 0; i < n; i++) bout[i] = tot[(size_t)i * LANES + b] < 0.0f;
                if (bc[b] == 0) { conv[b] = 1; its[b] = it; bad[b] = 0; }
                else            { its[b] = it; bad[b] = bc[b]; }
                done[b] = 1;
                ndone++;
            }
        }
        if (ndone == LANES || it == iters) break;
        /* one layered sweep over the checks */
        for (c = 0; c < nc; c++) {
            const int32_t e0 = cptr[c], e1 = cptr[c + 1];
            __m256 min1 = inf, min2 = inf;
            __m256i first = _mm256_set1_epi32(-1), par = _mm256_setzero_si256();
            for (e = e0; e < e1; e++) {
                const __m256 vv = _mm256_sub_ps(
                    _mm256_loadu_ps(tot + (size_t)cidx[e] * LANES),
                    _mm256_loadu_ps(msg + (size_t)e * LANES));
                v[e - e0] = vv;
                par = _mm256_xor_si256(par, _mm256_srli_epi32(_mm256_castps_si256(vv), 31));
                const __m256 mag = _mm256_andnot_ps(signmask, vv);
                const __m256 lt1 = _mm256_cmp_ps(mag, min1, _CMP_LT_OQ);
                const __m256 lt2 = _mm256_cmp_ps(mag, min2, _CMP_LT_OQ);
                /* else-if: a new min1 pushes the old one down; otherwise a
                 * value under min2 replaces min2 -- the scalar kernel's rule */
                min2 = _mm256_blendv_ps(_mm256_blendv_ps(min2, mag, lt2), min1, lt1);
                min1 = _mm256_blendv_ps(min1, mag, lt1);
                first = _mm256_blendv_epi8(first, _mm256_set1_epi32(e), _mm256_castps_si256(lt1));
            }
            /* E110: the magnitude every edge but the first-minimum one receives */
            __m256 others = min1;
            if (corr) {
                const __m256 c = _mm256_sub_ps(
                    _mm256_add_ps(min1, softplus_neg(_mm256_add_ps(min1, min2))),
                    softplus_neg(_mm256_sub_ps(min2, min1)));
                others = _mm256_max_ps(c, _mm256_setzero_ps());
            }
            for (e = e0; e < e1; e++) {
                const __m256 vv = v[e - e0];
                const __m256i isfirst = _mm256_cmpeq_epi32(first, _mm256_set1_epi32(e));
                __m256 val = _mm256_mul_ps(
                    _mm256_blendv_ps(others, min2, _mm256_castsi256_ps(isfirst)), alphav);
                /* sign: the check's parity without this edge */
                const __m256i sgn = _mm256_xor_si256(par, _mm256_srli_epi32(_mm256_castps_si256(vv), 31));
                val = _mm256_xor_ps(val, _mm256_castsi256_ps(_mm256_slli_epi32(_mm256_and_si256(sgn, one), 31)));
                _mm256_storeu_ps(msg + (size_t)e * LANES, val);
                _mm256_storeu_ps(tot + (size_t)cidx[e] * LANES, _mm256_add_ps(vv, val));
            }
        }
    }
    return 0;
}

/* E118 -- lane compaction.
 *
 * MEASURED FIRST (RF25 banked on the 32-core box, 378 batches, 11,633 blocks,
 * ATSC3_LDPC_DUMP): a vector runs until its SLOWEST lane is done, and on a
 * carrier where most blocks converge the slow lane is a near miss that runs
 * all `iters` sweeps while seven lanes that retired at sweep 8-16 ride along
 * doing nothing -- 37,323 vector-sweeps where the lanes' own sweep counts sum
 * to 27,674 lane-vectors: 26% of the kernel's time, exactly.
 *
 * So: the caller hands this entry point ANY number of blocks.  Eight start;
 * whenever a lane's block retires (converged, aborted or at the cap) the next
 * block is loaded into that lane -- its column of `tot` set to its LLRs, its
 * column of `msg` zeroed -- with its own sweep counter starting at 0.  Lanes
 * never interact, so a block sees exactly the sequence of operations the
 * eight-at-a-time kernel would have given it: the syndrome checked BEFORE its
 * first sweep (a clean block retires with its = 0, at load time, and the lane
 * takes the next one at once), the same sweep arithmetic, the same abort
 * rule on its own counter.  Bit-identical per block by construction; the
 * referee is lab/gate_e118_ldpc_pack.py (every bit, its and bad against the
 * old entry point on banked LLRs).  The old entry point stays for that
 * comparison and for any caller that wants it. */
__attribute__((target("avx2")))
static int32_t layered_pack_avx2(const float *llr, int32_t nblk, int32_t n,
                                 int32_t nc, int32_t nedges,
                                 const int32_t *cptr, const int32_t *cidx,
                                 int32_t iters, float alpha, int32_t corr,
                                 int32_t abort_k, int32_t abort_t, int32_t abort_stall,
                                 uint8_t *bits, uint8_t *conv, int32_t *its,
                                 int32_t *bad, float *work)
{
    float *tot = work;                      /* n      x LANES */
    float *msg = work + (size_t)n * LANES;  /* nedges x LANES */
    const __m256 signmask = _mm256_set1_ps(-0.0f);
    const __m256 inf = _mm256_set1_ps(3.402823466e+38f);
    const __m256 alphav = _mm256_set1_ps(alpha);
    const __m256i one = _mm256_set1_epi32(1);
    __m256 v[MAX_ROW_DEG];
    int32_t c, e, i, b;
    for (c = 0; c < nc; c++)
        if (cptr[c + 1] - cptr[c] > MAX_ROW_DEG) return -1;
    if (nblk <= 0) return -2;
    for (b = 0; b < nblk; b++) { conv[b] = 0; its[b] = iters; bad[b] = -1; }
    int32_t lane_blk[LANES], lane_it[LANES], last_bc[LANES], stall[LANES];
    int32_t next = 0, nactive = 0;
    for (b = 0; b < LANES; b++) { lane_blk[b] = -1; lane_it[b] = 0; last_bc[b] = 0x7fffffff; stall[b] = 0; }
    /* Columns of a lane that holds no block are never read for a verdict and
     * never influence another lane, so they are left as they are.  Only a
     * lane about to run needs its column written: tot = the block's LLRs,
     * msg = 0 -- done below for every refilled lane in ONE pass per array
     * (contiguous writes, a masked vector store), with the sweep-0 syndrome
     * taken by the same vector check the sweeps use.  A clean block retires
     * there with its = 0 and its lane is refilled again at once. */
    uint32_t pend = (1u << LANES) - 1;      /* lanes waiting for a block */
    int32_t bc[LANES];

    for (;;) {
        /* ---- refill: every lane in `pend`, as one event ---- */
        while (pend) {
            const float *src[LANES];
            uint32_t fresh = 0;
            for (b = 0; b < LANES; b++) {
                src[b] = NULL;
                if (!((pend >> b) & 1u)) continue;
                if (next < nblk) {
                    lane_blk[b] = next;
                    src[b] = llr + (size_t)next * n;
                    next++;
                    fresh |= 1u << b;
                } else {
                    lane_blk[b] = -1;          /* nothing left: lane parked */
                    pend &= ~(1u << b);
                }
            }
            if (!fresh) break;
            /* tot columns, i-outer so the writes are contiguous */
            for (i = 0; i < n; i++) {
                float *row = tot + (size_t)i * LANES;
                for (b = 0; b < LANES; b++)
                    if ((fresh >> b) & 1u) row[b] = src[b][i];
            }
            /* msg columns of the fresh lanes -> 0, one masked vector pass */
            {
                int32_t keep[LANES];
                for (b = 0; b < LANES; b++) keep[b] = ((fresh >> b) & 1u) ? 0 : -1;
                const __m256 km = _mm256_castsi256_ps(_mm256_loadu_si256((const __m256i *)keep));
                for (e = 0; e < nedges; e++) {
                    float *m = msg + (size_t)e * LANES;
                    _mm256_storeu_ps(m, _mm256_and_ps(_mm256_loadu_ps(m), km));
                }
            }
            /* the sweep-0 check of the fresh lanes: the vector syndrome */
            {
                __m256i badcnt = _mm256_setzero_si256();
                for (c = 0; c < nc; c++) {
                    __m256i par = _mm256_setzero_si256();
                    for (e = cptr[c]; e < cptr[c + 1]; e++) {
                        const __m256 t = _mm256_loadu_ps(tot + (size_t)cidx[e] * LANES);
                        par = _mm256_xor_si256(par, _mm256_srli_epi32(_mm256_castps_si256(t), 31));
                    }
                    badcnt = _mm256_add_epi32(badcnt, par);
                }
                _mm256_storeu_si256((__m256i *)bc, badcnt);
            }
            for (b = 0; b < LANES; b++) {
                if (!((fresh >> b) & 1u)) continue;
                const int32_t blk = lane_blk[b];
                if (bc[b] == 0) {      /* clean before any sweep: its = 0, as before */
                    uint8_t *bout = bits + (size_t)blk * n;
                    for (i = 0; i < n; i++) bout[i] = tot[(size_t)i * LANES + b] < 0.0f;
                    conv[blk] = 1; its[blk] = 0; bad[blk] = 0;
                    lane_blk[b] = -1;            /* stays in pend: refill again */
                } else {
                    lane_it[b] = 0; stall[b] = 0; last_bc[b] = bc[b];   /* sweep-0 bookkeeping */
                    nactive++;
                    pend &= ~(1u << b);
                }
            }
        }
        if (nactive == 0) break;

        /* ---- one layered sweep over the checks: the eight-lane kernel's, verbatim ---- */
        for (c = 0; c < nc; c++) {
            const int32_t e0 = cptr[c], e1 = cptr[c + 1];
            __m256 min1 = inf, min2 = inf;
            __m256i first = _mm256_set1_epi32(-1), par = _mm256_setzero_si256();
            for (e = e0; e < e1; e++) {
                const __m256 vv = _mm256_sub_ps(
                    _mm256_loadu_ps(tot + (size_t)cidx[e] * LANES),
                    _mm256_loadu_ps(msg + (size_t)e * LANES));
                v[e - e0] = vv;
                par = _mm256_xor_si256(par, _mm256_srli_epi32(_mm256_castps_si256(vv), 31));
                const __m256 mag = _mm256_andnot_ps(signmask, vv);
                const __m256 lt1 = _mm256_cmp_ps(mag, min1, _CMP_LT_OQ);
                const __m256 lt2 = _mm256_cmp_ps(mag, min2, _CMP_LT_OQ);
                min2 = _mm256_blendv_ps(_mm256_blendv_ps(min2, mag, lt2), min1, lt1);
                min1 = _mm256_blendv_ps(min1, mag, lt1);
                first = _mm256_blendv_epi8(first, _mm256_set1_epi32(e), _mm256_castps_si256(lt1));
            }
            __m256 others = min1;
            if (corr) {
                const __m256 cc = _mm256_sub_ps(
                    _mm256_add_ps(min1, softplus_neg(_mm256_add_ps(min1, min2))),
                    softplus_neg(_mm256_sub_ps(min2, min1)));
                others = _mm256_max_ps(cc, _mm256_setzero_ps());
            }
            for (e = e0; e < e1; e++) {
                const __m256 vv = v[e - e0];
                const __m256i isfirst = _mm256_cmpeq_epi32(first, _mm256_set1_epi32(e));
                __m256 val = _mm256_mul_ps(
                    _mm256_blendv_ps(others, min2, _mm256_castsi256_ps(isfirst)), alphav);
                const __m256i sgn = _mm256_xor_si256(par, _mm256_srli_epi32(_mm256_castps_si256(vv), 31));
                val = _mm256_xor_ps(val, _mm256_castsi256_ps(_mm256_slli_epi32(_mm256_and_si256(sgn, one), 31)));
                _mm256_storeu_ps(msg + (size_t)e * LANES, val);
                _mm256_storeu_ps(tot + (size_t)cidx[e] * LANES, _mm256_add_ps(vv, val));
            }
        }
        for (b = 0; b < LANES; b++) if (lane_blk[b] >= 0) lane_it[b]++;

        /* ---- syndrome of the current hard decisions, per lane ---- */
        {
            __m256i badcnt = _mm256_setzero_si256();
            for (c = 0; c < nc; c++) {
                __m256i par = _mm256_setzero_si256();
                for (e = cptr[c]; e < cptr[c + 1]; e++) {
                    const __m256 t = _mm256_loadu_ps(tot + (size_t)cidx[e] * LANES);
                    par = _mm256_xor_si256(par, _mm256_srli_epi32(_mm256_castps_si256(t), 31));
                }
                badcnt = _mm256_add_epi32(badcnt, par);
            }
            _mm256_storeu_si256((__m256i *)bc, badcnt);
        }
        for (b = 0; b < LANES; b++) {
            const int32_t blk = lane_blk[b];
            if (blk < 0) continue;
            const int32_t it = lane_it[b];
            int32_t give_up = 0;
            if (abort_k > 0 && it >= abort_k && bc[b] != 0) {
                if (abort_t > 0 && bc[b] >= abort_t) give_up = 1;
                if (abort_stall > 0 && stall[b] >= abort_stall) give_up = 1;
            }
            if (bc[b] >= last_bc[b]) stall[b]++; else stall[b] = 0;
            last_bc[b] = bc[b];
            if (bc[b] == 0 || it == iters || give_up) {
                uint8_t *bout = bits + (size_t)blk * n;
                for (i = 0; i < n; i++) bout[i] = tot[(size_t)i * LANES + b] < 0.0f;
                if (bc[b] == 0) { conv[blk] = 1; its[blk] = it; bad[blk] = 0; }
                else            { its[blk] = it; bad[blk] = bc[b]; }
                lane_blk[b] = -1;
                nactive--;
                pend |= 1u << b;
            }
        }
    }
    return 0;
}
#endif

/* Decode up to LANES codewords at once.
 *
 * llr    : nblk*n float32, C-contiguous          nblk  : 1..LANES
 * corr   : 0 = normalised min-sum; 1 = with the two-minimum correction (E110)
 * abort_k/abort_t/abort_stall : E112 -- from sweep abort_k on, give a block up
 *          when its syndrome weight is >= abort_t, or has not improved for
 *          abort_stall consecutive sweeps; 0 disables each.  All 0 = ABI 2.
 * cptr   : nc+1 edge offsets per check           cidx  : nedges variable index per edge
 * work   : (n + nedges) * LANES float32 scratch (one per concurrent call)
 * bits   : nblk*n uint8 out    conv/its/bad : nblk out (the flooding kernel's meaning)
 * returns 0; -1 a check degree over MAX_ROW_DEG; -2 nblk out of range;
 *         -3 this machine has no AVX2 (nothing written: keep the flooding kernel)
 */
EXPORT int32_t ldpc_layered_f32(
    const float *llr, int32_t nblk, int32_t n, int32_t nc, int32_t nedges,
    const int32_t *cptr, const int32_t *cidx,
    int32_t iters, float alpha, int32_t corr,
    int32_t abort_k, int32_t abort_t, int32_t abort_stall,
    uint8_t *bits, uint8_t *conv, int32_t *its, int32_t *bad,
    float *work)
{
#if defined(__x86_64__) || defined(_M_X64)
    if (__builtin_cpu_supports("avx2"))
        return layered_avx2(llr, nblk, n, nc, nedges, cptr, cidx, iters, alpha, corr,
                            abort_k, abort_t, abort_stall,
                            bits, conv, its, bad, work);
#endif
    (void)llr; (void)nblk; (void)n; (void)nc; (void)nedges; (void)cptr; (void)cidx;
    (void)iters; (void)alpha; (void)corr; (void)abort_k; (void)abort_t; (void)abort_stall; (void)bits; (void)conv; (void)its; (void)bad; (void)work;
    return -3;
}

/* E118: the same call for ANY nblk >= 1 (one call per thread's share of the
 * batch), lanes refilled as blocks retire.  Same outputs, same meanings,
 * bit-identical per block to ldpc_layered_f32 (gate_e118_ldpc_pack.py).
 * work is still (n + nedges) * LANES float32. */
EXPORT int32_t ldpc_layered_pack_f32(
    const float *llr, int32_t nblk, int32_t n, int32_t nc, int32_t nedges,
    const int32_t *cptr, const int32_t *cidx,
    int32_t iters, float alpha, int32_t corr,
    int32_t abort_k, int32_t abort_t, int32_t abort_stall,
    uint8_t *bits, uint8_t *conv, int32_t *its, int32_t *bad,
    float *work)
{
#if defined(__x86_64__) || defined(_M_X64)
    if (__builtin_cpu_supports("avx2"))
        return layered_pack_avx2(llr, nblk, n, nc, nedges, cptr, cidx, iters, alpha, corr,
                                 abort_k, abort_t, abort_stall,
                                 bits, conv, its, bad, work);
#endif
    (void)llr; (void)nblk; (void)n; (void)nc; (void)nedges; (void)cptr; (void)cidx;
    (void)iters; (void)alpha; (void)corr; (void)abort_k; (void)abort_t; (void)abort_stall; (void)bits; (void)conv; (void)its; (void)bad; (void)work;
    return -3;
}

/* E130 -- the layered schedule on SIXTEEN int8 lanes.
 *
 * MEASURED FIRST (13.29, RF33 bank, 3970X): this kernel's float32 working set
 * is 1.4 MB per 64800-bit block (tot n x 8 x 4 B + msg nedges x 8 x 4 B,
 * shared by eight), and a worker's PLP 1 LDPC runs 62 ms alone, 124 ms with
 * four workers on the box, 226 with eight, 356 with sixteen -- four workers
 * pinned one per core complex 74 ms, the same four packed into one complex
 * 134 ms.  Cache capacity, not CPU.  Sixteen int8 lanes hold a block in
 * 0.34 MB.  The idea of eight-bit lanes we learned from OnAir
 * (adswill/OnAir, GPL-3.0, atsc3_ldpc.cpp); the quantiser, the ranges, the
 * tie rule and the schedule here are written from the min-sum rule.
 *
 * ARITHMETIC.  Posteriors |tot| <= Q8_POST (96), messages |msg| <= Q8_MSG
 * (31), so tot - msg fits int8 without saturating.  LLRs enter as
 * round(llr * scale) clamped to +-Q8_POST (the caller picks `scale` from the
 * batch's own LLR magnitudes; normalised min-sum is scale-free except through
 * the clamps and the unit step).  alpha is applied in Q8 (x alpha*256 >> 8).
 * No two-minimum correction (it needs the LLRs' true scale); a block this
 * kernel does not converge is handed back for the float32 kernel to try
 * (the caller's hybrid), so the converged set can only grow.  Same E112
 * abort rule, same first-minimum tie rule, same syndrome-before-first-sweep.
 *
 * work: (n + nedges) * 16 int8.  returns 0; -1 degree > MAX_ROW_DEG;
 * -2 nblk out of 1..16; -3 no AVX2. */
#define Q8_LANES 16
#define Q8_POST 96
#define Q8_MSG 31

#if defined(__x86_64__) || defined(_M_X64)
__attribute__((target("avx2")))
static inline __m128i q8_scale(__m128i mag, __m128i alpha16)
{
    /* unsigned 8-bit magnitudes x alpha (Q8) >> 8, back to 8 bits */
    const __m128i lo = _mm_cvtepu8_epi16(mag);
    const __m128i hi = _mm_cvtepu8_epi16(_mm_srli_si128(mag, 8));
    return _mm_packus_epi16(_mm_srli_epi16(_mm_mullo_epi16(lo, alpha16), 8),
                            _mm_srli_epi16(_mm_mullo_epi16(hi, alpha16), 8));
}

__attribute__((target("avx2")))
static int32_t layered_q8_x86(const float *llr, int32_t nblk, int32_t n,
                              int32_t nc, int32_t nedges,
                              const int32_t *cptr, const int32_t *cidx,
                              int32_t iters, float alpha, float scale,
                              int32_t abort_k, int32_t abort_t, int32_t abort_stall,
                              uint8_t *bits, uint8_t *conv, int32_t *its,
                              int32_t *bad, int8_t *work)
{
    int8_t *tot = work;                          /* n      x 16 */
    int8_t *msg = work + (size_t)n * Q8_LANES;   /* nedges x 16 */
    const __m128i zero = _mm_setzero_si128();
    const __m128i one8 = _mm_set1_epi8(1);
    const __m128i ff = _mm_set1_epi8((char)0xFF);          /* 255 as unsigned: "infinity" */
    const __m128i msgmax = _mm_set1_epi8(Q8_MSG);
    const __m128i pmax = _mm_set1_epi8(Q8_POST), pmin = _mm_set1_epi8(-Q8_POST);
    int32_t a16 = (int32_t)(alpha * 256.0f + 0.5f);
    if (a16 < 1) a16 = 1;
    if (a16 > 256) a16 = 256;
    const __m128i alpha16 = _mm_set1_epi16((short)a16);
    __m128i v[MAX_ROW_DEG];
    int32_t c, e, i, b;
    for (c = 0; c < nc; c++)
        if (cptr[c + 1] - cptr[c] > MAX_ROW_DEG) return -1;
    if (nblk <= 0 || nblk > Q8_LANES) return -2;
    /* quantise into the lanes; unused lanes hold the all-zero word */
    for (i = 0; i < n; i++) {
        int8_t *row = tot + (size_t)i * Q8_LANES;
        for (b = 0; b < Q8_LANES; b++) {
            if (b >= nblk) { row[b] = 0; continue; }
            float x = llr[(size_t)b * n + i] * scale;
            int32_t q = (int32_t)(x >= 0.0f ? x + 0.5f : x - 0.5f);
            if (q > Q8_POST) q = Q8_POST;
            if (q < -Q8_POST) q = -Q8_POST;
            row[b] = (int8_t)q;
        }
    }
    memset(msg, 0, (size_t)nedges * Q8_LANES);
    uint8_t done[Q8_LANES];
    for (b = 0; b < Q8_LANES; b++) done[b] = b >= nblk;
    int32_t last_bc[Q8_LANES], stall[Q8_LANES];
    for (b = 0; b < Q8_LANES; b++) { last_bc[b] = 0x7fffffff; stall[b] = 0; }
    for (b = 0; b < nblk; b++) { conv[b] = 0; its[b] = iters; bad[b] = -1; }
    int32_t ndone = Q8_LANES - nblk;
    for (int32_t it = 0; it <= iters; it++) {
        /* syndrome of the current hard decisions, per lane (saturating int16 counts) */
        __m128i bad_lo = zero, bad_hi = zero;
        for (c = 0; c < nc; c++) {
            __m128i par = zero;
            for (e = cptr[c]; e < cptr[c + 1]; e++)
                par = _mm_xor_si128(par, _mm_cmplt_epi8(
                    _mm_loadu_si128((const __m128i *)(tot + (size_t)cidx[e] * Q8_LANES)), zero));
            bad_lo = _mm_subs_epi16(bad_lo, _mm_cvtepi8_epi16(par));
            bad_hi = _mm_subs_epi16(bad_hi, _mm_cvtepi8_epi16(_mm_srli_si128(par, 8)));
        }
        int16_t bc16[Q8_LANES];
        _mm_storeu_si128((__m128i *)bc16, bad_lo);
        _mm_storeu_si128((__m128i *)(bc16 + 8), bad_hi);
        for (b = 0; b < nblk; b++) {
            if (done[b]) continue;
            const int32_t bcb = bc16[b];
            int32_t give_up = 0;
            if (abort_k > 0 && it >= abort_k && bcb != 0) {
                if (abort_t > 0 && bcb >= abort_t) give_up = 1;
                if (abort_stall > 0 && stall[b] >= abort_stall) give_up = 1;
            }
            if (bcb >= last_bc[b]) stall[b]++; else stall[b] = 0;
            last_bc[b] = bcb;
            if (bcb == 0 || it == iters || give_up) {
                uint8_t *bout = bits + (size_t)b * n;
                for (i = 0; i < n; i++) bout[i] = tot[(size_t)i * Q8_LANES + b] < 0;
                if (bcb == 0) { conv[b] = 1; its[b] = it; bad[b] = 0; }
                else          { its[b] = it; bad[b] = bcb; }
                done[b] = 1;
                ndone++;
            }
        }
        if (ndone == Q8_LANES || it == iters) break;
        /* one layered sweep over the checks */
        for (c = 0; c < nc; c++) {
            const int32_t e0 = cptr[c], e1 = cptr[c + 1];
            __m128i min1 = ff, min2 = ff, first = ff, par = zero;
            for (e = e0; e < e1; e++) {
                const __m128i vv = _mm_subs_epi8(
                    _mm_loadu_si128((const __m128i *)(tot + (size_t)cidx[e] * Q8_LANES)),
                    _mm_loadu_si128((const __m128i *)(msg + (size_t)e * Q8_LANES)));
                v[e - e0] = vv;
                par = _mm_xor_si128(par, _mm_cmplt_epi8(vv, zero));
                const __m128i mag = _mm_abs_epi8(vv);
                /* strictly-less, unsigned: max(mag, m) == m and mag != m */
                const __m128i lt1 = _mm_andnot_si128(_mm_cmpeq_epi8(mag, min1),
                                                     _mm_cmpeq_epi8(_mm_max_epu8(mag, min1), min1));
                const __m128i lt2 = _mm_andnot_si128(_mm_cmpeq_epi8(mag, min2),
                                                     _mm_cmpeq_epi8(_mm_max_epu8(mag, min2), min2));
                min2 = _mm_blendv_epi8(_mm_blendv_epi8(min2, mag, lt2), min1, lt1);
                min1 = _mm_blendv_epi8(min1, mag, lt1);
                first = _mm_blendv_epi8(first, _mm_set1_epi8((char)(e - e0)), lt1);
            }
            const __m128i m_others = _mm_min_epu8(q8_scale(min1, alpha16), msgmax);
            const __m128i m_first = _mm_min_epu8(q8_scale(min2, alpha16), msgmax);
            for (e = e0; e < e1; e++) {
                const __m128i vv = v[e - e0];
                const __m128i isfirst = _mm_cmpeq_epi8(first, _mm_set1_epi8((char)(e - e0)));
                const __m128i mag = _mm_blendv_epi8(m_others, m_first, isfirst);
                /* sign: the check's parity without this edge; -1 where negative, +1 elsewhere */
                const __m128i sgn = _mm_or_si128(_mm_xor_si128(par, _mm_cmplt_epi8(vv, zero)), one8);
                const __m128i val = _mm_sign_epi8(mag, sgn);
                _mm_storeu_si128((__m128i *)(msg + (size_t)e * Q8_LANES), val);
                const __m128i nt = _mm_min_epi8(_mm_max_epi8(_mm_adds_epi8(vv, val), pmin), pmax);
                _mm_storeu_si128((__m128i *)(tot + (size_t)cidx[e] * Q8_LANES), nt);
            }
        }
    }
    return 0;
}
#endif

EXPORT int32_t ldpc_layered_q8(
    const float *llr, int32_t nblk, int32_t n, int32_t nc, int32_t nedges,
    const int32_t *cptr, const int32_t *cidx,
    int32_t iters, float alpha, float scale,
    int32_t abort_k, int32_t abort_t, int32_t abort_stall,
    uint8_t *bits, uint8_t *conv, int32_t *its, int32_t *bad,
    int8_t *work)
{
#if defined(__x86_64__) || defined(_M_X64)
    if (__builtin_cpu_supports("avx2"))
        return layered_q8_x86(llr, nblk, n, nc, nedges, cptr, cidx, iters, alpha, scale,
                              abort_k, abort_t, abort_stall, bits, conv, its, bad, work);
#endif
    (void)llr; (void)nblk; (void)n; (void)nc; (void)nedges; (void)cptr; (void)cidx;
    (void)iters; (void)alpha; (void)scale; (void)abort_k; (void)abort_t; (void)abort_stall;
    (void)bits; (void)conv; (void)its; (void)bad; (void)work;
    return -3;
}

/* E131 -- the layered schedule on SIXTEEN int16 lanes (one AVX2 vector), with
 * the two-minimum correction.
 *
 * E130's int8 lanes (above) quartered the working set and lost 11 % of the
 * Blocks at the cliff: no correction, and a quantiser step of 1/8 of the mean
 * LLR.  Here the step is 1/1024 of the mean (the caller scales the batch so
 * its mean |LLR| sits at Q16_MEAN units), the correction is kept (computed per
 * check on the two minima in float32 -- 16 values per check, not per edge --
 * and rounded back), and the footprint is half the float32 kernel's: tot
 * n x 16 x 2 B, msg nedges x 16 x 2 B, 0.69 MB per 64800-bit block.
 * Posteriors clamp at +-Q16_POST and messages at +-Q16_MSG so tot - msg
 * cannot overflow; all adds and subtracts saturate.  Same first-minimum tie
 * rule, syndrome before the first sweep, E112 abort rule.  Referee:
 * lab/gate_e130_q8.py --q16 (no block lost vs float32 at any level, bytes
 * identical on shared blocks).
 *
 * work: (n + nedges) * 16 int16.  returns 0; -1 degree > MAX_ROW_DEG;
 * -2 nblk out of 1..16; -3 no AVX2. */
#define Q16_LANES 16
#define Q16_POST 24576
#define Q16_MSG 8191

#if defined(__x86_64__) || defined(_M_X64)
__attribute__((target("avx2")))
static inline __m256i q16_scale(__m256i mag, __m256i alpha_q15)
{
    /* unsigned 16-bit magnitudes (<= Q16_POST + Q16_MSG < 32768) x alpha (Q15): mulhrs = (a*b + 2^14) >> 15 */
    return _mm256_mulhrs_epi16(mag, alpha_q15);
}

/* ln(1 + e^-t), t >= 0, as the largest of four tangents (the f32 kernel's); on 8 floats */
__attribute__((target("avx2,fma")))
static inline __m256 softplus_neg8(__m256 t)
{
    __m256 r = _mm256_fnmadd_ps(_mm256_set1_ps(0.4378235f), t, _mm256_set1_ps(0.8259394f));   /* tangent at 0.25 */
    r = _mm256_max_ps(r, _mm256_fnmadd_ps(_mm256_set1_ps(0.2689414f), t, _mm256_set1_ps(0.5822404f)));   /* 1.0 */
    r = _mm256_max_ps(r, _mm256_fnmadd_ps(_mm256_set1_ps(0.1192029f), t, _mm256_set1_ps(0.3652897f)));   /* 2.0 */
    r = _mm256_max_ps(r, _mm256_fnmadd_ps(_mm256_set1_ps(0.0293122f), t, _mm256_set1_ps(0.1323831f)));   /* 3.5 */
    return _mm256_max_ps(r, _mm256_setzero_ps());
}

/* the corrected magnitude for every edge but the first-minimum one:
 *   max(0, min1 + f(min1 + min2) - f(min2 - min1)),  f = ln(1 + e^-t),
 * in LLR units: the lanes hold mins in Q16 units, `inv` = 1 / units-per-LLR */
__attribute__((target("avx2,fma")))
static inline __m256i q16_boxplus(__m256i min1, __m256i min2, float inv, float units)
{
    const __m256 invv = _mm256_set1_ps(inv), unitsv = _mm256_set1_ps(units);
    __m256i out_lo, out_hi;
    {
        const __m256 a = _mm256_mul_ps(_mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm256_castsi256_si128(min1))), invv);
        const __m256 b = _mm256_mul_ps(_mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm256_castsi256_si128(min2))), invv);
        __m256 c = _mm256_sub_ps(_mm256_add_ps(a, softplus_neg8(_mm256_add_ps(a, b))), softplus_neg8(_mm256_sub_ps(b, a)));
        c = _mm256_max_ps(c, _mm256_setzero_ps());
        out_lo = _mm256_cvtps_epi32(_mm256_mul_ps(c, unitsv));
    }
    {
        const __m256 a = _mm256_mul_ps(_mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm256_extracti128_si256(min1, 1))), invv);
        const __m256 b = _mm256_mul_ps(_mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm256_extracti128_si256(min2, 1))), invv);
        __m256 c = _mm256_sub_ps(_mm256_add_ps(a, softplus_neg8(_mm256_add_ps(a, b))), softplus_neg8(_mm256_sub_ps(b, a)));
        c = _mm256_max_ps(c, _mm256_setzero_ps());
        out_hi = _mm256_cvtps_epi32(_mm256_mul_ps(c, unitsv));
    }
    /* int32 -> uint16 saturating, lanes back in order (packus works per 128-bit half) */
    const __m256i packed = _mm256_packus_epi32(out_lo, out_hi);
    return _mm256_permute4x64_epi64(packed, 0xD8);
}

__attribute__((target("avx2,fma")))
static int32_t layered_q16_avx2(const float *llr, int32_t nblk, int32_t n,
                                int32_t nc, int32_t nedges,
                                const int32_t *cptr, const int32_t *cidx,
                                int32_t iters, float alpha, float scale, int32_t corr,
                                int32_t abort_k, int32_t abort_t, int32_t abort_stall,
                                uint8_t *bits, uint8_t *conv, int32_t *its,
                                int32_t *bad, int16_t *work)
{
    int16_t *tot = work;                           /* n      x 16 */
    int16_t *msg = work + (size_t)n * Q16_LANES;   /* nedges x 16 */
    const __m256i zero = _mm256_setzero_si256();
    const __m256i one16 = _mm256_set1_epi16(1);
    const __m256i big = _mm256_set1_epi16((short)0x7FFF);
    const __m256i msgmax = _mm256_set1_epi16(Q16_MSG);
    const __m256i pmax = _mm256_set1_epi16(Q16_POST), pmin = _mm256_set1_epi16(-Q16_POST);
    int32_t aq = (int32_t)(alpha * 32768.0f + 0.5f);
    if (aq < 1) aq = 1;
    if (aq > 32767) aq = 32767;
    const __m256i alpha_q15 = _mm256_set1_epi16((short)aq);
    const float inv = 1.0f / scale, units = scale;
    __m256i v[MAX_ROW_DEG];
    int32_t c, e, i, b;
    for (c = 0; c < nc; c++)
        if (cptr[c + 1] - cptr[c] > MAX_ROW_DEG) return -1;
    if (nblk <= 0 || nblk > Q16_LANES) return -2;
    for (i = 0; i < n; i++) {
        int16_t *row = tot + (size_t)i * Q16_LANES;
        for (b = 0; b < Q16_LANES; b++) {
            if (b >= nblk) { row[b] = 0; continue; }
            float x = llr[(size_t)b * n + i] * scale;
            int32_t q = (int32_t)(x >= 0.0f ? x + 0.5f : x - 0.5f);
            if (q > Q16_POST) q = Q16_POST;
            if (q < -Q16_POST) q = -Q16_POST;
            row[b] = (int16_t)q;
        }
    }
    memset(msg, 0, (size_t)nedges * Q16_LANES * sizeof(int16_t));
    uint8_t done[Q16_LANES];
    for (b = 0; b < Q16_LANES; b++) done[b] = b >= nblk;
    int32_t last_bc[Q16_LANES], stall[Q16_LANES];
    for (b = 0; b < Q16_LANES; b++) { last_bc[b] = 0x7fffffff; stall[b] = 0; }
    for (b = 0; b < nblk; b++) { conv[b] = 0; its[b] = iters; bad[b] = -1; }
    int32_t ndone = Q16_LANES - nblk;
    for (int32_t it = 0; it <= iters; it++) {
        /* syndrome per lane: saturating int16 counts (a count past 32767 is hopeless anyway) */
        __m256i badcnt = zero;
        for (c = 0; c < nc; c++) {
            __m256i par = zero;
            for (e = cptr[c]; e < cptr[c + 1]; e++)
                par = _mm256_xor_si256(par, _mm256_cmpgt_epi16(zero,
                    _mm256_loadu_si256((const __m256i *)(tot + (size_t)cidx[e] * Q16_LANES))));
            badcnt = _mm256_subs_epi16(badcnt, par);
        }
        int16_t bc16[Q16_LANES];
        _mm256_storeu_si256((__m256i *)bc16, badcnt);
        for (b = 0; b < nblk; b++) {
            if (done[b]) continue;
            const int32_t bcb = bc16[b];
            int32_t give_up = 0;
            if (abort_k > 0 && it >= abort_k && bcb != 0) {
                if (abort_t > 0 && bcb >= abort_t) give_up = 1;
                if (abort_stall > 0 && stall[b] >= abort_stall) give_up = 1;
            }
            if (bcb >= last_bc[b]) stall[b]++; else stall[b] = 0;
            last_bc[b] = bcb;
            if (bcb == 0 || it == iters || give_up) {
                uint8_t *bout = bits + (size_t)b * n;
                for (i = 0; i < n; i++) bout[i] = tot[(size_t)i * Q16_LANES + b] < 0;
                if (bcb == 0) { conv[b] = 1; its[b] = it; bad[b] = 0; }
                else          { its[b] = it; bad[b] = bcb; }
                done[b] = 1;
                ndone++;
            }
        }
        if (ndone == Q16_LANES || it == iters) break;
        for (c = 0; c < nc; c++) {
            const int32_t e0 = cptr[c], e1 = cptr[c + 1];
            __m256i min1 = big, min2 = big, first = _mm256_set1_epi16(-1), par = zero;
            for (e = e0; e < e1; e++) {
                const __m256i vv = _mm256_subs_epi16(
                    _mm256_loadu_si256((const __m256i *)(tot + (size_t)cidx[e] * Q16_LANES)),
                    _mm256_loadu_si256((const __m256i *)(msg + (size_t)e * Q16_LANES)));
                v[e - e0] = vv;
                par = _mm256_xor_si256(par, _mm256_cmpgt_epi16(zero, vv));
                const __m256i mag = _mm256_abs_epi16(vv);          /* <= 32767: |vv| <= Q16_POST + Q16_MSG */
                const __m256i lt1 = _mm256_cmpgt_epi16(min1, mag);  /* strictly less, signed is fine: all < 32768 */
                const __m256i lt2 = _mm256_cmpgt_epi16(min2, mag);
                min2 = _mm256_blendv_epi8(_mm256_blendv_epi8(min2, mag, lt2), min1, lt1);
                min1 = _mm256_blendv_epi8(min1, mag, lt1);
                first = _mm256_blendv_epi8(first, _mm256_set1_epi16((short)(e - e0)), lt1);
            }
            __m256i others = min1;
            if (corr) others = q16_boxplus(min1, min2, inv, units);
            const __m256i m_others = _mm256_min_epu16(q16_scale(others, alpha_q15), msgmax);
            const __m256i m_first = _mm256_min_epu16(q16_scale(min2, alpha_q15), msgmax);
            for (e = e0; e < e1; e++) {
                const __m256i vv = v[e - e0];
                const __m256i isfirst = _mm256_cmpeq_epi16(first, _mm256_set1_epi16((short)(e - e0)));
                const __m256i mag = _mm256_blendv_epi8(m_others, m_first, isfirst);
                const __m256i sgn = _mm256_or_si256(_mm256_xor_si256(par, _mm256_cmpgt_epi16(zero, vv)), one16);
                const __m256i val = _mm256_sign_epi16(mag, sgn);
                _mm256_storeu_si256((__m256i *)(msg + (size_t)e * Q16_LANES), val);
                const __m256i nt = _mm256_min_epi16(_mm256_max_epi16(_mm256_adds_epi16(vv, val), pmin), pmax);
                _mm256_storeu_si256((__m256i *)(tot + (size_t)cidx[e] * Q16_LANES), nt);
            }
        }
    }
    return 0;
}
#endif

EXPORT int32_t ldpc_layered_q16(
    const float *llr, int32_t nblk, int32_t n, int32_t nc, int32_t nedges,
    const int32_t *cptr, const int32_t *cidx,
    int32_t iters, float alpha, float scale, int32_t corr,
    int32_t abort_k, int32_t abort_t, int32_t abort_stall,
    uint8_t *bits, uint8_t *conv, int32_t *its, int32_t *bad,
    int16_t *work)
{
#if defined(__x86_64__) || defined(_M_X64)
    if (__builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma"))
        return layered_q16_avx2(llr, nblk, n, nc, nedges, cptr, cidx, iters, alpha, scale, corr,
                                abort_k, abort_t, abort_stall, bits, conv, its, bad, work);
#endif
    (void)llr; (void)nblk; (void)n; (void)nc; (void)nedges; (void)cptr; (void)cidx;
    (void)iters; (void)alpha; (void)scale; (void)corr; (void)abort_k; (void)abort_t; (void)abort_stall;
    (void)bits; (void)conv; (void)its; (void)bad; (void)work;
    return -3;
}
