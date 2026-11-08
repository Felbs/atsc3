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

EXPORT int32_t ldpc_layered_abi(void) { return 4; }   /* 4: + ldpc_layered_pack_f32 (E118) */

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
