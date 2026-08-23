/* qmf_kernel.c -- E109: the AC-4 A-SPX QMF analysis/synthesis bank, in C.
 *
 * WHAT THIS IS
 * ------------
 * m33_qmf.analyse and .synthesise are Python loops over TIMESLOTS -- 2880 of
 * them for 3.84 s of audio -- and each iteration issues a handful of tiny
 * NumPy calls: a 640-tap window multiply, a five-fold add, a 64x128 matvec,
 * a shift of a 640- or 1280-element register.  Profiled on a Raspberry Pi 5,
 * the two of them are 71% of A-SPX high-frequency regeneration, which is in
 * turn the largest single cost in the audio worker.  Like E53's LDPC, this is
 * dispatch-bound rather than arithmetic-bound: the work per call is far
 * smaller than the cost of making it.
 *
 * SEMANTICS (mirrored from m33_qmf; the gate is the referee):
 *   analyse, per timeslot ts:
 *     filt shifts down by 64; filt[0:64] = the next 64 input samples REVERSED
 *     z = filt * w                        (640 taps)
 *     u = sum of the five 128-sample blocks of z
 *     out[:, ts] = M @ u                  (M is 64x128 complex, u real)
 *   synthesise, per timeslot ts:
 *     filt shifts down by 128; filt[0:128] = real(N @ Q[:, ts])
 *     g gathers 64 from each 256-block and 64 from that block + 192
 *     ww = g * w; out[ts*64:(ts+1)*64] = the ten 64-sample rows of ww summed
 *
 * The matrices are passed in already built (they are pure functions of the
 * geometry and NumPy builds them once), so this kernel owns only the loop.
 * Accumulation is in the natural order; it is NOT claimed bit-identical to
 * NumPy's BLAS gemv, and the gate measures the difference rather than
 * assuming it -- see lab/gate_e109_qmf.py, which also re-runs the bank's own
 * round-trip SNR gate.
 *
 * Build: python lab/build_qmf_kernel.py
 */

#include <stddef.h>
#include <stdint.h>

#ifdef _WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif

#define NSB  64      /* num_qmf_subbands */
#define NWIN 640     /* num_qmf_win_coef */

EXPORT int32_t qmf_kernel_abi(void) { return 1; }

/* pcm : nts*NSB doubles
 * w   : NWIN doubles
 * M   : NSB x (2*NSB) complex128 = 64 x 128, row-major, as 2 doubles each
 * out : NSB x nts complex128, row-major (matches np.empty((64, nts)))
 */
EXPORT int32_t qmf_analyse_f64(const double *pcm, int32_t nts,
                               const double *w, const double *M,
                               double *out)
{
    if (nts < 0) return 1;
    double filt[NWIN];
    double z[NWIN];
    double u[2 * NSB];
    for (int32_t i = 0; i < NWIN; ++i) filt[i] = 0.0;

    for (int32_t ts = 0; ts < nts; ++ts) {
        /* shift down by NSB, then load the next 64 samples reversed */
        for (int32_t i = NWIN - 1; i >= NSB; --i) filt[i] = filt[i - NSB];
        const double *src = pcm + (size_t)ts * NSB;
        for (int32_t i = 0; i < NSB; ++i) filt[i] = src[NSB - 1 - i];

        for (int32_t i = 0; i < NWIN; ++i) z[i] = filt[i] * w[i];
        for (int32_t i = 0; i < 2 * NSB; ++i) u[i] = z[i];
        for (int32_t k = 1; k < 5; ++k)
            for (int32_t i = 0; i < 2 * NSB; ++i)
                u[i] += z[k * 2 * NSB + i];

        for (int32_t sb = 0; sb < NSB; ++sb) {
            const double *Mr = M + (size_t)sb * (2 * NSB) * 2;
            double ar = 0.0, ai = 0.0;
            for (int32_t i = 0; i < 2 * NSB; ++i) {
                ar += Mr[2 * i] * u[i];
                ai += Mr[2 * i + 1] * u[i];
            }
            const size_t o = ((size_t)sb * nts + ts) * 2;
            out[o] = ar;
            out[o + 1] = ai;
        }
    }
    return 0;
}

/* Q   : NSB x nts complex128, row-major
 * w   : NWIN doubles
 * N   : (2*NSB) x NSB complex128 = 128 x 64, row-major
 * out : nts*NSB doubles (real PCM)
 */
EXPORT int32_t qmf_synthesise_f64(const double *Q, int32_t nts,
                                  const double *w, const double *N,
                                  double *out)
{
    if (nts < 0) return 1;
    double filt[10 * 2 * NSB];      /* 1280 */
    double g[NWIN];
    for (int32_t i = 0; i < 10 * 2 * NSB; ++i) filt[i] = 0.0;

    for (int32_t ts = 0; ts < nts; ++ts) {
        for (int32_t i = 10 * 2 * NSB - 1; i >= 2 * NSB; --i)
            filt[i] = filt[i - 2 * NSB];
        /* filt[0:128] = real(N @ Q[:, ts]) */
        for (int32_t n = 0; n < 2 * NSB; ++n) {
            const double *Nr = N + (size_t)n * NSB * 2;
            double ar = 0.0;
            for (int32_t sb = 0; sb < NSB; ++sb) {
                const size_t q = ((size_t)sb * nts + ts) * 2;
                /* real part of (Nr[sb] * Q[sb, ts]) */
                ar += Nr[2 * sb] * Q[q] - Nr[2 * sb + 1] * Q[q + 1];
            }
            filt[n] = ar;
        }
        for (int32_t n = 0; n < 5; ++n) {
            for (int32_t i = 0; i < NSB; ++i)
                g[128 * n + i] = filt[256 * n + i];
            for (int32_t i = 0; i < NSB; ++i)
                g[128 * n + NSB + i] = filt[256 * n + 192 + i];
        }
        double *o = out + (size_t)ts * NSB;
        for (int32_t i = 0; i < NSB; ++i) o[i] = g[i] * w[i];
        for (int32_t r = 1; r < 10; ++r)
            for (int32_t i = 0; i < NSB; ++i)
                o[i] += g[r * NSB + i] * w[r * NSB + i];
    }
    return 0;
}
