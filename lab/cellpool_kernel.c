/* cellpool_kernel.c -- E103: compiled per-symbol equalise+scatter, complex128.
 *
 * WHAT THIS IS
 * ------------
 * The per-class inner work of m9_fast.FrameDecoder.cell_pool_fast, written as
 * plain C.  MEASURED before it was written (x86, RF33, 8 reps): cell_pool is
 * 20.1 ms/Frame, of which the scipy FFT is 4.7 ms (23%) -- untouchable, it is
 * already pocketfft -- and the per-class work is 14.6 ms (73%).  That work is
 * ~240,000 complex128 divisions per Frame plus four fancy-index passes
 * (pilot gather, interpolation, data gather, scatter), each materialising a
 * fresh array.  This kernel fuses all of it into one sweep per symbol: the
 * channel is interpolated ON THE FLY per data cell, so H is never built, and
 * the only array written is the de-interleaved output.
 *
 * SEMANTICS (mirrored from cell_pool_fast; the gate is the referee):
 *   hp[k] = Yrow[bins_pk[k]] / refpk[k]
 *   H(c)  = hp[jj[c]] * (1 - wt[c]) + hp[jj[c]+1] * wt[c]
 *   z[i]  = Yrow[bins_d[i]] / H(dmap[i])
 *   xg[hmat[i]] = z[i]
 * The interpolation is evaluated in exactly the numpy path's association --
 * hp[lo]*(1-w) + hp[hi]*w -- which is itself cpu_fast's standing
 * re-association of np.interp, decoded-bytes gated.
 *
 * COMPLEX DIVISION uses Smith's algorithm (scale by the smaller part), the
 * same overflow-safe form the C runtime and numpy use, rather than the naive
 * (a*conj(b))/|b|^2.  The naive form is faster but changes both the rounding
 * and the overflow behaviour, and a channel estimate can legitimately be
 * tiny; a demapper that divides by a deep fade must not produce an infinity
 * where numpy produced a number.
 *
 * DELIBERATELY NOT HERE: threads (cell_pool_fast's ThreadPoolExecutor already
 * maps classes across cores and ctypes drops the GIL), Python.h, and
 * -ffast-math (the build pins -ffp-contract=off).
 *
 * Build: python lab/build_cellpool_kernel.py
 */

#include <stddef.h>
#include <stdint.h>

#ifdef _WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif

EXPORT int32_t cellpool_kernel_abi(void) { return 1; }

/* Smith's algorithm: (ar + i*ai) / (br + i*bi) -> (*qr, *qi). */
static inline void cdiv(double ar, double ai, double br, double bi,
                        double *qr, double *qi)
{
    double r, den;
    double babs_r = br < 0.0 ? -br : br;
    double babs_i = bi < 0.0 ? -bi : bi;
    /* NumPy's nc_quot multiplies by a RECIPROCAL rather than dividing twice
     * (scl = 1/den, then two multiplies).  That is a different rounding, and
     * the gate measured it: dividing twice left the pool ~1 ulp off real air.
     * Mirror numpy exactly so the pool comes back bit-identical. */
    if (babs_r >= babs_i) {
        r = bi / br;
        den = 1.0 / (br + bi * r);
        *qr = (ar + ai * r) * den;
        *qi = (ai - ar * r) * den;
    } else {
        r = br / bi;
        den = 1.0 / (br * r + bi);
        *qr = (ar * r + ai) * den;
        *qi = (ai * r - ar) * den;
    }
}

/* Equalise and de-interleave ONE OFDM symbol.
 *
 * Yrow    : nfft complex128 (2*nfft doubles) -- one FFT-shifted symbol
 * bins_pk : npk int32  -- pilot bin indices into Yrow
 * refpk   : npk complex128 -- the reference pilot values
 * jj, wt  : ncar each -- lower pilot index and interpolation weight per carrier
 * bins_d  : nd int32   -- data bin indices into Yrow
 * dmap    : nd int32   -- carrier index of each data cell (into jj/wt)
 * hmat    : nd int32   -- frequency-de-interleave target slot of each data cell
 * hp      : npk complex128 scratch (caller-owned, reused across symbols)
 * xg      : nd complex128 output, de-interleaved
 *
 * Returns 0, or non-zero if an index would leave its array (the caller then
 * falls back to numpy rather than reading out of bounds).
 */
EXPORT int32_t cellpool_row_f64(const double *Yrow, int32_t nfft,
                                const int32_t *bins_pk, const double *refpk,
                                int32_t npk,
                                const int32_t *jj, const double *wt,
                                int32_t ncar,
                                const int32_t *bins_d, const int32_t *dmap,
                                const int32_t *hmat, int32_t nd,
                                double *hp, double *xg)
{
    if (npk < 2 || nd < 0 || ncar < 0) return 1;

    for (int32_t k = 0; k < npk; ++k) {
        const int32_t b = bins_pk[k];
        if (b < 0 || b >= nfft) return 2;
        cdiv(Yrow[2 * b], Yrow[2 * b + 1],
             refpk[2 * k], refpk[2 * k + 1],
             &hp[2 * k], &hp[2 * k + 1]);
    }

    for (int32_t i = 0; i < nd; ++i) {
        const int32_t c = dmap[i];
        if (c < 0 || c >= ncar) return 3;
        const int32_t lo = jj[c];
        if (lo < 0 || lo + 1 >= npk) return 4;
        const double w = wt[c], w1 = 1.0 - w;
        /* H = hp[lo]*(1-w) + hp[lo+1]*w  -- the numpy association */
        const double hr = hp[2 * lo] * w1 + hp[2 * (lo + 1)] * w;
        const double hi = hp[2 * lo + 1] * w1 + hp[2 * (lo + 1) + 1] * w;

        const int32_t b = bins_d[i];
        if (b < 0 || b >= nfft) return 5;
        const int32_t t = hmat[i];
        if (t < 0 || t >= nd) return 6;
        cdiv(Yrow[2 * b], Yrow[2 * b + 1], hr, hi, &xg[2 * t], &xg[2 * t + 1]);
    }
    return 0;
}

/* Equalise and de-interleave ONE symbol against an ALREADY-BUILT channel.
 *
 * cell_pool_fast_sm (E58/E60 smoothed channel estimation) computes H itself
 * -- from the boxcar-smoothed pilot grid, or from the symbol's own pilots
 * when the change detector sends that symbol down the per-symbol path -- and
 * then does exactly the same tail as cell_pool_fast: gather bins_d, divide by
 * H[dmap], scatter through hmat.  That tail is where the ~208k complex
 * divisions per Frame live, so it gets the same treatment.  H stays numpy's
 * business: the smoothing, the per-symbol gain and the detector are policy,
 * not arithmetic to be hidden in C.
 *
 * H  : ncar complex128 -- channel per carrier, already interpolated
 * xg : nd complex128 output, de-interleaved
 */
EXPORT int32_t cellpool_div_scatter_f64(const double *Yrow, int32_t nfft,
                                        const double *H, int32_t ncar,
                                        const int32_t *bins_d,
                                        const int32_t *dmap,
                                        const int32_t *hmat, int32_t nd,
                                        double *xg)
{
    if (nd < 0 || ncar < 0) return 1;
    for (int32_t i = 0; i < nd; ++i) {
        const int32_t c = dmap[i];
        if (c < 0 || c >= ncar) return 3;
        const int32_t b = bins_d[i];
        if (b < 0 || b >= nfft) return 5;
        const int32_t t = hmat[i];
        if (t < 0 || t >= nd) return 6;
        cdiv(Yrow[2 * b], Yrow[2 * b + 1], H[2 * c], H[2 * c + 1],
             &xg[2 * t], &xg[2 * t + 1]);
    }
    return 0;
}
