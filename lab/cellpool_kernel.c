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

/* ABI 2 adds cellpool_row_sm_f64 (the smoothed-CE row). A stale ABI-1
 * build is refused so the loader falls back to numpy rather than
 * missing a symbol at argtypes time. */
/* ABI 3 adds derot_f64 (the front end's fused de-rotation). */
EXPORT int32_t cellpool_kernel_abi(void) { return 3; }

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

/* Complex multiply in NumPy's association: (ar*br - ai*bi, ar*bi + ai*br). */
static inline void cmul(double ar, double ai, double br, double bi,
                        double *pr, double *pi)
{
    *pr = ar * br - ai * bi;
    *pi = ar * bi + ai * br;
}

/* E106: equalise ONE symbol whose channel comes from the SMOOTHED pilot grid.
 *
 * cell_pool_fast_sm built the channel as a full (nrows, ncar) complex128
 * array -- two fancy-index gathers off the grid, an interpolation, and a
 * gain multiply, ~3.4 MB per gather -- and only then divided by it. Measured
 * on x86: that H-building was 13.0 ms/Frame, MORE than the FFT it sits next
 * to. Nothing needs H to exist as an array: each data cell needs exactly one
 * value of it, so the interpolation happens here, in registers, one cell at
 * a time, and the array is never allocated.
 *
 * The arithmetic is the numpy path's, association for association:
 *   Hc = Hg[jg[c]] * (1 - wtg[c]) + Hg[jg[c]+1] * wtg[c]
 *   H  = Hc * cl                       (complex multiply, numpy's ordering)
 *   z  = Yrow[bins_d[i]] / H
 *
 * Hg : ngrid complex128 -- this symbol's column of the smoothed pilot grid
 * cl : the per-symbol complex gain that rides OUTSIDE the boxcar (E60)
 */
EXPORT int32_t cellpool_row_sm_f64(const double *Yrow, int32_t nfft,
                                   const double *Hg, int32_t ngrid,
                                   const int32_t *jg, const double *wtg,
                                   int32_t ncar,
                                   double cl_re, double cl_im,
                                   const int32_t *bins_d, const int32_t *dmap,
                                   const int32_t *hmat, int32_t nd,
                                   double *xg)
{
    if (nd < 0 || ncar < 0 || ngrid < 2) return 1;
    for (int32_t i = 0; i < nd; ++i) {
        const int32_t c = dmap[i];
        if (c < 0 || c >= ncar) return 3;
        const int32_t lo = jg[c];
        if (lo < 0 || lo + 1 >= ngrid) return 4;
        const double w = wtg[c], w1 = 1.0 - w;
        const double cr = Hg[2 * lo] * w1 + Hg[2 * (lo + 1)] * w;
        const double ci = Hg[2 * lo + 1] * w1 + Hg[2 * (lo + 1) + 1] * w;
        double hr, hi;
        cmul(cr, ci, cl_re, cl_im, &hr, &hi);

        const int32_t b = bins_d[i];
        if (b < 0 || b >= nfft) return 5;
        const int32_t t = hmat[i];
        if (t < 0 || t >= nd) return 6;
        cdiv(Yrow[2 * b], Yrow[2 * b + 1], hr, hi, &xg[2 * t], &xg[2 * t + 1]);
    }
    return 0;
}

/* E107: the front end's de-rotation, fused into ONE pass.
 *
 * m11_stream applied the cached CFO ramp and the per-chunk phasor as two
 * separate NumPy multiplies:  seg *= ramp;  seg *= scalar.  The stream is
 * complex128, so an 854k-sample block is 13.7 MB and each pass reads and
 * rewrites all of it -- measured on a Raspberry Pi 5 at 15-20 ms per block,
 * two blocks per Frame, which is most of that box's 58 ms de-rotation.
 *
 * Per element the operations and their order are the reference's --
 *   tmp = seg[i] * ramp[i];  seg[i] = tmp * scalar
 * -- and tmp stays in a register instead of being written out and read back,
 * which is the whole point: one streaming pass instead of two.
 *
 * It is NOT bit-identical to NumPy, and that was MEASURED rather than
 * assumed: NumPy's complex128 multiply reproduces neither the textbook
 * (ar*br - ai*bi, ar*bi + ai*br) nor an FMA-contracted form of it -- both
 * were implemented and compared, and both differ by ~1e-16 relative.  So the
 * bar here is the same one cpu_fast itself runs under: a tight relative
 * bound plus DECODED-BYTES identity.  On 14-bit ADC samples a 1e-15 relative
 * perturbation of the de-rotation is nine orders below the float32 the
 * demapper immediately converts to.
 *
 * y     : n complex128, modified IN PLACE
 * ramp  : n complex128 (at least n entries)
 * sr/si : the per-chunk scalar phasor
 */
EXPORT int32_t derot_f64(double *y, const double *ramp, int32_t n,
                         double sr, double si)
{
    if (n < 0) return 1;
    for (int32_t i = 0; i < n; ++i) {
        const double ar = y[2 * i], ai = y[2 * i + 1];
        const double br = ramp[2 * i], bi = ramp[2 * i + 1];
        double tr, ti;
        cmul(ar, ai, br, bi, &tr, &ti);       /* seg * ramp */
        cmul(tr, ti, sr, si, &y[2 * i], &y[2 * i + 1]);   /* * scalar */
    }
    return 0;
}
