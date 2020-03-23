
/* constants, must be set before any of packed_get1, packed_set1, packed_get2,
*      packed_set2 are used.
*  log2pixelsperbyte is raw[1] if "raw" is the input packed state vector
*      (must be 0, 1, 2, or 3)
*  log2bitsperpixel is 3 - log2pixelsperbyte
*  mask_remainder is (1 << log2pixelsperbyte) - 1
*  mask_pixel is (1 << (1 << log2bitsperpixel)) - 1
*  nrow is raw[2-5] considered as unsigned int
*  ncol is raw[6-9] considered as unsigned int
*  x is raw+10
*
*  table of values
*  log2pixelsperbyte    mask_remainder    log2bitsperpixel    mask_pixel
*              0          0000 binary               3       11111111 binary
*              1          0001 binary               2       00001111 binary
*              2          0011 binary               1       00000011 binary
*              3          0111 binary               0       00000001 binary
*/

static unsigned int log2pixelsperbyte = 0;
static unsigned int log2bitsperpixel = 0;
static unsigned int mask_remainder = 0;
static unsigned int mask_pixel = 0;
static unsigned int nrow = 0;
static unsigned int ncol = 0;
static unsigned char *x = 0;

static unsigned int packed_get1(unsigned int idx);
static void packed_set1(unsigned int idx, unsigned int val);
#ifdef OOPSIE /* clang says unused */
static inline unsigned int packed_get2(unsigned int idx_row,
    unsigned int idx_col);
static inline void packed_set2(unsigned int idx_row, unsigned int idx_col,
    unsigned int val);
#endif /* OOPSIE */

/* ripped off from src/main/saveload.c in R core distribution */
#ifdef ASSERT
/* The line below requires an ANSI C preprocessor (stringify operator) */
#define R_assert(e) ((e) ? (void) 0 : error("assertion `%s' failed: file `%s', line %d\n", #e, __FILE__, __LINE__))
#else
#define R_assert(e) ((void) 0)
#endif /* ASSERT */

static unsigned int packed_get1(unsigned int idx)
{
    register unsigned int idx_mod = idx >> log2pixelsperbyte;
    register unsigned int idx_rem = idx & mask_remainder;
    register unsigned int b = x[idx_mod];
    return (b >> (idx_rem << log2bitsperpixel)) & mask_pixel;
}

static void packed_set1(unsigned int idx, unsigned int val)
{
    register unsigned int idx_mod = idx >> log2pixelsperbyte;
    register unsigned int idx_rem = idx & mask_remainder;
    register unsigned int idx_rem_shift = idx_rem << log2bitsperpixel;
    register unsigned int m = mask_pixel << idx_rem_shift;
    register unsigned int b = (val & mask_pixel) << idx_rem_shift;
    x[idx_mod] = (x[idx_mod] & (~m)) | b;
}

static inline unsigned int packed_get2(unsigned int idx_row,
    unsigned int idx_col)
{
    return packed_get1(idx_row + nrow * idx_col);
}

static inline void packed_set2(unsigned int idx_row,
    unsigned int idx_col, unsigned int val)
{
    packed_set1(idx_row + nrow * idx_col, val);
}

#include <R_ext/Error.h>
#include "potts.h"

#ifdef BLATHER
#include <stdio.h>
#endif /* BLATHER */

void packPotts(int *myx, int *nrowin, int *ncolin, int *ncolorin,
    int *lenrawin, unsigned char *raw)
{
    nrow = nrowin[0];
    ncol = ncolin[0];
    unsigned int ncolor = ncolorin[0];
    unsigned int lenraw = lenrawin[0];

    if (nrow < 1)
        error("nrow < 1");
    if (ncol < 1)
        error("ncol < 1");
    if (ncolor < 2 || ncolor > 256)
        error("ncolor < 2 || ncolor > 256");

    log2bitsperpixel = 0;
    mask_pixel = 1;
#ifdef BLATHER
    fprintf(stderr, "ncolor = %4d\n", ncolor);
    fprintf(stderr, "(ncolor - 1) & mask_pixel = %2x (hexadecimal), ncolor - 1 = %2x (hexadecimal)\n", (ncolor - 1) & mask_pixel, ncolor - 1);
#endif /* BLATHER */
    while (((ncolor - 1) & mask_pixel) != (ncolor - 1)) {
        log2bitsperpixel++;
        mask_pixel = (1 << (1 << log2bitsperpixel)) - 1;
#ifdef BLATHER
        fprintf(stderr, "log2bitsperpixel = %4d (decimal), mask_pixel = %4x (hexadecimal)\n", log2bitsperpixel, mask_pixel);
        fprintf(stderr, "(ncolor - 1) & mask_pixel = %2x (hexadecimal), ncolor - 1 = %2x (hexadecimal)\n", (ncolor - 1) & mask_pixel, ncolor - 1);
#endif /* BLATHER */
    }
    log2pixelsperbyte = 3 - log2bitsperpixel;
    mask_remainder = (1 << log2pixelsperbyte) - 1;
    int pixelsperbyte = 1 << log2pixelsperbyte;

    if ((lenraw - 10) * pixelsperbyte < nrow * ncol)
        error("can't happen; not enough space allocated for output vector");

    raw[0] = ncolor - 1;
    raw[1] = log2pixelsperbyte;
    raw[2] = nrow >> 24;
    raw[3] = nrow >> 16;
    raw[4] = nrow >> 8;
    raw[5] = nrow;
    raw[6] = ncol >> 24;
    raw[7] = ncol >> 16;
    raw[8] = ncol >> 8;
    raw[9] = ncol;
    x = raw + 10;

#ifdef BLATHER
    printf("nrow = %4d\n", nrow);
    printf("ncol = %4d\n", ncol);
    printf("ncolor = %4d\n", ncolor);
    printf("log2pixelsperbyte = %4d\n", log2pixelsperbyte);
#endif /* BLATHER */

    for (unsigned int i = 0; i < nrow * ncol; i++)
        packed_set1(i, myx[i] - 1);
}

void inspectPotts(unsigned char *raw, int *ncolorout,
    int *nrowout, int *ncolout)
{
    ncolorout[0] = raw[0] + 1;
    int foo = 0;
    foo |= raw[2];
    foo <<= 8;
    foo |= raw[3];
    foo <<= 8;
    foo |= raw[4];
    foo <<= 8;
    foo |= raw[5];
    nrowout[0] = foo;
    foo = 0;
    foo |= raw[6];
    foo <<= 8;
    foo |= raw[7];
    foo <<= 8;
    foo |= raw[8];
    foo <<= 8;
    foo |= raw[9];
    ncolout[0] = foo;
}

void unpackPotts(unsigned char *raw, int *lenrawin, int *ncolorin,
    int *nrowin, int *ncolin, int *myx)
{
    int nrow = nrowin[0];
    int ncol = ncolin[0];
    int ncolor = ncolorin[0];
    int lenraw = lenrawin[0];

    int nrow_raw, ncol_raw, ncolor_raw;
    inspectPotts(raw, &ncolor_raw, &nrow_raw, &ncol_raw);
    if (nrow != nrow_raw)
        error("nrow for raw vector does not match argument nrow");
    if (ncol != ncol_raw)
        error("ncol for raw vector does not match argument ncol");
    if (ncolor != ncolor_raw)
        error("ncolor for raw vector does not match argument ncolor");

    log2pixelsperbyte = raw[1];
    log2bitsperpixel = 3 - log2pixelsperbyte;
    mask_remainder = (1 << log2pixelsperbyte) - 1;
    mask_pixel = (1 << (1 << log2bitsperpixel)) - 1;
    x = raw + 10;

    int pixelsperbyte = 1 << log2pixelsperbyte;

    if ((lenraw - 10) * pixelsperbyte < nrow * ncol)
        error("can't happen; input vector too short");

    for (int i = 0; i < nrow * ncol; i++)
        myx[i] = packed_get1(i) + 1;
}

#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>

#define BDRY_TORUS 1
#define BDRY_FREE 2
#define BDRY_CONDITION 3

static void equiv_init(unsigned int *ff, int nvert);
static void equiv_update(unsigned int *ff, int p0, int q0);
static void equiv_finish(unsigned int *ff, int nvert);
#ifdef OOPSIE /* gcc says unused */
static void order_alpha(double *alpha, int ncolor, int *aorder);
#endif /* OOPSIE */

static SEXP out_func = NULL;
static SEXP out_envir = NULL;

/* yes, this is a kludge, but too painful to rewrite everything now */
SEXP outfun_setup(SEXP func, SEXP envir)
{
    if (! isFunction(func))
        error("outfun_setup: argument func must be function");
    if (! isEnvironment(envir))
        error("outfun_setup: argument envir must be environment");
    out_func = func;
    out_envir = envir;
    return R_NilValue;
}

void outfun_shutdown(void)
{
    out_func = NULL;
    out_envir = NULL;
}

static SEXP outfun_do(double *tt, int len_tt)
{
    SEXP call, result, tstat, foo;

    if (! out_func)
        error("outfun: not setup yet");
    PROTECT(tstat = allocVector(REALSXP, len_tt));
    for (int i = 0; i < len_tt; i++)
        REAL(tstat)[i] = tt[i];
    PROTECT(call = lang2(out_func, tstat));
    PROTECT(result = eval(call, out_envir));
    if (! isVectorAtomic(result))
        error("outfun: result of function call must be atomic");
    if (! isNumeric(result))
        error("outfun: result of function call must be numeric");
    PROTECT(foo = coerceVector(result, REALSXP));
    int len_foo = LENGTH(foo);
    for (int i = 0; i < len_foo; i++)
        if (! R_finite(REAL(foo)[i]))
            error("outfun returned vector with non-finite element");
    UNPROTECT(4);
    return foo;
}

static int outfun_len(double *tt, int len_tt)
{
    return LENGTH(outfun_do(tt, len_tt));
}

static void outfun_val(double *tt, int len_tt, double *buff, int len_buff)
{
    SEXP foo;

    PROTECT(foo = outfun_do(tt, len_tt));
    if (LENGTH(foo) != len_buff)
        error("outfun return vector length changed from initial");
    for (int i = 0; i < len_buff; i++)
        buff[i] = REAL(foo)[i];
    UNPROTECT(1);
}

static void compute_canonical(unsigned char *xx, double *tt, int len_tt,
    int code)
{
    int ncolor = len_tt - 1;
    int nvert = nrow * ncol;

    for (int k = 0; k < ncolor; k++)
        tt[k] = 0;
    int tstar = 0;

    if (code == BDRY_TORUS) {

        for (int i = 0; i < nvert; i++)
            tt[xx[i]]++;

        for (unsigned int i = 0; i < nrow; i++)
            for (unsigned int j = 0; j < ncol; j++) {
                int idx1 = i + nrow * j;
                int idx2 = ((i + 1) % nrow) + nrow * j;
                if (xx[idx1] == xx[idx2])
                    tstar++;
            }

        for (unsigned int i = 0; i < nrow; i++)
            for (unsigned int j = 0; j < ncol; j++) {
                int idx1 = i + nrow * j;
                int idx2 = i + nrow * ((j + 1) % ncol);
                if (xx[idx1] == xx[idx2])
                    tstar++;
            }

    } else if (code == BDRY_FREE) {

        for (int i = 0; i < nvert; i++)
            tt[xx[i]]++;

        for (unsigned int i = 0; i < nrow - 1; i++)
            for (unsigned int j = 0; j < ncol; j++) {
                int idx1 = i + nrow * j;
                int idx2 = (i + 1) + nrow * j;
                if (xx[idx1] == xx[idx2])
                    tstar++;
            }

        for (unsigned int i = 0; i < nrow; i++)
            for (unsigned int j = 0; j < ncol - 1; j++) {
                int idx1 = i + nrow * j;
                int idx2 = i + nrow * (j + 1);
                if (xx[idx1] == xx[idx2])
                    tstar++;
            }

    } else if (code == BDRY_CONDITION) {

        for (unsigned int i = 1; i < nrow - 1; i++)
            for (unsigned int j = 1; j < ncol - 1; j++)
                tt[xx[i + nrow * j]]++;

        for (unsigned int i = 0; i < nrow - 1; i++)
            for (unsigned int j = 1; j < ncol - 1; j++) {
                int idx1 = i + nrow * j;
                int idx2 = (i + 1) + nrow * j;
                if (xx[idx1] == xx[idx2])
                    tstar++;
            }

        for (unsigned int i = 1; i < nrow - 1; i++)
            for (unsigned int j = 0; j < ncol - 1; j++) {
                int idx1 = i + nrow * j;
                int idx2 = i + nrow * (j + 1);
                if (xx[idx1] == xx[idx2])
                    tstar++;
            }

    } else /* Can't happen */ {
        error("impossible value of code");
    }

    tt[ncolor] = tstar;
}

void outfun_len_init(unsigned char *raw, int *codein, int *nout)
{
    int code = codein[0];
    int nrow_raw, ncol_raw, ncolor_raw;
    inspectPotts(raw, &ncolor_raw, &nrow_raw, &ncol_raw);
    int nvert = nrow_raw * ncol_raw;
    int nparm = ncolor_raw + 1;

    unsigned char *xx = (unsigned char *) R_alloc(nvert, sizeof(unsigned char));
    double *tt = (double *) R_alloc(nparm, sizeof(double));

    for (int i = 0; i < nvert; i++)
        xx[i] = packed_get1(i);
    compute_canonical(xx, tt, nparm, code);
    nout[0] = outfun_len(tt, nparm);
}

void potts(unsigned char *raw, double *theta, int *nbatchin, int *blenin,
    int *nspacin, int *codein, double *batch, int *debugin, int *pstate,
    int *hstate, int *vstate, int *patch, double *hunif, double *vunif,
    double *punif)
{
    int nbatch = nbatchin[0];
    int blen = blenin[0];
    int nspac = nspacin[0];
    int code = codein[0];
    int is_debug = debugin[0];

    if (! (code == BDRY_TORUS || code == BDRY_FREE || code == BDRY_CONDITION))
        error("Can't happen: integer code for boundary conditions bad\n");

    int nrow_raw, ncol_raw, ncolor_raw;
    inspectPotts(raw, &ncolor_raw, &nrow_raw, &ncol_raw);
    log2pixelsperbyte = raw[1];
    log2bitsperpixel = 3 - log2pixelsperbyte;
    mask_remainder = (1 << log2pixelsperbyte) - 1;
    mask_pixel = (1 << (1 << log2bitsperpixel)) - 1;
    nrow = nrow_raw;
    ncol = ncol_raw;
    x = raw + 10;
    int ncolor = ncolor_raw;

    double *alpha = theta;
    double beta = theta[ncolor];
    if (beta < 0.0)
        error("Swendsen-Wang algorithm can only handle nonnegative theta[ncolor + 1]");

    double bprob = (- expm1(- beta));    /* conditional probability of bond */

    int nvert = nrow * ncol;
    unsigned int *ff = (unsigned int *) R_alloc(nvert, sizeof(unsigned int));
    unsigned int *gg = (unsigned int *) R_alloc(nvert, sizeof(unsigned int));
    unsigned char *xx = (unsigned char *) R_alloc(nvert, sizeof(unsigned char));
    unsigned char *ss = (unsigned char *) R_alloc(nvert, sizeof(unsigned char));
    unsigned char *cc = (unsigned char *) R_alloc(nvert, sizeof(unsigned char));
#ifdef OOPSIE /* gcc says unused */
    int *aorder = (int *) R_alloc(ncolor, sizeof(int));
#endif /* OOPSIE */
    double *pp = (double *) R_alloc(ncolor, sizeof(double));
    double alpha_max;
    int alpha_is_zero;

    for (int i = 0; i < nvert; i++)
        xx[i] = packed_get1(i);

    GetRNGstate();

    alpha_is_zero = 1;
    for (int k = 0; k < ncolor; k++)
        alpha_is_zero &= (alpha[k] == 0.0);

    alpha_max = alpha[0];
    for (int k = 0; k < ncolor; k++)
        if (alpha_max < alpha[k])
            alpha_max = alpha[k];

    // need buffers for batch means, canonical statistic vector, and outfun

    int nparm = ncolor + 1;
    double *tt = (double *) R_alloc(nparm, sizeof(double));
    compute_canonical(xx, tt, nparm, code);

    int nout = nparm;
    if (out_func)
        nout = outfun_len(tt, nparm);
    double *batch_buff = (double *) R_alloc(nout, sizeof(double));
    double *outfun_buff = (double *) R_alloc(nout, sizeof(double));

    int niter = nbatch * blen * nspac;

    for (int ibatch = 0, iiter = 0; ibatch < nbatch; ibatch++) {

    for (int i = 0; i < nout; i++)
        batch_buff[i] = 0.0;

    for (int jbatch = 0; jbatch < blen; jbatch++) {

    for (int ispac = 0; ispac < nspac; ispac++, iiter++) {

        if (is_debug)
            for (unsigned int i = 0; i < nrow; i++)
                for (unsigned int j = 0; j < ncol; j++)
                    pstate[iiter + niter * (i + nrow * j)] =
                        xx[i + nrow * j] + 1;

        /* do bonds */

        equiv_init(ff, nvert);

        if (code == BDRY_TORUS) {

            for (unsigned int i = 0; i < nrow; i++)
                for (unsigned int j = 0; j < ncol; j++) {
                    int idx1 = i + nrow * j;
                    int idx2 = ((i + 1) % nrow) + nrow * j;
                    R_assert(0 <= idx1 && idx1 < nvert);
                    R_assert(0 <= idx2 && idx1 < nvert);
                    if (xx[idx1] == xx[idx2]) {
                        double u = unif_rand();
                        if (u < bprob) {
                            equiv_update(ff, idx1, idx2);
                            if (is_debug)
                                hstate[iiter + niter * (i + nrow * j)] = 1;
                        }
                        if (is_debug)
                            hunif[iiter + niter * (i + nrow * j)] = u;
                    }
                }

            for (unsigned int i = 0; i < nrow; i++)
                for (unsigned int j = 0; j < ncol; j++) {
                    int idx1 = i + nrow * j;
                    int idx2 = i + nrow * ((j + 1) % ncol);
                    R_assert(0 <= idx1 && idx1 < nvert);
                    R_assert(0 <= idx2 && idx1 < nvert);
                    if (xx[idx1] == xx[idx2]) {
                        double u = unif_rand();
                        if (u < bprob) {
                            equiv_update(ff, idx1, idx2);
                            if (is_debug)
                                vstate[iiter + niter * (i + nrow * j)] = 1;
                        }
                        if (is_debug)
                            vunif[iiter + niter * (i + nrow * j)] = u;
                    }
                }

        } else if (code == BDRY_FREE) {

            for (unsigned int i = 0; i < nrow - 1; i++)
                for (unsigned int j = 0; j < ncol; j++) {
                    int idx1 = i + nrow * j;
                    int idx2 = (i + 1) + nrow * j;
                    R_assert(0 <= idx1 && idx1 < nvert);
                    R_assert(0 <= idx2 && idx1 < nvert);
                    if (xx[idx1] == xx[idx2]) {
                        double u = unif_rand();
                        if (u < bprob) {
                            equiv_update(ff, idx1, idx2);
                            if (is_debug)
                                hstate[iiter + niter * (i + nrow * j)] = 1;
                        }
                        if (is_debug)
                            hunif[iiter + niter * (i + nrow * j)] = u;
                    }
                }

            for (unsigned int i = 0; i < nrow; i++)
                for (unsigned int j = 0; j < ncol - 1; j++) {
                    int idx1 = i + nrow * j;
                    int idx2 = i + nrow * (j + 1);
                    R_assert(0 <= idx1 && idx1 < nvert);
                    R_assert(0 <= idx2 && idx1 < nvert);
                    if (xx[idx1] == xx[idx2]) {
                        double u = unif_rand();
                        if (u < bprob) {
                            equiv_update(ff, idx1, idx2);
                            if (is_debug)
                                vstate[iiter + niter * (i + nrow * j)] = 1;
                        }
                        if (is_debug)
                            vunif[iiter + niter * (i + nrow * j)] = u;
                    }
                }

        } else if (code == BDRY_CONDITION) {

            for (unsigned int i = 0; i < nrow - 1; i++)
                for (unsigned int j = 1; j < ncol - 1; j++) {
                    int idx1 = i + nrow * j;
                    int idx2 = (i + 1) + nrow * j;
                    R_assert(0 <= idx1 && idx1 < nvert);
                    R_assert(0 <= idx2 && idx1 < nvert);
                    if (xx[idx1] == xx[idx2]) {
                        double u = unif_rand();
                        if (u < bprob) {
                            equiv_update(ff, idx1, idx2);
                            if (is_debug)
                                hstate[iiter + niter * (i + nrow * j)] = 1;
                        }
                        if (is_debug)
                            hunif[iiter + niter * (i + nrow * j)] = u;
                    }
                }

            for (unsigned int i = 1; i < nrow - 1; i++)
                for (unsigned int j = 0; j < ncol - 1; j++) {
                    int idx1 = i + nrow * j;
                    int idx2 = i + nrow * (j + 1);
                    R_assert(0 <= idx1 && idx1 < nvert);
                    R_assert(0 <= idx2 && idx1 < nvert);
                    if (xx[idx1] == xx[idx2]) {
                        double u = unif_rand();
                        if (u < bprob) {
                            equiv_update(ff, idx1, idx2);
                            if (is_debug)
                                vstate[iiter + niter * (i + nrow * j)] = 1;
                        }
                        if (is_debug)
                            vunif[iiter + niter * (i + nrow * j)] = u;
                    }
                }

        } else /* Can't happen */ {
            error("impossible value of code");
        }

        equiv_finish(ff, nvert);

        /* Now vertices are in same patch if and only if have same ff[j] */

        if (is_debug)
            for (unsigned int i = 0; i < nrow; i++)
                for (unsigned int j = 0; j < ncol; j++)
                    patch[iiter + niter * (i + nrow * j)] =
                        ff[i + nrow * j] + 1;

        /* Count pixels in patches */

        for (int i = 0; i < nvert; i++)
            gg[i] = 0;
        for (int i = 0; i < nvert; i++)
            gg[ff[i]]++;

        /* Now patch ff[j] has size gg[j] */

        /* Determine which patches are fixed by conditioning on boundary */

        for (int i = 0; i < nvert; i++)
            ss[i] = 0;
        if (code == BDRY_CONDITION) {

            for (unsigned int i = 0; i < nrow; i++) {
                int idx1 = i + nrow * 0;
                int idx2 = i + nrow * (ncol - 1);
                R_assert(0 <= idx1 && idx1 < nvert);
                R_assert(0 <= idx2 && idx1 < nvert);
                ss[ff[idx1]] = 1;
                ss[ff[idx2]] = 1;
            }
            for (unsigned int j = 0; j < ncol; j++) {
                int idx1 = 0 + nrow * j;
                int idx2 = (nrow - 1) + nrow * j;
                R_assert(0 <= idx1 && idx1 < nvert);
                R_assert(0 <= idx2 && idx1 < nvert);
                ss[ff[idx1]] = 1;
                ss[ff[idx2]] = 1;
            }
        }

        /* Now patch ff[j] is fixed if and only if ss[j] == 1 */

        /* Determine new patch colors */

        for (int i = 0; i < nvert; i++)
            if (gg[i] > 0) /* there is patch number i */ {
                int newcolor;
                double u = unif_rand();
                if (alpha_is_zero) {
                    newcolor = u * ncolor;
                    if (newcolor >= ncolor)
                        newcolor--;
                } else /* alphas not all zero */ {
                    for (int k = 0; k < ncolor; k++)
                        pp[k] = exp(gg[i] * (alpha[k] - alpha_max));
                    for (int k = 1; k < ncolor; k++)
                        pp[k] += pp[k - 1];
                    double psum = pp[ncolor - 1];
                    for (int k = 0; k < ncolor; k++) {
                        newcolor = k;
                        if (u < pp[k] / psum)
                            break;
                    }
                }
                cc[i] = newcolor;
                if (is_debug)
                    punif[iiter + niter * i] = u;
            }

        /* now cc[j] is the new color of the patch numbered j */

        for (int i = 0; i < nvert; i++)
            if (ss[ff[i]] == 0)
                xx[i] = cc[ff[i]];

        /* compute canonical statistics */

        compute_canonical(xx, tt, nparm, code);

        R_CheckUserInterrupt();

    } /* end of inner loop (one iteration) */

    if (out_func) {
        outfun_val(tt, nparm, outfun_buff, nout);
        for (int i = 0; i < nout; i++)
            batch_buff[i] += outfun_buff[i];
    } else {
        for (int i = 0; i < nparm; i++)
            batch_buff[i] += tt[i];
    }

    } /* end of middle loop (one batch) */

    for (int i = 0; i < nout; i++)
        batch[i + nout * ibatch] = batch_buff[i] / blen;

    } /* end of outer loop */

    PutRNGstate();

    for (int i = 0; i < nvert; i++)
        packed_set1(i, xx[i]);
}

static void equiv_init(unsigned int *ff, int nvert)
{
    for (int i = 0; i < nvert; i++)
        ff[i] = i;
}

static void equiv_update(unsigned int *ff, int p0, int q0)
{
    int p1, q1;

    p1 = ff[p0];
    q1 = ff[q0];

    while (p1 != q1)
        if (q1 < p1) {
            ff[p0] = q1;
            p0 = p1;
            p1 = ff[p1];
        } else /* q1 > p1 */ {
            ff[q0] = p1;
            q0 = q1;
            q1 = ff[q1];
        }
}

static void equiv_finish(unsigned int *ff, int nvert)
{
    for (int i = 0; i < nvert; i++)
        ff[i] = ff[ff[i]];
}

#ifdef OOPSIE /* gcc says unused */
static void order_alpha(double *alpha, int ncolor, int *aorder)
{
    for (int i = 0; i < ncolor; i++)
        aorder[i] = i;

    for (int i = 0; i < ncolor; i++) {
        int tmp = aorder[i];
        double amin = alpha[aorder[i]];
        int imin = i;
        for (int j = i + 1; j < ncolor; j++)
            if (alpha[aorder[j]] < amin) {
                amin = alpha[aorder[j]];
                imin = j;
            }
            aorder[i] = aorder[imin];
            aorder[imin] = tmp;
    }
}
#endif /* OOPSIE */

