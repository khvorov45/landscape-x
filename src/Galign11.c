#include "mltaln.h"

#define DEBUG 0
#define XXXXXXX 0
#define USE_PENALTY_EX 1

#define TERMGAPFAC 0.0
#define TERMGAPFAC_EX 0.0

static void
match_calc_mtx(double** mtx, double* match, char** s1, char** s2, int i1, int lgth2) {
    char*   seq2 = s2[0];
    double* doubleptr = mtx[(unsigned char)s1[0][i1]];

    while (lgth2--)
        *match++ = doubleptr[(unsigned char)*seq2++];
}

static double
Atracking(aln_Opts opts, Context* ctx, double* lasthorizontalw, double* lastverticalw, char** seq1, char** seq2, char** mseq1, char** mseq2, int** ijp, int tailgp, int* warpis, int* warpjs, int warpbase) {
    int i, j, l, iin, jin, ifi, jfi, lgth1, lgth2, k, limk;
    //	char gap[] = "-";
    char* gap;
    gap = ctx->newgapstr;
    lgth1 = strlen(seq1[0]);
    lgth2 = strlen(seq2[0]);
    double wm, g;
    double fpenalty = (double)opts.penalty;
    double fpenalty_ex = (double)ctx->penalty_ex;

#if 0
	for( i=0; i<lgth1; i++ ) 
	{
		fprintf( stderr, "lastverticalw[%d] = %f\n", i, lastverticalw[i] );
	}
#endif

    for (i = 0; i < lgth1 + 1; i++) {
        ijp[i][0] = i + 1;
    }
    for (j = 0; j < lgth2 + 1; j++) {
        ijp[0][j] = -(j + 1);
    }

    //	if( tailgp == 1 || ijp[lgth1][lgth2] >= warpbase )
    if (tailgp == 1)
        ;
    else {
#if 1
        //		reporterr( "lastverticalw[lgth1-1] = %f\n", lastverticalw[lgth1-1] );
        //		reporterr( "lasthorizontalw[lgth2-1] = %f\n", lasthorizontalw[lgth2-1] );
        wm = lasthorizontalw[lgth2 - 1] - 1.0;  // lasthorizontalw[lgth2-1] yori kanarazu chiisai.
        for (j = lgth2 - 2; j >= 0; j--) {
            if ((g = lasthorizontalw[j] + (fpenalty * TERMGAPFAC + fpenalty_ex * (lgth2 - 1 - j) * TERMGAPFAC_EX)) > wm) {
                wm = g;
                iin = lgth1 - 1;
                jin = j;
                ijp[lgth1][lgth2] = -(lgth2 - j);
            }
        }
        for (i = lgth1 - 2; i >= 0; i--) {
            if ((g = lastverticalw[i] + (fpenalty * TERMGAPFAC + fpenalty_ex * (lgth1 - 1 - i) * TERMGAPFAC_EX)) > wm) {
                wm = g;
                iin = i;
                jin = lgth2 - 1;
                ijp[lgth1][lgth2] = +(lgth1 - i);
            }
        }
        if (lasthorizontalw[lgth2 - 1] > wm)  // score ga onaji baai erabarenai
        {
            wm = lasthorizontalw[lgth2 - 1];
            iin = lgth1 - 1;
            jin = lgth2 - 1;
            ijp[lgth1][lgth2] = 0;
        }
#else
        wm = lastverticalw[0];
        for (i = 0; i < lgth1; i++) {
            if (lastverticalw[i] >= wm) {
                wm = lastverticalw[i];
                iin = i;
                jin = lgth2 - 1;
                ijp[lgth1][lgth2] = +(lgth1 - i);
            }
        }
        for (j = 0; j < lgth2; j++) {
            if (lasthorizontalw[j] >= wm) {
                wm = lasthorizontalw[j];
                iin = lgth1 - 1;
                jin = j;
                ijp[lgth1][lgth2] = -(lgth2 - j);
            }
        }
#endif
    }

    mseq1[0] += lgth1 + lgth2;
    *mseq1[0] = 0;
    mseq2[0] += lgth1 + lgth2;
    *mseq2[0] = 0;

    iin = lgth1;
    jin = lgth2;
    limk = lgth1 + lgth2 + 1;
    for (k = 0; k < limk; k++) {
        if (ijp[iin][jin] >= warpbase) {
            //			fprintf( stderr, "WARP!\n" );
            ifi = warpis[ijp[iin][jin] - warpbase];
            jfi = warpjs[ijp[iin][jin] - warpbase];
        } else if (ijp[iin][jin] < 0) {
            ifi = iin - 1;
            jfi = jin + ijp[iin][jin];
        } else if (ijp[iin][jin] > 0) {
            ifi = iin - ijp[iin][jin];
            jfi = jin - 1;
        } else {
            ifi = iin - 1;
            jfi = jin - 1;
        }

        if (ifi == -warpbase && jfi == -warpbase) {
            l = iin;
            while (--l >= 0) {
                *--mseq1[0] = seq1[0][l];
                *--mseq2[0] = *gap;
                k++;
            }
            l = jin;
            while (--l >= 0) {
                *--mseq1[0] = *gap;
                *--mseq2[0] = seq2[0][l];
                k++;
            }
            break;
        } else {
            l = iin - ifi;
            while (--l > 0) {
                *--mseq1[0] = seq1[0][ifi + l];
                *--mseq2[0] = *gap;
                k++;
            }
            l = jin - jfi;
            while (--l > 0) {
                *--mseq1[0] = *gap;
                *--mseq2[0] = seq2[0][jfi + l];
                k++;
            }
        }
        if (iin <= 0 || jin <= 0)
            break;
        *--mseq1[0] = seq1[0][ifi];
        *--mseq2[0] = seq2[0][jfi];
        k++;
        iin = ifi;
        jin = jfi;
    }

    //	fprintf( stderr, "%s\n", mseq1[0] );
    //	fprintf( stderr, "%s\n", mseq2[0] );
    return (wm);
}

double
G__align11(aln_Opts opts, Context* ctx, double** n_dynamicmtx, char** seq1, char** seq2, int alloclen, int headgp, int tailgp) {
    //	int k;
    register int i, j;
    int          lasti; /* outgap == 0 -> lgth1, outgap == 1 -> lgth1+1 */
    int          lastj;
    int          lgth1, lgth2;
    int          resultlen;
    double       wm, wmo; /* int ?????? */
    double       g;
    double *     currentw, *previousw;
    double       fpenalty = (double)opts.penalty;
#if USE_PENALTY_EX
    double fpenalty_ex = (double)ctx->penalty_ex;
    double fpenalty_ex_i;
#endif
#if 1
    double* wtmp;
    int*    ijppt;
    double *mjpt, *prept, *curpt;
    int*    mpjpt;
#endif
    static double   mi = 0.0;
    static double*  m = NULL;
    static int**    ijp = NULL;
    static int      mpi = 0;
    static int*     mp = NULL;
    static double*  w1 = NULL;
    static double*  w2 = NULL;
    static double*  match = NULL;
    static double*  initverticalw = NULL; /* kufuu sureba iranai */
    static double*  lastverticalw = NULL; /* kufuu sureba iranai */
    static char**   mseq1 = NULL;
    static char**   mseq2 = NULL;
    static char**   mseq = NULL;
    static int**    intwork = NULL;
    static double** doublework = NULL;
    static int      orlgth1 = 0, orlgth2 = 0;
    static double** amino_dynamicmtx = NULL;  // ??

    int*    warpis = NULL;
    int*    warpjs = NULL;
    int     warpbase;

    if (seq1 == NULL) {
        if (orlgth1 > 0 && orlgth2 > 0) {
            orlgth1 = 0;
            orlgth2 = 0;
            if (mseq1)
                free(mseq1);
            mseq1 = NULL;
            if (mseq2)
                free(mseq2);
            mseq2 = NULL;
            if (w1)
                FreeFloatVec(w1);
            w1 = NULL;
            if (w2)
                FreeFloatVec(w2);
            w2 = NULL;
            if (match)
                FreeFloatVec(match);
            match = NULL;
            if (initverticalw)
                FreeFloatVec(initverticalw);
            initverticalw = NULL;
            if (lastverticalw)
                FreeFloatVec(lastverticalw);
            lastverticalw = NULL;

            if (m)
                FreeFloatVec(m);
            m = NULL;
            if (mp)
                FreeIntVec(mp);
            mp = NULL;

            if (mseq)
                FreeCharMtx(mseq);
            mseq = NULL;

            if (doublework)
                FreeFloatMtx(doublework);
            doublework = NULL;
            if (intwork)
                FreeIntMtx(intwork);
            intwork = NULL;

            if (amino_dynamicmtx)
                FreeDoubleMtx(amino_dynamicmtx);
            amino_dynamicmtx = NULL;
        }
        orlgth1 = 0;
        orlgth2 = 0;
        return (0.0);
    }

    lgth1 = strlen(seq1[0]);
    lgth2 = strlen(seq2[0]);

    warpbase = lgth1 + lgth2;
    warpis = NULL;
    warpjs = NULL;

    if (lgth1 == 0 && lgth2 == 0)
        return (0.0);

    if (lgth1 == 0) {
        seq1[0][lgth2] = 0;
        while (lgth2)
            seq1[0][--lgth2] = *ctx->newgapstr;
        //		fprintf( stderr, "seq1[0] = %s\n", seq1[0] );
        return (0.0);
    }

    if (lgth2 == 0) {
        seq2[0][lgth1] = 0;
        while (lgth1)
            seq2[0][--lgth1] = *ctx->newgapstr;
        //		fprintf( stderr, "seq2[0] = %s\n", seq2[0] );
        return (0.0);
    }

    wm = 0.0;

    if (orlgth1 == 0) {
        mseq1 = AllocateCharMtx(2, 0);  // 2020/Apr
        mseq2 = AllocateCharMtx(2, 0);  // 2020/Apr
    }

    if (lgth1 > orlgth1 || lgth2 > orlgth2) {
        int ll1, ll2;

        if (orlgth1 > 0 && orlgth2 > 0) {
            FreeFloatVec(w1);
            FreeFloatVec(w2);
            FreeFloatVec(match);
            FreeFloatVec(initverticalw);
            FreeFloatVec(lastverticalw);

            FreeFloatVec(m);
            FreeIntVec(mp);

            FreeCharMtx(mseq);

            FreeFloatMtx(doublework);
            FreeIntMtx(intwork);
            FreeDoubleMtx(amino_dynamicmtx);
        }

        ll1 = MAX((int)(1.3 * lgth1), orlgth1) + 100;
        ll2 = MAX((int)(1.3 * lgth2), orlgth2) + 100;

#if DEBUG
        fprintf(stderr, "\ntrying to allocate (%d+%d)xn matrices ... ", ll1, ll2);
#endif

        w1 = AllocateFloatVec(ll2 + 2);
        w2 = AllocateFloatVec(ll2 + 2);
        match = AllocateFloatVec(ll2 + 2);

        initverticalw = AllocateFloatVec(ll1 + 2);
        lastverticalw = AllocateFloatVec(ll1 + 2);

        m = AllocateFloatVec(ll2 + 2);
        mp = AllocateIntVec(ll2 + 2);

        mseq = AllocateCharMtx(2, ll1 + ll2);  // 2020/Apr

        doublework = AllocateFloatMtx(ctx->nalphabets, MAX(ll1, ll2) + 2);
        intwork = AllocateIntMtx(ctx->nalphabets, MAX(ll1, ll2) + 2);
        amino_dynamicmtx = AllocateDoubleMtx(0x100, 0x100);

#if DEBUG
        fprintf(stderr, "succeeded\n");
#endif

        orlgth1 = ll1 - 100;
        orlgth2 = ll2 - 100;
    }
    for (i = 0; i < ctx->nalphabets; i++)
        for (j = 0; j < ctx->nalphabets; j++)
            amino_dynamicmtx[(uint8_t)ctx->amino[i]][(uint8_t)ctx->amino[j]] = (double)n_dynamicmtx[i][j];

    mseq1[0] = mseq[0];
    mseq2[0] = mseq[1];

    if (orlgth1 > ctx->commonAlloc1 || orlgth2 > ctx->commonAlloc2) {
        int ll1, ll2;

        if (ctx->commonAlloc1 && ctx->commonAlloc2) {
            FreeIntMtx(ctx->commonIP);
        }

        ll1 = MAX(orlgth1, ctx->commonAlloc1);
        ll2 = MAX(orlgth2, ctx->commonAlloc2);

#if DEBUG
        fprintf(stderr, "\n\ntrying to allocate %dx%d matrices ... ", ll1 + 1, ll2 + 1);
#endif

        ctx->commonIP = AllocateIntMtx(ll1 + 10, ll2 + 10);

#if DEBUG
        fprintf(stderr, "succeeded\n\n");
#endif

        ctx->commonAlloc1 = ll1;
        ctx->commonAlloc2 = ll2;
    }
    ijp = ctx->commonIP;

#if 0
	for( i=0; i<lgth1; i++ ) 
		fprintf( stderr, "ogcp1[%d]=%f\n", i, ogcp1[i] );
#endif

    currentw = w1;
    previousw = w2;

    match_calc_mtx(amino_dynamicmtx, initverticalw, seq2, seq1, 0, lgth1);
    match_calc_mtx(amino_dynamicmtx, currentw, seq1, seq2, 0, lgth2);

    if (headgp == 1) {
        for (i = 1; i < lgth1 + 1; i++) {
            initverticalw[i] += fpenalty;
#if USE_PENALTY_EX
//			initverticalw[i] += fpenalty_ex * i; // ato de fukkatsu
//			reporterr( "added _ex\n" );
#endif
        }
        for (j = 1; j < lgth2 + 1; j++) {
            currentw[j] += fpenalty;
#if USE_PENALTY_EX
//			currentw[j] += fpenalty_ex * j; // ato de fukkatsu
//			reporterr( "added _ex\n" );
#endif
        }
    } else {
        for (i = 1; i < lgth1 + 1; i++) {
            initverticalw[i] += fpenalty * TERMGAPFAC;
#if USE_PENALTY_EX  // 2018/Apr/22
            initverticalw[i] += fpenalty_ex * i * TERMGAPFAC_EX;
//			reporterr( "added _ex\n" );
#endif
        }
        for (j = 1; j < lgth2 + 1; j++) {
            currentw[j] += fpenalty * TERMGAPFAC;
#if USE_PENALTY_EX  // 2018/Apr/22
            currentw[j] += fpenalty_ex * j * TERMGAPFAC_EX;
//			reporterr( "added _ex\n" );
#endif
        }
    }

    for (j = 1; j < lgth2 + 1; ++j) {
        m[j] = currentw[j - 1];
        mp[j] = 0;
    }

    if (lgth2 == 0)
        lastverticalw[0] = 0.0;  // lgth2==0 no toki error
    else
        lastverticalw[0] = currentw[lgth2 - 1];  // lgth2==0 no toki error

    if (tailgp)
        lasti = lgth1 + 1;
    else
        lasti = lgth1;
    lastj = lgth2 + 1;

#if XXXXXXX
    fprintf(stderr, "currentw = \n");
    for (i = 0; i < lgth1 + 1; i++) {
        fprintf(stderr, "%5.2f ", currentw[i]);
    }
    fprintf(stderr, "\n");
    fprintf(stderr, "initverticalw = \n");
    for (i = 0; i < lgth2 + 1; i++) {
        fprintf(stderr, "%5.2f ", initverticalw[i]);
    }
    fprintf(stderr, "\n");
#endif

    for (i = 1; i < lasti; i++) {
        wtmp = previousw;
        previousw = currentw;
        currentw = wtmp;

        previousw[0] = initverticalw[i - 1];

        match_calc_mtx(amino_dynamicmtx, currentw, seq1, seq2, i, lgth2);
#if XXXXXXX
        fprintf(stderr, "\n");
        fprintf(stderr, "i=%d\n", i);
        fprintf(stderr, "currentw = \n");
        for (j = 0; j < lgth2; j++) {
            fprintf(stderr, "%5.2f ", currentw[j]);
        }
        fprintf(stderr, "\n");
#endif
#if XXXXXXX
        fprintf(stderr, "\n");
        fprintf(stderr, "i=%d\n", i);
        fprintf(stderr, "currentw = \n");
        for (j = 0; j < lgth2; j++) {
            fprintf(stderr, "%5.2f ", currentw[j]);
        }
        fprintf(stderr, "\n");
#endif
        currentw[0] = initverticalw[i];

        mi = previousw[0];
        mpi = 0;

        ijppt = ijp[i] + 1;
        mjpt = m + 1;
        prept = previousw;
        curpt = currentw + 1;
        mpjpt = mp + 1;

        if (i < lgth1)
            fpenalty_ex_i = fpenalty_ex;
        else
            fpenalty_ex_i = 0.0;  // 2018/May/11

        for (j = 1; j < lastj; j++) {
            wm = *prept;
            *ijppt = 0;

#if 0
			fprintf( stderr, "%5.0f->", wm );
#endif
#if 0
			fprintf( stderr, "%5.0f?", g );
#endif
            if ((g = mi + fpenalty) > wm) {
                wm = g;
                *ijppt = -(j - mpi);
            }
            if ((g = *prept) >= mi)
            //			if( (g=*prept) > mi )
            {
                mi = g;
                mpi = j - 1;
            }
#if USE_PENALTY_EX
            mi += fpenalty_ex_i;
//			mi += fpenalty_ex;
#endif

            if ((g = *mjpt + fpenalty) > wm) {
                wm = g;
                *ijppt = +(i - *mpjpt);
            }
            if ((g = *prept) >= *mjpt)
            {
                *mjpt = g;
                *mpjpt = i - 1;
            }
#if USE_PENALTY_EX
            if (j < lgth2)
                m[j] += fpenalty_ex;
#endif

            *curpt++ += wm;
            ijppt++;
            mjpt++;
            prept++;
            mpjpt++;
        }
        lastverticalw[i] = currentw[lgth2 - 1];  // lgth2==0 no toki error
    }

    wmo = Atracking(opts, ctx, currentw, lastverticalw, seq1, seq2, mseq1, mseq2, ijp, tailgp, warpis, warpjs, warpbase);
    if (!tailgp)
        wm = wmo;

    //	reporterr( "wm (after tracking) = %f\n", wm );
    if (warpis)
        free(warpis);
    if (warpjs)
        free(warpjs);

    resultlen = strlen(mseq1[0]);
    aln_assert(!(alloclen < resultlen || resultlen > N));

    strcpy(seq1[0], mseq1[0]);
    strcpy(seq2[0], mseq2[0]);
#if 0
	fprintf( stderr, "\n" );
	fprintf( stderr, ">\n%s\n", mseq1[0] );
	fprintf( stderr, ">\n%s\n", mseq2[0] );
	fprintf( stderr, "wm = %f\n", wm );
#endif

    return (wm);
}

double
G__align11_noalign(Context* ctx, double** n_dynamicmtx, int penal, int penal_ex, char** seq1, char** seq2)
/* warp mitaiou */
{
    //	int k;
    register int i, j;
    int          lasti; /* outgap == 0 -> lgth1, outgap == 1 -> lgth1+1 */
    int          lgth1, lgth2;
    //	int resultlen;
    double  wm; /* int ?????? */
    double  g;
    double *currentw, *previousw;
    double  fpenalty = (double)penal;
#if USE_PENALTY_EX
    double fpenalty_ex = (double)penal_ex;
    double fpenalty_ex_i;
#endif
#if 1
    double* wtmp;
    double *mjpt, *prept, *curpt;
//	int *mpjpt;
#endif
    static double   mi, *m;
    static double * w1, *w2;
    static double*  match;
    static double*  initverticalw; /* kufuu sureba iranai */
    static double*  lastverticalw; /* kufuu sureba iranai */
    static int**    intwork;
    static double** doublework;
    static int      orlgth1 = 0, orlgth2 = 0;
    static double** amino_dynamicmtx;

    if (seq1 == NULL) {
        if (orlgth1 > 0 && orlgth2 > 0) {
            orlgth1 = 0;
            orlgth2 = 0;
            FreeFloatVec(w1);
            FreeFloatVec(w2);
            FreeFloatVec(match);
            FreeFloatVec(initverticalw);
            FreeFloatVec(lastverticalw);
            free(m);

            FreeFloatMtx(doublework);
            FreeIntMtx(intwork);
            FreeDoubleMtx(amino_dynamicmtx);
        }
        return (0.0);
    }

    wm = 0.0;

    lgth1 = strlen(seq1[0]);
    lgth2 = strlen(seq2[0]);

#if 0
	if( lgth1 <= 0 || lgth2 <= 0 )
	{
		fprintf( stderr, "WARNING (g11): lgth1=%d, lgth2=%d\n", lgth1, lgth2 );
	}
#endif

    if (lgth1 > orlgth1 || lgth2 > orlgth2) {
        int ll1, ll2;

        if (orlgth1 > 0 && orlgth2 > 0) {
            FreeFloatVec(w1);
            FreeFloatVec(w2);
            FreeFloatVec(match);
            FreeFloatVec(initverticalw);
            FreeFloatVec(lastverticalw);

            FreeFloatVec(m);

            FreeFloatMtx(doublework);
            FreeIntMtx(intwork);

            FreeDoubleMtx(amino_dynamicmtx);
        }

        ll1 = MAX((int)(1.3 * lgth1), orlgth1) + 100;
        ll2 = MAX((int)(1.3 * lgth2), orlgth2) + 100;

#if DEBUG
        fprintf(stderr, "\ntrying to allocate (%d+%d)xn matrices ... ", ll1, ll2);
#endif

        w1 = AllocateFloatVec(ll2 + 2);
        w2 = AllocateFloatVec(ll2 + 2);
        match = AllocateFloatVec(ll2 + 2);

        initverticalw = AllocateFloatVec(ll1 + 2);
        lastverticalw = AllocateFloatVec(ll1 + 2);

        m = AllocateFloatVec(ll2 + 2);

        doublework = AllocateFloatMtx(ctx->nalphabets, MAX(ll1, ll2) + 2);
        intwork = AllocateIntMtx(ctx->nalphabets, MAX(ll1, ll2) + 2);

        amino_dynamicmtx = AllocateDoubleMtx(0x100, 0x100);
#if DEBUG
        fprintf(stderr, "succeeded\n");
#endif

        orlgth1 = ll1 - 100;
        orlgth2 = ll2 - 100;
    }

    for (i = 0; i < ctx->nalphabets; i++)
        for (j = 0; j < ctx->nalphabets; j++)
            amino_dynamicmtx[(int)ctx->amino[i]][(int)ctx->amino[j]] = (double)n_dynamicmtx[i][j];

#if 0
	for( i=0; i<lgth1; i++ ) 
		fprintf( stderr, "ogcp1[%d]=%f\n", i, ogcp1[i] );
#endif

    currentw = w1;
    previousw = w2;

    match_calc_mtx(amino_dynamicmtx, initverticalw, seq2, seq1, 0, lgth1);

    match_calc_mtx(amino_dynamicmtx, currentw, seq1, seq2, 0, lgth2);

    if (1)  // tsuneni outgap-1
    {
        for (i = 1; i < lgth1 + 1; i++) {
            initverticalw[i] += fpenalty;
#if USE_PENALTY_EX  // 2018/Apr/23
//			initverticalw[i] += fpenalty_ex * i; // ato de fukkatsu
//			reporterr( "added _ex\n" );
#endif
        }
        for (j = 1; j < lgth2 + 1; j++) {
            currentw[j] += fpenalty;
#if USE_PENALTY_EX  // 2018/Apr/23
//			currentw[j] += fpenalty_ex * j; // ato de fukkatsu
//			reporterr( "added _ex\n" );
#endif
        }
    }

    for (j = 1; j < lgth2 + 1; ++j) {
        m[j] = currentw[j - 1];
    }

    if (lgth2 == 0)
        lastverticalw[0] = 0.0;  // lgth2==0 no toki error
    else
        lastverticalw[0] = currentw[lgth2 - 1];  // lgth2==0 no toki error

    if (1)
        lasti = lgth1 + 1;
    else
        lasti = lgth1;  // tsuneni outgap-1

#if XXXXXXX
    fprintf(stderr, "currentw = \n");
    for (i = 0; i < lgth1 + 1; i++) {
        fprintf(stderr, "%5.2f ", currentw[i]);
    }
    fprintf(stderr, "\n");
    fprintf(stderr, "initverticalw = \n");
    for (i = 0; i < lgth2 + 1; i++) {
        fprintf(stderr, "%5.2f ", initverticalw[i]);
    }
    fprintf(stderr, "\n");
#endif

    for (i = 1; i < lasti; i++) {
        wtmp = previousw;
        previousw = currentw;
        currentw = wtmp;

        previousw[0] = initverticalw[i - 1];

        match_calc_mtx(amino_dynamicmtx, currentw, seq1, seq2, i, lgth2);
#if XXXXXXX
        fprintf(stderr, "\n");
        fprintf(stderr, "i=%d\n", i);
        fprintf(stderr, "currentw = \n");
        for (j = 0; j < lgth2; j++) {
            fprintf(stderr, "%5.2f ", currentw[j]);
        }
        fprintf(stderr, "\n");
#endif
#if XXXXXXX
        fprintf(stderr, "\n");
        fprintf(stderr, "i=%d\n", i);
        fprintf(stderr, "currentw = \n");
        for (j = 0; j < lgth2; j++) {
            fprintf(stderr, "%5.2f ", currentw[j]);
        }
        fprintf(stderr, "\n");
#endif
        currentw[0] = initverticalw[i];

        mi = previousw[0];

        mjpt = m + 1;
        prept = previousw;
        curpt = currentw + 1;

        if (i < lgth1)
            fpenalty_ex_i = fpenalty_ex;
        else
            fpenalty_ex_i = 0.0;  // 2018/May/11

        for (j = 1; j < lgth2 + 1; j++) {
            wm = *prept;

#if 0
			fprintf( stderr, "%5.0f->", wm );
#endif
#if 0
			fprintf( stderr, "%5.0f?", g );
#endif
            if ((g = mi + fpenalty) > wm) {
                wm = g;
            }
            //			if( (g=*prept) >= mi )
            if ((g = *prept) > mi)  // onaji hazu
            {
                mi = g;
            }
#if USE_PENALTY_EX
            mi += fpenalty_ex_i;
#endif

#if 0 
			fprintf( stderr, "%5.0f?", g );
#endif
            if ((g = *mjpt + fpenalty) > wm) {
                wm = g;
            }
            //			if( (g=*prept) >= *mjpt )
            if ((g = *prept) > *mjpt)  // onaji hazu
            {
                *mjpt = g;
            }
#if USE_PENALTY_EX
            if (j < lgth2)  // 2018/May/11
                m[j] += fpenalty_ex;
#endif

#if 0
			fprintf( stderr, "%5.0f ", wm );
#endif
            *curpt++ += wm;
            mjpt++;
            prept++;
        }
        lastverticalw[i] = currentw[lgth2 - 1];  // lgth2==0 no toki error
    }

#if 0
	fprintf( stderr, "\n" );
	fprintf( stderr, ">\n%s\n", mseq1[0] );
	fprintf( stderr, ">\n%s\n", mseq2[0] );
	fprintf( stderr, "wm (noalign) = %f\n", wm );
#endif

    return (wm);
}
