#include "mltaln.h"

#define MACHIGAI 0
#define DEBUG 0
#define XXXXXXX 0
#define USE_PENALTY_EX 1
#define SLOW 0

#define TERMGAPFAC 0.0
#define TERMGAPFAC_EX 0.0

static double** impmtx = NULL;
static int      impalloclen = 0;

static void
imp_match_out_vead(double* imp, int i1, int lgth2) {
    double* pt = impmtx[i1];
    while (lgth2--)
        *imp++ += *pt++;
}
static void
imp_match_out_vead_tate(double* imp, int j1, int lgth1) {
    int i;
    for (i = 0; i < lgth1; i++)
        *imp++ += impmtx[i][j1];
}

void
imp_match_init_strict(aln_Context* ctx, int clus1, int clus2, int lgth1, int lgth2, char** seq1, char** seq2, double* eff1, double* eff2, double* eff1_kozo, double* eff2_kozo, aln_LocalHom*** localhom, char* swaplist, int* orinum1, int* orinum2) {
    if (seq1 == NULL) {
        if (impmtx)
            FreeFloatMtx(impmtx);
        impmtx = NULL;

        return;
    }

    if (impalloclen < lgth1 + 2 || impalloclen < lgth2 + 2) {
        if (impmtx)
            FreeFloatMtx(impmtx);
        impalloclen = MAX(lgth1, lgth2) + 2;
        impmtx = AllocateFloatMtx(impalloclen, impalloclen);
    }

    fillimp(ctx, impmtx, clus1, clus2, lgth1, lgth2, seq1, seq2, eff1, eff2, eff1_kozo, eff2_kozo, localhom, swaplist, orinum1, orinum2);  // uselh -> target -> localhomtable. seedinlh12 -> localhom ni haitteiru.
}

static void
match_calc(aln_Context* ctx, double** n_dynamicmtx, double* match, double** cpmx1, double** cpmx2, int i1, int lgth2, double** doublework, int** intwork, int initialize) {
    int      j, l;
    double** cpmxpd = doublework;
    int**    cpmxpdn = intwork;
    double * matchpt, *cpmxpdpt, **cpmxpdptpt;
    int *    cpmxpdnpt, **cpmxpdnptpt;
    double*  scarr;
    scarr = calloc(ctx->nalphabets, sizeof(double));
    if (initialize) {
        int count = 0;
        for (j = 0; j < lgth2; j++) {
            count = 0;
            for (l = 0; l < ctx->nalphabets; l++) {
                if (cpmx2[l][j]) {
                    cpmxpd[j][count] = cpmx2[l][j];
                    cpmxpdn[j][count] = l;
                    count++;
                }
            }
            cpmxpdn[j][count] = -1;
        }
    }

    {
        for (l = 0; l < ctx->nalphabets; l++) {
            scarr[l] = 0.0;
            for (j = 0; j < ctx->nalphabets; j++)
                scarr[l] += n_dynamicmtx[j][l] * cpmx1[j][i1];
        }
        matchpt = match;
        cpmxpdnptpt = cpmxpdn;
        cpmxpdptpt = cpmxpd;
        while (lgth2--) {
            *matchpt = 0.0;
            cpmxpdnpt = *cpmxpdnptpt++;
            cpmxpdpt = *cpmxpdptpt++;
            while (*cpmxpdnpt > -1)
                *matchpt += scarr[*cpmxpdnpt++] * *cpmxpdpt++;
            matchpt++;
        }
    }
    free(scarr);
}

static void
Atracking_localhom(double* impwmpt, char** seq1, char** seq2, char** mseq1, char** mseq2, int** ijp, int icyc, int jcyc, int* warpis, int* warpjs, int warpbase, int* ngap1, int* ngap2, int reuseprofiles, char** gt1, char** gt2) {
    int   i, j, l, iin, jin, ifi, jfi, lgth1, lgth2, k, limk;
    char *gaptable1, *gt1bk;
    char *gaptable2, *gt2bk;
    lgth1 = strlen(seq1[0]);
    lgth2 = strlen(seq2[0]);

    if (gt1 == NULL) {
        gt1bk = AllocateCharVec(lgth1 + lgth2 + 1);
        gt2bk = AllocateCharVec(lgth1 + lgth2 + 1);
    } else {
        gt1bk = *gt1;
        gt2bk = *gt2;
    }

    for (i = 0; i < lgth1 + 1; i++) {
        ijp[i][0] = i + 1;
    }
    for (j = 0; j < lgth2 + 1; j++) {
        ijp[0][j] = -(j + 1);
    }

    gaptable1 = gt1bk + lgth1 + lgth2;
    *gaptable1 = 0;
    gaptable2 = gt2bk + lgth1 + lgth2;
    *gaptable2 = 0;
    //*ngap1 = *ngap2 = 0; // shita de keisan

    iin = lgth1;
    jin = lgth2;
    limk = lgth1 + lgth2 + 1;
    *impwmpt = 0.0;
    for (k = 0; k < limk; k++) {
        if (ijp[iin][jin] >= warpbase) {
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
                *--gaptable1 = 'o';
                *--gaptable2 = '-';
                k++;
                //*ngap2 = 1; // shita de keisan
            }
            l = jin;
            while (--l >= 0) {
                *--gaptable1 = '-';
                *--gaptable2 = 'o';
                k++;
                //*ngap1 = 1;
            }
            break;
        } else {
            l = iin - ifi;
            while (--l) {
                *--gaptable1 = 'o';
                *--gaptable2 = '-';
                k++;
                //*ngap2 = 1;
            }
            l = jin - jfi;
            while (--l) {
                *--gaptable1 = '-';
                *--gaptable2 = 'o';
                k++;
                //*ngap1 = 1;
            }
        }
        if (iin == lgth1 || jin == lgth2)
            ;
        else {
            *impwmpt += (double)impmtx[iin][jin];

            //		fprintf( stderr, "impwm = %f (iin=%d, jin=%d) seq1=%c, seq2=%c\n", *impwmpt, iin, jin, seq1[0][iin], seq2[0][jin] );
        }
        if (iin <= 0 || jin <= 0)
            break;
        *--gaptable1 = 'o';
        *--gaptable2 = 'o';
        k++;
        iin = ifi;
        jin = jfi;
    }

    if (strchr(gaptable1, '-'))
        *ngap1 = 1;
    else
        *ngap1 = 0;
    if (strchr(gaptable2, '-'))
        *ngap2 = 1;
    else
        *ngap2 = 0;

    if (*ngap1 == 0 && reuseprofiles)
        ;
    else if (*ngap1 == 0) {
        limk = gt1bk + lgth1 + lgth2 - gaptable1;
        for (i = 0; i < icyc; i++) {
            strncpy(mseq1[i], seq1[i], limk);
            mseq1[i][limk] = 0;
        }
    } else {
        for (i = 0; i < icyc; i++)
            gapireru(mseq1[i], seq1[i], gaptable1);
    }

    if (*ngap2 == 0 && reuseprofiles)
        ;
    else if (*ngap2 == 0) {
        limk = gt2bk + lgth1 + lgth2 - gaptable2;  // nakutemo
        for (j = 0; j < jcyc; j++) {
            strncpy(mseq2[j], seq2[j], limk);
            mseq2[j][limk] = 0;
        }
    } else {
        for (j = 0; j < jcyc; j++)
            gapireru(mseq2[j], seq2[j], gaptable2);
    }

    if (gt1 == NULL) {
        free(gt1bk);
        free(gt2bk);
    } else {
        *gt1 = gaptable1;
        *gt2 = gaptable2;
    }
}

static void
createcpmxresult(aln_Context* ctx, double** cpmxresult, int limk, double eff1, double eff2, double*** cpmx1, double*** cpmx2, char* gaptable1, char* gaptable2, int usehist1, int usehist2) {
    int i, j, p;
    int alen = strlen(gaptable1);

#if 1  // sukoshi osoi
    for (i = 0; i < ctx->nalphabets; i++) {
        cpmxresult[i] = calloc(limk + 1, sizeof(double));
        for (j = 0; j < alen; j++)
            cpmxresult[i][j] = 0.0;
        for (j = 0, p = 0; j < alen; j++) {
            if (gaptable1[j] == '-')
                ;
            //if( amino_n['-'] == i ) cpmxresult[amino_n['-']][j] = eff1; // tsukawanai
            else
                cpmxresult[i][j] += (*cpmx1)[i][p++] * eff1;
        }

        for (j = 0, p = 0; j < alen; j++) {
            if (gaptable2[j] == '-')
                ;
            //if( amino_n['-'] == i ) cpmxresult[amino_n['-']][j] = eff2; // tsukawanai
            else
                cpmxresult[i][j] += (*cpmx2)[i][p++] * eff2;
        }

        if (usehist1) {
            free((*cpmx1)[i]);
            (*cpmx1)[i] = NULL;
        }
        if (usehist2) {
            free((*cpmx2)[i]);
            (*cpmx2)[i] = NULL;
        }
    }
#if 1
    if (usehist1) {
        free(*cpmx1);
        *cpmx1 = NULL;
    }
    if (usehist2) {
        free(*cpmx2);
        *cpmx2 = NULL;
    }
#endif
#endif
}

static void
creategapfreqresult(double** gapfresult, int limk, double eff1, double eff2, double* gapf1, double* gapf2, char* gaptable1, char* gaptable2)  // allocate nomi
{
    int j, p;
    int alen = strlen(gaptable1);
    (*gapfresult) = calloc(limk + 1, sizeof(double));

#if 1  // sukoshi osoi
    for (j = 0; j < alen + 1; j++)
        (*gapfresult)[j] = 0.0;
    for (j = 0, p = 0; j < alen + 1; j++) {
        if (gaptable1[j] == '-')
            ;
        else
            (*gapfresult)[j] += gapf1[p++] * eff1;
    }

    for (j = 0, p = 0; j < alen; j++) {
        if (gaptable2[j] == '-')
            ;
        else
            (*gapfresult)[j] += gapf2[p++] * eff2;
    }
    (*gapfresult)[j] = 1.0;

#endif
}

static void
createogresult(double** gapfresult, int limk, double eff1, double eff2, double* ori1, double* ori2, double* gf1, double* gf2, char* gaptable1, char* gaptable2)  // allocate nomi
{
    int j, p;
    int alen = strlen(gaptable1);

    *gapfresult = calloc(limk + 1, sizeof(double));
#if 1
    for (j = 0; j < alen; j++)
        (*gapfresult)[j] = 0.0;
    for (j = 0, p = 0; j < alen; j++) {
        if (gaptable1[j] == '-') {
            if (j == 0) {
                (*gapfresult)[j] += 1.0 * eff1;
            } else if (j && gaptable1[j - 1] != '-') {
                (*gapfresult)[j] += (gf1[p - 1]) * eff1;  // p-1 daijoubu?
            }
        } else {
            if (j == 0 || (j && gaptable1[j - 1] != '-')) {
                (*gapfresult)[j] += ori1[p] * eff1;
            }
            p++;
        }
    }

    for (j = 0, p = 0; j < alen; j++) {
        if (gaptable2[j] == '-') {
            if (j == 0) {
                (*gapfresult)[j] += 1.0 * eff2;
            } else if (j && gaptable2[j - 1] != '-') {
                (*gapfresult)[j] += (gf2[p - 1]) * eff2;  // p-1 daijoubu?
            }
        } else {
            if (j == 0 || (j && gaptable2[j - 1] != '-')) {
                (*gapfresult)[j] += ori2[p] * eff2;
            }
            p++;
        }
    }

#endif
}

static void
createfgresult(double** gapfresult, int limk, double eff1, double eff2, double* ori1, double* ori2, double* gf1, double* gf2, char* gaptable1, char* gaptable2)  // allocate nomi
{
    int j, p;
    int alen = strlen(gaptable1);
    (*gapfresult) = calloc(limk + 1, sizeof(double));

    for (j = 0; j < alen; j++)
        (*gapfresult)[j] = 0.0;
    for (j = 0, p = 0; j < alen; j++) {
        if (gaptable1[j] == '-') {
            if (j == alen - 1) {
                (*gapfresult)[j] += eff1;
            } else if (gaptable1[j + 1] != '-') {
                (*gapfresult)[j] += gf1[p] * eff1;
            }
        } else {
            if (gaptable1[j + 1] != '-')
                (*gapfresult)[j] += ori1[p] * eff1;
            p++;
        }
    }

    for (j = 0, p = 0; j < alen; j++) {
        if (gaptable2[j] == '-') {
            if (j == alen - 1) {
                (*gapfresult)[j] += eff2;
            } else if (gaptable2[j + 1] != '-') {
                (*gapfresult)[j] += gf2[p] * eff2;
            }
        } else {
            if (gaptable2[j + 1] != '-')
                (*gapfresult)[j] += ori2[p] * eff2;
            p++;
        }
    }
}

static double
Atracking(double* lasthorizontalw, double* lastverticalw, double fpenalty, double fpenalty_ex, char** seq1, char** seq2, char** mseq1, char** mseq2, int** ijp, int icyc, int jcyc, int tailgp, int* warpis, int* warpjs, int warpbase, int* ngap1, int* ngap2, int reuseprofiles, char** gt1, char** gt2) {
    int    i, j, l, iin, jin, ifi, jfi, lgth1, lgth2, k, limk;
    double wm;
    char * gaptable1, *gt1bk;
    char * gaptable2, *gt2bk;
    lgth1 = strlen(seq1[0]);
    lgth2 = strlen(seq2[0]);

    if (gt1 == NULL) {
        gt1bk = AllocateCharVec(lgth1 + lgth2 + 1);
        gt2bk = AllocateCharVec(lgth1 + lgth2 + 1);
        //		gaptable1 = gt1bk + lgth1+lgth2;
        //		gaptable2 = gt2bk + lgth1+lgth2;
    } else {
        //		gaptable1 = *gt1 + lgth1+lgth2;
        //		gaptable2 = *gt2 + lgth1+lgth2;
        gt1bk = *gt1;
        gt2bk = *gt2;
    }
    //	*gaptable1 = 0;
    //	*gaptable2 = 0;

#if 0
	for( i=0; i<lgth1; i++ ) 
	{
		fprintf( stderr, "lastverticalw[%d] = %f\n", i, lastverticalw[i] );
	}
#endif

    if (tailgp == 1)
        ;
    else {
#if 1
        double g;
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

    for (i = 0; i < lgth1 + 1; i++) {
        ijp[i][0] = i + 1;
    }
    for (j = 0; j < lgth2 + 1; j++) {
        ijp[0][j] = -(j + 1);
    }

    gaptable1 = gt1bk + lgth1 + lgth2;
    *gaptable1 = 0;
    gaptable2 = gt2bk + lgth1 + lgth2;
    *gaptable2 = 0;

    iin = lgth1;
    jin = lgth2;
    limk = lgth1 + lgth2 + 1;
    for (k = 0; k < limk; k++) {
        if (ijp[iin][jin] >= warpbase) {
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
                *--gaptable1 = 'o';
                *--gaptable2 = '-';
                k++;
            }
            l = jin;
            while (--l >= 0) {
                *--gaptable1 = '-';
                *--gaptable2 = 'o';
                k++;
            }
            break;
        } else {
            l = iin - ifi;
            while (--l) {
                *--gaptable1 = 'o';
                *--gaptable2 = '-';
                k++;
            }
            l = jin - jfi;
            while (--l) {
                *--gaptable1 = '-';
                *--gaptable2 = 'o';
                k++;
            }
        }
        if (iin <= 0 || jin <= 0)
            break;
        *--gaptable1 = 'o';
        *--gaptable2 = 'o';
        k++;
        iin = ifi;
        jin = jfi;
    }

    if (strchr(gaptable1, '-'))
        *ngap1 = 1;
    else
        *ngap1 = 0;
    if (strchr(gaptable2, '-'))
        *ngap2 = 1;
    else
        *ngap2 = 0;
    if (*ngap1 == 0 && reuseprofiles)
        ;
    else if (*ngap1 == 0) {
        limk = gt1bk + lgth1 + lgth2 - gaptable1;
        for (i = 0; i < icyc; i++) {
            strncpy(mseq1[i], seq1[i], limk);
            mseq1[i][limk] = 0;
        }
    } else {
        for (i = 0; i < icyc; i++)
            gapireru(mseq1[i], seq1[i], gaptable1);
    }

    if (*ngap2 == 0 && reuseprofiles)
        ;
    else if (*ngap2 == 0) {
        limk = gt2bk + lgth1 + lgth2 - gaptable2;  // nakutemo
        for (j = 0; j < jcyc; j++) {
            strncpy(mseq2[j], seq2[j], limk);
            mseq2[j][limk] = 0;
        }
    } else {
        for (j = 0; j < jcyc; j++)
            gapireru(mseq2[j], seq2[j], gaptable2);
    }

    if (gt1 == NULL) {
        free(gt1bk);
        free(gt2bk);
    } else {
        *gt1 = gaptable1;
        *gt2 = gaptable2;
    }

    return (wm);
}

double
A__align(aln_Context* ctx, double** n_dynamicmtx, int penalty_l, int penalty_ex_l, char** seq1, char** seq2, double* eff1, double* eff2, int icyc, int jcyc, int alloclen, int constraint, double* impmatch, char* sgap1, char* sgap2, char* egap1, char* egap2, int headgp, int tailgp, int firstmem, int calledbyfulltreebase, double*** cpmxchild0, double*** cpmxchild1, double*** cpmxresult, double orieff1, double orieff2) {
    int          reuseprofiles;
    static int   previousfirstlen;
    static int   previousicyc;
    static int   previousfirstmem;
    static int   previouscall;
    int          ngap1, ngap2;
    register int i, j;
    int          lasti, lastj;
    int          lgth1, lgth2;
    int          resultlen;
    double       wm = 0.0; /* int ?????? */
    double       wmo = 0.0;
    double       g;
    double *     currentw, *previousw;
#if USE_PENALTY_EX
    double fpenalty_ex = (double)penalty_ex_l;
#endif
#if 1
    double* wtmp;
    int*    ijppt;
    double *mjpt, *prept, *curpt;
    int*    mpjpt;
#endif
    static double   mi, *m;
    static int**    ijp;
    static int      mpi, *mp;
    static double * w1, *w2;
    static double*  match;
    static double*  initverticalw; /* kufuu sureba iranai */
    static double*  lastverticalw; /* kufuu sureba iranai */
    char**          mseq1;
    char**          mseq2;
    char**          mseq;
    static double * ogcp1, *ogcp1o;
    static double * ogcp2, *ogcp2o;
    static double * fgcp1, *fgcp1o;
    static double * fgcp2, *fgcp2o;
    double *        ogcp1opt, *ogcp2opt, *fgcp1opt, *fgcp2opt;
    static double** cpmx1;
    double***       cpmx1pt = NULL;
    static double** cpmx2;
    double***       cpmx2pt = NULL;
    static int**    intwork;
    static double** doublework;
    static int      orlgth1 = 0, orlgth2 = 0;
    static double*  gapfreq1;
    double*         gapfreq1pt;
    static double*  gapfreq2;
    double*         gapfreq2pt;
    double          fpenalty = (double)penalty_l;
    double*         fgcp2pt;
    double*         ogcp2pt;
    double          fgcp1va;
    double          ogcp1va;
    double*         gf2pt;
    double*         gf2ptpre;
    double          gf1va;
    double          gf1vapre;
    double          headgapfreq1;
    double          headgapfreq2;

    int* warpis = NULL;
    int* warpjs = NULL;

    int   warpbase;
    char *gt1, *gt2, *gt1bk, *gt2bk;

    if (seq1 == NULL) {
        if (orlgth1) {
            orlgth1 = 0;
            orlgth2 = 0;

            imp_match_init_strict(ctx, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

            FreeFloatVec(w1);
            FreeFloatVec(w2);
            FreeFloatVec(match);
            FreeFloatVec(initverticalw);
            FreeFloatVec(lastverticalw);

            FreeFloatVec(m);
            FreeIntVec(mp);

            FreeFloatVec(ogcp1);
            FreeFloatVec(ogcp1o);
            FreeFloatVec(ogcp2);
            FreeFloatVec(ogcp2o);
            FreeFloatVec(fgcp1);
            FreeFloatVec(fgcp1o);
            FreeFloatVec(fgcp2);
            FreeFloatVec(fgcp2o);

            FreeFloatMtx(cpmx1);
            FreeFloatMtx(cpmx2);

            FreeFloatVec(gapfreq1);
            FreeFloatVec(gapfreq2);

            FreeFloatMtx(doublework);
            FreeIntMtx(intwork);
        }
        return (0.0);
    }

    lgth1 = strlen(seq1[0]);
    lgth2 = strlen(seq2[0]);

    if (lgth1 == 0 && lgth2 == 0)
        return (0.0);

    if (lgth1 == 0) {
        for (i = 0; i < icyc; i++) {
            j = lgth2;
            seq1[i][j] = 0;
            while (j)
                seq1[i][--j] = '-';
        }
        return (0.0);
    }

    if (lgth2 == 0) {
        for (i = 0; i < jcyc; i++) {
            j = lgth1;
            seq2[i][j] = 0;
            while (j)
                seq2[i][--j] = '-';
        }
        return (0.0);
    }

    warpbase = lgth1 + lgth2;
    warpis = NULL;
    warpjs = NULL;

    mseq1 = AllocateCharMtx(icyc, 0);
    mseq2 = AllocateCharMtx(jcyc, 0);
    mseq = AllocateCharMtx(icyc + jcyc, lgth1 + lgth2 + 100);

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

            FreeFloatVec(ogcp1);
            FreeFloatVec(ogcp1o);
            FreeFloatVec(ogcp2);
            FreeFloatVec(ogcp2o);
            FreeFloatVec(fgcp1);
            FreeFloatVec(fgcp1o);
            FreeFloatVec(fgcp2);
            FreeFloatVec(fgcp2o);

            FreeFloatMtx(cpmx1);
            FreeFloatMtx(cpmx2);

            FreeFloatVec(gapfreq1);
            FreeFloatVec(gapfreq2);

            FreeFloatMtx(doublework);
            FreeIntMtx(intwork);
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

        ogcp1 = AllocateFloatVec(ll1 + 2);
        ogcp1o = AllocateFloatVec(ll1 + 2);
        ogcp2 = AllocateFloatVec(ll2 + 2);
        ogcp2o = AllocateFloatVec(ll2 + 2);
        fgcp1 = AllocateFloatVec(ll1 + 2);
        fgcp1o = AllocateFloatVec(ll1 + 2);
        fgcp2 = AllocateFloatVec(ll2 + 2);
        fgcp2o = AllocateFloatVec(ll2 + 2);

        cpmx1 = AllocateFloatMtx(ctx->nalphabets, ll1 + 2);
        cpmx2 = AllocateFloatMtx(ctx->nalphabets, ll2 + 2);
        previousfirstlen = -1;
        previousicyc = -1;

        gapfreq1 = AllocateFloatVec(ll1 + 2);
        gapfreq2 = AllocateFloatVec(ll2 + 2);

        doublework = AllocateFloatMtx(MAX(ll1, ll2) + 2, ctx->nalphabets);
        intwork = AllocateIntMtx(MAX(ll1, ll2) + 2, ctx->nalphabets + 1);

#if DEBUG
        fprintf(stderr, "succeeded\n");
#endif

        orlgth1 = ll1 - 100;
        orlgth2 = ll2 - 100;
    }

    for (i = 0; i < icyc; i++) {
        mseq1[i] = mseq[i];
        seq1[i][lgth1] = 0;
    }
    for (j = 0; j < jcyc; j++) {
        mseq2[j] = mseq[icyc + j];
        seq2[j][lgth2] = 0;
    }

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

    if (calledbyfulltreebase == 1 && previouscall && firstmem >= 0 && firstmem == previousfirstmem && lgth1 == previousfirstlen && icyc == previousicyc + 1)
        reuseprofiles = 1;
    else
        reuseprofiles = 0;

    if (ctx->n_dis[0][ctx->amino_n['-']] != 0) {
        exit(1);
    }

    if (cpmxresult) {
        if (sgap1) {
            exit(1);
        }

        if (cpmxchild0 && *cpmxchild0) {
            cpmx1pt = (cpmxchild0);
            gapfreq1pt = (*cpmxchild0)[ctx->nalphabets];
            ogcp1opt = (*cpmxchild0)[ctx->nalphabets + 1];
            fgcp1opt = (*cpmxchild0)[ctx->nalphabets + 2];
        } else {
            cpmx1pt = &cpmx1;
            cpmx_calc_new(ctx, seq1, *cpmx1pt, eff1, lgth1, icyc);

            gapfreq1pt = gapfreq1;
            gapcountf(gapfreq1pt, seq1, icyc, eff1, lgth1);
            for (i = 0; i < lgth1; i++)
                gapfreq1pt[i] = 1.0 - gapfreq1pt[i];

            ogcp1opt = ogcp1o;
            fgcp1opt = fgcp1o;
            st_OpeningGapCount(ogcp1opt, icyc, seq1, eff1, lgth1);
            st_FinalGapCount(fgcp1opt, icyc, seq1, eff1, lgth1);
        }

        if (cpmxchild1 && *cpmxchild1) {
            cpmx2pt = (cpmxchild1);
            gapfreq2pt = (*cpmxchild1)[ctx->nalphabets];
            ogcp2opt = (*cpmxchild1)[ctx->nalphabets + 1];
            fgcp2opt = (*cpmxchild1)[ctx->nalphabets + 2];
        } else {
            cpmx2pt = &cpmx2;
            cpmx_calc_new(ctx, seq2, *cpmx2pt, eff2, lgth2, jcyc);

            gapfreq2pt = gapfreq2;
            gapcountf(gapfreq2pt, seq2, jcyc, eff2, lgth2);
            for (i = 0; i < lgth2; i++)
                gapfreq2pt[i] = 1.0 - gapfreq2pt[i];

            ogcp2opt = ogcp2o;
            fgcp2opt = fgcp2o;
            st_OpeningGapCount(ogcp2opt, jcyc, seq2, eff2, lgth2);
            st_FinalGapCount(fgcp2opt, jcyc, seq2, eff2, lgth2);
        }

        //sgap1 == 1 ni taiou shite iani node, legacygap niyorazu, headgap to tailgap ha nai to minasu.

        headgapfreq1 = 1.0;
        headgapfreq2 = 1.0;
        gapfreq1pt[lgth1] = 1.0;
        gapfreq2pt[lgth2] = 1.0;
    } else {
        cpmx1pt = &cpmx1;
        cpmx2pt = &cpmx2;
        gapfreq1pt = gapfreq1;
        gapfreq2pt = gapfreq2;
        ogcp1opt = ogcp1o;
        ogcp2opt = ogcp2o;
        fgcp1opt = fgcp1o;
        fgcp2opt = fgcp2o;

        if (reuseprofiles) {
            cpmx_calc_add(ctx, seq1, *cpmx1pt, eff1, lgth1, icyc);
        } else {
            cpmx_calc_new(ctx, seq1, *cpmx1pt, eff1, lgth1, icyc);
        }
        cpmx_calc_new(ctx, seq2, *cpmx2pt, eff2, lgth2, jcyc);

        if (sgap1) {
            new_OpeningGapCount(ogcp1opt, icyc, seq1, eff1, lgth1, sgap1);
            new_FinalGapCount(fgcp1opt, icyc, seq1, eff1, lgth1, egap1);

            new_OpeningGapCount(ogcp2opt, jcyc, seq2, eff2, lgth2, sgap2);
            new_FinalGapCount(fgcp2opt, jcyc, seq2, eff2, lgth2, egap2);

            outgapcount(&headgapfreq1, icyc, sgap1, eff1);
            outgapcount(&headgapfreq2, jcyc, sgap2, eff2);
            outgapcount(gapfreq1pt + lgth1, icyc, egap1, eff1);
            outgapcount(gapfreq2pt + lgth2, jcyc, egap2, eff2);
        } else {
            if (reuseprofiles) {
                st_OpeningGapAdd(ogcp1opt, icyc, seq1, eff1, lgth1);
                st_FinalGapAdd(fgcp1opt, icyc, seq1, eff1, lgth1);
            } else {
                st_OpeningGapCount(ogcp1opt, icyc, seq1, eff1, lgth1);
                st_FinalGapCount(fgcp1opt, icyc, seq1, eff1, lgth1);
            }

            st_OpeningGapCount(ogcp2opt, jcyc, seq2, eff2, lgth2);
            st_FinalGapCount(fgcp2opt, jcyc, seq2, eff2, lgth2);

            headgapfreq1 = 0.0;
            headgapfreq2 = 0.0;
            gapfreq1pt[lgth1] = 0.0;
            gapfreq2pt[lgth2] = 0.0;
        }

        {
            if (reuseprofiles)
                gapcountadd(gapfreq1pt, seq1, icyc, eff1, lgth1);
            else
                gapcountf(gapfreq1pt, seq1, icyc, eff1, lgth1);
            gapcountf(gapfreq2pt, seq2, jcyc, eff2, lgth2);
            for (i = 0; i < lgth1 + 1; i++)
                gapfreq1pt[i] = 1.0 - gapfreq1pt[i];
            for (i = 0; i < lgth2 + 1; i++)
                gapfreq2pt[i] = 1.0 - gapfreq2pt[i];
            headgapfreq1 = 1.0 - headgapfreq1;
            headgapfreq2 = 1.0 - headgapfreq2;
        }
    }

    for (i = 0; i < lgth1; i++) {
        ogcp1[i] = 0.5 * (1.0 - ogcp1opt[i]) * fpenalty * (gapfreq1pt[i]);
        fgcp1[i] = 0.5 * (1.0 - fgcp1opt[i]) * fpenalty * (gapfreq1pt[i]);
    }

    for (i = 0; i < lgth2; i++) {
        ogcp2[i] = 0.5 * (1.0 - ogcp2opt[i]) * fpenalty * (gapfreq2pt[i]);
        fgcp2[i] = 0.5 * (1.0 - fgcp2opt[i]) * fpenalty * (gapfreq2pt[i]);
    }
#if 0
	for( i=0; i<lgth1; i++ ) 
		fprintf( stderr, "ogcp1[%d]=%f\n", i, ogcp1[i] );
#endif

    currentw = w1;
    previousw = w2;

    match_calc(ctx, n_dynamicmtx, initverticalw, *cpmx2pt, *cpmx1pt, 0, lgth1, doublework, intwork, 1);
    if (constraint)
        imp_match_out_vead_tate(initverticalw, 0, lgth1);  // 060306

    match_calc(ctx, n_dynamicmtx, currentw, *cpmx1pt, *cpmx2pt, 0, lgth2, doublework, intwork, 1);
    if (constraint)
        imp_match_out_vead(currentw, 0, lgth2);  // 060306

    if (headgp == 1) {
        for (i = 1; i < lgth1 + 1; i++) {
            initverticalw[i] += (ogcp1[0] * headgapfreq2 + fgcp1[i - 1] * gapfreq2pt[0]);
#if USE_PENALTY_EX  // 2018/Apr/23
            initverticalw[i] += fpenalty_ex * i;
#endif
        }
        for (j = 1; j < lgth2 + 1; j++) {
            currentw[j] += (ogcp2[0] * headgapfreq1 + fgcp2[j - 1] * gapfreq1pt[0]);
#if USE_PENALTY_EX  // 2018/Apr/23
            currentw[j] += fpenalty_ex * j;
#endif
        }
    } else {
        for (i = 1; i < lgth1 + 1; i++) {
            initverticalw[i] += fpenalty * TERMGAPFAC;  // motto tanjun de yoi?
#if USE_PENALTY_EX  // 2018/Apr/22
            initverticalw[i] += fpenalty_ex * i * TERMGAPFAC_EX;
#endif
        }
        for (j = 1; j < lgth2 + 1; j++) {
            currentw[j] += fpenalty * TERMGAPFAC;  // motto tanjun de yoi?
#if USE_PENALTY_EX  // 2018/Apr/22
            currentw[j] += fpenalty_ex * j * TERMGAPFAC_EX;
#endif
        }
    }

    for (j = 1; j < lgth2 + 1; ++j) {
        //		m[j] = currentw[j-1] + ogcp1[1]; mp[j] = 0;
        m[j] = currentw[j - 1] + ogcp1[1] * gapfreq2pt[j - 1];
        mp[j] = 0;
        ;
    }
    if (lgth2 == 0)
        lastverticalw[0] = 0.0;  // Falign kara yobaretatoki kounarukanousei ari
    else
        lastverticalw[0] = currentw[lgth2 - 1];

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
    fprintf(stderr, "fcgp\n");
    for (i = 0; i < lgth1; i++)
        fprintf(stderr, "fgcp1[%d]=%f\n", i, ogcp1[i]);
    for (i = 0; i < lgth2; i++)
        fprintf(stderr, "fgcp2[%d]=%f\n", i, ogcp2[i]);
#endif

    for (i = 1; i < lasti; i++) {
#ifdef enablemultithread
        //		fprintf( stderr, "chudan = %d, %d\n", *chudanpt, chudanref );
        if (chudanpt && *chudanpt != chudanref) {
            //			fprintf( stderr, "\n\n## CHUUDAN!!! S\n" );
            *chudanres = 1;
            free(mseq1);  // 2021/Sep
            free(mseq2);  // 2021/Sep
            FreeCharMtx(mseq);  // 2021/Sep
            return (-1.0);
        }
#endif
        wtmp = previousw;
        previousw = currentw;
        currentw = wtmp;

        previousw[0] = initverticalw[i - 1];

        match_calc(ctx, n_dynamicmtx, currentw, *cpmx1pt, *cpmx2pt, i, lgth2, doublework, intwork, 0);
#if XXXXXXX
        fprintf(stderr, "\n");
        fprintf(stderr, "i=%d\n", i);
        fprintf(stderr, "currentw = \n");
        for (j = 0; j < lgth2; j++) {
            fprintf(stderr, "%5.2f ", currentw[j]);
        }
        fprintf(stderr, "\n");
#endif
        if (constraint) {
            imp_match_out_vead(currentw, i, lgth2);
        }
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

#if 0
		fprintf( stderr, "%c ", seq1[0][i] );
		for( j=0; j<lgth2+1; j++ )
		{
			fprintf( stderr, "%5.0f ", currentw[j] );
		}
		fprintf( stderr, "\n"  );
#endif

        //		mi = previousw[0] + ogcp2[1]; mpi = 0;
        mi = previousw[0] + ogcp2[1] * gapfreq1pt[i - 1];
        mpi = 0;
        ijppt = ijp[i] + 1;
        mjpt = m + 1;
        prept = previousw;
        curpt = currentw + 1;
        mpjpt = mp + 1;
        fgcp2pt = fgcp2;
        ogcp2pt = ogcp2 + 1;
        fgcp1va = fgcp1[i - 1];
        ogcp1va = ogcp1[i];
        gf1va = gapfreq1pt[i];
        gf1vapre = gapfreq1pt[i - 1];
        gf2pt = gapfreq2pt + 1;
        gf2ptpre = gapfreq2pt;

        for (j = 1; j < lastj; j++) {
#ifdef xxxenablemultithread
            //			fprintf( stderr, "chudan = %d, %d\n", *chudanpt, chudanref );
            if (chudanpt && *chudanpt != chudanref) {
                //				fprintf( stderr, "\n\n## CHUUDAN!!! S\n" );
                free(mseq1);  // 2021/Sep
                free(mseq2);  // 2021/Sep
                FreeCharMtx(mseq);  // 2021/Sep
                *chudanres = 1;
                return (-1.0);
            }
#endif
            wm = *prept;
            *ijppt = 0;

#if 0
			fprintf( stderr, "\n i=%d, j=%d %c, %c", i, j, seq1[0][i], seq2[0][j] );
			fprintf( stderr, "%5.0f->", wm );
			fprintf( stderr, "%5.0f? (penal=%5.2f)", g=mi+*fgcp2pt*(1.0-gapfreq1pt[i]), *fgcp2pt*(1.0-gapfreq1pt[i]) );
#endif
            if ((g = mi + *fgcp2pt * gf1va) > wm) {
                wm = g;
                *ijppt = -(j - mpi);
                //				fprintf( stderr, "Jump to %d (%c)!", mpi, seq2[0][mpi] );
            }
            if ((g = *prept + *ogcp2pt * gf1vapre) >= mi)
            //			if( (g=*prept+*ogcp2pt*gf1vapre) > mi )
            {
                mi = g;
                mpi = j - 1;
            }
#if USE_PENALTY_EX
            mi += fpenalty_ex;
#endif

#if 0 
			fprintf( stderr, "%5.0f->", wm );
			fprintf( stderr, "%5.0f? (penal=%5.2f)", g=*mjpt+fgcp1va*(1.0-gapfreq2pt[j]), fgcp1va*(1.0-gapfreq2pt[j]) );
#endif
            if ((g = *mjpt + fgcp1va * *gf2pt) > wm) {
                wm = g;
                *ijppt = +(i - *mpjpt);
                //				fprintf( stderr, "Jump to %d (%c)!", *mpjpt, seq1[0][*mpjpt] );
            }
            if ((g = *prept + ogcp1va * *gf2ptpre) >= *mjpt)
            //			if( (g=*prept+ ogcp1va* *gf2ptpre) > *mjpt )
            {
                *mjpt = g;
                *mpjpt = i - 1;
            }
#if USE_PENALTY_EX
            m[j] += fpenalty_ex;
#endif

            *curpt++ += wm;
            ijppt++;
            mjpt++;
            prept++;
            mpjpt++;
            fgcp2pt++;
            ogcp2pt++;
            gf2ptpre++;
            gf2pt++;
        }
        lastverticalw[i] = currentw[lgth2 - 1];
    }

    gt1 = gt1bk = AllocateCharVec(lgth1 + lgth2 + 1);
    gt2 = gt2bk = AllocateCharVec(lgth1 + lgth2 + 1);

    if (constraint) {
        Atracking_localhom(impmatch, seq1, seq2, mseq1, mseq2, ijp, icyc, jcyc, warpis, warpjs, warpbase, &ngap1, &ngap2, reuseprofiles, &gt1, &gt2);
    } else {
        wmo = Atracking(currentw, lastverticalw, fpenalty, fpenalty_ex, seq1, seq2, mseq1, mseq2, ijp, icyc, jcyc, tailgp, warpis, warpjs, warpbase, &ngap1, &ngap2, reuseprofiles, &gt1, &gt2);
        if (!tailgp)
            wm = wmo;
    }

#if 1
    if (cpmxresult) {
        if (icyc + jcyc > 20)
        {
            int limk = gt1bk + lgth1 + lgth2 - gt1;
            limk = (int)(((limk + 1000) / 1000) * 1000);  // kiriage

            double totaleff1 = 0.0;
            double totaleff2 = 0.0;
            for (i = 0; i < icyc; i++)
                totaleff1 += eff1[i];
            for (j = 0; j < jcyc; j++)
                totaleff2 += eff2[j];

            if (fabs(totaleff1 - 1.0) > 0.001 || fabs(totaleff2 - 1.0) > 0.001) {
                exit(1);
            }
            totaleff1 = totaleff1 * orieff1 / (orieff1 + orieff2);
            totaleff2 = totaleff2 * orieff2 / (orieff1 + orieff2);

            *cpmxresult = AllocateDoubleMtx(ctx->nalphabets + 3, 0);  // gapcount, opg, fng no bun
            createcpmxresult(ctx, *cpmxresult, limk, totaleff1, totaleff2, cpmx1pt, cpmx2pt, gt1, gt2, (cpmx1 != *cpmx1pt), (cpmx2 != *cpmx2pt));  // naka de free
            creategapfreqresult(*cpmxresult + ctx->nalphabets, limk, totaleff1, totaleff2, gapfreq1pt, gapfreq2pt, gt1, gt2);  // naka deha free shinai
            // gapfreq1, gapfreq2 ha mada tsukau
            createogresult(*cpmxresult + ctx->nalphabets + 1, limk, totaleff1, totaleff2, ogcp1opt, ogcp2opt, gapfreq1pt, gapfreq2pt, gt1, gt2);  // naka deha free shinai
            if (cpmx1 != *cpmx1pt)
                free(ogcp1opt);
            if (cpmx2 != *cpmx2pt)
                free(ogcp2opt);
            createfgresult(*cpmxresult + ctx->nalphabets + 2, limk, totaleff1, totaleff2, fgcp1opt, fgcp2opt, gapfreq1pt, gapfreq2pt, gt1, gt2);  // naka deha free shinai
            if (cpmx1 != *cpmx1pt)
                free(fgcp1opt);
            if (cpmx2 != *cpmx2pt)
                free(fgcp2opt);

            if (cpmx1 != *cpmx1pt)
                free(gapfreq1pt);
            if (cpmx2 != *cpmx2pt)
                free(gapfreq2pt);
        } else
            *cpmxresult = NULL;
    }
#endif

    free(gt1bk);
    free(gt2bk);

    if (warpis)
        free(warpis);
    if (warpjs)
        free(warpjs);

    resultlen = strlen(mseq1[0]);
    aln_assert(!(alloclen < resultlen || resultlen > N));

    if (ngap1 || !reuseprofiles)
        for (i = 0; i < icyc; i++)
            strcpy(seq1[i], mseq1[i]);
    if (ngap2 || !reuseprofiles)
        for (j = 0; j < jcyc; j++)
            strcpy(seq2[j], mseq2[j]);

#if 0
	fprintf( stderr, "\n" );
	for( i=0; i<icyc; i++ ) fprintf( stderr, "%s\n", mseq1[i] );
	fprintf( stderr, "#####\n" );
	for( j=0; j<jcyc; j++ ) fprintf( stderr, "%s\n", mseq2[j] );
//exit( 1 );
#endif

    //	fprintf( stdout, "firstmem=%d, icyc=%d, jcyc=%d, wm = %f\n", firstmem, icyc, jcyc, wm );

    //	fprintf( stderr, "lgth1 = %d\n", lgth1 );
    //	fprintf( stderr, "->      %d\n", strlen( seq1[0] ) );
    previousfirstlen = lgth1;
    previousfirstmem = firstmem;
    previousicyc = icyc;
    previouscall = calledbyfulltreebase;

    free(mseq1);
    free(mseq2);
    FreeCharMtx(mseq);

    return (wm);
}
