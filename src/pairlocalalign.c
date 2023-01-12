#include "mltaln.h"

#define DEBUG 0
#define IODEBUG 0

#define NODIST -9999

static char* laraparams;
static char  foldalignopt[1000];
static int   stdout_align;
static int   stdout_dist;
static int   store_localhom;
static int   store_dist;
static int   laste;
static int   lastm;
static int   lastsubopt;
static int   lastonce;

typedef struct _lastres {
    int   score;
    int   start1;
    int   start2;
    char* aln1;
    char* aln2;
} Lastres;

typedef struct _reg {
    int start;
    int end;
} Reg;

typedef struct _aln {
    int  nreg;
    Reg* reg1;
    Reg* reg2;
} Aln;

typedef struct _lastresx {
    int  score;
    int  naln;
    Aln* aln;
} Lastresx;

static void
putlocalhom_last(Context* ctx, char* s1, char* s2, LocalHom* localhompt, Lastresx* lastresx) {
    char *    pt1, *pt2;
    int       naln, nreg;
    int       iscore;
    int       isumscore;
    int       sumoverlap;
    LocalHom* tmppt = localhompt;
    LocalHom* tmppt2;
    LocalHom* localhompt0;
    Reg *     rpt1, *rpt2;
    Aln*      apt;
    int       nlocalhom = 0;
    int       len;

    naln = lastresx->naln;
    apt = lastresx->aln;

    if (naln == 0)
        return;
    while (naln--) {
        rpt1 = apt->reg1;
        rpt2 = apt->reg2;
        nreg = apt->nreg;
        isumscore = 0;
        sumoverlap = 0;
        while (nreg--) {
            if (nlocalhom++ > 0) {
                tmppt->next = (LocalHom*)calloc(1, sizeof(LocalHom));
                tmppt = tmppt->next;
                tmppt->next = NULL;
            }
            tmppt->start1 = rpt1->start;
            tmppt->start2 = rpt2->start;
            tmppt->end1 = rpt1->end;
            tmppt->end2 = rpt2->end;
            tmppt->korh = 'h';
            if (rpt1 == apt->reg1)
                localhompt0 = tmppt;  // ?

            len = tmppt->end1 - tmppt->start1 + 1;

            pt1 = s1 + tmppt->start1;
            pt2 = s2 + tmppt->start2;
            iscore = 0;
            while (len--) {
                iscore += ctx->n_dis[(int)ctx->amino_n[(unsigned char)*pt1++]][(int)ctx->amino_n[(unsigned char)*pt2++]];  // - offset はいらないかも
            }

            if (ctx->divpairscore) {
                tmppt->overlapaa = tmppt->end2 - tmppt->start2 + 1;
                tmppt->opt = (double)iscore / tmppt->overlapaa * 5.8 / 600;
            } else {
                isumscore += iscore;
                sumoverlap += tmppt->end2 - tmppt->start2 + 1;
            }
            rpt1++;
            rpt2++;
        }

        if (!ctx->divpairscore) {
            for (tmppt2 = localhompt0; tmppt2; tmppt2 = tmppt2->next) {
                tmppt2->overlapaa = sumoverlap;
                tmppt2->opt = (double)isumscore * 5.8 / (600 * sumoverlap);
            }
        }
        apt++;
    }
}

int
countamino(char* s, int end) {
    int val = 0;
    while (end--)
        if (*s++ != '-')
            val++;
    return (val);
}

static double
score2dist(double pscore, double selfscore1, double selfscore2) {
    double val;
    double bunbo;
    //	fprintf( stderr, "In score2dist\n" );

    if ((bunbo = MIN(selfscore1, selfscore2)) == 0.0)
        val = 2.0;
    else if (bunbo < pscore)  // mondai ari
        val = 0.0;
    else
        val = (1.0 - pscore / bunbo) * 2.0;
    return (val);
}

static void
pairalign(aln_Opts opts, Context* ctx, const char* const* name, char** seq, char** aseq, char** dseq, int* thereisxineachseq, char** mseq1, char** mseq2, int alloclen, Lastresx** lastresx, double** distancemtx, LocalHom** localhomtable, int ngui) {
    int       i, j, ilim, jst, jj;
    int       off1, off2, thereisx;
    double    pscore = 0.0;
    FILE *    hat2p, *hat3p;
    double*   selfscore;
    double*   effarr1;
    double*   effarr2;
    char*     hat2file = "hat2";
    LocalHom* tmpptr;
    int       intdum;
    char***   bpp = NULL;
    char **   distseq1, **distseq2;
    int *     targetmap, *targetmapr;

    targetmap = calloc(ctx->njob, sizeof(int));
    targetmapr = calloc(ctx->njob, sizeof(int));
    for (int i = 0; i < ctx->njob; i++) {
        targetmap[i] = targetmapr[i] = i;
    }

    char** dumseq1 = 0;
    char** dumseq2 = 0;

    distseq1 = AllocateCharMtx(1, 0);  // muda
    distseq2 = AllocateCharMtx(1, 0);  // muda

    selfscore = AllocateDoubleVec(ctx->njob);
    effarr1 = AllocateDoubleVec(ctx->njob);
    effarr2 = AllocateDoubleVec(ctx->njob);

    for (int i = 0; i < ctx->njob; i++) {
        pscore = 0.0;
        for (char* pt = seq[i]; *pt; pt++) {
            pscore += ctx->amino_dis[(unsigned char)*pt][(unsigned char)*pt];
        }
        selfscore[i] = pscore;
    }

    ilim = ctx->njob - 1;
    for (int i = 0; i < ilim; i++) {
        jst = i + 1;
        for (j = jst; j < ctx->njob; j++) {

            strcpy(aseq[0], seq[i]);
            strcpy(aseq[1], seq[j]);

            effarr1[0] = 1.0;
            effarr2[0] = 1.0;
            mseq1[0] = aseq[0];
            mseq2[0] = aseq[1];

            thereisx = thereisxineachseq[i] + thereisxineachseq[j];
            distseq1[0] = dseq[i];
            distseq2[0] = dseq[j];

            if (opts.use_fft) {
                pscore = Falign(opts, ctx, NULL, NULL, ctx->n_dis_consweight_multi, mseq1, mseq2, effarr1, effarr2, NULL, NULL, 1, 1, alloclen, &intdum);
                off1 = off2 = 0;
            } else {
                if (ctx->usenaivescoreinsteadofalignmentscore) {
                    G__align11(ctx, ctx->n_dis_consweight_multi, mseq1, mseq2, alloclen, ctx->outgap, ctx->outgap);
                    pscore = (double)naivepairscore11(ctx, mseq1[0], mseq2[0], 0.0);  // uwagaki
                } else {
                    if (store_localhom && (targetmap[i] != -1 || targetmap[j] != -1)) {
                        pscore = G__align11(ctx, ctx->n_dis_consweight_multi, mseq1, mseq2, alloclen, ctx->outgap, ctx->outgap);
                        if (thereisx)
                            pscore = G__align11_noalign(ctx, ctx->n_dis_consweight_multi, ctx->penalty, ctx->penalty_ex, distseq1, distseq2);
                    } else
                        pscore = G__align11_noalign(ctx, ctx->n_dis_consweight_multi, ctx->penalty, ctx->penalty_ex, distseq1, distseq2);
                }
                off1 = off2 = 0;
                break;
            }

            if (opts.alg == 't' || (mseq1[0][0] != 0 && mseq2[0][0] != 0))  // 't' no jouken ha iranai to omou. if( ( mseq1[0][0] != 0 && mseq2[0][0] != 0  ) )
            {
                if ((ctx->nadd == 0 || (opts.alg != 'Y' && opts.alg != 'r') || (i < ctx->njob - ctx->nadd && ctx->njob - ctx->nadd <= j))) {
                    if (!store_localhom)
                        ;
                    else if (opts.alg == 'R')
                        putlocalhom_last(ctx, mseq1[0], mseq2[0], localhomtable[i] + j, lastresx[i] + j);
                    else if (opts.alg == 'r')
                        putlocalhom_last(ctx, mseq1[0], mseq2[0], localhomtable[i] + j - (ctx->njob - ctx->nadd), lastresx[i] + j - (ctx->njob - ctx->nadd));
                    else if (opts.alg == 'H')
                        putlocalhom_ext(ctx, mseq1[0], mseq2[0], localhomtable[i] + j, off1, off2, 'h');
                    else if (opts.alg == 'Y')
                        putlocalhom2(ctx, mseq1[0], mseq2[0], localhomtable[i] + j - (ctx->njob - ctx->nadd), off1, off2, 'h');
                    else if (opts.alg != 'S' && opts.alg != 'V')
                        putlocalhom2(ctx, mseq1[0], mseq2[0], localhomtable[i] + j - i, off1, off2, 'h');
                    else {
                        if (targetmap[i] != -1 && targetmap[j] != -1) {
                            putlocalhom2(ctx, mseq1[0], mseq2[0], localhomtable[targetmap[i]] + j, off1, off2, 'h');
                            putlocalhom2(ctx, mseq2[0], mseq1[0], localhomtable[targetmap[j]] + i, off2, off1, 'h');  // sukoshi muda.
                        } else if (targetmap[j] != -1)
                            putlocalhom2(ctx, mseq2[0], mseq1[0], localhomtable[targetmap[j]] + i, off2, off1, 'h');
                        else if (targetmap[i] != -1)
                            putlocalhom2(ctx, mseq1[0], mseq2[0], localhomtable[targetmap[i]] + j, off1, off2, 'h');
                        else {
                            reporterr("okashii\n");
                            exit(1);
                        }
                    }
                }

                pscore = score2dist(pscore, selfscore[i], selfscore[j]);
            } else {
                pscore = 2.0;
            }

            if (stdout_align) {
                if (opts.alg != 't') {
                    fprintf(stdout, "sequence %d - sequence %d, pairwise distance = %10.5f\n", i + 1, j + 1, pscore);
                    fprintf(stdout, ">%s\n", name[i]);
                    write1seq(stdout, mseq1[0]);
                    fprintf(stdout, ">%s\n", name[j]);
                    write1seq(stdout, mseq2[0]);
                    fprintf(stdout, "\n");
                }
            }
            if (stdout_dist)
                fprintf(stdout, "%d %d d=%.3f\n", i + 1, j + 1, pscore);
            if (store_dist) {
                if (opts.alg == 'Y' || opts.alg == 'r')
                    distancemtx[i][j - (ctx->njob - ctx->nadd)] = pscore;
                else
                    distancemtx[i][j - i] = pscore;
            }
        }
    }

    if (store_dist && ngui == 0) {
        hat2p = fopen(hat2file, "w");
        if (!hat2p)
            ErrorExit("Cannot open hat2.");
        if (opts.alg == 'Y' || opts.alg == 'r')
            WriteHat2_part_pointer(hat2p, ctx->njob, ctx->nadd, name, distancemtx);
        else
            WriteFloatHat2_pointer_halfmtx(ctx, hat2p, ctx->njob, name, distancemtx);  // jissiha double
        fclose(hat2p);
    }

    hat3p = fopen("hat3", "w");
    if (!hat3p)
        ErrorExit("Cannot open hat3.");
    if (store_localhom && ngui == 0) {
        fprintf(stderr, "\n\n##### writing hat3\n");
        if (opts.alg == 'Y' || opts.alg == 'r')
            ilim = ctx->njob - ctx->nadd;
        else
            ilim = ctx->njob - 1;
        for (i = 0; i < ilim; i++) {
            if (opts.alg == 'Y' || opts.alg == 'r') {
                jst = ctx->njob - ctx->nadd;
                jj = 0;
            } else {
                jst = i;
                jj = 0;
            }
            for (j = jst; j < ctx->njob; j++, jj++) {
                for (tmpptr = localhomtable[i] + jj; tmpptr; tmpptr = tmpptr->next) {
                    //					fprintf( stderr, "j=%d, jj=%d\n", j, jj );
                    if (tmpptr->opt == -1.0)
                        continue;
                    // tmptmptmptmptmp
                    //					if( alg == 'B' || alg == 'T' )
                    //						fprintf( hat3p, "%d %d %d %7.5f %d %d %d %d %p\n", i, j, tmpptr->overlapaa, 1.0, tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, (void *)tmpptr->next );
                    //					else
                    if (targetmap[j] == -1 || targetmap[i] < targetmap[j])
                        fprintf(hat3p, "%d %d %d %7.5f %d %d %d %d h\n", targetmapr[i], j, tmpptr->overlapaa, tmpptr->opt, tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2);
                    //						fprintf( hat3p, "%d %d %d %7.5f %d %d %d %d h\n", i, j, tmpptr->overlapaa, tmpptr->opt, tmpptr->start1, tmpptr->end1, tmpptr->start2+1, tmpptr->end2+1 ); // zettai dame!!!!
                }
            }
        }
//		if( ngui == 0 )
//		{
#if DEBUG
        fprintf(stderr, "calling FreeLocalHomTable\n");
#endif
        if (opts.alg == 'Y' || opts.alg == 'r')
            FreeLocalHomTable_part(localhomtable, (ctx->njob - ctx->nadd), ctx->nadd);
        else
            FreeLocalHomTable_half(localhomtable, ctx->njob);
#if DEBUG
        fprintf(stderr, "done. FreeLocalHomTable\n");
#endif
        //		}
    }
    fclose(hat3p);

    if (opts.alg == 's') {
        char** ptpt;
        for (i = 0; i < ctx->njob; i++) {
            ptpt = bpp[i];
            while (1) {
                if (*ptpt)
                    free(*ptpt);
                else
                    break;
                ptpt++;
            }
            free(bpp[i]);
        }
        free(bpp);
    }
    free(selfscore);
    free(effarr1);
    free(effarr2);
    if (opts.alg == 'N') {
        FreeCharMtx(dumseq1);
        FreeCharMtx(dumseq2);
    }
    free(distseq1);
    free(distseq2);
    if (store_dist && ngui == 0) {
        if (opts.alg == 'Y' || opts.alg == 'r')
            FreeDoubleMtx(distancemtx);  // 2020/Oct/23
        else
            FreeDoubleHalfMtx(distancemtx, ctx->njob);
    }

    free(targetmap);
    free(targetmapr);
}

int
pairlocalalign(
    aln_Opts           pairLocalAlignOpts,
    Context*           ctx,
    const char* const* name,
    char**             seq,
    double**           iscore,
    LocalHom**         localhomtable
) {
    aln_assert(ctx->njob >= 2);

    laste = 5000;
    lastm = 3;
    ctx->nadd = 0;
    lastsubopt = 0;
    lastonce = 0;
    foldalignopt[0] = 0;
    laraparams = NULL;
    ctx->fftkeika = 0;
    ctx->pslocal = -1000.0;
    ctx->fmodel = 0;
    ctx->fftscore = 1;
    ctx->fftRepeatStop = 0;
    ctx->fftNoAnchStop = 0;
    ctx->weight = 3;
    ctx->tbutree = 1;
    ctx->disp = 0;
    ctx->outgap = 1;
    ctx->mix = 0;
    ctx->tbitr = 0;
    ctx->tbweight = 0;
    ctx->tbrweight = 3;
    ctx->checkC = 0;
    ctx->kobetsubunkatsu = 0;
    ctx->divpairscore = 0;
    stdout_align = 0;
    stdout_dist = 0;
    store_dist = 1;
    store_localhom = 1;
    ctx->ppenalty_OP = NOTSPECIFIED;
    ctx->ppenalty_EX = NOTSPECIFIED;
    ctx->kimuraR = NOTSPECIFIED;
    ctx->pamN = NOTSPECIFIED;
    ctx->geta2 = GETA2;
    ctx->fftWinSize = NOTSPECIFIED;
    ctx->fftThreshold = NOTSPECIFIED;
    ctx->RNAppenalty = NOTSPECIFIED;
    ctx->RNApthr = NOTSPECIFIED;
    ctx->usenaivescoreinsteadofalignmentscore = 0;
    ctx->nwildcard = 0;

    constants(pairLocalAlignOpts, ctx, ctx->njob, seq);
    initSignalSM(ctx);
    initFiles(ctx);

    // TODO(sen) Sequence verification?

    int    alloclen = ctx->maxInputSeqLen * 2;
    char** bseq = AllocateCharMtx(ctx->njob, alloclen + 10);
    char** dseq = AllocateCharMtx(ctx->njob, alloclen + 10);
    int*   countsOfXs = AllocateIntVec(ctx->njob);
    for (int i = 0; i < ctx->njob; i++) {
        copyWithNoGaps(bseq[i], seq[i]);
        {
            char* d = dseq[i];
            char* m = bseq[i];
            int   countOfXs = 0;
            while (*m != 0) {
                if (*m == 'X' || *m == 'x') {
                    m++;
                    countOfXs++;
                } else {
                    *d++ = *m++;
                }
            }
            *d = 0;
            countsOfXs[i] = countOfXs;
        }
    }

    char** aseq = AllocateCharMtx(2, alloclen + 10);
    char** mseq1 = AllocateCharMtx(ctx->njob, 0);
    char** mseq2 = AllocateCharMtx(ctx->njob, 0);
    pairalign(
        pairLocalAlignOpts,
        ctx,
        name,
        bseq,
        aseq,
        dseq,
        countsOfXs,
        mseq1,
        mseq2,
        alloclen,
        0,
        iscore,
        localhomtable,
        ctx->njob
    );

    return 0;
}
