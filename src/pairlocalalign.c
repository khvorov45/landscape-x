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

static void
pairalign(aln_Opts opts, Context* ctx, const char* const* name, char** seq, char** aseq, char** dseq, int* thereisxineachseq, char** mseq1, char** mseq2, int alloclen, double** distancemtx) {
    int     j, ilim, jst;
    int     thereisx;
    double  pscore = 0.0;
    FILE*   hat2p;
    double* selfscore;
    double* effarr1;
    double* effarr2;
    char*   hat2file = "hat2";
    char ** distseq1, **distseq2;
    int *   targetmap, *targetmapr;

    targetmap = calloc(ctx->njob, sizeof(int));
    targetmapr = calloc(ctx->njob, sizeof(int));
    for (int i = 0; i < ctx->njob; i++) {
        targetmap[i] = targetmapr[i] = i;
    }

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

            if (store_localhom && (targetmap[i] != -1 || targetmap[j] != -1)) {
                pscore = G__align11(ctx, ctx->n_dis_consweight_multi, mseq1, mseq2, alloclen, ctx->outgap, ctx->outgap);
                if (thereisx) {
                    pscore = G__align11_noalign(ctx, ctx->n_dis_consweight_multi, ctx->penalty, ctx->penalty_ex, distseq1, distseq2);
                }
            } else {
                pscore = G__align11_noalign(ctx, ctx->n_dis_consweight_multi, ctx->penalty, ctx->penalty_ex, distseq1, distseq2);
            }
        }
    }

    if (store_dist && ctx->njob == 0) {
        hat2p = fopen(hat2file, "w");
        if (!hat2p)
            ErrorExit("Cannot open hat2.");
        if (opts.alg == 'Y' || opts.alg == 'r')
            WriteHat2_part_pointer(hat2p, ctx->njob, ctx->nadd, name, distancemtx);
        else
            WriteFloatHat2_pointer_halfmtx(ctx, hat2p, ctx->njob, name, distancemtx);  // jissiha double
        fclose(hat2p);
    }
}

int
pairlocalalign(
    aln_Opts           pairLocalAlignOpts,
    Context*           ctx,
    const char* const* name,
    char**             seq,
    double**           iscore
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
        iscore
    );

    return 0;
}
