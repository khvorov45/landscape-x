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

int
pairlocalalign(aln_Opts pairLocalAlignOpts, Context* ctx, char** seq) {
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
    {
        int* targetmap = calloc(ctx->njob, sizeof(int));
        int* targetmapr = calloc(ctx->njob, sizeof(int));
        for (int i = 0; i < ctx->njob; i++) {
            targetmap[i] = targetmapr[i] = i;
        }

        char** distseq1 = AllocateCharMtx(1, 0);  // muda
        char** distseq2 = AllocateCharMtx(1, 0);  // muda

        double* selfscore = AllocateDoubleVec(ctx->njob);
        double* effarr1 = AllocateDoubleVec(ctx->njob);
        double* effarr2 = AllocateDoubleVec(ctx->njob);

        double pscore = 0.0;
        for (int i = 0; i < ctx->njob; i++) {
            pscore = 0.0;
            for (char* pt = bseq[i]; *pt; pt++) {
                pscore += ctx->amino_dis[(unsigned char)*pt][(unsigned char)*pt];
            }
            selfscore[i] = pscore;
        }

        int ilim = ctx->njob - 1;
        for (int i = 0; i < ilim; i++) {
            for (int j = i + 1; j < ctx->njob; j++) {
                strcpy(aseq[0], bseq[i]);
                strcpy(aseq[1], bseq[j]);

                effarr1[0] = 1.0;
                effarr2[0] = 1.0;
                mseq1[0] = aseq[0];
                mseq2[0] = aseq[1];

                distseq1[0] = dseq[i];
                distseq2[0] = dseq[j];

                pscore = G__align11(ctx, ctx->n_dis_consweight_multi, mseq1, mseq2, alloclen, ctx->outgap, ctx->outgap);
                int thereisx = countsOfXs[i] + countsOfXs[j];
                if (thereisx) {
                    pscore = G__align11_noalign(ctx, ctx->n_dis_consweight_multi, ctx->penalty, ctx->penalty_ex, distseq1, distseq2);
                }
            }
        }
    }

    return 0;
}
