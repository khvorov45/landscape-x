#include <stdint.h>
#include <stdbool.h>

#define aln_IMPLEMENTATION
#include "mltaln.h"

aln_AlignResult
tbfast_main(aln_Str* strings, intptr_t stringsCount, void* out, intptr_t outBytes, aln_Opts opts) {
    aln_Arena permArena_ = {.base = out, .size = outBytes / 4};
    aln_Arena tempArena_ = {.base = (uint8_t*)out + permArena_.size, .size = outBytes - permArena_.size};

    aln_Arena* permArena = &permArena_;
    aln_Arena* tempArena = &tempArena_;

    Context* ctx = aln_arenaAllocStruct(tempArena, Context);

    ctx->RNAscoremtx = 'n';
    ctx->parallelizationstrategy = BAATARI1;
    ctx->newgapstr = "-";
    ctx->nalphabets = 26;
    ctx->ndistclass = 10;
    ctx->maxdistclass = -1;
    ctx->lhlimit = INT_MAX;
    ctx->njob = stringsCount;
    ctx->rnaprediction = 'm';
    ctx->addprofile = 1;
    ctx->outgap = 1;
    ctx->kimuraR = NOTSPECIFIED;
    ctx->pamN = NOTSPECIFIED;
    ctx->fftWinSize = NOTSPECIFIED;
    ctx->fftThreshold = NOTSPECIFIED;
    ctx->RNAppenalty = NOTSPECIFIED;
    ctx->RNAppenalty_ex = NOTSPECIFIED;
    ctx->RNApthr = NOTSPECIFIED;
    ctx->TMorJTT = JTT;
    ctx->consweight_multi = 1.0;

    ctx->ppenalty_OP = NOTSPECIFIED;
    ctx->ppenalty_EX = NOTSPECIFIED;

    // TODO(sen) What do we need the name array for?
    // TODO(sen) Can we use the strings directly as aos?
    // NOTE(sen) Process input
    const char* const* name = aln_arenaAllocArray(tempArena, const char*, stringsCount);
    int*               nlen = aln_arenaAllocArray(tempArena, int, stringsCount);
    for (int32_t strIndex = 0; strIndex < stringsCount; strIndex++) {
        aln_Str str = strings[strIndex];
        ctx->maxInputSeqLen = aln_max(ctx->maxInputSeqLen, str.len);
        ((const char**)name)[strIndex] = "name";
        nlen[strIndex] = str.len;
    }

    // TODO(sen) Remove the requirement for null-terminated array
    char** seq = aln_arenaAllocArray(tempArena, char*, stringsCount + 1);
    for (int32_t strIndex = 0; strIndex < stringsCount; strIndex++) {
        aln_Str str = strings[strIndex];
        // TODO(sen) Remove the requirement for null-terminated str
        char* thisSeq = aln_arenaAllocArray(tempArena, char, ctx->maxInputSeqLen + 1);
        for (int32_t charInd = 0; charInd < str.len; charInd++) {
            thisSeq[charInd] = str.ptr[charInd];
        }
        seq[strIndex] = thisSeq;
    }

    int* targetmap = calloc(ctx->njob, sizeof(int));
    int* targetmapr = calloc(ctx->njob, sizeof(int));
    for (int i = 0; i < ctx->njob; i++) {
        targetmap[i] = targetmapr[i] = i;
    }

    LocalHom** localhomtable = (LocalHom**)calloc(ctx->njob, sizeof(LocalHom*));
    {
        int ilim = ctx->njob;
        for (int i = 0; i < ctx->njob; i++) {
            localhomtable[i] = (LocalHom*)calloc(ilim, sizeof(LocalHom));
            for (int j = 0; j < ilim; j++) {
                localhomtable[i][j].start1 = -1;
                localhomtable[i][j].end1 = -1;
                localhomtable[i][j].start2 = -1;
                localhomtable[i][j].end2 = -1;
                localhomtable[i][j].overlapaa = -1.0;
                localhomtable[i][j].opt = -1.0;
                localhomtable[i][j].importance = -1.0;
                localhomtable[i][j].next = NULL;
                localhomtable[i][j].nokori = 0;
                localhomtable[i][j].extended = -1;
                localhomtable[i][j].last = localhomtable[i] + j;
                localhomtable[i][j].korh = 'h';
            }
            ilim--;
        }
    }

    constants(opts, ctx);

    {
        aln_assert(ctx->njob >= 2);


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

            char** distseq1 = AllocateCharMtx(1, 0);
            char** distseq2 = AllocateCharMtx(1, 0);

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

                    pscore = G__align11(opts, ctx, ctx->n_dis_consweight_multi, mseq1, mseq2, alloclen, ctx->outgap, ctx->outgap);
                    int thereisx = countsOfXs[i] + countsOfXs[j];
                    if (thereisx) {
                        pscore = G__align11_noalign(ctx, ctx->n_dis_consweight_multi, opts.penalty, ctx->penalty_ex, distseq1, distseq2);
                    }
                }
            }
        }
    }

    {
        int ilim = ctx->njob;
        for (int i = 0; i < ctx->njob; i++) {
            for (int j = 0; j < ilim; j++) {
                for (LocalHom* tmpptr = localhomtable[i] + j; tmpptr; tmpptr = tmpptr->next) {
                    if (tmpptr->opt == -1.0)
                        continue;
                    tmpptr->opt = (tmpptr->opt) / 5.8 * 600;
                }
            }
            ilim--;
        }
    }

    // TODO(sen) Sequence verification?

    int***   topol = AllocateIntCub(ctx->njob, 2, 0);
    double** len = AllocateFloatMtx(ctx->njob, 2);
    Treedep* dep = (Treedep*)calloc(ctx->njob, sizeof(Treedep));
    {
        double** iscore = AllocateFloatHalfMtx(ctx->njob);
        fixed_musclesupg_double_realloc_nobk_halfmtx_memsave(opts, ctx, ctx->njob, iscore, topol, len, dep, 1, 1);
    }

    int** localmem = AllocateIntMtx(2, ctx->njob + 1);
    localmem[0][0] = -1;
    topolorderz(localmem[0], topol, dep, ctx->njob - 2, 2);

    double* eff = AllocateDoubleVec(ctx->njob);
    counteff_simple_double_nostatic_memsave(ctx->njob, topol, len, dep, eff);
    for (int i = ctx->njob; i < ctx->njob; i++) {
        eff[i] /= (double)100;
    }

    int    alloclen = ctx->maxInputSeqLen * 2 + 1;
    char** bseq = AllocateCharMtx(ctx->njob, alloclen);

    for (int i = 0; i < ctx->njob; i++) {
        copyWithNoGaps(bseq[i], seq[i]);
    }

    {
        int       m1, m2;
        double    dumdb = 0.0;
        double*** cpmxchild0 = NULL;
        double*** cpmxchild1 = NULL;
        double    orieff1, orieff2;

        char* swaplist = 0;

        int*    fftlog = AllocateIntVec(ctx->njob);
        double* effarr1 = AllocateDoubleVec(ctx->njob);
        double* effarr2 = AllocateDoubleVec(ctx->njob);
        char*   indication1 = AllocateCharVec(150);
        char*   indication2 = AllocateCharVec(150);

        double**  dynamicmtx = AllocateDoubleMtx(ctx->nalphabets, ctx->nalphabets);
        int**     localmem = calloc(sizeof(int*), 2);
        double*** cpmxhist = (double***)calloc(ctx->njob - 1, sizeof(double**));
        int**     memhist = (int**)calloc(ctx->njob - 1, sizeof(int*));
        double*   effarr1_kozo = AllocateDoubleVec(ctx->njob);
        double*   effarr2_kozo = AllocateDoubleVec(ctx->njob);

        LocalHom*** localhomshrink = (LocalHom***)calloc(ctx->njob, sizeof(LocalHom**));
        for (int i = 0; i < ctx->njob; i++) {
            localhomshrink[i] = (LocalHom**)calloc(ctx->njob, sizeof(LocalHom*));
        }

        int* alreadyaligned = AllocateIntVec(ctx->njob);
        for (int i = 0; i < ctx->njob; i++)
            alreadyaligned[i] = 1;
        for (int i = ctx->njob; i < ctx->njob; i++)
            alreadyaligned[i] = 0;

        for (int l = 0; l < ctx->njob; l++)
            fftlog[l] = 1;

        calcimportance_half(ctx, ctx->njob, eff, bseq, localhomtable, alloclen);

        char** mseq1 = AllocateCharMtx(ctx->njob, 0);
        char** mseq2 = AllocateCharMtx(ctx->njob, 0);

        int clus1 = 0;
        int clus2 = 0;

        for (int l = 0; l < ctx->njob - 1; l++) {
            m1 = topol[l][0][0];
            m2 = topol[l][1][0];

            if (dep[l].child0 == -1) {
                cpmxchild0 = NULL;
            } else {
                cpmxchild0 = cpmxhist + dep[l].child0;
            }

            if (dep[l].child1 == -1) {
                cpmxchild1 = NULL;
            } else {
                cpmxchild1 = cpmxhist + dep[l].child1;
            }

            if (dep[l].child0 == -1) {
                localmem[0] = calloc(sizeof(int), 2);
                localmem[0][0] = m1;
                localmem[0][1] = -1;
                clus1 = 1;
            } else {
                localmem[0] = memhist[dep[l].child0];
                clus1 = intlen(localmem[0]);
            }
            if (dep[l].child1 == -1) {
                localmem[1] = calloc(sizeof(int), 2);
                localmem[1][0] = m2;
                localmem[1][1] = -1;
                clus2 = 1;
            } else {
                localmem[1] = memhist[dep[l].child1];
                clus2 = intlen(localmem[1]);
            }

            if (l != ctx->njob - 2) {
                memhist[l] = calloc(sizeof(int), clus1 + clus2 + 1);
                intcpy(memhist[l], localmem[0]);
                intcpy(memhist[l] + clus1, localmem[1]);
                memhist[l][clus1 + clus2] = -1;
            }

            makedynamicmtx(ctx, dynamicmtx, ctx->n_dis_consweight_multi, dep[l].distfromtip);

            int len1 = strlen(bseq[m1]);
            int len2 = strlen(bseq[m2]);

            // TODO(sen) Out of memory error?
            aln_assert(alloclen >= len1 + len2);

            clus1 = fastconjuction_noname(localmem[0], bseq, mseq1, effarr1, eff, indication1, opts.minimumweight, &orieff1);
            clus2 = fastconjuction_noname(localmem[1], bseq, mseq2, effarr2, eff, indication2, opts.minimumweight, &orieff2);

            fastshrinklocalhom_half(localmem[0], localmem[1], localhomtable, localhomshrink);

            imp_match_init_strict(
                opts,
                clus1,
                clus2,
                strlen(mseq1[0]),
                strlen(mseq2[0]),
                mseq1,
                mseq2,
                effarr1,
                effarr2,
                effarr1_kozo,
                effarr2_kozo,
                localhomshrink,
                swaplist,
                localmem[0],
                localmem[1]
            );

            A__align(
                opts,
                ctx,
                dynamicmtx,
                opts.penalty,
                ctx->penalty_ex,
                mseq1,
                mseq2,
                effarr1,
                effarr2,
                clus1,
                clus2,
                alloclen,
                opts.constraint,
                &dumdb,
                NULL,
                NULL,
                NULL,
                NULL,
                ctx->outgap,
                ctx->outgap,
                localmem[0][0],
                1,
                cpmxchild0,
                cpmxchild1,
                cpmxhist + l,
                orieff1,
                orieff2
            );

            nlen[m1] = 0.5 * (nlen[m1] + nlen[m2]);
        }
    }

    aln_Str* alignedSeqs = aln_arenaAllocArray(permArena, aln_Str, stringsCount);
    for (int32_t strIndex = 0; strIndex < stringsCount; strIndex++) {
        aln_Str alignedSeq = {bseq[strIndex], strlen(bseq[strIndex])};
        char*   permCopy = aln_arenaAllocArray(permArena, char, alignedSeq.len);
        for (int32_t charIndex = 0; charIndex < alignedSeq.len; charIndex++) {
            permCopy[charIndex] = alignedSeq.ptr[charIndex];
        }
        alignedSeqs[strIndex] = (aln_Str) {permCopy, alignedSeq.len};
    }

    aln_AlignResult result = {.seqs = alignedSeqs, .seqCount = stringsCount, .bytesWritten = permArena->used};
    return result;
}
