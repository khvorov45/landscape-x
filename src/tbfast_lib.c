#include <stdint.h>
#include <stdbool.h>

#define aln_IMPLEMENTATION
#include "mltaln.h"

#define DEBUG 0
#define IODEBUG 0

static void
WriteOptions(aln_Opts opts, Context* ctx, FILE* fp) {
    if (ctx->dorp == 'd')
        fprintf(fp, "DNA\n");
    else {
        if (opts.scoremtx == 0)
            fprintf(fp, "JTT %dPAM\n", ctx->pamN);
        else if (opts.scoremtx == 1)
            fprintf(fp, "BLOSUM %d\n", opts.nblosum);
        else if (opts.scoremtx == 2)
            fprintf(fp, "M-Y\n");
    }
    fprintf(stderr, "Gap Penalty = %+5.2f, %+5.2f, %+5.2f\n", (double)opts.ppenalty / 1000, (double)opts.ppenalty_ex / 1000, (double)opts.poffset / 1000);
    if (opts.use_fft)
        fprintf(fp, "FFT on\n");

    fprintf(fp, "tree-base method\n");
    if (ctx->tbrweight == 0)
        fprintf(fp, "unweighted\n");
    else if (ctx->tbrweight == 3)
        fprintf(fp, "clustalw-like weighting\n");
    if (ctx->tbitr || ctx->tbweight) {
        fprintf(fp, "iterate at each step\n");
        if (ctx->tbitr && ctx->tbrweight == 0)
            fprintf(fp, "  unweighted\n");
        if (ctx->tbitr && ctx->tbrweight == 3)
            fprintf(fp, "  reversely weighted\n");
        if (ctx->tbweight)
            fprintf(fp, "  weighted\n");
        fprintf(fp, "\n");
    }

    fprintf(fp, "Gap Penalty = %+5.2f, %+5.2f, %+5.2f\n", (double)opts.ppenalty / 1000, (double)opts.ppenalty_ex / 1000, (double)opts.poffset / 1000);

    if (opts.alg == 'a')
        fprintf(fp, "Algorithm A\n");
    else if (opts.alg == 'A')
        fprintf(fp, "Algorithm A+\n");
    else if (opts.alg == 'C')
        fprintf(fp, "Apgorithm A+/C\n");
    else
        fprintf(fp, "Unknown algorithm\n");

    if (opts.treemethod == 'X')
        fprintf(fp, "Tree = UPGMA (mix).\n");
    else if (opts.treemethod == 'E')
        fprintf(fp, "Tree = UPGMA (average).\n");
    else if (opts.treemethod == 'q')
        fprintf(fp, "Tree = Minimum linkage.\n");
    else
        fprintf(fp, "Unknown tree.\n");

    if (opts.use_fft) {
        fprintf(fp, "FFT on\n");
        if (ctx->dorp == 'd')
            fprintf(fp, "Basis : 4 nucleotides\n");
        else {
            if (ctx->fftscore)
                fprintf(fp, "Basis : Polarity and Volume\n");
            else
                fprintf(fp, "Basis : 20 amino acids\n");
        }
        fprintf(fp, "Threshold   of anchors = %d%%\n", ctx->fftThreshold);
        fprintf(fp, "window size of anchors = %dsites\n", ctx->fftWinSize);
    } else
        fprintf(fp, "FFT off\n");
    fflush(fp);
}

aln_AlignResult
tbfast_main(aln_Str* strings, intptr_t stringsCount, void* out, intptr_t outBytes, aln_Opts opts) {
    aln_Arena permArena_ = {.base = out, .size = outBytes / 4};
    aln_Arena tempArena_ = {.base = (uint8_t*)out + permArena_.size, .size = outBytes - permArena_.size};

    aln_Arena* permArena = &permArena_;
    aln_Arena* tempArena = &tempArena_;

    aln_Opts pairLocalAlignOpts = opts;
    pairLocalAlignOpts.ppenalty = -2000;
    pairLocalAlignOpts.poffset = 100;
    pairLocalAlignOpts.use_fft = 0;
    pairLocalAlignOpts.constraint = 0;

    Context* ctx = aln_arenaAllocStruct(tempArena, Context);

    // TODO(sen) This is 'dna or protein'. Figure out what to do with this
    // and why this even matters
    ctx->dorp = 'p';

    ctx->RNAscoremtx = 'n';
    ctx->parallelizationstrategy = BAATARI1;
    ctx->newgapstr = "-";
    ctx->nalphabets = 26;
    ctx->nscoredalphabets = 20;
    ctx->ndistclass = 10;
    ctx->maxdistclass = -1;
    ctx->lhlimit = INT_MAX;
    ctx->njob = stringsCount;
    ctx->rnaprediction = 'm';
    ctx->addprofile = 1;
    ctx->fftscore = 1;
    ctx->weight = 3;
    ctx->tbutree = 1;
    ctx->outgap = 1;
    ctx->tbrweight = 3;
    ctx->kimuraR = NOTSPECIFIED;
    ctx->pamN = NOTSPECIFIED;
    ctx->geta2 = GETA2;
    ctx->fftWinSize = NOTSPECIFIED;
    ctx->fftThreshold = NOTSPECIFIED;
    ctx->RNAppenalty = NOTSPECIFIED;
    ctx->RNAppenalty_ex = NOTSPECIFIED;
    ctx->RNApthr = NOTSPECIFIED;
    ctx->TMorJTT = JTT;
    ctx->consweight_multi = 1.0;

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

    char**   bseq = NULL;
    double** iscore = NULL;
    double** expdist = NULL;
    Treedep* dep = NULL;
    FILE*    prep = NULL;
    int**    localmem = NULL;

    int       alloclen = 0;
    int       ntarget = 0;
    int*      targetmap = NULL;
    int*      targetmapr = NULL;
    int*      uselh = NULL;

    char**   mseq1 = AllocateCharMtx(ctx->njob, 0);
    char**   mseq2 = AllocateCharMtx(ctx->njob, 0);
    int***   topol = AllocateIntCub(ctx->njob, 2, 0);
    double** len = AllocateFloatMtx(ctx->njob, 2);
    double*  eff = AllocateDoubleVec(ctx->njob);

    dep = (Treedep*)calloc(ctx->njob, sizeof(Treedep));
    localmem = AllocateIntMtx(2, ctx->njob + 1);

    iscore = AllocateFloatHalfMtx(ctx->njob);

    ntarget = ctx->njob;
    targetmap = calloc(ctx->njob, sizeof(int));
    targetmapr = calloc(ctx->njob, sizeof(int));
    for (int i = 0; i < ctx->njob; i++) {
        targetmap[i] = targetmapr[i] = i;
    }

    LocalHom** localhomtable = (LocalHom**)calloc(ntarget, sizeof(LocalHom*));
    {
        int ilim = ctx->njob;
        for (int i = 0; i < ntarget; i++) {
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

        {
            pairlocalalign(pairLocalAlignOpts, ctx, ctx->njob, name, seq, iscore, localhomtable, expdist);
            if (expdist)
                FreeDoubleMtx(expdist);
            expdist = NULL;

            {
                ilim = ctx->njob;
                for (int i = 0; i < ntarget; i++) {
                    for (int j = 0; j < ilim; j++) {
                        for (LocalHom* tmpptr = localhomtable[i] + j; tmpptr; tmpptr = tmpptr->next) {
                            if (tmpptr->opt == -1.0)
                                continue;
                            tmpptr->opt = (tmpptr->opt) / 5.8 * 600;
                        }
                    }
                    ilim--;
                }

                {
                    prep = fopen("hat3", "w");
                    if (!prep)
                        ErrorExit("Cannot open hat3 to write.");

                    fprintf(stderr, "Writing hat3 for iterative refinement\n");
                    ilim = ctx->njob - 1;
                    for (int i = 0; i < ilim; i++) {
                        int jst = i;
                        int jj = 0;
                        for (int j = jst; j < ctx->njob; j++, jj++) {
                            for (LocalHom* tmpptr = localhomtable[i] + jj; tmpptr; tmpptr = tmpptr->next) {
                                if (tmpptr->opt == -1.0)
                                    continue;
                                if (targetmap[j] == -1 || targetmap[i] < targetmap[j])
                                    fprintf(prep, "%d %d %d %7.5f %d %d %d %d %c\n", targetmapr[i], j, tmpptr->overlapaa, tmpptr->opt / 600 * 5.8, tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, tmpptr->korh);
                            }
                        }
                    }
                    fclose(prep);

                    prep = fopen("hat2", "w");
                    WriteFloatHat2_pointer_halfmtx(ctx, prep, ctx->njob, name, iscore);
                    fclose(prep);
                }
            }
        }
    }

    constants(opts, ctx, ctx->njob, seq);

    initSignalSM(ctx);
    initFiles(ctx);
    WriteOptions(opts, ctx, ctx->trap_g);

    // TODO(sen) Sequence verification?

    fixed_musclesupg_double_realloc_nobk_halfmtx_memsave(opts, ctx, ctx->njob, iscore, topol, len, dep, 1, 1);

    localmem[0][0] = -1;
    topolorderz(localmem[0], topol, dep, ctx->njob - 2, 2);

    ctx->weight = 3;
    counteff_simple_double_nostatic_memsave(ctx->njob, topol, len, dep, eff);
    for (int i = ctx->njob - ctx->nadd; i < ctx->njob; i++) {
        eff[i] /= (double)100;
    }

    alloclen = ctx->maxInputSeqLen * 2 + 1;
    bseq = AllocateCharMtx(ctx->njob, alloclen);

    for (int i = 0; i < ctx->njob; i++) {
        copyWithNoGaps(bseq[i], seq[i]);
    }

    {
        int       len1, len2;
        int       clus1, clus2;
        int*      seedinlh1 = NULL;
        int*      seedinlh2 = NULL;
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
        for (int i = 0; i < ctx->njob - ctx->nadd; i++)
            alreadyaligned[i] = 1;
        for (int i = ctx->njob - ctx->nadd; i < ctx->njob; i++)
            alreadyaligned[i] = 0;

        for (int l = 0; l < ctx->njob; l++)
            fftlog[l] = 1;

        calcimportance_half(ctx, ctx->njob, eff, bseq, localhomtable, alloclen);

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

            makedynamicmtx(opts, ctx, dynamicmtx, ctx->n_dis_consweight_multi, dep[l].distfromtip);

            len1 = strlen(bseq[m1]);
            len2 = strlen(bseq[m2]);

            // TODO(sen) Out of memory error?
            aln_assert(alloclen >= len1 + len2);

            clus1 = fastconjuction_noname(localmem[0], bseq, mseq1, effarr1, eff, indication1, opts.minimumweight, &orieff1);  // orieff tsukau!
            clus2 = fastconjuction_noname(localmem[1], bseq, mseq2, effarr2, eff, indication2, opts.minimumweight, &orieff2);  // orieff tsukau!

            fastshrinklocalhom_half(localmem[0], localmem[1], localhomtable, localhomshrink);

            imp_match_init_strict(opts, ctx, clus1, clus2, strlen(mseq1[0]), strlen(mseq2[0]), mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, swaplist, localmem[0], localmem[1], uselh, seedinlh1, seedinlh2, (ctx->compacttree == 3) ? l : -1, 0);
            A__align(opts, ctx, dynamicmtx, ctx->penalty, ctx->penalty_ex, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, opts.constraint, &dumdb, NULL, NULL, NULL, NULL, ctx->outgap, ctx->outgap, localmem[0][0], 1, cpmxchild0, cpmxchild1, cpmxhist + l, orieff1, orieff2);

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
