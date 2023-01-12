#include <stdint.h>
#include <stdbool.h>

#define aln_IMPLEMENTATION
#include "mltaln.h"

#define DEBUG 0
#define IODEBUG 0

typedef struct msacompactdistmtxthread_arg {
    Context* ctx;
    int      njob;
    int      thread_no;
    int*     selfscore;
    double** partmtx;
    char**   seq;
    int**    skiptable;
    double*  mindist;
    int*     mindistfrom;
    int*     jobpospt;
} msacompactdistmtxthread_arg_t;

typedef struct TbfastOpts {
    int32_t treein;
    int32_t treeout;
    int32_t distout;
    int32_t noalign;
    int32_t multidist;
    int32_t subalignment;
    int32_t subalignmentoffset;
    int32_t keeplength;
    int32_t ndeleted;
    int32_t mapout;
    int32_t smoothing;
    int32_t callpairlocalalign;
} TbfastOpts;

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

    TbfastOpts  tempOpts_ = {};
    TbfastOpts* tempOpts = &tempOpts_;
    tempOpts->callpairlocalalign = 1;

    char**   bseq = NULL;
    double** iscore = NULL;
    double * eff = NULL, *eff_kozo_mapped = NULL;
    int      i, j;
    int***   topol = NULL;
    double** expdist = NULL;
    Treedep* dep = NULL;
    double** len = NULL;
    FILE*    prep = NULL;
    char*    mergeoralign = NULL;
    int      nsubalignments, maxmem;
    int**    subtable;
    int*     insubtable;
    int*     preservegaps;
    int**    localmem = NULL;

    int          alloclen = 0;
    LocalHom*    tmpptr;
    static char* kozoarivec = NULL;
    int          ntarget = 0;
    int*         targetmap = NULL;
    int*         targetmapr = NULL;
    int          jst, jj;

    int* uselh = NULL;
    int  nseed = 0;
    int* nfilesfornode = NULL;

    int nkozo = 0;

    if (tempOpts->subalignment) {
        readsubalignmentstable(ctx->njob, NULL, NULL, &nsubalignments, &maxmem);
        fprintf(stderr, "nsubalignments = %d\n", nsubalignments);
        fprintf(stderr, "maxmem = %d\n", maxmem);
        subtable = AllocateIntMtx(nsubalignments, maxmem + 1);
        insubtable = AllocateIntVec(ctx->njob);
        for (i = 0; i < ctx->njob; i++)
            insubtable[i] = 0;
        preservegaps = AllocateIntVec(ctx->njob);
        for (i = 0; i < ctx->njob; i++)
            preservegaps[i] = 0;
        readsubalignmentstable(ctx->njob, subtable, preservegaps, NULL, NULL);
    }

    char** mseq1 = AllocateCharMtx(ctx->njob, 0);
    char** mseq2 = AllocateCharMtx(ctx->njob, 0);

    topol = AllocateIntCub(ctx->njob, 2, 0);
    len = AllocateFloatMtx(ctx->njob, 2);
    eff = AllocateDoubleVec(ctx->njob);
    kozoarivec = AllocateCharVec(ctx->njob);

    mergeoralign = AllocateCharVec(ctx->njob);

    dep = (Treedep*)calloc(ctx->njob, sizeof(Treedep));
    localmem = AllocateIntMtx(2, ctx->njob + 1);

    iscore = AllocateFloatHalfMtx(ctx->njob);

    tempOpts->ndeleted = 0;

    ntarget = ctx->njob;
    targetmap = calloc(ctx->njob, sizeof(int));
    targetmapr = calloc(ctx->njob, sizeof(int));
    for (i = 0; i < ctx->njob; i++) {
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
            tempOpts->callpairlocalalign = 1;
            if (expdist)
                FreeDoubleMtx(expdist);
            expdist = NULL;

            {
                for (ilim = ctx->njob, i = 0; i < ntarget; i++) {
                    for (j = 0; j < ilim; j++) {
                        for (tmpptr = localhomtable[i] + j; tmpptr; tmpptr = tmpptr->next) {
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
                    for (i = 0; i < ilim; i++) {
                        jst = i;
                        jj = 0;
                        for (j = jst; j < ctx->njob; j++, jj++) {
                            for (tmpptr = localhomtable[i] + jj; tmpptr; tmpptr = tmpptr->next) {
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

        nkozo = 0;
        for (i = 0; i < ctx->njob; i++) {
            //			fprintf( stderr, "kozoarivec[%d] = %d\n", i, kozoarivec[i] );
            if (kozoarivec[i])
                nkozo++;
        }
        if (nkozo) {
            eff_kozo_mapped = AllocateDoubleVec(ctx->njob);
        }
    }

    constants(opts, ctx, ctx->njob, seq);

    initSignalSM(ctx);
    initFiles(ctx);
    WriteOptions(opts, ctx, ctx->trap_g);

    aln_assert(!(tempOpts->distout && !tempOpts->treeout && tempOpts->noalign));

    // TODO(sen) Sequence verification?

    fixed_musclesupg_double_realloc_nobk_halfmtx_memsave(opts, ctx, ctx->njob, iscore, topol, len, dep, 1, 1);

    localmem[0][0] = -1;
    topolorderz(localmem[0], topol, dep, ctx->njob - 2, 2);

    aln_assert(!(tempOpts->treeout && tempOpts->noalign));

    ctx->weight = 3;
    counteff_simple_double_nostatic_memsave(ctx->njob, topol, len, dep, eff);
    for (i = ctx->njob - ctx->nadd; i < ctx->njob; i++) {
        eff[i] /= (double)100;
    }

    alloclen = ctx->maxInputSeqLen * 2 + 1;
    bseq = AllocateCharMtx(ctx->njob, alloclen);

    for (i = 0; i < ctx->njob; i++) {
        copyWithNoGaps(bseq[i], seq[i]);
    }

    for (i = 0; i < ctx->njob - 1; i++) {
        mergeoralign[i] = 'a';
    }

    {
        int             i, l, m;
        int             len1nocommongap;
        int             len1, len2;
        int             clus1, clus2;
        double          pscore, tscore;
        char *          indication1, *indication2;
        double*         effarr1 = NULL;
        double*         effarr2 = NULL;
        double*         effarr1_kozo = NULL;
        double*         effarr2_kozo = NULL;
        LocalHom***     localhomshrink = NULL;
        int*            seedinlh1 = NULL;
        int*            seedinlh2 = NULL;
        char*           swaplist = NULL;
        int*            fftlog;
        int             m1, m2;
        int*            gaplen;
        int*            gapmap;
        int*            alreadyaligned;
        double          dumdb = 0.0;
        int             ffttry;
        RNApair ***     grouprna1 = NULL, ***grouprna2 = NULL;
        static double** dynamicmtx;
        int             gapmaplen;
        int**           localmem = NULL;
        int             nfiles;
        double***       cpmxhist = NULL;
        int**           memhist = NULL;
        double***       cpmxchild0 = NULL;
        double***       cpmxchild1 = NULL;
        double          orieff1, orieff2;

        if (ctx->rnakozo && ctx->rnaprediction == 'm') {
            grouprna1 = (RNApair***)calloc(ctx->njob, sizeof(RNApair**));
            grouprna2 = (RNApair***)calloc(ctx->njob, sizeof(RNApair**));
        } else {
            grouprna1 = grouprna2 = NULL;
        }

        fftlog = AllocateIntVec(ctx->njob);
        effarr1 = AllocateDoubleVec(ctx->njob);
        effarr2 = AllocateDoubleVec(ctx->njob);
        indication1 = AllocateCharVec(150);
        indication2 = AllocateCharVec(150);
        gaplen = AllocateIntVec(alloclen + 10);
        gapmap = AllocateIntVec(alloclen + 10);
        alreadyaligned = AllocateIntVec(ctx->njob);
        dynamicmtx = AllocateDoubleMtx(ctx->nalphabets, ctx->nalphabets);
        localmem = calloc(sizeof(int*), 2);
        cpmxhist = (double***)calloc(ctx->njob - 1, sizeof(double**));
        for (i = 0; i < ctx->njob - 1; i++)
            cpmxhist[i] = NULL;
        memhist = (int**)calloc(ctx->njob - 1, sizeof(int*));
        for (i = 0; i < ctx->njob - 1; i++)
            memhist[i] = NULL;

        swaplist = NULL;
        if (opts.constraint && ctx->compacttree != 3) {
            localhomshrink = (LocalHom***)calloc(ctx->njob, sizeof(LocalHom**));
            for (i = 0; i < ctx->njob; i++) {
                localhomshrink[i] = (LocalHom**)calloc(ctx->njob, sizeof(LocalHom*));
            }
        } else if (opts.constraint && nseed) {
            localhomshrink = (LocalHom***)calloc(nseed, sizeof(LocalHom**));
            for (i = 0; i < nseed; i++)
                localhomshrink[i] = (LocalHom**)calloc(nseed, sizeof(LocalHom*));

            seedinlh1 = calloc(ctx->njob, sizeof(int));
            seedinlh2 = calloc(ctx->njob, sizeof(int));
        } else if (opts.constraint && nseed == 0) {
            seedinlh1 = NULL;
            seedinlh2 = NULL;
            localhomshrink = NULL;
        }

        effarr1_kozo = AllocateDoubleVec(ctx->njob);
        effarr2_kozo = AllocateDoubleVec(ctx->njob);
        for (i = 0; i < ctx->njob; i++)
            effarr1_kozo[i] = 0.0;
        for (i = 0; i < ctx->njob; i++)
            effarr2_kozo[i] = 0.0;

        for (i = 0; i < ctx->njob - ctx->nadd; i++)
            alreadyaligned[i] = 1;
        for (i = ctx->njob - ctx->nadd; i < ctx->njob; i++)
            alreadyaligned[i] = 0;

        for (l = 0; l < ctx->njob; l++)
            fftlog[l] = 1;

        if (opts.constraint && ctx->compacttree != 3) {
            calcimportance_half(ctx, ctx->njob, eff, bseq, localhomtable, alloclen);
        } else if (opts.constraint && nseed) {
            dontcalcimportance_half(ctx, nseed, bseq, localhomtable);
        }

        tscore = 0.0;
        for (l = 0; l < ctx->njob - 1; l++) {
            m1 = topol[l][0][0];
            m2 = topol[l][1][0];

            if (eff_kozo_mapped) {
                cpmxchild0 = NULL;
                cpmxchild1 = NULL;
            } else {
                if (dep[l].child0 == -1)
                    cpmxchild0 = NULL;
                else
                    cpmxchild0 = cpmxhist + dep[l].child0;
                if (dep[l].child1 == -1)
                    cpmxchild1 = NULL;
                else
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

            if (mergeoralign[l] == 'n') {
                free(localmem[0]);
                free(localmem[1]);
                continue;
            }

            makedynamicmtx(opts, ctx, dynamicmtx, ctx->n_dis_consweight_multi, dep[l].distfromtip);

            len1 = strlen(bseq[m1]);
            len2 = strlen(bseq[m2]);
            if (alloclen < len1 + len2) {
                fprintf(stderr, "\nReallocating..");
                alloclen = (len1 + len2) + 1000;
                ReallocateCharMtx(bseq, ctx->njob, alloclen + 10);
                gaplen = realloc(gaplen, (alloclen + 10) * sizeof(int));
                if (gaplen == NULL) {
                    fprintf(stderr, "Cannot realloc gaplen\n");
                    exit(1);
                }
                gapmap = realloc(gapmap, (alloclen + 10) * sizeof(int));
                if (gapmap == NULL) {
                    fprintf(stderr, "Cannot realloc gapmap\n");
                    exit(1);
                }
                fprintf(stderr, "done. alloclen = %d\n", alloclen);
            }

            if (eff_kozo_mapped) {
                clus1 = fastconjuction_noname_kozo(localmem[0], bseq, mseq1, effarr1, eff, effarr1_kozo, eff_kozo_mapped, indication1);
                clus2 = fastconjuction_noname_kozo(localmem[1], bseq, mseq2, effarr2, eff, effarr2_kozo, eff_kozo_mapped, indication2);
            } else {
                clus1 = fastconjuction_noname(localmem[0], bseq, mseq1, effarr1, eff, indication1, opts.minimumweight, &orieff1);  // orieff tsukau!
                clus2 = fastconjuction_noname(localmem[1], bseq, mseq2, effarr2, eff, indication2, opts.minimumweight, &orieff2);  // orieff tsukau!
            }

            if (mergeoralign[l] == '1' || mergeoralign[l] == '2') {
                ctx->newgapstr = "=";
            } else
                ctx->newgapstr = "-";

            len1nocommongap = len1;
            if (mergeoralign[l] == '1')  // nai
            {
                findcommongaps(clus2, mseq2, gapmap);
                commongappick(clus2, mseq2);
            } else if (mergeoralign[l] == '2') {
                findcommongaps(clus1, mseq1, gapmap);
                commongappick(clus1, mseq1);
                len1nocommongap = strlen(mseq1[0]);
            }

            if (ctx->compacttree == 3)
                nfiles = nfilesfornode[l];
            else
                nfiles = 0;

            if (l < 1000 || l % 100 == 0)
                fprintf(stderr, "\rSTEP % 5d /%d ", l + 1, ctx->njob - 1);

            if (opts.constraint && ctx->compacttree != 3) {
                fastshrinklocalhom_half(localmem[0], localmem[1], localhomtable, localhomshrink);
            } else if (opts.constraint && nseed) {
                fastshrinklocalhom_half_seed(localmem[0], localmem[1], nseed, seedinlh1, seedinlh2, localhomtable, localhomshrink);
                for (i = 0; i < ctx->njob; i++)
                    reporterr("seedinlh1[%d]=%d\n", i, seedinlh1[i]);
                for (i = 0; i < ctx->njob; i++)
                    reporterr("seedinlh2[%d]=%d\n", i, seedinlh2[i]);
            }

            if (ctx->rnakozo && ctx->rnaprediction == 'm') {
                makegrouprna(grouprna1, NULL, localmem[0]);
                makegrouprna(grouprna2, NULL, localmem[1]);
            }

            if (!ctx->nevermemsave && (opts.constraint != 2 && opts.alg != 'M' && (len1 > 30000 || len2 > 30000))) {
                aln_assert(!"should not execute");
            }

            if (fftlog[m1] && fftlog[m2])
                ffttry = (nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 1000 && clus2 < 1000);
            else
                ffttry = 0;
            if (opts.constraint == 2) {
                if (opts.alg == 'M') {
                    fprintf(stderr, "\n\nMemory saving mode is not supported.\n\n");
                    exit(1);
                }
                if (opts.alg == 'A') {
                    imp_match_init_strict(opts, ctx, clus1, clus2, strlen(mseq1[0]), strlen(mseq2[0]), mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, swaplist, localmem[0], localmem[1], uselh, seedinlh1, seedinlh2, (ctx->compacttree == 3) ? l : -1, nfiles);
                    if (ctx->rnakozo)
                        imp_rna(ctx, clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2);
                    pscore = A__align(opts, ctx, dynamicmtx, ctx->penalty, ctx->penalty_ex, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, opts.constraint, &dumdb, NULL, NULL, NULL, NULL, ctx->outgap, ctx->outgap, localmem[0][0], 1, cpmxchild0, cpmxchild1, cpmxhist + l, orieff1, orieff2);
                }
                if (opts.alg == 'd') {
                    imp_match_init_strictD(opts, ctx, clus1, clus2, strlen(mseq1[0]), strlen(mseq2[0]), mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, swaplist, localmem[0], localmem[1], uselh, seedinlh1, seedinlh2, (ctx->compacttree == 3) ? l : -1, nfiles);
                    if (ctx->rnakozo)
                        imp_rnaD(ctx, clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2);
                    pscore = D__align(opts, ctx, dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, opts.constraint, &dumdb, ctx->outgap, ctx->outgap);
                } else if (opts.alg == 'Q') {
                    aln_assert(!"not supported");
                }
            } else if (ctx->force_fft || (opts.use_fft && ffttry)) {
                fprintf(stderr, " f\b\b");
                if (opts.alg == 'M') {
                    fprintf(stderr, "m");
                    pscore = Falign_udpari_long(opts, ctx, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, alloclen, fftlog + m1);
                } else
                    pscore = Falign(opts, ctx, NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, alloclen, fftlog + m1);
            } else {
                fprintf(stderr, " d\b\b");
                fftlog[m1] = 0;
                switch (opts.alg) {
                    case ('a'):
                        pscore = Aalign(ctx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen);
                        break;
                    case ('M'):
                        fprintf(stderr, "m");
                        pscore = MSalignmm(opts, ctx, dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, NULL, NULL, NULL, NULL, ctx->outgap, ctx->outgap, cpmxchild0, cpmxchild1, cpmxhist + l, orieff1, orieff2);
                        break;
                    case ('A'):
                        pscore = A__align(opts, ctx, dynamicmtx, ctx->penalty, ctx->penalty_ex, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, 0, &dumdb, NULL, NULL, NULL, NULL, ctx->outgap, ctx->outgap, localmem[0][0], 1, cpmxchild0, cpmxchild1, cpmxhist + l, orieff1, orieff2);
                        break;
                    case ('d'):
                        pscore = D__align(opts, ctx, dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, 0, &dumdb, ctx->outgap, ctx->outgap);
                        break;
                    default:
                        ErrorExit("ERROR IN SOURCE FILE");
                }
            }

            nlen[m1] = 0.5 * (nlen[m1] + nlen[m2]);

            tscore += pscore;

            if (ctx->disp)
                display(ctx, bseq, ctx->njob);

            if (mergeoralign[l] == '1') {
                reporterr("Check source!!\n");
                exit(1);
            }
            if (mergeoralign[l] == '2') {
                gapmaplen = strlen(mseq1[0]) - len1nocommongap + len1;
                adjustgapmap(gapmaplen, gapmap, mseq1[0]);
                if (tempOpts->smoothing) {
                    restorecommongapssmoothly(ctx->njob, ctx->njob - (clus1 + clus2), bseq, localmem[0], localmem[1], gapmap, alloclen, '-');
                    findnewgaps(ctx, 0, mseq1, gaplen);
                    insertnewgaps_bothorders(opts, ctx, ctx->njob, alreadyaligned, bseq, localmem[0], localmem[1], gaplen, gapmap, gapmaplen, alloclen, opts.alg, '-');
                } else {
                    restorecommongaps(ctx->njob, ctx->njob - (clus1 + clus2), bseq, localmem[0], localmem[1], gapmap, alloclen, '-');
                    findnewgaps(ctx, 0, mseq1, gaplen);
                    insertnewgaps(opts, ctx, ctx->njob, alreadyaligned, bseq, localmem[0], localmem[1], gaplen, gapmap, alloclen, opts.alg, '-');
                }
                eq2dashmatometehayaku(mseq1, clus1);
                eq2dashmatometehayaku(mseq2, clus2);
                for (i = 0; (m = localmem[1][i]) > -1; i++)
                    alreadyaligned[m] = 1;
            }

            free(localmem[0]);
            free(localmem[1]);
        }

        if (cpmxhist) {
            for (i = 0; i < ctx->njob - 1; i++) {
                if (cpmxhist[i]) {
                    FreeDoubleMtx(cpmxhist[i]);
                    cpmxhist[i] = NULL;
                }
            }
            free(cpmxhist);
            cpmxhist = NULL;
        }

        free(memhist);
        memhist = NULL;

        bool scoreout = false;
        if (scoreout) {
            fprintf(stderr, "totalscore = %10.2f\n\n", tscore);
        }

        if (ctx->rnakozo && ctx->rnaprediction == 'm') {
            if (grouprna1)
                free(grouprna1);  // nakami ha?
            if (grouprna2)
                free(grouprna2);  // nakami ha?
            grouprna1 = grouprna2 = NULL;
        }

        if (opts.constraint) {
            if (localhomshrink)  // nen no tame
            {
                if (ctx->compacttree == 3)
                    m = nseed;
                else
                    m = ctx->njob;
                for (i = 0; i < m; i++) {
                    free(localhomshrink[i]);
                    localhomshrink[i] = NULL;
                }
                free(localhomshrink);
                localhomshrink = NULL;
            }
            if (seedinlh1)
                free(seedinlh1);
            if (seedinlh2)
                free(seedinlh2);
        }

        Falign(opts, ctx, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL);
        Falign_udpari_long(opts, ctx, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL);
        D__align(opts, ctx, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, 0, 0);
        A__align(opts, ctx, NULL, 0, 0, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, 0, 0, -1, -1, NULL, NULL, NULL, 0.0, 0.0);
        imp_match_init_strictD(opts, ctx, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0);
        imp_match_init_strict(opts, ctx, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0);
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
