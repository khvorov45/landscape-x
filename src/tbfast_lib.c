#include <stdint.h>
#include <stdbool.h>

#include "mltaln.h"

#define DEBUG 0
#define IODEBUG 0
#define SCOREOUT 0
#define SHISHAGONYU 0  // for debug

#define REPORTCOSTS 0

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
    int32_t outputhat23;
    int32_t nthreadtb;
} TbfastOpts;

static void
arguments(Context* ctx, TbfastOpts* opts, int argc, char* argv[], int* pac, char** pav, int* tac, char** tav)  // 2 kai yobaremasu.
{
    int c;
    int i;

    ctx->nthread = 1;
    opts->nthreadtb = 1;
    ctx->nthreadpair = 1;
    ctx->outnumber = 0;
    ctx->scoreout = 0;
    ctx->spscoreout = 0;
    ctx->rnaprediction = 'm';
    ctx->rnakozo = 0;
    ctx->nevermemsave = 0;
    ctx->inputfile = NULL;
    ctx->addfile = NULL;
    ctx->addprofile = 1;
    ctx->fftkeika = 0;
    ctx->constraint = 0;
    ctx->nblosum = 62;
    ctx->fmodel = 0;
    ctx->use_fft = 0;  // chuui
    ctx->force_fft = 0;
    ctx->fftscore = 1;
    ctx->fftRepeatStop = 0;
    ctx->fftNoAnchStop = 0;
    ctx->weight = 3;
    ctx->tbutree = 1;
    ctx->disp = 0;
    ctx->outgap = 1;
    ctx->alg = 'A';
    ctx->mix = 0;
    ctx->tbitr = 0;
    ctx->tbweight = 0;
    ctx->tbrweight = 3;
    ctx->checkC = 0;
    ctx->treemethod = 'X';
    ctx->sueff_global = 0.1;
    ctx->scoremtx = 1;
    ctx->kobetsubunkatsu = 0;
    ctx->ppenalty_dist = NOTSPECIFIED;
    ctx->ppenalty = NOTSPECIFIED;
    ctx->penalty_shift_factor = 1000.0;
    ctx->ppenalty_ex = NOTSPECIFIED;
    ctx->poffset = NOTSPECIFIED;
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
    ctx->consweight_rna = 0.0;
    ctx->legacygapcost = 0;
    ctx->specificityconsideration = 0.0;
    specifictarget = 0;
    ctx->nwildcard = 0;
    nadd = 0;

    if (pac) {
        pav[0] = "tbfast-pair";
        *pac = 1;
        tav[0] = "tbfast";
        *tac = 1;

        for (i = 0; i < argc; i++) {
            if (argv[i][0] == '_') {
                opts->callpairlocalalign = 1;

                for (i++; i < argc; i++) {
                    if (argv[i][0] == '_') {
                        argc -= 2;
                        argv += 2;

                        goto pavend;
                    }
                    pav[*pac] = argv[i];
                    *pac += 1;
                }
            }
        }

    pavend:

        for (i -= 1; i < argc; i++) {
            tav[*tac] = argv[i];
            *tac += 1;
        }

        argc -= *pac - 1;
        argv += *pac - 1;
    }

    while (--argc > 0 && (*++argv)[0] == '-') {
        while ((c = *++argv[0])) {
            switch (c) {
                case 'i':
                    ctx->inputfile = *++argv;
                    --argc;
                    goto nextoption;
                case 'I':
                    nadd = myatoi(*++argv);
                    --argc;
                    goto nextoption;
                case 'e':
                    ctx->RNApthr = (int)(atof(*++argv) * 1000 - 0.5);
                    --argc;
                    goto nextoption;
                case 'o':
                    ctx->RNAppenalty = (int)(atof(*++argv) * 1000 - 0.5);
                    --argc;
                    goto nextoption;
                case 'V':
                    ctx->ppenalty_dist = (int)(atof(*++argv) * 1000 - 0.5);
                    --argc;
                    goto nextoption;
                case 'f':
                    ctx->ppenalty = (int)(atof(*++argv) * 1000 - 0.5);
                    --argc;
                    goto nextoption;
                case 'Q':
                    ctx->penalty_shift_factor = atof(*++argv);
                    --argc;
                    goto nextoption;
                case 'g':
                    ctx->ppenalty_ex = (int)(atof(*++argv) * 1000 - 0.5);
                    --argc;
                    goto nextoption;
                case 'h':
                    ctx->poffset = (int)(atof(*++argv) * 1000 - 0.5);
                    //					fprintf( stderr, "poffset = %d\n", poffset );
                    --argc;
                    goto nextoption;
                case 'k':
                    ctx->kimuraR = myatoi(*++argv);
                    //					fprintf( stderr, "kappa = %d\n", kimuraR );
                    --argc;
                    goto nextoption;
                case 'b':
                    ctx->nblosum = myatoi(*++argv);
                    ctx->scoremtx = 1;
                    //					fprintf( stderr, "blosum %d / kimura 200\n", nblosum );
                    --argc;
                    goto nextoption;
                case 'j':
                    ctx->pamN = myatoi(*++argv);
                    ctx->scoremtx = 0;
                    ctx->TMorJTT = JTT;
                    //					fprintf( stderr, "jtt/kimura %d\n", pamN );
                    --argc;
                    goto nextoption;
                case 'm':
                    ctx->pamN = myatoi(*++argv);
                    ctx->scoremtx = 0;
                    ctx->TMorJTT = TM;
                    //					fprintf( stderr, "tm %d\n", pamN );
                    --argc;
                    goto nextoption;
                case 'l':
                    ctx->fastathreshold = atof(*++argv);
                    ctx->constraint = 2;
                    --argc;
                    goto nextoption;
                case 'r':
                    ctx->consweight_rna = atof(*++argv);
                    ctx->rnakozo = 1;
                    --argc;
                    goto nextoption;
                case 'c':
                    ctx->consweight_multi = atof(*++argv);
                    --argc;
                    goto nextoption;
                case 'C':
                    ctx->nthreadpair = ctx->nthread = myatoi(*++argv);
                    --argc;
                    ctx->nthread = 0;
                    goto nextoption;
                case 's':
                    ctx->specificityconsideration = (double)myatof(*++argv);
                    //					fprintf( stderr, "specificityconsideration = %f\n", specificityconsideration );
                    --argc;
                    goto nextoption;
                case 'R':
                    ctx->rnaprediction = 'r';
#if 1
                case 'a':
                    ctx->fmodel = 1;
                    break;
#endif
                case 'K':
                    ctx->addprofile = 0;
                    break;
                case 'y':
                    opts->distout = 1;
                    break;
                case 't':
                    opts->treeout = 1;
                    break;
                case '^':
                    opts->treeout = 2;
                    break;
                case 'T':
                    opts->noalign = 1;
                    break;
                case 'D':
                    ctx->dorp = 'd';
                    break;
                case 'P':
                    ctx->dorp = 'p';
                    break;
                case 'L':
                    ctx->legacygapcost = 1;
                    break;
#if 1
                case 'O':
                    ctx->outgap = 0;
                    break;
#else
                case 'O':
                    fftNoAnchStop = 1;
                    break;
#endif
#if 0
				case 'S' :
					scoreout = 1; // for checking parallel calculation
					break;
#else
                case 'S':
                    ctx->spscoreout = 1;  // 2014/Dec/30, sp score
                    break;
#endif
                case 'H':
                    opts->subalignment = 1;
                    opts->subalignmentoffset = myatoi(*++argv);
                    --argc;
                    goto nextoption;
#if 0
				case 'e':
					fftscore = 0;
					break;
				case 'r':
					fmodel = -1;
					break;
				case 'R':
					fftRepeatStop = 1;
					break;
				case 's':
					treemethod = 's';
					break;
#endif
                case 'X':
                    ctx->treemethod = 'X';
                    ctx->sueff_global = atof(*++argv);
                    //					fprintf( stderr, "sueff_global = %f\n", sueff_global );
                    --argc;
                    goto nextoption;
                case 'E':
                    ctx->treemethod = 'E';
                    break;
                case 'q':
                    ctx->treemethod = 'q';
                    break;
                case 'n':
                    ctx->outnumber = 1;
                    break;
#if 0
				case 'a':
					alg = 'a';
					break;
				case 'H':
					alg = 'H';
					break;
				case 'Q':
					alg = 'Q';
					break;
#endif
                case '@':
                    ctx->alg = 'd';
                    break;
                case 'A':
                    ctx->alg = 'A';
                    break;
                case 'M':
                    ctx->alg = 'M';
                    break;
                case 'N':
                    ctx->nevermemsave = 1;
                    break;
                case 'B':  // hitsuyou! memopt -M -B no tame
                    break;
                case 'F':
                    ctx->use_fft = 1;
                    break;
                case 'G':
                    ctx->force_fft = 1;
                    ctx->use_fft = 1;
                    break;
                case 'U':
                    opts->treein = 1;
                    break;
                case 'u':
                    ctx->tbrweight = 0;
                    ctx->weight = 0;
                    break;
                case 'v':
                    ctx->tbrweight = 3;
                    break;
                case 'd':
                    opts->multidist = 1;
                    break;
#if 0
				case 'd':
					disp = 1;
					break;
#endif
                    /* Modified 01/08/27, default: user tree */
                case 'J':
                    ctx->tbutree = 0;
                    break;
/* modification end. */
#if 0
				case 'z':
					fftThreshold = myatoi( *++argv );
					--argc; 
					goto nextoption;
#endif
                case 'w':
                    ctx->fftWinSize = myatoi(*++argv);
                    --argc;
                    goto nextoption;
                case 'W':
                    ctx->minimumweight = atof(*++argv);
                    //					fprintf( stderr, "minimumweight = %f\n", minimumweight );
                    --argc;
                    goto nextoption;
                case 'Y':
                    opts->keeplength = 1;
                    break;
                case 'z':
                    opts->mapout = 2;
                    break;
                case 'Z':
                    opts->mapout = 1;
                    break;
                case 'p':
                    opts->smoothing = 1;
                    break;
                case '=':
                    specifictarget = 1;
                    break;
                case ':':
                    ctx->nwildcard = 1;
                    break;
                case '+':
                    opts->outputhat23 = myatoi(*++argv);
                    reporterr("outputhat23=%d\n", opts->outputhat23);
                    --argc;
                    goto nextoption;
                default:
                    // TODO(sen) This needs to be fixed
                    fprintf(stderr, "illegal option %c\n", c);
                    argc = 0;
                    break;
            }
        }
    nextoption:
        ;
    }

    if (argc == 1) {
        argc--;
    }
    if (argc != 0) {
        fprintf(stderr, "argc=%d, tbfast options: Check source file !\n", argc);
        exit(1);
    }
    if (ctx->tbitr == 1 && ctx->outgap == 0) {
        fprintf(stderr, "conflicting options : o, m or u\n");
        exit(1);
    }
    if (ctx->alg == 'C' && ctx->outgap == 0) {
        fprintf(stderr, "conflicting options : C, o\n");
        exit(1);
    }
}

static double
preferenceval(int ori, int pos, int max) {
    pos -= ori;
    if (pos < 0)
        pos += max;
    return (0.00000000000001 * pos);
}

static void*
msacompactdisthalfmtxthread(msacompactdistmtxthread_arg_t* targ) {
    Context* ctx = targ->ctx;
    int      njob = targ->njob;
    int      thread_no = targ->thread_no;
    int*     selfscore = targ->selfscore;
    double** partmtx = targ->partmtx;
    char**   seq = targ->seq;
    int**    skiptable = targ->skiptable;
    double*  mindist = targ->mindist;
    int*     mindistfrom = targ->mindistfrom;
    int*     jobpospt = targ->jobpospt;
    double   tmpdist, preference, tmpdistx, tmpdisty;
    int      i, j;

    while (1) {
        {
            i = *jobpospt;
            if (i == njob - 1) {
                return (NULL);
            }
            *jobpospt = i + 1;
        }

        if (i % 100 == 0) {
            if (ctx->nthreadpair)
                fprintf(stderr, "\r% 5d / %d (thread %4d)", i, njob, thread_no);
            else
                fprintf(stderr, "\r% 5d / %d", i, njob);
        }

        for (j = i + 1; j < njob; j++) {
            tmpdist = distcompact_msa(ctx, seq[i], seq[j], skiptable[i], skiptable[j], selfscore[i], selfscore[j]);  // osoikedo,

            preference = preferenceval(i, j, njob);
            tmpdistx = tmpdist + preference;
            if (tmpdistx < mindist[i]) {
                mindist[i] = tmpdistx;
                mindistfrom[i] = j;
            }

            preference = preferenceval(j, i, njob);
            tmpdisty = tmpdist + preference;
            if (tmpdisty < mindist[j]) {
                mindist[j] = tmpdisty;
                mindistfrom[j] = i;
            }
            if (partmtx[i])
                partmtx[i][j] = tmpdist;
            if (partmtx[j])
                partmtx[j][i] = tmpdist;
        }
    }
}

void
treebase(Context* ctx, TbfastOpts* opts, int* nlen, char** aseq, int nadd, char* mergeoralign, char** mseq1, char** mseq2, int*** topol, Treedep* dep, double* effarr, int* alloclen, LocalHom** localhomtable, RNApair*** singlerna, double* effarr_kozo, int* targetmap, int* targetmapr, int ntarget, int* uselh, int nseed, int* nfilesfornode) {
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
#if REPORTCOSTS
    time_t starttime, startclock;
    starttime = time(NULL);
    startclock = clock();
#endif

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
    gaplen = AllocateIntVec(*alloclen + 10);
    gapmap = AllocateIntVec(*alloclen + 10);
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
    if (ctx->constraint && compacttree != 3) {
        if (specifictarget)
            swaplist = calloc(ctx->njob, sizeof(char));
        localhomshrink = (LocalHom***)calloc(ctx->njob, sizeof(LocalHom**));
        for (i = 0; i < ctx->njob; i++) {
            localhomshrink[i] = (LocalHom**)calloc(ctx->njob, sizeof(LocalHom*));
        }
    } else if (ctx->constraint && nseed) {
        localhomshrink = (LocalHom***)calloc(nseed, sizeof(LocalHom**));
        for (i = 0; i < nseed; i++)
            localhomshrink[i] = (LocalHom**)calloc(nseed, sizeof(LocalHom*));

        seedinlh1 = calloc(ctx->njob, sizeof(int));
        seedinlh2 = calloc(ctx->njob, sizeof(int));
    } else if (ctx->constraint && nseed == 0) {
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

    for (i = 0; i < ctx->njob - nadd; i++)
        alreadyaligned[i] = 1;
    for (i = ctx->njob - nadd; i < ctx->njob; i++)
        alreadyaligned[i] = 0;

    for (l = 0; l < ctx->njob; l++)
        fftlog[l] = 1;

    if (ctx->constraint && compacttree != 3) {
        if (specifictarget)
            calcimportance_target(ctx, ctx->njob, ntarget, effarr, aseq, localhomtable, targetmap, targetmapr, *alloclen);
        else
            calcimportance_half(ctx, ctx->njob, effarr, aseq, localhomtable, *alloclen);
    } else if (ctx->constraint && nseed) {
        dontcalcimportance_half(ctx, nseed, aseq, localhomtable);
    }

    tscore = 0.0;
    for (l = 0; l < ctx->njob - 1; l++) {
        m1 = topol[l][0][0];
        m2 = topol[l][1][0];

        if (effarr_kozo) {
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

        makedynamicmtx(ctx, dynamicmtx, ctx->n_dis_consweight_multi, dep[l].distfromtip);

        len1 = strlen(aseq[m1]);
        len2 = strlen(aseq[m2]);
        if (*alloclen < len1 + len2) {
            fprintf(stderr, "\nReallocating..");
            *alloclen = (len1 + len2) + 1000;
            ReallocateCharMtx(aseq, ctx->njob, *alloclen + 10);
            gaplen = realloc(gaplen, (*alloclen + 10) * sizeof(int));
            if (gaplen == NULL) {
                fprintf(stderr, "Cannot realloc gaplen\n");
                exit(1);
            }
            gapmap = realloc(gapmap, (*alloclen + 10) * sizeof(int));
            if (gapmap == NULL) {
                fprintf(stderr, "Cannot realloc gapmap\n");
                exit(1);
            }
            fprintf(stderr, "done. *alloclen = %d\n", *alloclen);
        }

        if (effarr_kozo) {
            clus1 = fastconjuction_noname_kozo(localmem[0], aseq, mseq1, effarr1, effarr, effarr1_kozo, effarr_kozo, indication1);
            clus2 = fastconjuction_noname_kozo(localmem[1], aseq, mseq2, effarr2, effarr, effarr2_kozo, effarr_kozo, indication2);
        } else {
            clus1 = fastconjuction_noname(localmem[0], aseq, mseq1, effarr1, effarr, indication1, ctx->minimumweight, &orieff1);  // orieff tsukau!
            clus2 = fastconjuction_noname(localmem[1], aseq, mseq2, effarr2, effarr, indication2, ctx->minimumweight, &orieff2);  // orieff tsukau!
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

        if (compacttree == 3)
            nfiles = nfilesfornode[l];
        else
            nfiles = 0;

        if (l < 1000 || l % 100 == 0)
            fprintf(stderr, "\rSTEP % 5d /%d ", l + 1, ctx->njob - 1);

#if REPORTCOSTS
        if (l < 1000 || l % 100 == 0)
            reporterr("\nclus1=%d, clus2=%d\n", clus1, clus2);
#endif

        if (ctx->constraint && compacttree != 3) {
            if (specifictarget)
                fastshrinklocalhom_target(localmem[0], localmem[1], localhomtable, localhomshrink, swaplist, targetmap);
            else
                fastshrinklocalhom_half(localmem[0], localmem[1], localhomtable, localhomshrink);
        } else if (ctx->constraint && nseed) {
            fastshrinklocalhom_half_seed(localmem[0], localmem[1], nseed, seedinlh1, seedinlh2, localhomtable, localhomshrink);
            for (i = 0; i < ctx->njob; i++)
                reporterr("seedinlh1[%d]=%d\n", i, seedinlh1[i]);
            for (i = 0; i < ctx->njob; i++)
                reporterr("seedinlh2[%d]=%d\n", i, seedinlh2[i]);
        }

        if (ctx->rnakozo && ctx->rnaprediction == 'm') {
            makegrouprna(grouprna1, singlerna, localmem[0]);
            makegrouprna(grouprna2, singlerna, localmem[1]);
        }

        if (!ctx->nevermemsave && (ctx->constraint != 2 && ctx->alg != 'M' && (len1 > 30000 || len2 > 30000))) {
            fprintf(stderr, "\nlen1=%d, len2=%d, Switching to the memsave mode.\n", len1, len2);
            ctx->alg = 'M';
            if (ctx->commonIP)
                FreeIntMtx(ctx->commonIP);
            ctx->commonIP = NULL;
            ctx->commonAlloc1 = 0;
            ctx->commonAlloc2 = 0;
        }

        if (fftlog[m1] && fftlog[m2])
            ffttry = (nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 1000 && clus2 < 1000);
        else
            ffttry = 0;
        if (ctx->constraint == 2) {
            if (ctx->alg == 'M') {
                fprintf(stderr, "\n\nMemory saving mode is not supported.\n\n");
                exit(1);
            }
            if (ctx->alg == 'A') {
                imp_match_init_strict(ctx, clus1, clus2, strlen(mseq1[0]), strlen(mseq2[0]), mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, swaplist, localmem[0], localmem[1], uselh, seedinlh1, seedinlh2, (compacttree == 3) ? l : -1, nfiles);
                if (ctx->rnakozo)
                    imp_rna(ctx, clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2);
#if REPORTCOSTS
//				reporterr(       "\n\n %d - %d (%d x %d) : \n", topol[l][0][0], topol[l][1][0], clus1, clus2 );
#endif
                pscore = A__align(ctx, dynamicmtx, ctx->penalty, ctx->penalty_ex, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, ctx->constraint, &dumdb, NULL, NULL, NULL, NULL, ctx->outgap, ctx->outgap, localmem[0][0], 1, cpmxchild0, cpmxchild1, cpmxhist + l, orieff1, orieff2);
            }
            if (ctx->alg == 'd') {
                imp_match_init_strictD(ctx, clus1, clus2, strlen(mseq1[0]), strlen(mseq2[0]), mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, swaplist, localmem[0], localmem[1], uselh, seedinlh1, seedinlh2, (compacttree == 3) ? l : -1, nfiles);
                if (ctx->rnakozo)
                    imp_rnaD(ctx, clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2);
                pscore = D__align(ctx, dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, ctx->constraint, &dumdb, ctx->outgap, ctx->outgap);
            } else if (ctx->alg == 'Q') {
                fprintf(stderr, "Not supported\n");
                exit(1);
            }
        } else if (ctx->force_fft || (ctx->use_fft && ffttry)) {
            fprintf(stderr, " f\b\b");
            if (ctx->alg == 'M') {
                fprintf(stderr, "m");
                pscore = Falign_udpari_long(ctx, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog + m1);
            } else
                pscore = Falign(ctx, NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog + m1);
        } else {
            fprintf(stderr, " d\b\b");
            fftlog[m1] = 0;
            switch (ctx->alg) {
                case ('a'):
                    pscore = Aalign(ctx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen);
                    break;
                case ('M'):
                    fprintf(stderr, "m");
                    pscore = MSalignmm(ctx, dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, NULL, NULL, NULL, ctx->outgap, ctx->outgap, cpmxchild0, cpmxchild1, cpmxhist + l, orieff1, orieff2);
                    break;
                case ('A'):
                    pscore = A__align(ctx, dynamicmtx, ctx->penalty, ctx->penalty_ex, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, 0, &dumdb, NULL, NULL, NULL, NULL, ctx->outgap, ctx->outgap, localmem[0][0], 1, cpmxchild0, cpmxchild1, cpmxhist + l, orieff1, orieff2);
                    break;
                case ('d'):
                    pscore = D__align(ctx, dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, 0, &dumdb, ctx->outgap, ctx->outgap);
                    break;
                default:
                    ErrorExit("ERROR IN SOURCE FILE");
            }
        }

        nlen[m1] = 0.5 * (nlen[m1] + nlen[m2]);

#if SCOREOUT
        fprintf(stderr, "score = %10.2f\n", pscore);
#endif
        tscore += pscore;

        if (ctx->disp)
            display(ctx, aseq, ctx->njob);

        if (mergeoralign[l] == '1') {
            reporterr("Check source!!\n");
            exit(1);
        }
        if (mergeoralign[l] == '2') {
            gapmaplen = strlen(mseq1[0]) - len1nocommongap + len1;
            adjustgapmap(gapmaplen, gapmap, mseq1[0]);
            if (opts->smoothing) {
                restorecommongapssmoothly(ctx->njob, ctx->njob - (clus1 + clus2), aseq, localmem[0], localmem[1], gapmap, *alloclen, '-');
                findnewgaps(ctx, 0, mseq1, gaplen);
                insertnewgaps_bothorders(ctx, ctx->njob, alreadyaligned, aseq, localmem[0], localmem[1], gaplen, gapmap, gapmaplen, *alloclen, ctx->alg, '-');
            } else {
                restorecommongaps(ctx->njob, ctx->njob - (clus1 + clus2), aseq, localmem[0], localmem[1], gapmap, *alloclen, '-');
                findnewgaps(ctx, 0, mseq1, gaplen);
                insertnewgaps(ctx, ctx->njob, alreadyaligned, aseq, localmem[0], localmem[1], gaplen, gapmap, *alloclen, ctx->alg, '-');
            }
            eq2dashmatometehayaku(mseq1, clus1);
            eq2dashmatometehayaku(mseq2, clus2);
            for (i = 0; (m = localmem[1][i]) > -1; i++)
                alreadyaligned[m] = 1;
        }

        free(localmem[0]);
        free(localmem[1]);
#if REPORTCOSTS
        if (l < 1000 || l % 100 == 0) {
            use_getrusage();
            reporterr("real = %f min\n", (float)(time(NULL) - starttime) / 60.0);
            reporterr("user = %f min\n", (float)(clock() - startclock) / CLOCKS_PER_SEC / 60);
        }
#endif
    }
#if REPORTCOSTS
    use_getrusage();
    reporterr("real = %f min\n", (float)(time(NULL) - starttime) / 60.0);
    reporterr("user = %f min\n", (float)(clock() - startclock) / CLOCKS_PER_SEC / 60);
#endif

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

    if (ctx->constraint) {
        if (localhomshrink)  // nen no tame
        {
            if (compacttree == 3)
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
        if (specifictarget)
            free(swaplist);
    }

    free(fftlog);
    free(effarr1);
    free(effarr2);
    free(indication1);
    free(indication2);
    free(gaplen);
    free(gapmap);
    free(alreadyaligned);
    FreeDoubleMtx(dynamicmtx);
    free(localmem);
    free(effarr1_kozo);
    free(effarr2_kozo);
    Falign(ctx, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL);
    Falign_udpari_long(ctx, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL);
    D__align(ctx, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, 0, 0);
    A__align(ctx, NULL, 0, 0, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, 0, 0, -1, -1, NULL, NULL, NULL, 0.0, 0.0);
    imp_match_init_strictD(ctx, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0);
    imp_match_init_strict(ctx, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0);
    FreeCommonIP(ctx);
}

static void
WriteOptions(Context* ctx, FILE* fp) {
    if (ctx->dorp == 'd')
        fprintf(fp, "DNA\n");
    else {
        if (ctx->scoremtx == 0)
            fprintf(fp, "JTT %dPAM\n", ctx->pamN);
        else if (ctx->scoremtx == 1)
            fprintf(fp, "BLOSUM %d\n", ctx->nblosum);
        else if (ctx->scoremtx == 2)
            fprintf(fp, "M-Y\n");
    }
    fprintf(stderr, "Gap Penalty = %+5.2f, %+5.2f, %+5.2f\n", (double)ctx->ppenalty / 1000, (double)ctx->ppenalty_ex / 1000, (double)ctx->poffset / 1000);
    if (ctx->use_fft)
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

    fprintf(fp, "Gap Penalty = %+5.2f, %+5.2f, %+5.2f\n", (double)ctx->ppenalty / 1000, (double)ctx->ppenalty_ex / 1000, (double)ctx->poffset / 1000);

    if (ctx->alg == 'a')
        fprintf(fp, "Algorithm A\n");
    else if (ctx->alg == 'A')
        fprintf(fp, "Algorithm A+\n");
    else if (ctx->alg == 'C')
        fprintf(fp, "Apgorithm A+/C\n");
    else
        fprintf(fp, "Unknown algorithm\n");

    if (ctx->treemethod == 'X')
        fprintf(fp, "Tree = UPGMA (mix).\n");
    else if (ctx->treemethod == 'E')
        fprintf(fp, "Tree = UPGMA (average).\n");
    else if (ctx->treemethod == 'q')
        fprintf(fp, "Tree = Minimum linkage.\n");
    else
        fprintf(fp, "Unknown tree.\n");

    if (ctx->use_fft) {
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

static double**
preparepartmtx(int nseq) {
    double** val = (double**)calloc(nseq, sizeof(double*));
    return val;
}

int
tbfast_main(int argc, char* argv[]) {
    Context* ctx = calloc(sizeof(Context), 1);
    ctx->dorp = NOTSPECIFIED;
    ctx->penalty_shift_factor = 100.0;
    ctx->outgap = 1;
    ctx->addprofile = 1;
    ctx->consweight_multi = 1.0;
    ctx->RNAscoremtx = 'n';
    ctx->nthread = 1;
    ctx->nthreadpair = 1;
    ctx->parallelizationstrategy = BAATARI1;
    ctx->minimumweight = 0.0005;
    ctx->newgapstr = "-";
    ctx->nalphabets = 26;
    ctx->nscoredalphabets = 20;
    ctx->specificityconsideration = 0.0;
    ctx->ndistclass = 10;
    ctx->maxdistclass = -1;
    ctx->sueff_global = SUEFF;

    TbfastOpts opts_ = {};
    TbfastOpts* opts = &opts_;

    int*     nlen = NULL;
    int*     selfscore = NULL;
    int      nogaplen;
    char **  name = NULL, **seq = NULL;
    char **  mseq1 = NULL, **mseq2 = NULL;
    char**   bseq = NULL;
    double **iscore = NULL, **iscore_kozo = NULL;
    int**    skiptable;
    double * eff = NULL, *eff_kozo = NULL, *eff_kozo_mapped = NULL;
    int      i, j, k, ien, ik, jk;
    int ***  topol = NULL, ***topol_kozo = NULL;
    double** expdist = NULL;
    int*     addmem;
    Treedep* dep = NULL;
    double **len = NULL, **len_kozo = NULL;
    FILE*    prep = NULL;
    FILE*    infp = NULL;
    FILE*    orderfp = NULL;
    FILE*    hat2p = NULL;
    double   unweightedspscore;
    size_t   alignmentlength;
    char*    mergeoralign = NULL;
    int      foundthebranch;
    int      nsubalignments, maxmem;
    int**    subtable;
    int*     insubtable;
    int*     preservegaps;
    char***  subalnpt;
    char*    originalgaps = NULL;
    char**   addbk = NULL;
    GapPos** deletelist = NULL;
    FILE*    dlf = NULL;
    int**    localmem = NULL;
    int      includememberres0, includememberres1;
    int*     mindistfrom = NULL;
    double*  mindist = NULL;
    double** partmtx = NULL;

    char         c;
    int          alloclen = 0;
    LocalHom**   localhomtable = NULL;
    LocalHom*    tmpptr;
    RNApair***   singlerna = NULL;
    double       ssi, ssj, bunbo;
    static char* kozoarivec = NULL;
    int          nkozo;
    int          ntarget;
    int *        targetmap = NULL, *targetmapr = NULL;
    int          ilim, jst, jj;

    FILE* fp;
    int*  uselh = NULL;
    int   nseed = 0;
    int*  nfilesfornode = NULL;

    char** pav = calloc(argc, sizeof(char*));
    char** tav = calloc(argc, sizeof(char*));

    int        pac = 0;
    int        tac = 0;
    arguments(ctx, opts, argc, argv, &pac, pav, &tac, tav);

    if (opts->treein) {
        int    dumx, dumy;
        double dumz;
        opts->treein = check_guidetreefile(&dumx, &dumy, &dumz);
        if (opts->treein == 'C') {
            compacttree = 2;
            opts->treein = 0;
            ctx->use_fft = 0;
        } else if (opts->treein == 'n') {
            compacttree = 3;
            opts->treein = 0;
            ctx->use_fft = 0;
        }
    }

    reporterr("treein = %d\n", opts->treein);
    reporterr("compacttree = %d\n", compacttree);

    if (ctx->fastathreshold < 0.0001)
        ctx->constraint = 0;

    if (ctx->inputfile) {
        infp = fopen(ctx->inputfile, "rb");
        if (!infp) {
            fprintf(stderr, "Cannot open %s\n", ctx->inputfile);
            exit(1);
        }
    } else
        infp = stdin;

    getnumlen(ctx, infp);
    rewind(infp);

    nkozo = 0;

#if !defined(mingw) && !defined(_MSC_VER)
    setstacksize(200 * ctx->njob);
#endif

    if (opts->subalignment) {
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
        subalnpt = AllocateCharCub(nsubalignments, maxmem, 0);
        readsubalignmentstable(ctx->njob, subtable, preservegaps, NULL, NULL);
    }

    seq = AllocateCharMtx(ctx->njob, ctx->nlenmax + 1);
    mseq1 = AllocateCharMtx(ctx->njob, 0);
    mseq2 = AllocateCharMtx(ctx->njob, 0);

    name = AllocateCharMtx(ctx->njob, B + 1);
    nlen = AllocateIntVec(ctx->njob);
    selfscore = AllocateIntVec(ctx->njob);

    topol = AllocateIntCub(ctx->njob, 2, 0);
    len = AllocateFloatMtx(ctx->njob, 2);
    eff = AllocateDoubleVec(ctx->njob);
    kozoarivec = AllocateCharVec(ctx->njob);

    mergeoralign = AllocateCharVec(ctx->njob);

    dep = (Treedep*)calloc(ctx->njob, sizeof(Treedep));
    if (nadd)
        addmem = AllocateIntVec(nadd + 1);
    localmem = AllocateIntMtx(2, ctx->njob + 1);

    if (compacttree == 3)
        nfilesfornode = calloc(sizeof(int), ctx->njob - 1);

#if REPORTCOSTS
    reporterr("before allocating iscore\n");
    use_getrusage();
#endif

    if (ctx->tbutree && compacttree != 3)
        iscore = AllocateFloatHalfMtx(ctx->njob);

    opts->ndeleted = 0;

    readData_pointer(ctx, infp, name, nlen, seq);
    fclose(infp);

    if (opts->treein) {
        loadtree(ctx, ctx->njob, topol, len, name, dep, opts->treeout);
        fprintf(stderr, "\ndone.\n\n");
        if (opts->callpairlocalalign && ctx->specificityconsideration > 0.0) {
            int* mem0 = calloc(sizeof(int), ctx->njob);
            int* mem1 = calloc(sizeof(int), ctx->njob);
            expdist = AllocateDoubleMtx(ctx->njob, ctx->njob);
            for (i = 0; i < ctx->njob - 1; i++) {
                topolorderz(mem0, topol, dep, i, 0);
                topolorderz(mem1, topol, dep, i, 1);

                for (j = 0; mem0[j] != -1; j++)
                    for (k = 0; mem1[k] != -1; k++) {
                        expdist[mem0[j]][mem1[k]] += (len[i][0] + len[i][1]);
                        expdist[mem1[k]][mem0[j]] += (len[i][0] + len[i][1]);
                    }
            }
            free(mem0);
            free(mem1);
        }
    }

    if (specifictarget && compacttree != 3) {
        targetmap = calloc(ctx->njob, sizeof(int));
        ntarget = 0;
        for (i = 0; i < ctx->njob; i++) {
            targetmap[i] = -1;
            if (!strncmp(name[i] + 1, "_focus_", 7))
                targetmap[i] = ntarget++;
        }
        targetmapr = calloc(ntarget, sizeof(int));
        for (i = 0; i < ctx->njob; i++)
            if (targetmap[i] != -1)
                targetmapr[targetmap[i]] = i;

    } else if (compacttree != 3) {
        ntarget = ctx->njob;
        targetmap = calloc(ctx->njob, sizeof(int));
        targetmapr = calloc(ctx->njob, sizeof(int));
        for (i = 0; i < ctx->njob; i++)
            targetmap[i] = targetmapr[i] = i;
    }

    if (ctx->constraint && compacttree != 3) {
        ilim = ctx->njob;
        localhomtable = (LocalHom**)calloc(ntarget, sizeof(LocalHom*));
        for (i = 0; i < ntarget; i++) {
            localhomtable[i] = (LocalHom*)calloc(ilim, sizeof(LocalHom));
            for (j = 0; j < ilim; j++) {
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
            if (!specifictarget)
                ilim--;
        }

        if (opts->callpairlocalalign) {
            pairlocalalign(ctx, ctx->njob, name, seq, iscore, localhomtable, pac, pav, expdist);
            arguments(ctx, opts, tac, tav, NULL, NULL, NULL, NULL);
            opts->callpairlocalalign = 1;
            if (expdist)
                FreeDoubleMtx(expdist);
            expdist = NULL;
            if (ctx->fastathreshold < 0.0001)
                ctx->constraint = 0;
            if (compacttree != 3) {
                for (ilim = ctx->njob, i = 0; i < ntarget; i++) {
                    for (j = 0; j < ilim; j++) {
                        for (tmpptr = localhomtable[i] + j; tmpptr; tmpptr = tmpptr->next) {
                            if (tmpptr->opt == -1.0)
                                continue;
#if SHISHAGONYU
                            char buff[100];
                            sprintf(buff, "%10.5f", tmpptr->opt);
                            tmpptr->opt = 0.0;
                            sscanf(buff, "%lf", &(tmpptr->opt));
#endif
                            tmpptr->opt = (tmpptr->opt) / 5.8 * 600;
                        }
                    }
                    if (!specifictarget)
                        ilim--;
                }

                prep = fopen("hat3.seed", "r");
                if (prep) {
                    fprintf(stderr, "Loading 'hat3.seed' ... ");
                    if (specifictarget)
                        readlocalhomtable2_target(prep, localhomtable, kozoarivec, targetmap);
                    else
                        readlocalhomtable2_half(prep, ctx->njob, localhomtable, kozoarivec);
                    fclose(prep);
                    fprintf(stderr, "\ndone.\n");
                } else
                    fprintf(stderr, "No hat3.seed. No problem.\n");

                if (opts->outputhat23) {
                    prep = fopen("hat3", "w");
                    if (!prep)
                        ErrorExit("Cannot open hat3 to write.");

                    fprintf(stderr, "Writing hat3 for iterative refinement\n");
                    if (specifictarget)
                        ilim = ntarget;
                    else
                        ilim = ctx->njob - 1;
                    for (i = 0; i < ilim; i++) {
                        if (specifictarget) {
                            jst = 0;
                            jj = 0;
                        } else {
                            jst = i;
                            jj = 0;
                        }
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
                } else if (opts->distout) {
                    prep = fopen("hat2", "w");
                    WriteFloatHat2_pointer_halfmtx(ctx, prep, ctx->njob, name, iscore);
                    fclose(prep);
                }
            } else {
                prep = fopen("hat3.seed", "r");
                if (prep) {
                    char r;
                    r = fgetc(prep);
                    if (isalnum(r) || r == ' ') {
                        reporterr("Structural alignment is not yet supported in the --memsavepair mode. Try normal mode,\n");
                        exit(1);
                    }
                    fclose(prep);
                }
            }
        } else {
            fprintf(stderr, "Loading 'hat3' ... ");
            prep = fopen("hat3", "r");
            if (prep == NULL)
                ErrorExit("Make hat3.");
            if (specifictarget)
                readlocalhomtable2_target(prep, localhomtable, kozoarivec, targetmap);
            else
                readlocalhomtable2_half(prep, ctx->njob, localhomtable, kozoarivec);
            fclose(prep);
            fprintf(stderr, "\ndone.\n");
        }

        nkozo = 0;
        for (i = 0; i < ctx->njob; i++) {
            //			fprintf( stderr, "kozoarivec[%d] = %d\n", i, kozoarivec[i] );
            if (kozoarivec[i])
                nkozo++;
        }
        if (nkozo) {
            topol_kozo = AllocateIntCub(nkozo, 2, 0);
            len_kozo = AllocateFloatMtx(nkozo, 2);
            iscore_kozo = AllocateFloatHalfMtx(nkozo);
            eff_kozo = AllocateDoubleVec(nkozo);
            eff_kozo_mapped = AllocateDoubleVec(ctx->njob);
        }
    } else if (compacttree != 3) {
        if (opts->callpairlocalalign) {
            pairlocalalign(ctx, ctx->njob, name, seq, iscore, NULL, pac, pav, expdist);
            arguments(ctx, opts, tac, tav, NULL, NULL, NULL, NULL);
            opts->callpairlocalalign = 1;
            if (expdist)
                FreeDoubleMtx(expdist);
            expdist = NULL;
            if (ctx->fastathreshold < 0.0001)
                ctx->constraint = 0;
            fprintf(stderr, "blosum %d / kimura 200\n", ctx->nblosum);
            fprintf(stderr, "scoremtx=%d\n", ctx->scoremtx);
            fprintf(stderr, "fastathreshold=%f\n", ctx->fastathreshold);
        }
        if (opts->distout || opts->outputhat23) {
            reporterr("\nwriting hat2 (1)\n");
            prep = fopen("hat2", "w");
            WriteFloatHat2_pointer_halfmtx(ctx, prep, ctx->njob, name, iscore);
            fclose(prep);
        }
    } else  // ie, conpacttree == 3 // ntarget ha tsukawanai. uselh <- nodepair kara
    {
        specifictarget = 0;  // ichiou uwagaki. '-=' option ha mushi sareru.
        //		reporterr( "compacttree=3, ntarget=%d\n", ntarget );
        // hat3 <- hat3dir
        // hat3.seed

        prep = fopen("hat3.seed", "r");
        if (prep) {
            char r;
            fprintf(stderr, "Checking 'hat3.seed' ... ");
            r = fgetc(prep);
            if (isalnum(r) || r == ' ')
                nkozo = 1;
            else {
                nkozo = 0;
                fclose(prep);
            }
        } else {
            nkozo = 0;
            fprintf(stderr, "No hat3.seed. No problem.\n");
            localhomtable = NULL;
        }

        nseed = 0;
        if (nkozo) {
            for (i = 0; i < ctx->njob; i++)
                if (strncmp(name[i] + 1, "_seed", 4))
                    break;  // konran!!!
            nseed = i;

            topol_kozo = AllocateIntCub(nseed, 2, 0);
            len_kozo = AllocateFloatMtx(nseed, 2);
            iscore_kozo = AllocateFloatHalfMtx(nseed);
            eff_kozo = AllocateDoubleVec(nseed);
            eff_kozo_mapped = AllocateDoubleVec(ctx->njob);

            if (localhomtable) {
                reporterr("bug. localhomtable is already allocated?\n");
                exit(1);
            }

            ilim = nseed;
            localhomtable = (LocalHom**)calloc(nseed, sizeof(LocalHom*));
            for (i = 0; i < nseed; i++) {
                localhomtable[i] = (LocalHom*)calloc(ilim, sizeof(LocalHom));
                for (j = 0; j < ilim; j++) {
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
                    localhomtable[i][j].korh = 'k';  // dochirademo
                }
                ilim--;
            }
            readlocalhomtable2_half(prep, nseed, localhomtable, kozoarivec);
            fclose(prep);

            nkozo = 0;  // kazoenaoshi
            for (i = 0; i < ctx->njob; i++) {
                fprintf(stderr, "kozoarivec[%d] = %d\n", i, kozoarivec[i]);
                if (kozoarivec[i])
                    nkozo++;
            }
            if (nkozo != nseed) {
                reporterr("problem in input file?  nkozo != nseed\n");
                exit(1);
            }
        }
    }

    free(tav);
    free(pav);
    constants(ctx, ctx->njob, seq);

    if (ctx->sueff_global < 0.0001 || compacttree == 3) {
        ctx->nthreadreadlh = ctx->nthread;
        if (ctx->nthreadreadlh == 0)
            ctx->nthreadreadlh = 1;

        opts->nthreadtb = 0;
        ctx->nthread = 0;
    } else {
        ctx->nthreadreadlh = 1;
        opts->nthreadtb = ctx->nthread;
    }

    initSignalSM(ctx);
    initFiles(ctx);
    WriteOptions(ctx, ctx->trap_g);

    if (opts->distout && !opts->treeout && opts->noalign) {
        writeData_pointer(ctx->prep_g, ctx->njob, name, seq);
        fprintf(stderr, "\n");
        goto chudan;
    }

    c = seqcheck(ctx, seq);
    if (c) {
        fprintf(stderr, "Illegal character %c\n", c);
        exit(1);
    }

    if (nadd && opts->keeplength) {
        originalgaps = (char*)calloc(ctx->nlenmax + 1, sizeof(char));
        recordoriginalgaps(originalgaps, ctx->njob - nadd, seq);

        if (opts->mapout) {
            addbk = (char**)calloc(nadd + 1, sizeof(char*));
            for (i = 0; i < nadd; i++) {
                ien = strlen(seq[ctx->njob - nadd + i]);
                addbk[i] = (char*)calloc(ien + 1, sizeof(char));
                gappick0(addbk[i], seq[ctx->njob - nadd + i]);
            }
            addbk[nadd] = NULL;
        } else
            addbk = NULL;
    } else {
        originalgaps = NULL;
        addbk = NULL;
    }

    if (!opts->treein) {
        reporterr("tbutree = %d, compacttree = %d\n", ctx->tbutree, compacttree);
        if (compacttree == 3) {
            iscore = NULL;
        } else if (ctx->tbutree == 0 && compacttree) {
            iscore = NULL;
            reporterr("Making a compact tree from msa, step 1.. \n");
            skiptable = AllocateIntMtx(ctx->njob, 0);
            makeskiptable(ctx->njob, skiptable, seq);
            mindistfrom = (int*)calloc(ctx->njob, sizeof(int));
            mindist = (double*)calloc(ctx->njob, sizeof(double));
            partmtx = preparepartmtx(ctx->njob);

            for (i = 0; i < ctx->njob; i++) {
                selfscore[i] = (int)naivepairscorefast(ctx, seq[i], seq[i], skiptable[i], skiptable[i], ctx->penalty_dist);
            }

            {
                int jobpos = 0;

                {
                    for (j = 0; j < ctx->njob; j++) {
                        mindist[j] = 999.9;
                        mindistfrom[j] = -1;
                    }

                    msacompactdistmtxthread_arg_t targ = {
                        .ctx = ctx,
                        .thread_no = 0,
                        .njob = ctx->njob,
                        .selfscore = selfscore,
                        .partmtx = partmtx,
                        .seq = seq,
                        .skiptable = skiptable,
                        .jobpospt = &jobpos,
                        .mindistfrom = mindistfrom,
                        .mindist = mindist,
                    };

                    msacompactdisthalfmtxthread(&targ);
                }

                for (i = 0; i < ctx->njob; i++)
                    mindist[i] -= preferenceval(i, mindistfrom[i], ctx->njob);
            }
            reporterr("\rdone.                                          \n");
        } else if (ctx->tbutree == 0 && compacttree == 0) {
            reporterr("Making a distance matrix from msa .. \n");
            iscore = AllocateFloatHalfMtx(ctx->njob);

            for (i = 1; i < ctx->njob; i++) {
                if (nlen[i] != nlen[0]) {
                    fprintf(stderr, "Input pre-aligned seqences or make hat2.\n");
                    exit(1);
                }
            }

            skiptable = AllocateIntMtx(ctx->njob, 0);
            makeskiptable(ctx->njob, skiptable, seq);
            ien = ctx->njob - 1;
            for (i = 0; i < ctx->njob; i++) {
                selfscore[i] = (int)naivepairscorefast(ctx, seq[i], seq[i], skiptable[i], skiptable[i], ctx->penalty_dist);
            }

            {
                for (i = 0; i < ien; i++) {
                    if (i % 10 == 0) {
                        fprintf(stderr, "\r% 5d / %d", i, ien);
                    }
                    ssi = selfscore[i];
                    for (j = i + 1; j < ctx->njob; j++) {
                        ssj = selfscore[j];
                        bunbo = MIN(ssi, ssj);
                        if (bunbo == 0.0)
                            iscore[i][j - i] = 2.0;
                        else
                            iscore[i][j - i] = (1.0 - naivepairscorefast(ctx, seq[i], seq[j], skiptable[i], skiptable[j], ctx->penalty_dist) / bunbo) * 2.0;  // 2014/Aug/15 fast
                        if (iscore[i][j - i] > 10)
                            iscore[i][j - i] = 10.0;
                    }
                }
            }
            //			fprintf( stderr, "\ndone.\n\n" );
            FreeIntMtx(skiptable);
            //			fflush( stderr );
            reporterr("\rdone.                                           \n");

        } else {
            if (opts->callpairlocalalign) {
                if (opts->multidist) {
                    reporterr("Bug in v7.290.  Please email katoh@ifrec.osaka-u.ac.jp\n");
                    exit(1);
                }
            } else {
                if (opts->multidist) {
                    fprintf(stderr, "Loading 'hat2n' (aligned sequences - new sequences) ... ");
                    prep = fopen("hat2n", "r");
                    if (prep == NULL)
                        ErrorExit("Make hat2.");
                    readhat2_doublehalf_pointer(prep, ctx->njob, iscore);
                    fclose(prep);
                    fprintf(stderr, "done.\n");

                    fprintf(stderr, "Loading 'hat2i' (aligned sequences) ... ");
                    prep = fopen("hat2i", "r");
                    if (prep == NULL)
                        ErrorExit("Make hat2i.");
                    readhat2_doublehalf_pointer(prep, ctx->njob - nadd, iscore);
                    fclose(prep);
                    fprintf(stderr, "done.\n");
                } else {
                    fprintf(stderr, "Loading 'hat2' ... ");
                    prep = fopen("hat2", "r");
                    if (prep == NULL)
                        ErrorExit("Make hat2.");
                    readhat2_doublehalf_pointer(prep, ctx->njob, iscore);
                    fclose(prep);
                    fprintf(stderr, "done.\n");
                }

                if (opts->distout) {
                    reporterr("\nwriting hat2 (2)\n");
                    hat2p = fopen("hat2", "w");
                    WriteFloatHat2_pointer_halfmtx(ctx, hat2p, ctx->njob, name, iscore);
                    fclose(hat2p);
                }
            }
        }

        if (nkozo && compacttree != 3) {
            ien = ctx->njob - 1;
            ik = 0;
            for (i = 0; i < ien; i++) {
                jk = ik + 1;
                for (j = i + 1; j < ctx->njob; j++) {
                    if (kozoarivec[i] && kozoarivec[j]) {
                        iscore_kozo[ik][jk - ik] = iscore[i][j - i];
                    }
                    if (kozoarivec[j])
                        jk++;
                }
                if (kozoarivec[i])
                    ik++;
            }
        } else if (nkozo && compacttree == 3)  // kozo ha saisho ni atsumarunode ik, jk ha hontouha iranai.
        {
            for (i = 0; i < ctx->njob; i++) {
                if (kozoarivec[i])
                    selfscore[i] = naivepairscore11(ctx, seq[i], seq[i], 0.0);
                else
                    selfscore[i] = -1;
            }
            ien = ctx->njob - 1;
            ik = 0;
            for (i = 0; i < ien; i++) {
                jk = ik + 1;
                for (j = i + 1; j < ctx->njob; j++) {
                    if (kozoarivec[i] && kozoarivec[j]) {
                        reporterr("seq0=%s\n", seq[i]);
                        reporterr("seq1=%s\n", seq[j]);
                        reporterr("selfscore0=%d\n", selfscore[0]);
                        reporterr("selfscore1=%d\n", selfscore[1]);
                        iscore_kozo[ik][jk - ik] = distdp_noalign(ctx, seq[i], seq[j], (double)selfscore[i], (double)selfscore[j], alloclen);
                        reporterr("iscore_kozo[%d][%d]=%f\n", ik, jk, iscore_kozo[ik][jk - ik]);
                    }
                    if (kozoarivec[j])
                        jk++;
                }
                if (kozoarivec[i])
                    ik++;
                G__align11_noalign(ctx, NULL, 0, 0, NULL, NULL);
            }
        }

        if (opts->subalignment) {
            fprintf(stderr, "Constructing a UPGMA tree ... ");
            fixed_supg_double_realloc_nobk_halfmtx_treeout_constrained(ctx, ctx->njob, iscore, topol, len, name, dep, nsubalignments, subtable, 1);
        } else if (compacttree == 3) {
            fp = fopen("hat3dir/tree", "rb");  // windows no tame rb
            treein_bin(fp, ctx->njob, topol, len, dep, nfilesfornode);
            fclose(fp);

            if (ctx->constraint) {
                uselh = AllocateIntVec(ctx->njob);
                fp = fopen("hat3dir/uselh", "rb");
                if (uselhin(fp, ctx->njob, uselh)) {
                    free(uselh);
                    uselh = NULL;
                }
                fclose(fp);
                //				for( i=0; i<ctx->njob; i++ ) reporterr( "uselh[%d]=%d\n", i, uselh[i] );
            }
        } else if (ctx->tbutree == 0 && compacttree)  // tbutree != 0 no toki (aln->mtx) ha, 6merdistance -> disttbfast.c; dp distance -> muzukashii
        {
            reporterr("Constructing a tree ... nthread=%d", ctx->nthread);
            compacttree_memsaveselectable(ctx, ctx->njob, partmtx, mindistfrom, mindist, NULL, selfscore, seq, skiptable, topol, len, name, NULL, dep, opts->treeout, compacttree, 1);

            if (mindistfrom)
                free(mindistfrom);
            mindistfrom = NULL;
            if (mindist)
                free(mindist);
            mindist = NULL;
            if (skiptable)
                FreeIntMtx(skiptable);
            skiptable = NULL;
            free(partmtx);
        } else if (opts->treeout) {
            fprintf(stderr, "Constructing a UPGMA tree ... ");
            fixed_musclesupg_double_realloc_nobk_halfmtx_treeout_memsave(ctx, ctx->njob, iscore, topol, len, name, dep, 1, opts->treeout);
        } else {
            fprintf(stderr, "Constructing a UPGMA tree ... ");
            fixed_musclesupg_double_realloc_nobk_halfmtx_memsave(ctx, ctx->njob, iscore, topol, len, dep, 1, 1);
        }

        if (nkozo) {
            fixed_musclesupg_double_realloc_nobk_halfmtx(ctx, nkozo, iscore_kozo, topol_kozo, len_kozo, NULL, 1, 1);
        }
        fprintf(stderr, "\ndone.\n\n");
    }

    localmem[0][0] = -1;
    topolorderz(localmem[0], topol, dep, ctx->njob - 2, 2);

    orderfp = fopen("order", "w");
    if (!orderfp) {
        fprintf(stderr, "Cannot open 'order'\n");
        exit(1);
    }

    for (i = 0; i < ctx->njob; i++)
        fprintf(orderfp, "%d\n", localmem[0][i]);

    fclose(orderfp);

    if (opts->treeout && opts->noalign) {
        writeData_pointer(ctx->prep_g, ctx->njob, name, seq);
        fprintf(stderr, "\n");
        goto chudan;  // 2016Jul31
    }

    if (ctx->tbrweight) {
        ctx->weight = 3;
        counteff_simple_double_nostatic_memsave(ctx->njob, topol, len, dep, eff);
        for (i = ctx->njob - nadd; i < ctx->njob; i++)
            eff[i] /= (double)100;
        if (nkozo) {
            for (i = 0, j = 0; i < ctx->njob; i++) {
                if (kozoarivec[i]) {
                    eff_kozo_mapped[i] = eff[i];  // single weight
                    j++;
                } else
                    eff_kozo_mapped[i] = 0.0;
            }
        }
    } else {
        for (i = 0; i < ctx->njob; i++)
            eff[i] = 1.0;
        if (nkozo) {
            for (i = 0; i < ctx->njob; i++) {
                if (kozoarivec[i])
                    eff_kozo_mapped[i] = 1.0;
                else
                    eff_kozo_mapped[i] = 0.0;
            }
        }
    }

    if (iscore)
        FreeFloatHalfMtx(iscore, ctx->njob);
    iscore = NULL;
    if (iscore_kozo)
        FreeFloatHalfMtx(iscore_kozo, nkozo);
    iscore = NULL;
    if (topol_kozo)
        FreeIntCub(topol_kozo);
    topol_kozo = NULL;
    if (len_kozo)
        FreeFloatMtx(len_kozo);
    len_kozo = NULL;
    if (eff_kozo)
        free(eff_kozo);
    eff_kozo = NULL;
    FreeFloatMtx(len);

    alloclen = ctx->nlenmax * 2 + 1;  //chuui!
    bseq = AllocateCharMtx(ctx->njob, alloclen);

    if (nadd) {
        alignmentlength = strlen(seq[0]);
        for (i = 0; i < ctx->njob - nadd; i++) {
            if (alignmentlength != strlen(seq[i])) {
                fprintf(stderr, "#################################################################################\n");
                fprintf(stderr, "# ERROR!                                                                        #\n");
                fprintf(stderr, "# The original%4d sequences must be aligned                                    #\n", ctx->njob - nadd);
                fprintf(stderr, "#################################################################################\n");
                exit(1);
            }
        }
        if (ctx->addprofile) {
            alignmentlength = strlen(seq[ctx->njob - nadd]);
            for (i = ctx->njob - nadd; i < ctx->njob; i++) {
                if (alignmentlength != strlen(seq[i])) {
                    fprintf(stderr, "###############################################################################\n");
                    fprintf(stderr, "# ERROR!                                                                      #\n");
                    fprintf(stderr, "# The%4d additional sequences must be aligned                                #\n", nadd);
                    fprintf(stderr, "# Otherwise, try the '--add' option, instead of '--addprofile' option.        #\n");
                    fprintf(stderr, "###############################################################################\n");
                    exit(1);
                }
            }
            for (i = 0; i < nadd; i++)
                addmem[i] = ctx->njob - nadd + i;
            addmem[nadd] = -1;
            foundthebranch = 0;
            for (i = 0; i < ctx->njob - 1; i++) {
                localmem[0][0] = -1;
                topolorderz(localmem[0], topol, dep, i, 0);
                localmem[1][0] = -1;
                topolorderz(localmem[1], topol, dep, i, 1);

                if (samemember(localmem[0], addmem)) {
                    mergeoralign[i] = '1';
                    foundthebranch = 1;
                } else if (samemember(localmem[1], addmem)) {
                    mergeoralign[i] = '2';
                    foundthebranch = 1;
                } else {
                    mergeoralign[i] = 'n';
                }
            }
            if (!foundthebranch) {
                fprintf(stderr, "###############################################################################\n");
                fprintf(stderr, "# ERROR!                                                                      #\n");
                fprintf(stderr, "# There is no appropriate position to add the%4d sequences in the guide tree.#\n", nadd);
                fprintf(stderr, "# Check whether the%4d sequences form a monophyletic cluster.                #\n", nadd);
                fprintf(stderr, "# If not, try the '--add' option, instead of the '--addprofile' option.       #\n");
                fprintf(stderr, "############################################################################### \n");
                exit(1);
            }
            commongappick(nadd, seq + ctx->njob - nadd);
            for (i = ctx->njob - nadd; i < ctx->njob; i++)
                strcpy(bseq[i], seq[i]);
        } else {
            for (i = 0; i < ctx->njob - 1; i++)
                mergeoralign[i] = 'n';
            for (i = 0; i < nadd; i++)
                addmem[i] = ctx->njob - nadd + i;
            addmem[nadd] = -1;
            for (i = 0; i < ctx->njob - 1; i++) {
                localmem[0][0] = -1;
                topolorderz(localmem[0], topol, dep, i, 0);
                localmem[1][0] = -1;
                topolorderz(localmem[1], topol, dep, i, 1);

                includememberres0 = includemember(localmem[0], addmem);
                includememberres1 = includemember(localmem[1], addmem);
                if (includememberres0 && includememberres1) {
                    mergeoralign[i] = 'w';
                } else if (includememberres0) {
                    mergeoralign[i] = '1';
                } else if (includememberres1) {
                    mergeoralign[i] = '2';
                }
            }

            for (i = ctx->njob - nadd; i < ctx->njob; i++)
                gappick0(bseq[i], seq[i]);
        }

        commongappick(ctx->njob - nadd, seq);
        for (i = 0; i < ctx->njob - nadd; i++)
            strcpy(bseq[i], seq[i]);
    } else if (opts->subalignment) {
        for (i = 0; i < ctx->njob - 1; i++)
            mergeoralign[i] = 'a';
        for (i = 0; i < nsubalignments; i++) {
            fprintf(stderr, "Checking subalignment %d:\n", i + 1);
            alignmentlength = strlen(seq[subtable[i][0]]);
            //			for( j=0; subtable[i][j]!=-1; j++ )
            //				fprintf( stderr, " %d. %-30.30s\n", subtable[i][j]+1, name[subtable[i][j]]+1 );
            for (j = 0; subtable[i][j] != -1; j++) {
                if (subtable[i][j] >= ctx->njob) {
                    fprintf(stderr, "No such sequence, %d.\n", subtable[i][j] + 1);
                    exit(1);
                }
                if (alignmentlength != strlen(seq[subtable[i][j]])) {
                    fprintf(stderr, "\n");
                    fprintf(stderr, "###############################################################################\n");
                    fprintf(stderr, "# ERROR!\n");
                    fprintf(stderr, "# Subalignment %d must be aligned.\n", i + 1);
                    fprintf(stderr, "# Please check the alignment lengths of following sequences.\n");
                    fprintf(stderr, "#\n");
                    fprintf(stderr, "# %d. %-10.10s -> %zu letters (including gaps)\n", subtable[i][0] + 1, name[subtable[i][0]] + 1, alignmentlength);
                    fprintf(stderr, "# %d. %-10.10s -> %d letters (including gaps)\n", subtable[i][j] + 1, name[subtable[i][j]] + 1, (int)strlen(seq[subtable[i][j]]));
                    fprintf(stderr, "#\n");
                    fprintf(stderr, "# See http://mafft.cbrc.jp/alignment/software/merge.html for details.\n");
                    if (opts->subalignmentoffset) {
                        fprintf(stderr, "#\n");
                        fprintf(stderr, "# You specified seed alignment(s) consisting of %d sequences.\n", opts->subalignmentoffset);
                        fprintf(stderr, "# In this case, the rule of numbering is:\n");
                        fprintf(stderr, "#   The aligned seed sequences are numbered as 1 .. %d\n", opts->subalignmentoffset);
                        fprintf(stderr, "#   The input sequences to be aligned are numbered as %d .. %d\n", opts->subalignmentoffset + 1, opts->subalignmentoffset + ctx->njob);
                    }
                    fprintf(stderr, "###############################################################################\n");
                    fprintf(stderr, "\n");
                    exit(1);
                }
                insubtable[subtable[i][j]] = 1;
            }
            for (j = 0; j < ctx->njob - 1; j++) {
                if (includemember(topol[j][0], subtable[i]) && includemember(topol[j][1], subtable[i])) {
                    mergeoralign[j] = 'n';
                }
            }
            foundthebranch = 0;
            for (j = 0; j < ctx->njob - 1; j++) {
                if (samemember(topol[j][0], subtable[i]) || samemember(topol[j][1], subtable[i])) {
                    foundthebranch = 1;
                    fprintf(stderr, " -> OK\n");
                    break;
                }
            }
            if (!foundthebranch) {
                system("cp infile.tree GuideTree");  // tekitou
                fprintf(stderr, "\n");
                fprintf(stderr, "###############################################################################\n");
                fprintf(stderr, "# ERROR!\n");
                fprintf(stderr, "# Subalignment %d does not form a monophyletic cluster\n", i + 1);
                fprintf(stderr, "# in the guide tree ('GuideTree' in this directory) internally computed.\n");
                fprintf(stderr, "# If you really want to use this subalignment, pelase give a tree with --treein \n");
                fprintf(stderr, "# http://mafft.cbrc.jp/alignment/software/treein.html\n");
                fprintf(stderr, "# http://mafft.cbrc.jp/alignment/software/merge.html\n");
                if (opts->subalignmentoffset) {
                    fprintf(stderr, "#\n");
                    fprintf(stderr, "# You specified seed alignment(s) consisting of %d sequences.\n", opts->subalignmentoffset);
                    fprintf(stderr, "# In this case, the rule of numbering is:\n");
                    fprintf(stderr, "#   The aligned seed sequences are numbered as 1 .. %d\n", opts->subalignmentoffset);
                    fprintf(stderr, "#   The input sequences to be aligned are numbered as %d .. %d\n", opts->subalignmentoffset + 1, opts->subalignmentoffset + ctx->njob);
                }
                fprintf(stderr, "############################################################################### \n");
                fprintf(stderr, "\n");
                exit(1);
            }
        }

        for (i = 0; i < ctx->njob; i++) {
            if (insubtable[i])
                strcpy(bseq[i], seq[i]);
            else
                gappick0(bseq[i], seq[i]);
        }

        for (i = 0; i < nsubalignments; i++) {
            for (j = 0; subtable[i][j] != -1; j++)
                subalnpt[i][j] = bseq[subtable[i][j]];
            if (!preservegaps[i])
                commongappick(j, subalnpt[i]);
        }

        FreeIntMtx(subtable);
        free(insubtable);
        for (i = 0; i < nsubalignments; i++)
            free(subalnpt[i]);
        free(subalnpt);
        free(preservegaps);
    } else {
        for (i = 0; i < ctx->njob; i++)
            gappick0(bseq[i], seq[i]);
        for (i = 0; i < ctx->njob - 1; i++)
            mergeoralign[i] = 'a';
    }

    if (ctx->rnakozo && ctx->rnaprediction == 'm') {
        singlerna = (RNApair***)calloc(ctx->njob, sizeof(RNApair**));
        prep = fopen("hat4", "r");
        if (prep == NULL)
            ErrorExit("Make hat4 using mccaskill.");
        fprintf(stderr, "Loading 'hat4' ... ");
        for (i = 0; i < ctx->njob; i++) {
            nogaplen = strlen(bseq[i]);
            singlerna[i] = (RNApair**)calloc(nogaplen + 1, sizeof(RNApair*));
            for (j = 0; j < nogaplen; j++) {
                singlerna[i][j] = (RNApair*)calloc(1, sizeof(RNApair));
                singlerna[i][j][0].bestpos = -1;
                singlerna[i][j][0].bestscore = -1.0;
            }
            singlerna[i][nogaplen] = NULL;
            readmccaskill(prep, singlerna[i], nogaplen);
        }
        fclose(prep);
        fprintf(stderr, "\ndone.\n");
    } else
        singlerna = NULL;

    fprintf(stderr, "Progressive alignment ... \n");

    treebase(ctx, opts, nlen, bseq, nadd, mergeoralign, mseq1, mseq2, topol, dep, eff, &alloclen, localhomtable, singlerna, eff_kozo_mapped, targetmap, targetmapr, ntarget, uselh, nseed, nfilesfornode);

    fprintf(stderr, "\ndone.\n");

    if (opts->keeplength) {
        dlf = fopen("_deletelist", "w");
        deletelist = (GapPos**)calloc(nadd + 1, sizeof(GapPos*));
        for (i = 0; i < nadd; i++) {
            deletelist[i] = calloc(1, sizeof(GapPos));
            deletelist[i][0].pos = -1;
            deletelist[i][0].len = 0;
        }
        deletelist[nadd] = NULL;
        opts->ndeleted = deletenewinsertions_whole(ctx->njob - nadd, nadd, bseq, bseq + ctx->njob - nadd, deletelist);

        for (i = 0; i < nadd; i++) {
            if (deletelist[i])
                for (j = 0; deletelist[i][j].pos != -1; j++)
                    fprintf(dlf, "%d %d %d\n", ctx->njob - nadd + i, deletelist[i][j].pos, deletelist[i][j].len);  // 0origin
        }
        fclose(dlf);

        restoreoriginalgaps(ctx->njob, bseq, originalgaps);
        free(originalgaps);
        originalgaps = NULL;

        if (opts->mapout) {
            dlf = fopen("_deletemap", "w");
            if (opts->mapout == 1)
                reconstructdeletemap(ctx, nadd, addbk, deletelist, bseq + ctx->njob - nadd, dlf, name + ctx->njob - nadd);
            else
                reconstructdeletemap_compact(ctx, nadd, addbk, deletelist, seq + ctx->njob - nadd, dlf, name + ctx->njob - nadd);
            FreeCharMtx(addbk);
            addbk = NULL;
            fclose(dlf);
        }

        for (i = 0; deletelist[i] != NULL; i++)
            free(deletelist[i]);
        free(deletelist);
        deletelist = NULL;
    }

    if (ctx->scoreout) {
        unweightedspscore = plainscore(ctx, ctx->njob, bseq);
        fprintf(stderr, "\nSCORE %s = %.0f, ", "(treebase)", unweightedspscore);
        fprintf(stderr, "SCORE / residue = %f", unweightedspscore / (ctx->njob * strlen(bseq[0])));
        fprintf(stderr, "\n\n");
    }

    fprintf(ctx->trap_g, "done.\n");
    free(mergeoralign);
    freeconstants(ctx);

    if (ctx->rnakozo && ctx->rnaprediction == 'm') {
        if (singlerna)  // nen no tame
        {
            for (i = 0; i < ctx->njob; i++) {
                for (j = 0; singlerna[i][j] != NULL; j++) {
                    if (singlerna[i][j])
                        free(singlerna[i][j]);
                }
                if (singlerna[i])
                    free(singlerna[i]);
            }
            free(singlerna);
            singlerna = NULL;
        }
    }

    writeData_pointer(ctx->prep_g, ctx->njob, name, bseq);

#if IODEBUG
    fprintf(stderr, "OSHIMAI\n");
#endif

    if (ctx->constraint && compacttree != 3) {
        if (specifictarget)
            FreeLocalHomTable_part(localhomtable, ntarget, ctx->njob);
        else
            FreeLocalHomTable_half(localhomtable, ctx->njob);
    } else if (ctx->constraint && nkozo) {
        FreeLocalHomTable_half(localhomtable, nkozo);
    }

    if (compacttree != 3) {
        free(targetmap);
        free(targetmapr);
    }

    if (ctx->constraint && compacttree == 3 && uselh)
        free(uselh);
    uselh = NULL;

    if (ctx->spscoreout)
        reporterr("Unweighted sum-of-pairs score = %10.5f\n", sumofpairsscore(ctx, ctx->njob, bseq));
    ctx->nthread = MAX(opts->nthreadtb, ctx->nthreadreadlh);  // toriaezu
    if (opts->ndeleted > 0) {
        reporterr("\nTo keep the alignment length, %d letters were DELETED.\n", opts->ndeleted);
        if (opts->mapout)
            reporterr("The deleted letters are shown in the (filename).map file.\n");
        else
            reporterr("To know the positions of deleted letters, rerun the same command with the --mapout option.\n");
    }

    free(kozoarivec);
    FreeCharMtx(seq);
    FreeCharMtx(bseq);
    free(mseq1);
    free(mseq2);

    FreeCharMtx(name);
    free(nlen);
    free(selfscore);

    FreeIntCub(topol);
    topol = NULL;
    free(eff);
    free(dep);
    if (nfilesfornode)
        free(nfilesfornode);
    nfilesfornode = NULL;
    closeFiles(ctx);
    if (nadd)
        free(addmem);
    FreeIntMtx(localmem);
    if (eff_kozo_mapped)
        free(eff_kozo_mapped);
    eff_kozo_mapped = NULL;

    return (0);

chudan:
    if (seq)
        FreeCharMtx(seq);
    seq = NULL;
    if (mseq1)
        free(mseq1);
    mseq1 = NULL;
    if (mseq2)
        free(mseq2);
    mseq2 = NULL;

    if (name)
        FreeCharMtx(name);
    name = NULL;
    if (nlen)
        free(nlen);
    nlen = NULL;
    if (selfscore)
        free(selfscore);
    selfscore = NULL;
    if (mergeoralign)
        free(mergeoralign);
    mergeoralign = NULL;

    if (localhomtable) {
        reporterr("freeing localhomtable\n");
        if (specifictarget)
            FreeLocalHomTable_part(localhomtable, ntarget, ctx->njob);
        else
            FreeLocalHomTable_half(localhomtable, ctx->njob);
    }
    localhomtable = NULL;
    if (targetmap)
        free(targetmap);
    targetmap = NULL;
    if (targetmapr)
        free(targetmapr);
    targetmapr = NULL;
    if (uselh)
        free(uselh);
    uselh = NULL;

    if (kozoarivec)
        free(kozoarivec);
    kozoarivec = NULL;
    if (eff_kozo_mapped)
        free(eff_kozo_mapped);
    eff_kozo_mapped = NULL;
    if (eff_kozo)
        free(eff_kozo);
    eff_kozo = NULL;

    if (topol)
        FreeIntCub(topol);
    topol = NULL;
    if (topol_kozo)
        FreeIntCub(topol_kozo);
    topol_kozo = NULL;
    if (len)
        FreeFloatMtx(len);
    len = NULL;
    if (len_kozo)
        FreeFloatMtx(len_kozo);
    len_kozo = NULL;
    if (iscore)
        FreeFloatHalfMtx(iscore, ctx->njob);
    iscore = NULL;
    if (iscore_kozo && nseed)
        FreeFloatHalfMtx(iscore_kozo, nseed);
    iscore = NULL;  // ? nseed?
    if (eff)
        free(eff);
    eff = NULL;
    if (dep)
        free(dep);
    dep = NULL;
    if (nfilesfornode)
        free(nfilesfornode);
    nfilesfornode = NULL;

    if (addmem)
        free(addmem);
    addmem = NULL;
    if (localmem)
        FreeIntMtx(localmem);
    localmem = NULL;

    freeconstants(ctx);
    closeFiles(ctx);
    FreeCommonIP(ctx);
    return (0);
}
