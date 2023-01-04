#include <stdint.h>
#include <stdbool.h>

#include "mltaln.h"

#define DEBUG 0
#define IODEBUG 0
#define SCOREOUT 0
#define SHISHAGONYU 0  // for debug

#define REPORTCOSTS 0

typedef struct msacompactdistmtxthread_arg {
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

static TbfastOpts
arguments(int argc, char* argv[], int* pac, char** pav, int* tac, char** tav)  // 2 kai yobaremasu.
{
    TbfastOpts opts = {};

    int c;
    int i;

    nthread = 1;
    opts.nthreadtb = 1;
    nthreadpair = 1;
    outnumber = 0;
    scoreout = 0;
    spscoreout = 0;
    rnaprediction = 'm';
    rnakozo = 0;
    nevermemsave = 0;
    inputfile = NULL;
    addfile = NULL;
    addprofile = 1;
    fftkeika = 0;
    constraint = 0;
    nblosum = 62;
    fmodel = 0;
    calledByXced = 0;
    devide = 0;
    use_fft = 0;  // chuui
    force_fft = 0;
    fftscore = 1;
    fftRepeatStop = 0;
    fftNoAnchStop = 0;
    weight = 3;
    utree = 1;
    tbutree = 1;
    refine = 0;
    check = 1;
    cut = 0.0;
    disp = 0;
    outgap = 1;
    alg = 'A';
    mix = 0;
    tbitr = 0;
    scmtd = 5;
    tbweight = 0;
    tbrweight = 3;
    checkC = 0;
    treemethod = 'X';
    sueff_global = 0.1;
    contin = 0;
    scoremtx = 1;
    kobetsubunkatsu = 0;
    ppenalty_dist = NOTSPECIFIED;
    ppenalty = NOTSPECIFIED;
    penalty_shift_factor = 1000.0;
    ppenalty_ex = NOTSPECIFIED;
    poffset = NOTSPECIFIED;
    kimuraR = NOTSPECIFIED;
    pamN = NOTSPECIFIED;
    geta2 = GETA2;
    fftWinSize = NOTSPECIFIED;
    fftThreshold = NOTSPECIFIED;
    RNAppenalty = NOTSPECIFIED;
    RNAppenalty_ex = NOTSPECIFIED;
    RNApthr = NOTSPECIFIED;
    TMorJTT = JTT;
    consweight_multi = 1.0;
    consweight_rna = 0.0;
    legacygapcost = 0;
    specificityconsideration = 0.0;
    specifictarget = 0;
    nwildcard = 0;
    nadd = 0;

    if (pac) {
        pav[0] = "tbfast-pair";
        *pac = 1;
        tav[0] = "tbfast";
        *tac = 1;

        for (i = 0; i < argc; i++) {
            if (argv[i][0] == '_') {
                opts.callpairlocalalign = 1;

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
                    inputfile = *++argv;
                    --argc;
                    goto nextoption;
                case 'I':
                    nadd = myatoi(*++argv);
                    --argc;
                    goto nextoption;
                case 'e':
                    RNApthr = (int)(atof(*++argv) * 1000 - 0.5);
                    --argc;
                    goto nextoption;
                case 'o':
                    RNAppenalty = (int)(atof(*++argv) * 1000 - 0.5);
                    --argc;
                    goto nextoption;
                case 'V':
                    ppenalty_dist = (int)(atof(*++argv) * 1000 - 0.5);
                    --argc;
                    goto nextoption;
                case 'f':
                    ppenalty = (int)(atof(*++argv) * 1000 - 0.5);
                    --argc;
                    goto nextoption;
                case 'Q':
                    penalty_shift_factor = atof(*++argv);
                    --argc;
                    goto nextoption;
                case 'g':
                    ppenalty_ex = (int)(atof(*++argv) * 1000 - 0.5);
                    --argc;
                    goto nextoption;
                case 'h':
                    poffset = (int)(atof(*++argv) * 1000 - 0.5);
                    //					fprintf( stderr, "poffset = %d\n", poffset );
                    --argc;
                    goto nextoption;
                case 'k':
                    kimuraR = myatoi(*++argv);
                    //					fprintf( stderr, "kappa = %d\n", kimuraR );
                    --argc;
                    goto nextoption;
                case 'b':
                    nblosum = myatoi(*++argv);
                    scoremtx = 1;
                    //					fprintf( stderr, "blosum %d / kimura 200\n", nblosum );
                    --argc;
                    goto nextoption;
                case 'j':
                    pamN = myatoi(*++argv);
                    scoremtx = 0;
                    TMorJTT = JTT;
                    //					fprintf( stderr, "jtt/kimura %d\n", pamN );
                    --argc;
                    goto nextoption;
                case 'm':
                    pamN = myatoi(*++argv);
                    scoremtx = 0;
                    TMorJTT = TM;
                    //					fprintf( stderr, "tm %d\n", pamN );
                    --argc;
                    goto nextoption;
                case 'l':
                    fastathreshold = atof(*++argv);
                    constraint = 2;
                    --argc;
                    goto nextoption;
                case 'r':
                    consweight_rna = atof(*++argv);
                    rnakozo = 1;
                    --argc;
                    goto nextoption;
                case 'c':
                    consweight_multi = atof(*++argv);
                    --argc;
                    goto nextoption;
                case 'C':
                    nthreadpair = nthread = myatoi(*++argv);
                    --argc;
                    nthread = 0;
                    goto nextoption;
                case 's':
                    specificityconsideration = (double)myatof(*++argv);
                    //					fprintf( stderr, "specificityconsideration = %f\n", specificityconsideration );
                    --argc;
                    goto nextoption;
                case 'R':
                    rnaprediction = 'r';
#if 1
                case 'a':
                    fmodel = 1;
                    break;
#endif
                case 'K':
                    addprofile = 0;
                    break;
                case 'y':
                    opts.distout = 1;
                    break;
                case 't':
                    opts.treeout = 1;
                    break;
                case '^':
                    opts.treeout = 2;
                    break;
                case 'T':
                    opts.noalign = 1;
                    break;
                case 'D':
                    dorp = 'd';
                    break;
                case 'P':
                    dorp = 'p';
                    break;
                case 'L':
                    legacygapcost = 1;
                    break;
#if 1
                case 'O':
                    outgap = 0;
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
                    spscoreout = 1;  // 2014/Dec/30, sp score
                    break;
#endif
                case 'H':
                    opts.subalignment = 1;
                    opts.subalignmentoffset = myatoi(*++argv);
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
                    treemethod = 'X';
                    sueff_global = atof(*++argv);
                    //					fprintf( stderr, "sueff_global = %f\n", sueff_global );
                    --argc;
                    goto nextoption;
                case 'E':
                    treemethod = 'E';
                    break;
                case 'q':
                    treemethod = 'q';
                    break;
                case 'n':
                    outnumber = 1;
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
                    alg = 'd';
                    break;
                case 'A':
                    alg = 'A';
                    break;
                case 'M':
                    alg = 'M';
                    break;
                case 'N':
                    nevermemsave = 1;
                    break;
                case 'B':  // hitsuyou! memopt -M -B no tame
                    break;
                case 'F':
                    use_fft = 1;
                    break;
                case 'G':
                    force_fft = 1;
                    use_fft = 1;
                    break;
                case 'U':
                    opts.treein = 1;
                    break;
                case 'u':
                    tbrweight = 0;
                    weight = 0;
                    break;
                case 'v':
                    tbrweight = 3;
                    break;
                case 'd':
                    opts.multidist = 1;
                    break;
#if 0
				case 'd':
					disp = 1;
					break;
#endif
                    /* Modified 01/08/27, default: user tree */
                case 'J':
                    tbutree = 0;
                    break;
/* modification end. */
#if 0
				case 'z':
					fftThreshold = myatoi( *++argv );
					--argc; 
					goto nextoption;
#endif
                case 'w':
                    fftWinSize = myatoi(*++argv);
                    --argc;
                    goto nextoption;
                case 'W':
                    minimumweight = atof(*++argv);
                    //					fprintf( stderr, "minimumweight = %f\n", minimumweight );
                    --argc;
                    goto nextoption;
                case 'Y':
                    opts.keeplength = 1;
                    break;
                case 'z':
                    opts.mapout = 2;
                    break;
                case 'Z':
                    opts.mapout = 1;
                    break;
                case 'p':
                    opts.smoothing = 1;
                    break;
                case '=':
                    specifictarget = 1;
                    break;
                case ':':
                    nwildcard = 1;
                    break;
                case '+':
                    opts.outputhat23 = myatoi(*++argv);
                    reporterr("outputhat23=%d\n", opts.outputhat23);
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
        cut = atof((*argv));
        argc--;
    }
    if (argc != 0) {
        fprintf(stderr, "argc=%d, tbfast options: Check source file !\n", argc);
        exit(1);
    }
    if (tbitr == 1 && outgap == 0) {
        fprintf(stderr, "conflicting options : o, m or u\n");
        exit(1);
    }
    if (alg == 'C' && outgap == 0) {
        fprintf(stderr, "conflicting options : C, o\n");
        exit(1);
    }

    return opts;
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
            if (nthreadpair)
                fprintf(stderr, "\r% 5d / %d (thread %4d)", i, njob, thread_no);
            else
                fprintf(stderr, "\r% 5d / %d", i, njob);
        }

        for (j = i + 1; j < njob; j++) {
            tmpdist = distcompact_msa(seq[i], seq[j], skiptable[i], skiptable[j], selfscore[i], selfscore[j]);  // osoikedo,

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
treebase(TbfastOpts opts, int* nlen, char** aseq, int nadd, char* mergeoralign, char** mseq1, char** mseq2, int*** topol, Treedep* dep, double* effarr, int* alloclen, LocalHom** localhomtable, RNApair*** singlerna, double* effarr_kozo, int* targetmap, int* targetmapr, int ntarget, int* uselh, int nseed, int* nfilesfornode) {
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

    if (rnakozo && rnaprediction == 'm') {
        grouprna1 = (RNApair***)calloc(njob, sizeof(RNApair**));
        grouprna2 = (RNApair***)calloc(njob, sizeof(RNApair**));
    } else {
        grouprna1 = grouprna2 = NULL;
    }

    fftlog = AllocateIntVec(njob);
    effarr1 = AllocateDoubleVec(njob);
    effarr2 = AllocateDoubleVec(njob);
    indication1 = AllocateCharVec(150);
    indication2 = AllocateCharVec(150);
    gaplen = AllocateIntVec(*alloclen + 10);
    gapmap = AllocateIntVec(*alloclen + 10);
    alreadyaligned = AllocateIntVec(njob);
    dynamicmtx = AllocateDoubleMtx(nalphabets, nalphabets);
    localmem = calloc(sizeof(int*), 2);
    cpmxhist = (double***)calloc(njob - 1, sizeof(double**));
    for (i = 0; i < njob - 1; i++)
        cpmxhist[i] = NULL;
    memhist = (int**)calloc(njob - 1, sizeof(int*));
    for (i = 0; i < njob - 1; i++)
        memhist[i] = NULL;

    swaplist = NULL;
    if (constraint && compacttree != 3) {
        if (specifictarget)
            swaplist = calloc(njob, sizeof(char));
        localhomshrink = (LocalHom***)calloc(njob, sizeof(LocalHom**));
        for (i = 0; i < njob; i++) {
            localhomshrink[i] = (LocalHom**)calloc(njob, sizeof(LocalHom*));
        }
    } else if (constraint && nseed) {
        localhomshrink = (LocalHom***)calloc(nseed, sizeof(LocalHom**));
        for (i = 0; i < nseed; i++)
            localhomshrink[i] = (LocalHom**)calloc(nseed, sizeof(LocalHom*));

        seedinlh1 = calloc(njob, sizeof(int));
        seedinlh2 = calloc(njob, sizeof(int));
    } else if (constraint && nseed == 0) {
        seedinlh1 = NULL;
        seedinlh2 = NULL;
        localhomshrink = NULL;
    }

    effarr1_kozo = AllocateDoubleVec(njob);
    effarr2_kozo = AllocateDoubleVec(njob);
    for (i = 0; i < njob; i++)
        effarr1_kozo[i] = 0.0;
    for (i = 0; i < njob; i++)
        effarr2_kozo[i] = 0.0;

    for (i = 0; i < njob - nadd; i++)
        alreadyaligned[i] = 1;
    for (i = njob - nadd; i < njob; i++)
        alreadyaligned[i] = 0;

    for (l = 0; l < njob; l++)
        fftlog[l] = 1;

    if (constraint && compacttree != 3) {
        if (specifictarget)
            calcimportance_target(njob, ntarget, effarr, aseq, localhomtable, targetmap, targetmapr, *alloclen);
        else
            calcimportance_half(njob, effarr, aseq, localhomtable, *alloclen);
    } else if (constraint && nseed)  // ie, compacttree == 3 && constraint > 0
    {
        dontcalcimportance_half(nseed, aseq, localhomtable);  //CHUUI
    }

    tscore = 0.0;
    for (l = 0; l < njob - 1; l++) {
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

        if (l != njob - 2) {
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

        makedynamicmtx(dynamicmtx, n_dis_consweight_multi, dep[l].distfromtip);

        len1 = strlen(aseq[m1]);
        len2 = strlen(aseq[m2]);
        if (*alloclen < len1 + len2) {
            fprintf(stderr, "\nReallocating..");
            *alloclen = (len1 + len2) + 1000;
            ReallocateCharMtx(aseq, njob, *alloclen + 10);
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
            clus1 = fastconjuction_noname(localmem[0], aseq, mseq1, effarr1, effarr, indication1, minimumweight, &orieff1);  // orieff tsukau!
            clus2 = fastconjuction_noname(localmem[1], aseq, mseq2, effarr2, effarr, indication2, minimumweight, &orieff2);  // orieff tsukau!
        }

        if (mergeoralign[l] == '1' || mergeoralign[l] == '2') {
            newgapstr = "=";
        } else
            newgapstr = "-";

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
            fprintf(stderr, "\rSTEP % 5d /%d ", l + 1, njob - 1);

#if REPORTCOSTS
        if (l < 1000 || l % 100 == 0)
            reporterr("\nclus1=%d, clus2=%d\n", clus1, clus2);
#endif

        if (constraint && compacttree != 3) {
            if (specifictarget)
                fastshrinklocalhom_target(localmem[0], localmem[1], localhomtable, localhomshrink, swaplist, targetmap);
            else
                fastshrinklocalhom_half(localmem[0], localmem[1], localhomtable, localhomshrink);
        } else if (constraint && nseed) {
            fastshrinklocalhom_half_seed(localmem[0], localmem[1], nseed, seedinlh1, seedinlh2, localhomtable, localhomshrink);
            for (i = 0; i < njob; i++)
                reporterr("seedinlh1[%d]=%d\n", i, seedinlh1[i]);
            for (i = 0; i < njob; i++)
                reporterr("seedinlh2[%d]=%d\n", i, seedinlh2[i]);
        }

        if (rnakozo && rnaprediction == 'm') {
            makegrouprna(grouprna1, singlerna, localmem[0]);
            makegrouprna(grouprna2, singlerna, localmem[1]);
        }

        if (!nevermemsave && (constraint != 2 && alg != 'M' && (len1 > 30000 || len2 > 30000))) {
            fprintf(stderr, "\nlen1=%d, len2=%d, Switching to the memsave mode.\n", len1, len2);
            alg = 'M';
            if (commonIP)
                FreeIntMtx(commonIP);
            commonIP = NULL;
            commonAlloc1 = 0;
            commonAlloc2 = 0;
        }

        if (fftlog[m1] && fftlog[m2])
            ffttry = (nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 1000 && clus2 < 1000);
        else
            ffttry = 0;
        //		ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 5000 && clus2 < 5000 ); // v6.708
        //		fprintf( stderr, "f=%d, len1/fftlog[m1]=%f, clus1=%d, len2/fftlog[m2]=%f, clus2=%d\n", ffttry, (double)len1/fftlog[m1], clus1, (double)len2/fftlog[m2], clus2 );
        //		fprintf( stderr, "f=%d, clus1=%d, fftlog[m1]=%d, clus2=%d, fftlog[m2]=%d\n", ffttry, clus1, fftlog[m1], clus2, fftlog[m2] );
        if (constraint == 2) {
            if (alg == 'M') {
                fprintf(stderr, "\n\nMemory saving mode is not supported.\n\n");
                exit(1);
            }
            //			fprintf( stderr, "c" );
            if (alg == 'A') {
                imp_match_init_strict(clus1, clus2, strlen(mseq1[0]), strlen(mseq2[0]), mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, swaplist, localmem[0], localmem[1], uselh, seedinlh1, seedinlh2, (compacttree == 3) ? l : -1, nfiles);
                if (rnakozo)
                    imp_rna(clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2);
#if REPORTCOSTS
//				reporterr(       "\n\n %d - %d (%d x %d) : \n", topol[l][0][0], topol[l][1][0], clus1, clus2 );
#endif
                pscore = A__align(dynamicmtx, penalty, penalty_ex, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, constraint, &dumdb, NULL, NULL, NULL, NULL, outgap, outgap, localmem[0][0], 1, cpmxchild0, cpmxchild1, cpmxhist + l, orieff1, orieff2);
            }
            if (alg == 'd') {
                imp_match_init_strictD(clus1, clus2, strlen(mseq1[0]), strlen(mseq2[0]), mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, swaplist, localmem[0], localmem[1], uselh, seedinlh1, seedinlh2, (compacttree == 3) ? l : -1, nfiles);
                if (rnakozo)
                    imp_rnaD(clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2);
                pscore = D__align(dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, constraint, &dumdb, outgap, outgap);
            } else if (alg == 'Q') {
                fprintf(stderr, "Not supported\n");
                exit(1);
            }
        } else if (force_fft || (use_fft && ffttry)) {
            fprintf(stderr, " f\b\b");
            if (alg == 'M') {
                fprintf(stderr, "m");
                pscore = Falign_udpari_long(NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog + m1);
            } else
                pscore = Falign(NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog + m1);
        } else {
            fprintf(stderr, " d\b\b");
            fftlog[m1] = 0;
            switch (alg) {
                case ('a'):
                    pscore = Aalign(mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen);
                    break;
                case ('M'):
                    fprintf(stderr, "m");
                    pscore = MSalignmm(dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, NULL, NULL, NULL, outgap, outgap, cpmxchild0, cpmxchild1, cpmxhist + l, orieff1, orieff2);
                    break;
                case ('A'):
                    pscore = A__align(dynamicmtx, penalty, penalty_ex, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, 0, &dumdb, NULL, NULL, NULL, NULL, outgap, outgap, localmem[0][0], 1, cpmxchild0, cpmxchild1, cpmxhist + l, orieff1, orieff2);
                    break;
                case ('d'):
                    pscore = D__align(dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, 0, &dumdb, outgap, outgap);
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

        if (disp)
            display(aseq, njob);

        if (mergeoralign[l] == '1') {
            reporterr("Check source!!\n");
            exit(1);
        }
        if (mergeoralign[l] == '2') {
            gapmaplen = strlen(mseq1[0]) - len1nocommongap + len1;
            adjustgapmap(gapmaplen, gapmap, mseq1[0]);
            if (opts.smoothing) {
                restorecommongapssmoothly(njob, njob - (clus1 + clus2), aseq, localmem[0], localmem[1], gapmap, *alloclen, '-');
                findnewgaps(0, mseq1, gaplen);
                insertnewgaps_bothorders(njob, alreadyaligned, aseq, localmem[0], localmem[1], gaplen, gapmap, gapmaplen, *alloclen, alg, '-');
            } else {
                restorecommongaps(njob, njob - (clus1 + clus2), aseq, localmem[0], localmem[1], gapmap, *alloclen, '-');
                findnewgaps(0, mseq1, gaplen);
                insertnewgaps(njob, alreadyaligned, aseq, localmem[0], localmem[1], gaplen, gapmap, *alloclen, alg, '-');
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
        for (i = 0; i < njob - 1; i++) {
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

    if (rnakozo && rnaprediction == 'm') {
        if (grouprna1)
            free(grouprna1);  // nakami ha?
        if (grouprna2)
            free(grouprna2);  // nakami ha?
        grouprna1 = grouprna2 = NULL;
    }

    if (constraint) {
        if (localhomshrink)  // nen no tame
        {
            if (compacttree == 3)
                m = nseed;
            else
                m = njob;
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
    Falign(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL);
    Falign_udpari_long(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL);
    D__align(NULL, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, 0, 0);
    A__align(NULL, 0, 0, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, 0, 0, -1, -1, NULL, NULL, NULL, 0.0, 0.0);
    imp_match_init_strictD(0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0);
    imp_match_init_strict(0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0);
    FreeCommonIP();
}

static void
WriteOptions(FILE* fp) {
    if (dorp == 'd')
        fprintf(fp, "DNA\n");
    else {
        if (scoremtx == 0)
            fprintf(fp, "JTT %dPAM\n", pamN);
        else if (scoremtx == 1)
            fprintf(fp, "BLOSUM %d\n", nblosum);
        else if (scoremtx == 2)
            fprintf(fp, "M-Y\n");
    }
    fprintf(stderr, "Gap Penalty = %+5.2f, %+5.2f, %+5.2f\n", (double)ppenalty / 1000, (double)ppenalty_ex / 1000, (double)poffset / 1000);
    if (use_fft)
        fprintf(fp, "FFT on\n");

    fprintf(fp, "tree-base method\n");
    if (tbrweight == 0)
        fprintf(fp, "unweighted\n");
    else if (tbrweight == 3)
        fprintf(fp, "clustalw-like weighting\n");
    if (tbitr || tbweight) {
        fprintf(fp, "iterate at each step\n");
        if (tbitr && tbrweight == 0)
            fprintf(fp, "  unweighted\n");
        if (tbitr && tbrweight == 3)
            fprintf(fp, "  reversely weighted\n");
        if (tbweight)
            fprintf(fp, "  weighted\n");
        fprintf(fp, "\n");
    }

    fprintf(fp, "Gap Penalty = %+5.2f, %+5.2f, %+5.2f\n", (double)ppenalty / 1000, (double)ppenalty_ex / 1000, (double)poffset / 1000);

    if (alg == 'a')
        fprintf(fp, "Algorithm A\n");
    else if (alg == 'A')
        fprintf(fp, "Algorithm A+\n");
    else if (alg == 'C')
        fprintf(fp, "Apgorithm A+/C\n");
    else
        fprintf(fp, "Unknown algorithm\n");

    if (treemethod == 'X')
        fprintf(fp, "Tree = UPGMA (mix).\n");
    else if (treemethod == 'E')
        fprintf(fp, "Tree = UPGMA (average).\n");
    else if (treemethod == 'q')
        fprintf(fp, "Tree = Minimum linkage.\n");
    else
        fprintf(fp, "Unknown tree.\n");

    if (use_fft) {
        fprintf(fp, "FFT on\n");
        if (dorp == 'd')
            fprintf(fp, "Basis : 4 nucleotides\n");
        else {
            if (fftscore)
                fprintf(fp, "Basis : Polarity and Volume\n");
            else
                fprintf(fp, "Basis : 20 amino acids\n");
        }
        fprintf(fp, "Threshold   of anchors = %d%%\n", fftThreshold);
        fprintf(fp, "window size of anchors = %dsites\n", fftWinSize);
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
    TbfastOpts opts = arguments(argc, argv, &pac, pav, &tac, tav);

    if (opts.treein) {
        int    dumx, dumy;
        double dumz;
        opts.treein = check_guidetreefile(&dumx, &dumy, &dumz);
        if (opts.treein == 'C') {
            compacttree = 2;
            opts.treein = 0;
            use_fft = 0;
        } else if (opts.treein == 'n') {
            compacttree = 3;
            opts.treein = 0;
            use_fft = 0;
        }
    }

    reporterr("treein = %d\n", opts.treein);
    reporterr("compacttree = %d\n", compacttree);

    if (fastathreshold < 0.0001)
        constraint = 0;

    if (inputfile) {
        infp = fopen(inputfile, "rb");
        if (!infp) {
            fprintf(stderr, "Cannot open %s\n", inputfile);
            exit(1);
        }
    } else
        infp = stdin;

    getnumlen(infp);
    rewind(infp);

    nkozo = 0;

#if !defined(mingw) && !defined(_MSC_VER)
    setstacksize(200 * njob);
#endif

    if (opts.subalignment) {
        readsubalignmentstable(njob, NULL, NULL, &nsubalignments, &maxmem);
        fprintf(stderr, "nsubalignments = %d\n", nsubalignments);
        fprintf(stderr, "maxmem = %d\n", maxmem);
        subtable = AllocateIntMtx(nsubalignments, maxmem + 1);
        insubtable = AllocateIntVec(njob);
        for (i = 0; i < njob; i++)
            insubtable[i] = 0;
        preservegaps = AllocateIntVec(njob);
        for (i = 0; i < njob; i++)
            preservegaps[i] = 0;
        subalnpt = AllocateCharCub(nsubalignments, maxmem, 0);
        readsubalignmentstable(njob, subtable, preservegaps, NULL, NULL);
    }

    seq = AllocateCharMtx(njob, nlenmax + 1);
    mseq1 = AllocateCharMtx(njob, 0);
    mseq2 = AllocateCharMtx(njob, 0);

    name = AllocateCharMtx(njob, B + 1);
    nlen = AllocateIntVec(njob);
    selfscore = AllocateIntVec(njob);

    topol = AllocateIntCub(njob, 2, 0);
    len = AllocateFloatMtx(njob, 2);
    eff = AllocateDoubleVec(njob);
    kozoarivec = AllocateCharVec(njob);

    mergeoralign = AllocateCharVec(njob);

    dep = (Treedep*)calloc(njob, sizeof(Treedep));
    if (nadd)
        addmem = AllocateIntVec(nadd + 1);
    localmem = AllocateIntMtx(2, njob + 1);

    if (compacttree == 3)
        nfilesfornode = calloc(sizeof(int), njob - 1);

#if REPORTCOSTS
    reporterr("before allocating iscore\n");
    use_getrusage();
#endif

    if (tbutree && compacttree != 3)
        iscore = AllocateFloatHalfMtx(njob);

    opts.ndeleted = 0;

    readData_pointer(infp, name, nlen, seq);
    fclose(infp);

    if (opts.treein) {
        loadtree(njob, topol, len, name, dep, opts.treeout);
        fprintf(stderr, "\ndone.\n\n");
        if (opts.callpairlocalalign && specificityconsideration > 0.0) {
            int* mem0 = calloc(sizeof(int), njob);
            int* mem1 = calloc(sizeof(int), njob);
            expdist = AllocateDoubleMtx(njob, njob);
            for (i = 0; i < njob - 1; i++) {
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
        targetmap = calloc(njob, sizeof(int));
        ntarget = 0;
        for (i = 0; i < njob; i++) {
            targetmap[i] = -1;
            if (!strncmp(name[i] + 1, "_focus_", 7))
                targetmap[i] = ntarget++;
        }
        targetmapr = calloc(ntarget, sizeof(int));
        for (i = 0; i < njob; i++)
            if (targetmap[i] != -1)
                targetmapr[targetmap[i]] = i;

    } else if (compacttree != 3) {
        ntarget = njob;
        targetmap = calloc(njob, sizeof(int));
        targetmapr = calloc(njob, sizeof(int));
        for (i = 0; i < njob; i++)
            targetmap[i] = targetmapr[i] = i;
    }

    if (constraint && compacttree != 3) {
        ilim = njob;
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

        if (opts.callpairlocalalign) {
            pairlocalalign(ctx, njob, name, seq, iscore, localhomtable, pac, pav, expdist);
            arguments(tac, tav, NULL, NULL, NULL, NULL);
            opts.callpairlocalalign = 1;
            if (expdist)
                FreeDoubleMtx(expdist);
            expdist = NULL;
            if (fastathreshold < 0.0001)
                constraint = 0;
            if (compacttree != 3) {
                for (ilim = njob, i = 0; i < ntarget; i++) {
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
                        readlocalhomtable2_half(prep, njob, localhomtable, kozoarivec);
                    fclose(prep);
                    fprintf(stderr, "\ndone.\n");
                } else
                    fprintf(stderr, "No hat3.seed. No problem.\n");

                if (opts.outputhat23) {
                    prep = fopen("hat3", "w");
                    if (!prep)
                        ErrorExit("Cannot open hat3 to write.");

                    fprintf(stderr, "Writing hat3 for iterative refinement\n");
                    if (specifictarget)
                        ilim = ntarget;
                    else
                        ilim = njob - 1;
                    for (i = 0; i < ilim; i++) {
                        if (specifictarget) {
                            jst = 0;
                            jj = 0;
                        } else {
                            jst = i;
                            jj = 0;
                        }
                        for (j = jst; j < njob; j++, jj++) {
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
                    WriteFloatHat2_pointer_halfmtx(prep, njob, name, iscore);
                    fclose(prep);
                } else if (opts.distout) {
                    prep = fopen("hat2", "w");
                    WriteFloatHat2_pointer_halfmtx(prep, njob, name, iscore);
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
                readlocalhomtable2_half(prep, njob, localhomtable, kozoarivec);
            fclose(prep);
            fprintf(stderr, "\ndone.\n");
        }

        nkozo = 0;
        for (i = 0; i < njob; i++) {
            //			fprintf( stderr, "kozoarivec[%d] = %d\n", i, kozoarivec[i] );
            if (kozoarivec[i])
                nkozo++;
        }
        if (nkozo) {
            topol_kozo = AllocateIntCub(nkozo, 2, 0);
            len_kozo = AllocateFloatMtx(nkozo, 2);
            iscore_kozo = AllocateFloatHalfMtx(nkozo);
            eff_kozo = AllocateDoubleVec(nkozo);
            eff_kozo_mapped = AllocateDoubleVec(njob);
        }
    } else if (compacttree != 3) {
        if (opts.callpairlocalalign) {
            pairlocalalign(ctx, njob, name, seq, iscore, NULL, pac, pav, expdist);
            arguments(tac, tav, NULL, NULL, NULL, NULL);
            opts.callpairlocalalign = 1;
            if (expdist)
                FreeDoubleMtx(expdist);
            expdist = NULL;
            if (fastathreshold < 0.0001)
                constraint = 0;
            fprintf(stderr, "blosum %d / kimura 200\n", nblosum);
            fprintf(stderr, "scoremtx=%d\n", scoremtx);
            fprintf(stderr, "fastathreshold=%f\n", fastathreshold);
        }
        if (opts.distout || opts.outputhat23) {
            reporterr("\nwriting hat2 (1)\n");
            prep = fopen("hat2", "w");
            WriteFloatHat2_pointer_halfmtx(prep, njob, name, iscore);
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
            for (i = 0; i < njob; i++)
                if (strncmp(name[i] + 1, "_seed", 4))
                    break;  // konran!!!
            nseed = i;

            topol_kozo = AllocateIntCub(nseed, 2, 0);
            len_kozo = AllocateFloatMtx(nseed, 2);
            iscore_kozo = AllocateFloatHalfMtx(nseed);
            eff_kozo = AllocateDoubleVec(nseed);
            eff_kozo_mapped = AllocateDoubleVec(njob);

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
            for (i = 0; i < njob; i++) {
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
    constants(ctx, njob, seq);

    if (sueff_global < 0.0001 || compacttree == 3) {
        nthreadreadlh = nthread;
        if (nthreadreadlh == 0)
            nthreadreadlh = 1;

        opts.nthreadtb = 0;
        nthread = 0;
    } else {
        nthreadreadlh = 1;
        opts.nthreadtb = nthread;
    }

    initSignalSM();
    initFiles();
    WriteOptions(trap_g);

    if (opts.distout && !opts.treeout && opts.noalign) {
        writeData_pointer(prep_g, njob, name, seq);
        fprintf(stderr, "\n");
        goto chudan;
    }

    c = seqcheck(seq);
    if (c) {
        fprintf(stderr, "Illegal character %c\n", c);
        exit(1);
    }

    if (nadd && opts.keeplength) {
        originalgaps = (char*)calloc(nlenmax + 1, sizeof(char));
        recordoriginalgaps(originalgaps, njob - nadd, seq);

        if (opts.mapout) {
            addbk = (char**)calloc(nadd + 1, sizeof(char*));
            for (i = 0; i < nadd; i++) {
                ien = strlen(seq[njob - nadd + i]);
                addbk[i] = (char*)calloc(ien + 1, sizeof(char));
                gappick0(addbk[i], seq[njob - nadd + i]);
            }
            addbk[nadd] = NULL;
        } else
            addbk = NULL;
    } else {
        originalgaps = NULL;
        addbk = NULL;
    }

    if (!opts.treein) {
        reporterr("tbutree = %d, compacttree = %d\n", tbutree, compacttree);
        if (compacttree == 3) {
            iscore = NULL;
        } else if (tbutree == 0 && compacttree) {
            iscore = NULL;
            reporterr("Making a compact tree from msa, step 1.. \n");
            skiptable = AllocateIntMtx(njob, 0);
            makeskiptable(njob, skiptable, seq);
            mindistfrom = (int*)calloc(njob, sizeof(int));
            mindist = (double*)calloc(njob, sizeof(double));
            partmtx = preparepartmtx(njob);

            for (i = 0; i < njob; i++) {
                selfscore[i] = (int)naivepairscorefast(seq[i], seq[i], skiptable[i], skiptable[i], penalty_dist);
            }

            {
                int jobpos = 0;

                {
                    for (j = 0; j < njob; j++) {
                        mindist[j] = 999.9;
                        mindistfrom[j] = -1;
                    }

                    msacompactdistmtxthread_arg_t targ = {
                        .thread_no = 0,
                        .njob = njob,
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

                for (i = 0; i < njob; i++)
                    mindist[i] -= preferenceval(i, mindistfrom[i], njob);
            }
            reporterr("\rdone.                                          \n");
        } else if (tbutree == 0 && compacttree == 0) {
            reporterr("Making a distance matrix from msa .. \n");
            iscore = AllocateFloatHalfMtx(njob);

            for (i = 1; i < njob; i++) {
                if (nlen[i] != nlen[0]) {
                    fprintf(stderr, "Input pre-aligned seqences or make hat2.\n");
                    exit(1);
                }
            }

            skiptable = AllocateIntMtx(njob, 0);
            makeskiptable(njob, skiptable, seq);
            ien = njob - 1;
            for (i = 0; i < njob; i++) {
                selfscore[i] = (int)naivepairscorefast(seq[i], seq[i], skiptable[i], skiptable[i], penalty_dist);
            }

            {
                for (i = 0; i < ien; i++) {
                    if (i % 10 == 0) {
                        fprintf(stderr, "\r% 5d / %d", i, ien);
                    }
                    ssi = selfscore[i];
                    for (j = i + 1; j < njob; j++) {
                        ssj = selfscore[j];
                        bunbo = MIN(ssi, ssj);
                        if (bunbo == 0.0)
                            iscore[i][j - i] = 2.0;
                        else
                            iscore[i][j - i] = (1.0 - naivepairscorefast(seq[i], seq[j], skiptable[i], skiptable[j], penalty_dist) / bunbo) * 2.0;  // 2014/Aug/15 fast
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
            if (opts.callpairlocalalign) {
                if (opts.multidist) {
                    reporterr("Bug in v7.290.  Please email katoh@ifrec.osaka-u.ac.jp\n");
                    exit(1);
                }
            } else {
                if (opts.multidist) {
                    fprintf(stderr, "Loading 'hat2n' (aligned sequences - new sequences) ... ");
                    prep = fopen("hat2n", "r");
                    if (prep == NULL)
                        ErrorExit("Make hat2.");
                    readhat2_doublehalf_pointer(prep, njob, iscore);
                    fclose(prep);
                    fprintf(stderr, "done.\n");

                    fprintf(stderr, "Loading 'hat2i' (aligned sequences) ... ");
                    prep = fopen("hat2i", "r");
                    if (prep == NULL)
                        ErrorExit("Make hat2i.");
                    readhat2_doublehalf_pointer(prep, njob - nadd, iscore);
                    fclose(prep);
                    fprintf(stderr, "done.\n");
                } else {
                    fprintf(stderr, "Loading 'hat2' ... ");
                    prep = fopen("hat2", "r");
                    if (prep == NULL)
                        ErrorExit("Make hat2.");
                    readhat2_doublehalf_pointer(prep, njob, iscore);
                    fclose(prep);
                    fprintf(stderr, "done.\n");
                }

                if (opts.distout) {
                    reporterr("\nwriting hat2 (2)\n");
                    hat2p = fopen("hat2", "w");
                    WriteFloatHat2_pointer_halfmtx(hat2p, njob, name, iscore);
                    fclose(hat2p);
                }
            }
        }

        if (nkozo && compacttree != 3) {
            ien = njob - 1;
            ik = 0;
            for (i = 0; i < ien; i++) {
                jk = ik + 1;
                for (j = i + 1; j < njob; j++) {
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
            for (i = 0; i < njob; i++) {
                if (kozoarivec[i])
                    selfscore[i] = naivepairscore11(seq[i], seq[i], 0.0);
                else
                    selfscore[i] = -1;
            }
            ien = njob - 1;
            ik = 0;
            for (i = 0; i < ien; i++) {
                jk = ik + 1;
                for (j = i + 1; j < njob; j++) {
                    if (kozoarivec[i] && kozoarivec[j]) {
                        reporterr("seq0=%s\n", seq[i]);
                        reporterr("seq1=%s\n", seq[j]);
                        reporterr("selfscore0=%d\n", selfscore[0]);
                        reporterr("selfscore1=%d\n", selfscore[1]);
                        iscore_kozo[ik][jk - ik] = distdp_noalign(seq[i], seq[j], (double)selfscore[i], (double)selfscore[j], alloclen);
                        reporterr("iscore_kozo[%d][%d]=%f\n", ik, jk, iscore_kozo[ik][jk - ik]);
                    }
                    if (kozoarivec[j])
                        jk++;
                }
                if (kozoarivec[i])
                    ik++;
                G__align11_noalign(NULL, 0, 0, NULL, NULL);
            }
        }

        if (opts.subalignment) {
            fprintf(stderr, "Constructing a UPGMA tree ... ");
            fixed_supg_double_realloc_nobk_halfmtx_treeout_constrained(njob, iscore, topol, len, name, dep, nsubalignments, subtable, 1);
        } else if (compacttree == 3) {
            fp = fopen("hat3dir/tree", "rb");  // windows no tame rb
            treein_bin(fp, njob, topol, len, dep, nfilesfornode);
            fclose(fp);

            if (constraint) {
                uselh = AllocateIntVec(njob);
                fp = fopen("hat3dir/uselh", "rb");
                if (uselhin(fp, njob, uselh)) {
                    free(uselh);
                    uselh = NULL;
                }
                fclose(fp);
                //				for( i=0; i<njob; i++ ) reporterr( "uselh[%d]=%d\n", i, uselh[i] );
            }
        } else if (tbutree == 0 && compacttree)  // tbutree != 0 no toki (aln->mtx) ha, 6merdistance -> disttbfast.c; dp distance -> muzukashii
        {
            reporterr("Constructing a tree ... nthread=%d", nthread);
            compacttree_memsaveselectable(njob, partmtx, mindistfrom, mindist, NULL, selfscore, seq, skiptable, topol, len, name, NULL, dep, opts.treeout, compacttree, 1);

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
        } else if (opts.treeout) {
            fprintf(stderr, "Constructing a UPGMA tree ... ");
            fixed_musclesupg_double_realloc_nobk_halfmtx_treeout_memsave(njob, iscore, topol, len, name, dep, 1, opts.treeout);
        } else {
            fprintf(stderr, "Constructing a UPGMA tree ... ");
            fixed_musclesupg_double_realloc_nobk_halfmtx_memsave(njob, iscore, topol, len, dep, 1, 1);
        }

        if (nkozo) {
            fixed_musclesupg_double_realloc_nobk_halfmtx(nkozo, iscore_kozo, topol_kozo, len_kozo, NULL, 1, 1);
        }
        fprintf(stderr, "\ndone.\n\n");
    }

    localmem[0][0] = -1;
    topolorderz(localmem[0], topol, dep, njob - 2, 2);

    orderfp = fopen("order", "w");
    if (!orderfp) {
        fprintf(stderr, "Cannot open 'order'\n");
        exit(1);
    }

    for (i = 0; i < njob; i++)
        fprintf(orderfp, "%d\n", localmem[0][i]);

    fclose(orderfp);

    if (opts.treeout && opts.noalign) {
        writeData_pointer(prep_g, njob, name, seq);
        fprintf(stderr, "\n");
        goto chudan;  // 2016Jul31
    }

    if (tbrweight) {
        weight = 3;
        counteff_simple_double_nostatic_memsave(njob, topol, len, dep, eff);
        for (i = njob - nadd; i < njob; i++)
            eff[i] /= (double)100;
        if (nkozo) {
            for (i = 0, j = 0; i < njob; i++) {
                if (kozoarivec[i]) {
                    eff_kozo_mapped[i] = eff[i];  // single weight
                    j++;
                } else
                    eff_kozo_mapped[i] = 0.0;
            }
        }
    } else {
        for (i = 0; i < njob; i++)
            eff[i] = 1.0;
        if (nkozo) {
            for (i = 0; i < njob; i++) {
                if (kozoarivec[i])
                    eff_kozo_mapped[i] = 1.0;
                else
                    eff_kozo_mapped[i] = 0.0;
            }
        }
    }

    if (iscore)
        FreeFloatHalfMtx(iscore, njob);
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

    alloclen = nlenmax * 2 + 1;  //chuui!
    bseq = AllocateCharMtx(njob, alloclen);

    if (nadd) {
        alignmentlength = strlen(seq[0]);
        for (i = 0; i < njob - nadd; i++) {
            if (alignmentlength != strlen(seq[i])) {
                fprintf(stderr, "#################################################################################\n");
                fprintf(stderr, "# ERROR!                                                                        #\n");
                fprintf(stderr, "# The original%4d sequences must be aligned                                    #\n", njob - nadd);
                fprintf(stderr, "#################################################################################\n");
                exit(1);
            }
        }
        if (addprofile) {
            alignmentlength = strlen(seq[njob - nadd]);
            for (i = njob - nadd; i < njob; i++) {
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
                addmem[i] = njob - nadd + i;
            addmem[nadd] = -1;
            foundthebranch = 0;
            for (i = 0; i < njob - 1; i++) {
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
            commongappick(nadd, seq + njob - nadd);
            for (i = njob - nadd; i < njob; i++)
                strcpy(bseq[i], seq[i]);
        } else {
            for (i = 0; i < njob - 1; i++)
                mergeoralign[i] = 'n';
            for (i = 0; i < nadd; i++)
                addmem[i] = njob - nadd + i;
            addmem[nadd] = -1;
            for (i = 0; i < njob - 1; i++) {
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

            for (i = njob - nadd; i < njob; i++)
                gappick0(bseq[i], seq[i]);
        }

        commongappick(njob - nadd, seq);
        for (i = 0; i < njob - nadd; i++)
            strcpy(bseq[i], seq[i]);
    } else if (opts.subalignment) {
        for (i = 0; i < njob - 1; i++)
            mergeoralign[i] = 'a';
        for (i = 0; i < nsubalignments; i++) {
            fprintf(stderr, "Checking subalignment %d:\n", i + 1);
            alignmentlength = strlen(seq[subtable[i][0]]);
            //			for( j=0; subtable[i][j]!=-1; j++ )
            //				fprintf( stderr, " %d. %-30.30s\n", subtable[i][j]+1, name[subtable[i][j]]+1 );
            for (j = 0; subtable[i][j] != -1; j++) {
                if (subtable[i][j] >= njob) {
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
                    if (opts.subalignmentoffset) {
                        fprintf(stderr, "#\n");
                        fprintf(stderr, "# You specified seed alignment(s) consisting of %d sequences.\n", opts.subalignmentoffset);
                        fprintf(stderr, "# In this case, the rule of numbering is:\n");
                        fprintf(stderr, "#   The aligned seed sequences are numbered as 1 .. %d\n", opts.subalignmentoffset);
                        fprintf(stderr, "#   The input sequences to be aligned are numbered as %d .. %d\n", opts.subalignmentoffset + 1, opts.subalignmentoffset + njob);
                    }
                    fprintf(stderr, "###############################################################################\n");
                    fprintf(stderr, "\n");
                    exit(1);
                }
                insubtable[subtable[i][j]] = 1;
            }
            for (j = 0; j < njob - 1; j++) {
                if (includemember(topol[j][0], subtable[i]) && includemember(topol[j][1], subtable[i])) {
                    mergeoralign[j] = 'n';
                }
            }
            foundthebranch = 0;
            for (j = 0; j < njob - 1; j++) {
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
                if (opts.subalignmentoffset) {
                    fprintf(stderr, "#\n");
                    fprintf(stderr, "# You specified seed alignment(s) consisting of %d sequences.\n", opts.subalignmentoffset);
                    fprintf(stderr, "# In this case, the rule of numbering is:\n");
                    fprintf(stderr, "#   The aligned seed sequences are numbered as 1 .. %d\n", opts.subalignmentoffset);
                    fprintf(stderr, "#   The input sequences to be aligned are numbered as %d .. %d\n", opts.subalignmentoffset + 1, opts.subalignmentoffset + njob);
                }
                fprintf(stderr, "############################################################################### \n");
                fprintf(stderr, "\n");
                exit(1);
            }
        }

        for (i = 0; i < njob; i++) {
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
        for (i = 0; i < njob; i++)
            gappick0(bseq[i], seq[i]);
        for (i = 0; i < njob - 1; i++)
            mergeoralign[i] = 'a';
    }

    if (rnakozo && rnaprediction == 'm') {
        singlerna = (RNApair***)calloc(njob, sizeof(RNApair**));
        prep = fopen("hat4", "r");
        if (prep == NULL)
            ErrorExit("Make hat4 using mccaskill.");
        fprintf(stderr, "Loading 'hat4' ... ");
        for (i = 0; i < njob; i++) {
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

    treebase(opts, nlen, bseq, nadd, mergeoralign, mseq1, mseq2, topol, dep, eff, &alloclen, localhomtable, singlerna, eff_kozo_mapped, targetmap, targetmapr, ntarget, uselh, nseed, nfilesfornode);

    fprintf(stderr, "\ndone.\n");

    if (opts.keeplength) {
        dlf = fopen("_deletelist", "w");
        deletelist = (GapPos**)calloc(nadd + 1, sizeof(GapPos*));
        for (i = 0; i < nadd; i++) {
            deletelist[i] = calloc(1, sizeof(GapPos));
            deletelist[i][0].pos = -1;
            deletelist[i][0].len = 0;
        }
        deletelist[nadd] = NULL;
        opts.ndeleted = deletenewinsertions_whole(njob - nadd, nadd, bseq, bseq + njob - nadd, deletelist);

        for (i = 0; i < nadd; i++) {
            if (deletelist[i])
                for (j = 0; deletelist[i][j].pos != -1; j++)
                    fprintf(dlf, "%d %d %d\n", njob - nadd + i, deletelist[i][j].pos, deletelist[i][j].len);  // 0origin
        }
        fclose(dlf);

        restoreoriginalgaps(njob, bseq, originalgaps);
        free(originalgaps);
        originalgaps = NULL;

        if (opts.mapout) {
            dlf = fopen("_deletemap", "w");
            if (opts.mapout == 1)
                reconstructdeletemap(nadd, addbk, deletelist, bseq + njob - nadd, dlf, name + njob - nadd);
            else
                reconstructdeletemap_compact(nadd, addbk, deletelist, seq + njob - nadd, dlf, name + njob - nadd);
            FreeCharMtx(addbk);
            addbk = NULL;
            fclose(dlf);
        }

        for (i = 0; deletelist[i] != NULL; i++)
            free(deletelist[i]);
        free(deletelist);
        deletelist = NULL;
    }

    if (scoreout) {
        unweightedspscore = plainscore(njob, bseq);
        fprintf(stderr, "\nSCORE %s = %.0f, ", "(treebase)", unweightedspscore);
        fprintf(stderr, "SCORE / residue = %f", unweightedspscore / (njob * strlen(bseq[0])));
        fprintf(stderr, "\n\n");
    }

    fprintf(trap_g, "done.\n");
    free(mergeoralign);
    freeconstants();

    if (rnakozo && rnaprediction == 'm') {
        if (singlerna)  // nen no tame
        {
            for (i = 0; i < njob; i++) {
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

    writeData_pointer(prep_g, njob, name, bseq);

#if IODEBUG
    fprintf(stderr, "OSHIMAI\n");
#endif

    if (constraint && compacttree != 3) {
        if (specifictarget)
            FreeLocalHomTable_part(localhomtable, ntarget, njob);
        else
            FreeLocalHomTable_half(localhomtable, njob);
    } else if (constraint && nkozo) {
        FreeLocalHomTable_half(localhomtable, nkozo);
    }

    if (compacttree != 3) {
        free(targetmap);
        free(targetmapr);
    }

    if (constraint && compacttree == 3 && uselh)
        free(uselh);
    uselh = NULL;

    if (spscoreout)
        reporterr("Unweighted sum-of-pairs score = %10.5f\n", sumofpairsscore(njob, bseq));
    nthread = MAX(opts.nthreadtb, nthreadreadlh);  // toriaezu
    if (opts.ndeleted > 0) {
        reporterr("\nTo keep the alignment length, %d letters were DELETED.\n", opts.ndeleted);
        if (opts.mapout)
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
    closeFiles();
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
            FreeLocalHomTable_part(localhomtable, ntarget, njob);
        else
            FreeLocalHomTable_half(localhomtable, njob);
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
        FreeFloatHalfMtx(iscore, njob);
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

    freeconstants();
    closeFiles();
    FreeCommonIP();
    return (0);
}
