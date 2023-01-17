#include "mltaln.h"

#define DEBUG 0
#define CANONICALTREEFORMAT 1
#define MEMSAVE 1

#define HAT3SORTED 0
#define DISPPAIRID 0  // tbfast ha ugokanakunaru

#define LHBLOCKFACTOR 2
#define MINBLOCKLEN2 1000000000  // 100000 pairs * 100 sites * 100 sites

#define N0LOOPFIRST 0
#define YOUNGER0TREE 1  // --add ni hitsuyou

#define RECURSIVETOP 0

#define TREE7325 0

int
seqlen(Context* ctx, char* seq) {
    int val = 0;
    if (*ctx->newgapstr == '-') {
        while (*seq)
            if (*seq++ != '-')
                val++;
    } else {
        while (*seq) {
            if (*seq != '-' && *seq != *ctx->newgapstr)
                val++;
            seq++;
        }
    }
    return (val);
}

int
intlen(int* num) {
    int* numbk = num;
    while (*num++ != -1)
        ;
    return (num - numbk - 1);
}

void
intcat(int* s1, int* s2) {
    while (*s1 != -1)
        s1++;
    while (*s2 != -1) {
        //		reporterr(       "copying %d\n", *s2 );
        *s1++ = *s2++;
    }
    *s1 = -1;
}

void
intcpy(int* s1, int* s2) {
    while (*s2 != -1) {
        //		reporterr(       "copying %d\n", *s2 );
        *s1++ = *s2++;
    }
    *s1 = -1;
}

void
intncpy(int* s1, int* s2, int n) {
    while (n--)
        *s1++ = *s2++;
}

void
fltncpy(double* s1, double* s2, int n) {
    while (n--)
        *s1++ = *s2++;
}

#define END_OF_VEC -1

#define BLOCKSIZE 100
#define LARGEBLOCKSIZE 100

static void
setnearest(Bchain* acpt, double** eff, double* mindisfrompt, int* nearestpt, int pos) {
    int    j;
    double tmpdouble;
    double mindisfrom;
    int    nearest;
    //	double **effptpt;
    Bchain* acptj;

    mindisfrom = 999.9;
    nearest = -1;

    //	printf( "[%d], %f, dist=%d ->", pos, *mindisfrompt, *nearestpt );

    //	if( (acpt+pos)->next ) effpt = eff[pos]+(acpt+pos)->next->pos-pos;

    //	for( j=pos+1; j<nseq; j++ )
    for (acptj = (acpt + pos)->next; acptj != NULL; acptj = acptj->next) {
        j = acptj->pos;
        //		if( (tmpdouble=*effpt++) < *mindisfrompt )
        if ((tmpdouble = eff[pos][j - pos]) < mindisfrom) {
            mindisfrom = tmpdouble;
            nearest = j;
        }
    }
    //	effptpt = eff;
    //	for( j=0; j<pos; j++ )
    for (acptj = acpt; (acptj && acptj->pos != pos); acptj = acptj->next) {
        j = acptj->pos;
        //		if( (tmpdouble=(*effptpt++)[pos-j]) < *mindisfrompt )
        if ((tmpdouble = eff[j][pos - j]) < mindisfrom) {
            mindisfrom = tmpdouble;
            nearest = j;
        }
    }

    *mindisfrompt = mindisfrom;
    *nearestpt = nearest;
    //	printf( "%f, %d \n", pos, *mindisfrompt, *nearestpt );
}

static void
shufflelennum(Lennum* ary, int size) {
    int i;
    for (i = 0; i < size; i++) {
        int    j = rand() % size;
        Lennum t = ary[i];  // kouzoutai dainyu?
        ary[i] = ary[j];
        ary[j] = t;
    }
}

void
stringshuffle(int* ary, int size) {
    int i;
    for (i = 0; i < size; i++) {
        int j = rand() % size;
        int t = ary[i];
        ary[i] = ary[j];
        ary[j] = t;
    }
}

static int
compfunc(const void* p, const void* q) {
    return (((Lennum*)q)->len - ((Lennum*)p)->len);
}

void
limitlh(int* uselh, Lennum* in, int size, int limit) {
    int i;
    //	for(i=0;i<size;i++) reporterr( "i=%d, onum=%d, len=%d\n", i, in[i].num, in[i].len );

    qsort(in, size, sizeof(Lennum), compfunc);
    shufflelennum(in, size / 2);
    shufflelennum(in + size / 2, size - size / 2);

    if (limit > size)
        limit = size;
    //	reporterr( "numpairs=%llu, ULLONG_MAX=%llu, nn=%lld, INT_MAX=%d, n=%d\n", numpairs, ULLONG_MAX, nn, INT_MAX, n );

    for (i = 0; i < limit; i++)
        uselh[in[i].num] = 1;
    for (i = limit; i < size; i++)
        uselh[in[i].num] = 0;
}

void
sortbylength(int* uselh, Lennum* in, int size, unsigned long long numpairs) {
    int                i;
    unsigned long long nn;
    int                n;
    //	for(i=0;i<size;i++) reporterr( "i=%d, onum=%d, len=%d\n", i, in[i].num, in[i].len );

    qsort(in, size, sizeof(Lennum), compfunc);

    nn = ((unsigned long long)sqrt(1 + 8 * numpairs) + 1) / 2;
    if (nn > INT_MAX)
        nn = INT_MAX;

    n = (int)nn;
    if (n > size)
        n = size;
    //	reporterr( "numpairs=%llu, ULLONG_MAX=%llu, nn=%lld, INT_MAX=%d, n=%d\n", numpairs, ULLONG_MAX, nn, INT_MAX, n );

    for (i = 0; i < n; i++)
        uselh[in[i].num] = 1;
    for (i = n; i < size; i++)
        uselh[in[i].num] = 0;
}

typedef struct _TopDep {
    Treedep* dep;
    int***   topol;
} TopDep;

static TopDep* tdpglobal = NULL;
static int*
topolorder_lessargs(int* order, int pos) {
    if ((tdpglobal->dep)[pos].child0 == -1) {
        *order++ = (tdpglobal->topol)[pos][0][0];
        *order = -1;
    } else {
        order = topolorder_lessargs(order, (tdpglobal->dep)[pos].child0);
    }

    if ((tdpglobal->dep)[pos].child1 == -1) {
        *order++ = (tdpglobal->topol)[pos][1][0];
        *order = -1;
    } else {
        order = topolorder_lessargs(order, (tdpglobal->dep)[pos].child1);
    }

    return (order);
}

int*
topolorderz(int* order, int*** topol, Treedep* dep, int pos, int nchild) {
    tdpglobal = (TopDep*)calloc(sizeof(TopDep), 1);
    tdpglobal->topol = topol;
    tdpglobal->dep = dep;

    int child;

    if (nchild == 0 || nchild == 2) {
        if ((child = (dep)[pos].child0) == -1) {
            *order++ = (topol)[pos][0][0];
            *order = -1;
        } else {
            order = topolorder_lessargs(order, child);
        }
    }
    if (nchild == 1 || nchild == 2) {
        if ((child = (dep)[pos].child1) == -1) {
            *order++ = (topol)[pos][1][0];
            *order = -1;
        } else {
            order = topolorder_lessargs(order, child);
        }
    }

#if 1
    free(tdpglobal);
    tdpglobal = NULL;
#endif

    return (order);
}

static double sueff1, sueff05;

static double
cluster_mix_double(double d1, double d2) {
    return (MIN(d1, d2) * sueff1 + (d1 + d2) * sueff05);
}

static double
cluster_average_double(double d1, double d2) {
    return ((d1 + d2) * 0.5);
}

static double
cluster_minimum_double(double d1, double d2) {
    return (MIN(d1, d2));
}

void
fixed_musclesupg_double_realloc_nobk_halfmtx_memsave(aln_Opts opts, Context* ctx, int nseq, double** eff, int*** topol, double** len, Treedep* dep, int progressout, int efffree) {
    int     i, j, k, miniim, maxiim, minijm, maxijm;
    int*    intpt;
    double  tmpdouble;
    double  eff1, eff0;
    double* tmptmplen = NULL;  // static -> local, 2012/02/25
    int*    hist = NULL;  // static -> local, 2012/02/25
    Bchain* ac = NULL;  // static -> local, 2012/02/25
    int     im = -1, jm = -1;
    Bchain *acjmnext, *acjmprev;
    int     prevnode;
    Bchain* acpti;
    int *   pt1, *pt2, *pt11;
    int*    nmemar;  // static -> local, 2012/02/25
    int     nmemim, nmemjm;
    double  minscore;
    int*    nearest = NULL;  // by Mathog, a guess
    double* mindisfrom = NULL;  // by Mathog, a guess
    double (*clusterfuncpt[1])(double, double);

    sueff1 = 1 - (double)opts.sueff_global;
    sueff05 = (double)opts.sueff_global * 0.5;
    if (opts.treemethod == 'X')
        clusterfuncpt[0] = cluster_mix_double;
    else if (opts.treemethod == 'E')
        clusterfuncpt[0] = cluster_average_double;
    else if (opts.treemethod == 'q')
        clusterfuncpt[0] = cluster_minimum_double;
    else {
        reporterr("Unknown treemethod, %c\n", opts.treemethod);
        exit(1);
    }

    if (!hist) {
        hist = AllocateIntVec(ctx->njob);
        tmptmplen = AllocateFloatVec(ctx->njob);
        ac = (Bchain*)malloc(ctx->njob * sizeof(Bchain));
        nmemar = AllocateIntVec(ctx->njob);
        mindisfrom = AllocateFloatVec(ctx->njob);
        nearest = AllocateIntVec(ctx->njob);
    }

    for (i = 0; i < nseq; i++) {
        ac[i].next = ac + i + 1;
        ac[i].prev = ac + i - 1;
        ac[i].pos = i;
    }
    ac[nseq - 1].next = NULL;

    for (i = 0; i < nseq; i++)
        setnearest(ac, eff, mindisfrom + i, nearest + i, i);  // muscle

    for (i = 0; i < nseq; i++)
        tmptmplen[i] = 0.0;
    for (i = 0; i < nseq; i++) {
        hist[i] = -1;
        nmemar[i] = 1;
    }

    if (progressout)
        reporterr("\n");
    for (k = 0; k < nseq - 1; k++) {
        if (progressout && k % 10 == 0)
            reporterr("\r% 5d / %d", k, nseq);

        minscore = 999.9;
        for (acpti = ac; acpti->next != NULL; acpti = acpti->next) {
            i = acpti->pos;
            //			reporterr(       "k=%d i=%d\n", k, i );
            if (mindisfrom[i] < minscore)  // muscle
            {
                im = i;
                minscore = mindisfrom[i];
            }
        }
        jm = nearest[im];
        if (jm < im) {
            j = jm;
            jm = im;
            im = j;
        }

        prevnode = hist[im];
        if (dep)
            dep[k].child0 = prevnode;
        nmemim = nmemar[im];
        intpt = topol[k][0] = (int*)realloc(topol[k][0], (2) * sizeof(int));  // memsave
        //		intpt = topol[k][0] = (int *)realloc( topol[k][0], ( nmemim + 1 ) * sizeof( int ) );
        if (prevnode == -1) {
            *intpt++ = im;
            *intpt = -1;
        } else {
            pt1 = topol[prevnode][0];
            pt2 = topol[prevnode][1];
            if (*pt1 > *pt2) {
                pt11 = pt2;
                //				pt22 = pt1;
            } else {
                pt11 = pt1;
                //				pt22 = pt2;
            }
#if 1
            *intpt++ = *pt11;
            *intpt = -1;
#else
            for (intpt2 = pt11; *intpt2 != -1;)
                *intpt++ = *intpt2++;
            for (intpt2 = pt22; *intpt2 != -1;)
                *intpt++ = *intpt2++;
            *intpt = -1;
#endif
        }

        prevnode = hist[jm];
        if (dep)
            dep[k].child1 = prevnode;
        nmemjm = nmemar[jm];
        intpt = topol[k][1] = (int*)realloc(topol[k][1], (2) * sizeof(int));
        //		intpt = topol[k][1] = (int *)realloc( topol[k][1], ( nmemjm + 1 ) * sizeof( int ) );
        if (!intpt) {
            reporterr("Cannot reallocate topol\n");
            exit(1);
        }
        if (prevnode == -1) {
            *intpt++ = jm;
            *intpt = -1;
        } else {
            pt1 = topol[prevnode][0];
            pt2 = topol[prevnode][1];
            if (*pt1 > *pt2) {
                pt11 = pt2;
                //				pt22 = pt1;
            } else {
                pt11 = pt1;
                //				pt22 = pt2;
            }
#if 1
            *intpt++ = *pt11;
            *intpt = -1;
#else
            for (intpt2 = pt11; *intpt2 != -1;)
                *intpt++ = *intpt2++;
            for (intpt2 = pt22; *intpt2 != -1;)
                *intpt++ = *intpt2++;
            *intpt = -1;
#endif
        }

        minscore *= 0.5;

        len[k][0] = (minscore - tmptmplen[im]);
        len[k][1] = (minscore - tmptmplen[jm]);

        if (dep)
            dep[k].distfromtip = minscore;

        tmptmplen[im] = minscore;

        hist[im] = k;
        nmemar[im] = nmemim + nmemjm;

        mindisfrom[im] = 999.9;
        for (acpti = ac; acpti != NULL; acpti = acpti->next) {
            i = acpti->pos;
            if (i != im && i != jm) {
                if (i < im) {
                    miniim = i;
                    maxiim = im;
                    minijm = i;
                    maxijm = jm;
                } else if (i < jm) {
                    miniim = im;
                    maxiim = i;
                    minijm = i;
                    maxijm = jm;
                } else {
                    miniim = im;
                    maxiim = i;
                    minijm = jm;
                    maxijm = i;
                }
                eff0 = eff[miniim][maxiim - miniim];
                eff1 = eff[minijm][maxijm - minijm];
                tmpdouble = eff[miniim][maxiim - miniim] =
#if 0
				MIN( eff0, eff1 ) * sueff1 + ( eff0 + eff1 ) * sueff05;
#else
                    (clusterfuncpt[0])(eff0, eff1);
#endif
                if (tmpdouble < mindisfrom[i]) {
                    mindisfrom[i] = tmpdouble;
                    nearest[i] = im;
                }
                if (tmpdouble < mindisfrom[im]) {
                    mindisfrom[im] = tmpdouble;
                    nearest[im] = i;
                }
                if (nearest[i] == jm) {
                    nearest[i] = im;
                }
            }
        }

        //		reporterr(       "im,jm=%d,%d\n", im, jm );
        acjmprev = ac[jm].prev;
        acjmnext = ac[jm].next;
        acjmprev->next = acjmnext;
        if (acjmnext != NULL)
            acjmnext->prev = acjmprev;
        if (efffree) {
            free((void*)eff[jm]);
            eff[jm] = NULL;
        }

#if 1  // muscle seems to miss this.
        for (acpti = ac; acpti != NULL; acpti = acpti->next) {
            i = acpti->pos;
            if (nearest[i] == im) {
                if (i < im) {
                    miniim = i;
                    maxiim = im;
                } else {
                    miniim = im;
                    maxiim = i;
                }
                if (eff[miniim][maxiim - miniim] > mindisfrom[i])
                    setnearest(ac, eff, mindisfrom + i, nearest + i, i);
            }
        }
#endif

#if 0
        fprintf( stdout, "vSTEP-%03d:\n", k+1 );
		fprintf( stdout, "len0 = %f\n", len[k][0] );
        for( i=0; topol[k][0][i]>-1; i++ ) fprintf( stdout, " %03d", topol[k][0][i]+1 );
        fprintf( stdout, "\n" );
		fprintf( stdout, "len1 = %f\n", len[k][1] );
        for( i=0; topol[k][1][i]>-1; i++ ) fprintf( stdout, " %03d", topol[k][1][i]+1 );
        fprintf( stdout, "\n" );
#endif
    }
    free((void*)tmptmplen);
    tmptmplen = NULL;
    free(hist);
    hist = NULL;
    free((char*)ac);
    ac = NULL;
    free((void*)nmemar);
    nmemar = NULL;
    free(mindisfrom);
    free(nearest);
}

void
counteff_simple_double_nostatic_memsave(int nseq, int*** topol, double** len, Treedep* dep, double* node) {
    int     i, j, s1, s2;
    double  total;
    double* rootnode;
    double* eff;
    int**   localmem;
    int**   memhist;

    rootnode = AllocateDoubleVec(nseq);
    eff = AllocateDoubleVec(nseq);
    localmem = AllocateIntMtx(2, 0);
    memhist = AllocateIntMtx(nseq - 1, 0);
    for (i = 0; i < nseq - 1; i++) {
        memhist[i] = NULL;
    }

    for (i = 0; i < nseq; i++) {
        if (len[i][0] < 0.0) {
            reporterr("WARNING: negative branch length %f, step %d-0\n", len[i][0], i);
            len[i][0] = 0.0;
        }
        if (len[i][1] < 0.0) {
            reporterr("WARNING: negative branch length %f, step %d-1\n", len[i][1], i);
            len[i][1] = 0.0;
        }
    }

    for (i = 0; i < nseq; i++) {
        rootnode[i] = 0.0;
        eff[i] = 1.0;
    }

    for (i = 0; i < nseq - 1; i++) {
        if (dep[i].child0 == -1) {
            localmem[0] = calloc(sizeof(int), 2);
            localmem[0][0] = topol[i][0][0];
            localmem[0][1] = -1;
            s1 = 1;
        } else {
            localmem[0] = memhist[dep[i].child0];
            s1 = intlen(localmem[0]);
        }
        if (dep[i].child1 == -1) {
            localmem[1] = calloc(sizeof(int), 2);
            localmem[1][0] = topol[i][1][0];
            localmem[1][1] = -1;
            s2 = 1;
        } else {
            localmem[1] = memhist[dep[i].child1];
            s2 = intlen(localmem[1]);
        }

        memhist[i] = calloc(sizeof(int), s1 + s2 + 1);
        intcpy(memhist[i], localmem[0]);
        intcpy(memhist[i] + s1, localmem[1]);
        memhist[i][s1 + s2] = -1;

        for (j = 0; (s1 = localmem[0][j]) > -1; j++) {
            rootnode[s1] += (double)len[i][0] * eff[s1];
            eff[s1] *= 0.5;
            /*
           	rootnode[s1] *= 0.5;
*/
        }
        for (j = 0; (s2 = localmem[1][j]) > -1; j++) {
            rootnode[s2] += (double)len[i][1] * eff[s2];
            eff[s2] *= 0.5;
            /*
           	rootnode[s2] *= 0.5;
*/
        }
        free(localmem[0]);
        free(localmem[1]);
    }
    free(localmem);
    free(memhist[nseq - 2]);
    free(memhist);

    for (i = 0; i < nseq; i++) {
#if 1 /* 97.9.29 */
        rootnode[i] += GETA3;
#endif
#if 0
		reporterr(       "### rootnode for %d = %f\n", i, rootnode[i] );
#endif
    }
#if 1
    total = 0.0;
    for (i = 0; i < nseq; i++) {
        total += rootnode[i];
    }
#else
    total = 1.0;
#endif

    for (i = 0; i < nseq; i++) {
        node[i] = rootnode[i] / total;
    }

#if 0
	reporterr(       "weight array in counteff_simple\n" );
	for( i=0; i<nseq; i++ )
		reporterr(       "%f\n", node[i] );
	printf( "\n" );
	exit( 1 );
#endif
    free(rootnode);
    free(eff);
}

void
counteff_simple_double_nostatic(int nseq, int*** topol, double** len, double* node) {
    int     i, j, s1, s2;
    double  total;
    double* rootnode;
    double* eff;

    rootnode = AllocateDoubleVec(nseq);
    eff = AllocateDoubleVec(nseq);

    for (i = 0; i < nseq; i++)  // 2014/06/07, fu no eff wo sakeru.
    {
        if (len[i][0] < 0.0) {
            reporterr("WARNING: negative branch length %f, step %d-0\n", len[i][0], i);
            len[i][0] = 0.0;
        }
        if (len[i][1] < 0.0) {
            reporterr("WARNING: negative branch length %f, step %d-1\n", len[i][1], i);
            len[i][1] = 0.0;
        }
    }
#if DEBUG
    for (i = 0; i < nseq - 1; i++) {
        reporterr("\nstep %d, group 0\n", i);
        for (j = 0; topol[i][0][j] != -1; j++)
            reporterr("%3d ", topol[i][0][j]);
        reporterr("\n", i);
        reporterr("step %d, group 1\n", i);
        for (j = 0; topol[i][1][j] != -1; j++)
            reporterr("%3d ", topol[i][1][j]);
        reporterr("\n", i);
        reporterr("len0 = %f\n", len[i][0]);
        reporterr("len1 = %f\n", len[i][1]);
    }
#endif
    for (i = 0; i < nseq; i++) {
        rootnode[i] = 0.0;
        eff[i] = 1.0;
        /*
		rootnode[i] = 1.0;
*/
    }
    for (i = 0; i < nseq - 1; i++) {
        for (j = 0; (s1 = topol[i][0][j]) > -1; j++) {
            rootnode[s1] += (double)len[i][0] * eff[s1];
            eff[s1] *= 0.5;
            /*
           	rootnode[s1] *= 0.5;
*/
        }
        for (j = 0; (s2 = topol[i][1][j]) > -1; j++) {
            rootnode[s2] += (double)len[i][1] * eff[s2];
            eff[s2] *= 0.5;
            /*
           	rootnode[s2] *= 0.5;
*/
        }
    }
    for (i = 0; i < nseq; i++) {
#if 1 /* 97.9.29 */
        rootnode[i] += GETA3;
#endif
#if 0
		reporterr(       "### rootnode for %d = %f\n", i, rootnode[i] );
#endif
    }
#if 1
    total = 0.0;
    for (i = 0; i < nseq; i++) {
        total += rootnode[i];
    }
#else
    total = 1.0;
#endif

    for (i = 0; i < nseq; i++) {
        node[i] = rootnode[i] / total;
    }

#if 0
	reporterr(       "weight array in counteff_simple\n" );
	for( i=0; i<nseq; i++ )
		reporterr(       "%f\n", node[i] );
	printf( "\n" );
	exit( 1 );
#endif
    free(rootnode);
    free(eff);
}

void
counteff_simple(int nseq, int*** topol, double** len, double* node) {
    int    i, j, s1, s2;
    double total;

    double* rootnode;
    double* eff;
    rootnode = AllocateDoubleVec(nseq);
    eff = AllocateDoubleVec(nseq);

#if DEBUG
    for (i = 0; i < nseq; i++) {
        reporterr("len0 = %f\n", len[i][0]);
        reporterr("len1 = %f\n", len[i][1]);
    }
#endif
    for (i = 0; i < nseq; i++) {
        rootnode[i] = 0.0;
        eff[i] = 1.0;
        /*
		rootnode[i] = 1.0;
*/
    }
    for (i = 0; i < nseq - 1; i++) {
        for (j = 0; (s1 = topol[i][0][j]) > -1; j++) {
            rootnode[s1] += len[i][0] * eff[s1];
            eff[s1] *= 0.5;
            /*
           	rootnode[s1] *= 0.5;
*/
        }
        for (j = 0; (s2 = topol[i][1][j]) > -1; j++) {
            rootnode[s2] += len[i][1] * eff[s2];
            eff[s2] *= 0.5;
            /*
           	rootnode[s2] *= 0.5;
*/
        }
    }
    for (i = 0; i < nseq; i++) {
#if 1 /* 97.9.29 */
        rootnode[i] += GETA3;
#endif
#if 0
		reporterr(       "### rootnode for %d = %f\n", i, rootnode[i] );
#endif
    }
#if 1
    total = 0.0;
    for (i = 0; i < nseq; i++) {
        total += rootnode[i];
    }
#else
    total = 1.0;
#endif

    for (i = 0; i < nseq; i++) {
        node[i] = rootnode[i] / total;
    }

#if 0
	reporterr(       "weight array in counteff_simple\n" );
	for( i=0; i<nseq; i++ )
		reporterr(       "%f\n", node[i] );
	printf( "\n" );
	exit( 1 );
#endif
#if 1
    free(rootnode);
    free(eff);
#endif
}

void
copyWithNoGaps(char* aseq, const char* seq) {
    for (; *seq != 0; seq++) {
        if (*seq != '-')
            *aseq++ = *seq;
    }
    *aseq = 0;
}

int
isallgap(char* seq) {
    for (; *seq != 0; seq++) {
        if (*seq != '-')
            return (0);
    }
    return (1);
}

void
gappick(int nseq, int s, char** aseq, char** mseq2, double** eff, double* effarr) {
    int i, j, count, countjob, len, allgap;
    len = strlen(aseq[0]);
    for (i = 0, count = 0; i < len; i++) {
        allgap = 1;
        for (j = 0; j < nseq; j++)
            if (j != s)
                allgap *= (aseq[j][i] == '-');
        if (allgap == 0) {
            for (j = 0, countjob = 0; j < nseq; j++) {
                if (j != s) {
                    mseq2[countjob][count] = aseq[j][i];
                    countjob++;
                }
            }
            count++;
        }
    }
    for (i = 0; i < nseq - 1; i++)
        mseq2[i][count] = 0;

    for (i = 0, countjob = 0; i < nseq; i++) {
        if (i != s) {
            effarr[countjob] = eff[s][i];
            countjob++;
        }
    }
    /*
fprintf( stdout, "effarr in gappick s = %d\n", s+1 );
for( i=0; i<countjob; i++ ) 
	fprintf( stdout, " %f", effarr[i] );
printf( "\n" );
*/
}

void
commongappick_record(int nseq, char** seq, int* map) {
    int i, j, count;
    int len = strlen(seq[0]);

    for (i = 0, count = 0; i <= len; i++) {
        /*
		allgap = 1;
		for( j=0; j<nseq; j++ ) 
			allgap *= ( seq[j][i] == '-' );
		if( !allgap )
	*/
        for (j = 0; j < nseq; j++)
            if (seq[j][i] != '-')
                break;
        if (j != nseq) {
            for (j = 0; j < nseq; j++) {
                seq[j][count] = seq[j][i];
            }
            map[count] = i;
            count++;
        }
    }
}

void
doublencpy(double* vec1, double* vec2, int len) {
    while (len--)
        *vec1++ = *vec2++;
}

#define SEGMENTSIZE 150

void
dontcalcimportance_half(Context* ctx, int nseq, char** seq, LocalHom** localhom) {
    int       i, j;
    LocalHom* ptr;
    int*      nogaplen;

    nogaplen = AllocateIntVec(nseq);

    for (i = 0; i < nseq; i++) {
        nogaplen[i] = seqlen(ctx, seq[i]);
    }

    for (i = 0; i < nseq; i++) {
        for (j = 0; j < nseq; j++) {
            if (i >= j)
                continue;
            for (ptr = localhom[i] + j - i; ptr; ptr = ptr->next) {
//				reporterr(       "i,j=%d,%d,ptr=%p\n", i, j, ptr );
#if 1
                //				ptr->importance = ptr->opt / ptr->overlapaa;
                ptr->importance = ptr->opt;
//				ptr->fimportance = (double)ptr->importance;
#else
                ptr->importance = ptr->opt / MIN(nogaplen[i], nogaplen[j]);
#endif
            }
        }
    }
    free(nogaplen);
}

void
dontcalcimportance_firstone(int nseq, LocalHom** localhom) {
    int       i, j, nseq1;
    LocalHom* ptr;

    nseq1 = nseq - 1;
    for (i = 0; i < nseq1; i++) {
        j = 0;
        {
            for (ptr = localhom[i] + j; ptr; ptr = ptr->next) {
//				reporterr(       "i,j=%d,%d,ptr=%p\n", i, j, ptr );
#if 1
                //				ptr->importance = ptr->opt / ptr->overlapaa;
                ptr->importance = ptr->opt * 0.5;  // tekitou
//				ptr->fimportance = (double)ptr->importance;
//				reporterr(       "i=%d, j=%d, importance = %f, opt=%f\n", i, j, ptr->fimportance, ptr->opt  );
#else
                ptr->importance = ptr->opt / MIN(nogaplen[i], nogaplen[j]);
#endif
            }
        }
    }
#if 1
#else
    free(nogaplen);
#endif
}

void
calcimportance_target(Context* ctx, int nseq, int ntarget, double* eff, char** seq, LocalHom** localhom, int* targetmap, int* targetmapr, int alloclen) {
    int       i, j, pos, len, ti, tj;
    double*   importance;  // static -> local, 2012/02/25
    double    tmpdouble;
    double *  ieff, totaleff;  // counteff_simple_double ni utsusu kamo
    int*      nogaplen;  // static -> local, 2012/02/25
    LocalHom* tmpptr;

    importance = AllocateDoubleVec(alloclen);
    nogaplen = AllocateIntVec(nseq);
    ieff = AllocateDoubleVec(nseq);

    totaleff = 0.0;
    for (i = 0; i < nseq; i++) {
        nogaplen[i] = seqlen(ctx, seq[i]);
        if (nogaplen[i] == 0)
            ieff[i] = 0.0;
        else
            ieff[i] = eff[i];
        totaleff += ieff[i];
    }
    for (i = 0; i < nseq; i++)
        ieff[i] /= totaleff;
    for (i = 0; i < nseq; i++)
        printf("eff[%d] = %30.25f\n", i, ieff[i]);

#if 0
	for( i=0; i<nseq; i++ ) for( j=0; j<nseq; j++ )
	{
		tmpptr = localhom[i]+j;
		reporterr(       "%d-%d\n", i, j );
		do
		{
			reporterr(       "reg1=%d-%d, reg2=%d-%d, opt=%f\n", tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, tmpptr->opt );
		} while( tmpptr=tmpptr->next );
	}
#endif

    //	for( i=0; i<nseq; i++ )
    for (ti = 0; ti < ntarget; ti++) {
        i = targetmapr[ti];
        //		reporterr(       "i = %d\n", i );
        for (pos = 0; pos < alloclen; pos++)
            importance[pos] = 0.0;
        for (j = 0; j < nseq; j++) {
            if (i == j)
                continue;
            //			tmpptr = localhom[ti]+j;
            for (tmpptr = localhom[ti] + j; tmpptr; tmpptr = tmpptr->next) {
                if (tmpptr->opt == -1)
                    continue;
                for (pos = tmpptr->start1; pos <= tmpptr->end1; pos++) {
#if 1
                    //					if( pos == 0 ) reporterr( "hit! i=%d, j=%d, pos=%d\n", i, j, pos );
                    importance[pos] += ieff[j];
#else
                    importance[pos] += ieff[j] * tmpptr->opt / MIN(nogaplen[i], nogaplen[j]);
                    importance[pos] += ieff[j] * tmpptr->opt / tmpptr->overlapaa;
#endif
                }
            }
        }
#if 0
		reporterr(       "position specific importance of seq %d:\n", i );
		for( pos=0; pos<alloclen; pos++ )
			reporterr(       "%d: %f\n", pos, importance[pos] );
		reporterr(       "\n" );
#endif
        for (j = 0; j < nseq; j++) {
            //			reporterr(       "i=%d, j=%d\n", i, j );
            if (i == j)
                continue;
            if (localhom[ti][j].opt == -1.0)
                continue;
#if 1
            for (tmpptr = localhom[ti] + j; tmpptr; tmpptr = tmpptr->next) {
                if (tmpptr->opt == -1.0)
                    continue;
                tmpdouble = 0.0;
                len = 0;
                for (pos = tmpptr->start1; pos <= tmpptr->end1; pos++) {
                    tmpdouble += importance[pos];
                    len++;
                }

                tmpdouble /= (double)len;

                tmpptr->importance = tmpdouble * tmpptr->opt;
                //				tmpptr->fimportance = (double)tmpptr->importance;
            }
#else
            tmpdouble = 0.0;
            len = 0;
            for (tmpptr = localhom[ti] + j; tmpptr; tmpptr = tmpptr->next) {
                if (tmpptr->opt == -1.0)
                    continue;
                for (pos = tmpptr->start1; pos <= tmpptr->end1; pos++) {
                    tmpdouble += importance[pos];
                    len++;
                }
            }

            tmpdouble /= (double)len;

            for (tmpptr = localhom[ti] + j; tmpptr; tmpptr = tmpptr->next) {
                if (tmpptr->opt == -1.0)
                    continue;
                tmpptr->importance = tmpdouble * tmpptr->opt;
                //				tmpptr->importance = tmpptr->opt / tmpptr->overlapaa; //$B$J$+$C$?$3$H$K$9$k(B
            }
#endif

            //			reporterr(       "importance of match between %d - %d = %f\n", i, j, tmpdouble );
        }
    }

#if 0
	printf(       "before averaging:\n" );

	for( ti=0; ti<ntarget; ti++ ) for( j=0; j<nseq; j++ )
	{
		i = targetmapr[ti];
		if( i == j ) continue;
		printf(       "%d-%d\n", i, j );
		for( tmpptr = localhom[ti]+j; tmpptr; tmpptr=tmpptr->next )
		{
			printf(       "reg1=%d-%d, reg2=%d-%d, imp=%f -> %f opt=%30.25f\n", tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, tmpptr->opt / tmpptr->overlapaa, eff[i] * tmpptr->importance, tmpptr->opt );
		}
	}
#endif

#if 1
    //	reporterr(       "average?\n" );
    //	for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ )
    for (ti = 0; ti < ntarget; ti++)
        for (tj = ti + 1; tj < ntarget; tj++) {
            double    imp;
            LocalHom *tmpptr1, *tmpptr2;

            i = targetmapr[ti];
            j = targetmapr[tj];
            //		if( i == j ) continue;

            //		reporterr(       "i=%d, j=%d\n", i, j );

            tmpptr1 = localhom[ti] + j;
            tmpptr2 = localhom[tj] + i;
            for (; tmpptr1 && tmpptr2; tmpptr1 = tmpptr1->next, tmpptr2 = tmpptr2->next) {
                if (tmpptr1->opt == -1.0 || tmpptr2->opt == -1.0) {
                    //				reporterr(       "WARNING: i=%d, j=%d, tmpptr1->opt=%f, tmpptr2->opt=%f\n", i, j, tmpptr1->opt, tmpptr2->opt );
                    continue;
                }
                //			reporterr(       "## importances = %f, %f\n", tmpptr1->importance, tmpptr2->importance );
                imp = 0.5 * (tmpptr1->importance + tmpptr2->importance);
                tmpptr1->importance = tmpptr2->importance = imp;
                //			tmpptr1->fimportance = tmpptr2->fimportance = (double)imp;

                //			reporterr(       "## importance = %f\n", tmpptr1->importance );
            }

#if 0  // commented out, 2012/02/10
		if( ( tmpptr1 && !tmpptr2 ) || ( !tmpptr1 && tmpptr2 ) )
		{
			reporterr(       "ERROR: i=%d, j=%d\n", i, j );
			exit( 1 );
		}
#endif
        }

    for (ti = 0; ti < ntarget; ti++)
        for (j = 0; j < nseq; j++) {
            double    imp;
            LocalHom* tmpptr1;

            i = targetmapr[ti];
            if (i == j)
                continue;
            if (targetmap[j] != -1)
                continue;

            //		reporterr(       "i=%d, j=%d\n", i, j );

            tmpptr1 = localhom[ti] + j;
            for (; tmpptr1; tmpptr1 = tmpptr1->next) {
                if (tmpptr1->opt == -1.0) {
                    //				reporterr(       "WARNING: i=%d, j=%d, tmpptr1->opt=%f, tmpptr2->opt=%f\n", i, j, tmpptr1->opt, tmpptr2->opt );
                    continue;
                }
                //			reporterr(       "## importances = %f, %f\n", tmpptr1->importance, tmpptr2->importance );
                imp = 0.5 * (tmpptr1->importance);
                //			imp = 1.0 * ( tmpptr1->importance );
                tmpptr1->importance = imp;
                //			tmpptr1->fimportance = (double)imp;

                //			reporterr(       "## importance = %f\n", tmpptr1->importance );
            }

#if 0  // commented out, 2012/02/10
		if( ( tmpptr1 && !tmpptr2 ) || ( !tmpptr1 && tmpptr2 ) )
		{
			reporterr(       "ERROR: i=%d, j=%d\n", i, j );
			exit( 1 );
		}
#endif
        }
#endif
#if 0
	printf(       "after averaging:\n" );

	for( ti=0; ti<ntarget; ti++ ) for( j=0; j<nseq; j++ )
	{
		i = targetmapr[ti];
		if( i == j ) continue;
		printf(       "%d-%d\n", i, j );
		for( tmpptr = localhom[ti]+j; tmpptr; tmpptr=tmpptr->next )
		{
			if( tmpptr->end1 )
				printf(       "%d-%d, reg1=%d-%d, reg2=%d-%d, imp=%f -> %f opt=%f\n", i, j, tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, tmpptr->opt / tmpptr->overlapaa, tmpptr->importance, tmpptr->opt );
		}
	}
//exit( 1 );
#endif
    free(importance);
    free(nogaplen);
    free(ieff);
}

void
calcimportance_half(Context* ctx, int nseq, double* eff, char** seq, LocalHom** localhom, int alloclen) {
    int       i, j, pos, len;
    double*   importance;  // static -> local, 2012/02/25
    double    tmpdouble;
    double *  ieff, totaleff;  // counteff_simple_double ni utsusu kamo
    int*      nogaplen;  // static -> local, 2012/02/25
    LocalHom* tmpptr;

    importance = AllocateDoubleVec(alloclen);
    //	reporterr("alloclen=%d, nlenmax=%d\n", alloclen, nlenmax );
    nogaplen = AllocateIntVec(nseq);
    ieff = AllocateDoubleVec(nseq);

    totaleff = 0.0;
    for (i = 0; i < nseq; i++) {
        nogaplen[i] = seqlen(ctx, seq[i]);
        //		reporterr(       "nogaplen[] = %d\n", nogaplen[i] );
        if (nogaplen[i] == 0)
            ieff[i] = 0.0;
        else
            ieff[i] = eff[i];
        totaleff += ieff[i];
    }
    for (i = 0; i < nseq; i++)
        ieff[i] /= totaleff;
        //	for( i=0; i<nseq; i++ ) reporterr(       "eff[%d] = %f\n", i, ieff[i] );

#if 0
	for( i=0; i<nseq; i++ ) for( j=0; j<nseq; j++ )
	{
		tmpptr = localhom[i]+j;
		reporterr(       "%d-%d\n", i, j );
		do
		{
			reporterr(       "reg1=%d-%d, reg2=%d-%d, opt=%f\n", tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, tmpptr->opt );
		} while( tmpptr=tmpptr->next );
	}
#endif

    for (i = 0; i < nseq; i++) {
        //		reporterr(       "i = %d\n", i );
        for (pos = 0; pos < alloclen; pos++) {
            importance[pos] = 0.0;
        }
        for (j = 0; j < nseq; j++) {
            if (i == j)
                continue;

            else if (i < j) {
                for (tmpptr = localhom[i] + j - i; tmpptr; tmpptr = tmpptr->next) {
                    //					reporterr( "pos=%d, alloclen=%d\n", pos, alloclen );
                    if (tmpptr->opt == -1)
                        continue;
                    for (pos = tmpptr->start1; pos <= tmpptr->end1; pos++) {
#if 1
                        //						if( pos == 0 ) reporterr( "hit! i=%d, j=%d, pos=%d\n", i, j, pos );
                        importance[pos] += ieff[j];
#else
                        importance[pos] += ieff[j] * tmpptr->opt / MIN(nogaplen[i], nogaplen[j]);
                        importance[pos] += ieff[j] * tmpptr->opt / tmpptr->overlapaa;
#endif
                    }
                }
            } else {
                for (tmpptr = localhom[j] + i - j; tmpptr; tmpptr = tmpptr->next) {
                    if (tmpptr->opt == -1)
                        continue;
                    for (pos = tmpptr->start2; pos <= tmpptr->end2; pos++) {
#if 1
                        //						if( pos == 0 ) reporterr( "hit! i=%d, j=%d, pos=%d\n", i, j, pos );
                        importance[pos] += ieff[j];
#else
                        importance[pos] += ieff[j] * tmpptr->opt / MIN(nogaplen[i], nogaplen[j]);
                        importance[pos] += ieff[j] * tmpptr->opt / tmpptr->overlapaa;
#endif
                    }
                }
            }
        }
#if 0
		reporterr(       "position specific importance of seq %d:\n", i );
		for( pos=0; pos<nlenmax; pos++ )
			reporterr(       "%d: %f\n", pos, importance[pos] );
		reporterr(       "\n" );
#endif
        for (j = 0; j < nseq; j++) {
            //			reporterr(       "i=%d, j=%d\n", i, j );
            if (i == j)
                continue;

            else if (i < j) {
                if (localhom[i][j - i].opt == -1.0)
                    continue;

                for (tmpptr = localhom[i] + j - i; tmpptr; tmpptr = tmpptr->next) {
                    if (tmpptr->opt == -1.0)
                        continue;
                    tmpdouble = 0.0;
                    len = 0;
                    for (pos = tmpptr->start1; pos <= tmpptr->end1; pos++) {
                        tmpdouble += importance[pos];
                        len++;
                    }

                    tmpdouble /= (double)len;

                    tmpptr->importance = tmpdouble * tmpptr->opt;
                    //					tmpptr->fimportance = (double)tmpptr->importance;
                }
            } else {
                if (localhom[j][i - j].opt == -1.0)
                    continue;

                for (tmpptr = localhom[j] + i - j; tmpptr; tmpptr = tmpptr->next) {
                    if (tmpptr->opt == -1.0)
                        continue;
                    tmpdouble = 0.0;
                    len = 0;
                    for (pos = tmpptr->start2; pos <= tmpptr->end2; pos++) {
                        tmpdouble += importance[pos];
                        len++;
                    }

                    tmpdouble /= (double)len;

                    tmpptr->rimportance = tmpdouble * tmpptr->opt;
                    //					tmpptr->fimportance = (double)tmpptr->importance;
                }
            }

            //			reporterr(       "importance of match between %d - %d = %f\n", i, j, tmpdouble );
        }
    }

#if 0
	printf(       "before averaging:\n" );

	for( i=0; i<nseq; i++ ) for( j=0; j<nseq; j++ )
	{
		if( i == j ) continue;

		else if( i < j ) 
		{
			printf(       "%d-%d\n", i, j );
			for( tmpptr = localhom[i]+j-i; tmpptr; tmpptr=tmpptr->next )
			{
				printf(       "reg1=%d-%d, reg2=%d-%d, imp=%f -> %f opt=%f\n", tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, tmpptr->opt / tmpptr->overlapaa, eff[i] * tmpptr->importance, tmpptr->opt );
			}
		}
		else
		{
			printf(       "%d-%d\n", i, j );
			for( tmpptr = localhom[j]+i-j; tmpptr; tmpptr=tmpptr->next )
			{
				printf(       "reg1=%d-%d, reg2=%d-%d, imp=%f -> %f opt=%f\n", tmpptr->start2, tmpptr->end2, tmpptr->start1, tmpptr->end1, tmpptr->opt / tmpptr->overlapaa, eff[i] * tmpptr->rimportance, tmpptr->opt );
			}
		}
	}
#endif

#if 1
    //	reporterr(       "average?\n" );
    for (i = 0; i < nseq - 1; i++)
        for (j = i + 1; j < nseq; j++) {
            double    imp;
            LocalHom* tmpptr1;

            //		reporterr(       "i=%d, j=%d\n", i, j );

            tmpptr1 = localhom[i] + j - i;
            for (; tmpptr1; tmpptr1 = tmpptr1->next) {
                if (tmpptr1->opt == -1.0) {
                    //				reporterr(       "WARNING: i=%d, j=%d, tmpptr1->opt=%f, tmpptr2->opt=%f\n", i, j, tmpptr1->opt, tmpptr2->opt );
                    continue;
                }
                //			reporterr(       "## importances = %f, %f\n", tmpptr1->importance, tmpptr2->importance );
                imp = 0.5 * (tmpptr1->importance + tmpptr1->rimportance);
                tmpptr1->importance = tmpptr1->rimportance = imp;
                //			tmpptr1->fimportance = tmpptr2->fimportance = (double)imp;

                //			reporterr(       "## importance = %f\n", tmpptr1->importance );
            }

#if 0  // commented out, 2012/02/10
		if( ( tmpptr1 && !tmpptr2 ) || ( !tmpptr1 && tmpptr2 ) )
		{
			reporterr(       "ERROR: i=%d, j=%d\n", i, j );
			exit( 1 );
		}
#endif
        }
#endif
#if 0
	printf(       "after averaging:\n" );

	for( i=0; i<nseq; i++ ) for( j=0; j<nseq; j++ )
	{
		if( i < j ) for( tmpptr = localhom[i]+j-i; tmpptr; tmpptr=tmpptr->next )
		{
			if( tmpptr->end1 && tmpptr->start1 != -1 )
				printf(       "%d-%d, reg1=%d-%d, reg2=%d-%d, imp=%f -> %f opt=%f\n", i, j, tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, tmpptr->opt / tmpptr->overlapaa, tmpptr->importance, tmpptr->opt );
		}
		else for( tmpptr = localhom[j]+i-j; tmpptr; tmpptr=tmpptr->next )
		{
			if( tmpptr->end2 && tmpptr->start2 != -1 )
				printf(       "%d-%d, reg1=%d-%d, reg2=%d-%d, imp=%f -> %f opt=%f\n", i, j, tmpptr->start2, tmpptr->end2, tmpptr->start1, tmpptr->end1, tmpptr->opt / tmpptr->overlapaa, tmpptr->importance, tmpptr->opt );
		}
	}
exit( 1 );
#endif
    free(importance);
    free(nogaplen);
    free(ieff);
}

void
gapireru(Context* ctx, char* res, char* ori, char* gt) {
    char g;
    char gapchar = *ctx->newgapstr;
    while ((g = *gt++)) {
        if (g == '-') {
            *res++ = gapchar;
        } else {
            *res++ = *ori++;
        }
    }
    *res = 0;
}

void
getkyokaigap(char* g, char** s, int pos, int n) {
    //	char *bk = g;
    //	while( n-- ) *g++ = '-';
    while (n--)
        *g++ = (*s++)[pos];

    //	reporterr(       "bk = %s\n", bk );
}

void
new_OpeningGapCount(double* ogcp, int clus, char** seq, double* eff, int len, char* sgappat)
#if 0
{
	int i, j, gc, gb; 
	double feff;

	
	for( i=0; i<len+1; i++ ) ogcp[i] = 0.0;
	for( j=0; j<clus; j++ ) 
	{
		feff = (double)eff[j];
		gc = ( sgappat[j] == '-' );
		for( i=0; i<len; i++ ) 
		{
			gb = gc;
			gc = ( seq[j][i] == '-' );
			if( !gb *  gc ) ogcp[i] += feff;
		}
	}
}
#else
{
    int     i, j, gc, gb;
    double  feff;
    double* fpt;
    char*   spt;

    fpt = ogcp;
    i = len;
    while (i--)
        *fpt++ = 0.0;
    for (j = 0; j < clus; j++) {
        feff = (double)eff[j];
        spt = seq[j];
        fpt = ogcp;
        gc = (sgappat[j] == '-');
        i = len;
        while (i--) {
            gb = gc;
            gc = (*spt++ == '-');
            {
                if (!gb * gc)
                    *fpt += feff;
                fpt++;
            }
        }
    }
}
#endif
void
new_OpeningGapCount_zure(double* ogcp, int clus, char** seq, double* eff, int len, char* sgappat, char* egappat)
#if 0
{
	int i, j, gc, gb; 
	double feff;

	
	for( i=0; i<len+1; i++ ) ogcp[i] = 0.0;
	for( j=0; j<clus; j++ ) 
	{
		feff = (double)eff[j];
		gc = ( sgappat[j] == '-' );
		for( i=0; i<len; i++ ) 
		{
			gb = gc;
			gc = ( seq[j][i] == '-' );
			if( !gb *  gc ) ogcp[i] += feff;
		}
		{
			gb = gc;
			gc = ( egappat[j] == '-' );
			if( !gb *  gc ) ogcp[i] += feff;
		}
	}
}
#else
{
    int     i, j, gc, gb;
    double  feff;
    double* fpt;
    char*   spt;

    fpt = ogcp;
    i = len + 2;
    while (i--)
        *fpt++ = 0.0;
    for (j = 0; j < clus; j++) {
        feff = (double)eff[j];
        spt = seq[j];
        fpt = ogcp;
        gc = (sgappat[j] == '-');
        i = len;
        while (i--) {
            gb = gc;
            gc = (*spt++ == '-');
            {
                if (!gb * gc)
                    *fpt += feff;
                fpt++;
            }
        }
        {
            gb = gc;
            gc = (egappat[j] == '-');
            if (!gb * gc)
                *fpt += feff;
        }
    }
}
#endif

void
new_FinalGapCount_zure(double* fgcp, int clus, char** seq, double* eff, int len, char* sgappat, char* egappat)
#if 0
{
	int i, j, gc, gb; 
	double feff;
	
	for( i=0; i<len+1; i++ ) fgcp[i] = 0.0;
	for( j=0; j<clus; j++ ) 
	{
		feff = (double)eff[j];
		gc = ( sgappat[j] == '-' );
		for( i=0; i<len; i++ ) 
		{
			gb = gc;
			gc = ( seq[j][i] == '-' );
			{
				if( gb * !gc ) fgcp[i] += feff;
			}
		}
		{
			gb = gc;
			gc = ( egappat[j] == '-' );
			{
				if( gb * !gc ) fgcp[len] += feff;
			}
		}
	}
}
#else
{
    int     i, j, gc, gb;
    double  feff;
    double* fpt;
    char*   spt;

    fpt = fgcp;
    i = len + 2;
    while (i--)
        *fpt++ = 0.0;
    for (j = 0; j < clus; j++) {
        feff = (double)eff[j];
        fpt = fgcp;
        spt = seq[j];
        gc = (sgappat[j] == '-');
        i = len;
        while (i--) {
            gb = gc;
            gc = (*spt++ == '-');
            {
                if (gb * !gc)
                    *fpt += feff;
                fpt++;
            }
        }
        {
            gb = gc;
            gc = (egappat[j] == '-');
            {
                if (gb * !gc)
                    *fpt += feff;
            }
        }
    }
}
#endif
void
new_FinalGapCount(double* fgcp, int clus, char** seq, double* eff, int len, char* egappat)
#if 0
{
	int i, j, gc, gb; 
	double feff;
	
	for( i=0; i<len; i++ ) fgcp[i] = 0.0;
	for( j=0; j<clus; j++ ) 
	{
		feff = (double)eff[j];
		gc = ( seq[j][0] == '-' );
		for( i=1; i<len; i++ ) 
		{
			gb = gc;
			gc = ( seq[j][i] == '-' );
			{
				if( gb * !gc ) fgcp[i-1] += feff;
			}
		}
		{
			gb = gc;
			gc = ( egappat[j] == '-' );
			{
				if( gb * !gc ) fgcp[len-1] += feff;
			}
		}
	}
}
#else
{
    int     i, j, gc, gb;
    double  feff;
    double* fpt;
    char*   spt;

    fpt = fgcp;
    i = len;
    while (i--)
        *fpt++ = 0.0;
    for (j = 0; j < clus; j++) {
        feff = (double)eff[j];
        fpt = fgcp;
        spt = seq[j];
        gc = (*spt == '-');
        i = len;
        while (i--) {
            gb = gc;
            gc = (*++spt == '-');
            {
                if (gb * !gc)
                    *fpt += feff;
                fpt++;
            }
        }
        {
            gb = gc;
            gc = (egappat[j] == '-');
            {
                if (gb * !gc)
                    *fpt += feff;
            }
        }
    }
}
#endif

void
st_OpeningGapAdd(double* ogcp, int clus, char** seq, double* eff, int len) {
    int     i, j, gc, gb;
    double* fpt;
    char*   spt;
    int     newmem = clus - 1;
    double  neweff = eff[newmem];
    double  orieff = 1.0 - neweff;
    double  feff;

    //	fpt = ogcp;
    //	i = len;
    //	while( i-- ) *fpt++ = 0.0;

    j = clus - 1;
    //	for( j=0; j<clus; j++ )
    {
        feff = (double)eff[j];
        spt = seq[j];
        fpt = ogcp;
        i = len;
        gc = 0;
        while (i--) {
            gb = gc;
            gc = (*spt++ == '-');
            *fpt *= orieff;
            if (!gb * gc)
                *fpt += feff;
            fpt++;
        }
    }
    ogcp[len] = 0.0;

#if 0
	for( i=0; i<len; i++ )
		reporterr( "ogcp[%d]=%f\n", i, ogcp[i] );
	for( i=0; i<clus; i++ )
		reporterr( "%s\n", seq[i] );
	exit( 1 );
#endif
}

void
st_OpeningGapCount(double* ogcp, int clus, char** seq, double* eff, int len) {
    int     i, j, gc, gb;
    double  feff;
    double* fpt;
    char*   spt;

    fpt = ogcp;
    i = len;
    while (i--)
        *fpt++ = 0.0;
    for (j = 0; j < clus; j++) {
        feff = (double)eff[j];
        spt = seq[j];
        fpt = ogcp;
        gc = 0;
        //		gc = 1;
        i = len;
        while (i--) {
            gb = gc;
            gc = (*spt++ == '-');
            {
                if (!gb * gc)
                    *fpt += feff;
                fpt++;
            }
        }
    }
    ogcp[len] = 0.0;
}

void
st_FinalGapCount_zure(double* fgcp, int clus, char** seq, double* eff, int len) {
    int     i, j, gc, gb;
    double  feff;
    double* fpt;
    char*   spt;

    fpt = fgcp;
    i = len + 1;
    while (i--)
        *fpt++ = 0.0;
    for (j = 0; j < clus; j++) {
        feff = (double)eff[j];
        fpt = fgcp + 1;
        spt = seq[j];
        gc = (*spt == '-');
        i = len;
        //		for( i=1; i<len; i++ )
        while (i--) {
            gb = gc;
            gc = (*++spt == '-');
            {
                if (gb * !gc)
                    *fpt += feff;
                fpt++;
            }
        }
        {
            gb = gc;
            gc = 0;
            //			gc = 1;
            {
                if (gb * !gc)
                    *fpt += feff;
            }
        }
    }
}

void
st_FinalGapAdd(double* fgcp, int clus, char** seq, double* eff, int len) {
    int     i, j, gc, gb;
    double* fpt;
    char*   spt;
    int     newmem = clus - 1;
    double  neweff = eff[newmem];
    double  orieff = 1.0 - neweff;
    double  feff;

    //	fpt = fgcp;
    //	i = len;
    //	while( i-- ) *fpt++ = 0.0;

    j = clus - 1;
    //	for( j=0; j<clus; j++ )
    {
        feff = (double)eff[j];
        fpt = fgcp;
        spt = seq[j];
        gc = (*spt == '-');
        i = len;
        //		for( i=1; i<len; i++ )
        while (i--) {
            *fpt *= orieff;
            gb = gc;
            gc = (*++spt == '-');
            {
                if (gb * !gc)
                    *fpt += feff;
                fpt++;
            }
        }
        {
            *fpt *= orieff;
            gb = gc;
            gc = 0;
            //			gc = 1;
            {
                if (gb * !gc)
                    *fpt += feff;
            }
        }
    }
}

void
st_FinalGapCount(double* fgcp, int clus, char** seq, double* eff, int len) {
    int     i, j, gc, gb;
    double  feff;
    double* fpt;
    char*   spt;

    fpt = fgcp;
    i = len;
    while (i--)
        *fpt++ = 0.0;
    for (j = 0; j < clus; j++) {
        feff = (double)eff[j];
        fpt = fgcp;
        spt = seq[j];
        gc = (*spt == '-');
        i = len;
        //		for( i=1; i<len; i++ )
        while (i--) {
            gb = gc;
            gc = (*++spt == '-');
            {
                if (gb * !gc)
                    *fpt += feff;
                fpt++;
            }
        }
        {
            gb = gc;
            gc = 0;
            //			gc = 1;
            {
                if (gb * !gc)
                    *fpt += feff;
            }
        }
    }
}

void
getGapPattern(double* fgcp, int clus, char** seq, double* eff, int len) {
    int     i, j, gc, gb;
    double  feff;
    double* fpt;
    char*   spt;

    fpt = fgcp;
    i = len + 1;
    while (i--)
        *fpt++ = 0.0;
    for (j = 0; j < clus; j++) {
        feff = (double)eff[j];
        fpt = fgcp;
        spt = seq[j];
        gc = (*spt == '-');
        i = len + 1;
        while (i--) {
            gb = gc;
            gc = (*++spt == '-');
            {
                if (gb * !gc)
                    *fpt += feff;
                fpt++;
            }
        }
#if 0
		{
			gb = gc;
			gc = ( egappat[j] == '-' );
			{
				if( gb * !gc ) *fpt += feff;
			}
		}
#endif
    }
    for (j = 0; j < len; j++) {
        reporterr("%d, %f\n", j, fgcp[j]);
    }
}

void
getdigapfreq_st(double* freq, int clus, char** seq, double* eff, int len) {
    int    i, j;
    double feff;
    for (i = 0; i < len + 1; i++)
        freq[i] = 0.0;
    for (i = 0; i < clus; i++) {
        feff = eff[i];
        if (0 && seq[i][0] == '-')  // machigai kamo
            freq[0] += feff;
        for (j = 1; j < len; j++) {
            if (seq[i][j] == '-' && seq[i][j - 1] == '-')
                freq[j] += feff;
        }
        if (0 && seq[i][len - 1] == '-')
            freq[len] += feff;
    }
    //	reporterr(       "\ndigapf = \n" );
    //	for( i=0; i<len+1; i++ ) reporterr(       "%5.3f ", freq[i] );
}

void
getdiaminofreq_x(double* freq, int clus, char** seq, double* eff, int len) {
    int    i, j;
    double feff;
    for (i = 0; i < len + 2; i++)
        freq[i] = 0.0;
    for (i = 0; i < clus; i++) {
        feff = eff[i];
        if (seq[i][0] != '-')  // tadashii
            freq[0] += feff;
        for (j = 1; j < len; j++) {
            if (seq[i][j] != '-' && seq[i][j - 1] != '-')
                freq[j] += feff;
        }
        if (1 && seq[i][len - 1] != '-')  // xxx wo tsukawanaitoki [len-1] nomi
            freq[len] += feff;
    }
    //	reporterr(       "\ndiaaf = \n" );
    //	for( i=0; i<len+1; i++ ) reporterr(       "%5.3f ", freq[i] );
}

void
getdiaminofreq_st(double* freq, int clus, char** seq, double* eff, int len) {
    int    i, j;
    double feff;
    for (i = 0; i < len + 1; i++)
        freq[i] = 0.0;
    for (i = 0; i < clus; i++) {
        feff = eff[i];
        if (seq[i][0] != '-')
            freq[0] += feff;
        for (j = 1; j < len; j++) {
            if (seq[i][j] != '-' && seq[i][j - 1] != '-')
                freq[j] += feff;
        }
        //		if( 1 && seq[i][len-1] != '-' ) // xxx wo tsukawanaitoki [len-1] nomi
        freq[len] += feff;
    }
    //	reporterr(       "\ndiaaf = \n" );
    //	for( i=0; i<len+1; i++ ) reporterr(       "%5.3f ", freq[i] );
}

void
getdigapfreq_part(double* freq, int clus, char** seq, double* eff, int len, char* sgappat, char* egappat) {
    int    i, j;
    double feff;
    for (i = 0; i < len + 2; i++)
        freq[i] = 0.0;
    for (i = 0; i < clus; i++) {
        feff = eff[i];
        //		if( seq[i][0] == '-' )
        if (seq[i][0] == '-' && sgappat[i] == '-')
            freq[0] += feff;
        for (j = 1; j < len; j++) {
            if (seq[i][j] == '-' && seq[i][j - 1] == '-')
                freq[j] += feff;
        }
        //		if( seq[i][len] == '-' && seq[i][len-1] == '-' ) // xxx wo tsukawanaitoki arienai
        if (egappat[i] == '-' && seq[i][len - 1] == '-')
            freq[len] += feff;
    }
    //	reporterr(       "\ndigapf = \n" );
    //	for( i=0; i<len+1; i++ ) reporterr(       "%5.3f ", freq[i] );
}

void
getdiaminofreq_part(double* freq, int clus, char** seq, double* eff, int len, char* sgappat, char* egappat) {
    int    i, j;
    double feff;
    for (i = 0; i < len + 2; i++)
        freq[i] = 0.0;
    for (i = 0; i < clus; i++) {
        feff = eff[i];
        if (seq[i][0] != '-' && sgappat[i] != '-')
            freq[0] += feff;
        for (j = 1; j < len; j++) {
            if (seq[i][j] != '-' && seq[i][j - 1] != '-')
                freq[j] += feff;
        }
        //		if( 1 && seq[i][len-1] != '-' ) // xxx wo tsukawanaitoki [len-1] nomi
        if (egappat[i] != '-' && seq[i][len - 1] != '-')  // xxx wo tsukawanaitoki [len-1] nomi
            freq[len] += feff;
    }
    //	reporterr(       "\ndiaaf = \n" );
    //	for( i=0; i<len+1; i++ ) reporterr(       "%5.3f ", freq[i] );
}

void
getgapfreq_zure_part(double* freq, int clus, char** seq, double* eff, int len, char* sgap) {
    int    i, j;
    double feff;
    for (i = 0; i < len + 2; i++)
        freq[i] = 0.0;
    for (i = 0; i < clus; i++) {
        feff = eff[i];
        if (sgap[i] == '-')
            freq[0] += feff;
        for (j = 0; j < len; j++) {
            if (seq[i][j] == '-')
                freq[j + 1] += feff;
        }
        //		if( egap[i] == '-' )
        //			freq[len+1] += feff;
    }
    //	reporterr(       "\ngapf = \n" );
    //	for( i=0; i<len+1; i++ ) reporterr(       "%5.3f ", freq[i] );
}

void
getgapfreq_zure(double* freq, int clus, char** seq, double* eff, int len) {
    int    i, j;
    double feff;
    for (i = 0; i < len + 1; i++)
        freq[i] = 0.0;
    for (i = 0; i < clus; i++) {
        feff = eff[i];
        for (j = 0; j < len; j++) {
            if (seq[i][j] == '-')
                freq[j + 1] += feff;
        }
    }
    freq[len + 1] = 0.0;
    //	reporterr(       "\ngapf = \n" );
    //	for( i=0; i<len+1; i++ ) reporterr(       "%5.3f ", freq[i] );
}

void
getgapfreq(double* freq, int clus, char** seq, double* eff, int len) {
    int    i, j;
    double feff;
    for (i = 0; i < len + 1; i++)
        freq[i] = 0.0;
    for (i = 0; i < clus; i++) {
        feff = eff[i];
        for (j = 0; j < len; j++) {
            if (seq[i][j] == '-')
                freq[j] += feff;
        }
    }
    freq[len] = 0.0;
    //	reporterr(       "\ngapf = \n" );
    //	for( i=0; i<len+1; i++ ) reporterr(       "%5.3f ", freq[i] );
}

void
st_getGapPattern(Gappat** pat, int clus, char** seq, double* eff, int len) {
    int      i, j, k, gb, gc;
    int      known;
    double   feff;
    Gappat** fpt;
    char*    spt;
    int      gaplen;

    fpt = pat;
    i = len + 1;
    while (i--) {
        if (*fpt)
            free(*fpt);
        *fpt++ = NULL;
    }

    for (j = 0; j < clus; j++) {
        //		reporterr(       "seq[%d] = %s\n", j, seq[j] );
        feff = (double)eff[j];

        fpt = pat;
        *fpt = NULL;  // Falign.c kara yobareru tokiha chigau.
        spt = seq[j];
        gc = 0;
        gaplen = 0;

        for (i = 0; i < len + 1; i++)
        //		while( i-- )
        {
            //			reporterr(       "i=%d, gaplen = %d\n", i, gaplen );
            gb = gc;
            gc = (i != len && *spt++ == '-');
            if (gc)
                gaplen++;
            else {
                if (gb && gaplen) {
                    k = 1;
                    known = 0;
                    if (*fpt)
                        for (; (*fpt)[k].len != -1; k++) {
                            if ((*fpt)[k].len == gaplen) {
                                //							reporterr(       "known\n" );
                                known = 1;
                                break;
                            }
                        }

                    if (known == 0) {
                        *fpt = (Gappat*)realloc(*fpt, (k + 3) * sizeof(Gappat));  // mae1 (total), ato2 (len0), term
                        if (!*fpt) {
                            reporterr("Cannot allocate gappattern!'n");
                            reporterr("Use an approximate method, with the --mafft5 option.\n");
                            exit(1);
                        }
                        (*fpt)[k].freq = 0.0;
                        (*fpt)[k].len = gaplen;
                        (*fpt)[k + 1].len = -1;
                        (*fpt)[k + 1].freq = 0.0;  // iranai
                        //						reporterr(       "gaplen=%d, Unknown, %f\n", gaplen, (*fpt)[k].freq );
                    }

                    //					reporterr(       "adding pos %d, len=%d, k=%d, freq=%f->", i, gaplen, k, (*fpt)[k].freq );
                    (*fpt)[k].freq += feff;
                    //					reporterr(       "%f\n", (*fpt)[k].freq );
                    gaplen = 0;
                }
            }
            fpt++;
        }
    }
#if 1
    for (j = 0; j < len + 1; j++) {
        if (pat[j]) {
            //			reporterr(       "j=%d\n", j );
            //			for( i=1; pat[j][i].len!=-1; i++ )
            //				reporterr(       "pos=%d, i=%d, len=%d, freq=%f\n", j, i, pat[j][i].len, pat[j][i].freq );

            pat[j][0].len = 0;  // iminashi
            pat[j][0].freq = 0.0;
            for (i = 1; pat[j][i].len != -1; i++) {
                pat[j][0].freq += pat[j][i].freq;
                //				reporterr(       "totaling, i=%d, result = %f\n", i, pat[j][0].freq );
            }
            //			reporterr(       "totaled, result = %f\n", pat[j][0].freq );

            pat[j][i].freq = 1.0 - pat[j][0].freq;
            pat[j][i].len = 0;  // imiari
            pat[j][i + 1].len = -1;
        } else {
            pat[j] = (Gappat*)calloc(3, sizeof(Gappat));
            pat[j][0].freq = 0.0;
            pat[j][0].len = 0;  // iminashi

            pat[j][1].freq = 1.0 - pat[j][0].freq;
            pat[j][1].len = 0;  // imiari
            pat[j][2].len = -1;
        }
    }
#endif
}

static int
minimum(int i1, int i2) {
    return MIN(i1, i2);
}

static void
commongappickpairfast(char* r1, char* r2, const char* i1, const char* i2, int* skip1, int* skip2) {
    int skip, skipped1, skipped2;
    skipped1 = skipped2 = 0;
    while (1) {
        skip = minimum(*skip1 - skipped1, *skip2 - skipped2);
        i1 += skip;
        i2 += skip;
        skipped1 += skip;
        skipped2 += skip;
        //		fprintf( stderr, "i1 pos =%d, nlenmax=%d\n", (int)(i1- i1bk), nlenmax );
        if (!*i1)
            break;
        //		reporterr( "%d, %c-%c\n", i1-i1bk, *i1, *i2 );
        //		if( *i1 == '-' && *i2 == '-' )  // iranai?
        //		{
        //			reporterr( "Error in commongappickpairfast" );
        //			exit( 1 );
        //			i1++;
        //			i2++;
        //		}
        if (*i1 != '-') {
            skipped1 = 0;
            skip1++;
        } else
            skipped1++;

        if (*i2 != '-') {
            skipped2 = 0;
            skip2++;
        } else
            skipped2++;

        *r1++ = *i1++;
        *r2++ = *i2++;
    }
    *r1 = 0;
    *r2 = 0;
}

static void
commongappickpair(char* r1, char* r2, const char* i1, const char* i2) {
    while (*i1) {
        if (*i1 == '-' && *i2 == '-') {
            i1++;
            i2++;
        } else {
            *r1++ = *i1++;
            *r2++ = *i2++;
        }
    }
    *r1 = 0;
    *r2 = 0;
}

double
naivepairscorefast(Context* ctx, const char* seq1, const char* seq2, int* skip1, int* skip2, int penal) {
    double vali;
    int    len = strlen(seq1);
    char * s1, *s2;
    char * p1, *p2;

    s1 = calloc(len + 1, sizeof(char));
    s2 = calloc(len + 1, sizeof(char));
    {
        vali = 0.0;
        commongappickpairfast(s1, s2, seq1, seq2, skip1, skip2);

        p1 = s1;
        p2 = s2;
        while (*p1) {
            if (*p1 == '-') {
                //				reporterr(       "Penal! %c-%c in %d-%d, %f\n", *(p1-1), *(p2-1), i, j, feff );
                vali += (double)penal;
                //				while( *p1 == '-' || *p2 == '-' )
                while (*p1 == '-')  // SP
                {
                    p1++;
                    p2++;
                }
                continue;
            }
            if (*p2 == '-') {
                //				reporterr(       "Penal! %c-%c in %d-%d, %f\n", *(p1-1), *(p2-1), i, j, feff );
                vali += (double)penal;
                //				while( *p2 == '-' || *p1 == '-' )
                while (*p2 == '-')  // SP
                {
                    p1++;
                    p2++;
                }
                continue;
            }
            //			reporterr(       "adding %c-%c, %d\n", *p1, *p2, amino_dis[*p1][*p2] );
            vali += (double)ctx->amino_dis[(unsigned char)*p1++][(unsigned char)*p2++];
        }
    }
    free(s1);
    free(s2);
    //	reporterr(       "###vali = %d\n", vali );
    return (vali);
}

double
naivepairscore11_dynmtx(Context* ctx, double** mtx, char* seq1, char* seq2, int penal) {
    double vali;
    int    len = strlen(seq1);
    char * s1, *s2, *p1, *p2;
    int    c1, c2;

    s1 = calloc(len + 1, sizeof(char));
    s2 = calloc(len + 1, sizeof(char));
    {
        vali = 0.0;
        commongappickpair(s1, s2, seq1, seq2);
        //		reporterr(       "###i1 = %s\n", s1 );
        //		reporterr(       "###i2 = %s\n", s2 );
        //		reporterr(       "###penal = %d\n", penal );

        p1 = s1;
        p2 = s2;
        while (*p1) {
            if (*p1 == '-') {
                //				reporterr(       "Penal! %c-%c in %d-%d, %f\n", *(p1-1), *(p2-1), i, j, feff );
                vali += (double)penal;
                //				while( *p1 == '-' || *p2 == '-' )
                while (*p1 == '-')  // SP
                {
                    p1++;
                    p2++;
                }
                continue;
            }
            if (*p2 == '-') {
                //				reporterr(       "Penal! %c-%c in %d-%d, %f\n", *(p1-1), *(p2-1), i, j, feff );
                vali += (double)penal;
                //				while( *p2 == '-' || *p1 == '-' )
                while (*p2 == '-')  // SP
                {
                    p1++;
                    p2++;
                }
                continue;
            }
            c1 = ctx->amino_n[(unsigned char)*p1++];
            c2 = ctx->amino_n[(unsigned char)*p2++];
            vali += (double)mtx[c1][c2];
        }
    }
    free(s1);
    free(s2);
    //	reporterr(       "###vali = %d\n", vali );
    return (vali);
}

double
naivepairscore11(Context* ctx, const char* seq1, const char* seq2, int penal) {
    double vali;
    int    len = strlen(seq1);
    char * s1, *s2, *p1, *p2;

    s1 = calloc(len + 1, sizeof(char));
    s2 = calloc(len + 1, sizeof(char));
    {
        vali = 0.0;
        commongappickpair(s1, s2, seq1, seq2);
        //		reporterr(       "###i1 = %s\n", s1 );
        //		reporterr(       "###i2 = %s\n", s2 );
        //		reporterr(       "###penal = %d\n", penal );

        p1 = s1;
        p2 = s2;
        while (*p1) {
            if (*p1 == '-') {
                //				reporterr(       "Penal! %c-%c in %d-%d, %f\n", *(p1-1), *(p2-1), i, j, feff );
                vali += (double)penal;
                //				while( *p1 == '-' || *p2 == '-' )
                while (*p1 == '-')  // SP
                {
                    p1++;
                    p2++;
                }
                continue;
            }
            if (*p2 == '-') {
                //				reporterr(       "Penal! %c-%c in %d-%d, %f\n", *(p1-1), *(p2-1), i, j, feff );
                vali += (double)penal;
                //				while( *p2 == '-' || *p1 == '-' )
                while (*p2 == '-')  // SP
                {
                    p1++;
                    p2++;
                }
                continue;
            }
            vali += (double)ctx->amino_dis[(unsigned char)*p1++][(unsigned char)*p2++];
        }
    }
    free(s1);
    free(s2);
    //	reporterr(       "###vali = %d\n", vali );
    return (vali);
}

double
plainscore(Context* ctx, int nseq, char** s) {
    int    i, j, ilim;
    double v = 0.0;

    ilim = nseq - 1;
    for (i = 0; i < ilim; i++)
        for (j = i + 1; j < nseq; j++) {
            v += (double)naivepairscore11(ctx, s[i], s[j], ctx->penalty);
        }

    reporterr("penalty = %d\n", ctx->penalty);

    return (v);
}

int
samemember(int* mem, int* cand) {
    int i, j;
    int nm, nc;

    nm = 0;
    for (i = 0; mem[i] > -1; i++)
        nm++;
    nc = 0;
    for (i = 0; cand[i] > -1; i++)
        nc++;

    if (nm != nc)
        return (0);

    for (i = 0; mem[i] > -1; i++) {
        for (j = 0; cand[j] > -1; j++)
            if (mem[i] == cand[j])
                break;
        if (cand[j] == -1)
            return (0);
    }

    if (mem[i] == -1) {
        return (1);
    } else {
        return (0);
    }
}

int
samemembern(int* mem, int* cand, int nc) {
    int i, j;
    int nm;

    nm = 0;
    for (i = 0; mem[i] > -1; i++) {
        nm++;
        if (nm > nc)
            return (0);
    }

    if (nm != nc)
        return (0);

    for (i = 0; mem[i] > -1; i++) {
        for (j = 0; j < nc; j++)
            if (mem[i] == cand[j])
                break;
        if (j == nc)
            return (0);
    }

    if (mem[i] == -1) {
#if 0
		reporterr(       "mem = " );
		for( i=0; mem[i]>-1; i++ )	reporterr(       "%d ", mem[i] );
		reporterr(       "\n" );
	
		reporterr(       "cand = " );
		for( i=0; cand[i]>-1; i++ )	reporterr(       "%d ", cand[i] );
		reporterr(       "\n" );
#endif
        return (1);
    } else {
        return (0);
    }
}

int
includemember(int* mem, int* cand)  // mem in cand
{
    int i, j;

#if 0
	reporterr(       "mem = " );
	for( i=0; mem[i]>-1; i++ )	reporterr(       "%d ", mem[i] );
	reporterr(       "\n" );

	reporterr(       "cand = " );
	for( i=0; cand[i]>-1; i++ )	reporterr(       "%d ", cand[i] );
	reporterr(       "\n" );
#endif

    for (i = 0; mem[i] > -1; i++) {
        for (j = 0; cand[j] > -1; j++)
            if (mem[i] == cand[j])
                break;
        if (cand[j] == -1)
            return (0);
    }
    //	reporterr(       "INCLUDED! mem[0]=%d\n", mem[0] );
    return (1);
}

int
overlapmember(int* mem1, int* mem2) {
    int i, j;

    for (i = 0; mem1[i] > -1; i++)
        for (j = 0; mem2[j] > -1; j++)
            if (mem1[i] == mem2[j])
                return (1);
    return (0);
}

void
gapcount(double* freq, char** seq, int nseq, double* eff, int lgth) {
    int    i, j;
    double fr;

    //	for( i=0; i<lgth; i++ ) freq[i] = 0.0;
    //	return;

    for (i = 0; i < lgth; i++) {
        fr = 0.0;
        for (j = 0; j < nseq; j++) {
            if (seq[j][i] == '-')
                fr += eff[j];
        }
        freq[i] = fr;
        //		reporterr(       "freq[%d] = %f\n", i, freq[i] );
    }
    //	reporterr(       "\n" );
    return;
}

void
gapcountadd(double* freq, char** seq, int nseq, double* eff, int lgth) {
    int    i;
    int    j = nseq - 1;
    double newfr = eff[j];
    double orifr = 1.0 - newfr;

    //	for( i=0; i<lgth; i++ ) freq[i] = 0.0;
    //	return;
    //	for( i=0; i<nseq; i++ )
    //		reporterr( "%s\n", seq[i] );

    for (i = 0; i < lgth; i++) {
        //		reporterr(       "freq[%d] = %f", i, freq[i] );
        freq[i] = 1.0 - freq[i];  // modosu
        freq[i] *= orifr;

        if (seq[j][i] == '-')
            freq[i] += newfr;
        //		reporterr(       "->         %f\n", i, freq[i] );
    }
    //	reporterr(       "\n" );
    return;
}
void
gapcountf(double* freq, char** seq, int nseq, double* eff, int lgth) {
    int    i, j;
    double fr;

    //	for( i=0; i<lgth; i++ ) freq[i] = 0.0;
    //	return;

    for (i = 0; i < lgth; i++) {
        fr = 0.0;
        for (j = 0; j < nseq; j++) {
            if (seq[j][i] == '-')
                fr += eff[j];
        }
        freq[i] = fr;
        //		reporterr(       "in gapcountf, freq[%d] = %f\n", i, freq[i] );
    }
    //	reporterr(       "\n" );
    return;
}

void
outgapcount(double* freq, int nseq, char* gappat, double* eff) {
    int    j;
    double fr;

    fr = 0.0;
    for (j = 0; j < nseq; j++) {
        if (gappat[j] == '-')
            fr += eff[j];
    }
    *freq = fr;
    return;
}

static double
dist2offset(double dist) {
    double val = dist * 0.5;
    if (val > 0.0)
        val = 0.0;
    return val;
}

void
makedynamicmtx(Context* ctx, double** out, double** in, double offset) {
    int    i, j, ii, jj;
    double av;

    offset = dist2offset(offset * 2.0);  // offset 0..1 -> 0..2

    //	if( offset > 0.0 ) offset = 0.0;
    //	reporterr(       "dynamic offset = %f\n", offset );

    for (i = 0; i < ctx->nalphabets; i++)
        for (j = 0; j < ctx->nalphabets; j++) {
            out[i][j] = in[i][j];
        }
    if (offset == 0.0)
        return;

    for (i = 0; i < ctx->nalphabets; i++) {
        ii = (int)ctx->amino[i];
        if (ii == '-')
            continue;  // text no toki arieru
        for (j = 0; j < ctx->nalphabets; j++) {
            jj = (int)ctx->amino[j];
            if (jj == '-')
                continue;  // text no toki arieru
            out[i][j] = in[i][j] + offset * 600;
            //			reporterr(       "%c-%c: %f\n", ii, jj, out[i][j] );
        }
    }

    //	reporterr(       "offset = %f\n", offset );
    //	reporterr(       "out[W][W] = %f\n", out[amino_n['W']][amino_n['W']] );
    //	reporterr(       "out[A][A] = %f\n", out[amino_n['A']][amino_n['A']] );

    return;

    // Taikaku youso no heikin ga 600 ni naruyouni re-scale.
    // Hitaikaku youso ga ookiku narisugi.

    av = 0.0;
    for (i = 0; i < ctx->nalphabets; i++) {
        if (ii == '-')
            continue;  // text no toki arieru
        av += out[i][i];
    }
    av /= (double)ctx->nalphabets;

    for (i = 0; i < ctx->nalphabets; i++) {
        if (ctx->amino[i] == '-')
            continue;  // text no toki arieru
        for (j = 0; j < ctx->nalphabets; j++) {
            if (ctx->amino[j] == '-')
                continue;  // text no toki arieru
            out[i][j] = out[i][j] * 600 / av;
            reporterr("%c-%c: %f\n", ctx->amino[i], ctx->amino[j], out[i][j]);
        }
    }
}

void
makeskiptable(int n, int** skip, char** seq) {
    char* nogapseq;
    int   nogaplen, alnlen;
    int   i, j, posinseq;

    nogapseq = calloc(strlen(seq[0]) + 1, sizeof(char));
    for (i = 0; i < n; i++) {
        copyWithNoGaps(nogapseq, seq[i]);
        nogaplen = strlen(nogapseq);
        alnlen = strlen(seq[i]);
        skip[i] = calloc(nogaplen + 1, sizeof(int));

        //		reporterr( "%s\n", nogapseq );

        posinseq = 0;
        for (j = 0; j < alnlen; j++) {
            if (seq[i][j] == '-') {
                skip[i][posinseq]++;
            } else {
                posinseq++;
            }
        }
        //		for( j=0; j<nogaplen+1; j++ )
        //			reporterr( "%d ", skip[i][j] );
        //		reporterr( "\n" );
        //		exit( 1 );
    }
    free(nogapseq);
}

int
generatesubalignmentstable(int nseq, int*** tablept, int* nsubpt, int* maxmempt, int*** topol, double** len, double threshold) {
    int     i, j, rep0, rep1, nmem, mem;
    double  distfromtip0, distfromtip1;
    double* distfromtip;
    reporterr("\n\n\n");

    *maxmempt = 0;
    *nsubpt = 0;

    distfromtip = calloc(nseq, sizeof(double));
    for (i = 0; i < nseq - 1; i++) {
#if 0
		reporterr( "STEP %d\n", i );
		for( j=0; topol[i][0][j]!=-1; j++ )
			reporterr( "%3d ", topol[i][0][j] );
		reporterr( "\n" );
		reporterr( "len=%f\n", len[i][0] );
#endif

        rep0 = topol[i][0][0];
        distfromtip0 = distfromtip[rep0];
        distfromtip[rep0] += len[i][0];
        //		reporterr( "distfromtip[%d] = %f->%f\n", rep0, distfromtip0, distfromtip[rep0] );

#if 0
		for( j=0; topol[i][1][j]!=-1; j++ )
			reporterr( "%3d ", topol[i][1][j] );
		reporterr( "\n" );
		reporterr( "len=%f\n", len[i][1] );
#endif

        rep1 = topol[i][1][0];
        distfromtip1 = distfromtip[rep1];
        distfromtip[rep1] += len[i][1];
        //		reporterr( "distfromtip[%d] = %f->%f\n", rep1, distfromtip1, distfromtip[rep1] );

        if (topol[i][0][1] != -1 && distfromtip0 <= threshold && threshold < distfromtip[rep0]) {
            //			reporterr( "HIT 0!\n" );
            *tablept = realloc(*tablept, sizeof(char*) * (*nsubpt + 2));
            for (j = 0, nmem = 0; (mem = topol[i][0][j]) != -1; j++)
                nmem++;
            //			reporterr( "allocating %d\n", nmem+1 );
            (*tablept)[*nsubpt] = calloc(nmem + 1, sizeof(int));
            (*tablept)[*nsubpt + 1] = NULL;
            intcpy((*tablept)[*nsubpt], topol[i][0]);
            if (*maxmempt < nmem)
                *maxmempt = nmem;
            *nsubpt += 1;
        }

        if (topol[i][1][1] != -1 && distfromtip1 <= threshold && threshold < distfromtip[rep1]) {
            //			reporterr( "HIT 1!\n" );
            *tablept = realloc(*tablept, sizeof(char*) * (*nsubpt + 2));
            for (j = 0, nmem = 0; (mem = topol[i][1][j]) != -1; j++)
                nmem++;
            //			reporterr( "allocating %d\n", nmem+1 );
            (*tablept)[*nsubpt] = calloc(nmem + 1, sizeof(int));
            (*tablept)[*nsubpt + 1] = NULL;
            intcpy((*tablept)[*nsubpt], topol[i][1]);
            if (*maxmempt < nmem)
                *maxmempt = nmem;
            *nsubpt += 1;
        }
    }

    if (distfromtip[0] <= threshold) {
        free(distfromtip);
        return (1);
    }

    free(distfromtip);
    return (0);
}

double
sumofpairsscore(Context* ctx, int nseq, char** seq) {
    double v = 0;
    int    i, j;
    for (i = 1; i < nseq; i++) {
        for (j = 0; j < i; j++) {
            v += naivepairscore11(ctx, seq[i], seq[j], ctx->penalty) / 600;
        }
    }
    //	v /= ( (nseq-1) * nseq ) / 2;
    return (v);
}

static void
movereg(char* seq1, char* seq2, LocalHom* tmpptr, int* start1pt, int* start2pt, int* end1pt, int* end2pt) {
    char* pt;
    int   tmpint;

    pt = seq1;
    tmpint = -1;
    while (*pt != 0) {
        if (*pt++ != '-')
            tmpint++;
        if (tmpint == tmpptr->start1)
            break;
    }
    *start1pt = (int)(pt - seq1) - 1;

    if (tmpptr->start1 == tmpptr->end1)
        *end1pt = *start1pt;
    else {
        while (*pt != 0) {
            //			fprintf( stderr, "tmpint = %d, end1 = %d pos = %d\n", tmpint, tmpptr->end1, pt-seq1[i] );
            if (*pt++ != '-')
                tmpint++;
            if (tmpint == tmpptr->end1)
                break;
        }
        *end1pt = (int)(pt - seq1) - 1;
    }

    pt = seq2;
    tmpint = -1;
    while (*pt != 0) {
        if (*pt++ != '-')
            tmpint++;
        if (tmpint == tmpptr->start2)
            break;
    }
    *start2pt = (int)(pt - seq2) - 1;
    if (tmpptr->start2 == tmpptr->end2)
        *end2pt = *start2pt;
    else {
        while (*pt != 0) {
            if (*pt++ != '-')
                tmpint++;
            if (tmpint == tmpptr->end2)
                break;
        }
        *end2pt = (int)(pt - seq2) - 1;
    }
}

static void
movereg_swap(char* seq1, char* seq2, LocalHom* tmpptr, int* start1pt, int* start2pt, int* end1pt, int* end2pt) {
    char* pt;
    int   tmpint;

    pt = seq1;
    tmpint = -1;
    while (*pt != 0) {
        if (*pt++ != '-')
            tmpint++;
        if (tmpint == tmpptr->start2)
            break;
    }
    *start1pt = (int)(pt - seq1) - 1;

    if (tmpptr->start2 == tmpptr->end2)
        *end1pt = *start1pt;
    else {
        while (*pt != 0) {
            //			fprintf( stderr, "tmpint = %d, end1 = %d pos = %d\n", tmpint, tmpptr->end1, pt-seq1[i] );
            if (*pt++ != '-')
                tmpint++;
            if (tmpint == tmpptr->end2)
                break;
        }
        *end1pt = (int)(pt - seq1) - 1;
    }

    pt = seq2;
    tmpint = -1;
    while (*pt != 0) {
        if (*pt++ != '-')
            tmpint++;
        if (tmpint == tmpptr->start1)
            break;
    }
    *start2pt = (int)(pt - seq2) - 1;
    if (tmpptr->start1 == tmpptr->end1)
        *end2pt = *start2pt;
    else {
        while (*pt != 0) {
            if (*pt++ != '-')
                tmpint++;
            if (tmpint == tmpptr->end1)
                break;
        }
        *end2pt = (int)(pt - seq2) - 1;
    }
}

void
fillimp(aln_Opts opts, double** impmtx, int clus1, int clus2, int lgth1, int lgth2, char** seq1, char** seq2, double* eff1, double* eff2, double* eff1_kozo, double* eff2_kozo, LocalHom*** localhom, char* swaplist, int* orinum1, int* orinum2) {
    int       i, j, k1, k2, start1, start2, end1, end2;
    double    effij, effijx, effij_kozo;
    char *    pt1, *pt2;
    LocalHom* tmpptr;
    void (*movefunc)(char*, char*, LocalHom*, int*, int*, int*, int*);

    for (i = 0; i < lgth1; i++)
        for (j = 0; j < lgth2; j++)
            impmtx[i][j] = 0.0;
    effijx = 1.0 * opts.fastathreshold;
    for (i = 0; i < clus1; i++) {
        if (swaplist && swaplist[i])
            movefunc = movereg_swap;
        else
            movefunc = movereg;
        for (j = 0; j < clus2; j++) {
            if (swaplist == NULL && orinum1 && orinum2)  // muda.
            {
                if (orinum1[i] > orinum2[j])
                    movefunc = movereg_swap;
                else
                    movefunc = movereg;
            }

            //			effij = eff1[i] * eff2[j] * effijx;
            effij = eff1[i] * eff2[j] * effijx;
            effij_kozo = eff1_kozo[i] * eff2_kozo[j] * effijx;
            tmpptr = localhom[i][j];
            while (tmpptr) {
                //				fprintf( stderr, "start1 = %d\n", tmpptr->start1 );
                //				fprintf( stderr, "end1   = %d\n", tmpptr->end1   );
                //				fprintf( stderr, "i = %d, seq1 = \n%s\n", i, seq1[i] );
                //				fprintf( stderr, "j = %d, seq2 = \n%s\n", j, seq2[j] );

                movefunc(seq1[i], seq2[j], tmpptr, &start1, &start2, &end1, &end2);

                //				fprintf( stderr, "start1 = %d (%c), end1 = %d (%c), start2 = %d (%c), end2 = %d (%c)\n", start1, seq1[i][start1], end1, seq1[i][end1], start2, seq2[j][start2], end2, seq2[j][end2] );
                //				fprintf( stderr, "step 0\n" );
                if (end1 - start1 != end2 - start2) {
                    //					fprintf( stderr, "CHUUI!!, start1 = %d, end1 = %d, start2 = %d, end2 = %d\n", start1, end1, start2, end2 );
                }

                k1 = start1;
                k2 = start2;
                pt1 = seq1[i] + k1;
                pt2 = seq2[j] + k2;
                while (*pt1 && *pt2) {
                    if (*pt1 != '-' && *pt2 != '-') {
                        // �Ťߤ���Ťˤ����ʤ��褦�����դ��Ʋ�������
                        //						impmtx[k1][k2] += tmpptr->wimportance * fastathreshold;
                        //						impmtx[k1][k2] += tmpptr->importance * effij;
                        //						impmtx[k1][k2] += tmpptr->fimportance * effij;
                        if (tmpptr->korh == 'k')
                            impmtx[k1][k2] += tmpptr->importance * effij_kozo;
                        else
                            impmtx[k1][k2] += tmpptr->importance * effij;
                        //						fprintf( stderr, "k1=%d, k2=%d, impalloclen=%d\n", k1, k2, impalloclen );
                        //						fprintf( stderr, "mark, %d (%c) - %d (%c) \n", k1, *pt1, k2, *pt2 );
                        k1++;
                        k2++;
                        pt1++;
                        pt2++;
                    } else if (*pt1 != '-' && *pt2 == '-') {
                        //						fprintf( stderr, "skip, %d (%c) - %d (%c) \n", k1, *pt1, k2, *pt2 );
                        k2++;
                        pt2++;
                    } else if (*pt1 == '-' && *pt2 != '-') {
                        //						fprintf( stderr, "skip, %d (%c) - %d (%c) \n", k1, *pt1, k2, *pt2 );
                        k1++;
                        pt1++;
                    } else if (*pt1 == '-' && *pt2 == '-') {
                        //						fprintf( stderr, "skip, %d (%c) - %d (%c) \n", k1, *pt1, k2, *pt2 );
                        k1++;
                        pt1++;
                        k2++;
                        pt2++;
                    }
                    if (k1 > end1 || k2 > end2)
                        break;
                }
                tmpptr = tmpptr->next;
            }
        }
    }
#if 0
	printf( "orinum1=%d, orinum2=%d\n", *orinum1, *orinum2 );
	if( *orinum1 == 0 )
	{
		fprintf( stdout, "impmtx = \n" );
		for( k2=0; k2<lgth2; k2++ )
			fprintf( stdout, "%6.3f ", (double)k2 );
		fprintf( stdout, "\n" );
		for( k1=0; k1<lgth1; k1++ )
		{
			fprintf( stdout, "%d", k1 );
			for( k2=0; k2<lgth2; k2++ )
				fprintf( stdout, "%2.1f ", impmtx[k1][k2] );
			fprintf( stdout, "\n" );
		}
		exit( 1 );
	}
#endif
}

static void
readlocalhomtable2_single_bin_noseek(FILE* fp, LocalHom* localhomtable)  // pos ha tsukawanai
{
    double    opt;
    int       size;
    LocalHom* tmpptr1;
    int       st1, st2, len;
    char      c;

    void* m;
    int*  p;

    fread(&size, sizeof(int), 1, fp);
    fread(&opt, sizeof(double), 1, fp);

    m = malloc(sizeof(int) * 3 * size);
    fread(m, size * sizeof(int), 3, fp);
    p = (int*)m;
    while (size--) {
        st1 = *(p++);
        st2 = *(p++);
        len = *(p++);

        if (localhomtable->nokori++ > 0) {
            tmpptr1 = localhomtable->last;
            tmpptr1->next = (LocalHom*)calloc(1, sizeof(LocalHom));
            tmpptr1 = tmpptr1->next;
            tmpptr1->extended = -1;
            tmpptr1->next = NULL;
            localhomtable->last = tmpptr1;
        } else {
            tmpptr1 = localhomtable;
        }

        tmpptr1->start1 = st1;
        tmpptr1->start2 = st2;
        tmpptr1->end1 = st1 + len;
        tmpptr1->end2 = st2 + len;
        //		tmpptr1->opt = ( opt / overlapaa + 0.00 ) / 5.8  * 600;
        //		tmpptr1->opt = opt;
        tmpptr1->opt = ((double)opt + 0.00) / 5.8 * 600;
        tmpptr1->importance = ((double)opt + 0.00) / 5.8 * 600;  // C0 to itchi shinai
        tmpptr1->overlapaa = len;  // tsukau toki ha chuui
        tmpptr1->korh = 'h';

        //		fprintf( stderr, " %f %d-%d %d-%d \n",  tmpptr1->opt, tmpptr1->start1, tmpptr1->end1, tmpptr1->start2, tmpptr1->end2 );
    }
    free(m);
    fread(&c, sizeof(char), 1, fp);
    if (c != '\n') {
        reporterr("\n\nError in binary hat3  \n");
        exit(1);
    }
}

static int
readlocalhomfromfile_autofid(LocalHom* lhpt, FILE* fp, int o1, int o2)  // for hat3node
{
    int swap;

    lhpt->start1 = -1;
    lhpt->end1 = -1;
    lhpt->start2 = -1;
    lhpt->end2 = -1;
    lhpt->overlapaa = -1.0;
    lhpt->opt = -1.0;
    lhpt->importance = -1.0;
    lhpt->next = NULL;
    lhpt->nokori = 0;
    lhpt->extended = -1;
    lhpt->last = lhpt;
    lhpt->korh = 'h';

    if (o2 > o1) {
        swap = 0;
    } else {
        swap = 1;
    }

    if (fp) {
        readlocalhomtable2_single_bin_noseek(fp, lhpt);
    }
    return (swap);
}

static int
whichpair(int* ipt, int* jpt, FILE* fp) {
    if (fread(ipt, sizeof(int), 1, fp) < 1)
        return (1);
    if (fread(jpt, sizeof(int), 1, fp) < 1)
        return (1);  // <1 ha nai
    return (0);
}

typedef struct _readloopthread_arg {
    aln_Opts opts;
    int                 nodeid;
    int                 nfiles;
    double**            impmtx;
    char**              seq1;
    char**              seq2;
    int*                orinum1;
    int*                orinum2;
    double*             eff1;
    double*             eff2;
    unsigned long long* ndone;
    int*                subidpt;
} readloopthread_arg_t;

static void*
readloopthread(void* arg) {
    readloopthread_arg_t* targ = (readloopthread_arg_t*)arg;
    aln_Opts opts = targ->opts;
    int                   nodeid = targ->nodeid;
    //	int thread_no = targ->thread_no;
    double**            impmtx = targ->impmtx;
    char**              seq1 = targ->seq1;
    char**              seq2 = targ->seq2;
    int*                orinum1 = targ->orinum1;
    int*                orinum2 = targ->orinum2;
    double*             eff1 = targ->eff1;
    double*             eff2 = targ->eff2;
    unsigned long long* ndone = targ->ndone;
    int*                subidpt = targ->subidpt;
    int                 nfiles = targ->nfiles;
    int                 subid = -1;
    int       i, j, k1, k2, start1, start2, end1, end2;
    double    effij, effijx;
    char *    pt1, *pt2;
    LocalHom* tmpptr;
    FILE*     fp = NULL;
    char*     fn;
    LocalHom  lhsingle;
    int       res;
    void (*movefunc)(char*, char*, LocalHom*, int*, int*, int*, int*);
    initlocalhom1(&lhsingle);
    effijx = 1.0 * opts.fastathreshold;

    while (1) {
        if (subid == -1 || whichpair(&i, &j, fp)) {
            while (1) {
                if (fp)
                    fclose(fp);

                subid = (*subidpt)++;

                if (subid >= nfiles) {
                    return (NULL);
                }

                fn = calloc(100, sizeof(char));
                sprintf(fn, "hat3dir/%d-/hat3node-%d-%d", (int)(nodeid / HAT3NODEBLOCK) * HAT3NODEBLOCK, nodeid, subid);
                fp = fopen(fn, "rb");
                if (fp == NULL) {
                    reporterr("Cannot open %s\n", fn);
                    exit(1);
                }
                free(fn);
                setvbuf(fp, NULL, _IOFBF, MYBUFSIZE);
                if (!whichpair(&i, &j, fp))
                    break;
            }
        }
        (*ndone)++;

        {
            //			effij = eff1[i] * eff2[j] * effijx;
            effij = eff1[i] * eff2[j] * effijx;
            //			effij_kozo = eff1_kozo[i] * eff2_kozo[j] * effijx;

            res = readlocalhomfromfile_autofid(&lhsingle, fp, orinum1[i], orinum2[j]);
            if (res == -1)
                tmpptr = NULL;

            else if (res == 1) {
                movefunc = movereg_swap;  // h3i ga arutoki swaplist ha mushi
                tmpptr = &lhsingle;
            } else {
                movefunc = movereg;  // h3i ga arutoki swaplist ha mushi
                tmpptr = &lhsingle;
            }

            while (tmpptr) {
                //				fprintf( stderr, "start1 = %d\n", tmpptr->start1 );
                //				fprintf( stderr, "end1   = %d\n", tmpptr->end1   );
                //				fprintf( stderr, "i = %d, seq1 = \n%s\n", i, seq1[i] );
                //				fprintf( stderr, "j = %d, seq2 = \n%s\n", j, seq2[j] );

                movefunc(seq1[i], seq2[j], tmpptr, &start1, &start2, &end1, &end2);

                //				fprintf( stderr, "start1 = %d (%c), end1 = %d (%c), start2 = %d (%c), end2 = %d (%c)\n", start1, seq1[i][start1], end1, seq1[i][end1], start2, seq2[j][start2], end2, seq2[j][end2] );
                //				fprintf( stderr, "step 0\n" );
                //				if( end1 - start1 != end2 - start2 )
                //					fprintf( stderr, "CHUUI!!, start1 = %d, end1 = %d, start2 = %d, end2 = %d\n", start1, end1, start2, end2 );

                k1 = start1;
                k2 = start2;
                pt1 = seq1[i] + k1;
                pt2 = seq2[j] + k2;
                while (*pt1 && *pt2) {
                    if (*pt1 != '-' && *pt2 != '-') {
                        impmtx[k1][k2] += tmpptr->importance * effij;
                        k1++;
                        k2++;
                        pt1++;
                        pt2++;
                    } else if (*pt1 != '-' && *pt2 == '-') {
                        //						fprintf( stderr, "skip, %d (%c) - %d (%c) \n", k1, *pt1, k2, *pt2 );
                        k2++;
                        pt2++;
                    } else if (*pt1 == '-' && *pt2 != '-') {
                        //						fprintf( stderr, "skip, %d (%c) - %d (%c) \n", k1, *pt1, k2, *pt2 );
                        k1++;
                        pt1++;
                    } else if (*pt1 == '-' && *pt2 == '-') {
                        //						fprintf( stderr, "skip, %d (%c) - %d (%c) \n", k1, *pt1, k2, *pt2 );
                        k1++;
                        pt1++;
                        k2++;
                        pt2++;
                    }
                    if (k1 > end1 || k2 > end2)
                        break;
                }
                tmpptr = tmpptr->next;
            }
            freelocalhom1(&lhsingle);
        }
    }
}

void
fillimp_file(aln_Opts opts, Context* ctx, double** impmtx, int clus1, int clus2, int lgth1, int lgth2, char** seq1, char** seq2, double* eff1, double* eff2, double* eff1_kozo, double* eff2_kozo, LocalHom*** localhom, int* orinum1, int* orinum2, int* uselh, int* seedinlh1, int* seedinlh2, int nodeid, int nfiles) {
    int                i, j, k1, k2, start1, start2, end1, end2, m0, m1, m2;
    double             effijx, effij_kozo;
    char *             pt1, *pt2;
    LocalHom*          tmpptr;
    unsigned long long npairs;
    //	LocalHom lhsingle;
    //	FILE *fp = NULL;
    //	char *fn;
    //	int subid, res;
    void (*movefunc)(char*, char*, LocalHom*, int*, int*, int*, int*);
    unsigned long long  ndone;
    int                 subid;

#if 0
	fprintf( stderr, "eff1 in _init_strict = \n" );
	for( i=0; i<clus1; i++ )
		fprintf( stderr, "eff1[] = %f\n", eff1[i] );
	for( i=0; i<clus2; i++ )
		fprintf( stderr, "eff2[] = %f\n", eff2[i] );
#endif

    for (i = 0; i < lgth1; i++)
        for (j = 0; j < lgth2; j++)
            impmtx[i][j] = 0.0;
    effijx = 1.0 * opts.fastathreshold;

    if (ctx->nadd) {
        npairs = 0;
        for (i = 0; i < clus1; i++)
            for (j = 0; j < clus2; j++) {
                m1 = orinum1[i];
                m2 = orinum2[j];
                if (m1 > m2) {
                    m0 = m1;
                    m1 = m2;
                    m2 = m0;
                }
                if (m2 >= ctx->njob - ctx->nadd && (uselh == NULL || uselh[m1] || uselh[m2]))  // saikentou
                {
                    //				reporterr( "%d x %d\n", m1, m2 );
                    npairs++;
                }
            }
    } else if (uselh) {
        //		npairs = (unsigned long long)clus1 * clus2;
        npairs = 0;
        for (i = 0; i < clus1; i++)
            for (j = 0; j < clus2; j++)  // sukoshi muda. lh == NULL nara zenbu tsukau toka.
            {
                m1 = orinum1[i];
                m2 = orinum2[j];
                if (uselh[m1] || uselh[m2]) {
                    //				reporterr( "%d x %d\n", m1, m2 );
                    npairs++;
                }
            }
    } else {
        npairs = (unsigned long long)clus1 * clus2;
    }

    if (localhom)  // seed yurai
    {
        for (i = 0; i < clus1; i++) {
            if (seedinlh1[i] == -1)
                continue;
            for (j = 0; j < clus2; j++) {
                if (seedinlh2[j] == -1)
                    continue;

                //				reporterr( "adding %d-%d\n", orinum1[i], orinum2[j] );

                effij_kozo = eff1_kozo[i] * eff2_kozo[j] * effijx;
                tmpptr = localhom[seedinlh1[i]][seedinlh2[j]];

                if (orinum1[i] > orinum2[j])
                    movefunc = movereg_swap;
                else
                    movefunc = movereg;

                while (tmpptr) {
                    movefunc(seq1[i], seq2[j], tmpptr, &start1, &start2, &end1, &end2);

                    k1 = start1;
                    k2 = start2;
                    pt1 = seq1[i] + k1;
                    pt2 = seq2[j] + k2;
                    while (*pt1 && *pt2) {
                        if (*pt1 != '-' && *pt2 != '-') {
                            if (tmpptr->korh == 'k')
                                impmtx[k1][k2] += tmpptr->importance * effij_kozo;
                            else  // naihazu
                            {
                                reporterr("okashii\n");
                                exit(1);
                            }
                            //							fprintf( stderr, "k1=%d, k2=%d, impalloclen=%d\n", k1, k2, impalloclen );
                            //							fprintf( stderr, "mark, %d (%c) - %d (%c) \n", k1, *pt1, k2, *pt2 );
                            k1++;
                            k2++;
                            pt1++;
                            pt2++;
                        } else if (*pt1 != '-' && *pt2 == '-') {
                            //							fprintf( stderr, "skip, %d (%c) - %d (%c) \n", k1, *pt1, k2, *pt2 );
                            k2++;
                            pt2++;
                        } else if (*pt1 == '-' && *pt2 != '-') {
                            //							fprintf( stderr, "skip, %d (%c) - %d (%c) \n", k1, *pt1, k2, *pt2 );
                            k1++;
                            pt1++;
                        } else if (*pt1 == '-' && *pt2 == '-') {
                            //							fprintf( stderr, "skip, %d (%c) - %d (%c) \n", k1, *pt1, k2, *pt2 );
                            k1++;
                            pt1++;
                            k2++;
                            pt2++;
                        }
                        if (k1 > end1 || k2 > end2)
                            break;
                    }
                    tmpptr = tmpptr->next;
                }
            }
        }
    }

    {
        subid = 0;
        ndone = 0;

        readloopthread_arg_t targ = {
            .opts = opts,
            .nodeid = nodeid,
            .seq1 = seq1,
            .seq2 = seq2,
            .orinum1 = orinum1,
            .orinum2 = orinum2,
            .eff1 = eff1,
            .eff2 = eff2,
            .subidpt = &subid,
            .nfiles = nfiles,
            .ndone = &ndone,
            .impmtx = impmtx,
        }; 
        readloopthread(&targ);
        npairs -= ndone;
    }

    if (npairs != 0) {
        reporterr("okashii. npairs = %d\n", npairs);
        exit(1);
    }

#if 0
	reporterr( "\n" );
	printf( "orinum1=%d, orinum2=%d\n", *orinum1, *orinum2 );
	if( *orinum1 == 0 )
	{
		fprintf( stdout, "impmtx = \n" );
		for( k2=0; k2<lgth2; k2++ )
			fprintf( stdout, "%6.3f ", (double)k2 );
		fprintf( stdout, "\n" );
		for( k1=0; k1<lgth1; k1++ )
		{
			fprintf( stdout, "%d", k1 );
			for( k2=0; k2<lgth2; k2++ )
				fprintf( stdout, "%2.1f ", impmtx[k1][k2] );
			fprintf( stdout, "\n" );
		}
		exit( 1 );
	}
#endif
}
