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
seqlen(char* seq, char gap) {
    int val = 0;
    while (*seq) {
        if (*seq++ != gap) {
            val++;
        }
    }
    return val;
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

    free(rootnode);
    free(eff);
}

void
copyWithNoGaps(char* aseq, const char* seq) {
    for (; *seq != 0; seq++) {
        if (*seq != '-')
            *aseq++ = *seq;
    }
    *aseq = 0;
}

#define SEGMENTSIZE 150

void
calcimportance_half(aln_Opts opts, int nseq, double* eff, char** seq, LocalHom** localhom, int alloclen) {
    int       i, j, pos, len;
    double*   importance;
    double    tmpdouble;
    double *  ieff, totaleff;
    int*      nogaplen;
    LocalHom* tmpptr;

    importance = AllocateDoubleVec(alloclen);
    //	reporterr("alloclen=%d, nlenmax=%d\n", alloclen, nlenmax );
    nogaplen = AllocateIntVec(nseq);
    ieff = AllocateDoubleVec(nseq);

    totaleff = 0.0;
    for (i = 0; i < nseq; i++) {
        nogaplen[i] = seqlen(seq[i], opts.gap);
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
gapireru(aln_Opts opts, char* res, char* ori, char* gt) {
    char g = 0;
    char gapchar = opts.gap;
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

void
gapcountadd(double* freq, char** seq, int nseq, double* eff, int lgth) {
    int    i;
    int    j = nseq - 1;
    double newfr = eff[j];
    double orifr = 1.0 - newfr;
    for (i = 0; i < lgth; i++) {
        freq[i] = 1.0 - freq[i];
        freq[i] *= orifr;
        if (seq[j][i] == '-')
            freq[i] += newfr;
    }
    return;
}

void
gapcountf(double* freq, char** seq, int nseq, double* eff, int lgth) {
    int    i, j;
    double fr;
    for (i = 0; i < lgth; i++) {
        fr = 0.0;
        for (j = 0; j < nseq; j++) {
            if (seq[j][i] == '-')
                fr += eff[j];
        }
        freq[i] = fr;
    }
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
        }
    }

    return;

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
}
