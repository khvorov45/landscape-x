#include "mltaln.h"

#define DEBUG 0
#define USEDISTONTREE 1

void
cpmx_calc(Context* ctx, char** seq, double** cpmx, double* eff, int lgth, int clus) {
    int    i, j, k;
    double totaleff = 0.0;

    for (i = 0; i < clus; i++)
        totaleff += eff[i];
    for (i = 0; i < ctx->nalphabets; i++)
        for (j = 0; j < lgth; j++)
            cpmx[i][j] = 0.0;
    for (j = 0; j < lgth; j++)
        for (k = 0; k < clus; k++)
            cpmx[(int)ctx->amino_n[(unsigned char)seq[k][j]]][j] += (double)eff[k] / totaleff;
}

void
cpmx_calc_add(Context* ctx, char** seq, double** cpmx, double* eff, int lgth, int clus) {
    double neweff, orieff;
    int    newmem, i, j;

    newmem = clus - 1;
    neweff = eff[clus - 1];
    orieff = 1.0 - neweff;
    for (j = 0; j < lgth; j++) {
        for (i = 0; i < ctx->nalphabets; i++)
            cpmx[i][j] *= orieff;
        cpmx[(unsigned char)ctx->amino_n[(unsigned char)seq[newmem][j]]][j] += neweff;
    }
}

void
cpmx_calc_new(Context* ctx, char** seq, double** cpmx, double* eff, int lgth, int clus) {
    int     i, j, k;
    double  feff;
    double *cpmxpt, **cpmxptpt;
    char*   seqpt;

    j = ctx->nalphabets;
    cpmxptpt = cpmx;
    while (j--) {
        cpmxpt = *cpmxptpt++;
        i = lgth;
        while (i--)
            *cpmxpt++ = 0.0;
    }
    for (k = 0; k < clus; k++) {
        feff = (double)eff[k];
        seqpt = seq[k];
        //		fprintf( stderr, "seqpt = %s, lgth=%d\n", seqpt, lgth );
        for (j = 0; j < lgth; j++) {
            cpmx[(unsigned char)ctx->amino_n[(unsigned char)*seqpt++]][j] += feff;
        }
    }
}

void
MScpmx_calc_new(Context* ctx, char** seq, double** cpmx, double* eff, int lgth, int clus) {
    int      i, j, k;
    double   feff;
    double **cpmxptpt, *cpmxpt;
    char*    seqpt;

    j = lgth;
    cpmxptpt = cpmx;
    while (j--) {
        cpmxpt = *cpmxptpt++;
        i = ctx->nalphabets;
        while (i--)
            *cpmxpt++ = 0.0;
    }
    for (k = 0; k < clus; k++) {
        feff = (double)eff[k];
        seqpt = seq[k];
        cpmxptpt = cpmx;
        j = lgth;
        while (j--)
            (*cpmxptpt++)[(int)ctx->amino_n[(unsigned char)*seqpt++]] += feff;
    }
}

void
cpmx_ribosum(Context* ctx, char** seq, char** seqr, char* dir, double** cpmx, double* eff, int lgth, int clus) {
    int      i, j, k;
    double   feff;
    double **cpmxptpt, *cpmxpt;
    char *   seqpt, *seqrpt, *dirpt;
    int      code, code1, code2;

    j = lgth;
    cpmxptpt = cpmx;
    while (j--) {
        cpmxpt = *cpmxptpt++;
        i = 37;
        while (i--)
            *cpmxpt++ = 0.0;
    }
    for (k = 0; k < clus; k++) {
        feff = (double)eff[k];
        seqpt = seq[k];
        seqrpt = seqr[k];
        dirpt = dir;
        cpmxptpt = cpmx;
        j = lgth;
        while (j--) {
#if 0
			code1 = amino_n[(int)*seqpt];
			if( code1 > 3 ) code = 36;
			else
				code = code1;
#else
            code1 = ctx->amino_n[(unsigned char)*seqpt];
            code2 = ctx->amino_n[(unsigned char)*seqrpt];
            if (code1 > 3) {
                code = 36;
            } else if (code2 > 3) {
                code = code1;
            } else if (*dirpt == '5') {
                code = 4 + code2 * 4 + code1;
            } else if (*dirpt == '3') {
                code = 20 + code2 * 4 + code1;
            } else  // if( *dirpt == 'o' ) // nai
            {
                code = code1;
            }
#endif

            seqpt++;
            seqrpt++;
            dirpt++;

            (*cpmxptpt++)[code] += feff;
        }
    }
}

int
conjuctionforgaln(int s0, int s1, char** seq, char** aseq, double* peff, double* eff, char* d) {
    int    m, k;
    char   b[B];
    double total;

#if 0
	fprintf( stderr, "s0 = %d, s1 = %d\n", s0, s1 );
#endif

    total = 0.0;
    d[0] = 0;
    for (m = s0, k = 0; m < s1; m++) {
        {
            sprintf(b, " %d", m + 1);
#if 1
            if (strlen(d) < 100)
                strcat(d, b);
#else
            strcat(d, b);
#endif
            aseq[k] = seq[m];
            peff[k] = eff[m];
            total += peff[k];
#if 0
			strcpy( aname[k], name[m] );
#endif
            k++;
        }
    }
#if 1
    for (m = 0; m < k; m++) {
        peff[m] /= total;
        //		fprintf( stderr, "peff[%d] = %f\n", m, peff[m] );
    }
#endif
    return (k);
}

void
makegrouprna(RNApair*** group, RNApair*** all, int* memlist) {
    int k, m;
    for (k = 0; (m = *memlist) != -1; memlist++, k++) {
        group[k] = all[m];
    }
}

int
fastconjuction_noweight(int* memlist, char** seq, char** aseq, double* peff, char* d) {
    int    m, k, dln;
    char   b[B];
    double total;

#if DEBUG
    fprintf(stderr, "s = %d\n", s);
#endif

    total = 0.0;
    d[0] = 0;
    dln = 0;
    for (k = 0; *memlist != -1; memlist++, k++) {
        m = *memlist;
        dln += sprintf(b, " %d", m + 1);
#if 1
        if (dln < 100)
            strcat(d, b);
#else
        strcat(d, b);
#endif
        aseq[k] = seq[m];
        peff[k] = 1.0;
        total += peff[k];
    }
#if 1
    for (m = 0; m < k; m++)
        peff[m] /= total;
#endif
    return (k);
}

int
fastconjuction_noname_kozo(int* memlist, char** seq, char** aseq, double* peff, double* eff, double* peff_kozo, double* eff_kozo, char* d) {
    int    m, k, dln;
    char   b[B];
    double total;
    double total_kozo;

#if DEBUG
    fprintf(stderr, "s = %d\n", s);
#endif

    total = 0.0;
    total_kozo = 0.0;
    d[0] = 0;
    dln = 0;
    for (k = 0; *memlist != -1; memlist++, k++) {
        m = *memlist;
        dln += sprintf(b, " %d", m + 1);
#if 1
        if (dln < 100)
            strcat(d, b);
#else
        strcat(d, b);
#endif
        aseq[k] = seq[m];
        peff[k] = eff[m];
        peff_kozo[k] = eff_kozo[m];
        total += peff[k];
        total_kozo += peff_kozo[k];
    }
#if 1
    for (m = 0; m < k; m++) {
        //		fprintf( stderr, "peff[%d] = %f\n", m, peff[m] );
        peff[m] /= total;
    }
    if (total_kozo) {
        for (m = 0; m < k; m++) {
            //			fprintf( stderr, "peff[%d] = %f\n", m, peff[m] );
            peff_kozo[m] /= total_kozo;
            if (peff_kozo[m] > 0.0)
                peff_kozo[m] += peff[m];
        }
    } else  //iranai
    {
        for (m = 0; m < k; m++) {
            peff_kozo[m] = 0.0;
        }
    }
#endif

    //	fprintf( stderr, "\n\ntbfast_total_kozo = %f\n\n", total_kozo );
    return (k);
}

#if 0
int fastconjuction_target( int *memlist, char **seq, char **aseq, double *peff, double *eff, char *d, double mineff, int *targetmap )
{
	int m, k, dln;
	char b[B];
	double total;
	int *memlistbk = memlist;

#if DEBUG
	fprintf( stderr, "s = %d\n", s );
#endif

	total = 0.0;
	d[0] = 0;
	dln = 0;
	for( k=0; *memlist!=-1; memlist++, k++ )
	{
		m = *memlist;
		dln += sprintf( b, " %d", m+1 );
#if 1
		if( dln < 100 ) strcat( d, b );
#else
		strcat( d, b );
#endif
		aseq[k] = seq[m];
		if( eff[m] < mineff )
			peff[k] = mineff;
		else
			peff[k] = eff[m];

		total += peff[k];
	}
#if 1
	for( m=0; m<k; m++ )
	{
//		fprintf( stderr, "Apr17   peff[%d] = %f\n", m, peff[m] );
		peff[m] /= total;
	}
#endif

	return( k );
}
#endif

int
fastconjuction_noname(int* memlist, char** seq, char** aseq, double* peff, double* eff, char* d, double mineff, double* oritotal) {
    int    m, k, dln;
    char   b[B];
    double total;

#if DEBUG
    fprintf(stderr, "s = %d\n", s);
#endif

    total = 0.0;
    d[0] = 0;
    dln = 0;
    for (k = 0; *memlist != -1; memlist++, k++) {
        m = *memlist;
        dln += sprintf(b, " %d", m + 1);
#if 1
        if (dln < 100)
            strcat(d, b);
#else
        strcat(d, b);
#endif
        aseq[k] = seq[m];
        if (eff[m] < mineff)
            peff[k] = mineff;
        else
            peff[k] = eff[m];

        total += peff[k];
    }
    if (oritotal)
        *oritotal = total;
#if 1
    for (m = 0; m < k; m++) {
        //		fprintf( stderr, "Apr17   peff[%d] = %20.10f\n", m, peff[m] );
        peff[m] /= total;
    }
#endif
    return (k);
}

int
fastconjuction(int* memlist, char** seq, char** aseq, double* peff, double* eff, char* d) {
    int    m, k, dln;
    char   b[B];
    double total;

#if DEBUG
    fprintf(stderr, "s = %d\n", s);
#endif

    total = 0.0;
    d[0] = 0;
    dln = 0;
    for (k = 0; *memlist != -1; memlist++, k++) {
        m = *memlist;
        dln += sprintf(b, " %d", m + 1);
#if 1
        if (dln < 100)
            strcat(d, b);
#else
        strcat(d, b);
#endif
        aseq[k] = seq[m];
        peff[k] = eff[m];
        total += peff[k];
#if 0
			strcpy( aname[k], name[m] );
#endif
    }
#if 1
    for (m = 0; m < k; m++)
        peff[m] /= total;
#endif
    return (k);
}

void
chardelete(char* seq, int d) {
    char b[N];

    strcpy(b, seq + d + 1);
    strcpy(seq + d, b);
}

void
nodeFromABranch(int nseq, int* result, int** pairwisenode, int*** topol, int step, int num) {
    int  i, s, count;
    int* innergroup;
    int* outergroup1;
#if 0
	int outergroup2[nseq];
	int table[nseq];
#else
    static int* outergroup2 = NULL;
    static int* table = NULL;
    if (outergroup2 == NULL) {
        outergroup2 = AllocateIntVec(nseq);
        table = AllocateIntVec(nseq);
    }
#endif
    innergroup = topol[step][num];
    outergroup1 = topol[step][!num];
    for (i = 0; i < nseq; i++)
        table[i] = 1;
    for (i = 0; (s = innergroup[i]) > -1; i++)
        table[s] = 0;
    for (i = 0; (s = outergroup1[i]) > -1; i++)
        table[s] = 0;
    for (i = 0, count = 0; i < nseq; i++) {
        if (table[i]) {
            outergroup2[count] = i;
            count++;
        }
    }
    outergroup2[count] = -1;

    for (i = 0; (s = innergroup[i]) > -1; i++) {
        result[s] = pairwisenode[s][outergroup1[0]]
            + pairwisenode[s][outergroup2[0]]
            - pairwisenode[outergroup1[0]][outergroup2[0]] - 1;
        result[s] /= 2;
    }
    for (i = 0; (s = outergroup1[i]) > -1; i++) {
        result[s] = pairwisenode[s][outergroup2[0]]
            + pairwisenode[s][innergroup[0]]
            - pairwisenode[innergroup[0]][outergroup2[0]] + 1;
        result[s] /= 2;
    }
    for (i = 0; (s = outergroup2[i]) > -1; i++) {
        result[s] = pairwisenode[s][outergroup1[0]]
            + pairwisenode[s][innergroup[0]]
            - pairwisenode[innergroup[0]][outergroup1[0]] + 1;
        result[s] /= 2;
    }

#if 0
	for( i=0; i<nseq; i++ ) 
		fprintf( stderr, "result[%d] = %d\n", i+1, result[i] );
#endif
}

void
OneClusterAndTheOther_fast(int locnjob, int* memlist1, int* memlist2, int* s1, int* s2, char* pair, int*** topol, int step, int branch, double** smalldistmtx, double* distontree) {
    int i, k, j;
    int r1;
    //	char *pair;

    //	pair = calloc( locnjob, sizeof( char ) );

    for (i = 0; i < locnjob; i++)
        pair[i] = 0;
    for (i = 0, k = 0; (r1 = topol[step][branch][i]) > -1; i++) {
        pair[r1] = 1;
        memlist1[k++] = r1;
    }
    memlist1[k] = -1;

    for (i = 0, k = 0; i < locnjob; i++) {
        if (!pair[i]) {
            memlist2[k++] = i;
        }
    }
    memlist2[k] = -1;

    *s1 = memlist1[0];
    *s2 = memlist2[0];

    if (smalldistmtx) {
        int im, jm;
#if USEDISTONTREE
        for (i = 0; (im = memlist1[i]) != -1; i++)
            for (j = 0; (jm = memlist2[j]) != -1; j++)
                smalldistmtx[i][j] = distontree[im] + distontree[jm];
#else
        for (i = 0; (im = memlist1[i]) != -1; i++)
            for (j = 0; (jm = memlist2[j]) != -1; j++)
                smalldistmtx[i][j] = distmtx[MIN(im, jm)][MAX(im, jm)];
#endif

#if 0
		reporterr( "\n" );
		for( i=0; (im=memlist1[i])!=-1; i++ ) for( j=0; (jm=memlist2[j])!=-1; j++ )
			reporterr( "smalldistmtx[%d][%d] = %f\n", i, j, smalldistmtx[i][j] );


		for( i=0; (im=memlist1[i])!=-1; i++ ) for( j=0; (jm=memlist2[j])!=-1; j++ )
			smalldistmtx[i][j] = distmtx[MIN(im,jm)][MAX(im,jm)];

		for( i=0; (im=memlist1[i])!=-1; i++ ) for( j=0; (jm=memlist2[j])!=-1; j++ )
			reporterr( "old smalldistmtx[%d][%d] = %f\n", i, j, smalldistmtx[i][j] );
if( im > 10 && jm > 10 ) exit( 1 );
#endif
    }
}

void
makeEffMtx(int nseq, double** mtx, double* vec) {
    int i, j;
    for (i = 0; i < nseq; i++)
        for (j = 0; j < nseq; j++)
            mtx[i][j] = vec[i] * vec[j];
}

int
msshrinklocalhom_fast_target(int* memlist1, int* memlist2, LocalHom** localhom, LocalHom*** localhomshrink, char* swaplist, int* targetmap) {
    int m1, k1, m2, k2, t1, i2;

    for (k1 = 0; (m1 = memlist1[k1]) != -1; k1++) {
        if (targetmap[m1] == -1) {
            swaplist[k1] = 1;
            //			swaplist[k1] = 0; // DAME!!!
            for (k2 = 0; (m2 = memlist2[k2]) != -1; k2++) {
                if (targetmap[m2] == -1) {
                    localhomshrink[k1][k2] = NULL;
                    continue;
                }

                t1 = targetmap[m2];  // start1 <-> start2, end1 <-> end2
                i2 = m1;

                if (localhom[t1][i2].opt == -1)
                    localhomshrink[k1][k2] = NULL;
                else
                    localhomshrink[k1][k2] = localhom[t1] + i2;
            }
        } else {
            swaplist[k1] = 0;
            for (k2 = 0; (m2 = memlist2[k2]) != -1; k2++) {
                t1 = targetmap[m1];
                i2 = m2;

                if (localhom[t1][i2].opt == -1)
                    localhomshrink[k1][k2] = NULL;
                else
                    localhomshrink[k1][k2] = localhom[t1] + i2;
            }
        }
    }
    return (0);
}

int
msshrinklocalhom_fast_half(int* memlist1, int* memlist2, LocalHom** localhom, LocalHom*** localhomshrink) {
    int m1, k1, m2, k2;

    for (k1 = 0; (m1 = memlist1[k1]) != -1; k1++) {
        for (k2 = 0; (m2 = memlist2[k2]) != -1; k2++) {
            if (m1 < m2) {
                if (localhom[m1][m2 - m1].opt == -1)
                    localhomshrink[k1][k2] = NULL;
                else
                    localhomshrink[k1][k2] = localhom[m1] + m2 - m1;
            } else {
                if (localhom[m2][m1 - m2].opt == -1)
                    localhomshrink[k1][k2] = NULL;
                else
                    localhomshrink[k1][k2] = localhom[m2] + m1 - m2;
            }
        }
    }
    return (0);
}

int
msshrinklocalhom_fast(int* memlist1, int* memlist2, LocalHom** localhom, LocalHom*** localhomshrink) {
    int m1, k1, m2, k2;

    for (k1 = 0; (m1 = memlist1[k1]) != -1; k1++) {
        for (k2 = 0; (m2 = memlist2[k2]) != -1; k2++) {
            if (localhom[m1][m2].opt == -1)
                localhomshrink[k1][k2] = NULL;
            else
                localhomshrink[k1][k2] = localhom[m1] + m2;
        }
    }
    return (0);
}
int
fastshrinklocalhom_one(int* mem1, int* mem2, int norg, LocalHom** localhom, LocalHom*** localhomshrink) {
    int  k1, k2;
    int *intpt1, *intpt2;

    for (intpt1 = mem1, k1 = 0; *intpt1 != -1; intpt1++, k1++) {
        for (intpt2 = mem2, k2 = 0; *intpt2 != -1; intpt2++, k2++) {
            if (*intpt2 != norg) {
                fprintf(stderr, "ERROR! *intpt2 = %d\n", *intpt2);
                exit(1);
            }
            if (localhom[*intpt1][0].opt == -1)
                localhomshrink[k1][k2] = NULL;
            else
                localhomshrink[k1][k2] = localhom[*intpt1];
        }
    }
    return (0);
}

int
fastshrinklocalhom(int* mem1, int* mem2, LocalHom** localhom, LocalHom*** localhomshrink) {
    int  k1, k2;
    int *intpt1, *intpt2;

    for (intpt1 = mem1, k1 = 0; *intpt1 != -1; intpt1++, k1++) {
        for (intpt2 = mem2, k2 = 0; *intpt2 != -1; intpt2++, k2++) {
            if (localhom[*intpt1][*intpt2].opt == -1)
                localhomshrink[k1][k2] = NULL;
            else
                localhomshrink[k1][k2] = localhom[*intpt1] + *intpt2;

            //			if( localhomshrink[k1][k2] != NULL )
            //				printf( "ori localhomshrink[%d][%d].opt = %f\n", k1, k2, localhomshrink[k1][k2]->opt );
        }
    }
    return (0);
}

int
fastshrinklocalhom_half_seed(int* mem1, int* mem2, int nseed, int* posinlsh1, int* posinlsh2, LocalHom** localhom, LocalHom*** localhomshrink) {
    int  k1, k2, sk1, sk2;
    int *intpt1, *intpt2;

    for (intpt1 = mem1, k1 = 0, sk1 = 0; *intpt1 != -1; intpt1++, k1++) {
        if (*intpt1 >= nseed)
            posinlsh1[k1] = -1;
        else
            posinlsh1[k1] = sk1++;
    }
    for (intpt2 = mem2, k2 = 0, sk2 = 0; *intpt2 != -1; intpt2++, k2++) {
        if (*intpt2 >= nseed)
            posinlsh2[k2] = -1;
        else
            posinlsh2[k2] = sk2++;
    }

    for (intpt1 = mem1, sk1 = 0; *intpt1 != -1; intpt1++) {
        if (*intpt1 >= nseed)
            continue;
        for (intpt2 = mem2, sk2 = 0; *intpt2 != -1; intpt2++) {
            if (*intpt2 >= nseed)
                continue;
            if (*intpt1 < *intpt2) {
                if (localhom[*intpt1][*intpt2 - *intpt1].opt == -1)
                    localhomshrink[sk1][sk2] = NULL;
                else
                    localhomshrink[sk1][sk2] = localhom[*intpt1] + *intpt2 - *intpt1;
            } else {
                if (localhom[*intpt2][*intpt1 - *intpt2].opt == -1)
                    localhomshrink[sk1][sk2] = NULL;
                else
                    localhomshrink[sk1][sk2] = localhom[*intpt2] + *intpt1 - *intpt2;
            }

            //			if( localhomshrink[k1][k2] != NULL )
            //				printf( "ori localhomshrink[%d][%d].opt = %f, .importance=%f\n", k1, k2, localhomshrink[k1][k2]->opt, localhomshrink[k1][k2]->importance );
            sk2++;
        }
        sk1++;
    }
    return (0);
}

int
fastshrinklocalhom_half(int* mem1, int* mem2, LocalHom** localhom, LocalHom*** localhomshrink) {
    int  k1, k2;
    int *intpt1, *intpt2;

    for (intpt1 = mem1, k1 = 0; *intpt1 != -1; intpt1++, k1++) {
        for (intpt2 = mem2, k2 = 0; *intpt2 != -1; intpt2++, k2++) {
            if (*intpt1 < *intpt2) {
                if (localhom[*intpt1][*intpt2 - *intpt1].opt == -1)
                    localhomshrink[k1][k2] = NULL;
                else
                    localhomshrink[k1][k2] = localhom[*intpt1] + *intpt2 - *intpt1;
            } else {
                if (localhom[*intpt2][*intpt1 - *intpt2].opt == -1)
                    localhomshrink[k1][k2] = NULL;
                else
                    localhomshrink[k1][k2] = localhom[*intpt2] + *intpt1 - *intpt2;
            }

            //			if( localhomshrink[k1][k2] != NULL )
            //				printf( "ori localhomshrink[%d][%d].opt = %f, .importance=%f\n", k1, k2, localhomshrink[k1][k2]->opt, localhomshrink[k1][k2]->importance );
        }
    }
    return (0);
}

int
fastshrinklocalhom_target(int* mem1, int* mem2, LocalHom** localhom, LocalHom*** localhomshrink, char* swaplist, int* targetmap) {
    int  k1, k2;
    int *intpt1, *intpt2;
    int  t1, i2;

    for (intpt1 = mem1, k1 = 0; *intpt1 != -1; intpt1++, k1++) {
        if (targetmap[*intpt1] == -1) {
            swaplist[k1] = 1;
            //			swaplist[k1] = 0; // DAME!!!
            for (intpt2 = mem2, k2 = 0; *intpt2 != -1; intpt2++, k2++) {
                if (targetmap[*intpt2] == -1) {
                    localhomshrink[k1][k2] = NULL;
                    continue;
                }

                t1 = targetmap[*intpt2];  // end1<->end2, start1<->start2
                i2 = *intpt1;

                if (localhom[t1][i2].opt == -1)
                    localhomshrink[k1][k2] = NULL;
                else
                    localhomshrink[k1][k2] = localhom[t1] + i2;

                //				if( localhomshrink[k1][k2] != NULL )
                //					printf( "localhomshrink[%d][%d].opt = %f\n", k1, k2, localhomshrink[k1][k2]->opt );
                //				else
                //					printf( "localhomshrink[%d][%d] = NULL\n", k1, k2 );
            }
        } else {
            swaplist[k1] = 0;
            for (intpt2 = mem2, k2 = 0; *intpt2 != -1; intpt2++, k2++) {
                t1 = targetmap[*intpt1];
                i2 = *intpt2;

                if (localhom[t1][i2].opt == -1)
                    localhomshrink[k1][k2] = NULL;
                else
                    localhomshrink[k1][k2] = localhom[t1] + i2;

                //				if( localhomshrink[k1][k2] != NULL )
                //					printf( "localhomshrink[%d][%d].opt = %f\n", k1, k2, localhomshrink[k1][k2]->opt );
                //				else
                //					printf( "localhomshrink[%d][%d] = NULL\n", k1, k2 );
            }
        }
    }
    return (0);
}

int
msfastshrinklocalhom(int* mem1, int* mem2, LocalHom** localhom, LocalHom*** localhomshrink) {
    int  k1, k2;
    int *intpt1, *intpt2;
    int  m1, m2;

    for (intpt1 = mem1, k1 = 0; *intpt1 != -1; intpt1++, k1++) {
        for (intpt2 = mem2, k2 = 0; *intpt2 != -1; intpt2++, k2++) {
            m1 = MIN(*intpt1, *intpt2);
            m2 = MAX(*intpt1, *intpt2);
            if (localhom[m1][m2].opt == -1)
                localhomshrink[k1][k2] = NULL;
            else
                localhomshrink[k1][k2] = localhom[m1] + m2;
        }
    }
    return (0);
}
