#include "mltaln.h"

#define DEBUG 0
#define IODEBUG 0

#define NODIST -9999

static char* whereispairalign;
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

typedef struct _lastres {
    int   score;
    int   start1;
    int   start2;
    char* aln1;
    char* aln2;
} Lastres;

typedef struct _reg {
    int start;
    int end;
} Reg;

typedef struct _aln {
    int  nreg;
    Reg* reg1;
    Reg* reg2;
} Aln;

typedef struct _lastresx {
    int  score;
    int  naln;
    Aln* aln;
} Lastresx;

typedef struct lastcallthread_arg_t {
    aln_Opts   opts;
    Context*   ctx;
    int        nq, nd;
    char**     dseq;
    char**     qseq;
    Lastresx** lastresx;
} lastcallthread_arg_t;

static void
t2u(char* seq) {
    while (*seq) {
        if (*seq == 'A')
            *seq = 'a';
        else if (*seq == 'a')
            *seq = 'a';
        else if (*seq == 'T')
            *seq = 'u';
        else if (*seq == 't')
            *seq = 'u';
        else if (*seq == 'U')
            *seq = 'u';
        else if (*seq == 'u')
            *seq = 'u';
        else if (*seq == 'G')
            *seq = 'g';
        else if (*seq == 'g')
            *seq = 'g';
        else if (*seq == 'C')
            *seq = 'c';
        else if (*seq == 'c')
            *seq = 'c';
        else
            *seq = 'n';
        seq++;
    }
}

static int
removex(char* d, char* m) {
    int val = 0;
    while (*m != 0) {
        if (*m == 'X' || *m == 'x') {
            m++;
            val++;
        } else {
            *d++ = *m++;
        }
    }
    *d = 0;
    return (val);
}

static void
putlocalhom_last(Context* ctx, char* s1, char* s2, LocalHom* localhompt, Lastresx* lastresx) {
    char *    pt1, *pt2;
    int       naln, nreg;
    int       iscore;
    int       isumscore;
    int       sumoverlap;
    LocalHom* tmppt = localhompt;
    LocalHom* tmppt2;
    LocalHom* localhompt0;
    Reg *     rpt1, *rpt2;
    Aln*      apt;
    int       nlocalhom = 0;
    int       len;

    naln = lastresx->naln;
    apt = lastresx->aln;

    if (naln == 0)
        return;
    while (naln--) {
        rpt1 = apt->reg1;
        rpt2 = apt->reg2;
        nreg = apt->nreg;
        isumscore = 0;
        sumoverlap = 0;
        while (nreg--) {
            if (nlocalhom++ > 0) {
                tmppt->next = (LocalHom*)calloc(1, sizeof(LocalHom));
                tmppt = tmppt->next;
                tmppt->next = NULL;
            }
            tmppt->start1 = rpt1->start;
            tmppt->start2 = rpt2->start;
            tmppt->end1 = rpt1->end;
            tmppt->end2 = rpt2->end;
            tmppt->korh = 'h';
            if (rpt1 == apt->reg1)
                localhompt0 = tmppt;  // ?

            len = tmppt->end1 - tmppt->start1 + 1;

            pt1 = s1 + tmppt->start1;
            pt2 = s2 + tmppt->start2;
            iscore = 0;
            while (len--) {
                iscore += ctx->n_dis[(int)ctx->amino_n[(unsigned char)*pt1++]][(int)ctx->amino_n[(unsigned char)*pt2++]];  // - offset はいらないかも
            }

            if (ctx->divpairscore) {
                tmppt->overlapaa = tmppt->end2 - tmppt->start2 + 1;
                tmppt->opt = (double)iscore / tmppt->overlapaa * 5.8 / 600;
            } else {
                isumscore += iscore;
                sumoverlap += tmppt->end2 - tmppt->start2 + 1;
            }
            rpt1++;
            rpt2++;
        }

        if (!ctx->divpairscore) {
            for (tmppt2 = localhompt0; tmppt2; tmppt2 = tmppt2->next) {
                tmppt2->overlapaa = sumoverlap;
                tmppt2->opt = (double)isumscore * 5.8 / (600 * sumoverlap);
            }
        }
        apt++;
    }
}

static int
countcomma(char* s) {
    int v = 0;
    while (*s)
        if (*s++ == ',')
            v++;
    return (v);
}

static double
recallpairfoldalign(Context* ctx, char** mseq1, char** mseq2, int m1, int m2, int* of1pt, int* of2pt, int alloclen) {
    static FILE* fp = NULL;
    double       value;
    char*        aln1;
    char*        aln2;
    int          of1tmp, of2tmp;

    if (fp == NULL) {
        fp = fopen("_foldalignout", "r");
        if (fp == NULL) {
            fprintf(stderr, "Cannot open _foldalignout\n");
            exit(1);
        }
    }

    aln1 = calloc(alloclen, sizeof(char));
    aln2 = calloc(alloclen, sizeof(char));

    readpairfoldalign(fp, *mseq1, *mseq2, aln1, aln2, m1, m2, &of1tmp, &of2tmp, alloclen);

    if (strstr(foldalignopt, "-global")) {
        fprintf(stderr, "Calling G__align11\n");
        value = G__align11(ctx, ctx->n_dis_consweight_multi, mseq1, mseq2, alloclen, ctx->outgap, ctx->outgap);
        *of1pt = 0;
        *of2pt = 0;
    } else {
        fprintf(stderr, "Calling L__align11\n");
        value = L__align11(ctx, ctx->n_dis_consweight_multi, 0.0, mseq1, mseq2, alloclen, of1pt, of2pt);
    }

    if (aln1[0] == 0) {
        fprintf(stderr, "FOLDALIGN returned no alignment between %d and %d.  Sequence alignment is used instead.\n", m1 + 1, m2 + 1);
    } else {
        strcpy(*mseq1, aln1);
        strcpy(*mseq2, aln2);
        *of1pt = of1tmp;
        *of2pt = of2tmp;
    }

    free(aln1);
    free(aln2);

    return (value);
}

static void
block2reg(char* block, Reg* reg1, Reg* reg2, int start1, int start2) {
    Reg*  rpt1;
    char *tpt, *npt;
    int   pos1, pos2;
    int   len, glen1, glen2;
    pos1 = start1;
    pos2 = start2;
    rpt1 = reg1;
    while (block) {
        block++;
        //		fprintf( stderr, "block = %s\n", block );
        tpt = strchr(block, ':');
        npt = strchr(block, ',');
        if (!tpt || tpt > npt) {
            len = atoi(block);
            reg1->start = pos1;
            reg2->start = pos2;
            pos1 += len - 1;
            pos2 += len - 1;
            reg1->end = pos1;
            reg2->end = pos2;
            //			fprintf( stderr, "in loop reg1: %d-%d\n", reg1->start, reg1->end );
            //			fprintf( stderr, "in loop reg2: %d-%d\n", reg2->start, reg2->end );
            reg1++;
            reg2++;
        } else {
            sscanf(block, "%d:%d", &glen1, &glen2);
            pos1 += glen1 + 1;
            pos2 += glen2 + 1;
        }
        block = npt;
    }
    reg1->start = reg1->end = reg2->start = reg2->end = -1;

    while (rpt1->start != -1) {
        //		fprintf( stderr, "reg1: %d-%d\n", rpt1->start, rpt1->end );
        //		fprintf( stderr, "reg2: %d-%d\n", rpt2->start, rpt2->end );
        rpt1++;
    }
    //	*apt1 = *apt2 = 0;
    //	fprintf( stderr, "aln1 = %s\n", aln1 );
    //	fprintf( stderr, "aln2 = %s\n", aln2 );
}

static void
readlastresx_singleq(aln_Opts opts, Context* ctx, FILE* fp, int nameq, Lastresx** lastresx) {
    char* gett;
    Aln*  tmpaln;
    int   prevnaln, naln, nreg;
#if 0
	int i, pstart, pend, end1, end2;
#endif
    int  score, name1, start1, alnSize1, seqSize1;
    int  name2, start2, alnSize2, seqSize2;
    char strand1, strand2;
    int  includeintoscore;
    gett = calloc(10000, sizeof(char));

    //	fprintf( stderr, "seq2[0] = %s\n", seq2[0] );
    //	fprintf( stderr, "seq1[0] = %s\n", seq1[0] );

    while (1) {
        fgets(gett, 9999, fp);
        if (feof(fp))
            break;
        if (gett[0] == '#')
            continue;
        //		fprintf( stdout, "gett = %s\n", gett );
        if (gett[strlen(gett) - 1] != '\n') {
            fprintf(stderr, "Too long line?\n");
            exit(1);
        }

        sscanf(gett, "%d %d %d %d %c %d %d %d %d %c %d", &score, &name1, &start1, &alnSize1, &strand1, &seqSize1, &name2, &start2, &alnSize2, &strand2, &seqSize2);

        if (opts.alg == 'R' && name2 <= name1)
            continue;
        if (name2 != nameq) {
            fprintf(stderr, "BUG!!!\n");
            exit(1);
        }

        //		if( lastresx[name1][name2].score ) continue; // dame!!!!

        prevnaln = lastresx[name1][name2].naln;
#if 0
		for( i=0; i<prevnaln; i++ )
		{
			nreg = lastresx[name1][name2].aln[i].nreg;

			pstart = lastresx[name1][name2].aln[i].reg1[0].start + 0;
			pend   = lastresx[name1][name2].aln[i].reg1[nreg-1].end - 0;
			end1 = start1 + alnSize1;
//			fprintf( stderr, "pstart = %d, pend = %d\n", pstart, pend );
			if( pstart <= start1 && start1 <= pend && pend - start1 > 1 ) break;
			if( pstart <= end1   && end1   <= pend && end1 - pstart > 1 ) break;

			pstart = lastresx[name1][name2].aln[i].reg2[0].start + 0;
			pend   = lastresx[name1][name2].aln[i].reg2[nreg-1].end - 0;
			end2 = start2 + alnSize2;
//			fprintf( stderr, "pstart = %d, pend = %d\n", pstart, pend );
			if( pstart <= start2 && start2 <= pend && pend - start2 > 1 ) break;
			if( pstart <= end2   && end2   <= pend && end2 - pstart > 1 ) break;
		}
		includeintoscore = ( i == prevnaln );
#else
        if (prevnaln)
            includeintoscore = 0;
        else
            includeintoscore = 1;
#endif
        if (!includeintoscore && !lastsubopt)
            continue;

        naln = prevnaln + 1;
        lastresx[name1][name2].naln = naln;
        //		fprintf( stderr, "OK! add this alignment to hat3, %d-%d, naln = %d->%d\n", name1, name2, prevnaln, naln );

        if ((tmpaln = (Aln*)realloc(lastresx[name1][name2].aln, (naln) * sizeof(Aln))) == NULL)  // yoyu nashi
        {
            fprintf(stderr, "Cannot reallocate lastresx[][].aln\n");
            exit(1);
        } else
            lastresx[name1][name2].aln = tmpaln;

        nreg = countcomma(gett) / 2 + 1;
        lastresx[name1][name2].aln[prevnaln].nreg = nreg;
        //		lastresx[name1][name2].aln[naln].nreg = -1;
        //		lastresx[name1][name2].aln[naln].reg1 = NULL;
        //		lastresx[name1][name2].aln[naln].reg2 = NULL;
        //		fprintf( stderr, "name1=%d, name2=%d, nreg=%d, prevnaln=%d\n", name1, name2, nreg, prevnaln );

        if ((lastresx[name1][name2].aln[prevnaln].reg1 = (Reg*)calloc(nreg + 1, sizeof(Reg))) == NULL)  // yoyu nashi
        {
            fprintf(stderr, "Cannot reallocate lastresx[][].reg2\n");
            exit(1);
        }

        if ((lastresx[name1][name2].aln[prevnaln].reg2 = (Reg*)calloc(nreg + 1, sizeof(Reg))) == NULL)  // yoyu nashi
        {
            fprintf(stderr, "Cannot reallocate lastresx[][].reg2\n");
            exit(1);
        }

        //		lastresx[name1][name2].aln[prevnaln].reg1[0].start = -1; // iranai?
        //		lastresx[name1][name2].aln[prevnaln].reg2[0].start = -1; // iranai?
        block2reg(strrchr(gett, '\t'), lastresx[name1][name2].aln[prevnaln].reg1, lastresx[name1][name2].aln[prevnaln].reg2, start1, start2);

        if (includeintoscore) {
            if (lastresx[name1][name2].score)
                score += ctx->penalty;
            lastresx[name1][name2].score += score;
        }

        //		fprintf( stderr, "score(%d,%d) = %d\n", name1, name2, lastresx[name1][name2].score );
    }
    free(gett);
}

#ifdef enablemultithread
#if 0
static void readlastresx_group( FILE *fp, Lastresx **lastresx )
{
	char *gett;
	Aln *tmpaln;
	int prevnaln, naln, nreg;
#if 0
	int i, pstart, pend, end1, end2;
#endif
	int score, name1, start1, alnSize1, seqSize1;
	int        name2, start2, alnSize2, seqSize2;
	char strand1, strand2;
	int includeintoscore;
	gett = calloc( 10000, sizeof( char ) );

//	fprintf( stderr, "seq2[0] = %s\n", seq2[0] );
//	fprintf( stderr, "seq1[0] = %s\n", seq1[0] );

	while( 1 )
	{
		fgets( gett, 9999, fp );
		if( feof( fp ) ) break;
		if( gett[0] == '#' ) continue;
//		fprintf( stdout, "gett = %s\n", gett );
		if( gett[strlen(gett)-1] != '\n' )
		{
			fprintf( stderr, "Too long line?\n" );
			exit( 1 );
		}

		sscanf( gett, "%d %d %d %d %c %d %d %d %d %c %d", 
					&score, &name1, &start1, &alnSize1, &strand1, &seqSize1,
					        &name2, &start2, &alnSize2, &strand2, &seqSize2 );

		if( alg == 'R' && name2 <= name1 ) continue;

//		if( lastresx[name1][name2].score ) continue; // dame!!!!

		prevnaln = lastresx[name1][name2].naln;
#if 0
		for( i=0; i<prevnaln; i++ )
		{
			nreg = lastresx[name1][name2].aln[i].nreg;

			pstart = lastresx[name1][name2].aln[i].reg1[0].start + 0;
			pend   = lastresx[name1][name2].aln[i].reg1[nreg-1].end - 0;
			end1 = start1 + alnSize1;
//			fprintf( stderr, "pstart = %d, pend = %d\n", pstart, pend );
			if( pstart <= start1 && start1 <= pend && pend - start1 > 3 ) break;
			if( pstart <= end1   && end1   <= pend && end1 - pstart > 3 ) break;

			pstart = lastresx[name1][name2].aln[i].reg2[0].start + 0;
			pend   = lastresx[name1][name2].aln[i].reg2[nreg-1].end - 0;
			end2 = start2 + alnSize2;
//			fprintf( stderr, "pstart = %d, pend = %d\n", pstart, pend );
			if( pstart <= start2 && start2 <= pend && pend - start2 > 3 ) break;
			if( pstart <= end2   && end2   <= pend && end2 - pstart > 3 ) break;
		}
		includeintoscore = ( i == prevnaln );
#else
		if( prevnaln ) includeintoscore = 0;
		else includeintoscore = 1;
#endif
		if( !includeintoscore && !lastsubopt )
			continue;

		naln = prevnaln + 1;
		lastresx[name1][name2].naln = naln;
//		fprintf( stderr, "OK! add this alignment to hat3, %d-%d, naln = %d->%d\n", name1, name2, prevnaln, naln );

		if( ( tmpaln = (Aln *)realloc( lastresx[name1][name2].aln, (naln) * sizeof( Aln ) ) ) == NULL ) // yoyu nashi
		{
			fprintf( stderr, "Cannot reallocate lastresx[][].aln\n" );
			exit( 1 );
		}
		else
			lastresx[name1][name2].aln = tmpaln;



		nreg = countcomma( gett )/2 + 1;
		lastresx[name1][name2].aln[prevnaln].nreg = nreg;
//		lastresx[name1][name2].aln[naln].nreg = -1;
//		lastresx[name1][name2].aln[naln].reg1 = NULL;
//		lastresx[name1][name2].aln[naln].reg2 = NULL;
//		fprintf( stderr, "name1=%d, name2=%d, nreg=%d, prevnaln=%d\n", name1, name2, nreg, prevnaln );

		if( ( lastresx[name1][name2].aln[prevnaln].reg1 = (Reg *)calloc( nreg+1, sizeof( Reg ) ) ) == NULL ) // yoyu nashi
		{
			fprintf( stderr, "Cannot reallocate lastresx[][].reg2\n" );
			exit( 1 );
		}

		if( ( lastresx[name1][name2].aln[prevnaln].reg2 = (Reg *)calloc( nreg+1, sizeof( Reg ) ) ) == NULL ) // yoyu nashi
		{
			fprintf( stderr, "Cannot reallocate lastresx[][].reg2\n" );
			exit( 1 );
		}

//		lastresx[name1][name2].aln[prevnaln].reg1[0].start = -1; // iranai?
//		lastresx[name1][name2].aln[prevnaln].reg2[0].start = -1; // iranai?
		block2reg( strrchr( gett, '\t' ), lastresx[name1][name2].aln[prevnaln].reg1, lastresx[name1][name2].aln[prevnaln].reg2, start1, start2 );

		if( includeintoscore )
		{
			if( lastresx[name1][name2].score ) score += penalty;
			lastresx[name1][name2].score += score;
		}

//		fprintf( stderr, "score(%d,%d) = %d\n", name1, name2, lastresx[name1][name2].score );
	}
	free( gett );
}
#endif
#endif

static void
readlastresx(aln_Opts opts, Context* ctx, FILE* fp, Lastresx** lastresx) {
    char* gett;
    Aln*  tmpaln;
    int   prevnaln, naln, nreg;
#if 0
	int i, pstart, pend, end1, end2;
#endif
    int  score, name1, start1, alnSize1, seqSize1;
    int  name2, start2, alnSize2, seqSize2;
    char strand1, strand2;
    int  includeintoscore;
    gett = calloc(10000, sizeof(char));

    //	fprintf( stderr, "seq2[0] = %s\n", seq2[0] );
    //	fprintf( stderr, "seq1[0] = %s\n", seq1[0] );

    while (1) {
        fgets(gett, 9999, fp);
        if (feof(fp))
            break;
        if (gett[0] == '#')
            continue;
        //		fprintf( stdout, "gett = %s\n", gett );
        if (gett[strlen(gett) - 1] != '\n') {
            fprintf(stderr, "Too long line?\n");
            exit(1);
        }

        sscanf(gett, "%d %d %d %d %c %d %d %d %d %c %d", &score, &name1, &start1, &alnSize1, &strand1, &seqSize1, &name2, &start2, &alnSize2, &strand2, &seqSize2);

        if (opts.alg == 'R' && name2 <= name1)
            continue;

        //		if( lastresx[name1][name2].score ) continue; // dame!!!!

        prevnaln = lastresx[name1][name2].naln;
#if 0
		for( i=0; i<prevnaln; i++ )
		{
			nreg = lastresx[name1][name2].aln[i].nreg;

			pstart = lastresx[name1][name2].aln[i].reg1[0].start + 0;
			pend   = lastresx[name1][name2].aln[i].reg1[nreg-1].end - 0;
			end1 = start1 + alnSize1;
//			fprintf( stderr, "pstart = %d, pend = %d\n", pstart, pend );
			if( pstart <= start1 && start1 <= pend && pend - start1 > 3 ) break;
			if( pstart <= end1   && end1   <= pend && end1 - pstart > 3 ) break;

			pstart = lastresx[name1][name2].aln[i].reg2[0].start + 0;
			pend   = lastresx[name1][name2].aln[i].reg2[nreg-1].end - 0;
			end2 = start2 + alnSize2;
//			fprintf( stderr, "pstart = %d, pend = %d\n", pstart, pend );
			if( pstart <= start2 && start2 <= pend && pend - start2 > 3 ) break;
			if( pstart <= end2   && end2   <= pend && end2 - pstart > 3 ) break;
		}
		includeintoscore = ( i == prevnaln );
#else
        if (prevnaln)
            includeintoscore = 0;
        else
            includeintoscore = 1;
#endif
        if (!includeintoscore && !lastsubopt)
            continue;

        naln = prevnaln + 1;
        lastresx[name1][name2].naln = naln;
        //		fprintf( stderr, "OK! add this alignment to hat3, %d-%d, naln = %d->%d\n", name1, name2, prevnaln, naln );

        if ((tmpaln = (Aln*)realloc(lastresx[name1][name2].aln, (naln) * sizeof(Aln))) == NULL)  // yoyu nashi
        {
            fprintf(stderr, "Cannot reallocate lastresx[][].aln\n");
            exit(1);
        } else
            lastresx[name1][name2].aln = tmpaln;

        nreg = countcomma(gett) / 2 + 1;
        lastresx[name1][name2].aln[prevnaln].nreg = nreg;
        //		lastresx[name1][name2].aln[naln].nreg = -1;
        //		lastresx[name1][name2].aln[naln].reg1 = NULL;
        //		lastresx[name1][name2].aln[naln].reg2 = NULL;
        //		fprintf( stderr, "name1=%d, name2=%d, nreg=%d, prevnaln=%d\n", name1, name2, nreg, prevnaln );

        if ((lastresx[name1][name2].aln[prevnaln].reg1 = (Reg*)calloc(nreg + 1, sizeof(Reg))) == NULL)  // yoyu nashi
        {
            fprintf(stderr, "Cannot reallocate lastresx[][].reg2\n");
            exit(1);
        }

        if ((lastresx[name1][name2].aln[prevnaln].reg2 = (Reg*)calloc(nreg + 1, sizeof(Reg))) == NULL)  // yoyu nashi
        {
            fprintf(stderr, "Cannot reallocate lastresx[][].reg2\n");
            exit(1);
        }

        //		lastresx[name1][name2].aln[prevnaln].reg1[0].start = -1; // iranai?
        //		lastresx[name1][name2].aln[prevnaln].reg2[0].start = -1; // iranai?
        block2reg(strrchr(gett, '\t'), lastresx[name1][name2].aln[prevnaln].reg1, lastresx[name1][name2].aln[prevnaln].reg2, start1, start2);

        if (includeintoscore) {
            if (lastresx[name1][name2].score)
                score += ctx->penalty;
            lastresx[name1][name2].score += score;
        }

        //		fprintf( stderr, "score(%d,%d) = %d\n", name1, name2, lastresx[name1][name2].score );
    }
    free(gett);
}

static void*
lastcallthread(lastcallthread_arg_t* targ) {
    Context*   ctx = targ->ctx;
    aln_Opts   opts = targ->opts;
    int        k, i;
    int        nq = targ->nq;
    int        nd = targ->nd;
    Lastresx** lastresx = targ->lastresx;
    char**     dseq = targ->dseq;
    char**     qseq = targ->qseq;
    char       command[5000];
    FILE*      lfp;
    int        msize;
    int        klim;
    char       kd[1000];

    k = -1;
    while (1) {
        {
            k++;
            if (k == nq)
                break;
            fprintf(stderr, "\r%d / %d                    \r", k, nq);
        }

        if (opts.alg == 'R') {
            klim = MIN(k, ctx->njob - ctx->nadd);
            if (klim == k) {
                sprintf(command, "_db%dd", k);
                lfp = fopen(command, "w");
                if (!lfp) {
                    fprintf(stderr, "Cannot open _db.");
                    exit(1);
                }
                for (i = 0; i < klim; i++)
                    fprintf(lfp, ">%d\n%s\n", i, dseq[i]);
                fclose(lfp);

                //				sprintf( command, "md5sum _db%dd > /dev/tty", k );
                //				system( command );

                if (ctx->dorp == 'd')
                    sprintf(command, "%s/lastdb _db%dd _db%dd", whereispairalign, k, k);
                else
                    sprintf(command, "%s/lastdb -p _db%dd _db%dd", whereispairalign, k, k);
                system(command);
                sprintf(kd, "%d", k);
            } else {
                kd[0] = 0;
            }
        } else {
            kd[0] = 0;
        }

        sprintf(command, "_q%d", k);
        lfp = fopen(command, "w");
        if (!lfp) {
            fprintf(stderr, "Cannot open %s", command);
            exit(1);
        }
        fprintf(lfp, ">%d\n%s\n", k, qseq[k]);
        fclose(lfp);

        //		if( alg == 'R' ) msize = MAX(10,k+nq);
        //			else msize = MAX(10,nd+nq);
        if (opts.alg == 'R')
            msize = MAX(10, k * lastm);
        else
            msize = MAX(10, nd * lastm);

        //		fprintf( stderr, "Calling lastal from lastcallthread, msize = %d, k=%d\n", msize, k );
        //		sprintf( command, "grep '>' _db%sd", kd );
        //		system( command );
        sprintf(command, "%s/lastal -m %d -e %d -f 0 -s 1 -p _scoringmatrixforlast -a %d -b %d _db%sd _q%d > _lastres%d", whereispairalign, msize, laste, -ctx->penalty, -ctx->penalty_ex, kd, k, k);
        if (system(command))
            exit(1);

        sprintf(command, "_lastres%d", k);
        lfp = fopen(command, "r");
        if (!lfp) {
            fprintf(stderr, "Cannot read _lastres%d", k);
            exit(1);
        }
        readlastresx_singleq(opts, ctx, lfp, k, lastresx);
        fclose(lfp);
    }
    return (NULL);
}

static void
calllast_fast(aln_Opts opts, Context* ctx, int nd, char** dseq, int nq, char** qseq, Lastresx** lastresx) {
    int   i, j;
    FILE* lfp;
    char  command[1000];

    lfp = fopen("_scoringmatrixforlast", "w");
    if (!lfp) {
        fprintf(stderr, "Cannot open _scoringmatrixforlast");
        exit(1);
    }
    if (ctx->dorp == 'd') {
        fprintf(lfp, "      ");
        for (j = 0; j < 4; j++)
            fprintf(lfp, " %c ", ctx->amino[j]);
        fprintf(lfp, "\n");
        for (i = 0; i < 4; i++) {
            fprintf(lfp, "%c ", ctx->amino[i]);
            for (j = 0; j < 4; j++)
                fprintf(lfp, " %d ", ctx->n_dis[i][j]);
            fprintf(lfp, "\n");
        }
    } else {
        fprintf(lfp, "      ");
        for (j = 0; j < 20; j++)
            fprintf(lfp, " %c ", ctx->amino[j]);
        fprintf(lfp, "\n");
        for (i = 0; i < 20; i++) {
            fprintf(lfp, "%c ", ctx->amino[i]);
            for (j = 0; j < 20; j++)
                fprintf(lfp, " %d ", ctx->n_dis[i][j]);
            fprintf(lfp, "\n");
        }
    }
    fclose(lfp);

    //	if( alg == 'r' ) // if 'R' -> lastcallthread, kokonoha nadd>0 no toki nomi shiyou
    {
        sprintf(command, "_dbd");
        lfp = fopen(command, "w");
        if (!lfp) {
            fprintf(stderr, "Cannot open _dbd");
            exit(1);
        }
        if (opts.alg == 'R')
            j = ctx->njob - ctx->nadd;
        else
            j = nd;
        for (i = 0; i < j; i++)
            fprintf(lfp, ">%d\n%s\n", i, dseq[i]);

        fclose(lfp);
        if (ctx->dorp == 'd')
            sprintf(command, "%s/lastdb _dbd _dbd", whereispairalign);
        else
            sprintf(command, "%s/lastdb -p _dbd _dbd", whereispairalign);
        system(command);
    }

    {
        lastcallthread_arg_t targ = {
            .ctx = ctx,
            .nq = nq,
            .nd = nd,
            .dseq = dseq,
            .qseq = qseq,
            .lastresx = lastresx,
        };
        lastcallthread(&targ);
    }
}

static void
calllast_once(aln_Opts opts, Context* ctx, int nd, char** dseq, int nq, char** qseq, Lastresx** lastresx) {
    int   i, j;
    char  command[5000];
    FILE* lfp;
    int   msize;
    int   res;

    fprintf(stderr, "nq=%d\n", nq);

    lfp = fopen("_db", "w");
    if (!lfp) {
        fprintf(stderr, "Cannot open _db");
        exit(1);
    }
    for (i = 0; i < nd; i++)
        fprintf(lfp, ">%d\n%s\n", i, dseq[i]);
    fclose(lfp);

    if (ctx->dorp == 'd') {
        sprintf(command, "%s/lastdb _db _db", whereispairalign);
        system(command);
        lfp = fopen("_scoringmatrixforlast", "w");
        if (!lfp) {
            fprintf(stderr, "Cannot open _scoringmatrixforlast");
            exit(1);
        }
        fprintf(lfp, "      ");
        for (j = 0; j < 4; j++)
            fprintf(lfp, " %c ", ctx->amino[j]);
        fprintf(lfp, "\n");
        for (i = 0; i < 4; i++) {
            fprintf(lfp, "%c ", ctx->amino[i]);
            for (j = 0; j < 4; j++)
                fprintf(lfp, " %d ", ctx->n_dis[i][j]);
            fprintf(lfp, "\n");
        }
        fclose(lfp);
    } else {
        sprintf(command, "%s/lastdb -p _db _db", whereispairalign);
        system(command);
        lfp = fopen("_scoringmatrixforlast", "w");
        if (!lfp) {
            fprintf(stderr, "Cannot open _scoringmatrixforlast");
            exit(1);
        }
        fprintf(lfp, "      ");
        for (j = 0; j < 20; j++)
            fprintf(lfp, " %c ", ctx->amino[j]);
        fprintf(lfp, "\n");
        for (i = 0; i < 20; i++) {
            fprintf(lfp, "%c ", ctx->amino[i]);
            for (j = 0; j < 20; j++)
                fprintf(lfp, " %d ", ctx->n_dis[i][j]);
            fprintf(lfp, "\n");
        }
        fclose(lfp);
        //		fprintf( stderr, "Not written yet\n" );
    }

    lfp = fopen("_q", "w");
    if (!lfp) {
        fprintf(stderr, "Cannot open _q");
        exit(1);
    }
    for (i = 0; i < nq; i++) {
        fprintf(lfp, ">%d\n%s\n", i, qseq[i]);
    }
    fclose(lfp);

    msize = MAX(10, nd * lastm);

    sprintf(command, "%s/lastal -v -m %d -e %d -f 0 -s 1 -p _scoringmatrixforlast -a %d -b %d _db _q > _lastres", whereispairalign, msize, laste, -ctx->penalty, -ctx->penalty_ex);
    res = system(command);
    if (res) {
        fprintf(stderr, "LAST aborted\n");
        exit(1);
    }

    lfp = fopen("_lastres", "r");
    if (!lfp) {
        fprintf(stderr, "Cannot read _lastres");
        exit(1);
    }
    //	readlastres( lfp, nd, nq, lastres, dseq, qseq );
    fprintf(stderr, "Reading lastres\n");
    readlastresx(opts, ctx, lfp, lastresx);
    fclose(lfp);
}

static void
callfoldalign(int nseq, char** mseq) {
    FILE*       fp;
    int         i;
    int         res;
    static char com[10000];

    for (i = 0; i < nseq; i++)
        t2u(mseq[i]);

    fp = fopen("_foldalignin", "w");
    if (!fp) {
        fprintf(stderr, "Cannot open _foldalignin\n");
        exit(1);
    }
    for (i = 0; i < nseq; i++) {
        fprintf(fp, ">%d\n", i + 1);
        fprintf(fp, "%s\n", mseq[i]);
    }
    fclose(fp);

    sprintf(com, "env PATH=%s  foldalign210 %s _foldalignin > _foldalignout ", whereispairalign, foldalignopt);
    res = system(com);
    if (res) {
        fprintf(stderr, "Error in foldalign\n");
        exit(1);
    }
}

static void
calllara(int nseq, char** mseq, char* laraarg) {
    FILE*       fp;
    int         i;
    int         res;
    static char com[10000];

    //	for( i=0; i<nseq; i++ )

    fp = fopen("_larain", "w");
    if (!fp) {
        fprintf(stderr, "Cannot open _larain\n");
        exit(1);
    }
    for (i = 0; i < nseq; i++) {
        fprintf(fp, ">%d\n", i + 1);
        fprintf(fp, "%s\n", mseq[i]);
    }
    fclose(fp);

    //	fprintf( stderr, "calling LaRA\n" );
    sprintf(com, "env PATH=%s:/bin:/usr/bin mafft_lara -i _larain -w _laraout -o _lara.params %s", whereispairalign, laraarg);
    res = system(com);
    if (res) {
        fprintf(stderr, "Error in lara\n");
        exit(1);
    }
}

static double
recalllara(Context* ctx, char** mseq1, char** mseq2, int alloclen) {
    static FILE* fp = NULL;
    static char* ungap1;
    static char* ungap2;
    static char* ori1;
    static char* ori2;
    //	int res;
    static char com[10000];
    double      value;

    if (fp == NULL) {
        fp = fopen("_laraout", "r");
        if (fp == NULL) {
            fprintf(stderr, "Cannot open _laraout\n");
            exit(1);
        }
        ungap1 = AllocateCharVec(alloclen);
        ungap2 = AllocateCharVec(alloclen);
        ori1 = AllocateCharVec(alloclen);
        ori2 = AllocateCharVec(alloclen);
    }

    strcpy(ori1, *mseq1);
    strcpy(ori2, *mseq2);

    fgets(com, 999, fp);
    myfgets(com, 9999, fp);
    strcpy(*mseq1, com);
    myfgets(com, 9999, fp);
    strcpy(*mseq2, com);

    copyWithNoGaps(ungap1, *mseq1);
    copyWithNoGaps(ungap2, *mseq2);
    t2u(ungap1);
    t2u(ungap2);
    t2u(ori1);
    t2u(ori2);

    if (strcmp(ungap1, ori1) || strcmp(ungap2, ori2)) {
        fprintf(stderr, "SEQUENCE CHANGED!!\n");
        fprintf(stderr, "*mseq1  = %s\n", *mseq1);
        fprintf(stderr, "ungap1  = %s\n", ungap1);
        fprintf(stderr, "ori1    = %s\n", ori1);
        fprintf(stderr, "*mseq2  = %s\n", *mseq2);
        fprintf(stderr, "ungap2  = %s\n", ungap2);
        fprintf(stderr, "ori2    = %s\n", ori2);
        exit(1);
    }

    value = (double)naivepairscore11(ctx, *mseq1, *mseq2, ctx->penalty);

    //	fclose( fp ); // saigo dake yatta houga yoi.

    return (value);
}

static double
calldafs_giving_bpp(Context* ctx, char** mseq1, char** mseq2, char** bpp1, char** bpp2, int i, int j) {
    FILE*  fp;
    int    res;
    char*  com;
    double value;
    char*  dirname;

    dirname = calloc(100, sizeof(char));
    com = calloc(1000, sizeof(char));
    sprintf(dirname, "_%d-%d", i, j);
    sprintf(com, "rm -rf %s", dirname);
    system(com);
    sprintf(com, "mkdir %s", dirname);
    system(com);

    sprintf(com, "%s/_bpporg", dirname);
    fp = fopen(com, "w");
    if (!fp) {
        fprintf(stderr, "Cannot write to %s/_bpporg\n", dirname);
        exit(1);
    }
    fprintf(fp, ">a\n");
    while (*bpp1)
        fprintf(fp, "%s", *bpp1++);

    fprintf(fp, ">b\n");
    while (*bpp2)
        fprintf(fp, "%s", *bpp2++);
    fclose(fp);

    sprintf(com, "tr -d '\\r' < %s/_bpporg > %s/_bpp", dirname, dirname);
    system(com);  // for cygwin, wakaran

    t2u(*mseq1);
    t2u(*mseq2);

    sprintf(com, "%s/_dafsinorg", dirname);
    fp = fopen(com, "w");
    if (!fp) {
        fprintf(stderr, "Cannot open %s/_dafsinorg\n", dirname);
        exit(1);
    }
    fprintf(fp, ">1\n");
    //	fprintf( fp, "%s\n", *mseq1 );
    write1seq(fp, *mseq1);
    fprintf(fp, ">2\n");
    //	fprintf( fp, "%s\n", *mseq2 );
    write1seq(fp, *mseq2);
    fclose(fp);

    sprintf(com, "tr -d '\\r' < %s/_dafsinorg > %s/_dafsin", dirname, dirname);
    system(com);  // for cygwin, wakaran

    sprintf(com, "_dafssh%s", dirname);
    fp = fopen(com, "w");
    fprintf(fp, "cd %s\n", dirname);
    fprintf(fp, "%s/dafs --mafft-in _bpp _dafsin > _dafsout 2>_dum\n", whereispairalign);
    fprintf(fp, "exit $tatus\n");
    fclose(fp);

    sprintf(com, "tr -d '\\r' < _dafssh%s > _dafssh%s.unix", dirname, dirname);
    system(com);  // for cygwin, wakaran

    sprintf(com, "sh _dafssh%s.unix 2>_dum%s", dirname, dirname);
    res = system(com);
    if (res) {
        fprintf(stderr, "Error in dafs\n");
        exit(1);
    }

    sprintf(com, "%s/_dafsout", dirname);

    fp = fopen(com, "r");
    if (!fp) {
        fprintf(stderr, "Cannot open %s/_dafsout\n", dirname);
        exit(1);
    }

    myfgets(com, 999, fp);  // nagai kanousei ga arunode
    fgets(com, 999, fp);
    myfgets(com, 999, fp);  // nagai kanousei ga arunode
    fgets(com, 999, fp);
    load1SeqWithoutName_new(ctx, fp, *mseq1);
    fgets(com, 999, fp);
    load1SeqWithoutName_new(ctx, fp, *mseq2);

    fclose(fp);

    //	fprintf( stderr, "*mseq1 = %s\n", *mseq1 );
    //	fprintf( stderr, "*mseq2 = %s\n", *mseq2 );

    value = (double)naivepairscore11(ctx, *mseq1, *mseq2, ctx->penalty);

#if 0
	sprintf( com, "rm -rf %s > /dev/null 2>&1", dirname );
	if( system( com ) )
	{
		fprintf( stderr, "retrying to rmdir\n" );
		usleep( 2000 );
		system( com );
	}
#endif

    free(dirname);
    free(com);

    return (value);
}

static double
callmxscarna_giving_bpp(Context* ctx, char** mseq1, char** mseq2, char** bpp1, char** bpp2, int i, int j) {
    FILE*  fp;
    int    res;
    char*  com;
    double value;
    char*  dirname;

    dirname = calloc(100, sizeof(char));
    com = calloc(1000, sizeof(char));
    sprintf(dirname, "_%d-%d", i, j);
    sprintf(com, "rm -rf %s", dirname);
    system(com);
    sprintf(com, "mkdir %s", dirname);
    system(com);

    sprintf(com, "%s/_bpporg", dirname);
    fp = fopen(com, "w");
    if (!fp) {
        fprintf(stderr, "Cannot write to %s/_bpporg\n", dirname);
        exit(1);
    }
    fprintf(fp, ">a\n");
    while (*bpp1)
        fprintf(fp, "%s", *bpp1++);

    fprintf(fp, ">b\n");
    while (*bpp2)
        fprintf(fp, "%s", *bpp2++);
    fclose(fp);

    sprintf(com, "tr -d '\\r' < %s/_bpporg > %s/_bpp", dirname, dirname);
    system(com);  // for cygwin, wakaran

    t2u(*mseq1);
    t2u(*mseq2);

    sprintf(com, "%s/_mxscarnainorg", dirname);
    fp = fopen(com, "w");
    if (!fp) {
        fprintf(stderr, "Cannot open %s/_mxscarnainorg\n", dirname);
        exit(1);
    }
    fprintf(fp, ">1\n");
    //	fprintf( fp, "%s\n", *mseq1 );
    write1seq(fp, *mseq1);
    fprintf(fp, ">2\n");
    //	fprintf( fp, "%s\n", *mseq2 );
    write1seq(fp, *mseq2);
    fclose(fp);

    sprintf(com, "tr -d '\\r' < %s/_mxscarnainorg > %s/_mxscarnain", dirname, dirname);
    system(com);  // for cygwin, wakaran

#if 0
	sprintf( com, "cd %s; %s/mxscarnamod -readbpp _mxscarnain > _mxscarnaout 2>_dum", dirname, whereispairalign );
#else
    sprintf(com, "_mxscarnash%s", dirname);
    fp = fopen(com, "w");
    fprintf(fp, "cd %s\n", dirname);
    fprintf(fp, "%s/mxscarnamod -readbpp _mxscarnain > _mxscarnaout 2>_dum\n", whereispairalign);
    fprintf(fp, "exit $tatus\n");
    fclose(fp);
    //sleep( 10000 );

    sprintf(com, "tr -d '\\r' < _mxscarnash%s > _mxscarnash%s.unix", dirname, dirname);
    system(com);  // for cygwin, wakaran

    sprintf(com, "sh _mxscarnash%s.unix 2>_dum%s", dirname, dirname);
#endif
    res = system(com);
    if (res) {
        fprintf(stderr, "Error in mxscarna\n");
        exit(1);
    }

    sprintf(com, "%s/_mxscarnaout", dirname);

    fp = fopen(com, "r");
    if (!fp) {
        fprintf(stderr, "Cannot open %s/_mxscarnaout\n", dirname);
        exit(1);
    }

    fgets(com, 999, fp);
    load1SeqWithoutName_new(ctx, fp, *mseq1);
    fgets(com, 999, fp);
    load1SeqWithoutName_new(ctx, fp, *mseq2);

    fclose(fp);

    //	fprintf( stderr, "*mseq1 = %s\n", *mseq1 );
    //	fprintf( stderr, "*mseq2 = %s\n", *mseq2 );

    value = (double)naivepairscore11(ctx, *mseq1, *mseq2, ctx->penalty);

#if 0
	sprintf( com, "rm -rf %s > /dev/null 2>&1", dirname );
	if( system( com ) )
	{
		fprintf( stderr, "retrying to rmdir\n" );
		usleep( 2000 );
		system( com );
	}
#endif

    free(dirname);
    free(com);

    return (value);
}

static void
readhat4(FILE* fp, char*** bpp) {
    char oneline[1000];
    int  bppsize;
    int  onechar;
    //	double prob;
    //	int posi, posj;

    bppsize = 0;
    //	fprintf( stderr, "reading hat4\n" );
    onechar = getc(fp);
    //	fprintf( stderr, "onechar = %c\n", onechar );
    if (onechar != '>') {
        fprintf(stderr, "Format error\n");
        exit(1);
    }
    ungetc(onechar, fp);
    fgets(oneline, 999, fp);
    while (1) {
        onechar = getc(fp);
        ungetc(onechar, fp);
        if (onechar == '>' || onechar == EOF) {
            //			fprintf( stderr, "Next\n" );
            *bpp = realloc(*bpp, (bppsize + 2) * sizeof(char*));
            (*bpp)[bppsize] = NULL;
            break;
        }
        fgets(oneline, 999, fp);
        //		fprintf( stderr, "oneline=%s\n", oneline );
        //		sscanf( oneline, "%d %d %lf", &posi, &posj, &prob );
        //		fprintf( stderr, "%d %d -> %f\n", posi, posj, prob );
        *bpp = realloc(*bpp, (bppsize + 2) * sizeof(char*));
        (*bpp)[bppsize] = calloc(100, sizeof(char));
        strcpy((*bpp)[bppsize], oneline);
        bppsize++;
    }
}

static void
preparebpp(int nseq, char*** bpp) {
    FILE* fp;
    int   i;

    fp = fopen("hat4", "r");
    if (!fp) {
        fprintf(stderr, "Cannot open hat4\n");
        exit(1);
    }
    for (i = 0; i < nseq; i++)
        readhat4(fp, bpp + i);
    fclose(fp);
}

int
countamino(char* s, int end) {
    int val = 0;
    while (end--)
        if (*s++ != '-')
            val++;
    return (val);
}

static double
score2dist(double pscore, double selfscore1, double selfscore2) {
    double val;
    double bunbo;
    //	fprintf( stderr, "In score2dist\n" );

    if ((bunbo = MIN(selfscore1, selfscore2)) == 0.0)
        val = 2.0;
    else if (bunbo < pscore)  // mondai ari
        val = 0.0;
    else
        val = (1.0 - pscore / bunbo) * 2.0;
    return (val);
}

static void
pairalign(aln_Opts opts, Context* ctx, const char* const* name, char** seq, char** aseq, char** dseq, int* thereisxineachseq, char** mseq1, char** mseq2, int alloclen, Lastresx** lastresx, double** distancemtx, LocalHom** localhomtable, double** expdist, int ngui) {
    int    i, j, ilim, jst, jj;
    int    off1, off2, dum1, dum2, thereisx;
    double pscore = 0.0;  // by D.Mathog
    FILE * hat2p, *hat3p;
    //	double **distancemtx;
    double*     selfscore;
    double*     effarr1;
    double*     effarr2;
    const char* pt;
    char*       hat2file = "hat2";
    //	LocalHom **localhomtable = NULL,
    LocalHom* tmpptr;
    int       intdum;
    char***   bpp = NULL;  // mxscarna no toki dake
    char **   distseq1, **distseq2;
    char **   dumseq1, **dumseq2;
    double    dist;
    double    scoreoffset;
    int       ntarget;
    int *     targetmap, *targetmapr;

    {
        ntarget = ctx->njob;
        targetmap = calloc(ctx->njob, sizeof(int));
        targetmapr = calloc(ctx->njob, sizeof(int));
        for (i = 0; i < ctx->njob; i++)
            targetmap[i] = targetmapr[i] = i;
    }

    if (store_localhom && localhomtable == NULL) {
        if (opts.alg == 'Y' || opts.alg == 'r') {
            ilim = ctx->njob - ctx->nadd;
            jst = ctx->nadd;
        } else {
            ilim = ntarget;
            jst = ctx->njob;
        }
        localhomtable = (LocalHom**)calloc(ilim, sizeof(LocalHom*));
        for (i = 0; i < ilim; i++) {
            localhomtable[i] = (LocalHom*)calloc(jst, sizeof(LocalHom));
            for (j = 0; j < jst; j++) {
                localhomtable[i][j].start1 = -1;
                localhomtable[i][j].end1 = -1;
                localhomtable[i][j].start2 = -1;
                localhomtable[i][j].end2 = -1;
                localhomtable[i][j].opt = -1.0;
                localhomtable[i][j].next = NULL;
                localhomtable[i][j].nokori = 0;
                localhomtable[i][j].extended = -1;
                localhomtable[i][j].last = localhomtable[i] + j;
                localhomtable[i][j].korh = 'h';
            }
            if (opts.alg != 'Y' && opts.alg != 'r')
                jst--;
        }
    }

    if (store_dist) {
        if (ngui == 0) {
            if (opts.alg == 'Y' || opts.alg == 'r')
                distancemtx = AllocateDoubleMtx(ctx->njob - ctx->nadd, ctx->nadd);  // 2020/Oct/23
            else
                distancemtx = AllocateDoubleHalfMtx(ctx->njob);
        }
    } else
        distancemtx = NULL;

    if (opts.alg == 'N') {
        dumseq1 = AllocateCharMtx(1, alloclen + 10);
        dumseq2 = AllocateCharMtx(1, alloclen + 10);
    }
    distseq1 = AllocateCharMtx(1, 0);  // muda
    distseq2 = AllocateCharMtx(1, 0);  // muda

    selfscore = AllocateDoubleVec(ctx->njob);
    effarr1 = AllocateDoubleVec(ctx->njob);
    effarr2 = AllocateDoubleVec(ctx->njob);

    reporterr("All-to-all alignment.\n");
    if (opts.alg == 'R') {
        fprintf(stderr, "Calling last (http://last.cbrc.jp/)\n");
        if (lastonce)
            calllast_once(opts, ctx, ctx->njob, seq, ctx->njob, seq, lastresx);
        else
            calllast_fast(opts, ctx, ctx->njob, seq, ctx->njob, seq, lastresx);
        fprintf(stderr, "done.\n");
    }

    if (opts.alg == 'r') {
        fprintf(stderr, "Calling last (http://last.cbrc.jp/)\n");
        fprintf(stderr, "nadd=%d\n", ctx->nadd);
        if (lastonce)
            calllast_once(opts, ctx, ctx->njob - ctx->nadd, seq, ctx->nadd, seq + ctx->njob - ctx->nadd, lastresx);
        else
            calllast_fast(opts, ctx, ctx->njob - ctx->nadd, seq, ctx->nadd, seq + ctx->njob - ctx->nadd, lastresx);

        fprintf(stderr, "nadd=%d\n", ctx->nadd);
        fprintf(stderr, "done.\n");
    }

    if (opts.alg == 'H') {
        fprintf(stderr, "Calling FOLDALIGN with option '%s'\n", foldalignopt);
        callfoldalign(ctx->njob, seq);
        fprintf(stderr, "done.\n");
    }
    if (opts.alg == 'B') {
        fprintf(stderr, "Running LARA (Bauer et al. http://www.planet-lisa.net/)\n");
        calllara(ctx->njob, seq, "");
        fprintf(stderr, "done.\n");
    }
    if (opts.alg == 'T') {
        fprintf(stderr, "Running SLARA (Bauer et al. http://www.planet-lisa.net/)\n");
        calllara(ctx->njob, seq, "-s");
        fprintf(stderr, "done.\n");
    }
    if (opts.alg == 's') {
        fprintf(stderr, "Preparing bpp\n");
        //		bpp = AllocateCharCub( ctx->njob, nlenmax, 0 );
        bpp = calloc(ctx->njob, sizeof(char**));
        preparebpp(ctx->njob, bpp);
        fprintf(stderr, "done.\n");
        fprintf(stderr, "Running MXSCARNA (Tabei et al. http://www.ncrna.org/software/mxscarna)\n");
    }
    if (opts.alg == 'G') {
        fprintf(stderr, "Preparing bpp\n");
        //		bpp = AllocateCharCub( ctx->njob, nlenmax, 0 );
        bpp = calloc(ctx->njob, sizeof(char**));
        preparebpp(ctx->njob, bpp);
        fprintf(stderr, "done.\n");
        fprintf(stderr, "Running DAFS (Sato et al. http://www.ncrna.org/)\n");
    }

    for (i = 0; i < ctx->njob; i++) {
        pscore = 0.0;
        for (pt = seq[i]; *pt; pt++)
            pscore += ctx->amino_dis[(unsigned char)*pt][(unsigned char)*pt];
        selfscore[i] = pscore;
    }

    {
        double** dynamicmtx = NULL;
        if (opts.specificityconsideration > 0.0)
            dynamicmtx = AllocateDoubleMtx(ctx->nalphabets, ctx->nalphabets);

        if (opts.alg == 'Y' || opts.alg == 'r')
            ilim = ctx->njob - ctx->nadd;
        else
            ilim = ctx->njob - 1;
        for (i = 0; i < ilim; i++) {
            if (stdout_dist)
                fprintf(stdout, "%d %d d=%.3f\n", i + 1, i + 1, 0.0);
            fprintf(stderr, "% 5d / %d\r", i, ctx->njob - ctx->nadd);
            fflush(stderr);

            if (opts.alg == 'Y' || opts.alg == 'r')
                jst = ctx->njob - ctx->nadd;
            else
                jst = i + 1;
            for (j = jst; j < ctx->njob; j++) {
                if (strlen(seq[i]) == 0 || strlen(seq[j]) == 0) {
                    if (store_dist) {
                        if (opts.alg == 'Y' || opts.alg == 'r')
                            distancemtx[i][j - (ctx->njob - ctx->nadd)] = 3.0;
                        else
                            distancemtx[i][j - i] = 3.0;
                    }
                    if (stdout_dist)
                        fprintf(stdout, "%d %d d=%.3f\n", i + 1, j + 1, 3.0);
                    continue;
                }

                strcpy(aseq[0], seq[i]);
                strcpy(aseq[1], seq[j]);

                effarr1[0] = 1.0;
                effarr2[0] = 1.0;
                mseq1[0] = aseq[0];
                mseq2[0] = aseq[1];

                thereisx = thereisxineachseq[i] + thereisxineachseq[j];
                distseq1[0] = dseq[i];
                distseq2[0] = dseq[j];

                if (opts.use_fft) {
                    pscore = Falign(opts, ctx, NULL, NULL, ctx->n_dis_consweight_multi, mseq1, mseq2, effarr1, effarr2, NULL, NULL, 1, 1, alloclen, &intdum);
                    off1 = off2 = 0;
                } else {
                    switch (opts.alg) {
                        case ('t'):
                            pscore = G__align11_noalign(ctx, ctx->n_dis_consweight_multi, ctx->penalty, ctx->penalty_ex, distseq1, distseq2);
                            off1 = off2 = 0;
                            break;
                        case ('A'):
                            if (ctx->usenaivescoreinsteadofalignmentscore) {
                                G__align11(ctx, ctx->n_dis_consweight_multi, mseq1, mseq2, alloclen, ctx->outgap, ctx->outgap);
                                pscore = (double)naivepairscore11(ctx, mseq1[0], mseq2[0], 0.0);  // uwagaki
                            } else {
                                //								if( store_localhom )
                                if (store_localhom && (targetmap[i] != -1 || targetmap[j] != -1)) {
                                    pscore = G__align11(ctx, ctx->n_dis_consweight_multi, mseq1, mseq2, alloclen, ctx->outgap, ctx->outgap);
                                    if (thereisx)
                                        pscore = G__align11_noalign(ctx, ctx->n_dis_consweight_multi, ctx->penalty, ctx->penalty_ex, distseq1, distseq2);
                                    if (opts.specificityconsideration > 0.0) {
                                        if (expdist)
                                            dist = expdist[i][j];
                                        else
                                            dist = score2dist(pscore, selfscore[i], selfscore[j]);
                                        if (dist2offset(opts, dist) < 0.0) {
                                            makedynamicmtx(opts, ctx, dynamicmtx, ctx->n_dis_consweight_multi, 0.5 * dist);  // upgma ni awaseru.
                                            strcpy(mseq1[0], seq[i]);
                                            strcpy(mseq2[0], seq[j]);
                                            G__align11(ctx, dynamicmtx, mseq1, mseq2, alloclen, ctx->outgap, ctx->outgap);
                                        }
                                    }
                                } else
                                    pscore = G__align11_noalign(ctx, ctx->n_dis_consweight_multi, ctx->penalty, ctx->penalty_ex, distseq1, distseq2);
                            }
                            off1 = off2 = 0;
                            break;
                        case ('N'):
                            if (ctx->usenaivescoreinsteadofalignmentscore) {
                                genL__align11(ctx, ctx->n_dis_consweight_multi, mseq1, mseq2, alloclen, &off1, &off2);
                                pscore = (double)naivepairscore11(ctx, mseq1[0], mseq2[0], 0.0);  // uwagaki
                            } else {
                                pscore = genL__align11(ctx, ctx->n_dis_consweight_multi, mseq1, mseq2, alloclen, &off1, &off2);
                                if (thereisx) {
                                    strcpy(dumseq1[0], distseq1[0]);
                                    strcpy(dumseq2[0], distseq2[0]);
                                    pscore = genL__align11(ctx, ctx->n_dis_consweight_multi, dumseq1, dumseq2, alloclen, &dum1, &dum2);  // uwagaki
                                }
                                if (opts.specificityconsideration > 0.0) {
                                    //									fprintf( stderr, "dist = %f\n", score2dist( pscore, selfscore[i], selfscore[j] ) );
                                    if (expdist)
                                        dist = expdist[i][j];
                                    else
                                        dist = score2dist(pscore, selfscore[i], selfscore[j]);
                                    if (dist2offset(opts, dist) < 0.0) {
                                        makedynamicmtx(opts, ctx, dynamicmtx, ctx->n_dis_consweight_multi, 0.5 * dist);  // upgma ni awaseru.
                                        strcpy(mseq1[0], seq[i]);
                                        strcpy(mseq2[0], seq[j]);
                                        genL__align11(ctx, dynamicmtx, mseq1, mseq2, alloclen, &off1, &off2);
                                    }
                                }
                            }
                            break;
                        case ('R'):
                            if (ctx->nadd && ctx->njob - ctx->nadd <= j && ctx->njob - ctx->nadd <= i)  // new sequence doushi ha mushi
                                pscore = 0.0;
                            else
                                pscore = (double)lastresx[i][j].score;  // all pair
                            break;
                        case ('r'):
                            if (ctx->nadd == 0 || (i < ctx->njob - ctx->nadd && ctx->njob - ctx->nadd <= j))
                                pscore = (double)lastresx[i][j - (ctx->njob - ctx->nadd)].score;
                            else
                                pscore = 0.0;
                            break;
                        case ('L'):
                            if (ctx->nadd && ctx->njob - ctx->nadd <= j && ctx->njob - ctx->nadd <= i)  // new sequence doushi ha mushi
                                pscore = 0.0;
                            else {
                                if (ctx->usenaivescoreinsteadofalignmentscore) {
                                    L__align11(ctx, ctx->n_dis_consweight_multi, 0.0, mseq1, mseq2, alloclen, &off1, &off2);
                                    pscore = (double)naivepairscore11(ctx, mseq1[0], mseq2[0], 0.0);  // uwagaki
                                } else {
                                    //									if( store_localhom )
                                    if (store_localhom && (targetmap[i] != -1 || targetmap[j] != -1)) {
                                        pscore = L__align11(ctx, ctx->n_dis_consweight_multi, 0.0, mseq1, mseq2, alloclen, &off1, &off2);  // all pair
                                        if (thereisx)
                                            pscore = L__align11_noalign(ctx, ctx->n_dis_consweight_multi, distseq1, distseq2);  // all pair
                                        if (opts.specificityconsideration > 0.0) {
                                            if (expdist)
                                                dist = expdist[i][j];
                                            else
                                                dist = score2dist(pscore, selfscore[i], selfscore[j]);
                                            if ((scoreoffset = dist2offset(opts, dist)) < 0.0) {
                                                makedynamicmtx(opts, ctx, dynamicmtx, ctx->n_dis_consweight_multi, 0.5 * dist);  // upgma ni awaseru.
                                                strcpy(mseq1[0], seq[i]);
                                                strcpy(mseq2[0], seq[j]);
                                                L__align11(ctx, dynamicmtx, scoreoffset, mseq1, mseq2, alloclen, &off1, &off2);
                                            }
                                        }
                                    } else
                                        pscore = L__align11_noalign(ctx, ctx->n_dis_consweight_multi, distseq1, distseq2);  // all pair
                                }
                            }
                            break;
                        case ('Y'):
                            if (ctx->nadd == 0 || (i < ctx->njob - ctx->nadd && ctx->njob - ctx->nadd <= j))  // new sequence vs exiting sequence nomi keisan
                            {
                                if (ctx->usenaivescoreinsteadofalignmentscore) {
                                    L__align11(ctx, ctx->n_dis_consweight_multi, 0.0, mseq1, mseq2, alloclen, &off1, &off2);
                                    pscore = (double)naivepairscore11(ctx, mseq1[0], mseq2[0], 0.0);  // uwagaki
                                } else {
                                    if (store_localhom) {
                                        pscore = L__align11(ctx, ctx->n_dis_consweight_multi, 0.0, mseq1, mseq2, alloclen, &off1, &off2);
                                        if (thereisx)
                                            pscore = L__align11_noalign(ctx, ctx->n_dis_consweight_multi, distseq1, distseq2);  // uwagaki
                                    } else
                                        pscore = L__align11_noalign(ctx, ctx->n_dis_consweight_multi, distseq1, distseq2);
                                }
                            } else
                                pscore = 0.0;
                            break;
                        case ('a'):
                            pscore = Aalign(ctx, mseq1, mseq2, effarr1, effarr2, 1, 1, alloclen);
                            off1 = off2 = 0;
                            break;
#if 0
						case( 'K' ):
							pscore = genG__align11( mseq1, mseq2, alloclen );
							off1 = off2 = 0;
							break;
#endif
                        case ('H'):
                            pscore = recallpairfoldalign(ctx, mseq1, mseq2, i, j, &off1, &off2, alloclen);
                            break;
                        case ('B'):
                        case ('T'):
                            pscore = recalllara(ctx, mseq1, mseq2, alloclen);
                            off1 = off2 = 0;
                            break;
                        case ('s'):
                            pscore = callmxscarna_giving_bpp(ctx, mseq1, mseq2, bpp[i], bpp[j], i, j);
                            off1 = off2 = 0;
                            break;
                        case ('G'):
                            pscore = calldafs_giving_bpp(ctx, mseq1, mseq2, bpp[i], bpp[j], i, j);
                            off1 = off2 = 0;
                            break;
                        case ('M'):
                            pscore = MSalign11(ctx, mseq1, mseq2, alloclen);
                            break;
                        default:
                            ErrorExit("ERROR IN SOURCE FILE");
                    }
                }

                if (opts.alg == 't' || (mseq1[0][0] != 0 && mseq2[0][0] != 0))  // 't' no jouken ha iranai to omou. if( ( mseq1[0][0] != 0 && mseq2[0][0] != 0  ) )
                {
                    if ((ctx->nadd == 0 || (opts.alg != 'Y' && opts.alg != 'r') || (i < ctx->njob - ctx->nadd && ctx->njob - ctx->nadd <= j))) {
                        if (!store_localhom)
                            ;
                        else if (opts.alg == 'R')
                            putlocalhom_last(ctx, mseq1[0], mseq2[0], localhomtable[i] + j, lastresx[i] + j);
                        else if (opts.alg == 'r')
                            putlocalhom_last(ctx, mseq1[0], mseq2[0], localhomtable[i] + j - (ctx->njob - ctx->nadd), lastresx[i] + j - (ctx->njob - ctx->nadd));
                        else if (opts.alg == 'H')
                            putlocalhom_ext(ctx, mseq1[0], mseq2[0], localhomtable[i] + j, off1, off2, 'h');
                        else if (opts.alg == 'Y')
                            putlocalhom2(ctx, mseq1[0], mseq2[0], localhomtable[i] + j - (ctx->njob - ctx->nadd), off1, off2, 'h');
                        else if (opts.alg != 'S' && opts.alg != 'V')
                            putlocalhom2(ctx, mseq1[0], mseq2[0], localhomtable[i] + j - i, off1, off2, 'h');
                        else {
                            if (targetmap[i] != -1 && targetmap[j] != -1) {
                                putlocalhom2(ctx, mseq1[0], mseq2[0], localhomtable[targetmap[i]] + j, off1, off2, 'h');
                                putlocalhom2(ctx, mseq2[0], mseq1[0], localhomtable[targetmap[j]] + i, off2, off1, 'h');  // sukoshi muda.
                            } else if (targetmap[j] != -1)
                                putlocalhom2(ctx, mseq2[0], mseq1[0], localhomtable[targetmap[j]] + i, off2, off1, 'h');
                            else if (targetmap[i] != -1)
                                putlocalhom2(ctx, mseq1[0], mseq2[0], localhomtable[targetmap[i]] + j, off1, off2, 'h');
                            else {
                                reporterr("okashii\n");
                                exit(1);
                            }
                        }
                    }

                    pscore = score2dist(pscore, selfscore[i], selfscore[j]);
                } else {
                    pscore = 2.0;
                }

                if (stdout_align) {
                    if (opts.alg != 't') {
                        fprintf(stdout, "sequence %d - sequence %d, pairwise distance = %10.5f\n", i + 1, j + 1, pscore);
                        fprintf(stdout, ">%s\n", name[i]);
                        write1seq(stdout, mseq1[0]);
                        fprintf(stdout, ">%s\n", name[j]);
                        write1seq(stdout, mseq2[0]);
                        fprintf(stdout, "\n");
                    }
                }
                if (stdout_dist)
                    fprintf(stdout, "%d %d d=%.3f\n", i + 1, j + 1, pscore);
                if (store_dist) {
                    if (opts.alg == 'Y' || opts.alg == 'r')
                        distancemtx[i][j - (ctx->njob - ctx->nadd)] = pscore;
                    else
                        distancemtx[i][j - i] = pscore;
                }
            }
        }
        if (dynamicmtx)
            FreeDoubleMtx(dynamicmtx);
    }

    if (store_dist && ngui == 0) {
        hat2p = fopen(hat2file, "w");
        if (!hat2p)
            ErrorExit("Cannot open hat2.");
        if (opts.alg == 'Y' || opts.alg == 'r')
            WriteHat2_part_pointer(hat2p, ctx->njob, ctx->nadd, name, distancemtx);
        else
            WriteFloatHat2_pointer_halfmtx(ctx, hat2p, ctx->njob, name, distancemtx);  // jissiha double
        fclose(hat2p);
    }

    hat3p = fopen("hat3", "w");
    if (!hat3p)
        ErrorExit("Cannot open hat3.");
    if (store_localhom && ngui == 0) {
        fprintf(stderr, "\n\n##### writing hat3\n");
        if (opts.alg == 'Y' || opts.alg == 'r')
            ilim = ctx->njob - ctx->nadd;
        else
            ilim = ctx->njob - 1;
        for (i = 0; i < ilim; i++) {
            if (opts.alg == 'Y' || opts.alg == 'r') {
                jst = ctx->njob - ctx->nadd;
                jj = 0;
            } else {
                jst = i;
                jj = 0;
            }
            for (j = jst; j < ctx->njob; j++, jj++) {
                for (tmpptr = localhomtable[i] + jj; tmpptr; tmpptr = tmpptr->next) {
                    //					fprintf( stderr, "j=%d, jj=%d\n", j, jj );
                    if (tmpptr->opt == -1.0)
                        continue;
                    // tmptmptmptmptmp
                    //					if( alg == 'B' || alg == 'T' )
                    //						fprintf( hat3p, "%d %d %d %7.5f %d %d %d %d %p\n", i, j, tmpptr->overlapaa, 1.0, tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, (void *)tmpptr->next );
                    //					else
                    if (targetmap[j] == -1 || targetmap[i] < targetmap[j])
                        fprintf(hat3p, "%d %d %d %7.5f %d %d %d %d h\n", targetmapr[i], j, tmpptr->overlapaa, tmpptr->opt, tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2);
                    //						fprintf( hat3p, "%d %d %d %7.5f %d %d %d %d h\n", i, j, tmpptr->overlapaa, tmpptr->opt, tmpptr->start1, tmpptr->end1, tmpptr->start2+1, tmpptr->end2+1 ); // zettai dame!!!!
                }
            }
        }
//		if( ngui == 0 )
//		{
#if DEBUG
        fprintf(stderr, "calling FreeLocalHomTable\n");
#endif
        if (opts.alg == 'Y' || opts.alg == 'r')
            FreeLocalHomTable_part(localhomtable, (ctx->njob - ctx->nadd), ctx->nadd);
        else
            FreeLocalHomTable_half(localhomtable, ctx->njob);
#if DEBUG
        fprintf(stderr, "done. FreeLocalHomTable\n");
#endif
        //		}
    }
    fclose(hat3p);

    if (opts.alg == 's') {
        char** ptpt;
        for (i = 0; i < ctx->njob; i++) {
            ptpt = bpp[i];
            while (1) {
                if (*ptpt)
                    free(*ptpt);
                else
                    break;
                ptpt++;
            }
            free(bpp[i]);
        }
        free(bpp);
    }
    free(selfscore);
    free(effarr1);
    free(effarr2);
    if (opts.alg == 'N') {
        FreeCharMtx(dumseq1);
        FreeCharMtx(dumseq2);
    }
    free(distseq1);
    free(distseq2);
    if (store_dist && ngui == 0) {
        if (opts.alg == 'Y' || opts.alg == 'r')
            FreeDoubleMtx(distancemtx);  // 2020/Oct/23
        else
            FreeDoubleHalfMtx(distancemtx, ctx->njob);
    }

    free(targetmap);
    free(targetmapr);
}

int
pairlocalalign(
    aln_Opts           pairLocalAlignOpts,
    Context*           ctx,
    const char* const* name,
    char**             seq,
    double**           iscore,
    LocalHom**         localhomtable
) {
    aln_assert(ctx->njob >= 2);

    int i;

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
    ctx->usenaivescoreinsteadofalignmentscore = 0;
    ctx->nwildcard = 0;

    int alloclen = ctx->maxInputSeqLen * 2;

    char** aseq = AllocateCharMtx(2, alloclen + 10);
    char** bseq = AllocateCharMtx(ctx->njob, alloclen + 10);
    char** dseq = AllocateCharMtx(ctx->njob, alloclen + 10);
    char** mseq1 = AllocateCharMtx(ctx->njob, 0);
    char** mseq2 = AllocateCharMtx(ctx->njob, 0);
    int*   thereisxineachseq = AllocateIntVec(ctx->njob);

    Lastresx** lastresx = NULL;

    constants(pairLocalAlignOpts, ctx, ctx->njob, seq);
    initSignalSM(ctx);
    initFiles(ctx);

    // TODO(sen) Sequence verification?

    if (ctx->dorp == 'p' && pairLocalAlignOpts.scoremtx == 1 && pairLocalAlignOpts.nblosum > 0) {
        for (i = 0; i < ctx->njob; i++) {
            copyWithNoGaps(bseq[i], seq[i]);
            thereisxineachseq[i] = removex(dseq[i], bseq[i]);
        }
    } else {
        for (i = 0; i < ctx->njob; i++) {
            copyWithNoGaps(bseq[i], seq[i]);
            strcpy(dseq[i], bseq[i]);
            thereisxineachseq[i] = 0;
        }
    }

    pairalign(
        pairLocalAlignOpts,
        ctx,
        name,
        bseq,
        aseq,
        dseq,
        thereisxineachseq,
        mseq1,
        mseq2,
        alloclen,
        lastresx,
        iscore,
        localhomtable,
        0,
        ctx->njob
    );

    fprintf(ctx->trap_g, "done.\n");
    fclose(ctx->trap_g);
    fclose(ctx->prep_g);

    return 0;
}
