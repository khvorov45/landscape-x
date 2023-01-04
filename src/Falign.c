#include "mltaln.h"

#if 0
static FILE *fftfp;
#endif
static int n20or4or2;

#define KEIKA 0
#define RND 0
#define DEBUG 0

#if RND  // by D.Mathog
static void
generateRndSeq(char* seq, int len) {
    while (len--)
#if 1
        *seq++ = (int)(rnd() * n20or4or2);
#else
        *seq++ = (int)1;
#endif
}
#endif

static void
vec_init(Fukusosuu* result, int nlen) {
    while (nlen--) {
        result->R = result->I = 0.0;
        result++;
    }
}

#if 0  // by D.Mathog
static void vec_init2( Fukusosuu **result, char *seq, double eff, int st, int ed )
{
	int i;
	for( i=st; i<ed; i++ )
		result[(int)*seq++][i].R += eff;
}
#endif

static void
seq_vec_3(Context* ctx, Fukusosuu** result, double incr, char* seq) {
    int i;
    int n;
    for (i = 0; *seq; i++) {
        n = ctx->amino_n[(int)*seq++];
        if (n < n20or4or2 && n >= 0)
            result[n][i].R += incr;
    }
}

static void
seq_vec_5(Context* ctx, Fukusosuu* result, double* score1, double* score2, double incr, char* seq) {
    int n;
    for (; *seq; result++) {
        n = ctx->amino_n[(int)*seq++];
        if (n > 20)
            continue;
        result->R += incr * score1[n];
        result->I += incr * score2[n];
#if 0
		fprintf( stderr, "n=%d, score=%f, inc=%f R=%f\n",n,  score[n], incr * score[n], result->R );
#endif
    }
}

static void
seq_vec_4(Fukusosuu* result, double incr, char* seq) {
    char s;
    for (; *seq; result++) {
        s = *seq++;
        if (s == 'a')
            result->R += incr;
        else if (s == 't')
            result->R -= incr;
        else if (s == 'g')
            result->I += incr;
        else if (s == 'c')
            result->I -= incr;
    }
}

#if 0  // by D.Mathog
static void seq_vec( Fukusosuu *result, char query, double incr, char *seq )
{
#if 0
	int bk = nlen;
#endif
	while( *seq )
	{
		if( *seq++ == query ) result->R += incr;
		result++;
#if 0
fprintf( stderr, "i = %d result->R = %f\n", bk-nlen, (result-1)->R );
#endif
	}
}

static int checkRepeat( int num, int *cutpos )
{
	int tmp, buf;

	buf = *cutpos;
	while( num-- )
	{
		if( ( tmp = *cutpos++ ) < buf ) return( 1 );
		buf = tmp;
	}
	return( 0 );
}

static int segcmp( void *ptr1, void *ptr2 )
{
	int diff;
	Segment **seg1 = (Segment **)ptr1;
	Segment **seg2 = (Segment **)ptr2;
#if 0
	return( (*seg1)->center - (*seg2)->center );
#else
	diff = (*seg1)->center - (*seg2)->center;
	if( diff ) return( diff );

	diff = (*seg1)->start - (*seg2)->start;
	if( diff ) return( diff );

	diff = (*seg1)->end - (*seg2)->end;
	if( diff ) return( diff );

	fprintf( stderr, "USE STABLE SORT !!\n" );
	exit( 1 );
	return( 0 );
#endif
}
#endif

static void
mymergesort(int first, int last, Segment** seg) {
    int              middle;
    static int       i, j, k, p;
    static int       allo = 0;
    static Segment** work = NULL;

    if (seg == NULL) {
        if (work)
            free(work);
        work = NULL;
        allo = 0;
        return;
    }

    if (last > allo) {
        allo = last;
        if (work)
            free(work);
        work = (Segment**)calloc(allo / 2 + 1, sizeof(Segment*));
    }

    if (first < last) {
        middle = (first + last) / 2;
        mymergesort(first, middle, seg);
        mymergesort(middle + 1, last, seg);
        p = 0;
        for (i = first; i <= middle; i++)
            work[p++] = seg[i];
        i = middle + 1;
        j = 0;
        k = first;
        while (i <= last && j < p) {
            if (work[j]->center <= seg[i]->center)
                seg[k++] = work[j++];
            else
                seg[k++] = seg[i++];
        }
        while (j < p)
            seg[k++] = work[j++];
    }
}

double
Falign(Context* ctx, int** whichmtx, double*** scoringmatrices, double** n_dynamicmtx, char** seq1, char** seq2, double* eff1, double* eff2, double** eff1s, double** eff2s, int clus1, int clus2, int alloclen, int* fftlog) {
    int        i, j, k, l, m, maxk;
    int        nlen, nlen2;
    static int crossscoresize = 0;
    char**     tmpseq1 = NULL;
    char**     tmpseq2 = NULL;
    char**     tmpptr1 = NULL;
    char**     tmpptr2 = NULL;
    char**     tmpres1 = NULL;
    char**     tmpres2 = NULL;
    char**     result1 = NULL;
    char**     result2 = NULL;
#if RND
    char** rndseq1 = NULL;
    char** rndseq2 = NULL;
#endif
    static Fukusosuu** seqVector1 = NULL;
    static Fukusosuu** seqVector2 = NULL;
    static Fukusosuu** naiseki = NULL;
    static Fukusosuu*  naisekiNoWa = NULL;
    static double*     soukan = NULL;
    static double**    crossscore = NULL;
    int                nlentmp;
    static int*        kouho = NULL;
    static Segment*    segment = NULL;
    static Segment*    segment1 = NULL;
    static Segment*    segment2 = NULL;
    static Segment**   sortedseg1 = NULL;
    static Segment**   sortedseg2 = NULL;
    static int*        cut1 = NULL;
    static int*        cut2 = NULL;
    char *             sgap1, *egap1, *sgap2, *egap2;
    static int         localalloclen = 0;
    int                lag;
    int                tmpint;
    int                count, count0;
    int                len1, len2;
    int                totallen;
    double             totalscore;
    double             dumdb = 0.0;
    int                headgp, tailgp;
    static double*     gstart = NULL;
    static double*     gend = NULL;
    static double**    codonscoremtx;

    if (seq1 == NULL) {
        if (kouho) {
            //			fprintf( stderr, "Freeing localarrays in Falign\n" );
            localalloclen = 0;
            crossscoresize = 0;
            mymergesort(0, 0, NULL);
            alignableReagion(ctx, 0, 0, NULL, NULL, NULL, NULL, NULL);
            fft(0, NULL, 1);
            A__align(ctx, NULL, 0, 0, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, 0, 0, -1, -1, NULL, NULL, NULL, 0.0, 0.0);
            D__align(ctx, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, 0, 0);
            A__align_variousdist(ctx, NULL, NULL, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, 0, 0);
            D__align_variousdist(ctx, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, 0, 0);
            G__align11(ctx, NULL, NULL, NULL, 0, 0, 0);
            G__align11psg(ctx, NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL);
            blockAlign2(NULL, NULL, NULL, NULL, NULL, NULL);
            if (crossscore)
                FreeDoubleMtx(crossscore);
            crossscore = NULL;
            free(kouho);
            kouho = NULL;
            free(cut1);
            free(cut2);
            free(segment);
            free(segment1);
            free(segment2);
            free(sortedseg1);
            free(sortedseg2);
            if (gstart)
                free(gstart);
            if (gend)
                free(gend);
            if (codonscoremtx)
                FreeDoubleMtx(codonscoremtx);
            if (!kobetsubunkatsu) {
                FreeFukusosuuMtx(seqVector1);
                FreeFukusosuuMtx(seqVector2);
                FreeFukusosuuVec(naisekiNoWa);
                FreeFukusosuuMtx(naiseki);
                FreeDoubleVec(soukan);
            }
        } else {
            //			fprintf( stderr, "Did not allocate localarrays in Falign\n" );
        }

        return (0.0);
    }

    len1 = strlen(seq1[0]);
    len2 = strlen(seq2[0]);
    nlentmp = MAX(len1, len2);

    nlen = 1;
    while (nlentmp >= nlen)
        nlen <<= 1;

    nlen2 = nlen / 2;

#if DEBUG
    fprintf(stderr, "len1 = %d, len2 = %d\n", len1, len2);
    fprintf(stderr, "nlentmp = %d, nlen = %d\n", nlentmp, nlen);
#endif

    result1 = AllocateCharMtx(clus1, alloclen);
    result2 = AllocateCharMtx(clus2, alloclen);
    tmpres1 = AllocateCharMtx(clus1, alloclen);
    tmpres2 = AllocateCharMtx(clus2, alloclen);
    sgap1 = AllocateCharVec(clus1);
    egap1 = AllocateCharVec(clus1);
    sgap2 = AllocateCharVec(clus2);
    egap2 = AllocateCharVec(clus2);
    tmpptr1 = calloc(clus1, sizeof(char*));
    tmpptr2 = calloc(clus2, sizeof(char*));
    tmpseq1 = AllocateCharMtx(clus1, nlen);
    tmpseq2 = AllocateCharMtx(clus2, nlen);
#if RND
    rndseq1 = AllocateCharMtx(clus1, nlen);
    rndseq2 = AllocateCharMtx(clus2, nlen);
    for (i = 0; i < clus1; i++)
        generateRndSeq(rndseq1[i], nlen);
    for (i = 0; i < clus2; i++)
        generateRndSeq(rndseq2[i], nlen);
#endif

    if (!localalloclen) {
        kouho = AllocateIntVec(NKOUHO);
        cut1 = AllocateIntVec(MAXSEG);
        cut2 = AllocateIntVec(MAXSEG);
        //		crossscore = AllocateDoubleMtx( MAXSEG, MAXSEG );
        segment = (Segment*)calloc(MAXSEG, sizeof(Segment));
        segment1 = (Segment*)calloc(MAXSEG, sizeof(Segment));
        segment2 = (Segment*)calloc(MAXSEG, sizeof(Segment));
        sortedseg1 = (Segment**)calloc(MAXSEG, sizeof(Segment*));
        sortedseg2 = (Segment**)calloc(MAXSEG, sizeof(Segment*));
        if (!(segment && segment1 && segment2 && sortedseg1 && sortedseg2))
            ErrorExit("Allocation error\n");

        if (scoremtx == -1)
            n20or4or2 = 1;
        else if (fftscore)
            n20or4or2 = 1;
        else
            n20or4or2 = 20;

        gstart = NULL;
        gend = NULL;
        if (codonpos) {
            FILE* cfp;
            char* buf = calloc(sizeof(char), 1000);

            if (dorp != 'd') {
                reporterr("\n\nThe --codonpos and --codonscore options are only for DNA data.\n\n");
                exit(1);
            }

            //			reporterr( "\nIn Falign, loading position-specific gap costs\n" );
            gstart = calloc(sizeof(double), len1 + 1);
            for (i = 0; i < len1 + 1; i++)
                gstart[i] = 1.0 * 0.5;  // init
            gend = calloc(sizeof(double), len1 + 1);
            for (i = 0; i < len1 + 1; i++)
                gend[i] = 1.0 * 0.5;  // init
#define STRONG 3.0

            i = 0;
            cfp = fopen("_codonpos", "r");
            if (cfp == NULL) {
                reporterr("Cannot open _codonpos file\n");
                exit(1);
            }
            while (fgets(buf, 1000, cfp)) {
                if (i == len1) {
                    reporterr("\n\nNumber of lines in the codonposition file must be the same as the length of the first sequences (%d).\n\n", len1);
                    exit(1);
                }

                if (buf[0] == '#') {
                    continue;
                } else if (buf[0] == '0') {
                    gstart[i] = 1.0 * 0.5;  // noncoding
                    gend[i] = 1.0 * 0.5;  // noncoding
                } else if (buf[0] == '1') {
                    gstart[i] = STRONG * 0.5;
                    gend[i] = 1.0 * 0.5;
                } else if (buf[0] == '2') {
                    gstart[i] = STRONG * 0.5;
                    gend[i] = STRONG * 0.5;
                } else if (buf[0] == '3') {
                    gstart[i] = 1.0 * 0.5;
                    gend[i] = STRONG * 0.5;
                } else {
                    reporterr("In the codonposition file, 1st letter in a line must be either of:\n");
                    reporterr("	0 (noncoding)\n");
                    reporterr("	1 (1st position in codon)\n");
                    reporterr("	2 (2nd position in codon)\n");
                    reporterr("	3 (3rd position in codon)\n");
                    reporterr("When mutliple difference frames are used, set 2.\n");
                    exit(1);
                }
                i++;
            }
            fclose(cfp);
            free(buf);
            if (i < len1) {
                reporterr("\n\nNumber of lines in the codonposition file (%d) is less than the length of the first sequences (%d).\n\n", i, len1);
                exit(1);
            }
        }

        if (codonscore) {
            int     i, j;
            FILE*   cfp;
            double* codonfreq = calloc(sizeof(double), 64);
            double  totalcount, codonscore0av, codonscore1av;
            if (!codonpos) {
                reporterr("\n\n --codonpos is necessary for --codonscore\n\n");
                exit(1);
            }
            codonscoremtx = AllocateDoubleMtx(64, 64);
            cfp = fopen("_codonscore", "r");
            loadcodonscore(cfp, codonscoremtx);
            fclose(cfp);
            for (i = 0; i < 64; i++)
                codonfreq[i] = 0.0;
            totalcount = 0.0;
            for (i = 3; i < len1; i++) {
                //				reporterr( "i=%d\n", i );
                if (gstart[i - 2] == STRONG * 0.5 && gend[i - 2] == 1.0 * 0.5 && gstart[i - 1] == STRONG * 0.5 && gend[i - 1] == STRONG * 0.5 && gstart[i - 0] == 1.0 * 0.5 && gend[i - 0] == STRONG * 0.5) {
                    //					reporterr( "codon=%.3s, id=%d\n", seq1[0]+i-2, codon2id(seq1[0]+i-2) );
                    codonfreq[codon2id(seq1[0] + i - 2)] += 1.0;
                    totalcount += 1.0;
                }
            }
            for (i = 0; i < 64; i++)
                codonfreq[i] /= totalcount;
            //			for( i=0; i<64; i++ ) reporterr( "%d, %f\n", i, codonfreq[i] );

            codonscore0av = 0.0;
            for (i = 0; i < 64; i++)
                for (j = 0; j < 64; j++)
                    codonscore0av += codonfreq[i] * codonfreq[j] * codonscoremtx[i][j];
            //			reporterr( "0av=%f\n", codonscore0av );
            for (i = 0; i < 64; i++)
                for (j = 0; j < 64; j++)
                    codonscoremtx[i][j] -= codonscore0av;

            codonscore1av = 0.0;
            for (i = 0; i < 64; i++)
                codonscore1av += codonfreq[i] * codonscoremtx[i][i];
            //			reporterr( "1av=%f\n", codonscore1av );
            for (i = 0; i < 64; i++)
                for (j = 0; j < 64; j++)
                    codonscoremtx[i][j] /= codonscore1av;

#if 0
			codonscore1av = 0.0;
			for( i=0; i<64; i++ ) codonscore1av += codonfreq[i] * codonscoremtx[i][i];
			reporterr( "1av=%f\n", codonscore1av );

			codonscore0av = 0.0;
			for( i=0; i<64; i++ ) for( j=0; j<64; j++ ) codonscore0av += codonfreq[i]*codonfreq[j] * codonscoremtx[i][j];
			reporterr( "0av=%f\n", codonscore0av );
#endif

            free(codonfreq);
        } else
            codonscoremtx = NULL;
    }

    if (localalloclen < nlen) {
        if (localalloclen) {
#if 1
            if (!kobetsubunkatsu) {
                FreeFukusosuuMtx(seqVector1);
                FreeFukusosuuMtx(seqVector2);
                FreeFukusosuuVec(naisekiNoWa);
                FreeFukusosuuMtx(naiseki);
                FreeDoubleVec(soukan);
            }
#endif
        }

        if (!kobetsubunkatsu) {
            naisekiNoWa = AllocateFukusosuuVec(nlen);
            naiseki = AllocateFukusosuuMtx(n20or4or2, nlen);
            seqVector1 = AllocateFukusosuuMtx(n20or4or2 + 1, nlen + 1);
            seqVector2 = AllocateFukusosuuMtx(n20or4or2 + 1, nlen + 1);
            soukan = AllocateDoubleVec(nlen + 1);
        }
        localalloclen = nlen;
    }

    for (j = 0; j < clus1; j++)
        strcpy(tmpseq1[j], seq1[j]);
    for (j = 0; j < clus2; j++)
        strcpy(tmpseq2[j], seq2[j]);

#if 0
fftfp = fopen( "input_of_Falign", "w" );
fprintf( fftfp, "nlen = %d\n", nlen );
fprintf( fftfp, "seq1: ( %d sequences ) \n", clus1 );
for( i=0; i<clus1; i++ )
	fprintf( fftfp, "%s\n", seq1[i] );
fprintf( fftfp, "seq2: ( %d sequences ) \n", clus2 );
for( i=0; i<clus2; i++ )
	fprintf( fftfp, "%s\n", seq2[i] );
fclose( fftfp );
system( "less input_of_Falign < /dev/tty > /dev/tty" );
#endif
    if (!kobetsubunkatsu) {
        if (fftkeika)
            fprintf(stderr, " FFT ... ");

        for (j = 0; j < n20or4or2; j++)
            vec_init(seqVector1[j], nlen);
        if (fftscore && scoremtx != -1) {
            for (i = 0; i < clus1; i++) {
#if 1
                seq_vec_5(ctx, seqVector1[0], polarity, volume, eff1[i], tmpseq1[i]);
#else
                seq_vec_2(seqVector1[0], polarity, eff1[i], tmpseq1[i]);
                seq_vec_2(seqVector1[1], volume, eff1[i], tmpseq1[i]);
#endif
            }
        } else {
#if 0
			for( i=0; i<clus1; i++ ) for( j=0; j<n20or4or2; j++ ) 
				seq_vec( seqVector1[j], amino[j], eff1[i], tmpseq1[i] );
#else
            for (i = 0; i < clus1; i++)
                seq_vec_3(ctx, seqVector1, eff1[i], tmpseq1[i]);
#endif
        }
#if RND
        for (i = 0; i < clus1; i++) {
            vec_init2(seqVector1, rndseq1[i], eff1[i], len1, nlen);
        }
#endif
#if 0
fftfp = fopen( "seqVec", "w" );
fprintf( fftfp, "before transform\n" );
for( k=0; k<n20or4or2; k++ ) 
{
   fprintf( fftfp, "nlen=%d\n", nlen );
   fprintf( fftfp, "%c\n", amino[k] );
   for( l=0; l<nlen; l++ )
   fprintf( fftfp, "%f %f\n", seqVector1[k][l].R, seqVector1[k][l].I );
}
fclose( fftfp );
system( "less seqVec < /dev/tty > /dev/tty" );
#endif

        for (j = 0; j < n20or4or2; j++)
            vec_init(seqVector2[j], nlen);
        if (fftscore && scoremtx != -1) {
            for (i = 0; i < clus2; i++) {
#if 1
                seq_vec_5(ctx, seqVector2[0], polarity, volume, eff2[i], tmpseq2[i]);
#else
                seq_vec_2(seqVector2[0], polarity, eff2[i], tmpseq2[i]);
                seq_vec_2(seqVector2[1], volume, eff2[i], tmpseq2[i]);
#endif
            }
        } else {
#if 0
			for( i=0; i<clus2; i++ ) for( j=0; j<n20or4or2; j++ ) 
				seq_vec( seqVector2[j], amino[j], eff2[i], tmpseq2[i] );
#else
            for (i = 0; i < clus2; i++)
                seq_vec_3(ctx, seqVector2, eff2[i], tmpseq2[i]);
#endif
        }
#if RND
        for (i = 0; i < clus2; i++) {
            vec_init2(seqVector2, rndseq2[i], eff2[i], len2, nlen);
        }
#endif

#if 0
fftfp = fopen( "seqVec2", "w" );
fprintf( fftfp, "before fft\n" );
for( k=0; k<n20or4or2; k++ ) 
{
   fprintf( fftfp, "%c\n", amino[k] );
   for( l=0; l<nlen; l++ )
   fprintf( fftfp, "%f %f\n", seqVector2[k][l].R, seqVector2[k][l].I );
}
fclose( fftfp );
system( "less seqVec2 < /dev/tty > /dev/tty" );
#endif

        for (j = 0; j < n20or4or2; j++) {
            fft(nlen, seqVector2[j], 0);
            fft(nlen, seqVector1[j], 0);
        }
#if 0
fftfp = fopen( "seqVec2", "w" );
fprintf( fftfp, "#after fft\n" );
for( k=0; k<n20or4or2; k++ ) 
{
   fprintf( fftfp, "#%c\n", amino[k] );
   for( l=0; l<nlen; l++ )
	   fprintf( fftfp, "%f %f\n", seqVector2[k][l].R, seqVector2[k][l].I );
}
fclose( fftfp );
system( "less seqVec2 < /dev/tty > /dev/tty" );
#endif

        for (k = 0; k < n20or4or2; k++) {
            for (l = 0; l < nlen; l++)
                calcNaiseki(naiseki[k] + l, seqVector1[k] + l, seqVector2[k] + l);
        }
        for (l = 0; l < nlen; l++) {
            naisekiNoWa[l].R = 0.0;
            naisekiNoWa[l].I = 0.0;
            for (k = 0; k < n20or4or2; k++) {
                naisekiNoWa[l].R += naiseki[k][l].R;
                naisekiNoWa[l].I += naiseki[k][l].I;
            }
        }

#if 0
	fftfp = fopen( "naisekiNoWa", "w" );
	fprintf( fftfp, "#Before fft\n" );
	for( l=0; l<nlen; l++ )
		fprintf( fftfp, "%d  %f %f\n", l, naisekiNoWa[l].R, naisekiNoWa[l].I ); 
	fclose( fftfp );
	system( "less naisekiNoWa < /dev/tty > /dev/tty " );
#endif

        fft(-nlen, naisekiNoWa, 0);

        for (m = 0; m <= nlen2; m++)
            soukan[m] = naisekiNoWa[nlen2 - m].R;
        for (m = nlen2 + 1; m < nlen; m++)
            soukan[m] = naisekiNoWa[nlen + nlen2 - m].R;

#if 0
	fftfp = fopen( "naisekiNoWa", "w" );
	fprintf( fftfp, "#After fft\n" );
	for( l=0; l<nlen; l++ )
		fprintf( fftfp, "%d  %f\n", l, naisekiNoWa[l].R ); 
	fclose( fftfp );
	fftfp = fopen( "list.plot", "w"  );
	fprintf( fftfp, "plot 'naisekiNoWa'\npause -1" );
	fclose( fftfp );
	system( "/usr/bin/gnuplot list.plot &" );
#endif
#if 0
	fprintf( stderr, "soukan\n" );
	for( l=0; l<nlen; l++ )
		fprintf( stderr, "%d  %f\n", l-nlen2, soukan[l] );
#if 0
	fftfp = fopen( "list.plot", "w"  );
	fprintf( fftfp, "plot 'frt'\n pause +1" );
	fclose( fftfp );
	system( "/usr/bin/gnuplot list.plot" );
#endif
#endif

        getKouho(kouho, NKOUHO, soukan, nlen);

#if 0
		for( i=0; i<NKOUHO; i++ )
		{
			fprintf( stderr, "kouho[%d] = %d\n", i, kouho[i] );
		}
#endif
    }

#if KEIKA
    fprintf(stderr, "Searching anchors ... ");
#endif
    count = 0;

#define CAND 0
#if CAND
    fftfp = fopen("cand", "w");
    fclose(fftfp);
#endif
    if (kobetsubunkatsu) {
        maxk = 1;
        kouho[0] = 0;
    } else {
        maxk = NKOUHO;
    }

    for (k = 0; k < maxk; k++) {
        lag = kouho[k];
        if (lag <= -len1 || len2 <= lag)
            continue;
        zurasu2(lag, clus1, clus2, seq1, seq2, tmpptr1, tmpptr2);
#if CAND
        fftfp = fopen("cand", "a");
        fprintf(fftfp, ">Candidate No.%d lag = %d\n", k + 1, lag);
        fprintf(fftfp, "%s\n", tmpptr1[0]);
        fprintf(fftfp, ">Candidate No.%d lag = %d\n", k + 1, lag);
        fprintf(fftfp, "%s\n", tmpptr2[0]);
        fprintf(fftfp, ">\n", k + 1, lag);
        fclose(fftfp);
#endif

        //		fprintf( stderr, "lag = %d\n", lag );
        tmpint = alignableReagion(ctx, clus1, clus2, tmpptr1, tmpptr2, eff1, eff2, segment + count);

        //		if( lag == -50 ) exit( 1 );

        if (count + tmpint > MAXSEG - 3)
            ErrorExit("TOO MANY SEGMENTS.\n");

        if (tmpint == 0)
            break;  // 060430 iinoka ?
        while (tmpint-- > 0) {
#if 0
			if( segment[count].end - segment[count].start < fftWinSize )
			{
				count++;
				continue;
			}
#endif
            if (lag > 0) {
                segment1[count].start = segment[count].start;
                segment1[count].end = segment[count].end;
                segment1[count].center = segment[count].center;
                segment1[count].score = segment[count].score;

                segment2[count].start = segment[count].start + lag;
                segment2[count].end = segment[count].end + lag;
                segment2[count].center = segment[count].center + lag;
                segment2[count].score = segment[count].score;
            } else {
                segment1[count].start = segment[count].start - lag;
                segment1[count].end = segment[count].end - lag;
                segment1[count].center = segment[count].center - lag;
                segment1[count].score = segment[count].score;

                segment2[count].start = segment[count].start;
                segment2[count].end = segment[count].end;
                segment2[count].center = segment[count].center;
                segment2[count].score = segment[count].score;
            }
#if 0
			fprintf( stderr, "in 1 %d\n", segment1[count].center );
			fprintf( stderr, "in 2 %d\n", segment2[count].center );
#endif
            segment1[count].pair = &segment2[count];
            segment2[count].pair = &segment1[count];
            count++;
        }
    }
#if 0
	if( !kobetsubunkatsu && fftkeika )
		fprintf( stderr, "%d anchors found\r", count );
#endif
    if (!count && fftNoAnchStop)
        ErrorExit("Cannot detect anchor!");
#if 0
	fprintf( stderr, "RESULT before sort:\n" );
	for( l=0; l<count+1; l++ )
	{
		fprintf( stderr, "cut[%d]=%d, ", l, segment1[l].center );
		fprintf( stderr, "%d score = %f\n", segment2[l].center, segment1[l].score );
	}
#endif

#if KEIKA
    fprintf(stderr, "done. (%d anchors)\n", count);
    fprintf(stderr, "Aligning anchors ... ");
#endif
    for (i = 0; i < count; i++) {
        sortedseg1[i] = &segment1[i];
        sortedseg2[i] = &segment2[i];
    }
#if 0
	tmpsort( count, sortedseg1 ); 
	tmpsort( count, sortedseg2 ); 
	qsort( sortedseg1, count, sizeof( Segment * ), segcmp );
	qsort( sortedseg2, count, sizeof( Segment * ), segcmp );
#else
    mymergesort(0, count - 1, sortedseg1);
    mymergesort(0, count - 1, sortedseg2);
#endif
    for (i = 0; i < count; i++)
        sortedseg1[i]->number = i;
    for (i = 0; i < count; i++)
        sortedseg2[i]->number = i;

    if (kobetsubunkatsu) {
        for (i = 0; i < count; i++) {
            cut1[i + 1] = sortedseg1[i]->center;
            cut2[i + 1] = sortedseg2[i]->center;
        }
        cut1[0] = 0;
        cut2[0] = 0;
        cut1[count + 1] = len1;
        cut2[count + 1] = len2;
        count += 2;
    } else {
        if (crossscoresize < count + 2) {
            crossscoresize = count + 2;
#if 1
            if (fftkeika)
                fprintf(stderr, "######allocating crossscore, size = %d\n", crossscoresize);
#endif
            if (crossscore)
                FreeDoubleMtx(crossscore);
            crossscore = AllocateDoubleMtx(crossscoresize, crossscoresize);
        }
        for (i = 0; i < count + 2; i++)
            for (j = 0; j < count + 2; j++)
                crossscore[i][j] = 0.0;
        for (i = 0; i < count; i++) {
            crossscore[segment1[i].number + 1][segment1[i].pair->number + 1] = segment1[i].score;
            cut1[i + 1] = sortedseg1[i]->center;
            cut2[i + 1] = sortedseg2[i]->center;
        }

#if 0
		fprintf( stderr, "AFTER SORT\n" );
		for( i=0; i<count+1; i++ ) fprintf( stderr, "%d, %d\n", cut1[i], cut2[i] );
		fprintf( stderr, "crossscore = \n" );
		for( i=0; i<count+1; i++ )
		{
			for( j=0; j<count+1; j++ )
				fprintf( stderr, "%.0f ", crossscore[i][j] );
			fprintf( stderr, "\n" );
		}
#endif

        crossscore[0][0] = 10000000.0;
        cut1[0] = 0;
        cut2[0] = 0;
        crossscore[count + 1][count + 1] = 10000000.0;
        cut1[count + 1] = len1;
        cut2[count + 1] = len2;
        count += 2;
        count0 = count;

        blockAlign2(cut1, cut2, sortedseg1, sortedseg2, crossscore, &count);

        //		if( count-count0 )
        //			fprintf( stderr, "%d unused anchors\n", count0-count );

        if (!kobetsubunkatsu && fftkeika)
            fprintf(stderr, "%d anchors found\n", count);
        if (fftkeika) {
            if (count0 > count) {
#if 0
				fprintf( stderr, "\7 REPEAT!? \n" );
#else
                fprintf(stderr, "REPEAT!? \n");
#endif
                if (fftRepeatStop)
                    exit(1);
            }
#if KEIKA
            else
                fprintf(stderr, "done\n");
#endif
        }
    }

#if 0
	fftfp = fopen( "fft", "a" );
	fprintf( fftfp, "RESULT after sort:\n" );
	for( l=0; l<count; l++ )
	{
		fprintf( fftfp, "cut[%d]=%d, ", l, segment1[l].center );
		fprintf( fftfp, "%d\n", segment2[l].center );
	}
	fclose( fftfp );
#endif

#if 0
	fprintf( stderr, "RESULT after blckalign:\n" );
	for( l=0; l<count+1; l++ )
	{
		fprintf( stderr, "cut : %d %d\n", cut1[l], cut2[l] );
	}
#endif

#if 0
	fprintf( trap_g, "Devided to %d segments\n", count-1 );
	fprintf( trap_g, "%d  %d forg\n", MIN( clus1, clus2 ), count-1 );
#endif

    totallen = 0;
    for (j = 0; j < clus1; j++)
        result1[j][0] = 0;
    for (j = 0; j < clus2; j++)
        result2[j][0] = 0;
    totalscore = 0.0;
    *fftlog = -1;
    for (i = 0; i < count - 1; i++) {
        *fftlog += 1;
        if (i == 0)
            headgp = outgap;
        else
            headgp = 1;
        if (i == count - 2)
            tailgp = outgap;
        else
            tailgp = 1;

#if 0
		if( cut1[i] ) // chuui
		{
//			getkyokaigap( sgap1, tmpres1, nlen-1, clus1 );
//			getkyokaigap( sgap2, tmpres2, nlen-1, clus2 );
			getkyokaigap( sgap1, tmpres1, nlen-1, clus1 );
			getkyokaigap( sgap2, tmpres2, nlen-1, clus2 );
		}
		else
		{
			for( j=0; j<clus1; j++ ) sgap1[j] = 'o';
			for( j=0; j<clus2; j++ ) sgap2[j] = 'o';
		}
		if( cut1[i+1] != len1 ) // chuui
		{       
			getkyokaigap( egap1, seq1, cut1[i+1], clus1 );
			getkyokaigap( egap2, seq2, cut2[i+1], clus2 );
		}       
		else    
		{       
			for( j=0; j<clus1; j++ ) egap1[j] = 'o';
			for( j=0; j<clus2; j++ ) egap2[j] = 'o';
		}
#else
        if (cut1[i])
            getkyokaigap(sgap1, seq1, cut1[i] - 1, clus1);
        else
            for (j = 0; j < clus1; j++)
                sgap1[j] = 'o';
        if (cut2[i])
            getkyokaigap(sgap2, seq2, cut2[i] - 1, clus2);
        else
            for (j = 0; j < clus2; j++)
                sgap2[j] = 'o';

        if (cut1[i + 1] != len1)
            getkyokaigap(egap1, seq1, cut1[i + 1], clus1);
        else
            for (j = 0; j < clus1; j++)
                egap1[j] = 'o';
        if (cut2[i + 1] != len2)
            getkyokaigap(egap2, seq2, cut2[i + 1], clus2);
        else
            for (j = 0; j < clus2; j++)
                egap2[j] = 'o';
#endif
#if 0
		{
			fprintf( stderr, "kyokkaigap1(%d)=", cut1[i]-1 );
			for( j=0; j<clus1; j++ )
				fprintf( stderr, "%c", sgap1[j] );
			fprintf( stderr, "=kyokkaigap1-start\n" );
		}
		{
			fprintf( stderr, "kyokkaigap2(%d)=", cut2[i]-1 );
			for( j=0; j<clus2; j++ )
				fprintf( stderr, "%c", sgap2[j] );
			fprintf( stderr, "=kyokkaigap2-start\n" );
		}
		{
			fprintf( stderr, "kyokkaigap1(%d)=", cut1[i]-1 );
			for( j=0; j<clus1; j++ )
				fprintf( stderr, "%c", egap1[j] );
			fprintf( stderr, "=kyokkaigap1-end\n" );
		}
		{
			fprintf( stderr, "kyokkaigap2(%d)=", cut2[i]-1 );
			for( j=0; j<clus2; j++ )
				fprintf( stderr, "%c", egap2[j] );
			fprintf( stderr, "=kyokkaigap2-end\n" );
		}
#endif

#if DEBUG
        fprintf(stderr, "DP %03d / %03d %4d to ", i + 1, count - 1, totallen);
#else
#if KEIKA
        fprintf(stderr, "DP %03d / %03d\r", i + 1, count - 1);
#endif
#endif

        //reporterr( "cut1[] = %d\n", cut1[i] );
        //reporterr( "cut2[] = %d\n", cut2[i] );

        for (j = 0; j < clus1; j++) {
            strncpy(tmpres1[j], seq1[j] + cut1[i], cut1[i + 1] - cut1[i]);
            tmpres1[j][cut1[i + 1] - cut1[i]] = 0;
        }
        if (kobetsubunkatsu && fftkeika)
            commongappick(clus1, tmpres1);  //dvtditr $B$K8F$P$l$?$H$-(B fftkeika=1
        //		if( kobetsubunkatsu ) commongappick( clus1, tmpres1 );
        for (j = 0; j < clus2; j++) {
            strncpy(tmpres2[j], seq2[j] + cut2[i], cut2[i + 1] - cut2[i]);
            tmpres2[j][cut2[i + 1] - cut2[i]] = 0;
        }
        if (kobetsubunkatsu && fftkeika)
            commongappick(clus2, tmpres2);  //dvtditr $B$K8F$P$l$?$H$-(B fftkeika=1
        //		if( kobetsubunkatsu ) commongappick( clus2, tmpres2 );

        if (constraint) {
            fprintf(stderr, "Not supported\n");
            exit(1);
        }

        switch (alg) {
            case ('a'):
                totalscore += Aalign(ctx, tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen);
                break;
            case ('M'):
                if (scoringmatrices)
                    totalscore += MSalignmm_variousdist(ctx, scoringmatrices, tmpres1, tmpres2, eff1, eff2, eff1s, eff2s, clus1, clus2, alloclen, sgap1, sgap2, egap1, egap2, headgp, tailgp);
                else
                    totalscore += MSalignmm(ctx, n_dynamicmtx, tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen, sgap1, sgap2, egap1, egap2, headgp, tailgp, NULL, NULL, NULL, 0.0, 0.0);
                break;
            case ('d'):
                if (clus1 == 1 && clus2 == 1) {
                    totalscore += G__align11(ctx, n_dynamicmtx, tmpres1, tmpres2, alloclen, headgp, tailgp);
                } else {
                    if (scoringmatrices)  // called by tditeration.c
                    {
                        totalscore += D__align_variousdist(ctx, whichmtx, scoringmatrices, tmpres1, tmpres2, eff1, eff2, eff1s, eff2s, clus1, clus2, alloclen, 0, &dumdb, headgp, tailgp);
                    } else
                        totalscore += D__align(ctx, n_dynamicmtx, tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen, 0, &dumdb, headgp, tailgp);
                }
                break;
            case ('A'):
                if (clus1 == 1 && clus2 == 1) {
                    if (codonpos || codonscore) {
                        totalscore += G__align11psg(ctx, codonscoremtx, n_dynamicmtx, tmpres1, tmpres2, alloclen, headgp, tailgp, gstart + cut1[i], gend + cut1[i]);
                    } else
                        totalscore += G__align11(ctx, n_dynamicmtx, tmpres1, tmpres2, alloclen, headgp, tailgp);
                } else {
                    if (codonpos) {
                        reporterr("\n\ncodonpos will be soon supported for a reference MSA. For now, use a single sequence as reference.\n\n\n");
                        exit(1);
                    }
                    if (scoringmatrices)  // called by tditeration.c
                    {
                        totalscore += A__align_variousdist(ctx, whichmtx, scoringmatrices, penalty, penalty_ex, tmpres1, tmpres2, eff1, eff2, eff1s, eff2s, clus1, clus2, alloclen, 0, &dumdb, sgap1, sgap2, egap1, egap2, headgp, tailgp);
                    } else
                        totalscore += A__align(ctx, n_dynamicmtx, penalty, penalty_ex, tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen, 0, &dumdb, sgap1, sgap2, egap1, egap2, headgp, tailgp, -1, -1, NULL, NULL, NULL, 0.0, 0.0);
                }
                break;
            default:
                fprintf(stderr, "alg = %c\n", alg);
                ErrorExit("ERROR IN SOURCE FILE Falign.c");
                break;
        }

#ifdef enablemultithread
        if (chudanres && *chudanres) {
            //			fprintf( stderr, "\n\n## CHUUDAN!!! at Falign_localhom\n" );
            //			Added 2021/Jul/25.
            FreeCharMtx(result1);
            FreeCharMtx(result2);
            FreeCharMtx(tmpres1);
            FreeCharMtx(tmpres2);
            FreeCharMtx(tmpseq1);
            FreeCharMtx(tmpseq2);
            free(sgap1);
            free(egap1);
            free(sgap2);
            free(egap2);
            free(tmpptr1);
            free(tmpptr2);
#if RND
            FreeCharMtx(rndseq1);
            FreeCharMtx(rndseq2);
#endif
            return (-1.0);
        }
#endif

        nlen = strlen(tmpres1[0]);
        if (totallen + nlen > alloclen) {
            fprintf(stderr, "totallen=%d +  nlen=%d > alloclen = %d\n", totallen, nlen, alloclen);
            ErrorExit("LENGTH OVER in Falign\n ");
        }
        for (j = 0; j < clus1; j++)
            strcat(result1[j], tmpres1[j]);
        for (j = 0; j < clus2; j++)
            strcat(result2[j], tmpres2[j]);
        totallen += nlen;
#if 0
		fprintf( stderr, "$#####$$$$ i=%d", i );
		fprintf( stderr, "%4d\n", totallen );
		fprintf( stderr, "\n\n" );
		for( j=0; j<clus1; j++ ) 
		{
			fprintf( stderr, "%s\n", tmpres1[j] );
		}
		fprintf( stderr, "-------\n" );
		for( j=0; j<clus2; j++ ) 
		{
			fprintf( stderr, "%s\n", tmpres2[j] );
		}
#endif
    }

#if KEIKA
    fprintf(stderr, "DP ... done   \n");
#endif

    for (j = 0; j < clus1; j++)
        strcpy(seq1[j], result1[j]);
    for (j = 0; j < clus2; j++)
        strcpy(seq2[j], result2[j]);
#if 0
	for( j=0; j<clus1; j++ ) 
	{
		fprintf( stderr, "in Falign, %s\n", result1[j] );
	}
	fprintf( stderr, "- - - - - - - - - - -\n" );
	for( j=0; j<clus2; j++ ) 
	{
		fprintf( stderr, "in Falign, %s\n", result2[j] );
	}
#endif

    FreeCharMtx(result1);
    FreeCharMtx(result2);
    FreeCharMtx(tmpres1);
    FreeCharMtx(tmpres2);
    FreeCharMtx(tmpseq1);
    FreeCharMtx(tmpseq2);
    free(sgap1);
    free(egap1);
    free(sgap2);
    free(egap2);
    free(tmpptr1);
    free(tmpptr2);
#if RND
    FreeCharMtx(rndseq1);
    FreeCharMtx(rndseq2);
#endif

    return (totalscore);
}

double
Falign_udpari_long(
    Context*  ctx,
    double*** scoringmatrices,
    double**  n_dynamicmtx,
    char**    seq1,
    char**    seq2,
    double*   eff1,
    double*   eff2,
    double**  eff1s,
    double**  eff2s,
    int       clus1,
    int       clus2,
    int       alloclen,
    int*      fftlog
) {
    int        i, j, k, l, m, maxk;
    int        nlen, nlen2;
    static int crossscoresize = 0;
    char**     tmpseq1 = NULL;
    char**     tmpseq2 = NULL;
    char**     tmpptr1 = NULL;
    char**     tmpptr2 = NULL;
    char**     tmpres1 = NULL;
    char**     tmpres2 = NULL;
    char**     result1 = NULL;
    char**     result2 = NULL;
#if RND
    char** rndseq1 = NULL;
    char** rndseq2 = NULL;
#endif
    static Fukusosuu** seqVector1 = NULL;
    static Fukusosuu** seqVector2 = NULL;
    static Fukusosuu** naiseki = NULL;
    static Fukusosuu*  naisekiNoWa = NULL;
    static double*     soukan = NULL;
    static double**    crossscore = NULL;
    int                nlentmp;
    static int*        kouho = NULL;
    static Segment*    segment = NULL;
    static Segment*    segment1 = NULL;
    static Segment*    segment2 = NULL;
    static Segment**   sortedseg1 = NULL;
    static Segment**   sortedseg2 = NULL;
    static int*        cut1 = NULL;
    static int*        cut2 = NULL;
    char *             sgap1, *egap1, *sgap2, *egap2;
    static int         localalloclen = 0;
    int                lag;
    int                tmpint;
    int                count, count0;
    int                len1, len2;
    int                totallen;
    double             totalscore;
    int                nkouho = 0;
    int                headgp, tailgp;
    //	double dumfl = 0.0;

    if (seq1 == NULL) {
        if (kouho) {
            //			fprintf( stderr, "### Freeing localarrays in Falign\n" );
            localalloclen = 0;
            crossscoresize = 0;
            mymergesort(0, 0, NULL);
            alignableReagion(ctx, 0, 0, NULL, NULL, NULL, NULL, NULL);
            fft(0, NULL, 1);
            A__align(ctx, NULL, 0, 0, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, 0, 0, -1, -1, NULL, NULL, NULL, 0.0, 0.0);
            A__align_variousdist(ctx, NULL, NULL, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, 0, 0);
            D__align_variousdist(ctx, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, 0, NULL, 0, 0);
            G__align11(ctx, NULL, NULL, NULL, 0, 0, 0);
            blockAlign2(NULL, NULL, NULL, NULL, NULL, NULL);
            if (crossscore)
                FreeDoubleMtx(crossscore);
            crossscore = NULL;  // reallocate sareru kanousei ga arunode.
            free(kouho);
            kouho = NULL;
            free(cut1);
            free(cut2);
            free(segment);
            free(segment1);
            free(segment2);
            free(sortedseg1);
            free(sortedseg2);
            if (!kobetsubunkatsu) {
                FreeFukusosuuMtx(seqVector1);
                FreeFukusosuuMtx(seqVector2);
                FreeFukusosuuVec(naisekiNoWa);
                FreeFukusosuuMtx(naiseki);
                FreeDoubleVec(soukan);
            }
        } else {
            //			fprintf( stderr, "Did not allocate localarrays in Falign\n" );
        }

        return (0.0);
    }

    len1 = strlen(seq1[0]);
    len2 = strlen(seq2[0]);
    nlentmp = MAX(len1, len2);

    nlen = 1;
    while (nlentmp >= nlen)
        nlen <<= 1;
#if 0
	fprintf( stderr, "###   nlen    = %d\n", nlen );
#endif

    nlen2 = nlen / 2;

#if 0
	fprintf( stderr, "len1 = %d, len2 = %d\n", len1, len2 );
	fprintf( stderr, "nlentmp = %d, nlen = %d\n", nlentmp, nlen );
#endif

    result1 = AllocateCharMtx(clus1, alloclen);
    result2 = AllocateCharMtx(clus2, alloclen);
    tmpres1 = AllocateCharMtx(clus1, alloclen);
    tmpres2 = AllocateCharMtx(clus2, alloclen);
    sgap1 = AllocateCharVec(clus1);
    egap1 = AllocateCharVec(clus1);
    sgap2 = AllocateCharVec(clus2);
    egap2 = AllocateCharVec(clus2);

    tmpseq1 = AllocateCharMtx(clus1, nlen);
    tmpseq2 = AllocateCharMtx(clus2, nlen);
    tmpptr1 = calloc(clus1, sizeof(char*));
    tmpptr2 = calloc(clus2, sizeof(char*));

#if RND
    rndseq1 = AllocateCharMtx(clus1, nlen);
    rndseq2 = AllocateCharMtx(clus2, nlen);
    for (i = 0; i < clus1; i++)
        generateRndSeq(rndseq1[i], nlen);
    for (i = 0; i < clus2; i++)
        generateRndSeq(rndseq2[i], nlen);
#endif

    if (!localalloclen) {
        kouho = AllocateIntVec(NKOUHO_LONG);
        cut1 = AllocateIntVec(MAXSEG);
        cut2 = AllocateIntVec(MAXSEG);
        segment = (Segment*)calloc(MAXSEG, sizeof(Segment));
        segment1 = (Segment*)calloc(MAXSEG, sizeof(Segment));
        segment2 = (Segment*)calloc(MAXSEG, sizeof(Segment));
        sortedseg1 = (Segment**)calloc(MAXSEG, sizeof(Segment*));
        sortedseg2 = (Segment**)calloc(MAXSEG, sizeof(Segment*));
        if (!(segment && segment1 && segment2 && sortedseg1 && sortedseg2))
            ErrorExit("Allocation error\n");

        if (scoremtx == -1)
            n20or4or2 = 1;
        else if (fftscore)
            n20or4or2 = 1;
        else
            n20or4or2 = 20;
    }

    if (localalloclen < nlen) {
        if (localalloclen) {
#if 1
            if (!kobetsubunkatsu) {
                FreeFukusosuuMtx(seqVector1);
                FreeFukusosuuMtx(seqVector2);
                FreeFukusosuuVec(naisekiNoWa);
                FreeFukusosuuMtx(naiseki);
                FreeDoubleVec(soukan);
            }
#endif
        }

        if (!kobetsubunkatsu) {
            naisekiNoWa = AllocateFukusosuuVec(nlen);
            naiseki = AllocateFukusosuuMtx(n20or4or2, nlen);
            seqVector1 = AllocateFukusosuuMtx(n20or4or2, nlen + 1);
            seqVector2 = AllocateFukusosuuMtx(n20or4or2, nlen + 1);
            soukan = AllocateDoubleVec(nlen + 1);
        }
        localalloclen = nlen;
    }

    for (j = 0; j < clus1; j++)
        strcpy(tmpseq1[j], seq1[j]);
    for (j = 0; j < clus2; j++)
        strcpy(tmpseq2[j], seq2[j]);

#if 0
fftfp = fopen( "input_of_Falign", "w" );
fprintf( fftfp, "nlen = %d\n", nlen );
fprintf( fftfp, "seq1: ( %d sequences ) \n", clus1 );
for( i=0; i<clus1; i++ )
	fprintf( fftfp, "%s\n", seq1[i] );
fprintf( fftfp, "seq2: ( %d sequences ) \n", clus2 );
for( i=0; i<clus2; i++ )
	fprintf( fftfp, "%s\n", seq2[i] );
fclose( fftfp );
system( "less input_of_Falign < /dev/tty > /dev/tty" );
#endif
    if (!kobetsubunkatsu) {
        if (fftkeika)
            fprintf(stderr, " FFT ... ");

        for (j = 0; j < n20or4or2; j++)
            vec_init(seqVector1[j], nlen);
        if (scoremtx == -1) {
            for (i = 0; i < clus1; i++)
                seq_vec_4(seqVector1[0], eff1[i], tmpseq1[i]);
        } else if (fftscore) {
            for (i = 0; i < clus1; i++) {
#if 0
				seq_vec_2( seqVector1[0], polarity, eff1[i], tmpseq1[i] );
				seq_vec_2( seqVector1[1], volume,   eff1[i], tmpseq1[i] );
#else
                seq_vec_5(ctx, seqVector1[0], polarity, volume, eff1[i], tmpseq1[i]);
#endif
            }
        } else {
            for (i = 0; i < clus1; i++)
                seq_vec_3(ctx, seqVector1, eff1[i], tmpseq1[i]);
        }
#if RND
        for (i = 0; i < clus1; i++) {
            vec_init2(seqVector1, rndseq1[i], eff1[i], len1, nlen);
        }
#endif
#if 0
fftfp = fopen( "seqVec", "w" );
fprintf( fftfp, "before transform\n" );
for( k=0; k<n20or4or2; k++ ) 
{
   fprintf( fftfp, "nlen=%d\n", nlen );
   fprintf( fftfp, "%c\n", amino[k] );
   for( l=0; l<nlen; l++ )
   fprintf( fftfp, "%f %f\n", seqVector1[k][l].R, seqVector1[k][l].I );
}
fclose( fftfp );
system( "less seqVec < /dev/tty > /dev/tty" );
#endif

        for (j = 0; j < n20or4or2; j++)
            vec_init(seqVector2[j], nlen);
        if (scoremtx == -1) {
            for (i = 0; i < clus2; i++)
                seq_vec_4(seqVector2[0], eff2[i], tmpseq2[i]);
        } else if (fftscore) {
            for (i = 0; i < clus2; i++) {
#if 0
				seq_vec_2( seqVector2[0], polarity, eff2[i], tmpseq2[i] );
				seq_vec_2( seqVector2[1], volume,   eff2[i], tmpseq2[i] );
#else
                seq_vec_5(ctx, seqVector2[0], polarity, volume, eff2[i], tmpseq2[i]);
#endif
            }
        } else {
            for (i = 0; i < clus2; i++)
                seq_vec_3(ctx, seqVector2, eff2[i], tmpseq2[i]);
        }
#if RND
        for (i = 0; i < clus2; i++) {
            vec_init2(seqVector2, rndseq2[i], eff2[i], len2, nlen);
        }
#endif

#if 0
fftfp = fopen( "seqVec2", "w" );
fprintf( fftfp, "before fft\n" );
for( k=0; k<n20or4or2; k++ ) 
{
   fprintf( fftfp, "%c\n", amino[k] );
   for( l=0; l<nlen; l++ )
   fprintf( fftfp, "%f %f\n", seqVector2[k][l].R, seqVector2[k][l].I );
}
fclose( fftfp );
system( "less seqVec2 < /dev/tty > /dev/tty" );
#endif

        for (j = 0; j < n20or4or2; j++) {
            fft(nlen, seqVector2[j], 0);
            fft(nlen, seqVector1[j], 0);
        }
#if 0
fftfp = fopen( "seqVec2", "w" );
fprintf( fftfp, "#after fft\n" );
for( k=0; k<n20or4or2; k++ ) 
{
   fprintf( fftfp, "#%c\n", amino[k] );
   for( l=0; l<nlen; l++ )
	   fprintf( fftfp, "%f %f\n", seqVector2[k][l].R, seqVector2[k][l].I );
}
fclose( fftfp );
system( "less seqVec2 < /dev/tty > /dev/tty" );
#endif

        for (k = 0; k < n20or4or2; k++) {
            for (l = 0; l < nlen; l++)
                calcNaiseki(naiseki[k] + l, seqVector1[k] + l, seqVector2[k] + l);
        }
        for (l = 0; l < nlen; l++) {
            naisekiNoWa[l].R = 0.0;
            naisekiNoWa[l].I = 0.0;
            for (k = 0; k < n20or4or2; k++) {
                naisekiNoWa[l].R += naiseki[k][l].R;
                naisekiNoWa[l].I += naiseki[k][l].I;
            }
        }

#if 0
	fftfp = fopen( "naisekiNoWa", "w" );
	fprintf( fftfp, "#Before fft\n" );
	for( l=0; l<nlen; l++ )
		fprintf( fftfp, "%d  %f %f\n", l, naisekiNoWa[l].R, naisekiNoWa[l].I ); 
	fclose( fftfp );
	system( "less naisekiNoWa < /dev/tty > /dev/tty " );
#endif

        fft(-nlen, naisekiNoWa, 0);

        for (m = 0; m <= nlen2; m++)
            soukan[m] = naisekiNoWa[nlen2 - m].R;
        for (m = nlen2 + 1; m < nlen; m++)
            soukan[m] = naisekiNoWa[nlen + nlen2 - m].R;

#if 0
	fftfp = fopen( "naisekiNoWa", "w" );
	fprintf( fftfp, "#After fft\n" );
	for( l=0; l<nlen; l++ )
		fprintf( fftfp, "%d  %f\n", l, naisekiNoWa[l].R ); 
	fclose( fftfp );
	fftfp = fopen( "list.plot", "w"  );
	fprintf( fftfp, "plot 'naisekiNoWa'\npause -1" );
	fclose( fftfp );
	system( "/usr/bin/gnuplot list.plot &" );
#endif
#if 0
	fprintf( stderr, "soukan\n" );
	for( l=0; l<nlen; l++ )
		fprintf( stderr, "%d  %f\n", l-nlen2, soukan[l] );
#if 0
	fftfp = fopen( "list.plot", "w"  );
	fprintf( fftfp, "plot 'frt'\n pause +1" );
	fclose( fftfp );
	system( "/usr/bin/gnuplot list.plot" );
#endif
#endif

        nkouho = getKouho(kouho, NKOUHO_LONG, soukan, nlen);

#if 0
		for( i=0; i<nkouho; i++ )
		{
			fprintf( stderr, "kouho[%d] = %d\n", i, kouho[i] );
		}
#endif
    }

#if KEIKA
    fprintf(stderr, "Searching anchors ... ");
#endif
    count = 0;

#define CAND 0
#if CAND
    fftfp = fopen("cand", "w");
    fclose(fftfp);
#endif
    if (kobetsubunkatsu) {
        maxk = 1;
        kouho[0] = 0;
    } else {
        maxk = nkouho;
    }

    for (k = 0; k < maxk; k++) {
        lag = kouho[k];
        if (lag <= -len1 || len2 <= lag)
            continue;
        //		fprintf( stderr, "k=%d, lag=%d\n", k, lag );
        zurasu2(lag, clus1, clus2, seq1, seq2, tmpptr1, tmpptr2);
#if CAND
        fftfp = fopen("cand", "a");
        fprintf(fftfp, ">Candidate No.%d lag = %d\n", k + 1, lag);
        fprintf(fftfp, "%s\n", tmpptr1[0]);
        fprintf(fftfp, ">Candidate No.%d lag = %d\n", k + 1, lag);
        fprintf(fftfp, "%s\n", tmpptr2[0]);
        fprintf(fftfp, ">\n", k + 1, lag);
        fclose(fftfp);
#endif

        //		fprintf( stderr, "lag = %d\n", lag );
        tmpint = alignableReagion(ctx, clus1, clus2, tmpptr1, tmpptr2, eff1, eff2, segment + count);
        //		fprintf( stderr, "lag = %d, %d found\n", lag, tmpint );

        //		if( lag == -50 ) exit( 1 );

        if (count + tmpint > MAXSEG - 3)
            ErrorExit("TOO MANY SEGMENTS.\n");

        //		fprintf( stderr, "##### k=%d / %d\n", k, maxk );
        //		if( tmpint == 0 ) break; // 060430 iinoka ? // 090530 yameta
        while (tmpint-- > 0) {
#if 0
			if( segment[count].end - segment[count].start < fftWinSize )
			{
				count++;
				continue;
			}
#endif
            if (lag > 0) {
                segment1[count].start = segment[count].start;
                segment1[count].end = segment[count].end;
                segment1[count].center = segment[count].center;
                segment1[count].score = segment[count].score;

                segment2[count].start = segment[count].start + lag;
                segment2[count].end = segment[count].end + lag;
                segment2[count].center = segment[count].center + lag;
                segment2[count].score = segment[count].score;
            } else {
                segment1[count].start = segment[count].start - lag;
                segment1[count].end = segment[count].end - lag;
                segment1[count].center = segment[count].center - lag;
                segment1[count].score = segment[count].score;

                segment2[count].start = segment[count].start;
                segment2[count].end = segment[count].end;
                segment2[count].center = segment[count].center;
                segment2[count].score = segment[count].score;
            }
#if 0
			fprintf( stderr, "##### k=%d / %d\n", k, maxk );
			fprintf( stderr, "anchor %d, score = %f\n", count, segment1[count].score );
			fprintf( stderr, "in 1 %d\n", segment1[count].center );
			fprintf( stderr, "in 2 %d\n", segment2[count].center );
#endif
            segment1[count].pair = &segment2[count];
            segment2[count].pair = &segment1[count];
            count++;
#if 0
			fprintf( stderr, "count=%d\n", count );
#endif
        }
    }
#if 1
    if (!kobetsubunkatsu)
        if (fftkeika)
            fprintf(stderr, "done. (%d anchors) ", count);
#endif
    if (!count && fftNoAnchStop)
        ErrorExit("Cannot detect anchor!");
#if 0
	fprintf( stderr, "RESULT before sort:\n" );
	for( l=0; l<count+1; l++ )
	{
		fprintf( stderr, "cut[%d]=%d, ", l, segment1[l].center );
		fprintf( stderr, "%d score = %f\n", segment2[l].center, segment1[l].score );
	}
#endif

    for (i = 0; i < count; i++) {
        sortedseg1[i] = &segment1[i];
        sortedseg2[i] = &segment2[i];
    }
#if 0
	tmpsort( count, sortedseg1 ); 
	tmpsort( count, sortedseg2 ); 
	qsort( sortedseg1, count, sizeof( Segment * ), segcmp );
	qsort( sortedseg2, count, sizeof( Segment * ), segcmp );
#else
    mymergesort(0, count - 1, sortedseg1);
    mymergesort(0, count - 1, sortedseg2);
#endif
    for (i = 0; i < count; i++)
        sortedseg1[i]->number = i;
    for (i = 0; i < count; i++)
        sortedseg2[i]->number = i;

    if (kobetsubunkatsu) {
        for (i = 0; i < count; i++) {
            cut1[i + 1] = sortedseg1[i]->center;
            cut2[i + 1] = sortedseg2[i]->center;
        }
        cut1[0] = 0;
        cut2[0] = 0;
        cut1[count + 1] = len1;
        cut2[count + 1] = len2;
        count += 2;
    }

    else {
        if (count < 5000) {
            if (crossscoresize < count + 2) {
                crossscoresize = count + 2;
#if 1
                if (fftkeika)
                    fprintf(stderr, "######allocating crossscore, size = %d\n", crossscoresize);
#endif
                if (crossscore)
                    FreeDoubleMtx(crossscore);
                crossscore = AllocateDoubleMtx(crossscoresize, crossscoresize);
            }
            for (i = 0; i < count + 2; i++)
                for (j = 0; j < count + 2; j++)
                    crossscore[i][j] = 0.0;
            for (i = 0; i < count; i++) {
                crossscore[segment1[i].number + 1][segment1[i].pair->number + 1] = segment1[i].score;
                cut1[i + 1] = sortedseg1[i]->center;
                cut2[i + 1] = sortedseg2[i]->center;
            }

#if 0
			fprintf( stderr, "AFTER SORT\n" );
			for( i=0; i<count+1; i++ ) fprintf( stderr, "%d, %d\n", cut1[i], cut2[i] );
			fprintf( stderr, "crossscore = \n" );
			for( i=0; i<count+1; i++ )
			{
				for( j=0; j<count+1; j++ )
					fprintf( stderr, "%.0f ", crossscore[i][j] );
				fprintf( stderr, "\n" );
			}
#endif

            crossscore[0][0] = 10000000.0;
            cut1[0] = 0;
            cut2[0] = 0;
            crossscore[count + 1][count + 1] = 10000000.0;
            cut1[count + 1] = len1;
            cut2[count + 1] = len2;
            count += 2;
            count0 = count;

            //			fprintf( stderr, "\n\n\ncalling blockAlign2\n\n\n\n" );
            blockAlign2(cut1, cut2, sortedseg1, sortedseg2, crossscore, &count);

            //			if( count-count0 )
            //				fprintf( stderr, "%d unused anchors\n", count0-count );

            if (!kobetsubunkatsu && fftkeika)
                fprintf(stderr, "%d anchors found\n", count);
            if (fftkeika) {
                if (count0 > count) {
#if 0
					fprintf( stderr, "\7 REPEAT!? \n" );
#else
                    fprintf(stderr, "REPEAT!? \n");
#endif
                    if (fftRepeatStop)
                        exit(1);
                }
#if KEIKA
                else
                    fprintf(stderr, "done\n");
#endif
            }
        }

        else {
            fprintf(stderr, "\nMany anchors were found. The upper-level DP is skipped.\n\n");

            cut1[0] = 0;
            cut2[0] = 0;
            count0 = 0;
            for (i = 0; i < count; i++) {
                //				fprintf( stderr, "i=%d, %d-%d ?\n", i, sortedseg1[i]->center, sortedseg1[i]->pair->center );
                if (sortedseg1[i]->center > cut1[count0]
                    && sortedseg1[i]->pair->center > cut2[count0]) {
                    count0++;
                    cut1[count0] = sortedseg1[i]->center;
                    cut2[count0] = sortedseg1[i]->pair->center;
                } else {
                    if (i && sortedseg1[i]->score > sortedseg1[i - 1]->score) {
                        if (sortedseg1[i]->center > cut1[count0 - 1]
                            && sortedseg1[i]->pair->center > cut2[count0 - 1]) {
                            cut1[count0] = sortedseg1[i]->center;
                            cut2[count0] = sortedseg1[i]->pair->center;
                        } else {
                            //							count0--;
                        }
                    }
                }
            }
            //			if( count-count0 )
            //				fprintf( stderr, "%d anchors unused\n", count-count0 );
            cut1[count0 + 1] = len1;
            cut2[count0 + 1] = len2;
            count = count0 + 2;
            count0 = count;
        }
    }

    //	exit( 0 );

#if 0
	fftfp = fopen( "fft", "a" );
	fprintf( fftfp, "RESULT after sort:\n" );
	for( l=0; l<count; l++ )
	{
		fprintf( fftfp, "cut[%d]=%d, ", l, segment1[l].center );
		fprintf( fftfp, "%d\n", segment2[l].center );
	}
	fclose( fftfp );
#endif

#if 0
	fprintf( stderr, "RESULT after blckalign:\n" );
	for( l=0; l<count+1; l++ )
	{
		fprintf( stderr, "cut : %d %d\n", cut1[l], cut2[l] );
	}
#endif

#if 0
	fprintf( trap_g, "Devided to %d segments\n", count-1 );
	fprintf( trap_g, "%d  %d forg\n", MIN( clus1, clus2 ), count-1 );
#endif

    totallen = 0;
    for (j = 0; j < clus1; j++)
        result1[j][0] = 0;
    for (j = 0; j < clus2; j++)
        result2[j][0] = 0;
    totalscore = 0.0;
    *fftlog = -1;
    //	reporterr( "\nin Falign_udpari(), *fftlog = %d\n", *fftlog );
    for (i = 0; i < count - 1; i++) {
        *fftlog += 1;
        if (i == 0)
            headgp = outgap;
        else
            headgp = 1;
        if (i == count - 2)
            tailgp = outgap;
        else
            tailgp = 1;

#if 0
		if( cut1[i] )
		{
//			getkyokaigap( sgap1, seq1, cut1[i]-1, clus1 );
//			getkyokaigap( sgap2, seq2, cut2[i]-1, clus2 );
			getkyokaigap( sgap1, tmpres1, nlen-1, clus1 );
			getkyokaigap( sgap2, tmpres2, nlen-1, clus2 );
		}
		else
		{
			for( j=0; j<clus1; j++ ) sgap1[j] = 'o';
			for( j=0; j<clus2; j++ ) sgap2[j] = 'o';
		}
		if( cut1[i+1] != len1 )
		{       
			getkyokaigap( egap1, seq1, cut1[i+1], clus1 );
			getkyokaigap( egap2, seq2, cut2[i+1], clus2 );
		}       
		else    
		{       
			for( j=0; j<clus1; j++ ) egap1[j] = 'o';
			for( j=0; j<clus2; j++ ) egap2[j] = 'o';
		}
#else
        if (cut1[i])
            getkyokaigap(sgap1, seq1, cut1[i] - 1, clus1);
        else
            for (j = 0; j < clus1; j++)
                sgap1[j] = 'o';
        if (cut2[i])
            getkyokaigap(sgap2, seq2, cut2[i] - 1, clus2);
        else
            for (j = 0; j < clus2; j++)
                sgap2[j] = 'o';

        if (cut1[i + 1] != len1)
            getkyokaigap(egap1, seq1, cut1[i + 1], clus1);
        else
            for (j = 0; j < clus1; j++)
                egap1[j] = 'o';
        if (cut2[i + 1] != len2)
            getkyokaigap(egap2, seq2, cut2[i + 1], clus2);
        else
            for (j = 0; j < clus2; j++)
                egap2[j] = 'o';
#endif

#if DEBUG
        fprintf(stderr, "DP %03d / %03d %4d to ", i + 1, count - 1, totallen);
#else
#if 1
        if (1 || fftkeika)
            fprintf(stderr, "DP %05d / %05d \b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", i + 1, count - 1);
#endif
#endif
        for (j = 0; j < clus1; j++) {
            strncpy(tmpres1[j], seq1[j] + cut1[i], cut1[i + 1] - cut1[i]);
            tmpres1[j][cut1[i + 1] - cut1[i]] = 0;
        }
        if (kobetsubunkatsu && fftkeika)
            commongappick(clus1, tmpres1);  //dvtditr $B$K8F$P$l$?$H$-(B fftkeika=1
        //		if( kobetsubunkatsu ) commongappick( clus1, tmpres1 );
        for (j = 0; j < clus2; j++) {
            //			fprintf( stderr, "### cut2[i+1]-cut2[i] = %d\n", cut2[i+1]-cut2[i] );
            if (cut2[i + 1] - cut2[i] <= 0)
                fprintf(stderr, "### cut2[i+1]=%d, cut2[i]=%d\n", cut2[i + 1], cut2[i]);
            strncpy(tmpres2[j], seq2[j] + cut2[i], cut2[i + 1] - cut2[i]);
            tmpres2[j][cut2[i + 1] - cut2[i]] = 0;
        }
        if (kobetsubunkatsu && fftkeika)
            commongappick(clus2, tmpres2);  //dvtditr $B$K8F$P$l$?$H$-(B fftkeika=1
        //		if( kobetsubunkatsu ) commongappick( clus2, tmpres2 );

        if (constraint) {
            fprintf(stderr, "Not supported\n");
            exit(1);
        }
#if 0
		fprintf( stderr, "i=%d, before alignment", i );
		fprintf( stderr, "%4d\n", totallen );
		fprintf( stderr, "\n\n" );
		for( j=0; j<clus1; j++ ) 
		{
			fprintf( stderr, "%s\n", tmpres1[j] );
		}
		fprintf( stderr, "-------\n" );
		for( j=0; j<clus2; j++ ) 
		{
			fprintf( stderr, "%s\n", tmpres2[j] );
		}
#endif

#if 0
		fprintf( stdout, "writing input\n" );
		for( j=0; j<clus1; j++ )
		{
			fprintf( stdout, ">%d of GROUP1\n", j );
			fprintf( stdout, "%s\n", tmpres1[j] );
		}
		for( j=0; j<clus2; j++ )
		{
			fprintf( stdout, ">%d of GROUP2\n", j );
			fprintf( stdout, "%s\n", tmpres2[j] );
		}
		fflush( stdout );
#endif
        switch (alg) {
            case ('M'):
                if (scoringmatrices)  // called by tditeration.c
                    totalscore += MSalignmm_variousdist(ctx, scoringmatrices, tmpres1, tmpres2, eff1, eff2, eff1s, eff2s, clus1, clus2, alloclen, sgap1, sgap2, egap1, egap2, headgp, tailgp);
                else
                    totalscore += MSalignmm(ctx, n_dynamicmtx, tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen, sgap1, sgap2, egap1, egap2, headgp, tailgp, NULL, NULL, NULL, 0.0, 0.0);
                break;
            default:
                fprintf(stderr, "alg = %c\n", alg);
                ErrorExit("ERROR IN SOURCE FILE Falign.c");
                break;
        }

        nlen = strlen(tmpres1[0]);
        if (totallen + nlen > alloclen) {
            fprintf(stderr, "totallen=%d +  nlen=%d > alloclen = %d\n", totallen, nlen, alloclen);
            ErrorExit("LENGTH OVER in Falign\n ");
        }
        for (j = 0; j < clus1; j++)
            strcat(result1[j], tmpres1[j]);
        for (j = 0; j < clus2; j++)
            strcat(result2[j], tmpres2[j]);
        totallen += nlen;
#if 0
		fprintf( stderr, "i=%d", i );
		fprintf( stderr, "%4d\n", totallen );
		fprintf( stderr, "\n\n" );
		for( j=0; j<clus1; j++ ) 
		{
			fprintf( stderr, "%s\n", tmpres1[j] );
		}
		fprintf( stderr, "-------\n" );
		for( j=0; j<clus2; j++ ) 
		{
			fprintf( stderr, "%s\n", tmpres2[j] );
		}
#endif
    }

    //	reporterr( "\nafter Falign_udpari(), *fftlog = %d\n", *fftlog );

#if KEIKA
    fprintf(stderr, "DP ... done   \n");
#endif

    for (j = 0; j < clus1; j++)
        strcpy(seq1[j], result1[j]);
    for (j = 0; j < clus2; j++)
        strcpy(seq2[j], result2[j]);
#if 0
	for( j=0; j<clus1; j++ ) 
	{
		fprintf( stderr, "%s\n", result1[j] );
	}
	fprintf( stderr, "- - - - - - - - - - -\n" );
	for( j=0; j<clus2; j++ ) 
	{
		fprintf( stderr, "%s\n", result2[j] );
	}
#endif

    FreeCharMtx(result1);
    FreeCharMtx(result2);
    FreeCharMtx(tmpres1);
    FreeCharMtx(tmpres2);
    free(sgap1);
    free(egap1);
    free(sgap2);
    free(egap2);
    FreeCharMtx(tmpseq1);
    FreeCharMtx(tmpseq2);
    free(tmpptr1);
    free(tmpptr2);
#if RND
    FreeCharMtx(rndseq1);
    FreeCharMtx(rndseq2);
#endif
    return (totalscore);
}
