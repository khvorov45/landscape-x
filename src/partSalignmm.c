#include "mltaln.h"

#define MACHIGAI 0
#define OUTGAP0TRY 1
#define DEBUG 0
#define XXXXXXX 0
#define USE_PENALTY_EX 0
#define FASTMATCHCALC 1

#if 0
static void st_OpeningGapCount( double *ogcp, int clus, char **seq, double *eff, int len )
{
	int i, j, gc, gb; 
	double feff;
	
	for( i=0; i<len; i++ ) ogcp[i] = 0.0;
	for( j=0; j<clus; j++ ) 
	{
		feff = (double)eff[j];
		gc = 0;
		for( i=0; i<len; i++ ) 
		{
			gb = gc;
			gc = ( seq[j][i] == '-' );
			{
				if( !gb *  gc ) ogcp[i] += feff;
			}
		}
	}
}

static void st_FinalGapCount( double *fgcp, int clus, char **seq, double *eff, int len )
{
	int i, j, gc, gb; 
	double feff;
	
	for( i=0; i<len; i++ ) fgcp[i] = 0.0;
	for( j=0; j<clus; j++ ) 
	{
		feff = (double)eff[j];
		gc = ( seq[j][0] == '-' );
		for( i=1; i<len+1; i++ ) 
		{
			gb = gc;
			gc = ( seq[j][i] == '-' );
			{
				if( gb * !gc ) fgcp[i-1] += feff;
			}
		}
	}
}
#endif

static int      impalloclen = 0;
static double** impmtx = NULL;
double
part_imp_match_out_sc(int i1, int j1) {
    //	fprintf( stderr, "impalloclen = %d\n", impalloclen );
    //	fprintf( stderr, "i1,j1=%d,%d -> impmtx=%f\n", i1, j1, impmtx[i1][j1] );
    return (impmtx[i1][j1]);
#if 0
	if( i1 == l1 || j1 == l2 ) return( 0.0 );
	return( impmtx[i1+start1][j1+start2] );
#endif
}

#if 1
void
part_imp_match_init_strict(int clus1, int clus2, int lgth1, int lgth2, char** seq1, char** seq2, double* eff1, double* eff2, double* eff1_kozo, double* eff2_kozo, LocalHom*** localhom, char* swaplist, int* orinum1, int* orinum2) {
    //	int i, j, k1, k2, tmpint, start1, start2, end1, end2;
    //	double effij;
    //	double effij_kozo;
    //	double effijx;
    //	char *pt, *pt1, *pt2;
    //	static char *nocount1 = NULL;
    //	static char *nocount2 = NULL;
    //	LocalHom *tmpptr;

    if (seq1 == NULL) {
        if (impmtx)
            FreeFloatMtx(impmtx);
        impmtx = NULL;
        //		if( nocount1 ) free( nocount1 );
        //		nocount1 = NULL;
        //		if( nocount2 ) free( nocount2 );
        //		nocount2 = NULL;

        return;
    }

    if (impalloclen < lgth1 + 2 || impalloclen < lgth2 + 2) {
        if (impmtx)
            FreeFloatMtx(impmtx);
        //		if( nocount1 ) free( nocount1 );
        //		if( nocount2 ) free( nocount2 );
        impalloclen = MAX(lgth1, lgth2) + 2;
        impmtx = AllocateFloatMtx(impalloclen, impalloclen);
        //		nocount1 = AllocateCharVec( impalloclen );
        //		nocount2 = AllocateCharVec( impalloclen );
    }

    fillimp(impmtx, clus1, clus2, lgth1, lgth2, seq1, seq2, eff1, eff2, eff1_kozo, eff2_kozo, localhom, swaplist, orinum1, orinum2);
}
#else
#endif
