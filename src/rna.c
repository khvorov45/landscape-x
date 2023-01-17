#include "mltaln.h"

#define MEMSAVE 1

#define DEBUG 1
#define USE_PENALTY_EX 1
#define STOREWM 1

#if 0
static double singleribosumscore( int n1, int n2, char **s1, char **s2, double *eff1, double *eff2, int p1, int p2 )
{
	double val;
	int i, j;
	int code1, code2;

	val = 0.0;
	for( i=0; i<n1; i++ ) for( j=0; j<n2; j++ )
	{
		code1 = amino_n[(int)s1[i][p1]];
		if( code1 > 3 ) code1 = 36;
		code2 = amino_n[(int)s2[j][p2]];
		if( code2 > 3 ) code2 = 36;

//		fprintf( stderr, "'l'%c-%c: %f\n", s1[i][p1], s2[j][p2], (double)ribosumdis[code1][code2] );

		val += (double)ribosumdis[code1][code2] * eff1[i] * eff2[j];
	}
	return( val );
}
static double pairedribosumscore53( int n1, int n2, char **s1, char **s2, double *eff1, double *eff2, int p1, int p2, int c1, int c2 )
{
	double val;
	int i, j;
	int code1o, code1u, code2o, code2u, code1, code2;

	val = 0.0;
	for( i=0; i<n1; i++ ) for( j=0; j<n2; j++ )
	{
		code1o = amino_n[(int)s1[i][p1]];
		code1u = amino_n[(int)s1[i][c1]];
		if( code1o > 3 ) code1 = code1o = 36;
		else if( code1u > 3 ) code1 = 36;
		else code1 = 4 + code1o * 4 + code1u;

		code2o = amino_n[(int)s2[j][p2]];
		code2u = amino_n[(int)s2[j][c2]];
		if( code2o > 3 ) code2 = code1o = 36;
		else if( code2u > 3 ) code2 = 36;
		else code2 = 4 + code2o * 4 + code2u;


//		fprintf( stderr, "%c%c-%c%c: %f\n", s1[i][p1], s1[i][c1], s2[j][p2], s2[j][c2], (double)ribosumdis[code1][code2] );

		if( code1 == 36 || code2 == 36 )
			val += (double)n_dis[code1o][code2o] * eff1[i] * eff2[j];
		else
			val += (double)ribosumdis[code1][code2] * eff1[i] * eff2[j];
	}
	return( val );
}

static double pairedribosumscore35( int n1, int n2, char **s1, char **s2, double *eff1, double *eff2, int p1, int p2, int c1, int c2 )
{
	double val;
	int i, j;
	int code1o, code1u, code2o, code2u, code1, code2;

	val = 0.0;
	for( i=0; i<n1; i++ ) for( j=0; j<n2; j++ )
	{
		code1o = amino_n[(int)s1[i][p1]];
		code1u = amino_n[(int)s1[i][c1]];
		if( code1o > 3 ) code1 = code1o = 36;
		else if( code1u > 3 ) code1 = 36;
		else code1 = 4 + code1u * 4 + code1o;

		code2o = amino_n[(int)s2[j][p2]];
		code2u = amino_n[(int)s2[j][c2]];
		if( code2o > 3 ) code2 = code1o = 36;
		else if( code2u > 3 ) code2 = 36;
		else code2 = 4 + code2u * 4 + code2o;


//		fprintf( stderr, "%c%c-%c%c: %f\n", s1[i][p1], s1[i][c1], s2[j][p2], s2[j][c2], (double)ribosumdis[code1][code2] );

		if( code1 == 36 || code2 == 36 )
			val += (double)n_dis[code1o][code2o] * eff1[i] * eff2[j];
		else
			val += (double)ribosumdis[code1][code2] * eff1[i] * eff2[j];
	}
	return( val );
}
#endif

