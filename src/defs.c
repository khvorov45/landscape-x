#include "mltaln.h"

int    scoreout = 0;
int    spscoreout = 0;
int    outnumber = 0;
int    legacygapcost = 0;
double minimumweight = 0.0005;
int    nwildcard = 0;

char* newgapstr = "-";

int nalphabets = 26;
int nscoredalphabets = 20;

double specificityconsideration = 0.0;
int    ndistclass = 10;
int    maxdistclass = -1;

int gmsg = 0;

double sueff_global = SUEFF;

double lenfaca, lenfacb, lenfacc, lenfacd;
int    maxl, tsize;

char codonpos = 0;
char codonscore = 0;

// for usetmpfile
int compacttree = 0;
int lhlimit = INT_MAX;
int specifictarget = 0;
int nadd = 0;  // <- static in tbfast.c, pairlocalalign.c
int usenaivescoreinsteadofalignmentscore = 0;
int nthreadreadlh = 1;
int LineLengthInFASTA = -1;
