#include "mltaln.h"

// for usetmpfile
int compacttree = 0;
int lhlimit = INT_MAX;
int specifictarget = 0;
int nadd = 0;  // <- static in tbfast.c, pairlocalalign.c
int usenaivescoreinsteadofalignmentscore = 0;
int nthreadreadlh = 1;
int LineLengthInFASTA = -1;
