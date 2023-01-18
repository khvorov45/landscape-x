#include <stddef.h>
#include <stdarg.h>
#include <stdint.h>

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>
#include <unistd.h>

#include <sys/shm.h>  // shared memory
#include <sys/mman.h>  // shm_open

#include "align.h"

#define NKOUHO 20
#define NKOUHO_LONG 500

#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

#define FFT_THRESHOLD 80
#define FFT_WINSIZE_P 20
#define FFT_WINSIZE_D 100
#define DISPSEQF 60
#define DISPSITEI 0
#define MAXITERATION 500
#define N 5000000 /* nlen no saidaiti */
#define MAXSEG 100000
#define B 256
#define C 60 /*  1 gyou no mojisuu */
#define D 6
#define DFORMAT "%#6.3f"
#define rnd() ((1.0 / (RAND_MAX + 1.0)) * rand())
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define G(X) (((X) > (0)) ? (X) : (0))
#define BEFF 1.0 /* 0.6 ni suruto zureru */
#define WIN 3
#define SGAPP -1000
#define GETA2 0.001
#define GETA3 0.001
#define NOTSPECIFIED 100009
#define SUEFF 0.1 /* upg/(spg+upg)  -> sueff.sed */
#define DIVLOCAL 0
#define INTMTXSCALE 1000000.0
#define JTT 201
#define TM 202

#define BESTFIRST 0
#define BAATARI0 1
#define BAATARI1 2
#define BAATARI2 3

#define HAT3NODEBLOCK 500
#define MYBUFSIZE 1000 * 1000 * 100

typedef struct NodeInCub {
    int step;
    int LorR;
} NodeInCub;

typedef struct Node {
    struct Node* children[3];
    int          tmpChildren[3];
    double       length[3];
    double*      weightptr[3];
    int          top[3];
    int*         members[3];
} Node;

typedef struct Segment {
    int             start;
    int             end;
    int             center;
    double          score;
    int             skipForeward;
    int             skipBackward;
    struct Segment* pair;
    int             number;
} Segment;

typedef struct Segments {
    Segment group1;
    Segment group2;
    int     number1;
    int     number2;
} Segments;

typedef struct Bchain {
    struct Bchain* next;
    struct Bchain* prev;
    int            pos;
} Bchain;

typedef struct Achain {
    int next;
    int prev;
} Achain;

typedef struct Fukusosuu {
    double R;
    double I;
} Fukusosuu;

typedef struct Gappattern {
    int    len;
    double freq;
} Gappat;

typedef struct RNApair {
    int    uppos;
    double upscore;
    int    downpos;
    double downscore;
    int    bestpos;
    double bestscore;
} RNApair;

typedef struct Treedep {
    int    child0;
    int    child1;
    int    done;
    double distfromtip;
} Treedep;

typedef struct Addtree {
    int    nearest;
    double dist1;
    char*  neighbors;
    double dist2;
} Addtree;

typedef struct Lennum {
    int len;
    int num;
} Lennum;

typedef struct Pairnum {
    unsigned long long npairs;
    int                num;
    int                n0;
    int                n1;
} Pairnum;

typedef struct ExtAnch {
    int i;
    int j;
    int starti;
    int endi;
    int startj;
    int endj;
    int score;
} ExtAnch;

typedef struct GapPos {
    int pos;
    int len;
} GapPos;

void      MtxuntDouble(double**, int);
void      MtxmltDouble(double**, double**, int);
char*     AllocateCharVec(int);
char**    AllocateCharMtx(int, int);
void      ReallocateCharMtx(char**, int, int);
void      FreeCharMtx(char**);
double*   AllocateFloatVec(int);
void      FreeFloatVec(double*);
double**  AllocateFloatHalfMtx(int);
double**  AllocateFloatMtx(int, int);
void      FreeFloatHalfMtx(double**, int);
void      FreeFloatMtx(double**);
void      FreeFloatTri(double**);
int*      AllocateIntVec(int);
int*      AllocateIntVecLarge(unsigned long long);
void      FreeIntVec(int*);
int**     AllocateIntMtx(int, int);
void      FreeIntMtx(int**);
char***   AllocateCharCub(int, int, int);
void      FreeCharCub(char***);
int***    AllocateIntCub(int, int, int);
void      FreeIntCub(int***);
double*   AllocateDoubleVec(int);
void      FreeDoubleVec(double*);
double**  AllocateDoubleHalfMtx(int);
double**  AllocateDoubleMtx(int, int);
void      FreeDoubleHalfMtx(double**, int);
void      FreeDoubleMtx(double**);
double*** AllocateDoubleCub(int, int, int);
void      FreeDoubleCub(double***);
double*** AllocateFloatCub(int, int, int);
void      FreeFloatCub(double***);
short*    AllocateShortVec(int);
void      FreeShortVec(short*);
short**   AllocateShortMtx(int, int);
void      FreeShortMtx(short**);
void      freeintmtx(int**, int);

extern double A__align(aln_Context* ctx, double** scoringmtx, int penalty, int penalty_ex, char** seq1, char** seq2, double* eff1, double* eff2, int icyc, int jcyc, int alloclen, int constraint, double* impmatch, char* gs1, char* gs2, char* ge1, char* ge2, int headgp, int tailgp, int firstmem, int calledby, double*** cpmxchild0, double*** cpmxchild1, double*** cpmxresult, double orieff1, double orieff2);
extern double G__align11(aln_Context* ctx, char** seq1, char** seq2, int alloclen);
extern double G__align11_noalign(aln_Context* ctx, char** seq1, char** seq2);
extern void   cpmx_calc(aln_Context* ctx, char** seq, double** cpmx, double* eff, int lgth, int clus);
extern void   intergroup_score(char**, char**, double*, double*, int, int, int, double*);
extern int    fastconjuction(int* memlist, char** seq, char** aseq, double* peff, double* eff, char* d);

extern int       intlen(int* num);
extern void      veryfastsupg_double(int nseq, double** oeff, int*** topol, double** len);
extern void      counteff_simple_double_nostatic_memsave(int nseq, int*** topol, double** len, Treedep* dep, double* node);
extern void      calcimportance_half(int nseq, double* eff, char** seq, aln_LocalHom** localhom, int alloclen);
extern void      cpmx_calc_new(aln_Context* ctx, char** seq, double** cpmx, double* eff, int lgth, int clus);
extern void      cpmx_calc_add(aln_Context* ctx, char** seq, double** cpmx, double* eff, int lgth, int clus);
extern void      MScpmx_calc_new(aln_Context* ctx, char** seq, double** cpmx, double* eff, int lgth, int clus);
extern int       conjuctionforgaln(int s0, int s1, char** seq, char** aseq, double* peff, double* eff, char* d);
extern int       fastconjuction(int* memlist, char** seq, char** aseq, double* peff, double* eff, char* d);
extern int       fastconjuction_noname_kozo(int* memlist, char** seq, char** aseq, double* peff, double* eff, double* peff_kozo, double* eff_kozo, char* d);
extern int       fastconjuction_noname(int* memlist, char** seq, char** aseq, double* peff, double* eff, char* d, double mineff, double* oritotal);
extern int       fastconjuction_target(int* memlist, char** seq, char** aseq, double* peff, double* eff, char* d, double mineff, int* targetmap);
extern int       fastconjuction_noweight(int* memlist, char** seq, char** aseq, double* peff, char* d);
extern void      chardelete(char* seq, int d);
extern void      nodeFromABranch(int nseq, int* result, int** pairwisenode, int*** topol, int step, int num);
extern void      OneClusterAndTheOther_fast(int locnjob, int* memlist1, int* memlist2, int* s1, int* s2, char* pairbuf, int*** topol, int step, int branch, double** smalldistmtx, double* distontree);
extern void      makeEffMtx(int nseq, double** mtx, double* vec);
extern int       msshrinklocalhom_fast(int* memlist1, int* memlist2, aln_LocalHom** localhom, aln_LocalHom*** localhomshrink);
extern int       msshrinklocalhom_fast_half(int* memlist1, int* memlist2, aln_LocalHom** localhom, aln_LocalHom*** localhomshrink);
extern int       msshrinklocalhom_fast_target(int* memlist1, int* memlist2, aln_LocalHom** localhom, aln_LocalHom*** localhomshrink, char* swaplist, int* targetmap);
extern int       fastshrinklocalhom(int* mem1, int* mem2, aln_LocalHom** localhom, aln_LocalHom*** localhomshrink);
extern int       fastshrinklocalhom_half(int* mem1, int* mem2, aln_LocalHom** localhom, aln_LocalHom*** localhomshrink);
extern int       fastshrinklocalhom_target(int* mem1, int* mem2, aln_LocalHom** localhom, aln_LocalHom*** localhomshrink, char* swaplist, int* targetmap);
extern int       fastshrinklocalhom_one(int* mem1, int* mem2, int norg, aln_LocalHom** localhom, aln_LocalHom*** localhomshrink);
extern int       msfastshrinklocalhom(int* mem1, int* mem2, aln_LocalHom** localhom, aln_LocalHom*** localhomshrink);
extern int       fastshrinklocalhom_half_seed(int* mem1, int* mem2, int nseed, int* posinlsh1, int* posinlsh2, aln_LocalHom** localhom, aln_LocalHom*** localhomshrink);
extern void      checkMinusLength(int nseq, double** len);
extern void      negativeMember2(int* mem, int* query, int locnseq);
extern int*      negativeMember(int* query, int locnseq);
extern int       IntExistsInVec(int query, int* vector);
extern NodeInCub searchParent(int top, int*** topol, int Start, int End);
extern void      stopolInit(int n, Node* stopol);
extern void      treeCnv(Node* stopol, int locnseq, int*** topol, double** len, double** bw);
extern int       isLeaf(Node node);
extern double    syntheticLength(Node* ob, Node* oppositeNode);
extern double    calcW(Node* ob, Node* op);
extern void      branchWeightToPairWeight(int locnseq, int*** topol, double** pw, double** bw);
extern void      weightFromABranch(int nseq, double* result, Node* stopol, int*** topol, int step, int LorR);
extern void      distFromABranch(int nseq, double* result, Node* stopol, int*** topol, double** len, int step, int LorR);

extern void gapireru(char* res, char* ori, char* gt);
extern int  seqlen(char* seq, char gap);
extern void st_FinalGapCount(double* fgcp, int clus, char** seq, double* eff, int len);
extern void st_FinalGapAdd(double* fgcp, int clus, char** seq, double* eff, int len);
extern void st_OpeningGapCount(double* ogcp, int clus, char** seq, double* eff, int len);
extern void st_OpeningGapAdd(double* ogcp, int clus, char** seq, double* eff, int len);
extern void new_FinalGapCount(double* fgcp, int clus, char** seq, double* eff, int len, char* g);
extern void new_OpeningGapCount(double* ogcp, int clus, char** seq, double* eff, int len, char* g);
extern void getGapPattern(double* fgcp, int clus, char** seq, double* eff, int len);
extern void getkyokaigap(char* g, char** s, int pos, int n);
void        makegrouprna(RNApair*** group, RNApair*** all, int* memlist);
extern void fixed_musclesupg_double_realloc_nobk_halfmtx_memsave(aln_Context* ctx, int nseq, double** eff, int*** topol, double** len, Treedep*, int efffree);
extern void imp_match_init_strict(aln_Context* ctx, int clus1, int clus2, int lgth1, int lgth2, char** seq1, char** seq2, double* eff1, double* eff2, aln_LocalHom*** localhom, char* swaplist, int* memlist1, int* memlist2);
extern void cpmx_ribosum(aln_Context* ctx, char** seq, char** seqr, char* dir, double** cpmx, double* eff, int lgth, int clus);
extern void assignstrweight(int nseq, double* strweight, Node* stopol, int*** topol, int step, int LorR, char* kozoari, double* seqweight);
extern void intcpy(int* s1, int* s2);
extern void intncpy(int* s1, int* s2, int n);
extern void fltncpy(double* s1, double* s2, int n);
extern void intcat(int* s1, int* s2);
extern void gapcountf(double* freq, char** seq, int nseq, double* eff, int lgth);
extern void gapcountadd(double* freq, char** seq, int nseq, double* eff, int lgth);
extern void outgapcount(double* freq, int nseq, char* gappat, double* eff);
extern void makedynamicmtx(aln_Context* ctx, double** out, double offset);

extern int* topolorderz(int* order, int*** topol, Treedep* dep, int pos, int nchild);
extern void fillimp(aln_Context* ctx, double** impmtx, int clus1, int clus2, int lgth1, int lgth2, char** seq1, char** seq2, double* eff1, double* eff2, aln_LocalHom*** localhom, char* swaplist, int* orinum1, int* orinum2);
