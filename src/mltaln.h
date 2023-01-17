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

extern void reporterr(const char* str, ...);
#define aln_assertAction() \
    do { \
        reporterr("assertion failure at %s:%d\n", __FILE__, __LINE__); \
        aln_debugbreak(); \
        _exit(1); \
    } while (0)

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

typedef struct LocalHom {
    struct LocalHom* next;
    struct LocalHom* last;
    int              start1;
    int              end1;
    int              start2;
    int              end2;
    double           opt;
    int              overlapaa;
    int              extended;
    double           importance;
    double           rimportance;
    char             korh;
    int              nokori;
} LocalHom;

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

typedef struct Context {
    char     modelname[500];
    int      njob;
    int      maxInputSeqLen;
    int      amino_n[0x100];
    char     amino_grp[0x100];
    int**    amino_dis;
    double** amino_dis_consweight_multi;
    int**    n_dis;
    double** n_disLN;
    int**    n_disFFT;
    double** n_dis_consweight_multi;
    uint8_t  amino[0x100];
    double   polarity[0x100];
    double   volume[0x100];
    int      ribosumdis[37][37];
    int      ppid;
    int      pslocal;
    int      divpairscore;
    int      fmodel;
    int      kobetsubunkatsu;
    int      dorp;
    int      weight;
    int      utree;
    int      tbutree;
    int      trywarp;
    int      penalty;
    int      penaltyLN;
    int      penalty_dist;
    int      RNApenalty;
    int      RNAppenalty;
    int      RNApenalty_ex;
    int      RNAppenalty_ex;
    int      penalty_ex;
    int      ppenalty_ex;
    int      penalty_exLN;
    int      penalty_EX;
    int      ppenalty_EX;
    int      penalty_OP;
    int      ppenalty_OP;
    int      penalty_shift;
    int      offset;
    int      offsetLN;
    int      offsetFFT;
    int      RNAthr;
    int      RNApthr;
    int      TMorJTT;
    char     force_fft;
    int      nevermemsave;
    int      fftscore;
    int      fftWinSize;
    int      fftThreshold;
    int      fftRepeatStop;
    int      fftNoAnchStop;
    int      divWinSize;
    int      divThreshold;
    int      disp;
    int      outgap;
    char     alg;
    int      cnst;
    int      mix;
    int      tbitr;
    int      tbweight;
    int      tbrweight;
    int      disopt;
    int      pamN;
    int      checkC;
    double   geta2;
    int      kimuraR;
    char*    swopt;
    int      fftkeika;
    int      score_check;
    char*    addfile;
    int      addprofile;
    double   consweight_multi;
    double   consweight_rna;
    char     RNAscoremtx;
    char*    signalSM;
    char**   seq_g;
    char**   res_g;
    int      rnakozo;
    char     rnaprediction;
    int      commonAlloc1;
    int      commonAlloc2;
    int**    commonIP;
    int**    commonJP;
    int      randomseed;
    int      parallelizationstrategy;
    int      scoreout;
    int      spscoreout;
    int      outnumber;
    int      legacygapcost;
    int      nwildcard;
    char*    newgapstr;
    int      nalphabets;
    int      nscoredalphabets;
    int      ndistclass;
    int      maxdistclass;
    double   lenfaca;
    double   lenfacb;
    double   lenfacc;
    double   lenfacd;
    int      maxl;
    int      tsize;
    char     codonpos;
    char     codonscore;
    int      compacttree;
    int      lhlimit;
    int      nadd;
} Context;

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
int**     AllocateIntMtxLarge(unsigned long long, unsigned long long);
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

extern void   constants(aln_Opts opts, Context* ctx, int nseq, char** seq);
extern double A__align(aln_Opts opts, Context* ctx, double** scoringmtx, int penalty, int penalty_ex, char** seq1, char** seq2, double* eff1, double* eff2, int icyc, int jcyc, int alloclen, int constraint, double* impmatch, char* gs1, char* gs2, char* ge1, char* ge2, int headgp, int tailgp, int firstmem, int calledby, double*** cpmxchild0, double*** cpmxchild1, double*** cpmxresult, double orieff1, double orieff2);
extern double G__align11(Context* ctx, double** scoringmtx, char** seq1, char** seq2, int alloclen, int headgp, int tailgp);
extern double imp_match_out_sc(int, int);
extern void   ErrorExit(char* message);
extern void   cpmx_calc(Context* ctx, char** seq, double** cpmx, double* eff, int lgth, int clus);
extern void   intergroup_score(char**, char**, double*, double*, int, int, int, double*);
extern int    fastconjuction(int* memlist, char** seq, char** aseq, double* peff, double* eff, char* d);

extern int       intlen(int* num);
extern void      veryfastsupg_double(int nseq, double** oeff, int*** topol, double** len);
extern void      counteff_simple(int nseq, int*** topol, double** len, double* node);
extern void      counteff_simple_double_nostatic(int nseq, int*** topol, double** len, double* node);
extern void      counteff_simple_double_nostatic_memsave(int nseq, int*** topol, double** len, Treedep* dep, double* node);
extern void      copyWithNoGaps(char* aseq, const char* seq);
extern void      gappick(int nseq, int s, char** aseq, char** mseq2, double** eff, double* effarr);
extern void      commongappick_record(int nseq, char** seq, int* map);
extern void      doublencpy(double* vec1, double* vec2, int len);
extern void      calcimportance_target(Context* ctx, int nseq, int ntarget, double* eff, char** seq, LocalHom** localhom, int* targetmap, int* targetmapr, int alloclen);
extern void      dontcalcimportance_firstone(int nseq, LocalHom** localhom);
extern void      dontcalcimportance_half(Context* ctx, int nseq, char** seq, LocalHom** localhom);
extern void      calcimportance_half(Context* ctx, int nseq, double* eff, char** seq, LocalHom** localhom, int alloclen);
extern void      cpmx_calc_new(Context* ctx, char** seq, double** cpmx, double* eff, int lgth, int clus);
extern void      cpmx_calc_add(Context* ctx, char** seq, double** cpmx, double* eff, int lgth, int clus);
extern void      MScpmx_calc_new(Context* ctx, char** seq, double** cpmx, double* eff, int lgth, int clus);
extern int       conjuctionforgaln(int s0, int s1, char** seq, char** aseq, double* peff, double* eff, char* d);
extern int       fastconjuction(int* memlist, char** seq, char** aseq, double* peff, double* eff, char* d);
extern int       fastconjuction_noname_kozo(int* memlist, char** seq, char** aseq, double* peff, double* eff, double* peff_kozo, double* eff_kozo, char* d);
extern int       fastconjuction_noname(int* memlist, char** seq, char** aseq, double* peff, double* eff, char* d, double mineff, double* oritotal);
extern int       fastconjuction_target(int* memlist, char** seq, char** aseq, double* peff, double* eff, char* d, double mineff, int* targetmap);
extern int       fastconjuction_noweight(int* memlist, char** seq, char** aseq, double* peff, char* d);
extern void      chardelete(char* seq, int d);
extern int       RootBranchNode(int nseq, int*** topol, int step, int branch);
extern void      BranchLeafNode(int nseq, int*** topol, int* node, int step, int branch);
extern void      RootLeafNode(int nseq, int*** topol, int* node);
extern void      nodeFromABranch(int nseq, int* result, int** pairwisenode, int*** topol, int step, int num);
extern void      OneClusterAndTheOther_fast(int locnjob, int* memlist1, int* memlist2, int* s1, int* s2, char* pairbuf, int*** topol, int step, int branch, double** smalldistmtx, double* distontree);
extern void      makeEffMtx(int nseq, double** mtx, double* vec);
extern int       msshrinklocalhom_fast(int* memlist1, int* memlist2, LocalHom** localhom, LocalHom*** localhomshrink);
extern int       msshrinklocalhom_fast_half(int* memlist1, int* memlist2, LocalHom** localhom, LocalHom*** localhomshrink);
extern int       msshrinklocalhom_fast_target(int* memlist1, int* memlist2, LocalHom** localhom, LocalHom*** localhomshrink, char* swaplist, int* targetmap);
extern int       fastshrinklocalhom(int* mem1, int* mem2, LocalHom** localhom, LocalHom*** localhomshrink);
extern int       fastshrinklocalhom_half(int* mem1, int* mem2, LocalHom** localhom, LocalHom*** localhomshrink);
extern int       fastshrinklocalhom_target(int* mem1, int* mem2, LocalHom** localhom, LocalHom*** localhomshrink, char* swaplist, int* targetmap);
extern int       fastshrinklocalhom_one(int* mem1, int* mem2, int norg, LocalHom** localhom, LocalHom*** localhomshrink);
extern int       msfastshrinklocalhom(int* mem1, int* mem2, LocalHom** localhom, LocalHom*** localhomshrink);
extern int       fastshrinklocalhom_half_seed(int* mem1, int* mem2, int nseed, int* posinlsh1, int* posinlsh2, LocalHom** localhom, LocalHom*** localhomshrink);
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
extern void      calcBranchWeight(double** bw, int locnseq, Node* stopol, int*** topol);
extern void      branchWeightToPairWeight(int locnseq, int*** topol, double** pw, double** bw);
extern void      weightFromABranch(int nseq, double* result, Node* stopol, int*** topol, int step, int LorR);
extern void      distFromABranch(int nseq, double* result, Node* stopol, int*** topol, double** len, int step, int LorR);
extern double    imp_match_out_scD(int i1, int j1);
extern double    MSalignmm(Context* ctx, double** n_dynamicmtx, char** seq1, char** seq2, double* eff1, double* eff2, int icyc, int jcyc, int alloclen, char*, char*, char*, char*, int headgp, int tailgp, double*** cpmxchild0, double*** cpmxchild1, double*** cpmxresult, double orieff1, double orieff2);
extern double    MSalignmm_variousdist(Context* ctx, double*** matrices, char** seq1, char** seq2, double* eff1, double* eff2, double** eff1s, double** eff2s, int icyc, int jcyc, int alloclen, char*, char*, char*, char*, int headgp, int tailgp);
extern double    Lalignmm_hmout(Context* ctx, char** seq1, char** seq2, double* eff1, double* eff2, int icyc, int jcyc, char*, char*, char*, double** map);
extern double    MSalign11(Context* ctx, char** seq1, char** seq2, int alloclen);
extern double    G__align11_noalign(Context* ctx, double** scoringmtx, int penal, int penal_ex, char** seq1, char** seq2);
extern void      JTTmtx(double** rsr, double* freq, unsigned char locamino[0x80], char locgrp[0x80], int isTM);
extern void      BLOSUMmtx(int n, double** matrix, double* freq, unsigned char* amino, char* amino_grp);
extern void      putlocalhom(char* al1, char* al2, LocalHom* localhompt, int off1, int off2, int opt, int overlapaa, char korh);
extern void      ErrorExit(char* message);
extern void      initSignalSM(Context* ctx);
extern void      readlocalhomtable_half(FILE* fp, int njob, LocalHom** localhomtable, char* kozoarivec);
extern void      readlocalhomtable_target(FILE* fp, int nt, int njob, LocalHom** localhomtable, char* kozoarivec, int* targetmap);
extern void      freelocalhom1(LocalHom* lh);
extern void      initlocalhom1(LocalHom* lh);

extern void   gapireru(Context* ctx, char* res, char* ori, char* gt);
extern int    seqlen(Context* ctx, char* seq);
extern void   st_FinalGapCount(double* fgcp, int clus, char** seq, double* eff, int len);
extern void   st_FinalGapAdd(double* fgcp, int clus, char** seq, double* eff, int len);
extern void   st_OpeningGapCount(double* ogcp, int clus, char** seq, double* eff, int len);
extern void   st_OpeningGapAdd(double* ogcp, int clus, char** seq, double* eff, int len);
extern void   st_FinalGapCount_zure(double* fgcp, int clus, char** seq, double* eff, int len);
extern void   getdiaminofreq_x(double* freq, int clus, char** seq, double* eff, int len);
extern void   new_FinalGapCount_zure(double* fgcp, int clus, char** seq, double* eff, int len, char* s, char* e);
extern void   new_FinalGapCount(double* fgcp, int clus, char** seq, double* eff, int len, char* g);
extern void   new_OpeningGapCount(double* ogcp, int clus, char** seq, double* eff, int len, char* g);
extern void   new_OpeningGapCount_zure(double* ogcp, int clus, char** seq, double* eff, int len, char* s, char* e);
extern void   getGapPattern(double* fgcp, int clus, char** seq, double* eff, int len);
extern void   getgapfreq(double* freq, int clus, char** seq, double* eff, int len);
extern void   getgapfreq_zure(double* freq, int clus, char** seq, double* eff, int len);
extern void   getgapfreq_zure_part(double* freq, int clus, char** seq, double* eff, int len, char* s);
extern void   getdiaminofreq_part(double* freq, int clus, char** seq, double* eff, int len, char* s, char* e);
extern void   getdigapfreq_part(double* freq, int clus, char** seq, double* eff, int len, char* s, char* e);
extern void   getdiaminofreq_st(double* freq, int clus, char** seq, double* eff, int len);
extern void   getdigapfreq_st(double* freq, int clus, char** seq, double* eff, int len);
extern void   st_getGapPattern(Gappat** gpat, int clus, char** seq, double* eff, int len);
extern void   getkyokaigap(char* g, char** s, int pos, int n);
extern double naivepairscore11(Context* ctx, const char* seq1, const char* seq2, int penal);
extern double naivepairscore11_dynmtx(Context* ctx, double**, char* seq1, char* seq2, int penal);
extern double naivepairscorefast(Context* ctx, const char* seq1, const char* seq2, int* skip1, int* skip2, int penal);
void          makegrouprna(RNApair*** group, RNApair*** all, int* memlist);
extern void   fixed_musclesupg_double_realloc_nobk_halfmtx_memsave(aln_Opts opts, Context* ctx, int nseq, double** eff, int*** topol, double** len, Treedep*, int progressout, int efffree);
extern void   imp_match_init_strict(aln_Opts opts, Context* ctx, int clus1, int clus2, int lgth1, int lgth2, char** seq1, char** seq2, double* eff1, double* eff2, double* eff1kozo, double* eff2kozo, LocalHom*** localhom, char* swaplist, int* memlist1, int* memlist2, int* uselh, int* seedinlh1, int* seedinlh2, int nodeid, int nfiles);
extern void   cpmx_ribosum(Context* ctx, char** seq, char** seqr, char* dir, double** cpmx, double* eff, int lgth, int clus);
extern void   assignstrweight(int nseq, double* strweight, Node* stopol, int*** topol, int step, int LorR, char* kozoari, double* seqweight);
extern double plainscore(Context* ctx, int nseq, char** seq);
extern void   eq2dash(char* s);
extern void   eq2dashmatometehayaku(char** s, int n);
extern void   findnewgaps(Context* ctx, int rep, char** seq, int* gaplen);
extern void   findcommongaps(int, char**, int*);
extern void   adjustgapmap(int, int*, char*);
extern void   restorecommongaps(int n, int n0, char** seq, int* top0, int* top1, int* gaplen, int alloclen, char gapchar);
extern void   restorecommongapssmoothly(int n, int n0, char** seq, int* top0, int* top1, int* gaplen, int alloclen, char gapchar);
extern int    samemember(int* mem, int* cand);
extern int    samemembern(int* mem, int* cand, int candn);
extern int    includemember(int* mem, int* cand);
extern int    overlapmember(int* mem1, int* mem2);
extern void   intcpy(int* s1, int* s2);
extern void   intncpy(int* s1, int* s2, int n);
extern void   fltncpy(double* s1, double* s2, int n);
extern void   intcat(int* s1, int* s2);
extern void   gapcount(double* freq, char** seq, int nseq, double* eff, int lgth);
extern void   gapcountf(double* freq, char** seq, int nseq, double* eff, int lgth);
extern void   gapcountadd(double* freq, char** seq, int nseq, double* eff, int lgth);
extern void   outgapcount(double* freq, int nseq, char* gappat, double* eff);
extern void   makedynamicmtx(Context* ctx, double** out, double** in, double offset);
extern void   freeconstants(Context* ctx);
extern void   makeskiptable(int n, int** skip, char** seq);
extern int    generatesubalignmentstable(int nseq, int*** tablept, int* nsubpt, int* maxmempt, int*** topol, double** len, double threshold);
extern double sumofpairsscore(Context* ctx, int nseq, char** seq);

extern int    isallgap(char*);
extern int    deletenewinsertions_whole(int on, int an, char** oseq, char** aseq, GapPos** deletelist);
extern int    deletenewinsertions_whole_eq(int on, int an, char** oseq, char** aseq, GapPos** deletelist);
extern int    deletenewinsertions_difflist(int on, int an, char** oseq, char** aseq, GapPos** difflist);
extern int    recordoriginalgaps(char* originallygapped, int n, char** s);
extern void   restoreoriginalgaps(int n, char** seq, char* originalgaps);
extern void   reconstructdeletemap(Context* ctx, int nadd, char** addbk, GapPos** deletelist, char** realn, FILE* fp, const char* const* name);
extern void   reconstructdeletemap_compact(Context* ctx, int nadd, char** addbk, GapPos** deletelist, char** realn, FILE* fp, const char* const* name);
extern void   topolorder(int* order, int* posinorder, int*** topol, Treedep* dep, int pos, int child);
extern int*   topolorderz(int* order, int*** topol, Treedep* dep, int pos, int nchild);
extern int*   topolordery(int* order, int*** topol, Treedep* dep, int pos, int nchild);
extern void   fillimp(aln_Opts opts, double** impmtx, int clus1, int clus2, int lgth1, int lgth2, char** seq1, char** seq2, double* eff1, double* eff2, double* eff1_kozo, double* eff2_kozo, LocalHom*** localhom, char* swaplist, int* orinum1, int* orinum2);
extern void   fillimp_file(aln_Opts opts, Context* ctx, double** impmtx, int clus1, int clus2, int lgth1, int lgth2, char** seq1, char** seq2, double* eff1, double* eff2, double* eff1_kozo, double* eff2_kozo, LocalHom*** localhom, int* orinum1, int* orinum2, int* uselh, int* seedinlh1, int* seedinlh2, int nodeid, int nfiles);
extern int    pairlocalalign(aln_Opts opts, Context* ctx, char** seqgui);
extern void   use_getrusage(void);
extern void   sortbylength(int* uselh, Lennum* in, int size, unsigned long long numpairs);
extern void   limitlh(int* uselh, Lennum* in, int size, int limit);
