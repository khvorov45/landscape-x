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

#include <sys/resource.h>  // for setstacksize, 2016/Jun
#include <sys/shm.h>  // shared memory
#include <sys/mman.h>  // shm_open

// TODO(sen) Figure out asserts in this context
// #define aln_assertAction() \
//     do { \
//         reporterr("assertion failure at %s:%d\n", __FILE__, __LINE__); \
//         aln_debugbreak(); \
//         _exit(1); \
//     } while (0)

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
    FILE*    prep_g;
    FILE*    trap_g;
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
    int      specifictarget;
    int      nadd;
    int      usenaivescoreinsteadofalignmentscore;
    int      LineLengthInFASTA;
} Context;

void      MtxuntDouble(double**, int);
void      MtxmltDouble(double**, double**, int);
char*     AllocateCharVec(int);
void      FreeCharVec(char*);
char**    AllocateCharMtx(int, int);
void      ReallocateCharMtx(char**, int, int);
void      FreeCharMtx(char**);
double*   AllocateFloatVec(int);
void      FreeFloatVec(double*);
double**  AllocateFloatHalfMtx(int);
double**  AllocateFloatMtx(int, int);
void      FreeFloatHalfMtx(double**, int);
void      FreeFloatMtx(double**);
double**  AlocateFloatTri(int);
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
extern char** Calignm1(double* wm, char** aseq, char* seq, double* effarr, int icyc, int ex);
extern char** Dalignm1();
extern char** align0();
extern double Cscore_m_1(char**, int, int, double**);
extern double score_m_1(char**, int, int, double**);
extern char   seqcheck(Context* ctx, char**);
extern double ipower(double, int);
extern double translate_and_Calign(char** mseq1, char** mseq2, double* effarr1, double* effarr2, int clus1, int clus2, int alloclen);
extern double A__align(aln_Opts opts, Context* ctx, double** scoringmtx, int penalty, int penalty_ex, char** seq1, char** seq2, double* eff1, double* eff2, int icyc, int jcyc, int alloclen, int constraint, double* impmatch, char* gs1, char* gs2, char* ge1, char* ge2, int headgp, int tailgp, int firstmem, int calledby, double*** cpmxchild0, double*** cpmxchild1, double*** cpmxresult, double orieff1, double orieff2);
extern double A__align11();
extern double A__align_gapmap(void);
extern double L__align11(Context* ctx, double** scoringmtx, double scoreoffset, char** seq1, char** seq2, int alloclen, int* off1pt, int* off2pt);
extern double G__align11(Context* ctx, double** scoringmtx, char** seq1, char** seq2, int alloclen, int headgp, int tailgp);
extern double Falign(aln_Opts opts, Context* ctx, int** whichmtx, double*** scoringmatrices, double** scoreingmtx, char** seq1, char** seq2, double* eff1, double* eff2, double** eff1s, double** eff2s, int clus1, int clus2, int alloclen, int* fftlog);
extern double Conalign();
extern double Aalign(Context* ctx, char** seq1, char** seq2, double* eff1, double* eff2, int icyc, int jcyc, int alloclen);
extern double imp_match_out_sc(int, int);
extern void   ErrorExit(char* message);
extern void   cpmx_calc(Context* ctx, char** seq, double** cpmx, double* eff, int lgth, int clus);
extern void   intergroup_score(char**, char**, double*, double*, int, int, int, double*);
extern int    conjuctionfortbfast();
extern int    fastconjuction(int* memlist, char** seq, char** aseq, double* peff, double* eff, char* d);

extern int          intlen(int* num);
extern void         exitall(char arr[]);
extern void         display(Context* ctx, char** seq, int nseq);
extern void         intergroup_score_new(char** seq1, char** seq2, double* eff1, double* eff2, int clus1, int clus2, int len, double* value);
extern void         veryfastsupg_int_realloc_nobk(int njob, int** mtx, int*** topol, double** len);
extern void         veryfastsupg_double(int nseq, double** oeff, int*** topol, double** len);
extern void         veryfastsupg_int(int nseq, int** oeff, int*** topol, double** len);
extern double       ipower(double x, int n);
extern void         counteff_simple(int nseq, int*** topol, double** len, double* node);
extern void         counteff_simple_double_nostatic(int nseq, int*** topol, double** len, double* node);
extern void         counteff_simple_double_nostatic_memsave(int nseq, int*** topol, double** len, Treedep* dep, double* node);
extern void         gappick_samestring(char* aseq);
extern void         gappick0(char* aseq, const char* seq);
extern void         gappick(int nseq, int s, char** aseq, char** mseq2, double** eff, double* effarr);
extern void         commongappick_record(int nseq, char** seq, int* map);
extern void         commongappick(int nseq, char** seq);
extern int          commongapcount(int, int, char**, char**);
extern void         strins(char* str1, char* str2);
extern int          isaligned(int nseq, char** seq);
extern void         doublencpy(double* vec1, double* vec2, int len);
extern char*        progName(char* str);
extern void         calcimportance_target(Context* ctx, int nseq, int ntarget, double* eff, char** seq, LocalHom** localhom, int* targetmap, int* targetmapr, int alloclen);
extern void         dontcalcimportance_lastone(int nseq, double* eff, char** seq, LocalHom** localhom);
extern void         dontcalcimportance_firstone(int nseq, LocalHom** localhom);
extern void         dontcalcimportance_half(Context* ctx, int nseq, char** seq, LocalHom** localhom);
extern void         calcimportance_half(Context* ctx, int nseq, double* eff, char** seq, LocalHom** localhom, int alloclen);
extern void         weightimportance2(int nseq, double* eff, LocalHom** localhom);
extern void         weightimportance4(int clus1, int clus2, double* eff1, double* eff2, LocalHom*** localhom);
extern void         extendlocalhom(int nseq, LocalHom** localhom);
extern void         cpmx_calc_new(Context* ctx, char** seq, double** cpmx, double* eff, int lgth, int clus);
extern void         cpmx_calc_add(Context* ctx, char** seq, double** cpmx, double* eff, int lgth, int clus);
extern void         MScpmx_calc_new(Context* ctx, char** seq, double** cpmx, double* eff, int lgth, int clus);
extern void         strnbcat(char* s1, char* s2, int m);
extern int          conjuctionforgaln(int s0, int s1, char** seq, char** aseq, double* peff, double* eff, char* d);
extern int          fastconjuction(int* memlist, char** seq, char** aseq, double* peff, double* eff, char* d);
extern int          fastconjuction_noname_kozo(int* memlist, char** seq, char** aseq, double* peff, double* eff, double* peff_kozo, double* eff_kozo, char* d);
extern int          fastconjuction_noname(int* memlist, char** seq, char** aseq, double* peff, double* eff, char* d, double mineff, double* oritotal);
extern int          fastconjuction_target(int* memlist, char** seq, char** aseq, double* peff, double* eff, char* d, double mineff, int* targetmap);
extern int          fastconjuction_noweight(int* memlist, char** seq, char** aseq, double* peff, char* d);
extern void         chardelete(char* seq, int d);
extern int          RootBranchNode(int nseq, int*** topol, int step, int branch);
extern void         BranchLeafNode(int nseq, int*** topol, int* node, int step, int branch);
extern void         RootLeafNode(int nseq, int*** topol, int* node);
extern void         nodeFromABranch(int nseq, int* result, int** pairwisenode, int*** topol, int step, int num);
extern void         OneClusterAndTheOther_fast(int locnjob, int* memlist1, int* memlist2, int* s1, int* s2, char* pairbuf, int*** topol, int step, int branch, double** smalldistmtx, double* distontree);
extern void         makeEffMtx(int nseq, double** mtx, double* vec);
extern int          msshrinklocalhom_fast(int* memlist1, int* memlist2, LocalHom** localhom, LocalHom*** localhomshrink);
extern int          msshrinklocalhom_fast_half(int* memlist1, int* memlist2, LocalHom** localhom, LocalHom*** localhomshrink);
extern int          msshrinklocalhom_fast_target(int* memlist1, int* memlist2, LocalHom** localhom, LocalHom*** localhomshrink, char* swaplist, int* targetmap);
extern int          fastshrinklocalhom(int* mem1, int* mem2, LocalHom** localhom, LocalHom*** localhomshrink);
extern int          fastshrinklocalhom_half(int* mem1, int* mem2, LocalHom** localhom, LocalHom*** localhomshrink);
extern int          fastshrinklocalhom_target(int* mem1, int* mem2, LocalHom** localhom, LocalHom*** localhomshrink, char* swaplist, int* targetmap);
extern int          fastshrinklocalhom_one(int* mem1, int* mem2, int norg, LocalHom** localhom, LocalHom*** localhomshrink);
extern int          msfastshrinklocalhom(int* mem1, int* mem2, LocalHom** localhom, LocalHom*** localhomshrink);
extern int          fastshrinklocalhom_half_seed(int* mem1, int* mem2, int nseed, int* posinlsh1, int* posinlsh2, LocalHom** localhom, LocalHom*** localhomshrink);
extern void         checkMinusLength(int nseq, double** len);
extern void         negativeMember2(int* mem, int* query, int locnseq);
extern int*         negativeMember(int* query, int locnseq);
extern int          IntExistsInVec(int query, int* vector);
extern NodeInCub    searchParent(int top, int*** topol, int Start, int End);
extern void         stopolInit(int n, Node* stopol);
extern void         treeCnv(Node* stopol, int locnseq, int*** topol, double** len, double** bw);
extern int          isLeaf(Node node);
extern double       syntheticLength(Node* ob, Node* oppositeNode);
extern double       calcW(Node* ob, Node* op);
extern void         calcBranchWeight(double** bw, int locnseq, Node* stopol, int*** topol);
extern void         branchWeightToPairWeight(int locnseq, int*** topol, double** pw, double** bw);
extern void         weightFromABranch(int nseq, double* result, Node* stopol, int*** topol, int step, int LorR);
extern void         distFromABranch(int nseq, double* result, Node* stopol, int*** topol, double** len, int step, int LorR);
extern void         keika(char* str, int current, int all);
extern double       maxItch(double* soukan, int size);
extern void         calcNaiseki(Fukusosuu* value, Fukusosuu* x, Fukusosuu* y);
extern Fukusosuu*   AllocateFukusosuuVec(int l1);
extern Fukusosuu**  AllocateFukusosuuMtx(int l1, int l2);
extern Fukusosuu*** AllocateFukusosuuCub(int l1, int l2, int l3);
extern void         FreeFukusosuuVec(Fukusosuu* vec);
extern void         FreeFukusosuuMtx(Fukusosuu** mtx);
extern int          getKouho(int* kouho, int nkouho, double* soukan, int nlen2);
extern void         zurasu2(int lag, int clus1, int clus2, char** seq1, char** seq2, char** aseq1, char** aseq2);
extern void         zurasu(int lag, int clus1, int clus2, char** seq1, char** seq2, char** aseq1, char** aseq2);
extern int          alignableReagion(Context* ctx, int clus1, int clus2, char** seq1, char** seq2, double* eff1, double* eff2, Segment* seg);
extern void         blockAlign(int* cut1, int* cut2, double** ocrossscore, int* ncut);
extern void         blockAlign2(Context* ctx, int* cut1, int* cut2, Segment** seg1, Segment** seg2, double** ocrossscore, int* ncut);
extern double       imp_match_out_scD(int i1, int j1);
extern void         imp_match_init_strictD(aln_Opts opts, Context* ctx, int clus1, int clus2, int lgth1, int lgth2, char** seq1, char** seq2, double* eff1, double* eff2, double* eff1kozo, double* eff2kozo, LocalHom*** localhom, char* swaplist, int* memlist1, int* memlist2, int* uselh, int* seedinlh1, int* seedinlh2, int nodeid, int nfiles);
extern double       MSalignmm(aln_Opts opts, Context* ctx, double** n_dynamicmtx, char** seq1, char** seq2, double* eff1, double* eff2, int icyc, int jcyc, int alloclen, char*, char*, char*, char*, int headgp, int tailgp, double*** cpmxchild0, double*** cpmxchild1, double*** cpmxresult, double orieff1, double orieff2);
extern double       MSalignmm_variousdist(Context* ctx, double*** matrices, char** seq1, char** seq2, double* eff1, double* eff2, double** eff1s, double** eff2s, int icyc, int jcyc, int alloclen, char*, char*, char*, char*, int headgp, int tailgp);
extern double       Lalignmm_hmout(Context* ctx, char** seq1, char** seq2, double* eff1, double* eff2, int icyc, int jcyc, char*, char*, char*, double** map);
extern double       MSalign11(Context* ctx, char** seq1, char** seq2, int alloclen);
extern double       A__align_variousdist(aln_Opts opts, Context* ctx, int** which, double*** scoringmatrices, int penalty, int penalty_ex, char** seq1, char** seq2, double* eff1, double* eff2, double** eff1s, double** eff2s, int icyc, int jcyc, int alloclen, int constraint, double* impmatch, char* gs1, char* gs2, char* ge1, char* ge2, int headgp, int tailgp);
extern double       A__align_gapmap(void);
extern double       translate_and_Calign(char** mseq1, char** mseq2, double* effarr1, double* effarr2, int clus1, int clus2, int alloclen);
extern double       Falign_udpari_long(aln_Opts opts, Context* ctx, double*** scoringmatrices, double** scoringmtx, char** seq1, char** seq2, double* eff1, double* eff2, double** eff1s, double** eff2s, int clus1, int clus2, int alloclen, int* fftlog);
extern void         part_imp_match_init(double* imp, int clus1, int clus2, int lgth1, int lgth2, char** seq1, char** seq2, double* eff1, double* eff2, LocalHom*** localhom);
extern double       G__align11psg(Context* ctx, double** codonmtx, double** scoringmtx, char** seq1, char** seq2, int alloclen, int headgp, int tailgp, double* gstart, double* gend);
extern double       G__align11_noalign(Context* ctx, double** scoringmtx, int penal, int penal_ex, char** seq1, char** seq2);
extern double       L__align11_noalign(Context* ctx, double** scoringmtx, char** seq1, char** seq2);
extern double       genL__align11(Context* ctx, double** scoringmtx, char** seq1, char** seq2, int alloclen, int* off1pt, int* off2pt);
extern double       genG__align11(char** seq1, char** seq2, int alloclen);
extern double       VAalign11(char** seq1, char** seq2, int alloclen, int* off1pt, int* off2pt, LocalHom* lhmpt);
extern int          fft(int n, Fukusosuu* x, int dum);
extern void         JTTmtx(double** rsr, double* freq, unsigned char locamino[0x80], char locgrp[0x80], int isTM);
extern void         BLOSUMmtx(Context* ctx, int n, double** matrix, double* freq, unsigned char* amino, char* amino_grp, int* rescale);
extern int          extendedmtx(Context* ctx, double** matrix, double* freq, unsigned char* amino, char* amino_grp);
extern void         putlocalhom2(Context* ctx, char* al1, char* al2, LocalHom* localhompt, int off1, int off2, char korh);
extern void         putlocalhom_str(char* al1, char* al2, double* equiv, double scale, LocalHom* localhompt, int off1, int off2, char korh);
extern void         putlocalhom_ext(Context* ctx, char* al1, char* al2, LocalHom* localhompt, int off1, int off2, char korh);
extern void         putlocalhom(char* al1, char* al2, LocalHom* localhompt, int off1, int off2, int opt, int overlapaa, char korh);
extern char*        cutal(char* al, int al_display_start, int start, int end);
extern void         ErrorExit(char* message);
extern void         seqUpper(int nseq, char** seq);
extern void         seqLower(int nseq, char** seq);
extern int          getaline_fp_eof(char* s, int l, FILE* fp);
extern int          getaline_fp_eof_new(char s[], int l, FILE* fp);
extern int          myfgets(char s[], int l, FILE* fp);
extern double       input_new(FILE* fp, int d);
extern int          allSpace(char* str);
extern void         kake2hiku(char* str);
extern int          countATGC(char* s, int* total);
extern void         writeDataforgaln(FILE* fp, int locnjob, char** name, char** aseq);
extern void         writeData_pointer(Context* ctx, FILE* fp, int locnjob, const char* const* name, char** aseq);
extern void         readhat2_doublehalf_pointer(FILE* fp, int nseq, double** mtx);
extern void         readhat2_double(FILE* fp, int nseq, double** mtx);
extern void         readhat2_int(FILE* fp, int nseq, int** mtx);
extern void         readhat2_pointer(FILE* fp, int nseq, double** mtx);
extern void         readhat2(FILE* fp, int nseq, double** mtx);
extern void         WriteFloatHat2_pointer_halfmtx(Context* ctx, FILE* hat2p, int locnjob, const char* const* name, double** mtx);
extern void         WriteHat2_pointer(FILE* hat2p, int locnjob, char** name, double** mtx);
extern void         WriteHat2_part_pointer(FILE* hat2p, int locnjob, int nadd, const char* const* name, double** mtx);
extern int          ReadBlastm7_scoreonly(FILE* fp, double* dis, int nin);
extern int          ReadBlastm7_avscore(FILE* fp, double* dis, int nin);
extern void         initSignalSM(Context* ctx);
extern void         initFiles(Context* ctx);
extern void         readlocalhomtable(FILE* fp, int njob, LocalHom** localhomtable, char* kozoarivec);
extern void         readlocalhomtable_half(FILE* fp, int njob, LocalHom** localhomtable, char* kozoarivec);
extern void         readlocalhomtable_target(FILE* fp, int nt, int njob, LocalHom** localhomtable, char* kozoarivec, int* targetmap);
extern void         readlocalhomtable2(FILE* fp, LocalHom** localhomtable, char* kozoarivec);
extern void         readlocalhomtable2_half(FILE* fp, int njob, LocalHom** localhomtable, char* kozoarivec);
extern void         readlocalhomtable2_target(FILE* fp, LocalHom** localhomtable, char* kozoarivec, int* targetmap);
extern void         readlocalhomtable_part(FILE* fp, int njob, int nadd, LocalHom** localhomtable, char* kozoarivec);
extern void         readlocalhomtable_two(FILE* fp, int njob, int nadd, LocalHom** localhomtable, LocalHom** localhomtablex);
extern void         readlocalhomtable_one(FILE* fp, int njob, int nadd, LocalHom** localhomtable);
extern void         outlocalhom(LocalHom** localhom, int nseq);
extern void         outlocalhom_part(LocalHom** localhom, int norg, int nadd);
extern void         outlocalhompt(LocalHom*** localhom, int n1, int n2);
extern void         FreeLocalHomTable_half(LocalHom** localhomtable, int n);
extern void         FreeLocalHomTable(LocalHom** localhomtable, int n);
extern void         FreeLocalHomTable_part(LocalHom** localhomtable, int n, int m);
extern void         FreeLocalHomTable_two(LocalHom** localhomtable, int n, int m);
extern void         FreeLocalHomTable_one(LocalHom** localhomtable, int n, int m);
extern void         freelocalhom1(LocalHom* lh);
extern void         initlocalhom1(LocalHom* lh);
extern void         clustalout_pointer(FILE* fp, int nseq, int maxlen, char** seq, char** name, char* mark, char* comment, int* order, int namelen);
extern void         phylipout_pointer(FILE* fp, int nseq, int maxlen, char** seq, char** name, int* order, int namelen);
extern void         writeData_reorder(FILE* fp, int locnjob, char name[][B], char** aseq, int* order);

extern int                load1SeqWithoutName_new(Context* ctx, FILE* fpp, char* cbuf);
extern char*              load1SeqWithoutName_realloc_casepreserve(FILE* fpp);
extern void               searchKUorWA(FILE* fp);
extern void               gapireru(Context* ctx, char* res, char* ori, char* gt);
extern int                seqlen(Context* ctx, char* seq);
extern void               st_FinalGapCount(double* fgcp, int clus, char** seq, double* eff, int len);
extern void               st_FinalGapAdd(double* fgcp, int clus, char** seq, double* eff, int len);
extern void               st_OpeningGapCount(double* ogcp, int clus, char** seq, double* eff, int len);
extern void               st_OpeningGapAdd(double* ogcp, int clus, char** seq, double* eff, int len);
extern void               st_FinalGapCount_zure(double* fgcp, int clus, char** seq, double* eff, int len);
extern void               getdiaminofreq_x(double* freq, int clus, char** seq, double* eff, int len);
extern void               new_FinalGapCount_zure(double* fgcp, int clus, char** seq, double* eff, int len, char* s, char* e);
extern void               new_FinalGapCount(double* fgcp, int clus, char** seq, double* eff, int len, char* g);
extern void               new_OpeningGapCount(double* ogcp, int clus, char** seq, double* eff, int len, char* g);
extern void               new_OpeningGapCount_zure(double* ogcp, int clus, char** seq, double* eff, int len, char* s, char* e);
extern void               getGapPattern(double* fgcp, int clus, char** seq, double* eff, int len);
extern void               getgapfreq(double* freq, int clus, char** seq, double* eff, int len);
extern void               getgapfreq_zure(double* freq, int clus, char** seq, double* eff, int len);
extern void               getgapfreq_zure_part(double* freq, int clus, char** seq, double* eff, int len, char* s);
extern void               getdiaminofreq_part(double* freq, int clus, char** seq, double* eff, int len, char* s, char* e);
extern void               getdigapfreq_part(double* freq, int clus, char** seq, double* eff, int len, char* s, char* e);
extern void               getdiaminofreq_st(double* freq, int clus, char** seq, double* eff, int len);
extern void               getdigapfreq_st(double* freq, int clus, char** seq, double* eff, int len);
extern void               st_getGapPattern(Gappat** gpat, int clus, char** seq, double* eff, int len);
extern void               getkyokaigap(char* g, char** s, int pos, int n);
extern double*            loadaamtx(Context* ctx, int* rescalept);
extern double             naivepairscore11(Context* ctx, const char* seq1, const char* seq2, int penal);
extern double             naivepairscore11_dynmtx(Context* ctx, double**, char* seq1, char* seq2, int penal);
extern double             naivepairscorefast(Context* ctx, const char* seq1, const char* seq2, int* skip1, int* skip2, int penal);
extern void               foldrna(Context* ctx, int nseq1, int nseq2, char** seq1, char** seq2, double* eff1, double* eff2, RNApair*** gr1, RNApair*** gr2, double** impmtx);
extern void               foldrna_gappick(int nseq1, int nseq2, char** seq1, char** seq2, double* eff1, double* eff2, RNApair*** gr1, RNApair*** gr2, double** impmtx, int* gapmap1, int* gapmap2, RNApair* pair);
extern void               imp_rna(Context* ctx, int nseq1, int nseq2, char** seq1, char** seq2, double* eff1, double* eff2, RNApair*** gr1, RNApair*** gr2);
extern void               imp_rnaD(Context* ctx, int nseq1, int nseq2, char** seq1, char** seq2, double* eff1, double* eff2, RNApair*** gr1, RNApair*** gr2);
extern void               foldalignedrna(int clus1, int clus2, char** mseq1, char** mseq2, double* effarr1, double* effarr2, RNApair* rnapairboth);
void                      readmccaskill(FILE* fp, RNApair** pairprob, int length);
void                      makegrouprna(RNApair*** group, RNApair*** all, int* memlist);
void                      makegrouprnait(RNApair*** group, RNApair*** all, char* pair, int s);
extern void               fixed_musclesupg_double_realloc_nobk_halfmtx(aln_Opts opts, Context* ctx, int nseq, double** eff, int*** topol, double** len, Treedep*, int progressout, int efffree);
extern void               fixed_musclesupg_double_realloc_nobk_halfmtx_memsave(aln_Opts opts, Context* ctx, int nseq, double** eff, int*** topol, double** len, Treedep*, int progressout, int efffree);
extern void               loadtree(Context* ctx, int nseq, int*** topol, double** len, const char* const* name, Treedep*, int treeout);
extern int                check_guidetreefile(int* seed, int* npick, double* limitram);
extern void               fixed_musclesupg_double_realloc_nobk_halfmtx_treeout_memsave(aln_Opts opts, Context* ctx, int nseq, double** eff, int*** topol, double** len, const char* const* name, Treedep*, int efffree, int treeout);
extern void               fixed_supg_double_realloc_nobk_halfmtx_treeout_constrained(aln_Opts opts, Context* ctx, int nseq, double** eff, int*** topol, double** len, const char* const* name, Treedep*, int ncons, int** constraints, int efffree);
extern void               imp_match_init_strict(aln_Opts opts, Context* ctx, int clus1, int clus2, int lgth1, int lgth2, char** seq1, char** seq2, double* eff1, double* eff2, double* eff1kozo, double* eff2kozo, LocalHom*** localhom, char* swaplist, int* memlist1, int* memlist2, int* uselh, int* seedinlh1, int* seedinlh2, int nodeid, int nfiles);
extern void               miyataout_reorder_pointer(FILE* fp, int locnjob, int nlenmax, char** name, int* nlen, char** aseq, int* order);
extern void               cpmx_ribosum(Context* ctx, char** seq, char** seqr, char* dir, double** cpmx, double* eff, int lgth, int clus);
extern void               rnaalifoldcall(Context* ctx, char** seq, int nseq, RNApair** pairprob);
extern void               readpairfoldalign(FILE* fp, char* seq1, char* seq2, char* aln1, char* aln2, int q1, int q2, int* of1, int* of2, int sumlen);
extern void               write1seq(FILE* fp, char* aseq);
extern void               assignstrweight(int nseq, double* strweight, Node* stopol, int*** topol, int step, int LorR, char* kozoari, double* seqweight);
extern double             plainscore(Context* ctx, int nseq, char** seq);
extern void               eq2dash(char* s);
extern void               eq2dashmatometehayaku(char** s, int n);
extern void               findnewgaps(Context* ctx, int rep, char** seq, int* gaplen);
extern void               findcommongaps(int, char**, int*);
extern void               adjustgapmap(int, int*, char*);
extern void               insertnewgaps_bothorders(aln_Opts opts, Context* ctx, int njob, int* alreadyaligned, char** seq, int* ex1, int* ex2, int* gaplen, int* gapmap, int gapmaplen, int alloclen, char alg, char gapchar);
extern void               insertnewgaps(aln_Opts opts, Context* ctx, int njob, int* alreadyaligned, char** seq, int* ex1, int* ex2, int* gaplen, int* gapmap, int alloclen, char alg, char gapchar);
extern void               restorecommongaps(int n, int n0, char** seq, int* top0, int* top1, int* gaplen, int alloclen, char gapchar);
extern void               restorecommongapssmoothly(int n, int n0, char** seq, int* top0, int* top1, int* gaplen, int alloclen, char gapchar);
extern int                samemember(int* mem, int* cand);
extern int                samemembern(int* mem, int* cand, int candn);
extern int                includemember(int* mem, int* cand);
extern int                overlapmember(int* mem1, int* mem2);
extern void               sreverse(char* r, char* s);
extern void               intcpy(int* s1, int* s2);
extern void               intncpy(int* s1, int* s2, int n);
extern void               fltncpy(double* s1, double* s2, int n);
extern void               intcat(int* s1, int* s2);
extern void               readsubalignmentstable(int n, int** table, int* preservegaps, int* nsubpt, int* maxmempt);
extern int                myatoi(char*);
extern unsigned long long myatoll(char*);
extern double             myatof(char*);
extern void               gapcount(double* freq, char** seq, int nseq, double* eff, int lgth);
extern void               gapcountf(double* freq, char** seq, int nseq, double* eff, int lgth);
extern void               gapcountadd(double* freq, char** seq, int nseq, double* eff, int lgth);
extern void               outgapcount(double* freq, int nseq, char* gappat, double* eff);
extern void               makedynamicmtx(aln_Opts opts, Context* ctx, double** out, double** in, double offset);
extern double             dist2offset(aln_Opts opts, double dist);
extern void               reporterr(const char* str, ...);
extern void               freeconstants(Context* ctx);
extern void               closeFiles(Context* ctx);
extern void               FreeCommonIP(Context* ctx);
extern void               makeskiptable(int n, int** skip, char** seq);
extern int                generatesubalignmentstable(int nseq, int*** tablept, int* nsubpt, int* maxmempt, int*** topol, double** len, double threshold);
extern double             sumofpairsscore(Context* ctx, int nseq, char** seq);

extern int    isallgap(char*);
extern int    deletenewinsertions_whole(int on, int an, char** oseq, char** aseq, GapPos** deletelist);
extern int    deletenewinsertions_whole_eq(int on, int an, char** oseq, char** aseq, GapPos** deletelist);
extern int    deletenewinsertions_difflist(int on, int an, char** oseq, char** aseq, GapPos** difflist);
extern int    recordoriginalgaps(char* originallygapped, int n, char** s);
extern void   restoreoriginalgaps(int n, char** seq, char* originalgaps);
extern void   reconstructdeletemap(Context* ctx, int nadd, char** addbk, GapPos** deletelist, char** realn, FILE* fp, const char* const* name);
extern void   reconstructdeletemap_compact(Context* ctx, int nadd, char** addbk, GapPos** deletelist, char** realn, FILE* fp, const char* const* name);
extern double D__align(aln_Opts opts, Context* ctx, double** n_dynamicmtx, char** seq1, char** seq2, double* eff1, double* eff2, int icyc, int jcyc, int alloclen, int constraint, double* impmatch, int headgp, int tailgp);
extern double D__align_variousdist(aln_Opts opts, Context* ctx, int** whichmtx, double*** matrices, char** seq1, char** seq2, double* eff1, double* eff2, double** eff1s, double** eff2s, int icyc, int jcyc, int alloclen, int constraint, double* impmatch, int headgp, int tailgp);
extern void   stringshuffle(int* ary, int size);
extern void   topolorder(int* order, int* posinorder, int*** topol, Treedep* dep, int pos, int child);
extern int*   topolorderz(int* order, int*** topol, Treedep* dep, int pos, int nchild);
extern int*   topolordery(int* order, int*** topol, Treedep* dep, int pos, int nchild);
extern void   compacttree_memsaveselectable(aln_Opts opts, Context* ctx, int nseq, double** partmtx, int* nearest, double* mindist, int** pointt, int* selfscore, char** seq, int** skiptable, int*** topol, double** len, const char* const* name, int* nlen, Treedep* dep, int treeout, int howcompact, int memsave);
extern double distcompact(Context* ctx, int len1, int len2, int* table1, int* point2, int ss1, int ss2);
extern double distcompact_msa(Context* ctx, const char* seq1, const char* seq2, int* skiptable1, int* skiptable2, int ss1, int ss2);
extern void   fillimp(aln_Opts opts, double** impmtx, int clus1, int clus2, int lgth1, int lgth2, char** seq1, char** seq2, double* eff1, double* eff2, double* eff1_kozo, double* eff2_kozo, LocalHom*** localhom, char* swaplist, int* orinum1, int* orinum2);
extern void   fillimp_file(aln_Opts opts, Context* ctx, double** impmtx, int clus1, int clus2, int lgth1, int lgth2, char** seq1, char** seq2, double* eff1, double* eff2, double* eff1_kozo, double* eff2_kozo, LocalHom*** localhom, int* orinum1, int* orinum2, int* uselh, int* seedinlh1, int* seedinlh2, int nodeid, int nfiles);
extern int    pairlocalalign(aln_Opts opts, Context* ctx, int ngui, const char* const* namegui, char** seqgui, double** distancemtx, LocalHom** localhomtable, int argc, char** argv, double** expdist);
extern char   creverse(char f);
extern void   setstacksize(rlim_t);
extern void   use_getrusage(void);
extern void   treeout_bin(FILE* treefp, int n, int*** topol, double** len, Treedep* dep, int* nfilesfornode);
extern void   treein_bin(FILE* treefp, int n, int*** topol, double** len, Treedep* dep, int* nfilesfornode);
extern void   uselhout(FILE*, int n, int*);
extern int    uselhin(FILE*, int n, int*);
extern void   sortbylength(int* uselh, Lennum* in, int size, unsigned long long numpairs);
extern void   limitlh(int* uselh, Lennum* in, int size, int limit);

extern double   distdp_noalign(Context* ctx, char* s1, char* s2, double selfscore1, double selfscore2, int alloclen);  // tbfast.c kara yobareru
extern void     getweightfromname(int n, double* w, char** name);
extern void     readexternalanchors(ExtAnch** extanch, int nseq, int* nogaplen);
extern double** loadcodonscore(FILE*, double** mtx);
extern int      codon2id(char*);
