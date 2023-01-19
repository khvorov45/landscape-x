#include "mltaln.h"

int
intlen(int* num) {
    int* numbk = num;
    while (*num++ != -1)
        ;
    return (num - numbk - 1);
}

void
intcpy(int* s1, int* s2) {
    while (*s2 != -1) {
        *s1++ = *s2++;
    }
    *s1 = -1;
}

typedef struct _TopDep {
    Treedep* dep;
    int***   topol;
} TopDep;

static TopDep* tdpglobal = NULL;
static int*
topolorder_lessargs(int* order, int pos) {
    if ((tdpglobal->dep)[pos].child0 == -1) {
        *order++ = (tdpglobal->topol)[pos][0][0];
        *order = -1;
    } else {
        order = topolorder_lessargs(order, (tdpglobal->dep)[pos].child0);
    }

    if ((tdpglobal->dep)[pos].child1 == -1) {
        *order++ = (tdpglobal->topol)[pos][1][0];
        *order = -1;
    } else {
        order = topolorder_lessargs(order, (tdpglobal->dep)[pos].child1);
    }

    return (order);
}

int*
topolorderz(int* order, int*** topol, Treedep* dep, int pos, int nchild) {
    tdpglobal = (TopDep*)calloc(sizeof(TopDep), 1);
    tdpglobal->topol = topol;
    tdpglobal->dep = dep;

    int child;

    if (nchild == 0 || nchild == 2) {
        if ((child = (dep)[pos].child0) == -1) {
            *order++ = (topol)[pos][0][0];
            *order = -1;
        } else {
            order = topolorder_lessargs(order, child);
        }
    }
    if (nchild == 1 || nchild == 2) {
        if ((child = (dep)[pos].child1) == -1) {
            *order++ = (topol)[pos][1][0];
            *order = -1;
        } else {
            order = topolorder_lessargs(order, child);
        }
    }

#if 1
    free(tdpglobal);
    tdpglobal = NULL;
#endif

    return (order);
}

static void
setnearest(Bchain* acpt, double** iscore, double* mindisfrompt, int* nearestpt, int pos) {
    double mindisfrom = 999.9;
    int    nearest = -1;

    for (Bchain* acptj = (acpt + pos)->next; acptj; acptj = acptj->next) {
        int    j = acptj->pos;
        double tmpdouble = iscore[pos][j - pos];
        if (tmpdouble < mindisfrom) {
            mindisfrom = tmpdouble;
            nearest = j;
        }
    }

    for (Bchain* acptj = acpt; (acptj && acptj->pos != pos); acptj = acptj->next) {
        int    j = acptj->pos;
        double tmpdouble = iscore[j][pos - j];
        if (tmpdouble < mindisfrom) {
            mindisfrom = tmpdouble;
            nearest = j;
        }
    }

    *mindisfrompt = mindisfrom;
    *nearestpt = nearest;
}

void
fixed_musclesupg_double_realloc_nobk_halfmtx_memsave(aln_Context* ctx, double** iscore, int*** topol, double** len, Treedep* dep) {
    int     i, j, miniim, maxiim, minijm, maxijm;
    int*    intpt;
    double  tmpdouble;
    double  eff1, eff0;
    int     im = -1, jm = -1;
    Bchain *acjmnext, *acjmprev;
    int     prevnode;
    Bchain* acpti;
    int *   pt1, *pt2, *pt11;
    int     nmemim, nmemjm;
    double  minscore;

    double sueff1 = 1 - (double)ctx->sueff_global;
    double sueff05 = (double)ctx->sueff_global * 0.5;

    int*    hist = AllocateIntVec(ctx->njob);
    double* tmptmplen = AllocateFloatVec(ctx->njob);
    Bchain* ac = (Bchain*)malloc(ctx->njob * sizeof(Bchain));
    int*    nmemar = AllocateIntVec(ctx->njob);
    double* mindisfrom = AllocateFloatVec(ctx->njob);
    int*    nearest = AllocateIntVec(ctx->njob);

    for (int32_t i = 0; i < ctx->njob; i++) {
        ac[i].next = ac + i + 1;
        ac[i].prev = ac + i - 1;
        ac[i].pos = i;
    }
    ac[ctx->njob - 1].next = 0;

    for (int32_t i = 0; i < ctx->njob; i++) {
        setnearest(ac, iscore, mindisfrom + i, nearest + i, i);
    }

    for (int32_t i = 0; i < ctx->njob; i++)
        tmptmplen[i] = 0.0;
    for (int32_t i = 0; i < ctx->njob; i++) {
        hist[i] = -1;
        nmemar[i] = 1;
    }

    for (int32_t k = 0; k < ctx->njob - 1; k++) {
        minscore = 999.9;
        for (acpti = ac; acpti->next != NULL; acpti = acpti->next) {
            i = acpti->pos;
            if (mindisfrom[i] < minscore) {
                im = i;
                minscore = mindisfrom[i];
            }
        }
        jm = nearest[im];
        if (jm < im) {
            j = jm;
            jm = im;
            im = j;
        }

        prevnode = hist[im];
        if (dep)
            dep[k].child0 = prevnode;
        nmemim = nmemar[im];
        intpt = topol[k][0] = (int*)realloc(topol[k][0], (2) * sizeof(int));
        if (prevnode == -1) {
            *intpt++ = im;
            *intpt = -1;
        } else {
            pt1 = topol[prevnode][0];
            pt2 = topol[prevnode][1];
            if (*pt1 > *pt2) {
                pt11 = pt2;
            } else {
                pt11 = pt1;
            }
            *intpt++ = *pt11;
            *intpt = -1;
        }

        prevnode = hist[jm];
        if (dep)
            dep[k].child1 = prevnode;
        nmemjm = nmemar[jm];
        intpt = topol[k][1] = (int*)realloc(topol[k][1], (2) * sizeof(int));
        if (!intpt) {
            exit(1);
        }
        if (prevnode == -1) {
            *intpt++ = jm;
            *intpt = -1;
        } else {
            pt1 = topol[prevnode][0];
            pt2 = topol[prevnode][1];
            if (*pt1 > *pt2) {
                pt11 = pt2;
            } else {
                pt11 = pt1;
            }

            *intpt++ = *pt11;
            *intpt = -1;
        }

        minscore *= 0.5;

        len[k][0] = (minscore - tmptmplen[im]);
        len[k][1] = (minscore - tmptmplen[jm]);

        if (dep)
            dep[k].distfromtip = minscore;

        tmptmplen[im] = minscore;

        hist[im] = k;
        nmemar[im] = nmemim + nmemjm;

        mindisfrom[im] = 999.9;
        for (acpti = ac; acpti != NULL; acpti = acpti->next) {
            i = acpti->pos;
            if (i != im && i != jm) {
                if (i < im) {
                    miniim = i;
                    maxiim = im;
                    minijm = i;
                    maxijm = jm;
                } else if (i < jm) {
                    miniim = im;
                    maxiim = i;
                    minijm = i;
                    maxijm = jm;
                } else {
                    miniim = im;
                    maxiim = i;
                    minijm = jm;
                    maxijm = i;
                }
                eff0 = iscore[miniim][maxiim - miniim];
                eff1 = iscore[minijm][maxijm - minijm];
                tmpdouble = iscore[miniim][maxiim - miniim] = MIN(eff0, eff1) * sueff1 + (eff0 + eff1) * sueff05;
                if (tmpdouble < mindisfrom[i]) {
                    mindisfrom[i] = tmpdouble;
                    nearest[i] = im;
                }
                if (tmpdouble < mindisfrom[im]) {
                    mindisfrom[im] = tmpdouble;
                    nearest[im] = i;
                }
                if (nearest[i] == jm) {
                    nearest[i] = im;
                }
            }
        }

        acjmprev = ac[jm].prev;
        acjmnext = ac[jm].next;
        acjmprev->next = acjmnext;
        if (acjmnext != NULL) {
            acjmnext->prev = acjmprev;
        }

        iscore[jm] = NULL;

        for (acpti = ac; acpti != NULL; acpti = acpti->next) {
            i = acpti->pos;
            if (nearest[i] == im) {
                if (i < im) {
                    miniim = i;
                    maxiim = im;
                } else {
                    miniim = im;
                    maxiim = i;
                }
                if (iscore[miniim][maxiim - miniim] > mindisfrom[i]) {
                    setnearest(ac, iscore, mindisfrom + i, nearest + i, i);
                }
            }
        }
    }
}

void
counteff_simple_double_nostatic_memsave(int nseq, int*** topol, double** len, Treedep* dep, double* node) {
    int     j, s1, s2;
    double  total;
    double* rootnode;
    double* eff;
    int**   localmem;
    int**   memhist;

    rootnode = AllocateDoubleVec(nseq);
    eff = AllocateDoubleVec(nseq);
    localmem = AllocateIntMtx(2, 0);
    memhist = AllocateIntMtx(nseq - 1, 0);
    for (int32_t i = 0; i < nseq - 1; i++) {
        memhist[i] = NULL;
    }

    for (int32_t i = 0; i < nseq; i++) {
        if (len[i][0] < 0.0) {
            len[i][0] = 0.0;
        }
        if (len[i][1] < 0.0) {
            len[i][1] = 0.0;
        }
    }

    for (int32_t i = 0; i < nseq; i++) {
        rootnode[i] = 0.0;
        eff[i] = 1.0;
    }

    for (int32_t i = 0; i < nseq - 1; i++) {
        if (dep[i].child0 == -1) {
            localmem[0] = calloc(sizeof(int), 2);
            localmem[0][0] = topol[i][0][0];
            localmem[0][1] = -1;
            s1 = 1;
        } else {
            localmem[0] = memhist[dep[i].child0];
            s1 = intlen(localmem[0]);
        }
        if (dep[i].child1 == -1) {
            localmem[1] = calloc(sizeof(int), 2);
            localmem[1][0] = topol[i][1][0];
            localmem[1][1] = -1;
            s2 = 1;
        } else {
            localmem[1] = memhist[dep[i].child1];
            s2 = intlen(localmem[1]);
        }

        memhist[i] = calloc(sizeof(int), s1 + s2 + 1);
        intcpy(memhist[i], localmem[0]);
        intcpy(memhist[i] + s1, localmem[1]);
        memhist[i][s1 + s2] = -1;

        for (j = 0; (s1 = localmem[0][j]) > -1; j++) {
            rootnode[s1] += (double)len[i][0] * eff[s1];
            eff[s1] *= 0.5;
            /*
           	rootnode[s1] *= 0.5;
*/
        }
        for (j = 0; (s2 = localmem[1][j]) > -1; j++) {
            rootnode[s2] += (double)len[i][1] * eff[s2];
            eff[s2] *= 0.5;
            /*
           	rootnode[s2] *= 0.5;
*/
        }
        free(localmem[0]);
        free(localmem[1]);
    }
    free(localmem);
    free(memhist[nseq - 2]);
    free(memhist);

    for (int32_t i = 0; i < nseq; i++) {
        rootnode[i] += GETA3;
    }
    total = 0.0;
    for (int32_t i = 0; i < nseq; i++) {
        total += rootnode[i];
    }

    for (int32_t i = 0; i < nseq; i++) {
        node[i] = rootnode[i] / total;
    }

    free(rootnode);
    free(eff);
}

#define SEGMENTSIZE 150

static int
seqlen(char* seq) {
    int val = 0;
    while (*seq) {
        if (*seq++ != '-') {
            val++;
        }
    }
    return val;
}

void
calcimportance_half(int nseq, double* eff, char** seq, aln_LocalHom** localhom, int alloclen) {
    int           j, pos, len;
    double*       importance;
    double        tmpdouble;
    double *      ieff, totaleff;
    int*          nogaplen;
    aln_LocalHom* tmpptr;

    importance = AllocateDoubleVec(alloclen);
    nogaplen = AllocateIntVec(nseq);
    ieff = AllocateDoubleVec(nseq);

    totaleff = 0.0;
    for (int32_t i = 0; i < nseq; i++) {
        nogaplen[i] = seqlen(seq[i]);
        if (nogaplen[i] == 0)
            ieff[i] = 0.0;
        else
            ieff[i] = eff[i];
        totaleff += ieff[i];
    }
    for (int32_t i = 0; i < nseq; i++)
        ieff[i] /= totaleff;

    for (int32_t i = 0; i < nseq; i++) {
        for (pos = 0; pos < alloclen; pos++) {
            importance[pos] = 0.0;
        }
        for (j = 0; j < nseq; j++) {
            if (i == j)
                continue;

            else if (i < j) {
                for (tmpptr = localhom[i] + j - i; tmpptr; tmpptr = tmpptr->next) {
                    if (tmpptr->opt == -1)
                        continue;
                    for (pos = tmpptr->start1; pos <= tmpptr->end1; pos++) {
                        importance[pos] += ieff[j];
                    }
                }
            } else {
                for (tmpptr = localhom[j] + i - j; tmpptr; tmpptr = tmpptr->next) {
                    if (tmpptr->opt == -1)
                        continue;
                    for (pos = tmpptr->start2; pos <= tmpptr->end2; pos++) {
                        importance[pos] += ieff[j];
                    }
                }
            }
        }

        for (j = 0; j < nseq; j++) {
            if (i == j)
                continue;

            else if (i < j) {
                if (localhom[i][j - i].opt == -1.0)
                    continue;

                for (tmpptr = localhom[i] + j - i; tmpptr; tmpptr = tmpptr->next) {
                    if (tmpptr->opt == -1.0)
                        continue;
                    tmpdouble = 0.0;
                    len = 0;
                    for (pos = tmpptr->start1; pos <= tmpptr->end1; pos++) {
                        tmpdouble += importance[pos];
                        len++;
                    }

                    tmpdouble /= (double)len;
                    tmpptr->importance = tmpdouble * tmpptr->opt;
                }
            } else {
                if (localhom[j][i - j].opt == -1.0)
                    continue;

                for (tmpptr = localhom[j] + i - j; tmpptr; tmpptr = tmpptr->next) {
                    if (tmpptr->opt == -1.0)
                        continue;
                    tmpdouble = 0.0;
                    len = 0;
                    for (pos = tmpptr->start2; pos <= tmpptr->end2; pos++) {
                        tmpdouble += importance[pos];
                        len++;
                    }
                    tmpdouble /= (double)len;
                    tmpptr->rimportance = tmpdouble * tmpptr->opt;
                }
            }
        }
    }

    for (int32_t i = 0; i < nseq - 1; i++) {
        for (j = i + 1; j < nseq; j++) {
            double        imp;
            aln_LocalHom* tmpptr1;
            tmpptr1 = localhom[i] + j - i;
            for (; tmpptr1; tmpptr1 = tmpptr1->next) {
                if (tmpptr1->opt == -1.0) {
                    continue;
                }
                imp = 0.5 * (tmpptr1->importance + tmpptr1->rimportance);
                tmpptr1->importance = tmpptr1->rimportance = imp;
            }
        }
    }

    free(importance);
    free(nogaplen);
    free(ieff);
}

void
gapireru(char* res, char* ori, char* gt) {
    char g = 0;
    char gapchar = '-';
    while ((g = *gt++)) {
        if (g == '-') {
            *res++ = gapchar;
        } else {
            *res++ = *ori++;
        }
    }
    *res = 0;
}

void
getkyokaigap(char* g, char** s, int pos, int n) {
    while (n--)
        *g++ = (*s++)[pos];
}

void
new_OpeningGapCount(double* ogcp, int clus, char** seq, double* eff, int len, char* sgappat) {
    int     i, j, gc, gb;
    double  feff;
    double* fpt;
    char*   spt;

    fpt = ogcp;
    i = len;
    while (i--)
        *fpt++ = 0.0;
    for (j = 0; j < clus; j++) {
        feff = (double)eff[j];
        spt = seq[j];
        fpt = ogcp;
        gc = (sgappat[j] == '-');
        i = len;
        while (i--) {
            gb = gc;
            gc = (*spt++ == '-');
            {
                if (!gb * gc)
                    *fpt += feff;
                fpt++;
            }
        }
    }
}

void
new_FinalGapCount(double* fgcp, int clus, char** seq, double* eff, int len, char* egappat) {
    int     i, j, gc, gb;
    double  feff;
    double* fpt;
    char*   spt;

    fpt = fgcp;
    i = len;
    while (i--)
        *fpt++ = 0.0;
    for (j = 0; j < clus; j++) {
        feff = (double)eff[j];
        fpt = fgcp;
        spt = seq[j];
        gc = (*spt == '-');
        i = len;
        while (i--) {
            gb = gc;
            gc = (*++spt == '-');
            {
                if (gb * !gc)
                    *fpt += feff;
                fpt++;
            }
        }
        {
            gb = gc;
            gc = (egappat[j] == '-');
            {
                if (gb * !gc)
                    *fpt += feff;
            }
        }
    }
}

void
st_OpeningGapAdd(double* ogcp, int clus, char** seq, double* eff, int len) {
    int     i, j, gc, gb;
    double* fpt;
    char*   spt;
    int     newmem = clus - 1;
    double  neweff = eff[newmem];
    double  orieff = 1.0 - neweff;
    double  feff;

    j = clus - 1;
    {
        feff = (double)eff[j];
        spt = seq[j];
        fpt = ogcp;
        i = len;
        gc = 0;
        while (i--) {
            gb = gc;
            gc = (*spt++ == '-');
            *fpt *= orieff;
            if (!gb * gc)
                *fpt += feff;
            fpt++;
        }
    }
    ogcp[len] = 0.0;
}

void
st_OpeningGapCount(double* ogcp, int clus, char** seq, double* eff, int len) {
    int     i, j, gc, gb;
    double  feff;
    double* fpt;
    char*   spt;

    fpt = ogcp;
    i = len;
    while (i--)
        *fpt++ = 0.0;
    for (j = 0; j < clus; j++) {
        feff = (double)eff[j];
        spt = seq[j];
        fpt = ogcp;
        gc = 0;
        //		gc = 1;
        i = len;
        while (i--) {
            gb = gc;
            gc = (*spt++ == '-');
            {
                if (!gb * gc)
                    *fpt += feff;
                fpt++;
            }
        }
    }
    ogcp[len] = 0.0;
}

void
st_FinalGapAdd(double* fgcp, int clus, char** seq, double* eff, int len) {
    int     i, j, gc, gb;
    double* fpt;
    char*   spt;
    int     newmem = clus - 1;
    double  neweff = eff[newmem];
    double  orieff = 1.0 - neweff;
    double  feff;

    j = clus - 1;
    //	for( j=0; j<clus; j++ )
    {
        feff = (double)eff[j];
        fpt = fgcp;
        spt = seq[j];
        gc = (*spt == '-');
        i = len;
        //		for( i=1; i<len; i++ )
        while (i--) {
            *fpt *= orieff;
            gb = gc;
            gc = (*++spt == '-');
            {
                if (gb * !gc)
                    *fpt += feff;
                fpt++;
            }
        }
        {
            *fpt *= orieff;
            gb = gc;
            gc = 0;
            //			gc = 1;
            {
                if (gb * !gc)
                    *fpt += feff;
            }
        }
    }
}

void
st_FinalGapCount(double* fgcp, int clus, char** seq, double* eff, int len) {
    int     i, j, gc, gb;
    double  feff;
    double* fpt;
    char*   spt;

    fpt = fgcp;
    i = len;
    while (i--)
        *fpt++ = 0.0;
    for (j = 0; j < clus; j++) {
        feff = (double)eff[j];
        fpt = fgcp;
        spt = seq[j];
        gc = (*spt == '-');
        i = len;
        //		for( i=1; i<len; i++ )
        while (i--) {
            gb = gc;
            gc = (*++spt == '-');
            {
                if (gb * !gc)
                    *fpt += feff;
                fpt++;
            }
        }
        {
            gb = gc;
            gc = 0;
            //			gc = 1;
            {
                if (gb * !gc)
                    *fpt += feff;
            }
        }
    }
}

void
getGapPattern(double* fgcp, int clus, char** seq, double* eff, int len) {
    int     i, j, gc, gb;
    double  feff;
    double* fpt;
    char*   spt;

    fpt = fgcp;
    i = len + 1;
    while (i--)
        *fpt++ = 0.0;
    for (j = 0; j < clus; j++) {
        feff = (double)eff[j];
        fpt = fgcp;
        spt = seq[j];
        gc = (*spt == '-');
        i = len + 1;
        while (i--) {
            gb = gc;
            gc = (*++spt == '-');
            {
                if (gb * !gc)
                    *fpt += feff;
                fpt++;
            }
        }
    }
}

void
gapcountadd(double* freq, char** seq, int nseq, double* eff, int lgth) {
    int    j = nseq - 1;
    double newfr = eff[j];
    double orifr = 1.0 - newfr;
    for (int32_t i = 0; i < lgth; i++) {
        freq[i] = 1.0 - freq[i];
        freq[i] *= orifr;
        if (seq[j][i] == '-')
            freq[i] += newfr;
    }
    return;
}

void
gapcountf(double* freq, char** seq, int nseq, double* eff, int lgth) {
    int    j;
    double fr;
    for (int32_t i = 0; i < lgth; i++) {
        fr = 0.0;
        for (j = 0; j < nseq; j++) {
            if (seq[j][i] == '-')
                fr += eff[j];
        }
        freq[i] = fr;
    }
    return;
}

void
outgapcount(double* freq, int nseq, char* gappat, double* eff) {
    int    j;
    double fr;

    fr = 0.0;
    for (j = 0; j < nseq; j++) {
        if (gappat[j] == '-')
            fr += eff[j];
    }
    *freq = fr;
    return;
}

static double
dist2offset(double dist) {
    double val = dist * 0.5;
    if (val > 0.0)
        val = 0.0;
    return val;
}

void
makedynamicmtx(aln_Context* ctx, double** out, double offset) {
    int j, ii, jj;

    offset = dist2offset(offset * 2.0);  // offset 0..1 -> 0..2

    for (int32_t i = 0; i < ctx->nalphabets; i++)
        for (j = 0; j < ctx->nalphabets; j++) {
            out[i][j] = aln_matrix2get(ctx->n_dis_consweight_multi, i, j);
        }

    if (offset == 0.0)
        return;

    for (int32_t i = 0; i < ctx->nalphabets; i++) {
        ii = (int)ctx->amino[i];
        if (ii == '-')
            continue;  // text no toki arieru
        for (j = 0; j < ctx->nalphabets; j++) {
            jj = (int)ctx->amino[j];
            if (jj == '-')
                continue;  // text no toki arieru
            out[i][j] = aln_matrix2get(ctx->n_dis_consweight_multi, i, j) + offset * 600;
        }
    }

    return;
}

static void
movereg(char* seq1, char* seq2, aln_LocalHom* tmpptr, int* start1pt, int* start2pt, int* end1pt, int* end2pt) {
    char* pt;
    int   tmpint;

    pt = seq1;
    tmpint = -1;
    while (*pt != 0) {
        if (*pt++ != '-')
            tmpint++;
        if (tmpint == tmpptr->start1)
            break;
    }
    *start1pt = (int)(pt - seq1) - 1;

    if (tmpptr->start1 == tmpptr->end1)
        *end1pt = *start1pt;
    else {
        while (*pt != 0) {
            if (*pt++ != '-')
                tmpint++;
            if (tmpint == tmpptr->end1)
                break;
        }
        *end1pt = (int)(pt - seq1) - 1;
    }

    pt = seq2;
    tmpint = -1;
    while (*pt != 0) {
        if (*pt++ != '-')
            tmpint++;
        if (tmpint == tmpptr->start2)
            break;
    }
    *start2pt = (int)(pt - seq2) - 1;
    if (tmpptr->start2 == tmpptr->end2)
        *end2pt = *start2pt;
    else {
        while (*pt != 0) {
            if (*pt++ != '-')
                tmpint++;
            if (tmpint == tmpptr->end2)
                break;
        }
        *end2pt = (int)(pt - seq2) - 1;
    }
}

static void
movereg_swap(char* seq1, char* seq2, aln_LocalHom* tmpptr, int* start1pt, int* start2pt, int* end1pt, int* end2pt) {
    char* pt;
    int   tmpint;

    pt = seq1;
    tmpint = -1;
    while (*pt != 0) {
        if (*pt++ != '-')
            tmpint++;
        if (tmpint == tmpptr->start2)
            break;
    }
    *start1pt = (int)(pt - seq1) - 1;

    if (tmpptr->start2 == tmpptr->end2)
        *end1pt = *start1pt;
    else {
        while (*pt != 0) {
            //			fprintf( stderr, "tmpint = %d, end1 = %d pos = %d\n", tmpint, tmpptr->end1, pt-seq1[i] );
            if (*pt++ != '-')
                tmpint++;
            if (tmpint == tmpptr->end2)
                break;
        }
        *end1pt = (int)(pt - seq1) - 1;
    }

    pt = seq2;
    tmpint = -1;
    while (*pt != 0) {
        if (*pt++ != '-')
            tmpint++;
        if (tmpint == tmpptr->start1)
            break;
    }
    *start2pt = (int)(pt - seq2) - 1;
    if (tmpptr->start1 == tmpptr->end1)
        *end2pt = *start2pt;
    else {
        while (*pt != 0) {
            if (*pt++ != '-')
                tmpint++;
            if (tmpint == tmpptr->end1)
                break;
        }
        *end2pt = (int)(pt - seq2) - 1;
    }
}

void
fillimp(aln_Context* ctx, double** impmtx, int clus1, int clus2, int lgth1, int lgth2, char** seq1, char** seq2, double* eff1, double* eff2, aln_LocalHom*** localhom, char* swaplist, int* orinum1, int* orinum2) {
    int           j, k1, k2, start1, start2, end1, end2;
    double        effij, effijx;
    char *        pt1, *pt2;
    aln_LocalHom* tmpptr;
    void (*movefunc)(char*, char*, aln_LocalHom*, int*, int*, int*, int*);

    for (int32_t i = 0; i < lgth1; i++)
        for (j = 0; j < lgth2; j++)
            impmtx[i][j] = 0.0;
    effijx = 1.0 * ctx->fastathreshold;
    for (int32_t i = 0; i < clus1; i++) {
        if (swaplist && swaplist[i])
            movefunc = movereg_swap;
        else
            movefunc = movereg;
        for (j = 0; j < clus2; j++) {
            if (swaplist == NULL && orinum1 && orinum2)  // muda.
            {
                if (orinum1[i] > orinum2[j])
                    movefunc = movereg_swap;
                else
                    movefunc = movereg;
            }

            effij = eff1[i] * eff2[j] * effijx;
            tmpptr = localhom[i][j];
            while (tmpptr) {
                movefunc(seq1[i], seq2[j], tmpptr, &start1, &start2, &end1, &end2);

                k1 = start1;
                k2 = start2;
                pt1 = seq1[i] + k1;
                pt2 = seq2[j] + k2;
                while (*pt1 && *pt2) {
                    if (*pt1 != '-' && *pt2 != '-') {
                        impmtx[k1][k2] += tmpptr->importance * effij;
                        k1++;
                        k2++;
                        pt1++;
                        pt2++;
                    } else if (*pt1 != '-' && *pt2 == '-') {
                        k2++;
                        pt2++;
                    } else if (*pt1 == '-' && *pt2 != '-') {
                        k1++;
                        pt1++;
                    } else if (*pt1 == '-' && *pt2 == '-') {
                        k1++;
                        pt1++;
                        k2++;
                        pt2++;
                    }
                    if (k1 > end1 || k2 > end2)
                        break;
                }
                tmpptr = tmpptr->next;
            }
        }
    }
}
