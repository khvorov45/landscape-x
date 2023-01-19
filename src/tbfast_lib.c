#include <stdint.h>
#include <stdbool.h>

#define aln_IMPLEMENTATION
#include "mltaln.h"

void
copyWithNoGaps(char* aseq, const char* seq) {
    for (; *seq != 0; seq++) {
        if (*seq != '-')
            *aseq++ = *seq;
    }
    *aseq = 0;
}

aln_AlignResult
tbfast_main(aln_Str* strings, intptr_t stringsCount, void* out, intptr_t outBytes) {
    aln_Arena permArena_ = {.base = out, .size = outBytes / 4};
    aln_Arena tempArena_ = {.base = (uint8_t*)out + permArena_.size, .size = outBytes - permArena_.size};

    aln_Arena* permArena = &permArena_;
    aln_Arena* tempArena = &tempArena_;

    aln_Context* ctx = aln_arenaAllocStruct(tempArena, aln_Context);
    ctx->penalty = -918;
    ctx->penalty_dist = 918;
    ctx->offset = 60;
    ctx->constraint = 2;
    ctx->penalty_ex = -60;
    ctx->minimumweight = 0.00001;
    ctx->fastathreshold = 2.7;
    ctx->sueff_global = 0.1;
    ctx->nalphabets = 26;
    ctx->njob = stringsCount;
    ctx->outgap = 1;
    ctx->consweight_multi = 1.0;

    {
        char* localamino = "ARNDCQEGHILKMFPSTWYVBZX.-J";
        for (int i = 0; i < 128; i++) {
            ctx->amino_n[i] = -1;
        }
        for (int i = 0; i < 26; i++) {
            char amino = localamino[i];
            ctx->amino[i] = amino;
            ctx->amino_n[(int32_t)amino] = i;
        }
    }

    // TODO(sen) Can we use the strings directly as aos?
    // NOTE(sen) Process input
    int maxInputSeqLen = 0;
    for (int32_t strIndex = 0; strIndex < stringsCount; strIndex++) {
        aln_Str str = strings[strIndex];
        maxInputSeqLen = aln_max(maxInputSeqLen, str.len);
    }

    // TODO(sen) Remove the requirement for null-terminated array
    char** seq = aln_arenaAllocArray(tempArena, char*, stringsCount + 1);
    for (int32_t strIndex = 0; strIndex < stringsCount; strIndex++) {
        aln_Str str = strings[strIndex];
        // TODO(sen) Remove the requirement for null-terminated str
        char* thisSeq = aln_arenaAllocArray(tempArena, char, maxInputSeqLen + 1);
        for (int32_t charInd = 0; charInd < str.len; charInd++) {
            thisSeq[charInd] = str.ptr[charInd];
        }
        seq[strIndex] = thisSeq;
    }

    aln_LocalHom** localhomtable = aln_arenaAllocArray(tempArena, aln_LocalHom*, ctx->njob);
    {
        int32_t ilim = ctx->njob;
        for (int32_t i = 0; i < ctx->njob; i++) {
            localhomtable[i] = aln_arenaAllocArray(tempArena, aln_LocalHom, ilim);
            for (int32_t j = 0; j < ilim; j++) {
                localhomtable[i][j].start1 = -1;
                localhomtable[i][j].end1 = -1;
                localhomtable[i][j].start2 = -1;
                localhomtable[i][j].end2 = -1;
                localhomtable[i][j].opt = -1.0;
                localhomtable[i][j].importance = -1.0;
                localhomtable[i][j].next = 0;
            }
            ilim--;
        }
    }

    {
        aln_Matrix2F32 n_distmp = aln_arenaAllocMatrix2F32(tempArena, 20, 20);
        {
            // clang-format off

            // NOTE(sen) BLOSUM62
            double tmpmtx[] = { 
                5.893685,
                -2.120252,  8.210189,
                -2.296072, -0.659672,  8.479856,
                -2.630151, -2.408668,  1.907550,  8.661363,
                -0.612761, -5.083814, -3.989626, -5.189966, 12.873172,
                -1.206025,  1.474162,  0.002529, -0.470069, -4.352838,  7.927704,
                -1.295821, -0.173087, -0.402015,  2.265459, -5.418729,  2.781955,  7.354247,
                0.239392, -3.456163, -0.634136, -1.970281, -3.750621, -2.677743, -3.165266,  8.344902,
                -2.437724, -0.374792,  0.867735, -1.678363, -4.481724,  0.672051, -0.176497, -3.061315, 11.266586,
                -1.982718, -4.485360, -4.825558, -4.681732, -1.841495, -4.154454, -4.791538, -5.587336, -4.847345,  5.997760,
                -2.196882, -3.231860, -5.068375, -5.408471, -1.916207, -3.200863, -4.269723, -5.440437, -4.180099,  2.282412,  5.774148,
                -1.101017,  3.163105, -0.268534, -1.052724, -4.554510,  1.908859,  1.163010, -2.291924, -1.081539, -4.005209, -3.670219,  6.756827,
                -1.402897, -2.050705, -3.226290, -4.587785, -2.129758, -0.631437, -2.997038, -4.014898, -2.326896,  1.690191,  2.987638, -2.032119,  8.088951,
                -3.315080, -4.179521, -4.491005, -5.225795, -3.563219, -4.746598, -4.788639, -4.661029, -1.851231, -0.241317,  0.622170, -4.618016,  0.018880,  9.069126,
                -1.221394, -3.162863, -3.000581, -2.220163, -4.192770, -1.922917, -1.674258, -3.200320, -3.241363, -4.135001, -4.290107, -1.520445, -3.714633, -5.395930, 11.046892,
                1.673639, -1.147170,  0.901353, -0.391548, -1.312485, -0.151708, -0.220375, -0.438748, -1.322366, -3.522266, -3.663923, -0.305170, -2.221304, -3.553533, -1.213470,  5.826527,
                -0.068042, -1.683495, -0.069138, -1.576054, -1.299983, -1.012997, -1.294878, -2.363065, -2.528844, -1.076382, -1.796229, -1.004336, -0.999449, -3.161436, -1.612919,  2.071710,  6.817956,
                -3.790328, -4.019108, -5.543911, -6.321502, -3.456164, -2.919725, -4.253197, -3.737232, -3.513238, -3.870811, -2.447829, -4.434676, -2.137255,  1.376341, -5.481260, -4.127804, -3.643382, 15.756041,
                -2.646022, -2.540799, -3.122641, -4.597428, -3.610671, -2.131601, -3.030688, -4.559647,  2.538948, -1.997058, -1.593097, -2.730047, -1.492308,  4.408690, -4.379667, -2.528713, -2.408996,  3.231335,  9.892544,
                -0.284140, -3.753871, -4.314525, -4.713963, -1.211518, -3.297575, -3.663425, -4.708118, -4.676220,  3.820569,  1.182672, -3.393535,  1.030861, -1.273542, -3.523054, -2.469318, -0.083276, -4.251392, -1.811267,  5.653391
            };

            // clang-format on

            int32_t count = 0;
            for (int32_t i = 0; i < 20; i++) {
                for (int32_t j = 0; j <= i; j++) {
                    aln_matrix2get(n_distmp, i, j) = aln_matrix2get(n_distmp, j, i) = (double)tmpmtx[count++];
                }
            }
        }

        double freq1[20] = {0.077, 0.051, 0.043, 0.052, 0.020, 0.041, 0.062, 0.074, 0.023, 0.052, 0.091, 0.059, 0.024, 0.040, 0.051, 0.069, 0.059, 0.014, 0.032, 0.066};

        {
            double average = 0.0;
            for (int32_t i = 0; i < 20; i++) {
                for (int32_t j = 0; j < 20; j++) {
                    average += aln_matrix2get(n_distmp, i, j) * freq1[i] * freq1[j];
                }
            }
            for (int32_t i = 0; i < 20; i++) {
                for (int32_t j = 0; j < 20; j++) {
                    aln_matrix2get(n_distmp, i, j) -= average;
                }
            }
        }

        {
            double average = 0.0;
            for (int32_t i = 0; i < 20; i++)
                average += aln_matrix2get(n_distmp, i, i) * freq1[i];
            for (int32_t i = 0; i < 20; i++) {
                for (int32_t j = 0; j < 20; j++) {
                    aln_matrix2get(n_distmp, i, j) *= 600.0 / average;
                }
            }
        }

        for (int32_t i = 0; i < 20; i++) {
            for (int32_t j = 0; j < 20; j++) {
                aln_matrix2get(n_distmp, i, j) -= ctx->offset;
            }
        }

        for (int32_t i = 0; i < 20; i++) {
            for (int32_t j = 0; j < 20; j++) {
                int32_t out = 0;
                {
                    int32_t in = aln_matrix2get(n_distmp, i, j);
                    if (in > 0.0)
                        out = ((int32_t)(in + 0.5));
                    else if (in == 0.0)
                        out = (0);
                    else if (in < 0.0)
                        out = ((int32_t)(in - 0.5));
                    else
                        out = 0;
                }
                aln_matrix2get(n_distmp, i, j) = out;
            }
        }

        ctx->n_dis = aln_arenaAllocMatrix2I32(tempArena, ctx->nalphabets, ctx->nalphabets);
        for (int32_t i = 0; i < 20; i++) {
            for (int32_t j = 0; j < 20; j++) {
                aln_matrix2get(ctx->n_dis, i, j) = (int32_t)aln_matrix2get(n_distmp, i, j);
            }
        }
    }

    {
        int32_t charsize = 128;
        for (int32_t i = 0; i < charsize; i++) {
            ctx->amino_n[i] = -1;
        }
        for (int32_t i = 0; i < ctx->nalphabets; i++) {
            ctx->amino_n[(int32_t)ctx->amino[i]] = i;
        }
        ctx->n_dis_consweight_multi = aln_arenaAllocMatrix2F32(tempArena, ctx->nalphabets, ctx->nalphabets);
        for (int32_t i = 0; i < ctx->nalphabets; i++) {
            for (int32_t j = 0; j < ctx->nalphabets; j++) {
                aln_matrix2get(ctx->n_dis_consweight_multi, i, j) = (float)aln_matrix2get(ctx->n_dis, i, j) * ctx->consweight_multi;
            }
        }
    }

    // TODO(sen) Sequence verification?

    int***   topol = AllocateIntCub(ctx->njob, 2, 0);
    double** len = AllocateFloatMtx(ctx->njob, 2);
    Treedep* dep = (Treedep*)calloc(ctx->njob, sizeof(Treedep));
    {
        double** iscore = AllocateFloatHalfMtx(ctx->njob);
        fixed_musclesupg_double_realloc_nobk_halfmtx_memsave(ctx, iscore, topol, len, dep);
    }

    int** localmem = AllocateIntMtx(2, ctx->njob + 1);
    localmem[0][0] = -1;
    topolorderz(localmem[0], topol, dep, ctx->njob - 2, 2);

    double* eff = AllocateDoubleVec(ctx->njob);
    counteff_simple_double_nostatic_memsave(ctx->njob, topol, len, dep, eff);
    for (int i = ctx->njob; i < ctx->njob; i++) {
        eff[i] /= (double)100;
    }

    // TODO(sen) Figure out what to do when running out here
    int    alloclen = maxInputSeqLen * 4 + 1;
    char** bseq = AllocateCharMtx(ctx->njob, alloclen);

    for (int i = 0; i < ctx->njob; i++) {
        copyWithNoGaps(bseq[i], seq[i]);
    }

    {
        int       m1, m2;
        double    dumdb = 0.0;
        double*** cpmxchild0 = NULL;
        double*** cpmxchild1 = NULL;
        double    orieff1, orieff2;

        char* swaplist = 0;

        int*    fftlog = AllocateIntVec(ctx->njob);
        double* effarr1 = AllocateDoubleVec(ctx->njob);
        double* effarr2 = AllocateDoubleVec(ctx->njob);
        char*   indication1 = AllocateCharVec(150);
        char*   indication2 = AllocateCharVec(150);

        double**  dynamicmtx = AllocateDoubleMtx(ctx->nalphabets, ctx->nalphabets);
        int**     localmem = calloc(sizeof(int*), 2);
        double*** cpmxhist = (double***)calloc(ctx->njob - 1, sizeof(double**));
        int**     memhist = (int**)calloc(ctx->njob - 1, sizeof(int*));

        aln_LocalHom*** localhomshrink = (aln_LocalHom***)calloc(ctx->njob, sizeof(aln_LocalHom**));
        for (int i = 0; i < ctx->njob; i++) {
            localhomshrink[i] = (aln_LocalHom**)calloc(ctx->njob, sizeof(aln_LocalHom*));
        }

        int* alreadyaligned = AllocateIntVec(ctx->njob);
        for (int i = 0; i < ctx->njob; i++)
            alreadyaligned[i] = 1;
        for (int i = ctx->njob; i < ctx->njob; i++)
            alreadyaligned[i] = 0;

        for (int l = 0; l < ctx->njob; l++)
            fftlog[l] = 1;

        calcimportance_half(ctx->njob, eff, bseq, localhomtable, alloclen);

        char** mseq1 = AllocateCharMtx(ctx->njob, 0);
        char** mseq2 = AllocateCharMtx(ctx->njob, 0);

        int clus1 = 0;
        int clus2 = 0;

        for (int l = 0; l < ctx->njob - 1; l++) {
            m1 = topol[l][0][0];
            m2 = topol[l][1][0];

            if (dep[l].child0 == -1) {
                cpmxchild0 = NULL;
            } else {
                cpmxchild0 = cpmxhist + dep[l].child0;
            }

            if (dep[l].child1 == -1) {
                cpmxchild1 = NULL;
            } else {
                cpmxchild1 = cpmxhist + dep[l].child1;
            }

            if (dep[l].child0 == -1) {
                localmem[0] = calloc(sizeof(int), 2);
                localmem[0][0] = m1;
                localmem[0][1] = -1;
                clus1 = 1;
            } else {
                localmem[0] = memhist[dep[l].child0];
                clus1 = intlen(localmem[0]);
            }
            if (dep[l].child1 == -1) {
                localmem[1] = calloc(sizeof(int), 2);
                localmem[1][0] = m2;
                localmem[1][1] = -1;
                clus2 = 1;
            } else {
                localmem[1] = memhist[dep[l].child1];
                clus2 = intlen(localmem[1]);
            }

            if (l != ctx->njob - 2) {
                memhist[l] = calloc(sizeof(int), clus1 + clus2 + 1);
                intcpy(memhist[l], localmem[0]);
                intcpy(memhist[l] + clus1, localmem[1]);
                memhist[l][clus1 + clus2] = -1;
            }

            makedynamicmtx(ctx, dynamicmtx, dep[l].distfromtip);

            int len1 = strlen(bseq[m1]);
            int len2 = strlen(bseq[m2]);

            // TODO(sen) Out of memory error?
            aln_assert(alloclen >= len1 + len2);

            clus1 = fastconjuction_noname(localmem[0], bseq, mseq1, effarr1, eff, indication1, ctx->minimumweight, &orieff1);
            clus2 = fastconjuction_noname(localmem[1], bseq, mseq2, effarr2, eff, indication2, ctx->minimumweight, &orieff2);

            fastshrinklocalhom_half(localmem[0], localmem[1], localhomtable, localhomshrink);

            imp_match_init_strict(
                ctx,
                clus1,
                clus2,
                strlen(mseq1[0]),
                strlen(mseq2[0]),
                mseq1,
                mseq2,
                effarr1,
                effarr2,
                localhomshrink,
                swaplist,
                localmem[0],
                localmem[1]
            );

            A__align(
                ctx,
                dynamicmtx,
                ctx->penalty,
                ctx->penalty_ex,
                mseq1,
                mseq2,
                effarr1,
                effarr2,
                clus1,
                clus2,
                alloclen,
                ctx->constraint,
                &dumdb,
                NULL,
                NULL,
                NULL,
                NULL,
                ctx->outgap,
                ctx->outgap,
                localmem[0][0],
                1,
                cpmxchild0,
                cpmxchild1,
                cpmxhist + l,
                orieff1,
                orieff2
            );
        }
    }

    aln_Str* alignedSeqs = aln_arenaAllocArray(permArena, aln_Str, stringsCount);
    for (int32_t strIndex = 0; strIndex < stringsCount; strIndex++) {
        aln_Str alignedSeq = {bseq[strIndex], strlen(bseq[strIndex])};
        char*   permCopy = aln_arenaAllocArray(permArena, char, alignedSeq.len);
        for (int32_t charIndex = 0; charIndex < alignedSeq.len; charIndex++) {
            permCopy[charIndex] = alignedSeq.ptr[charIndex];
        }
        alignedSeqs[strIndex] = (aln_Str) {permCopy, alignedSeq.len};
    }

    aln_AlignResult result = {.seqs = alignedSeqs, .seqCount = stringsCount, .bytesWritten = permArena->used};
    return result;
}
