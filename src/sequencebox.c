#include <R.h>
#include <Rinternals.h>
#include <R_ext/Random.h>

// clang-format off
#define aln_IMPLEMENTATION
#define aln_PUBLICAPI static
#define aln_assert(cond) do {if (cond) {} else {error("assertion failure at %s:%d", __FILE__, __LINE__);}} while (0)
#include "align.h"
// clang-format on

static aln_Str*
rStrArrayToAlnStrArray(aln_Arena* arena, SEXP strs, int strsCount) {
    aln_Str* alnStrings = aln_arenaAllocArray(arena, aln_Str, strsCount);
    for (int strIndex = 0; strIndex < strsCount; strIndex++) {
        SEXP rSeq = STRING_ELT(strs, strIndex);
        alnStrings[strIndex] = (aln_Str) {(char*)CHAR(rSeq), LENGTH(rSeq)};
    }
    return alnStrings;
}

static SEXP
createRMatrix(int rows, int cols, int type) {
    SEXP matrix = PROTECT(allocVector(type, rows * cols));
    SEXP dim = PROTECT(allocVector(INTSXP, 2));
    SET_INTEGER_ELT(dim, 0, rows);
    SET_INTEGER_ELT(dim, 1, cols);
    setAttrib(matrix, R_DimSymbol, dim);
    UNPROTECT(2);
    return matrix;
}

static void
assertLen(char* name, SEXP obj, int len) {
    if (LENGTH(obj) != len) {
        error("%s should be length 1, not length %d", name, LENGTH(obj));
    }
}

static void
assertType(char* name, SEXP obj, int type, char* typename) {
    if (TYPEOF(obj) != type) {
        error("%s should be of type %s", name, typename);
    }
}

static void
assertLenAndType(char* name, SEXP obj, int len, int type, char* typename) {
    assertLen(name, obj, len);
    assertType(name, obj, type, typename);
}

static aln_Rng
alnRngFromR(void) {
    GetRNGstate();
    aln_Rng rng = aln_createRng((uint32_t)(unif_rand() * (double)UINT32_MAX));
    PutRNGstate();
    return rng;
}

static aln_Str
alnStrFromRstrArray(SEXP rstr) {
    SEXP    rstr0 = STRING_ELT(rstr, 0);
    aln_Str result = {(char*)CHAR(rstr0), LENGTH(rstr0)};
    return result;
}

SEXP
align_c(SEXP references, SEXP sequences, SEXP mode, SEXP matrices) {
    assertType("references", references, STRSXP, "string");
    assertType("sequences", sequences, STRSXP, "string");
    assertLenAndType("mode", mode, 1, STRSXP, "string");
    assertLenAndType("matrices", matrices, 1, LGLSXP, "boolean");

    int referenceCount = LENGTH(references);
    int sequencesCount = LENGTH(sequences);

    bool reconstructToCommon = false;
    {
        SEXP    mode0 = STRING_ELT(mode, 0);
        aln_Str modechar = {(char*)CHAR(mode0), LENGTH(mode0)};
        bool    indiv = aln_streq(modechar, aln_STR("individual"));
        reconstructToCommon = aln_streq(modechar, aln_STR("common"));
        if (!indiv && !reconstructToCommon) {
            error("mode should be either 'individual' or 'common'");
        }
    }

    if (reconstructToCommon) {
        assertLen("references", references, 1);
    }

    if (referenceCount != 1 && referenceCount != sequencesCount) {
        error("reference count (%d) should be either 1 or the same as sequence count (%d)", referenceCount, sequencesCount);
    }

    bool returnMatrices = LOGICAL_ELT(matrices, 0);

    SEXP resultSeqs = PROTECT(allocVector(STRSXP, sequencesCount));
    SEXP resultRefs = 0;
    if (reconstructToCommon) {
        resultRefs = PROTECT(allocVector(STRSXP, 1));
    } else {
        resultRefs = PROTECT(allocVector(STRSXP, sequencesCount));
    }

    SEXP resultScores = 0;
    SEXP resultDirections = 0;
    if (returnMatrices) {
        resultScores = PROTECT(allocVector(VECSXP, sequencesCount));
        resultDirections = PROTECT(allocVector(VECSXP, sequencesCount));
    }

    int  resultEntryCount = 2 + ((int)returnMatrices * 2);
    SEXP result = PROTECT(allocVector(VECSXP, resultEntryCount));
    SET_VECTOR_ELT(result, 0, resultSeqs);
    SET_VECTOR_ELT(result, 1, resultRefs);

    SEXP resultNames = PROTECT(allocVector(STRSXP, resultEntryCount));
    SET_STRING_ELT(resultNames, 0, mkChar("sequences"));
    SET_STRING_ELT(resultNames, 1, mkChar("references"));
    if (returnMatrices) {
        SET_VECTOR_ELT(result, 2, resultScores);
        SET_VECTOR_ELT(result, 3, resultDirections);
        SET_STRING_ELT(resultNames, 2, mkChar("scores"));
        SET_STRING_ELT(resultNames, 3, mkChar("directions"));
    }
    setAttrib(result, R_NamesSymbol, resultNames);

    intptr_t   totalMemoryBytes = 20 * 1024 * 1024;
    aln_Memory alnMem = aln_createMemory(R_alloc(totalMemoryBytes, 1), totalMemoryBytes, totalMemoryBytes / 4);
    if (alnMem.perm.base) {
        aln_Str* alnRefs = rStrArrayToAlnStrArray(&alnMem.perm, references, referenceCount);
        aln_Str* alnStrings = rStrArrayToAlnStrArray(&alnMem.perm, sequences, sequencesCount);

        aln_AlignResult alnResult = aln_align(
            (aln_StrArray) {alnRefs, referenceCount},
            (aln_StrArray) {alnStrings, sequencesCount},
            (aln_Config) {.storeFinalMatrices = returnMatrices},
            &alnMem
        );

        if (alnResult.alignments.len == sequencesCount) {
            if (reconstructToCommon) {
                aln_ReconstructToCommonRefResult reconstructResult = aln_reconstructToCommonRef(alnResult.alignments, alnRefs[0], (aln_StrArray) {alnStrings, sequencesCount}, &alnMem);
                SET_STRING_ELT(resultRefs, 0, mkCharLen(reconstructResult.commonRef.ptr, (int)reconstructResult.commonRef.len));
                aln_assert(reconstructResult.alignedStrs.len == sequencesCount);
                for (int seqIndex = 0; seqIndex < sequencesCount; seqIndex++) {
                    aln_Str alnStr = reconstructResult.alignedStrs.ptr[seqIndex];
                    SET_STRING_ELT(resultSeqs, seqIndex, mkCharLen(alnStr.ptr, (int)alnStr.len));
                }
            } else {
                for (int seqIndex = 0; seqIndex < alnResult.alignments.len; seqIndex++) {
                    aln_Alignment alignment = alnResult.alignments.ptr[seqIndex];
                    aln_Str       thisRef = alnRefs[0];
                    if (referenceCount > 1) {
                        thisRef = alnRefs[seqIndex];
                    }
                    aln_Str alnStr = aln_reconstruct(alignment, aln_Reconstruct_Str, thisRef, alnStrings[seqIndex], &alnMem.perm);
                    SET_STRING_ELT(resultSeqs, seqIndex, mkCharLen(alnStr.ptr, (int)alnStr.len));
                    aln_Str alnRef = aln_reconstruct(alignment, aln_Reconstruct_Ref, thisRef, alnStrings[seqIndex], &alnMem.perm);
                    SET_STRING_ELT(resultRefs, seqIndex, mkCharLen(alnRef.ptr, (int)alnRef.len));
                }
            }

            if (returnMatrices) {
                aln_assert(alnResult.alignments.len == alnResult.matrices.len);
                for (int matIndex = 0; matIndex < alnResult.matrices.len; matIndex++) {
                    aln_Matrix2NW mat = alnResult.matrices.ptr[matIndex];
                    int           rwidth = mat.height;
                    int           rheight = mat.width;
                    SEXP          rmatScores = PROTECT(createRMatrix(rwidth, rheight, REALSXP));
                    SEXP          rmatDirs = PROTECT(createRMatrix(rwidth, rheight, INTSXP));
                    for (int row = 0; row < mat.height; row++) {
                        for (int col = 0; col < mat.width; col++) {
                            aln_NWEntry entry = aln_matrix2get(mat, row, col);
                            int         rrow = col;
                            int         rcol = row;
                            int         rindex = rrow * rwidth + rcol;
                            SET_REAL_ELT(rmatScores, rindex, entry.score);
                            SET_INTEGER_ELT(rmatDirs, rindex, (int)entry.cameFromDir);
                        }
                    }
                    SET_VECTOR_ELT(resultScores, matIndex, rmatScores);
                    SET_VECTOR_ELT(resultDirections, matIndex, rmatDirs);
                    UNPROTECT(2);
                }
            }
        } else {
            error("unexpected alignment result");
        }
    } else {
        error("could not allocate memory");
    }

    UNPROTECT(2 + resultEntryCount);
    return result;
}

SEXP
generate_random_sequence_c(SEXP src, SEXP len) {
    assertLenAndType("src", src, 1, STRSXP, "string");
    assertLen("len", len, 1);

    int lenInt = 0;
    switch (TYPEOF(len)) {
        case REALSXP: {
            double lenDouble = REAL_ELT(len, 0);
            lenInt = (int)lenDouble;
        } break;

        case INTSXP: {
            lenInt = INTEGER_ELT(len, 0);
        } break;

        default: error("length should be a number"); break;
    }

    aln_Str alnSrc = alnStrFromRstrArray(src);

    void*     mem = R_alloc(lenInt, 1);
    aln_Arena arena = {mem, lenInt};

    aln_Rng rng = alnRngFromR();

    aln_Str alnStr = aln_randomString(&rng, alnSrc, lenInt, &arena);
    SEXP    result = allocVector(STRSXP, 1);
    SET_STRING_ELT(result, 0, mkCharLen(alnStr.ptr, (int)alnStr.len));
    return result;
}

SEXP
random_sequence_mod_c(
    SEXP src,
    SEXP trim_start_max,
    SEXP trim_end_max,
    SEXP mutation,
    SEXP deletion,
    SEXP insertion,
    SEXP insertion_src
) {
    assertLenAndType("src", src, 1, STRSXP, "string");
    assertLenAndType("trim_start_max", trim_start_max, 1, REALSXP, "real number");
    assertLenAndType("trim_end_max", trim_end_max, 1, REALSXP, "real number");
    assertLenAndType("mutation", mutation, 1, REALSXP, "real number");
    assertLenAndType("deletion", deletion, 1, REALSXP, "real number");
    assertLenAndType("insertion", insertion, 1, REALSXP, "real number");
    assertLenAndType("insertion_src", insertion_src, 1, STRSXP, "string");

    aln_Rng rng = alnRngFromR();

    aln_Str alnSrc = alnStrFromRstrArray(src);
    aln_Str alnInsertionSrc = alnStrFromRstrArray(insertion_src);

    intptr_t  memsize = alnSrc.len * 2;
    void*     mem = R_alloc(memsize, 1);
    aln_Arena arena = {mem, memsize};

    aln_StrModSpec spec = {
        .trimStartMaxProp = REAL_ELT(trim_start_max, 0),
        .trimEndMaxProp = REAL_ELT(trim_end_max, 0),
        .mutationProb = REAL_ELT(mutation, 0),
        .deletionProb = REAL_ELT(deletion, 0),
        .insertionProb = REAL_ELT(insertion, 0),
        .insertionSrc = alnInsertionSrc,
    };

    aln_Str alnMod = aln_randomStringMod(&rng, alnSrc, spec, &arena);
    SEXP    result = allocVector(STRSXP, 1);
    SET_STRING_ELT(result, 0, mkCharLen(alnMod.ptr, (int)alnMod.len));
    return result;
}
