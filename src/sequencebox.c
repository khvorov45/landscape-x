#include <R.h>
#include <Rinternals.h>
#include <R_ext/Random.h>

// clang-format off
#define aln_IMPLEMENTATION
#define aln_PUBLICAPI static
#define aln_assert(cond) do {if (cond) {} else {error("assertion failure at %s:%d", __FILE__, __LINE__);}} while (0)
#include "align.h"
// clang-format on

static aln_StrArray
rStrArrayToAlnStrArray(SEXP strs) {
    intptr_t strsCount = LENGTH(strs);
    aln_Str* alnStrings = (aln_Str*)R_alloc(strsCount, sizeof(aln_Str));
    for (int strIndex = 0; strIndex < strsCount; strIndex++) {
        SEXP rSeq = STRING_ELT(strs, strIndex);
        alnStrings[strIndex] = (aln_Str) {(char*)CHAR(rSeq), LENGTH(rSeq)};
    }
    aln_StrArray result = {alnStrings, strsCount};
    return result;
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
alnStrFromRstrArray(SEXP rstr, intptr_t ind) {
    SEXP    rstr0 = STRING_ELT(rstr, ind);
    aln_Str result = {(char*)CHAR(rstr0), LENGTH(rstr0)};
    return result;
}

static aln_Arena
arenaFromR(intptr_t size) {
    void* mem = R_alloc(size, 1);
    aln_assert(mem);
    aln_Arena arena = {mem, size};
    return arena;
}

static intptr_t
intFromRNumber(char* name, SEXP rnumber) {
    assertLen(name, rnumber, 1);
    intptr_t result = 0;
    switch (TYPEOF(rnumber)) {
        case REALSXP: {
            double dbl = REAL_ELT(rnumber, 0);
            result = (intptr_t)dbl;
        } break;
        case INTSXP: {
            result = INTEGER_ELT(rnumber, 0);
        } break;
        default: error("%s should be a number", name); break;
    }
    return result;
}

SEXP
align_sequences_c(SEXP references, SEXP sequences, SEXP mode, SEXP matrices) {
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

    aln_StrArray alnRefs = rStrArrayToAlnStrArray(references);
    aln_StrArray alnStrings = rStrArrayToAlnStrArray(sequences);

    aln_AlignConfig alignConfig = {0};
    aln_fillAlignConfigLens(alnRefs, alnStrings, &alignConfig);
    alignConfig.maxGrid.ptr = R_alloc(alignConfig.maxGrid.len, 1);
    alignConfig.alignedSeqs.arr.ptr = (aln_Alignment*)R_alloc(alignConfig.alignedSeqs.arr.len, sizeof(aln_Alignment));
    alignConfig.alignedSeqs.data.ptr = R_alloc(alignConfig.alignedSeqs.data.len, 1);

    if (returnMatrices) {
        alignConfig.storedMatrices.arr.ptr = (aln_Matrix2NW*)R_alloc(alignConfig.storedMatrices.arr.len, sizeof(aln_Matrix2NW));
        alignConfig.storedMatrices.data.ptr = R_alloc(alignConfig.storedMatrices.data.len, 1);
    }

    aln_align(alnRefs, alnStrings, &alignConfig);

    if (reconstructToCommon) {
        aln_ReconstructToCommonConfig reconstructConfig = {0};
        aln_Str commonRef = alnRefs.ptr[0];
        aln_ReferenceGaps refGaps = {(intptr_t*)R_alloc(commonRef.len + 1, sizeof(intptr_t)), commonRef.len + 1};
        aln_fillReconstructToCommonRefConfigLens(alignConfig.alignedSeqs.arr, commonRef, alnStrings, refGaps, &reconstructConfig);
        reconstructConfig.seqs.ptr = (aln_Str*)R_alloc(reconstructConfig.seqs.len, sizeof(aln_Str));
        reconstructConfig.data.ptr = R_alloc(reconstructConfig.data.len, 1);

        aln_reconstructToCommonRef(alignConfig.alignedSeqs.arr, commonRef, alnStrings, refGaps, &reconstructConfig);

        SET_STRING_ELT(resultRefs, 0, mkCharLen(reconstructConfig.commonRef.ptr, (int)reconstructConfig.commonRef.len));
        for (int seqIndex = 0; seqIndex < sequencesCount; seqIndex++) {
            aln_Str alnStr = reconstructConfig.seqs.ptr[seqIndex];
            SET_STRING_ELT(resultSeqs, seqIndex, mkCharLen(alnStr.ptr, (int)alnStr.len));
        }
    } else {
        aln_Arena reconstructArena = {.size = 10 * 1024 * 1024};
        reconstructArena.base = R_alloc(reconstructArena.size, 1);
        for (int seqIndex = 0; seqIndex < alignConfig.alignedSeqs.arr.len; seqIndex++) {
            aln_Alignment alignment = alignConfig.alignedSeqs.arr.ptr[seqIndex];
            aln_Str ogstr = alnStrings.ptr[seqIndex];
            aln_Str thisRef = alnRefs.ptr[referenceCount > 1 ? seqIndex : 0];

            aln_Str alnStr = aln_reconstruct(alignment, aln_Reconstruct_Str, thisRef, ogstr, &reconstructArena);
            SET_STRING_ELT(resultSeqs, seqIndex, mkCharLen(alnStr.ptr, (int)alnStr.len));
            reconstructArena.used = 0;

            aln_Str alnRef = aln_reconstruct(alignment, aln_Reconstruct_Ref, thisRef, ogstr, &reconstructArena);
            SET_STRING_ELT(resultRefs, seqIndex, mkCharLen(alnRef.ptr, (int)alnRef.len));
            reconstructArena.used = 0;
        }
    }

    if (returnMatrices) {
        for (int matIndex = 0; matIndex < alignConfig.storedMatrices.arr.len; matIndex++) {
            aln_Matrix2NW mat = alignConfig.storedMatrices.arr.ptr[matIndex];
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

    UNPROTECT(2 + resultEntryCount);
    return result;
}

SEXP
generate_random_sequence_c(SEXP src, SEXP len) {
    assertLenAndType("src", src, 1, STRSXP, "string");

    int       lenInt = intFromRNumber("length", len);
    aln_Str   alnSrc = alnStrFromRstrArray(src, 0);
    aln_Arena arena = arenaFromR(lenInt);
    aln_Rng   rng = alnRngFromR();

    aln_Str alnStr = aln_randomString(&rng, alnSrc, lenInt, &arena);
    SEXP    result = allocVector(STRSXP, 1);
    SET_STRING_ELT(result, 0, mkCharLen(alnStr.ptr, (int)alnStr.len));
    return result;
}

SEXP
random_sequence_mod_c(
    SEXP src,
    SEXP n_mods_per_seq,
    SEXP trim_start_max,
    SEXP trim_end_max,
    SEXP mutation,
    SEXP deletion,
    SEXP insertion,
    SEXP insertion_src
) {
    assertType("src", src, STRSXP, "string");
    assertLen("n_mods_per_seq", n_mods_per_seq, 1);
    assertLenAndType("trim_start_max", trim_start_max, 1, REALSXP, "real number");
    assertLenAndType("trim_end_max", trim_end_max, 1, REALSXP, "real number");
    assertLenAndType("mutation", mutation, 1, REALSXP, "real number");
    assertLenAndType("deletion", deletion, 1, REALSXP, "real number");
    assertLenAndType("insertion", insertion, 1, REALSXP, "real number");
    assertLenAndType("insertion_src", insertion_src, 1, STRSXP, "string");

    // TODO(sen) Better memory determination
    aln_Rng      rng = alnRngFromR();
    aln_Arena    arena = arenaFromR(20 * 1024 * 1024);
    aln_StrArray alnStrings = rStrArrayToAlnStrArray(src);
    aln_Str      alnInsertionSrc = alnStrFromRstrArray(insertion_src, 0);

    aln_StrModSpec spec = {
        .trimStartMaxProp = REAL_ELT(trim_start_max, 0),
        .trimEndMaxProp = REAL_ELT(trim_end_max, 0),
        .mutationProb = REAL_ELT(mutation, 0),
        .deletionProb = REAL_ELT(deletion, 0),
        .insertionProb = REAL_ELT(insertion, 0),
        .insertionSrc = alnInsertionSrc,
    };

    intptr_t modsPerSeq = intFromRNumber("n_mods_per_seq", n_mods_per_seq);

    SEXP result = PROTECT(allocVector(STRSXP, alnStrings.len * modsPerSeq));
    for (intptr_t ind = 0; ind < alnStrings.len; ind++) {
        for (intptr_t modInd = 0; modInd < modsPerSeq; modInd++) {
            aln_Str alnSrc = alnStrings.ptr[ind];
            aln_Str alnMod = aln_randomStringMod(&rng, alnSrc, spec, &arena);
            SET_STRING_ELT(result, ind * modsPerSeq + modInd, mkCharLen(alnMod.ptr, (int)alnMod.len));
        }
    }

    UNPROTECT(1);
    return result;
}
