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
rStrArrayToAlnStrArray(aln_Arena* arena, SEXP strs) {
    intptr_t strsCount = LENGTH(strs);
    aln_Str* alnStrings = aln_arenaAllocArray(arena, aln_Str, strsCount);
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
setRDfRownames(SEXP df, int len, aln_Arena* arena) {
    SEXP names = allocVector(STRSXP, len);
    setAttrib(df, R_RowNamesSymbol, names);
    for (intptr_t ind = 0; ind < len; ind++) {
        aln_TempMemory temp = aln_beginTempMemory(arena);

        char* ptr = aln_arenaFreePtr(arena);
        int   len = snprintf(ptr, aln_arenaFreeSize(arena), "%ld", ind + 1);
        SET_STRING_ELT(names, ind, mkCharLen(ptr, len));

        aln_endTempMemory(temp);
    }
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

static aln_Arena
arenaFromR(intptr_t size) {
    void* mem = R_alloc(size, 1);
    aln_assert(mem);
    aln_Arena arena = {mem, size};
    return arena;
}

static aln_Memory
memoryFromR(intptr_t size) {
    void* ptr = R_alloc(size, 1);
    aln_assert(ptr);
    aln_Memory mem = aln_createMemory(ptr, size, size / 4);
    return mem;
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

    // TODO(sen) Figure out how to work out memory size in a more robust way
    intptr_t     totalMemoryBytes = 20 * 1024 * 1024;
    aln_Memory   alnMem = memoryFromR(totalMemoryBytes);
    aln_StrArray alnRefs = rStrArrayToAlnStrArray(&alnMem.perm, references);
    aln_StrArray alnStrings = rStrArrayToAlnStrArray(&alnMem.perm, sequences);

    aln_AlignResult alnResult = aln_align(
        alnRefs,
        alnStrings,
        (aln_Config) {.storeFinalMatrices = returnMatrices},
        &alnMem
    );
    aln_assert(alnResult.alignments.len == sequencesCount);

    if (reconstructToCommon) {
        aln_ReconstructToCommonRefResult reconstructResult = aln_reconstructToCommonRef(alnResult.alignments, alnRefs.ptr[0], alnStrings, &alnMem);
        SET_STRING_ELT(resultRefs, 0, mkCharLen(reconstructResult.commonRef.ptr, (int)reconstructResult.commonRef.len));
        aln_assert(reconstructResult.alignedStrs.len == sequencesCount);
        for (int seqIndex = 0; seqIndex < sequencesCount; seqIndex++) {
            aln_Str alnStr = reconstructResult.alignedStrs.ptr[seqIndex];
            SET_STRING_ELT(resultSeqs, seqIndex, mkCharLen(alnStr.ptr, (int)alnStr.len));
        }
    } else {
        for (int seqIndex = 0; seqIndex < alnResult.alignments.len; seqIndex++) {
            aln_Alignment alignment = alnResult.alignments.ptr[seqIndex];
            aln_Str       thisRef = alnRefs.ptr[0];
            if (referenceCount > 1) {
                thisRef = alnRefs.ptr[seqIndex];
            }
            aln_Str alnStr = aln_reconstruct(alignment, aln_Reconstruct_Str, thisRef, alnStrings.ptr[seqIndex], &alnMem.perm);
            SET_STRING_ELT(resultSeqs, seqIndex, mkCharLen(alnStr.ptr, (int)alnStr.len));
            aln_Str alnRef = aln_reconstruct(alignment, aln_Reconstruct_Ref, thisRef, alnStrings.ptr[seqIndex], &alnMem.perm);
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

    UNPROTECT(2 + resultEntryCount);
    return result;
}

SEXP
generate_random_sequence_c(SEXP src, SEXP len) {
    assertLenAndType("src", src, 1, STRSXP, "string");

    int       lenInt = intFromRNumber("length", len);
    aln_Str   alnSrc = alnStrFromRstrArray(src);
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

    aln_Rng      rng = alnRngFromR();
    aln_Arena    arena = arenaFromR(20 * 1024 * 1024);
    aln_StrArray alnStrings = rStrArrayToAlnStrArray(&arena, src);
    aln_Str      alnInsertionSrc = alnStrFromRstrArray(insertion_src);

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

SEXP
create_tree_c(SEXP seqs) {
    assertType("seqs", seqs, STRSXP, "string");

    aln_Memory   mem = memoryFromR(1024 * 1024);
    aln_StrArray alnStrs = rStrArrayToAlnStrArray(&mem.perm, seqs);
    aln_Tree     alnTree = aln_createTree(alnStrs, &mem);

    SEXP result = PROTECT(allocVector(VECSXP, 2));
    SEXP nodeDf = allocVector(VECSXP, 2);
    SET_VECTOR_ELT(result, 0, nodeDf);
    SEXP branchDf = allocVector(VECSXP, 3);
    SET_VECTOR_ELT(result, 1, branchDf);

    {
        SEXP names = allocVector(STRSXP, 2);
        setAttrib(result, R_NamesSymbol, names);
        SET_STRING_ELT(names, 0, mkChar("node"));
        SET_STRING_ELT(names, 1, mkChar("branch"));
    }

    {
        SEXP names = allocVector(STRSXP, 1);
        classgets(nodeDf, names);
        SET_STRING_ELT(names, 0, mkChar("data.frame"));
        classgets(branchDf, names);
    }

    setRDfRownames(nodeDf, alnTree.nodes.len, &mem.temp);
    setRDfRownames(branchDf, alnTree.branches.len, &mem.temp);

    {
        SEXP names = allocVector(STRSXP, 2);
        setAttrib(nodeDf, R_NamesSymbol, names);
        SET_STRING_ELT(names, 0, mkChar("node"));
        SET_STRING_ELT(names, 1, mkChar("internal"));
    }

    {
        SEXP names = allocVector(STRSXP, 3);
        setAttrib(branchDf, R_NamesSymbol, names);
        SET_STRING_ELT(names, 0, mkChar("node1Index"));
        SET_STRING_ELT(names, 1, mkChar("node2Index"));
        SET_STRING_ELT(names, 2, mkChar("len"));
    }

    {
        SEXP nodeIndices = allocVector(INTSXP, alnTree.nodes.len);
        SET_VECTOR_ELT(nodeDf, 0, nodeIndices);
        SEXP nodeInternal = allocVector(LGLSXP, alnTree.nodes.len);
        SET_VECTOR_ELT(nodeDf, 1, nodeInternal);

        for (intptr_t ind = 0; ind < alnTree.nodes.len; ind++) {
            aln_TreeNode node = alnTree.nodes.ptr[ind];
            SET_INTEGER_ELT(nodeIndices, ind, ind);
            SET_LOGICAL_ELT(nodeInternal, ind, node.isInternal);
        }

        SEXP branchNode1Indices = allocVector(INTSXP, alnTree.branches.len);
        SET_VECTOR_ELT(branchDf, 0, branchNode1Indices);
        SEXP branchNode2Indices = allocVector(INTSXP, alnTree.branches.len);
        SET_VECTOR_ELT(branchDf, 1, branchNode2Indices);
        SEXP branchLens = allocVector(REALSXP, alnTree.branches.len);
        SET_VECTOR_ELT(branchDf, 2, branchLens);

        for (intptr_t ind = 0; ind < alnTree.branches.len; ind++) {
            aln_TreeBranch branch = alnTree.branches.ptr[ind];
            SET_INTEGER_ELT(branchNode1Indices, ind, branch.node1Index);
            SET_INTEGER_ELT(branchNode2Indices, ind, branch.node2Index);
            SET_REAL_ELT(branchLens, ind, branch.len);
        }
    }

    UNPROTECT(1);
    return result;
}
