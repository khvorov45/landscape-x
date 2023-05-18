#include <R.h>
#include <Rinternals.h>

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

SEXP
align_c(SEXP references, SEXP sequences, SEXP mode) {
    bool reconstructToCommon = false;
    {
        if (LENGTH(mode) != 1) {
            error("mode should be of length 1, not length %d", LENGTH(mode));
        }
        SEXP    mode0 = STRING_ELT(mode, 0);
        aln_Str modechar = {(char*)CHAR(mode0), LENGTH(mode0)};
        bool    indiv = aln_streq(modechar, aln_STR("individual"));
        reconstructToCommon = aln_streq(modechar, aln_STR("common"));
        if (!indiv && !reconstructToCommon) {
            error("mode should be either 'individual' or 'common'");
        }
    }

    int referenceCount = LENGTH(references);
    int sequencesCount = LENGTH(sequences);

    if (reconstructToCommon && referenceCount != 1) {
        error("reference count should be 1, not %d", referenceCount);
    }

    if (referenceCount != 1 && referenceCount != sequencesCount) {
        error("reference count (%d) should be either 1 or the same as sequence count (%d)", referenceCount, sequencesCount);
    }

    SEXP resultSeqs = PROTECT(allocVector(STRSXP, sequencesCount));
    SEXP resultRefs = 0;
    if (reconstructToCommon) {
        resultRefs = PROTECT(allocVector(STRSXP, 1));
    } else {
        resultRefs = PROTECT(allocVector(STRSXP, sequencesCount));
    }

    SEXP result = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(result, 0, resultSeqs);
    SET_VECTOR_ELT(result, 1, resultRefs);

    SEXP resultNames = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(resultNames, 0, mkChar("sequences"));
    SET_STRING_ELT(resultNames, 1, mkChar("references"));
    setAttrib(result, R_NamesSymbol, resultNames);

    intptr_t   totalMemoryBytes = 20 * 1024 * 1024;
    aln_Memory alnMem = aln_createMemory(R_alloc(totalMemoryBytes, 1), totalMemoryBytes, totalMemoryBytes / 4);
    if (alnMem.perm.base) {
        aln_Str* alnRefs = rStrArrayToAlnStrArray(&alnMem.perm, references, referenceCount);
        aln_Str* alnStrings = rStrArrayToAlnStrArray(&alnMem.perm, sequences, sequencesCount);

        aln_AlignResult alnResult = aln_align(
            (aln_StrArray) {alnRefs, referenceCount},
            (aln_StrArray) {alnStrings, sequencesCount},
            (aln_Config) {.storeFinalMatrices = false},
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
        } else {
            error("unexpected alignment result");
        }
    } else {
        error("could not allocate memory");
    }

    UNPROTECT(4);
    return result;
}
