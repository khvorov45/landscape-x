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
align_c(SEXP references, SEXP sequences) {
    int referenceCount = LENGTH(references);
    int sequencesCount = LENGTH(sequences);
    if (referenceCount != 1 && referenceCount != sequencesCount) {
        error("reference count (%d) should be either 1 or the same as sequence count (%d)", referenceCount, sequencesCount);
    }

    SEXP result = PROTECT(allocVector(STRSXP, sequencesCount));

    intptr_t  totalMemoryBytes = 20 * 1024 * 1024;
    aln_Memory alnMem = aln_createMemory(R_alloc(totalMemoryBytes, 1), totalMemoryBytes, totalMemoryBytes / 4);
    if (alnMem.perm.base) {
        aln_Str* alnRefs = rStrArrayToAlnStrArray(&alnMem.perm, references, referenceCount);
        aln_Str* alnStrings = rStrArrayToAlnStrArray(&alnMem.perm, sequences, sequencesCount);

        aln_AlignResult alnResult = aln_align(
            (aln_StrArray) {alnRefs, referenceCount},
            (aln_StrArray) {alnStrings, sequencesCount},
            (aln_Config) { .storeFinalMatrices = false},
            &alnMem
        );

        if (alnResult.alignments.len == sequencesCount) {
            for (int seqIndex = 0; seqIndex < alnResult.alignments.len; seqIndex++) {
                aln_Alignment alignment = alnResult.alignments.ptr[seqIndex];
                aln_Str thisRef = alnRefs[0];
                if (referenceCount > 1) {
                    thisRef = alnRefs[seqIndex];
                }
                aln_Str       alnStr = aln_reconstruct(alignment, aln_Reconstruct_Str, thisRef, alnStrings[seqIndex], &alnMem.perm);
                SET_STRING_ELT(result, seqIndex, mkCharLen(alnStr.ptr, (int)alnStr.len));
            }
        } else {
            error("unexpected alignment result");
        }
    } else {
        error("could not allocate memory");
    }

    UNPROTECT(1);
    return result;
}
