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
    aln_Arena arena = {.base = R_alloc(totalMemoryBytes, 1), .size = totalMemoryBytes};
    if (arena.base) {
        aln_Str* alnRefs = rStrArrayToAlnStrArray(&arena, references, referenceCount);
        aln_Str* alnStrings = rStrArrayToAlnStrArray(&arena, sequences, sequencesCount);

        aln_Arena       alnOutput = aln_createArenaFromArena(&arena, aln_arenaFreeSize(&arena) / 4);
        aln_AlignResult alnResult = aln_align(
            alnRefs,
            referenceCount,
            alnStrings,
            sequencesCount,
            (aln_Config) {
                .outmem = alnOutput.base,
                .outmemBytes = alnOutput.size,
                .tempmem = aln_arenaFreePtr(&arena),
                .tempmemBytes = aln_arenaFreeSize(&arena),
            }
        );
        aln_arenaChangeUsed(&arena, alnOutput.size - alnResult.bytesWrittenToOutput);

        if (alnResult.strCount == sequencesCount) {
            for (int seqIndex = 0; seqIndex < alnResult.strCount; seqIndex++) {
                aln_Alignment alignment = alnResult.strs[seqIndex];
                aln_Str thisRef = alnRefs[0];
                if (referenceCount > 1) {
                    thisRef = alnRefs[seqIndex];
                }
                aln_Str       alnStr = aln_reconstruct(alignment, aln_Reconstruct_Str, thisRef, alnStrings[seqIndex], aln_arenaFreePtr(&arena), aln_arenaFreeSize(&arena));
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
