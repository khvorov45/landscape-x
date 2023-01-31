#include "cbuild.h"

#define aln_assert(condition) prb_assert(condition)
#define aln_PUBLICAPI static
#define aln_IMPLEMENTATION
#include "../src/align.h"

#define function static

typedef uint64_t u64;
typedef uint32_t u32;
typedef int32_t  i32;
typedef uint8_t  u8;

typedef struct GenerateSequencesResult {
    prb_Str  full;
    prb_Str* seqs;
    i32      seqCount;
} GenerateSequencesResult;

function GenerateSequencesResult
generateSequences(prb_Arena* arena, prb_Rng* rng, char* choices, i32 choicesCount, i32 lengthFullSeq, i32 seqCount) {
    char* fullSeqBuf = prb_arenaAllocArray(arena, char, lengthFullSeq + 1);
    for (i32 aaIndex = 0; aaIndex < lengthFullSeq; aaIndex++) {
        u32  choiceIndex = prb_randomU32Bound(rng, choicesCount);
        char choice = choices[choiceIndex];
        fullSeqBuf[aaIndex] = choice;
    }
    fullSeqBuf[lengthFullSeq] = '\0';
    prb_Str fullSeq = {fullSeqBuf, lengthFullSeq};

    i32  maxTrimFromEnds = lengthFullSeq / 10;
    i32  maxMutations = lengthFullSeq / 10;
    i32* mutationBuffer = prb_arenaAllocArray(arena, i32, maxMutations);

    prb_Str* seqs = prb_arenaAllocArray(arena, prb_Str, seqCount);
    for (i32 seqIndex = 0; seqIndex < seqCount; seqIndex++) {
        i32 startIndex = prb_randomU32Bound(rng, maxTrimFromEnds + 1);
        i32 onePastEndIndex = lengthFullSeq - prb_randomU32Bound(rng, maxTrimFromEnds + 1);
        i32 seqLen = onePastEndIndex - startIndex;
        i32 mutationCount = prb_randomU32Bound(rng, maxMutations + 1);
        for (i32 mutationIndex = 0; mutationIndex < mutationCount; mutationIndex++) {
            mutationBuffer[mutationIndex] = prb_randomU32Bound(rng, seqLen);
        }

        char* seqBuf = prb_arenaAllocArray(arena, char, onePastEndIndex - startIndex + 1);
        for (i32 seqIndex = 0; seqIndex < seqLen; seqIndex++) {
            bool isMutated = false;
            for (i32 mutationIndex = 0; mutationIndex < mutationCount && !isMutated; mutationIndex++) {
                if (seqIndex == mutationBuffer[mutationIndex]) {
                    isMutated = true;
                }
            }

            if (isMutated) {
                u32  choiceIndex = prb_randomU32Bound(rng, choicesCount);
                char choice = choices[choiceIndex];
                seqBuf[seqIndex] = choice;
            } else {
                seqBuf[seqIndex] = fullSeqBuf[seqIndex + startIndex];
            }
        }
        seqBuf[seqLen] = '\0';

        i32 xLocation = prb_randomU32Bound(rng, seqLen);
        seqBuf[xLocation] = 'X';

        seqs[seqIndex] = (prb_Str) {seqBuf, seqLen};
    }

    GenerateSequencesResult result = {.full = fullSeq, .seqs = seqs, .seqCount = seqCount};
    return result;
}

int
main() {
    prb_TimeStart testsStart = prb_timeStart();
    prb_Arena     arena_ = prb_createArenaFromVmem(1 * prb_GIGABYTE);
    prb_Arena*    arena = &arena_;
    prb_Rng       rng_ = prb_createRng(1);
    prb_Rng*      rng = &rng_;

    GenerateSequencesResult genSeq = {};
    {
        char aminoAcids[] = {'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'};
        i32  aminoAcidsCount = prb_arrayCount(aminoAcids);
        prb_assert(aminoAcidsCount == 20);
        genSeq = generateSequences(arena, rng, aminoAcids, aminoAcidsCount, 100, 8);
        prb_writelnToStdout(arena, genSeq.full);
        for (i32 seqIndex = 0; seqIndex < genSeq.seqCount; seqIndex++) {
            prb_writelnToStdout(arena, genSeq.seqs[seqIndex]);
        }
    }

    aln_AlignResult alnResult = aln_align((aln_Str*)genSeq.seqs, genSeq.seqCount, prb_arenaFreePtr(arena), 20 * prb_MEGABYTE);
    prb_arenaChangeUsed(arena, alnResult.bytesWritten);

    prb_writelnToStdout(arena, prb_fmt(arena, "tests took %.2fms", prb_getMsFrom(testsStart)));
    return 0;
}
