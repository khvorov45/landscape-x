#include "cbuild.h"

// clang-format off
#define aln_IMPLEMENTATION
#define aln_PUBLICAPI static
#define aln_assert(cond) prb_assert(cond)
#include "../src/align.h"
// clang-format on

#define function static

typedef int32_t  i32;
typedef uint32_t u32;

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
    prb_Arena  arena_ = prb_createArenaFromVmem(1 * prb_GIGABYTE);
    prb_Arena* arena = &arena_;

    prb_Rng  rng_ = prb_createRng(1);
    prb_Rng* rng = &rng_;

    GenerateSequencesResult genSeq = {};
    {
        char aminoAcids[] = {'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'};
        i32  aminoAcidsCount = prb_arrayCount(aminoAcids);
        prb_assert(aminoAcidsCount == 20);
        genSeq = generateSequences(arena, rng, aminoAcids, aminoAcidsCount, 1000, 100);
        bool displayGenSeq = false;
        if (displayGenSeq) {
            prb_writelnToStdout(arena, genSeq.full);
            for (i32 seqIndex = 0; seqIndex < genSeq.seqCount; seqIndex++) {
                prb_writelnToStdout(arena, genSeq.seqs[seqIndex]);
            }
        }
    }

    prb_Str benchDir = prb_getParentDir(arena, prb_STR(__FILE__));

    prb_Str fastaOutputPath = prb_pathJoin(arena, benchDir, prb_STR("testseqs.fasta"));
    {
        prb_GrowingStr gstr = prb_beginStr(arena);
        for (i32 seqIndex = 0; seqIndex < genSeq.seqCount; seqIndex++) {
            prb_addStrSegment(&gstr, ">virus%d\n%.*s\n", seqIndex, prb_LIT(genSeq.seqs[seqIndex]));
        }
        prb_Str fastaContent = prb_endStr(&gstr);
        prb_writeEntireFile(arena, fastaOutputPath, fastaContent.ptr, fastaContent.len);
    }

    {
        prb_Str       mafftOutputPath = prb_pathJoin(arena, benchDir, prb_STR("mafft.out"));
        prb_Str       mafftCmd = prb_fmt(arena, "mafft --auto %.*s", prb_LIT(fastaOutputPath));
        prb_Process   proc = prb_createProcess(mafftCmd, (prb_ProcessSpec) {.redirectStdout = true, .redirectStderr = true, .stderrFilepath = mafftOutputPath, .stdoutFilepath = mafftOutputPath});
        prb_TimeStart mafftStart = prb_timeStart();
        prb_assert(prb_launchProcesses(arena, &proc, 1, prb_Background_No));
        prb_writeToStdout(prb_fmt(arena, "mafft total time: %.2fms\n", prb_getMsFrom(mafftStart)));
    }

    {
        prb_Arena     alnOut = prb_createArenaFromArena(arena, 20 * prb_MEGABYTE);
        prb_TimeStart start = prb_timeStart();
        aln_align(
            (aln_Str) {(char*)genSeq.full.ptr, genSeq.full.len},
            (aln_Str*)genSeq.seqs,
            genSeq.seqCount,
            (aln_Config) {
                .outmem = alnOut.base,
                .outmemBytes = alnOut.size,
                .tempmem = prb_arenaFreePtr(arena),
                .tempmemBytes = prb_arenaFreeSize(arena),
            }
        );
        prb_writeToStdout(prb_fmt(arena, "my total time: %.2fms\n", prb_getMsFrom(start)));
    }

    return 0;
}