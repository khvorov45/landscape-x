#include "../programmable_build.h"

#define function static

typedef uint64_t u64;
typedef uint32_t u32;
typedef int32_t  i32;
typedef uint8_t  u8;

typedef struct Rng {
    u64 state;
    u64 inc;
} Rng;

function u32
getRandomU32(Rng* rng) {
    u64 oldstate = rng->state;
    rng->state = oldstate * 6364136223846793005ULL + rng->inc;
    u32 xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    u32 rot = oldstate >> 59u;
    u32 result = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
    return result;
}

function Rng
createRng(u64 initstate, u64 initseq) {
    Rng rng = {.state = 0, .inc = (initseq << 1u) | 1u};
    getRandomU32(&rng);
    rng.state += initstate;
    getRandomU32(&rng);
    return rng;
}

function u32
getRandomU32BiasedBound(Rng* rng, i32 bound) {
    u32 randU32 = getRandomU32(rng);
    u32 result = randU32 % bound;
    return result;
}

typedef struct GenerateSequencesResult {
    prb_String  full;
    prb_String* seqs;
    i32         seqCount;
} GenerateSequencesResult;

function GenerateSequencesResult
generateSequences(prb_Arena* arena, Rng* rng, char* choices, i32 choicesCount, i32 lengthFullSeq, i32 seqCount) {
    char* fullSeqBuf = prb_arenaAllocArray(arena, char, lengthFullSeq + 1);
    for (i32 aaIndex = 0; aaIndex < lengthFullSeq; aaIndex++) {
        u32  choiceIndex = getRandomU32BiasedBound(rng, choicesCount);
        char choice = choices[choiceIndex];
        fullSeqBuf[aaIndex] = choice;
    }
    fullSeqBuf[lengthFullSeq] = '\0';
    prb_String fullSeq = {fullSeqBuf, lengthFullSeq};

    i32  maxTrimFromEnds = lengthFullSeq / 10;
    i32  maxMutations = lengthFullSeq / 10;
    i32* mutationBuffer = prb_arenaAllocArray(arena, i32, maxMutations);

    prb_String* seqs = prb_arenaAllocArray(arena, prb_String, seqCount);
    for (i32 seqIndex = 0; seqIndex < seqCount; seqIndex++) {
        i32 startIndex = getRandomU32BiasedBound(rng, maxTrimFromEnds + 1);
        i32 onePastEndIndex = lengthFullSeq - getRandomU32BiasedBound(rng, maxTrimFromEnds + 1);
        i32 seqLen = onePastEndIndex - startIndex;
        i32 mutationCount = getRandomU32BiasedBound(rng, maxMutations + 1);
        for (i32 mutationIndex = 0; mutationIndex < mutationCount; mutationIndex++) {
            mutationBuffer[mutationIndex] = getRandomU32BiasedBound(rng, seqLen);
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
                u32  choiceIndex = getRandomU32BiasedBound(rng, choicesCount);
                char choice = choices[choiceIndex];
                seqBuf[seqIndex] = choice;
            } else {
                seqBuf[seqIndex] = fullSeqBuf[seqIndex + startIndex];
            }
        }
        seqBuf[seqLen] = '\0';
        seqs[seqIndex] = (prb_String) {seqBuf, seqLen};
    }

    GenerateSequencesResult result = {.full = fullSeq, .seqs = seqs, .seqCount = seqCount};
    return result;
}

function prb_String*
alignWithMafft(prb_Arena* arena, prb_String mafftExe, prb_String inputPath, prb_String mafftOuputPath) {
    {
        prb_String cmd = prb_fmt(arena, "%.*s --globalpair --maxiterate 1000 %.*s", prb_LIT(mafftExe), prb_LIT(inputPath));
        prb_writelnToStdout(cmd);
        prb_ProcessHandle proc = prb_execCmd(arena, cmd, prb_ProcessFlag_RedirectStdout, mafftOuputPath);
        prb_assert(proc.status == prb_ProcessStatus_CompletedSuccess);
    }

    prb_String* mafftAlignedSeqs = 0;
    {
        prb_ReadEntireFileResult mafftOutput = prb_readEntireFile(arena, mafftOuputPath);
        prb_assert(mafftOutput.success);
        prb_String       mafftContentLeft = prb_strFromBytes(mafftOutput.content);
        prb_LineIterator lineIter = prb_createLineIter(mafftContentLeft);
        while (prb_lineIterNext(&lineIter) == prb_Success) {
            prb_assert(lineIter.curLine.ptr[0] == '>');
            prb_assert(prb_lineIterNext(&lineIter) == prb_Success);
            arrput(mafftAlignedSeqs, lineIter.curLine);
        }
    }

    return mafftAlignedSeqs;
}

int
main() {
    prb_TimeStart testsStart = prb_timeStart();
    prb_Arena     arena_ = prb_createArenaFromVmem(1 * prb_GIGABYTE);
    prb_Arena*    arena = &arena_;
    Rng           rng_ = createRng(1, 3);
    Rng*          rng = &rng_;
    prb_String    testsDir = prb_getParentDir(arena, prb_STR(__FILE__));
    prb_String    rootDir = prb_getParentDir(arena, testsDir);

    GenerateSequencesResult genSeq = {};
    {
        char aminoAcids[] = {'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'};
        i32  aminoAcidsCount = prb_arrayLength(aminoAcids);
        prb_assert(aminoAcidsCount == 20);
        genSeq = generateSequences(arena, rng, aminoAcids, aminoAcidsCount, 10, 3);
        prb_writelnToStdout(genSeq.full);
        for (i32 seqIndex = 0; seqIndex < genSeq.seqCount; seqIndex++) {
            prb_writelnToStdout(genSeq.seqs[seqIndex]);
        }
    }

    prb_String fastaOutputPath = prb_pathJoin(arena, testsDir, prb_STR("testseqs.fasta"));
    {
        prb_GrowingString gstr = prb_beginString(arena);
        for (i32 seqIndex = 0; seqIndex < genSeq.seqCount; seqIndex++) {
            prb_addStringSegment(&gstr, ">virus%d\n%.*s\n", seqIndex, prb_LIT(genSeq.seqs[seqIndex]));
        }
        prb_String fastaContent = prb_endString(&gstr);
        prb_writeEntireFile(arena, fastaOutputPath, fastaContent.ptr, fastaContent.len);
    }

// TODO(sen) Need to set different env variables for global/local maffts
#if 0
    prb_String  mafftOuputPath = prb_pathJoin(arena, testsDir, prb_STR("mafft-testseqs.fasta"));
    prb_String* mafftAlignedSeqs = alignWithMafft(arena, prb_STR("mafft"), fastaOutputPath, mafftOuputPath);
    prb_assert(arrlen(mafftAlignedSeqs) == genSeq.seqCount);

    for (i32 seqIndex = 0; seqIndex < arrlen(mafftAlignedSeqs); seqIndex++) {
        prb_writelnToStdout(mafftAlignedSeqs[seqIndex]);
    }
#endif

    prb_String  localMafftExe = prb_pathJoin(arena, rootDir, prb_STR("mafft/core/mafft.tmpl"));
    prb_String  localMafftOuputPath = prb_pathJoin(arena, testsDir, prb_STR("localmafft-testseqs.fasta"));
    prb_String* localMafftAlignedSeqs = alignWithMafft(arena, localMafftExe, fastaOutputPath, localMafftOuputPath);
    prb_assert(arrlen(localMafftAlignedSeqs) == genSeq.seqCount);

    prb_writelnToStdout(prb_fmt(arena, "tests took %.2fms", prb_getMsFrom(testsStart)));
    return 0;
}
