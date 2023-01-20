#include "../cbuild.h"

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

#if 0
function prb_Str*
getSeqsFromFile(prb_Arena* arena, prb_Str filepath) {
    prb_Str*                 seqs = 0;
    prb_ReadEntireFileResult readRes = prb_readEntireFile(arena, filepath);
    prb_assert(readRes.success);
    prb_Str        contentLeft = prb_strFromBytes(readRes.content);
    prb_StrScanner lineIter = prb_createStrScanner(contentLeft);
    for (;;) {
        prb_StrFindSpec linebr = {.mode = prb_StrFindMode_LineBreak};
        if (prb_strScannerMove(&lineIter, linebr, prb_StrScannerSide_AfterMatch) == prb_Failure) {
            break;
        }
        prb_assert(lineIter.betweenLastMatches.ptr[0] == '>');
        prb_GrowingStr gstr = prb_beginStr(arena);
        for (;;) {
            prb_StrScanner lineIterCopy = lineIter;
            if (prb_strScannerMove(&lineIterCopy, linebr, prb_StrScannerSide_AfterMatch) == prb_Failure) {
                break;
            }
            if (lineIterCopy.betweenLastMatches.ptr[0] == '>') {
                break;
            }
            lineIter = lineIterCopy;
            prb_addStrSegment(&gstr, "%.*s", prb_LIT(lineIter.betweenLastMatches));
        }
        prb_Str seq = prb_endStr(&gstr);
        arrput(seqs, seq);
    }
    return seqs;
}

function prb_Str*
alignWithMafft(prb_Arena* arena, prb_Str mafftExe, prb_Str inputPath, prb_Str mafftOuputPath) {
    {
        prb_Str cmd = prb_fmt(arena, "%.*s --globalpair --maxiterate 1000 %.*s", prb_LIT(mafftExe), prb_LIT(inputPath));
        prb_writelnToStdout(arena, cmd);
        prb_Process proc = prb_createProcess(cmd, (prb_ProcessSpec) {.redirectStdout = true, .stdoutFilepath = mafftOuputPath});
        prb_assert(prb_launchProcesses(arena, &proc, 1, prb_Background_No));
    }
    prb_Str* mafftAlignedSeqs = getSeqsFromFile(arena, mafftOuputPath);
    return mafftAlignedSeqs;
}

aln_AlignResult tbfast_main(aln_Str* strings, intptr_t stringsCount, void* out, intptr_t outBytes);
#endif 

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

#if 0
    prb_Str       testsDir = prb_getParentDir(arena, prb_STR(__FILE__));
    prb_Str       rootDir = prb_getParentDir(arena, testsDir);

    prb_Str fastaOutputPath = prb_pathJoin(arena, testsDir, prb_STR("testseqs.fasta"));
    {
        prb_GrowingStr gstr = prb_beginStr(arena);
        for (i32 seqIndex = 0; seqIndex < genSeq.seqCount; seqIndex++) {
            prb_addStrSegment(&gstr, ">virus%d\n%.*s\n", seqIndex, prb_LIT(genSeq.seqs[seqIndex]));
        }
        prb_Str fastaContent = prb_endStr(&gstr);
        prb_writeEntireFile(arena, fastaOutputPath, fastaContent.ptr, fastaContent.len);
    }

    prb_assert(prb_removePathIfExists(arena, prb_pathJoin(arena, rootDir, prb_STR("mafft/core/logfile.txt"))) == prb_Success);

    prb_Str          mafftBinEnvName = prb_STR("MAFFT_BINARIES");
    prb_GetenvResult oldMafftBin = prb_getenv(arena, mafftBinEnvName);
    prb_assert(oldMafftBin.found);
    prb_assert(prb_unsetenv(arena, mafftBinEnvName));
    prb_Str  mafftOuputPath = prb_pathJoin(arena, testsDir, prb_STR("mafft-testseqs.fasta"));
    prb_Str* mafftAlignedSeqs = alignWithMafft(arena, prb_STR("mafft"), fastaOutputPath, mafftOuputPath);
    prb_assert(arrlen(mafftAlignedSeqs) == genSeq.seqCount);
    prb_assert(prb_setenv(arena, mafftBinEnvName, oldMafftBin.str));

    prb_Str  localMafftExe = prb_pathJoin(arena, rootDir, prb_STR("mafft/core/mafft.tmpl"));
    prb_Str  localMafftOuputPath = prb_pathJoin(arena, testsDir, prb_STR("localmafft-testseqs.fasta"));
    prb_Str* localMafftAlignedSeqs = alignWithMafft(arena, localMafftExe, fastaOutputPath, localMafftOuputPath);
    prb_assert(arrlen(localMafftAlignedSeqs) == genSeq.seqCount);

    for (i32 seqIndex = 0; seqIndex < genSeq.seqCount; seqIndex++) {
        prb_assert(prb_streq(mafftAlignedSeqs[seqIndex], localMafftAlignedSeqs[seqIndex]));
    }

    // NOTE(sen) Call what I pulled out directly
    prb_Str tempDir = prb_pathJoin(arena, testsDir, prb_STR("tmp"));
    prb_assert(prb_clearDir(arena, tempDir));
    prb_Str cwd = prb_getWorkingDir(arena);
    prb_assert(prb_setWorkingDir(arena, tempDir));
    intptr_t        outBytes = prb_arenaFreeSize(arena);
    void*           outBuf = prb_arenaFreePtr(arena);
    aln_AlignResult myAlignResult = tbfast_main((aln_Str*)genSeq.seqs, genSeq.seqCount, outBuf, outBytes);
    prb_arenaChangeUsed(arena, myAlignResult.bytesWritten);
    prb_assert(prb_setWorkingDir(arena, cwd));

    prb_assert(myAlignResult.seqCount == genSeq.seqCount);
    for (i32 seqIndex = 0; seqIndex < genSeq.seqCount; seqIndex++) {
        aln_Str myStr = myAlignResult.seqs[seqIndex];
        prb_Str myStrPrb = {myStr.ptr, (int32_t)myStr.len};
        prb_writelnToStdout(arena, myStrPrb);
        // prb_writelnToStdout(arena, mafftAlignedSeqs[seqIndex]);
        // prb_assert(prb_streq(mafftAlignedSeqs[seqIndex], myStrPrb));
    }
#endif

    aln_AlignResult alnResult = aln_align((aln_Str*)genSeq.seqs, genSeq.seqCount, prb_arenaFreePtr(arena), 20 * prb_MEGABYTE);
    prb_arenaChangeUsed(arena, alnResult.bytesWritten);

    prb_writelnToStdout(arena, prb_fmt(arena, "tests took %.2fms", prb_getMsFrom(testsStart)));
    return 0;
}
