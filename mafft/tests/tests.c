#include "../programmable_build.h"

#define mafft_IMPLEMENTATION
#define mafft_assert(x) prb_assert(x)
#include "../rewrite/mafft.h"

#define function static

typedef int32_t i32;
typedef uint8_t u8;

function mafft_Arena
mafftArenaFromPrbArena(prb_Arena* arena) {
    prb_Arena   mafftMem = prb_createArenaFromArena(arena, 1 * prb_MEGABYTE);
    mafft_Arena mafftArena = {};
    mafft_initArena(&mafftArena, prb_arenaFreePtr(&mafftMem), prb_arenaFreeSize(&mafftMem));
    return mafftArena;
}

function void
test_readFasta(prb_Arena* arena) {
    prb_TempMemory        temp = prb_beginTempMemory(arena);
    prb_String            fastaStr = prb_STR(">Virus1\nQDFPGND\r\n>Virus2\rQDLPGNDNSTAT\n>Virus3\nQDLPGN\nDNSTAT");
    mafft_Arena           mafftArena = mafftArenaFromPrbArena(arena);
    mafft_ReadFastaResult fastaResult = mafft_readFasta(&mafftArena, fastaStr.ptr, fastaStr.len);
    prb_assert(fastaResult.success);
    prb_assert(fastaResult.fasta.entryCount == 3);
    prb_assert(prb_streq(prb_STR(fastaResult.fasta.entries[0].name.ptr), prb_STR("Virus1")));
    prb_assert(prb_streq(prb_STR(fastaResult.fasta.entries[0].seq.ptr), prb_STR("QDFPGND")));
    prb_assert(prb_streq(prb_STR(fastaResult.fasta.entries[1].name.ptr), prb_STR("Virus2")));
    prb_assert(prb_streq(prb_STR(fastaResult.fasta.entries[1].seq.ptr), prb_STR("QDLPGNDNSTAT")));
    prb_assert(prb_streq(prb_STR(fastaResult.fasta.entries[2].name.ptr), prb_STR("Virus3")));
    prb_assert(prb_streq(prb_STR(fastaResult.fasta.entries[2].seq.ptr), prb_STR("QDLPGNDNSTAT")));
    prb_assert(mafftArena.tempCount == 0);
    prb_assert(mafftArena.lockedForString == 0);
    prb_endTempMemory(temp);
}

function void
test_alignSeq(prb_Arena* arena) {
    prb_TempMemory temp = prb_beginTempMemory(arena);

    mafft_String seqs[] = {
        mafft_STR("QDFPGNDNSTATLCLGHHAVPNGTLVKTITNDQIEVTNATELVQSSSTGKICNNPHRILD"),
        mafft_STR("QDFPGNDNSTATLCLGHHAVPNGTLVKTITNDQIEVTNATELVQSSSTGKICNNPHRIL"),
        mafft_STR("DFPGNDNSTATLCLGHHAVPNGTLVKTITNDQIEVTNATELVQSSSTGKICNNPHRILD"),
    };

    i32   alignedOutputBufSize = mafft_getAlignOutputBufferSize(seqs, prb_arrayLength(seqs));
    void* alignedOutputBuf = prb_arenaAllocAndZero(arena, alignedOutputBufSize, 1);

    mafft_alignSeq(seqs, prb_arrayLength(seqs), alignedOutputBuf);
    prb_String seqsAligned[] = {
        prb_STR("QDFPGNDNSTATLCLGHHAVPNGTLVKTITNDQIEVTNATELVQSSSTGKICNNPHRILD"),
        prb_STR("QDFPGNDNSTATLCLGHHAVPNGTLVKTITNDQIEVTNATELVQSSSTGKICNNPHRIL-"),
        prb_STR("-DFPGNDNSTATLCLGHHAVPNGTLVKTITNDQIEVTNATELVQSSSTGKICNNPHRILD"),
    };

    u8* outBufLeft = (u8*)alignedOutputBuf;
    for (i32 seqIndex = 0; seqIndex < prb_arrayLength(seqs); seqIndex++) {
        prb_String seqAlignedMafft = prb_STR((const char*)outBufLeft);
        prb_String seqAlignedRef = seqsAligned[seqIndex];
        prb_assert(prb_streq(seqAlignedMafft, seqAlignedRef));
    }

    prb_endTempMemory(temp);
}

int
main() {
    prb_TimeStart testsStart = prb_timeStart();
    prb_Arena     arena_ = prb_createArenaFromVmem(1 * prb_GIGABYTE);
    prb_Arena*    arena = &arena_;

    test_readFasta(arena);
    test_alignSeq(arena);

    prb_writelnToStdout(prb_fmt(arena, "tests took %.2fms", prb_getMsFrom(testsStart)));
    return 0;
}
