#include "cbuild.h"

// clang-format off
#define aln_IMPLEMENTATION
#define aln_PUBLICAPI static
#define aln_assert(cond) prb_assert(cond)
#include "../src/align.h"
// clang-format on

typedef intptr_t isize;

int
main() {
    prb_Arena  arena_ = prb_createArenaFromVmem(1 * prb_GIGABYTE);
    prb_Arena* arena = &arena_;

    aln_Str reference = aln_STR("GTCCG");
    aln_Str seqs[] = {aln_STR("TCC"), aln_STR("GTCC")};

    prb_Arena       alnOut = prb_createArenaFromArena(arena, 20 * prb_MEGABYTE);
    aln_AlignResult alignResult = aln_align(
        reference,
        seqs,
        prb_arrayCount(seqs),
        (aln_Config) {
            .outmem = alnOut.base,
            .outmemBytes = alnOut.size,
            .tempmem = prb_arenaFreePtr(arena),
            .tempmemBytes = prb_arenaFreeSize(arena),
            .storeFinalMatrices = true,
        }
    );

    for (isize seqInd = 0; seqInd < prb_arrayCount(seqs); seqInd++) {
        aln_Str       ogstr = seqs[seqInd];
        aln_Matrix2NW mat = alignResult.matrices[seqInd];

        prb_GrowingStr matStrBuilder = prb_beginStr(arena);

        // NOTE(sen) First row with referece seq
        prb_addStrSegment(&matStrBuilder, "     ");
        for (isize colIndex = 1; colIndex < mat.width; colIndex++) {
            prb_addStrSegment(&matStrBuilder, " %c ", reference.ptr[colIndex - 1]);
        }
        prb_addStrSegment(&matStrBuilder, "\n");

        for (isize rowIndex = 0; rowIndex < mat.height; rowIndex++) {
            if (rowIndex > 0) {
                prb_addStrSegment(&matStrBuilder, "    ");
                for (isize colIndex = 1; colIndex < mat.width; colIndex++) {
                    aln_CameFromDir cameFrom = aln_matrix2get(mat, rowIndex, colIndex).cameFromDir;
                    switch (cameFrom) {
                        case aln_CameFromDir_TopLeft: prb_addStrSegment(&matStrBuilder, "↖  "); break;
                        case aln_CameFromDir_Top: prb_addStrSegment(&matStrBuilder, "  ↑"); break;
                        case aln_CameFromDir_Left: prb_addStrSegment(&matStrBuilder, "   "); break;
                    }
                }
                prb_addStrSegment(&matStrBuilder, "\n%c", ogstr.ptr[rowIndex - 1]);
            } else {
                prb_addStrSegment(&matStrBuilder, " ");
            }

            for (isize colIndex = 0; colIndex < mat.width; colIndex++) {
                aln_NWEntry     entry = aln_matrix2get(mat, rowIndex, colIndex);
                aln_CameFromDir cameFrom = entry.cameFromDir;
                if (cameFrom == aln_CameFromDir_Left && colIndex > 0 && rowIndex > 0) {
                    prb_addStrSegment(&matStrBuilder, "←");
                } else {
                    prb_addStrSegment(&matStrBuilder, " ");
                }

                if (entry.score >= 0) {
                    prb_addStrSegment(&matStrBuilder, " ");
                }
                prb_addStrSegment(&matStrBuilder, "%.0f", entry.score);
            }

            prb_addStrSegment(&matStrBuilder, "\n");
        }  // for row
        prb_addStrSegment(&matStrBuilder, "\n");

        prb_Str matStr = prb_endStr(&matStrBuilder);
        prb_writeToStdout(matStr);

        aln_Str alignedStr = alignResult.strs[seqInd];
        prb_writelnToStdout(arena, (prb_Str) {alignedStr.ptr, alignedStr.len});
    }

    return 0;
}