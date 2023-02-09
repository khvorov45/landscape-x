#include "cbuild.h"

// clang-format off
#define aln_IMPLEMENTATION
#define aln_PUBLICAPI static
#define aln_assert(cond) prb_assert(cond)
#include "../src/align.h"
// clang-format on

#define function static

typedef intptr_t isize;

typedef struct MatrixStr {
    prb_Str str;
    isize   charPitch;
    isize   width;
    isize   height;
    isize   cellWidth;
    isize   cellHeight;
} MatrixStr;

function MatrixStr
createMatrixStr(prb_Arena* arena, isize width, isize height) {
    isize cellWidth = 3;
    isize cellHeight = 2;
    isize charCols = width * cellWidth;
    isize charRows = height * cellHeight;
    isize charPitch = charCols + 1;  // NOTE(sen) +1 for newlines
    isize bytes = charPitch * charRows;
    char* buf = prb_arenaFreePtr(arena);
    prb_arenaChangeUsed(arena, bytes);
    for (isize byteIndex = 0; byteIndex < bytes; byteIndex++) {
        buf[byteIndex] = ' ';
    }
    for (isize rowEndIndex = charCols; rowEndIndex < bytes; rowEndIndex += charCols + 1) {
        buf[rowEndIndex] = '\n';
    }
    MatrixStr result = {{buf, bytes}, charPitch, width, height, cellWidth, cellHeight};
    return result;
}

function isize
getTopleftBufInd(MatrixStr* matStr, isize row, isize col) {
    prb_assert(row < matStr->height);
    prb_assert(col < matStr->width);
    isize cellLeft = col * matStr->cellWidth;
    isize cellTop = row * matStr->cellHeight;
    isize topleftBufInd = cellTop * matStr->charPitch + cellLeft;
    return topleftBufInd;
}

function void
addCellBottomRight(prb_Arena* arena, MatrixStr* matStr, isize row, isize col, float number) {
    isize   topleftBufInd = getTopleftBufInd(matStr, row, col);
    isize   bottomrightBufInd = topleftBufInd + (matStr->cellWidth - 1) + ((matStr->cellHeight - 1) * matStr->charPitch);
    prb_Str numStr = prb_fmt(arena, "%.0f", number);
    isize   bottomrightOffset = 0;
    for (isize strind = prb_min(1, numStr.len - 1); strind >= 0; strind--) {
        ((char*)matStr->str.ptr)[bottomrightBufInd + bottomrightOffset] = numStr.ptr[strind];
        bottomrightOffset -= 1;
    }
}

function void
addCellTopRight(MatrixStr* matStr, isize row, isize col, char ch) {
    isize topleftBufInd = getTopleftBufInd(matStr, row, col);
    isize topcenterBufInd = topleftBufInd + matStr->cellWidth - 1;
    ((char*)matStr->str.ptr)[topcenterBufInd] = ch;
}

function void
addCellBottomLeft(MatrixStr* matStr, isize row, isize col, char ch) {
    isize topleftBufInd = getTopleftBufInd(matStr, row, col);
    isize bottomleftBufInd = topleftBufInd + (matStr->cellHeight - 1) * matStr->charPitch;
    ((char*)matStr->str.ptr)[bottomleftBufInd] = ch;
}

function void
addCellTopleft(MatrixStr* matStr, isize row, isize col, char ch) {
    isize topleftBufInd = getTopleftBufInd(matStr, row, col);
    ((char*)matStr->str.ptr)[topleftBufInd] = ch;
}

function void
printMatrix(prb_Arena* arena, aln_Matrix2NW mat, aln_Str reference, aln_Str ogstr) {
    prb_TempMemory temp = prb_beginTempMemory(arena);

    MatrixStr matStr = createMatrixStr(arena, reference.len + 1, ogstr.len + 1);
    addCellBottomRight(arena, &matStr, 0, 0, 0.0f);

    for (isize refInd = 0; refInd < reference.len; refInd++) {
        addCellTopRight(&matStr, 0, refInd + 1, reference.ptr[refInd]);
        addCellBottomRight(arena, &matStr, 0, refInd + 1, -(float)(refInd + 1));
    }

    for (isize ogInd = 0; ogInd < ogstr.len; ogInd++) {
        addCellBottomLeft(&matStr, ogInd + 1, 0, ogstr.ptr[ogInd]);
        addCellBottomRight(arena, &matStr, ogInd + 1, 0, -(float)(ogInd + 1));
    }

    for (isize rowIndex = 1; rowIndex < mat.height; rowIndex++) {
        for (isize colIndex = 1; colIndex < mat.width; colIndex++) {
            aln_NWEntry entry = aln_matrix2get(mat, rowIndex, colIndex);
            addCellBottomRight(arena, &matStr, rowIndex, colIndex, entry.score);
            switch (entry.cameFromDir) {
                case aln_CameFromDir_TopLeft: addCellTopleft(&matStr, rowIndex, colIndex, '\\'); break;
                case aln_CameFromDir_Top: addCellTopRight(&matStr, rowIndex, colIndex, '^'); break;
                case aln_CameFromDir_Left: addCellBottomLeft(&matStr, rowIndex, colIndex, '<'); break;
            }
        }
    }

    prb_writeToStdout(matStr.str);
    prb_endTempMemory(temp);
}

function void
streq(prb_Arena* arena, aln_Str str1, aln_Str str2) {
    bool result = false;
    if (str1.len == str2.len) {
        result = prb_memeq(str1.ptr, str2.ptr, str1.len);
    }
    if (!result) {
        prb_writeToStdout(prb_fmt(arena, "expected: %.*s\nactual: %.*s\n", prb_LIT(str1), prb_LIT(str2)));
        prb_assert(!"failed");
    }
}

function void
alignAndReconstruct(
    prb_Arena* arena,
    aln_Str*   refs,
    intptr_t   refCount,
    aln_Str*   strs,
    intptr_t   strsCount,
    aln_Str*   expectedRefs,
    aln_Str*   expectedStrs
) {
    prb_TempMemory temp = prb_beginTempMemory(arena);

    prb_Arena       alnOut = prb_createArenaFromArena(arena, 20 * prb_MEGABYTE);
    aln_AlignResult alignResult = aln_align(
        refs,
        refCount,
        strs,
        strsCount,
        (aln_Config) {
            .outmem = alnOut.base,
            .outmemBytes = alnOut.size,
            .tempmem = prb_arenaFreePtr(arena),
            .tempmemBytes = prb_arenaFreeSize(arena),
            .storeFinalMatrices = true,
        }
    );

    for (isize seqInd = 0; seqInd < strsCount; seqInd++) {
        aln_Str ogstr = strs[seqInd];
        aln_Str reference = refs[0];
        if (refCount > 1) {
            reference = refs[seqInd];
        }
        aln_Alignment alignedStr = alignResult.strs[seqInd];

        bool printMats = false;
        if (printMats) {
            printMatrix(arena, alignResult.matrices[seqInd], reference, ogstr);
        }

        aln_Str refReconstructed = aln_reconstruct(alignedStr, aln_Reconstruct_Ref, reference, ogstr, prb_arenaFreePtr(arena), prb_arenaFreeSize(arena));
        prb_arenaChangeUsed(arena, refReconstructed.len);
        streq(arena, refReconstructed, expectedRefs[seqInd]);

        aln_Str strReconstructed = aln_reconstruct(alignedStr, aln_Reconstruct_Str, reference, ogstr, prb_arenaFreePtr(arena), prb_arenaFreeSize(arena));
        prb_arenaChangeUsed(arena, strReconstructed.len);
        streq(arena, strReconstructed, expectedStrs[seqInd]);
    }

    prb_endTempMemory(temp);
}

function void
test_alignAndReconstruct(prb_Arena* arena) {
    {
        aln_Str reference = aln_STR("ABC");
        aln_Str seqs[] = {aln_STR("ABC"), aln_STR("BC"), aln_STR("AB"), aln_STR("B"), aln_STR("DABC"), aln_STR("ABCD"), aln_STR("DABCD")};

        aln_Str expectedRefs[] = {aln_STR("ABC"), aln_STR("ABC"), aln_STR("ABC"), aln_STR("ABC"), aln_STR("-ABC"), aln_STR("ABC-"), aln_STR("-ABC-")};
        aln_Str expectedSeqs[] = {aln_STR("ABC"), aln_STR("-BC"), aln_STR("AB-"), aln_STR("-B-"), aln_STR("DABC"), aln_STR("ABCD"), aln_STR("DABCD")};

        alignAndReconstruct(
            arena,
            &reference,
            1,
            seqs,
            prb_arrayCount(seqs),
            expectedRefs,
            expectedSeqs
        );
    }

    {
        aln_Str references[] = {aln_STR("ABC"), aln_STR("DEFABCTYU"), aln_STR("ABCDEFGH")};
        aln_Str seqs[] = {aln_STR("ABC"), aln_STR("FABSTY"), aln_STR("ABCFGH")};

        aln_Str expectedRefs[] = {aln_STR("ABC"), aln_STR("DEFABCTYU"), aln_STR("ABCDEFGH")};
        aln_Str expectedSeqs[] = {aln_STR("ABC"), aln_STR("--FABSTY-"), aln_STR("ABC--FGH")};

        alignAndReconstruct(
            arena,
            references,
            prb_arrayCount(references),
            seqs,
            prb_arrayCount(seqs),
            expectedRefs,
            expectedSeqs
        );
    }
}

int
main() {
    prb_Arena  arena_ = prb_createArenaFromVmem(1 * prb_GIGABYTE);
    prb_Arena* arena = &arena_;

    test_alignAndReconstruct(arena);

    return 0;
}
