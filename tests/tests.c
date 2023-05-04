#include "cbuild.h"

// clang-format off
#define aln_IMPLEMENTATION
#define aln_PUBLICAPI static
#define aln_assert(cond) prb_assert(cond)
#include "../src/align.h"

#define aln_STR(x) (aln_Str) {x, sizeof(x) - 1}
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
        prb_writeToStdout(prb_fmt(arena, "expected: `%.*s`\nactual: `%.*s`\n", prb_LIT(str1), prb_LIT(str2)));
        prb_assert(!"failed");
    }
}

function void
alignAndReconstruct(
    prb_Arena*   arena,
    aln_StrArray refs,
    aln_StrArray strs,
    aln_StrArray expectedRefs,
    aln_StrArray expectedStrs
) {
    prb_assert(refs.len == 1 || refs.len == strs.len);
    prb_assert(strs.len == expectedRefs.len);
    prb_assert(expectedRefs.len == expectedStrs.len);

    aln_Memory alnMem = aln_createMemory(prb_arenaFreePtr(arena), prb_arenaFreeSize(arena), 20 * prb_MEGABYTE);
    aln_AlignResult alignResult = aln_align(refs, strs, (aln_Config) {.storeFinalMatrices = true}, &alnMem);

    prb_assert(alignResult.alignments.len == strs.len);
    prb_assert(alignResult.matrices.len == strs.len);

    for (isize seqInd = 0; seqInd < strs.len; seqInd++) {
        aln_Str ogstr = strs.ptr[seqInd];
        aln_Str reference = refs.ptr[0];
        if (refs.len > 1) {
            reference = refs.ptr[seqInd];
        }
        aln_Alignment alignedStr = alignResult.alignments.ptr[seqInd];

        bool printMats = false;
        if (printMats) {
            printMatrix(arena, alignResult.matrices.ptr[seqInd], reference, ogstr);
        }

        aln_Str refReconstructed = aln_reconstruct(alignedStr, aln_Reconstruct_Ref, reference, ogstr, &alnMem.perm);
        streq(arena, expectedRefs.ptr[seqInd], refReconstructed);

        aln_Str strReconstructed = aln_reconstruct(alignedStr, aln_Reconstruct_Str, reference, ogstr, &alnMem.perm);
        streq(arena, expectedStrs.ptr[seqInd], strReconstructed);
    }
}

function void
alignAndReconstructToCommon(
    prb_Arena*   arena,
    aln_Str      ref,
    aln_StrArray strs,
    aln_Str      expectedRef,
    aln_StrArray expectedStrs
) {
    prb_assert(strs.len == expectedStrs.len);

    aln_Memory alnMem = aln_createMemory(prb_arenaFreePtr(arena), prb_arenaFreeSize(arena), 20 * prb_MEGABYTE);
    aln_AlignResult alignResult = aln_align((aln_StrArray) {&ref, 1}, strs, (aln_Config) {.storeFinalMatrices = true}, &alnMem);

    aln_ReconstructToCommonRefResult reconstruction = aln_reconstructToCommonRef(alignResult.alignments, ref, strs, &alnMem);
    streq(arena, expectedRef, reconstruction.commonRef);
    prb_assert(reconstruction.alignedStrs.len == strs.len);
    for (isize strInd = 0; strInd < strs.len; strInd++) {
        aln_Str str = reconstruction.alignedStrs.ptr[strInd];
        aln_Str strExpected = expectedStrs.ptr[strInd];

        bool printMats = false;
        if (printMats) {
            printMatrix(arena, alignResult.matrices.ptr[strInd], ref, strs.ptr[strInd]);
        }

        streq(arena, strExpected, str);
    }
}

function void
test_alignAndReconstruct(prb_Arena* arena) {
    {
        aln_Str reference = aln_STR("ABC");
        aln_Str seqs[] = {aln_STR("ABC"), aln_STR("BC"), aln_STR("AB"), aln_STR("B"), aln_STR("DDABC"), aln_STR("ABCD"), aln_STR("DABCD"), aln_STR("AB12C")};

        aln_Str expectedRefs[] = {aln_STR("ABC"), aln_STR("ABC"), aln_STR("ABC"), aln_STR("ABC"), aln_STR("--ABC"), aln_STR("ABC-"), aln_STR("-ABC-"), aln_STR("AB--C")};
        aln_Str expectedSeqs[] = {aln_STR("ABC"), aln_STR("-BC"), aln_STR("AB-"), aln_STR("-B-"), aln_STR("DDABC"), aln_STR("ABCD"), aln_STR("DABCD"), aln_STR("AB12C")};

        alignAndReconstruct(
            arena,
            (aln_StrArray) {&reference, 1},
            (aln_StrArray) {seqs, prb_arrayCount(seqs)},
            (aln_StrArray) {expectedRefs, prb_arrayCount(expectedRefs)},
            (aln_StrArray) {expectedSeqs, prb_arrayCount(expectedSeqs)}
        );

        aln_Str expectedCommonRef = aln_STR("--AB--C-");
        aln_Str expectedCommonSeqs[] = {aln_STR("--AB--C-"), aln_STR("---B--C-"), aln_STR("--AB----"), aln_STR("---B----"), aln_STR("DDAB--C-"), aln_STR("--AB--CD"), aln_STR("-DAB--CD"), aln_STR("--AB12C-")};

        alignAndReconstructToCommon(
            arena,
            reference,
            (aln_StrArray) {seqs, prb_arrayCount(seqs)},
            expectedCommonRef,
            (aln_StrArray) {expectedCommonSeqs, prb_arrayCount(expectedCommonSeqs)}
        );
    }

    {
        aln_Str references[] = {aln_STR("ABC"), aln_STR("DEFABCTYU"), aln_STR("ABCDEFGH")};
        aln_Str seqs[] = {aln_STR("ABC"), aln_STR("FABSTY"), aln_STR("ABCFGH")};

        aln_Str expectedRefs[] = {aln_STR("ABC"), aln_STR("DEFABCTYU"), aln_STR("ABCDEFGH")};
        aln_Str expectedSeqs[] = {aln_STR("ABC"), aln_STR("--FABSTY-"), aln_STR("ABC--FGH")};

        alignAndReconstruct(
            arena,
            (aln_StrArray) {references, prb_arrayCount(references)},
            (aln_StrArray) {seqs, prb_arrayCount(seqs)},
            (aln_StrArray) {expectedRefs, prb_arrayCount(expectedRefs)},
            (aln_StrArray) {expectedSeqs, prb_arrayCount(expectedSeqs)}
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
