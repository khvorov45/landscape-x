#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
#endif

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#endif

#ifndef aln_HEADER
#define aln_HEADER

#include <stdint.h>
#include <stdbool.h>
#include <stdalign.h>

#ifndef aln_PUBLICAPI
#define aln_PUBLICAPI
#endif

#define aln_matrix2get(matrix, row, col) matrix.ptr[aln_matrix2index(matrix.width, matrix.height, row, col)]

typedef enum aln_Status {
    aln_Failure,
    aln_Success,
} aln_Status;

typedef struct aln_Arena {
    void*    base;
    intptr_t size;
    intptr_t used;
    intptr_t tempCount;
    bool     lockedForStr;
} aln_Arena;

typedef struct aln_Memory {
    aln_Arena perm;
    aln_Arena temp;
} aln_Memory;

typedef struct aln_Str {
    char*    ptr;
    intptr_t len;
} aln_Str;

typedef struct aln_StrArray {
    aln_Str* ptr;
    intptr_t len;
} aln_StrArray;

typedef struct aln_Config {
    bool storeFinalMatrices;
} aln_Config;

typedef enum aln_CameFromDir {
    aln_CameFromDir_TopLeft,
    aln_CameFromDir_Top,
    aln_CameFromDir_Left,
} aln_CameFromDir;

// Entry in a Needleman–Wunsch alignment matrix
typedef struct aln_NWEntry {
    float           score;
    aln_CameFromDir cameFromDir;
} aln_NWEntry;

typedef struct aln_Matrix2NW {
    aln_NWEntry* ptr;
    intptr_t     width;
    intptr_t     height;
} aln_Matrix2NW;

typedef struct aln_Matrix2NWArray {
    aln_Matrix2NW* ptr;
    intptr_t       len;
} aln_Matrix2NWArray;

typedef enum aln_AlignAction {
    aln_AlignAction_Match,
    aln_AlignAction_GapStr,
    aln_AlignAction_GapRef,
} aln_AlignAction;

typedef struct aln_Alignment {
    aln_AlignAction* actions;
    intptr_t         actionCount;
} aln_Alignment;

typedef struct aln_AlignmentArray {
    aln_Alignment* ptr;
    intptr_t       len;
} aln_AlignmentArray;

typedef struct aln_AlignResult {
    aln_AlignmentArray alignments;
    aln_Matrix2NWArray matrices;
} aln_AlignResult;

typedef enum aln_Reconstruct {
    aln_Reconstruct_Ref,
    aln_Reconstruct_Str,
} aln_Reconstruct;

typedef struct aln_ReconstructToCommonRefResult {
    aln_Str      commonRef;
    aln_StrArray alignedStrs;
} aln_ReconstructToCommonRefResult;

aln_PUBLICAPI aln_Memory                       aln_createMemory(void* base, intptr_t totalSize, intptr_t permSize);
aln_PUBLICAPI aln_Memory                       aln_createMemory2(void* permBase, intptr_t permSize, void* tempBase, intptr_t tempSize);
aln_PUBLICAPI aln_AlignResult                  aln_align(aln_StrArray references, aln_StrArray strings, aln_Config config, aln_Memory* memory);
aln_PUBLICAPI aln_Str                          aln_reconstruct(aln_Alignment aligned, aln_Reconstruct which, aln_Str reference, aln_Str ogstr, aln_Arena* arena);
aln_PUBLICAPI aln_ReconstructToCommonRefResult aln_reconstructToCommonRef(aln_AlignmentArray alignments, aln_Str reference, aln_StrArray strings, aln_Memory* memory);
aln_PUBLICAPI intptr_t                         aln_matrix2index(intptr_t matrixWidth, intptr_t matrixHeight, intptr_t row, intptr_t col);

#endif  // aln_HEADER

#ifdef aln_IMPLEMENTATION

// clang-format off
#define aln_max(a, b) (((a) > (b)) ? (a) : (b))
#define aln_min(a, b) (((a) < (b)) ? (a) : (b))
#define aln_abs(a) (((a) < (0)) ? (-(a)) : (a))
#define aln_clamp(x, a, b) (((x) < (a)) ? (a) : (((x) > (b)) ? (b) : (x)))
#define aln_arrayCount(arr) (intptr_t)(sizeof(arr) / sizeof(arr[0]))
#define aln_arenaAllocArray(arena, type, len) (type*)aln_arenaAllocAndZero(arena, (len) * (intptr_t)sizeof(type), alignof(type))
#define aln_arenaAllocStruct(arena, type) (type*)aln_arenaAllocAndZero(arena, sizeof(type), alignof(type))
#define aln_arenaAllocMatrix2(matrixType, dataType, arena, width_, height_) (matrixType) {.ptr = aln_arenaAllocArray(arena, dataType, width_ * height_), .width = width_, .height = height_}
#define aln_isPowerOf2(x) (((x) > 0) && (((x) & ((x)-1)) == 0))
#define aln_unused(x) ((x) = (x))
#define aln_STR(x) (aln_Str) {x, sizeof(x) - 1}

#ifndef aln_assert
#define aln_assert(condition)
#endif
// clang-format on

#ifndef aln_PRIVATEAPI
#define aln_PRIVATEAPI static
#endif

aln_PRIVATEAPI bool
aln_memeq(void* ptr1, void* ptr2, intptr_t len) {
    bool result = true;
    for (intptr_t ind = 0; ind < len; ind++) {
        uint8_t b1 = ((uint8_t*)ptr1)[ind];
        uint8_t b2 = ((uint8_t*)ptr2)[ind];
        if (b1 != b2) {
            result = false;
            break;
        }
    }
    return result;
}

aln_PRIVATEAPI bool
aln_streq(aln_Str str1, aln_Str str2) {
    bool result = false;
    if (str1.len == str2.len) {
        result = aln_memeq(str1.ptr, str2.ptr, str1.len);
    }
    return result;
}

aln_PRIVATEAPI intptr_t
aln_getOffsetForAlignment(void* ptr, intptr_t align) {
    aln_assert(aln_isPowerOf2(align));
    uintptr_t ptrAligned = (uintptr_t)((uint8_t*)ptr + (align - 1)) & (uintptr_t)(~(align - 1));
    aln_assert(ptrAligned >= (uintptr_t)ptr);
    intptr_t diff = (intptr_t)(ptrAligned - (uintptr_t)ptr);
    aln_assert(diff < align && diff >= 0);
    return (intptr_t)diff;
}

typedef struct aln_TempMemory {
    aln_Arena* arena;
    intptr_t   usedAtBegin;
    intptr_t   tempCountAtBegin;
} aln_TempMemory;

aln_PRIVATEAPI void*
aln_arenaFreePtr(aln_Arena* arena) {
    void* result = (uint8_t*)arena->base + arena->used;
    return result;
}

aln_PRIVATEAPI intptr_t
aln_arenaFreeSize(aln_Arena* arena) {
    intptr_t result = arena->size - arena->used;
    return result;
}

aln_PRIVATEAPI void
aln_arenaChangeUsed(aln_Arena* arena, intptr_t byteDelta) {
    aln_assert(aln_arenaFreeSize(arena) >= byteDelta);
    arena->used += byteDelta;
}

aln_PRIVATEAPI void
aln_arenaAlignFreePtr(aln_Arena* arena, intptr_t align) {
    intptr_t offset = aln_getOffsetForAlignment(aln_arenaFreePtr(arena), align);
    aln_arenaChangeUsed(arena, offset);
}

aln_PRIVATEAPI aln_Arena
aln_createArenaFromArena(aln_Arena* parent, intptr_t bytes) {
    aln_Arena arena = {.base = aln_arenaFreePtr(parent), .size = bytes};
    aln_arenaChangeUsed(parent, bytes);
    return arena;
}

aln_PRIVATEAPI void*
aln_arenaAllocAndZero(aln_Arena* arena, intptr_t size, intptr_t align) {
    aln_arenaAlignFreePtr(arena, align);
    void* result = aln_arenaFreePtr(arena);
    aln_arenaChangeUsed(arena, size);
    for (intptr_t ind = 0; ind < size; ind++) {
        ((uint8_t*)result)[ind] = 0;
    }
    return result;
}

aln_PRIVATEAPI aln_TempMemory
aln_beginTempMemory(aln_Arena* arena) {
    aln_TempMemory temp = {.arena = arena, .usedAtBegin = arena->used, .tempCountAtBegin = arena->tempCount};
    arena->tempCount += 1;
    return temp;
}

aln_PRIVATEAPI void
aln_endTempMemory(aln_TempMemory temp) {
    aln_assert(temp.arena->tempCount == temp.tempCountAtBegin + 1);
    temp.arena->used = temp.usedAtBegin;
    temp.arena->tempCount -= 1;
}

aln_PRIVATEAPI char*
aln_strGetNullTerminated(aln_Arena* arena, aln_Str str) {
    char* buf = aln_arenaAllocArray(arena, char, str.len + 1);
    for (intptr_t ind = 0; ind < str.len; ind++) {
        buf[ind] = str.ptr[ind];
    }
    return buf;
}

aln_PRIVATEAPI aln_Arena
aln_createArena(void* base, intptr_t size) {
    aln_Arena arena = {.base = base, .size = size};
    return arena;
}

aln_PUBLICAPI aln_Memory
aln_createMemory(void* base, intptr_t totalSize, intptr_t permSize) {
    aln_Arena  total = aln_createArena(base, totalSize);
    aln_Arena  perm = aln_createArenaFromArena(&total, permSize);
    aln_Arena  temp = aln_createArena(aln_arenaFreePtr(&total), aln_arenaFreeSize(&total));
    aln_Memory memory = {perm, temp};
    return memory;
}

aln_PUBLICAPI aln_Memory
aln_createMemory2(void* permBase, intptr_t permSize, void* tempBase, intptr_t tempSize) {
    aln_Arena  perm = aln_createArena(permBase, permSize);
    aln_Arena  temp = aln_createArena(tempBase, tempSize);
    aln_Memory memory = {perm, temp};
    return memory;
}

aln_PUBLICAPI aln_AlignResult
aln_align(aln_StrArray references, aln_StrArray strings, aln_Config config, aln_Memory* memory) {
    // NOTE(sen) Input validation
    {
        aln_assert(strings.ptr);
        aln_assert(strings.len > 0);
        aln_assert(references.ptr);
        aln_assert(references.len == 1 || references.len == strings.len);
        for (intptr_t strInd = 0; strInd < strings.len; strInd++) {
            aln_Str ogstr = strings.ptr[strInd];
            aln_assert(ogstr.ptr && ogstr.len > 0);
        }
        for (intptr_t refInd = 0; refInd < references.len; refInd++) {
            aln_Str ref = references.ptr[refInd];
            aln_assert(ref.ptr && ref.len > 0);
        }
    }

    aln_TempMemory temp = aln_beginTempMemory(&memory->temp);

    aln_Matrix2NW maxGrid = {.ptr = 0, .width = 0, .height = 0};
    {
        intptr_t longestInputLen = 0;
        for (intptr_t strInd = 0; strInd < strings.len; strInd++) {
            longestInputLen = aln_max(longestInputLen, strings.ptr[strInd].len);
        }

        intptr_t longestRefLen = 0;
        for (intptr_t refInd = 0; refInd < references.len; refInd++) {
            longestRefLen = aln_max(longestRefLen, references.ptr[refInd].len);
        }

        maxGrid = aln_arenaAllocMatrix2(aln_Matrix2NW, aln_NWEntry, &memory->temp, longestRefLen + 1, longestInputLen + 1);
    }

    aln_Matrix2NW* storedMatrices = 0;
    if (config.storeFinalMatrices) {
        storedMatrices = aln_arenaAllocArray(&memory->perm, aln_Matrix2NW, strings.len);
    }

    aln_Alignment* alignedSeqs = aln_arenaAllocArray(&memory->perm, aln_Alignment, strings.len);
    for (intptr_t strInd = 0; strInd < strings.len; strInd++) {
        aln_Str ogstr = strings.ptr[strInd];
        aln_Str thisRef = references.ptr[0];
        if (references.len > 1) {
            thisRef = references.ptr[strInd];
        }

        aln_Matrix2NW thisGrid = maxGrid;
        thisGrid.width = thisRef.len + 1;
        thisGrid.height = ogstr.len + 1;
        aln_assert(thisGrid.width <= maxGrid.width);
        aln_assert(thisGrid.height <= maxGrid.height);

        aln_matrix2get(thisGrid, 0, 0).score = 0;
        float edgeIndelScore = -0.5f;
        for (intptr_t colIndex = 1; colIndex < thisGrid.width; colIndex++) {
            aln_matrix2get(thisGrid, 0, colIndex) = (aln_NWEntry) {(float)(colIndex)*edgeIndelScore, aln_CameFromDir_Left};
        }
        for (intptr_t rowIndex = 1; rowIndex < thisGrid.height; rowIndex++) {
            aln_matrix2get(thisGrid, rowIndex, 0) = (aln_NWEntry) {(float)(rowIndex)*edgeIndelScore, aln_CameFromDir_Top};
        }

// TODO(khvorov) Diagonal iteration (start of)
#if 0
        for (intptr_t diagonalIndex = 2; diagonalIndex < thisGrid.width + thisGrid.height - 1; diagonalIndex++) {
            intptr_t entriesInDiagonal = 0;
            if (diagonalIndex < thisGrid.height) {
                entriesInDiagonal = aln_min(diagonalIndex + 1, thisGrid.width);
            } else {
                entriesInDiagonal = aln_min(thisGrid.width - (diagonalIndex - thisGrid.height + 1), thisGrid.height);
            }
            aln_unused(entriesInDiagonal);
        }
#endif

        // TODO(khvorov) This would need to go by diagonal rather than row/column
        for (intptr_t rowIndex = 1; rowIndex < thisGrid.height; rowIndex++) {
            for (intptr_t colIndex = 1; colIndex < thisGrid.width; colIndex++) {
                intptr_t topleftIndex = aln_matrix2index(thisGrid.width, thisGrid.height, rowIndex - 1, colIndex - 1);
                intptr_t topIndex = aln_matrix2index(thisGrid.width, thisGrid.height, rowIndex - 1, colIndex);
                intptr_t leftIndex = aln_matrix2index(thisGrid.width, thisGrid.height, rowIndex, colIndex - 1);

                float alignScore = 1.0f;
                {
                    char chRef = thisRef.ptr[colIndex - 1];
                    char chStr = ogstr.ptr[rowIndex - 1];
                    if (chRef != chStr) {
                        alignScore = -1.0f;
                    }
                }

                float indelScore = -2.0f;
                if (rowIndex == thisGrid.height - 1 || colIndex == thisGrid.width - 1) {
                    indelScore = edgeIndelScore;
                }

                float scoreTopleft = thisGrid.ptr[topleftIndex].score + alignScore;
                float scoreTop = thisGrid.ptr[topIndex].score + indelScore;
                float scoreLeft = thisGrid.ptr[leftIndex].score + indelScore;

                float           maxScore = scoreTopleft;
                aln_CameFromDir cameFrom = aln_CameFromDir_TopLeft;
                if (maxScore < scoreTop) {
                    maxScore = scoreTop;
                    cameFrom = aln_CameFromDir_Top;
                }
                if (maxScore < scoreLeft) {
                    maxScore = scoreLeft;
                    cameFrom = aln_CameFromDir_Left;
                }

                aln_matrix2get(thisGrid, rowIndex, colIndex) = (aln_NWEntry) {maxScore, cameFrom};
            }  // for mat col
        }  // for mat row

        if (config.storeFinalMatrices) {
            aln_Matrix2NW matCopy = aln_arenaAllocMatrix2(aln_Matrix2NW, aln_NWEntry, &memory->perm, thisGrid.width, thisGrid.height);
            for (intptr_t matIndex = 0; matIndex < matCopy.width * matCopy.height; matIndex++) {
                matCopy.ptr[matIndex] = thisGrid.ptr[matIndex];
            }
            storedMatrices[strInd] = matCopy;
        }

        // TODO(khvorov) This reconstruction would have to go by diagonal as well
        intptr_t      maxActionCount = (thisGrid.width - 1) + (thisGrid.height - 1);
        aln_Alignment alignedStr = {aln_arenaAllocArray(&memory->perm, aln_AlignAction, maxActionCount)};

        intptr_t matInd = thisGrid.width * thisGrid.height - 1;
        for (intptr_t actionIndex = maxActionCount - 1; matInd > 0; actionIndex--) {
            aln_assert(actionIndex >= 0);
            aln_NWEntry entry = thisGrid.ptr[matInd];

            switch (entry.cameFromDir) {
                case aln_CameFromDir_TopLeft: {
                    alignedStr.actions[actionIndex] = aln_AlignAction_Match;
                    matInd -= (thisGrid.width + 1);
                } break;
                case aln_CameFromDir_Top: {
                    alignedStr.actions[actionIndex] = aln_AlignAction_GapRef;
                    matInd -= thisGrid.width;
                } break;
                case aln_CameFromDir_Left: {
                    alignedStr.actions[actionIndex] = aln_AlignAction_GapStr;
                    matInd -= 1;
                } break;
            }

            alignedStr.actionCount += 1;
        }

        {
            intptr_t actionsEmpty = maxActionCount - alignedStr.actionCount;
            aln_assert(actionsEmpty);
            alignedStr.actions += actionsEmpty;
        }

        alignedSeqs[strInd] = alignedStr;
    }  // for str

    aln_endTempMemory(temp);
    aln_AlignResult result = {
        .alignments = {alignedSeqs, strings.len},
        .matrices = {storedMatrices, config.storeFinalMatrices ? strings.len : 0},
    };
    return result;
}

typedef struct aln_StrBuilder {
    aln_Str    str;
    aln_Arena* arena;
} aln_StrBuilder;

aln_PRIVATEAPI aln_StrBuilder
aln_strBegin(aln_Arena* arena) {
    aln_assert(!arena->lockedForStr);
    arena->lockedForStr = true;
    aln_StrBuilder builder = {{aln_arenaFreePtr(arena)}, arena};
    return builder;
}

aln_PRIVATEAPI void
aln_strBuilderAddChar(aln_StrBuilder* builder, char ch) {
    aln_arenaChangeUsed(builder->arena, 1);
    builder->str.ptr[builder->str.len++] = ch;
}

aln_PRIVATEAPI aln_Str
aln_strEnd(aln_StrBuilder* builder) {
    aln_assert(builder->arena->lockedForStr);
    builder->arena->lockedForStr = false;
    return builder->str;
}

aln_PUBLICAPI aln_Str
aln_reconstruct(aln_Alignment aligned, aln_Reconstruct which, aln_Str reference, aln_Str ogstr, aln_Arena* arena) {
    aln_StrBuilder builder = aln_strBegin(arena);

    aln_Str target = which == aln_Reconstruct_Ref ? reference : ogstr;

    {
        intptr_t cur = 0;
        for (intptr_t actionIndex = 0; actionIndex < aligned.actionCount; actionIndex++) {
            aln_AlignAction targetGapAction = which == aln_Reconstruct_Ref ? aln_AlignAction_GapRef : aln_AlignAction_GapStr;
            aln_AlignAction nontargetGapAction = which == aln_Reconstruct_Ref ? aln_AlignAction_GapStr : aln_AlignAction_GapRef;
            aln_AlignAction thisAction = aligned.actions[actionIndex];
            if (thisAction == aln_AlignAction_Match || thisAction == nontargetGapAction) {
                aln_strBuilderAddChar(&builder, target.ptr[cur++]);
            } else if (thisAction == targetGapAction) {
                aln_strBuilderAddChar(&builder, '-');
            }
        }

        while (cur < target.len) {
            aln_assert(cur < target.len);
            aln_strBuilderAddChar(&builder, target.ptr[cur]);
            cur += 1;
        }
    }

    aln_Str result = aln_strEnd(&builder);
    return result;
}

aln_PUBLICAPI aln_ReconstructToCommonRefResult
aln_reconstructToCommonRef(aln_AlignmentArray alignments, aln_Str reference, aln_StrArray strings, aln_Memory* memory) {
    aln_assert(alignments.len == strings.len);
    aln_TempMemory temp = aln_beginTempMemory(&memory->temp);

    intptr_t* refGaps = aln_arenaAllocArray(&memory->temp, intptr_t, reference.len + 1);

    for (intptr_t alnIndex = 0; alnIndex < alignments.len; alnIndex++) {
        aln_Alignment aligned = alignments.ptr[alnIndex];

        {
            intptr_t cur = 0;
            intptr_t curGaps = 0;
            for (intptr_t actionIndex = 0; actionIndex < aligned.actionCount + 1; actionIndex++) {
                if (actionIndex == aligned.actionCount || aligned.actions[actionIndex] != aln_AlignAction_GapRef) {
                    aln_assert(cur < reference.len + 1);
                    refGaps[cur] = aln_max(refGaps[cur], curGaps);
                    cur += 1;
                    curGaps = 0;
                } else {
                    curGaps += 1;
                }
            }
        }
    }

    aln_Str* reconStrs = aln_arenaAllocArray(&memory->perm, aln_Str, strings.len);
    for (intptr_t alnIndex = 0; alnIndex < alignments.len; alnIndex++) {
        aln_Alignment aligned = alignments.ptr[alnIndex];
        aln_Str       thisStr = strings.ptr[alnIndex];

        aln_StrBuilder builder = aln_strBegin(&memory->perm);

        {
            intptr_t curStr = 0;
            intptr_t curRef = 0;
            intptr_t curGapsInRef = 0;
            for (intptr_t actionIndex = 0; actionIndex < aligned.actionCount + 1; actionIndex++) {
                if (actionIndex == aligned.actionCount || aligned.actions[actionIndex] != aln_AlignAction_GapRef) {
                    intptr_t extraGaps = refGaps[curRef] - curGapsInRef;

                    if (actionIndex < aligned.actionCount) {
                        while (extraGaps > 0) {
                            aln_strBuilderAddChar(&builder, '-');
                            extraGaps -= 1;
                        }
                    }

                    while (curGapsInRef > 0) {
                        aln_assert(curStr < thisStr.len);
                        aln_strBuilderAddChar(&builder, thisStr.ptr[curStr++]);
                        curGapsInRef -= 1;
                    }

                    if (actionIndex == aligned.actionCount) {
                        while (extraGaps > 0) {
                            aln_strBuilderAddChar(&builder, '-');
                            extraGaps -= 1;
                        }
                    }

                    curRef += 1;

                    if (actionIndex < aligned.actionCount) {
                        char ch = '-';
                        if (aligned.actions[actionIndex] == aln_AlignAction_Match) {
                            aln_assert(curStr < thisStr.len);
                            ch = thisStr.ptr[curStr++];
                        }
                        aln_strBuilderAddChar(&builder, ch);
                    }
                } else {
                    curGapsInRef += 1;
                }
            }
        }

        reconStrs[alnIndex] = aln_strEnd(&builder);
    }

    aln_Str commonRef = {.ptr = 0, .len = 0};
    {
        aln_StrBuilder builder = aln_strBegin(&memory->perm);

        for (intptr_t refGapIndex = 0; refGapIndex < reference.len + 1; refGapIndex++) {
            while (refGaps[refGapIndex]) {
                aln_strBuilderAddChar(&builder, '-');
                refGaps[refGapIndex] -= 1;
            }
            if (refGapIndex < reference.len) {
                aln_strBuilderAddChar(&builder, reference.ptr[refGapIndex]);
            }
        }

        commonRef = aln_strEnd(&builder);
    }

    aln_endTempMemory(temp);
    aln_ReconstructToCommonRefResult result = {.commonRef = commonRef, .alignedStrs = {reconStrs, strings.len}};
    return result;
}

aln_PUBLICAPI intptr_t
aln_matrix2index(intptr_t matrixWidth, intptr_t matrixHeight, intptr_t row, intptr_t col) {
// TODO(sen) Enable
#if 0 
    // NOTE(sen) Diagonal-first storage. Diagonals start in the top-left corner and are oriented bottomleft to topright.
    intptr_t topleftRect = (row + 1) * col;
    intptr_t toprightTriangleColCount = aln_min(row, matrixWidth - col - 1);
    intptr_t toprightTriangle = toprightTriangleColCount * (row + (row - toprightTriangleColCount + 1)) / 2;
    intptr_t bottomleftTriangleColCount = aln_min(col, matrixHeight - row - 1);
    intptr_t bottomleftTriangle = bottomleftTriangleColCount * (col + (col - bottomleftTriangleColCount + 1)) / 2;
    intptr_t result = topleftRect + toprightTriangle + bottomleftTriangle;
#else
    aln_unused(matrixHeight);
    intptr_t result = row * matrixWidth + col;
#endif
    return result;
}

#endif  // aln_IMPLEMENTATION

#ifdef __clang__
#pragma clang diagnostic pop
#endif

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
