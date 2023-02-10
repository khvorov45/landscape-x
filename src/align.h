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

typedef struct aln_Str {
    char*    ptr;
    intptr_t len;
} aln_Str;

typedef struct aln_Config {
    void*    outmem;
    intptr_t outmemBytes;
    void*    tempmem;
    intptr_t tempmemBytes;
    bool     storeFinalMatrices;
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

typedef enum aln_AlignAction {
    aln_AlignAction_Match,
    aln_AlignAction_GapStr,
    aln_AlignAction_GapRef,
} aln_AlignAction;

typedef struct aln_Alignment {
    aln_AlignAction* actions;
    intptr_t         actionCount;
} aln_Alignment;

typedef struct aln_AlignResult {
    aln_Alignment* strs;
    intptr_t       strCount;
    aln_Matrix2NW* matrices;
    intptr_t       bytesWrittenToOutput;
} aln_AlignResult;

typedef enum aln_Reconstruct {
    aln_Reconstruct_Ref,
    aln_Reconstruct_Str,
} aln_Reconstruct;

typedef struct aln_ReconstructToCommonRefResult {
    aln_Str  commonRef;
    aln_Str* strs;
    intptr_t strCount;
    intptr_t bytesWritten;
} aln_ReconstructToCommonRefResult;

// TODO(sen) StringArray?
// TODO(sen) Pass arenas instead of buffers and have some helpers to create the arenas?

aln_PUBLICAPI aln_AlignResult                  aln_align(aln_Str* references, intptr_t refCount, aln_Str* strings, intptr_t stringCount, aln_Config config);
aln_PUBLICAPI intptr_t                         aln_matrix2index(intptr_t matrixWidth, intptr_t matrixHeight, intptr_t row, intptr_t col);
aln_PUBLICAPI aln_Str                          aln_reconstruct(aln_Alignment aligned, aln_Reconstruct which, aln_Str reference, aln_Str ogstr, void* buf, intptr_t bufBytes);
aln_PUBLICAPI aln_ReconstructToCommonRefResult aln_reconstructToCommonRef(aln_Alignment* alignments, intptr_t alignmentsCount, aln_Str reference, aln_Str* strings, intptr_t stringCount, void* buf, intptr_t bufBytes);

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
#define aln_STR(x) (aln_Str) {x, aln_strlen(x)}

#ifndef aln_assert
#define aln_assert(condition)
#endif
// clang-format on

#ifndef aln_PRIVATEAPI
#define aln_PRIVATEAPI static
#endif

aln_PRIVATEAPI intptr_t
aln_strlen(char* str) {
    intptr_t result = 0;
    while (str[result] != '\0') {
        result += 1;
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

// TODO(sen) Don't assert on out of memory?
typedef struct aln_Arena {
    void*    base;
    intptr_t size;
    intptr_t used;
    intptr_t tempCount;
} aln_Arena;

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

aln_PUBLICAPI aln_AlignResult
aln_align(aln_Str* references, intptr_t refCount, aln_Str* strings, intptr_t stringCount, aln_Config config) {
    // NOTE(sen) Input validation
    {
        aln_assert(strings);
        aln_assert(stringCount > 0);
        aln_assert(references);
        aln_assert(refCount == 1 || refCount == stringCount);
        aln_assert(config.outmem);
        aln_assert(config.outmemBytes > 0);
        aln_assert(config.tempmem);
        aln_assert(config.tempmemBytes > 0);
        for (intptr_t strInd = 0; strInd < stringCount; strInd++) {
            aln_Str ogstr = strings[strInd];
            aln_assert(ogstr.ptr && ogstr.len > 0);
        }
        for (intptr_t refInd = 0; refInd < refCount; refInd++) {
            aln_Str ref = references[refInd];
            aln_assert(ref.ptr && ref.len > 0);
        }
    }

    aln_Arena  outputArena_ = {.base = config.outmem, .size = config.outmemBytes};
    aln_Arena* outputArena = &outputArena_;
    aln_Arena  arena_ = {.base = config.tempmem, .size = config.tempmemBytes};
    aln_Arena* arena = &arena_;

    aln_Matrix2NW maxGrid = {.ptr = 0, .width = 0, .height = 0};
    {
        intptr_t longestInputLen = 0;
        for (intptr_t strInd = 0; strInd < stringCount; strInd++) {
            longestInputLen = aln_max(longestInputLen, strings[strInd].len);
        }

        intptr_t longestRefLen = 0;
        for (intptr_t refInd = 0; refInd < refCount; refInd++) {
            longestRefLen = aln_max(longestRefLen, references[refInd].len);
        }

        maxGrid = aln_arenaAllocMatrix2(aln_Matrix2NW, aln_NWEntry, arena, longestRefLen + 1, longestInputLen + 1);
    }

    aln_Matrix2NW* storedMatrices = 0;
    if (config.storeFinalMatrices) {
        storedMatrices = aln_arenaAllocArray(outputArena, aln_Matrix2NW, stringCount);
    }

    aln_Alignment* alignedSeqs = aln_arenaAllocArray(outputArena, aln_Alignment, stringCount);
    for (intptr_t strInd = 0; strInd < stringCount; strInd++) {
        aln_Str ogstr = strings[strInd];
        aln_Str thisRef = references[0];
        if (refCount > 1) {
            thisRef = references[strInd];
        }

        aln_Matrix2NW thisGrid = maxGrid;
        thisGrid.width = thisRef.len + 1;
        thisGrid.height = ogstr.len + 1;
        aln_assert(thisGrid.width <= maxGrid.width);
        aln_assert(thisGrid.height <= maxGrid.height);

        aln_matrix2get(thisGrid, 0, 0).score = 0;
        for (intptr_t colIndex = 1; colIndex < thisGrid.width; colIndex++) {
            aln_matrix2get(thisGrid, 0, colIndex) = (aln_NWEntry) {-(float)colIndex, aln_CameFromDir_Left};
        }
        for (intptr_t rowIndex = 1; rowIndex < thisGrid.height; rowIndex++) {
            aln_matrix2get(thisGrid, rowIndex, 0) = (aln_NWEntry) {-(float)rowIndex, aln_CameFromDir_Top};
        }

        // TODO(khvorov) Diagonal iteration
        for (intptr_t diagonalIndex = 2; diagonalIndex < thisGrid.width + thisGrid.height - 1; diagonalIndex++) {
            intptr_t entriesInDiagonal = 0;
            if (diagonalIndex < thisGrid.height) {
                entriesInDiagonal = aln_min(diagonalIndex + 1, thisGrid.width);
            } else {
                entriesInDiagonal = aln_min(thisGrid.width - (diagonalIndex - thisGrid.height + 1), thisGrid.height);
            }
            aln_unused(entriesInDiagonal);
        }

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

                float indelScore = -1.0f;

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
            aln_Matrix2NW matCopy = aln_arenaAllocMatrix2(aln_Matrix2NW, aln_NWEntry, outputArena, thisGrid.width, thisGrid.height);
            for (intptr_t matIndex = 0; matIndex < matCopy.width * matCopy.height; matIndex++) {
                matCopy.ptr[matIndex] = thisGrid.ptr[matIndex];
            }
            storedMatrices[strInd] = matCopy;
        }

        // TODO(khvorov) This reconstruction would have to go by diagonal as well
        intptr_t      maxActionCount = (thisGrid.width - 1) + (thisGrid.height - 1);
        aln_Alignment alignedStr = {aln_arenaAllocArray(outputArena, aln_AlignAction, maxActionCount)};

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

    aln_AlignResult result = {.bytesWrittenToOutput = outputArena->used, .strs = alignedSeqs, .strCount = stringCount, .matrices = storedMatrices};
    return result;
}

typedef struct aln_StrBuilder {
    aln_Str  str;
    intptr_t capacity;
} aln_StrBuilder;

aln_PRIVATEAPI void
aln_strBuilderAddChar(aln_StrBuilder* builder, char ch) {
    aln_assert(builder->capacity > builder->str.len);
    builder->str.ptr[builder->str.len++] = ch;
}

aln_PUBLICAPI aln_Str
aln_reconstruct(aln_Alignment aligned, aln_Reconstruct which, aln_Str reference, aln_Str ogstr, void* buf, intptr_t bufBytes) {
    aln_StrBuilder builder = {{buf}, bufBytes};

    aln_Str target = which == aln_Reconstruct_Ref ? reference : ogstr;

    {
        intptr_t cur = 0;
        for (intptr_t actionIndex = 0; actionIndex < aligned.actionCount; actionIndex++) {
            aln_AlignAction targetAction = which == aln_Reconstruct_Ref ? aln_AlignAction_GapRef : aln_AlignAction_GapStr;
            aln_AlignAction thisAction = aligned.actions[actionIndex];
            if (thisAction == aln_AlignAction_Match) {
                aln_strBuilderAddChar(&builder, target.ptr[cur++]);
            } else if (thisAction == targetAction) {
                aln_strBuilderAddChar(&builder, '-');
            }
        }

        while (cur < target.len) {
            aln_assert(cur < target.len);
            aln_strBuilderAddChar(&builder, target.ptr[cur]);
            cur += 1;
        }
    }

    return builder.str;
}

aln_PUBLICAPI aln_ReconstructToCommonRefResult
aln_reconstructToCommonRef(aln_Alignment* alignments, intptr_t alignmentsCount, aln_Str reference, aln_Str* strings, intptr_t stringCount, void* buf, intptr_t bufBytes) {
    aln_assert(alignmentsCount == stringCount);

    aln_Arena  arena_ = {.base = buf, .size = bufBytes};
    aln_Arena* arena = &arena_;

    // TODO(sen) Temp arena
    intptr_t* refGaps = aln_arenaAllocArray(arena, intptr_t, reference.len + 1);

    for (intptr_t alnIndex = 0; alnIndex < alignmentsCount; alnIndex++) {
        aln_Alignment aligned = alignments[alnIndex];

        {
            intptr_t cur = 0;
            intptr_t curGaps = 0;
            for (intptr_t actionIndex = 0; actionIndex < aligned.actionCount; actionIndex++) {
                aln_AlignAction thisAction = aligned.actions[actionIndex];
                switch (thisAction) {
                    case aln_AlignAction_Match:
                    case aln_AlignAction_GapStr:
                        refGaps[cur] = aln_max(refGaps[cur], curGaps);
                        cur += 1;
                        curGaps = 0;
                        break;
                    case aln_AlignAction_GapRef:
                        curGaps += 1;
                        break;
                }
            }
            // TODO(sen) Do something about this repetition
            refGaps[cur] = aln_max(refGaps[cur], curGaps);
        }
    }

    aln_Str* reconStrs = aln_arenaAllocArray(arena, aln_Str, stringCount);
    for (intptr_t alnIndex = 0; alnIndex < alignmentsCount; alnIndex++) {
        aln_Alignment aligned = alignments[alnIndex];
        aln_Str       thisStr = strings[alnIndex];

        aln_StrBuilder builder = {.str.ptr = aln_arenaFreePtr(arena), .capacity = aln_arenaFreeSize(arena)};

        {
            intptr_t curStr = 0;
            intptr_t curRef = 0;
            intptr_t curGapsInRef = 0;
            for (intptr_t actionIndex = 0; actionIndex < aligned.actionCount; actionIndex++) {
                aln_AlignAction thisAction = aligned.actions[actionIndex];
                switch (thisAction) {
                    case aln_AlignAction_Match:
                    case aln_AlignAction_GapStr: {
                        intptr_t extraGaps = refGaps[curRef] - curGapsInRef;
                        while (extraGaps > 0) {
                            aln_strBuilderAddChar(&builder, '-');
                            extraGaps -= 1;
                        }
                        while (curGapsInRef > 0) {
                            aln_assert(curStr < thisStr.len);
                            aln_strBuilderAddChar(&builder, thisStr.ptr[curStr++]);
                            curGapsInRef -= 1;
                        }
                        curRef += 1;

                        char ch = '-';
                        if (thisAction == aln_AlignAction_Match) {
                            aln_assert(curStr < thisStr.len);
                            ch = thisStr.ptr[curStr++];
                        }
                        aln_strBuilderAddChar(&builder, ch);
                    } break;
                    case aln_AlignAction_GapRef:
                        curGapsInRef += 1;
                        break;
                }
            }
            // TODO(sen) Do something about this repetition
            intptr_t extraGaps = refGaps[curRef] - curGapsInRef;
            while (extraGaps > 0) {
                aln_strBuilderAddChar(&builder, '-');
                extraGaps -= 1;
            }
            while (curGapsInRef > 0) {
                aln_assert(curStr < thisStr.len);
                aln_strBuilderAddChar(&builder, thisStr.ptr[curStr++]);
                curGapsInRef -= 1;
            }
        }

        aln_arenaChangeUsed(arena, builder.str.len);
        reconStrs[alnIndex] = builder.str;
    }

    aln_Str commonRef = {.ptr = 0, .len = 0};
    {
        aln_StrBuilder builder = {.str.ptr = aln_arenaFreePtr(arena), .capacity = aln_arenaFreeSize(arena)};

        for (intptr_t refGapIndex = 0; refGapIndex < reference.len + 1; refGapIndex++) {
            while (refGaps[refGapIndex]) {
                aln_strBuilderAddChar(&builder, '-');
                refGaps[refGapIndex] -= 1;
            }
            if (refGapIndex < reference.len) {
                aln_strBuilderAddChar(&builder, reference.ptr[refGapIndex]);
            }
        }

        aln_arenaChangeUsed(arena, builder.str.len);
        commonRef = builder.str;
    }

    aln_ReconstructToCommonRefResult result = {.bytesWritten = arena->used, .commonRef = commonRef, .strs = reconStrs, .strCount = stringCount};
    return result;
}

#endif  // aln_IMPLEMENTATION

#ifdef __clang__
#pragma clang diagnostic pop
#endif

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
