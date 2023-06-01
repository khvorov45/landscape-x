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
#include <float.h>

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

typedef struct aln_Matrix2f {
    float*   ptr;
    intptr_t width;
    intptr_t height;
} aln_Matrix2f;

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

typedef struct aln_Rng {
    uint64_t state;
    // Must be odd
    uint64_t inc;
} aln_Rng;

typedef struct aln_StrModSpec {
    float   trimStartMaxProp;
    float   trimEndMaxProp;
    float   mutationProb;
    float   deletionProb;
    float   insertionProb;
    aln_Str insertionSrc;
} aln_StrModSpec;

typedef struct aln_TreeNode {
    bool isInternal;
} aln_TreeNode;

typedef struct aln_TreeNodeArray {
    aln_TreeNode* ptr;
    intptr_t      len;
} aln_TreeNodeArray;

typedef struct aln_TreeBranch {
    int32_t node1Index;
    int32_t node2Index;
    float   len;
} aln_TreeBranch;

typedef struct aln_TreeBranchArray {
    aln_TreeBranch* ptr;
    intptr_t        len;
} aln_TreeBranchArray;

typedef struct aln_Tree {
    aln_TreeNodeArray   nodes;
    aln_TreeBranchArray branches;
} aln_Tree;

aln_PUBLICAPI aln_Memory aln_createMemory(void* base, intptr_t totalSize, intptr_t permSize);
aln_PUBLICAPI aln_Memory aln_createMemory2(void* permBase, intptr_t permSize, void* tempBase, intptr_t tempSize);

aln_PUBLICAPI aln_AlignResult                  aln_align(aln_StrArray references, aln_StrArray strings, aln_Config config, aln_Memory* memory);
aln_PUBLICAPI aln_Str                          aln_reconstruct(aln_Alignment aligned, aln_Reconstruct which, aln_Str reference, aln_Str ogstr, aln_Arena* arena);
aln_PUBLICAPI aln_ReconstructToCommonRefResult aln_reconstructToCommonRef(aln_AlignmentArray alignments, aln_Str reference, aln_StrArray strings, aln_Memory* memory);

aln_PUBLICAPI aln_Rng  aln_createRng(uint32_t seed);
aln_PUBLICAPI uint32_t aln_randomU32(aln_Rng* rng);
aln_PUBLICAPI uint32_t aln_randomU32Bound(aln_Rng* rng, uint32_t max);
aln_PUBLICAPI float    aln_randomF3201(aln_Rng* rng);
aln_PUBLICAPI aln_Str  aln_randomString(aln_Rng* rng, aln_Str src, intptr_t len, aln_Arena* arena);
aln_PUBLICAPI aln_Str  aln_randomStringMod(aln_Rng* rng, aln_Str src, aln_StrModSpec spec, aln_Arena* arena);

aln_PUBLICAPI aln_Tree aln_createTree(aln_StrArray strings, aln_Memory* memory);

aln_PUBLICAPI intptr_t aln_matrix2index(intptr_t matrixWidth, intptr_t matrixHeight, intptr_t row, intptr_t col);

#endif  // aln_HEADER

#ifdef aln_IMPLEMENTATION

// clang-format off
#define aln_max(a, b) (((a) > (b)) ? (a) : (b))
#define aln_min(a, b) (((a) < (b)) ? (a) : (b))
#define aln_abs(a) (((a) < (0)) ? (-(a)) : (a))
#define aln_clamp(x, a, b) (((x) < (a)) ? (a) : (((x) > (b)) ? (b) : (x)))
#define aln_arrayCount(arr) (intptr_t)(sizeof(arr) / sizeof(arr[0]))
#define aln_arenaAllocArray(arena, type, len) (type*)aln_arenaAllocAndZero(arena, (len) * (intptr_t)sizeof(type), alignof(type))
#define aln_arenaAllocArrayNoZero(arena, type, len) (type*)aln_arenaAllocNoZero(arena, (len) * (intptr_t)sizeof(type), alignof(type))
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
aln_arenaAllocNoZero(aln_Arena* arena, intptr_t size, intptr_t align) {
    aln_arenaAlignFreePtr(arena, align);
    void* result = aln_arenaFreePtr(arena);
    aln_arenaChangeUsed(arena, size);
    return result;
}

aln_PRIVATEAPI void*
aln_arenaAllocAndZero(aln_Arena* arena, intptr_t size, intptr_t align) {
    void* result = aln_arenaAllocNoZero(arena, size, align);
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
    buf[str.len] = 0;
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
    aln_assert(row >= 0 && row < matrixHeight);
    aln_assert(col >= 0 && col < matrixWidth);
    intptr_t result = row * matrixWidth + col;
    return result;
}

aln_PUBLICAPI aln_Rng
aln_createRng(uint32_t seed) {
    aln_Rng rng = {.state = seed, .inc = seed | 1};
    // NOTE(khvorov) When seed is 0 the first 2 numbers are always 0 which is probably not what we want
    aln_randomU32(&rng);
    aln_randomU32(&rng);
    return rng;
}

// PCG-XSH-RR
// state_new = a * state_old + b
// output = rotate32((state ^ (state >> 18)) >> 27, state >> 59)
// as per `PCG: A Family of Simple Fast Space-Efficient Statistically Good Algorithms for Random Number Generation`
aln_PUBLICAPI uint32_t
aln_randomU32(aln_Rng* rng) {
    uint64_t state = rng->state;
    uint64_t xorWith = state >> 18u;
    uint64_t xored = state ^ xorWith;
    uint64_t shifted64 = xored >> 27u;
    uint32_t shifted32 = (uint32_t)shifted64;
    uint32_t rotateBy = state >> 59u;
    uint32_t shiftRightBy = rotateBy;
    uint32_t resultRight = shifted32 >> shiftRightBy;
    // NOTE(khvorov) This is `32 - rotateBy` but produces 0 when rotateBy is 0
    // Shifting a 32 bit value by 32 is apparently UB and the compiler is free to remove that code
    // I guess, so we are avoiding it by doing this weird bit hackery
    uint32_t shiftLeftBy = (-rotateBy) & 0x1Fu;
    uint32_t resultLeft = shifted32 << shiftLeftBy;
    uint32_t result = resultRight | resultLeft;
    // NOTE(khvorov) This is just one of those magic LCG constants "in common use"
    // https://en.wikipedia.org/wiki/Linear_congruential_generator#Parameters_in_common_use
    rng->state = 6364136223846793005ULL * state + rng->inc;
    return result;
}

aln_PUBLICAPI uint32_t
aln_randomU32Bound(aln_Rng* rng, uint32_t max) {
    // NOTE(khvorov) This is equivalent to (UINT32_MAX + 1) % max;
    uint32_t threshold = -max % max;
    uint32_t unbound = aln_randomU32(rng);
    while (unbound < threshold) {
        unbound = aln_randomU32(rng);
    }
    uint32_t result = unbound % max;
    return result;
}

aln_PUBLICAPI float
aln_randomF3201(aln_Rng* rng) {
    uint32_t randomU32 = aln_randomU32(rng);
    float    randomF32 = (float)randomU32;
    float    onePastMaxRandomU32 = (float)(1ULL << 32ULL);
    float    result = randomF32 / onePastMaxRandomU32;
    return result;
}

aln_PUBLICAPI aln_Str
aln_randomString(aln_Rng* rng, aln_Str src, intptr_t len, aln_Arena* arena) {
    char* ptr = aln_arenaAllocArrayNoZero(arena, char, len);
    for (intptr_t ind = 0; ind < len; ind++) {
        uint32_t index = aln_randomU32Bound(rng, src.len);
        ptr[ind] = src.ptr[index];
    }
    aln_Str result = {ptr, len};
    return result;
}

aln_PUBLICAPI aln_Str
aln_randomStringMod(aln_Rng* rng, aln_Str src, aln_StrModSpec spec, aln_Arena* arena) {
    intptr_t trimStartMax = (intptr_t)(spec.trimStartMaxProp * (float)src.len + 0.5f);
    intptr_t trimStart = aln_randomU32Bound(rng, trimStartMax + 1);

    intptr_t trimEndMax = (intptr_t)(spec.trimEndMaxProp * (float)src.len + 0.5f);
    intptr_t trimEnd = aln_randomU32Bound(rng, trimEndMax + 1);

    aln_StrBuilder builder = aln_strBegin(arena);

    for (intptr_t ind = trimStart; ind < src.len - trimEnd; ind++) {
        float r01 = aln_randomF3201(rng);

        bool mutate = false;
        bool delete = false;
        bool insert = false;
        if (r01 < spec.mutationProb) {
            mutate = true;
        } else if (r01 >= spec.mutationProb && r01 < spec.mutationProb + spec.deletionProb) {
            delete = true;
        } else if (r01 >= spec.mutationProb + spec.deletionProb && r01 < spec.mutationProb + spec.deletionProb + spec.insertionProb) {
            insert = true;
        }

        if (mutate || insert) {
            uint32_t insertionCharInd = aln_randomU32Bound(rng, spec.insertionSrc.len);
            aln_strBuilderAddChar(&builder, spec.insertionSrc.ptr[insertionCharInd]);
        } else if (!delete) {
            aln_strBuilderAddChar(&builder, src.ptr[ind]);
        }
    }

    aln_Str result = aln_strEnd(&builder);
    return result;
}

aln_PRIVATEAPI aln_Tree
aln_createTreeFromDistanceMatrix(aln_Matrix2f curDistMat, aln_Memory* memory) {
    aln_assert(curDistMat.height == curDistMat.width);
    aln_TempMemory temp = aln_beginTempMemory(&memory->temp);

    intptr_t ogmatside = curDistMat.width;

    aln_Tree tree = {.branches = {.ptr = 0, .len = 0}, .nodes = {.ptr = 0, .len = 0}};
    tree.nodes.len = ogmatside + ogmatside - 2;
    tree.nodes.ptr = aln_arenaAllocArray(&memory->perm, aln_TreeNode, tree.nodes.len);
    for (intptr_t ind = ogmatside; ind < tree.nodes.len; ind++) {
        tree.nodes.ptr[ind].isInternal = true;
    }
    tree.branches.len = (ogmatside - 2) * 2 + 1;
    tree.branches.ptr = aln_arenaAllocArray(&memory->perm, aln_TreeBranch, tree.branches.len);

    aln_Matrix2f nextDistMat = aln_arenaAllocMatrix2(aln_Matrix2f, float, &memory->temp, ogmatside, ogmatside);

    int32_t* curDistMatNodeIndices = aln_arenaAllocArray(&memory->temp, int32_t, curDistMat.width);
    for (intptr_t ind = 0; ind < curDistMat.width; ind++) {
        curDistMatNodeIndices[ind] = ind;
    }

    for (intptr_t joinIndex = 0; joinIndex < ogmatside - 2; joinIndex++) {
        aln_TempMemory joinTemp = aln_beginTempMemory(&memory->temp);

        float* distanceRowSum = aln_arenaAllocArray(&memory->temp, float, curDistMat.height);
        for (intptr_t row = 0; row < curDistMat.height; row++) {
            float sum = 0;
            for (intptr_t col = 0; col < curDistMat.width; col++) {
                sum += aln_matrix2get(curDistMat, row, col);
            }
            distanceRowSum[row] = sum;
        }

        bool     minIsSet = false;
        float    qabMin = FLT_MAX;
        intptr_t rowQabMin = 0;
        intptr_t colQabMin = 0;
        for (intptr_t row = 0; row < curDistMat.height - 1; row++) {
            for (intptr_t col = row + 1; col < curDistMat.width; col++) {
                float dab = aln_matrix2get(curDistMat, row, col);
                float t1 = (float)(curDistMat.width - 2) * dab;
                float t2 = distanceRowSum[row];
                float t3 = distanceRowSum[col];
                float qab = t1 - t2 - t3;
                if (qab < qabMin) {
                    qabMin = qab;
                    rowQabMin = row;
                    colQabMin = col;
                    minIsSet = true;
                }
            }
        }
        aln_assert(minIsSet);
        aln_assert(rowQabMin < colQabMin);

        {
            intptr_t rowQabMinNodeIndex = curDistMatNodeIndices[rowQabMin];
            intptr_t colQabMinNodeIndex = curDistMatNodeIndices[colQabMin];
            intptr_t newNodeIndex = ogmatside + joinIndex;

            aln_TreeBranch branch1 = {.node1Index = rowQabMinNodeIndex, .node2Index = newNodeIndex, .len = 0};
            aln_TreeBranch branch2 = {.node1Index = colQabMinNodeIndex, .node2Index = newNodeIndex, .len = 0};

            {
                float dfg = aln_matrix2get(curDistMat, rowQabMin, colQabMin);
                float t1 = 0.5 * dfg;
                float t2 = 1.0f / (float)(2 * (curDistMat.width - 2));
                float t3 = 0;
                float t4 = 0;
                for (intptr_t ind = 0; ind < curDistMat.width; ind++) {
                    t3 += aln_matrix2get(curDistMat, rowQabMin, ind);
                    t4 += aln_matrix2get(curDistMat, colQabMin, ind);
                }
                branch1.len = t1 + t2 * (t3 - t4);
                branch2.len = dfg - branch1.len;
            }

            intptr_t firstBranchIndex = joinIndex * 2;
            tree.branches.ptr[firstBranchIndex] = branch1;
            tree.branches.ptr[firstBranchIndex + 1] = branch2;

            if (joinIndex == ogmatside - 3) {
                aln_assert(curDistMat.width == 3 && curDistMat.height == 3);
                intptr_t leftoverQab = 0;
                if (rowQabMin == 0 || colQabMin == 0) {
                    leftoverQab = 1;
                }
                if (rowQabMin == 1 || colQabMin == 1) {
                    leftoverQab = 2;
                }
                aln_assert(rowQabMin != leftoverQab && colQabMin != leftoverQab);

                intptr_t leftoverQabNodeIndex = curDistMatNodeIndices[leftoverQab];
                float    dfg = aln_matrix2get(curDistMat, rowQabMin, leftoverQab);

                aln_TreeBranch branch3 = {.node1Index = leftoverQabNodeIndex, .node2Index = newNodeIndex, .len = dfg - branch1.len};
                aln_assert(firstBranchIndex + 2 < tree.branches.len);
                tree.branches.ptr[firstBranchIndex + 2] = branch3;
            }
        }

        if (curDistMat.height > 3) {
            nextDistMat.height = curDistMat.height - 1;
            nextDistMat.width = curDistMat.width - 1;

            intptr_t jstart = aln_min(rowQabMin, colQabMin);
            intptr_t jend = aln_max(rowQabMin, colQabMin) - 1;
            for (intptr_t row = 0; row < nextDistMat.height - 1; row++) {
                for (intptr_t col = row; col < nextDistMat.width; col++) {
                    intptr_t curRow = row + (intptr_t)(row >= jstart) + (intptr_t)(row >= jend);
                    intptr_t curCol = col + (intptr_t)(col >= jstart) + (intptr_t)(col >= jend);

                    if (curCol == curDistMat.width) {
                        float d1 = aln_matrix2get(curDistMat, rowQabMin, curRow);
                        float d2 = aln_matrix2get(curDistMat, colQabMin, curRow);
                        float d3 = aln_matrix2get(curDistMat, rowQabMin, colQabMin);
                        float val = 0.5 * (d1 + d2 - d3);
                        aln_matrix2get(nextDistMat, row, col) = val;
                        aln_matrix2get(nextDistMat, col, row) = val;
                    } else {
                        float val = aln_matrix2get(curDistMat, curRow, curCol);
                        aln_matrix2get(nextDistMat, row, col) = val;
                        aln_matrix2get(nextDistMat, col, row) = val;
                    }
                }
            }
            aln_matrix2get(nextDistMat, nextDistMat.height - 1, nextDistMat.width - 1) = 0;

            {
                aln_Matrix2f temp = curDistMat;
                curDistMat = nextDistMat;
                nextDistMat = temp;
            }

            for (intptr_t ind = jstart; ind < jend; ind++) {
                curDistMatNodeIndices[ind] = curDistMatNodeIndices[ind + 1];
            }
            for (intptr_t ind = jend; ind < curDistMat.width - 1; ind++) {
                curDistMatNodeIndices[ind] = curDistMatNodeIndices[ind + 2];
            }
            curDistMatNodeIndices[curDistMat.width - 1] = ogmatside + joinIndex;
        }

        aln_endTempMemory(joinTemp);
    }

    aln_endTempMemory(temp);
    return tree;
}

aln_PUBLICAPI aln_Tree
aln_createTree(aln_StrArray strings, aln_Memory* memory) {
    aln_TempMemory temp = aln_beginTempMemory(&memory->temp);

    aln_Matrix2f curDistMat = aln_arenaAllocMatrix2(aln_Matrix2f, float, &memory->temp, strings.len, strings.len);
    for (intptr_t row = 0; row < curDistMat.height; row++) {
        for (intptr_t col = row; col < curDistMat.width; col++) {
            aln_Str s1 = strings.ptr[row];
            aln_Str s2 = strings.ptr[col];
            aln_assert(s1.len == s2.len);
            float distance = 0.0f;
            for (intptr_t ind = 0; ind < s1.len; ind++) {
                char c1 = s1.ptr[ind];
                char c2 = s2.ptr[ind];
                if (c1 != c2) {
                    distance += 1.0f;
                }
            }
            aln_matrix2get(curDistMat, row, col) = distance;
            aln_matrix2get(curDistMat, col, row) = distance;
        }
    }

    aln_Tree tree = aln_createTreeFromDistanceMatrix(curDistMat, memory);

    aln_endTempMemory(temp);
    return tree;
}

#endif  // aln_IMPLEMENTATION

#ifdef __clang__
#pragma clang diagnostic pop
#endif

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
