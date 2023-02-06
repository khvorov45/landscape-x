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

#define aln_matrix2index(matrix, row, col) (((row)*matrix.width) + (col))
#define aln_matrix2get(matrix, row, col) matrix.ptr[aln_matrix2index(matrix, row, col)]

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

typedef struct aln_AlignResult {
    aln_Str*       strs;
    intptr_t       strCount;
    aln_Matrix2NW* matrices;
    intptr_t       bytesWrittenToOutput;
} aln_AlignResult;

aln_PUBLICAPI aln_AlignResult aln_align(aln_Str reference, aln_Str* strings, intptr_t stringCount, aln_Config config);

#endif  // aln_HEADER

#ifdef aln_IMPLEMENTATION

// clang-format off
#define aln_max(a, b) (((a) > (b)) ? (a) : (b))
#define aln_min(a, b) (((a) < (b)) ? (a) : (b))
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

aln_PUBLICAPI aln_AlignResult
aln_align(aln_Str reference, aln_Str* strings, intptr_t stringCount, aln_Config config) {
    aln_assert(reference.ptr && reference.len > 0);
    aln_assert(strings);
    aln_assert(stringCount > 0);
    aln_assert(config.outmem);
    aln_assert(config.outmemBytes > 0);
    aln_assert(config.tempmem);
    aln_assert(config.tempmemBytes > 0);
    for (intptr_t strInd = 0; strInd < stringCount; strInd++) {
        aln_Str ogstr = strings[strInd];
        aln_assert(ogstr.ptr && ogstr.len > 0);
    }

    aln_Arena  outputArena_ = {.base = config.outmem, .size = config.outmemBytes};
    aln_Arena* outputArena = &outputArena_;
    aln_Arena  arena_ = {.base = config.tempmem, .size = config.tempmemBytes};
    aln_Arena* arena = &arena_;

    intptr_t longestInputLen = reference.len;
    for (intptr_t strInd = 0; strInd < stringCount; strInd++) {
        longestInputLen = aln_max(longestInputLen, strings[strInd].len);
    }

    aln_Matrix2NW maxGrid = aln_arenaAllocMatrix2(aln_Matrix2NW, aln_NWEntry, arena, reference.len + 1, longestInputLen + 1);

    aln_Matrix2NW* storedMatrices = 0;
    if (config.storeFinalMatrices) {
        storedMatrices = aln_arenaAllocArray(outputArena, aln_Matrix2NW, stringCount);
    }

    aln_Str* alignedSeqs = aln_arenaAllocArray(outputArena, aln_Str, stringCount);
    for (intptr_t strInd = 0; strInd < stringCount; strInd++) {
        aln_Str ogstr = strings[strInd];

        aln_Matrix2NW thisGrid = maxGrid;
        thisGrid.height = ogstr.len + 1;
        aln_assert(thisGrid.height <= maxGrid.height);

        aln_matrix2get(thisGrid, 0, 0).score = 0;
        for (intptr_t colIndex = 1; colIndex < thisGrid.width; colIndex++) {
            aln_matrix2get(thisGrid, 0, colIndex) = (aln_NWEntry) {-(float)colIndex, aln_CameFromDir_Left};
        }
        for (intptr_t rowIndex = 1; rowIndex < thisGrid.height; rowIndex++) {
            aln_matrix2get(thisGrid, rowIndex, 0) = (aln_NWEntry) {-(float)rowIndex, aln_CameFromDir_Top};
        }

        intptr_t thisGridMaxScoreIndex = 0;
        float    thisGridMaxScore = 0;
        for (intptr_t rowIndex = 1; rowIndex < thisGrid.height; rowIndex++) {
            for (intptr_t colIndex = 1; colIndex < thisGrid.width; colIndex++) {
                intptr_t topleftIndex = aln_matrix2index(thisGrid, rowIndex - 1, colIndex - 1);
                intptr_t topIndex = aln_matrix2index(thisGrid, rowIndex - 1, colIndex);
                intptr_t leftIndex = aln_matrix2index(thisGrid, rowIndex, colIndex - 1);

                float alignScore = 1.0f;
                {
                    char chRef = reference.ptr[colIndex - 1];
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

                if (thisGridMaxScore < maxScore) {
                    thisGridMaxScore = maxScore;
                    thisGridMaxScoreIndex = aln_matrix2index(thisGrid, rowIndex, colIndex);
                }
            }  // for mat col
        }  // for mat row

        if (config.storeFinalMatrices) {
            aln_Matrix2NW matCopy = aln_arenaAllocMatrix2(aln_Matrix2NW, aln_NWEntry, outputArena, thisGrid.width, thisGrid.height);
            for (intptr_t matIndex = 0; matIndex < matCopy.width * matCopy.height; matIndex++) {
                matCopy.ptr[matIndex] = thisGrid.ptr[matIndex];
            }
            storedMatrices[strInd] = matCopy;
        }

        aln_Str  alignedStr = {aln_arenaAllocArray(outputArena, char, reference.len), reference.len};
        intptr_t alignedStrIndex = (thisGridMaxScoreIndex % thisGrid.width) - 1;

        for (intptr_t padIndex = alignedStrIndex + 1; padIndex < reference.len; padIndex++) {
            alignedStr.ptr[padIndex] = '-';
        }

        for (intptr_t matInd = thisGridMaxScoreIndex; matInd > 0 && alignedStrIndex >= 0;) {
            aln_NWEntry entry = thisGrid.ptr[matInd];

            switch (entry.cameFromDir) {
                case aln_CameFromDir_TopLeft: {
                    intptr_t row = matInd / thisGrid.width;
                    alignedStr.ptr[alignedStrIndex] = ogstr.ptr[row - 1];
                    matInd -= (thisGrid.width + 1);
                } break;
                case aln_CameFromDir_Top: {
                    alignedStr.ptr[alignedStrIndex] = '!';
                    matInd -= thisGrid.width;
                } break;
                case aln_CameFromDir_Left: {
                    alignedStr.ptr[alignedStrIndex] = '-';
                    matInd -= 1;
                } break;
            }

            alignedStrIndex -= 1;
        }

        alignedSeqs[strInd] = alignedStr;
    }  // for str

    aln_AlignResult result = {.bytesWrittenToOutput = outputArena->used, .strs = alignedSeqs, .strCount = stringCount, .matrices = storedMatrices};
    return result;
}

#endif  // aln_IMPLEMENTATION

#ifdef __clang__
#pragma clang diagnostic pop
#endif

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
