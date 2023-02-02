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

typedef enum aln_Status {
    aln_Failure,
    aln_Success,
} aln_Status;

typedef struct aln_Str {
    char*    ptr;
    intptr_t len;
} aln_Str;

typedef struct aln_AlignResult {
    aln_Str* seqs;
    intptr_t seqCount;
    intptr_t bytesWrittenToOutput;
} aln_AlignResult;

aln_PUBLICAPI aln_AlignResult aln_align(aln_Str* strings, intptr_t stringCount, void* outmem, intptr_t outmemBytes, void* tempmem, intptr_t tempmemBytes);

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
#define aln_arenaAllocMatrix2I32(arena, width, height) aln_arenaAllocMatrix2(aln_Matrix2I32, intptr_t, arena, width, height)
#define aln_arenaAllocMatrix2F32(arena, width, height) aln_arenaAllocMatrix2(aln_Matrix2F32, float, arena, width, height)
#define aln_arenaAllocMatrix2U8(arena, width, height) aln_arenaAllocMatrix2(aln_Matrix2U8, uint8_t, arena, width, height)
#define aln_matrix2get(matrix, row, col) matrix.ptr[((row) * matrix.width) + (col)] 
#define aln_isPowerOf2(x) (((x) > 0) && (((x) & ((x)-1)) == 0))
#define aln_unused(x) ((x) = (x))

#ifndef aln_assert
#define aln_assert(condition)
#endif
// clang-format on

#ifndef aln_PRIVATEAPI
#define aln_PRIVATEAPI static
#endif

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
aln_align(aln_Str* strings, intptr_t stringCount, void* outmem, intptr_t outmemBytes, void* tempmem, intptr_t tempmemBytes) {
    aln_Arena  arena_ = {.base = tempmem, .size = tempmemBytes};
    aln_Arena* arena = &arena_;

    aln_unused(arena);

    aln_Arena outputArena = {.base = outmem, .size = outmemBytes};
    aln_Str*  alignedSeqs = aln_arenaAllocArray(&outputArena, aln_Str, stringCount);
    for (intptr_t strInd = 0; strInd < stringCount; strInd++) {
        aln_Str ogstr = strings[strInd];
        char* alignedSeq = aln_arenaAllocArray(&outputArena, char, ogstr.len);
        for (intptr_t charInd = 0; charInd < ogstr.len; charInd++) {
            alignedSeq[charInd] = ogstr.ptr[charInd];
        }
        alignedSeqs[strInd] = (aln_Str) {alignedSeq, ogstr.len};
    }

    aln_AlignResult result = {.bytesWrittenToOutput = outputArena.used, .seqs = alignedSeqs, .seqCount = stringCount};
    return result;
}

#endif  // aln_IMPLEMENTATION

#ifdef __clang__
#pragma clang diagnostic pop
#endif

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
