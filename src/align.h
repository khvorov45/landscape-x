#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
#endif

#ifndef aln_HEADER
#define aln_HEADER

//
// SECTION Alignment
//

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
    intptr_t bytesWritten;
} aln_AlignResult;

// TODO(sen) Remove private forward decls once we pull everything into 1 TU

//
// SECTION Private
//

// clang-format off
#define aln_max(a, b) (((a) > (b)) ? (a) : (b))
#define aln_min(a, b) (((a) < (b)) ? (a) : (b))
#define aln_clamp(x, a, b) (((x) < (a)) ? (a) : (((x) > (b)) ? (b) : (x)))
#define aln_arrayCount(arr) (intptr_t)(sizeof(arr) / sizeof(arr[0]))
#define aln_arenaAllocArray(arena, type, len) (type*)aln_arenaAllocAndZero(arena, (len) * (int32_t)sizeof(type), alignof(type))
#define aln_arenaAllocStruct(arena, type) (type*)aln_arenaAllocAndZero(arena, sizeof(type), alignof(type))
#define aln_arenaAllocMatrix2(matrixType, dataType, arena, width_, height_) (matrixType) {.ptr = aln_arenaAllocArray(arena, dataType, width_ * height_), .width = width_, .height = height_}
#define aln_arenaAllocMatrix2I32(arena, width, height) aln_arenaAllocMatrix2(aln_Matrix2I32, int32_t, arena, width, height)
#define aln_arenaAllocMatrix2F32(arena, width, height) aln_arenaAllocMatrix2(aln_Matrix2F32, float, arena, width, height)
#define aln_arenaAllocMatrix2U8(arena, width, height) aln_arenaAllocMatrix2(aln_Matrix2U8, uint8_t, arena, width, height)
#define aln_matrix2get(matrix, row, col) matrix.ptr[((row) * matrix.width) + (col)] 
#define aln_isPowerOf2(x) (((x) > 0) && (((x) & ((x)-1)) == 0))
#define aln_unused(x) ((x) = (x))

#ifndef aln_assert
#define aln_assert(condition) do {if (condition) {} else {__builtin_debugtrap();}} while(0)
#endif
// clang-format on

#ifndef aln_PRIVATEAPI
#define aln_PRIVATEAPI static
#endif

typedef struct aln_Arena {
    void*    base;
    intptr_t size;
    intptr_t used;
    int32_t  tempCount;
} aln_Arena;

typedef struct aln_TempMemory {
    aln_Arena* arena;
    intptr_t   usedAtBegin;
    int32_t    tempCountAtBegin;
} aln_TempMemory;

typedef struct aln_Matrix2I32 {
    int32_t* ptr;
    int32_t  width;
    int32_t  height;
} aln_Matrix2I32;

typedef struct aln_Matrix2F32 {
    float*  ptr;
    int32_t width;
    int32_t height;
} aln_Matrix2F32;

typedef struct aln_Matrix2U8 {
    uint8_t* ptr;
    int32_t  width;
    int32_t  height;
} aln_Matrix2U8;

typedef struct aln_Context {
    int32_t        penalty;
    int32_t        offset;
    int32_t        constraint;
    float          minimumweight;
    float          fastathreshold;
    float          sueff_global;
    char           treemethod;
    int32_t        njob;
    int32_t        amino_n[256];
    aln_Matrix2I32 n_dis;
    aln_Matrix2F32 n_dis_consweight_multi;
    uint8_t        amino[256];
    int32_t        penalty_dist;
    int32_t        penalty_ex;
    int32_t        outgap;
    float          consweight_multi;
    int32_t        commonAlloc1;
    int32_t        commonAlloc2;
    int32_t**      commonIP;
    int32_t        nalphabets;
} aln_Context;

typedef struct aln_LocalHom {
    struct aln_LocalHom* next;
    int32_t              start1;
    int32_t              end1;
    int32_t              start2;
    int32_t              end2;
    double               opt;
    double               importance;
    double               rimportance;
} aln_LocalHom;

// TODO(sen) Don't assert on out of memory?

aln_PRIVATEAPI int32_t        aln_getOffsetForAlignment(void* ptr, int32_t align);
aln_PRIVATEAPI aln_Arena      aln_createArenaFromArena(aln_Arena* arena, intptr_t bytes);
aln_PRIVATEAPI void*          aln_arenaAllocAndZero(aln_Arena* arena, int32_t size, int32_t align);
aln_PRIVATEAPI void           aln_arenaAlignFreePtr(aln_Arena* arena, int32_t align);
aln_PRIVATEAPI void*          aln_arenaFreePtr(aln_Arena* arena);
aln_PRIVATEAPI intptr_t       aln_arenaFreeSize(aln_Arena* arena);
aln_PRIVATEAPI void           aln_arenaChangeUsed(aln_Arena* arena, intptr_t byteDelta);
aln_PRIVATEAPI aln_TempMemory aln_beginTempMemory(aln_Arena* arena);
aln_PRIVATEAPI void           aln_endTempMemory(aln_TempMemory temp);

#endif  // aln_HEADER

#ifdef aln_IMPLEMENTATION

aln_PRIVATEAPI int32_t
aln_getOffsetForAlignment(void* ptr, int32_t align) {
    aln_assert(aln_isPowerOf2(align));
    uintptr_t ptrAligned = (uintptr_t)((uint8_t*)ptr + (align - 1)) & (uintptr_t)(~(align - 1));
    aln_assert(ptrAligned >= (uintptr_t)ptr);
    intptr_t diff = (intptr_t)(ptrAligned - (uintptr_t)ptr);
    aln_assert(diff < align && diff >= 0);
    return (int32_t)diff;
}

aln_PRIVATEAPI aln_Arena
aln_createArenaFromArena(aln_Arena* parent, intptr_t bytes) {
    aln_Arena arena = {.base = aln_arenaFreePtr(parent), .size = bytes};
    aln_arenaChangeUsed(parent, bytes);
    return arena;
}

aln_PRIVATEAPI void*
aln_arenaAllocAndZero(aln_Arena* arena, int32_t size, int32_t align) {
    aln_arenaAlignFreePtr(arena, align);
    void* result = aln_arenaFreePtr(arena);
    aln_arenaChangeUsed(arena, size);
    for (int32_t ind = 0; ind < size; ind++) {
        ((uint8_t*)result)[ind] = 0;
    }
    return result;
}

aln_PRIVATEAPI void
aln_arenaAlignFreePtr(aln_Arena* arena, int32_t align) {
    int32_t offset = aln_getOffsetForAlignment(aln_arenaFreePtr(arena), align);
    aln_arenaChangeUsed(arena, offset);
}

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
    for (int32_t ind = 0; ind < str.len; ind++) {
        buf[ind] = str.ptr[ind];
    }
    return buf;
}

#endif  // aln_IMPLEMENTATION

#ifdef __clang__
#pragma clang diagnostic pop
#endif
