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
#include <stdalign.h>

#ifndef aln_PUBLICAPI
#define aln_PUBLICAPI
#endif

typedef enum aln_Status {
    aln_Failure,
    aln_Success,
} aln_Status;

typedef struct aln_Str {
    char*   ptr;
    int32_t len;
} aln_Str;

typedef struct aln_Opts {
    int32_t outputhat23;
} aln_Opts;

aln_PUBLICAPI aln_Opts aln_defaultOpts(void);

// TODO(sen) Remove private forward decls once we pull everything into 1 TU

//
// SECTION Private
//

#define aln_max(a, b) (((a) > (b)) ? (a) : (b))
#define aln_min(a, b) (((a) < (b)) ? (a) : (b))
#define aln_clamp(x, a, b) (((x) < (a)) ? (a) : (((x) > (b)) ? (b) : (x)))
#define aln_arrayCount(arr) (int32_t)(sizeof(arr) / sizeof(arr[0]))
#define aln_arenaAllocArray(arena, type, len) (type*)aln_arenaAllocAndZero(arena, (len) * (int32_t)sizeof(type), alignof(type))
#define aln_arenaAllocStruct(arena, type) (type*)aln_arenaAllocAndZero(arena, sizeof(type), alignof(type))
#define aln_isPowerOf2(x) (((x) > 0) && (((x) & ((x)-1)) == 0))
#define aln_unused(x) ((x) = (x))

// clang-format off
// Taken from portable snippets
// https://github.com/nemequ/portable-snippets/blob/master/debug-trap/debug-trap.h
#if defined(__has_builtin) && !defined(__ibmxl__)
#if __has_builtin(__builtin_debugtrap)
#define aln_debugbreak() __builtin_debugtrap()
#elif __has_builtin(__debugbreak)
#define aln_debugbreak() __debugbreak()
#endif
#endif
#if !defined(aln_debugbreak)
#if defined(_MSC_VER) || defined(__INTEL_COMPILER)
#define aln_debugbreak() __debugbreak()
#elif defined(__ARMCC_VERSION)
#define aln_debugbreak() __breakpoint(42)
#elif defined(__ibmxl__) || defined(__xlC__)
#include <builtins.h>
#define aln_debugbreak() __trap(42)
#elif defined(__DMC__) && defined(_M_IX86)
#define aln_debugbreak() __asm int 3h;
#elif defined(__i386__) || defined(__x86_64__)
#define aln_debugbreak() __asm__ __volatile__("int3")
#elif defined(__thumb__)
#define aln_debugbreak() __asm__ __volatile__(".inst 0xde01")
#elif defined(__aarch64__)
#define aln_debugbreak() __asm__ __volatile__(".inst 0xd4200000")
#elif defined(__arm__)
#define aln_debugbreak() __asm__ __volatile__(".inst 0xe7f001f0")
#elif defined (__alpha__) && !defined(__osf__)
#define aln_debugbreak() __asm__ __volatile__("bpt")
#elif defined(_54_)
#define aln_debugbreak() __asm__ __volatile__("ESTOP")
#elif defined(_55_)
#define aln_debugbreak() __asm__ __volatile__(";\n .if (.MNEMONIC)\n ESTOP_1\n .else\n ESTOP_1()\n .endif\n NOP")
#elif defined(_64P_)
#define aln_debugbreak() __asm__ __volatile__("SWBP 0")
#elif defined(_6x_)
#define aln_debugbreak() __asm__ __volatile__("NOP\n .word 0x10000000")
#elif defined(__STDC_HOSTED__) && (__STDC_HOSTED__ == 0) && defined(__GNUC__)
#define aln_debugbreak() __builtin_trap()
#else
#include <signal.h>
#if defined(SIGTRAP)
#define aln_debugbreak() raise(SIGTRAP)
#else
#define aln_debugbreak() raise(SIGABRT)
#endif
#endif
#endif

#ifndef aln_assertAction
#define aln_assertAction() do {aln_debugbreak();} while (0)
#endif

#ifndef aln_assert
#define aln_assert(condition) do { if (condition) {} else { aln_assertAction(); } } while (0)
#endif
// clang-format on

#ifndef aln_PRIVATEAPI
#define aln_PRIVATEAPI static
#endif

// TODO(sen) Don't assert on out of memory?

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

aln_PRIVATEAPI int32_t        aln_getOffsetForAlignment(void* ptr, int32_t align);
aln_PRIVATEAPI aln_Arena      aln_createArenaFromArena(aln_Arena* arena, intptr_t bytes);
aln_PRIVATEAPI void*          aln_arenaAllocAndZero(aln_Arena* arena, int32_t size, int32_t align);
aln_PRIVATEAPI void           aln_arenaAlignFreePtr(aln_Arena* arena, int32_t align);
aln_PRIVATEAPI void*          aln_arenaFreePtr(aln_Arena* arena);
aln_PRIVATEAPI intptr_t       aln_arenaFreeSize(aln_Arena* arena);
aln_PRIVATEAPI void           aln_arenaChangeUsed(aln_Arena* arena, intptr_t byteDelta);
aln_PRIVATEAPI aln_TempMemory aln_beginTempMemory(aln_Arena* arena);
aln_PRIVATEAPI void           aln_endTempMemory(aln_TempMemory temp);

#endif // aln_HEADER

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

aln_PUBLICAPI aln_Opts 
aln_defaultOpts(void) {
    aln_Opts opts = {
        .outputhat23 = 16,
    };
    return opts;
}

#endif  // aln_IMPLEMENTATION

#ifdef __clang__
#pragma clang diagnostic pop
#endif
