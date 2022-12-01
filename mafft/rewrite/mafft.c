#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>

#define STB_SPRINTF_IMPLEMENTATION
#include "stb_sprintf.h"

#if defined(WIN32) || defined(_WIN32)
#define PLATFORM_WINDOWS 1
#elif (defined(linux) || defined(__linux) || defined(__linux__))
#define PLATFORM_LINUX 1
#else
#error unrecognized platform
#endif

#if PLATFORM_WINDOWS

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <shellapi.h>
#pragma comment(lib, "Shell32.lib")

#elif PLATFORM_LINUX

#include <unistd.h>

#endif

#define BYTE 1
#define KILOBYTE 1024 * BYTE
#define MEGABYTE 1024 * KILOBYTE
#define GIGABYTE 1024 * MEGABYTE

#define unused(x) ((x) = (x))
#define function static
#define isPowerOf2(x) (((x) > 0) && (((x) & ((x)-1)) == 0))
#define arenaAllocArray(arena, type, len) (type*)arenaAllocAndZero(arena, (len) * sizeof(type), alignof(type))

// clang-format off
// Debug break taken from SDL
// https://github.com/libsdl-org/SDL/blob/main/include/SDL_assert.h
#ifdef __has_builtin
    #define _HAS_BUILTIN(x) __has_builtin(x)
#else
    #define _HAS_BUILTIN(x) 0
#endif

#if defined(_MSC_VER)
    /* Don't include intrin.h here because it contains C++ code */
    extern void __cdecl __debugbreak(void);
    #define debugbreak() __debugbreak()
#elif _HAS_BUILTIN(__builtin_debugtrap)
    #define debugbreak() __builtin_debugtrap()
#elif ((!defined(__NACL__)) && ((defined(__GNUC__) || defined(__clang__)) && (defined(__i386__) || defined(__x86_64__))))
    #define debugbreak() __asm__ __volatile__("int $3\n\t")
#elif (defined(__APPLE__) && (defined(__arm64__) || defined(__aarch64__))) /* this might work on other ARM targets, but this is a known quantity... */
    #define debugbreak() __asm__ __volatile__("brk #22\n\t")
#elif defined(__APPLE__) && defined(__arm__)
    #define debugbreak() __asm__ __volatile__("bkpt #22\n\t")
#elif defined(__386__) && defined(__WATCOMC__)
    #define debugbreak() { _asm { int 0x03 } }
#elif defined(HAVE_SIGNAL_H) && !defined(__WATCOMC__)
    #include <signal.h>
    #define debugbreak() raise(SIGTRAP)
#else
    /* How do we trigger breakpoints on this platform? */
    #define debugbreak()
#endif

#ifndef assertAction
#define assertAction() do {\
    debugbreak();\
    terminate(1);\
} while (0)
#endif

// NOTE(khvorov) Assign condition to a variable to catch assert(x = y) instead of assert(x == y)
#define assert(condition) do { bool assertbool = condition; if (!(assertbool)) { assertAction(); } } while (0)

#define STR(x) (String) {x, cstringLength(x)}
#define LIT(x) (x).len, (x).ptr
// clang-format on

typedef uint8_t   u8;
typedef int32_t   i32;
typedef uintptr_t usize;
typedef intptr_t  isize;

typedef enum Status {
    Failure,
    Success,
} Status;

function void
terminate(i32 code) {
#if PLATFORM_WINDOWS

#error unimplemented

#elif PLATFORM_LINUX

    _exit(code);

#else
#error unimplemented
#endif
}

//
// SECTION Memory
//

#if PLATFORM_LINUX
#include <sys/mman.h>
#endif

typedef struct Arena {
    void* base;
    i32   size;
    i32   used;
    i32   tempCount;
} Arena;

typedef struct TempMemory {
    Arena* arena;
    i32    usedAtBegin;
    i32    tempCountAtBegin;
} TempMemory;

function i32
getOffsetForAlignment(void* ptr, i32 align) {
    assert(isPowerOf2(align));
    usize ptrAligned = (usize)((uint8_t*)ptr + (align - 1)) & (usize)(~(align - 1));
    assert(ptrAligned >= (usize)ptr);
    isize diff = ptrAligned - (usize)ptr;
    assert(diff < align && diff >= 0);
    return (i32)diff;
}

function void*
arenaFreePtr(Arena* arena) {
    void* result = (uint8_t*)arena->base + arena->used;
    return result;
}

function i32
arenaFreeSize(Arena* arena) {
    i32 result = arena->size - arena->used;
    return result;
}

function void
arenaChangeUsed(Arena* arena, i32 byteDelta) {
    assert(arenaFreeSize(arena) >= byteDelta);
    arena->used += byteDelta;
}

function void
arenaAlignFreePtr(Arena* arena, i32 align) {
    i32 offset = getOffsetForAlignment(arenaFreePtr(arena), align);
    arenaChangeUsed(arena, offset);
}

function void*
arenaAllocAndZero(Arena* arena, i32 size, i32 align) {
    arenaAlignFreePtr(arena, align);
    void* result = arenaFreePtr(arena);
    arenaChangeUsed(arena, size);
    for (i32 index = 0; index < size; index++) {
        ((u8*)result)[index] = 0;
    }
    return result;
}

function TempMemory
beginTempMemory(Arena* arena) {
    TempMemory temp = {.arena = arena, .usedAtBegin = arena->used, .tempCountAtBegin = arena->tempCount};
    arena->tempCount += 1;
    return temp;
}

function void
endTempMemory(TempMemory temp) {
    assert(temp.arena->tempCount == temp.tempCountAtBegin + 1);
    temp.arena->used = temp.usedAtBegin;
    temp.arena->tempCount -= 1;
}

typedef struct VmemAllocResult {
    bool  success;
    void* ptr;
} VmemAllocResult;

function VmemAllocResult
vmemAlloc(i32 bytes) {
    VmemAllocResult result = {};

#if PLATFORM_WINDOWS

    void* ptr = VirtualAlloc(0, bytes, MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE);

#elif PLATFORM_LINUX

    result.ptr = mmap(0, bytes, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    result.success = result.ptr != MAP_FAILED;

#else
#error unimplemented
#endif

    return result;
}

typedef struct CreateArenaFromVmemResult {
    bool  success;
    Arena arena;
} CreateArenaFromVmemResult;

function CreateArenaFromVmemResult
createArenaFromVmem(i32 bytes) {
    CreateArenaFromVmemResult result = {};
    VmemAllocResult           vmem = vmemAlloc(bytes);
    if (vmem.success) {
        result.success = true;
        result.arena.base = vmem.ptr;
        result.arena.size = bytes;
    }
    return result;
}

//
// SECTION Strings
//

typedef struct String {
    const char* ptr;
    i32         len;
} String;

function i32
cstringLength(const char* cstring) {
    i32 len = 0;
    while (cstring[len]) {
        len += 1;
    }
    return len;
}

function String
strSliceForward(String str, int32_t bytes) {
    assert(bytes <= str.len);
    String result = {str.ptr + bytes, str.len - bytes};
    return result;
}

function String
vfmtCustomBuffer(void* buf, i32 bufSize, const char* fmt, va_list args) {
    i32    len = stbsp_vsnprintf((char*)buf, bufSize, fmt, args);
    String result = {(const char*)buf, len};
    return result;
}

function String
fmt(Arena* arena, const char* fmt, ...) {
    va_list args;
    va_start(args, fmt);
    String result = vfmtCustomBuffer(arenaFreePtr(arena), arenaFreeSize(arena), fmt, args);
    arenaChangeUsed(arena, result.len + 1);  // NOTE(khvorov) Null terminator
    va_end(args);
    return result;
}

function const char*
strGetNullTerminated(Arena* arena, String str) {
    const char* result = fmt(arena, "%.*s", LIT(str)).ptr;
    return result;
}

function Status
writeToStdout(String str) {
    Status status = Failure;

#if PLATFORM_WINDOWS

    HANDLE out = GetStdHandle(STD_OUTPUT_HANDLE);
    WriteFile(out, msg.ptr, msg.len, 0, 0);

#elif PLATFORM_LINUX

    ssize_t writeResult = write(STDOUT_FILENO, str.ptr, str.len);
    status = writeResult == str.len ? Success : Failure;

#else
#error unimplemented
#endif

    return status;
}

function void
writelnToStdout(String str) {
    writeToStdout(str);
    writeToStdout(STR("\n"));
}

typedef struct StringFindResult {
    bool found;
    i32  matchByteIndex;
} StringFindResult;

function StringFindResult
strFindAny(String input, String chars) {
    StringFindResult result = {};

    if (chars.len > 0) {
        for (i32 inputIndex = 0; inputIndex < input.len && !result.found; inputIndex++) {
            for (i32 charsIndex = 0; charsIndex < chars.len && !result.found; charsIndex++) {
                if (input.ptr[inputIndex] == chars.ptr[charsIndex]) {
                    result.found = true;
                    result.matchByteIndex = inputIndex;
                }
            }
        }
    }

    return result;
}

typedef struct LineIterator {
    String ogstr;
    i32    curLineCount;
    i32    curByteOffset;
    String curLine;
    i32    curLineEndLen;
} LineIterator;

function LineIterator
createLineIter(String str) {
    LineIterator iter = {
        .ogstr = str,
        .curLineCount = 0,
        .curByteOffset = 0,
        .curLine = (String) {},
        .curLineEndLen = 0,
    };
    return iter;
}

function Status
lineIterNext(LineIterator* iter) {
    Status result = Failure;

    iter->curByteOffset += iter->curLine.len + iter->curLineEndLen;
    iter->curLine = (String) {};
    iter->curLineEndLen = 0;

    if (iter->curByteOffset < iter->ogstr.len) {
        iter->curLine = strSliceForward(iter->ogstr, iter->curByteOffset);
        StringFindResult lineEndResult = strFindAny(iter->curLine, STR("\r\n"));
        if (lineEndResult.found) {
            iter->curLine.len = lineEndResult.matchByteIndex;
            iter->curLineEndLen = 1;
            if (iter->curLine.ptr[iter->curLine.len] == '\r'
                && iter->curByteOffset + iter->curLine.len + 1 < iter->ogstr.len
                && iter->curLine.ptr[iter->curLine.len + 1] == '\n') {
                iter->curLineEndLen += 1;
            }
        }
        iter->curLineCount += 1;
        result = Success;
    }

    return result;
}

//
// SECTION Filesystem
//

#if PLATFORM_LINUX

#include <sys/stat.h>
#include <fcntl.h>

typedef struct linux_GetFileStatResult {
    bool        success;
    struct stat stat;
} linux_GetFileStatResult;

typedef struct linux_OpenResult {
    bool success;
    int  handle;
} linux_OpenResult;

function linux_OpenResult
linux_open(Arena* arena, String path, int oflags, mode_t mode) {
    TempMemory       temp = beginTempMemory(arena);
    linux_OpenResult result = {};
    const char*      pathNull = strGetNullTerminated(arena, path);
    result.handle = open(pathNull, oflags, mode);
    result.success = result.handle != -1;
    endTempMemory(temp);
    return result;
}

#endif

typedef struct Bytes {
    uint8_t* data;
    i32      len;
} Bytes;

typedef struct ReadEntireFileResult {
    bool  success;
    Bytes bytes;
} ReadEntireFileResult;

function ReadEntireFileResult
readEntireFile(Arena* arena, String path) {
    ReadEntireFileResult result = {};

#if PLATFORM_WINDOWS

#error unimplemented

#elif PLATFORM_LINUX

    linux_OpenResult handleResult = linux_open(arena, path, O_RDONLY, 0);
    if (handleResult.success) {
        int         handle = handleResult.handle;
        struct stat statBuf = {};
        if (fstat(handle, &statBuf) == 0) {
            uint8_t* buf = (uint8_t*)arenaAllocAndZero(arena, statBuf.st_size + 1, 1);  // NOTE(sen) Null terminator just in case
            i32      readResult = read(handle, buf, statBuf.st_size);
            if (readResult == statBuf.st_size) {
                result.success = true;
                result.bytes = (Bytes) {buf, readResult};
            }
            close(handle);
        }
    }

#else
#error unimplemented
#endif

    return result;
}

//
// SECTION Main
//

typedef struct FastaEntry {
    String name;
    String seq;
} FastaEntry;

typedef struct Fasta {
    FastaEntry* entry;
    i32         entryCount;
} Fasta;

typedef struct ReadFastaEntryResult {
    bool       success;
    i32        bytesRead;
    FastaEntry entry;
} ReadFastaEntryResult;

function ReadFastaEntryResult
readFastaEntry(String input) {
    ReadFastaEntryResult result = {};

    // NOTE(sen) Skip empty lines
    LineIterator lineIter = createLineIter(input);
    while (lineIterNext(&lineIter) == Success && lineIter.curLine.len == 0) {}

    if (lineIter.curLine.len > 0 && lineIter.curLine.ptr[0] == '>') {
        result.entry.name = lineIter.curLine;
        while (lineIterNext(&lineIter) == Success) {
            writelnToStdout(lineIter.curLine);
        }
    }
    return result;
}

typedef struct ReadFastaResult {
    bool  success;
    Fasta fasta;
} ReadFastaResult;

function ReadFastaResult
readFasta(Bytes fileContent) {
    ReadFastaResult result = {.success = true};
    String          str = (String) {(const char*)fileContent.data, fileContent.len};
    while (str.len > 0 && result.success) {
        ReadFastaEntryResult entryRead = readFastaEntry(str);
        if (entryRead.success) {
            str = strSliceForward(str, entryRead.bytesRead);
        } else {
            result.success = false;
        }
    }
    return result;
}

int
main() {
    CreateArenaFromVmemResult arenaResult = createArenaFromVmem(1 * GIGABYTE);
    if (arenaResult.success) {
        Arena*               arena = &arenaResult.arena;
        String               inputpath = STR("mafft/temp.fs");
        ReadEntireFileResult readResult = readEntireFile(arena, inputpath);
        if (readResult.success) {
            Bytes           inputContent = readResult.bytes;
            ReadFastaResult inputFastaResult = readFasta(inputContent);
            if (inputFastaResult.success) {
                Fasta inputFasta = inputFasta;
                unused(inputFasta);
            } else {
                writelnToStdout(STR("Failed to parse input file"));
            }
        } else {
            writelnToStdout(STR("Failed to read input file"));
        }
    } else {
        writelnToStdout(STR("Failed to allocate memory"));
    }
    return 0;
}
