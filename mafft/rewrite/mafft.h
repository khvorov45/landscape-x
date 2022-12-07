#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdalign.h>

#ifndef mafft_PUBLICDEC
#define mafft_PUBLICDEC
#endif

#ifndef mafft_PUBLICDEF
#define mafft_PUBLICDEF
#endif

#ifndef static
#define static
#endif

#ifndef mafft_assert
#define mafft_assert(condition)
#endif

typedef struct mafft_Arena {
    void*   base;
    int32_t size;
    int32_t used;
    int32_t tempCount;
    bool    lockedForString;
} mafft_Arena;

typedef struct mafft_String {
    const char* ptr;
    int32_t     len;
} mafft_String;

typedef struct mafft_FastaEntry {
    mafft_String name;
    mafft_String seq;
} mafft_FastaEntry;

typedef struct mafft_Fasta {
    mafft_FastaEntry* entries;
    int32_t           entryCount;
} mafft_Fasta;

typedef struct mafft_ReadFastaResult {
    bool        success;
    mafft_Fasta fasta;
} mafft_ReadFastaResult;

mafft_PUBLICDEC void                  mafft_initArena(mafft_Arena* arena, void* ptr, int32_t len);
mafft_PUBLICDEC mafft_ReadFastaResult mafft_readFasta(mafft_Arena* arena, const void* ptr, int32_t len);
mafft_PUBLICDEC void                  mafft_alignSeq(mafft_String* seqs, int32_t seqCount, void* output);

#ifdef mafft_IMPLEMENTATION

#define mafft_BYTE 1
#define mafft_KILOBYTE 1024 * mafft_BYTE
#define mafft_MEGABYTE 1024 * mafft_KILOBYTE
#define mafft_GIGABYTE 1024 * mafft_MEGABYTE

#define mafft_max(a, b) (((a) > (b)) ? (a) : (b))
#define mafft_min(a, b) (((a) < (b)) ? (a) : (b))
#define mafft_clamp(x, a, b) (((x) < (a)) ? (a) : (((x) > (b)) ? (b) : (x)))
#define mafft_arrayLength(arr) (int32_t)(sizeof(arr) / sizeof(arr[0]))
#define mafft_unused(x) ((x) = (x))
#define mafft_isPowerOf2(x) (((x) > 0) && (((x) & ((x)-1)) == 0))

// clang-format off
#define mafft_STR(x) (mafft_String){x, mafft_cstringLength(x)}
// clang-format on

static int32_t
mafft_getOffsetForAlignment(void* ptr, int32_t align) {
    mafft_assert(mafft_isPowerOf2(align));
    uintptr_t ptrAligned = (uintptr_t)((uint8_t*)ptr + (align - 1)) & (uintptr_t)(~(align - 1));
    mafft_assert(ptrAligned >= (uintptr_t)ptr);
    intptr_t diff = ptrAligned - (uintptr_t)ptr;
    mafft_assert(diff < align && diff >= 0);
    return (int32_t)diff;
}

static void*
mafft_arenaFreePtr(mafft_Arena* arena) {
    void* result = (uint8_t*)arena->base + arena->used;
    return result;
}

static int32_t
mafft_arenaFreeSize(mafft_Arena* arena) {
    int32_t result = arena->size - arena->used;
    return result;
}

typedef enum mafft_Status {
    mafft_Failure,
    mafft_Success,
} mafft_Status;

#ifdef __GNUC__
#define mafft_MUST_USE __attribute__((warn_unused_result))
#elif
#define mafft_MUST_USE
#endif

mafft_MUST_USE static mafft_Status
mafft_arenaChangeUsed(mafft_Arena* arena, int32_t byteDelta) {
    bool canChange = byteDelta == 0
        || (byteDelta > 0 && mafft_arenaFreeSize(arena) >= byteDelta)
        || (byteDelta < 0 && arena->used >= -byteDelta);
    mafft_Status result = mafft_Failure;
    if (canChange) {
        arena->used += byteDelta;
        result = mafft_Success;
    }
    return result;
}

mafft_MUST_USE static mafft_Status
mafft_arenaAlignFreePtr(mafft_Arena* arena, int32_t align) {
    int32_t      offset = mafft_getOffsetForAlignment(mafft_arenaFreePtr(arena), align);
    mafft_Status result = mafft_arenaChangeUsed(arena, offset);
    return result;
}

typedef struct mafft_ArenaAllocAndZeroResult {
    bool  success;
    void* ptr;
} mafft_ArenaAllocAndZeroResult;

mafft_MUST_USE static mafft_ArenaAllocAndZeroResult
mafft_arenaAllocAndZero(mafft_Arena* arena, int32_t size, int32_t align) {
    mafft_ArenaAllocAndZeroResult result = {};
    if (mafft_arenaAlignFreePtr(arena, align) == mafft_Success) {
        void* tentativeResult = mafft_arenaFreePtr(arena);
        if (mafft_arenaChangeUsed(arena, size) == mafft_Success) {
            result.ptr = tentativeResult;
            result.success = true;
            for (int32_t index = 0; index < size; index++) {
                ((uint8_t*)result.ptr)[index] = 0;
            }
        }
    }
    return result;
}

typedef struct mafft_TempMemory {
    mafft_Arena* arena;
    int32_t      usedAtBegin;
    int32_t      tempCountAtBegin;
} mafft_TempMemory;

static mafft_TempMemory
mafft_beginTempMemory(mafft_Arena* arena) {
    mafft_TempMemory temp = {.arena = arena, .usedAtBegin = arena->used, .tempCountAtBegin = arena->tempCount};
    arena->tempCount += 1;
    return temp;
}

static void
mafft_endTempMemory(mafft_TempMemory temp) {
    mafft_assert(temp.arena->tempCount == temp.tempCountAtBegin + 1);
    temp.arena->used = temp.usedAtBegin;
    temp.arena->tempCount -= 1;
}

static void
mafft_keepTempMemory(mafft_TempMemory temp) {
    mafft_assert(temp.arena->tempCount == temp.tempCountAtBegin + 1);
    temp.arena->tempCount -= 1;
}

static int32_t
mafft_cstringLength(const char* cstring) {
    int32_t len = 0;
    while (cstring[len]) {
        len += 1;
    }
    return len;
}

static mafft_String
mafft_strSliceForward(mafft_String str, int32_t bytes) {
    mafft_assert(bytes <= str.len);
    mafft_String result = {str.ptr + bytes, str.len - bytes};
    return result;
}

typedef struct mafft_StringFindResult {
    bool    found;
    int32_t matchByteIndex;
} mafft_StringFindResult;

static mafft_StringFindResult
mafft_strFindAny(mafft_String input, mafft_String chars) {
    mafft_StringFindResult result = {};

    if (chars.len > 0) {
        for (int32_t inputIndex = 0; inputIndex < input.len && !result.found; inputIndex++) {
            for (int32_t charsIndex = 0; charsIndex < chars.len && !result.found; charsIndex++) {
                if (input.ptr[inputIndex] == chars.ptr[charsIndex]) {
                    result.found = true;
                    result.matchByteIndex = inputIndex;
                }
            }
        }
    }

    return result;
}

typedef struct mafft_LineIterator {
    mafft_String ogstr;
    int32_t      curLineCount;
    int32_t      curByteOffset;
    mafft_String curLine;
    int32_t      curLineEndLen;
} mafft_LineIterator;

static mafft_LineIterator
mafft_createLineIter(mafft_String str) {
    mafft_LineIterator iter = {
        .ogstr = str,
        .curLineCount = 0,
        .curByteOffset = 0,
        .curLine = (mafft_String) {},
        .curLineEndLen = 0,
    };
    return iter;
}

static mafft_Status
mafft_lineIterNext(mafft_LineIterator* iter) {
    mafft_Status result = mafft_Failure;

    iter->curByteOffset += iter->curLine.len + iter->curLineEndLen;
    iter->curLine = (mafft_String) {};
    iter->curLineEndLen = 0;

    if (iter->curByteOffset < iter->ogstr.len) {
        iter->curLine = mafft_strSliceForward(iter->ogstr, iter->curByteOffset);
        mafft_StringFindResult lineEndResult = mafft_strFindAny(iter->curLine, mafft_STR("\r\n"));
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
        result = mafft_Success;
    }

    return result;
}

typedef struct mafft_GrowingString {
    mafft_Arena* arena;
    mafft_String string;
} mafft_GrowingString;

static mafft_GrowingString
mafft_beginString(mafft_Arena* arena) {
    mafft_assert(!arena->lockedForString);
    arena->lockedForString = true;
    mafft_String        str = {(const char*)mafft_arenaFreePtr(arena), 0};
    mafft_GrowingString result = {arena, str};
    return result;
}

mafft_MUST_USE static mafft_Status
mafft_addStringSegment(mafft_GrowingString* gstr, mafft_String seg) {
    mafft_assert(gstr->arena->lockedForString);
    mafft_Status result = mafft_arenaChangeUsed(gstr->arena, seg.len);
    if (result == mafft_Success) {
        for (int32_t segIndex = 0; segIndex < seg.len; segIndex++) {
            ((char*)gstr->string.ptr)[gstr->string.len] = seg.ptr[segIndex];
            gstr->string.len += 1;
        }
    }
    return result;
}

typedef struct mafft_EndStringResult {
    bool         success;
    mafft_String string;
} mafft_EndStringResult;

mafft_MUST_USE static mafft_EndStringResult
mafft_endString(mafft_GrowingString* gstr) {
    mafft_assert(gstr->arena->lockedForString);
    gstr->arena->lockedForString = false;
    mafft_EndStringResult result = {};
    // NOTE(khvorov) Null terminator
    mafft_ArenaAllocAndZeroResult res = mafft_arenaAllocAndZero(gstr->arena, 1, 1);
    if (res.success) {
        result.string = gstr->string;
        result.success = true;
    }
    *gstr = (mafft_GrowingString) {};
    return result;
}

mafft_MUST_USE static mafft_EndStringResult
mafft_strCopy(mafft_Arena* arena, mafft_String str) {
    mafft_GrowingString   gstr = mafft_beginString(arena);
    mafft_EndStringResult result = {};
    if (mafft_addStringSegment(&gstr, str) == prb_Success) {
        result = mafft_endString(&gstr);
    }
    return result;
}

typedef struct mafft_ReadFastaEntryResult {
    bool             success;
    int32_t          bytesRead;
    mafft_FastaEntry entry;
} mafft_ReadFastaEntryResult;

static void
mafft_keepIfTrue(mafft_TempMemory temp, bool success) {
    if (success) {
        mafft_keepTempMemory(temp);
    } else {
        mafft_endTempMemory(temp);
    }
}

static mafft_ReadFastaEntryResult
mafft_readFastaEntry(mafft_Arena* arena, mafft_String input) {
    mafft_ReadFastaEntryResult result = {};
    mafft_TempMemory           temp = mafft_beginTempMemory(arena);

    mafft_LineIterator lineIter = mafft_createLineIter(input);

    // NOTE(sen) Skip empty lines
    while (mafft_lineIterNext(&lineIter) == mafft_Success && lineIter.curLine.len == 0) {}

    if (lineIter.curLine.len > 0 && lineIter.curLine.ptr[0] == '>') {
        result.success = true;
        mafft_String          name = mafft_strSliceForward(lineIter.curLine, 1);
        mafft_EndStringResult copyRes = mafft_strCopy(arena, name);
        if (copyRes.success) {
            result.entry.name = copyRes.string;
            mafft_GrowingString gstr = mafft_beginString(arena);
            bool                addingSuccessfully = true;
            while (mafft_lineIterNext(&lineIter) == mafft_Success && addingSuccessfully) {
                if (lineIter.curLine.len > 0 && lineIter.curLine.ptr[0] == '>') {
                    break;
                } else {
                    addingSuccessfully = mafft_addStringSegment(&gstr, lineIter.curLine) == prb_Success;
                }
            }
            if (addingSuccessfully) {
                mafft_EndStringResult endRes = mafft_endString(&gstr);
                if (endRes.success) {
                    result.success = true;
                    result.entry.seq = endRes.string;
                }
            }
        }
    }

    mafft_keepIfTrue(temp, result.success);
    result.bytesRead = lineIter.curByteOffset;
    return result;
}

mafft_PUBLICDEF void
mafft_initArena(mafft_Arena* arena, void* ptr, int32_t len) {
    *arena = (mafft_Arena) {};
    arena->base = ptr;
    arena->size = len;
}

mafft_PUBLICDEF mafft_ReadFastaResult
mafft_readFasta(mafft_Arena* arena, const void* ptr, int32_t len) {
    mafft_TempMemory      temp = mafft_beginTempMemory(arena);
    mafft_ReadFastaResult result = {.success = true};
    mafft_String          str = {(const char*)ptr, len};

    int32_t entryCount = 0;
    for (mafft_String strLeft = str;;) {
        mafft_StringFindResult headFind = mafft_strFindAny(strLeft, mafft_STR(">"));
        if (headFind.found) {
            strLeft = mafft_strSliceForward(strLeft, headFind.matchByteIndex + 1);
            entryCount += 1;
        } else {
            break;
        }
    }

    if (entryCount > 0) {
        mafft_ArenaAllocAndZeroResult arrAllocRes = mafft_arenaAllocAndZero(arena, sizeof(mafft_FastaEntry) * entryCount, alignof(mafft_FastaEntry));
        if (arrAllocRes.success) {
            result.fasta.entries = (mafft_FastaEntry*)arrAllocRes.ptr;
            for (mafft_String strLeft = str; strLeft.len > 0 && result.success;) {
                mafft_ReadFastaEntryResult entryRead = mafft_readFastaEntry(arena, strLeft);
                if (entryRead.success) {
                    strLeft = mafft_strSliceForward(strLeft, entryRead.bytesRead);
                    mafft_assert(result.fasta.entryCount < entryCount);
                    result.fasta.entries[result.fasta.entryCount] = entryRead.entry;
                    result.fasta.entryCount += 1;
                } else {
                    result.success = false;
                }
            }
        }
    }

    mafft_keepIfTrue(temp, result.success);
    return result;
}

mafft_PUBLICDEF int32_t
mafft_getAlignOutputBufferSize(mafft_String* seqs, int32_t seqCount) {
    int32_t maxLen = 0;
    for (int32_t seqIndex = 0; seqIndex < seqCount; seqIndex++) {
        mafft_String seq = seqs[seqIndex];
        maxLen = mafft_max(maxLen, seq.len);
    }
    // NOTE(sen) +1 for the null terminator
    int32_t size = (maxLen + 1) * seqCount;
    return size;
}

mafft_PUBLICDEF void
mafft_alignSeq(mafft_String* seqs, int32_t seqCount, void* output) {
    
}

#endif  // mafft_IMPLEMENTATION
