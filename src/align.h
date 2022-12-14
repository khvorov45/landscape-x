#include <stdint.h>

typedef enum aln_Status {
    aln_Failure,
    aln_Success,
} aln_Status;

typedef struct aln_String {
    char*   ptr;
    int32_t len;
} aln_String;

aln_Status aln_alignStrings(aln_String* strings, int32_t stringsCount, void* out, int32_t outBytes);

#ifdef aln_IMPLEMENTATION

aln_Status
aln_alignStrings(aln_String* strings, int32_t stringsCount, void* out, int32_t outBytes) {
    aln_Status result = aln_Failure;
    return result;
}

#endif  // aln_IMPLEMENTATION
