#include "programmable_build.h"

static int
strcmpSort(const void* a, const void* b) {
    const char* astr = ((prb_String*)a)->ptr;
    const char* bstr = ((prb_String*)b)->ptr;
    int         result = strcmp(astr, bstr);
    return result;
}

static void
printObjList(prb_Arena* arena, const char* title, prb_String* objs) {
    int32_t objsLen = arrlen(objs);
    qsort(objs, objsLen, sizeof(*objs), strcmpSort);
    prb_fmtAndPrintln(arena, "\n%s: %d", title, objsLen);
    for (int32_t index = 0; index < objsLen; index++) {
        prb_writelnToStdout(objs[index]);
    }
}

int
main() {
    prb_Arena  arena = prb_createArenaFromVmem(1 * prb_GIGABYTE);
    prb_String rootDir = prb_getParentDir(&arena, prb_STR(__FILE__));
    prb_String coreDir = prb_pathJoin(&arena, rootDir, prb_STR("core"));
    prb_String outDir = prb_pathJoin(&arena, rootDir, prb_STR("build-debug"));
    prb_createDirIfNotExists(&arena, outDir);
    prb_clearDirectory(&arena, outDir);

    // NOTE(sen) Objs
    prb_String* objsWithMain = 0;
    prb_String* objsWithoutMain = 0;
    prb_String  objDir = prb_pathJoin(&arena, outDir, prb_STR("objs"));
    prb_createDirIfNotExists(&arena, objDir);
    prb_PathFindIterator sourceIter = prb_createPathFindIter((prb_PathFindSpec) {.arena = &arena, .dir = coreDir, .mode = prb_PathFindMode_Glob, .glob.pattern = prb_STR("*.c")});
    prb_ProcessHandle*   objProccesses = 0;
    while (prb_pathFindIterNext(&sourceIter)) {
        prb_String inname = prb_getLastEntryInPath(sourceIter.curPath);
        if (!prb_streq(inname, prb_STR("JTT.c")) && !prb_streq(inname, prb_STR("blosum.c"))) {
            prb_String        outname = prb_replaceExt(&arena, inname, prb_STR("obj"));
            prb_String        outpath = prb_pathJoin(&arena, objDir, outname);
            prb_LastModResult lastModIn = prb_getLastModifiedFromPath(&arena, sourceIter.curPath);
            prb_assert(lastModIn.success);
            prb_LastModResult lastModOut = prb_getLastModifiedFromPath(&arena, outpath);
            if (!lastModOut.success || lastModOut.timestamp < lastModIn.timestamp) {
                prb_String cmd = prb_fmtAndPrintln(
                    &arena,
                    "clang -Denablemultithread -Werror %.*s -c -o %.*s",
                    prb_LIT(sourceIter.curPath),
                    prb_LIT(outpath)
                );
                prb_ProcessHandle cmdProcess = prb_execCmd(&arena, cmd, prb_ProcessFlag_DontWait, (prb_String) {});
                prb_assert(cmdProcess.status == prb_ProcessStatus_Launched);
                arrput(objProccesses, cmdProcess);
            }

            prb_Bytes            fileContent = prb_readEntireFile(&arena, sourceIter.curPath);
            prb_String           fileContentStr = (prb_String) {(const char*)fileContent.data, fileContent.len};
            prb_StringFindResult mainFound = prb_strFind((prb_StringFindSpec) {.str = fileContentStr, .pattern = prb_STR("main("), .mode = prb_StringFindMode_Exact});
            if (mainFound.found) {
                arrput(objsWithMain, outpath);
            } else {
                arrput(objsWithoutMain, outpath);
            }
        }
    }
    prb_destroyPathFindIter(&sourceIter);
    prb_assert(prb_waitForProcesses(objProccesses, arrlen(objProccesses)) == prb_Success);

    printObjList(&arena, "With main", objsWithMain);
    printObjList(&arena, "Without main", objsWithoutMain);

    return 0;
}
