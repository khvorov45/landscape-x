#include "programmable_build.h"

typedef int32_t i32;

static int
strcmpSort(const void* a, const void* b) {
    const char* astr = ((prb_String*)a)->ptr;
    const char* bstr = ((prb_String*)b)->ptr;
    int         result = strcmp(astr, bstr);
    return result;
}

static void
printObjList(prb_Arena* arena, const char* title, prb_String* objs) {
    i32 objsLen = arrlen(objs);
    qsort(objs, objsLen, sizeof(*objs), strcmpSort);
    prb_writelnToStdout(prb_fmt(arena, "\n%s: %d", title, objsLen));
    for (i32 index = 0; index < objsLen; index++) {
        prb_writelnToStdout(objs[index]);
    }
}

static bool
notIn(prb_String val, prb_String* arr, i32 arrLen) {
    bool result = true;
    for (i32 arrIndex = 0; arrIndex < arrLen && result; arrIndex++) {
        prb_String arrVal = arr[arrIndex];
        if (prb_streq(val, arrVal)) {
            result = false;
        }
    }
    return result;
}

int
main() {
    prb_TimeStart scriptStart = prb_timeStart();
    prb_Arena     arena = prb_createArenaFromVmem(1 * prb_GIGABYTE);
    prb_String    rootDir = prb_getParentDir(&arena, prb_STR(__FILE__));
    prb_String    coreDir = prb_pathJoin(&arena, rootDir, prb_STR("core"));
    prb_String    outDir = prb_pathJoin(&arena, rootDir, prb_STR("build-debug"));
    prb_createDirIfNotExists(&arena, outDir);
    // prb_clearDirectory(&arena, outDir);

    // NOTE(sen) Objs
    prb_String srcFilesNoObj[] = {
        prb_STR("JTT.c"),
        prb_STR("blosum.c"),
        prb_STR("iteration.c")};
    prb_String* objsWithMain = 0;
    prb_String* objsWithoutMain = 0;
    prb_String  objDir = prb_pathJoin(&arena, outDir, prb_STR("objs"));
    prb_createDirIfNotExists(&arena, objDir);
    prb_PathFindIterator sourceIter = prb_createPathFindIter((prb_PathFindSpec) {.arena = &arena, .dir = coreDir, .mode = prb_PathFindMode_Glob, .glob.pattern = prb_STR("*.c")});
    prb_ProcessHandle*   objProccesses = 0;
    prb_Multitime        objMultitime = prb_createMultitime();
    while (prb_pathFindIterNext(&sourceIter)) {
        prb_String inname = prb_getLastEntryInPath(sourceIter.curPath);
        if (notIn(inname, srcFilesNoObj, prb_arrayLength(srcFilesNoObj))) {
            prb_String        outname = prb_replaceExt(&arena, inname, prb_STR("obj"));
            prb_String        outpath = prb_pathJoin(&arena, objDir, outname);
            prb_FileTimestamp lastModIn = prb_getLastModified(&arena, sourceIter.curPath);
            prb_assert(lastModIn.valid);
            prb_FileTimestamp lastModOut = prb_getLastModified(&arena, outpath);
            prb_multitimeAdd(&objMultitime, lastModOut);

            if (!lastModOut.valid || lastModOut.timestamp < lastModIn.timestamp) {
                prb_String cmd = prb_fmt(
                    &arena,
                    "clang -g -Denablemultithread -Werror %.*s -c -o %.*s",
                    prb_LIT(sourceIter.curPath),
                    prb_LIT(outpath)
                );
                prb_writelnToStdout(cmd);
                prb_ProcessHandle cmdProcess = prb_execCmd(&arena, cmd, prb_ProcessFlag_DontWait, (prb_String) {});
                prb_assert(cmdProcess.status == prb_ProcessStatus_Launched);
                arrput(objProccesses, cmdProcess);
            }

            prb_ReadEntireFileResult fileRead = prb_readEntireFile(&arena, sourceIter.curPath);
            prb_assert(fileRead.success);
            prb_String           fileContentStr = prb_strFromBytes(fileRead.content);
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

    // NOTE(sen) Executables
    bool printObjs = false;
    if (printObjs) {
        printObjList(&arena, "With main", objsWithMain);
        printObjList(&arena, "Without main", objsWithoutMain);
    }
    prb_ProcessHandle* exeProcesses = 0;
    prb_String         exeDir = prb_pathJoin(&arena, outDir, prb_STR("exes"));
    prb_createDirIfNotExists(&arena, exeDir);
    // prb_clearDirectory(&arena, exeDir);
    for (i32 withMainIndex = 0; withMainIndex < arrlen(objsWithMain); withMainIndex++) {
        prb_String objWithMain = objsWithMain[withMainIndex];
        prb_String objWithMainName = prb_getLastEntryInPath(objWithMain);
        if (!prb_streq(objWithMainName, prb_STR("interface.obj"))) {
            prb_String           exeName = objWithMainName;
            prb_StringFindResult dotFind = prb_strFind((prb_StringFindSpec) {
                .str = exeName,
                .pattern = prb_STR("."),
                .mode = prb_StringFindMode_AnyChar,
                .direction = prb_StringDirection_FromEnd,
            });
            prb_assert(dotFind.found);
            exeName.len = dotFind.matchByteIndex;
            prb_String exePath = prb_pathJoin(&arena, exeDir, exeName);

            bool shouldRecompile = true;
            if (objMultitime.invalidAddedTimestampsCount == 0) {
                prb_FileTimestamp exeLastMod = prb_getLastModified(&arena, exePath);
                if (exeLastMod.valid) {
                    shouldRecompile = exeLastMod.timestamp < objMultitime.timeLatest;
                }
            }

            if (shouldRecompile) {
                prb_GrowingString gstr = prb_beginString(&arena);
                prb_addStringSegment(&gstr, "clang %.*s", prb_LIT(objWithMain));
                for (i32 withoutMainIndex = 0; withoutMainIndex < arrlen(objsWithoutMain); withoutMainIndex++) {
                    prb_String objWithoutMain = objsWithoutMain[withoutMainIndex];
                    prb_String objWithoutMainName = prb_getLastEntryInPath(objWithoutMain);

                    bool addThisObj = true;
                    if (prb_streq(objWithoutMainName, prb_STR("pairlocalalign.obj"))) {
                        prb_String exesToNotCouple[] = {prb_STR("pairash.obj"), prb_STR("rnatest.obj"), prb_STR("multi2hat3s.obj")};
                        addThisObj = notIn(objWithMainName, exesToNotCouple, prb_arrayLength(exesToNotCouple));
                    }

                    if (addThisObj) {
                        prb_addStringSegment(&gstr, " %.*s", prb_LIT(objWithoutMain));
                    }
                }
                prb_addStringSegment(&gstr, " -o %.*s -lpthread -lm", prb_LIT(exePath));
                prb_String cmd = prb_endString(&gstr);

                prb_writelnToStdout(cmd);
                prb_ProcessHandle exeHandle = prb_execCmd(&arena, cmd, prb_ProcessFlag_DontWait, (prb_String) {});
                prb_assert(exeHandle.status == prb_ProcessStatus_Launched);
                arrput(exeProcesses, exeHandle);
            }
        }
    }
    prb_assert(prb_waitForProcesses(exeProcesses, arrlen(exeProcesses)) == prb_Success);

    // NOTE(sen) My mafft
    prb_String mymafftPath = prb_pathJoin(&arena, rootDir, prb_STR("rewrite/mafft_cli.c"));
    prb_String mymafftExe = prb_pathJoin(&arena, outDir, prb_STR("mafft_cli.bin"));
    prb_String mymafftCompile = prb_fmt(&arena, "clang -g -Wall -Wextra -Werror %.*s -o %.*s", prb_LIT(mymafftPath), prb_LIT(mymafftExe));
    prb_writelnToStdout(mymafftCompile);
    prb_assert(prb_execCmd(&arena, mymafftCompile, 0, (prb_String) {}).status == prb_ProcessStatus_CompletedSuccess);

    prb_writelnToStdout(prb_fmt(&arena, "total: %.2fms", prb_getMsFrom(scriptStart)));
    return 0;
}
