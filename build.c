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

static void
copyFile(prb_Arena* arena, prb_String path, prb_String dir) {
    prb_String               copyPath = prb_pathJoin(arena, dir, prb_getLastEntryInPath(path));
    prb_ReadEntireFileResult fileContent = prb_readEntireFile(arena, path);
    prb_assert(fileContent.success);
    prb_assert(prb_writeEntireFile(arena, copyPath, fileContent.content.data, fileContent.content.len) == prb_Success);
}

int
main() {
    prb_TimeStart scriptStart = prb_timeStart();
    prb_Arena     arena_ = prb_createArenaFromVmem(1 * prb_GIGABYTE);
    prb_Arena*    arena = &arena_;
    prb_String    rootDir = prb_getParentDir(arena, prb_STR(__FILE__));
    prb_String    srcDir = prb_pathJoin(arena, rootDir, prb_STR("src"));
    prb_String    outDir = prb_pathJoin(arena, rootDir, prb_STR("build-debug"));
    prb_assert(prb_createDirIfNotExists(arena, outDir) == prb_Success);
    // prb_assert(prb_clearDirectory(arena, outDir) == prb_Success);

    // NOTE(sen) Get mafft
    prb_String mafftDir = prb_pathJoin(arena, rootDir, prb_STR("mafft"));
    if (!prb_isDirectory(arena, mafftDir) || prb_directoryIsEmpty(arena, mafftDir)) {
        prb_String cmd = prb_fmt(arena, "git clone https://gitlab.com/sysimm/mafft %.*s", prb_LIT(mafftDir));
        prb_writelnToStdout(cmd);
        prb_execCmd(arena, cmd, 0, (prb_String) {});
    }

    // NOTE(sen) Compile mafft
    prb_String coreDir = prb_pathJoin(arena, mafftDir, prb_STR("core"));
    prb_String mafftSuppressedWarnings[] = {
        prb_STR("-Wno-deprecated-non-prototype"),
        prb_STR("-Wno-unknown-directives"),
        prb_STR("-Wno-fortify-source"),
        prb_STR("-Wno-parentheses"),
        prb_STR("-Wno-format"),
        prb_STR("-Wno-incompatible-pointer-types"),
    };
    prb_String mafftSuppressedWarningsStr = prb_stringsJoin(arena, mafftSuppressedWarnings, prb_arrayLength(mafftSuppressedWarnings), prb_STR(" "));

    prb_String logfilePath = prb_pathJoin(arena, coreDir, prb_STR("logfile.txt"));
    prb_removeFileIfExists(arena, logfilePath);

    // NOTE(sen) Objs
    prb_String srcFilesNoObj[] = {
        prb_STR("JTT.c"),
        prb_STR("blosum.c"),
        prb_STR("iteration.c"),
    };
    prb_String* objsWithMain = 0;
    prb_String* objsWithoutMain = 0;
    prb_String  mafftOutDir = prb_pathJoin(arena, outDir, prb_STR("mafft"));
    prb_assert(prb_createDirIfNotExists(arena, mafftOutDir) == prb_Success);
    prb_String mafftObjDir = prb_pathJoin(arena, mafftOutDir, prb_STR("objs"));
    prb_assert(prb_createDirIfNotExists(arena, mafftObjDir) == prb_Success);
    // prb_assert(prb_clearDirectory(arena, mafftObjDir) == prb_Success);
    prb_ProcessHandle*   objProccesses = 0;
    prb_PathFindSpec     srcIterSpec = {.arena = arena, .dir = coreDir, .mode = prb_PathFindMode_Glob, .pattern = prb_STR("*.c")};
    prb_PathFindIterator sourceIter = prb_createPathFindIter(srcIterSpec);
    bool                 anyObjRecompiled = false;
    while (prb_pathFindIterNext(&sourceIter)) {
        prb_String inname = prb_getLastEntryInPath(sourceIter.curPath);
        if (notIn(inname, srcFilesNoObj, prb_arrayLength(srcFilesNoObj))) {
            prb_String outname = prb_replaceExt(arena, inname, prb_STR("obj"));
            prb_String outpath = prb_pathJoin(arena, mafftObjDir, outname);

            prb_ReadEntireFileResult fileRead = prb_readEntireFile(arena, sourceIter.curPath);
            prb_assert(fileRead.success);
            prb_String           fileContentStr = prb_strFromBytes(fileRead.content);
            prb_StringFindResult mainFound = prb_strFind((prb_StringFindSpec) {.str = fileContentStr, .pattern = prb_STR("main("), .mode = prb_StringFindMode_Exact});
            prb_String           srcPath = sourceIter.curPath;
            if (mainFound.found) {
                arrput(objsWithMain, outpath);

                bool       mainTakesArgs = false;
                prb_String mainArgcName = {};
                prb_String mainArgvName = {};
                {
                    prb_String       fromMain = prb_strSliceForward(fileContentStr, mainFound.matchByteIndex);
                    prb_LineIterator iter = prb_createLineIter(fromMain);
                    prb_assert(prb_lineIterNext(&iter) == prb_Success);
                    prb_StringFindResult intFound = prb_strFind((prb_StringFindSpec) {.str = iter.curLine, .pattern = prb_STR("int"), .mode = prb_StringFindMode_Exact});
                    mainTakesArgs = intFound.found;
                    if (intFound.found) {
                        prb_String           postInt = prb_strSliceForward(iter.curLine, intFound.matchByteIndex + intFound.matchLen);
                        prb_StringFindResult commaFound = prb_strFind((prb_StringFindSpec) {.str = postInt, .pattern = prb_STR(","), .mode = prb_StringFindMode_AnyChar});
                        prb_assert(commaFound.found);
                        mainArgcName = postInt;
                        mainArgcName.len = commaFound.matchByteIndex;
                        mainArgcName = prb_strTrim(mainArgcName);

                        bool insideAlpha = false;
                        i32  lastAlpha = 0;
                        i32  firstAlpha = 0;
                        for (i32 lineIndex = iter.curLine.len - 1; lineIndex >= 0; lineIndex--) {
                            char ch = iter.curLine.ptr[lineIndex];
                            bool isAlpha = (ch >= 'a' && ch <= 'z') || ((ch >= 'A' && ch <= 'Z'));
                            if (isAlpha && !insideAlpha) {
                                lastAlpha = lineIndex;
                                insideAlpha = true;
                            } else if (!isAlpha && insideAlpha) {
                                firstAlpha = lineIndex + 1;
                                break;
                            }
                        }

                        mainArgvName = prb_strSliceBetween(iter.curLine, firstAlpha, lastAlpha + 1);
                    }
                }

                if (mainTakesArgs) {
                    // NOTE(sen) Generate alternative file where main is modified
                    i32 mainLineStart = 0;
                    {
                        prb_String toMain = fileContentStr;
                        toMain.len = mainFound.matchByteIndex;
                        prb_StringFindResult res = prb_strFind((prb_StringFindSpec) {.str = toMain, .pattern = prb_STR("\r\n"), .direction = prb_StringDirection_FromEnd, .mode = prb_StringFindMode_AnyChar});
                        prb_assert(res.found);
                        mainLineStart = res.matchByteIndex;
                    }

                    i32 mainLineEnd = 0;
                    {
                        prb_String           fromMain = prb_strSliceForward(fileContentStr, mainFound.matchByteIndex);
                        prb_StringFindResult res = prb_strFind((prb_StringFindSpec) {.str = fromMain, .pattern = prb_STR("{"), .direction = prb_StringDirection_FromStart, .mode = prb_StringFindMode_AnyChar});
                        prb_assert(res.found);
                        mainLineEnd = res.matchByteIndex + mainFound.matchByteIndex + 1;
                    }

                    if (false) {
                        prb_writelnToStdout(prb_strSliceBetween(fileContentStr, mainLineStart, mainLineEnd));
                        prb_writelnToStdout(mainArgcName);
                        prb_writelnToStdout(mainArgvName);
                    }

                    prb_GrowingString gstr = prb_beginString(arena);
                    prb_String        fromMainEnd = prb_strSliceForward(fileContentStr, mainLineEnd);
                    prb_addStringSegment(&gstr, "%.*s\n", mainLineEnd, fileContentStr.ptr);
                    prb_addStringSegment(&gstr, "{\n");
                    prb_addStringSegment(&gstr, "FILE* logfile = fopen(\"%.*s\", \"a\");\n", prb_LIT(logfilePath));
                    prb_addStringSegment(&gstr, "for (int argIndex = 0; argIndex < %.*s; argIndex++) {\n", prb_LIT(mainArgcName));
                    prb_addStringSegment(&gstr, "char* arg = %.*s[argIndex];\n", prb_LIT(mainArgvName));
                    prb_addStringSegment(&gstr, "fwrite(arg, strlen(arg), 1, logfile);\n");
                    prb_addStringSegment(&gstr, "fwrite(\" \", 1, 1, logfile);\n");
                    prb_addStringSegment(&gstr, "}\n");
                    prb_addStringSegment(&gstr, "fwrite(\"\\n\", 1, 1, logfile);\n");
                    prb_addStringSegment(&gstr, "fclose(logfile);\n");
                    prb_addStringSegment(&gstr, "}\n");
                    prb_addStringSegment(&gstr, "%.*s", prb_LIT(fromMainEnd));

                    prb_String newMainContent = prb_endString(&gstr);
                    prb_String newMainPath = prb_fmt(arena, "%.*s.new", prb_LIT(sourceIter.curPath));
                    prb_assert(prb_writeEntireFile(arena, newMainPath, newMainContent.ptr, newMainContent.len));

                    srcPath = newMainPath;
                }

            } else {
                arrput(objsWithoutMain, outpath);
            }

            if (!prb_isFile(arena, outpath)) {
                anyObjRecompiled = true;
                prb_String cmd = prb_fmt(
                    arena,
                    "clang -g -Denablemultithread -Werror -Wfatal-errors %.*s -x c %.*s -c -o %.*s",
                    prb_LIT(mafftSuppressedWarningsStr),
                    prb_LIT(srcPath),
                    prb_LIT(outpath)
                );
                prb_writelnToStdout(cmd);
                prb_ProcessHandle cmdProcess = prb_execCmd(arena, cmd, prb_ProcessFlag_DontWait, (prb_String) {});
                prb_assert(cmdProcess.status == prb_ProcessStatus_Launched);
                arrput(objProccesses, cmdProcess);
            }
        }
    }
    prb_destroyPathFindIter(&sourceIter);
    prb_assert(prb_waitForProcesses(objProccesses, arrlen(objProccesses)) == prb_Success);

    // NOTE(sen) Executables
    bool printObjs = false;
    if (printObjs) {
        printObjList(arena, "With main", objsWithMain);
        printObjList(arena, "Without main", objsWithoutMain);
    }
    prb_ProcessHandle* exeProcesses = 0;
    prb_String         exeDir = prb_pathJoin(arena, mafftOutDir, prb_STR("exes"));
    prb_assert(prb_createDirIfNotExists(arena, exeDir) == prb_Success);
    // prb_assert(prb_clearDirectory(arena, exeDir) == prb_Success);
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
            prb_String exePath = prb_pathJoin(arena, exeDir, exeName);

            if (anyObjRecompiled || !prb_isFile(arena, exePath)) {
                prb_GrowingString gstr = prb_beginString(arena);
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
                prb_ProcessHandle exeHandle = prb_execCmd(arena, cmd, prb_ProcessFlag_DontWait, (prb_String) {});
                prb_assert(exeHandle.status == prb_ProcessStatus_Launched);
                arrput(exeProcesses, exeHandle);
            }
        }
    }
    prb_assert(prb_waitForProcesses(exeProcesses, arrlen(exeProcesses)) == prb_Success);

    // NOTE(sen) Copy from mafft if necessary
    if (false) {
        for (i32 index = 0; index < arrlen(objsWithoutMain); index++) {
            prb_String ogsrcname = prb_replaceExt(arena, prb_getLastEntryInPath(objsWithoutMain[index]), prb_STR("c"));
            prb_String ogsrcpath = prb_pathJoin(arena, coreDir, ogsrcname);
            copyFile(arena, ogsrcpath, srcDir);
        }

        prb_PathFindSpec headerSpec = srcIterSpec;
        headerSpec.pattern = prb_STR("*.h");
        prb_PathFindIterator headerIter = prb_createPathFindIter(headerSpec);
        while (prb_pathFindIterNext(&headerIter)) {
            copyFile(arena, headerIter.curPath, srcDir);
        }

        for (i32 index = 0; index < prb_arrayLength(srcFilesNoObj); index++) {
            copyFile(arena, prb_pathJoin(arena, coreDir, srcFilesNoObj[index]), srcDir);
        }
    }

    // NOTE(sen) Compile what I pulled out
    prb_String srcOutDir = prb_pathJoin(arena, outDir, prb_STR("src"));
    prb_assert(prb_createDirIfNotExists(arena, srcOutDir));
    prb_assert(prb_clearDirectory(arena, srcOutDir));
    srcIterSpec.dir = srcDir;
    sourceIter = prb_createPathFindIter(srcIterSpec);
    prb_ProcessHandle* srcCompileProcs = 0;
    prb_String*        srcObjPaths = 0;
    while (prb_pathFindIterNext(&sourceIter)) {
        prb_String inname = prb_getLastEntryInPath(sourceIter.curPath);
        if (notIn(inname, srcFilesNoObj, prb_arrayLength(srcFilesNoObj))) {
            prb_String outname = prb_replaceExt(arena, inname, prb_STR("obj"));
            prb_String outpath = prb_pathJoin(arena, srcOutDir, outname);
            arrput(srcObjPaths, outpath);
            prb_String cmd = prb_fmt(arena, "clang -g -Wall -Werror -Wfatal-errors -Denablemultithread -c %.*s -o %.*s", prb_LIT(sourceIter.curPath), prb_LIT(outpath));
            prb_writelnToStdout(cmd);
            prb_ProcessHandle proc = prb_execCmd(arena, cmd, 0, (prb_String) {});
            prb_assert(proc.status == prb_ProcessStatus_CompletedSuccess);
            arrput(srcCompileProcs, proc);
        }
    }
    prb_assert(prb_waitForProcesses(srcCompileProcs, arrlen(srcCompileProcs)));

    {
        prb_String srcLibPath = prb_pathJoin(arena, srcOutDir, prb_STR("align.a"));
        prb_String srcObjsStr = prb_stringsJoin(arena, srcObjPaths, arrlen(srcObjPaths), prb_STR(" "));
        prb_String cmd = prb_fmt(arena, "ar rcs %.*s %.*s", prb_LIT(srcLibPath), prb_LIT(srcObjsStr));
        prb_writelnToStdout(cmd);
        prb_ProcessHandle proc = prb_execCmd(arena, cmd, 0, (prb_String) {});
        prb_assert(proc.status == prb_ProcessStatus_CompletedSuccess);
    }

    prb_writelnToStdout(prb_fmt(arena, "total: %.2fms", prb_getMsFrom(scriptStart)));
    return 0;
}
