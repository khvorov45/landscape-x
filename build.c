#include "cbuild.h"

typedef int32_t i32;
static int
strcmpSort(const void* a, const void* b) {
    const char* astr = ((prb_Str*)a)->ptr;
    const char* bstr = ((prb_Str*)b)->ptr;
    int         result = strcmp(astr, bstr);
    return result;
}

static void
printObjList(prb_Arena* arena, const char* title, prb_Str* objs) {
    i32 objsLen = arrlen(objs);
    qsort(objs, objsLen, sizeof(*objs), strcmpSort);
    prb_writeToStdout(prb_fmt(arena, "\n%s: %d\n", title, objsLen));
    for (i32 index = 0; index < objsLen; index++) {
        prb_writelnToStdout(arena, objs[index]);
    }
}

static bool
notIn(prb_Str val, prb_Str* arr, i32 arrLen) {
    bool result = true;
    for (i32 arrIndex = 0; arrIndex < arrLen && result; arrIndex++) {
        prb_Str arrVal = arr[arrIndex];
        if (prb_streq(val, arrVal)) {
            result = false;
        }
    }
    return result;
}

static void
copyFile(prb_Arena* arena, prb_Str path, prb_Str dir) {
    prb_Str                  copyPath = prb_pathJoin(arena, dir, prb_getLastEntryInPath(path));
    prb_ReadEntireFileResult fileContent = prb_readEntireFile(arena, path);
    prb_assert(fileContent.success);
    prb_assert(prb_writeEntireFile(arena, copyPath, fileContent.content.data, fileContent.content.len) == prb_Success);
}

int
main() {
    prb_TimeStart scriptStart = prb_timeStart();
    prb_Arena     arena_ = prb_createArenaFromVmem(1 * prb_GIGABYTE);
    prb_Arena*    arena = &arena_;
    prb_Str       rootDir = prb_getParentDir(arena, prb_STR(__FILE__));
    prb_Str       srcDir = prb_pathJoin(arena, rootDir, prb_STR("src"));
    prb_Str       outDir = prb_pathJoin(arena, rootDir, prb_STR("build-debug"));
    prb_assert(prb_createDirIfNotExists(arena, outDir) == prb_Success);
    // prb_assert(prb_clearDir(arena, outDir) == prb_Success);

    // NOTE(sen) Get mafft
    prb_Str mafftDir = prb_pathJoin(arena, rootDir, prb_STR("mafft"));
    if (!prb_isDir(arena, mafftDir) || prb_dirIsEmpty(arena, mafftDir)) {
        prb_Str cmd = prb_fmt(arena, "git clone https://gitlab.com/sysimm/mafft %.*s", prb_LIT(mafftDir));
        prb_writelnToStdout(arena, cmd);
        prb_Process proc = prb_createProcess(cmd, (prb_ProcessSpec) {});
        prb_assert(prb_launchProcesses(arena, &proc, 1, prb_Background_No));
    }

    prb_Str srcFilesNoObj[] = {
        prb_STR("JTT.c"),
        prb_STR("blosum.c"),
        prb_STR("iteration.c"),
    };

    // NOTE(sen) Compile mafft
    if (false) {
        prb_Str coreDir = prb_pathJoin(arena, mafftDir, prb_STR("core"));
        prb_Str mafftSuppressedWarnings[] = {
            prb_STR("-Wno-deprecated-non-prototype"),
            prb_STR("-Wno-unknown-directives"),
            prb_STR("-Wno-fortify-source"),
            prb_STR("-Wno-parentheses"),
            prb_STR("-Wno-format"),
            prb_STR("-Wno-incompatible-pointer-types"),
        };
        prb_Str mafftSuppressedWarningsStr = prb_stringsJoin(arena, mafftSuppressedWarnings, prb_arrayCount(mafftSuppressedWarnings), prb_STR(" "));

        prb_Str logfilePath = prb_pathJoin(arena, coreDir, prb_STR("logfile.txt"));
        prb_removePathIfExists(arena, logfilePath);

        // NOTE(sen) Find all source files
        prb_Str* allFilesInCore = prb_getAllDirEntries(arena, coreDir, prb_Recursive_No);

        // NOTE(sen) Objs

        prb_Str* objsWithMain = 0;
        prb_Str* objsWithoutMain = 0;
        prb_Str  mafftOutDir = prb_pathJoin(arena, outDir, prb_STR("mafft"));
        prb_assert(prb_createDirIfNotExists(arena, mafftOutDir) == prb_Success);
        prb_Str mafftObjDir = prb_pathJoin(arena, mafftOutDir, prb_STR("objs"));
        prb_assert(prb_createDirIfNotExists(arena, mafftObjDir) == prb_Success);

        // NOTE(sen) Force clean mafft
        // prb_assert(prb_clearDir(arena, mafftObjDir) == prb_Success);

        prb_Process* objProccesses = 0;
        bool         anyObjRecompiled = false;
        for (i32 fileIndex = 0; fileIndex < arrlen(allFilesInCore); fileIndex++) {
            prb_Str thisFile = allFilesInCore[fileIndex];
            if (prb_strEndsWith(thisFile, prb_STR(".c"))) {
                prb_Str inname = prb_getLastEntryInPath(thisFile);
                if (notIn(inname, srcFilesNoObj, prb_arrayCount(srcFilesNoObj))) {
                    prb_Str outname = prb_replaceExt(arena, inname, prb_STR("obj"));
                    prb_Str outpath = prb_pathJoin(arena, mafftObjDir, outname);

                    prb_ReadEntireFileResult fileRead = prb_readEntireFile(arena, thisFile);
                    prb_assert(fileRead.success);
                    prb_Str        srcPath = thisFile;
                    prb_StrScanner fileContentScanner = prb_createStrScanner(prb_strFromBytes(fileRead.content));
                    if (prb_strScannerMove(&fileContentScanner, (prb_StrFindSpec) {.pattern = prb_STR("main(")}, prb_StrScannerSide_AfterMatch)) {
                        arrput(objsWithMain, outpath);

                        prb_assert(prb_strScannerMove(&fileContentScanner, (prb_StrFindSpec) {.mode = prb_StrFindMode_LineBreak, .direction = prb_StrDirection_FromEnd}, prb_StrScannerSide_BeforeMatch));
                        prb_assert(prb_strScannerMove(&fileContentScanner, (prb_StrFindSpec) {.pattern = prb_STR("{")}, prb_StrScannerSide_AfterMatch));
                        prb_Str mainLine = fileContentScanner.betweenLastMatches;

                        bool    mainTakesArgs = false;
                        prb_Str mainArgcName = {};
                        prb_Str mainArgvName = {};
                        {
                            prb_StrScanner scanner = prb_createStrScanner(mainLine);
                            prb_assert(prb_strScannerMove(&scanner, (prb_StrFindSpec) {.pattern = prb_STR("(")}, prb_StrScannerSide_AfterMatch));
                            if (prb_strScannerMove(&scanner, (prb_StrFindSpec) {.pattern = prb_STR("int")}, prb_StrScannerSide_AfterMatch)) {
                                mainTakesArgs = true;

                                prb_assert(prb_strScannerMove(&scanner, (prb_StrFindSpec) {.pattern = prb_STR(",")}, prb_StrScannerSide_AfterMatch));
                                mainArgcName = prb_strTrim(scanner.betweenLastMatches);

                                bool insideAlpha = false;
                                i32  lastAlpha = 0;
                                i32  firstAlpha = 0;
                                for (i32 lineIndex = mainLine.len - 1; lineIndex >= 0; lineIndex--) {
                                    char ch = mainLine.ptr[lineIndex];
                                    bool isAlpha = (ch >= 'a' && ch <= 'z') || ((ch >= 'A' && ch <= 'Z'));
                                    if (isAlpha && !insideAlpha) {
                                        lastAlpha = lineIndex;
                                        insideAlpha = true;
                                    } else if (!isAlpha && insideAlpha) {
                                        firstAlpha = lineIndex + 1;
                                        break;
                                    }
                                }

                                mainArgvName = prb_strSlice(mainLine, firstAlpha, lastAlpha + 1);
                            }
                        }

                        if (mainTakesArgs) {
                            // NOTE(sen) Generate alternative file where main is modified

                            if (false) {
                                prb_writelnToStdout(arena, mainLine);
                                prb_writelnToStdout(arena, mainArgcName);
                                prb_writelnToStdout(arena, mainArgvName);
                            }

                            prb_GrowingStr gstr = prb_beginStr(arena);
                            prb_addStrSegment(&gstr, "%.*s", prb_LIT(fileContentScanner.beforeMatch));
                            prb_addStrSegment(&gstr, "%.*s\n", prb_LIT(fileContentScanner.match));
                            prb_addStrSegment(&gstr, "{\n");
                            prb_addStrSegment(&gstr, "FILE* logfile = fopen(\"%.*s\", \"a\");\n", prb_LIT(logfilePath));
                            prb_addStrSegment(&gstr, "for (int argIndex = 0; argIndex < %.*s; argIndex++) {\n", prb_LIT(mainArgcName));
                            prb_addStrSegment(&gstr, "char* arg = %.*s[argIndex];\n", prb_LIT(mainArgvName));
                            prb_addStrSegment(&gstr, "fwrite(arg, strlen(arg), 1, logfile);\n");
                            prb_addStrSegment(&gstr, "fwrite(\" \", 1, 1, logfile);\n");
                            prb_addStrSegment(&gstr, "}\n");
                            prb_addStrSegment(&gstr, "fwrite(\"\\n\", 1, 1, logfile);\n");
                            prb_addStrSegment(&gstr, "fclose(logfile);\n");
                            prb_addStrSegment(&gstr, "}\n");
                            prb_addStrSegment(&gstr, "%.*s", prb_LIT(fileContentScanner.afterMatch));

                            prb_Str newMainContent = prb_endStr(&gstr);
                            prb_Str newMainPath = prb_fmt(arena, "%.*s.new", prb_LIT(thisFile));
                            prb_assert(prb_writeEntireFile(arena, newMainPath, newMainContent.ptr, newMainContent.len));

                            srcPath = newMainPath;
                        }

                    } else {
                        arrput(objsWithoutMain, outpath);
                    }

                    if (!prb_isFile(arena, outpath)) {
                        anyObjRecompiled = true;
                        prb_Str cmd = prb_fmt(
                            arena,
                            "clang -g -Denablemultithread -Werror -Wfatal-errors %.*s -x c %.*s -c -o %.*s",
                            prb_LIT(mafftSuppressedWarningsStr),
                            prb_LIT(srcPath),
                            prb_LIT(outpath)
                        );
                        prb_writelnToStdout(arena, cmd);
                        prb_Process cmdProcess = prb_createProcess(cmd, (prb_ProcessSpec) {});
                        arrput(objProccesses, cmdProcess);
                    }
                }
            }
        }

        prb_assert(prb_launchProcesses(arena, objProccesses, arrlen(objProccesses), prb_Background_Yes));
        prb_assert(prb_waitForProcesses(objProccesses, arrlen(objProccesses)));

        // NOTE(sen) Executables
        bool printObjs = false;
        if (printObjs) {
            printObjList(arena, "With main", objsWithMain);
            printObjList(arena, "Without main", objsWithoutMain);
        }
        prb_Process* exeProcesses = 0;
        prb_Str      exeDir = prb_pathJoin(arena, mafftOutDir, prb_STR("exes"));
        prb_assert(prb_createDirIfNotExists(arena, exeDir) == prb_Success);

        // NOTE(khvorov) Force clean executables
        // prb_assert(prb_clearDir(arena, exeDir) == prb_Success);
        for (i32 withMainIndex = 0; withMainIndex < arrlen(objsWithMain); withMainIndex++) {
            prb_Str objWithMain = objsWithMain[withMainIndex];
            prb_Str objWithMainName = prb_getLastEntryInPath(objWithMain);
            if (!prb_streq(objWithMainName, prb_STR("interface.obj"))) {
                prb_StrFindResult dotFind = prb_strFind(objWithMainName, (prb_StrFindSpec) {
                                                                             .pattern = prb_STR("."),
                                                                             .direction = prb_StrDirection_FromEnd,
                                                                         });
                prb_assert(dotFind.found);
                prb_Str exeName = dotFind.beforeMatch;
                prb_Str exePath = prb_pathJoin(arena, exeDir, exeName);

                if (anyObjRecompiled || !prb_isFile(arena, exePath)) {
                    prb_GrowingStr gstr = prb_beginStr(arena);
                    prb_addStrSegment(&gstr, "clang %.*s", prb_LIT(objWithMain));
                    for (i32 withoutMainIndex = 0; withoutMainIndex < arrlen(objsWithoutMain); withoutMainIndex++) {
                        prb_Str objWithoutMain = objsWithoutMain[withoutMainIndex];
                        prb_Str objWithoutMainName = prb_getLastEntryInPath(objWithoutMain);

                        bool addThisObj = true;
                        if (prb_streq(objWithoutMainName, prb_STR("pairlocalalign.obj"))) {
                            prb_Str exesToNotCouple[] = {prb_STR("pairash.obj"), prb_STR("rnatest.obj"), prb_STR("multi2hat3s.obj")};
                            addThisObj = notIn(objWithMainName, exesToNotCouple, prb_arrayCount(exesToNotCouple));
                        }

                        if (addThisObj) {
                            prb_addStrSegment(&gstr, " %.*s", prb_LIT(objWithoutMain));
                        }
                    }
                    prb_addStrSegment(&gstr, " -o %.*s -lpthread -lm", prb_LIT(exePath));
                    prb_Str cmd = prb_endStr(&gstr);

                    prb_writelnToStdout(arena, cmd);
                    prb_Process exeHandle = prb_createProcess(cmd, (prb_ProcessSpec) {});
                    arrput(exeProcesses, exeHandle);
                }
            }
        }

        prb_assert(prb_launchProcesses(arena, exeProcesses, arrlen(exeProcesses), prb_Background_Yes));
        prb_assert(prb_waitForProcesses(exeProcesses, arrlen(exeProcesses)));

        // NOTE(sen) Copy from mafft if necessary
        if (false) {
            for (i32 index = 0; index < arrlen(objsWithoutMain); index++) {
                prb_Str ogsrcname = prb_replaceExt(arena, prb_getLastEntryInPath(objsWithoutMain[index]), prb_STR("c"));
                prb_Str ogsrcpath = prb_pathJoin(arena, coreDir, ogsrcname);
                copyFile(arena, ogsrcpath, srcDir);
            }

            for (i32 fileIndex = 0; fileIndex < arrlen(allFilesInCore); fileIndex++) {
                prb_Str thisFile = allFilesInCore[fileIndex];
                if (prb_strEndsWith(thisFile, prb_STR(".h"))) {
                    copyFile(arena, thisFile, srcDir);
                }
            }

            for (i32 index = 0; index < prb_arrayCount(srcFilesNoObj); index++) {
                copyFile(arena, prb_pathJoin(arena, coreDir, srcFilesNoObj[index]), srcDir);
            }
        }
    }

    // NOTE(sen) Compile what I pulled out
    prb_Str srcOutDir = prb_pathJoin(arena, outDir, prb_STR("src"));
    prb_assert(prb_createDirIfNotExists(arena, srcOutDir));
    prb_assert(prb_clearDir(arena, srcOutDir));
    prb_Str*     allFilesInMysrc = prb_getAllDirEntries(arena, srcDir, prb_Recursive_No);
    prb_Process* srcCompileProcs = 0;
    prb_Str*     srcObjPaths = 0;
    prb_Str*     srcObjPathsLog = 0;
    for (i32 fileIndex = 0; fileIndex < arrlen(allFilesInMysrc); fileIndex++) {
        prb_Str thisFile = allFilesInMysrc[fileIndex];
        if (prb_strEndsWith(thisFile, prb_STR(".c"))) {
            prb_Str inname = prb_getLastEntryInPath(thisFile);
            if (notIn(inname, srcFilesNoObj, prb_arrayCount(srcFilesNoObj))) {
                prb_Str outname = prb_replaceExt(arena, inname, prb_STR("obj"));
                prb_Str outpath = prb_pathJoin(arena, srcOutDir, outname);
                arrput(srcObjPaths, outpath);
                prb_Str outlog = prb_replaceExt(arena, outpath, prb_STR("log"));
                arrput(srcObjPathsLog, outlog);
                prb_Str cmd = prb_fmt(arena, "clang -g -Wall -Werror -Wextra -fno-caret-diagnostics -c %.*s -o %.*s", prb_LIT(thisFile), prb_LIT(outpath));
                prb_writelnToStdout(arena, cmd);
                prb_Process proc = prb_createProcess(cmd, (prb_ProcessSpec) {.redirectStderr = true, .stderrFilepath = outlog});
                arrput(srcCompileProcs, proc);
            }
        }
    }

    prb_assert(prb_launchProcesses(arena, srcCompileProcs, arrlen(srcCompileProcs), prb_Background_Yes));
    if (prb_waitForProcesses(srcCompileProcs, arrlen(srcCompileProcs)) == prb_Failure) {
        for (i32 logIndex = 0; logIndex < arrlen(srcObjPathsLog); logIndex++) {
            prb_Str                  logfile = srcObjPathsLog[logIndex];
            prb_ReadEntireFileResult readRes = prb_readEntireFile(arena, logfile);
            prb_assert(readRes.success);
            prb_writeToStdout(prb_strFromBytes(readRes.content));
        }
        prb_assert(!"compilation failed");
    }

    {
        prb_Str srcLibPath = prb_pathJoin(arena, srcOutDir, prb_STR("align.a"));
        prb_Str srcObjsStr = prb_stringsJoin(arena, srcObjPaths, arrlen(srcObjPaths), prb_STR(" "));
        prb_Str cmd = prb_fmt(arena, "ar rcs %.*s %.*s", prb_LIT(srcLibPath), prb_LIT(srcObjsStr));
        prb_writelnToStdout(arena, cmd);
        prb_Process proc = prb_createProcess(cmd, (prb_ProcessSpec) {});
        prb_assert(prb_launchProcesses(arena, &proc, 1, prb_Background_No));
    }

    prb_writelnToStdout(arena, prb_fmt(arena, "total: %.2fms", prb_getMsFrom(scriptStart)));
    return 0;
}
