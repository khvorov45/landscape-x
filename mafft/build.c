#include "programmable_build.h"

int
main() {
    prb_Arena arena = prb_createArenaFromVmem(1 * prb_GIGABYTE);
    prb_String rootDir = prb_getParentDir(&arena, prb_STR(__FILE__));
    prb_String coreDir = prb_pathJoin(&arena, rootDir, prb_STR("core"));
    prb_String outDir = prb_pathJoin(&arena, rootDir, prb_STR("build-debug"));
    prb_createDirIfNotExists(&arena, outDir);
    prb_clearDirectory(&arena, outDir);
    prb_PathFindIterator iter = prb_createPathFindIter((prb_PathFindSpec){.arena = &arena, .dir = coreDir, .mode = prb_PathFindMode_Glob, .glob.pattern = prb_STR("*.c")});
    prb_ProcessHandle* processes = 0;
    while (prb_pathFindIterNext(&iter)) {
        prb_String inname = prb_getLastEntryInPath(iter.curPath);
        prb_Bytes fileContent = prb_readEntireFile(&arena, iter.curPath);
        prb_String fileContentStr = (prb_String){(const char*)fileContent.data, fileContent.len};
        prb_StringFindResult mainFound = prb_strFind((prb_StringFindSpec){.str = fileContentStr, .pattern = prb_STR("main("), .mode = prb_StringFindMode_Exact});
        if (mainFound.found) {
            //prb_writelnToStdout(iter.curPath);
        }
        if (!prb_streq(inname, prb_STR("JTT.c")) && !prb_streq(inname, prb_STR("blosum.c"))) {
            prb_String outname = prb_replaceExt(&arena, inname, prb_STR("obj"));
            prb_String outpath = prb_pathJoin(&arena, outDir, outname);
            prb_LastModResult lastModIn = prb_getLastModifiedFromPath(&arena, iter.curPath);
            prb_assert(lastModIn.success);
            prb_LastModResult lastModOut = prb_getLastModifiedFromPath(&arena, outpath);
            if (!lastModOut.success || lastModOut.timestamp < lastModIn.timestamp) {
                prb_String cmd = prb_fmtAndPrintln(&arena, "clang -Denablemultithread -Werror %.*s -c -o %.*s", prb_LIT(iter.curPath), prb_LIT(outpath));
                prb_ProcessHandle cmdProcess = prb_execCmd(&arena, cmd, 0, (prb_String){});
                prb_assert(cmdProcess.status == prb_ProcessStatus_CompletedSuccess);
                arrput(processes, cmdProcess);
            }
        }
    }
    prb_destroyPathFindIter(&iter);
    prb_assert(prb_waitForProcesses(processes, arrlen(processes)) == prb_Success);
    return 0;
}
