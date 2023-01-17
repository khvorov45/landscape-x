#include "mltaln.h"

#define DEBUG 0
#define IODEBUG 0

void
ErrorExit(char* message) {
    fprintf(stderr, "%s\n", message);
    exit(1);
}

void
WriteFloatHat2_pointer_halfmtx(Context* ctx, FILE* hat2p, int locnjob, const char* const* name, double** mtx) {
    int    i, j, ijsa;
    double max = 0.0;
    for (i = 0; i < locnjob - 1; i++)
        for (j = 1; j < locnjob - i; j++)
            if (mtx[i][j] > max)
                max = mtx[i][j];

    fprintf(hat2p, "%5d\n", 1);
    fprintf(hat2p, "%5d\n", locnjob);
    fprintf(hat2p, " %#6.3f\n", max * 2.5);

    for (i = 0; i < locnjob; i++)
        fprintf(hat2p, "%4d. %s\n", i + 1, name[i]);
    for (i = 0; i < locnjob; i++) {
        for (j = i + 1; j < ctx->njob; j++) {
            fprintf(hat2p, DFORMAT, mtx[i][j - i]);
            ijsa = j - i;
            if (ijsa % 12 == 0 || ijsa == locnjob - i - 1)
                fprintf(hat2p, "\n");
        }
    }
}

void
WriteHat2_part_pointer(FILE* hat2p, int locnjob, int nadd, const char* const* name, double** mtx) {
    int    i, j;
    int    norg = locnjob - nadd;
    double max = 0.0;

    fprintf(hat2p, "%5d\n", 1);
    fprintf(hat2p, "%5d\n", locnjob);
    fprintf(hat2p, " %#6.3f\n", max * 2.5);

    for (i = 0; i < locnjob; i++)
        fprintf(hat2p, "%4d. %s\n", i + 1, name[i]);
    for (i = 0; i < norg; i++) {
        for (j = 0; j < nadd; j++) {
            fprintf(hat2p, DFORMAT, mtx[i][j]);
            if ((j + 1) % 12 == 0 || j == nadd - 1)
                fprintf(hat2p, "\n");
        }
    }
}

void
initSignalSM(Context* ctx) {
    if (!ctx->ppid) {
        ctx->signalSM = NULL;
    }
}

void
initFiles(Context* ctx) {
    char pname[100];
    if (ctx->ppid)
        sprintf(pname, "/tmp/pre.%d", ctx->ppid);
    else
        sprintf(pname, "pre");
    ctx->prep_g = fopen(pname, "w");
    if (!ctx->prep_g)
        ErrorExit("Cannot open pre");
    ctx->trap_g = fopen("trace", "w");
    if (!ctx->trap_g)
        ErrorExit("cannot open trace");
    fprintf(ctx->trap_g, "PID = %d\n", getpid());
    fflush(ctx->trap_g);
}

void
initlocalhom1(LocalHom* lh) {
    lh->start1 = -1;
    lh->end1 = -1;
    lh->start2 = -1;
    lh->end2 = -1;
    lh->opt = -1.0;
    lh->next = NULL;
    lh->nokori = 0;
    lh->extended = -1;
    lh->last = lh;
    lh->korh = 'h';
}

static void
showaamtxexample() {
    fprintf(stderr, "Format error in aa matrix\n");
    fprintf(stderr, "# Example:\n");
    fprintf(stderr, "# comment\n");
    fprintf(stderr, "   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V\n");
    fprintf(stderr, "A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0\n");
    fprintf(stderr, "R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3\n");
    fprintf(stderr, "...\n");
    fprintf(stderr, "V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4\n");
    fprintf(stderr, "frequency 0.07 0.05 0.04 0.05 0.02 .. \n");
    fprintf(stderr, "# Example end\n");
    fprintf(stderr, "Only the lower half is loaded\n");
    fprintf(stderr, "The last line (frequency) is optional.\n");
    exit(1);
}

double*
loadaamtx(Context* ctx, int* rescalept) {
    int      i, j, k, ii, jj;
    double*  val;
    double** raw;
    int*     map;
    char*    aaorder = "ARNDCQEGHILKMFPSTWYV";
    char*    inorder;
    char*    line;
    char*    ptr1;
    char*    ptr2;
    char*    mtxfname = "_aamtx";
    FILE*    mf;
    char     key[1000];

    raw = AllocateDoubleMtx(21, 20);
    val = AllocateDoubleVec(420);
    map = AllocateIntVec(20);

    if (ctx->dorp != 'p') {
        fprintf(stderr, "User-defined matrix is not supported for DNA\n");
        exit(1);
    }

    mf = fopen(mtxfname, "r");
    if (mf == NULL) {
        fprintf(stderr, "Cannot open the _aamtx file\n");
        exit(1);
    }

    inorder = calloc(1000, sizeof(char));
    line = calloc(1000, sizeof(char));

    while (!feof(mf)) {
        fgets(inorder, 999, mf);
        if (inorder[0] != '#')
            break;
    }
    ptr1 = ptr2 = inorder;
    while (*ptr2) {
        if (isalpha(*ptr2)) {
            *ptr1 = toupper(*ptr2);
            ptr1++;
        }
        ptr2++;
    }
    inorder[20] = 0;

    for (i = 0; i < 20; i++) {
        ptr2 = strchr(inorder, aaorder[i]);
        if (ptr2 == NULL) {
            fprintf(stderr, "%c: not found in the first 20 letters.\n", aaorder[i]);
            showaamtxexample();
        } else {
            map[i] = ptr2 - inorder;
        }
    }

    i = 0;
    while (!feof(mf)) {
        fgets(line, 999, mf);
        //		fprintf( stderr, "line = %s\n", line );
        if (line[0] == '#')
            continue;
        ptr1 = line;
        //		fprintf( stderr, "line = %s\n", line );
        for (j = 0; j <= i; j++) {
            while (!isdigit(*ptr1) && *ptr1 != '-' && *ptr1 != '.')
                ptr1++;

            raw[i][j] = atof(ptr1);
            //			fprintf( stderr, "raw[][]=%f, %c-%c %d-%d\n", raw[i][j], inorder[i], inorder[j], i, j );
            ptr1 = strchr(ptr1, ' ');
            if (ptr1 == NULL && j < i)
                showaamtxexample();
        }
        i++;
        if (i > 19)
            break;
    }

    *rescalept = 1;
    for (i = 0; i < 20; i++)
        raw[20][i] = -1.0;
    while (!feof(mf)) {
        fgets(line, 999, mf);

        sscanf(line, "%s", key);

        if (!strcmp(key, "norescale")) {
            reporterr("no rescale\n");
            *rescalept = 0;
            break;
        }
        //		else if( line[0] == 'f' )
        else if (!strcmp(key, "frequency")) {
            //			fprintf( stderr, "found! line = %s\n", line );
            ptr1 = line;
            for (j = 0; j < 20; j++) {
                while (!isdigit(*ptr1) && *ptr1 != '-' && *ptr1 != '.')
                    ptr1++;

                raw[20][j] = atof(ptr1);
                //				fprintf( stderr, "raw[20][]=%f, %c %d\n", raw[20][j], inorder[i], j );
                ptr1 = strchr(ptr1, ' ');
                if (ptr1 == NULL && j < 19)
                    showaamtxexample();
            }
            break;
        }
    }

    k = 0;
    for (i = 0; i < 20; i++) {
        for (j = 0; j <= i; j++) {
            if (i != j) {
                ii = MAX(map[i], map[j]);
                jj = MIN(map[i], map[j]);
            } else
                ii = jj = map[i];
            val[k++] = raw[ii][jj];
            //			fprintf( stderr, "%c-%c, %f\n", aaorder[i], aaorder[j], val[k-1] );
        }
    }
    for (i = 0; i < 20; i++)
        val[400 + i] = raw[20][map[i]];

    fprintf(stderr, "inorder = %s\n", inorder);
    fclose(mf);
    free(inorder);
    free(line);
    FreeDoubleMtx(raw);
    free(map);
    return (val);
}

void
freelocalhom1(LocalHom* lh) {
    if (lh == NULL)
        return;
    LocalHom* tmpptr = lh;
    LocalHom* ppp;
    for (; tmpptr; tmpptr = ppp) {
        ppp = tmpptr->next;
        if (tmpptr != lh) {
            free(tmpptr);
            continue;
        }
        tmpptr->start1 = -1;
        tmpptr->end1 = -1;
        tmpptr->start2 = -1;
        tmpptr->end2 = -1;
        tmpptr->opt = -1.0;
        tmpptr->next = NULL;
        tmpptr->nokori = 0;
        tmpptr->extended = -1;
        tmpptr->last = tmpptr;
        tmpptr->korh = 'h';
    }
}
void
reporterr(const char* str, ...) {
    va_list args;
    va_start(args, str);
    vfprintf(stderr, str, args);
    va_end(args);
    return;
}
