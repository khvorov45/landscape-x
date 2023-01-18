#include "mltaln.h"

void
ErrorExit(char* message) {
    fprintf(stderr, "%s\n", message);
    exit(1);
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
