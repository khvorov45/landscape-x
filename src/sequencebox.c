#include <R.h>
#include <Rinternals.h>

SEXP
sequencebox_c(void) {
    SEXP result = PROTECT(allocVector(REALSXP, 1));
    (REAL(result))[0] = 0;
    UNPROTECT(1);
    return result;
}