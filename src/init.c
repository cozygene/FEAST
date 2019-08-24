#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _FEAST_schur(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_FEAST_schur", (DL_FUNC) &_FEAST_schur, 2},
    {NULL, NULL, 0}
};

void R_init_FEAST(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}