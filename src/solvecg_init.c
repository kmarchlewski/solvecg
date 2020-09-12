#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>

extern SEXP invert_cg(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP lehmer(SEXP);
extern SEXP solve_cg(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"invert_cg", (DL_FUNC) &invert_cg, 5},
    {"lehmer",    (DL_FUNC) &lehmer,    1},
    {"solve_cg",  (DL_FUNC) &solve_cg,  6},
    {NULL, NULL, 0}
};

void R_init_solvecg(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
