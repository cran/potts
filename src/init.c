
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include "potts.h"

static R_NativePrimitiveArgType pack_types[6] =
    {INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, RAWSXP};

static R_NativePrimitiveArgType inspect_types[4] =
    {RAWSXP, INTSXP, INTSXP, INTSXP};

static R_NativePrimitiveArgType unpack_types[6] =
    {RAWSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP};

static R_NativePrimitiveArgType potts_types[15] =
    {RAWSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, LGLSXP,
    INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP};

static R_NativePrimitiveArgType out_types[3] = {RAWSXP, INTSXP, INTSXP};

static R_CMethodDef cMethods[] = {
    {"packPotts", (DL_FUNC) &packPotts, 6, pack_types},
    {"inspectPotts", (DL_FUNC) &inspectPotts, 4, inspect_types},
    {"unpackPotts", (DL_FUNC) &unpackPotts, 6, unpack_types},
    {"potts", (DL_FUNC) &potts, 15, potts_types},
    {"outfun_shutdown", (DL_FUNC) &outfun_shutdown, 0, NULL},
    {"outfun_len_init", (DL_FUNC) &outfun_len_init, 3, out_types},
    {NULL, NULL, 0, NULL}
};
 
static R_CallMethodDef callMethods[]  = {
    {"outfun_setup", (DL_FUNC) &outfun_setup, 2},
    {NULL, NULL, 0}
};

void attribute_visible R_init_potts(DllInfo *info)
{
    R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}

