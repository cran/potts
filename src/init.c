
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "potts.h"

static R_NativePrimitiveArgType pack_types[6] =
    {INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, RAWSXP};

static R_NativeArgStyle pack_styles[6] =
    {R_ARG_IN, R_ARG_IN, R_ARG_IN, R_ARG_IN, R_ARG_IN, R_ARG_OUT};

static R_NativePrimitiveArgType inspect_types[4] =
    {RAWSXP, INTSXP, INTSXP, INTSXP};

static R_NativeArgStyle inspect_styles[4] =
    {R_ARG_IN, R_ARG_OUT, R_ARG_OUT, R_ARG_OUT};

static R_NativePrimitiveArgType unpack_types[6] =
    {RAWSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP};

static R_NativeArgStyle unpack_styles[6] =
    {R_ARG_IN, R_ARG_IN, R_ARG_IN, R_ARG_IN, R_ARG_IN, R_ARG_OUT};

static R_NativePrimitiveArgType potts_types[15] =
    {RAWSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP,
    INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP};

static R_NativeArgStyle potts_styles[15] =
    {R_ARG_IN_OUT, R_ARG_IN, R_ARG_IN, R_ARG_IN, R_ARG_IN, R_ARG_IN,
    R_ARG_OUT, R_ARG_IN, R_ARG_OUT, R_ARG_OUT, R_ARG_OUT, R_ARG_OUT,
    R_ARG_OUT, R_ARG_OUT, R_ARG_OUT};

static R_CMethodDef cMethods[] = {
    {"packPotts", (DL_FUNC) &packPotts, 6, pack_types, pack_styles},
    {"inspectPotts", (DL_FUNC) &inspectPotts, 4, inspect_types, inspect_styles},
    {"unpackPotts", (DL_FUNC) &unpackPotts, 6, unpack_types, unpack_styles},
    {"potts", (DL_FUNC) &potts, 15, potts_types, potts_styles},
    {NULL, NULL, 0, NULL, NULL}
};
 
static R_CallMethodDef callMethods[]  = {
    {NULL, NULL, 0}
};

void R_init_aster(DllInfo *info)
{
    R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
}

