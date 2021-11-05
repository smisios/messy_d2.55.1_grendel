/*
 * Activate or deactivate generation and trapping of floating point
 * exceptions.  Wrapper for access via fortran.
 */
#if defined(__xlC__)
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <fptrap.h>
#include "cfortran.h"

FCALLSCFUN0(LOGICAL,fp_any_enable,F_FP_ANY_ENABLE,f_fp_any_enable)

FCALLSCFUN1(LOGICAL,fp_is_enabled,F_FP_IS_ENABLED,f_fp_is_enabled,\
            INT)

FCALLSCSUB0(fp_enable_all,F_FP_ENABLE_all, f_fp_enable_all)

FCALLSCSUB1(fp_enable,F_FP_ENABLE,f_fp_enable,\
            INT)

FCALLSCSUB0(fp_disable_all,F_FP_DISABLE_all, f_fp_disable_all)

FCALLSCSUB1(fp_disable,F_FP_DISABLE,f_fp_disable,\
            INT)

FCALLSCFUN1(INT,fp_trap,F_FP_TRAP,f_fp_trap,INT)
#endif
