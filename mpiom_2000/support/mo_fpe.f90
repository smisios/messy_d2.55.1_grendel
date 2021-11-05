MODULE mo_fpe
#if defined(__xlC__)
# include "aix_fptrap.inc"
# include <fexcp.h>
#endif
  IMPLICIT NONE
  PRIVATE
#if defined(__xlC__)
  INTEGER(FPSCR_KIND) :: oldmode
#endif
  LOGICAL :: tracing_active = .FALSE.
  PUBLIC :: enable_fpe_tracing, disable_fpe_tracing
CONTAINS
  SUBROUTINE enable_fpe_tracing
    IF (.NOT. tracing_active) THEN
#if defined(__xlC__)
      CALL SIGNAL(SIGFPE,xl__trce)
      CALL f_fp_enable(IOR(IOR(TRP_INVALID,TRP_DIV_BY_ZERO),TRP_OVERFLOW))
      oldmode = f_fp_trap(FP_TRAP_FASTMODE)
#endif
      tracing_active = .TRUE.
    END IF
  END SUBROUTINE enable_fpe_tracing
  SUBROUTINE disable_fpe_tracing
    IF (tracing_active) THEN
#if defined(__xlC__)
      CALL f_fp_disable_all
      oldmode = f_fp_trap(oldmode)
#endif
      tracing_active = .FALSE.
    END IF
  END SUBROUTINE disable_fpe_tracing
END MODULE mo_fpe
