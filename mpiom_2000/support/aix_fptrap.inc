#if defined(__xlC__)
! Fortran interface of C wrapper functions for fptrap.h defined interface
! (c) Thomas Jahns, 2009
! (c) DKRZ GmbH, 2009
!  EXTERNAL f_fp_any_enabled, f_fp_is_enabled, f_fp_enable, f_fp_disable, &
!       f_fp_enable_all, f_fp_disable_all
! replace 4 with FPSCR_KIND
  USE xlf_fp_util
INTERFACE
  FUNCTION f_fp_any_enabled()
    LOGICAL :: f_fp_any_enabled
  END FUNCTION f_fp_any_enabled
END INTERFACE
INTERFACE
  FUNCTION f_fp_is_enabled(mask)
    LOGICAL :: f_fp_is_enabled
    INTEGER(4), INTENT(in) :: mask
  END FUNCTION f_fp_is_enabled
END INTERFACE
INTERFACE
  SUBROUTINE f_fp_enable_all()
  END SUBROUTINE f_fp_enable_all
END INTERFACE
INTERFACE
  SUBROUTINE f_fp_enable(mask)
    INTEGER(4), INTENT(in) :: mask
  END SUBROUTINE f_fp_enable
END INTERFACE
INTERFACE
  SUBROUTINE f_fp_disable_all()
  END SUBROUTINE f_fp_disable_all
END INTERFACE
INTERFACE
  SUBROUTINE f_fp_disable(mask)
    INTEGER(4), INTENT(in) :: mask
  END SUBROUTINE f_fp_disable
END INTERFACE
INTERFACE
  FUNCTION f_fp_trap(mode)
    INTEGER(4) :: f_fp_trap
    INTEGER(4), INTENT(in) :: mode
  END FUNCTION f_fp_trap
END INTERFACE
!!$  INTEGER(4) TRP_INVALID,TRP_DIV_BY_ZERO,TRP_OVERFLOW,TRP_UNDERFLOW,&
!!$       TRP_INEXACT
!!$  PARAMETER (TRP_INVALID=128)
!!$  PARAMETER (TRP_OVERFLOW=64)
!!$  PARAMETER (TRP_UNDERFLOW=32)
!!$  PARAMETER (TRP_DIV_BY_ZERO=16)
!!$  PARAMETER (TRP_INEXACT=8)
#endif
