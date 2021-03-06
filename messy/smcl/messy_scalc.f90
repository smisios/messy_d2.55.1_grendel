! **********************************************************************
MODULE messy_scalc
! **********************************************************************

  USE messy_main_constants_mem, ONLY: DP

  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modstr = 'scalc'
  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modver = '1.9'

  INTEGER,         PARAMETER, PUBLIC :: STRLEN = 1000

  PRIVATE

  PUBLIC :: DP

  INTEGER, PARAMETER, PUBLIC :: IS_2D   = 1
  INTEGER, PARAMETER, PUBLIC :: IS_3D_Z = 2
  INTEGER, PARAMETER, PUBLIC :: IS_3D_N = 3

  PUBLIC :: delta2frac

CONTAINS

  ELEMENTAL SUBROUTINE delta2frac(frac, delta, Rstd, nr)

    IMPLICIT NONE

    ! I/O
    REAL(DP), INTENT(OUT) :: frac   ! fration of rare (minor)
    REAL(DP), INTENT(IN)  :: delta  ! permil
    REAL(DP), INTENT(IN)  :: Rstd   ! ratio of standard
    REAL(DP), INTENT(IN)  :: nr     ! number of (isotope) atoms / molecule

    ! LOCAL
    REAL(DP) :: ratio

    ratio = Rstd * (delta / 1000.0_dp + 1.0_dp)
    frac  = nr * ratio / (1.0_dp + ratio)

  END SUBROUTINE delta2frac

! **********************************************************************
END MODULE messy_scalc
! **********************************************************************
