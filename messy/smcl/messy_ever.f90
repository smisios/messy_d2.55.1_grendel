! **********************************************************************
MODULE messy_ever
! **********************************************************************

  ! EXPLOSIVE VOLCANIC ERUPTIONS (EVER)
  !
  ! CORE MODULE (MESSy/SMCL)
  !
  ! Author: Matthias Kohl, MPICH, February 2020
  !

  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: TINY

  PUBLIC :: DP
  PUBLIC :: XERF

  ! ----------- <

  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modstr = 'ever'
  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modver = '1.0'
  REAL(DP),        PARAMETER, PUBLIC :: ZERO_EPS = TINY(0.0_DP)    ! zero



CONTAINS

  ! ----------------------------------------------------------------------
  REAL(dp) FUNCTION XERF(X)

     ! ***** ERROR FUNCTION WITH RATIONAL APPROXIMATION
     ! ***** M.ABRAMOWITZ, I.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS
     ! ***** DOVER, 10th PRINTING 1972, p.299, # 7.1.26

     IMPLICIT NONE

     ! function parameters

     REAL(dp), INTENT(in) :: X

     ! parameters

     REAL(dp), PARAMETER :: P  = 0.3275911_dp
     REAL(dp), PARAMETER :: A1 = 0.254829592_dp
     REAL(dp), PARAMETER :: A2 =-0.284496736_dp
     REAL(dp), PARAMETER :: A3 = 1.421413741_dp
     REAL(dp), PARAMETER :: A4 =-1.453152027_dp
     REAL(dp), PARAMETER :: A5 = 1.061405429_dp

     ! local variables

     REAL(dp) :: T,ARG

     T    = 1.0_dp / ( 1.0_dp + P * ABS(X) )
     ARG  = MIN( X**2, 75.0_dp )
     XERF = 1.0_dp - ( A1*T + A2*T**2 + A3*T**3 + A4*T**4 + A5*T**5 ) &
            * EXP( -ARG )
     XERF = SIGN( XERF, X )

     RETURN

     END FUNCTION XERF

  ! -------------------------------------------------------------------------

! **********************************************************************
END MODULE messy_ever
! **********************************************************************
