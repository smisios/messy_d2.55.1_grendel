! ************************************************************************
MODULE messy_tnudge
! ************************************************************************

  ! MODULE FOR TRACER NUDGING  (MESSy-SMCL)
  !
  ! Authors:
  !    Patrick Joeckel, MPICH, December 2003
  !       - complete namelist control setup for
  !         arbitrary tracers/channel objects
  !

  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: TINY

  PUBLIC :: DP

  ! ----------- <

  ! GLOBAL PARAMETERS
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'tnudge'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '3.1'
  REAL(DP),         PARAMETER, PUBLIC :: EPS = TINY(0.0_DP) ! zero

! ************************************************************************
END MODULE messy_tnudge
! ************************************************************************
