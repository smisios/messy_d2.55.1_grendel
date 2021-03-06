! **********************************************************************
!
! MODULE FOR TESTING MASS CONSERVATION AND MONOTONICITY OF
! EULERIAN TRANSPORT ALGORITHMS WITH PROGNOSTIC TRACERS
!
! CORE MODULE
!
! Author : Patrick Joeckel, MPICH, July  2002
!
! References:
!
! **********************************************************************

! **********************************************************************
MODULE messy_ptrac
! **********************************************************************

  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: DP

  ! ----------- <

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'ptrac'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '3.1'

! **********************************************************************
END MODULE messy_ptrac
! **********************************************************************
