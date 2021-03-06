! **************************************************************************
MODULE messy_mmd2way
! **************************************************************************

  ! MODULE FOR exchanging data with another model of finer resolution
  !
  ! Author: Klaus Ketelsen, MPICH, Oct 2008
  !
  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP, I4, I8

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: DP, I4, I8

  ! ----------- <

  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modstr = 'mmd2way'
  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modver = '1.1'

! **************************************************************************
END MODULE messy_mmd2way
! **************************************************************************
