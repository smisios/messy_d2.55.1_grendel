! **************************************************************************
MODULE messy_mmd2way_parent
! **************************************************************************

  ! MODULE FOR exchanging data with another model of finer resolution
  !
  ! Author: Klaus Ketelsen, MPICH, Oct 2008
  !
  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: DP

  ! ----------- <

  CHARACTER(LEN=*),PARAMETER, PUBLIC :: submodstr = 'mmd2way_parent'
  CHARACTER(LEN=*),PARAMETER, PUBLIC :: submodver = '1.0'

! **************************************************************************
END MODULE messy_mmd2way_parent
! **************************************************************************
