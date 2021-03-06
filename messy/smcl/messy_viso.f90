! **************************************************************************
MODULE messy_viso
! **************************************************************************

  ! MODULE FOR VALUES ON (HORIZONTAL) ISOSURFACES
  !
  ! Author: Patrick Joeckel, MPICH, Feb 2004
  !
  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: DP

  ! ----------- <

  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modstr = 'viso'
  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modver = '2.3'

! **************************************************************************
END MODULE messy_viso
! **************************************************************************
