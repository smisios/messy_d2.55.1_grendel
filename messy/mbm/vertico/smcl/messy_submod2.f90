! **************************************************************************
MODULE messy_submod2
! **************************************************************************

  ! MODULE FOR ILLUSTRATING THE USE OF CHANNEL / OBJECTS
  !
  ! Author: Patrick Joeckel, MPICH, Jan 2009
  !
  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: DP

  ! ----------- <

  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modstr = 'submod2'
  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modver = '1.0'


! **************************************************************************
END MODULE messy_submod2
! **************************************************************************
