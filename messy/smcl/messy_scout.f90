! **********************************************************************
MODULE messy_scout
! **********************************************************************

  ! Selectable Column OUTput
  ! MODULE FOR HF-OUTPUT OF TRACERS/CHANNEL OBJECTS (COLUMNS)
  ! AT SELECTED GEOGRAPHIC LOCATIONS
  !
  ! CORE MODULE (MESSy/SMCL)
  !
  ! Author: Patrick Joeckel, MPICH, December 2003
  !

  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: DP

  ! ----------- <

  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modstr = 'scout'
  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modver = '2.1'

! **********************************************************************
END MODULE messy_scout
! **********************************************************************
