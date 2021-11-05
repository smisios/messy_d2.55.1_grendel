! -*- f90 -*-

  USE messy_main_constants_mem, ONLY: dp

  IMPLICIT NONE
  PUBLIC
  SAVE

  INTEGER, PARAMETER :: nlev = 90
  INTEGER, PARAMETER :: nlon =  1
  INTEGER, PARAMETER :: nlat =  1

  INTEGER :: nproma  = nlon
  INTEGER :: npromz  = nlon
  INTEGER :: ngpblks = nlat

  INTEGER :: jrow
  INTEGER :: kproma  = nlon

  REAL(dp) :: apzero         = 1._dp
