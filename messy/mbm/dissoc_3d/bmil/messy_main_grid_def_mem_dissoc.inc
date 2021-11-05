! -*- f90 -*-
  USE messy_main_constants_mem, ONLY: dp

  IMPLICIT NONE

  PUBLIC
  SAVE

  INTEGER, PARAMETER :: nlev = 10
  INTEGER, PARAMETER :: nlon = 36
  INTEGER, PARAMETER :: nlat = 18

  INTEGER :: nproma  = nlon   ! vector length
  INTEGER :: ngpblks = nlat   ! number of vector rows
  INTEGER :: npromz  = nlon   ! vector length of last row

  INTEGER :: kproma  = nlon   ! vector length of current row
  INTEGER :: jrow             ! row loop index
!  INTEGER :: je               !
