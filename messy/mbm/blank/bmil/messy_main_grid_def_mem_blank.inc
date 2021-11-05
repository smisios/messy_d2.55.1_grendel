! -*- f90 -*-
  ! GRID DEFINED BY BASEMODEL
  USE messy_main_constants_mem, ONLY: dp

  IMPLICIT NONE
  PUBLIC
  SAVE

  INTEGER, PARAMETER :: nlev =  2
  INTEGER, PARAMETER :: nlon = 36
  INTEGER, PARAMETER :: nlat = 18

  INTEGER :: nproma           ! vector length
  INTEGER :: ngpblks          ! number of vector rows
  INTEGER :: npromz           ! vector length of last row

  INTEGER :: kproma           ! vector length of current row
  INTEGER :: jrow             ! row loop index
  INTEGER :: je               !

#ifdef MESSYDWARF
  INTEGER :: nlevp1 = nlev + 1
  INTEGER :: nlevm1 = nlev - 1
  INTEGER :: nvclev = nlev + 1
  INTEGER :: nn = -50

  ! imitate ECHAM5
  LOGICAL :: lmidatm = .FALSE.
#endif
