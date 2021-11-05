MODULE messy_main_grid_def_mem_bi

  USE messy_main_constants_mem, ONLY: dp
  INTEGER, PARAMETER :: nlon = 1
  INTEGER, PARAMETER :: nlat = 1
  INTEGER, PARAMETER :: ngl = 1

  INTEGER, PARAMETER :: nn = 42
  
  INTEGER, PARAMETER :: jrow = 1, kproma = 1, nproma=1
  INTEGER, PARAMETER :: ngpblks = 1

  INTEGER, PARAMETER :: nlev = 90          ! number of vertical levels
  LOGICAL            :: lmidatm = .TRUE.
  INTEGER, PARAMETER :: nlevp1 = nlev+1
  INTEGER, PARAMETER :: nvclev = nlevp1
  REAL(dp) :: vct(nvclev*2)
  REAL(dp) :: apzero = 1013.*1e2

#include "grid_def_mem_L90.inc"

 ! additional info see grid_bi:
 INTEGER, ALLOCATABLE, DIMENSION(:) :: BASEGRID_ID

END MODULE messy_main_grid_def_mem_bi
