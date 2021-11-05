! **********************************************************************+
MODULE messy_main_grid_def_mem_bi
! **********************************************************************+

  ! THIS MODULE PROVIDES THE LINK BETWEEN THE BASMODEL (BML) AND THE BMIL

  USE messy_main_constants_mem, ONLY: dp

  IMPLICIT NONE
  PUBLIC
 
  ! THESE ARE THE TRACER SET DATA DIMENSIONS, E.G. THE GRID
  INTEGER, PARAMETER :: nproma  = 5
  INTEGER, PARAMETER :: npromz  = 5
  INTEGER, PARAMETER :: nlev    = 2
  INTEGER, PARAMETER :: ngpblks = 10
  !
  INTEGER, PARAMETER :: kproma  = nproma

  ! THIS IS THE LOOP VARIABLE FOR THE 'LOCAL LOOP' ALONG THE 3rd DIMENSION
  INTEGER, SAVE :: jrow

! **********************************************************************+
END MODULE messy_main_grid_def_mem_bi
! **********************************************************************+
