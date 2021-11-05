! **********************************************************************+
MODULE messy_main_grid_def_bi
! **********************************************************************+

  ! THIS MODULE PROVIDES THE LINK BETWEEN THE BASMODEL (BML) AND THE BMIL

  USE messy_main_grid_def_mem_bi, ONLY: kproma, ngpblks, nlev
  USE messy_main_constants_mem,   ONLY: dp

  IMPLICIT NONE
  PUBLIC
 
  REAL(DP), DIMENSION(kproma, ngpblks),         SAVE :: gboxarea_2d  ! [m^2]
  REAL(DP), DIMENSION(kproma, nlev, ngpblks),   SAVE :: grmass       ! [kg]

! **********************************************************************+
END MODULE messy_main_grid_def_bi
! **********************************************************************+
