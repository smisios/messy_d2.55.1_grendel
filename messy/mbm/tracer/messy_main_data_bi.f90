! **********************************************************************+
MODULE messy_main_data_bi
! **********************************************************************+

  ! THIS MODULE PROVIDES THE LINK BETWEEN THE BASMODEL (BML) AND THE BMIL

  USE messy_main_grid_def_mem_bi, ONLY: kproma, ngpblks, nlev
  USE messy_main_constants_mem,   ONLY: dp

  IMPLICIT NONE
  PUBLIC
 
  ! LENGTH OF THE INTEGRATION TIME STEP
  REAL(DP), PARAMETER :: time_step_len = 1.0_dp ! [s]

  ! THIS IS THE FLAG FOR THE RESTART MODE, WHICH IS NOT IMPLEMENTED IN
  ! THIS BOXMODEL
  LOGICAL, PARAMETER :: lresume = .FALSE.

  ! PHYSICAL DATA TO BE PROVIDED BY THE BASE MODEL (FOR TRACER_PDEF)
  REAL(DP), DIMENSION(kproma, nlev+1, ngpblks), SAVE :: pressi_3d    ! [Pa]

! **********************************************************************+
END MODULE messy_main_data_bi
! **********************************************************************+
