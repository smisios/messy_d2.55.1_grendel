MODULE messy_mbm_qbo_bl

  USE messy_main_constants_mem,   ONLY: dp
  USE messy_main_grid_def_mem_bi, ONLY: nlev

  IMPLICIT NONE

  PUBLIC
  SAVE

  CHARACTER(LEN=*), PARAMETER :: modstr = 'QBO'
  CHARACTER(LEN=*), PARAMETER :: modver = '1.1'

  REAL(dp) :: um1(1,1,nlev)    = 0._dp
  REAL(dp) :: vom_3d(1,1,nlev) = 0._dp
  ! SPECIFIC FOR QBO -------------------------------------------------------

!CONTAINS


END MODULE messy_mbm_qbo_bl
