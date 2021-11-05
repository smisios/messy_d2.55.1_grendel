! **********************************************************************
MODULE channel_mem_bml
! **********************************************************************

  USE messy_main_constants_mem, ONLY: dp
  IMPLICIT NONE
  SAVE

  ! NAME AND VERSION OF THE BASEMODEL
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'CHANNEL'
  CHARACTER(len=*), PARAMETER, PUBLIC :: modver = '1.0'

  ! LENGTH OF TIME STEP
  REAL(DP), PARAMETER :: time_step_len = 1.0_dp ! hours

  ! RESTART FREQUENCY [time steps]
  INTEGER, PARAMETER :: i_restfreq = 10

  ! NUMBER OF TIME STEPS TO RUN
  INTEGER, PARAMETER :: NH_max = 9999

  ! TIME STEP COUNTER
  INTEGER :: NH = 0

  ! RESTART CYCLE COUNTER
  INTEGER          :: cy = 0
  CHARACTER(LEN=4) :: cystr = '0000'

  ! RESTART MODE ?
  LOGICAL :: lresume = .FALSE.

! **********************************************************************
END MODULE channel_mem_bml
! **********************************************************************
