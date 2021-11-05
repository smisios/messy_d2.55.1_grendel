! -*- f90 -*-

  USE messy_main_constants_mem, ONLY: dp
  USE messy_main_blather,       ONLY: start_message, end_message

  IMPLICIT NONE
  PUBLIC
  SAVE

  ! NAME AND VERSION OF THE BASEMODEL
  CHARACTER(LEN=*), PARAMETER :: modstr = 'MBM_MPIOM'
  CHARACTER(LEN=*), PARAMETER :: modver = '1.0'

  LOGICAL, PARAMETER :: l2tls = .FALSE.
  ! LOCALIZED PARAMETERS
  ! - GRID CONTROL HERE FIXED VALUES normally determined by grid definition

  logical :: lcouple
  ! dummy
  LOGICAL, PARAMETER :: L_IS_CHILD = .FALSE.

  ! defined constants
  REAL(dp) :: eps  = 0.1_dp

  ! exemplary variables for standard basemodel
  REAL(DP), POINTER :: yr,mo,dy,hr,mi,se,ms

CONTAINS

  !----------------------------------------------------------------------------
  SUBROUTINE main_data_initialize

    IMPLICIT NONE

  END SUBROUTINE main_data_initialize
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE main_data_init_memory

    ! MESSy
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         , new_attribute
    USE messy_main_channel_repr,     ONLY: get_representation_id

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_data_init_memory'
    INTEGER :: status
    INTEGER :: reprid

    ! create new channel
    CALL start_message(modstr,'CREATE STANDARD CHANNELS/OBJECTS',substr)

    CALL get_representation_id(status, 'SCALAR', reprid)
    CALL channel_halt(substr, status)

    CALL new_channel(status, modstr, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'channel_info' &
         , c = 'standard basemodel channel' )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'yr', yr)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'yr' &
         , 'long_name', c = 'year')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'mo', mo)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'mo' &
         , 'long_name', c = 'month')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'dy', dy)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dy' &
         , 'long_name', c = 'day')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'hr', hr)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'hr' &
         , 'long_name', c = 'hour')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'mi', mi)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'mi' &
         , 'long_name', c = 'minute')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'se', se)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'se' &
         , 'long_name', c = 'second')
    CALL channel_halt(substr, status)

    ! FORCE RESTART FILE CREATION
    CALL new_channel_object(status, modstr, 'ms', ms, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ms' &
         , 'long_name', c = 'millisecond')
    CALL channel_halt(substr, status)

    CALL end_message(modstr,'CREATE STANDARD CHANNELS/OBJECTS',substr)

  END SUBROUTINE main_data_init_memory
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE main_data_global_start

    IMPLICIT NONE

    INTRINSIC :: DATE_AND_TIME, TRIM, ADJUSTL

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_data_global_start'
    CHARACTER (8)   :: ydate
    CHARACTER (10)  :: ytime

    INTEGER :: iyr, imo, idy, ihr, imi, ise, ims

    CALL DATE_AND_TIME(ydate, ytime)

    READ(ydate,'(i4)') iyr
    READ(ydate,'(4x,i2)') imo
    READ(ydate,'(6x,i2)') idy

    READ(ytime,'(i2)') ihr
    READ(ytime,'(2x,i2)') imi
    READ(ytime,'(4x,i2)') ise
    READ(ytime,'(7x,i3)') ims

    yr = REAL(iyr,dp)
    mo = REAL(imo,dp)
    dy = REAL(idy,dp)
    hr = REAL(ihr,dp)
    mi = REAL(imi,dp)
    se = REAL(ise,dp)
    ms = REAL(ims,dp)

  END SUBROUTINE main_data_global_start
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE main_data_local_start

    IMPLICIT NONE

    ! jrow is set in main program loop (region loop)

  END SUBROUTINE main_data_local_start
  !----------------------------------------------------------------------------
