! -*- f90 -*-

  USE messy_main_grid_def_mem_bi ! USE without ONLY
  USE messy_main_grid_def,      ONLY: submodstr
  USE messy_main_constants_mem, ONLY: dp
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi

  IMPLICIT NONE
  PUBLIC
  SAVE

CONTAINS

  !----------------------------------------------------------------------------
  SUBROUTINE main_grid_def_initialize

    IMPLICIT NONE

    nproma  = nlon
    ngpblks = nlat
    npromz  = nproma

  END SUBROUTINE main_grid_def_initialize
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE main_grid_def_init_memory

    ! MESSy
!!$    USE messy_main_channel_error_bi, ONLY: channel_halt
!!$    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
!!$                                         , new_attribute
!!$    USE messy_main_channel_repr,     ONLY: get_representation_id

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_grid_def_init_memory'

    ! create new channel
    CALL start_message_bi(submodstr,'CREATE STANDARD CHANNELS/OBJECTS',substr)

    CALL end_message_bi(submodstr,'CREATE STANDARD CHANNELS/OBJECTS',substr)

  END SUBROUTINE main_grid_def_init_memory
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE main_grid_def_read_restart

 
    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_grid_def_read_restart'

  END SUBROUTINE main_grid_def_read_restart
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE main_grid_def_global_start

 
    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_grid_def_global_start'

  END SUBROUTINE main_grid_def_global_start
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE main_grid_def_local_start

    IMPLICIT NONE

    ! jrow is set in main program loop (region loop)

    IF ( jrow == ngpblks ) THEN
       kproma = npromz
    ELSE
       kproma = nproma
    END IF

  END SUBROUTINE main_grid_def_local_start
  !----------------------------------------------------------------------------
 
 !----------------------------------------------------------------------------
  SUBROUTINE main_grid_def_free_memory

 
    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_grid_def_free_memory'

  END SUBROUTINE main_grid_def_free_memory
  !----------------------------------------------------------------------------

