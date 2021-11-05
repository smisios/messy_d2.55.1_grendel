! -*- f90 -*-
 USE messy_main_constants_mem, ONLY: dp
 USE messy_main_grid_def,      ONLY: submodstr

  IMPLICIT NONE

  PUBLIC
  SAVE
 
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: longitude
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: latitude
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: level
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: hyai

CONTAINS
  ! --------------------------------------------------------------------------
  SUBROUTINE main_grid_def_setup

    USE messy_main_grid_def_mem_bi, ONLY: nlon, nlat, nlev

    IMPLICIT NONE

    INTEGER :: i

    ! DEFINE GRID
    ! LONGITUDE
    ALLOCATE(longitude(nlon))
    DO i=1, nlon
       longitude(i) = REAL((i-1)*10, DP) + 5.0_dp  ! degrees_east
    END DO
    ! LATITUDE
    ALLOCATE(latitude(nlat))
    DO i=1, nlat
       latitude(i) = -85.0_dp + REAL((i-1)*10, DP) ! degrees_north
    END DO
    ! LEVELS
    ALLOCATE(level(nlev))
    DO i=1, nlev
       level(i) = REAL(i, DP)
    END DO
    ALLOCATE(hyai(nlev+1))
    hyai(nlev+1) = 101325.0
    DO i= nlev,1 ,-1
       hyai(i) = hyai(i+1) / 2. 
    END DO

  END SUBROUTINE main_grid_def_setup
  ! --------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE main_grid_def_initialize

    IMPLICIT NONE

  END SUBROUTINE main_grid_def_initialize
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE main_grid_def_init_memory

    ! MESSy
!!$    USE messy_main_channel_error_bi, ONLY: channel_halt
!!$    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
!!$                                         , new_attribute
!!$    USE messy_main_channel_repr,     ONLY: get_representation_id
!!$    USE messy_main_constants_mem,    ONLY: DTR
!!$    USE messy_main_blather,          ONLY: start_message, end_message

    IMPLICIT NONE
 !!$   INTRINSIC :: COS, SIN, REAL

    ! LOCAL
!!$    CHARACTER(LEN=*), PARAMETER :: substr = 'main_grid_def_init_memory'
!!$    INTEGER :: status
!!$    INTEGER :: reprid

    ! create new channel
!!$    CALL start_message(submodstr,'CREATE STANDARD CHANNELS/OBJECTS',substr)
!!$
!!$    CALL get_representation_id(status, 'SCALAR', reprid)
!!$    CALL channel_halt(substr, status)
!!$
!!$    CALL new_channel(status, submodstr, reprid=reprid)
!!$    CALL channel_halt(substr, status)
!!$    CALL new_attribute(status, submodstr, 'channel_info' &
!!$         , c = 'standard basemodel channel' )
!!$    CALL channel_halt(substr, status)

   ! ------------------------------------------------------------------------

!!$    CALL end_message(submodstr,'CREATE STANDARD CHANNELS/OBJECTS',substr)

  END SUBROUTINE main_grid_def_init_memory
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE main_grid_def_read_restart

    IMPLICIT NONE

    ! NOTHING TO DO

  END SUBROUTINE main_grid_def_read_restart
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE main_grid_def_global_start

    IMPLICIT NONE

  END SUBROUTINE main_grid_def_global_start
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE main_grid_def_local_start

    IMPLICIT NONE

    ! jrow is set in main program loop (region loop)

  END SUBROUTINE main_grid_def_local_start
  !----------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE main_grid_def_free_memory

    IMPLICIT NONE

    DEALLOCATE(latitude)
    DEALLOCATE(longitude)
    DEALLOCATE(level)

  END SUBROUTINE main_grid_def_free_memory
  ! --------------------------------------------------------------------------
