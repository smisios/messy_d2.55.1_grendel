! -*- f90 -*-

  ! GRID DEFINED BY BASEMODEL
  USE messy_main_grid_def_mem_bi ! use without only
  USE messy_main_grid_def,      ONLY: submodstr
  USE messy_main_constants_mem, ONLY: dp
  USE messy_main_blather,       ONLY: start_message, end_message

  IMPLICIT NONE
  PUBLIC
  SAVE

  !  LOCALIZED PARAMETERS
  ! - GRID CONTROL HERE FIXED VALUES normally determined by grid definition

  REAL(DP), ALLOCATABLE, DIMENSION(:) :: longitude
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: latitude
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: level
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: hyai
  !
  REAL(DP), DIMENSION(:,:), POINTER :: sinlat_2d => NULL()
  REAL(DP), DIMENSION(:,:), POINTER :: coslat_2d => NULL()
  REAL(DP), DIMENSION(:,:), POINTER :: sinlon_2d => NULL()
  REAL(DP), DIMENSION(:,:), POINTER :: coslon_2d => NULL()
  ! 
  REAL(dp), POINTER, DIMENSION(:)    :: ceta       => NULL()
  ! 
CONTAINS

  ! --------------------------------------------------------------------------
  SUBROUTINE main_grid_def_setup

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
    DO i=nlev,1,-1
       hyai(i) = hyai(i+1) /2.
    END DO

  END SUBROUTINE main_grid_def_setup
  ! --------------------------------------------------------------------------

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
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         , new_attribute
    USE messy_main_channel_repr,     ONLY: get_representation_id
    USE messy_main_constants_mem,    ONLY: DTR

    IMPLICIT NONE
    INTRINSIC :: COS, SIN

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_grid_def_init_memory'
    INTEGER :: status
    INTEGER :: reprid
    INTEGER :: i, j

    ! create new channel
    CALL start_message(submodstr,'CREATE STANDARD CHANNELS/OBJECTS',substr)

    CALL get_representation_id(status, 'SCALAR', reprid)
    CALL channel_halt(substr, status)

    CALL new_channel(status, submodstr, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, submodstr, 'channel_info' &
         , c = 'standard basemodel channel' )
    CALL channel_halt(substr, status)

    ! ------------------------------------------------------------------------
    ! ------------------------------------------------------------------------
    CALL get_representation_id(status, 'GP_2D_HORIZONTAL', reprid)
    ! ------------------------------------------------------------------------

    ! ------------------------------------------------------------------------
    ! COS(longitude)
    CALL new_channel_object(status, submodstr,  'coslon', &
         p2=coslon_2d, reprid=reprid, lstatic=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, submodstr, 'coslon', &
         'long_name', c='cos(longitude)')
    CALL channel_halt(substr, status)
    DO i = 1, nproma
       DO j = 1, ngpblks
          coslon_2d(i,j) = COS( longitude(i) * DTR )
       ENDDO
    ENDDO
    ! ------------------------------------------------------------------------
    ! ------------------------------------------------------------------------
    ! SIN(longitude)
    CALL new_channel_object(status, submodstr,  'sinlon', &
         p2=sinlon_2d, reprid=reprid, lstatic=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, submodstr, 'sinlon', &
         'long_name', c='sin(longitude)')
    CALL channel_halt(substr, status)
    DO i = 1, nproma
       DO j = 1, ngpblks
          sinlon_2d(i,j) = SIN( longitude(i) * DTR )
       ENDDO
    ENDDO
    ! ------------------------------------------------------------------------
    ! ------------------------------------------------------------------------
    ! COS(latitude)
    CALL new_channel_object(status, submodstr,  'coslat', &
         p2=coslat_2d, reprid=reprid, lstatic=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, submodstr, 'coslat', &
         'long_name', c='cos(latitude)')
    CALL channel_halt(substr, status)
    DO i = 1, nproma
       DO j = 1, ngpblks
          coslat_2d(i,j) = COS( latitude(j) * DTR )
       ENDDO
    ENDDO
    ! ------------------------------------------------------------------------
    ! ------------------------------------------------------------------------
    ! SIN(latitude)
    CALL new_channel_object(status, submodstr,  'sinlat', &
         p2=sinlat_2d, reprid=reprid, lstatic=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, submodstr, 'sinlat', &
         'long_name', c='sin(latitude)')
    CALL channel_halt(substr, status)
    DO i = 1, nproma
       DO j = 1, ngpblks
          sinlat_2d(i,j) = SIN( latitude(j) * DTR )
       ENDDO
    ENDDO
    ! ------------------------------------------------------------------------

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

    DEALLOCATE(latitude)
    DEALLOCATE(longitude)
    DEALLOCATE(level)
    DEALLOCATE(hyai)

  END SUBROUTINE main_grid_def_free_memory
  !----------------------------------------------------------------------------
