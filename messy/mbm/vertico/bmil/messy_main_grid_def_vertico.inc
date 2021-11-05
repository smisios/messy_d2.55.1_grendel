! -*- f90 -*-

  USE messy_main_grid_def_mem_bi ! use without only as all parameters
                                 ! are set here
  USE messy_main_grid_def,      ONLY: submodstr
  USE messy_main_constants_mem, ONLY: dp
  USE messy_main_blather,       ONLY: start_message, end_message

  IMPLICIT NONE
  PUBLIC
  SAVE
  
  ! base model variables
  REAL(DP), POINTER, DIMENSION(:,:,:) :: grmass => NULL() ! grid mass
  REAL(DP), POINTER, DIMENSION(:,:,:) :: grmassdry => NULL() ! dry grid mass
  REAL(DP), POINTER, DIMENSION(:,:,:) :: grvol => NULL() ! grid volume
 
  REAL(DP), POINTER, DIMENSION(:,:,:) :: deltaz => NULL() ! layerthickness
  REAL(DP), POINTER, DIMENSION(:,:,:) :: altitude_gnd => NULL() ! height above ground
  REAL(DP), POINTER, DIMENSION(:,:,:) :: altitudei_gnd => NULL() ! height above ground
  ! variables for the onemis submodel
  REAL(DP), POINTER, DIMENSION(:,:)   :: philat_2d => NULL()     ! geographical latitude
  REAL(DP), POINTER, DIMENSION(:,:)   :: gboxarea_2d => NULL()   ! gridbox area
  REAL(DP), POINTER, DIMENSION(:,:)   :: philon_2d => NULL()     ! geographical longitude

CONTAINS

  !----------------------------------------------------------------------------
  SUBROUTINE main_grid_def_initialize

    IMPLICIT NONE

    nproma  = nlon
    ngpblks = nlat
    npromz  = nproma

    nlevp1 = nlev+1
    nvclev = nlevp1     ! number of levels with vertical coefficients
    ALLOCATE(vct(nvclev*2))
    vct(:) = 0.0_dp

  END SUBROUTINE main_grid_def_initialize
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE main_grid_def_init_memory

    ! MESSy
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         , new_attribute
    USE messy_main_channel_repr,     ONLY: get_representation_id

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_grid_def_init_memory'
    INTEGER :: status
    INTEGER :: reprid_scalar
    INTEGER :: reprid_int3d
    INTEGER :: reprid_mid3d
    INTEGER :: reprid_2d

    ! create new channel
    CALL start_message(submodstr,'CREATE STANDARD CHANNELS/OBJECTS',substr)
    ! get representations ids
    !
    CALL get_representation_id(status, 'GP_3D_MID', reprid_mid3d)
    CALL channel_halt(substr, status)
    !
    CALL get_representation_id(status, 'GP_3D_INT', reprid_int3d)
    CALL channel_halt(substr, status)
    !
    CALL get_representation_id(status, 'SCALAR', reprid_scalar)
    CALL channel_halt(substr, status)

    CALL get_representation_id(status, 'GP_2D_HORIZONTAL', reprid_2d)
    CALL channel_halt(substr, status)

    CALL new_channel(status, submodstr, reprid=reprid_scalar)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, submodstr, 'channel_info' &
         , c = 'standard basemodel channel' )
    CALL channel_halt(substr, status)


!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    !! PHYSICAL QUANTITIES
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! ------------------------------------------------------------------------
    !  grid mass in kg
    CALL new_channel_object(status, submodstr,  'grmass', &
         p3=grmass, reprid=reprid_mid3d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, submodstr, 'grmass', 'long_name', c='mass of gridbox')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, submodstr, 'grmass', 'units', c='kg')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    !  mass of dry air in kg
    CALL new_channel_object(status, submodstr,  'grmassdry', &
         p3=grmassdry, reprid=reprid_mid3d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, submodstr, 'grmassdry', 'long_name' &
         , c='mass of dry air')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, submodstr, 'grmassdry', 'units', c='kg')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    !  volume of grid box in m**3
    CALL new_channel_object(status, submodstr,  'grvol', &
         p3=grvol, reprid=reprid_mid3d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, submodstr, 'grvol', 'long_name', c='volume of gridbox')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, submodstr, 'grvol', 'units', c='m**3')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, submodstr,  'deltaz', p3=deltaz &
         , reprid=reprid_mid3d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, submodstr, 'deltaz', 'long_name', c='layerthickness')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, submodstr, 'deltaz', 'units', c='m')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, submodstr,  'altitude_gnd', p3=altitude_gnd &
         , reprid=reprid_mid3d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, submodstr, 'altitude_gnd', 'long_name', c='altitude above ground')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, submodstr, 'altitude_gnd', 'units', c='m')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, submodstr,  'altitudei_gnd', p3=altitudei_gnd &
         , reprid=reprid_mid3d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, submodstr, 'altitudei_gnd', 'long_name', c='altitude above ground at interfaces')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, submodstr, 'altitudei_gnd', 'units', c='m')
    CALL channel_halt(substr, status)

  END SUBROUTINE main_grid_def_init_memory
  !----------------------------------------------------------------------------

   SUBROUTINE main_grid_def_read_restart

     IMPLICIT NONE

     ! NOTHING TO DO HERE, AS VETICO IS NOT PARALLISED

   END SUBROUTINE main_grid_def_read_restart

  !----------------------------------------------------------------------------
  SUBROUTINE main_grid_def_global_start

    ! BML/MESSy
 
    IMPLICIT NONE

    ! LOCAL

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

  END SUBROUTINE main_grid_def_free_memory
  !----------------------------------------------------------------------------
