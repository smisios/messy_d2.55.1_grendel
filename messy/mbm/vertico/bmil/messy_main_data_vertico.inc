! -*- f90 -*-

  USE messy_main_constants_mem, ONLY: dp
  USE messy_main_blather,       ONLY: start_message, end_message

  IMPLICIT NONE
  PUBLIC
  SAVE

  ! NAME AND VERSION OF THE BASEMODEL
  CHARACTER(LEN=*), PARAMETER :: modstr = 'VERTICO'
  CHARACTER(LEN=*), PARAMETER :: modver = '1.0'

  ! dummy
  LOGICAL, PARAMETER :: l2tls = .FALSE.
  LOGICAL, PARAMETER :: L_IS_CHILD = .FALSE.
  LOGICAL            :: ledith  = .FALSE.

  ! defined constants
  REAL(dp) :: eps  = 0.1_dp

  ! exemplary variables for standard basemodel
  REAL(DP), POINTER :: yr,mo,dy,hr,mi,se,ms

   ! base model variables
   REAL(DP), POINTER, DIMENSION(:,:,:) :: pressi_3d => NULL() ! press at interface
   REAL(DP), POINTER, DIMENSION(:,:,:) :: press_3d => NULL() ! press at midpoint
   REAL(DP), POINTER, DIMENSION(:,:,:) :: geopoti_3d => NULL() ! geopotential at interface
   REAL(DP), POINTER, DIMENSION(:,:,:) :: geopot_3d => NULL() ! geopotential

   REAL(DP), POINTER, DIMENSION(:,:)   :: az0 => NULL()       ! roughness length

   ! variables for the ddep submodel
   REAL(DP), POINTER, DIMENSION(:,:,:) :: tm1 => NULL()     ! temperature
   REAL(DP), POINTER, DIMENSION(:,:)   :: cvs => NULL()       ! snow cover
   REAL(DP), POINTER, DIMENSION(:,:)   :: cvw => NULL()       ! wet skin fraction
   REAL(DP), POINTER, DIMENSION(:,:)   :: alake => NULL()     ! lake fraction
   REAL(DP), POINTER, DIMENSION(:,:)   :: slm => NULL()       ! land mask
   REAL(DP), POINTER, DIMENSION(:,:)   :: seaice => NULL()    ! sea ice fraction
   REAL(DP), POINTER, DIMENSION(:,:)   :: vgrat => NULL()     ! vegetation fraction
   REAL(DP), POINTER, DIMENSION(:,:)   :: cdnl => NULL()      ! neutral drag coeff., land
   REAL(DP), POINTER, DIMENSION(:,:)   :: cfncl => NULL()     ! exchange parameter, land
   REAL(DP), POINTER, DIMENSION(:,:)   :: cfml => NULL()      ! momentum drag coeff., land
   REAL(DP), POINTER, DIMENSION(:,:,:) :: um1 => NULL()       ! horizontal wind velocity
   REAL(DP), POINTER, DIMENSION(:,:,:) :: vm1 => NULL()       ! horizontal wind velocity
   REAL(DP), POINTER, DIMENSION(:,:)   :: tvir => NULL()      ! surface virtual temperature
   REAL(DP), POINTER, DIMENSION(:,:)   :: tvl => NULL()       ! surface virtual temperature (land)
   REAL(DP), POINTER, DIMENSION(:,:)   :: ril => NULL()       ! Richardson number (land)
   REAL(DP), POINTER, DIMENSION(:,:)   :: cdnw => NULL()      ! neutral drag coeff., water
   REAL(DP), POINTER, DIMENSION(:,:)   :: cfmw => NULL()      ! momentum drag coeff., water
   REAL(DP), POINTER, DIMENSION(:,:)   :: cfncw => NULL()     ! exchange parameter, water
   REAL(DP), POINTER, DIMENSION(:,:)   :: riw => NULL()       ! Richardson number (water)
   REAL(DP), POINTER, DIMENSION(:,:)   :: tvw => NULL()       ! surface virtal temperature (water)
   REAL(DP), POINTER, DIMENSION(:,:)   :: cdni => NULL()      ! neutral drag coeff., ice
   REAL(DP), POINTER, DIMENSION(:,:)   :: cfmi => NULL()      ! momentum drag coeff., ice
   REAL(DP), POINTER, DIMENSION(:,:)   :: cfnci => NULL()     ! exchange parameter, ice
   REAL(DP), POINTER, DIMENSION(:,:)   :: rii => NULL()       ! Richardson number (ice)
   REAL(DP), POINTER, DIMENSION(:,:)   :: tvi => NULL()       ! surface virtual temperature (ice)
   REAL(DP), POINTER, DIMENSION(:,:)   :: tslm1 => NULL()     ! ground surface temperature
   REAL(DP), POINTER, DIMENSION(:,:)   :: rh_2m => NULL()     ! relative humidity at 2m
   REAL(DP), POINTER, DIMENSION(:,:)   :: rco_leaf => NULL()  ! leaf stomatal resistance
   REAL(DP), POINTER, DIMENSION(:,:)   :: fws => NULL()       ! soil moisture stress function
   REAL(DP), POINTER, DIMENSION(:,:)   :: u10 => NULL()       ! 10m u-velocity
   REAL(DP), POINTER, DIMENSION(:,:)   :: v10 => NULL()       ! 10m v-velocity

   !mz_ho_20160412+
   ! variables for the vdiff/vertex submodel
   REAL(DP), POINTER, DIMENSION(:,:,:) :: xlm1 => NULL()      ! cloud water
   REAL(DP), POINTER, DIMENSION(:,:,:) :: xim1 => NULL()      ! cloud ice
   !mz_ho_20160412-

   ! variables for the megan submodel
   REAL(DP), POINTER, DIMENSION(:,:)   :: slf => NULL()     ! sea-land fraction

   ! variables for the onemis submodel
   REAL(DP), POINTER, DIMENSION(:,:)   :: tsw => NULL()           ! surface temperature of water
   REAL(DP), POINTER, DIMENSION(:,:,:) :: rho_air_dry_3d => NULL()! density of dry air
   REAL(DP), POINTER, DIMENSION(:,:)   :: wind10_2d => NULL()     ! 10 m wind speed
   REAL(DP), POINTER, DIMENSION(:,:)   :: srfl => NULL()          ! net surface radiative flux
   REAL(DP), POINTER, DIMENSION(:,:)   :: prc => NULL()           ! convective precipitation
   REAL(DP), POINTER, DIMENSION(:,:)   :: prl => NULL()           ! large-scale precipitation
   REAL(DP), POINTER, DIMENSION(:,:)   :: wsoil => NULL()            ! soil moisture
   REAL(DP), POINTER, DIMENSION(:,:,:) :: tsoil => NULL()         ! deep soil temperatures
   REAL(DP), POINTER, DIMENSION(:,:)   :: glac => NULL()          ! fraction of land covered by glaciers
   REAL(DP), POINTER, DIMENSION(:,:)   :: wsmx => NULL()          ! soil moisture at field capacity
   REAL(DP), POINTER, DIMENSION(:,:)   :: zust_2d => NULL()       ! surface friction velocity
   REAL(DP), POINTER, DIMENSION(:,:,:) :: tm1_3d => NULL()        ! air temperature
   REAL(DP), POINTER, DIMENSION(:,:,:) :: qm1_3d => NULL()        ! specific humidity
   REAL(DP), POINTER, DIMENSION(:,:) :: tsurf_2d => NULL()        ! surface temperature

   ! variables for the jval submodel
   REAL(DP), POINTER, DIMENSION(:,:,:) :: tte_3d   => NULL()     ! air temperature tendency
   REAL(DP), POINTER, DIMENSION(:,:,:) :: rhum_3d  => NULL()     ! rel. hum.
   REAL(DP), POINTER, DIMENSION(:,:,:) :: xlm1_3d  => NULL()     ! cloud water
   REAL(DP), POINTER, DIMENSION(:,:,:) :: xlte_3d  => NULL()     ! cloud water tendency
   REAL(DP), POINTER, DIMENSION(:,:,:) :: xim1_3d  => NULL()     ! cloud ice
   REAL(DP), POINTER, DIMENSION(:,:,:) :: xite_3d  => NULL()     ! cloud ice tendency
   REAL(DP), POINTER, DIMENSION(:,:)   :: albedo   => NULL()     ! surface albedo
   REAL(DP), POINTER, DIMENSION(:,:,:) :: aclc     => NULL()     ! 3D(!?) cloud cover

   REAL(DP), POINTER, DIMENSION(:,:,:) :: o3h   => NULL()        ! column relative ozone
   REAL(DP), POINTER, DIMENSION(:,:,:) :: v3h   => NULL()        ! vertical ozone column
   REAL(DP), POINTER, DIMENSION(:,:,:) :: pressh => NULL()       ! column pressure

   ! variables for the mecca submodel
   REAL(DP), POINTER, DIMENSION(:,:,:) :: qte_3d   => NULL()     ! spec. humidity  tendency
   REAL(DP), POINTER, DIMENSION(:,:,:) :: qm1      => NULL()     ! spec. humidity

CONTAINS

  !----------------------------------------------------------------------------
  SUBROUTINE main_data_initialize

    USE messy_main_timer,        ONLY: timer_set_time_step_len

    IMPLICIT NONE

    CALL timer_set_time_step_len(l2tls)

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
    INTEGER :: reprid_scalar
    INTEGER :: reprid_int3d
    INTEGER :: reprid_mid3d
    INTEGER :: reprid_2d

    ! create new channel
    CALL start_message(modstr,'CREATE STANDARD CHANNELS/OBJECTS',substr)
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

    CALL new_channel(status, modstr, reprid=reprid_scalar)
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

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    !! PHYSICAL QUANTITIES
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! pressure at middle of box ("full level pressure")
    CALL new_channel_object(status, modstr,  'press', &
         p3=press_3d, reprid= reprid_mid3d, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'press', &
         'long_name', c='pressure')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'press', 'units', c='Pa')
    CALL channel_halt(substr, status)

    ! pressure at level interfaces ("half level pressure")
    CALL new_channel_object(status, modstr,  'pressi', &
         p3=pressi_3d, reprid= reprid_int3d, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'pressi', &
         'long_name', c='pressure at interface')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'pressi', 'units', c='Pa')
    CALL channel_halt(substr, status)

    ! geopotential at middle of box ("full level geopotential")
    CALL new_channel_object(status, modstr,  'geopot', &
         p3=geopot_3d, reprid= reprid_mid3d, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'geopot', &
         'long_name', c='geopotential')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'geopot', 'units', c='m2 s-2')
    CALL channel_halt(substr, status)

    ! geopotential at level interface ("half level geopotential")
    CALL new_channel_object(status, modstr,  'geopoti', &
         p3=geopoti_3d, reprid= reprid_int3d, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'geopoti', &
         'long_name', c='geopotential at interface')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'geopoti', 'units', c='m2 s-2')
    CALL channel_halt(substr, status)

    ! ------------------------------------------------------------------------
    CALL new_channel_object(status, modstr,  'az0', &
         p2=az0, reprid=reprid_2d, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'az0', &
         'long_name', c='roughness length')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'az0', 'units', c='m')
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
    READ(ydate,'(3x,i2)') imo
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

