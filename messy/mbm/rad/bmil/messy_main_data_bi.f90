!*****************************************************************************
#include "messy_main_ppd_bi.inc"

MODULE messy_main_data_bi
!*****************************************************************************

  USE messy_main_grid_def_mem_bi, ONLY: nlon, nlat, nlev, nlevp1 &
                                      , lv_echam, ladd_tte
  USE messy_main_constants_mem, ONLY : dp &
                                     , dtr    ! Degrees to radians
  USE messy_main_blather,       ONLY : start_message, end_message
  USE messy_main_timer,         ONLY : time_step_len, delta_time

  IMPLICIT NONE
  PUBLIC
  SAVE

  ! NAME AND VERSION OF THE BASEMODEL
  CHARACTER(LEN=*), PARAMETER :: modstr = 'MBM_RAD'
  CHARACTER(LEN=*), PARAMETER :: modver = '2.0'

  LOGICAL, PARAMETER :: l2tls = .FALSE.      
  ! LOCALIZED PARAMETERS
  ! - GRID CONTROL HERE FIXED VALUES normally determined by grid definition
  REAL(dp), PARAMETER :: eps    = 0.1_dp
  ! 
  ! - 2D - fields boundary conditions
  !
  REAL(dp), POINTER, DIMENSION(:,:) :: slf        => NULL() ! sea land fraction
  REAL(dp), POINTER, DIMENSION(:,:) :: alb        => NULL() ! surface background albedo
  REAL(dp), POINTER, DIMENSION(:,:) :: forest     => NULL() ! forest fraction
  REAL(dp), POINTER, DIMENSION(:,:) :: glac       => NULL() ! fraction of land covered by glaciers
  REAL(dp), POINTER, DIMENSION(:,:) :: slm        => NULL() ! land mask
  REAL(dp), POINTER, DIMENSION(:,:) :: geosp      => NULL() ! surface geopotential
  REAL(dp), POINTER, DIMENSION(:,:) :: vlt        => NULL() ! leaf area index (LAI)
  REAL(dp), POINTER, DIMENSION(:,:) :: seaice     => NULL() ! ice cover (fraction of 1-SLM)
  REAL(dp), POINTER, DIMENSION(:,:) :: icecov     => NULL() ! ice cover fraction
  REAL(dp), POINTER, DIMENSION(:,:) :: seacov     => NULL() ! sea cover fraction
  REAL(dp), POINTER, DIMENSION(:,:) :: sni        => NULL() ! water equivalent of snow on ice
  REAL(dp), POINTER, DIMENSION(:,:) :: aps        => NULL() ! surface pressure
  REAL(dp), POINTER, DIMENSION(:,:) :: cvs        => NULL() ! snow cover
  REAL(dp), POINTER, DIMENSION(:,:) :: cvsc       => NULL() ! snow covered canopy
  REAL(dp), POINTER, DIMENSION(:,:) :: tslm1      => NULL() ! land surface temperature
  REAL(dp), POINTER, DIMENSION(:,:) :: tslnew     => NULL() ! land surface temperature for sensible heat flux
  REAL(dp), POINTER, DIMENSION(:,:) :: tsi        => NULL() ! surface temperature of ice
  REAL(dp), POINTER, DIMENSION(:,:) :: tsw        => NULL() ! surface temperature of water
  ! 
  ! - 2D - fields surface diagnostics
  !
  REAL(dp), POINTER, DIMENSION(:,:) :: radflxw_2d => NULL() ! radiation flux (lw + sw) over water
  REAL(dp), POINTER, DIMENSION(:,:) :: srfl_2d    => NULL() ! 
  !
  ! - 3D - fields used for radiation coupling
  !
  REAL(dp), POINTER, DIMENSION(:,:,:) :: tm1      => NULL() ! temperature
  REAL(dp), POINTER, DIMENSION(:,:,:) :: tm1_imp  => NULL() ! temperature (imported)
  REAL(dp), POINTER, DIMENSION(:,:,:) :: qm1      => NULL() ! specific humidity
  ! 
  ! - 3D - fields
  !
  REAL(dp), POINTER, DIMENSION(:,:,:) :: press_3d => NULL() ! pressure at full model levels
  REAL(dp), POINTER, DIMENSION(:,:,:) :: pressi_3d=> NULL() ! pressure at half model levels
  REAL(dp), POINTER, DIMENSION(:,:,:) :: geopoti_3d=> NULL() ! interface geopotential in m2 s-2'
  REAL(dp), POINTER, DIMENSION(:,:,:) :: tvirt    => NULL() ! virtual temperature
  
  REAL(dp), ALLOCATABLE, DIMENSION(:,:)   :: aphm1  ! pressure at half level
  REAL(dp), ALLOCATABLE, DIMENSION(:,:)   :: apm1   ! pressure at full level
  REAL(dp), ALLOCATABLE, DIMENSION(:,:)   :: app1   ! pressure at full level
  REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: tte_3d ! temperature tendency
  
  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: loland_2d ! logical lan
  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: loglac_2d ! logical glacier
  
  LOGICAL, PARAMETER :: L_IS_CHILD = .FALSE.
  LOGICAL, PARAMETER :: lcouple = .FALSE.

  ! extension in ionosphere/thermosphere ?
  LOGICAL :: ledith  = .FALSE.

CONTAINS

  !----------------------------------------------------------------------------
  SUBROUTINE main_data_initialize

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_data_initialize'

    CALL start_message(modstr,'INITIALISATION',substr)
    
    ! initialise 2D arrays
    ! 
    ALLOCATE( &
           loland_2d(nlon,nlat), loglac_2d(nlon,nlat) &
         , aphm1(nlon,nlevp1), apm1(nlon,nlev), app1(nlon,nlev) &
         )
    loland_2d(:,:) = .FALSE.
    loglac_2d(:,:) = .FALSE.
    apm1(:,:)  = 0._dp
    aphm1(:,:) = 0._dp
    app1(:,:)  = 0._dp
    
    ! initialise 3D arrays
    ! 
    ALLOCATE( &
          tte_3d(_RI_XYZ__(nlon,nlat,nlev)) &
         )
    tte_3d(:,:,:) = 0._dp
    !
    CALL end_message(modstr,'INITIALISATION',substr)
    !
  END SUBROUTINE main_data_initialize
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE main_data_init_memory

    ! MESSy/BMIL
    USE messy_main_channel_error_bi, ONLY: channel_halt
    ! MESSy/SMCL
    USE messy_main_channel,       ONLY: new_channel, new_channel_object &
                                      , new_attribute
    USE messy_main_channel_repr,  ONLY: get_representation_id

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_data_init_memory'
    INTEGER :: status
    INTEGER :: GP_3D_MID_ID, GP_3D_INT_ID, GP_2D_HORIZONTAL_ID
    
    ! create new channel
    CALL start_message(modstr,'CREATE STANDARD CHANNELS/OBJECTS',substr)

    CALL get_representation_id(status, 'GP_3D_MID', reprid=GP_3D_MID_ID)
    CALL channel_halt(substr, status)
    
    CALL get_representation_id(status, 'GP_3D_INT', reprid=GP_3D_INT_ID)
    CALL channel_halt(substr, status)
    
    CALL get_representation_id(status, 'GP_2D_HORIZONTAL', reprid=GP_2D_HORIZONTAL_ID)
    CALL channel_halt(substr, status)
    
    CALL new_channel(status, modstr)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'channel_info', c = 'standard basemodel channel' )
    CALL channel_halt(substr, status)
    !
    ! 2D fields
    !
    CALL new_channel_object(status, modstr, 'radflxw_2d', p2=radflxw_2d, reprid=GP_2D_HORIZONTAL_ID)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'radflxw_2d', 'long_name', c = 'radiation flux (lw + sw) over water')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'radflxw_2d', 'units', c = 'W/m^2')
    CALL channel_halt(substr, status)
            
    CALL new_channel_object(status, modstr, 'srfl_2d', p2=srfl_2d, reprid=GP_2D_HORIZONTAL_ID)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'srfl_2d', 'long_name', c = 'surface radiative flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'srfl_2d', 'units', c = 'W m-2')
    CALL channel_halt(substr, status)
        
    CALL new_channel_object(status, modstr, 'tslnew', p2=tslnew, reprid=GP_2D_HORIZONTAL_ID)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tslnew', 'long_name', c = 'land surface temperature for sensible heat flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tslnew', 'units', c = '')
    CALL channel_halt(substr, status)
    !
    ! 3D fields
    !
    CALL new_channel_object(status, modstr, 'tm1', p3=tm1, reprid=GP_3D_MID_ID)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tm1' &
         , 'long_name', c = 'temperature')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tm1' &
         , 'units', c = 'K')
    CALL channel_halt(substr, status)
    
    CALL new_channel_object(status, modstr, 'tvirt', p3=tvirt, reprid=GP_3D_MID_ID)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tvirt' &
         , 'long_name', c = 'virtual temperature')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tvirt' &
         , 'units', c = 'K')
    CALL channel_halt(substr, status)
    
    CALL new_channel_object(status, modstr, 'press', p3=press_3d, reprid=GP_3D_MID_ID)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'press' &
         , 'long_name', c = 'pressure')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'press' &
         , 'units', c = 'pa')
    CALL channel_halt(substr, status)
    
    CALL new_channel_object(status, modstr, 'pressi', p3=pressi_3d, reprid=GP_3D_INT_ID)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'pressi' &
         , 'long_name', c = 'pressure')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'pressi' &
         , 'units', c = 'pa')
    CALL channel_halt(substr, status)
    
    ! geopotential height at interface level
    CALL new_channel_object(status, modstr,  'geopoti', &
         p3=geopoti_3d, reprid=GP_3D_INT_ID)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'geopoti', &
         'long_name', c='interface geopotential')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'geopoti', 'units', c='m2 s-2')
    CALL channel_halt(substr, status)
    !
    ! constant start temperature
    !
    tm1        = 280._dp
    radflxw_2d = 0._dp
    srfl_2d    = 0._dp
    tslnew     = 0._dp
    
    CALL end_message(modstr,'CREATE STANDARD CHANNELS/OBJECTS',substr)

  END SUBROUTINE main_data_init_memory
  
  !----------------------------------------------------------------------------
  
  SUBROUTINE main_data_init_coupling
  
    ! ------------------------------------------------------------------
    ! This soubroutine is used to set pointers
    ! (channel objects and/or tracers) for coupling to the
    ! basemodel and to other submodes.
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_blather_bi,       ONLY: start_message_bi, end_message_bi
    ! MESSy SUBMODEL CORE (SMCL)
    USE messy_main_channel,          ONLY: get_channel_object
    
    ! LOCAL
    CHARACTER(LEN=*),  PARAMETER             :: substr = 'main_data_init_coupling'
    INTEGER                                  :: status
    
    CALL start_message_bi(modstr,'COUPLING',substr)  ! log-output
    
    CALL get_channel_object(status, 'import_grid', 'INP_tm1_st', p3=tm1_imp)
    CALL channel_halt(substr//': object INP_tm1_st in channel import_grid not found!',&
         status)
         
    CALL get_channel_object(status, 'import_grid', 'INP_qm1_q', p3=qm1)
    CALL channel_halt(substr//': object INP_qm1_q in channel import_grid not found!',&
         status)
    
    CALL get_channel_object(status, 'import_grid', 'INP_tsl_tsl', p2=tslm1)
    CALL channel_halt(substr//': object INP_tsl_tsl in channel import_grid not found!',&
         status)
         
    CALL get_channel_object(status, 'import_grid', 'INP_tsw_tsw', p2=tsw)
    CALL channel_halt(substr//': object INP_tsw_tsw in channel import_grid not found!',&
         status)
         
    CALL get_channel_object(status, 'import_grid', 'INP_tsi_tsi', p2=tsi)
    CALL channel_halt(substr//': object INP_tsi_tsi in channel import_grid not found!',&
         status)
    
    CALL get_channel_object(status, 'import_grid', 'INP_alb_alb', p2=alb)
    CALL channel_halt(substr//': object INP_alb_alb in channel import_grid not found!',&
         status)

    CALL get_channel_object(status, 'import_grid', 'INP_cvs_cvs', p2=cvs)
    CALL channel_halt(substr//': object INP_cvs_cvs in channel import_grid not found!',&
         status)
         
    CALL get_channel_object(status, 'import_grid', 'INP_cvsc_cvsc', p2=cvsc)
    CALL channel_halt(substr//': object INP_cvsc_cvsc in channel import_grid not found!',&
         status)

    CALL get_channel_object(status, 'import_grid', 'INP_sni_sni', p2=sni)
    CALL channel_halt(substr//': object INP_sni_sni in channel import_grid not found!',&
         status)

    CALL get_channel_object(status, 'import_grid', 'INP_glac_glac', p2=glac)
    CALL channel_halt(substr//': object INP_glac_glac in channel import_grid not found!',&
         status)

    CALL get_channel_object(status, 'import_grid', 'INP_slm_slm', p2=slm)
    CALL channel_halt(substr//': object INP_slm_slm in channel import_grid not found!',&
         status)

    CALL get_channel_object(status, 'import_grid', 'INP_forest_forest', p2=forest)
    CALL channel_halt(substr//': object INP_forest_forest in channel import_grid not found!',&
         status)

    CALL get_channel_object(status, 'import_grid', 'INP_icecov_icecov', p2=icecov)
    CALL channel_halt(substr//': object INP_icecov_icecov in channel import_grid not found!',&
         status)
         
    CALL get_channel_object(status, 'import_grid', 'INP_seacov_seacov', p2=seacov)
    CALL channel_halt(substr//': object INP_seacov_seacov in channel import_grid not found!',&
         status)

    CALL get_channel_object(status, 'import_grid', 'INP_seaice_seaice', p2=seaice)
    CALL channel_halt(substr//': object INP_seaice_seaice in channel import_grid not found!',&
         status)

    CALL get_channel_object(status, 'import_grid', 'INP_vlt_vlt', p2=vlt)
    CALL channel_halt(substr//': object INP_vlt_vlt in channel import_grid not found!',&
         status)

    CALL get_channel_object(status, 'import_grid', 'INP_aps_aps', p2=aps)
    CALL channel_halt(substr//': object INP_aps_aps in channel import_grid not found!',&
         status)

    CALL get_channel_object(status, 'import_grid', 'INP_geosp_geosp', p2=geosp)
    CALL channel_halt(substr//': object INP_geosp_geosp in channel import_grid not found!',&
         status)

    CALL end_message_bi(modstr,'COUPLING',substr)  ! log-output

  END SUBROUTINE main_data_init_coupling
  
  !----------------------------------------------------------------------------
  
  SUBROUTINE main_data_local_start
  
    USE messy_main_grid_def_mem_bi, ONLY: jrow, ngpblks, npromz, nproma, kproma

    USE messy_rad_fubrad_mem,  ONLY: nswlev
    USE messy_main_timer,      ONLY: lstart

    
    IMPLICIT NONE
    
    CHARACTER(LEN=*),  PARAMETER             :: substr = 'main_data_local_start'
    INTEGER :: jl
    
    ! temperature changes only at FUBRad levels
    !
    IF (ladd_tte) THEN
       IF (lstart) THEN
          tm1(_RI_XYZ__(:,jrow,:)) = tm1_imp(_RI_XYZ__(:,jrow,:))
       ELSE
          !
          ! temperature changes only at FUBRad levels
          !
          tm1(_RI_XYZ__(:,jrow,1:nswlev)) = tm1(_RI_XYZ__(:,jrow,1:nswlev)) + &
                                         tte_3d(_RI_XYZ__(:,jrow,1:nswlev)) * delta_time
       END IF
    ELSE
       tm1(_RI_XYZ__(:,jrow,:)) = tm1_imp(_RI_XYZ__(:,jrow,:))
    END IF
    apm1  (:,:)    = press_3d (_RI_XYZ__(:,jrow,:))
    app1  (:,:)    = press_3d (_RI_XYZ__(:,jrow,:))
    aphm1 (:,:)    = pressi_3d(_RI_XYZ__(:,jrow,:))
    tslnew(:,jrow) = tslm1(:,jrow)

    IF ( jrow == ngpblks ) THEN
       kproma = npromz
    ELSE
       kproma = nproma
    END IF

    DO jl=1,kproma
       loland_2d(jl,jrow) = slm(jl,jrow) > 0._dp
       loglac_2d(jl,jrow) = loland_2d(jl,jrow) .AND. glac(jl,jrow) > 0._dp
    END DO
    
  END SUBROUTINE main_data_local_start
 
  !----------------------------------------------------------------------------
  
  SUBROUTINE main_data_global_start

    USE messy_main_grid_def_mem_bi, ONLY: nlev, nlevp1, ngpblks, npromz, nproma
    USE messy_main_grid_def_bi,     ONLY: gboxarea_2d, grvol, grmass, ph, pf &
                                        , deltaz
    ! BML/MESSy
    USE messy_main_timer,         ONLY: time_step_len, lstart
    USE messy_main_constants_mem, ONLY: g, M_air, R_gas, rd

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_data_global_start'
    INTEGER :: i, zjrow, zkproma, jk, ikp
    !
    IF (lv_echam) THEN
       ! PRESSURE AT LAYER INTERFACES
       DO i=1, nlevp1
          pressi_3d(_RI_XYZ__(:,:,i)) = ph(i) * aps(:,:) ! [Pa]
       END DO
       ! PRESSURE AT LAYER MID
       DO i=1, nlev
          press_3d(_RI_XYZ__(:,:,i)) = pf(i) * aps(:,:) ! [Pa]
       END DO
    ELSE
      !
      DO i = 1, nlevp1
         pressi_3d(_RI_XYZ__(:,:,i)) = ph(i)
      END DO
      !
      DO i = 1, nlev
         press_3d(_RI_XYZ__(:,:,i))  = pf(i)
      END DO
    END IF
    !
    ! AIR MASS IN GRID BOX
    grmass(_RI_XYZ__(:,:,:)) = ((pressi_3d(_RI_XYZ__(:,:,2:nlevp1)) - &
                                 pressi_3d(_RI_XYZ__(:,:,1:nlev))) / g) * &
                               SPREAD(gboxarea_2d,_IZ_XYZ__,nlev)

    IF (lstart) THEN
       tm1(:,:,:) = tm1_imp(:,:,:)
    END IF

    DO zjrow=1, nlat
       IF ( zjrow == ngpblks ) THEN
          zkproma = npromz
       ELSE
          zkproma = nproma
       END IF
       grvol(_RI_XYZ__(1:zkproma,zjrow,:)) = grmass(_RI_XYZ__(1:zkproma,zjrow,:)) / &
            ( press_3d(_RI_XYZ__(1:zkproma,zjrow,:)) * (1.0E-03_dp * M_air) &
            / (tm1(_RI_XYZ__(1:zkproma,zjrow,:)) * R_gas) )
       !
       ! calculate the virtual temperature
       !
       tvirt(_RI_XYZ__(1:zkproma,zjrow,:)) = tm1(_RI_XYZ__(1:zkproma,zjrow,:)) *  &
           (1._dp + 0.607717_dp * qm1(_RI_XYZ__(1:zkproma,zjrow,:)))
       !
       !  calculate the geopotential at interface levels
       !  --> Integrate hydrostatic equation
       !
       geopoti_3d(_RI_XYZ__(1:zkproma,zjrow,nlevp1)) = geosp(1:zkproma,zjrow)
       DO jk = nlevp1-1, 2, -1
          ikp = jk + 1
          geopoti_3d(_RI_XYZ__(1:zkproma,zjrow,jk)) =  &
                 rd*LOG(pressi_3d(_RI_XYZ__(1:zkproma,zjrow,ikp))/ &
                        pressi_3d(_RI_XYZ__(1:zkproma,zjrow,jk))) * &
                            tvirt(_RI_XYZ__(1:zkproma,zjrow,jk)) + &
                       geopoti_3d(_RI_XYZ__(1:zkproma,zjrow,ikp))
       END DO
       geopoti_3d(_RI_XYZ__(1:zkproma,zjrow,1)) = &
              rd*LOG(pressi_3d(_RI_XYZ__(1:zkproma,zjrow,3))/ &
                     pressi_3d(_RI_XYZ__(1:zkproma,zjrow,2))) * &
                         tvirt(_RI_XYZ__(1:zkproma,zjrow,1)) + &
                    geopoti_3d(_RI_XYZ__(1:zkproma,zjrow,2))

       DO jk = 1, nlev
          deltaz(_RI_XYZ__(1:zkproma,zjrow,jk)) =  &
               (geopoti_3d(_RI_XYZ__(1:zkproma,zjrow,jk)) - &
               geopoti_3d(_RI_XYZ__(1:zkproma,zjrow,jk+1))) / g
       END DO
    END DO
    
  END SUBROUTINE main_data_global_start
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE main_data_free_memory

    IMPLICIT NONE
    
    DEALLOCATE(loland_2d, loglac_2d, aphm1, apm1, app1)
    DEALLOCATE(tte_3d)

  END SUBROUTINE main_data_free_memory
  !----------------------------------------------------------------------------

!*****************************************************************************
END MODULE messy_main_data_bi
!*****************************************************************************
