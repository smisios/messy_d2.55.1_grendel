!========================================================================
!SUBMODEL interface layer (SMIL) riutines for messy submodel ACCF
!Author: Feijia Yin 
! 
!Reference: see messy_accf.f90
!
!========================================================================
!
!****************************************************
MODULE messy_accf_e5
!****************************************************
   USE messy_main_blather_bi,     ONLY: start_message_bi, end_message_bi,&
                                        error_bi, info_bi, warning_bi
   USE messy_main_constants_mem,  ONLY: DP
  
   USE messy_accf

   IMPLICIT NONE

   ! Pointer for output
   REAL(DP), DIMENSION(:,:,:), POINTER  :: gp_atr20_o3 => NULL()
   REAL(DP), DIMENSION(:,:,:), POINTER  :: gp_atr20_ch4 => NULL()
   REAL(DP), DIMENSION(:,:,:), POINTER  :: gp_atr20_h2o => NULL()
   REAL(DP), DIMENSION(:,:,:), POINTER  :: gp_atr20_contrail=> NULL()
   REAL(DP), DIMENSION(:,:,:), POINTER  :: gp_atr20_co2=> NULL()
   
   !pointer for coupled channel objects
   REAL(DP), DIMENSION(:,:,:), POINTER  :: OLR
   REAL(DP), DIMENSION(:,:,:), POINTER  :: potcov
   REAL(DP), DIMENSION(:,:,:), POINTER  :: PV
   REAL(DP), DIMENSION(:,:),   POINTER  :: cossza
   REAL(DP), POINTER                    :: dec_sun
   !
   CHARACTER(len=20),DIMENSION(2)   :: c_LR_rad
   CHARACTER(len=20),DIMENSION(2)   :: c_pv_tropop
   CHARACTER(len=20),DIMENSION(2)   :: c_cossza_orbit
   CHARACTER(len=20),DIMENSION(2)   :: c_dec_orbit
   CHARACTER(len=20),DIMENSION(2)   :: c_potcov_contrail

 !Public subroutine
   PUBLIC :: accf_initialize
   PUBLIC :: accf_init_memory
   PUBLIC :: accf_init_coupling
   PUBLIC :: accf_physc
  ! private sub routine 
   PRIVATE :: accf_read_nml_cpl

CONTAINS
  !############################################################
  !PUBLIC SUBROUTINE
  !###########################################################
!====================================================================
 SUBROUTINE accf_initialize

  !###########################################################
  !This subroutine is used to 
  !- read (and broadcast) the CTRL-namelist.
  !- read (and broadcast) the CPL-namelist.
  !- perform the basic setup of the submodel.
  !##########################################################
   USE messy_main_mpi_bi,  ONLY : p_parallel_io, p_io, p_bcast
   USE messy_main_tools,   ONLY : find_next_free_unit

   IMPLICIT NONE

   CHARACTER(LEN=*),   PARAMETER  :: substr = 'accf_initialize'
   INTEGER                        :: status
   INTEGER                        :: iou

   CALL start_message_bi(modstr, 'INITIALIZATION', substr) 

   IF (p_parallel_io) THEN
     iou = find_next_free_unit(100,200)
     CALL accf_read_nml_ctrl(status, iou)
     IF (status/=0) CALL error_bi('Error in reading CTRL namelist',substr)
   END IF 
  ! Broadcast the CTRL namelist entries from I/O-PE to all other PEs 
   CALL p_bcast(l_cpl_rad,     p_io)
   CALL p_bcast(l_cpl_tropop,  p_io)
   CALL p_bcast(l_cpl_orbit,   p_io)
   CALL p_bcast(l_cpl_contrail,p_io)
   
   IF (p_parallel_io) THEN
     iou = find_next_free_unit(100,200)
     CALL accf_read_nml_cpl(status, iou)
     IF (status/=0) CALL error_bi('Error in reading CPL namelist',substr)
   END IF 
   
  ! Broadcast the CPL namelist entries from I/O-PE to all other PEs 
   CALL p_bcast(c_LR_rad,         p_io)
   CALL p_bcast(c_pv_tropop,      p_io)
   CALL p_bcast(c_cossza_orbit,   p_io)
   CALL p_bcast(c_potcov_contrail,p_io)
   CALL p_bcast(c_dec_orbit,      p_io)


   CALL end_message_bi(modstr, 'INITIALIZATION', substr)

 END SUBROUTINE accf_initialize
!==========================================================================
 
!=========================================================================  
 SUBROUTINE accf_init_memory

   !############################################################
   ! This subroutine is used to request memory for the submodel.
   ! The preferable method is to use "channel objects"
   ! Allocate your own memory, only if absolutely required
   !#############################################################

   USE messy_main_channel_error_bi, ONLY : channel_halt
   USE messy_main_channel_bi,       ONLY : GP_3D_MID
   USE messy_main_channel,          ONLY : new_channel, new_channel_object, &
                                           new_attribute
   
   IMPLICIT NONE
   
   CHARACTER(LEN=*),            PARAMETER   :: substr = 'accf_init_memory'
   INTEGER                                  :: status

   CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)
     CALL new_channel(status, modstr//'_gp', reprid=GP_3D_MID)
     CALL channel_halt(substr,status)
     
     CALL new_channel_object(status,modstr//'_gp','atr20_o3', p3=gp_atr20_o3)
     CALL channel_halt(substr,status)
     CALL new_attribute(status, modstr//'_gp','atr20_o3','long_name',&
                        c='ATR20 O3 by NO2')
     CALL channel_halt(substr,status)
     CALL new_attribute(status,modstr//'_gp', 'atr20_o3', 'units', &
                        c = 'K/kg(NO2)')
     CALL channel_halt(substr,status)
     CALL info_bi('channel/object'//modstr//'_gp/'//'atr20_o3'//&
                 'was created', substr)

     CALL new_channel_object(status,modstr//'_gp','atr20_ch4', p3=gp_atr20_ch4)
     CALL channel_halt(substr,status)
     CALL new_attribute(status, modstr//'_gp','atr20_ch4','long_name',&
                        c='ATR20 CH4 by NO2')
     CALL channel_halt(substr,status)
     CALL new_attribute(status,modstr//'_gp', 'atr20_ch4', 'units', &
                        c = 'K/kg(NO2)')
     CALL channel_halt(substr,status)
     CALL info_bi('channel/object'//modstr//'_gp/'//'atr20_ch4'//&
                 'was created', substr)

     CALL new_channel_object(status,modstr//'_gp','atr20_h2o', p3=gp_atr20_h2o)
     CALL channel_halt(substr,status)
     CALL new_attribute(status, modstr//'_gp','atr20_h2o','long_name',&
                        c='ATR20 H2O by NO2')
     CALL channel_halt(substr,status)
     CALL new_attribute(status,modstr//'_gp', 'atr20_h2o', 'units', &
                        c = 'K/kg(fuel)')
     CALL channel_halt(substr,status)
     CALL info_bi('channel/object'//modstr//'_gp/'//'atr20_h2o'//&
                 'was created', substr)

     CALL new_channel_object(status,modstr//'_gp','atr20_contrail', p3=gp_atr20_contrail)
     CALL channel_halt(substr,status)
     CALL new_attribute(status, modstr//'_gp','atr20_contrail','long_name',&
                        c='ATR20 contrail')
     CALL channel_halt(substr,status)
     CALL new_attribute(status,modstr//'_gp', 'atr20_contrail', 'units', &
                        c = 'K/km')
     CALL channel_halt(substr,status)
     CALL info_bi('channel/object'//modstr//'_gp/'//'atr20_contrail'//&
                 'was created', substr)
     
     CALL new_channel_object(status,modstr//'_gp','atr20_co2', p3=gp_atr20_co2)
     CALL channel_halt(substr,status)
     CALL new_attribute(status, modstr//'_gp','atr20_co2','long_name',&
                        c='ATR20 co2')
     CALL channel_halt(substr,status)
     CALL new_attribute(status,modstr//'_gp', 'atr20_co2', 'units', &
                        c = 'K/kg(fuel)')
     CALL channel_halt(substr,status)
     CALL info_bi('channel/object'//modstr//'_gp/'//'atr20_co2'//&
                 'was created', substr)

   CALL end_message_bi(modstr,'CHANNEL DEFINITION',substr)
 END SUBROUTINE accf_init_memory
!=======================================================================

!=======================================================================
SUBROUTINE accf_init_coupling

   USE messy_main_channel,          ONLY : get_channel_object
   USE messy_main_channel_error_bi, ONLY : channel_halt
   
   IMPLICIT NONE

   CHARACTER(LEN=*), PARAMETER  :: substr='accf_init_coupling'
   INTEGER                      :: status

   CALL start_message_bi(modstr,'COUPLING INITIALIZATION',substr)   
   IF(l_cpl_rad)THEN 
      CALL get_channel_object(status, TRIM(c_LR_rad(1)),TRIM(c_LR_rad(2)),p3=OLR) 
      CALL channel_halt(substr,status) 
   ENDIF
   IF (l_cpl_tropop) THEN
      CALL get_channel_object(status, TRIM(c_pv_tropop(1)),TRIM(c_pv_tropop(2)),p3=PV) 
      CALL channel_halt(substr,status) 
   ENDIF
   IF (l_cpl_orbit) THEN
      CALL get_channel_object(status,TRIM(c_cossza_orbit(1)),TRIM(c_cossza_orbit(2)),p2=cossza) 
      CALL channel_halt(substr,status)
      CALL get_channel_object(status,TRIM(c_dec_orbit(1)),TRIM(c_dec_orbit(2)),p0=dec_sun) 
      CALL channel_halt(substr,status) 
   ENDIF
   IF (l_cpl_contrail)THEN
      CALL get_channel_object(status, TRIM(c_potcov_contrail(1)),TRIM(c_potcov_contrail(2)),p3=potcov) 
      CALL channel_halt(substr,status)
   ENDIF
   CALL end_message_bi(modstr,'COUPLING INITIALIZATION',substr)
   
END SUBROUTINE accf_init_coupling
!======================================================================

!=====================================================================
SUBROUTINE accf_physc
   
   USE messy_main_data_bi,        ONLY : tm1_3d,geopot_3d,geosp,&
                                         rhum_3d
   USE messy_main_grid_def_bi,     ONLY: philat_2d,philon_2d
   USE messy_main_grid_def_mem_bi, ONLY: kproma, nlev, jrow
! op_pj_20180622: non-MESSy-conform USEage of ECHAM5 routine ...
!!$   USE mo_time_control,     ONLY:  current_date,get_year_day
   USE messy_main_timer,    ONLY : current_date,get_year_day=>YearDay

   IMPLICIT NONE
   
   INTEGER   :: DAYOFYEAR

   DAYOFYEAR = NINT(get_year_day(current_date))
    
   CALL accf_atr20_cal(nlev,kproma,jrow,&
                        geopot_3d,geosp,tm1_3d,rhum_3d,PV,OLR,cossza,&
                        dec_sun,philat_2d,philon_2d,DAYOFYEAR,potcov,&
                        gp_atr20_o3,gp_atr20_ch4,gp_atr20_h2o,&
                        gp_atr20_contrail,gp_atr20_co2)
                                   

END SUBROUTINE accf_physc
!=====================================================================
!=====================================================================
SUBROUTINE accf_read_nml_cpl(status,iou)
  ! ------------------------------------------------------------------
    ! This subroutine is used to read the CPL-namelist of the submodel.
    ! This subroutine is executed when starting very first calculations,
    ! also restarting.
    ! ------------------------------------------------------------------

    ! MESSy
    USE messy_main_tools,  ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /CPL/ c_LR_rad,c_pv_tropop,c_cossza_orbit,c_dec_orbit,c_potcov_contrail



    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='accf_read_nml_cpl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status

    ! INITIALIZE
    status = 1 !ERROR


    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist
   
    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

END SUBROUTINE accf_read_nml_cpl
!=====================================================================
!******************************************************
END MODULE messy_accf_e5
!******************************************************
