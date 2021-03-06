! **********************************************************************
! 
! SUBMODEL CORE LAYER (SMCL) ROUTINES FOR MESSy SUBMODEL airtraf
!
! THIS SUBMODEL IS USED TO CALCULATE AIR TRAFFIC EMISSIONS
!
! Author : Hiroshi Yamashita, DLR-IPA, December 2011
!          Volker Grewe, DLR-IPA, December 2011
!
! References:
!
! * None yet
!
! see SMIL for more details
! **********************************************************************

! **********************************************************************
MODULE messy_airtraf
! **********************************************************************

  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP, OneDay, g, atm2Pa   !20140724-2
!!$  USE messy_main_timer,         ONLY: Julian_date_start
!We should use ONLY sentence below to select using variables. 
  USE messy_airtraf_gc 
! USE messy_airtraf_wind        !20131214-1, 20140417-5
  USE messy_airtraf_tools_ga,                 ONLY: armoga  !20140417-5
!  USE messy_main_tools,         ONLY: nn_index  !see20130909-6,20130812-5
    
  IMPLICIT NONE
  PRIVATE

! PUBLIC :: DP  20140619-1

  ! ----------- <

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'airtraf'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '1.1'

  ! Parameters for emission calculation
  REAL(DP), PARAMETER, PUBLIC :: H2O_INDEX               = 1230.0_dp  ! [g(H2O)/kg(fuel)], see20121204-4, 
! CTRL-NAMELIST PARAMETERS

  INTEGER, PUBLIC :: nwaypoints, ngroutes         
  INTEGER, PUBLIC :: nwaypoints_out, ngroutes_out   
  INTEGER, PUBLIC :: option_traj_calc             ! 0:great circle  1: time optimal 2: fuel optimal, etc.
  INTEGER, PUBLIC :: option_output                ! 0:standard, 1: ac locations in addition
  INTEGER, PUBLIC :: option_max_ac_reached        ! 0:Stop 1:ignore flight 
  LOGICAL, PUBLIC :: lupdate_traj                  
  LOGICAL, PUBLIC :: ldaily_fp                    ! repeat flight next day

  ! Parameter settings for properties
! INTEGER, PARAMETER, PUBLIC                  :: nprops=9   !5
  !Yin_20170423
  INTEGER, PARAMETER, PUBLIC                  :: nprops=16  !5
  !Yin_20170423

  ! on individual PEs
  INTEGER, PUBLIC            :: nlroutes, nlroutes_out ! currently not used 
  INTEGER, PUBLIC            :: nl_flightplan, ng_flightplan ! currently not used

  ! GLOBAL PARAMETERS
  INTEGER, PARAMETER                          :: len_city_code = 3     ! City Code: MUC=Munich, JFK=NY
  INTEGER, PARAMETER, PUBLIC                  :: NOx  = 1
  INTEGER, PARAMETER, PUBLIC                  :: H2O  = 2
  INTEGER, PARAMETER, PUBLIC                  :: DIST = 3
  INTEGER, PARAMETER, PUBLIC                  :: FUEL = 4
  ! These four parameters below should be moved in another file, in order to 'USE' by other programs, 20170216-5
  !Yin_20170423
  INTEGER, PARAMETER, PUBLIC                  :: C_PC  = 5
  INTEGER, PARAMETER, PUBLIC                  :: ATR20O3 = 6
  INTEGER, PARAMETER, PUBLIC                  :: ATR20CH4 = 7
  INTEGER, PARAMETER, PUBLIC                  :: ATR20H2O = 8
  INTEGER, PARAMETER, PUBLIC                  :: ATR20CPC = 9
  INTEGER, PARAMETER, PUBLIC                  :: ATR20CO2 = 10
  INTEGER, PARAMETER, PUBLIC                  :: ATR20TOT = 11
  !Yin_20170423

  INTEGER, PARAMETER, PUBLIC                  :: GC = 0
  INTEGER, PARAMETER, PUBLIC                  :: Wind_opt = 1
  INTEGER, PARAMETER, PUBLIC                  :: Fuel_opt = 2          ! 20160301-1 
  INTEGER, PARAMETER, PUBLIC                  :: NOx_opt  = 3          ! 20170215-2 
  INTEGER, PARAMETER, PUBLIC                  :: H2O_opt  = 4          ! 20170222-1
  !Yin_20170423
  INTEGER, PARAMETER, PUBLIC                  :: ContrailPC_opt = 5    ! 20170430-1
  !Yin_20170423
  INTEGER, PARAMETER, PUBLIC                  :: Cost_opt = 6          ! 20170224-1
  INTEGER, PARAMETER, PUBLIC                  :: COC_opt  = 7          ! 20170316-1
  INTEGER, PARAMETER, PUBLIC                  :: ATR20_opt = 8           ! 20170801 
  INTEGER, PARAMETER, PUBLIC                  :: COSTCLIM_opt = 9           ! 20170801 
  INTEGER, PARAMETER, PUBLIC                  :: COSTCPC_opt = 10           ! 20170801 

  ! These parameters should be utilized from SMIL !!!
  INTEGER, PARAMETER, PRIVATE  :: AC_POS=10 , AC_POS_OLD=11 
  INTEGER, PARAMETER, PRIVATE  :: add_lon  = 1
  INTEGER, PARAMETER, PRIVATE  :: add_lat  = 2
  INTEGER, PARAMETER, PRIVATE  :: add_alt  = 3 
  INTEGER, PARAMETER, PRIVATE  :: add_emis = 4

! DOUBLE PRECISION, PRIVATE    :: total_flight_time_s2 ![s],see20131003-3,20131114-3
  REAL(DP)                     :: gc_opt_fl_time         !20140527-3
  REAL(DP)                     :: gc_opt_bada_fuel_nom   !20140601-3
  !Yin_20170423
  REAL(DP)                     :: gc_opt_cpc             ![km]20161212
  REAL(DP)                     :: contrail_opt_fl_time   !20140524-4
  REAL(DP)                     :: contrail_opt_cpc       ![km]20161207FY
  !Yin_20170423
  REAL(DP)                     :: wind_opt_fl_time  !, fuel_opt_fl_time   !20140524-4,20160302-1,20160606-7
  REAL(DP), DIMENSION(:), ALLOCATABLE  :: gc_opt_vtas_along_traj    !(nwaypoints) !20140527-3 
  REAL(DP), DIMENSION(:), ALLOCATABLE  :: wind_opt_vtas_along_traj !,fuel_opt_vtas_along_traj !(nwaypoints)![m/s] 20140524-6,20160302-1
!                                                                                             !20160606-7
  !Yin_20170423
  REAL(DP), DIMENSION(:), ALLOCATABLE  :: contrail_opt_vtas_along_traj !(nwaypoints)  ![m/s] 201612076
  !Yin_20170423

  REAL(DP), DIMENSION(:,:), ALLOCATABLE  :: wind_opt_vtas !,fuel_opt_vtas  !(nwaypoints,101)  ![m/s] 20150313-2,20160302-1,20160606-7
  !Yin_20170423
  REAL(DP), DIMENSION(:,:), ALLOCATABLE  :: contrail_opt_vtas          !(nwaypoints,101)  ![m/s] 20161207
  !Yin_20170423

  REAL(DP), DIMENSION(:,:), ALLOCATABLE  :: wind_opt_vground !,fuel_opt_vground !(nwaypoints,101)  ![m/s] 20150313-2,20160302-1
  REAL(DP), DIMENSION(:,:), ALLOCATABLE  :: wind_opt_potcov  !(nwaypoints,101)  ![-] 20170608-2
  !Yin_20170423
  REAL(DP), DIMENSION(:,:), ALLOCATABLE  :: contrail_opt_vground       !(nwaypoints,101)  ![m/s] 20161207
  !Yin_20170423

  REAL(DP), DIMENSION(:,:), ALLOCATABLE  :: wind_opt_rho !,fuel_opt_vground !(nwaypoints,101)  ![kg/m^3] 20160623-2
!                                                                                            
  INTEGER                      :: fl_direction, fl_distance   !20141215-1,20141216-1
! Module data 


  ! ----------- >
  !DUMMY VARIABLES
   
!  REAL(DP), DIMENSION(:,:,:,1), PUBLIC :: d_ac_routes
!  REAL(DP), DIMENSION(:,:), PUBLIC     :: d_p2_ac_routes_desc

!  REAL(DP), DIMENSION(:,:,:), PUBLIC   :: d_glob_NOx_emis 
!  REAL(DP), PUBLIC                     :: d_NOx

  ! ----------- <

! Part 1 - Logical units
!
! Flightplan:
! ===========
!

  ! PUBLIC SUBROUTINES (to be called from messy_airtraf_e5.f90)
  PUBLIC :: airtraf_read_nml_ctrl  
! PUBLIC :: prepare_trajectory_and_flightplan     ! Reorganise Flightplan etc. in case of departure
  PUBLIC :: calculate_trajectory
  PUBLIC :: calculate_emissions_along_trajectory
! PUBLIC :: update_trajectory
  PUBLIC :: fly_aircraft
! PUBLIC :: arrival
  PUBLIC :: compress_traj
! PUBLIC :: append_traj
  PUBLIC :: fill_emission_field
  ! PRIVATE SUBROUTINES
  ! ### add your own private subroutines here
!  PRIVATE ::  

CONTAINS
  ! =========================================================================
!  SUBROUTINE prepare_trajectory_and_flightplan
  
!  IMPLICIT NONE

!  TYPE(t_flightplan), pointer            :: a_traj, i_traj
  ! I/O
! TYPE(t_flightplan), INTENT(INOUT)      :: flightplan,flightplan_last
! TYPE(t_trajectory), INTENT(INOUT)      :: activ_traj,inactiv_traj

  ! Prepare new Trajectory
  !    -go to last active 
  !    -check if inactive_traj available
  !     yes:  add first inactive to last activ
  !     no: if option_max_ac_reached=yes: do return else stop
  !    -move first inactive to last active
  !    -put flightplan data into new trajectory
  !    -if daily flightplan then move to the end and adapt time (+24h)
  !    -delete flightplan item
  
  

!  END SUBROUTINE prepare_trajectory_and_flightplan
  ! =========================================================================

  ! =========================================================================
  SUBROUTINE calculate_trajectory(d_ac_routes, d_p2_ac_routes_desc,    &
                                  d_philon, d_philat, d_zgl_geopot_3d, &
                                  d_uwind_g, d_vwind_g, d_v_z_g, d_t_scb_g, d_fl_direction, d_cpc_g, &
                                  d_rho_air_dry_3d_g, d_press_3d_g,d_ATR20O3_g,&
                                  d_ATR20CH4_g,d_ATR20H2O_g,d_ATR20CPC_g,d_ATR20CO2_g) 
    ! ------------------------------------------------------------------
    ! This routine calculates waypoints regarding options. This routine  
    ! determins flightplans.
    !
    ! Outputs are waypoints info (lon,lat,alt,time,ac_speed,dist) added into ac_routes properties.
    ! If an option is selected, the option is applied for all flights. 
    !
    ! Option 0   : Great circle at fixed altitude.
    !
    ! Option 1   : Time-optimal trajectories are calculated by Optimization modules. 
    !              
    ! Option 2   : Fuel_optinal. d_rho_air_dry_3d_g(OPTIONAL) is used for this option.20160301-1
    ! Option 3   : NOx_optinal.  d_press_3d_g(OPTIONAL) is used for this option.20170215-2
    ! Option 4   : H2O_optinal.  H2O_INDEX(const) is multiplied to Fuel use. 20170222-1
    !                            The code is almost similar to the Fuel_opt. 
    ! Option 5   : Contrail_optinal. Yin_20140423 contrail potential coverage in added. 
    !              Contrail potential coverage optimal trajectories are clculated (this might be changing later on 20161207) 
    ! Option 6   : Cost_optimal. Time- and fuel-related costs are included in obj function. 20170223-5
    !                            The code is almost similar to the Fuel_opt. 
    ! Option 7   : Cash operating cost optimal.
    !
    ! Option 8   : ATR20_optional,Climate cost function will be used. 
    ! Option 9   : COSTCLIM_optional,Climate cost function will be used. 
    ! Option 1-9 : Treajectory Optimization will be performed within this subroutine. 
    !  
    !
    ! NOTE: This routine is executed only once, when an aircraft reaches its departure time.
    !
    ! ------------------------------------------------------------------
    ! TIPS:
    !     1. d_p2_ac_routes_desc(D_time) is defined in Julian_date, provided by SMIL. 
    !
    !     2. If we use Option-0(GC), we can set 'd_p2_ac_routes_desc' as INTENT(IN)
    !
    !     3. See 20170531-2

     IMPLICIT NONE 
!    REAL(DP), DIMENSION(:,:), INTENT(OUT)                :: d_ac_routes
     REAL(DP), DIMENSION(:,:), INTENT(INOUT)              :: d_ac_routes
!    INTEGER, DIMENSION(:), INTENT(INOUT)                 :: d_p2_ac_routes_desc
!GC  REAL(DP), DIMENSION(:), INTENT(IN)                   :: d_p2_ac_routes_desc
     REAL(DP), DIMENSION(:), INTENT(INOUT)                :: d_p2_ac_routes_desc
     REAL(DP), DIMENSION(:), INTENT(INOUT)                :: d_philon          !global field,20140212-1,20140527-3
     REAL(DP), DIMENSION(:), INTENT(INOUT)                :: d_philat          !global field
     REAL(DP), DIMENSION(:,:,:), INTENT(INOUT)            :: d_zgl_geopot_3d   !global field 
     REAL(DP), DIMENSION(:,:,:), INTENT(INOUT)            :: d_uwind_g         !global field 
     REAL(DP), DIMENSION(:,:,:), INTENT(INOUT)            :: d_vwind_g         !global field
     REAL(DP), DIMENSION(:,:,:), INTENT(INOUT)            :: d_v_z_g           !global field 
     REAL(DP), DIMENSION(:,:,:), INTENT(INOUT)            :: d_t_scb_g         !global field  20140522-2 
     ! Yin_20140423
     REAL(DP), DIMENSION(:,:,:), INTENT(INOUT)            :: d_cpc_g           !global field  20161207
     REAL(DP), DIMENSION(:,:,:), INTENT(IN), OPTIONAL     :: d_ATR20O3_g       !global field  20170801
     REAL(DP), DIMENSION(:,:,:), INTENT(IN), OPTIONAL     :: d_ATR20CH4_g      !global field  20170801
     REAL(DP), DIMENSION(:,:,:), INTENT(IN), OPTIONAL     :: d_ATR20H2O_g           !global field  20170801
     REAL(DP), DIMENSION(:,:,:), INTENT(IN), OPTIONAL     :: d_ATR20CPC_g           !global field  20170801
     REAL(DP), DIMENSION(:,:,:), INTENT(IN), OPTIONAL     :: d_ATR20CO2_g           !global field  20170801
     ! Yin_20140423
     INTEGER, INTENT(OUT)                                 :: d_fl_direction    !20141217-1                                            
     REAL(DP), DIMENSION(:,:,:), INTENT(IN), OPTIONAL     :: d_rho_air_dry_3d_g !global field  20160301-2,20170203-1 
     REAL(DP), DIMENSION(:,:,:), INTENT(IN), OPTIONAL     :: d_press_3d_g       !global field  20170215-2 
     REAL(DP), DIMENSION(:,:,:), ALLOCATABLE              :: temp_rho_air_dry_3d_g !global field  20170203-1 
     REAL(DP), DIMENSION(:,:,:), ALLOCATABLE              :: temp_press_3d_g    !global field  20170216-4
     REAL(DP), DIMENSION(:,:,:), ALLOCATABLE              :: temp_ATR20O3_3d_g    !global field  20170801
     REAL(DP), DIMENSION(:,:,:), ALLOCATABLE              :: temp_ATR20CH4_3d_g    !global field  20170801
     REAL(DP), DIMENSION(:,:,:), ALLOCATABLE              :: temp_ATR20H2O_3d_g    !global field  20170801
     REAL(DP), DIMENSION(:,:,:), ALLOCATABLE              :: temp_ATR20CPC_3d_g    !global field  20170801
     REAL(DP), DIMENSION(:,:,:), ALLOCATABLE              :: temp_ATR20CO2_3d_g    !global field  20170801
     INTEGER                                              :: wind_opt_gen  !,fuel_opt_gen  !HY20140420-3,20160302-1,20160606-7
     !Yin_20170423
     INTEGER                                              :: contrail_opt_gen       !FY 20161207
     !Yin_20170423
!    INTEGDER                                             :: i_lat, j_lon, i_alt   !20130812-5,20130905-3,20130909-6 
     INTEGER                                              :: nlon_rho, nlev_rho, ngl_rho   !20170203-1     

     !Check: global field,20130916-4
!    write(*,*)'d_philon',shape(d_philon),d_philon(:)
!    write(*,*)'d_philat',shape(d_philat),d_philat(:)
!    write(*,*)'d_zgl_geopot_3d',shape(d_zgl_geopot_3d),d_zgl_geopot_3d(:,:,:)
!    write(*,*)'d_uwind_g',shape(d_uwind_g),d_uwind_g(:,:,:)
!    write(*,*)'d_vwind_g',shape(d_vwind_g),d_vwind_g(:,:,:)
!    write(*,*)'d_v_z_g',shape(d_v_z_g),d_v_z_g(:,:,:)
!    write(*,*)'airtraff90_d_cpc_g:',shape(d_cpc_g),d_cpc_g(:,:,:)
!    write(*,*)'airtraff90_d_t_scb_g:',shape(d_t_scb_g),d_t_scb_g(:,:,:)
!    write(*,*)'d_press_3d_g:',shape(d_press_3d_g),d_press_3d_g(:,:,:)
     
     !Allocation of vtas arrays, see20140721-1,2 
     ALLOCATE(gc_opt_vtas_along_traj(nwaypoints))
     ALLOCATE(wind_opt_vtas_along_traj(nwaypoints))
     ALLOCATE(wind_opt_vtas(nwaypoints,101))      !20150313-2
     ALLOCATE(wind_opt_vground(nwaypoints,101))   !20150313-2
     ALLOCATE(wind_opt_rho(nwaypoints,101))   !20160623-2,20170215-1
     ALLOCATE(wind_opt_potcov(nwaypoints,101))   !20170608-2


     !Yin_20170423
     ALLOCATE(contrail_opt_vtas_along_traj(nwaypoints))
     ALLOCATE(contrail_opt_vtas(nwaypoints,101))      !20150313-2
     ALLOCATE(contrail_opt_vground(nwaypoints,101))   !20150313-2
     !Yin_20170423

     IF(PRESENT(d_rho_air_dry_3d_g))THEN          !implemented for fuel_opt, NOx_opt, H2O_opt, Cost_opt.
        nlon_rho=SIZE(d_rho_air_dry_3d_g,dim=1)   !lon,20170203-1          
        nlev_rho=SIZE(d_rho_air_dry_3d_g,dim=2)   !lev,20170203-1
        ngl_rho =SIZE(d_rho_air_dry_3d_g,dim=3)   !lat,20170203-1 
        ALLOCATE(temp_rho_air_dry_3d_g(nlon_rho, nlev_rho, ngl_rho))     
        temp_rho_air_dry_3d_g=d_rho_air_dry_3d_g  !20170203-1
!        write(*,*)"nlon_rho,nlev_rho,ngl_rho",nlon_rho,nlev_rho,ngl_rho
     ENDIF
     IF(PRESENT(d_press_3d_g))THEN          !implemented for NOx_opt. 20170216-4
!       nlon_rho=SIZE(d_press_3d_g,dim=1)   !lon,20170203-1          
!       nlev_rho=SIZE(d_press_3d_g,dim=2)   !lev,20170203-1
!       ngl_rho =SIZE(d_press_3d_g,dim=3)   !lat,20170203-1 
        ALLOCATE(temp_press_3d_g(nlon_rho, nlev_rho, ngl_rho))     
        temp_press_3d_g=d_press_3d_g        !20170216-4
!       write(*,*)"nlon_rho,nlev_rho,ngl_rho",nlon_rho,nlev_rho,ngl_rho
     ENDIF
     IF(PRESENT(d_ATR20O3_g))THEN          !implemented for ATR20_opt,COSTCLIM_opt. 20170801
        ALLOCATE(temp_ATR20O3_3d_g(nlon_rho, nlev_rho, ngl_rho))     
        temp_ATR20O3_3d_g=d_ATR20O3_g        !20170801
     ENDIF
     IF(PRESENT(d_ATR20CH4_g))THEN          !implemented for ATR20_opt,COSTCLIM_opt. 20170801
        ALLOCATE(temp_ATR20CH4_3d_g(nlon_rho, nlev_rho, ngl_rho))     
        temp_ATR20CH4_3d_g=d_ATR20CH4_g        !20170801
     ENDIF
     IF(PRESENT(d_ATR20H2O_g))THEN          !implemented for ATR20_opt,COSTCLIM_opt. 20170801
        ALLOCATE(temp_ATR20H2O_3d_g(nlon_rho, nlev_rho, ngl_rho))     
        temp_ATR20H2O_3d_g=d_ATR20H2O_g        !20170801
     ENDIF
     IF(PRESENT(d_ATR20CPC_g))THEN          !implemented for ATR20_opt,COSTCLIM_opt. 20170801
        ALLOCATE(temp_ATR20CPC_3d_g(nlon_rho, nlev_rho, ngl_rho))     
        temp_ATR20CPC_3d_g=d_ATR20CPC_g        !20170801
     ENDIF
     IF(PRESENT(d_ATR20CO2_g))THEN          !implemented for ATR20_opt,COSTCLIM_opt. 20170801
        ALLOCATE(temp_ATR20CO2_3d_g(nlon_rho, nlev_rho, ngl_rho))     
        temp_ATR20CO2_3d_g=d_ATR20CO2_g        !20170801
     ENDIF


     SELECT CASE(option_traj_calc)
!------------------------------------------------------------------
!    Routing option: GC=0
!------------------------------------------------------------------
     CASE(GC)  
!        write(*,*)'SMCL:option GC'
        CALL gc_calculate(nwaypoints, d_ac_routes(:,:),d_p2_ac_routes_desc(:),d_philon(:),d_philat(:),          & 
                          d_zgl_geopot_3d(:,:,:),d_uwind_g(:,:,:),d_vwind_g(:,:,:),d_v_z_g(:,:,:),              &
                          d_t_scb_g(:,:,:), d_cpc_g(:,:,:), gc_opt_fl_time, gc_opt_bada_fuel_nom, gc_opt_vtas_along_traj,  &
                          fl_direction, fl_distance, gc_opt_cpc) !HY20140527-3,20141215-1,20141216-1
        contrail_opt_cpc = gc_opt_cpc  ![km]
!       write(*,*)'props_emis_h2o=',props_emis_h2o
!       write(*,*)'total_flight_time_s',total_flight_time_s
!       total_flight_time_s2=total_flight_time_s    ![s],20131003-3, This value from sub-sub SMCL.20131114-3 

!------------------------------------------------------------------
!    Routing option: Wind_opt=1
!------------------------------------------------------------------
     CASE(Wind_opt) 
     !  CALL wind_opt_calculate(d_ac_routes(:,:), d_p2_ac_routes_desc(:))
     !20130904-2 
     !IF(PRESENT(d_philon).and.PRESENT(d_philat).and.PRESENT(d_zgl_geopot_3d).and. &
     !   PRESENT(d_uwind_g).and.PRESENT(d_vwind_g).and.PRESENT(d_v_z_g))THEN
!        write(*,*)'SMCL:option Wind_opt',Wind_opt
!20140417-2        CALL wind_opt_calculate(d_ac_routes(:,:),d_p2_ac_routes_desc(:),d_philon(:),d_philat(:),  &
!                                          d_zgl_geopot_3d(:,:,:),d_uwind_g(:,:,:),d_vwind_g(:,:,:),d_v_z_g(:,:,:)) !20140212-1 
!20140417-5
        CALL armoga(option_traj_calc, nwaypoints, d_ac_routes(:,:),d_p2_ac_routes_desc(:),d_philon(:),d_philat(:), &
                    d_zgl_geopot_3d(:,:,:),d_uwind_g(:,:,:),d_vwind_g(:,:,:),d_v_z_g(:,:,:),              &
                    d_t_scb_g(:,:,:),d_cpc_g(:,:,:), wind_opt_fl_time, wind_opt_gen, wind_opt_vtas_along_traj,           &
                    fl_direction, wind_opt_vtas, wind_opt_vground, wind_opt_rho, wind_opt_potcov, & !20140420-3,20141216-3,20150313-1 
                    contrail_opt_cpc)                                                               !20170608-2 
     ! The nn_index can be used in SMCL for index search. 20130812-5
     ! In SMCL, nearest points in each lon, alt, lat direction are searched by nn_index 
     ! within optimization routine: internal index search (global) with nn_index.
     !    - refer wind values at each waypoint. 
     !
     ! Then i_lat and j_lon should be determined.
     ! d_philon and d_philat are geographical coordinates of global field.
     !    d_philon(64(=nlon)): 0..->>..+360[deg]
     !    d_philat(32(=ngl)) : +90..->>..-90[deg]
     !  
     ! Calculate global horizontal nearest indices (i_lat, j_lon).20130813-2
     !! CALL nn_index(d_philat(:), lat value on waypoints(my value), i_lat)
     !! CALL nn_index(d_philon(:), lon value on waypoints(my value), j_lon)
     !
     ! Calculate global vertical index (i_alt)
     !coding on zalt using geopotential altitude.20130905-3
     !! CALL nn_index( , my value, i_alt)
     !
     ! Use global wind field:20130813-2,20130905-3
     !!     d_uwind_g(j_lon, i_alt, i_lat)
     !!     d_vwind_g(j_lon, i_alt, i_lat)
     !!       d_v_z_g(j_lon, i_alt, i_lat)
  
!------------------------------------------------------------------
!    Routing option: Fuel_opt=2, H2O_opt=4, Cost_opt=6, COC_opt=7
!------------------------------------------------------------------
     CASE(Fuel_opt, H2O_opt, Cost_opt, COC_opt,COSTCPC_opt)           !20170222-1,20170224-1,20170316-1  
!       IF(PRESENT(d_rho_air_dry_3d_g))THEN               !20170203-1
!          temp_rho_air_dry_3d_g=d_rho_air_dry_3d_g       !copy global field to use it in 'CALL armoga'  
!           IF(option_traj_calc==Fuel_opt)write(*,*)'SMCL:option Fuel_opt',Fuel_opt
!           IF(option_traj_calc==H2O_opt)write(*,*)'SMCL:option H2O_opt',H2O_opt      !20170222-1
!           IF(option_traj_calc==Cost_opt)write(*,*)'SMCL:option Cost_opt',Cost_opt   !20170224-1
!           IF(option_traj_calc==COC_opt)write(*,*)'SMCL:option COC_opt',COC_opt      !20170316-1
!           IF(option_traj_calc==COSTCPC_opt)write(*,*)'SMCL:option COSTCPC_opt',COSTCPC_opt      !20170316-1
           CALL armoga(option_traj_calc, nwaypoints, d_ac_routes(:,:),d_p2_ac_routes_desc(:),d_philon(:),d_philat(:), &
                       d_zgl_geopot_3d(:,:,:),d_uwind_g(:,:,:),d_vwind_g(:,:,:),d_v_z_g(:,:,:),              &
                       d_t_scb_g(:,:,:),d_cpc_g(:,:,:), wind_opt_fl_time, wind_opt_gen, wind_opt_vtas_along_traj,           &
                       fl_direction, wind_opt_vtas, wind_opt_vground, wind_opt_rho, wind_opt_potcov, & !20140420-3,20141216-3,20150313-1  
                       contrail_opt_cpc, temp_rho_air_dry_3d_g(:,:,:)) !,    & !20160301-2,20160302-1,20160606-6,7,20170203-1
                      !temp_rho_air_dry_3d_g(:,:,:),                    & !20160301-2,20160302-1,20160606-6,7,20170203-1
                      !f_opt_fuel_use)                                    !20160623-2,20170213-1,20170222-1,20170224-1
           !write(*,*)"f_opt_fuel_use_check:",f_opt_fuel_use,"[kg]"       !20170213-1,20170608-2
           !After optmization, final waypoints info are added into ac_routes(nwaypoints,1:6) and are backed to SMIL.
!       ENDIF
           
     !ELSE
     !   write(*,*)'ERROR:subroutine calc_traj_wind_opt'
     !ENDIF

!------------------------------------------------------------------
!    Routing option: NOx_opt=3
!------------------------------------------------------------------
     CASE(NOx_opt)  
!       IF(PRESENT(d_rho_air_dry_3d_g))THEN               !20170203-1,20170215-2
!          temp_rho_air_dry_3d_g=d_rho_air_dry_3d_g       !copy global field to use it in 'CALL armoga'  
!          temp_press_3d_g=
!           write(*,*)'SMCL:option NOx_opt',NOx_opt
           CALL armoga(option_traj_calc, nwaypoints, d_ac_routes(:,:),d_p2_ac_routes_desc(:),d_philon(:),d_philat(:), &
                       d_zgl_geopot_3d(:,:,:),d_uwind_g(:,:,:),d_vwind_g(:,:,:),d_v_z_g(:,:,:),              &
                       d_t_scb_g(:,:,:), d_cpc_g(:,:,:), wind_opt_fl_time, wind_opt_gen, wind_opt_vtas_along_traj,           &
                       fl_direction, wind_opt_vtas, wind_opt_vground, wind_opt_rho, wind_opt_potcov, & !20140420-3,20141216-3,20150313-1  
                       contrail_opt_cpc, temp_rho_air_dry_3d_g(:,:,:), & !20160301-2,20160302-1,20160606-6,7,20170203-1
                       temp_press_3d_g(:,:,:))    !,                 & !20170215-3,20170608-2
                      !f_opt_fuel_use)                                    !20160623-2,20170213-1
           !write(*,*)"f_opt_fuel_use_check:",f_opt_fuel_use,"[kg]"       !20170213-1
           !After optmization, final waypoints info are added into ac_routes(nwaypoints,1:6) and are backed to SMIL.
!       ENDIF
           
     !ELSE
     !   write(*,*)'ERROR:subroutine calc_traj_wind_opt'
     !ENDIF

!Yin_20170423
!------------------------------------------------------------------
!   Routing option: ContrailPC_opt=5
!------------------------------------------------------------------
     !If combine ContrailPC_opt into Wind_opt, see HY20170608-4 
     CASE(ContrailPC_opt)
!       write(*,*)'SMCL:option ContrailPC_opt'
        CALL armoga(option_traj_calc, nwaypoints, d_ac_routes(:,:),d_p2_ac_routes_desc(:),d_philon(:),d_philat(:),     &
                    d_zgl_geopot_3d(:,:,:),d_uwind_g(:,:,:),d_vwind_g(:,:,:),d_v_z_g(:,:,:),        &
                    d_t_scb_g(:,:,:), d_cpc_g(:,:,:), contrail_opt_fl_time, contrail_opt_gen,       &
                    contrail_opt_vtas_along_traj,fl_direction, contrail_opt_vtas,                   &
                    contrail_opt_vground, wind_opt_rho, wind_opt_potcov, contrail_opt_cpc)      !20170608-2  
        !option_traj_calc & wind_opt_rho are added by HY just in case. wind_opt_rho should be changed as contrail_opt_rho
!------------------------------------------------------------------
!  Routing option: ATR20_opt = 8, COSTCLIM=9
!------------------------------------------------------------------
     CASE(ATR20_opt,COSTCLIM_opt)
!        IF(option_traj_calc==ATR20_opt)write(*,*)'SMCL:option ATR20_opt',ATR20_opt
!        IF(option_traj_calc==COSTCLIM_opt)write(*,*)'SMCL:option COSTCLIM_opt',COSTCLIM_opt
        CALL armoga(option_traj_calc, nwaypoints,d_ac_routes(:,:),d_p2_ac_routes_desc(:),d_philon(:),d_philat(:),  &
                    d_zgl_geopot_3d(:,:,:),d_uwind_g(:,:,:),d_vwind_g(:,:,:),d_v_z_g(:,:,:),  &
                    d_t_scb_g(:,:,:),d_cpc_g(:,:,:), wind_opt_fl_time, wind_opt_gen,        &
                    wind_opt_vtas_along_traj,fl_direction, wind_opt_vtas,                   &
                    wind_opt_vground,wind_opt_rho,wind_opt_potcov,contrail_opt_cpc,          &
                    temp_rho_air_dry_3d_g(:,:,:),temp_press_3d_g(:,:,:),                    &              
                    temp_ATR20O3_3d_g(:,:,:),temp_ATR20CH4_3d_g(:,:,:),&
                    temp_ATR20H2O_3d_g(:,:,:),temp_ATR20CPC_3d_g(:,:,:),temp_ATR20CO2_3d_g(:,:,:))
!Yin_20170423
!------------------------------------------------------------------
!  DUMMY 
!------------------------------------------------------------------
     CASE DEFAULT 
        d_ac_routes(:,:)       = 1.0_dp
        d_p2_ac_routes_desc(:) = 1.0_dp
!        write(*,*)'dummy_calculate_traj'
     ENDSELECT
!     write(*,*)'CPCpropcheck_GA_2:',d_ac_routes(:,props_potcov)
!     write(*,*)'CPCpropcheck_GA_3:',shape(d_cpc_g),d_cpc_g(:,:,:)
!     write(*,*)'CPCpropcheck_tscb_GA_3:',shape(d_t_scb_g),d_t_scb_g(:,:,:)
     !add fl_direction into d_fl_direction,20141217-1
     d_fl_direction = fl_direction
     !USAGE CONFIRMATION
!     write(*,*)"calculate_trajectory subroutine is done"
     !For fuel_opt, NOx_opt, H2O_opt, Cost_opt. from line 256 
     IF(ALLOCATED(temp_rho_air_dry_3d_g))DEALLOCATE(temp_rho_air_dry_3d_g)   !20170203-1,20170212-3
     !For NOx_opt. from line 264 
     IF(ALLOCATED(temp_press_3d_g))DEALLOCATE(temp_press_3d_g)               !20170216-4
     IF(ALLOCATED(temp_ATR20O3_3d_g))DEALLOCATE(temp_ATR20O3_3d_g)               !20170801
     IF(ALLOCATED(temp_ATR20CH4_3d_g))DEALLOCATE(temp_ATR20CH4_3d_g)               !20170801
     IF(ALLOCATED(temp_ATR20H2O_3d_g))DEALLOCATE(temp_ATR20H2O_3d_g)               !20170801
     IF(ALLOCATED(temp_ATR20CPC_3d_g))DEALLOCATE(temp_ATR20CPC_3d_g)               !20170801
     IF(ALLOCATED(temp_ATR20CO2_3d_g))DEALLOCATE(temp_ATR20CO2_3d_g)               !20170801

  END SUBROUTINE calculate_trajectory
  ! =========================================================================

  ! =========================================================================
!HY20120404-1
  SUBROUTINE calculate_emissions_along_trajectory(d_ac_routes, d_pass_phyval_along_traj, d_i_traj_dep, d_p_pe,&
                                                  d_philon,d_philat,d_zgl_geopot_3d,&         ! for index calculation
                                                  d_ATR20O3_g,d_ATR20CH4_g,d_ATR20H2O_g,d_ATR20CPC_g,&
                                                  d_ATR20CO2_g)

    ! ------------------------------------------------------------------
    ! This routine calculates expected emissions along a trajectory.
    ! Emission indeces are used now. In furute, emission will be calculated 
    ! by considering meteorological conditions at departure time, with a 
    ! emission caluculation module.
    !
    ! The emissions are calculated by using: 
    !    ac_routes properties (lon,lat,alt,time,ac_speed,dist).
    ! The emission values are added as  
    ! one of the properties in ac_routes.
    !
    ! Regarding ac_routes(dist) property:see opt20130131-5B
    !    (GC)The 'dist' corresponds to the true_distance between waypoints.
    !    (Wind_opt)The 'dist' means (true_distance - distance flying by winds) 
    !              that is, 'dist' <= true_distance
    !              -> IS THIS IDEA STILL CORRECT, when applying total energy model ??
    !
    ! Another simulation code (or climate cost function) will be adapted for the 
    ! emission calculations.
    ! See note 20120404-1,20120719-4,20120814-3,20120820-3
    !
    ! NOTE: This routine is executed only once, when departure time is assigned to an aircraft.
    !
    ! Note: - d_pass_phyval_along_traj(nwaypoints-1,1): pressure, [Pa]
    !       - d_pass_phyval_along_traj(nwaypoints-1,2): temperature, [K]
    !       - d_pass_phyval_along_traj(nwaypoints-1,3): rho, [kg/m^3]
    !
    ! ------------------------------------------------------------------
    ! TIPS:
    !     1. d_p2_ac_routes_desc(D_time) is defined in Julian_date, provided by SMIL. 
    !
     USE messy_airtraf_wind, ONLY:coc_calc, cost_calc
     USE messy_main_tools,   ONLY: nn_index
     USE messy_main_constants_mem, ONLY: g
     IMPLICIT NONE 
     REAL(DP), DIMENSION(:,:), INTENT(INOUT)   :: d_ac_routes
     REAL(DP), DIMENSION(:,:), INTENT(IN)      :: d_pass_phyval_along_traj 
     REAL(DP), DIMENSION(:),   INTENT(IN)      :: d_philon 
     REAL(DP), DIMENSION(:),   INTENT(IN)      :: d_philat 
     REAL(DP), DIMENSION(:,:,:),   INTENT(IN)      :: d_zgl_geopot_3d 
     INTEGER, INTENT(IN)                       :: d_i_traj_dep, d_p_pe  !20131115-1,20131125-1 
!    INTEGER, DIMENSION(:), INTENT(INOUT)      :: d_p2_ac_routes_desc
!    REAL(DP), DIMENSION(:), INTENT(INOUT)     :: d_p2_ac_routes_desc
     REAL(DP)                                  :: arr_ac_mass,total_nox_emission,total_h2o_emission,total_fl_distance  !20140926-3,4 
     REAL(DP)                                  :: total_atr20_tot  !20180525-2 
     REAL(DP)                                  :: ave_ac_mass      !20141216-4
     INTEGER                                   :: i,j,k
     REAL(DP), DIMENSION(:,:,:), INTENT(IN)    :: d_ATR20O3_g       !global field  20170801
     REAL(DP), DIMENSION(:,:,:), INTENT(IN)    :: d_ATR20CH4_g      !global field  20170801
     REAL(DP), DIMENSION(:,:,:), INTENT(IN)    :: d_ATR20H2O_g           !global field  20170801
     REAL(DP), DIMENSION(:,:,:), INTENT(IN)    :: d_ATR20CPC_g           !global field  20170801
     REAL(DP), DIMENSION(:,:,:), INTENT(IN)    :: d_ATR20CO2_g           !global field  20170801
     REAL(DP), DIMENSION(:),   ALLOCATABLE     :: d_zalt            ! geopotentia 
     !Declaration for Total Energy Model
     REAL(DP), PARAMETER                       :: extra_fuel_rate = 0.03_dp             !Percentage of extra fuel,20140724-2
     REAL(DP)                                  :: ita,Cfcr,CD2,S_area,ac_mass,new_ac_mass,CD0,Cf1,Cf2
                                                                                        !,g see20130717-3  !,rho see20130902-1 
     REAL(DP)                                  :: ac_oew, max_payload,load_factor       !20131114-3 
     REAL(DP)                                  :: aa,bb,cc,sq,de,x1,x2                  !,CL,CD,DRAG_FORCE,THRUST_FORCE 
     REAL(DP)                                  :: wind_opt_bada_fuel_nom                !Fuel from BADA A333_ptf,[kg/s] 20140524-7
     !Yin_20170423
     REAL(DP)                                  :: contrail_opt_bada_fuel_nom            !Fuel from BADA A333_ptf,[kg/s] 20140524-7
     !Yin_20170423

     REAL(DP)                                  :: wind_opt_ave_alt                      !Average altitude of wind opt traj. 20140524-7 
     !Yin_20170423
     REAL(DP)                                  :: contrail_opt_ave_alt                  !Average altitude of CPC opt traj.  
     !Yin_20170423

     REAL(DP)                                  :: temp_wind_opt_ave_alt                 !20140524-7
     !Yin_20170423
     REAL(DP)                                  :: temp_contrail_opt_ave_alt             !20140524-7
     !Yin_20170423

     !REAL(DP), DIMENSION(:,:)                 :: BADA_DATA(8,5) !this should be given by Aircraft table.20130803-1,20131108-1,20170430-1
     !Yin_20170423
     REAL(DP)                  :: BADA_DATA(8,6) !this should be given by Aircraft table.20170531-3
     !Yin_20170423
     REAL(DP)                    :: VTAS(nwaypoints)    !20140602-2
     INTEGER                                   :: j1,nlev_d_zalt,d_j_lon,d_j_lat,d_j_alt
    
     REAL(DP), DIMENSION(:,:,:)                :: traj_for_nnindex_d(nwaypoints,3,1)    !20180511-3
     !Declaration for DLR Fuel Flow Method
     REAL(DP), PARAMETER                       :: temp_sealv    = 288.15_dp             !Temperature at sea level[K],20140724-2
     INTEGER,PARAMETER                         :: num_icao_data = 4
     INTEGER,PARAMETER                         :: order_app_equ = 2
     REAL(DP),DIMENSION(:)                     :: EINOx(nwaypoints), Wfuel_ref(nwaypoints), EINOx_ref(nwaypoints),Wfuel(nwaypoints) 
     REAL(DP),DIMENSION(:)                     :: coefficient(order_app_equ+1),z_array(order_app_equ+1) 
     REAL(DP),DIMENSION(:,:)                   :: w_array(order_app_equ+1,order_app_equ+1) 
     REAL(DP),DIMENSION(:,:)                   :: ICAO_DATA(2,4) !this should be given by Aircraft table.20130803-1,20140723-1 
     REAL(DP)                                  :: pressure_total, temperature_total, num_engine, theta_total, delta_total
     REAL(DP)                                  :: hmd_corr_factor, specific_hmd !,pressure_actual,temperature_actual 
     REAL(DP)                                  :: s_chol, temp_XXN_d
     INTEGER                                   :: j2

     !Cost(cost_opt=6) and COC(coc_opt=7) calculations for all routing options
     REAL(DP)                                  :: coc_total_flight_time                 ![s],20170428-2
     REAL(DP)                                  :: coc_total_fuel_use                    ![kg(fuel)],20170428-2
     REAL(DP)                                  :: coc_value                             ![US Dollar],20170428-2
     REAL(DP)                                  :: cost_value                            ![US Dollar],20170428-6

     nlev_d_zalt = SIZE(d_zgl_geopot_3d,dim=2)   !ALT: in 2nd dimention.
     ALLOCATE(d_zalt(nlev_d_zalt))               !20140429-2

     !COMMON DATA
     !ICAO DATA input,20140723-1
     DATA ((ICAO_DATA(i,j),i=1,2),j=1,4) /0.228_dp,4.88_dp,0.724_dp,12.66_dp,2.245_dp,22.01_dp,2.767_dp,28.72_dp/     

     !BADA DATA input,A333_ptf,BADA_ISA_TABLE, 20131108-1, 20140524-7
     DATA ((BADA_DATA(i,j),i=1,8),j=1,3) /280.0_dp,290.0_dp,310.0_dp,330.0_dp,350.0_dp,370.0_dp,390.0_dp,410.0_dp,         &  !FL[ft]   
                                          305.79_dp,304.48_dp,301.86_dp,299.21_dp,296.54_dp,295.07_dp,295.07_dp,295.07_dp, &  !a[m/s]  
                                          112.4_dp,108.7_dp,101.7_dp,95.5_dp,90.0_dp,85.5_dp,81.9_dp,79.0_dp/      !Fuel[nom, kg/min]
     do i=1,8
        BADA_DATA(i,4)=BADA_DATA(i,1)*100.0_dp*0.30480_dp   !Altitude [m],20131108-1
        BADA_DATA(i,5)=BADA_DATA(i,3)/60.0_dp               !Fuel [NOM case, kg/s]
        !Yin_20170423
        BADA_DATA(i,6)=BADA_DATA(i,5)*0.75_dp               !Fuel cost [NOM case, dollar/s]
        !Yin_20170423
        !write(*,*)'BADA_DATA_check',(BADA_DATA(i,j),j=1,5)
     enddo

     !FROM BADA&ICAO date (from d_ac_routes_desc? or other new Aircraft table structure?) 
     !Specific aircraft/engine data should be included into new Aircraft table structure.
     ac_oew          = 125100.0_dp       !=Mmin from BADA A333_opf,[kg],20131114-3
     max_payload     = 47900.0_dp        !=Mpyld from BADA A333_opf,[kg],20131114-3
     load_factor     = 0.62_dp         !=0.609,20130426-1  ; change to 0.62(ICAO 2008),20131126-1B
     S_area          = 361.6_dp        !from BADA A333_opf,[m^2]

     CD0             = 0.019805_dp     ![-]
     CD2             = 0.031875_dp     ![-]
     Cf1             = 0.61503_dp      ![kg/min/kN]
     Cf2             = 919.03_dp       ![Kt]
     Cfcr            = 0.93655_dp      ![-]
     num_engine      = 2.0_dp          !A330 2_engines

     !Initialization for Total Energy Model
!    FUEL_U(1)       = 0.0_dp    !see20130711-1
!    g               = 9.81_dp   ![m/s^2],see20130717-3
!    CL              = 0.0_dp
!    CD              = 0.0_dp
     ita             = 0.0_dp 
     aa              = 0.0_dp
     bb              = 0.0_dp
     cc              = 0.0_dp
     de              = 0.0_dp
     sq              = 0.0_dp
     x1              = 0.0_dp
     x2              = 0.0_dp
     wind_opt_bada_fuel_nom     = 0.0_dp   !20140524-7
     !Yin_20170423
     contrail_opt_bada_fuel_nom = 0.0_dp   !20161215
     !Yin_20170423

     wind_opt_ave_alt           = 0.0_dp   !20140524-7
     !Yin_20170423
     contrail_opt_ave_alt       = 0.0_dp   !20161215
     !Yin_20170423

     arr_ac_mass     = 0.0_dp          !20140926-3
     ave_ac_mass     = 0.0_dp          !20141216-4
     total_nox_emission  = 0.0_dp          !20140926-3
     total_h2o_emission  = 0.0_dp          !20140926-3 
     total_fl_distance   = 0.0_dp          !20140926-4
     total_atr20_tot    = 0.0_dp          !20180525-2
     !Initialization for DLR Fuel Flow Method 
     EINOx           = 0.0_dp
     EINOx_ref       = 0.0_dp
     Wfuel_ref       = 0.0_dp
     Wfuel           = 0.0_dp
     delta_total     = 0.0_dp
     theta_total     = 0.0_dp
     pressure_total  = 0.0_dp
     temperature_total = 0.0_dp
     hmd_corr_factor = 0.0_dp
     specific_hmd    = 0.0_dp
     coefficient     = 0.0_dp
!    z_array         = 0.0_dp
!    w_array         = 0.0_dp
!    s_chol          = 0.0_dp

     !Initialization for Cost/COC calculations
     coc_total_flight_time = 0.0_dp                ![s],20170428-2
     coc_total_fuel_use    = 0.0_dp                ![kg(fuel)],20170428-2
     coc_value             = 0.0_dp                ![US Dollar],20170428-2
     cost_value            = 0.0_dp                ![US Dollar],20170428-6

     !FROM EMAC
!    rho               = 0.4_dp        ![kg/m^3] at 10,000m ,see20130902-1
!    pressure_actual   = 26437.0_dp    ![Pa] at 10,000m ,see20130902-2
!    temperature_actual= 223.15_dp     ![K]  at 10,000m ,see20130902-2 

     !FROM BADA&ICAO date (from d_ac_routes_desc? or other new Aircraft table structure?) 
     !Specific aircraft/engine data should be included into new Aircraft table structure.
     !ac_oew          = 125100_dp       !=Mmin from BADA A333_opf,[kg],20131114-3
     !max_payload     = 47900_dp        !=Mpyld from BADA A333_opf,[kg],20131114-3
     !load_factor     = 0.62_dp         !=0.609,20130426-1  ; change to 0.62(ICAO 2008),20131126-1B
     !S_area          = 361.6_dp        !from BADA A333_opf,[m^2]

     !Initial AC mass 1st rough estimation, by using BADA data. 
     SELECT CASE(option_traj_calc)     !20140524-5
     CASE(GC)                          !20140602-1
        !gc_opt_bada_fuel_nom,[kg/s] 20140601-3 
!       ac_mass         = 174000_dp ![kg]  "m_ref" value from A333_ptf, A333_opf 
!       ac_mass         = total_flight_time_s2*fuel_consumption*0.03_dp + 125100_dp + 47900_dp*0.609_dp ![kg],20130725-1,20131025-1
!       ac_mass         = gc_opt_fl_time*BADA_FUEL_NOM*0.03_dp+ac_oew+max_payload*load_factor ![kg],20131114-2,-3
        ac_mass         = gc_opt_fl_time*gc_opt_bada_fuel_nom*extra_fuel_rate  &     !20140724-2
                         + ac_oew + max_payload*load_factor                          ![kg]20140601-3,20140602-1
        VTAS            = gc_opt_vtas_along_traj          !See 20140524-9   Volker:250.0_dp,245.55966_dp [m/s],M=0.82
                                                          !20130719-1,20131114-2,20131216-3,20140602-2,20150306-2
!        write(*,*)"(GC)gc_opt_fl_time=",gc_opt_fl_time,"[s]" 
!        write(*,*)"(GC)gc_opt_vtas_along_traj=",gc_opt_vtas_along_traj,"[m/s]"
!        write(*,*)'(GC)initial_ac_mass=',ac_mass,"[kg]"
!        write(*,*)'(GC)fuel_nom=',gc_opt_bada_fuel_nom,"[kg/s]"
!        write(*,*)'(GC)VTAS=',VTAS,"[m/s]"
        !Yin_20170423
!        write(*,*)'(GC)Contrail_potcov=',gc_opt_cpc,"[km]"
        !Yin_20170423
       
        !Output to compare flight time,20140612-2,20141215-2,4,20141216-5,20150310-4
!        IF(d_p_pe.eq.0)THEN
!        IF((d_p_pe.ge.0).and.(d_p_pe.le.5))THEN
!           IF(fl_direction.eq.0)THEN       !Eastbound
!              OPEN(81,file='ft_gc_pe0_eb.dat',position='append')
!                 write(81,*)d_i_traj_dep, gc_opt_fl_time, d_ac_routes(1,props_alt), fl_direction, "Eastbound"      
!                 !  write(10,'(a)')">"
!              CLOSE(81)
!              OPEN(11,file='vtas_alt_pe0_eb.dat',position='append')   !20150310-4
!                 do i=1,nwaypoints
!                    write(11,*)d_i_traj_dep,VTAS(i),d_ac_routes(i,props_alt)
!                 enddo
!                 write(11,*)
!                 write(11,*)
!              CLOSE(11)
!              OPEN(13,file='vground_alt_pe0_eb.dat',position='append')   !20150310-4
!                 do i=1,nwaypoints-1
!                    write(13,*)d_i_traj_dep,d_ac_routes(i,props_ac_speed)/3.6_dp,d_ac_routes(i,props_alt)
!                 enddo
!                 write(13,*)
!                 write(13,*)
!              CLOSE(13)
!           ELSEIF(fl_direction.eq.2)THEN   !Westbound
!              OPEN(82,file='ft_gc_pe0_wb.dat',position='append')
!                 write(82,*)d_i_traj_dep, gc_opt_fl_time, d_ac_routes(1,props_alt), fl_direction, "Westbound"     
!                 !  write(10,'(a)')">"
!              CLOSE(82)
!              OPEN(12,file='vtas_alt_pe0_wb.dat',position='append')   !20150310-4
!                 do i=1,nwaypoints
!                    write(12,*)d_i_traj_dep,VTAS(i),d_ac_routes(i,props_alt)
!                 enddo
!                 write(12,*)
!                 write(12,*)
!              CLOSE(12)
!              OPEN(14,file='vground_alt_pe0_wb.dat',position='append')   !20150310-4
!                 do i=1,nwaypoints-1
!                    write(14,*)d_i_traj_dep,d_ac_routes(i,props_ac_speed)/3.6_dp,d_ac_routes(i,props_alt)
!                 enddo
!                 write(14,*)
!                 write(14,*)
!              CLOSE(14)
!           ENDIF
!!          OPEN(10,file='ft_gc_001.dat',position='append')
!!             write(10,*)gc_opt_fl_time
!          !!  write(10,'(a)')">"
!!          CLOSE(10)
!!          OPEN(31,file='arr_ac_mass_001.dat',position='append')
!!             write(31,*)ac_mass        !20140926-1,20141215-2
!!          CLOSE(31)
!        ENDIF
!!        IF(d_p_pe.eq.1)THEN
!       IF((d_p_pe.ge.6).and.(d_p_pe.le.11))THEN
!           IF(fl_direction.eq.0)THEN       !Eastbound
!              OPEN(83,file='ft_gc_pe1_eb.dat',position='append')
!                 write(83,*)d_i_traj_dep, gc_opt_fl_time, d_ac_routes(1,props_alt), fl_direction, "Eastbound"      
!                 !  write(10,'(a)')">"
!              CLOSE(83)
!              OPEN(15,file='vtas_alt_pe1_eb.dat',position='append')   !20150310-4
!                 do i=1,nwaypoints
!                    write(15,*)d_i_traj_dep,VTAS(i),d_ac_routes(i,props_alt)
!                 enddo
!                 write(15,*)
!                 write(15,*)
!              CLOSE(15)
!              OPEN(17,file='vground_alt_pe1_eb.dat',position='append')   !20150310-4
!                 do i=1,nwaypoints-1
!                    write(17,*)d_i_traj_dep,d_ac_routes(i,props_ac_speed)/3.6_dp,d_ac_routes(i,props_alt)
!                 enddo
!                 write(17,*)
!                 write(17,*)
!              CLOSE(17)
!           ELSEIF(fl_direction.eq.2)THEN   !Westbound
!              OPEN(84,file='ft_gc_pe1_wb.dat',position='append')
!                 write(84,*)d_i_traj_dep, gc_opt_fl_time, d_ac_routes(1,props_alt), fl_direction, "Westbound"     
!                 !  write(10,'(a)')">"
!              CLOSE(84)
!              OPEN(16,file='vtas_alt_pe1_wb.dat',position='append')   !20150310-4
!                 do i=1,nwaypoints
!                    write(16,*)d_i_traj_dep,VTAS(i),d_ac_routes(i,props_alt)
!                 enddo
!                 write(16,*)
!                 write(16,*)
!              CLOSE(16)
!              OPEN(18,file='vground_alt_pe1_wb.dat',position='append')   !20150310-4
!                 do i=1,nwaypoints-1
!                    write(18,*)d_i_traj_dep,d_ac_routes(i,props_ac_speed)/3.6_dp,d_ac_routes(i,props_alt)
!                 enddo
!                 write(18,*)
!                 write(18,*)
!              CLOSE(18)
!           ENDIF
!!          OPEN(20,file='ft_gc_002.dat',position='append')
!!             write(20,*)gc_opt_fl_time
!          !!  write(20,'(a)')">"
!!          CLOSE(20)
!!          OPEN(32,file='arr_ac_mass_002.dat',position='append')
!!             write(32,*)ac_mass        !20140926-1
!!          CLOSE(32)
!        ENDIF
!!        IF(d_p_pe.eq.2)THEN
!        IF((d_p_pe.ge.12).and.(d_p_pe.le.17))THEN
!           IF(fl_direction.eq.0)THEN       !Eastbound
!              OPEN(85,file='ft_gc_pe2_eb.dat',position='append')
!                 write(85,*)d_i_traj_dep, gc_opt_fl_time, d_ac_routes(1,props_alt), fl_direction, "Eastbound"      
!                 !  write(10,'(a)')">"
!              CLOSE(85)
!              OPEN(21,file='vtas_alt_pe2_eb.dat',position='append')   !20150310-4
!                 do i=1,nwaypoints
!                    write(21,*)d_i_traj_dep,VTAS(i),d_ac_routes(i,props_alt)
!                 enddo
!                 write(21,*)
!                 write(21,*)
!              CLOSE(21)
!              OPEN(23,file='vground_alt_pe2_eb.dat',position='append')   !20150310-4
!                 do i=1,nwaypoints-1
!                    write(23,*)d_i_traj_dep,d_ac_routes(i,props_ac_speed)/3.6_dp,d_ac_routes(i,props_alt)
!                 enddo
!                 write(23,*)
!                 write(23,*)
!              CLOSE(23)
!           ELSEIF(fl_direction.eq.2)THEN   !Westbound
!              OPEN(86,file='ft_gc_pe2_wb.dat',position='append')
!                 write(86,*)d_i_traj_dep, gc_opt_fl_time, d_ac_routes(1,props_alt), fl_direction, "Westbound"     
!                 !  write(10,'(a)')">"
!              CLOSE(86)
!              OPEN(22,file='vtas_alt_pe2_wb.dat',position='append')   !20150310-4
!                 do i=1,nwaypoints
!                    write(22,*)d_i_traj_dep,VTAS(i),d_ac_routes(i,props_alt)
!                 enddo
!                 write(22,*)
!                 write(22,*)
!              CLOSE(22)
!              OPEN(24,file='vground_alt_pe2_wb.dat',position='append')   !20150310-4
!                 do i=1,nwaypoints-1
!                    write(24,*)d_i_traj_dep,d_ac_routes(i,props_ac_speed)/3.6_dp,d_ac_routes(i,props_alt)
!                 enddo
!                 write(24,*)
!                 write(24,*)
!              CLOSE(24)
!           ENDIF
!!          OPEN(30,file='ft_gc_003.dat',position='append')
!!             write(30,*)gc_opt_fl_time
!          !!  write(30,'(a)')">"
!!          CLOSE(30)
!!          OPEN(33,file='arr_ac_mass_003.dat',position='append')
!!             write(33,*)ac_mass        !20140926-1
!!          CLOSE(33)
!        ENDIF
!!        IF(d_p_pe.eq.3)THEN
!        IF((d_p_pe.ge.18).and.(d_p_pe.le.23))THEN
!           IF(fl_direction.eq.0)THEN       !Eastbound
!              OPEN(87,file='ft_gc_pe3_eb.dat',position='append')
!                 write(87,*)d_i_traj_dep, gc_opt_fl_time, d_ac_routes(1,props_alt), fl_direction, "Eastbound"      
!                 !  write(10,'(a)')">"
!              CLOSE(87)
!              OPEN(25,file='vtas_alt_pe3_eb.dat',position='append')   !20150310-4
!                 do i=1,nwaypoints
!                    write(25,*)d_i_traj_dep,VTAS(i),d_ac_routes(i,props_alt)
!                 enddo
!                 write(25,*)
!                 write(25,*)
!              CLOSE(25)
!              OPEN(27,file='vground_alt_pe3_eb.dat',position='append')   !20150310-4
!                 do i=1,nwaypoints-1
!                    write(27,*)d_i_traj_dep,d_ac_routes(i,props_ac_speed)/3.6_dp,d_ac_routes(i,props_alt)
!                 enddo
!                 write(27,*)
!                 write(27,*)
!              CLOSE(27)
!           ELSEIF(fl_direction.eq.2)THEN   !Westbound
!              OPEN(88,file='ft_gc_pe3_wb.dat',position='append')
!                 write(88,*)d_i_traj_dep, gc_opt_fl_time, d_ac_routes(1,props_alt), fl_direction, "Westbound"     
!                 !  write(10,'(a)')">"
!              CLOSE(88)
!              OPEN(26,file='vtas_alt_pe3_wb.dat',position='append')   !20150310-4
!                 do i=1,nwaypoints
!                    write(26,*)d_i_traj_dep,VTAS(i),d_ac_routes(i,props_alt)
!                 enddo
!                 write(26,*)
!                 write(26,*)
!              CLOSE(26)
!              OPEN(28,file='vground_alt_pe3_wb.dat',position='append')   !20150310-4
!                 do i=1,nwaypoints-1
!                    write(28,*)d_i_traj_dep,d_ac_routes(i,props_ac_speed)/3.6_dp,d_ac_routes(i,props_alt)
!                 enddo
!                 write(28,*)
!                 write(28,*)
!              CLOSE(28)
!           ENDIF
!!          OPEN(40,file='ft_gc_004.dat',position='append')
!!             write(40,*)gc_opt_fl_time
!          !!  write(40,'(a)')">"
!!          CLOSE(40)
!!          OPEN(34,file='arr_ac_mass_004.dat',position='append')
!!             write(34,*)ac_mass        !20140926-1
!!          CLOSE(34)
!        ENDIF

        !Flight_check output, 20141216-1
        !This is implemented only for GC case, because fl_distance varies in wind_opt case. 
!        OPEN(39,file="flight_check.dat",position='APPEND')
!           write(39,*)d_i_traj_dep, d_p_pe, fl_direction, fl_distance,                              &
!              "D_city:",d_ac_routes(1,props_lon),d_ac_routes(1,props_lat),d_ac_routes(1,props_alt), &
!              "A_city:",d_ac_routes(nwaypoints,props_lon),d_ac_routes(nwaypoints,props_lat),d_ac_routes(nwaypoints,props_alt)
!        CLOSE(39)

     CASE(Wind_opt, Fuel_opt, NOx_opt, H2O_opt, Cost_opt,COC_opt,ATR20_opt,COSTCLIM_opt,COSTCPC_opt)
        !wind_opt_bada_fuel_nom,[kg/s] 20140524-7,8
        !Calc wind_opt_ave_alt. 
        temp_wind_opt_ave_alt = 0.0_dp 
        do i=1,nwaypoints
           temp_wind_opt_ave_alt = temp_wind_opt_ave_alt + d_ac_routes(i,props_alt)  !sum altitudes 20140524-8 
        enddo
        wind_opt_ave_alt = temp_wind_opt_ave_alt/nwaypoints  !average altitude [m]    
!        write(*,*)"wind_opt_ave_alt=",wind_opt_ave_alt,"[m]"        
        
        !Calculate average fuel consumption at averaged wind-opt-traj altitude, by interplating BADA data.
        !This calcuation is moved from messy_airtraf_wind.f90, 20140524-7,8 
        do i=1,7
           if((BADA_DATA(i,4)<=wind_opt_ave_alt).and.(wind_opt_ave_alt<BADA_DATA(i+1,4)))then
              wind_opt_bada_fuel_nom=(BADA_DATA(i,5)*(wind_opt_ave_alt-BADA_DATA(i+1,4))   & 
                                     -BADA_DATA(i+1,5)*(wind_opt_ave_alt-BADA_DATA(i,4)))/ &
                                     (BADA_DATA(i,4)-BADA_DATA(i+1,4))      ![kg/s]
!              write(*,*)"if_bada_fuel_1",i,wind_opt_ave_alt,BADA_DATA(i+1,4)
              exit
           elseif(wind_opt_ave_alt.lt.BADA_DATA(1,4))then
              wind_opt_bada_fuel_nom=BADA_DATA(1,5)
!              write(*,*)"if_bada_fuel_2",i,wind_opt_ave_alt,BADA_DATA(1,4)
              exit
           elseif(wind_opt_ave_alt.ge.BADA_DATA(8,4))then
              wind_opt_bada_fuel_nom=BADA_DATA(8,5)
!              write(*,*)"if_bada_fuel_3",i,wind_opt_ave_alt,BADA_DATA(8,4)
              exit  
           endif
        enddo
!        write(*,*)'(Wind)wind_opt_bada_fuel_nom=',wind_opt_bada_fuel_nom,"[kg/s]"
       
        !wind_opt_fl_time shows total flight time [s] of wind-opt traj.
        !wind_opt_vtas_along_traj(nwaypoints)[m/s] shows vtas of wind-opt traj. 
        !Use wind_opt_fl_time instead of total_flight_time_s of GC   !20140524-4,5 
        !ac_mass        = wind_opt_fl_time*BADA_FUEL_NOM*0.03_dp+ac_oew+max_payload*load_factor
        ac_mass         = wind_opt_fl_time*wind_opt_bada_fuel_nom*extra_fuel_rate  &   !20140724-2
                         + ac_oew + max_payload*load_factor
        VTAS            = wind_opt_vtas_along_traj                 !20140524-9,20140602-2,20150306-2  
!        write(*,*)"(Wind)wind_opt_fl_time=",wind_opt_fl_time,"[s]" 
!        write(*,*)"(Wind)wind_opt_vtas_along_traj=",wind_opt_vtas_along_traj,"[m/s]"  !20140524-6
!        write(*,*)'(Wind)itnitial_ac_mass=',ac_mass,"[kg]"
!        write(*,*)'(Wind)fuel_nom=',wind_opt_bada_fuel_nom,"[kg/s]"
!        write(*,*)'(Wind)VTAS=',VTAS,"[m/s]" 
        !Yin_20170423
!        write(*,*)'(Wind)contrail_potcov=',contrail_opt_cpc,"[km]"                     !20170531-2
        !Yin_20170423

        !Output to compare flight time,20140612-2,20141215-4,20141216-5 
!        IF(d_p_pe.eq.0)THEN
!        IF((d_p_pe.ge.0).and.(d_p_pe.le.5))THEN
!           IF(fl_direction.eq.0)THEN       !Eastbound
!              OPEN(91,file='ft_wind_pe0_eb.dat',position='append')
!                 write(91,*)d_i_traj_dep, wind_opt_fl_time, wind_opt_ave_alt, fl_direction, "Eastbound"      
!                 !  write(10,'(a)')">"
!              CLOSE(91)
!              OPEN(31,file='vtas_alt_pe0_eb.dat',position='append')   !20150310-4
!                 do i=1,nwaypoints
!                    write(31,*)d_i_traj_dep,VTAS(i),d_ac_routes(i,props_alt)
!                 enddo
!                 write(31,*)
!                 write(31,*)
!              CLOSE(31)
!              OPEN(33,file='vground_alt_pe0_eb.dat',position='append')   !20150310-4
!                 do i=1,nwaypoints-1
!                    write(33,*)d_i_traj_dep,d_ac_routes(i,props_ac_speed)/3.6_dp,d_ac_routes(i,props_alt)
!                 enddo
!                 write(33,*)
!                 write(33,*)
!              CLOSE(33)
!           ELSEIF(fl_direction.eq.2)THEN   !Westbound
!              OPEN(92,file='ft_wind_pe0_wb.dat',position='append')
!                 write(92,*)d_i_traj_dep, wind_opt_fl_time, wind_opt_ave_alt, fl_direction, "Westbound"     
!                 !  write(10,'(a)')">"
!              CLOSE(92)
!              OPEN(32,file='vtas_alt_pe0_wb.dat',position='append')   !20150310-4
!                 do i=1,nwaypoints
!                    write(32,*)d_i_traj_dep,VTAS(i),d_ac_routes(i,props_alt)
!                 enddo
!                 write(32,*)
!                 write(32,*)
!              CLOSE(32)
!              OPEN(34,file='vground_alt_pe0_wb.dat',position='append')   !20150310-4
!                 do i=1,nwaypoints-1
!                    write(34,*)d_i_traj_dep,d_ac_routes(i,props_ac_speed)/3.6_dp,d_ac_routes(i,props_alt)
!                 enddo
!                 write(34,*)
!                 write(34,*)
!              CLOSE(34)
!           ENDIF
!!          OPEN(50,file='ft_wind_001.dat',position='append')
!!             write(50,*)wind_opt_fl_time
!          !!  write(50,'(a)')">"
!!          CLOSE(50)
!!          OPEN(35,file='arr_ac_mass_001.dat',position='append')
!!             write(35,*)ac_mass        !20140926-2,20141215-3
!!          CLOSE(35)
!        ENDIF
!!        IF(d_p_pe.eq.1)THEN
!        IF((d_p_pe.ge.6).and.(d_p_pe.le.11))THEN
!           IF(fl_direction.eq.0)THEN       !Eastbound
!              OPEN(93,file='ft_wind_pe1_eb.dat',position='append')
!                 write(93,*)d_i_traj_dep, wind_opt_fl_time, wind_opt_ave_alt, fl_direction, "Eastbound"      
!                 !  write(10,'(a)')">"
!              CLOSE(93)
!              OPEN(35,file='vtas_alt_pe1_eb.dat',position='append')   !20150310-4
!                 do i=1,nwaypoints
!                    write(35,*)d_i_traj_dep,VTAS(i),d_ac_routes(i,props_alt)
!                 enddo
!                 write(35,*)
!                 write(35,*)
!              CLOSE(35)
!              OPEN(37,file='vground_alt_pe1_eb.dat',position='append')   !20150310-4
!                 do i=1,nwaypoints-1
!                    write(37,*)d_i_traj_dep,d_ac_routes(i,props_ac_speed)/3.6_dp,d_ac_routes(i,props_alt)
!                 enddo
!                 write(37,*)
!                 write(37,*)
!              CLOSE(37)
!           ELSEIF(fl_direction.eq.2)THEN   !Westbound
!              OPEN(94,file='ft_wind_pe1_wb.dat',position='append')
!                 write(94,*)d_i_traj_dep, wind_opt_fl_time, wind_opt_ave_alt, fl_direction, "Westbound"     
!                 !  write(10,'(a)')">"
!              CLOSE(94)
!              OPEN(36,file='vtas_alt_pe1_wb.dat',position='append')   !20150310-4
!                 do i=1,nwaypoints
!                    write(36,*)d_i_traj_dep,VTAS(i),d_ac_routes(i,props_alt)
!                 enddo
!                 write(36,*)
!                 write(36,*)
!              CLOSE(36)
!              OPEN(38,file='vground_alt_pe1_wb.dat',position='append')   !20150310-4
!                 do i=1,nwaypoints-1
!                    write(38,*)d_i_traj_dep,d_ac_routes(i,props_ac_speed)/3.6_dp,d_ac_routes(i,props_alt)
!                 enddo
!                 write(38,*)
!                 write(38,*)
!              CLOSE(38)
!           ENDIF
!!          OPEN(60,file='ft_wind_002.dat',position='append')
!!             write(60,*)wind_opt_fl_time
!          !!  write(60,'(a)')">"
!!          CLOSE(60)
!!          OPEN(36,file='arr_ac_mass_002.dat',position='append')
!!             write(36,*)ac_mass        !20140926-2
!!          CLOSE(36)
!        ENDIF
!!        IF(d_p_pe.eq.2)THEN
!        IF((d_p_pe.ge.12).and.(d_p_pe.le.17))THEN
!           IF(fl_direction.eq.0)THEN       !Eastbound
!              OPEN(95,file='ft_wind_pe2_eb.dat',position='append')
!                 write(95,*)d_i_traj_dep, wind_opt_fl_time, wind_opt_ave_alt, fl_direction, "Eastbound"      
!                 !  write(10,'(a)')">"
!              CLOSE(95)
!              OPEN(41,file='vtas_alt_pe2_eb.dat',position='append')   !20150310-4
!                 do i=1,nwaypoints
!                    write(41,*)d_i_traj_dep,VTAS(i),d_ac_routes(i,props_alt)
!                 enddo
!                 write(41,*)
!                 write(41,*)
!              CLOSE(41)
!              OPEN(43,file='vground_alt_pe2_eb.dat',position='append')   !20150310-4
!                 do i=1,nwaypoints-1
!                    write(43,*)d_i_traj_dep,d_ac_routes(i,props_ac_speed)/3.6_dp,d_ac_routes(i,props_alt)
!                 enddo
!                 write(43,*)
!                 write(43,*)
!              CLOSE(43)
!           ELSEIF(fl_direction.eq.2)THEN   !Westbound
!              OPEN(96,file='ft_wind_pe2_wb.dat',position='append')
!                 write(96,*)d_i_traj_dep, wind_opt_fl_time, wind_opt_ave_alt, fl_direction, "Westbound"     
!                 !  write(10,'(a)')">"
!              CLOSE(96)
!              OPEN(42,file='vtas_alt_pe2_wb.dat',position='append')   !20150310-4
!                 do i=1,nwaypoints
!                    write(42,*)d_i_traj_dep,VTAS(i),d_ac_routes(i,props_alt)
!                 enddo
!                 write(42,*)
!                 write(42,*)
!              CLOSE(42)
!              OPEN(44,file='vground_alt_pe2_wb.dat',position='append')   !20150310-4
!                 do i=1,nwaypoints-1
!                    write(44,*)d_i_traj_dep,d_ac_routes(i,props_ac_speed)/3.6_dp,d_ac_routes(i,props_alt)
!                 enddo
!                 write(44,*)
!                 write(44,*)
!              CLOSE(44)
!           ENDIF
!!          OPEN(70,file='ft_wind_003.dat',position='append')
!!             write(70,*)wind_opt_fl_time
!          !!  write(70,'(a)')">"
!!          CLOSE(70)
!!          OPEN(37,file='arr_ac_mass_003.dat',position='append')
!!             write(37,*)ac_mass        !20140926-2
!!          CLOSE(37)
!        ENDIF
!!        IF(d_p_pe.eq.3)THEN
!        IF((d_p_pe.ge.18).and.(d_p_pe.le.23))THEN
!           IF(fl_direction.eq.0)THEN       !Eastbound
!              OPEN(97,file='ft_wind_pe3_eb.dat',position='append')
!                 write(97,*)d_i_traj_dep, wind_opt_fl_time, wind_opt_ave_alt, fl_direction, "Eastbound"      
!                 !  write(10,'(a)')">"
!              CLOSE(97)
!              OPEN(45,file='vtas_alt_pe3_eb.dat',position='append')   !20150310-4
!                 do i=1,nwaypoints
!                    write(45,*)d_i_traj_dep,VTAS(i),d_ac_routes(i,props_alt)
!                 enddo
!                 write(45,*)
!                 write(45,*)
!              CLOSE(45)
!              OPEN(47,file='vground_alt_pe3_eb.dat',position='append')   !20150310-4
!                 do i=1,nwaypoints-1
!                    write(47,*)d_i_traj_dep,d_ac_routes(i,props_ac_speed)/3.6_dp,d_ac_routes(i,props_alt)
!                 enddo
!                 write(47,*)
!                 write(47,*)
!              CLOSE(47)
!           ELSEIF(fl_direction.eq.2)THEN   !Westbound
!              OPEN(98,file='ft_wind_pe3_wb.dat',position='append')
!                 write(98,*)d_i_traj_dep, wind_opt_fl_time, wind_opt_ave_alt, fl_direction, "Westbound"     
!                 !  write(10,'(a)')">"
!              CLOSE(98)
!              OPEN(46,file='vtas_alt_pe3_wb.dat',position='append')   !20150310-4
!                 do i=1,nwaypoints
!                    write(46,*)d_i_traj_dep,VTAS(i),d_ac_routes(i,props_alt)
!                 enddo
!                 write(46,*)
!                 write(46,*)
!              CLOSE(46)
!              OPEN(48,file='vground_alt_pe3_wb.dat',position='append')   !20150310-4
!                 do i=1,nwaypoints-1
!                    write(48,*)d_i_traj_dep,d_ac_routes(i,props_ac_speed)/3.6_dp,d_ac_routes(i,props_alt)
!                 enddo
!                 write(48,*)
!                 write(48,*)
!              CLOSE(48)
!           ENDIF
!!          OPEN(80,file='ft_wind_004.dat',position='append')
!!             write(80,*)wind_opt_fl_time
!          !!  write(80,'(a)')">"
!!          CLOSE(80)
!!          OPEN(38,file='arr_ac_mass_004.dat',position='append')
!!             write(38,*)ac_mass        !20140926-2
!!          CLOSE(38)
!        ENDIF
        
        !Output Vground/Vtas for 3 city-pairs, 20150313-2
        !MUC >> JFK
!        IF(d_p_pe.eq.0 .and. d_i_traj_dep.eq.20)THEN 
!           OPEN(51,file='vtas_plane_muc_jfk.dat')
!              do j=1,101  
!                 do i=1,nwaypoints    !20150313-1,2 
!                    write(51,*)(d_ac_routes(i,props_time)-d_ac_routes(1,props_time))*OneDay, 0.0_dp+150.0_dp*(j-1),  &
!                                wind_opt_vtas(i,j) 
!                 enddo
!              enddo
!           CLOSE(51)
!           OPEN(52,file='vground_plane_muc_jfk.dat')
!              do j=1,101  
!                 do i=1,nwaypoints-1  !20150313-1, the value at A_city is 0.  
!                    write(52,*)(d_ac_routes(i,props_time)-d_ac_routes(1,props_time))*OneDay, 0.0_dp+150.0_dp*(j-1),  &
!                                0.5_dp*wind_opt_rho(i,j)*wind_opt_vtas(i,j)*wind_opt_vtas(i,j),  &   !20160623-1
!                                wind_opt_vground(i,j), wind_opt_vground(i,j)/wind_opt_vtas(i,j), wind_opt_potcov(i,j)
!                 enddo
!              enddo
!           CLOSE(52)
!           OPEN(31,file='traj_on_plane_muc_jfk.dat')
!                 do i=1,nwaypoints    !20150313-2
!                    write(31,*)(d_ac_routes(i,props_time)-d_ac_routes(1,props_time))*OneDay, d_ac_routes(i,props_alt) 
!                 enddo
!           CLOSE(31)
!        ENDIF
!        !JFK >> MUC 
!        IF(d_p_pe.eq.3 .and. d_i_traj_dep.eq.1)THEN 
!           OPEN(53,file='vtas_plane_jfk_muc.dat')
!              do j=1,101  
!                 do i=1,nwaypoints    !20150313-2 
!                    write(53,*)(d_ac_routes(i,props_time)-d_ac_routes(1,props_time))*OneDay, 0.0_dp+150.0_dp*(j-1),  &
!                                wind_opt_vtas(i,j) 
!                 enddo
!              enddo
!           CLOSE(53)
!           OPEN(54,file='vground_plane_jfk_muc.dat')
!              do j=1,101  
!                 do i=1,nwaypoints-1  !20150313-1, the value at A_city is 0.  
!                    write(54,*)(d_ac_routes(i,props_time)-d_ac_routes(1,props_time))*OneDay, 0.0_dp+150.0_dp*(j-1),  &
!                                0.5_dp*wind_opt_rho(i,j)*wind_opt_vtas(i,j)*wind_opt_vtas(i,j),  &   !20160623-1
!                                wind_opt_vground(i,j), wind_opt_vground(i,j)/wind_opt_vtas(i,j), wind_opt_potcov(i,j)
!                 enddo
!              enddo
!           CLOSE(54)
!           OPEN(32,file='traj_on_plane_jfk_muc.dat')
!                 do i=1,nwaypoints    !20150313-2
!                    write(32,*)(d_ac_routes(i,props_time)-d_ac_routes(1,props_time))*OneDay, d_ac_routes(i,props_alt) 
!                 enddo
!           CLOSE(32)
!        ENDIF
!        !EHAM >> KMSP
!        IF(d_p_pe.eq.2 .and. d_i_traj_dep.eq.14)THEN 
!           OPEN(55,file='vtas_plane_eham_kmsp.dat')
!              do j=1,101  
!                 do i=1,nwaypoints    !20150313-2 
!                    write(55,*)(d_ac_routes(i,props_time)-d_ac_routes(1,props_time))*OneDay, 0.0_dp+150.0_dp*(j-1),  &
!                                wind_opt_vtas(i,j) 
!                 enddo
!              enddo
!           CLOSE(55)
!           OPEN(56,file='vground_plane_eham_kmsp.dat')
!              do j=1,101  
!                 do i=1,nwaypoints-1  !20150313-1, the value at A_city is 0.  
!                    write(56,*)(d_ac_routes(i,props_time)-d_ac_routes(1,props_time))*OneDay, 0.0_dp+150.0_dp*(j-1),  &
!                                0.5_dp*wind_opt_rho(i,j)*wind_opt_vtas(i,j)*wind_opt_vtas(i,j),  &   !20160623-1
!                                wind_opt_vground(i,j), wind_opt_vground(i,j)/wind_opt_vtas(i,j), wind_opt_potcov(i,j)
!                 enddo
!              enddo
!           CLOSE(56)
!           OPEN(33,file='traj_on_plane_eham_kmsp.dat')
!                 do i=1,nwaypoints    !20150313-2
!                    write(33,*)(d_ac_routes(i,props_time)-d_ac_routes(1,props_time))*OneDay, d_ac_routes(i,props_alt) 
!                 enddo
!           CLOSE(33)
!        ENDIF
!        !KMSP >> EHAM
!        IF(d_p_pe.eq.1 .and. d_i_traj_dep.eq.24)THEN 
!           OPEN(57,file='vtas_plane_kmsp_eham.dat')
!              do j=1,101  
!                 do i=1,nwaypoints    !20150313-2 
!                    write(57,*)(d_ac_routes(i,props_time)-d_ac_routes(1,props_time))*OneDay, 0.0_dp+150.0_dp*(j-1),  &
!                                wind_opt_vtas(i,j) 
!                 enddo
!              enddo
!           CLOSE(57)
!           OPEN(58,file='vground_plane_kmsp_eham.dat')
!              do j=1,101  
!                 do i=1,nwaypoints-1  !20150313-1, the value at A_city is 0.  
!                    write(58,*)(d_ac_routes(i,props_time)-d_ac_routes(1,props_time))*OneDay, 0.0_dp+150.0_dp*(j-1),  &
!                                0.5_dp*wind_opt_rho(i,j)*wind_opt_vtas(i,j)*wind_opt_vtas(i,j),  &   !20160623-1
!                                wind_opt_vground(i,j), wind_opt_vground(i,j)/wind_opt_vtas(i,j), wind_opt_potcov(i,j)
!                 enddo
!              enddo
!           CLOSE(58)
!           OPEN(34,file='traj_on_plane_kmsp_eham.dat')
!                 do i=1,nwaypoints    !20150313-2
!                    write(34,*)(d_ac_routes(i,props_time)-d_ac_routes(1,props_time))*OneDay, d_ac_routes(i,props_alt) 
!                 enddo
!           CLOSE(34)
!        ENDIF
!        !EHAM >> KSEA
!        IF(d_p_pe.eq.1 .and. d_i_traj_dep.eq.17)THEN 
!           OPEN(59,file='vtas_plane_eham_ksea.dat')
!              do j=1,101  
!                 do i=1,nwaypoints    !20150313-2 
!                    write(59,*)(d_ac_routes(i,props_time)-d_ac_routes(1,props_time))*OneDay, 0.0_dp+150.0_dp*(j-1),  &
!                                wind_opt_vtas(i,j) 
!                 enddo
!              enddo
!           CLOSE(59)
!           OPEN(60,file='vground_plane_eham_ksea.dat')
!              do j=1,101  
!                 do i=1,nwaypoints-1  !20150313-1, the value at A_city is 0.  
!                    write(60,*)(d_ac_routes(i,props_time)-d_ac_routes(1,props_time))*OneDay, 0.0_dp+150.0_dp*(j-1),  &
!                                0.5_dp*wind_opt_rho(i,j)*wind_opt_vtas(i,j)*wind_opt_vtas(i,j),  &   !20160623-1
!                                wind_opt_vground(i,j), wind_opt_vground(i,j)/wind_opt_vtas(i,j), wind_opt_potcov(i,j) 
!                 enddo
!              enddo
!           CLOSE(60)
!           OPEN(35,file='traj_on_plane_eham_ksea.dat')
!                 do i=1,nwaypoints    !20150313-2
!                    write(35,*)(d_ac_routes(i,props_time)-d_ac_routes(1,props_time))*OneDay, d_ac_routes(i,props_alt) 
!                 enddo
!           CLOSE(35)
!        ENDIF
!        !KSEA >> EHAM
!        IF(d_p_pe.eq.1 .and. d_i_traj_dep.eq.21)THEN 
!           OPEN(61,file='vtas_plane_ksea_eham.dat')
!              do j=1,101  
!                 do i=1,nwaypoints    !20150313-2 
!                    write(61,*)(d_ac_routes(i,props_time)-d_ac_routes(1,props_time))*OneDay, 0.0_dp+150.0_dp*(j-1),  &
!                                wind_opt_vtas(i,j) 
!                 enddo
!              enddo
!           CLOSE(61)
!           OPEN(62,file='vground_plane_ksea_eham.dat')
!              do j=1,101  
!                 do i=1,nwaypoints-1  !20150313-1, the value at A_city is 0.  
!                    write(62,*)(d_ac_routes(i,props_time)-d_ac_routes(1,props_time))*OneDay, 0.0_dp+150.0_dp*(j-1),  &
!                                0.5_dp*wind_opt_rho(i,j)*wind_opt_vtas(i,j)*wind_opt_vtas(i,j),  &   !20160623-1
!                                wind_opt_vground(i,j), wind_opt_vground(i,j)/wind_opt_vtas(i,j), wind_opt_potcov(i,j)
!                 enddo
!              enddo
!           CLOSE(62)
!           OPEN(36,file='traj_on_plane_ksea_eham.dat')
!                 do i=1,nwaypoints    !20150313-2
!                    write(36,*)(d_ac_routes(i,props_time)-d_ac_routes(1,props_time))*OneDay, d_ac_routes(i,props_alt) 
!                 enddo
!           CLOSE(36)
!        ENDIF

    !Yin_20170423
     CASE(ContrailPC_opt)
        temp_contrail_opt_ave_alt = 0.0_dp 
        do i=1,nwaypoints
           temp_contrail_opt_ave_alt = temp_contrail_opt_ave_alt + d_ac_routes(i,props_alt)  !sum altitudes 20140524-8 
        enddo
        contrail_opt_ave_alt = temp_contrail_opt_ave_alt/nwaypoints  !average altitude [m]    
!        write(*,*)"contrail_opt_ave_alt=",contrail_opt_ave_alt,"[m]"        
        
        !Calculate average fuel consumption at averaged wind-opt-traj altitude, by interplating BADA data.
        !This calcuation is moved from messy_airtraf_wind.f90, 20140524-7,8 
        do i=1,7
           if((BADA_DATA(i,4)<=contrail_opt_ave_alt).and.(contrail_opt_ave_alt<BADA_DATA(i+1,4)))then
              contrail_opt_bada_fuel_nom=(BADA_DATA(i,5)*(contrail_opt_ave_alt-BADA_DATA(i+1,4))   & 
                                     -BADA_DATA(i+1,5)*(contrail_opt_ave_alt-BADA_DATA(i,4)))/ &
                                     (BADA_DATA(i,4)-BADA_DATA(i+1,4))      ![kg/s]
!              write(*,*)"if_bada_fuel_1",i,contrail_opt_ave_alt,BADA_DATA(i+1,4)
              exit
           elseif(contrail_opt_ave_alt.lt.BADA_DATA(1,4))then
              contrail_opt_bada_fuel_nom=BADA_DATA(1,5)
!              write(*,*)"if_bada_fuel_2",i,contrail_opt_ave_alt,BADA_DATA(1,4)
              exit
           elseif(contrail_opt_ave_alt.ge.BADA_DATA(8,4))then
              contrail_opt_bada_fuel_nom=BADA_DATA(8,5)
!              write(*,*)"if_bada_fuel_3",i,contrail_opt_ave_alt,BADA_DATA(8,4)
              exit  
           endif
        enddo
!        write(*,*)'(Contrail)contrail_opt_bada_fuel_nom=',contrail_opt_bada_fuel_nom,"[kg/s]"
       
        !wind_opt_fl_time shows total flight time [s] of wind-opt traj.
        !wind_opt_vtas_along_traj(nwaypoints)[m/s] shows vtas of wind-opt traj. 
        !Use wind_opt_fl_time instead of total_flight_time_s of GC   !20140524-4,5 
        !ac_mass        = wind_opt_fl_time*BADA_FUEL_NOM*0.03_dp+ac_oew+max_payload*load_factor
        ac_mass         = contrail_opt_fl_time*contrail_opt_bada_fuel_nom*extra_fuel_rate  &   !20140724-2
                         + ac_oew + max_payload*load_factor
        VTAS            = contrail_opt_vtas_along_traj                 !20140524-9,20140602-2,20150306-2  
!        write(*,*)"(Contrail)contrail_opt_fl_time=",contrail_opt_fl_time,"[s]" 
!        write(*,*)"(Contrail)contrail_opt_vtas_along_traj=",contrail_opt_vtas_along_traj,"[m/s]"  !20140524-6
!        write(*,*)'(Contrail)itnitial_ac_mass=',ac_mass,"[kg]"
!        write(*,*)'(Contrail)fuel_nom=',contrail_opt_bada_fuel_nom,"[kg/s]"
!        write(*,*)'(Contrail)VTAS=',VTAS,"[m/s]" 
!        write(*,*)'(Contrail)contrail_potcov=',contrail_opt_cpc,"[km]" 
        
     ENDSELECT    
!--------------------------------------------------
!    Total Energy Model: calculate fuel flow 
!--------------------------------------------------
     !   - d_ac_routes(:,props_fuel_use) is calculated here.
     !
     !   - d_ac_routes(:,props_time) shows global passing time in Julian_date.  
     !     This value is affected by wind velocity in messy_airtraf_gc.f90.20140602-1
     !     NOTE: Time(Julian_date) * OneDay(86400.0) = Time(sec)
     !   
     ! AC_MACH is used from MODULE messy_airtraf_gc
     ! VTAS is passed from MODULE messy_airtraf_gc(GC option) or messy_airtraf_wind(Wind option).
     !
     ! The 'ita' is recalculated at each waypoint. 
     ! Because VTAS changes at waypoints. see 20130521-3,20130718-1,20131104-1
     ! Thus ita is included following DO structure.20140602-3   
     
     !!see 20131031-1,20130228,20130115-1,20130314-1

     !ita = Cf1*(1.0_dp + VTAS/0.514443004_dp/Cf2)/60.0_dp/1000.0_dp  ![kg/s/N] 20130521-3,20140602-3
!     open(71,file="ac_mass_pe0_eb.dat",position='APPEND')   !20131115-1,20141215-1,2
!     open(72,file="ac_mass_pe0_wb.dat",position='APPEND')   !
!     open(73,file="ac_mass_pe1_eb.dat",position='APPEND')   !fl_direction can be used. see 20141215-1 
!     open(74,file="ac_mass_pe1_wb.dat",position='APPEND')   !
!     open(75,file="ac_mass_pe2_eb.dat",position='APPEND')   !Eastbound: fl_direction = 0 
!     open(76,file="ac_mass_pe2_wb.dat",position='APPEND')   !Westbound: fl_direction = 2
!     open(77,file="ac_mass_pe3_eb.dat",position='APPEND')   
!     open(78,file="ac_mass_pe3_wb.dat",position='APPEND')   
!        if(d_p_pe.eq.0 .and. fl_direction.eq.0)write(71,*)d_i_traj_dep, ac_mass, fl_direction, "Eastbound"
!        if(d_p_pe.eq.0 .and. fl_direction.eq.2)write(72,*)d_i_traj_dep, ac_mass, fl_direction, "Westbound"
!        if(d_p_pe.eq.1 .and. fl_direction.eq.0)write(73,*)d_i_traj_dep, ac_mass, fl_direction, "Eastbound"
!        if(d_p_pe.eq.1 .and. fl_direction.eq.2)write(74,*)d_i_traj_dep, ac_mass, fl_direction, "Westbound"
!        if(d_p_pe.eq.2 .and. fl_direction.eq.0)write(75,*)d_i_traj_dep, ac_mass, fl_direction, "Eastbound"
!        if(d_p_pe.eq.2 .and. fl_direction.eq.2)write(76,*)d_i_traj_dep, ac_mass, fl_direction, "Westbound"
!        if(d_p_pe.eq.3 .and. fl_direction.eq.0)write(77,*)d_i_traj_dep, ac_mass, fl_direction, "Eastbound"
!        if(d_p_pe.eq.3 .and. fl_direction.eq.2)write(78,*)d_i_traj_dep, ac_mass, fl_direction, "Westbound"
!        if(d_p_pe.ge.0 .and. d_p_pe.le.5 .and. fl_direction.eq.0)write(71,*)d_i_traj_dep, ac_mass, fl_direction, "Eastbound"
!        if(d_p_pe.ge.0 .and. d_p_pe.le.5 .and. fl_direction.eq.2)write(72,*)d_i_traj_dep, ac_mass, fl_direction, "Westbound"
!        if(d_p_pe.ge.6 .and. d_p_pe.le.11 .and. fl_direction.eq.0)write(73,*)d_i_traj_dep, ac_mass, fl_direction, "Eastbound"
!        if(d_p_pe.ge.6 .and. d_p_pe.le.11 .and. fl_direction.eq.2)write(74,*)d_i_traj_dep, ac_mass, fl_direction, "Westbound"
!        if(d_p_pe.ge.12 .and. d_p_pe.le.17 .and. fl_direction.eq.0)write(75,*)d_i_traj_dep, ac_mass, fl_direction, "Eastbound"
!        if(d_p_pe.ge.12 .and. d_p_pe.le.17 .and. fl_direction.eq.2)write(76,*)d_i_traj_dep, ac_mass, fl_direction, "Westbound"
!        if(d_p_pe.ge.18 .and. d_p_pe.le.23 .and. fl_direction.eq.0)write(77,*)d_i_traj_dep, ac_mass, fl_direction, "Eastbound"
!        if(d_p_pe.ge.18 .and. d_p_pe.le.23 .and. fl_direction.eq.2)write(78,*)d_i_traj_dep, ac_mass, fl_direction, "Westbound"
     !20140926-3
     arr_ac_mass = ac_mass     !researve aircraft (initial) mass at A_city
     ave_ac_mass = ac_mass     !add ac_mass at A_city,20141216-4   
     !Backward calculation on ac_mass from A_city to D_city.
     DO j1=nwaypoints,2,-1
        !CL= 2.0_dp*g*ac_mass/rho/VTAS**2/S_area  
        !CD= CD0 + CD2*CL**2
        !DRAG_FORCE=0.5_dp*rho*VTAS**2*S_area*CD
        !THRUST_FORCE=DRAG_FORCE*0.001_dp
        ita = Cf1*(1.0_dp + VTAS(j1-1)/0.514443004_dp/Cf2)/60.0_dp*0.001_dp  ![kg/s/N] 20130521-3,20140602-3,20140620-8
        
        !write(*,*)'CL=',CL,'CD=',CD,'AC_V*1000=',AC_V*1000.0_dp,'VTAS=',VTAS,'DRAG_FORCE=',DRAG_FORCE,'THUST_FORCE=',THRUST_FORCE
!        write(*,*)'j1=',j1,'flight_time[s]=',(d_ac_routes(j1,props_time)-d_ac_routes(j1-1,props_time))*OneDay, &
!                  'VTAS=',VTAS(j1-1),'ita=',ita,'ac_mass=',ac_mass
     
        !20130711-1,20130522-1
       !aa = ita*Cfcr*((d_ac_routes(j1,props_time)-d_ac_routes(j1-1,props_time))*OneDay)*CD2*2.0_dp*(g**2)/VTAS**2/rho/S_area
        aa = ita*Cfcr*((d_ac_routes(j1,props_time)-d_ac_routes(j1-1,props_time))*OneDay)* & 
             CD2*2.0_dp*(g**2)/VTAS(j1-1)**2/d_pass_phyval_along_traj(j1-1,3)/S_area   !20130902-1,20140602-3 
        bb = -1.0_dp
       !cc = ac_mass + ita*Cfcr*((d_ac_routes(j1,props_time)-d_ac_routes(j1-1,props_time))*OneDay)*rho*(VTAS**2)*S_area*CD0*0.5_dp
        cc = ac_mass + ita*Cfcr*((d_ac_routes(j1,props_time)-d_ac_routes(j1-1,props_time))*OneDay)* & 
             d_pass_phyval_along_traj(j1-1,3)*(VTAS(j1-1)**2)*S_area*CD0*0.5_dp        !20130902-1,20140602-3
        !write(*,*)'aa=',aa,'bb=',bb,'cc=',cc

        sq = dsqrt(bb*bb-4.0_dp*aa*cc)
        de = 2.0_dp*aa
        x1 = (-bb+sq)/de
        x2 = -(bb+sq)/de
        !write(*,*)'x1=',x1,'x2=',x2

        !Determine a solution (x1 or x2) as new_ac_mass.
        IF(sq.gt.0.0_dp)THEN    !D>0: Two solution exist
           IF(x1.ge.0.0_dp .and. x2.ge.0.0_dp)THEN
              IF(x1.ge.ac_mass .and. x2.ge.ac_mass)THEN
                 IF(dabs(x1-ac_mass) .gt. dabs(x2-ac_mass))new_ac_mass = x2
                 IF(dabs(x1-ac_mass) .lt. dabs(x2-ac_mass))new_ac_mass = x1
              ELSE
                 IF(x1.lt.ac_mass .and. x2.ge.ac_mass)new_ac_mass = x2
                 IF(x2.lt.ac_mass .and. x1.ge.ac_mass)new_ac_mass = x1
              ENDIF
           ELSEIF(x1.gt.0.0_dp .and. x2.lt.0.0_dp)THEN
              new_ac_mass = x1
           ELSEIF(x1.lt.0.0_dp .and. x2.gt.0.0_dp)THEN
              new_ac_mass = x2
           ELSE
              write(*,*)'ERROR: Total_energy_model_calculation(D>0)'
           ENDIF

        ELSEIF(sq.eq.0.0_dp)THEN  !D=0: One solution exists(x1=x2)
           IF(x1.gt.0.0_dp)new_ac_mass = x1
           IF(x1.le.0.0_dp)write(*,*)'ERROR: Total_energy_model_calculation(D=0)'
        ENDIF

        !Researve Fuel_use value on current segment and exchange mass values 
        !for calculation at next waypoint.
        d_ac_routes(j1,props_fuel_use) = new_ac_mass - ac_mass  ![kg],see20130711-1
        ac_mass     = new_ac_mass                               ![kg]
        ave_ac_mass = ave_ac_mass + ac_mass                     !20141216-4 calculate average value of ac_mass [kg]
!        if(d_p_pe.eq.0 .and. fl_direction.eq.0)write(71,*)d_i_traj_dep, ac_mass, fl_direction, "Eastbound"   !20141215-2
!        if(d_p_pe.eq.0 .and. fl_direction.eq.2)write(72,*)d_i_traj_dep, ac_mass, fl_direction, "Westbound"
!        if(d_p_pe.eq.1 .and. fl_direction.eq.0)write(73,*)d_i_traj_dep, ac_mass, fl_direction, "Eastbound"
!        if(d_p_pe.eq.1 .and. fl_direction.eq.2)write(74,*)d_i_traj_dep, ac_mass, fl_direction, "Westbound"
!        if(d_p_pe.eq.2 .and. fl_direction.eq.0)write(75,*)d_i_traj_dep, ac_mass, fl_direction, "Eastbound"
!        if(d_p_pe.eq.2 .and. fl_direction.eq.2)write(76,*)d_i_traj_dep, ac_mass, fl_direction, "Westbound"
!        if(d_p_pe.eq.3 .and. fl_direction.eq.0)write(77,*)d_i_traj_dep, ac_mass, fl_direction, "Eastbound"
!        if(d_p_pe.eq.3 .and. fl_direction.eq.2)write(78,*)d_i_traj_dep, ac_mass, fl_direction, "Westbound"
!        if(d_p_pe.ge.0 .and. d_p_pe.le.5 .and. fl_direction.eq.0)write(71,*)d_i_traj_dep, ac_mass, fl_direction, "Eastbound"
!        if(d_p_pe.ge.0 .and. d_p_pe.le.5 .and. fl_direction.eq.2)write(72,*)d_i_traj_dep, ac_mass, fl_direction, "Westbound"
!        if(d_p_pe.ge.6 .and. d_p_pe.le.11 .and. fl_direction.eq.0)write(73,*)d_i_traj_dep, ac_mass, fl_direction, "Eastbound"
!        if(d_p_pe.ge.6 .and. d_p_pe.le.11 .and. fl_direction.eq.2)write(74,*)d_i_traj_dep, ac_mass, fl_direction, "Westbound"
!        if(d_p_pe.ge.12 .and. d_p_pe.le.17 .and. fl_direction.eq.0)write(75,*)d_i_traj_dep, ac_mass, fl_direction, "Eastbound"
!        if(d_p_pe.ge.12 .and. d_p_pe.le.17 .and. fl_direction.eq.2)write(76,*)d_i_traj_dep, ac_mass, fl_direction, "Westbound"
!        if(d_p_pe.ge.18 .and. d_p_pe.le.23 .and. fl_direction.eq.0)write(77,*)d_i_traj_dep, ac_mass, fl_direction, "Eastbound"
!        if(d_p_pe.ge.18 .and. d_p_pe.le.23 .and. fl_direction.eq.2)write(78,*)d_i_traj_dep, ac_mass, fl_direction, "Westbound"
!       if(d_p_pe.eq.0)write(10,*)d_i_traj_dep, ac_mass
!       if(d_p_pe.eq.1)write(20,*)d_i_traj_dep+460, ac_mass
!       if(d_p_pe.eq.2)write(30,*)d_i_traj_dep+920, ac_mass
!       if(d_p_pe.eq.3)write(40,*)d_i_traj_dep+1380, ac_mass
     ENDDO     !Backward calculation loop end
     !After this DO loop, ac_mass shows the as_mass at D_city. 
!     if(d_p_pe.eq.0 .and. fl_direction.eq.0)write(71,*)          !20131127-1,20141215-2
!     if(d_p_pe.eq.0 .and. fl_direction.eq.0)write(71,*)    
!     if(d_p_pe.eq.0 .and. fl_direction.eq.2)write(72,*)    
!     if(d_p_pe.eq.0 .and. fl_direction.eq.2)write(72,*)    
!      if(d_p_pe.ge.0 .and. d_p_pe.le.5 .and. fl_direction.eq.0)write(71,*)
!      if(d_p_pe.ge.0 .and. d_p_pe.le.5 .and. fl_direction.eq.0)write(71,*)
!      if(d_p_pe.ge.0 .and. d_p_pe.le.5 .and. fl_direction.eq.2)write(72,*)
!      if(d_p_pe.ge.0 .and. d_p_pe.le.5 .and. fl_direction.eq.2)write(72,*)
!     if(d_p_pe.eq.1 .and. fl_direction.eq.0)write(73,*)      
!     if(d_p_pe.eq.1 .and. fl_direction.eq.0)write(73,*)    
!     if(d_p_pe.eq.1 .and. fl_direction.eq.2)write(74,*)    
!     if(d_p_pe.eq.1 .and. fl_direction.eq.2)write(74,*)    
!      if(d_p_pe.ge.6 .and. d_p_pe.le.11 .and. fl_direction.eq.0)write(73,*)
!      if(d_p_pe.ge.6 .and. d_p_pe.le.11 .and. fl_direction.eq.0)write(73,*)
!      if(d_p_pe.ge.6 .and. d_p_pe.le.11 .and. fl_direction.eq.2)write(74,*)
!      if(d_p_pe.ge.6 .and. d_p_pe.le.11 .and. fl_direction.eq.2)write(74,*)
!     if(d_p_pe.eq.2 .and. fl_direction.eq.0)write(75,*)      
!     if(d_p_pe.eq.2 .and. fl_direction.eq.0)write(75,*)    
!     if(d_p_pe.eq.2 .and. fl_direction.eq.2)write(76,*)    
!     if(d_p_pe.eq.2 .and. fl_direction.eq.2)write(76,*)    
!      if(d_p_pe.ge.12 .and. d_p_pe.le.17 .and. fl_direction.eq.0)write(75,*)
!      if(d_p_pe.ge.12 .and. d_p_pe.le.17 .and. fl_direction.eq.0)write(75,*)
!      if(d_p_pe.ge.12 .and. d_p_pe.le.17 .and. fl_direction.eq.2)write(76,*)
!      if(d_p_pe.ge.12 .and. d_p_pe.le.17 .and. fl_direction.eq.2)write(76,*)
!     if(d_p_pe.eq.3 .and. fl_direction.eq.0)write(77,*)      
!     if(d_p_pe.eq.3 .and. fl_direction.eq.0)write(77,*)    
!     if(d_p_pe.eq.3 .and. fl_direction.eq.2)write(78,*)    
!     if(d_p_pe.eq.3 .and. fl_direction.eq.2)write(78,*)    
!      if(d_p_pe.ge.18 .and. d_p_pe.le.23 .and. fl_direction.eq.0)write(77,*)
!      if(d_p_pe.ge.18 .and. d_p_pe.le.23 .and. fl_direction.eq.0)write(77,*)
!      if(d_p_pe.ge.18 .and. d_p_pe.le.23 .and. fl_direction.eq.2)write(78,*)
!      if(d_p_pe.ge.18 .and. d_p_pe.le.23 .and. fl_direction.eq.2)write(78,*)
!     close(71)   !20131115-1,20141215-2 
!     close(72)   
!     close(73)  
!     close(74)  
!     close(75)  
!     close(76)  
!     close(77)  
!     close(78)  

     !statistics of ac_mass output, 20140926-2,3, 20141215-3,20141215-5,20141216-3,4
!     IF(d_p_pe.eq.0)THEN
!     IF((d_p_pe.ge.0).and.(d_p_pe.le.5))THEN
!        IF(fl_direction.eq.0)THEN       !Eastbound
!           OPEN(41,file='sta_ac_mass_pe0_eb.dat',position='append')                           !20141215-2,20141216-3,4,20160606-7
!              write(41,*)d_i_traj_dep, ac_mass, arr_ac_mass, ave_ac_mass/nwaypoints, fl_direction, "Eastbound" !ac_mass: ac_mass at D_city
!           CLOSE(41)                                                                          !arr_ac_mass: ac_mass at A_city
!           OPEN(81,file='total_fuel_use_pe0_eb.dat',position='append')           !20141215-5  !ave_ac_mass: ave. of ac_mass 
!              write(81,*)d_i_traj_dep, ac_mass-arr_ac_mass, fl_direction, "Eastbound"         
!           CLOSE(81)
!        ELSEIF(fl_direction.eq.2)THEN   !Westbound
!           OPEN(42,file='sta_ac_mass_pe0_wb.dat',position='append')
!              write(42,*)d_i_traj_dep, ac_mass, arr_ac_mass, ave_ac_mass/nwaypoints, fl_direction, "Westbound"      
!           CLOSE(42)
!           OPEN(82,file='total_fuel_use_pe0_wb.dat',position='append')           !20141215-5
!              write(82,*)d_i_traj_dep, ac_mass-arr_ac_mass, fl_direction, "Westbound"         
!           CLOSE(82)
!        ENDIF
!     ENDIF
!!     IF(d_p_pe.eq.1)THEN
!     IF((d_p_pe.ge.6).and.(d_p_pe.le.11))THEN
!        IF(fl_direction.eq.0)THEN       !Eastbound
!           OPEN(43,file='sta_ac_mass_pe1_eb.dat',position='append')
!              write(43,*)d_i_traj_dep, ac_mass, arr_ac_mass, ave_ac_mass/nwaypoints, fl_direction, "Eastbound"      
!           CLOSE(43)
!           OPEN(83,file='total_fuel_use_pe1_eb.dat',position='append')           !20141215-5
!              write(83,*)d_i_traj_dep, ac_mass-arr_ac_mass, fl_direction, "Eastbound"         
!           CLOSE(83)
!        ELSEIF(fl_direction.eq.2)THEN   !Westbound
!           OPEN(44,file='sta_ac_mass_pe1_wb.dat',position='append')
!              write(44,*)d_i_traj_dep, ac_mass, arr_ac_mass, ave_ac_mass/nwaypoints, fl_direction, "Westbound"    
!           CLOSE(44)
!           OPEN(84,file='total_fuel_use_pe1_wb.dat',position='append')           !20141215-5
!              write(84,*)d_i_traj_dep, ac_mass-arr_ac_mass, fl_direction, "Westbound"         
!           CLOSE(84)
!        ENDIF
!     ENDIF
!!     IF(d_p_pe.eq.2)THEN
!     IF((d_p_pe.ge.12).and.(d_p_pe.le.17))THEN
!        IF(fl_direction.eq.0)THEN       !Eastbound
!           OPEN(45,file='sta_ac_mass_pe2_eb.dat',position='append')
!              write(45,*)d_i_traj_dep, ac_mass, arr_ac_mass, ave_ac_mass/nwaypoints, fl_direction, "Eastbound"      
!           CLOSE(45)
!           OPEN(85,file='total_fuel_use_pe2_eb.dat',position='append')           !20141215-5
!              write(85,*)d_i_traj_dep, ac_mass-arr_ac_mass, fl_direction, "Eastbound"         
!           CLOSE(85)
!        ELSEIF(fl_direction.eq.2)THEN   !Westbound
!           OPEN(46,file='sta_ac_mass_pe2_wb.dat',position='append')
!              write(46,*)d_i_traj_dep, ac_mass, arr_ac_mass, ave_ac_mass/nwaypoints, fl_direction, "Westbound"      
!           CLOSE(46)
!           OPEN(86,file='total_fuel_use_pe2_wb.dat',position='append')           !20141215-5
!              write(86,*)d_i_traj_dep, ac_mass-arr_ac_mass, fl_direction, "Westbound"         
!           CLOSE(86)
!        ENDIF
!     ENDIF
!!     IF(d_p_pe.eq.3)THEN
!     IF((d_p_pe.ge.18).and.(d_p_pe.le.23))THEN
!        IF(fl_direction.eq.0)THEN       !Eastbound
!           OPEN(47,file='sta_ac_mass_pe3_eb.dat',position='append')
!              write(47,*)d_i_traj_dep, ac_mass, arr_ac_mass, ave_ac_mass/nwaypoints, fl_direction, "Eastbound"       
!           CLOSE(47)
!           OPEN(87,file='total_fuel_use_pe3_eb.dat',position='append')           !20141215-5
!              write(87,*)d_i_traj_dep, ac_mass-arr_ac_mass, fl_direction, "Eastbound"         
!           CLOSE(87)
!        ELSEIF(fl_direction.eq.2)THEN   !Westbound
!           OPEN(48,file='sta_ac_mass_pe3_wb.dat',position='append')
!              write(48,*)d_i_traj_dep, ac_mass, arr_ac_mass, ave_ac_mass/nwaypoints, fl_direction, "Westbound"      
!           CLOSE(48)
!           OPEN(88,file='total_fuel_use_pe3_wb.dat',position='append')           !20141215-5
!              write(88,*)d_i_traj_dep, ac_mass-arr_ac_mass, fl_direction, "Westbound"         
!           CLOSE(88)
!        ENDIF
!     ENDIF


!--------------------------------------------------------
!    DLR Fuel Flow Method: correct Emission Index
!--------------------------------------------------------
     ! EINOx is calculated by DLR Fuel Flow Method for each segment. 
     ! see 20130524-4,20130527-3,20130727-1,20130802-1 
     ! 
     !   - 1. Approximation fuction is calculated by the Least Squares method (2nd order) 
     !        + Cholesky method (see20130803-2), with ICAO DATA (4 points). The definition of 
     !        the function is:
     !           EINOx,ref = coefficient(1) + coefficient(2)*Wfuel_flow,ref + coefficient(3)*Wfuel_flow,ref**2
     !           see20130803-2,20130805-1 
     !        NOTE: ICAO data(WFuel_flow,ref and EINOx,ref) should be reffered from Aircraft table.
     !   - 2. EINOx is calculated by the DLR fuel flow method. 
     !
     ! NOTE: The approximation function calculation is performed only once (at first usage).
     !       This function is specific for each engine. This function is NOT affected by meteorological
     !       conditions, and difference in position along waypoints.see20130805-4


     ! START: The Least Squares method (2nd order),     see20130803-2
     ! ------W_array and z_array calculations
!     write(*,*)'ICAO_DATA=',ICAO_DATA(:,:)
     w_array(1,1)=num_icao_data
     DO 10 j=2,order_app_equ+1
        w_array(j,1)=0.0_dp
        DO 20 k=1,num_icao_data
           w_array(j,1)=w_array(j,1)+ICAO_DATA(1,k)**(j-1)
        20 continue      
        w_array(order_app_equ+1,j)=0.0_dp
        DO 30 k=1,num_icao_data
           w_array(order_app_equ+1,j)=w_array(order_app_equ+1,j)+ICAO_DATA(1,k)**(j+order_app_equ-1)
        30 continue
        DO 40 i=1,j-1
           w_array(i,j-i+1)=w_array(j,1)
        40 continue      
        DO 50 i=1,order_app_equ+1-j
           w_array(order_app_equ+1-i,j+i)=w_array(order_app_equ+1,j)
        50 continue
     10 continue
     DO 60 j=1,order_app_equ+1
        z_array(j)=0.0_dp
        DO 70 k=1,num_icao_data
           IF(j.eq.1)THEN
              z_array(j)=z_array(j)+ICAO_DATA(2,k)
           ELSE
              z_array(j)=z_array(j)+(ICAO_DATA(1,k)**(j-1))*ICAO_DATA(2,k)
           ENDIF
        70 continue
     60 continue

     ! ------Wa=z calculation by Cholesky method
     ! START: The Cholesky method
     ! ------Trans(R)DR calculation
     !write(*,*)'coefficient_check0=',coefficient(:)
     DO 80 i=1,order_app_equ+1
        s_chol = 0.0_dp
        DO 90 k=1,i-1
           s_chol=s_chol+w_array(k,i)*w_array(k,i)*w_array(k,k)
        90 continue
        w_array(i,i)=w_array(i,i)-s_chol
        DO 110 j=i+1,order_app_equ+1
           s_chol=0.0_dp
           DO 100 k=1,i-1
              s_chol=s_chol+w_array(k,i)*w_array(k,j)*w_array(k,k)
           100 continue
           w_array(i,j)=(w_array(i,j)-s_chol)/w_array(i,i)
        110 continue
     80 continue 
     ! ------Trans(R)y=b calculation
     DO 130 i=2,order_app_equ+1
        s_chol=0.0_dp
        DO 120 j=1,i-1
           s_chol=s_chol+w_array(j,i)*z_array(j)
        120 continue       
        z_array(i)=z_array(i)-s_chol
     130 continue
     ! ------DRx=y calculation
     DO 150 i=order_app_equ+1,1,-1
        !write(*,*)'coefficient_check1=',i,coefficient(:)
        s_chol=0.0_dp
        DO 140 j=i+1,order_app_equ+1
           !write(*,*)'coefficient_check2=',coefficient(j)
           s_chol=s_chol+w_array(i,j)*coefficient(j)
        140 continue
        coefficient(i)=z_array(i)/w_array(i,i)-s_chol
     150 continue       
     !write(*,*)'coefficient=',coefficient(:)
     ! END: The Cholsky method 
     ! END: The Least Squares method (2nd order)


     ! NOTE: Pressure and temperature values should be expressed by dimension(:) along waypoints.see20130727-1
     !       The EINOx is considered as a value of one engine. However, consequently the EINOx value of one aircraft is       
     !       identical to that of one engine. 20130802-2,20131130-2
     ! Wfuel        : fuel consumption for one engine, [kg(fuel)/s].
     ! Wfuel_ref    : fuel consumption (reference) for one engine, [kg(fuel)/s].
     ! EINOx        : NOx emission index for one engine, [g(NOx)/kg(fuel)].
     ! EINOx_ref    : NOx emission index (reference) for one engine, [g(NOx)/kg(fuel)].
!     open(11,file="fuel_cons_pe0_eb.dat",position='APPEND')   !20131125-1,20141215-5
!     open(12,file="fuel_cons_pe0_wb.dat",position='APPEND')   
!     open(13,file="fuel_cons_pe1_eb.dat",position='APPEND')   
!     open(14,file="fuel_cons_pe1_wb.dat",position='APPEND')   
!     open(15,file="fuel_cons_pe2_eb.dat",position='APPEND')   
!     open(16,file="fuel_cons_pe2_wb.dat",position='APPEND')   
!     open(17,file="fuel_cons_pe3_eb.dat",position='APPEND')   
!     open(18,file="fuel_cons_pe3_wb.dat",position='APPEND')   

!     open(21,file="nox_pe0_eb.dat",position='APPEND')   !20131130-3,20141215-6
!     open(22,file="nox_pe0_wb.dat",position='APPEND')   
!     open(23,file="nox_pe1_eb.dat",position='APPEND')   
!     open(24,file="nox_pe1_wb.dat",position='APPEND')   
!     open(25,file="nox_pe2_eb.dat",position='APPEND')   
!     open(26,file="nox_pe2_wb.dat",position='APPEND')   
!     open(27,file="nox_pe3_eb.dat",position='APPEND')   
!     open(28,file="nox_pe3_wb.dat",position='APPEND')   
     DO j2=1,nwaypoints-1     !20140603-1
        Wfuel(j2)         = d_ac_routes(j2+1,props_fuel_use) / &
             ((d_ac_routes(j2+1,props_time)-d_ac_routes(j2,props_time)) &
             *OneDay)/num_engine 
!       if(d_p_pe.eq.0)write(11,*)d_i_traj_dep, Wfuel(j2)*num_engine*60.0_dp      !Fuel consumption of one aircraft [kg(fuel)/min]
!       if(d_p_pe.eq.1)write(12,*)d_i_traj_dep+460, Wfuel(j2)*num_engine*60.0_dp  !20131125-1  
!       if(d_p_pe.eq.2)write(13,*)d_i_traj_dep+920, Wfuel(j2)*num_engine*60.0_dp  
!       if(d_p_pe.eq.3)write(14,*)d_i_traj_dep+1380, Wfuel(j2)*num_engine*60.0_dp
        !20141215-5
!        if(d_p_pe.eq.0 .and. fl_direction.eq.0)write(11,*)d_i_traj_dep, Wfuel(j2)*num_engine*60.0_dp, fl_direction, "Eastbound" 
!        if(d_p_pe.eq.0 .and. fl_direction.eq.2)write(12,*)d_i_traj_dep, Wfuel(j2)*num_engine*60.0_dp, fl_direction, "Westbound"
!        if(d_p_pe.eq.1 .and. fl_direction.eq.0)write(13,*)d_i_traj_dep, Wfuel(j2)*num_engine*60.0_dp, fl_direction, "Eastbound"
!        if(d_p_pe.eq.1 .and. fl_direction.eq.2)write(14,*)d_i_traj_dep, Wfuel(j2)*num_engine*60.0_dp, fl_direction, "Westbound"
!        if(d_p_pe.eq.2 .and. fl_direction.eq.0)write(15,*)d_i_traj_dep, Wfuel(j2)*num_engine*60.0_dp, fl_direction, "Eastbound"
!        if(d_p_pe.eq.2 .and. fl_direction.eq.2)write(16,*)d_i_traj_dep, Wfuel(j2)*num_engine*60.0_dp, fl_direction, "Westbound"
!        if(d_p_pe.eq.3 .and. fl_direction.eq.0)write(17,*)d_i_traj_dep, Wfuel(j2)*num_engine*60.0_dp, fl_direction, "Eastbound"
!        if(d_p_pe.eq.3 .and. fl_direction.eq.2)write(18,*)d_i_traj_dep, Wfuel(j2)*num_engine*60.0_dp, fl_direction, "Westbound"
!    if(d_p_pe.ge.0 .and. d_p_pe.le.5 .and. fl_direction.eq.0)write(11,*)d_i_traj_dep,Wfuel(j2)*num_engine*60.0_dp,fl_direction,"Eastbound" 
!    if(d_p_pe.ge.0 .and. d_p_pe.le.5 .and. fl_direction.eq.2)write(12,*)d_i_traj_dep,Wfuel(j2)*num_engine*60.0_dp,fl_direction,"Westbound"
!    if(d_p_pe.ge.6 .and. d_p_pe.le.11 .and. fl_direction.eq.0)write(13,*)d_i_traj_dep,Wfuel(j2)*num_engine*60.0_dp,fl_direction,"Eastbound"
!    if(d_p_pe.ge.6 .and. d_p_pe.le.11 .and. fl_direction.eq.2)write(14,*)d_i_traj_dep,Wfuel(j2)*num_engine*60.0_dp,fl_direction,"Westbound"
!   if(d_p_pe.ge.12 .and. d_p_pe.le.17 .and. fl_direction.eq.0)write(15,*)d_i_traj_dep,Wfuel(j2)*num_engine*60.0_dp,fl_direction,"Eastbound"
!   if(d_p_pe.ge.12 .and. d_p_pe.le.17 .and. fl_direction.eq.2)write(16,*)d_i_traj_dep,Wfuel(j2)*num_engine*60.0_dp,fl_direction,"Westbound"
!   if(d_p_pe.ge.18 .and. d_p_pe.le.23 .and. fl_direction.eq.0)write(17,*)d_i_traj_dep,Wfuel(j2)*num_engine*60.0_dp,fl_direction,"Eastbound"
!   if(d_p_pe.ge.18 .and. d_p_pe.le.23 .and. fl_direction.eq.2)write(18,*)d_i_traj_dep,Wfuel(j2)*num_engine*60.0_dp,fl_direction,"Westbound"

        pressure_total    = d_pass_phyval_along_traj(j2,1)*(1.0_dp+0.2_dp*(AC_MACH**2))**3.5 !AC_MACH should be from Aircraft table 
        temperature_total = d_pass_phyval_along_traj(j2,2)*(1.0_dp+0.2_dp*(AC_MACH**2))      !AC_MACH should be from Aircraft table
        delta_total       = pressure_total/atm2Pa
        theta_total       = temperature_total/temp_sealv                         !20140724-2
        Wfuel_ref(j2)     = Wfuel(j2)*(delta_total**(-1))*(theta_total**(-0.5))
!        write(*,*)'num_engine=',num_engine,'AC_MACH=',AC_MACH,'atm2Pa=',atm2Pa,'temp_sealv=',temp_sealv
!        write(*,*)'Wfuel_ref=',Wfuel_ref(j2),'Wfuel=',Wfuel(j2)

        EINOx_ref(j2)   = coefficient(1) + coefficient(2)*Wfuel_ref(j2) + coefficient(3)*(Wfuel_ref(j2)**2)
        specific_hmd    = (10.0_dp**(-3))*DEXP(-0.0001426_dp*(d_ac_routes(j2,props_alt)*3.280839895013_dp-12900.0_dp))
        hmd_corr_factor = -19.0_dp*(specific_hmd-0.00634_dp)
        EINOx(j2)       = EINOx_ref(j2)*(delta_total**0.4)*(theta_total**3)*DEXP(hmd_corr_factor) !*num_engine,see20130802-2
!        write(*,*)'EINOx_ref=',EINOx_ref(j2),'EINOx=',EINOx(j2),'specific_hmd:',specific_hmd      !see20131130-2 
!        if(d_p_pe.eq.0 .and. fl_direction.eq.0)write(21,*)d_i_traj_dep, EINOx(j2), fl_direction, "Eastbound"   !20141215-6 
!        if(d_p_pe.eq.0 .and. fl_direction.eq.2)write(22,*)d_i_traj_dep, EINOx(j2), fl_direction, "Westbound"
!        if(d_p_pe.eq.1 .and. fl_direction.eq.0)write(23,*)d_i_traj_dep, EINOx(j2), fl_direction, "Eastbound"
!        if(d_p_pe.eq.1 .and. fl_direction.eq.2)write(24,*)d_i_traj_dep, EINOx(j2), fl_direction, "Westbound"
!        if(d_p_pe.eq.2 .and. fl_direction.eq.0)write(25,*)d_i_traj_dep, EINOx(j2), fl_direction, "Eastbound"
!        if(d_p_pe.eq.2 .and. fl_direction.eq.2)write(26,*)d_i_traj_dep, EINOx(j2), fl_direction, "Westbound"
!        if(d_p_pe.eq.3 .and. fl_direction.eq.0)write(27,*)d_i_traj_dep, EINOx(j2), fl_direction, "Eastbound"
!        if(d_p_pe.eq.3 .and. fl_direction.eq.2)write(28,*)d_i_traj_dep, EINOx(j2), fl_direction, "Westbound"
!    if(d_p_pe.ge.0 .and. d_p_pe.le.5 .and. fl_direction.eq.0)write(21,*)d_i_traj_dep, EINOx(j2), fl_direction, "Eastbound"   !20141215-6 
!    if(d_p_pe.ge.0 .and. d_p_pe.le.5 .and. fl_direction.eq.2)write(22,*)d_i_traj_dep, EINOx(j2), fl_direction, "Westbound"
!    if(d_p_pe.ge.6 .and. d_p_pe.le.11 .and. fl_direction.eq.0)write(23,*)d_i_traj_dep, EINOx(j2), fl_direction, "Eastbound"
!    if(d_p_pe.ge.6 .and. d_p_pe.le.11 .and. fl_direction.eq.2)write(24,*)d_i_traj_dep, EINOx(j2), fl_direction, "Westbound"
!    if(d_p_pe.ge.12 .and. d_p_pe.le.17 .and. fl_direction.eq.0)write(25,*)d_i_traj_dep, EINOx(j2), fl_direction, "Eastbound"
!    if(d_p_pe.ge.12 .and. d_p_pe.le.17 .and. fl_direction.eq.2)write(26,*)d_i_traj_dep, EINOx(j2), fl_direction, "Westbound"
!    if(d_p_pe.ge.18 .and. d_p_pe.le.23 .and. fl_direction.eq.0)write(27,*)d_i_traj_dep, EINOx(j2), fl_direction, "Eastbound"
!    if(d_p_pe.ge.18 .and. d_p_pe.le.23 .and. fl_direction.eq.2)write(28,*)d_i_traj_dep, EINOx(j2), fl_direction, "Westbound"
     ENDDO
!     if(d_p_pe.eq.0 .and. fl_direction.eq.0)write(11,*)        !20141215-5 
!     if(d_p_pe.eq.0 .and. fl_direction.eq.0)write(11,*)    
!     if(d_p_pe.eq.0 .and. fl_direction.eq.2)write(12,*)    
!     if(d_p_pe.eq.0 .and. fl_direction.eq.2)write(12,*)    
!     if(d_p_pe.eq.1 .and. fl_direction.eq.0)write(13,*)      
!     if(d_p_pe.eq.1 .and. fl_direction.eq.0)write(13,*)    
!     if(d_p_pe.eq.1 .and. fl_direction.eq.2)write(14,*)    
!     if(d_p_pe.eq.1 .and. fl_direction.eq.2)write(14,*)    
!     if(d_p_pe.eq.2 .and. fl_direction.eq.0)write(15,*)      
!     if(d_p_pe.eq.2 .and. fl_direction.eq.0)write(15,*)    
!     if(d_p_pe.eq.2 .and. fl_direction.eq.2)write(16,*)    
!     if(d_p_pe.eq.2 .and. fl_direction.eq.2)write(16,*)    
!     if(d_p_pe.eq.3 .and. fl_direction.eq.0)write(17,*)      
!     if(d_p_pe.eq.3 .and. fl_direction.eq.0)write(17,*)    
!     if(d_p_pe.eq.3 .and. fl_direction.eq.2)write(18,*)    
!     if(d_p_pe.eq.3 .and. fl_direction.eq.2)write(18,*)    
!      if(d_p_pe.ge.0 .and. d_p_pe.le.5 .and. fl_direction.eq.0)write(11,*)
!      if(d_p_pe.ge.0 .and. d_p_pe.le.5 .and. fl_direction.eq.0)write(11,*)
!      if(d_p_pe.ge.0 .and. d_p_pe.le.5 .and. fl_direction.eq.2)write(12,*)
!      if(d_p_pe.ge.0 .and. d_p_pe.le.5 .and. fl_direction.eq.2)write(12,*)
!      if(d_p_pe.ge.6 .and. d_p_pe.le.11 .and. fl_direction.eq.0)write(13,*)
!      if(d_p_pe.ge.6 .and. d_p_pe.le.11 .and. fl_direction.eq.0)write(13,*)
!      if(d_p_pe.ge.6 .and. d_p_pe.le.11 .and. fl_direction.eq.2)write(14,*)
!      if(d_p_pe.ge.6 .and. d_p_pe.le.11 .and. fl_direction.eq.2)write(14,*)
!      if(d_p_pe.ge.12 .and. d_p_pe.le.17 .and. fl_direction.eq.0)write(15,*)
!      if(d_p_pe.ge.12 .and. d_p_pe.le.17 .and. fl_direction.eq.0)write(15,*)
!      if(d_p_pe.ge.12 .and. d_p_pe.le.17 .and. fl_direction.eq.2)write(16,*)
!      if(d_p_pe.ge.12 .and. d_p_pe.le.17 .and. fl_direction.eq.2)write(16,*)
!      if(d_p_pe.ge.18 .and. d_p_pe.le.23 .and. fl_direction.eq.0)write(17,*)
!      if(d_p_pe.ge.18 .and. d_p_pe.le.23 .and. fl_direction.eq.0)write(17,*)
!      if(d_p_pe.ge.18 .and. d_p_pe.le.23 .and. fl_direction.eq.2)write(18,*)
!      if(d_p_pe.ge.18 .and. d_p_pe.le.23 .and. fl_direction.eq.2)write(18,*)

!     if(d_p_pe.eq.0)write(21,*)          !20131130-3,20141215-6
!     if(d_p_pe.eq.0)write(21,*)
!     if(d_p_pe.eq.1)write(22,*)
!     if(d_p_pe.eq.1)write(22,*)
!     if(d_p_pe.eq.2)write(23,*)
!     if(d_p_pe.eq.2)write(23,*)
!     if(d_p_pe.eq.3)write(24,*)
!     if(d_p_pe.eq.3)write(24,*)
!     if(d_p_pe.eq.0 .and. fl_direction.eq.0)write(21,*)        !20141215-5 
!     if(d_p_pe.eq.0 .and. fl_direction.eq.0)write(21,*)    
!     if(d_p_pe.eq.0 .and. fl_direction.eq.2)write(22,*)    
!     if(d_p_pe.eq.0 .and. fl_direction.eq.2)write(22,*)    
!     if(d_p_pe.eq.1 .and. fl_direction.eq.0)write(23,*)      
!     if(d_p_pe.eq.1 .and. fl_direction.eq.0)write(23,*)    
!     if(d_p_pe.eq.1 .and. fl_direction.eq.2)write(24,*)    
!     if(d_p_pe.eq.1 .and. fl_direction.eq.2)write(24,*)    
!     if(d_p_pe.eq.2 .and. fl_direction.eq.0)write(25,*)      
!     if(d_p_pe.eq.2 .and. fl_direction.eq.0)write(25,*)    
!     if(d_p_pe.eq.2 .and. fl_direction.eq.2)write(26,*)    
!     if(d_p_pe.eq.2 .and. fl_direction.eq.2)write(26,*)    
!     if(d_p_pe.eq.3 .and. fl_direction.eq.0)write(27,*)      
!     if(d_p_pe.eq.3 .and. fl_direction.eq.0)write(27,*)    
!     if(d_p_pe.eq.3 .and. fl_direction.eq.2)write(28,*)    
!     if(d_p_pe.eq.3 .and. fl_direction.eq.2)write(28,*)    
!      if(d_p_pe.ge.0 .and. d_p_pe.le.5 .and. fl_direction.eq.0)write(21,*)
!      if(d_p_pe.ge.0 .and. d_p_pe.le.5 .and. fl_direction.eq.0)write(21,*)
!      if(d_p_pe.ge.0 .and. d_p_pe.le.5 .and. fl_direction.eq.2)write(22,*)
!      if(d_p_pe.ge.0 .and. d_p_pe.le.5 .and. fl_direction.eq.2)write(22,*)
!      if(d_p_pe.ge.6 .and. d_p_pe.le.11 .and. fl_direction.eq.0)write(23,*)
!      if(d_p_pe.ge.6 .and. d_p_pe.le.11 .and. fl_direction.eq.0)write(23,*)
!      if(d_p_pe.ge.6 .and. d_p_pe.le.11 .and. fl_direction.eq.2)write(24,*)
!      if(d_p_pe.ge.6 .and. d_p_pe.le.11 .and. fl_direction.eq.2)write(24,*)
!      if(d_p_pe.ge.12 .and. d_p_pe.le.17 .and. fl_direction.eq.0)write(25,*)
!      if(d_p_pe.ge.12 .and. d_p_pe.le.17 .and. fl_direction.eq.0)write(25,*)
!      if(d_p_pe.ge.12 .and. d_p_pe.le.17 .and. fl_direction.eq.2)write(26,*)
!      if(d_p_pe.ge.12 .and. d_p_pe.le.17 .and. fl_direction.eq.2)write(26,*)
!      if(d_p_pe.ge.18 .and. d_p_pe.le.23 .and. fl_direction.eq.0)write(27,*)
!      if(d_p_pe.ge.18 .and. d_p_pe.le.23 .and. fl_direction.eq.0)write(27,*)
!      if(d_p_pe.ge.18 .and. d_p_pe.le.23 .and. fl_direction.eq.2)write(28,*)
!      if(d_p_pe.ge.18 .and. d_p_pe.le.23 .and. fl_direction.eq.2)write(28,*)

!     close(11)   !20131125-1,20141215-6
!     close(12)    
!     close(13)    
!     close(14)    
!     close(15)    
!     close(16)    
!     close(17)    
!     close(18)    

!     close(21)   !20131130-3,20141215-6
!     close(22)    
!     close(23)    
!     close(24)    
!     close(25)    
!     close(26)    
!     close(27)    
!     close(28)    

!-----------------------------------------
!    Cost and COC calculation
!-----------------------------------------
     SELECT CASE(option_traj_calc)                           !20170428-2,-6
     CASE(GC)
        coc_total_flight_time = gc_opt_fl_time               ![s],20170428-2    
        coc_total_fuel_use    = ac_mass - arr_ac_mass        ![kg(fuel)],20170428-2  
        CALL coc_calc(coc_total_flight_time, coc_total_fuel_use, coc_value) 
        CALL cost_calc(coc_total_flight_time, coc_total_fuel_use, cost_value)     !20170428-6 
!        write(*,*)'RO_coc1=',option_traj_calc,'PE=',d_p_pe,'FD=',fl_direction,     & 
!                  'FT=',coc_total_flight_time,'Fuel=',coc_total_fuel_use,'COST=',cost_value,'COC=',coc_value

     CASE(Wind_opt, Fuel_opt, NOx_opt, H2O_opt, Cost_opt, COC_opt, ATR20_opt)     !20180511-3       
        coc_total_flight_time = wind_opt_fl_time             ![s],20170428-2    
        coc_total_fuel_use    = ac_mass - arr_ac_mass        ![kg(fuel)],20170428-2 
        CALL coc_calc(coc_total_flight_time, coc_total_fuel_use, coc_value) 
        CALL cost_calc(coc_total_flight_time, coc_total_fuel_use, cost_value)     !20170428-6 
!        write(*,*)'RO_coc2=',option_traj_calc,'PE=',d_p_pe,'FD=',fl_direction,     & 
!                  'FT=',coc_total_flight_time,'Fuel=',coc_total_fuel_use,'COST=',cost_value,'COC=',coc_value

     CASE(ContrailPC_opt)
        coc_total_flight_time = contrail_opt_fl_time   !wind_opt_fl_time             ![s],20170428-2    
        coc_total_fuel_use    = ac_mass - arr_ac_mass          ![kg(fuel)],20170428-2 
        CALL coc_calc(coc_total_flight_time, coc_total_fuel_use, coc_value) 
        CALL cost_calc(coc_total_flight_time, coc_total_fuel_use, cost_value)     !20170428-6 
!        write(*,*)'RO_coc3=',option_traj_calc,'PE=',d_p_pe,'FD=',fl_direction,     & 
!                  'FT=',coc_total_flight_time,'Fuel=',coc_total_fuel_use,'COST=',cost_value,'COC=',coc_value
     ENDSELECT

!-----------------------------------------
!    Add Emission values into ac_routes
!-----------------------------------------
     SELECT CASE(option_traj_calc)
     ! All options
     !Yin_20170423 
     CASE(GC, Wind_opt, Fuel_opt, NOx_opt, H2O_opt, ContrailPC_opt, Cost_opt, COC_opt,ATR20_opt,COSTCLIM_opt,COSTCPC_opt)   !20170501-1
!        DO i=1,nwaypoints   !20180511-3
           d_zalt = 0.0_dp
           DO i = 1,nwaypoints
              traj_for_nnindex_d(i,1,1) = d_ac_routes(i,props_lon) 
              traj_for_nnindex_d(i,2,1) = d_ac_routes(i,props_lat) 
              traj_for_nnindex_d(i,3,1) = d_ac_routes(i,props_alt) 
           ENDDO
           temp_XXN_d=0.0_dp   
           DO i = 1,nwaypoints
              IF(i.eq.1)j=i     !20180511-4
              IF(i.gt.1)j=i-1   !20180511-4
              IF((-180.0_dp.le.traj_for_nnindex_d(j,1,1)).and.(traj_for_nnindex_d(j,1,1).lt.0.0_dp))THEN
                 temp_XXN_d = traj_for_nnindex_d(j,1,1)+360.0_dp
                 CALL nn_index(d_philon(:),temp_XXN_d,d_j_lon)
              ELSE
                 temp_XXN_d = traj_for_nnindex_d(j,1,1)
                 CALL nn_index(d_philon(:),temp_XXN_d,d_j_lon)
              ENDIF
              CALL nn_index(d_philat(:), traj_for_nnindex_d(j,2,1), d_j_lat)
              d_zalt(:)= d_zgl_geopot_3d(d_j_lon,:,d_j_lat)/g
              CALL nn_index(d_zalt(:), traj_for_nnindex_d(j,3,1), d_j_alt) 
              IF(i.eq.1)THEN
                 d_ac_routes(i,props_emis_nox)=0.0_dp
                 d_ac_routes(i,props_emis_h2o)=0.0_dp
                 d_ac_routes(i,props_atr20o3)=0.0_dp
                 d_ac_routes(i,props_atr20ch4)=0.0_dp
                 d_ac_routes(i,props_atr20h2o)=0.0_dp
                 d_ac_routes(i,props_atr20co2)=0.0_dp
                 d_ac_routes(i,props_atr20cpc)=0.0_dp
                 d_ac_routes(i,props_atr20tot)=0.0_dp
              ELSEIF(i.gt.1)THEN
                 d_ac_routes(i,props_emis_nox) = d_ac_routes(i,props_fuel_use)*EINOx(i-1)  !d_ac_routes(i,props_emis_nox), [g(NO2)]
                 d_ac_routes(i,props_emis_h2o) = d_ac_routes(i,props_fuel_use)*H2O_INDEX   !d_ac_routes(i,props_emis_h2o), [g(H2O)]
                 d_ac_routes(i,props_atr20o3)  = d_ac_routes(i,props_emis_nox)*d_ATR20O3_g(d_j_lon,d_j_alt,d_j_lat)/1000._dp !d_ac_routes(i,props_atr20o3), [K]
                 d_ac_routes(i,props_atr20ch4) = d_ac_routes(i,props_emis_nox)*d_ATR20CH4_g(d_j_lon,d_j_alt,d_j_lat)/1000._dp
                 d_ac_routes(i,props_atr20h2o) = d_ac_routes(i,props_fuel_use)*d_ATR20H2O_g(d_j_lon,d_j_alt,d_j_lat)
                 d_ac_routes(i,props_atr20co2) = d_ac_routes(i,props_fuel_use)*d_ATR20CO2_g(d_j_lon,d_j_alt,d_j_lat)
                 d_ac_routes(i,props_atr20cpc) = d_ac_routes(i,props_potcov)*d_ATR20CPC_g(d_j_lon,d_j_alt,d_j_lat)
!                 write(*,*)'props_atr20_check:',d_ac_routes(i,props_potcov),d_ATR20CPC_g(d_j_lon,d_j_alt,d_j_lat),d_j_lon,d_j_alt,d_j_lat
                 d_ac_routes(i,props_atr20tot) = d_ac_routes(i,props_atr20o3)+d_ac_routes(i,props_atr20ch4)+  &
                                                 d_ac_routes(i,props_atr20h2o)+d_ac_routes(i,props_atr20cpc)+ &
                                                 d_ac_routes(i,props_atr20co2)   !20180511-4  
               !20140926-3,4
                 total_nox_emission = total_nox_emission + d_ac_routes(i,props_emis_nox)
                 total_h2o_emission = total_h2o_emission + d_ac_routes(i,props_emis_h2o)
                 total_fl_distance  = total_fl_distance  + d_ac_routes(i,props_dist)
                 total_atr20_tot  = total_atr20_tot  + d_ac_routes(i,props_atr20tot)  !20180525-2                       
              ENDIF
           ENDDO
!        ENDDO       !20180511-3         
        !DUMMY 
        !d_ac_routes(:,props_emis_nox)       = 1.0_dp   !props_emis_nox=7
        !d_ac_routes(:,props_emis_h2o)       = 1.0_dp   !props_emis_h2o=8
!        write(*,*)'calculate_emissions_along_traj'
        !write(*,*)'props_dist,props_fuel_use,props_emis_nox,props_emis_h2o',props_dist,props_fuel_use,props_emis_nox,props_emis_h2o
!        write(*,*)'H2O_INDEX=',H2O_INDEX

        !20140926-3,4, 20141215-7
!        IF(d_p_pe.eq.0)THEN
!        IF((d_p_pe.ge.0).and.(d_p_pe.le.5))THEN
!           IF(fl_direction.eq.0)THEN       !Eastbound
!              OPEN(51,file='total_nox_pe0_eb.dat',position='append')
!                 write(51,*)d_i_traj_dep, total_nox_emission, fl_direction, "Eastbound"
!              CLOSE(51)
!              OPEN(61,file='total_h2o_pe0_eb.dat',position='append')
!                 write(61,*)d_i_traj_dep, total_h2o_emission, fl_direction, "Eastbound"
!              CLOSE(61)
!              OPEN(71,file='total_fl_dist_pe0_eb.dat',position='append')
!                 write(71,*)d_i_traj_dep, total_fl_distance, fl_direction, "Eastbound"
!              CLOSE(71)
!              OPEN(81,file='total_cost_pe0_eb.dat',position='append')      !20170428-3,-7
!                 write(81,*)d_p_pe, d_ac_routes(1,props_time), d_i_traj_dep, coc_total_flight_time, coc_total_fuel_use, cost_value, &
!                            coc_value, total_atr20_tot, contrail_opt_cpc, fl_direction
!              CLOSE(81)
!           ELSEIF(fl_direction.eq.2)THEN   !Westbound
!              OPEN(52,file='total_nox_pe0_wb.dat',position='append')
!                 write(52,*)d_i_traj_dep, total_nox_emission, fl_direction, "Westbound"
!              CLOSE(52)
!              OPEN(62,file='total_h2o_pe0_wb.dat',position='append')
!                 write(62,*)d_i_traj_dep, total_h2o_emission, fl_direction, "Westbound"
!              CLOSE(62)
!              OPEN(72,file='total_fl_dist_pe0_wb.dat',position='append')
!                 write(72,*)d_i_traj_dep, total_fl_distance, fl_direction, "Westbound"
!              CLOSE(72)
!              OPEN(82,file='total_cost_pe0_wb.dat',position='append')      !20170428-3,-7
!                 write(82,*)d_p_pe, d_ac_routes(1,props_time), d_i_traj_dep, coc_total_flight_time, coc_total_fuel_use, cost_value, &
!                            coc_value, total_atr20_tot, contrail_opt_cpc, fl_direction
!              CLOSE(82)
!           ENDIF
!        ENDIF
!!        IF(d_p_pe.eq.1)THEN
!        IF((d_p_pe.ge.6).and.(d_p_pe.le.11))THEN
!           IF(fl_direction.eq.0)THEN       !Eastbound
!              OPEN(53,file='total_nox_pe1_eb.dat',position='append')
!                 write(53,*)d_i_traj_dep, total_nox_emission, fl_direction, "Eastbound"
!              CLOSE(53)
!              OPEN(63,file='total_h2o_pe1_eb.dat',position='append')
!                 write(63,*)d_i_traj_dep, total_h2o_emission, fl_direction, "Eastbound"
!              CLOSE(63)
!              OPEN(73,file='total_fl_dist_pe1_eb.dat',position='append')
!                 write(73,*)d_i_traj_dep, total_fl_distance, fl_direction, "Eastbound"
!              CLOSE(73)
!              OPEN(83,file='total_cost_pe1_eb.dat',position='append')      !20170428-3,-7
!                 write(83,*)d_p_pe, d_ac_routes(1,props_time), d_i_traj_dep, coc_total_flight_time, coc_total_fuel_use, cost_value, &
!                            coc_value, total_atr20_tot, contrail_opt_cpc, fl_direction
!              CLOSE(83)
!           ELSEIF(fl_direction.eq.2)THEN   !Westbound
!              OPEN(54,file='total_nox_pe1_wb.dat',position='append')
!                 write(54,*)d_i_traj_dep, total_nox_emission, fl_direction, "Westbound"
!              CLOSE(54)
!              OPEN(64,file='total_h2o_pe1_wb.dat',position='append')
!                 write(64,*)d_i_traj_dep, total_h2o_emission, fl_direction, "Westbound"
!              CLOSE(64)
!              OPEN(74,file='total_fl_dist_pe1_wb.dat',position='append')
!                 write(74,*)d_i_traj_dep, total_fl_distance, fl_direction, "Westbound"
!              CLOSE(74)
!              OPEN(84,file='total_cost_pe1_wb.dat',position='append')      !20170428-3,-7
!                 write(84,*)d_p_pe, d_ac_routes(1,props_time), d_i_traj_dep, coc_total_flight_time, coc_total_fuel_use, cost_value, &
!                            coc_value, total_atr20_tot, contrail_opt_cpc, fl_direction
!              CLOSE(84)
!           ENDIF
!        ENDIF
!!        IF(d_p_pe.eq.2)THEN
!        IF((d_p_pe.ge.12).and.(d_p_pe.le.17))THEN
!           IF(fl_direction.eq.0)THEN       !Eastbound
!              OPEN(55,file='total_nox_pe2_eb.dat',position='append')
!                 write(55,*)d_i_traj_dep, total_nox_emission, fl_direction, "Eastbound"
!              CLOSE(55)
!              OPEN(65,file='total_h2o_pe2_eb.dat',position='append')
!                 write(65,*)d_i_traj_dep, total_h2o_emission, fl_direction, "Eastbound"
!              CLOSE(65)
!              OPEN(75,file='total_fl_dist_pe2_eb.dat',position='append')
!                 write(75,*)d_i_traj_dep, total_fl_distance, fl_direction, "Eastbound"
!              CLOSE(75)
!              OPEN(85,file='total_cost_pe2_eb.dat',position='append')      !20170428-3,-7
!                 write(85,*)d_p_pe, d_ac_routes(1,props_time), d_i_traj_dep, coc_total_flight_time, coc_total_fuel_use, cost_value, &
!                            coc_value, total_atr20_tot, contrail_opt_cpc, fl_direction
!              CLOSE(85)
!           ELSEIF(fl_direction.eq.2)THEN   !Westbound
!              OPEN(56,file='total_nox_pe2_wb.dat',position='append')
!                 write(56,*)d_i_traj_dep, total_nox_emission, fl_direction, "Westbound"
!              CLOSE(56)
!              OPEN(66,file='total_h2o_pe2_wb.dat',position='append')
!                 write(66,*)d_i_traj_dep, total_h2o_emission, fl_direction, "Westbound"
!              CLOSE(66)
!              OPEN(76,file='total_fl_dist_pe2_wb.dat',position='append')
!                 write(76,*)d_i_traj_dep, total_fl_distance, fl_direction, "Westbound"
!              CLOSE(76)
!              OPEN(86,file='total_cost_pe2_wb.dat',position='append')      !20170428-3,-7
!                 write(86,*)d_p_pe, d_ac_routes(1,props_time), d_i_traj_dep, coc_total_flight_time, coc_total_fuel_use, cost_value, &
!                            coc_value, total_atr20_tot, contrail_opt_cpc, fl_direction
!              CLOSE(86)
!           ENDIF
!        ENDIF
!!        IF(d_p_pe.eq.3)THEN
!        IF((d_p_pe.ge.18).and.(d_p_pe.le.23))THEN
!           IF(fl_direction.eq.0)THEN       !Eastbound
!              OPEN(57,file='total_nox_pe3_eb.dat',position='append')
!                 write(57,*)d_i_traj_dep, total_nox_emission, fl_direction, "Eastbound"
!              CLOSE(57)
!              OPEN(67,file='total_h2o_pe3_eb.dat',position='append')
!                 write(67,*)d_i_traj_dep, total_h2o_emission, fl_direction, "Eastbound"
!              CLOSE(67)
!              OPEN(77,file='total_fl_dist_pe3_eb.dat',position='append')
!                 write(77,*)d_i_traj_dep, total_fl_distance, fl_direction, "Eastbound"
!              CLOSE(77)
!              OPEN(87,file='total_cost_pe3_eb.dat',position='append')      !20170428-3,-7
!                 write(87,*)d_p_pe, d_ac_routes(1,props_time), d_i_traj_dep, coc_total_flight_time, coc_total_fuel_use, cost_value, &
!                            coc_value, total_atr20_tot, contrail_opt_cpc, fl_direction
!              CLOSE(87)
!           ELSEIF(fl_direction.eq.2)THEN   !Westbound
!              OPEN(58,file='total_nox_pe3_wb.dat',position='append')
!                 write(58,*)d_i_traj_dep, total_nox_emission, fl_direction, "Westbound"
!              CLOSE(58)
!              OPEN(68,file='total_h2o_pe3_wb.dat',position='append')
!                 write(68,*)d_i_traj_dep, total_h2o_emission, fl_direction, "Westbound"
!              CLOSE(68)
!              OPEN(78,file='total_fl_dist_pe3_wb.dat',position='append')
!                 write(78,*)d_i_traj_dep, total_fl_distance, fl_direction, "Westbound"
!              CLOSE(78)
!              OPEN(88,file='total_cost_pe3_wb.dat',position='append')      !20170428-3,-7
!                 write(88,*)d_p_pe, d_ac_routes(1,props_time), d_i_traj_dep, coc_total_flight_time, coc_total_fuel_use, cost_value, &
!                            coc_value, total_atr20_tot, contrail_opt_cpc, fl_direction
!              CLOSE(88)
!           ENDIF
!           !original 20140926-3,4, 20141215-7
!!           OPEN(54,file='total_nox_004.dat',position='append')
!!              write(54,*)total_nox_emission
!!           CLOSE(54)
!!           OPEN(58,file='total_h2o_004.dat',position='append')
!!              write(58,*)total_h2o_emission
!!           CLOSE(58)
!!           OPEN(64,file='total_fl_dist_004.dat',position='append')
!!              write(64,*)total_fl_distance
!!           CLOSE(64)
!        ENDIF

     !DUMMY 
     CASE DEFAULT 
        d_ac_routes(:,:)       = 1.0_dp  
!       d_p2_ac_routes_desc(:) = 1.0_dp
        write(*,*)'dummy_calculate_emissions_along_traj'
     ENDSELECT

     !Deallocation of vtas array:20140721-2
     DEALLOCATE(gc_opt_vtas_along_traj)
     DEALLOCATE(wind_opt_vtas_along_traj)
!    DEALLOCATE(fuel_opt_vtas_along_traj) !20160302-1,20160606-7
     DEALLOCATE(wind_opt_vtas)             !20150313-2
!    DEALLOCATE(fuel_opt_vtas)            !20160302-1,20160606-7
     DEALLOCATE(wind_opt_vground)          !20150313-2
     DEALLOCATE(wind_opt_rho)              !20160623-2,20170215-1
     DEALLOCATE(wind_opt_potcov)           !20170608-3
!    DEALLOCATE(fuel_opt_vground)         !20160302-1,20160606-7

     !Yin_20170423
     DEALLOCATE(contrail_opt_vtas_along_traj)
     DEALLOCATE(contrail_opt_vtas)          !20150313-2
     DEALLOCATE(contrail_opt_vground)       !20150313-2
     DEALLOCATE(d_zalt)                     !20150313-2

     !Yin_20170423

     !USAGE CONFIRMATION
!     write(*,*)"calculate_emissions subroutine is done"

  END SUBROUTINE calculate_emissions_along_trajectory 
  ! =========================================================================

  ! =========================================================================
!  SUBROUTINE update_trajectory 

!  IMPLICIT NONE 


!  END SUBROUTINE update_trajectory
  ! =========================================================================

  ! =========================================================================
  SUBROUTINE fly_aircraft(d_now, d_ac_routes, d_p2_ac_routes_desc) 

    ! ------------------------------------------------------------------
    ! This routine calculates followings: 
    ! 
    !     - 1. The subroutine makes aircrafts fly along flight plan determined 
    !       by ac_routes(lon,lat,alt,time). 
    !     - 2. New(AC_POS) and old(AC_POS_OLD) aircraft positions are caulcated,
    !       and expressed by waypoint numbers, such as
    !             D_city=1.0, AC_POS_OLD=1.5, AC_POS=3.2, A_city=nwaypoints
    !
    !     - 3. The subroutine executes 'Arrival check.'
    !
    ! ------------------------------------------------------------------
    ! TIPS:
    !     1. We assume that all aircraft fly along the determined flight plan (ac_routes).
    !      
    !     2. d_now, and d_p2_ac_routes_desc(D_time) are defined in Julian_date, provided by SMIL. 
    !
    !     3. Departure check is performed in airtraf SMIL. When aircrafts come to this subroutine 
    !       for its first time by departure check, the aircraft is at D_city:
    !            AC_POS_OLD = 1.0
    !            AC_POS     = 1.0 
    !    
    !     4. Arrival check is performed in this subroutine. When aircrafts arrive at A_city,
    !       AC_POS becomes the value of nwaypoints, defined by CTRL-NAMELIST.
    !  
    !     5. Tips 3 and 4 treatments can calculate emissions at both D- and A-cities.
    !        See note 20120903-1. This treatment is necessary for very short flight.  

     IMPLICIT NONE 
!    INTEGER, PARAMETER, private              :: AC_POS=10 , AC_POS_OUT=11 

     REAL(DP), DIMENSION(:,:), INTENT(INOUT)  :: d_ac_routes
!    INTEGER, DIMENSION(:), INTENT(INOUT)     :: d_p2_ac_routes_desc
     REAL(DP), DIMENSION(:), INTENT(INOUT)    :: d_p2_ac_routes_desc
     REAL(DP), INTENT(IN)                     :: d_now
     INTEGER                                  :: i, upper_waypoint, lower_waypoint
     REAL(DP)                                 :: ratio_btw_waypoints

! Initilizations
     upper_waypoint      = 0
     lower_waypoint      = 0
     ratio_btw_waypoints = 0.0_dp

! 1. Calculate new AC_pos and AC_pos_old
! Note: 20120404-2,20120405-1,20120724-2,20120719-5
     d_p2_ac_routes_desc(AC_POS_OLD)=d_p2_ac_routes_desc(AC_POS) !AC_pos -> AC_pos_old

     ! For first step
     IF(d_now <= d_ac_routes(1,props_time))THEN  
        d_p2_ac_routes_desc(AC_POS)     = 1.0_dp 
        d_p2_ac_routes_desc(AC_POS_OLD) = 1.0_dp 
!        write(*,*)'check_no_case1',d_p2_ac_routes_desc(AC_POS),d_p2_ac_routes_desc(AC_POS_OLD)
!        write(*,*)d_now,'<=',d_ac_routes(1,props_time) 
     ENDIF


     ! Inflight step
     ! Values of current lon and lat can be also calculated in the same way, if necessary.
     ! Program on Arrival check! see note 20120719-5,-7
     ! When aircraft reaches A_city, p2_ac_routes_desc(i_traj,AC_pos) should be changed
     ! into the value of nwaypoints (D_city point array number)

!HY 20120903-2
!    IF(d_now > d_ac_routes(1,props_time))THEN  
     IF(d_now > d_ac_routes(1,props_time) .and. d_now < d_ac_routes(nwaypoints,props_time))THEN  
        DO i=1, nwaypoints-1
           IF(d_ac_routes(i,props_time)<d_now .and. d_ac_routes(i+1,props_time)>=d_now)THEN 
              lower_waypoint = i
              upper_waypoint = i+1 
              ratio_btw_waypoints = (d_now-d_ac_routes(lower_waypoint,props_time))/   &
                              (d_ac_routes(upper_waypoint,props_time)-d_ac_routes(lower_waypoint,props_time))
!              write(*,*)'check_btw1',(d_now-d_ac_routes(lower_waypoint,props_time))
!              write(*,*)'check_btw2',(d_ac_routes(upper_waypoint,props_time)-d_ac_routes(lower_waypoint,props_time))
              d_p2_ac_routes_desc(AC_POS)=ratio_btw_waypoints*(i+1) + (1.0_dp-ratio_btw_waypoints)*i
              exit  ! leave do loop
           ELSE
!              write(*,*)'check_no_case2 ','props_time=',props_time 
!              write(*,*)d_ac_routes(i,props_time),'<',d_now,'<=',d_ac_routes(i+1,props_time) 
           ENDIF
        ENDDO
     ELSE IF(d_now >= d_ac_routes(nwaypoints,props_time))THEN
        d_p2_ac_routes_desc(AC_POS)= REAL(nwaypoints,DP)   ! Aircrafts have arrived now!
!        write(*,*)'arival_check is done'
     ELSE 
!        write(*,*)'Aircaft just departured!'
     ENDIF

!     write(*,*)'fly_aircraft_check'
!     write(*,*)'d_now',d_now
!     write(*,*)'lower_waypoint',lower_waypoint,'upper_waypoint',upper_waypoint
!     write(*,*)'ratio_btw_waypoints',ratio_btw_waypoints
!     write(*,*)'AC_POS_OLD=',d_p2_ac_routes_desc(AC_POS_OLD),'AC_POS=',d_p2_ac_routes_desc(AC_POS)
!     write(*,*)'REAL=',REAL(nwaypoints,DP),'DBLE=',DBLE(nwaypoints)
    
     !DUMMY
!     d_ac_routes(:,:)       = 1.0_dp
!     d_p2_ac_routes_desc(:) = 1.0_dp
!HY   d_p2_ac_routes_desc(AC_POS) = 1.0_dp

     !USAGE CONFIRMATION
!     write(*,*)"fly_aircraft subroutine is done"

  END SUBROUTINE fly_aircraft
  ! =========================================================================

  ! =========================================================================
!  SUBROUTINE arrival 

!  IMPLICIT NONE 



!  END SUBROUTINE arrival 
  ! =========================================================================

  ! =========================================================================
  SUBROUTINE compress_traj(d_ac_routes, d_ac_routes_out, Julian_date_start) 

    ! ------------------------------------------------------------------
    ! This subroutine calculates followings: 
    ! 
    !     - The idea is that the output is too large for a long, like 10 year simulation.
    !       So 'compress' means that chunks are calculated an longer legs stored.
    !       For example if we have 100 waypoints per trajectory, we calculate the Sum over 
    !       20 gridpoints for emissions and take middle value for location, time, and other
    !       properties. We get a much smaller output. Whenever an aircraft lands, the 
    !       trajectory is taken out of the active array (ac_routes), compressed and stored
    !       in the compressed output(ac_routes_out). That means that the order changes,
    !       the active array(ac_routes) gets new trajs in the order of the departures;
    !       while the compress(ac_routes_out) gets trajs in the order of arrivals.
    !        
    !     - ac_routes should be compressed and transferred to ac_routes_out: 
    !       ac_routes -> ac_routes_out
    !       See note 20120906-1 
    !  
    !     - NOTE: Now, data are NOT compressed in 'nlroutes' dimension.
    !             Current code can compress data only in 'nwaypoints' dimension.
    !             See note 20120906-1, 20120911-2
    !
    ! ------------------------------------------------------------------
    ! TIPS:
    !     1. When this subroutine is applied to an aircraft, the aircraft has just got 
    !        arrival check. The order of arrival is expressed by i_traj_out.
    !       
    !     2. The extent of compression is defined by ngroutes_out and nwaypoints_out in airtraf.nml. 
    !        ngroutes_out is decomposed in SMCL, and nlroutes_out is calculated on each PE.
    !       
    !     3. Active array: ac_routes(nwaypoints, nprops, nlroutes,1)
    !        Output array: ac_routes_out(nwaypoints_out, nprops, nlroutes_out,1)
    !        NOTE: The properties of props_time in ac_routes is Julian_date.
    !              The properties of props_time in ac_routes_out is (Julian_date - Julian_date_start) 
    !
    !     4. nwaypoints should be even number, which can be devided by nwaypoints_out.
    !        I.e. nwaypoints_out is one of a measure of nwaypoints.
    !        E.g.) when nwaypoints = 12,
    !              nwaypoints_out can be selected among (1,2,3,4,6,12).
    !        

     IMPLICIT NONE 
     REAL(DP), DIMENSION(:,:), INTENT(IN)   :: d_ac_routes
     REAL(DP), DIMENSION(:,:), INTENT(OUT)  :: d_ac_routes_out
     REAL(DP),                 INTENT(IN)   :: Julian_date_start
     INTEGER                                :: i, j, compressed_waypoints
     
     ! With compressed 
     IF(nwaypoints_out < nwaypoints)THEN 
        compressed_waypoints = nwaypoints/nwaypoints_out
!        write(*,*)'compressed_waypoints',compressed_waypoints 
!        write(*,*)'calculate_compress_traj'
!        write(*,*)'nwaypoints=',nwaypoints,'nwaypoints_out=',nwaypoints_out
!        write(*,*)'ngroutes=',ngroutes,'ngroutes_out=',ngroutes_out
        DO j=1,nprops   !nprops=8
           DO i=1,nwaypoints_out
              d_ac_routes_out(i,j)=sum(d_ac_routes((i*compressed_waypoints-(compressed_waypoints-1)): &
                                                   (i*compressed_waypoints),j))/compressed_waypoints
              IF(j.eq.props_time)d_ac_routes_out(i,j)=d_ac_routes_out(i,j)-Julian_date_start

!              write(*,*)'nprops=',j,'nwaypoints_out=',i,'compressed_waypoints',compressed_waypoints
!              write(*,*)'d_ac_routes=',d_ac_routes((i*compressed_waypoints-(compressed_waypoints-1)): &
!                                                   (i*compressed_waypoints),j)
!              write(*,*)'sum(d_ac_routes)=',sum(d_ac_routes((i*compressed_waypoints-(compressed_waypoints-1)): &
!                                                            (i*compressed_waypoints),j))
!              write(*,*)'average_sum(d_ac_routes)=',sum(d_ac_routes((i*compressed_waypoints-(compressed_waypoints-1)): &
!                                                                    (i*compressed_waypoints),j))/compressed_waypoints
!              write(*,*)'j_1=',j,'d_ac_routes_out',d_ac_routes_out(i,j)
           ENDDO
        ENDDO
        
     ! Without compressed
     ELSEIF(nwaypoints_out==nwaypoints)THEN
        DO j=1,nprops
           DO i=1,nwaypoints_out
              d_ac_routes_out(i,j) = d_ac_routes(i,j)
              IF(j.eq.props_time)d_ac_routes_out(i,j)=d_ac_routes_out(i,j)-Julian_date_start
           ENDDO
!           write(*,*)'calculate_compress_traj(nwaypoints_out=nwaypoints)'
!           write(*,*)'d_ac_routes_out(lon)=',d_ac_routes_out(:,1)
!           write(*,*)'d_ac_routes(lon)=',d_ac_routes(:,1)
!           write(*,*)'j_2=',j,'d_ac_routes_out(time)=',d_ac_routes_out(:,4)
!           write(*,*)'j_2=',j,'d_ac_routes(time)=',d_ac_routes(:,4)
        ENDDO

     !DUMMY
     ELSE  
        d_ac_routes_out(:,:) = 1.0_dp
!        write(*,*)'error:calculate_compress_traj'
!        write(*,*)'nwaypoints_out should be: nwaypoints_out < nwaypoints'
     ENDIF

     !USAGE CONFIRMATION
!     write(*,*)"compress_traj subroutine is done"

  END SUBROUTINE compress_traj
  ! =========================================================================

  ! =========================================================================
!  SUBROUTINE append_traj 

!  IMPLICIT NONE 



!  END SUBROUTINE append_traj
  ! =========================================================================

  ! =========================================================================
  SUBROUTINE fill_emission_field(emistype, d_ac_routes, d_p2_ac_routes_desc_new, &
                                 d_p2_ac_routes_desc_old, d_add_glob_emis, d_i_add_emis) 

    ! ------------------------------------------------------------------
    ! This subroutine caluclates emissions for 1 time step by using 
    ! aircraft position information. 
    ! 
    !     - The d_add_glob_emis gathers info(lon,lat,alt,emiss) on a trajectory. 
    !       for 1 time step. Emissions will be calculated by using previous(AC_POS_OLD),
    !       current(AC_POS) and ac_routes properties. The emissions are added into 
    !       glob_(NOx)_emis in SMIL after this routine.
    !       
    !       See 20121114-4,20121115-2,20120911-6,20120801-2  
    !  
    !     - d_i_add_emis shows a number of elements of d_add_glob_emis.
    !
    ! ------------------------------------------------------------------
    ! TIPS:
    !     1. d_ac_rotes has aircraft properties.  
    !       
    !     2. d_p2_ac_routes_desc_new and d_p2_ac_routes_desc_old have aircraft position
    !        informations expressed by waypoints number(e.g., AC_POS_OLD =1.5, AC_POS=3.4). 
    !       
    !     3. When aircraft is assiged as departure time,
    !             AC_POS_OLD = 1.0
    !             AC_POS     = 1.0
    !        in this case, NO emission is estimated:
    !             d_add_glob_emis = 0.0
    !             d_i_add_emis    = 0.0
    !             

     IMPLICIT NONE
     REAL(DP), DIMENSION(:,:), INTENT(IN)                :: d_ac_routes
!    INTEGER, INTENT(INOUT)                              :: d_p2_ac_routes_desc_new
!    INTEGER, INTENT(INOUT)                              :: d_p2_ac_routes_desc_old
     REAL(DP), INTENT(INOUT)                             :: d_p2_ac_routes_desc_new
     REAL(DP), INTENT(INOUT)                             :: d_p2_ac_routes_desc_old
!    REAL(DP), DIMENSION(:,:,:), INTENT(INOUT)           :: d_glob_NOx_emis 
     REAL(DP), DIMENSION(:,:), INTENT(INOUT)             :: d_add_glob_emis
!    REAL(DP), DIMENSION(:,:,:), INTENT(INOUT)           :: d_glob_emis 
     INTEGER, INTENT(IN)                                 :: emistype
     INTEGER, INTENT(OUT)                                :: d_i_add_emis
     INTEGER                                             :: i


     !INITIALIZATION
     d_add_glob_emis(:,:) = 0.0_dp
     d_i_add_emis         = 0

     !write(*,*)'check_initialization on d_add_glob_emis =',d_add_glob_emis(:,:) 
!     write(*,*)"Emistype = ", emistype
!     write(*,*)"AC_POS= ", d_p2_ac_routes_desc_new, "AC_POS_OLD= ",d_p2_ac_routes_desc_old
!     write(*,*)'props_lon=',props_lon,'props_lat=',props_lat,'props_alt=',props_alt
!     write(*,*)'props_emis_nox=',props_emis_nox,'props_emis_h2o=',props_emis_h2o
!     write(*,*)'props_dist=',props_dist
     !Yin_20170423
!     write(*,*)'props_potcov=',props_potcov
     !Yin_20170423
     SELECT CASE(emistype)
!------------------------------------------
     CASE(NOx)
!------------------------------------------
!-------AIRCRAFT IS ASSIGNED DEPARTURE TIME
        IF((d_p2_ac_routes_desc_old == 1.0_dp).and. &
           (d_p2_ac_routes_desc_new == 1.0_dp))THEN
!           write(*,*)'NOx_Aircraft gets departure time.' 
           !Nothing: d_i_add_emis=0 

!-------AIRCRAFT LEAVES FROM D_CITY
        ELSEIF((d_p2_ac_routes_desc_old == 1.0_dp)                .and. &
               (d_p2_ac_routes_desc_old < d_p2_ac_routes_desc_new).and. &
               (d_p2_ac_routes_desc_new < REAL(nwaypoints,DP)))THEN 
!               write(*,*)'NOx_Aircraft leaves from D_city.' 
               DO i=1,nwaypoints
                  IF(d_p2_ac_routes_desc_old + REAL(i,DP) < d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_emis_nox)
                  ELSEIF(d_p2_ac_routes_desc_old + REAL(i,DP) >= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_emis_nox)*&
                                                            (d_p2_ac_routes_desc_new - REAL(i,DP))
                     exit  !
                  ELSE
                     write(*,*)'NOx_ERROR: Aircraft leaves from D_city.'
                  ENDIF
               ENDDO
!               write(*,*)'NOx_d_i_add_emis= ',d_i_add_emis

!-------AIRCRAFT IS INFLIGHT
        ELSEIF((d_p2_ac_routes_desc_old > 1.0_dp)                 .and. &
               (d_p2_ac_routes_desc_old < d_p2_ac_routes_desc_new).and. &
               (d_p2_ac_routes_desc_new < REAL(nwaypoints,DP)))THEN 
!               write(*,*)'NOx_Aircraft is inflight.' 
            IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) < d_p2_ac_routes_desc_new)THEN
               IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) /= d_p2_ac_routes_desc_old)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_emis_nox)*&
                                                            (DBLE(CEILING(d_p2_ac_routes_desc_old)) - d_p2_ac_routes_desc_old)
               ENDIF
               DO i=1,nwaypoints
                  IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) < d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_emis_nox)
                  ELSEIF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) >= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_emis_nox)*&
                                                            (d_p2_ac_routes_desc_new - &
                                                             DBLE(CEILING(d_p2_ac_routes_desc_old)) - REAL(i-1,DP))
                     exit  !
                  ELSE
                     write(*,*)'NOx_ERROR_1: Aircraft inflight.'
                  ENDIF
               ENDDO
            ELSEIF(DBLE(CEILING(d_p2_ac_routes_desc_old)) >= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_emis_nox)*&
                                                            (d_p2_ac_routes_desc_new - d_p2_ac_routes_desc_old)
            ELSE
                     write(*,*)'NOx_ERROR_2: Aircraft inflight.'
            ENDIF 
!            write(*,*)'NOx_d_i_add_emis= ',d_i_add_emis

!-------AIRCRAFT ARRIVES AT A_CITY
        ELSEIF((d_p2_ac_routes_desc_old > 1.0_dp)                 .and. &
               (d_p2_ac_routes_desc_old < d_p2_ac_routes_desc_new).and. &
               (d_p2_ac_routes_desc_new == REAL(nwaypoints,DP)))THEN 
!               write(*,*)'NOx_Aircraft arrives at A_city.' 
            IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) <= d_p2_ac_routes_desc_new)THEN
               IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) /= d_p2_ac_routes_desc_old)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_emis_nox)*&
                                                            (DBLE(CEILING(d_p2_ac_routes_desc_old)) - d_p2_ac_routes_desc_old)
               ENDIF
               DO i=1,nwaypoints
                  IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) <= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_emis_nox)
                  ENDIF
                  !See 20121115-2, case(III)
                  IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) >= d_p2_ac_routes_desc_new)exit  !
               ENDDO
            ELSE
               write(*,*)'NOx_ERROR: Aircraft arrives at A_city.'
            ENDIF 
!            write(*,*)'NOx_d_i_add_emis= ',d_i_add_emis

!-------AIRCRAFT FLY FROM D_CITY TO A_CITY WITHIN 1 TIME STEP
        ELSEIF((d_p2_ac_routes_desc_old == 1.0_dp).and. &
               (d_p2_ac_routes_desc_new == REAL(nwaypoints,DP)))THEN
!               write(*,*)'NOx_Aircraft fly from D_ to A_city within 1 time step.' 
               DO i=1,nwaypoints-1
                  d_i_add_emis = d_i_add_emis + 1
                  d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lon)
                  d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lat)
                  d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_alt)
                  d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_emis_nox)
               ENDDO
!               write(*,*)'NOx_d_i_add_emis= ',d_i_add_emis

!-------ERROR MESSAGE
        ELSE
           write(*,*)'NOx_ERROR: Fill_emission_field.'
        ENDIF

        !DUMMY
!       d_glob_emis(:,:,:) = 1.0_dp + d_glob_emis(:,:,:) 
!        write(*,*)"NOx_fill_emission_field."

!------------------------------------------
     CASE(H2O)
!------------------------------------------
!-------AIRCRAFT IS ASSIGNED DEPARTURE TIME
        IF((d_p2_ac_routes_desc_old == 1.0_dp).and. &
           (d_p2_ac_routes_desc_new == 1.0_dp))THEN
!           write(*,*)'H2O_Aircraft gets departure time.' 
           !Nothing: d_i_add_emis=0 

!-------AIRCRAFT LEAVES FROM D_CITY
        ELSEIF((d_p2_ac_routes_desc_old == 1.0_dp)                .and. &
               (d_p2_ac_routes_desc_old < d_p2_ac_routes_desc_new).and. &
               (d_p2_ac_routes_desc_new < REAL(nwaypoints,DP)))THEN 
!               write(*,*)'H2O_Aircraft leaves from D_city.' 
               DO i=1,nwaypoints
                  IF(d_p2_ac_routes_desc_old + REAL(i,DP) < d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_emis_h2o)
                  ELSEIF(d_p2_ac_routes_desc_old + REAL(i,DP) >= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_emis_h2o)*&
                                                            (d_p2_ac_routes_desc_new - REAL(i,DP))
                     exit  !
                  ELSE
                     write(*,*)'H2O_ERROR: Aircraft leaves from D_city.'
                  ENDIF
               ENDDO
!               write(*,*)'H2O_d_i_add_emis= ',d_i_add_emis

!-------AIRCRAFT IS INFLIGHT
        ELSEIF((d_p2_ac_routes_desc_old > 1.0_dp)                 .and. &
               (d_p2_ac_routes_desc_old < d_p2_ac_routes_desc_new).and. &
               (d_p2_ac_routes_desc_new < REAL(nwaypoints,DP)))THEN 
!               write(*,*)'H2O_Aircraft is inflight.' 
            IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) < d_p2_ac_routes_desc_new)THEN
               IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) /= d_p2_ac_routes_desc_old)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_emis_h2o)*&
                                                            (DBLE(CEILING(d_p2_ac_routes_desc_old)) - d_p2_ac_routes_desc_old)
               ENDIF
               DO i=1,nwaypoints
                  IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) < d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_emis_h2o)
                  ELSEIF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) >= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_emis_h2o)*&
                                                            (d_p2_ac_routes_desc_new - &
                                                             DBLE(CEILING(d_p2_ac_routes_desc_old)) - REAL(i-1,DP))
                     exit  !
                  ELSE
                     write(*,*)'H2O_ERROR_1: Aircraft inflight.'
                  ENDIF
               ENDDO
            ELSEIF(DBLE(CEILING(d_p2_ac_routes_desc_old)) >= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_emis_h2o)*&
                                                            (d_p2_ac_routes_desc_new - d_p2_ac_routes_desc_old)
            ELSE
                     write(*,*)'H2O_ERROR_2: Aircraft inflight.'
            ENDIF 
!            write(*,*)'H2O_d_i_add_emis= ',d_i_add_emis

!-------AIRCRAFT ARRIVES AT A_CITY
        ELSEIF((d_p2_ac_routes_desc_old > 1.0_dp)                 .and. &
               (d_p2_ac_routes_desc_old < d_p2_ac_routes_desc_new).and. &
               (d_p2_ac_routes_desc_new == REAL(nwaypoints,DP)))THEN 
!               write(*,*)'H2O_Aircraft arrives at A_city.' 
            IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) <= d_p2_ac_routes_desc_new)THEN
               IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) /= d_p2_ac_routes_desc_old)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_emis_h2o)*&
                                                            (DBLE(CEILING(d_p2_ac_routes_desc_old)) - d_p2_ac_routes_desc_old)
               ENDIF
               DO i=1,nwaypoints
                  IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) <= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_emis_h2o)
                  ENDIF
                  !See 20121115-2, case(III)
                  IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) >= d_p2_ac_routes_desc_new)exit  !
               ENDDO
            ELSE
               write(*,*)'H2O_ERROR: Aircraft arrives at A_city.'
            ENDIF 
!            write(*,*)'H2O_d_i_add_emis= ',d_i_add_emis

!-------AIRCRAFT FLY FROM D_CITY TO A_CITY WITHIN 1 TIME STEP
        ELSEIF((d_p2_ac_routes_desc_old == 1.0_dp).and. &
               (d_p2_ac_routes_desc_new == REAL(nwaypoints,DP)))THEN
!               write(*,*)'H2O_Aircraft fly from D_ to A_city within 1 time step.' 
               DO i=1,nwaypoints-1
                  d_i_add_emis = d_i_add_emis + 1
                  d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lon)
                  d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lat)
                  d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_alt)
                  d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_emis_h2o)
               ENDDO
!               write(*,*)'H2O_d_i_add_emis= ',d_i_add_emis

!-------ERROR MESSAGE
        ELSE
           write(*,*)'H2O_ERROR: Fill_emission_field.'
        ENDIF

        !DUMMY
!       d_add_glob_emis(:,:) = 1.0_dp + d_add_glob_emis(:,:)
!       d_glob_emis(:,:,:) = 1.0_dp + d_glob_emis(:,:,:)
!        write(*,*)"H2O_fill_emission_field."

!-----------------------------------------
     CASE(DIST)
!-----------------------------------------
!-------AIRCRAFT IS ASSIGNED DEPARTURE TIME
        IF((d_p2_ac_routes_desc_old == 1.0_dp).and. &
           (d_p2_ac_routes_desc_new == 1.0_dp))THEN
!           write(*,*)'DIST_Aircraft gets departure time.' 
           !Nothing: d_i_add_emis=0 

!-------AIRCRAFT LEAVES FROM D_CITY
        ELSEIF((d_p2_ac_routes_desc_old == 1.0_dp)                .and. &
               (d_p2_ac_routes_desc_old < d_p2_ac_routes_desc_new).and. &
               (d_p2_ac_routes_desc_new < REAL(nwaypoints,DP)))THEN 
!               write(*,*)'DIST_Aircraft leaves from D_city.' 
               DO i=1,nwaypoints
                  IF(d_p2_ac_routes_desc_old + REAL(i,DP) < d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_dist)
                  ELSEIF(d_p2_ac_routes_desc_old + REAL(i,DP) >= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_dist)*&
                                                            (d_p2_ac_routes_desc_new - REAL(i,DP))
                     exit  !
                  ELSE
                     write(*,*)'DIST_ERROR: Aircraft leaves from D_city.'
                  ENDIF
               ENDDO
!               write(*,*)'DIST_d_i_add_emis= ',d_i_add_emis

!-------AIRCRAFT IS INFLIGHT
        ELSEIF((d_p2_ac_routes_desc_old > 1.0_dp)                 .and. &
               (d_p2_ac_routes_desc_old < d_p2_ac_routes_desc_new).and. &
               (d_p2_ac_routes_desc_new < REAL(nwaypoints,DP)))THEN 
!               write(*,*)'DIST_Aircraft is inflight.' 
            IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) < d_p2_ac_routes_desc_new)THEN
               IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) /= d_p2_ac_routes_desc_old)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_dist)*&
                                                            (DBLE(CEILING(d_p2_ac_routes_desc_old)) - d_p2_ac_routes_desc_old)
               ENDIF
               DO i=1,nwaypoints
                  IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) < d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_dist)
                  ELSEIF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) >= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_dist)*&
                                                            (d_p2_ac_routes_desc_new - &
                                                             DBLE(CEILING(d_p2_ac_routes_desc_old)) - REAL(i-1,DP))
                     exit  !
                  ELSE
                     write(*,*)'DIST_ERROR_1: Aircraft inflight.'
                  ENDIF
               ENDDO
            ELSEIF(DBLE(CEILING(d_p2_ac_routes_desc_old)) >= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_dist)*&
                                                            (d_p2_ac_routes_desc_new - d_p2_ac_routes_desc_old)
            ELSE
                     write(*,*)'DIST_ERROR_2: Aircraft inflight.'
            ENDIF 
!            write(*,*)'DIST_d_i_add_emis= ',d_i_add_emis

!-------AIRCRAFT ARRIVES AT A_CITY
        ELSEIF((d_p2_ac_routes_desc_old > 1.0_dp)                 .and. &
               (d_p2_ac_routes_desc_old < d_p2_ac_routes_desc_new).and. &
               (d_p2_ac_routes_desc_new == REAL(nwaypoints,DP)))THEN 
!               write(*,*)'DIST_Aircraft arrives at A_city.' 
            IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) <= d_p2_ac_routes_desc_new)THEN
               IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) /= d_p2_ac_routes_desc_old)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_dist)*&
                                                            (DBLE(CEILING(d_p2_ac_routes_desc_old)) - d_p2_ac_routes_desc_old)
               ENDIF
               DO i=1,nwaypoints
                  IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) <= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_dist)
                  ENDIF
                  !See 20121115-2, case(III)
                  IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) >= d_p2_ac_routes_desc_new)exit  !
               ENDDO
            ELSE
               write(*,*)'DIST_ERROR: Aircraft arrives at A_city.'
            ENDIF 
!            write(*,*)'DIST_d_i_add_emis= ',d_i_add_emis

!-------AIRCRAFT FLY FROM D_CITY TO A_CITY WITHIN 1 TIME STEP
        ELSEIF((d_p2_ac_routes_desc_old == 1.0_dp).and. &
               (d_p2_ac_routes_desc_new == REAL(nwaypoints,DP)))THEN
!               write(*,*)'DIST_Aircraft fly from D_ to A_city within 1 time step.' 
               DO i=1,nwaypoints-1
                  d_i_add_emis = d_i_add_emis + 1
                  d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lon)
                  d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lat)
                  d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_alt)
                  d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_dist)
               ENDDO
!               write(*,*)'DIST_d_i_add_emis= ',d_i_add_emis

!-------ERROR MESSAGE
        ELSE
           write(*,*)'DIST_ERROR: Fill_emission_field.'
        ENDIF

        !DUMMY
!       d_add_glob_emis(:,:) = 1.0_dp + d_add_glob_emis(:,:)
!       d_glob_emis(:,:,:) = 1.0_dp + d_glob_emis(:,:,:)
!        write(*,*)"DIST_fill_emission_field."

!-----------------------------------------
     CASE(FUEL)
!-----------------------------------------
!-------AIRCRAFT IS ASSIGNED DEPARTURE TIME
        IF((d_p2_ac_routes_desc_old == 1.0_dp).and. &
           (d_p2_ac_routes_desc_new == 1.0_dp))THEN
!           write(*,*)'FUEL_Aircraft gets departure time.' 
           !Nothing: d_i_add_emis=0 

!-------AIRCRAFT LEAVES FROM D_CITY
        ELSEIF((d_p2_ac_routes_desc_old == 1.0_dp)                .and. &
               (d_p2_ac_routes_desc_old < d_p2_ac_routes_desc_new).and. &
               (d_p2_ac_routes_desc_new < REAL(nwaypoints,DP)))THEN 
!               write(*,*)'FUEL_Aircraft leaves from D_city.' 
               DO i=1,nwaypoints
                  IF(d_p2_ac_routes_desc_old + REAL(i,DP) < d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_fuel_use)
!                    d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_dist)*&
!                                                            FUEL_CONSUMPTION_RATE
                  ELSEIF(d_p2_ac_routes_desc_old + REAL(i,DP) >= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_fuel_use)*&
                                                            (d_p2_ac_routes_desc_new - REAL(i,DP))
!                    d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_dist)*&
!                                                           (d_p2_ac_routes_desc_new - REAL(i,DP))*FUEL_CONSUMPTION_RATE
                     exit  !
                  ELSE
                     write(*,*)'FUEL_ERROR: Aircraft leaves from D_city.'
                  ENDIF
               ENDDO
!               write(*,*)'FUEL_d_i_add_emis= ',d_i_add_emis

!-------AIRCRAFT IS INFLIGHT
        ELSEIF((d_p2_ac_routes_desc_old > 1.0_dp)                 .and. &
               (d_p2_ac_routes_desc_old < d_p2_ac_routes_desc_new).and. &
               (d_p2_ac_routes_desc_new < REAL(nwaypoints,DP)))THEN 
!               write(*,*)'FUEL_Aircraft is inflight.' 
            IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) < d_p2_ac_routes_desc_new)THEN
               IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) /= d_p2_ac_routes_desc_old)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_fuel_use)*&
                                                            (DBLE(CEILING(d_p2_ac_routes_desc_old)) - d_p2_ac_routes_desc_old)
!                    d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_dist)*&
!                                                           (DBLE(CEILING(d_p2_ac_routes_desc_old)) - d_p2_ac_routes_desc_old)*&
!                                                            FUEL_CONSUMPTION_RATE
               ENDIF
               DO i=1,nwaypoints
                  IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) < d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_fuel_use)
!                    d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_dist)*&
!                                                            FUEL_CONSUMPTION_RATE
                  ELSEIF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) >= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_fuel_use)*&
                                                            (d_p2_ac_routes_desc_new - &
                                                             DBLE(CEILING(d_p2_ac_routes_desc_old)) - REAL(i-1,DP))
!                    d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_dist)*&
!                                                           (d_p2_ac_routes_desc_new - &
!                                                            DBLE(CEILING(d_p2_ac_routes_desc_old)) - REAL(i-1,DP))*&
!                                                            FUEL_CONSUMPTION_RATE
                     exit  !
                  ELSE
                     write(*,*)'FUEL_ERROR_1: Aircraft inflight.'
                  ENDIF
               ENDDO
            ELSEIF(DBLE(CEILING(d_p2_ac_routes_desc_old)) >= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_fuel_use)*&
                                                            (d_p2_ac_routes_desc_new - d_p2_ac_routes_desc_old)
!                    d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_dist)*&
!                                                           (d_p2_ac_routes_desc_new - d_p2_ac_routes_desc_old)*&
!                                                            FUEL_CONSUMPTION_RATE
            ELSE
                     write(*,*)'FUEL_ERROR_2: Aircraft inflight.'
            ENDIF 
!            write(*,*)'FUEL_d_i_add_emis= ',d_i_add_emis

!-------AIRCRAFT ARRIVES AT A_CITY
        ELSEIF((d_p2_ac_routes_desc_old > 1.0_dp)                 .and. &
               (d_p2_ac_routes_desc_old < d_p2_ac_routes_desc_new).and. &
               (d_p2_ac_routes_desc_new == REAL(nwaypoints,DP)))THEN 
!               write(*,*)'FUEL_Aircraft arrives at A_city.' 
            IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) <= d_p2_ac_routes_desc_new)THEN
               IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) /= d_p2_ac_routes_desc_old)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_fuel_use)*&
                                                            (DBLE(CEILING(d_p2_ac_routes_desc_old)) - d_p2_ac_routes_desc_old)
!                    d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_dist)*&
!                                                           (DBLE(CEILING(d_p2_ac_routes_desc_old)) - d_p2_ac_routes_desc_old)*&
!                                                            FUEL_CONSUMPTION_RATE
               ENDIF
               DO i=1,nwaypoints
                  IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) <= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_fuel_use)
!                    d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_dist)*&
!                                                            FUEL_CONSUMPTION_RATE
                  ENDIF
                  !See 20121115-2, case(III)
                  IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) >= d_p2_ac_routes_desc_new)exit  !
               ENDDO
            ELSE
               write(*,*)'FUEL_ERROR: Aircraft arrives at A_city.'
            ENDIF 
!            write(*,*)'FUEL_d_i_add_emis= ',d_i_add_emis

!-------AIRCRAFT FLY FROM D_CITY TO A_CITY WITHIN 1 TIME STEP
        ELSEIF((d_p2_ac_routes_desc_old == 1.0_dp).and. &
               (d_p2_ac_routes_desc_new == REAL(nwaypoints,DP)))THEN
!               write(*,*)'FUEL_Aircraft fly from D_ to A_city within 1 time step.' 
               DO i=1,nwaypoints-1
                  d_i_add_emis = d_i_add_emis + 1
                  d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lon)
                  d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lat)
                  d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_alt)
                  d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_fuel_use)
!                 d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_dist)*&
!                                                         FUEL_CONSUMPTION_RATE
               ENDDO
!               write(*,*)'FUEL_d_i_add_emis= ',d_i_add_emis

!-------ERROR MESSAGE
        ELSE
           write(*,*)'FUEL_ERROR: Fill_emission_field.'
        ENDIF

        !DUMMY
!       d_add_glob_emis(:,:) = 1.0_dp + d_add_glob_emis(:,:)
!       d_glob_emis(:,:,:) = 1.0_dp + d_glob_emis(:,:,:)
!        write(*,*)"FUEL_fill_emission_field."

!Yin_20170423
!------------------------------------------
     CASE(C_PC)
!------------------------------------------
!-------AIRCRAFT IS ASSIGNED DEPARTURE TIME
        IF((d_p2_ac_routes_desc_old == 1.0_dp).and. &
           (d_p2_ac_routes_desc_new == 1.0_dp))THEN
!           write(*,*)'CPC_Aircraft gets departure time.' 
           !Nothing: d_i_add_emis=0 

!-------AIRCRAFT LEAVES FROM D_CITY
        ELSEIF((d_p2_ac_routes_desc_old == 1.0_dp)                .and. &
               (d_p2_ac_routes_desc_old < d_p2_ac_routes_desc_new).and. &
               (d_p2_ac_routes_desc_new < REAL(nwaypoints,DP)))THEN 
!               write(*,*)'CPC_Aircraft leaves from D_city.' 
               DO i=1,nwaypoints
                  IF(d_p2_ac_routes_desc_old + REAL(i,DP) < d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_potcov)
                  ELSEIF(d_p2_ac_routes_desc_old + REAL(i,DP) >= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_potcov)*&
                                                            (d_p2_ac_routes_desc_new - REAL(i,DP))
                     exit  !
                  ELSE
                     write(*,*)'CPC_ERROR: Aircraft leaves from D_city.'
                  ENDIF
               ENDDO
!               write(*,*)'CPC_d_i_add_emis= ',d_i_add_emis

!-------AIRCRAFT IS INFLIGHT
        ELSEIF((d_p2_ac_routes_desc_old > 1.0_dp)                 .and. &
               (d_p2_ac_routes_desc_old < d_p2_ac_routes_desc_new).and. &
               (d_p2_ac_routes_desc_new < REAL(nwaypoints,DP)))THEN 
!               write(*,*)'CPC_Aircraft is inflight.' 
            IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) < d_p2_ac_routes_desc_new)THEN
               IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) /= d_p2_ac_routes_desc_old)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_potcov)*&
                                                            (DBLE(CEILING(d_p2_ac_routes_desc_old)) - d_p2_ac_routes_desc_old)
               ENDIF
               DO i=1,nwaypoints
                  IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) < d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_potcov)
                  ELSEIF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) >= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_potcov)*&
                                                            (d_p2_ac_routes_desc_new - &
                                                             DBLE(CEILING(d_p2_ac_routes_desc_old)) - REAL(i-1,DP))
                     exit  !
                  ELSE
                     write(*,*)'CPC_ERROR_1: Aircraft inflight.'
                  ENDIF
               ENDDO
            ELSEIF(DBLE(CEILING(d_p2_ac_routes_desc_old)) >= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_potcov)*&
                                                            (d_p2_ac_routes_desc_new - d_p2_ac_routes_desc_old)
            ELSE
                     write(*,*)'CPC_ERROR_2: Aircraft inflight.'
            ENDIF 
!            write(*,*)'CPC_d_i_add_emis= ',d_i_add_emis

!-------AIRCRAFT ARRIVES AT A_CITY
        ELSEIF((d_p2_ac_routes_desc_old > 1.0_dp)                 .and. &
               (d_p2_ac_routes_desc_old < d_p2_ac_routes_desc_new).and. &
               (d_p2_ac_routes_desc_new == REAL(nwaypoints,DP)))THEN 
!               write(*,*)'CPC_Aircraft arrives at A_city.' 
            IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) <= d_p2_ac_routes_desc_new)THEN
               IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) /= d_p2_ac_routes_desc_old)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_potcov)*&
                                                            (DBLE(CEILING(d_p2_ac_routes_desc_old)) - d_p2_ac_routes_desc_old)
               ENDIF
               DO i=1,nwaypoints
                  IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) <= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_potcov)
                  ENDIF
                  !See 20121115-2, case(III)
                  IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) >= d_p2_ac_routes_desc_new)exit  !
               ENDDO
            ELSE
               write(*,*)'CPC_ERROR: Aircraft arrives at A_city.'
            ENDIF 
!            write(*,*)'CPC_d_i_add_emis= ',d_i_add_emis

!-------AIRCRAFT FLY FROM D_CITY TO A_CITY WITHIN 1 TIME STEP
        ELSEIF((d_p2_ac_routes_desc_old == 1.0_dp).and. &
               (d_p2_ac_routes_desc_new == REAL(nwaypoints,DP)))THEN
!               write(*,*)'CPC_Aircraft fly from D_ to A_city within 1 time step.' 
               DO i=1,nwaypoints-1
                  d_i_add_emis = d_i_add_emis + 1
                  d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lon)
                  d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lat)
                  d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_alt)
                  d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_potcov)
               ENDDO
!               write(*,*)'CPC_d_i_add_emis= ',d_i_add_emis

!-------ERROR MESSAGE
        ELSE
           write(*,*)'CPC_ERROR: Fill_emission_field.'
        ENDIF

        !DUMMY
!       d_add_glob_emis(:,:) = 1.0_dp + d_add_glob_emis(:,:)
!       d_glob_emis(:,:,:) = 1.0_dp + d_glob_emis(:,:,:)
!        write(*,*)"CPC_fill_emission_field."
!------------------------------------------
     CASE(ATR20O3)
!------------------------------------------
!-------AIRCRAFT IS ASSIGNED DEPARTURE TIME
        IF((d_p2_ac_routes_desc_old == 1.0_dp).and. &
           (d_p2_ac_routes_desc_new == 1.0_dp))THEN
!           write(*,*)'ATR20O3_Aircraft gets departure time.' 
           !Nothing: d_i_add_emis=0 

!-------AIRCRAFT LEAVES FROM D_CITY
        ELSEIF((d_p2_ac_routes_desc_old == 1.0_dp)                .and. &
               (d_p2_ac_routes_desc_old < d_p2_ac_routes_desc_new).and. &
               (d_p2_ac_routes_desc_new < REAL(nwaypoints,DP)))THEN 
!               write(*,*)'ATR20O3_Aircraft leaves from D_city.' 
               DO i=1,nwaypoints
                  IF(d_p2_ac_routes_desc_old + REAL(i,DP) < d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_atr20o3)
                  ELSEIF(d_p2_ac_routes_desc_old + REAL(i,DP) >= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_atr20o3)*&
                                                            (d_p2_ac_routes_desc_new - REAL(i,DP))
                     exit  !
                  ELSE
                     write(*,*)'ATR20O3_ERROR: Aircraft leaves from D_city.'
                  ENDIF
               ENDDO
!               write(*,*)'ATR20O3_d_i_add_emis= ',d_i_add_emis

!-------AIRCRAFT IS INFLIGHT
        ELSEIF((d_p2_ac_routes_desc_old > 1.0_dp)                 .and. &
               (d_p2_ac_routes_desc_old < d_p2_ac_routes_desc_new).and. &
               (d_p2_ac_routes_desc_new < REAL(nwaypoints,DP)))THEN 
!               write(*,*)'ATR20O3_Aircraft is inflight.' 
            IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) < d_p2_ac_routes_desc_new)THEN
               IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) /= d_p2_ac_routes_desc_old)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_atr20o3)*&
                                                            (DBLE(CEILING(d_p2_ac_routes_desc_old)) - d_p2_ac_routes_desc_old)
               ENDIF
               DO i=1,nwaypoints
                  IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) < d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_atr20o3)
                  ELSEIF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) >= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_atr20o3)*&
                                                            (d_p2_ac_routes_desc_new - &
                                                             DBLE(CEILING(d_p2_ac_routes_desc_old)) - REAL(i-1,DP))
                     exit  !
                  ELSE
                     write(*,*)'ATR20O3_ERROR_1: Aircraft inflight.'
                  ENDIF
               ENDDO
            ELSEIF(DBLE(CEILING(d_p2_ac_routes_desc_old)) >= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_atr20o3)*&
                                                            (d_p2_ac_routes_desc_new - d_p2_ac_routes_desc_old)
            ELSE
                     write(*,*)'ATR20O3_ERROR_2: Aircraft inflight.'
            ENDIF 
!            write(*,*)'ATR20O3_d_i_add_emis= ',d_i_add_emis

!-------AIRCRAFT ARRIVES AT A_CITY
        ELSEIF((d_p2_ac_routes_desc_old > 1.0_dp)                 .and. &
               (d_p2_ac_routes_desc_old < d_p2_ac_routes_desc_new).and. &
               (d_p2_ac_routes_desc_new == REAL(nwaypoints,DP)))THEN 
!               write(*,*)'ATR20O3_Aircraft arrives at A_city.' 
            IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) <= d_p2_ac_routes_desc_new)THEN
               IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) /= d_p2_ac_routes_desc_old)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_atr20o3)*&
                                                            (DBLE(CEILING(d_p2_ac_routes_desc_old)) - d_p2_ac_routes_desc_old)
               ENDIF
               DO i=1,nwaypoints
                  IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) <= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_atr20o3)
                  ENDIF
                  !See 20121115-2, case(III)
                  IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) >= d_p2_ac_routes_desc_new)exit  !
               ENDDO
            ELSE
               write(*,*)'ATR20O3_ERROR: Aircraft arrives at A_city.'
            ENDIF 
!            write(*,*)'ATR20O3_d_i_add_emis= ',d_i_add_emis

!-------AIRCRAFT FLY FROM D_CITY TO A_CITY WITHIN 1 TIME STEP
        ELSEIF((d_p2_ac_routes_desc_old == 1.0_dp).and. &
               (d_p2_ac_routes_desc_new == REAL(nwaypoints,DP)))THEN
!               write(*,*)'ATR20O3_Aircraft fly from D_ to A_city within 1 time step.' 
               DO i=1,nwaypoints-1
                  d_i_add_emis = d_i_add_emis + 1
                  d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lon)
                  d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lat)
                  d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_alt)
                  d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_atr20o3)
               ENDDO
!               write(*,*)'ATR20O3_d_i_add_emis= ',d_i_add_emis

!-------ERROR MESSAGE
        ELSE
           write(*,*)'ATR20O3_ERROR: Fill_emission_field.'
        ENDIF

        !DUMMY
!        write(*,*)"ATR20O3_fill_emission_field."
!------------------------------------------
     CASE(ATR20CH4)
!------------------------------------------
!-------AIRCRAFT IS ASSIGNED DEPARTURE TIME
        IF((d_p2_ac_routes_desc_old == 1.0_dp).and. &
           (d_p2_ac_routes_desc_new == 1.0_dp))THEN
!           write(*,*)'ATR20CH4_Aircraft gets departure time.' 

!-------AIRCRAFT LEAVES FROM D_CITY
        ELSEIF((d_p2_ac_routes_desc_old == 1.0_dp)                .and. &
               (d_p2_ac_routes_desc_old < d_p2_ac_routes_desc_new).and. &
               (d_p2_ac_routes_desc_new < REAL(nwaypoints,DP)))THEN 
!               write(*,*)'ATR20CH4_Aircraft leaves from D_city.' 
               DO i=1,nwaypoints
                  IF(d_p2_ac_routes_desc_old + REAL(i,DP) < d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_atr20ch4)
                  ELSEIF(d_p2_ac_routes_desc_old + REAL(i,DP) >= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_atr20ch4)*&
                                                            (d_p2_ac_routes_desc_new - REAL(i,DP))
                     exit  !
                  ELSE
                     write(*,*)'ATR20CH4_ERROR: Aircraft leaves from D_city.'
                  ENDIF
               ENDDO
!               write(*,*)'ATR20CH4_d_i_add_emis= ',d_i_add_emis

!-------AIRCRAFT IS INFLIGHT
        ELSEIF((d_p2_ac_routes_desc_old > 1.0_dp)                 .and. &
               (d_p2_ac_routes_desc_old < d_p2_ac_routes_desc_new).and. &
               (d_p2_ac_routes_desc_new < REAL(nwaypoints,DP)))THEN 
!               write(*,*)'ATR20CH4_Aircraft is inflight.' 
            IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) < d_p2_ac_routes_desc_new)THEN
               IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) /= d_p2_ac_routes_desc_old)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_atr20ch4)*&
                                                            (DBLE(CEILING(d_p2_ac_routes_desc_old)) - d_p2_ac_routes_desc_old)
               ENDIF
               DO i=1,nwaypoints
                  IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) < d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_atr20ch4)
                  ELSEIF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) >= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_atr20ch4)*&
                                                            (d_p2_ac_routes_desc_new - &
                                                             DBLE(CEILING(d_p2_ac_routes_desc_old)) - REAL(i-1,DP))
                     exit  !
                  ELSE
                     write(*,*)'ATR20CH4_ERROR_1: Aircraft inflight.'
                  ENDIF
               ENDDO
            ELSEIF(DBLE(CEILING(d_p2_ac_routes_desc_old)) >= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_atr20ch4)*&
                                                            (d_p2_ac_routes_desc_new - d_p2_ac_routes_desc_old)
            ELSE
                     write(*,*)'ATR20CH4_ERROR_2: Aircraft inflight.'
            ENDIF 
!            write(*,*)'ATR20CH4_d_i_add_emis= ',d_i_add_emis

!-------AIRCRAFT ARRIVES AT A_CITY
        ELSEIF((d_p2_ac_routes_desc_old > 1.0_dp)                 .and. &
               (d_p2_ac_routes_desc_old < d_p2_ac_routes_desc_new).and. &
               (d_p2_ac_routes_desc_new == REAL(nwaypoints,DP)))THEN 
!               write(*,*)'ATR20CH4_Aircraft arrives at A_city.' 
            IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) <= d_p2_ac_routes_desc_new)THEN
               IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) /= d_p2_ac_routes_desc_old)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_atr20ch4)*&
                                                            (DBLE(CEILING(d_p2_ac_routes_desc_old)) - d_p2_ac_routes_desc_old)
               ENDIF
               DO i=1,nwaypoints
                  IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) <= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_atr20ch4)
                  ENDIF
                  !See 20121115-2, case(III)
                  IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) >= d_p2_ac_routes_desc_new)exit  !
               ENDDO
            ELSE
               write(*,*)'ATR20CH4_ERROR: Aircraft arrives at A_city.'
            ENDIF 
!            write(*,*)'ATR20CH4_d_i_add_emis= ',d_i_add_emis

!-------AIRCRAFT FLY FROM D_CITY TO A_CITY WITHIN 1 TIME STEP
        ELSEIF((d_p2_ac_routes_desc_old == 1.0_dp).and. &
               (d_p2_ac_routes_desc_new == REAL(nwaypoints,DP)))THEN
!               write(*,*)'ATR20CH4_Aircraft fly from D_ to A_city within 1 time step.' 
               DO i=1,nwaypoints-1
                  d_i_add_emis = d_i_add_emis + 1
                  d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lon)
                  d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lat)
                  d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_alt)
                  d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_atr20ch4)
               ENDDO
!               write(*,*)'ATR20CH4_d_i_add_emis= ',d_i_add_emis

!-------ERROR MESSAGE
        ELSE
           write(*,*)'ATR20CH4_ERROR: Fill_emission_field.'
        ENDIF

        !DUMMY
!        write(*,*)"ATR20CH4_fill_emission_field."
!------------------------------------------
     CASE(ATR20H2O)
!------------------------------------------
!-------AIRCRAFT IS ASSIGNED DEPARTURE TIME
        IF((d_p2_ac_routes_desc_old == 1.0_dp).and. &
           (d_p2_ac_routes_desc_new == 1.0_dp))THEN
!           write(*,*)'ATR20H2O_Aircraft gets departure time.' 

!-------AIRCRAFT LEAVES FROM D_CITY
        ELSEIF((d_p2_ac_routes_desc_old == 1.0_dp)                .and. &
               (d_p2_ac_routes_desc_old < d_p2_ac_routes_desc_new).and. &
               (d_p2_ac_routes_desc_new < REAL(nwaypoints,DP)))THEN 
!               write(*,*)'ATR20H2O_Aircraft leaves from D_city.' 
               DO i=1,nwaypoints
                  IF(d_p2_ac_routes_desc_old + REAL(i,DP) < d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_atr20h2o)
                  ELSEIF(d_p2_ac_routes_desc_old + REAL(i,DP) >= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_atr20h2o)*&
                                                            (d_p2_ac_routes_desc_new - REAL(i,DP))
                     exit  !
                  ELSE
                     write(*,*)'ATR20H2O_ERROR: Aircraft leaves from D_city.'
                  ENDIF
               ENDDO
!               write(*,*)'ATR20H2O_d_i_add_emis= ',d_i_add_emis

!-------AIRCRAFT IS INFLIGHT
        ELSEIF((d_p2_ac_routes_desc_old > 1.0_dp)                 .and. &
               (d_p2_ac_routes_desc_old < d_p2_ac_routes_desc_new).and. &
               (d_p2_ac_routes_desc_new < REAL(nwaypoints,DP)))THEN 
!               write(*,*)'ATR20H2O_Aircraft is inflight.' 
            IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) < d_p2_ac_routes_desc_new)THEN
               IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) /= d_p2_ac_routes_desc_old)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_atr20h2o)*&
                                                            (DBLE(CEILING(d_p2_ac_routes_desc_old)) - d_p2_ac_routes_desc_old)
               ENDIF
               DO i=1,nwaypoints
                  IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) < d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_atr20h2o)
                  ELSEIF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) >= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_atr20h2o)*&
                                                            (d_p2_ac_routes_desc_new - &
                                                             DBLE(CEILING(d_p2_ac_routes_desc_old)) - REAL(i-1,DP))
                     exit  !
                  ELSE
                     write(*,*)'ATR20H2O_ERROR_1: Aircraft inflight.'
                  ENDIF
               ENDDO
            ELSEIF(DBLE(CEILING(d_p2_ac_routes_desc_old)) >= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_atr20h2o)*&
                                                            (d_p2_ac_routes_desc_new - d_p2_ac_routes_desc_old)
            ELSE
                     write(*,*)'ATR20H2O_ERROR_2: Aircraft inflight.'
            ENDIF 
!            write(*,*)'ATR20H2O_d_i_add_emis= ',d_i_add_emis

!-------AIRCRAFT ARRIVES AT A_CITY
        ELSEIF((d_p2_ac_routes_desc_old > 1.0_dp)                 .and. &
               (d_p2_ac_routes_desc_old < d_p2_ac_routes_desc_new).and. &
               (d_p2_ac_routes_desc_new == REAL(nwaypoints,DP)))THEN 
!               write(*,*)'ATR20H2O_Aircraft arrives at A_city.' 
            IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) <= d_p2_ac_routes_desc_new)THEN
               IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) /= d_p2_ac_routes_desc_old)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_atr20h2o)*&
                                                            (DBLE(CEILING(d_p2_ac_routes_desc_old)) - d_p2_ac_routes_desc_old)
               ENDIF
               DO i=1,nwaypoints
                  IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) <= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_atr20h2o)
                  ENDIF
                  !See 20121115-2, case(III)
                  IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) >= d_p2_ac_routes_desc_new)exit  !
               ENDDO
            ELSE
               write(*,*)'ATR20H2O_ERROR: Aircraft arrives at A_city.'
            ENDIF 
!            write(*,*)'ATR20H2O_d_i_add_emis= ',d_i_add_emis

!-------AIRCRAFT FLY FROM D_CITY TO A_CITY WITHIN 1 TIME STEP
        ELSEIF((d_p2_ac_routes_desc_old == 1.0_dp).and. &
               (d_p2_ac_routes_desc_new == REAL(nwaypoints,DP)))THEN
!               write(*,*)'ATR20H2O_Aircraft fly from D_ to A_city within 1 time step.' 
               DO i=1,nwaypoints-1
                  d_i_add_emis = d_i_add_emis + 1
                  d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lon)
                  d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lat)
                  d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_alt)
                  d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_atr20h2o)
               ENDDO
!               write(*,*)'ATR20H2O_d_i_add_emis= ',d_i_add_emis

!-------ERROR MESSAGE
        ELSE
           write(*,*)'ATR20H2O_ERROR: Fill_emission_field.'
        ENDIF

        !DUMMY
!        write(*,*)"ATR20H2O_fill_emission_field."
!------------------------------------------
!------------------------------------------
     CASE(ATR20CO2)
!------------------------------------------
!-------AIRCRAFT IS ASSIGNED DEPARTURE TIME
        IF((d_p2_ac_routes_desc_old == 1.0_dp).and. &
           (d_p2_ac_routes_desc_new == 1.0_dp))THEN
!           write(*,*)'ATR20CO2_Aircraft gets departure time.' 

!-------AIRCRAFT LEAVES FROM D_CITY
        ELSEIF((d_p2_ac_routes_desc_old == 1.0_dp)                .and. &
               (d_p2_ac_routes_desc_old < d_p2_ac_routes_desc_new).and. &
               (d_p2_ac_routes_desc_new < REAL(nwaypoints,DP)))THEN 
!               write(*,*)'ATR20CO2_Aircraft leaves from D_city.' 
               DO i=1,nwaypoints
                  IF(d_p2_ac_routes_desc_old + REAL(i,DP) < d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_atr20co2)
                  ELSEIF(d_p2_ac_routes_desc_old + REAL(i,DP) >= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_atr20co2)*&
                                                            (d_p2_ac_routes_desc_new - REAL(i,DP))
                     exit  !
                  ELSE
                     write(*,*)'ATR20O2_ERROR: Aircraft leaves from D_city.'
                  ENDIF
               ENDDO
!               write(*,*)'ATR20CH4_d_i_add_emis= ',d_i_add_emis

!-------AIRCRAFT IS INFLIGHT
        ELSEIF((d_p2_ac_routes_desc_old > 1.0_dp)                 .and. &
               (d_p2_ac_routes_desc_old < d_p2_ac_routes_desc_new).and. &
               (d_p2_ac_routes_desc_new < REAL(nwaypoints,DP)))THEN 
!               write(*,*)'ATR20CO2_Aircraft is inflight.' 
            IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) < d_p2_ac_routes_desc_new)THEN
               IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) /= d_p2_ac_routes_desc_old)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_atr20co2)*&
                                                            (DBLE(CEILING(d_p2_ac_routes_desc_old)) - d_p2_ac_routes_desc_old)
               ENDIF
               DO i=1,nwaypoints
                  IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) < d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_atr20co2)
                  ELSEIF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) >= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_atr20co2)*&
                                                            (d_p2_ac_routes_desc_new - &
                                                             DBLE(CEILING(d_p2_ac_routes_desc_old)) - REAL(i-1,DP))
                     exit  !
                  ELSE
                     write(*,*)'ATR20CO2_ERROR_1: Aircraft inflight.'
                  ENDIF
               ENDDO
            ELSEIF(DBLE(CEILING(d_p2_ac_routes_desc_old)) >= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_atr20co2)*&
                                                            (d_p2_ac_routes_desc_new - d_p2_ac_routes_desc_old)
            ELSE
                     write(*,*)'ATR20CO2_ERROR_2: Aircraft inflight.'
            ENDIF 
!            write(*,*)'ATR20CO2_d_i_add_emis= ',d_i_add_emis

!-------AIRCRAFT ARRIVES AT A_CITY
        ELSEIF((d_p2_ac_routes_desc_old > 1.0_dp)                 .and. &
               (d_p2_ac_routes_desc_old < d_p2_ac_routes_desc_new).and. &
               (d_p2_ac_routes_desc_new == REAL(nwaypoints,DP)))THEN 
!               write(*,*)'ATR20CO2_Aircraft arrives at A_city.' 
            IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) <= d_p2_ac_routes_desc_new)THEN
               IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) /= d_p2_ac_routes_desc_old)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_atr20co2)*&
                                                            (DBLE(CEILING(d_p2_ac_routes_desc_old)) - d_p2_ac_routes_desc_old)
               ENDIF
               DO i=1,nwaypoints
                  IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) <= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_atr20co2)
                  ENDIF
                  !See 20121115-2, case(III)
                  IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) >= d_p2_ac_routes_desc_new)exit  !
               ENDDO
            ELSE
               write(*,*)'ATR20CO2_ERROR: Aircraft arrives at A_city.'
            ENDIF 
!            write(*,*)'ATR20CO2_d_i_add_emis= ',d_i_add_emis

!-------AIRCRAFT FLY FROM D_CITY TO A_CITY WITHIN 1 TIME STEP
        ELSEIF((d_p2_ac_routes_desc_old == 1.0_dp).and. &
               (d_p2_ac_routes_desc_new == REAL(nwaypoints,DP)))THEN
!               write(*,*)'ATR20CO2_Aircraft fly from D_ to A_city within 1 time step.' 
               DO i=1,nwaypoints-1
                  d_i_add_emis = d_i_add_emis + 1
                  d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lon)
                  d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lat)
                  d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_alt)
                  d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_atr20co2)
               ENDDO
!               write(*,*)'ATR20CO2_d_i_add_emis= ',d_i_add_emis

!-------ERROR MESSAGE
        ELSE
           write(*,*)'ATR20CO2_ERROR: Fill_emission_field.'
        ENDIF

        !DUMMY
!        write(*,*)"ATR20CO2_fill_emission_field."
!------------------------------------------


!------------------------------------------
     CASE(ATR20CPC)
!------------------------------------------
!-------AIRCRAFT IS ASSIGNED DEPARTURE TIME
        IF((d_p2_ac_routes_desc_old == 1.0_dp).and. &
           (d_p2_ac_routes_desc_new == 1.0_dp))THEN
!           write(*,*)'ATR20CPC_Aircraft gets departure time.' 

!-------AIRCRAFT LEAVES FROM D_CITY
        ELSEIF((d_p2_ac_routes_desc_old == 1.0_dp)                .and. &
               (d_p2_ac_routes_desc_old < d_p2_ac_routes_desc_new).and. &
               (d_p2_ac_routes_desc_new < REAL(nwaypoints,DP)))THEN 
!               write(*,*)'ATR20CPC_Aircraft leaves from D_city.' 
               DO i=1,nwaypoints
                  IF(d_p2_ac_routes_desc_old + REAL(i,DP) < d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_atr20cpc)
                  ELSEIF(d_p2_ac_routes_desc_old + REAL(i,DP) >= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_atr20cpc)*&
                                                            (d_p2_ac_routes_desc_new - REAL(i,DP))
                     exit  !
                  ELSE
                     write(*,*)'ATR20CPC_ERROR: Aircraft leaves from D_city.'
                  ENDIF
               ENDDO
!               write(*,*)'ATR20CPC_d_i_add_emis= ',d_i_add_emis

!-------AIRCRAFT IS INFLIGHT
        ELSEIF((d_p2_ac_routes_desc_old > 1.0_dp)                 .and. &
               (d_p2_ac_routes_desc_old < d_p2_ac_routes_desc_new).and. &
               (d_p2_ac_routes_desc_new < REAL(nwaypoints,DP)))THEN 
!               write(*,*)'ATR20CPC_Aircraft is inflight.' 
            IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) < d_p2_ac_routes_desc_new)THEN
               IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) /= d_p2_ac_routes_desc_old)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_atr20cpc)*&
                                                            (DBLE(CEILING(d_p2_ac_routes_desc_old)) - d_p2_ac_routes_desc_old)
               ENDIF
               DO i=1,nwaypoints
                  IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) < d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_atr20cpc)
                  ELSEIF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) >= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_atr20cpc)*&
                                                            (d_p2_ac_routes_desc_new - &
                                                             DBLE(CEILING(d_p2_ac_routes_desc_old)) - REAL(i-1,DP))
                     exit  !
                  ELSE
                     write(*,*)'ATR20CPC_ERROR_1: Aircraft inflight.'
                  ENDIF
               ENDDO
            ELSEIF(DBLE(CEILING(d_p2_ac_routes_desc_old)) >= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_atr20cpc)*&
                                                            (d_p2_ac_routes_desc_new - d_p2_ac_routes_desc_old)
            ELSE
                     write(*,*)'ATR20CPC_ERROR_2: Aircraft inflight.'
            ENDIF 
!            write(*,*)'ATR20CPC_d_i_add_emis= ',d_i_add_emis

!-------AIRCRAFT ARRIVES AT A_CITY
        ELSEIF((d_p2_ac_routes_desc_old > 1.0_dp)                 .and. &
               (d_p2_ac_routes_desc_old < d_p2_ac_routes_desc_new).and. &
               (d_p2_ac_routes_desc_new == REAL(nwaypoints,DP)))THEN 
!               write(*,*)'ATR20CPC_Aircraft arrives at A_city.' 
            IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) <= d_p2_ac_routes_desc_new)THEN
               IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) /= d_p2_ac_routes_desc_old)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_atr20cpc)*&
                                                            (DBLE(CEILING(d_p2_ac_routes_desc_old)) - d_p2_ac_routes_desc_old)
               ENDIF
               DO i=1,nwaypoints
                  IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) <= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_atr20cpc)
                  ENDIF
                  !See 20121115-2, case(III)
                  IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) >= d_p2_ac_routes_desc_new)exit  !
               ENDDO
            ELSE
               write(*,*)'ATR20CPC_ERROR: Aircraft arrives at A_city.'
            ENDIF 
!            write(*,*)'ATR20CPC_d_i_add_emis= ',d_i_add_emis

!-------AIRCRAFT FLY FROM D_CITY TO A_CITY WITHIN 1 TIME STEP
        ELSEIF((d_p2_ac_routes_desc_old == 1.0_dp).and. &
               (d_p2_ac_routes_desc_new == REAL(nwaypoints,DP)))THEN
!               write(*,*)'ATR20CH4_Aircraft fly from D_ to A_city within 1 time step.' 
               DO i=1,nwaypoints-1
                  d_i_add_emis = d_i_add_emis + 1
                  d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lon)
                  d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lat)
                  d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_alt)
                  d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_atr20cpc)
               ENDDO
!               write(*,*)'ATR20CPC_d_i_add_emis= ',d_i_add_emis

!-------ERROR MESSAGE
        ELSE
           write(*,*)'ATR20CPC_ERROR: Fill_emission_field.'
        ENDIF

        !DUMMY
!        write(*,*)"ATR20CPC_fill_emission_field."

!--------------------------------------------------------------
     CASE(ATR20TOT)
!--------------------------------------------------------------
!-------AIRCRAFT IS ASSIGNED DEPARTURE TIME
        IF((d_p2_ac_routes_desc_old == 1.0_dp).and. &
           (d_p2_ac_routes_desc_new == 1.0_dp))THEN
!           write(*,*)'ATR2020TOT_Aircraft gets departure time.' 

!-------AIRCRAFT LEAVES FROM D_CITY
        ELSEIF((d_p2_ac_routes_desc_old == 1.0_dp)                .and. &
               (d_p2_ac_routes_desc_old < d_p2_ac_routes_desc_new).and. &
               (d_p2_ac_routes_desc_new < REAL(nwaypoints,DP)))THEN 
!               write(*,*)'ATR20TOT_Aircraft leaves from D_city.' 
               DO i=1,nwaypoints
                  IF(d_p2_ac_routes_desc_old + REAL(i,DP) < d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_atr20tot)
                  ELSEIF(d_p2_ac_routes_desc_old + REAL(i,DP) >= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_atr20tot)*&
                                                            (d_p2_ac_routes_desc_new - REAL(i,DP))
                     exit  !
                  ELSE
                     write(*,*)'ATR20TOT_ERROR: Aircraft leaves from D_city.'
                  ENDIF
               ENDDO
!               write(*,*)'ATR20TOT_d_i_add_emis= ',d_i_add_emis

!-------AIRCRAFT IS INFLIGHT
        ELSEIF((d_p2_ac_routes_desc_old > 1.0_dp)                 .and. &
               (d_p2_ac_routes_desc_old < d_p2_ac_routes_desc_new).and. &
               (d_p2_ac_routes_desc_new < REAL(nwaypoints,DP)))THEN 
!               write(*,*)'ATR20TOT_Aircraft is inflight.' 
            IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) < d_p2_ac_routes_desc_new)THEN
               IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) /= d_p2_ac_routes_desc_old)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_atr20tot)*&
                                                            (DBLE(CEILING(d_p2_ac_routes_desc_old)) - d_p2_ac_routes_desc_old)
               ENDIF
               DO i=1,nwaypoints
                  IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) < d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_atr20tot)
                  ELSEIF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) >= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_atr20tot)*&
                                                            (d_p2_ac_routes_desc_new - &
                                                             DBLE(CEILING(d_p2_ac_routes_desc_old)) - REAL(i-1,DP))
                     exit  !
                  ELSE
                     write(*,*)'ATR20TOT_ERROR_1: Aircraft inflight.'
                  ENDIF
               ENDDO
            ELSEIF(DBLE(CEILING(d_p2_ac_routes_desc_old)) >= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_atr20tot)*&
                                                            (d_p2_ac_routes_desc_new - d_p2_ac_routes_desc_old)
            ELSE
                     write(*,*)'ATR20TOT_ERROR_2: Aircraft inflight.'
            ENDIF 
!            write(*,*)'ATR20TOT_d_i_add_emis= ',d_i_add_emis

!-------AIRCRAFT ARRIVES AT A_CITY
        ELSEIF((d_p2_ac_routes_desc_old > 1.0_dp)                 .and. &
               (d_p2_ac_routes_desc_old < d_p2_ac_routes_desc_new).and. &
               (d_p2_ac_routes_desc_new == REAL(nwaypoints,DP)))THEN 
!               write(*,*)'ATR20TOT_Aircraft arrives at A_city.' 
            IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) <= d_p2_ac_routes_desc_new)THEN
               IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) /= d_p2_ac_routes_desc_old)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(CEILING(d_p2_ac_routes_desc_old),props_atr20tot)*&
                                                            (DBLE(CEILING(d_p2_ac_routes_desc_old)) - d_p2_ac_routes_desc_old)
               ENDIF
               DO i=1,nwaypoints
                  IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) <= d_p2_ac_routes_desc_new)THEN
                     d_i_add_emis = d_i_add_emis + 1
                     d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lon)
                     d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_lat)
                     d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_alt)
                     d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + CEILING(d_p2_ac_routes_desc_old),props_atr20tot)
                  ENDIF
                  !See 20121115-2, case(III)
                  IF(DBLE(CEILING(d_p2_ac_routes_desc_old)) + REAL(i,DP) >= d_p2_ac_routes_desc_new)exit  !
               ENDDO
            ELSE
               write(*,*)'ATR20TOT_ERROR: Aircraft arrives at A_city.'
            ENDIF 
!            write(*,*)'ATR20TOT_d_i_add_emis= ',d_i_add_emis

!-------AIRCRAFT FLY FROM D_CITY TO A_CITY WITHIN 1 TIME STEP
        ELSEIF((d_p2_ac_routes_desc_old == 1.0_dp).and. &
               (d_p2_ac_routes_desc_new == REAL(nwaypoints,DP)))THEN
!               write(*,*)'ATR20TOT_Aircraft fly from D_ to A_city within 1 time step.' 
               DO i=1,nwaypoints-1
                  d_i_add_emis = d_i_add_emis + 1
                  d_add_glob_emis(d_i_add_emis,add_lon) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lon)
                  d_add_glob_emis(d_i_add_emis,add_lat) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_lat)
                  d_add_glob_emis(d_i_add_emis,add_alt) = d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_alt)
                  d_add_glob_emis(d_i_add_emis,add_emis)= d_ac_routes(i + IDINT(d_p2_ac_routes_desc_old),props_atr20tot)
               ENDDO
!               write(*,*)'ATR20TOT_d_i_add_emis= ',d_i_add_emis

!-------ERROR MESSAGE
        ELSE
           write(*,*)'ATR20TOT_ERROR: Fill_emission_field.'
        ENDIF

        !DUMMY
!        write(*,*)"ATR20TOT_fill_emission_field."

!Yin_20170423
!-----------------------------------------
     ENDSELECT

     !USAGE CONFIRMATION
!     write(*,*)"fill_emission_field subroutine is done"

  END SUBROUTINE fill_emission_field
  ! =========================================================================

  ! =========================================================================
  SUBROUTINE airtraf_read_nml_ctrl(status, iou)

    ! ------------------------------------------------------------------
    ! This routine is used to read the CTRL-namelist of the submodel.
    ! This subroutine is executed when starting very first calculations,
    ! also restarting.
    ! ------------------------------------------------------------------

    ! MESSy INTERFACE
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    NAMELIST /CTRL/ nwaypoints, ngroutes, nwaypoints_out, & 
                    ngroutes_out, option_traj_calc,       &
                    option_output, option_max_ac_reached, &
                    lupdate_traj, ldaily_fp

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER       :: substr='airtraf_read_nml_ctrl'
    LOGICAL                           :: lex          ! file exists ?
    INTEGER                           :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR
    
    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist
    
    ! ### ADD HERE DIAGNOSTIC OUPUT FOR LOG-FILE
    !refer to the /CPL part of smil messy_airtraf_e5.f90
    WRITE(*,*) 'nwaypoints     = ',nwaypoints
    WRITE(*,*) 'ngroutes       = ',ngroutes
    WRITE(*,*) 'nwaypoints_out = ',nwaypoints_out
    WRITE(*,*) 'ngroutes_out   = ',ngroutes_out
    SELECT CASE(option_traj_calc) 
    CASE(0)
       WRITE(*,*) 'All aircrafts fly along Great circle trajectories.'
    CASE(1)
       WRITE(*,*) 'ALL aircrafts fly along Wind optimal trajectories'
    CASE(2)
       WRITE(*,*) 'ALL aircrafts fly along Fuel optimal trajectories'
    CASE(3)
       WRITE(*,*) 'ALL aircrafts fly along NOx optimal trajectories'
    CASE(4)
       WRITE(*,*) 'ALL aircrafts fly along H2O optimal trajectories'
    !Yin_20170423
    CASE(5)
       WRITE(*,*) 'ALL aircrafts fly along Contrail optimal trajectories'    
    !Yin_20170423
    CASE(6)
       WRITE(*,*) 'ALL aircrafts fly along Cost optimal trajectories'
    CASE(7)
       WRITE(*,*) 'ALL aircrafts fly along COC optimal trajectories'
!   !Yin_20170423
!   CASE(8)
!      WRITE(*,*) 'ALL aircrafts fly along multiobjective optimal trajectories'
!   !Yin_20170423 
    CASE DEFAULT
       WRITE(*,*) 'option_traj_calc is NOT defined.'
    ENDSELECT
    SELECT CASE(option_output) 
    CASE(0)
       WRITE(*,*) 'option_output         =0: Standatd output.'
    CASE(1)
       WRITE(*,*) 'option_output         =1: AC locations in addition'
    CASE DEFAULT
       WRITE(*,*) 'option_output is NOT defined.'
    ENDSELECT
    SELECT CASE(option_max_ac_reached) 
    CASE(0)
       WRITE(*,*) 'option_max_ac_reached =0: Stop.'
    CASE(1)
       WRITE(*,*) 'option_max_ac_reached =1: Ignore ac flight.'
    CASE DEFAULT
       WRITE(*,*) 'option_max_ac_reached is NOT defined.'
    ENDSELECT
    if (lupdate_traj) then
       WRITE(*,*) "lupdate_traj is .TRUE."
    else
       WRITE(*,*) "lupdate_traj is .FALSE."
    endif
    if (ldaily_fp) then
       WRITE(*,*) "Repeat flight next day"
    else
       WRITE(*,*) "NOT repeat flight next day"
    endif
    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR
    
  END SUBROUTINE airtraf_read_nml_ctrl
  ! =========================================================================

! **********************************************************************
END MODULE messy_airtraf
! **********************************************************************

