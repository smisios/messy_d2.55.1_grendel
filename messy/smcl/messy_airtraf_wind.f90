! **********************************************************************
! 
! SUB_SUBMODEL ROUTINES FOR MESSy SUBMODEL airtraf
!
! THIS SUB_SUBMODEL IS USED TO CALCULATE OPTIMAL TRAJECTORIES
!
! Authors : Hiroshi Yamashita, DLR-IPA, March 2012
!           Volker Grewe, DLR-IPA, March 2012
!           Feijia Yin, TUDelft, Dec. 2016
!
! References:
!
! * None yet
!
! see SMCL for more details
! **********************************************************************
! Objective: Simulation of an optimal trajectory, which has minimal objective function value, 
!            on arbitrary city pairs.  
!            - Flight routing options: 
!                    (6) minimal contrail potential coverage trajectories, FY
! Input:    Flighplan information
! Output:   Flight Time table including waypoints
! Method: 
!
! - Input ranges  : -90 =< gc_D_lat =< +90, -180 =< gc_D_lon =< +180 [deg]
!                   -90 =< gc_A_lat =< +90, -180 =< gc_A_lon =< +180 [deg]
!
! - Output ranges : Waypoints output is,
!                   - OUT_SWITCH = 0 (default) -> -180 =< lon =< +180 [deg] (ALL FL_DIR)
!                   - OUT_SWITCH = 1 (option)  ->    0 =< lon =< +360 [deg] (FL_DIR=1,3) 
!
! - These lon & lat are in decimal degrees.(NOT DD:MM:SS expressions)  
!
! - Flight time is accounted and integrated by 'FLIGHT_TIME(:)' at every waypoint in [sec]. 
!   The 'T_FLIGHT_TIME(DP)' indicates a total flight time from D_city to A_city in [sec].
!   Flight time and flight distance are calculated in spherical coordinates.
!
! - T_FLIGHT_DISTANCE  : Total flight time from D_ to A_city, [km]
!   FLIGHT_DISTANCE(:) : Flight distance between waypoints, [km]  
!
! - Aircrafts fly along GC:
!     - strady flight (AC_MACH = const)
!     - constant altitude [m] (ALT = const) 
!     - constant sound speed [m/s] (SOUND_V = const) 
!     - no winds 
!      
!     
! - Included subroutines:
!     - (gc_distance) not used in default
!     - gc_waypoints
!     - gc_ac_traj
!
!
!
!=========================================================================
!15 MARCH 2012, H. Yamashita, Oberpfaffenhofen
!14 AUGUST 2012, updated
!05 Dec. 2016, Feijia Yin, TUDelft 
!=========================================================================
!
! **********************************************************************
MODULE messy_airtraf_wind
! **********************************************************************

!USE ONLY
 USE messy_main_constants_mem, ONLY: DP, PI, DTR, OneDay, radius_earth, g, atm2Pa
 USE messy_main_tools,         ONLY: nn_index !20131213-1
 USE messy_airtraf_gc,         ONLY: AC_MACH &  !delete ALT 20140522-1,20140524-10, BADA_FUEL_NOM 20140524-7, SOUND_V 20140524-9 
                                    ,props_lon,props_lat,props_alt,props_time        &  !,total_flight_time_s 20140524-4
                                    ,props_ac_speed,props_dist,props_fuel_use,props_emis_nox &
                                    ,props_emis_h2o   &   !20131216-1   ,DIV_WAY 20140721-4
                                    ,props_potcov ,props_atr20o3,props_atr20ch4&
                                    ,props_atr20h2o,props_atr20cpc,props_atr20co2,props_atr20tot   !Yin_20170423
! USE messy_airtraf_e5,  ONLY:D_lon, D_lat, D_time, A_lon, A_lat 

! USE messy_airtraf,  ONLY: NOx, H2O, DIST, FUEL
! USE messy_airtraf,  ONLY: nwaypoints
 USE messy_airtraf_tools_ga_armogaset,       ONLY: Lower_alt, Upper_alt    !HY20140428-3

  IMPLICIT NONE
  ! GLOBAL PARAMETERS 
  !Comments:If we can use the information on aircraft type for each route, 
  !         average flight speed of each aircraft type can be utilized as AC_MACH. 
!  DOUBLE PRECISION, PARAMETER, PUBLIC :: ALT     = 10058.40_dp !8839.2d0 !Flight altitude(=geopotential altitude,20131105-1), [m] 
!  DOUBLE PRECISION, PARAMETER, PUBLIC :: AC_MACH = 0.82_dp  !0.8d0       !Flight Mach number(=Mcr from A333_ptf), [-] see20130719-1
                                                                          !AC_MACH should be from Aircraft table(BADA A333_ptf) 
                                                                          !0.8d0(opt,20130131-1)
!  DOUBLE PRECISION, PARAMETER, PUBLIC :: SOUND_V = 299.463_dp !340.3d0   !Speed of sound from BADA_ISA,[m/s] see20130719-1 
!  DOUBLE PRECISION, PUBLIC            :: SOUND_V                         !Speed of sound,[m/s] see20131112-1 
!  DOUBLE PRECISION, PUBLIC            :: BADA_FUEL_NOM                   !Fuel from BADA A333_ptf,[kg/s] see20131112-1 
!  DOUBLE PRECISION, PUBLIC            :: total_flight_time_s             ![s],20131003-3,20131114-2
 
  ! The parameter should be used from SMIL
!  INTEGER, PARAMETER                  :: DIV_WAY = 11    !(DIV_WAY = nwaypoints - 1)
                                                 ! this variable should be determined automatically,see20130712-1.
  INTEGER, PRIVATE                     :: DIV_WAY         !20140721-5
  INTEGER, PARAMETER, PRIVATE          :: D_lon=1, D_lat=2, D_time=3, A_lon=4, A_lat=5
  REAL(DP), PARAMETER, PRIVATE         :: adiabatic_index_air = 1.4_dp         ![-]20140724-3 
  REAL(DP), PARAMETER, PRIVATE         :: gas_constant_air    = 287.05287_dp   ![m^2/(k*s^2)]20140724-3
  INTEGER, PARAMETER                   :: NUM_CP_ALT   = 5
  INTEGER, PARAMETER                   :: NUM_CP_LONLAT= 3

  !Parameters for ac_routes properties
!  INTEGER, PARAMETER, PUBLIC                  :: props_lon       = 1
!  INTEGER, PARAMETER, PUBLIC                  :: props_lat       = 2
!  INTEGER, PARAMETER, PUBLIC                  :: props_alt       = 3
!  INTEGER, PARAMETER, PUBLIC                  :: props_time      = 4
!  INTEGER, PARAMETER, PUBLIC                  :: props_ac_speed  = 5
!  INTEGER, PARAMETER, PUBLIC                  :: props_dist      = 6   !subroutine calculate_emissions_along_trajectory
!  INTEGER, PARAMETER, PUBLIC                  :: props_fuel_use  = 7   !subroutine calculate_emissions_along_trajectory
!  INTEGER, PARAMETER, PUBLIC                  :: props_emis_nox  = 8   !subroutine calculate_emissions_along_trajectory
!  INTEGER, PARAMETER, PUBLIC                  :: props_emis_h2o  = 9   !subroutine calculate_emissions_along_trajectory

  ! PRIVATE PARAMETERS
  INTEGER, PARAMETER, PRIVATE :: ISWITCH = 2          !0: Spherical law of cosines
                                                      !1: Haversine formula
                                                      !2: Vincenty formula
  INTEGER, PARAMETER, PRIVATE :: OUT_SWITCH = 0       !0: -180 =< lon =< +180 (Default) 
                                                      !1:    0 =< lon =< +360 (only cases FL_DIR = 1 and 3)
  REAL(DP), PARAMETER, PUBLIC :: obj_dummy = 1.0d10   !HY20140428-3,20160606-3,20170215-6,20170217-5,20170221-1,20170222-5
  REAL(DP), PARAMETER, PUBLIC :: Fuel_price   = 1.545_dp   !Average fuel Price in March 2017(US Doller/US Gallon)20170224-5
  REAL(DP), PARAMETER, PUBLIC :: Fuel_density = 6.71_dp    !Fuel_density[Pounds/US Gallon]
! REAL(DP), PARAMETER, PUBLIC :: Ct        = 2710.0_dp     !Unit time costs[US Dollar/h], 20170428-4
! REAL(DP), PARAMETER, PUBLIC :: Co        = 0.0_dp        !Other costs, 20170428-4
! REAL(DP)                    :: Cf                        !Unit fuel costs[Cents/lbs], 20170428-4
 
  ! PRIVATE
  !REAL(DP), DIMENSION(:), PRIVATE :: THETA_OUT(DIV_WAY+1), PHI_OUT(DIV_WAY+1)  !20140721-4
!HY REAL(DP), PRIVATE   :: AC_V  20140522-4 
  REAL(DP), PRIVATE           :: RARAD1,RARAD2,DCRAD1,DCRAD2 
  REAL(DP), PRIVATE           :: DELDC2,DELRA,DELRA2,SINDIS,SINDIS1,SINDIS2
  REAL(DP), PRIVATE           :: THETAR, PHIR
! REAL(DP), PRIVATE           :: gc_D_time
  INTEGER, PRIVATE            :: FL_DIR        !see 20111219-1,20130514-2-left

  REAL(DP)                       :: gc_D_lat        !latitude of departure city 
  REAL(DP)                       :: gc_D_lon        !longitude of departure city
  REAL(DP), PRIVATE              :: gc_D_time       !20140728-1
  REAL(DP)                       :: gc_A_lat        !latitude of arrivar city 
  REAL(DP)                       :: gc_A_lon        !longitude of arrival city

! REAL(DP), DIMENSION(:) :: FLIGHT_TIME(DIV_WAY+1), FLIGHT_TIME_JULIAN(DIV_WAY+1) !20140721-4 
! REAL(DP), DIMENSION(:) :: FLIGHT_DISTANCE(DIV_WAY+1), FLIGHT_SPEED(DIV_WAY+1)   !opt20130214-4,opt20130222-9,20140721-4
! REAL(DP)               :: T_FLIGHT_TIME, A_time, T_FLIGHT_DISTANCE !A_time in Julian date,20140728-1
  REAL(DP)               :: temp_gc_A_lon        !For FL_DIR=3, 180 =< temp_gc_A_lon =< 540
  REAL(DP)               :: temp_gc_D_lon        !For FL_DIR=1, 180 =< temp_gc_D_lon =< 540

  ! Parameter for Single Opbjective Optimization Problem on 'MINIMUM FLIGHT TIME': Wind opt option. 
  ! The DV(:) is used to express arbitray trajectories.
! REAL(DP), DIMENSION(:) :: XXN(DIV_WAY+1),YYN(DIV_WAY+1),ZZN(DIV_WAY+1)   !20140721-5
  REAL(DP)               :: DV(11)  !HY20140415-6
  REAL(DP)               :: lon_dv7, lon_dv8, lon_dv9, lon_dv10, lon_dv11

!HY20140415-7,20140416-1
! temp_best_traj     : best traj propaties in current generation
! temp_best_fl_time  : best(minimum) flight time in current generation
! temp_best_gen      : generation number including the current best traj
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: temp_best_traj          !(DIV_WAY+1,7) 20140721-5,20170607-1
  REAL(DP)                              :: temp_best_fl_time
  REAL(DP)                              :: temp_best_fuel_use      !20160606-5 
  REAL(DP)                              :: temp_best_nox_emis      !20170216-2 
  REAL(DP)                              :: temp_best_h2o_emis      !20170222-4 
  !Yin_20170423
  REAL(DP)                              :: temp_best_potcov        !Total CPC [km] 20170608-1
  REAL(DP)                              :: temp_best_atr20         !Total ATR20 [K] 20170801
  REAL(DP)                              :: temp_best_costclim      !20170801
  REAL(DP)                              :: temp_best_costcpc      !20170801
  !Yin_20170423
  REAL(DP)                              :: temp_best_cost          !20170224-4 
  REAL(DP)                              :: temp_best_coc           !20170316-4 
  INTEGER                               :: temp_best_gen
  REAL(DP), DIMENSION(:), ALLOCATABLE   :: temp_vtas_along_traj    !(DIV_WAY+1)    ![m/s] 20140524-5,20140721-5 
  !Yin_20170423
  !REAL(DP), DIMENSION(:), ALLOCATABLE  :: temp_cpc_along_traj     !(DIV_WAY+1)    ![km] 20170606-4,20170607-2
  !Yin_20170423
  INTEGER                               :: w_fl_direction          !20141216-2

!Output Vground/Vtas for 3 pair-cities,20150312-1
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: temp_best_vtas          !(DIV_WAY+1,101) 
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: temp_best_vground       !(DIV_WAY+1,101) 
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: temp_best_rho           !(DIV_WAY+1,101) 20160623-2 
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: temp_best_cpc           !(DIV_WAY+1,101) 20170608-1

!information of routing option 20170216-5
  INTEGER                               :: w_option_traj_calc_g   

  ! PUBLIC SUBROUTINES (to be called from messy_airtraf_tools_ga.f90)
! PUBLIC :: gc_calculate
! PUBLIC :: wind_opt_calculate
!20140328-3
  PUBLIC :: calcobj   !20140328-4,20140417-2
  PUBLIC :: cost_calc !20170428-6
  PUBLIC :: coc_calc  !20170428-1
!Yin_20170801+  
  PRIVATE :: atr20_calc  
  PRIVATE :: cost_atr20_calc  
  PRIVATE :: contrail_cost_calc  
!Yin_20170801-
  ! PRIVATE SUBROUTINES
! PRIVATE :: gc_distance 
! PRIVATE :: gc_waypoints 
  PRIVATE :: determine_lon_dvs  !20140417-2 
  PRIVATE :: arbitrary_traj 
  PRIVATE :: gc_ac_traj 
  PRIVATE :: fuel_calc          !20160531-2 
!HY  PRIVATE :: inner_point     !20140508-3 
  PRIVATE :: wind_effect        !20140509-1,20140721-1
  CONTAINS
!----------------------------------------------------------------------
!20140417-2,20140328-3
!HY20140420-2  SUBROUTINE calcobj(gc_ac_routes, gc_p2_ac_routes_desc, w_philon, w_philat,  &
!                                 w_zgl_geopot_3d, w_uwind_g, w_vwind_g, w_v_z_g,          &
!                                 dvs,num_dvs,obj_val,num_obj,const_val,num_const,num_gen,num_pop) !20140417-2
  SUBROUTINE calcobj(dvs,num_dvs,obj_val,num_obj,const_val,num_const,num_gen,num_pop,  &
                     w_nwaypoints, w_p2_ac_routes_desc, w_philon, w_philat, w_zgl_geopot_3d,         &
                     w_uwind_g, w_vwind_g, w_v_z_g, w_t_scb_g, w_cpc_g, &  !20140420-2,20140522-3,20140721-7 
                     w_option_traj_calc, w_rho_air_dry_3d_g,   &  !20160302-3,20160525-1,20170607-3     
                     w_press_3d_g,w_ATR20O3_3d_g,w_ATR20CH4_3d_g,&
                     w_ATR20H2O_3d_g,w_ATR20CPC_3d_g,w_ATR20CO2_3d_g)             !20170801
!----------------------------------------------------------------------
!
! This subroutine contains the subroutines on wind_opt calculate and controls I/O.
! This SUBROUTINE is called by messy_airtraf_tools_ga.f90.
!
  IMPLICIT NONE
!HY20140420-2  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: gc_ac_routes
  INTEGER, INTENT(IN)                   :: w_nwaypoints         !20140721-7
  INTEGER, INTENT(IN)                   :: w_option_traj_calc   !20160525-1, routing option
  REAL(DP), DIMENSION(:),   INTENT(IN)  :: w_p2_ac_routes_desc  !HY20140420-2
  REAL(DP), DIMENSION(:),   INTENT(IN)  :: w_philon             !20140212-2
  REAL(DP), DIMENSION(:),   INTENT(IN)  :: w_philat
  REAL(DP), DIMENSION(:,:,:), INTENT(IN):: w_zgl_geopot_3d   !global field [m^2/s^2] 
  REAL(DP), DIMENSION(:,:,:), INTENT(IN):: w_uwind_g         !global field [m/s] 
  REAL(DP), DIMENSION(:,:,:), INTENT(IN):: w_vwind_g         !global field [m/s]
  REAL(DP), DIMENSION(:,:,:), INTENT(IN):: w_v_z_g           !global field [m/s]
  REAL(DP), DIMENSION(:,:,:), INTENT(IN):: w_t_scb_g         !global field [K]
  !Yin_20170423
  REAL(DP), DIMENSION(:,:,:), INTENT(IN):: w_cpc_g           !global field [-]
  !Yin_20170423
  REAL(DP), DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: w_rho_air_dry_3d_g    !global field [kg/m^3] 20160302-3,20170607-3
  REAL(DP), DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: w_press_3d_g          !global field [Pa] 20170215-4,20170607-3
  !Yin_20170801+
  REAL(DP), DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: w_ATR20O3_3d_g          
  REAL(DP), DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: w_ATR20CH4_3d_g         
  REAL(DP), DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: w_ATR20H2O_3d_g         
  REAL(DP), DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: w_ATR20CO2_3d_g         
  REAL(DP), DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: w_ATR20CPC_3d_g         
  !Yin_20170810-
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: w_zalt            ! geopotential altitude[m],20131105-1
  REAL(DP)                              :: wind_along_traj(w_nwaypoints,3)  ![m/s] 20140429-2,20140721-7 
  REAL(DP)                              :: vtas_along_traj(w_nwaypoints)    ![m/s] 20140522-4,20140721-7 
  !Yin_20170423+
  REAL(DP)                              :: cpc_along_traj(w_nwaypoints)     ![-] fraction 20161205,20161215 
  REAL(DP)                              :: ATR20O3_along_traj(w_nwaypoints)     ![K/kg(NO2)] ATR20 ozone,20170801 
  REAL(DP)                              :: ATR20CH4_along_traj(w_nwaypoints)     ![K/kg(NO2)] ATR20 CH4,20170801 
  REAL(DP)                              :: ATR20H2O_along_traj(w_nwaypoints)     ![K/kg(NO2)]ATR20 H2O, 20170801  
  REAL(DP)                              :: ATR20CPC_along_traj(w_nwaypoints)     ![K/km(contrail)]ATR20 contrail,20170801 
  REAL(DP)                              :: ATR20CO2_along_traj(w_nwaypoints)     ![K/km(contrail)]ATR20 contrail,20170801 
  !Yin_20170423-
  REAL(DP)                              :: rho_along_traj(w_nwaypoints)     ![kg/m^3] 20160531-1 
  REAL(DP)                              :: press_along_traj(w_nwaypoints)   ![Pa] 20170215-5 
  REAL(DP)                              :: t_scb_along_traj(w_nwaypoints)   ![K] 20170215-5 

  REAL(DP)                              :: XXN(w_nwaypoints),YYN(w_nwaypoints),ZZN(w_nwaypoints)   !20140721-5,20140721-7
!HY20140524-9  REAL(DP), DIMENSION(:,:) :: BADA_DATA(8,5) !this should be given by Aircraft table.20130803-1,20131108-1
  REAL(DP)                              :: FLIGHT_TIME(w_nwaypoints), FLIGHT_TIME_JULIAN(w_nwaypoints) !20140721-7
  REAL(DP)                              :: FLIGHT_DISTANCE(w_nwaypoints), FLIGHT_SPEED(w_nwaypoints)  !,FUEL_USE(DIV_WAY+1)
  REAL(DP)                              :: contrail_coverage(w_nwaypoints)   ![km] 20170607-1
  REAL(DP)                              :: T_FLIGHT_TIME, A_time, T_FLIGHT_DISTANCE !A_time in Julian date
  !Yin_20170423 T_CPC is the total travel distance in contrail formation region, 20170607-2
  REAL(DP)                              :: T_CPC      ![km]   
  !Yin_20170423
  REAL(DP)                              :: temp_XXN   !20140605-1
!Objective function values
  REAL(DP)                              :: Total_fuel_use   !Obj value for fuel_opt routing option, 20160606-3
  REAL(DP)                              :: Total_nox_emis   !Obj value for NOx_opt routing option, 20170215-6
  REAL(DP)                              :: Total_h2o_emis   !Obj value for H2O_opt routing option, 20170222-3
  REAL(DP)                              :: Total_cost       !Obj value for Cost_opt routing option, 20170224-3
  REAL(DP)                              :: Total_coc        !Obj value for COC_opt routing option, 20170316-3
  !Yin_20170801+
  REAL(DP)                              :: Total_atr20      !Obj value for ATR20_opt routing option, 20170801
  REAL(DP)                              :: Total_atr20o3    !Obj value for ATR20_opt routing option, 20170801
  REAL(DP)                              :: Total_atr20ch4   !Obj value for ATR20_opt routing option, 20170801
  REAL(DP)                              :: Total_atr20h2o     !Obj value for ATR20_opt routing option, 20170801
  REAL(DP)                              :: Total_atr20cpc    !Obj value for ATR20_opt routing option, 20170801
  REAL(DP)                              :: Total_atr20co2    !Obj value for ATR20_opt routing option, 20170801
  REAL(DP)                              :: Total_costclim    !Obj value for Multiobj routing option, 20170801
  REAL(DP)                              :: Total_costcpc    !Obj value for Multiobj routing option, 20170801
  !Yin_20170801-
  INTEGER  :: i,j, w_i_lat, w_j_lon, w_i_alt   !20140212-2 
  INTEGER  :: nlev_w_zalt  !20140429-2
  !For optimizations:   !20140417-2
  INTEGER                       :: num_dvs,num_obj,num_const,num_pop !i,j, !20140328-3
  INTEGER                       :: nsum1,nsum2,ipop,idv,iof               !HY20140415-5
  INTEGER                       :: num_gen                                !HY20140415-2
  INTEGER                       :: wp_const                               !HY20140428-3 
  REAL(DP)              :: dvs(num_pop*num_dvs),obj_val(num_pop*num_obj),const_val(num_pop*num_const+1)
  !Output Vground/Vtas for 3 city-pairs,20150311-2
  REAL(DP)                      :: traj_for_nnindex(w_nwaypoints,3,101)    !20150311-2
  INTEGER                       :: i_add_nnindex     !0:NO execute   ,1:Execute 

  !COMMON DATA FOR ALL POPPULATION
  !BADA DATA input,A333_ptf,BADA_ISA_TABLE, 20131108-1, see20140524-9
!  DATA ((BADA_DATA(i,j),i=1,8),j=1,3) /280.0_dp,290.0_dp,310.0_dp,330.0_dp,350.0_dp,370.0_dp,390.0_dp,410.0_dp,         &  !FL[ft]
!                                       305.79_dp,304.48_dp,301.86_dp,299.21_dp,296.54_dp,295.07_dp,295.07_dp,295.07_dp, &  !a[m/s]
!                                       112.4_dp,108.7_dp,101.7_dp,95.5_dp,90.0_dp,85.5_dp,81.9_dp,79.0_dp/      !Fuel[nom, kg/min]
!  do i=1,8 
!     BADA_DATA(i,4)=BADA_DATA(i,1)*100.0_dp*0.30480_dp   ![ft] to [m],20131108-1
!     BADA_DATA(i,5)=BADA_DATA(i,3)/60.0_dp               ![kg/min] to [kg/s]
!     write(*,*)'BADA_DATA_check',(BADA_DATA(i,j),j=1,5)
!  enddo

  !Allocation:20140721-5
  !Deallocation is performed in messy_airtraf_tools_ga.f90 at line 170,20140722-1
  IF(.not.ALLOCATED(temp_best_traj))ALLOCATE(temp_best_traj(w_nwaypoints,10))       !20170607-1
  IF(.not.ALLOCATED(temp_vtas_along_traj))ALLOCATE(temp_vtas_along_traj(w_nwaypoints))

  !output Vground/Vtas for 3 city-pairs,20150312-2
  IF(.not.ALLOCATED(temp_best_vtas))ALLOCATE(temp_best_vtas(w_nwaypoints,101))
  IF(.not.ALLOCATED(temp_best_vground))ALLOCATE(temp_best_vground(w_nwaypoints,101))
  IF(.not.ALLOCATED(temp_best_rho))ALLOCATE(temp_best_rho(w_nwaypoints,101))
  IF(.not.ALLOCATED(temp_best_cpc))ALLOCATE(temp_best_cpc(w_nwaypoints,101))     !20170608-2

  nlev_w_zalt = SIZE(w_zgl_geopot_3d,dim=2)   !ALT: in 2nd dimention.
  ALLOCATE(w_zalt(nlev_w_zalt))               !20140429-2
!  write(*,*)'windf90_tscb_check:',shape(w_t_scb_g),w_t_scb_g(:,:,:)
!  write(*,*)'windf90_cpc_check:',shape(w_cpc_g),w_cpc_g(:,:,:)
  !Determine:DIV_WAY 20140721-6,7
  DIV_WAY = w_nwaypoints-1
! write(*,*)"w_nwaypoints_check",w_nwaypoints
! write(*,*)"nlev_w_zalt=",nlev_w_zalt
  !CALCULATE FOR EACH IPOP(ARBITRARY TRAJECOTRY ON 1 CITY-PAIR)
  !Loop for ipop, HY20140415-5 
  nsum1=1
  nsum2=1
  temp_best_fl_time  = obj_dummy  !=1.0d8 20140416-1,20170217-5
  temp_best_fuel_use = obj_dummy  !=1.0d8 20160606-5,20170217-5
  temp_best_nox_emis = obj_dummy  !=1.0d8 20170216-2,20170217-5
  temp_best_h2o_emis = obj_dummy  !=1.0d8 20170216-2,20170222-4
  !Yin_20170423
  temp_best_potcov   = obj_dummy  ![km] FY 20161215
  temp_best_atr20    = obj_dummy  ![K] FY 20170801
  temp_best_costclim = obj_dummy  ![dollar] FY 20170801
  temp_best_costcpc  = obj_dummy  ![dollar] FY 20170801
  !Yin_20170423
  temp_best_cost     = obj_dummy  !=1.0d8 20170224-4
  temp_best_coc      = obj_dummy  !=1.0d8 20170316-4

! add info of routing option to global variable 20170216-5,20170221-1
  w_option_traj_calc_g = w_option_traj_calc 
! write(*,*)"Check:routing_option",w_option_traj_calc,w_option_traj_calc_g

  !Output Vground/Vtas for 3 city-pairs, 20150311-2
  traj_for_nnindex=0.0_dp  !Initialization,20150311-2
  i_add_nnindex   =0

  !From line 241 to 303 were moved outside of Population loop. 20140725-1  

  !Definition of D_city and A_city
  !Flightplan information from SMCL 
  gc_D_lon  = w_p2_ac_routes_desc(D_lon)      !lon of departure-city   
  gc_D_lat  = w_p2_ac_routes_desc(D_lat)      !lat of departure-city   
  gc_D_time = w_p2_ac_routes_desc(D_time)     !departure time in Julian date 
  gc_A_lon  = w_p2_ac_routes_desc(A_lon)      !lon of arrival-city
  gc_A_lat  = w_p2_ac_routes_desc(A_lat)      !lat of arrival-city  
  temp_gc_A_lon = gc_A_lon + 360.0_dp
  temp_gc_D_lon = gc_D_lon + 360.0_dp

!  OPEN(40,file="north_check.dat",position='APPEND')  !20131203-3
!     if(gc_D_lat.lt.0.0_dp)write(40,*)"1",gc_D_lat,gc_D_lon,gc_A_lat,gc_A_lon
!     if(gc_A_lat.lt.0.0_dp)write(40,*)"2",gc_D_lat,gc_D_lon,gc_A_lat,gc_A_lon
!     if(gc_A_lat.lt.0.0_dp.and.gc_D_lat.lt.0.0_dp)write(40,*)"3",gc_D_lat,gc_D_lon,gc_A_lat,gc_A_lon
!  CLOSE(40)

  !Definition of the city being placed at higher lon. 
  IF(gc_D_lon.lt.gc_A_lon)THEN   !20111219-1
     RARAD1= gc_A_lon*DTR        !lon of higher city [rad]
     RARAD2= gc_D_lon*DTR        !lon of lower  city [rad]
     DCRAD1= gc_A_lat*DTR        !lat of higher city [rad]    
     DCRAD2= gc_D_lat*DTR        !lat of lower  city [rad]
     IF((RARAD1-RARAD2).le.PI)THEN
        FL_DIR= 0
     ELSE
        FL_DIR= 1
     ENDIF
  ELSE IF(gc_D_lon.gt.gc_A_lon)THEN
     RARAD1= gc_D_lon*DTR     
     RARAD2= gc_A_lon*DTR
     DCRAD1= gc_D_lat*DTR     
     DCRAD2= gc_A_lat*DTR     
     IF((RARAD1-RARAD2).le.PI)THEN 
        FL_DIR= 2
     ELSE
        FL_DIR= 3
     ENDIF
  ELSE     !(gc_D_lon.eq.gc_A_lon)
     IF(gc_D_lat.gt.gc_A_lat)THEN
        RARAD1= gc_D_lon*DTR     
        RARAD2= gc_A_lon*DTR
        DCRAD1= gc_D_lat*DTR     
        DCRAD2= gc_A_lat*DTR     
        FL_DIR= 4
     ELSE IF(gc_D_lat.lt.gc_A_lat)THEN
        RARAD1= gc_A_lon*DTR 
        RARAD2= gc_D_lon*DTR
        DCRAD1= gc_A_lat*DTR     
        DCRAD2= gc_D_lat*DTR 
        FL_DIR= 5
     ELSE     !(gc_D_lat.eq.gc_A_lat)
        write(*,*)'City Selection Error: They must be identical city!'
     ENDIF
  ENDIF
  !20141216-2
  w_fl_direction = FL_DIR       

  !opt20130205-3
  CALL determine_lon_dvs(lon_dv7, lon_dv8, lon_dv9, lon_dv10, lon_dv11)
  !write(*,*)'lon_dvs_check=: ',lon_dv7, lon_dv8, lon_dv9, lon_dv10, lon_dv11
! After this CALL, the range of variables are:
!    -180 =< gc_D_lon & gc_A_lon =< 180 [For all FL_DIRs]
!    -180 =< lon_dv7 to lon_dv11 =< 180 [FL_DIR=0,2]
!       0 =< lon_dv7 to lon_dv11 =< 360 [FL_DIR=1,3]
!     -90 =< lon_dv7 to lon_dv11 =< +90 [FL_DIR=4,5](<- those variables show latitude values)

!------------------------------------------
!Population loop for one generation(START)
!------------------------------------------
  DO ipop=1,num_pop  
     wp_const        = 0        !HY20140428-3
     w_zalt          = 0.0_dp   !HY20140430-1
     wind_along_traj = 0.0_dp   !HY20140430-1
     vtas_along_traj = 0.0_dp   !HY20140522-4
     rho_along_traj  = 0.0_dp   !HY20160531-1
     press_along_traj= 0.0_dp   !HY20170215-5
     t_scb_along_traj= 0.0_dp   !HY20170215-5
     !Yin_20170423
     cpc_along_traj  = 0.0_dp
     ATR20O3_along_traj  = 0.0_dp
     ATR20CH4_along_traj  = 0.0_dp
     ATR20H2O_along_traj  = 0.0_dp
     ATR20CPC_along_traj  = 0.0_dp
     ATR20CO2_along_traj  = 0.0_dp
     !Yin_20170423

  !Input DVs
  !HY20140327-4,20140328-4  open(10,file='mot.input')
     DO idv=1,num_dvs  !HY20140415-5
        DV(idv)=dvs(nsum1)  !HY20140415-6
        nsum1=nsum1+1
     ENDDO
! Ranges of DVs:
!    -180 =< DV(1),DV(3),DV(5) =< 180 [FL_DIR=0,2,4,5]
!       0 =< DV(1),DV(3),DV(5) =< 360 [FL_DIR=1,3]



  !IMPORTANT!!!20131216-3
  ! SOUND_V (i.e. VTAS) should be calculated at each waypoints during optimization.
  ! the following "ALT" is changed by each altitude!!!!!! delete USE -, ONlY:ALT from gc_option above.

  !Calculate SOUND_V and BADA_FUEL_NOM at ALT by using BADA_DATA array.
  !maybe this calculation is performed after generating waypoints for each trajectory......instead of using "ALT"! 
  
  !see 20140524-9
  !SOUND_V,[m/s]
!  do i=1,7
!     if((BADA_DATA(i,4)<=ALT).and.(ALT<BADA_DATA(i+1,4)))then
!        SOUND_V=(BADA_DATA(i,2)*(ALT-BADA_DATA(i+1,4))-BADA_DATA(i+1,2)*(ALT-BADA_DATA(i,4)))/ &
!                (BADA_DATA(i,4)-BADA_DATA(i+1,4))
!        write(*,*)"if_soundv_1",i,ALT,BADA_DATA(i+1,4)
!        exit
!     elseif(ALT.lt.BADA_DATA(1,4))then
!        SOUND_V=BADA_DATA(1,2)
!        write(*,*)"if_soundv_2",i,ALT,BADA_DATA(1,4)
!        exit
!     elseif(ALT.ge.BADA_DATA(8,4))then
!        SOUND_V=BADA_DATA(8,2)
!        write(*,*)"if_soundv_3",i,ALT,BADA_DATA(8,4)
!        exit  
!     endif
!  enddo
!  write(*,*)'SOUND_V_check',SOUND_V
  
  !This calculation is moved to messy_airtraf.f90.  see20140524-7 
  !BADA_FUEL_NOM,[kg/s]
!  do i=1,7
!     if((BADA_DATA(i,4)<=ALT).and.(ALT<BADA_DATA(i+1,4)))then
!        BADA_FUEL_NOM=(BADA_DATA(i,5)*(ALT-BADA_DATA(i+1,4))-BADA_DATA(i+1,5)*(ALT-BADA_DATA(i,4)))/ &
!                (BADA_DATA(i,4)-BADA_DATA(i+1,4))
!        write(*,*)"if_bada_fuel_1",i,ALT,BADA_DATA(i+1,4)
!        exit
!     elseif(ALT.lt.BADA_DATA(1,4))then
!        BADA_FUEL_NOM=BADA_DATA(1,5)
!        write(*,*)"if_bada_fuel_2",i,ALT,BADA_DATA(1,4)
!        exit
!     elseif(ALT.ge.BADA_DATA(8,4))then
!        BADA_FUEL_NOM=BADA_DATA(8,5)
!        write(*,*)"if_bada_fuel_3",i,ALT,BADA_DATA(8,4)
!        exit  
!     endif
!  enddo
!  write(*,*)'BADA_FUEL_NOM_check',BADA_FUEL_NOM

  !Definition of D_city and A_city
  !Flightplan information from SMCL 
!  gc_D_lon  = w_p2_ac_routes_desc(D_lon)      !lon of departure-city   
!  gc_D_lat  = w_p2_ac_routes_desc(D_lat)      !lat of departure-city   
!  gc_D_time = w_p2_ac_routes_desc(D_time)     !departure time in Julian date 
!  gc_A_lon  = w_p2_ac_routes_desc(A_lon)      !lon of arrival-city
!  gc_A_lat  = w_p2_ac_routes_desc(A_lat)      !lat of arrival-city  
!  temp_gc_A_lon = gc_A_lon + 360.0_dp
!  temp_gc_D_lon = gc_D_lon + 360.0_dp

!  OPEN(40,file="north_check.dat",position='APPEND')  !20131203-3
!     if(gc_D_lat.lt.0.0_dp)write(40,*)"1",gc_D_lat,gc_D_lon,gc_A_lat,gc_A_lon
!     if(gc_A_lat.lt.0.0_dp)write(40,*)"2",gc_D_lat,gc_D_lon,gc_A_lat,gc_A_lon
!     if(gc_A_lat.lt.0.0_dp.and.gc_D_lat.lt.0.0_dp)write(40,*)"3",gc_D_lat,gc_D_lon,gc_A_lat,gc_A_lon
!  CLOSE(40)
  !write(*,*)'SUBSUB_check',gc_D_lon,gc_D_lat,gc_D_time,gc_A_lon,gc_A_lat
  !write(*,*)'SUBSUB_check',D_lon,D_lat,D_time,A_lon,A_lat
  !write(*,*)'SUBSUB_check',PI, DTR,OneDay,radius_earth

  !Initialization: Obj-function, 20170609-2
  T_FLIGHT_TIME  = obj_dummy   !=0.0_dp
  Total_fuel_use = obj_dummy   !=0.0_dp,20160606-3
  Total_nox_emis = obj_dummy   !=0.0_dp,20170215-6
  Total_h2o_emis = obj_dummy   !=0.0_dp,20170222-3
  !Yin_20170423
  T_CPC          = obj_dummy   !=0.0_dp,20170607-2
  Total_atr20    = obj_dummy   !=0.0_dp,20170801
  Total_costclim = obj_dummy   !=0.0_dp,20170801
  Total_costcpc = obj_dummy   !=0.0_dp,20170801
  !Yin_20170423
  Total_cost     = obj_dummy   !=0.0_dp,20170224-3
  Total_coc      = obj_dummy   !=0.0_dp,20170316-3
!HY  AC_V   = AC_MACH*SOUND_V/1000.0_dp    !AC speed [km/s], 20130719-1,20131031-1,-2,20131104-1,20131114-2,20130315-2
                                           !20140522-4  

  !Definition of the city being placed at higher lon. 
!  IF(gc_D_lon.lt.gc_A_lon)THEN   !20111219-1
!     RARAD1= gc_A_lon*DTR        !lon of higher city [rad]
!     RARAD2= gc_D_lon*DTR        !lon of lower  city [rad]
!     DCRAD1= gc_A_lat*DTR        !lat of higher city [rad]    
!     DCRAD2= gc_D_lat*DTR        !lat of lower  city [rad]
!     IF((RARAD1-RARAD2).le.PI)THEN
!        FL_DIR= 0
!     ELSE
!        FL_DIR= 1
!     ENDIF
!  ELSE IF(gc_D_lon.gt.gc_A_lon)THEN
!     RARAD1= gc_D_lon*DTR     
!     RARAD2= gc_A_lon*DTR
!     DCRAD1= gc_D_lat*DTR     
!     DCRAD2= gc_A_lat*DTR     
!     IF((RARAD1-RARAD2).le.PI)THEN 
!        FL_DIR= 2
!     ELSE
!        FL_DIR= 3
!     ENDIF
!  ELSE     !(gc_D_lon.eq.gc_A_lon)
!     IF(gc_D_lat.gt.gc_A_lat)THEN
!        RARAD1= gc_D_lon*DTR     
!        RARAD2= gc_A_lon*DTR
!        DCRAD1= gc_D_lat*DTR     
!        DCRAD2= gc_A_lat*DTR     
!        FL_DIR= 4
!     ELSE IF(gc_D_lat.lt.gc_A_lat)THEN
!        RARAD1= gc_A_lon*DTR 
!        RARAD2= gc_D_lon*DTR
!        DCRAD1= gc_A_lat*DTR     
!        DCRAD2= gc_D_lat*DTR 
!        FL_DIR= 5
!     ELSE     !(gc_D_lat.eq.gc_A_lat)
!        write(*,*)'City Selection Error: They must be identical city!'
!     ENDIF
!  ENDIF

  !opt20130205-3
!  CALL determine_lon_dvs(lon_dv7, lon_dv8, lon_dv9, lon_dv10, lon_dv11)
  !write(*,*)'lon_dvs_check=: ',lon_dv7, lon_dv8, lon_dv9, lon_dv10, lon_dv11
! After this CALL, the range of variables are:
!    -180 =< gc_D_lon & gc_A_lon =< 180 [For all FL_DIRs]
!    -180 =< lon_dv7 to lon_dv11 =< 180 [FL_DIR=0,2]
!       0 =< lon_dv7 to lon_dv11 =< 360 [FL_DIR=1,3]
!     -90 =< lon_dv7 to lon_dv11 =< +90 [FL_DIR=4,5](<- those variables show latitude values)

  !Calculation of GC distance
  !HY 20120810
  !CALL gc_distance(COSIS, T_DIS, T_DISN)

  CALL arbitrary_traj(XXN, YYN, ZZN)
! After this CALL, waypoints are detemined. The range of variables are:see opt20130215-4,20140512-4
!    -180 =< XXN(LON) =< 180 [FL_DIR=0,2,4,5]
!       0 =< XXN(LON) =< 360 [FL_DIR=1,3]
! write(*,*)'CHECK_(LON, LAT, ALT) TRAJ'
!  do i=1, DIV_WAY+1
!     write(*,*)XXN(i),YYN(i),ZZN(i)
!  enddo
 
  !Constraints check on waypoints  HY20140428-3
  DO i=1, DIV_WAY+1
     !write(*,*)"ZZN_Lower_alt",ZZN(i),Lower_alt,"Upper_alt_ZZN",Upper_alt,ZZN(i)
     IF(ZZN(i).lt.Lower_alt .OR. Upper_alt.lt.ZZN(i))THEN
        wp_const = 1
        exit
     ENDIF
  ENDDO

  !Calculation of waypoints on GC
  !CALL gc_waypoints(THETA_OUT, PHI_OUT)


  !wp_const=0: waypoints generation well done -> extract wind-values by nn_index for each wp. 
  !                                           -> calculate obj value by CALL gc_ac_traj. 
  !wp_const=1: waypoints contain odd points.  -> dummy_obj is given as obj value.
  !
  !Detailed description of the following part: 20170609-2,3
  IF(wp_const.eq.0)THEN
     temp_XXN=0.0_dp     !20140605-1
     DO i=1,DIV_WAY+1    !?? DIV_WAY??
        !Extract wind values by nn_index. 
        !nn_index check  20140212-2,20130905-3
        !nn_index can find nearest index.
        
        !If [-180=< XXN(i) < 0], add +360 to XXN(i) and use temp_XXN for nn_index search.20140605-2
        IF((-180.0_dp.le.XXN(i)).and.(XXN(i).lt.0.0_dp))THEN
           temp_XXN = XXN(i)+360.0_dp
           CALL nn_index(w_philon(:), temp_XXN, w_j_lon) !lon value on waypoints,20140605-1
        ELSE
           temp_XXN = XXN(i)
           CALL nn_index(w_philon(:), temp_XXN, w_j_lon)   !lon value on waypoints
        ENDIF
       !CALL nn_index(w_philon(:), XXN(i), w_j_lon)  !lon value on waypoints
        CALL nn_index(w_philat(:), YYN(i), w_i_lat)  !lat value on waypoints
        w_zalt(:)=w_zgl_geopot_3d(w_j_lon, :, w_i_lat)/g   ![m] geopotential altitude,20131105-1,20140430-1
        !write(*,*)'(wind)g=',g,'(wind)nlev_w_zalt=',nlev_w_zalt
        !write(*,*)'(wind)w_zgl_geopot_3d=',w_zgl_geopot_3d(w_j_lon,:,w_i_lat)/g
        !write(*,*)'(wind)w_zalt=',w_zalt(:)
        CALL nn_index(w_zalt(:), ZZN(i), w_i_alt)  
        !write(*,*)'(wind)alt=', ZZN(i)
        !write(*,*)"Nearest indices=", w_i_lat, w_j_lon, w_i_alt
        !write(*,*)"nn_index_lat_check", w_i_lat, 'YYN(i)=',YYN(i)   !20140604-1
        !write(*,*)"nn_index_lon_check", w_j_lon, 'XXN(i)=',XXN(i),'temp_XXN=',temp_XXN   !20140604-1,20140605-2
        !write(*,*)"nn_index_alt_check", w_i_alt, 'ZZN(i)=',ZZN(i)   !20140604-1
        !Extract wind values from global wind fields on each wp by using resulting indices.20130813-2,20130905-3
        !If sound speed is also needed, see20140212-3.
        !Reserve wind values along waypoints; they are passed to CALL gc_ac_traj.
        !wind_along_traj(:,1) :u [m/s] ,20140430-1
        !               (:,2) :v [m/s]
        !               (:,3) :w [m/s]
        !vtas_along_traj(:)   :AC_MACH*speed of sound [m/s]
        wind_along_traj(i,1) = w_uwind_g(w_j_lon, w_i_alt, w_i_lat)    ![m/s]
        wind_along_traj(i,2) = w_vwind_g(w_j_lon, w_i_alt, w_i_lat)    ![m/s]
        wind_along_traj(i,3) =   w_v_z_g(w_j_lon, w_i_alt, w_i_lat)    !global vertical wind velocity field [m/s]
        !Yin_20170423
        cpc_along_traj(i)    =   w_cpc_g(w_j_lon, w_i_alt, w_i_lat)    !global potcov field [-]
        !Yin_20170423
       ! write(*,*)'ATR20CPC_along_traj:',SIZE(ATR20CPC_along_traj)
       ! write(*,*)'ATR20CPC_along_traj_check',i,ATR20CPC_along_traj(i)  
!       write(*,*)'wind_cpc_along_traj:',SIZE(cpc_along_traj)
!       write(*,*)'wind_cpc_along_traj_check',i,cpc_along_traj(i)
        !       vtas_along_traj(i)   = AC_MACH*SQRT(w_t_scb_g(w_j_lon, w_i_alt, w_i_lat)*1.4_dp*287.05287_dp)  ![m/s] 20140522-4 
        vtas_along_traj(i)   = AC_MACH*SQRT(w_t_scb_g(w_j_lon, w_i_alt, w_i_lat)*   &
                               adiabatic_index_air*gas_constant_air)   ![m/s] 20140724-3
        !write(*,*)"wind_along_traj=",(wind_along_traj(i,j),j=1,3) 
        !write(*,*)"vtas_along_traj=",vtas_along_traj(i),"[m/s]" 
        if(wind_along_traj(i,3).gt.1.0_dp)write(*,*)"Check_this_case"  !20140429-2

        !Extract rho values from global rho fields on each wp by using resulting indices.
        !Reserve rho values along waypoints.
        IF(w_option_traj_calc_g==2 .or. w_option_traj_calc_g==4 .or.  &   !Fuel_opt=2,20170221-1,H2O_opt=4,20170222-3
           w_option_traj_calc_g==6 .or. w_option_traj_calc_g==7.or.&
           w_option_traj_calc_g==10)THEN      !Cost_opt=6,20170224-2,COC_opt=7,20170316-2
           rho_along_traj(i) = w_rho_air_dry_3d_g(w_j_lon, w_i_alt, w_i_lat)    !20170607-3
!          write(*,*)"rho_along_traj_1=",rho_along_traj(i),"[kg/m^3]"  !20170223-4              
        ELSEIF(w_option_traj_calc_g==3)THEN   !NOx_opt=3 20170215-5,20170221-1   
           rho_along_traj(i) = w_rho_air_dry_3d_g(w_j_lon, w_i_alt, w_i_lat)    !20170607-3
           press_along_traj(i) = w_press_3d_g(w_j_lon, w_i_alt, w_i_lat)        !20170607-3 
           t_scb_along_traj(i) = w_t_scb_g(w_j_lon, w_i_alt, w_i_lat)
        !Yin_20170801+
        ELSEIF(w_option_traj_calc_g==8.or.w_option_traj_calc_g==9)THEN   !ATR20_opt=8 , COSTCLIM_opt=9,20170801   
           rho_along_traj(i) = w_rho_air_dry_3d_g(w_j_lon, w_i_alt, w_i_lat)    !20170607-3
           press_along_traj(i) = w_press_3d_g(w_j_lon, w_i_alt, w_i_lat)        !20170607-3 
           t_scb_along_traj(i) = w_t_scb_g(w_j_lon, w_i_alt, w_i_lat) 
           ATR20O3_along_traj(i)=   w_ATR20O3_3d_g(w_j_lon, w_i_alt, w_i_lat)    !global ATR20O3 field [K/kg(NO2)]
           ATR20CH4_along_traj(i)    =   w_ATR20CH4_3d_g(w_j_lon, w_i_alt, w_i_lat)    !global ATR20CH4 field [K/kg(NO2)]
           ATR20H2O_along_traj(i)    =   w_ATR20H2O_3d_g(w_j_lon, w_i_alt, w_i_lat)    !global ATR20H2O field [K/kg(fuel)]
           ATR20CPC_along_traj(i)    =   w_ATR20CPC_3d_g(w_j_lon, w_i_alt, w_i_lat)    !global ATR20CPC field [K/km(contrail)]
           ATR20CO2_along_traj(i)    =   w_ATR20CO2_3d_g(w_j_lon, w_i_alt, w_i_lat)    !global ATR20CO2 field [K/kg(fuel)]
          ! write(*,*)'ATR20O3_along_traj:',SIZE(ATR20O3_along_traj)
          ! write(*,*)'ATR20CH4_along_traj:',SIZE(ATR20CH4_along_traj)
          ! write(*,*)'ATR20H2O_along_traj:',SIZE(ATR20H2O_along_traj)
          ! write(*,*)'ATR20O3_along_traj_check',i,ATR20O3_along_traj(i)
          ! write(*,*)'ATR20CH4_along_traj_check',i,ATR20CH4_along_traj(i)
          ! write(*,*)'ATR20H2O_along_traj_check',i,ATR20H2O_along_traj(i) 
          ! write(*,*)'ATR20CO2_along_traj_check',i,ATR20CO2_along_traj(i) 
        !Yin_20170801-
        ENDIF
     ENDDO
     
!     IF(w_option_traj_calc_g==3)write(*,*)'w_press_3d_g[Pa]=',&
     !press_3d_g(w_j_lon, w_i_alt, w_i_lat)  !20170215-4 [Pa],20170607-3

     !Calculation of Obj (time_opt) and ac_routes properties (time_opt, fuel_opt, NOx_opt, H2O_opt, Cost_opt, COC_opt).
     CALL gc_ac_traj(XXN, YYN, ZZN, cpc_along_traj, wind_along_traj, vtas_along_traj,    &   !vtas_along_traj 20140522-4
                     FLIGHT_TIME, FLIGHT_TIME_JULIAN, T_FLIGHT_TIME, &   !20170607-1
                     A_time, FLIGHT_DISTANCE, T_FLIGHT_DISTANCE, FLIGHT_SPEED, contrail_coverage, T_CPC, ATR20CPC_along_traj, &
                     Total_atr20cpc)   !FLIGHT_SPEED = ac_gs[m/s]
     !After this CALL, the range of variables are: see opt 20130215-5

     !Calculation of obj functions: Total_fuel_use (fuel_opt) by Total Energy model 20160531-2
     !                              Total_nox_emis (NOx_opt) by DLR Fuel Flow method 20170215-6 
     !                              Total_h2o_emis (H2O_opt) 20170222-3 
     !                              Total_cost (Cost_opt) 20170224-2
     !                              Total_coc (COC_opt) by NASA Liebeck model 20170316-3
     SELECT CASE(w_option_traj_calc_g)   !20170221-1
     CASE(2,4,6,7,10) !Fuel_opt=2, H2O_opt=4, Cost_opt=6, COC_opt=7  !20170222-3,20170224-3,20170316-3 
        CALL fuel_calc(ZZN(:),T_FLIGHT_TIME, vtas_along_traj(:),FLIGHT_TIME_JULIAN(:),   &
                       rho_along_traj(:), Total_fuel_use)   !20170216-5,fuel_opt,H2O_opt,Cost_opt,COC_opt Total_nox_emis=dummy 
!       Total_h2o_emis = Total_fuel_use*1230.0_dp  ! [g(H2O)/kg(fuel)] 20170222-3
!       Cf = Fuel_price*100.0_dp/Fuel_density      !Unit fuel costs[Cents/lbs] 20170224-5,20170309-3,20170428-4        
!       Total_cost = (Cf*0.01_dp)*Total_fuel_use*2.2046226218_dp + Ct*(T_FLIGHT_TIME/3600.0_dp)+Co ![US Dollar]20170224-3,20170224-5

        !H2O emission calc. for H2O_opt=4, 20170428-5
        IF(w_option_traj_calc_g==4)THEN
           Total_h2o_emis = Total_fuel_use*1230.0_dp                  ![g(H2O)/kg(fuel)] 20170222-3, 20170428-5 
!           write(*,*)"CHECK_H2O:",Total_h2o_emis,"[g(H2O)/kg(fuel)] routing_option:",w_option_traj_calc_g
        ENDIF 
 
        !Cost calc. for Cost_opt=6, 20170428-4
        IF(w_option_traj_calc_g==6)THEN
           CALL cost_calc(T_FLIGHT_TIME, Total_fuel_use, Total_cost) 
!           write(*,*)"CHECK_COST:",Total_cost,"[US Dollar] routing_option:",w_option_traj_calc_g
        ENDIF 

        !COC calc. for COC_opt=7, 20170316-7
        IF(w_option_traj_calc_g==7)THEN
           CALL coc_calc(T_FLIGHT_TIME, Total_fuel_use, Total_coc) 
!           write(*,*)"CHECK_COC:",Total_coc,"[US Dollar] routing_option:",w_option_traj_calc_g
        ENDIF 
        IF(w_option_traj_calc_g==10)THEN
           CALL contrail_cost_calc(T_FLIGHT_TIME,T_CPC, Total_fuel_use, Total_costcpc) 
!           write(*,*)"CHECK_COC:",Total_costcpc,"[US Dollar] routing_option:",w_option_traj_calc_g
        ENDIF 

!       write(*,*)"rho_along_traj_2=",rho_along_traj,"[kg/m^3]"   !20170223-4
!       write(*,*)"Total_fuel_use(Obj)=",Total_fuel_use ,'w_nwaypoints=',w_nwaypoints   ![kg] 
!       write(*,*)"Total_nox_emis(Obj)=",Total_nox_emis,'w_nwaypoints=',w_nwaypoints    ![g(NO2)] 
!       write(*,*)"Total_h2o_emis(Obj)=",Total_h2o_emis,'w_nwaypoints=',w_nwaypoints    ![g(H2O)] 
!       write(*,*)"Total_cost(Obj)=",Total_cost,"[US Dollar], Cf=",Cf,"[Cents/Pound]"   ![US Dollar]20170309-3
!       write(*,*)"Total_coc(Obj)=",Total_coc,"[US Dollar]"   ![US Dollar]20170428-4
!       write(*,*)"ZZN",ZZN
!       write(*,*)"T_FLIGHT_TIME",T_FLIGHT_TIME
!       write(*,*)"vtas_along_traj",vtas_along_traj
!       write(*,*)"FLIGHT_TIME_JULIAN",FLIGHT_TIME_JULIAN
     CASE(3,8,9) !NOx_opt=3 , ATR20_opt=8,COSTCLIM_opt=9
        CALL fuel_calc(ZZN(:),T_FLIGHT_TIME, vtas_along_traj(:),FLIGHT_TIME_JULIAN(:),   &
                       rho_along_traj(:), Total_fuel_use,                                & !20160603-1
                       Total_nox_emis,&
                       press_along_traj(:), t_scb_along_traj(:),ATR20O3_along_traj(:),   &
                       ATR20CH4_along_traj(:),ATR20H2O_along_traj(:),ATR20CO2_along_traj(:),&
                       Total_atr20o3,&
                       Total_atr20ch4,Total_atr20h2o,Total_atr20co2)           !20170215-6,20170216-5 
!        write(*,*)"Total_fuel_use(Obj)=",Total_fuel_use,'w_nwaypoints=',w_nwaypoints    ![kg] 
!        write(*,*)"Total_nox_emis(Obj)=",Total_nox_emis,'w_nwaypoints=',w_nwaypoints    ![g(NO2)] 
!       write(*,*)"Total_h2o_emis(Obj)=",Total_h2o_emis,'w_nwaypoints=',w_nwaypoints    ![g(H2O)] 
!       write(*,*)"Total_cost(Obj)=",Total_cost,"[US Dollar], Cf=",Cf,"[Cents/Pound]"   ![US Dollar]20170309-3,20170428-5
!       write(*,*)"Total_coc(Obj)=",Total_coc,"[US Dollar]"   ![US Dollar]201703
        IF(w_option_traj_calc_g==8) THEN !ATR20_opt=8
             CALL atr20_calc(Total_atr20o3,Total_atr20ch4,Total_atr20h2o,   &
                             Total_atr20cpc,Total_atr20co2,Total_atr20)           !20170801 
                                                                                                   !20170309-3,20170428-5 
!             write(*,*)"Total_atr20(Obj)=",Total_atr20,'w_nwaypoints=',w_nwaypoints    ![K] 
        ELSEIF(w_option_traj_calc_g==9) THEN !ATR20_opt=8,COSTCLIM_opt=9
             CALL cost_atr20_calc(T_FLIGHT_TIME,Total_fuel_use,&
                                  Total_atr20o3,Total_atr20ch4,Total_atr20h2o,   &
                                  Total_atr20cpc,Total_atr20co2,Total_costclim)           !20170801 
                                                                                                   !20170309-3,20170428-5 
!             write(*,*)"Total_costclim(Obj)=",Total_costclim,'w_nwaypoints=',w_nwaypoints    ![K] 
        ENDIF
     CASE DEFAULT !Time_opt=1, CPC_opt=5,20170609-2,3
!       write(*,*)'CHECK:w_option_traj_calc_g', w_option_traj_calc_g                   !20170221-1 
     ENDSELECT 
  ELSEIF(wp_const.eq.1)THEN
     !Add obj_dummy. 20170217-5
     T_FLIGHT_TIME  = obj_dummy  !20140428-3
     Total_fuel_use = obj_dummy  !20160606-3 
     Total_nox_emis = obj_dummy  !20170215-6
     Total_h2o_emis = obj_dummy  !20170222-3
     !Yin_20170423
     T_CPC          = obj_dummy  !20170609-2,3 
     Total_atr20    = obj_dummy  !20170801 
     Total_costclim = obj_dummy  !20170801 
     Total_costcpc  = obj_dummy  !20170801 
     !Yin_20170423
     Total_cost     = obj_dummy  !20170224-3
     Total_coc      = obj_dummy  !20170316-3
     write(*,*)"wp_const_error",wp_const
  ENDIF

  !Calculation of fuel_use for each segment using Total Energy model 
  !CALL total_energy_model_calculation(FLIGHT_TIME, FUEL_USE) 

!HY  total_flight_time_s = T_FLIGHT_TIME   !20131003-3,20131114-2,20140524-4
!HY  write(*,*)'total_flight_time_s=',total_flight_time_s 

  !Outputs of flightplan information
!  write(*,'(x,2a,2x,f10.4,2x,a,f10.4,2x,a,f8.1)')'Departure city     :','lat',gc_D_lat,'lon',gc_D_lon,'Alt',ALT
!  write(*,'(x,2a,2x,f10.4,2x,a,f10.4,2x,a,f8.1)')'Arrival city       :','lat',gc_A_lat,'lon',gc_A_lon,'Alt',ALT
!  write(*,*)'Departure time[Julian date]  :',gc_D_time
!  write(*,*)'Arrival time[Julian date]    :',A_time
!  write(*,*)'Total flight time    :',T_FLIGHT_TIME,'sec'
!  write(*,'(x,a,f9.4,a,f6.4,a)')'Flight Speed(Mach)   :',AC_V*3600.0_dp,'km/h(',AC_MACH,')'
!  write(*,*)'ISWITCH              :',ISWITCH
!  write(*,*)'OUT_ISWITCH          :',OUT_SWITCH
!  write(*,*)'Flight direction     :',FL_DIR           !FL_DIR: see 20111219-1
!! write(*,'(x,a,f7.4,x,a)')'Central angle between D_city and A-city:',COSIS,'radians'
!! write(*,'(x,a,2(f10.4,a))')'Flight distance (GC) between D_city and A-city:',T_DIS,'km;',T_DISN,'nm'
!  write(*,'(x,a,2(f10.4,a))')'Flight distance (GC) between D_city and A-city:',T_FLIGHT_DISTANCE,'km;'&
!                             ,T_FLIGHT_DISTANCE*0.539957_dp,'nm'
!  write(*,*)'Flight distance (wind) between D_city and A-city:',T_FLIGHT_DISTANCE,'km;'
!  write(*,'(/,a)')'>>>>>>>>>>>>>>>>> Flight PLAN on GC >>>>>>>>>>>>>>>>>>>'
!  write(*,'(a,3x,a,3(5x,a))')'Time[Julian]','Time [sec]','Lat','Lon','Alt'

  !Arrange the outputs range
!  SELECT CASE(OUT_SWITCH)
!  CASE(0)  
     ! Only for FL_DIR = 1,3, the following DO structure will be adapted. 
     ! FL_DIR of (MUC to JFK) is 'FL_DIR = 2'.
!     DO i=1,DIV_WAY+1
!        IF(PHI_OUT(i).gt.180.0_dp)PHI_OUT(i) = PHI_OUT(i) - 360.0_dp    !opt 20130215-5
!        IF(XXN(i).gt.180.0d0)XXN(i) = XXN(i) - 360.0d0
!     ENDDO
!     DO j=1, DIV_WAY+1
!        write(*,'(f15.6,2x,i6,2x,f10.4,x,f10.4,2x,f7.1)')&
!        write(*,'(f15.6,2x,f15.6,2x,f10.4,x,f10.4,2x,f7.1)')&
!               FLIGHT_TIME_JULIAN(j),INT(FLIGHT_TIME(j)),THETA_OUT(j),PHI_OUT(j),ALT
!               FLIGHT_TIME_JULIAN(j),FLIGHT_TIME(j),YYN(j),XXN(j),ZZN(j)
!     ENDDO

!  CASE(1)
!     DO j=1, DIV_WAY+1
!        write(*,'(f15.6,2x,i6,2x,f10.4,x,f10.4,2x,f7.1)')&
!               FLIGHT_TIME_JULIAN(j),INT(FLIGHT_TIME(j)),THETA_OUT(j),PHI_OUT(j),ALT
!HY               FLIGHT_TIME_JULIAN(j),FLIGHT_TIME(j),YYN(j),XXN(j),ZZN(j)   
!     ENDDO

!  CASE DEFAULT
!     write(*,*) 'OUT_SWITCH error: Unknown OUT_SWITCH on GC output!'
!  ENDSELECT

  !Outputs array for upper SMCL(ac_routes properties): wind optimal trajectory info.
  !This should be done in Subroutine armoga of messy_airtraf_tools_ga.f90.
!  DO j=1, DIV_WAY+1
!     gc_ac_routes(j,props_lon)           = PHI_OUT(j)             !lon
!     gc_ac_routes(j,props_lat)           = THETA_OUT(j)           !lat
!     gc_ac_routes(j,props_alt)           = ALT                    !alt [m]
!     gc_ac_routes(j,props_time)          = FLIGHT_TIME_JULIAN(j)  !Passing time at each waypoint, [Julian date]
!     gc_ac_routes(j,props_ac_speed)      = AC_V*3600.0_dp         !Aircraft_speed, [km/h] !should be FLIGHT_SPEED array!!!! 
!     gc_ac_routes(j,props_dist)          = FLIGHT_DISTANCE(j)     !Distance for each segment, [km]  
!!    gc_ac_routes(j,props_fuel_use)      = FUEL_USE(j)            !Fuel consumption for each segment, [kg]  
!!    gc_ac_routes(j,props_emis_nox)      = 1.0_dp                 !Dummy for nox_emis 
!!    gc_ac_routes(j,props_emis_h2o)      = 1.0_dp                 !Dummy for h2o_emis
!!    write(10,'(f10.4,x,f10.4,2x,2(8x,a))')PHI_OUT(j),THETA_OUT(j),'1.000','0.01'
!  ENDDO

  !Outputs for GNUPLOT.
!  OPEN(10,file='gnuplot.dat') 
!    DO i=1,DIV_WAY+1
!       IF(PHI_OUT(i).gt.180.0d0)PHI_OUT(i) = PHI_OUT(i) - 360.0d0
!    ENDDO
!     DO j=1, DIV_WAY+1
!       write(10,'(f7.4,x,f9.4,2x,f7.1)')THETA_OUT(j),PHI_OUT(j),ALT
!        write(10,'(f15.6,x,f10.4,x,f10.4,2x,2(8x,a))')FLIGHT_TIME_JULIAN(j),PHI_OUT(j),THETA_OUT(j),'1.000','0.01'
!!       write(10,'(f15.6,x,f10.4,x,f10.4,2x,2(8x,a))')FLIGHT_TIME(j),XXN(j),YYN(j),'1.000','0.01'   !
!        write(10,'(f15.6,x,f10.4,x,f10.4,x,f10.4,2x,2(8x,a))')FLIGHT_TIME(j),XXN(j),YYN(j),ZZN(j),'1.000','0.01'   !
!     ENDDO
!  CLOSE(10)

  !Outputs for GMT. see 20130118-1
!  OPEN(20,file='gmt.dat',position='append')
!     DO j=1, DIV_WAY+1
!        write(20,*)PHI_OUT(j),THETA_OUT(j)
!     ENDDO
!     write(20,'(a)')">"
!  CLOSE(20)

  !Output for tecplot, 20140704-1,20141218-2,20150305-1,2
  !MUC >> JFK  
!  IF((INT(gc_D_lon).eq.11) .and. (INT(gc_D_lat).eq.48) .and. (INT(gc_A_lon).eq.-73) .and. (INT(gc_A_lat).eq.40))THEN
!!  IF((INT(gc_D_lon).eq.-83) .and. (INT(gc_D_lat).eq.42) .and. (INT(gc_A_lon).eq.8) .and. (INT(gc_A_lat).eq.50))THEN
!!    write(*,*)"output_tecplot_muc_jfk",INT(11.786),INT(48.35),INT(-73.77),INT(40.63)
!     i_add_nnindex = 1     !20150312-1
!     OPEN(20,file='tecplot_muc_jfk.dat',position='append')
!        DO j=1, DIV_WAY+1
!           write(20,*)XXN(j),YYN(j),ZZN(j)
!        ENDDO
!        write(20,*)
!        write(20,*)
!     CLOSE(20)
!     OPEN(22,file='tecplot_muc_jfk_gmt.dat',position='append')
!        DO j=1, DIV_WAY+1
!           write(22,*)XXN(j),YYN(j),ZZN(j)
!        ENDDO
!        write(22,'(a)')">" 
!     CLOSE(22)
!  ENDIF
!  !JFK >> MUC  
!  IF((INT(gc_D_lon).eq.-73) .and. (INT(gc_D_lat).eq.40) .and. (INT(gc_A_lon).eq.11) .and. (INT(gc_A_lat).eq.48))THEN
!!    write(*,*)"output_tecplot_jfk_muc",INT(11.786),INT(48.35),INT(-73.77),INT(40.63)
!     i_add_nnindex = 1     !20150312-1
!     OPEN(21,file='tecplot_jfk_muc.dat',position='append')
!        DO j=1, DIV_WAY+1
!           write(21,*)XXN(j),YYN(j),ZZN(j)
!        ENDDO
!        write(21,*)
!        write(21,*)
!     CLOSE(21)
!     OPEN(23,file='tecplot_jfk_muc_gmt.dat',position='append')
!        DO j=1, DIV_WAY+1
!           write(23,*)XXN(j),YYN(j),ZZN(j)
!        ENDDO
!        write(23,'(a)')">" 
!     CLOSE(23)
!  ENDIF
!  !EHAM >> KSEA
!  IF((INT(gc_D_lon).eq.4) .and. (INT(gc_D_lat).eq.52) .and. (INT(gc_A_lon).eq.-122) .and. (INT(gc_A_lat).eq.47))THEN
!     i_add_nnindex = 1     !20150312-1
!     OPEN(24,file='tecplot_eham_ksea.dat',position='append')
!        DO j=1, DIV_WAY+1
!           write(24,*)XXN(j),YYN(j),ZZN(j)
!        ENDDO
!        write(24,*)
!        write(24,*)
!     CLOSE(24)
!     OPEN(26,file='tecplot_eham_ksea_gmt.dat',position='append')
!        DO j=1, DIV_WAY+1
!           write(26,*)XXN(j),YYN(j),ZZN(j)
!        ENDDO
!        write(26,'(a)')">" 
!     CLOSE(26)
!  ENDIF
!  !KSEA >> EHAM
!  IF((INT(gc_D_lon).eq.-122) .and. (INT(gc_D_lat).eq.47) .and. (INT(gc_A_lon).eq.4) .and. (INT(gc_A_lat).eq.52))THEN
!     i_add_nnindex = 1     !20150312-1
!     OPEN(25,file='tecplot_ksea_eham.dat',position='append')
!        DO j=1, DIV_WAY+1
!           write(25,*)XXN(j),YYN(j),ZZN(j)
!        ENDDO
!        write(25,*)
!        write(25,*)
!     CLOSE(25)
!     OPEN(27,file='tecplot_ksea_eham_gmt.dat',position='append')
!        DO j=1, DIV_WAY+1
!           write(27,*)XXN(j),YYN(j),ZZN(j)
!        ENDDO
!        write(27,'(a)')">" 
!     CLOSE(27)
!  ENDIF
!  !EHAM >> KMSP
!  IF((INT(gc_D_lon).eq.4) .and. (INT(gc_D_lat).eq.52) .and. (INT(gc_A_lon).eq.-93) .and. (INT(gc_A_lat).eq.44) .and. &
!     (INT(gc_D_time).ge.2443510))THEN
!     i_add_nnindex = 1     !20150312-1
!     OPEN(28,file='tecplot_eham_kmsp.dat',position='append')
!        DO j=1, DIV_WAY+1
!           write(28,*)XXN(j),YYN(j),ZZN(j)
!        ENDDO
!        write(28,*)
!        write(28,*)
!     CLOSE(28)
!     OPEN(30,file='tecplot_eham_kmsp_gmt.dat',position='append')
!        DO j=1, DIV_WAY+1
!           write(30,*)XXN(j),YYN(j),ZZN(j)
!        ENDDO
!        write(30,'(a)')">" 
!     CLOSE(30)
!  ENDIF
!  !KMSP >> EHAM
!  IF((INT(gc_D_lon).eq.-93) .and. (INT(gc_D_lat).eq.44) .and. (INT(gc_A_lon).eq.4) .and. (INT(gc_A_lat).eq.52) .and. &
!     (INT(gc_D_time).ge.2443510))THEN
!     i_add_nnindex = 1     !20150312-1
!     OPEN(29,file='tecplot_kmsp_eham.dat',position='append')
!        DO j=1, DIV_WAY+1
!           write(29,*)XXN(j),YYN(j),ZZN(j)
!        ENDDO
!        write(29,*)
!        write(29,*)
!     CLOSE(29)
!     OPEN(31,file='tecplot_kmsp_eham_gmt.dat',position='append')
!        DO j=1, DIV_WAY+1
!           write(31,*)XXN(j),YYN(j),ZZN(j)
!        ENDDO
!        write(31,'(a)')">" 
!     CLOSE(31)
!  ENDIF

  !Flight_check,20131203-1,20141214-1
!  OPEN(30,file="flight_check.dat",position='APPEND') 
!     if(T_FLIGHT_DISTANCE.ge.0.0_dp .and. T_FLIGHT_DISTANCE.lt.1000.0_dp)write(30,*)"0"
!     if(T_FLIGHT_DISTANCE.ge.1000.0_dp .and. T_FLIGHT_DISTANCE.lt.2000.0_dp)write(30,*)"1"
!     if(T_FLIGHT_DISTANCE.ge.2000.0_dp .and. T_FLIGHT_DISTANCE.lt.3000.0_dp)write(30,*)"2"
!     if(T_FLIGHT_DISTANCE.ge.3000.0_dp .and. T_FLIGHT_DISTANCE.lt.4000.0_dp)write(30,*)"3"
!     if(T_FLIGHT_DISTANCE.ge.4000.0_dp .and. T_FLIGHT_DISTANCE.lt.5000.0_dp)write(30,*)"4"
!     if(T_FLIGHT_DISTANCE.ge.5000.0_dp .and. T_FLIGHT_DISTANCE.lt.6000.0_dp)write(30,*)"5"
!     if(T_FLIGHT_DISTANCE.ge.6000.0_dp .and. T_FLIGHT_DISTANCE.lt.7000.0_dp)write(30,*)"6"
!     if(T_FLIGHT_DISTANCE.ge.7000.0_dp .and. T_FLIGHT_DISTANCE.lt.8000.0_dp)write(30,*)"7"
!     if(T_FLIGHT_DISTANCE.ge.8000.0_dp .and. T_FLIGHT_DISTANCE.lt.9000.0_dp)write(30,*)"8"
!     if(T_FLIGHT_DISTANCE.ge.9000.0_dp)write(30,*)"9"
!  CLOSE(30)
!C test 
!  call cprogram
!  call system('sh /data/yama_hi/EMAC/AirTraf/messy_2.41/messy/smcl/ctext.sh')

  !Outputs for Objective functions
  !Objective function is flight time in [s](not [Julian]). see20140420-3
  !
!HY20140327-5  OPEN(20,file='mot.output1')
!HY20140328-4  OPEN(20,file='mot.output')
!HY20140328-4     write(20,*)T_FLIGHT_TIME  !Without penalty
!HY     write(20,*)T_FLIGHT_TIME  !Without penalty
     !write(*,*)"check_lower_alt,upper_alt",Lower_alt, Upper_alt, obj_dummy, wp_const 
     DO iof=1,num_obj
!HY20140415-5        obj_val(1)=T_FLIGHT_TIME  !maybe obj_val(1)=REAL(T_FLIGHT_TIME)
!HY20140415-5        obj_val(2)=T_FLIGHT_TIME  !20140328-4
        SELECT CASE(w_option_traj_calc_g)      !20160606-3,20170221-1
        CASE(1)   !Time_opt=1
           obj_val(nsum2)=T_FLIGHT_TIME  !HY20140415-5 [s]
        CASE(2)   !Fuel_opt=2
           obj_val(nsum2)=Total_fuel_use  !20160606-3 [kg(fuel)]
        CASE(3)   !NOx_opt=3
           obj_val(nsum2)=Total_nox_emis  !20170216-1 [g(NO2)]
        CASE(4)   !H2O_opt=4
           obj_val(nsum2)=Total_h2o_emis  !20170222-4 [g(H2O)]
       !Yin_20170423
        CASE(5)   !Contrail_opt=5
           obj_val(nsum2)=T_CPC           !20170608-1 [km]
       !Yin_20170423
        CASE(6)   !Cost_opt=6
           obj_val(nsum2)=Total_cost      !20170224-3 [US Dollar]
        CASE(7)   !COC_opt=7
           obj_val(nsum2)=Total_coc       !20170316-3 [US Dollar]
       !Yin_20170801+
        CASE(8) 
           obj_val(nsum2)= Total_atr20 !20170801
        CASE(9)
           obj_val(nsum2)= Total_costclim 
        CASE(10)
           obj_val(nsum2)= Total_costcpc 
      
       !Yin_20170801-
        CASE DEFAULT
           obj_val(nsum2)=obj_dummy       !20160606-3 [N/A]
           write(*,*)"Error:obj_value"
        ENDSELECT
        nsum2=nsum2+1
     ENDDO
!    write(20,*)T_FLIGHT_TIME + (10**6)*(max(0.0d0,T_FLIGHT_TIME-23838.75493755843d0*1.1d0))  !Penalty term included
!    write(20,*)abs(23839.39-T_FLIGHT_TIME)  !'REAL' for ARMOGA format.
!HY20140328-4  CLOSE(20)
!HY20140327-5  OPEN(30,file='mot.output2')
!HY     write(30,*)T_FLIGHT_TIME  !Without penalty
!    write(30,*)T_FLIGHT_TIME + (10**6)*(max(0.0d0,T_FLIGHT_TIME-23838.75493755843d0*1.1d0))  !penalty term included 
!    write(30,*)abs(23839.39-T_FLIGHT_TIME)  !'REAL' for ARMOGA format.
!HY  CLOSE(30)

      !Find temporary best solution through the current generation.        
      !Yin_20170423
      !- Program for contrail potential coverage (Contrail_opt=5) is added
      !- The contrail potential coverage is collected for output
      !Yin_20170423
      SELECT CASE(w_option_traj_calc_g)  !20160606-4,20170221-1
      CASE(1)   !Time_opt=1 
         IF(T_FLIGHT_TIME.lt.temp_best_fl_time)THEN
            temp_best_fl_time = T_FLIGHT_TIME
            temp_best_gen     = num_gen
            !Yin_20170423
            temp_best_potcov  = T_CPC                                  !Total CPC [km] 20170607-2
            !Yin_20170423
            !Outputs for upper SMCL(ac_routes properties)
            DO i=1,DIV_WAY+1
               temp_best_traj(i,props_lon)     =XXN(i)                 !Lon
               temp_best_traj(i,props_lat)     =YYN(i)                 !Lat
               temp_best_traj(i,props_alt)     =ZZN(i)                 !Alt [m]
               temp_best_traj(i,props_time)    =FLIGHT_TIME_JULIAN(i)  !Passing time at each waypoint, [Julian date]
               temp_best_traj(i,props_ac_speed)=FLIGHT_SPEED(i)*3.6_dp !AC_ground_speed[km/h],20140423-2,FLIGHT_SPEED[m/s]20140523-2 
               temp_best_traj(i,props_dist)    =FLIGHT_DISTANCE(i)     !Distance for each segment, [km]
               temp_best_traj(i,props_potcov)  =contrail_coverage(i)   !potcov for each segment, [km] 20170607-1
               temp_vtas_along_traj(i)         =vtas_along_traj(i)     !vtas [m/s] 20140524-5 
               !Yin_20170423
               !temp_cpc_along_traj(i)         =cpc_along_traj(i)     !CPC [km] 20170606-4,20170607-2
               !Yin_20170423
            ENDDO
         ENDIF
      CASE(2)   !Fuel_opt=2 20160606-4
         IF(Total_fuel_use.lt.temp_best_fuel_use)THEN
            temp_best_fl_time = T_FLIGHT_TIME
            temp_best_gen     = num_gen
            temp_best_fuel_use= Total_fuel_use 
            !Yin_20170423
            temp_best_potcov  = T_CPC                                  !Total CPC [km] 20170607-2
            !Yin_20170423
            !Outputs for upper SMCL(ac_routes properties)
            DO i=1,DIV_WAY+1
               temp_best_traj(i,props_lon)     =XXN(i)                 !Lon
               temp_best_traj(i,props_lat)     =YYN(i)                 !Lat
               temp_best_traj(i,props_alt)     =ZZN(i)                 !Alt [m]
               temp_best_traj(i,props_time)    =FLIGHT_TIME_JULIAN(i)  !Passing time at each waypoint, [Julian date]
               temp_best_traj(i,props_ac_speed)=FLIGHT_SPEED(i)*3.6_dp !AC_ground_speed[km/h],20140423-2,FLIGHT_SPEED[m/s]20140523-2 
               temp_best_traj(i,props_dist)    =FLIGHT_DISTANCE(i)     !Distance for each segment, [km]
               temp_best_traj(i,props_potcov)  =contrail_coverage(i)   !potcov for each segment, [km] 20170607-1
               temp_vtas_along_traj(i)         =vtas_along_traj(i)     !vtas [m/s] 20140524-5 
               !Yin_20170423
               !temp_cpc_along_traj(i)         =cpc_along_traj(i)     !CPC [km] 20170606-4,20170607-2
               !Yin_20170423
               !temp_best_traj(i,props_fuel_use)=                       !IMPORTANT:20160606-4,20170216-1 
            ENDDO
         ENDIF
      CASE(3)   !NOx_opt=3 20170216-1
         IF(Total_nox_emis.lt.temp_best_nox_emis)THEN
            temp_best_fl_time = T_FLIGHT_TIME
            temp_best_gen     = num_gen
            temp_best_nox_emis= Total_nox_emis                         ![g(NO2)] 
            !Yin_20170423
            temp_best_potcov  = T_CPC                                  !Total CPC [km] 20170607-2
            !Yin_20170423
            !Outputs for upper SMCL(ac_routes properties)
            DO i=1,DIV_WAY+1
               temp_best_traj(i,props_lon)     =XXN(i)                 !Lon
               temp_best_traj(i,props_lat)     =YYN(i)                 !Lat
               temp_best_traj(i,props_alt)     =ZZN(i)                 !Alt [m]
               temp_best_traj(i,props_time)    =FLIGHT_TIME_JULIAN(i)  !Passing time at each waypoint, [Julian date]
               temp_best_traj(i,props_ac_speed)=FLIGHT_SPEED(i)*3.6_dp !AC_ground_speed[km/h],20140423-2,FLIGHT_SPEED[m/s]20140523-2 
               temp_best_traj(i,props_dist)    =FLIGHT_DISTANCE(i)     !Distance for each segment, [km]
               temp_best_traj(i,props_potcov)  =contrail_coverage(i)   !potcov for each segment, [km] 20170607-1
               temp_vtas_along_traj(i)         =vtas_along_traj(i)     !vtas [m/s] 20140524-5 
               !Yin_20170423
               !temp_cpc_along_traj(i)         =cpc_along_traj(i)     !CPC [km] 20170606-4,20170607-2
               !Yin_20170423
               !temp_best_traj(i,props_fuel)    =                       !IMPORTANT:20160606-4 
               !temp_best_traj(i,props_emis_nox)=                       !IMPORTANT:20170216-1
            ENDDO
         ENDIF
      CASE(4)   !H2O_opt=4 20170222-4
         IF(Total_h2o_emis.lt.temp_best_h2o_emis)THEN
            temp_best_fl_time = T_FLIGHT_TIME
            temp_best_gen     = num_gen
            temp_best_h2o_emis= Total_h2o_emis 
            !Yin_20170423
            temp_best_potcov  = T_CPC                                  !Total CPC [km] 20170607-2
            !Yin_20170423
            !Outputs for upper SMCL(ac_routes properties)
            DO i=1,DIV_WAY+1
               temp_best_traj(i,props_lon)     =XXN(i)                 !Lon
               temp_best_traj(i,props_lat)     =YYN(i)                 !Lat
               temp_best_traj(i,props_alt)     =ZZN(i)                 !Alt [m]
               temp_best_traj(i,props_time)    =FLIGHT_TIME_JULIAN(i)  !Passing time at each waypoint, [Julian date]
               temp_best_traj(i,props_ac_speed)=FLIGHT_SPEED(i)*3.6_dp !AC_ground_speed[km/h],20140423-2,FLIGHT_SPEED[m/s]20140523-2 
               temp_best_traj(i,props_dist)    =FLIGHT_DISTANCE(i)     !Distance for each segment, [km]
               temp_best_traj(i,props_potcov)  =contrail_coverage(i)   !potcov for each segment, [km] 20170607-1
               temp_vtas_along_traj(i)         =vtas_along_traj(i)     !vtas [m/s] 20140524-5 
               !Yin_20170423
               !temp_cpc_along_traj(i)         =cpc_along_traj(i)     !CPC [km] 20170606-4,20170607-2
               !Yin_20170423
               !temp_best_traj(i,props_fuel_use)=                       !IMPORTANT:20160606-4,20170216-1 
            ENDDO
         ENDIF
      CASE(5)   !Contrail_opt=5 20170608-1
         !Yin_20170423 
         IF(T_CPC.lt.temp_best_potcov)THEN
!           write(*,*)"optimization works"
            temp_best_fl_time = T_FLIGHT_TIME
            temp_best_gen     = num_gen
            temp_best_potcov  = T_CPC                                  !Total CPC [km] 20170607-2,20170608-1
            !Outputs for upper SMCL(ac_routes properties)
            DO i=1,DIV_WAY+1
               temp_best_traj(i,props_lon)     =XXN(i)                 !Lon
               temp_best_traj(i,props_lat)     =YYN(i)                 !Lat
               temp_best_traj(i,props_alt)     =ZZN(i)                 !Alt [m]
               temp_best_traj(i,props_time)    =FLIGHT_TIME_JULIAN(i)  !Passing time at each waypoint, [Julian date]
               temp_best_traj(i,props_ac_speed)=FLIGHT_SPEED(i)*3.6_dp !AC_ground_speed[km/h],20140423-2,FLIGHT_SPEED[m/s]20140523-2 
               temp_best_traj(i,props_dist)    =FLIGHT_DISTANCE(i)     !Distance for each segment, [km]
               temp_best_traj(i,props_potcov)  =contrail_coverage(i)   !potcov for each segment, [km] 20170607-1,20170608-1
               temp_vtas_along_traj(i)         =vtas_along_traj(i)     !vtas [m/s] 20140524-5 
               !temp_cpc_along_traj(i)         =cpc_along_traj(i)      !CPC [km] 20170606-4,20170607-2 
            ENDDO
         ENDIF
         !Yin_20170423
      CASE(6)   !Cost_opt=6 20170224-3
         IF(Total_cost.lt.temp_best_cost)THEN
            temp_best_fl_time = T_FLIGHT_TIME
            temp_best_gen     = num_gen
            temp_best_cost    = Total_cost 
            !Yin_20170423
            temp_best_potcov  = T_CPC                                  !Total CPC [km] 20170607-2
            !Yin_20170423
            !Outputs for upper SMCL(ac_routes properties)
            DO i=1,DIV_WAY+1
               temp_best_traj(i,props_lon)     =XXN(i)                 !Lon
               temp_best_traj(i,props_lat)     =YYN(i)                 !Lat
               temp_best_traj(i,props_alt)     =ZZN(i)                 !Alt [m]
               temp_best_traj(i,props_time)    =FLIGHT_TIME_JULIAN(i)  !Passing time at each waypoint, [Julian date]
               temp_best_traj(i,props_ac_speed)=FLIGHT_SPEED(i)*3.6_dp !AC_ground_speed[km/h],20140423-2,FLIGHT_SPEED[m/s]20140523-2 
               temp_best_traj(i,props_dist)    =FLIGHT_DISTANCE(i)     !Distance for each segment, [km]
               temp_best_traj(i,props_potcov)  =contrail_coverage(i)   !potcov for each segment, [km] 20170607-1
               temp_vtas_along_traj(i)         =vtas_along_traj(i)     !vtas [m/s] 20140524-5 
               !Yin_20170423
               !temp_cpc_along_traj(i)         =cpc_along_traj(i)     !CPC [km] 20170606-4,20170607-2
               !Yin_20170423
               !temp_best_traj(i,props_fuel_use)=                       !IMPORTANT:20160606-4,20170216-1 
            ENDDO
         ENDIF
      CASE(7)   !COC_opt=7 20170316-4
         IF(Total_coc.lt.temp_best_coc)THEN
            temp_best_fl_time = T_FLIGHT_TIME
            temp_best_gen     = num_gen
            temp_best_coc     = Total_coc 
            !Yin_20170423
            temp_best_potcov  = T_CPC                                  !Total CPC [km] 20170607-2
            !Yin_20170423
            !Outputs for upper SMCL(ac_routes properties)
            DO i=1,DIV_WAY+1
               temp_best_traj(i,props_lon)     =XXN(i)                 !Lon
               temp_best_traj(i,props_lat)     =YYN(i)                 !Lat
               temp_best_traj(i,props_alt)     =ZZN(i)                 !Alt [m]
               temp_best_traj(i,props_time)    =FLIGHT_TIME_JULIAN(i)  !Passing time at each waypoint, [Julian date]
               temp_best_traj(i,props_ac_speed)=FLIGHT_SPEED(i)*3.6_dp !AC_ground_speed[km/h],20140423-2,FLIGHT_SPEED[m/s]20140523-2 
               temp_best_traj(i,props_dist)    =FLIGHT_DISTANCE(i)     !Distance for each segment, [km]
               temp_best_traj(i,props_potcov)  =contrail_coverage(i)   !potcov for each segment, [km] 20170607-1
               temp_vtas_along_traj(i)         =vtas_along_traj(i)     !vtas [m/s] 20140524-5 
            ENDDO
         ENDIF
      CASE(8)   !ATR20_opt=8 20170801
         IF(Total_atr20.lt.temp_best_atr20)THEN
            temp_best_fl_time = T_FLIGHT_TIME
            temp_best_gen     = num_gen
            temp_best_atr20   = Total_atr20 
            !Yin_20170423
            temp_best_potcov  = T_CPC                                  !Total CPC [km] 20170607-2
            !Yin_20170423
            !Outputs for upper SMCL(ac_routes properties)
            DO i=1,DIV_WAY+1
               temp_best_traj(i,props_lon)     =XXN(i)                 !Lon
               temp_best_traj(i,props_lat)     =YYN(i)                 !Lat
               temp_best_traj(i,props_alt)     =ZZN(i)                 !Alt [m]
               temp_best_traj(i,props_time)    =FLIGHT_TIME_JULIAN(i)  !Passing time at each waypoint, [Julian date]
               temp_best_traj(i,props_ac_speed)=FLIGHT_SPEED(i)*3.6_dp !AC_ground_speed[km/h],20140423-2,FLIGHT_SPEED[m/s]20140523-2 
               temp_best_traj(i,props_dist)    =FLIGHT_DISTANCE(i)     !Distance for each segment, [km]
               temp_best_traj(i,props_potcov)  =contrail_coverage(i)   !potcov for each segment, [km] 20170607-1
               temp_vtas_along_traj(i)         =vtas_along_traj(i)     !vtas [m/s] 20140524-5 
            ENDDO
         ENDIF
      CASE(9)   !COSTCLIM_opt=9 20170801
         IF(Total_costclim.lt.temp_best_costclim)THEN
            temp_best_fl_time = T_FLIGHT_TIME
            temp_best_gen     = num_gen
            temp_best_costclim= Total_costclim 
            !Yin_20170423
            temp_best_potcov  = T_CPC                                  !Total CPC [km] 20170607-2
            !Yin_20170423
            !Outputs for upper SMCL(ac_routes properties)
            DO i=1,DIV_WAY+1
               temp_best_traj(i,props_lon)     =XXN(i)                 !Lon
               temp_best_traj(i,props_lat)     =YYN(i)                 !Lat
               temp_best_traj(i,props_alt)     =ZZN(i)                 !Alt [m]
               temp_best_traj(i,props_time)    =FLIGHT_TIME_JULIAN(i)  !Passing time at each waypoint, [Julian date]
               temp_best_traj(i,props_ac_speed)=FLIGHT_SPEED(i)*3.6_dp !AC_ground_speed[km/h],20140423-2,FLIGHT_SPEED[m/s]20140523-2 
               temp_best_traj(i,props_dist)    =FLIGHT_DISTANCE(i)     !Distance for each segment, [km]
               temp_best_traj(i,props_potcov)  =contrail_coverage(i)   !potcov for each segment, [km] 20170607-1
               temp_vtas_along_traj(i)         =vtas_along_traj(i)     !vtas [m/s] 20140524-5 
            ENDDO
         ENDIF
      CASE(10)   !CPCCOST_opt=10 20170801
         IF(Total_costcpc.lt.temp_best_costcpc)THEN
            temp_best_fl_time = T_FLIGHT_TIME
            temp_best_gen     = num_gen
            temp_best_costcpc= Total_costcpc 
            !Yin_20170423
            temp_best_potcov  = T_CPC                                  !Total CPC [km] 20170607-2
            !Yin_20170423
            !Outputs for upper SMCL(ac_routes properties)
            DO i=1,DIV_WAY+1
               temp_best_traj(i,props_lon)     =XXN(i)                 !Lon
               temp_best_traj(i,props_lat)     =YYN(i)                 !Lat
               temp_best_traj(i,props_alt)     =ZZN(i)                 !Alt [m]
               temp_best_traj(i,props_time)    =FLIGHT_TIME_JULIAN(i)  !Passing time at each waypoint, [Julian date]
               temp_best_traj(i,props_ac_speed)=FLIGHT_SPEED(i)*3.6_dp !AC_ground_speed[km/h],20140423-2,FLIGHT_SPEED[m/s]20140523-2 
               temp_best_traj(i,props_dist)    =FLIGHT_DISTANCE(i)     !Distance for each segment, [km]
               temp_best_traj(i,props_potcov)  =contrail_coverage(i)   !potcov for each segment, [km] 20170607-1
               temp_vtas_along_traj(i)         =vtas_along_traj(i)     !vtas [m/s] 20140524-5 
            ENDDO
         ENDIF

      CASE DEFAULT
         write(*,*)"Error:Find temporary best solution"
      ENDSELECT
      
  ENDDO   !HY20140415-5
!------------------------------------------
!Population loop for one generation(END)
!------------------------------------------
!  DEALLOCATE(w_zalt)                       !Array, 20150312-5 

  !output for congergence history, 20140704-2,20140706-2
  !IF((INT(gc_D_lon).eq.11).and.(INT(gc_D_lat).eq.48).and.(INT(gc_A_lon).eq.-74).and.(INT(gc_A_lat).eq.40))THEN
  !OPEN(21,file='converge.dat',position='append')
  !      write(21,*)temp_best_gen,temp_best_fl_time
  !CLOSE(21)
  !ENDIF

!Output: Vgound/Vtas for 3 city-pairs, 20150311-2
!Execute followings, if current case is one of the 3 city-pairs.20150312-1,2, 20160606-5, 20170214-1
  IF(i_add_nnindex.eq.1)THEN   !for the 3 city-pairs
     DO j=1,101            !101:num of traj. 
        w_zalt = 0.0_dp         
        wind_along_traj=0.0_dp  
        vtas_along_traj=0.0_dp  
        rho_along_traj=0.0_dp
        !Yin_20170423
        cpc_along_traj=0.0_dp
        ATR20CPC_along_traj=0.0_dp
        !Yin_20170423

        !Determine waypoints of candidated-trajectory.This is instead of "CALL arbitrary_traj" 
        DO i=1,DIV_WAY+1   !num of wp
           traj_for_nnindex(i,1,j) = temp_best_traj(i,props_lon)   !Lon
           traj_for_nnindex(i,2,j) = temp_best_traj(i,props_lat)   !Lat
           traj_for_nnindex(i,3,j) = 0.0_dp + 150.0_dp*(j-1)       !Alt, [m]
        ENDDO
        !arrange of longitude
        temp_XXN=0.0_dp   
        DO i=1,DIV_WAY+1
           !If [-180=< XXN(i) < 0], add +360 to XXN(i) and use temp_XXN for nn_index search.20140605-2
           IF((-180.0_dp.le.traj_for_nnindex(i,1,j)).and.(traj_for_nnindex(i,1,j).lt.0.0_dp))THEN
              temp_XXN = traj_for_nnindex(i,1,j)+360.0_dp
              CALL nn_index(w_philon(:), temp_XXN, w_j_lon) !lon value on waypoints,20140605-1
           ELSE
              temp_XXN = traj_for_nnindex(i,1,j)
              CALL nn_index(w_philon(:), temp_XXN, w_j_lon) !lon value on waypoints
           ENDIF

           CALL nn_index(w_philat(:), traj_for_nnindex(i,2,j), w_i_lat)  !lat value on waypoints
           w_zalt(:)=w_zgl_geopot_3d(w_j_lon, :, w_i_lat)/g   ![m] geopotential altitude,20131105-1,20140430-1
           CALL nn_index(w_zalt(:), traj_for_nnindex(i,3,j), w_i_alt)  
           wind_along_traj(i,1) = w_uwind_g(w_j_lon, w_i_alt, w_i_lat)    ![m/s]
           wind_along_traj(i,2) = w_vwind_g(w_j_lon, w_i_alt, w_i_lat)    ![m/s]
           wind_along_traj(i,3) =   w_v_z_g(w_j_lon, w_i_alt, w_i_lat)    !global vertical wind velocity field [m/s]
           !Yin_20170423
           cpc_along_traj(i)    =   w_cpc_g(w_j_lon, w_i_alt, w_i_lat)    !global potcov field [-]
           IF(w_option_traj_calc_g==8)THEN   !20180525 
           ATR20CPC_along_traj(i)    =   w_ATR20CPC_3d_g(w_j_lon, w_i_alt, w_i_lat)    !global potcov field [-]
           ENDIF 
           !Yin_20170423
           vtas_along_traj(i)   = AC_MACH*SQRT(w_t_scb_g(w_j_lon, w_i_alt, w_i_lat)*   &
                                  adiabatic_index_air*gas_constant_air)   ![m/s] 20140724-3
           !For fuel_opt, NOx_opt, H2O_opt, Cost_opt, COC_opt.
           !For time_opt, CPC_opt, rho_along_traj(i)=0 by line 1188.
           IF(w_option_traj_calc_g==2 .or. w_option_traj_calc_g==3 .or. w_option_traj_calc_g==4 .or.  &
              w_option_traj_calc_g==6 .or. w_option_traj_calc_g==7 .or. w_option_traj_calc_g==8.or. &
              w_option_traj_calc_g==9 .or. w_option_traj_calc_g==10)THEN 
              !20170216-2,20170221-1,20170222-4,20170224-4,20170316-4,20170608-1
              rho_along_traj(i) = w_rho_air_dry_3d_g(w_j_lon, w_i_alt, w_i_lat) !20160623-2,20170214-2,20170607-3 
           ENDIF
        ENDDO
        CALL gc_ac_traj(traj_for_nnindex(:,1,j), traj_for_nnindex(:,2,j), traj_for_nnindex(:,3,j),    &
                        cpc_along_traj, wind_along_traj, vtas_along_traj,  &     !vtas_along_traj 20140522-4
                        FLIGHT_TIME, FLIGHT_TIME_JULIAN, T_FLIGHT_TIME,         &
                        A_time, FLIGHT_DISTANCE, T_FLIGHT_DISTANCE, FLIGHT_SPEED, contrail_coverage, &
                        T_CPC, ATR20CPC_along_traj, Total_atr20cpc) !FLIGHT_SPEED = ac_gs[m/s]
        DO i=1,DIV_WAY+1   !num of wp
           !Regarding values of temp_best_vtas, temp_best_vground, temp_best_rho for each routing option,           
           !they are listed in 20170216-4
           !These arrays are passed to messy_airtraf_tools_ga.f90
           temp_best_vtas(i,j)    = vtas_along_traj(i)      !Vtas,[m/s] 
           temp_best_vground(i,j) = FLIGHT_SPEED(i)         !Vground,[m/s] 
           !For time_opt, CPC_opt, rho_along_traj(i)=0.0  20170214-2,20170608-1
           temp_best_rho(i,j) = rho_along_traj(i)           !rho [kg/m^3] 20160623-2,20170216-4
           temp_best_cpc(i,j) = cpc_along_traj(i)           !potcov [-] 20170608-2
        ENDDO
     ENDDO
  ELSEIF(i_add_nnindex.eq.0)THEN     ! for other city-pairs
     !These arrays are passed to messy_airtraf_tools_ga.f90
     temp_best_vtas    =0.0_dp                           !Vtas,[m/s] 
     temp_best_vground =0.0_dp                           !Vground,[m/s] 
     temp_best_rho     =0.0_dp                           !rho,[kg/m^3] 20160623-2 
     temp_best_cpc     =0.0_dp                           !potcov,[-] 20170608-2
  ENDIF

  DEALLOCATE(w_zalt)                       !Array, 20150312-5 

  END SUBROUTINE calcobj 
!----------------------------------------------------------------------

!----------------------------------------------------------------------
  SUBROUTINE gc_distance(COSIS, T_DIS, T_DISN)
!----------------------------------------------------------------------
!
! This subroutine calcuate great circre distance.
!
  IMPLICIT NONE
  REAL(DP), INTENT(OUT) :: COSIS, T_DIS, T_DISN 

  SELECT CASE (ISWITCH)  
  CASE(0)
     DELRA  = RARAD1-RARAD2
     COSIS  = ACOS(SIN(DCRAD1)*SIN(DCRAD2) + COS(DCRAD1)*COS(DCRAD2)*&
                   COS(DELRA))               !in [rad]

  CASE(1) 
     DELDC2 = (DCRAD1-DCRAD2)/2.0_dp            !delta lat
     DELRA2 = (RARAD1-RARAD2)/2.0_dp            !delta lon
     SINDIS = SQRT(SIN(DELDC2)*SIN(DELDC2) + COS(DCRAD1)*COS(DCRAD2)*&
                   SIN(DELRA2)*SIN(DELRA2))
     COSIS  = 2.0_dp*DASIN(SINDIS)              !in [rad]

  CASE(2)
     DELRA2 = RARAD1-RARAD2                  !delta lon
     SINDIS1= SQRT((COS(DCRAD2)*SIN(DELRA2))**2 + (COS(DCRAD1)*SIN(DCRAD2)-&
                    SIN(DCRAD1)*COS(DCRAD2)*COS(DELRA2))**2)
     SINDIS2= SIN(DCRAD1)*SIN(DCRAD2) + COS(DCRAD1)*COS(DCRAD2)*COS(DELRA2) 
     COSIS  = DATAN2(SINDIS1,SINDIS2)        !in [rad]

  CASE DEFAULT
     write(*,*) 'ISWITCH error: Unknown ISWITCH on GC distance calculation!'
  ENDSELECT

! T_DIS    = RADIUS*COSIS                    !GC distance = radius * radians [km]
! T_DIS    = (RADIUS + ALT/1000.0)*COSIS     !GC distance = radius * radians [km] 
!20140524-10  T_DIS    = (radius_earth*0.001_dp + ALT/1000.0_dp)*COSIS     !GC distance = radius * radians [km] 
!20140524-10  T_DISN   = T_DIS*0.539957_dp                  !convert T_DIS into nautical miles [nm]

  END SUBROUTINE gc_distance
!----------------------------------------------------------------------

!----------------------------------------------------------------------
  SUBROUTINE arbitrary_traj(XX, YY, ZZ)
!----------------------------------------------------------------------
!  HY20130206-4, 20140423-4, opt20121023-10
!  Ranges of dvs:
!    -180 =< dv1,dv3,dv5 =< 180 [FL_DIR=0,2,4,5]
!       0 =< dv1,dv3,dv5 =< 360 [FL_DIR=1,3]
!
!    -180 =< lon_dv7 to lon_dv11 =< 180 [FL_DIR=0,2]
!       0 =< lon_dv7 to lon_dv11 =< 360 [FL_DIR=1,3]
!     -90 =< lon_dv7 to lon_dv11 =< +90 [FL_DIR=4,5](<- those variables mean positions in latitudes)
!
!     -90 =< gc_D_lat =< +90, -180 =< gc_D_lon =< 180 
!     -90 =< gc_A_lat =< +90, -180 =< gc_A_lon =< 180 
!
!     180 =< temp_gc_D_lon =< 540 [FL_DIR=1]
!     180 =< temp_gc_A_lon =< 540 [FL_DIR=3]

 IMPLICIT NONE
! When we use b_spline.f, use these declerations below.
 INTEGER, PARAMETER                  :: lonlat_case = 1  !HY20140424-1
 INTEGER, PARAMETER                  :: lonalt_case = 2  !HY20140424-1
!HY INTEGER, PARAMETER                  :: internal_lonlat_case = 3  !
 REAL(DP)                            :: lon_cp(NUM_CP_LONLAT+2), lat_cp(NUM_CP_LONLAT+2)
 REAL(DP)                            :: lon4alt_cp(NUM_CP_ALT+2), alt_cp(NUM_CP_ALT+2)
 REAL(DP)                            :: lat4alt_cp(NUM_CP_ALT+2)
!HY20140424-5 REAL(DP), DIMENSION(:) :: new_lon(DIV_WAY+1), new_lat(DIV_WAY+1)
 REAL(DP), DIMENSION(:), ALLOCATABLE :: new_lon, new_lat
!HY20140424-5 REAL(DP), DIMENSION(:) :: new_alt(DIV_WAY+1), new_lon4alt(DIV_WAY+1)
 REAL(DP), DIMENSION(:), ALLOCATABLE :: new_alt, new_lon4alt
!HY20140424-5 REAL(DP), DIMENSION(:) :: new_lat4alt(DIV_WAY+1)
 REAL(DP), DIMENSION(:), ALLOCATABLE :: new_lat4alt
 REAL(DP)              , INTENT(OUT) :: XX(DIV_WAY+1), YY(DIV_WAY+1), ZZ(DIV_WAY+1)
 INTEGER                             :: i, j, new_out

 !For wind components calculation in gc_ac_traj. 
!HY REAL(DP), DIMENSION(:), ALLOCATABLE :: internal_new_lon, internal_new_lat
!HY INTEGER                             :: internal_new_out



! When we use gcontrol.f, use these declerations below.
! REAL, DIMENSION(:) :: lon_cp(NUM_CP_LONLAT+2), lat_cp(NUM_CP_LONLAT+2), alt_cp(NUM_CP_ALT+2)
! REAL, DIMENSION(:) :: new_lon(DIV_WAY+1),new_lat(DIV_WAY+1),new_alt(DIV_WAY+1),SN(DIV_WAY+1)
! DOUBLE PRECISION, DIMENSION(:),INTENT(OUT) :: XX(DIV_WAY+1),YY(DIV_WAY+1),ZZ(DIV_WAY+1)
! INTEGER :: i,j

 REAL(DP)          :: ratio_btw_waypoints
 INTEGER           :: lower_waypoint, upper_waypoint

  !Initilizations
  lon_cp      = 0.0_dp
  lat_cp      = 0.0_dp
  alt_cp      = 0.0_dp
  lon4alt_cp  = 0.0_dp
  lat4alt_cp  = 0.0_dp
!HY20140424-5  new_lon     = 0.0_dp
!HY20140424-5  new_lat     = 0.0_dp
!HY20140424-5  new_alt     = 0.0_dp
!HY20140424-5  new_lon4alt = 0.0_dp
!HY20140424-5  new_lat4alt = 0.0_dp
! SN          = 0.0_dp     !For grcv3d.
  XX          = 0.0_dp
  YY          = 0.0_dp
  ZZ          = 0.0_dp

  !CALCULATE (LON, LAT) TRAJ USING B_SPLINE 
  SELECT CASE(FL_DIR)
  CASE(0,1,2,3)
     !LONGITUDE
     IF(FL_DIR.eq.1)THEN
          lon_cp(1)  = temp_gc_D_lon !fixed
       ELSEIF(FL_DIR.ne.1)THEN   !For FL_DIR=0,2,3
          lon_cp(1)  = gc_D_lon      !fixed
     ENDIF
!HY20140415-6     lon_cp(2)  = dv1
!HY20140415-6     lon_cp(3)  = dv3
!HY20140415-6     lon_cp(4)  = dv5
     lon_cp(2)  = DV(1)
     lon_cp(3)  = DV(3)
     lon_cp(4)  = DV(5)
     IF(FL_DIR.eq.3)THEN
          lon_cp(5)  = temp_gc_A_lon !fixed
       ELSEIF(FL_DIR.ne.3)THEN   !For FL_DIR=0,1,2
          lon_cp(5)  = gc_A_lon      !fixed
     ENDIF
     !LATITUDE
     lat_cp(1)  = gc_D_lat      !fixed
!HY20140415-6     lat_cp(2)  = dv2
!HY20140415-6     lat_cp(3)  = dv4
!HY20140415-6     lat_cp(4)  = dv6
     lat_cp(2)  = DV(2)
     lat_cp(3)  = DV(4)
     lat_cp(4)  = DV(6)
     lat_cp(5)  = gc_A_lat      !fixed

  CASE(4,5)
     !LATITUDE
     lat_cp(1)  = gc_D_lat      !fixed
!HY20140415-6     lat_cp(2)  = dv2
!HY20140415-6     lat_cp(3)  = dv4
!HY20140415-6     lat_cp(4)  = dv6
     lat_cp(2)  = DV(2)
     lat_cp(3)  = DV(4)
     lat_cp(4)  = DV(6)
     lat_cp(5)  = gc_A_lat      !fixed

     !LONGITUDE
     lon_cp(1)  = gc_D_lon      !fixed
!HY20140415-6     lon_cp(2)  = dv1
!HY20140415-6     lon_cp(3)  = dv3
!HY20140415-6     lon_cp(4)  = dv5
     lon_cp(2)  = DV(1)
     lon_cp(3)  = DV(3)
     lon_cp(4)  = DV(5)
     lon_cp(5)  = gc_A_lon      !fixed

  END SELECT

  !ALTITUDE (For grcv3d)
  ! alt_cp(1)  = dv7  !ALT
  ! alt_cp(2)  = dv8  !ALT
  ! alt_cp(3)  = dv9  !ALT
  ! alt_cp(4)  = dv10 !ALT
  ! alt_cp(5)  = dv11 !ALT

!  do i=1, NUM_CP_LONLAT+2
!     write(*,*)'Control points for (LON, LAT) TRAJ:',lon_cp(i),lat_cp(i)
!  enddo

  SELECT CASE(FL_DIR)
  CASE(0,1,2,3)
     CALL B_SPLINE(lonlat_case, DIV_WAY+1, NUM_CP_LONLAT+2, lon_cp, lat_cp, new_out, new_lon, new_lat)  !HY20140423-4
!    CALL B_SPLINE(1000,NUM_CP,lon_cp,lat_cp,new_out,new_lon,new_lat)
!    CALL grcv3d(5,lon_cp,lat_cp,alt_cp,DIV_WAY+1,1,1.0,1.0,new_lon,new_lat,new_alt,SN)
     
     !For wind components calculation !HY20140506-3
!HY     CALL B_SPLINE(internal_lonlat_case, DIV_WAY+1, NUM_CP_LONLAT+2, lon_cp, lat_cp, internal_new_out, &
!HY                   internal_new_lon, internal_new_lat)
  CASE(4,5)
     CALL B_SPLINE(lonlat_case, DIV_WAY+1, NUM_CP_LONLAT+2, lat_cp, lon_cp, new_out, new_lat, new_lon)  !HY20140423-4
     
     !For wind components calculation !HY20140506-3
!HY     CALL B_SPLINE(internal_lonlat_case, DIV_WAY+1, NUM_CP_LONLAT+2, lat_cp, lon_cp, internal_new_out, &
!HY                   internal_new_lat, internal_new_lon) 
  END SELECT

  !new_out = DIV_WAY+1(nwaypoints) basically 
  !XX(DIV_WAY+1),YY(DIV_WAY+1)
  !XX and YY do not pass through the exact control points.
  !write(*,*)"XX(j),YY(j),ZZ(j)="
  do j=1, DIV_WAY+1
     XX(j) = new_lon(j)   !determined waypoints(DIV_WAY) in lon
     YY(j) = new_lat(j)   !determined waypoints(DIV_WAY) in lat
!    ZZ(j) = ALT 
!    write(*,*)XX(j),YY(j),ZZ(j)
  enddo
 !write(*,*)'new_out=',new_out
!HY write(*,*)"inter_XX(j),YY(j),ZZ(j)="   !HY20140506-3
!HY do j=1, internal_new_out
!     XX(j) = new_lon(j)   !determined waypoints(DIV_WAY) in lon
!     YY(j) = new_lat(j)   !determined waypoints(DIV_WAY) in lat
!    ZZ(j) = ALT 
!HY    write(*,*)internal_new_lon(j), internal_new_lat(j)
!HY  enddo
!HY write(*,*)'inter_new_out=',internal_new_out


  !----------For grcv3d
!  write(*,*)'after_grcv3d'
!  do j=1, DIV_WAY+1
!     XX(j) = new_lon(j)*1.0d0
!     YY(j) = new_lat(j)*1.0d0
!     ZZ(j) = new_alt(j)*1.0d0 !ALT 
!     write(*,*)XX(j),YY(j),ZZ(j)
!  enddo
!-----------------------------

  !CALCULATE (LON, ALT) TRAJ USING B_SPLINE 
  SELECT CASE(FL_DIR)
  CASE(0,1,2,3)
     !LONGITUDE
     IF(FL_DIR.eq.1)THEN
          lon4alt_cp(1) = temp_gc_D_lon !fixed
       ELSEIF(FL_DIR.ne.1)THEN  !For FL_DIR= 0,2,3
          lon4alt_cp(1) = gc_D_lon      !fixed 
     ENDIF
     lon4alt_cp(2) = lon_dv7
     lon4alt_cp(3) = lon_dv8
     lon4alt_cp(4) = lon_dv9
     lon4alt_cp(5) = lon_dv10
     lon4alt_cp(6) = lon_dv11
     IF(FL_DIR.eq.3)THEN
          lon4alt_cp(7) = temp_gc_A_lon !fixed
       ELSEIF(FL_DIR.ne.3)THEN !For FL_DIR= 0,1,2
          lon4alt_cp(7) = gc_A_lon      !fixed
     ENDIF
     !ALTITUDE
     alt_cp(1)     = Lower_alt  !ALT20140522-1          !fixed
!HY20140415-6     alt_cp(2)     = dv7
!HY20140415-6     alt_cp(3)     = dv8
!HY20140415-6     alt_cp(4)     = dv9
!HY20140415-6     alt_cp(5)     = dv10
!HY20140415-6     alt_cp(6)     = dv11
     alt_cp(2)     = DV(7)
     alt_cp(3)     = DV(8)
     alt_cp(4)     = DV(9)
     alt_cp(5)     = DV(10)
     alt_cp(6)     = DV(11)
     alt_cp(7)     = Lower_alt   !ALT20140522-1          !fixed

  !CALCULATE (LAT, ALT) TRAJ USING B_SPLINE 
  CASE(4,5)
     !LATITUDE
     lat4alt_cp(1) = gc_D_lat
     lat4alt_cp(2) = lon_dv7
     lat4alt_cp(3) = lon_dv8
     lat4alt_cp(4) = lon_dv9
     lat4alt_cp(5) = lon_dv10
     lat4alt_cp(6) = lon_dv11
     lat4alt_cp(7) = gc_A_lat

     !ALTITUDE
     alt_cp(1)     = Lower_alt   !ALT20140522-1          !fixed
!HY20140415-6     alt_cp(2)     = dv7
!HY20140415-6     alt_cp(3)     = dv8
!HY20140415-6     alt_cp(4)     = dv9
!HY20140415-6     alt_cp(5)     = dv10
!HY20140415-6     alt_cp(6)     = dv11
     alt_cp(2)     = DV(7)
     alt_cp(3)     = DV(8)
     alt_cp(4)     = DV(9)
     alt_cp(5)     = DV(10)
     alt_cp(6)     = DV(11)
     alt_cp(7)     = Lower_alt   !ALT20140522-1          !fixed

  END SELECT

!  do i=1, NUM_CP_ALT+2
!     write(*,*)'Control points for (LON, ALT) TRAJ:',lon4alt_cp(i),alt_cp(i)
!  enddo

  SELECT CASE(FL_DIR)
  CASE(0,1,2,3)
     CALL B_SPLINE(lonalt_case, DIV_WAY+1, NUM_CP_ALT+2, lon4alt_cp, alt_cp, new_out, new_lon4alt, new_alt)  !HY20140424-1
  CASE(4,5)
     CALL B_SPLINE(lonalt_case, DIV_WAY+1, NUM_CP_ALT+2, lat4alt_cp, alt_cp, new_out, new_lat4alt, new_alt)  !HY20140424-1
  END SELECT

  !NOTE: new_out =\ DIV_WAY+1(nwaypoints), see 20121023-10
  !new_lon4alt(new_out), new_alt(new_out)
  !new_lon4alt and new_alt do not pass through the exact control points.
!  write(*,*)'CHECK_(LON, ALT) TRAJ. new_lon4alt, new_alt='
!  do i=1, new_out  
!     write(*,*)new_lon4alt(i), new_alt(i)
!  enddo
!  write(*,*)"new_out=",new_out 

  !Calculate ZZ(j) corresponding to 'XX(j)' generated by (LON, LAT) traj B-spline interpolation.  
  !ZZ(DIV_WAY+1)
  ZZ(1)         = new_alt(1)
  ZZ(DIV_WAY+1) = new_alt(new_out)

  !NOTE: new_out =\ DIV_WAY+1, see 20121023-10
  SELECT CASE(FL_DIR)
  CASE(0,3)   !opt20130208-2, opt20130208-4, 20140426-3, 20141211-1
  OUTER1: DO j=2, DIV_WAY   !Loop for XX
     INNER1: DO i=1, new_out-1   !Loop for new_lon4alt
       IF(new_lon4alt(i)<=XX(j) .and. new_lon4alt(i+1)>XX(j))THEN  !FL_DIR=0 
           lower_waypoint = i
           upper_waypoint = i+1
           ratio_btw_waypoints = (new_lon4alt(upper_waypoint)-XX(j))/   &
                                 (new_lon4alt(upper_waypoint)-new_lon4alt(lower_waypoint))
           ZZ(j)=ratio_btw_waypoints*new_alt(lower_waypoint) + (1.0_dp-ratio_btw_waypoints)*new_alt(upper_waypoint)
           exit INNER1
       ENDIF
     ENDDO INNER1
  ENDDO OUTER1

  CASE(1,2)  !opt20130207-2, opt20130208-1 
  OUTER2: DO j=2, DIV_WAY
     INNER2: DO i=1, new_out-1
       IF(new_lon4alt(i)>=XX(j) .and. new_lon4alt(i+1)<XX(j))THEN  !FL_DIR=2, MUC to JFK
           lower_waypoint = i
           upper_waypoint = i+1
           ratio_btw_waypoints = (XX(j)-new_lon4alt(lower_waypoint))/   &
                                 (new_lon4alt(upper_waypoint)-new_lon4alt(lower_waypoint))
           ZZ(j)=ratio_btw_waypoints*new_alt(upper_waypoint) + (1.0_dp-ratio_btw_waypoints)*new_alt(lower_waypoint)
           exit INNER2
       ENDIF
     ENDDO INNER2
  ENDDO OUTER2

  CASE(4)  !opt20130213-4
  OUTER3: DO j=2, DIV_WAY
     INNER3: DO i=1, new_out-1
       IF(new_lat4alt(i)>=YY(j) .and. new_lat4alt(i+1)<YY(j))THEN
           lower_waypoint = i
           upper_waypoint = i+1
           ratio_btw_waypoints = (YY(j)-new_lat4alt(lower_waypoint))/   &
                                 (new_lat4alt(upper_waypoint)-new_lat4alt(lower_waypoint))
           ZZ(j)=ratio_btw_waypoints*new_alt(upper_waypoint) + (1.0_dp-ratio_btw_waypoints)*new_alt(lower_waypoint)
           exit INNER3
       ENDIF
     ENDDO INNER3
  ENDDO OUTER3

  CASE(5)  !opt20130214-3
  OUTER4: DO j=2, DIV_WAY
     INNER4: DO i=1, new_out-1
       IF(new_lat4alt(i)<=YY(j) .and. new_lat4alt(i+1)>YY(j))THEN
           lower_waypoint = i
           upper_waypoint = i+1
           ratio_btw_waypoints = (new_lat4alt(upper_waypoint)-YY(j))/   &
                                 (new_lat4alt(upper_waypoint)-new_lat4alt(lower_waypoint))
           ZZ(j)=ratio_btw_waypoints*new_alt(lower_waypoint) + (1.0_dp-ratio_btw_waypoints)*new_alt(upper_waypoint)
           exit INNER4
       ENDIF
     ENDDO INNER4
  ENDDO OUTER4

  CASE DEFAULT
      write(*,*) 'subroutine arbitrary_traj error: this aircraft is not assinged of any FL_DIR!'
  END SELECT

  END SUBROUTINE arbitrary_traj
!----------------------------------------------------------------------

!----------------------------------------------------------------------
  SUBROUTINE gc_waypoints(THETA_OUT, PHI_OUT)
!----------------------------------------------------------------------
!
! This subroutine calculate waypoints along GC.
! Equations below shoulb be checked on precision!!!
! see20111114-2,20140507-4,5
  IMPLICIT NONE
  REAL(DP) :: THETA(DIV_WAY+1), PHI(DIV_WAY+1) 
  REAL(DP), INTENT(OUT) :: THETA_OUT(DIV_WAY+1), PHI_OUT(DIV_WAY+1) 
  REAL(DP) :: DELDC3, DELRA3
  INTEGER          :: i

  !Initialization
  DELDC3 = 0.0_dp
  DELRA3 = 0.0_dp
  THETA  = 0.0_dp
  PHI    = 0.0_dp

  DELRA3 = RARAD1 - RARAD2     !Delta lon between gc_D_lon & gc_A_lon
  DELDC3 = DCRAD1 - DCRAD2     !Delta lat between gc_D_lat & gc_A_lat

  IF((DELRA3.gt.0.0_dp).and.(DELRA3.le.PI))THEN
     DO i=1, DIV_WAY+1
        PHIR = RARAD2 + (i-1)*(RARAD1-RARAD2)/DIV_WAY              !in [rad]
        THETAR = atan(sin(DCRAD1)*sin(PHIR-RARAD2)/(cos(DCRAD1)*sin(RARAD1-RARAD2))+&
                      sin(DCRAD2)*sin(RARAD1-PHIR)/(cos(DCRAD2)*sin(RARAD1-RARAD2)))   !in[rad]
        PHI(i) = PHIR/DTR         !lon along wp [deg]
        THETA(i) = THETAR/DTR     !lat along wp [deg]
     ENDDO 
  ELSE IF(DELRA3.gt.PI)THEN
     DO i=1,DIV_WAY+1
        PHIR = RARAD1 + (i-1)*((RARAD2+2.0_dp*PI) - RARAD1)/DIV_WAY     !lon [rad]
        THETAR = atan(sin(DCRAD1)*sin(PHIR-RARAD2)/(cos(DCRAD1)*sin(RARAD1-RARAD2))+&
                      sin(DCRAD2)*sin(RARAD1-PHIR)/(cos(DCRAD2)*sin(RARAD1-RARAD2)))   !in[rad]
        PHI(i) = PHIR/DTR         !lon along wp [deg]
        THETA(i) = THETAR/DTR     !lat along wp [deg]
     ENDDO
  ELSE                            !(DELRA3.eq.0.0)
     DO i=1, DIV_WAY+1
        THETAR = DCRAD2 + (i-1)*DELDC3/DIV_WAY  !lon [rad]
        PHI(i) = RARAD1/DTR       !lon along wp [deg]
        THETA(i) = THETAR/DTR     !lat along wp [deg]   !in[rad]
     ENDDO 
  ENDIF

  !Waypoint output be in order
  !Cases FL_DIR = 0,3,5 are in same coding
  SELECT CASE (FL_DIR)  
  CASE(0,3,5)
     DO i=1,DIV_WAY+1   
       PHI_OUT(i) = PHI(i)
       THETA_OUT(i) = THETA(i)
     ENDDO
!  CASE(3) !Output is converted from 0 =<lon=< +360 into -180 =<lon=< +180 
!     DO i=1,DIV_WAY+1   
!       PHI_OUT(i) = PHI(i)           !lon along wp [deg]
!       THETA_OUT(i) = THETA(i)       !lat along wp [deg]
!     ENDDO
!     DO i=1,DIV_WAY+1
!       IF(PHI_OUT(i).gt.180.0d0)PHI_OUT(i) = PHI_OUT(i) - 360.0d0
!     ENDDO
   
  !Cases FL_DIR = 1,2,4 are in same coding
  CASE(1,2,4)
     DO i=1,DIV_WAY+1   
       PHI_OUT(i) = PHI(DIV_WAY+2-i)
       THETA_OUT(i) = THETA(DIV_WAY+2-i)
     ENDDO
!  CASE(1) !Output is converted from 0 =<lon=< +360 into -180 =<lon=< +180
!     DO i=1,DIV_WAY+1   
!       PHI_OUT(i) = PHI(DIV_WAY+2-i)      !lon along wp [deg]
!       THETA_OUT(i) = THETA(DIV_WAY+2-i)  !lat along wp [deg]
!     ENDDO
!     DO i=1,DIV_WAY+1
!       IF(PHI_OUT(i).gt.180.0d0)PHI_OUT(i) = PHI_OUT(i) - 360.0d0
!     ENDDO
  CASE DEFAULT 
      write(*,*) 'Output error: Unknown Waypoints on GC!'
  END SELECT

  END SUBROUTINE gc_waypoints
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!20140417-4  SUBROUTINE gc_ac_traj(T_LAT, T_LON, FL_T, FL_T_Julian, T_FL_T, T_FL_T_Julian, FL_DIST, T_FL_DIST) 
  SUBROUTINE gc_ac_traj(T_LON, T_LAT, T_ALT, cpc_value, wind_value, ac_vtas,      &
                        FL_T, FL_T_Julian, T_FL_T,                     &
                        T_FL_T_Julian, FL_DIST, T_FL_DIST, ac_gs, ac_cpc, &
                        T_potcov, atr20cpc_value, T_atr20cpc)   !FL_SPEED 20140523-1
!----------------------------------------------------------------------
!
! This subroutine calculate aircraft flight trajectories.
! Equations below should be checked on precision!!!
  IMPLICIT NONE
!HY20140430-2  REAL(DP), DIMENSION(:), INTENT(IN) :: T_LON(DIV_WAY+1), T_LAT(DIV_WAY+1), T_ALT(DIV_WAY+1)  
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: wind_value   ![m/s] 20140430-2 
  REAL(DP), DIMENSION(:), INTENT(IN)   :: ac_vtas      ![m/s] 20140522-4 
  !Yin_20170423
  REAL(DP), DIMENSION(:), INTENT(IN)   :: cpc_value    ![-] 20161215,20170607-1,3
  REAL(DP), DIMENSION(:), INTENT(IN),OPTIONAL   :: atr20cpc_value    ![-] 20170801
  !Yin_20170423
  REAL(DP), DIMENSION(:), INTENT(IN)   :: T_LON, T_LAT, T_ALT
  REAL(DP),               INTENT(OUT)  :: FL_T(DIV_WAY+1), FL_T_Julian(DIV_WAY+1), FL_DIST(DIV_WAY+1)  !, FL_SPEED(DIV_WAY+1)
  REAL(DP),               INTENT(OUT)  :: ac_gs(DIV_WAY+1)  !20140523-1 ac ground speed including wind effects [m/s]
  REAL(DP),               INTENT(OUT)  :: ac_cpc(DIV_WAY+1)  !CPC [km] 20170602-1
  REAL(DP), INTENT(OUT)                :: T_FL_T, T_FL_T_Julian, T_FL_DIST
  !Yin_20170423
  REAL(DP), INTENT(OUT)                :: T_potcov           ![km] 20170607-1
  REAL(DP), INTENT(OUT),OPTIONAL       :: T_atr20cpc           ![K] 20170801
  !Yin_20170423
  REAL(DP), DIMENSION(:),ALLOCATABLE   :: temp_atr20cpc_value
  REAL(DP)                             :: ac_atr20cpc(DIV_WAY+1)  !ATR20CPC [K] 20170801
  REAL(DP)                             :: W1_lat, W1_lon, W1_alt, W2_lat, W2_lon, W2_alt 
  REAL(DP)                             :: RARAD1_traj, RARAD2_traj, DCRAD1_traj, DCRAD2_traj   !20140726-2
  REAL(DP)                             :: DELDC1_traj, DELRA1_traj, DELRA2_traj, COSIS_traj    
  REAL(DP)                             :: SINDIS_traj, SINDIS1_traj, SINDIS2_traj, DIS !,DISN,COSIS_DEGREE 
  REAL(DP)                             :: DELT_T, DELT_T_Julian 
! REAL(DP)                             :: near_lon, near_lat  !20140508-7
  INTEGER                              :: j1,nlen
  !Yin_20170423
  !INTEGER                             :: j2
  !Yin_20170423

  IF(PRESENT(atr20cpc_value))THEN
    nlen=SIZE(atr20cpc_value)
    ALLOCATE(temp_atr20cpc_value(nlen))
    temp_atr20cpc_value = atr20cpc_value
    T_atr20cpc      = 0.0_dp     !20170810[K]
  ENDIF
 !Initilizations
! AC_V            = 0.0_dp 
  FL_DIST(1)      = 0.0_dp 
  ac_cpc(1)       = 0.0_dp   !20170607-1
  ac_atr20cpc(1)  = 0.0_dp   !20170801
! FL_SPEED(DIV_WAY+1)  = 0.0_dp
  ac_gs(DIV_WAY+1)= 0.0_dp     !20140523-1
  T_FL_DIST       = 0.0_dp 
  DELT_T          = 0.0_dp
  DELT_T_Julian   = 0.0_dp     ![Julian date]
  T_FL_T          = 0.0_dp   
  T_FL_T_Julian   = 0.0_dp     ![Julian date]
  DIS             = 0.0_dp
  FL_T(1)         = 0.0_dp 
  FL_T_Julian(1)  = gc_D_time  ![Julian date]
  !Yin_20170423
  T_potcov        = 0.0_dp     !20161212,20170607-1 [km]
  !Yin_20170423

  !write(*,*)'gc_D_time_check:',gc_D_time,'DIV_WAY:',DIV_WAY,'DTR:',DTR,'ISWITCH',ISWITCH
  !write(*,*)'radius_earth:',radius_earth,'OneDay',OneDay  !'ALT',ALT 20140524-9
!  write(*,*)"wind_value_check:",wind_value
!  write(*,*)"ac_vtas_check",ac_vtas,"[m/s]"
  
  !Calculate flight details from D_city to A_city along the waypoints
  !Waypoints are already in order from D to A city.
  !The range of variables are:see opt20130215-4,20140512-4
  !    -180 =< T_LON =< 180 [FL_DIR=0,2,4,5]
  !       0 =< T_LON =< 360 [FL_DIR=1,3]
  !20120104-1
  DO j1=1, DIV_WAY   !Loop for waypoint
     W1_lat = T_LAT(j1)
     W1_lon = T_LON(j1)
     W1_alt = T_ALT(j1)
     W2_lat = T_LAT(j1+1)
     W2_lon = T_LON(j1+1)
     W2_alt = T_ALT(j1+1)
     !Yin_20170423
!     T_CPC  = T_CPC+cpc_value(j1+1)
     !Yin_20170423

     !Definition of the waypoint being placed at higher lon between W1 and W2
     IF(W1_lon.lt.W2_lon)THEN
        RARAD1_traj= W2_lon*DTR   !lon of higher waypoint [rad]
        RARAD2_traj= W1_lon*DTR   !lon of lower  waypoint [rad]
        DCRAD1_traj= W2_lat*DTR   !lat of higher waypoint [rad]    
        DCRAD2_traj= W1_lat*DTR   !lat of lower  waypoint [rad]
!       WP_FL_DIR= 0         !20140508-2
     ELSE IF(W1_lon.gt.W2_lon)THEN
        RARAD1_traj= W1_lon*DTR     
        RARAD2_traj= W2_lon*DTR
        DCRAD1_traj= W1_lat*DTR     
        DCRAD2_traj= W2_lat*DTR     
!       WP_FL_DIR= 1         !20140508-2
     ELSE                    !(W1_lon.eq.W2_lon)
        IF(W1_lat.gt.W2_lat)THEN
           RARAD1_traj= W1_lon*DTR     
           RARAD2_traj= W2_lon*DTR
           DCRAD1_traj= W1_lat*DTR     
           DCRAD2_traj= W2_lat*DTR     
!          WP_FL_DIR= 2      !20140508-2
        ELSE IF(W1_lat.lt.W2_lat)THEN
           RARAD1_traj= W2_lon*DTR 
           RARAD2_traj= W1_lon*DTR
           DCRAD1_traj= W2_lat*DTR     
           DCRAD2_traj= W1_lat*DTR 
!          WP_FL_DIR= 3      !20140508-2
        ELSE                 !(W1_lat.eq.W2_lat)
           write(*,*)'Waypoints Error: They must be identical waypoints!'
        ENDIF
     ENDIF

     !Calculate central angle between W1 and W2.
     !COSIS: central angle in [rad]
     SELECT CASE (ISWITCH)  
     CASE(0)
        DELRA1_traj  = RARAD1_traj-RARAD2_traj
        COSIS_traj  = ACOS(SIN(DCRAD1_traj)*SIN(DCRAD2_traj) + COS(DCRAD1_traj)*COS(DCRAD2_traj)*&
                      COS(DELRA1_traj))               !in [rad]
     CASE(1) 
        DELDC1_traj = (DCRAD1_traj-DCRAD2_traj)/2.0_dp            !delta lat
        DELRA2_traj = (RARAD1_traj-RARAD2_traj)/2.0_dp            !delta lon
        SINDIS_traj = SQRT(SIN(DELDC1_traj)*SIN(DELDC1_traj) + COS(DCRAD1_traj)*COS(DCRAD2_traj)*&
                      SIN(DELRA2_traj)*SIN(DELRA2_traj))
        COSIS_traj  = 2.0_dp*DASIN(SINDIS_traj)              !in [rad]
     CASE(2)
        DELRA2_traj = RARAD1_traj-RARAD2_traj                  !delta lon
        SINDIS1_traj= SQRT((COS(DCRAD2_traj)*SIN(DELRA2_traj))**2 + (COS(DCRAD1_traj)*SIN(DCRAD2_traj)-&
                       SIN(DCRAD1_traj)*COS(DCRAD2_traj)*COS(DELRA2_traj))**2)
        SINDIS2_traj= SIN(DCRAD1_traj)*SIN(DCRAD2_traj) + COS(DCRAD1_traj)*COS(DCRAD2_traj)*COS(DELRA2_traj) 
        COSIS_traj  = DATAN2(SINDIS1_traj,SINDIS2_traj)        !in [rad]
     CASE DEFAULT
        write(*,*) 'ISWITCH error: Unknown ISWITCH on WP distance calculation!'
     END SELECT

     !Calculate flight distance between W1 and W2, opt20121016-2 
     !DIS does not change with or without wind effect. Wind affects only ac ground speed. 

     !Same altitude at W1 and W2
!    DIS     = (radius_earth*0.001_dp + ALT/1000.0_dp)*COSIS      !WP distance = radius * radians [km] 
!    DIS     = (RADIUS + ALT/1000.0_dp)*COSIS                     !WP distance = radius * radians [km] 
     
     !Different altitude at W1 and W2     
     DIS = SQRT(((radius_earth*0.001_dp+W1_alt*0.001_dp)**2) + ((radius_earth*0.001_dp+W2_alt*0.001_dp)**2)  &
                -2.0_dp*(radius_earth*0.001_dp+W1_alt*0.001_dp)*(radius_earth*0.001_dp+W2_alt*0.001_dp)*COS(COSIS_traj))  ![km]
!    DISN    = DIS*0.539957_dp                                    !Unit conversion DIS into nautical miles [nm]
!    write(*,*)'pi=',pi,'W1_alt=',W1_alt,'W2_alt=',W2_alt,'COSIS=',COSIS,"COS(COSIS)=",COS(COSIS)  ! 'COSIS_DEGREE=',COSIS_DEGREE
!    write(*,*)'radius_earth*0.001d0+W1_alt/1000=',radius_earth*0.001d0+W1_alt/1000.0d0
!    write(*,*)'radius_earth*0.001d0+W2_alt/1000=',radius_earth*0.001d0+W2_alt/1000.0d0

     T_FL_DIST = T_FL_DIST + DIS                                  !Total flight distance from D_city to A_city [km]
     FL_DIST(j1+1) = DIS                                          !Flight distance between W1 and W2 [km]

!    write(*,*)'Distance |W',j1,'-W',j1+1,'| =',DIS,'[km]',DISN,'[km]',' Total flight distance =',T_FL_DIST,'[km]'
!    write(*,*)'ISWITCH =',ISWITCH,'FL_DIR =',FL_DIR

     !Calculate ac ground speed considering wind effects at each waypoints.
     !   ac_vtas    :AC true air speed at waypoint1 W1 [m/s]
     !   wind_value :wind values(u,v,w) at waypoint1 W1 [m/s]
     !   ac_gs      :AC ground speed at waypoint1 W1 [m/s]
     !CALL wind_component_lonlat(W1_lat, W1_lon, W1_alt, W2_lat, W2_lon, W2_alt, DIS, wind_value(j1,:), AC_V, new_AC_V) 
     CALL wind_effect(W1_lat, W1_lon, W1_alt, W2_lat, W2_lon, W2_alt, DIS, wind_value(j1,:), ac_vtas(j1), ac_gs(j1)) 
     !write(*,*)"AC_SPEED_CHECK","ac_vtas[m/s]=",ac_vtas(j1),"ac_gs[m/s]=",ac_gs(j1) 
 
     !AC_V   = AC_MACH*SOUND_V/1000.0_dp                        !Aircraft speed [km/s] 
!    FL_SPEED(j1) = AC_V  ! <<<<<<-- input wind fields,  AC_V[km/s], i.e., FL_SPEED[km/s] 
!HY  FL_SPEED(j1) = ac_gs*0.001_dp  ! FL_SPEED[km/s]
                                 
     !Calculate flight time between W1 and W2
     !DELT_T  : Flight time between W1 and W2 [s]  20140523-1
     !FL_DIST : Flight distance between W1 and W2 [km]
     !ac_gs   : AC ground speed[m/s] = ac_vtas + wind_speed 20140523-1
!HY  DELT_T = FL_DIST(j1+1)/FL_SPEED(j1)                 
     DELT_T = FL_DIST(j1+1)/(ac_gs(j1)*0.001_dp)           !20140523-1 
     FL_T(j1+1)  = FL_T(j1) + DELT_T                       !Aircraft passing time at each waypoint [s],see20130513-1
     if(j1.eq.DIV_WAY)T_FL_T = FL_T(j1+1)                  !Total flight time from D_city to A_city [s] 
!    write(*,*)'DELT_T_CHECK=',DELT_T,'j1=',j1,"FL_SPEED=",FL_SPEED(j1),"FL_DIST=",FL_DIST(j1+1),"FL_T=",FL_T(j1+1)
     !write(*,*)'j1=',j1,'W1=',W1_alt,'W2=',W2_alt,"FL_DIST=",FL_DIST(j1+1),'DELT_T=',DELT_T,"FL_T=",FL_T(j1+1) 

     !Flight time [Julian date]
     DELT_T_Julian = DELT_T/OneDay                         !Flight time between W1 and W2 [Julian date]
     FL_T_Julian(j1+1) = FL_T_Julian(j1) + DELT_T_Julian   !Aircraft passing time at each waypoint [Julian date]
     if(j1.eq.DIV_WAY)T_FL_T_Julian = FL_T_Julian(j1+1)    !Arrival time at A_city [Julian date] 

     !Fuel_use calculation see20130521-3
     ! FUEL_U(j1+1)=   *DELT_T                             !DELT_T[s]  
     ! FUEL_U(j1+1)=    DELT_T                             !DELT_T[s]  

     !Calculate Total CPC, 20170602-1,20170607-1
     !cpc_value  :potential contrail coverage at waypoint [-]
     !ac_cpc     :potential contrail coverage for each segment [km]
     !T_potcov   :total potential contrail coverage [km]
     ac_cpc(j1+1) = cpc_value(j1)*FL_DIST(j1+1)            !potcov between W1 and W2 [km]
     T_potcov = T_potcov + ac_cpc(j1+1)                    ![km]
     IF(PRESENT(atr20cpc_value))THEN
       ac_atr20cpc(j1+1) = ac_cpc(j1+1)*atr20cpc_value(j1)   !ATR20 of potcov between W1 and W2 [K]
       T_atr20cpc = T_atr20cpc + ac_atr20cpc(j1+1)           ![K]
     ENDIF
     !Initialization for next waypoints calculation
     DIS             = 0.0_dp
     DELT_T          = 0.0_dp   
     DELT_T_Julian   = 0.0_dp   
  ENDDO   !End loop for waypoint

  !Yin_20170423
!  DO j2=1,DIV_WAY+1
!     IF(j2.le.DIV_WAY)THEN
!        cpc_value_temp(j2) = cpc_value(j2)*(FL_DIST(j2)+FL_DIST(j2+1))/2
!     ELSEIF(j2.eq.DIV_WAY+1)THEN
!        cpc_value_temp(j2) = cpc_value(j2)*FL_DIST(j2)/2
!     ENDIF
!     T_CPC = T_CPC+cpc_value_temp(j2)
!  ENDDO
  !Yin_20170423

!HY 20170607-1
!  DO j2=1,DIV_WAY
!     cpc_value_temp(j2) = cpc_value(j2)*FL_DIST(j2+1)
!     T_CPC = T_CPC+cpc_value_temp(j2)
!  ENDDO
    IF (ALLOCATED(temp_atr20cpc_value)) DEALLOCATE(temp_atr20cpc_value)
  END SUBROUTINE gc_ac_traj
!----------------------------------------------------------------------

!----------------------------------------------------------------------
  SUBROUTINE fuel_calc(fo_ZZN, fo_wind_opt_fl_time, fo_wind_opt_vtas_along_traj,  &
                       fo_flight_time_julian, fo_rho_along_traj, fo_total_fuel_use, &
                       fo_total_nox_emis, &
                       fo_press_along_traj, fo_t_scb_along_traj,&
                       fo_ATR20O3_along_traj,fo_ATR20CH4_along_traj,fo_ATR20H2O_along_traj,&
                       fo_ATR20CO2_along_traj,&
                       fo_total_atr20o3,fo_total_atr20ch4,fo_total_atr20h2o,fo_total_atr20co2)
!----------------------------------------------------------------------
! This subroutine calculates Obj for fuel_opt routing option: total fuel usage for a flight
! using Total Energy model.20160603-1
!
!
  IMPLICIT NONE     !20160603-1
  REAL(DP), PARAMETER                       :: fo_extra_fuel_rate = 0.03_dp    !Percentage of extra fuel,20140724-2,20160605-1
 
  REAL(DP)                                  :: fo_wind_opt_ave_alt             !Average altitude of fuel opt traj. 20160603-1 
  REAl(DP)                                  :: fo_temp_wind_opt_ave_alt     !20160603-1
  REAL(DP)                                  :: fo_wind_opt_bada_fuel_nom                !Fuel from BADA A333_ptf,[kg/s] 20140524-7
  REAL(DP), DIMENSION(:), INTENT(IN)        :: fo_ZZN                    !altitude(= ac_routes(i,props_alt)), [m]
  REAL(DP), DIMENSION(:), INTENT(IN)        :: fo_flight_time_julian     !Passing time at each waypoint,[Julian]20160605-2
  REAL(DP), DIMENSION(:), INTENT(IN)        :: fo_rho_along_traj         !rho, [kg/m^3] 
  REAL(DP), DIMENSION(:), INTENT(IN),OPTIONAL :: fo_press_along_traj       !pressure, [Pa] 20170215-6,20170217-1
  REAL(DP), DIMENSION(:), INTENT(IN),OPTIONAL :: fo_t_scb_along_traj       !temperature, [K] 20170215-6,20170217-1
!Yin_20170801+
  REAL(DP), DIMENSION(:), INTENT(IN),OPTIONAL :: fo_ATR20O3_along_traj       !ATR20 zone, [K/kg(NO2)]
  REAL(DP), DIMENSION(:), INTENT(IN),OPTIONAL :: fo_ATR20CH4_along_traj       !ATR20 methane, [K/kg(NO2)]
  REAL(DP), DIMENSION(:), INTENT(IN),OPTIONAL :: fo_ATR20H2O_along_traj       !ATR20 h2o, [K/kg(fuel)]
  REAL(DP), DIMENSION(:), INTENT(IN),OPTIONAL :: fo_ATR20CO2_along_traj       !ATR20 co2, [K/kg(fuel)]
!Yin_20170801-
  REAL(DP)                                  :: fo_BADA_DATA(8,5) !this should be given by Aircraft table.20160603-1
  REAL(DP)    :: fo_VTAS(DIV_WAY+1)    !20140602-2,20160605-1,20160606-7
  REAL(DP)    :: fo_ac_mass,fo_ita,fo_Cf1,fo_Cf2,fo_Cfcr,fo_CD2,fo_S_area,fo_CD0,fo_new_ac_mass  !20160605-1 
  REAL(DP)                                  :: fo_arr_ac_mass, fo_ave_ac_mass       !20160605-1 
  REAL(DP), INTENT(OUT)                     :: fo_total_fuel_use             ![kg]20160606-1   
  REAL(DP), INTENT(OUT),OPTIONAL            :: fo_total_nox_emis             ![g(NO2)]20170215-6,20170217-1,4   
  REAL(DP), INTENT(OUT),OPTIONAL            :: fo_total_atr20o3              ![K] 20170801
  REAL(DP), INTENT(OUT),OPTIONAL            :: fo_total_atr20ch4             ![K] 20170801
  REAL(DP), INTENT(OUT),OPTIONAL            :: fo_total_atr20h2o             ![K] 20170801 
  REAL(DP), INTENT(OUT),OPTIONAL            :: fo_total_atr20co2             ![K] 20170801 
  REAL(DP), INTENT(IN)                      :: fo_wind_opt_fl_time          !20160605-1
  REAL(DP)                                  :: fo_ac_oew, fo_max_payload, fo_load_factor    !20160605-1
  REAL(DP)                                  :: fo_aa,fo_bb,fo_cc,fo_sq,fo_de,fo_x1,fo_x2
  REAL(DP), DIMENSION(:), INTENT(IN)        :: fo_wind_opt_vtas_along_traj      ![m/s] 20140522-4,20160605-1

  INTEGER                                   :: i,j,j1,k,j2,nlen_press,nlen_atr20
 
  !DLR Fuel flow method for NOx_opt 20170217-1,2
  REAL(DP)                                  :: fo_d_ac_routes(DIV_WAY+1)   !20170217-2,-3,20170221-2
  !Yin_20170801
  REAL(DP),DIMENSION(:),ALLOCATABLE         :: fo_temp_press_along_traj    !20170217-2,-3,20170221-2
  REAL(DP),DIMENSION(:),ALLOCATABLE         :: fo_temp_t_scb_along_traj    !20170217-2,-3,20170221-2
  REAL(DP),DIMENSION(:),ALLOCATABLE         :: fo_temp_ATR20O3_along_traj    !20170217-2,-3,20170221-2
  REAL(DP),DIMENSION(:),ALLOCATABLE         :: fo_temp_ATR20CH4_along_traj    !20170217-2,-3,20170221-2
  REAL(DP),DIMENSION(:),ALLOCATABLE         :: fo_temp_ATR20H2O_along_traj    !20170217-2,-3,20170221-2
  REAL(DP),DIMENSION(:),ALLOCATABLE         :: fo_temp_ATR20CO2_along_traj    !20170217-2,-3,20170221-2
  !Yin_20170801
  REAL(DP), PARAMETER                       :: fo_temp_sealv    = 288.15_dp             !Temperature at sea level[K],20140724-2
  INTEGER,PARAMETER                         :: fo_num_icao_data = 4
  INTEGER,PARAMETER                         :: fo_order_app_equ = 2
  REAL(DP)                                  :: fo_EINOx(DIV_WAY+1) &
       ,fo_Wfuel_ref(DIV_WAY+1),fo_EINOx_ref(DIV_WAY+1),fo_Wfuel(DIV_WAY+1)
  REAL(DP)                                  :: fo_ATR20O3(DIV_WAY+1) &
       ,fo_ATR20CH4(DIV_WAY+1),fo_ATR20H2O(DIV_WAY+1),fo_ATR20CO2(DIV_WAY+1)
  REAL(DP)                                  :: fo_coefficient(fo_order_app_equ+1),fo_z_array(fo_order_app_equ+1)
  REAL(DP)                                  :: fo_w_array(fo_order_app_equ+1,fo_order_app_equ+1)
  REAL(DP)                                  :: fo_ICAO_DATA(2,4) !this should be given by Aircraft table.20130803-1,20140723-1
  REAL(DP)                                  :: fo_pressure_total,fo_temperature_total,fo_num_engine,fo_theta_total,fo_delta_total
  REAL(DP)                                  :: fo_hmd_corr_factor,fo_specific_hmd !,pressure_actual,temperature_actual
  REAL(DP)                                  :: fo_s_chol

 
  !COMMON DATA 20170217-2
  !ICAO DATA input,20140723-1
  DATA ((fo_ICAO_DATA(i,j),i=1,2),j=1,4) /0.228_dp,4.88_dp,0.724_dp,12.66_dp,2.245_dp,22.01_dp,2.767_dp,28.72_dp/

  !BADA DATA input,A333_ptf,BADA_ISA_TABLE, 20131108-1, 20140524-7,20160603-1
  DATA ((fo_BADA_DATA(i,j),i=1,8),j=1,3) /280.0_dp,290.0_dp,310.0_dp,330.0_dp,350.0_dp,370.0_dp,390.0_dp,410.0_dp,         &  !FL[ft]   
                                          305.79_dp,304.48_dp,301.86_dp,299.21_dp,296.54_dp,295.07_dp,295.07_dp,295.07_dp, &  !a[m/s]  
                                          112.4_dp,108.7_dp,101.7_dp,95.5_dp,90.0_dp,85.5_dp,81.9_dp,79.0_dp/      !Fuel[nom, kg/min]
  do i=1,8   !20160603-1
     fo_BADA_DATA(i,4)=fo_BADA_DATA(i,1)*100.0_dp*0.30480_dp   !Altitude [m],20131108-1
     fo_BADA_DATA(i,5)=fo_BADA_DATA(i,3)/60.0_dp               !Fuel [NOM case, kg/s]
     !write(*,*)'fo_BADA_DATA_check',(fo_BADA_DATA(i,j),j=1,5)
  enddo

!FROM BADA&ICAO date (from d_ac_routes_desc? or other new Aircraft table structure?) 
!Specific aircraft/engine data should be included into new Aircraft table structure.
 fo_ac_oew          = 125100.0_dp       !=Mmin from BADA A333_opf,[kg],20131114-3,20160605-1
 fo_max_payload     = 47900.0_dp        !=Mpyld from BADA A333_opf,[kg],20131114-3,20160605-1
 fo_load_factor     = 0.62_dp         !=0.609,20130426-1  ; change to 0.62(ICAO 2008),20131126-1B,20160605-1
 fo_Cf1             = 0.61503_dp      ![kg/min/kN]
 fo_Cf2             = 919.03_dp       ![Kt]
 fo_CD2             = 0.031875_dp     ![-]
 fo_Cfcr            = 0.93655_dp      ![-]
 fo_num_engine      = 2.0_dp          !A330 2_engines 
 fo_S_area          = 361.6_dp        !from BADA A333_opf,[m^2]
 fo_CD0             = 0.019805_dp     ![-]


! Initialiyation 20160603-1
 fo_total_fuel_use         = 0.0_dp   !20160606-1, Obj value for fuel_opt routing option 
 fo_wind_opt_ave_alt       = 0.0_dp
 fo_wind_opt_bada_fuel_nom = 0.0_dp   !20140524-7
 fo_arr_ac_mass            = 0.0_dp   !20140926-3
 fo_ave_ac_mass            = 0.0_dp   !20141216-4
 fo_ita                    = 0.0_dp
 fo_aa                     = 0.0_dp
 fo_bb                     = 0.0_dp
 fo_cc                     = 0.0_dp
 fo_sq                     = 0.0_dp
 fo_de                     = 0.0_dp
 fo_x1                     = 0.0_dp
 fo_x2                     = 0.0_dp
 !Initialization for DLR Fuel Flow Method 20170217-2
 fo_Wfuel           = 0.0_dp
 fo_coefficient     = 0.0_dp         
 fo_d_ac_routes     = 0.0_dp
 fo_pressure_total  = 0.0_dp
 fo_temperature_total = 0.0_dp
 fo_hmd_corr_factor = 0.0_dp
 fo_specific_hmd    = 0.0_dp
 !ATR20 calculation initial 20170801
  !Yin_20170801+
  IF (present(fo_press_along_traj)) THEN
     nlen_press = SIZE(fo_press_along_traj)
     ALLOCATE(fo_temp_press_along_traj(nlen_press))
     fo_temp_press_along_traj = fo_press_along_traj
  ENDIF
  IF (present(fo_t_scb_along_traj)) THEN
     ALLOCATE(fo_temp_t_scb_along_traj(nlen_press))
     fo_temp_t_scb_along_traj = fo_t_scb_along_traj
     fo_EINOx           = 0.0_dp
     fo_EINOx_ref       = 0.0_dp
     fo_Wfuel_ref       = 0.0_dp
     fo_delta_total     = 0.0_dp
     fo_theta_total     = 0.0_dp
     fo_total_nox_emis  = 0.0_dp   !20170215-6,20170217-1 Obj value for NOx_opt routing option
  ENDIF

  IF (present(fo_ATR20O3_along_traj)) THEN
     nlen_atr20 = SIZE(fo_ATR20O3_along_traj)
     ALLOCATE(fo_temp_ATR20O3_along_traj(nlen_atr20))
     fo_temp_ATR20O3_along_traj = fo_ATR20O3_along_traj
     fo_ATR20O3         = 0.0_dp
     fo_total_atr20o3   = 0.0_dp
  ENDIF
  
  IF (present(fo_ATR20CH4_along_traj)) THEN
     ALLOCATE(fo_temp_ATR20CH4_along_traj(nlen_atr20))
     fo_temp_ATR20CH4_along_traj = fo_ATR20CH4_along_traj
     fo_ATR20CH4        = 0.0_dp
     fo_total_atr20ch4  = 0.0_dp
  ENDIF
  IF (present(fo_ATR20H2O_along_traj)) THEN
     ALLOCATE(fo_temp_ATR20H2O_along_traj(nlen_atr20))
     fo_temp_ATR20H2O_along_traj = fo_ATR20H2O_along_traj
     fo_ATR20H2O        = 0.0_dp
     fo_total_atr20h2o  = 0.0_dp
  ENDIF
  IF (present(fo_ATR20CO2_along_traj)) THEN
     ALLOCATE(fo_temp_ATR20CO2_along_traj(nlen_atr20))
     fo_temp_ATR20CO2_along_traj = fo_ATR20CO2_along_traj
     fo_ATR20CO2        = 0.0_dp
     fo_total_atr20co2  = 0.0_dp
  ENDIF
 
!Yin_20170801
 
 !Initial AC mass 1st rough estimation, by using BADA data.20160603-1
 !fo_wind_opt_bada_fuel_nom,[kg/s] 20140524-7,8
 !Calc fo_wind_opt_ave_alt.
 fo_temp_wind_opt_ave_alt = 0.0_dp
 do i=1, DIV_WAY+1             !w_nwaypoints,20160606-7
    fo_temp_wind_opt_ave_alt = fo_temp_wind_opt_ave_alt + fo_ZZN(i)  !sum altitudes 20140524-8 
 enddo
 fo_wind_opt_ave_alt = fo_temp_wind_opt_ave_alt/(DIV_WAY+1)   !average altitude [m],20160606-7 
! write(*,*)"fo_wind_opt_ave_alt=",fo_wind_opt_ave_alt,"[m]",DIV_WAY+1
!write(*,*)'fo_ZZN',fo_ZZN   !20170223-4

 !Calculate average fuel consumption at averaged fuel-opt-traj altitude, by interplating BADA data.20160603-1,2
 do i=1,7
    if((fo_BADA_DATA(i,4)<=fo_wind_opt_ave_alt).and.(fo_wind_opt_ave_alt<fo_BADA_DATA(i+1,4)))then
       fo_wind_opt_bada_fuel_nom=(fo_BADA_DATA(i,5)*(fo_wind_opt_ave_alt-fo_BADA_DATA(i+1,4))   &
                                 -fo_BADA_DATA(i+1,5)*(fo_wind_opt_ave_alt-fo_BADA_DATA(i,4)))/ &
                                 (fo_BADA_DATA(i,4)-fo_BADA_DATA(i+1,4))      ![kg/s]
!       write(*,*)"if_bada_fuel_4",i,fo_wind_opt_ave_alt,fo_BADA_DATA(i+1,4)
       exit
    elseif(fo_wind_opt_ave_alt.lt.fo_BADA_DATA(1,4))then
       fo_wind_opt_bada_fuel_nom=fo_BADA_DATA(1,5)
!       write(*,*)"if_bada_fuel_5",i,fo_wind_opt_ave_alt,fo_BADA_DATA(1,4)
       exit
    elseif(fo_wind_opt_ave_alt.ge.fo_BADA_DATA(8,4))then
       fo_wind_opt_bada_fuel_nom=fo_BADA_DATA(8,5)
!       write(*,*)"if_bada_fuel_6",i,fo_wind_opt_ave_alt,fo_BADA_DATA(8,4)
       exit
    endif
 enddo
! write(*,*)'(Fuel)fo_wind_opt_bada_fuel_nom=',fo_wind_opt_bada_fuel_nom,"[kg/s]"


 !fo_wind_opt_fl_time shows total flight time [s] of a current design candidate.
 !fo_wind_opt_vtas_along_traj(nwaypoints)[m/s] shows vtas of a current design candidate. 
 fo_ac_mass         = fo_wind_opt_fl_time*fo_wind_opt_bada_fuel_nom*fo_extra_fuel_rate  &   !20140724-2
                    + fo_ac_oew + fo_max_payload*fo_load_factor
 fo_VTAS            = fo_wind_opt_vtas_along_traj                 !20140524-9,20140602-2,20150306-2
! write(*,*)"(Fuel)fo_wind_opt_fl_time=",fo_wind_opt_fl_time,"[s]"
! write(*,*)"(Fuel)fo_wind_opt_vtas_along_traj=",fo_wind_opt_vtas_along_traj,"[m/s]"  !20140524-6,20170223-4
! write(*,*)'(Fuel)fo_itnitial_ac_mass=',fo_ac_mass,"[kg]"
! write(*,*)'(Fuel)fo_fuel_nom=',fo_wind_opt_bada_fuel_nom,"[kg/s]"
! write(*,*)'(Fuel)fo_VTAS=',fo_VTAS,"[m/s]"   !20170223-4

!--------------------------------------------------
!    Total Energy Model: calculate fuel flow 
!--------------------------------------------------
     !   - d_as_routes(:,props_fuel_use) is calculated here.
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

     !!see 20131031-1,20130228,20130115-1,20130314-1,20160605-1

 !20140926-3,20160605-1
 fo_arr_ac_mass = fo_ac_mass     !researve aircraft (initial) mass at A_city
 fo_ave_ac_mass = fo_ac_mass     !add ac_mass at A_city,20141216-4   
 !Backward calculation on ac_mass from A_city to D_city.
 DO j1=DIV_WAY+1,2,-1     !20160606-7
    fo_ita = fo_Cf1*(1.0_dp + fo_VTAS(j1-1)/0.514443004_dp/fo_Cf2)/60.0_dp*0.001_dp ![kg/s/N] 20130521-3,20140602-3,20140620-8,20160605-1
!   write(*,*)'j1=',j1,'fo_VTAS=',fo_VTAS(j1-1),'fo_ita=',fo_ita,'fo_ac_mass=',fo_ac_mass   !20170223-4
!   write(*,*)'fo_flight_time_julian(j1)',fo_flight_time_julian(j1),'fo_rho_along_traj(j1)',fo_rho_along_traj(j1)   !20170223-4
    !20130711-1,20130522-1
    fo_aa = fo_ita*fo_Cfcr*((fo_flight_time_julian(j1)-fo_flight_time_julian(j1-1))*OneDay)* &
            fo_CD2*2.0_dp*(g**2)/fo_VTAS(j1-1)**2/fo_rho_along_traj(j1-1)/fo_S_area   !20130902-1,20140602-3 
    fo_bb = -1.0_dp
    fo_cc = fo_ac_mass + fo_ita*fo_Cfcr*((fo_flight_time_julian(j1)-fo_flight_time_julian(j1-1))*OneDay)* &
            fo_rho_along_traj(j1-1)*(fo_VTAS(j1-1)**2)*fo_S_area*fo_CD0*0.5_dp        !20130902-1,20140602-3

    !20160605-2
    fo_sq = dsqrt(fo_bb*fo_bb-4.0_dp*fo_aa*fo_cc)
    fo_de = 2.0_dp*fo_aa
    fo_x1 = (-fo_bb+fo_sq)/fo_de
    fo_x2 = -(fo_bb+fo_sq)/fo_de
    !write(*,*)'fo_x1=',fo_x1,'fo_x2=',fo_x2

    !Determine a solution (fo_x1 or fo_x2) as fo_new_ac_mass.20160606-1
    IF(fo_sq.gt.0.0_dp)THEN    !D>0: Two solution exist
       IF(fo_x1.ge.0.0_dp .and. fo_x2.ge.0.0_dp)THEN
          IF(fo_x1.ge.fo_ac_mass .and. fo_x2.ge.fo_ac_mass)THEN
             IF(dabs(fo_x1-fo_ac_mass) .gt. dabs(fo_x2-fo_ac_mass))fo_new_ac_mass = fo_x2
             IF(dabs(fo_x1-fo_ac_mass) .lt. dabs(fo_x2-fo_ac_mass))fo_new_ac_mass = fo_x1
          ELSE
             IF(fo_x1.lt.fo_ac_mass .and. fo_x2.ge.fo_ac_mass)fo_new_ac_mass = fo_x2
             IF(fo_x2.lt.fo_ac_mass .and. fo_x1.ge.fo_ac_mass)fo_new_ac_mass = fo_x1
          ENDIF
       ELSEIF(fo_x1.gt.0.0_dp .and. fo_x2.lt.0.0_dp)THEN
          fo_new_ac_mass = fo_x1
       ELSEIF(fo_x1.lt.0.0_dp .and. fo_x2.gt.0.0_dp)THEN
          fo_new_ac_mass = fo_x2
       ELSE
          write(*,*)'ERROR: Total_energy_model_calculation(D>0)'
       ENDIF

    ELSEIF(fo_sq.eq.0.0_dp)THEN  !D=0: One solution exists(fo_x1=fo_x2)
       IF(fo_x1.gt.0.0_dp)fo_new_ac_mass = fo_x1
       IF(fo_x1.le.0.0_dp)write(*,*)'ERROR: Total_energy_model_calculation(D=0)'
    ENDIF

    !Researve Fuel_use value on current segment and exchange mass values 
    !for calculation at next waypoint.
    !d_ac_routes(j1,props_fuel_use) = fo_new_ac_mass - fo_ac_mass  ![kg],see20130711-1
    fo_d_ac_routes(j1) = fo_new_ac_mass - fo_ac_mass  ![kg],20170217-1,20170221-2
    fo_ac_mass     = fo_new_ac_mass                               ![kg]
    fo_ave_ac_mass = fo_ave_ac_mass + fo_ac_mass !see20160607-1   !20141216-4 calculate average value of ac_mass [kg]
!    write(*,*)'check_fo_ave_ac_mass',fo_ave_ac_mass, fo_ac_mass, fo_new_ac_mass
 ENDDO     !Backward calculation loop end
 !After this DO loop, fo_ac_mass shows the ac_mass at D_city.
 !Calculate total fuel use for fuel_opt, NOx_opt, H2O_opt and Cost_opt. 
 fo_total_fuel_use = fo_ac_mass - fo_arr_ac_mass   !Total fuel usage of a design candidate(a trajectory)[kg]20160607-1
! write(*,*)'fo_total_fuel_use', fo_ac_mass-fo_arr_ac_mass,'=', fo_total_fuel_use   !Total fuel usage[kg]20160606-1

! ------- added for NOx_opt 20170215-6 -------
! DIV_WAY+1 = w_nwaypoints (line 263) should be used. 20160606-7
! DLR fuel flow method described below is NOT applied for fuel_opt, H2O_opt and Cost_opt. 20170215-6
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

 IF((w_option_traj_calc_g.eq.3).or.(w_option_traj_calc_g.eq.8).or.(w_option_traj_calc_g.eq.9))THEN   !NOx_opt or ATR20_opt
!   write(*,*)"Check:routing_option_in_fuel_calc_1",w_option_traj_calc_g     !20180512-2           

     ! START: The Least Squares method (2nd order),     see20130803-2
     ! ------W_array and z_array calculations
!     write(*,*)'fo_ICAO_DATA=',fo_ICAO_DATA(:,:)
     fo_w_array(1,1)=fo_num_icao_data
     DO 10 j=2,fo_order_app_equ+1
        fo_w_array(j,1)=0.0_dp
        DO 20 k=1,fo_num_icao_data
           fo_w_array(j,1)=fo_w_array(j,1)+fo_ICAO_DATA(1,k)**(j-1)
        20 continue      
        fo_w_array(fo_order_app_equ+1,j)=0.0_dp
        DO 30 k=1,fo_num_icao_data
           fo_w_array(fo_order_app_equ+1,j)=fo_w_array(fo_order_app_equ+1,j)+fo_ICAO_DATA(1,k)**(j+fo_order_app_equ-1)
        30 continue
        DO 40 i=1,j-1
           fo_w_array(i,j-i+1)=fo_w_array(j,1)
        40 continue      
        DO 50 i=1,fo_order_app_equ+1-j
           fo_w_array(fo_order_app_equ+1-i,j+i)=fo_w_array(fo_order_app_equ+1,j)
        50 continue
     10 continue
     DO 60 j=1,fo_order_app_equ+1
        fo_z_array(j)=0.0_dp
        DO 70 k=1,fo_num_icao_data
           IF(j.eq.1)THEN
              fo_z_array(j)=fo_z_array(j)+fo_ICAO_DATA(2,k)
           ELSE
              fo_z_array(j)=fo_z_array(j)+(fo_ICAO_DATA(1,k)**(j-1))*fo_ICAO_DATA(2,k)
           ENDIF
        70 continue
     60 continue

     ! ------Wa=z calculation by Cholesky method
     ! START: The Cholesky method
     ! ------Trans(R)DR calculation
     !write(*,*)'coefficient_check0=',coefficient(:)
     DO 80 i=1,fo_order_app_equ+1
        fo_s_chol = 0.0_dp
        DO 90 k=1,i-1
           fo_s_chol=fo_s_chol+fo_w_array(k,i)*fo_w_array(k,i)*fo_w_array(k,k)
        90 continue
        fo_w_array(i,i)=fo_w_array(i,i)-fo_s_chol
        DO 110 j=i+1,fo_order_app_equ+1
           fo_s_chol=0.0_dp
           DO 100 k=1,i-1
              fo_s_chol=fo_s_chol+fo_w_array(k,i)*fo_w_array(k,j)*fo_w_array(k,k)
           100 continue
           fo_w_array(i,j)=(fo_w_array(i,j)-fo_s_chol)/fo_w_array(i,i)
        110 continue
     80 continue 
     ! ------Trans(R)y=b calculation
     DO 130 i=2,fo_order_app_equ+1
        fo_s_chol=0.0_dp
        DO 120 j=1,i-1
           fo_s_chol=fo_s_chol+fo_w_array(j,i)*fo_z_array(j)
        120 continue       
        fo_z_array(i)=fo_z_array(i)-fo_s_chol
     130 continue
     ! ------DRx=y calculation
     DO 150 i=fo_order_app_equ+1,1,-1
        !write(*,*)'coefficient_check1=',i,coefficient(:)
        fo_s_chol=0.0_dp
        DO 140 j=i+1,fo_order_app_equ+1
           !write(*,*)'coefficient_check2=',coefficient(j)
           fo_s_chol=fo_s_chol+fo_w_array(i,j)*fo_coefficient(j)
        140 continue
        fo_coefficient(i)=fo_z_array(i)/fo_w_array(i,i)-fo_s_chol
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
     DO j2=1,DIV_WAY     !20140603-1
        fo_Wfuel(j2) = fo_d_ac_routes(j2+1)/((fo_flight_time_julian(j2+1)-fo_flight_time_julian(j2))*OneDay)/fo_num_engine
        fo_specific_hmd    = (10.0_dp**(-3))*DEXP(-0.0001426_dp*(fo_ZZN(j2)*3.280839895013_dp-12900.0_dp))
        fo_hmd_corr_factor = -19.0_dp*(fo_specific_hmd-0.00634_dp)
        IF(present(fo_press_along_traj))then
          fo_pressure_total    = fo_temp_press_along_traj(j2)*(1.0_dp+0.2_dp*(AC_MACH**2))**3.5_dp !AC_MACH should be from Aircraft table 
          fo_temperature_total = fo_temp_t_scb_along_traj(j2)*(1.0_dp+0.2_dp*(AC_MACH**2))      !AC_MACH should be from Aircraft table
          fo_delta_total       = fo_pressure_total/atm2Pa
          fo_theta_total       = fo_temperature_total/fo_temp_sealv                         !20140724-2
          fo_Wfuel_ref(j2)     = fo_Wfuel(j2)*(fo_delta_total**(-1))*(fo_theta_total**(-0.5_dp))
          !EINOx_ref is calcalculated by using approximation fuction and Wfuel_ref.
          !EINOx_ref(j2)   = 15.0_dp  !dummy
          fo_EINOx_ref(j2)     = fo_coefficient(1) + &
               fo_coefficient(2)*fo_Wfuel_ref(j2) + fo_coefficient(3)*(fo_Wfuel_ref(j2)**2._dp)
          fo_EINOx(j2)       = fo_EINOx_ref(j2)*(fo_delta_total**0.4_dp)*(fo_theta_total**3._dp)*DEXP(fo_hmd_corr_factor) !*num_engine,20130802-2
!          write(*,*) 'fo_Wfuel_ref=',fo_Wfuel_ref(j2)      
!          write(*,*)'fo_EINOx_ref=',fo_EINOx_ref(j2),'fo_EINOx=',fo_EINOx(j2)     !see20131130-2
        ENDIF
!        write(*,*)'fo_num_engine=',fo_num_engine,'AC_MACH=',AC_MACH,'atm2Pa=',atm2Pa,'fo_temp_sealv=',fo_temp_sealv
!        write(*,*)'fo_Wfuel=',fo_Wfuel(j2),'fo_specific_hmd:',fo_specific_hmd 
        IF(present(fo_ATR20O3_along_traj))THEN
          fo_ATR20H2O(j2+1)     = fo_ATR20H2O_along_traj(j2)*fo_d_ac_routes(j2+1)        !20180511-2       
          fo_ATR20CO2(j2+1)     = fo_ATR20CO2_along_traj(j2)*fo_d_ac_routes(j2+1)        !20180511-2 
          fo_ATR20O3(j2+1)     = fo_ATR20O3_along_traj(j2)*fo_d_ac_routes(j2+1)*fo_EINOx(j2)/1000._dp      !20180511-2 
          fo_ATR20CH4(j2+1)     = fo_ATR20CH4_along_traj(j2)*fo_d_ac_routes(j2+1)*fo_EINOx(j2)/1000._dp      !20180511-2 
        ENDIF
!       if(d_p_pe.eq.0)write(21,*)d_i_traj_dep,EINOx(j2)      ![g(NOx)/kg(fuel)],20131130-3,20141215-6
     ENDDO
 
! IF((w_option_traj_calc_g.eq.3).or.(w_option_traj_calc_g.eq.8).or.(w_option_traj_calc_g.eq.9))THEN   !NOx_opt or ATR20_opt
!    write(*,*)"Check:routing_option_in_fuel_calc_1",w_option_traj_calc_g   !20180512-2         

    !Add Obj value for NOx_opt routing option.  
    DO i=1,DIV_WAY+1
       IF(i.eq.1)THEN
!          write(*,*)'Total_nox_emis at dep time(wp=1)',fo_total_nox_emis
!          write(*,*)'Total_atr20o3 at dep time(wp=1)',fo_total_atr20o3
!          write(*,*)'Total_atr20ch4 at dep time(wp=1)',fo_total_atr20ch4
!          write(*,*)'Total_atr20h2o at dep time(wp=1)',fo_total_atr20h2o
!          write(*,*)'Total_atr20co2 at dep time(wp=1)',fo_total_atr20co2
       ELSEIF(i.gt.1)THEN
          fo_total_nox_emis = fo_total_nox_emis + fo_d_ac_routes(i)*fo_EINOx(i-1)  !d_ac_routes(i,props_emis_nox),[g(NO2)],20180511-2
          fo_total_atr20o3  = fo_total_atr20o3 + fo_ATR20O3(i)  ! total atr20 ozone, [K]
          fo_total_atr20ch4 = fo_total_atr20ch4 + fo_ATR20CH4(i)  ! total atr20 methane, [K]
          fo_total_atr20h2o  = fo_total_atr20h2o + fo_ATR20H2O(i)  ! total atr20 h2o, [K]
          fo_total_atr20co2  = fo_total_atr20co2 + fo_ATR20CO2(i)  ! total atr20 co2, [K]
       ENDIF
    ENDDO

 !ELSEIF(w_option_traj_calc_g.eq.2 .or. w_option_traj_calc_g.eq.4 .or.  &
 !       w_option_traj_calc_g.eq.6 .or. w_option_traj_calc_g.eq.7 .or.  &
 !       w_option_traj_calc_g .eq. 10)THEN 
    !Fuel_opt=2, H2O_opt=4(20170222-5), Cost_opt=6(20170224-4) or COC_opt=7(20170316-5)
    !Add dummy value for fuel_opt, H2O_opt, Cost_opt and COC_opt.  
 !   fo_total_nox_emis = obj_dummy   !20170216-5
 !   fo_total_atr20o3 = obj_dummy   !20170216-5
 !   fo_total_atr20ch4 = obj_dummy   !20170216-5
 !   fo_total_atr20h2o = obj_dummy   !20170216-5
 !   write(*,*)"Check:routing_option_in_fuel_calc_2",w_option_traj_calc_g
 !   write(*,*)'fuel_calc:fo_total_nox_emis=',fo_total_nox_emis 
 ENDIF

  END SUBROUTINE fuel_calc
!----------------------------------------------------------------------

!----------------------------------------------------------------------
  SUBROUTINE atr20_calc(atr20_total_atr20o3,atr20_total_atr20ch4, atr20_total_atr20h2o,  &
                        atr20_total_atr20cpc,atr20_total_atr20co2,atr20_total)
!----------------------------------------------------------------------
! This subroutine calculates Obj for atr20_opt routing option: total fuel usage for a flight
! using Total Energy model.20160603-1
!
!
  IMPLICIT NONE     
  REAL(DP), INTENT(IN)      :: atr20_total_atr20o3
  REAL(DP), INTENT(IN)      :: atr20_total_atr20ch4
  REAL(DP), INTENT(IN)      :: atr20_total_atr20h2o
  REAL(DP), INTENT(IN)      :: atr20_total_atr20cpc
  REAL(DP), INTENT(IN)      :: atr20_total_atr20co2
  REAL(DP), INTENT(OUT)     :: atr20_total
 
  atr20_total = 0._dp
 ! atr20_total = atr20_total_atr20o3       !20180511-2
  atr20_total = atr20_total_atr20o3+atr20_total_atr20ch4+atr20_total_atr20h2o+   &
                atr20_total_atr20cpc+atr20_total_atr20co2
!  WRITE(*,*)'ATR20O3_tot=',atr20_total_atr20o3
!  WRITE(*,*)'ATR20CH4_tot=',atr20_total_atr20ch4
!  WRITE(*,*)'ATR20H2O_tot=',atr20_total_atr20h2o
!  WRITE(*,*)'ATR20CPC_tot=',atr20_total_atr20cpc
!  WRITE(*,*)'ATR20CO2_tot=',atr20_total_atr20co2
  
  END SUBROUTINE atr20_calc
!----------------------------------------------------------------------
!----------------------------------------------------------------------
  SUBROUTINE cost_atr20_calc(atr20_T_FLIGHT_TIME,atr20_Total_fuel_use,&
                            atr20_total_atr20o3,atr20_total_atr20ch4, atr20_total_atr20h2o,  &
                            atr20_total_atr20cpc,atr20_total_atr20co2,atr20_total_costclim)
!----------------------------------------------------------------------
! This subroutine calculates Obj for atr20_opt routing option: total fuel usage for a flight
! using Total Energy model.20160603-1
!
!
  IMPLICIT NONE    
   
  REAL(DP), INTENT(IN)      :: atr20_total_atr20o3
  REAL(DP), INTENT(IN)      :: atr20_total_atr20ch4
  REAL(DP), INTENT(IN)      :: atr20_total_atr20h2o
  REAL(DP), INTENT(IN)      :: atr20_total_atr20co2
  REAL(DP), INTENT(IN)      :: atr20_total_atr20cpc
  REAL(DP), INTENT(IN)      :: atr20_T_FLIGHT_TIME             ![s]
  REAL(DP), INTENT(IN)      :: atr20_Total_fuel_use            ![kg(fuel)]   
  REAL(DP), INTENT(OUT)     :: atr20_total_costclim
  REAL(DP), PARAMETER       :: Ct = 2710.0_dp                 !Unit time costs[US Dollar/h], 20170428-4
  REAL(DP), PARAMETER       :: Co = 0.0_dp                    !Other costs, 20170428-4
  REAL(DP), PARAMETER       :: alpha = 0.8_dp                 !weight, 20170505-1
  REAL(DP), PARAMETER       :: coeff_costclim = 100._dp        !coefficient, 20170505-1
  REAL(DP)                  :: Cf                             !Unit fuel costs[Cents/lbs], 20170428-4
  REAL(DP)                  :: atr20_total,atr20_Total_cost   !Total ATR20, Total cost

  
  atr20_total_costclim = 0._dp 

  atr20_total = atr20_total_atr20o3+atr20_total_atr20ch4+atr20_total_atr20h2o+&
                atr20_total_atr20cpc+atr20_total_atr20co2

  atr20_total = atr20_total_atr20o3

!  WRITE(*,*)'ATR20O3_tot=',atr20_total_atr20o3
!  WRITE(*,*)'ATR20CH4_tot=',atr20_total_atr20ch4
!  WRITE(*,*)'ATR20H2O_tot=',atr20_total_atr20h2o
!  WRITE(*,*)'ATR20CPC_tot=',atr20_total_atr20cpc
!  WRITE(*,*)'ATR20CO2_tot=',atr20_total_atr20co2
  
  Cf = Fuel_price*100.0_dp/Fuel_density                                       !Unit fuel costs[Cents/lbs] 20170224-5,20170309-3
  atr20_Total_cost = (Cf*0.01_dp)*atr20_Total_fuel_use*2.2046226218_dp + &      ![US Dollar]20170224-3,20170224-5,20170309-3
                     Ct*(atr20_T_FLIGHT_TIME/3600.0_dp) + Co 

!  write(*,*)'Fuel_price=',Fuel_price,'Fuel_density=',Fuel_density,'Cf=',Cf,'[Cents/Pound]'
!  write(*,*)'atr20_T_FLIGHT_TIME=',atr20_T_FLIGHT_TIME,'atr20_Total_fuel_use=',atr20_Total_fuel_use
!  write(*,*)'atr20_Total_cost=',atr20_Total_cost
!  write(*,*)"Cost cauclation (END)"

  atr20_total_costclim = atr20_Total_cost*(alpha)+(1._dp-alpha)*coeff_costclim*atr20_total
  END SUBROUTINE cost_atr20_calc
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
  SUBROUTINE contrail_cost_calc(costcpc_T_FLIGHT_TIME,costcpc_CPC,costcpc_Total_fuel_use,&
                                      costcpc_total_costcpc)
!----------------------------------------------------------------------
! This subroutine calculates Obj for atr20_opt routing option: total fuel usage for a flight
! using Total Energy model.20160603-1
!
!
  IMPLICIT NONE    
   
  REAL(DP), INTENT(IN)      :: costcpc_T_FLIGHT_TIME             ![s]
  REAL(DP), INTENT(IN)      :: costcpc_CPC
  REAL(DP), INTENT(IN)      :: costcpc_Total_fuel_use            ![kg(fuel)]   
  REAL(DP), INTENT(OUT)     :: costcpc_total_costcpc
  REAL(DP), PARAMETER       :: Ct = 2710.0_dp                 !Unit time costs[US Dollar/h], 20170428-4
  REAL(DP), PARAMETER       :: Co = 0.0_dp                    !Other costs, 20170428-4
  REAL(DP), PARAMETER       :: alpha = 0.8_dp                 !weight, 20170505-1
  REAL(DP), PARAMETER       :: coeff_costcpc = 100._dp        !coefficient, 20170505-1
  REAL(DP)                  :: Cf                             !Unit fuel costs[Cents/lbs], 20170428-4
  REAL(DP)                  :: costcpc_cost                   !Total ATR20, Total cost

  
  costcpc_total_costcpc = 0._dp 
  
  Cf = Fuel_price*100.0_dp/Fuel_density                                       !Unit fuel costs[Cents/lbs] 20170224-5,20170309-3
  costcpc_cost = (Cf*0.01_dp)*costcpc_Total_fuel_use*2.2046226218_dp + &      ![US Dollar]20170224-3,20170224-5,20170309-3
                     Ct*(costcpc_T_FLIGHT_TIME/3600.0_dp) + Co 

!  write(*,*)'Fuel_price=',Fuel_price,'Fuel_density=',Fuel_density,'Cf=',Cf,'[Cents/Pound]'
!  write(*,*)'costcpc_T_FLIGHT_TIME=',costcpc_T_FLIGHT_TIME,'costcpc_Total_fuel_use=',costcpc_Total_fuel_use
!  write(*,*)'costcpc_cost=',costcpc_cost
!  write(*,*)"Cost cauclation (END)"

  costcpc_total_costcpc = costcpc_cost*(alpha)+(1._dp-alpha)*coeff_costcpc*costcpc_CPC
  END SUBROUTINE contrail_cost_calc
!----------------------------------------------------------------------

!----------------------------------------------------------------------
  SUBROUTINE cost_calc(cost_T_FLIGHT_TIME, cost_Total_fuel_use, cost_Total_cost)
!----------------------------------------------------------------------
! This subroutine calculates Obj for cost_opt routing option. 20170428-4
! The idea is based on the Cost index.
!
!
  IMPLICIT NONE  
  REAL(DP), PARAMETER                       :: Ct = 2710.0_dp                 !Unit time costs[US Dollar/h], 20170428-4
  REAL(DP), PARAMETER                       :: Co = 0.0_dp                    !Other costs, 20170428-4
  REAL(DP), PARAMETER                       :: alpha = 1.0_dp                 !weight, 20170505-1
  REAL(DP)                                  :: Cf                             !Unit fuel costs[Cents/lbs], 20170428-4

  REAL(DP), INTENT(IN)                      :: cost_T_FLIGHT_TIME             ![s]
  REAL(DP), INTENT(IN)                      :: cost_Total_fuel_use            ![kg(fuel)]   
  REAL(DP), INTENT(OUT)                     :: cost_Total_cost                ![US Dollar], Obj for SOC_opt   

  !Initialization
  cost_Total_cost    = 0.0_dp

  Cf = Fuel_price*100.0_dp/Fuel_density                                       !Unit fuel costs[Cents/lbs] 20170224-5,20170309-3
  cost_Total_cost = (Cf*0.01_dp)*cost_Total_fuel_use*2.2046226218_dp + &      ![US Dollar]20170224-3,20170224-5,20170309-3
                     Ct*(cost_T_FLIGHT_TIME/3600.0_dp) + Co 

! This is for only Abt seminar, 20170505-1
!  cost_Total_cost = (1.0_dp - alpha)*((Cf*0.01_dp)*cost_Total_fuel_use*2.2046226218_dp) + & ![US Dollar]
!                    alpha*(Ct*(cost_T_FLIGHT_TIME/3600.0_dp)) + Co 
 
!  write(*,*)'Fuel_price=',Fuel_price,'Fuel_density=',Fuel_density,'Cf=',Cf,'[Cents/Pound]'
!  write(*,*)'cost_T_FLIGHT_TIME=',cost_T_FLIGHT_TIME,'cost_Total_fuel_use=',cost_Total_fuel_use
!  write(*,*)'cost_Total_cost=',cost_Total_cost
!  write(*,*)"Cost cauclation (END)"

  END SUBROUTINE cost_calc
!----------------------------------------------------------------------

!----------------------------------------------------------------------
  SUBROUTINE coc_calc(coc_T_FLIGHT_TIME, coc_Total_fuel_use, coc_Total_coc)
!----------------------------------------------------------------------
! This subroutine calculates Obj for coc_opt routing option. Total coc for a
! flight is caluclated using NASA Liebeck (1995) model. 20170316-7,20170302-1
!
! NOTE: Costs for Depreciation, Insurance, and Interest are not included for 
! the coc_Total_coc (Obj).
!
!
  IMPLICIT NONE                                                      !20170316-8
  REAL(DP), PARAMETER                       :: MTOGW = 467379.0_dp   ![lb] Maximum take-off weight, A330-301
  REAL(DP), PARAMETER                       :: AFW   = 227238.19_dp  ![lb] Estimated airframe weight, A330-301 
  REAL(DP), PARAMETER                       :: SLST  = 60400.0_dp    ![lbf] Thrust per engine, A330-301
! REAL(DP), PARAMETER                       :: MLGW = 383603.0_dp    ![lb] Maximum landing gross weight, A330-300 
  REAL(DP), PARAMETER                       :: num_seats = 295.0_dp  !Number of seats, A330-301,20180514-4
  REAL(DP), PARAMETER                       :: num_eng   = 2.0_dp    !Number of engines, A330-301 

  REAL(DP), PARAMETER                       :: inflation_rate   = 2.28_dp     !Ave. US inflation rate [%] (1994-2015),20170315-3  
  INTEGER, PARAMETER                        :: cur_year   = 2015              !Current year
  INTEGER, PARAMETER                        :: ref_year   = 1993              !Reference year by Liebeck method 

  REAL(DP), INTENT(IN)                      :: coc_T_FLIGHT_TIME              ![s]
  REAL(DP), INTENT(IN)                      :: coc_Total_fuel_use             ![kg(fuel)]   
  REAL(DP), INTENT(OUT)                     :: coc_Total_coc                  ![US Dollar], Obj for COC_opt   

  REAL(DP)                                  :: cost_FC, cost_CC,  cost_LF,  cost_NF
  REAL(DP)                                  :: cost_MA, cost_AML, cost_AMM, cost_AAMB
  REAL(DP)                                  :: cost_ME, cost_EML, cost_EMM, cost_EAMB
  REAL(DP)                                  :: cost_fuel
  REAL(DP)                                  :: temp_coc_Total_coc 
  REAL(DP)                                  :: AFLAB_MMH_FH, AFLAB_MMH_FC, AFLAB_MMH_TRIP 
  REAL(DP)                                  :: AFMAT_MAT_FH, AFMAT_MAT_FC, AFMAT_MAT_TRIP 
  REAL(DP)                                  :: ENGLAB_MMH_TRIP, ENGLAB_MAT_TRIP 

 !Initialization,20170316-8 
 cost_FC            = 0.0_dp
 cost_CC            = 0.0_dp
 cost_LF            = 0.0_dp
 cost_NF            = 0.0_dp
 cost_MA            = 0.0_dp
 cost_AML           = 0.0_dp
 cost_AMM           = 0.0_dp
 cost_AAMB          = 0.0_dp
 cost_ME            = 0.0_dp
 cost_EML           = 0.0_dp
 cost_EMM           = 0.0_dp
 cost_EAMB          = 0.0_dp
 cost_fuel          = 0.0_dp
 coc_Total_coc      = 0.0_dp          !Obj value for coc_opt routing option 
 temp_coc_Total_coc = 0.0_dp      

 !Flight Crew(cost_FC), 20170316-8
 !Domestic
 !cost_FC = (coc_T_FLIGHT_TIME/3600.0_dp)*(440_dp + 0.532_dp*(MTOGW/1000.0_dp))
 !Intetnational
 cost_FC = (coc_T_FLIGHT_TIME/3600.0_dp)*(482_dp + 0.590_dp*(MTOGW/1000.0_dp))

 !Cabin Crew(cost_CC), 20170316-8
 !Domestic
 !cost_CC = (coc_T_FLIGHT_TIME/3600.0_dp)*(num_seats/35.0_dp)*60.0_dp 
 !Intetnational
 cost_CC = (coc_T_FLIGHT_TIME/3600.0_dp)*(num_seats/30.0_dp)*78.0_dp 

 !Landing Fees(cost_LF), 20170316-8
 !Domestic
 !cost_LF = 1.50_dp*(MLGW/1000.0_dp) 
 !Intetnational
 cost_LF = 4.25_dp*(MTOGW/1000.0_dp)

 !Navigation Fees(cost_NF), 20170316-8
 !Domestic
 !cost_NF = 0.0_dp 
 !Intetnational
 cost_NF = 0.136_dp*500.0_dp*SQRT(MTOGW/1000.0_dp) 

 !Maintenance - Airframe(cost_MA),20170316-9
 !(a)Airframe Maintenance Labor(cost_AML)
 AFLAB_MMH_FH = 1.260_dp + (1.774_dp*AFW/100000.0_dp) - 0.1071_dp*((AFW/100000.0_dp)**2.0)
 AFLAB_MMH_FC = 1.614_dp + (0.7227_dp*AFW/100000.0_dp) + 0.1024_dp*((AFW/100000.0_dp)**2.0) 
 AFLAB_MMH_TRIP = AFLAB_MMH_FH*(coc_T_FLIGHT_TIME/3600.0_dp) + AFLAB_MMH_FC
 cost_AML = AFLAB_MMH_TRIP*25.0_dp

 !(b)Airframe Maintenance Materials(cost_AMM)
 AFMAT_MAT_FH = 12.39_dp + (29.80_dp*AFW/100000.0_dp) + 0.1806_dp*((AFW/100000.0_dp)**2.0)
 AFMAT_MAT_FC = 15.20_dp + (97.33_dp*AFW/100000.0_dp) - 2.862_dp*((AFW/100000.0_dp)**2.0)
 AFMAT_MAT_TRIP = AFMAT_MAT_FH*(coc_T_FLIGHT_TIME/3600.0_dp) + AFMAT_MAT_FC
 cost_AMM = AFMAT_MAT_TRIP

 !(c)Airframe Applied Maintenance Burden(cost_AAMB)
 cost_AAMB = 2.0_dp*cost_AML

 !sum of (a)-(c):
 cost_MA = cost_AML + cost_AMM + cost_AAMB

 !Maintenance - Engines(cost_ME),20170316-9
 !(a)Engine Maintenance Labor(cost_EML)
 ENGLAB_MMH_TRIP = (0.645_dp + (0.05*SLST/10000.0_dp))*(0.566_dp + 0.434_dp/(coc_T_FLIGHT_TIME/3600.0_dp))* &
                   (coc_T_FLIGHT_TIME/3600.0_dp)*num_eng
 cost_EML = ENGLAB_MMH_TRIP*25.0_dp 

 !(b)Engine Maintenance Materials(cost_EMM)
 ENGLAB_MAT_TRIP = (25.0_dp + (18.0_dp*SLST/10000.0_dp))*(0.62_dp + 0.38_dp/(coc_T_FLIGHT_TIME/3600.0_dp))* & 
                   (coc_T_FLIGHT_TIME/3600.0_dp)*num_eng
 cost_EMM = ENGLAB_MAT_TRIP

 !(c)Engine Applied Maintenance Burden(cost_EAMB)
 cost_EAMB = 2.0_dp*cost_EML

 !sum of (a)-(c):
 cost_ME = cost_EML + cost_EMM + cost_EAMB 

 !Fuel(cost_fuel):Liebeck methods is not used. The current fuel price is directly used to calculate cost_fuel. 20170316-10
 !Domestic
 !cost_fuel = (0.65_dp/6.7_dp)*coc_Total_fuel_use*2.2046226218_dp
 !Intetnational
 !cost_fuel = (0.70_dp/6.7_dp)*coc_Total_fuel_use*2.2046226218_dp  
 !
 cost_fuel = (Fuel_price/Fuel_density)*coc_Total_fuel_use*2.2046226218_dp   ![US Dollar] 20170224-5,20170309-3        
! write(*,*)"Fuel_price=", Fuel_price, "Fuel_density=", Fuel_density, "coc_Total_fuel_use=", coc_Total_fuel_use

 !The costs (except for fuel_cost) are calculated based on the price in 1993. 20170316-10 
 !Thus, the costs are scaled by inflation rate. cost_fuel is already calculated based on a present price.
 !Total COC costs:
 temp_coc_Total_coc = (cost_FC + cost_CC + cost_LF + cost_NF + cost_MA + cost_ME)*  &   !20170310-2,20170310-4,20170315-1,3
                      ((1.0_dp + inflation_rate*0.01_dp)**(cur_year-ref_year)) 
 coc_Total_coc = temp_coc_Total_coc + cost_fuel 
! write(*,*)"COC cauclation (END)"

  END SUBROUTINE coc_calc
!----------------------------------------------------------------------

!----------------------------------------------------------------------
  SUBROUTINE determine_lon_dvs(pos_dv7, pos_dv8, pos_dv9, pos_dv10, pos_dv11)
!----------------------------------------------------------------------
!
! This subroutine calculate positions of dv7 to dv11 in longitude.
! In the case of FL_DIR=4 and 5, this subroutine calculates positions in latitude
! of dv7 to dv11. 
!
! NOTE: Inputs from AirTraf,  -180 =< gc_D_lon(gc_A_lon) =< 180 [For all FL_DIRs]
!
!       Output: -180 =< pos_dv7 to pos_dv11 =< 180 [FL_DIR=0,2]
!                  0 =< pos_dv7 to pos_dv11 =< 360 [FL_DIR=1,3]
!                -90 =< pos_dv7 to pos_dv11 =< +90 [FL_DIR=4,5] (those variables mean position in latitudes)
!
!       After this subroutine, the range of gc_D_lon(gc_A_lon) is not changed:
!               -180 =< gc_D_lon(gc_A_lon) =< 180 [For all FL_DIRs]
!                     
  IMPLICIT NONE
! REAL(DP), INTENT(OUT)         :: lon_dv7, lon_dv8, lon_dv9, lon_dv10, lon_dv11  !20140724-4
  REAL(DP), INTENT(OUT)         :: pos_dv7, pos_dv8, pos_dv9, pos_dv10, pos_dv11
  REAL(DP)                      :: DIV_DV_ALT
! DOUBLE PRECISION              :: temp_gc_A_lon        !For FL_DIR=3, 0 =< temp_gc_A_lon =< 360
! DOUBLE PRECISION              :: temp_gc_D_lon        !For FL_DIR=1, 0 =< temp_gc_D_lon =< 360

  !Initializations
  DIV_DV_ALT    = 0.0_dp
  pos_dv7       = 0.0_dp
  pos_dv8       = 0.0_dp
  pos_dv9       = 0.0_dp
  pos_dv10      = 0.0_dp
  pos_dv11      = 0.0_dp
! temp_gc_A_lon = 0.0_dp
! temp_gc_D_lon = 0.0_dp
! write(*,*)'check FL_DIR=',FL_DIR,NUM_CP_ALT,gc_D_lon,gc_A_lon

  !Definitions: pos_dv7, pos_dv8, pos_dv9, pos_dv10, pos_dv11  20140724-4
  SELECT CASE(FL_DIR)
  CASE(0,3) !opt20130206-2
     IF(FL_DIR.eq.3)THEN
!       temp_gc_A_lon = gc_A_lon + 360.0d0 
        DIV_DV_ALT    = (temp_gc_A_lon - gc_D_lon)/(NUM_CP_ALT+1)
     ELSEIF(FL_DIR.eq.0)THEN
        DIV_DV_ALT    = (gc_A_lon - gc_D_lon)/(NUM_CP_ALT+1)
     ENDIF
     pos_dv7    =  gc_D_lon +  DIV_DV_ALT
     pos_dv8    =  gc_D_lon + (DIV_DV_ALT*2.0_dp)
     pos_dv9    =  gc_D_lon + (DIV_DV_ALT*3.0_dp)
     pos_dv10   =  gc_D_lon + (DIV_DV_ALT*4.0_dp)
     pos_dv11   =  gc_D_lon + (DIV_DV_ALT*5.0_dp)
  CASE(1) !opt20130206-2
!    temp_gc_D_lon = gc_D_lon + 360.0d0 
     DIV_DV_ALT    = (temp_gc_D_lon - gc_A_lon)/(NUM_CP_ALT+1)
     pos_dv7    =  temp_gc_D_lon -  DIV_DV_ALT
     pos_dv8    =  temp_gc_D_lon - (DIV_DV_ALT*2.0_dp)
     pos_dv9    =  temp_gc_D_lon - (DIV_DV_ALT*3.0_dp)
     pos_dv10   =  temp_gc_D_lon - (DIV_DV_ALT*4.0_dp)
     pos_dv11   =  temp_gc_D_lon - (DIV_DV_ALT*5.0_dp)
  CASE(2) !opt20130206-1
     DIV_DV_ALT = (gc_D_lon - gc_A_lon)/(NUM_CP_ALT+1)
     pos_dv7    =  gc_D_lon -  DIV_DV_ALT
     pos_dv8    =  gc_D_lon - (DIV_DV_ALT*2.0_dp)
     pos_dv9    =  gc_D_lon - (DIV_DV_ALT*3.0_dp)
     pos_dv10   =  gc_D_lon - (DIV_DV_ALT*4.0_dp)
     pos_dv11   =  gc_D_lon - (DIV_DV_ALT*5.0_dp)
  CASE(4) !opt20130206-3,20140724-3
     DIV_DV_ALT = (gc_D_lat - gc_A_lat)/(NUM_CP_ALT+1)
     pos_dv7    =  gc_D_lat -  DIV_DV_ALT
     pos_dv8    =  gc_D_lat - (DIV_DV_ALT*2.0_dp)
     pos_dv9    =  gc_D_lat - (DIV_DV_ALT*3.0_dp)
     pos_dv10   =  gc_D_lat - (DIV_DV_ALT*4.0_dp)
     pos_dv11   =  gc_D_lat - (DIV_DV_ALT*5.0_dp)
  CASE(5) !opt20130206-3,20140724-3 
     DIV_DV_ALT = (gc_A_lat - gc_D_lat)/(NUM_CP_ALT+1)
     pos_dv7    =  gc_D_lat +  DIV_DV_ALT
     pos_dv8    =  gc_D_lat + (DIV_DV_ALT*2.0_dp)
     pos_dv9    =  gc_D_lat + (DIV_DV_ALT*3.0_dp)
     pos_dv10   =  gc_D_lat + (DIV_DV_ALT*4.0_dp)
     pos_dv11   =  gc_D_lat + (DIV_DV_ALT*5.0_dp)

  CASE DEFAULT
     write(*,*) 'Dv_calc error1: this aircraft is not assinged of any FL_DIR!'
  END SELECT
! write(*,*)'FL_DIR=',FL_DIR,'pos_dv7=',pos_dv7,'pos_dv8=',pos_dv8,'pos_dv9=',pos_dv9,'pos_dv10=',pos_dv10,'pos_dv11=',pos_dv11

  END SUBROUTINE determine_lon_dvs
!----------------------------------------------------------------------

!---------------------------------------------------------------
 SUBROUTINE B_SPLINE(i_case, nd, nin, xin, yin, nout, xout, yout)  !HY20140424-1
!---------------------------------------------------------------
 IMPLICIT NONE
!INTEGER, PARAMETER                              :: nn=50
 INTEGER, INTENT(IN)                             :: nd, nin, i_case
 INTEGER, INTENT(OUT)                            :: nout
 REAL(DP),               INTENT(IN)              :: xin(nin), yin(nin)
!HY20140424-2 REAL(DP), DIMENSION(:), INTENT(OUT):: xout(nd),yout(nd)
 REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT):: xout, yout
 REAL(DP)                                        :: xp(0:nin+1),yp(0:nin+1)
 REAL(DP)                                        :: ds, s, q0, q1, q2, q3, x, y
 INTEGER                                         :: i, j, k, nn, nd2

!HY20140424-1,2, opt20121019-4
!Interpolation in (lon,lat).
 IF(i_case.eq.1)THEN            
!HY SELECT CASE(i_case)
!HY CASE(1)
    nn = INT((nd-1)/(nin-1))        !"nout" becomes DIV_WAY+1(nwaypoints)
    ALLOCATE(xout(nd),yout(nd))     !HY20140424-2

!Interpolation in (lon,alt). Instead of using "nd", "nd2" is calculated. 
!The nd2 is the number of data points on (lon,alt) interplated Bspline line.
!Currently, "nn" is given by designer.see20140424-2,20140424-6 
 ELSEIF(i_case.eq.2)THEN         
!HY CASE(2)
!   nn = 2*INT((nd-1)/(nin-1))      !"nout" is independent from DIV_WAY+1(nwaypoints), 20140424-2
    nn = 20  !10!2   !HY20140424-6,20140706-2,20141211-1,2 !"nout" becomes 20(nn)*6(NUM_CP_ALT-1) +1 =121, see line1311 below. 
    nd2= (nin-1)*nn+1               !Thus, we allocate "xout" and "yout" as same size.
    ALLOCATE(xout(nd2),yout(nd2)) 
 
 !For wind components calcuation
!HY CASE(3)  
!HY    nn = 5*INT((nd-1)/(nin-1))      !HY20140506-3
!HY    nd2= (nin-1)*nn+1               !Thus, we allocate "xout" and "yout" as same size.
!HY    ALLOCATE(xout(nd2),yout(nd2))   !
 
 ENDIF
!HY ENDSELECT

 !write(*,*)'nn=',nn
 !write(*,*)'nin=',nin
 !write(*,*)'nd=',nd,'nd2=',nd2
 !write(*,*)'xin=',xin,'yin=',yin

 do i=1,nin
    xp(i) = xin(i)
    yp(i) = yin(i)
 enddo
! write(*,*)'xp(0),yp(0)',xp(0),yp(0)
! write(*,*)'xp(nin+1),yp(nin+1)',xp(nin+1),yp(nin+1)

 xp(0)     = 2.0_dp*xp(1)-xp(2)
 yp(0)     = 2.0_dp*yp(1)-yp(2)
 xp(nin+1) = 2.0_dp*xp(nin)-xp(nin-1)
 yp(nin+1) = 2.0_dp*yp(nin)-yp(nin-1)

! write(*,*)'xp(0),yp(0)',xp(0),yp(0)
! write(*,*)'xp(nin+1),yp(nin+1)',xp(nin+1),yp(nin+1)

 ds=1.0_dp/DBLE(nn)
 do i=1, nin-1
    do j=1, nn
       s  = ds*DBLE(j-1)
       q0 = (s**3)/6.0_dp
       q1 = (-3.0_dp*s**3 +3.0_dp*s**2 + 3.0_dp*s + 1.0_dp)/6.0_dp
       q2 = ( 3.0_dp*s**3 -6.0_dp*s**2            + 4.0_dp)/6.0_dp
       q3 = (-1.0_dp*s**3 +3.0_dp*s**2 - 3.0_dp*s + 1.0_dp)/6.0_dp
       x  = xp(i+2)*q0+xp(i+1)*q1+xp(i)*q2+xp(i-1)*q3
       y  = yp(i+2)*q0+yp(i+1)*q1+yp(i)*q2+yp(i-1)*q3
       k  = (i-1)*nn+j

       xout(k) = x
       yout(k) = y
!      write(*,*)k,xout(k),yout(k)
    enddo
 enddo
 nout = (nin-1)*nn+1
 xout(nout) = xp(nin)
 yout(nout) = yp(nin)

! write(*,*)'-------------------------'
! do k=1,nout-1
!    write(*,*)k,xout(k),yout(k) 
! enddo
! write(*,*)nout,xout(nout),yout(nout)

 RETURN
!---------------------------------------------------------------
 END SUBROUTINE B_SPLINE
!---------------------------------------------------------------

!----------------------------------------------------------------------
  SUBROUTINE wind_effect(W1_lat_in, W1_lon_in, W1_alt_in, W2_lat_in, W2_lon_in, W2_alt_in, DIS_in, &
                                   wind_value_in, AC_V_in, new_AC_V_in)
!----------------------------------------------------------------------
!
! This subroutine calculate waypoints along GC.
! 20140508-3 
  IMPLICIT NONE
  INTEGER, PARAMETER                  :: DIV_WAY_in = 20  !20140508-6
  REAL(DP), DIMENSION(:),INTENT(IN)   :: wind_value_in  !20140509-1[m/s] 
  REAL(DP), INTENT(IN)                :: AC_V_in        !ac_vtas at W1 [m/s], 20140509-2,20140522-4   
  REAL(DP), INTENT(OUT)               :: new_AC_V_in    !ac_ground speed including wind effect [m/s], 20140509-2   
  REAL(DP)                            :: THETA_in(DIV_WAY_in+1), PHI_in(DIV_WAY_in+1) 
  REAL(DP)                            :: THETA_OUT_in(DIV_WAY_in+1), PHI_OUT_in(DIV_WAY_in+1) 
! REAL(DP), INTENT(OUT)               :: near_lon_in, near_lat_in   !20140509-1 
  REAL(DP)                            :: THETAR_in, PHIR_in
  REAL(DP)                            :: DELDC3_in, DELRA3_in
! REAL(DP)                            :: org_W1_lon, org_W2_lon 
  REAL(DP)                            :: RARAD1_in, RARAD2_in, DCRAD1_in, DCRAD2_in
! REAL(DP), INTENT(IN)                :: W1_lat_in, W1_lon_in, W2_lat_in, W2_lon_in
  REAL(DP)                            :: W1_lat_in, W1_lon_in, W1_alt_in, W2_lat_in, W2_lon_in, W2_alt_in
  REAL(DP)                            :: theta_wind, ux_wind_component, vy_wind_component, AC_V_temp  !,theta_wind2 20140512-2
  REAL(DP)                            :: delta_wind, DIS_in, wz_wind_component    !20140514-1
  INTEGER                             :: i, WP_FL_DIR

! The range of variables are:see opt20130215-4,20140512-4
!    -180 =< W1_lon_in,W2_lon_in =< 180 [FL_DIR=0,2,4,5]
!       0 =< W1_lon_in,W2_lon_in =< 360 [FL_DIR=1,3]

!Calculate nearest innter point of (W1_lon, W1_lat) on GC 
  !Initialization
  DELDC3_in = 0.0_dp
  DELRA3_in = 0.0_dp
  THETA_in  = 0.0_dp
  PHI_in    = 0.0_dp

  !To use original subroutine gc_waypoints, input range of LON should be as [-180 =< LON =< +180].20140508-5 
  !FL_DIR=1,3, these if sentences should be used.20140512-4
  IF(W1_lon_in.gt.180.0_dp)W1_lon_in = W1_lon_in-360.0_dp   ![deg]
  IF(W2_lon_in.gt.180.0_dp)W2_lon_in = W2_lon_in-360.0_dp   ![deg]

! HERE!! The range of variables are:see opt20130215-4,20140512-4
!    -180 =< W1_lon_in,W2_lon_in =< 180 [FL_DIR=ALL]

  !Definition of the waypoint being placed at higher lon. 
  IF(W1_lon_in.lt.W2_lon_in)THEN   
     RARAD1_in= W2_lon_in*DTR        !lon of higher waypoint [rad]
     RARAD2_in= W1_lon_in*DTR        !lon of lower  waypoint [rad]
     DCRAD1_in= W2_lat_in*DTR        !lat of higher waypoint [rad]    
     DCRAD2_in= W1_lat_in*DTR        !lat of lower  waypoint [rad]
     IF((RARAD1_in-RARAD2_in).le.PI)THEN
        WP_FL_DIR= 0
     ELSE
        WP_FL_DIR= 1
     ENDIF
  ELSE IF(W1_lon_in.gt.W2_lon_in)THEN
     RARAD1_in= W1_lon_in*DTR     
     RARAD2_in= W2_lon_in*DTR
     DCRAD1_in= W1_lat_in*DTR     
     DCRAD2_in= W2_lat_in*DTR     
     IF((RARAD1_in-RARAD2_in).le.PI)THEN 
        WP_FL_DIR= 2
     ELSE
        WP_FL_DIR= 3
     ENDIF
  ELSE     !(W1_lon_in.eq.W2_lon_in)
     IF(W1_lat_in.gt.W2_lat_in)THEN
        RARAD1_in= W1_lon_in*DTR     
        RARAD2_in= W2_lon_in*DTR
        DCRAD1_in= W1_lat_in*DTR     
        DCRAD2_in= W2_lat_in*DTR     
        WP_FL_DIR= 4
     ELSE IF(W1_lat_in.lt.W2_lat_in)THEN
        RARAD1_in= W2_lon_in*DTR 
        RARAD2_in= W1_lon_in*DTR
        DCRAD1_in= W2_lat_in*DTR     
        DCRAD2_in= W1_lat_in*DTR 
        WP_FL_DIR= 5
     ELSE     !(W1_lat_in.eq.W2_lat_in)
        write(*,*)'Waypoints Error: They must be identical waypoint!'
     ENDIF
  ENDIF

  !The following code is identical to the subroutine "gc_waypoints"
  DELRA3_in = RARAD1_in - RARAD2_in     !Delta lon between W1_lon & W2_lon
  DELDC3_in = DCRAD1_in - DCRAD2_in     !Delta lat between W1_lat & W2_lat

  !WP_FL_DIR=0,2  Regarding WP_FL_DIR, see20140508-6
  !On calculation, see20140507-5
  IF((DELRA3_in.gt.0.0_dp).and.(DELRA3_in.le.PI))THEN
     DO i=1, DIV_WAY_in+1
        PHIR_in = RARAD2_in + (i-1)*(RARAD1_in-RARAD2_in)/DIV_WAY_in              !in [rad]
        THETAR_in = atan(sin(DCRAD1_in)*sin(PHIR_in-RARAD2_in)/(cos(DCRAD1_in)*sin(RARAD1_in-RARAD2_in))+&
                         sin(DCRAD2_in)*sin(RARAD1_in-PHIR_in)/(cos(DCRAD2_in)*sin(RARAD1_in-RARAD2_in)))
        PHI_in(i) = PHIR_in/DTR         !lon along wp [deg]
        THETA_in(i) = THETAR_in/DTR     !lat along wp [deg]
     ENDDO 
  !WP_FL_DIR=1,3
  ELSE IF(DELRA3_in.gt.PI)THEN
     DO i=1, DIV_WAY_in+1
        PHIR_in = RARAD1_in + (i-1)*((RARAD2_in+2.0_dp*PI) - RARAD1_in)/DIV_WAY_in     !lon [rad]
        THETAR_in = atan(sin(DCRAD1_in)*sin(PHIR_in-RARAD2_in)/(cos(DCRAD1_in)*sin(RARAD1_in-RARAD2_in))+&
                         sin(DCRAD2_in)*sin(RARAD1_in-PHIR_in)/(cos(DCRAD2_in)*sin(RARAD1_in-RARAD2_in)))
        PHI_in(i) = PHIR_in/DTR         !lon along wp [deg]
        THETA_in(i) = THETAR_in/DTR     !lat along wp [deg]
     ENDDO
  !WP_FL_DIR=4,5 
  ELSE     !(DELRA3_in.eq.0.0)
     DO i=1, DIV_WAY_in+1
        THETAR_in = DCRAD2_in + (i-1)*DELDC3_in/DIV_WAY_in  !lon [rad]
        PHI_in(i) = RARAD1_in/DTR       !lon along wp [deg]
        THETA_in(i) = THETAR_in/DTR     !lat along wp [deg]
     ENDDO 
  ENDIF

  !Waypoint output should be in order
  !WP_FL_DIR = 0,3,5: current array order and output order are same
  SELECT CASE (WP_FL_DIR)  
  CASE(0,3,5)
     DO i=1,DIV_WAY_in+1   
       PHI_OUT_in(i)   = PHI_in(i)
       THETA_OUT_in(i) = THETA_in(i)
     ENDDO
!  CASE(3) !Output is converted from 0 =<lon=< +360 into -180 =<lon=< +180 
!     DO i=1,DIV_WAY+1   
!       PHI_OUT(i) = PHI(i)           !lon along wp [deg]
!       THETA_OUT(i) = THETA(i)       !lat along wp [deg]
!     ENDDO
!     DO i=1,DIV_WAY+1
!       IF(PHI_OUT(i).gt.180.0d0)PHI_OUT(i) = PHI_OUT(i) - 360.0d0
!     ENDDO
   
  !WP_FL_DIR = 1,2,4: current array order and output order are opposite
  CASE(1,2,4)
     DO i=1,DIV_WAY_in+1   
       PHI_OUT_in(i)   = PHI_in(DIV_WAY_in+2-i)
       THETA_OUT_in(i) = THETA_in(DIV_WAY_in+2-i)
     ENDDO
!  CASE(1) !Output is converted from 0 =<lon=< +360 into -180 =<lon=< +180
!     DO i=1,DIV_WAY+1   
!       PHI_OUT(i) = PHI(DIV_WAY+2-i)      !lon along wp [deg]
!       THETA_OUT(i) = THETA(DIV_WAY+2-i)  !lat along wp [deg]
!     ENDDO
!     DO i=1,DIV_WAY+1
!       IF(PHI_OUT(i).gt.180.0d0)PHI_OUT(i) = PHI_OUT(i) - 360.0d0
!     ENDDO
  CASE DEFAULT 
      write(*,*) 'Output error: Unknown Waypoints on GC!'
  END SELECT


! After above calculation, the range of variables are:see opt20130215-4,20140512-4
!    -180 =< PHI_OUT_in =< 180 [FL_DIR=0,2,4,5]
!       0 =< PHI_OUT_in =< 360 [FL_DIR=1,3]
! NOTE!! 
!    -180 =< W1_lon_in,W2_lon_in =< 180 [FL_DIR=ALL]
!
 !Current(start) waypoint1  :(PHI_OUT_in(1)(=W1_lon_in), THETA_OUT_in(1)(=W1_lat_in))
 !Nearest point of Waypoint1:(PHI_OUT_in(2)            , THETA_OUT_in(2)            )
 !see20140508-7
 !HY20140509-1 near_lon_in=PHI_OUT_in(2)       !LON
 !HY20140509-1 near_lat_in=THETA_OUT_in(2)     !LAT

  !Calculate tail- or head-wind component along flight direction in (LON,LAT).
  !Theta_wind is angle between flight direction and Ux(wind component).
  !wind_value_in(1):u [m/s] ,20140509-1
  !             (2):v [m/s]
  !             (3):w [m/s]
  
!  IF(W1_lon_in.ne.PHI_OUT_in(1))write(*,*)"W1_lon_in=/PHI_OUT_in(1)",W1_lon_in,PHI_OUT_in(1)
!  IF(W1_lat_in.ne.THETA_OUT_in(1))write(*,*)"W1_lat_in=/THETA_OUT_in(1)",W1_lat_in,THETA_OUT_in(1) 
  
  !CASE:A  20140509-1,20140508-7,20140510-1,2, 20140512-3 
! IF(W1_lon_in.lt.PHI_OUT_in(2) .and. W1_lat_in.ne.THETA_OUT_in(2))THEN   20140512-5
  IF(PHI_OUT_in(1).lt.PHI_OUT_in(2) .and. THETA_OUT_in(1).ne.THETA_OUT_in(2))THEN
     theta_wind        = ATAN2(ABS(THETA_OUT_in(2)-THETA_OUT_in(1)), ABS(PHI_OUT_in(2)-PHI_OUT_in(1)))  !theta_wind is [rad]
!    theta_wind        = ATAN2(ABS(THETA_OUT_in(2)-W1_lat_in), ABS(PHI_OUT_in(2)-W1_lon_in))  !theta_wind is [rad]  20140512-5
!    theta_wind2       = atan2(ABS((THETA_OUT_in(2)-W1_lat_in)*DTR), ABS((PHI_OUT_in(2)-W1_lon_in)*DTR))  !20140512-2
     delta_wind        = ASIN(ABS(W2_alt_in-W1_alt_in)/(DIS_in*1000.0_dp)) !delta_wind[rad], DIS_in*1000[m] 20140514-1
                                                                         !W1_alt_in[m], W2_alt_in[m] 
     ux_wind_component = ABS(wind_value_in(1))*COS(theta_wind)           
     ux_wind_component = ux_wind_component*COS(delta_wind)               ![m/s] 20140514-1,20140523-1 
     vy_wind_component = ABS(wind_value_in(2))*SIN(theta_wind)
     vy_wind_component = vy_wind_component*COS(delta_wind)               ![m/s] 20140514-1
     wz_wind_component = ABS(wind_value_in(3))*SIN(delta_wind)           ![m/s] 20140519-1

     !Consider wind direct to flight direct.
     !Add both wind components and clulcate new AC speed along flight dir in (LON,LAT) plane.  
     !CASE:A-1
     IF(THETA_OUT_in(1).lt.THETA_OUT_in(2))THEN
!    IF(W1_lat_in.lt.THETA_OUT_in(2))THEN  20140512-5
        !Add wind component of Ux, which is along to flight dir, to AC speed. 
        IF(wind_value_in(1).ge.0.0_dp)AC_V_temp   = AC_V_in   + ux_wind_component   
        IF(wind_value_in(1).lt.0.0_dp)AC_V_temp   = AC_V_in   - ux_wind_component   
        !Add wind component of Vy, which is along to flight dir, to AC speed.
        IF(wind_value_in(2).ge.0.0_dp)AC_V_temp = AC_V_temp + vy_wind_component  !20140519-1 
        IF(wind_value_in(2).lt.0.0_dp)AC_V_temp = AC_V_temp - vy_wind_component  
        !Add wind component of Wz, which is along to flight dir, to AC speed.
        IF(W1_alt_in.le.W2_alt_in)THEN
           IF(wind_value_in(3).ge.0.0_dp)new_AC_V_in = AC_V_temp + wz_wind_component  !20140519-1 
           IF(wind_value_in(3).lt.0.0_dp)new_AC_V_in = AC_V_temp - wz_wind_component   
        ELSEIF(W1_alt_in.gt.W2_alt_in)THEN
           IF(wind_value_in(3).ge.0.0_dp)new_AC_V_in = AC_V_temp - wz_wind_component  !20140519-1 
           IF(wind_value_in(3).lt.0.0_dp)new_AC_V_in = AC_V_temp + wz_wind_component   
        ENDIF
     
     !CASE:A-2
     ELSEIF(THETA_OUT_in(1).gt.THETA_OUT_in(2))THEN
!    ELSEIF(W1_lat_in.gt.THETA_OUT_in(2))THEN  20140512-5
        !Add wind component of Ux, which is along to flight dir, to AC speed.
        IF(wind_value_in(1).ge.0.0_dp)AC_V_temp   = AC_V_in   + ux_wind_component   
        IF(wind_value_in(1).lt.0.0_dp)AC_V_temp   = AC_V_in   - ux_wind_component   
        !Add wind component of Vy, which is along to flight dir, to AC speed.
        IF(wind_value_in(2).ge.0.0_dp)AC_V_temp = AC_V_temp - vy_wind_component   !20140519-2
        IF(wind_value_in(2).lt.0.0_dp)AC_V_temp = AC_V_temp + vy_wind_component   
        !Add wind component of Wz, which is along to flight dir, to AC speed.
        IF(W1_alt_in.le.W2_alt_in)THEN
           IF(wind_value_in(3).ge.0.0_dp)new_AC_V_in = AC_V_temp + wz_wind_component  !20140519-2 
           IF(wind_value_in(3).lt.0.0_dp)new_AC_V_in = AC_V_temp - wz_wind_component   
        ELSEIF(W1_alt_in.gt.W2_alt_in)THEN
           IF(wind_value_in(3).ge.0.0_dp)new_AC_V_in = AC_V_temp - wz_wind_component  !20140519-2 
           IF(wind_value_in(3).lt.0.0_dp)new_AC_V_in = AC_V_temp + wz_wind_component   
        ENDIF
     ENDIF
  
  !CASE:B 20140508-7, 20140510-3, 20140512-3
  ELSEIF(PHI_OUT_in(1).gt.PHI_OUT_in(2) .and. THETA_OUT_in(1).ne.THETA_OUT_in(2))THEN
! ELSEIF(W1_lon_in.gt.PHI_OUT_in(2) .and. W1_lat_in.ne.THETA_OUT_in(2))THEN  20140512-5
     theta_wind        = ATAN2(ABS(THETA_OUT_in(2)-THETA_OUT_in(1)), ABS(PHI_OUT_in(1)-PHI_OUT_in(2)))   !theta_wind is [rad]
!    theta_wind        = ATAN2(ABS(THETA_OUT_in(2)-W1_lat_in), ABS(W1_lon_in-PHI_OUT_in(2)))   !theta_wind is [rad]  20140512-5
!    theta_wind2       = atan2(ABS((THETA_OUT_in(2)-W1_lat_in)*DTR), ABS((W1_lon_in-PHI_OUT_in(2))*DTR))  !20140512-2
     delta_wind        = ASIN(ABS(W2_alt_in-W1_alt_in)/(DIS_in*1000.0_dp)) !delta_wind[rad], DIS_in*1000[m] 20140516-1
                                                                         !W1_alt_in[m], W2_alt_in[m] 
     ux_wind_component = ABS(wind_value_in(1))*COS(theta_wind)           ![m/s]
     ux_wind_component = ux_wind_component*COS(delta_wind)               ![m/s] 20140516-1,20140523-1 
     vy_wind_component = ABS(wind_value_in(2))*SIN(theta_wind)
     vy_wind_component = vy_wind_component*COS(delta_wind)               ![m/s] 20140516-1
     wz_wind_component = ABS(wind_value_in(3))*SIN(delta_wind)           ![m/s] 20140519-1
     
     !CASE:B-3
     IF(THETA_OUT_in(1).lt.THETA_OUT_in(2))THEN
!    IF(W1_lat_in.lt.THETA_OUT_in(2))THEN  20140512-5
        !Add wind component of Ux, which is along to flight dir, to AC speed.
        IF(wind_value_in(1).ge.0.0_dp)AC_V_temp   = AC_V_in   - ux_wind_component   
        IF(wind_value_in(1).lt.0.0_dp)AC_V_temp   = AC_V_in   + ux_wind_component   
        !Add wind component of Vy, which is along to flight dir, to AC speed.
        IF(wind_value_in(2).ge.0.0_dp)AC_V_temp = AC_V_temp + vy_wind_component  !20140519-3 
        IF(wind_value_in(2).lt.0.0_dp)AC_V_temp = AC_V_temp - vy_wind_component   
        !Add wind component of Wz, which is along to flight dir, to AC speed.
        IF(W1_alt_in.le.W2_alt_in)THEN
           IF(wind_value_in(3).ge.0.0_dp)new_AC_V_in = AC_V_temp + wz_wind_component  !20140519-3 
           IF(wind_value_in(3).lt.0.0_dp)new_AC_V_in = AC_V_temp - wz_wind_component   
        ELSEIF(W1_alt_in.gt.W2_alt_in)THEN
           IF(wind_value_in(3).ge.0.0_dp)new_AC_V_in = AC_V_temp - wz_wind_component  !20140519-3 
           IF(wind_value_in(3).lt.0.0_dp)new_AC_V_in = AC_V_temp + wz_wind_component   
        ENDIF
     
     !CASE:B-4
     ELSEIF(THETA_OUT_in(1).gt.THETA_OUT_in(2))THEN  
!    ELSEIF(W1_lat_in.gt.THETA_OUT_in(2))THEN  20140512-5
        !Add wind component of Ux, which is along to flight dir, to AC speed.
        IF(wind_value_in(1).ge.0.0_dp)AC_V_temp   = AC_V_in   - ux_wind_component   
        IF(wind_value_in(1).lt.0.0_dp)AC_V_temp   = AC_V_in   + ux_wind_component   
        !Add wind component of Vy, which is along to flight dir, to AC speed.
        IF(wind_value_in(2).ge.0.0_dp)AC_V_temp = AC_V_temp - vy_wind_component   !20140520-1 
        IF(wind_value_in(2).lt.0.0_dp)AC_V_temp = AC_V_temp + vy_wind_component   
        !Add wind component of Wz, which is along to flight dir, to AC speed.
        IF(W1_alt_in.le.W2_alt_in)THEN
           IF(wind_value_in(3).ge.0.0_dp)new_AC_V_in = AC_V_temp + wz_wind_component  !20140520-1 
           IF(wind_value_in(3).lt.0.0_dp)new_AC_V_in = AC_V_temp - wz_wind_component   
        ELSEIF(W1_alt_in.gt.W2_alt_in)THEN
           IF(wind_value_in(3).ge.0.0_dp)new_AC_V_in = AC_V_temp - wz_wind_component  !20140520-1 
           IF(wind_value_in(3).lt.0.0_dp)new_AC_V_in = AC_V_temp + wz_wind_component   
        ENDIF
     ENDIF
  ENDIF

  !CASE:C  use Vy(=wind_value_in(2)) directly. 20140508-7, 20140510-5, 20140512-3
  IF(PHI_OUT_in(1).eq.PHI_OUT_in(2) .and. THETA_OUT_in(1).lt.THETA_OUT_in(2))THEN
! IF(W1_lon_in.eq.PHI_OUT_in(2) .and. W1_lat_in.lt.THETA_OUT_in(2))THEN  20140512-5
     delta_wind        = ASIN(ABS(W2_alt_in-W1_alt_in)/(DIS_in*1000.0_dp))    !delta_wind[rad], DIS_in*1000[m] 20140516-2
                                                                              !W1_alt_in[m], W2_alt_in[m] 
     vy_wind_component = ABS(wind_value_in(2))                                ![m/s] 20140523-2
     vy_wind_component = vy_wind_component*COS(delta_wind)                    ![m/s] 20140516-2,20140514-3
     wz_wind_component = ABS(wind_value_in(3))*SIN(delta_wind)                ![m/s] 20140520-2
     
     !Add wind component of Vy, which is along to flight dir, to AC speed.
     IF(wind_value_in(2).ge.0.0_dp)AC_V_temp = AC_V_in + vy_wind_component    ![m/s]20140520-2 
     IF(wind_value_in(2).lt.0.0_dp)AC_V_temp = AC_V_in - vy_wind_component    !   
     !Add wind component of Wz, which is along to flight dir, to AC speed.
     IF(W1_alt_in.le.W2_alt_in)THEN
        IF(wind_value_in(3).ge.0.0_dp)new_AC_V_in = AC_V_temp + wz_wind_component  !20140520-2 
        IF(wind_value_in(3).lt.0.0_dp)new_AC_V_in = AC_V_temp - wz_wind_component   
     ELSEIF(W1_alt_in.gt.W2_alt_in)THEN
        IF(wind_value_in(3).ge.0.0_dp)new_AC_V_in = AC_V_temp - wz_wind_component  !20140520-2 
        IF(wind_value_in(3).lt.0.0_dp)new_AC_V_in = AC_V_temp + wz_wind_component   
     ENDIF

  !CASE:D  use Vy(=wind_value_in(2)) directly. 20140508-7, 20140510-5, 20140512-3
  ELSEIF(PHI_OUT_in(1).eq.PHI_OUT_in(2) .and. THETA_OUT_in(1).gt.THETA_OUT_in(2))THEN
! ELSEIF(W1_lon_in.eq.PHI_OUT_in(2) .and. W1_lat_in.gt.THETA_OUT_in(2))THEN  20140512-5
     delta_wind        = ASIN(ABS(W2_alt_in-W1_alt_in)/(DIS_in*1000.0_dp))    !delta_wind[rad], DIS_in*1000[m] 20140516-2
                                                                              !W1_alt_in[m], W2_alt_in[m] 
     vy_wind_component = ABS(wind_value_in(2))                                ![m/s] 20140523-2
     vy_wind_component = vy_wind_component*COS(delta_wind)                    ![m/s] 20140516-2,20140514-4
     wz_wind_component = ABS(wind_value_in(3))*SIN(delta_wind)                ![m/s] 20140520-3
     
     !Add wind component of Vy, which is along to flight dir, to AC speed.
     IF(wind_value_in(2).ge.0.0_dp)AC_V_temp = AC_V_in - vy_wind_component  ![m/s] 20140520-3 
     IF(wind_value_in(2).lt.0.0_dp)AC_V_temp = AC_V_in + vy_wind_component  !
     !Add wind component of Wz, which is along to flight dir, to AC speed.
     IF(W1_alt_in.le.W2_alt_in)THEN
        IF(wind_value_in(3).ge.0.0_dp)new_AC_V_in = AC_V_temp + wz_wind_component  !20140520-3 
        IF(wind_value_in(3).lt.0.0_dp)new_AC_V_in = AC_V_temp - wz_wind_component   
     ELSEIF(W1_alt_in.gt.W2_alt_in)THEN
        IF(wind_value_in(3).ge.0.0_dp)new_AC_V_in = AC_V_temp - wz_wind_component  !20140520-3 
        IF(wind_value_in(3).lt.0.0_dp)new_AC_V_in = AC_V_temp + wz_wind_component   
     ENDIF
  ENDIF

  !CASE:E  use Ux(=wind_value_in(1)) directly. 20140508-7, 20140510-5, 20140512-3
  IF(THETA_OUT_in(1).eq.THETA_OUT_in(2) .and. PHI_OUT_in(1).lt.PHI_OUT_in(2))THEN
! IF(W1_lat_in.eq.THETA_OUT_in(2) .and. W1_lon_in.lt.PHI_OUT_in(2))THEN  20140512-5
     delta_wind        = ASIN(ABS(W2_alt_in-W1_alt_in)/(DIS_in*1000.0_dp))    !delta_wind[rad], DIS_in*1000[m] 20140516-2
                                                                              !W1_alt_in[m], W2_alt_in[m] 
     ux_wind_component = ABS(wind_value_in(1))                                ![m/s] 20140523-2
     ux_wind_component = ux_wind_component*COS(delta_wind)                    ![m/s] 20140516-3,20140514-5 
     wz_wind_component = ABS(wind_value_in(3))*SIN(delta_wind)                ![m/s] 20140520-4
     
     !Add wind component of Ux, which is along to flight dir, to AC speed.
     IF(wind_value_in(1).ge.0.0_dp)AC_V_temp = AC_V_in + ux_wind_component  ![m/s] 20140520-4 
     IF(wind_value_in(1).lt.0.0_dp)AC_V_temp = AC_V_in - ux_wind_component  !
     !Add wind component of Wz, which is along to flight dir, to AC speed.
     IF(W1_alt_in.le.W2_alt_in)THEN
        IF(wind_value_in(3).ge.0.0_dp)new_AC_V_in = AC_V_temp + wz_wind_component  !20140520-4 
        IF(wind_value_in(3).lt.0.0_dp)new_AC_V_in = AC_V_temp - wz_wind_component   
     ELSEIF(W1_alt_in.gt.W2_alt_in)THEN
        IF(wind_value_in(3).ge.0.0_dp)new_AC_V_in = AC_V_temp - wz_wind_component  !20140520-4 
        IF(wind_value_in(3).lt.0.0_dp)new_AC_V_in = AC_V_temp + wz_wind_component   
     ENDIF

  !CASE:F  use Ux(=wind_value_in(1)) directly. 20140508-7, 20140510-5, 20140512-3
  ELSEIF(THETA_OUT_in(1).eq.THETA_OUT_in(2) .and. PHI_OUT_in(1).gt.PHI_OUT_in(2))THEN
! ELSEIF(W1_lat_in.eq.THETA_OUT_in(2) .and. W1_lon_in.gt.PHI_OUT_in(2))THEN  20140512-5
     delta_wind        = ASIN(ABS(W2_alt_in-W1_alt_in)/(DIS_in*1000.0_dp))    !delta_wind[rad], DIS_in*1000[m] 20140516-2
                                                                              !W1_alt_in[m], W2_alt_in[m] 
     ux_wind_component = ABS(wind_value_in(1))                                ![m/s] 20140523-2
     ux_wind_component = ux_wind_component*COS(delta_wind)                    ![m/s] 20140516-3,20140514-6 
     wz_wind_component = ABS(wind_value_in(3))*SIN(delta_wind)                ![m/s] 20140520-5
     
     !Add wind component of Ux, which is along to flight dir, to AC speed.
     IF(wind_value_in(1).ge.0.0_dp)AC_V_temp = AC_V_in - ux_wind_component  ![m/s] 20140520-5  
     IF(wind_value_in(1).lt.0.0_dp)AC_V_temp = AC_V_in + ux_wind_component  !   
     !Add wind component of Wz, which is along to flight dir, to AC speed.
     IF(W1_alt_in.le.W2_alt_in)THEN
        IF(wind_value_in(3).ge.0.0_dp)new_AC_V_in = AC_V_temp + wz_wind_component  !20140520-5 
        IF(wind_value_in(3).lt.0.0_dp)new_AC_V_in = AC_V_temp - wz_wind_component   
     ELSEIF(W1_alt_in.gt.W2_alt_in)THEN
        IF(wind_value_in(3).ge.0.0_dp)new_AC_V_in = AC_V_temp - wz_wind_component  !20140520-5 
        IF(wind_value_in(3).lt.0.0_dp)new_AC_V_in = AC_V_temp + wz_wind_component   
     ENDIF
  ENDIF
!  write(*,*)"(Wind)wind_components=",ux_wind_component,"[m/s]",vy_wind_component,"[m/s]",wz_wind_component,"[m/s]"  

  END SUBROUTINE wind_effect 
!----------------------------------------------------------------------

END MODULE messy_airtraf_wind

