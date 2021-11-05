! **********************************************************************
! 
! SUB_SUBMODEL ROUTINES FOR MESSy SUBMODEL airtraf
!
! THIS SUB_SUBMODEL IS USED TO CALCULATE GREAT CIRCLE(GC)
!
! Author : Hiroshi Yamashita, DLR-IPA, March 2012
!          Volker Grewe, DLR-IPA, March 2012
!
! References:
!
! * None yet
!
! see SMCL for more details
! **********************************************************************
! Objective: Simulation of great circle on arbitrary city pairs
! Input:     Flighplan information
! Output:    Flight Time table including waypoints
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
! 
!=========================================================================
!
! **********************************************************************
MODULE messy_airtraf_gc
! **********************************************************************

!USE ONLY
 USE messy_main_constants_mem, ONLY: DP, PI, DTR, OneDay, radius_earth, g

 USE messy_main_tools,         ONLY: nn_index !20131214-1,20140530-1

! USE messy_airtraf_e5,  ONLY:D_lon, D_lat, D_time, A_lon, A_lat 

! USE messy_airtraf,  ONLY: NOx, H2O, DIST, FUEL
! USE messy_airtraf,  ONLY: nwaypoints


  IMPLICIT NONE
  ! GLOBAL PARAMETERS
  !Comments:If we can use the information on aircraft type for each route, 
  !         average flight speed of each aircraft type can be utilized as AC_MACH. 
  REAL(DP), PARAMETER, PUBLIC :: ALT     = 10668._dp !12496.8_dp  !8839.2_dp !Flight altitude(=geopotential altitude,20131105-1), [m] 
  REAL(DP), PARAMETER, PUBLIC :: AC_MACH = 0.82_dp          !Flight Mach number(=Mcr from A333_ptf), [-] see20130719-1
                                                                       !AC_MACH should be from Aircraft table(BADA A333_ptf) 
! DOUBLE PRECISION, PARAMETER, PUBLIC :: SOUND_V = 299.463_dp          !Speed of sound from BADA_ISA,[m/s] see20130719-1 
 !REAL(DP), PUBLIC            :: SOUND_V                       !Speed of sound,[m/s] see20131112-1,20140601-2 
 !REAL(DP), PUBLIC            :: BADA_FUEL_NOM                 !Fuel from BADA A333_ptf,[kg/s] see20131112-1,20140601-3 
 !REAL(DP), PUBLIC           :: total_flight_time_s           ![s],20131003-3,20131114-2,20140601-2
 
  ! The parameter should be used from SMIL
! INTEGER, PARAMETER, PUBLIC   :: DIV_WAY =60 !DIV_WAY = nwaypoints-1,this variable should be determined automatically,see20130712-1.
  INTEGER, PRIVATE             :: DIV_WAY     !20140721-2,20140722-1,3
  INTEGER, PARAMETER, PRIVATE  :: D_lon=1, D_lat=2, D_time=3, A_lon=4, A_lat=5
  REAL(DP), PARAMETER, PRIVATE :: adiabatic_index_air = 1.4_dp         ![-]20140724-2 
  REAL(DP), PARAMETER, PRIVATE :: gas_constant_air    = 287.05287_dp   ![m^2/(k*s^2)]20140724-2 

  !Parameters for ac_routes properties
  INTEGER, PARAMETER, PUBLIC                  :: props_lon       = 1
  INTEGER, PARAMETER, PUBLIC                  :: props_lat       = 2
  INTEGER, PARAMETER, PUBLIC                  :: props_alt       = 3
  INTEGER, PARAMETER, PUBLIC                  :: props_time      = 4
  INTEGER, PARAMETER, PUBLIC                  :: props_ac_speed  = 5
  INTEGER, PARAMETER, PUBLIC                  :: props_dist      = 6   !subroutine calculate_emissions_along_trajectory
  INTEGER, PARAMETER, PUBLIC                  :: props_fuel_use  = 7   !subroutine calculate_emissions_along_trajectory
  INTEGER, PARAMETER, PUBLIC                  :: props_emis_nox  = 8   !subroutine calculate_emissions_along_trajectory
  INTEGER, PARAMETER, PUBLIC                  :: props_emis_h2o  = 9   !subroutine calculate_emissions_along_trajectory
  !Yin_20170423
  INTEGER, PARAMETER, PUBLIC                  :: props_potcov    = 10  !potcov along trajectory
  INTEGER, PARAMETER, PUBLIC                  :: props_atr20o3   = 11  !atr20o3 along trajectory
  INTEGER, PARAMETER, PUBLIC                  :: props_atr20ch4  = 12  !atr20ch4 along trajectory
  INTEGER, PARAMETER, PUBLIC                  :: props_atr20h2o  = 13  !atr20h2o along trajectory
  INTEGER, PARAMETER, PUBLIC                  :: props_atr20cpc  = 14  !atr20cpc along trajectory
  INTEGER, PARAMETER, PUBLIC                  :: props_atr20co2  = 15  !atr20cpc along trajectory
  INTEGER, PARAMETER, PUBLIC                  :: props_atr20tot  = 16  !atr20tot along trajectory
  !Yin_20170423

  ! PRIVATE PARAMETERS
  INTEGER, PARAMETER, PRIVATE :: ISWITCH = 2          !0: Spherical law of cosines
                                                      !1: Haversine formula
                                                      !2: Vincenty formula
  INTEGER, PARAMETER, PRIVATE :: OUT_SWITCH = 0       !0: -180 =< lon =< +180 (Default) 
                                                      !1:    0 =< lon =< +360 (only cases FL_DIR = 1 and 3)

  INTEGER, PARAMETER, PRIVATE :: VTAS_SWITCH = 0      !0: Calculate vtas by using sound speed from ECHAM5 data (Default)20140528-1 
                                                      !1: Calculate vtas by using sound speed from interpolating BADA(ISA atm.) data.  
  
  INTEGER, PARAMETER, PRIVATE :: WIND_SWITCH = 1      !0: No wind 20140531-1 
                                                      !1: Include wind effects(u,v,w)(Default)  
  
  ! PRIVATE
  !REAL(DP), DIMENSION(:), PRIVATE :: THETA_OUT(DIV_WAY+1), PHI_OUT(DIV_WAY+1)   !20140528-2
! REAL(DP), DIMENSION(:), PRIVATE  :: XXN(DIV_WAY+1),YYN(DIV_WAY+1),ZZN(DIV_WAY+1)  !20140721-3
 !REAL(DP), PRIVATE   :: AC_V   !20140601-2
  REAL(DP), PRIVATE   :: RARAD1,RARAD2,DCRAD1,DCRAD2 
  REAL(DP), PRIVATE   :: DELDC2,DELRA,DELRA2,SINDIS,SINDIS1,SINDIS2
  REAL(DP), PRIVATE   :: THETAR, PHIR
  REAL(DP), PRIVATE   :: gc_D_time
  INTEGER, PRIVATE    :: FL_DIR        !see 20111219-1,20130514-2-left

  ! PUBLIC SUBROUTINES (to be called from messy_airtraf.f90)
  PUBLIC :: gc_calculate

  ! PRIVATE SUBROUTINES
 !PRIVATE :: gc_distance 
  PRIVATE :: gc_waypoints 
  PRIVATE :: gc_ac_traj 
  PRIVATE :: wind_effect   !20140531-2 

  CONTAINS
!----------------------------------------------------------------------
  SUBROUTINE gc_calculate(gc_nwaypoints, gc_ac_routes, gc_p2_ac_routes_desc, gc_philon, gc_philat,    &
                          gc_zgl_geopot_3d, gc_uwind_g, gc_vwind_g, gc_v_z_g,                         &
                          gc_t_scb_g, gc_cpc_g, gc_fl_time, gc_bada_fuel_nom, gc_vtas_along_traj,     &
                          gc_fl_direction, gc_fl_distance, gc_cpc)  !20140527-4,20140721-2,20141215-1,20141216-1
!----------------------------------------------------------------------
!
! This subroutine contains the subroutines on GC and controls I/O.
! This SUBROUTINE is called by messy_airtraf(SMCL).
! Yin_20170423  This subroutine was modified to add contrail potential coverage
!
  IMPLICIT NONE
  INTEGER, INTENT(IN)                              :: gc_nwaypoints      !20140721-2 
  REAL(DP), DIMENSION(:,:), INTENT(OUT)            :: gc_ac_routes
  REAL(DP), DIMENSION(:),   INTENT(IN)             :: gc_p2_ac_routes_desc
  REAL(DP), DIMENSION(:),   INTENT(IN)             :: gc_philon          !added 20140527-4
  REAL(DP), DIMENSION(:),   INTENT(IN)             :: gc_philat
  REAL(DP), DIMENSION(:,:,:), INTENT(IN)           :: gc_zgl_geopot_3d   !global field  [m^2/s^2] 
  REAL(DP), DIMENSION(:,:,:), INTENT(IN)           :: gc_uwind_g         !global field  [m/s]
  REAL(DP), DIMENSION(:,:,:), INTENT(IN)           :: gc_vwind_g         !global field  [m/s]
  REAL(DP), DIMENSION(:,:,:), INTENT(IN)           :: gc_v_z_g           !global field  [m/s]
  REAL(DP), DIMENSION(:,:,:), INTENT(IN)           :: gc_t_scb_g         !global field  [K] 20140522-3
  !Yin_20170423
  REAL(DP), DIMENSION(:,:,:), INTENT(IN)           :: gc_cpc_g           !global field  [-] 20161212
  !Yin_20170423
! REAL(DP), DIMENSION(:),   INTENT(IN), OPTIONAL   :: gc_philat  !20131213-1,20131214-1
  REAL(DP), INTENT(OUT)                            :: gc_fl_time         !20140527-4
  REAL(DP), INTENT(OUT)                            :: gc_bada_fuel_nom   !20140527-4,20140601-2  [kg/s]
  !Yin_20170423
  REAL(DP), INTENT(OUT)                            :: gc_cpc             !20161212
  !Yin_20170423
  REAL(DP), DIMENSION(:),   ALLOCATABLE            :: gc_zalt            ! geopotential altitude[m],20131105-1,20140527-4
  REAL(DP)                                         :: gc_wind_along_traj(gc_nwaypoints,3)  ![m/s] 20140429-2 
  REAL(DP),INTENT(INOUT)                           :: gc_vtas_along_traj(gc_nwaypoints)    ![m/s] 20140527-4
  INTEGER, INTENT(OUT)                             :: gc_fl_direction    !20141215-1 
  INTEGER, INTENT(OUT)                             :: gc_fl_distance     !20141216-1

  REAL(DP) :: XXN(gc_nwaypoints),YYN(gc_nwaypoints),ZZN(gc_nwaypoints)  !20140721-3
  !Yin_20170423
  REAL(DP) :: gc_cpc_along_traj(gc_nwaypoints)    ![-] 20161212
  !Yin_20170423
  REAL(DP) :: BADA_DATA(8,5) !this should be given by Aircraft table.20130803-1,20131108-1
  REAL(DP) :: FLIGHT_TIME(gc_nwaypoints), FLIGHT_TIME_JULIAN(gc_nwaypoints)
  REAL(DP) :: FLIGHT_DISTANCE(gc_nwaypoints), FLIGHT_SPEED(gc_nwaypoints)   !,FUEL_USE(DIV_WAY+1)
  REAL(DP) :: contrail_coverage(gc_nwaypoints)   ![km] 20170602-1
  REAL(DP) :: T_FLIGHT_TIME, A_time, T_FLIGHT_DISTANCE !A_time in Julian date
  REAL(DP) :: gc_D_lat        !latitude of departure city 
  REAL(DP) :: gc_D_lon        !longitude of departure city
  REAL(DP) :: gc_A_lat        !latitude of arrivar city 
  REAL(DP) :: gc_A_lon        !longitude of arrival city
  REAL(DP) :: SOUND_V         !speed of sound [m/s] in standard armosphere
  REAL(DP) :: temp_XXN        !20140605-1
  INTEGER  :: i, j, gc_i_lat, gc_j_lon, gc_i_alt   !20131213-1,20131214-1,20140530-1
  INTEGER  :: nlev_gc_zalt

  !BADA DATA input,A333_ptf,BADA_ISA_TABLE, 20131108-1
  DATA ((BADA_DATA(i,j),i=1,8),j=1,3) /280.0_dp,290.0_dp,310.0_dp,330.0_dp,350.0_dp,370.0_dp,390.0_dp,410.0_dp,         &  !FL[ft]
                                       305.79_dp,304.48_dp,301.86_dp,299.21_dp,296.54_dp,295.07_dp,295.07_dp,295.07_dp, &  !a[m/s]
                                       112.4_dp,108.7_dp,101.7_dp,95.5_dp,90.0_dp,85.5_dp,81.9_dp,79.0_dp/      !Fuel[nom, kg/min]
 
  !Determine:DIV_WAY,20140721-2,20140722-1
  DIV_WAY = gc_nwaypoints-1

!  write(*,*)"gc_nwaypoints_check",gc_nwaypoints 
  !Allocation for gc_zalt, 20140527-4 
  nlev_gc_zalt = SIZE(gc_zgl_geopot_3d,dim=2)   !ALT: in 2nd dimention.
  ALLOCATE(gc_zalt(nlev_gc_zalt)) 
  !write(*,*)"nlev_gc_zalt=",nlev_gc_zalt   !20140620-8
 
  !Initialization 20140527-4
  gc_fl_time        =0.0_dp
  gc_bada_fuel_nom  =0.0_dp  !20140601-2
  gc_wind_along_traj=0.0_dp
  gc_vtas_along_traj=0.0_dp
  !Yin_20170423
  gc_cpc_along_traj =0.0_dp
  !Yin_20170423
  gc_zalt           =0.0_dp
  SOUND_V           =0.0_dp
 !BADA_FUEL_NOM     =0.0_dp  !20140602-3 
  gc_fl_direction   =0       !20141215-1
  gc_fl_distance    =0       !20141216-1

  do i=1,8 
     BADA_DATA(i,4)=BADA_DATA(i,1)*100.0_dp*0.30480_dp   ![ft] to [m],20131108-1,20131031-1
     BADA_DATA(i,5)=BADA_DATA(i,3)/60.0_dp               ![kg/min] to [kg/s]
     !write(*,*)'BADA_DATA_check',(BADA_DATA(i,j),j=1,5) !20140620-8 
  enddo
  
  IF(VTAS_SWITCH.eq.1)THEN  !20140528-1
  !Calculate SOUND_V at ALT by using BADA_DATA array.
  !SOUND_V,[m/s]
     do i=1,7
        if((BADA_DATA(i,4)<=ALT).and.(ALT<BADA_DATA(i+1,4)))then
           SOUND_V=(BADA_DATA(i,2)*(ALT-BADA_DATA(i+1,4))-BADA_DATA(i+1,2)*(ALT-BADA_DATA(i,4)))/ &
                   (BADA_DATA(i,4)-BADA_DATA(i+1,4))
!           write(*,*)"if_soundv_1",i,ALT,BADA_DATA(i+1,4)
           exit
        elseif(ALT.lt.BADA_DATA(1,4))then
           SOUND_V=BADA_DATA(1,2)
!           write(*,*)"if_soundv_2",i,ALT,BADA_DATA(1,4)
           exit
        elseif(ALT.ge.BADA_DATA(8,4))then
           SOUND_V=BADA_DATA(8,2)
!           write(*,*)"if_soundv_3",i,ALT,BADA_DATA(8,4)
           exit  
        endif
     enddo
!     write(*,*)'SOUND_V_check',SOUND_V
  ENDIF 
  !Calculate gc_bada_fuel_nom at ALT by using BADA_DATA array.20140528-1
  !gc_bada_fuel_nom,[kg/s]
  do i=1,7
     if((BADA_DATA(i,4)<=ALT).and.(ALT<BADA_DATA(i+1,4)))then
        gc_bada_fuel_nom=(BADA_DATA(i,5)*(ALT-BADA_DATA(i+1,4))-BADA_DATA(i+1,5)*(ALT-BADA_DATA(i,4)))/ &
                         (BADA_DATA(i,4)-BADA_DATA(i+1,4))
!        write(*,*)"if_bada_fuel_1",i,ALT,BADA_DATA(i+1,4)
        exit
     elseif(ALT.lt.BADA_DATA(1,4))then
        gc_bada_fuel_nom=BADA_DATA(1,5)
!        write(*,*)"if_bada_fuel_2",i,ALT,BADA_DATA(1,4)
        exit
     elseif(ALT.ge.BADA_DATA(8,4))then
        gc_bada_fuel_nom=BADA_DATA(8,5)
!        write(*,*)"if_bada_fuel_3",i,ALT,BADA_DATA(8,4)
        exit  
     endif
  enddo
!  write(*,*)'gc_bada_fuel_nom_check',gc_bada_fuel_nom

  !Flightplan information from SMCL 
  gc_D_lon  = gc_p2_ac_routes_desc(D_lon)        
  gc_D_lat  = gc_p2_ac_routes_desc(D_lat)         
  gc_D_time = gc_p2_ac_routes_desc(D_time)    !Julian 
  gc_A_lon  = gc_p2_ac_routes_desc(A_lon)      
  gc_A_lat  = gc_p2_ac_routes_desc(A_lat)        
!  OPEN(40,file="north_check.dat",position='APPEND')  !20131203-3
!     if(gc_D_lat.lt.0.0_dp)write(40,*)"1",gc_D_lat,gc_D_lon,gc_A_lat,gc_A_lon
!     if(gc_A_lat.lt.0.0_dp)write(40,*)"2",gc_D_lat,gc_D_lon,gc_A_lat,gc_A_lon
!     if(gc_A_lat.lt.0.0_dp.and.gc_D_lat.lt.0.0_dp)write(40,*)"3",gc_D_lat,gc_D_lon,gc_A_lat,gc_A_lon
!  CLOSE(40)
!  write(*,*)'SUBSUB_check',gc_D_lon,gc_D_lat,gc_D_time,gc_A_lon,gc_A_lat
!  write(*,*)'SUBSUB_check',D_lon,D_lat,D_time,A_lon,A_lat
!  write(*,*)'SUBSUB_check',PI, DTR,OneDay,radius_earth

  !Initialization
  T_FLIGHT_TIME = 0.0_dp
!HY  AC_V   = AC_MACH*SOUND_V/1000.0_dp    !Aircraft speed [km/s], 20130719-1,20131031-1,-2,20131104-1,20131114-2,20130315-2 
                                           !20140528-1
  !Definition of the city being placed at higher lon. 
  IF(gc_D_lon.lt.gc_A_lon)THEN
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
  gc_fl_direction = FL_DIR   !20141215-1

  !nn_index check 20131213-1,20131214-1
!  CALL nn_index(gc_philat(:), 50.0_dp, gc_i_lat)
!  write(*,*)"gc_i_lat_SMCL=",gc_i_lat

  !Calculation of GC distance
  !HY 20120810
  !CALL gc_distance(COSIS, T_DIS, T_DISN)

  !Calculation of waypoints on GC
  CALL gc_waypoints(XXN, YYN)   
  !Detremine altitude(constant) of all waypoints.
  ZZN = ALT     ![m] 20140528-2
  ! After this CALL, waypoints are detemined. The range of variables are:see opt20130215-4,20140512-4
  !    -180 =< XXN(LON) =< 180 [FL_DIR=0,2,4,5]
  !       0 =< XXN(LON) =< 360 [FL_DIR=1,3]
!  write(*,*)'CHECK_(LON, LAT, ALT) TRAJ'
!  do i=1, DIV_WAY+1  
!     write(*,*)XXN(i),YYN(i),ZZN(i)
!  enddo

     !see20140530-1
     temp_XXN=0.0_dp     !20140605-1
     DO i=1,DIV_WAY+1    !?? DIV_WAY??
        !Extract wind values by nn_index for each wp. 
        !nn_index check  20140212-2,20130905-3
        !nn_index can find nearest index.
        
        !If [-180=< XXN(i) < 0], add +360 to XXN(i) and use temp_XXN for nn_index search.20140605-1,20140604-3 
        IF((-180.0_dp.le.XXN(i)).and.(XXN(i).lt.0.0_dp))THEN
           temp_XXN = XXN(i)+360.0_dp
           CALL nn_index(gc_philon(:), temp_XXN, gc_j_lon) !lon value on waypoints,20140605-1
        ELSE   
           temp_XXN = XXN(i)
           CALL nn_index(gc_philon(:), temp_XXN, gc_j_lon)   !lon value on waypoints
        ENDIF   
       !CALL nn_index(gc_philon(:), XXN(i), gc_j_lon)  !lon value on waypoints
        CALL nn_index(gc_philat(:), YYN(i), gc_i_lat)  !lat value on waypoints
        gc_zalt(:)=gc_zgl_geopot_3d(gc_j_lon, :, gc_i_lat)/g   ![m] geopotential altitude,20131105-1,20140430-1
!        write(*,*)'(GC)g=',g,'(GC)nlev_gc_zalt=',nlev_gc_zalt
!        write(*,*)'(GC)gc_zgl_geopot_3d=',gc_zgl_geopot_3d(gc_j_lon,:,gc_i_lat)/g
!        write(*,*)'(GC)gc_zalt=',gc_zalt(:)
        CALL nn_index(gc_zalt(:), ZZN(i), gc_i_alt)
!        write(*,*)'(GC)alt=', ZZN(i)
!        write(*,*)"Nearest indices=", gc_i_lat, gc_j_lon, gc_i_alt
!        write(*,*)"nn_index_lat_check", gc_i_lat, 'YYN(i)=',YYN(i)   !20140604-1
!        write(*,*)"nn_index_lon_check", gc_j_lon, 'XXN(i)=',XXN(i),'temp_XXN=',temp_XXN   !20140604-1,20140605-1
!        write(*,*)"nn_index_alt_check", gc_i_alt, 'ZZN(i)=',ZZN(i)   !20140604-1
        !Extract wind values from global wind fields on each wp by using resulting indices.20130813-2,20130905-3
        !If sound speed is also needed, see20140212-3.
        !Calculate speed of sound at waypoints.
        !Reserve wind values along waypoints; they are passed to CALL gc_ac_traj.
        !wind_along_traj(:,1) :u [m/s] ,20140430-1
        !               (:,2) :v [m/s]
        !               (:,3) :w [m/s]
        !vtas_along_traj(:)   :AC_MACH*speed of sound [m/s]
        gc_wind_along_traj(i,1) = gc_uwind_g(gc_j_lon, gc_i_alt, gc_i_lat)    ![m/s]
        gc_wind_along_traj(i,2) = gc_vwind_g(gc_j_lon, gc_i_alt, gc_i_lat)    ![m/s]
        gc_wind_along_traj(i,3) =   gc_v_z_g(gc_j_lon, gc_i_alt, gc_i_lat)    !global vertical wind velocity field [m/s]
        !Yin_20170423
        gc_cpc_along_traj(i)    =   gc_cpc_g(gc_j_lon, gc_i_alt, gc_i_lat)    !global contrail potential coverage [-]
        !Yin_20170423        
!        write(*,*)'gc_cpc_along_traj:',SIZE(gc_cpc_along_traj),gc_cpc_g(gc_j_lon, gc_i_alt, gc_i_lat)

        IF(VTAS_SWITCH.eq.0)THEN   !20140602-1
           gc_vtas_along_traj(i)   = AC_MACH*SQRT(gc_t_scb_g(gc_j_lon, gc_i_alt, gc_i_lat)*   &
                                     adiabatic_index_air*gas_constant_air)  ![m/s] 20140522-4,20140724-2 
        ELSEIF(VTAS_SWITCH.eq.1)THEN
           !AC_V   = AC_MACH*SOUND_V/1000.0_dp    !Aircraft speed [km/s], 20130719-1,20131031-1,-2,20131104-1,20131114-2,20130315-2 
           gc_vtas_along_traj(i)   = AC_MACH*SOUND_V   ![m/s]
        ENDIF
!        write(*,*)"gc_wind_along_traj=",(gc_wind_along_traj(i,j),j=1,3)
!        write(*,*)"gc_vtas_along_traj=",gc_vtas_along_traj(i),"[m/s]"
        !Yin_20170423
!        write(*,*)"gc_cpc_along_traj=",gc_cpc_along_traj(i),"[-]"
        !Yin_20170423
        ! write(*,*)gc_cpc_along_traj(i),"[-]"
        if(gc_wind_along_traj(i,3).gt.1.0_dp)write(*,*)"Check_this_case"  !20140429-2
     ENDDO
     
     !Calculation of objective functions along flight trajectories by CALL gc_ac_traj.
     !CALL gc_ac_traj(THETA_OUT, PHI_OUT, FLIGHT_TIME, FLIGHT_TIME_JULIAN, T_FLIGHT_TIME, A_time, &
     !                FLIGHT_DISTANCE, T_FLIGHT_DISTANCE) 
     !Original-- Calculation of aircraft flight trajectory
     !CALL gc_ac_traj(XXN, YYN, FLIGHT_TIME, FLIGHT_TIME_JULIAN, T_FLIGHT_TIME, A_time, &
     !                FLIGHT_DISTANCE, T_FLIGHT_DISTANCE) 
     !Yin_20170423 gc_cpc_along_traj, gc_cpc are added
     CALL gc_ac_traj(XXN, YYN, ZZN, gc_cpc_along_traj, gc_wind_along_traj, gc_vtas_along_traj,     &     !gc_vtas_along_traj 20140530-1
                     FLIGHT_TIME, FLIGHT_TIME_JULIAN, T_FLIGHT_TIME,  &
                     A_time, FLIGHT_DISTANCE, T_FLIGHT_DISTANCE, FLIGHT_SPEED, &
                     contrail_coverage, gc_cpc)  !FLIGHT_SPEED = ac_gs[m/s]
     !After this CALL, the range of variables are: see opt 20130215-5                                      !CPC 20170602-1

  !Calculation of fuel_use for each segment using Total Energy model 
!  CALL total_energy_model_calculation(FLIGHT_TIME, FUEL_USE) 
  !total_flight_time_s = T_FLIGHT_TIME   !20131003-3,20131114-2,20140601-1
  !write(*,*)'total_flight_time_s=',total_flight_time_s

  !Outputs of flightplan information
!  write(*,'(1x,2a,2x,f10.4,2x,a,f10.4,2x,a,f8.1)')'Departure city     :','lat',gc_D_lat,'lon',gc_D_lon,'Alt',ALT
!  write(*,'(1x,2a,2x,f10.4,2x,a,f10.4,2x,a,f8.1)')'Arrival city       :','lat',gc_A_lat,'lon',gc_A_lon,'Alt',ALT
!  write(*,*)'Departure time[Julian date]  :',gc_D_time
!  write(*,*)'Arrival time[Julian date]    :',A_time
!  write(*,*)'Total flight time    :',INT(T_FLIGHT_TIME),'sec'
!  write(*,'(1x,a,f9.4,a,f6.4,a)')'Flight Speed(Mach)   :',AC_MACH
!  write(*,*)'ISWITCH              :',ISWITCH
!  write(*,*)'OUT_ISWITCH          :',OUT_SWITCH
!  write(*,*)'Flight direction     :',FL_DIR           !FL_DIR: see 20111219-1
! write(*,'(x,a,f7.4,x,a)')'Central angle between D_city and A-city:',COSIS,'radians'
! write(*,'(x,a,2(f10.4,a))')'Flight distance (GC) between D_city and A-city:',T_DIS,'km;',T_DISN,'nm'
  write(*,'(1x,a,2(f10.4,a))')'Flight distance (GC) between D_city and A-city:',T_FLIGHT_DISTANCE,'km;'&
                             ,T_FLIGHT_DISTANCE*0.539957_dp,'nm'
  write(*,'(/,a)')'>>>>>>>>>>>>>>>>> Flight PLAN on GC >>>>>>>>>>>>>>>>>>>'
  write(*,'(a,3x,a,3(5x,a))')'Time[Julian]','Time [sec]','Lat','Lon','Alt'

  !Arrange the outputs range
  SELECT CASE(OUT_SWITCH)
  CASE(0)  
     ! Only for FL_DIR = 1,3, the following DO structure will be adapted. 
     DO i=1,DIV_WAY+1
     !  IF(PHI_OUT(i).gt.180.0_dp)PHI_OUT(i) = PHI_OUT(i) - 360.0_dp
        IF(XXN(i).gt.180.0_dp)XXN(i) = XXN(i) - 360.0_dp
     ENDDO
!     DO j=1,DIV_WAY+1
!        write(*,'(f15.6,2x,i6,2x,f10.4,1x,f10.4,2x,f7.1)')&
!               FLIGHT_TIME_JULIAN(j),INT(FLIGHT_TIME(j)),YYN(j),XXN(j),ZZN(j)
!     ENDDO

  CASE(1)
!     DO j=1,DIV_WAY+1
!        write(*,'(f15.6,2x,i6,2x,f10.4,1x,f10.4,2x,f7.1)')&
!               FLIGHT_TIME_JULIAN(j),INT(FLIGHT_TIME(j)),YYN(j),XXN(j),ZZN(j)
!     ENDDO

  CASE DEFAULT
     write(*,*) 'OUT_SWITCH error: Unknown OUT_SWITCH on GC output!'
  ENDSELECT

  !Outputs for upper SMCL(ac_routes properties)
  gc_fl_time = T_FLIGHT_TIME                                      ![s] 20140601-2 
  DO j=1,DIV_WAY+1
     gc_ac_routes(j,props_lon)           = XXN(j)                 !lon
     gc_ac_routes(j,props_lat)           = YYN(j)                 !lat
     gc_ac_routes(j,props_alt)           = ZZN(j)                 !alt [m]  20140601-3
     gc_ac_routes(j,props_time)          = FLIGHT_TIME_JULIAN(j)  !Passing time at each waypoint, [Julian date]
     gc_ac_routes(j,props_ac_speed)      = FLIGHT_SPEED(j)*3.6_dp !AC ground speed, [km/h] 20140601-2 
     gc_ac_routes(j,props_dist)          = FLIGHT_DISTANCE(j)     !Distance for each segment, [km]  
     gc_ac_routes(j,props_potcov)        = contrail_coverage(j)   !potcov for each segment, [km] 20170602-1  
     !Yin_20170423
!    gc_ac_routes(j,props_potcov)        = gc_cpc_along_traj(j)   !potcov for each segment, [-]  
     !Yin_20170423
!    gc_ac_routes(j,props_fuel_use)      = FUEL_USE(j)            !Fuel consumption for each segment, [kg]  
!    gc_ac_routes(j,props_emis_nox)      = 1.0_dp                 !Dummy for nox_emis 
!    gc_ac_routes(j,props_emis_h2o)      = 1.0_dp                 !Dummy for h2o_emis
!    write(10,'(f10.4,x,f10.4,2x,2(8x,a))')XXN(j),YYN(j),'1.000','0.01'
  ENDDO

  !Outputs for GNUPLOT.
!  OPEN(10,file='gnuplot.dat') 
!    DO i=1,DIV_WAY+1
!       IF(XXN(i).gt.180.0d0)XXN(i) = XXN(i) - 360.0d0
!    ENDDO
!     DO j=1,DIV_WAY+1
!       write(10,'(f7.4,x,f9.4,2x,f7.1)')YYN(j),XXN(j),ALT
!        write(10,'(f15.6,x,f10.4,x,f10.4,2x,2(8x,a))')FLIGHT_TIME_JULIAN(j),XXN(j),YYN(j),'1.000','0.01'
!     ENDDO
!  CLOSE(10)

  !Outputs for GMT. see 20130118-1
!  OPEN(20,file='gmt.dat',position='append')
!     DO j=1, DIV_WAY+1
!        write(20,*)XXN(j),YYN(j)
!     ENDDO
!     write(20,'(a)')">"
!  CLOSE(20)

  !Flight distance check,20131203-1,20141216-1
!  OPEN(30,file="flight_check.dat",position='APPEND') 
     if(T_FLIGHT_DISTANCE.ge.0.0_dp .and. T_FLIGHT_DISTANCE.lt.1000.0_dp)gc_fl_distance = 0
     if(T_FLIGHT_DISTANCE.ge.1000.0_dp .and. T_FLIGHT_DISTANCE.lt.2000.0_dp)gc_fl_distance = 1
     if(T_FLIGHT_DISTANCE.ge.2000.0_dp .and. T_FLIGHT_DISTANCE.lt.3000.0_dp)gc_fl_distance = 2
     if(T_FLIGHT_DISTANCE.ge.3000.0_dp .and. T_FLIGHT_DISTANCE.lt.4000.0_dp)gc_fl_distance = 3
     if(T_FLIGHT_DISTANCE.ge.4000.0_dp .and. T_FLIGHT_DISTANCE.lt.5000.0_dp)gc_fl_distance = 4
     if(T_FLIGHT_DISTANCE.ge.5000.0_dp .and. T_FLIGHT_DISTANCE.lt.6000.0_dp)gc_fl_distance = 5
     if(T_FLIGHT_DISTANCE.ge.6000.0_dp .and. T_FLIGHT_DISTANCE.lt.7000.0_dp)gc_fl_distance = 6
     if(T_FLIGHT_DISTANCE.ge.7000.0_dp .and. T_FLIGHT_DISTANCE.lt.8000.0_dp)gc_fl_distance = 7
     if(T_FLIGHT_DISTANCE.ge.8000.0_dp .and. T_FLIGHT_DISTANCE.lt.9000.0_dp)gc_fl_distance = 8
     if(T_FLIGHT_DISTANCE.ge.9000.0_dp)gc_fl_distance = 9
!  CLOSE(30)
!C test 
!  call cprogram
!  call system('sh /data/yama_hi/EMAC/AirTraf/messy_2.41/messy/smcl/ctext.sh')

  DEALLOCATE(gc_zalt)                       !Array 20140601-3

  END SUBROUTINE gc_calculate
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
     SINDIS1= SQRT((COS(DCRAD2)*SIN(DELRA2))**2.0 + (COS(DCRAD1)*SIN(DCRAD2)-&
                    SIN(DCRAD1)*COS(DCRAD2)*COS(DELRA2))**2.0)
     SINDIS2= SIN(DCRAD1)*SIN(DCRAD2) + COS(DCRAD1)*COS(DCRAD2)*COS(DELRA2) 
     COSIS  = DATAN2(SINDIS1,SINDIS2)        !in [rad]

  CASE DEFAULT
     write(*,*) 'ISWITCH error: Unknown ISWITCH on GC distance calculation!'
  ENDSELECT

! T_DIS    = RADIUS*COSIS                    !GC distance = radius * radians [km]
! T_DIS    = (RADIUS + ALT/1000.0)*COSIS     !GC distance = radius * radians [km] 
  T_DIS    = (radius_earth*0.001_dp + ALT/1000.0_dp)*COSIS     !GC distance = radius * radians [km] 
  T_DISN   = T_DIS*0.539957_dp                  !convert T_DIS into nautical miles [nm]

  END SUBROUTINE gc_distance
!----------------------------------------------------------------------

!----------------------------------------------------------------------
  SUBROUTINE gc_waypoints(PHI_OUT, THETA_OUT)
!----------------------------------------------------------------------
!
! This subroutine calculate waypoints along GC.
! Equations below shoulb be checked on precision!!!
  IMPLICIT NONE
  REAL(DP) :: THETA(DIV_WAY+1), PHI(DIV_WAY+1) 
  REAL(DP), INTENT(OUT) :: THETA_OUT(DIV_WAY+1), PHI_OUT(DIV_WAY+1) 
  REAL(DP) :: DELDC3, DELRA3
  INTEGER i

  !Initialization
  DELDC3 = 0.0_dp
  DELRA3 = 0.0_dp
  THETA  = 0.0_dp
  PHI    = 0.0_dp

  DELRA3 = RARAD1 - RARAD2     !Delta lon between gc_D_lon & gc_A_lon
  DELDC3 = DCRAD1 - DCRAD2     !Delta lat between gc_D_lat & gc_A_lat

  IF((DELRA3.gt.0.0_dp).and.(DELRA3.le.PI))THEN
     DO i=1,DIV_WAY+1
        PHIR = RARAD2 + (i-1)*(RARAD1-RARAD2)/DIV_WAY              !in [rad]  20140722-2
        THETAR = atan(sin(DCRAD1)*sin(PHIR-RARAD2)/(cos(DCRAD1)*sin(RARAD1-RARAD2))+&
                      sin(DCRAD2)*sin(RARAD1-PHIR)/(cos(DCRAD2)*sin(RARAD1-RARAD2)))
        PHI(i) = PHIR/DTR         !lon along wp [deg]
        THETA(i) = THETAR/DTR     !lat along wp [deg]
     ENDDO 
  ELSE IF(DELRA3.gt.PI)THEN
     DO i=1,DIV_WAY+1
        PHIR = RARAD1 + (i-1)*((RARAD2+2.0_dp*PI) - RARAD1)/DIV_WAY               !lon [rad]  20140722-2
        THETAR = atan(sin(DCRAD1)*sin(PHIR-RARAD2)/(cos(DCRAD1)*sin(RARAD1-RARAD2))+&
                      sin(DCRAD2)*sin(RARAD1-PHIR)/(cos(DCRAD2)*sin(RARAD1-RARAD2)))
        PHI(i) = PHIR/DTR         !lon along wp [deg]
        THETA(i) = THETAR/DTR     !lat along wp [deg]
     ENDDO
  ELSE                            !(DELRA3.eq.0.0)
     DO i=1,DIV_WAY+1
        THETAR = DCRAD2 + (i-1)*DELDC3/DIV_WAY            !lon [rad]  20140722-2
        PHI(i) = RARAD1/DTR       !lon along wp [deg]
        THETA(i) = THETAR/DTR     !lat along wp [deg]
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
  SUBROUTINE gc_ac_traj(T_LON, T_LAT, T_ALT, cpc_value, wind_value, ac_vtas,   &
                        FL_T, FL_T_Julian, T_FL_T,                  &
                        T_FL_T_Julian, FL_DIST, T_FL_DIST, ac_gs, ac_cpc, T_potcov) 
!----------------------------------------------------------------------
!
! This subroutine calculate aircraft flight trajectories.
! Equations below shoulc be checked on precision!!!

  IMPLICIT NONE
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: wind_value   ![m/s] 20140531-1
  REAL(DP), DIMENSION(:), INTENT(IN)   :: ac_vtas      ![m/s] 20140531-1
  !Yin_20170423
  REAL(DP), DIMENSION(:), INTENT(IN)   :: cpc_value    ![-] 20170602-1
  !Yin_20170423
  REAL(DP), INTENT(IN)   :: T_LON(DIV_WAY+1), T_LAT(DIV_WAY+1), T_ALT(DIV_WAY+1)
  REAL(DP), INTENT(OUT)  :: FL_T(DIV_WAY+1), FL_T_Julian(DIV_WAY+1), FL_DIST(DIV_WAY+1)
  REAL(DP), INTENT(OUT)  :: ac_gs(DIV_WAY+1)  !20140531-1 ac ground speed including wind effects [m/s]
  REAL(DP), INTENT(OUT)  :: ac_cpc(DIV_WAY+1)  !CPC [km] 20170602-1
  REAL(DP), INTENT(OUT)                :: T_FL_T, T_FL_T_Julian, T_FL_DIST
  !Yin_20170423
  REAL(DP), INTENT(OUT)                :: T_potcov          ![km] 20170602-1  
  !Yin_20170423

  REAL(DP)                             :: W1_lat, W1_lon, W1_alt, W2_lat, W2_lon, W2_alt
  REAL(DP)                             :: RARAD1, RARAD2, DCRAD1, DCRAD2
  REAL(DP)                             :: COSIS, DIS !,DISN 
  REAL(DP)                             :: DELT_T, DELT_T_Julian 
  INTEGER                              :: j1

 !Initilizations
! AC_V            = 0.0_dp 
  FL_DIST(1)      = 0.0_dp 
  ac_cpc(1)       = 0.0_dp   !20170602-1
  ac_gs(DIV_WAY+1)= 0.0_dp    !20140531-1
  T_FL_DIST       = 0.0_dp 
  DELT_T          = 0.0_dp
  DELT_T_Julian   = 0.0_dp     ![Julian date]
  T_FL_T          = 0.0_dp   
  T_FL_T_Julian   = 0.0_dp     ![Julian date]
  DIS             = 0.0_dp
  FL_T(1)         = 0.0_dp 
  FL_T_Julian(1)  = gc_D_time  ![Julian date]
  !Yin_20170423
  T_potcov        = 0.0_dp  ! 20161212, [km]
  !Yin_20170423

  !Yin_20170423
!  write(*,*)'gc_D_time_check:',gc_D_time,'DTR:',DTR,'ISWITCH',ISWITCH
!  write(*,*)'radius_earth:',radius_earth,'ALT',ALT,'OneDay',OneDay
!  write(*,*)"wind_value_check:",wind_value
!  write(*,*)"ac_vtas_check",ac_vtas,"[m/s]"
  
  !Calculate flight time tables from D_city to A_city along the waypoints
  !Waypoints are already in order from D to A city.
  !The range of variables are:see opt20130215-4,20140512-4
  !    -180 =< T_LON =< 180 [FL_DIR=0,2,4,5]
  !       0 =< T_LON =< 360 [FL_DIR=1,3]
  !20140531-1 
  DO j1=1, DIV_WAY
     W1_lat = T_LAT(j1)
     W1_lon = T_LON(j1)
     W1_alt = T_ALT(j1)    !20140531-1
     W2_lat = T_LAT(j1+1)
     W2_lon = T_LON(j1+1)
     W2_alt = T_ALT(j1+1)  !20140531-1
     !Yin_20170423
!     T_potcov = T_potcov+T_CPC(j1+1)
     !Yin_20170423
     !Definition of the waypoint being placed at higher lon between W1 and W2
     IF(W1_lon.lt.W2_lon)THEN
        RARAD1= W2_lon*DTR   !lon of higher waypoint [rad]
        RARAD2= W1_lon*DTR   !lon of lower  waypoint [rad]
        DCRAD1= W2_lat*DTR   !lat of higher waypoint [rad]    
        DCRAD2= W1_lat*DTR   !lat of lower  waypoint [rad]
     ELSE IF(W1_lon.gt.W2_lon)THEN
        RARAD1= W1_lon*DTR     
        RARAD2= W2_lon*DTR
        DCRAD1= W1_lat*DTR     
        DCRAD2= W2_lat*DTR     
     ELSE                    !(W1_lon.eq.W2_lon)
        IF(W1_lat.gt.W2_lat)THEN
           RARAD1= W1_lon*DTR     
           RARAD2= W2_lon*DTR
           DCRAD1= W1_lat*DTR     
           DCRAD2= W2_lat*DTR     
        ELSE IF(W1_lat.lt.W2_lat)THEN
           RARAD1= W2_lon*DTR 
           RARAD2= W1_lon*DTR
           DCRAD1= W2_lat*DTR     
           DCRAD2= W1_lat*DTR 
        ELSE                 !(W1_lat.eq.W2_lat)
           write(*,*)'Waypoints Error: They must be identical waypoints!'
        ENDIF
     ENDIF

     !Calculate central angle between W1 and W2.
     !COSIS: central angle in [rad]
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
        write(*,*) 'ISWITCH error: Unknown ISWITCH on WP distance calculation!'
     END SELECT

     !Calculate flight distance between W1 and W2,
     !DIS does not change with or without wind effect. Wind affects only ac ground speed.
     !add following DIS eq. 20141205-1
     DIS = SQRT(((radius_earth*0.001_dp+W1_alt*0.001_dp)**2) + ((radius_earth*0.001_dp+W2_alt*0.001_dp)**2)  &
                -2.0_dp*(radius_earth*0.001_dp+W1_alt*0.001_dp)*(radius_earth*0.001_dp+W2_alt*0.001_dp)*COS(COSIS))  ![km]
     !Same altitude at W1 and W2
!     write(*,*)"(GC)same alt check",W1_alt, W2_alt, ALT
!    DIS     = (radius_earth*0.001_dp + ALT*0.001_dp)*COSIS      !WP distance = radius * radians [km],20141205-1 
!    DIS     = (RADIUS + ALT/1000.0_dp)*COSIS                     !WP distance = radius * radians [km] 
!    DISN    = DIS*0.539957_dp                                    !Unit conversion DIS into nautical miles [nm]
     
     T_FL_DIST = T_FL_DIST + DIS                                  !Total flight distance from D_city to A_city [km]
     FL_DIST(j1+1) = DIS                                          !Flight distance between W1 and W2 [km]

!    write(*,*)'Distance |W',j1,'-W',j1+1,'| =',DIS,'[km]',DISN,'[km]',' Total flight distance =',T_FL_DIST,'[km]'
!    write(*,*)'ISWITCH =',ISWITCH,'FL_DIR =',FL_DIR
    
     IF(WIND_SWITCH.eq.1)THEN    !20140531-1
        !Calculate ac ground speed considering wind effects at each waypoints.
        !   ac_vtas    :AC true air speed at waypoint1 W1 [m/s]
        !   wind_value :wind values(u,v,w) at waypoint1 W1 [m/s]
        !   ac_gs      :AC ground speed at waypoint1 W1 [m/s]
        CALL wind_effect(W1_lat, W1_lon, W1_alt, W2_lat, W2_lon, W2_alt, DIS, wind_value(j1,:), ac_vtas(j1), ac_gs(j1))
!        write(*,*)"AC_SPEED_CHECK(Wind on)","ac_vtas[m/s]=",ac_vtas(j1),"ac_gs[m/s]=",ac_gs(j1)
     ELSEIF(WIND_SWITCH.eq.0)THEN
        !No wind effects.
        ac_gs(j1) = ac_vtas(j1)    ![m/s] 20140531-1,2
!        write(*,*)"AC_SPEED_CHECK(Wind off)","ac_vtas[m/s]=",ac_vtas(j1),"ac_gs[m/s]=",ac_gs(j1)
        !AC_V   = AC_MACH*SOUND_V/1000.0_dp                        !Aircraft speed [km/s] 
!       FL_SPEED(j1) = AC_V  ! <<<<<<-- input wind fields,  AC_V[km/s], i.e., FL_SPEED[km/s] 
!HY     FL_SPEED(j1) = ac_gs*0.001_dp  ! FL_SPEED[km/s]
     ENDIF
     
     !Calculate flight time between W1 and W2
     !DELT_T  : Flight time between W1 and W2 [s]  20140523-1
     !FL_DIST : Flight distance between W1 and W2 [km]
     !ac_gs   : AC ground speed[m/s] = ac_vtas + wind_speed 20140523-1
     !AC_V   = AC_MACH*SOUND_V/1000.0_dp                   !Aircraft speed [km/s] 
     !DELT_T = DIS/AC_V                                    !Flight time between W1 and W2 [s]
     DELT_T = FL_DIST(j1+1)/(ac_gs(j1)*0.001_dp)           !20140601-1
     FL_T(j1+1)  = FL_T(j1) + DELT_T                       !Aircraft passing time at each waypoint [s],see20130513-1
     if(j1.eq.DIV_WAY)T_FL_T = FL_T(j1+1)                  !Total flight time from D_city to A_city [s] 
     !write(*,*)'DELT_T_CHECK=',DELT_T,'j1=',j1
!     write(*,*)'j1=',j1,'W1=',W1_alt,'W2=',W2_alt,"FL_DIST=",FL_DIST(j1+1),'DELT_T=',DELT_T,"FL_T=",FL_T(j1+1)

     !Flight time [Julian date]
     DELT_T_Julian = DELT_T/OneDay                         !Flight time between W1 and W2 [Julian date]
     FL_T_Julian(j1+1) = FL_T_Julian(j1) + DELT_T_Julian   !Aircraft passing time at each waypoint [Julian date]
     if(j1.eq.DIV_WAY)T_FL_T_Julian = FL_T_Julian(j1+1)    !Arrival time at A_city [Julian date] 

     !Fuel_use calculation see20130521-3
     ! FUEL_U(j1+1)=   *DELT_T                             !DELT_T[s]  
     ! FUEL_U(j1+1)=    DELT_T                             !DELT_T[s]  

     !Calculate Total CPC, 20170602-1
     !cpc_value  :potential contrail coverage at waypoint [-]
     !ac_cpc     :potential contrail coverage for each segment [km]
     !T_potcov   :total potential contrail coverage [km]
     ac_cpc(j1+1) = cpc_value(j1)*FL_DIST(j1+1)                !potcov between W1 and W2 [km]
     T_potcov = T_potcov + ac_cpc(j1+1)                    ![km] 

     !Initialization for next waypoints calculation
     DIS             = 0.0_dp
     DELT_T          = 0.0_dp   
     DELT_T_Julian   = 0.0_dp   
  ENDDO   !End loop for waypoint

  END SUBROUTINE gc_ac_traj
!----------------------------------------------------------------------

!----------------------------------------------------------------------
  SUBROUTINE total_energy_model_calculation(flight_time_table, FUEL_U)
!----------------------------------------------------------------------
! This subroutine calculates as_routes(fuel_use) using Total Energy model.
!
!
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: flight_time_table(DIV_WAY+1) 
  REAL(DP), INTENT(OUT):: FUEL_U(DIV_WAY+1)

  REAL(DP) :: ita,Cfcr,CD2,g,rho,S_area,ac_mass,new_ac_mass,CD0,VTAS,Cf1,Cf2
  REAL(DP) :: aa,bb,cc,sq,de,x1,x2     !,CL,CD,DRAG_FORCE,THRUST_FORCE 
  INTEGER j1


! Initialization
 FUEL_U(1)       = 0.0_dp 
 ita             = 0.0_dp !0.00001566666  !0.94_dp
 Cfcr            = 0.93655_dp ![-]
 g               = 9.81_dp   ![m/s^2]
 rho             = 0.4_dp    ![kg/m^3]
 S_area          = 361.6_dp  ![m^2]
 ac_mass         = 174000.0_dp ![kg]
 CD0             = 0.019805_dp ![-]
 CD2             = 0.031875_dp ![-]
! CL             = 0.0_dp
! CD             = 0.0_dp
 VTAS            = 250.0     ![m/s]
 Cf1             = 0.61503_dp ![kg/min/kN]
 Cf2             = 919.03_dp  ![Kt]

 aa              = 0.0_dp
 bb              = 0.0_dp
 cc              = 0.0_dp
 de              = 0.0_dp
 sq              = 0.0_dp
 x1              = 0.0_dp
 x2              = 0.0_dp

! The below 'ita' is constant during flight with GC option. 
! If wind case, this ita should be calculated within DO sentence below.
! Because VTAS changes due to wind velocity. see 20130521-3
  ita = Cf1*(1.0_dp + VTAS/0.514443004_dp/Cf2)/60.0_dp/1000.0_dp  ![kg/s/N]


! This DO-calculation is executed backward from A- to D-city.
  DO j1=DIV_WAY+1,2,-1
!    CL= 2.0_dp*g*ac_mass/rho/VTAS**2/S_area  
!    CD= CD0 + CD2*CL**2
!    DRAG_FORCE=0.5_dp*rho*VTAS**2*S_area*CD
!    THRUST_FORCE=DRAG_FORCE*0.001_dp
!    write(*,*)'CL=',CL,'CD=',CD,'AC_V*1000=',AC_V*1000.0_dp,'VTAS=',VTAS,'DRAG_FORCE=',DRAG_FORCE,'THUST_FORCE=',THRUST_FORCE
     
!     write(*,*)'ac_mass=',ac_mass,'ita=',ita,'j1=',j1,'time=',flight_time_table(j1)-flight_time_table(j1-1) 
     aa = ita*Cfcr*(flight_time_table(j1)-flight_time_table(j1-1))*CD2*2.0_dp*(g**2)/VTAS**2/rho/S_area
     bb = -1.0_dp
     cc = ac_mass + ita*Cfcr*(flight_time_table(j1)-flight_time_table(j1-1))*rho*(VTAS**2)*S_area*CD0*0.5_dp
!     write(*,*)'aa=',aa,'bb=',bb,'cc=',cc

     sq = dsqrt(bb*bb-4.0_dp*aa*cc)
     de = 2.0_dp*aa
     x1 = (-bb+sq)/de
     x2 = -(bb+sq)/de
!     write(*,*)'x1=',x1,'x2=',x2

     !Determine a solution as new_ac_mass, x1 or x2 ?
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
     !for next waypoint calculation
     FUEL_U(j1) = new_ac_mass - ac_mass  ![kg]
     ac_mass    = new_ac_mass            ![kg]
     
  ENDDO

  END SUBROUTINE total_energy_model_calculation
!----------------------------------------------------------------------

!----------------------------------------------------------------------
  SUBROUTINE wind_effect(W1_lat_in, W1_lon_in, W1_alt_in, W2_lat_in, W2_lon_in, W2_alt_in, DIS_in, &
                         wind_value_in, AC_V_in, new_AC_V_in)   !20140531-2
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
!  write(*,*)"(GC)wind_components=",ux_wind_component,"[m/s]",vy_wind_component,"[m/s]",wz_wind_component,"[m/s]"  

  END SUBROUTINE wind_effect 
!----------------------------------------------------------------------

END MODULE messy_airtraf_gc

