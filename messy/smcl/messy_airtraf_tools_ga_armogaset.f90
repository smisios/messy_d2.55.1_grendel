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
! - DVs Output ranges : -180 =< dv1,dv3,dv5 =< +180 [deg] (FL_DIR=0,2,4,5)
!                          0 =< dv1,dv3,dv5 =< +360 [deg] (FL_DIR=1,3) 
!
! - Flight time is accounted and integrated by 'FLIGHT_TIME(:)' at every waypoint in [sec]. 
!   The 'T_FLIGHT_TIME(DP)' indicates a total flight time from D_city to A_city in [sec].
!   Flight time and flight distance are calculated in spherical coordinates.
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
!HY20140401-2 Program messy_airtraf_gc
MODULE messy_airtraf_tools_ga_armogaset
! **********************************************************************
 USE messy_main_constants_mem, ONLY: DP, PI, DTR   !HY20140422-2,20140722-4

  IMPLICIT NONE
  ! GLOBAL PARAMETERS
  ! Comments:If we can use the information on aircraft type for each route, 
  !         average flight speed of each aircraft type can be utilized as AC_MACH. 
! REAL(DP), PARAMETER :: pi           = 3.14159265358979323846_dp
! REAL(DP), PARAMETER :: OneDay       = 86400.0_dp   ! one day [s]
! REAL(DP), PARAMETER :: radius_earth = 6371000.0_dp   !Radius of earth in km 
!HY  REAL(DP), PARAMETER :: ALT          = 0.0_dp !10000.0d0   !Altitude under steady flight, [m] 
! REAL(DP), PARAMETER :: AC_MACH      = 0.82_dp !0.710157704d0 !0.8d0(opt,20130131-1)  !Mach number in steady flight, [-]  
! REAL(DP), PARAMETER :: SOUND_V      = 340.3_dp     !Speed of sound [m/s] 
! REAL(DP), PARAMETER :: DTR          = pi/180.0_dp    ! Degrees to radians 
!HY20140428-1  REAL(DP), PARAMETER :: delta_x      = 0.05_dp  !5.0_dp    ! Degrees to radians 
!HY20140428-1  REAL(DP), PARAMETER :: delta_y      = 2.0_dp  !15.0_dp    ! Degrees to radians 
  REAL(DP), PARAMETER :: percent_x_region      = 0.1_dp  ! DVs region in lon 
  REAL(DP), PARAMETER :: percent_y_region      = 0.3_dp  ! DVs region in lat 
  REAL(DP), PARAMETER :: Lower_alt    = 8839.2_dp    ! corresponding to 29000ft, in [m] 
  REAL(DP), PARAMETER :: Upper_alt    = 12496.8_dp    ! corresponding to 41000ft, in [m] 

  ! The parameter should be used from SMIL
  INTEGER, PARAMETER, PRIVATE :: D_lon=1, D_lat=2, D_time=3, A_lon=4, A_lat=5  !HY20140418-3
  INTEGER, PARAMETER, PRIVATE :: DIV_WAY = 100  !HY20140425-1,this DIV_WAY is independent from another DIV_WAY.20140722-3 
  INTEGER, PARAMETER          :: NUM_CP_LONLAT  = 3

  ! PRIVATE PARAMETERS
! INTEGER, PARAMETER :: ISWITCH = 2                   !0: Spherical law of cosines
                                                      !1: Haversine formula
                                                      !2: Vincenty formula
! INTEGER, PARAMETER :: OUT_SWITCH = 0                !0: -180 =< lon =< +180 (Default) 
                                                      !1:    0 =< lon =< +360 (only cases FL_DIR = 1 and 3)

  ! FROM MUN TO JFK
!HY20140418-2  DOUBLE PRECISION, PARAMETER :: D_lat = 48.35388889d0       !lat of departure-city 
!HY20140418-2  DOUBLE PRECISION, PARAMETER :: D_lon = 11.78611111d0      !lon of departure-city 
!HY20140418-2  DOUBLE PRECISION, PARAMETER :: A_lat = 40.63972222d0       !lat of arrival-city
!HY20140418-2  DOUBLE PRECISION, PARAMETER :: A_lon = -73.77888889d0      !lon of arrival-city

  !To determine DVS region
  REAL(DP)                       :: delta_x, delta_y  !HY20140428-1

  ! PRIVATE
  REAL(DP)                       :: THETA_OUT(DIV_WAY+1), PHI_OUT(DIV_WAY+1)
! REAL(DP)                       :: AC_V 
  REAL(DP)                       :: RARAD1,RARAD2,DCRAD1,DCRAD2 
  REAL(DP)                       :: DELDC2,DELRA,DELRA2,SINDIS,SINDIS1,SINDIS2
  REAL(DP)                       :: THETAR, PHIR
  INTEGER                        :: FL_DIR
! CHARACTER(len=3)               :: C1,C2,C3,C4,C5,C6


!HY20140425-2  REAL(DP), DIMENSION(:) :: FLIGHT_TIME(DIV_WAY+1), FLIGHT_TIME_JULIAN(DIV_WAY+1)
!HY20140425-2  REAL(DP)               :: T_FLIGHT_TIME, A_time, FLIGHT_DISTANCE !A_time in Julian date
  REAL(DP)               :: gc_D_lat        !latitude of departure city 
  REAL(DP)               :: gc_D_lon        !longitude of departure city
  REAL(DP)               :: gc_A_lat        !latitude of arrivar city 
  REAL(DP)               :: gc_A_lon        !longitude of arrival city


! REAL(DP), DIMENSION(:) :: XXN(DIV_WAY+1),YYN(DIV_WAY+1),ZZN(DIV_WAY+1)

  ! Parameter for Single Opbjective Optimization Problem on 'MINIMUM FLIGHT TIME'. 
  ! These dvs are used to express arbitray trajectories.
!HY20140401-3  DOUBLE PRECISION :: L_dv1, U_dv1, L_dv2, U_dv2, L_dv3, U_dv3, L_dv4, U_dv4, L_dv5, U_dv5, L_dv6, U_dv6
!HY20140401-3  DOUBLE PRECISION :: L_dv7, U_dv7, L_dv8, U_dv8, L_dv9, U_dv9, L_dv10, U_dv10, L_dv11, U_dv11
  REAL(DP) :: dv1,dv2,dv3,dv4,dv5,dv6
! DOUBLE PRECISION :: DIV_DV 

 
  PUBLIC :: calc_boundary  !HY20140401-3 
  !HY20140401-2
  CONTAINS 
!-----------------------------------------------------------------
  SUBROUTINE calc_boundary(ndv_val, phlb_val, phub_val, calc_p2_ac_routes_desc) !HY20140418-2
!-----------------------------------------------------------------
  IMPLICIT NONE
  INTEGER, INTENT(IN)                         :: ndv_val
  REAL(DP),               INTENT(OUT)         :: phlb_val(ndv_val),phub_val(ndv_val)
  REAL(DP), DIMENSION(:), INTENT(IN)          :: calc_p2_ac_routes_desc  !HY20140418-2
  !READ armoga.set
!HY20140401-3  open(10,file='armoga.set')
!HY     read(10,*)
!HY     read(10,*)
!HY     read(10,*)
!HY     read(10,*)L_dv1, U_dv1
!HY     read(10,*)L_dv2, U_dv2
!HY     read(10,*)L_dv3, U_dv3
!HY     read(10,*)L_dv4, U_dv4
!HY     read(10,*)L_dv5, U_dv5
!HY     read(10,*)L_dv6, U_dv6
!HY     read(10,*)L_dv7, U_dv7
!HY     read(10,*)L_dv8, U_dv8
!HY     read(10,*)L_dv9, U_dv9
!HY     read(10,*)L_dv10, U_dv10
!HY     read(10,*)L_dv11, U_dv11
!HY     read(10,*)
!HY  close(10)

!  write(*,*)'armoga.set_check'
!  write(*,*)
!  write(*,*)
!  write(*,*)
!  write(*,*)L_dv1, U_dv1
!  write(*,*)L_dv2, U_dv2
!  write(*,*)L_dv3, U_dv3
!  write(*,*)L_dv4, U_dv4
!  write(*,*)L_dv5, U_dv5
!  write(*,*)L_dv6, U_dv6
!  write(*,*)L_dv7, U_dv7
!  write(*,*)L_dv8, U_dv8
!  write(*,*)L_dv9, U_dv9
!  write(*,*)L_dv10, U_dv10
!  write(*,*)L_dv11, U_dv11
!  write(*,*)
!stop
  !Initialization
!HY20140425-2  T_FLIGHT_TIME = 0.0_dp
! AC_V          = AC_MACH*SOUND_V/1000.0_dp    !Aircraft speed [km/s] 

  !Definition of D_city and A_city
  !HY20140418-2
  gc_D_lon  = calc_p2_ac_routes_desc(D_lon)        !D_lon        
  gc_D_lat  = calc_p2_ac_routes_desc(D_lat)        !D_lat         
  gc_A_lon  = calc_p2_ac_routes_desc(A_lon)        !A_lon      
  gc_A_lat  = calc_p2_ac_routes_desc(A_lat)        !A_lat        

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

  !Definitions of dv1, dv3, dv5
!  DIV_DV = (gc_D_lon-gc_A_lon)/(NUM_CP_LONLAT+1)
!  dv1    =  gc_D_lon- DIV_DV
!  dv3    =  gc_D_lon-(DIV_DV*2.0d0)
!  dv5    =  gc_D_lon-(DIV_DV*3.0d0)
!  write(*,*)dv1,dv3,dv5

  !Calculation of GC distance
  !HY 20120810
  !CALL gc_distance(COSIS, T_DIS, T_DISN)

!  CALL arbitrary_traj(XXN, YYN, ZZN) 

  !Calculation of waypoints on GC
  CALL gc_waypoints(THETA_OUT, PHI_OUT)
  !See opt20130204-3
  !After this call, in cases of FL_DIR=0,2,4,5, -180 =< PHI_OUT(LON) =< 180
  !                 in cases of FL_DIR=1,3    ,    0 =< PHI_OUT(LON) =< 360

  !Calculation of objective functions along flight trajectories
!GC CALL gc_ac_traj(THETA_OUT, PHI_OUT, FLIGHT_TIME, FLIGHT_TIME_JULIAN, T_FLIGHT_TIME, A_time, FLIGHT_DISTANCE) 
!  CALL gc_ac_traj(YYN, XXN, FLIGHT_TIME, FLIGHT_TIME_JULIAN, T_FLIGHT_TIME, A_time, FLIGHT_DISTANCE) 

  !Calclations of dv2, dv4, dv6
!  CALL determine_domains(THETA_OUT, PHI_OUT, dv1, dv3, dv5, dv2, dv4, dv6)
  CALL determine_domains(THETA_OUT, PHI_OUT, dv1, dv2, dv3, dv4, dv5, dv6)

  !Outputs of flightplan information
!  write(*,'(x,2a,2x,f10.4,2x,a,f10.4,2x,a,f8.1)')'Departure city     :','lat',gc_D_lat,'lon',gc_D_lon,'Alt',ALT
!  write(*,'(x,2a,2x,f10.4,2x,a,f10.4,2x,a,f8.1)')'Arrival city       :','lat',gc_A_lat,'lon',gc_A_lon,'Alt',ALT
!  write(*,*)'Total flight time    :',T_FLIGHT_TIME,'sec'
!  write(*,'(x,a,f9.4,a,f6.4,a)')'Flight Speed(Mach)   :',AC_V*3600.0d0,'km/h(',AC_MACH,')'
!  write(*,*)'ISWITCH              :',ISWITCH
!  write(*,*)'OUT_ISWITCH          :',OUT_SWITCH
!  write(*,*)'Flight direction     :',FL_DIR
!! write(*,'(x,a,f7.4,x,a)')'Central angle between D_city and A-city:',COSIS,'radians'
!! write(*,'(x,a,2(f10.4,a))')'Flight distance (GC) between D_city and A-city:',T_DIS,'km;',T_DISN,'nm'
!  write(*,'(x,a,f10.4,a)')'Flight distance (GC) between D_city and A-city:',FLIGHT_DISTANCE,'km'
!!GC  write(*,'(/,a)')'>>>>>>>>>>>>>>>>> Flight PLAN on GC >>>>>>>>>>>>>>>>>>>'
!  write(*,'(/,a)')'>>>>>>>>>>>>>>>>> Flight PLAN on traj set by current dvs >>>>>>>>>>>>>>>>>>>'
!  write(*,'(a,3x,a,3(5x,a))')'Time[Julian]','Time [sec]','Lat','Lon','Alt'

  !Arrange the outputs range
!  SELECT CASE(OUT_SWITCH)
!  CASE(0)  
     ! Only for FL_DIR = 1,3, the following DO structure will be adapted. 
!     DO i=1,DIV_WAY+1
!        IF(PHI_OUT(i).gt.180.0_dp)PHI_OUT(i) = PHI_OUT(i) - 360.0_dp     !
!        write(*,*)PHI_OUT(i),THETA_OUT(i)
!     ENDDO
!     DO j=1, DIV_WAY+1
!        write(*,'(f15.6,2x,i6,2x,f10.4,x,f10.4,2x,f7.1)')&
!            FLIGHT_TIME_JULIAN(j),INT(FLIGHT_TIME(j)),THETA_OUT(j),PHI_OUT(j),ALT     !
!note!               FLIGHT_TIME_JULIAN(j),INT(FLIGHT_TIME(j)),YYN(j),XXN(j),ZZN(j)     !
!         write(*,*)PHI_OUT(j),THETA_OUT(j)
!     ENDDO

!  CASE(1)
!     DO j=1, DIV_WAY+1
!        write(*,'(f15.6,2x,i6,2x,f10.4,x,f10.4,2x,f7.1)')&
!            FLIGHT_TIME_JULIAN(j),INT(FLIGHT_TIME(j)),THETA_OUT(j),PHI_OUT(j),ALT     !
!note!               FLIGHT_TIME_JULIAN(j),INT(FLIGHT_TIME(j)),YYN(j),XXN(j),ZZN(j)     !
!         write(*,*)PHI_OUT(j),THETA_OUT(j)
!     ENDDO
!  CASE DEFAULT
!     write(*,*) 'OUT_SWITCH error: Unknown OUT_SWITCH on GC output!'
!  ENDSELECT

  !Outputs array for upper SMCL
!  DO j=1, DIV_WAY+1
!     gc_ac_routes(j,D_lon) = PHI_OUT(j)             !lon
!     gc_ac_routes(j,D_lat) = THETA_OUT(j)           !lat
!     gc_ac_routes(j,D_time) = FLIGHT_TIME_JULIAN(j)  !time
!     gc_ac_routes(j,4) = 1.0d0                    !Dummy
!     gc_ac_routes(j,5) = 1.0d0                    !Dummy
!    write(10,'(f10.4,x,f10.4,2x,2(8x,a))')PHI_OUT(j),THETA_OUT(j),'1.000','0.01'
!  ENDDO

  !Outputs for GNUPLOT and GMT.
!  OPEN(10,file='gnuplot.dat') 
!!    DO i=1,DIV_WAY+1
!!       IF(PHI_OUT(i).gt.180.0d0)PHI_OUT(i) = PHI_OUT(i) - 360.0d0
!!    ENDDO
!     DO j=1, DIV_WAY+1
!!       write(10,'(f7.4,x,f9.4,2x,f7.1)')THETA_OUT(j),PHI_OUT(j),ALT
!!GC     write(10,'(f15.6,x,f10.4,x,f10.4,2x,2(8x,a))')FLIGHT_TIME_JULIAN(j),PHI_OUT(j),THETA_OUT(j),'1.000','0.01'   !
!        write(10,'(f15.6,x,f10.4,x,f10.4,2x,2(8x,a))')FLIGHT_TIME(j),XXN(j),YYN(j),'1.000','0.01'   !
!     ENDDO
!  CLOSE(10)
 
  ! Calulate all dvs domains   
  !HY20140401-3,4, 20140426-1,20140428-2
  !lon and lat dir.
!  SELECT CASE(FL_DIR)
!  CASE(0,1,2,3)
     phlb_val(1) = dv1 - delta_x 
     phub_val(1) = dv1 + delta_x
     phlb_val(2) = dv2 - delta_y
     phub_val(2) = dv2 + delta_y
     phlb_val(3) = dv3 - delta_x
     phub_val(3) = dv3 + delta_x
     phlb_val(4) = dv4 - delta_y
     phub_val(4) = dv4 + delta_y
     phlb_val(5) = dv5 - delta_x
     phub_val(5) = dv5 + delta_x
     phlb_val(6) = dv6 - delta_y
     phub_val(6) = dv6 + delta_y
!  CASE(4,5)
!     phlb_val(1) = dv1 - delta_y   !In cases FI_DIR=4 & 5, delta_x and delta_y are used in Lat and Lon dir., respectively. 
!     phub_val(1) = dv1 + delta_y   !HY20140426-1, 20140428-2
!     phlb_val(2) = dv2 - delta_x
!     phub_val(2) = dv2 + delta_x
!     phlb_val(3) = dv3 - delta_y
!     phub_val(3) = dv3 + delta_y
!     phlb_val(4) = dv4 - delta_x
!     phub_val(4) = dv4 + delta_x
!     phlb_val(5) = dv5 - delta_y
!     phub_val(5) = dv5 + delta_y
!     phlb_val(6) = dv6 - delta_x
!     phub_val(6) = dv6 + delta_x
!  END SELECT
  !Altitude dir.
  phlb_val(7) = Lower_alt    
  phub_val(7) = Upper_alt   
  phlb_val(8) = Lower_alt  
  phub_val(8) = Upper_alt 
  phlb_val(9) = Lower_alt  
  phub_val(9) = Upper_alt  
  phlb_val(10)= Lower_alt  
  phub_val(10)= Upper_alt  
  phlb_val(11)= Lower_alt  
  phub_val(11)= Upper_alt  

  !CHECK DVS regions
  !write(*,*)"CHECK_DVS_regions"
  !do i=1,6
  !   write(*,*)phlb_val(i),phub_val(i)
  !enddo

  !WRITE atmoga.set
!HY20140401-3  open(20,file='armoga.set')
!HY     write(20,*)
!HY     write(20,*)
!HY     write(20,*)
!HY     write(20,*)L_dv1, U_dv1
!HY     write(20,*)L_dv2, U_dv2
!HY     write(20,*)L_dv3, U_dv3
!HY     write(20,*)L_dv4, U_dv4
!HY     write(20,*)L_dv5, U_dv5
!HY     write(20,*)L_dv6, U_dv6
!HY     write(20,*)L_dv7, U_dv7
!HY     write(20,*)L_dv8, U_dv8
!HY     write(20,*)L_dv9, U_dv9
!HY     write(20,*)L_dv10, U_dv10
!HY     write(20,*)L_dv11, U_dv11
!HY     write(20,*)
!HY  close(20)

  !check dvs
!  write(*,*)'-------------------'
!  write(*,*)'Center of DVs are:'
!  write(*,*)gc_D_lon,gc_D_lat
!  write(*,*)dv1,dv2
!  write(*,*)dv3,dv4
!  write(*,*)dv5,dv6
!  write(*,*)gc_A_lon,gc_A_lat

!HY20140401-2  CONTAINS
 END SUBROUTINE calc_boundary 
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!  SUBROUTINE gc_distance(COSIS, T_DIS, T_DISN)
!----------------------------------------------------------------------
!
! This subroutine calcuate great circre distance.
!
!  IMPLICIT NONE
!!HY20140522-2  REAL(DP), INTENT(OUT) :: COSIS, T_DIS, T_DISN 
!  REAL(DP)  :: COSIS, T_DIS, T_DISN 

!  SELECT CASE (ISWITCH)  
!  CASE(0)
!     DELRA  = RARAD1-RARAD2
!     COSIS  = ACOS(SIN(DCRAD1)*SIN(DCRAD2) + COS(DCRAD1)*COS(DCRAD2)*&
!                   COS(DELRA))               !in [rad]

!  CASE(1) 
!     DELDC2 = (DCRAD1-DCRAD2)/2.0_dp            !delta lat
!     DELRA2 = (RARAD1-RARAD2)/2.0_dp            !delta lon
!     SINDIS = SQRT(SIN(DELDC2)*SIN(DELDC2) + COS(DCRAD1)*COS(DCRAD2)*&
!                   SIN(DELRA2)*SIN(DELRA2))
!     COSIS  = 2.0_dp*DASIN(SINDIS)              !in [rad]

!  CASE(2)
!     DELRA2 = RARAD1-RARAD2                  !delta lon
!     SINDIS1= SQRT((COS(DCRAD2)*SIN(DELRA2))**2 + (COS(DCRAD1)*SIN(DCRAD2)-&
!                    SIN(DCRAD1)*COS(DCRAD2)*COS(DELRA2))**2)
!     SINDIS2= SIN(DCRAD1)*SIN(DCRAD2) + COS(DCRAD1)*COS(DCRAD2)*COS(DELRA2) 
!     COSIS  = DATAN2(SINDIS1,SINDIS2)        !in [rad]

!  CASE DEFAULT
!     write(*,*) 'ISWITCH error: Unknown ISWITCH on GC distance calculation!'
!  ENDSELECT

!! T_DIS    = RADIUS*COSIS                    !GC distance = radius * radians [km]
!! T_DIS    = (RADIUS + ALT/1000.0)*COSIS     !GC distance = radius * radians [km] 
!!HY20140522-1  T_DIS    = (radius_earth*0.001_dp + ALT/1000.0_dp)*COSIS     !GC distance = radius * radians [km] 
!  T_DISN   = T_DIS*0.539957_dp                  !convert T_DIS into nautical miles [nm]

!  END SUBROUTINE gc_distance
!----------------------------------------------------------------------

!----------------------------------------------------------------------
  SUBROUTINE arbitrary_traj(XX, YY, ZZ) 
!----------------------------------------------------------------------

 IMPLICIT NONE
 REAL(DP) :: XG(NUM_CP_LONLAT+2), YG(NUM_CP_LONLAT+2), ZG(NUM_CP_LONLAT+2)
 REAL(DP) :: XN(DIV_WAY+1),YN(DIV_WAY+1)!,ZN(DIV_WAY+1),SN(DIV_WAY+1)
 REAL(DP),INTENT(OUT) :: XX(DIV_WAY+1),YY(DIV_WAY+1),ZZ(DIV_WAY+1)
 INTEGER :: i,j

  ! Initilization for GRCV3D
  XG = 0.0_dp
  YG = 0.0_dp
  ZG = 0.0_dp
  XN = 0.0_dp
  YN = 0.0_dp
  !ZN = 0.0_dp
  !SN = 0.0_dp
  XX = 0.0_dp
  YY = 0.0_dp
  ZZ = 0.0_dp

  XG(1)=D_lon
  XG(2)=dv1
  XG(3)=dv3
  XG(4)=dv5
  XG(5)=A_lon

  YG(1)=D_lat
  YG(2)=dv2
  YG(3)=dv4
  YG(4)=dv6
  YG(5)=A_lat

!20140522-1  ZG(1)=ALT
!  ZG(2)=ALT
!  ZG(3)=ALT
!  ZG(4)=ALT
!  ZG(5)=ALT

  do i=1, NUM_CP_LONLAT+2
     write(*,*)'arbitrary_subroutine',XG(i),YG(i),ZG(i),NUM_CP_LONLAT+2
  enddo

!  CALL grcv3d(5,XG,YG,ZG,DIV_WAY+1,1,1.0,1.0,XN,YN,ZN,SN)
!  CALL grcv3d(5,XG,YG,ZG,81,1,1.0,1.0,XN,YN,ZN,SN)

  write(*,*)'after_grcv3d'
  do j=1, DIV_WAY+1
     XX(j) = XN(j)*1.0_dp
     YY(j) = YN(j)*1.0_dp
!HY20140522-2     ZZ(j) = ALT !ZN(j)
     write(*,*)XX(j),YY(j),ZZ(j)
  enddo

  END SUBROUTINE arbitrary_traj 
!----------------------------------------------------------------------

!----------------------------------------------------------------------
  SUBROUTINE gc_waypoints(THETA_OUT, PHI_OUT)
!----------------------------------------------------------------------
!
! This subroutine calculate waypoints along the great circle.
!
  IMPLICIT NONE
  REAL(DP)              :: THETA(DIV_WAY+1), PHI(DIV_WAY+1) 
  REAL(DP), INTENT(OUT) :: THETA_OUT(DIV_WAY+1), PHI_OUT(DIV_WAY+1) 
  REAL(DP)              :: DELDC3, DELRA3
  INTEGER               :: i

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
        THETAR = DATAN(DSIN(DCRAD1)*DSIN(PHIR-RARAD2)/(DCOS(DCRAD1)*DSIN(RARAD1-RARAD2))+&
                       DSIN(DCRAD2)*DSIN(RARAD1-PHIR)/(DCOS(DCRAD2)*DSIN(RARAD1-RARAD2)))
        PHI(i) = PHIR/DTR         !lon along wp [deg]
        THETA(i) = THETAR/DTR     !lat along wp [deg]
     ENDDO 
  ELSE IF(DELRA3.gt.PI)THEN
     DO i=1,DIV_WAY+1
        PHIR = RARAD1 + (i-1)*((RARAD2+2.0_dp*PI) - RARAD1)/DIV_WAY     !lon [rad]
        THETAR = DATAN(DSIN(DCRAD1)*DSIN(PHIR-RARAD2)/(DCOS(DCRAD1)*DSIN(RARAD1-RARAD2))+&
                       DSIN(DCRAD2)*DSIN(RARAD1-PHIR)/(DCOS(DCRAD2)*DSIN(RARAD1-RARAD2)))
        PHI(i) = PHIR/DTR         !lon along wp [deg]
        THETA(i) = THETAR/DTR     !lat along wp [deg]
     ENDDO
  ELSE                            !(DELRA3.eq.0.0)
     DO i=1, DIV_WAY+1
        THETAR = DCRAD2 + (i-1)*DELDC3/DIV_WAY  !lon [rad]
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
!  SUBROUTINE gc_ac_traj(T_LAT, T_LON, FL_T, FL_T_Julian, T_FL_T, T_FL_T_Julian, FL_DIST) 
!----------------------------------------------------------------------
!
! This subroutine calculate aircraft flight trajectories.
!
!  IMPLICIT NONE
!  REAL(DP), DIMENSION(:), INTENT(IN) :: T_LAT(DIV_WAY+1), T_LON(DIV_WAY+1) 
!!  REAL, DIMENSION(:), INTENT(IN) :: T_LAT(DIV_WAY+1), T_LON(DIV_WAY+1) 
!  REAL(DP), DIMENSION(:), INTENT(OUT):: FL_T(DIV_WAY+1), FL_T_Julian(DIV_WAY+1)
!  REAL(DP), INTENT(OUT)              :: T_FL_T, T_FL_T_Julian, FL_DIST

!  REAL(DP) :: W1_lat, W1_lon, W2_lat, W2_lon
!  REAL(DP) :: RARAD1, RARAD2, DCRAD1, DCRAD2
!  REAL(DP) :: COSIS, DIS !,FL_DIST, DISN 
!  REAL(DP) :: DELT_T, DELT_T_Julian 
!  INTEGER  :: j1

 !Initilizations
!  FL_DIST         = 0.0_dp 
!  DELT_T          = 0.0_dp
!  DELT_T_Julian   = 0.0_dp     !Julian date
!  T_FL_T          = 0.0_dp   
!  T_FL_T_Julian   = 0.0_dp     !Julian date
!  DIS             = 0.0_dp
!  FL_T(1)         = 0.0_dp
!  FL_T_Julian(1)  = 0.0_dp      !Julian date


  !Calculate flight time tables from D_city to A_city along the waypoints
!  DO j1=1, DIV_WAY
!     W1_lat = T_LAT(j1)
!     W1_lon = T_LON(j1)
!     W2_lat = T_LAT(j1+1)
!     W2_lon = T_LON(j1+1)
     !Definition of the city being placed at higher lon between W1 and W2
!     IF(W1_lon.lt.W2_lon)THEN
!        RARAD1= W2_lon*DTR   !lon of higher waypoint [rad]
!        RARAD2= W1_lon*DTR   !lon of lower  waypoint [rad]
!        DCRAD1= W2_lat*DTR   !lat of higher waypoint [rad]    
!        DCRAD2= W1_lat*DTR   !lat of lower  waypoint [rad]
!     ELSE IF(W1_lon.gt.W2_lon)THEN
!        RARAD1= W1_lon*DTR     
!        RARAD2= W2_lon*DTR
!        DCRAD1= W1_lat*DTR     
!        DCRAD2= W2_lat*DTR     
!     ELSE                    !(W1_lon.eq.W2_lon)
!        IF(W1_lat.gt.W2_lat)THEN
!           RARAD1= W1_lon*DTR     
!           RARAD2= W2_lon*DTR
!           DCRAD1= W1_lat*DTR     
!           DCRAD2= W2_lat*DTR     
!        ELSE IF(W1_lat.lt.W2_lat)THEN
!           RARAD1= W2_lon*DTR 
!           RARAD2= W1_lon*DTR
!           DCRAD1= W2_lat*DTR     
!           DCRAD2= W1_lat*DTR 
!        ELSE                 !(W1_lat.eq.W2_lat)
!           write(*,*)'Waypoints Error: They must be identical waypoints!'
!        ENDIF
!     ENDIF

     !Distance calcluations between W1 and W2
!     SELECT CASE (ISWITCH)  
!     CASE(0)
!        DELRA  = RARAD1-RARAD2
!        COSIS  = ACOS(SIN(DCRAD1)*SIN(DCRAD2) + COS(DCRAD1)*COS(DCRAD2)*&
!                      COS(DELRA))               !in [rad]
!     CASE(1) 
!        DELDC2 = (DCRAD1-DCRAD2)/2.0_dp            !delta lat
!        DELRA2 = (RARAD1-RARAD2)/2.0_dp            !delta lon
!        SINDIS = SQRT(SIN(DELDC2)*SIN(DELDC2) + COS(DCRAD1)*COS(DCRAD2)*&
!                      SIN(DELRA2)*SIN(DELRA2))
!        COSIS  = 2.0_dp*DASIN(SINDIS)              !in [rad]
!     CASE(2)
!        DELRA2 = RARAD1-RARAD2                  !delta lon
!        SINDIS1= SQRT((COS(DCRAD2)*SIN(DELRA2))**2 + (COS(DCRAD1)*SIN(DCRAD2)-&
!                       SIN(DCRAD1)*COS(DCRAD2)*COS(DELRA2))**2)
!        SINDIS2= SIN(DCRAD1)*SIN(DCRAD2) + COS(DCRAD1)*COS(DCRAD2)*COS(DELRA2) 
!        COSIS  = DATAN2(SINDIS1,SINDIS2)        !in [rad]
!     CASE DEFAULT
!        write(*,*) 'ISWITCH error: Unknown ISWITCH on WP distance calculation!'
!     END SELECT

     !Flight distance calculation between W1 and W2
!!HY20140522-1     DIS     = (radius_earth*0.001_dp + ALT/1000.0_dp)*COSIS      !WP distance = radius * radians [km] 
!!    DIS     = (RADIUS + ALT/1000.0d0)*COSIS      !WP distance = radius * radians [km] 
!!    DISN    = DIS*0.539957d0                     !Unit conversion DIS into nautical miles [nm]
!     FL_DIST = FL_DIST + DIS                    !Total flight distance from D_city to A_city [km]

!!    write(*,*)'Distance |W',j1,'-W',j1+1,'| =',DIS,'[km]',DISN,'[km]',' Total flight distance =',FL_DIST,'[km]'
!!    write(*,*)'ISWITCH =',ISWITCH,'FL_DIR =',FL_DIR


     !!Flight time calculation between W1 and W2
     !!If we re-calculate aircraft speed at every waypoint, AC_V should be calculated here. 
     !!AC_V   = AC_MACH*SOUND_V/1000.0d0                        !Aircraft speed [km/s] 
     !!in [sec]
     !DELT_T = DIS/AC_V                                     !Flight time between W1 and W2 [s]
     !FL_T(j1+1)  = FL_T(j1) + DELT_T                       !Aircraft passing time at each waypoint [s]
     !if(j1.eq.DIV_WAY)T_FL_T = FL_T(j1+1)                  !Total flight time from D_city to A_city [s] 
     !!in [Julian date]
     !DELT_T_Julian = DELT_T/OneDay                         !FLight time between W1 and W2 [Julian date]
     !FL_T_Julian(j1+1) = FL_T_Julian(j1) + DELT_T_Julian   !Aircraft passing time at each waypoint [Julian date]
     !if(j1.eq.DIV_WAY)T_FL_T_Julian = FL_T_Julian(j1+1)    !Arrival time at A_city [Julian date] 


     !!Initialization for next waypoints calculation
     !DIS             = 0.0_dp
     !DELT_T          = 0.0_dp   
     !DELT_T_Julian   = 0.0_dp   
  !ENDDO

!  END SUBROUTINE gc_ac_traj
!----------------------------------------------------------------------

!----------------------------------------------------------------------
  SUBROUTINE determine_domains(THETA_OUT, PHI_OUT, dv1, dv2, dv3, dv4, dv5, dv6)
!----------------------------------------------------------------------
!
! This subroutine calculate waypoints along the great circle.
!
  IMPLICIT NONE
!  REAL(DP)                           :: THETA(DIV_WAY+1), PHI(DIV_WAY+1) 
  REAL(DP), INTENT(IN)               :: THETA_OUT(DIV_WAY+1), PHI_OUT(DIV_WAY+1) 
  REAL(DP), INTENT(OUT)              :: dv1, dv2, dv3, dv4, dv5, dv6 
  REAL(DP)                           :: ratio_btw_waypoints, DIV_DV 
  INTEGER                            :: i, lower_waypoint, upper_waypoint

  !Initializations
  DIV_DV = 0.0_dp
  dv1    = 0.0_dp
  dv2    = 0.0_dp
  dv3    = 0.0_dp
  dv4    = 0.0_dp
  dv5    = 0.0_dp
  dv6    = 0.0_dp
! write(*,*)'check FL_DIR=',FL_DIR

  !Definitions of the centers of dv1, dv3, dv5
  !Definitions of delta_x, delta_y  HY20140428-1
  SELECT CASE(FL_DIR)  
  CASE(0,3)
     DIV_DV = (PHI_OUT(DIV_WAY+1)-PHI_OUT(1))/(NUM_CP_LONLAT+1)
     dv1    =  PHI_OUT(1)+ DIV_DV
     dv3    =  PHI_OUT(1)+(DIV_DV*2.0_dp)
     dv5    =  PHI_OUT(1)+(DIV_DV*3.0_dp)
     delta_x = (PHI_OUT(DIV_WAY+1)-PHI_OUT(1))*percent_x_region*0.5_dp   !HY20140428-1,20141211-1
     delta_y = (PHI_OUT(DIV_WAY+1)-PHI_OUT(1))*percent_y_region*0.5_dp
  CASE(1,2)
     DIV_DV = (PHI_OUT(1)-PHI_OUT(DIV_WAY+1))/(NUM_CP_LONLAT+1)
     dv1    =  PHI_OUT(1)- DIV_DV
     dv3    =  PHI_OUT(1)-(DIV_DV*2.0_dp)
     dv5    =  PHI_OUT(1)-(DIV_DV*3.0_dp)
     delta_x = (PHI_OUT(1)-PHI_OUT(DIV_WAY+1))*percent_x_region*0.5_dp   !HY20140428-1
     delta_y = (PHI_OUT(1)-PHI_OUT(DIV_WAY+1))*percent_y_region*0.5_dp
  CASE(4)
     dv1    =  PHI_OUT(1)
     dv3    =  PHI_OUT(1)
     dv5    =  PHI_OUT(1)
     delta_x = (THETA_OUT(1)-THETA_OUT(DIV_WAY+1))*percent_y_region*0.5_dp   !HY20140428-1,2
     delta_y = (THETA_OUT(1)-THETA_OUT(DIV_WAY+1))*percent_x_region*0.5_dp
  CASE(5)
     dv1    =  PHI_OUT(1)
     dv3    =  PHI_OUT(1)
     dv5    =  PHI_OUT(1)
     delta_x = (THETA_OUT(DIV_WAY+1)-THETA_OUT(1))*percent_y_region*0.5_dp   !HY20140428-1,2
     delta_y = (THETA_OUT(DIV_WAY+1)-THETA_OUT(1))*percent_x_region*0.5_dp
  CASE DEFAULT 
     write(*,*) 'Dv_calc error1: this aircraft is not assinged of any FL_DIR!'
  END SELECT
! write(*,*)'FL_DIR=',FL_DIR,'dv1=',dv1,'dv3=',dv3,'dv5=',dv5


  !Definitions of the centers of dv2, dv4, dv6
  SELECT CASE(FL_DIR)
  CASE(0,3) 
     DO i=1, DIV_WAY
        IF(PHI_OUT(i)<=dv1 .and. PHI_OUT(i+1)>dv1)THEN  
           lower_waypoint = i
           upper_waypoint = i+1
           ratio_btw_waypoints = (PHI_OUT(upper_waypoint)-dv1)/   &
                                 (PHI_OUT(upper_waypoint)-PHI_OUT(lower_waypoint))
           dv2=ratio_btw_waypoints*THETA_OUT(lower_waypoint) + (1.0_dp-ratio_btw_waypoints)*THETA_OUT(upper_waypoint)
           exit  ! leave do loop
        ENDIF
     ENDDO
     DO i=1, DIV_WAY
        IF(PHI_OUT(i)<=dv3 .and. PHI_OUT(i+1)>dv3)THEN
           lower_waypoint = i
           upper_waypoint = i+1
           ratio_btw_waypoints = (PHI_OUT(upper_waypoint)-dv3)/   &
                                 (PHI_OUT(upper_waypoint)-PHI_OUT(lower_waypoint))
           dv4=ratio_btw_waypoints*THETA_OUT(lower_waypoint) + (1.0_dp-ratio_btw_waypoints)*THETA_OUT(upper_waypoint)
           exit  ! leave do loop
        ENDIF
     ENDDO
     DO i=1, DIV_WAY
        IF(PHI_OUT(i)<=dv5 .and. PHI_OUT(i+1)>dv5)THEN
           lower_waypoint = i
           upper_waypoint = i+1
           ratio_btw_waypoints = (PHI_OUT(upper_waypoint)-dv5)/   &
                                 (PHI_OUT(upper_waypoint)-PHI_OUT(lower_waypoint))
           dv6=ratio_btw_waypoints*THETA_OUT(lower_waypoint) + (1.0_dp-ratio_btw_waypoints)*THETA_OUT(upper_waypoint)
           exit  ! leave do loop
        ENDIF
     ENDDO

  CASE(1,2)
     DO i=1, DIV_WAY
       IF(PHI_OUT(i)>=dv1 .and. PHI_OUT(i+1)<dv1)THEN   !see opt20121023-2,20130130-5, MUC to JFK
           lower_waypoint = i
           upper_waypoint = i+1
           ratio_btw_waypoints = (dv1-PHI_OUT(lower_waypoint))/   &
                                 (PHI_OUT(upper_waypoint)-PHI_OUT(lower_waypoint))
           dv2=ratio_btw_waypoints*THETA_OUT(upper_waypoint) + (1.0_dp-ratio_btw_waypoints)*THETA_OUT(lower_waypoint)
           exit  ! leave do loop
        ENDIF
     ENDDO
     DO i=1, DIV_WAY
       IF(PHI_OUT(i)>=dv3 .and. PHI_OUT(i+1)<dv3)THEN
           lower_waypoint = i
           upper_waypoint = i+1
           ratio_btw_waypoints = (dv3-PHI_OUT(lower_waypoint))/   &
                                 (PHI_OUT(upper_waypoint)-PHI_OUT(lower_waypoint))
           dv4=ratio_btw_waypoints*THETA_OUT(upper_waypoint) + (1.0_dp-ratio_btw_waypoints)*THETA_OUT(lower_waypoint)
           exit  ! leave do loop
        ENDIF
     ENDDO
     DO i=1, DIV_WAY
       IF(PHI_OUT(i)>=dv5 .and. PHI_OUT(i+1)<dv5)THEN
           lower_waypoint = i
           upper_waypoint = i+1
           ratio_btw_waypoints = (dv5-PHI_OUT(lower_waypoint))/   &
                                 (PHI_OUT(upper_waypoint)-PHI_OUT(lower_waypoint))
           dv6=ratio_btw_waypoints*THETA_OUT(upper_waypoint) + (1.0_dp-ratio_btw_waypoints)*THETA_OUT(lower_waypoint)
           exit  ! leave do loop
        ENDIF
     ENDDO

  CASE(4)
     DIV_DV = (THETA_OUT(1)-THETA_OUT(DIV_WAY+1))/(NUM_CP_LONLAT+1)
     dv2    =  THETA_OUT(1)- DIV_DV
     dv4    =  THETA_OUT(1)-(DIV_DV*2.0_dp)
     dv6    =  THETA_OUT(1)-(DIV_DV*3.0_dp)

  CASE(5)
     DIV_DV = (THETA_OUT(DIV_WAY+1)-THETA_OUT(1))/(NUM_CP_LONLAT+1)
     dv2    =  THETA_OUT(1)+ DIV_DV
     dv4    =  THETA_OUT(1)+(DIV_DV*2.0_dp)
     dv6    =  THETA_OUT(1)+(DIV_DV*3.0_dp)

  CASE DEFAULT 
      write(*,*) 'Dv_calc error2: this aircraft is not assinged of any FL_DIR!'
  END SELECT

  END SUBROUTINE determine_domains 
!----------------------------------------------------------------------

END MODULE messy_airtraf_tools_ga_armogaset
!HY20140401-2 END PROGRAM messy_airtraf_gc

