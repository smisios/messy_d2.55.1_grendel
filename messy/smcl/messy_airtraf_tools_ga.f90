!-----------------------------------------------
!Reserved from original mainf.f
!-----------------------------------------------
!      program ARMOGAMAIN      
!      integer inistat(2)
!
!      inistat(1)=0
!      call armoga(inistat)
!
!      stop
!      end
!
!-----------------------------------------------
!Reserved from messy_airtraf_tools_ga_mainf.f90, 20140417-6
!-----------------------------------------------
! program ARMOGAMAIN
!      USE messy_airtraf_tools_ga,    ONLY: armoga
!      IMPLICIT NONE
!!      INTEGER,DIMENSION(:) :: inistat(2)
!      DOUBLE PRECISION, DIMENSION(:,:) :: best_traj(41,3)
!      DOUBLE PRECISION :: best_time
!      INTEGER :: i,best_gen
!
!!      inistat(1)=0
!!      call armoga(inistat, best_traj, best_time, best_gen)
!      call armoga(best_traj, best_time, best_gen)
!      write(*,*)"Final opt solution"
!      write(*,*)"Final generation =",best_gen
!      write(*,*)"Final flight time=",best_time
!      DO i=1,41
!         write(*,*)best_traj(i,1),best_traj(i,2),best_traj(i,3)
!      ENDDO
!!      stop
! end program ARMOGAMAIN
!
!-----------------------------------------------
!PROGRAM messy_airtraf_tools_ga
MODULE messy_airtraf_tools_ga
!      USE messy_airtraf_tools_ga_parameter
 USE messy_main_constants_mem, ONLY: DP   !HY20140422-2

IMPLICIT NONE

! Global variables and arrays:

INTEGER :: ga_option_traj_calc_g    !20160302-1,20160304-1
INTEGER :: nlon_ga, nlev_ga, ngl_ga !20170203-1

PUBLIC :: armoga

CONTAINS
!-----------------------------------------------
      !opt_traj >>> ac_routes.  Arguments of subroutine should be corrected!!!!!
      subroutine armoga(ga_option_traj_calc, ga_nwaypoints, ga_ac_routes, ga_p2_ac_routes_desc, ga_philon, ga_philat, &
                        ga_zgl_geopot_3d, ga_uwind_g, ga_vwind_g, ga_v_z_g,                      &
                        ga_t_scb_g, ga_cpc_g, opt_fl_time, opt_gen, opt_vtas_along_traj,         & !20140420-2,opt_traj, 20140524-6
                        opt_fl_direction, opt_vtas, opt_vground, opt_rho, opt_potcov, & !20141216-2,20150312-2,20160301-2,20170605-1  
                        opt_cpc, ga_rho_air_dry_3d_g, ga_press_3d_g,&
                        ga_ATR20O3_g,ga_ATR20CH4_g,ga_ATR20H2O_g,ga_ATR20CPC_g,ga_ATR20CO2_g) ! 20170801  
!HY20140416-4      subroutine armoga(inistat, opt_traj, opt_fl_time, opt_gen)
!Yin_20170423 contrail potential coverage is added
!-----------------------------------------------
!CMAIN  Program EA for MO problems
!   programmed by Daisuke Sasaki
!      Since January 13, 2002
!      Revised April 25, 2003

!  ##### MAIN ROUTINE #####
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!     USE messy_airtraf_gc_w_alt,   ONLY:temp_best_traj,       &  !20140417-6
      USE messy_airtraf_wind,       ONLY:temp_best_traj,       &
                                         temp_best_fl_time,    &
                                         temp_best_fuel_use,   &  !20160606-5
                                         temp_best_nox_emis,   &  !20170216-2 
                                         temp_best_h2o_emis,   &  !20170222-4
                                         temp_best_potcov,     &  !Yin_20170423
                                         temp_best_cost,       &  !20170224-4
                                         temp_best_coc,        &  !20170316-4
                                         temp_best_gen,        &  !20140416-1
                                         temp_vtas_along_traj, &  !20140524-5
                                         w_fl_direction,       &  !20141216-2
                                         temp_best_vtas,       &  !20150312-2
                                         temp_best_vground,    &    !20150312-2
                                         temp_best_rho,        &    !20160623-2 
                                         temp_best_cpc,        &  !20170608-2
                                         temp_best_atr20,      &  !20170801
                                         temp_best_costcpc,    &  !20170801
                                         temp_best_costclim,   &  !20170801
                                         obj_dummy                  !20170221-1
     USE messy_airtraf_gc,          ONLY:props_lon,props_lat,props_alt,props_time,  &  
                                         props_ac_speed,props_dist, & !HY20140420-3  !DIV_WAY 20140721-3
                                         props_potcov                 !Yin_20170423

  IMPLICIT NONE
!HY20140422-3  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: opt_traj(41,6)  !HY20140416-5
  INTEGER, INTENT(INOUT)                        :: ga_nwaypoints !20140721-3,6
  REAL(DP)                                      :: opt_traj(ga_nwaypoints,10)  !HY20140416-5,20140422-3,20140721-3,20170607-2
  INTEGER                                       :: inistat(2)
  REAL(DP), INTENT(OUT)                         :: opt_fl_time
  REAL(DP)                                      :: opt_fuel_use       !20160606-5,20170213-1
  REAL(DP)                                      :: opt_nox_emis       !20170216-3
  REAL(DP)                                      :: opt_h2o_emis       !20170222-5
  !Yin_20170423
  REAL(DP), INTENT(OUT)                         :: opt_cpc            ![km] opt_obj, 20170605-1
  REAL(DP)                                      :: opt_atr20          ![K] opt_obj, 20170801
  REAL(DP)                                      :: opt_costclim       ![K] opt_obj, 20170801
  REAL(DP)                                      :: opt_costcpc        ![K] opt_obj, 20170801
  !Yin_20170423
  REAL(DP)                                      :: opt_cost           !20170224-4
  REAL(DP)                                      :: opt_coc            !20170316-4
  INTEGER, INTENT(OUT)                          :: opt_gen 
  INTEGER, INTENT(OUT)                          :: opt_fl_direction   !20141216-2
  REAL(DP),               INTENT(OUT)           :: opt_vtas_along_traj(ga_nwaypoints)  ![m/s] 20140524-5,20140721-3 
  INTEGER                                       :: iccset,ilnum,i_wp, ga_option_traj_calc   !20160302-1
  LOGICAL                                       :: condition
! INTEGER inistat(2)

! Arguments from messy_airtraf.f90, CALL armoga:
  REAL(DP), DIMENSION(:,:), INTENT(OUT)            :: ga_ac_routes
  REAL(DP), DIMENSION(:),   INTENT(INOUT)          :: ga_p2_ac_routes_desc
  REAL(DP), DIMENSION(:),   INTENT(INOUT)          :: ga_philon   
  REAL(DP), DIMENSION(:),   INTENT(INOUT)          :: ga_philat
  REAL(DP), DIMENSION(:,:,:), INTENT(INOUT)        :: ga_zgl_geopot_3d   !global field 
  REAL(DP), DIMENSION(:,:,:), INTENT(INOUT)        :: ga_uwind_g         !global field 
  REAL(DP), DIMENSION(:,:,:), INTENT(INOUT)        :: ga_vwind_g         !global field
  REAL(DP), DIMENSION(:,:,:), INTENT(INOUT)        :: ga_v_z_g           !global field 
  REAL(DP), DIMENSION(:,:,:), INTENT(INOUT)        :: ga_t_scb_g         !global field  20140522-3
  !Yin_20170423
  REAL(DP), DIMENSION(:,:,:), INTENT(INOUT)        :: ga_cpc_g           !global field  20161207
  REAL(DP), DIMENSION(:,:,:), INTENT(IN),OPTIONAL  :: ga_ATR20O3_g       !global field  20170801
  REAL(DP), DIMENSION(:,:,:), INTENT(IN),OPTIONAL  :: ga_ATR20CH4_g      !global field  20170801
  REAL(DP), DIMENSION(:,:,:), INTENT(IN),OPTIONAL  :: ga_ATR20H2O_g      !global field  20170801
  REAL(DP), DIMENSION(:,:,:), INTENT(IN),OPTIONAL  :: ga_ATR20CPC_g      !global field  20170801
  REAL(DP), DIMENSION(:,:,:), INTENT(IN),OPTIONAL  :: ga_ATR20CO2_g      !global field  20170801
  !Yin_20170423
  REAL(DP), DIMENSION(:,:,:), INTENT(IN),OPTIONAL  :: ga_rho_air_dry_3d_g  !global field  20160301-2,20170203-1 
  REAL(DP), DIMENSION(:,:,:), INTENT(IN),OPTIONAL  :: ga_press_3d_g      !global field  20170215-3 
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: temp_ga_rho_air_dry_3d_g  !global field  20170203-1 
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: temp_ga_press_3d_g  !global field  20170215-3 
!Yin_20170801+
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: temp_ga_ATR20O3_3d_g  !global field  20170801 
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: temp_ga_ATR20CH4_3d_g  !global field  20170801 
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: temp_ga_ATR20H2O_3d_g  !global field  20170801 
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: temp_ga_ATR20CPC_3d_g  !global field  20170801 
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: temp_ga_ATR20CO2_3d_g  !global field  20170801 
!Yin_20170801-
!Outpur Vground/Vtas for 3 city-pairs,20150312-2
  REAL(DP), INTENT(OUT)         :: opt_vtas(ga_nwaypoints,101)     ![m/s] 20150312-2 
  REAL(DP), INTENT(OUT)         :: opt_vground(ga_nwaypoints,101)  ![m/s] 20150312-2 
  REAL(DP), INTENT(OUT)         :: opt_rho(ga_nwaypoints,101)      ![kg/m^3] 20160623-2 
  REAL(DP), INTENT(OUT)         :: opt_potcov(ga_nwaypoints,101)   ![-] 20160608-2 
      
!Add dummy for objective functions, HY20140416-1,20170221-1   
      opt_fl_time  = obj_dummy !=1.0d8 Initial value of total flight time in [s].
      opt_fuel_use = obj_dummy !=1.0d8 Initial value of total fuel use in [kg].20160606-5
      opt_nox_emis = obj_dummy !=1.0d8 Initial value of total nox emission in [g(NO2)].20170216-3
      opt_h2o_emis = obj_dummy !=1.0d8 Initial value of total h2o emission in [g(H2O)].20170222-5
      !Yin_20170423
      opt_cpc      = obj_dummy !1.0d8  Initial value of total potcov [km]. opt_obj, 20170605-1
      opt_atr20    = obj_dummy !1.0d8  Initial value of total ATR20 [K]. opt_obj, 20170801
      opt_costclim = obj_dummy !1.0d8  Initial value of total ATR20 [K]. opt_obj, 20170801
      opt_costcpc  = obj_dummy !1.0d8  Initial value of total ATR20 [K]. opt_obj, 20170801
      !Yin_20170423
      opt_cost     = obj_dummy !=1.0d8 Initial value of total cost in [US Dollar].20170224-4
      opt_coc      = obj_dummy !=1.0d8 Initial value of total cost in [US Dollar].20170316-4
!      write(*,*)'check:obj_dummy in GA',obj_dummy
      ga_option_traj_calc_g = ga_option_traj_calc   !20160304-1
      IF(PRESENT(ga_rho_air_dry_3d_g))THEN         !for fuel_opt, NOx_opt, H2O_opt, Cost_opt, COC_opt
         nlon_ga=SIZE(ga_rho_air_dry_3d_g,dim=1)   !lon,20170203-1
         nlev_ga=SIZE(ga_rho_air_dry_3d_g,dim=2)   !lev,20170203-1
         ngl_ga =SIZE(ga_rho_air_dry_3d_g,dim=3)   !lat,20170203-1         
         ALLOCATE(temp_ga_rho_air_dry_3d_g(nlon_ga, nlev_ga, ngl_ga))
         temp_ga_rho_air_dry_3d_g=ga_rho_air_dry_3d_g   !20170203-1
!         write(*,*)"allocation1:nlon_ga,nlev_ga,ngl_ga",nlon_ga,nlev_ga,ngl_ga
      ENDIF
      IF(PRESENT(ga_press_3d_g))THEN               !for NOx_opt 20170215-3
         ALLOCATE(temp_ga_press_3d_g(nlon_ga, nlev_ga, ngl_ga))
         temp_ga_press_3d_g=ga_press_3d_g 
      ENDIF
      IF(PRESENT(ga_ATR20O3_g))THEN               !for ATR20_opt 20170801
         ALLOCATE(temp_ga_ATR20O3_3d_g(nlon_ga, nlev_ga, ngl_ga))
         temp_ga_ATR20O3_3d_g=ga_ATR20O3_g 
      ENDIF
      IF(PRESENT(ga_ATR20CH4_g))THEN               !for ATR20_opt 20170801
         ALLOCATE(temp_ga_ATR20CH4_3d_g(nlon_ga, nlev_ga, ngl_ga))
         temp_ga_ATR20CH4_3d_g=ga_ATR20CH4_g 
      ENDIF
      IF(PRESENT(ga_ATR20H2O_g))THEN               !for ATR20_opt 20170801
         ALLOCATE(temp_ga_ATR20H2O_3d_g(nlon_ga, nlev_ga, ngl_ga))
         temp_ga_ATR20H2O_3d_g=ga_ATR20H2O_g 
      ENDIF
      IF(PRESENT(ga_ATR20CPC_g))THEN               !for ATR20_opt 20170801
         ALLOCATE(temp_ga_ATR20CPC_3d_g(nlon_ga, nlev_ga, ngl_ga))
         temp_ga_ATR20CPC_3d_g=ga_ATR20CPC_g 
      ENDIF
      IF(PRESENT(ga_ATR20CO2_g))THEN               !for ATR20_opt 20170801
         ALLOCATE(temp_ga_ATR20CO2_3d_g(nlon_ga, nlev_ga, ngl_ga))
         temp_ga_ATR20CO2_3d_g=ga_ATR20CO2_g 
      ENDIF

      inistat(1)  =0    !HY20140416-4  
      iccset      =inistat(1)      
!CDS      write(mtout,*)'STATUS',iccset
      write(mtout,*)'##### ARMOGA Ver ',version,' 1.01'
!      write(*,*)'check_ga_option_traj_calc_1', ga_option_traj_calc_g, ga_option_traj_calc         !20160304-1
      call Initialisation(condition, iccset, ga_p2_ac_routes_desc) !HY20140418-2
 
!cf90      do while(condition.eq..true.)
      do while(condition.eqv..true.)
!C1.0.1        write(mtout,*)'Generation:',igen
!C1.2.0        write(mtout,*)'Generation:',igen+1
        write(mtout,*)'Generation:',iter

        SELECT CASE(ga_option_traj_calc_g)   !20160302-2
        CASE(1,5) !Wind_opt, CPC_opt Yin_20170423, HY20170608-4
           call Evaluation(ga_nwaypoints, ga_p2_ac_routes_desc, ga_philon, ga_philat, ga_zgl_geopot_3d, &
                           ga_uwind_g, ga_vwind_g, ga_v_z_g, ga_t_scb_g, ga_cpc_g)   !HY20140418-3, 20140522-3,20140721-6
        CASE(2,4,6,7,10) !Fuel_opt, H2O_opt, Cost_opt or COC_opt   !20170222-2,20170224-1,20170316-2
           call Evaluation(ga_nwaypoints, ga_p2_ac_routes_desc, ga_philon, ga_philat, ga_zgl_geopot_3d, &
                           ga_uwind_g, ga_vwind_g, ga_v_z_g, ga_t_scb_g, ga_cpc_g, temp_ga_rho_air_dry_3d_g) 
                !HY20140418-3,20140522-3,20140721-6
                                                                                              !20160302-2,20170203-1
!           write(*,*)'check_ga_option_traj_calc_2', ga_option_traj_calc_g, ga_option_traj_calc  !20160304-1
        CASE(3) !NOx_opt
           call Evaluation(ga_nwaypoints, ga_p2_ac_routes_desc, ga_philon, ga_philat, ga_zgl_geopot_3d, &
                           ga_uwind_g, ga_vwind_g, ga_v_z_g, ga_t_scb_g, ga_cpc_g, temp_ga_rho_air_dry_3d_g, &
                           temp_ga_press_3d_g)                      !HY20140418-3,20140522-3,20140721-6
                                                              !20160302-2,20170203-1,20170215-3
!           write(*,*)'check_ga_option_traj_calc_2', ga_option_traj_calc_g, ga_option_traj_calc  !20160304-1
        !Yin_20170801
        CASE(8,9)  !ATR20_opt
!           write(*,*)'tscb_in_ga01:',shape(ga_t_scb_g),ga_t_scb_g(:,:,:)
!           write(*,*)'CPC_in_ga01:',shape(ga_cpc_g),ga_cpc_g(:,:,:)
           call Evaluation(ga_nwaypoints, ga_p2_ac_routes_desc, ga_philon, ga_philat, ga_zgl_geopot_3d, &
                           ga_uwind_g, ga_vwind_g, ga_v_z_g, ga_t_scb_g, ga_cpc_g,temp_ga_rho_air_dry_3d_g,&
                           temp_ga_press_3d_g,temp_ga_ATR20O3_3d_g,temp_ga_ATR20CH4_3d_g,&
                           temp_ga_ATR20H2O_3d_g,temp_ga_ATR20CPC_3d_g,temp_ga_ATR20CO2_3d_g)          
        !Yin_20170801-
        CASE DEFAULT
           write(*,*)'ERROR:ga_option_traj_calc on Evaluation',ga_option_traj_calc_g
        ENDSELECT

        if(ncon.gt.0)call Constraint
        do ilnum=1,nisland
          call UpdateData(ilnum)
          if(istat(ilnum,1).eq.1)then
            call RangeAdapt(ilnum)
          else
            call AssignFitness(ilnum)
            call Reproduction(ilnum)
            call NextPop(ilnum)
          endif
        enddo
        call Judgement(condition)

        !HY20140416-2  Reserve optimal solution  
        ! DO i=1,DIV_WAY+1
        !    temp_best_traj(i,props_lon)=XXN(i)                  !Lon
        !    temp_best_traj(i,props_lat)=YYN(i)                  !Lat
        !    temp_best_traj(i,props_alt)=ZZN(i)                  !Alt [m]
        !    temp_best_traj(i,props_time)=FLIGHT_TIME_JULIAN(i)  !Passing time at each waypoint, [Julian]
        !    temp_best_traj(i,props_ac_speed)=FLIGHT_SPEED(i)*3600.0_dp    !Aircraft_speed, [km/h] !should be FLIGHT_SPEED array!!!!
        !    temp_best_traj(i,props_dist)=FLIGHT_DISTANCE(i)     !Distance for each segment, [km]
        !ENDDO
        SELECT CASE(ga_option_traj_calc_g)   !20160606-5        
        ! Regargin values of temp_best_vtas, temp_best_vground, temp_best_rho for each routing option, 
        ! they are listed in 20170216-4
        ! NOTE: opt_fl_time is passed to messy_airtraf.f90; opt_fuel_use and opt_nox_emis are not passed to messy_airtraf.f90.
        ! Because fuel and NOx calculations are performed again to the optimal traj. in messy_airtraf.f90. 20170216-3
        CASE(1)   !Time_opt=1 
           if(temp_best_fl_time.lt.opt_fl_time)then
              opt_fl_time=temp_best_fl_time               !Total flight time, [s]  !HY20140421-2
              opt_traj   =temp_best_traj                  !temp_best_traj(waypoints,props_time) is Passing time at each waypoint, [Julian]
                                                          !temp_best_traj(waypoints,props_ac_speed) is AC ground speed in [km/h]
                                                          !temp_best_traj(waypoints,props_potcov) is CPC in [km] 20170607-1
              opt_gen    =temp_best_gen
              opt_vtas_along_traj  =temp_vtas_along_traj  !vtas along temp_best traj [m/s] !20140524-5
              opt_vtas   =temp_best_vtas                  !vtas in a plane along waypoints,20150312-2
              opt_vground=temp_best_vground               !vground in a plane along waypoints,20150312-2 
              opt_rho    =temp_best_rho                   !rho in a plane along waypoints,20160623-2,20170214-1 
              opt_potcov          = temp_best_cpc         !Potcov in a plane along waypoints [-]20170606-4,20170608-2
              !Yin_20170423
              opt_cpc    =temp_best_potcov                !Total CPC [km], 20170605-1,2,20170607-3
              !opt_cpc_along_traj   =temp_cpc_along_traj  !CPC along temp_best traj [km] !20170605-2,20170607-2              
              !Yin_20170423
           endif
        CASE(2)   !Fuel_opt  20160606-5
           if(temp_best_fuel_use.lt.opt_fuel_use)then
              opt_fuel_use=temp_best_fuel_use             !Total fuel use, [kg]  !20160606-5
              opt_fl_time=temp_best_fl_time               !Total flight time, [s]  !HY20140421-2
              opt_traj   =temp_best_traj                  !temp_best_traj(waypoints,props_time) is Passing time at each waypoint, [Julian]
                                                          !temp_best_traj(waypoints,props_ac_speed) is AC ground speed in [km/h]
                                                          !temp_best_traj(waypoints,props_potcov) is CPC in [km] 20170607-1
              opt_gen    =temp_best_gen
              opt_vtas_along_traj  =temp_vtas_along_traj  !vtas along temp_best traj [m/s] !20140524-5
              opt_vtas   =temp_best_vtas                  !vtas in a plane along waypoints,20150312-2
              opt_vground=temp_best_vground               !vground in a plane along waypoints,20150312-2 
              opt_rho    =temp_best_rho                   !rho in a plane along waypoints,20160623-2,20140214-1 
              opt_potcov          = temp_best_cpc         !Potcov in a plane along waypoints [-]20170606-4,20170608-2
              !Yin_20170423
              opt_cpc    =temp_best_potcov                !Total CPC [km], 20170605-1,2,20170607-3
              !opt_cpc_along_traj   =temp_cpc_along_traj  !CPC along temp_best traj [km] !20170605-2,20170607-2              
              !Yin_20170423
           endif
        CASE(3)   !NOx_opt  20170216-3
           if(temp_best_nox_emis.lt.opt_nox_emis)then
              opt_nox_emis=temp_best_nox_emis             !Total NOx emission, [g(NO2)] 20170216-3 
              opt_fl_time=temp_best_fl_time               !Total flight time, [s]  !HY20140421-2
              opt_traj   =temp_best_traj                  !temp_best_traj(waypoints,props_time) is Passing time at each waypoint, [Julian]
                                                          !temp_best_traj(waypoints,props_ac_speed) is AC ground speed in [km/h]
                                                          !temp_best_traj(waypoints,props_potcov) is CPC in [km] 20170607-1
              opt_gen    =temp_best_gen
              opt_vtas_along_traj  =temp_vtas_along_traj  !vtas along temp_best traj [m/s] !20140524-5
              opt_vtas   =temp_best_vtas                  !vtas in a plane along waypoints,20150312-2
              opt_vground=temp_best_vground               !vground in a plane along waypoints,20150312-2 
              opt_rho    =temp_best_rho                   !rho in a plane along waypoints,20160623-2,20140214-1 
              opt_potcov          = temp_best_cpc         !Potcov in a plane along waypoints [-]20170606-4,20170608-2
              !Yin_20170423
              opt_cpc    =temp_best_potcov                !Total CPC [km], 20170605-1,2,20170607-3 
              !opt_cpc_along_traj   =temp_cpc_along_traj  !CPC along temp_best traj [km] !20170605-2,20170607-2              
              !Yin_20170423
           endif
        CASE(4)   !H2O_opt  20170222-4
           if(temp_best_h2o_emis.lt.opt_h2o_emis)then
              opt_h2o_emis=temp_best_h2o_emis             !Total H2O emission, [g(H2O)]  !20170222-4
              opt_fl_time=temp_best_fl_time               !Total flight time, [s]  !HY20140421-2
              opt_traj   =temp_best_traj                  !temp_best_traj(waypoints,props_time) is Passing time at each waypoint, [Julian]
                                                          !temp_best_traj(waypoints,props_ac_speed) is AC ground speed in [km/h]
                                                          !temp_best_traj(waypoints,props_potcov) is CPC in [km] 20170607-1
              opt_gen    =temp_best_gen
              opt_vtas_along_traj  =temp_vtas_along_traj  !vtas along temp_best traj [m/s] !20140524-5
              opt_vtas   =temp_best_vtas                  !vtas in a plane along waypoints,20150312-2
              opt_vground=temp_best_vground               !vground in a plane along waypoints,20150312-2 
              opt_rho    =temp_best_rho                   !rho in a plane along waypoints,20160623-2,20140214-1 
              opt_potcov          = temp_best_cpc         !Potcov in a plane along waypoints [-]20170606-4,20170608-2
              !Yin_20170423
              opt_cpc    =temp_best_potcov                !Total CPC [km], 20170605-1,2,20170607-3 
              !opt_cpc_along_traj   =temp_cpc_along_traj  !CPC along temp_best traj [km] !20170605-2,20170607-2              
              !Yin_20170423
           endif
        !Yin_20170423
        CASE(5)   !ContrailPC_opt
!           write(*,*)'SMCL: option ContrailPC_opt'
           if(temp_best_potcov.lt.opt_cpc)then          !20170605-1,20170608-3
              opt_cpc             = temp_best_potcov    !Total CPC [km], 20170605-1,2,20170607-3
              opt_fl_time         = temp_best_fl_time
              opt_traj            = temp_best_traj
              opt_gen             = temp_best_gen
              opt_vtas_along_traj = temp_vtas_along_traj
              opt_vtas            = temp_best_vtas
              opt_vground         = temp_best_vground
              opt_rho             = temp_best_rho         !rho in a plane along waypoints,20160623-2,20140214-1,20180511-3 
              opt_potcov          = temp_best_cpc         !Potcov in a plane along waypoints [-]20170606-4,20170608-2
              !Yin_20170423
              !opt_cpc_along_traj  =temp_cpc_along_traj_c !CPC along temp_best traj [km] !20170605-2,20470607-2
              !Yin_20170423
           endif          
        !Yin_20170423
        CASE(6)   !Cost_opt  20170224-4
           if(temp_best_cost.lt.opt_cost)then
              opt_cost   =temp_best_cost                  !Total cost, [US Dollar]  !20170224-4
              opt_fl_time=temp_best_fl_time               !Total flight time, [s]  !HY20140421-2
              opt_traj   =temp_best_traj                  !temp_best_traj(waypoints,props_time) is Passing time at each waypoint, [Julian]
                                                          !temp_best_traj(waypoints,props_ac_speed) is AC ground speed in [km/h]
                                                          !temp_best_traj(waypoints,props_potcov) is CPC in [km] 20170607-1
              opt_gen    =temp_best_gen
              opt_vtas_along_traj  =temp_vtas_along_traj  !vtas along temp_best traj [m/s] !20140524-5
              opt_vtas   =temp_best_vtas                  !vtas in a plane along waypoints,20150312-2
              opt_vground=temp_best_vground               !vground in a plane along waypoints,20150312-2 
              opt_rho    =temp_best_rho                   !rho in a plane along waypoints,20160623-2,20140214-1 
              opt_potcov          = temp_best_cpc         !Potcov in a plane along waypoints [-]20170606-4,20170608-2
              !Yin_20170423
              opt_cpc    =temp_best_potcov                !Total CPC [km], 20170605-1,2,20170607-3 
              !opt_cpc_along_traj   =temp_cpc_along_traj  !CPC along temp_best traj [km] !20170605-2,20470607-2              
              !Yin_20170423
           endif
        CASE(7)   !COC_opt  20170316-4
           if(temp_best_coc.lt.opt_coc)then
              opt_coc    =temp_best_coc                   !Total coc, [US Dollar]  !20170316-4
              opt_fl_time=temp_best_fl_time               !Total flight time, [s]  !HY20140421-2
              opt_traj   =temp_best_traj                  !temp_best_traj(waypoints,props_time) is Passing time at each waypoint, [Julian]
                                                          !temp_best_traj(waypoints,props_ac_speed) is AC ground speed in [km/h]
                                                          !temp_best_traj(waypoints,props_potcov) is CPC in [km] 20170607-1
              opt_gen    =temp_best_gen
              opt_vtas_along_traj  =temp_vtas_along_traj  !vtas along temp_best traj [m/s] !20140524-5
              opt_vtas   =temp_best_vtas                  !vtas in a plane along waypoints,20150312-2
              opt_vground=temp_best_vground               !vground in a plane along waypoints,20150312-2 
              opt_rho    =temp_best_rho                   !rho in a plane along waypoints,20160623-2,20140214-1 
              opt_potcov          = temp_best_cpc         !Potcov in a plane along waypoints [-]20170606-4,20170608-2
              !Yin_20170423
              opt_cpc    =temp_best_potcov                !Total CPC [km], 20170605-1,2,20170607-3 
              !opt_cpc_along_traj   =temp_cpc_along_traj  !CPC along temp_best traj [km] !20170605-2,20470607-2              
              !Yin_20170423
           endif
        !Yin_20170801
        CASE(8)   !ATR20_opt
!           write(*,*)'SMCL: option ATR20_opt'
           if(temp_best_atr20.lt.opt_atr20)then
              opt_atr20           = temp_best_atr20
              opt_fl_time         = temp_best_fl_time
              opt_traj            = temp_best_traj
              opt_gen             = temp_best_gen
              opt_vtas_along_traj = temp_vtas_along_traj
              opt_vtas            = temp_best_vtas
              opt_vground         = temp_best_vground
              opt_rho             = temp_best_rho
              opt_potcov          = temp_best_cpc         !20180511-3
              opt_cpc             = temp_best_potcov      !Total CPC [km], 20170605-1,2,20170607-3,20180511-3 
           endif
        CASE(9)   !COSTCLIM_opt
!           write(*,*)'SMCL: option COSTCLIM_opt'
           if(temp_best_costclim.lt.opt_costclim)then
              opt_costclim        = temp_best_costclim
              opt_fl_time         = temp_best_fl_time
              opt_traj            = temp_best_traj
              opt_gen             = temp_best_gen
              opt_vtas_along_traj = temp_vtas_along_traj
              opt_vtas            = temp_best_vtas
              opt_vground         = temp_best_vground
              opt_rho             = temp_best_rho
              opt_potcov          = temp_best_cpc         !20180511-3
              opt_cpc             = temp_best_potcov      !Total CPC [km], 20170605-1,2,20170607-3,20180511-3 
           endif
        CASE(10)   !COSTCPC_opt
!           write(*,*)'SMCL: option COSTCPC_opt'
           if(temp_best_costcpc.lt.opt_costcpc)then
              opt_costcpc         = temp_best_costcpc
              opt_fl_time         = temp_best_fl_time
              opt_traj            = temp_best_traj
              opt_gen             = temp_best_gen
              opt_vtas_along_traj = temp_vtas_along_traj
              opt_vtas            = temp_best_vtas
              opt_vground         = temp_best_vground
              opt_rho             = temp_best_rho
              opt_potcov          = temp_best_cpc         !20180511-3
              opt_cpc             = temp_best_potcov      !Total CPC [km], 20170605-1,2,20170607-3,20180511-3 
           endif

        !Yin_20170423
        CASE DEFAULT
           write(*,*)"Error:Reserve optimal solution"
        ENDSELECT

        !check current optimal solution
        !write(*,*)"Current opt solution" 
        !write(*,*)"Generation =",opt_gen
        !write(*,*)"Flight time=",opt_fl_time
        !write(*,*)"opt_vtas=",opt_vtas_along_traj

        !output: history of convergence, 20140706-2,20141218-3,20150305-3
        !MUC >> JFK
!        IF((INT(ga_p2_ac_routes_desc(1)).eq.11).and.(INT(ga_p2_ac_routes_desc(2)).eq.48).and.  &
!           (INT(ga_p2_ac_routes_desc(4)).eq.-73).and.(INT(ga_p2_ac_routes_desc(5)).eq.40))THEN
!        IF((INT(ga_p2_ac_routes_desc(1)).eq.-83) .and. (INT(ga_p2_ac_routes_desc(2)).eq.42).and.  &
!           (INT(ga_p2_ac_routes_desc(4)).eq.8) .and. (INT(ga_p2_ac_routes_desc(5)).eq.50))THEN
!           write(*,*)"history_muc_jfk",INT(ga_p2_ac_routes_desc(1)),INT(ga_p2_ac_routes_desc(2)),INT(ga_p2_ac_routes_desc(4)), &
!                     INT(ga_p2_ac_routes_desc(5))
!           OPEN(23,file='history_muc_jfk.dat',position='append')
!              if(ga_option_traj_calc_g.eq.1)write(23,*)opt_gen,opt_fl_time    !Time_opt
!              if(ga_option_traj_calc_g.eq.2)write(23,*)opt_gen,opt_fuel_use   !Fuel_opt,20160606-6 
!              if(ga_option_traj_calc_g.eq.3)write(23,*)opt_gen,opt_nox_emis   !NOx_opt,20170216-3 
!              if(ga_option_traj_calc_g.eq.4)write(23,*)opt_gen,opt_h2o_emis   !H2O_opt,20170222-5 
!              if(ga_option_traj_calc_g.eq.5)write(23,*)opt_gen,opt_cpc        !CPC_opt,20170605-1
!              if(ga_option_traj_calc_g.eq.6)write(23,*)opt_gen,opt_cost       !Cost_opt,20170224-4 
!              if(ga_option_traj_calc_g.eq.7)write(23,*)opt_gen,opt_coc        !COC_opt,20170316-5 
!              if(ga_option_traj_calc_g.eq.8)write(23,*)opt_gen,opt_atr20      !ATR_opt,20180511-2 
!           CLOSE(23)
!        ENDIF
        !JFK >> MUC 20141218-3
!        IF((INT(ga_p2_ac_routes_desc(1)).eq.-73).and.(INT(ga_p2_ac_routes_desc(2)).eq.40).and.  &
!           (INT(ga_p2_ac_routes_desc(4)).eq.11).and.(INT(ga_p2_ac_routes_desc(5)).eq.48))THEN
!          write(*,*)"history_jfk_muc" 
!           OPEN(24,file='history_jfk_muc.dat',position='append')
!              if(ga_option_traj_calc_g.eq.1)write(24,*)opt_gen,opt_fl_time    !Time_opt
!              if(ga_option_traj_calc_g.eq.2)write(24,*)opt_gen,opt_fuel_use   !Fuel_opt,20160606-6    
!              if(ga_option_traj_calc_g.eq.3)write(24,*)opt_gen,opt_nox_emis   !NOx_opt,20170216-3   
!              if(ga_option_traj_calc_g.eq.4)write(24,*)opt_gen,opt_h2o_emis   !H2O_opt,20170222-5 
!              if(ga_option_traj_calc_g.eq.5)write(24,*)opt_gen,opt_cpc        !CPC_opt,20170605-1
!              if(ga_option_traj_calc_g.eq.6)write(24,*)opt_gen,opt_cost       !Cost_opt,20170224-4 
!              if(ga_option_traj_calc_g.eq.7)write(24,*)opt_gen,opt_coc        !COC_opt,20170316-5 
!              if(ga_option_traj_calc_g.eq.8)write(24,*)opt_gen,opt_atr20      !ATR_opt,20180511-2 
!           CLOSE(24)
!        ENDIF
        !EHAM >> KSEA
!        IF((INT(ga_p2_ac_routes_desc(1)).eq.4).and.(INT(ga_p2_ac_routes_desc(2)).eq.52).and.  &
!           (INT(ga_p2_ac_routes_desc(4)).eq.-122).and.(INT(ga_p2_ac_routes_desc(5)).eq.47))THEN
!           OPEN(25,file='history_eham_ksea.dat',position='append')
!              if(ga_option_traj_calc_g.eq.1)write(25,*)opt_gen,opt_fl_time    !Time_opt
!              if(ga_option_traj_calc_g.eq.2)write(25,*)opt_gen,opt_fuel_use   !Fuel_opt,20160606-6
!              if(ga_option_traj_calc_g.eq.3)write(25,*)opt_gen,opt_nox_emis   !NOx_opt,20170216-3
!              if(ga_option_traj_calc_g.eq.4)write(25,*)opt_gen,opt_h2o_emis   !H2O_opt,20170222-5 
!              if(ga_option_traj_calc_g.eq.5)write(25,*)opt_gen,opt_cpc        !CPC_opt,20170605-1
!              if(ga_option_traj_calc_g.eq.6)write(25,*)opt_gen,opt_cost       !Cost_opt,20170224-4 
!              if(ga_option_traj_calc_g.eq.7)write(25,*)opt_gen,opt_coc        !COC_opt,20170316-5 
!              if(ga_option_traj_calc_g.eq.8)write(25,*)opt_gen,opt_atr20      !ATR_opt,20180511-2 
!           CLOSE(25)
!        ENDIF
        !KSEA >> EHAM
!        IF((INT(ga_p2_ac_routes_desc(1)).eq.-122).and.(INT(ga_p2_ac_routes_desc(2)).eq.47).and.  &
!           (INT(ga_p2_ac_routes_desc(4)).eq.4).and.(INT(ga_p2_ac_routes_desc(5)).eq.52))THEN
!           OPEN(26,file='history_ksea_eham.dat',position='append')
!              if(ga_option_traj_calc_g.eq.1)write(26,*)opt_gen,opt_fl_time    !Time_opt
!              if(ga_option_traj_calc_g.eq.2)write(26,*)opt_gen,opt_fuel_use   !Fuel_opt,20160606-6   
!              if(ga_option_traj_calc_g.eq.3)write(26,*)opt_gen,opt_nox_emis   !NOx_opt,20170216-3   
!              if(ga_option_traj_calc_g.eq.4)write(26,*)opt_gen,opt_h2o_emis   !H2O_opt,20170222-5 
!              if(ga_option_traj_calc_g.eq.5)write(26,*)opt_gen,opt_cpc        !CPC_opt,20170605-1
!              if(ga_option_traj_calc_g.eq.6)write(26,*)opt_gen,opt_cost       !Cost_opt,20170224-4 
!              if(ga_option_traj_calc_g.eq.7)write(26,*)opt_gen,opt_coc        !COC_opt,20170316-5 
!              if(ga_option_traj_calc_g.eq.8)write(26,*)opt_gen,opt_atr20      !ATR_opt,20180511-2 
!           CLOSE(26)
!        ENDIF
        !EHAM >> KMSP
!        IF((INT(ga_p2_ac_routes_desc(1)).eq.4).and.(INT(ga_p2_ac_routes_desc(2)).eq.52).and.  &
!           (INT(ga_p2_ac_routes_desc(4)).eq.-93).and.(INT(ga_p2_ac_routes_desc(5)).eq.44).and.  &
!           (INT(ga_p2_ac_routes_desc(3)).ge.2443510))THEN
!           OPEN(27,file='history_eham_kmsp.dat',position='append')
!              if(ga_option_traj_calc_g.eq.1)write(27,*)opt_gen,opt_fl_time    !Time_opt
!              if(ga_option_traj_calc_g.eq.2)write(27,*)opt_gen,opt_fuel_use   !Fuel_opt,20160606-6
!              if(ga_option_traj_calc_g.eq.3)write(27,*)opt_gen,opt_nox_emis   !NOx_opt,20170216-3
!              if(ga_option_traj_calc_g.eq.4)write(27,*)opt_gen,opt_h2o_emis   !H2O_opt,20170222-5 
!              if(ga_option_traj_calc_g.eq.5)write(27,*)opt_gen,opt_cpc        !CPC_opt,20170605-1
!              if(ga_option_traj_calc_g.eq.6)write(27,*)opt_gen,opt_cost       !Cost_opt,20170224-4 
!              if(ga_option_traj_calc_g.eq.7)write(27,*)opt_gen,opt_coc        !COC_opt,20170316-5 
!              if(ga_option_traj_calc_g.eq.8)write(27,*)opt_gen,opt_atr20      !ATR_opt,20180511-2 
!           CLOSE(27)
!        ENDIF
        !KMSP >> EHAM
!        IF((INT(ga_p2_ac_routes_desc(1)).eq.-93).and.(INT(ga_p2_ac_routes_desc(2)).eq.44).and.  &
!           (INT(ga_p2_ac_routes_desc(4)).eq.4).and.(INT(ga_p2_ac_routes_desc(5)).eq.52).and.  &
!           (INT(ga_p2_ac_routes_desc(3)).ge.2443510))THEN
!           OPEN(28,file='history_kmsp_eham.dat',position='append')
!              if(ga_option_traj_calc_g.eq.1)write(28,*)opt_gen,opt_fl_time    !Time_opt
!              if(ga_option_traj_calc_g.eq.2)write(28,*)opt_gen,opt_fuel_use   !Fuel_opt,20160606-6 
!              if(ga_option_traj_calc_g.eq.3)write(28,*)opt_gen,opt_nox_emis   !NOx_opt,20170216-3 
!              if(ga_option_traj_calc_g.eq.4)write(28,*)opt_gen,opt_h2o_emis   !H2O_opt,20170222-5 
!              if(ga_option_traj_calc_g.eq.5)write(28,*)opt_gen,opt_cpc        !CPC_opt,20170605-1
!              if(ga_option_traj_calc_g.eq.6)write(28,*)opt_gen,opt_cost       !Cost_opt,20170224-4 
!              if(ga_option_traj_calc_g.eq.7)write(28,*)opt_gen,opt_coc        !COC_opt,20170316-5 
!              if(ga_option_traj_calc_g.eq.8)write(28,*)opt_gen,opt_atr20      !ATR_opt,20180511-2 
!           CLOSE(28)
!        ENDIF

      enddo
      DO i_wp=1,ga_nwaypoints   !HY20140420-2,20140421-2
         ga_ac_routes(i_wp,props_lon)     =opt_traj(i_wp,props_lon)  !XXN(i)                !Lon
         ga_ac_routes(i_wp,props_lat)     =opt_traj(i_wp,props_lat)  !YYN(i)                !Lat
         ga_ac_routes(i_wp,props_alt)     =opt_traj(i_wp,props_alt)  !ZZN(i)                !Alt [m]
         ga_ac_routes(i_wp,props_time)    =opt_traj(i_wp,props_time) !FLIGHT_TIME_JULIAN(i) !Passing time at each waypoint,[Julian]
         ga_ac_routes(i_wp,props_ac_speed)=opt_traj(i_wp,props_ac_speed)  !FLIGHT_SPEED(i)  !AC ground speed, [km/h] 
         ga_ac_routes(i_wp,props_dist)    =opt_traj(i_wp,props_dist)  !FLIGHT_DISTANCE(i)   !Distance for each segment, [km]
         ga_ac_routes(i_wp,props_potcov)  =opt_traj(i_wp,props_potcov) !CPC for each segment, [km]20170605-2,20170606-4,20170607-2
!         write(*,*)'CPCpropcheck_GA:',ga_ac_routes(i_wp,props_potcov)
      ENDDO
      !20141216-2

!HY
!     IF(ga_option_traj_calc_g/=5)THEN   !20170607-3
      opt_fl_direction = w_fl_direction
      !deallocate:20140722-1
      DEALLOCATE(temp_best_traj)
      DEALLOCATE(temp_vtas_along_traj)
      !Yin_20170423
      !DEALLOCATE(temp_cpc_along_traj)   !20170607-2
      !Yin_20170423
      !20150312-2 
      DEALLOCATE(temp_best_vtas)
      DEALLOCATE(temp_best_vground)
      DEALLOCATE(temp_best_rho)   !20160623-2,20170214-1
      DEALLOCATE(temp_best_cpc)   !20180511-3
!     ENDIF   !20170607-3

      !Yin_20170423
!     IF(option_traj_calc==ContrailPC_opt)THEN
!     IF(ga_option_traj_calc_g==ContrailPC_opt)THEN
     ! IF(ga_option_traj_calc_g==5)THEN
     !   opt_fl_direction = w_fl_direction         !contrail_fl_direction
     !   DEALLOCATE(temp_best_traj)
     !   DEALLOCATE(temp_vtas_along_traj)
       !DEALLOCATE(temp_cpc_along_traj_c)         !20170608-3
     !   DEALLOCATE(temp_best_vtas)
     !   DEALLOCATE(temp_best_vground)
!     ELSEIF(option_traj_calc==ECOCLIM)THEN
!     ELSEIF(ga_option_traj_calc_g==ECOCLIM)THEN
!        opt_fl_direction = multiObj_fl_direction
!        DEALLOCATE(temp_best_traj_mo)
!        DEALLOCATE(temp_vtas_along_traj_mo)
!        DEALLOCATE(temp_cpc_along_traj_mo)
!        DEALLOCATE(temp_best_vtas_mo)
!        DEALLOCATE(temp_best_vground_mo)
     ! ENDIF
      !Yin_20170423

      !For fuel_opt, NOx_opt, H2O_opt, Cost_opt, COC_opt.
      IF(ALLOCATED(temp_ga_rho_air_dry_3d_g))DEALLOCATE(temp_ga_rho_air_dry_3d_g)   !20170203-1,20170212-3,20170216-3
      !For NOx_opt. 
      IF(ALLOCATED(temp_ga_press_3d_g))DEALLOCATE(temp_ga_press_3d_g)   !20170215-3
      !Yin_20170801+
      IF(ALLOCATED(temp_ga_ATR20O3_3d_g))DEALLOCATE(temp_ga_ATR20O3_3d_g)   !20170801
      IF(ALLOCATED(temp_ga_ATR20CH4_3d_g))DEALLOCATE(temp_ga_ATR20CH4_3d_g)   !20170801
      IF(ALLOCATED(temp_ga_ATR20H2O_3d_g))DEALLOCATE(temp_ga_ATR20H2O_3d_g)   !20170801
      IF(ALLOCATED(temp_ga_ATR20CPC_3d_g))DEALLOCATE(temp_ga_ATR20CPC_3d_g)   !20170801
      IF(ALLOCATED(temp_ga_ATR20CO2_3d_g))DEALLOCATE(temp_ga_ATR20CO2_3d_g)   !20170801
      !Yin_20170801-
!C1.2.0      if(igen.gt.0)call PostProcess
      if(iter.gt.1)call PostProcess

!CMAIN      stop
      return
      end subroutine armoga

! *****************************************
      subroutine Initialisation(condition, iccset, init_p2_ac_routes_desc)  !HY20140418-2
! *****************************************
!      USE fort_1,      ONLY:testp   !HY20140331-2
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6      USE messy_airtraf_tools_ga_initial,     ONLY: CommonSetting, InputOperator, CheckDesign, &
!HY                                      SetParam, InitialParam, InitialPop, OutputInitial       
!HY      USE messy_airtraf_tools_ga_testfunction,    ONLY: TestDesign
!HY      USE messy_airtraf_tools_ga_testdirect,     ONLY: SetDesign 
!HY      USE messy_airtraf_tools_ga_iofile,     ONLY: InputGen
      IMPLICIT NONE
      REAL(DP), DIMENSION(:),   INTENT(INOUT)  :: init_p2_ac_routes_desc
      INTEGER :: iccset,icheck,ncont,ndvt,nobjt
      LOGICAL :: condition
      condition = .true.
      icheck = 1
   
! --------- INPUT PARAMETERS ---------
      call CommonSetting(nobjt,ncont,ndvt,iccset)

! --------- DESIGN INFORMATION ---------
      if(testp.ne.'OUTE')call TestDesign
      if(isetdes.eq.1)call SetDesign(init_p2_ac_routes_desc)  !HY20140418-2

!      if(testp.eq.'OUTE')then
!        call SetDesign
!      else
!        call TestDesign
!      endif

! --------- INPUT EA OPERATORS ---------
      call InputOperator

! --------- CHECK DESIGN ---------
      call CheckDesign(nobjt,ncont,ndvt)

! --------- SET INITIAL PARAMETER ---------
      call SetParam

!C1.2.0      if(igen.lt.0)then
      if(iter.lt.1)then
        write(mtout,*)'Initial Generation'
        call InitialParam
        call InitialPop
        call OutputInitial
        if(iman.eq.1)then
          write(mtout,*) 'Manual Mode Optimization:'
          write(mtout,*) '  0th Generation has been created'
          condition = .false.
        endif
      else
!C1.0.1        write(mtout,*)'Restart Generation:',igen
!C1.2.0        write(mtout,*)'Restart Generation:',igen+1
        write(mtout,*)'Restart Generation:',iter
!CDS-IMPORTANT        if((iflag(8).gt.0).and.(igen.gt.0))
!C1.2.0        if((iflag(8).gt.0).and.(igen.ge.0))
!HY20140327-2     if((iflag(8).gt.0).and.(iter.ge.1))
!HY20140327-2  &      call InputAllData(mtall,icheck)
        call InputGen(icheck)
!        if(icheck.eq.3)call ArcCheck(0)
      endif

! --------- OUTPUT PARAMETERS ---------
!HY20140327-1   call OutputOperator

      return
      end subroutine Initialisation

! *******************************
      subroutine Evaluation(eval_nwaypoints, eval_p2_ac_routes_desc, eval_philon, eval_philat, eval_zgl_geopot_3d, &
                            eval_uwind_g, eval_vwind_g, eval_v_z_g, eval_t_scb_g, eval_cpc_g, &  !HY20140418-3,20140522-3,20140721-6 
                            eval_rho_air_dry_3d_g, eval_press_3d_g,eval_ATR20O3_3d_g,&
                            eval_ATR20CH4_3d_g,eval_ATR20H2O_3d_g,eval_ATR20CPC_3d_g,eval_ATR20CO2_3d_g) !20170801 
! *******************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6      USE messy_airtraf_tools_ga_testfunction,    ONLY: EvaPop
  IMPLICIT NONE
  INTEGER, INTENT(INOUT)                           :: eval_nwaypoints         !20140721-6
  REAL(DP), DIMENSION(:),   INTENT(INOUT)          :: eval_p2_ac_routes_desc  !HY20140418-3
  REAL(DP), DIMENSION(:),   INTENT(INOUT)          :: eval_philon
  REAL(DP), DIMENSION(:),   INTENT(INOUT)          :: eval_philat
  REAL(DP), DIMENSION(:,:,:), INTENT(INOUT)        :: eval_zgl_geopot_3d   !global field 
  REAL(DP), DIMENSION(:,:,:), INTENT(INOUT)        :: eval_uwind_g         !global field 
  REAL(DP), DIMENSION(:,:,:), INTENT(INOUT)        :: eval_vwind_g         !global field
  REAL(DP), DIMENSION(:,:,:), INTENT(INOUT)        :: eval_v_z_g           !global field 
  REAL(DP), DIMENSION(:,:,:), INTENT(INOUT)        :: eval_t_scb_g         !global field  20140522-3 
  !Yin_20170423
  REAL(DP), DIMENSION(:,:,:), INTENT(INOUT)        :: eval_cpc_g           !global field  20161215 
  !Yin_20170423
  REAL(DP), DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: eval_rho_air_dry_3d_g  !global field  20160302-2,20170206-1 
  REAL(DP), DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: eval_press_3d_g      !global field 20170215-3 
  REAL(DP), DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: eval_ATR20O3_3d_g      !global field 20170801 
  REAL(DP), DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: eval_ATR20CH4_3d_g      !global field 20170801 
  REAL(DP), DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: eval_ATR20H2O_3d_g      !global field 20170801 
  REAL(DP), DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: eval_ATR20CPC_3d_g      !global field 20170801 
  REAL(DP), DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: eval_ATR20CO2_3d_g      !global field 20170801 
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: temp_eval_rho_air_dry_3d_g  !global field  20170206-1 
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: temp_eval_press_3d_g  !global field  20170215-3
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: temp_eval_ATR20O3_3d_g  !global field  20170801
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: temp_eval_ATR20CH4_3d_g  !global field  20170801
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: temp_eval_ATR20H2O_3d_g  !global field  20170801
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: temp_eval_ATR20CPC_3d_g  !global field  20170801
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: temp_eval_ATR20CO2_3d_g  !global field  20170801
  INTEGER                                          :: ist,ied,ilnum,ip,ipop

      IF(PRESENT(eval_rho_air_dry_3d_g))THEN   !for fuel_opt, NOx_opt, H2O_opt, Cost_opt, COC_opt. 20170206-1
         ALLOCATE(temp_eval_rho_air_dry_3d_g(nlon_ga, nlev_ga, ngl_ga))
         temp_eval_rho_air_dry_3d_g=eval_rho_air_dry_3d_g
!         write(*,*)"allocation2:nlon_ga,nlev_ga,ngl_ga",nlon_ga,nlev_ga,ngl_ga
      ENDIF
      IF(PRESENT(eval_press_3d_g))THEN   !for NOx_opt. 20170215-3
         ALLOCATE(temp_eval_press_3d_g(nlon_ga, nlev_ga, ngl_ga))
         temp_eval_press_3d_g=eval_press_3d_g
      ENDIF
      IF(PRESENT(eval_ATR20O3_3d_g))THEN   !for ATR20_opt. 20170801
         ALLOCATE(temp_eval_ATR20O3_3d_g(nlon_ga, nlev_ga, ngl_ga))
         temp_eval_ATR20O3_3d_g=eval_ATR20O3_3d_g
      ENDIF
      IF(PRESENT(eval_ATR20CH4_3d_g))THEN   !for ATR20_opt. 20170801
         ALLOCATE(temp_eval_ATR20CH4_3d_g(nlon_ga, nlev_ga, ngl_ga))
         temp_eval_ATR20CH4_3d_g=eval_ATR20CH4_3d_g
      ENDIF
      IF(PRESENT(eval_ATR20H2O_3d_g))THEN   !for ATR20_opt. 20170801
         ALLOCATE(temp_eval_ATR20H2O_3d_g(nlon_ga, nlev_ga, ngl_ga))
         temp_eval_ATR20H2O_3d_g=eval_ATR20H2O_3d_g
      ENDIF
      IF(PRESENT(eval_ATR20CPC_3d_g))THEN   !for ATR20_opt. 20170801
         ALLOCATE(temp_eval_ATR20CPC_3d_g(nlon_ga, nlev_ga, ngl_ga))
         temp_eval_ATR20CPC_3d_g=eval_ATR20CPC_3d_g
      ENDIF
      IF(PRESENT(eval_ATR20CO2_3d_g))THEN   !for ATR20_opt. 20170801
         ALLOCATE(temp_eval_ATR20CO2_3d_g(nlon_ga, nlev_ga, ngl_ga))
         temp_eval_ATR20CO2_3d_g=eval_ATR20CO2_3d_g
      ENDIF

      if(indveva.eq.1)then
        do ilnum=1,nisland
!HY20140415-1          do ip=1,istat(ilnum,11)
!HY20140415-1            ipop=ifitpop(ilnum,ip)
!HY20140414-7          write(*,*)"from",ist,"to",ied
!HY20140414-7          write(*,*)"from",ipop,"to",ipop
!HY20140414-7          write(*,*)"from",ifitpop(ilnum,1),"to",ifitpop(ilnum,istat(ilnum,11))
            ist=nallpop-npop+1
            ied=nallpop
!HY20140415-1            call EvaPop(ipop,ipop)
            SELECT CASE(ga_option_traj_calc_g)  !20160302-2
            !F.Yin: contrail potential coverage is added "eval_cpc_g"
            !HY   : eval_cpc_g is added to other routing options in EvaPop, 20170606-1
            !
            CASE(1,5)  !Wind_opt, CPC_opt
            call EvaPop(ist, ied, eval_nwaypoints, eval_p2_ac_routes_desc, eval_philon, eval_philat, eval_zgl_geopot_3d, &
                        eval_uwind_g, eval_vwind_g, eval_v_z_g, eval_t_scb_g, eval_cpc_g)   !HY20140420-1,20140522-3,20140721-6
!HY20140415-1          enddo
            CASE(2,4,6,7,10)  !Fuel_opt, H2O_opt, Cost_opt or COC_opt   !20170222-2,20170224-2,20170316-2
           !IF(PRESENT(eval_rho_air_dry_3d_g))temp_eval_rho_air_dry_3d_g=eval_rho_air_dry_3d_g   !20170206-1 
            call EvaPop(ist, ied, eval_nwaypoints, eval_p2_ac_routes_desc, eval_philon, eval_philat, eval_zgl_geopot_3d, &
                        eval_uwind_g, eval_vwind_g, eval_v_z_g, eval_t_scb_g, eval_cpc_g, &   !HY20140420-1,20140522-3,20140721-6
                        temp_eval_rho_air_dry_3d_g)                               !20160302-2,20170206-1
!            write(*,*)'check_ga_option_traj_calc_3', ga_option_traj_calc_g        !20160304-1
            CASE(3)  !NOx_opt
           !IF(PRESENT(eval_rho_air_dry_3d_g))temp_eval_rho_air_dry_3d_g=eval_rho_air_dry_3d_g   !20170206-1 
            call EvaPop(ist, ied, eval_nwaypoints, eval_p2_ac_routes_desc, eval_philon, eval_philat, eval_zgl_geopot_3d, &
                        eval_uwind_g, eval_vwind_g, eval_v_z_g, eval_t_scb_g, eval_cpc_g, &   !HY20140420-1,20140522-3,20140721-6
                        temp_eval_rho_air_dry_3d_g, temp_eval_press_3d_g)         !20160302-2,20170206-1,20170215-3
!            write(*,*)'check_ga_option_traj_calc_3', ga_option_traj_calc_g        !20160304-1

            !Yin_20170801
            CASE(8,9)  !ATR20_opt,COSTCLIM_opt
            call EvaPop(ist, ied, eval_nwaypoints, eval_p2_ac_routes_desc, eval_philon, eval_philat, eval_zgl_geopot_3d, &
                        eval_uwind_g, eval_vwind_g, eval_v_z_g, eval_t_scb_g, eval_cpc_g,&
                        temp_eval_rho_air_dry_3d_g,temp_eval_press_3d_g,temp_eval_ATR20O3_3d_g,&
                        temp_eval_ATR20CH4_3d_g,temp_eval_ATR20H2O_3d_g,temp_eval_ATR20CPC_3d_g,&
                        temp_eval_ATR20CO2_3d_g)   
            !Yin_20170801
            CASE DEFAULT
               write(*,*)'ERROR:ga_option_traj_calc on EvaPop',ga_option_traj_calc_g
            ENDSELECT
        enddo
      else
        ist=nallpop-npop+1
        ied=nallpop
        call EvaPop(ist,ied)
      endif
      !for fuel_opt, NOx_opt, H2O_opt, Cost_opt, COC_opt
      IF(ALLOCATED(temp_eval_rho_air_dry_3d_g))DEALLOCATE(temp_eval_rho_air_dry_3d_g)   !20170206-1,20170212-3
      !for NOx_opt
      IF(ALLOCATED(temp_eval_press_3d_g))DEALLOCATE(temp_eval_press_3d_g)               !20170215-3
      !Yin_20170810+
      IF(ALLOCATED(temp_eval_ATR20O3_3d_g))DEALLOCATE(temp_eval_ATR20O3_3d_g)               !20170801
      IF(ALLOCATED(temp_eval_ATR20CH4_3d_g))DEALLOCATE(temp_eval_ATR20CH4_3d_g)               !20170801
      IF(ALLOCATED(temp_eval_ATR20H2O_3d_g))DEALLOCATE(temp_eval_ATR20H2O_3d_g)               !20170801
      IF(ALLOCATED(temp_eval_ATR20CPC_3d_g))DEALLOCATE(temp_eval_ATR20CPC_3d_g)               !20170801
      IF(ALLOCATED(temp_eval_ATR20CO2_3d_g))DEALLOCATE(temp_eval_ATR20CO2_3d_g)               !20170801
      !Yin_20170801-
      return
      end subroutine Evaluation

! *******************************
      subroutine Constraint
! *******************************
!      USE fort_1,   ONLY:testp   !HY20140331-2
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6      USE messy_airtraf_tools_ga_testfunction,     ONLY: GvalPop
      IMPLICIT NONE
      INTEGER :: ist,ied,ilnum,ip,ipop 
      if(testp.eq.'OUTE')return
            
      if(indvcnr.eq.1)then
        do ilnum=1,nisland
          do ip=1,istat(ilnum,11)
            ipop=ifitpop(ilnum,ip)
            call GvalPop(ipop,ipop,1,ncon)
          enddo
        enddo
      else
        ist=nallpop-npop+1
        ied=nallpop
        call GvalPop(ist,ied,1,ncon)
      endif

      return
      end subroutine Constraint


! **********************************
      subroutine UpdateData(ilnum)
! **********************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6       USE messy_airtraf_tools_ga_update,     ONLY: StatSetting
!HY20140404-6       USE messy_airtraf_tools_ga_archive,    ONLY: ConstrRankAdd
      IMPLICIT NONE
      INTEGER :: ilnum
      call StatSetting(ilnum)
      call ConstrRankAdd(ilnum)
      istat(ilnum,9)=1

      return
      end subroutine UpdateData


! **********************************************************
      subroutine RangeAdapt(ilnum)
! **********************************************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6       USE messy_airtraf_tools_ga_adapt,      ONLY: ArcAdapt, UpdateStatistic, GenotypeRange, RandomGeneration
!HY20140404-6       USE messy_airtraf_tools_ga_archive,    ONLY: UpdateArcInfo
!HY20140404-6       USE messy_airtraf_tools_ga_update,     ONLY: NextGen
      IMPLICIT NONE
      INTEGER :: idra(mtpop)
      INTEGER :: ilnum,nselect
      call ArcAdapt(ilnum,idra,nselect)
      call UpdateStatistic(ilnum,idra,nselect)
      call UpdateArcInfo(ilnum,idra,nselect)
      call GenotypeRange(ilnum)
      call RandomGeneration(ilnum,0)
      istat(ilnum,12)=istat(ilnum,11)
      call NextGen(ilnum)

      return
      end subroutine RangeAdapt

! **********************************************************
      subroutine AssignFitness(ilnum)
! **********************************************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6       USE messy_airtraf_tools_ga_afit,      ONLY: FitnessValue, SFitnessValue
      IMPLICIT NONE
      REAL(DP) :: fdec(maxpop,mobj)
!      dimension gdec(maxpop,mcon),murank(maxpop)
      INTEGER :: ilnum,nfitpop
      nfitpop=istat(ilnum,12)
 
      call FitnessValue(ilnum,nfitpop,fdec,mfiteval)


!cMAY      call MaxminObj(ilnum,nfitpop,fdec)
!cMAY      if(idebug.eq.1)call Debug_Obj(ilnum,istat(ilnum,12),mtdeb)
!cMAY
!cMAY      if(ncon.gt.0)then
!cMAY        call MaxminGval(ilnum,nfitpop,gdec)
!cMAY        call ConstrHandling(ilnum,nfitpop,gdec,icnr)
!cMAY        print *,'Pareto Ranking'
!cMAY        if(mprank.eq.1)
!cMAY     &      call ConstrParetoRanking(ilnum,nfitpop,fdec,gdec,murank,icnr)
!cMAY        if(mprank.eq.2)
!cMAY     &      call ConstrFrontRanking(ilnum,nfitpop,fdec,gdec,murank,icnr)
!cMAY        if(idebug.eq.1)call Debug_ConstrRank(ilnum,gdec,mtdeb)
!cMAY      else
!cMAY        if(mprank.eq.1)call ParetoRanking(ilnum,fdec,murank)
!cMAY        if(mprank.eq.2)call FrontRanking(ilnum,fdec,murank)
!cMAY      endif
!CHECK
!cMAY      if(idebug.eq.1)call Debug_Fitness(ilnum,mtdeb)

!cMAY      if(mfiteval.eq.0)call RankBasedFitness(ilnum)
!cMAY      if(mfiteval.eq.1)call AverageFitness(ilnum,murank) 
!CHECK
!cMAY      if(idebug.eq.1)call Debug_Fitness(ilnum,mtdeb)

      if(ishare.gt.0)then
        call SFitnessValue(ilnum,nfitpop,fdec)

!cMAY        numpop=istat(ilnum,12)
!cMAY        if(isig2n.eq.0)numpop=istat(ilnum,11)
!cMAY        call NormSigmaShare(numpop,nobj,sshare)
!cMAY        if(msdist.eq.1)call NormNicheCount(ilnum,fdec,sshare)
!cMAY        if(msdist.eq.2)call NormNicheCountDV(ilnum,0)
!cMAY        if(msdist.eq.3)call NormNicheCountDV(ilnum,1)
!cMAY        if(msdist.eq.4)call NormNicheCountDV(ilnum,2)
!cMAY        call SharedFitness(ilnum)
!cMAY        if(idebug.eq.1)call Debug_Niche(ilnum,sshare,mtdeb) 
      endif
!HY20140328-1      if(idebug.eq.1)call Debug_Fitness(ilnum,nfitpop,mtdeb)
 
      return
      end subroutine AssignFitness

! ********************************************
      subroutine Reproduction(ilnum)
! ********************************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6       USE messy_airtraf_tools_ga_afit,     ONLY: PopSorting
!HY20140404-6       USE messy_airtraf_tools_ga_archive,      ONLY: ArcParent, ArcNBest
!HY20140404-6       USE messy_airtraf_tools_ga_selection,    ONLY: Selection, MatingPool
      IMPLICIT NONE
      INTEGER, DIMENSION(:) :: iorder(maxpop)
      INTEGER :: ilnum,newpop,nfitpop
      nfitpop=istat(ilnum,12)
 
      call PopSorting(ilnum,nfitpop,iorder)

!cMAY      if(iarchiv.ge.2)then
      if(iarchiv.gt.0)then
!cMAY        if( (iarchiv.eq.2) .or. (istat(ilnum,1).eq.2) )then
        if( (iarchiv.eq.1) .or. (istat(ilnum,1).eq.2) )then
!C1.0.1          print *,'ARCHIVING',igen,istat(ilnum,1)
!C1.2.0          print *,'ARCHIVING',igen+1,istat(ilnum,1)
          print *,'ARCHIVING',iter,istat(ilnum,1)
          call ArcParent(ilnum)
          call AssignFitness(ilnum)
          call PopSorting(ilnum,nfitpop,iorder)
        endif
      endif
 
      call Selection(ilnum)

      newpop=istat(ilnum,11)
!HY20140328-1      if(idebug.eq.1)call Debug_MatingPool(ilnum,newpop,mtdeb)
 
      call MatingPool(ilnum,newpop)
!HY20140328-1      if(idebug.eq.1)call Debug_Obj(ilnum,newpop,mtdeb)

!      if( (iarchiv.eq.1) .and. (istat(ilnum,2).gt.0) )then
!cMAY      if(iarchiv.eq.1)call ArcNBest(ilnum)
      if(ibnarc.eq.1)call ArcNBest(ilnum)
!      endif
 
      return
      end subroutine Reproduction

! ******************************************************
      subroutine NextPop(ilnum)
! ******************************************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6       USE messy_airtraf_tools_ga_selection,     ONLY: PairRandom, PairOrder
!HY20140404-6       USE messy_airtraf_tools_ga_cross,      ONLY: Crossover
!HY20140404-6       USE messy_airtraf_tools_ga_mutation,     ONLY: Mutation
!HY20140404-6       USE messy_airtraf_tools_ga_adapt,      ONLY: DecodePop
!HY20140404-6       USE messy_airtraf_tools_ga_testfunction,     ONLY: GvalPop
!HY20140404-6       USE messy_airtraf_tools_ga_update,     ONLY: NextGen, NextFPop
      IMPLICIT NONE
      INTEGER  :: matepair(maxpop,3)
      INTEGER  :: istpop,ncount,numloop,ilnum,icn,igst,iged,inum,ioff1,ioff2,mate1,mate2,maxl
      INTEGER  :: newpop,numpair
      newpop  = istat(ilnum,11)
      istpop  = istat(ilnum,10)
      numpair = newpop/2
      ncount  = numpair+1
!      ncount  = 1

      igst    = imoga(ilnum,7)
      iged    = imoga(ilnum,8)

!      nsc = 1
!      if(mrpair.gt.0)then
!        call PairRandom(ilnum,matepair,numpair,newpop)
!      else
!        call PairOrder(ilnum,matepair,numpair,newpop,mrandp)
!      endif 


      inum = 1

      do while(inum.le.newpop)
        numloop=0
 10     continue

        if((numloop.gt.0).and.(mrpair.eq.2))ncount=ncount+1 
        if(numloop.gt.loopmax)ncount=ncount+1
        if(ncount.gt.numpair)then
          ncount=1
          if(mrpair.gt.0)then
            call PairRandom(ilnum,matepair,numpair,newpop)
          else
            maxl=int(rnpair*newpop)+1
            if(numloop.gt.loopmax)maxl=max(maxl,numloop)
            call PairOrder(ilnum,matepair,numpair,newpop,maxl)
          endif
        endif

!CD        if(mrpair.eq.2)then
!CD          if(nsc.eq.0)ncount=ncount+1 
!CD          if(ncount.gt.numpair)then
!CD            ncount=1
!CD            call PairRandom(ilnum,matepair,numpair,newpop)
!CD          endif
!CD        endif

        ioff1 = istpop + inum -1
        ioff2 = ioff1 + 1
        mate1 = matepair(ncount,1) 
        mate2 = matepair(ncount,2) 
!CHECK
!        write(*,*)'PopNo.',ioff1,ioff2,mate1,mate2,ncount
       
        if(idebug.eq.1)write(mtdeb,*)'PopNo.',ioff1,mate1,mate2

        numloop=numloop+1
!        nsc=0

        call Crossover(ilnum,ioff1,ioff2,mate1,mate2)

        call Mutation(ilnum,ioff1)
        call DecodePop(ilnum,ioff1,ioff1,1,ndv)
!cDS        if(icnrgeo.eq.1)then
        if(igst.gt.0)then
          call GvalPop(ioff1,ioff1,igst,iged)
          do icn=igst,iged
            if(gvall(ioff1,icn).gt.gvaltol)goto 10
          enddo
        endif

        call Mutation(ilnum,ioff2)
        call DecodePop(ilnum,ioff2,ioff2,1,ndv)
!cDS        if(icnrgeo.eq.1)then
        if(igst.gt.0)then
          call GvalPop(ioff2,ioff2,igst,iged)
          do icn=igst,iged
            if(gvall(ioff2,icn).gt.gvaltol)goto 10
          enddo
        endif

        istat(ilnum,13)=istat(ilnum,13)+1
        idarc(ilnum,istat(ilnum,13))=ioff1
        istat(ilnum,13)=istat(ilnum,13)+1
        idarc(ilnum,istat(ilnum,13))=ioff2
         
!        nsc=1
        ncount=ncount+1
        inum = inum + 2
      enddo

      
      if(ibestn.gt.0)call NextFPop(ilnum)
      call NextGen(ilnum)
      
      return
      end subroutine NextPop
  

! ************************************
      subroutine Judgement(condition)
! ************************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      LOGICAL :: condition
  
!C1.0.1      if((igen.ge.ngend).or.(iman.eq.1))condition=.false.
!C1.2.0      if((igen+1.ge.ngend).or.(iman.eq.1))condition=.false.
      if((iter.ge.ngend).or.(iman.eq.1))condition=.false.

      igen = igen + 1
!C1.2.0
      iter = iter + 1
      nallpop = nallpop + npop

      return
      end subroutine Judgement

! *********************************
      subroutine PostProcess
! *********************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6       USE messy_airtraf_tools_ga_iofile,     ONLY: OutputObj
      IMPLICIT NONE
      INTEGER :: ilnum,ipst,iped
      ipst=nallpop-npop+1
      iped=nallpop
!HY20140327-2   call OutputGPtype(ipst,iped,mtgen,0) 
!HY20170327-1   call OutputGPtype(ipst,iped,mteva,1) 

!HY20140327-1   call OutputSystemPtype(ipst,iped,mtsys,1) 

!      ist = iflag(9) * (1-imtall) * (igen-1)
!      ied = iflag(9) * (igen-1)


! ---THINK
      do ilnum=1,nisland
        ipst=istat(ilnum,13)-istat(ilnum,11)*2+1
        iped=istat(ilnum,13)-istat(ilnum,11)
!HY20140327-1     call OutputObj(ilnum,ipst,iped,1,mtobj)

        ipst=1
        iped=istat(ilnum,13)-istat(ilnum,11)
!HY20140327-2     call OutputObj(ilnum,ipst,iped,0,mtprt)

        if(idebug.eq.1) call OutputObj(ilnum,ipst,iped,1,72)
      enddo

! ---THINK
      ipst=1
      iped=nallpop-npop
!HY20140327-2   call OutputAllData(ipst,iped,mtall)

      return
      end subroutine PostProcess

!------------------------------------------------------------
!
! ORIGINAL: messy_airtraf_tools_ga_initial.f90 
!
! ***********************************************
      subroutine CommonSetting(nobjt,ncont,ndvt,iccset)
! ***********************************************
      USE messy_airtraf_tools_ga_input             !HY20140331-5,20140402-2 
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      INTEGER :: iccset,ndvt,ncont,nobjt
!HY20140401-7      include 'ga.par'
!HY20140331-5      NAMELIST/input/testp, imaxmin, ngend, &
!HY       iman, idman, indveva, indvcnr, idebug, & 
!HY       icnrgeo, icnrrnk, icnrvio, &
!HY       mtout, mtparam, mtgop, mtgen, mteva, ipout, mtobj, &
!HY       mtprt, mtall, mtdeb, imtall, imemory, irank1, infopt, &
!HY       mfiteval,  ishare, msdist, iallsh, shalpha, irksp, mprank, &
!HY       ictrlelt, rkcoef, mselection, ibestn, isig2n, alpdmpar, &
!HY       ialpnorm, iarchiv, arcratio, arclimit, aratio, icrobd, &
!HY       ismpop,xratio,iarbase,arpopratio,ararcratio,mrpair,rnpair, &
!HY       loopmax, isetdes, isubar, mfitsel, clr, iadinit, gvaltol, &
!HY       mtsys, ibnarc, idvsm

! .....................................................
! ..... General Information ...........................
!  >>>> testp   ... Name of Test Problem
!  >>>> imaxmin  1 ... Maximization Problem
!  >>>> imaxmin  0 ... Minimization Problem
!  >>>> ngend   ... Number of Final Generation 
!  >>>> iman     0 ... Auto Mode Optimization 
!  >>>> iman     1 ... Manual Mode Optimization
!  >>>> idman    0 ... [[ Under Construction ]] 
!  >>>> idman    1 ... Manual Operator Input
!  >>>> indveva  0 ... All evaluations are performed simultaneously.
!  >>>> indveva  1 ... Each evaluation is performed separately.
!  >>>> indvcnr  0 ... All constraints are calculated simultaneously.
!  >>>> indvcnr  1 ... Each constraint is calcualted separately.
!  >>>  idebug   0 ... Standard Mode
!  >>>> idebug   1 ... Debug Mode 
      testp   =  testp_input_d   !HY20140331-5
      imaxmin =  imaxmin_input_d
      ngend   =  ngend_input_d
      iman    =  iman_input_d
      idman   =  idman_input_d
      indveva =  indveva_input_d
      indvcnr =  indvcnr_input_d
      idebug  =  idebug_input_d
      isetdes =  isetdes_input_d
      isubar  =  isubar_input_d
      iadinit =  iadinit_input_d
! ......................................................
! ..... I/O Data .......................................
!  >>>> mtout    ... Normal Output
!  >>>>>>> infopt   0  ... No Optimization Information
!  >>>>>>> infopt   1  ... Optimization Information Output
!  >>>> mtparam  ... Parameter File
!  >>>> mtgop    ... Output Operator Data
!  >>>> mtgen    ... Latest Population Data (Gtype)
!  >>>> mteva    ... Latest Population Data (Ptype)
!  >>>>>>> ipout 0  ... Single Phenotype Data
!  >>>>>>> ipout 1  ... Separate Phenotype Data
!  >>>> mtobj    ... Latest Objective Function
!  >>>> mtprt    ... Latest Pareto Optimal
!  >>>> mtall    ... Archive File (imem=1)
!  >>>>>>> imtall   0 ... Final Data is preserved.
!  >>>>>>> imtall   1 ... All Data is preserved.
!  >>>>>>> imemory  0 ... Present Population Data is stored.
!  >>>>>>> imemory  1 ... All Population Data is stored.[*1]
!  >>>>>>> irank1   0 ... Ranking for All Population
!  >>>>>>> irank1   1 ... Ranking Only for Rank1 Population [*2]
!  >>>> mtdeb    ... Debug Output
      mtout   = mtout_input_d     !HY20140331-5 
      mtparam = mtparam_input_d 
      mtgop   = mtgop_input_d 
      mtgen   = mtgen_input_d
      mteva   = mteva_input_d
      ipout   = ipout_input_d 
      mtobj   = mtobj_input_d
      mtprt   = mtprt_input_d
      mtall   = mtall_input_d
      mtsys   = mtsys_input_d
      imtall  = imtall_input_d 
      imemory = imemory_input_d 
      irank1  = irank1_input_d 
      infopt  = infopt_input_d 
      mtdeb   = mtdeb_input_d
! ......................................................
! ..... Genetic Operator ...............................
!  [ Constraint Handling ]
!  >>>> icnrgeo  0 ... No Geometry Constraints Treatment
!  >>>> icnrgeo  1 ... Prevention of Infeasible Region 
!  >>>>>>> gvaltol = 0. ... All designs are always feasible
!  >>>>>>> gvaltol > 0. ... Tolerance for Infeasible Region
!  >>>>>>>                   (ngeocon > 0 is required)
!  >>>> icnrrnk  1 ... Constrained Domination
!  >>>> icnrrnk  2 ... Constrained Domination (1-gvalue)
!  >>>> icnrvio  0 ... Violated POP is not considered in Adaptaion
!  >>>> icnrvio  1 ... Min. violated POP is chosen in Adaptation
!  [ Selection ]
!  >>>> mfiteval 1 ... Fitness Based on Rank (=1/rank)
!  >>>> mfiteval 2 ... Average Fitness
!  >>>> mfiteval 3 ... Fitness Based on Rank (by Michel...)
!  >>>> ishare   0 ... No Sharing Function
!  >>>> ishare   1 ... Sharing Function [for (Rt)]
!  >>>> ishare   2 ... Sharing Function [for (Rt)&(Pt+1)]
!  >>>>>>> msdist   1 ... Sharing Distance by Objective Function Space
!  >>>>>>> msdist   2 ... Sharing Distance by Present Design Space
!  >>>>>>> msdist   3 ... Sharing Distance by Original Design Space
!  >>>>>>> iallsh   0 ... Sharing for Same Rank Population
!  >>>>>>> iallsh   1 ... Sharing for All Population
!  >>>>>>> shalpha (r>0)  ... Sharing Coefficients, ALPHA(=0.5/1.0/2.0)
!  >>>>>>> irksp    0 ... Order by Shared Fitness Only
!  >>>>>>> irksp    1 ... Order by Rank with Sharing 
!  >>>> mprank   1 ... Fonseca's Pareto Ranking
!  >>>> mprank   2 ... Goldberg's Pareto Ranking (Front)
!  >>>>>>> ictrlelt 0 ... No Control (mprank=2 is required)
!  >>>>>>> ictrlelt 1 ... Controled Elitism (Deb) [*3]
!  >>>>>>> rkcoef (0<r<1) ... Controled Coefficient
!  >>>> mselection  0 ... No Selection (Mating Pool = Pt+1)
!  >>>> mselection  1 ... SUS Selection 
!  >>>> mselection  2 ... Roulette Wheel Selection (RWS)
!  >>>> ibestn   0 ... No Elitism (Always N set)
!  >>>> ibestn   1 ... Best N Strategy
!  >>>>>>> isig2n   0 ... SigmaShare for N Pop 
!  >>>>>>> isig2n   1 ... SigmaShare for 2N Pop 
!  [ ALPHA-Domination ]
!  >>>> alpdmpar (r=0) ... Non Dominated Strategy
!  >>>> alpdmpar (r>0) ... Alpha-Domination Strategy
!  >>>>>>> ialpnorm 0 ... No Normalization (Manual Usage)
!  >>>>>>> ialpnorm 1 ... Normalization by All Max-Min
!  >>>>>>> ialpnorm 2 ... Normalization by Present Max-Min
!  [ Archiving ] 
!  >>>> iarchiv     1  ... Archiving for Next Parent Set
!  >>>> iarchiv     2  ... Archiving for Actual Parent Set (MATING POOL)
!  >>>> iarchiv     3  ... Archiving for Actual Parent Set (Next Adaptation)
!  >>>> arcratio (r>0) ... Archiving Ratio for NPOP
!  >>>> arclimit (r>=0) ... Archiving Size = (arclimit)*NPOP
!  [ Adaptive Range Algorithms ]
!  >>>> aratio   (r=0) ... Only Rank 1 POPs are sampled
!  >>>> aratio   (r>0) ... Sampling population Ratio
!  [ Crossover ]
!  >>>> icrobd   0 ... Both Children have to exist design region 
!  >>>> icrobd   1 ... If new pop is over the bound, then pop is bound 
! -------------------------------------------------------
      icnrgeo = icnrgeo_input_d    !HY20140331-5
      gvaltol = gvaltol_input_d
      icnrrnk = icnrrnk_input_d
      icnrvio = icnrvio_input_d  
! --------------------
      mfiteval = mfiteval_input_d  
      ishare   = ishare_input_d 
      msdist   = msdist_input_d  
      iallsh   = iallsh_input_d  
      shalpha  = shalpha_input_d
      irksp    = irksp_input_d   
      mprank   = mprank_input_d  
      ictrlelt = ictrlelt_input_d  
      rkcoef   = rkcoef_input_d
      mselection = mselection_input_d
      ibestn   = ibestn_input_d  
      isig2n   = isig2n_input_d  
! --------------------
!      mfiteval =   1
      alpdmpar = alpdmpar_input_d
      ialpnorm = ialpnorm_input_d  
      iarchiv  = iarchiv_input_d  
      ibnarc   = ibnarc_input_d  
      arcratio = arcratio_input_d
      arclimit = arclimit_input_d
      aratio   = aratio_input_d
      icrobd   = icrobd_input_d  
! ----------------------------------------
!  NEW NEW NEW
      ismpop   = ismpop_input_d
      xratio   = xratio_input_d
      iarbase  = iarbase_input_d
      arpopratio = arpopratio_input_d
      ararcratio = ararcratio_input_d
      mrpair   = mrpair_input_d
      rnpair   = rnpair_input_d
      loopmax  = loopmax_input_d
      mfitsel  = mfitsel_input_d
      clr      = clr_input_d
      idvsm    = idvsm_input_d 
! ----------------------------------------

      if(iccset.eq.0)then
!HY20140331-5        read(5,input)
!HY20140331-5        if(idebug.eq.1)write(6,input)
      else
        testp='OUTE'
        ngend=8
        iman=1
        isetdes=1
        mtgop=0
        infopt=1
      endif

      if(infopt.eq.1)write(mtout,'(A)')testp
      testcode(1:3) = testp 

      if(imem.eq.0)then
        imtall=0
        imemory=0
      endif

      if(mfitsel.eq.0)mfitsel=mfiteval

      if(testp.eq.'OUTE')isetdes=1

!  ********* Input Optimisation Information from (mtgen) ******
!C1.2.0      read(mtgen,*)igen, npop, nisland
!HY20140329-2      read(mtgen,*)iter, npop, nisland
!HY20140329-2      read(mtgen,*)nobj, ncon, ndv
!HY20140329-2      read(mtgen,*)idum, iy
      iter = iter_fort_42 
      npop = npop_fort_42
      nisland = nisland_fort_42
      nobj = nobj_fort_42
      ncon = ncon_fort_42
      ndv = ndv_fort_42
      idum = idum_fort_42
      iy = iy_fort_42
!C1.2.0 -ref-
      igen=iter-1

      nobjt = nobj
      ncont = ncon
      ndvt  = ndv

!  ********* IFLAG **************
!      if((imtall.eq.1).and.(igen.gt.0))then
!        iflag(7)=1
!      else
!        iflag(7)=0
!      endif
      iflag(8)=imtall
      iflag(9)=imemory

      return
      end subroutine CommonSetting

! ********************************************
      subroutine CheckDesign(nobjt,ncont,ndvt)
! ********************************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!C1.2.0      if(igen.ge.0)then
      IMPLICIT NONE
      INTEGER :: nobjt,ncont,ndvt 
      if(iter.ge.1)then
        if(nobjt.ne.nobj) stop 'Check! Wrong Generation File [ nobj ]' 
        if(ncont.ne.ncon) stop 'Check! Wrong Generation File [ ncon ]' 
        if( ndvt.ne. ndv) stop 'Check! Wrong Generation File [ ndv ]' 
      endif
   
      if(ngend.ge.mgend) stop 'Check! [ ngend / mgend ]' 
      if( npop.gt. mpop) stop 'Check! [ npop / mpop ]'
      if(  ndv.gt.  mdv) stop 'Check! [ ndv / mdv ]'
      if( nobj.gt. mobj) stop 'Check! [ nobj / mobj ]'
      if( ncon.gt. mcon) stop 'Check! [ ncon / mcon ]'
      if(nisland.gt.misland) stop 'Check! [ nisland / misland ]'

!      if(imtall.gt.imemory) stop 'Check! [ imtall / imem ]'

      return
      end subroutine CheckDesign

! ************************************
      subroutine InputOperator
! ************************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6       USE messy_airtraf_tools_ga_iofile,     ONLY: InputByHand
      IMPLICIT NONE
      if(idman.eq.1)then
        call InputByHand
      else
!check
!       call InputByParam
       stop
      endif

      return
      end subroutine InputOperator

! ************************************
      subroutine SetParam
! ************************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6       USE messy_airtraf_tools_ga_testfunction,     ONLY: DesignMAXMIN
      IMPLICIT NONE
      INTEGER :: ilnum,iofi,iofj,idv
!  ********* CONSTRAINTS **************
      if(ncon.lt.1)then
        icnrrnk = 0
        icnrvio = 0
        nlinear = 0
        ngeocon = 0
      endif
      if(ngeocon.lt.1)then
        icnrgeo = 0
      endif  

!  *** ELITIST .. FUTURE WORK FOR SINGLE OBJ  *****
!HY20140407-3      sumn = 0   
      nelite = 0
      do ilnum=1,nisland
        if(nelt(ilnum).gt.0)nelite = 1
      enddo

!  ******* ALPHA Domination ************
      do iofi=1,nobj
        do iofj=1,nobj
          alpdm(iofi,iofj) = alpdmpar
        enddo
        alpdm(iofi,iofi) = 1.0_dp  !20140619-1
      enddo 

!  ******** DESIGN SPACE, INITIAL VALUES  **********
      do ilnum=1,nisland
        do idv=1,ndv
          geub(ilnum,idv) = 1.0_dp  !20140619-1
          gelb(ilnum,idv) = 0.0_dp  !20140619-1
        enddo
        call DesignMAXMIN(ilnum)
      enddo

!  ********* NALLPOP **************
      nallpop = 0

!  ********* ISTAT **************
      do ilnum=1,nisland
        istat(ilnum,1)=0
        istat(ilnum,2)=0
        istat(ilnum,3)=0

        istat(ilnum,4)=imoga(ilnum,1)
        istat(ilnum,5)=imoga(ilnum,2)
!        if(icnrgeo.eq.1)then
!          istat(ilnum,6)=imoga(ilnum,7)
!          istat(ilnum,7)=imoga(ilnum,8)
!        else
          istat(ilnum,6)=imoga(ilnum,5)
          istat(ilnum,7)=imoga(ilnum,6)
!        endif

        istat(ilnum,9)=0
        istat(ilnum,10)=0
        istat(ilnum,11)=0
        istat(ilnum,12)=0
        istat(ilnum,13)=0
        istat(ilnum,14)=0
      enddo
 
!  ********* IMOGA **************
!      do ilnum=1,nisland
!        imoga(ilnum,1)=1
!        imoga(ilnum,2)=nobj
!        imoga(ilnum,3)=1
!        imoga(ilnum,4)=nobj
!        imoga(ilnum,9)=0
!        imoga(ilnum,10)=0
!        if(ncon.gt.0)then
!          imoga(ilnum,9)=1
!          imoga(ilnum,10)=ncon
!        endif
!      enddo  

      return
      end subroutine SetParam

! ************************************
      subroutine InitialParam
! ************************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      REAL(DP) :: r
      INTEGER  :: irm,i,idv,ilnum,ip    !HY20140407-3
!HY20140404-6       USE messy_airtraf_tools_ga_function,     ONLY: ran1
! << Initialized Parameters >>
      igen = 0
!C1.2.0
      iter = 1

! << Initialization of Random Number >>
      irm  = idum
      idum = 0
!      igen = 0
      if(irm.gt.1000000)irm=irm/100000
      do i=1,irm
        r = ran1(idum)
      enddo

! << Island Information >>
      do ilnum=1,nisland
        ip = npop/nisland
        istat(ilnum,10)=(ilnum-1)*ip+1
        istat(ilnum,11)=ip
        istat(ilnum,12)=0
        istat(ilnum,13)=0
! << Parameter for Design Variables >>
        do idv=1,ndv
          pave(ilnum,idv)=0.5_dp*(phlb(idv)+phub(idv))
          stdv(ilnum,idv)=pave(ilnum,idv)-phlb(idv)
          flat(ilnum,idv)=stdv(ilnum,idv)
          pout(ilnum,idv)=0.0_dp
        enddo
      enddo

      return
      end subroutine InitialParam

! **********************************
      subroutine InitialPop
! **********************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6       USE messy_airtraf_tools_ga_adapt,     ONLY: RandomGeneration
!HY20140404-6       USE messy_airtraf_tools_ga_update,    ONLY: InitialObj
      IMPLICIT NONE
      INTEGER :: ist,ied,ilnum,ipop
      do ilnum=1,nisland
! << Restriction of Initial Boundary >>
!        call AdjustDV(ilnum,0)

        if(iadinit.eq.1)then
          call AdjustRandomGeneration(ilnum)
        else
          call RandomGeneration(ilnum,0)
        endif

        ist = istat(ilnum,10)
        ied = istat(ilnum,10)+istat(ilnum,11)-1

!        do while(ipop.le.ied)
!          do idv=1,ndv
!            gtype = ran1(idum)
!            if(gtype.lt.0.)gtype=0.
!            if(gtype.gt.1.)gtype=1.
!            ptype = stdv(ilnum,idv) * (2.*(gtype-0.5)) + pave(ilnum,idv)
!            if(ptype.lt.phlb(idv))ptype=phlb(idv)
!            if(ptype.gt.phub(idv))ptype=phub(idv)
!            ptall(ipop,idv) = ptype
!            gtall(ipop,idv) = gtype
!          enddo
!        enddo

        do ipop=ist,ied
!C1.2.0          ipinfo(ipop,1)=igen
          ipinfo(ipop,1)=iter
          ipinfo(ipop,2)=ipop
          ipinfo(ipop,3)=ilnum
          ipinfo(ipop,4)=0  
          istat(ilnum,12)=istat(ilnum,12)+1
          ifitpop(ilnum,istat(ilnum,12))=ipop
        enddo

        call InitialObj(ist,ied)
        istat(ilnum,10)=istat(ilnum,10)+npop

!        call AdjustDV(ilnum,1)

      enddo
      nallpop=npop

          
      return
      end subroutine InitialPop

! ********************************
      subroutine OutputInitial
! ********************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      INTEGER :: ipst,iped
      ipst=nallpop-npop+1
      iped=nallpop

!HY20140327-2   call OutputGPtype(ipst,iped,mtgen,0)
!HY20140327-1   call OutputGPtype(ipst,iped,mteva,1)
!HY20140327-1   call OutputSystemPtype(ipst,iped,mtsys,1)

      return
      end subroutine OutputInitial

! **********************************************
      subroutine AdjustRandomGeneration(ilnum)
! **********************************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6       USE messy_airtraf_tools_ga_adapt,     ONLY: RandomGeneration
      IMPLICIT NONE
      REAL(DP) :: prob(mdv),plo(mdv),pup(mdv),press(mdv)
      REAL(DP) :: cen(mdv),wid(mdv)
      REAL(DP) :: ptmp1,ptmp2,ptmp3
      INTEGER :: il,numtype,ilnum,num,idv,ist,ied,igst,iged,inst,ined,itype,numpop
      ist=istat(ilnum,10)
      ied=istat(ilnum,10) + istat(ilnum,11) - 1

      igst=imoga(ilnum,7)
      iged=imoga(ilnum,8)


      open(10,file='armoga.ini')
      if(ilnum.gt.1)then
        do il=1,ilnum-1
          read(10,*)numtype
          do num=1,numtype
            read(10,*)itype,numpop
            if(itype.eq.1)then
              do idv=1,ndv
                read(10,*)ptmp1,ptmp2,ptmp3
              enddo
            elseif(itype.eq.2)then
              do idv=1,ndv
                read(10,*)ptmp1
              enddo
            endif
          enddo
        enddo
      endif

      read(10,*)numtype
      do num=1,numtype
        read(10,*)itype,numpop
        inst=ist
        ined=ist+numpop-1
        if(itype.eq.1)then
          do idv=1,ndv
            read(10,*)prob(idv),plo(idv),pup(idv)
          enddo
          call DisturbGen(ilnum,inst,ined,prob,plo,pup)
        elseif(itype.eq.2)then
          do idv=1,ndv
            read(10,*)press(idv)
            plo(idv)=pave(ilnum,idv)-flat(ilnum,idv)
            pup(idv)=pave(ilnum,idv)+flat(ilnum,idv)
!            print *,plo(idv),pup(idv)
          enddo
          call FeasibleGen(ilnum,inst,ined,press,plo,pup)
        elseif(itype.eq.3)then
          do idv=1,ndv
            read(10,*)prob(idv),cen(idv),wid(idv)
          enddo
          call NormalGen(ilnum,inst,ined,prob,cen,wid)
        endif
        ist=ined+1
      enddo
!c 1.1.3
      if(ined.lt.ied) call RandomGeneration(ilnum,ist)

      close(10)

      return
      end subroutine AdjustRandomGeneration


      subroutine DisturbGen(ilnum,ist,ied,prob,plo,pup)
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6       USE messy_airtraf_tools_ga_testfunction,     ONLY: GvalPop
!HY20140404-6       USE messy_airtraf_tools_ga_iofile,     ONLY: ErrorPtype
!HY20140404-6       USE messy_airtraf_tools_ga_adapt,     ONLY: Encode
!HY20140404-6       USE messy_airtraf_tools_ga_function,     ONLY: ran1
      IMPLICIT NONE
      REAL(DP) :: prob(mdv),plo(mdv),pup(mdv)
      REAL(DP) :: flo,fup,gtype,pot,ptype,sdv
      INTEGER  :: ist,ied,ilnum,icn,idv,ipop,igst,iged
      igst=imoga(ilnum,7)
      iged=imoga(ilnum,8)

      do ipop=ist,ied
        do idv=1,ndv
          ptall(ipop,idv)=0.0_dp   !20140619-1
        enddo
      enddo

      do ipop=ist,ied
 20     continue
        do idv=1,ndv
          ptype=0.0_dp   !20140619-1
          if(ran1(idum).le.prob(idv))then
            ptype = (pup(idv)-plo(idv))*ran1(idum) + plo(idv)
!CAUTION
          endif
          ptall(ipop,idv)=ptype
        enddo
        if(igst.gt.0)then
          call GvalPop(ipop,ipop,igst,iged)
          do icn=igst,iged
            print *,'DIST',ipop,gvall(ipop,icn)
            if(gvall(ipop,icn).gt.gvaltol)goto 20
          enddo
        endif
        do idv=1,ndv
          ptype=ptall(ipop,idv)
          pot=pout(ilnum,idv)
          sdv=stdv(ilnum,idv)
          flo=pave(ilnum,idv)-flat(ilnum,idv)
          fup=pave(ilnum,idv)+flat(ilnum,idv)
          gtype=Encode(ptype,pot,sdv,flo,fup)
          gtall(ipop,idv) = gtype
          if((gtype.lt.gelb(ilnum,idv)).or.(gtype.gt.geub(ilnum,idv))) &
           call ErrorPtype(ilnum,ipop,idv)
        enddo
        istat(ilnum,13)=istat(ilnum,13)+1
        idarc(ilnum,istat(ilnum,13))=ipop 
      enddo

      return
      end subroutine DisturbGen



      subroutine FeasibleGen(ilnum,ist,ied,press,plo,pup)
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6       USE messy_airtraf_tools_ga_testfunction,     ONLY: GvalPop
!HY20140404-6       USE messy_airtraf_tools_ga_iofile,      ONLY: ErrorPtype
!HY20140404-6       USE messy_airtraf_tools_ga_adapt,    ONLY: Encode
!HY20140404-6       USE messy_airtraf_tools_ga_function,     ONLY: ran1
      IMPLICIT NONE
      REAL(DP) :: press(mdv),plo(mdv),pup(mdv)
      REAL(DP) :: flo,fup,gtype,pot,ptype,sdv
      INTEGER  :: ilnum,ist,ied,icn,idv,ipop,igst,iged
      igst=imoga(ilnum,7)
      iged=imoga(ilnum,8)

      do ipop=ist,ied
        do idv=1,ndv
          ptall(ipop,idv)=0.0_dp   !20140619-1
        enddo
      enddo

      do ipop=ist,ied
        do idv=1,ndv
          ptype = (pup(idv)-plo(idv))*ran1(idum) + plo(idv)
          ptall(ipop,idv)=ptype
        enddo
 20     continue
        do idv=1,ndv
          ptall(ipop,idv)=ptall(ipop,idv)*press(idv)
        enddo
        if(igst.gt.0)then
          call GvalPop(ipop,ipop,igst,iged)
          do icn=igst,iged
            print *,'FEAS',ipop,gvall(ipop,icn)
            if(gvall(ipop,icn).gt.gvaltol)goto 20
          enddo
        endif
        do idv=1,ndv
          ptype=ptall(ipop,idv)
          pot=pout(ilnum,idv)
          sdv=stdv(ilnum,idv)
          flo=pave(ilnum,idv)-flat(ilnum,idv)
          fup=pave(ilnum,idv)+flat(ilnum,idv)
          gtype=Encode(ptype,pot,sdv,flo,fup)
          gtall(ipop,idv) = gtype
          if((gtype.lt.gelb(ilnum,idv)).or.(gtype.gt.geub(ilnum,idv))) &
           call ErrorPtype(ilnum,ipop,idv)
        enddo
        istat(ilnum,13)=istat(ilnum,13)+1
        idarc(ilnum,istat(ilnum,13))=ipop 
      enddo

      return
      end subroutine FeasibleGen

      subroutine NormalGen(ilnum,ist,ied,prob,cen,wid)
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6       USE messy_airtraf_tools_ga_testfunction,     ONLY: GvalPop
!HY20140404-6       USE messy_airtraf_tools_ga_iofile,     ONLY: ErrorPtype
!HY20140404-6       USE messy_airtraf_tools_ga_adapt,     ONLY: Encode
!HY20140404-6       USE messy_airtraf_tools_ga_function,     ONLY: ran1, xmin1
      IMPLICIT NONE
      REAL(DP) :: prob(mdv),cen(mdv),wid(mdv)
      REAL(DP) :: flo, fup,u,gtype,fnor,pot,ptype,sdv
      INTEGER  :: ilnum,ist,ied,igst,iged,icn,idv,ipop
      igst=imoga(ilnum,7)
      iged=imoga(ilnum,8)

      do ipop=ist,ied
        do idv=1,ndv
          ptall(ipop,idv)=0.0_dp   !20140619-1
        enddo
      enddo

      do ipop=ist,ied
 20     continue
        do idv=1,ndv
 21       continue
          ptype=cen(mdv)
          if(ran1(idum).le.prob(idv))then
            u = xmin1(ran1(idum),fnor)/xminpar  !xmin1: function 
            ptype = u*wid(idv)+cen(idv)
            flo=pave(ilnum,idv)-flat(ilnum,idv)
            fup=pave(ilnum,idv)+flat(ilnum,idv)
            if (ptype.lt.flo.or.ptype.gt.fup)goto 21
!CAUTION
          endif
          ptall(ipop,idv)=ptype
        enddo
        if(igst.gt.0)then
          call GvalPop(ipop,ipop,igst,iged)
          do icn=igst,iged
            print *,'DIST',ipop,gvall(ipop,icn)
            if(gvall(ipop,icn).gt.gvaltol)goto 20
          enddo
        endif
        do idv=1,ndv
          ptype=ptall(ipop,idv)
          pot=pout(ilnum,idv)
          sdv=stdv(ilnum,idv)
          flo=pave(ilnum,idv)-flat(ilnum,idv)
          fup=pave(ilnum,idv)+flat(ilnum,idv)
          gtype=Encode(ptype,pot,sdv,flo,fup)
          gtall(ipop,idv) = gtype
          if((gtype.lt.gelb(ilnum,idv)).or.(gtype.gt.geub(ilnum,idv))) &
           call ErrorPtype(ilnum,ipop,idv)
        enddo
        istat(ilnum,13)=istat(ilnum,13)+1
        idarc(ilnum,istat(ilnum,13))=ipop 
      enddo

      return
      end subroutine NormalGen
!-----------------------------------------------------------------
!
!
! ORIGINAL: messy_airtraf_tools_ga_testfunction.f90
!
!
! *******************************************************
!     Mutliobjective Optimization
! *******************************************************
! *********** Zitzler-Deb-Thiele's (ZDT) Test Problems
! *********** 2-Objective Minimization
! *********** ZDT1 ...
! *********** ZDT2 ...
! *********** ZDT3 ...
! *********** ZDT4 ...
! *********** ZDT6 ...
! ********************************************************
! *********** Schaffer's (SCH) Test Problems
! *********** 2-Objective Minimization
! *********** SCH1 ... Range A is changeable (10 to 10e5)
! *********** SCH2 ...
! ********************************************************
! *********** Fonseca and Fleming's (FON) Test Problem
! *********** 2-Objective Minimization
! *********** FON1 ... # of DV is changeable 
! ********************************************************
! *********** Kursawe's (KUR) Test Problem
! *********** 2-Objective Minimization
! *********** KUR1 ...
! ********************************************************
! *********** Poloni et al. (POL) Test Problem
! *********** 2-Objective Minimization
! *********** POL1 ...
! ********************************************************
! *********** Viennet (VNT) Test Problem
! *********** 3-Objective Minimization
! *********** VNT1 ...
! ********************************************************
! *********** EX-MIN (EXN) Test Problems by Textbook
! *********** 2-Objective Minimization
! *********** EXN1 ...
! ********************************************************
! *********** EX-MAX (EXX) Test Problems by Textbook
! *********** 2-Objective Maximization
! *********** EXX1 ...
! ********************************************************
! *********** Test Problems at Rolls Royce 
! *********** 2-Objective Minimisation
! *********** 2-Design Variables 
! *********** RRT0 ... RRT0(2DVs):Convex Case
! *********** RRT1 ... FON1(2DVs):Convex Case(bd=4.)
! *********** RRT2 ... KUR1(3DVs):Discontinuous
! *********** RRT3 ... SCH1(2DVs):Concave Case
! *********** RRT4 ... FON1(20DVs):Convex Case(bd=1000.)
! *********** RRT5 ... VNT1(2DVs):3 Objectives
! *********** RRT6 ... ZDT1(30DVs):Convex Case
! *********** RRT7 ... ZDT2(30DVs):Concave Case
! *********** RRT8 ... ZDT3(30DVs):Discontinuous Case
! *********** RRT9 ... RRT9(2DVs):Multimodal Case
! *********** RRC1 ... CEX1(2DVs):Constrained Case
! *********** RRC2 ... TNK1(2DVs):Constrained Case
! *********** RRC3 ... WATER :Real Problem Case
! ******************** WRITTEN BY C PROGRAM (rrobj.c)
! ********************************************************
! *********** MT21 ... SCHF1(2OBJs(1),3DVs)
! *********** MT22 ... SCHF2(2DVs)
! *********** MT23 ... SCHF2(2DVs)
! *********** MT24 ... SCHF2(2DVs)
! *********** MT25 ... INT5(16DVs,integer)
! ********************************************************
! ********************************************************
! *********** Optimization for EXECUTABLE FILE
! *********** OUTE ... 
! *********** [REQURIED FILES]
! ***********   armoga.set  ... Settings of Optimisation 
! ***********   armoga.db   ... Design Information File
! ***********   obj.sh      ... Script for Evaluation
! ***********   cgeo.sh     ... Script for Geometry Constraints
! ***********                    (if icnrgeo=1)
! *********** [OPTIONAL] if iajinit=1
! ***********   armoga.ini  ... Settings of Optimisation 
! ********************************************************
! **************** OLD VERSION ***************************
! ***********   armoga.set  ... Settings of Optimisation 
! ***********   armoga.dv   ... Design Variables 
! ***********   armoga.obj  ... Objective Functions
! ***********   armoga.cstr ... Constraints 
! ***********   obj.sh      ... Script for Evaluation
! ***********   cgeo.sh     ... Script for Geometry Constraints
! ***********                    (if icnrgeo=1)
! *******************************************************
! *******************************************************
!     Constrained Mutliobjective Optimization
! *******************************************************
! *********** Constrained Test Problem (CTP)
! *********** 2-Objective Minimization with 1 Constraint
! *********** CTP2 ...
! *********** CTP6 ...
! *********** CTP7 ...
! ********************************************************
! *********** Constr-EX (CEX) Test Problems
! *********** 2-Objective Minimization with 2 Constraints
! *********** CEX1 ...
! ********************************************************
! *********** Tanaka (TNK) Test Problem
! *********** 2-Objective Minimization with 2 Constrains
! *********** TNK1 ...
! ********************************************************
! *********** Osyczka and Kundu (OSY) Test Problem
! *********** 2-Objective Minimization with 6 Constrains
! *********** OSY1 ...
! ********************************************************
! *********** Srinivas and Deb (SRN) Test Problem
! *********** 2-Objective Minimization with 2 Constrains
! *********** SRN1 ...
! ********************************************************
! *********** Ray (RAY) Test Problem: 'WATER'
! *********** 5-Objective Minimization with 7 Constrains
! *********** RAY1 ...
! ********************************************************
! ------------------------------------
      subroutine TestDesign
! ------------------------------------
!      USE fort_1,   ONLY:testp
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
!HY20140405-1      external dvall, dvexn1, dvosy1, dvray1
!HY20140405-1      external dvrrt0, dvrrt9
!HY20140405-1      external dvzdt4
!HY20140405-1      external debugdv
      REAL(DP) :: bd

!     phub ... maximum design parameter value constraint
!     phlb ... minimum design parameter value constraint
!     initial population distributes between
!     (pave+sd) and (pave-sd)
      bd   = 1000.0_dp   !20140619-1

!HY20140403-9      if(testp.eq.'ZDT1')call TestDV(-2,0,0,30,0.,1., dvall)
!HY      if(testp.eq.'ZDT2')call TestDV(-2,0,0,30,0.,1., dvall)
!HY      if(testp.eq.'ZDT3')call TestDV(-2,0,0,30,0.,1., dvall)
!HY      if(testp.eq.'ZDT4')call TestDV(-2,0,0,10,0.,0.,dvzdt4)
!HY      if(testp.eq.'ZDT6')call TestDV(-2,0,0,10,0.,1., dvall)

!HY      if(testcode.eq.'CTP')call TestDV(-2,-1,0,10,0.,0.,dvzdt4)

!HY      if(idebug.eq.2)call TestDV(-2,-ncon,0,ndv,0.,0.,debugdv)
 
!HY      if(testp.eq.'SCH1')call TestDV(-2,0,0,  1, -bd, bd, dvall)
!HY      if(testp.eq.'SCH2')call TestDV(-2,0,0,  1, -5.,10., dvall)
!HY      if(testp.eq.'FON1')call TestDV(-2,0,0,ndv, -4., 4., dvall)
!HY      if(testp.eq.'KUR1')call TestDV(-2,0,0,  3, -5., 5., dvall)
!HY      if(testp.eq.'POL1')call TestDV(-2,0,0,  2, -pi, pi, dvall)
!HY      if(testp.eq.'EXN1')call TestDV(-2,0,0,  2,  0., 0.,dvexn1)
!HY      if(testp.eq.'EXX1')call TestDV( 2,0,0,  2,  0., 0.,dvexn1)

!HY      if(testp.eq.'VNT1')call TestDV(-3,0,0,  2, -3., 3., dvall)
   
!HY      if(testp.eq.'CEX1')call TestDV(-2, 2,2, 2,  0., 0.,dvexn1)
!HY      if(testp.eq.'TNK1')call TestDV(-2, 2,0, 2,  0., pi, dvall)
!HY      if(testp.eq.'OSY1')call TestDV(-2, 6,4, 6,  0., 0.,dvosy1)
!HY      if(testp.eq.'SRN1')call TestDV(-2, 2,1, 2,-20.,20., dvall)
!      if(testp.eq.'RAY1')call TestDV(-5, 7,0, 3,  0., 0.,dvray1)
!HY      if(testp.eq.'RAY1') stop 'RRC3'

!HY      if(testp.eq.'RRT0')call TestDV(-2,0,0,2, 0.,0.,dvrrt0)
!HY      if(testp.eq.'RRT1')call TestDV(-2,0,0,2,-4.,4.,dvall)
!      if(testp.eq.'RRT1')call TestDV(-2,0,0,2,-160.,160.,dvall)
!HY      if(testp.eq.'RRT2')call TestDV(-2,0,0,3,-5.,5.,dvall)
!HY      if(testp.eq.'RRT3')call TestDV(-2,0,0,2,-4.,4.,dvall)
!      if(testp.eq.'RRT4')call TestDV(-2,0,0,2,-1000.,1000.,dvall)
!HY      if(testp.eq.'RRT4')call TestDV(-2,0,0,20,-4.,4.,dvall)
!HY      if(testp.eq.'RRT5')call TestDV(-3,0,0,2,-3.,3.,dvall)
!HY      if(testp.eq.'RRT6')call TestDV(-2,0,0,30, 0.,1.,dvall)
!HY      if(testp.eq.'RRT7')call TestDV(-2,0,0,30, 0.,1.,dvall)
!HY      if(testp.eq.'RRT8')call TestDV(-2,0,0,30, 0.,1.,dvall)
!HY      if(testp.eq.'RRT9')call TestDV(-2,0,0,2, 0.,0.,dvrrt9)
!HY      if(testp.eq.'RRC1')call TestDV(-2,2,2,2, 0.,0.,dvexn1) 
!HY      if(testp.eq.'RRC2')call TestDV(-2,2,0,2, 0.,pi, dvall)
!HY      if(testp.eq.'RRC3')call TestDV(-5,7,0,3, 0.,0.,dvray1)
!
!HY      if(testp.eq.'MT21')call TestDV(-2,0,0, 3,-10.,10.,dvall)
!HY      if(testp.eq.'MT22')call TestDV(-2,0,0, 1, -6., 6.,dvall)
!HY      if(testp.eq.'MT23')call TestDV(-2,0,0, 1, -2., 8.,dvall)
!HY      if(testp.eq.'MT24')call TestDV(-2,0,0, 2, -2., 2.,dvall)
!HY20140403-9      if(testp.eq.'MT25')call TestDV(-2,0,0,16,  0., 1.,dvall)
!      ngeocon = ncon
 
      return
      end subroutine TestDesign

! --------------------------------
      subroutine EvaPop(ist, ied, evap_nwaypoints, evap_p2_ac_routes_desc, evap_philon, evap_philat, evap_zgl_geopot_3d, &
                        evap_uwind_g, evap_vwind_g, evap_v_z_g, evap_t_scb_g, evap_cpc_g, evap_rho_air_dry_3d_g, & !HY20140420-1,20140522-3
                        evap_press_3d_g,evap_ATR20O3_3d_g,evap_ATR20CH4_3d_g,evap_ATR20H2O_3d_g,evap_ATR20CPC_3d_g,&
                        evap_ATR20CO2_3d_g) 
! --------------------------------
!### evaluation
!###  if indveva=1 then 
!       Evaluation will be performed for each individual
!      USE fort_1,   ONLY:testp
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6       USE messy_airtraf_tools_ga_testdirect,     ONLY: SystemExec
  IMPLICIT NONE
  INTEGER, INTENT(INOUT), OPTIONAL                          :: evap_nwaypoints  !20140721-6
  REAL(DP), DIMENSION(:),   INTENT(INOUT),OPTIONAL          :: evap_p2_ac_routes_desc  !HY20140420-1,20140422-3
  REAL(DP), DIMENSION(:),   INTENT(INOUT),OPTIONAL          :: evap_philon
  REAL(DP), DIMENSION(:),   INTENT(INOUT),OPTIONAL          :: evap_philat
  REAL(DP), DIMENSION(:,:,:), INTENT(INOUT),OPTIONAL        :: evap_zgl_geopot_3d   !global field 
  REAL(DP), DIMENSION(:,:,:), INTENT(INOUT),OPTIONAL        :: evap_uwind_g         !global field 
  REAL(DP), DIMENSION(:,:,:), INTENT(INOUT),OPTIONAL        :: evap_vwind_g         !global field
  REAL(DP), DIMENSION(:,:,:), INTENT(INOUT),OPTIONAL        :: evap_v_z_g           !global field 
  REAL(DP), DIMENSION(:,:,:), INTENT(INOUT),OPTIONAL        :: evap_t_scb_g         !global field  20140522-3 
  !Yin_20170423
  REAL(DP), DIMENSION(:,:,:), INTENT(INOUT),OPTIONAL        :: evap_cpc_g         !global field  20140522-3 
  !Yin_20170423
  REAL(DP), DIMENSION(:,:,:), INTENT(IN),OPTIONAL           :: evap_rho_air_dry_3d_g !global field  20160302-2,20170206-1 
  REAL(DP), DIMENSION(:,:,:), INTENT(IN),OPTIONAL           :: evap_press_3d_g      !global field  20170215-4 
  !Yin_20170801+
  REAL(DP), DIMENSION(:,:,:), INTENT(IN),OPTIONAL           :: evap_ATR20O3_3d_g      !global field  20170801 
  REAL(DP), DIMENSION(:,:,:), INTENT(IN),OPTIONAL           :: evap_ATR20CH4_3d_g      !global field  20170801 
  REAL(DP), DIMENSION(:,:,:), INTENT(IN),OPTIONAL           :: evap_ATR20H2O_3d_g      !global field  20170801 
  REAL(DP), DIMENSION(:,:,:), INTENT(IN),OPTIONAL           :: evap_ATR20CPC_3d_g      !global field  20170801 
  REAL(DP), DIMENSION(:,:,:), INTENT(IN),OPTIONAL           :: evap_ATR20CO2_3d_g      !global field  20170801 
  !Yin_20170801-
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                   :: temp_evap_rho_air_dry_3d_g !global field  20170206-1 
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                   :: temp_evap_press_3d_g !global field  20170215-4
  !Yin_20170801+
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                   :: temp_evap_ATR20O3_3d_g !global field  20170801
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                   :: temp_evap_ATR20CH4_3d_g !global field  20170801
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                   :: temp_evap_ATR20H2O_3d_g !global field  20170801
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                   :: temp_evap_ATR20CPC_3d_g !global field  20170801
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                   :: temp_evap_ATR20CO2_3d_g !global field  20170801
  !Yin_20170801-
  REAL(DP) :: x(mdv), objf(mobj)
  REAL(DP) :: xm(mpop*mdv), objfm(mpop*mobj),gvm(mpop*mcon)
  INTEGER :: ist,ied,nsum,icn,iof,idv,ipop,nsum1,nsum2,nsum3,itype,nump
!HY20140405-1      external fzdt1, fzdt6, gzdt1, gzdt4, gzdt6, hzdt1, hzdt2, hzdt3
!HY20140405-1      external fsch1, fsch2, ffon1, fkur1, fpol1, fexn1, fexx1
!HY20140405-1      external fvnt1, ftnk1, fosy1, fsrn1
 
      IF(PRESENT(evap_rho_air_dry_3d_g))THEN   !20170206-1 for fuel_opt, NOx_opt, H2O_opt, Cost-opt, COC_opt.
         ALLOCATE(temp_evap_rho_air_dry_3d_g(nlon_ga, nlev_ga, ngl_ga))
         temp_evap_rho_air_dry_3d_g=evap_rho_air_dry_3d_g
!         write(*,*)"allocation3:nlon_ga,nlev_ga,ngl_ga",nlon_ga,nlev_ga,ngl_ga
      ENDIF
      IF(PRESENT(evap_press_3d_g))THEN   !20170215-4 for NOx_opt
         ALLOCATE(temp_evap_press_3d_g(nlon_ga, nlev_ga, ngl_ga))
         temp_evap_press_3d_g=evap_press_3d_g
      ENDIF
!Yin_20170801+
      IF(PRESENT(evap_ATR20O3_3d_g))THEN   !20170801 for ATR20_opt
         ALLOCATE(temp_evap_ATR20O3_3d_g(nlon_ga, nlev_ga, ngl_ga))
         temp_evap_ATR20O3_3d_g=evap_ATR20O3_3d_g
      ENDIF
      IF(PRESENT(evap_ATR20CH4_3d_g))THEN   !20170801 for ATR20_opt
         ALLOCATE(temp_evap_ATR20CH4_3d_g(nlon_ga, nlev_ga, ngl_ga))
         temp_evap_ATR20CH4_3d_g=evap_ATR20CH4_3d_g
      ENDIF
      IF(PRESENT(evap_ATR20H2O_3d_g))THEN   !20170801 for ATR20_opt
         ALLOCATE(temp_evap_ATR20H2O_3d_g(nlon_ga, nlev_ga, ngl_ga))
         temp_evap_ATR20H2O_3d_g=evap_ATR20H2O_3d_g
      ENDIF
      IF(PRESENT(evap_ATR20CPC_3d_g))THEN   !20170801 for ATR20_opt
         ALLOCATE(temp_evap_ATR20CPC_3d_g(nlon_ga, nlev_ga, ngl_ga))
         temp_evap_ATR20CPC_3d_g=evap_ATR20CPC_3d_g
      ENDIF
      IF(PRESENT(evap_ATR20CO2_3d_g))THEN   !20170801 for ATR20_opt
         ALLOCATE(temp_evap_ATR20CO2_3d_g(nlon_ga, nlev_ga, ngl_ga))
         temp_evap_ATR20CO2_3d_g=evap_ATR20CO2_3d_g
      ENDIF

!Yin_20170801-
      if(testp.eq.'OUTE')then
        itype = 2
        nump  = ied-ist+1
        nsum1 = 1
        nsum2 = 1
        nsum3 = 1
        do ipop=ist,ied
          do idv=1,ndv
            xm(nsum1) = ptall(ipop,idv)
            nsum1 = nsum1 + 1
          enddo
          do iof=1,nobj
            objfm(nsum2) = objinit
            nsum2 = nsum2 + 1
          enddo
          if(ncon.gt.0)then
            do icn=1,ncon
              gvm(nsum3) = gvalinit
              nsum3 = nsum3 + 1
            enddo
          endif
        enddo
        if((iman.eq.1).and.(indveva.eq.0))itype=0 
!C1.2.0        call SystemExec(xm,ndv,objfm,nobj,gvm,ncon,igen,nump,itype)
!HY20140420-1  call SystemExec(xm,ndv,objfm,nobj,gvm,ncon,iter,nump,itype)
        SELECT CASE(ga_option_traj_calc_g)   !20160302-2
        CASE(1,5)  !Wind_opt, CPC_opt Yin_20170423, HY 20170609-1
        call SystemExec(xm,ndv,objfm,nobj,gvm,ncon,iter,nump,itype, &
                        evap_nwaypoints, evap_p2_ac_routes_desc, evap_philon, evap_philat, evap_zgl_geopot_3d, &
                        evap_uwind_g, evap_vwind_g, evap_v_z_g, evap_t_scb_g, evap_cpc_g)  !HY20140420-1,20140522-3,20140721-6 
        CASE(2,4,6,7,10)  !Fuel_opt, H2O_opt, Cost_opt or COC_opt   !20170222-2,20170224-2,20170316-2
        !IF(PRESENT(evap_rho_air_dry_3d_g))temp_evap_rho_air_dry_3d_g=evap_rho_air_dry_3d_g   !20170206-1
        call SystemExec(xm,ndv,objfm,nobj,gvm,ncon,iter,nump,itype, &
                        evap_nwaypoints, evap_p2_ac_routes_desc, evap_philon, evap_philat, evap_zgl_geopot_3d, &
                        evap_uwind_g, evap_vwind_g, evap_v_z_g, evap_t_scb_g, evap_cpc_g, & !HY20140420-1,20140522-3,20140721-6 
                        temp_evap_rho_air_dry_3d_g)                             !20160302-2,20170206-1
!        write(*,*)'check_ga_option_traj_calc_4', ga_option_traj_calc_g          !20160304-1
        CASE(3)  !NOx_opt
        !IF(PRESENT(evap_rho_air_dry_3d_g))temp_evap_rho_air_dry_3d_g=evap_rho_air_dry_3d_g   !20170206-1
        call SystemExec(xm,ndv,objfm,nobj,gvm,ncon,iter,nump,itype, &
                        evap_nwaypoints, evap_p2_ac_routes_desc, evap_philon, evap_philat, evap_zgl_geopot_3d, &
                        evap_uwind_g, evap_vwind_g, evap_v_z_g, evap_t_scb_g, evap_cpc_g, & !HY20140420-1,20140522-3,20140721-6 
                        temp_evap_rho_air_dry_3d_g, temp_evap_press_3d_g)       !20160302-2,20170206-1,20170215-4
!        write(*,*)'check_ga_option_traj_calc_4', ga_option_traj_calc_g          !20160304-1
        !Yin_20170801+
        CASE(8,9)  !ATR20_opt
        call SystemExec(xm,ndv,objfm,nobj,gvm,ncon,iter,nump,itype, &
                        evap_nwaypoints, evap_p2_ac_routes_desc, evap_philon, evap_philat, evap_zgl_geopot_3d, &
                        evap_uwind_g, evap_vwind_g, evap_v_z_g, evap_t_scb_g, evap_cpc_g,&
                        temp_evap_rho_air_dry_3d_g,temp_evap_press_3d_g,&
                        temp_evap_ATR20O3_3d_g,temp_evap_ATR20CH4_3d_g,temp_evap_ATR20H2O_3d_g,&
                        temp_evap_ATR20CPC_3d_g,temp_evap_ATR20CO2_3d_g)  
        !Yin_20170801-

        CASE DEFAULT
           write(*,*)'ERROR:ga_option_traj_calc on SystemExec',ga_option_traj_calc_g
        ENDSELECT
        nsum = 1
!        write(*,*)"ist,ied",ist,ied
        do ipop=ist,ied
          do iof=1,nobj
            fvall(ipop,iof) = objfm(nsum)
            nsum = nsum + 1
          enddo
          ipinfo(ipop,4) = 1
        enddo
        if(ncon.gt.0)then
!          call SystemCstrAll(gvm,ncon,nump)
          nsum = 1
          do ipop=ist,ied
            do icn=1,ncon 
              gvall(ipop,icn) = gvm(nsum)
              nsum = nsum + 1
            enddo
          enddo
        endif
        return
      endif

      do ipop=ist,ied
        do idv=1,ndv
          x(idv) = ptall(ipop,idv)
        enddo
!HY20140403-9        if(testp.eq.'ZDT1') &
!HY             call Tunable2(objf,nobj,x,ndv,fzdt1,gzdt1,hzdt1)
!HY        if(testp.eq.'ZDT2') &
!HY             call Tunable2(objf,nobj,x,ndv,fzdt1,gzdt1,hzdt2)
!HY        if(testp.eq.'ZDT3') &
!HY             call Tunable2(objf,nobj,x,ndv,fzdt1,gzdt1,hzdt3)
!HY        if((testp.eq.'ZDT4').or.(testcode.eq.'CTP')) &
!HY             call Tunable2(objf,nobj,x,ndv,fzdt1,gzdt4,hzdt1)
!HY        if(testp.eq.'ZDT6') &
!HY             call Tunable2(objf,nobj,x,ndv,fzdt6,gzdt6,hzdt2)
!HY        if(testp.eq.'SCH1') &
!HY             call Multiobj(objf,nobj,x,ndv,fsch1)
!HY        if(testp.eq.'SCH2') &
!HY             call Multiobj(objf,nobj,x,ndv,fsch2)
!HY        if(testp.eq.'FON1') &
!HY             call Multiobj(objf,nobj,x,ndv,ffon1)
!HY        if(testp.eq.'KUR1') &
!HY             call Multiobj(objf,nobj,x,ndv,fkur1)
!HY        if(testp.eq.'POL1') &
!HY             call Multiobj(objf,nobj,x,ndv,fpol1)
!HY        if(testp.eq.'EXN1') &
!HY             call Multiobj(objf,nobj,x,ndv,fexn1)
!HY        if(testp.eq.'EXX1') &
!HY             call Multiobj(objf,nobj,x,ndv,fexx1)
!HY        if(testp.eq.'VNT1') &
!HY             call Multiobj(objf,nobj,x,ndv,fvnt1)
!HY        if(testp.eq.'CEX1') &
!HY             call Multiobj(objf,nobj,x,ndv,fexn1)
!HY        if(testp.eq.'TNK1') &
!HY             call Multiobj(objf,nobj,x,ndv,ftnk1)
!HY        if(testp.eq.'OSY1') &
!HY             call Multiobj(objf,nobj,x,ndv,fosy1)
!HY        if(testp.eq.'SRN1') &
!HY20140403-9             call Multiobj(objf,nobj,x,ndv,fsrn1)
!HY20140327-3     if((testcode.eq.'RRT').or.(testcode.eq.'RRC')
!HY20140327-3  &       .or.(testcode.eq.'MT2'))
!HY20140327-3  &        call DirectCF(objf,nobj,x,ndv,testp)

        do iof=1,nobj
!HY20140404-1          fvall(ipop,iof) = objf(iof)
          ipinfo(ipop,4) = 1
        enddo
      enddo
      !for fuel_opt, NOx_opt, H2O_opt, Cost_opt, COC_opt 
      IF(ALLOCATED(temp_evap_rho_air_dry_3d_g))DEALLOCATE(temp_evap_rho_air_dry_3d_g)   !20170206-1,20170212-3
      !for NO_opt  
      IF(ALLOCATED(temp_evap_press_3d_g))DEALLOCATE(temp_evap_press_3d_g)               !20170215-4
!Yin_20170801+
      IF(ALLOCATED(temp_evap_ATR20O3_3d_g))DEALLOCATE(temp_evap_ATR20O3_3d_g)               !20170801
      IF(ALLOCATED(temp_evap_ATR20CH4_3d_g))DEALLOCATE(temp_evap_ATR20CH4_3d_g)               !20170801
      IF(ALLOCATED(temp_evap_ATR20H2O_3d_g))DEALLOCATE(temp_evap_ATR20H2O_3d_g)               !20170801
      IF(ALLOCATED(temp_evap_ATR20CPC_3d_g))DEALLOCATE(temp_evap_ATR20CPC_3d_g)               !20170801
      IF(ALLOCATED(temp_evap_ATR20CO2_3d_g))DEALLOCATE(temp_evap_ATR20CO2_3d_g)               !20170801
!Yin_20170801-
      return
      end subroutine EvaPop

! -----------------------------------------
      subroutine GvalPop(ist,ied,ngst,nged)
! -----------------------------------------
!### evaluation of constraint value (Gvalue)
!###  if indvcnr=1 then
!       Evaluation will be performed for each individual
!      USE fort_1,   ONLY:testp
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6       USE messy_airtraf_tools_ga_testdirect,     ONLY: SystemExec
      IMPLICIT NONE
      REAL(DP) :: x(mdv), f(mobj), gv(mcon)
      REAL(DP) :: xm(mpop*mdv), objfm(mpop*mobj),gvm(mpop*mcon)
      INTEGER :: ngst,nged,ist,ied,nsum,icn,iof,ipop,idv,nsum1,nsum2,nsum3,nump
!HY20140405-1      external ccex1, ctnk1, cosy1, csrn1

      if((ngst.gt.nged).or.(nged.lt.1))then
!        write(*,*)ngst,nged,testp
        write(mtout,*)'GvalPop has problem.'
        return
      endif

      if(testp.eq.'OUTE')then
        nump  = ied-ist+1
        nsum1 = 1
        nsum2 = 1
        nsum3 = 1
        do ipop=ist,ied
          do idv=1,ndv
            xm(nsum1) = ptall(ipop,idv)
            nsum1 = nsum1 + 1
          enddo
          do iof=1,nobj
            objfm(nsum2) = objinit
            nsum2 = nsum2 + 1
          enddo
          do icn=1,ncon
            gvm(nsum3) = gvalinit
            nsum3 = nsum3 + 1
          enddo
        enddo
!C1.2.0        call SystemExec(xm,ndv,objfm,nobj,gvm,ncon,igen,nump,1)
        call SystemExec(xm,ndv,objfm,nobj,gvm,ncon,iter,nump,1)
        nsum = 1
        do ipop=ist,ied
          do icn=ngst,nged
            gvall(ipop,icn) = gvm(nsum)
            nsum = nsum + 1
          enddo
        enddo
        return
      endif


      do ipop=ist,ied
        do idv=1,ndv
          x(idv) = ptall(ipop,idv)
        enddo
        do iof=1,nobj
          f(iof) = fvall(ipop,iof)
        enddo

!HY20140403-7        if(testp.eq.'CEX1') &
!HY             call Gcalc(gv,ncon,x,ndv,ccex1)
!HY        if(testp.eq.'TNK1') &
!HY             call Gcalc(gv,ncon,x,ndv,ctnk1)
!HY        if(testp.eq.'OSY1') &
!HY             call Gcalc(gv,ncon,x,ndv,cosy1)
!HY        if(testp.eq.'SRN1') &
!HY20140403-7             call Gcalc(gv,ncon,x,ndv,csrn1)
        if(testp.eq.'CTP2') &
         call Ctunable6(gv,ncon,f,nobj,x,ndv,-0.2_dp*pi,0.2_dp,10.0_dp,1,6,1.0_dp)
        if(testp.eq.'CTP6') &
         call Ctunable6(gv,ncon,f,nobj,x,ndv,0.1_dp*pi,40.0_dp,0.5_dp,1,2,-2.0_dp)
        if(testp.eq.'CTP7') &
         call Ctunable6(gv,ncon,f,nobj,x,ndv,-0.05_dp*pi,40.0_dp,5.0_dp,1,6,0.0_dp)

!HY20140327-3     if(testcode.eq.'RRC')
!HY20140327-3  &        call DirectCG(gv,nged,f,nobj,x,ndv,testp)

        do icn=ngst,nged
          gvall(ipop,icn) = gv(icn)
        enddo
       enddo

       return
       end subroutine GvalPop

! ------------------------------------
      subroutine DesignMAXMIN(ilnum)
! ------------------------------------
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      INTEGER :: ilnum,icn,iof
      objinit  = (-1.0_dp) * dble(imaxmin) * 1.0d5 !HY20140414-5 
      gvalinit = 1.0d5                            !HY20140414-5

      do iof=1,nobj
        objmax(ilnum,iof)    = (-1.0_dp) * abs(objinit)
        objmin(ilnum,iof)    =            abs(objinit)
        objallmax(ilnum,iof) = objmax(ilnum,iof)
        objallmin(ilnum,iof) = objmin(ilnum,iof)
      enddo
      if(ncon.gt.0)then
        do icn=1,ncon
          gvalmax(ilnum,icn) = (-1.0_dp)*gvalinit
          gvalmin(ilnum,icn) =       gvalinit
          gvalallmax(ilnum,icn) = gvalmax(ilnum,icn)
          gvalallmin(ilnum,icn) = gvalmin(ilnum,icn)
        enddo
      endif

      return
      end subroutine DesignMAXMIN

! --------------------------
      subroutine AdjustDV(ilnum,itype)
! --------------------------
      IMPLICIT NONE
      INTEGER :: ilnum,itype
      return
      end subroutine AdjustDV

! ---------------------------------------------------------
      subroutine TestDV(numobj,numcon,nlcon,numdv,bl,bu,dv)
! ---------------------------------------------------------
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      REAL(DP), EXTERNAL :: dv
      REAL(DP) :: bu, bl
      INTEGER  :: nlcon,numdv,numobj,numcon,idv,ilnum
      if(numobj.lt.0)then
        imaxmin = -1
      else
        imaxmin =  1
      endif
      nobj    = abs(numobj)

      if(numcon.lt.0)then
        igwtobj = 1
      else
        igwtobj = 0
      endif
      ncon    = abs(numcon)
      
      ngeocon = ncon
      nlinear = nlcon
      ndv     = numdv
      if(icnrgeo.eq.0)ngeocon=0

      do idv=1,ndv 
        phub(idv) = dv(idv, bu,  1) 
        phlb(idv) = dv(idv, bl, -1)
      enddo


!  ********* IMOGA **************  
      do ilnum=1,nisland
        imoga(ilnum,1)=1
        imoga(ilnum,2)=nobj
        imoga(ilnum,3)=1
        imoga(ilnum,4)=nobj
        imoga(ilnum,5)=1
        imoga(ilnum,6)=ncon
        imoga(ilnum,7)=1
        imoga(ilnum,8)=ngeocon
        if(ncon.eq.0)imoga(ilnum,5)=0
        if(ngeocon.eq.0)imoga(ilnum,7)=0
      enddo

      return
      end subroutine TestDV

! --------------------------------------------
      subroutine Tunable2(f12,nobj,x,ndv,f,g,h)
! --------------------------------------------
      IMPLICIT NONE
      INTEGER  :: nobj,ndv
      REAL(DP), EXTERNAL :: f, g, h
      REAL(DP) :: fx, gx, hfg
      REAL(DP) :: f12(nobj), x(ndv)
      fx  = f(x(1))
      gx  = g(x,ndv)
      hfg = h(fx,gx)

      f12(1) = fx
      f12(2) = gx * hfg

!     if(f12(1).lt.bv)f12(1)=0.
      
      return 
      end subroutine Tunable2

! --------------------------------------------
      subroutine Multiobj(fmobj,nobj,x,ndv,f)
! --------------------------------------------
      IMPLICIT NONE
      INTEGER :: nobj,ndv,iof 
      REAL(DP) :: fmobj(nobj), x(ndv)
      REAL(DP), EXTERNAL :: f

! ---- check ----
!      write(*,*)(x(ndv),ndv=1,5)
!      stop
! --------------
      do iof=1,nobj
        fmobj(iof) = f(x,ndv,iof)
      enddo

! ---- check ----
!      write(*,*)(fmobj(iof),iof=1,5)
!      stop
! --------------
      return 
      end subroutine Multiobj


! --------------------------------------------
      subroutine Gcalc(gv,ncon,x,ndv,c)
! --------------------------------------------
      IMPLICIT NONE
      INTEGER :: ncon,ndv,icn
      REAL(DP) :: gv(ncon), x(ndv)
      REAL(DP), EXTERNAL :: c
      do icn=1,ncon
        gv(icn) = c(x,ndv,icn)
      enddo

      return 
      end subroutine Gcalc

! ---------------------------------------------------------------
      subroutine Ctunable6(gv,ncon,f,nobj,x,ndv,theta,a,b,ic,id,e)
! ---------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ic,id,ncon,nobj,ndv
      REAL(DP) :: gv(ncon), f(nobj), x(ndv)
      REAL(DP) :: pi, f1, f2, theta, a, b, e
      REAL(DP) :: eqlhs, ein, eqrhs, cx  
      pi = 4.0_dp*atan(1.0_dp)   !20140619-1

      f1 = f(1)
      f2 = f(2)

      eqlhs = cos(theta)*(f2-e) - sin(theta)*f1
      ein   = sin(theta)*(f2-e) + cos(theta)*f1
      eqrhs = a * (abs( sin(b * pi * ein**ic) ))**id
      cx    = eqlhs - eqrhs
      gv(1) = (-1.0_dp) * cx   !20140619-1

      return 
      end subroutine Ctunable6

! -----------------------------------
      function dvall(idv, bd, idir )
! -----------------------------------
      IMPLICIT NONE
      REAL(DP) :: dvall
      REAL(DP) :: bd
      INTEGER  :: idv, idir
      dvall = bd

      return
      end function dvall
       
! -----------------------------------
      function dvzdt4(idv, tmp, idir)
! -----------------------------------
      IMPLICIT NONE
      REAL(DP) :: dvzdt4 
      REAL(DP) :: tmp
      INTEGER  :: idv, idir
      if(idir.gt.0)then
        dvzdt4 =  5.0_dp   !20140619-1
        if(idv.eq.1)dvzdt4 = 1.0_dp   !20140619-1
      else
        dvzdt4 = -5.0_dp   !20140619-1
        if(idv.eq.1)dvzdt4 = 0.0_dp   !20140619-1
      endif

      return
      end function dvzdt4

! -----------------------------------
      function dvexn1(idv, tmp, idir)
! -----------------------------------
      IMPLICIT NONE
      REAL(DP) :: dvexn1 
      REAL(DP) :: tmp
      INTEGER  :: idv, idir
      if(idir.gt.0)then
        dvexn1 = 5.0_dp   !20140619-1
        if(idv.eq.1)dvexn1 = 1.0_dp   !20140619-1
      else
        dvexn1 = 0.0_dp   !20140619-1
        if(idv.eq.1)dvexn1 = 0.1_dp   !20140619-1
      endif

      return
      end function dvexn1

! -----------------------------------
      function dvosy1(idv, tmp, idir)
! -----------------------------------
      IMPLICIT NONE
      REAL(DP) :: dvosy1 
      REAL(DP) :: tmp
      INTEGER  :: idv, idir
      if(idir.gt.0)then
        dvosy1 = 10.0_dp   !20140619-2
        if((idv.eq.3).or.(idv.eq.5))dvosy1 = 5.0_dp   !20140619-2
        if(idv.eq.4)dvosy1 = 6.0_dp   !20140619-2
      else
        dvosy1 = 0.0_dp   !20140619-2
        if((idv.eq.3).or.(idv.eq.5))dvosy1 = 1.0_dp   !20140619-2
      endif

      return
      end function dvosy1

! -----------------------------------
      function dvray1(idv, tmp, idir)
! -----------------------------------
      IMPLICIT NONE
      REAL(DP) :: dvray1 
      REAL(DP) :: tmp
      INTEGER  :: idv, idir
      if(idir.gt.0)then
        dvray1=0.1_dp   !20140619-2
        if(idv.eq.1)dvray1=0.45_dp   !20140619-2
      else
        dvray1 = 0.01_dp   !20140619-2
      endif

      return
      end function dvray1

! -----------------------------------
      function dvrrt0(idv, tmp, idir)
! -----------------------------------
      IMPLICIT NONE
      REAL(DP) :: dvrrt0 
      REAL(DP) :: tmp
      INTEGER  :: idv, idir
      if(idir.gt.0)then
        dvrrt0 = 1.0_dp   !20140619-2
      else
        dvrrt0 = -1.0_dp   !20140619-2
        if(idv.eq.1)dvrrt0=0.0_dp   !20140619-2
      endif

      return
      end function dvrrt0

! -----------------------------------
      function dvrrt9(idv, tmp, idir)
! -----------------------------------
      IMPLICIT NONE
      REAL(DP) :: dvrrt9 
      REAL(DP) :: tmp
      INTEGER  :: idv, idir
      if(idir.gt.0)then
        dvrrt9 = 1.0_dp   !20140619-2
      else
        dvrrt9 = 0.0_dp   !20140619-2
        if(idv.eq.1)dvrrt9=0.1_dp   !20140619-2
      endif

      return
      end function dvrrt9

! -----------------------------------
      function fzdt1(x1)
! -----------------------------------
      IMPLICIT NONE
      REAL(DP) :: fzdt1 
      REAL(DP) :: x1
      fzdt1 = x1
      return
      end function fzdt1

! -----------------------------------
      function fzdt6(x1)
! -----------------------------------
      IMPLICIT NONE
      REAL(DP) :: fzdt6 
      REAL(DP) :: pi, x1
      pi = 4.0_dp*atan(1.0_dp)   !20140619-2
      fzdt6 = 1.0_dp - exp(-4.0_dp*x1) * (sin(6.0_dp*pi*x1))**6    !20140619-2
      return
      end function fzdt6

! ----------------------------
      function gzdt1(x,ndv)
! ----------------------------
      IMPLICIT NONE
      REAL(DP) :: gzdt1 
      INTEGER :: ndv,idv
      REAL(DP) :: x(ndv)
      gzdt1 = 0.0_dp   !20140619-2
      do idv=2,ndv
        gzdt1 = gzdt1 + x(idv)
      enddo
      gzdt1 = 1.0_dp + 9.0_dp/DBLE(ndv-1)*gzdt1   !HY20140414-5,20140619-2

      return
      end function gzdt1

! ----------------------------
      function gzdt4(x,ndv)
! ----------------------------
      IMPLICIT NONE
      INTEGER  :: ndv,idv
      REAL(DP) :: gzdt4 
      REAL(DP) :: x(ndv)
      REAL(DP) :: pi
      pi = 4.0_dp*atan(1.0_dp)   !20140619-2

      gzdt4 = 0.0_dp   !20140619-2
      do idv=2,ndv
        gzdt4 = gzdt4 + (x(idv)**2 -10.0_dp*cos(4.0_dp*pi*x(idv)))   !20140619-2
      enddo
      gzdt4 = 1.0_dp + 10.0_dp*DBLE(ndv-1)+gzdt4   !HY20140414-5,20140619-2

      return
      end function gzdt4

! ----------------------------
      function gzdt6(x,ndv)
! ----------------------------
      IMPLICIT NONE
      INTEGER :: ndv,idv
      REAL(DP) :: gzdt6
      REAL(DP) :: x(ndv)
      gzdt6 = 0.0_dp   !20140619-2
      do idv=2,ndv
        gzdt6 = gzdt6 + x(idv)
      enddo
      gzdt6 = 1.0_dp + 9.0_dp*(gzdt6/9.0_dp)**0.25   !20140619-2

      return
      end function gzdt6

! ----------------------------
      function hzdt1(f1,g)
! ----------------------------
      IMPLICIT NONE
      REAL(DP) :: hzdt1 
      REAL(DP) :: g, f1
      hzdt1 = 1.0_dp - sqrt(f1/g)   !20140619-2
      return
      end function hzdt1

! ----------------------------
      function hzdt2(f1,g)
! ----------------------------
      IMPLICIT NONE
      REAL(DP) :: hzdt2 
      REAL(DP) :: g, f1
      hzdt2 = 1.0_dp - (f1/g)**2   !20140619-2
      return
      end function hzdt2

! ----------------------------
      function hzdt3(f1,g)
! ----------------------------
      IMPLICIT NONE
      REAL(DP) :: hzdt3 
      REAL(DP) :: pi,g, f1 
      pi = 4.0_dp*atan(1.0_dp)   !20140619-2
      hzdt3 = 1.0_dp - sqrt(f1/g) - (f1/g)*sin(10.0_dp*pi*f1)   !20140619-2
      return
      end function hzdt3

! -----------------------------------
      function debugdv(idv, tmp, idir)
! -----------------------------------
      IMPLICIT NONE
      REAL(DP) :: debugdv 
      REAL(DP) :: tmp
      INTEGER  :: idv, idir
      if(idir.gt.0)then
        debugdv =   0.10_dp      !20140619-2
        if(idv.eq.1)debugdv = 1.0_dp   !20140619-2
      else
        debugdv =  -0.10_dp   !20140619-2
        if(idv.eq.1)debugdv = 0.0_dp   !20140619-2
      endif

      return
      end function debugdv

! ----------------------------
      function fsch1(x,ndv,iof)
! ----------------------------
      IMPLICIT NONE
      INTEGER :: iof, ndv,i
      REAL(DP) :: fsch1 
      REAL(DP) :: x(ndv)
      fsch1 = 0.0_dp   !20140619-2
      if(iof.eq.1)then
        do i=1,ndv
          fsch1 = fsch1 + x(i)**2
        enddo
      else
        do i=1,ndv
          fsch1 = fsch1 + (x(i)-2.0_dp)**2   !20140619-2
        enddo
      endif
      fsch1 = fsch1 / DBLE(ndv)   !HY20140414-5

      return
      end function fsch1

! ----------------------------
      function fsch2(x,ndv,iof)
! ----------------------------
      IMPLICIT NONE
      INTEGER  :: iof, ndv
      REAL(DP) :: fsch2 
      REAL(DP) :: x(ndv)
      REAL(DP) :: x1 
      x1 = x(1)

      if(iof.eq.1)then
        if(x1.le.1.0_dp)then   !20140619-3
          fsch2 = -1.0_dp * x1   !20140619-3
        elseif(x1.le.3.0_dp)then   !20140619-3
          fsch2 = x1 - 2.0_dp   !20140619-3
        elseif(x1.le.4.0_dp)then   !20140619-3
          fsch2 = 4.0_dp - x1   !20140619-3
        else
          fsch2 = x1 - 4.0_dp   !20140619-3
        endif
      else
        fsch2 = (x1-5.0_dp)**2   !20140619-3
      endif

      return
      end function fsch2

! ------------------------------
      function ffon1(x,ndv,iof)
! ------------------------------
      IMPLICIT NONE
      INTEGER  :: iof, ndv,i
      REAL(DP) :: ffon1 
      REAL(DP) :: x(ndv)
      REAL(DP) :: fsum
      fsum = 0.0_dp   !20140619-3

      if(iof.eq.1)then
        do i=1,ndv
          fsum = fsum + (x(i)-1.0_dp/sqrt(DBLE(ndv)))**2  !HY20140414-5,20140619-3
        enddo
        ffon1 = 1.0_dp - exp(-1.0_dp*fsum)   !20140619-3
      else
        do i=1,ndv
          fsum = fsum + (x(i)+1.0_dp/sqrt(DBLE(ndv)))**2  !HY20140414-5,20140619-3
        enddo
        ffon1 = 1.0_dp - exp(-1.0_dp*fsum)   !20140619-3
      endif

      return
      end function ffon1

! ------------------------------
      function fkur1(x,ndv,iof)
! ------------------------------
      IMPLICIT NONE
      INTEGER  :: iof, ndv,i
      REAL(DP) :: fkur1 
      REAL(DP) :: x(ndv)
      REAL(DP) :: fsum
      fsum = 0.0_dp   !20140619-3

      if(iof.eq.1)then
        do i=1,2
          fsum = fsum + & 
            (-10.0_dp*exp(-0.2_dp*sqrt(x(i)**2+x(i+1)**2)))   !20140619-3
        enddo
        fkur1 = fsum
      else
        do i=1,ndv
          fsum = fsum + (abs(x(i))**0.8 + 5.0_dp*sin(x(i)**3))   !20140619-3
        enddo
        fkur1 = fsum 
      endif

      return
      end function fkur1

! ------------------------------
      function fpol1(x,ndv,iof)
! ------------------------------
      IMPLICIT NONE
      INTEGER  :: iof, ndv
      REAL(DP) :: fpol1 
      REAL(DP) :: x(ndv)
      REAL(DP) :: a1, a2, b1, b2,x1,x2
      x1 = x(1)
      x2 = x(2)

      if(iof.eq.1)then
        a1 = 0.5_dp*sin(1.0_dp) - 2.0_dp*cos(1.0_dp) +        sin(2.0_dp) -1.5_dp*cos(2.0_dp)   !20140619-3
        a2 = 1.5_dp*sin(1.0_dp) -        cos(1.0_dp) + 2.0_dp*sin(2.0_dp) -0.5_dp*cos(2.0_dp)   !20140619-3
        b1 = 0.5_dp*sin(x1)     -        cos(x1)     +        sin(x2)     -1.5_dp*cos(x2)       !20140619-3
        b2 = 1.5_dp*sin(x1)     -        cos(x1)     + 2.0_dp*sin(x2)     -0.5_dp*cos(x2)       !20140619-3
        fpol1 = 1.0_dp + (a1-b1)**2 + (a2-b2)**2   !20140619-3
      else
        fpol1 = (x1+3.0_dp)**2 + (x2+1.0_dp)**2    !20140619-3
      endif

      return
      end function fpol1

! ------------------------------
      function fvnt1(x,ndv,iof)
! ------------------------------
      IMPLICIT NONE
      INTEGER  :: iof, ndv
      REAL(DP) :: fvnt1 
      REAL(DP) :: x(ndv)
      REAL(DP) :: x1,x2
      x1 = x(1)
      x2 = x(2)

      if(iof.eq.1)then
        fvnt1 = 0.5_dp * (x1**2 + x2**2) + sin(x1**2 + x2**2)   !20140619-3
      elseif(iof.eq.2)then
        fvnt1 = (3.0_dp*x1-2.0_dp*x2+4.0_dp)**2/8.0_dp + (x1-x2+1.0_dp)**2/27.0_dp + 15.0_dp   !20140619-3
      else
        fvnt1 = 1.0_dp/(x1**2+x2**2+1.0_dp)-1.1_dp*exp(-1.0_dp*(x1**2+x2**2))   !20140619-3
      endif

      return
      end function fvnt1

! ------------------------------
      function fexn1(x,ndv,iof)
! ------------------------------
      IMPLICIT NONE
      INTEGER  :: iof, ndv
      REAL(DP) :: fexn1 
      REAL(DP) :: x(ndv)
      REAL(DP) :: x1,x2
      x1 = x(1)
      x2 = x(2)

      if(iof.eq.1)then
        fexn1 = x1
      else
        fexn1 = (1.0_dp+x2)/x1   !20140619-3
      endif

      return
      end function fexn1

! ------------------------------
      function fexx1(x,ndv,iof)
! ------------------------------
      IMPLICIT NONE
      INTEGER  :: iof, ndv
      REAL(DP) :: fexx1 
      REAL(DP) :: x(ndv)
      REAL(DP) :: x1,x2
      x1 = x(1)
      x2 = x(2)

      if(iof.eq.1)then
        fexx1 = 1.1_dp - x1   !20140619-3
      else
        fexx1 = 60.0_dp - (1.0_dp+x2)/x1   !20140619-3
      endif

      return
      end function fexx1

! ------------------------------
      function ftnk1(x,ndv,iof)
! ------------------------------
      IMPLICIT NONE
      INTEGER :: iof, ndv
      REAL(DP) :: ftnk1 
      REAL(DP) :: x(ndv)
      ftnk1 = x(iof)

      return
      end function ftnk1

! ------------------------------
      function fosy1(x,ndv,iof)
! ------------------------------
      IMPLICIT NONE
      INTEGER :: iof, ndv,idv
      REAL(DP) :: fosy1 
      REAL(DP) :: x(ndv)

      if(iof.eq.1)then
        fosy1 = -1.0_dp * ( 25.0_dp*(x(1)-2.0_dp)**2 + (x(2)-2.0_dp)**2 + &
           (x(3)-1.0_dp)**2 + (x(4)-4.0_dp)**2 + (x(5)-1.0_dp)**2 )   !20140619-3
      else
        fosy1 = 0.0_dp   !20140619-3
        do idv=1,ndv
          fosy1 = fosy1 + x(idv)**2
        enddo
      endif

      return
      end function fosy1

! ------------------------------
      function fsrn1(x,ndv,iof)
! ------------------------------
      IMPLICIT NONE
      INTEGER  :: iof, ndv
      REAL(DP) :: fsrn1 
      REAL(DP) :: x(ndv)
      REAL(DP) :: x1,x2

      x1 = x(1)
      x2 = x(2)

      if(iof.eq.1)then
        fsrn1 = 2.0_dp + (x1-2.0_dp)**2 + (x2-1.0_dp)**2   !20140620-1
      else
        fsrn1 = 9.0_dp*x1 - (x2-1.0_dp)**2   !20140620-1
      endif

      return
      end function fsrn1

! -------------------------------
      function ccex1(x,ndv,icn)
! -------------------------------
      IMPLICIT NONE
      INTEGER  :: icn, ndv
      REAL(DP) :: ccex1 
      REAL(DP) :: x(ndv)
      REAL(DP) :: x1,x2
      x1 = x(1)
      x2 = x(2)

      if(icn.eq.1)then
        ccex1 =  x2 + 9.0_dp*x1 - 6.0_dp   !20140620-1
      else
        ccex1 = -x2 + 9.0_dp*x1 - 1.0_dp   !20140620-1 
      endif
      ccex1 = (-1.0_dp) * ccex1   !20140620-1
      
      return
      end function ccex1

! -------------------------------
      function ctnk1(x,ndv,icn)
! -------------------------------
      IMPLICIT NONE
      INTEGER  :: icn, ndv
      REAL(DP) :: ctnk1 
      REAL(DP) :: x(ndv)
      REAL(DP) :: x1,x2
      x1 = x(1)
      x2 = x(2)

      if(icn.eq.1)then
        ctnk1 = x1**2 + x2**2 -1.0_dp -0.1_dp*cos(16.0_dp*atan(x1/x2))   !20140620-1
        ctnk1 = (-1.0_dp) * ctnk1   !20140620-1
      else
        ctnk1 = (x1-0.5_dp)**2 + (x2-0.5_dp)**2 - 0.5_dp   !20140620-1
      endif
      
      return
      end function ctnk1

! -------------------------------
      function cosy1(x,ndv,icn)
! -------------------------------
      IMPLICIT NONE
      INTEGER :: icn, ndv
      REAL(DP) :: cosy1 
      REAL(DP) :: x(ndv)

      if(icn.eq.1)then
        cosy1 = x(1) + x(2) - 2.0_dp   !20140620-1 
      elseif(icn.eq.2)then
        cosy1 = 6.0_dp - x(1) - x(2)   !20140620-1
      elseif(icn.eq.3)then
        cosy1 = 2.0_dp - x(2) + x(1)   !20140620-1
      elseif(icn.eq.4)then
        cosy1 = 2.0_dp - x(1) + 3.0_dp*x(2)   !20140620-1
      elseif(icn.eq.5)then
        cosy1 = 4.0_dp - (x(3)-3.0_dp)**2 - x(4)   !20140620-1
      else
        cosy1 = (x(5)-3.0_dp)**2 + x(6) -4.0_dp   !20140620-1
      endif
      cosy1 = (-1.0_dp) * cosy1   !20140620-1
      
      return
      end function cosy1

! -------------------------------
      function csrn1(x,ndv,icn)
! -------------------------------
      IMPLICIT NONE
      INTEGER  :: icn, ndv
      REAL(DP) :: csrn1 
      REAL(DP) :: x(ndv)
      REAL(DP) :: x1,x2
      x1 = x(1)
      x2 = x(2)

      if(icn.eq.1)then
        csrn1 = x1 - 3.0_dp*x2 + 10.0_dp   !20140620-1
      else
        csrn1 = x1**2 + x2**2 - 225.0_dp   !20140620-1
      endif
      
      return
      end function csrn1

!------------------------------------------------------------------
!
!
! ORIGINAL: messy_airtraf_tools_ga_testdirect.f90
!
!
! *******************************************************
!     Subroutines for connecting with C program 
! *******************************************************
!       YOU MAY CHANGE NAME OF SUBROUINTINE  
!         ACCORDING TO SYSTEM AND COMPILER
! -------------------------------------------------------
! --------------------------------------------------
!HY20140327-3      subroutine DirectCF(fmobj,nobj,x,ndv,testp)
! --------------------------------------------------
!HY      dimension fmobj(nobj), x(ndv)
!HY      character*4 testp
      
!HY      if(testp.eq.'RRT0')call testp_rrmot0(x,fmobj)
!HY      if(testp.eq.'RRT1')call testp_rrmot1(x,fmobj)
!HY      if(testp.eq.'RRT2')call testp_rrmot2(x,fmobj)
!HY      if(testp.eq.'RRT3')call testp_rrmot3(x,fmobj)
!HY      if(testp.eq.'RRT4')call testp_rrmot1(x,fmobj)
!HY      if(testp.eq.'RRT5')call testp_rrmot5(x,fmobj)
!HY      if(testp.eq.'RRT6')call testp_rrmot6(x,fmobj)
!HY      if(testp.eq.'RRT7')call testp_rrmot7(x,fmobj)
!HY      if(testp.eq.'RRT8')call testp_rrmot8(x,fmobj)
!HY      if(testp.eq.'RRT9')call testp_rrmot9(x,fmobj)
!HY      if(testp.eq.'RRC1')call testcf_rrmot1(x,fmobj)
!HY      if(testp.eq.'RRC2')call testcf_rrmot2(x,fmobj)
!HY      if(testp.eq.'RRC3')call testcf_rrmot3(x,fmobj)

!HY      if(testp.eq.'MT21')call testp_mot21(x,fmobj)
!HY      if(testp.eq.'MT22')call testp_mot22(x,fmobj)
!HY      if(testp.eq.'MT23')call testp_mot23(x,fmobj)
!HY      if(testp.eq.'MT24')call testp_mot24(x,fmobj)
!HY      if(testp.eq.'MT25')call testp_mot25(ndv,x,fmobj)

!C2      if(testp.eq.'RRT0')call testp_rrmot0_(x,fmobj)
!C2      if(testp.eq.'RRT1')call testp_rrmot1_(x,fmobj)
!C2      if(testp.eq.'RRT2')call testp_rrmot2_(x,fmobj)
!C2      if(testp.eq.'RRT3')call testp_rrmot3_(x,fmobj)
!C2      if(testp.eq.'RRT4')call testp_rrmot4_(ndv,x,fmobj)
!C2      if(testp.eq.'RRT5')call testp_rrmot5_(x,fmobj)
!C2      if(testp.eq.'RRT6')call testp_rrmot6_(x,fmobj)
!C2      if(testp.eq.'RRT7')call testp_rrmot7_(x,fmobj)
!C2      if(testp.eq.'RRT8')call testp_rrmot8_(x,fmobj)
!C2      if(testp.eq.'RRT9')call testp_rrmot9_(x,fmobj)
!C2      if(testp.eq.'RRC1')call testcf_rrmot1_(x,fmobj)
!C2      if(testp.eq.'RRC2')call testcf_rrmot2_(x,fmobj)
!C2      if(testp.eq.'RRC3')call testcf_rrmot3_(x,fmobj)
   
!HY      return
!HY      end

! -------------------------------------------------------
!HY      subroutine DirectCG(gv,ncon,fmobj,nobj,x,ndv,testp)
! -------------------------------------------------------
!HY      dimension gv(ncon), fmobj(nobj), x(ndv)
!HY      character*4 testp
      
!HY      if(testp.eq.'RRC1')call testcg_rrmot1(x,fmobj,gv)
!HY      if(testp.eq.'RRC2')call testcg_rrmot2(x,fmobj,gv)
!HY      if(testp.eq.'RRC3')call testcg_rrmot3(x,fmobj,gv)

!CLINUX      if(testp.eq.'RRC1')call testcg_rrmot1_(x,fmobj,gv)
!CLINUX      if(testp.eq.'RRC2')call testcg_rrmot2_(x,fmobj,gv)
!CLINUX      if(testp.eq.'RRC3')call testcg_rrmot3_(x,fmobj,gv)
   
!HY      return
!HY20140327-3      end

! --------------------------------------------
      subroutine SetDesign(set_p2_ac_routes_desc)  !HY20140418-2
! --------------------------------------------
      USE messy_airtraf_tools_ga_input,      ONLY:nobj_fort_42,ncon_fort_42,ndv_fort_42,ndetail_armogaset !HY20140402-2
      USE messy_airtraf_tools_ga_armogaset,  ONLY:calc_boundary
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
!HY20140401-7      include 'ga.par'
!      dimension nobjst(misland),nobjed(misland)
!      dimension nsubst(misland),nsubed(misland)
!      dimension nconst(misland),nconed(misland)
!      dimension ngcnst(misland),ngcned(misland)
      REAL(DP), DIMENSION(:),   INTENT(INOUT)  :: set_p2_ac_routes_desc   !HY20140418-2
      INTEGER :: i,ilnum,num,numcon,numobj
      nlinear=0
      ngeocon=0
      do ilnum=1,nisland
        do i=1,8
          imoga(ilnum,i)=0
        enddo
      enddo

!HY20140401-2      open(10,file='armoga.set')
!HY      read(10,*)numobj
!HY      read(10,*)numcon
!HY      read(10,*)ndv
      numobj = nobj_fort_42  !HY20140401-2
      numcon = ncon_fort_42
      ndv    = ndv_fort_42
!HY20140401-2      do idv=1,ndv
!HY        read(10,*)phlb(idv),phub(idv)
!HY      enddo
      call calc_boundary(ndv, phlb, phub, set_p2_ac_routes_desc)   !HY20140401-3,20140418-2   
!HY20140401-3      read(10,*)ndetail
      ndetail = ndetail_armogaset   !HY20140401-3
      if(ndetail.eq.1)then
!HY20140401-3        do i=1,4
!HY          read(10,*)(imoga(ilnum,i),ilnum=1,nisland)
!HY        enddo
!HY        if(numcon.ne.0)then
!HY          do i=5,6
!HY            read(10,*)(imoga(ilnum,i),ilnum=1,nisland)
!HY          enddo
!HY          if(icnrgeo.eq.1)then
!HY            do i=7,8
!HY              read(10,*)(imoga(ilnum,i),ilnum=1,nisland)
!HY            enddo
!HY          endif
!HY        endif
      else
        do ilnum=1,nisland
          imoga(ilnum,1)=1
          imoga(ilnum,2)=abs(numobj)
          imoga(ilnum,3)=1
          imoga(ilnum,4)=abs(numobj)
          if(numcon.ne.0)imoga(ilnum,5)=1
          imoga(ilnum,6)=abs(numcon)
        enddo
      endif
!HY20140401-3      close(10)
  
      if(numobj.lt.0)then
        imaxmin = -1
      else
        imaxmin =  1
      endif
      nobj    = abs(numobj)

      if(numcon.lt.0)then
        igwtobj = 1
      else
        igwtobj = 0
      endif
      ncon    = abs(numcon)

      do ilnum=1,nisland
        if(imoga(ilnum,7).gt.0)then
          num=imoga(ilnum,8)-imoga(ilnum,7)+1
          if(num.gt.ngeocon)ngeocon=num
        endif
      enddo

!cPREV      open(10,file='armoga.set')
!cPREV      read(10,*)numobj
!cPREV      read(10,*)numcon
!cPREV      read(10,*)ndv
!cPREV      do idv=1,ndv
!cPREV        read(10,*)phlb(idv),phub(idv)
!cPREV      enddo
!cPREV      read(10,*)ndetail
!cPREV      if(ndetail.eq.1)then
!cPREV        read(10,*)nlinear
!cPREV        read(10,*)ngeocon
!cPREV      else
!cPREV        nlinear = 0
!cPREV        ngeocon = 0
!cPREV      endif 
      
      return
      end subroutine SetDesign

! --------------------------------------------
      subroutine SystemExe(fmobj,nobj,x,ndv)
! --------------------------------------------
      IMPLICIT NONE
      INTEGER :: nobj,ndv,idv,iof 
      REAL(DP) :: fmobj(nobj), x(ndv)

      open(11,file='armoga.dv')
      do idv=1,ndv
        write(11,*)x(idv)
      enddo
      close(11)

!HY20140328-5      call cs
!CLINUX      call cs_
   
      open(12,file='armoga.obj')
      do iof=1,nobj
        read(12,*)fmobj(iof)
      enddo
      close(12)

      return 
      end subroutine SystemExe

! --------------------------------------------
      subroutine SystemCstr(gv,ncon)
! --------------------------------------------
      IMPLICIT NONE
      INTEGER :: ncon,icn
      REAL(DP) :: gv(ncon)
!HY20140328-5      call cs_cstr
!CLINUX      Call cs_cstr_
   
!HY20140328-5      open(13,file='armoga.cstr')
      do icn=1,ncon
        read(13,*)gv(icn)
      enddo
      close(13)

      return
      end subroutine SystemCstr

! --------------------------------------------
      subroutine SystemCstrAll(gm,ncon,npop)
! --------------------------------------------
      IMPLICIT NONE
      INTEGER :: ncon,npop,i
      REAL(DP) :: gm(npop*ncon)
!HY20140328-5      call cs_cstr
!CLINUX      call cs_cstr_
   
!HY20140328-5      open(13,file='armoga.cstr')
      do i=1,npop*ncon
        read(13,*)gm(i)
      enddo
      close(13)

      return 
      end subroutine SystemCstrAll

! ----------------------------------------------------------------
      subroutine SystemExec(x,ndv,objf,nobj,gv,ncon,ig,npop,itype, &
                            sys_nwaypoints, sys_p2_ac_routes_desc, sys_philon, sys_philat, sys_zgl_geopot_3d, &
                            sys_uwind_g, sys_vwind_g, sys_v_z_g, sys_t_scb_g,sys_cpc_g, & !20140420-1,20140522-3,20160302-2
                            sys_rho_air_dry_3d_g, sys_press_3d_g,sys_ATR20O3_3d_g,&
                            sys_ATR20CH4_3d_g,sys_ATR20H2O_3d_g,sys_ATR20CPC_3d_g,sys_ATR20CO2_3d_g) 
! ----------------------------------------------------------------
      USE messy_airtraf_wind,      ONLY:calcobj   !20140417-6
      !Yin_20170423      
!     USE messy_airtraf_MO,        ONLY:calcobjMO   !20140417-6
      !Yin_20170423      
      IMPLICIT NONE
      INTEGER :: ig,itype,npop,ndv,nobj,ncon,i,ipop,nsum1,nsum2,nsum3
      REAL(DP)               :: x(npop*ndv)
      REAL(DP)               :: objf(npop*nobj)
      REAL(DP)               :: gv(npop*ncon+1)
  INTEGER, INTENT(INOUT),OPTIONAL                           :: sys_nwaypoints  !20140721-7
  REAL(DP), DIMENSION(:),   INTENT(INOUT),OPTIONAL          :: sys_p2_ac_routes_desc  !HY20140420-1,20140422-4
  REAL(DP), DIMENSION(:),   INTENT(INOUT),OPTIONAL          :: sys_philon
  REAL(DP), DIMENSION(:),   INTENT(INOUT),OPTIONAL          :: sys_philat
  REAL(DP), DIMENSION(:,:,:), INTENT(INOUT),OPTIONAL        :: sys_zgl_geopot_3d   !global field 
  REAL(DP), DIMENSION(:,:,:), INTENT(INOUT),OPTIONAL        :: sys_uwind_g         !global field 
  REAL(DP), DIMENSION(:,:,:), INTENT(INOUT),OPTIONAL        :: sys_vwind_g         !global field
  REAL(DP), DIMENSION(:,:,:), INTENT(INOUT),OPTIONAL        :: sys_v_z_g           !global field 
  REAL(DP), DIMENSION(:,:,:), INTENT(INOUT),OPTIONAL        :: sys_t_scb_g         !global field  20140522-3 
  !Yin_20170423
  REAL(DP), DIMENSION(:,:,:), INTENT(INOUT),OPTIONAL        :: sys_cpc_g         !global field  20140522-3 
  !Yin_20170423
  REAL(DP), DIMENSION(:,:,:), INTENT(IN),OPTIONAL           :: sys_rho_air_dry_3d_g  !global field  20160302-2,20170206-1 
  REAL(DP), DIMENSION(:,:,:), INTENT(IN),OPTIONAL           :: sys_press_3d_g      !global field  20170215-4 
  !Yin_20170801+
  REAL(DP), DIMENSION(:,:,:), INTENT(IN),OPTIONAL           :: sys_ATR20O3_3d_g      !global field  20170801 
  REAL(DP), DIMENSION(:,:,:), INTENT(IN),OPTIONAL           :: sys_ATR20CH4_3d_g      !global field  20170801 
  REAL(DP), DIMENSION(:,:,:), INTENT(IN),OPTIONAL           :: sys_ATR20H2O_3d_g      !global field  20170801 
  REAL(DP), DIMENSION(:,:,:), INTENT(IN),OPTIONAL           :: sys_ATR20CPC_3d_g      !global field  20170801 
  REAL(DP), DIMENSION(:,:,:), INTENT(IN),OPTIONAL           :: sys_ATR20CO2_3d_g      !global field  20170801 
  !Yin_20170801-
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                   :: temp_sys_rho_air_dry_3d_g  !global field  20170206-1 
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                   :: temp_sys_press_3d_g  !global field  20170215-4 
  !Yin_20170801+
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                   :: temp_sys_ATR20O3_3d_g  !global field  20170801 
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                   :: temp_sys_ATR20CH4_3d_g  !global field  20170801 
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                   :: temp_sys_ATR20H2O_3d_g  !global field  20170801 
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                   :: temp_sys_ATR20CPC_3d_g  !global field  20170801
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                   :: temp_sys_ATR20CO2_3d_g  !global field  20170801
  !Yin_20170801- 
      nsum1=1
      nsum2=1
      nsum3=1
 
      IF(PRESENT(sys_rho_air_dry_3d_g))THEN   !20170206-1 for Fuel_opt, NOx_opt, H2O_opt, Cost_opt, COC_opt
         ALLOCATE(temp_sys_rho_air_dry_3d_g(nlon_ga, nlev_ga, ngl_ga))
         temp_sys_rho_air_dry_3d_g=sys_rho_air_dry_3d_g
!         write(*,*)"allocation4:nlon_ga,nlev_ga,ngl_ga",nlon_ga,nlev_ga,ngl_ga
      ENDIF
      IF(PRESENT(sys_press_3d_g))THEN         !20170215-4 for NOx_opt
         ALLOCATE(temp_sys_press_3d_g(nlon_ga, nlev_ga, ngl_ga))
         temp_sys_press_3d_g=sys_press_3d_g
      ENDIF
      !Yin_20170801+
      IF(PRESENT(sys_ATR20O3_3d_g))THEN         !20170801 for ATR20_opt
         ALLOCATE(temp_sys_ATR20O3_3d_g(nlon_ga, nlev_ga, ngl_ga))
         temp_sys_ATR20O3_3d_g=sys_ATR20O3_3d_g
      ENDIF
      IF(PRESENT(sys_ATR20CH4_3d_g))THEN         !20170801 for ATR20_opt
         ALLOCATE(temp_sys_ATR20CH4_3d_g(nlon_ga, nlev_ga, ngl_ga))
         temp_sys_ATR20CH4_3d_g=sys_ATR20CH4_3d_g
      ENDIF
      IF(PRESENT(sys_ATR20H2O_3d_g))THEN         !20170801 for ATR20_opt
         ALLOCATE(temp_sys_ATR20H2O_3d_g(nlon_ga, nlev_ga, ngl_ga))
         temp_sys_ATR20H2O_3d_g=sys_ATR20H2O_3d_g
      ENDIF
      IF(PRESENT(sys_ATR20CPC_3d_g))THEN         !20170801 for ATR20_opt
         ALLOCATE(temp_sys_ATR20CPC_3d_g(nlon_ga, nlev_ga, ngl_ga))
         temp_sys_ATR20CPC_3d_g=sys_ATR20CPC_3d_g
      ENDIF
      IF(PRESENT(sys_ATR20CO2_3d_g))THEN         !20170801 for ATR20_opt
         ALLOCATE(temp_sys_ATR20CO2_3d_g(nlon_ga, nlev_ga, ngl_ga))
         temp_sys_ATR20CO2_3d_g=sys_ATR20CO2_3d_g
      ENDIF

      !Yin_20170801-
      if(itype.gt.0)then
!HY20140327-4,20140328-3
!HY     open(11,file='armoga.db')
!HY     open(11,file='mot.input')
!HY     do ipop=1,npop
!HY       write(11,*)ig,ipop,1
!HY       write(11,*)nobj,ncon,ndv
!HY          do i=1,nobj
!HY            write(11,*)objf(nsum2)
!HY            nsum2 = nsum2 + 1
!HY          enddo
!HY          if(ncon.gt.0)then
!HY            do i=1,ncon
!HY              write(11,*)gv(nsum3)
!HY              nsum3 = nsum3 + 1
!HY            enddo
!HY          endif
!HY          do i=1,ndv
!HY            write(11,*)x(nsum1)
!HY            nsum1 = nsum1 + 1
!HY          enddo
!HY     enddo
!HY     close(11)

        if(itype.eq.2)then
!HY20140328-3       call cs
!HY20140420-1       call calcobj(x,ndv,objf,nobj,gv,ncon,ig,npop)
          SELECT CASE(ga_option_traj_calc_g)   !20160302-3
          CASE(1,5)   !Wind_opt, CPC_opt Yin_20170423, HY 20170609-1
          call calcobj(x,ndv,objf,nobj,gv,ncon,ig,npop,  &
                       sys_nwaypoints, sys_p2_ac_routes_desc, sys_philon, sys_philat, sys_zgl_geopot_3d, &
                       sys_uwind_g, sys_vwind_g, sys_v_z_g, sys_t_scb_g, sys_cpc_g, &  !HY20140420-1, 20140522-3,20140721-7 
                       ga_option_traj_calc_g)                               !20160525-1 
          CASE(2,4,6,7,10)   !Fuel_opt, H2O_opt, Cost_opt or COC_opt   !20170222-2,20170224-2,20170316-2  
          !IF(PRESENT(sys_rho_air_dry_3d_g))temp_sys_rho_air_dry_3d_g=sys_rho_air_dry_3d_g   !20170206-1
          call calcobj(x,ndv,objf,nobj,gv,ncon,ig,npop,  &
                       sys_nwaypoints, sys_p2_ac_routes_desc, sys_philon, sys_philat, sys_zgl_geopot_3d, &
                       sys_uwind_g, sys_vwind_g, sys_v_z_g, sys_t_scb_g, sys_cpc_g, &  !HY20140420-1, 20140522-3,20140721-7 
                       ga_option_traj_calc_g, temp_sys_rho_air_dry_3d_g)         !20160302-3,20160525-1,20170206-1
!          write(*,*)'check_ga_option_traj_calc_5', ga_option_traj_calc_g    !20160304-1
          CASE(3)   !NOx_opt  
          !IF(PRESENT(sys_rho_air_dry_3d_g))temp_sys_rho_air_dry_3d_g=sys_rho_air_dry_3d_g   !20170206-1
          call calcobj(x,ndv,objf,nobj,gv,ncon,ig,npop,  &
                       sys_nwaypoints, sys_p2_ac_routes_desc, sys_philon, sys_philat, sys_zgl_geopot_3d, &
                       sys_uwind_g, sys_vwind_g, sys_v_z_g, sys_t_scb_g, sys_cpc_g, &  !HY20140420-1, 20140522-3,20140721-7 
                       ga_option_traj_calc_g, temp_sys_rho_air_dry_3d_g, &  !20160302-3,20160525-1,20170206-1
                       temp_sys_press_3d_g)                                 !20170215-4 
!          write(*,*)'check_ga_option_traj_calc_5', ga_option_traj_calc_g    !20160304-1
          !Yin_20170801+
          CASE(8,9)   !ATR20_opt,COSTCLIM_opt
          call calcobj(x,ndv,objf,nobj,gv,ncon,ig,npop,  &
                      sys_nwaypoints, sys_p2_ac_routes_desc, sys_philon, sys_philat, sys_zgl_geopot_3d, &
                      sys_uwind_g, sys_vwind_g, sys_v_z_g, sys_t_scb_g,sys_cpc_g,&
                      ga_option_traj_calc_g,temp_sys_rho_air_dry_3d_g, &
                      temp_sys_press_3d_g,temp_sys_ATR20O3_3d_g,temp_sys_ATR20CH4_3d_g,&
                      temp_sys_ATR20H2O_3d_g,temp_sys_ATR20CPC_3d_g,temp_sys_ATR20CO2_3d_g)
          !Yin_20170801-
          CASE DEFAULT
             write(*,*)'ERROR:ga_option_traj_calc on SystemExec',ga_option_traj_calc_g
          ENDSELECT
!CLINUX          call cs_
        else
!HY20140328-5          call cs_geo
!CLINUX          call cs_geo_
        endif

      endif

      nsum1=1
      nsum2=1
      nsum3=1

!HY20140328-3   open(12,file='armoga.db')
!HY   open(12,file='mot.output')
      do ipop=1,npop
!HY     read(12,*)id1,id2,id3
!HY     read(12,*)id1,id2,id3
        do i=1,nobj
!HY          read(12,*)objf(nsum2)
          !print *,'GA,F:',ipop,i,objf(nsum2)   !To check optimization process on display.
          nsum2 = nsum2 + 1
        enddo
!HY        if(ncon.gt.0)then
!HY          do i=1,ncon
!HY            read(12,*)gv(nsum3)
!HY            print *,'GA,G:',ipop,i,gv(nsum3)
!HY            nsum3 = nsum3 + 1
!HY          enddo
!HY        endif
!HY     do i=1,ndv
!HY       read(12,*)xd
!HY     enddo
      enddo
!HY20140327-4,20140328-3   close(12)
      !for fuel_opt, NOx_opt, H2O_opt, Cost_opt, COC_opt     
      IF(ALLOCATED(temp_sys_rho_air_dry_3d_g))DEALLOCATE(temp_sys_rho_air_dry_3d_g)   !20170206-1,20170212-3
      !for NOx_opt      
      IF(ALLOCATED(temp_sys_press_3d_g))DEALLOCATE(temp_sys_press_3d_g)               !20170215-4
      !Yin_20170801+
      IF(ALLOCATED(temp_sys_ATR20O3_3d_g))DEALLOCATE(temp_sys_ATR20O3_3d_g)               !20170801
      IF(ALLOCATED(temp_sys_ATR20CH4_3d_g))DEALLOCATE(temp_sys_ATR20CH4_3d_g)               !20170801
      IF(ALLOCATED(temp_sys_ATR20H2O_3d_g))DEALLOCATE(temp_sys_ATR20H2O_3d_g)               !20170801
      IF(ALLOCATED(temp_sys_ATR20CPC_3d_g))DEALLOCATE(temp_sys_ATR20CPC_3d_g)               !20170801
      IF(ALLOCATED(temp_sys_ATR20CO2_3d_g))DEALLOCATE(temp_sys_ATR20CO2_3d_g)               !20170801
      !Yin_20170801-
      return 
      end subroutine SystemExec

! --------------------------------------------
      subroutine SystemGeo(gv,ncon,x,ndv)
! --------------------------------------------
      IMPLICIT NONE
      INTEGER :: ncon, ndv,icn,idv
      REAL(DP) :: gv(ncon), x(ndv)
      open(11,file='armoga.dv')
      do idv=1,ndv
        write(11,*)x(idv)
      enddo
      close(11)

!HY20140328-5      call cs_geo
!CLINUX      call cs_geo_
   
!HY20140328-5      open(13,file='armoga.cstr')
      do icn=1,ncon
        read(13,*)gv(icn)
      enddo
      close(13)

      return 
      end subroutine SystemGeo

!-----------------------------------------------------------
!
!
! ORIGINAL: messy_airtraf_tools_ga_update.f90 
!
!
! *******************************
      subroutine NextGen(ilnum)
! *******************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      INTEGER :: ilnum,ist,ied,ip,ipop,new,newpop
      newpop=istat(ilnum,11)
      do ip=1,newpop
        new=istat(ilnum,13)-istat(ilnum,11)+ip
        ipop=idarc(ilnum,new)
        ifitpop(ilnum,ip)=ipop
!C1.2.0        ipinfo(ipop,1)=igen+1
        ipinfo(ipop,1)=iter
        ipinfo(ipop,2)=istat(ilnum,10)+ip-1
        ipinfo(ipop,3)=ilnum
        ipinfo(ipop,4)=0
        ifrank(ilnum,ip)=0
        ffit(ilnum,ip)=0.0_dp   !20140620-1
      enddo  

      ist=istat(ilnum,13)-istat(ilnum,11)+1
      ied=istat(ilnum,13)
      call InitialObj(ist,ied)
      
      istat(ilnum,10)=istat(ilnum,10)+npop

      return
      end subroutine NextGen

! *******************************
      subroutine NextFPop(ilnum)
! *******************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      INTEGER :: ilnum,inum,new,newpop,nmax
      newpop=istat(ilnum,11)
      nmax = (ibestn+1)*newpop
      inum = istat(ilnum,12)
      do while (inum.gt.0)
        new = inum + newpop
        if(new.le.nmax)then
          ifitpop(ilnum,new)=ifitpop(ilnum,inum)
          if(new.gt.istat(ilnum,12))istat(ilnum,12)=new
        endif
        inum=inum-1
      enddo
!CHECK
      istat(ilnum,2)=ibestn   

      return
      end subroutine NextFPop

! ************************************************
      subroutine StatSetting(ilnum)
! ************************************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      INTEGER :: istart, inter,ilnum,i1,i2,i3,igs
      i1=istat(ilnum,11)
      i2=istat(ilnum,12)
      i3=istat(ilnum,13)
 
! ------  ISTAT(1)  <ADAPTIVE RANGE>  ------
      istat(ilnum,1)=0
      istart = iast(ilnum)
!C1.0.1      igs=igen
      igs=igen+1
!C1.1.4      if(igs.ge.istart)then
!C1.2.0      if(igen.ge.istart)then
      if(iter.ge.istart)then
        inter = iaint(ilnum)
!C1.1.4        istat(ilnum,1)=mod(igs-istart,inter)+1
!C1.2.0        istat(ilnum,1)=mod(igen-istart,inter)+1
        istat(ilnum,1)=mod(iter-istart,inter)+1
      endif
 
! ------  ISTAT(2)  <FITNESS NUMBER (BEST N)>  ------
      istat(ilnum,2)=0
      if((i2-i1).gt.0)istat(ilnum,2)=ibestn
 
! ------  ISTAT(3)  <ARCHIVING>  ------
      istat(ilnum,3)=0
!cMAY      if(i3-i1.gt.0)istat(ilnum,3)=iarchiv
      if(i3-i1.gt.0)istat(ilnum,3)=ibnarc

! ------  ISTAT(4&5)  <ACTIVE OBJECTIVES>  ------
      istat(ilnum,4)=imoga(ilnum,1)
      istat(ilnum,5)=imoga(ilnum,2)
      if((isubar.eq.1).and.(istat(ilnum,1).eq.1))then
        istat(ilnum,4)=imoga(ilnum,3)
        istat(ilnum,5)=imoga(ilnum,4)
      endif

! ------  ISTAT(6&7)  <ACTIVE CONSTRAINTS>  ------
! ------         BEFORE RANKING             ------
      istat(ilnum,6)=imoga(ilnum,5)
      istat(ilnum,7)=imoga(ilnum,6)
 
! ------  ISTAT(9)  <RANKING METHOD FOR ARCHIVING>  ------
! ------     -1 ... Only Rank 1 is Marked (Efficient)
! ------      0 ... New Ranking
! ------      1 ... Ranking Based on previous Ranking
      istat(ilnum,9)=0
!CMAY
      if(istat(ilnum,2).eq.0)istat(ilnum,9)=0
      if(istat(ilnum,1).eq.2)istat(ilnum,9)=0
      if((isubar.eq.1).and.(istat(ilnum,1).eq.1))istat(ilnum,9)=1
      if(irank1.eq.1)istat(ilnum,9)=-1
 
      return
      end subroutine StatSetting   

! **********************************
      subroutine InitialObj(ist,ied)
! **********************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      INTEGER :: ist,ied,icn,iof,ipop
      do ipop=ist,ied
        do iof=1,nobj
          fvall(ipop,iof) = objinit
        enddo
        if(ncon.gt.0)then
          do icn=1,ncon
            gvall(ipop,icn) = gvalinit
          enddo
        endif
      enddo

      return
      end subroutine InitialObj

! ****************************************************
      subroutine MaxminObj(ilnum,numpop,fdec) 
! ****************************************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      REAL(DP) :: fdec(maxpop,mobj)
      REAL(DP) :: objsum,value 
      INTEGER  :: ilnum,numpop,istobj,iedobj,iof,ip,ipop
!      numpop = istat(ilnum,12)
 
      istobj = istat(ilnum,4)
      iedobj = istat(ilnum,5)
 
      ipop=ifitpop(ilnum,1)
      do iof=istobj,iedobj
        objmax(ilnum,iof)=fvall(ipop,iof)
        objmin(ilnum,iof)=fvall(ipop,iof)
      enddo
 
      do iof=istobj,iedobj
        objsum = 0.0_dp   !20140620-1
        do ip=1,numpop
          ipop=ifitpop(ilnum,ip)
          value=fvall(ipop,iof)
          fdec(ip,iof)=value
          if(value.gt.objmax(ilnum,iof))objmax(ilnum,iof)=value
          if(value.lt.objmin(ilnum,iof))objmin(ilnum,iof)=value
          objsum = objsum + value
        enddo
        objave(ilnum,iof) = objsum / DBLE(numpop)     !HY20140414-5
        if(objmax(ilnum,iof).gt.objallmax(ilnum,iof)) &
           objallmax(ilnum,iof) = objmax(ilnum,iof)
        if(objmin(ilnum,iof).lt.objallmin(ilnum,iof)) &
           objallmin(ilnum,iof) = objmin(ilnum,iof)
      enddo
   
!      if(infopt.eq.1)then
!        do iof=istobj,iedobj
!          write(mtout,*)'OBJ:', iof, &
!             objmin(ilnum,iof),objmax(ilnum,iof)
!        enddo
!      endif
 
      return
      end subroutine MaxminObj  

! ****************************************************
      subroutine MaxminGval(ilnum,numpop,gdec)
! ****************************************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      REAL(DP) :: gdec(maxpop,mcon)
      REAL(DP) :: gvalsum,value 
      INTEGER  :: ilnum,numpop,icn,istcon,iedcon,ip,ipop
!      numpop = istat(ilnum,12)
      istcon = istat(ilnum,6)
      iedcon = istat(ilnum,7)
 
      ipop=ifitpop(ilnum,1)
      do icn=istcon,iedcon
        gvalmax(ilnum,icn)=gvall(ipop,icn)
        gvalmin(ilnum,icn)=gvall(ipop,icn)
      enddo
 
      do icn=istcon,iedcon
        gvalsum = 0.0_dp   !20140620-1
        do ip=1,numpop
          ipop = ifitpop(ilnum,ip)
          value = gvall(ipop,icn)
          gdec(ip,icn) = value
          if(value.gt.gvalmax(ilnum,icn)) &
             gvalmax(ilnum,icn)=value
          if(value.lt.gvalmin(ilnum,icn)) &
             gvalmin(ilnum,icn)=value
          gvalsum = gvalsum + value
          if(gvalmax(ilnum,icn).gt.gvalallmax(ilnum,icn)) &
             gvalallmax(ilnum,icn) = gvalmax(ilnum,icn)
          if(gvalmin(ilnum,icn).lt.gvalallmin(ilnum,icn)) &
             gvalallmin(ilnum,icn) = gvalmin(ilnum,icn)
        enddo
        gvalave(ilnum,icn)= gvalsum / DBLE(numpop)  !HY20140414-5   
      enddo
 
      if(infopt.eq.1)then
        do icn=istcon,iedcon
         write(mtout,*)'CON:',icn, &
          gvalmin(ilnum,icn), gvalmax(ilnum,icn)
        enddo
      endif
 
      return
      end subroutine MaxminGval    

!------------------------------------------------------------------
!
!
! ORIGINAL: messy_airtraf_tools_ga_archive.f90
!
!
! **************************************************************
      subroutine ConstrRankAdd(ilnum)
! **************************************************************
!  ..... calc {n*(n-1)}/2  (?)  ...
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      REAL(DP) :: gn(mcon),gf(mcon)
      REAL(DP) :: fj,gj,bd
      INTEGER  :: jump,ncpop,npfin,nn,ifin,itotal,nf,ilnum,icn,istcon,iedcon,istobj,iedobj
      INTEGER  :: inew,iof,ip,numobj
      jump=0
      bd = 1.0_dp   !20140620-1
      if(imaxmin.lt.0)bd = -1.0_dp   !20140620-1

      istobj=istat(ilnum,4)
      iedobj=istat(ilnum,5)
      istcon=istat(ilnum,6)
      iedcon=istat(ilnum,7)
      numobj=iedobj-istobj+1


      if(istat(ilnum,9).eq.1)then
        ncpop = istat(ilnum,11)
        npfin = istat(ilnum,13)-istat(ilnum,11)
      else
        ncpop = istat(ilnum,13)
        npfin = 0
        istat(ilnum,14)=0
      endif
      if(istat(ilnum,9).eq.-1)jump=1


      do ip=1,ncpop
        inew=npfin+ip
        irank(ilnum,inew)=1
        if(ncon.gt.0)then
          nn=0
          do icn=istcon,iedcon
            if(gvall(idarc(ilnum,inew),icn).gt.0)nn=nn+1
          enddo
          if(nn.eq.0)istat(ilnum,14)=istat(ilnum,14)+1
        endif
      enddo

      do ip=1,ncpop
        inew = npfin + ip
        ifin = 1
        itotal = npfin + ip - 1 

        do while(ifin.le.itotal)
          nn = 0
          nf = 0
          if(istcon.gt.0)then
            do icn=istcon,iedcon
              gn(icn) = max(0.0_dp,gvall(idarc(ilnum,inew),icn))    !20140620-1
              gf(icn) = max(0.0_dp,gvall(idarc(ilnum,ifin),icn))    !20140620-1
              if(gn(icn).gt.0.0_dp)nn = nn + 1   !20140620-1
              if(gf(icn).gt.0.0_dp)nf = nf + 1   !20140620-1
            enddo
          endif
          if((nn.eq.0).and.(nf.eq.0))then
            do iof=istobj,iedobj
              fj = bd * & 
            (fvall(idarc(ilnum,inew),iof)-fvall(idarc(ilnum,ifin),iof))
              if(fj.le.0.0_dp)nn = nn + 1   !20140620-1
              if(fj.ge.0.0_dp)nf = nf + 1   !20140620-1
            enddo
            if((nn.eq.numobj).and.(nf.eq.numobj))goto 10
            if(nn.eq.numobj)then
              irank(ilnum,inew)=irank(ilnum,inew)+1
              if(jump.eq.1)goto 20
            endif
            if(nf.eq.numobj)irank(ilnum,ifin)=irank(ilnum,ifin)+1
          else
!CDS--ST
            if(nn.eq.nf)then
              do icn=istcon,iedcon
                gj = gn(icn) - gf(icn)
                if(gj.le.0.0_dp)nf = nf + 1   !20140620-2
                if(gj.ge.0.0_dp)nn = nn + 1   !20140620-2
              enddo
            endif
!CDS--ED
            if(nn.gt.nf)then
              irank(ilnum,inew)=irank(ilnum,inew)+1
              if(jump.eq.1)goto 20
            endif
            if(nn.lt.nf)irank(ilnum,ifin)=irank(ilnum,ifin)+1
          endif
  10      ifin = ifin + 1
        enddo
  20    continue
      enddo

      return
      end subroutine ConstrRankAdd

! ********************************************************
      subroutine DrawArchive(ilnum,inst,ined,iarcst,iarced,narcin)
! ********************************************************
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6       USE messy_airtraf_tools_ga_adapt   !,        ONLY: Encode
!HY20140404-6       USE messy_airtraf_tools_ga_function,     ONLY: ran1
      IMPLICIT NONE
      INTEGER :: id(mtpop),nsum(mtpop)
      REAL(DP) :: gtype,ptype
      INTEGER  :: insum,nummax,iarank,iarcpop,narcin,iarcst,iarced,inst,ined,ilnum,iarc,ip
      INTEGER  :: idv,ipop,nselect
      insum = ined - inst + 1
      nummax = 1
      iarank = 1
      nselect = 0
      call SelectArchive(ilnum,iarank,iarank,iarced,id,nselect)

      do ip=1,nselect
        nsum(ip)=0
      enddo

      narcin=nselect
      if(ismpop.ne.1)then
        do iarc=1,nselect
          iarcpop=idarc(ilnum,id(iarc))
          do ip=1,inst
            ipop=ifitpop(ilnum,ip)
            if(iarcpop.eq.ipop)then
              narcin=narcin-1
              nsum(iarcpop)=nsum(iarcpop)+1
              goto 11
            endif
          enddo
 11       continue
        enddo
      endif
     
      if(narcin.lt.insum)ined=ined-narcin+1
      if(narcin.le.0)return

      ip=inst
      do while(ip.le.ined)
        iarc = int(ran1(idum)*DBLE(nselect))+1   !HY20140414-5
        if(iarc.gt.nselect)iarc=nselect
        if(nsum(iarc).lt.nummax)then
          ipop=idarc(ilnum,id(iarc))
          ifitpop(ilnum,ip)=ipop
          if(ipinfo(ipop,4).eq.2)then
            do idv=1,ndv
              ptype=ptall(ipop,idv)
              gtype=Encode( ptype, pout(ilnum,idv), stdv(ilnum,idv),& 
                pave(ilnum,idv)-flat(ilnum,idv), &
                pave(ilnum,idv)+flat(ilnum,idv) )
              gtall(ipop,idv)=gtype
            enddo
          endif
          if(idebug.eq.1)write(mtdeb,*)'ARC: ',ip,iarc,ipop
          ip=ip+1
          nsum(iarc)=nsum(iarc)+1
        endif
      enddo


      narcin=ined-inst+1

      return
      end subroutine DrawArchive

       
! *****************************************************************
      subroutine SelectArchive(ilnum,irkfr,irkto,ied,id,nselect)
! *****************************************************************
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      INTEGER :: id(mtpop)
      INTEGER :: irkfr, irkto,nselect,ied,ilnum,ip,ist 
      ist = 1

      do ip=ist,ied
        if( (irank(ilnum,ip).ge.irkfr) .and.&
           (irank(ilnum,ip).le.irkto) ) then
          nselect = nselect + 1
          id(nselect) = ip
        endif
      enddo
          
      return
      end subroutine SelectArchive


! *****************************************************************
      subroutine SelectViolation(ilnum,idra,nselect)
! *****************************************************************
! -------------------- FIX MUST BE REQUIRED ------------------
! -------------------- FIX MUST BE REQUIRED ------------------
! -------------------- FIX MUST BE REQUIRED ------------------
! -------------------- FIX MUST BE REQUIRED ------------------
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      INTEGER :: idra(mtpop)
      REAL(DP) :: gvmin, gv 
      INTEGER  :: iselect,nselect,ilnum,icheck,icn,ist,ied,ipop,ip 
      ist=1
      ied=istat(ilnum,13)

      do icn=1,ncon
        icheck=0
        gvmin = gvalallmax(ilnum,icn)
!        do ip=ist+1,ied
        do ip=ist,ied
          ipop=idarc(ilnum,ip)
          gv=gvalallmax(ilnum,icn)
          if((gv.gt.0.0_dp).and.(gv.le.gvmin))then   !20140620-2
            gvmin = gv
            iselect = ipop
          endif
        enddo
        icheck=nselect
        do while(icheck.gt.0)
          if(idra(icheck).eq.iselect)icheck=icheck-nselect
          icheck=icheck-1
        enddo
        if(icheck.eq.0)then
          nselect=nselect+1
          idra(nselect)=iselect
        endif
      enddo

      return
      end subroutine SelectViolation

! *******************************************************
      subroutine ArcParent(ilnum)
! *******************************************************
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6       USE messy_airtraf_tools_ga_function,     ONLY: ran1
      IMPLICIT NONE
      INTEGER :: idra(mtpop), id(mtpop)
      INTEGER :: ilnum,i,iarc,iarcst,iarced,ifp,inst,ined,inew,inum,ip 
      INTEGER :: ipop,irkst,irked,ism,narcin,nselect
      nselect=0
      iarcst=1
      iarced=istat(ilnum,13)-istat(ilnum,11)

      if(iarced.le.0)return

!      narcin = int(arcratio*istat(ilnum,11))
      narcin = int(arcratio*istat(ilnum,12))
      if(narcin.lt.1)narcin=1

      irkst=1
      irked=0
      do while(nselect.lt.narcin)
        irked=irked+1
        nselect=0
        call SelectArchive(ilnum,irkst,irked,iarced,idra,nselect)
      enddo
  
      call DuplRemove(idra,nselect)
      call UnfRemove(idra,nselect)

!      if(nselect.lt.narcin)narcin=nselect 

      do i=1,narcin
        iarc = int(ran1(idum)*DBLE(nselect))+1   !HY20140414-5
        if(iarc.gt.nselect)iarc=nselect   
        id(i)=idra(iarc)
      enddo

      call DuplRemove(id,narcin)

      inst=istat(ilnum,12)-narcin+1
      ined=istat(ilnum,12)

      inum=0
      do ip=inst,ined
        inum=inum+1
        idra(inum)=ifitpop(ilnum,ip)
      enddo
      do ip=1,narcin
        inum=inum+1
        idra(inum)=id(ip)
      enddo
      nselect=inum 

      call DuplRemove(idra,nselect)

      inew=ined
      ip=nselect
      do while (inew.ge.inst)
        ipop=idra(ip)
        ism=0
        do ifp=1,inst-1
          if(ifitpop(ilnum,ifp).eq.ipop)ism=ism+1
        enddo
        if(ism.eq.0)then
          ifitpop(ilnum,inew)=ipop
          inew=inew-1
        endif
        ip=ip-1
      enddo    

!      call DrawArchive(ilnum,inst,ined,iarcst,iarced,narcin)

!      if(narcin.gt.0)then
!        fmax=ffit(ilnum,ined)
!        do ip=1,inst-1
!          if(ffit(ilnum,ip).gt.fmax)fmax=ffit(ilnum,ip)
!        enddo
!        do ip=inst,ined
!          ffit(ilnum,ip)=fmax
!        enddo
!      endif

      return
      end subroutine ArcParent

! *******************************************************
      subroutine ArcNBest(ilnum)
! *******************************************************
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6       USE messy_airtraf_tools_ga_function,      ONLY: ran1
      IMPLICIT NONE
      INTEGER :: idra(mtpop),id(mtpop)
      INTEGER :: imax,ilnum,i,iarc,iarcst,iarced,ifp,inst,ined,inew,inum,ip,ipop
      INTEGER :: irkst,irked,ism,narcin,nmax,nselect
      nselect=0
      iarcst=1
      iarced=istat(ilnum,13)
      !write(*,*)'iarced',iarced,'istat(ilnum,11)',istat(ilnum,11)
      if(iarced-istat(ilnum,11).le.0)return
      imax=istat(ilnum,11)*ibestn
      nmax=istat(ilnum,11)*(ibestn+1)

!      narcin = int((istat(ilnum,12)-istat(ilnum,11))*xratio)
      narcin = int(imax*xratio)
      if(narcin.lt.1)narcin=1

      irkst=1
      irked=0
      do while(nselect.lt.narcin)
        irked=irked+1
        nselect=0
        call SelectArchive(ilnum,irkst,irked,iarced,idra,nselect)
      enddo
     
      call DuplRemove(idra,nselect)

      call UnfRemove(idra,nselect)   

      do i=1,narcin
        iarc = int(ran1(idum)*DBLE(nselect))+1   !HY20140414-5
        if(iarc.gt.nselect)iarc=nselect
        id(i)=idra(iarc)
      enddo
 
      call DuplRemove(id,narcin)

! ------ THINK! ------
      if(istat(ilnum,12).lt.nmax)then
        inst=istat(ilnum,12)-narcin+1
        ined=istat(ilnum,12)
      else
        inst=istat(ilnum,12)-istat(ilnum,11)-narcin+1
        ined=istat(ilnum,12)-istat(ilnum,11)
      endif
! ------ THINK! ------

      inum=0
      do ip=inst,ined
        inum=inum+1
        idra(inum)=ifitpop(ilnum,ip)
      enddo
      do ip=1,narcin
        inum=inum+1
        idra(inum)=id(ip)
      enddo
      nselect=inum
 
      call DuplRemove(idra,nselect)

      inew=ined
      ip=nselect
      do while (inew.ge.inst)
        !write(*,*)'idra_check',idra(ip),'ip=',ip,'inew=',inew,'inst=',inst   !HY20140606-2
        ipop=idra(ip)
        ism=0
        do ifp=1,inst-1
          if(ifitpop(ilnum,ifp).eq.ipop)ism=ism+1
        enddo
        if(ism.eq.0)then
          ifitpop(ilnum,inew)=ipop
          inew=inew-1
        endif
        ip=ip-1
      enddo


!      call DrawArchive(ilnum,inst,ined,iarcst,iarced,narcin)
!      print *,'ARC1',narcin

!HY20140328-1      if(idebug.eq.1)call Debug_Obj(ilnum,istat(ilnum,12),mtdeb)  

      return
      end subroutine ArcNBest
      
! *******************************************************
      subroutine UpdateArcInfo(ilnum,idra,nselect)
! *******************************************************
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6       USE messy_airtraf_tools_ga_adapt,     ONLY: EncodeCheck
      IMPLICIT NONE
      INTEGER :: idra(mtpop)
      INTEGER :: nselect,ilnum,idv,ist,ied,ip,ipop
      ist=1
      ied=istat(ilnum,13)
      do ipop=ist,ied
        ipinfo(ipop,4)=-1
        do idv=1,ndv
          gtall(ipop,idv)=0.0_dp   !20140620-2
        enddo
      enddo
      ist=1
      ied=nselect
      do ip=ist,ied
        ipop=idra(ip)
        call EncodeCheck(ilnum,ipop,ipop,1,ndv,1)
!        ipinfo(ipop,4)=1
      enddo

      return
      end subroutine UpdateArcInfo

! *******************************************************
      subroutine ArcCheck(ilinf)
! *******************************************************
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      INTEGER :: ilinf,idv,ist,ied,ilnum,info,ipop,numpop,ip
      if(ilinf.eq.0)then
        ist=1
        ied=nisland
      else
        ist=ilinf
        ied=ilinf
      endif

      do ilnum=ist,ied
        numpop=istat(ilnum,13)
        do ip=1,numpop
          ipop=idarc(ip,ilnum)
          info=ipinfo(ipop,4)
          if(info.eq.2)then
            idv = 1
            do while(idv.le.ndv)
              if((gtall(ipop,idv).lt.gelb(ilnum,idv)).or.&
                (gtall(ipop,idv).gt.geub(ilnum,idv))) &
                idv=idv+ndv
              idv=idv+1
            enddo
            idv=idv-1
            if(idv.eq.ndv)ipinfo(ipop,4)=1
          endif
        enddo
      enddo
   
      return
      end subroutine ArcCheck


! *******************************************************
      subroutine DuplRemove(id,num)
! *******************************************************
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      INTEGER :: id(mtpop)
      INTEGER :: itmp(mtpop),nsum(mtpop)
      INTEGER :: num,i,inum,j
      do i=1,num
        itmp(i)=id(i)
        nsum(i)=0
      enddo

      do i=1,num-1
        do j=i+1,num
          if(id(i).eq.id(j))nsum(i)=nsum(i)+1
        enddo
      enddo

      inum=0
      do i=1,num
        if(nsum(i).eq.0)then
          inum=inum+1
          id(inum)=itmp(i)
        endif
      enddo

      num=inum
      !write(*,*)'num_check_dupl',num   !HY20140606-4
      return
      end subroutine DuplRemove

! *******************************************************
      subroutine UnfRemove(id,num)
! *******************************************************
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6       USE messy_airtraf_tools_ga_adapt,     ONLY: EncodeCheck
      IMPLICIT NONE
      INTEGER :: id(mtpop)
      INTEGER :: itmp(mtpop)
      INTEGER :: num,i,ip,ilnum,inum,ipop
      do i=1,num 
        itmp(i)=id(i)
      enddo

      inum=0
      do ip=1,num
        ipop=id(ip)
        if((ipinfo(ipop,4).eq.-1).or.(ipinfo(ipop,4).eq.2))then
          ilnum=ipinfo(ipop,3)
          call EncodeCheck(ilnum,ipop,ipop,1,ndv,1)
        endif
        if(ipinfo(ipop,4).ne.3)then
          inum=inum+1
          id(inum)=ipop
        endif
      enddo

      num=inum

      return
      end subroutine UnfRemove

!------------------------------------------------------------
!
!
! ORIGINAL: messy_airtraf_tools_ga_afit.f90
!
!
! **************************************************
      subroutine FitnessValue(ilnum,nfitpop,fdec,mtyp)
! ***************************************************
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6       USE messy_airtraf_tools_ga_update,     ONLY: MaxminObj, MaxminGval
      IMPLICIT NONE
      REAL(DP) :: fdec(maxpop,mobj),gdec(maxpop,mcon)
      INTEGER :: murank(maxpop)
      INTEGER :: mtyp,nfitpop,ilnum,icnr
 
      call MaxminObj(ilnum,nfitpop,fdec)
!HY20140328-1      if(idebug.eq.1)call Debug_Obj(ilnum,nfitpop,mtdeb)

      if(ncon.gt.0)then
        call MaxminGval(ilnum,nfitpop,gdec)
        call ConstrHandling(ilnum,nfitpop,gdec,icnr)
        if(mprank.eq.1)&
         call ConstrParetoRanking(ilnum,nfitpop,fdec,gdec,murank,icnr)
        if(mprank.eq.2)&
         call ConstrFrontRanking(ilnum,nfitpop,fdec,gdec,murank,icnr)
!HY20140328-1        if(idebug.eq.1)call Debug_ConstrRank(ilnum,nfitpop,gdec,mtdeb)
      else
        if(mprank.eq.1)&
         call ParetoRanking(ilnum,nfitpop,fdec,murank)
        if(mprank.eq.2)&
         call FrontRanking(ilnum,nfitpop,fdec,murank)
      endif

      if(mtyp.eq.1)call RankBasedFitness(ilnum,nfitpop)
      if(mtyp.eq.2)call AverageFitness(ilnum,nfitpop,murank) 
      if(mtyp.eq.3)call RankPowFitness(ilnum,nfitpop) 

!HY20140328-1      if(idebug.eq.1)call Debug_Fitness(ilnum,nfitpop,mtdeb) 
 
      return
      end subroutine FitnessValue

! **************************************************
      subroutine SFitnessValue(ilnum,nfitpop,fdec)
! ***************************************************
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      REAL(DP) :: fdec(maxpop,mobj)
      REAL(DP) :: sshare
      INTEGER  :: numspop,ilnum,nfitpop
!      dimension gdec(maxpop,mobj)

      numspop=nfitpop
      if(isig2n.eq.0)numspop=istat(ilnum,11)

      call NormSigmaShare(numspop,nobj,sshare)

      if(msdist.eq.1)call NormNicheCount(ilnum,nfitpop,fdec,sshare)
      if(msdist.eq.2)call NormNicheCountDV(ilnum,nfitpop,sshare,0)
      if(msdist.eq.3)call NormNicheCountDV(ilnum,nfitpop,sshare,1)
      if(msdist.eq.4)call NormNicheCountDV(ilnum,nfitpop,sshare,2)

      call SharedFitness(ilnum,nfitpop)

!HY20140328-1      if(idebug.eq.1)call Debug_Niche(ilnum,nfitpop,sshare,mtdeb) 

      return
      end subroutine SFitnessValue



! **************************************************
      subroutine PopSorting(ilnum,nfitpop,iorder)
! ***************************************************
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      INTEGER :: iorder(maxpop)
      INTEGER :: ilnum,nfitpop
!      nfitpop = istat(ilnum,12)
      call OrderByFit(ilnum,nfitpop,iorder)
 
!HY20140328-1      if(idebug.eq.1)call Debug_Order(nfitpop,iorder,mtdeb)
 
      if(ibestn.gt.0)then
         if((irksp.eq.1).or.(ictrlelt.eq.1))&
           call OrderByRankwF(ilnum,nfitpop,iorder)
         if(ictrlelt.eq.1)&
           call ControlElitist(ilnum,nfitpop,iorder,npop)
      endif
 
      call SortByOrder(ilnum,nfitpop,iorder)   

      return
      end subroutine PopSorting

! *******************************************************
      subroutine PreRank(ilnum,numpop,istobj,iedobj,&
                              nrk,murank,bdval)
! *******************************************************
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      INTEGER :: murank(maxpop),nrk(maxpop)
      REAL(DP) :: bdval(mobj)
      INTEGER :: istobj,iedobj,ilnum,numpop,iof,ip
      do ip=1,numpop
        nrk(ip)    = 1
        murank(ip) = 0
      enddo

      do iof=istobj,iedobj
        if(ialpnorm.eq.1)bdval(iof)&
         = objallmax(ilnum,iof) - objallmin(ilnum,iof)
        if(ialpnorm.eq.2)bdval(iof)&
         = objmax(ilnum,iof) - objmin(ilnum,iof)
        if((bdval(iof).lt.bv).or.(ialpnorm.eq.0))bdval(iof) =  1.0_dp   !20140620-2
        if(imaxmin.lt.0) bdval(iof) = -1.0_dp * bdval(iof)   !20140620-2
      enddo   

      return
      end subroutine PreRank

! *******************************************************
      subroutine ParetoRanking(ilnum,numpop,fdec,murank)
! *******************************************************
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      REAL(DP) :: fdec(maxpop,mobj)
      INTEGER :: murank(maxpop), nrk(maxpop)
      REAL(DP) :: fj(mobj),bdval(mobj)
      REAL(DP) :: gj,gs
      INTEGER  :: ilnum,numpop,iedobj,istobj,iof,iofi,iofj,ip,is
!      numpop = istat(ilnum,12)
      istobj = istat(ilnum,4)
      iedobj = istat(ilnum,5)

      call PreRank(ilnum,numpop,istobj,iedobj,nrk,murank,bdval)

      do ip = 1, numpop
        do is = 1, numpop
          gs = 0.0_dp   !20140620-2
          do iof=istobj,iedobj
            fj(iof)=(fdec(ip,iof)-fdec(is,iof))/bdval(iof)
          enddo
          do iofi=istobj,iedobj
            gj = 0.0_dp   !20140620-2
            do iofj=istobj,iedobj
              gj = gj + alpdm(iofi,iofj)*fj(iofj)
            enddo
            if(gj.gt.0.0_dp)goto 10   !20140620-2
            gs = gs + gj 
          enddo
          if(gs.ge.0.0_dp)goto 10   !20140620-2
          nrk(ip) = nrk(ip) + 1 
 10       continue
        enddo
        ifrank(ilnum,ip) = nrk(ip)
      enddo

      do ip=1,numpop
        murank(ifrank(ilnum,ip)) = murank(ifrank(ilnum,ip)) + 1
      enddo

      return
      end subroutine ParetoRanking

! ******************************************************
      subroutine FrontRanking(ilnum,numpop,fdec,murank)
! ******************************************************
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      REAL(DP) :: fdec(maxpop,mobj)
      INTEGER :: murank(maxpop),nrk(maxpop)
      REAL(DP) :: fj(mobj),bdval(mobj)
      REAL(DP) :: gj,gs
      INTEGER  :: ilnum,numpop,iedobj,istobj,ifinish,inum,iof,iofi,iofj,ip,is
!      numpop = istat(ilnum,12)
      istobj = istat(ilnum,4)
      iedobj = istat(ilnum,5)

      inum = 1
      ifinish = 0

      call PreRank(ilnum,numpop,istobj,iedobj,nrk,murank,bdval)

      do while(inum.le.numpop)

        do ip = 1, numpop
          if(nrk(ip).lt.0)goto 11
          nrk(ip) = 1
          do is = 1, numpop
            if(nrk(is).lt.0)goto 10
            gs = 0.0_dp   !20140620-2
            do iof=istobj,iedobj
              fj(iof)=(fdec(ip,iof)-fdec(is,iof))/bdval(iof)
            enddo
            do iofi=istobj,iedobj
              gj = 0.0_dp   !20140620-2
              do iofj=istobj,iedobj
                gj = gj + alpdm(iofi,iofj)*fj(iofj)
              enddo
              if(gj.gt.0.0_dp)goto 10    !20140620-2
              gs = gs + gj
            enddo
            if(gs.ge.0.0_dp)goto 10    !20140620-2
            nrk(ip) = nrk(ip) + 1
 10         continue
          enddo
 11       continue
        enddo

        ifinish = ifinish + 1
        do ip=1,numpop
          if(nrk(ip).eq.1)then
            nrk(ip) = -1
            inum = inum + 1
            ifrank(ilnum,ip) = ifinish
            murank(ifinish) = murank(ifinish) + 1 
          endif
        enddo

      enddo

      return
      end subroutine FrontRanking

! ****************************************************
      subroutine AverageFitness(ilnum,numpop,murank)
! ****************************************************
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      INTEGER :: murank(maxpop)
      INTEGER :: musum,nrk,ilnum,numpop,ip,k
!      numpop = istat(ilnum,12)

      do ip=1,numpop
        musum = 0
        nrk   = ifrank(ilnum,ip)
        do k=1,nrk-1
          musum = musum + murank(k)
        enddo
        ffit(ilnum,ip) = DBLE(numpop-musum)-0.5_dp*DBLE(murank(nrk)-1)  !HY20140414-5
      enddo

      return
      end subroutine AverageFitness

! ****************************************************
      subroutine RankBasedFitness(ilnum,numpop)
! ****************************************************
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      INTEGER :: nrk,ilnum,numpop,ip
!      numpop=istat(ilnum,12)

      do ip=1,numpop
        nrk   = ifrank(ilnum,ip)
        ffit(ilnum,ip) = 1.0_dp/DBLE(nrk)   !HY20140414-5
      enddo

      return
      end subroutine RankBasedFitness

! ****************************************************
      subroutine RankPowFitness(ilnum,numpop)
! ****************************************************
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      INTEGER :: nrk,ilnum,numpop,ip
      do ip=1,numpop
        nrk   = ifrank(ilnum,ip)
        ffit(ilnum,ip) = clr*(1.0_dp-clr)**(nrk-1)   !20140620-2
      enddo

      return
      end subroutine RankPowFitness


! ***********************************************
      subroutine NormSigmaShare(n,m,x)
! ***********************************************
! #### Normalized Sigma Share 
! #### [Ref.] K., Deb., "Mulit-Objective Optimization using 
! ####            Evolutionary Algorithms," Wiley, 2001, pp.197-199.
! #### n ... Number of Population
! #### m ... Number of Objective Functions
! #### x ... Normalized Sigma Share 
! #### Sigma Share is calculatedd by Newton-Raphson Algorithm
      IMPLICIT NONE
      REAL(DP) :: x,g,f 
      INTEGER  :: i,n,m
 
      if(n.lt.2)return

      i=1
!     i=0
      g=0.0_dp   !20140620-2
      x=0.0_dp   !20140620-2

      if(m.eq.2)then
        x = 2.0_dp / DBLE(n-1)   !HY20140414-6,20140620-2
        return
      endif

      do while(x.le.1.e-4)
        x = 2.0_dp / DBLE(n-1) * 10.0_dp**i   !HY20140414-6,20140620-2
        f = 1.0_dp   !20140620-2  
        do while(abs(f).gt.1.e-3) 
          f = (1.0_dp+x)**m -DBLE(n)*x**m -1.0_dp   !HY20140414-6,20140620-2
          g = DBLE(m)*(1.0_dp+x)**(m-1) -DBLE(n)*m*x**(m-1) -1.0_dp   !HY20140414-6,20140620-2
          x = x - f/g 
        enddo
        i=i+1
      enddo

      return
      end subroutine NormSigmaShare

! ****************************************************
      subroutine NormNicheCount(ilnum,numpop,fdec,sshare)
! **************************************************** 
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      REAL(DP) :: fdec(maxpop,mobj)
      REAL(DP) :: objd(mobj) 
      REAL(DP) :: objn,sshare,dijnorm,dini
      INTEGER  :: ilnum,numpop,istobj,iedobj,iof,ip,is,nrkp,nrks
!      numpop = istat(ilnum,12)
      istobj = istat(ilnum,4)
      iedobj = istat(ilnum,5)

      dini = 0.0_dp   !20140620-2
      do iof=1,nobj
        objd(iof) = objmax(ilnum,iof) - objmin(ilnum,iof)
        if(objd(iof).lt.bv)then
          objd(iof) = 1.0_dp   !20140620-2
          dini = 1.0_dp   !20140620-2
        endif
      enddo  

      do ip = 1, numpop 
        nrkp = ifrank(ilnum,ip)
        sniche(ilnum,ip) = 0.0_dp   !20140620-2
        do is = 1, numpop
          nrks = ifrank(ilnum,is)
!          dijnorm = 0.
          dijnorm = dini
          if((nrkp.eq.nrks).or.(iallsh.eq.1))then
            do iof=istobj,iedobj
              objn    = fdec(ip,iof) - fdec(is,iof)
              dijnorm = dijnorm + (objn/objd(iof))**2 
            enddo
            dijnorm = sqrt(dijnorm)
            sniche(ilnum,ip) = sniche(ilnum,ip) + & 
            ( 1.0_dp - (min(dijnorm/sshare,1.0_dp))**shalpha )   !20140620-2
          endif
        enddo
        if(sniche(ilnum,ip).lt.1.0_dp)sniche(ilnum,ip)=1.0_dp   !20140620-2
      enddo
 
      return
      end subroutine NormNicheCount

! ****************************************************
      subroutine NormNicheCountDV(ilnum,numpop,sshare,inorm)
! **************************************************** 
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      REAL(DP) :: dvd(mdv) 
      REAL(DP) :: dvmax(mdv),dvmin(mdv)
      REAL(DP) :: dvn,sshare,dijnorm,dini,value
      INTEGER  :: inorm,ilnum,numpop,idv,ip,ipop,is,isub,nrkp,nrks
      if(inorm.eq.2)then
        do idv=1,ndv
          dvmax(idv)=pave(ilnum,idv)+flat(ilnum,idv)
          dvmin(idv)=pave(ilnum,idv)-flat(ilnum,idv)
        enddo
      else
        do idv=1,ndv
          dvmax(idv)=phub(idv)
          dvmin(idv)=phlb(idv)
        enddo
        if(inorm.eq.0)then
          do idv=1,ndv
            value      = dvmax(idv)
            dvmax(idv) = dvmin(idv)
            dvmin(idv) = value
            do ip=1,numpop
              ipop  = ifitpop(ilnum,ip)
              value = ptall(ipop,idv)  
              if(value.gt.dvmax(idv))dvmax(idv)=value
              if(value.lt.dvmin(idv))dvmin(idv)=value
            enddo
          enddo
        endif
      endif

      dini = 0.0_dp   !20140620-3
      do idv=1,ndv
        dvd(idv) = dvmax(idv)-dvmin(idv)
        if(dvd(idv).lt.bv)then
          dvd(idv) = 1.0_dp   !20140620-3
          dini = 1.0_dp   !20140620-3
        endif
      enddo  

      do ip = 1, numpop 
        sniche(ilnum,ip) = 0.0_dp   !20140620-3
        nrkp = ifrank(ilnum,ip)
        ipop=ifitpop(ilnum,ip)
        do is = 1, numpop
          nrks = ifrank(ilnum,is)
          isub=ifitpop(ilnum,is)
!          dijnorm = 0.
          dijnorm = dini
          if((nrkp.eq.nrks).or.(iallsh.eq.1))then
            do idv=1,ndv
              dvn = ptall(ipop,idv) - ptall(isub,idv)
              dijnorm = dijnorm + (dvn/dvd(idv))**2 
            enddo
            dijnorm = sqrt(dijnorm)
            sniche(ilnum,ip) = sniche(ilnum,ip) + & 
            ( 1.0_dp - (min(dijnorm/sshare,1.0_dp))**shalpha )   !20140620-3
          endif
        enddo
        if(sniche(ilnum,ip).lt.1.0_dp)sniche(ilnum,ip)=1.0_dp   !20140620-3
      enddo
 
      return
      end subroutine NormNicheCountDV

! **********************************************
      subroutine SharedFitness(ilnum,numpop)
! **********************************************
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      INTEGER :: ilnum,numpop,ip
!      numpop = istat(ilnum,12)

      

      do ip=1,numpop
        ffit(ilnum,ip) = ffit(ilnum,ip) / sniche(ilnum,ip)
      enddo

      return
      end subroutine SharedFitness

! ***********************************************
      subroutine OrderByFit(ilnum,numpop,iorder)
! ***********************************************
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      REAL(DP) :: fit(maxpop)
      INTEGER :: iorder(maxpop), ir(maxpop)
      REAL(DP) :: ftmp
      INTEGER  :: itmp,ilnum,numpop,i1,i2,ip
!      numpop = istat(ilnum,12)

      do ip=1,numpop
        fit(ip) = ffit(ilnum,ip)
        ir(ip) = ip
      enddo

      do i1 = 1, numpop-1
        do i2 = i1+1, numpop
          if(fit(i1).lt.fit(i2))then
             ftmp    = fit(i1)
             fit(i1) = fit(i2)
             fit(i2) = ftmp
             itmp    = ir(i1)
             ir(i1)  = ir(i2)
             ir(i2)  = itmp
          endif
        enddo
      enddo

      do ip=1,numpop
        iorder(ir(ip)) = ip
      enddo

      return
      end subroutine OrderByFit


! ***********************************************
      subroutine SortByOrder(ilnum,numpop,iorder)
! ***********************************************
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      INTEGER :: iorder(maxpop)
      REAL(DP) :: ftmp
      INTEGER  :: itmp,ilnum,numpop,i1,i2
!      numpop=istat(ilnum,12)

      do i1=1,numpop-1
        do i2=i1+1,numpop
          if(iorder(i2).lt.iorder(i1))then
            itmp = ifitpop(ilnum,i1)
            ifitpop(ilnum,i1) = ifitpop(ilnum,i2)
            ifitpop(ilnum,i2) = itmp
            itmp = iorder(i1)
            iorder(i1) = iorder(i2)
            iorder(i2) = itmp
            ftmp = ffit(ilnum,i1)
            ffit(ilnum,i1)=ffit(ilnum,i2)
            ffit(ilnum,i2)=ftmp
          endif
        enddo
      enddo

      return
      end subroutine SortByOrder

! ********************************************************
      subroutine ControlElitist(ilnum,numpop,iorder,nmax)
! ********************************************************
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      INTEGER :: iorder(maxpop) 
      INTEGER :: ir(maxpop),murank(maxpop),munum(maxpop)
      REAL(DP) :: ratio
      INTEGER  :: irkend, nro, ni, nit, jnum, nrk,nmax,ilnum,numpop,i,inum,ipop,irk
!      numpop = istat(ilnum,12)

      inum = 0
      irk = 1
      irkend = 1
      nro = 0
   
      do i=1,numpop
        murank(i) = 0
        murank(ifrank(ilnum,i)) = murank(ifrank(ilnum,i)) + 1
        if(ifrank(ilnum,i).gt.irkend)irkend=ifrank(ilnum,i)
      enddo
      do i=1,irkend
        munum(i) = murank(i)
      enddo

      do while(inum.lt.nmax)
        ratio=(1.0_dp-rkcoef)/(1.0_dp-rkcoef**irkend)*rkcoef**(irk-1)   !20140620-3
        ni  = nint(nmax*ratio) + nro
        nit = murank(irk)
        if(nit.gt.ni)then
          if((inum+ni).gt.nmax)then
            ni = nmax - inum
            murank(irk) = murank(irk) - ni
            nro  = 0
            inum = inum + ni
            if(idebug.eq.1)write(mtdeb,*)'case1',irk,inum,ni 
          else
            murank(irk) = murank(irk) - ni
            nro  = 0
            inum = inum + ni
            if(idebug.eq.1)write(mtdeb,*)'case2',irk,inum,ni 
          endif
        else
          if((inum+nit).gt.nmax)then
            ni = nmax - inum
            murank(irk) = murank(irk) - ni
            nro = 0
            inum = inum + ni
            if(idebug.eq.1)write(mtdeb,*)'case3',irk,inum,ni
          else
            murank(irk) = murank(irk) - nit
            nro = ni - nit
            inum = inum + nit
            if(idebug.eq.1)write(mtdeb,*)'case4',irk,inum,nit
          endif
        endif
        irk = irk + 1
        if(irk.gt.irkend)irk = 1
      enddo

      do i=1,irkend
        munum(i) = munum(i) - murank(i)
        if(idebug.eq.1)write(mtdeb,*)i,munum(i),murank(i)
      enddo

      do i=1,numpop 
        ir(iorder(i))=i
      enddo

      inum = 0
      jnum = nmax
      do ipop=1,numpop
        nrk = ifrank(ilnum,ipop) 
        if(munum(nrk).gt.0)then
          inum = inum + 1
          iorder(ir(ipop)) = inum 
          munum(nrk) = munum(nrk) - 1
        else
          jnum = jnum + 1
          iorder(ir(ipop)) = jnum
        endif  
      enddo

      return
      end subroutine ControlElitist

! ************************************************
      subroutine OrderByRankwF(ilnum,numpop,iorder)
! ************************************************
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      INTEGER :: iorder(maxpop)
      INTEGER :: nrk(maxpop),ir(maxpop)
      INTEGER :: itmp,ilnum,numpop,i1,i2,ip,ntmp
!      numpop = istat(ilnum,12)

      do ip=1,numpop
        ir(iorder(ip))=ip
        nrk(iorder(ip))=ifrank(ilnum,ip)
      enddo

      do i1=1,numpop-1
        do i2=i1+1,numpop
           if(nrk(i1).gt.nrk(i2))then
             ntmp    = nrk(i1)
             nrk(i1) = nrk(i2)
             nrk(i2) = ntmp
             itmp    = ir(i1)
             ir(i1)  = ir(i2)
             ir(i2)  = itmp
           endif
         enddo
       enddo

      do ip=1,numpop
        iorder(ir(ip)) = ip
      enddo

      return
      end subroutine OrderByRankwF

! ************************************************************
      subroutine ConstrParetoRanking(ilnum,numpop,fdec,gdec,murank,icnr)
! ************************************************************
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      REAL(DP) :: fdec(maxpop,mobj),gdec(maxpop,mobj)
      INTEGER :: murank(maxpop),nrk(maxpop)
      REAL(DP) :: bdval(mobj),fj(mobj)
      REAL(DP) :: gj,gs,gpop,gsub
      INTEGER  :: ilnum,numpop,icnr,icn,istcon,iedcon,istobj,iedobj,ifinish,igp,igs,inum,iof,iofi,iofj
      INTEGER  :: ipop,isub
!      numpop = istat(ilnum,12)
      istobj = istat(ilnum,4)
      iedobj = istat(ilnum,5)
      istcon = istat(ilnum,6)
      iedcon = istat(ilnum,7)

      call PreRank(ilnum,numpop,istobj,iedobj,nrk,murank,bdval)
      inum = 0
      ifinish = 0

      do ipop = 1, numpop
        do isub = 1, numpop
          gs = 0.0_dp   !20140620-3
          igp = 0
          igs = 0
          do icn=istcon,iedcon
            gpop = max(0.0_dp,gdec(ipop,icn))   !20140620-3
            gsub = max(0.0_dp,gdec(isub,icn))   !20140620-3
            if(gpop.gt.0.0_dp)igp=igp+1   !20140620-3
            if(gsub.gt.0.0_dp)igs=igs+1   !20140620-3
            if((icnr.eq.1).and.(gpop.gt.gsub).and.(gsub.gt.0.0_dp))igp=igp+1   !20140620-3
          enddo
          if((igp.gt.0).or.(igs.gt.0))then
            if(igp.gt.igs)nrk(ipop)=nrk(ipop)+1
            goto 10
          else
            gs = 0.0_dp   !20140620-3
          endif
          do iof=istobj,iedobj
            fj(iof)=(fdec(ipop,iof)-fdec(isub,iof))/bdval(iof)
          enddo
          do iofi=istobj,iedobj
            gj = 0.0_dp   !20140620-3
            do iofj=istobj,iedobj
              gj = gj + alpdm(iofi,iofj)*fj(iofj)
            enddo
            if(gj.gt.0.0_dp)goto 10   !20140620-3
            gs = gs + gj
          enddo
          if(gs.ge.0.0_dp)goto 10   !20140620-3
          nrk(ipop) = nrk(ipop) + 1
 10       continue
        enddo
        ifrank(ilnum,ipop) = nrk(ipop)
      enddo

      do ipop=1,numpop
        murank(ifrank(ilnum,ipop)) = murank(ifrank(ilnum,ipop))+1
      enddo

      return
      end subroutine ConstrParetoRanking

! *********************************************************
      subroutine ConstrFrontRanking(ilnum,numpop,fdec,gdec,murank,icnr)
! *********************************************************
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      REAL(DP) :: fdec(maxpop,mobj),gdec(maxpop,mcon)
      INTEGER :: murank(maxpop),nrk(maxpop)
      REAL(DP) :: bdval(mobj),fj(mobj)
      REAL(DP) :: gj,gpop,gs,gsub
      INTEGER  :: ilnum,numpop,icnr,icn,iedcon,iedobj,ifinish,igp,igs,inum,iof,iofi,iofj,ipop
      INTEGER  :: istcon, istobj,isub
!      numpop = istat(ilnum,12)
      istobj = istat(ilnum,4)
      iedobj = istat(ilnum,5)
      istcon = istat(ilnum,6)
      iedcon = istat(ilnum,7)

      inum = 1
      ifinish = 0
      call PreRank(ilnum,numpop,istobj,iedobj,nrk,murank,bdval)

      do while(inum.le.numpop)

        do ipop = 1, numpop
          if(nrk(ipop).lt.0)goto 11
          nrk(ipop) = 1
          do isub = 1, numpop
            if(nrk(isub).lt.0)goto 10
            gs  = 0.0_dp   !20140620-3
            igp = 0
            igs = 0
            do icn=istcon,iedcon
              gpop = max(0.0_dp,gdec(ipop,icn))   !20140620-3
              gsub = max(0.0_dp,gdec(isub,icn))   !20140620-3
              if(gpop.gt.0.0_dp)igp=igp+1   !20140620-3
              if(gsub.gt.0.0_dp)igs=igs+1   !20140620-3
            if((icnr.eq.1).and.(gpop.gt.gsub).and.(gsub.gt.0.0_dp))igp=igp+1   !20140620-3
            enddo
            if((igp.gt.0).or.(igs.gt.0))then
              if(igp.gt.igs)nrk(ipop)=nrk(ipop)+1
              goto 10
            else
              gs = 0.0_dp   !20140620-3
            endif
            do iof=istobj,iedobj
              fj(iof)=(fdec(ipop,iof)-fdec(isub,iof))/bdval(iof)
            enddo
            do iofi=istobj,iedobj
              gj = 0.0_dp   !20140620-3
              do iofj=istobj,iedobj
                gj = gj + alpdm(iofi,iofj)*fj(iofj)
              enddo
              if(gj.gt.0.0_dp)goto 10   !20140620-3
              gs = gs + gj
            enddo
            if(gs.ge.0.0_dp)goto 10   !20140620-3
            nrk(ipop) = nrk(ipop) + 1
 10         continue
          enddo
 11       continue
        enddo

        ifinish = ifinish + 1
        do ipop=1,numpop
          if(nrk(ipop).eq.1)then
            nrk(ipop) = -1
            inum = inum + 1
            ifrank(ilnum,ipop) = ifinish
            murank(ifinish) = murank(ifinish) + 1
          endif
        enddo

      enddo

      return
      end subroutine ConstrFrontRanking

! *************************************************
      subroutine ConstrHandling(ilnum,numpop,gdec,icnr)
! ************************************************
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      REAL(DP) :: gdec(maxpop,mcon)
      INTEGER :: ilnum,numpop,icnr,icn,ip
!      numpop = istat(ilnum,12)
      
      if(icnrrnk.gt.0)then
        if(icnrrnk.eq.1)then
          icnr = 0
        else     
          icnr = 1
        endif
      else
        icnr = 0
        do ip=1,numpop
          do icn=1,ncon
            gdec(ip,icn)=0.0_dp   !20140620-3
          enddo
        enddo
      endif
          
      return
      end subroutine ConstrHandling

!-----------------------------------------------------
!
!
! ORIGINAL: messy_airtraf_tools_ga_selection.f90
!
!
! ************************************************
      subroutine Selection(ilnum)
! ************************************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6       USE messy_airtraf_tools_ga_afit,     ONLY: FitnessValue, SFitnessValue
      IMPLICIT NONE
      REAL(DP) :: fdec(maxpop,mobj)
      INTEGER :: noffpop,ilnum,nfitpop
!      dimension gdec(maxpop,mcon)
!      dimension murank(maxpop)

      noffpop =istat(ilnum,11)
      nfitpop=istat(ilnum,11)
 
      if(mselection.eq.0)then
        call EQS(ilnum,noffpop)
      else
        call FitnessValue(ilnum,nfitpop,fdec,mfitsel)

!       if(ishare.eq.2)then
!         call MaxminObj(ilnum,nfitpop,fdec)
!         if(ncon.gt.0)then
!           call MaxminGval(ilnum,nfitpop,gdec)
!           call ConstrHandling(ilnum,nfitpop,gdec,icnr)
!           if(mprank.eq.1)
!     &       call ConstrParetoRanking(ilnum,nfitpop,fdec,gdec,murank,icnr)
!           if(mprank.eq.2)
!     &       call ConstrFrontRanking(ilnum,nfitpop,fdec,gdec,murank,icnr)
!         else
!           if(mprank.eq.1)
!     &        call ParetoRanking(ilnum,nfitpop,fdec,murank)
!           if(mprank.eq.2)
!     &        call FrontRanking(ilnum,fdec,murank)
!         endif
!         if(mfiteval.eq.0)call RankBasedFitness(ilnum)
!         if(mfiteval.eq.1)call AverageFitness(ilnum,murank)
         
       if(ishare.eq.2)then
         call SFitnessValue(ilnum,nfitpop,fdec)
!         numpop=istat(ilnum,12)
!         if(isig2n.eq.0)numpop=istat(ilnum,11)

!         call NormSigmaShare(numpop,nobj,sshare)
!         if(msdist.eq.1)call NormNicheCount(ilnum,fdec,sshare)
!         if(msdist.eq.2)call NormNicheCountDV(ilnum,0)
!         if(msdist.eq.3)call NormNicheCountDV(ilnum,1)
!         if(msdist.eq.4)call NormNicheCountDV(ilnum,2)  
!         call SharedFitness(ilnum)
       endif

       if(mselection.eq.1)call SUS(ilnum,nfitpop,noffpop)
       if(mselection.eq.2)call RWS(ilnum,nfitpop,noffpop)
      endif
 
      return
      end subroutine Selection

! ***********************************************
      subroutine MatingPool(ilnum,numpop)
! ***********************************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      INTEGER :: ilnum,numpop,inum,ip
      inum = 1
      ip = 1
      do while(ip.le.numpop)
        if(nunit(ilnum,inum).gt.0)then
          idparent(ilnum,ip)=ifitpop(ilnum,inum)
          nunit(ilnum,inum) = nunit(ilnum,inum) -1
          ip = ip + 1
        else
          inum = inum + 1
        endif
      enddo
 
      return
      end subroutine MatingPool

! ****************************************
      subroutine SUS(ilnum,numpop,newpop)
! ****************************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6       USE messy_airtraf_tools_ga_function,     ONLY: ran1
      IMPLICIT NONE
      REAL(DP) :: fitsum(maxpop)
      REAL(DP) :: point,r
      INTEGER  :: ilnum,numpop,newpop,i,inum,ip,num
!      newpop = istat(ilnum,11)

      do ip=1,numpop
        nunit(ilnum,ip) = 0     !Corrected: HY20140414-1
      enddo

      fitsum(1) = ffit(ilnum,1)
      do ip=2,numpop
        fitsum(ip)=fitsum(ip-1) + ffit(ilnum,ip)
      enddo

      do i=1,numpop
        fitsum(i) = fitsum(i) / fitsum(numpop)
      enddo

      num=0
      ip=1
      r=ran1(idum)/DBLE(newpop)   !HY20140414-6
      do inum=1,newpop
  10    point = r + DBLE(inum-1)/DBLE(newpop)   !HY20140414-6
        if(point.le.fitsum(ip))then
          num  = num + 1
          nunit(ilnum,ip) = nunit(ilnum,ip) + 1
        else
          ip = ip + 1
          goto 10
        endif
      enddo
!check
      if(num.ne.newpop)then
        nunit(ilnum,newpop) = nunit(ilnum,newpop) + 1
        write(6,*)'Wrong nunit',num, newpop
      endif

      return
      end subroutine SUS

! ****************************************
      subroutine RWS(ilnum,numpop,newpop)
! ****************************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6       USE messy_airtraf_tools_ga_function,     ONLY: ran1
      IMPLICIT NONE
      REAL(DP) :: fitsum(maxpop)
      REAL(DP) :: point
      INTEGER  :: ilnum,numpop,newpop,i,ip,num 
!      newpop = istat(ilnum,11)

      do i=1,numpop
        nunit(ilnum,i) = 0
      enddo

      fitsum(1) = ffit(ilnum,1)
      do ip=2,numpop
        fitsum(ip)=fitsum(ip-1) + ffit(ilnum,ip)
      enddo
      do i=1,numpop
        fitsum(i) = fitsum(i) / fitsum(numpop)
      enddo   

      do num=1,newpop
  10    point = ran1(idum)
        ip = 1
        do while(point.ge.fitsum(ip)) 
          ip = ip + 1
        enddo
        if(ip.gt.newpop)ip=newpop
        nunit(ilnum,ip) = nunit(ilnum,ip) + 1
      enddo

      return
      end subroutine RWS

! ********************************************
      subroutine EQS(ilnum,numpop)
! ********************************************
! ##### No Selection 
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      INTEGER :: ilnum,numpop,ip
!      numpop = istat(ilnum,11)

      do ip=1,numpop
        nunit(ilnum,ip)=1
      enddo

      return
      end subroutine EQS

! *******************************************************
      subroutine Pair(numpop,mate1,mate2)
! *******************************************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6       USE messy_airtraf_tools_ga_function,     ONLY: ran1
      IMPLICIT NONE
      INTEGER :: numpop,mate1,mate2
  11  mate1 = int( DBLE(numpop)*ran1(idum) ) + 1   !HY20140414-6
      if(mate1.eq.(numpop+1))mate1=numpop

  12  mate2 = int( DBLE(numpop)*ran1(idum) ) + 1   !HY20140414-6
      if(mate2.eq.(numpop+1))mate2=numpop

      if(mate2.eq.mate1)goto 12

      return
      end subroutine Pair

! ************************************************************
      subroutine PairRandom(ilnum,matepair,numpair,numparent)
! ************************************************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      INTEGER :: matepair(maxpop,3)
      INTEGER :: loop,numpair,numparent,ilnum,i,j,inum,m1,m2,mate1,mate2
      loop=1

      inum=1
16    continue
      do while(inum.le.numpair)
        call Pair(numparent,m1,m2)
        mate1=idparent(ilnum,m1)
        mate2=idparent(ilnum,m2)
        if(mate1.ne.mate2)then
          matepair(inum,1)=mate1
          matepair(inum,2)=mate2
         matepair(inum,3)=matepair(inum,1)*nallpop &
             +matepair(inum,2)
          inum=inum+1
        endif
      enddo

      inum=1
      do i=1,numpair-1
        do j=i+1,numpair
          if(matepair(i,3).eq.matepair(j,3))then
            loop=loop+1
            if(loop.lt.10*numparent**2)goto 16
          endif
        enddo
      enddo
        

      return
      end subroutine PairRandom

! ************************************************************
      subroutine PairOrder(ilnum,matepair,numpair,numparent,maxl)
! ************************************************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6       USE messy_airtraf_tools_ga_function,     ONLY: ran1
      IMPLICIT NONE
      INTEGER :: matepair(maxpop,3)
      INTEGER :: matetmp(maxpop)
      INTEGER :: nuni(maxpop)
      INTEGER :: ndipl,n1,mtmp,numpair,numparent,ilnum,maxl,i,j,m1,m2,ntmp
!      maxl = mrandp

      do i=1,numparent
        matetmp(i)=idparent(ilnum,i)
        nuni(i)=1
      enddo

      ndipl=0
      do i=1,numparent-1
        do j=i+1,numparent
          if((nuni(i).eq.1).and.(nuni(j).eq.1))then
            m1=matetmp(i)
            m2=matetmp(j)
            if(m1.eq.m2)then
              nuni(i)=nuni(i)+1
              nuni(j)=0
              if(nuni(i).gt.ndipl)then
                n1=i
                ndipl=nuni(i)
              endif
            endif
          endif
        enddo
      enddo
      ntmp=0
      if(ndipl.gt.numparent/4)then
        do i=1,numparent
          if(nuni(i).gt.0)then
            if((nuni(i).gt.ntmp).and.(i.ne.n1))then
              ntmp=nuni(i)
            endif
          endif
        enddo
        ndipl=ndipl+ntmp
        if(ndipl.le.numparent/2)ndipl=0
      else
        ndipl=0
      endif
          

 17   continue
      do i=1,maxl
!        call Pair(numparent,m1,m2)
        m1=int(DBLE(numparent-1)*ran1(idum))+1   !HY20140414-6
        m2=m1+1
        mtmp = matetmp(m1)
        matetmp(m1) = matetmp(m2) 
        matetmp(m2) = mtmp
      enddo   

      maxl=1
      do i=1,numpair
        m1 = 2*i-1
        m2 = 2*i
        if(matetmp(m1).eq.matetmp(m2))goto 17
        matepair(i,1)=matetmp(m1)
        matepair(i,2)=matetmp(m2)
        matepair(i,3)=matepair(i,1)*nallpop+matepair(i,2)
      enddo

      if(ndipl.eq.0)then
        do i=1,numpair-1
          do j=i+1,numpair
            if(matepair(i,3).eq.matepair(j,3))goto 17
          enddo
        enddo
      endif
        
      return
      end subroutine PairOrder

!-----------------------------------------------------------------
!
!
! ORIGINAL: messy_airtraf_tools_ga_cross.f90
!
!
! *******************************************************
      subroutine Crossover(ilnum,ioff1,ioff2,mate1,mate2)
! *******************************************************
!     mcros   ... Crossover Method
!       0 ARITHMETIC CROSSOVER (FIXED CROSSOVER)
!       1 BLENDED CROSSOVER
!       2 SIMULATED BINARY CROSSOVER (FIXED CROSSOVER)
!       3 SIMULATED BINARY
!       4 NEW SIMULATED BINARY CROSSOVER        
!
!  NEW DEFINITION BY 2004 ---CAUTION---
! *******************************************************
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      INTEGER :: mate1,mate2,ilnum,ioff1,ioff2,method
      method = mcros(ilnum)

      if(method.eq.0)call CrossBLX(ilnum,ioff1,ioff2,mate1,mate2,1)
      if(method.eq.1)call CrossBLX(ilnum,ioff1,ioff2,mate1,mate2,0)
      if(method.eq.2)call CrossSBX(ilnum,ioff1,ioff2,mate1,mate2,1)
      if(method.eq.3)call CrossSBX(ilnum,ioff1,ioff2,mate1,mate2,0)
      if(method.eq.4)call CrossNEWSBX(ilnum,ioff1,ioff2,mate1,mate2,0)

      return
      end subroutine Crossover

! --------------------------------------------------
      subroutine CrossBLX(ilnum,j1,j2,m1,m2,ifix)
! --------------------------------------------------
! ####### IFIX=0 BLENDED CROSSOVER BY ESHELMAN AND SCHAFFER
! ####### IFIX=1 ARITHMETIC CROSSOVER BY MICHALEWICZ AND ..
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      REAL(DP) :: rand0, rand1,alpha,bl,bu,pm1,pm2,pmemp1,pmemp2,rate
      INTEGER  :: ilnum,j1,j2,m1,m2,ifix,idv
!HY20140404-6       USE messy_airtraf_tools_ga_function,    ONLY: ran1
      rate = rcros(ilnum)
      alpha = pcros(ilnum)

      if(ifix.eq.0)then

!  < BLENDED CROSSOVER >
        do idv=1,ndv
          if(ran1(idum).le.rate)then
!cAPR            pmemp1 = pmemp(m1,idv)
!cAPR            pmemp2 = pmemp(m2,idv)
            pmemp1 = gtall(m1,idv)
            pmemp2 = gtall(m2,idv)
            bl     = gelb(ilnum,idv)
            bu     = geub(ilnum,idv)
            
  11        rand0  = (2.0_dp*alpha+1.0_dp)*ran1(idum)-alpha   !20140620-3
            rand1 = 1.0_dp - rand0   !20140620-3

            pm1 = rand0*pmemp1 + rand1*pmemp2
            pm2 = rand1*pmemp1 + rand0*pmemp2

            if(icrobd.eq.1)then
              if(pm1.lt.bl)pm1=bl
              if(pm1.gt.bu)pm1=bu
              if(pm2.lt.bl)pm2=bl
              if(pm2.gt.bu)pm2=bu
            else
              if((pm1.lt.bl).or.(pm1.gt.bu))goto 11
              if((pm2.lt.bl).or.(pm2.gt.bu))goto 11
            endif

!cAPR            pmem(j1,idv)=pm1
!cAPR            pmem(j2,idv)=pm2
            gtall(j1,idv)=pm1
            gtall(j2,idv)=pm2
          else
!cAPR            pmem(j1,idv)=pmemp(m1,idv)
!cAPR            pmem(j2,idv)=pmemp(m2,idv)
            gtall(j1,idv)=gtall(m1,idv)
            gtall(j2,idv)=gtall(m2,idv)
          endif
         enddo
  
      else
      
!  < ARITHMETIC CROSSOVER >
        if(ran1(idum).le.rate)then

  12      rand0 = (2.0_dp*alpha+1.0_dp)*ran1(idum)-alpha   !20140620-3
          rand1 = 1.0_dp - rand0    !20140620-3

          do idv=1,ndv
!cAPR            pmemp1 = pmemp(m1,idv)
!cAPR            pmemp2 = pmemp(m2,idv)
            pmemp1 = gtall(m1,idv)
            pmemp2 = gtall(m2,idv)
            bl     = gelb(ilnum,idv)
            bu     = geub(ilnum,idv)

            pm1 = rand0*pmemp1 + rand1*pmemp2
            pm2 = rand1*pmemp1 + rand0*pmemp2

            if(icrobd.eq.1)then
              if(pm1.lt.bl)pm1=bl
              if(pm1.gt.bu)pm1=bu
              if(pm2.lt.bl)pm2=bl
              if(pm2.gt.bu)pm2=bu
            else
              if((pm1.lt.bl).or.(pm1.gt.bu))goto 12
              if((pm2.lt.bl).or.(pm2.gt.bu))goto 12
            endif
!cAPR            pmem(j1,idv)=pm1
!cAPR            pmem(j2,idv)=pm2
            gtall(j1,idv)=pm1
            gtall(j2,idv)=pm2
          enddo
        else
          do idv=1,ndv
!cAPR            pmem(j1,idv)=pmemp(m1,idv)
!cAPR            pmem(j2,idv)=pmemp(m2,idv)
            gtall(j1,idv)=gtall(m1,idv)
            gtall(j2,idv)=gtall(m2,idv)
          enddo
        endif

      endif

      return 
      end subroutine CrossBLX

! --------------------------------------------------
      subroutine CrossSBX(ilnum,j1,j2,m1,m2,ifix)
! --------------------------------------------------
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      REAL(DP) :: sbxeta,betaq,bl,bu,pm1,pm2,pmemp1,pmemp2,rate,u
      INTEGER  :: ilnum,j1,j2,m1,m2,ifix,idv
!HY20140404-6       USE messy_airtraf_tools_ga_function,     ONLY: ran1
      rate = rcros(ilnum)
      sbxeta = pcros(ilnum)

      if(ifix.eq.0)then

!  < SIMULATED BINARY CROSSOVER >
        do idv=1,ndv
          if(ran1(idum).le.rate)then 
            pmemp1 = gtall(m1,idv)
            pmemp2 = gtall(m2,idv)
            bl     = gelb(ilnum,idv)
            bu     = geub(ilnum,idv)

  11        u = ran1(idum)
            if(u.le.0.5_dp)then   !20140620-4
              betaq = (2.0_dp*u)**(1.0_dp/(sbxeta+1.0_dp))   !20140620-4
            elseif(u.lt.1.0_dp)then   !20140620-4
              betaq = (1.0_dp/(2.0_dp*(1.0_dp-u)))**(1.0_dp/(sbxeta+1.0_dp))   !20140620-4
            else
              goto 11
            endif

            pm1 = 0.5_dp*((1.0_dp+betaq)*pmemp1 + (1.0_dp-betaq)*pmemp2)   !20140620-4
            pm2 = 0.5_dp*((1.0_dp-betaq)*pmemp1 + (1.0_dp+betaq)*pmemp2)   !20140620-4

            if(icrobd.eq.1)then
              if(pm1.lt.bl)pm1=bl
              if(pm1.gt.bu)pm1=bu
              if(pm2.lt.bl)pm2=bl
              if(pm2.gt.bu)pm2=bu
            else
              if((pm1.lt.bl).or.(pm1.gt.bu))goto 11
              if((pm2.lt.bl).or.(pm2.gt.bu))goto 11
            endif

            gtall(j1,idv)=pm1
            gtall(j2,idv)=pm2
          else
            gtall(j1,idv)=gtall(m1,idv)
            gtall(j2,idv)=gtall(m2,idv)
          endif
        enddo

      else

!  < SIMULATED BINARY CROSSOVER (FIXED) >
        if(ran1(idum).le.rate)then
          do idv=1,ndv
            pmemp1 = gtall(m1,idv)
            pmemp2 = gtall(m2,idv)
            bl     = gelb(ilnum,idv)
            bu     = geub(ilnum,idv)

  12        u = ran1(idum)
            if(u.le.0.5_dp)then   !20140620-4
              betaq = (2.0_dp*u)**(1.0_dp/(sbxeta+1.0_dp))   !20140620-4
            elseif(u.lt.1.0_dp)then   !20140620-4
              betaq = (1.0_dp/(2.0_dp*(1.0_dp-u)))**(1.0_dp/(sbxeta+1.0_dp))   !20140620-4
            else
              goto 12
            endif

            pm1 = 0.5_dp*((1.0_dp+betaq)*pmemp1 + (1.0_dp-betaq)*pmemp2)   !20140620-4
            pm2 = 0.5_dp*((1.0_dp-betaq)*pmemp1 + (1.0_dp+betaq)*pmemp2)   !20140620-4

            if(icrobd.eq.1)then
              if(pm1.lt.bl)pm1=bl
              if(pm1.gt.bu)pm1=bu
              if(pm2.lt.bl)pm2=bl
              if(pm2.gt.bu)pm2=bu
            else
              if((pm1.lt.bl).or.(pm1.gt.bu))goto 12
              if((pm2.lt.bl).or.(pm2.gt.bu))goto 12
            endif
            gtall(j1,idv)=pm1
            gtall(j2,idv)=pm2
          enddo
        else
          do idv=1,ndv
            gtall(j1,idv)=gtall(m1,idv)
            gtall(j2,idv)=gtall(m2,idv)
          enddo
        endif

      endif

      return
      end subroutine CrossSBX

! --------------------------------------------------
      subroutine CrossNEWSBX(ilnum,j1,j2,m1,m2,ifix)
! --------------------------------------------------
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      REAL(DP) :: csize, beta, expv, rand,alpha,betaq,bl,bu,eta,pm1,pm2
      REAL(DP) :: pmemp1,pmemp2,rate,u 
      INTEGER  :: idveq,idv
      INTEGER  :: ilnum,j1,j2,m1,m2,ifix
!HY20140404-6       USE messy_airtraf_tools_ga_function,     ONLY: ran1
      rate = rcros(ilnum)
      eta  = pcros(ilnum)
      csize = 0.000001_dp   !20140620-4

!  < NEW SIMULATED BINARY CROSSOVER >

11    continue
      idveq=0

!cMAY      if(ran1(idum).le.rate)then 
        do idv=1,ndv
          pmemp1 = gtall(m1,idv)
          pmemp2 = gtall(m2,idv)
          bl     = gelb(ilnum,idv)
          bu     = geub(ilnum,idv)

          u = ran1(idum)
!cMAY          if(u.le.0.5)then
          if(u.le.rate)then
            if(abs(pmemp1-pmemp2).gt.csize)then
              if(pmemp2.gt.pmemp1)then
                 pm1 = pmemp1
                 pm2 = pmemp2
              else
                 pm1 = pmemp2
                 pm2 = pmemp1
              endif
        
              if((pm1-bl).gt.(bu-pm2))then
                beta = 1.0_dp + (2.0_dp*(bu-pm2)/(pm2-pm1))   !20140620-4
              else
                beta = 1.0_dp + (2.0_dp*(pm1-bl)/(pm2-pm1))   !20140620-4
              endif

!     ------- FIND ALPHA -----
              expv = eta + 1.0_dp   !20140620-4
              beta = 1.0_dp / beta   !20140620-4
              alpha = 2.0_dp - beta**expv   !20140620-4
  
              rand=ran1(idum)
              if(rand.le.(1.0_dp/alpha))then   !20140620-4
                alpha = alpha * rand
                expv  = 1.0_dp/(eta+1.0_dp)   !20140620-4
                betaq = alpha**expv
              else
                alpha = alpha * rand
                alpha = 1.0_dp/(2.0_dp-alpha)   !20140620-4
                expv  = 1.0_dp/(eta+1.0_dp)   !20140620-4 
                if(alpha.lt.0.0_dp) stop 'crossover'   !20140620-4
                betaq = alpha**expv
              endif

              gtall(j1,idv)=0.5_dp*((1.0_dp+betaq)*pm1+(1.0_dp-betaq)*pm2)   !20140620-4
              gtall(j2,idv)=0.5_dp*((1.0_dp-betaq)*pm1+(1.0_dp+betaq)*pm2)   !20140620-4

            else 

              betaq = 1.0_dp   !20140620-4
              pm1 = pmemp1
              pm2 = pmemp2

              gtall(j1,idv)=0.5_dp*((1.0_dp+betaq)*pm1+(1.0_dp-betaq)*pm2)   !20140620-4
              gtall(j2,idv)=0.5_dp*((1.0_dp-betaq)*pm1+(1.0_dp+betaq)*pm2)   !20140620-4

            endif

            if(gtall(j1,idv).lt.bl)gtall(j1,idv)=bl
            if(gtall(j1,idv).gt.bu)gtall(j1,idv)=bu
            if(gtall(j2,idv).lt.bl)gtall(j2,idv)=bl
            if(gtall(j2,idv).gt.bu)gtall(j2,idv)=bu

          else
            gtall(j1,idv)=pmemp1
            gtall(j2,idv)=pmemp2
            idveq=idveq+1
          endif  
        enddo
        if((idveq.eq.ndv).and.(idvsm.ne.1))goto 11
!cMAY      endif

      return
      end subroutine CrossNEWSBX


!--------------------------------------------------------------------
!
!
! ORIGINAL: messy_airtraf_tools_ga_mutation.f90
!
!
! *****************************************************
      subroutine Mutation(ilnum,ioff)
! *****************************************************
!      mmute = Mutation Method
!         0 RANDOM MUTATION
!         1 UNIFORM MUTATION
!         2 POLYNOMIAL MUTATION
!         3 NEW POLYNOMIAL MUTATION
!        (4) UNDER CONSTRUCTION (NON UNIFORM MUTATION)
!        (5) UNDER CONSTRUCTION (NORMAL MUTATION)
! *****************************************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      INTEGER :: ioff,ilnum,method
      method = mmute(ilnum)

      if(method.eq.0)call RandomMutation(ilnum,ioff)
      if(method.eq.1)call UniformMutation(ilnum,ioff)
      if(method.eq.2)call PolynomialMutation(ilnum,ioff)
      if(method.eq.3)call NewPolynomialMutation(ilnum,ioff)
!      call NonUnifMutation(ilnum,ioff)

      return
      end subroutine Mutation

! --------------------------------------------
      subroutine UniformMutation(ilnum,ioff)
! --------------------------------------------
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6       USE messy_airtraf_tools_ga_function,     ONLY: ran1
      IMPLICIT NONE
      REAL(DP) :: bl,bu,parm,pm,r,rate,rmu
      INTEGER  :: ioff,ilnum,idv
      rate = rmute(ilnum)
      parm = pmute(ilnum)

      do idv=1,ndv
        if(ran1(idum).le.rate)then
!cAPR          pm = pmem(ioff,idv)
          pm = gtall(ioff,idv)
          bl = gelb(ilnum,idv)
          bu = geub(ilnum,idv)

          if(pm.gt.bl)then
            rmu = (ran1(idum)-0.5_dp)*parm   !20140620-5
!cAPR  10        pm = pmem(ioff,idv) + rmu
  10        pm = gtall(ioff,idv) + rmu
            if((pm.lt.bl).or.(pm.gt.bu))then
              rmu = rmu * 0.9_dp   !20140620-5
              goto 10 
            endif
          else
            r = ran1(idum)
            pm = r * (bu-bl) + bl 
          endif
          if(pm.lt.bl)pm=bl
          if(pm.gt.bu)pm=bu
!cAPR          pmem(ioff,idv) = pm
          gtall(ioff,idv) = pm
        endif
      enddo

      return
      end subroutine UniformMutation

! --------------------------------------------
      subroutine RandomMutation(ilnum,ioff)
! --------------------------------------------
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      REAL(DP) :: pmu,rate
      INTEGER  :: ilnum,ioff,idv
!HY20140404-6       USE messy_airtraf_tools_ga_function,     ONLY: ran1
      rate = rmute(ilnum)

      do idv=1,ndv
        if(ran1(idum).le.rate)then
  10      pmu = ( geub(ilnum,idv) - gelb(ilnum,idv) )&
            * ran1(idum) + gelb(ilnum,idv)
          if(pmu.lt.gelb(ilnum,idv))pmu=gelb(ilnum,idv)
          if(pmu.gt.geub(ilnum,idv))pmu=geub(ilnum,idv)
!cAPR          pmem(ioff,idv) = pmu
          gtall(ioff,idv) = pmu
        endif
      enddo

      return
      end subroutine RandomMutation

! --------------------------------------------------
      subroutine PolynomialMutation(ilnum,ioff)
! --------------------------------------------------
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      REAL(DP) :: poleta,bl,bu,delta,pm,r,rate,rmu 
      INTEGER  :: ilnum,ioff,idv
!HY20140404-6       USE messy_airtraf_tools_ga_function,     ONLY: ran1
      rate = rmute(ilnum)
      poleta = pmute(ilnum)

      do idv=1,ndv
        if(ran1(idum).le.rate)then
!cAPR          pm = pmem(ioff,idv)
          pm = gtall(ioff,idv)
          bl = gelb(ilnum,idv)
          bu = geub(ilnum,idv)

          if(pm.gt.bl)then
            r = ran1(idum)
            if(r.lt.0.5_dp)then   !20140620-4
              delta = (2.0_dp*r)**(1.0_dp/(poleta+1.0_dp))-1    !20140620-5
            else
              delta = 1.0_dp-(2.0_dp*(1.0_dp-r))**(1.0_dp/(poleta+1.0_dp))   !20140620-5
            endif
            rmu = (bu-bl)*delta
!cAPR  10        pm = pmem(ioff,idv) + rmu
  10        pm = gtall(ioff,idv) + rmu
            if((pm.lt.bl).or.(pm.gt.bu))then
              rmu = rmu * 0.9_dp   !20140620-5
              goto 10 
            endif
          else
            r = ran1(idum)
            pm = r * (bu-bl) + bl 
          endif
          if(pm.lt.bl)pm=bl
          if(pm.gt.bu)pm=bu
!cAPR          pmem(ioff,idv) = pm
          gtall(ioff,idv) = pm
        endif
      enddo

      return
      end subroutine PolynomialMutation

! --------------------------------------------------
      subroutine NewPolynomialMutation(ilnum,ioff)
! --------------------------------------------------
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      REAL(DP) :: etap, pol, delq,bl,bu,delta,eta,pm,r,rate
      INTEGER  :: ilnum,ioff,idv
!HY20140404-6       USE messy_airtraf_tools_ga_function,     ONLY: ran1
      rate = rmute(ilnum)
      eta  = pmute(ilnum)

      do idv=1,ndv
        if(ran1(idum).le.rate)then
!cAPR          pm = pmem(ioff,idv)
          pm = gtall(ioff,idv)
          bl = gelb(ilnum,idv)
          bu = geub(ilnum,idv)

          if(pm.gt.bl)then
            if((pm-bl).lt.(bu-pm))then
              delta = (pm-bl)/(bu-bl)
            else
              delta = (bu-pm)/(bu-bl)
            endif
 
            r = ran1(idum)
            etap = 1.0_dp/(eta+1.0_dp)   !20140620-5
 
            if(r.le.0.5_dp)then   !20140620-5
              pol  = 2.0_dp*r+(1.0_dp-2.0_dp*r)*(1.0_dp-delta)**(eta+1.0_dp)   !20140620-5
              delq = pol**etap - 1.0_dp   !20140620-5
            else
              pol  = 2.0_dp*(1.0_dp-r)+2.0_dp*(r-0.5_dp)*(1.0_dp-delta)**(eta+1.0_dp)   !20140620-5
              delq = 1.0_dp - pol**etap   !20140620-5
            endif

            pm = pm + delq * (bu-bl) 
          else
            r = ran1(idum)
            pm = r * (bu-bl) + bl 
          endif
          if(pm.lt.bl)pm=bl
          if(pm.gt.bu)pm=bu
!cAPR          pmem(ioff,idv) = pm
          gtall(ioff,idv) = pm
        endif
      enddo

      return
      end subroutine NewPolynomialMutation

!-----------------------------------------------------------------
!
!
! ORIGINAL: messy_airtraf_tools_ga_adapt.f90
!
!
! *****************************************
      subroutine DecodeGen
! *****************************************
!     decode from genotype to phenotype
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter 
      IMPLICIT NONE
      INTEGER :: ist,ied,ilnum
      do ilnum=1,nisland 
       ist = istat(ilnum,13)-istat(ilnum,11)+1
       ied = istat(ilnum,13)
       call DecodePop(ilnum,ist,ied,1,ndv)
      enddo

      return
      end subroutine DecodeGen

! *********************************************************
      subroutine DecodePop(ilnum,ippst,ipped,idvst,idved)
! *********************************************************
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6       USE messy_airtraf_tools_ga_iofile,      ONLY: ErrorPtype
      IMPLICIT NONE
      REAL(DP) :: flo,fup,gtype,ptype,pot,sdv
      INTEGER  :: idvst,idved,ippst,ipped,ilnum,idv,ipop
 
!     reset 
      do ipop=ippst,ipped
        do idv=idvst,idved
          ptall(ipop,idv) = 0.0_dp   !20140620-5
        enddo
      enddo

      do idv=idvst,idved
        pot = pout(ilnum,idv)
        sdv = stdv(ilnum,idv)
        flo = pave(ilnum,idv) - flat(ilnum,idv)
        fup = pave(ilnum,idv) + flat(ilnum,idv)

        do ipop=ippst,ipped
          gtype = gtall(ipop,idv)
          ptype = Decode(gtype,pot,sdv,flo,fup)
          ptall(ipop,idv) = ptype
          if((ptype.lt.phlb(idv)).or.(ptype.gt.phub(idv)))& 
           call ErrorPtype(ilnum,ipop,idv)
        enddo
      enddo

      return
      end subroutine DecodePop

! -----------------------------------------------
      function Decode(gtype,pout,stdv,phlo,phup)
! -----------------------------------------------
!HY20140404-6       USE messy_airtraf_tools_ga_function,    ONLY: xmin1
      IMPLICIT NONE
      REAL(DP) :: Decode
      REAL(DP) :: g1,phup,phlo,stdv,gtype,pout,fnor,g0
      if(gtype.lt.0.5_dp*pout)then   !20140620-5
        g0 = gtype/pout
        Decode = stdv * xmin1(g0,fnor) + phlo
      elseif(gtype.gt.(1.0_dp-0.5_dp*pout))then   !20140620-5
        g1 = 1.0_dp-0.5_dp*pout   !20140620-5
        g0 = 0.5_dp*((gtype-g1)/(1.0_dp-g1))+0.5_dp   !20140620-5
        Decode = stdv * xmin1(g0,fnor) + phup
      else
        g0 = (gtype-0.5_dp*pout)/(1.0_dp-pout)    !20140620-5
        Decode = (phup-phlo) * g0 + phlo
      endif

      return
      end function Decode

! -----------------------------------------------
      function Encode(ptype,pout,stdv,phlo,phup)
! -----------------------------------------------
!HY20140404-6       USE messy_airtraf_tools_ga_function,     ONLY: erfc1
      IMPLICIT NONE
      REAL(DP) :: Encode
      REAL(DP) :: p0,phlo,phup,stdv,pout,ptype
 
      if(ptype.lt.phlo)then
        p0 = (ptype-phlo)/stdv
        Encode = pout - erfc1(p0) * pout
      elseif(ptype.gt.phup)then
        p0 = (ptype-phup)/stdv
        Encode =   1.0_dp - erfc1(p0) * pout   !20140620-5
      else
        p0 = (ptype-phlo)/(phup-phlo)
        Encode = p0 * (1.0_dp-pout) + 0.5_dp * pout   !20140620-5
      endif

      return
      end function Encode

! *****************************************************
      subroutine GenotypeRange(ilnum)
! *****************************************************
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6       USE messy_airtraf_tools_ga_function,     ONLY: erfc1
      IMPLICIT NONE
!check
!C1.0.1      write(6,*)'Range Adaptation:',igen
!C1.2.0      write(6,*)'Range Adaptation:',igen+1
      REAL(DP) :: flup, fllo, pup, plo,g0
      INTEGER  :: ilnum,idv
      write(6,*)'Range Adaptation:',iter

      do idv=1,ndv
        flup = pave(ilnum,idv) + flat(ilnum,idv)
        fllo = pave(ilnum,idv) - flat(ilnum,idv)

! .. Upper Boundary
        if( phub(idv).lt.fllo )then
          g0 = (phub(idv)-fllo)/stdv(ilnum,idv)
          geub(ilnum,idv) = pout(ilnum,idv) - erfc1(g0)*pout(ilnum,idv)
        elseif( phub(idv).gt.flup )then
          g0 = (phub(idv)-flup)/stdv(ilnum,idv)
          geub(ilnum,idv) = 1.0_dp              - erfc1(g0)*pout(ilnum,idv)   !20140620-5
        else
          g0 = (phub(idv)-fllo)/(flup-fllo)
          geub(ilnum,idv) = g0*(1.0_dp-pout(ilnum,idv))+0.5_dp*pout(ilnum,idv)   !20140620-5
        endif
        geub(ilnum,idv) = min ( geub(ilnum,idv), 1.0_dp )   !20140620-5

! .. Lower Boundary
        if( phlb(idv).lt.fllo )then
          g0 = (phlb(idv)-fllo)/stdv(ilnum,idv)
          gelb(ilnum,idv) = pout(ilnum,idv) - erfc1(g0)*pout(ilnum,idv)
        elseif( phlb(idv).gt.flup )then
          g0 = (phlb(idv)-flup)/stdv(ilnum,idv)
          gelb(ilnum,idv) = 1.0_dp              - erfc1(g0)*pout(ilnum,idv)   !20140620-5
        else
          g0 = (phlb(idv)-fllo)/(flup-fllo)
          gelb(ilnum,idv) = g0*(1.0_dp-pout(ilnum,idv))+0.5_dp*pout(ilnum,idv)   !20140620-5
        endif
        gelb(ilnum,idv) = max ( gelb(ilnum,idv), 0.0_dp )   !20140620-5

        pup = Decode(geub(ilnum,idv),pout(ilnum,idv),stdv(ilnum,idv),&
               fllo,flup)
        plo = Decode(gelb(ilnum,idv),pout(ilnum,idv),stdv(ilnum,idv),&
               fllo,flup)
        if(infopt.eq.1)write(mtout,11)idv,gelb(ilnum,idv),&
           geub(ilnum,idv),'|',plo,fllo,flup,pup
!        if(infopt.eq.1)write(mtout,10)idv,gelb(ilnum,idv),
!     &      geub(ilnum,idv),'|',plo,pup
!        if(infopt.eq.1)write(mtout,10)idv,plo,fllo,'|',flup,pup
      enddo
   10 format(i4,2f12.6,a3,2f12.6)
   11 format(i4,2f8.3,a3,4f10.4)
      return
      end subroutine GenotypeRange

! **************************************************
      subroutine RandomGeneration(ilnum,info)
! ************************************************
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6       USE messy_airtraf_tools_ga_testfunction,   ONLY: GvalPop
!HY20140404-6       USE messy_airtraf_tools_ga_function,     ONLY: ran1 
      IMPLICIT NONE
      REAL(DP) :: gtype
      INTEGER  :: ilnum,info,icn,idv,ied,iged,igst,ipop,ist
      if(info.le.0)then
        ist = istat(ilnum,10)
      else
        ist = info
      endif
      ied = istat(ilnum,10) + istat(ilnum,11) - 1

      igst = imoga(ilnum,7)
      iged = imoga(ilnum,8)

      do ipop=ist,ied
      !write(*,*)"c1",ist,ied,ipop
  20    continue
        do idv=1,ndv
  10      gtype = ( geub(ilnum,idv) - gelb(ilnum,idv) )& 
           * ran1(idum) + gelb(ilnum,idv)
          if( (gtype.lt.gelb(ilnum,idv)) .or.& 
             (gtype.gt.geub(ilnum,idv)) )goto 10
          gtall(ipop,idv) = gtype
          !write(*,*)"c2",ist,ied,ipop

        enddo
        !write(*,*)"c3",ilnum,ipop,ipop,ndv,ist,ied
        call DecodePop(ilnum,ipop,ipop,1,ndv)
        !write(*,*)"c4",ilnum,ipop,ipop,ndv,ist,ied

!cDS        if(icnrgeo.eq.1)then
        !write(*,*)"c5",igst,ist,ied,ipop
        if(igst.gt.0)then
          !write(*,*)"c6",igst
          call GvalPop(ipop,ipop,igst,iged)
          do icn=igst,iged
            if(gvall(ipop,icn).gt.gvaltol)goto 20
          enddo
        endif 
!HY        write(*,*)"c7",igst,iged,ist,ied,istat(ilnum,13),ipop,idarc(ilnum,istat(ilnum,13))
        istat(ilnum,13)=istat(ilnum,13)+1
        idarc(ilnum,istat(ilnum,13))=ipop
!HY        write(*,*)"c8",igst,iged,ist,ied,istat(ilnum,13),ipop,idarc(ilnum,istat(ilnum,13))
      enddo
       
      return
      end subroutine RandomGeneration

! *******************************************************
      subroutine ArcAdapt(ilnum,idra,nselect)
! *******************************************************
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6       USE messy_airtraf_tools_ga,          ONLY: AssignFitness
!HY20140404-6       USE messy_airtraf_tools_ga_afit,     ONLY: PopSorting
!HY20140404-6       USE messy_airtraf_tools_ga_update,   ONLY: MaxminObj
!HY20140404-6       USE messy_airtraf_tools_ga_archive,  ONLY: SelectArchive
      IMPLICIT NONE
      INTEGER :: idra(mtpop),idtmp(maxpop),iorder(maxpop)
!HY20140405-1      dimension fdec(maxpop)
      REAL(DP) :: fdec(maxpop,mobj)       
      INTEGER :: iminpop, ir, k,ilnum,nselect,i,iarced,ip,irk,is,nfitpop,nmax,ntmp
 
      ntmp=0
      nselect=0
      nfitpop=istat(ilnum,12)
      if(arpopratio.gt.0.0_dp)then   !20140620-6
        iminpop = int(arpopratio*istat(ilnum,12))
        if(iminpop.lt.2)iminpop=2
        call AssignFitness(ilnum)
        call PopSorting(ilnum,nfitpop,iorder)
        do ip=1,iminpop
          idtmp(ip)=ifitpop(ilnum,ip)
        enddo
        ntmp=iminpop
      else
!        nfitpop=istat(ilnum,12)
        call MaxminObj(ilnum,nfitpop,fdec)
      endif
 
!CHECK
      if(idebug.eq.1)then
        do i=1,ntmp
          write(mtout,*)'ARPOP:',i,idtmp(i)
        enddo
      endif
 
      if(ararcratio.gt.0.0_dp)then   !20140620-6
        iminpop=int(ararcratio*istat(ilnum,11))
        if(iarbase.eq.1)iminpop=int(ararcratio*istat(ilnum,13))
        if(iminpop.lt.2)iminpop=2
        iarced  = istat(ilnum,13)
        irk = 1
        do while(nselect.lt.iminpop)
          call SelectArchive(ilnum,1,irk,iarced,idra,nselect)
!          write(mtout,*)'ARARC',igen,irk,nselect
          irk = irk + 1
        enddo
        nmax=nselect
        if(ntmp.gt.0)then
          do ip=1,ntmp
            ir=0
            do is=1,nmax
              if(idtmp(ip).eq.idra(is))ir=1
            enddo
            if(ir.eq.0)then
              nselect=nselect+1
              idra(nselect)=idtmp(ip)
            endif
          enddo
        endif
      else
        nselect=ntmp
        do ip=1,nselect
          idra(ip)=idtmp(ip)
        enddo
      endif

!      write(mtout,*)'AFTER ARPOP:',igen,nselect
    
!CHECK
      if(idebug.eq.1)then
        do i=1,nselect
          write(mtout,*)'ARALL:',i,idra(i)
        enddo
      endif
 
      if( (ncon.gt.0).and.(icnrvio.eq.1) )then
        irk=istat(ilnum,14)+1
        if(irk.gt.nselect)then
          iarced=istat(ilnum,13)
          call SelectArchive(ilnum,irk,irk,iarced,idra,nselect)
        endif
!        call SelectViolation(ilnum,idra,nselect)
!C1.0.1        write(mtout,*)'AFTER SV:',igen,nselect
!C1.2.0        write(mtout,*)'AFTER SV:',igen+1,nselect
        write(mtout,*)'AFTER SV:',iter,nselect
        write(mtout,*)'VIO',nselect-1,idra(nselect-1),&
          (gvall(idra(nselect-1),k),k=1,ncon)
        write(mtout,*)'VIO',nselect,idra(nselect),&
          (gvall(idra(nselect),k),k=1,ncon)
      endif

!      do ip=1,ilnum(istat,13)
!        ipop=idarc(ilnum,ip)
!        ipinfo(ipop,4)=-1
!      enddo
!      do ip=1,nselect
!        ipop=idra(ip)
!        ipinfo(ipop,4)=1
!      enddo
 
      return
      end subroutine ArcAdapt  

! *******************************************************
      subroutine UpdateStatistic(ilnum,idra,nselect)
! *******************************************************
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      INTEGER :: idra(mtpop)
      REAL(DP) :: pup, plo, pave0, stdv0
      INTEGER  :: ilnum,nselect,idv,ip,ipop
      do idv=1,ndv
        pup = pave(ilnum,idv) - flat(ilnum,idv)
        plo = pave(ilnum,idv) + flat(ilnum,idv)
        pave0 = pave(ilnum,idv)
        stdv0 = stdv(ilnum,idv)
        pave(ilnum,idv) = 0.0_dp   !20140620-6
        stdv(ilnum,idv) = 0.0_dp   !20140620-6
 
        do ip=1,nselect
          ipop = idra(ip)
          pave(ilnum,idv) = pave(ilnum,idv) + ptall(ipop,idv)
          pup = max(pup,ptall(ipop,idv))
          plo = min(plo,ptall(ipop,idv))
        enddo
        pave(ilnum,idv) = pave(ilnum,idv) / DBLE(nselect)   !HY20140414-6

        if ( (pave(ilnum,idv).lt.phlb(idv)) .or.&
            (pave(ilnum,idv).gt.phub(idv)) ) stop 'Statistic Error'

 
        do ip=1,nselect
          ipop = idra(ip)
          stdv(ilnum,idv) = stdv(ilnum,idv) + & 
          ( ptall(ipop,idv) - pave(ilnum,idv) )**2
        enddo
        stdv(ilnum,idv) = sqrt( stdv(ilnum,idv) / DBLE(nselect) )   !HY20140414-6
 
        pave(ilnum,idv) = pave0 + rave(ilnum) * (pave(ilnum,idv)-pave0)
!cMAY        stdv(ilnum,idv) = stdv0 + rstd(ilnum) * (stdv(ilnum,idv)-stdv0)

        if(rstd(ilnum).lt.0.0_dp)then   !20140620-6
          stdv(ilnum,idv) = stdv0 + & 
            abs(rstd(ilnum))*(stdv(ilnum,idv)-stdv0)
        else
         if(stdv(ilnum,idv).lt.stdv0) &
         stdv(ilnum,idv) = stdv0 + rstd(ilnum)*(stdv(ilnum,idv)-stdv0)
        endif
 
        if(pup.gt.phub(idv))pup = phub(idv)
        if(plo.lt.phlb(idv))plo = phlb(idv)
        pave(ilnum,idv) = 0.5_dp * (pup + plo)   !20140620-6
        flat(ilnum,idv) = 0.5_dp * (pup - plo)   !20140620-6
        pout(ilnum,idv) = rout(ilnum)
      enddo
 
      return
      end subroutine UpdateStatistic 

! *********************************************************
      subroutine EncodeCheck(ilnum,ippst,ipped,idvst,idved,ityp)
! *********************************************************
!HY20140401-6      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6       USE messy_airtraf_tools_ga_iofile,    ONLY: ErrorPtype
      IMPLICIT NONE
!     reset
      REAL(DP) :: flo, fup,gtype,pot,ptype,sdv 
      INTEGER  :: ityp,idvst,idved,ippst,ipped,ilnum,idv,ipop
      do ipop=ippst,ipped
        do idv=idvst,idved
          gtall(ipop,idv) = 0.0_dp   !20140620-6
        enddo
      enddo

      do ipop=ippst,ipped
        ipinfo(ipop,4)=-1 
        do idv=idvst,idved
          pot = pout(ilnum,idv)
          sdv = stdv(ilnum,idv)
          flo = pave(ilnum,idv) - flat(ilnum,idv)
          fup = pave(ilnum,idv) + flat(ilnum,idv)
 
          ptype = ptall(ipop,idv)
          if((ptype.lt.phlb(idv)).or.(ptype.gt.phub(idv)))then
            ipinfo(ipop,4)=3
            goto 12
          endif
          if(ityp.eq.1)then
            gtype=Encode(ptype,pot,sdv,flo,fup)
            gtall(ipop,idv) = gtype
!cMAY            if((gtype.lt.phlb(idv)).or.(gtype.gt.phub(idv)))
            if((gtype.lt.gelb(ilnum,idv)).or.(gtype.gt.geub(ilnum,idv))) &
             call ErrorPtype(ilnum,ipop,idv)
          endif
        enddo
        ipinfo(ipop,4)=1
 12     continue
      enddo
 
      return
      end subroutine EncodeCheck

!-------------------------------------------------------------
!
!
! ORIGINAL: messy_airtraf_tools_ga_iofile.f90
!
!
! *********************************
      subroutine InputGen(icheck)
! *********************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
!HY20140404-6       USE messy_airtraf_tools_ga_adapt,     ONLY: GenotypeRange, DecodePop
      IMPLICIT NONE
      INTEGER :: nd, ick,icheck,icn,idv,ii,ilnum,iof,ip,ipop,j,k,mt,numpop
      mt = mtgen

      if(mt.le.0)return

! ...  < Read Restart File >  ...
      do j=1,8
        read(mt,*)(iv(4*(j-1)+k),k=1,4)
      enddo
      do idv=1,ndv
        read(mt,*)phub(idv),phlb(idv)
        do ii=1,nisland
          read(mt,*)pave(ii,idv),stdv(ii,idv)
          read(mt,*)flat(ii,idv),pout(ii,idv)
        enddo
      enddo

      do ip=1,npop
        ipop=nallpop+1
        nallpop=ipop
        read(mt,*)nd,ilnum
        do iof=1,nobj
          read(mt,*)fvall(ipop,iof)
        enddo
        if(ncon.gt.0)then
          do icn=1,ncon
            read(mt,*)gvall(ipop,icn)
          enddo
        endif
        do idv=1,ndv
          read(mt,*)gtall(ipop,idv)
        enddo
!C1.2.0        ipinfo(ipop,1)=igen
        ipinfo(ipop,1)=iter
        ipinfo(ipop,2)=ip
        ipinfo(ipop,3)=ilnum
        ipinfo(ipop,4)=0
        istat(ilnum,11)=istat(ilnum,11)+1
        istat(ilnum,12)=istat(ilnum,12)+1
        istat(ilnum,13)=istat(ilnum,13)+1
        ifitpop(ilnum,istat(ilnum,12))=ipop
        idarc(ilnum,istat(ilnum,13))=ipop
      enddo

      if(ibestn.gt.0)then
        do ip=1,npop*ibestn
          read(mt,*,end=10)ipop,ilnum
          do iof=1,nobj
            read(mt,*)fvall(ipop,iof)
          enddo
          if(ncon.gt.0)then
            do icn=1,ncon
              read(mt,*)gvall(ipop,icn)
            enddo
          endif
          do idv=1,ndv
            read(mt,*)gtall(ipop,idv)
          enddo  
          if(iflag(9).eq.0)then
            ipinfo(ipop,1)=0
            ipinfo(ipop,2)=ip
          endif
          ipinfo(ipop,3)=ilnum
          ipinfo(ipop,4)=1
          istat(ilnum,12)=istat(ilnum,12)+1
          ifitpop(ilnum,istat(ilnum,12))=ipop
        enddo
 10     continue      
      endif

      istat(1,10)=nallpop+1
      if(nisland.gt.1)then
        do ilnum=1,nisland
          istat(ilnum,10)=istat(ilnum-1,10)+istat(ilnum-1,11)
        enddo
      endif

! ...  < Genotype Range > ...
      ick=0
      do ilnum=1,nisland
        if(pout(ilnum,1).gt.0.0_dp)then   !20140620-6
          call GenotypeRange(ilnum)
          ick=ick+1
        endif
      enddo
      if(ick.gt.0)icheck=icheck+1
      
! ...  < Decode from Genotype to Phenotype > ...
      do ilnum=1,nisland
        numpop=istat(ilnum,12)
        do ip=1,numpop
          ipop=ifitpop(ilnum,ip)
          call DecodePop(ilnum,ipop,ipop,1,ndv)
        enddo
      enddo

      return
      end subroutine InputGen

! *********************************
      subroutine InputByHand
! *********************************
      USE messy_airtraf_tools_ga_input   !HY20140402-2 
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
!HY20140401-7      include 'ga.par'
! >>> iast    ... Start generation of range adaptation
! >>>           if iast > ngen then no range adaptation is performed
! >>> iaint   ... Interval of range adaptation
! >>> rave    ... Relaxation coefficient for average (range adaptation)
! >>>  (0 < rave < 1) 1.0 is better
! >>> rstd    ... Relaxation coefficient for standard deviation (range adaptation)
! >>>  (0 < rstd < 1) 0.5 is better
! >>> rout    ... Outer ratio of variable distribution (range adaptation) 
! >>>  (0 < rout < 1) 
! >>> mcros   ... Crossover method
! >>> rcros   ... Crossover rate (0<1)
! >>> pcros   ... Parameter for Crossover
! >>> mmute   ... Mutation method
! >>> rmute   ... Mutation rate (0<1)
! >>> pmute   ... Parameter for Mutation
! >>> nelt    ... Number of keep individual (Elite)
! >>>             (It should be set to 0 for present MOGA)

! ....... Range Adaptation Start .......
!HY20140328-6      call ReadInteger(mtparam,iast)
! ....... Interval of Range Adaptation .......
!HY20140328-6      call ReadInteger(mtparam,iaint)
! ....... Relaxation Coef. of Average for Range Adaptation .......
!HY20140328-6      call ReadReal(mtparam,rave)
! ....... Relaxation Coef. of St.Deviation for Range Adaptation ...
!HY20140328-6      call ReadReal(mtparam,rstd)
! ....... Outer Rate of Population by ARMOGA ......
!HY20140328-6      call ReadReal(mtparam,rout)
! ....... Crossover Method .......
!HY20140328-6      call ReadInteger(mtparam,mcros)
! ....... Crossover Rate .......
!HY20140328-6      call ReadReal(mtparam,rcros)
! ....... Crossover Parameter .......
!HY20140328-6      call ReadReal(mtparam,pcros)
! ....... Mutation Method .......
!HY20140328-6      call ReadInteger(mtparam,mmute)
! ....... Mutation Rate .......
!HY20140328-6      call ReadReal(mtparam,rmute)
! ....... Mutation Parameter .......
!HY20140328-6      call ReadReal(mtparam,pmute)
! ....... Elitist .......
!HY20140328-6      call ReadInteger(mtparam,nelt)
! ....... Number of Parents for selection .......
!check call ReadInteger(mtparam,noffs)
! ....... Number of Parents for selection .......
!check call ReadInteger(mtparam,nselp)
      iast(nisland)  = iast_fort_1          !HY20140329-1,20140410-2
      iaint(nisland) = iaint_fort_1         !These parameteds are read by fort.1 originally.
      rave(nisland)  = rave_fort_1          !If misland.ne.1 (or nisland.ne.1), correction is needed.
      rstd(nisland)  = rstd_fort_1          !If AirTraf uses island model, see 20140328-6
      rout(nisland)  = rout_fort_1          
      mcros(nisland) = mcros_fort_1 
      rcros(nisland) = rcros_fort_1
      pcros(nisland) = pcros_fort_1
      mmute(nisland) = mmute_fort_1
      rmute(nisland) = rmute_fort_1
      pmute(nisland) = pmute_fort_1
      nelt(nisland)  = idummy_fort_1
!HY20140329-1 write(*,*)iast,iaint,rave,rstd,rout,mcros,rcros,pcros,mmute,rmute,pmute,nelt_
      return     
      end subroutine InputByHand

! ***************************************
      subroutine ReadInteger(iread,itmp)
! ***************************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      INTEGER :: itmp(misland)
      INTEGER :: iread,id,ii
      if(nisland.eq.1)then
       read(iread,*)itmp(1)
      else
       read(iread,*)id
       if(id.ne.0)then
        read(iread,*)itmp(1)
        do ii=2,nisland
         itmp(ii)=itmp(1)
        enddo
       else
        read(iread,*)(itmp(ii),ii=1,nisland)
       endif
      endif

      return
      end subroutine ReadInteger

! ***************************************
      subroutine ReadReal(iread,rtmp)
! ***************************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      REAL(DP) :: rtmp(misland)
      INTEGER :: iread,id,ii
      if(nisland.eq.1)then
       read(iread,*)rtmp(1)
      else
       read(iread,*)id
       if(id.ne.0)then
        read(iread,*)rtmp(1)
        do ii=2,nisland
         rtmp(ii)=rtmp(1)
        enddo
       else
        read(iread,*)(rtmp(ii),ii=1,nisland)
       endif
      endif

      return
      end subroutine ReadReal

! ***************************************************
      subroutine OutputSystemPtype(ist,ied,mt,itype)
! ***************************************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter 
      IMPLICIT NONE
      INTEGER :: ist,ied,mt,itype,icn,idv,ilnum,iof,ip
      INTEGER :: ipst,iped,ipop 
      if(mt.le.0)return
 
! ...  < Write Restart File >  ...
      open(mt)
      rewind(mt)
      do ipop=ist,ied
!C1.2.0        write(mt,*)igen, npop, nisland
        write(mt,*)iter, npop, nisland
        write(mt,*)nobj, ncon, ndv
        do iof=1,nobj
          write(mt,*)fvall(ipop,iof)
        enddo
        if(ncon.gt.0)then
          do icn=1,ncon
            write(mt,*)gvall(ipop,icn)
          enddo
        endif
        do idv=1,ndv
          if(itype.eq.0)write(mt,*)gtall(ipop,idv)
          if(itype.eq.1)write(mt,*)ptall(ipop,idv)
        enddo
      enddo
 
      do ilnum=1,nisland
        if(istat(ilnum,12).gt.istat(ilnum,11))then
          ipst=istat(ilnum,11)+1
          iped=istat(ilnum,12)
          do ip=ipst,iped
            ipop=ifitpop(ilnum,ip)
            write(mt,*)ipinfo(ipop,1),ipinfo(ipop,2), &
              ipinfo(ipop,3) 
            write(mt,*)nobj, ncon, ndv
            do iof=1,nobj
              write(mt,*)fvall(ipop,iof)
            enddo
            if(ncon.gt.0)then
              do icn=1,ncon
                write(mt,*)gvall(ipop,icn)
              enddo
            endif
            do idv=1,ndv
              if(itype.eq.0)write(mt,*)gtall(ipop,idv)
              if(itype.eq.1)write(mt,*)ptall(ipop,idv)
            enddo
          enddo
        endif
      enddo
 
      return
      end subroutine OutputSystemPtype 

! ***************************************************
      subroutine OutputGPtype(ist,ied,mt,itype)
! ***************************************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      INTEGER :: ist,ied,mt,itype,icn,idv,ilnum,iof,ip   
      INTEGER :: ipst,iped,ipop,j,k
      if(mt.le.0)return

! ...  < Write Restart File >  ...
      open(mt)
      rewind(mt)
!C1.2.0      write(mt,*)igen, npop, nisland
      write(mt,*)iter, npop, nisland
      write(mt,*)nobj, ncon, ndv
      write(mt,*)idum, iy
      do j=1,8
        write(mt,*)(iv(4*(j-1)+k),k=1,4)
      enddo
      do idv=1,ndv
        write(mt,*)phub(idv),phlb(idv)
        do ilnum=1,nisland
          write(mt,*)pave(ilnum,idv),stdv(ilnum,idv)
          write(mt,*)flat(ilnum,idv),pout(ilnum,idv)
        enddo
      enddo

      do ipop=ist,ied
        write(mt,*)ipop,ipinfo(ipop,3)
        do iof=1,nobj
          write(mt,*)fvall(ipop,iof)
        enddo
       if(ncon.gt.0)then
         do icn=1,ncon
           write(mt,*)gvall(ipop,icn)
         enddo
       endif
       do idv=1,ndv
        if(itype.eq.0)write(mt,*)gtall(ipop,idv)
        if(itype.eq.1)write(mt,*)ptall(ipop,idv)
       enddo
      enddo

      do ilnum=1,nisland
        if(istat(ilnum,12).gt.istat(ilnum,11))then
          ipst=istat(ilnum,11)+1
          iped=istat(ilnum,12)
          do ip=ipst,iped
            ipop=ifitpop(ilnum,ip)
            write(mt,*)ipop,ilnum
            do iof=1,nobj
              write(mt,*)fvall(ipop,iof)
            enddo
            if(ncon.gt.0)then
              do icn=1,ncon
                write(mt,*)gvall(ipop,icn)
              enddo
            endif
            do idv=1,ndv
              if(itype.eq.0)write(mt,*)gtall(ipop,idv)
              if(itype.eq.1)write(mt,*)ptall(ipop,idv)
            enddo
          enddo
        endif
      enddo

      return
      end subroutine OutputGPtype

! ***************************************
      subroutine OutputOperator
! ***************************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      REAL(DP) :: parm
      INTEGER  :: narld,ii,method
!  ------- explanation of EA ----------

      narld=0

! <<<  Optimization Problem  >>>
      write(mtgop,*)'##################################################'
      write(mtgop,10)'Optimization Name: ',testp
      if(imaxmin.gt.0)then
        write(mtgop,*)nobj,'-Objective Maximization Problem'
      else
        write(mtgop,*)nobj,'-Objective Minimization Problem'
      endif
      if(ncon.gt.0)then
        write(mtgop,*)' with',ncon,' Constraint(s)'
      endif
      write(mtgop,*)'Number of Design Variables',ndv

! <<<  Optimization Mode  >>>
      if(idebug.eq.1)then
        write(mtgop,*)'******** DEBUG MODE OPTIMIZATION ********'
        if(mtdeb.gt.0)write(mtgop,*)'Output File No.',mtdeb
      endif

! <<<  Optimization Method  >>>
      if(iman.eq.1)then
        write(mtgop,*)'Manual Optimization'
!C1.2.0        write(mtgop,*)'Present Generation',igen
        write(mtgop,*)'Present Generation',iter
      else
        write(mtgop,*)'Automatic Optimization'
      endif

! <<<  Genetic Operator Input  >>>
      if(idman.eq.1)then
        write(mtgop,*)'Manual Genetic Operator Input'
        write(mtgop,*)'Input File No.',mtparam
      else
        write(mtgop,*)'Semi-Auto Genetic Operator Input'
      endif

! <<<  Evaluation Process  >>>
      if(indveva.eq.1)then
        write(mtgop,*)'Each Evaluation is Performed Separately.'
       else
        write(mtgop,*)'All Evaluations are Performed Simultaneously.'
       endif

! <<<  Constraint Process  >>>
      if(ncon.gt.0)then
        if(indvcnr.eq.1)then
          write(mtgop,*)'Each Constraint is Calculated Separately.'
         else
          write(mtgop,*)'All Constraints are Calculated Simultaneously.'
         endif
       endif

! <<<  Data Preservation  >>>
      if(iflag(9).eq.1)then
        write(mtgop,*)'All Population Data is Stored.'
      else
        write(mtgop,*)'Only Latest Population Data is Stored.'
      endif
      if(imtall.eq.1)then
        write(mtgop,*)'All Data is Preserved.'
      else
        write(mtgop,*)'Only Final Data is Preserved.'
      endif
      if(irank1.eq.1)then
        write(mtgop,*)'Only Rank1 Population is Ranked. [*]'
      else
        write(mtgop,*)'All Population is Ranked Correctly.'
      endif
      if(mtall.gt.0)write(mtgop,*)'Output File No.',mtall

! <<<  I/O Setting  >>>
      if(mtout.gt.0)write(mtgop,*)'Normal Output File No.',mtout
      if(mtgen.gt.0)write(mtgop,*)'Latest Population File (P) No.',mtgen
      if(mteva.gt.0)write(mtgop,*)'Latest Population File (G) No.',mteva
      if(ipout.eq.1)then
        write(mtgop,*)' Separate Phenotype Data' 
      else
        write(mtgop,*)' Single Phenotype Data' 
      endif
      if(mtobj.gt.0)write(mtgop,*)'Latest Objective Functions No.',mtobj
      if(mtprt.gt.0)write(mtgop,*)'Latest Pareto Optimal No.',mtprt

! <<<  Operator Data  >>>
      if(ncon.gt.0)then
        if(icnrgeo.eq.0)then
          write(mtgop,*)'No Geometry Constraints Treatment'
        else
          write(mtgop,*)'Prevention of Infeasible Region'
        endif
        if(icnrrnk.eq.2)then
          write(mtgop,*)'Constrained Domination (single constraint)'
        else
          write(mtgop,*)'Constrained Domination'
        endif
      endif
         
      write(mtgop,*)'Number of Total Population',npop
      write(mtgop,*)'Number of Generations',ngend

      if(nisland.gt.1)then
        write(mtgop,*)'Distributed Model' 
!HY20140404-1        if(imocir.eq.1)then
!HY          write(mtgop,*)'Stepping Stone Model'
!HY        else
!HY          write(mtgop,*)'Random Migration Model'
!HY        endif
        write(mtgop,*)'Number of Islands',nisland
!HY20140404-1        write(mtgop,*)'Interval of Migration',iminter
!HY        write(mtgop,*)'Migration Rate',DBLE(imrate*nisland)/DBLE(npop),&  !HY20140414-6
!HY                       '(',imrate,')'
      endif
        
      do ii=1,nisland

        write(mtgop,'(A)')
        if(nisland.gt.1)write(mtgop,*)'<   Island Number',ii,'   >'
        if(ibestn.eq.1)then
          write(mtgop,*)'Elitist MOEA'
          if(ictrlelt.eq.1)then
            write(mtgop,*)'Controled Elitist Model'
            write(mtgop,*)'Control Coefficient:',rkcoef
          endif
        else
          write(mtgop,*)'Non-Elitist MOEA'
        endif
        if(iast(ii).lt.ngend)then
          write(mtgop,*)'Adaptive Range MOEA'
          write(mtgop,*)'First Adaptation:',iast(ii)
          write(mtgop,*)'Adaptation Interval:',iaint(ii)
          write(mtgop,*)'Population Ratio for Sampling:',aratio,&
            '(',int(aratio*npop),')' 
          write(mtgop,*)'Relaxation Coefficient for Standard Deviation:'&
            ,rstd(ii)
        endif

        if(alpdmpar.gt.0.0_dp)then   !20140620-6
          write(mtgop,*)'ALPHA-Domination'
          if(ialpnorm.eq.0)then
            write(mtgop,*)'No Normalization of Alpha'
          elseif(ialpnorm.eq.1)then
            write(mtgop,*)'Normalization of Alpha by All Max-Min'
          else
            write(mtgop,*)'Normalization of Alpha by Present Max-Min'
          endif
          write(mtgop,*)'D-ALPHA(i,j):',alpdmpar
        endif

        if(mprank.eq.1)then 
          write(mtgop,*)'Pareto Ranking Method (Fonseca)'
        else
          write(mtgop,*)'Front Ranking Method (Goldberg)'
        endif

        if(mfiteval.eq.2)then
          write(mtgop,*)'Average Fitness' 
        elseif(mfiteval.eq.1)then
          write(mtgop,*)'Rank-Based Fitness (F=1/rank)' 
        endif
        if(irksp.eq.1)write(mtgop,*)'Rank Order is Surely Preserved.'
        if(ishare.gt.0)then
          if(ibestn.eq.0)then
            write(mtgop,*)'Fitness Sharing for Selection'
          else
            write(mtgop,*)'Fitness Sharing for 2N Population Set'
            if(ishare.eq.2)write(mtgop,*)' and for Selection'
          endif
          if(iallsh.eq.1)then
            write(mtgop,*)'Sharing for All Population'
          else
            write(mtgop,*)'Sharing for Only Same Rank Population'
          endif
          write(mtgop,*)'Sharing-ALPHA:',shalpha
        endif

        if(mselection.eq.0)then
          write(mtgop,*)'Best-N Selection (No Selection)'
        elseif(mselection.eq.1)then
          write(mtgop,*)'Stochastic Universal Sampling (SUS)'
        elseif(mselection.eq.2)then
          write(mtgop,*)'Roulette Wheel Selection (RWS)'
        endif

        if(iarchiv.gt.0)then
          if(iarchiv.eq.1)then
            write(mtgop,*) &
              'Archiving is Included in Next Parents Set. (2N)'
          elseif(iarchiv.ge.2)then
            write(mtgop,*)&
              'Archiving is Included in Actual Parents Set. (N)'
            if(iarchiv.eq.3)&
              write(mtgop,*)'Only Next Generation of Adaptation' 
          endif
          write(mtgop,*)'Population Ratio for Archiving',aratio,&
            '(',int(arcratio*npop),')'
        endif

        if(rcros(ii).gt.0.0_dp)then    !20140620-6
          method=mcros(ii)
          parm=pcros(ii)
          write(mtgop,*)'Crossover Rate:',rcros(ii)
          if(method.eq.0)then 
            write(mtgop,*)'Arithmetic Crossover'
          elseif(method.eq.1)then 
            write(mtgop,*)'Blended Crossover (BLX)'
            write(mtgop,*)'Range ALPHA:',parm
          elseif(method.eq.2)then
            write(mtgop,*)'Simulated Binary Crossover (SBX): DV FIXED'
            write(mtgop,*)'ETAc:',parm
          elseif(method.eq.3)then
            write(mtgop,*)'Simulated Binary Crossover (SBX)'
            write(mtgop,*)'ETAc:',parm
          endif
        else
            write(mtgop,*)'NO CROSSOVER'
        endif
           
        if(rmute(ii).gt.0.0_dp)then   !20140620-6
          method=mmute(ii)
          parm=pmute(ii)
          write(mtgop,*)'Mutation Rate:',rmute(ii)
          if(method.eq.0)then
            write(mtgop,*)'Random Mutation'
          elseif(method.eq.1)then
            write(mtgop,*)'Uniform Mutation'
            write(mtgop,*)'Width of Disturbance:',parm
          elseif(method.eq.2)then
            write(mtgop,*)'Polynomial Mutation'
            write(mtgop,*)'ETAm:',parm
          endif
        endif

      enddo

      if(nelite.gt.0)then
       write(mtgop,*)'Elitism is used'
!HY20140404-1       if(ieltbst.eq.1)then
!HY        write(mtgop,*)'Only the best is preserved'
!HY        write(mtgop,*)'Number of Preserved Population',nelite
!HY       endif
      endif

!       if(nelite.gt.0.and.ieltbst.ne.1)
!     &  write(mtgop,*)'Number of Preserved Population',nelt(ii)
!       write(mtgop,*)'Selection Constant for Ranking',clr
!check  write(mtgop,*)'Maximum Selection Number as Parent(s)',nselp(ii)
!check  if(nparent.eq.1.or.cross(ii).lt.0.001)then
!check   write(mtgop,*)1,'  parent to',1,'  child'
!check   write(mtgop,*)'Only Mutation (No Crossover)'
!check  else
!check   write(mtgop,*)nparent,'  parent(s) to',noffs(ii),'  child(ren)'
!check   write(mtgop,*)'Crossover Rate',cross(ii)
!check  endif
!check   if(isptr(ii).gt.0)then
!check    write(mtgop,*)'Constrained Mutation'
!check    if(isptr(ii).eq.1)then
!check     write(mtgop,*)'G=0 Mutation is performed'
!check    else
!check     write(mtgop,*)'G<0 Mutation is performed'
!check    endif
!check   else
!check   endif
!      enddo

!check if(nisland.gt.1.and.narld.gt.0)then
!check  write(mtgop,'(A)')
!check  write(mtgop,*)'In case of large difference between range and DV'
!check  write(mtgop,*)'at immigration to adaptve island,'
!check  if(irres.eq.0)write(mtgop,*)
!check&  'adaptation is performed for all DVs'
!check  if(irres.eq.1)write(mtgop,*)
!check&  'adaptation is done only for its DV independently'
!check  if(irres.eq.2)write(mtgop,*)
!check&  'adaptation is done for its DV with same population statistic'
!check endif
  10  format(a20,a4)

      return
      end subroutine OutputOperator

! ****************************************************
      subroutine ErrorPtype(ilnum,ipop,idv)
! ****************************************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      REAL(DP) :: pup, plo, gup, glo,gtype,ptype
      INTEGER  :: ilnum,ipop,idv
      ptype = ptall(ipop,idv)
      gtype = gtall(ipop,idv)
      pup   = phub(idv)
      plo   = phlb(idv)
      gup   = geub(ilnum,idv)
      glo   = gelb(ilnum,idv)

      write(mtout,*)'******* ErrorPtype *********'
      write(mtout,*)'POP:',ipop,'DV:',idv
      write(mtout,*)'GTYPE:',gtype,' PTYPE:',ptype
      write(mtout,*)'(G) Upper Boundary:',gup
      write(mtout,*)'Genotype Value:',gtype
      write(mtout,*)'(G) Lower Boundary:',glo
      write(mtout,*)'(P) Upper Boundary:',pup
      write(mtout,*)'Phenotype Value:',ptype
      write(mtout,*)'(P) Lower Boundary:',plo

      if(ptype.gt.pup)then
        ptall(ipop,idv) = pup
        gtall(ipop,idv) = gup
      else
        ptall(ipop,idv) = plo
        gtall(ipop,idv) = glo
      endif

      if(istop.ne.0) stop

      return
      end subroutine ErrorPtype

! ****************************************************
      subroutine OutputCObj(ilnum,ist,ied,icst,iced,iall,mt)
! ****************************************************
! ##### iall 0 ...  Feasible Pareto Solution
! ##### iall 1 ...  All Feasible Solution
! ##### iall 2 ...  All Solution
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      REAL(DP) :: gs
      INTEGER  :: icst,iced,mt,iall,ist,ied,ilnum,icn,iof,ip,ipop
      if(mt.le.0)return

      rewind(mt)

      do ip=ist,ied
        gs = 0.0_dp   !20140620-6
        if(iced.gt.0)then
          ipop=idarc(ilnum,ip)
          do icn=icst,iced
            gs = gs + max(0.0_dp,gvall(ipop,icn))   !20140620-6
          enddo
          if( ((irank(ilnum,ip).eq.1).and.(gs.le.0.0_dp)) .or.&   !20140620-6
            ((gs.le.0.0_dp).and.(iall.eq.1)) .or. (iall.eq.2) )&   !20140620-6
             write(mt,*)(fvall(ipop,iof),iof=1,nobj)
        endif
      enddo

      return
      end subroutine OutputCObj

! *************************************************
      subroutine OutputObj(ilnum,ist,ied,iall,mt)
! *************************************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      INTEGER :: mt,iall,ilnum,ist,ied,iof,ip,ipop
      if(mt.le.0)return

      rewind(mt)
      do ip=ist,ied
        if((irank(ilnum,ip).eq.1).or.(iall.eq.1))then
          ipop=idarc(ilnum,ip)
          write(mt,*)(fvall(ipop,iof),iof=1,nobj)
        endif
      enddo

      return
      end subroutine OutputObj

! **********************************************
      subroutine OutputAllData(ist,ied,mt)
! **********************************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      INTEGER :: mt,ist,ied,icn,idv,igst,iged,iof,ipop,num
      if(mt.le.0)return

      if(iflag(9).eq.0)then
!C1.2.0        igst = igen-1
!C1.2.0        iged = igen-1
        igst = iter-2
        iged = iter-2
      else
        igst=ipinfo(ist,1)
        iged=ipinfo(ied,1)
      endif

      num=ied-ist+1
      rewind(mt)
      do ipop=ist,ied
        
        write(mt,*)ipinfo(ipop,1),ipinfo(ipop,2),ipinfo(ipop,3)
        write(mt,*)nobj,ncon,ndv
        do iof=1,nobj
          write(mt,*)fvall(ipop,iof)
        enddo
        if(ncon.gt.0)then
          do icn=1,ncon
            write(mt,*)gvall(ipop,icn)
          enddo
        endif
        do idv=1,ndv
          write(mt,*)ptall(ipop,idv)
        enddo
      enddo
    
      return
      end subroutine OutputAllData


! **************************************************
      subroutine InputAllData(mt,icheck)
! **************************************************
!HY20140401-7      include 'ga.par'
      USE messy_airtraf_tools_ga_parameter
      IMPLICIT NONE
      INTEGER :: ig, in,mt,icheck,iobj,icon,icn,idv,igp,ilnum,inum
      INTEGER :: iof,ip,ipop
      if(mt.le.0)return

      rewind(mt)
      inum=0
      igp=-1
      do ipop=1,mtpop
        read(mt,*,end=10)ig,ip,ilnum
        read(mt,*)iobj,icon,idv
        if(icheck.eq.1)then
          if((idv.ne.ndv).or.(iobj.ne.nobj).or.(icon.ne.ncon))&
            stop 'Failure of Read All Data (mtall)'
        endif 
        do iof=1,nobj
          read(mt,*)fvall(ipop,iof)
          if(fvall(ipop,iof).gt.objallmax(ilnum,iof))&
            objallmax(ilnum,iof)=fvall(ipop,iof)
          if(fvall(ipop,iof).lt.objallmin(ilnum,iof))&
            objallmin(ilnum,iof)=fvall(ipop,iof)
        enddo
        if(ncon.gt.0)then
          do icn=1,ncon
            read(mt,*)gvall(ipop,icn)
            if(gvall(ipop,iof).gt.gvalallmax(ilnum,iof))&
              gvalallmax(ilnum,iof)=gvall(ipop,iof)
            if(gvall(ipop,iof).lt.gvalallmax(ilnum,iof))&
              gvalallmax(ilnum,iof)=gvall(ipop,iof)
          enddo
        endif
        do idv=1,ndv
          read(mt,*)ptall(ipop,idv)
          gtall(ipop,idv)=0.0_dp   !20140620-6
        enddo
        if(ig.ne.igp)then
          igp=ig 
          in=0
        endif
        in=in+1
        ipinfo(ipop,1)=ig
        ipinfo(ipop,2)=in
        ipinfo(ipop,3)=ilnum
        ipinfo(ipop,4)=-1
        istat(ilnum,13)=istat(ilnum,13)+1
        idarc(ilnum,istat(ilnum,13))=ipop
        inum=inum+1
      enddo
 10   continue      

      nallpop=inum
      icheck=2

      return
      end subroutine InputAllData

!--------------------------------------------------------------
!
!
! ORIGINAL: messy_airtraf_tools_ga_function.f90
!
!
! ------------------------------------
      function ran1(idum)
! ------------------------------------
      USE messy_airtraf_tools_ga_parameter,  ONLY:iy,iv  !HY20140401-5,6
      IMPLICIT NONE
      REAL(DP)           :: ran1
      INTEGER, PARAMETER :: im=2147483647
      REAL(DP), PARAMETER :: am=1.0_dp/im, eps=1.2e-7, rnmx=1.0_dp-eps 
      INTEGER, PARAMETER :: ia=16807, iq=127773, ir=2836, &
                            ntab=32, ndiv=1+(im-1)/ntab
      INTEGER :: i0,idum,j,k 
!HY20140401-5      common/ranc/iy, iv(32)

      if(idum.eq.0)then
      do i0=1,ntab
      iv(i0)=0
      enddo
      iy = 0
      endif

      if (idum.le.0 .or. iy.eq.0 )then
         idum=max(-idum,1)
         do j=ntab+8,1,-1
            k=idum/iq
            idum=ia*(idum-k*iq)-ir*k
            if (idum.lt.0) idum=idum+im
            if (j.le.ntab) iv(j)=idum
         enddo
         iy=iv(1)
      endif
      k=idum/iq
      idum=ia*(idum-k*iq) -ir*k
      if (idum.lt.0) idum=idum+im
      j=1+iy/ndiv
      iy=iv(j)
      iv(j)=idum
      ran1=min(am*iy,rnmx)
      return
      end function ran1

! -----------------------------------------
      function func(x,s)
! -----------------------------------------
      IMPLICIT NONE
      REAL(DP) :: func 
      REAL(DP) :: s,x
      func = abs( erfc1(x) - (1.0_dp-s) )   !20140620-6
      return
      end function func


! -----------------------------------------
      function xmin1(s,fnor)
! -----------------------------------------
!     routine for golden section search
      IMPLICIT NONE
      REAL(DP) :: xmin1 
      REAL(DP), PARAMETER :: r=0.61803399_dp, c=1.0_dp-r   !20140620-6
      REAL(DP), PARAMETER :: tol=1.e-6
      REAL(DP) :: ax, bx, cx, x0, x3,s,fnor,f1,f2,x1,x2
      if( (s.gt.0.0_dp).and.(s.lt.1.0_dp) )then   !20140620-6
!     ok
      elseif(s.lt.0.5_dp)then   !20140620-6
      xmin1 = -3.0_dp   !20140620-6
      elseif(s.gt.0.5_dp)then   !20140620-6
      xmin1 =  3.0_dp   !20140620-6
      else
      write(6,*)'bad s value in xmin1 ',s
      stop
      endif
      
!     99% included
      ax = -2.6_dp   !20140620-6
      cx =  2.6_dp   !20140620-6

      bx = 0.5_dp*(ax+cx)   !20140620-6
      x0 = ax
      x3 = cx

      if(abs(cx-bx).gt.abs(bx-ax))then
       x1 = bx
       x2 = bx+c*(cx-bx)
      else
       x2 = bx
       x1 = bx-c*(bx-ax)
      endif
      f1 = func(x1,s)
      f2 = func(x2,s)
 1    if(abs(x3-x0).gt.tol*(abs(x1)+abs(x2)))then
       if(f2.lt.f1)then
        x0 = x1
        x1 = x2
        x2 = r*x1+c*x3
        f1 = f2
        f2 = func(x2,s)
        else
        x3 = x2
        x2 = x1
        x1 = r*x2+c*x0
        f2 = f1
        f1 = func(x1,s)
       endif
       goto 1
      endif
      if(f1.lt.f2)then
       fnor  = f1
       xmin1 = x1
       else
       fnor  = f2
       xmin1 = x2
      endif

      return
      end function xmin1
!
      function erf0(x)
!*****************************************
!*                2.      x   -t*t       *
!*  erf0(x) = ---------- S   e      dt   *
!*             sqrt(pi)   0              *
!*****************************************
      IMPLICIT NONE
      REAL(DP) :: erf0 
      REAL(DP) :: x 
      if(x.lt.0.0_dp)then
       erf0=-gammp(0.5_dp,x**2)
       else
       erf0= gammp(0.5_dp,x**2)
      endif
      return
      end function erf0
!
      function erf1(x)
!********************************************
!*                 1.      x   -t*t/2       *
!*   erf1(x) = ----------- S   e       dt   *
!*            sqrt(2.*pi)  0                *
!********************************************
      IMPLICIT NONE
      REAL(DP) :: erf1 
      REAL(DP) :: x
      if(x.lt.0.0_dp)then
       erf1=0.5_dp*(-gammp(0.5_dp,0.5_dp*x**2))   !20140620-7
       else
       erf1=0.5_dp*( gammp(0.5_dp,0.5_dp*x**2))   !20140620-7
      endif
      return
      end function erf1
!
      function erfc0(x)
!*****************************************
!*                 2.      00   -t*t     *
!*  erfc0(x) = ---------- S    e      dt *
!*              sqrt(pi)   x             *
!*****************************************
      IMPLICIT NONE
      REAL(DP) :: erfc0 
      REAL(DP) :: x
      if(x.lt.0.0_dp)then
       erfc0=1.0_dp+gammp(0.5_dp,x**2)   !20140620-7
       else
       erfc0=       gammq(0.5_dp,x**2)   !20140620-7
      endif
      return
      end function erfc0
!
      function erfc1(x)
!*********************************************
!*                  1.       00   -t*t/2     *
!*   erfc1(x) = ------------ S    e       dt *
!*             sqrt(2.*pi)   x               *
!*********************************************
      IMPLICIT NONE
      REAL(DP) :: erfc1 
      REAL(DP) :: x
      if(x.lt.0.0_dp)then
       erfc1=0.5_dp*(1.0_dp+gammp(0.5_dp,0.5_dp*x**2))
       else
       erfc1=0.5_dp*(       gammq(0.5_dp,0.5_dp*x**2))
      endif
      return
      end function erfc1
!
      function gammp(a,x)
      IMPLICIT NONE
      REAL(DP) :: gammp
      REAL(DP) :: a,x,gammcf,gamser,gln
      if(x.lt.0.0_dp)then   !20140620-7
         write(6,*)'bad x value in gammp ',x
         stop
      endif
      if(a.le.0.0_dp)then   !20140620-7
         write(6,*)'bad a value in gammp ',a
         stop
      endif
      if(x.lt.a+1.0_dp)then   !20140620-7
         call gser(gamser,a,x,gln)
         gammp=gamser
      else
         call gcf(gammcf,a,x,gln)
         gammp=1.0_dp-gammcf   !20140620-7
      endif
      return
      end function gammp
!
      function gammq(a,x)
      IMPLICIT NONE
      REAL(DP) :: gammq
      REAL(DP) :: x,a,gammcf,gamser,gln
      if(x.lt.0.0_dp)then   !20140620-7
         write(6,*)'bad x value in gammq ',x
         stop
      endif
      if(a.le.0.0_dp)then   !20140620-7
         write(6,*)'bad a value in gammq ',a
         stop
      endif
      if(x.lt.a+1.0_dp)then   !20140620-7
         call gser(gamser,a,x,gln)
         gammq=1.0_dp-gamser   !20140620-7
      else
         call gcf(gammcf,a,x,gln)
         gammq=gammcf
      endif
      return
      end function gammq
!
      subroutine gser(gamser,a,x,gln)
      IMPLICIT NONE
      INTEGER, PARAMETER :: itmax=100
      REAL(DP), PARAMETER :: eps=3.e-7
      REAL(DP) :: ap, sum,gln,a,x,gamser,del
      INTEGER  :: n
      gln=gammln(a)
      if(x.le.0.0_dp)then   !20140620-7
         if(x.lt.0.0_dp)then   !20140620-7
            write(6,*)'x < 0. in gser ',x
            stop
         endif
         gamser=0.0_dp   !20140620-7
         return
      endif
      ap=a
      sum=1.0_dp/a   !20140620-7
      del=sum
      do n=1,itmax
         ap=ap+1.0_dp   !20140620-7
         del=del*x/ap
         sum=sum+del
         if(abs(del).lt.abs(sum)*eps)goto 1
      enddo
      write(6,*)'a is too large, itmax must be greater gse'
 1    gamser=sum*exp(-x+a*log(x)-gln)
      return
      end subroutine gser
!
      subroutine gcf(gammcf,a,x,gln)
      IMPLICIT NONE
      INTEGER, PARAMETER :: itmax=100
      REAL(DP), PARAMETER :: eps=3.e-7,fpmin=1.e-30
      REAL(DP) :: b, c, d, h, an,a,x,gln,gammcf,del
      INTEGER  :: i
      gln=gammln(a)
      b=x+1.0_dp-a   !20140620-7
      c=1.0_dp/fpmin   !20140620-7
      d=1.0_dp/b   !20140620-7
      h=d
      do i=1,itmax
         an=-i*(i-a)
         b=b+2.0_dp   !20140620-7
         d=an*d+b
         if(abs(d).lt.fpmin)d=fpmin
         c=b+an/c
         if(abs(c).lt.fpmin)c=fpmin
         d=1.0_dp/d   !20140620-7
         del=d*c
         h=h*del 
         if(abs(del-1.0_dp).lt.eps) goto 1   !20140620-7
      enddo
      write(6,*)'a is too large, itmax must be greater gcf'
 1    gammcf=exp(-x+a*log(x)-gln)*h
      return
      end subroutine gcf
!
!     returns the value ln[gamma(xx)] for xx > 0.
      function gammln(xx)
      IMPLICIT NONE
      REAL(DP) :: gammln
      REAL(DP), DIMENSION (:) :: cof(6)
      REAL(DP) :: x,xx, stp, ser, y,tmp
      INTEGER  :: j
      cof(1) =  76.18009172947146_dp   !20140620-7
      cof(2) = -86.50532032941677_dp   !20140620-7
      cof(3) =  24.01409824083091_dp   !20140620-7
      cof(4) = -1.231739572450155_dp   !20140620-7
      cof(5) =  0.1208650973866179e-2
      cof(6) = -0.5395239384953e-5
      stp    =  2.5066282746310005_dp   !20140620-7
      ser    =  1.000000000190015_dp   !20140620-7
      x = xx
      y = x
      tmp = x+5.5_dp   !20140620-7
      tmp = (x+0.5_dp)*log(tmp)-tmp   !20140620-7
      do j=1,6
         y=y+1.0_dp   !20140620-7
         ser = ser+cof(j)/y
      enddo
      gammln=tmp+log(stp*ser/x)
      return
      end function gammln

!---------------------------------------------------------

END MODULE messy_airtraf_tools_ga
!END PROGRAM messy_airtraf_tools_ga

