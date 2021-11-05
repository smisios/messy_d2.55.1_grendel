! **********************************************************************
!
! SUBMODEL INTERFACE LAYER (SMIL) ROUTINES FOR MESSy SUBMODEL AirTraf
!
! Author : Volker Grewe, DLR-IPA, December 2011    (SMIL)
!          Hiroshi Yamashita, DLR-IPA, December 2011   (SMCL)
!
! References: see messy_airtraf.f90
!
! **********************************************************************
! Objectiv: Simulation of air traffic in EMAC
! Input:    Flightplan Citypairs, timetable, aircraft
! Output:   Flight trajectories and emissions
! Method:   
! 
!                    SMIL                                     SMCL
!
!
! |-airtraf_initialize 
! |          |-----------------------------------------airtraf_read_ctrl
! |          |                                                |- ngroutes, nwaypoints, nwaypoints_out
! |          |                                                |  option_traj_calc: 0=great circle
! |          |                                                |  output options 0: standard (details via channel) 
! |          |                                                |                 1: ac locations
! |          |                                                |  options max_ac_reached:  0: Stop  1:ignore ac
! |          |                                                |  options max_ac_out: 0:Stop 1:ignore 2:add to last 
! |          |                                                |  option: update_trajectory: updated trajectory every timeste
! |          |                                                |  option: lfeedback_NO, lfeedback_H2O
! |          |                                                |  ldaily_fp: usage of flightplan daily airtraf?
! |          |-airtraf_read_cpl
! |          |    flightpath_file=" ", var_names
! |          |    
! |          | Decomposition of trajectories
! |
! |-airtraf_init_memory
! |          | Allocate global emission fields on all CPUs
! |          | if (p_parallel_io)
! |          |    Read Flightplan 
! |          |    Allocate flightplan
! |          |    Resort optimal for decomposition
! |          |    Adapt time if necessary
! |          |  endif  
! |          |  broadcast ng_flightplan
! |          |  Decompose flightplan
! |          |------------------------------------------number_of_ac
! |          |  Allocate trajectories
! |          |  New channels output, output ac_locations and rerun; _gp for emissions
! |          |  New representations: 
! |          |       REPR_AIRTRAF_3D_AC_OUT, REPR_AIRTRAF_3D_AC_RERUN, REPR_AIRTRAF_1D_AC_FP
! |          |  New channel objects 
! |          |      - active = rerun
! |          |      - output+fp=output
! |          |      - Emission field NO, H2O, fuel use, distance
! |          |      - Ac location for film
! |
! |-airtraf_init_coupling
! |          | If (lfeedback_NO) get NO tracers
! |          | If (lfeedback_H2O) get H2O tracers
! |
! |-airtraf_global_end (originally, global_start, see20130902-5)
! |          | Get global wind field on each cpu:
! |          |    
! |          |-Do 
! |          |  check_departure
! |          |           
! |          |  if (no_dep) exit
! |          |    --------------------------------------prepare_trajectory
! |          |                                            - check if active ac is available
! |          |                                            - copy flightplan to active 
! |          |                                            - delete flight from flightplan or append in the end
! |          |    --------------------------------------calculate_trajectory
! |          |                                            - heart of everything
! |          |                                            - calculate emissions along flighttrac
! |          |-Enddo
! |          |- if (lupdate_traj) ----------------------update_trajectory
! |          |------------------------------------------fly_aircraft
! |          |                                            - simple temporal interpolation of precalculated values
! |          |                                            - calculate current emissions add to gl_emissions
! |          |-Do --------------------------------------check_arrival
! |          |                                          check_out_trajectory
! |          |-enddo
! |          |------------------------------------------fill_emission_field (Fill global 3D field with emissions from trajectories)
! |          | Decompose emission fields
!{|-airtraf_physc
!{|          |-update emissions to all NO tracers -> NO and H2O see lnox_si
! |-airtraf_global_end
! |          |- collect array information from cpus
! |          |- print maximum number of ac in air, ac in output.
! |-airtraf_free_memory
!===========================
!3 march 2011, V. Grewe, Oberpfaffenhofen
!21 nov 2011, updated
!02.12.2011, revision
!13.12.2011, ...
! 
!=========================================================================
!
! **********************************************************************
MODULE messy_airtraf_e5
! **********************************************************************

  ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
  USE messy_main_blather_bi,    ONLY: start_message_bi,end_message_bi,&
                                      error_bi,info_bi


  ! SMCL
  USE messy_airtraf                 ! Describe which data
 
  USE messy_main_constants_mem, ONLY: DP

  USE messy_airtraf_gc  !, ONLY: props_time,   Describe which data


!Yin_20170321+
#ifdef MESSYTENDENCY
  USE messy_main_tendency_bi,   ONLY: mtend_get_handle, &
                                      mtend_add_l,      &
                                      mtend_register
#endif
!Yin_20170321-



  IMPLICIT NONE
  INTRINSIC :: NULL
  PRIVATE

!Yin_20170321+
#ifdef MESSYTENDENCY
  INTEGER :: my_handle
#endif
!Yin_20170321-

  ! GLOBAL PARAMETERS
  INTEGER, PARAMETER                          :: len_city_code=3       ! City Code: MUC=Munich, JFK=NY
  INTEGER, PARAMETER                          :: len_desc=10           ! Not clear whether needed at all??

  !Yin_20170321+
  INTEGER                          :: nairtraf_nox_trac_gp 
  INTEGER , DIMENSION(:), POINTER  :: idt_list_teNOx => NULL()
  !Yin_20170321-

  ! CPL-NAMELIST PARAMETERS
  
  ! Module data 


!      
! 4 Data storages: 1) Data storage for active trajectories = currently flying aircraft 
!                        and inactive trajectories, basically depending on the
!                        number of available aircraft
!                  2) Data storage for completed trajectories = aircraft is landed
!                        for collecting data untill output flashes memory 
!                        Possibility to compress data to save memory
!                        basically depending on number of daily flights and length of output interval
!                  3) Data storage for flightplan.
!                  4) Data providing information on flightplan and aircraft for individual trajectories
!                        
! Part - 1 Physical Units 
!
! Active routes: (output only in rerun files)
! ac_routes    decomposed on CPUs: waypoints, properties, nlroutes,1  allocated in init_memory
! p3_ac_routes  pointer to 3D-field-ac_routes with new_channel_obj = organises output.
!
! Information on active routes: (output only in rerun files)
! ac_routes_desc Description of status, decomposed: 1, 1, nlroutes, number of status variables  allocated in init_memory
! p2_ac_routes_desc  pointer to 2D-Field on status information with new_channel_obj 
!                           D_lon=1, D_lat=2, D_time=3, A_lon=4, A_lat=5, A_time=6, AC_Type=7, Flight_num=8, 
!                           AC_status=9: notused=-1, inflight=1, arrived=2
!                           AC_pos=10, AC_pos_old=11

  REAL(dp), DIMENSION(:,:,:,:), POINTER                    :: ac_routes=>NULL()    &
!                                                                                   & ! Allocated in init_memory 
!                                                                                   & ! (waypoints, properties, ng_routes, 1)
                                                             ,ac_routes_out=>NULL()& ! (waypoints_out, properties, ng_routes_out, 1)
                                                             ,ac_fp_info_r=>NULL()   ! Information on flightplan and aircraft
  REAL(dp), DIMENSION(:,:,:,:), POINTER                    :: ac_fp_info_i=>NULL()   ! for individual trajectories
!HY 20120717  CHARACTER(len=len_city_code),DIMENSION(:,:,:,:), POINTER :: ac_fp_info_city=>NULL()
!HY 20120717  CHARACTER(len=len_desc),DIMENSION(:,:,:,:), POINTER      :: ac_fp_info_c=>NULL()

! 
  REAL(dp), DIMENSION(:,:,:,:), POINTER :: ac_routes_desc                ! status of ac_routes
 
! Indices for flightplan and information on active and output routes 
! INTEGER, PARAMETER      ::  flightplan_integer=2, AC_type=1, Flight_num=2
  INTEGER, PARAMETER      ::  flightplan_integer=0 !, AC_type=1, Flight_num=2
! INTEGER, PARAMETER      ::  flightplan_real=6,      & 
  INTEGER, PARAMETER      ::  flightplan_real=5,      & 
!                             D_lon=3, D_lat=4, D_time=5, A_lon=6, A_lat=7, A_time=8
                              D_lon=1, D_lat=2, D_time=3, A_lon=4, A_lat=5 !, A_time=8
! INTEGER, PARAMETER      ::  flightplan_char=2 , D_city=9, A_city=10
  INTEGER, PARAMETER      ::  flightplan_char=0 , total_arrival=6, corr_i_traj_out=7 ,corr_ifp=8    ! D_city=6, A_city=7
  INTEGER, PARAMETER      ::  n_routes_descr=11, ac_status=9, ac_pos=10, ac_pos_old=11
  REAL(DP), PARAMETER     ::  notused=-1._dp, inflight=1._dp, arrived=2._dp  , initial_arrival=0._dp ! 1= status 2=AC-type
  INTEGER, PARAMETER      ::  total_fp_input=flightplan_real+flightplan_integer+flightplan_char
!
! Flightplan_input file  -> via coupling namelist
  CHARACTER(len=128)                               :: flightplan_filename
  CHARACTER(len=10),dimension(total_fp_input)      :: flightplan_varname
  CHARACTER(len=5),dimension(total_fp_input), PARAMETER :: flightplan_vardescription=(/              &
                                               "D_lon", "D_lat", "Dtime", "A_lon", "A_lat" &
                                                /)
!                                              "ACtyp","Flnum",                                       &
!                                              "D_lon", "D_lat", "Dtime", "A_lon", "A_lat", "Atime" , &
!                                              "Dcity", "Acity" &
!                                               /)

!Yin_20170801+
  CHARACTER(len=20), DIMENSION(2) :: c_potcov
  CHARACTER(len=20), DIMENSION(2) :: c_atr20_o3
  CHARACTER(len=20), DIMENSION(2) :: c_atr20_ch4
  CHARACTER(len=20), DIMENSION(2) :: c_atr20_h2o
  CHARACTER(len=20), DIMENSION(2) :: c_atr20_contrail
  CHARACTER(len=20), DIMENSION(2) :: c_atr20_co2
!Yin_20170801-



  LOGICAL      :: lfeedback_NO
  LOGICAL      :: lfeedback_H2O
! 
  TYPE FP_TYPE
     REAL(dp), dimension(:), pointer :: r
     REAL(dp), dimension(:), pointer :: i
     CHARACTER(len=len_desc), dimension(:), pointer :: c
  END TYPE
! Part - 2 Pointer to physical units - decomposed fields allocated with new_chnnel_obj
  REAL(dp), DIMENSION(:,:,:),POINTER                   :: p3_ac_routes=>NULL(), p3_ac_routes_out=>NULL()  ! used in new_channel_obj
  REAL(dp), DIMENSION(:,:),POINTER                     :: p2_ac_fp_info_r=>NULL()
  REAL(dp), DIMENSION(:,:),POINTER                     :: p2_ac_fp_info_i=>NULL()
!HY 20120717  CHARACTER(len=len_city_code),DIMENSION(:,:),POINTER  :: p2_ac_fp_info_city=>NULL()
!HY 20120717  CHARACTER(len=len_desc),DIMENSION(:,:), POINTER      :: p2_ac_fp_info_c=>NULL()
  TYPE(FP_TYPE), DIMENSION(total_fp_input)             :: p1_ac_flightplan
  REAL(DP), DIMENSION(:,:),POINTER                     :: p2_ac_routes_desc=>NULL()

 !Yin_20170423
 ! CPL-NAMELIST PARAMETER
 ! LOGICAL                                              :: AC = .TRUE.
 ! TYPE(t_chaobj_cpl)                                   :: AC_SUBM
 !Yin_20170423  POINTER FOR COUPLED CHANNEL OBJECTS (contrail) 
  REAL(DP), DIMENSION(:,:,:), POINTER                  :: cpc => NULL()
 !Yin_20170423
!Yin_20170801+
  REAL(DP), DIMENSION(:,:,:), POINTER                  :: ATR20_o3 => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER                  :: ATR20_ch4 => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER                  :: ATR20_h2o => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER                  :: ATR20_contrail => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER                  :: ATR20_co2 => NULL()
!Yin_20170801-
! ========================================================================================
! 
! The 3 data storages need to have a new representation for decomposition
  INTEGER,SAVE    ::  REPR_AIRTRAF_3D_AC_OUT, REPR_AIRTRAF_3D_AC_RERUN, REPR_AIRTRAF_1D_AC_FP, &!REPR_AIRTRAF_2D_AC_FP, &
                      REPR_AIRTRAF_2D_AC_ic_RE,REPR_AIRTRAF_2D_AC_ir_RE,REPR_AIRTRAF_2D_AC_ii_RE,                     &
                     ! REPR_AIRTRAF_2D_AC_info_r_RERUN,REPR_AIRTRAF_2D_AC_info_i_RERUN,REPR_AIRTRAF_2D_AC_info_c_RERUN,&
!                     REPR_AIRTRAF_2D_AC_ROUTES_DESC 
!                     REPR_AIRTRAF_3D_AC_info_r_RERUN,REPR_AIRTRAF_3D_AC_info_i_RERUN,REPR_AIRTRAF_3D_AC_info_c_RERUN 
                      REPR_AIRTRAF_2D_AC_RDESC
! Local
!
 ! GLOBALIZED FIELDS VISIBLE ON ANY CPU
 ! VIA messy_main_data_bi; to be globalized
!  REAL, DIMENSION(:,:,:), POINTER :: uwind       => NULL() ! uwind: global field
!  REAL, DIMENSION(:,:,:), POINTER :: vwind       => NULL() ! vwind: global field
 ! Decomposed(local) fields -> Allocated via channel
!  REAL, DIMENSION(:,:,:), POINTER :: NOx_emis    => NULL() ! NOx emission in kg(N)/box/s
!  REAL, DIMENSION(:,:,:), POINTER :: H2O_emis    => NULL() ! H2O emission in kg(H2O)/box/s
!  REAL, DIMENSION(:,:,:), POINTER :: distance    => NULL() ! Flown distance in km/box
!  REAL, DIMENSION(:,:,:), POINTER :: fuel_use    => NULL() ! Used fuel in kg/box/s
  REAL(DP), DIMENSION(:,:,:), POINTER :: NOx_emis    => NULL() ! NOx emission in kg(N)/box/s
  REAL(DP), DIMENSION(:,:,:), POINTER :: H2O_emis    => NULL() ! H2O emission in kg(H2O)/box/s
  REAL(DP), DIMENSION(:,:,:), POINTER :: distance    => NULL() ! Flown distance in km/box
  REAL(DP), DIMENSION(:,:,:), POINTER :: fuel_use    => NULL() ! Used fuel in kg/box/s
  !Yin_20170423
  REAL(DP), DIMENSION(:,:,:), POINTER :: CPC_traj    => NULL() ! contrail distance
  REAL(DP), DIMENSION(:,:,:), POINTER :: ATR20_O3_traj    => NULL() ! ATR20 by ozone
  REAL(DP), DIMENSION(:,:,:), POINTER :: ATR20_CH4_traj    => NULL() ! ATR20 by CH4
  REAL(DP), DIMENSION(:,:,:), POINTER :: ATR20_H2O_traj    => NULL() ! ATR20 by H2O
  REAL(DP), DIMENSION(:,:,:), POINTER :: ATR20_CPC_traj    => NULL() ! ATR20 by contrail distance
  REAL(DP), DIMENSION(:,:,:), POINTER :: ATR20_CO2_traj    => NULL() ! ATR20 by contrail distance
  REAL(DP), DIMENSION(:,:,:), POINTER :: ATR20_TOT_traj    => NULL() ! total climate impact ATR20
  !Yin_20170423
  REAL(DP), DIMENSION(:,:,:), POINTER :: vervel_3d   => NULL() ! vertical velocity [Pa/s], 20130910-2,20130911-1,20130916-3
!Yin_20170321+
  REAL(DP), DIMENSION(:,:,:), POINTER :: teNOx  => NULL() ! NOx concenrtration tendency
!Yin_20170321-

  PUBLIC :: airtraf_initialize    ! initialize submodel
  PUBLIC :: airtraf_init_memory   ! request memory
!!$  PUBLIC :: airtraf_new_tracer ! define new tracers
  PUBLIC :: airtraf_init_coupling ! set pointers for coupling to BM and other SMs
!!$  PUBLIC :: airtraf_init_tracer! initilize tracers
  PUBLIC :: airtraf_global_end    ! entry point in time loop (all vectors)
  PUBLIC :: airtraf_physc         ! entry point in time loop (current vector)
  PUBLIC :: airtraf_write_output  ! 
  PUBLIC :: airtraf_free_memory   ! free allocated memory

  ! PRIVATE SUBROTINES
  PRIVATE :: airtraf_read_nml_cpl

CONTAINS

  ! ####################################################################
  ! PUBLIC SUBROUTINES
  ! ####################################################################

  ! ====================================================================
  SUBROUTINE airtraf_initialize

    ! ------------------------------------------------------------------
    ! This subroutine is used to
    ! - read (and broadcast) the CTRL-namelist,
    ! - read (and broadcast) the CPL-namelist,
    ! - perform the basic setup of the submodel.
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
   


 
    USE messy_main_mpi_bi,    ONLY: p_parallel_io, p_io, p_bcast, p_pe
    USE messy_main_tools,     ONLY: find_next_free_unit
    USE messy_main_transform_bi,   ONLY: get_dc_index

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'airtraf_initialize'
    INTEGER                     :: status ! error status
    INTEGER                     :: iou    ! I/O unit

    INTEGER, DIMENSION(:,:),   POINTER :: IDX => NULL()


    CALL start_message_bi(modstr,'INITIALISATION',substr)  ! log-output

    ! READ CTRL namelist
    IF (p_parallel_io) THEN                  ! read only on I/O-PE
       iou = find_next_free_unit(100,200)    ! find free I/O unit
       CALL airtraf_read_nml_ctrl(status, iou)  ! read CTRL-namelist
       ! terminate if error
       IF (status /= 0) CALL error_bi('Error in reading CTRL namelist',substr)
    END IF
    ! BROADCAST CTRL namleist entries from I/O-PE(probably PE=0) to ALL OTHER PEs
    !HY
    CALL p_bcast(nwaypoints,            p_io)   ! Please include all variabe from NML
    CALL p_bcast(ngroutes,              p_io)    
    CALL p_bcast(nwaypoints_out,        p_io)    
    CALL p_bcast(ngroutes_out,          p_io)  
    CALL p_bcast(option_traj_calc,      p_io) 
    CALL p_bcast(option_output,         p_io)  
    CALL p_bcast(option_max_ac_reached, p_io)   
    CALL p_bcast(lupdate_traj,          p_io)   
!   CALL p_bcast(lfeedback_NO,          p_io)  
!   CALL p_bcast(lfeedback_H2O,         p_io)   
    CALL p_bcast(ldaily_fp,             p_io)   


    ! READ CPL namelist
    IF (p_parallel_io) THEN                  ! read only on I/O-PE
       iou = find_next_free_unit(100,200)    ! find next free I/O unit
       CALL airtraf_read_nml_cpl(status, iou)  ! read CPL-namelist
       ! terminate if error
       IF (status /= 0) CALL error_bi('Error in reading CPL namelist',substr)
    END IF
    ! BROADCAST CPL namleist entries from I/O-PE to ALL OTHER PEs
    !HY
    CALL p_bcast(flightplan_filename,   p_io)   ! Please include all variabe from NML
    CALL p_bcast(flightplan_varname,    p_io)   
    CALL p_bcast(lfeedback_NO,          p_io)   
    CALL p_bcast(lfeedback_H2O,         p_io)   
!Yin_20170801+
    CALL p_bcast(c_potcov, p_io)
    CALL p_bcast(c_atr20_o3, p_io)
    CALL p_bcast(c_atr20_ch4, p_io)
    CALL p_bcast(c_atr20_h2o, p_io)
    CALL p_bcast(c_atr20_contrail, p_io)
    CALL p_bcast(c_atr20_co2, p_io)
!Yin_20170801-
    ! ### PERFORM INITIAL SETUP (CALL RESPECTIVE SMCL ROUTINE(S)) HERE
    ! 
    ! Define decomosition of trajectories
    ! for calculation of fligth trajectory
    call get_dc_index(ngroutes,idx)
    nlroutes=IDX(p_pe,2) - IDX(p_pe,1) + 1

    ! for calculation of output fligth trajectory

    call get_dc_index(ngroutes_out,idx)
    nlroutes_out=IDX(p_pe,2) - IDX(p_pe,1) + 1

    write(*,*) "PE=",p_pe," has routes/outputroutes",nlroutes,nlroutes_out

!Yin_20170321+
#ifdef MESSYTENDENCY
    my_handle=mtend_get_handle(modstr)
#endif
!Yin_20170321-

    CALL end_message_bi(modstr,'INITIALISATION',substr)  ! log-output

  END SUBROUTINE airtraf_initialize
  ! ====================================================================

!!$  ! ====================================================================
!!$  SUBROUTINE airtraf_new_tracer
!!$
!!$    ! ------------------------------------------------------------------
!!$    ! This subroutine is used to define new tracers. See
!!$    ! http://www.atmos-chem-phys.net/8/1677   (including supplement !)
!!$    ! for full documentation.
!!$    ! ------------------------------------------------------------------
!!$
!!$    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
!!$    USE messy_main_tracer_mem_bi, ONLY: GPTRSTR
!!$    USE messy_main_tracer_bi,     ONLY: tracer_halt
!!$    ! MESSy
!!$    USE messy_main_tracer,        ONLY: new_tracer, set_tracer &
!!$                                      , R_molarmass ! ,ON, OFF
!!$    USE messy_main_constants_mem, ONLY: MO
!!$
!!$    IMPLICIT NONE
!!$
!!$    ! LOCAL
!!$    CHARACTER(LEN=*), PARAMETER :: substr = 'airtraf_new_tracer'
!!$    INTEGER                     :: status
!!$    INTEGER                     :: idt_myO3   ! tracer index (identifier)
!!$
!!$    CALL start_message_bi(modstr,'TRACER DEFINITION',substr)  ! log-output
!!$
!!$    ! ### define new tracers here
!!$    CALL new_tracer(status, GPTRSTR, 'myO3', modstr, idx=idt_myO3)
!!$    CALL tracer_halt(substr, status)   ! terminate if error
!!$    CALL set_tracer(status, GPTRSTR, idt_myO3 &
!!$         , R_molarmass, r=3.0_dp*MO)
!!$    CALL tracer_halt(substr, status)   ! terminate if error
!!$
!!$    CALL end_message_bi(modstr,'TRACER DEFINITION',substr)  ! log-output
!!$
!!$  END SUBROUTINE airtraf_new_tracer
!!$  ! ====================================================================

!!$  ! ====================================================================
  SUBROUTINE airtraf_init_memory
  
    ! ------------------------------------------------------------------
    ! This subroutine is used to request memory for the submodel.
    ! The preferable method is to use "channel objects".
    ! Allocate your own memory, only if absolutely required.
    ! ------------------------------------------------------------------

    ! BMIL
    USE messy_main_channel_error_bi,   ONLY: channel_halt
    USE messy_main_channel_bi,         ONLY: DC_IX, GP_3D_MID     !, SCALAR
    USE messy_main_channel,            ONLY: new_channel, new_channel_object, &
                                             new_attribute
    USE messy_main_mpi_bi,             ONLY: p_parallel_io, p_pe, p_nprocs &
                                           , p_bcast, p_io
    USE messy_main_channel_repr,       ONLY: new_representation, AUTO
    USE messy_main_channel_dimensions, ONLY: new_dimension

    USE messy_main_grid_netcdf,        ONLY: t_ncvar, import_ncvar
     
    USE messy_main_grid_def_mem_bi,    ONLY: nlev, nproma, ngpblks 
 
    USE messy_main_transform_bi,  ONLY: get_dc_index, scatter_glix  !_1d
    USE messy_main_timer,         ONLY: lstart

!Yin_20170321+
    USE messy_main_tracer_mem_bi,   ONLY: GPTRSTR
    USE messy_main_tracer_tools_bi, ONLY: tracer_halt
    USE messy_main_tracer,          ONLY: get_tracer_list
    USE messy_main_constants_mem,   ONLY: STRLEN_MEDIUM
!Yin_20170321-

    IMPLICIT NONE

    TYPE (t_ncvar),DIMENSION(total_fp_input) :: myvar

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'airtraf_init_memory'
   ! CHARACTER(LEN=8)            :: ch
    INTEGER                     :: status

    INTEGER, DIMENSION(:,:),   POINTER :: IDX => NULL()

! Yin_20170321+
    INTEGER                                 :: jt
    CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(:), POINTER :: subnames => NULL()
! Yin_20170321-

  !--- Declare dimensions for the channels: -------------------------------------!

  ! LOCAL
     INTEGER                               :: i
    
   ! LOGICAL                               :: ldiag, lpost

   ! CHARACTER(LEN=1)                      :: char1
   ! CHARACTER(LEN=2)                      :: char2



    ! auxiliary pointer for data managment
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: mem => NULL()
    INTEGER                               :: DIMID_AC_routes,DIMID_AC_wayp,DIMID_AC_props,      &
                                             DIMID_AC_routes_out,DIMID_AC_wayp_out,             &
                                             DIMID_AC_routes_fp,DIMID_AC_INFO_c,DIMID_AC_INFO_r,&
                                             DIMID_AC_INFO_i,DIMID_AC_routes_desc
    INTEGER                               :: istart, j, ifp, ifp_id 


    
    TYPE (FP_TYPE), dimension(total_fp_input) :: p1_ac_fp_h
 
    CHARACTER(LEN=*), PARAMETER      :: channel_name=modstr//'_ac'

    ! Here, 'p_parallel_io' is .true.
    ! Parallel simulation is executed here
    IF(p_parallel_io) &
    CALL start_message_bi(modstr,'Channel definition',substr)
    

    CALL start_message_bi(modstr,'MEMORY ALLOCATION',substr)  ! log-output
    ! ### ALLOCATE OWN MEMORY HERE, BUT ONLY IF ABSOLUTELY REQUIRED!
    ! These allocation are probably performed for docmposed (local) array.
    ! Dimensions in waypoints and properties mustnot be decomposed; dimension in routes will be decomposed on each PE.
     
    ALLOCATE(ac_routes(nwaypoints,nprops,nlroutes,1))  
    ALLOCATE(ac_routes_out(nwaypoints_out,nprops,nlroutes_out,1))   
    ALLOCATE(ac_fp_info_r(1,nprops,nlroutes,1))
    ALLOCATE(ac_fp_info_i(1,nprops,nlroutes,1))
!HY 20120717    ALLOCATE(ac_fp_info_city(1,nprops,nlroutes,1))
!HY 20120717    ALLOCATE(ac_fp_info_c(1,nprops,nlroutes,1))
    ALLOCATE(ac_routes_desc(1,1,nlroutes,n_routes_descr))           !n_routes_descr = 11
    ALLOCATE(vervel_3d(nproma, nlev, ngpblks)) !20130911-1,20130916-3
    vervel_3d(:,:,:)=0.0_dp 
!    write(*,*)'PE=',p_pe,'check_vervel_3d_init_memory',shape(vervel_3d),vervel_3d(:,nlev,ngpblks)
! output check
!    write(*,*)'allocation check on PE=',p_pe
!    write(*,*)'nwaypoints=',nwaypoints,'nprops=',nprops,'nlroutes=',nlroutes
!    write(*,*)'nwaypoints_out=',nwaypoints_out,'nlroutes_out=',nlroutes_out
!    write(*,*)'n_routes_descr=',n_routes_descr


    CALL end_message_bi(modstr,'MEMORY ALLOCATION',substr)  ! log-output

    CALL start_message_bi(modstr,'READ DATA - SHOULD BE DONE IN IMPORT ...',substr)  ! log-output
   
    ! Here, 'p_parallel_io' is .true.
    ! i.e. parallel simulation is executed now, however, this if structure (below) is performed by
    ! only one PE. "IO_P" shows the PE. 
    if (p_parallel_io) then
!      write(*,*) "IO_P=",P_PE
      do ifp_id=1,total_fp_input
!      CALL IMPORT_NCVAR(myvar(i),var_name=TRIM(flightplan_varname(i)) &
!         ,file='/athome/pa2m/_data/AirTraf/FlightPlan/flightplan_1day.nc')
       CALL IMPORT_NCVAR(myvar(ifp_id),varname=TRIM(flightplan_varname(ifp_id)),file=flightplan_filename)


!       write(*,*) "   Variable read for ",flightplan_varname(ifp_id),": ",myvar(ifp_id)%name
!       write(*,*) "   Variable dimension:",myvar(ifp_id)%dat%dim(1)
!      write(*,*) "   First 10 data ",myvar(ifp_id)%dat%vr(1:10)

       ! HY
       if (ifp_id.le.flightplan_real)then
          write(*,*) "   First 10 data ",myvar(ifp_id)%dat%vr(1:10)
       else
          write(*,*) "   First 10 data ",myvar(ifp_id)%dat%vc(1:10)
       endif

       ng_flightplan=myvar(ifp_id)%dat%dim(1)
       WRITE(*,*) ""

       ! Make sure that data are read in the right order!
       if (ifp_id.le.flightplan_real) then
             ALLOCATE(p1_ac_fp_h(ifp_id)%r(ng_flightplan))
       elseif (ifp_id.le.flightplan_integer+flightplan_real) then
!            This elseif is not used now.
             ALLOCATE(p1_ac_fp_h(ifp_id)%i(ng_flightplan))
       else
             ALLOCATE(p1_ac_fp_h(ifp_id)%c(ng_flightplan))
       endif

       WRITE(*,*) "----> Interim flight plan allocated"

       !Re-organize the flightplan so that all CPU have similar work load
       istart=1
       j=1
       do ifp=1,ng_flightplan
          if (ifp_id<=flightplan_real) then
             p1_ac_fp_h(ifp_id)%r(j)=myvar(ifp_id)%dat%vr(ifp)       
          elseif (ifp_id<=flightplan_real+flightplan_integer) then
!            This elseif is not used now. 
             p1_ac_fp_h(ifp_id)%i(j)=real(myvar(ifp_id)%dat%vi(ifp),dp)
          else
             p1_ac_fp_h(ifp_id)%c(j)=myvar(ifp_id)%dat%vc(ifp)
          endif
          j=j+p_nprocs
          if (j>ng_flightplan) then 
            istart=istart+1
            j=istart
          endif
       enddo

      write(*,*) "ifp_id check do inner,PE",ifp_id,p_pe

      enddo
      write(*,*) "All data resorted for optimal decomposition"
      write(*,*) "ifp_id check do outer,PE",ifp_id,p_pe

    endif
    write(*,*)'p_nprocs=',p_nprocs
    write(*,*)'total_fp_input,ifp_id,ng_flightplan,PE are=',total_fp_input,ifp_id,ng_flightplan,p_pe   

    ! Parallel simulation is executed here, however, the if structure (above) is performed by
    ! only one 'PE'. The 'PE' only has information of ng_flightplan=1840 and ifp_id=6. 
    ! Thus ng_flightplan is broadcasted by 'p_bcast'.
    ! Simulation are decomposed from here !!!!
    CALL p_bcast(ng_flightplan, p_io)

    call get_dc_index(ng_flightplan,idx)
    nl_flightplan=IDX(p_pe,2) - IDX(p_pe,1) + 1
   
    write(*,*) "PE=",p_pe,"nl_flightplan=",nl_flightplan, IDX(p_pe,1),IDX(p_pe,2)
!   write(*,*) "PE=",p_pe,"pr_flightplan (total)=",p1_ac_fp_h(D_lon)%r(1:10)
!   write(*,*) "PE=",p_pe,"pr_flightplan (on pe)=",pr_flightplan(IDX(p_pe,1):IDX(p_pe,1)+10,A_lat,1,1)

    CALL end_message_bi(modstr,'READ DATA - Should be done via import ...',substr)  ! log-output

    ! CHANNEL AND CHANNEL OBJECTS
    CALL start_message_bi(modstr,'CHANNEL DEFINITION',substr)  ! log-output

    !--- 1) Construct the airtraf channels: --------------------------------------!
    !
    ! New dimension is defined here
    ! CALL NEW_DIMENSION (status, id, name, len)
    ! id    : dimension identifer
    ! name  : name of dimension
    ! len   : dimension length
    ! 
    ! for output
    CALL new_dimension(status, DIMID_AC_routes_out, 'AirTraf_routes_out', ngroutes_out)
    CALL channel_halt(substr, status)
    CALL new_dimension(status, DIMID_AC_wayp_out, 'AirTraf_waypoints_out', nwaypoints_out)
    CALL channel_halt(substr, status)

    ! for Flight Plan = fp
    CALL new_dimension(status, DIMID_AC_routes_fp, 'AirTraf_routes_fp',ng_flightplan)
    CALL channel_halt(substr, status)

    ! currently active routes
    CALL new_dimension(status, DIMID_AC_routes, 'AirTraf_routes', ngroutes)
    CALL channel_halt(substr, status)
    CALL new_dimension(status, DIMID_AC_wayp, 'AirTraf_waypoints', nwaypoints)
    CALL channel_halt(substr, status)
    CALL new_dimension(status, DIMID_AC_props, 'AirTraf_properties', nprops)
    CALL channel_halt(substr, status)
    CALL new_dimension(status, DIMID_AC_INFO_r, 'flightplan_info_real', flightplan_real)
    CALL channel_halt(substr, status)
    CALL new_dimension(status, DIMID_AC_INFO_i, 'flightplan_info_integer', flightplan_integer)
    CALL channel_halt(substr, status)
    CALL new_dimension(status, DIMID_AC_INFO_c, 'flightplan_info_char', flightplan_char)
    CALL channel_halt(substr, status)

    ! added from the modificatoin on 17.04.2012
    CALL new_dimension(status, DIMID_AC_routes_desc, 'AirTraf_routes_desc', n_routes_descr)
    CALL channel_halt(substr, status)
    write(*,*)'new_dimension is done PE=',p_pe

! output check
!    write(*,*)'ngroutes_out=',ngroutes_out,'nwaypoints_out=',nwaypoints_out,'ng_flightplan=',ng_flightplan
!    write(*,*)'ngroutes=',ngroutes,'nwaypoints=',nwaypoints,'nprops=',nprops
!    write(*,*)'flightplan_real=',flightplan_real,'flightplan_char=',flightplan_char
!    write(*,*)'flightplan_integer=',flightplan_integer,'n_routes_descr=',n_routes_descr




    ! Data structures are difined here 
    ! CALL NEW_REPRESENTATIONS (status, id, name, rank, link, dctype, dimension_ids, axis, 
    !                           ldimlen, dcindex)
    ! id           : representation identifer
    ! name         : name of representation
    ! rank         : dimensions of using data
    !                dimensions corresponds to the pointers of 'CALL NER_CHANNEL_OBJECT'
    ! link         : rank mapping
    !                ('x')used or ('-') unused 
    ! dctype       : type of parallel decomposition
    ! dimension_ids: vector of dimension identifers
    !                determined by 'CALL NEW_DIMENSION', which means max(global) dimension number of the rank
    ! axis         : geometry information
    ! ldimlen      : local(i.e. parallel decomposed)dimension length
    !                'AUTO' means no decomposition along the coressponding dimension
    ! dcindex      : dcindex defined which index should be decomposed in the using data.
    !                The value of dcindex is based on rank(1-4). 
    !                added from the modification on 20120509-7,8 

    CALL new_representation(status, REPR_AIRTRAF_3D_AC_OUT, 'REPR_AIRTRAF_3D_AC_OUT'       &
         , rank = 3, link = 'xxx-', dctype = DC_IX                                         &
!        , dimension_ids = (/ DIMID_AC_routes_out, DIMID_AC_wayp_out, DIMID_AC_props /)    &
         , dimension_ids = (/ DIMID_AC_wayp_out, DIMID_AC_props, DIMID_AC_routes_out /)    &
         , axis = 'NNN-'                                                                   &
         , ldimlen       = (/ AUTO, AUTO, nlroutes_out /)                                  & 
         , dcindex = 3                                                                     &
         )
!        , axis = '-NN-'
!        , output_order  = (/ 3,1,4,2 /)
    CALL channel_halt(substr, status)

    CALL new_representation(status, REPR_AIRTRAF_3D_AC_RERUN, 'REPR_AIRTRAF_3D_AC_RERUN'   &
         , rank = 3, link = 'xxx-', dctype = DC_IX                                         &
!        , dimension_ids = (/ DIMID_AC_routes, DIMID_AC_wayp, DIMID_AC_props /)            &
         , dimension_ids = (/ DIMID_AC_wayp, DIMID_AC_props, DIMID_AC_routes /)            &
         , axis = 'NNN-'                                                                   &
         , ldimlen       = (/ AUTO, AUTO, nlroutes /)                                      &
         , dcindex = 3                                                                     &
         )
!        , output_order  = (/ 3,1,4,2 /)
    CALL channel_halt(substr, status)

    ! HY define representations for _r _i _c
!   CALL new_representation(status, REPR_AIRTRAF_3D_AC_info_r_RERUN, 'REPR_AIRTRAF_3D_AC_INFO_r_RERUN'  &
!   CALL new_representation(status, REPR_AIRTRAF_2D_AC_info_r_RERUN, 'REPR_AIRTRAF_2D_AC_INFO_r_RERUN'  &
    CALL new_representation(status, REPR_AIRTRAF_2D_AC_ir_RE, 'REPR_AIRTRAF_2D_AC_Ir_RE'  &
!        , rank = 3, link = 'xxx-', dctype = DC_IX                                                      &
         , rank = 2, link = '-xx-', dctype = DC_IX                                                      &
!        , dimension_ids = (/ DIMID_AC_props, DIMID_AC_INFO_r, DIMID_AC_routes /)                       &
         , dimension_ids = (/ DIMID_AC_props, DIMID_AC_routes /)                                        &
!        , axis = 'NNN-'                                                                                &
         , axis = '-NN-'                                                                                &
!        , ldimlen       = (/ AUTO, AUTO, nlroutes /)                                                   &
         , ldimlen       = (/ AUTO, nlroutes /)                                                         &
         , dcindex = 3                                                                                  &
         )
!        , output_order  = (/ 3,1,4,2 /)
    CALL channel_halt(substr, status)

!   CALL new_representation(status, REPR_AIRTRAF_3D_AC_info_i_RERUN, 'REPR_AIRTRAF_3D_AC_INFO_i_RERUN'  &
!   CALL new_representation(status, REPR_AIRTRAF_2D_AC_info_i_RERUN, 'REPR_AIRTRAF_2D_AC_INFO_i_RERUN'  &
    CALL new_representation(status, REPR_AIRTRAF_2D_AC_ii_RE, 'REPR_AIRTRAF_2D_AC_Ii_RE'  &
!        , rank = 3, link = 'xxx-', dctype = DC_IX                                                      &
         , rank = 2, link = '-xx-', dctype = DC_IX                                                      &
!        , dimension_ids = (/ DIMID_AC_props, DIMID_AC_INFO_i, DIMID_AC_routes /)                       &
         , dimension_ids = (/ DIMID_AC_props, DIMID_AC_routes /)                                        &
!        , axis = 'NNN-'                                                                                &
         , axis = '-NN-'                                                                                &
!        , ldimlen       = (/ AUTO, AUTO, nlroutes /)                                                   &
         , ldimlen       = (/ AUTO, nlroutes /)                                                         &
         , dcindex = 3                                                                                  &
         )
!        , output_order  = (/ 3,1,4,2 /)
    CALL channel_halt(substr, status)

!   CALL new_representation(status, REPR_AIRTRAF_3D_AC_info_c_RERUN, 'REPR_AIRTRAF_3D_AC_INFO_c_RERUN'  &
!   CALL new_representation(status, REPR_AIRTRAF_2D_AC_info_c_RERUN, 'REPR_AIRTRAF_2D_AC_INFO_c_RERUN'  &
    CALL new_representation(status, REPR_AIRTRAF_2D_AC_ic_RE, 'REPR_AIRTRAF_2D_AC_Ic_RE'  &
!        , rank = 3, link = 'xxx-', dctype = DC_IX                                                      &
         , rank = 2, link = '-xx-', dctype = DC_IX                                                      &
!        , dimension_ids = (/ DIMID_AC_props, DIMID_AC_INFO_c, DIMID_AC_routes /)                       &
         , dimension_ids = (/ DIMID_AC_props, DIMID_AC_routes /)                                        &
!        , axis = 'NNN-'                                                                                &
         , axis = '-NN-'                                                                                &
!        , ldimlen       = (/ AUTO, AUTO, nlroutes /)                                                   &
         , ldimlen       = (/ AUTO, nlroutes /)                                                         &
         , dcindex = 3                                                                                  &
         )
!        , output_order  = (/ 3,1,4,2 /)
    CALL channel_halt(substr, status)

    CALL new_representation(status, REPR_AIRTRAF_1D_AC_FP, 'REPR_AIRTRAF_1D_AC_FP'         &
         , rank = 1, link = 'x---', dctype = DC_IX                                         &
         , dimension_ids = (/ DIMID_AC_routes_fp /)                                        &
         , axis = 'N---'                                                                   &
         , ldimlen       = (/ nl_flightplan /)                                             &
         , dcindex = 1                                                                     &
         )
!        , output_order  = (/ 3,1,4,2 /)
    CALL channel_halt(substr, status)

!   CALL new_representation(status, REPR_AIRTRAF_2D_AC_ROUTES_DESC, 'REPR_AIRTRAF_2D_AC_ROUTES_DESC'         &
    CALL new_representation(status, REPR_AIRTRAF_2D_AC_RDESC, 'REPR_AIRTRAF_2D_AC_RDESC'    &
         , rank = 2, link = '--xx', dctype = DC_IX                                          &
         , dimension_ids = (/ DIMID_AC_routes, DIMID_AC_routes_desc /)                      &
         , axis = '--NN'                                                                    &
!        , ldimlen       = (/ nl_flightplan, AUTO /)                                        &
         , ldimlen       = (/ nlroutes, AUTO /)                                             &
         , dcindex = 3                                                                      &
         )
!        , output_order  = (/ 3,1,4,2 /)
    CALL channel_halt(substr, status)
    write(*,*)'new_representation is done'


    ! New channels (such as large data box) are difined with the unique name 'cname' 
    ! CALL NEW_CHANNEL (status, cname, reprid) 
    ! cname        : name of channels 
    ! reprid       : identifer of default representation(data structure)
    !                the default representation(data structure) for all channel objects of this channel
    !                can be pre-set. However, this can be overwritten by CALL NEW_CHANNEL_OBJECT

    CALL new_channel(status, modstr//'_ac', reprid=REPR_AIRTRAF_3D_AC_OUT, lrestreq=.true.)
    CALL channel_halt(substr, status)
    
!   CALL new_channel(status, modstr//'_ac_rerun', reprid=REPR_AIRTRAF_3D_AC_RERUN)
    CALL new_channel(status, modstr//'_ac_rr', reprid=REPR_AIRTRAF_3D_AC_RERUN, lrestreq=.true.)
    CALL channel_halt(substr, status)
    write(*,*)'new_channels for output and rerun are defined'

! Here, simulations have been already decomposed now.

     ! HY
!    CALL new_channel(status, modstr//'_ac', reprid=REPR_AIRTRAF_2D_AC_info_r_RERUN)
!    CALL channel_halt(substr, status)

!    CALL new_channel(status, modstr//'_ac', reprid=REPR_AIRTRAF_2D_AC_info_r_RERUN)
!    CALL channel_halt(substr, status)

!    CALL new_channel(status, modstr//'_ac', reprid=REPR_AIRTRAF_2D_AC_ROUTES_DESC)
!    CALL channel_halt(substr, status)


    ! New channel objects (such as small data box) are difined with the unique name 'cname' in an existing channel.   
    ! CALL NEW_CHANNEL_OBJECT (status, cname, oname, p0(-p4), reprid, mem) 
    ! cname        : name of channel 
    ! oname        : name of channel object
    ! p0 to p4     : pointer to data array
    ! reprid       : identifer of representation(data structure)
    ! mem          : external data array
    !                mem is an optional pointer. With the mem, pre-allocated memory
    !                for the channel object can be specified.

    mem => ac_routes_out(:,:,1:nlroutes_out,1:1)
    CALL new_channel_object(status, modstr//'_ac', 'routes_out'                     &
            , p3=p3_ac_routes_out, reprid=REPR_AIRTRAF_3D_AC_OUT, mem=mem)
    CALL channel_halt(substr, status)

    mem => ac_routes(:,:,1:nlroutes,1:1)
!   CALL new_channel_object(status, modstr//'_ac_rerun', 'routes_rerun'             &
    CALL new_channel_object(status, modstr//'_ac_rr', 'routes_rerun'             &
            , p3=p3_ac_routes, reprid=REPR_AIRTRAF_3D_AC_RERUN, mem=mem)
    CALL channel_halt(substr, status)

    ! HY ADD p2_...r ...i ...c to rerunfile
    ! Define REPR_AIRTRAF_2D_AC_info_r_RERUN
    mem => ac_fp_info_r(1:1,:,1:nlroutes,1:1)
!   CALL new_channel_object(status, modstr//'_ac_rerun', 'routes_rerun_info_r'      &
    CALL new_channel_object(status, modstr//'_ac_rr', 'routes_rerun_info_r'      &
!          ,p2=p2_ac_fp_info_r, reprid=REPR_AIRTRAF_2D_AC_info_r_RERUN,mem=mem)
           ,p2=p2_ac_fp_info_r, reprid=REPR_AIRTRAF_2D_AC_ir_RE, mem=mem)
    CALL channel_halt(substr, status)

    p2_ac_fp_info_r(:,:)=-1._dp

    mem => ac_fp_info_i(1:1,:,1:nlroutes,1:1)
!   CALL new_channel_object(status, modstr//'_ac_rerun', 'routes_rerun_info_i'      &
    CALL new_channel_object(status, modstr//'_ac_rr', 'routes_rerun_info_i'      &
!          ,p2=p2_ac_fp_info_i, reprid=REPR_AIRTRAF_2D_AC_info_i_RERUN,mem=mem)
           ,p2=p2_ac_fp_info_i, reprid=REPR_AIRTRAF_2D_AC_ii_RE, mem=mem)
    CALL channel_halt(substr, status)

    p2_ac_fp_info_i(:,:)=-1._dp

!    mem => ac_fp_info_c(1:1,:,1:nlroutes,1:1)
!    CALL new_channel_object(status, modstr//'_ac_rerun', 'routes_rerun_info_c'      &
!          ,p2=p2_ac_fp_info_c, reprid=REPR_AIRTRAF_2D_AC_info_c_RERUN,mem=mem)
!    CALL channel_halt(substr, status)

    mem => ac_routes_desc(1:1,1:1,1:nlroutes,:)
!   CALL new_channel_object(status, modstr//'_ac_rerun', 'routes_rerun_desc'        &
    CALL new_channel_object(status, modstr//'_ac_rr', 'routes_rerun_desc'        &
!          ,p2=p2_ac_routes_desc, reprid=REPR_AIRTRAF_2D_AC_info_r_RERUN,mem=mem)
!          ,p2=p2_ac_routes_desc, reprid=REPR_AIRTRAF_2D_AC_ROUTES_DESC,mem=mem)
           ,p2=p2_ac_routes_desc, reprid=REPR_AIRTRAF_2D_AC_RDESC, mem=mem)
    CALL channel_halt(substr, status)


    do ifp_id=1,flightplan_real
      CALL new_channel_object(status, modstr//'_ac', 'routes_flightplan_'//flightplan_vardescription(ifp_id) &
           ,p1=p1_ac_flightplan(ifp_id)%r, reprid=REPR_AIRTRAF_1D_AC_FP)
      CALL channel_halt(substr, status)
!     call scatter_glix(p1_ac_flightplan(ifp)%r,p1_ac_flightplan(ifp)%r)
      CALL scatter_glix(p1_ac_fp_h(ifp_id)%r,p1_ac_flightplan(ifp_id)%r)  !(sendbuf, recevbuf)
!     DEALLOCATE(p1_ac_flightplan(ifp)%r) ; NULLIFY(p1_ac_flightplan(ifp)%r)
!
! Here, simulations have been already decomposed now.
!     write(*,*)">>>>>check0 PE=",p_pe,"p1_ac_flightplan(fast 10)=",p1_ac_flightplan(ifp_id)%r(1:10)
      write(*,*)"check0_PE=", p_pe, "p1_ac_flightplan_fast10=", p1_ac_flightplan(ifp_id)%r(1:MIN(SIZE( p1_ac_flightplan(ifp_id)%r),4))
      write(*,*)'check1 ifp_id and PE are =',ifp_id,p_pe,ASSOCIATED(p1_ac_fp_h(ifp_id)%r)
!
      if(p_pe.eq.0)then 
         write(*,*)'check2 ifp_id and PE are =',ifp_id,p_pe,ASSOCIATED(p1_ac_fp_h(ifp_id)%r)
         DEALLOCATE(p1_ac_fp_h(ifp_id)%r) ; NULLIFY(p1_ac_fp_h(ifp_id)%r)
      endif
    enddo

!HY20120920-2
!    CALL new_channel_object(status, modstr//'_ac_rr', 'i_traj_out' &
!         ,p0=i_traj_out, reprid=SCALAR)
!    CALL channel_halt(substr, status)
!    CALL new_channel_object(status, modstr//'_ac_rr', 'i_traj_dep' &
!         ,p0=i_traj_dep, reprid=SCALAR)
!    CALL channel_halt(substr, status)

! comment out '!'on 2012.0709, See note!
! '!!' should be commented out if '!' is delated

!    do ifp_id=1+flightplan_real, flightplan_integer+flightplan_real
!      CALL new_channel_object(status, modstr//'_ac', 'routes_flightplan_'//flightplan_vardescription(ifp_id) &
!           ,p1=p1_ac_flightplan(ifp_id)%i, reprid=REPR_AIRTRAF_1D_AC_FP)
!      CALL channel_halt(substr, status)
!!     call scatter_glix(p1_ac_flightplan(ifp)%i,p1_ac_flightplan(ifp)%i) 
!      call scatter_glix(p1_ac_fp_h(ifp_id)%i,p1_ac_flightplan(ifp_id)%i)  !(sendbuf, recevbuf)
!!     DEALLOCATE(p1_ac_flightplan(ifp)%i) ; NULLIFY(p1_ac_flightplan(ifp)%i)
!      if (p_pe.eq.0)then
!!        write(*,*)'ifp_id and PE are =',ifp_id,p_pe,ASSOCIATED(p1_ac_fp_h(ifp_id)%r)
!         DEALLOCATE(p1_ac_fp_h(ifp_id)%i) ; NULLIFY(p1_ac_fp_h(ifp_id)%i)
!      endif
!    enddo


!    do ifp=1+flightplan_integer+flightplan_real,flightplan_char+flightplan_integer+flightplan_real
!      CALL new_channel_object(status, modstr//'_ac', 'routes_flightplan_'//flightplan_vardescription(ifp) &
!            , p1=p1_ac_flightplan(ifp_id)%c, reprid=REPR_AIRTRAF_1D_AC_FP)
!      CALL channel_halt(substr, status)
!!     call scatter_glix(p1_ac_flightplan(ifp)%c,p1_ac_flightplan(ifp)%c)
!      call scatter_glix(p1_ac_fp_h(ifp)%c,p1_ac_flightplan(ifp)%c)
!!     DEALLOCATE(p1_ac_flightplan(ifp)%c) ; NULLIFY(p1_ac_flightplan(ifp)%c)
!      DEALLOCATE(p1_ac_fp_h(ifp)%c) ; NULLIFY(p1_ac_fp_h(ifp)%c)
!    enddo
!
!!   write(*,*) "PE=",p_pe,"p1_ac_flightplan (D_lon)=",p1_ac_flightplan(D_lon)%r(1:nwaypoints)
!    write(*,*) "PE=",p_pe,"p1_ac_flightplan_D_lon=",p1_ac_flightplan(D_lon)%r(1:nwaypoints)
!    write(*,*) "PE=",p_pe,"p1_ac_flightplan (D_lat)=",p1_ac_flightplan(D_lat)%r(1:nwaypoints)
!    write(*,*) "PE=",p_pe,"p1_ac_flightplan (D_time)=",p1_ac_flightplan(D_time)%r(1:nwaypoints)
!    write(*,*) "PE=",p_pe,"p1_ac_flightplan (A_lon)=",p1_ac_flightplan(A_lon)%r(1:nwaypoints)
!    write(*,*) "PE=",p_pe,"p1_ac_flightplan (A_lat)=",p1_ac_flightplan(A_lat)%r(1:nwaypoints)


    !--- 2) Construct channels for emissions, fuel, and distance: --------------------------------------!
    !
    ! New channels (such as large data box) are difined with the unique name 'cname' 
    ! CALL NEW_CHANNEL (status, cname, reprid) 
    ! cname        : name of channels 
    ! reprid       : identifer of default representation(data structure)
    !                the default representation(data structure) for all channel objects of this channel
    !                can be pre-set. However, this can be overwritten by CALL NEW_CHANNEL_OBJECT

!   CALL new_channel(status, modstr//'_gp', reprid=GP_3D_MID, lrestreq=.true.)
    CALL new_channel(status, modstr//'_gp', reprid=GP_3D_MID)
    CALL channel_halt(substr, status)
    write(*,*)'new_channel _gp is defined'


    ! New channel objects (such as small data box) are difined with the unique name 'cname' in an existing channel.   
    ! CALL NEW_CHANNEL_OBJECT (status, cname, oname, p0(-p4), reprid, mem) 
    ! cname        : name of channel 
    ! oname        : name of channel object
    ! p0 to p4     : pointer to data array
    ! reprid       : which is not needed to use now. The CALL NEW_CHANNEL(airtraf_gp) has already 'reprrid', 
    !                which is applied for all NEW_CHANNEL_OBJECTS in the CHANNEL(airtraf_gp)
    !
    ! New attributes for 'channel' and 'channel object' are defined.   
    ! CALL NEW_ATTRIBUTE (status, cname, oname, oaname, c) 
    ! cname        : name of channel 
    ! oname        : name of channel object
    ! oaname       : name of new channel object attribute
    ! c            : string value of attribute 

!   CALL new_channel_object(status, modstr//'_gp', 'NOx-Emissions', p3=NOx_emis)
    CALL new_channel_object(status, modstr//'_gp', 'NOx_Emissions', p3=NOx_emis)
    CALL channel_halt(substr, status)
!   CALL new_attribute(status,modstr//'_gp','long_name', c='Air Traffic NOx emissions')
    CALL new_attribute(status,modstr//'_gp','NOx_Emissions','long_name', c='Air Traffic NOx emissions')
    CALL channel_halt(substr, status)
!   CALL new_attribute(status,modstr//'_gp','units', c='kg(N)/box/s ')
    CALL new_attribute(status,modstr//'_gp','NOx_Emissions','units', c='kg(N)/box/s ')
    CALL channel_halt(substr, status)
    CALL info_bi('channel/object' //modstr//'_gp/'// 'ac_NOx_emis'//' was created')
    write(*,*)'new_channel_object NOx_Emissions is done'


!   CALL new_channel_object(status, modstr//'_gp', 'H2O-Emissions', p3=H2O_emis)
    CALL new_channel_object(status, modstr//'_gp', 'H2O_Emissions', p3=H2O_emis)
    CALL channel_halt(substr, status)
!   CALL new_attribute(status,modstr//'_gp','long_name', c='Air Traffic H2O emissions')
    CALL new_attribute(status,modstr//'_gp','H2O_Emissions','long_name', c='Air Traffic H2O emissions')
    CALL channel_halt(substr, status)
!   CALL new_attribute(status,modstr//'_gp','units', c='kg(H2O)/box/s ')
    CALL new_attribute(status,modstr//'_gp','H2O_Emissions','units', c='kg(H2O)/box/s ')
    CALL channel_halt(substr, status)
    CALL info_bi('channel/object' //modstr//'_gp/'// 'ac_H2O_emis'//' was created')
    write(*,*)'new_channel_object H2O_Emissions is done'


!   CALL new_channel_object(status, modstr//'_gp', 'Flown distance', p3=distance)
    CALL new_channel_object(status, modstr//'_gp', 'Flown_Distance', p3=distance)
    CALL channel_halt(substr, status)
!   CALL new_attribute(status,modstr//'_gp','long_name', c='Air Traffic flown distance')
    CALL new_attribute(status,modstr//'_gp','Flown_Distance','long_name', c='Air Traffic flown distance')
    CALL channel_halt(substr, status)
!   CALL new_attribute(status,modstr//'_gp','units', c='km/box/s ')
    CALL new_attribute(status,modstr//'_gp','Flown_Distance','units', c='km/box/s ')
    CALL channel_halt(substr, status)
    CALL info_bi('channel/object' //modstr//'_gp/'// 'ac_dist'//' was created')
    write(*,*)'new_channel_object Flown_Distance is done'


!   CALL new_channel_object(status, modstr//'_gp', 'Fuel Consumption', p3=fuel_use)
    CALL new_channel_object(status, modstr//'_gp', 'Fuel_Consumption', p3=fuel_use)
    CALL channel_halt(substr, status)
!   CALL new_attribute(status,modstr//'_gp','long_name', c='Air Traffic fuel consumption')
    CALL new_attribute(status,modstr//'_gp','Fuel_Consumption','long_name', c='Air Traffic fuel consumption')
    CALL channel_halt(substr, status)
!   CALL new_attribute(status,modstr//'_gp','units', c='kg(fuel)/box/s ')
    CALL new_attribute(status,modstr//'_gp','Fuel_Consumption','units', c='kg(fuel)/box/s ')
    CALL channel_halt(substr, status)
    CALL info_bi('channel/object' //modstr//'_gp/'// 'ac_fuel_use'//' was created')
    write(*,*)'new_channel_object Fuel_Consumption is done'

    !Yin_20170423
    CALL new_channel_object(status, modstr//'_gp', 'contrail_potcov',p3=CPC_traj)
    CALL channel_halt(substr, status)
    CALL new_attribute(status,modstr//'_gp','contrail_potcov','long_name', c='contrail potential coverage by airtraf')
    CALL channel_halt(substr, status)
    CALL new_attribute(status,modstr//'_gp','contrail_potcov','units', c='fraction')
    CALL channel_halt(substr, status)
    CALL info_bi('channel/object' //modstr//'_gp/'// 'ac_contrail_pc'//' was created')
    write(*,*)'new_channel_object contrail_potcov is done'
    !Yin_20170423

    !Yin_20170801+
    CALL new_channel_object(status, modstr//'_gp', 'ATR20O3_air',p3=ATR20_O3_traj)
    CALL channel_halt(substr, status)
    CALL new_attribute(status,modstr//'_gp','ATR20O3_air','long_name', c='ATR20 ozone by airtraf NOx')
    CALL channel_halt(substr, status)
    CALL new_attribute(status,modstr//'_gp','ATR20O3_air','units', c='K/box/s')
    CALL channel_halt(substr, status)
    CALL info_bi('channel/object' //modstr//'_gp/'// 'ac_atr20_o3_pc'//' was created')
    write(*,*)'new_channel_object ATR20O3_air is done'
    
    CALL new_channel_object(status, modstr//'_gp', 'ATR20CH4_air',p3=ATR20_CH4_traj)
    CALL channel_halt(substr, status)
    CALL new_attribute(status,modstr//'_gp','ATR20CH4_air','long_name', c='ATR20 methane by airtraf NOx')
    CALL channel_halt(substr, status)
    CALL new_attribute(status,modstr//'_gp','ATR20CH4_air','units', c='K/box/s')
    CALL channel_halt(substr, status)
    CALL info_bi('channel/object' //modstr//'_gp/'// 'ac_atr20_ch4_pc'//' was created')
    write(*,*)'new_channel_object ATR20CH4_air is done'

    CALL new_channel_object(status, modstr//'_gp', 'ATR20H2O_air',p3=ATR20_H2O_traj)
    CALL channel_halt(substr, status)
    CALL new_attribute(status,modstr//'_gp','ATR20H2O_air','long_name', c='ATR20 H2O by fuel burn')
    CALL channel_halt(substr, status)
    CALL new_attribute(status,modstr//'_gp','ATR20H2O_air','units', c='K/box/s')
    CALL channel_halt(substr, status)
    CALL info_bi('channel/object' //modstr//'_gp/'// 'ac_atr20_h2o_pc'//' was created')
    write(*,*)'new_channel_object ATR20H2O_air is done'

    CALL new_channel_object(status, modstr//'_gp', 'ATR20CPC_air',p3=ATR20_CPC_traj)
    CALL channel_halt(substr, status)
    CALL new_attribute(status,modstr//'_gp','ATR20CPC_air','long_name', c='ATR20 CPC by airtraf contrail')
    CALL channel_halt(substr, status)
    CALL new_attribute(status,modstr//'_gp','ATR20CPC_air','units', c='K/box/s')
    CALL channel_halt(substr, status)
    CALL info_bi('channel/object' //modstr//'_gp/'// 'ac_atr20_contrail_pc'//' was created')
    write(*,*)'new_channel_object ATR20CPC_air is done'
    
    CALL new_channel_object(status, modstr//'_gp', 'ATR20CO2_air',p3=ATR20_CO2_traj)
    CALL channel_halt(substr, status)
    CALL new_attribute(status,modstr//'_gp','ATR20CO2_air','long_name', c='ATR20 CO2 by airtraf contrail')
    CALL channel_halt(substr, status)
    CALL new_attribute(status,modstr//'_gp','ATR20CO2_air','units', c='K/box/s')
    CALL channel_halt(substr, status)
    CALL info_bi('channel/object' //modstr//'_gp/'//'ac_atr20_co2_pc'//' was created')
    write(*,*)'new_channel_object ATR20CO2_air is done'
    
    CALL new_channel_object(status, modstr//'_gp', 'ATR20TOT_air',p3=ATR20_TOT_traj)
    CALL channel_halt(substr, status)
    CALL new_attribute(status,modstr//'_gp','ATR20TOT_air','long_name', c='total ATR20 by airtraf')
    CALL channel_halt(substr, status)
    CALL new_attribute(status,modstr//'_gp','ATR20TOT_air','units', c='K/box/s')
    CALL channel_halt(substr, status)
    CALL info_bi('channel/object' //modstr//'_gp/'// 'ac_atr20_tot_pc'//' was created')
    write(*,*)'new_channel_object ATR20TOT_air is done'

    !Yin_20170801-
!Yin_20170321+
    CALL new_channel_object(status, modstr//'_gp', 'teNOx', p3=teNOx)
    CALL channel_halt(substr, status)
    CALL new_attribute(status,modstr//'_gp','teNOx','long_name', c='Tendency air Traffic NOx concerntration')
    CALL channel_halt(substr, status)
    CALL new_attribute(status,modstr//'_gp','teNOx','units', c='mol/mol/s ')
    CALL channel_halt(substr, status)
    CALL info_bi('channel/object' //modstr//'_gp/'// 'ac_teNOx'//' was created')
    write(*,*)'new_channel_object teNOx is done'
!Yin_20170321-

    CALL end_message_bi(modstr,'MEMORY ALLOCATION',substr)  ! log-output

!   nseg = gp_nseg
!   ALLOCATE(start(nseg,IRANK))
!   ALLOCATE(cnt(nseg,IRANK))
!   ALLOCATE(meml(nseg,IRANK))
!   ALLOCATE(memu(nseg,IRANK))
!   
!   start(:,:) = gp_start(:,:)
!   cnt(:,:) = gp_cnt(:,:)
!   meml(:,:) = gp_meml(:,:)
!   memu(:,:) = gp_memu(:,:)
!   
!   cnt(:,3) = nmod
!   memu(:,3) = nmod
!   
!   CALL set_representation_decomp(status, REPR_GMXE_4D_NMOD &
!        , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
!   CALL channel_halt(substr, status)
    
!   DEALLOCATE(start) ; NULLIFY(start)
!   DEALLOCATE(cnt)   ; NULLIFY(cnt)
!   DEALLOCATE(meml)  ; NULLIFY(meml)
!   DEALLOCATE(memu)  ; NULLIFY(memu)


   ! ### ADD MORE CHANNEL OBJECTS HERE
   if (lstart) then
      do i=1,nlroutes
!HY 20120831-2
!        p3_ac_routes(:,:,i)=real(p_pe,kind=dp)         ! real(p_pe*nlroutes*100+i,kind=dp)
         p3_ac_routes(:,:,i)=0.0_dp                     ! 
                                                        ! Add PE number as initial values
      enddo
!      do i=1,nlroutes_out
!        p3_ac_routes_out(i,:,:)=real(p_pe*nlroutes*100+i,kind=dp)
!      enddo
         p3_ac_routes_out(:,:,:)=0.0_dp                  ! Add zero as initial values. See 2012.5.9-4 memo
         write (*,*) "VGVG P3 initialised"
   endif
!Yin_20170321+
    CALL start_message_bi(modstr,'LOOKING FOR TRACERS',substr)

    CALL get_tracer_list(status,GPTRSTR,'NO',idt_list_teNOx,subnames)
    CALL tracer_halt(substr,status)
   ! write(*,*)'tracer list',idt_list_teNOx
    nairtraf_nox_trac_gp=SIZE(idt_list_teNOx)
   ! write(*,*)'tracer size',nairtraf_nox_trac_gp 
    IF(p_parallel_io) WRITE(*,*)'GP_TRACERS:'
    DO jt=1,nairtraf_nox_trac_gp
      IF (p_parallel_io) THEN
         IF(TRIM(subnames(jt))=='')THEN
           WRITE(*,*)'... NO'
         ELSE
           WRITE(*,*)'... NO_'//TRIM(subnames(jt))
         ENDIF
      ENDIF
#ifdef MESSYTENDENCY
!      CALL mtend_register(my_handle,mtend_id_tracer,idt=idt_list_teNOx(jt))
      CALL mtend_register(my_handle,idt_list_teNOx(jt))
#endif
    ENDDO

    IF(ASSOCIATED(subnames))DEALLOCATE(subnames)
    NULLIFY(subnames)

    CALL end_message_bi(modstr,'LOOKING FOR TRACER', substr)
!Yin_20170321-
  
   CALL end_message_bi(modstr,'CHANNEL DEFINITION',substr)  ! log-output

  END SUBROUTINE airtraf_init_memory
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE airtraf_init_coupling

    ! ------------------------------------------------------------------
    ! This soubroutine is used to set pointers
    ! (channel objects and/or tracers) for coupling to the 
    ! basemodel and to other submodels.
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel,          ONLY: get_channel_object
    ! Yin_20170423    
    !USE messy_main_mpi_bi,     ONLY: p_parallel_io     
    ! Yin_20170423
!!$    USE messy_main_tracer_bi, ONLY: tracer_halt
!!$    USE messy_main_tracer,    ONLY: get_tracer

!HY
!    USE messy_main_mpi_bi,          ONLY: p_pe
    IMPLICIT NONE
    ! Yin_20170423
    INTRINSIC :: TRIM
    ! Yin_20170423
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'airtraf_init_coupling'
    INTEGER                     :: status

    CALL start_message_bi(modstr,'COUPLING',substr)  ! log-output


    if (lfeedback_NO) then
        WRITE(*,*) "NO Feedback-Option not yet supported" 
    endif
    if (lfeedback_H2O) then
        WRITE(*,*) "H2O Feedback-Option not yet supported" 
    endif
    ! ### set pointers to channel objects here

!Yin_20170801+
    CALL get_channel_object(status,TRIM(c_potcov(1)),TRIM(c_potcov(2)),p3=cpc)
    CALL channel_halt(substr,status)

    CALL get_channel_object(status,TRIM(c_atr20_o3(1)),TRIM(c_atr20_o3(2)),p3=ATR20_o3)
    CALL channel_halt(substr,status)

    CALL get_channel_object(status,TRIM(c_atr20_ch4(1)),TRIM(c_atr20_ch4(2)),p3=ATR20_ch4)
    CALL channel_halt(substr,status)

    CALL get_channel_object(status,TRIM(c_atr20_h2o(1)),TRIM(c_atr20_h2o(2)),p3=ATR20_h2o)
    CALL channel_halt(substr,status)

    CALL get_channel_object(status,TRIM(c_atr20_contrail(1)),TRIM(c_atr20_contrail(2)),p3=ATR20_contrail)
    CALL channel_halt(substr,status)
    
    CALL get_channel_object(status,TRIM(c_atr20_co2(1)),TRIM(c_atr20_co2(2)),p3=ATR20_co2)
    CALL channel_halt(substr,status)

!Yin_20170801-

!   write(*,*)'init_coupling_PE=',p_pe,'ishape_cpc:',shape(cpc),cpc(:,:,:)
!   write(*,*)'init_coupling_PE=',p_pe,'ishape_ATR20_contrail:',shape(ATR20_contrail),ATR20_contrail(:,:,:)
!   write(*,*)'init_coupling_PE=',p_pe,'ishape_ATR20_co2:',shape(ATR20_co2),ATR20_co2(:,:,:)

    ! ### set pointers to tracers here
!!$    CALL get_tracer(status, ...)
!!$    CALL tracer_halt(substr, status)    ! terminate on error

    CALL end_message_bi(modstr,'COUPLING',substr)  ! log-output

     write(*,*)'calculation check, before end_init_coupling'

  END SUBROUTINE airtraf_init_coupling
  ! ====================================================================

!!$  ! ====================================================================
!!$  SUBROUTINE airtraf_init_tracer
!!$
!!$    ! ------------------------------------------------------------------
!!$    ! This subroutine is used to initialise tracers (via NCREGRID)
!!$    ! according to the REGRID-namelists in airtraf_t.nml.
!!$    ! For a full documention of NCREGRID see
!!$    ! http://www.atmos-chem-phys.net/6/3557  (including supplement !)
!!$    ! ------------------------------------------------------------------
!!$
!!$    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
!!$    USE messy_main_tracer_bi, ONLY: tracer_init
!!$
!!$    IMPLICIT NONE
!!$    
!!$    CALL tracer_init(modstr) ! initialise tracers via NCREGRID
!!$
!!$  END SUBROUTINE airtraf_init_tracer
!!$  ! ====================================================================

  ! ====================================================================
  SUBROUTINE airtraf_global_end

    USE messy_main_timer,           ONLY: delta_time, YEAR, MONTH, DAY, HOUR, MINUTE, SECOND, lstart, lresume
    USE messy_main_grid_def_mem_bi, ONLY: nproma, npromz, ngpblks, nlev &   ! 20130814-4,20130909-3
                                        , nlon, ngl
    USE messy_main_grid_def_bi,     ONLY: philon, philat, grvol              !philon_2d, philat_2d
    USE messy_main_data_bi,         ONLY: u_scb, v_scb,               &   ! scan_buffer
                                          geopot_3d, press_3d         &
                                        , t_scb , rho_air_dry_3d      &
                                        , geopoti_3d, pressi_3d   !vervel
                                          !Basically, physical values reffered by messy_main_data_bi are "decomposed field".

    USE messy_main_transform_bi   !,ONLY: trp_gpdc_gpgl, locate_in_decomp !This locate_in_decomp is used only in SMIL.20130512-4

    USE messy_main_mpi_bi,          ONLY: p_pe!, p_io

    USE messy_main_timer,           ONLY: gregor2julian, Julian_date_start
    USE messy_main_tools,           ONLY: nn_index  
    USE messy_main_constants_mem,   ONLY: g, OneDay,M_air,MN,MO!,R_gas,N_A     !Earth acceleration [m/s^2],M_air [g/mol]
 
    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'airtraf_global_end'
    INTEGER :: status
    ! BELOW!! need '=>NULL?' and 'REAL(DP)' see 20121121-10,20130813-1
   ! REAL    , DIMENSION(:,:,:), POINTER :: zu_scb, zv_scb             ! Decomposed wind fields

    !Global fields allocated in init_memory -> make decomposed(local) fields out of it.
!     REAL, DIMENSION(:,:,:), POINTER :: glob_NOx_emis    => NULL() ! NOx emission in kg(N)/box/s
!     REAL, DIMENSION(:,:,:), POINTER :: glob_H2O_emis    => NULL() ! H2O emission in kg(H2O)/box/s
!     REAL, DIMENSION(:,:,:), POINTER :: glob_distance    => NULL() ! flown distance in km/box
!     REAL, DIMENSION(:,:,:), POINTER :: glob_fuel_use    => NULL() ! Used fuel in kg/box/s
    REAL(DP), DIMENSION(:,:,:), POINTER   :: glob_NOx_emis    => NULL() ! NOx emission in g(N)/box/s
    REAL(DP), DIMENSION(:,:,:), POINTER   :: glob_H2O_emis    => NULL() ! H2O emission in g(H2O)/box/s
    REAL(DP), DIMENSION(:,:,:), POINTER   :: glob_distance    => NULL() ! flown distance in km/box
    REAL(DP), DIMENSION(:,:,:), POINTER   :: glob_fuel_use    => NULL() ! Used fuel in kg/box/s
   !Yin_20170321+
    REAL(DP), DIMENSION(:,:,:), POINTER   :: glob_teNOx_emis    => NULL() ! Used fuel in kg/box/s
   !Yin_20170321-

    !Yin_20170801+
    REAL(DP), DIMENSION(:,:,:), POINTER   :: glob_CPC         => NULL() ! Potcov in fraction
    REAL(DP), DIMENSION(:,:,:), POINTER   :: glob_ATR20O3     => NULL() ! ATR20 ozone in K/box/s
    REAL(DP), DIMENSION(:,:,:), POINTER   :: glob_ATR20CH4    => NULL() ! ATR20 methane in K/box/s
    REAL(DP), DIMENSION(:,:,:), POINTER   :: glob_ATR20H2O    => NULL() ! ATR20 h2o in K/box/s
    REAL(DP), DIMENSION(:,:,:), POINTER   :: glob_ATR20CPC    => NULL() ! ATR20 contrail in K/box/s
    REAL(DP), DIMENSION(:,:,:), POINTER   :: glob_ATR20CO2    => NULL() ! ATR20 CO2 in K/box/s
    REAL(DP), DIMENSION(:,:,:), POINTER   :: glob_ATR20TOT    => NULL() ! ATR20 total in K/box/s
    !Yin_20170801-
    REAL(DP), DIMENSION(:,:,:), POINTER   :: zgl_geopot_3d    => NULL() ! 
    REAL(DP), DIMENSION(:,:,:), POINTER   :: press_3d_g       => NULL() ! Pressure in Pa 
    REAL(DP), DIMENSION(:,:,:), POINTER   :: t_scb_g          => NULL() ! Temperature in K 
    REAL(DP), DIMENSION(:,:,:), POINTER   :: rho_air_dry_3d_g => NULL() ! Rho in kg/m3 
    REAL(DP), DIMENSION(:,:,:), POINTER   :: uwind_g          => NULL() ! uwind in m/s 
    REAL(DP), DIMENSION(:,:,:), POINTER   :: vwind_g          => NULL() ! vwind in m/s
    REAL(DP), DIMENSION(:,:,:), POINTER   :: v_z_g            => NULL() ! vertical wind in m/s 
    REAL(DP), DIMENSION(:,:,:), POINTER   :: v_z              => NULL() ! (decomposed)vertical wind in m/s 
  !Yin_20170801+
    REAL(DP), DIMENSION(:,:,:), POINTER   :: grvol_g          => NULL() ! grid volume, 20180526 
    REAL(DP), DIMENSION(:,:,:), POINTER   :: cpc_g            => NULL() ! Contrail potential coverage, 20170531-1 
    REAL(DP), DIMENSION(:,:,:), POINTER   :: ATR20_o3_g       => NULL() ! ATR20 ozone 
    REAL(DP), DIMENSION(:,:,:), POINTER   :: ATR20_ch4_g       => NULL() ! ATR20 methane 
    REAL(DP), DIMENSION(:,:,:), POINTER   :: ATR20_h2o_g       => NULL() ! ATR20 h2o 
    REAL(DP), DIMENSION(:,:,:), POINTER   :: ATR20_contrail_g  => NULL() ! ATR20 contrail distance 
    REAL(DP), DIMENSION(:,:,:), POINTER   :: ATR20_co2_g       => NULL() ! ATR20 CO2 
  !Yin_20170801-
!HY

! Local 
    LOGICAL                               :: departure
    REAL(DP)                              :: now
    REAL(DP), SAVE                        :: delta_time_julian          !,time_step_len_julian, 20160119-1  
!   REAL(DP)                              :: my_press  
    REAL(DP), DIMENSION(:), ALLOCATABLE   :: d_time_julian                   ! to check departure 
    INTEGER                               :: i_traj, jk, ifp, i_add, i_traj_notused, i_traj_inflight !, ifp_id
    INTEGER, SAVE                         :: i_traj_dep, i_traj_out     !20160119-1
    INTEGER                               :: i_traj_out_reuse, i_wp, kkproma, jjrow
    INTEGER                               :: d_pe, d_jp, d_jrow, d_jgx, d_jgy, d_idx 
    INTEGER                               :: i_add_NOx_emis      !see 20121117-1,2
    INTEGER                               :: i_add_H2O_emis
    INTEGER                               :: i_add_distance
    INTEGER                               :: i_add_fuel_use
    !Yin_20170423
    INTEGER                               :: i_add_CPC
    INTEGER                               :: i_add_ATR20O3
    INTEGER                               :: i_add_ATR20CH4
    INTEGER                               :: i_add_ATR20H2O
    INTEGER                               :: i_add_ATR20CPC
    INTEGER                               :: i_add_ATR20CO2
    INTEGER                               :: i_add_ATR20TOT
    !Yin_20170423
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: add_glob_NOx_emis   ! To add NOx emission into glob_NOx_emis
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: add_glob_H2O_emis   ! To add H2O emission into glob_H2O_emis
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: add_glob_distance   ! To add Flown distance into glob_distance
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: add_glob_fuel_use   ! To add Used fuel into glob_fuel_use
    !Yin_20170423
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: add_glob_CPC        ! To add CPC into glob_CPC
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: add_glob_ATR20O3   ! To add ATR20O3 into glob_ATR20O3
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: add_glob_ATR20CH4  ! To add ATR20CH4 into glob_ATR20CH4
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: add_glob_ATR20H2O  ! To add ATR20H2O into glob_ATR20H2O
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: add_glob_ATR20CPC  ! To add ATR20CPC into glob_ATR20CPC
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: add_glob_ATR20CO2  ! To add ATR20CO2 into glob_ATR20CO2
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: add_glob_ATR20TOT  ! To add ATR20TOT into glob_ATR20TOT
    !Yin_20170423
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: pass_phyval_along_traj   ! To pass physical values to SMCL
    REAL(DP), DIMENSION(:),   ALLOCATABLE :: zalt                ! geopotential altitude,20131105-1 
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: zdz               ! differences in (geopotential) altitude  
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: zdp               ! differences in pressure
    INTEGER                               :: ac_fl_direction     ! FL_DIR value, 20141217-1
!Viaualization of flight routes.20140610-2
!    INTEGER                               :: idx_vis, j_vis!,i_fl !20140610-2,20140611-1 
    INTEGER                               :: j_vis!,i_fl !20140610-2,20140611-1 
    REAL(DP)                              :: dmin_vis            !20140610-2
!Vialalization of wind field, 20140703-2     
    !INTEGER                               :: i,j,k
!    write(*,*)'calculation check, after global_end',"PE=",p_pe,"nlroutes=",nlroutes
!     do i_traj=1,10 !nlroutes
!             write(*,*) "restart_check_S_PE=",p_pe,"ac_routes(lon)=",ac_routes(:,1,i_traj,1)
!             write(*,*) "restart_check_S_PE=",p_pe,"ac_routes(lat)=",ac_routes(:,2,i_traj,1)
!             write(*,*) "restart_check_S_PE=",p_pe,"ac_routes(alt)=",ac_routes(:,3,i_traj,1)
!             write(*,*) "restart_check_S_PE=",p_pe,"ac_routes(time)=",ac_routes(:,4,i_traj,1)
!             write(*,*) "restart_check_S_PE=",p_pe,"ac_routes_out(lon)=",ac_routes_out(:,1,i_traj,1)
!             write(*,*) "restart_check_S_PE=",p_pe,"ac_routes_out(lat)=",ac_routes_out(:,2,i_traj,1)
!             write(*,*) "restart_check_S_PE=",p_pe,"ac_routes_out(alt)=",ac_routes_out(:,3,i_traj,1)
!             write(*,*) "restart_check_S_PE=",p_pe,"ac_routes_out(time)=",ac_routes_out(:,4,i_traj,1)
!     enddo
!    write(*,*)'shape_vervel_0',shape(vervel_3d)
    
    !GMT and GNUPLOT: global flight routes visualization.20140610-4
    !IF(p_pe.eq.0)OPEN(20,file='gmt.dat',position='append')

    ! ----------------------------------------------------------------------
    ! PART I: Get global winds fields
    ! ----------------------------------------------------------------------
    ! - Do this only once per time step.
    
    ! uwind and vwind
    ! This calculation is performed in GC, Wind_opt and Fuel_opt options.
    !IF(option_traj_calc==Wind_opt)THEN    !see20140524-10
!    write(*,*)"calculation check in PART1-1, PE=",p_pe

!    zu_scb => u_scb   !u_wind
!    zv_scb => v_scb   !v_wind
    
    ! See note 20121102-2,20120329-2,20120410-2
    ! Each PE has to call trp_gpdc_gpgl.
    ! CALL trp_gpdc_gpgl(sign, lf, gf, method)
    !    sign = 1         : decomposed field -> global field 
    !    lf               : (local)decomposed field = zu_scb, zv_scb 
    !    gf               : global field            = uwind, vwind 
    !    method(optional) :   
!    CALL trp_gpdc_gpgl(1, zu_scb, uwind)  
!    CALL trp_gpdc_gpgl(1, zv_scb, vwind)  
    ! After this CALL, each PE can refer uwind and vwind (global field).

    ! DIVIDE (U,V)*cos(lat) BY cos(lat) TO GET HORIZONTAL WIND

    ! Comment: to look at messy_gwave_si.f90 

!    zrcst(:,:) = 1./ sqcst(:,:)
!    DO jk = 1, nlev
!      uwind(:,jk,:) = uwind(:,jk,:) * zrcst(:,:)
!      vwind(:,jk,:) = vwind(:,jk,:) * zrcst(:,:)
!    END DO

    !Decomposed field: u_scb  :[m/s]
    !                  v_scb  :[m/s]
    !Global field    : uwind_g, vwind_g
    CALL trp_gpdc_gpgl(1, u_scb, uwind_g)  
    CALL trp_gpdc_gpgl(1, v_scb, vwind_g)  

    !Parallel decomposed:
    ! vertical wind
    !Decomposed field: zdz
    !                  v_z
    !                  zdp: difference in pressure 
    ALLOCATE(zdz(nproma, nlev, ngpblks))   
    ALLOCATE(v_z(nproma, nlev, ngpblks))    
    ALLOCATE(zdp(nproma, nlev, ngpblks))    
    !Initilize decomposed fields
    zdz(:,:,:) = 0.0_dp 
    v_z(:,:,:) = 0.0_dp 
    zdp(:,:,:) = 0.0_dp 

    !Calculate vertical depth of grid-box from difference in geopotential height.
    DO jk=1, nlev
       !geopoti_3d: geopotential at box interface.20130813-1
       !pressi_3d : pressure at box interface 20130814-3 
       !zdz       : decomposed geopotential altitude [m],20131105-1,20130814-3
       !zdp       : decomposed altitude difference [Pa],20130814-3
       zdz(:,jk,:)=(geopoti_3d(:,jk,:)-geopoti_3d(:,jk+1,:))/g  
       zdp(:,jk,:)=pressi_3d(:,jk+1,:)-pressi_3d(:,jk,:)        
    ENDDO
!    write(*,*)'PE=',p_pe,'vervel_3dcheck_1',vervel_3d(:,:,:)
    DO jjrow=1,ngpblks
       IF(jjrow==ngpblks)THEN
          kkproma=npromz
       ELSE
          kkproma=nproma
       ENDIF
!       write(*,*)'before_v_z'
!       write(*,*)'after_v_z'
       !Calculate vertical wind velocity [m/s]
       !vervel_3d: decomposed vertical wind velocity [Pa/s] at box center.20130814-1,3
       !v_z      : decomposed vertical wind [m/s]; positive(air dir.); negative(ground dir.)20130814-3
       v_z(1:kkproma,:,jjrow)=vervel_3d(1:kkproma,:,jjrow)/zdp(1:kkproma,:,jjrow) &
                              *zdz(1:kkproma,:,jjrow) 
       !still parallel decomposed!
    ENDDO

    !Make results globally available on all PEs.
    !v_z_g: global vertical wind velocity field [m/s]
    !       All PEs can refer wind values from v_z_g. 
    CALL trp_gpdc_gpgl(1, v_z, v_z_g)
    vervel_3d(:,:,:)=0.0_dp    !20130911-1,20130916-3
!    write(*,*)'check_vervel_3d_global_end'
    DEALLOCATE(zdz)
    DEALLOCATE(zdp)
    !ENDIF    !see20140524-10

    ! ----------------------------------------------------------------------
    ! PART II: Get global fields geopot_3d, press_3d, t_scb and rho_air_dry_3d
    ! ----------------------------------------------------------------------
    ! - See20121121-5,9,10
    !
    ! - Each PE has to call trp_gpdc_gpgl. Do this only once per time step.
    !
    ! - The Pointers, used by CALL trp_gpdc_gpgl, do not need ALLOCATION.    
    !   However, the Pointers need DEALLOCATE and NULLIFY. see20130813-1 
    !
    ! - Decomposed fields: geopot_3d       :[m^2/s^2]
    !                      press_3d        :[Pa]
    !                      t_scb           :[K]
    !                      rho_air_dry_3d  :[kg/m^3]
    !                      cpc             
    !                      grvol
    ! - Global fields    : zgl_geopot_3d
    !                      press_3d_g
    !                      t_scb_g
    !                      rho_air_dry_3d_g
    !                      cpc_g
    !                      grvol_g
      CALL trp_gpdc_gpgl(1, geopot_3d, zgl_geopot_3d)  !geopot_3d shows the value at box center. 
      CALL trp_gpdc_gpgl(1, press_3d, press_3d_g)
      CALL trp_gpdc_gpgl(1, t_scb, t_scb_g)
      CALL trp_gpdc_gpgl(1, rho_air_dry_3d, rho_air_dry_3d_g)
      CALL trp_gpdc_gpgl(1, grvol, grvol_g)            !20180526
    ! Yin_20170423 
      CALL trp_gpdc_gpgl(1, cpc, cpc_g)
      CALL trp_gpdc_gpgl(1, ATR20_o3, ATR20_o3_g)
      CALL trp_gpdc_gpgl(1, ATR20_ch4,ATR20_ch4_g)
      CALL trp_gpdc_gpgl(1, ATR20_h2o,ATR20_h2o_g)
      CALL trp_gpdc_gpgl(1, ATR20_contrail,ATR20_contrail_g)
      CALL trp_gpdc_gpgl(1, ATR20_co2,ATR20_co2_g)
    ! Yin_20170423
    ! After this CALL, each PE can refer global fields. 

    ! ----------------------------------------------------------------------
    ! PART III: Check Departures
    ! ----------------------------------------------------------------------
!    write(*,*)"calculation check in PART1-2, PE=",p_pe
!   write(*,*)"Check_hy:i_traj_dep=",i_traj_dep,"delta_time_julian=",delta_time_julian
    now = gregor2julian(YEAR, MONTH, DAY, HOUR, MINUTE, SECOND)  ! 'now' is defined in Julian_date

!    write(*,*)"now_partIII=",now,"PE=",p_pe
!    write(*,*)"julian_date_start=",Julian_date_start,"PE=",p_pe
!    write(*,*)"calculation check in PART1-3, PE=",p_pe

!    if(p_pe==3)write(*,*)'PE=',p_pe,now,'check_partII_geopot_3d_1',geopot_3d(1:nproma,1,ngpblks)
!    if(p_pe==3)write(*,*)'PE=',p_pe,now,'check_partII_geopot_3d_2',geopot_3d(1:npromz,1,ngpblks)
!    write(*,*)'npromz=',npromz,'nproma=',nproma
!    write(*,*)'shape_geopot_3d',shape(geopot_3d)
!    write(*,*)'shape_zgl_geopot_3d',shape(zgl_geopot_3d)
!    write(*,*)'PE=',p_pe,now,'shape_vervel_3d',shape(vervel_3d)
!    write(*,*)'PE=',p_pe,now,'shape_t_scb',shape(t_scb),t_scb(:,:,:)
!    write(*,*)'PE=',p_pe,now,'shape_t_scb_g',shape(t_scb_g),t_scb_g(:,:,:)   !HY20140522-2
!    write(*,*)'PE=',p_pe,now,'shape_cpc:',shape(cpc),cpc(:,:,:)
!    write(*,*)'PE=',p_pe,now,'shape_cpc_g:',shape(cpc_g),cpc_g(:,:,:)
!    write(*,*)'PE=',p_pe,now,'shape_ATR20_o3:',shape(ATR20_o3),ATR20_o3(:,:,:)
!    write(*,*)'PE=',p_pe,now,'shape_ATR20_o3_g:',shape(ATR20_o3_g),ATR20_o3_g(:,:,:)
!    write(*,*)'PE=',p_pe,now,'shape_ATR20_ch4:',shape(ATR20_ch4),ATR20_ch4(:,:,:)
!    write(*,*)'PE=',p_pe,now,'shape_ATR20_ch4_g:',shape(ATR20_ch4_g),ATR20_ch4_g(:,:,:)
!    write(*,*)'PE=',p_pe,now,'shape_ATR20_h2o:',shape(ATR20_h2o),ATR20_h2o(:,:,:)
!    write(*,*)'PE=',p_pe,now,'shape_ATR20_h2o_g:',shape(ATR20_h2o_g),ATR20_h2o_g(:,:,:)
!    write(*,*)'PE=',p_pe,now,'shape_ATR20_contrail:',shape(ATR20_contrail),ATR20_contrail(:,:,:)
!    write(*,*)'PE=',p_pe,now,'shape_ATR20_contrail_g:',shape(ATR20_contrail_g),ATR20_contrail_g(:,:,:)
!    write(*,*)'PE=',p_pe,now,'shape_ATR20_co2_g:',shape(ATR20_co2_g),ATR20_co2_g(:,:,:)
!    if(p_pe==3)write(*,*)'PE=',p_pe,now,'check_partII_zgl_geopot_3d',zgl_geopot_3d(:,:,:)
!    if(p_pe==3)write(*,*)'PE=',p_pe,now,'check_partII_press_3d',press_3d(:,:,:)
!    if(p_pe==3)write(*,*)'PE=',p_pe,now,'check_partII_g_press_3d',press_3d_g(:,:,:)
!    if(p_pe==3)write(*,*)'PE=',p_pe,now,'check_partII_t_scb',t_scb(:,:,:)
!    if(p_pe==3)write(*,*)'PE=',p_pe,now,'check_partII_global_t_scb_g',t_scb_g(:,:,:)
!    if(p_pe==3)write(*,*)'PE=',p_pe,now,'check_partII_rho',rho_air_dry_3d(:,:,:)
!    if(p_pe==3)write(*,*)'PE=',p_pe,now,'check_partII_global_rho_g',rho_air_dry_3d_g(:,:,:)

!    write(*,*) "AAtscb_PE=",p_pe,shape(t_scb_g),"tscb_check_t_scb_g_1=",t_scb_g(:,:,:)
!    write(*,*) "AArho_PE=",p_pe,shape(rho_air_dry_3d_g),"rho_check_g_1=",rho_air_dry_3d_g(:,:,:)
!    write(*,*) "AAuwind_PE=",p_pe,shape(uwind_g),"check_uwind_g_1=",uwind_g(:,:,:)
!    write(*,*) "AAATR20_contrail_PE=",p_pe,shape(ATR20_contrail_g),"ATR20_contrail_1=",ATR20_contrail_g(:,:,:)
!    write(*,*)'PE=',p_pe,now,'shape_grvol:',shape(grvol),'shape_grvol_g:',shape(grvol_g)

    IF((now-Julian_date_start).ge.1.0_dp)THEN
       write(*,*)"Julian_date_start=",Julian_date_start,"PE=",p_pe
       write(*,*)"now-Julian_date_start=",now-Julian_date_start,"PE=",p_pe
    ENDIF

! Initializations at very first step.
    if (lstart) then
       ! Initilizazion on p2_ac_routes_desc(i_traj,AC_status)
       do i_traj=1,nlroutes 
!         p2_ac_routes_desc(i_traj,AC_status)=notused
          p2_ac_routes_desc(i_traj,:)=notused 
       enddo

       ! Convert all p1_ac_flightplan(D_time)%r(ifp) to D_time_julian(ifp)(Julian_date)
       ALLOCATE(d_time_julian(nl_flightplan))
       do ifp=1, nl_flightplan
          d_time_julian(ifp) = p1_ac_flightplan(D_time)%r(ifp)/OneDay     !86400.0_DP 20140722-5
!          write(*,*)"PE=",p_pe
!          write(*,*)">> d_time_julian=",d_time_julian(ifp),"OLD p1_ac_flightplan(D_time)=",p1_ac_flightplan(D_time)%r(ifp)
          p1_ac_flightplan(D_time)%r(ifp) = now + d_time_julian(ifp)
!          write(*,*)">> now=",now,"NEW p1_ac_flightplan(D_time)=",p1_ac_flightplan(D_time)%r(ifp)
       enddo
       DEALLOCATE(d_time_julian)

       ! Convert delta_time to delta_time_julian (Julian_date)
       ! delta_time has always an identical value defined by timer.nml during global_end calculation.
!       delta_time_julian = REAL(delta_time,DP)/(60.0_DP * 60.0_DP * 24.0_DP)
!       write(*,*)"calculation check in PART1-4, PE=",p_pe
!       write(*,*)"time_step_len (sec)=",time_step_len
!       write(*,*)"delta_time (sec)=",delta_time,"delta_time_julian=",delta_time_julian
!!    endif

       ! Initilization of counting variables. 
       ! If performing more 1-day calu., we need modifications. See note 20120724-5.
!      i_traj_dep = 0
!      i_traj_out = 0
    endif

! Initializations at very first time step and at every restart.
    if (lstart .or. lresume) then
       ! Initialize delta_time_julian at very first step and at every restart.
       ! Convert delta_time to delta_time_julian (Julian_date)
       ! delta_time(= real(dp)) has always an identical value defined by timer.nml during global_end calculation.
       delta_time_julian = delta_time/OneDay     !86400.0_DP  20140722-5
!       write(*,*)"calculation check in PART1-4, PE=",p_pe
!       write(*,*)"time_step_len (sec)=",time_step_len
!       write(*,*)"delta_time (sec)=",delta_time,"delta_time_julian=",delta_time_julian

       ! Initilization of counting variables. See 20120724-5 
       ! First, i_traj_dep and i_traj_out become 0.
       ! If these values are needed, thay are calculated below line 1257
       i_traj_dep = 0
       i_traj_out = 0

       ! Check number of notused and inflight aircrafts.
       ! There is no aircraft in 'Arrival' status. When an aircraft arrived, 'Arrived' status is assinged and 
       ! after that 'notused' status is assigned in that arrived time-step.
       i_traj_notused  = 0
       i_traj_inflight = 0
       do i_traj=1,nlroutes
          if(p2_ac_routes_desc(i_traj,AC_status)==notused) i_traj_notused  = i_traj_notused  + 1
          if(p2_ac_routes_desc(i_traj,AC_status)==inflight)i_traj_inflight = i_traj_inflight + 1
       enddo
!       write(*,*)now,"Number of aircrafts(notused) (PE=",p_pe,"):","i_traj_notused =",i_traj_notused
!       write(*,*)now,"Number of aircrafts(inflight)(PE=",p_pe,"):","i_traj_inflight=",i_traj_inflight
    endif

! Initialize for CALL LOCATE_IN_DECOMP to obtain local and global fields indeces.
! Each PE has to do this initilization process.20121127-1,20130828-1,20140606-1
    if (lstart .or. lresume) then
       CALL locate_in_decomp(status, 0.0_dp, 0.0_dp, d_pe, d_jp, d_jrow) 
    endif

! Re-initialize and re-calculate number of departured/arrived aircrafs at every restart.
!   if (lresume) then
!   if (lresume .and. ((now-Julian_date_start).lt.2.0_dp)) then
!      write(*,*)'lresume_check',(lresume .and. ((now-Julian_date_start).le.2.0_dp))
!      i_traj_dep = 0
!      i_traj_out = 0

! Re-caucluate counting numbers of Departured/Arrived aircrafts on each PE.
! Julian_date_start (=real(dp))
    !Departured aircrafts(i_traj_dep): needed for only 1 day
    IF(lresume .and. ((now-Julian_date_start).lt.1.0_dp))THEN
!       write(*,*)'lresume_i_traj_dep',(lresume .and. ((now-Julian_date_start).lt.1.0_dp))
       i_traj_dep = 0
       do i_traj=1,nlroutes 
          if(p2_ac_routes_desc(i_traj,corr_ifp).gt.-1.0_dp)then
!         if(0.0_dp < ac_routes(1,4,i_traj,1) .and. ac_routes(1,4,i_traj,1) < now)then
             i_traj_dep = i_traj_dep + 1
          endif
!          write(*,*)"PE=",p_pe,"i_traj=",i_traj,"ac_routes(dep): ",ac_routes(1,4,i_traj,1),"< ",now
!          write(*,*)"PE=",p_pe,"i_traj=",i_traj,"p3_ac_routes(dep): ",p3_ac_routes(1,4,i_traj),"< ",now
       enddo
    ENDIF

    !Arrived aircrafts(i_traj_out): needed for only 2 days
    IF(lresume .and. ((now-Julian_date_start).lt.2.0_dp))THEN
!       write(*,*)'lresume_i_traj_out',(lresume .and. ((now-Julian_date_start).lt.2.0_dp))
       i_traj_out = 0
       do i_traj=1,nlroutes 
!         if(0.0_dp < ac_routes(nwaypoints,4,i_traj,1) .and. &
!            ac_routes(nwaypoints,4,i_traj,1) <= (now - delta_time_julian))then
          if(p2_ac_routes_desc(i_traj,corr_i_traj_out).gt.-1.0_dp)then
             i_traj_out = i_traj_out + 1
          endif
!          write(*,*)"PE=",p_pe,"i_traj=",i_traj,"ac_routes(arr):",ac_routes(nwaypoints,4,i_traj,1),"<= ",now-delta_time_julian
!          write(*,*)"PE=",p_pe,"i_traj=",i_traj,"p3_ac_routes(arr): ",p3_ac_routes(nwaypoints,4,i_traj),"<= ",now-delta_time_julian
       enddo
    ENDIF
!    write(*,*)'now=',now,'corr_ifp=',corr_ifp
!    write(*,*)now,"check:i_traj_dep=",i_traj_dep," i_traj_out=",i_traj_out," PE=",p_pe

! Dirty search all -> Change to pointer at current. ....
! HY20120724
!    i_traj_out=0
!    i_traj=0

!    now = gregor2julian(YEAR, MONTH, DAY, HOUR, MINUTE, SECOND)  ! 'now' is defined in Julian_date

!    write(*,*)"now=",now,"PE=",p_pe
!    write(*,*)"calculation check in PART1-4, PE=",p_pe

!! Convert all p1_ac_flightplan(D_time)%r(ifp) to D_time_julian(ifp) at first step (Julian_date)
!    if (lstart) then
!    ALLOCATE(d_time_julian(nl_flightplan))
!       do ifp=1, nl_flightplan
!          d_time_julian(ifp) = REAL(p1_ac_flightplan(D_time)%r(ifp),DP)/(60.0_DP * 60.0_DP * 24.0_DP) 
!          write(*,*)"PE=",p_pe
!          write(*,*)">> d_time_julian=",d_time_julian(ifp),"OLD p1_ac_flightplan(D_time)=",p1_ac_flightplan(D_time)%r(ifp)
!          p1_ac_flightplan(D_time)%r(ifp) = now + d_time_julian(ifp)
!          write(*,*)">> now=",now,"NEW p1_ac_flightplan(D_time)=",p1_ac_flightplan(D_time)%r(ifp)
!       enddo
!    DEALLOCATE(d_time_julian)
!    endif

! Convert delta_time to delta_time_julian (Julian_date)
!    time_step_len_julian = REAL(time_step_len,DP)/(60.0_DP * 60.0_DP * 24.0_DP)
!    delta_time_julian = REAL(delta_time,DP)/(60.0_DP * 60.0_DP * 24.0_DP)

!       write(*,*)"calculation check in PART1-4, PE=",p_pe
!       write(*,*)"time_step_len (sec)=",time_step_len
!       write(*,*)"delta_time (sec)=",delta_time,"delta_time_julian=",delta_time_julian

   !write(*,*)"restart_check:p2_ac_routes_desc=",p2_ac_routes_desc(1,1:11),"PE=",p_pe
   !write(*,*)"restart_check:p2_ac_routes_desc=",p2_ac_routes_desc(2,1:11),"PE=",p_pe
   !write(*,*)"restart_check:p2_ac_routes_desc=",p2_ac_routes_desc(3,1:11),"PE=",p_pe
   !write(*,*)"restart_check:p2_ac_routes_desc=",p2_ac_routes_desc(4,1:11),"PE=",p_pe
   !write(*,*)"restart_check:p2_ac_routes_desc=",p2_ac_routes_desc(5,1:11),"PE=",p_pe
   !write(*,*)"restart_check:p2_ac_routes_desc=",p2_ac_routes_desc(231,AC_status),"PE=",p_pe


    ALLOCATE(zalt(nlev))                             !Array
    ALLOCATE(pass_phyval_along_traj(nwaypoints-1,3)) !Array: 1:pressure, 2:temperature, 3:rho 
    DO ifp=1,nl_flightplan

!       write(*,*)"calculation check in PART1-5, PE=",p_pe
!       write(*,*)"nl_flightplan=",nl_flightplan

! Departure check !
! The departure check is applied for each aircraft only once for one day simulation.
!       departure=(       p1_ac_flightplan(ifp_id)%r(D_time)>=now    &
!                   .and. p1_ac_flightplan(ifp_id)%r(D_time)<now+time_step_len)                                  
! HY 20120718       
!       departure=(       p1_ac_flightplan(D_time)%r(ifp)>=now    &
!                   .and. p1_ac_flightplan(D_time)%r(ifp)<now+time_step_len)                                  
! HY 20120720
!       departure=(       p1_ac_flightplan(D_time)%r(ifp)>=now    &
!                   .and. p1_ac_flightplan(D_time)%r(ifp)<now+time_step_len_julian)                                  

       departure=(       p1_ac_flightplan(D_time)%r(ifp)>=now    &
                   .and. p1_ac_flightplan(D_time)%r(ifp)<now+delta_time_julian)                                  
!       write(*,*)"DEPARTURE_CHECK: ",departure,"PE= ",p_pe,"ifp=",ifp
!       write(*,*)"now<= p1_ac_flightplan(D_time)< now+delta_time_julian"
!       write(*,*)now,"<= ",p1_ac_flightplan(D_time)%r(ifp),"< ",now+delta_time_julian
!       write(*,*)"delta_time_julian=",delta_time_julian

       if (.not.departure) cycle   ! should I add the if condition, like ".and. =.not.inflight" ?

!       write(*,*)"calculation check in PART1-6, PE=",p_pe
       outer: do 
!HY20120724
!         i_traj=i_traj+1
!         i_traj_dep = i_traj_dep + 1

          if((now-Julian_date_start).lt.1.0_dp)then
             i_traj_dep = i_traj_dep + 1
!             write(*,*)'i_traj_dep_check1',i_traj_dep, 'now=',now,'PE=',p_pe
          elseif((now-Julian_date_start).ge.1.0_dp)then
             inner: do i_traj=1,nlroutes 
               !if(IDINT(p2_ac_routes_desc(i_traj,corr_ifp))==ifp)then  !20170212-2
                if(INT(p2_ac_routes_desc(i_traj,corr_ifp))==ifp)then
                   i_traj_dep = i_traj 
!                   write(*,*)'i_traj_dep_check2',i_traj_dep,'i_traj=',i_traj, 'now=',now,'PE=',p_pe
!                   write(*,*)'idint_check',INT(p2_ac_routes_desc(i_traj,corr_ifp)),'=',ifp
                  !write(*,*)'idint_check',IDINT(p2_ac_routes_desc(i_traj,corr_ifp)),'=',ifp   !20170212-2
                   exit inner !leave loop
                endif 
             enddo inner
          endif
!          write(*,*)"Number of departured aircrafts(PE=",p_pe,"):","i_traj_dep=",i_traj_dep
!HY20120724
!         if (i_traj.gt.nlroutes) then
          if (i_traj_dep.gt.nlroutes) then
                 ! Do whatever necessary to solve this problem. 
                 ! Depending on parameters: Stop run
                 ! cancel start
                 ! see note 2012.03.30-2 
          endif
!          write(*,*)"nlroutes=",nlroutes
!          if (p2_ac_routes_desc(i_traj,AC_status)==notused) then
!              p2_ac_routes_desc(i_traj,AC_status)=inflight        ! free trajectory found
!              p2_ac_routes_desc(i_traj,D_time)=p1_ac_flightplan(ifp_id)%r(D_time)
!              p2_ac_routes_desc(i_traj,D_lon)=p1_ac_flightplan(ifp_id)%r(D_lon)
!              p2_ac_routes_desc(i_traj,D_lat)=p1_ac_flightplan(ifp_id)%r(D_lat)              
!              p2_ac_routes_desc(i_traj,A_lon)=p1_ac_flightplan(ifp_id)%r(A_lon)              
!              p2_ac_routes_desc(i_traj,A_lat)=p1_ac_flightplan(ifp_id)%r(A_lat)
!              p2_ac_routes_desc(i_traj,AC_type)=real(p1_ac_flightplan(ifp_id)%i(AC_type),dp)
!              exit     ! leave loop 
!           endif   
!          write(*,*)"AC_status=",AC_status,"D_time",D_time,"D_lon",D_lon,"D_lat",D_lat,"A_lon",A_lon,"A_lat",A_lat
!          write(*,*)"notused=",notused,"inflight=",inflight
!          write(*,*)"p2_ac_routes_desc=",p2_ac_routes_desc(i_traj_dep,AC_status)
! HY20120724-6: In if structure below, all "i_traj" were changed into "i_traj_dep".
          if (p2_ac_routes_desc(i_traj_dep,AC_status)==notused) then
              p2_ac_routes_desc(i_traj_dep,AC_status)    = inflight           ! free trajectory found
!             p2_ac_routes_desc(i_traj_dep,total_arrival)= initial_arrival    ! used for first calc. 
!             p2_ac_routes_desc(i_traj_dep,corr_ifp)     = ifp                ! used for first calc.
              p2_ac_routes_desc(i_traj_dep,D_lon)        = p1_ac_flightplan(D_lon)%r(ifp)
              p2_ac_routes_desc(i_traj_dep,D_lat)        = p1_ac_flightplan(D_lat)%r(ifp)              
              p2_ac_routes_desc(i_traj_dep,D_time)       = p1_ac_flightplan(D_time)%r(ifp)! Julian_data
              p2_ac_routes_desc(i_traj_dep,A_lon)        = p1_ac_flightplan(A_lon)%r(ifp)              
              p2_ac_routes_desc(i_traj_dep,A_lat)        = p1_ac_flightplan(A_lat)%r(ifp)
!             p2_ac_routes_desc(i_traj,AC_type)=p1_ac_flightplan(AC_type)%i(ifp)
!             p2_ac_routes_desc(i_traj_dep,AC_POS)=1.0_dp       !20120822-2, 20120829-1
!             p2_ac_routes_desc(i_traj_dep,AC_POS_OLD)=1.0_dp

              !Used in first calc.
              IF(p2_ac_routes_desc(i_traj_dep,total_arrival)== -1.0_dp)THEN
                 p2_ac_routes_desc(i_traj_dep,total_arrival)= initial_arrival   
                 p2_ac_routes_desc(i_traj_dep,corr_ifp)     = REAL(ifp,DP)              
!              write(*,*)'ifp_check=',REAL(ifp,DP),'now=',now,'p_pe=',p_pe
              ENDIF 
              !Add 1 day to flightplan(D_time) for tomorrow's flight, see 20121220-4
              p1_ac_flightplan(D_time)%r(ifp) = p1_ac_flightplan(D_time)%r(ifp) + 1.0_dp

              exit outer    ! leave loop 
          endif   
       enddo outer
!VG
! Call prepare_trajectory_and_flightplan.

    !option_traj_calc check!
!    write(*,*)'0_PE=',p_pe,'option_traj_calc=',option_traj_calc
!    write(*,*)'0_PE=',p_pe,'GC=',GC,'Wind_opt=',Wind_opt,'Fuel_opt=',Fuel_opt &
       !,'NOx_opt=',NOx_opt,'H2O_opt=',H2O_opt,'Cost_opt=',Cost_opt
!    write(*,*)'0_PE=',p_pe,'COC_opt=',COC_opt
    !Yin_20170423+
!    write(*,*)'0_PE=',p_pe,'ContrailPC_opt=',ContrailPC_opt,'ATR20_opt='&
       !,ATR20_opt,'COSTCLIM_opt=',COSTCLIM_opt,'COSTCPC_opt=',COSTCPC_opt
    !Yin_20170423-

    ! ----------------------------------------------------------------------
    ! This subroutine determines projected flight routes, e.g., waypoints.
    ! p2_ac_routes_desc(:,:) expresses an identical flight(aircraft), i.e.
    ! ac_routes(). 
    !
    ! Note: - p2_ac_routes_desc(:,D_time) is defined in Julian_date.
    !       - i_traj_dep indicates an order of departure. 20120103-1
    !       - For Fuel_opt option, fuel use calculation is performed in calculate_trajectory, 20160301-1 
    ! ----------------------------------------------------------------------
!HY 20120724-6
!      call calculate_trajectory(ac_routes(:,:,i_traj,1),   &        ! 
!                                 p2_ac_routes_desc(i_traj,:))       ! 
!!     call calculate_trajectory (ac_routes(:,:,:,1),Wind, ...)    :for example
       !write(*,*)'PE=',p_pe,'philon=',philon
       !write(*,*)'PE=',p_pe,'philat=',philat
       IF(option_traj_calc==GC)THEN                                     ! 20130904-3
          ! Is the ac_routes maybe 'p3_ac_routes'?
          CALL calculate_trajectory(ac_routes(:,:,i_traj_dep,1),    &   
               !Route description, note 20120724-2,2012.4.3,20111129-3
                                    p2_ac_routes_desc(i_traj_dep,:),&   
                                    philon(:),philat(:),zgl_geopot_3d(:,:,:),&    !added more arrays, see20140524-10 
                                    uwind_g(:,:,:),vwind_g(:,:,:),v_z_g(:,:,:), &
                                    t_scb_g(:,:,:),ac_fl_direction,cpc_g(:,:,:))     ! ac_fl_direction20141217-1
          ! The GC trajectory is found.
       ELSEIF(option_traj_calc==Wind_opt .or. option_traj_calc==ContrailPC_opt)THEN  !20170608-3
          ! Pass global wind fields to optimization routine in SMCL.
          ! philon(:) and plilat(:) are geographical coordinates of global field.20130812-7
          !    philon(64(=nlon)): 0..->>..+360[deg]
          !    philat(32(=ngl)) : +90..->>..-90[deg]
          ! v_z_g(:,:,:)        : global vertical wind velocity field [m/s] 
          ! Yin_20170423 ContrailPC option is added
          CALL calculate_trajectory(ac_routes(:,:,i_traj_dep,1),    &   ! 20130904-3,20130909-2 
                                    p2_ac_routes_desc(i_traj_dep,:),& 
                                    philon(:),philat(:),zgl_geopot_3d(:,:,:),&   
                                    uwind_g(:,:,:),vwind_g(:,:,:),v_z_g(:,:,:), &
                                    t_scb_g(:,:,:),ac_fl_direction,cpc_g(:,:,:))     !t_scb_g 20140522-2, ac_fl_direction20141217-1 
          ! The time/CPC optimal trajectory is found.                                        
          ! tecplot output 20140704-1, 20141218-3
!!         IF((p_pe.eq.0).and.(i_traj_dep.eq.291))THEN 
!          IF((p_pe.eq.0).and.(i_traj_dep.eq.20))THEN       !MUC >> JFK 
!             call system('mv /data/yama_hi/EMAC/AirTraf/messy_2.41/workdir/tecplot_muc_jfk.dat alltraj_muc_jfk.dat')
!             call system('mv /data/yama_hi/EMAC/AirTraf/messy_2.41/workdir/history_muc_jfk.dat converge_muc_jfk.dat')
!          ELSEIF((p_pe.eq.3).and.(i_traj_dep.eq.1))THEN    !JFK >> MUC
!             call system('mv /data/yama_hi/EMAC/AirTraf/messy_2.41/workdir/tecplot_jfk_muc.dat alltraj_jfk_muc.dat')
!             call system('mv /data/yama_hi/EMAC/AirTraf/messy_2.41/workdir/history_jfk_muc.dat converge_jfk_muc.dat')
!!         ELSE     !20141218-3
!!            call system('rm /data/yama_hi/EMAC/AirTraf/messy_2.41/workdir/tecplot.dat')
!!            call system('rm /data/yama_hi/EMAC/AirTraf/messy_2.41/workdir/converge.dat')
!          ENDIF
       ELSEIF(option_traj_calc==Fuel_opt .or. option_traj_calc==H2O_opt .or.  &   !20160301-1,20170222-1 
              option_traj_calc==Cost_opt .or. option_traj_calc==COC_opt .or. &
              option_traj_calc==COSTCPC_opt)THEN      !20170224-1,20170316-1
          ! Pass global wind fields to optimization routine in SMCL. 
          ! philon(:) and plilat(:) are geographical coordinates of global field.20130812-7
          !    philon(64(=nlon)): 0..->>..+360[deg]
          !    philat(32(=ngl)) : +90..->>..-90[deg]
          ! v_z_g(:,:,:)        : global vertical wind velocity field [m/s] 
          CALL calculate_trajectory(ac_routes(:,:,i_traj_dep,1),    &   ! 20130904-3,20130909-2 
                                    p2_ac_routes_desc(i_traj_dep,:),& 
                                    philon(:),philat(:),zgl_geopot_3d(:,:,:),&   
                                    uwind_g(:,:,:),vwind_g(:,:,:),v_z_g(:,:,:), &
                                    t_scb_g(:,:,:),ac_fl_direction,cpc_g(:,:,:), &     !t_scb_g 20140522-2, ac_fl_direction20141217-1
                                    rho_air_dry_3d_g(:,:,:))              !20160301-1   
          ! The fuel/H2O_emis/Cost/COC optimal trajectory is found.                !20170222-1,20170224-1                        
       
       ELSEIF(option_traj_calc==NOx_opt)THEN     !20170215-2
          CALL calculate_trajectory(ac_routes(:,:,i_traj_dep,1),    &   ! 20130904-3,20130909-2 
                                    p2_ac_routes_desc(i_traj_dep,:),& 
                                    philon(:),philat(:),zgl_geopot_3d(:,:,:),&   
                                    uwind_g(:,:,:),vwind_g(:,:,:),v_z_g(:,:,:), &
                                    t_scb_g(:,:,:),ac_fl_direction,cpc_g(:,:,:), &     !t_scb_g 20140522-2, ac_fl_direction20141217-1
                                    rho_air_dry_3d_g(:,:,:), press_3d_g(:,:,:))     !20160301-1, press_3d_g 20170215-2   
          ! The NOx_emis optimal trajectory is found.                                        
       ELSEIF(option_traj_calc==ATR20_opt.or.option_traj_calc==COSTCLIM_opt)THEN     !20170801
!          write(*,*) "XXX1_PE=",p_pe,now,shape(cpc_g),"contrail_check_cpc_g_1=",cpc_g(:,:,:)
!          write(*,*) "XXXtscb_PE=",p_pe,now,shape(t_scb_g),"tscb_check_t_scb_g_1=",t_scb_g(:,:,:)
!          write(*,*) "XXXrho_PE=",p_pe,now,shape(rho_air_dry_3d_g),"rho_check_g_1=",rho_air_dry_3d_g(:,:,:)
!          write(*,*) "XXXuwind_PE=",p_pe,now,shape(uwind_g),"check_uwind_g_1=",uwind_g(:,:,:)
!          write(*,*) "XXXATR20_contrail_PE=",p_pe,now,shape(ATR20_contrail_g),"ATR20_contrail_1=",ATR20_contrail_g(:,:,:)
!          write(*,*) "XXXATR20_o3_PE=",p_pe,now,shape(ATR20_o3_g),"ATR20_o3_1=",ATR20_o3_g(:,:,:)
!          write(*,*) "XXXATR20_ch4_PE=",p_pe,now,shape(ATR20_ch4_g),"ATR20_ch4_1=",ATR20_ch4_g(:,:,:)
!          write(*,*) "XXXATR20_h2o_PE=",p_pe,now,shape(ATR20_h2o_g),"ATR20_h2o_1=",ATR20_h2o_g(:,:,:)
!          write(*,*) "XXXATR20_co2_PE=",p_pe,now,shape(ATR20_co2_g),"ATR20_co2_1=",ATR20_co2_g(:,:,:)
          CALL calculate_trajectory(ac_routes(:,:,i_traj_dep,1),    &    
                                    p2_ac_routes_desc(i_traj_dep,:),& 
                                    philon(:),philat(:),zgl_geopot_3d(:,:,:),&   
                                    uwind_g(:,:,:),vwind_g(:,:,:),v_z_g(:,:,:), &
                                    t_scb_g(:,:,:),ac_fl_direction,cpc_g(:,:,:), &   
                                    rho_air_dry_3d_g(:,:,:), press_3d_g(:,:,:),&
                                    ATR20_o3_g(:,:,:),ATR20_ch4_g(:,:,:),ATR20_h2o_g(:,:,:),&
                                    ATR20_contrail_g(:,:,:),ATR20_co2_g(:,:,:))     !20170801 atr20_g 
      
       ELSE
!          write(*,*)'ERROR:SMIL_call_calcuate_trajectory'
       ENDIF
!       write(*,*) "1_PE=",p_pe,"GC_ac_routes(lon)     =", ac_routes(:,1,i_traj_dep,1)
!       write(*,*) "1_PE=",p_pe,"GC_ac_routes(lat)     =", ac_routes(:,2,i_traj_dep,1)
!       write(*,*) "1_PE=",p_pe,"GC_ac_routes(alt)     =", ac_routes(:,3,i_traj_dep,1)
!       write(*,*) "1_PE=",p_pe,"GC_ac_routes(time)    =", ac_routes(:,4,i_traj_dep,1)
!       write(*,*) "1_PE=",p_pe,"GC_ac_routes(speed)   =", ac_routes(:,5,i_traj_dep,1)   !AC ground speed [km/h]
!       write(*,*) "1_PE=",p_pe,"GC_ac_routes(dist)    =", ac_routes(:,6,i_traj_dep,1)
!       write(*,*) "1_PE=",p_pe,"GC_ac_routes(fuel_use)=", ac_routes(:,7,i_traj_dep,1)
       !Yin_20170423
!       write(*,*) "1_PE=",p_pe,"GC_ac_routes(contrail_pc) =", ac_routes(:,10,i_traj_dep,1)
       !Yin_20170423
!       write(*,*) "XXX1_PE=",p_pe,"contrail_check_cpc_g_2=",cpc_g(:,:,:)
!       write(*,*) "XXXtscb__PE=",p_pe,"tscb_check_t_scb_g_1=",t_scb_g(:,:,:)
!          write(*,*) "BBB1_PE=",p_pe,now,shape(cpc_g),"contrail_check_cpc_g_1=",cpc_g(:,:,:)
!          write(*,*) "BBBtscb_PE=",p_pe,now,shape(t_scb_g),"tscb_check_t_scb_g_1=",t_scb_g(:,:,:)
!          write(*,*) "BBBrho_PE=",p_pe,now,shape(rho_air_dry_3d_g),"rho_check_g_1=",rho_air_dry_3d_g(:,:,:)
!          write(*,*) "BBBuwind_PE=",p_pe,now,shape(uwind_g),"check_uwind_g_1=",uwind_g(:,:,:)
!          write(*,*) "BBBATR20_contrail_PE=",p_pe,now,shape(ATR20_contrail_g),"ATR20_contrail_1=",ATR20_contrail_g(:,:,:)
!          write(*,*) "BBBATR20_o3_PE=",p_pe,now,shape(ATR20_o3_g),"ATR20_o3_1=",ATR20_o3_g(:,:,:)
!          write(*,*) "BBBATR20_ch4_PE=",p_pe,now,shape(ATR20_ch4_g),"ATR20_ch4_1=",ATR20_ch4_g(:,:,:)
!          write(*,*) "BBBATR20_h2o_PE=",p_pe,now,shape(ATR20_h2o_g),"ATR20_h2o_1=",ATR20_h2o_g(:,:,:)
!          write(*,*) "BBBATR20_co2_PE=",p_pe,now,shape(ATR20_co2_g),"ATR20_co2_1=",ATR20_co2_g(:,:,:)

       !GMT and GNUPLOT: global flight routes visualization.20140611-1,20141216-6,20150304-1,2
!       write(*,*)"ac_fl_direction_check:",ac_fl_direction
!       IF(p_pe.eq.0)THEN
!       IF((p_pe.ge.0).and.(p_pe.le.5))THEN
!          IF(ac_fl_direction.eq.0)THEN       !Eastbound
!             OPEN(50,file='fr_pe0_eb.dat',position='append')
!                do i_fl=1,nwaypoints     
!                   write(50,*)ac_routes(i_fl,1,i_traj_dep,1),ac_routes(i_fl,2,i_traj_dep,1),ac_routes(i_fl,3,i_traj_dep,1)
!                enddo
!                write(50,*)                  !for gnuplot
!                write(50,*)
!             CLOSE(50)
!             OPEN(52,file='fr_pe0_eb_gmt.dat',position='append')
!!                do i_fl=1,nwaypoints     
!!                   write(52,*)ac_routes(i_fl,1,i_traj_dep,1),ac_routes(i_fl,2,i_traj_dep,1),ac_routes(i_fl,3,i_traj_dep,1)
!!                enddo
!                write(52,'(a)')">"           !for GMT 
!             CLOSE(52)
!          ELSEIF(ac_fl_direction.eq.2)THEN   !Westbound
!             OPEN(51,file='fr_pe0_wb.dat',position='append')
!                do i_fl=1,nwaypoints     
!                   write(51,*)ac_routes(i_fl,1,i_traj_dep,1),ac_routes(i_fl,2,i_traj_dep,1),ac_routes(i_fl,3,i_traj_dep,1)
!                enddo
!                write(51,*)                  !for gnuplot
!                write(51,*)
!             CLOSE(51)
!             OPEN(53,file='fr_pe0_wb_gmt.dat',position='append')
!                do i_fl=1,nwaypoints     
!                   write(53,*)ac_routes(i_fl,1,i_traj_dep,1),ac_routes(i_fl,2,i_traj_dep,1),ac_routes(i_fl,3,i_traj_dep,1)
!                enddo
!                write(53,'(a)')">"           !for GMT
!             CLOSE(53)
!          ENDIF
!       ENDIF
!       IF(p_pe.eq.1)THEN
!       IF((p_pe.ge.6).and.(p_pe.le.11))THEN
!          IF(ac_fl_direction.eq.0)THEN       !Eastbound
!             OPEN(60,file='fr_pe1_eb.dat',position='append')
!                do i_fl=1,nwaypoints     
!                   write(60,*)ac_routes(i_fl,1,i_traj_dep,1),ac_routes(i_fl,2,i_traj_dep,1),ac_routes(i_fl,3,i_traj_dep,1)
!                enddo
!                write(60,*)
!                write(60,*)
!             CLOSE(60)
!             OPEN(62,file='fr_pe1_eb_gmt.dat',position='append')
!                do i_fl=1,nwaypoints     
!                   write(62,*)ac_routes(i_fl,1,i_traj_dep,1),ac_routes(i_fl,2,i_traj_dep,1),ac_routes(i_fl,3,i_traj_dep,1)
!                enddo
!                write(62,'(a)')">"
!             CLOSE(62)
!          ELSEIF(ac_fl_direction.eq.2)THEN   !Westbound
!             OPEN(61,file='fr_pe1_wb.dat',position='append')
!                do i_fl=1,nwaypoints     
!                   write(61,*)ac_routes(i_fl,1,i_traj_dep,1),ac_routes(i_fl,2,i_traj_dep,1),ac_routes(i_fl,3,i_traj_dep,1)
!                enddo
!                write(61,*)
!                write(61,*)
!             CLOSE(61)
!             OPEN(63,file='fr_pe1_wb_gmt.dat',position='append')
!                do i_fl=1,nwaypoints     
!                   write(63,*)ac_routes(i_fl,1,i_traj_dep,1),ac_routes(i_fl,2,i_traj_dep,1),ac_routes(i_fl,3,i_traj_dep,1)
!                enddo
!                write(63,'(a)')">"
!             CLOSE(63)
!          ENDIF
!       ENDIF
!       IF(p_pe.eq.2)THEN
!       IF((p_pe.ge.12).and.(p_pe.le.17))THEN
!          IF(ac_fl_direction.eq.0)THEN       !Eastbound
!             OPEN(70,file='fr_pe2_eb.dat',position='append')
!                do i_fl=1,nwaypoints     
!                   write(70,*)ac_routes(i_fl,1,i_traj_dep,1),ac_routes(i_fl,2,i_traj_dep,1),ac_routes(i_fl,3,i_traj_dep,1)
!                enddo
!                write(70,*)
!                write(70,*)
!             CLOSE(70)
!             OPEN(72,file='fr_pe2_eb_gmt.dat',position='append')
!                do i_fl=1,nwaypoints     
!                   write(72,*)ac_routes(i_fl,1,i_traj_dep,1),ac_routes(i_fl,2,i_traj_dep,1),ac_routes(i_fl,3,i_traj_dep,1)
!                enddo
!                write(72,'(a)')">"
!             CLOSE(72)
!          ELSEIF(ac_fl_direction.eq.2)THEN   !Westbound
!             OPEN(71,file='fr_pe2_wb.dat',position='append')
!                do i_fl=1,nwaypoints     
!                   write(71,*)ac_routes(i_fl,1,i_traj_dep,1),ac_routes(i_fl,2,i_traj_dep,1),ac_routes(i_fl,3,i_traj_dep,1)
!                enddo
!                write(71,*)
!                write(71,*)
!             CLOSE(71)
!             OPEN(73,file='fr_pe2_wb_gmt.dat',position='append')
!                do i_fl=1,nwaypoints     
!                   write(73,*)ac_routes(i_fl,1,i_traj_dep,1),ac_routes(i_fl,2,i_traj_dep,1),ac_routes(i_fl,3,i_traj_dep,1)
!                enddo
!                write(73,'(a)')">"
!             CLOSE(73)
!          ENDIF
!       ENDIF
!!       IF(p_pe.eq.3)THEN
!       IF((p_pe.ge.18).and.(p_pe.le.23))THEN
!          IF(ac_fl_direction.eq.0)THEN       !Eastbound
!             OPEN(80,file='fr_pe3_eb.dat',position='append')
!                do i_fl=1,nwaypoints     
!                   write(80,*)ac_routes(i_fl,1,i_traj_dep,1),ac_routes(i_fl,2,i_traj_dep,1),ac_routes(i_fl,3,i_traj_dep,1)
!                enddo
!                write(80,*)
!                write(80,*)
!             CLOSE(80)
!             OPEN(82,file='fr_pe3_eb_gmt.dat',position='append')
!                do i_fl=1,nwaypoints     
!                   write(82,*)ac_routes(i_fl,1,i_traj_dep,1),ac_routes(i_fl,2,i_traj_dep,1),ac_routes(i_fl,3,i_traj_dep,1)
!                enddo
!                write(82,'(a)')">"
!             CLOSE(82)
!          ELSEIF(ac_fl_direction.eq.2)THEN   !Westbound
!             OPEN(81,file='fr_pe3_wb.dat',position='append')
!                do i_fl=1,nwaypoints     
!                   write(81,*)ac_routes(i_fl,1,i_traj_dep,1),ac_routes(i_fl,2,i_traj_dep,1),ac_routes(i_fl,3,i_traj_dep,1)
!                enddo
!                write(81,*)
!                write(81,*)
!             CLOSE(81)
!             OPEN(83,file='fr_pe3_wb_gmt.dat',position='append')
!                do i_fl=1,nwaypoints     
!                   write(83,*)ac_routes(i_fl,1,i_traj_dep,1),ac_routes(i_fl,2,i_traj_dep,1),ac_routes(i_fl,3,i_traj_dep,1)
!                enddo
!                write(83,'(a)')">"
!             CLOSE(83)
!          ENDIF
!       ENDIF

!see 20130125-1       
!       if((-150.0_dp<=ac_routes(1,1,i_traj_dep,1) .and. ac_routes(1,1,i_traj_dep,1)<= -60.0_dp) .and.&
!          (0.0_dp   <=ac_routes(nwaypoints,1,i_traj_dep,1) .and.  ac_routes(nwaypoints,1,i_traj_dep,1)<= 30.0_dp))then
!          write(*,*)'FROM_USA_TO_FRA',i_traj_dep,'PE=',p_pe 
!       endif
       !see20131007-2
!       if((-75.0_dp<=ac_routes(nwaypoints,1,i_traj_dep,1) .and.  ac_routes(nwaypoints,1,i_traj_dep,1)<= -71.0_dp) .and.&
!          (7.0_dp   <=ac_routes(1,1,i_traj_dep,1) .and. ac_routes(1,1,i_traj_dep,1)<= 9.0_dp))then
!          write(*,*)'FROM_FRA_TO_JFK',i_traj_dep,'PE=',p_pe 
!       endif
!       if((-75.0_dp<=ac_routes(nwaypoints,1,i_traj_dep,1) .and.  ac_routes(nwaypoints,1,i_traj_dep,1)<= -71.0_dp) .and.&
!          (10.0_dp   <=ac_routes(1,1,i_traj_dep,1) .and. ac_routes(1,1,i_traj_dep,1)<= 12.0_dp))then
!          write(*,*)'FROM_MUC_TO_JFK',i_traj_dep,'PE=',p_pe 
!       endif
!see20140224-4
!       if((-88.0_dp<=ac_routes(nwaypoints,1,i_traj_dep,1) .and.  ac_routes(nwaypoints,1,i_traj_dep,1)<= -86.0_dp) .and.&
!          (7.0_dp   <=ac_routes(1,1,i_traj_dep,1) .and. ac_routes(1,1,i_traj_dep,1)<= 9.0_dp))then
!          write(*,*)'FROM_FRA_TO_ORD',i_traj_dep,ac_routes(1,1,i_traj_dep,1),ac_routes(1,2,i_traj_dep,1),&
!                                      ac_routes(nwaypoints,1,i_traj_dep,1),ac_routes(nwaypoints,2,i_traj_dep,1) 
!       endif
!       if((-78.0_dp<=ac_routes(nwaypoints,1,i_traj_dep,1) .and.  ac_routes(nwaypoints,1,i_traj_dep,1)<= -76.0_dp) .and.&
!          (7.0_dp   <=ac_routes(1,1,i_traj_dep,1) .and. ac_routes(1,1,i_traj_dep,1)<= 9.0_dp))then
!          write(*,*)'FROM_FRA_TO_IAD',i_traj_dep,ac_routes(1,1,i_traj_dep,1),ac_routes(1,2,i_traj_dep,1),&
!                                      ac_routes(nwaypoints,1,i_traj_dep,1),ac_routes(nwaypoints,2,i_traj_dep,1)
!       endif

       !Pick out pressure, temperature, rho at waypoints using the indices (d_jgx,d_idx,d_jgy) from global fields.
       !Global field index    : d_jgx, d_idx, d_jgy
       !Decomposed field index: d_jp, (d_idx: common to global firld), d_jrow, d_pe(local PE)
       !By using locate_in_decomp and nn_index, information on global field index, local PE(d_pe), decomposed field index
       !are obtained. Thus, they are utilized if necessary. 
       !Loop over nwaypoints-1
       pass_phyval_along_traj = 0.0_dp
       DO i_wp=1,nwaypoints-1
          !HORIZONTAL INDEX 
          !write(*,*)'(press)now=',now
          !write(*,*)'(press)check_before_locate','p_pe=',p_pe,'p_io=',p_io
          !write(*,*)'(press)i_wp =',i_wp,'lon=',ac_routes(i_wp,1,i_traj_dep,1),'lat=',ac_routes(i_wp,2,i_traj_dep,1)
          CALL locate_in_decomp(status, ac_routes(i_wp,1,i_traj_dep,1), ac_routes(i_wp,2,i_traj_dep,1), &
                                d_pe, d_jp, d_jrow, d_jgx, d_jgy) 
          !VERTICAL INDEX 
          zalt(:)=zgl_geopot_3d(d_jgx,:,d_jgy)/g   !zalt:geopotential altitude,20131105-1
          !write(*,*)'locate_in_check', ' lon=',ac_routes(i_wp,1,i_traj_dep,1),'index=',d_jgx   !20140606-1
          !write(*,*)'(press)g=',g,'(press)nlev=',nlev 
          !write(*,*)'(press)zgl_geopot_3d=',zgl_geopot_3d(d_jgx,:,d_jgy)/g 
          !write(*,*)'(press)zalt=',zalt(:)
          CALL nn_index(zalt(:), ac_routes(i_wp,3,i_traj_dep,1), d_idx)
          pass_phyval_along_traj(i_wp,1)=press_3d_g(d_jgx,d_idx,d_jgy) 
          pass_phyval_along_traj(i_wp,2)=t_scb_g(d_jgx,d_idx,d_jgy) 
          pass_phyval_along_traj(i_wp,3)=rho_air_dry_3d_g(d_jgx,d_idx,d_jgy) 

!         write(*,*)'(press)alt=',ac_routes(i_wp,3,i_traj_dep,1),'d_idx=',d_idx
!         write(*,*)'(press)d_jgx, d_idx, d_jgy, g =',d_jgx,d_idx,d_jgy,g
!         write(*,*)'(press)press_3d_g =',press_3d_g(d_jgx,d_idx,d_jgy)
!         write(*,*)'(press)pass_phyval_along_traj=',pass_phyval_along_traj(i_wp,1)
!         write(*,*)'(temp)t_scb_g =',t_scb_g(d_jgx,d_idx,d_jgy)
!         write(*,*)'(temp)pass_phyval_along_traj=',pass_phyval_along_traj(i_wp,2)
!         write(*,*)'(rho)rho_air_dra_3d_g =',rho_air_dry_3d_g(d_jgx,d_idx,d_jgy)
!         write(*,*)'(rho)pass_phyval_along_traj=',pass_phyval_along_traj(i_wp,3)
!         write(*,*)'---------------------------------------------------------------'
       ENDDO

    ! ----------------------------------------------------------------------
    ! This subroutine calculates fuel flow and emissions along the projected
    ! flight routes. Pass physical values along waypoints to SMCL:
    ! 
    ! Note: - pass_phyval_along_traj(nwaypoints-1,1): pressure, [Pa]
    !       - pass_phyval_along_traj(nwaypoints-1,2): temperature, [K]
    !       - pass_phyval_along_traj(nwaypoints-1,3): rho, [kg/m^3]
    ! ----------------------------------------------------------------------
!HY20120724-6,20120404-1
!      call calculate_emissions_along_trajectory(ac_routes(:,:,i_traj,1),   &    ! One Route
!                                 p2_ac_routes_desc(i_traj,:))    ! Route description,Is the ac_routes maybe 'p3_ac_routes'?
!HY20120820-1
!      call calculate_emissions_along_trajectory(ac_routes(:,:,i_traj_dep,1),   &    ! One Route
!                                 p2_ac_routes_desc(i_traj_dep,:))    ! Route description,Is the ac_routes maybe 'p3_ac_routes'?
!       CALL calculate_emissions_along_trajectory(ac_routes(:,:,i_traj_dep,1), pass_phyval_along_traj(:,:))
       CALL calculate_emissions_along_trajectory(ac_routes(:,:,i_traj_dep,1), pass_phyval_along_traj(:,:),i_traj_dep,p_pe,&
                                                 philon(:),philat(:),zgl_geopot_3d(:,:,:),&
                                                 ATR20_o3_g,ATR20_ch4_g,ATR20_h2o_g,ATR20_contrail_g,ATR20_co2_g)
!       write(*,*) "2_PE=",p_pe,"GC_ac_routes(fuel_use2)=", ac_routes(:,7,i_traj_dep,1)
!       write(*,*) "2_PE=",p_pe,"GC_ac_routes(emis_nox) =", ac_routes(:,8,i_traj_dep,1)
!       write(*,*) "2_PE=",p_pe,"GC_ac_routes(emis_h2o) =", ac_routes(:,9,i_traj_dep,1)
!       write(*,*) "2_PE=",p_pe,"GC_ac_routes(potcov) =", ac_routes(:,10,i_traj_dep,1)
   !Yin_20170801+
!       write(*,*) "1_PE=",p_pe,"GC_ac_routes(atr20o3_pc) =", ac_routes(:,11,i_traj_dep,1)
!       write(*,*) "1_PE=",p_pe,"GC_ac_routes(atr20ch4_pc) =", ac_routes(:,12,i_traj_dep,1)
!       write(*,*) "1_PE=",p_pe,"GC_ac_routes(atr20h2o_pc) =", ac_routes(:,13,i_traj_dep,1)
!       write(*,*) "1_PE=",p_pe,"GC_ac_routes(atr20cpc_pc) =", ac_routes(:,14,i_traj_dep,1)
!       write(*,*) "1_PE=",p_pe,"GC_ac_routes(atr20co2_pc) =", ac_routes(:,15,i_traj_dep,1)
!       write(*,*) "1_PE=",p_pe,"GC_ac_routes(atr20tot_pc) =", ac_routes(:,16,i_traj_dep,1)
!          write(*,*) "CCC1_PE=",p_pe,now,shape(cpc_g),"contrail_check_cpc_g_1=",cpc_g(:,:,:)
!          write(*,*) "CCCtscb_PE=",p_pe,now,shape(t_scb_g),"tscb_check_t_scb_g_1=",t_scb_g(:,:,:)
!          write(*,*) "CCCrho_PE=",p_pe,now,shape(rho_air_dry_3d_g),"rho_check_g_1=",rho_air_dry_3d_g(:,:,:)
!          write(*,*) "CCCuwind_PE=",p_pe,now,shape(uwind_g),"check_uwind_g_1=",uwind_g(:,:,:)
!          write(*,*) "CCCATR20_contrail_PE=",p_pe,now,shape(ATR20_contrail_g),"ATR20_contrail_1=",ATR20_contrail_g(:,:,:)
!          write(*,*) "CCCATR20_o3_PE=",p_pe,now,shape(ATR20_o3_g),"ATR20_o3_1=",ATR20_o3_g(:,:,:)
!          write(*,*) "CCCATR20_ch4_PE=",p_pe,now,shape(ATR20_ch4_g),"ATR20_ch4_1=",ATR20_ch4_g(:,:,:)
!          write(*,*) "CCCATR20_h2o_PE=",p_pe,now,shape(ATR20_h2o_g),"ATR20_h2o_1=",ATR20_h2o_g(:,:,:)
!          write(*,*) "CCCATR20_co2_PE=",p_pe,now,shape(ATR20_co2_g),"ATR20_co2_1=",ATR20_co2_g(:,:,:)
   !Yin_20170801-

    !Visualization for 4th EMAC, 20140703-1,20141218-4,20150306-1
    !MUC >> JFK
!    IF(p_pe.eq.0 .and. i_traj_dep.eq.20)THEN 
!       OPEN(61,file='vis_windfield_muc_jfk.dat')
!          do i=1,128   !64 longitude
!             do k=1,31   !19 altitude
!                do j=1,64   !32 latitude    
!                   write(61,*)philon(i),philat(j),zalt(k),uwind_g(i,k,j),vwind_g(i,k,j),v_z_g(i,k,j)
!                enddo
!             enddo
!          enddo
!       CLOSE(61)
!    ENDIF
!    !JFK >> MUC  20141218-4
!    IF(p_pe.eq.3 .and. i_traj_dep.eq.1)THEN 
!       OPEN(62,file='vis_windfield_jfk_muc.dat')
!          do i=1,128   !64 longitude
!             do k=1,31   !19 altitude
!                do j=1,64   !32 latitude    
!                   write(62,*)philon(i),philat(j),zalt(k),uwind_g(i,k,j),vwind_g(i,k,j),v_z_g(i,k,j)
!                enddo
!             enddo
!          enddo
!       CLOSE(62)
!    ENDIF
!    !EHAM >> KMSP
!    IF(p_pe.eq.2 .and. i_traj_dep.eq.14)THEN 
!       OPEN(63,file='vis_windfield_eham_kmsp.dat')
!          do i=1,128   !64 longitude
!             do k=1,31   !19 altitude
!                do j=1,64   !32 latitude    
!                   write(63,*)philon(i),philat(j),zalt(k),uwind_g(i,k,j),vwind_g(i,k,j),v_z_g(i,k,j)
!                 !  write(63,*)philon(i),philat(j),zalt(k),uwind_g(i,k,j),vwind_g(i,k,j),cpc_g(i,k,j)   !20170508-1
!                enddo
!             enddo
!          enddo
!       CLOSE(63)
!    ENDIF
!    !KMSP >> EHAM
!    IF(p_pe.eq.1 .and. i_traj_dep.eq.24)THEN 
!       OPEN(64,file='vis_windfield_kmsp_eham.dat')
!          do i=1,128   !64 longitude
!             do k=1,31   !19 altitude
!                do j=1,64   !32 latitude    
!                   write(64,*)philon(i),philat(j),zalt(k),uwind_g(i,k,j),vwind_g(i,k,j),v_z_g(i,k,j)
!                 !  write(64,*)philon(i),philat(j),zalt(k),uwind_g(i,k,j),vwind_g(i,k,j),cpc_g(i,k,j)   !20170508-1
!                enddo
!             enddo
!          enddo
!       CLOSE(64)
!    ENDIF
!    !EHAM >> KSEA
!    IF(p_pe.eq.1 .and. i_traj_dep.eq.17)THEN 
!       OPEN(65,file='vis_windfield_eham_ksea.dat')
!          do i=1,128   !64 longitude
!             do k=1,31   !19 altitude
!                do j=1,64   !32 latitude    
!                   write(65,*)philon(i),philat(j),zalt(k),uwind_g(i,k,j),vwind_g(i,k,j),v_z_g(i,k,j)
!                enddo
!             enddo
!          enddo
!       CLOSE(65)
!    ENDIF
!    !KSEA >> EHAM
!    IF(p_pe.eq.1 .and. i_traj_dep.eq.21)THEN 
!       OPEN(66,file='vis_windfield_ksea_eham.dat')
!          do i=1,128   !64 longitude
!             do k=1,31   !19 altitude
!                do j=1,64   !32 latitude    
!                   write(66,*)philon(i),philat(j),zalt(k),uwind_g(i,k,j),vwind_g(i,k,j),v_z_g(i,k,j)
!                enddo
!             enddo
!          enddo
!       CLOSE(66)
!    ENDIF


    ENDDO
    DEALLOCATE(zalt)                       !Array 
    DEALLOCATE(pass_phyval_along_traj)     !Array
    
!   if (lupdate_traj) call update_trajectory

    DO i_traj=1,nlroutes    !20120719
         if (p2_ac_routes_desc(i_traj,AC_status)==inflight)then      !&    !inflight=1.0
!            p2_ac_routes_desc(i_traj,AC_POS_OLD)=p2_ac_routes_desc(i_traj,AC_POS) !AC_pos -> AC_pos_old            
             !write(*,*)'d_now',now,'ac_routes(time)',ac_routes(:,4,i_traj,1) 
             !write(*,*)'aircraft_number=',i_traj 
             !write(*,*)'inflight_check',p2_ac_routes_desc(i_traj,AC_status) 
             !write(*,*)'p2_ac_routes_desc(AC_POS_OLD,AC_POS)1',p2_ac_routes_desc(i_traj,AC_POS_OLD),p2_ac_routes_desc(i_traj,AC_POS) 
             call fly_aircraft(now, ac_routes(:,:,i_traj,1), p2_ac_routes_desc(i_traj,:))
             !write(*,*)'p2_ac_routes_desc(AC_POS_OLD,AC_POS)2',p2_ac_routes_desc(i_traj,AC_POS_OLD),p2_ac_routes_desc(i_traj,AC_POS) 
         
             !GMT and GNUPLOT: global flight routes visualization for time series.20140610-2,4
             !This structure is the same of subroutine nn_index of messy_main_tools.f90
             !see line 1955 
             !idx_vis  = 0
             dmin_vis = ABS(ac_routes(1,4,i_traj,1)-ac_routes(nwaypoints,4,i_traj,1))
             !write(*,*)'dmin_vis_init=',dmin_vis
             DO j_vis=1,nwaypoints
                IF(ABS(ac_routes(j_vis,4,i_traj,1)-now) <= dmin_vis)THEN  
                   dmin_vis = ABS(ac_routes(j_vis,4,i_traj,1)-now)
                   !idx_vis  = j_vis
                   !write(*,*)'dmin_vis=',dmin_vis,'idx_vis',idx_vis 
                ENDIF
             ENDDO       
             !write(*,*)'now=',now,ac_routes(idx_vis,1,i_traj,1),ac_routes(idx_vis,2,i_traj,1),ac_routes(idx_vis,3,i_traj,1)
!             IF(p_pe.eq.0)THEN
!             IF((p_pe.ge.0).and.(p_pe.le.5))THEN
!                OPEN(10,file='gnu001.dat',position='append')
!                   write(10,*)ac_routes(idx_vis,1,i_traj,1),ac_routes(idx_vis,2,i_traj,1),ac_routes(idx_vis,3,i_traj,1)
                !  write(10,*)
                !  write(10,'(a)')">"
!                CLOSE(10)
!             ENDIF
!             IF(p_pe.eq.1)THEN
!             IF((p_pe.ge.6).and.(p_pe.le.11))THEN
!                OPEN(20,file='gnu002.dat',position='append')
!                   write(20,*)ac_routes(idx_vis,1,i_traj,1),ac_routes(idx_vis,2,i_traj,1),ac_routes(idx_vis,3,i_traj,1)
                !  write(20,*)
                !  write(20,'(a)')">"
!                CLOSE(20)
!             ENDIF
!             IF(p_pe.eq.2)THEN
!             IF((p_pe.ge.12).and.(p_pe.le.17))THEN
!                OPEN(30,file='gnu003.dat',position='append')
!                   write(30,*)ac_routes(idx_vis,1,i_traj,1),ac_routes(idx_vis,2,i_traj,1),ac_routes(idx_vis,3,i_traj,1)
                !  write(30,*)
                !  write(30,'(a)')">"
!                CLOSE(30)
!             ENDIF
!             IF(p_pe.eq.3)THEN
!             IF((p_pe.ge.18).and.(p_pe.le.23))THEN
!                OPEN(40,file='gnu004.dat',position='append')
!                   write(40,*)ac_routes(idx_vis,1,i_traj,1),ac_routes(idx_vis,2,i_traj,1),ac_routes(idx_vis,3,i_traj,1)
                !  write(40,*)
                !  write(40,'(a)')">"
!                CLOSE(40)
!             ENDIF
         endif

! After emission calculation at last A_city waypoint, Ac_pos and AC_pos_old should be changed into initial value.
! That is, Emission values at final waypoint(A_city) should be calculated.
!HY      if (p2_ac_routes_desc(i_traj,AC_pos)==1._dp) then     
           ! Now 1.0 is added soon by SMCL as dummy subroutine. Note 20120724-5. 
         ! The dummy 1.0_dp should be changed. See 20120903-2,20120827-2  
         if (p2_ac_routes_desc(i_traj,AC_pos)==REAL(nwaypoints,DP)) then 
             p2_ac_routes_desc(i_traj,AC_status)=arrived        ! Arrived=2.0 
             ! Counting arrivals 
             p2_ac_routes_desc(i_traj,total_arrival)=p2_ac_routes_desc(i_traj,total_arrival) + 1.0_dp  
!             write(*,*)'aircraft_arrived!'
!             write(*,*)'p2_ac_routes_desc(i_traj,AC_pos)=',p2_ac_routes_desc(i_traj,AC_pos)
!             write(*,*)'p2_ac_routes_desc(i_traj,AC_pos_old)=',p2_ac_routes_desc(i_traj,AC_pos_old)
!             write(*,*)'p2_ac_routes_desc(i_traj,AC_status)=',p2_ac_routes_desc(i_traj,AC_status)
             IF(p2_ac_routes_desc(i_traj,total_arrival)== 1.0_dp)THEN
                i_traj_out = i_traj_out + 1
                p2_ac_routes_desc(i_traj,corr_i_traj_out)=REAL(i_traj_out,DP)
                call compress_traj(ac_routes(:,:,i_traj,1), ac_routes_out(:,:,i_traj_out,1), Julian_date_start)
!                write(*,*)'i_traj_out=',i_traj_out,"PE=",p_pe,'now=',now
             ELSEIF(p2_ac_routes_desc(i_traj,total_arrival)>= 2.0_dp)THEN
               !i_traj_out_reuse = IDINT(p2_ac_routes_desc(i_traj,corr_i_traj_out))   !20170212-2 
                i_traj_out_reuse = INT(p2_ac_routes_desc(i_traj,corr_i_traj_out)) 
                call compress_traj(ac_routes(:,:,i_traj,1), ac_routes_out(:,:,i_traj_out_reuse,1), Julian_date_start)
!                write(*,*)'i_traj_out_reuse=',i_traj_out_reuse,"PE=",p_pe,'now=',now
             ENDIF
             if (i_traj_out > nlroutes_out) then 
                 ! find a solution according to parameters
                 ! - just ignore output or
                 ! - stop simulation
             endif
!            call compress_traj(ac_routes(:,:,i_traj,1), ac_routes_out(:,:,i_traj_out,1), Julian_date_start)
             ! See 20120906-1,-4,20120905-3,20120404-2,20120724-3
!             write(*,*) "3_PE=",p_pe,"ac_routes(lon)=",ac_routes(:,1,i_traj,1)
!             write(*,*) "3_PE=",p_pe,"ac_routes_out(lon)=",ac_routes_out(:,1,i_traj_out,1)
!             write(*,*) "4_PE=",p_pe,"ac_routes(time)=",ac_routes(:,4,i_traj,1)
!             write(*,*) "4_PE=",p_pe,"ac_routes_out(time)=",ac_routes_out(:,4,i_traj_out,1)
         endif
    ENDDO
!    write(*,*)"Number of arrived aircrafts(PE=",p_pe,"):",i_traj_out

! Collect emissions one after another, i.e. do not allocate global fields in parallel decomposition.
! Fill emission field for all 'inflight' and 'arrived' trajectories     
! Calculate sum over all PEs and decompose it
! When you create 'call fill_emission_field', check notes 20120329-3, 20120410-1, 20120713-3,-4, 20111215-8
! See 20120404-2, should we calculate emissions along traj here? or within call fly_aircraft? 
!-------------------------------------------
!   NOx
!-------------------------------------------
    ALLOCATE(glob_NOx_emis(nlon, nlev, ngl))    !Pointer: 
    !Yin_20170321+
    ALLOCATE(glob_teNOx_emis(nlon, nlev, ngl))    !Pointer:
    !Yin_20170321-
    ALLOCATE(add_glob_NOx_emis(nwaypoints-1,4)) !Array: 1:lon, 2:lat, 3:alt, 4:Emissions 
    ALLOCATE(zalt(nlev))                        !Array:
    glob_NOx_emis = 0.0_dp   
    glob_teNOx_emis = 0.0_dp                    !HY_20180524
    !write(*,*)"NOx_check_used_PEs_1=",p_pe,"p_io=",p_io
    do i_traj=1,nlroutes 
    !write(*,*)"NOx_check_used_PEs_2=",p_pe,"p_io=",p_io
        if (p2_ac_routes_desc(i_traj,AC_status)==arrived .or.     &
            p2_ac_routes_desc(i_traj,AC_status)==inflight)then   !&
            call fill_emission_field(NOx,ac_routes(:,:,i_traj,1),                    &
                                     p2_ac_routes_desc(i_traj,AC_pos),               &
                                     p2_ac_routes_desc(i_traj,AC_pos_old),add_glob_NOx_emis(:,:),i_add_NOx_emis)
!                                    p2_ac_routes_desc(i_traj,AC_pos_old),glob_NOx_emis)
!           write(*,*)'i_add_NOx_emis(SMIL)=',i_add_NOx_emis,'p_pe=',p_pe,'p_io=',p_io
!!         my_press = 500
           IF(i_add_NOx_emis/=0)THEN 
              DO i_add=1,i_add_NOx_emis
                 !HORIZONTAL INDEX 
!                 write(*,*)'(NOx)check_before_locate','p_pe=',p_pe,'p_io=',p_io
!                 write(*,*)'add_glob_NOx_emis(SMIL)1 =',add_glob_NOx_emis(i_add,:)
                 ! See 20121127-1,20121117-2
                 CALL locate_in_decomp(status, add_glob_NOx_emis(i_add,1), add_glob_NOx_emis(i_add,2), &
                                       d_pe, d_jp, d_jrow, d_jgx, d_jgy)   !20140606-1
                 ! If use pressure, see 20121108-10
!!               write(*,*)'press_3d=',press_3d(d_jp,:,d_jrow)
!!               IF(d_pe==p_pe)THEN

                 !VERTICAL INDEX 
                 zalt(:)=zgl_geopot_3d(d_jgx,:,d_jgy)/g   !zalt:geopotential altitude,20131105-1
!                 write(*,*)'(NOx)zalt=',zalt(:)
                 CALL nn_index(zalt(:), add_glob_NOx_emis(i_add,3), d_idx)
!!               CALL nn_index(press_3d(d_jp,:,d_jrow), my_press, d_idx)
                 glob_NOx_emis(d_jgx,d_idx,d_jgy)=add_glob_NOx_emis(i_add,4)
                 !Yin_20170321+
                 glob_teNOx_emis(d_jgx,d_idx,d_jgy)=glob_NOx_emis(d_jgx,d_idx,d_jgy)/(MN+MO)/&
                                                    (rho_air_dry_3d_g(d_jgx,d_idx,d_jgy)*grvol_g(d_jgx,d_idx,d_jgy)*1000._dp/M_air)
!HY20180526-0                                       (rho_air_dry_3d(d_jgx,d_idx,d_jgy)*grvol(d_jgx,d_idx,d_jgy)*1000._dp/M_air)
                 !Yin_20170321- 
!                 write(*,*)'(NOx)d_jgx, d_idx, d_jgy, g =',d_jgx,d_idx,d_jgy,g
!                 write(*,*)'(NOx)i_add_NOx_emis(SMIL)=',i_add_NOx_emis
!                 write(*,*)'(NOx)add_glob_NOx_emis(SMIL)2 =',add_glob_NOx_emis(i_add,:)
!                 write(*,*)'(NOx)glob_NOx_emis(SMIL) =',glob_NOx_emis(d_jgx,d_idx,d_jgy)
!!               ENDIF
              ENDDO
           ENDIF
        endif
    enddo 
!!  write(*,*)'press_3d=',press_3d(:,:,:)
!   write(*,*)'before_trp_N: NOx_emis=',NOx_emis(:,:,:),'glob_NOx_emis=',glob_NOx_emis(:,:,:),'PE=',p_pe,'M_SUM=',M_SUM
!    write(*,*)'(NOx)nlon,nlev,ngl=',nlon, nlev, ngl
!    write(*,*)"(NOx)check_used_PEs3=",p_pe,"p_io=",p_io
    ! See 20121102-2,20120329-2,20120410-2
    ! CALL trp_gpdc_gpgl(sign, lf, gf, method)
    !    sign = -1        : global field -> decomposed field 
    !    lf               : (local)decomposed field = NOx_emis
    !    gf               : global field            = glob_NOx_emis
    !    method(optional) : M_SUM=0 which means 'Sum over all PEs'.
    !                       M_SUM comes from 'messy_main_transform_bi.f90'.
    ! Decomposed fiels 'NOx_emis' can be used for other calculations after CALL trp_gpdc_gpgl.
    ! E.g., the field is converted into other climate impact metrics.
    CALL trp_gpdc_gpgl(-1, NOx_emis, glob_NOx_emis, M_SUM) 
    !Yin_20170321+
    CALL trp_gpdc_gpgl(-1, teNOx, glob_teNOx_emis, M_SUM)
    !Yin_20170321-

!   write(*,*)'after_trp_N: NOx_emis=',NOx_emis(:,:,:),'glob_NOx_emis=',glob_NOx_emis(:,:,:),'PE=',p_pe,'M_SUM=',M_SUM
    DEALLOCATE(glob_NOx_emis); NULLIFY(glob_NOx_emis)   !Pointer 
    DEALLOCATE(glob_teNOx_emis); NULLIFY(glob_teNOx_emis)   !Pointer
    DEALLOCATE(add_glob_NOx_emis)                       !Array
    DEALLOCATE(zalt)                                    !Array

!-------------------------------------------
!   H2O
!-------------------------------------------
    ALLOCATE(glob_H2O_emis(nlon, nlev, ngl))  
    ALLOCATE(add_glob_H2O_emis(nwaypoints-1,4)) !1:lon, 2:lat, 3:alt, 4:Emissions
    ALLOCATE(zalt(nlev))
    glob_H2O_emis = 0.0_dp   
    !write(*,*)"H2O_check_used_PEs_1=",p_pe,"p_io=",p_io
    do i_traj=1,nlroutes 
    !write(*,*)"H2O_check_used_PEs_2=",p_pe,"p_io=",p_io
        if (p2_ac_routes_desc(i_traj,AC_status)==arrived .or.     &
            p2_ac_routes_desc(i_traj,AC_status)==inflight)then   !&
            call fill_emission_field(H2O,ac_routes(:,:,i_traj,1),                    &
                                     p2_ac_routes_desc(i_traj,AC_pos),               &
                                     p2_ac_routes_desc(i_traj,AC_pos_old),add_glob_H2O_emis(:,:),i_add_H2O_emis)
!                                    p2_ac_routes_desc(i_traj,AC_pos_old),glob_H2O_emis)
!           write(*,*)'i_add_H2O_emis(SMIL)=',i_add_H2O_emis,'p_pe=',p_pe,'p_io=',p_io
           IF(i_add_H2O_emis/=0)THEN 
              DO i_add=1,i_add_H2O_emis
                 !HORIZONTAL INDEX 
!                 write(*,*)'(H2O)check_before_locate','p_pe=',p_pe,'p_io=',p_io
!                 write(*,*)'add_glob_H2O_emis(SMIL)1 =',add_glob_H2O_emis(i_add,:)
                 CALL locate_in_decomp(status, add_glob_H2O_emis(i_add,1), add_glob_H2O_emis(i_add,2), &
                                       d_pe, d_jp, d_jrow, d_jgx, d_jgy)   !20140606-1
                 !VERTICAL INDEX 
                 zalt(:)=zgl_geopot_3d(d_jgx,:,d_jgy)/g   !zalt:geopotential altitude,20131105-1
!                 write(*,*)'(H2O)zalt=',zalt(:)
                 CALL nn_index(zalt(:), add_glob_H2O_emis(i_add,3), d_idx)
                 glob_H2O_emis(d_jgx,d_idx,d_jgy)=add_glob_H2O_emis(i_add,4) 
!                 write(*,*)'(H2O)d_jgx, d_idx, d_jgy, g =',d_jgx,d_idx,d_jgy,g
!                 write(*,*)'(H2O)i_add_H2O_emis(SMIL)=',i_add_H2O_emis
!                 write(*,*)'(H2O)add_glob_H2O_emis(SMIL)2 =',add_glob_H2O_emis(i_add,:)
!                 write(*,*)'(H2O)glob_H2O_emis(SMIL) =',glob_H2O_emis(d_jgx,d_idx,d_jgy)
              ENDDO
           ENDIF
        endif
    enddo 
!    write(*,*)'(H2O)nlon,nlev,ngl=',nlon, nlev, ngl
!    write(*,*)"(H2O)check_used_PEs3=",p_pe,"p_io=",p_io
    CALL trp_gpdc_gpgl(-1, H2O_emis, glob_H2O_emis, M_SUM)
    DEALLOCATE(glob_H2O_emis); NULLIFY(glob_H2O_emis)
    DEALLOCATE(add_glob_H2O_emis)
    DEALLOCATE(zalt) 

!-------------------------------------------
!   DISTANCE
!-------------------------------------------
    ALLOCATE(glob_distance(nlon, nlev, ngl))  
    ALLOCATE(add_glob_distance(nwaypoints-1,4)) !1:lon, 2:lat, 3:alt, 4:Emissions
    ALLOCATE(zalt(nlev))
    glob_distance = 0.0_dp   
    !write(*,*)"DIST_check_used_PEs_1=",p_pe,"p_io=",p_io
    do i_traj=1,nlroutes 
    !write(*,*)"DIST_check_used_PEs_2=",p_pe,"p_io=",p_io
        if (p2_ac_routes_desc(i_traj,AC_status)==arrived .or.     &
            p2_ac_routes_desc(i_traj,AC_status)==inflight)then   !&
            call fill_emission_field(DIST,ac_routes(:,:,i_traj,1),                    &
                                     p2_ac_routes_desc(i_traj,AC_pos),               &
                                     p2_ac_routes_desc(i_traj,AC_pos_old),add_glob_distance(:,:),i_add_distance)
!                                    p2_ac_routes_desc(i_traj,AC_pos_old),glob_distance)
!           write(*,*)'i_add_distance(SMIL)=',i_add_distance,'p_pe=',p_pe,'p_io=',p_io
           IF(i_add_distance/=0)THEN 
              DO i_add=1,i_add_distance
                 !HORIZONTAL INDEX 
!                 write(*,*)'(DIST)check_before_locate','p_pe=',p_pe,'p_io=',p_io
!                 write(*,*)'add_glob_distance(SMIL)1 =',add_glob_distance(i_add,:)
                 CALL locate_in_decomp(status, add_glob_distance(i_add,1), add_glob_distance(i_add,2), &
                                       d_pe, d_jp, d_jrow, d_jgx, d_jgy)   !20140606-1
                 !VERTICAL INDEX 
                 zalt(:)=zgl_geopot_3d(d_jgx,:,d_jgy)/g  !zalt:geopotential altitutde,20131105-1
!                 write(*,*)'(DIST)zalt=',zalt(:)
                 CALL nn_index(zalt(:), add_glob_distance(i_add,3), d_idx)
                 glob_distance(d_jgx,d_idx,d_jgy)=add_glob_distance(i_add,4) 
!                 write(*,*)'(DIST)d_jgx, d_idx, d_jgy, g =',d_jgx,d_idx,d_jgy,g
!                 write(*,*)'(DIST)i_add_distance(SMIL)=',i_add_distance
!                 write(*,*)'(DIST)add_glob_distance(SMIL)2 =',add_glob_distance(i_add,:)
!                 write(*,*)'(DIST)glob_distance(SMIL) =',glob_distance(d_jgx,d_idx,d_jgy)
              ENDDO
           ENDIF
        endif
    enddo     
!    write(*,*)'(DIST)nlon,nlev,ngl=',nlon, nlev, ngl
!    write(*,*)"(DIST)check_used_PEs3=",p_pe,"p_io=",p_io
    CALL trp_gpdc_gpgl(-1, distance, glob_distance, M_SUM)
    DEALLOCATE(glob_distance); NULLIFY(glob_distance)
    DEALLOCATE(add_glob_distance)
    DEALLOCATE(zalt) 

!-------------------------------------------
!   FUEL_USE
!-------------------------------------------
    ALLOCATE(glob_fuel_use(nlon, nlev, ngl))  
    ALLOCATE(add_glob_fuel_use(nwaypoints-1,4)) !1:lon, 2:lat, 3:alt, 4:Emissions
    ALLOCATE(zalt(nlev))
    glob_fuel_use = 0.0_dp   
    !write(*,*)"FUEL_check_used_PEs_1=",p_pe,"p_io=",p_io
    do i_traj=1,nlroutes 
    !write(*,*)"FUEL_check_used_PEs_2=",p_pe,"p_io=",p_io
        if (p2_ac_routes_desc(i_traj,AC_status)==arrived .or.     &
            p2_ac_routes_desc(i_traj,AC_status)==inflight)then   !&
            call fill_emission_field(FUEL,ac_routes(:,:,i_traj,1),                   &
                                     p2_ac_routes_desc(i_traj,AC_pos),               &
                                     p2_ac_routes_desc(i_traj,AC_pos_old),add_glob_fuel_use(:,:),i_add_fuel_use)
!                                    p2_ac_routes_desc(i_traj,AC_pos_old),glob_fuel_use)
!           write(*,*)'i_add_fuel_use(SMIL)=',i_add_fuel_use,'p_pe=',p_pe,'p_io=',p_io
           IF(i_add_fuel_use/=0)THEN 
              DO i_add=1,i_add_fuel_use
                 !HORIZONTAL INDEX 
!                 write(*,*)'(FUEL)check_before_locate','p_pe=',p_pe,'p_io=',p_io
!                 write(*,*)'add_glob_fuel_use(SMIL)1 =',add_glob_fuel_use(i_add,:)
                 CALL locate_in_decomp(status, add_glob_fuel_use(i_add,1), add_glob_fuel_use(i_add,2), &
                                       d_pe, d_jp, d_jrow, d_jgx, d_jgy)   !20140606-1
                 !VERTICAL INDEX 
                 zalt(:)=zgl_geopot_3d(d_jgx,:,d_jgy)/g   !zalt:geopotential altitude,20131105-1
!                 write(*,*)'(FUEL)zalt=',zalt(:)
                 CALL nn_index(zalt(:), add_glob_fuel_use(i_add,3), d_idx)
                 glob_fuel_use(d_jgx,d_idx,d_jgy)=add_glob_fuel_use(i_add,4) 
!                 write(*,*)'(FUEL)d_jgx, d_idx, d_jgy, g =',d_jgx,d_idx,d_jgy,g
!                 write(*,*)'(FUEL)i_add_fuel_use(SMIL)=',i_add_fuel_use
!                 write(*,*)'(FUEL)add_glob_fuel_use(SMIL)2 =',add_glob_fuel_use(i_add,:)
!                 write(*,*)'(FUEL)glob_fuel_use(SMIL) =',glob_fuel_use(d_jgx,d_idx,d_jgy)
              ENDDO
           ENDIF
        endif
    enddo     
!    write(*,*)'(FUEL)nlon,nlev,ngl=',nlon, nlev, ngl
!    write(*,*)"(FUEL)check_used_PEs3=",p_pe,"p_io=",p_io
    CALL trp_gpdc_gpgl(-1, fuel_use, glob_fuel_use, M_SUM)
    DEALLOCATE(glob_fuel_use); NULLIFY(glob_fuel_use)
    DEALLOCATE(add_glob_fuel_use)
    DEALLOCATE(zalt) 

!Yin_20170423 
!-------------------------------------------
!  contrail potential coverage
!-------------------------------------------
    ALLOCATE(glob_CPC(nlon, nlev, ngl))
    ALLOCATE(add_glob_CPC(nwaypoints-1,4)) !1:lon, 2:lat, 3:alt, 4:Emissions
    ALLOCATE(zalt(nlev))
    glob_CPC = 0.0_dp
    !write(*,*)"FUEL_check_used_PEs_1=",p_pe,"p_io=",p_io
    do i_traj=1,nlroutes
    !write(*,*)"FUEL_check_used_PEs_2=",p_pe,"p_io=",p_io
        if (p2_ac_routes_desc(i_traj,AC_status)==arrived .or.     &
            p2_ac_routes_desc(i_traj,AC_status)==inflight)then   !&
            call fill_emission_field(C_PC,ac_routes(:,:,i_traj,1),                   &
                                     p2_ac_routes_desc(i_traj,AC_pos),               &
                                     p2_ac_routes_desc(i_traj,AC_pos_old),add_glob_CPC(:,:),i_add_CPC)
!                                    p2_ac_routes_desc(i_traj,AC_pos_old),glob_fuel_use)
!           write(*,*)'i_add_CPC(SMIL)=',i_add_CPC,'p_pe=',p_pe,'p_io=',p_io
           IF(i_add_CPC/=0)THEN
              DO i_add=1,i_add_CPC
                 !HORIZONTAL INDEX 
!                 write(*,*)'(CPC)check_before_locate','p_pe=',p_pe,'p_io=',p_io
!                 write(*,*)'add_glob_CPC(SMIL)1 =',add_glob_CPC(i_add,:)
                 CALL locate_in_decomp(status, add_glob_CPC(i_add,1), add_glob_CPC(i_add,2), &
                                       d_pe, d_jp, d_jrow, d_jgx, d_jgy)   !20140606-1
                 !VERTICAL INDEX 
                 zalt(:)=zgl_geopot_3d(d_jgx,:,d_jgy)/g   !zalt:geopotential altitude,20131105-1
!                 write(*,*)'(CPC)zalt=',zalt(:)
                 CALL nn_index(zalt(:), add_glob_CPC(i_add,3), d_idx)
                 glob_CPC(d_jgx,d_idx,d_jgy)=add_glob_CPC(i_add,4)
!                 write(*,*)'(CPC)d_jgx, d_idx, d_jgy, g =',d_jgx,d_idx,d_jgy,g
!                 write(*,*)'(CPC)i_add_CPC(SMIL)=',i_add_CPC
!                 write(*,*)'(CPC)add_glob_CPC(SMIL)2 =',add_glob_CPC(i_add,:)
!                 write(*,*)'(CPC)glob_CPC(SMIL) =',glob_CPC(d_jgx,d_idx,d_jgy)
              ENDDO
           ENDIF
        endif
    enddo
!    write(*,*)'(CPC)nlon,nlev,ngl=',nlon, nlev, ngl
!    write(*,*)"(CPC)check_used_PEs3=",p_pe,"p_io=",p_io
    CALL trp_gpdc_gpgl(-1, CPC_traj, glob_CPC, M_SUM)
    DEALLOCATE(glob_CPC); NULLIFY(glob_CPC)
    DEALLOCATE(add_glob_CPC)
    DEALLOCATE(zalt)
!Yin_20170423

!Yin_20170801+ 
!-------------------------------------------
!  ATR20 ozone
!-------------------------------------------
    ALLOCATE(glob_ATR20O3(nlon, nlev, ngl))
    ALLOCATE(add_glob_ATR20O3(nwaypoints-1,4)) !1:lon, 2:lat, 3:alt, 4:Emissions
    ALLOCATE(zalt(nlev))
    glob_ATR20O3 = 0.0_dp
    !write(*,*)"FUEL_check_used_PEs_1=",p_pe,"p_io=",p_io
    do i_traj=1,nlroutes
    !write(*,*)"FUEL_check_used_PEs_2=",p_pe,"p_io=",p_io
        if (p2_ac_routes_desc(i_traj,AC_status)==arrived .or.     &
            p2_ac_routes_desc(i_traj,AC_status)==inflight)then   !&
            call fill_emission_field(ATR20O3,ac_routes(:,:,i_traj,1),                   &
                                     p2_ac_routes_desc(i_traj,AC_pos),               &
                                     p2_ac_routes_desc(i_traj,AC_pos_old),add_glob_ATR20O3(:,:),i_add_ATR20O3)
!                                    p2_ac_routes_desc(i_traj,AC_pos_old),glob_fuel_use)
!           write(*,*)'i_add_ATR20O3(SMIL)=',i_add_ATR20O3,'p_pe=',p_pe,'p_io=',p_io
           IF(i_add_ATR20O3/=0)THEN
              DO i_add=1,i_add_ATR20O3
                 !HORIZONTAL INDEX 
!                 write(*,*)'(ATR20O3)check_before_locate','p_pe=',p_pe,'p_io=',p_io
!                 write(*,*)'add_glob_ATR20O3(SMIL)1 =',add_glob_ATR20O3(i_add,:)
                 CALL locate_in_decomp(status, add_glob_ATR20O3(i_add,1), add_glob_ATR20O3(i_add,2), &
                                       d_pe, d_jp, d_jrow, d_jgx, d_jgy)   !20140606-1
                 !VERTICAL INDEX 
                 zalt(:)=zgl_geopot_3d(d_jgx,:,d_jgy)/g   !zalt:geopotential altitude,20131105-1
!                 write(*,*)'(ATR20O3)zalt=',zalt(:)
                 CALL nn_index(zalt(:), add_glob_ATR20O3(i_add,3), d_idx)
                 glob_ATR20O3(d_jgx,d_idx,d_jgy)=add_glob_ATR20O3(i_add,4)
!                 write(*,*)'(ATR20O3)d_jgx, d_idx, d_jgy, g =',d_jgx,d_idx,d_jgy,g
!                 write(*,*)'(ATR20O3)i_add_ATR20O3(SMIL)=',i_add_ATR20O3
!                 write(*,*)'(ATR20O3)add_glob_ATR20O3(SMIL) =',add_glob_ATR20O3(i_add,:)
!                 write(*,*)'(ATR20O3)glob_ATR20O3(SMIL) =',glob_ATR20O3(d_jgx,d_idx,d_jgy)
              ENDDO
           ENDIF
        endif
    enddo
!    write(*,*)'(ATR20O3)nlon,nlev,ngl=',nlon, nlev, ngl
!    write(*,*)"(ATR20O3)check_used_PEs3=",p_pe,"p_io=",p_io
    CALL trp_gpdc_gpgl(-1, ATR20_O3_traj, glob_ATR20O3, M_SUM)
    DEALLOCATE(glob_ATR20O3); NULLIFY(glob_ATR20O3)
    DEALLOCATE(add_glob_ATR20O3)
    DEALLOCATE(zalt)
!-------------------------------------------
!  ATR20 methane
!-------------------------------------------
    ALLOCATE(glob_ATR20CH4(nlon, nlev, ngl))
    ALLOCATE(add_glob_ATR20CH4(nwaypoints-1,4)) !1:lon, 2:lat, 3:alt, 4:Emissions
    ALLOCATE(zalt(nlev))
    glob_ATR20CH4 = 0.0_dp
    !write(*,*)"FUEL_check_used_PEs_1=",p_pe,"p_io=",p_io
    do i_traj=1,nlroutes
    !write(*,*)"FUEL_check_used_PEs_2=",p_pe,"p_io=",p_io
        if (p2_ac_routes_desc(i_traj,AC_status)==arrived .or.     &
            p2_ac_routes_desc(i_traj,AC_status)==inflight)then   !&
            call fill_emission_field(ATR20CH4,ac_routes(:,:,i_traj,1),                   &
                                     p2_ac_routes_desc(i_traj,AC_pos),               &
                                     p2_ac_routes_desc(i_traj,AC_pos_old),add_glob_ATR20CH4(:,:),i_add_ATR20CH4)
!           write(*,*)'i_add_ATR20CH4(SMIL)=',i_add_ATR20CH4,'p_pe=',p_pe,'p_io=',p_io
           IF(i_add_ATR20CH4/=0)THEN
              DO i_add=1,i_add_ATR20CH4
                 !HORIZONTAL INDEX 
!                 write(*,*)'(ATR20CH4)check_before_locate','p_pe=',p_pe,'p_io=',p_io
!                 write(*,*)'add_glob_ATR20CH4(SMIL)1 =',add_glob_ATR20CH4(i_add,:)
                 CALL locate_in_decomp(status, add_glob_ATR20CH4(i_add,1), add_glob_ATR20CH4(i_add,2), &
                                       d_pe, d_jp, d_jrow, d_jgx, d_jgy)   !20140606-1
                 !VERTICAL INDEX 
                 zalt(:)=zgl_geopot_3d(d_jgx,:,d_jgy)/g   !zalt:geopotential altitude,20131105-1
!                 write(*,*)'(ATR20CH4)zalt=',zalt(:)
                 CALL nn_index(zalt(:), add_glob_ATR20CH4(i_add,3), d_idx)
                 glob_ATR20CH4(d_jgx,d_idx,d_jgy)=add_glob_ATR20CH4(i_add,4)
!                 write(*,*)'(ATR20CH4)d_jgx, d_idx, d_jgy, g =',d_jgx,d_idx,d_jgy,g
!                 write(*,*)'(ATR20CH4)i_add_ATR20CH4(SMIL)=',i_add_ATR20CH4
!                 write(*,*)'(ATR20CH4)add_glob_ATR20CH4(SMIL) =',add_glob_ATR20CH4(i_add,:)
!                 write(*,*)'(ATR20CH4)glob_ATR20CH4(SMIL) =',glob_ATR20CH4(d_jgx,d_idx,d_jgy)
              ENDDO
           ENDIF
        endif
    enddo
!    write(*,*)'(ATR20CH4)nlon,nlev,ngl=',nlon, nlev, ngl
!    write(*,*)"(ATR20CH4)check_used_PEs3=",p_pe,"p_io=",p_io
    CALL trp_gpdc_gpgl(-1, ATR20_CH4_traj, glob_ATR20CH4, M_SUM)
    DEALLOCATE(glob_ATR20CH4); NULLIFY(glob_ATR20CH4)
    DEALLOCATE(add_glob_ATR20CH4)
    DEALLOCATE(zalt)
!-------------------------------------------
!  ATR20 h2o
!-------------------------------------------
    ALLOCATE(glob_ATR20H2O(nlon, nlev, ngl))
    ALLOCATE(add_glob_ATR20H2O(nwaypoints-1,4)) !1:lon, 2:lat, 3:alt, 4:Emissions
    ALLOCATE(zalt(nlev))
    glob_ATR20H2O = 0.0_dp
    !write(*,*)"FUEL_check_used_PEs_1=",p_pe,"p_io=",p_io
    do i_traj=1,nlroutes
    !write(*,*)"FUEL_check_used_PEs_2=",p_pe,"p_io=",p_io
        if (p2_ac_routes_desc(i_traj,AC_status)==arrived .or.     &
            p2_ac_routes_desc(i_traj,AC_status)==inflight)then   !&
            call fill_emission_field(ATR20H2O,ac_routes(:,:,i_traj,1),                   &
                                     p2_ac_routes_desc(i_traj,AC_pos),               &
                                     p2_ac_routes_desc(i_traj,AC_pos_old),add_glob_ATR20H2O(:,:),i_add_ATR20H2O)
!           write(*,*)'i_add_ATR20H2O(SMIL)=',i_add_ATR20H2O,'p_pe=',p_pe,'p_io=',p_io
           IF(i_add_ATR20H2O/=0)THEN
              DO i_add=1,i_add_ATR20H2O
                 !HORIZONTAL INDEX 
!                 write(*,*)'(ATR20H2O)check_before_locate','p_pe=',p_pe,'p_io=',p_io
!                 write(*,*)'add_glob_ATR20H2O(SMIL) =',add_glob_ATR20H2O(i_add,:)
                 CALL locate_in_decomp(status, add_glob_ATR20H2O(i_add,1), add_glob_ATR20H2O(i_add,2), &
                                       d_pe, d_jp, d_jrow, d_jgx, d_jgy)   !20140606-1
                 !VERTICAL INDEX 
                 zalt(:)=zgl_geopot_3d(d_jgx,:,d_jgy)/g   !zalt:geopotential altitude,20131105-1
!                 write(*,*)'(ATR20H2O)zalt=',zalt(:)
                 CALL nn_index(zalt(:), add_glob_ATR20H2O(i_add,3), d_idx)
                 glob_ATR20H2O(d_jgx,d_idx,d_jgy)=add_glob_ATR20H2O(i_add,4)
!                 write(*,*)'(ATR20H2O)d_jgx, d_idx, d_jgy, g =',d_jgx,d_idx,d_jgy,g
!                 write(*,*)'(ATR20H2O)i_add_ATR20H2O(SMIL)=',i_add_ATR20H2O
!                 write(*,*)'(ATR20H2O)add_glob_ATR20H2O(SMIL) =',add_glob_ATR20H2O(i_add,:)
!                 write(*,*)'(ATR20H2O)glob_ATR20H2O(SMIL) =',glob_ATR20H2O(d_jgx,d_idx,d_jgy)
              ENDDO
           ENDIF
        endif
    enddo
!    write(*,*)'(ATR20H2O)nlon,nlev,ngl=',nlon, nlev, ngl
!    write(*,*)"(ATR20H2O)check_used_PEs3=",p_pe,"p_io=",p_io
    CALL trp_gpdc_gpgl(-1, ATR20_H2O_traj, glob_ATR20H2O, M_SUM)
    DEALLOCATE(glob_ATR20H2O); NULLIFY(glob_ATR20H2O)
    DEALLOCATE(add_glob_ATR20H2O)
    DEALLOCATE(zalt)
!-------------------------------------------
!  ATR20 CPC
!-------------------------------------------
    ALLOCATE(glob_ATR20CPC(nlon, nlev, ngl))
    ALLOCATE(add_glob_ATR20CPC(nwaypoints-1,4)) !1:lon, 2:lat, 3:alt, 4:Emissions
    ALLOCATE(zalt(nlev))
    glob_ATR20CPC = 0.0_dp
    !write(*,*)"FUEL_check_used_PEs_1=",p_pe,"p_io=",p_io
    do i_traj=1,nlroutes
    !write(*,*)"FUEL_check_used_PEs_2=",p_pe,"p_io=",p_io
        if (p2_ac_routes_desc(i_traj,AC_status)==arrived .or.     &
            p2_ac_routes_desc(i_traj,AC_status)==inflight)then   !&
            call fill_emission_field(ATR20CPC,ac_routes(:,:,i_traj,1),                   &
                                     p2_ac_routes_desc(i_traj,AC_pos),               &
                                     p2_ac_routes_desc(i_traj,AC_pos_old),add_glob_ATR20CPC(:,:),i_add_ATR20CPC)
!           write(*,*)'i_add_ATR20CPC(SMIL)=',i_add_ATR20CPC,'p_pe=',p_pe,'p_io=',p_io
           IF(i_add_ATR20CPC/=0)THEN
              DO i_add=1,i_add_ATR20CPC
                 !HORIZONTAL INDEX 
!                 write(*,*)'(ATR20CPC)check_before_locate','p_pe=',p_pe,'p_io=',p_io
!                 write(*,*)'add_glob_ATR20CPC(SMIL) =',add_glob_ATR20CPC(i_add,:)
                 CALL locate_in_decomp(status, add_glob_ATR20CPC(i_add,1), add_glob_ATR20CPC(i_add,2), &
                                       d_pe, d_jp, d_jrow, d_jgx, d_jgy)   !20140606-1
                 !VERTICAL INDEX 
                 zalt(:)=zgl_geopot_3d(d_jgx,:,d_jgy)/g   !zalt:geopotential altitude,20131105-1
!                 write(*,*)'(ATR20CPC)zalt=',zalt(:)
                 CALL nn_index(zalt(:), add_glob_ATR20CPC(i_add,3), d_idx)
                 glob_ATR20CPC(d_jgx,d_idx,d_jgy)=add_glob_ATR20CPC(i_add,4)
!                 write(*,*)'(ATR20CPC)d_jgx, d_idx, d_jgy, g =',d_jgx,d_idx,d_jgy,g
!                 write(*,*)'(ATR20CPC)i_add_ATR20CPC(SMIL)=',i_add_ATR20CPC
!                 write(*,*)'(ATR20CPC)add_glob_ATR20CPC(SMIL) =',add_glob_ATR20CPC(i_add,:)
!                 write(*,*)'(ATR20CPC)glob_ATR20CPC(SMIL) =',glob_ATR20CPC(d_jgx,d_idx,d_jgy)
              ENDDO
           ENDIF
        endif
    enddo
!    write(*,*)'(ATR20CPC)nlon,nlev,ngl=',nlon, nlev, ngl
!    write(*,*)"(ATR20CPC)check_used_PEs3=",p_pe,"p_io=",p_io
    CALL trp_gpdc_gpgl(-1, ATR20_CPC_traj, glob_ATR20CPC, M_SUM)
    DEALLOCATE(glob_ATR20CPC); NULLIFY(glob_ATR20CPC)
    DEALLOCATE(add_glob_ATR20CPC)
    DEALLOCATE(zalt)
!-------------------------------------------
!  ATR20 CO2
!-------------------------------------------
    ALLOCATE(glob_ATR20CO2(nlon, nlev, ngl))
    ALLOCATE(add_glob_ATR20CO2(nwaypoints-1,4)) !1:lon, 2:lat, 3:alt, 4:Emissions
    ALLOCATE(zalt(nlev))
    glob_ATR20CO2 = 0.0_dp
    !write(*,*)"FUEL_check_used_PEs_1=",p_pe,"p_io=",p_io
    do i_traj=1,nlroutes
    !write(*,*)"FUEL_check_used_PEs_2=",p_pe,"p_io=",p_io
        if (p2_ac_routes_desc(i_traj,AC_status)==arrived .or.     &
            p2_ac_routes_desc(i_traj,AC_status)==inflight)then   !&
            call fill_emission_field(ATR20CO2,ac_routes(:,:,i_traj,1),                   &
                                     p2_ac_routes_desc(i_traj,AC_pos),               &
                                     p2_ac_routes_desc(i_traj,AC_pos_old),add_glob_ATR20CO2(:,:),i_add_ATR20CO2)
!           write(*,*)'i_add_ATR20CO2(SMIL)=',i_add_ATR20CO2,'p_pe=',p_pe,'p_io=',p_io
           IF(i_add_ATR20CO2/=0)THEN
              DO i_add=1,i_add_ATR20CO2
                 !HORIZONTAL INDEX 
!                 write(*,*)'(ATR20CO2)check_before_locate','p_pe=',p_pe,'p_io=',p_io
!                 write(*,*)'add_glob_ATR20CO2(SMIL) =',add_glob_ATR20CO2(i_add,:)
                 CALL locate_in_decomp(status, add_glob_ATR20CO2(i_add,1), add_glob_ATR20CO2(i_add,2), &
                                       d_pe, d_jp, d_jrow, d_jgx, d_jgy)   !20140606-1
                 !VERTICAL INDEX 
                 zalt(:)=zgl_geopot_3d(d_jgx,:,d_jgy)/g   !zalt:geopotential altitude,20131105-1
!                 write(*,*)'(ATR20CO2)zalt=',zalt(:)
                 CALL nn_index(zalt(:), add_glob_ATR20CO2(i_add,3), d_idx)
                 glob_ATR20CO2(d_jgx,d_idx,d_jgy)=add_glob_ATR20CO2(i_add,4)
!                 write(*,*)'(ATR20CO2)d_jgx, d_idx, d_jgy, g =',d_jgx,d_idx,d_jgy,g
!                 write(*,*)'(ATR20CO2)i_add_ATR20CO2(SMIL)=',i_add_ATR20CO2
!                 write(*,*)'(ATR20CO2)add_glob_ATR20CO2(SMIL) =',add_glob_ATR20CO2(i_add,:)
!                 write(*,*)'(ATR20CO2)glob_ATR20CO2(SMIL) =',glob_ATR20CO2(d_jgx,d_idx,d_jgy)
              ENDDO
           ENDIF
        endif
    enddo
!    write(*,*)'(ATR20CO2)nlon,nlev,ngl=',nlon, nlev, ngl
!    write(*,*)"(ATR20CO2)check_used_PEs3=",p_pe,"p_io=",p_io
    CALL trp_gpdc_gpgl(-1, ATR20_CO2_traj, glob_ATR20CO2, M_SUM)
    DEALLOCATE(glob_ATR20CO2); NULLIFY(glob_ATR20CO2)
    DEALLOCATE(add_glob_ATR20CO2)
    DEALLOCATE(zalt)

!-------------------------------------------
!  ATR20 TOT
!-------------------------------------------
    ALLOCATE(glob_ATR20TOT(nlon, nlev, ngl))
    ALLOCATE(add_glob_ATR20TOT(nwaypoints-1,4)) !1:lon, 2:lat, 3:alt, 4:Emissions
    ALLOCATE(zalt(nlev))
    glob_ATR20TOT = 0.0_dp
    !write(*,*)"FUEL_check_used_PEs_1=",p_pe,"p_io=",p_io
    do i_traj=1,nlroutes
    !write(*,*)"FUEL_check_used_PEs_2=",p_pe,"p_io=",p_io
        if (p2_ac_routes_desc(i_traj,AC_status)==arrived .or.     &
            p2_ac_routes_desc(i_traj,AC_status)==inflight)then   !&
            call fill_emission_field(ATR20TOT,ac_routes(:,:,i_traj,1),                   &
                                     p2_ac_routes_desc(i_traj,AC_pos),               &
                                     p2_ac_routes_desc(i_traj,AC_pos_old),add_glob_ATR20TOT(:,:),i_add_ATR20TOT)
!           write(*,*)'i_add_ATR20TOT(SMIL)=',i_add_ATR20TOT,'p_pe=',p_pe,'p_io=',p_io
           IF(i_add_ATR20TOT/=0)THEN
              DO i_add=1,i_add_ATR20TOT
                 !HORIZONTAL INDEX 
!                 write(*,*)'(ATR20TOT)check_before_locate','p_pe=',p_pe,'p_io=',p_io
!                 write(*,*)'add_glob_ATR20TOT(SMIL) =',add_glob_ATR20TOT(i_add,:)
                 CALL locate_in_decomp(status, add_glob_ATR20TOT(i_add,1), add_glob_ATR20TOT(i_add,2), &
                                       d_pe, d_jp, d_jrow, d_jgx, d_jgy)   !20140606-1
                 !VERTICAL INDEX 
                 zalt(:)=zgl_geopot_3d(d_jgx,:,d_jgy)/g   !zalt:geopotential altitude,20131105-1
!                 write(*,*)'(ATR20TOT)zalt=',zalt(:)
                 CALL nn_index(zalt(:), add_glob_ATR20TOT(i_add,3), d_idx)
                 glob_ATR20TOT(d_jgx,d_idx,d_jgy)=add_glob_ATR20TOT(i_add,4)
!                 write(*,*)'(ATR20TOT)d_jgx, d_idx, d_jgy, g =',d_jgx,d_idx,d_jgy,g
!                 write(*,*)'(ATR20TOT)i_add_ATR20TOT(SMIL)=',i_add_ATR20TOT
!                 write(*,*)'(ATR20TOT)add_glob_ATR20TOT(SMIL) =',add_glob_ATR20TOT(i_add,:)
!                 write(*,*)'(ATR20TOT)glob_ATR20TOT(SMIL) =',glob_ATR20TOT(d_jgx,d_idx,d_jgy)
              ENDDO
           ENDIF
        endif
    enddo
!    write(*,*)'(ATR20TOT)nlon,nlev,ngl=',nlon, nlev, ngl
!    write(*,*)"(ATR20TOT)check_used_PEs3=",p_pe,"p_io=",p_io
    CALL trp_gpdc_gpgl(-1, ATR20_TOT_traj, glob_ATR20TOT, M_SUM)
    DEALLOCATE(glob_ATR20TOT); NULLIFY(glob_ATR20TOT)
    DEALLOCATE(add_glob_ATR20TOT)
    DEALLOCATE(zalt)

!Yin_20170801-


    !GMT and GNUPLOT: global flight routes visualization for time series.20140610-3
    !See line 1700 
    !IF(p_pe.eq.0)THEN
    !   OPEN(20,file='gmt.dat',position='append')
    !      write(20,*)
    !     write(20,*)
    !     write(20,'(a)')">"
    !   CLOSE(20)
    !ENDIF
!    IF(p_pe.eq.0)THEN
!       OPEN(10,file='gnu001.dat',position='append')
!          write(10,*)
       !  write(10,'(a)')">"
!       CLOSE(10)
!    ENDIF
!    IF(p_pe.eq.1)THEN
!       OPEN(20,file='gnu002.dat',position='append')
!          write(20,*)
       !  write(20,'(a)')">"
!       CLOSE(20)
!    ENDIF
!    IF(p_pe.eq.2)THEN
!       OPEN(30,file='gnu003.dat',position='append')
!          write(30,*)
       !  write(30,'(a)')">"
!       CLOSE(30)
!    ENDIF
!    IF(p_pe.eq.3)THEN
!       OPEN(40,file='gnu004.dat',position='append')
!          write(40,*)
       !  write(40,'(a)')">"
!       CLOSE(40)
!    ENDIF

! Release arrived trajectories.   
! After emission calculation is finished at A_city, AC_status, AC_pos, AC_pos_old 
! D_lon, D_lat, D_time, A_lon, and A_lat are changed into notused(= -1.0dp).
! The below 'do' structure is used every time step (12 min) on each PE. 
! See: 20120911-2,20121220-4 
    do i_traj=1,nlroutes 
        if(p2_ac_routes_desc(i_traj,AC_status)==arrived)then
!            write(*,*)'aircraft_number=',i_traj
!            write(*,*)'arrive_to_notused'
!            write(*,*)'B_p2_ac_routes_desc(i_traj,AC_status)=',p2_ac_routes_desc(i_traj,AC_status)
!            write(*,*)'B_p2_ac_routes_desc(i_traj,AC_pos)=',p2_ac_routes_desc(i_traj,AC_pos)
!            write(*,*)'B_p2_ac_routes_desc(i_traj,AC_pos_old)=',p2_ac_routes_desc(i_traj,AC_pos_old)
            !Initilization of p2_ac_routes_desc
            p2_ac_routes_desc(i_traj,AC_status) = notused
            p2_ac_routes_desc(i_traj,AC_pos)    = notused
            p2_ac_routes_desc(i_traj,AC_pos_old)= notused
            p2_ac_routes_desc(i_traj,D_lon)     = notused
            p2_ac_routes_desc(i_traj,D_lat)     = notused              
            p2_ac_routes_desc(i_traj,D_time)    = notused 
            p2_ac_routes_desc(i_traj,A_lon)     = notused              
            p2_ac_routes_desc(i_traj,A_lat)     = notused
            
            !Initilization of ac_routes
            ac_routes(:,:,i_traj,1) = 0.0_dp            

!            write(*,*)'A_p2_ac_routes_desc(i_traj,AC_status)=',p2_ac_routes_desc(i_traj,AC_status)
!            write(*,*)'A_p2_ac_routes_desc(i_traj,AC_pos)=',p2_ac_routes_desc(i_traj,AC_pos)
!            write(*,*)'A_p2_ac_routes_desc(i_traj,AC_pos_old)=',p2_ac_routes_desc(i_traj,AC_pos_old)
!            if(p_pe.eq.0)write(*,*)'used_here','d_now',now,'ac_routes(time)',ac_routes(:,4,i_traj,1) 
        endif
    enddo

    ! Deacllocate/nullify global fields: release memory
    ! see20121121-1
    DEALLOCATE(zgl_geopot_3d); NULLIFY(zgl_geopot_3d)
    DEALLOCATE(press_3d_g); NULLIFY(press_3d_g)
    DEALLOCATE(t_scb_g); NULLIFY(t_scb_g)
    DEALLOCATE(rho_air_dry_3d_g); NULLIFY(rho_air_dry_3d_g)
    DEALLOCATE(grvol_g); NULLIFY(grvol_g) !20180526
    !IF(option_traj_calc==Wind_opt)THEN   !20140524-10
    DEALLOCATE(uwind_g); NULLIFY(uwind_g)
    DEALLOCATE(vwind_g); NULLIFY(vwind_g)
    DEALLOCATE(v_z_g); NULLIFY(v_z_g)
    DEALLOCATE(v_z); NULLIFY(v_z)
    !ENDIF     !see20140524-10
    DEALLOCATE(cpc_g); NULLIFY(cpc_g)     !20170531-1
    DEALLOCATE(ATR20_o3_g); NULLIFY(ATR20_o3_g)     !20170531-1
    DEALLOCATE(ATR20_ch4_g); NULLIFY(ATR20_ch4_g)     !20170531-1
    DEALLOCATE(ATR20_h2o_g); NULLIFY(ATR20_h2o_g)     !20170531-1
    DEALLOCATE(ATR20_contrail_g); NULLIFY(ATR20_contrail_g)     !20170531-1
    DEALLOCATE(ATR20_co2_g); NULLIFY(ATR20_co2_g)     !20170531-1
! ----------------------------------------------------------------------
! This part is used every time step (e.g.12 min) on each PE.  
!     write(*,*)now,'calculation check, before end_global_end'
  END SUBROUTINE airtraf_global_end
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE airtraf_physc

    ! ------------------------------------------------------------------
    ! This subroutine is called within the time loop.
    ! Therefore, this subroutine is used every time step (12 min) by each PE.
    !
    ! It constitutes the main entry point for additional processes 
    ! or diagnostics.
    ! Here, only the current vector of the grid-point-fields is
    ! accessible.
    !
    ! ilroutes     : Total number of active routes(aircraft), related to
    !                ac_routes(nwaypoints, nprops, nlroutes, 1)
    !
    ! ilroutes_out : Total number of output routes(aircraft), related to
    !                ac_routes_out(nwaypoints, nprops, nlroutes_out, 1) 
    !
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    !USE messy_main_timer,         ONLY: lstart, time_step_len  !lresume
    USE messy_main_grid_def_mem_bi, ONLY: kproma, jrow !, nlev, philon_2d,philat_2d,vervel
    USE messy_main_data_bi,         ONLY: vervel=>vervel_3d ! vervel
    !USE messy_main_mpi_bi,        ONLY: p_parallel_io,p_io,p_bcast, p_pe
!Yin_20170321+
#ifndef MESSYTENDENCY
    USE messy_main_tracer_mem_bi, ONLY: pxtte => qxtte
#endif
!Yin_20170321-

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER      :: substr = 'airtraf_physc'
    INTEGER                          :: i
!    REAL(DP), DIMENSION(:), POINTER  :: my_tte => NULL() ! perturbation tendency

    INTEGER                          :: ilroutes, ilroutes_out,idt,jt

!20130911-1,20130916-3
    !Calculate vervel_3d from vervel.
    !The vervel_3d is used in airtraf_global_end.
    !vervel: decomposed vertical wind velocity [Pa/s] at box center.
    !        vervel(nproma,nlev) is 2D array,20130911-1
    !        vervel value in positive(air direction); nevative (ground direction)20130814-3 
    !verval_3d(nproma,nlev,ngpblks): decomposed vertical wind velocity [Pa/s] in 3D. 
!   write(*,*)'PE=',p_pe,'physc_shape_vervel_3d',shape(vervel_3d),vervel_3d(:,nlev,:)
    vervel_3d(1:kproma,:,jrow)=vervel(1:kproma,:,jrow)
!   write(*,*)'PE=',p_pe,'jrow=',jrow,'kproma=',kproma,'nlev=',nlev,'physic_vervel_3dcheck'

! see note 20120904-5, 20120905-1,-2, 20120830-5,-6
!   write(*,*)'calculation check, after airtraf_physc,(PE=)',p_pe
!   write(*,*)'nlroutes,nlroutes_out',nlroutes,nlroutes_out
! this initilization is correct? see 20120830-5
    ilroutes_out = 0
    ilroutes     = 0
! note 2012.05.09-4, 20120724-3
    do i=1,nlroutes_out  !see note 20120830-5
!!!   if (p3_ac_routes_out(i,1,1).gt.0._dp)ilroutes_out=ilroutes_out+1
!HY See 20120831-2, this if structure needs modification!
!     if (p3_ac_routes_out(1,1,i).gt.0._dp)ilroutes_out=ilroutes_out+1
      if (p3_ac_routes_out(1,4,i).gt.0._dp)ilroutes_out=ilroutes_out+1
    enddo

    do i=1,nlroutes !see note 20120830-5
!!!   if (p3_ac_routes(i,1,1).gt.0._dp)ilroutes=ilroutes+1
!HY 20120831-2
!     if (p3_ac_routes(1,1,i).gt.0._dp)ilroutes=ilroutes+1
      if (p3_ac_routes(1,4,i).gt.0._dp)ilroutes=ilroutes+1
    enddo
!!! write(*,*) "Total number of active/output routes ILROUTES=",ilroutes,ilroutes_out,p3_ac_routes(1:nlroutes,1,1)
!   write(*,*) "PE=",p_pe,"Total number of active/output routes ILROUTES=",ilroutes,ilroutes_out,p3_ac_routes(1,1,1:nlroutes)
!   write(*,*) "PE=",p_pe,'ac_routes(time)_physc',ac_routes(:,4,1:5,1) 

    ! This is dummy for arrival check!!! See note 20120509-4
    do i=1,nlroutes   
!!!    p3_ac_routes(i,:,:)=p3_ac_routes(i,:,:)+i+real(p_pe,kind=dp)
!HY 20120831-3       p3_ac_routes(:,:,i)=p3_ac_routes(:,:,i)+real(i,kind=dp)+real(p_pe,kind=dp)
!       p3_ac_routes(:,:,i)=p3_ac_routes(:,:,i)
!!!    if (p3_ac_routes(i,1,1).gt.100._dp) then 
!HY 20120831-3
!      if (p3_ac_routes(1,1,i).gt.100._dp)then   ! note 20120509-1,-4, if .gt.100, aircrafts are changed into arrivali status 
       if (p3_ac_routes(1,4,i).gt.0.0_dp)then    !  
          ilroutes_out = ilroutes_out + 1         ! Counting the additional arrival aircrafts(maybe dummy).
!          write(*,*) "Nextoutput ilroutes_out,ilroute",ilroutes_out, i
          if (ilroutes_out.gt.nlroutes_out)then
             write(*,*) "ILROUTES_OUT,nlroutes_out=",ilroutes_out, nlroutes_out ! Here, ilroutes_ous is 461 now. why??
             ilroutes_out = nlroutes_out
!!!          STOP "ilroutes_out.gt.nlroutes_out"
          endif
!!!       p3_ac_routes_out(ilroutes_out,1:nwaypoints_out,:)=p3_ac_routes(i,1:nwaypoints_out,:)
!!!       p3_ac_routes(i,:,:)=0.
        !!! Now, nwaypoints = nwaypoints was set. so no error occurs.
        !!! If nwaypoints=/ nwaypoints, error occurs.
        !!! See 2012.07.03-1,4 memo
!!!       p3_ac_routes_out(1:nwaypoints_out,:,ilroutes_out)=p3_ac_routes(1:nwaypoints_out,:,i)
!20120910-5,          p3_ac_routes_out(1:nwaypoints_out,:,ilroutes_out)=p3_ac_routes(1:nwaypoints,:,i)
!HY 20120831-3
!         p3_ac_routes(:,:,i)=0.
       endif
    enddo
!!! p3_ac_routes(2,:,:)=22._dp
!Yin_20170321+
    IF(nairtraf_nox_trac_gp>0)THEN
       DO jt=1,nairtraf_nox_trac_gp
          idt = idt_list_teNOx(jt)
         ! write(*,*)'idt=',idt
#ifndef MESSYTENDENCY
        !  pxtte(1:kproma,:,jrow) = pxtte(1:kproma,:,jrow)+teNOx_emis(1:kproma,:,jrow)
          pxtte(1:kproma,:,idt_list_teNOx(jt)) = pxtte(1:kproma,:,idt_list_teNOx(jt))+teNOx(1:kproma,:,jrow)
#else
!          CALL mtend_add_l(my_handle,mtend_id_tracer,&
!                           px=teNOx(1:kproma,:,jrow),idt=idt)
          CALL mtend_add_l(my_handle,idt, px=teNOx(1:kproma,:,jrow))
!HY_20180524+              px=teNOx(1:krpoma,:,jrow),idt=idt)
#endif
       ENDDO
    ENDIF
!Yin_20170321-

!   write(*,*)'calculation check, before end_airtraf_physc,(PE=)',p_pe

  END SUBROUTINE airtraf_physc
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE airtraf_free_memory

    ! ------------------------------------------------------------------
    ! This subroutine is used to deallocate the memory, which has
    ! been "manually" allocated in airtraf_init_memory.
    ! Note: channel object memory must not be deallocated! This is
    !       performed centrally.
    !
    ! This subroutine is executed after finishing all calculations of
    ! airtraf_globa_start and airtraf_physc. Also, this is executed 
    ! after outputting outputs- and restarts-files, when restart is performed.
    ! ------------------------------------------------------------------

    USE messy_main_mpi_bi,          ONLY: p_pe

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'airtraf_free_memory'

    DEALLOCATE(ac_routes)       ; NULLIFY(ac_routes) 
    DEALLOCATE(ac_routes_out)   ; NULLIFY(ac_routes_out)
!!   DEALLOCATE(flightplan)      ; NULLIFY(flightplan)

    ! Added by HY
    DEALLOCATE(ac_fp_info_r)    ; NULLIFY(ac_fp_info_r)
    DEALLOCATE(ac_fp_info_i)    ; NULLIFY(ac_fp_info_i)
!!HY 20120717    DEALLOCATE(ac_fp_info_city) ; NULLIFY(ac_fp_info_city)
!!HY 20120717    DEALLOCATE(ac_fp_info_c)    ; NULLIFY(ac_fp_info_c)
    DEALLOCATE(ac_routes_desc)  ; NULLIFY(ac_routes_desc)
!!   DEALLOCATE(p1_ac_fp_h)      ; NULLIFY(p1_ac_fp_h)  !maybe we should use this sentence.
    DEALLOCATE(vervel_3d)       ; NULLIFY(vervel_3d) !20130911-1 

    write(*,*)'airtraf_free_memory_check(PE=)',p_pe
    !Yin_20170321+
    IF(ASSOCIATED(idt_list_teNOx)) DEALLOCATE(idt_list_teNOx)
!HY_20180524
!    IF(ASSOCIATED(teNOx))  DEALLOCATE(teNOx)
    !Yin_20170321-

  END SUBROUTINE airtraf_free_memory
  ! ====================================================================

  ! ####################################################################
  ! PRIVATE SUBROUTINES
  ! ####################################################################

  ! ====================================================================
  SUBROUTINE airtraf_read_nml_cpl(status, iou)
   
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

    NAMELIST /CPL/ flightplan_filename, &
                   flightplan_varname,  &
                   lfeedback_NO,        &
                   lfeedback_H2O,       &
                   c_potcov,            &
                   c_atr20_o3,          &
                   c_atr20_ch4,         &
                   c_atr20_h2o,         &
                   c_atr20_contrail,    &
                   c_atr20_co2      



    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='airtraf_read_nml_cpl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status
    INTEGER                     :: i

    ! INITIALIZE
    status = 1 !ERROR

    ! Set default values
    lfeedback_NO  = .false.
    lfeedback_H2O = .false.

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! ### ADD HERE DIAGNOSTIC OUTPUT FOR LOG-FILE
    WRITE(*,*) 'Filename for Flightplan: ',TRIM(flightplan_filename)
    do i=1,total_fp_input
       WRITE(*,*) "Required variable :",flightplan_vardescription(i),"  ---> ",TRIM(flightplan_varname(i))
    enddo
    if (lfeedback_NO) then
       WRITE(*,*) "Emissions are fed back to other modules for NO"
    else 
       WRITE(*,*) "Emissions are NOT fed back to other modules for NO"
    endif
    if (lfeedback_H2O) then
       WRITE(*,*) "Emissions are fed back to other modules for H2O"
    else 
       WRITE(*,*) "Emissions are NOT fed back to other modules for H2O"
    endif
    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE airtraf_read_nml_cpl
  ! ====================================================================

  SUBROUTINE airtraf_write_output
     
  END SUBROUTINE airtraf_write_output

! **********************************************************************
END MODULE messy_airtraf_e5
! **********************************************************************
