MODULE messy_mpiom_e5

  !
  !  MESSy- submodel interface for MPIOM .
  !
  !  AUTHOR:  Pozzer Andrea, MPICH, May 2007
  !
  !  MESSy Interface for the Submodel kernel
  !

  ! MESSy
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi &
                                      , error_bi ! op_pj_20161104
  USE messy_main_tools,         ONLY: PTR_2D_ARRAY,PTR_1D_ARRAY
  USE messy_main_channel,       ONLY: STRLEN_CHANNEL, STRLEN_OBJECT, &
                                      DIMID_UNDEF
  USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM
  USE messy_main_timer_bi,       ONLY: p_bcast_event, timer_event_init
  USE messy_main_timer_event,    ONLY: io_time_event, TRIG_FIRST            &
                                     , time_event                           &
                                     , TIME_INC_DAYS, TIME_INC_YEARS        &
                                     , TIME_INC_MONTHS, TIME_INC_SECONDS
  USE messy_mpiom

  !MPIOM-core
  ! mz_bk_20110220+
#if defined (MPIOM_13B)
  USE messy_mpiom_mem_e5
#elif defined (MPIOM_2000)
  USE messy_mpiom_mem_e5,        ONLY: ie, je, ie_g, je_g, ke, dt, dz       &
                                     , imean, nprocx, nprocy                &
                                     , lgmdiag, lforcediag, lhfldiag        &
                                     , lconvdiag, ldiffdiag, lcalcdifi      &
                                     , lisopyc, ladpo, ibolk, iocad         &
                                     , nfixYearLen, cah00, aus, bofric      &
                                     , rayfric, av0, dv0, cstabeps, crelsal &
                                     , creltem, numriv, numglac             &
                                     , iaufr, iaufw, istart, i3drest        &
                                     , ioasisflux, cdvocon, cavocon         &
                                     , h0, hmin, rleadclose, armin, armax   &
                                     , hsntoice, sicthmin, sice, isnflg     &
                                     , ioconv, icontro, lnonblock, ltidal   &
                                     , ltidal_diag, lswr_jerlov             &
                                     , lfb_bgc_oce, jerlov_bluefrac         &
                                     , lwith_one_layer_shelfs, lsaoclose    &
                                     , iter_sor, iter_sor_hack, rtsorpar    &
                                     , ladfs, ibbl_transport, iocaduv       &
                                     , lmpitype, jerlov_atten               &
                                     , rtsorpar_hack, saf, taf, dzw, dti    &
                                     , spongezone, ndtday, aulapts, aulapuv &
                                     , ah00, tiestw                         &
                                     , luse_river_runoff_stations           &
                                     , luse_glac_calv, almzer, conn, stabn  &
                                     , fclou, cono, fprec, jto, km, kbb     &
                                     , kbm, kbot, weto, ielimi, gila, giph  &
                                     , api, dlxp, dlyp, zo, ddpo, lyears    &
                                     , lmonts, ldays, lyear1, lmont1        &
                                     , istart_new_topo_update               &
                                     , istart_new_run, istart_new_topo      &
                                     , tho, sictho, sicomo, tmelt, sicsno   &
                                     , io_stdout, sao, uoo, voe, dvo, avo   &
                                     , wo, z1o, sicuo, sicve, hibzeto       &
                                     , hibeto, hibzete, hibete, tice, kep   &
                                     , mmccdt, rhoo, fu10, swsum, swrab     &
                                     , txo, tye, rhoref_water, stabio       &
                                     , lwith_barotropic_stokes_drift        &
                                     , uaccel, vaccel, lbounds_exch_tp      &
                                     , have_g_js, po, eminpo, rivrun
  USE messy_mpiom_mem_e5,        ONLY: ltstranspose, ltswrite               &
                                     , spzndamp_time                        &
                                     , setunits, p_deco, contro, bounds_exch &
                                     , alloc_mem_bounds_exch_halo           &
                                     , set_param1, alloc_mem_commo1         &
                                     , alloc_mem_octdiff, alloc_mem_forcing &
                                     , alloc_mem_commoau2, alloc_mem_dilcor &
                                     , alloc_mem_commoau3, alloc_mem_elicom &
                                     , alloc_mem_para2, init_levitus_2d     &
                                     , init_levitus_3d, ocice, octher       &
                                     , setup_ocean_vertical_mixing          &
                                     , alloc_mem_stokes_drift, alloc_mem_adpo &
                                     , alloc_mem_commobbl, alloc_mem_gmbolus &
                                     , setup_grid, boden, coriol, findalfa  &
                                     , itprep, trotest2, glac_calv_ini      &
                                     , alloc_mem_swr_absorption             &
                                     , alloc_mem_tidal, foreph_ini          &
                                     , river_runoff_stations_ini            &
                                     , levitus_read_3d_stratification       &
                                     , levitus_horizontal_stratification    &
                                     , levitus_read_3d_restore, calc_rinum  &
                                     , levitus_per_month_setup, calc_dens   &
                                     , tipouv, itprep, troneu, ocvtro       &
                                     , calc_icecutoff, cell_thickness       &
                                     , slopetrans, ocadpo_base, ocadpo_trf  &
                                     , update_zo, octdiff_base, octdiff_trf &
                                     , relax_ts, ocvisc, correct_zo         &
                                     , river_runoff_omip_ini, foreph        &
                                     , convection, rhoicwa, rhosnwa         &
                                     ! mz_bk_20110315+
                                     , calc_global_mean, global_sum_salinity &
                                     , global_sum_temperature, global_mass  &
                                     , global_volume, global_salt_content   &
                                     ! mz_bk_20110315-
                                     ! op_mk_20180108+
                                     , dduo, ddue, dlxu, dlyu, dlxv, dlyv,  &
                                     , uko, vke
                                     ! _20180108-

  ! mz_bk_20110220-
#endif

  IMPLICIT NONE

  PRIVATE


  ! TIME STEP
  INTEGER, PUBLIC :: DT_steps
  INTEGER         :: DT_default

  ! NEEDED FOR TIME LOOP
  LOGICAL, SAVE ::  l_trig_mpiom_year = .TRUE.
  TYPE(io_time_event), PUBLIC, SAVE :: trig_mpiom_year
  TYPE(time_event),    PUBLIC, SAVE :: ev_trig_mpiom_year
  LOGICAL, SAVE ::  l_trig_mpiom_month = .TRUE.
  TYPE(io_time_event), PUBLIC, SAVE :: trig_mpiom_month
  TYPE(time_event),    PUBLIC, SAVE :: ev_trig_mpiom_month
  LOGICAL, SAVE ::  l_trig_mpiom_day = .TRUE.
  TYPE(io_time_event), PUBLIC, SAVE :: trig_mpiom_day
  TYPE(time_event),    PUBLIC, SAVE :: ev_trig_mpiom_day
  LOGICAL, SAVE, PUBLIC ::  l_trig_mpiom = .TRUE.
  TYPE(io_time_event), PUBLIC, SAVE :: trig_mpiom
  TYPE(time_event),    PUBLIC, SAVE :: ev_trig_mpiom


  ! OUTPUT VARIABLES

  ! mz_bk_20110315+
#if defined (MPIOM_2000)
  ! diagnostics
  REAL(dp)                  , POINTER :: out_global_sum_temperature => NULL()
  REAL(dp)                  , POINTER :: out_global_sum_salinity => NULL()
  REAL(dp)                  , POINTER :: out_global_mass => NULL()
  REAL(dp)                  , POINTER :: out_global_volume => NULL()
  REAL(dp)                  , POINTER :: out_global_salt_content => NULL()

  LOGICAL :: ldiag
#endif
  ! mz_bk_20110315-

  ! GRID RELATED
  REAL(dp), DIMENSION(:,:,:), POINTER :: decomp_mpiom  => NULL()
  REAL(dp), DIMENSION(:,:)  , POINTER :: lat_mpiom  => NULL()
  REAL(dp), DIMENSION(:,:)  , POINTER :: lon_mpiom  => NULL()
  REAL(dp), DIMENSION(:,:)  , POINTER :: out_dlxp  => NULL()
  REAL(dp), DIMENSION(:,:)  , POINTER :: out_dlyp  => NULL()
  ! op_mk_20180108+
  REAL(dp), DIMENSION(:,:)  , POINTER :: out_dlxu  => NULL()
  REAL(dp), DIMENSION(:,:)  , POINTER :: out_dlyu  => NULL()
  REAL(dp), DIMENSION(:,:)  , POINTER :: out_dlxv  => NULL()
  REAL(dp), DIMENSION(:,:)  , POINTER :: out_dlyv  => NULL()
  ! op_mk_20180108-

  ! 3D-variables
  REAL(dp), DIMENSION(:,:,:), POINTER :: OUT_RHOO => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: OUT_THO  => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: OUT_THO_K  => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: OUT_SAO  => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: OUT_WO   => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: OUT_PO   => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: OUT_VOE  => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: OUT_UOO  => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: OUT_DVO  => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: OUT_AVO  => NULL()
  ! op_mk_20180108+
  REAL(dp), DIMENSION(:,:,:), POINTER :: OUT_UKO  => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: OUT_VKE  => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: OUT_STABIO => NULL()
  ! op_mk_20180108-
  !PUBLIC 3D
  REAL(dp), DIMENSION(:,:,:), PUBLIC, POINTER :: OMMASS => NULL()
  REAL(dp), DIMENSION(:,:,:), PUBLIC, POINTER :: OMVOL => NULL()
  !2D-variables
  REAL(dp), DIMENSION(:,:), POINTER    :: OUT_ZO => NULL()
  REAL(dp), DIMENSION(:,:), POINTER    :: OUT_Z1O=> NULL()
  ! mz_bk_20111209+
  REAL(dp), DIMENSION(:,:,:), POINTER  :: OUT_DDPO=> NULL()
  ! mz_bk_20111209-
  ! op_mk_20180108+
  REAL(dp), DIMENSION(:,:,:), POINTER  :: OUT_DDUO=> NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER  :: OUT_DDUE=> NULL()
  ! op_mk_20180108-

  REAL(dp), DIMENSION(:,:), POINTER    :: OUT_SICTHO  => NULL()
  REAL(dp), DIMENSION(:,:), POINTER    :: OUT_TICE    => NULL()
  REAL(dp), DIMENSION(:,:), POINTER    :: OUT_SICOMO  => NULL()
  REAL(dp), DIMENSION(:,:), POINTER    :: OUT_SICUO   => NULL()
  REAL(dp), DIMENSION(:,:), POINTER    :: OUT_SICVE   => NULL()
  REAL(dp), DIMENSION(:,:), POINTER    :: OUT_HIBZETO => NULL()
  REAL(dp), DIMENSION(:,:), POINTER    :: OUT_HIBETO  => NULL()
  REAL(dp), DIMENSION(:,:), POINTER    :: OUT_HIBZETE => NULL()
  REAL(dp), DIMENSION(:,:), POINTER    :: OUT_HIBETE  => NULL()
  REAL(dp), DIMENSION(:,:), POINTER    :: OUT_SICSNO  => NULL()

  REAL(dp), DIMENSION(:,:), POINTER    :: OUT_SOCU  => NULL()
  REAL(dp), DIMENSION(:,:), POINTER    :: OUT_SOCV  => NULL()

  REAL(dp),                 POINTER    :: OUT_mmccdt   => NULL() !
  REAL(dp), DIMENSION(:,:), POINTER    :: OUT_RIVRUN   => NULL() !
  REAL(dp), DIMENSION(:,:), POINTER    :: OUT_EMINPO   => NULL() !

  ! COUPLING !!!!!!!!
  REAL(dp),DIMENSION(:,:),POINTER :: AOFLTXWO(:,:) !Zonal wind stress on water
  REAL(dp),DIMENSION(:,:),POINTER :: AOFLTYWO(:,:) !Meridional wind stress on water
  REAL(dp),DIMENSION(:,:),POINTER :: AOFLTXIO(:,:) !Zonal wind stress on snow/ice
  REAL(dp),DIMENSION(:,:),POINTER :: AOFLTYIO(:,:) !Meridional wind stress on snow/ice
  REAL(dp),DIMENSION(:,:),POINTER :: AOFLFRIO(:,:) !solid freshwater flux (over ice only)
  REAL(dp),DIMENSION(:,:),POINTER :: AOFLFRWO(:,:) !liquid freshwater flux (over water and ice)
  REAL(dp),DIMENSION(:,:),POINTER :: AOFLRHIO(:,:) !solid freshwater flux (over ice only)
  REAL(dp),DIMENSION(:,:),POINTER :: AOFLCHIO(:,:) !conductive heat flux through ice
  REAL(dp),DIMENSION(:,:),POINTER :: AOFLNHWO(:,:) !net heat flux over water
  REAL(dp),DIMENSION(:,:),POINTER :: AOFLSHWO(:,:) !downwelling solar radiation
  REAL(dp),DIMENSION(:,:),POINTER :: AOFLWSVO(:,:) !wind stress velocity

  REAL(dp),DIMENSION(:,:,:),POINTER :: abs_oce(:,:,:) !ocean absorption
  REAL(dp) :: opendep

!  ! HD-related
!  LOGICAL :: L_HD

  TYPE, PUBLIC :: TYP_LINK
     CHARACTER(LEN=STRLEN_CHANNEL) :: channel  = ''
     CHARACTER(LEN=STRLEN_OBJECT ) :: object   = ''
  END TYPE TYP_LINK

  ! GRID MASK
  REAL(dp), DIMENSION (:,:,:), POINTER :: mask  => NULL()

  ! location of initialization file
  CHARACTER(LEN=200) :: mpiom_init_file_name  = ''
  REAL(DP), ALLOCATABLE :: input_data3d(:,:,:)
  REAL(DP), ALLOCATABLE :: input_data2d(:,:)

  ! PUBLIC ECHAM-5 INTERFACE ROUTINES TO 'messy_SUBMODELS'
  PUBLIC :: mpiom_initialize
  PUBLIC :: mpiom_init_memory
  PUBLIC :: mpiom_init_coupling
  PUBLIC :: mpiom_global_start
  PUBLIC :: mpiom_global_end
  PUBLIC :: mpiom_free_memory

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      GENERAL DESCRIPTION OF MPIOM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!                       MP MPIOM
!                      SBR BELEG
!                          BODEN
!                          CORIOL
!                          ITPREP
!                            !
!                 !  --->  THERMODYNAMIC FORCING
!                 !        WIND FORCING
!      TIME       !        DECOMPOSITION INTO BAROTROPIC AND BAROCLINIC FIELD
!      STEPPING   !        BAROTROPIC SYSTEM
!                 !        BAROCLINIC SYSTEM
!                 !        MOMENTUM ADVECTION
!                 !        TRACER ADVECTION
!                 !        TRACER EDDY DIFFUSION
!                 !        TRACER DIFFUSION
!                 !        MOMENTUM DIFFUSION
!                 ! <---     !
!                            !
!                          OUTPUT-ROUTINES

!          ------------------------------
! NUMBER OF VERTICAL LAYERS IS KE.
!
!************************************************************
! PARAMETER :  IE   NUMBER OF GRID POINTS IN X
!              JE                            Y
!              KE                            Z
!              KBB=NMAX,ILL=MATR MUST BE SET BY USER!!!
!
! SOME IMPORTANT VARIABLES AND FIELDS :
!             DT        TIME STEP
!             TIESTU    DEPTH OF HORIZONTAL VELOCITY POINTS
!             TIESTW    DEPTH OF VERTICAL VELOCITY POINTS
!             ZO        SEA SURFACE ELEVATION
!             AMSUE/O   LAND/SEA-MASK FOR VECTORFIELD (LAND=0/SEA=1)
!             WETO                    FOR SCALARFIELD       "
!             UKO       ZONAL VELOCITY COMPONENT
!             VKE       MERIDIONAL "       "
!             WO        VERTICAL VELOCITY COMPONENT
!             THO       TEMPERATURE
!             SAO       SALINITY
!             PO        PRESSURE
!             TXO       ZONAL WIND STRESS
!             TYE       MERIDIONAL WIND STRESS
!             AVO       VARIABLE VERTICAL EDDY VISCOSITY
!                       (DEFINED ON SCALAR POINTS!!)
!             DVO       VARIABLE VERTICAL EDDY DIFFUSIVITY
!
!
! E=MERIDIONAL VELOCITY POINTS, O=ZONAL VELOCITY POINTS/SCALAR POINTS
!***********************************************************************

CONTAINS

  ! ---------------------------------------------------------------------------

  SUBROUTINE mpiom_initialize

       ! MPIOM MODULE ROUTINE (ECHAM-5 INTERFACE)
       !
       ! INITIALIZATION OF GLOBAL VARIABLES FROM NAMELIST
       ! IN PARALLEL ENVIRONMENT
       !
       ! Author: Pozzer Andrea, MPICH, May 2007
       ! ECHAM5/MESSy
#ifdef MPIOM_13B
       USE messy_main_mpi_bi,            ONLY: p_parallel_io, p_io, p_bcast !!$, finish
#elif defined (MPIOM_2000)
       USE messy_main_mpi_bi,            ONLY: p_parallel_io, p_io, p_pe, &
                                               p_bcast,                   &
                                               p_abort, init_mpi_datatypes
#endif
       USE messy_main_tools,             ONLY: find_next_free_unit
       USE messy_main_grid_def_mem_bi,   ONLY: jrow, nlev, nlevp1
       USE messy_main_timer,             ONLY: delta_time, time_step_len
       USE messy_main_channel_error_bi,  ONLY: channel_halt
       USE messy_main_channel_bi,        ONLY: GP_2D_HORIZONTAL, DC_GP_MPIOM   &
                                             , GP_3D_MPIOM, GP_3D_MPIOM_INT    &
                                             , GP_2D_MPIOM, SCALAR             &
                                             , DIMID_LEV_MPIOM                 &
                                             , DIMID_LEV_MPIOM_INT             &
                                             , DIMID_DEPTH_MPIOM               &
                                             , DIMID_DEPTH_MPIOM_INT &
                                             , DIMID_LON_MPIOM, DIMID_LAT_MPIOM
       USE messy_main_channel_repr,       ONLY: new_representation, AUTO
       USE messy_main_channel_dimensions, ONLY: new_dimension,    &
                                                add_dimension_variable_att, &
                                                add_dimension_variable


       IMPLICIT NONE

       INTRINSIC INT,FLOOR

       ! DEFINE NCREGRID EVENT TRIGGERs
       ! LOCAL
       CHARACTER(LEN=*), PARAMETER :: substr='mpiom_initialize'
       INTEGER                 :: status
       INTEGER                 :: iou    ! I/O unit
       REAL(DP), ALLOCATABLE, DIMENSION(:) :: array
       REAL(DP), ALLOCATABLE, DIMENSION(:) :: depthw

       CALL start_message_bi(modstr, 'INITIALIZATION', substr)

       ! set input-output unit

       CALL SETUNITS

!       IF (io_stdout > 0) THEN
!         CALL OPEN_STDOUT(io_stdout,'oceout')
!       ENDIF

!-----------------------------------------------------------------------
!                 CTRL NAMELIST
!-----------------------------------------------------------------------

       ! INITIALIZE CTRL
       IF (p_parallel_io) THEN
          iou = find_next_free_unit(100,200)
          ! *** CALL CORE ROUTINE:
          CALL mpiom_read_nml_ctrl(status, iou)
!!$       IF (status /= 0) CALL finish(substr)     ! op_pj_20161104
          IF (status /= 0) CALL error_bi('',substr)   ! op_pj_20161104
       END IF


     ! BROADCAST RESULTS
     ! GLOBAL SETTINGS
     CALL p_bcast(GRID              , p_io)
     CALL p_bcast(VLEVELS           , p_io)
     CALL p_bcast(delta_time_mpiom  , p_io)
     CALL p_bcast(NPROCA            , p_io)
     CALL p_bcast(NPROCB            , p_io)
     CALL p_bcast(L_TIDES           , p_io)
     CALL p_bcast(L_COUPLING        , p_io)
     CALL p_bcast(L_HAMOCC_COUPLING , p_io)
     CALL p_bcast(L_CHECK           , p_io)

!-----------------------------------------------------------------------
!                 GRID AND TIME STEP
!-----------------------------------------------------------------------

     SELECT CASE (GRID)
        CASE ("GR60")
           IE_G = 60
           JE_G = 50
           DT_default   = 10800
        CASE ("GR30")
           IE_G = 122
           JE_G = 101
           DT_default   = 8640
        CASE ("GR15")
           IE_G = 256
           JE_G = 220
           DT_default   = 4800
        CASE ("TP04")
           IE_G = 802
           JE_G = 401
           DT_default   = 3600
        CASE ("TP10")
           IE_G = 362
           JE_G = 192
           DT_default   = 4800
        CASE ("TP40")
           IE_G = 82
           JE_G = 41
           DT_default   = 8640
        CASE DEFAULT
           WRITE(*,*) "MPIOM GRID NOT FROM DEFAULT!"
           WRITE(*,*) "GRID: ",GRID
           IE_G = IE_G_nml
           JE_G = JE_G_nml
           !CALL finish(substr)
     END SELECT

     SELECT CASE (VLEVELS)
        CASE("L20")
           cdzw(1:20)=(/20.,20., 20., 30.,40.,50.,70.        &
                ,90.,120.,150.,180.,210.,250.,300.           &
                ,400.,500.,600.,700.,900.,1400./)
           KE=20
        CASE("L40")
           cdzw(1:40)=(/12.,10.,10.,10.,10.,10.,13.,15.,20.,25.,    &
                 30.,35.,40.,45.,50.,55.,60.,70.,80.,90.,           &
                 100.,110.,120.,130.,140.,150.,170.,180.,190.,200., &
                 220.,250.,270.,300.,350.,400.,450.,500.,500.,600./)
           KE=40
        CASE("L3")
           cdzw(1:3)=(/12.,10.,5000./)
           KE=3
        CASE DEFAULT
! op_pj_20161104+
!!$        WRITE(*,*) "WRONG SLECTION OF MPIOM VERTICAL LEVELS!"
!!$        CALL finish(substr) ! op_pj_20161104
           CALL error_bi('WRONG SLECTION OF MPIOM VERTICAL LEVELS!',substr) ! op_pj_20161104
! op_pj_20161104-
     END SELECT

     ! TIME STEP definition

      IF (delta_time_mpiom == 0 ) THEN
        DT = DT_default
      ELSE
        DT = delta_time_mpiom
      ENDIF

      DT_steps=FLOOR(DT/delta_time)
      DT = DT_steps*delta_time

      IF (p_pe==p_io) THEN
         IF (L_TIDES) WRITE(*,*) "Tides calculation = ON"
         WRITE (*,*) "DT_steps = ", DT_steps
         WRITE (*,*) "DT = ", DT
      ENDIF


      CALL p_bcast(IE_G        , p_io)
      CALL p_bcast(JE_G        , p_io)
      CALL p_bcast(DT          , p_io)
      CALL p_bcast(cdzw        , p_io)

!-----------------------------------------------------------------------
!                 OCECTL NAMELIST
!-----------------------------------------------------------------------

       IF (p_parallel_io) THEN
          iou = find_next_free_unit(100,200)
          CALL mpiom_read_nml_ocectl(status, iou)
! op_pj_20161104+
!!$       IF (status /= 0) CALL finish(substr)
          IF (status /= 0) CALL error_bi('',substr)
! op_pj_20161104+
       END IF

#if defined (MPIOM_13B)
      CALL p_bcast(CAULAPTS    , p_io)
      CALL p_bcast(CAULAPUV    , p_io)
      CALL p_bcast(CAH00       , p_io)
      CALL p_bcast(AUS         , p_io)
      CALL p_bcast(AV0         , p_io)
      CALL p_bcast(DV0         , p_io)
      CALL p_bcast(CWT         , p_io)
      CALL p_bcast(CSTABEPS    , p_io)
      CALL p_bcast(DBACK       , p_io)
      CALL p_bcast(ABACK       , p_io)
      CALL p_bcast(CRELSAL     , p_io)
      CALL p_bcast(CRELTEM     , p_io)
      CALL p_bcast(ICONVA      , p_io)
      CALL p_bcast(IWIRB       , p_io)
      CALL p_bcast(CDVOCON     , p_io)
      CALL p_bcast(CAVOCON     , p_io)
      CALL p_bcast(IOCAD       , p_io)
      CALL p_bcast(isnflg      , p_io)
      CALL p_bcast(H0          , p_io)
      CALL p_bcast(HMIN        , p_io)
      CALL p_bcast(ARMIN       , p_io)
      CALL p_bcast(ARMAX       , p_io)
      CALL p_bcast(HSNTOICE    , p_io)
      CALL p_bcast(SICTHMIN    , p_io)
      CALL p_bcast(SICE        , p_io)
      CALL p_bcast(D3          , p_io)
      CALL p_bcast(ISTART      , p_io)
      CALL p_bcast(I3DREST     , p_io)
      CALL p_bcast(LY_END      , p_io)
      CALL p_bcast(LY_START    , p_io)
      CALL p_bcast(LM_START    , p_io)
      CALL p_bcast(ICONTRO     , p_io)

      IF (p_parallel_io) THEN
         WRITE(*,*) " IE_G     = ", IE_G
         WRITE(*,*) " JE_G     = ", JE_G
         WRITE(*,*) " KE       = ", KE
         WRITE(*,*) " CAULAPTS = ", CAULAPTS
         WRITE(*,*) " CAULAPUV = ", CAULAPUV
         WRITE(*,*) " CAH00    = ", CAH00
         WRITE(*,*) " AUS      = ", AUS
         WRITE(*,*) " AV0      = ", AV0
         WRITE(*,*) " DV0      = ", DV0
         WRITE(*,*) " CWT      = ", CWT
         WRITE(*,*) " CSTABEPS = ", CSTABEPS
         WRITE(*,*) " DBACK    = ", DBACK
         WRITE(*,*) " ABACK    = ", ABACK
         WRITE(*,*) " CRELSAL  = ", CRELSAL
         WRITE(*,*) " CRELTEM  = ", CRELTEM
         WRITE(*,*) " ICONVA   = ", ICONVA
         WRITE(*,*) " IWIRB    = ", IWIRB
         WRITE(*,*) " CDVOCON  = ", CDVOCON
         WRITE(*,*) " CAVOCON  = ", CAVOCON
         WRITE(*,*) " IOCAD    = ", IOCAD
         WRITE(*,*) " isnflg   = ", isnflg
         WRITE(*,*) " H0       = ", H0
         WRITE(*,*) " HMIN     = ", HMIN
         WRITE(*,*) " ARMIN    = ", ARMIN
         WRITE(*,*) " ARMAX    = ", ARMAX
         WRITE(*,*) " HSNTOICE = ", HSNTOICE
         WRITE(*,*) " SICTHMIN = ", SICTHMIN
         WRITE(*,*) " SICE     = ", SICE
         WRITE(*,*) " D3       = ", D3
         WRITE(*,*) " ISTART   = ", ISTART
         WRITE(*,*) " I3DREST  = ", I3DREST
         WRITE(*,*) " LY_END   = ", LY_END
         WRITE(*,*) " LY_START = ", LY_START
         WRITE(*,*) " LM_START = ", LM_START
         WRITE(*,*) " ICONTRO  = ", ICONTRO
      ENDIF

#elif defined (MPIOM_2000)
    IF (p_parallel_io) THEN
    write(*,*) 'iocad=',iocad
      IF (imean .EQ. 0) THEN
        IF ( LGMDIAG ) WRITE(*,*)'imean = 0 => LGMDIAG disabled'
        LGMDIAG = .FALSE.
        IF ( LFORCEDIAG ) WRITE(*,*)'imean = 0 => LFORCEDIAG disabled'
        LFORCEDIAG = .FALSE.
        IF ( LHFLDIAG ) WRITE(*,*)'imean = 0 => LHFLDIAG disabled'
        LHFLDIAG = .FALSE.
        !todo check LHFLDIAG and then remove
#ifdef __coupled
        IF ( LHFLDIAG ) WRITE(*,*)'coupled => LHFLDIAG disabled'
        LHFLDIAG = .FALSE.
#endif
        IF ( LCONVDIAG ) WRITE(*,*)'imean = 0 => LCONVDIAG disabled'
        LCONVDIAG = .FALSE.
        IF ( LDIFFDIAG ) WRITE(*,*)'imean = 0 => LDIFFDIAG disabled'
        LDIFFDIAG = .FALSE.
        IF ( LGRIDINFO ) WRITE(*,*)'imean = 0 => LGRIDINFO disabled'
        LGRIDINFO = .FALSE.
        IF ( LCALCDIFI ) WRITE(*,*)'imean = 0 => LCALCDIFI disabled'
        LCALCDIFI = .FALSE.
      END IF

      IF ( IBOLK .EQ. 0 )  LGMDIAG = .FALSE.
      IF ( IBOLK .EQ. 0 ) WRITE(*,*)'ibolk = 0 => GM scheme disabled'
      IF ( IBOLK .GT. 0 ) WRITE(*,*)'ibolk > 0 => GM scheme enabled'

      IF ( .NOT. LISOPYC ) THEN
        IF ( IBOLK .GT. 0 ) WRITE(*,*)                                  &
             'ATTN => GM scheme enabled with horizontal ',              &
             ' instead of isopycnal diffusion'
      ENDIF
      IF ( IBOLK .LT. 0 ) WRITE(*,*)                                  &
           'ibolk < 0 => GM + Visbek scheme enabled'
      IF (.NOT. LISOPYC) THEN
        IF ( IBOLK .LT. 0 ) WRITE(*,*)                                  &
             'ATTN => GM+Visbek scheme enabled with horizontal ',       &
             ' instead of isopycnal diffusion'
      ENDIF

      IF (iocad == 4) THEN
        WRITE(*,'(/,a,/,a,/)') ' iocad = 4 is obsolete => use iocad = 3 '// &
             '(TVD advection scheme) ', '  and ibbl_transport = 1 '//  &
             '(BBL transport to neutral density level) instead'
        iocad = 3
        ibbl_transport = 1
      END IF

      SELECT CASE (iocad)
      CASE (1:3, 8)
        ladpo = .TRUE.    ! use tracer advection routines from module mo_adpo
        ladfs = .FALSE.
      CASE (5:7)
        ladpo = .FALSE.
        ladfs = .TRUE.    ! use tracer advection routine ocadfs
        IF (ibbl_transport == 1) THEN
          WRITE (*,*) 'iocad =', iocad, '=> BBL transport scheme disabled'
          ibbl_transport = 0
        END IF
      CASE DEFAULT
        WRITE (*,*) 'ERROR: The specified iocad value', iocad, 'is not valid!'
        CALL p_abort
      END SELECT

      IF (iocaduv /= 3 .AND. iocaduv /= 8) THEN
        WRITE (*,*) 'ERROR: The specified iocaduv value', iocaduv, &
             'is not valid!'
        CALL p_abort
      ENDIF

    ENDIF ! (p_pe==p_io)

    CALL p_bcast(CAULAPTS,p_io)
    CALL p_bcast(CAULAPUV,p_io)
    CALL p_bcast(CAH00,p_io)

    CALL p_bcast(nfixYearLen,p_io)
    CALL p_bcast(AUS,p_io)
    CALL p_bcast(bofric,p_io)
    CALL p_bcast(rayfric,p_io)

    CALL p_bcast(AV0,p_io)
    CALL p_bcast(DV0,p_io)
    CALL p_bcast(CWT,p_io)

    CALL p_bcast(CWA,p_io)
    CALL p_bcast(CSTABEPS,p_io)
    CALL p_bcast(DBACK,p_io)
    CALL p_bcast(ABACK,p_io)
    CALL p_bcast(CRELSAL,p_io)
    CALL p_bcast(CRELTEM,p_io)

    CALL p_bcast(NUMRIV,p_io)
    CALL p_bcast(NUMGLAC,p_io)

    CALL p_bcast(IAUFR,p_io)
    CALL p_bcast(IAUFW,p_io)
    CALL p_bcast(NYEARS,p_io)
    CALL p_bcast(NMONTS,p_io)
    CALL p_bcast(NDAYS,p_io)
    CALL p_bcast(IMEAN,p_io)
    CALL p_bcast(ISTART,p_io)
    CALL p_bcast(I3DREST,p_io)
    CALL p_bcast(IOASISFLUX,p_io)

    CALL p_bcast(RRELSAL,p_io)
    CALL p_bcast(RRELTEM,p_io)
    CALL p_bcast(CDVOCON,p_io)
    CALL p_bcast(CAVOCON,p_io)

    CALL p_bcast(rleadclose,p_io)
    CALL p_bcast(H0,p_io)
    CALL p_bcast(HMIN,p_io)
    CALL p_bcast(ARMIN,p_io)
    CALL p_bcast(ARMAX,p_io)
    CALL p_bcast(HSNTOICE,p_io)
    CALL p_bcast(SICTHMIN,p_io)
    CALL p_bcast(SICE,p_io)

    CALL p_bcast(ibbl_transport,p_io)

    CALL p_bcast(IBOLK,p_io)
    CALL p_bcast(ISNFLG,p_io)
    CALL p_bcast(IOCAD,p_io)
    CALL p_bcast(IOCADUV,p_io)
    CALL p_bcast(IOCONV,p_io)
    CALL p_bcast(ICONTRO,p_io)

!    CALL p_bcast(ihalo,p_io)
    CALL p_bcast(imod,p_io)

    CALL p_bcast(LISOPYC,p_io)
    CALL p_bcast(lmpitype,p_io)
    CALL p_bcast(lnonblock,p_io)
    CALL p_bcast(ltidal,p_io)
    CALL p_bcast(ltidal_diag,p_io)

    CALL p_bcast(lswr_jerlov,p_io)
    CALL p_bcast(lfb_bgc_oce,p_io)
    CALL p_bcast(jerlov_atten,p_io)
    CALL p_bcast(jerlov_bluefrac,p_io)

    CALL p_bcast(lwith_one_layer_shelfs,p_io)

    CALL p_bcast(LFORCEDIAG,p_io)
    CALL p_bcast(LHFLDIAG,p_io)
    CALL p_bcast(LCONVDIAG,p_io)
    CALL p_bcast(LDIFFDIAG,p_io)
    CALL p_bcast(LGMDIAG,p_io)
    CALL p_bcast(LGRIDINFO,p_io)
    CALL p_bcast(LCALCDIFI,p_io)

    CALL p_bcast(ladpo,p_io)
    CALL p_bcast(ladfs,p_io)
    CALL p_bcast(lsaoclose,p_io)
    CALL p_bcast(lundelayed_momentum_advection,p_io)
    CALL p_bcast(fp_tracing_enabled, p_io)


    CALL p_bcast(iter_sor,p_io)
    CALL p_bcast(iter_sor_hack,p_io)
    CALL p_bcast(rtsorpar,p_io)
    CALL p_bcast(rtsorpar_hack,p_io)

    ! mz_bk_20110317+
    CALL p_bcast(ldiag,p_io)
    ! mz_bk_20110317-

    IF (p_parallel_io) THEN
      WRITE(*,*) "CAULAPTS               =", CAULAPTS
      WRITE(*,*) "CAULAPUV               =", CAULAPUV
      WRITE(*,*) "CAH00                  =", CAH00
      WRITE(*,*) "nfixYearLen            =", nfixYearLen
      WRITE(*,*) "AUS                    =", AUS
      WRITE(*,*) "bofric                 =", bofric
      WRITE(*,*) "rayfric                =", rayfric
      WRITE(*,*) "AV0                    =", AV0
      WRITE(*,*) "DV0                    =", DV0
      WRITE(*,*) "CWT                    =", CWT
      WRITE(*,*) "CWA                    =", CWA
      WRITE(*,*) "CSTABEPS               =", CSTABEPS
      WRITE(*,*) "DBACK                  =", DBACK
      WRITE(*,*) "ABACK                  =", ABACK
      WRITE(*,*) "CRELSAL                =", CRELSAL
      WRITE(*,*) "CRELTEM                =", CRELTEM
      WRITE(*,*) "NUMRIV                 =", NUMRIV
      WRITE(*,*) "NUMGLAC                =", NUMGLAC
      WRITE(*,*) "IAUFR                  =", IAUFR
      WRITE(*,*) "IAUFW                  =", IAUFW
      WRITE(*,*) "NYEARS                 =", NYEARS
      WRITE(*,*) "NMONTS                 =", NMONTS
      WRITE(*,*) "NDAYS                  =", NDAYS
      WRITE(*,*) "IMEAN                  =", IMEAN
      WRITE(*,*) "ISTART                 =", ISTART
      WRITE(*,*) "I3DREST                =", I3DREST
      WRITE(*,*) "IOASISFLUX             =", IOASISFLUX
      WRITE(*,*) "RRELSAL                =", RRELSAL
      WRITE(*,*) "RRELTEM                =", RRELTEM
      WRITE(*,*) "CDVOCON                =", CDVOCON
      WRITE(*,*) "CAVOCON                =", CAVOCON
      WRITE(*,*) "rleadclose             =", rleadclose
      WRITE(*,*) "H0                     =", H0
      WRITE(*,*) "HMIN                   =", HMIN
      WRITE(*,*) "ARMIN                  =", ARMIN
      WRITE(*,*) "ARMAX                  =", ARMAX
      WRITE(*,*) "HSNTOICE               =", HSNTOICE
      WRITE(*,*) "SICTHMIN               =", SICTHMIN
      WRITE(*,*) "SICE                   =", SICE
      WRITE(*,*) "ibbl_transport         =", ibbl_transport
      WRITE(*,*) "IBOLK                  =", IBOLK
      WRITE(*,*) "ISNFLG                 =", ISNFLG
      WRITE(*,*) "IOCAD                  =", IOCAD
      WRITE(*,*) "IOCADUV                =", IOCADUV
      WRITE(*,*) "IOCONV                 =", IOCONV
      WRITE(*,*) "ICONTRO                =", ICONTRO
!      WRITE(*,*) "ihalo                  =", ihalo
      WRITE(*,*) "imod                   =", imod
      WRITE(*,*) "LISOPYC                =", LISOPYC
      WRITE(*,*) "lmpitype               =", lmpitype
      WRITE(*,*) "lnonblock              =", lnonblock
      WRITE(*,*) "ltidal                 =", ltidal
      WRITE(*,*) "ltidal_diag            =", ltidal_diag
      WRITE(*,*) "lswr_jerlov            =", lswr_jerlov
      WRITE(*,*) "lfb_bgc_oce            =", lfb_bgc_oce
      WRITE(*,*) "jerlov_atten           =", jerlov_atten
      WRITE(*,*) "jerlov_bluefrac        =", jerlov_bluefrac
      WRITE(*,*) "lwith_one_layer_shelfs =", lwith_one_layer_shelfs
      WRITE(*,*) "LFORCEDIAG             =", LFORCEDIAG
      WRITE(*,*) "LHFLDIAG               =", LHFLDIAG
      WRITE(*,*) "LCONVDIAG              =", LCONVDIAG
      WRITE(*,*) "LDIFFDIAG              =", LDIFFDIAG
      WRITE(*,*) "LGMDIAG                =", LGMDIAG
      WRITE(*,*) "LGRIDINFO              =", LGRIDINFO
      WRITE(*,*) "LCALCDIFI              =", LCALCDIFI
      WRITE(*,*) "ladpo                  =", ladpo
      WRITE(*,*) "ladfs                  =", ladfs
      WRITE(*,*) "lsaoclose              =", lsaoclose
      WRITE(*,*) "fp_tracing_enabled     =", fp_tracing_enabled
      WRITE(*,*) "iter_sor               =", iter_sor
      WRITE(*,*) "iter_sor_hack          =", iter_sor_hack
      WRITE(*,*) "rtsorpar               =", rtsorpar
      WRITE(*,*) "rtsorpar_hack          =", rtsorpar_hack
      WRITE(*,*) "lundelayed_momentum_advection =", lundelayed_momentum_advection
    ENDIF
#endif

!-----------------------------------------------------------------------
!                 INIT NAMELIST
!-----------------------------------------------------------------------
    IF (ISTART.EQ.4) THEN
      IF (p_parallel_io) THEN
         iou = find_next_free_unit(100,200)
         CALL mpiom_read_nml_init(status, iou)
! op_pj_20161104+
!!$      IF (status /= 0) CALL finish(substr)
         IF (status /= 0) CALL error_bi('',substr)
! op_pj_20161104+
      END IF
    CALL p_bcast(mpiom_init_file_name,p_io)
    ENDIF


!-----------------------------------------------------------------------
!              SURFACE RELAXATION EQUAL TO ZERO IF COUPLED
!-----------------------------------------------------------------------
      IF (L_COUPLING) THEN
        CRELTEM = 0.0_dp !no surface relaxation
        CRELSAL = 0.0_dp !no surface relaxation
      ENDIF
!-----------------------------------------------------------------------
!                          EVENT INITIALIZATION
!-----------------------------------------------------------------------

  ! initialize mpiom event
!  trig_mpiom = io_time_event ( INT(DT_steps), 'steps',TRIG_FIRST, 0)
  trig_mpiom = io_time_event ( DT, TIME_INC_SECONDS,TRIG_FIRST, 0)
  trig_mpiom_day = io_time_event ( 1, TIME_INC_DAYS,TRIG_FIRST,0)
  trig_mpiom_month = io_time_event ( 1, TIME_INC_MONTHS,TRIG_FIRST,0)
  trig_mpiom_year = io_time_event ( 1, TIME_INC_YEARS,TRIG_FIRST,0)


  CALL p_bcast_event (trig_mpiom, p_io)
  CALL p_bcast_event (trig_mpiom_day, p_io)
  CALL p_bcast_event (trig_mpiom_month, p_io)
  CALL p_bcast_event (trig_mpiom_year, p_io)

  CALL timer_event_init(ev_trig_mpiom, trig_mpiom, &
       'mpiom computation', 'present')
  CALL timer_event_init(ev_trig_mpiom_day, trig_mpiom_day, &
       'mpiom computation_day', 'present')
  CALL timer_event_init(ev_trig_mpiom_month, trig_mpiom_month, &
       'mpiom computation_month', 'present')
  CALL timer_event_init(ev_trig_mpiom_year, trig_mpiom_year, &
       'mpiom computation_year', 'present')

  IF (p_pe==p_io) THEN
     WRITE (*,*) 'trig_mpiom: '      ,trig_mpiom
     WRITE (*,*) 'trig_mpiom_day: '  ,trig_mpiom_day
     WRITE (*,*) 'trig_mpiom_month: ',trig_mpiom_month
     WRITE (*,*) 'trig_mpiom_year: ' ,trig_mpiom_year
  ENDIF


!-----------------------------------------------------------------------
!                          DECOMPOSITION
!-----------------------------------------------------------------------

      nprocx=NPROCA
      nprocy=NPROCB

      IF (p_parallel_io) THEN
         WRITE(*,*) " ##### MPIOM DECOMPOSITION ########"
         WRITE(*,*) " nprocx    = ",nprocx
         WRITE(*,*) " nprocy    = ",nprocy
      ENDIF

      ! initialize decomposition
      call p_deco

#if defined (MPIOM_2000)
    CALL alloc_mem_bounds_exch_halo
    CALL init_MPI_datatypes(IE,JE,KE)
#endif

! um_ak_20130705+  moved to mpiom_init_memory
!!$! um_ak_20130620+ ! mz_ap_20130620+
!!$!-----------------------------------------------------------------------
!!$!                      GRID DEFINTION
!!$!-----------------------------------------------------------------------
!!$    ! define ocean grid for use in other model part e.g. IMPORT-GRID
!!$    CALL define_mpiom_grid
!!$! um_ak_20130620- ! mz_ap_20130620-
! um_ak_20130705-
!-----------------------------------------------------------------------
!                      NEW REPRESENTATION
!-----------------------------------------------------------------------

    ! (1b) NEW DIMENSIONS
    !-----------------_LEVELS-----------------------------------
    ALLOCATE(array(ke))
    DO i=1,ke
       array(i) = REAL(i, DP)
    END DO

    CALL new_dimension(status, DIMID_LEV_MPIOM, 'mpiom_level', ke)
    CALL channel_halt(substr, status)

    CALL add_dimension_variable(status, 'mpiom_level', 'mpiom_level',array )
    CALL channel_halt(substr, status)

    CALL add_dimension_variable_att(status, 'mpiom_level', 'mpiom_level', &
         'long_name', c='ocean mpiom levels ')
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'mpiom_level', 'mpiom_level', &
         'units', c='level')
    CALL channel_halt(substr, status)

    DEALLOCATE(array)
    !-----------------_INTERFACE LEVELS-----------------------------------
    ALLOCATE(array(ke+1))
    DO i=1,ke+1
       array(i) = REAL(i, DP)
    END DO

    CALL new_dimension(status, DIMID_LEV_MPIOM_INT, 'mpiom_interface_level', ke+1)
    CALL channel_halt(substr, status)

    CALL add_dimension_variable(status, 'mpiom_interface_level', 'mpiom_interface_level',array )
    CALL channel_halt(substr, status)

    CALL add_dimension_variable_att(status, 'mpiom_interface_level', 'mpiom_interface_level', &
         'long_name', c='ocean mpiom levels ')
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'mpiom_interface_level', 'mpiom_interface_level', &
         'units', c='level')
    CALL channel_halt(substr, status)

    DEALLOCATE(array)
    !-----------------_LEVELS DEPTHS-----------------------------------
    ALLOCATE(array(ke))
    ALLOCATE(depthw(ke+1))

    depthw(1) = 0
    DO i=1,KE
      depthw(i+1)=depthw(i)+cdzw(i)
      array(i)=0.5*(depthw(i+1)+depthw(i))
    ENDDO

    CALL new_dimension(status, DIMID_DEPTH_MPIOM, 'mpiom_depth', ke)
    CALL channel_halt(substr, status)

    CALL add_dimension_variable(status, 'mpiom_depth', 'mpiom_depth',array )
    CALL channel_halt(substr, status)

    CALL add_dimension_variable_att(status, 'mpiom_depth', 'mpiom_depth', &
         'long_name', c='ocean mpiom_depth ')
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'mpiom_depth', 'mpiom_depth', &
         'units', c='m')
    CALL channel_halt(substr, status)

    DEALLOCATE(array)
    DEALLOCATE(depthw)
    !-----------------_INTERFACE LEVELS DEPTHS-----------------------------------
    ALLOCATE(depthw(ke+1))

    depthw(1) = 0
    DO i=1,KE
      depthw(i+1)=depthw(i)+cdzw(i)
    ENDDO

    CALL new_dimension(status, DIMID_DEPTH_MPIOM_INT, 'mpiom_interface_depth', ke+1)
    CALL channel_halt(substr, status)

    CALL add_dimension_variable(status, 'mpiom_interface_depth', 'mpiom_interface_depth',depthw )
    CALL channel_halt(substr, status)

    CALL add_dimension_variable_att(status, 'mpiom_interface_depth', 'mpiom_interface_depth', &
         'long_name', c='ocean mpiom_interface_depth ')
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'mpiom_interface_depth', 'mpiom_interface_depth', &
         'units', c='m')
    CALL channel_halt(substr, status)

    DEALLOCATE(depthw)
    !-----------------_LON-----------------------------------

    ALLOCATE(array(ie_g))
    DO i=1, ie_g
       array(i) = REAL(i, DP)
    END DO

    CALL new_dimension(status, DIMID_LON_MPIOM, 'mpiom_lon', ie_g)
    CALL channel_halt(substr, status)

    CALL add_dimension_variable(status, 'mpiom_lon', 'mpiom_lon',array )
    CALL channel_halt(substr, status)

    CALL add_dimension_variable_att(status, 'mpiom_lon', 'mpiom_lon', &
         'long_name', c='ocean grid cell longitude index')
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'mpiom_lon', 'mpiom_lon', &
         'units', c='')
    CALL channel_halt(substr, status)

    DEALLOCATE(array)
    !-----------------_LAT-----------------------------------
    ALLOCATE(array(je_g))
    DO i=1, je_g
       array(i) = REAL(i, DP)
    END DO


    CALL new_dimension(status, DIMID_LAT_MPIOM, 'mpiom_lat', je_g)
    CALL channel_halt(substr, status)

    CALL add_dimension_variable(status, 'mpiom_lat', 'mpiom_lat',array )
    CALL channel_halt(substr, status)

    CALL add_dimension_variable_att(status, 'mpiom_lat', 'mpiom_lat', &
         'long_name', c='ocean grid cell latitude index')
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'mpiom_lat', 'mpiom_lat', &
         'units', c='')
    CALL channel_halt(substr, status)

    DEALLOCATE(array)

    !-------------------------------------------------

    ! (1c) NEW REPRESENTATION
    CALL new_representation(status, GP_3D_MPIOM, 'GP_3D_MPIOM' &
         , rank = 3, link = 'xxx-', dctype = DC_GP_MPIOM       &
         , dimension_ids = (/ DIMID_LON_MPIOM, DIMID_LAT_MPIOM, DIMID_DEPTH_MPIOM /) &
         , ldimlen       = (/ ie , je, AUTO   /) &
         , output_order  = (/ 1,2,3 /)                           &
         , axis = 'XYZ-'                                         &
         )
    CALL channel_halt(substr, status)

    ! (1c) NEW REPRESENTATION
    CALL new_representation(status, GP_3D_MPIOM_INT, 'GP_3D_MPIOM_INT' &
         , rank = 3, link = 'xxx-', dctype = DC_GP_MPIOM     &
         , dimension_ids = (/ DIMID_LON_MPIOM, DIMID_LAT_MPIOM, DIMID_DEPTH_MPIOM_INT /) &
         , ldimlen       = (/ ie , je, AUTO   /) &
         , output_order  = (/ 1,2,3 /)                           &
         , axis = 'XYZ-'                                         &
         )
    CALL channel_halt(substr, status)

    CALL new_representation(status, GP_2D_MPIOM, 'GP_2D_MPIOM' &
         , rank = 2, link = 'xx--', dctype = DC_GP_MPIOM       &
         , dimension_ids = (/ DIMID_LON_MPIOM, DIMID_LAT_MPIOM /) &
         , ldimlen       = (/ ie , je  /) &
         , output_order  = (/ 1,2 /)                           &
         , axis = 'XY--'                                       &
         )
    CALL channel_halt(substr, status)


  END SUBROUTINE mpiom_initialize

  ! ---------------------------------------------------------------------------

  SUBROUTINE mpiom_init_memory

    USE messy_main_mpi_bi,           ONLY: p_parallel_io, p_io, p_pe, p_bcast
    USE messy_main_channel_dimensions
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_2D_HORIZONTAL, DC_GP_MPIOM   &
                                         , GP_3D_MPIOM, GP_3D_MPIOM_INT    &
                                         , GP_2D_MPIOM, SCALAR
    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         , new_attribute
    USE messy_mpiom_tools_e5,        ONLY: nullify_borders

    IMPLICIT NONE

    ! CHANNEL MANAGEMENT
    INTEGER                 :: status
    INTEGER                 :: iou    ! I/O unit
    INTRINSIC INT, TRIM

    ! LOCAL
    CHARACTER(len=*), PARAMETER  :: substr = 'mpiom_init_memory'


    CALL start_message_bi(modstr, 'MEMORY INITIALIZATION', substr)


    ! (2) INITIALIZE MEMORY

    ! (2b) new channel
    CALL new_channel(status, modstr, reprid=GP_3D_MPIOM, lrestreq=.TRUE.)
    CALL channel_halt(substr,status)

! CHANNEL OBJECT

    ! (2c) new channel objects

!-------------------------------------------------------------

!--------------------   decomposition

    CALL new_channel_object(status, modstr, 'pe' &
         , p3 = decomp_mpiom, lstatic=.TRUE.)
    CALL channel_halt(substr, status)

!--------------------   latitude

    CALL new_channel_object(status, modstr, 'latitude' &
         , p2 = lat_mpiom, reprid = GP_2D_MPIOM, lstatic=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'latitude' &
                , 'units', c='degrees_north')
    CALL channel_halt(substr, status)

!--------------------   longitude

    CALL new_channel_object(status, modstr, 'longitude' &
         , p2 = lon_mpiom, reprid = GP_2D_MPIOM, lstatic=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'longitude' &
                , 'units', c='degrees_east')
    CALL channel_halt(substr, status)

!--------------------   x-lenght

    CALL new_channel_object(status, modstr, 'dlxp' &
         , p2 = out_dlxp, reprid = GP_2D_MPIOM, lstatic=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'dlxp' &
                , 'long_name', c='x length grid box' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dlxp' &
                , 'units', c='m')
    CALL channel_halt(substr, status)

!--------------------   y-lenght

    CALL new_channel_object(status, modstr, 'dlyp' &
         , p2 = out_dlyp, reprid = GP_2D_MPIOM, lstatic=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'dlyp' &
                , 'long_name', c='y length grid box' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dlyp' &
                , 'units', c='m')
    CALL channel_halt(substr, status)

    ! op_mk_20180108+
!--------------------   grid distance x (vector u)

    CALL new_channel_object(status, modstr, 'dlxu' &
         , p2 = out_dlxu, reprid = GP_2D_MPIOM, lstatic=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'dlxu' &
                , 'long_name', c='x length grid box vector u' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dlxu', 'units', c='m')
    CALL channel_halt(substr, status)

!--------------------   grid distance y (vector u)

    CALL new_channel_object(status, modstr, 'dlyu' &
         , p2 = out_dlyu, reprid = GP_2D_MPIOM, lstatic=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'dlyu' &
                , 'long_name', c='y length grid box vector u' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dlyu', 'units', c='m')
    CALL channel_halt(substr, status)

!--------------------   grid distance x (vector v)

    CALL new_channel_object(status, modstr, 'dlxv' &
         , p2 = out_dlxv, reprid = GP_2D_MPIOM, lstatic=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'dlxv' &
                , 'long_name', c='x length grid box vector v' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dlxv', 'units', c='m')
    CALL channel_halt(substr, status)

!--------------------   grid distance y (vector v)

    CALL new_channel_object(status, modstr, 'dlyv' &
         , p2 = out_dlyv, reprid = GP_2D_MPIOM, lstatic=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'dlyv' &
                , 'long_name', c='y length grid box vector v' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dlyv', 'units', c='m')
    CALL channel_halt(substr, status)
    ! op_mk_20180108-

!--------------------   land-sea mask

    CALL new_channel_object(status, modstr, 'lsm' &
         , p3 = mask, reprid = GP_3D_MPIOM, lstatic=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'lsm' &
                , 'long_name', c='land-sea mask' )
    CALL channel_halt(substr, status)

!--------------------   temperature

    CALL new_channel_object(status, modstr, 'tho' &
         , p3 = out_tho)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'tho' &
                , 'long_name', c='temperature' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tho' &
                , 'units', c='C')
    CALL channel_halt(substr, status)

!--------------------   temperature in K

    CALL new_channel_object(status, modstr, 'tho_K' &
         , p3 = out_tho_K)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'tho_K' &
                , 'long_name', c='temperature' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tho_K' &
                , 'units', c='K')
    CALL channel_halt(substr, status)

!--------------------   salinity

    CALL new_channel_object(status, modstr, 'sao' &
         , p3 = out_sao)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'sao' &
                , 'long_name', c='salinity' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'sao' &
                , 'units', c='psu')
    CALL channel_halt(substr, status)

! mz_bk_20110317+
!--------------------   ocean diagnostics
#if defined (MPIOM_2000)
    IF (ldiag) THEN
       CALL new_channel_object(status, modstr, 'global_sum_temperature' &
            , p0 = out_global_sum_temperature, reprid = SCALAR)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr , 'global_sum_temperature' &
            , 'long_name', c='global temperature sum' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'global_sum_temperature' &
            , 'units', c='degrees C')
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr, 'global_sum_salinity' &
            , p0 = out_global_sum_salinity, reprid = SCALAR)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr , 'global_sum_salinity' &
            , 'long_name', c='global salinity sum' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'global_sum_salinity' &
            , 'units', c='PSU')
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr, 'global_mass' &
            , p0 = out_global_mass, reprid = SCALAR)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr , 'global_mass' &
            , 'long_name', c='global mass' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'global_mass' &
            , 'units', c='kg')
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr, 'global_volume' &
            , p0 = out_global_volume, reprid = SCALAR)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr , 'global_volume' &
            , 'long_name', c='global volume' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'global_volume' &
            , 'units', c='m**3')
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr, 'global_salt_content' &
            , p0 = out_global_salt_content, reprid = SCALAR)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr , 'global_salt_content' &
            , 'long_name', c='global salt content' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'global_salt_content' &
            , 'units', c='PSU')
       CALL channel_halt(substr, status)
    ENDIF
#endif
! mz_bk_20110317-

!--------------------   z-velocity

    CALL new_channel_object(status, modstr, 'wo' &
         , p3 = out_wo, reprid=GP_3D_MPIOM_INT)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'wo' &
                , 'long_name', c='z-velocity' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wo' &
                , 'units', c='m/s')
    CALL channel_halt(substr, status)

!--------------------    total zonal velocity

    CALL new_channel_object(status, modstr, 'uoo' &
         , p3 = out_uoo)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'uoo' &
                , 'long_name', c='total zonal velocity' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'uoo' &
                , 'units', c='m/s')
    CALL channel_halt(substr, status)

!--------------------    total meridional velocity

    CALL new_channel_object(status, modstr, 'voe' &
         , p3 = out_voe)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'voe' &
                , 'long_name', c='total meridional velocity' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'voe' &
                , 'units', c='m/s')
    CALL channel_halt(substr, status)

    ! op_mk_20180108+
!--------------------    x velocity

    CALL new_channel_object(status, modstr, 'uko' &
         , p3 = out_uko)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'uko' &
                , 'long_name', c='x velocity' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'uko' &
                , 'units', c='m/s')
    CALL channel_halt(substr, status)

!--------------------    y velocity

    CALL new_channel_object(status, modstr, 'vke' &
         , p3 = out_vke)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'vke' &
                , 'long_name', c='y velocity' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vke' &
                , 'units', c='m/s')
    CALL channel_halt(substr, status)

!--------------------   negative of vertical density gradient (stability)

    CALL new_channel_object(status, modstr, 'stabio' &
         , p3 = out_stabio)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'stabio' &
                , 'long_name', c='negative of vertical density gradient' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'stabio' &
                , 'units', c='kg m-4')
    CALL channel_halt(substr, status)
    ! op_mk_20180108-

!!--------------------   x-velocity divergence free
!
!    CALL new_channel_object(status, modstr, 'ukomfl' &
!         , p3 = ukomfl_out)
!    CALL channel_halt(substr, status)
!    CALL new_attribute(status, modstr , 'ukomfl' &
!                , 'long_name', c='x-velocity divergence free' )
!    CALL channel_halt(substr, status)
!    CALL new_attribute(status, modstr, 'ukomfl' &
!                , 'units', c='m/s')
!    CALL channel_halt(substr, status)
!
!!--------------------   y-velocity divergence free
!
!    CALL new_channel_object(status, modstr, 'vkemfl' &
!         , p3 = vkemfl_out)
!    CALL channel_halt(substr, status)
!    CALL new_attribute(status, modstr , 'vkemfl' &
!                , 'long_name', c='x-velocity divergence free' )
!    CALL channel_halt(substr, status)
!    CALL new_attribute(status, modstr, 'vkemfl' &
!                , 'units', c='m/s')
!    CALL channel_halt(substr, status)
!
!--------------------   volume

    CALL new_channel_object(status, modstr, 'volume' &
         , p3 = omvol)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'volume' &
                , 'long_name', c='Box volume' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'volume' &
                , 'units', c='m^3')
    CALL channel_halt(substr, status)
!--------------------   mass

    CALL new_channel_object(status, modstr, 'mass' &
         , p3 = ommass)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'mass' &
                , 'long_name', c='water mass' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'mass' &
                , 'units', c='Kg')
    CALL channel_halt(substr, status)
!--------------------   insitu density

    CALL new_channel_object(status, modstr, 'rhoo' &
         , p3 = out_rhoo)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'rhoo' &
                , 'long_name', c='density' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rhoo' &
                , 'units', c='Kg/m^3')
    CALL channel_halt(substr, status)

!--------------------   pressure

    CALL new_channel_object(status, modstr, 'po' &
         , p3 = out_po)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'po' &
                , 'long_name', c='pressure' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'po' &
                , 'units', c='Pa')
    CALL channel_halt(substr, status)

!--------------------   vertical T,S diffusion

    CALL new_channel_object(status, modstr, 'dvo' &
         , p3 = out_dvo, reprid=GP_3D_MPIOM_INT)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'dvo' &
                , 'long_name', c='vertical diffusion (T,S)' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dvo' &
                , 'units', c='m^2 /s')
    CALL channel_halt(substr, status)

!--------------------   vertical momentum diffusion

    CALL new_channel_object(status, modstr, 'avo' &
         , p3 = out_avo, reprid=GP_3D_MPIOM_INT)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'avo' &
                , 'long_name', c='vertical momentum diffusion' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'avo' &
                , 'units', c='m^2 /s')
    CALL channel_halt(substr, status)

!--------------------   sea-level

    CALL new_channel_object(status, modstr, 'zo' &
         , p2 = out_zo, REPRID= GP_2D_MPIOM)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'zo' &
                , 'long_name', c='sea-level' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'zo' &
                , 'units', c='m')
    CALL channel_halt(substr, status)

    ! mz_bk_20111209+
!--------------------   ddpo

    CALL new_channel_object(status, modstr, 'ddpo' &
         , p3 = out_ddpo &
         , lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'ddpo' &
                , 'long_name', c='layer thickness at scalar points' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ddpo' &
                , 'units', c='m')
    CALL channel_halt(substr, status)
    ! mz_bk_20111209-
    ! op_mk_20180108+
!--------------------   dduo

    CALL new_channel_object(status, modstr, 'dduo' &
         , p3 = out_dduo &
         , lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'dduo' &
                , 'long_name', c='layer thickness at vector points' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dduo' &
                , 'units', c='m')
    CALL channel_halt(substr, status)

!--------------------   ddue

    CALL new_channel_object(status, modstr, 'ddue' &
         , p3 = out_ddue &
         , lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'ddue' &
                , 'long_name', c='layer thickness at vector points' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ddue' &
                , 'units', c='m')
    CALL channel_halt(substr, status)
    ! op_mk_20180108-

!--------------------   sea-level-change

    CALL new_channel_object(status, modstr, 'z1o' &
         , p2 = out_z1o, REPRID= GP_2D_MPIOM)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'z1o' &
                , 'long_name', c='sea-level change' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'z1o' &
                , 'units', c='m')
    CALL channel_halt(substr, status)

!--------------------   ice thickness

    CALL new_channel_object(status, modstr, 'sictho' &
         , p2 = out_sictho, REPRID= GP_2D_MPIOM)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'sictho' &
                , 'long_name', c='ice thickness' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'sictho' &
                , 'units', c='m')
    CALL channel_halt(substr, status)

    IF (.not.L_COUPLING) THEN !ice temperature calculated in E5
!--------------------   ice temperature

    CALL new_channel_object(status, modstr, 'tice' &
         , p2 = out_tice, REPRID= GP_2D_MPIOM)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'tice' &
                , 'long_name', c='ice temperature' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tice' &
                , 'units', c='m')
    CALL channel_halt(substr, status)

     ENDIF

!--------------------   ice compactness

    CALL new_channel_object(status, modstr, 'sicomo' &
         , p2 = out_sicomo, REPRID= GP_2D_MPIOM)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'sicomo' &
                , 'long_name', c='ice compactness' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'sicomo' &
                , 'units', c=' frac.')
    CALL channel_halt(substr, status)

!--------------------  x ice velocity

    CALL new_channel_object(status, modstr, 'sicuo' &
         , p2 = out_sicuo, REPRID= GP_2D_MPIOM)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'sicuo' &
                , 'long_name', c='x ice velocity' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'sicuo' &
                , 'units', c=' m/s')
    CALL channel_halt(substr, status)

!--------------------  y ice velocity

    CALL new_channel_object(status, modstr, 'sicve' &
         , p2 = out_sicve, REPRID= GP_2D_MPIOM)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'sicve' &
                , 'long_name', c='y ice velocity' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'sicve' &
                , 'units', c=' m/s')
    CALL channel_halt(substr, status)

!--------------------  non linear shear viscosity of ice

    CALL new_channel_object(status, modstr, 'hibete' &
         , p2 = out_hibete, REPRID= GP_2D_MPIOM)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'hibete' &
                , 'long_name', c='non linear shear viscosity of ice' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'hibete' &
                , 'units', c=' ...')
    CALL channel_halt(substr, status)

!--------------------  non linear shear viscosity of ice

    CALL new_channel_object(status, modstr, 'hibeto' &
         , p2 = out_hibeto, REPRID= GP_2D_MPIOM)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'hibeto' &
                , 'long_name', c='non linear shear viscosity of ice' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'hibeto' &
                , 'units', c=' ...')
    CALL channel_halt(substr, status)

!--------------------  non linear bulk viscosity of ice

    CALL new_channel_object(status, modstr, 'hibzeto' &
         , p2 = out_hibzeto, REPRID= GP_2D_MPIOM)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'hibzeto' &
                , 'long_name', c='non linear bulk viscosity of ice' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'hibzeto' &
                , 'units', c=' ...')
    CALL channel_halt(substr, status)


!--------------------  non linear bul viscosity of ice

    CALL new_channel_object(status, modstr, 'hibzete' &
         , p2 = out_hibzete, REPRID= GP_2D_MPIOM)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'hibzete' &
                , 'long_name', c='non linear bulk viscosity of ice' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'hibzete' &
                , 'units', c=' ...')
    CALL channel_halt(substr, status)

!--------------------   ocean westward velocity

    CALL new_channel_object(status, modstr, 'socu' &
         , p2 = out_socu, REPRID= GP_2D_MPIOM)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'socu' &
                , 'long_name', c='ocean westward velocity' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'socu' &
                , 'units', c='m s-1')
    CALL channel_halt(substr, status)

!--------------------   ocean northward velocity

    CALL new_channel_object(status, modstr, 'socv' &
         , p2 = out_socv, REPRID= GP_2D_MPIOM)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'socv' &
                , 'long_name', c='ocean northward velocity' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'socv' &
                , 'units', c='m s-1')
    CALL channel_halt(substr, status)



!--------------------   snow thickness

    CALL new_channel_object(status, modstr, 'sicsno' &
         , p2 = out_sicsno, REPRID= GP_2D_MPIOM)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'sicsno' &
                , 'long_name', c='snow thickness' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'sicsno' &
                , 'units', c='m')
    CALL channel_halt(substr, status)

    IF (.not.L_COUPLING) THEN

!--------------------  freshwater flux by restoring

    CALL new_channel_object(status, modstr, 'eminpo' &
         , p2 = out_eminpo, REPRID= GP_2D_MPIOM)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'eminpo' &
                , 'long_name', c='freshwater flux by restoring' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'eminpo' &
                , 'units', c='m/s')
    CALL channel_halt(substr, status)

!--------------------   river runoff
    CALL new_channel_object(status, modstr, 'riv_run' &
         , p2 = out_rivrun, REPRID= GP_2D_MPIOM)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr , 'riv_run' &
                , 'long_name', c='river runoff' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'riv_run' &
                , 'units', c='m/s')
    CALL channel_halt(substr, status)
    ENDIF

#if defined (MPIOM_13B)
    IF (L_TIDES) THEN
#elif defined (MPIOM_2000)
    IF (ltidal) THEN
#endif
    !--------------------   tidal phase, needed for rerun

        CALL new_channel_object(status, modstr, 'mmccdt' &
             , p0 = OUT_mmccdt, REPRID= SCALAR)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr , 'mmccdt' &
                    , 'long_name', c='daily tidal phase' )
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr, 'mmccdt' &
                    , 'units', c='')
        CALL channel_halt(substr, status)

    ENDIF

!-----------------------------------------------------------------------
!                        NULLIFY POINTERS
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!                          FIELD ALLOCATION
!-----------------------------------------------------------------------

    ! allocate decomposed field

    CALL set_param1
    CALL alloc_mem_commo1
    CALL alloc_mem_octdiff
    CALL alloc_mem_commoau2
    CALL alloc_mem_commoau3
    !CALL alloc_mem_diag
#if defined (MPIOM_2000)
    CALL alloc_mem_forcing
#endif


    CALL alloc_mem_dilcor

    CALL alloc_mem_elicom
    CALL alloc_mem_para2

    ! define as zero all the variables and parameters not allocatable
#ifdef MPIOM_13B
    ZERO=0.0
#endif
    CALL BELEG_ZERO

    ! Before namelist
    !As MPIOM... but obsolete. To be cleaned

#ifdef MPIOM_13B
     NYEARS=0
     NMONTS=1


     NANF=0
     NNNDT = 0
     NDAYS = 0
     DH = 0.
     AH = 0.
     FOUR=4.0
     TWO=2.
     ONE=1.0
     HALF=0.5
     FOURTH=0.25
     EIGHTH=0.125
     TENM4=1.E-4
     ALMZER=1.E-19

     monlen(:) = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

#ifdef YEAR360
         monlen(:) = 30
#endif /*YEAR360*/

#ifdef DEBUG_ONEDAY
         monlen(:) = 1
#endif /*DEBUG_ONEDAY*/


     WINTUR(1)=0.
     WINTUR(2)=1.E-3
     DO K=3,KEP
        WINTUR(K)=0.3*WINTUR(K-1)
     ENDDO
#endif

   ! DZwW(K) DECOFTES THE THICKNESS OF LAYER K IN METERS
   ! SUM OF DZW(1-KE) IS THE BOTTOM OF THE MODEL, I.E. SOLID GROUND !!!
   ! DZ(K) DECOFTES THE DISTANCE BETWEEN VECTOR POINTS BUT WILL BE
   ! COMPUTED IN SBR BODEN FOR CONSISTENCY WITH W-POINTS
     DO K=1,KE
        !mz_ap_20070731+
        ! change from the original value... these are more realistic values
        ! constant for all vertical layers
        SAF(K)=34.8
        TAF(K)=17.0
        !SAF(K)=34.8
        !TAF(K)=1.
        !mz_ap_20070731-
        DZW(K)=0.
     ENDDO

     DO K=1,KE
             DZW(k)=CDZW(k)
     ENDDO

#ifdef MPIOM_2000
     spongezone=(/1,1,4,ie_g,je_g,ke/)
      IF (cwa .LT. 0._wp) cwa = cwt

!done after
!      !HH   CHECK LAYER THICKNESS
!      IF (p_parallel_io) THEN
!      DO k=1, ke
!        ! FIXME: shouldn't this also check less than zero?
!        IF (dzw(k) .EQ. 0._wp) THEN
!          WRITE(*, *) ' LAYER: ', k, ' THICKNESS IS ZERO !!!'
!          CALL stop_all('check of layer thickness failed => run aborted')
!        ELSE
!          WRITE(*, *) k, ' LAYERTHICKNESS    (DZW): ', dzw(k)
!        END IF
!      END DO
!      END IF
#endif
     DTI=1./DT
     NDTDAY=NINT(86400./DT)

     IF (ISTART .LT. 0 .OR. ISTART .GT. 4)  THEN
! op_pj_20161104+
!!$      IF (p_pe==p_io) WRITE(*,*)'ISTART NOT SUPPORTED!!!!!'
!!$      CALL FINISH(substr)
         CALL error_bi('ISTART NOT SUPPORTED!',substr)
! op_pj_20161104-
     ENDIF
     IF (ISTART .LT. 1)  THEN
! op_pj_20161104+
!!$         IF (p_pe==p_io) THEN
!!$             WRITE(*,*)'ISTART=0 not desired in MESSy'
!!$             WRITE(*,*)'Please use MPIOM alone for this option'
!!$             WRITE(*,*)'(rewrite topography)'
!!$         ENDIF
!!$      CALL FINISH(substr)
         CALL error_bi('ISTART=0 not desired in MESSy; Please use MPIOM alone for this option (rewrite topography)',substr)
! op_pj_20161104-
     ENDIF
     IF (I3DREST .LT. 0 .OR. I3DREST .GT. 1) THEN
! op_pj_20161104+
!!$      IF (p_pe==p_io) WRITE(*,*)'I3DREST NOT SUPPORTED!!!!!'
!!$      CALL FINISH(substr)
         CALL error_bi('I3DREST NOT SUPPORTED!',substr)
! op_pj_20161104-
     ENDIF

#ifdef MPIOM_13B
     IF (I3DREST .GT. 0 .OR. ISTART .LT. 3) THEN
      CALL init_levitus(ie,je,ke)
     ENDIF
#elif defined (MPIOM_2000)
     write (*,*) "init levitus 2d"
     CALL init_levitus_2d
     IF (I3DREST .GT. 0 .OR. ISTART .LT. 3) THEN
     write (*,*) "init levitus 3d"
       CALL init_levitus_3d
     !TODO: why finish?????
     !ELSE
     !  CALL FINISH(substr)
     ENDIF
#endif

     AULAPTS=CAULAPTS*DT/3600.
     AULAPUV=CAULAPUV*DT/3600.
     AH00=CAH00/4.E5
#if defined (MPIOM_13B)
     WT=CWT/(6.**3)
#endif

   !HH   CHECK LAYER THICKNESS
     DO K=1,KE
       IF(DZW(K).EQ.0.) THEN
         IF (p_pe==p_io) THEN
           WRITE(*,*)' LAYER: ',K,' THICKNESS IS ZERO !!!'
         ENDIF
! op_pj_20161104+
!!$      CALL FINISH(substr)
         CALL error_bi('LAYER THICKNESS IS ZERO!',substr)
! op_pj_20161104-
       ELSE
         IF (p_pe==p_io) THEN
           WRITE(*,*) K, ' LAYERTHICKNESS    (DZW): ',DZW(K)
         ENDIF
       ENDIF
     ENDDO

#ifdef MPIOM_13B
   !HH   BACKGROUNDDIFFUSION
     DO K=1,KEP
       ABACKV(K) = ABACK
       DBACKV(K) = DBACK
     ENDDO
#endif

     TIESTW(1) = 0.
     DO K=1,KE
     TIESTW(K+1)  = TIESTW(K) + DZW(K)
     ENDDO

#if defined (MPIOM_13B)
     PI=4.*ATAN(1.)
     PIBOG=180./PI

#ifdef DBACKPROFIL
         IDBACK=0
         DO K=1,KEP
            ABACKV(K) = ABACK
            DBACKV(K) = DBACK + (1-IDBACK) * 1.E-4 *                       &
        &               (0.8 + 1.05/PI*ATAN(4.5*1.E-3*(TIESTW(K)-2500.)))* &
        &               (0.5+SIGN(0.5,TIESTW(K)-500.)) *                   &
        &               SQRT(ABS((TIESTW(K)-500.)/(3500.-500.)))
         ENDDO
#endif

         DO K=1,KEP
            GFDL_DIFF = 1.E-4 *                                            &
        &            (0.8 + 1.05/PI*ATAN(4.5*1.E-3*(TIESTW(K)-2500.)))
#ifdef DBACKGFDL
            DBACKV(K) = GFDL_DIFF
#endif
          ENDDO
#ifdef DBACKGFDL2
            DBACKV(K) = 0.5*(DBACK+GFDL_DIFF)
#endif
          DO K=1,KEP
            IF (p_pe==p_io) THEN
                WRITE(*,6002)'BACKGROUND DIFFUSIVITY AT '              &
            &                        ,INT(TIESTW(K))                           &
            &         ,'M : HOPE : ',DBACKV(K),' GFDL : ',GFDL_DIFF
                WRITE(*,6002)'BACKGROUND VISCOSITY AT ',INT(TIESTW(K)) &
            &         ,'M : HOPE : ',ABACKV(K)
            ENDIF
          ENDDO

    6002  FORMAT(1X,A27,I5,A10,E12.3,A8,E12.3)
#elif defined (MPIOM_2000)

    CALL setup_ocean_vertical_mixing(cwt, cwa, aback, dback)

    !  CHECK RIVER RUNOFF SCHEME
    IF ( numriv > 0 ) THEN
      luse_river_runoff_stations=.TRUE.
      IF (p_pe==p_io) WRITE(*,*)'river_runoff_stations enabled -> omip rivers disabled'
    ENDIF

    !  CHECK GLACIER CALVING SCHEME
    IF ( numglac > 0 ) THEN
      luse_glac_calv=.TRUE.
      IF (p_pe==p_io) WRITE(*,*)'glacier calving enabled'
    ENDIF

#endif


#if defined (MPIOM_13B)
   !HH    RELAXATION TIME SALINITY
   !      RELSAL = CRELSAL*20./DZW(1)
         RELSAL = CRELSAL
         IF (RELSAL.GT. ALMZER) THEN
   !         WRITE(IO_STDOUT,26668)1./(RELSAL*24.*3600.)
            if (p_pe==p_io) then
             WRITE(*,*) ' SURFACE SALINITY : '
             WRITE(*,*) ' RELAXATION DZW(1) [m]  =',DZW(1)
             WRITE(*,*) ' RELAXATION TIME [DAYS] =',1./(RELSAL*24.*3600.)
             WRITE(*,*) ' PISTON VELOCITY [m/s]  =',DZW(1)*(RELSAL)
             WRITE(*,*) ' RELAXATION TIME (relative to 20m) [DAYS] =', &
                          1./(RELSAL*24.*3600.*DZW(1)/20.)
          endif
         ELSE
            if (p_pe==p_io) then
               WRITE(*,*) 'SSS relaxation switched off !!'
              !mz_ap_20070801+
              RELSAL=0.0
              !mz_ap_20070801-
            endif
         ENDIF

   !HH    RELAXATION TIME TEMPERATURE
   !      RELTEM = CRELTEM*20./DZW(1)
         RELTEM = CRELTEM
         IF (RELTEM.GT. ALMZER) THEN
   !         WRITE(IO_STDOUT,26669)1./(RELTEM*24.*3600.)
          if (p_pe==p_io) then
             WRITE(*,*) 'SURFACE TEMPERATURE : '
             WRITE(*,*) 'RELAXATION DZW(1) [m]  =',DZW(1)
             WRITE(*,*)' RELAXATION TIME [DAYS] =',1./(RELTEM*24.*3600.)
             WRITE(*,*)' PISTON VELOCITY [m/s]  =',DZW(1)*(RELTEM)
             !WRITE(*,*)'  RELAXATION TIME TEMPERATURE COUPLING : '&
             !    ,DZW(1),' m /',(RELTEM*24.*3600.),' DAYS'
             WRITE(*,*)' RELAXATION TIME (relative to 20m) [DAYS] =', &
                         1./(RELTEM*24.*3600.*DZW(1)/20.)
          endif
        ELSE
           if (p_pe==p_io) then
              WRITE(*,*) 'SST relaxation switched off !!'
              !mz_ap_20070801+
              RELTEM=0.0
              !mz_ap_20070801-
           endif
         ENDIF


#elif defined (MPIOM_2000)
!RELSAL AND RELTEM OBSOLETE
   !HH    RELAXATION TIME SALINITY
         IF (CRELSAL.GT. ALMZER) THEN
   !         WRITE(IO_STDOUT,26668)1./(RELSAL*24.*3600.)
            if (p_pe==p_io) then
             WRITE(*,*) ' SURFACE SALINITY : '
             WRITE(*,*) ' RELAXATION DZW(1) [m]  =',DZW(1)
             WRITE(*,*) ' RELAXATION TIME [DAYS] =',1./(CRELSAL*24.*3600.)
             WRITE(*,*) ' PISTON VELOCITY [m/s]  =',DZW(1)*(CRELSAL)
             WRITE(*,*) ' RELAXATION TIME (relative to 20m) [DAYS] =', &
                          1./(CRELSAL*24.*3600.*DZW(1)/20.)
          endif
         ELSE
            if (p_pe==p_io) then
               WRITE(*,*) 'SSS relaxation switched off !!'
              !mz_ap_20070801+
              CRELSAL=0.0
              !mz_ap_20070801-
            endif
         ENDIF

   !HH    RELAXATION TIME TEMPERATURE
   !      RELTEM = CRELTEM*20./DZW(1)
         IF (CRELTEM.GT. ALMZER) THEN
   !         WRITE(IO_STDOUT,26669)1./(RELTEM*24.*3600.)
          if (p_pe==p_io) then
             WRITE(*,*) 'SURFACE TEMPERATURE : '
             WRITE(*,*) 'RELAXATION DZW(1) [m]  =',DZW(1)
             WRITE(*,*)' RELAXATION TIME [DAYS] =',1./(CRELTEM*24.*3600.)
             WRITE(*,*)' PISTON VELOCITY [m/s]  =',DZW(1)*(CRELTEM)
             !WRITE(*,*)'  RELAXATION TIME TEMPERATURE COUPLING : '&
             !    ,DZW(1),' m /',(CRELTEM*24.*3600.),' DAYS'
             WRITE(*,*)' RELAXATION TIME (relative to 20m) [DAYS] =', &
                         1./(CRELTEM*24.*3600.*DZW(1)/20.)
          endif
        ELSE
           if (p_pe==p_io) then
              WRITE(*,*) 'SST relaxation switched off !!'
              !mz_ap_20070801+
              CRELTEM=0.0
              !mz_ap_20070801-
           endif
         ENDIF

    IF ( lwith_barotropic_stokes_drift ) THEN
      CALL alloc_mem_stokes_drift
    ENDIF

#endif


#if defined (MPIOM_13B)
   !-----------------------------------------------------------------------

         IF(IOCAD.EQ.3.OR.IOCAD.EQ.4)THEN
            CALL alloc_mem_adpo
         ENDIF

         IF(IOCAD.EQ.4)THEN
            CALL alloc_mem_commobbl
         ENDIF

   !----------------------------------------------------------------------
   ! OPEN FILES
#ifndef LITTLE_ENDIAN
         IF(p_pe==p_io) THEN
           OPEN(IO_IN_ARCG,FILE='arcgri',                                  &
        &          ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
           OPEN(IO_IN_INIT,FILE='INITEM',STATUS='UNKNOWN',&
        &          ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
           OPEN(IO_IN_INIS,FILE='INISAL',STATUS='UNKNOWN',&
        &          ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
           OPEN(IO_IN_SURS,FILE='SURSAL',STATUS='UNKNOWN',&
        &          ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
!           IF (RELTEM .GT. ALMZER) THEN
             OPEN(IO_IN_SURT,FILE='SURTEM',STATUS='UNKNOWN',               &
        &            ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
!           ENDIF

#else
         IF(p_pe==p_io) THEN
#ifndef NOENDIANCONVERT
           OPEN(IO_IN_ARCG,FILE='arcgri',                                  &
        &          ACCESS='SEQUENTIAL',FORM='UNFORMATTED', CONVERT='BIG_ENDIAN')
           OPEN(IO_IN_INIT,FILE='INITEM',STATUS='UNKNOWN',                 &
        &          ACCESS='SEQUENTIAL',FORM='UNFORMATTED', CONVERT='BIG_ENDIAN')
           OPEN(IO_IN_INIS,FILE='INISAL',STATUS='UNKNOWN',                 &
        &          ACCESS='SEQUENTIAL',FORM='UNFORMATTED', CONVERT='BIG_ENDIAN')
           OPEN(IO_IN_SURS,FILE='SURSAL',STATUS='UNKNOWN',                 &
        &          ACCESS='SEQUENTIAL',FORM='UNFORMATTED', CONVERT='BIG_ENDIAN')
!           IF (RELTEM .GT. ALMZER) THEN
             OPEN(IO_IN_SURT,FILE='SURTEM',STATUS='UNKNOWN',               &
        &            ACCESS='SEQUENTIAL',FORM='UNFORMATTED', CONVERT='BIG_ENDIAN')
!           ENDIF

#else
    ! ERROR: compiler does not support convert='big_endian'
#endif
#endif

        IF (.not.L_COUPLING) THEN
            if (p_pe==p_io) WRITE(*,*) 'open omip'
           CALL open_omip
        ENDIF


        ENDIF ! p_pe == p_io

   !-------------------------------------------------------------------

   !     DU/DT - F*V = -G*( STABN* (DZ/DX)(NEW) + STABO* (DZ/DX)(OLD) )

   !     DZ/DT = H*( CONN* (DU/DX)(NEW) + CONO* (DU/DX)(OLD) + ... )

    STABN=0.6
    STABO=1.-STABN

    CONN=0.50
    CONO=1.-CONN

    IF (p_pe==p_io) WRITE(*,*)' STABN = ',STABN,' CONN = ',CONN

   !
   !-------------------------------------------------------------------
   !  OMEGA  :  ANGULAR VELOCITY OF OUR NICE BLUE PLANET
   !  RADIUS :  RADIUS OF THE ABOVE MENTIONED PLANET
   !  G      :  GRAVITATIONAL ACCELERATION ( CONSTANT 9.81 M/S**2
   !                TO ACCOUNT FOR THE EXISTENCE OF MINI BLACK HOLES )
   !  ROCP   :  HEAT CAPACITY PER CUBICMETER
   !  ROCD   :  WIND STRESS DRAG COEFFICIENT
   !------------------------------------------------------------------
         PI=4.*ATAN(1.)
         OMEGA=7.292E-5
         RADIUS=6371.E+3
         G=9.81
         GHN=G*STABN
         GHO=G*STABO
         ROCP=4.E06
         ROCD=1.2*1.5E-3
   !-------------------------------------------------------------------
   ! PARAMETER ICE MODEL

         ACLO(:,:)=0.7
         PAO(:,:)=101300.
         RPRECO(:,:)=0.0
         FRSE(:,:)=0.0

         ISNFLG=1

         ALBI=0.75
         ALBM=0.66
         ALBW=0.10
         ALBSN=0.85
         ALBSNM=0.75
#ifdef ALBMSN07
         ALBSN=0.82
         ALBSNM=0.7
         ALBM=0.63
#endif /*ALBMSN07*/

#ifdef ALBOMIP
   !SJM TUNE ALBEDOS FOR OMIP WITH BERYLIAND
         ALBI=0.75    !ICE < 0
         ALBM=0.70    !ICE > 0
         ALBW=0.10    !WATER
         ALBSN=0.85   !SNOW < 0
         ALBSNM=0.70  !SNOW > 0
#endif /*ALBOMIP*/
#ifdef ALBNCEP
   !HH   TUNE ALBEDOS FOR NCEP WITH KOCH
         ALBI=0.76
         ALBM=0.71
         ALBW=0.1
         ALBSN=0.86
         ALBSNM=0.72
#endif
#ifdef ALBMELTHI
         ALBSNM=0.8
         ALBSN=0.87
         ALBI=0.78
         ALBM=0.7
#endif /*ALBMELTHI*/
         TMELT=273.16
         TFREZ=-1.9
         CC=4.2E6
         CW=0.0045
         CLO=3.02E8
         CLB=2.70E8
         RHOAIR=1.3E+00
         RHOWAT=1.025E+03
         RHOICE=0.91E+03
         RHOSNO=0.33E+03
         RHOICWA=RHOICE/RHOWAT
         RHOSNWA=RHOSNO/RHOWAT
         RHOSNIC=RHOSNO/RHOICE
         CON=2.1656
         CONSN=0.31
         H0=0.5
         ARMIN=0.15
         ARMAX=1.0
         HMIN=0.05
   !UWE  MAXIMUM ICE THICKNESS FOR SNOW-->ICE CONVERSION
   !      HSNTOICE = 17.
          HSNTOICE = 0.45 * DZW(1)
   !UWE  MINIMUM ICE THICKNESS IN NEW ICE GROWTH
         SICTHMIN=0.5
   !
         VAPL=2.5E6
         SUBL=2.834E6
         SICE=5.0
         D3=5.5E-08
   !
   !JJ OLD VALUES UP TO HOPS 62!
         D1=RHOAIR*1004.*1.75E-3
         D2I=RHOAIR*SUBL*1.75E-3
         D2W=RHOAIR*VAPL*1.75E-3

#ifdef DRAGGILL
   !HH   VALUES FROM GILL ATMOSPHERE OCEAN DYNAMICS
         D1=RHOAIR*1004.*1.1E-3
         D2I=RHOAIR*SUBL*1.5E-3
         D2W=RHOAIR*VAPL*1.5E-3
#endif
#ifdef DASILVA
         D1=RHOAIR*1004.
         D2I=RHOAIR*SUBL
         D2W=RHOAIR*VAPL
#endif
#ifdef BULK_KARA
         D1=RHOAIR*1004.67
         D2I=RHOAIR*SUBL
         D2W=RHOAIR*VAPL
#endif

         TFREEZ=TFREZ
         ENTMEL=320.E6
         SICHEC=2.
         HICCE=2.
         HICCP=20.
         HTH00=0.5
         STEBOL=5.67E-8

         N=JTO
         N1=N-1
         N2=N-2
         N3=N-3
         N4=N-4

         M=IE
         M1=M-1
         M2=M-2
         M3=M-3
         M4=M-4

         KB=KBB
         KM=KB+1
         KBM=KB+KM

         DEUTO(:,:)=0.
         DEUTE(:,:)=0.
         DEPTO(:,:)=0.

         IF(p_pe==p_io) READ(IO_IN_ARCG) IBLA
         CALL read_slice(IO_IN_ARCG,DEUTO)

         IF (ISTART .EQ. 0) THEN

             CALL bounds_exch('p',DEUTO,'mpiom 10')

                DO J=1,JE
                   DO I=1,IE
                      DEPTO(I,J)=DEUTO(I,J)
                   ENDDO
                ENDDO


                DO I=1,IE
#ifndef bounds_exch_tp
                   IF(have_g_js) THEN
                      DEPTO(I,1)=0.
                      DEPTO(I,2)=0.
                   ENDIF
#endif
                   IF(have_g_je) THEN
                      DEPTO(I,JE)=0.
                      DEPTO(I,JE1)=0.
                   ENDIF
                ENDDO


              DO 8261 J=2,JE1
              DO 8261 I=2,IE1
                 DEUTO(I,J)=DEPTO(I,J)
                 IF(DEPTO(I,J).LT.1.) GO TO 8261
                 IF(DEPTO(I-1,J).LT.1..AND.DEPTO(I+1,J).LT.1.                  &
            &        .AND.DEPTO(I,J-1).LT.1..AND.DEPTO(I,J+1).LT.1.)           &
            &        DEUTO(I,J)=0.
       8261   CONTINUE
              CALL bounds_exch('u+',DEUTO,'mpiom 11')
              DO 8262 J=1,JE
              DO 8262 I=1,IE
                 DEPTO(I,J)=DEUTO(I,J)
       8262   CONTINUE

          ENDIF ! ISTART=0

          IF(p_pe==p_io) READ(IO_IN_ARCG) IBLA
          CALL read_slice(IO_IN_ARCG, DLXP)
          IF(p_pe==p_io) READ(IO_IN_ARCG) IBLA
          CALL read_slice(IO_IN_ARCG, DLXU)
          IF(p_pe==p_io) READ(IO_IN_ARCG) IBLA
          CALL read_slice(IO_IN_ARCG, DLXV)
          IF(p_pe==p_io) READ(IO_IN_ARCG) IBLA
          CALL read_slice(IO_IN_ARCG, DLYP)
          IF(p_pe==p_io) READ(IO_IN_ARCG) IBLA
          CALL read_slice(IO_IN_ARCG, DLYU)
          IF(p_pe==p_io) READ(IO_IN_ARCG) IBLA
          CALL read_slice(IO_IN_ARCG, DLYV)
          IF(p_pe==p_io) READ(IO_IN_ARCG) IBLA
          CALL read_slice(IO_IN_ARCG, FTWOU)
          IF(p_pe==p_io) READ(IO_IN_ARCG) IBLA
          CALL read_slice(IO_IN_ARCG, FTWOV)
          IF(p_pe==p_io) CLOSE(IO_IN_ARCG)


          DLXP(:,:)=MAX(1.,DLXP(:,:))
          DLXU(:,:)=MAX(1.,DLXU(:,:))
          DLXV(:,:)=MAX(1.,DLXV(:,:))
          DLYP(:,:)=MAX(1.,DLYP(:,:))
          DLYU(:,:)=MAX(1.,DLYU(:,:))
          DLYV(:,:)=MAX(1.,DLYV(:,:))

          CALL bounds_exch('v+',DLXV,'mpiom 12')
          CALL bounds_exch('p',DLXP,'mpiom 13')
          CALL bounds_exch('u+',DLXU,'mpiom 14')
          CALL bounds_exch('p',DLYP,'mpiom 15')
          CALL bounds_exch('u+',DLYU,'mpiom 16')
          CALL bounds_exch('v+',DLYV,'mpiom 17')

   !SL
   !SL GRID DEFORMATION
   !SL
          DO J=2,JE1
            DO I=2,IE1
                CURVAV(I,J)=(DLXP(I,J+1)-DLXP(I,J))/(DLYV(I,J)*DLXV(I,J))
            ENDDO
          ENDDO

          CALL bounds_exch('v+',CURVAV,'mpiom 18')


          DO 13788 I=1,IE
             IF(have_g_je) THEN
                CURVAV(I,JE)=CURVAV(I,JE-1)
                DLYP(I,JE)=DLYP(I,JE-1)
                DLYV(I,JE)=DLYV(I,JE-1)
                DLXV(I,JE)=DLXV(I,JE-1)
                DLXP(I,JE)=DLXP(I,JE-1)
                DLXU(I,JE)=DLXU(I,JE-1)
                DLYU(I,JE)=DLYU(I,JE-1)
             ENDIF
   13788  CONTINUE
   !

#elif defined (MPIOM_2000)

    IF (ladpo) THEN
      CALL alloc_mem_adpo
    END IF
    IF (ibbl_transport .EQ. 1) THEN
      CALL alloc_mem_commobbl
    END IF

    IF (IBOLK .NE. 0) THEN
      CALL alloc_mem_gmbolus
    END IF

    ! stabn defined in mo_param1.f90

    !-------------------------------------------------------------------

    !     DU/DT - F*V = -G*( STABN* (DZ/DX)(NEW) + STABO* (DZ/DX)(OLD) )

    !     DZ/DT = H*( CONN* (DU/DX)(NEW) + CONO* (DU/DX)(OLD) + ... )

    conn = 0.50_wp
    cono = 1._wp - conn

    IF (p_pe==p_io) WRITE(*,*)' STABN = ',STABN,' CONN = ',CONN
    !
    !-------------------------------------------------------------------
    ! PARAMETER ICE MODEL

    fclou(:,:) = 0.7_wp

    fprec(:,:) = 0.0_wp

    ! FIXME: isnflg is read in OCECTL but overwritten here
    ISNFLG = 1

    !
    HSNTOICE = HSNTOICE * DZW(1) ! fraction of the upper layer [m]

    N = JTO
    M=IE

    KM=KBB+1
    KBM=KBB+KM


    CALL setup_grid

#endif

         CALL BELEG

         CALL BODEN

   !
   !-------------------------------------------------------------------
   ! SPECIFY CORIOLIS PARAMETER ETC.
   !
         CALL CORIOL
   !
   !-------------------------------------------------------------------
#ifdef MPIOM_13B

         IF (p_pe==p_io) THEN
           WRITE(IO_STDOUT,*) 'DZ ', DZ,DI,DZW,TIESTU
           WRITE(IO_STDOUT,*)'WETO:'
         ENDIF
         DO JMM=0,JMMM
            JMANF=1+JMM*120
            JMEND=MIN((JMM+1)*120,JE_G)
            IF (p_pe==p_io) THEN
              WRITE(IO_STDOUT,*)'JMM ',JMM,JMANF,JMEND
              WRITE(IO_STDOUT,6061)0,(MOD(J,10),J=JMEND,JMANF,-1)
              WRITE(IO_STDOUT,*)'        J  <=== '
            ENDIF
   !         DO I=1,IE_G
   !            WRITE(IO_STDOUT,6061)I,(NINT(WETO_G(I,J,1))*MOD(J,10)       &
   !                 -10*NINT(WETO_G(I,J,1)-1.),J=JMEND,JMANF,-1)
   !         ENDDO
         ENDDO

         DO J=1,JE_G
            ICOU=0
            DO I=2,IE_G-1
               IF(WETO_G(I,J,1).GT.0.5)ICOU=ICOU+1
            ENDDO
            IF (p_pe==p_io) THEN
               WRITE(IO_STDOUT,*)'ZAEHL : ',J,ICOU
            ENDIF
            ENDDO
    6061       FORMAT(I4,1X,120I1)
   !-------------------------------------------------------------------
#endif

#ifdef MPIOM_13B
        IF(IOCAD.EQ.4)THEN
   ! PREPARE FOR BOTTOM BOUNDARY PARAMETRIZATIONS
           CALL FINDBOT(KBOT,WETO,IE,JE,KE)
           CALL FINDALFA
   !      PRINT*,'nach SLOPECON_ADPO'
        ENDIF
#elif defined  (MPIOM_2000)
    IF (ibbl_transport .EQ. 1) THEN
      ! PREPARE FOR BOTTOM BOUNDARY PARAMETRIZATIONS
      CALL FINDBOT(KBOT,WETO,IE,JE,KE)
      CALL FINDALFA
      !    PRINT*,'Preparation for bottom boundary parametrizations completed'
    END IF
#endif

   !-------------------------------------------------------------------
   !  provide matrix for barotropic mode
   !     uwe  add iteration of matrix
   !
         IF(ielimi.GE.1) THEN
   !!???? !mz_ap_20070604 always lower than 1 for us
   !! See mo_param1.f90
            CALL trian
         ELSE
            CALL itprep
            CALL trotest2
         ENDIF

#if defined (MPIOM_13B)
   !------------------------------------------------------------
   ! Tidal mode
      IF (L_TIDES) THEN
         CALL alloc_mem_tidal
      ENDIF


#elif defined (MPIOM_2000)

    CALL alloc_mem_swr_absorption

! TO DO FORCING (need CDI)
!    CALL initialize_surface_forcing

    IF ( luse_glac_calv ) THEN
      CALL glac_calv_ini
    ENDIF

    IF ( ltidal ) THEN
      CALL alloc_mem_tidal
    END IF

#endif

   !-------------------------------------------------------------
   !  GRID AND CPU INFORMATION (OUTPUTTED)
   !-------------------------------------------------------------
    DO I=1,IE
      DO J=1,JE
#ifdef MPIOM_13B
       lon_mpiom(I,J) = gila(2*I,2*J)*180._dp/pi
       lat_mpiom(I,J) = giph(2*I,2*J)*180._dp/pi
#elif defined (MPIOM_2000)
       lon_mpiom(I,J) = gila(2*I,2*J)*180._dp/api
       lat_mpiom(I,J) = giph(2*I,2*J)*180._dp/api
#endif
       out_dlxp(I,J)  = dlxp(I,J)
       out_dlyp(I,J)  = dlyp(I,J)
       ! op_mk_20180108+
       out_dlxu(I,J)  = dlxu(I,J)
       out_dlyu(I,J)  = dlyu(I,J)
       out_dlxv(I,J)  = dlxv(I,J)
       out_dlyv(I,J)  = dlyv(I,J)
       ! op_mk_20180108-
       DO K=1,KE
         decomp_mpiom(I,J,K)=p_pe
         MASK(I,J,K) = WETO(I,J,K)
         IF (MASK(I,J,K) == 0.0_dp) MASK(I,J,K) = (-1.E34_dp)
         IF (K==1) THEN
            ! total volume (without ice and snow)
            OMVOL(I,J,K) = dlyp(I,J)*dlxp(I,J)* &
                           (zo(I,J)+ddpo(I,J,1) &
                           ! mz_bk_20110311+
                           -sictho(i,j)*rhoicwa-sicsno(i,j)*rhosnwa)
                           ! mz_bk_20110311-
         ELSE
             OMVOL(I,J,K) = dlyp(I,J)*dlxp(I,J)* &
                            ddpo(I,J,K)
         ENDIF
       END DO
      END DO
    END DO

    CALL nullify_borders(status,OMVOL(:,:,:))


    CALL end_message_bi(modstr, 'MEMORY INITIALIZATION', substr)

! um_ak_20130705+ moved here from above
! um_ak_20130620+ ! mz_ap_20130620+
!-----------------------------------------------------------------------
!                      GRID DEFINTION
!-----------------------------------------------------------------------
    CALL start_message_bi(modstr, 'MPIOM GRID INITIALIZATION', substr)
    ! define ocean grid for use in other model part e.g. IMPORT-GRID
    CALL define_mpiom_grid
    CALL end_message_bi(modstr, 'MPIOM GRID INITIALIZATION', substr)
! um_ak_20130620- ! mz_ap_20130620-
! um_ak_20130705- moved here from above

  END SUBROUTINE mpiom_init_memory

  ! ---------------------------------------------------------------------------

  SUBROUTINE mpiom_init_coupling

    ! MESSy
    USE messy_main_data_bi,          ONLY: lcouple
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_2D_MPIOM
    USE messy_main_channel,          ONLY: get_channel_object, get_channel_info &
                                         , get_channel_object_info
#ifdef MPIOM_2000
    USE messy_main_mpi_bi,           ONLY: p_parallel_io, p_io, p_pe
#else
    USE messy_main_mpi_bi,           ONLY: p_parallel_io
#endif
    ! mz_bk_20110216-
    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'mpiom_init_coupling'
    !I/O
    INTEGER :: status

    CALL start_message_bi(modstr, 'COUPLING', substr)

   !-------------------------------------------------------------------
   !              FORCING ATMOSPHERE ?
   !-------------------------------------------------------------------

   ! COUPLING OCEAN --> ATMOSPHERE

      IF(lcouple) THEN
        CALL get_channel_info(status,cname='a2o')
        IF (status /= 0) THEN
! op_pj_20161104+
!!$       WRITE(*,*) "A2O submodel OFF, coupling OCEAN-->ATMOSPHERE not possible"
!!$       CALL finish(substr)
          CALL error_bi(&
               'A2O submodel OFF, coupling OCEAN-->ATMOSPHERE not possible"' &
               , substr)
! op_pj_20161104-
        ENDIF
      ENDIF

   ! COUPLING ATMOSPHERE --> OCEAN

      IF(L_COUPLING) THEN
        CALL get_channel_info(status,cname='a2o')
        IF (status /= 0) THEN
! op_pj_20161104+
!!$       WRITE(*,*) "A2O submodel OFF, coupling ATMOSPHERE-->OCEAN not possible"
!!$       CALL finish(substr)
          CALL error_bi( &
               'A2O submodel OFF, coupling ATMOSPHERE-->OCEAN not possible' &
               , substr)
! op_pj_20161104-
        ENDIF
      ENDIF

   !-------------------------------------------------------------------
   !   RUN OFF CALCULATION (FRESH WATER) <-- HD ?
   !-------------------------------------------------------------------

   ! HD model present??

      IF(L_COUPLING) THEN
        CALL get_channel_info(status,cname='hd')
        IF (status /= 0) THEN
          IF (p_pe==p_io) write(*,*) 'l_coupling active but no HD model active'
!         CALL FINISH(substr)
        ENDIF
      ENDIF

      IF(.not.L_COUPLING) THEN
#if defined (MPIOM_13B)
        call river_runoff_ini ! not needed <-- from read_omip
#elif defined (MPIOM_2000)
        IF ( luse_river_runoff_stations ) THEN
          CALL river_runoff_stations_ini
        ELSE
          CALL river_runoff_omip_ini
        ENDIF
#endif
!        MPIOM compiled with RIVER_GIRIV option
      ENDIF
   !-------------------------------------------------------------------

   !FIRST FOUR DIVIDED BY rhowat! Done later...

!      AOFLTXWO(:,:) !Zonal wind stress on water
      IF (L_COUPLING) THEN
         CALL get_channel_object(status          &
               , oname='AOFLTXWO', cname='a2o'  &
               , p2=AOFLTXWO )
         IF (status /= 0) THEN
! op_pj_20161104+
!!$         WRITE(*,*) "MISSING AOFLTXWO from a2o"
!!$      CALL finish(substr)
         CALL error_bi('MISSING AOFLTXWO from a2o',substr)
! op_pj_20161104-
         ENDIF
!  AOFLTYWO(:,:) !Meridional wind stress on water
          CALL get_channel_object(status          &
               , oname='AOFLTYWO', cname='a2o'  &
               , p2=AOFLTYWO )
         IF (status /= 0) THEN
! op_pj_20161104+
!!$         WRITE(*,*) "MISSING AOFLTYWO from a2o"
!!$      CALL finish(substr)
         CALL error_bi('MISSING AOFLTYWO from a2o',substr)
! op_pj_20161104+
         ENDIF
!  AOFLTXIO(:,:) !Zonal wind stress on snow/ice
          CALL get_channel_object(status          &
               , oname='AOFLTXIO', cname='a2o'  &
               , p2=AOFLTXIO )
         IF (status /= 0) THEN
! op_pj_20161104+
!!$      WRITE(*,*) "MISSING AOFLTXIO from a2o"
!!$      CALL finish(substr)
         CALL error_bi('MISSING AOFLTXIO from a2o',substr)
! op_pj_20161104-
         ENDIF
!  AOFLTYIO(:,:) !Meridional wind stress on snow/ice
          CALL get_channel_object(status          &
               , oname='AOFLTYIO', cname='a2o'  &
                , p2=AOFLTYIO)
         IF (status /= 0) THEN
! op_pj_20161104+
!!$      WRITE(*,*) "MISSING AOFLTYIO from a2o"
!!$      CALL finish(substr)
         CALL error_bi('MISSING AOFLTYIO from a2o',substr)
! op_pj_20161104-
         ENDIF
!  AOFLFRIO(:,:) !solid freshwater flux (over ice only)
          CALL get_channel_object(status          &
               , oname='AOFLFRIO', cname='a2o'  &
               , p2=AOFLFRIO)
         IF (status /= 0) THEN
! op_pj_20161104+
!!$      WRITE(*,*) "MISSING AOFLFRIO from a2o"
!!$      CALL finish(substr)
         CALL error_bi('MISSING AOFLFRIO from a2o',substr)
! op_pj_20161104-
         ENDIF
!  AOFLFRWO(:,:) !liquid freshwater flux (over water and ice)
          CALL get_channel_object(status          &
               , oname='AOFLFRWO', cname='a2o'  &
               , p2=AOFLFRWO)
         IF (status /= 0) THEN
! op_pj_20161104+
!!$      WRITE(*,*) "MISSING AOFLFRWO from a2o"
!!$      CALL finish(substr)
         CALL error_bi('MISSING AOFLFRWO from a2o',substr)
! op_pj_20161104-
         ENDIF
!  AOFLRHIO(:,:) !residual heat flux (sea-ice topmelt heat flux)
          CALL get_channel_object(status          &
               , oname='AOFLRHIO', cname='a2o'  &
               ,p2=AOFLRHIO)
         IF (status /= 0) THEN
! op_pj_20161104+
!!$      WRITE(*,*) "MISSING AOFLRHIO from a2o"
!!$      CALL finish(substr)
         CALL error_bi('MISSING AOFLRHIO from a2o',substr)
! op_pj_20161104-
         ENDIF
!  AOFLCHIO(:,:) !conductive heat flux through ice
          CALL get_channel_object(status          &
               , oname='AOFLCHIO', cname='a2o'  &
               ,p2=AOFLCHIO)
         IF (status /= 0) THEN
! op_pj_20161104+
!!$      WRITE(*,*) "MISSING AOFLCHIO from a2o"
!!$      CALL finish(substr)
         CALL error_bi('MISSING AOFLCHIO from a2o',substr)
! op_pj_20161104-
         ENDIF
!  AOFLNHWO(:,:) !net heat flux over water
          CALL get_channel_object(status          &
               , oname='AOFLNHWO', cname='a2o'  &
               ,p2=AOFLNHWO)
         IF (status /= 0) THEN
! op_pj_20161104+
!!$      WRITE(*,*) "MISSING AOFLNHWO from a2o"
!!$      CALL finish(substr)
         CALL error_bi('MISSING AOFLNHWO from a2o',substr)
! op_pj_20161104-
         ENDIF
!  AOFLSHWO(:,:) !downwelling solar radiation
          CALL get_channel_object(status          &
               , oname='AOFLSHWO', cname='a2o'  &
               ,p2=AOFLSHWO)
         IF (status /= 0) THEN
! op_pj_20161104+
!!$      WRITE(*,*) "MISSING AOFLSHWO from a2o"
!!$      CALL finish(substr)
         CALL error_bi('MISSING AOFLSHWO from a2o',substr)
! op_pj_20161104-
         ENDIF
!  AOFLWSVO(:,:) !wind stress velocity
          CALL get_channel_object(status          &
               , oname='AOFLWSVO', cname='a2o'  &
               ,p2=AOFLWSVO)
         IF (status /= 0) THEN
! op_pj_20161104+
!!$      WRITE(*,*) "MISSING AOFLWSVO from a2o"
!!$      CALL finish(substr)
         CALL error_bi('MISSING AOFLWSVO from a2o',substr)
! op_pj_20161104-
         ENDIF
      ENDIF
   !-------------------------------------------------------------------
   !              HAMOCC COUPLING
   !-------------------------------------------------------------------

    IF (L_HAMOCC_COUPLING) THEN
      CALL get_channel_object(status          &
            , oname='abs_oce', cname='hamocc'  &
            , p3=abs_oce )
      IF (status /= 0) THEN
         IF (p_pe==p_io) THEN
           WRITE(*,*) "NO abs_oce in HAMOCC"
           WRITE(*,*) "L_HAMOCC_COUPLING = FALSE"
         ENDIF
         L_HAMOCC_COUPLING = .FALSE.
         ALLOCATE(abs_oce(ie,je,ke))
         abs_oce(:,:,:) = 0.0_dp
      ELSE
         IF (p_pe==p_io) WRITE(*,*) "L_HAMOCC_COUPLING = TRUE"
      ENDIF
    ELSE
       ALLOCATE(abs_oce(ie,je,ke))
       abs_oce(:,:,:) = 0.0_dp
       L_HAMOCC_COUPLING = .FALSE.
       IF (p_pe==p_io) WRITE(*,*) "L_HAMOCC_COUPLING = FALSE"
    ENDIF

    CALL end_message_bi(modstr, 'COUPLING', substr)

  END SUBROUTINE mpiom_init_coupling

  ! ---------------------------------------------------------------------------

  SUBROUTINE mpiom_global_start


!!$  USE messy_main_mpi_bi,    ONLY: finish
#ifdef MPIOM_2000
  USE messy_main_mpi_bi,           ONLY: p_io, p_pe
#endif
  USE messy_main_grid_def_mem_bi,  ONLY: jrow, ngpblks, kproma
  USE messy_main_timer,     ONLY: time_step_len &
                                , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND &
                                , lstart, lresume
  ! um_ak_20080414+
!  USE messy_main_bmluse_bi, ONLY: get_date_components, current_date &
!                                , time_days, event_state
  USE messy_main_timer_bi,  ONLY: event_state
  USE messy_main_timer,     ONLY: current_date
  ! um_ak_20080414-
  USE messy_mpiom_tools_e5, ONLY: nullify_borders, read_file_3d, read_file_2d


  IMPLICIT NONE

  INTEGER :: IM1,JM1
  INTEGER :: status

  INTRINSIC INT

  CHARACTER(LEN=*), PARAMETER :: substr = 'mpiom_global_start'

!  ocean timestep= DT
!  messy/echam timestep = time_step_len
! l_trig_mpiom = TRUE  --> do the time loop
!                FALSE --> skip this time step of E5/M1

  l_trig_mpiom       = event_state(ev_trig_mpiom, current_date)
  l_trig_mpiom_day   = event_state(ev_trig_mpiom_day, current_date)
  l_trig_mpiom_month = event_state(ev_trig_mpiom_month, current_date)
  l_trig_mpiom_year  = event_state(ev_trig_mpiom_year, current_date)


! -------------------------------------------------------------
!  UPDATE THE TIME MANAGER (WE ARE IN E5/M1 LOOP!)
!------------------------------------------------------------------------

   lyears = YEAR
   lmonts = MONTH
   ldays  = DAY
   lyear1 = YEAR
   lmont1 = MONTH

!------------------------------------------------------------------------
! NEW START --> INITIAL CONDITIONS
!------------------------------------------------------------------------
!
!     START FROM STATUS LEVITUS OR HORIZONTAL STRATIFICATION
!     CALL LEVIRE WITH ARGUMENT =  0 : HORIZONTAL STRATIFICATION
!     CALL LEVIRE WITH ARGUMENT = 13 : 3D STRATIFICATION

      IF (lstart) THEN

#if defined (MPIOM_13B)
        ! rewind and read levitus climatology directly in sao and tho
        ! no corrections for ice
        IF (ISTART .EQ. 0) THEN
          ! loop to obtain the actual month!
          IF(p_pe==p_io) THEN
            REWIND(IO_IN_INIT)
            REWIND(IO_IN_INIS)
          ENDIF
          DO IREADC=1,LMONTS
            CALL LEVIRE(13)
          ENDDO
        ENDIF
        ! constant sao and tho (= saf, taf)
        IF (ISTART .EQ. 1) CALL LEVIRE(0)
        ! rewind and read levitus climatology directly in sao and tho
        ! corrections for ice
        IF (ISTART .EQ. 2) THEN
          ! loop to obtain the actual month
          IF(p_pe==p_io) THEN
            REWIND(IO_IN_INIT)
            REWIND(IO_IN_INIS)
          ENDIF
          DO IREADC=1,LMONTS
            CALL LEVIRE(13)
          ENDDO
          DO J=1,JE
           DO I=1,IE
            IF (THO(I,J,1)*WETO(I,J,1).LT.-0.5) THEN
              SICTHO(I,J)=MAX(0.,3.*(-0.5-THO(I,J,1)))
              SICOMO(I,J)=1.
            ENDIF
           ENDDO
          ENDDO
        ENDIF ! ISTART
#elif defined (MPIOM_2000)
        !istart_new_topo_update =0 => start from climatology
        !istart_new_run = 1        => start from unifirm ts-profile
        !istart_new_topo =2        => start from climatology
        !istart_new_continueo =3   => start from rerun
        IF (ISTART .EQ. istart_new_topo_update) THEN
          DO IREADC=1,LMONTS
            CALL levitus_read_3d_stratification
          ENDDO
        ENDIF
        IF (ISTART .EQ. istart_new_run) CALL levitus_horizontal_stratification
        IF (ISTART .EQ. istart_new_topo) THEN
          DO IREADC=1,LMONTS
            CALL levitus_read_3d_stratification
          ENDDO
          DO J=1,JE
           DO I=1,IE
            IF (THO(I,J,1)*WETO(I,J,1).LT.-0.5) THEN
              SICTHO(I,J)=MAX(0.,3.*(-0.5-THO(I,J,1)))
              SICOMO(I,J)=1.
            ENDIF
           ENDDO
          ENDDO
        ENDIF ! ISTART
#endif
!mz_ap_20160929+
        IF (ISTART .EQ. 4) THEN
           IF (mpiom_init_file_name=='') THEN
! op_pj_20161104+
!!$          write(*,*) 'missing netcdf file for initial conditions'
!!$          CALL finish(substr)
             CALL error_bi('missing netcdf file for initial conditions',substr)
! op_pj_20161104-
           ENDIF
           allocate(input_data3d(ie_g,je_g,ke))
           allocate(input_data2d(ie_g,je_g))
              IF (p_parallel_io) call read_file_3d(mpiom_init_file_name, 'tho', input_data3d)
              do K=1,ke
                call scatter_arr(input_data3d(:,:,k),tho(:,:,k), p_io)
              enddo
              IF (p_parallel_io) call read_file_3d(mpiom_init_file_name, 'sao', input_data3d)
              do K=1,ke
                call scatter_arr(input_data3d(:,:,k),sao(:,:,k), p_io)
              enddo
              IF (p_parallel_io) call read_file_2d(mpiom_init_file_name, 'z1o',    input_data2d)
              call scatter_arr(input_data2d(:,:),           z1o(:,:), p_io)
              IF (p_parallel_io) call read_file_2d(mpiom_init_file_name, 'sictho', input_data2d)
              call scatter_arr(input_data2d(:,:),           sictho(:,:), p_io)
              IF (p_parallel_io) call read_file_2d(mpiom_init_file_name, 'sicomo', input_data2d)
              call scatter_arr(input_data2d(:,:),           sicomo(:,:), p_io)
              IF (p_parallel_io) call read_file_2d(mpiom_init_file_name, 'sicuo',  input_data2d)
              call scatter_arr(input_data2d(:,:),           sicuo(:,:), p_io)
              IF (p_parallel_io) call read_file_2d(mpiom_init_file_name, 'sicve',  input_data2d)
              call scatter_arr(input_data2d(:,:),           sicve(:,:), p_io)
              IF (p_parallel_io) call read_file_2d(mpiom_init_file_name, 'sicsno', input_data2d)
              call scatter_arr(input_data2d(:,:),           sicsno(:,:), p_io)
              IF (p_parallel_io) call read_file_2d(mpiom_init_file_name, 'hibete', input_data2d)
              call scatter_arr(input_data2d(:,:),           hibete(:,:), p_io)
              IF (p_parallel_io) call read_file_2d(mpiom_init_file_name, 'hibeto', input_data2d)
              call scatter_arr(input_data2d(:,:),           hibeto(:,:), p_io)
              IF (p_parallel_io) call read_file_2d(mpiom_init_file_name, 'hibzete',input_data2d)
              call scatter_arr(input_data2d(:,:),           hibzete(:,:), p_io)
              IF (p_parallel_io) call read_file_2d(mpiom_init_file_name, 'hibzeto',input_data2d)
              call scatter_arr(input_data2d(:,:),           hibzeto(:,:), p_io)
              IF (p_parallel_io) call read_file_3d(mpiom_init_file_name, 'dvo',    input_data3d)
              do K=1,ke
                call scatter_arr(input_data3d(:,:,k),dvo(:,:,k), p_io)
              enddo
              IF (p_parallel_io) call read_file_3d(mpiom_init_file_name, 'avo', input_data3d)
              do K=1,ke
                call scatter_arr(input_data3d(:,:,k),avo(:,:,k), p_io)
              enddo
              call read_file_3d(mpiom_init_file_name, 'wo', input_data3d)
              do K=1,ke
                call scatter_arr(input_data3d(:,:,k),wo(:,:,k), p_io)
              enddo
           deallocate(input_data3d)
           deallocate(input_data2d)
        ENDIF
!mz_ap_20160929-

        !------------------------------------------------------------
        ! MASK VALUES OVER LAND!!!!!!!!!!
        !------------------------------------------------------------
        ! here...

        !------------------------------------------------------------
        !  CHECK GLOBAL SALT CONTENT
        IF ( ICONTRO .ne. 0 ) THEN
           CALL CONTRO(24)
        ENDIF
        !------------------------------------------------------------
        ! initial values:
        IF (L_CHECK) CALL mpiom_check("start of simulation")
        !------------------------------------------------------------
        ! Tidal mode
        ! see mo_tidal.f90:
        ! write (*,*) "YEAR,MONTH",YEAR,MONTH
#if defined (MPIOM_13B)
        IF (L_TIDES) THEN
#elif defined (MPIOM_2000)
        IF (ltidal) THEN
#endif
           CALL foreph_ini
        ENDIF
        !------------------------------------------------------------
        ! initial condition: needed at lstart to force the atmosphere
        ! mz_bk_20110315+
#if defined (MPIOM_2000)
        IF (ldiag) THEN
!          CALL calc_global_mean
           out_global_sum_temperature = global_sum_temperature
           out_global_sum_salinity = global_sum_salinity
           out_global_mass = global_mass
           out_global_volume = global_volume
           out_global_salt_content = global_salt_content
        ENDIF
#endif
        ! mz_bk_20110315-

        DO I=1,IE
          DO J=1,JE
           DO k=1,KE
              OUT_THO_K(I,J,K) = (THO(I,J,K)+tmelt)!
              OUT_THO(I,J,K)   = (THO(I,J,K))
              OUT_SAO(I,J,K)  =  SAO(I,J,K)
           END DO
           OUT_SICTHO(I,J)  = SICTHO(I,J)
           OUT_SICOMO(I,J)  = SICOMO(I,J)
           OUT_SICSNO(I,J)  = SICSNO(I,J)
          END DO
        END DO

      ENDIF !(lstart)


!------------------------------------------------------------------------
! INPUT FILES  FOR RESTORING BOUNDARY CONDITIONS (SURFACE AND 3D)
!------------------------------------------------------------------------

      IF (l_trig_mpiom_month.or.lstart.or.lresume) THEN
#if defined (MPIOM_13B)
         !------------------------------------------------------------------------
         ! 3D - RESTORING INPUT FILES (CONTROLLED BY I3DREST IN THE NAMELIST)
         !------------------------------------------------------------------------
         IF (I3DREST .EQ. 2) THEN
            !------------------ RESTORING ANNUAL CLIMATOLOGY -------------------
            LMONTS=0
            CALL LEVIRE(-1)
            LMONTS=MONTH
         ENDIF
         IF (I3DREST .EQ. 1) THEN
            !------------------ RESTORING MONTHLY CLIMATOLOGY -------------------
            !:: REWIND INISAL/INITEM TO BEGIN OF YEAR
            IF(p_pe==p_io) THEN
              REWIND(IO_IN_INIT)
              REWIND(IO_IN_INIS)
            ENDIF
            !:: READ APPROPRIATE INISAL/INITEM FOR MONTHLY RESTORING
            IF (p_pe==p_io) THEN
              WRITE(IO_STDOUT,*)                                               &
            'READING MONTHLY LEVITUS FIELDS IN MONTH ',LMONTS
            ENDIF
            DO IREADC=1,LMONTS
             CALL LEVIRE(-1)
            ENDDO
         ENDIF
         !------------------------------------------------------------------------
         ! SURFACE RESTORING INPUT FILES
         ! (CONTROLLED BY CRELTEM AND CRELSAL IN THE NAMELIST)
         !------------------------------------------------------------------------
          IF (RELSAL .GT. ALMZER) THEN
             IF (p_pe==p_io) THEN
               WRITE(IO_STDOUT,*) " READ SALINITY SURFACE RESTORING VALUE"
             ENDIF
             IF (p_pe==p_io) REWIND(IO_IN_SURS)   !mz_ap_2007062+
             DO IREADC=1,LMONTS
               CALL LEVIRE(-2) !mz_ap_20070629+
             ENDDO
          ENDIF
          IF (RELTEM .GT. ALMZER) THEN
             IF (p_pe==p_io) THEN
               WRITE(IO_STDOUT,*) " READ TEMPERATURE SURFACE RESTORING VALUE"
             ENDIF
             IF (p_pe==p_io) REWIND(IO_IN_SURT)
             DO IREADC=1,LMONTS
               CALL LEVIRE(-3)
             ENDDO
          ENDIF

#elif defined (MPIOM_2000)
         !------------------------------------------------------------------------
         ! 3D - RESTORING INPUT FILES (CONTROLLED BY I3DREST IN THE NAMELIST)
         ! I3DREST OPTIONS FOR 3-D RESTORING
         ! I3DREST == 0: NO RESTORINg (DEFAULT)
         !------------------------------------------------------------------------
         IF (I3DREST .GT. 0) THEN
           CALL levitus_read_3d_restore
         ENDIF
         CALL levitus_per_month_setup(lmont)
#endif

      ENDIF

!------------------------------------------------------------------------------
!                       FORCING FILES
!------------------------------------------------------------------------------

#if defined (MPIOM_13B)
      IF (.not.L_COUPLING) THEN
         IF (l_trig_mpiom_month.or.lstart.or.lresume) THEN
             ! we keep l_trig_mpiom_month in case of 29 days!
             !---------------------REWIND FORCING FILES ----------------------------
             CALL rewind_omip
             !---------------------SPOOL FORCING FILES ----------------------------
             IF (p_pe==p_io) THEN
               WRITE(IO_STDOUT,*) 'call spool omip'
             ENDIF
             CALL spool_omip     ! reach the correct month
             IF (lstart.or.lresume) THEN
                 CALL spool_omip_day ! reach the correct day
             ENDIF
         ENDIF
         IF (l_trig_mpiom_day.or.lstart.or.lresume) THEN
             IF (p_pe==p_io) THEN
               WRITE(IO_STDOUT,*)'READ SURFACE FORCING AT TIMESTEP '  &
                                  ,LDAYS,LMONTS,LYEARS
               WRITE(IO_STDOUT,*) 'call read omip'
             ENDIF
             CALL read_omip
         ENDIF
      ENDIF
#elif defined (MPIOM_2000)
  ! TODO OMIP forcing reading
  !see mo_forcing!
  ! need CDI
#endif

!------------------------------------------------------------------------
! START FROM RESTART
!------------------------------------------------------------------------

      IF (lresume) THEN
        IF (p_pe==p_io) THEN
           WRITE(IO_STDOUT,*)'FIRST MPIOM TIME STEP AFTER RERUN (FROM PREVIOUS CALCULATION)'
        ENDIF
        ! mz_bk_20110315+
#if defined (MPIOM_2000)
        IF (ldiag) THEN
           global_sum_temperature = out_global_sum_temperature
           global_sum_salinity = out_global_sum_salinity
           global_mass = out_global_mass
           global_volume = out_global_volume
           global_salt_content = out_global_salt_content
        ENDIF
#endif
        ! mz_bk_20110315-
        DO I=1,IE
          DO J=1,JE
           DO k=1,KE
             ! see later
             THO(I,J,K)  =  OUT_THO(I,J,K)
             SAO(I,J,K)  =  OUT_SAO(I,J,K)
! cy_ap_20091022+
             UOO(I,J,K)  =  OUT_UOO(I,J,K)
             VOE(I,J,K)  =  OUT_VOE(I,J,K)
! cy_ap_20091022-
             ! op_mk_20180108+
             UKO(I,J,K)  =  OUT_UKO(I,J,K)
             VKE(I,J,K)  =  OUT_VKE(I,J,K)
             DDUO(I,J,K) =  OUT_DDUO(I,J,K)
             DDUE(I,J,K) =  OUT_DDUE(I,J,K)
             STABIO(I,J,K) = OUT_STABIO(I,J,K)
             ! op_mk_20180108-
             ! mz_bk_20111209+
             DDPO(I,J,K) = OUT_DDPO(I,J,K)
             ! mz_bk_20111209-
           END DO
           DO k=1,KE+1
             DVO(I,J,K)  =  OUT_DVO(I,J,K)
             AVO(I,J,K)  =  OUT_AVO(I,J,K)
             WO(I,J,K)   =  OUT_WO(I,J,K)
           END DO
           ZO(I,J)      = OUT_ZO(I,J)
           Z1O(I,J)     = OUT_Z1O(I,J)
           SICTHO(I,J)  = OUT_SICTHO(I,J)
           SICOMO(I,J)  = OUT_SICOMO(I,J)
           SICUO(I,J)   = OUT_SICUO(I,J)
           SICVE(I,J)   = OUT_SICVE(I,J)
           HIBZETO(I,J) = OUT_HIBZETO(I,J)
           HIBETO(I,J)  = OUT_HIBETO(I,J)
           HIBZETE(I,J) = OUT_HIBZETE(I,J)
           HIBETE(I,J)  = OUT_HIBETE(I,J)
           SICSNO(I,J)  = OUT_SICSNO(I,J)
           IF (.not.L_COUPLING) TICE(I,J)    = OUT_TICE(I,J)
          END DO
        END DO
        !NOTE
        ! needed for rerun:
        IF (p_pe==p_io) write (*,*) "MPIOM RERUN INITIALIZATION"
#if defined (MPIOM_13B)
        CALL bounds_exch('u',UOO,'aufr 1')
        CALL bounds_exch('v',VOE,'aufr 2')
        !:: UPDATE VELOCITY FIELDS
        CALL OCTIMF
        CALL bounds_exch('p',THO,'aufr 3')
        CALL bounds_exch('p',SAO,'aufr 4')
        CALL bounds_exch('p',Z1O,'aufr 6')
        CALL bounds_exch('p',SICTHO,'aufr 7')
        CALL bounds_exch('p',SICOMO,'aufr 8')
        CALL bounds_exch('u',SICUO,'aufr 9')
        CALL bounds_exch('v',SICVE,'aufr 10')
        CALL bounds_exch('p',SICSNO,'aufr 11')
        CALL bounds_exch('s',hibete,'aufr 12')
        CALL bounds_exch('p',hibeto,'aufr 13')
        CALL bounds_exch('s',hibzete,'aufr 14')
        CALL bounds_exch('p',hibzeto,'aufr 15')
        CALL bounds_exch('p',dvo,'aufr 16')
        CALL bounds_exch('p',avo,'aufr 17')
        CALL bounds_exch('p',wo,'aufr 18')
#elif defined (MPIOM_2000)
        CALL bounds_exch(1,'u',UOO,'aufr 1')
        CALL bounds_exch(1,'v',VOE,'aufr 2')
        !:: UPDATE VELOCITY FIELDS
        CALL OCTIMF
        CALL bounds_exch(1,'p',THO,'aufr 3')
        CALL bounds_exch(1,'p',SAO,'aufr 4')
        CALL bounds_exch(1,'p',Z1O,'aufr 6')
        CALL bounds_exch(1,'p',SICTHO,'aufr 7')
        CALL bounds_exch(1,'p',SICOMO,'aufr 8')
        CALL bounds_exch(1,'u',SICUO,'aufr 9')
        CALL bounds_exch(1,'v',SICVE,'aufr 10')
        CALL bounds_exch(1,'p',SICSNO,'aufr 11')
        CALL bounds_exch(1,'s',hibete,'aufr 12')
        CALL bounds_exch(1,'p',hibeto,'aufr 13')
        CALL bounds_exch(1,'s',hibzete,'aufr 14')
        CALL bounds_exch(1,'p',hibzeto,'aufr 15')
        CALL bounds_exch(1,'p',dvo,'aufr 16')
        CALL bounds_exch(1,'p',avo,'aufr 17')
        CALL bounds_exch(1,'p',wo,'aufr 18')
#endif
        DO K=1,KEP
        DO J=1,JE
         DO I=1,IE
          DVO(I,J,1)=0.
          DVO(I,J,KEP)=0.
          AVO(I,J,1)=0.
          AVO(I,J,KEP)=0.
         ENDDO
        ENDDO
        ENDDO
#if defined (MPIOM_13B)
        CALL bounds_exch('p',tice,'aufr 19')
#elif defined (MPIOM_2000)
        CALL bounds_exch(1,'p',tice,'aufr 19')
#endif
        DO K=1,KE
        DO J=1,JE
        DO I=1,IE
          SAO(I,J,K)=MIN(SAO(I,J,K),70.)
          THO(I,J,K)=MIN(THO(I,J,K),70.)
        ENDDO
        ENDDO
        ENDDO
        !  CHECK GLOBAL SALT CONTENT
             IF ( ICONTRO .ne. 0 ) THEN
              CALL CONTRO(13)
             ENDIF

      IF (L_CHECK) CALL mpiom_check("start of rerun")

      !------------------------------------------------------------
      ! Tidal mode
      ! see mo_tidal.f90:
      ! write (*,*) "YEAR,MONTH",YEAR,MONTH
#if defined (MPIOM_13B)
      IF (L_TIDES) THEN
#elif defined (MPIOM_2000)
      IF (ltidal) THEN
#endif
         CALL foreph_ini
         mmccdt = NINT(OUT_mmccdt)
      ENDIF

      ENDIF !lresume


!    CALCULATION OF SOME VARIABLES NEEDED IN OTHER SUBMODELS OR LATER


     DO I=1,IE
      DO J=1,JE
       DO K=1,KE
         IF (K==1) THEN
            ! total volume (without ice and snow)
            OMVOL(I,J,K) = dlyp(I,J)*dlxp(I,J)* &
                           (zo(I,J)+ddpo(I,J,1) &
                           ! mz_bk_20110311+
                           -sictho(i,j)*rhoicwa-sicsno(i,j)*rhosnwa)
                           ! mz_bk_20110311-
         ELSE
            OMVOL(I,J,K) = dlyp(I,J)*dlxp(I,J)* &
                           ddpo(I,J,K)
         ENDIF
         OMMASS(I,J,K) = omvol(I,J,K)*rhoo(i,j,k)
       ENDDO
      ENDDO
     ENDDO

    CALL nullify_borders(status,OMMASS(:,:,:))
    CALL nullify_borders(status,OMVOL(:,:,:))

  END SUBROUTINE mpiom_global_start

! ---------------------------------------------------------------------------

  SUBROUTINE mpiom_global_end

#ifdef MPIOM_2000
  USE messy_main_mpi_bi,          ONLY: p_io, p_pe
#endif
  USE messy_main_grid_def_mem_bi, ONLY: jrow, ngpblks, kproma
  USE messy_main_timer,     ONLY: time_step_len &
                                , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND &
                                , lstart, lresume
  USE messy_mpiom_tools_e5, ONLY: nullify_borders
  USE messy_main_constants_mem, ONLY: iouerr

  IMPLICIT NONE

  INTEGER :: status
  INTEGER :: IM1,JM1


!#############################################################################
!                               TIME LOOP !!!!!
!#############################################################################

   IF (l_trig_mpiom) THEN !mpiom_time_step

    IF (L_COUPLING) THEN
      FU10(1:IE,1:JE)=AOFLWSVO(1:IE,1:JE)
    ENDIF

    IF (L_CHECK) CALL mpiom_check("beginning of time-step")

    IF (p_pe==p_io) THEN
      WRITE(*,'(a,i4.4,''-'',i2.2,''-'',i2.2,''    '',i2.2,'':'',i2.2)') &
           'MPIOM: current date: ', YEAR,MONTH,DAY, HOUR,MINUTE
    ENDIF

#if defined (MPIOM_13B)
    IF (L_TIDES) THEN
#elif defined (MPIOM_2000)
    IF (ltidal) THEN
#endif
       CALL foreph
    ENDIF


!--------------------------------------------------------------------
! HAMOCC COUPLING
#if defined (MPIOM_13B)
    IF (.NOT.L_HAMOCC_COUPLING) THEN
      opendep=11.
       DO k=1,ke
        DO j=1,je
          DO i=1,ie
          abs_oce(i,j,k) = EXP(-tiestw(k)/opendep)
          END DO
        END DO
       END DO
    ENDIF
#elif defined (MPIOM_2000)
!NB in MPIOM_2000 abs_oce => swr_frac (check mo_swr_absorption.f90)!!!
   IF (.NOT.L_HAMOCC_COUPLING) THEN
     opendep=11.
     DO k=1,ke
      DO j=1,je
        DO i=1,ie
          IF ( .NOT. lswr_jerlov ) THEN !old absorption
           abs_oce(i,j,k) = EXP(-tiestw(k)/opendep)
          ELSE                          !new absorption
           abs_oce(i,j,k) = jerlov_bluefrac*EXP(-tiestw(k)*jerlov_atten)
          ENDIF
        END DO
      END DO
     END DO
   ENDIF
#endif
    !     length scale for penetration of sw-radiation [m]
    !$OMP DO
    DO j=1,je
      DO i=1,ie
        swsum(i,j)=abs_oce(i,j,2)
!        swsumi(i,j)=1./abs_oce(i,j,2)
      END DO
    END DO
    DO k=1,ke-1
      DO j=1,je
        DO i=1,ie
!          swrab(i,j,k)=swsumi(i,j)*(abs_oce(i,j,k)-abs_oce(i,j,k+1))
          swrab(i,j,k)=(1./swsum(i,j))*(abs_oce(i,j,k)-abs_oce(i,j,k+1))
        END DO
      END DO
    END DO
    k=ke
    DO j=1,je
      DO i=1,ie
        swrab(i,j,k)=(1./swsum(i,j))*abs_oce(i,j,k)
      END DO
    END DO


!--------------------------------------------------------------------
! OCICE : THERMODYNAMIC FORCING

    IF (L_COUPLING) THEN
      DO I=1,IE
      DO J=1,JE
#if defined (MPIOM_13B)
      TXO(I,J) = AOFLTXIO(I,J)/rhowat
      TYE(I,J) = AOFLTYIO(I,J)/rhowat
#elif defined (MPIOM_2000)
      TXO(I,J) = AOFLTXIO(I,J)/rhoref_water
      TYE(I,J) = AOFLTYIO(I,J)/rhoref_water
#endif
!      IF (aoflrhio(I,J) < 0.0 .AND. weto(I,J,1) > 0.5) THEN
!         WRITE(*,*) 'Modify residual heat flux at :', I,J,aoflrhio(I,J),aoflchio(I,J)
!         aoflchio(I,J) = aoflchio(I,J) + aoflrhio(I,J)
!         aoflrhio(I,J) = 0.0
      ENDDO
      ENDDO
      CALL OCICE(AOFLFRIO,AOFLFRWO,AOFLRHIO,AOFLCHIO,AOFLNHWO,AOFLSHWO, &
               AOFLWSVO)

    ELSE
      CALL OCICE
    ENDIF


    IF (L_CHECK) CALL mpiom_check("after OCICE")

!--------------------------------------------------------------------
! OCTHER :  relaxation (some not needed in the uncoupled model)

    IF (L_COUPLING) THEN
         ! CALL OCTHER is equal to:
         !* CALL relax_surf
         !    surface relaxing term accordingly to
         !    RELTEM and RELSAL (see mpiom.nml)
         !    not needed in coupled mode
         !* CALL river_runoff
         !    HD present??? if not you need this call!
         !    it modifies SAO and ZO acordingly to input files
         !* (if hamocc) CALL dilcor_ptrf2
         !    update concentration based on new ZO values
         !    only in uncoupled model
#if defined (MPIOM_13B)
         CALL convection
         CALL calc_rinum
         CALL calc_dens
#elif defined (MPIOM_2000)
         CALL convection
         CALL calc_rinum
         IF  (ioconv .NE. 1 ) THEN
           CALL calc_dens
         ENDIF
#endif

    ELSE
         CALL OCTHER
    ENDIF

    IF (L_CHECK) CALL mpiom_check("after OCTHER")

    IF ( ICONTRO .ne. 0 ) THEN
       CALL CONTRO(1)
    ENDIF

!--------------------------------------------------------------------
! OCWIND :

    IF (L_COUPLING) THEN
      DO I=1,IE
      DO J=1,JE
#if defined (MPIOM_13B)
      TXO(I,J) = AOFLTXWO(I,J)/rhowat
      TYE(I,J) = AOFLTYWO(I,J)/rhowat
#elif defined (MPIOM_2000)
      TXO(I,J) = AOFLTXWO(I,J)/rhoref_water
      TYE(I,J) = AOFLTYWO(I,J)/rhoref_water
#endif
      ENDDO
      ENDDO
    ENDIF

#if defined (MPIOM_2000)
    IF (IBOLK .LT. 0) THEN
      IF (l_trig_mpiom_day) CALL CALCGMVIS
    ENDIF
#endif

    CALL OCWIND

    IF (L_CHECK) CALL mpiom_check("after OCWIND")

!--------------------------------------------------------------------
! TDES calculation
#if defined (MPIOM_13B)
    IF (L_TIDES) THEN
#elif defined (MPIOM_2000)
    IF (ltidal) THEN
#endif
       CALL tipouv
    ENDIF

    IF (L_CHECK) CALL mpiom_check("after tipouv")

!--------------------------------------------------------------------
!  OCTIMF : UPDATE VELOCITY

    CALL OCTIMF

    IF (L_CHECK) CALL mpiom_check("after OCTIMF")

!------------------------------------------------------------------
!  OCMODMOM : DECOMPOSITION INTO BAROTROPIC AND BAROCLINIC FIELD

    CALL OCMODMOM

    IF (L_CHECK) CALL mpiom_check("after OCMODMOM")

!------------------------------------------------------------------

#if defined (MPIOM_13B)
    CALL OCBARP
#elif defined (MPIOM_2000)
    IF ( lwith_barotropic_stokes_drift ) THEN
      CALL itprep
    ENDIF
    IF ( .NOT. lwith_barotropic_stokes_drift ) THEN
      CALL OCBARP
    ENDIF
#endif

    IF (L_CHECK) CALL mpiom_check("after OCBARP")

!------------------------------------------------------------------

#if defined (MPIOM_2000)
    IF ( lwith_barotropic_stokes_drift ) THEN
      CALL troneu
    ENDIF
#endif

    CALL OCCLIT

    IF (L_CHECK) CALL mpiom_check("after OCCLIT")
!------------------------------------------------------------------

    IF ( ICONTRO .ne. 0 ) THEN
        CALL contro(38)
    ENDIF
#if defined (MPIOM_13B)
    IF(IELIMI.GE.1) THEN
      CALL BARTIM
    ELSE
      CALL TRONEU
    ENDIF

#elif defined (MPIOM_2000)
    IF(IELIMI.GE.1) THEN
      CALL BARTIM
    ELSE
      IF ( .NOT. lwith_barotropic_stokes_drift ) THEN
!TRONEU2 NOT WORKING!!!!
                 CALL TRONEU
         !CALL TRONEU2(IHALO_SOR,IMOD)
      ENDIF
    ENDIF
#endif

    IF ( ICONTRO .ne. 0 ) THEN
      CALL contro(39)
    ENDIF

    IF (L_CHECK) CALL mpiom_check("after TRONEU")

!------------------------------------------------------------------

    CALL OCVTRO

    IF (L_CHECK) CALL mpiom_check("after OCVTRO")

#if defined (MPIOM_13B)
    ! done after advection in MPIOM_2000
    CALL update_zo
#endif
    ! TO DO : HERE TRACER DILUTION!

    IF (L_CHECK) CALL mpiom_check("after update_zo")

!--------------------------------------------------------------------

    IF ( ICONTRO .ne. 0 ) THEN
      CALL contro(40)
    ENDIf
    CALL OCVTOT
!   PRINT*,'nach ocvtot'
    IF ( ICONTRO .ne. 0 ) THEN
      CALL contro(41)
    ENDIF

    !UWE   RESET STABIO TO DRHODZ
    DO K=1,KE
     DO J=1,JE
      DO I=1,IE
       STABIO(I,J,K)=(1000./DZ(K))*STABIO(I,J,K)
      ENDDO
     ENDDO
    ENDDO
    IF ( ICONTRO .ne. 0 ) THEN
      CALL contro(42)
    ENDIF

#ifdef MPIOM_13B
    IF ( ICONTRO .ne. 0 ) THEN
       call contro(43)
    ENDIF
    CALL OCUAD(UOO)
    CALL OCVAD(VOE)
    IF ( ICONTRO .ne. 0 ) THEN
       call contro(44)
    ENDIF
#elif defined (MPIOM_2000)
    CALL calc_icecutoff  ! this should go into mo_ocice

    CALL cell_thickness  ! update cell thickness

    IF ( lundelayed_momentum_advection ) THEN
      CALL ocuad(uoo)
      CALL ocvad(voe)
      uaccel(:,:,:) = 0._wp
      vaccel(:,:,:) = 0._wp
    ELSE
      ! new "delayed" momentum advection
      IF ( lbounds_exch_tp ) THEN
        !maybe not needed
        CALL bounds_exch(1,'u',uoo,'mpiom 29')
        CALL bounds_exch(1,'v',voe,'mpiom 29')

        IF ( have_g_js ) THEN
          uoo(:,1,:) = 0._wp
        ENDIF
      ENDIF

      uaccel(:,:,:)=uoo(:,:,:)
      vaccel(:,:,:)=voe(:,:,:)

      CALL ocuad(uaccel)
      CALL ocvad(vaccel)

      uaccel(:,:,:)=uaccel(:,:,:)-uoo(:,:,:)
      vaccel(:,:,:)=vaccel(:,:,:)-voe(:,:,:)

      IF ( ICONTRO .NE. 0 ) THEN
        WRITE(iouerr,*) 'uaccel nach ocuad',MAXVAL(uaccel),MAXLOC(uaccel)
        WRITE(iouerr,*) 'vaccle nach ocvad',MAXVAL(vaccel),MAXLOC(vaccel)
        CALL contro(44)
      END IF

    ENDIF

#endif

!--------------------------------------------------------------------
!  ADVECTION

#ifdef MPIOM_13B
    IF(IOCAD.EQ.4)THEN               ! SLOPECON_ADPO
       CALL SLOPETRANS
       IF ( ICONTRO .ne. 0 ) THEN
          CALL CONTRO(2)
       ENDIF
    ENDIF

    IF(IOCAD.EQ.3.OR.IOCAD.EQ.4)THEN
       CALL OCADPO(THO)
       ! SALT ADVECTION
       CALL OCADPO(SAO)
    ENDIF

    IF(IOCAD.EQ.5) CALL OCADFS

#elif defined (MPIOM_2000)
    IF(ibbl_transport.EQ.1)THEN               ! slope convection
      CALL SLOPETRANS
      IF ( ICONTRO .NE. 0 ) THEN
        CALL CONTRO(2)
      END IF
    END IF

    IF (ladfs) CALL OCADFS

    !tracer advection
    CALL bounds_exch(1,'p',rhoo,'mpiom 29')

    !  Preparation for advection
    IF (ladpo) CALL ocadpo_base

    IF (ladpo) CALL ocadpo_trf(tho)
    CALL bounds_exch(1,'p',tho,'mpiom 30')
    IF (ladpo) CALL ocadpo_trf(sao)
    CALL bounds_exch(1,'p',sao,'mpiom 31')

    !--------------------------------------------------------------------
    !  UPDATE_ZO : UPDATE OF SEALEVEL
    CALL update_zo
#endif

#if defined (MPIOM_13B)
!--------------------------------------------------------------------
! OCJITR : ISOPYCNAL DIFFUSION (W/O GM PARAMETERISATION)

    CALL OCJITR

    IF ( ICONTRO .ne. 0 ) THEN
       CALL CONTRO(2713)
    ENDIF

!--------------------------------------------------------------------
!  TEMPERATURE AND SALINITY  DIFFUSION

    CALL bounds_exch('p',rhoo,'mpiom 29')
    CALL bounds_exch('p',tho,'mpiom 30')
    CALL bounds_exch('p',sao,'mpiom 31')

    DO IMAL=1,1
       CALL OCTDIFF_BASE
       CALL OCTDIFF_TRF(tho)
       CALL OCTDIFF_TRF(sao)
    ENDDO

#elif defined (MPIOM_2000)
!--------------------------------------------------------------------
! OCJITR : ISOPYCNAL DIFFUSION (W/O GM PARAMETERISATION)

     IF ( IBOLK .NE. 0 ) CALL ocjitr_base
     IF ( IBOLK .NE. 0 ) CALL ocjitr_trf(tho)
     IF ( IBOLK .NE. 0 ) CALL ocjitr_trf(sao)

!--------------------------------------------------------------------
!  TEMPERATURE AND SALINITY  DIFFUSION

    CALL octdiff_base
    CALL octdiff_trf(tho)
    CALL octdiff_trf(sao)

#endif

#ifdef MPIOM_2000
    ! mz_bk_20110315+
    CALL calc_global_mean
    ! mz_bk_20110315-
#endif

!--------------------------------------------------------------------
!  SALINITY and TEMPERATURE RELAXATION tau=(1/(4*months))
!  only if I3DREST GT 0
    CALL RELAX_TS
    IF ( ICONTRO .ne. 0 ) THEN
       CALL CONTRO(4)
    ENDIF

#ifdef MPIOM_2000
!NUDGING OF MPIOM!!!!
!    CALL READ_ECCO
!    CALL NUDGE_T(tho, rreltem)
!    CALL NUDGE_S(sao, rrelsal)
#endif

!--------------------------------------------------------------------
!  OCEAN VELOCITY UPDATE

    CALL OCTIMF
    IF ( ICONTRO .ne. 0 ) THEN
       CALL CONTRO(5)
    ENDIF

!--------------------------------------------------------------------
!  MOMENTUM DIFFUSION

    IF (AUS .GT. 1.E-12) CALL OCSCHEP

!--------------------------------------------------------------------
!  BIHARMONIC MOMENTUM DIFFUSION, BOTTOM FRICTION,
!  SHEAR DEPENDENT DIFFUSION

    CALL OCVISC

    IF ( ICONTRO .ne. 0 ) THEN
       CALL CONTRO(6)
    ENDIF

    CALL OCTIMF


!--------------------------------------------------------------------
!     correct ZO to a global mean of zero (once per day in the original code)
#if defined (MPIOM_13B)
      call correct_zo
#elif defined (MPIOM_2000)
      IF (l_trig_mpiom_day.or.lstart.or.lresume) THEN
       call correct_zo
      ENDIF
#endif

#if defined (MPIOM_13B)
    IF (L_TIDES) THEN
#elif defined (MPIOM_2000)
    IF (ltidal) THEN
#endif
      OUT_mmccdt = REAL(mmccdt,DP)
    ENDIF

    IF (L_CHECK) CALL mpiom_check("end of time-step")


!--------------------------------------------------------------------
!    FINISH WITH  OUTPUT ASSOCIATION

    ! mz_bk_20110315+
#if defined (MPIOM_2000)
    IF (ldiag) THEN
       out_global_sum_temperature = global_sum_temperature
       out_global_sum_salinity = global_sum_salinity
       out_global_mass = global_mass
       out_global_volume = global_volume
       out_global_salt_content = global_salt_content
    ENDIF
#endif
    ! mz_bk_20110315-

    DO I=1,IE
      DO J=1,JE
       DO k=1,KE
         OUT_RHOO(I,J,K) =  RHOO(I,J,K)
         OUT_THO(I,J,K)  =  THO(I,J,K)
         OUT_SAO(I,J,K)  =  SAO(I,J,K)
         OUT_PO(I,J,K)   =  PO(I,J,K)
         OUT_UOO(I,J,K)  =  UOO(I,J,K)
         OUT_VOE(I,J,K)  =  VOE(I,J,K)
         ! op_mk_20180108+
         OUT_UKO(I,J,K)  = UKO(I,J,K)
         OUT_VKE(I,J,K)  = VKE(I,J,K)
         OUT_DDUO(I,J,K) = DDUO(I,J,K)
         OUT_DDUE(I,J,K) = DDUE(I,J,K)
         OUT_STABIO(I,J,K) = STABIO(I,J,K)
         ! op_mk_20180108-
         ! needed for forcing atmosphere
         OUT_THO_K(I,J,K)=  (THO(I,J,K)+tmelt)
         ! general output
         ! mz_bk_20111209+
         OUT_DDPO(I,J,K) = DDPO(I,J,K)
         ! mz_bk_20111209-
         IF (K==1) THEN
            ! total volume (without ice and snow)
            OMVOL(I,J,K) = dlyp(I,J)*dlxp(I,J)* &
                           (zo(I,J)+ddpo(I,J,1) &
                           ! mz_bk_20110311+
                           -sictho(i,j)*rhoicwa-sicsno(i,j)*rhosnwa)
                           ! mz_bk_20110311-
         ELSE
            OMVOL(I,J,K) = dlyp(I,J)*dlxp(I,J)* &
                           ddpo(I,J,K)
         ENDIF
         OMMASS(I,J,K) = omvol(I,J,K)*rhoo(i,j,k)
       END DO
       DO k=1,KE+1
         OUT_DVO(I,J,K)  =  DVO(I,J,K)
         OUT_AVO(I,J,K)  =  AVO(I,J,K)
         OUT_WO(I,J,K)   =  WO(I,J,K)
       END DO
       OUT_ZO(I,J)      = ZO(I,J)
       OUT_Z1O(I,J)     = Z1O(I,J)
       OUT_SICTHO(I,J)  = SICTHO(I,J)
       OUT_SICOMO(I,J)  = SICOMO(I,J)
       OUT_SICUO(I,J)   = SICUO(I,J)
       OUT_SICVE(I,J)   = SICVE(I,J)
       OUT_HIBZETO(I,J) = HIBZETO(I,J)
       OUT_HIBETO(I,J)  = HIBETO(I,J)
       OUT_HIBZETE(I,J) = HIBZETE(I,J)
       OUT_HIBETE(I,J)  = HIBETE(I,J)
       OUT_SICSNO(I,J)  = SICSNO(I,J)
       IF (.not.L_COUPLING) OUT_TICE(I,J)    = TICE(I,J)
       IF (.not.L_COUPLING) THEN
        OUT_EMINPO(I,J)  = EMINPO(I,J)
        OUT_RIVRUN(I,J) = RIVRUN(I,J)
       ENDIF
      END DO
    END DO

    CALL nullify_borders(status,OMMASS(:,:,:))
    CALL nullify_borders(status,OMVOL(:,:,:))

    ! needed for forcing atmosphere
    DO I=2,IE-1
      DO J=2,JE-1
       IM1=I-1
       JM1=J-1
       IF (JM1 == 0) JM1=1
       IF (IM1 == 0) IM1=1
       OUT_SOCU(I,J) = (1.-SICOMO(I,J))*(UOO(I,J,1)+UOO(IM1,J,1))*0.5+ &
                           SICOMO(I,J)*(SICUO(I,J)+SICUO(IM1,J))*0.5
       OUT_SOCV(I,J) = (1.-SICOMO(I,J))*(VOE(I,J,1)+VOE(I,JM1,1))*0.5 &
                          +SICOMO(I,J)*(SICVE(I,J)+SICVE(I,JM1))*0.5
      END DO
    END DO
#if defined (MPIOM_13B)
    CALL bounds_exch('p',OUT_SOCU  ,'start 3')
    CALL bounds_exch('p',OUT_SOCV  ,'start 4')
#elif defined (MPIOM_2000)
    CALL bounds_exch(1,'p',OUT_SOCU  ,'start 3')
    CALL bounds_exch(1,'p',OUT_SOCV  ,'start 4')
#endif

!    CALL bounds_exch('p',OUT_RHOO   ,'output 1')
!    CALL bounds_exch('u',OUT_UKO    ,'output 2')
!    CALL bounds_exch('v',OUT_VKE    ,'output 3')
!    CALL bounds_exch('p',OUT_THO    ,'output 4')
!    CALL bounds_exch('p',OUT_SAO    ,'output 5')
!    CALL bounds_exch('p',OUT_WO     ,'output 6')
!    CALL bounds_exch('p',OUT_PO     ,'output 7')
!    CALL bounds_exch('u',OUT_UOO    ,'output 8')
!    CALL bounds_exch('v',OUT_VOE    ,'output 9')
!    CALL bounds_exch('p',OUT_DVO    ,'output 10')
!    CALL bounds_exch('p',OUT_AVO    ,'output 11')
!    CALL bounds_exch('p',OUT_THO_K  ,'output 12')
!    CALL bounds_exch('p',OUT_ZO     ,'output 13')
!    CALL bounds_exch('p',OUT_Z1O    ,'output 14')
!    CALL bounds_exch('p',OUT_SICTHO ,'output 15')
!    CALL bounds_exch('p',OUT_TICE   ,'output 16')
!    CALL bounds_exch('p',OUT_SICOMO ,'output 17')
!    CALL bounds_exch('u',OUT_SICUO  ,'output 18')
!    CALL bounds_exch('v',OUT_SICVE  ,'output 19')
!    CALL bounds_exch('p',OUT_HIBZETO,'output 20')
!    CALL bounds_exch('p',OUT_HIBETO ,'output 21')
!    CALL bounds_exch('s',OUT_HIBZETE,'output 22')
!    CALL bounds_exch('s',OUT_HIBETE ,'output 23')
!    CALL bounds_exch('p',OUT_SICSNO ,'output 24')


!#############################################################################
!                             END OF TIME LOOP !!!!!
!#############################################################################

   END IF !mpiom_time_step

  END SUBROUTINE mpiom_global_end

  ! ---------------------------------------------------------------------------

  SUBROUTINE mpiom_free_memory

    IMPLICIT NONE
    INTRINSIC ALLOCATED, ASSOCIATED

    IF (.not.L_HAMOCC_COUPLING) THEN
        DEALLOCATE(abs_oce)
    ENDIF


  END SUBROUTINE mpiom_free_memory

  ! ---------------------------------------------------------------------------

! ===========================================================================
! PRIVATE MPIOM INTERFACE ROUTINES
! ===========================================================================

  SUBROUTINE mpiom_read_nml_ocectl(status, iou)
    ! MESSy
    USE messy_main_tools,     ONLY: read_nml_open, read_nml_check, read_nml_close
    USE messy_main_timer,     ONLY: lresume

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

#if defined (MPIOM_13B)
    NAMELIST /ocectl/ &
              CAULAPTS, CAULAPUV, CAH00, AUS                           &
             ,AV0,DV0,CWT,CSTABEPS,DBACK,ABACK,CRELSAL,CRELTEM         &
             ,ICONVA, IWIRB                                            &
             ,CDVOCON,CAVOCON                                          &
             ,IOCAD,isnflg                                             &
             ,H0, HMIN, ARMIN, ARMAX, HSNTOICE, SICTHMIN, SICE         &
             ,D3                                                       &
             ,ISTART,I3DREST                                           &
             ,ICONTRO
#elif defined (MPIOM_2000)
    NAMELIST /ocectl/                                                  &
              CAULAPTS, CAULAPUV, CAH00, AUS, bofric, rayfric          &
             ,AV0,DV0,CWT,CWA,CSTABEPS,DBACK,ABACK,CRELSAL,CRELTEM     &
             ,ltidal,ltidal_diag,lswr_jerlov,jerlov_atten,jerlov_bluefrac &
             ,lwith_one_layer_shelfs,lfb_bgc_oce                               &
             ,SPONGEZONE, spzndamp_time, RRELTEM, RRELSAL              &
             ,CDVOCON,CAVOCON,IBOLK,lisopyc                            &
             ,IOCAD,iocaduv,ioconv,isnflg,numriv,numglac               &
             ,H0, HMIN, ARMIN, ARMAX, HSNTOICE, SICTHMIN, SICE         &
             ,rleadclose                                               &
             ,IAUFR, IAUFW                                             &
             ,ISTART,I3DREST,nfixYearLen                               &
             ,ICONTRO,IOASISFLUX                                       &
             ,LFORCEDIAG,LHFLDIAG,LCONVDIAG,LDIFFDIAG,LGMDIAG          &
             ,LGRIDINFO,LCALCDIFI,lmpitype,lnonblock                   &
!             ,ITSDIAG,LTSTRANSPOSE,ihalo_sor,imod                      &
             ,ihalo_sor,imod                                           &
             ,fp_tracing_enabled,ibbl_transport, lsaoclose             &
             ,iter_sor, rtsorpar, iter_sor_hack,rtsorpar_hack          &
             ,lundelayed_momentum_advection                            &
             ,lwith_barotropic_stokes_drift                            &
             ! mz_bk_20110317+
             ,ldiag
             ! mz_bk_20110317-

#endif

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'mpiom_read_nml_ocectl'
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1
!----------------------------------------------------------------------
! DEFAULT PARAMETER SETTINGS - CAN BE OVERWRITEN BY THE NAMELIST
!----------------------------------------------------------------------

#if defined (MPIOM_13B)
!  TRACER ADVECTION ROUTINES : SEVERAL OPTIONS
!  IOCAD ==  1: UPWIND
!  IOCAD ==  2: not used
!  IOCAD ==  3: ADPO
!  IOCAD ==  4: ADPO + SLOPECON_ADPO
!  IOCAD ==  5: ADFS
!  IOCAD ==  6: not yet --> QUICK
!  IOCAD ==  7: mot yet --> QUICK2
        IOCAD=4

        ICONTRO=0

!UWE   CONSTANTS FOR BIHARMONIC DIFFUSION
         CAULAPTS=0.0002

!UWE   CONSTANTS FOR HARMONIC DIFFUSION
         CAH00=0.

#ifdef ISOPYK
         CAULAPTS=0.
#endif

#ifdef ISOPYK
         CAH00=1000.
#endif

!HH    CONSTANT FOR WINDMIXING IN OCTHER
        CWT=0.5E-3
        CSTABEPS=0.05

!HH    BACKGROUNDDIFFUSION
        ABACK=5.E-5
        DBACK=5.E-5

!----------------------------------------------------------------------
!     ICONVA : 1   WITH CONVECTIVE ADJUSTMENT
!              0    NO      "          "
      ICONVA=1
!-------------------------------------------------------
!     IWIRB  : 1   COMPUTE VARIABLE EDDY VISCOSITY
!              0   CONSTANT EDDY VISCOSITY
      IWIRB=1
!----------------------------------------------------------------------
! SWITCH FOR RESPECTIVE RESTART FILES
!
      IFLAG=1

!----------------------------------------------------------------------
! PARAMETER ICE MODEL

      ISNFLG=1


      H0=0.5
      HMIN=0.05
      ARMIN=0.15
      ARMAX=1.0

!UWE  MAXIMUM ICE THICKNESS FOR SNOW-->ICE CONVERSION
!      HSNTOICE = 17. !to be corrected! mz_ap
     HSNTOICE = 0.45 * CDZW(1)
!UWE  MINIMUM ICE THICKNESS IN NEW ICE GROWTH
      SICTHMIN=0.5
      SICE=5.0
      D3=5.5E-08
! ISTART :: START OPTIONS
! ISTART == 0: COMPLETELY NEW SETUP,
!              topography read from anta and written to topo (!!! WARNING !!!)
!              start from climatology
! ISTART == 1: new run, topography read from topo
!              start from horizontally uniform ts-profile (taf,saf)
! ISTART == 2: new run, topography read from topo
!              start from climatology
! ISTART == 3: continuing run (default)
         ISTART =3

#elif defined (MPIOM_2000)

    ! I3DREST OPTIONS FOR 3-D RESTORING
    ! I3DREST == 0: NO RESTORINg (DEFAULT)
    ! I3DREST == 1: RESTORING to annual climatology
    ! I3DREST == 2: RESTORING to monthly climatology
    I3DREST = 0
    !--------------------------------------------------------------------
    !  TRACER ADVECTION ROUTINES : SEVERAL OPTIONS
    !  IOCAD ==  1: UPWIND SCHEME
    !  IOCAD ==  2: CENTRAL-DIFFERENCE SCHEME (not tested)
    !  IOCAD ==  3: ADPO   (TVD scheme after Sweby, 1984)
    !  IOCAD ==  4: obsolete, use IOCAD=3 and IBBL_TRANSPORT=1 for ADPO + SLOPECON_ADPO
    !  IOCAD ==  5: ADFS   (predictor-corrector advection scheme)
    !  IOCAD ==  6: QUICK  (quick-scheme after Farrow and Stevens, 1995)
    !  IOCAD ==  7: QUICK2 (quick-scheme with modified boundary treatment)
    !  IOCAD ==  8: as IOCAD=3 but with splitted vertical tracer transport
    !               (may improve numerical stability for high resolution setups)
    IOCAD = 3
    !
    !  IMPULSE ADVECTION ROUTINES
    !  IOCADUV = 1 UPWIND SCHEME (NOT IMPLEMENTED YET!!!)
    !  IOCADUV = 2 CENTRAL-DIFFERENCE SCHEME (NOT IMPLEMENTED YET!!!)
    !  IOCADUV = 3 ADPO (TVD scheme after Sweby, 1984)
    !  IOCADUV = 8 as IOCAD=3 but with splitted vertical impulse transport
    !               (may improve numerical stability for high resolution setups)
    !
    !  ATTENTION: valid iocaduv values are 3 (default) and 8 currently
    !
    IOCADUV = 3
    !
    !--------------------------------------------------------------------
    !  BOTTOM BOUNDARY LAYER ADVECTIVE TRANSPORT SCHEME (SLOPE CONVECTION)
    !  IBBL_TRANSPORT = 1: BBL advective transport to neutral buoyancy level
    !                      can not be used with IOCAD=5,6,7

    IBBL_TRANSPORT = 1

    !--------------------------------------------------------------------
    !! Parametrization schemes for convective adjustment (Attion: only option 1 works with HAMOCC).
    !! ioconv = 1 : enhancement of vertical diffusivity (default)
    !! ioconv = 2 : complete mixing
    !! ioconv = 3 : interchange of upper and lower box
    !! ioconv = 4 : plume convection model

    IOCONV = 1

    ICONTRO = 0

    !   CONSTANTS FOR BIHARMONIC DIFFUSION
    caulapts = 0.0_wp
    !   CONSTANTS FOR HARMONIC DIFFUSION
    cah00 = 1000._wp

    !  SWITCH FOR COMMUNICATION WITH MPI DATA TYPES
    lmpitype = .FALSE.
    lnonblock = .FALSE.

    !  SWITCH TO ENABLE TIDAL SUBMODEL
    ltidal = .FALSE.
    ltidal_diag = .FALSE.

    !  SWITCH TO ENABLE JERLOV_SWR (default is Jerlov Type II water)
    lswr_jerlov = .true.
    jerlov_atten = 0.12_wp
    jerlov_bluefrac = 0.28_wp


    !  SWITCH TO ENABLE ISOPYCNAL DIFFUSION
    lisopyc = .TRUE.

    !SOR DEFAULT : CALCULATE SOR WITH A HALO OF TWO; IMOD SHOULD BE SET TO IHALO/2
    ihalo_sor = 2
    imod = 1

    !BOFRIC  FRICTION COEFFICENT AT BOTTOM [M**2/S]
    bofric = 1.E-3_wp
    !COEFFICIENT FOR QUADRATIC BOTTOM FRICTION []
    rayfric = 3.E-3_wp

    IBOLK = 1000

    !HH    CONSTANT FOR WINDMIXING IN OCTHER
    cwt = 0.5E-3_wp
    cwa = -9e33_wp   ! cwa is set to cwt if the special value
    !  is not overwritten in the namelist
    cstabeps = 0.05_wp

    !HH    BACKGROUNDDIFFUSION
    aback = 5.E-5_wp
    dback = 5.E-5_wp

    !HH    DEFAULT MEAN OUTPUT
    !IMEAN = 2
!mz_ap_20100923 NO OUTPUT!
    IMEAN = 0

    !HH    SPONGE ZONE RELAXATION TIME
    spzndamp_time = 1._wp /(86400._wp * 30._wp * 4._wp)   ! default is 4 month

    ! LEADCLOSING parameters (default) can be overwritten in the namelist
    rleadclose(1) = 0.5_wp
    rleadclose(2) = 5._wp
    rleadclose(3) = 4._wp

    ! sea ice parameters (namelist)
    h0 = 0.5_wp
    armin = 0.15_wp
    armax = 1.0_wp
    hmin = 0.05_wp
    !UWE  MAXIMUM ICE THICKNESS FOR SNOW-->ICE CONVERSION
    hsntoice = 0.45_wp ! percentage of the upper layer
    !UWE  MINIMUM ICE THICKNESS IN NEW ICE GROWTH
    sicthmin = 0.5_wp

    ! salinity of sea ice
    sice = 5.0_wp

    ! number of observed river runoff stations (if used OMIP rivers are disabled)
    numriv=0
    luse_river_runoff_stations=.false.

    ! number of glacier calcing stations
    numglac=0
    luse_glac_calv=.false.


    lwith_barotropic_stokes_drift = .FALSE.

    ! mz_bk_20110317+
    ldiag = .FALSE.
    ! mz_bk_20110317-
    !  SWITCH FOR WRITING OASIS COUPLED FLUXES

    IOASISFLUX = 0


    LGMDIAG    = .FALSE.
    LFORCEDIAG = .FALSE.
    LHFLDIAG   = .FALSE.
    LCONVDIAG  = .FALSE.
    LDIFFDIAG  = .FALSE.
    LCALCDIFI  = .FALSE.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!FOR I/O not really used!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ltstranspose = .TRUE.
    ltswrite = .TRUE.     !  set to .FALSE. after first call of diagnosis

#endif

!UWE   CONSTANTS FOR BIHARMONIC FRICTION
         CAULAPUV=0.0045

!UWE   CONSTANTS FOR HARMONIC FRICTION
         AUS=3.E-6

!HH    CONSTANT FOR DIFFUSION IN OCTHER
        DV0=0.5E-2
        AV0=0.5E-2
        CDVOCON=20.
        CAVOCON=0.

!HH    RELAXATION TIME SALINITY
        CRELSAL = 3.E-7

!HH    RELAXATION TIME TEMERATURE
        CRELTEM = 0.
!
! I3DREST OPTIONS FOR 3-D RESTORING
! I3DREST == 0: NO RESTORINg (DEFAULT)
! I3DREST == 1: RESTORING to annual climatology
! I3DREST == 2: RESTORING to monthly climatology
         I3DREST=0

!--------------------------------------------------------------------

    CALL read_nml_open(lex, substr, iou, 'OCECTL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=OCECTL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'OCECTL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

    IF (lresume) ISTART=3

  END SUBROUTINE mpiom_read_nml_ocectl

!-------------------------------------------------------------------------------

  SUBROUTINE mpiom_check(text)

    USE messy_main_mpi_bi,    ONLY: p_parallel_io, p_io, p_pe, p_bcast !!$, finish

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: text


      IF (p_pe==p_io) THEN
        WRITE(*,*) "MPIOM at: "
        WRITE(*,*) text
        WRITE(*,*) "------------------------------------------"
        WRITE(*,*) "MAXVAL(UOO    )->",MAXVAL(UOO    *WETO)
        WRITE(*,*) "MAXVAL(VOE    )->",MAXVAL(VOE    *WETO)
        WRITE(*,*) "MAXVAL(THO    )->",MAXVAL(THO    *WETO)
        WRITE(*,*) "MAXVAL(SAO    )->",MAXVAL(SAO    *WETO)
        WRITE(*,*) "MAXVAL(ZO     )->",MAXVAL(ZO     *WETO(:,:,1))
        WRITE(*,*) "MAXVAL(Z1O    )->",MAXVAL(Z1O    *WETO(:,:,1))
        WRITE(*,*) "MAXVAL(SICTHO )->",MAXVAL(SICTHO *WETO(:,:,1))
        WRITE(*,*) "MAXVAL(SICOMO )->",MAXVAL(SICOMO *WETO(:,:,1))
        WRITE(*,*) "MAXVAL(SICUO  )->",MAXVAL(SICUO  *WETO(:,:,1))
        WRITE(*,*) "MAXVAL(SICVE  )->",MAXVAL(SICVE  *WETO(:,:,1))
        WRITE(*,*) "MAXVAL(SICSNO )->",MAXVAL(SICSNO *WETO(:,:,1))
        WRITE(*,*) "MAXVAL(hibete )->",MAXVAL(hibete *WETO(:,:,1))
        WRITE(*,*) "MAXVAL(hibeto )->",MAXVAL(hibeto *WETO(:,:,1))
        WRITE(*,*) "MAXVAL(hibzete)->",MAXVAL(hibzete*WETO(:,:,1))
        WRITE(*,*) "MAXVAL(hibzeto)->",MAXVAL(hibzeto*WETO(:,:,1))
        WRITE(*,*) "MAXVAL(dvo    )->",MAXVAL(dvo    )
        WRITE(*,*) "MAXVAL(avo    )->",MAXVAL(avo    )
        WRITE(*,*) "MAXVAL(wo     )->",MAXVAL(wo     )
        WRITE(*,*) "------------------------------------------"
!        WRITE(*,*) "MINVAL(UOO    )->",MINVAL(UOO    *WETO)
!        WRITE(*,*) "MINVAL(VOE    )->",MINVAL(VOE    *WETO)
!        WRITE(*,*) "MINVAL(THO    )->",MINVAL(THO    *WETO)
!        WRITE(*,*) "MINVAL(SAO    )->",MINVAL(SAO    *WETO)
!        WRITE(*,*) "MINVAL(ZO     )->",MINVAL(ZO     *WETO(:,:,1))
!        WRITE(*,*) "MINVAL(Z1O    )->",MINVAL(Z1O    *WETO(:,:,1))
!        WRITE(*,*) "MINVAL(SICTHO )->",MINVAL(SICTHO *WETO(:,:,1))
!        WRITE(*,*) "MINVAL(SICOMO )->",MINVAL(SICOMO *WETO(:,:,1))
!        WRITE(*,*) "MINVAL(SICUO  )->",MINVAL(SICUO  *WETO(:,:,1))
!        WRITE(*,*) "MINVAL(SICVE  )->",MINVAL(SICVE  *WETO(:,:,1))
!        WRITE(*,*) "MINVAL(SICSNO )->",MINVAL(SICSNO *WETO(:,:,1))
!        WRITE(*,*) "MINVAL(hibete )->",MINVAL(hibete *WETO(:,:,1))
!        WRITE(*,*) "MINVAL(hibeto )->",MINVAL(hibeto *WETO(:,:,1))
!        WRITE(*,*) "MINVAL(hibzete)->",MINVAL(hibzete*WETO(:,:,1))
!        WRITE(*,*) "MINVAL(hibzeto)->",MINVAL(hibzeto*WETO(:,:,1))
!        WRITE(*,*) "MINVAL(dvo    )->",MINVAL(dvo    *WETO)
!        WRITE(*,*) "MINVAL(avo    )->",MINVAL(avo    *WETO)
!        WRITE(*,*) "MINVAL(wo     )->",MINVAL(wo     *WETO)
!        WRITE(*,*) "------------------------------------------"
      ENDIF

  END SUBROUTINE mpiom_check

  ! um_ak_20130620+ ! mz_ap_20130620+
  SUBROUTINE define_mpiom_grid

!#ifdef FALSE
  ! MESSy/BMIL
!!$  USE messy_main_blather_bi,   ONLY: error_bi ! op_pj_20161104 (global USE)

  ! MESSy/SMCL
  USE messy_main_timer,        ONLY: YEAR_START, MONTH_START, DAY_START &
                                   , HOUR_START, MINUTE_START, SECOND_START &
                                   , YEAR, MONTH, DAY &
                                   , HOUR, MINUTE, SECOND

  ! MESSy
  USE messy_main_constants_mem, ONLY: api=>pi, sp
  USE messy_main_timer,         ONLY: time_span_d

  ! GRID MODULES
  USE MESSY_MAIN_GRID_TRAFO, ONLY: RGEMPTY, COMPLETE_GEOHYBGRID
  USE MESSY_MAIN_GRID_NETCDF, ONLY: ERRMSG  &
                                  , NULL_DIMID, NULL_VARID &
                                  , VTYPE_REAL, VTYPE_DOUBLE &
                                  , NF90_FLOAT, NF90_DOUBLE  &
                                  , POSITION, INIT_NARRAY, ADD_NCATT
  USE MESSY_MAIN_GRID,        ONLY: T_GEOHYBGRID,   INIT_GEOHYBGRID &
                                  , new_geohybgrid, locate_geohybgrid &
                                  , grid_error

  IMPLICIT NONE

  !LOCAL
  ! Global mpiom grid
  TYPE(t_geohybgrid) :: ggrid
  INTEGER            :: GGRID_ID
  ! Local mpiom grid
  TYPE(t_geohybgrid) :: lgrid
  INTEGER            :: lGRID_ID

  REAL(DP) :: dts
  INTEGER  :: i,j, n
  INTEGER  :: status
  CHARACTER(LEN=100)       :: tunit
  INTEGER  :: idx, jdx
  REAL(dp) :: gloni(ie_g+1,je_g+1)
  REAL(dp) :: glati(ie_g+1,je_g+1)

  ! --------------------------------
  ! DEFINE GLOBAL MPIOM
  ! --------------------------------

  ! INIT
  CALL INIT_GEOHYBGRID(ggrid)

  ggrid%name = 'MPIOMg'

  ggrid%file = 'MPIOMg-_GEO-HYBRID-GRID'       ! Filename
  ggrid%t    = 0                              ! time step

  ggrid%clonc = .TRUE. ! Longitude Axis is circular
  ! Curvilinear LONGITUDE (MID) ...
  ggrid%clonm%name  = 'lon'
  ggrid%clonm%id    = NULL_VARID
  ggrid%clonm%xtype = NF90_DOUBLE
  ! ... dimensions
  ggrid%clonm%ndims = 2
  ALLOCATE(ggrid%clonm%dim(ggrid%clonm%ndims), STAT=status)
  CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
  ggrid%clonm%dim(1)%name  = 'lon'
  ggrid%clonm%dim(1)%id    = NULL_DIMID
  ggrid%clonm%dim(1)%len   = ie_g
  ggrid%clonm%dim(1)%fuid  = .false.
  ggrid%clonm%dim(1)%varid = NULL_VARID
  ggrid%clonm%dim(2)%name  = 'lat'
  ggrid%clonm%dim(2)%id    = NULL_DIMID
  ggrid%clonm%dim(2)%len   = je_g
  ggrid%clonm%dim(2)%fuid  = .false.
  ggrid%clonm%dim(2)%varid = NULL_VARID
  ! ... data
  CALL INIT_NARRAY(ggrid%clonm%dat, ggrid%clonm%ndims &
       , (/ggrid%clonm%dim(1)%len,ggrid%clonm%dim(2)%len/) &
       ,VTYPE_DOUBLE)
  DO i=1, ie_g
     DO j = 1 , je_g
        n = (j-1) * ggrid%clonm%dim(1)%len + i
        ggrid%clonm%dat%vd(n) = gila_g(i*2,j*2)*180._dp/pi
     END DO
  END DO

  ! ... attributes
  CALL ADD_NCATT(ggrid%clonm, 'long_name', vs='longitude')
  CALL ADD_NCATT(ggrid%clonm, 'units', vs='degrees_east')
  ggrid%ranges(1,1) = RGEMPTY  ! equidistant; no correction of range required
  ggrid%ranges(1,2) = RGEMPTY

  ! LATITUDE (MID) ...
  ggrid%clatm%name  = 'lat'
  ggrid%clatm%id    = NULL_VARID
  ggrid%clatm%xtype = NF90_DOUBLE
  ! ... dimensions
  ggrid%clatm%ndims = 2
  ALLOCATE(ggrid%clatm%dim(ggrid%clatm%ndims), STAT=status)
  CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
  ggrid%clatm%dim(1)%name  = 'lon'
  ggrid%clatm%dim(1)%id    = NULL_DIMID
  ggrid%clatm%dim(1)%len   = ie_g
  ggrid%clatm%dim(1)%fuid  = .false.
  ggrid%clatm%dim(1)%varid = NULL_VARID
  ggrid%clatm%dim(2)%name  = 'lat'
  ggrid%clatm%dim(2)%id    = NULL_DIMID
  ggrid%clatm%dim(2)%len   = je_g
  ggrid%clatm%dim(2)%fuid  = .false.
  ggrid%clatm%dim(2)%varid = NULL_VARID
  ! ... data
  CALL INIT_NARRAY(ggrid%clatm%dat, ggrid%clatm%ndims &
       , (/ggrid%clatm%dim(1)%len, ggrid%clatm%dim(2)%len/) &
       ,VTYPE_DOUBLE)
  DO i = 1, ie
     DO j=1,je
        n = (j-1) * ggrid%clatm%dim(1)%len + i
        ggrid%clatm%dat%vd(n) = giph_g(i*2,j*2)*180._dp/pi
     END DO
  END DO

  ! ... attributes
  CALL ADD_NCATT(ggrid%clatm, 'long_name', vs='mid latitude')
  CALL ADD_NCATT(ggrid%clatm, 'units',     vs='degrees_north')

  ggrid%ranges(2,1) = RGEMPTY !-90.0_dp
  ggrid%ranges(2,2) = RGEMPTY

  ! ----------------------------------------------------------
  ! define interfaces:
  ! Curvilinear LONGITUDE (INTERFACES) ...
  ggrid%cloni%name  = 'lon_I'
  ggrid%cloni%id    = NULL_VARID
  ggrid%cloni%xtype = NF90_DOUBLE
  ! ... dimensions
  ggrid%cloni%ndims = 2
  ALLOCATE(ggrid%cloni%dim(ggrid%cloni%ndims), STAT=status)
  CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
  ggrid%cloni%dim(1)%name  = 'lon_I'
  ggrid%cloni%dim(1)%id    = NULL_DIMID
  ggrid%cloni%dim(1)%len   = ie_g+1
  ggrid%cloni%dim(1)%fuid  = .false.
  ggrid%cloni%dim(1)%varid = NULL_VARID
  ggrid%cloni%dim(2)%name  = 'lat_I'
  ggrid%cloni%dim(2)%id    = NULL_DIMID
  ggrid%cloni%dim(2)%len   = je_g+1
  ggrid%cloni%dim(2)%fuid  = .false.
  ggrid%cloni%dim(2)%varid = NULL_VARID
  ! ... data
  CALL INIT_NARRAY(ggrid%cloni%dat, ggrid%cloni%ndims &
       , (/ggrid%cloni%dim(1)%len,ggrid%cloni%dim(2)%len/) &
       ,VTYPE_DOUBLE)
  DO i=1, ie_g
     DO j = 1 , je_g
        n = (j-1) * ggrid%cloni%dim(1)%len + i
         ggrid%cloni%dat%vd(n) = gila_g(i*2-1,j*2)
         gloni(i,j) = ggrid%cloni%dat%vd(n)
     END DO
  END DO
  ! for ie_g
  i = 1  ! Value of ie_g+1 is equal to value at place 1
  DO j = 1 , je_g
     n = (j-1) * ggrid%cloni%dim(1)%len + ie_g+1
     ggrid%cloni%dat%vd(n) = gila_g(i*2-1,j*2)
     gloni(ie_g+1,j) = ggrid%cloni%dat%vd(n)
  END DO
  ! for je_g
  ! interface at je_g+1 does not exist
  DO i = 1 , ie_g
     n = (je_g+1-1) * ggrid%cloni%dim(1)%len + i
     ggrid%cloni%dat%vd(n) = gila_g(i*2-1,je_g*2)
     gloni(i,je_g+1) = ggrid%cloni%dat%vd(n)
  END DO
  ! ie_g+1, je_g+1
  n = (je_g+1-1) * ggrid%cloni%dim(1)%len + ie_g+1
  ggrid%cloni%dat%vd(n) = gila_g(1*2-1,je_g*2)
  gloni(ie_g+1,je_g+1) = ggrid%cloni%dat%vd(n)

  ! ... attributes
  CALL ADD_NCATT(ggrid%cloni, 'long_name', vs='interface longitude')
  CALL ADD_NCATT(ggrid%cloni, 'units',     vs='degrees_east')

  ggrid%ranges(1,1) = RGEMPTY  ! equidistant; no correction of range required
  ggrid%ranges(1,2) = RGEMPTY

  ! LATITUDE (MID) ...
  ggrid%clati%name  = 'lat_I'
  ggrid%clati%id    = NULL_VARID
  ggrid%clati%xtype = NF90_DOUBLE
  ! ... dimensions
  ggrid%clati%ndims = 2
  ALLOCATE(ggrid%clati%dim(ggrid%clati%ndims), STAT=status)
  CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
  ggrid%clati%dim(1)%name  = 'lon_I'
  ggrid%clati%dim(1)%id    = NULL_DIMID
  ggrid%clati%dim(1)%len   = ie_g+1
  ggrid%clati%dim(1)%fuid  = .false.
  ggrid%clati%dim(1)%varid = NULL_VARID
  ggrid%clati%dim(2)%name  = 'lat_I'
  ggrid%clati%dim(2)%id    = NULL_DIMID
  ggrid%clati%dim(2)%len   = je_g+1
  ggrid%clati%dim(2)%fuid  = .false.
  ggrid%clati%dim(2)%varid = NULL_VARID
  ! ... data
  CALL INIT_NARRAY(ggrid%clati%dat, ggrid%clati%ndims &
       , (/ggrid%clati%dim(1)%len, ggrid%clati%dim(2)%len/) &
       ,VTYPE_DOUBLE)
  DO i=1, ie_g
     DO j = 1 , je_g
        n = (j-1) * ggrid%clati%dim(1)%len + i
        ggrid%clati%dat%vd(n) = giph_g(i*2,j*2-1)
        glati(i,j) = ggrid%clati%dat%vd(n)
     END DO
  END DO
  ! for ie_g
  i = 1  ! Value of ie_g+1 is equal to value at place 1
  DO j = 1 , je_g
     n = (j-1) * ggrid%cloni%dim(1)%len + ie_g+1
     ggrid%clati%dat%vd(n) = gila_g(i*2,j*2-1)
     glati(ie_g+1,j) = ggrid%clati%dat%vd(n)
  END DO
  ! for je_g
  ! interface at je_g+1 does not exist
  DO i = 1 , ie_g
     n = (je_g+1-1) * ggrid%cloni%dim(1)%len + i
     ggrid%clati%dat%vd(n) = gila_g(i*2,je_g*2)
     glati(i,je_g+1) = ggrid%clati%dat%vd(n)
  END DO
  ! ie_g+1, je_g+1
  n = (je_g+1-1) * ggrid%cloni%dim(1)%len + ie_g+1
  ggrid%clati%dat%vd(n) = gila_g(1*2,je_g*2)
  glati(ie_g+1,je_g+1)  = ggrid%clati%dat%vd(n)

  ! ... attributes
  CALL ADD_NCATT(ggrid%clati, 'long_name', vs='interface latitude')
  CALL ADD_NCATT(ggrid%clati, 'units', vs='degrees_north')

  ggrid%ranges(2,1) = RGEMPTY !-90.0_sp
  ggrid%ranges(2,2) = RGEMPTY

  ! TIME (MID) ...
  ggrid%timem%name  = 'time'
  ggrid%timem%id    = NULL_VARID
  ggrid%timem%xtype = NF90_DOUBLE
  ! ... dimensions
  ggrid%timem%ndims = 1
  ALLOCATE(ggrid%timem%dim(ggrid%timem%ndims), STAT=status)
  CALL ERRMSG('INTERFACE_GEOHYBGRID',status,6)
  ggrid%timem%dim(1)%name  = 'time'
  ggrid%timem%dim(1)%id    = NULL_DIMID
  ggrid%timem%dim(1)%len   = 1
  ggrid%timem%dim(1)%fuid  = .true.
  ggrid%timem%dim(1)%varid = NULL_VARID
  ! ... data
  CALL INIT_NARRAY(ggrid%timem%dat, ggrid%timem%ndims, (/ggrid%timem%dim(1)%len/) &
                   ,VTYPE_DOUBLE)
  ! TIME: SECONDS SINCE MODEL START
  CALL time_span_d(dts   &
       , YEAR_START, MONTH_START, DAY_START &
       , HOUR_START, MINUTE_START, SECOND_START  &
       , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND)
  ggrid%timem%dat%vd(1) = dts * 86400.0_dp  ! days -> seconds

  ! ... attributes
  CALL ADD_NCATT(ggrid%timem, 'long_name'                              &
               ,vs='time in seconds since model start')
  WRITE(tunit, &
       '("sec since ",I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
       YEAR_START, MONTH_START, DAY_START &
       , HOUR_START, MINUTE_START, SECOND_START

  CALL ADD_NCATT(ggrid%timem, 'units', vs=TRIM(tunit))

  ! CALCULATE INTs from MIDs
  CALL COMPLETE_GEOHYBGRID(ggrid)

  CALL new_geohybgrid(status, GGRID_ID, ggrid)
  IF (status /= 0 .AND. status /= 01) &
       CALL error_bi(grid_error(status), 'mpiom grid definition')

  ! --------------------------------
  ! DEFINE LOCAL MPIOM
  ! --------------------------------

  ! INIT
  CALL INIT_GEOHYBGRID(lgrid)

  lgrid%name = 'MPIOMl'

  lgrid%file = 'MPIOMl-_GEO-HYBRID-GRID'       ! Filename
  lgrid%t    = 0                              ! time step

  lgrid%clonc = .FALSE. ! Longitude Axis is circular
  ! Curvilinear LONGITUDE (MID) ...
  lgrid%clonm%name  = 'lon'
  lgrid%clonm%id    = NULL_VARID
  lgrid%clonm%xtype = NF90_DOUBLE
  ! ... dimensions
  lgrid%clonm%ndims = 2
  ALLOCATE(lgrid%clonm%dim(lgrid%clonm%ndims), STAT=status)
  CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
  lgrid%clonm%dim(1)%name  = 'lon'
  lgrid%clonm%dim(1)%id    = NULL_DIMID
  lgrid%clonm%dim(1)%len   = ie
  lgrid%clonm%dim(1)%fuid  = .false.
  lgrid%clonm%dim(1)%varid = NULL_VARID
  lgrid%clonm%dim(2)%name  = 'lat'
  lgrid%clonm%dim(2)%id    = NULL_DIMID
  lgrid%clonm%dim(2)%len   = je
  lgrid%clonm%dim(2)%fuid  = .false.
  lgrid%clonm%dim(2)%varid = NULL_VARID
  ! ... data
  CALL INIT_NARRAY(lgrid%clonm%dat, lgrid%clonm%ndims &
       , (/lgrid%clonm%dim(1)%len,lgrid%clonm%dim(2)%len/) &
       ,VTYPE_DOUBLE)
  DO i=1, ie
     DO j = 1 , je
        n = (j-1) * lgrid%clonm%dim(1)%len + i
        idx= i + p_ioff
        jdx= j + p_joff
        lgrid%clonm%dat%vd(n) = gila_g(idx*2,jdx*2)*180._dp/pi
     END DO
  END DO

  ! ... attributes
  CALL ADD_NCATT(lgrid%clonm, 'long_name', vs='longitude')
  CALL ADD_NCATT(lgrid%clonm, 'units', vs='degrees_east')
  lgrid%ranges(1,1) = RGEMPTY  ! equidistant; no correction of range required
  lgrid%ranges(1,2) = RGEMPTY

  ! LATITUDE (MID) ...
  lgrid%clatm%name  = 'lat'
  lgrid%clatm%id    = NULL_VARID
  lgrid%clatm%xtype = NF90_DOUBLE
  ! ... dimensions
  lgrid%clatm%ndims = 2
  ALLOCATE(lgrid%clatm%dim(lgrid%clatm%ndims), STAT=status)
  CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
  lgrid%clatm%dim(1)%name  = 'lon'
  lgrid%clatm%dim(1)%id    = NULL_DIMID
  lgrid%clatm%dim(1)%len   = ie
  lgrid%clatm%dim(1)%fuid  = .false.
  lgrid%clatm%dim(1)%varid = NULL_VARID
  lgrid%clatm%dim(2)%name  = 'lat'
  lgrid%clatm%dim(2)%id    = NULL_DIMID
  lgrid%clatm%dim(2)%len   = je
  lgrid%clatm%dim(2)%fuid  = .false.
  lgrid%clatm%dim(2)%varid = NULL_VARID
  ! ... data
  CALL INIT_NARRAY(lgrid%clatm%dat, lgrid%clatm%ndims &
       , (/lgrid%clatm%dim(1)%len, lgrid%clatm%dim(2)%len/) &
       ,VTYPE_DOUBLE)
  DO i = 1, ie
     DO j=1,je
        n = (j-1) * lgrid%clatm%dim(1)%len + i
        idx= i + p_ioff
        jdx= j + p_joff
        lgrid%clatm%dat%vd(n) = giph_g(idx*2,jdx*2)*180._dp/pi
     END DO
  END DO

  ! ... attributes
  CALL ADD_NCATT(lgrid%clatm, 'long_name', vs='mid latitude')
  CALL ADD_NCATT(lgrid%clatm, 'units',     vs='degrees_north')

  lgrid%ranges(2,1) = RGEMPTY !-90.0_dp
  lgrid%ranges(2,2) = RGEMPTY

  ! ----------------------------------------------------------
  ! define interfaces:
  ! Curvilinear LONGITUDE (INTERFACES) ...
  lgrid%cloni%name  = 'lon_I'
  lgrid%cloni%id    = NULL_VARID
  lgrid%cloni%xtype = NF90_DOUBLE
  ! ... dimensions
  lgrid%cloni%ndims = 2
  ALLOCATE(lgrid%cloni%dim(lgrid%cloni%ndims), STAT=status)
  CALL ERRMSG('INTERFACE_GEOHYBGRID',status,1)
  lgrid%cloni%dim(1)%name  = 'lon_I'
  lgrid%cloni%dim(1)%id    = NULL_DIMID
! um_ak_20130705+
!  lgrid%cloni%dim(1)%len   = ie_g+1
  lgrid%cloni%dim(1)%len   = ie+1
! um_ak_20130705-
  lgrid%cloni%dim(1)%fuid  = .false.
  lgrid%cloni%dim(1)%varid = NULL_VARID
  lgrid%cloni%dim(2)%name  = 'lat_I'
  lgrid%cloni%dim(2)%id    = NULL_DIMID
! um_ak_20130705+
  !lgrid%cloni%dim(2)%len   = je_g+1
  lgrid%cloni%dim(2)%len   = je+1
! um_ak_20130705-
  lgrid%cloni%dim(2)%fuid  = .false.
  lgrid%cloni%dim(2)%varid = NULL_VARID
  ! ... data
  CALL INIT_NARRAY(lgrid%cloni%dat, lgrid%cloni%ndims &
       , (/lgrid%cloni%dim(1)%len,lgrid%cloni%dim(2)%len/) &
       ,VTYPE_DOUBLE)

! um_ak_20130705+
!  DO i=1, ie_g+1
!     DO j = 1 , je_g+1
  DO i=1, ie+1
     DO j = 1 , je+1
! um_ak_20130705-
        n = (j-1) * lgrid%cloni%dim(1)%len + i
        idx= i + p_ioff
        jdx= j + p_joff
        lgrid%cloni%dat%vd(n) = gloni(idx,jdx)
     END DO
  END DO

  ! ... attributes
  CALL ADD_NCATT(lgrid%cloni, 'long_name', vs='interface longitude')
  CALL ADD_NCATT(lgrid%cloni, 'units',     vs='degrees_east')

  lgrid%ranges(1,1) = RGEMPTY  ! equidistant; no correction of range required
  lgrid%ranges(1,2) = RGEMPTY

  ! LATITUDE (MID) ...
  lgrid%clati%name  = 'lat_I'
  lgrid%clati%id    = NULL_VARID
  lgrid%clati%xtype = NF90_DOUBLE
  ! ... dimensions
  lgrid%clati%ndims = 2
  ALLOCATE(lgrid%clati%dim(lgrid%clati%ndims), STAT=status)
  CALL ERRMSG('INTERFACE_GEOHYBGRID',status,2)
  lgrid%clati%dim(1)%name  = 'lon_I'
  lgrid%clati%dim(1)%id    = NULL_DIMID
  lgrid%clati%dim(1)%len   = ie+1
  lgrid%clati%dim(1)%fuid  = .false.
  lgrid%clati%dim(1)%varid = NULL_VARID
  lgrid%clati%dim(2)%name  = 'lat_I'
  lgrid%clati%dim(2)%id    = NULL_DIMID
  lgrid%clati%dim(2)%len   = je+1
  lgrid%clati%dim(2)%fuid  = .false.
  lgrid%clati%dim(2)%varid = NULL_VARID
  ! ... data
  CALL INIT_NARRAY(lgrid%clati%dat, lgrid%clati%ndims &
       , (/lgrid%clati%dim(1)%len, lgrid%clati%dim(2)%len/) &
       ,VTYPE_DOUBLE)
! um_ak_20130705+
!  DO i=1, ie_g+1
!     DO j = 1 , je_g+1
  DO i=1, ie+1
     DO j = 1 , je+1
! um_ak_20130705-
        n = (j-1) * lgrid%clati%dim(1)%len + i
        idx= i + p_ioff
        jdx= j + p_joff
        lgrid%clati%dat%vd(n) = glati(idx,jdx)
     END DO
  END DO

  ! ... attributes
  CALL ADD_NCATT(lgrid%clati, 'long_name', vs='interface latitude')
  CALL ADD_NCATT(lgrid%clati, 'units', vs='degrees_north')

  lgrid%ranges(2,1) = RGEMPTY !-90.0_sp
  lgrid%ranges(2,2) = RGEMPTY

  ! TIME (MID) ...
  lgrid%timem%name  = 'time'
  lgrid%timem%id    = NULL_VARID
  lgrid%timem%xtype = NF90_DOUBLE
  ! ... dimensions
  lgrid%timem%ndims = 1
  ALLOCATE(lgrid%timem%dim(lgrid%timem%ndims), STAT=status)
  CALL ERRMSG('INTERFACE_GEOHYBGRID',status,6)
  lgrid%timem%dim(1)%name  = 'time'
  lgrid%timem%dim(1)%id    = NULL_DIMID
  lgrid%timem%dim(1)%len   = 1
  lgrid%timem%dim(1)%fuid  = .true.
  lgrid%timem%dim(1)%varid = NULL_VARID
  ! ... data
  CALL INIT_NARRAY(lgrid%timem%dat, lgrid%timem%ndims, (/lgrid%timem%dim(1)%len/) &
                   ,VTYPE_DOUBLE)
  ! TIME: SECONDS SINCE MODEL START
  CALL time_span_d(dts   &
       , YEAR_START, MONTH_START, DAY_START &
       , HOUR_START, MINUTE_START, SECOND_START  &
       , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND)
  lgrid%timem%dat%vd(1) = dts * 86400.0_dp  ! days -> seconds

  ! ... attributes
  CALL ADD_NCATT(lgrid%timem, 'long_name'                              &
               ,vs='time in seconds since model start')
  WRITE(tunit, &
       '("sec since ",I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
       YEAR_START, MONTH_START, DAY_START &
       , HOUR_START, MINUTE_START, SECOND_START

  CALL ADD_NCATT(lgrid%timem, 'units', vs=TRIM(tunit))

  ! CALCULATE INTs from MIDs
  CALL COMPLETE_GEOHYBGRID(lgrid)

  CALL new_geohybgrid(status, LGRID_ID, lgrid)
  IF (status /= 0 .AND. status /= 01) &
       CALL error_bi(grid_error(status), 'mpiom grid definition')
!#endif

  END SUBROUTINE define_mpiom_grid
  ! um_ak_20130620- ! mz_ap_20130620-

  SUBROUTINE mpiom_read_nml_init(status, iou)

    ! read namelist for 'coupling' to ECHAM5
    !
    ! Author: Andrea Pozzer, MPICH, Aug 2003

    ! MESSy
    USE messy_main_tools,         ONLY: read_nml_open, read_nml_check, read_nml_close
!!$    USE messy_main_mpi_bi,        ONLY: finish

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /INIT/ mpiom_init_file_name


    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='mpiom_read_nml_init'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status
    INTEGER                     :: jt

    status = 1

    CALL read_nml_open(lex, substr, iou, 'INIT', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=INIT, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'INIT', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

   WRITE(*,*) '.........................................................'
   WRITE(*,*) '       MPIOM INITILIZATION FROM NETCDF FILE              '
   WRITE(*,*) ' FILE=', mpiom_init_file_name
   WRITE(*,*) '.........................................................'

  END SUBROUTINE mpiom_read_nml_init

! ===========================================================================

END MODULE messy_mpiom_e5
