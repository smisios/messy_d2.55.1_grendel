PROGRAM MPIOM
!
!    MPIOM
!
!     Version:
!     ---------
!     $URL: http://svn.zmaw.de/svn/mpiom/branches/mpiom-latest/src/mpiom.f90 $
!     $Rev: 1999 $
!
!    VERSION OF THE HOPE OGCM ON C-GRID
!
!    DEVELOPED BY MAIER-REIMER TILL 1997
!
!    VERSION OF UWE MIKOLAJEWICZ,JOHANN JUNGCLAUS AND HELMUTH HAAK, 12/99
!    INCLUDES CONFORMAL MAPPING
!
!    MODIFIED :
!    ---------
!   HOPS65 : CREATED JULI 25, 2000  H. HAAK, J. JUNGCLAUS
!   HOPS67 : CREATED  NOV 09, 2000  H. HAAK, J. JUNGCLAUS, U.MIKOLAJEWICZ
!   HOPS68 : CREATED JUNE 20, 2001  H. HAAK, J. JUNGCLAUS, U.MIKOLAJEWICZ
!                        Nov. 2001  O. Boeringer  - TRANSFER TO FORTRAN90
!                                   S. Legutke    - interface to HAMOCC5
!   HOPS69 : CREATED   JAN 5, 2002  H. HAAK, J. JUNGCLAUS, U.MIKOLAJEWICZ
!   MPIOM  :           JAN 14 2003  S. Legutke    - created mo_couple.F90
!   MPIOM  :           JUN 03 2003  J. Jungclaus  - update GIRIV, RELTEM, ZOCORR
!   MPIOM  :           AUG 15 2007  S. Lorenz     - omp-parallel tracer loops
!*************************************************************************
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
!                                                                      *
!UWE  AKTUALISIEREN!!!!!!!!!!!!!!!!
!
!     DATA SETS AND FORTRAN CHANNEL NUMBER CODING                      *
!                                                                      *
!     UNIT NO.   NAME     DESCRIPTION                      LOCATION    *
!
!HH    IO_IN_INIT INITEM  TEMPERATURE LEVITUS          MPIOM.F90,SBR LEVIRE
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!
  USE MO_PARAM1
  USE mo_constants, ONLY: api
  USE MO_MPI
  USE MO_PARALLEL
#ifndef MESSY
  USE mo_parallel_diags, ONLY: timing_diagnostics
#endif
  USE mo_boundsexch, ONLY : bounds_exch, print_stats &
                           ,alloc_mem_bounds_exch_halo

  USE mo_restart, ONLY: setup_restart_defaults, setup_restart, &
#ifdef OLD_IO
       aufr, aufw, &
#endif
       lm_start, ly_start, ly_end
  USE MO_COMMO1

#ifndef NO_NEW_IO
  USE mo_iolist, ONLY : iolist_init, iolist_fini, iolist_poll_write_avglist, &
                        iolist_write_initial, iolist_write_final, readRestartCDI, &
                        iolist_read_config, iolist_create_file_list, &
                        iolist_close_file_list, ocean_iolist, &
                        iolist_accumulate_staggered, iolist_accumulate
  USE mo_varlist, ONLY: ocean_varlist, build_ocean_varlist, &
       varlist_post_process_ocean_data
#endif
  USE mo_contro,  ONLY: contro
  USE mo_levitus
  USE MO_COMMOAU1
  USE MO_COMMOAU2
  USE MO_COMMOAU3
  USE MO_DIAGNOSIS
  USE MO_DIFFUSION, ONLY : alloc_mem_octdiff, octdiff_base, octdiff_trf, aulapts, ah00
  USE mo_ocean_vertical_mixing, ONLY: av0, dv0, setup_ocean_vertical_mixing
  USE mo_grid, ONLY: boden, coriol, setup_grid, wrte_gridinfo, lwith_one_layer_shelfs, cell_thickness
  USE MO_OCTHER
  USE MO_OCICE, ONLY: lsaoclose, ocice
  USE mo_swr_absorption, ONLY : alloc_mem_swr_absorption, &
       jerlov_swr_absorption, old_swr_absorption &
       ,jerlov_atten, jerlov_bluefrac,lfb_bgc_oce
  USE MO_TRO
  USE MO_EDDYDIAG
  USE mo_runoff, ONLY : river_runoff_omip_ini, river_runoff_stations_ini, glac_calv_ini, &
         luse_river_runoff_stations,numriv,luse_glac_calv,numglac

  USE mo_ocvisc, ONLY : bofric, rayfric, ocvisc

#ifdef _PROFILE
  USE mo_profile,      ONLY: trace_init, trace_start, trace_stop, &
                             trace_finalize
#endif
#ifndef MESSY
  USE mo_fpe, ONLY: enable_fpe_tracing
#endif

#ifdef CORE
  USE MO_NCAR_OCEAN_FLUXES
#else
  USE MO_OMIP
#endif

  USE mo_forcing, ONLY: read_namelist_forcctl, initialize_surface_forcing, &
       update_surface_forcing, finalize_surface_forcing, forcing_frequency
#ifndef __coupled
  USE mo_forcing, ONLY: spool_forcing
#endif

#ifdef NUDGE_ECCO
  USE MO_NUDGE_TS
#endif

  USE mo_tidal, ONLY : alloc_mem_tidal,foreph_ini,foreph,tipouv &
       ,init_tide_timeseries,p_tide_timeseries

#ifdef FLUXCORRECT
  USE MO_COMMO_FLUXCORR
#endif /*FLUXCORRECT*/

  USE MO_ADPO
  USE MO_COMMOBBL

  USE MO_MEAN

  USE MO_ELICOM
  USE MO_PARA2

  USE MO_UNITS

  USE mo_basin_masks, only: load_basin_masks
#ifndef MESSY
  USE MO_ISO_C_KINDS
#else
  USE mo_kind, ONLY : c_int64_t => i8, c_float => sp
#endif

  USE mo_model_time, ONLY: time_desc, monlen, inc_day, &
                           format_model_time, write_model_date, &
                           broadcast_time_desc, model_time_init, &
#ifdef PBGC
                           seconds_between_two_times, &
#endif
                           operator(+), operator(-)

#if defined (__coupled)
  USE mo_fluxes1, ONLY : &
#ifdef PBGC
       aoflwsvo, aoflshwo, &
#endif
       alloc_mem_fluxes1

  USE mo_couple, ONLY: couple_prep, couple_init, &
                       couple_put_o2a, couple_calendar , &
                       couple_end
#endif

#ifdef PBGC
  USE mo_bgc_diagnostic, ONLY : rate_tracer, store_tracer
#ifndef NO_NEW_IO
  USE mo_bgc_varlist, ONLY : bgc_varlist
  USE mo_bgc_iolist, ONLY : bgc_iolist
#endif

  USE mo_carbch
  USE mo_control_bgc
  USE mo_biomod
  USE mo_sedmnt
  USE mo_param1_bgc
#endif /*PBGC*/

  IMPLICIT NONE
#ifdef PBGC
!      REAL(wp), POINTER :: layer1_bgc(:,:)
!      REAL(wp), POINTER :: layer1_new(:,:)
  REAL(wp), ALLOCATABLE :: bgcddpo(:,:,:)
  REAL(wp), ALLOCATABLE :: bgcdpio(:,:,:)
  INTEGER ndtrun
#endif /*PBGC*/


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     DECLARATIONS
#ifdef __coupled
  INTEGER NACTYEAR
#endif
  REAL(wp) ABACK,CAULAPTS,CAULAPUV,CWT,CWA,DBACK,GFDL_DIFF
  REAL(wp) RRELSAL,RRELTEM

  INTEGER I,J,K,M,N,NANF,ENDDAY
  INTEGER LDTRUN, LDTDAY, LMON1, LMON2
  !
  !> time step variables
  INTEGER LDTYEAR, LDTMONTH, ldt
  INTEGER :: NYEARS, NMONTS, NDAYS

  INTEGER :: L, ll, lmont, ntracerloop
  INTEGER :: ihalo_sor, imod
  REAL(dp) :: ttts, tttr, ttt
  ! write grid information back to file
  LOGICAL :: lgridinfo = .FALSE.
  ! check floating-point operations for exceptions
  LOGICAL :: fp_tracing_enabled = .FALSE.
  ! apply sea level correction?
  LOGICAL :: lzo_correct = .TRUE.

  ! print timing diagnostics (0 = none, 1 = aggregate, 2 = every task,
  !                           3 = individual routines)
  INTEGER :: time_verbosity
#ifdef SORTEST
  REAL(wp), ALLOCATABLE :: z1(:,:),z2(:,:)
#endif

#ifdef CLOCK
  REAL(wp) :: zwtime,zutime,zstime, zrtime

!  External functions
  REAL(wp), EXTERNAL :: util_walltime
  INTEGER, EXTERNAL :: util_cputime
#endif

  LOGICAL ::lundelayed_momentum_advection = .FALSE.

  TYPE(time_desc) :: model_start_time, model_time, &
                     run_start_time, run_end_time, next_run_start_time

  INTEGER, PARAMETER :: tag_subr_surf_forcing=1, tag_subr_bgc=2, &
       tag_subr_octher=3, tag_subr_ocmodmom=4, tag_subr_occlit=5, &
       tag_subr_solver=6, tag_subr_ocvtot=7, tag_subr_ocuvad=8, &
       tag_subr_ocadfs=9, tag_region_adv=10, tag_subr_ocschep=11, &
       tag_subr_ocvisc=12, &
       tag_timestep=13, num_timing_tags=tag_timestep
  CHARACTER(len=26), PARAMETER :: region_tags(num_timing_tags) = (/ &
       'initialize surface forcing', &
       'BGC                       ', &
       'octher                    ', &
       'ocmodmom                  ', &
       'occlit                    ', &
       'solver                    ', &
       'ocvtot                    ', &
       'ocuad+ocvad               ', &
       'ocadfs                    ', &
       'advection/diffusion/gm    ', &
       'ocschep                   ', &
       'ocvisc                    ', &
       'complete timestep         ' /)
  REAL(dp) :: time_diffs_subr(num_timing_tags)
#ifndef MESSY
  TYPE(min_mean_max_dp) :: timestep_time_agg(num_timing_tags)
#endif


  CALL setup_defaults
  CALL setup_initial_rte

  CALL read_user_setup

  CALL reset_inconsistent_user_setup

  CALL broadcast_user_setup

  IF (cwa .LT. 0._wp) cwa = cwt

  dti = 1._wp / dt
  ! FIXME: why not put 86400 as solar_day in mo_planetary constants?
  ndtday = NINT(86400._wp/dt)


  CALL logwrite_user_setup


#ifndef MESSY
  IF (fp_tracing_enabled) CALL enable_fpe_tracing
#endif


  CALL setup_time_invariants


  PER_YEAR: DO LYEAR=LYEAR1,LYEAR2
    CALL per_year_setup
    MONTH_OF_YEAR: DO LMONT=LMON1,LMON2
      CALL per_month_setup
      DAY_OF_MONTH: DO LDAY=LDAY1,ENDDAY
        LDAYS=LDAY
        TIME_STEP_OF_DAY: DO LDTDAY = 1,NDTDAY
          CALL per_timestep_computation
        END DO TIME_STEP_OF_DAY
        ! END OF ONE DAY

        ! correct ZO to a global mean of zero once per day
        IF (lzo_correct) CALL correct_zo
      END DO DAY_OF_MONTH
!#ifdef PBGC
!      CALL AVRG_BGCMEAN(ie,je,ke)
!#endif /*PBGC*/
      IF ( LCALCDIFI) THEN
        CALL CALC_DIFI((LYEARS*10000)+(LMONTS*100)+LDAYS)
      END IF
      ! END OF ONE MONTH
#ifdef FLUXCORRECT
      CALL FLUX_CORRECT
#endif /*FLUXCORRECT*/
    END DO MONTH_OF_YEAR
    ! END OF ONE YEAR
  END DO PER_YEAR

  CALL finalize_surface_forcing

  !----------------------------------------------------------------------
  ! END OF TIME STEPPING

  CALL write_restart_file

#ifdef KONVDIAG
  CALL WRTE_KONVDIAG(LYEAR2,LMONT2,MONLEN(LMONT2,lyear2))
#endif /*KONVDIAG*/

  CALL WRTE_GRIDINFO(LYEAR2,LMONT2,MONLEN(LMONT2,lyear2))

#ifdef FLUXCORRECT
  CALL FLUXC_WRTE
#endif /*FLUXCORRECT*/

#ifdef PBGC
!----------------------------------------------------------------------
! Finish cleanly with marine bgc
!
  CALL END_BGC(ie,je,ke,ddpo,dlxp,dlyp,gila,giph,tiestu     &
       ,lyears,lmonts,ldays,ldt)
#endif /*PBGC*/

#if defined (__coupled)

! Finish cleanly

  CALL couple_end
#endif

#ifndef NO_NEW_IO
  call iolist_close_file_list(ocean_iolist)
#ifdef PBGC
  call iolist_close_file_list(bgc_iolist)
#endif
  CALL iolist_fini
#endif

!     Branch target for finishing cleanly
!99999 CONTINUE
  CALL print_stats

  WRITE(nerr,*) 'NORMAL END OF MPIOM'

#ifdef CLOCK
  IF (util_cputime(zutime, zstime) == -1) THEN
    WRITE(nerr,*)'Cannot determine used CPU time'
  ELSE
    zwtime = util_walltime()
    zrtime = (zutime+zstime)/zwtime

    IF(p_pe == p_io) THEN
      WRITE (nerr,'(a,f10.2,a)') ' Wallclock        : ', zwtime, ' s'
      WRITE (nerr,'(a,f10.2,a)') ' CPU-time (user)  : ', zutime, ' s'
      WRITE (nerr,'(a,f10.2,a)') ' CPU-time (system): ', zstime, ' s'
      WRITE (nerr,'(a,f10.2,a)') ' Ratio            : ', 100._wp * zrtime, ' %'
    END IF
  END IF
#endif

  IF (p_pe == p_io) CALL write_model_date(run_start_time, run_end_time, &
                    next_run_start_time)

#ifdef _PROFILE
  CALL trace_stop ('all', 0)
  CALL trace_finalize (p_pe)
#endif


#ifdef bounds_exch_put
  CALL close_put_window
#endif
  CALL p_stop



CONTAINS


  !> initialize variables that can be changed by the user
  !! and customize the run-time behaviour
  SUBROUTINE setup_defaults
#ifdef __IFC /* Intel Fortran Compiler */
    INTEGER ieee_handler
    EXTERNAL handler
    !
    n = ieee_handler('set','division',handler)
    n = ieee_handler('set','overflow',handler)
    n = ieee_handler('set','invalid',handler)
#endif


#ifdef CLOCK
    ! Initialize wallclock timer
    zwtime = util_walltime()
#endif


    ! initialize grid (default => false for bipolar grids)
    lbounds_exch_tp = .FALSE.
    ! set dimensions to illegal values so we can
    ! identify missing user setup values
    ie_g = -1
    je_g = -1
    ke = -1

    ! initialize MPI
    CALL p_start

#ifdef _PROFILE
    CALL trace_init ('mpiom', p_pe)
    CALL trace_start ('all', 0)
#endif

    ! do only print aggregate time information by default
    time_verbosity = 1

    !     SPECIFY LOGICAL I/O UNITS

    CALL SETUNITS
    !  OPEN STD OUT FILE
    !  This is needed for the MPI version since we do not want
    !  the output of all processors intermixed (and basically
    !  we want only the ouput of processor 0)

    IF (io_stdout > 0) THEN
      CALL OPEN_STDOUT(io_stdout,'oceout')
    END IF

    ! start of integration
    model_start_time = time_desc(0, 1, 1, 0, 0, 0)

    ! length of integration
    nyears = 0
    nmonts = 1
    nanf = 0
    ndays = 0

    !HH    DEFAULT TIME STEPS
    dt = 1800._wp

    nfixYearLen = -1    ! use real years : other options are fixed 365 days or fixed 360 days

    !
    !----------------------------------------------------------------------
    ! DEFAULT PARAMETER SETTINGS - CAN BE OVERWRITEN BY THE NAMELIST
    istart = istart_run_continue
    !
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

    !UWE   CONSTANTS FOR BIHARMONIC FRICTION
    caulapuv = 0.0045_wp

    !UWE   CONSTANTS FOR HARMONIC FRICTION
    aus = 3.E-6_wp

    !BOFRIC  FRICTION COEFFICENT AT BOTTOM [M**2/S]
    bofric = 1.E-3_wp
    !COEFFICIENT FOR QUADRATIC BOTTOM FRICTION []
    rayfric = 3.E-3_wp

    !HH    CONSTANT FOR DIFFUSION IN OCTHER
    dv0 = 0.5E-2_wp
    av0 = 0.5E-2_wp
    cdvocon = 20._wp
    cavocon = 0._wp

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
    IMEAN = 2

    !HH    SURFACE RELAXATION TIME SALINITY
    crelsal = 3.E-7_wp

    !HH    SURFACE RELAXATION TIME TEMERATURE
    creltem = 0._wp

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

    !  SWITCH FOR WRITING OASIS COUPLED FLUXES

    IOASISFLUX = 0


    LGMDIAG    = .FALSE.
    LFORCEDIAG = .FALSE.
    LHFLDIAG   = .FALSE.
    LCONVDIAG  = .FALSE.
    LDIFFDIAG  = .FALSE.
    LCALCDIFI  = .FALSE.

    !----------------------------------------------------------------------
    !FR   Initialize timeseries writing
    ! itsdiag=0     :  No Output
    ! itsdiag=1     :  one snapshot per day
    ! itsdiag=2     :  monthly averaged snapshots
    ! itsdiag=3     :  yearly  averaged snapshots
    ! itsdiag=4     :  output every timestep
    ! itsdiag=5     :  daily   average
    ! itsdiag=7     :  monthly average of daily means
    ! itsdiag=7     :  yearly  average of daily means
    ! ltstranspose=T:  one value for each code in diagnosis (extra-format)
    ! ltstranspose=F:  all codes in one array  in diagnosis (pseudo-extra)

    itsdiag = 5
    ltstranspose = .TRUE.
    ltswrite = .TRUE.     !  set to .FALSE. after first call of diagnosis

    CALL setup_restart_defaults
    !  reading and writing of restart files
    IAUFR = 1
    IAUFW = 1

  END SUBROUTINE setup_defaults

  !> initialize variables that are modified during the simulation
  SUBROUTINE setup_initial_rte
    ttts = 0.0_dp
    tttr = 0.0_dp
    ttt  = 0.0_dp


  END SUBROUTINE setup_initial_rte

  !>  read ocean namelist
  SUBROUTINE read_user_setup
    INTEGER :: ierror
    REAL(wp) :: cdzw(max_ke)
    ! ---------------------------------------------------------------------
    !
    !*    *NAMELIST* *OCECTL*   - Defines namelist parameters.
    !
    !*     VARIABLE  TYPE        PURPOSE.
    !      --------  ----        --------
    !      *DT*      *REAL*      Ocean time step in sec
    !
    ! ---------------------------------------------------------------------
    !
    ! declared in mo_param1/mo_commo1
    NAMELIST /OCEDIM/ ie_g, je_g, ke, lbounds_exch_tp

    ! declared in mo_parallel
    NAMELIST /NPROCS/ nprocx, nprocy

    NAMELIST /OCECTL/                                              &
              DT                                                       &
             ,CAULAPTS, CAULAPUV, CAH00, AUS, bofric, rayfric          &
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
             ,NYEARS, NMONTS, NDAYS, IMEAN, LY_END, LY_START,LM_START  &
             ,ICONTRO,IOASISFLUX                                       &
             ,LFORCEDIAG,LHFLDIAG,LCONVDIAG,LDIFFDIAG,LGMDIAG          &
             ,LGRIDINFO,LCALCDIFI,lmpitype,lnonblock                   &
             ,ITSDIAG,LTSTRANSPOSE,ihalo_sor,imod                      &
             ,fp_tracing_enabled,ibbl_transport, lsaoclose             &
             ,iter_sor, rtsorpar, iter_sor_hack,rtsorpar_hack          &
             ,lzo_correct,lundelayed_momentum_advection                &
             ,lwith_barotropic_stokes_drift, model_start_time          &
             ,time_verbosity

    NAMELIST /OCEDZW/ CDZW

    ierror = 0

    IF (p_pe == p_io) THEN
      !     Open parameter input file and ...
      OPEN(io_in_octl, file='OCECTL', status='UNKNOWN', &
           access='SEQUENTIAL', form='FORMATTED')
      ! ... read grid dimensions and grid type (regular/tripolar) ...
      READ(io_in_octl, ocedim, iostat=ierror)
      IF (ierror /= 0) THEN
        WRITE(0, *)'read OCEDIM : error=', ierror
        CALL stop_all('error in reading namelist OCEDIM - run aborted')
      END IF
      IF (ie_g < 1 .OR. je_g < 1 .OR. ke < 1) THEN
        WRITE(0, *) 'illegal dimension, must be positive integer: ', &
             'ie_g:', ie_g, 'je_g:', je_g, 'ke:', ke
        CALL stop_all('error while reading namelist - run aborted')
      END IF
      ! ... read the number of processors along x- and y-direction
      READ(io_in_octl, nprocs, iostat=ierror)
      IF (ierror /= 0) THEN
        WRITE(0, *)'read NPROCS : error=', ierror
        CALL stop_all('error in reading namelist NPROCS - run aborted')
      END IF

    END IF

    CALL setup_dimensional_data

    IF (p_pe == p_io) THEN

      READ(io_in_octl, ocectl, iostat=ierror)
      IF (ierror /= 0) THEN
        WRITE(0, *) 'read OCECTL : error=', ierror
        CALL stop_all('error in reading namelist OCECTL - run aborted')
      END IF


      READ(io_in_octl, ocedzw, iostat=ierror)
      IF (ierror /= 0) THEN
        WRITE(0, *) 'read OCEDZW : error=', ierror
        CALL stop_all('error in reading namelist OCEDZW - run aborted')
      END IF
      DZW(1:ke) = CDZW(1:ke)
    END IF

#ifndef NO_NEW_IO
    CALL iolist_read_config(io_in_octl, ierror)
#endif

    IF (p_pe == p_io .AND. ierror /= 0) THEN
      WRITE(0, *) 'read IOCTL : error=', ierror
      CALL stop_all('error in reading namelist IOCTL - run aborted')
    END IF

    CALL read_namelist_forcctl(model_start_time, io_in_octl, ierror)
    IF (p_pe == p_io .AND. ierror /= 0) THEN
      WRITE(0, *) 'read FORCCTL : error=', ierror
      CALL stop_all('error in reading namelist FORCCTL - run aborted')
    END IF

    IF (p_pe == p_io) THEN
      CLOSE(io_in_octl)
      IF ( ierror .NE. 0 ) THEN
        CALL stop_all('error in reading namelist on close - run aborted')
      ELSE
        WRITE(0,*) 'read OCECTL successfully'
      END IF
    END IF
  END SUBROUTINE read_user_setup


  SUBROUTINE setup_dimensional_data
    CALL p_bcast(ie_g,p_io)
    CALL p_bcast(je_g,p_io)
    CALL p_bcast(ke,p_io)
    CALL p_bcast(lbounds_exch_tp,p_io)
    !     Domain decomposition, determines local IE and JE
    CALL p_deco
!     WRITE(0,*) 'nach deco', p_pe, ie_g, je_g, ke, ie, je, p_ioff, p_joff

    CALL alloc_mem_bounds_exch_halo
#ifdef MESSY
    CALL init_MPI_datatypes(IE,JE,KE)
#elif
    CALL init_MPI_datatypes
#endif

    !     Set some dependent parameters and allocate arrays
    CALL set_param1
    CALL alloc_mem_commo1
    CALL alloc_mem_forcing
    CALL alloc_mem_commoau2
    CALL alloc_mem_commoau3
    CALL alloc_mem_diag
    CALL alloc_eddydiag

    CALL alloc_mem_dilcor

#ifndef NO_NEW_IO
    CALL iolist_init
#endif

#ifdef __coupled
    CALL alloc_mem_fluxes1
#endif /* __coupled */

    !#ifdef CONVDIAG
    !      CALL alloc_mem_commconv
    !#endif /*CONVDIAG*/

#ifdef PBGC
    !      ALLOCATE( layer1_bgc(ie,je) )
    !      ALLOCATE( layer1_new(ie,je) )
    ALLOCATE( bgcddpo(ie,je,ke) )
    ALLOCATE( bgcdpio(ie,je,ke) )
#endif /*PBGC*/

    CALL alloc_mem_elicom
!

#ifdef bounds_exch_put
    CALL set_put_window
#endif
#ifdef CORE
    WRITE(IO_STDOUT,*) 'allocate core'
    CALL alloc_mem_core
#endif

    CALL BELEG_ZERO



    !     DZwW(K) DECOFTES THE THICKNESS OF LAYER K IN METERS
    !     SUM OF DZW(1-KE) IS THE BOTTOM OF THE MODEL, I.E. SOLID GROUND !!!
    !     DZ(K) DECOFTES THE DISTANCE BETWEEN VECTOR POINTS BUT WILL BE
    !     COMPUTED IN SBR BODEN FOR CONSISTENCY WITH W-POINTS
    DO K=1,KE
      saf(k) = 34.8_wp
      taf(k) = 1._wp
      dzw(k) = 0._wp
    END DO

    spongezone=(/1,1,4,ie_g,je_g,ke/)     ! default is global TS restoring

  END SUBROUTINE setup_dimensional_data




  SUBROUTINE reset_inconsistent_user_setup
    IF (p_pe == p_io) THEN
      IF (imean .EQ. 0) THEN
        IF ( LGMDIAG ) WRITE(0,*)'imean = 0 => LGMDIAG disabled'
        LGMDIAG = .FALSE.
        IF ( LFORCEDIAG ) WRITE(0,*)'imean = 0 => LFORCEDIAG disabled'
        LFORCEDIAG = .FALSE.
        IF ( LHFLDIAG ) WRITE(0,*)'imean = 0 => LHFLDIAG disabled'
        LHFLDIAG = .FALSE.
#ifdef __coupled
        IF ( LHFLDIAG ) WRITE(0,*)'coupled => LHFLDIAG disabled'
        LHFLDIAG = .FALSE.
#endif
        IF ( LCONVDIAG ) WRITE(0,*)'imean = 0 => LCONVDIAG disabled'
        LCONVDIAG = .FALSE.
        IF ( LDIFFDIAG ) WRITE(0,*)'imean = 0 => LDIFFDIAG disabled'
        LDIFFDIAG = .FALSE.
        IF ( LGRIDINFO ) WRITE(0,*)'imean = 0 => LGRIDINFO disabled'
        LGRIDINFO = .FALSE.
        IF ( LCALCDIFI ) WRITE(0,*)'imean = 0 => LCALCDIFI disabled'
        LCALCDIFI = .FALSE.
      END IF

      IF ( IBOLK .EQ. 0 )  LGMDIAG = .FALSE.
      IF ( IBOLK .EQ. 0 ) WRITE(0,*)'ibolk = 0 => GM scheme disabled'
      IF ( IBOLK .GT. 0 ) WRITE(0,*)'ibolk > 0 => GM scheme enabled'

      IF ( .NOT. LISOPYC ) THEN
        IF ( IBOLK .GT. 0 ) WRITE(0,*)                                  &
             'ATTN => GM scheme enabled with horizontal ',              &
             ' instead of isopycnal diffusion'
      ENDIF
      IF ( IBOLK .LT. 0 ) WRITE(0,*)                                  &
           'ibolk < 0 => GM + Visbek scheme enabled'
      IF (.NOT. LISOPYC) THEN
        IF ( IBOLK .LT. 0 ) WRITE(0,*)                                  &
             'ATTN => GM+Visbek scheme enabled with horizontal ',       &
             ' instead of isopycnal diffusion'
      ENDIF

      IF (iocad == 4) THEN
        WRITE(0,'(/,a,/,a,/)') ' iocad = 4 is obsolete => use iocad = 3 '// &
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
          WRITE (0,*) 'iocad =', iocad, '=> BBL transport scheme disabled'
          ibbl_transport = 0
        END IF
      CASE DEFAULT
        WRITE (0,*) 'ERROR: The specified iocad value', iocad, 'is not valid!'
        CALL p_abort
      END SELECT

      IF (iocaduv /= 3 .AND. iocaduv /= 8) THEN
        WRITE (0,*) 'ERROR: The specified iocaduv value', iocaduv, &
             'is not valid!'
        CALL p_abort
      ENDIF

      !HH   CHECK LAYER THICKNESS
      DO k=1, ke
        ! FIXME: shouldn't this also check less than zero?
        IF (dzw(k) .EQ. 0._wp) THEN
          WRITE(io_stdout, *) ' LAYER: ', k, ' THICKNESS IS ZERO !!!'
          CALL stop_all('check of layer thickness failed => run aborted')
        ELSE
          WRITE(io_stdout, *) k, ' LAYERTHICKNESS    (DZW): ', dzw(k)
        END IF
      END DO


    ENDIF ! (p_pe==p_io)


  END SUBROUTINE reset_inconsistent_user_setup




  SUBROUTINE broadcast_user_setup
#ifndef NOMPI
#ifndef MESSY
    USE mo_mpi, ONLY :  p_all_comm, p_real, p_int
#endif

    IMPLICIT NONE

    include "mpif.h"

    INTEGER (kind=MPI_ADDRESS_KIND)  :: position

    INTEGER :: ierr
    INTEGER, PARAMETER :: nbytes=50000
    CHARACTER :: scratch(nbytes)


    position = 0_MPI_ADDRESS_KIND

    IF ( nprocxy > 1 .AND. p_pe /= nprocxy ) THEN

      IF ( p_pe == p_io ) THEN

#ifndef MESSY
        CALL mpi_pack(dt,1,p_real,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(caulapts,1,p_real,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(caulapuv,1,p_real,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(cah00,1,p_real,scratch,nbytes,position,p_all_comm,ierr)

        CALL mpi_pack(nfixyearlen,1,p_int,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(aus,1,p_real,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(bofric,1,p_real,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(rayfric,1,p_real,scratch,nbytes,position,p_all_comm,ierr)

        CALL mpi_pack(av0,1,p_real,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(dv0,1,p_real,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(cwt,1,p_real,scratch,nbytes,position,p_all_comm,ierr)

        CALL mpi_pack(cwa,1,p_real,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(cstabeps,1,p_real,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(dback,1,p_real,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(aback,1,p_real,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(crelsal,1,p_real,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(creltem,1,p_real,scratch,nbytes,position,p_all_comm,ierr)

        CALL mpi_pack(numriv,1,p_int,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(numglac,1,p_int,scratch,nbytes,position,p_all_comm,ierr)

        CALL mpi_pack(iaufr,1,p_int,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(iaufw,1,p_int,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(nyears,1,p_int,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(nmonts,1,p_int,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(ndays,1,p_int,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(imean,1,p_int,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(ly_end,1,p_int,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(ly_start,1,p_int,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(lm_start,1,p_int,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(istart,1,p_int,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(i3drest,1,p_int,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(spongezone,6,p_int,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(spzndamp_time,1,p_real,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(ioasisflux,1,p_int,scratch,nbytes,position,p_all_comm,ierr)

        CALL mpi_pack(rrelsal,1,p_real,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(rreltem,1,p_real,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(cdvocon,1,p_real,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(cavocon,1,p_real,scratch,nbytes,position,p_all_comm,ierr)

        CALL mpi_pack(rleadclose,3,p_real,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(h0,1,p_real,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(hmin,1,p_real,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(armin,1,p_real,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(armax,1,p_real,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(hsntoice,1,p_real,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(sicthmin,1,p_real,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(sice,1,p_real,scratch,nbytes,position,p_all_comm,ierr)

        CALL mpi_pack(ibbl_transport,1,p_int,scratch,nbytes,position,p_all_comm,ierr)

        CALL mpi_pack(ibolk,1,p_int,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(isnflg,1,p_int,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(iocad,1,p_int,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(iocaduv,1,p_int,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(ioconv,1,p_int,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(icontro,1,p_int,scratch,nbytes,position,p_all_comm,ierr)

        CALL mpi_pack(ihalo_sor,1,p_int,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(imod,1,p_int,scratch,nbytes,position,p_all_comm,ierr)

        CALL mpi_pack(lisopyc,1,mpi_logical,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(lmpitype,1,mpi_logical,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(lnonblock,1,mpi_logical,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(ltidal,1,mpi_logical,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(ltidal_diag,1,mpi_logical,scratch,nbytes,position,p_all_comm,ierr)

        CALL mpi_pack(lswr_jerlov,1,mpi_logical,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(lfb_bgc_oce,1,mpi_logical,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(jerlov_atten,1,p_real,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(jerlov_bluefrac,1,p_real,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(lwith_one_layer_shelfs,1,mpi_logical,scratch,nbytes,position,p_all_comm,ierr)

        CALL mpi_pack(lforcediag,1,mpi_logical,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(lhfldiag,1,mpi_logical,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(lconvdiag,1,mpi_logical,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(ldiffdiag,1,mpi_logical,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(lgmdiag,1,mpi_logical,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(lgridinfo,1,mpi_logical,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(lcalcdifi,1,mpi_logical,scratch,nbytes,position,p_all_comm,ierr)

        CALL mpi_pack(itsdiag,1,p_int,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(ltstranspose,1,mpi_logical,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(ltswrite,1,mpi_logical,scratch,nbytes,position,p_all_comm,ierr)

        CALL mpi_pack(ladpo,1,mpi_logical,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(ladfs,1,mpi_logical,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(lsaoclose,1,mpi_logical,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(lzo_correct,1,mpi_logical,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(lundelayed_momentum_advection,1,mpi_logical,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(lwith_barotropic_stokes_drift,1,mpi_logical,scratch,nbytes,position,p_all_comm,ierr)

        CALL mpi_pack(fp_tracing_enabled,1,mpi_logical,scratch,nbytes,position,p_all_comm,ierr)

        CALL mpi_pack(iter_sor,1,p_int,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(iter_sor_hack,1,p_int,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(rtsorpar,1,p_real,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(rtsorpar_hack,1,p_real,scratch,nbytes,position,p_all_comm,ierr)

        CALL mpi_pack(dzw,ke,p_real,scratch,nbytes,position,p_all_comm,ierr)
        CALL mpi_pack(time_verbosity,1,p_int,scratch,nbytes,position,&
             p_all_comm,ierr)

#endif

      ENDIF

      CALL mpi_bcast(scratch,nbytes,mpi_packed,p_io,p_all_comm,ierr)

      IF ( p_pe /= p_io ) THEN

        ! unpack needs to be in the same sequence as pack !!

        position = 0_MPI_ADDRESS_KIND
#ifndef MESSY
        CALL mpi_unpack(scratch,nbytes,position,dt,1,p_real,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,caulapts,1,p_real,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,caulapuv,1,p_real,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,cah00,1,p_real,p_all_comm,ierr)

        CALL mpi_unpack(scratch,nbytes,position,nfixyearlen,1,p_int,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,aus,1,p_real,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,bofric,1,p_real,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,rayfric,1,p_real,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,av0,1,p_real,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,dv0,1,p_real,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,cwt,1,p_real,p_all_comm,ierr)

        CALL mpi_unpack(scratch,nbytes,position,cwa,1,p_real,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,cstabeps,1,p_real,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,dback,1,p_real,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,aback,1,p_real,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,crelsal,1,p_real,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,creltem,1,p_real,p_all_comm,ierr)

        CALL mpi_unpack(scratch,nbytes,position,numriv,1,p_int,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,numglac,1,p_int,p_all_comm,ierr)

        CALL mpi_unpack(scratch,nbytes,position,iaufr,1,p_int,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,iaufw,1,p_int,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,nyears,1,p_int,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,nmonts,1,p_int,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,ndays,1,p_int,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,imean,1,p_int,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,ly_end,1,p_int,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,ly_start,1,p_int,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,lm_start,1,p_int,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,istart,1,p_int,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,i3drest,1,p_int,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,spongezone,6,p_int,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,spzndamp_time,1,p_real,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,ioasisflux,1,p_int,p_all_comm,ierr)


        CALL mpi_unpack(scratch,nbytes,position,rrelsal,1,p_real,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,rreltem,1,p_real,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,cdvocon,1,p_real,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,cavocon,1,p_real,p_all_comm,ierr)

        CALL mpi_unpack(scratch,nbytes,position,rleadclose,3,p_real,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,h0,1,p_real,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,hmin,1,p_real,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,armin,1,p_real,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,armax,1,p_real,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,hsntoice,1,p_real,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,sicthmin,1,p_real,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,sice,1,p_real,p_all_comm,ierr)

        CALL mpi_unpack(scratch,nbytes,position,ibbl_transport,1,p_int,p_all_comm,ierr)

        CALL mpi_unpack(scratch,nbytes,position,ibolk,1,p_int,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,isnflg,1,p_int,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,iocad,1,p_int,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,iocaduv,1,p_int,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,ioconv,1,p_int,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,icontro,1,p_int,p_all_comm,ierr)

        CALL mpi_unpack(scratch,nbytes,position,ihalo_sor,1,p_int,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,imod,1,p_int,p_all_comm,ierr)

        CALL mpi_unpack(scratch,nbytes,position,lisopyc,1,mpi_logical,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,lmpitype,1,mpi_logical,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,lnonblock,1,mpi_logical,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,ltidal,1,mpi_logical,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,ltidal_diag,1,mpi_logical,p_all_comm,ierr)

        CALL mpi_unpack(scratch,nbytes,position,lswr_jerlov,1,mpi_logical,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,lfb_bgc_oce,1,mpi_logical,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,jerlov_atten,1,p_real,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,jerlov_bluefrac,1,p_real,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,lwith_one_layer_shelfs,1,mpi_logical,p_all_comm,ierr)

        CALL mpi_unpack(scratch,nbytes,position,lforcediag,1,mpi_logical,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,lhfldiag,1,mpi_logical,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,lconvdiag,1,mpi_logical,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,ldiffdiag,1,mpi_logical,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,lgmdiag,1,mpi_logical,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,lgridinfo,1,mpi_logical,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,lcalcdifi,1,mpi_logical,p_all_comm,ierr)

        CALL mpi_unpack(scratch,nbytes,position,itsdiag,1,p_int,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,ltstranspose,1,mpi_logical,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,ltswrite,1,mpi_logical,p_all_comm,ierr)

        CALL mpi_unpack(scratch,nbytes,position,ladpo,1,mpi_logical,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,ladfs,1,mpi_logical,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,lsaoclose,1,mpi_logical,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,lzo_correct,1,mpi_logical,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,lundelayed_momentum_advection,1,mpi_logical,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,lwith_barotropic_stokes_drift,1,mpi_logical,p_all_comm,ierr)

        CALL mpi_unpack(scratch,nbytes,position,fp_tracing_enabled,1,mpi_logical,p_all_comm,ierr)

        CALL mpi_unpack(scratch,nbytes,position,iter_sor,1,p_int,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,iter_sor_hack,1,p_int,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,rtsorpar,1,p_real,p_all_comm,ierr)
        CALL mpi_unpack(scratch,nbytes,position,rtsorpar_hack,1,p_real,p_all_comm,ierr)

        CALL mpi_unpack(scratch,nbytes,position,dzw,ke,p_real,p_all_comm,ierr)

        CALL mpi_unpack(scratch,nbytes,position,time_verbosity,1,p_int,&
             p_all_comm,ierr)
#endif

      ENDIF
    ENDIF

    CALL broadcast_time_desc(model_start_time)


    !    CALL p_bcast(DT,p_io)
    !    CALL p_bcast(CAULAPTS,p_io)
    !    CALL p_bcast(CAULAPUV,p_io)
    !    CALL p_bcast(CAH00,p_io)

    !    CALL p_bcast(nfixYearLen,p_io)
    !    CALL p_bcast(AUS,p_io)
    !    CALL p_bcast(bofric,p_io)
    !    CALL p_bcast(rayfric,p_io)

    !    CALL p_bcast(AV0,p_io)
    !    CALL p_bcast(DV0,p_io)
    !    CALL p_bcast(CWT,p_io)

    !    CALL p_bcast(CWA,p_io)
    !    CALL p_bcast(CSTABEPS,p_io)
    !    CALL p_bcast(DBACK,p_io)
    !    CALL p_bcast(ABACK,p_io)
    !    CALL p_bcast(CRELSAL,p_io)
    !    CALL p_bcast(CRELTEM,p_io)

    !    CALL p_bcast(NUMRIV,p_io)
    !    CALL p_bcast(NUMGLAC,p_io)

    !    CALL p_bcast(IAUFR,p_io)
    !    CALL p_bcast(IAUFW,p_io)
    !    CALL p_bcast(NYEARS,p_io)
    !    CALL p_bcast(NMONTS,p_io)
    !    CALL p_bcast(NDAYS,p_io)
    !    CALL p_bcast(IMEAN,p_io)
    !    CALL p_bcast(LY_END,p_io)
    !    CALL p_bcast(LY_START,p_io)
    !    CALL p_bcast(LM_START,p_io)
    !    CALL p_bcast(ISTART,p_io)
    !    CALL p_bcast(I3DREST,p_io)
    !    CALL p_bcast(IOASISFLUX,p_io)


    !    CALL p_bcast(RRELSAL,p_io)
    !    CALL p_bcast(RRELTEM,p_io)
    !    CALL p_bcast(CDVOCON,p_io)
    !    CALL p_bcast(CAVOCON,p_io)

    !    CALL p_bcast(rleadclose,p_io)
    !    CALL p_bcast(H0,p_io)
    !    CALL p_bcast(HMIN,p_io)
    !    CALL p_bcast(ARMIN,p_io)
    !    CALL p_bcast(ARMAX,p_io)
    !    CALL p_bcast(HSNTOICE,p_io)
    !    CALL p_bcast(SICTHMIN,p_io)
    !    CALL p_bcast(SICE,p_io)

    !    CALL p_bcast(ibbl_transport,p_io)


    !    CALL p_bcast(IBOLK,p_io)
    !    CALL p_bcast(ISNFLG,p_io)
    !    CALL p_bcast(IOCAD,p_io)
    !    CALL p_bcast(IOCADUV,p_io)
    !    CALL p_bcast(IOCONV,p_io)
    !    CALL p_bcast(ICONTRO,p_io)

    !    CALL p_bcast(ihalo,p_io)
    !    CALL p_bcast(imod,p_io)

    !    CALL p_bcast(LISOPYC,p_io)
    !    CALL p_bcast(lmpitype,p_io)
    !    CALL p_bcast(lnonblock,p_io)
    !    CALL p_bcast(ltidal,p_io)
    !    CALL p_bcast(ltidal_diag,p_io)

    !    CALL p_bcast(lswr_jerlov,p_io)
    !    CALL p_bcast(lfb_bgc_oce,p_io)
    !    CALL p_bcast(jerlov_atten,p_io)
    !    CALL p_bcast(jerlov_bluefrac,p_io)

    !    CALL p_bcast(lwith_one_layer_shelfs,p_io)

    !    CALL p_bcast(LFORCEDIAG,p_io)
    !    CALL p_bcast(LHFLDIAG,p_io)
    !    CALL p_bcast(LCONVDIAG,p_io)
    !    CALL p_bcast(LDIFFDIAG,p_io)
    !    CALL p_bcast(LGMDIAG,p_io)
    !    CALL p_bcast(LGRIDINFO,p_io)
    !    CALL p_bcast(LCALCDIFI,p_io)

    !    CALL p_bcast(ITSDIAG,p_io)
    !    CALL p_bcast(LTSTRANSPOSE,p_io)
    !    CALL p_bcast(LTSWRITE,p_io)


    !    CALL p_bcast(ladpo,p_io)
    !    CALL p_bcast(ladfs,p_io)
    !    CALL p_bcast(lsaoclose,p_io)
    !    CALL p_bcast(lzo_correct,p_io)
    !    CALL p_bcast(lundelayed_momentum_advection,p_io)
    !    CALL p_bcast(fp_tracing_enabled, p_io)


    !    CALL p_bcast(iter_sor,p_io)
    !    CALL p_bcast(iter_sor_hack,p_io)
    !    CALL p_bcast(rtsorpar,p_io)
    !    CALL p_bcast(rtsorpar_hack,p_io)

    !    CALL p_bcast(dzw,p_io)


#endif


  END SUBROUTINE broadcast_user_setup




  SUBROUTINE logwrite_user_setup
    ! Write final namelist parameters
    !
    WRITE(IO_STDOUT,*)' SECONDS PER TIMESTEP   (DT): ',DT
    WRITE(IO_STDOUT,*)' TIME STEPS PER DAY (NDTDAY): ',NDTDAY
    WRITE(IO_STDOUT,*)'         (FORCING_FREQUENCY): ',FORCING_FREQUENCY
    WRITE(IO_STDOUT,*)'               (NFIXYEARLEN): ',nfixYearLen
    WRITE(IO_STDOUT,*)'                     (IOCAD): ',IOCAD
    WRITE(IO_STDOUT,*)'                   (IOCADUV): ',IOCADUV
    WRITE(IO_STDOUT,*)'            (IBBL_TRANSPORT): ',IBBL_TRANSPORT
    WRITE(IO_STDOUT,*)'                    (IOCONV): ',IOCONV
    WRITE(IO_STDOUT,*)'                  (CAULAPUV): ',CAULAPUV
    WRITE(IO_STDOUT,*)'                  (CAULAPTS): ',CAULAPTS
    WRITE(IO_STDOUT,*)'                       (AUS): ',AUS
    WRITE(IO_STDOUT,*)'                    (BOFRIC): ',BOFRIC
    WRITE(IO_STDOUT,*)'                   (RAYFRIC): ',RAYFRIC
    WRITE(IO_STDOUT,*)'                     (CAH00): ',CAH00
    WRITE(IO_STDOUT,*)'                       (AV0): ',AV0
    WRITE(IO_STDOUT,*)'                       (DV0): ',DV0
    WRITE(IO_STDOUT,*)'                       (CWT): ',CWT
    WRITE(IO_STDOUT,*)'                       (CWA): ',CWA
    WRITE(IO_STDOUT,*)'                  (CSTABEPS): ',CSTABEPS
    WRITE(IO_STDOUT,*)'                     (DBACK): ',DBACK
    WRITE(IO_STDOUT,*)'                     (ABACK): ',ABACK
    WRITE(IO_STDOUT,*)'                   (CRELSAL): ',CRELSAL
    WRITE(IO_STDOUT,*)'                   (CRELTEM): ',CRELTEM
    WRITE(IO_STDOUT,*)'                    (NUMRIV): ',NUMRIV
    WRITE(IO_STDOUT,*)'                   (NUMGLAC): ',NUMGLAC
    WRITE(IO_STDOUT,*)'                   (CDVOCON): ',CDVOCON
    WRITE(IO_STDOUT,*)'                   (CAVOCON): ',CAVOCON
    WRITE(IO_STDOUT,*)'                     (IBOLK): ',IBOLK
    WRITE(IO_STDOUT,*)'                   (LISOPYC): ',LISOPYC
    WRITE(IO_STDOUT,*)'                     (IMEAN): ',IMEAN
    WRITE(IO_STDOUT,*)'                  (LY_START): ',LY_START
    WRITE(IO_STDOUT,*)'                  (LM_START): ',LM_START
    WRITE(IO_STDOUT,*)'                    (LY_END): ',LY_END
    WRITE(IO_STDOUT,*)'                    (ISTART): ',ISTART
    WRITE(IO_STDOUT,*)'                   (I3DREST): ',I3DREST
    WRITE(IO_STDOUT,*)'                   (ICONTRO): ',ICONTRO
    WRITE(IO_STDOUT,*)'                (IOASISFLUX): ',IOASISFLUX
    WRITE(IO_STDOUT,*)'                   (LGMDIAG): ',LGMDIAG
    WRITE(IO_STDOUT,*)'                (LFORCEDIAG): ',LFORCEDIAG
    WRITE(IO_STDOUT,*)'                  (LHFLDIAG): ',LHFLDIAG
    WRITE(IO_STDOUT,*)'                 (LDIFFDIAG): ',LDIFFDIAG
    WRITE(IO_STDOUT,*)'                 (LCONVDIAG): ',LCONVDIAG
    WRITE(IO_STDOUT,*)'                 (LGRIDINFO): ',LGRIDINFO
    WRITE(IO_STDOUT,*)'                 (LCALCDIFI): ',LCALCDIFI
    WRITE(IO_STDOUT,*)'                   (ITSDIAG): ',ITSDIAG
    WRITE(IO_STDOUT,*)'              (LTSTRANSPOSE): ',LTSTRANSPOSE
    WRITE(IO_STDOUT,*)'                  (LMPITYPE): ',LMPITYPE
    WRITE(IO_STDOUT,*)'                 (LNONBLOCK): ',LNONBLOCK
    WRITE(IO_STDOUT,*)'           (LBOUNDS_EXCH_TP): ',LBOUNDS_EXCH_TP
    WRITE(IO_STDOUT,*)'                    (LTIDAL): ',LTIDAL
    WRITE(IO_STDOUT,*)'               (LTIDAL_DIAG): ',LTIDAL_DIAG

    WRITE(IO_STDOUT,*)'               (LSWR_JERLOV): ',LSWR_JERLOV
    WRITE(IO_STDOUT,*)'               (LFB_BGC_OCE): ',LFB_BGC_OCE
    WRITE(IO_STDOUT,*)'               (LSWR_ATTEN): ',JERLOV_atten
    WRITE(IO_STDOUT,*)'               (LSWR_BLUEFRAC): ',JERLOV_bluefrac

    WRITE(IO_STDOUT,*)'(LUNDELAYED_MOMENTUM_ADVECTION): ',LUNDELAYED_MOMENTUM_ADVECTION
    WRITE(IO_STDOUT,*)'(WITH_BAROTROPIC_STOKES_DRIFT): ',LWITH_BAROTROPIC_STOKES_DRIFT


#ifdef SOR
    WRITE(IO_STDOUT,*)'                 (IHALO_SOR): ',IHALO_SOR
    WRITE(IO_STDOUT,*)'                      (IMOD): ',IMOD
    WRITE(IO_STDOUT,*)'                  (RTSORPAR): ',RTSORPAR
    WRITE(IO_STDOUT,*)'          (RTSORPAR_HACK): ',RTSORPAR_HACK
    WRITE(IO_STDOUT,*)'                  (ITER_SOR): ',ITER_SOR
    WRITE(IO_STDOUT,*)'          (ITER_SOR_HACK): ',ITER_SOR_HACK
#endif
    IF (.NOT. istart_in_valid_range(istart))                              &
         WRITE(IO_STDOUT,*)'ISTART NOT SUPPORTED!!!!!'
    IF (I3DREST .LT. 0 .OR. I3DREST .GT. 2)                            &
         WRITE(IO_STDOUT,*)'I3DREST NOT SUPPORTED!!!!!'
    WRITE(IO_STDOUT,*)'THIS JOB WILL TRY TO INTEGRATE '               &
         ,NYEARS,' YEARS AND '                           &
         ,NMONTS,' MONTHS AND'                           &
         ,NDAYS,'DAYS'

  END SUBROUTINE logwrite_user_setup




  SUBROUTINE setup_time_invariants

    CALL init_levitus_2d

    IF (i3drest .GT. 0 .OR. istart .LT. istart_run_continue) THEN
      CALL init_levitus_3d
    END IF

    IF (istart .LT. istart_run_continue) iaufr = 0

    aulapts = caulapts * dt / 3600._wp
    aulapuv = caulapuv * dt / 3600._wp
    ah00 = cah00 / 4.e5_wp

    tiestw(1) = 0._wp
    DO k=1,ke
      tiestw(k+1) = tiestw(k) + dzw(k)
    END DO

    CALL setup_ocean_vertical_mixing(cwt, cwa, aback, dback)

    !  CHECK RIVER RUNOFF SCHEME
    IF ( numriv > 0 ) THEN
      luse_river_runoff_stations=.TRUE.
      WRITE(io_stdout,*)'river_runoff_stations enabled -> omip rivers disabled'
    ENDIF

    !  CHECK GLACIER CALVING SCHEME
    IF ( numglac > 0 ) THEN
      luse_glac_calv=.TRUE.
      WRITE(io_stdout,*)'glacier calving enabled'
    ENDIF

    !HH    RELAXATION TIME SALINITY
    !      RELSAL = CRELSAL*20./DZW(1)
    IF (crelsal .GT. almzer) THEN
      !         WRITE(IO_STDOUT,26668)1./(crelsal*24.*3600.)
      IF (p_pe == p_io) THEN
        WRITE(0,*) 'RELAXATION DZW(1) [m]  =', DZW(1)
        ! FIXME: put tsolar_day here
        WRITE(0,*) ' RELAXATION TIME [DAYS] =', 1._wp / (crelsal * 24._wp * 3600._wp)
        WRITE(0,*) ' PISTON VELOCITY [m/s]  =', dzw(1) * crelsal
        WRITE(0,*) ' RELAXATION TIME (relative to 20m) [DAYS] =', &
             1._wp / (crelsal * 24._wp * 3600._wp * dzw(1)/20._wp)
      END IF
    ELSE
      IF (p_pe == p_io) THEN
        WRITE(IO_STDOUT,*) 'SSS relaxation switched off !!'
      END IF
    END IF

    !26668 FORMAT('  RELAXATION TIME SALINITY COUPLING : ',F10.2,' DAYS')
    !HH    RELAXATION TIME TEMPERATURE
    !      RELTEM = CRELTEM*20./DZW(1)
    IF (creltem .GT. almzer) THEN
      !         WRITE(IO_STDOUT,26669)1./(creltem * 24.*3600.)
      IF (p_pe == p_io) THEN
        WRITE(IO_STDOUT,*)'  RELAXATION TIME TEMPEARTURE COUPLING : '&
             ,DZW(1),' m /', (creltem * 24._wp * 3600._wp),' DAYS'
      END IF
    ELSE
      IF (p_pe == p_io) THEN
        WRITE(IO_STDOUT,*) 'SST relaxation switched off !!'
      END IF
    END IF
    !26669 FORMAT('  RELAXATION TIME SST COUPLING : ',F10.2,' DAYS')


    !-----------------------------------------------------------------------

    CALL alloc_mem_para2

    IF ( lwith_barotropic_stokes_drift ) THEN
      CALL alloc_mem_stokes_drift
    ENDIF

    CALL alloc_mem_octdiff

    IF (ladpo) THEN
      CALL alloc_mem_adpo
    END IF
    IF (ibbl_transport .EQ. 1) THEN
      CALL alloc_mem_commobbl
    END IF

    IF (IBOLK .NE. 0) THEN
      CALL alloc_mem_gmbolus
    END IF

    !    write(0,*)'vor mean'

    IF (imean.NE.0) THEN
      CALL alloc_mem_mean
    END IF

#if defined (__coupled)
    !----------------------------------------------------------------------
    ! Prepare coupling
    !
    CALL couple_prep
#endif

    !-------------------------------------------------------------------

    !     DU/DT - F*V = -G*( STABN* (DZ/DX)(NEW) + STABO* (DZ/DX)(OLD) )

    !     DZ/DT = H*( CONN* (DU/DX)(NEW) + CONO* (DU/DX)(OLD) + ... )

    conn = 0.50_wp
    cono = 1._wp - conn

    WRITE(IO_STDOUT,*)' STABN = ',STABN,' CONN = ',CONN
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

    CALL BELEG
    !      PRINT*,'nach beleg'
    !
    CALL BODEN
    !      PRINT*,'nach boden'
    !
    !-------------------------------------------------------------------
    !OtB LMONTS used in LEVIRE but not defined yet
    LMONTS=0
    IF (I3DREST .GT. 0) THEN
      ! read 3d levitus data for use in 3d restoring
      CALL levitus_read_3d_restore
    END IF ! I3DREST
    !
    !
#ifdef NOLAND
    CALL find_wet_proc
#endif
    !-------------------------------------------------------------------
    ! SPECIFY CORIOLIS PARAMETER ETC.
    !
    CALL CORIOL
    !
    !
    !
    !-------------------------------------------------------------------
    ! DIAGNOSTICS : MASKS OF BASINS
    !               9 : GLOBAL
    !
    CALL load_basin_masks(istart)

    !-------------------------------------------------------------------
    !            Print*,'nach Zaehl'
    !-------------------------------------------------------------------
    IF (istart .LT. istart_run_continue) THEN
      CALL levitus_set_ts
    END IF !ISTART
    IF (ibbl_transport .EQ. 1) THEN
      ! PREPARE FOR BOTTOM BOUNDARY PARAMETRIZATIONS
      CALL FINDBOT(KBOT,WETO,IE,JE,KE)
      CALL FINDALFA
      !    PRINT*,'Preparation for bottom boundary parametrizations completed'
    END IF
    !-------------------------------------------------------------------
    !  provide matrix for barotropic mode
    !     uwe  add iteration of matrix
    !
    IF (ielimi .GE. 1) THEN
      CALL trian
    ELSE
      CALL itprep
      CALL trotest2
    END IF
    !      print*,'nach trotest'
    !-------------------------------------------------------------------
    IF ( luse_river_runoff_stations ) THEN
      CALL river_runoff_stations_ini
    ELSE
      CALL river_runoff_omip_ini
    ENDIF

    CALL alloc_mem_swr_absorption

    IF (time_verbosity > 2) tttr = p_time()
    ldt = -1
    CALL initialize_surface_forcing
    IF (time_verbosity > 2) &
         time_diffs_subr(tag_subr_surf_forcing) = p_time() - tttr

    IF ( luse_glac_calv ) THEN
      CALL glac_calv_ini
    ENDIF

    IF ( ltidal ) THEN
      CALL alloc_mem_tidal
      IF ( ltidal_diag ) THEN
        CALL init_tide_timeseries ! Initialize tide outdput positions
      END IF
    END IF

    ! INITIALIZE DIAGNOSTICS

    CALL DIAG_INI

    !   write(0,*)'nach diag_ini'

    ! Initialize time computations.
    call model_time_init(int(dt), nfixYearLen)

#ifndef NO_NEW_IO
    CALL build_ocean_varlist
    CALL iolist_create_file_list(ocean_iolist, ocean_varlist)
#endif

    CALL register_diagnostics

    ! INCLUDE WRITEOUT OF CPP OPTIONS AND PARAMETER SETTINGS
    CALL CPPOUTPUT
    CALL setup_time_stepping

    IF ( ltidal ) THEN
      CALL foreph_ini
    END IF

    IF ( .NOT. lswr_jerlov ) THEN
      CALL old_swr_absorption
    ELSE
      CALL jerlov_swr_absorption
    ENDIF

  END SUBROUTINE setup_time_invariants



  SUBROUTINE setup_time_stepping
    !-----------------------------------------------------------------------
    !
    !     TIMESTEPPING IS BASED ON 3 LOOPS
    !        PER_YEAR      : DO  LYEAR = LYEAR1, LYEAR2
    !        MONTH_OF_YEAR :   DO  LMONT = 1 OR LMONT1, LMONT2 OR 12
    !        DAY_OF_MONTH  :     DO  LDAY = 1,MONLEN(LMONTS,LYEARS)
    !        TIME_STEP_OF_DAY   :       DO  LDTDAY = 1,NDTDAY
    ! ADDITIONAL COUNTERS
    !        LYEARS = ACTUAL YEAR
    !        LMONTS = ACTUAL MONTH
    !        LDAYS  = ACTUAL DAY
    !
    !     TIME COUNTER ARE :
    !               LDT       : TIME STEP COUNTER FOR EXPERIMENT
    !               LDTRUN    : TIME STEP COUNTER FOR RUN
    !               LDTYEAR   : TIME STEP COUNTER FOR YEAR
    !               LDTMONTH  : TIME STEP COUNTER FOR MONTH
    !               NDTDAY  : TIME STEPS PER DAY

    TYPE(time_desc) :: dtime

    WRITE(IO_STDOUT,*)'START FROM PREVIOUS CALCULATION YES=1/NO=0 : ' &
         ,IAUFR
    WRITE(IO_STDOUT,*)'WRITE RESTART FILE              YES=1/NO=0 : ' &
         ,IAUFW


    !-----------------------------------------------------------------------
    !
    !     START FROM STATUS LEVITUS OR HORIZONTAL STRATIFICATION
    !
    CALL levitus_read_surface_salinity
    IF (creltem .GT. almzer) THEN
      CALL levitus_read_surface_temperature
    END IF

    CALL setup_restart(iaufr)

    IF (IAUFR .EQ. 0) THEN
      IF (istart .EQ. istart_new_topo_update) &
           CALL levitus_read_3d_stratification
      IF (istart .EQ. istart_new_run) CALL levitus_horizontal_stratification
      IF (istart .EQ. istart_new_topo) THEN

        CALL levitus_read_3d_stratification
        DO J=1,JE
          DO I=1,IE
            IF (tho(i, j, 1) * weto(i, j, 1) .LT. -0.5_wp) THEN
              sictho(i, j) = MAX(0._wp, 3._wp * (-0.5_wp - tho(i, j, 1)))
              sicomo(i, j) = 1._wp
            END IF
          END DO
        END DO
      END IF ! ISTART

      !  SET AGE OF LEVITUS DATA
      LYEARS = model_start_time%year
      LMONTS = model_start_time%month
      LDAYS = model_start_time%mday
      LDT = 0
      WRITE(IO_STDOUT,*)'TIME COUNTER SET TO ZERO '
#ifdef OLD_IO
      !     WRITE RESTART FILES WITH INITIAL STRATIFICATION
      !     MAKE SURE THAT BOTH RESTART FILES EXIST ON EXIT
      IF (IAUFW .EQ. 1) THEN
        CALL AUFW(ldt)
        CALL AUFW(ldt)
      END IF
#endif
    END IF

    !
    !-----------------------------------------------------------------------
    !     START FROM RESTART FILES Z37000 OR Z38000
    IF (IAUFR .EQ. 1) THEN
#ifdef _PROFILE
      CALL trace_start ('restartread', 10)
#endif
#ifdef OLD_IO
      ! Old I/O is overridden by new I/O, but we chose to run both as a test.
      CALL AUFR(ldt)
#endif
#ifndef NO_NEW_IO
      ! May override old I/O (see above).
      CALL readRestartCDI(ocean_iolist,ocean_varlist)
      CALL varlist_post_process_ocean_data
      WRITE(0,*)'read restart date',LYEARS,LMONTS,LDAYS
#endif
#ifdef _PROFILE
      CALL trace_stop ('restartread', 10)
#endif
    END IF

    !  CHECK GLOBAL SALT CONTENT
    IF ( ICONTRO .NE. 0 ) THEN
      CALL CONTRO(-999)
    END IF

    !-----------------------------------------------------------------------
    ! SET LIMITS OF YEAR/MONTH TIME STEPPING LOOP

    LDTRUN = 0
    LYEAR1 = LYEARS
    LMONT1 = LMONTS
    LDAY1  = LDAYS
    IF (iaufr .EQ. 1) CALL inc_day(lyear1, lmont1, lday1)

    run_start_time = time_desc(lyear1,lmont1,lday1,0,0,0)
    dtime = time_desc(nyears,nmonts,ndays,0,0,0)

    next_run_start_time = run_start_time + dtime
    run_end_time = next_run_start_time - INT(dt,i8)

    lyear2 = run_end_time%year
    lmont2 = run_end_time%month
    IF (lyear1 .EQ. lyear2 .AND. lmont1 .EQ. lmont2) THEN
      endday = run_end_time%mday
    ELSE
      endday = monlen(lmont1,lyear1)
    ENDIF

#ifdef PBGC
    !
    !--------------------------------------------------------------------
    ! Initialize bgc modules
    !
    !
    !  Length of run in model time steps:
    ndtrun=NINT(REAL(seconds_between_two_times(run_end_time, &
                     run_start_time),dp)/dt) + 1

    WRITE(IO_STDOUT,*)'Total no. of time steps this run: ',ndtrun

    CALL INI_BGC(iaufr,icycli,dt,ndtrun,ie,je,ke                    &
         ,ddpo,tho,sao,dlxp,dlyp,tiestu,tiestw         &
         ,lyears,lmonts,ldays,ldt,nyears,nmonts        &
         ,gila,giph)

    IF (ladfs) THEN
      WRITE(io_stdout,*) ' '
      WRITE(io_stdout,*) ' WARNING:'
      WRITE(io_stdout,*) ' iocad == ', iocad, ' i.e. advection routine OCADFS', &
           ' is chosen but BGC Model HAMOCC is active!'
      WRITE(io_stdout,*) ' There is NO ADVECTION of BGC tracers'
    END IF
#endif /*PBGC*/

#if defined (__coupled)
    !--------------------------------------------------------------------
    !  Initialise the coupled model :
    !     nactyear = first ocean year of run (from restart data) at entry
    !     nactyear = first year of actual date (from KONTROL.h)  at exit
    !
    nactyear = lyear1
    CALL couple_init(nactyear)
    !--------------------------------------------------------------------
#endif /*(__coupled)*/

    ! Write book-keeping files storing setup information for later use.

#ifndef NO_NEW_IO
    call iolist_write_initial(ocean_iolist, ocean_varlist, run_start_time)
#ifdef PBGC
    call iolist_write_initial(bgc_iolist, bgc_varlist, run_start_time)
#endif
#endif

    ! Print summary.

    IF (nfixYearLen.EQ.-1) THEN
      WRITE(IO_STDOUT,*)                                                &
           'THIS RUN ASSUMES REALISTIC YEAR WITH LEAP FEBRUARIES'
    ENDIF
    IF (nfixYearLen.EQ.360) THEN
      WRITE(IO_STDOUT,*)                                                &
           'THIS RUN ASSUMES IDEALIZED YEARS (12*30 DAYS EACH)'
    ENDIF
    IF (nfixYearLen.EQ.365) THEN
      WRITE(IO_STDOUT,*)                                                &
           'THIS RUN ASSUMES IDEALIZED YEARS WITH 365 DAYS'
    ENDIF

    WRITE(IO_STDOUT,*)                                                &
         ' INTEGRATION PERIOD IS YEAR ',LYEAR1,' TO ',LYEAR2

  END SUBROUTINE setup_time_stepping




  SUBROUTINE per_year_setup
    LYEARS = LYEAR
    LDTYEAR = 0

    WRITE(IO_STDOUT,*) ' TO BE CALCULATED NOW : YEAR = ',LYEARS

#ifndef __coupled
    CALL spool_forcing
#endif/*ndef __coupled */

    !--------------------------------------------------------------------
    !  SET MONTH LOOP LIMITS FOR EACH YEAR
    LMON1 = 1
    LMON2 = 12
    IF (LYEAR1 .EQ. LYEAR2) THEN
      LMON1 = LMONT1
      LMON2 = LMONT2
    ELSE
      IF (LYEARS .EQ. LYEAR2) THEN
        LMON1 = 1
        LMON2 = LMONT2
      END IF
      IF (LYEARS .EQ. LYEAR1) THEN
        LMON1 = LMONT1
        LMON2 = 12
      END IF
    END IF
    WRITE(IO_STDOUT,*) 'MONTHLY LOOP FROM MONTH ',LMON1,          &
         ' TO', LMON2
  END SUBROUTINE per_year_setup



  SUBROUTINE per_month_setup
    LMONTS = LMONT
    LDTMONTH = 0

    WRITE(IO_STDOUT, 18823) LYEARS, LMONTS
18823 FORMAT(' TO BE CALCULATED NOW : YEAR ',I4,' MONTH ',I4)

    CALL levitus_per_month_setup(lmont)

    IF (lyears .GT. lyear1 .OR. lmonts .GT. lmont1) THEN
      lday1 = 1
      endday = monlen(lmonts,lyears)
      IF (lmonts .EQ. lmont2 .AND. lyears .EQ. lyear2) &
                   endday=run_end_time%mday
    ENDIF

  END SUBROUTINE per_month_setup




  SUBROUTINE per_timestep_computation
#ifdef PBGC
    USE mo_planetary_constants, ONLY : rhoicwa, rhosnwa
#endif

    LDTDAYC = LDTDAY

    IF (time_verbosity > 0) &
         ttts = p_time()

    !*       6.  TIMESTEP THE OCEAN.
    !            -------- --- ------
    LDT = LDT + 1
    LDTRUN = LDTRUN + 1
    LDTYEAR = LDTYEAR + 1
    LDTMONTH =LDTMONTH + 1

    model_time%year = lyears
    model_time%month = lmonts
    model_time%mday = ldays
    model_time%hour = INT(REAL(ldtday - 1, dp) * dt / 3600.0_dp)
    model_time%minute = MODULO(INT(REAL(ldtday-1, dp) * dt) / 60, 60)
    model_time%second = MODULO(INT(REAL(ldtday-1, dp) * dt), 60)

    !  current date: yyyy-mm-dd hh:mm:ss (W3C XML Schema)
    IF (p_pe==p_io) THEN
      WRITE(0,'(a,i3,a,i7)') &
           'MPIOM: current date: '//TRIM(format_model_time(model_time))// &
           '  begin of timestep (day):', ldtday,'  (run):', ldtrun
    END IF

    CALL update_surface_forcing(model_time, model_start_time, ldtrun)

#ifdef _PROFILE
    CALL trace_start ('dt_loop', 7)
#endif

#ifdef PBGC
    IF (time_verbosity > 2) tttr = p_time()
    !--------------------------------------------------------------------

    CALL store_tracer(ddpo)

    ! CALL BIOGEOCHEMISTRY MODULES.
    DO J=1,JE
      DO I=1,IE
        bgcdpio(i, j, 1) = 0._wp
        bgcddpo(i, j, 1) = 0._wp
        !        LAYER1_BGC(I,J)=0.
        IF (ddpo(i, j, 1) .GT. 0._wp)THEN
          BGCDDPO(I,J,1)= DDPO(I,J,1)                                  &
               +ZO(I,J)-SICTHO(I,J)*RHOICWA-SICSNO(I,J)*RHOSNWA
          bgcdpio(i, j, 1) = 1._wp / bgcddpo(i, j, 1)
          !          LAYER1_BGC(I,J)=BGCDDPO(I,J,1)
        END IF
      END DO
    END DO
    DO K=2,KE
      DO J=1,JE
        DO I=1,IE
          BGCDDPO(I,J,K)= DDPO(I,J,K)
          BGCDPIO(I,J,K)= DPIO(I,J,K)
        END DO
      END DO
    END DO
#  ifdef __coupled
    CALL BGC(IE,JE,KE,                                               &
         AOFLSHWO,SICOMO,THO,SAO,BGCDDPO,DLXP,DLYP,              &
         TIESTU,BGCDPIO,AOFLWSVO,WO,fslp,LYEARS,LMONTS,LDAYS,     &
         MONLEN(LMONTS,lyears),LDTMONTH,LDTDAY)
#  else

    CALL bounds_exch(1,'p',FSWR,'mpiom 19')
    CALL bounds_exch(1,'p',SICOMO,'mpiom 20')
    CALL bounds_exch(1,'p',THO,'mpiom 20')
    CALL bounds_exch(1,'p',SAO,'mpiom 21')
    CALL bounds_exch(1,'p',BGCDDPO,'mpiom 22')
    CALL bounds_exch(1,'p',DLXP,'mpiom 23')
    CALL bounds_exch(1,'p',DLYP,'mpiom 24')
    CALL bounds_exch(1,'p',FU10,'mpiom 25')
    CALL bounds_exch(1,'p',WO,'mpiom 26')
    CALL bounds_exch(1,'p',fslp,'mpiom 27')
    CALL bounds_exch(1,'p',BGCDPIO,'mpiom 28')




    CALL BGC(IE,JE,KE,                                               &
         FSWR,SICOMO,THO,SAO,BGCDDPO,DLXP,DLYP,                  &
         TIESTU,BGCDPIO,FU10,WO,fslp,LYEARS,LMONTS,LDAYS,         &
         MONLEN(LMONTS,lyears),LDTMONTH,LDTDAY)
#  endif


    IF (time_verbosity > 2) time_diffs_subr(tag_subr_bgc) = p_time() - tttr

#endif /*PBGC*/

    !--------------------------------------------------------------------
    ! OCTHER : THERMODYNAMIC FORCING

    IF (time_verbosity > 2) tttr = p_time()

    IF ( ltidal ) THEN
      CALL FOREPH
    END IF

    IF ( ICONTRO .NE. 0 ) THEN
      CALL CONTRO(-10)
    END IF

    CALL OCICE


    IF ( ICONTRO .NE. 0 ) THEN
      CALL CONTRO(-11)
    END IF

    !      GOTO 999
    CALL OCTHER

    !      WRITE(0,*) 'WITHOUT OCTHER !!!'


    IF (time_verbosity > 2) time_diffs_subr(tag_subr_octher) = p_time() - tttr
123 FORMAT(A,F10.3)
    !--------------------------------------------------------------------


#ifdef CORE
    WRITE(IO_STDOUT,*) 'CALL NORM PEM'
    CALL NORMPEM
#endif



    IF (IBOLK .LT. 0) THEN
      ! UPDATE GM COEFFICIENTS ONCE/DAY !
      IF (LDTDAY .EQ. 1) CALL CALCGMVIS
    END IF


    IF ( ICONTRO .NE. 0 ) THEN
      CALL CONTRO(1)
    END IF
    !--------------------------------------------------------------------
    ! OCWIND :

    CALL OCWIND
    !      PRINT*,'NACH OCWIND'

    !--------------------------------------------------------------------
    IF (ltidal) THEN
      !  TIDES
      CALL TIPOUV
    END IF

    !--------------------------------------------------------------------
    !  OCTIMF : UPDATE VELOCITY
    !
    CALL OCTIMF
    !--------------------------------------------------------------------
    !  OCMODMOM : DECOMPOSITION INTO BAROTROPIC AND BAROCLINIC FIELD

    IF (time_verbosity > 2) tttr = p_time()
    CALL OCMODMOM

    IF ( lwith_barotropic_stokes_drift ) THEN
      CALL itprep
    ENDIF
    IF (time_verbosity > 2) time_diffs_subr(tag_subr_ocmodmom) = p_time() - tttr

    !--------------------------------------------------------------------
    IF ( .NOT. lwith_barotropic_stokes_drift ) THEN
      CALL OCBARP
    ENDIF
    !      PRINT*,'NACH OCBARP'
    !
    !--------------------------------------------------------------------

    IF (time_verbosity > 2) tttr = p_time()

    IF ( lwith_barotropic_stokes_drift ) THEN
      CALL troneu
    ENDIF

    CALL OCCLIT

    IF (time_verbosity > 2) time_diffs_subr(tag_subr_occlit) = p_time() - tttr

    !--------------------------------------------------------------------

    IF (time_verbosity > 2) tttr = p_time()

    IF ( ICONTRO .NE. 0 ) THEN
      CALL CONTRO(38)
    END IF

    IF(IELIMI.GE.1) THEN
      CALL BARTIM
    ELSE
#ifdef SORTEST
      Z1O(:,:)=0
#endif
      IF ( .NOT. lwith_barotropic_stokes_drift ) THEN
        !        CALL TRONEU
        CALL TRONEU2(IHALO_SOR,IMOD)
      ENDIF

#ifdef SORTEST
      ALLOCATE( z1(ie,je),z2(ie,je))
      z1(:,:)=z1o(:,:)
      z1o(:,:)=0.


      CALL TRONEU
      z2(:,:)=z1o(:,:)

      z1o(:,:)= z1(:,:)   !use z1o from troneu

      z1(:,:)=z1(:,:)-z2(:,:)
      WRITE(0,*)'z1',p_pe,MAXVAL(z1),MAXLOC(z1)

      CALL p_barrier

      DEALLOCATE(z1,z2)
      !        stop
#endif

    END IF
    IF ( ICONTRO .NE. 0 ) THEN
      CALL contro(39)
    END IF
    IF (time_verbosity > 2) time_diffs_subr(tag_subr_solver) = p_time() - tttr
    !
    !--------------------------------------------------------------------
    CALL OCVTRO
    !--------------------------------------------------------------------
    !
    IF (time_verbosity > 2) tttr = p_time()

    IF ( ICONTRO .NE. 0 ) THEN
      CALL contro(40)
    END IF
    CALL OCVTOT

    !      PRINT*,'nach ocvtot'
    IF ( ICONTRO .NE. 0 ) THEN
      CALL contro(41)
    END IF
    IF (time_verbosity > 2) time_diffs_subr(tag_subr_ocvtot) = p_time() - tttr

#ifdef OLD_IO
    IF(imean.NE.0)THEN
      CALL WRTE_MFL(LDAYS,LMONTS,LYEARS,IMEAN,NANF)
    END IF
#endif
#ifndef NO_NEW_IO
    CALL iolist_accumulate_staggered(ocean_iolist,ocean_varlist,model_time)
#endif

    IF (leddydiag) THEN
      CALL calc_eddydiag
    END IF
    !
    !UWE   RESET STABIO TO DRHODZ
    DO K=1,KE
      DO J=1,JE
        DO I=1,IE
          stabio(i, j, k) = (1000._wp / dz(k)) * stabio(i, j, k)
        END DO
      END DO
    END DO
    IF ( ICONTRO .NE. 0 ) THEN
      CALL contro(42)
    END IF


    CALL calc_icecutoff  ! this should go into mo_ocice

    CALL cell_thickness  ! update cell thickness

    CALL diagnosis
!    IF(itsdiag.GE.4) THEN
!    CALL diagnosis
!    ELSE
!      IF (MOD(ldtday,ndtday).EQ.ndtday/2) THEN
!        CALL diagnosis
!      END IF
!    END IF


    IF (time_verbosity > 2) tttr = p_time()
    IF ( ICONTRO .NE. 0 ) THEN
      CALL contro(43)
    END IF

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
        WRITE(0,*) 'uaccel nach ocuad',MAXVAL(uaccel),MAXLOC(uaccel)
        WRITE(0,*) 'vaccle nach ocvad',MAXVAL(vaccel),MAXLOC(vaccel)
        CALL contro(44)
      END IF

    ENDIF

    IF (time_verbosity > 2) time_diffs_subr(tag_subr_ocuvad) = p_time() - tttr

    IF(ibbl_transport.EQ.1)THEN               ! slope convection
      CALL SLOPETRANS
      IF ( ICONTRO .NE. 0 ) THEN
        CALL CONTRO(2)
      END IF
    END IF

    IF (time_verbosity > 2) tttr = p_time()

    IF (ladfs) CALL OCADFS
    IF (time_verbosity > 2) time_diffs_subr(tag_subr_ocadfs) = p_time() - tttr

    !--------------------------------------------------------------------
    !  ADVECTION | EDDY-DIFFUSION | DIFFUSION
    !--------------------------------------------------------------------
    IF (time_verbosity > 2) tttr = p_time()

#ifdef PBGC

    !  BGC tracer advection
    !
!!$! Dilution of bgc tracer due to fluxes of water.
!!$! BRINE contains ice volume change  (see SBR GROWTH).
!!$! PRECO contains atmospheric fluxes (see SBR GROWTH).
!!$! Tracer concentration is assumed to be zero in water added.
!!$
!!$        DO J=1,JE
!!$          DO I=1,IE
!!$            layer1_new(i,j)=0.
!!$            IF (ddpo(i,j,1).GT.0.)THEN
!!$              layer1_new(i,j)=DDPO(I,J,1)+ZO(I,J)                    &
!!$     &         -WO(I,J,1)*DT-SICTHO(I,J)*RHOICWA-SICSNO(I,J)*RHOSNWA
!!$            END IF
!!$          END DO
!!$        END DO
!!$
!!$      CALL DILUTE_BGC(IE,JE,layer1_bgc,layer1_new)


    !     IF(iocad .EQ. 3 .AND. ibbl_transport .EQ. 1)THEN
!!$       call maschk(ie,je,ke,44)
!!$ !$OMP PARALLEL
!!$ !$OMP DO
!!$         DO l=1,nocetra
!!$            DO K=1,KE
!!$               DO J=1,JE
!!$                  DO I=1,IE
!!$                     OCETRA(I,J,K,L)=MAX(0.,OCETRA(I,J,K,L))
!!$                  END DO
!!$               END DO
!!$            END DO
!!$ !$OMP END DO
!!$ !$OMP END PARALLEL
!!$         END DO
    !     END IF !  ADPO + SLOPECON_ADPO
#endif /*PBGC*/

    CALL bounds_exch(1,'p',rhoo,'mpiom 29')

    !  Preparation for advection
    IF (ladpo) CALL ocadpo_base

#ifdef PBGC
    !  Loop over tracers (S, T and advected BGC tracers)
    ntracerloop=ntraad+2
#else
    !  Loop over tracers (S and T)
    ntracerloop=2
#endif /*PBGC*/

#ifdef TRACER_OMP

    !  OpenMP routine-level parallel (outside routines) tracer calculations
    !   - advection, eddy-diffusion, and diffusion are calculated completely
    !     parallel for each tracer
    !   - use a maximum of ntracerloop OpenMP threads
    !   - without   HAMOCC: maximum of 2 OpenMP threads (for S and T)
    !   - including HAMOCC: maximum of ntraad+2 threads (COSMOS: 16)

!$OMP PARALLEL PRIVATE(L,LL)
!$OMP DO
#endif
    DO LL=1,ntracerloop
      L=LL-2
      IF (LL==1) THEN

        !  Advection of Temperature
        IF (ladpo) CALL ocadpo_trf(tho)
        CALL bounds_exch(1,'p',tho,'mpiom 30')

      ELSE IF (LL==2) THEN

        !  Advection of Salt
        IF (ladpo) CALL ocadpo_trf(sao)
        CALL bounds_exch(1,'p',sao,'mpiom 31')

#ifdef PBGC
      ELSE IF (LL>2) THEN

        ! Set tracer values at dry cells to bottom cell value as done
        ! for temperature and salinity in SBR OCTHER.

        DO K=2,KE
          DO J=1,JE
            DO I=1,IE
              ! FIXME: replace with lweto
              IF (weto(i, j, k) .LT. 0.5_wp) THEN
                OCETRA(I,J,K,L)=OCETRA(I,J,K-1,L)
              END IF
            END DO
          END DO
        END DO

        !  BGC tracer advection
        IF (ladpo) CALL ocadpo_trf(ocetra(1,1,1,L))

        ! Reset tracer values at dry cells to undefined.

#ifndef TRACER_OMP
!$OMP PARALLEL
!$OMP DO
#endif
        DO K=1,KE
          DO J=1,JE
            DO I=1,IE
              IF (weto(i, j, k) .LT. 0.5_wp) THEN
                OCETRA(I,J,K,L)=RMASKO
                !                     ELSE
                !                        OCETRA(I,J,K,L)=MAX(0.,OCETRA(I,J,K,L))
              END IF
            END DO
          END DO
        END DO
#ifndef TRACER_OMP
!$OMP END DO
!$OMP END PARALLEL
#endif

#endif /*PBGC*/

      END IF !  LL
    END DO ! LL
#ifdef TRACER_OMP
!$OMP END DO
!$OMP END PARALLEL
#endif


!    !--------------------------------------------------------------------
!    !  UPDATE_ZO : UPDATE OF SEALEVEL
    CALL update_zo

    !  Preparation for eddy-diffusion
    IF ( IBOLK .NE. 0 ) CALL ocjitr_base


    IF ( ICONTRO .NE. 0 ) CALL CONTRO(3)

#ifdef TRACER_OMP
!$OMP PARALLEL PRIVATE(L,LL)
!$OMP DO
#endif
    DO LL=1,ntracerloop
      L=LL-2
      IF (LL==1) THEN

        !  GM eddy-diffusion of Temperature
        IF ( IBOLK .NE. 0 ) CALL ocjitr_trf(tho)

      ELSE IF (LL==2) THEN

        !  GM eddy-diffusion of Salt
        IF ( IBOLK .NE. 0 ) CALL ocjitr_trf(sao)
        !           call maschk(ie,je,ke,47)

#ifdef PBGC
      ELSE IF (LL>2) THEN

        !  GM eddy-diffusion of BGC-tracer
        IF ( IBOLK .NE. 0 ) CALL ocjitr_trf(ocetra(1,1,1,L))

#endif /*PBGC*/

      END IF !  LL
    END DO ! LL
#ifdef TRACER_OMP
!$OMP END DO
!$OMP END PARALLEL
#endif


    IF ( ICONTRO .NE. 0 ) CALL CONTRO(133)

    !  Preparation for diffusion
    CALL octdiff_base

    IF (LCONVDIAG) CALL calc_potential_energy_release(1)

#ifdef TRACER_OMP
!$OMP PARALLEL PRIVATE(L,LL)
!$OMP DO
#endif
    DO LL=1,ntracerloop
      L = LL - 2
      IF (LL == 1) THEN

        !  Diffusion of Temperature
        CALL octdiff_trf(tho)

      ELSE IF (LL == 2) THEN

        !  Diffusion of Salt
        CALL octdiff_trf(sao)

#ifdef PBGC
      ELSE IF (LL > 2) THEN

        !  BGC tracer diffusion
        CALL octdiff_trf(ocetra(1, 1, 1, L))
#endif /*PBGC*/

      END IF !  LL
    END DO ! LL
#ifdef TRACER_OMP
!$OMP END DO
!$OMP END PARALLEL
#endif

    IF (LCONVDIAG) CALL calc_potential_energy_release(3)


    IF ( ICONTRO .NE. 0 ) CALL CONTRO(134)
    !           call maschk(ie,je,ke,48)
    !
    CALL RELAX_TS

    !      call maschk(ie,je,ke,51)
    !

#ifdef NUDGE_ECCO

    CALL READ_ECCO

    CALL NUDGE_T(tho, rreltem)

    CALL NUDGE_S(sao, rrelsal)
#endif

    IF ( ICONTRO .NE. 0 ) THEN
      CALL CONTRO(4)
    END IF

    CALL OCTIMF
    IF ( ICONTRO .NE. 0 ) THEN
      CALL CONTRO(5)
    END IF

    IF (time_verbosity > 2) time_diffs_subr(tag_region_adv) = p_time() - tttr
    !--------------------------------------------------------------------
    !  MOMENTUM DIFFUSION
    !
    IF (time_verbosity > 2) tttr = p_time()
    IF (AUS .GT. 1.e-12_wp) CALL OCSCHEP
    IF (time_verbosity > 2) time_diffs_subr(tag_subr_ocschep) = p_time() - tttr
    !--------------------------------------------------------------------
    !  BIHARMONIC MOMENTUM DIFFUSION, BOTTOM FRICTION,
    !  SHEAR DEPENDENT DIFFUSION

    IF (time_verbosity > 2) tttr = p_time()

    CALL OCVISC
    IF (time_verbosity > 2) time_diffs_subr(tag_subr_ocvisc) = p_time() - tttr
    IF ( ICONTRO .NE. 0 ) THEN
      CALL CONTRO(6)
    END IF

    CALL OCTIMF

#ifdef FLUXCORRECT
    !JJ     MONTHLY AVERAGE AND YEARLY AVERAGE FOR FLUX CORRECTION
    NANF = LDTDAY + ((LDAYS - 1)*NDTDAY)
    NEND = NDTDAY*MONLEN(LMONTS,LYEARS)
    ! CALCULATE END OF YEAR FOR 12 MONTH
    MYEND = monlen_sum(1,12,lyears)
    !:: CALCULATE ACTUAL DAY IN YEAR
    MYACT = monlen_sum(1,lmonts-1,lyears) + ldays -1
    MYACT = MYACT*NDTDAY+LDTDAY
    MYEND = MYEND*NDTDAY

    CALL FLUXC_MMEAN2D(LDAYS,LMONTS,LYEARS,                         &
         NANF,NEND,MYACT,MYEND)
#endif /*FLUXCORRECT*/

#ifdef OLD_IO
    IF (iMEAN .NE. 0) THEN
      CALL WRTE_MEAN(LDAYS,LMONTS,LYEARS,IMEAN,NANF)
    END IF
#endif
#ifndef NO_NEW_IO
    CALL iolist_accumulate(ocean_iolist,ocean_varlist,model_time)
    CALL iolist_poll_write_avglist(ocean_iolist, ocean_varlist, model_time)
#endif

#ifdef PBGC

    CALL rate_tracer(ddpo)

#ifdef OLD_IO
    CALL WRTE_MEANBGC
#endif
#ifndef NO_NEW_IO
    CALL iolist_accumulate(bgc_iolist, bgc_varlist,model_time)
    CALL iolist_poll_write_avglist(bgc_iolist, bgc_varlist, model_time)
#endif

#endif/*def PBGC */

    IF (ltidal .and. ltidal_diag) THEN
      CALL p_tide_timeseries
    ENDIF

#if defined (__coupled)

    ! Put data to coupler

    CALL couple_put_o2a(ldtrun)
#endif /*(__coupled)*/
    !
#if defined (__coupled)

    ! Update coupled run calendar

    CALL couple_calendar(dt,io_stdout)
#endif

#ifndef MESSY
    IF (time_verbosity > 0) THEN
      ttt = p_time()
      time_diffs_subr(tag_timestep) = ttt - ttts
      IF (time_verbosity > 2) &
        CALL timing_diagnostics(time_diffs_subr(1:tag_timestep - 1), &
             timestep_time_agg(1:tag_timestep - 1), &
             region_tags(1:tag_timestep - 1), time_verbosity > 3, ldtrun)
      IF (time_verbosity > 0) &
           CALL timing_diagnostics(time_diffs_subr(tag_timestep:tag_timestep), &
           timestep_time_agg(tag_timestep:tag_timestep), &
           region_tags(tag_timestep:tag_timestep), time_verbosity > 1, ldtrun)
      IF (time_verbosity > 1 .AND. p_pe == p_io) &
           WRITE(0,123) 'Time for 1 timestep', ttt - ttts
    END IF
#endif
#ifdef _PROFILE
    CALL trace_stop ('dt_loop', 7)
#endif
  END SUBROUTINE per_timestep_computation


  SUBROUTINE write_restart_file
    !  WRITE RESTART FILE AT END OF MONTH
    IF (IAUFW .EQ. 1) THEN
#ifdef _PROFILE
      CALL trace_start ('restartwrite', 11)
#endif
#ifdef OLD_IO
      ! Old I/O is overridden by new I/O, but we chose to run both as a test.
      CALL AUFW(ldt)
#endif
#ifndef NO_NEW_IO
      ! May override old I/O (see above).
      WRITE(0,*)'call iolist_write_final (ocean)'
      CALL iolist_write_final(ocean_iolist, ocean_varlist, model_time)
#endif
#ifdef _PROFILE
      CALL trace_stop ('restartwrite', 11)
#endif
#ifdef PBGC
#ifdef _PROFILE
      CALL trace_start ('restartwritebgc', 13)
#endif
      ! Old I/O is overridden by new I/O, but we chose to run both as a test.
      CALL AUFW_BGC(ie,je,ke,ddpo,gila,giph,tiestu                    &
           ,lyears,lmonts,ldays,ldt)
#ifndef NO_NEW_IO
      ! May override old I/O (see above).
      WRITE(0,*)'call iolist_write_final (bgc)'
      CALL iolist_write_final(bgc_iolist, bgc_varlist, model_time)
#endif
#ifdef _PROFILE
      CALL trace_stop ('restartwritebgc', 13)
#endif
#endif /*PBGC*/
    ELSE
      WRITE(IO_STDOUT,*)'STOPPED WITHOUT WRITING RESTART FILE, IAUFW= ',IAUFW
    END IF
  END SUBROUTINE write_restart_file

  SUBROUTINE register_diagnostics

#ifndef NO_NEW_IO
    USE mo_varlist, ONLY : ocean_varlist, varlist_is_code_registered  &
         ,TsSecOffset,TsCodeperSec,TsRegOffset,TsCodeperReg

    USE mo_diagnosis
    USE mo_eddydiag

    INTEGER :: n, ncode

    ! check if moc diagnostic is registered
    IF (varlist_is_code_registered(ocean_varlist,100) .OR. &
         varlist_is_code_registered(ocean_varlist,101) .OR. &
         varlist_is_code_registered(ocean_varlist,102)) THEN

      lcalc_moc = .TRUE.
      WRITE(io_stdout,*)'lcalc_moc=',lcalc_moc
      CALL ini_moc

    END IF

    IF (varlist_is_code_registered(ocean_varlist,27)) THEN
      lcalc_psi = .TRUE.
      WRITE(io_stdout,*)'lcalc_psi=',lcalc_psi
    END IF

    IF (varlist_is_code_registered(ocean_varlist,11)) THEN
      lcalc_zo_sqr = .TRUE.
      WRITE(io_stdout,*)'lcalc_zo_sqr=',lcalc_zo_sqr
    END IF

    IF (varlist_is_code_registered(ocean_varlist,12)) THEN
      lcalc_sst = .TRUE.
      WRITE(io_stdout,*)'lcalc_sst=',lcalc_sst
    END IF

    IF (varlist_is_code_registered(ocean_varlist,14)) THEN
      lcalc_sst_sqr = .TRUE.
      WRITE(io_stdout,*)'lcalc_sst_sqr=',lcalc_sst_sqr
    END IF

    IF (varlist_is_code_registered(ocean_varlist,16)) THEN
      lcalc_sss = .TRUE.
      WRITE(io_stdout,*)'lcalc_sss=',lcalc_sss
    END IF

    IF (varlist_is_code_registered(ocean_varlist,17)) THEN
      lcalc_bottom_pressure = .TRUE.
      WRITE(io_stdout,*)'lcalc_bottom_pressure=',lcalc_bottom_pressure
    END IF

    IF (varlist_is_code_registered(ocean_varlist,18)) THEN
      lcalc_rhopoto = .TRUE.
      WRITE(io_stdout,*)'lcalc_rhopoto=',lcalc_rhopoto
    END IF

    IF (varlist_is_code_registered(ocean_varlist,21) .OR.           &
      varlist_is_code_registered(ocean_varlist,22) ) THEN

      lcalc_upward_mass_transport = .TRUE.

      IF (varlist_is_code_registered(ocean_varlist,22)) THEN

        lcalc_upward_mass_transport_sqr = .TRUE.
        WRITE(io_stdout,*)'lcalc_upward_mass_transport_sqr='    &
             ,lcalc_upward_mass_transport_sqr
      ENDIF

      WRITE(io_stdout,*)'lcalc_upward_mass_transport='                &
           ,lcalc_upward_mass_transport

    END IF

    IF (varlist_is_code_registered(ocean_varlist,37) .OR.           &
         varlist_is_code_registered(ocean_varlist,38) ) THEN
      lcalc_ice_velocities=.TRUE.
      WRITE(io_stdout,*)'lcalc_ice_velocities=',lcalc_ice_velocities
    END IF


    IF (varlist_is_code_registered(ocean_varlist,23) .OR.           &
      varlist_is_code_registered(ocean_varlist,24) ) THEN

      lcalc_ocean_mass_transport = .TRUE.
      WRITE(io_stdout,*)'lcalc_ocean_mass_transport=',lcalc_ocean_mass_transport

    END IF

    IF (varlist_is_code_registered(ocean_varlist,138) .OR.           &
         varlist_is_code_registered(ocean_varlist,139) ) THEN
      lcalc_ice_transport=.TRUE.
      WRITE(io_stdout,*)'lcalc_ice_transport=',lcalc_ice_transport
    END IF

    IF (varlist_is_code_registered(ocean_varlist,142) .OR.           &
         varlist_is_code_registered(ocean_varlist,143) ) THEN
      lcalc_sictr=.TRUE.
      WRITE(io_stdout,*)'lcalc_sictr=',lcalc_sictr
    END IF


    IF (varlist_is_code_registered(ocean_varlist,181) .OR.    &
         varlist_is_code_registered(ocean_varlist,182) .OR.   &
         varlist_is_code_registered(ocean_varlist,183)) THEN
      lcalc_mixed_layer_thickness = .TRUE.
      WRITE(io_stdout,*)'lcalc_mixed_layer_thickness=',lcalc_mixed_layer_thickness
      IF (varlist_is_code_registered(ocean_varlist,181)) THEN
        lcalc_mixed_layer_thickness_sqr = .TRUE.
        WRITE(io_stdout,*)'lcalc_mixed_layer_thickness_sqr='        &
             ,lcalc_mixed_layer_thickness_sqr
      END IF
    END IF

    IF (varlist_is_code_registered(ocean_varlist,226)              &
         .OR. varlist_is_code_registered(ocean_varlist,227)        &
         .OR. varlist_is_code_registered(ocean_varlist,228)        &
         .OR. varlist_is_code_registered(ocean_varlist,229) ) THEN
      lcalc_ice_stress=.TRUE.
      WRITE(io_stdout,*)'lcalc_ice_stress=',lcalc_ice_stress
    END IF

    IF (varlist_is_code_registered(ocean_varlist,513)              &
         .OR. varlist_is_code_registered(ocean_varlist,514)        &
         .OR. varlist_is_code_registered(ocean_varlist,515)        &
         .OR. varlist_is_code_registered(ocean_varlist,516) ) THEN

      lcalc_global_mean=.TRUE.
      WRITE(io_stdout,*)'lcalc_global_mean=',lcalc_global_mean
    ENDIF

    IF (varlist_is_code_registered(ocean_varlist,518)) THEN
      lcalc_gmsl_diag=.TRUE.
      WRITE(io_stdout,*)'lcalc_global_mean=',lcalc_gmsl_diag
    ENDIF

    ! check if eddy diagnostic is registered
    DO ncode=230,244
      IF (varlist_is_code_registered(ocean_varlist,ncode)) THEN
        leddydiag=.TRUE.
      ENDIF
      WRITE(io_stdout,*)'leddydiag=',leddydiag
    ENDDO

    ! check if sections are registered

    DO n=1,SIZE(section)

      ncode=TsSecOffset+(n-1)*TsCodeperSec

      IF (varlist_is_code_registered(ocean_varlist,ncode)) THEN
        section(n)%register_netheat=.TRUE.
      END IF

      IF (varlist_is_code_registered(ocean_varlist,ncode+1)) THEN
        section(n)%register_netsalt=.TRUE.
      END IF

      IF (varlist_is_code_registered(ocean_varlist,ncode+2)) THEN
        section(n)%register_netwater=.TRUE.
      END IF

      IF (varlist_is_code_registered(ocean_varlist,ncode+3)) THEN
        section(n)%register_sice=.TRUE.
      END IF

      IF ( varlist_is_code_registered(ocean_varlist,ncode+4) .AND. Section(n)%SecLay1(2) /= 0  ) THEN
        section(n)%register_layer1=.TRUE.
      ELSE
        WRITE(io_stdout,*) 'Warning: register code ',ncode+4,'layer1transport_'//section(n)%SecName,  &
             'between L=',Section(n)%SecLay1(1),' and L=',Section(n)%SecLay1(2)               &
             ,' Layer is not defined and code is therefor skipped'
      END IF

      IF ( varlist_is_code_registered(ocean_varlist,ncode+5) .AND. Section(n)%SecLay2(2) /= 0  ) THEN
        section(n)%register_layer2=.TRUE.
      ELSE
        WRITE(io_stdout,*) 'Warning: register code ',ncode+5,'layer2transport_'//section(n)%SecName,  &
             'between L=',Section(n)%SecLay2(1),' and L=',Section(n)%SecLay2(2)               &
             ,' Layer is not defined and code is therefor skipped'
      END IF

    ENDDO

    ! check regions

    DO n=1,SIZE(region)

      ncode=TsRegOffset+(n-1)*TsCodeperReg

      IF (varlist_is_code_registered(ocean_varlist,ncode).OR.varlist_is_code_registered(ocean_varlist,ncode+1)) THEN
        lcalc_region_seaice=.TRUE.
        lcalc_region=.TRUE.
      END IF

      IF (varlist_is_code_registered(ocean_varlist,ncode+2).OR.varlist_is_code_registered(ocean_varlist,ncode+3)) THEN
        lcalc_region_fluxes=.TRUE.
        lcalc_region=.TRUE.
      END IF

      IF (varlist_is_code_registered(ocean_varlist,ncode+4).OR.varlist_is_code_registered(ocean_varlist,ncode+5) .OR.     &
           varlist_is_code_registered(ocean_varlist,ncode+6).OR.varlist_is_code_registered(ocean_varlist,ncode+7) ) THEN
        lcalc_region_temperature=.TRUE.
        lcalc_region=.TRUE.
      END IF

      IF (varlist_is_code_registered(ocean_varlist,ncode+8).OR.varlist_is_code_registered(ocean_varlist,ncode+9) .OR.     &
           varlist_is_code_registered(ocean_varlist,ncode+10).OR.varlist_is_code_registered(ocean_varlist,ncode+11) ) THEN
        lcalc_region_salinity=.TRUE.
        lcalc_region=.TRUE.
      END IF

      IF (varlist_is_code_registered(ocean_varlist,ncode+12).OR.varlist_is_code_registered(ocean_varlist,ncode+13) .OR.     &
           varlist_is_code_registered(ocean_varlist,ncode+14).OR.varlist_is_code_registered(ocean_varlist,ncode+15) ) THEN
        lcalc_region_kinetic_energy=.TRUE.
        lcalc_region=.TRUE.
      END IF

    END DO
#endif

  END SUBROUTINE register_diagnostics

END PROGRAM MPIOM

#ifdef __IFC /* Intel Fortran Compiler */
SUBROUTINE handler(isig,icode)

  WRITE(io_stdout,*) 'Handler: ',isig,icode

END SUBROUTINE handler
#endif
