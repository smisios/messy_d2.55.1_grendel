#if defined(__uxp__) || defined(__SX__) || defined(ES)
#define FAST_AND_DIRTY 1
#endif

SUBROUTINE control

  ! Description:
  !
  ! Control routine for the model.
  !
  ! Method:
  !
  ! This subroutine controls the running of the model.
  ! 
  ! *control* is called from the main program (*master*).
  !
  ! Externals:
  !   *initialize*    called to initialize modules.
  !   *stepon*        controls time integration
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, February 1982, original source
  ! R.G and M.J, ECMWF, December 1982, changed
  ! U. Schlese, MPI, August 1989, new structure
  ! U. Schlese, DKRZ, September 1994, interface *drive* added
  ! U. Schlese, DKRZ, January 1995, reading of optional files added
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! U. Schlese, DKRZ and M. Esch, MPI, July 1999, modifications for ECHAM5
  ! L. Kornblueh, MPI, June 1999, parallel version (MPI based)
  ! U. Schlese, DKRZ, November 1999, "drive" removed, "stepon" called directly
  ! U. Schlese, DKRZ, December 1999, modifications for coupling
  ! M.A. Giorgetta, MPI, May 2000, ozone initialization removed, see setrad
  ! S. Legutke, MPI M&D, July 00, modifications for coupling interface
  ! I. Kirchner, MPI, December 2000, date/time control/nudging
  ! U. Schulzweida, MPI, May 2002, blocking (nproma)
  ! I. Kirchner, MPI, Aug 2002, tendency diagnostics revision
  ! A. Rhodin, DWD, June 2002, call subroutines cleanup_...
  ! M. Esch, MPI, September 2002, modifications for mixed layer ocean
  ! L. Kornblueh, MPI, October 2003, added setup for AMIP2 global diagnostics
  ! 
  ! for more details see file AUTHORS
  ! 
#ifdef OBSOLETE
  USE mo_sst,             ONLY: readsst, readice, readflux
#else
  USE mo_sst,             ONLY: readsst, readice
#endif
#ifdef OBSOLETE
  USE mo_control,         ONLY: lcouple, lnudge, ltdiag, lhd, lmidatm, &
                                lnmi, lmlo, numfl1, numfl2, nlon, ngl, &
                                ldiagamip, lso4
#else
  USE mo_control,         ONLY: lcouple, lnudge, ltdiag, lmidatm, &
                                lnmi, numfl1, numfl2, nlon, ngl
#endif
  USE mo_machine,         ONLY: machine_setup
  USE mo_legendre,        ONLY: inileg, cleanup_legendre
#ifndef MESSY
  ! mz_rs_20040329+
  USE mo_tracer,          ONLY: reademi, ntrac, cleanstatr
  ! mz_rs_20040329-
#endif
  USE mo_nmi,             ONLY: NMI_Init, NMI_Close
#ifdef OBSOLETE
  USE mo_time_control,    ONLY: lstart, lfirst_cycle, lresume,         &
                                construct_events, lstop, l_diagrad
#else
  USE mo_time_control,    ONLY: lstart, lfirst_cycle, lresume,         &
                                construct_events, lstop
#endif
  USE mo_clim,            ONLY: readtslclim, readvltclim, readvgratclim
  USE mo_nudging_init,    ONLY: NudgingInit, NDG_CLEAN_MEM, NDG_CLOSE
  USE mo_diag_tendency,   ONLY: DIAG_Init, IDIAG_INIT, IDIAG_FREE
#ifdef OBSOLETE
  USE mo_diag_amip2,      ONLY: init_amip2_diag
  USE mo_couple,          ONLY: couple_init
  USE mo_column,          ONLY: resetcolumn
#endif
  USE mo_grib,            ONLY: init_grib, cleanup_grib
  USE mo_memory_streams,  ONLY: init_memory, free_memory
  USE mo_geoloc,          ONLY: init_geoloc, cleanup_geoloc
  USE mo_advection,       ONLY: iadvec, tpcore, semi_lagrangian,       &
                                spitfire
  USE mo_semi_lagrangian, ONLY: init_semi_lagrangian,                  &
                                cleanup_semi_lagrangian
  USE mo_spitfire,        ONLY: init_spitfire
  USE mo_tpcore,          ONLY: init_tpcore, cleanup_tpcore
  USE mo_timer,           ONLY: init_timer, cleanup_timer
#ifdef OBSOLETE
  USE mo_aero_tanre,      ONLY: init_aero, cleanup_aero
  USE mo_radiation,       ONLY: io3, iaero
#endif
  USE mo_decomposition,   ONLY: ldc=>local_decomposition,              &
                                cleanup_decomposition
  USE mo_exception,       ONLY: finish
#ifdef OBSOLETE
  USE mo_hydrology,       ONLY: init_hydrology, cleanup_hydrology
#endif
  USE mo_scan_buffer,     ONLY: cleanup_scanbuffer
#ifdef OBSOLETE
  USE mo_o3clim,          ONLY: cleanup_o3clim
#endif
  USE mo_clim,            ONLY: cleanup_clim
  USE mo_sst,             ONLY: cleanup_sst
  USE mo_tmp_buffer,      ONLY: cleanup_tmp_buffer
  USE mo_gaussgrid,       ONLY: cleanup_gaussgrid
#ifdef FFT991
  USE mo_fft991,          ONLY: cleanup_fft991
#else
  USE mo_fft992,          ONLY: cleanup_fft992
#endif
  USE m_alloc_mods,       ONLY: dealloc_mods
  USE mo_netcdf,          ONLY: cleanup_netcdf
  USE mo_io,              ONLY: cleanup_io
  USE mo_doctor,          ONLY: nerr,nout
  USE mo_field,           ONLY: field1, field2
  USE mo_mpi,             ONLY: p_parallel_io, p_parallel
#ifdef HAVE_YAXT
  USE mo_mpi,             ONLY: p_all_comm
#endif
  USE mo_memory_g3b,      ONLY: alake, slf
#ifdef OBSOLETE
  USE mo_so4,             ONLY: cleanup_so4
#endif

#ifdef FAST_AND_DIRTY
  USE mo_memory_f,        ONLY: resort_restart_read_memory_f
#endif
  USE mo_control,      ONLY: lyaxt_transposition
  USE mo_echam_yaxt,   ONLY: setup_yaxt_decomposition, generate_yaxt_redist

  IMPLICIT NONE

!  External subroutines 
#ifdef OBSOLETE
  EXTERNAL stepon, initialize, inipost,                                &
           labrun, readfld, inhysi, ioinitial, scan2
#else
  EXTERNAL stepon, initialize,                                         &
           labrun, readfld, inhysi, ioinitial, scan2
#endif

!  Executable statements 

!-- Default advection scheme

  iadvec = tpcore

!--  Print machine specific values

  CALL machine_setup

  CALL construct_events

  CALL init_timer

!--  Initialize modules and parallel decomposition

  CALL initialize

!-- Initialize time independent surface parameters 

  CALL init_geoloc

#ifdef OBSOLETE
  IF ( iaero == 2 .OR. iaero == 3 .OR. iaero == 4 ) THEN
    CALL init_aero
  END IF
#endif

  IF ( .NOT. ldc% lreg ) THEN
    IF ( iadvec .EQ. semi_lagrangian .OR.  &
#ifdef OBSOLETE
         iadvec .EQ. spitfire        .OR.  &
         l_diagrad                  ) THEN !.OR.  & ! mz_pj_20050615
#else
         iadvec .EQ. spitfire       ) THEN
#endif
!         lnudge                            ) THEN  ! mz_pj_20050615
      WRITE(nerr,*) "NPROMA does not work with:"
      WRITE(nerr,*) "  - semi_lagrangian advection"
      WRITE(nerr,*) "  - spitfire advection"
#ifdef OBSOLETE
      WRITE(nerr,*) "  - l_diagrad"
#endif
!      WRITE(nerr,*) "  - nudging" ! mz_pj_20050615
      CALL finish('control','NPROMA not implemented')
    END IF
  END IF

!-- Initialize tendency diagnostics
 
  IF (ltdiag) CALL DIAG_Init(IDIAG_INIT)

!--  Initialize memory

  CALL init_memory

#ifdef HAVE_YAXT
  IF (lyaxt_transposition) THEN
    CALL setup_yaxt_decomposition
    CALL generate_yaxt_redist(p_all_comm)
  ENDIF
#endif
  !-- Preset values needed in the advection scheme

  SELECT CASE (iadvec)
  CASE (semi_lagrangian)
    CALL init_semi_lagrangian
  CASE (spitfire)
    CALL init_spitfire
  CASE (tpcore)
    CALL init_tpcore
  END SELECT

!
!--  Initialize HD Model
!
#ifndef MESSY
!mz_ap_20070920+
  IF (lhd) CALL init_hydrology
!mz_ap_20070920-
#endif

!--  Compute *Legendre polynomials and parameters needed
!    for the *Legendre transforms.

  CALL inileg

#ifdef MESSY
! mz_pj_20040315+
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CALL messy_init_coupling
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! mz_pj_20040315-
#endif

  IF (lresume) THEN
     CALL iorestart
  ELSE IF (lstart) THEN
     CALL ioinitial
  END IF

  CALL inhysi

#ifdef FAST_AND_DIRTY  
  CALL resort_restart_read_memory_f
#endif
  IF (lstart) CALL scan2

  IF (lnmi .AND. lfirst_cycle) CALL NMI_Init(lnudge)


!-- Initialize postprocessing
#ifdef OBSOLETE
  CALL inipost

  IF (ldiagamip) CALL init_amip2_diag
#endif

  CALL init_grib

!-- Read optional sst-file if not in coupled mode

  IF (.NOT. lcouple) CALL readsst
  
! --  read optional sea ice if not in coupled mode

  IF (.NOT.lcouple) CALL readice
#ifdef __cpl_fluxes4
  IF (lcouple ) CALL readice
#endif

#ifdef OBSOLETE
! --  read optional flux correction if in mixed layer mode

  IF (lmlo) CALL readflux
#endif

!-- Read  climate land surface temperatures

  CALL readtslclim

!-- Read climate leaf-area index and climate vegetation ratio

  CALL readvltclim

  CALL readvgratclim

!-- Read optional fields

#ifndef MESSY
  ! mz_rs_20040329+
  CALL reademi
  ! mz_rs_20040329-
#endif

  CALL readfld

!
!---  control a coupled run
!        -initialise  data exchange with coupler
!        -check whether Echam control parameters are consistent  
!                                          with coupling system 

#ifndef MESSY 
!mz_ap_20070830+
! only used for OASIS coupler, not for MESSy
  IF (lcouple .AND. lfirst_cycle) THEN 
    CALL couple_init(nout,nerr)
  END IF
!mz_ap_20070830-
#endif

  CALL labrun

!-- Start time integration

  CALL stepon 

!-- Clean up

#ifdef OBSOLETE
  CALL resetcolumn
#endif
  IF (iadvec == tpcore) CALL cleanup_tpcore
  IF (iadvec == semi_lagrangian) CALL cleanup_semi_lagrangian

  CALL cleanup_timer

#ifndef MESSY
!mz_ap_20070920-
  IF(lhd) CALL cleanup_hydrology
!mz_ap_20070920-
#endif

  IF (numfl2 > 0) DEALLOCATE (field2)
  IF (numfl1 > 0) DEALLOCATE (field1)

  CALL free_memory
  CALL cleanup_scanbuffer
#ifdef OBSOLETE
  CALL cleanup_o3clim
#endif
  CALL cleanup_clim
  CALL cleanup_sst
  CALL cleanup_grib
  CALL cleanup_legendre
#ifdef OBSOLETE
  CALL cleanup_aero
#endif
  CALL cleanup_tmp_buffer
  CALL cleanup_gaussgrid
#ifdef FFT991
  CALL cleanup_fft991
#else
  CALL cleanup_fft992
#endif
  CALL cleanup_decomposition
  CALL dealloc_mods
  CALL cleanup_netcdf
  CALL cleanup_io
  CALL cleanup_geoloc
#ifndef MESSY
  ! mz_rs_20040329+
  CALL cleanstatr
  ! mz_rs_20040329-
#endif
#ifdef OBSOLETE
  IF(lso4) CALL cleanup_so4
#endif
  IF (lnudge) THEN
     CALL NudgingInit(NDG_CLEAN_MEM)
     IF (lstop) CALL NudgingInit(NDG_CLOSE)
  END IF
  
  IF (lnmi .AND. lstop) CALL NMI_Close

  IF (ltdiag)  CALL DIAG_Init(IDIAG_FREE)
 
END SUBROUTINE control
