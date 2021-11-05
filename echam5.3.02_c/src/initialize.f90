SUBROUTINE initialize

  ! Description:
  !
  ! Set up constants in various modules.
  !
  ! Method:
  !
  ! This subroutine initializes all the variables and arrays
  ! in modules.
  !
  ! *initialize* is called from *control*
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, February 1982, original source
  ! R.G AND M.J, ECMWF, December 1982, changed
  ! U. Schlese, DKRZ, in 1994, and 1995, changed
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_io,           ONLY: IO_init
#ifndef MESSY
  USE mo_tracer,       ONLY: initrac
#endif
  USE mo_time_control, ONLY: lfirst_cycle, init_manager, init_events, init_times
#ifdef OBSOLETE
  USE mo_column,       ONLY: setcolumn
#endif
  USE mo_nudging_init, ONLY: NudgingInit, NDG_INI_IO, NDG_INI_STREAM
#ifdef OBSOLETE
  USE mo_greenhouse_gases, ONLY: init_ghg
  USE mo_radiation,    ONLY: ighg
  USE mo_control,      ONLY: lso4
  USE mo_so4,          ONLY: read_so4nat, read_so4all
  USE mo_column,       ONLY: inicolumn
#endif
  USE m_alloc_mods,    ONLY: alloc_mods ! module subroutine
#ifdef OBSOLETE
  USE mo_control,      ONLY: lcolumn, lvctch, nlev, nlevp1, nvclev, vct
#else
  USE mo_control,      ONLY: lvctch, nlev, nlevp1, nvclev, vct
#endif
  USE mo_echam_yaxt,   ONLY: yaxt_initialize, add_yaxt_gp_nlevs, yaxt_init_gp_coords
  USE mo_mpi,          ONLY: p_all_comm
  USE mo_control,      ONLY: lyaxt_transposition, nlev
  IMPLICIT NONE

  !  External subroutines 
! op_pj_20130407+
#ifndef MESSY
! op_pj_20130407-
  EXTERNAL inidoc, inictl, setdyn, setphys, setrad, setgws
! op_pj_20130407+
#else
#ifdef OBSOLETE
  EXTERNAL inidoc, inictl, setdyn, setphys !!$, setgws ! op_pj_20160617
#else
  EXTERNAL inidoc, inictl, setdyn !!$, setphys !!$, setgws ! op_pj_20160617
#endif
#endif
! op_pj_20130407-

  !  Executable statements 


  !-- 1. Set control variables

  !-- 1.1 Set general control variables and time stepping
  !--     Set i/o units and buffer indices

  CALL inictl

  !-- 1.2 Initialize netCDF IO

#ifdef MESSY
  CALL messy_setup(1)
#endif
  CALL IO_init
#ifdef MESSY
  ! mz_pj_20080418+
  ! Do not move to IO_init! This will create circular dependencies
  ! (at least in E5301):
  ! mo_geoloc -> mo_o3clim -> mo_io -> messy_main_channel_bi
  ! op_pj_20140328+
  !!$CALL messy_channel_init_restart
  CALL messy_setup(2) ! -> main_channel_setup -> channel_init_restart_bi
  ! op_pj_20140328-
  ! mz_pj_20080418-
#endif

  !-- 1.2 Preset constants in *mo_doctor*

#ifdef OBSOLETE
  !-- 1.3 Initialize column model

  CALL inicolumn (lcolumn, lvctch, nlev, nlevp1, nvclev, vct)
#endif

  CALL alloc_mods

  !-- 1.4 Preset constants in *mo_doctor*

  CALL inidoc

  IF (lfirst_cycle) CALL init_manager

  CALL NudgingInit(NDG_INI_IO)

  CALL init_times

  !-- 2. Compute decomposition 

! mz_ab_20100503+
#ifndef MESSY
  CALL init_decomposition 
#else
! op_pj_20140328+
  CALL messy_setup(3) ! -> main_decomp_setup
! op_pj_20140328-
#endif
! mz_ab_20100503-

#ifdef HAVE_YAXT
  IF (lyaxt_transposition) THEN
    CALL yaxt_initialize(p_all_comm)
    CALL yaxt_init_gp_coords
    CALL add_yaxt_gp_nlevs((/1, nlev, nlev+1/))
  ENDIF
#endif

  CALL NudgingInit(NDG_INI_STREAM)

#ifdef OBSOLETE
  ! initialize column model

  CALL setcolumn
#endif

  !-- 3. Preset, modify and derive values needed in the
  !      dynamics and the initialisation and call helmo the first time.

  CALL setdyn

  !-- 4. Preset, modify and derive values needed in the physics.

#ifdef OBSOLETE
  CALL setphys
#endif

  !-- 5. Preset, modify and derive values needed in the radiation.

! op_pj_20130406+
#ifndef MESSY
! op_pj_20130406-
  CALL setrad
! op_pj_20130406+
#endif
! op_pj_20130406-

#ifdef OBSOLETE
  IF(lso4) THEN
    CALL read_so4nat
    CALL read_so4all
  ENDIF
#endif

  !-- 6. Preset, modify and derive values needed in the gwspectrum param.

! mz_ab_20100829
#ifndef MESSY
  CALL setgws
#endif

  !-- 7. Prepare greenhouse gas scenario

#ifdef OBSOLETE
  IF(ighg .NE. 0) CALL init_ghg
#endif

#ifdef MESSY
! um_ak_20090604+
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  CALL messy_initialize
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! um_ak_20090604-
#endif

  !-- final event evaluation

  CALL init_events

  !--    Initialize submodels

#ifndef MESSY
  CALL call_init_submodels
#endif

  !-- Preset values for tracer transport

#ifndef MESSY
  CALL initrac
#else
  ! mz_pj_20040329+
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  CALL messy_new_tracer
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! mz_pj_20040329-
#endif

END SUBROUTINE initialize
