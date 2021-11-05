!==============================================================================
!
! This file contains the calls to the main entry points of the MESSy submodels
! Authors:
!  Astrid Kerkweg,  UNI-MZ,        2010
!  Frank  Kalinka,  UNI-Frankfurt, 2010
!  Patrick Joeckel, DLR,           2010
!
!==============================================================================

SUBROUTINE messy_setup

  ! infrastructure
  USE messy_main_switch      !ONLY: USE_*
  USE messy_main_switch_bi,  ONLY: main_switch_setup
  USE messy_main_blather_bi, ONLY: main_blather_setup
  USE messy_main_qtimer_bi,  ONLY: main_qtimer_setup
  USE messy_main_timer_bi,   ONLY: main_timer_setup
!!$  USE messy_main_channel_bi, ONLY: messy_channel_init_restart
  USE messy_main_channel_bi, ONLY: main_channel_setup

  IMPLICIT NONE

  ! LOCAL
  INTEGER :: nclock = 1 ! ub_ak_20171109
  ! infrastructure
  CALL main_blather_setup
  ! TIMER SPLIT NEEDED HERE, TO ENABLE OTHER MODELS TO OVERWRITE NAMELIST
  ! ENTRY OF START DATE AND STOP DATE (e.g. MMD_CLIENT)
  ! ub_ak_20171109+
  !CALL main_timer_setup(1)
  CALL main_timer_setup(1, nclock)
  ! ub_ak_20171109-
  CALL main_channel_setup
  ! ub_ak_20171109+
  !CALL main_timer_setup(2) ! op_pj_20120307
  CALL main_timer_setup(2, nclock) ! op_pj_20120307
  ! ub_ak_20171109-
  ! do NOT change the following order of main_* subroutines
  CALL main_qtimer_setup
  CALL main_switch_setup

END SUBROUTINE messy_setup

!==============================================================================

SUBROUTINE messy_initialize

  USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io
  ! infrastructure
  USE messy_main_switch   ! ONLY: USE_*
  USE messy_main_channel_bi,ONLY: main_channel_initialize
  USE messy_main_timer_bi,  ONLY: main_timer_initialize
  USE messy_main_data_bi,   ONLY: main_data_initialize
  USE messy_main_grid_bi,   ONLY: main_grid_initialize
  USE messy_main_tracer_bi, ONLY: main_tracer_initialize
  USE messy_main_tools_bi,  ONLY: main_tools_initialize
  USE messy_main_import_bi, ONLY: main_import_initialize
  !submodels

  USE messy_aeropt_si,     ONLY:   aeropt_initialize
  USE messy_cloudopt_si,   ONLY: cloudopt_initialize
  USE messy_orbit_si,      ONLY:    orbit_initialize
  USE messy_rad_si,        ONLY:      rad_initialize, rad_mbm_initialize

  IMPLICIT NONE

  ! -> read in the MBM_RAD CPL_MBM namelist, the input
  !    is required for main_grid_initialize and main_data_initialize
  CALL rad_mbm_initialize
  ! infrastructure
!!! do NOT change the following order of main_* subroutines
  CALL main_timer_initialize
  CALL main_grid_initialize
  CALL main_data_initialize
  CALL main_channel_initialize
  CALL main_tracer_initialize
  CALL main_tools_initialize
  CALL main_import_initialize

  ! submodels
  IF (USE_AEROPT)   CALL   aeropt_initialize
  IF (USE_CLOUDOPT) CALL cloudopt_initialize
  IF (USE_ORBIT)    CALL    orbit_initialize
  IF (USE_RAD)      CALL      rad_initialize

END SUBROUTINE messy_initialize

!==============================================================================

 SUBROUTINE messy_read_restart

   USE messy_main_timer,        ONLY: lresume
   USE messy_main_channel_bi,   ONLY: main_channel_read_restart

   IF (.NOT. lresume) RETURN

   CALL main_channel_read_restart
   CALL main_channel_read_restart

 END SUBROUTINE messy_read_restart

!==============================================================================

 SUBROUTINE messy_new_tracer

   ! infrastructure
   USE messy_main_switch      !ONLY: USE_*
   USE messy_main_tracer_bi,  ONLY: main_tracer_new_tracer

   ! submodels
   ! USE messy_XXX_si,     ONLY: XXX_init_tracer

   IMPLICIT NONE

   CALL main_tracer_new_tracer(1)  ! define tracer set

   ! submodels
   !IF (USE_XXX) CALL XXX_new_tracer

   CALL main_tracer_new_tracer(2) ! define tracer families
   CALL main_tracer_new_tracer(3) ! diagnostic output

 END SUBROUTINE messy_new_tracer

!==============================================================================

 SUBROUTINE messy_init_memory

   ! infrastructure
   USE messy_main_switch      !ONLY: USE_*
   USE messy_main_qtimer_bi,  ONLY: main_qtimer_init_memory
   USE messy_main_data_bi,    ONLY: main_data_init_memory
   USE messy_main_grid_bi,    ONLY: main_grid_init_memory
   USE messy_main_tracer_bi,  ONLY: main_tracer_init_memory
   USE messy_main_channel_bi, ONLY: main_channel_init_memory
   USE messy_main_import_bi,  ONLY: main_import_init_memory

   ! submodels
   USE messy_aeropt_si,       ONLY:   aeropt_init_memory
   USE messy_cloudopt_si,     ONLY: cloudopt_init_memory
   USE messy_orbit_si,        ONLY:    orbit_init_memory
   USE messy_rad_si,          ONLY:      rad_init_memory

   IMPLICIT NONE

   ! infrastructure
   CALL main_qtimer_init_memory
   CALL main_channel_init_memory   ! setup tracer memory
   CALL main_tracer_init_memory(1) ! create data channel BML <-> SMIL
   CALL main_grid_init_memory(1)
   CALL main_data_init_memory
   CALL main_grid_init_memory(2)
   CALL main_import_init_memory

   ! submodels
   IF (USE_AEROPT)   CALL     aeropt_init_memory
   IF (USE_CLOUDOPT) CALL   cloudopt_init_memory
   IF (USE_ORBIT)    CALL      orbit_init_memory
   IF (USE_RAD)      CALL        rad_init_memory

   ! - associate tracer memory to MESSy channel(s)
   ! - setting meta information of family-members to fraction
   !   (for advection initialization)
   CALL main_tracer_init_memory(2)
 END SUBROUTINE messy_init_memory

 !==============================================================================

 SUBROUTINE messy_init_coupling

   ! infrastructure
   USE messy_main_switch      !ONLY: USE_*
   USE messy_main_channel_bi, ONLY: main_channel_init_coupling
   USE messy_main_tracer_bi,  ONLY: main_tracer_init_coupling
   USE messy_main_rnd_bi,     ONLY: main_rnd_init_coupling
   USE messy_main_data_bi,    ONLY: main_data_init_coupling

   ! submodels
   USE messy_aeropt_si,       ONLY:   aeropt_init_coupling
   USE messy_cloudopt_si,     ONLY: cloudopt_init_coupling
   USE messy_orbit_si,        ONLY:    orbit_init_coupling
   USE messy_rad_si,          ONLY:      rad_init_coupling

   IMPLICIT NONE

   ! infrastructure
   ! - resetting meta information of family-members (to tracers)
   !   (after advection initialisation)
   CALL main_tracer_init_coupling(1)
   ! - modify attributes according to ADD_ATT in CTRL namelist
   CALL main_channel_init_coupling(1)

   CALL main_data_init_coupling

   ! submodels
   IF (USE_AEROPT)   CALL   aeropt_init_coupling
   IF (USE_CLOUDOPT) CALL cloudopt_init_coupling
   IF (USE_ORBIT)    CALL    orbit_init_coupling ! new channel objects
   IF (USE_RAD)      CALL      rad_init_coupling

   CALL main_rnd_init_coupling     ! create channel objects for state vectors
   CALL main_channel_init_coupling(2)

 END SUBROUTINE messy_init_coupling

 !==============================================================================

 SUBROUTINE messy_init_tracer

   ! infrastructure
   USE messy_main_switch   ! ONLY: USE_*
   USE messy_main_tracer_bi, ONLY: main_tracer_init_tracer
   USE messy_main_rnd_bi,    ONLY: main_rnd_init_tracer
   USE messy_main_import_bi, ONLY: main_import_init_tracer

   ! submodels
   ! USE messy_XXX_si,     ONLY: XXX_init_tracer

   IMPLICIT NONE

   ! infrastructure
   CALL main_rnd_init_tracer       ! reset state vectors after restart
   CALL main_tracer_init_tracer(1) ! check tracer init from rerun
   CALL main_import_init_tracer

   CALL main_tracer_init_tracer(4) ! initialize via tracer_init (tracer.nml)

   ! submodels
   !IF (USE_XXX)  CALL XXX_init_tracer

   ! infrastructure
   CALL main_tracer_init_tracer(2) ! check tracer init by init_tracer
                                   ! set to constant if not yet initialized
   CALL main_tracer_init_tracer(3) ! diagnose tracer initialization

 END SUBROUTINE messy_init_tracer

 !==============================================================================

 SUBROUTINE messy_time(flag)

   USE messy_main_timer_bi, ONLY: main_timer_time

   INTEGER, INTENT(IN) :: flag

   CALL main_timer_time(flag)     ! 1: set, 2: reset

 END SUBROUTINE messy_time

 !==============================================================================

 SUBROUTINE messy_global_start

   ! infrastructure
   USE messy_main_switch      !ONLY: USE_*
   USE messy_main_data_bi,    ONLY: main_data_global_start
   USE messy_main_tracer_bi,  ONLY: main_tracer_global_start
   USE messy_main_import_bi,  ONLY: main_import_global_start
   USE messy_main_channel_bi, ONLY: main_channel_global_start

   ! submodels
   USE messy_orbit_si,        ONLY: orbit_global_start
   USE messy_rad_si,          ONLY:   rad_global_start

   IMPLICIT NONE

   ! infrastructure
   CALL main_tracer_global_start
   CALL main_channel_global_start
   CALL main_import_global_start
   CALL main_data_global_start

   ! submodels
   IF (USE_RAD)      CALL   rad_global_start(1)
   IF (USE_ORBIT)    CALL orbit_global_start
   IF (USE_RAD)      CALL   rad_global_start(2)

 END SUBROUTINE messy_global_start

!============================================================================

!============================================================================

 SUBROUTINE messy_local_start

   ! infrastructure
   USE messy_main_switch     !ONLY: USE_*
   USE messy_main_data_bi,   ONLY: main_data_local_start
   USE messy_main_tracer_bi, ONLY: main_tracer_local_start

   ! submodels
   ! USE messy_XXX_si,     ONLY: XXX_local_start

   IMPLICIT NONE

   CALL main_data_local_start
   CALL main_tracer_local_start

   ! submodels
   !IF (USE_XXX) CALL XXX_local_start

 END SUBROUTINE messy_local_start

!==============================================================================

!==============================================================================

SUBROUTINE messy_radiation

  USE messy_main_switch   ! ONLY: USE_*

  USE messy_aeropt_si,     ONLY:   aeropt_radiation
  USE messy_cloudopt_si,   ONLY: cloudopt_radiation
  USE messy_rad_si,        ONLY:      rad_radiation

  IMPLICIT NONE

  ! special requirements:
  ! - aeropt should be called directly before radiation
  !   such that no other process can influence the aerosol
  !   distribution after the optical properties have been determined
  IF (USE_AEROPT)   CALL   aeropt_radiation
  IF (USE_CLOUDOPT) CALL cloudopt_radiation
  IF (USE_RAD)      CALL      rad_radiation

END SUBROUTINE messy_radiation

!==============================================================================

!==============================================================================

SUBROUTINE messy_radheat

  USE messy_main_switch   ! ONLY: USE_*

  USE messy_rad_si,        ONLY:     rad_radheat

  IMPLICIT NONE

  ! special requirements:
  IF (USE_RAD)     CALL     rad_radheat

END SUBROUTINE messy_radheat

!==============================================================================


!==============================================================================
 SUBROUTINE messy_local_end

   ! infrastructure
   USE messy_main_switch     !ONLY: USE_*

   IMPLICIT NONE

   ! submodels
   ! IF (USE_XXX) CALL XXX_local_end

 END SUBROUTINE messy_local_end

!==============================================================================

 SUBROUTINE messy_global_end

   ! infrastructure
   USE messy_main_switch      !ONLY: USE_*
   USE messy_main_qtimer_bi,  ONLY: main_qtimer_global_end
   USE messy_main_tracer_bi,  ONLY: main_tracer_global_end
!!$   USE messy_main_channel_bi, ONLY: main_channel_global_end

   ! submodels
   !USE messy_XXX_si,      ONLY: submod1_XXX_end

   IMPLICIT NONE

   ! infrastructure

   ! submodels
   !IF (USE_XXX)   CALL  XXX_global_end

   ! infrastructure
   CALL main_tracer_global_end
!!$   CALL main_channel_global_end
   CALL main_qtimer_global_end

 END SUBROUTINE messy_global_end

 !==============================================================================

 SUBROUTINE messy_write_output

   ! infrastructure
   USE messy_main_switch       !ONLY: USE_*
   USE messy_main_timer,       ONLY: l_rerun
   USE messy_main_tracer_bi,   ONLY: main_tracer_write_output
   USE messy_main_channel_bi,  ONLY: main_channel_write_output
   USE messy_main_rnd_bi,      ONLY: main_rnd_write_output

   ! infrastructure
   CALL main_tracer_write_output(1)

   CALL main_rnd_write_output     ! save integer state vectors in channel objs
   CALL main_channel_write_output

 END SUBROUTINE messy_write_output

 !==============================================================================

 SUBROUTINE messy_write_restart

   USE messy_main_timer,        ONLY: l_rerun
   USE messy_main_channel_bi,   ONLY: main_channel_write_restart

   IF (l_rerun) THEN
      CALL main_channel_write_restart
   END IF

 END SUBROUTINE messy_write_restart

 !==============================================================================

 SUBROUTINE messy_free_memory

   ! infrastructure
   USE messy_main_switch        !ONLY: USE_*
   USE messy_main_qtimer_bi,    ONLY: main_qtimer_free_memory
   USE messy_main_tracer_bi,    ONLY: main_tracer_free_memory
   USE messy_main_channel_bi,   ONLY: main_channel_free_memory
   USE messy_main_grid_def_bi,  ONLY: main_grid_def_free_memory
   USE messy_main_data_bi,      ONLY: main_data_free_memory
   USE messy_main_import_bi,    ONLY: main_import_free_memory
!!$   USE messy_ncregrid_tools_bi, ONLY: messy_ncregrid_free_memory

   ! submodels
   USE messy_aeropt_si,      ONLY:     aeropt_free_memory
   USE messy_cloudopt_si,    ONLY:   cloudopt_free_memory

   IMPLICIT NONE

   ! submodels
   IF (USE_AEROPT)   CALL     aeropt_free_memory
   IF (USE_CLOUDOPT) CALL   cloudopt_free_memory

   ! infrastructure
   CALL main_tracer_free_memory
   CALL main_channel_free_memory
   CALL main_data_free_memory
   CALL main_grid_def_free_memory
!!$   CALL messy_ncregrid_free_memory
   CALL main_qtimer_free_memory

  END SUBROUTINE messy_free_memory

!==============================================================================

  SUBROUTINE messy_finalize

    USE messy_main_blather_bi, ONLY: main_blather_finalize

    IMPLICIT NONE

    CALL main_blather_finalize

  END SUBROUTINE messy_finalize

!==============================================================================
