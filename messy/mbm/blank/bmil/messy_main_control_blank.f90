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
  USE messy_main_blather_bi, ONLY: main_blather_setup
  USE messy_main_switch_bi,  ONLY: main_switch_setup
  USE messy_main_qtimer_bi,  ONLY: main_qtimer_setup
  USE messy_main_timer_bi,   ONLY: main_timer_setup
  USE messy_main_channel_bi, ONLY: main_channel_setup
  USE messy_main_grid_bi,    ONLY: main_grid_setup
  USE messy_main_mpi_bi,     ONLY: main_mpi_setup

  IMPLICIT NONE

  ! LOCAL
  INTEGER :: nclock = 1

  CALL main_grid_setup
  CALL main_mpi_setup
  ! infrastructure
  CALL main_blather_setup
  ! TIMER SPLIT NEEDED HERE, TO ENABLE OTHER MODELS TO OVERWRITE NAMELIST
  ! ENTRY OF START DATE AND STOP DATE (e.g. MMD_CLIENT)
  ! ub_ak_20171109+
  !CALL main_timer_setup(1)
  CALL main_timer_setup(1, nclock)
  ! ub_ak_20171109-
!!$  CALL main_timer_setup(2) ! op_pj_20120307
  CALL main_channel_setup
  ! ub_ak_20171109+
  !CALL main_timer_setup(2)    ! op_pj_20120307
  CALL main_timer_setup(2, nclock)    ! op_pj_20120307
  ! ub_ak_20171109-
  ! do NOT change the following order of main_* subroutines
  CALL main_qtimer_setup
  CALL main_switch_setup

END SUBROUTINE messy_setup

!==============================================================================

SUBROUTINE messy_initialize

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
#if defined(MBM_BLANK)
   USE messy_ptrac_si,       ONLY:   ptrac_initialize
   USE messy_submod3_si,     ONLY: submod3_initialize
   USE messy_orbit_si,       ONLY:   orbit_initialize
#endif
#if defined(MBM_QBO)
   USE messy_qbo_si,         ONLY:     qbo_initialize
#endif
#if defined(MBM_DISSOC)
   USE messy_ptrac_si,       ONLY:   ptrac_initialize
   USE messy_dissoc_si,      ONLY:  dissoc_initialize
#endif

   IMPLICIT NONE

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
#if defined (MBM_BLANK)
   IF (USE_PTRAC)   CALL   ptrac_initialize
   IF (USE_SUBMOD3) CALL submod3_initialize
   IF (USE_ORBIT)   CALL   orbit_initialize
#endif
#if defined (MBM_QBO)
   IF (USE_QBO)     CALL     qbo_initialize
#endif
#if defined (MBM_DISSOC)
   IF (USE_PTRAC)   CALL   ptrac_initialize
   IF (USE_DISSOC)  CALL  dissoc_initialize
#endif

 END SUBROUTINE messy_initialize

!==============================================================================

 SUBROUTINE messy_new_tracer

   ! infrastructure
   USE messy_main_switch      !ONLY: USE_*
   USE messy_main_tracer_bi,  ONLY: main_tracer_new_tracer

   ! submodels
#if defined (MBM_BLANK)
   USE messy_ptrac_si,        ONLY: ptrac_new_tracer
#endif
#if defined (MBM_DISSOC)
   USE messy_ptrac_si,        ONLY: ptrac_new_tracer
#endif

   IMPLICIT NONE

   CALL main_tracer_new_tracer(1)  ! define tracer set

   ! submodels
#if defined (MBM_BLANK)
   IF (USE_PTRAC) CALL ptrac_new_tracer
#endif
#if defined (MBM_DISSOC)
   IF (USE_PTRAC) CALL ptrac_new_tracer
#endif

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
#if defined (MBM_BLANK)
   USE messy_ptrac_si,        ONLY:    ptrac_init_memory
   USE messy_submod1_si,      ONLY:  submod1_init_memory
   USE messy_submod2_si,      ONLY:  submod2_init_memory
   USE messy_orbit_si,        ONLY:    orbit_init_memory
#endif
#if defined (MBM_QBO)
   USE messy_qbo_si,          ONLY:      qbo_init_memory
#endif
#if defined (MBM_DISSOC)
   USE messy_ptrac_si,        ONLY:    ptrac_init_memory
   USE messy_dissoc_si,       ONLY:   dissoc_init_memory
#endif

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
#if defined (MBM_BLANK)
   IF (USE_PTRAC)      CALL ptrac_init_memory
   IF (USE_SUBMOD1)    CALL submod1_init_memory
   IF (USE_SUBMOD2)    CALL submod2_init_memory
   IF (USE_ORBIT)      CALL   orbit_init_memory
#endif
#if defined (MBM_QBO)
   IF (USE_QBO)        CALL     qbo_init_memory
#endif
#if defined (MBM_DISSOC)
   IF (USE_PTRAC)      CALL   ptrac_init_memory
   IF (USE_DISSOC)     CALL  dissoc_init_memory
#endif

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

   ! submodels
#if defined(MBM_BLANK)
   USE messy_submod2_si,      ONLY: submod2_init_coupling
   USE messy_submod3_si,      ONLY: submod3_init_coupling
   USE messy_orbit_si,        ONLY:   orbit_init_coupling
   USE messy_testevent_si,    ONLY: testevent_init_coupling ! ub_ak_20190107
#endif
#if defined(MBM_QBO)
   USE messy_qbo_si,          ONLY:     qbo_init_coupling
#endif
#if defined(MBM_DISSOC)
   USE messy_dissoc_si,       ONLY:  dissoc_init_coupling
#endif

   IMPLICIT NONE

   ! infrastructure
   ! - resetting meta information of family-members (to tracers)
   !   (after advection initialisation)
   CALL main_tracer_init_coupling(1)
   ! - modify attributes according to ADD_ATT in CTRL namelist
   CALL main_channel_init_coupling(1)

   ! submodels
#if defined (MBM_BLANK)
   IF (USE_SUBMOD2)    CALL submod2_init_coupling
   IF (USE_SUBMOD3)    CALL submod3_init_coupling
   IF (USE_ORBIT)      CALL   orbit_init_coupling
   IF (USE_TESTEVENT)  CALL testevent_init_coupling ! ub_ak_20190107
#endif
#if defined (MBM_QBO)
   IF (USE_QBO)        CALL     qbo_init_coupling
#endif
#if defined(MBM_DISSOC)
   IF (USE_DISSOC)     CALL  dissoc_init_coupling
#endif

   CALL main_rnd_init_coupling     ! create channel objects for state vectors
   CALL main_channel_init_coupling(2)

 END SUBROUTINE messy_init_coupling

 !==============================================================================

 SUBROUTINE messy_read_restart

   USE messy_main_timer,        ONLY: lresume
   USE messy_main_channel_bi,   ONLY: main_channel_read_restart

   IF (.NOT. lresume) RETURN

   CALL main_channel_read_restart

 END SUBROUTINE messy_read_restart

!==============================================================================

 SUBROUTINE messy_init_tracer

   ! infrastructure
   USE messy_main_switch      ! ONLY: USE_*
   USE messy_main_tracer_bi,  ONLY: main_tracer_init_tracer
   USE messy_main_rnd_bi,     ONLY: main_rnd_init_tracer
   USE messy_main_import_bi,  ONLY: main_import_init_tracer

   ! submodels
#if defined (MBM_BLANK)
   USE messy_submod3_si,      ONLY: submod3_init_tracer
#endif

   IMPLICIT NONE

   ! infrastructure
   CALL main_rnd_init_tracer       ! reset state vectors after restart
   CALL main_tracer_init_tracer(1) ! check tracer init from rerun
   CALL main_import_init_tracer

   CALL main_tracer_init_tracer(4) ! initialize via tracer_init (tracer.nml)

   ! submodels
#if defined (MBM_BLANK)
   IF (USE_SUBMOD3)   CALL submod3_init_tracer
#endif

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
#if defined (MBM_BLANK)
   USE messy_submod3_si,      ONLY: submod3_global_start
   USE messy_orbit_si,        ONLY:   orbit_global_start
#endif
#if defined (MBM_DISSOC)
   USE messy_dissoc_si,       ONLY:  dissoc_global_start
#endif

   IMPLICIT NONE

   ! infrastructure
   CALL main_data_global_start
   CALL main_tracer_global_start
   CALL main_channel_global_start
   CALL main_import_global_start

   ! submodels
#if defined (MBM_BLANK)
   IF (USE_SUBMOD3) CALL submod3_global_start
   IF (USE_ORBIT)   CALL   orbit_global_start
#endif
#if defined (MBM_DISSOC)
   IF (USE_DISSOC)  CALL   dissoc_global_start
#endif

 END SUBROUTINE messy_global_start

!============================================================================

!============================================================================

 SUBROUTINE messy_local_start

   ! infrastructure
   USE messy_main_switch     !ONLY: USE_*
   USE messy_main_grid_bi,   ONLY: main_grid_local_start
   USE messy_main_tracer_bi, ONLY: main_tracer_local_start

   ! submodels
   ! USE messy_XXX_si,     ONLY: XXX_local_start

   IMPLICIT NONE

   ! infrastructure
   CALL main_grid_local_start
   CALL main_tracer_local_start

   ! submodels
   !IF (USE_XXX) CALL XXX_local_start

 END SUBROUTINE messy_local_start

 !==============================================================================

 !==============================================================================

 SUBROUTINE messy_physc

   ! infrastructure
   USE messy_main_switch     !ONLY: USE_*
   USE messy_main_grid_def_mem_bi, ONLY: jrow, nlat, ngpblks, kproma
   USE messy_main_mpi_bi,    ONLY: p_pe
   USE messy_main_timer,     ONLY: YEAR, MONTH, DAY, HOUR, MINUTE, SECOND
   USE messy_main_tracer_bi, ONLY: main_tracer_local_start

   ! submodels
#if defined (MBM_BLANK)
    USE messy_submod2_si,       ONLY:  submod2_physc
    USE messy_submod3_si,       ONLY:  submod3_physc
#endif
#if defined (MBM_QBO)
   USE messy_qbo_si,            ONLY:      qbo_physc
#endif
#if defined (MBM_DISSOC)
   USE messy_dissoc_si,         ONLY:   dissoc_physc
#endif

   IMPLICIT NONE

   ! infrastructure
   ! ONLY FOR DEBUGGING +
   IF (L_TIME_INFO) &
   WRITE(*, '(1x,i4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a3'//&
        &',a4,i4,a8,i5,a3,i5,a10,i5)') &
        YEAR,'-',MONTH,'-',DAY,' ',HOUR,':', MINUTE,':',SECOND,' # ' &
        ,' PE=',p_pe,' # JROW=',jrow,' / ',ngpblks,' # KPROMA=',kproma
   ! ONLY FOR DEBUGGING -

   CALL main_tracer_local_start ! must be called first

   ! submodels
#if defined (MBM_BLANK)
   IF (USE_SUBMOD2) CALL submod2_physc
   IF (USE_SUBMOD3) CALL submod3_physc
#endif
#if defined (MBM_QBO)
   IF (USE_QBO)     CALL     qbo_physc
#endif
#if defined (MBM_DISSOC)
   IF (USE_DISSOC) CALL   dissoc_physc
#endif

 END SUBROUTINE messy_physc
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

!==============================================================================

 SUBROUTINE messy_global_end

   ! infrastructure
   USE messy_main_switch      !ONLY: USE_*
   USE messy_main_qtimer_bi,  ONLY: main_qtimer_global_end
   USE messy_main_tracer_bi,  ONLY: main_tracer_global_end
!!$   USE messy_main_channel_bi, ONLY: main_channel_global_end

   ! submodels
#if defined (MBM_BLANK)
   USE messy_submod1_si,      ONLY: submod1_global_end
   USE messy_submod2_si,      ONLY: submod2_global_end
   USE messy_testevent_si,    ONLY: testevent_global_end ! ub_ak_20190107
#endif
#if defined (MBM_DISSOC)
   USE messy_dissoc_si,       ONLY:  dissoc_global_end
#endif

   IMPLICIT NONE

   ! infrastructure

   ! submodels
#if defined (MBM_BLANK)
   IF (USE_SUBMOD1)   CALL  submod1_global_end
   IF (USE_SUBMOD2)   CALL  submod2_global_end
   IF (USE_TESTEVENT) CALL  testevent_global_end ! ub_ak_20190107
#endif
#if defined (MBM_DISSOC)
   IF (USE_DISSOC)    CALL   dissoc_global_end
#endif

   ! infrastructure
   CALL main_tracer_global_end
!!$   CALL main_channel_global_end
   CALL main_qtimer_global_end

 END SUBROUTINE messy_global_end

 !==============================================================================

 SUBROUTINE messy_write_output

   ! infrastructure
   USE messy_main_switch       !ONLY: USE_*
   USE messy_main_tracer_bi,   ONLY: main_tracer_write_output
   USE messy_main_channel_bi,  ONLY: main_channel_write_output
   USE messy_main_rnd_bi,      ONLY: main_rnd_write_output

   ! infrastructure
   CALL main_tracer_write_output(1)

   ! submodels

   CALL main_rnd_write_output     ! save integer state vectors in channel objs
   CALL main_channel_write_output

 END SUBROUTINE messy_write_output

 !==============================================================================

 SUBROUTINE messy_write_restart

   USE messy_main_timer,       ONLY: l_rerun
   USE messy_main_channel_bi,  ONLY: main_channel_write_restart

   IF (l_rerun) CALL main_channel_write_restart

 END SUBROUTINE messy_write_restart

 !==============================================================================

 SUBROUTINE messy_free_memory

   ! infrastructure
   USE messy_main_switch        !ONLY: USE_*
   USE messy_main_qtimer_bi,    ONLY: main_qtimer_free_memory
   USE messy_main_tracer_bi,    ONLY: main_tracer_free_memory
   USE messy_main_channel_bi,   ONLY: main_channel_free_memory
   USE messy_main_data_bi,      ONLY: main_data_free_memory
   USE messy_main_grid_bi,      ONLY: main_grid_free_memory
   USE messy_main_import_bi,    ONLY: main_import_free_memory

   ! submodels
#if defined (MBM_BLANK)
   USE messy_submod1_si,        ONLY: submod1_free_memory
   USE messy_orbit_si,          ONLY:   orbit_free_memory
#endif
#if defined (MBM_QBO)
   USE messy_qbo_si,            ONLY:     qbo_free_memory
#endif
#if defined (MBM_DISSOC)
   USE messy_dissoc_si,         ONLY:  dissoc_free_memory
#endif

   IMPLICIT NONE

   ! submodels
#if defined (MBM_BLANK)
   IF (USE_SUBMOD1) CALL submod1_free_memory
   IF (USE_ORBIT)   CALL   orbit_free_memory
#endif
#if defined (MBM_QBO)
   IF (USE_QBO)     CALL     qbo_free_memory
#endif
#if defined (MBM_DISSOC)
   IF (USE_DISSOC)  CALL  dissoc_free_memory
#endif

   ! infrastructure
   CALL main_import_free_memory
   CALL main_tracer_free_memory
   CALL main_channel_free_memory
   CALL main_data_free_memory
   CALL main_grid_free_memory
   CALL main_qtimer_free_memory

  END SUBROUTINE messy_free_memory

!==============================================================================

  SUBROUTINE messy_finalize

    USE messy_main_blather_bi, ONLY: main_blather_finalize
 
    IMPLICIT NONE

    CALL main_blather_finalize

  END SUBROUTINE messy_finalize

!==============================================================================
