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
!!$  CALL main_timer_setup(2) ! op_pj_20120307
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
  USE messy_ptrac_si,       ONLY: ptrac_initialize
  USE messy_ddep_si,        ONLY: ddep_initialize
  USE messy_mxl_si,         ONLY: mxl_initialize
  USE messy_megan_si,       ONLY: megan_initialize
  USE messy_onemis_si,      ONLY: onemis_initialize
  USE messy_offemis_si,     ONLY: offemis_initialize
  USE messy_jval_si,        ONLY: jval_initialize
  USE messy_mecca_si,       ONLY: mecca_initialize
  USE messy_oracle_si,      ONLY: oracle_initialize
  USE messy_chemglue_si,    ONLY: chemglue_initialize

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
  IF (USE_PTRAC)    CALL ptrac_initialize
  IF (USE_DDEP)     CALL ddep_initialize
  IF (USE_MXL)      CALL mxl_initialize
  IF (USE_MEGAN)    CALL megan_initialize
  IF (USE_ONEMIS)   CALL onemis_initialize
  IF (USE_OFFEMIS)  CALL offemis_initialize
  IF (USE_JVAL)     CALL jval_initialize
  IF (USE_MECCA)    CALL mecca_initialize
  IF (USE_ORACLE)   CALL oracle_initialize
  IF (USE_CHEMGLUE) CALL chemglue_initialize

END SUBROUTINE messy_initialize

!==============================================================================

 SUBROUTINE messy_new_tracer

   ! infrastructure
   USE messy_main_switch      !ONLY: USE_*
   USE messy_main_tracer_bi,  ONLY: main_tracer_new_tracer

   ! submodels
   USE messy_ptrac_si,        ONLY: ptrac_new_tracer
   USE messy_mecca_si,        ONLY: mecca_new_tracer
   USE messy_oracle_si,       ONLY: oracle_new_tracer

   IMPLICIT NONE

   CALL main_tracer_new_tracer(1)  ! define tracer set

   ! submodels
   IF (USE_PTRAC)  CALL ptrac_new_tracer
   IF (USE_MECCA)  CALL mecca_new_tracer
   IF (USE_ORACLE) CALL oracle_new_tracer

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
   USE messy_ptrac_si,        ONLY:  ptrac_init_memory
   USE messy_ddep_si,         ONLY:  ddep_init_memory
   USE messy_mxl_si,          ONLY:  mxl_init_memory
   USE messy_megan_si,        ONLY:  megan_init_memory
   USE messy_onemis_si,       ONLY:  onemis_init_memory
   USE messy_jval_si,         ONLY:  jval_init_memory
   USE messy_mecca_si,        ONLY:  mecca_init_memory
   USE messy_oracle_si,       ONLY:  oracle_init_memory
   USE messy_chemglue_si,     ONLY:  chemglue_init_memory

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
   IF (USE_PTRAC)      CALL ptrac_init_memory
   IF (USE_DDEP)       CALL ddep_init_memory
   IF (USE_MXL)        CALL mxl_init_memory
   IF (USE_MEGAN)      CALL megan_init_memory
   IF (USE_ONEMIS)     CALL onemis_init_memory
   IF (USE_JVAL)       CALL jval_init_memory
   IF (USE_MECCA)      CALL mecca_init_memory
   IF (USE_ORACLE)     CALL oracle_init_memory
   IF (USE_CHEMGLUE)   CALL chemglue_init_memory

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
!  USE messy_main_rnd_bi,     ONLY: main_rnd_init_coupling

   ! submodels
   USE messy_ddep_si,         ONLY:  ddep_init_coupling
   USE messy_megan_si,        ONLY:  megan_init_coupling
   USE messy_onemis_si,       ONLY:  onemis_init_coupling
   USE messy_offemis_si,      ONLY:  offemis_init_coupling
   USE messy_jval_si,         ONLY:  jval_init_coupling
   USE messy_mecca_si,        ONLY:  mecca_init_coupling
   USE messy_oracle_si,       ONLY:  oracle_init_coupling
   !USE messy_chemglue_si,     ONLY: chemglue_init_coupling

   IMPLICIT NONE

   ! infrastructure
   ! - resetting meta information of family-members (to tracers)
   !   (after advection initialisation)
   CALL main_tracer_init_coupling(1)
   ! - modify attributes according to ADD_ATT in CTRL namelist
   CALL main_channel_init_coupling(1)

   ! submodels
   IF (USE_DDEP)     CALL ddep_init_coupling
   IF (USE_MEGAN)    CALL megan_init_coupling
   IF (USE_ONEMIS)   CALL onemis_init_coupling
   IF (USE_OFFEMIS)  CALL offemis_init_coupling
   IF (USE_JVAL)     CALL jval_init_coupling
   IF (USE_MECCA)    CALL mecca_init_coupling(1)
   IF (USE_MECCA)    CALL mecca_init_coupling(2)
   IF (USE_ORACLE)   CALL oracle_init_coupling
   !IF (USE_CHEMGLUE) CALL chemglue_init_coupling

!  CALL main_rnd_init_coupling     ! create channel objects for state vectors
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
   USE messy_main_switch   ! ONLY: USE_*
   USE messy_main_tracer_bi, ONLY: main_tracer_init_tracer
!  USE messy_main_rnd_bi,    ONLY: main_rnd_init_tracer
   USE messy_main_import_bi, ONLY: main_import_init_tracer

   ! submodels
   USE messy_mecca_si,      ONLY: mecca_init_tracer
   USE messy_mxl_si,        ONLY: mxl_init_tracer

   IMPLICIT NONE

   ! infrastructure
!  CALL main_rnd_init_tracer       ! reset state vectors after restart
   CALL main_tracer_init_tracer(1) ! check tracer init from rerun
   CALL main_import_init_tracer

   ! only for ECHAM5
   CALL main_tracer_init_tracer(4) ! initialize via tracer_init (tracer.nml)

   IF (USE_MXL)     CALL mxl_init_tracer
   ! submodels
   IF (USE_MECCA)   CALL mecca_init_tracer

   ! infrastructure
! ECHAM5 only; CALL main_tracer_init_tracer(2) ! check tracer init by init_tracer
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
   USE messy_ddep_si,      ONLY: ddep_global_start
   USE messy_mxl_si,       ONLY: mxl_global_start
   USE messy_onemis_si,    ONLY: onemis_global_start
   USE messy_jval_si,      ONLY: jval_global_start
   !USE messy_chemglue_si,   ONLY: chemglue_global_start

   IMPLICIT NONE

   ! infrastructure
   CALL main_data_global_start
   CALL main_tracer_global_start
   CALL main_channel_global_start
   CALL main_import_global_start

   ! submodels
   IF (USE_MXL)    CALL mxl_global_start
   if (USE_DDEP)   call ddep_global_start
   if (USE_ONEMIS) call onemis_global_start
   if (USE_JVAL)   call jval_global_start
   !IF (USE_CHEMGLUE) CALL chemglue_global_start

 END SUBROUTINE messy_global_start

!============================================================================

!============================================================================

 SUBROUTINE messy_local_start

   ! infrastructure
   USE messy_main_switch     !ONLY: USE_*
   USE messy_main_data_bi,   ONLY: main_data_local_start
   USE messy_main_grid_bi,   ONLY: main_grid_local_start
   USE messy_main_tracer_bi, ONLY: main_tracer_local_start

   ! submodels
   ! USE messy_XXX_si,     ONLY: XXX_local_start

   IMPLICIT NONE

   CALL main_grid_local_start
   CALL main_data_local_start
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
   !USE messy_ptrac_si,      ONLY:  ptrac_physc
   USE messy_jval_si,      ONLY:  jval_physc
   USE messy_mecca_si,     ONLY:  mecca_physc
   USE messy_oracle_si,    ONLY:  oracle_physc
   USE messy_chemglue_si,  ONLY: chemglue_physc

   IMPLICIT NONE

   ! infrastructure
   ! ONLY FOR DEBUGGING +
   WRITE(*, '(1x,i4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a3'//&
        &',a4,i4,a8,i5,a3,i5,a10,i5)') &
        YEAR,'-',MONTH,'-',DAY,' ',HOUR,':', MINUTE,':',SECOND,' # ' &
        ,' PE=',p_pe,' # JROW=',jrow,' / ',ngpblks,' # KPROMA=',kproma
   ! ONLY FOR DEBUGGING -

   ! special requirements:
   ! - chemglue_physc must be called before mecca_physc

   CALL main_tracer_local_start ! must be called first

   ! submodels
   IF (USE_JVAL)     CALL jval_physc
   IF (USE_CHEMGLUE) CALL chemglue_physc
   IF (USE_MECCA)    CALL mecca_physc
   IF (USE_ORACLE)   CALL oracle_physc

 END SUBROUTINE messy_physc

 !==============================================================================

 SUBROUTINE messy_vdiff
   ! infrastructure
   USE messy_main_switch      !ONLY: USE_*

   ! submodels
   USE messy_ddep_si,         ONLY:  ddep_vdiff
   USE messy_megan_si,        ONLY:  megan_vdiff
   USE messy_onemis_si,       ONLY:  onemis_vdiff
   USE messy_offemis_si,      ONLY:  offemis_vdiff
   USE messy_mecca_si,        ONLY:  mecca_vdiff
   USE messy_oracle_si,       ONLY:  oracle_vdiff

   IMPLICIT NONE
  ! special requirements:
  ! - call mecca_vdiff(1) after emissions (ONLEM, OFFLEM, OFFEMIS, AIRSEA,
  !   MEGAN)
  !   and before deposition (DRYDEP) to get the "real" new emitted aerosol
  ! - call mecca_vdiff(2) after DRYDEP
  ! - call m7_vdiff and gmxe_vdiff after emissions
  !   (ONLEM, OFFLEM, OFFEMIS, AIRSEA, MEGAN) to get the
  !   current emissions if m7/gmxe calculates its own emissions
  !   to get the
  !   current emissions if cam does not calculate its own emissions
  ! - tropop_vdiff must be called before viso_vdiff
  ! - viso_vdiff and tropop_vdiff
  !   must be called first, since they calulate diagnostic quantities
  !   needed in other submodels

   ! submodels
   IF (USE_MEGAN)    CALL megan_vdiff
   IF (USE_ONEMIS)   CALL onemis_vdiff
   IF (USE_OFFEMIS)  CALL offemis_vdiff
   IF (USE_MECCA)    CALL mecca_vdiff(1)
   IF (USE_DDEP)     CALL ddep_vdiff
   IF (USE_MECCA)    CALL mecca_vdiff(2)
   IF (USE_ORACLE)   CALL oracle_vdiff

 END SUBROUTINE messy_vdiff

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
   USE messy_megan_si,      ONLY: megan_global_end
   USE messy_onemis_si,     ONLY: onemis_global_end
   USE messy_offemis_si,    ONLY: offemis_global_end
   USE messy_jval_si,       ONLY: jval_global_end
   USE messy_mecca_si,      ONLY: mecca_global_end

   IMPLICIT NONE

   ! infrastructure

   ! submodels
   IF (USE_MEGAN)    CALL megan_global_end
   IF (USE_ONEMIS)   CALL onemis_global_end
   IF (USE_OFFEMIS)  CALL offemis_global_end
   IF (USE_JVAL)     CALL jval_global_end
   IF (USE_MECCA)    CALL mecca_global_end

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
!  USE messy_main_rnd_bi,      ONLY: main_rnd_write_output

   ! infrastructure
   CALL main_tracer_write_output(1)

!  CALL main_rnd_write_output     ! save integer state vectors in channel objs
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
!!$   USE messy_main_data_bi,      ONLY: main_data_free_memory
   USE messy_main_import_bi,    ONLY: main_import_free_memory

   ! submodels
   USE messy_ddep_si,           ONLY: ddep_free_memory
   USE messy_mxl_si,            ONLY: mxl_free_memory
   USE messy_megan_si,          ONLY: megan_free_memory
   USE messy_onemis_si,         ONLY: onemis_free_memory
   USE messy_offemis_si,        ONLY: offemis_free_memory
   USE messy_jval_si,           ONLY: jval_free_memory
   USE messy_mecca_si,          ONLY: mecca_free_memory
   USE messy_oracle_si,         ONLY: oracle_free_memory

   IMPLICIT NONE

   ! submodels
   IF (USE_MXL)     CALL mxl_free_memory
   IF (USE_DDEP)    CALL ddep_free_memory
   IF (USE_MEGAN)   CALL megan_free_memory
   IF (USE_ONEMIS)  CALL onemis_free_memory
   IF (USE_OFFEMIS) CALL offemis_free_memory
   IF (USE_JVAL)    CALL jval_free_memory
   IF (USE_MECCA)   CALL mecca_free_memory
   IF (USE_ORACLE)  CALL oracle_free_memory

   ! infrastructure
   CALL main_tracer_free_memory
   CALL main_channel_free_memory
   CALL main_channel_free_memory
!!$   CALL main_data_free_memory
   CALL main_import_free_memory
   CALL main_qtimer_free_memory

  END SUBROUTINE messy_free_memory

!==============================================================================

  SUBROUTINE messy_finalize

    USE messy_main_blather_bi, ONLY: main_blather_finalize

    IMPLICIT NONE

    CALL main_blather_finalize

  END SUBROUTINE messy_finalize

!==============================================================================
