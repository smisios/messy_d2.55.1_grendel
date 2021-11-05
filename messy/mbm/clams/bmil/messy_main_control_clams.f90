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

  ! submodels
  USE messy_clams_si,        ONLY: clams_setup
  USE messy_clamstraj_si,    ONLY: clamstraj_setup
  USE messy_clamssedi_si,    ONLY: clamssedi_setup
  USE messy_clamsmix_si,     ONLY: clamsmix_setup
  USE messy_clamsdeepconv_si,ONLY: clamsdeepconv_setup

  IMPLICIT NONE

  ! LOCAL
  INTEGER :: nclock = 1

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
  !CALL main_timer_setup(2)
  CALL main_timer_setup(2, nclock)
  ! ub_ak_20171109-
  ! do NOT change the following order of main_* subroutines
  CALL main_qtimer_setup
  CALL main_switch_setup

  CALL clams_setup
  if (USE_CLAMSTRAJ)     CALL clamstraj_setup
  if (USE_CLAMSSEDI)     CALL clamssedi_setup
  if (USE_CLAMSMIX)      CALL clamsmix_setup
  if (USE_CLAMSDEEPCONV) CALL clamsdeepconv_setup

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
   USE messy_clams_si,         ONLY: clams_initialize
   USE messy_clamstraj_si,     ONLY: clamstraj_initialize
   USE messy_clamsdeepconv_si, ONLY: clamsdeepconv_initialize
   USE messy_clamscirrus_si,   ONLY: clamscirrus_initialize
   USE messy_dissoc_si,        ONLY: dissoc_initialize
   USE messy_clamschem_si,     ONLY: clamschem_initialize
   USE messy_clamssedi_si,     ONLY: clamssedi_initialize
   USE messy_clamsmix_si,      ONLY: clamsmix_initialize
   USE messy_clamsbmix_si,     ONLY: clamsbmix_initialize
   USE messy_clamstracer_si,   ONLY: clamstracer_initialize

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
   CALL clams_initialize
   IF (USE_CLAMSTRAJ) CALL clamstraj_initialize
   IF (USE_CLAMSDEEPCONV) CALL clamsdeepconv_initialize
   IF (USE_DISSOC)    CALL dissoc_initialize
   IF (USE_CLAMSCHEM) CALL clamschem_initialize
   IF (USE_CLAMSCIRRUS) CALL clamscirrus_initialize
   IF (USE_CLAMSSEDI)  CALL clamssedi_initialize
   IF (USE_CLAMSMIX)  CALL clamsmix_initialize
   IF (USE_CLAMSBMIX) CALL clamsbmix_initialize
   IF (USE_CLAMSTRACER) CALL clamstracer_initialize

 END SUBROUTINE messy_initialize

!==============================================================================

 SUBROUTINE messy_new_tracer

   ! infrastructure
   USE messy_main_switch      !ONLY: USE_*
   USE messy_main_tracer_bi,  ONLY: main_tracer_new_tracer

   ! submodels

   IMPLICIT NONE

   CALL main_tracer_new_tracer(1)  ! define tracer set

   ! submodels

   CALL main_tracer_new_tracer(2) ! define tracer families
   CALL main_tracer_new_tracer(3) ! diagnostic output

 END SUBROUTINE messy_new_tracer

!==============================================================================

 SUBROUTINE messy_init_memory

   ! infrastructure
   USE messy_main_switch      !ONLY: USE_*
   USE messy_main_qtimer_bi,  ONLY: main_qtimer_init_memory
   USE messy_main_data_bi,    ONLY: main_data_init_memory
   USE messy_main_tracer_bi,  ONLY: main_tracer_init_memory
   USE messy_main_channel_bi, ONLY: main_channel_init_memory
   USE messy_main_import_bi,  ONLY: main_import_init_memory

   ! submodels
   USE messy_clams_si,        ONLY: clams_init_memory
   USE messy_clamstraj_si,    ONLY: clamstraj_init_memory
   USE messy_clamsdeepconv_si,ONLY: clamsdeepconv_init_memory
   USE messy_clamscirrus_si,  ONLY: clamscirrus_init_memory
   USE messy_dissoc_si,       ONLY: dissoc_init_memory
   USE messy_clamschem_si,    ONLY: clamschem_init_memory
   USE messy_clamssedi_si,    ONLY: clamssedi_init_memory
   USE messy_clamsmix_si,     ONLY: clamsmix_init_memory
   USE messy_clamsbmix_si,    ONLY: clamsbmix_init_memory
   USE messy_clamstracer_si,  ONLY: clamstracer_init_memory

   IMPLICIT NONE

   ! infrastructure
   CALL main_qtimer_init_memory
   CALL main_channel_init_memory   ! setup tracer memory
   CALL main_tracer_init_memory(1) ! create data channel BML <-> SMIL
   CALL main_data_init_memory
   CALL main_import_init_memory

   ! submodels
   ! write(*,*)'init_mem: ',USE_CLAMSTRAJ, USE_CLAMSCIRRUS, USE_DISSOC,&
   !      USE_CLAMSCHEM, USE_CLAMSMIX, USE_CLAMSBMIX
   CALL clams_init_memory
   IF (USE_CLAMSTRAJ)     CALL clamstraj_init_memory
   If (USE_CLAMSDEEPCONV) CALL clamsdeepconv_init_memory
   IF (USE_DISSOC)        CALL dissoc_init_memory
   IF (USE_CLAMSCHEM)     CALL clamschem_init_memory
   IF (USE_CLAMSCIRRUS)   CALL clamscirrus_init_memory
   IF (USE_CLAMSSEDI)     CALL clamssedi_init_memory
   IF (USE_CLAMSMIX)      CALL clamsmix_init_memory
   IF (USE_CLAMSBMIX)     CALL clamsbmix_init_memory
   IF (USE_CLAMSTRACER)   CALL clamstracer_init_memory

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

   ! submodels
   USE messy_clams_si,         ONLY: clams_init_coupling
   USE messy_clamstraj_si,     ONLY: clamstraj_init_coupling
   USE messy_clamsdeepconv_si, ONLY: clamsdeepconv_init_coupling
   USE messy_clamscirrus_si,   ONLY: clamscirrus_init_coupling
   USE messy_dissoc_si,        ONLY: dissoc_init_coupling
   USE messy_clamschem_si,     ONLY: clamschem_init_coupling
   USE messy_clamssedi_si,     ONLY: clamssedi_init_coupling
   USE messy_clamsmix_si,      ONLY: clamsmix_init_coupling
   USE messy_clamsbmix_si,     ONLY: clamsbmix_init_coupling
   USE messy_clamstracer_si,   ONLY: clamstracer_init_coupling

   IMPLICIT NONE

   ! infrastructure
   ! - resetting meta information of family-members (to tracers)
   !   (after advection initialisation)
   CALL main_tracer_init_coupling(1)
   ! - modify attributes according to ADD_ATT in CTRL namelist
   CALL main_channel_init_coupling(1)

   ! submodels
   IF (USE_CLAMS)         CALL clams_init_coupling(1)
   IF (USE_CLAMSTRAJ)     CALL clamstraj_init_coupling
   IF (USE_CLAMSDEEPCONV) CALL clamsdeepconv_init_coupling
   IF (USE_CLAMSCIRRUS)   CALL clamscirrus_init_coupling
   IF (USE_DISSOC)        CALL dissoc_init_coupling
   IF (USE_CLAMSCHEM)     CALL clamschem_init_coupling
   IF (USE_CLAMSSEDI)     CALL clamssedi_init_coupling
   IF (USE_CLAMSMIX)      CALL clamsmix_init_coupling
   IF (USE_CLAMSBMIX)     CALL clamsbmix_init_coupling
   IF (USE_CLAMSTRACER)   CALL clamstracer_init_coupling

   IF (USE_CLAMS)         CALL clams_init_coupling(2)

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
   USE messy_main_import_bi, ONLY: main_import_init_tracer

   ! submodels

   IMPLICIT NONE

   ! infrastructure
   CALL main_tracer_init_tracer(1) ! check tracer init from rerun
   CALL main_import_init_tracer

   CALL main_tracer_init_tracer(4) ! initialize via tracer_init (tracer.nml)

   ! submodels

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
   USE messy_clams_si,        ONLY: clams_global_start
   USE messy_clamstraj_si,    ONLY: clamstraj_global_start
   USE messy_dissoc_si,       ONLY: dissoc_global_start
   USE messy_clamssedi_si,    ONLY: clamssedi_global_start
   USE messy_clamsmix_si,     ONLY: clamsmix_global_start
   USE messy_clamstracer_si,  ONLY: clamstracer_global_start

   IMPLICIT NONE

   ! infrastructure
   CALL main_data_global_start
   CALL main_tracer_global_start
   CALL main_channel_global_start
   CALL main_import_global_start

   ! submodels
   CALL clams_global_start
   IF (USE_CLAMSSEDI)    CALL clamssedi_global_start
   IF (USE_CLAMSTRAJ)    CALL clamstraj_global_start
   IF (USE_DISSOC)       CALL dissoc_global_start
   IF (USE_CLAMSMIX)     CALL clamsmix_global_start
   IF (USE_CLAMSTRACER)  CALL clamstracer_global_start

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

   IMPLICIT NONE

   ! infrastructure
   ! ONLY FOR DEBUGGING +
   !if (p_pe==0) WRITE(*, '(1x,a,i4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a3'//&
   !     &',a4,i4,a8,i5,a3,i5,a10,i5)') &
   !     'messy_physc: ',YEAR,'-',MONTH,'-',DAY,' ',HOUR,':', MINUTE,':',SECOND,' # ' &
   !     ,' PE=',p_pe,' # JROW=',jrow,' / ',ngpblks,' # KPROMA=',kproma
   ! ONLY FOR DEBUGGING -

   CALL main_tracer_local_start ! must be called first

   ! submodels

 END SUBROUTINE messy_physc

!==============================================================================
 SUBROUTINE messy_local_end

   ! infrastructure
   USE messy_main_switch     !ONLY: USE_*

   IMPLICIT NONE

   ! submodels

 END SUBROUTINE messy_local_end

!==============================================================================

 SUBROUTINE messy_global_end

   ! infrastructure
   USE messy_main_switch      !ONLY: USE_*
   USE messy_main_qtimer_bi,  ONLY: main_qtimer_global_end
   USE messy_main_tracer_bi,  ONLY: main_tracer_global_end

   ! submodels
   USE messy_clamstraj_si,     ONLY: clamstraj_global_end
   USE messy_clamsdeepconv_si, ONLY: clamsdeepconv_global_end
   USE messy_clamscirrus_si,   ONLY: clamscirrus_global_end
   USE messy_dissoc_si,        ONLY: dissoc_global_end
   USE messy_clamschem_si,     ONLY: clamschem_global_end
   USE messy_clamssedi_si,     ONLY: clamssedi_global_end
   USE messy_clamsmix_si,      ONLY: clamsmix_global_end
   USE messy_clamsbmix_si,     ONLY: clamsbmix_global_end
   USE messy_clams_si,         ONLY: clams_global_end
   USE messy_clamstracer_si,   ONLY: clamstracer_global_end

   USE messy_clams_global,     ONLY: lchemevent

   IMPLICIT NONE

   ! infrastructure

   ! submodels
   IF (USE_CLAMS)         CALL clams_global_end(1)
   IF (USE_CLAMSDEEPCONV) CALL clamsdeepconv_global_end
   IF (USE_CLAMSSEDI)     CALL clamssedi_global_end
   IF (USE_CLAMSTRAJ)     CALL clamstraj_global_end
   IF (USE_CLAMSCIRRUS)   CALL clamscirrus_global_end
!!!!!
   IF (USE_DISSOC)        CALL dissoc_global_end
!   IF (USE_DISSOC .and. lchemevent) CALL dissoc_global_end
   IF (USE_CLAMSCHEM)     CALL clamschem_global_end
   IF (USE_CLAMSMIX)      CALL clamsmix_global_end
   IF (USE_CLAMSBMIX)     CALL clamsbmix_global_end
   IF (USE_CLAMSTRACER)   CALL clamstracer_global_end

   ! infrastructure
   CALL main_tracer_global_end
   CALL main_qtimer_global_end

 END SUBROUTINE messy_global_end

 !==============================================================================

 SUBROUTINE messy_write_output

   ! infrastructure
   USE messy_main_switch       !ONLY: USE_*
   USE messy_main_timer,       ONLY: l_rerun
   USE messy_main_tracer_bi,   ONLY: main_tracer_write_output
   USE messy_main_channel_bi,  ONLY: main_channel_write_output &
                                   , IOMODE_OUT, IOMODE_RST

   ! infrastructure
   CALL main_tracer_write_output(1)
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
   USE messy_clams_si,          ONLY: clams_free_memory
   USE messy_clamstraj_si,      ONLY: clamstraj_free_memory
   USE messy_clamsdeepconv_si,  ONLY: clamsdeepconv_free_memory
   USE messy_clamschem_si,      ONLY: clamschem_free_memory
   USE messy_clamssedi_si,      ONLY: clamssedi_free_memory
   USE messy_clamsmix_si,       ONLY: clamsmix_free_memory
   USE messy_clamsbmix_si,      ONLY: clamsbmix_free_memory
   USE messy_clamscirrus_si,    ONLY: clamscirrus_free_memory
   USE messy_clamstracer_si,    ONLY: clamstracer_free_memory
   USE messy_dissoc_si,         ONLY: dissoc_free_memory

   IMPLICIT NONE

   ! submodels
   CALL clams_free_memory
   IF (USE_CLAMSCIRRUS)   CALL clamscirrus_free_memory
   IF (USE_CLAMSMIX)      CALL clamsmix_free_memory
   IF (USE_CLAMSBMIX)     CALL clamsbmix_free_memory
   IF (USE_CLAMSCHEM)     CALL clamschem_free_memory
   IF (USE_CLAMSSEDI)     CALL clamssedi_free_memory
   IF (USE_CLAMSTRAJ)     CALL clamstraj_free_memory
   IF (USE_CLAMSDEEPCONV) CALL clamsdeepconv_free_memory
   IF (USE_CLAMSTRACER)   CALL clamstracer_free_memory
   IF (USE_DISSOC)        CALL dissoc_free_memory

   ! infrastructure
   CALL main_tracer_free_memory
   CALL main_channel_free_memory
!!$   CALL main_data_free_memory
   CALL main_qtimer_free_memory

  END SUBROUTINE messy_free_memory

!==============================================================================

  SUBROUTINE messy_finalize

    USE messy_main_blather_bi, ONLY: main_blather_finalize

    IMPLICIT NONE

    CALL main_blather_finalize

  END SUBROUTINE messy_finalize

!==============================================================================
