MODULE MESSY_MAIN_CONTROL

  ! CONTROL CORE
  !  written by Astrid Kerkweg, February 2018

  ! DATA used to identify MESSy entry points by name or index

  USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM

  IMPLICIT NONE
  PUBLIC
  SAVE


  ! DEFINE IDENTIFIER FOR location in MESSy entry points
  INTEGER, PARAMETER :: MEPS_BEGIN = 1
  INTEGER, PARAMETER :: MEPS_END   = 2
  INTEGER, PARAMETER :: MEPS_ONE   = 0

  ! TRIGGER EVERYWHERE IN ENTRY POINT (REQUIRED FOR FUNCTIONS, to enable calls
  ! from subroutines without case differentiation)
  INTEGER, PARAMETER :: MEPS_ALWAYS= -1

  INTEGER, PARAMETER :: MEPS_MIN   = -1
  INTEGER, PARAMETER :: MEPS_MAX   = 2

  ! DEFINE IDENTIFYER FOR MESSy entry points
  INTEGER, PARAMETER :: MEP_UNDEF = -99

  INTEGER, PARAMETER :: MEP_SETUP          = 1
  INTEGER, PARAMETER :: MEP_INITIALIZE     = 2
  INTEGER, PARAMETER :: MEP_NEW_TRACER     = 3
  INTEGER, PARAMETER :: MEP_TRACER_META    = 4
  INTEGER, PARAMETER :: MEP_INIT_MEMORY    = 5
  INTEGER, PARAMETER :: MEP_INIT_COUPLER   = 6
  INTEGER, PARAMETER :: MEP_INIT_COUPLING  = 7
  INTEGER, PARAMETER :: MEP_READ_RESTART   = 8
  INTEGER, PARAMETER :: MEP_INIT_TRACER    = 9
  INTEGER, PARAMETER :: MEP_INIT_LOOP      = 10
  INTEGER, PARAMETER :: MEP_TIME           = 11
  INTEGER, PARAMETER :: MEP_TENDENCY_RESET = 12
  INTEGER, PARAMETER :: MEP_GLOBAL_START   = 13
  INTEGER, PARAMETER :: MEP_LOCAL_START    = 14
  INTEGER, PARAMETER :: MEP_BEFOREADV      = 15
  INTEGER, PARAMETER :: MEP_AFTERADV       = 16
  INTEGER, PARAMETER :: MEP_RADIATION      = 17
  INTEGER, PARAMETER :: MEP_VDIFF          = 18
  INTEGER, PARAMETER :: MEP_TURBDIFF       = 19
  INTEGER, PARAMETER :: MEP_RADHEAT        = 20
  INTEGER, PARAMETER :: MEP_GWDRAG         = 21
  INTEGER, PARAMETER :: MEP_CONVEC         = 22
  INTEGER, PARAMETER :: MEP_MIXLO          = 23
  INTEGER, PARAMETER :: MEP_PHYSC          = 24
  INTEGER, PARAMETER :: MEP_LOCAL_END      = 25
  INTEGER, PARAMETER :: MEP_GLOBAL_END     = 26
  INTEGER, PARAMETER :: MEP_WRITE_OUTPUT   = 27
  INTEGER, PARAMETER :: MEP_WRITE_RESTART  = 28
  INTEGER, PARAMETER :: MEP_FREE_MEMORY    = 29
  INTEGER, PARAMETER :: MEP_FINALIZE       = 30

  INTEGER, PARAMETER :: MEP_ENTRY_MAX      = 30

  CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(MEP_ENTRY_MAX) :: entry_name = &
       (/ 'messy_setup             ','messy_initialize        ' &
         ,'messy_new_tracer        ','messy_tracer_meta       ' &
         ,'messy_init_memory       ','messy_init_coupler      ' &
         ,'messy_init_coupling     ','messy_read_restart      ' &
         ,'messy_init_tracer       ','messy_init_loop         ' &
         ,'messy_time              ','messy_tendency_reset    ' &
         ,'messy_global_start      ','messy_local_start       ' &
         ,'messy_beforeadv         ','messy_afteradv          ' &
         ,'messy_radiation         ','messy_vdiff             ' &
         ,'messy_trubdiff          ' &
         ,'messy_radheat           ','messy_gwdrag            ' &
         ,'messy_convec            ','messy_mixlo             ' &
         ,'messy_physc             ','messy_local_end         ' &
         ,'messy_global_end        ','messy_write_output      ' &
         ,'messy_write_restart     ','messy_free_memory       ' &
         ,'messy_finalize          ' /)

  ! which entry point?
  INTEGER :: entrypoint = MEP_UNDEF
  ! sub-entry point, e.g., in CALL messy_setup(0)
  !                        0 indicates the sub-entry point
  INTEGER :: subentry   = MEP_UNDEF

  ! SET DEBUGING OUTPUT
  LOGICAL :: ldebug = .FALSE.


CONTAINS

  !------------------------------------------------------

  SUBROUTINE main_control_set_context(pos, sp)

    USE messy_main_constants_mem, ONLY: iouerr

    IMPLICIT NONE

    INTEGER, INTENT(IN)           :: pos
    INTEGER, INTENT(IN), OPTIONAL :: sp

    entrypoint = pos

    IF (PRESENT(sp)) subentry = sp

    IF (ldebug) write(iouerr,*) 'CONTROL START ',entry_name(entrypoint)

  END SUBROUTINE main_control_set_context

  !------------------------------------------------------

  SUBROUTINE main_control_unset_context

    USE messy_main_constants_mem, ONLY: iouerr

    IMPLICIT NONE

    IF (ldebug) write(iouerr,*) 'CONTROL END   ',entry_name(entrypoint)

    entrypoint = MEP_UNDEF
    subentry   = MEP_UNDEF

  END SUBROUTINE main_control_unset_context

  !------------------------------------------------------

  SUBROUTINE control_get_idx_from_name(status, idx, name)

    IMPLICIT NONE

    INTEGER,                      INTENT(OUT) :: status
    INTEGER,                      INTENT(OUT) :: idx
    CHARACTER(LEN=STRLEN_MEDIUM), INTENT(IN)  :: name

    ! LOCAL
    INTEGER :: ii

    status = -1

    DO ii = 1, MEP_ENTRY_MAX
       IF (TRIM(entry_name(ii)) == ADJUSTL(TRIM(name))) THEN
          idx = ii
          status = 0
          RETURN
       END IF
    END DO

  END SUBROUTINE control_get_idx_from_name

  !------------------------------------------------------

  !------------------------------------------------------

  SUBROUTINE main_control_set_context_output(lout)

    IMPLICIT NONE

    LOGICAL, INTENT(IN) :: lout


    ldebug = lout

  END SUBROUTINE main_control_set_context_output

  !------------------------------------------------------

  !------------------------------------------------------
  LOGICAL FUNCTION main_control_is_context( ent, subent)

    IMPLICIT NONE

    INTEGER, INTENT(IN)           :: ent
    INTEGER, INTENT(IN), OPTIONAL :: subent

    main_control_is_context = .FALSE.

    IF (ent == entrypoint) THEN
       IF (PRESENT(subent)) THEN
          IF (subent == subentry .OR. subent == MEPS_ALWAYS) THEN
              main_control_is_context = .TRUE.
           END IF
       ELSE
          main_control_is_context = .TRUE.
       END IF
    END IF

  END FUNCTION main_control_is_context

  !------------------------------------------------------

  !------------------------------------------------------

  CHARACTER(STRLEN_MEDIUM) FUNCTION main_control_get_context_name()

    IMPLICIT NONE

    main_control_get_context_name = entry_name(entrypoint)

  END FUNCTION main_control_get_context_name

  !------------------------------------------------------

  !------------------------------------------------------

  INTEGER FUNCTION main_control_get_context_num()

    IMPLICIT NONE

    main_control_get_context_num = entrypoint

  END FUNCTION main_control_get_context_num

  !------------------------------------------------------

END MODULE MESSY_MAIN_CONTROL
