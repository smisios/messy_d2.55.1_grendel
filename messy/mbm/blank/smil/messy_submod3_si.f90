! **********************************************************************
MODULE messy_submod3_si
! **********************************************************************

  ! SMCL
  USE messy_main_tracer, ONLY: STRLEN_FNAME
  USE messy_submod3

  IMPLICIT NONE
  PRIVATE

  TYPE T_INIT_TRACER
     CHARACTER(LEN=STRLEN_FNAME) :: fullname
     REAL(DP)                    :: scale
     REAL(DP)                    :: decay
  END TYPE T_INIT_TRACER
  
  ! WORKSPACE
  INTEGER, PARAMETER                    :: NMAXI = 25
  TYPE(T_INIT_TRACER), DIMENSION(NMAXI) :: INIT_TRACER  ! CPL-NAMELIST
  INTEGER, DIMENSION(NMAXI)             :: IDX          ! tracer IDs

  PUBLIC :: submod3_initialize
  PUBLIC :: submod3_init_coupling
  PUBLIC :: submod3_init_tracer
  PUBLIC :: submod3_global_start
  PUBLIC :: submod3_physc
  !PRIVATE :: submod3_read_nml_cpl

CONTAINS

  ! --------------------------------------------------------------------
  SUBROUTINE submod3_initialize

    ! BMIL
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_bcast, p_io
    USE messy_main_blather_bi, ONLY: error_bi
    USE messy_main_tools,      ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'submod3_init_coupling'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status
    INTEGER                     :: i

    ! INITIALIZE CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL submod3_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF
    DO i=1, NMAXI
       CALL p_bcast(INIT_TRACER(i)%fullname, p_io)
       CALL p_bcast(INIT_TRACER(i)%scale, p_io)
       CALL p_bcast(INIT_TRACER(i)%decay, p_io)
    END DO

  END SUBROUTINE submod3_initialize
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  SUBROUTINE submod3_init_coupling

    ! BMIL
    USE messy_main_tracer_mem_bi,   ONLY: GPTRSTR
    USE messy_main_tracer_tools_bi, ONLY: tracer_halt
    ! SMCL
    USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM
    USE messy_main_tracer,        ONLY: get_tracer, full2base_sub

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'submod3_init_coupling'
    INTEGER :: status
    INTEGER :: i
    CHARACTER(LEN=STRLEN_MEDIUM) :: basename, subname

    IDX(:) = 0
    DO i=1, NMAXI
       
       CALL full2base_sub(status, INIT_TRACER(i)%fullname, basename, subname)
       CALL tracer_halt(substr, status)

       CALL get_tracer(status, GPTRSTR, basename, subname=subname  &
            , idx=IDX(i) )
       IF (status /= 0) IDX(i) = 0

    END DO

  END SUBROUTINE submod3_init_coupling
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  SUBROUTINE submod3_init_tracer

    ! BMIL
    USE messy_main_tracer_mem_bi,   ONLY: xt, GPTRSTR
    USE messy_main_tracer_tools_bi, ONLY: tracer_halt
    USE messy_main_timer,           ONLY: lresume
    ! SMCL
    USE messy_main_tracer,        ONLY: tracer_iniflag
    
    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'submod3_init_tracer'
    INTEGER                     :: i, jt
    INTEGER                     :: status

    IF (lresume) RETURN

    loop: DO i=1, NMAXI

       IF (IDX(i) == 0) CYCLE
       jt = IDX(i)

       CALL osm_fill(xt(:,:,jt,:), INIT_TRACER(i)%scale)

       CALL tracer_iniflag(status, GPTRSTR, jt, lset=.TRUE.)
       CALL tracer_halt(substr, status)

    END DO loop

  END SUBROUTINE submod3_init_tracer
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  SUBROUTINE submod3_global_start

    ! BMIL
    USE messy_main_tracer_mem_bi, ONLY: xt

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER         :: substr = 'submod3_global_start'
    INTEGER                             :: i, jt

    loop: DO i=1, NMAXI

       IF (IDX(i) == 0) CYCLE
       jt = IDX(i)

       CALL osm_diff(xt(:,:,jt,:))

    END DO loop

  END SUBROUTINE submod3_global_start
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  SUBROUTINE submod3_physc

    ! BMIL
    USE messy_main_tracer_mem_bi, ONLY: qxt

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER         :: substr = 'submod3_physc'
    INTEGER                             :: i, jt

    loop: DO i=1, NMAXI

       IF (IDX(i) == 0) CYCLE
       jt = IDX(i)

       CALL osm_decay(qxt(:,jt,:), INIT_TRACER(i)%decay)

    END DO loop

  END SUBROUTINE submod3_physc
  ! --------------------------------------------------------------------

! **********************************************************************
! PRIVATE INTERFACE ROUTINES
! **********************************************************************

! ----------------------------------------------------------------------
  SUBROUTINE submod3_read_nml_cpl(status, iou)

    ! SMCL
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'submod3_read_nml_cpl'

    NAMELIST /CPL/ INIT_TRACER

    ! LOCAL
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    !
    ! CHECK NAMELIST

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE submod3_read_nml_cpl
! ----------------------------------------------------------------------

! **********************************************************************
END MODULE messy_submod3_si
! **********************************************************************
