! ************************************************************************
MODULE messy_main_blather_bi
! ************************************************************************

  USE messy_main_blather

  IMPLICIT NONE

  PUBLIC :: main_blather_setup
  PUBLIC :: main_blather_finalize
  !
  PUBLIC :: start_message_bi    ! standard messages for start and end ...
  PUBLIC :: end_message_bi      ! .. of submodel-specific MESSy-routines
  PUBLIC :: info_bi
  PUBLIC :: warning_bi
  INTERFACE debug_bi
     MODULE PROCEDURE debug_s_bi
     MODULE PROCEDURE debug_r_bi
     MODULE PROCEDURE debug_i_bi
     MODULE PROCEDURE debug_r_3d2d3d_bi
  END INTERFACE debug_bi
  PUBLIC :: debug_bi
  PUBLIC :: error_bi
  PUBLIC :: messy_blather_endfile_bi

CONTAINS

  ! --------------------------------------------------------------------
  SUBROUTINE main_blather_setup

    USE messy_main_mpi_bi,           ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_tools,            ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_blather_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status

    ! DEFAULT FOR BASEMODELS, WHICH USE ENTIRE INFRASTRUCTURE
    ! need to be switched ON via $CTRL in blather.nml:
    ! DBG_SMS=1 : only SMs listed in sms-string
    ! DBG_ALL=2 : all SMs
    i_debug = DBG_OFF

    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL main_blather_read_nml_ctrl(status, iou)
       IF (status > 0) CALL error_bi( &
            'main_blather_read_nml_ctrl reported an error', substr)
    END IF
    CALL p_bcast(i_debug, p_io)
    CALL p_bcast(sms, p_io)

    CALL bl_init

  END SUBROUTINE main_blather_setup
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  SUBROUTINE main_blather_finalize

    IMPLICIT NONE

    CALL bl_finish

  END SUBROUTINE main_blather_finalize
  ! --------------------------------------------------------------------

  ! ####################################################################

  ! --------------------------------------------------------------------
  SUBROUTINE start_message_bi(modstr, str, substr)

    USE messy_main_mpi_bi,    ONLY: p_parallel_io

    IMPLICIT NONE

    ! I/O
    CHARACTER(LEN=*), INTENT(IN) :: modstr, str, substr

    CALL start_message(modstr, str, substr, p_parallel_io)

  END SUBROUTINE start_message_bi
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  SUBROUTINE end_message_bi(modstr, str, substr)

    USE messy_main_mpi_bi,    ONLY: p_parallel_io

    IMPLICIT NONE

    ! I/O
    CHARACTER(LEN=*), INTENT(IN) :: modstr, str, substr

    CALL end_message(modstr, str, substr, p_parallel_io)

  END SUBROUTINE end_message_bi
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  SUBROUTINE debug_s_bi(exmodstr, string, s, substr, l_print, iou)

    IMPLICIT NONE

    ! I/O
    CHARACTER(LEN=*), INTENT(IN)           :: exmodstr, string, substr
    CHARACTER(LEN=*), INTENT(IN)           :: s
    LOGICAL,          INTENT(IN), OPTIONAL :: l_print
    INTEGER,          INTENT(IN), OPTIONAL :: iou

    IF (i_debug == DBG_OFF) RETURN
    CALL debug(exmodstr, string, s, substr, l_print, iou)

  END SUBROUTINE debug_s_bi
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  SUBROUTINE debug_r_bi(exmodstr, string, r, substr, l_print, iou)

    USE messy_main_constants_mem, ONLY: dp

    IMPLICIT NONE

    ! I/O
    CHARACTER(LEN=*), INTENT(IN)           :: exmodstr, string, substr
    REAL(dp),         INTENT(IN)           :: r
    LOGICAL,          INTENT(IN), OPTIONAL :: l_print
    INTEGER,          INTENT(IN), OPTIONAL :: iou

    IF (i_debug == DBG_OFF) RETURN
    CALL debug(exmodstr, string, r, substr, l_print, iou)

  END SUBROUTINE debug_r_bi
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  SUBROUTINE debug_i_bi(exmodstr, string, i, substr, l_print, iou)

    IMPLICIT NONE

    ! I/O
    CHARACTER(LEN=*), INTENT(IN)           :: exmodstr, string, substr
    INTEGER,          INTENT(IN)           :: i
    LOGICAL,          INTENT(IN), OPTIONAL :: l_print
    INTEGER,          INTENT(IN), OPTIONAL :: iou

    IF (i_debug == DBG_OFF) RETURN
    CALL debug(exmodstr, string, i, substr, l_print, iou)

  END SUBROUTINE debug_i_bi
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  SUBROUTINE debug_r_3d2d3d_bi(exmodstr, string, substr &
       , var3d, var2d, yvar3d, l_print, iou)

    USE messy_main_constants_mem, ONLY: dp

    IMPLICIT NONE

    ! I/O
    CHARACTER(LEN=*),           INTENT(IN)           :: exmodstr, string, substr
    REAL(dp), DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: var3d
    REAL(dp), DIMENSION(:,:),   INTENT(IN), OPTIONAL :: var2d
    REAL(dp), DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: yvar3d
    LOGICAL,                    INTENT(IN), OPTIONAL :: l_print
    INTEGER,                    INTENT(IN), OPTIONAL :: iou

    IF (i_debug == DBG_OFF) RETURN
    CALL debug(exmodstr, string, substr, var3d, var2d, yvar3d, l_print, iou)

  END SUBROUTINE debug_r_3d2d3d_bi
  ! --------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE info_bi(string, substr, l_all)

    USE messy_main_mpi_bi, ONLY: p_parallel_io

    IMPLICIT NONE

    CHARACTER(LEN=*),           INTENT(IN) :: string
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: substr
    LOGICAL,          OPTIONAL, INTENT(IN) :: l_all

    LOGICAL :: zl_all

    IF (PRESENT(l_all)) THEN
       zl_all = l_all
    ELSE
       zl_all = p_parallel_io
    END IF

    CALL info(string, substr, l_print=zl_all)

  END SUBROUTINE info_bi
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
  SUBROUTINE warning_bi(string, substr)

    USE messy_main_mpi_bi, ONLY: p_parallel_io

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: string, substr

    CALL warning(string, substr, l_print=p_parallel_io)

  END SUBROUTINE warning_bi
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
  SUBROUTINE error_bi(string, substr)

    USE messy_main_mpi_bi,        ONLY: p_abort, p_parallel, p_pe
    USE messy_main_tools,         ONLY: find_next_free_unit
    USE messy_main_constants_mem, ONLY: iouerr

    IMPLICIT NONE
    INTRINSIC :: TRIM

    CHARACTER(LEN=*), INTENT(IN) :: string, substr
    INTEGER          :: iou
    LOGICAL          :: lex
    CHARACTER(LEN=7) :: efile ! e.g. "END0001"

    WRITE(iouerr,'(/,78("*"),/)')
    WRITE(iouerr,'("ERROR: p_pe   = ",I4)') p_pe
    WRITE(iouerr,'("ERROR: substr = ",A)')  TRIM(substr)
    WRITE(iouerr,'("ERROR: ",A)')           TRIM(string)
    WRITE(iouerr,'(/,78("*"),/)')

    ! write file 'END' for breaking rerun-chain in runscript xmessy
    iou = find_next_free_unit(100,300)
    WRITE(efile,'(a3,i4.4)') 'END',p_pe
    INQUIRE(FILE=efile, EXIST=lex) ! check if file exists
    IF (lex) THEN
      OPEN(iou, FILE=efile, STATUS='OLD', POSITION='APPEND') ! existing file
    ELSE
      OPEN(iou, FILE=efile, STATUS='UNKNOWN', POSITION='APPEND') ! new file
    ENDIF
    WRITE(iou,'("ERROR: p_pe   = ",I4)') p_pe
    WRITE(iou,'("ERROR: substr = ",A)')  TRIM(substr)
    WRITE(iou,'("ERROR: ",A)')           TRIM(string)
    CLOSE(iou)

#ifndef COSMO
    IF (p_parallel) THEN
       CALL p_abort
    ELSE
       CALL messy_blather_endfile_bi(string,substr)
       STOP
    ENDIF
#else
    IF (p_parallel) THEN
       CALL p_abort(string,substr)
    ELSE
       STOP
    ENDIF
#endif

  END SUBROUTINE error_bi
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE messy_blather_endfile_bi(string, psubstr)

    USE messy_main_mpi_bi, ONLY: p_pe
    USE messy_main_tools,  ONLY: find_next_free_unit
    USE messy_main_constants_mem, ONLY: STRLEN_ULONG

    IMPLICIT NONE

    INTRINSIC TRIM

    CHARACTER(LEN=*),           INTENT(IN) :: string
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: psubstr

    INTEGER          :: iou
    LOGICAL          :: lex
    CHARACTER(LEN=7) :: efile ! e.g. "END0001"
    CHARACTER(LEN=STRLEN_ULONG) :: substr = ''

    IF (PRESENT(psubstr)) THEN
       substr = TRIM(psubstr)
    ELSE
       substr = ''
    ENDIF

    CALL info_bi(string//" (endfile_bi) ", substr)

    ! write file 'END' for breaking rerun-chain in runscript xmessy
    iou = find_next_free_unit(100,300)
    WRITE(efile,'(a3,i4.4)') 'END',p_pe
    INQUIRE(FILE=efile, EXIST=lex) ! check if file exists
    IF (lex) THEN
      OPEN(iou, FILE=efile, STATUS='OLD', POSITION='APPEND') ! existing file
    ELSE
      OPEN(iou, FILE=efile, STATUS='UNKNOWN', POSITION='APPEND') ! new file
    ENDIF
    WRITE(iou,'("END: p_pe   = ",I4)') p_pe
    WRITE(iou,'("END: substr = ",A)')  TRIM(substr)
    WRITE(iou,'("END: ",A)')           TRIM(string)
    CLOSE(iou)

  END SUBROUTINE messy_blather_endfile_bi
!-------------------------------------------------------------------------

! ************************************************************************
END MODULE messy_main_blather_bi
! ************************************************************************
