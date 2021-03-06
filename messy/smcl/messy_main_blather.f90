! ************************************************************************
MODULE messy_main_blather
! ************************************************************************

  USE messy_main_constants_mem, ONLY: STRLEN_XLONG, STRLEN_MEDIUM

  IMPLICIT NONE
  INTRINSIC :: NULL
  PRIVATE
  SAVE

  ! GLOBAL PARAMETERS
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'blather'
  CHARACTER(len=*), PARAMETER, PUBLIC :: modver = '2.0'

  INTEGER, PARAMETER, PUBLIC :: DBG_OFF = 0 ! OFF
  INTEGER, PARAMETER, PUBLIC :: DBG_SMS = 1 ! only SMs from &CTRL
  INTEGER, PARAMETER, PUBLIC :: DBG_ALL = 2 ! ON, all SMs

  ! CTRL namelist parameters
  ! set default for MBMs, which do NOT use the BMIL & namelist read
  INTEGER, PUBLIC                     :: i_debug = DBG_ALL
  CHARACTER(LEN=STRLEN_XLONG), PUBLIC :: sms = ''

  ! GLOBAL WORKSPACE
  INTEGER :: nsms = 0
  CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(:), POINTER :: lsms => NULL()

  PUBLIC :: start_message       ! standard messages for start and end ...
  PUBLIC :: end_message         ! .. of submodel-specific MESSy-routines
  PUBLIC :: info
  PUBLIC :: warning
  INTERFACE debug
     MODULE PROCEDURE debug_s
     MODULE PROCEDURE debug_r
     MODULE PROCEDURE debug_i
     MODULE PROCEDURE debug_r_3d2d3d
  END INTERFACE debug
  PUBLIC :: debug
  PUBLIC :: main_blather_read_nml_ctrl
  PUBLIC :: bl_init
  PUBLIC :: bl_finish
  PUBLIC :: bl_check

CONTAINS

  ! --------------------------------------------------------------------
  SUBROUTINE start_message(exmodstr, str, substr, l_print)

    USE messy_main_constants_mem, ONLY: HLINE1
    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! I/O
    CHARACTER(LEN=*), INTENT(IN)  :: exmodstr, str, substr
    LOGICAL, INTENT(IN), OPTIONAL :: l_print

    ! LOCAL
    LOGICAL :: zl_print

    IF (PRESENT(l_print)) THEN
       zl_print = l_print
    ELSE
       zl_print = .TRUE.  ! DEFAULT
    END IF

    IF (.NOT. zl_print) RETURN

    WRITE(*,*) HLINE1
    WRITE(*,*) '*** START ',TRIM(exmodstr),': ',TRIM(str), &
         ' (',TRIM(substr),')'

  END SUBROUTINE start_message
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  SUBROUTINE end_message(exmodstr, str, substr, l_print)

    USE messy_main_constants_mem, ONLY: HLINE1
    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! I/O
    CHARACTER(LEN=*), INTENT(IN)  :: exmodstr, str, substr
    LOGICAL, INTENT(IN), OPTIONAL :: l_print

    ! LOCAL
    LOGICAL :: zl_print

    IF (PRESENT(l_print)) THEN
       zl_print = l_print
    ELSE
       zl_print = .TRUE.  ! DEFAULT
    END IF

    IF (.NOT. zl_print) RETURN

    WRITE(*,*) '*** END   ',TRIM(exmodstr),': ',TRIM(str), &
         ' (',TRIM(substr),')'
    WRITE(*,*) HLINE1

  END SUBROUTINE end_message
  ! --------------------------------------------------------------------

  ! ---------------------------------------------------------------------
  SUBROUTINE info(string, substr, l_print)

    IMPLICIT NONE
    INTRINSIC :: PRESENT, TRIM
    CHARACTER(LEN=*),           INTENT(IN) :: string
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: substr
    LOGICAL,          OPTIONAL, INTENT(IN) :: l_print

    ! LOCAL
    LOGICAL :: zl_print

    IF (PRESENT(l_print)) THEN
       zl_print = l_print
    ELSE
       zl_print = .TRUE.  ! DEFAULT
    END IF

    IF (.NOT. zl_print) RETURN

    IF (PRESENT(substr)) THEN
       WRITE (*,'(A," (",A,")")') TRIM(string), TRIM(substr)
    ELSE
       WRITE (*,'(A)') TRIM(string)
    ENDIF

  END SUBROUTINE info
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE warning(string, substr, l_print)

    IMPLICIT NONE
    INTRINSIC :: TRIM
    CHARACTER(LEN=*),           INTENT(IN) :: string, substr
    LOGICAL,          OPTIONAL, INTENT(IN) :: l_print

    ! LOCAL
    LOGICAL :: zl_print

    IF (PRESENT(l_print)) THEN
       zl_print = l_print
    ELSE
       zl_print = .TRUE.  ! DEFAULT
    END IF

    IF (.NOT. zl_print) RETURN

    WRITE (*,'("WARNING: ",A," (",A,")")') TRIM(string), TRIM(substr)

  END SUBROUTINE warning
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE debug_s(exmodstr, string, s, substr, l_print, iou)

    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! I/O
    CHARACTER(LEN=*), INTENT(IN)           :: exmodstr, string, substr
    CHARACTER(LEN=*), INTENT(IN)           :: s
    LOGICAL,          INTENT(IN), OPTIONAL :: l_print
    INTEGER,          INTENT(IN), OPTIONAL :: iou

    ! LOCAL
    LOGICAL :: zl_print
    INTEGER :: ziou

    IF (PRESENT(l_print)) THEN
       zl_print = l_print
    ELSE
       zl_print = .TRUE.  ! DEFAULT
    END IF

    IF (.NOT. zl_print) RETURN
    IF (.NOT. bl_check(exmodstr)) RETURN

    IF (PRESENT(iou)) THEN
       ziou = iou
    ELSE
       ziou = 0  ! DEFAULT
    END IF

    WRITE (ziou,'("DEBUG: ",A,": ",A," (",A,")")') &
         TRIM(string), TRIM(s), TRIM(substr)

  END SUBROUTINE debug_s
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE debug_r(exmodstr, string, r, substr, l_print, iou)

    USE messy_main_constants_mem, ONLY: dp

    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! I/O
    CHARACTER(LEN=*), INTENT(IN)           :: exmodstr, string, substr
    REAL(dp),         INTENT(IN)           :: r
    LOGICAL,          INTENT(IN), OPTIONAL :: l_print
    INTEGER,          INTENT(IN), OPTIONAL :: iou

    ! LOCAL
    LOGICAL :: zl_print
    INTEGER :: ziou

    IF (PRESENT(l_print)) THEN
       zl_print = l_print
    ELSE
       zl_print = .TRUE.  ! DEFAULT
    END IF

    IF (.NOT. zl_print) RETURN
    IF (.NOT. bl_check(exmodstr)) RETURN

    IF (PRESENT(iou)) THEN
       ziou = iou
    ELSE
       ziou = 0  ! DEFAULT
    END IF

    WRITE (ziou,*) 'DEBUG: ',TRIM(string),': ', r, '(', TRIM(substr), ')'

  END SUBROUTINE debug_r
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE debug_i(exmodstr, string, i, substr, l_print, iou)

    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! I/O
    CHARACTER(LEN=*), INTENT(IN)           :: exmodstr, string, substr
    INTEGER,          INTENT(IN)           :: i
    LOGICAL,          INTENT(IN), OPTIONAL :: l_print
    INTEGER,          INTENT(IN), OPTIONAL :: iou

    ! LOCAL
    LOGICAL :: zl_print
    INTEGER :: ziou

    IF (PRESENT(l_print)) THEN
       zl_print = l_print
    ELSE
       zl_print = .TRUE.  ! DEFAULT
    END IF

    IF (.NOT. zl_print) RETURN
    IF (.NOT. bl_check(exmodstr)) RETURN

    IF (PRESENT(iou)) THEN
       ziou = iou
    ELSE
       ziou = 0  ! DEFAULT
    END IF

    WRITE (ziou,*) 'DEBUG: ',TRIM(string),': ', i, '(', TRIM(substr), ')'

  END SUBROUTINE debug_i
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE debug_r_3d2d3d(exmodstr, string, substr &
       , var3d, var2d, yvar3d, l_print, iou)

    USE messy_main_constants_mem, ONLY: dp

    IMPLICIT NONE
    INTRINSIC :: MINVAL, MAXVAL, SUM

    ! I/O
    CHARACTER(LEN=*),           INTENT(IN)           :: exmodstr, string, substr
    REAL(dp), DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: var3d
    REAL(dp), DIMENSION(:,:),   INTENT(IN), OPTIONAL :: var2d
    REAL(dp), DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: yvar3d
    LOGICAL,                    INTENT(IN), OPTIONAL :: l_print
    INTEGER,                    INTENT(IN), OPTIONAL :: iou

    ! LOCAL
    LOGICAL :: zl_print
    INTEGER :: ziou

    IF (PRESENT(l_print)) THEN
       zl_print = l_print
    ELSE
       zl_print = .TRUE.  ! DEFAULT
    END IF

    IF (.NOT. zl_print) RETURN
    IF (.NOT. bl_check(exmodstr)) RETURN

    IF (PRESENT(iou)) THEN
       ziou = iou
    ELSE
       ziou = 0  ! DEFAULT
    END IF

    IF (PRESENT(var3d)) THEN
       WRITE (ziou,*) 'DEBUG: ',TRIM(string),'(', TRIM(substr), '): '  &
            , 'MIN: ',MINVAL(var3d), 'MAX ',MAXVAL(var3d)              &
            , 'SUM ', SUM(var3d)
    END IF

    IF (PRESENT(yvar3d)) THEN
       WRITE (ziou,*) 'DEBUG: ',TRIM(string),'(', TRIM(substr), '): '  &
            , 'MIN: ',MINVAL(yvar3d), 'MAX ',MAXVAL(yvar3d)            &
            , 'SUM ', SUM(yvar3d)
    END IF

    IF (PRESENT(var2d)) THEN
       WRITE (ziou,*) 'DEBUG: ',TRIM(string),'(', TRIM(substr), '): '  &
            , 'MIN: ',MINVAL(var2d), 'MAX ',MAXVAL(var2d)              &
            , 'SUM ', SUM(var2d)
    END IF

  END SUBROUTINE debug_r_3d2d3d
  !-------------------------------------------------------------------------

  ! ####################################################################

  !-------------------------------------------------------------------------
  SUBROUTINE main_blather_read_nml_ctrl(status, iou)

    ! READ NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    NAMELIST /CTRL/ i_debug, sms

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_blather_read_nml_ctrl'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.NOT.lex) THEN
       status = -1 ! special case only for BLATHER
       CALL end_message(TRIM(modstr), 'INITIALISATION', substr)
       RETURN      ! <modstr>.nml does not exist
    END IF

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE main_blather_read_nml_ctrl
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE bl_init

    USE messy_main_tools, ONLY: strcrack

    IMPLICIT NONE

    IF (i_debug == DBG_OFF) RETURN
    IF (i_debug == DBG_SMS) CALL strcrack(sms, ';', lsms, nsms)

  END SUBROUTINE bl_init
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE bl_finish

    IMPLICIT NONE

    IF (i_debug == DBG_OFF) RETURN

    IF (ASSOCIATED(lsms)) THEN
       DEALLOCATE(lsms)
       NULLIFY(lsms)
    END IF

  END SUBROUTINE bl_finish
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  LOGICAL FUNCTION bl_check(exmodstr)

    IMPLICIT NONE
    INTRINSIC :: TRIM, ADJUSTL

    ! I/O
    CHARACTER(LEN=*), INTENT(IN) :: exmodstr

    ! LOCAL
    INTEGER :: i

    bl_check = (i_debug > DBG_OFF)
    IF ((i_debug == DBG_OFF) .OR. (i_debug == DBG_ALL)) RETURN

    ! ONLY, IF i_debug == DBG_SMS
    bl_check = .FALSE.
    DO i=1, nsms
       IF ( TRIM(ADJUSTL(exmodstr)) == TRIM(ADJUSTL(lsms(i))) ) THEN
          bl_check = .TRUE.
          EXIT
       END IF
    END DO

  END FUNCTION bl_check
  ! -------------------------------------------------------------------

! ************************************************************************
END MODULE messy_main_blather
! ************************************************************************
