MODULE MESSY_CAT

  USE messy_main_constants_mem, ONLY: dp

  IMPLICIT NONE
  
  SAVE

  CHARACTER(len=*), PARAMETER :: MODSTR='cat'
  CHARACTER(LEN=*), PARAMETER :: modver='1.0'

  ! switch for the number of the cloud scheme
  INTEGER, PUBLIC  :: cat_param
  REAL(dp), PUBLIC :: nhours, DIV_C
  LOGICAL, PUBLIC  :: USE_DDVSI =.FALSE.
  LOGICAL, PUBLIC  :: l_tracmix
  PUBLIC :: cat_read_nml_ctrl

CONTAINS

!-------------------------------------------------------------------------------

  SUBROUTINE cat_read_nml_ctrl(iou, status)

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

   ! I/O
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit
    INTEGER, INTENT(OUT) :: status ! error status

    NAMELIST /CTRL/ cat_param, nhours, DIV_C, USE_DDVSI, l_tracmix

   ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'cat_read_nml_ctrl'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    ! INITIALIZE GLOBAL CONTROL VARIABLES
    ! -> DEFAULT VALUES ARE SET AT DECLARATION ABOVE

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)

    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR


  END SUBROUTINE cat_read_nml_ctrl


!------------------------------------------------------------------------------
  
END MODULE MESSY_CAT
