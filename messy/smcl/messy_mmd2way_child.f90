! **************************************************************************
MODULE messy_mmd2way_child
! **************************************************************************

  ! MODULE FOR exchanging data with another model of coarse resolution
  !
  ! Author: Klaus Ketelsen, MPICH, Oct 2008
  !
  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP, STRLEN_ULONG
  USE messy_main_timer_event,   ONLY: time_event, io_time_event
  USE messy_mmd2way

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: DP

  ! ----------- <

  CHARACTER(LEN=*),PARAMETER,  PUBLIC :: submodstr = 'mmd2way_child'

  LOGICAL,             SAVE,   PUBLIC :: lvertsnrgd=.FALSE.
  ! PUBLIC ROUTINE
  PUBLIC :: mmd2way_child_read_nml_ctrl

  CONTAINS

  !===========================================================================
  
  SUBROUTINE mmd2way_child_read_nml_ctrl(status, iou)

    ! MMD2WAY MODULE ROUTINE (CORE)
    !
    ! READ MMD2WAY_CHILD CTRL NAMELIST, CHECK IT, 
    ! AND INITIALIZE GLOBAL VARIABLES
    !
    ! Author: Patrick Joeckel, MPICH, Feb 2002
    ! Modified: Astrid Kerkweg, UniMz, Dec 2012

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status
    INTEGER, INTENT(IN)  :: iou   ! logical I/O unit

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER       :: substr='mmd2way_read_nml_ctrl'
    LOGICAL                           :: lex          ! file exists ?
    INTEGER                           :: fstat        ! file status

    NAMELIST /CTRL/  lvertsnrgd
         

    status = 1 ! ERROR ON RETURN

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    CALL read_nml_close(substr, iou, modstr)

    status = 0

  END SUBROUTINE mmd2way_child_read_nml_ctrl

! **************************************************************************
END MODULE messy_mmd2way_child
! **************************************************************************
