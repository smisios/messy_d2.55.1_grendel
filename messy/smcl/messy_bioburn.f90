MODULE messy_bioburn

  !----------------------------------------------------------------------------
  !  bioburn : Model for Biomass burning Emissions
  !
  !  AUTHOR:  David Cabrera, MPICH, Sept 2013
  !----------------------------------------------------------------------------

  USE messy_main_constants_mem,  ONLY: DP, SP, STRLEN_MEDIUM

  IMPLICIT NONE
  PRIVATE

  ! GLOBAL PARAMETERS
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'bioburn'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '0.1'

  PUBLIC :: dp, sp, STRLEN_MEDIUM

  ! switch for output
  LOGICAL, PUBLIC :: l_verbose = .true.   ! GLOBAL SWITCH

  PUBLIC :: bioburn_read_nml_ctrl

CONTAINS

  !----------------------------------------------------------------------------

  SUBROUTINE bioburn_read_nml_ctrl(status, iou)

    !  bioburn MODULE ROUTINE
    !
    ! read namelist
    !
    ! Author: Pozzer Andrea, MPICH, Oct 2004

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /CTRL/ l_verbose

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='bioburn_read_nml_ctrl'
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    ! INITIALIZE NAMELIST VARIABLES
    l_verbose = .false.

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0 ! NO ERROR

  END SUBROUTINE bioburn_read_nml_ctrl

END MODULE
