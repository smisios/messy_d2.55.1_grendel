MODULE messy_clamscheme5
! **************************************************************************
! MODULE FOR CLaMS CHEMISTRY MODULE 'CHEM' for ECHAM gridpoints
! **************************************************************************

  USE messy_main_constants_mem, ONLY: DP

  IMPLICIT NONE

  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modstr = 'clamscheme5'
  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modver = '1.0'

  REAL(DP) :: zeta_max   ! Max. Zeta Level for interpolation    
  INTEGER :: bound_up    ! Max. Level for interpolation(level-number)
  CHARACTER(160) :: e5chem_initfile

!----------- 
CONTAINS
!----------- 

  SUBROUTINE clamscheme5_read_nml(status, iou)
    !Set default values and read namelist

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE
    
    !I/O
    INTEGER, INTENT(OUT) ::   status
    INTEGER, INTENT(IN)  ::   iou
    !LOCAL
    CHARACTER(LEN=*),PARAMETER :: substr='clamscheme5_read_nml'
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    NAMELIST /CTRL/ zeta_max, e5chem_initfile

    status = 1 !ERROR

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist
    
    CALL read_nml_close(substr, iou, modstr)
    
    status = 0 !NO ERROR
    
  END SUBROUTINE clamscheme5_read_nml


END MODULE messy_clamscheme5
