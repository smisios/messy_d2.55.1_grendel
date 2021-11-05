MODULE MESSY_SATSIMS

! This is the core module for the Satellite Simulators
! At the moment only the ISCCP simulator v3.4 is implemented

  USE MESSY_MAIN_CONSTANTS_MEM,   ONLY: DP, STRLEN_ULONG

  IMPLICIT NONE
 
  PRIVATE
  CHARACTER(len=*), PUBLIC, PARAMETER :: MODSTR='satsims'
  CHARACTER(LEN=*), PUBLIC, PARAMETER :: modver='1.0'

  ! switch for the ISCCP simulator
  LOGICAL, PUBLIC :: l_isccp = .FALSE.
  CHARACTER(LEN=STRLEN_ULONG), PUBLIC :: taufile=''
  CHARACTER(LEN=STRLEN_ULONG), PUBLIC :: invfile=''
SAVE

  PUBLIC :: satsims_read_nml_ctrl, isccp_lookup_read

CONTAINS

!===============================================================================

  SUBROUTINE satsims_read_nml_ctrl(status,iou)

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

   ! I/O
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit
    INTEGER, INTENT(OUT) :: status ! error status

    NAMELIST /CTRL/ l_isccp

   ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'satsims_read_nml_ctrl'
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

  END SUBROUTINE satsims_read_nml_ctrl

!===============================================================================

  SUBROUTINE isccp_lookup_read(tautab, invtau, iunit, status)

    REAL(dp) :: tautab(0:255)      
    INTEGER  :: invtau(-20:45000) 
    INTEGER  :: status, iunit
    
    status = 1

    OPEN(iunit, FILE=TRIM(taufile), FORM='FORMATTED')
    READ(iunit,'(f30.20)') tautab
    CLOSE(iunit)
            
    OPEN(iunit, FILE=TRIM(invfile), FORM='FORMATTED')
    READ(iunit,'(i10)') invtau
    CLOSE(iunit)

    status = 0
    
  END SUBROUTINE isccp_lookup_read

!===============================================================================
END MODULE MESSY_SATSIMS
