! **********************************************************************
!
! SUBMODEL CORE LAYER (SMCL) ROUTINES FOR MESSy SUBMODEL TPULSE
!
! THIS SUBMODEL IS USED TO ILLUSTRATE THE BUTTERFLY EFFECT
! VIA INITIAL SURFACE PRESSURE AND/OR TEMPERATURE PERTURBATIONS.
!
! MOREOVER, THIS SUBMODEL SERVES AS A TEMPLATE FOR NEW SUBMODELS.
!
! Author : Patrick Joeckel, DLR-IPA, October  2009
!
! References:
!
! * P. Joeckel, R. Sander, A. Kerkweg, H. Tost, and J. Lelieveld,
!   Technical Note: The Modular Earth Submodel System (MESSy) - a new
!   approach towards Earth System Modeling,
!   Atmos. Chem. Phys., 5, 433-444, 2005.
!   http://www.atmos-chem-phys.net/5/433 
! * Patrick Joeckel, Technical note: Recursive rediscretisation of
!   geo-scientific data in the Modular Earth Submodel System (MESSy),
!   Atmos. Chem. Phys., 6, 3557-3562, 2006.
!   http://www.atmos-chem-phys.net/6/3557
! * P. Joeckel, A. Kerkweg, J. Buchholz-Dietsch, H. Tost, R. Sander, and
!   A. Pozzer, Technical Note: Coupling of chemical processes with the
!   Modular Earth Submodel System (MESSy) submodel TRACER,
!   Atmos. Chem. Phys., 8, 1677-1687, 2008.
!   http://www.atmos-chem-phys.net/8/1677 
!
! **********************************************************************

! **********************************************************************
MODULE messy_tpulse
! **********************************************************************

  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: DP

  ! ----------- <

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'tpulse'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '1.1'
 
  ! PUBLIC SUBROUTINES (to be called from messy_tpulse_e5.f90)
  PUBLIC :: tpulse_read_nml_ctrl
  ! ### add your own public subroutines here
!  PUBLIC :: perturb

  ! PRIVATE SUBROUTINES
  ! ### add your own private subroutines here

CONTAINS

  ! =========================================================================
  ! ### add your own public subroutines here
  ! =========================================================================

  ! =========================================================================

  ! =========================================================================
  SUBROUTINE tpulse_read_nml_ctrl(status, iou)

    ! ------------------------------------------------------------------
    ! This routine is used to read the CTRL-namelist of the submodel.
    ! ------------------------------------------------------------------

    ! MESSy INTERFACE
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

!   NAMELIST /CTRL/
    ! LOCAL
 !   CHARACTER(LEN=*), PARAMETER       :: substr='tpulse_read_nml_ctrl'
 !   LOGICAL                           :: lex          ! file exists ?
 !   INTEGER                           :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR
    
 !   CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
 !   IF (.not.lex) RETURN    ! <modstr>.nml does not exist
!
 !   READ(iou, NML=CTRL, IOSTAT=fstat)
 !   CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
 !   IF (fstat /= 0) RETURN  ! error while reading namelist
    
    ! ### ADD HERE DIAGNOSTIC OUPUT FOR LOG-FILE
 
 !   CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR
    
  END SUBROUTINE tpulse_read_nml_ctrl
  ! =========================================================================

  ! =========================================================================
  ! ### add your own private subroutines here
  ! =========================================================================

! **********************************************************************
END MODULE messy_tpulse
! **********************************************************************

