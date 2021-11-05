! **********************************************************************
!
! SUBMODEL CORE LAYER (SMCL) ROUTINES FOR MESSy SUBMODEL TBUDGET
!
! THIS SUBMODEL IS USED FOR CALCULATING CHEMICAL BUDGETS
!
! Author : Phoebe Graf, DLR-IPA, April 2013
!
! References:
!
! * P. Jöckel, R. Sander, A. Kerkweg, H. Tost, and J. Lelieveld,
!   Technical Note: The Modular Earth Submodel System (MESSy) - a new
!   approach towards Earth System Modeling,
!   Atmos. Chem. Phys., 5, 433-444, 2005.
!   http://www.atmos-chem-phys.net/5/433 
! * P. Jöckel, A. Kerkweg, J. Buchholz-Dietsch, H. Tost, R. Sander, and
!   A. Pozzer, Technical Note: Coupling of chemical processes with the
!   Modular Earth Submodel System (MESSy) submodel TRACER,
!   Atmos. Chem. Phys., 8, 1677-1687, 2008.
!   http://www.atmos-chem-phys.net/8/1677 
! * Jöckel, P., Kerkweg, A., Pozzer, A., Sander, R., Tost, H., Riede, 
!   H., Baumgaertner, A., Gromov, S., & Kern, B.: Development cycle 2 of 
!   the Modular Earth Submodel System (MESSy2), Geoscientific Model 
!   Development, 3, 717-752, doi: 10.5194/gmd-3-717-2010, 
!   URL http://www.geosci-model-dev.net/3/717/2010/ (2010) 
!
! **********************************************************************

! **********************************************************************
MODULE messy_tbudget
! **********************************************************************

  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP

  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC :: DP

  ! ----------- <

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'tbudget'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '1.1'
 
  ! CTRL-NAMELIST PARAMETERS

  ! PUBLIC SUBROUTINES (to be called from messy_tbudget_e5.f90)
!!$  PUBLIC :: tbudget_read_nml_ctrl
  PUBLIC :: dgl_step

  ! PRIVATE SUBROUTINES
  ! ### add your own private subroutines here

CONTAINS

  ! =========================================================================
  ! ### add your own public subroutines here
  ! =========================================================================

  ! =========================================================================
  ELEMENTAL SUBROUTINE dgl_step(dmrdt, prod, loss, mr, tot)

    IMPLICIT NONE
    INTRINSIC :: ABS, EPSILON

    ! I/O
    REAL(dp), INTENT(OUT) :: dmrdt  ! tendency of tracer [mol/mol/s]
    REAL(dp), INTENT(in)  :: prod   ! production rate    [mol/mol/s]
    REAL(dp), INTENT(in)  :: loss   ! loss rate          [mol/mol/s] (< 0 !!!)
    REAL(dp), INTENT(in)  :: mr     ! mixing ratio       [mol/mol]
    REAL(dp), INTENT(in)  :: tot    ! total abundance    [mol/mol]

    IF (ABS(tot) > EPSILON(0.0_dp)) THEN
       dmrdt = prod + loss * (mr / tot)     ! loss must be negative !!!
    ELSE
       dmrdt = 0.0_dp
    END IF

  END SUBROUTINE dgl_step
  ! =========================================================================

!!$  ! =========================================================================
!!$  SUBROUTINE tbudget_read_nml_ctrl(status, iou)
!!$
!!$    ! ------------------------------------------------------------------
!!$    ! This routine is used to read the CTRL-namelist of the submodel.
!!$    ! ------------------------------------------------------------------
!!$
!!$    ! MESSy INTERFACE
!!$    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
!!$
!!$    IMPLICIT NONE
!!$
!!$    ! I/O
!!$    INTEGER, INTENT(OUT) :: status ! error status
!!$    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit
!!$
!!$    NAMELIST /CTRL/ 
!!$
!!$    ! LOCAL
!!$    CHARACTER(LEN=*), PARAMETER       :: substr='tbudget_read_nml_ctrl'
!!$    LOGICAL                           :: lex          ! file exists ?
!!$    INTEGER                           :: fstat        ! file status
!!$
!!$    ! INITIALIZE
!!$    status = 1 ! ERROR
!!$    
!!$    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
!!$    IF (.not.lex) RETURN    ! <modstr>.nml does not exist
!!$
!!$    READ(iou, NML=CTRL, IOSTAT=fstat)
!!$    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
!!$    IF (fstat /= 0) RETURN  ! error while reading namelist
!!$    
!!$    CALL read_nml_close(substr, iou, modstr)
!!$    status = 0 ! NO ERROR
!!$    
!!$  END SUBROUTINE tbudget_read_nml_ctrl
!!$  ! =========================================================================

  ! =========================================================================
  ! ### add your own private subroutines here
  ! =========================================================================

! **********************************************************************
END MODULE messy_tbudget
! **********************************************************************

