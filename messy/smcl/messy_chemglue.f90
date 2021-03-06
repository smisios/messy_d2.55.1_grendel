!*****************************************************************************
!                Time-stamp: <2018-05-24 15:18:29 sander>
!*****************************************************************************

! submodel CHEMGLUE
! select one out of several MECCA mechansims
! Author: Rolf Sander, MPICH, Mainz, 2015-...

!*****************************************************************************

MODULE messy_chemglue

  USE messy_main_constants_mem, ONLY: DP
  IMPLICIT NONE
  PUBLIC
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'chemglue'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '1.0'
 
  INTEGER :: NMAXMECCA
  ! declare one variable for each batch file:
  ! (same as batch file name but "." and "-" replaced by "_")
  ! main:
  REAL(DP), SAVE :: ESMVal_butov_wet_bat = 0.0_dp
  REAL(DP), SAVE :: example_bat          = 0.0_dp
  REAL(DP), SAVE :: ff_bat               = 0.0_dp
  REAL(DP), SAVE :: iso_example_bat      = 0.0_dp
  REAL(DP), SAVE :: latex_bat            = 0.0_dp
  REAL(DP), SAVE :: mbl_bat              = 0.0_dp
  REAL(DP), SAVE :: mcfct_bat            = 0.0_dp
  REAL(DP), SAVE :: mom_bat              = 0.0_dp
  REAL(DP), SAVE :: mtchem_bat           = 0.0_dp
  REAL(DP), SAVE :: simple_bat           = 0.0_dp
  REAL(DP), SAVE :: simple_rxnrates_bat  = 0.0_dp
  REAL(DP), SAVE :: strato_bat           = 0.0_dp
  ! CASiMIR:
  REAL(DP), SAVE :: CASiMIR_05_bat       = 0.0_dp
  REAL(DP), SAVE :: CASiMIR_06_bat       = 0.0_dp
  REAL(DP), SAVE :: CASiMIR_07_bat       = 0.0_dp
  REAL(DP), SAVE :: CASiMIR_11_bat       = 0.0_dp
  ! CCMI:
  REAL(DP), SAVE :: CCMI_aero_02_bat     = 0.0_dp
  REAL(DP), SAVE :: CCMI_airtrac_02_bat  = 0.0_dp
  REAL(DP), SAVE :: CCMI_base_01_bat     = 0.0_dp
  REAL(DP), SAVE :: CCMI_base_01_tag_bat = 0.0_dp
  REAL(DP), SAVE :: CCMI_base_02_bat     = 0.0_dp
  REAL(DP), SAVE :: CCMI_base_02_polymeccatest_bat = 0.0_dp ! for testing only
  REAL(DP), SAVE :: CCMI_sens_01_bat     = 0.0_dp
  ! user_contributed:
  REAL(DP), SAVE :: chem_eval2_3_bat     = 0.0_dp
  REAL(DP), SAVE :: debug_bat            = 0.0_dp
  REAL(DP), SAVE :: e4chem_bat           = 0.0_dp
  REAL(DP), SAVE :: exb_example_bat      = 0.0_dp
  REAL(DP), SAVE :: lab_bat              = 0.0_dp
  REAL(DP), SAVE :: mbl_ocean_bat        = 0.0_dp
  REAL(DP), SAVE :: react4c_bat          = 0.0_dp
  REAL(DP), SAVE :: simple_MADE_bat      = 0.0_dp
  ! for CHEMGLUE and skeletal mechansims:
  REAL(DP), SAVE :: full_organic_bat     = 0.0_dp
  REAL(DP), SAVE :: simple_organic_bat   = 0.0_dp
  REAL(DP), SAVE :: skeleton_organic_bat = 0.0_dp
  REAL(DP), SAVE :: skeleton_lowterp_bat = 0.0_dp

  ! CTRL-NAMELIST PARAMETERS
  REAL(DP), PUBLIC :: ctrl_nml_dummy = 0.0_dp ! only dummy, currently not used

CONTAINS

  ! --------------------------------------------------------------------------

  SUBROUTINE assign_mecnum_names

    IMPLICIT NONE
    ! assign names to the numbers of the MECCA mechanisms:
    INCLUDE 'messy_chemglue.inc'

  END SUBROUTINE assign_mecnum_names

  ! --------------------------------------------------------------------------

  SUBROUTINE select_mechanism_from_pressure(mecnum, pressure)

    IMPLICIT NONE
    REAL(DP), INTENT(OUT) :: mecnum
    REAL(DP), INTENT(IN)  :: pressure

    ! as an example, the MECCA mechanism is based on the current pressure:
    IF (pressure > 2E4) THEN
      mecnum = CCMI_base_02_bat
    ELSE
      mecnum = CCMI_base_02_polymeccatest_bat
    ENDIF

  END SUBROUTINE select_mechanism_from_pressure

  ! --------------------------------------------------------------------------

  SUBROUTINE select_mechanism_from_mixrat(mecnum, y_C5H8, y_APINENE, y_TOLUENE)

    IMPLICIT NONE
    REAL(DP), INTENT(OUT) :: mecnum
    REAL(DP), INTENT(IN)  :: y_C5H8, y_APINENE, y_TOLUENE

    ! default is simple mechanism:
    mecnum = simple_organic_bat
    ! full mechanism when mixing ratios are high enough:
    IF (y_C5H8    > 1E-10) mecnum = full_organic_bat
    IF (y_APINENE > 1E-10) mecnum = full_organic_bat
    IF (y_TOLUENE > 1E-11) mecnum = full_organic_bat

  END SUBROUTINE select_mechanism_from_mixrat

  ! --------------------------------------------------------------------------

  SUBROUTINE select_mechanism_from_slm(mecnum, slm)

    IMPLICIT NONE
    REAL(DP), INTENT(OUT) :: mecnum
    REAL(DP), INTENT(IN)  :: slm

    ! as an example, the MECCA mechanism is based on the sea/land mask:
    IF (slm > 0.5) THEN
      mecnum = mom_bat
    ELSE
      mecnum = CCMI_base_02_bat
    ENDIF

  END SUBROUTINE select_mechanism_from_slm

  ! --------------------------------------------------------------------------

  SUBROUTINE select_mechanism_testing(mecnum, jp)

    IMPLICIT NONE
    REAL(DP), INTENT(OUT) :: mecnum
    INTEGER,  INTENT(IN)  :: jp

    ! for testing only:
    IF (MOD(jp,3)==0) THEN
      mecnum = CCMI_base_02_polymeccatest_bat
    ELSE
      mecnum = CCMI_base_02_bat
    ENDIF

  END SUBROUTINE select_mechanism_testing

  ! --------------------------------------------------------------------------

  SUBROUTINE chemglue_read_nml_ctrl(status, iou)

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit
    NAMELIST /CTRL/ ctrl_nml_dummy
    CHARACTER(LEN=*), PARAMETER       :: substr='chemglue_read_nml_ctrl'
    LOGICAL                           :: lex          ! file exists ?
    INTEGER                           :: fstat        ! file status

    status = 1 ! ERROR
    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN ! <modstr>.nml does not exist
    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN ! error while reading namelist
    WRITE(*,*) 'ctrl_nml_dummy = ', ctrl_nml_dummy
    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR
    
  END SUBROUTINE chemglue_read_nml_ctrl

!*****************************************************************************

END MODULE messy_chemglue

!*****************************************************************************

