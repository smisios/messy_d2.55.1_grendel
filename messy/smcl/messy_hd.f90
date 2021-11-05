! ************************************************************************
MODULE messy_hd
! ************************************************************************

  ! AUTHOR:
  !  Pozzer Andrea, MPICH, August 2007

  USE messy_main_constants_mem, ONLY: DP, SP, I4, I8

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: DP, SP, I4, I8

  CHARACTER(LEN=200),         PUBLIC :: inithd_files_path = '' ! CTRL namelsit
  LOGICAL, PUBLIC :: lhd_que     = .false.   ! GLOBAL SWITCH for diagnostic
  LOGICAL, PUBLIC :: lwater_corr = .false.   ! GLOBAL SWITCH for water correction (only coupled model)
  LOGICAL, PUBLIC :: lnew_glac   = .false.   ! mz_bk_20120928 ! GLOBAL SWITCH
  !                  ! for doing glacier calving as in ECHAM6 hydrology
  !                  ! false="old" ECHAM5 behavior / true="new" ECHAM6 version
  !                  ! default: ECHAM5 version
  ! ----------- <

  ! GLOBAL PARAMETERS
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'hd'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '2.1'
  ! NAMELIST

  PUBLIC :: hd_read_nml_ctrl

CONTAINS

  ! ---------------------------------------------------------------------------

  SUBROUTINE hd_read_nml_ctrl(status, iou)

    !  AIRSEA MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    !
    ! read namelist for 'coupling' to ECHAM5
    !
    ! Author: Pozzer Andrea, MPICH, Oct 2004


    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /CTRL/ inithd_files_path, lhd_que, lwater_corr                   &
         &        , lnew_glac ! mz_bk_20120928

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='hd_read_nml_ctrl'
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    ! INITIALIZE NAMELIST VARIABLES
    lhd_que     = .false.
    lwater_corr = .false.
    lnew_glac   = .false. ! mz_bk_20120928

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0 ! NO ERROR

  END SUBROUTINE hd_read_nml_ctrl

  ! ---------------------------------------------------------------------------
  

! ************************************************************************
END MODULE messy_hd
! ************************************************************************
