#define _DIAGONALEX_
! *************************************************************************
MODULE messy_main_decomp
! *************************************************************************

  ! DECOMPOSITION MODULE CORE

  USE messy_main_constants_mem, ONLY: dp

  IMPLICIT NONE
  PUBLIC

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'decomp'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '0.1c'

  ! GLOBAL DECOMPOSITION OPTIONS
  !

  INTEGER, SAVE :: NLON = 1
  INTEGER, SAVE :: NLAT = 1

  INTEGER, SAVE :: NPX = 1
  INTEGER, SAVE :: NPY = 1

  ! for boundary exchange
#ifndef _DIAGONALEX_
  INTEGER :: my_cart_neigh(4)
  INTEGER :: my_cart_neigh_pole(4)
#else
  INTEGER :: my_cart_neigh(8)
  INTEGER :: my_cart_neigh_pole(8)
#endif
  INTEGER :: isendbuflen
  REAL(dp), ALLOCATABLE  ::  &
    sendbuf(:,:)       ! sending buffer for boundary exchange

  NAMELIST /CTRL/ NPX, NPY

  !PUBLIC :: main_decomp_read_nml_ctrl

CONTAINS

  ! -------------------------------------------------------------
  SUBROUTINE main_decomp_read_nml_ctrl(status, iou)

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit
    
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'messy_main_read_nml_ctrl'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status

    ! INITIALIZE
    status = 1  ! DEFAULT: ERROR

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! NO ERROR

  END SUBROUTINE main_decomp_read_nml_ctrl
  ! -------------------------------------------------------------



! *************************************************************************
END MODULE messy_main_decomp
! *************************************************************************
