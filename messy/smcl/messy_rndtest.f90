MODULE MESSY_RNDTEST

  ! VERSION: see "modver" below
  ! AUTHOR(S):

  IMPLICIT NONE
  PRIVATE ! default for all
  SAVE

  INTEGER, PARAMETER                              :: NMAXDOM = 10
  INTEGER,  PUBLIC, DIMENSION(NMAXDOM)            :: NCELL
  INTEGER,  PUBLIC, DIMENSION(:), ALLOCATABLE     :: NCELL_LOC

  ! PUBLIC SUBROUTINES CALLED BY GCM/CTM
  ! ### SMIL: RNDTEST_INITIALIZE
  PUBLIC  :: RNDTEST_READ_NML_CTRL


  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'rndtest'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '1.0'

CONTAINS

! ###########################################################################
! ### PUBLIC SUBROUTINES
! ###########################################################################

! -------------------------------------------------------------------
  SUBROUTINE RNDTEST_READ_NML_CTRL(status, iou)

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    ! PURPOSE:
    !   READ RNDTEST NAMELIST
    !
    ! AUTHOR(S)
    !   AStrid Kerkweg, 2018/01, UNi Bonn

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    ! NAMELIST CTRL
    NAMELIST /CTRL/ NCELL

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER       :: substr = 'rndtest_read_nml_ctrl'
    LOGICAL                           :: lex          ! file exists ?
    INTEGER                           :: fstat        ! file status
    INTEGER                           :: i

    ! INITIALIZE
    status = 1 ! ERROR

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

 
    DO i = 1, SIZE(NCELL)
       WRITE(*,*) 'number of cells in block ',i, NCELL(i) 
    END DO

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE RNDTEST_READ_NML_CTRL
! ----------------------------------------------------------------------

END MODULE MESSY_RNDTEST
! ===================================================================
