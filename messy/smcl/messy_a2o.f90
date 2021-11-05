! ************************************************************************
MODULE messy_a2o
! ************************************************************************

  ! AUTHOR:
  !  Pozzer Andrea, MPICH, August 2007

  USE messy_main_constants_mem, ONLY: DP, SP, I4, I8

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: DP, SP, I4, I8

  ! ----------- <

  ! CTRL-NAMELIST PARAMETERS
  LOGICAL, PUBLIC :: l_verbose = .true.   ! GLOBAL SWITCH for scrip library 
                                          ! (verbose mode)

  ! GLOBAL PARAMETERS
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'a2o'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '2.0'

  PUBLIC :: check_max_min
  PUBLIC :: a2o_read_nml_ctrl

CONTAINS

  ! ---------------------------------------------------------------------------

  SUBROUTINE a2o_read_nml_ctrl(status, iou)

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

    NAMELIST /CTRL/ l_verbose      
                    

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='a2o_read_nml_ctrl'
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    ! INITIALIZE NAMELIST VARIABLES
    l_verbose = .true.

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0 ! NO ERROR

  END SUBROUTINE a2o_read_nml_ctrl

  ! ---------------------------------------------------------------------------

  SUBROUTINE check_max_min(zfld, zmax, zmin)

  IMPLICIT NONE 
  INTRINSIC SIZE

  REAL(DP), DIMENSION(:,:), INTENT (INOUT) :: zfld
  REAL(DP),                 INTENT(IN)     :: zmax
  REAL(DP),                 INTENT(IN)     :: zmin
  
  !PRIVATE
  INTEGER :: i,j

  DO i=1, SIZE(zfld(:,:),DIM=1)
     DO j=1, SIZE(zfld(:,:),DIM=2)
       IF (zfld(i,j).lt.zmin) zfld(i,j)=zmin
       IF (zfld(i,j).gt.zmax) zfld(i,j)=zmax
     ENDDO
  ENDDO
  
  END SUBROUTINE check_max_min
  

! ************************************************************************
END MODULE messy_a2o
! ************************************************************************
