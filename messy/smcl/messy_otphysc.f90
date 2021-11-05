MODULE messy_otphysc

! **********************************************************************
!
! Author : Andrea Pozzer, MPICH, September  2007
!
! **********************************************************************

  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP, SP, I4, I8, FLAGGED_BAD

  IMPLICIT NONE
  PRIVATE

  ! CTRL-NAMELIST PARAMETER
  REAL(DP), PUBLIC  :: MASK_VALUE
  PUBLIC :: DP, SP, I4, I8

  ! ----------- <

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'otphysc'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '2.1'

  PUBLIC :: otphysc_read_nml_ctrl

  CONTAINS
   ! ---------------------------------------------------------------------------
 
   SUBROUTINE otphysc_read_nml_ctrl(status, iou)
 
     !  MPIOM MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
     !
     ! read namelist for 'coupling' to ECHAM5
     !
     ! Author: Pozzer Andrea, MPICH, Oct 2004
 
 
     USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
 
     IMPLICIT NONE
 
     ! I/O
     INTEGER, INTENT(OUT) :: status      ! error status
     INTEGER, INTENT(IN)  :: iou         ! I/O unit
 
   NAMELIST /CTRL/ MASK_VALUE
 
     ! LOCAL
     CHARACTER(LEN=*), PARAMETER :: substr='otphysc_read_nml_ctrl'
     LOGICAL              :: lex      ! file exists ?
     INTEGER              :: fstat    ! file status
 
     status = 1
 
!----------------------------------------------------------------------
! DEFAULT PARAMETER SETTINGS - CAN BE OVERWRITEN BY THE NAMELIST
     MASK_VALUE = FLAGGED_BAD ! ==>  -1E34
!----------------------------------------------------------------------
     
     CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
     IF (.not.lex) RETURN    ! <modstr>.nml does not exist
 
     READ(iou, NML=CTRL, IOSTAT=fstat)

     CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
     IF (fstat /= 0) RETURN  ! error while reading namelist
 
     CALL read_nml_close(substr, iou, modstr)
 
     status = 0 ! NO ERROR
 
   END SUBROUTINE otphysc_read_nml_ctrl
 
   ! ---------------------------------------------------------------------------

! **********************************************************************
END MODULE messy_otphysc
! **********************************************************************
