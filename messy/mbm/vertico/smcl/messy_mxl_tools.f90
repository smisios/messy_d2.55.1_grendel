
MODULE messy_mxl_tools

  !-------------------------------------------------------------------------------------------------------
  !  mxl : MiXed Layer model for the diurnal dynamics of the convective boundary layer  
  !
  !  AUTHOR:  Ruud Janssen, MPIC, Sept 2013
  !-------------------------------------------------------------------------------------------------------


  ! MESSy
  USE messy_main_constants_mem,  ONLY: DP, SP, STRLEN_MEDIUM

  IMPLICIT NONE 
  PRIVATE

!  PUBLIC :: mxl_read_nml_ctrl

CONTAINS

! op_pj_20150206+
  ! some compilers complain about an empty CONTAINS-block ...
  SUBROUTINE dummy
  END SUBROUTINE dummy
! op_pj_20150206-

!  !----------------------------------------------------------------------------
!
!  SUBROUTINE mxl_read_nml_ctrl(status, iou)
!
!    ! copied from:  bioburn MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
!    !
!    ! read namelist for 'coupling' to ECHAM5
!    !
!    ! Author: Pozzer Andrea, MPICH, Oct 2004
!
!
!    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
!
!    IMPLICIT NONE
!
!    ! I/O
!    INTEGER, INTENT(OUT) :: status     ! error status
!    INTEGER, INTENT(IN)  :: iou        ! I/O unit
!
!    NAMELIST /CTRL/ lat, lon
!
!    ! LOCAL
!    CHARACTER(LEN=*), PARAMETER :: substr='mxl_read_nml_ctrl'
!    LOGICAL              :: lex      ! file exists ?
!    INTEGER              :: fstat    ! file status
!
!    status = 1
!
!    ! INITIALIZE NAMELIST VARIABLES
!    l_verbose = .false.
!
!    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
!    IF (.not.lex) RETURN    ! <modstr>.nml does not exist
!
!    READ(iou, NML=CTRL, IOSTAT=fstat)
!    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
!    IF (fstat /= 0) RETURN  ! error while reading namelist
!    
!    CALL read_nml_close(substr, iou, modstr)
!
!    status = 0 ! NO ERROR
!
!  END SUBROUTINE mxl_read_nml_ctrl


END MODULE
