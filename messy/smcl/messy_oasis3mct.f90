! **********************************************************************
!
! SUBMODEL CORE LAYER (SMCL) ROUTINES FOR MESSy SUBMODEL OASIS3-MCT
!
! THIS SUBMODEL INCLUDES THE OASIS3-MCT COUPLER SOFTWARE INTO THE MESSy
! INTERFACE
!
!
! Author : Astrid Kerkweg, Meteorological Institute, University of Bonn
!           January 2018
!
! **********************************************************************

! **********************************************************************
MODULE messy_oasis3mct
! **********************************************************************

  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP, STRLEN_VLONG

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: DP

  ! ----------- <

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'oasis3mct'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '0.1'
 
  ! STRLEN_VLONG (MESSY) == ic_lvar (OASIS) == 80
  CHARACTER(LEN=STRLEN_VLONG), PUBLIC :: compname = ''

  ! PUBLIC SUBROUTINES (to be called from messy_bufly_e5.f90)
  !PUBLIC :: oasis3mct_read_nml_ctrl
  ! PRIVATE SUBROUTINES

CONTAINS



  ! =========================================================================
!!$  SUBROUTINE oasis3mct_read_nml_ctrl(status, iou)
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
!!$    NAMELIST /CTRL/ compname
!!$
!!$    ! LOCAL
!!$    CHARACTER(LEN=*), PARAMETER       :: substr='oasis3mct_read_nml_ctrl'
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
!!$  END SUBROUTINE oasis3mct_read_nml_ctrl
  ! =========================================================================

  ! =========================================================================
  ! ### add your own private subroutines here
  ! =========================================================================
  ! ===========================================================================
CHARACTER(LEN=256) FUNCTION MESSy_OASIS_ERROR(status)

  ! This subroutine provides an error string for a given status

  IMPLICIT NONE
     
  INTEGER, INTENT(IN) :: status
    
  SELECT CASE(status)
     !
     ! GRID ERRORS
        !
  CASE(100)
     messy_oasis_error = &
          TRIM('location string must contain 1-3 colon seperated entries')
  CASE(101)
     messy_oasis_error = TRIM('parse_location: subroutine names do not much')
  CASE(102)
     messy_oasis_error = TRIM('parse_location: second chunk must equal B or E')
!!$  CASE(200)
!!$     messy_oasis_error = &
!!$          TRIM('postproc string must contain colon seperated entries')
!!$  CASE()
!!$     messy_oasis_error = &
!!$          TRIM('')
!!$  CASE()
!!$     messy_oasis_error = &
!!$          TRIM('')
!!$  CASE()
!!$     messy_oasis_error = &
!!$          TRIM('')
!!$  CASE()
!!$     messy_oasis_error = &
!!$          TRIM('')

  END SELECT
     
END FUNCTION MESSY_OASIS_ERROR
 ! ===========================================================================

! **********************************************************************
END MODULE messy_oasis3mct
! **********************************************************************

