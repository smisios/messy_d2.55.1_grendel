MODULE MESSY_CLOUD

  USE MESSY_CLOUD_MEM,        ONLY: NCDNC, NICNC
  USE messy_main_constants_mem, ONLY: dp
  USE messy_main_tools,         ONLY: t_reset_par

  IMPLICIT NONE
 
  PRIVATE
  CHARACTER(len=*), PUBLIC, PARAMETER :: MODSTR='cloud'
  CHARACTER(LEN=*), PUBLIC, PARAMETER :: modver='2.2'
  
  ! switch for the number of the cloud scheme
  INTEGER, PUBLIC :: cloud_param
  ! switch for the cover scheme
  LOGICAL, PUBLIC :: lcover = .FALSE.    

  ! switch to reset the tuning parameters
  TYPE(t_reset_par), PUBLIC :: rset_ccraut = t_reset_par(.FALSE.,0._dp)
  TYPE(t_reset_par), PUBLIC :: rset_ccsaut = t_reset_par(.FALSE.,0._dp)

  ! switch to reset the cloud parameters
  TYPE(t_reset_par), PUBLIC :: rset_cauloc = t_reset_par(.FALSE.,0._dp)
  TYPE(t_reset_par), PUBLIC :: rset_csatsc = t_reset_par(.FALSE.,0._dp)
  TYPE(t_reset_par), PUBLIC :: rset_crhsc = t_reset_par(.FALSE.,0._dp)
  ! switch for gridpoint clouds
  LOGICAL, PUBLIC :: lcloud_gp = .TRUE.
  ! switch for lagrangian clouds
  LOGICAL, PUBLIC :: lcloud_lg = .FALSE.

  SAVE

  PUBLIC :: cloud_read_nml_ctrl

CONTAINS

!===============================================================================

  SUBROUTINE cloud_read_nml_ctrl(iou, status)

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

   ! I/O
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit
    INTEGER, INTENT(OUT) :: status ! error status

    NAMELIST /CTRL/ cloud_param, ncdnc, nicnc, lcover &
         , rset_ccraut, rset_ccsaut &
         , rset_cauloc, rset_csatsc, rset_crhsc

   ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'cloud_read_nml_ctrl'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    ! INITIALIZE GLOBAL CONTROL VARIABLES
    ! -> DEFAULT VALUES ARE SET AT DECLARATION ABOVE

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)

    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

    SELECT CASE (cloud_param)
    CASE(1)
      ncdnc = 0
      nicnc = 0
    CASE(2)
      nicnc = 0     
    END SELECT


  END SUBROUTINE cloud_read_nml_ctrl

!===============================================================================


END MODULE MESSY_CLOUD
