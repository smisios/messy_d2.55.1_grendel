!*****************************************************************************
!                Time-stamp: <2014-02-10 14:42:29 steil>
!*****************************************************************************
MODULE messy_mtskip

  ! Module to rule them all.
  ! Routines to enable in messy the skipping of timesteps:
  ! mtskip = multi timestep skipping

  USE messy_main_constants_mem, ONLY: DP

  IMPLICIT NONE
  PRIVATE
  SAVE
  ! NAME OF SUBMODEL
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'mtskip'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '0.0' 
  
  !input-namelist
  INTEGER,  PUBLIC                    :: trigmtskip = 1 
  REAL(DP), PUBLIC                    :: addon_skip = -20._dp
  REAL(DP), PUBLIC                    :: phtl_offset = 0._dp
  INTEGER,  PUBLIC                    :: nstep_always = 10
  LOGICAL,  PUBLIC                    :: lfixneg = .FALSE.
  LOGICAL,  PUBLIC                    :: lbudneg = .FALSE.
  LOGICAL,  PUBLIC                    :: lmecca_mtskip = .FALSE.
  LOGICAL,  PUBLIC                    :: lgmxe_mtskip = .FALSE.
  !LOGICAL, PUBLIC                     :: lmsbm_mtskip = .FALSE.

  NAMELIST /CTRL/ trigmtskip, addon_skip, nstep_always, phtl_offset  &
                 ,lfixneg, lbudneg                                   &
                 ,lmecca_mtskip ,lgmxe_mtskip !, lmsbm_mtskip

  PUBLIC  :: DP
  PUBLIC  :: mtskip_read_nml_ctrl
  PUBLIC  :: mtskip_event_trig

  ! **************************************************************************
  CONTAINS
  ! **************************************************************************

    SUBROUTINE mtskip_read_nml_ctrl(status,iou)

      USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

      IMPLICIT NONE

      ! I/O
      INTEGER, INTENT(IN)  :: iou    ! logical I/O unit
      INTEGER, INTENT(OUT) :: status ! error status

      ! LOCAL
      CHARACTER(LEN=*), PARAMETER  :: substr = 'mtskip_read_nml_ctrl'
      LOGICAL                      :: lex          ! file exists?
      INTEGER                      :: fstat        ! file status

      status = 1            ! error

      
      CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
      IF (.not.lex) RETURN    ! <modstr>.nml does not exist

      READ(iou, NML=CTRL, IOSTAT=fstat)
      CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
      IF (fstat /= 0) RETURN  ! error while reading namelist

      ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES

      CALL read_nml_close(substr, iou, modstr)

      status = 0  ! no ERROR

    END SUBROUTINE mtskip_read_nml_ctrl

    ! ************************************************************************

    LOGICAL FUNCTION mtskip_event_trig(istep_trig, istep_modell)

      INTEGER, INTENT(IN) :: istep_trig, istep_modell

      mtskip_event_trig = .FALSE.
      IF ( MOD(istep_modell,istep_trig) .EQ. 0 ) mtskip_event_trig = .TRUE.  

    END FUNCTION mtskip_event_trig

END MODULE messy_mtskip
