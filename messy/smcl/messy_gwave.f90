! ***********************************************************************
MODULE messy_gwave
  ! ***********************************************************************

  USE messy_main_constants_mem, ONLY: DP
 
  IMPLICIT NONE
  PRIVATE
  SAVE

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'gwave'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '1.0'


  ! GLOBAL CTRL NAMELIST VARIABLES
  INTEGER, PUBLIC :: gwparam = 1
  LOGICAL, PUBLIC :: tf_updates = .FALSE.

  INTEGER, PARAMETER    :: numgwparam =5
  CHARACTER(LEN=200), DIMENSION(numgwparam), PUBLIC :: &
                                      gwparamname = (/'Hines (orig. ECHAM5     ' &
                                                     ,'Medvedev-Klassen        ' &
                                                     ,'Hybrid-Lindzen-Matsuno  ' &
                                                     ,'Alexander and Dunkerton ' &
                                                     ,'Yigit and Medvedev      ' /)

  ! POINTERS FOR CHANNEL OBJECTS
  REAL(dp), DIMENSION(:,:,:), POINTER, PUBLIC :: gwdrag_u   => NULL(), & 
                                                 gwdrag_v   => NULL(), & 
                                                 gwflux_u   => NULL(), & 
                                                 gwflux_v   => NULL(), &
                                                 gwheat     => NULL(), &
                                                 gweddy     => NULL() 

  REAL(dp), DIMENSION(:,:,:), POINTER, PUBLIC :: gwdrag_u_ec   => NULL(), & 
                                                 gwdrag_v_ec   => NULL(), & 
                                                 gwflux_u_ec   => NULL(), & 
                                                 gwflux_v_ec   => NULL(), &
                                                 gwheat_ec     => NULL(), &
                                                 gweddy_ec    => NULL() 

  REAL(dp), DIMENSION(:,:,:), POINTER, PUBLIC :: gwdrag_u_cm   => NULL(), & 
                                                 gwdrag_v_cm   => NULL(), & 
                                                 gwflux_u_cm   => NULL(), & 
                                                 gwflux_v_cm   => NULL(), &
                                                 gwheat_cm     => NULL(), &
                                                 gweddy_cm    => NULL() 

  ! for basemodel
  REAL(dp), DIMENSION(:,:,:), POINTER, PUBLIC :: gwdrag_u_cg_gl   => NULL(), & 
                                                 gwdrag_v_cg_gl   => NULL(), &
                                                 gwdrag_u_cg   => NULL(), &
                                                 gwdrag_v_cg   => NULL()
  PUBLIC :: gwave_read_nml_ctrl

CONTAINS

 
  ! =========================================================================
  SUBROUTINE gwave_read_nml_ctrl(status, iou)

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
    

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    NAMELIST /CTRL/ gwparam, tf_updates

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER       :: substr='gwave_read_nml_ctrl'
    LOGICAL                           :: lex          ! file exists ?
    INTEGER                           :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist
    
    IF (gwparam > numgwparam) THEN
       status = 2 ! ERROR: gwparam out of range
       RETURN 
    ENDIF

    ! DIAGNOSE NAMELIST
    WRITE(*,*) 'Selected gravity wave parametrisation: ',TRIM(gwparamname(gwparam))
    CALL read_nml_close(substr, iou, modstr)
 
    status = 0 ! NO ERROR

  END SUBROUTINE gwave_read_nml_ctrl
  ! =========================================================================


  ! ***********************************************************************
END MODULE messy_gwave
! ***********************************************************************
