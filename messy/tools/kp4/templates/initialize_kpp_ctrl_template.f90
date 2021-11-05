MODULE initialize_kpp_ctrl_template

  ! HEADER MODULE initialize_kpp_ctrl_template

  ! NOTES:
  ! - L_VECTOR is automatically defined by kp4
  ! - VL_DIM is automatically defined by kp4
  ! - I_LU_DI is automatically defined by kp4
  ! - WANTED is automatically defined by xmecca
  ! - ICNTRL RCNTRL are automatically defined by kpp
  ! - "USE messy_main_tools" is in Module_header of messy_mecca_kpp.f90
  ! - SAVE will be automatically added by kp4

  IMPLICIT NONE
  !SAVE

  ! FOR FIXED TIME STEP CONTROL
  ! ... max. number of fixed time steps (sum must be 1)
  INTEGER, PARAMETER         :: NMAXFIXSTEPS = 50
  ! ... switch for fixed time stepping
  LOGICAL, PUBLIC            :: l_fixed_step = .FALSE.
  INTEGER, PUBLIC            :: nfsteps = 1
  ! ... number of kpp control parameters
  INTEGER, PARAMETER, PUBLIC :: NKPPCTRL = 20
  !
  INTEGER,  DIMENSION(NKPPCTRL), PUBLIC     :: icntrl = 0
  REAL(DP), DIMENSION(NKPPCTRL), PUBLIC     :: rcntrl = 0.0_dp
  REAL(DP), DIMENSION(NMAXFIXSTEPS), PUBLIC :: t_steps = 0.0_dp

  ! END HEADER MODULE initialize_kpp_ctrl_template

CONTAINS

SUBROUTINE initialize_kpp_ctrl(status, iou, modstr)

  IMPLICIT NONE

  ! I/O
  INTEGER,          INTENT(OUT) :: status
  INTEGER,          INTENT(IN)  :: iou     ! logical I/O unit
  CHARACTER(LEN=*), INTENT(IN)  :: modstr  ! read <modstr>.nml

  ! LOCAL
  REAL(DP) :: tsum
  INTEGER  :: i

  CALL kpp_read_nml_ctrl(status, iou)
  IF (status /= 0) RETURN

  ! check fixed time steps
  tsum = 0.0_dp
  DO i=1, NMAXFIXSTEPS
     IF (t_steps(i) < TINY(0.0_DP)) EXIT
     tsum = tsum + t_steps(i)
  END DO

  nfsteps = i-1

  l_fixed_step = (nfsteps > 0) .AND. ( (tsum -1.0) < TINY(0.0_DP) )

  ! mz_pj_20070531+
  ! DIAGNOSTIC OUTPUT TO LOG-FILE
  WRITE(*,*) ' MECHANISM        : ',WANTED
  !
  IF (L_VECTOR) THEN
     WRITE(*,*) ' MODE             : VECTOR (LENGTH=',VL_DIM,')'
  ELSE
     WRITE(*,*) ' MODE             : SCALAR'
  END IF
  !
  WRITE(*,*) ' DE-INDEXING MODE :',I_LU_DI
  !
  WRITE(*,*) ' ICNTRL           : ',icntrl
  WRITE(*,*) ' RCNTRL           : ',rcntrl
  !
  ! NOTE: THIS IS ONLY MEANINGFUL FOR VECTORIZED (kp4) ROSENBROCK-METHODS
  IF (L_VECTOR) THEN
     IF (l_fixed_step) THEN
        WRITE(*,*) ' TIME STEPS       : FIXED (',t_steps(1:nfsteps),')'
     ELSE
        WRITE(*,*) ' TIME STEPS       : AUTOMATIC'
     END IF
  ELSE
     WRITE(*,*) ' TIME STEPS       : AUTOMATIC '//&
          &'(t_steps (CTRL_KPP) ignored in SCALAR MODE)'
  END IF
  ! mz_pj_20070531-

  status = 0

CONTAINS

  SUBROUTINE kpp_read_nml_ctrl(status, iou)

    ! READ MECCA NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! Author: Astrid Kerkweg,  MPICH, June 2007
    !         Patrick Joeckel, MPICH, June 2007

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    ! LOCAL
    LOGICAL :: lex   ! file exists?
    INTEGER :: fstat ! file status
    CHARACTER(LEN=*), PARAMETER :: substr = 'kpp_read_nml_ctrl'

    NAMELIST /CTRL_KPP/ icntrl, rcntrl, t_steps

    ! INITIALIZE
    status = 1 ! error

    ! INPUT NAMELIST
    CALL read_nml_open(lex, substr, iou, 'CTRL_KPP', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL_KPP, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL_KPP', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    WRITE(*,*) 'solver-specific method:      icntrl(3) = ', icntrl(3)
    WRITE(*,*) 'max. number of kpp-substeps: icntrl(4) = ', icntrl(4)

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! no error

  END SUBROUTINE kpp_read_nml_ctrl

END SUBROUTINE initialize_kpp_ctrl

SUBROUTINE error_output(C,ierr,PE)

  ! must be set in module header: USE messy_main_tools,      ONLY: str

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ierr
  INTEGER, INTENT(IN) :: PE
  REAL(dp), DIMENSION(:),INTENT(IN) :: C

  INTEGER,SAVE :: NUM =0
  INTEGER :: iou
  INTEGER :: i

  CHARACTER(LEN=250)  :: filename
  CHARACTER(LEN=1000) :: strnum
  CHARACTER(LEN=1000) :: strPE
  CHARACTER(256)      :: info

  LOGICAL             :: opened
  IF (ierr >= 0) RETURN

  NUM = NUM +1

  strnum=str(NUM)
  strPE=str(PE)

  WRITE(filename,*) 'mecca_PE'//TRIM(STRPE)//'_'//TRIM(STRNUM)//'.txt'

  iou = 0
  DO i=100,200
     INQUIRE(unit=i,opened=opened)
     IF (.NOT.opened) THEN
        iou = i
        EXIT
     END IF
  END DO

  IF (iou==0) THEN
     WRITE(info,'(a,i2.2,a,i2.2,a)') &
          'No unit in range < 100 : 200 > free.'
     RETURN
  ENDIF

  OPEN(IOU,FILE =TRIM(ADJUSTL(filename))&
       ,STATUS='NEW',ACTION= 'WRITE')

  ! ERROR STATUS
  WRITE(IOU,*) 'KPP ERRORSTATUS IERR:'
  WRITE(IOU,*) IERR

  ! SPECIES NAME
  WRITE(IOU,*) 'SELECTED MECHANISM'
  WRITE(IOU,*) WANTED
  WRITE(IOU,*)

  WRITE(IOU,*) 'NUMBER OF SPECIES:'
  WRITE(IOU,*) NSPEC
  WRITE(IOU,*)
  WRITE(IOU,*)  'NAMES OF SPECIES:'
  DO i=1,NSPEC
  WRITE(IOU,*)  SPC_NAMES(i)
  ENDDO
  WRITE(IOU,*)

  ! CONCENTRATIONS
  WRITE(IOU,*) 'concentrations (molec/cm^3(air) before KPP'
  DO i=1,NSPEC
     WRITE(IOU,*) C(i)
  ENDDO
  WRITE(IOU,*)

  ! rates
  WRITE(IOU,*) 'rate contants'
  DO i=1,NREACT
     WRITE(IOU,*) RCONST(i)
  ENDDO
  WRITE(IOU,*)
  CLOSE(IOU)

END SUBROUTINE error_output

END MODULE initialize_kpp_ctrl_template
