PROGRAM RANDOM

  USE messy_main_constants_mem, ONLY: DP
  USE messy_main_rnd,           ONLY: RND_F90, RND_MTW, RND_LUX &
                                    , RND_F90_GAUSS, RND_MTW_GAUSS &
                                    , RND_LUX_GAUSS &
                                    , rnd_init, rnd_number
  USE messy_main_tools,         ONLY: PTR_1D_ARRAY

  IMPLICIT NONE

  ! GLOBAL PARAMETERS
  CHARACTER(LEN=*), PARAMETER :: modstr = 'random'
  INTEGER, PARAMETER :: iou  = 17
  INTEGER, PARAMETER :: idat = 18

  INTEGER, PARAMETER :: NRPROC = 4
  INTEGER, DIMENSION(NRPROC) :: id

  ! NAMELIST PARAMETERS
  INTEGER, SAVE :: NCELL
  INTEGER, SAVE :: LOOPS
  INTEGER, DIMENSION(NRPROC) :: START_SEED = (/ 211169, 7249, 3372, 159941 /)
  INTEGER, PARAMETER, DIMENSION(NRPROC) :: RAND_MULT = (/ 1, 1, 1, 3 /)
  INTEGER, SAVE :: I_RANDOM_METHOD

  ! WORKSPACE
  INTEGER                               :: status
  TYPE(PTR_1D_ARRAY), DIMENSION(NRPROC) :: HARVEST
  INTEGER                               :: i, j

  ! INIT
  NCELL = 634880
  LOOPS = 10000

  CALL RANDOM_READ_NML_CPL(status, iou)
  IF (status /= 0) STOP 'ERROR in namelist'

  ! INIT MEMORY
  DO i=1, NRPROC
     ALLOCATE(HARVEST(i)%ptr(NCELL*RAND_MULT(i)))
  END DO

  ! INIT RANDOM GENERATORS
  DO i=1, NRPROC
     CALL rnd_init(status, id(i), I_RANDOM_METHOD, START_SEED(i))
     IF (status /= 0) THEN
        WRITE(*,*) 'ERROR in rnd_init (status = ',status,')'
        STOP
     END IF
  END DO

  OPEN(unit=idat, file="state.dat", status='unknown')

  DO j=1, LOOPS
     !WRITE(idat,*) RAND_STATE(:,1)
     DO i=1, NRPROC
        CALL rnd_number(id(i), HARVEST(i)%ptr(:))
     END DO
     WRITE(idat, *) HARVEST(1)%ptr(1),HARVEST(2)%ptr(1) &
          ,HARVEST(3)%ptr(1),HARVEST(4)%ptr(1)
  END DO  

  ! CLEAN MEMORY
  DO i=1, NRPROC
     DEALLOCATE(HARVEST(i)%ptr)
  END DO

  CLOSE(idat)

CONTAINS

! -------------------------------------------------------------------
  SUBROUTINE RANDOM_READ_NML_CPL(status, iou)

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    ! PURPOSE:
    !   READ ATTILA NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! AUTHOR(S)
    !   Michael Traub, MPICH, July 2003

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    ! NAMELIST CPL
    NAMELIST /CPL_BOX/ START_SEED, I_RANDOM_METHOD, NCELL, LOOPS

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER       :: substr = 'random_read_nml_cpl'
    LOGICAL                           :: lex          ! file exists ?
    INTEGER                           :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    CALL read_nml_open(lex, substr, iou, 'CPL_BOX', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL_BOX, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL_BOX', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    WRITE(*,*) 'GLOBAL SETTINGS'
    WRITE(*,*) '---------------'
    !
    SELECT CASE(I_RANDOM_METHOD)
    CASE(RND_F90)
       WRITE(*,*) 'RANDOM NUMBERS: F90-INTRINSIC (MACHINE DEPENDENT !)'

    CASE(RND_F90_GAUSS)
       WRITE(*,*) 'RANDOM NUMBERS: F90-INTRINSIC [GAUSS] (MACHINE DEPENDENT !)'

    CASE(RND_MTW)
       WRITE(*,*) 'RANDOM NUMBERS: MERSENNE TWISTER (MACHINE INDEPENDENT !)'

    CASE(RND_MTW_GAUSS)
       WRITE(*,*) 'RANDOM NUMBERS: MERSENNE TWISTER [GAUSS] (MACHINE INDEPENDENT !)'

    CASE(RND_LUX)
       WRITE(*,*) 'RANDOM NUMBERS: LUXURY (MACHINE INDEPENDENT !)'

    CASE(RND_LUX_GAUSS)
       WRITE(*,*) 'RANDOM NUMBERS: LUXURY [GAUSS] (MACHINE INDEPENDENT !)'

    CASE DEFAULT
       WRITE(*,*) '*** ERROR *** UNKNOW RANDOM NUMBER METHOD'
       ! ERROR
       RETURN
    END SELECT

    WRITE(*,*) 'START_SEED(:) : ', START_SEED(:)
    WRITE(*,*) 'NCELL         : ', NCELL
    WRITE(*,*) 'LOOPS         : ', LOOPS

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE RANDOM_READ_NML_CPL
! ----------------------------------------------------------------------

END PROGRAM RANDOM
