program rndj

  ! PROGRAM FOR TESTING THE 'JUMP AHEAD' FACILITY OF THE MERSENNE TWISTER
  ! (see messy_main_rnd_mtw_ja.f90) for efficient parallelisation.
  !
  ! Author: Patrick Joeckel, DLR, Oct 2012
  !
  ! The program runs one 'master' pseudo-random number stream and p_nprocs
  ! additional streams emulating potential parallel substreams.
  ! The sequence of random numbers produced from the master must be
  ! reproduced by the 'chunk-wise' concatenated sequence of the substreams.

  USE messy_main_constants_mem, ONLY: dp
  USE messy_main_rnd,           ONLY: rnd_init, rnd_number &
                                    , rnd_finish, rnd_jump &
                                    , RND_MTW !=> RND_LUX !or !=> RND_F90

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER :: modstr = 'rndj'

  ! CPLT NAMELIST PARAMETERS
  INTEGER :: iseed  = 211169   ! INITIAL SEED
  INTEGER :: NGCELL = 10340    ! NUMBER OF RANDOM NUMBERS
  INTEGER :: p_nprocs = 9      ! NUMBER OF PROCESSORS
  INTEGER :: nt = 5            ! NUMBER OF TIME STEPS
  INTEGER :: nshow = 10        ! RND NUMBERS TO SHOW

  ! INDEX BOUNDARIES NGCELL
  INTEGER, DIMENSION(:,:), POINTER     :: IDX => NULL()
  INTEGER, DIMENSION(:),   ALLOCATABLE :: NCELL

  ! RANDOM NUMBERS TO HARVEST
  ! in one step
  REAL(DP), DIMENSION(:), ALLOCATABLE :: harv
  ! in chunks
  REAL(DP), DIMENSION(:), ALLOCATABLE :: harv_c

  ! RANDOM NUMBER GENERATORS
  INTEGER                            :: id    ! master id
  INTEGER, DIMENSION(:), ALLOCATABLE :: id_c  ! ids for chunks

  ! AUXILIARY
  INTEGER :: status, i, j, it

  ! READ NAMELIST
  CALL read_nml_cpl(Status, 42)
  IF (status /= 0) THEN
     STOP 'ERROR IN NAMELIST'
  END IF
  
  ! MEMORY
  ALLOCATE(NCELL(0:p_nprocs-1))
  ALLOCATE(id_c(0:p_nprocs-1))
  ALLOCATE(harv(NGCELL))
  ALLOCATE(harv_c(NGCELL))

  ! DIVIDE INTO CHUNKS
  CALL get_dc_index(NGCELL, IDX)
  DO i=0, p_nprocs-1
     NCELL(i) = IDX(i,2) - IDX(i,1) + 1
  END DO

  ! CREATE MASTER RANDOM NUMBER STREAM AND GET STATE VECTOR
  CALL rnd_init(status, id, RND_MTW, iseed)
  IF (status /=0) STOP

  ! CREATE RANDOM NUMBER STREAMS FOR CHUNKS
  ! (for first chunk it is in the same state as the master)
  CALL rnd_init(status, id_c(0), RND_MTW, iseed)
  IF (status /=0) STOP
  DO i=1, p_nprocs-1
     CALL rnd_init(status, id_c(i), RND_MTW, iseed)
     IF (status /=0) STOP
     CALL rnd_jump(status, id_c(i), IDX(i-1,2))
     IF (status /=0) STOP
  END DO

  ! HARVEST RANDOM NUMBERS (INITIALISATION)
  CALL rnd_number(id, harv)
  DO i=0, p_nprocs-1
     CALL rnd_number(id_c(i), harv_c(IDX(i,1):IDX(i,2)))
  END DO
  !
  write(*,*) '-------------------------------------------------------------'
  write(*,*) 'RANDOM NUMBERS FOR INITIALISATION (ONE STREAM VS. '&
       ,p_nprocs,'CHUNKS)'
  DO i=1, MIN(NGCELL,nshow)
     write(*,*) i, harv(i), harv_c(i)
  END DO
  WRITE(*,*) '... ... ...'
  DO i= MAX(1,NGCELL-nshow), NGCELL
     write(*,*) i, harv(i), harv_c(i)
  END DO
  write(*,*) '-------------------------------------------------------------'

  DO it=1, nt ! TIME LOOP

     ! FROM NOW ON, JUMP BY NGCELL-NCELL(.) FOR ALL CHUNKS
     DO i=0, p_nprocs-1
        CALL rnd_jump(status, id_c(i), NGCELL-NCELL(i))
        IF (status /=0) STOP
     END DO

     ! HARVEST RANDOM NUMBERS
     CALL rnd_number(id, harv)
     DO i=0, p_nprocs-1
        CALL rnd_number(id_c(i), harv_c(IDX(i,1):IDX(i,2)))
     END DO
     !
     write(*,*) '-------------------------------------------------------------'
     write(*,*) 'RANDOM NUMBERS FOR STEP ',it,' (ONE STREAM VS. '&
       ,p_nprocs,'CHUNKS)'
     DO i=1, MIN(NGCELL,nshow)
        write(*,*) i, harv(i), harv_c(i)
     END DO
     WRITE(*,*) '... ... ...'
     DO i= MAX(1,NGCELL-nshow), NGCELL
        write(*,*) i, harv(i), harv_c(i)
     END DO
     write(*,*) '-------------------------------------------------------------'
     
  END DO

  ! FREE MEMORY
  CALL rnd_finish(id)
  DO i=0, p_nprocs-1
     CALL rnd_finish(id_c(i))
  END DO
  DEALLOCATE(IDX); NULLIFY(IDX)
  !
  DEALLOCATE(NCELL)
  DEALLOCATE(id_c)
  DEALLOCATE(harv)
  DEALLOCATE(harv_c)

CONTAINS

  ! =========================================================================
  SUBROUTINE read_nml_cpl(status, iou)

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    ! PURPOSE:
    !   READ NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    ! NAMELIST CTRL
    NAMELIST /CPL/ iseed &  ! INITIAL SEED
         , NGCELL        &  ! NUMBER OF RANDOM NUMBERS
         , p_nprocs      &  ! NUMBER OF PROCESSORS
         , nt            &  ! NUMBER OF TIME STEPS
         , nshow            ! RND NUMBERS TO SHOW

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER       :: substr = 'read_nml_cpl'
    LOGICAL                           :: lex          ! file exists ?
    INTEGER                           :: fstat        ! file status


    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    WRITE(*,*) 'GLOBAL SETTINGS'
    WRITE(*,*) '---------------'
    WRITE(*,*) 'INITIAL SEED                           : ',iseed
    WRITE(*,*) 'NUMBE OF RANDOM NUMBERS (per time step): ', NGCELL
    WRITE(*,*) 'NUMBER OF PARALLEL PRND STREAMS        : ',p_nprocs
    WRITE(*,*) 'NUMBER OF (TIME) STEPS                 : ',nt
    WRITE(*,*) 'SHOW FIRST/LAST ... OF RANDOM NUMBERS  : ',nshow

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE read_nml_cpl
 ! =========================================================================

 ! =========================================================================
 SUBROUTINE get_dc_index(N, IDX, LDO)

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, PRESENT

    ! I/O
    INTEGER,                                   INTENT(IN) :: N
    INTEGER, DIMENSION(:,:), POINTER                      :: IDX
    LOGICAL, DIMENSION(:),   POINTER, OPTIONAL            :: LDO

    ! LOCAL
    INTEGER :: DN, i

    IF (ASSOCIATED(IDX)) THEN
       DEALLOCATE(IDX)
       NULLIFY(IDX)
    END IF
    ALLOCATE(IDX(0:p_nprocs-1,2))
    ! START INDEX: IDX(p_pe,1)
    ! STOP  INDEX: IDX(p_pe,2)
    IF (PRESENT(LDO)) THEN
       IF (ASSOCIATED(LDO)) THEN
          DEALLOCATE(LDO)
          NULLIFY(LDO)
       END IF
       ALLOCATE(LDO(0:p_nprocs-1))
    END IF

    DN = N/p_nprocs
    DO i=0, p_nprocs-1
       IDX(i,1) = i*DN + 1
       IDX(i,2) = IDX(i,1) + DN - 1
    END DO
    ! CORRECT FOR ODD Ns
    IDX(p_nprocs-1,2) = N

    IF (PRESENT(LDO)) THEN
       DO i=0, p_nprocs-1
          IF (IDX(i,1) > IDX(i,2)) THEN
             LDO(i) = .false.
          ELSE
             LDO(i) = .true.
          END IF
       END DO
    END IF

  END SUBROUTINE get_dc_index
  ! =========================================================================

end program rndj
