! *********************************************************************
! SMIL INTERFACE FOR TEST OF RANDOM NUMBERS
! BASED ON ATTILA
! Version: see "modver" in messy_rndtest.f90
!
! Authors :  Patrick Joeckel, MPICH, DLR, 2003-2017 (RND part of ATTILA)
!            Astrid Kerkweg, Uni Bonn, 2018 (RNDTEST)
!
!
!
! NOTES:

! **********************************************************************
MODULE messy_rndtest_si
! **********************************************************************

  ! ECHAM5/MESSy
  USE messy_main_blather_bi,     ONLY: start_message_bi, end_message_bi &
                                     , error_bi
  USE messy_main_rnd_bi,         ONLY: RND_F90, RND_MTW, RND_LUX      &
                                     , RND_MP_PSJ                     &
                                     , rnd_init_bi, rnd_finish_bi     &
                                     , rnd_number_bi

  ! MESSy
  USE messy_main_constants_mem,  ONLY: dp
  USE messy_main_tools,          ONLY: PTR_2D_ARRAY
  USE messy_rndtest

  IMPLICIT NONE
  PRIVATE
  SAVE

  ! INTRINSIC PROCEDURES
  INTRINSIC :: ASSOCIATED, SIZE, NULL

  ! RANDOM NUMBER HANDLING --------------------------------------------------
  !
  INTEGER, ALLOCATABLE, DIMENSION(:) :: REPR_RNDTEST
  INTEGER, ALLOCATABLE, DIMENSION(:) :: DIMID_NCELL 
  ! NUMBER OF RANDOM NUMBER SERIES
  INTEGER, PARAMETER :: NRPROC     = 4 
  !
  ! RAND_MULT x NCELL RANDOM NUMBERS FOR EACH RANDOM PROCESS
  INTEGER, PARAMETER, DIMENSION(NRPROC) :: RAND_MULT = (/ 1, 1, 1, 3 /)
  !
  TYPE T_RND_NUMBERS
     ! START SEEDS
     INTEGER, DIMENSION(NRPROC) :: START_SEED = (/211169, 7249, 13372, 159941/)
     !
     ! WORKSPACE FOR RANDOM NUMBERS
     TYPE(PTR_2D_ARRAY), DIMENSION(NRPROC) :: HARVEST
     !
     ! PSEUDO RANDOM NUMBER STREAM IDs
     INTEGER, DIMENSION(NRPROC) :: rndid = (/ 0, 0, 0, 0 /)
     ! PSEUDO RANDOM NUMBER GENERATOR
     INTEGER                    :: RND_RNDTEST = RND_MTW  ! default
     ! TEST POINTER
     REAL(DP), DIMENSION(:), POINTER :: test_rnd => NULL()
     ! 
     ! INDEX BOUNDARIES OF NCELL IN NGCELL
     INTEGER, DIMENSION(:,:), POINTER :: IDX => NULL()
     INTEGER                          :: npx
     INTEGER, DIMENSION(:),   POINTER :: V_NCELL => NULL()
     !
  END type T_RND_NUMBERS

  TYPE(T_RND_NUMBERS), ALLOCATABLE, DIMENSION(:) ::  TRND
  !
  INTEGER, PARAMETER, PUBLIC :: I_RANDOM_F90 = 0  ! F90 INTRINSIC
  INTEGER, PARAMETER, PUBLIC :: I_RANDOM_MTW = 1  ! MERSENNE TWISTER
  INTEGER, PARAMETER, PUBLIC :: I_RANDOM_LUX = 2  ! LUXURY
  !
  ! ---------------------------------------------------------------------


  ! CPL-namleist VARIABLES
  LOGICAL :: L_INI_PARALLEL = .TRUE.        ! initialisation in paralell
  INTEGER :: I_RANDOM_METHOD = I_RANDOM_F90 ! WHICH RANDOM NUMBER GENERATOR
  !                                         ! (default: Fortran intrinsic)
  INTEGER :: I_RANDOM_PARALLEL = RND_MP_PSJ ! HOW TO PRALLELIZE
  !                                         ! (default: sync. by jump ahead)
  LOGICAL :: L_RANDOM_TEST = .TRUE.        ! create channel object
  !                                         ! with random numbers for testing
  !

  ! PUBLIC INTERFACE ROUTINES:
  ! ----------------------------------
  PUBLIC :: rndtest_initialize
  PUBLIC :: rndtest_init_memory
  PUBLIC :: rndtest_global_end
  PUBLIC :: rndtest_free_memory

  ! PRIVATE ROUTINES:
  ! -----------------------------------
  !PRIVATE rndtest_read_nml_cpl

CONTAINS

!###########################################################################
! ### PUBLIC SUBROUTINES ###################################################
!###########################################################################
  
!--------------------------------------------------------------------------
  SUBROUTINE rndtest_initialize

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,         ONLY: p_parallel_io, p_io, p_pe, p_bcast
    USE messy_main_transform_bi,   ONLY: get_dc_index
    USE messy_main_channel_bi,     ONLY: n_dom, DC_IX
    USE messy_main_channel_error_bi,   ONLY: channel_halt
    USE messy_main_channel_dimensions, ONLY: new_dimension, DIMID_UNDEF
    USE messy_main_channel_repr,       ONLY: new_representation   &
                                           , set_representation_decomp  &
                                           , IRANK, PIOTYPE_COL, REPR_UNDEF
    ! MESSy
    USE messy_main_tools,          ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'rndtest_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status
    INTEGER                     :: i
    ! PARALLEL DECOMPOSITION
    INTEGER                          :: nseg = 0
    INTEGER, DIMENSION(:,:), POINTER :: start => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()

    CALL start_message_bi(modstr, 'INITIALIZATION', substr)

    ALLOCATE(TRND(n_dom))

    ALLOCATE(NCELL_LOC(n_dom))

    ! INITIALIZE MAIN-CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL rndtest_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi('error in rndtest_read_nml_ctrl', substr)
    END IF

    DO i = 1, n_dom
       IF (p_parallel_io) WRITE (*,*) 'DOMAIN ',i,' got '&
            , NCELL(i),' CELLs', n_dom, SIZE(NCELL)
       CALL p_bcast(NCELL(i), p_io)

       CALL get_dc_index(NCELL(i), TRND(i)%IDX)
       TRND(i)%npx= SIZE(TRND(i)%IDX,1)
       ALLOCATE(TRND(i)%V_NCELL(0:TRND(i)%npx-1))
       TRND(i)%V_NCELL(0:TRND(i)%npx-1) = &
            TRND(i)%IDX(0:TRND(i)%npx-1,2) -TRND(i)% IDX(0:TRND(i)%npx-1,1) + 1
       NCELL_LOC(i) = TRND(i)%V_NCELL(p_pe)
       WRITE (*,*) 'DOMAIN ',i,' PE ', p_pe,' NUMBER OF CELLS: ', NCELL_LOC(i)
    END DO

    ! INITIALIZE CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL rndtest_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi('error in rndtest_read_nml_cpl', substr)
    END IF
    !
    ! CPL NAMELIST
    CALL p_bcast(L_INI_PARALLEL, p_io)
    CALL p_bcast(I_RANDOM_METHOD, p_io)
    CALL p_bcast(I_RANDOM_PARALLEL, p_io)
    CALL p_bcast(L_RANDOM_TEST, p_io)

    ! SET PSEUDO RANDOM NUMBER METHOD (FOR ALL PROCESSES)
    SELECT CASE(I_RANDOM_METHOD)
    CASE(I_RANDOM_F90)
       IF (p_parallel_io) &
            WRITE(*,*) 'RANDOM NUMBERS: F90-INTRINSIC    (MACHINE DEPENDENT!)'
       TRND(:)%RND_RNDTEST = RND_F90
    CASE(I_RANDOM_MTW)
       IF (p_parallel_io) &
            WRITE(*,*) 'RANDOM NUMBERS: MERSENNE TWISTER (MACHINE INDEPENDENT!)'
       TRND(:)%RND_RNDTEST = RND_MTW
    CASE(I_RANDOM_LUX)
       IF (p_parallel_io) &
            WRITE(*,*) 'RANDOM NUMBERS: LUXURY (MACHINE INDEPENDENT!)'
       TRND(:)%RND_RNDTEST = RND_LUX
    CASE DEFAULT
       CALL error_bi('UNKNOW RANDOM NUMBER METHOD', substr)
    END SELECT

    ! NEW DIMENSIONS
    ALLOCATE(DIMID_NCELL(n_dom))
    DIMID_NCELL(:)  = DIMID_UNDEF
    ALLOCATE(REPR_RNDTEST(n_dom))
    REPR_RNDTEST(:) = REPR_UNDEF
    DO i = 1, n_dom
       CALL new_dimension(status, DIMID_NCELL(i), 'NCELL', NCELL(i)    &
#ifdef ICON
       , dom_id=i &
#endif
       )
       CALL channel_halt(substr, status)

       ! NEW REPRESENTATIONS
       CALL new_representation(status, REPR_RNDTEST(i), 'REPR_RNDTEST' &
            , rank = 1, link = 'x---', dctype = DC_IX                  &
            , dimension_ids = (/ DIMID_NCELL(i) /)                     &
            , ldimlen       = (/ NCELL_LOC(i) /)                       &
            , axis = 'N---'                                            &
#ifdef ICON
            , dom_id=i                                                 &
#endif
            )
       CALL channel_halt(substr, status)

       nseg = 1
       ALLOCATE(start(nseg,IRANK))
       ALLOCATE(cnt(nseg,IRANK))
       ALLOCATE(meml(nseg,IRANK))
       ALLOCATE(memu(nseg,IRANK))
       
       start(:,:) = 1
       cnt(:,:) = 1
       meml(:,:) = 1
       memu(:,:) = 1
       
       start(:,1) = TRND(i)%IDX(p_pe,1)
       cnt(:,1)   = NCELL_LOC(i)
       meml(:,1)  = 1
       memu(:,1)  = NCELL_LOC(i)

       CALL set_representation_decomp(status, REPR_RNDTEST(i) &
            , start, cnt, memu, meml, piotype=PIOTYPE_COL)
       CALL channel_halt(substr, status)

       ! ----------------------------------------
       
       DEALLOCATE(start) ; NULLIFY(start)
       DEALLOCATE(cnt)   ; NULLIFY(cnt)
       DEALLOCATE(meml)  ; NULLIFY(meml)
       DEALLOCATE(memu)  ; NULLIFY(memu)

    ! ----------------------------------------
    END DO
    CALL end_message_bi(modstr, 'INITIALIZATION', substr)

  END SUBROUTINE rndtest_initialize
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
  SUBROUTINE rndtest_init_memory

    ! ECHAM5/MESSy
    USE messy_main_channel_error_bi, ONLY: channel_halt 
    ! MESSy
    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         , new_attribute
    USE messy_main_channel_mem,      ONLY: dom_unbound, dom_current

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr ='rndtest_init_memory'
    INTEGER                     :: status
    INTEGER                     :: i
    INTEGER                     :: RND_MP_TEST
    INTEGER                     :: np

    CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)

    IF (dom_current == dom_unbound) THEN
       np = 1
    ELSE
       np = dom_current
    END IF

    ! OPEN NEW OUTPUT CHANNEL
    CALL new_channel(status, modstr)
    CALL channel_halt(substr, status)

    IF (L_RANDOM_TEST) THEN
       CALL new_channel_object(status, modstr, 'test_rnd' &
            , p1=TRND(np)%test_rnd, reprid=REPR_RNDTEST(np))
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'test_rnd' &
            , 'long_name', c='test of random number sequence' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'test_rnd' &
            , 'units', c=' ' )
       CALL channel_halt(substr, status)
    END IF

    CALL end_message_bi(modstr, 'CHANNEL DEFINITION', substr)

    CALL start_message_bi(modstr, 'RANDOM NUMBER INITIALISATION', substr)

    ! ALLOCATE SPACE FOR RANDOM NUMBERS AND
    ! INITIALISE PSEUDO RANDOM NUMBER STREAMS.
    !
    ! NOTE: DO THIS ALWAYS FOR TURB, MC, CONV AND CAT, EVEN IF
    !       LLTBLTURB, LLTMOCA, LLCONV AND/OR LLCAT ARE .FALSE.; 
    !       THIS ALLOWS A MORE SIMPLER CODE FOR INTEGRATION AND
    !       SWITCHING ON/OFF PROCESSES AFTER RESTART WITHOUT
    !       MESSING UP THE DIFFERENT RANDOM NUMBER STREAMS.
    !
    ! ALLOCATE WORKSPACE FOR RANDOM NUMBERS
    DO i=1, NRPROC
       ALLOCATE(TRND(np)%HARVEST(i)%ptr(NCELL_LOC(np),RAND_MULT(i)))
       RND_MP_TEST = I_RANDOM_PARALLEL

       ! use default jump of p_pe*2^256 for RND_RNDTEST = RND_MP_PSJ
       CALL rnd_init_bi(TRND(np)%rndid(i), TRND(np)%RND_RNDTEST&
            , RND_MP_TEST,TRND(np)% START_SEED(i))
    END DO

    CALL end_message_bi(modstr, 'RANDOM NUMBER INITIALISATION', substr)

  END SUBROUTINE rndtest_init_memory
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
  SUBROUTINE rndtest_free_memory

    IMPLICIT NONE

    ! LOCAL
    INTEGER :: i, np

    ! RANDOM NUMBERS AND INDICES
    DO np = 1, SIZE(TRND)
       DO i=1, NRPROC
          IF (ASSOCIATED(TRND(np)%HARVEST(i)%ptr)) &
               DEALLOCATE(TRND(np)%HARVEST(i)%ptr)
          NULLIFY(TRND(np)%HARVEST(i)%ptr)
          CALL rnd_finish_bi(TRND(np)%rndid(i))
       END DO
       IF (ASSOCIATED(TRND(np)%IDX)) THEN
          DEALLOCATE(TRND(np)%IDX) ; NULLIFY(TRND(np)%IDX)
       END IF
       IF (ASSOCIATED(TRND(np)%V_NCELL)) THEN
          DEALLOCATE(TRND(np)%V_NCELL) ; NULLIFY(TRND(np)%V_NCELL)
       END IF
    END DO
    DEALLOCATE(TRND)
    DEALLOCATE(DIMID_NCELL)
    DEALLOCATE(REPR_RNDTEST)
    DEALLOCATE(NCELL_LOC)
    
  END SUBROUTINE rndtest_free_memory
!---------------------------------------------------------------------------

! ===========================================================================
! PRIVATE RNDTEST INTERFACE ROUTINES
! ===========================================================================

! ----------------------------------------------------------------------
  SUBROUTINE rndtest_read_nml_cpl(status, iou)

    ! RNDTEST MODULE ROUTINE 
    !
    ! read namelist for 'coupling' to BASEMODEL
    !
    ! Author: Patrick Joeckel, MPICH, Oct 2003

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'rndtest_read_nml_cpl'

    NAMELIST /CPL/ &
         L_INI_PARALLEL, I_RANDOM_METHOD, I_RANDOM_PARALLEL, L_RANDOM_TEST


    ! LOCAL
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    ! ERROR STATUS
    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)

    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    !
    ! CHECK NAMELIST

    CALL read_nml_close(substr, iou, modstr)


    status = 0  ! no ERROR

  END SUBROUTINE rndtest_read_nml_cpl
! ----------------------------------------------------------------------

!---------------------------------------------------------------------------
  SUBROUTINE rndtest_global_end

    ! MESSy
    USE messy_main_channel_mem,    ONLY: dom_unbound, dom_current

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'rndtest_global_end'
    INTEGER :: i, j, np

    IF (dom_current == dom_unbound) THEN
       np = 1
    ELSE
       np = dom_current
    END IF

    DO i=1, NRPROC
       DO j=1, RAND_MULT(i)
          CALL rnd_number_bi(TRND(np)%rndid(i), TRND(np)%HARVEST(i)%ptr(:,j)&
               , ng=NCELL(np), nl=TRND(np)%V_NCELL(0:TRND(np)%npx-1))
       END DO
       IF (L_RANDOM_TEST .AND. (i == 1)) THEN
          TRND(np)%test_rnd(:) = TRND(np)%HARVEST(i)%ptr(:,1)
       END IF
       !
    END DO

  END SUBROUTINE rndtest_global_end
!---------------------------------------------------------------------------


! ===========================================================================
END MODULE messy_rndtest_si
! ===========================================================================
