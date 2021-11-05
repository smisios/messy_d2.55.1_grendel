!*****************************************************************************

! Author: R. Sander, based on random number generator in
! messy_attila.f90 by P. Joeckel, M. Traub, and A. Pozzer

!*****************************************************************************

MODULE messy_main_rnd

  USE messy_main_constants_mem, ONLY: DP
  USE messy_main_tools,         ONLY: PTR_1D_ARRAY_INT
  USE messy_main_rnd_mtw_ja,    ONLY: nn

  IMPLICIT NONE
  SAVE
  PRIVATE

  PUBLIC  :: rnd_init   ! initialize
  PUBLIC  :: rnd_number ! generate array of random numbers
  PUBLIC  :: rnd_seed   ! handle state vector

  INTERFACE rnd_jump
     MODULE PROCEDURE rnd_jump_n  ! jump by arbitrary n
     MODULE PROCEDURE rnd_jump_2  ! jump by n * 2^p
  END INTERFACE
  PUBLIC  :: rnd_jump   ! jump ahead
  PUBLIC  :: rnd_finish ! deallocate

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'rnd'
  CHARACTER(len=*), PARAMETER, PUBLIC :: modver = '1.1'
  LOGICAL, PUBLIC                     :: L_VERBOSE = .TRUE.

  INTEGER, PARAMETER, PUBLIC :: &
    RND_F90       = 1, & ! Fortran90 intrinisc      (architecture   dependent)
    RND_MTW       = 2, & ! Mersenne Twister MT19937 (architecture independent)
    RND_LUX       = 3, & ! Luxury                   (architecture independent)
    RND_F90_GAUSS = 4, & ! Gauss normal distribution, based on F90
    RND_MTW_GAUSS = 5, & ! Gauss normal distribution, based on MTW
    RND_LUX_GAUSS = 6    ! Gauss normal distribution, based on LUX

  INTEGER, PARAMETER, PUBLIC :: RND_MAX_METHOD = 6
  CHARACTER(LEN=12), DIMENSION(RND_MAX_METHOD), PUBLIC :: RND_REPR_NAME = &
       (/ 'REPR_RND_F90', 'REPR_RND_MTW', 'REPR_RND_LUX', &
          'REPR_RND_F90', 'REPR_RND_MTW', 'REPR_RND_LUX'  /)
  CHARACTER(LEN=11), DIMENSION(RND_MAX_METHOD), PUBLIC :: RND_DIM_NAME = &
       (/ 'DIM_RND_F90', 'DIM_RND_MTW', 'DIM_RND_LUX', &
          'DIM_RND_F90', 'DIM_RND_MTW', 'DIM_RND_LUX'  /)

  ! The following variables are for different random number series,
  ! each with their own id:
  INTEGER, PARAMETER, PUBLIC :: ID_MAX  = 20 ! max number of series
  INTEGER, PUBLIC :: nstate(ID_MAX)     =  0 ! dimension of state vector
  INTEGER, PUBLIC :: rnd_method(ID_MAX)
  TYPE(PTR_1D_ARRAY_INT), DIMENSION(ID_MAX), PUBLIC :: state

  ! SAVE POLYNOM COEFFICIENTS UP TO 2^N_MAX_EXP ... TO AVOID RECALCULATION
  INTEGER, PARAMETER                  :: N_MAX_EXP = 20
  LOGICAL                             :: lpoly(N_MAX_EXP) = .FALSE.
  INTEGER, DIMENSION(0:nn, N_MAX_EXP) :: pcoeff ! coefficients

CONTAINS

  !***************************************************************************

  SUBROUTINE rnd_init(status, id, method, pseed)

    USE messy_main_rnd_mtw, ONLY: init_genrand, mt, mti, wi
    USE messy_main_rnd_lux, ONLY: rluxgo, isdext, rluxut

    IMPLICIT NONE

    INTEGER,           INTENT(OUT) :: status
    INTEGER,           INTENT(OUT) :: id
    INTEGER,           INTENT(IN)  :: method
    INTEGER, OPTIONAL, INTENT(IN)  :: pseed
    REAL :: randomseed
    INTEGER :: seed


    IF (L_VERBOSE) THEN
      WRITE(*,*)
      WRITE(*,*) 'Initializing random number generator:'
    ENDIF
    status = 0 ! ok

    ! define seed:
    IF (PRESENT(pseed)) THEN
      seed = pseed
      IF (L_VERBOSE) WRITE(*,*) 'seed   = ', seed
    ELSE
      CALL RANDOM_SEED
      CALL RANDOM_NUMBER(randomseed)
      seed = INT(HUGE(1)*randomseed)
      IF (L_VERBOSE) WRITE(*,*) 'seed   = ', seed, ' (computer-generated)'
    ENDIF

    ! assign an id to this series of random numbers:
    DO id=1,ID_MAX
      IF (nstate(id)==0) EXIT ! exit do loop
    ENDDO
    IF (L_VERBOSE) WRITE(*,*) 'id     = ', id
    IF (id>ID_MAX) THEN
      WRITE(*,*) 'rnd_init: ERROR: Too many series. Increase ID_MAX!'
      status = 1
      RETURN
    ENDIF

    ! allocate memory for the state array, and set the seed:
    rnd_method(id) = method
    IF (L_VERBOSE) WRITE(*,'(A)', ADVANCE='NO') ' method = '
    SELECT CASE(rnd_method(id))
    CASE(RND_F90,RND_F90_GAUSS)
      CALL RANDOM_SEED(put=SPREAD(seed,1,4))
      IF (L_VERBOSE) WRITE(*,'(A)', ADVANCE='NO') 'Fortran90 intrinisc'
      CALL RANDOM_SEED(size = nstate(id))
      ALLOCATE(state(id)%ptr(nstate(id)))
      state(id)%ptr(:) = seed
    CASE(RND_MTW,RND_MTW_GAUSS)
      IF (L_VERBOSE) WRITE(*,'(A)', ADVANCE='NO') 'Mersenne Twister MT19937'
      nstate(id) = SIZE(mt)+1
      ALLOCATE(state(id)%ptr(nstate(id)))
      CALL init_genrand(INT(seed,wi))
      mti = SIZE(mt)+1
      state(id)%ptr(1:624) = INT(mt(1:624))
      state(id)%ptr(625)   = INT(mti)
    CASE(RND_LUX,RND_LUX_GAUSS)
      IF (L_VERBOSE) WRITE(*,'(A)', ADVANCE='NO') 'Luxury'
      nstate(id) = SIZE(isdext)
      ALLOCATE(state(id)%ptr(nstate(id)))
      CALL rluxgo(4,seed,0,0)
      CALL rluxut
      state(id)%ptr(:) = INT(isdext(:))
    CASE default
      WRITE(*,*) 'rnd_init: ERROR: unknown method; id = ', id, &
           '; method = ', rnd_method(id)
      status = 2
    END SELECT

    IF (L_VERBOSE) THEN
      IF ((rnd_method(id)==RND_F90_GAUSS) .OR. &
        (rnd_method(id)==RND_MTW_GAUSS) .OR. &
        (rnd_method(id)==RND_LUX_GAUSS)) THEN
        WRITE(*,*) ' (Gauss normal distribution)'
      ELSE
        WRITE(*,*) ''
      ENDIF

      WRITE(*,*) 'nstate = ', nstate(id)
    ENDIF

  END SUBROUTINE rnd_init

  !***************************************************************************

  SUBROUTINE rnd_seed(id, size, put, get)

    IMPLICIT NONE

    INTEGER,           INTENT(IN)                :: id
    INTEGER, OPTIONAL, INTENT(OUT)               :: size
    INTEGER, OPTIONAL, INTENT(IN),  DIMENSION(:) :: put
    INTEGER, OPTIONAL, INTENT(OUT), DIMENSION(:) :: get

    IF (PRESENT(size)) THEN
      size = nstate(id)
    ENDIF

    IF (PRESENT(put)) THEN
      state(id)%ptr(:) = put(:)
    ENDIF

    IF (PRESENT(get)) THEN
      get(:) = state(id)%ptr(:)
    ENDIF

  END SUBROUTINE rnd_seed

  !***************************************************************************

  SUBROUTINE rnd_jump_n(status, id, n, get)

    USE messy_main_rnd_mtw_ja,    ONLY: mt_jump_ahead, nn

    IMPLICIT NONE
    INTRINSIC :: BIT_SIZE, BTEST

    ! I/O
    INTEGER, INTENT(OUT) :: status
    INTEGER, INTENT(IN)  :: id
    INTEGER, INTENT(IN)  :: n
    INTEGER, DIMENSION(:), INTENT(OUT), OPTIONAL :: get

    ! LOCAL
    REAL(DP), DIMENSION(:), ALLOCATABLE :: harvest
    REAL(DP), DIMENSION(1)              :: h1
    !
    INTEGER :: n_check
    INTEGER :: ja  ! JUMP AHEAD 2**ja RANDOM NUMBERS ...
    INTEGER,  DIMENSION(0:nn-1) :: v, w  ! state vectors

    ! NOTE: FOR THE FORTRAN INTRINSIC AND LUXURY RANDOM NUMBER GENERATORS
    !       NO FANCY 'JUMP AHEAD' FACILITY IS AVAILABLE AT THE TIME BEING.
    !       THIS ONLY THE 'POOR MAN'S' JUMP AHEAD (I.E., WASTING N RANDOM
    !       NUMBERS) IS IMPLEMENTED. THIS DOES NOT SCALE, OF COURSE, IF
    !       USED FOR PARALLELISATION!
    !       DEPENDING ON THE REQUESTED JUMP N, THE ARRAY SIZE COULD 
    !       BECOME TOO BIG. THE SINGE HARVEST COULD THEN POTENTIALLY
    !       REPLACED BY A LOOP.
    ! NOTE: FOR THE GAUSSIAN DISTRIBUTIONS JUMP AHEAD IS NOT POSSIBLE,
    !       BECAUSE A PRIORI IT IS UNKNOWN HOW MANY RNDOM NUMBERS ARE
    !       REQUIRED FOR THE Marsaglia polar method.

    status = 0

    IF (n == 0) THEN
       IF (PRESENT(get)) THEN
          get(:) = state(id)%ptr(:)
       ENDIF
       RETURN
    ENDIF

    SELECT CASE(rnd_method(id))

    CASE(RND_F90, RND_LUX)
       ALLOCATE(harvest(n))
       CALL rnd_number(id, harvest, get)
       DEALLOCATE(harvest)

    CASE(RND_F90_GAUSS, RND_LUX_GAUSS, RND_MTW_GAUSS)
       WRITE(*,*) 'rnd_jump: ERROR: jump ahead not possible for '//&
            &'GAUSSIAN DISTRIBUTION BASED ON Marsaglia polar method;'//&
            &' id = ', id, '; method = ', rnd_method(id)
       status = 1

    CASE(RND_MTW)

       ! COPY STATE VECTOR
       v(0:nn-1) = state(id)%ptr(1:nn)
       n_check = 0

       ! THE JUMP AHEAD CAN ONLY BE IN POWERS OF 2
       DO ja = BIT_SIZE(n)-1, 1, -1

          bit: IF (BTEST(n,ja)) THEN

             n_check = n_check + 2**ja
             !write(*,*) 'rnd_jump: ', 'id = ',id, '; n = ', n, '; ja = ', ja &
             !     , '; 2**ja = ', 2**ja, '; sum = ',n_check

             ! NOTE: The factor of 2 in the call of mt_jump_ahead is needed, 
             !       because generation of DP reals proceeds the state vector
             !       by 64-bit.

             ! ADVANCE STATE VECTOR, BUT ...
             ! ... RECALCULATE POLYNOMIAL COEFFICIENTS, ONLY IF NOT YET DONE ...
             IF (ja > N_MAX_EXP) THEN
                ! always, if out of range
                CALL mt_jump_ahead(v(0:nn-1), w(0:nn-1), ja, 2)
             ELSE
                IF (lpoly(ja)) THEN
                   ! already available
                   CALL mt_jump_ahead(v(0:nn-1), w(0:nn-1), ja, 2 &
                        , putp=pcoeff(0:nn,ja)) 
                ELSE
                   ! not yet available
                   CALL mt_jump_ahead(v(0:nn-1), w(0:nn-1), ja, 2 &
                        , getp=pcoeff(0:nn,ja))
                   lpoly(ja) = .TRUE. ! from now on available
                END IF
             END IF

             ! COPY RESULT FOR NEXT ITERATION
             v(0:nn-1) = w(0:nn-1)

          END IF bit

       END DO

       ! WRITE BACK STATE VECTOR
       state(id)%ptr(1:nn) = v(0:nn-1)

       ! IF REST IS ONE, ADVANCE STATE VECTOR BY HARVESTING ONE PSEUDO
       ! RANDOM NUMBER ...
       IF (BTEST(n, 0)) THEN
          n_check = n_check + 1
          !write(*,*) 'rnd_jump: ', 'id = ',id, '; n = ', n, '; ja = ', 0 &
          !     , '; 2**ja = ', 2**0, '; sum = ',n_check
          CALL rnd_number(id, h1, get)
       END IF

       IF (n_check /= n) THEN
          WRITE(*,*) 'rnd_jump: ERROR: decomposition of vector length failed', &
               ' ( n = ', n, '; n_check = ',n_check,' )'
          status = 2
          RETURN
       END IF

    CASE default
      WRITE(*,*) 'rnd_jump: ERROR: unknown method; id = ', id, &
           '; method = ', rnd_method(id)
      status = 3
    END SELECT

    ! RETURN STATE VECTOR
    IF (PRESENT(get)) THEN
       get(:) = state(id)%ptr(:)
    ENDIF

  END SUBROUTINE rnd_jump_n

  !***************************************************************************

  SUBROUTINE rnd_jump_2(status, id, n, p, get)

    USE messy_main_rnd_mtw_ja,    ONLY: mt_jump_ahead, nn

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status
    INTEGER, INTENT(IN)  :: id
    INTEGER, INTENT(IN)  :: n
    INTEGER, INTENT(IN)  :: p
    INTEGER, DIMENSION(:), INTENT(OUT), OPTIONAL :: get

    ! LOCAL
    INTEGER,  DIMENSION(0:nn-1) :: v, w  ! state vectors

    ! NOTE: FOR THE FORTRAN INTRINSIC AND LUXURY RANDOM NUMBER GENERATORS
    !       NO FANCY 'JUMP AHEAD' FACILITY IS AVAILABLE AT THE TIME BEING.
    !       SINCE THIS ROUTINE IS USED FOR LARGE JUMPS (TO ENABLE PARALLEL,
    !       NON-OVERLAPPING SUBSTREAMS), THE 'POOR MAN'S' JUMP AHEAD 
    !       (I.E., WASTING N RANDOM NUMBERS) IS NOT FEASIBLE HERE!

    status = 0

    SELECT CASE(rnd_method(id))

    CASE(RND_F90, RND_F90_GAUSS, RND_LUX, RND_LUX_GAUSS)
       WRITE(*,*) 'rnd_jump: ERROR: large jump ahead not feasible for '//&
            &'Fortran intrinsic or Luxury generators;'//&
            &' id = ', id, '; method = ', rnd_method(id)
       status = 1

    CASE(RND_MTW, RND_MTW_GAUSS)

       ! COPY STATE VECTOR
       v(0:nn-1) = state(id)%ptr(1:nn)

       ! NOTE: The factor of 2 in the call of mt_jump_ahead is needed, 
       !       because generation of DP reals proceeds the state vector
       !       by 64-bit.
       ! -> jump ahead 2*n * 2^p
       CALL mt_jump_ahead(v(0:nn-1), w(0:nn-1), p, 2*n)

       ! WRITE BACK STATE VECTOR
       state(id)%ptr(1:nn) = w(0:nn-1)

    CASE default
      WRITE(*,*) 'rnd_jump: ERROR: unknown method; id = ', id, &
           '; method = ', rnd_method(id)
      status = 3
    END SELECT

    ! RETURN STATE VECTOR
    IF (PRESENT(get)) THEN
       get(:) = state(id)%ptr(:)
    ENDIF

  END SUBROUTINE rnd_jump_2

  !***************************************************************************

  SUBROUTINE rnd_number(id, harvest, get)

    USE messy_main_rnd_mtw, ONLY: genrand_res53, mt, mti, WI
    USE messy_main_rnd_lux, ONLY: isdext, ranlux &
                                , rluxut, rluxin
    IMPLICIT NONE
    INTRINSIC :: random_seed, random_number

    INTEGER, INTENT(IN) :: id
    REAL(DP), DIMENSION(:), INTENT(OUT)           :: harvest
    INTEGER,  DIMENSION(:), INTENT(OUT), OPTIONAL :: get
    INTEGER :: j, n
    REAL, DIMENSION(:), ALLOCATABLE :: temp_array
    REAL(DP) :: p, q, y1, y2 ! for Marsaglia polar method

    n = SIZE(harvest)
    ! ------------------------------------------------------------------------
    SELECT CASE(rnd_method(id))
    ! ------------------------------------------------------------------------
    CASE(RND_F90)
      CALL RANDOM_SEED(put = state(id)%ptr(:))
      CALL RANDOM_NUMBER(harvest)
      CALL RANDOM_SEED(get = state(id)%ptr(:))
    ! ------------------------------------------------------------------------
    CASE(RND_MTW)
      mt(1:624) = INT(state(id)%ptr(1:624), WI)
      mti       = INT(state(id)%ptr(625),   WI)
      DO j = 1, n
        harvest(j) = REAL(genrand_res53(), DP)
      ENDDO
      state(id)%ptr(1:624) = INT(mt(1:624))
      state(id)%ptr(625)   = INT(mti)
    ! ------------------------------------------------------------------------
    CASE(RND_LUX)
      isdext(:) = INT(state(id)%ptr(:))
      CALL rluxin
      ALLOCATE(temp_array(n))
      CALL ranlux(temp_array, n)
      harvest(:) = REAL(temp_array(:), DP)
      DEALLOCATE(temp_array)
      CALL rluxut
      state(id)%ptr(:) = INT(isdext(:))
    ! ------------------------------------------------------------------------
    CASE(RND_F90_GAUSS)
      ! create normal distribution with Marsaglia polar method, see:
      ! http://de.wikipedia.org/wiki/Polar-Methode
      DO j=1,n
        DO
          CALL RANDOM_NUMBER (y1) ! produce a number between 0 and 1
          CALL RANDOM_NUMBER (y2) ! produce a number between 0 and 1
          q = (2.*y1-1)**2 + (2.*y2-1)**2
          IF (q<=1) EXIT
          ! if q>1, repeat calculation with new y1, y2
        ENDDO
        p = SQRT((-2.*LOG(q)/q))
        harvest(j)=(2.*y1-1.) * p
        !z2=(2.*y2-1.) * p ! z2 is not used here
      ENDDO
    ! ------------------------------------------------------------------------
   CASE(RND_LUX_GAUSS)
      ! create normal distribution with Marsaglia polar method, see:
      ! http://de.wikipedia.org/wiki/Polar-Methode
      isdext(:) = INT(state(id)%ptr(:))
      ALLOCATE(temp_array(2))
      DO j=1,n
        DO
          CALL ranlux(temp_array, 2)
          y1 = REAL(temp_array(1), DP)
          y2 = REAL(temp_array(2), DP)
          q = (2.*y1-1)**2 + (2.*y2-1)**2
          IF (q<=1) EXIT
          ! if q>1, repeat calculation with new y1, y2
        ENDDO
        p = SQRT((-2.*LOG(q)/q))
        harvest(j)=(2.*y1-1.) * p
        !z2=(2.*y2-1.) * p ! z2 is not used here
      ENDDO
      DEALLOCATE(temp_array)
      state(id)%ptr(:) = INT(isdext(:))
    ! ------------------------------------------------------------------------
    CASE(RND_MTW_GAUSS)
      ! create normal distribution with Marsaglia polar method, see:
      ! http://de.wikipedia.org/wiki/Polar-Methode
      mt(1:624) = INT(state(id)%ptr(1:624), WI)
      mti       = INT(state(id)%ptr(625),   WI)
      DO j = 1, n
        DO
          y1 = REAL(genrand_res53(), DP) ! produce a number between 0 and 1
          y2 = REAL(genrand_res53(), DP) ! produce a number between 0 and 1
          q = (2.*y1-1)**2 + (2.*y2-1)**2
          IF (q<=1) EXIT
          ! if q>1, repeat calculation with new y1, y2
        ENDDO
        p = SQRT((-2.*LOG(q)/q))
        harvest(j)=(2.*y1-1.) * p
        !z2=(2.*y2-1.) * p ! z2 is not used here
      ENDDO
      state(id)%ptr(1:624) = INT(mt(1:624))
      state(id)%ptr(625)   = INT(mti)
    CASE default
      WRITE(*,*) 'rnd_number: ERROR: unknown method; id = ', id, &
           '; method = ', rnd_method(id)
    ! ------------------------------------------------------------------------
    END SELECT
    ! ------------------------------------------------------------------------

    IF (PRESENT(get)) THEN
      get(:) = state(id)%ptr(:)
    ENDIF

  END SUBROUTINE rnd_number

  !***************************************************************************

  SUBROUTINE rnd_finish(id)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: id

    DEALLOCATE(state(id)%ptr)
    nstate(id) = 0 ! reset to make id available again

  END SUBROUTINE rnd_finish

  !***************************************************************************

END MODULE messy_main_rnd

!*****************************************************************************
