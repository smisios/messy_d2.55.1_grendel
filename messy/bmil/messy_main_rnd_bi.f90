!*****************************************************************************

! Author: P. Joeckel, DLR, May 2012
! BMIL file (for restart handling of state vectors)

!*****************************************************************************

MODULE messy_main_rnd_bi

  USE messy_main_constants_mem, ONLY: DP
  USE messy_main_blather_bi,    ONLY: error_bi, warning_bi, info_bi
  USE messy_main_tools,         ONLY: PTR_2D_ARRAY
  USE messy_main_rnd

  IMPLICIT NONE
  PRIVATE
  SAVE

  TYPE(PTR_2D_ARRAY), DIMENSION(ID_MAX) :: rstate
  LOGICAL                               :: lfirst_jump(ID_MAX) = .FALSE.

  ! METHODS FOR PARALLELISATION (see below)
  INTEGER, PARAMETER, PUBLIC :: RND_MP_UND = -1 ! undefined
  INTEGER, PARAMETER, PUBLIC :: RND_MP_SEQ = 1  ! sequential / scatter
  INTEGER, PARAMETER, PUBLIC :: RND_MP_PIN = 2  ! parallel individual
  INTEGER, PARAMETER, PUBLIC :: RND_MP_PSJ = 3  ! parallel sync. by jump ahead
  INTEGER, PARAMETER, PUBLIC :: RND_MP_PIJ = 4  ! paralle indep. by jump ahead
  INTEGER, DIMENSION(ID_MAX) :: MP = RND_MP_UND

  ! USE throught SMCL PARAMETERS
  PUBLIC :: RND_F90, RND_MTW, RND_LUX
  PUBLIC :: RND_F90_GAUSS, RND_MTW_GAUSS, RND_LUX_GAUSS

  ! CALLED BY SMIL
  PUBLIC :: rnd_init_bi     ! initialize pseudo random number stream 
  PUBLIC :: rnd_number_bi   ! harvest random numbers
  PUBLIC :: rnd_finish_bi   ! terminate stream

  ! CALLED BY BMIL
  PUBLIC :: main_rnd_init_coupling   ! create dimensions, representations 
  !                                  ! and channel objects
                                     ! for state vectors (save for restart)
  PUBLIC :: main_rnd_init_tracer     ! reset state vectors after restart
  PUBLIC :: main_rnd_write_output    ! save state vectors for restart

  ! PRIVATE SUBROUTINES
  !PRIVATE :: rnd_number_sja ! parallel harvest, synchronised by jump-ahead

  ! RECIPE:
  ! 1) CALL rnd_init_bi on all PEs to initialise state vector
  !    in initialisation phase of model
  ! 2) CALL rnd_number_bi on all PEs within time loop of model and
  !    harvest 'n' random numbers (see below)
  ! 3) CALL rnd_finish_bi on all PEs to finalise the random
  !    number generator in finalising phase of model
  !
  ! NOTE: Restart invariance is achieved, since the state vectors of
  !       all PEs are saved in the restart files.
  !
  ! CORRECT METHODS (depending on the application):
  ! 
  ! ==========================================================================
  ! METHOD DECOMPOSITION  RESTART     MEMORY          SCAL-      GENERATORS
  !        INDEPENDENT    INVARIANT   REQUIREMENTS    ABILITY
  ! ==========================================================================
  ! SEQ    YES            YES         NG+NL on PE=0   NO         F90(_GAUSS)
  !                                   NL    on PE>0              LUX(_GAUSS)
  !                                                              MTW(_GAUSS)
  ! --------------------------------------------------------------------------
  ! PIN    YES            YES         NG on all PEs   NO         F90(_GAUSS)
  !                                                              LUX(_GAUSS)
  !                                                              MTW(_GAUSS)
  ! --------------------------------------------------------------------------
  ! PSJ    YES            YES         NL on all PEs   ?          MTW
  !                                                              (F90,LUX
  !                                                               not useful)
  ! --------------------------------------------------------------------------
  ! PIJ    NO             YES         NL on all PEs   YES        MTW(_GAUSS)
  ! ==========================================================================
  !  NL = number of harvested random numbers per call on each PE
  !       (NL can be different on different PEs, but SUM_pe(NL) = NG
  !  NG = total number of harvested random numbers per call on all PEs
  !  np = number of PEs
  ! --------------------------------------------------------------------------

  ! SEQ) SEQUENTIAL (HARVEST ON ONE PE AND SCATTER): 
  !      Harvest NG random numbers on one dedicated PE (e.g., p_io) and
  !      scatter them (in chunks of NL) to individual PEs. 
  !      This delivers sequences, which are independent on the parallel 
  !      decomposition, however, the method does not scale with the number of 
  !      PEs.
  !
  !         CALL rnd_inint_bi(is, ... , RND_MP_SEQ, SEED)
  !         ...
  !         IF (p_pe == p_io) THEN
  !            ALLOCATE(g_harvest(NG))
  !         ENDIF
  !         CALL rnd_number_bi(id, harvest, p_io)
  !         ALLOCATE(l_harvest(NL))
  !         CALL scatter_glix(g_harvest, l_harvest, p_io)
  !         ...
  !         DEALLOCATE(...)
  !         CALL rnd_finish_bi(id)
  !
  ! PIN) PARALLEL INDIVIDUAL (HARVEST (THE SAME) SEQUENCE INDIVIDUALLY
  !      ON ALL PEs IN PARALLEL):
  !      Harvest NG random random numbers on each PE. All parallel streams 
  !      produce their own sequence. If the SEED is the same for all PEs,
  !      all sequences are equal. Note that even with different
  !      SEEDs an overlap between the sequences cannot be ruled out.
  !      This potentially causes undesired correlations!
  !
  !         CALL rnd_inint_bi(is, ... , RND_MP_PIN, SEED)
  !         ...
  !         ALLOCATE(harvest(NG))              ! NG equal on all PEs
  !         CALL rnd_number_bi(id, harvest)
  !         ...
  !         DEALLOCATE(...)
  !         CALL rnd_finish_bi(id)
  !
  ! PSJ) PARALLEL HARVEST SYNCHRONISED WITH JUMP AHEAD:
  !      Harvest NL random numbers on each PE with rnd_number_p_bi.
  !      Parallel pseudo random number streams are operated on all PEs,
  !      the jump ahead mechanism is used to guarantee the same sequence as
  !      with a sequential execution (SEQ).
  !      Further information on ng (=n_size) and nl (=n_chunk), 
  !      see subroutine rnd_number_sja below.
  !
  !         CALL rnd_inint_bi(is, ... , RND_MP_PSJ, SEED)
  !         ...
  !         ALLOCATE(harvest(NL))
  !         CALL rnd_number_bi(id, harvest, ng, nl(0:))
  !         ...
  !         DEALLOCATE(...)
  !         CALL rnd_finish_bi(id)
  !
  ! PIJ) PARALLEL INDEPENDENT STREAMS, NON OVERLAPPING DUE TO JUMP AHEAD:
  !      The sequences (of methot PIN) on different PEs can be made 
  !      non-overlapping by a sufficiently large, PE depending 
  !      jump-ahead right after the initialisation of the stream. 
  !      In this case, the sequences depend on the parallel decomposition.
  !
  !         ! jump ahead to sub-stream starting in p*2^256 steps
  !         CALL rnd_inint_bi(is, ... , RND_MP_PIJ, SEED, p_pe, n=256)
  !         ...
  !         ALLOCATE(harvest(NL))
  !         CALL rnd_number_bi(id, harvest)
  !         ...
  !         DEALLOCATE(...)
  !         CALL rnd_finish_bi(id)

CONTAINS

! =======================================================================
! ### SMIL ENTRY POINTS
! =======================================================================

! -----------------------------------------------------------------------
  SUBROUTINE rnd_init_bi(id, method, pmethod, pseed, n, p, get)

    USE messy_main_mpi_bi, ONLY: p_pe

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT)                         :: id
    INTEGER, INTENT(IN)                          :: method
    INTEGER, INTENT(IN)                          :: pmethod
    INTEGER, INTENT(IN),                OPTIONAL :: pseed
    INTEGER, INTENT(IN),                OPTIONAL :: n      ! jump ahead n*2^p
    INTEGER, INTENT(IN),                OPTIONAL :: p      ! ...
    INTEGER, DIMENSION(:), INTENT(OUT), OPTIONAL :: get

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'rnd_init_bi'
    INTEGER :: status
    INTEGER :: zn, zp
    CHARACTER(LEN=4) :: znstr, zpstr

    CALL rnd_init(status, id, method, pseed)
    if (status /= 0) &
         CALL error_bi('rnd_init reported an error',substr)
    
    MP(id) = pmethod
    SELECT CASE(MP(id))

       CASE(RND_MP_SEQ)
          CALL info_bi('sequential harvest and scatter' &
               , substr, .TRUE.)

       CASE(RND_MP_PIN)
          CALL info_bi( &
               'individual harvest of the same sequence in parallel' &
               , substr, .TRUE.)

       CASE(RND_MP_PSJ)
          CALL info_bi( &
               'shared harvest of one sequence in parallel; '//&
               &'synchronized by jump-ahead', substr, .TRUE.)

       CASE(RND_MP_PIJ)

          IF (PRESENT(n)) THEN
             zn = n
          ELSE
             zn = p_pe
          ENDIF
          IF (PRESENT(p)) THEN
             zp = p
          ELSE
             zp = 256
          END IF

          WRITE(znstr,*) zn
          WRITE(zpstr,*) zp
          CALL info_bi( &
               'individual harvest of subsequences in parallel; '//&
               &'non-overlapping by jump-ahead of'//&
               &znstr//' * 2^'//zpstr//' steps', substr, .TRUE.)

          CALL rnd_jump(status, id, zn, zp, get)
          IF (status /= 0) &
               CALL error_bi('rnd_jump reported an error',substr)

       CASE DEFAULT
          call error_bi('unknown parallelisaton method for rnd stream',substr)
    END SELECT

    IF (PRESENT(get)) THEN
       get(:) = state(id)%ptr(:)
    ENDIF

  END SUBROUTINE rnd_init_bi
! -----------------------------------------------------------------------

! -----------------------------------------------------------------------
  SUBROUTINE rnd_number_bi(id, harvest, pe, ng, nl, get)

    USE messy_main_mpi_bi,       ONLY: p_io, p_pe, p_nprocs

    IMPLICIT NONE
    
    ! I/O
    INTEGER,                           INTENT(IN)            :: id
    REAL(DP), DIMENSION(:), TARGET,    INTENT(OUT)           :: harvest
    INTEGER,                           INTENT(IN),  OPTIONAL :: pe
    INTEGER,                           INTENT(IN),  OPTIONAL :: ng
    INTEGER,  DIMENSION(0:p_nprocs-1), INTENT(IN),  OPTIONAL :: nl
    INTEGER,  DIMENSION(:),            INTENT(OUT), OPTIONAL :: get

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'rnd_number_bi'
    INTEGER :: zpe

    IF (id > ID_MAX) &
         CALL error_bi('pseudo random number stream id out of range',substr)

    IF (nstate(id) == 0) &
         CALL error_bi('pseudo random number stream id not initialised',substr)

    SELECT CASE(MP(id))

       CASE(RND_MP_SEQ) ! -----------------------------------------------------

          IF (.NOT. PRESENT(ng)) &
               CALL error_bi('ng must be present for RND_MP_SEQ',substr)

          IF (PRESENT(pe)) THEN
             zpe = pe
          ELSE
             zpe = p_io
          END IF

          IF (p_pe == zpe) THEN
             CALL rnd_number(id, harvest, get)
          END IF

       CASE(RND_MP_PIN) ! -----------------------------------------------------

          CALL rnd_number(id, harvest, get)

       CASE(RND_MP_PSJ) ! -----------------------------------------------------

          CALL rnd_number_sja(id, ng, nl, harvest, get)

       CASE(RND_MP_PIJ) ! -----------------------------------------------------

          CALL rnd_number(id, harvest, get)

       CASE DEFAULT ! ---------------------------------------------------------

          call error_bi('unknown parallelisaton method for rnd stream',substr)

    END SELECT

  END SUBROUTINE rnd_number_bi
! -----------------------------------------------------------------------

! -----------------------------------------------------------------------
  SUBROUTINE rnd_finish_bi(id)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: id

    MP(id) = RND_MP_UND
    CALL rnd_finish(id)

  END SUBROUTINE rnd_finish_bi
! -----------------------------------------------------------------------

! -----------------------------------------------------------------------
! --- PRIVATE SUBROUTINES
! -----------------------------------------------------------------------

! -----------------------------------------------------------------------
  SUBROUTINE rnd_number_sja(id, n_size, n_chunk, harvest, get)

    USE messy_main_mpi_bi, ONLY: p_pe, p_nprocs !!$, p_bcast

    ! Author: Patrick Joeckel, DLR, Oct 2012
    !
    ! THIS SUBROUTINE USES THE JUMP AHEAD FACILITY (efficiently ONLY
    ! implemented for the Mersenne Twister and NOT usable for GAUSSian
    ! distributions created with the Marsaglia polar method!)
    ! to harvest pseudo rundom number sequences, which are independent
    ! on the degree of parallelisation.
    !
    ! NOTES: - n_size is the total number (summed over all PEs) of random
    !          numbers to harvest.
    !        - n_cunk(:) contains the number of random numbers to harvest
    !          on each PE.
    !          This information is NOT calulated here, because this would
    !          require a lot of expensive calls to 'MPI_BCAST'. Usually this
    !          information (the decomposition) is already available from 
    !          the application code. 
    !       EXAMPLES:
    !       - ATTILA: n_size = NGCELL
    !                 n_chunk = (/NCELL(PE=0), NCELL(PE=1) ... /)
    !                 ( SUM(NCELL(0:p_nprocs-1)) = NGCELL )
    !       - ...

    IMPLICIT NONE
    INTRINSIC :: SUM

    ! I/O
    INTEGER,                           INTENT(IN)            :: id
    INTEGER,                           INTENT(IN)            :: n_size
    INTEGER,  DIMENSION(0:p_nprocs-1), INTENT(IN)            :: n_chunk
    REAL(DP), DIMENSION(:),            INTENT(OUT)           :: harvest
    INTEGER,  DIMENSION(:),            INTENT(OUT), OPTIONAL :: get

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'rnd_number_sja'
    INTEGER,  DIMENSION(0:p_nprocs-1) :: njump
    INTEGER :: np
    INTEGER :: status

    IF (SIZE(n_chunk) /= p_nprocs) &
         CALL error_bi('NUMBER OF CHUNKS DOES NOT MATCH NUMBER OF PEs',substr)

    IF (SIZE(harvest) /= n_chunk(p_pe)) &
         CALL error_bi('SIZE OF HARVEST DOES NOT MATCH CHUNK-SIZE',substr)

    IF (SUM(n_chunk(0:p_nprocs-1)) /= n_size) &
         CALL error_bi('SUM OF CHUNK SIZES DOES NOT MATCH TOTAL SIZE', substr)

    IF (SIZE(harvest) == 0) THEN
       IF (PRESENT(get)) THEN
          get(:) = state(id)%ptr(:)
       END IF
       RETURN ! nothing to do
    END IF

    njump(0:p_nprocs-1) = 0
    IF (lfirst_jump(id)) THEN
       ! njump(0) = 0 !!!
       ! JUMP AHEAD TO INITIAL POSITION ...
       DO np=1, p_nprocs-1 
          njump(np) = SUM(n_chunk(0:np-1))
       END DO
       lfirst_jump(id) = .FALSE.  ! ... ONLY ONCE
    ELSE
       ! JUMP AHEAD TO NEXT POSTION
       njump(0:p_nprocs-1) = n_size - n_chunk(0:p_nprocs-1)
    END IF

    IF (njump(p_pe) == 0) THEN
       ! NO JUMP AHEAD
       CALL rnd_number(id, harvest, get)
    ELSE
       ! JUMP AHEAD
       CALL rnd_jump(status, id, njump(p_pe))

       IF (status /= 0) &
            CALL error_bi('rnd_jump reported an error', substr)
       CALL rnd_number(id, harvest, get)
    ENDIF

  END SUBROUTINE rnd_number_sja
! -----------------------------------------------------------------------

! =======================================================================
! ### BMIL ENTRY POINTS
! =======================================================================

! -----------------------------------------------------------------------
  SUBROUTINE main_rnd_init_coupling

    USE messy_main_mpi_bi,             ONLY: p_nprocs
    USE messy_main_tools,              ONLY: str
    USE messy_main_channel_dimensions, ONLY: get_dimension_info, new_dimension
    USE messy_main_channel_repr,       ONLY: get_representation_info &
                                           , new_representation, AUTO   &
                                           , set_representation_decomp  &
                                           , IRANK, PIOTYPE_COL
    USE messy_main_channel_error_bi,   ONLY: channel_halt
    USE messy_main_channel_bi,         ONLY: DC_AG
    USE messy_main_channel,            ONLY: new_channel, new_channel_object &
                                           , new_attribute, get_channel_info &
                                           , AF_RST_CMP
    USE messy_main_timer,              ONLY: lstart

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_rnd_init_coupling'
    INTEGER :: id
    CHARACTER(LEN=2) :: istr = ''
    INTEGER :: status, status_d
    INTEGER :: reprid, dimid
    INTEGER :: dimlen
    ! PARALLEL DECOMPOSITION
    INTEGER                          :: nseg = 0
    INTEGER, DIMENSION(:,:), POINTER :: start => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()
    INTEGER :: DIMID_NCPUS

    ! INTITIALISE
    ! SET TRUE AT VERY FIRST TIME STEP (AND KEEP IT IN THIS STATE UNTIL
    ! THE FIRST JUMP-AHEAD) ... NOTE THAT lfirst_jump(:) = lstart
    ! WOULD BE WRONG ( = .FALSE. FROM THE 2nd TIME STEP ON, EVEN WITHOUT
    ! JUMP-AHEAD ...)
    IF (lstart) THEN
       lfirst_jump(:) = .TRUE.
    END IF
    
    ! shared between QTIMER AND RND ...
    CALL get_dimension_info(status,'NCPUS',DIMID_NCPUS,dimlen)
    IF (status == 0) THEN
       IF (dimlen /= p_nprocs) &
            CALL error_bi('dimension NCPUS exists with wrong length',substr)
    ELSE
       CALL new_dimension(status, DIMID_NCPUS, 'NCPUS', p_nprocs)
       CALL channel_halt(substr, status)
    END IF

    id_loop: DO id=1, id_max

       IF (nstate(id) == 0) CYCLE

       IF ( (rnd_method(id) < 0) .OR. &
            (rnd_method(id) > RND_MAX_METHOD ) ) THEN
          istr = str(id)
          CALL error_bi('unknown method id = '//istr, substr)
       END IF

       ! check, if representation exists
       CALL get_representation_info(status &
            , TRIM(RND_REPR_NAME(rnd_method(id))), reprid)
       repr: IF (status /= 0) THEN

          ! create new representation ...
          CALL get_dimension_info(status_d &
               , TRIM(RND_DIM_NAME(rnd_method(id))), dimid, dimlen)
          ! ... but check dimension first ...
          dim: IF (status_d /= 0) THEN
             ! create new dimension
             dimlen = nstate(id)
             CALL new_dimension(status &
                  , dimid, TRIM(RND_DIM_NAME(rnd_method(id))), dimlen)
             CALL channel_halt(substr, status)
          ELSE
             ! check length
             IF (dimlen /= nstate(id)) THEN
                CALL error_bi( &
                     'DIMENSION '//TRIM(RND_DIM_NAME(rnd_method(id)))//&
                     &' exists with wrong length', substr)
             END IF
          END IF dim

          ! create new representation
          CALL new_representation(status, REPRID &
               , TRIM(RND_REPR_NAME(rnd_method(id)))            &
               , rank = 2, link = 'x--x', dctype = DC_AG        &
               , dimension_ids = (/ dimid, DIMID_NCPUS /)       &
               , ldimlen       = (/ AUTO,  1 /)                 &
               , axis = 'N--N'                                  &
               )
          CALL channel_halt(substr, status)

          ! parallel I/O ...
          nseg = 1
          ALLOCATE(start(nseg,IRANK))
          ALLOCATE(cnt(nseg,IRANK))
          ALLOCATE(meml(nseg,IRANK))
          ALLOCATE(memu(nseg,IRANK))

          start(:,:) = 1
          cnt(:,:) = 1
          meml(:,:) = 1
          memu(:,:) = 1
          
          start(:,1) = 1
          cnt(:,1)   = dimlen
          meml(:,1)  = 1
          memu(:,1)  = dimlen

          start(:,4) = 1
          cnt(:,4)   = 1 ! p_nprocs
          meml(:,4)  = 1
          memu(:,4)  = 1 ! p_nprocs

          CALL set_representation_decomp(status, reprid &
               , start, cnt, memu, meml, .FALSE., piotype=PIOTYPE_COL)
          CALL channel_halt(substr, status)

          DEALLOCATE(start) ; NULLIFY(start)
          DEALLOCATE(cnt)   ; NULLIFY(cnt)
          DEALLOCATE(meml)  ; NULLIFY(meml)
          DEALLOCATE(memu)  ; NULLIFY(memu)
          ! ... parallel I/O

          !
       ENDIF repr

       ! CREATE NEW CHANNEL
       CALL get_channel_info(status, modstr)
       IF (status /= 0) THEN
          CALL new_channel(status, modstr, lrestreq=.TRUE.)
          CALL channel_halt(substr, status)
       END IF

       ! CREATE NEW CHANNEL OBJECT
       ! STATUS VECTOR OF RANDOM NUMBER GENERATOR
       istr = str(id,'(i2.2)')
       CALL new_channel_object(status, modstr, 'RAND_STATE_ID'//istr &
            , p2=rstate(id)%ptr, reprid=reprid)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'RAND_STATE_ID'//istr &
            , 'long_name', c='state vector of pseudo random number generator')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'RAND_STATE_ID'//istr &
            , 'parallel_method', i=MP(id), iflag=AF_RST_CMP)
       CALL channel_halt(substr, status)

    END DO id_loop

  END SUBROUTINE main_rnd_init_coupling
! -----------------------------------------------------------------------

! -----------------------------------------------------------------------
  SUBROUTINE main_rnd_init_tracer  

    USE messy_main_timer,         ONLY: lresume
    USE messy_main_tools,         ONLY: str

    IMPLICIT NONE
    INTRINSIC :: INT, SUM, ABS

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_rnd_init_tracer'
    INTEGER :: id
    CHARACTER(LEN=2) :: istr = ''
    INTEGER, DIMENSION(:), ALLOCATABLE :: tmp_state

    IF (.NOT. lresume) RETURN

    id_loop: DO id=1, id_max

       IF (nstate(id) == 0) CYCLE

       ! For each state vector a channel object (rstate) has been created
       ! in main_rnd_init_coupling (see above).
       ! The channel objects (rstate) are read from restart files and
       ! need to be transferred back to the integer state vector.
       ! 
       IF (ASSOCIATED(state(id)%ptr)) THEN
          ! check size
          IF ( SIZE(state(id)%ptr(:)) == SIZE(rstate(id)%ptr(:,1)) ) THEN

             ! Some care has to be taken, if the user modifies the
             ! random number stream(s) after restart, e.g., by
             ! 1) adding new stream(s)
             !   (These are added with new ID (> as exisiting IDs), the
             !    lrestreq flag will terminate the model during import of
             !    channel objects from restart file, unless the user
             !    forced to ignore this (channel.nml). In this case, the
             !    state vector rstate contains zero(s).)
             ! 2) switching off stream(s)
             !    (This will cause a shift of IDs and the channel objects
             !    (variable names and dimension lengths) in the restart files
             !    do not fit anymore. An error will occur when the channel
             !    objects are read from the restart files. If the user is
             !    smart enough, he/she will remove the restart file and force
             !    channel to ignore the lrestreq flag. In this case, the
             !    state vector rstate contains zero(s).)
             ! 3) changing the random number generator (method)
             !    (This is very similar to item 2) above ...)
             !
             ! In summary, rstate == 0. needs to be detected.
             !
             ALLOCATE(tmp_state(SIZE(rstate(id)%ptr(:,1))))
             tmp_state(:) = INT(rstate(id)%ptr(:,1))

             IF ( SUM(ABS(tmp_state(:))) == 0 ) THEN
                istr = str(id)
                CALL warning_bi('state vector is zero after restart (id = '//&
                     &istr//')',substr)
                ! IN THIS CASE DO NOTHING!
                ! rnd_init_bi MUST HAVE BEEN CALLED ALREADY, OTHERWISE
                ! THE STREAM WOULD NOT BE ACTIVE. THE CALL TO rnd_init_bi
                ! WAS EITHER WITH EXPLICIT SEED OR WITH INTERNAL DEFAULF SEED.
                ! IN BOTH CASES, THE STATE VECTOR IS NON-ZERO.
             ELSE
                ! RE-INITIALISE FROM RESTART FILE
                ! (OVERWRITE INITIAL SEED OF INI-PHASE)
                state(id)%ptr(:) = INT(rstate(id)%ptr(:,1))
             END IF
             
             DEALLOCATE(tmp_state)

          ELSE
             ! Since rstate is created from state (main_rnd_init_coupling)
             ! this might never be reached. If the user modified the
             ! state vector length (e.g. by changing the random number generator
             ! method) after a restart, an error occurs already when
             ! the channel objects are read from the restart files.
             istr = str(id)
             CALL error_bi('incompatible state vector length (id = '//&
                  &istr//' ) ',substr)
          ENDIF
       END IF

    END DO id_loop

  END SUBROUTINE main_rnd_init_tracer
! -----------------------------------------------------------------------

! -----------------------------------------------------------------------
  SUBROUTINE main_rnd_write_output 

    IMPLICIT NONE
    INTRINSIC :: REAL

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_rnd_write_output'
    INTEGER :: id

    id_loop: DO id=1, id_max

       IF (nstate(id) == 0) CYCLE   

       rstate(id)%ptr(:,1) = REAL(state(id)%ptr(:), dp)

    END DO id_loop

  END SUBROUTINE main_rnd_write_output
! -----------------------------------------------------------------------

END MODULE messy_main_rnd_bi
