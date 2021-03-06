! **********************************************************************
MODULE messy_main_qtimer_bi
! **********************************************************************

  USE messy_main_timer,          ONLY: L_TRIGGER_RESTART
  USE messy_main_channel_repr,   ONLY: REPR_UNDEF
  USE messy_main_constants_mem,  ONLY: STRLEN_MEDIUM
  USE messy_main_qtimer

  IMPLICIT NONE
  PRIVATE
  INTRINSIC :: NULL

  ! CLOCKS [s]
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE :: QC_TM1
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE :: QC_TSTEP
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE :: QC_TOTAL

  ! METHOD
  INTEGER, PARAMETER :: M_MAX = 1
  INTEGER, PARAMETER :: M_AVE = 2
  INTEGER, PARAMETER :: M_SUM = 3
  INTEGER, PARAMETER :: M_IND = 4

  ! CPL-NAMELIST PARAMETERS
  ! METHOD FOR PARALLEL PROCESING ('max', 'ave', 'sum', 'ind')
  CHARACTER(LEN=3) :: QMETHOD = 'ind'   ! 'max', 'ave', 'sum', 'ind'
  ! DIAGNOSTIC OUTPUT ?
  LOGICAL          :: L_DIAG = .FALSE.
  ! DEEP THOUGHT ?
  LOGICAL          :: L_DEEP_THOUGHT = .TRUE.

  ! GLOBAL VARIABLES
  LOGICAL,  SAVE :: LQTIMER = .FALSE.   ! USE QTIMER AT ALL?
  INTEGER,  SAVE :: QCTYPE  = 0         ! WHICH CLOCK TYPE
  INTEGER , SAVE :: PMETHOD = 0         ! WHICH PARALLEL METHOD
  REAL(DP), SAVE :: QMAXTIME = 0.0_DP   ! MAXIMUM QTIME [s]

  INTEGER, SAVE  :: REPR_QTIMER   = REPR_UNDEF
  REAL(DP), DIMENSION(:,:,:,:,:,:), POINTER, SAVE :: clock_mem => NULL()

  ! FOR ONLINE TIME MEASURMENTS
  INTEGER, PUBLIC, PARAMETER :: MEAS_ON  = 1
  INTEGER, PUBLIC, PARAMETER :: MEAS_OFF = -1
  INTEGER, PUBLIC, PARAMETER :: MEAS_ADD = 0

  TYPE t_tmeas
     CHARACTER(LEN=STRLEN_MEDIUM) :: name       = ' '
     CHARACTER(LEN=STRLEN_MEDIUM) :: submodname = ' '
     INTEGER                      :: istart     = -99
  END type t_tmeas

  INTEGER, PARAMETER                            :: NMAXMEAS = 200
  TYPE(t_tmeas), DIMENSION(NMAXMEAS),      SAVE :: timemeasure

  INTEGER,                                 SAVE :: maxmeas = 0
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER, SAVE :: measure_mem => NULL()

  ! MAIN ENTRY POINTS
  PUBLIC :: main_qtimer_setup
  PUBLIC :: main_qtimer_init_memory
  PUBLIC :: main_qtimer_tendency_reset
  PUBLIC :: main_qtimer_global_end
  PUBLIC :: main_qtimer_free_memory
  !
  PUBLIC :: main_qtimer_measure_init
  PUBLIC :: main_qtimer_measure

CONTAINS

  ! -----------------------------------------------------------------------
  ! PUBLIC SUBROUTINES
  ! -----------------------------------------------------------------------

  ! -----------------------------------------------------------------------
  SUBROUTINE main_qtimer_setup

    USE messy_main_mpi_bi,             ONLY: p_parallel_io, p_pe &
                                           , p_bcast, p_nprocs, p_io
    USE messy_main_blather_bi,         ONLY: start_message_bi &
                                           , end_message_bi, error_bi
    USE messy_main_tools,              ONLY: find_next_free_unit
    USE messy_main_channel_error_bi,   ONLY: channel_halt
    USE messy_main_channel_bi,         ONLY: DC_AG
    USE messy_main_channel_dimensions, ONLY: new_dimension, DIMID_UNDEF &
                                           , get_dimension_info
    USE messy_main_channel_repr,       ONLY: new_representation         &
                                           , set_representation_decomp  &
                                           , IRANK, PIOTYPE_COL
#ifdef ICON
    USE messy_main_channel_mem,        ONLY: dom_unbound
#endif

    IMPLICIT NONE
    INTRINSIC :: REAL, TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_qtimer_setup'
    INTEGER :: status
    INTEGER :: iou
    INTEGER :: i
    INTEGER :: hh(NCLOCKS)
    INTEGER :: mm(NCLOCKS)
    INTEGER :: ss(NCLOCKS)
    INTEGER :: ms(NCLOCKS)
    INTEGER :: DIMID_NCPUS   = DIMID_UNDEF
    INTEGER :: dimlen
    ! PARALLEL DECOMPOSITION
    INTEGER                          :: nseg = 0
    INTEGER, DIMENSION(:,:), POINTER :: start => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()

    CALL start_message_bi(modstr,'SETUP QTIMER',substr)

    ! READ CTRL-NAMELIST
    ! INITIALIZE CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL main_qtimer_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF
    ! BROADCAST RESULTS
    CALL p_bcast(QTIME, p_io)
    CALL p_bcast(QCLOCK, p_io)
    CALL p_bcast(QFRAC, p_io)

    ! SWITCH QTIMER
    LQTIMER = .NOT. &
         ( (QTIME(1) == 0) .AND. (QTIME(2) == 0) .AND. (QTIME(3) == 0) )

    ! INITIALIZE CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL main_qtimer_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF
    ! BROADCAST RESULTS
    CALL p_bcast(L_DIAG, p_io)
    CALL p_bcast(QMETHOD, p_io)
    CALL p_bcast(L_DEEP_THOUGHT, p_io)

    ! shared between QTIMER AND RND ...
    CALL get_dimension_info(status,'NCPUS', DIMID_NCPUS, dimlen)
    IF (status == 0) THEN
       IF (dimlen /= p_nprocs) &
            CALL error_bi('dimension NCPUS exists with wrong length',substr)
    ELSE
       CALL new_dimension(status, DIMID_NCPUS, 'NCPUS', p_nprocs)
       CALL channel_halt(substr, status)
    END IF

    ! NEW REPRESENTATION
    ! NOTE: DC_AG requires p_pe at 4th rank ...
    CALL new_representation(status, REPR_QTIMER, 'REPR_QTIMER' &
         , rank = 1, link = '---x', dctype = DC_AG             &
         , dimension_ids = (/ DIMID_NCPUS /)                   &
         , ldimlen       = (/ 1 /)                             &
         , axis          = '---N'                              &
#ifdef ICON
         , dom_id        = dom_unbound                         &
#endif
         )
    CALL channel_halt(substr, status)

    nseg = 1
    ALLOCATE(start(nseg,IRANK))
    ALLOCATE(cnt(nseg,IRANK))
    ALLOCATE(meml(nseg,IRANK))
    ALLOCATE(memu(nseg,IRANK))

    start(:,:) = 1
    cnt(:,:)   = 1
    meml(:,:)  = 1
    memu(:,:)  = 1

    start(:,1) = 1
    cnt(:,4)   = 1 ! p_nprocs
    meml(:,1)  = 1
    memu(:,4)  = 1 ! p_nprocs

    CALL set_representation_decomp(status, REPR_QTIMER &
         , start, cnt, memu, meml, piotype=PIOTYPE_COL)
    CALL channel_halt(substr, status)

    DEALLOCATE(start) ; NULLIFY(start)
    DEALLOCATE(cnt)   ; NULLIFY(cnt)
    DEALLOCATE(meml)  ; NULLIFY(meml)
    DEALLOCATE(memu)  ; NULLIFY(memu)

    ! SET UP CLOCK MEMORY
    ALLOCATE(QC_TM1  (NCLOCKS))
    QC_TM1(:)   = 0.0_DP
    ALLOCATE(QC_TSTEP(NCLOCKS))
    QC_TSTEP(:) = 0.0_DP
    ALLOCATE(QC_TOTAL(NCLOCKS))
    QC_TOTAL(:) = 0.0_DP

    ! INITIALIZE CLOCK FOR TM1 AT START ON LOCAL PE
    CALL qtimer_read_clock(QC_TM1(:))

    ! CALCULATE MAXIMUM QTIME
    QMAXTIME = REAL(QTIME(1), DP) * 3600.0_DP  & ! hours -> seconds
         + REAL(QTIME(2), DP) * 60.0_DP        & ! minutes -> seconds
         + REAL(QTIME(3), DP)                    ! seconds

    IF (p_parallel_io) THEN
       ! CTRL
       WRITE(*,'(a29,2(i2.2,a1),i2.2)') 'LENGTH OF QUEUE [hh:mm:ss]: ' &
            , qtime(1),':',qtime(2),':',qtime(3)
       IF (.NOT. LQTIMER) THEN
          WRITE(*,*) ' *** QTIMER SWITCHED OFF ***'
       ELSE
          WRITE(*,'(a29,f20.12)') '                       [s] : ', QMAXTIME
          WRITE(*,'(a29,a4)')     ' CLOCK LIMIT IS FOR        : ', QCLOCK
          WRITE(*,'(a29,f6.4)')   ' USABLE FRACTION IS        : ', QFRAC
          WRITE(*,'(a29,f20.12)') ' USABLE TIME IS    ->  [s] : ', &
               QMAXTIME*QFRAC
          ! CPL
          WRITE(*,'(a29,L1)')     ' DIAGNOSTIC LOG-OUTPUT IS  : ', L_DIAG
          WRITE(*,'(a29,a4)')     ' PARALLEL MODE IS          : ', QMETHOD
       END IF
    END IF
    ! CORRECT FOR FRACTION
    QMAXTIME = QMAXTIME * QFRAC

    IF (QMAXTIME < 0.0_DP) &
         CALL error_bi('ERROR: MAXIMUM QUEUE TIME < 0 s !',substr)

    ! SET CLOCK TYPE
    SELECT CASE(TRIM(QCLOCK))
    CASE('wall')
       QCTYPE = C_WALL
    CASE('cpu')
       QCTYPE = C_CPU
    CASE('user')
       QCTYPE = C_USER
    CASE('sys')
       QCTYPE = C_SYS
    CASE DEFAULT
       ! NEVER REACHED
    END SELECT

    ! SET PARALLEL MODE
    SELECT CASE(TRIM(QMETHOD))
    CASE('max')
       PMETHOD = M_MAX
    CASE('sum')
       PMETHOD = M_SUM
    CASE('ave')
       PMETHOD = M_AVE
    CASE('ind')
       PMETHOD = M_IND
    CASE DEFAULT
       ! NEVER REACHED
    END SELECT

    IF (L_DIAG) THEN
       CALL sec2hhmmssms(QC_TM1(:), hh(:), mm(:), ss(:), ms(:))
       !
       WRITE(*,'(a56)') &
            '========================================================'
       WRITE(*,'(a56)') &
            '                       START TIME                       '
       WRITE(*,'(a56)') &
            '     WALL         CPU          USER         SYS         '
       WRITE(*,'(a56)') &
            '--------------------------------------------------------'
       WRITE(*,'(a56)') &
            'P_PE hh:mm:se.ms  hh:mm:se.ms  hh:mm:se.ms  hh:mm:se.ms '
       WRITE(*,'(a56)') &
            '--------------------------------------------------------'
       WRITE(*,'(i4,4(1x,i2.2,a1,i2.2,a1,i2.2,a1,i3.3))') p_pe, &
            (hh(i),':', mm(i),':',ss(i),'.',ms(i), i=1,4)
       WRITE(*,'(a56)') &
            '========================================================'
    END IF

    CALL end_message_bi(modstr,'SETUP QTIMER',substr)

  END SUBROUTINE main_qtimer_setup
  ! -----------------------------------------------------------------------

  ! -----------------------------------------------------------------------
  SUBROUTINE main_qtimer_init_memory

    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         , new_attribute

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_qtimer_init_memory'
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: mem => NULL()
    INTEGER :: status
    INTEGER :: i

    ! Note: 2nd '1' can be replaced for additonal set of
    !       channel objects, eg. the formerly knowm _ACC accumulated fields
    ALLOCATE(clock_mem(NCLOCKS, 1, 1, 1, 1, 1))
    clock_mem(:,:,:,:,:,:) = 0.0_dp

    ! --------------------------------------------------------

    CALL new_channel(status, modstr)
    CALL channel_halt(substr, status)

    ! --------------------------------------------------------

    mem => clock_mem(1,1,:,:,:,:)
    CALL new_channel_object(status, modstr, 'WALL' &
         , mem=mem, reprid=REPR_QTIMER)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'WALL' &
         , 'long_name', c='wall clock time' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'WALL' &
         , 'units', c='s' )
    CALL channel_halt(substr, status)

    mem => clock_mem(2,1,:,:,:,:)
    CALL new_channel_object(status, modstr, 'CPU' &
         , mem=mem, reprid=REPR_QTIMER)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'CPU' &
         , 'long_name', c='cpu time' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'CPU' &
         , 'units', c='s' )
    CALL channel_halt(substr, status)

    mem => clock_mem(3,1,:,:,:,:)
    CALL new_channel_object(status, modstr, 'USER' &
         , mem=mem, reprid=REPR_QTIMER)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'USER' &
         , 'long_name', c='user time' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'USER' &
         , 'units', c='s' )
    CALL channel_halt(substr, status)

    mem => clock_mem(4,1,:,:,:,:)
    CALL new_channel_object(status, modstr, 'SYSTEM' &
         , mem=mem, reprid=REPR_QTIMER)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'SYSTEM' &
         , 'long_name', c='system time' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'SYSTEM' &
         , 'units', c='s' )
    CALL channel_halt(substr, status)

    mem => clock_mem(5,1,:,:,:,:)
    CALL new_channel_object(status, modstr, 'REAL' &
         , mem=mem, reprid=REPR_QTIMER)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'REAL' &
         , 'long_name', c='real time' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'REAL' &
          , 'units', c='day since 2015-01-01 00:00:00' )

    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'REAL' &
         , 'calendar', c='julian' )
    CALL channel_halt(substr, status)

    ! --------------------------------------------------------

    clock_mem(:,1,1,1,1,1) = QC_TM1(:)

    ! --------------------------------------------------------

    ! INIT MEMORY FOR TIME MEASURE
    ALLOCATE(measure_mem(maxmeas, 1, 1, 1, 1))
    measure_mem(:,:,:,:,:) = 0.0_dp

    DO i = 1, maxmeas
       mem => measure_mem(i,:,:,:,:)
       CALL new_channel_object(status, modstr, TRIM(timemeasure(i)%name) &
         , mem=mem, reprid=REPR_QTIMER &
         )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, TRIM(timemeasure(i)%name) &
            , 'long_name', c='measured time'  &
            )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, TRIM(timemeasure(i)%name) &
            , 'submodel', c=TRIM(timemeasure(i)%submodname)  &
            )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, TRIM(timemeasure(i)%name) &
            , 'units', c='s'  &
            )
       CALL channel_halt(substr, status)
    END DO

  END SUBROUTINE main_qtimer_init_memory
  ! -----------------------------------------------------------------------

  ! -----------------------------------------------------------------------
  SUBROUTINE main_qtimer_tendency_reset

    IMPLICIT NONE

    measure_mem(:,:,:,:,:) = 0._dp

  END SUBROUTINE main_qtimer_tendency_reset
  ! -----------------------------------------------------------------------

  ! -----------------------------------------------------------------------
  SUBROUTINE main_qtimer_global_end

    USE messy_main_mpi_bi,   ONLY: p_parallel_io, p_pe, p_nprocs &
                                 , p_sum, p_lor, p_max
    USE messy_main_timer,    ONLY: current_time_step
    USE messy_main_data_bi,  ONLY: L_IS_CHILD

    IMPLICIT NONE
    INTRINSIC :: REAL

    ! LOCAL
    INTEGER  :: hh(NCLOCKS, 2)  ! 1: TOTAL TIME, 2: LAST STEP
    INTEGER  :: mm(NCLOCKS, 2)  ! 1: TOTAL TIME, 2: LAST STEP
    INTEGER  :: ss(NCLOCKS, 2)  ! 1: TOTAL TIME, 2: LAST STEP
    INTEGER  :: ms(NCLOCKS, 2)  ! 1: TOTAL TIME, 2: LAST STEP
    INTEGER  :: i
    REAL(DP) :: tcheck
    LOGICAL  :: ltmp

    IF (L_DEEP_THOUGHT) THEN
       IF (p_parallel_io) THEN
          CALL deep_thought(current_time_step)
       END IF
    END IF

    IF (.NOT. LQTIMER) RETURN

    ! SET TOTAL ELAPSED TIME ON LOCAL PE
    CALL qtimer_read_clock(QC_TOTAL(:))

    ! CALCULATE TIME ELAPSED SINCE LAST STEP
    QC_TSTEP(:) = QC_TOTAL(:) - QC_TM1(:)

    ! SHIFT CLOCK TO TM1 FOR NEXT STEP
    QC_TM1(:) = QC_TOTAL(:)

    ! CALCULATE VALUE FOR PARALLEL PROCESSING
    ! Note: Here, once the all-to-all MPI communication is required!
    SELECT CASE(PMETHOD)
    CASE(M_MAX)
       tcheck = p_max(QC_TOTAL(QCTYPE))
    CASE(M_SUM)
       tcheck = p_sum(QC_TOTAL(QCTYPE))
    CASE(M_AVE)
       tcheck = p_sum(QC_TOTAL(QCTYPE)) / REAL(p_nprocs, DP)
    CASE(M_IND)
       tcheck = QC_TOTAL(QCTYPE)
    CASE DEFAULT
       ! NEVER REACHED
    END SELECT

    ! CHECK IF RESTART NEEDS TO BE TRIGGERED
    ! If the model is an MMD-client, it is not possible to trigger
    ! a restart, because the MMD-server model might still be waiting within
    ! MMD.
    !L_TRIGGER_RESTART = (tcheck >= QMAXTIME)
    L_TRIGGER_RESTART = (tcheck >= QMAXTIME) .AND. .NOT. L_IS_CHILD
    ! ... special for this method: make sure that all tasks trigger
    !     at the same time.
    IF (PMETHOD == M_IND) THEN
       ltmp = (tcheck >= QMAXTIME) ! local
       ! any taks and not child
       L_TRIGGER_RESTART = p_lor(ltmp) .AND. .NOT. L_IS_CHILD
    END IF

    ! DIAGNOSTIC OUTPUT
    IF (L_DIAG) THEN
       CALL sec2hhmmssms(QC_TOTAL(:), hh(:,1), mm(:,1), ss(:,1), ms(:,1))
       CALL sec2hhmmssms(QC_TSTEP(:), hh(:,2), mm(:,2), ss(:,2), ms(:,2))
       WRITE(*,'(a56)') &
            '========================================================'
       WRITE(*,'(a32,i10,a14)') &
            '                  ELAPSED TIME (',current_time_step, &
            ')             '
       WRITE(*,'(a56)') &
            '     WALL         CPU          USER         SYS         '
       WRITE(*,'(a56)') &
            '--------------------------------------------------------'
       WRITE(*,'(a56)') &
            'P_PE hh:mm:se.ms  hh:mm:se.ms  hh:mm:se.ms  hh:mm:se.ms '
       WRITE(*,'(a56)') &
            '--------------------------------------------------------'
       WRITE(*,'(i4,4(1x,i2.2,a1,i2.2,a1,i2.2,a1,i3.3))') p_pe, &
            (hh(i,1),':', mm(i,1),':',ss(i,1),'.',ms(i,1), i=1,4)
       WRITE(*,'(4x,4(1x,i2.2,a1,i2.2,a1,i2.2,a1,i3.3))') &
            (hh(i,2),':', mm(i,2),':',ss(i,2),'.',ms(i,2), i=1,4)

       WRITE(*,'(a56)') &
            '========================================================'
       WRITE(*,'(a40,f6.4)') ' FRACTION OF USABLE QUEUE TIME REACHED: ', &
            tcheck/QMAXTIME
       WRITE(*,'(a56)') &
            '========================================================'
    END IF

    IF (L_TRIGGER_RESTART) THEN
       IF (p_parallel_io) THEN
          WRITE(*,*) '###############################################'
          WRITE(*,*) '###                 QTIMER                  ###'
          WRITE(*,*) '### !!! USABLE QUEUE TIME LIMIT REACHED !!! ###'
          WRITE(*,*) '###         >>> TRIGGER RESTART <<<         ###'
          WRITE(*,*) '###############################################'
       END IF
    END IF

    ! ALWAYS SAVE p_pe-specific values on local task.
    ! This implies that also the statistics (ave, std, min, max, ...)
    ! are internally correct for each p_pe.
    ! For correct output to netCDF, however, p_io requires the
    ! information of all PEs ...
    clock_mem(:,1,1,1,1,1) = QC_TSTEP(:)
    ! ... This is achieved during output via MPI_ALLGATHER,
    ! (see representation with DC_AG).

  END SUBROUTINE main_qtimer_global_end
  ! -----------------------------------------------------------------------

  ! -----------------------------------------------------------------------
  SUBROUTINE main_qtimer_free_memory

    IMPLICIT NONE

    DEALLOCATE(QC_TM1)
    DEALLOCATE(QC_TOTAL)
    DEALLOCATE(QC_TSTEP)
    IF (ASSOCIATED(clock_mem)) THEN
       DEALLOCATE(clock_mem) ; NULLIFY(clock_mem)
    END IF
    IF (ASSOCIATED(measure_mem)) THEN
       DEALLOCATE(measure_mem) ; NULLIFY(measure_mem)
    END IF

  END SUBROUTINE main_qtimer_free_memory
  ! -----------------------------------------------------------------------

  ! -----------------------------------------------------------------------
  ! PRIVATE SUBROUTINES
  ! -----------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_qtimer_read_nml_cpl(status, iou)

    ! MODULE ROUTINE (SMIL)
    !
    ! READ NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! Author: Patrick Joeckel, MPICH, May 2005

    USE messy_main_tools,   ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    NAMELIST /CPL/ L_DIAG, QMETHOD, L_DEEP_THOUGHT

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_qtimer_read_nml_cpl'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! CHECK NAMELIST
    SELECT CASE(TRIM(QMETHOD))
    CASE('max', 'sum', 'ave', 'ind')
       ! OK
    CASE DEFAULT
       WRITE(*,*) &
            '*** ERROR: UNKNOWN METHOD FOR PARALLEL PROCESSING: '//QMETHOD
       status = 2
       RETURN
    END SELECT

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE main_qtimer_read_nml_cpl
  ! -------------------------------------------------------------------

  ! --------------------------------------------------------------------
  SUBROUTINE main_qtimer_measure_init(handle, name, submodname)

    USE messy_main_blather_bi,  ONLY: error_bi

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)  :: name
    CHARACTER(LEN=*), INTENT(IN)  :: submodname
    INTEGER,          INTENT(OUT) :: handle

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: SUBSTR = "main_qtimer_measure_init"
    INTEGER                     :: i

    IF (maxmeas > NMAXMEAS) &
         CALL error_bi(' MAX number of measure points exceeded', substr)

    DO i = 1, maxmeas
       IF (TRIM(ADJUSTL(name)) == TRIM(timemeasure(i)%name)) &
            CALL error_bi(' Name of time measure point exists already', substr)
    END DO
    maxmeas = maxmeas + 1

    timemeasure(maxmeas)%name       = ADJUSTL(name)
    timemeasure(maxmeas)%submodname = ADJUSTL(submodname)

    handle = maxmeas

  END SUBROUTINE main_qtimer_measure_init
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  SUBROUTINE main_qtimer_measure(handle,mode)

    USE messy_main_blather_bi,  ONLY: error_bi, warning_bi

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: handle
    INTEGER, INTENT(IN) :: mode

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_qtimer_measure'
    INTEGER  :: ir, im, icount
    REAL(dp) :: timediff

    SELECT CASE (mode)

    CASE (MEAS_ON)
       CALL SYSTEM_CLOCK( COUNT=timemeasure(handle)%istart &
                        , COUNT_RATE=ir, COUNT_MAX=im )
    CASE (MEAS_OFF, MEAS_ADD)
       IF (timemeasure(handle)%istart == -99) THEN
          CALL WARNING_BI('no start time available',substr)
          RETURN
       END IF

       CALL SYSTEM_CLOCK( COUNT=icount, COUNT_RATE=ir, COUNT_MAX=im )

       IF ( ir == 0 ) THEN
          CALL WARNING_BI('system clock is not present',substr)
          RETURN
       ELSE
          ! convert the clock counts to seconds
          IF ( icount >= timemeasure(handle)%istart ) THEN
             timediff = (REAL(icount - timemeasure(handle)%istart,dp)) &
                  / REAL(ir,dp)
          ELSE
             timediff = REAL(im- (timemeasure(handle)%istart-icount ),dp) &
                  / REAL(ir,dp)
          ENDIF
          timemeasure(handle)%istart = -99

          ! Store value in the appropriate entry:
          SELECT CASE (mode)
          CASE (MEAS_OFF)
             measure_mem(handle, 1,1,1, 1) = timediff
          CASE(MEAS_ADD)
             measure_mem(handle, 1,1,1, 1) = &
                  measure_mem(handle, 1,1,1, 1)  + timediff
          END SELECT
       ENDIF

    CASE DEFAULT
       CALL error_bi('unknown mode',substr)

    END SELECT

  END SUBROUTINE main_qtimer_measure
  ! --------------------------------------------------------------------

! **********************************************************************
END MODULE messy_main_qtimer_bi
! **********************************************************************
