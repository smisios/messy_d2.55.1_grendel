! ****************************************************************************
MODULE messy_main_timer_bi
! ****************************************************************************

!#define MESSYTIMERTEST

  !  Author: Astrid Kerkweg,  IPA, Uni-Mainz, 2008
  !
  ! MAIN-TIMER is based on the date and time management of ECHAM5
  ! For the use within MESSy the basic subroutines have been heavily
  ! simplified.
  ! MAIN_TIMER comprises 3 core modules:
  ! - messy_main_timer:
  !                       This module defines the basic date and
  !                       variables needed within a model run to coordinate
  !                       the start, stop, restarts and defined EVENTS.
  ! - messy_main_timer_manager:
  !                       It determines the "heart beat" of the  model.
  !                       Based on a given start date and a time step
  !                       messy_main_timer_manager will always provide
  !                       the current date and time step number.
  ! - messy_main_timer_event:
  !                       This module provides the interface routines to
  !                       determine and manage certain events (e.g. restarts,
  !                       data output or emissions).


  ! BM/MESSy
  USE messy_main_mpi_bi,       ONLY: p_parallel_io
  USE messy_main_blather_bi,   ONLY: error_bi, info_bi, warning_bi &
       , start_message_bi, end_message_bi

  ! MESSy
  USE messy_main_constants_mem, ONLY: dp, OneDay  &
       , STRLEN_ULONG, STRLEN_SHORT

  ! CORE
  USE messy_main_timer_manager
  USE messy_main_timer_event
  USE messy_main_timer

  IMPLICIT NONE
  PRIVATE
  SAVE

  ! PUBLIC attribute required for ECHAM timer adjustment ...
  TYPE(time_manager), PUBLIC, POINTER :: BM_time => NULL() !base model time axis

  TYPE t_BM_time
     !  base model time axis
     TYPE(time_manager), POINTER :: BM_time => NULL()
  END type t_BM_time

  TYPE(t_BM_time), DIMENSION(:), POINTER :: BMclocks => NULL()

#ifdef ICON
  LOGICAL, DIMENSION(:), ALLOCATABLE :: lfirstreset
#endif

  LOGICAL :: lprintnewdate=.FALSE.  ! print new current_date in timer_set_time

  ! 0: annual cycle; 1..12 : perpetual month experiment
  INTEGER, PUBLIC :: nmonth = 0

  ! CALLED FROM CONTROL
  PUBLIC :: main_timer_setup
  PUBLIC :: main_timer_initialize
  PUBLIC :: main_timer_time
  PUBLIC :: main_timer_free_memory
  !
  PUBLIC :: main_timer_set_domain
  ! CALLED FROM BMIL or SMIL
  PUBLIC :: messy_timer_init_manager
  !
  ! PRIVATE
  !PRIVATE :: main_timer_reset_time
  !PRIVATE :: main_timer_set_time
#ifdef COSMO
  PUBLIC :: messy_timer_COSMO_reinit_time
#endif
  ! PRIVATE :: timer_read_namelist
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  ! event_state
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  PUBLIC :: p_bcast_event
  PUBLIC :: event_state                     ! get the state of event
  INTERFACE event_state
     !  get interval
     MODULE PROCEDURE I_event_steps          ! int<-(I:event,I:delta) in steps
     MODULE PROCEDURE R_event_seconds        ! int<-(I:event)  in seconds
     MODULE PROCEDURE R_event_seconds_next   ! int<-(I:event,I:true) future (s)
     !MODULE PROCEDURE I_event_seconds        ! int<-(I:event)  in seconds
     !MODULE PROCEDURE I_event_seconds_next   ! int<-(I:event,I:true) future (s)
     !  evaluate next trigger, check date
     MODULE PROCEDURE L_event_trigger        ! log<-(I:event,I:time_days)
  END INTERFACE

  PUBLIC  :: event_eval        ! int<-(I:event,I:delta) fit delta into interval
  PUBLIC  :: event_print       ! (I:event[,I:format]) print event
  PUBLIC  :: timer_event_init

  ! PUBLIC HELPER ROUTINE
  PUBLIC :: timer_set_rerun_event
  PUBLIC :: timer_get_rerun_event
  PUBLIC :: timer_message
  PUBLIC :: timer_get_time_step

  PUBLIC :: write_date_debug
  ! PRIVATE :: write_date
  ! PRIVATE :: init_BM_manager
  ! PRIVATE :: reset_BM_manager

  TYPE(time_event) :: RERUN_EVENT

  TYPE(io_time_event)  :: IO_RERUN_EV = &
       io_time_event(1, TIME_INC_MONTHS,TRIG_EXACT,0._dp) ! DEFAULT

  INTEGER :: MODEL_START(6) = 0
  INTEGER :: MODEL_STOP(6)  = 0
  TYPE(time_days) :: dummy_date
  INTEGER :: ms

  ! DUMMY INTEGER (see no_cycles)
  INTEGER         :: ccount = 1
#ifdef ICON
  LOGICAL         :: loverwriteICONtimer =.FALSE.
#endif

CONTAINS

  ! -------------------------------------------------------------------------
  SUBROUTINE main_timer_setup(flag, numclock)

    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_bcast, p_io
    USE messy_main_blather_bi, ONLY: error_bi
    USE messy_main_tools,      ONLY: find_next_free_unit
#ifdef COSMO
    USE data_modelconfig,      ONLY: dt
#endif
#ifdef ICON
    USE mo_run_config,         ONLY: dtime
#endif

    IMPLICIT NONE

    INTEGER, INTENT(IN)    :: flag
    INTEGER, INTENT(INOUT) :: numclock

    CHARACTER(LEN=*), PARAMETER :: substr='main_timer_setup'
    INTEGER :: iou
    INTEGER :: status
    LOGICAL :: lsmaller
    LOGICAL :: linit
    CHARACTER(LEN=2) :: nmstr
    CHARACTER(LEN=STRLEN_ULONG) :: txtstart = ''
    CHARACTER(LEN=STRLEN_ULONG) :: txtstop  = ''
    INTEGER        :: ic
#ifdef ICON
    REAL(dp)       :: timestep
#endif

    SELECT CASE(flag)
    CASE(1)
#ifdef ICON
       IF (numclock < 0 ) THEN
          loverwriteICONtimer = .TRUE.
          CALL timer_ICON_get_ndom(numclock)
       END IF
#endif

       CALL timer_alloc_clock(numclock)
       CALL main_timer_set_clock(1)

       ALLOCATE(BMclocks(numclock))
       DO ic = 1,numclock
          ALLOCATE(BMclocks(ic)%BM_time)
       END DO
#ifdef ICON
       ALLOCATE(lfirstreset(numclock))
       lfirstreset = .TRUE.
#endif
       CALL main_timer_set_manager(1)

       ! INITIALIZE CPL
       IF (p_parallel_io) THEN
          iou = find_next_free_unit(100,200)
          CALL main_timer_read_nml_cpl(status, iou)
          IF (status /= 0) CALL error_bi(' ',substr)
       ! BROADCAST RESULTS
#ifndef ECHAM5
          IF (delta_time < 0) THEN
#ifndef ICON
             ! HERE, THE MESSy-TIMER is taking over control:
             ! time step lenght is set in timer.nml
             CALL info_bi('delta_time < 0 ',substr)
             CALL info_bi('delta_time must be given in namelist',substr)
             CALL error_bi('delta_time < 0 ',substr)
#else
             IF (.NOT. loverwriteICONtimer) THEN
                ! HERE, THE BASEMODEL (ICON) is taking over control, if
                ! time step lenght in timer.nml is < 0
                CALL warning_bi(substr, 'delta_time < 0 ')
                CALL warning_bi(substr, 'using dtime from the basemodel')
                delta_time = dtime
             ELSE
                ! HERE, THE MESSy-TIMER is taking over control:
                ! time step lenght is set in timer.nml
                CALL info_bi('delta_time < 0 ',substr)
                CALL info_bi('delta_time must be given in namelist',substr)
                CALL error_bi('delta_time < 0 ',substr)
             END IF
#endif
          ENDIF
#endif
       END IF

       CALL p_bcast(CAL_TYPE     , p_io)
       CALL p_bcast(nmonth       , p_io)
       CALL p_bcast(delta_time   , p_io)
       CALL p_bcast(lresume      , p_io)
       CALL p_bcast(NO_CYCLES    , p_io)
       CALL p_bcast(NO_DAYS      , p_io)
       CALL p_bcast(NO_STEPS     , p_io)
       CALL p_bcast(LABORT       , p_io)
       CALL p_bcast(lprintnewdate, p_io)
#ifdef COSMO
       IF (dt /= delta_time .AND. (dt /= 30._dp)) &
            CALL timer_message(3452, substr)
       dt = delta_time
#endif
#ifdef ICON
       IF (.NOT. loverwriteICONtimer) THEN
          IF (dtime /= delta_time) &
               & CALL timer_message(3453, substr)
       ENDIF
       dtime = delta_time

       timestep = delta_time
       time_step_len = timestep
       CALL main_timer_set_clock(1)
#endif
       DO ic = 1, numclock
          CALL timer_set_lfirst_cycle(.TRUE., ic)
          IF (.NOT. lresume) THEN
             CALL timer_set_lstart(.TRUE., ic)
          ELSE
             CALL timer_set_lstart(.FALSE., ic)
          ENDIF
       END DO

       ! (see ECHAM5 setrad.f90)
       SELECT CASE (nmonth)
       CASE(0)
          CALL info_bi('nmonth=0 --> annual cycle on', substr)
       CASE(1:12)
          WRITE(nmstr,'(i2)') nmonth
          CALL info_bi('nmonth='//nmstr//' --> perpetual month', substr)
       CASE default
          WRITE(nmstr,'(i2)') nmonth
          CALL error_bi(&
               'nmonth='//nmstr//' in timer namelist is not supported' &
               ,substr)
       END SELECT

       IF ((nmonth > 0) .AND. (CAL_TYPE /= CAL_JULIAN)) THEN
          call error_bi('perpetual month setup (nmonth > 0) can only'//&
               &' be used with Julian calender',substr)
       END IF

    CASE(2)

       CALL is_init(start_date, linit)

       IF (.NOT. linit) THEN

          CALL p_bcast(MODEL_START, p_io)

          IF (SUM(MODEL_START) == 0 ) THEN
             CALL error_bi(substr, 'MODEL START DATE NOT GIVEN IN NAMELIST')
          ELSE
             DO ic = 1,nclock
                CALL main_timer_set_clock(ic)
                YEAR_START   = MODEL_START(1)
                MONTH_START  = MODEL_START(2)
                DAY_START    = MODEL_START(3)
                HOUR_START   = MODEL_START(4)
                MINUTE_START = MODEL_START(5)
                SECOND_START = MODEL_START(6)
                MILLISECOND_START = 0
                CALL timer_set_date(status,'start' &
                     , MODEL_START(1), MODEL_START(2) &
                     , MODEL_START(3), MODEL_START(4) &
                     , MODEL_START(5), MODEL_START(6) &
                     , MILLISECOND_START)
                ! set current time first to start_date,
                ! will be overwritten later
                CALL timer_get_date(status,'start' &
                     ,YEAR,MONTH,DAY,HOUR,MINUTE,SECOND,MILLISECOND)
             END DO
             CALL main_timer_set_clock(1)
             CALL timer_set_date(status,'resume' &
                  ,YEAR,MONTH,DAY,HOUR,MINUTE,SECOND,MILLISECOND)
          ENDIF

       ENDIF

       CALL is_init(stop_date, linit)

       IF (.NOT. linit) THEN
          CALL p_bcast(MODEL_STOP , p_io)

          IF (SUM(MODEL_STOP) == 0 ) THEN
             CALL warning_bi(substr, 'MODEL STOP DATE NOT GIVEN IN NAMELIST')
             IF (NO_DAYS > 0) THEN
                CALL copy_date (resume_date, dummy_date)
                IF (lstart) THEN
                   CALL add_date  (no_days, 0._dp, dummy_date)
                ELSE
                   CALL add_date  (no_days, delta_time, dummy_date)
                ENDIF
                CALL warning_bi(substr, 'Using NO_DAYS for model stop.')
             ELSE IF (NO_STEPS > 0) THEN
                CALL copy_date (resume_date, dummy_date)
                CALL add_date  ( 0, REAL(no_steps,dp)*delta_time, dummy_date)
                CALL warning_bi(substr, 'Using NO_STEPS for model stop.')
             ELSE
                CALL error_bi(substr, 'MODEL STOP DATE NOT GIVEN IN NAMELIST')
             ENDIF
             CALL timer_get_date (status, dummy_date&
                  , MODEL_STOP(1), MODEL_STOP(2) &
                  , MODEL_STOP(3), MODEL_STOP(4) &
                  , MODEL_STOP(5), MODEL_STOP(6), ms)
          ENDIF

          CALL timer_set_date(status, 'stop', MODEL_STOP(1), MODEL_STOP(2) &
               , MODEL_STOP(3), MODEL_STOP(4), MODEL_STOP(5), MODEL_STOP(6) &
               , ms)

       ENDIF

#ifdef COSMO
       CALL messy_timer_COSMO_reinit_time
#endif

#ifdef ICON
       ! OVERWRITE namelist dates with ICON settings
       IF (.NOT. loverwriteICONtimer) CALL timer_initfromICONdates
#endif

       CALL if_less(stop_date,start_date,lsmaller,status)
       CALL timer_message(status,substr)

       IF (lsmaller) THEN
          CALL print_date_components(start_date, status, txtstart)
          CALL print_date_components(stop_date,  status, txtstop)
          CALL error_bi( &
               'stop_date smaller than start_date ('//TRIM(txtstart)&
               //' - '//TRIM(txtstop)//')', substr)
       END IF

       ! BROADCAST RERUN EVENT
       CALL p_bcast_event(IO_RERUN_EV, p_io)

#ifdef ICON
       ! OVERWRITE namelist dates with ICON settings
       IF (loverwriteICONtimer) CALL timer_overwriteICONdates
#else
       DO ic = 1, nclock
          CALL  main_timer_set_domain(ic)
          IF (lstart) &
               CALL messy_timer_init_manager(delta_time, INIT_STEP)
       END DO
       CALL  main_timer_set_domain(1)
#endif

    END SELECT

  END SUBROUTINE main_timer_setup
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE messy_timer_init_manager (timestep, step)

    USE messy_main_tools, ONLY: int2str

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! time manager initialization
    REAL(DP), INTENT(IN) :: timestep   ! delta time in seconds
    INTEGER, INTENT(IN)  :: step       ! model time step

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'messy_timer_init_manager'
    INTEGER :: status
    CHARACTER(LEN=STRLEN_ULONG) :: m_text
    INTEGER          :: nstep
    CHARACTER(LEN=2) :: cstr = '  '
#ifdef ICON
    TYPE(time_days)  :: nextdate, resdate
    LOGICAL         :: lfit
    REAL(dp)        :: dt_glob
#endif

    IF (current_clock > 1) THEN
       CALL int2str(cstr,current_clock,'0')
    END IF

    nstep = step

    IF (p_parallel_io) &
         CALL write_date(start_date, 'start_date '//substr)

    CALL manager_init(BM_time,'base model time '//TRIM(cstr) &
         , start_date, timestep, m_text)
    CALL timer_message (0, substr, m_text)

    IF (lresume) THEN
       !
       CALL write_date(resume_date,'Experiment resumed at: ')
       !
#ifdef ICON
       ! RECOVER STEP OF MANAGER
       ! The resume_date is valid for all ICON domains, however,
       ! nstep needs to be recalculated based on the timestep of the
       ! in dividual domains
       CALL manager_state (BM_time, resume_date, nstep, lfit, status)
       IF (status /= 0) CALL ERROR_BI(' ', substr)
#endif
       IF (lfirst_cycle) THEN
#ifndef ICON
          CALL manager_init(BM_time, nstep+1,status)
          CALL timer_message (status, substr)
          CALL info_bi(' Set manager to resumed time step + 1.', 'TIMER')
#else
          CALL manager_init(BM_time, nstep,status)
          CALL timer_message (status, substr)
          CALL info_bi(' Set manager to resumed time step ', 'TIMER')
#endif
       END IF
    ELSE
       CALL manager_state(BM_time,resume_date,nstep,status)
       CALL timer_message (status, substr)
    END IF

  END SUBROUTINE messy_timer_init_manager
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE main_timer_initialize

#ifdef MESSYTIMERTEST
    USE messy_main_mpi_bi,    ONLY: p_parallel_io
#endif

    IMPLICIT NONE
    INTRINSIC :: FLOOR

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='main_timer_initialize'
    INTEGER                     :: status
    INTEGER                     :: ic

    CALL start_message_bi(modstr, 'INITIALIZE TIMER',' ')

#ifdef ICON
    CALL  main_timer_set_clock(1)
    CALL timer_ICON_set_timestep(1, delta_time)

    DO ic = 1,nclock
       CALL main_timer_set_domain(ic)
       CALL messy_timer_init_manager(delta_time, INIT_STEP)
    END DO
#endif

    DO ic = 1, nclock
       CALL  main_timer_set_domain(ic)

       CALL init_BM_time

       ! SET START TIME
       CALL timer_get_date(status, 'start', YEAR_START,MONTH_START,DAY_START, &
            HOUR_START,MINUTE_START,SECOND_START,MILLISECOND_START)

       ! set current time first to start_date, will be overwritten later
       ! qqq needed here for NEW TIMER??

       IF (lstart) &
       CALL timer_get_date(status, 'start', YEAR,MONTH,DAY,HOUR,MINUTE,SECOND &
            ,MILLISECOND)

       IF (lresume) &
       CALL timer_get_date(status, 'resume', YEAR,MONTH,DAY,HOUR,MINUTE,SECOND &
            ,MILLISECOND)

       CALL timer_get_date(status, 'next', YEAR_NEXT,MONTH_NEXT,DAY_NEXT &
            ,HOUR_NEXT,MINUTE_NEXT,SECOND_NEXT,MILLISECOND_NEXT)

       ! SET CURRENT TIME STEP
       current_time_step = timer_get_time_step()

       DAYOFYEAR=FLOOR(YearDay(current_date))

       CALL timer_message(status, substr)
    END DO

    CALL  main_timer_set_domain(1)

    CALL end_message_bi(modstr, 'INITIALIZE TIMER',' ')
    ! --------------------------------------------------------

    ! --------------------------------------------------------
    CALL start_message_bi(modstr, 'INITIALIZE RERUN',' ')

    CALL timer_event_init(  RERUN_EVENT,  IO_RERUN_EV &
         , 'rerun_event', EV_TLEV_NEXT)

    CALL end_message_bi(modstr, 'INITIALIZE RERUN',' ')

  END SUBROUTINE main_timer_initialize
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE main_timer_time(flag, nclck)

#ifdef COSMO
    USE data_runcontrol, ONLY: ntstep, nstop, nfinalstop
#endif

    IMPLICIT NONE

    !I/O
    INTEGER, INTENT(IN)           :: flag
    INTEGER, INTENT(IN), OPTIONAL :: nclck

    INTEGER :: nc

    IF (PRESENT(nclck)) THEN
       nc = nclck
    ELSE
       nc = 1
    END IF

    SELECT CASE(flag)
    CASE(1)
       CALL main_timer_set_domain(nc)
       CALL main_timer_set_time
#ifdef MESSYTIMERTEST
       IF (p_parallel_io) THEN
          !write (*,*) 'FLAG=1: ', nclck, lbreak, lstop, l_rerun
          write (*,*) 'main_timer_time01:' &
            , YEAR,' ',MONTH,' ', DAY,' ', HOUR,' ' &
            , MINUTE,' ', SECOND,' ', MILLISECOND
       END IF
#endif

    CASE(2)
       CALL main_timer_set_domain(nc)
       CALL main_timer_reset_time
#ifdef MESSYTIMERTEST
       IF (p_parallel_io) THEN
          !write (*,*) 'FLAG=2: ', nclck, lbreak, lstop, l_rerun
          write (*,*) 'main_timer_time02:' &
            , YEAR,' ',MONTH,' ', DAY,' ', HOUR,' ' &
            , MINUTE,' ', SECOND,' ', MILLISECOND
       END IF
#endif
#ifdef COSMO
       IF (ntstep+1 == nstop .AND. nstop == nfinalstop) THEN
          L_TRIGGER_RESTART = .TRUE.
       END IF
#endif
     END SELECT

  END SUBROUTINE main_timer_time

  !----------------------------------------------------------------------------

  SUBROUTINE main_timer_free_memory

     IMPLICIT NONE

    INTEGER :: ic

    NULLIFY(BM_time)

    DO ic = 1,nclock
       DEALLOCATE(BMclocks(ic)%BM_time)
    END DO
    DEALLOCATE(BMclocks)

#ifdef ICON
    DEALLOCATE(lfirstreset)
#endif
    CALL timer_dealloc_clock

  END SUBROUTINE main_timer_free_memory

  !----------------------------------------------------------------------------

  SUBROUTINE  main_timer_set_manager(nc)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nc

     BM_time => BMclocks(nc)%BM_time

   END SUBROUTINE main_timer_set_manager

  !----------------------------------------------------------------------------

  SUBROUTINE  main_timer_set_domain(np)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: np

    CALL main_timer_set_manager(np)
    CALL main_timer_set_clock(np)

   END SUBROUTINE main_timer_set_domain

  !----------------------------------------------------------------------------
  SUBROUTINE main_timer_set_time

    ! MAKE TIME STEPPING FOR CLOCK AND UPDATE EVENTS

    ! evaluate events and date/time control elements at the
    ! beginning of the time loop
#ifdef ICON
    USE messy_main_tools, ONLY: str
#endif

    IMPLICIT NONE

    INTRINSIC :: FLOOR

    CHARACTER(LEN=*), PARAMETER :: substr='main_timer_set_time'
    INTEGER                     :: status
    LOGICAL                     :: lsmaller, lequal
    CHARACTER(LEN=STRLEN_SHORT) :: strg = ' '

    ! evaluate events
    l_rerun = .FALSE.
    IF (current_clock == 1) THEN
       IF (.NOT. lforcedtime) THEN
#ifdef MESSYTIMERTEST
          CALL event_print(RERUN_EVENT,status, .FALSE., p_parallel_io)
#endif
#if !defined(ICON)
          l_rerun  = event_state(RERUN_EVENT, next_date) .OR. L_TRIGGER_RESTART
#else
          l_rerun  = event_state(RERUN_EVENT, current_date) .OR. L_TRIGGER_RESTART
#endif
#ifdef MESSYTIMERTEST
          IF (p_parallel_io) write (*,*) 'RERUN event start', l_rerun
          CALL write_date(next_date,'RERUN EVAL DATE ')
#endif
       END IF
       ! evaluate the stop of the model
       CALL if_less(stop_date,next_date,lsmaller,status)
       CALL timer_message(status,substr)
       CALL if_equal(stop_date,next_date,lequal,status)
       CALL timer_message(status,substr)
       lstop = (lsmaller .OR. lequal)

       IF (lstop) THEN
          l_rerun = .TRUE.
          L_TRIGGER_RESTART  = .TRUE.
       ENDIF

       IF (l_rerun) THEN
          IF (current_clock == 1) THEN
             IF (no_cycles > 0) THEN
                lcycbreak  = (ccount >= no_cycles)   ! break rerun cycles
                ccount = ccount + 1                  ! count rerun intervals
             ELSE
                lcycbreak  = .FALSE.                 ! break rerun cycles
             END IF
             lbreak = lcycbreak .OR. L_TRIGGER_RESTART
          END IF
       END IF
#ifdef MESSYTIMERTEST
       IF (p_parallel_io) write (*,*) 'restart flags BSRT: ',HOUR,MINUTE,SECOND&
            , lbreak, lstop, l_rerun, L_TRIGGER_RESTART
#endif

       ! print settings
       IF (lstop) THEN
          CALL write_date(next_date &
               ,'Stop model, last prognostic date/time is: ')
       ELSE IF (lbreak) THEN
          CALL write_date(next_date &
               ,'Interrupt model, last prognostic date/time is: ')
       END IF
    ELSE
       IF (timer_get_l_rerun(1)) THEN
          ! set l_rerun=T in last iteration only
          CALL if_equal(next_date, timer_get_next_date(1), lequal)
          IF (lequal) THEN
             l_rerun = .TRUE.
          END IF
       END IF
       ! SET LSTOP
       IF (lstop) THEN
          ! set lstop=T in last iteration only
          CALL if_equal(next_date, timer_get_next_date(1), lequal)
          l_rerun = .TRUE.
          L_TRIGGER_RESTART  = .TRUE.
#ifdef COSMO
       CALL timer_COSMO_resetfinalstop
#endif
       ENDIF
    END IF
#ifdef MESSYTIMERTEST
    IF (p_parallel_io) write (*,*) substr, 'CCC FLAGS 1: ',current_clock &
         ,  lbreak, lstop, l_rerun, L_TRIGGER_RESTART
#endif
    ! calculate integration interval
    CALL manager_state(BM_time,delta_time,status)
    CALL timer_message(status,substr)

    ! SET CURRENT DATE / TIME
    CALL timer_get_date(status, 'current'     &
         ,YEAR,MONTH,DAY,HOUR,MINUTE,SECOND,MILLISECOND)
    CALL timer_message(status, substr)
    ! SET NEXT DATE / TIME
    CALL timer_get_date(status, 'next'        &
         ,YEAR_NEXT,MONTH_NEXT,DAY_NEXT,HOUR_NEXT,MINUTE_NEXT,SECOND_NEXT &
         ,MILLISECOND_NEXT)
    CALL timer_message(status, substr)

    ! SET CURRENT TIME STEP
    current_time_step = timer_get_time_step()

    ! GET DAY OF YEAR (e.g. 1 Feb = 32)
    DAYOFYEAR=FLOOR(YearDay(current_date))

    IF (lprintnewdate) THEN
       IF (p_parallel_io) THEN
#ifdef ICON
          strg='CLOCK '//TRIM(str(current_clock))
#endif
          CALL write_date(current_date, TRIM(strg)//' Current date is: ')
       END IF
    END IF

#ifdef MESSYTIMERTEST
     IF (p_parallel_io)  THEN
        write (*,*) ' timer_set_time CURRENT ' &
         , YEAR,MONTH,DAY,HOUR,MINUTE,SECOND
        write (*,*) ' timer_set_time NEXT ' &
         , YEAR_NEXT,MONTH_NEXT,DAY_NEXT,HOUR_NEXT,MINUTE_NEXT,SECOND_NEXT
     END IF
#endif

  END SUBROUTINE main_timer_set_time
  !----------------------------------------------------------------------------

  SUBROUTINE main_timer_reset_time

    ! evaluate actions at the end of the time step loop

    INTEGER :: istep, status
    CHARACTER(LEN=*), PARAMETER :: substr='main_timer_reset_time'

    ! reset time manager and dates
    ! 1. re-set previous date
    CALL manager_state(BM_time,previous_date,status)
    CALL timer_message(status, substr)

#ifdef MESSYTIMERTEST
    IF (p_parallel_io) CALL write_date(previous_date &
         ,'reset_time previous date: ')
#endif

    ! increment manager position by 1
    CALL manager_init (BM_time,1, status)
    CALL timer_message(status, substr)

    ! redefine time window
    ! get new manager position istep
    CALL manager_state(BM_time, istep,status)
    CALL timer_message(status, substr)
    ! re-define current_date
    CALL manager_state(BM_time, current_date,status)
    CALL timer_message(status, substr)
#ifdef MESSYTIMERTEST
   IF (p_parallel_io) CALL write_date(current_date &
         ,'reset_time current date: ')
#endif
    ! re-define next_date
    CALL manager_state(BM_time, next_date,istep+1,status)
    CALL timer_message(status, substr)
#ifdef MESSYTIMERTEST
    IF (p_parallel_io) THEN
       CALL write_date(next_date,'reset_time next date: ')

       CALL manager_print(BM_time, status)
       CALL timer_message(status,substr)
    ENDIF
#endif

    ! reset switches
#ifdef ICON
    IF (.NOT. lfirstreset(current_clock)) THEN
#endif
       lstart       = .FALSE.
       lresume      = .FALSE.
       lfirst_cycle = .FALSE.
#ifdef ICON
    ELSE
       lfirstreset(current_clock) = .FALSE.
    END IF
#endif
    !IF (p_parallel_io) write (*,*)'RESET TIME ', istep, INIT_STEP, lstart

  END SUBROUTINE main_timer_reset_time

  ! --------------------------------------------------------------------------

  FUNCTION timer_get_time_step() RESULT (istep)  !****************************

    ! provide time step from echam time manager

    IMPLICIT NONE

    INTEGER :: istep
    INTEGER :: status
    CHARACTER(LEN=*), PARAMETER :: substr='main_timer: get_time_step'

    CALL manager_state(BM_time,istep,status)
    CALL timer_message(status,substr)

  END FUNCTION timer_get_time_step

  ! ---------------------------------------------------------------------------

#ifdef ICON
  ! ---------------------------------------------------------------------------

  SUBROUTINE timer_ICON_get_ndom(nc)

    USE mo_grid_config,    ONLY: dynamics_grid_filename
    USE mo_impl_constants, ONLY: max_dom

    IMPLICIT NONE

    INTEGER, INTENT(INOUT) :: nc
    ! LOCAL
    LOGICAL  :: file_exists

    nc=1
    DO WHILE (dynamics_grid_filename(nc) /= "")
      IF (p_parallel_io) THEN
        INQUIRE (FILE=dynamics_grid_filename(nc), EXIST=file_exists)
        IF (.NOT. file_exists)   THEN
          CALL error_bi( &
               TRIM(dynamics_grid_filename(nc))//" file does not exist" &
               ,"timer_ICON_get_ndom")
        ENDIF
      ENDIF
      nc=nc+1
      IF (nc > max_dom) &
           CALL error_bi( " too many input grids", "timer_ICON_get_ndom")
    END DO
    nc = nc-1

  END SUBROUTINE timer_ICON_get_ndom

  ! ---------------------------------------------------------------------------

  SUBROUTINE timer_initfromICONdates

    USE mo_time_config,        ONLY: time_config
    USE mo_time_config,        ONLY: dt_restart
    USE mo_io_config,          ONLY: dt_checkpoint

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr="timer_initfromICONdates"
    INTEGER :: YEAR, MONTH, DAY, HOUR, MINUTE, SECOND, MS
    INTEGER :: status

    ! set start and end date to values from base model
    ! (ignores timer namelist settings)
    ! in mtime year is int64, for now convert to int32
    ! as underlying code in timer/event only supports int32 for now

    ! start date of thw whole experiment
    YEAR    = INT(time_config%tc_exp_startdate%date%year)
    MONTH   = time_config%tc_exp_startdate%date%month
    DAY     = time_config%tc_exp_startdate%date%day
    HOUR    = time_config%tc_exp_startdate%time%hour
    MINUTE  = time_config%tc_exp_startdate%time%minute
    SECOND  = time_config%tc_exp_startdate%time%second
    MS      = time_config%tc_exp_startdate%time%ms

    CALL timer_set_date(status,'start', YEAR, MONTH, DAY &
         , HOUR, MINUTE, SECOND, MS)

    !start date of the current restart chain
    YEAR    = INT(time_config%tc_current_date%date%year)
    MONTH   = time_config%tc_current_date%date%month
    DAY     = time_config%tc_current_date%date%day
    HOUR    = time_config%tc_current_date%time%hour
    MINUTE  = time_config%tc_current_date%time%minute
    SECOND  = time_config%tc_current_date%time%second
    MS      = time_config%tc_current_date%time%ms

    CALL timer_set_date(status,'resume',YEAR, MONTH, DAY &
         , HOUR, MINUTE, SECOND, MS)

    ! here it is the stop date of the whole experiment,
    ! we need to handle restart seperate!
    YEAR    = INT(time_config%tc_exp_stopdate%date%year)
    MONTH   = time_config%tc_exp_stopdate%date%month
    DAY     = time_config%tc_exp_stopdate%date%day
    HOUR    = time_config%tc_exp_stopdate%time%hour
    MINUTE  = time_config%tc_exp_stopdate%time%minute
    SECOND  = time_config%tc_exp_stopdate%time%second
    MS      = time_config%tc_exp_stopdate%time%ms

    CALL timer_set_date(status,'stop', YEAR, MONTH, DAY &
         , HOUR, MINUTE, SECOND, MS)

    ! SET RERUN event
    CALL timer_set_rerun_event(status, &
         NINT(dt_checkpoint),TIME_INC_SECONDS,TRIG_EXACT,0._dp)

    NO_CYCLES   = NINT(dt_restart/dt_checkpoint)

#ifdef MESSYTIMERTEST
    IF (p_parallel_io) WRITE (*,*) "timer_initfromICONdates NO CYCLES" &
         , NO_CYCLES &
         , dt_checkpoint,dt_restart
#endif

  END SUBROUTINE timer_initfromICONdates

  ! ---------------------------------------------------------------------------

  SUBROUTINE timer_overwriteICONdates

    use, intrinsic :: iso_c_binding, only: c_int, c_int32_t, c_int64_t

    USE mo_master_config,  ONLY: experimentStartDate, experimentStopDate &
         , experimentReferenceDate                 &
         , checkpointTimeIntval, restartTimeIntval &
         , calendar_str
    USE mo_time_config,    ONLY: restart_ini_datetime_string                &
         , restart_end_datetime_string                &
         , ini_datetime_string, end_datetime_string   &
         , set_is_relative_time, dt_restart, icalendar
    USE mo_io_config,      ONLY: dt_checkpoint
    USE mtime_timedelta

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = "timer_overwriteICONdates"
    INTEGER                     :: YYYY, MO, DD, HH, MM, SS, MS
    INTEGER                     :: IYYYY, IMO, IDD, IHH, IMM, ISS, IMS
    INTEGER                     :: status
    integer(c_int64_t)          :: isoCnt

    ! SET CALENDAR
    calendar_str = "proleptic gregorian"  ! qqq we use Julian in TIMER ???

    ! SET EXPERIMENT START DATE
    CALL timer_get_date(status,'start',IYYYY,IMO,IDD,IHH,IMM,ISS,IMS)

    experimentStartDate = ' '

    ! SET ONLY REFERENCE DATE DUE TO COMPLICATIONS WITH DATE CHECKING
    WRITE(experimentReferenceDate,&
         '(i4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1)') &
         IYYYY,'-',IMO,'-',IDD,'T',IHH,':',IMM,':',ISS,'Z'

    ! SET EXPERIMENT STOP DATE
    CALL timer_get_date(status,'stop',YYYY,MO,DD,HH,MM,SS,MS)
    WRITE(experimentStopDate &
         , '(i4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1)') &
         YYYY,'-',MO,'-',DD,'T',HH,':',MM,':',SS,'Z'

    ! SET (RE-)START DATE
    IF (.NOT. LRESUME) THEN
       ini_datetime_string = ' '
    ELSE
       CALL timer_get_date(status,'resume',IYYYY,IMO,IDD,IHH,IMM,ISS,IMS)
       WRITE(ini_datetime_string,&
            '(i4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1)') &
            IYYYY,'-',IMO,'-',IDD,'T',IHH,':',IMM,':',ISS,'Z'
       restart_ini_datetime_string = ini_datetime_string
    END IF

    ! SET STOP  DATE
    end_datetime_string = ""

    ! SET CHECKPOINT AND RESTART DATES
    SELECT CASE (TRIM(IO_RERUN_EV%unit))
    CASE (TIME_INC_SECONDS)
       isoCnt = IO_RERUN_EV%counter
       CALL getPTstringFromSeconds(isoCnt,checkpointTimeIntval, status)
       IF (status /=0) CALL ERROR_BI('error in getPTstring', substr)
       isoCnt= no_cycles * isoCnt
       CALL getPTstringFromSeconds(isoCnt,restartTimeIntval, status)
       IF (status /=0) CALL ERROR_BI('error in getPTstring', substr)
    CASE (TIME_INC_MINUTES)
       isoCnt = IO_RERUN_EV%counter
       CALL getPTstringFromMinutes(isoCnt,checkpointTimeIntval, status)
       IF (status /=0) CALL ERROR_BI('error in getPTstring', substr)
       isoCnt = no_cycles * isoCnt
       CALL getPTstringFromMinutes(isoCnt,restartTimeIntval, status)
       IF (status /=0) CALL ERROR_BI('error in getPTstring', substr)
    CASE (TIME_INC_HOURS)
       isoCnt = IO_RERUN_EV%counter
       CALL getPTstringFromHours(isoCnt,checkpointTimeIntval, status)
       IF (status /=0) CALL ERROR_BI('error in getPTstring', substr)
       isoCnt = no_cycles * isoCnt
       CALL getPTstringFromHours(isoCnt,restartTimeIntval, status)
       IF (status /=0) CALL ERROR_BI('error in getPTstring', substr)
    CASE (TIME_INC_DAYS)
       ! Function getPTstring only available up to hours
       isoCnt = 24 * IO_RERUN_EV%counter
       CALL getPTstringFromHours(isoCnt,checkpointTimeIntval, status)
       IF (status /=0) CALL ERROR_BI('error in getPTstring', substr)
       isoCnt = no_cycles * isoCnt
       CALL getPTstringFromHours(isoCnt,restartTimeIntval, status)
       IF (status /=0) CALL ERROR_BI('error in getPTstring', substr)
    CASE (TIME_INC_MONTHS)
       isoCnt = IO_RERUN_EV%counter * 24 * 30
       CALL getPTstringFromHours(isoCnt,checkpointTimeIntval, status)
       IF (status /=0) CALL ERROR_BI('error in getPTstring', substr)
       isoCnt = no_cycles * isoCnt
       CALL getPTstringFromHours(isoCnt,restartTimeIntval, status)
       IF (status /=0) CALL ERROR_BI('error in getPTstring', substr)
       CALL warning_bi &
            ('EXACT REPRESENTATION OF CHECKPOINT INTERVAL NOT POSSIBLE IN ICON'&
            , substr)
       CALL warning_bi &
            ('... BUT WILL BE OVERWRITTEN ANYWAY ... '&
            , substr)
    CASE (TIME_INC_YEARS)
       isoCnt = IO_RERUN_EV%counter * 24 * 365
       CALL getPTstringFromHours(isoCnt,checkpointTimeIntval, status)
       IF (status /=0) CALL ERROR_BI('error in getPTstring', substr)
       isoCnt = no_cycles * isoCnt
       CALL getPTstringFromHours(isoCnt,restartTimeIntval, status)
       IF (status /=0) CALL ERROR_BI('error in getPTstring', substr)
       CALL warning_bi &
            ('EXACT REPRESENTATION OF CHECKPOINT INTERVAL NOT POSSIBLE IN ICON'&
            , substr)
       CALL warning_bi &
            ('... BUT WILL BE OVERWRITTEN ANYWAY ... '&
            , substr)
    CASE DEFAULT
       CALL error_bi('UNKOWN RESTART UNIT', substr)
    END SELECT

    ! FORCE ICON TIMING TO COUNT EVENTS FROM EXPERIMENT START ON
    CALL set_is_relative_time(.FALSE.)
    icalendar     = -1
    dt_restart    = 0._dp
    dt_checkpoint = 0._dp

#ifdef MESSYTIMERTEST
    IF (p_parallel_io) THEN
       write (*,*) substr, 'experimentReferenceDate:' &
            , TRIM(experimentReferenceDate)
       write (*,*) substr, 'experimentStartDate:    ',TRIM(experimentStartDate)
       write (*,*) substr, 'experimentStopDate:     ',TRIM(experimentStopDate)
       write (*,*) substr, 'ini_datetime_string     ',TRIM(ini_datetime_string)
       write (*,*) substr, 'CheckPointTimeIntval:   ',TRIM(checkpointTimeIntval)
       write (*,*) substr, 'RestartTimeIntval:      ',TRIM(restartTimeIntval)
    ENDIF
#endif

  END SUBROUTINE timer_overwriteICONdates

  ! ---------------------------------------------------------------------------

  RECURSIVE SUBROUTINE timer_ICON_set_timestep(jg, dt_loc)

   ! unfortunately we have to use the original ICON module here,
   ! otherwise (USEing BMLUSE or DATA) a ECHAM5 develops a circular
   ! dependency of mo_time_control and timer_bi

   USE mo_model_domain,          ONLY: p_patch

    IMPLICIT NONE

    INTEGER,  INTENT(IN) :: jg
    REAL(dp), INTENT(IN) :: dt_loc

    INTEGER :: jn, jgc
    REAL(dp):: dt_sub

    CALL main_timer_set_clock(jg)
#ifdef MESSYTIMERTEST
    IF (p_parallel_io) &
         write (*,*) 'ICON_set_timestep ', dt_loc,' for patch ',jg
#endif
    delta_time    = dt_loc
    time_step_len = dt_loc

    dt_sub      = dt_loc/2._dp

    DO jn = 1, p_patch(jg)%n_childdom
       jgc = p_patch(jg)%child_id(jn)
       CALL timer_ICON_set_timestep(jgc,dt_sub)
    END DO

  END SUBROUTINE timer_ICON_set_timestep

  ! ---------------------------------------------------------------------------

#endif

#ifdef COSMO

  ! --------------------------------------------------------------------------

  SUBROUTINE messy_timer_COSMO_reinit_time

    ! COSMO
    USE data_runcontrol,    ONLY: hstart, nstart, nstop, hstop, nfinalstop &
         , hlastmxu,   hnextmxu,   hincmxu   &
         , hlastmxt,   hnextmxt,   hincmxt   &
         , nlastmxu,   nnextmxu,   nlastmxt  &
         , nnextmxt                          &
         , yakdat1,    yakdat2, itype_calendar

    USE data_modelconfig,   ONLY: dt
    USE data_io,            ONLY: ydate_ini
    USE utilities,          ONLY: get_utc_date
    ! MESSy
    USE messy_main_constants_mem, ONLY: sp!, iouerr

    IMPLICIT NONE
    INTRINSIC :: NINT, REAL

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'reinit_COSMO_time'
    REAL(dp) :: timediff ! time difference in second (resume-start date)
    INTEGER  :: ryr, rmo, rdy, rhr, rmn, rse  ! resume date components
    INTEGER  :: dummy1
    REAL(dp) :: dummy2
    INTEGER  :: status

    CALL timer_get_date(status, 'resume', ryr, rmo, rdy, rhr, rmn, rse)
    CALL timer_message(status, substr)
    CALL time_span_d( timediff , YEAR_START, MONTH_START,  DAY_START    &
         , HOUR_START, MINUTE_START, SECOND_START &
         , ryr, rmo,   rdy, rhr,     rmn, rse     )
    hstart = REAL(REAL(timediff * 24.0_dp,sp),dp) ! resume hour
    nstart = NINT(timediff*86400.0_dp /dt)
    !write(iouerr,*) ' HSTART ',YEAR_START, MONTH_START,  DAY_START    &
    !  , HOUR_START, MINUTE_START, SECOND_START, timediff

    CALL timer_get_date(status, 'stop', ryr, rmo, rdy, rhr, rmn, rse)
    CALL timer_message(status, substr)
    CALL time_span_d( timediff , YEAR_START, MONTH_START,  DAY_START    &
         , HOUR_START, MINUTE_START, SECOND_START &
         , ryr, rmo,   rdy, rhr,     rmn, rse     )
    hstop  =  timediff * 24.0_dp
    nfinalstop = NINT(timediff*86400.0_dp /dt) - 1

#ifdef CCLMORI
    IF (NO_STEPS > 0) THEN
       nstop = NO_STEPS
    ELSE IF (NO_DAYS > 0) THEN
       nstop = NO_DAYS * 86400.0_dp / dt
    END IF
#else
    nstop = nfinalstop
#endif

    ! RE-DO some CHECKS WHICH ARE PERFORMED ORIGIALLY IN ORGANIZE_SETUP
    ! endless loop for finding the last hour (for restart runs)

    hlastmxu = 0.0_dp
    endless_u: DO
       IF ( (hlastmxu <= hstart) .AND. (hstart <= hlastmxu + hincmxu) ) THEN
          EXIT endless_u
       ENDIF
       hlastmxu = hlastmxu + hincmxu
    ENDDO endless_u
    hnextmxu = hlastmxu + hincmxu
    nlastmxu = NINT (hlastmxu * 3600.0_dp / dt)
    nnextmxu = NINT (hnextmxu * 3600.0_dp / dt)

    hlastmxt = 0.0_dp
    endless_t: DO
       IF ( (hlastmxt <= hstart) .AND. (hstart <= hlastmxt + hincmxt) ) THEN
          EXIT endless_t
       ENDIF
       hlastmxt = hlastmxt + hincmxt
    ENDDO endless_t
    hnextmxt = hlastmxt + hincmxt
    nlastmxt = NINT (hlastmxt * 3600.0_dp / dt)
    nnextmxt = NINT (hnextmxt * 3600.0_dp / dt)

    ! compute the actual date
    CALL get_utc_date(nstart, ydate_ini, dt, itype_calendar, yakdat1 &
         , yakdat2, dummy1, dummy2)

  END SUBROUTINE messy_timer_COSMO_reinit_time

  ! --------------------------------------------------------------------------

  SUBROUTINE timer_COSMO_resetfinalstop

    USE data_runcontrol,   ONLY: nfinalstop

    IMPLICIT NONE

    nfinalstop=current_time_step+1


  END SUBROUTINE timer_COSMO_resetfinalstop

  ! --------------------------------------------------------------------------

#endif

  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  ! EVENT_MANAGER
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------

  SUBROUTINE timer_event_init (event, io_ev, ev_name, eval_date &
       , clock_id)

    USE messy_main_tools,         ONLY: int2str

    ! the initial state of an event is calculated
    ! it is dependent on the initial date of a run
    ! at rerun the event triggers will be recalculated from beginning

    IMPLICIT NONE
    INTRINSIC :: NINT, TRIM

    TYPE(time_event),    INTENT(inout) :: event     ! will be initialized
    TYPE(io_time_event), INTENT(inout) :: io_ev     ! pass external settings
    CHARACTER(len=*)                   :: ev_name   ! submit a name
    CHARACTER(len=*)                   :: eval_date ! the evaluation date
    INTEGER,  OPTIONAL,  INTENT(IN)    :: clock_id

    TYPE(time_days)             :: my_date
    REAL(dp)                    :: zdtime, zztime
    CHARACTER(len=STRLEN_SHORT) :: newunit
    INTEGER                     :: newcounter, idtime
    CHARACTER(LEN=*),PARAMETER  :: substr='timer_event_init'
    CHARACTER(len=STRLEN_ULONG) :: message
    INTEGER                     :: status
    LOGICAL                     :: lsm
    CHARACTER(LEN=3)            :: scid
    TYPE(time_manager), POINTER :: locBMtime => NULL() !base model time axis
    INTEGER                     :: locINIT_STEP
    INTEGER                     :: clockmax
#ifdef ICON
    LOGICAL                     :: leq
#endif

    IF (PRESENT(clock_id)) THEN
#ifdef ICON
       clockmax = MAX(1,clock_id)
#else
       clockmax= clock_id
#endif
       CALL INT2STR(scid,clock_id,'_')
       locBMtime => BMclocks(clockmax)%BM_time
       CALL timer_get_init_step(locINIT_STEP, clock_id=clockmax)
    ELSE
       scid = '   '
       locBMtime    => BM_time
       locINIT_STEP = INIT_STEP
    END IF

    IF (p_parallel_io) write (*,*) substr, ' EVENT NAME: ', ev_name//TRIM(scid)

    CALL manager_state(locBMtime,zztime,status)
    CALL timer_message(status,substr)
    idtime = NINT(zztime)
    IF (ABS(REAL(idtime,dp) - zztime) >  TINY(1._dp )) THEN
       zdtime = zztime
    ELSE
       zdtime = zztime - 0.001_dp
    END IF

    ! convert from steps into other units
    IF (TRIM(io_ev%unit) == TIME_INC_STEPS &
         .OR. TRIM(io_ev%unit) == TIME_INC_ALWAYS) THEN

       IF (TRIM(io_ev%unit) == TIME_INC_ALWAYS) THEN
          ! for debugging and testing only
             io_ev%counter = 1
       END IF

       CALL convert_steps2unit(io_ev%counter,zztime,newunit,newcounter &
            , status, message)
       CALL timer_message(status,substr, message)

       SELECT CASE(TRIM(newunit))
       CASE(TIME_INC_MILSECS)
          CALL timer_message(0, ' ' &
               ,'Convert event interval from STEPS into MILLISECONDS')
       CASE(TIME_INC_SECONDS)
          CALL timer_message(0, ' ' &
               ,'Convert event interval from STEPS into SECONDS')
       CASE(TIME_INC_MINUTES)
          CALL timer_message(0,' '&
               ,'Convert event interval from STEPS into MINUTES')
       END SELECT
    ELSE
       newunit    = io_ev%unit
       newcounter = io_ev%counter
    END IF

    ! initialize event with external settings
    CALL event_init(event, ev_name &
         , newcounter, newunit, io_ev%adjustment, zdtime &
         , io_ev%offset, status, message, tag=clock_id)
    CALL timer_message(status,substr,message)

#ifdef MESSYTIMERTEST
    CALL event_print(event,status, .TRUE., p_parallel_io)
#endif

    ! check time stepping against event intervals
    SELECT CASE(io_ev%adjustment)
    CASE(TRIG_FIRST,TRIG_LAST)

       SELECT CASE(newunit)
       CASE(TIME_INC_MILSECS,TIME_INC_SECONDS, TIME_INC_MINUTES &
            , TIME_INC_HOURS, TIME_INC_DAYS)

          IF (event_eval(event,zztime) < 0) THEN
#ifdef MESSYTIMERTEST
             IF (p_parallel_io) THEN
                write (*,*) substr, ' io_ev%unit', newunit
                write (*,*) substr, ' ev%unit', event_unit(event)
             ENDIF
             CALL event_print(event, status, .FALSE., p_parallel_io )
#endif
             CALL timer_message(3400, substr &
                  , 'Event counter mismatch with time stepping.')
          END IF

       CASE(TIME_INC_MONTHS, TIME_INC_YEARS)
          CALL timer_message(0,' ',&
               'No time stepping mismatch check defined for MONTHS and YEARS.')

       CASE default
          CALL timer_message(3400,substr,'Event unit not defined.')

       END SELECT
    END SELECT

    ! define initial date and trigger dates
    CALL manager_state(locBMtime,my_date,locINIT_STEP,status)
    CALL timer_message(status,substr)
#ifdef MESSYTIMERTEST
    IF (p_parallel_io) CALL write_date(my_date,'my_date: timer_event_init 1')
#endif

    CALL event_init(event,my_date,.FALSE.,status,message)
    CALL timer_message(status,substr, message)

    ! find next trigger date starting at initial date
    IF (lresume) THEN

       !================= preliminar
       ! the finding of the next possible trigger can take a lot of
       ! time if the present date is very fare from the start date
       !
       ! revision needed with a new I/O concept:
       !     event elements must be available in a rerun file

       DO
          CALL event_next_date(event,my_date)
#ifdef MESSYTIMERTEST
          IF (p_parallel_io) THEN
             CALL write_date(my_date, &
                  TRIM(event_name(event))//' next_event_date')
             CALL write_date(next_date, &
                  TRIM(event_name(event))//' next_date BM_time')
             CALL write_date(current_date, &
                  TRIM(event_name(event))//' current_date BM_time')
             write (*,*) substr, 'EVAL DATE: ',eval_date
          END IF
#endif

          SELECT CASE(eval_date)
          CASE(EV_TLEV_PRES)          ! check with current date
             CALL if_less(my_date,current_date,lsm,status)
             CALL timer_message(status,substr)
#ifdef ICON
             CALL if_equal(my_date,current_date, leq, status)
             CALL timer_message(status,substr)
             lsm = lsm .OR. leq
#endif
             IF (.NOT.lsm) EXIT

          CASE(EV_TLEV_NEXT)          ! check with next date
             CALL if_less(my_date,next_date,lsm,status)
             CALL timer_message(status,substr)
#ifdef ICON
             CALL if_equal(my_date,next_date, leq, status)
             CALL timer_message(status,substr)
             lsm = lsm .OR. leq
#endif
             IF (.NOT. lsm) EXIT

          CASE(EV_TLEV_PREV)          ! check with previous date
             CALL if_less(my_date,previous_date,lsm,status)
             CALL timer_message(status,substr)
#ifdef ICON
             CALL if_equal(my_date,previous_date, leq, status)
             CALL timer_message(status,substr)
             lsm = lsm .OR. leq
#endif
             IF (.NOT. lsm) EXIT

          END SELECT

          ! next date smaller as current date, rotate all dates
          CALL event_init(event,my_date,.TRUE.,status, message)
          CALL timer_message(status,substr, message)

       END DO
#ifdef MESSYTIMERTEST
       IF (p_parallel_io) CALL write_date(my_date, 'next_date 2')
#endif
    END IF

    !print out event settings
    CALL event_print(event, status, .TRUE., p_parallel_io )

  END SUBROUTINE timer_event_init

  !-----------------------------------------------------------------------------

  SUBROUTINE p_bcast_event (event, p_source, comm)

    ! distribute events to all nodes

    USE messy_main_mpi_bi,  ONLY: p_bcast
#ifdef COSMO
    USE messy_main_constants_mem, ONLY: i4
#endif

    IMPLICIT NONE
#ifdef COSMO
    INTRINSIC :: INT
#endif

    TYPE(io_time_event), INTENT(inout) :: event
    INTEGER,             INTENT(in)    :: p_source
    INTEGER, OPTIONAL,   INTENT(in)    :: comm

    ! LOCAL
    INTEGER :: isender
#ifdef COSMO
    INTEGER :: icomm
#endif

    isender = INT(p_source)
#ifdef COSMO
    IF (PRESENT(comm)) THEN
       icomm   = INT(comm,i4)
       CALL p_bcast(event%counter,    isender, icomm=icomm)
       CALL p_bcast(event%unit,       isender, icomm=icomm)
       CALL p_bcast(event%adjustment, isender, icomm=icomm)
       CALL p_bcast(event%offset,     isender, icomm=icomm)
    ELSE
#endif
#ifdef ICON
    IF (PRESENT(comm)) THEN
       CALL p_bcast(event%counter,    isender, comm)
       CALL p_bcast(event%unit,       isender, comm)
       CALL p_bcast(event%adjustment, isender, comm)
       CALL p_bcast(event%offset,     isender, comm)
    ELSE
#endif
       CALL p_bcast(event%counter,    isender)
       CALL p_bcast(event%unit,       isender)
       CALL p_bcast(event%adjustment, isender)
       CALL p_bcast(event%offset,     isender)
#if defined(COSMO) || defined(ICON)
    ENDIF
#endif
  END SUBROUTINE p_bcast_event

  ! ---------------------------------------------------------------------------

  FUNCTION I_event_steps (event, delta) RESULT (ix)  !************************

    ! get out interval between two event trigger points in steps

    IMPLICIT NONE
    INTRINSIC :: NINT

    TYPE (time_event), INTENT(in) :: event
    REAL(dp),          INTENT(in) :: delta   ! seconds of special unit
    INTEGER                       :: ix
    REAL(dp)                      :: rsecs

    rsecs = R_event_seconds(event)
    ix = NINT(rsecs/delta)             ! convert into steps

  END FUNCTION I_event_steps

  ! ---------------------------------------------------------------------------

  FUNCTION R_event_seconds (event) RESULT (rx) !*******************************

    ! get out interval between two event trigges in seconds

    IMPLICIT NONE
    INTRINSIC :: REAL

    TYPE (time_event), INTENT(in) :: event
    REAL(DP)                       :: rx

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='main_timer event_seconds'
    INTEGER           :: iday, isec, ims
    TYPE (time_days)  :: my_date, previous_date
    INTEGER           :: status
    REAL(DP)          :: secfrac

    ! get difference between last and present trigger date
    CALL event_previous_date(event, previous_date)
    CALL date_get(previous_date, iday, isec, ims, status)
    CALL timer_message(status,substr)

    secfrac = REAL(isec,dp) + REAL(ims,dp) / 1000._dp
    iday = - iday
    isec = - isec
    ims = - ims

    CALL event_current_date(event, my_date)
    CALL add_date (iday, secfrac, my_date, status)
    CALL timer_message(status,substr)

    ! check result
    CALL date_get(my_date,iday,isec,ims,status)
    CALL timer_message(status,substr)

    IF (iday < 0 .OR. isec < 0 .OR. ims < 0) &
         CALL timer_message(3400,substr,'Event interval < 0')

    ! convert into seconds
    rx = REAL(iday,dp) * OneDay + REAL(isec,dp) + REAL(ims,dp) / 1000._dp

    ! at initial time the interval is zero
    ! this (may be) is important for accumulation

  END FUNCTION R_event_seconds

  ! ---------------------------------------------------------------------------

  FUNCTION R_event_seconds_next (event, lnext) RESULT (rx)  !******************

    ! return distance between present and next trigger

    IMPLICIT NONE
    INTRINSIC :: REAL

    TYPE (time_event), INTENT(in) :: event
    LOGICAL,           INTENT(in) :: lnext
    REAL(DP)                       :: rx

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='main_timer event_seconds_next'
    INTEGER           :: iday, isec, ims
    TYPE (time_days)  :: my_date, current_date
    INTEGER           :: status
    REAL(DP)          :: secfrac

    rx = 0.0_dp
    IF (lnext) THEN
       CALL event_current_date(event, current_date)
       CALL date_get(current_date,iday,isec,ims,status)
       CALL timer_message(status,substr)

       CALL event_next_date(event, my_date)
       secfrac = REAL(isec,dp) + REAL(ims,dp) / 1000._dp
       CALL add_date (-iday, -secfrac, my_date,status)
       CALL timer_message(status,substr)

       CALL date_get(my_date,iday,isec,ims,status)
       CALL timer_message(status,substr)

       IF (iday < 0 .OR. isec < 0 .OR. ims < 0) &
            CALL timer_message(3400, substr, 'Event interval < 0')

       ! convert into seconds
       rx = REAL(iday,dp) * OneDay + REAL(isec,dp) + REAL(ims,dp) / 1000._dp

    END IF

  END FUNCTION R_event_seconds_next

  ! ---------------------------------------------------------------------------

  LOGICAL FUNCTION L_event_trigger (event, date)
    !
    IMPLICIT NONE
    INTRINSIC :: INT, REAL, TRIM

    TYPE (time_event), INTENT(inout) :: event
    TYPE (time_days),  INTENT(in)    :: date

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='main_timer L_event_trigger'
    REAL(dp)         :: zsec  ! length of event increment in seconds
    INTEGER          :: iday1, iday2
    REAL(dp)         :: rsec1, rsec2
    INTEGER          :: status
    TYPE (time_days) :: adj_start, adj_stop
    LOGICAL          :: lequal1, lequal2, learly1, learly2
    CHARACTER(LEN=STRLEN_ULONG) :: m_text

    L_event_trigger = .FALSE.
    m_text = ' '

    IF (.NOT. event_is_active(event))  RETURN
    ! check if the date will fit into the trigger interval
    zsec = event_delta(event)
    IF (TRIM(event_adjust(event)) == TRIM(TRIG_EXACT)) &
         zsec = event_half_delta(event)
    iday1 = NINT((zsec)/OneDay)
    rsec1 = zsec - REAL(iday1,dp)*OneDay

    IF (TRIM(event_adjust(event)) == TRIM(TRIG_EXACT)) THEN
       zsec  = event_half_delta(event)
       iday2 = INT((zsec-1._dp+0.00001_dp)/OneDay)
       rsec2 = zsec - REAL(iday2,dp)*OneDay
    END IF

    ! detect first or last second in given unit
    CALL event_next_date(event, adj_start)
    CALL event_next_date(event, adj_stop)

    SELECT CASE(TRIM(event_adjust(event)))
    CASE(TRIM(TRIG_FIRST))
       CALL add_date( iday1, rsec1, adj_stop,status)
       CALL timer_message(status,substr)
    CASE(TRIM(TRIG_LAST))
       CALL add_date(-iday1,-rsec1, adj_start,status)
       CALL timer_message(status,substr)
    CASE(TRIM(TRIG_EXACT))
       CALL add_date( iday1, rsec1, adj_stop,status)
       CALL timer_message(status,substr)
       CALL add_date(-iday2,-rsec2, adj_start,status)
       CALL timer_message(status,substr)
    END SELECT

#ifdef MESSYTIMERTEST
    IF (p_parallel_io) THEN
       CALL write_date(date,      ' date1: ' //TRIM(event_name(event)))
       CALL write_date(adj_start, ' date2: ' //TRIM(event_name(event)))
       CALL write_date(adj_stop,  ' date3: ' //TRIM(event_name(event)))
    ENDIF
#endif

    ! check the position of the present date
    CALL if_less(adj_stop, date, learly1,status)
    CALL timer_message(status,substr)

    IF (learly1) THEN
       m_text = 'Warning: event trigger not longer in future <'&
            // TRIM(event_name(event)) // '>'
       CALL timer_message(0,substr,m_text)

    ELSE    ! adj_start <= my_date <= adj_stop
       CALL if_equal(adj_start,date,lequal1,status)
       CALL timer_message(status,substr)
       CALL if_equal(adj_stop,date,lequal2,status)
       CALL timer_message(status,substr)
       CALL if_less(adj_start,date,learly1,status)
       CALL timer_message(status,substr)
       CALL if_less(date,adj_stop,learly2,status)
       CALL timer_message(status,substr)

       L_event_trigger = lequal1 .OR. lequal2 .OR. (learly1 .AND. learly2)

       IF (L_event_trigger)  THEN
          CALL event_init (event, date, .TRUE. , status, m_text)
          CALL timer_message(status,substr,m_text)
       END IF

    ENDIF

  END FUNCTION L_event_trigger

  ! ---------------------------------------------------------------------------

  FUNCTION event_eval (event, delta) RESULT (slen)

    ! check the multiple of delta in event interval

    IMPLICIT NONE
    INTRINSIC :: INT, MOD, PRESENT, TRIM

    TYPE(time_event),   INTENT(in) :: event
    REAL(dp), OPTIONAL, INTENT(in) :: delta

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='main_timer event_eval'
    INTEGER                     :: slen
    INTEGER                     :: fit
    CHARACTER(LEN=STRLEN_ULONG) :: m_text

    slen = -1

    IF (.NOT. event_is_init(event)) THEN
       WRITE (*,*) ('event not initialised, no unit given.')

    ELSE
       SELECT CASE(event_unit(event))
       CASE(TIME_INC_MILSECS);   slen = event_count(event)
       CASE(TIME_INC_SECONDS);   slen = event_count(event)
       CASE(TIME_INC_MINUTES);   slen = event_count(event) * 60
       CASE(TIME_INC_HOURS);     slen = event_count(event) * 60 * 60
       CASE(TIME_INC_DAYS);      slen = event_count(event) * 60 * 60 * 24
       CASE(TIME_INC_MONTHS, TIME_INC_YEARS)
          WRITE (*,*) substr,'Exact value undefined.'
          WRITE (*,*) substr &
               ,'Event trigger interval for months or years may varied.'
       CASE default
          m_text = 'Counter unit unknown ::'//TRIM(event_unit(event))
          CALL timer_message(3400, TRIM(substr), TRIM(m_text))
       END SELECT

       IF (slen > 0 .AND. PRESENT(delta)) THEN
          fit = MOD(INT(slen*1000._dp),INT(delta*1000._dp))
          IF (fit /= 0) THEN
             WRITE(m_text,*) 'Event <',TRIM(event_name(event))&
                  ,   '> interval not a multiple of ',INT(delta)
             CALL timer_message(1,substr,m_text)
             slen = -fit
          END IF
       END IF

    END IF

  END FUNCTION event_eval

  ! ---------------------------------------------------------------------------

  SUBROUTINE event_print (event, status, short, lwrite)

    ! print event contents
    ! CALL for IO PE only

    IMPLICIT NONE
    INTRINSIC :: INT, NINT, PRESENT, TRIM

    TYPE (time_event), INTENT(in) :: event
    INTEGER,           INTENT(OUT):: status
    LOGICAL, OPTIONAL, INTENT(in) :: short
    LOGICAL, OPTIONAL, INTENT(IN) :: lwrite

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='main_timer: event_print'
    TYPE (time_days)            :: previousdate, currentdate, nextdate
    LOGICAL                     :: lequal, llshort
    CHARACTER(LEN=STRLEN_ULONG) :: date_text

    status = 0
    IF (PRESENT(lwrite)) THEN
       IF (.NOT. lwrite) RETURN
    ENDIF

    IF (event_is_init(event)) THEN

       IF (.NOT. PRESENT(SHORT)) THEN
          llshort = .FALSE.
       ELSE
          llshort = short
       ENDIF
       IF (llshort) THEN
          WRITE(*,*) 'Event <',TRIM(event_name(event)),&
               '> : interval ',event_count(event),' '&
               ,TRIM(event_unit(event)),&
               ' : adjustment ',TRIM(event_adjust(event)) &
               ,' : offset[sec] ',event_offset(event)
       ELSE
          IF (.NOT. event_is_active(event)) THEN
             WRITE(*,*) 'Event <',TRIM(event_name(event)),'> ... not active'
          ELSE
             WRITE(*,*) ' '
             WRITE(*,*) 'State of event <',TRIM(event_name(event)) &
                  ,'> ... initialized'

             SELECT CASE(event_unit(event))
             CASE(TIME_INC_SECONDS,TIME_INC_MILSECS)
                WRITE(*,*) ' trigger each ',event_count(event),' seconds'

             CASE(TIME_INC_MINUTES)
                WRITE(*,*) ' trigger each ',event_count(event),' minutes'

             CASE(TIME_INC_HOURS)
                WRITE(*,*) ' trigger each ',event_count(event),' hours'

             CASE(TIME_INC_DAYS)
                WRITE(*,*) ' trigger each ',event_count(event),' days'

             CASE(TIME_INC_MONTHS)
                WRITE(*,*) ' trigger each ',event_count(event),' months'

             CASE(TIME_INC_YEARS)
                WRITE(*,*) ' trigger each ',event_count(event),' years'

             CASE default
                WRITE(*,*) 'Counter unit unknown ::',event_unit(event)
                STATUS = 3400
                RETURN
             END SELECT

             CALL event_previous_date(event, previousdate)
             CALL print_date_components(previousdate, status &
                  , mess=date_text)
             IF (status /= 0) RETURN
             WRITE (*,*) '  last event trigger date is ...    ' // TRIM(date_text)

             CALL event_current_date(event, currentdate)
             CALL print_date_components(currentdate, status, mess=date_text)
             IF (status /= 0) RETURN
             CALL if_equal(currentdate,previousdate,lequal, status)
             IF (status /= 0) RETURN
             IF (lequal) THEN
                WRITE(*,*) '  initial event date is ...    ' // TRIM(date_text)
             ELSE
                WRITE(*,*) '  time between two triggers: ',&
                     R_event_seconds(event),' seconds'
                WRITE(*,*) '  present trigger date is ...  ' // TRIM(date_text)
             END IF

             CALL event_next_date(event, nextdate)
             CALL print_date_components(nextdate,status, mess=date_text)
             IF (status /=0 ) RETURN
             WRITE(*,*)    '  next trigger date is ...         ',&
                  TRIM(date_text),' (offset of ',event_offset(event) &
                  ,' seconds included)'

             SELECT CASE(TRIM(event_adjust(event)))
             CASE(TRIM(TRIG_FIRST))
                WRITE(*,*) '  adjustment at beginning of unit,',&
                     ' interval: ( trigger : trigger + '&
                     , NINT(event_delta(event)),'s )'

             CASE(TRIM(TRIG_LAST))
                WRITE(*,*) '  adjustment at end of unit,',&
                     ' interval: ( trigger - ', NINT(event_delta(event))&
                     ,'s : trigger )'

             CASE(TRIM(TRIG_EXACT))
                WRITE(*,*) '  no adjustment,',&
                     ' interval: ( trigger - ',INT(event_half_delta(event)),&
                     's : trigger + ',INT(event_half_delta(event))-1,'s )'

             END SELECT

          END IF

       END IF

    ELSE
       STATUS = 1
       WRITE (*,*)  substr,'no printout of uninitialized event'

    END IF

  END SUBROUTINE Event_print

  ! ------------------------------------------------------------------------
  ! ------------------------------------------------------------------------
  ! PRIVATE ROUTINES
  ! ------------------------------------------------------------------------
  ! ------------------------------------------------------------------------

  SUBROUTINE init_BM_time

    ! initialize BASE model dates
    ! CALLED FROM messy_timer_initialize

    IMPLICIT NONE
    INTRINSIC :: ABS

    INTEGER           :: istep, status
    INTEGER           :: yr, mo, dy, hr, mi, se, ms
    REAL(dp)          :: zdtold
    LOGICAL           :: lsmaller, lequal
    TYPE(time_days)   :: date
    CHARACTER(LEN=STRLEN_ULONG) :: m_text
    CHARACTER(LEN=*), PARAMETER :: substr='main_timer: init_BM_time'

    CALL manager_state(BM_time,zdtold,status)
    CALL timer_message(status,substr)

    IF (ABS(zdtold-delta_time) > 0.0_dp) THEN

       IF (lresume) THEN
          CALL info_bi('Attention: length of timestep has been changed!', ' ')
          WRITE (m_text,*) ' Old timestep was ', zdtold, ' seconds'
          CALL info_bi(m_text, ' ')
          WRITE (m_text,*) ' New timestep is  ', delta_time, ' seconds'
          CALL info_bi(m_text,' ')
       END IF
       CALL manager_init(BM_time,delta_time,status)
       CALL timer_message(status,substr)
       ! correct to old start_date (was changed by change of delta_time)
       CALL manager_init(BM_time,start_date,status)
       CALL timer_message(status,substr)
       ! Correct the start_date components
       CALL timer_get_date(status,start_date, yr, mo, dy, hr, mi, se, ms)
       CALL timer_message(status,substr)
    END IF

    IF (.NOT. lresume) THEN

       ! set the time manager to the init_step now
       ! reset time step at beginning of experiment
       ! count backward
       ! qqq needed for COSMO ?
       istep = INIT_STEP - timer_get_time_step()
       CALL manager_init(BM_time,istep,status)
       CALL timer_message(status,substr)
    ELSE

       CALL manager_state(BM_time, date,status)
       CALL timer_message(status,substr)

#ifdef MESSYTIMERTEST
       IF (p_parallel_io) THEN
          CALL write_date(date, '  TTT BMtime date:       ')
          CALL write_date(start_date,  ' TTT BMtime start_date: ')
       END IF
#endif
       CALL if_less(date,start_date,lsmaller,status)
       CALL timer_message(status,substr)
       CALL if_equal(date,start_date,lequal,status)
       CALL timer_message(status,substr)
       IF (lsmaller) THEN
          CALL write_date(date,'Current date ...')
          CALL write_date(start_date,'Start date ...')
          CALL timer_message (3400, substr,'Start date in future')

       ELSE IF (lequal) THEN
          lresume = .FALSE.
          CALL timer_message(0,substr &
               ,'Set start date to resume date, force initial run')
       END IF

    END IF

    ! preliminary initializion of  previous date and next date
    istep = timer_get_time_step()
#ifdef MESSYTIMERTEST
    IF (p_parallel_io) write (*,*) 'TTT istep ', istep
#endif
    CALL manager_state(BM_time,previous_date,istep-1, status)
    CALL timer_message(status, substr)
    CALL manager_state(BM_time,next_date,    istep+1, status)
    CALL timer_message(status, substr)

    CALL timer_message(0, ' ','Time step and start date evaluation done.')

    IF (lfirst_cycle)  CALL manager_init(BM_time,.TRUE.,status,m_text)
    CALL timer_message(status, substr,m_text)

    IF (p_parallel_io) THEN
       CALL manager_print(BM_time, status)
       CALL timer_message(status,substr)
    ENDIF

    ! define time stepping
    istep  = timer_get_time_step()
    lstart = (istep == INIT_STEP)
    time_step_len = delta_time

    ! start date of manager can be changed final definition of start date here
    CALL manager_state(BM_time,current_date, status)
    CALL timer_message(status,substr)

    CALL manager_state(BM_time,start_date,INIT_STEP,status)
    CALL timer_message(status,substr)

    ! stop of experiment evaluated only during first rerun cycle
    ! evaluation in the following order (highest priority left)
    ! NO_STEPS or NO_DAYS or DT_STOP or default

    IF (.NOT. lfirst_cycle) THEN
       CALL timer_message(0,' ','No MESSy evaluation of model stop date.')
    END IF

    ! check date order
    CALL if_less(start_date,stop_date, lsmaller,status)
    CALL timer_message(status,substr)
    IF (.NOT. (lsmaller)) THEN
       CALL timer_message(3441, substr)
    ELSE
       CALL if_less(stop_date,current_date, lsmaller,status)
       CALL timer_message(status,substr)
       IF (lsmaller) THEN
          CALL write_date(stop_date,'Stop date: ')
          CALL write_date(current_date,'smaller than current_date: ')
          CALL write_date(resume_date,'resume_date: ')
          CALL timer_message(3442, substr)
       END IF
    END IF

    CALL write_date(stop_date,'Stop experiment at: ')

    ! initialize previous date and next date
    CALL manager_state(BM_time,previous_date,istep-1,status)
    CALL timer_message(status,substr)
    CALL write_date   (previous_date,'Previous date: ')

    CALL manager_state(BM_time,next_date,    istep+1,status)
    CALL timer_message(status,substr)
    CALL write_date   (next_date,    'Next date    : ')

    CALL timer_message (0,' ',' ')

  END SUBROUTINE init_BM_time

  ! ---------------------------------------------------------------------------

#ifdef COSMO

  ! ---------------------------------------------------------------------------

  SUBROUTINE check_rerun_event(event)

    IMPLICIT NONE

    TYPE(time_event), INTENT(INOUT) :: event

    ! LOCAL
    LOGICAL :: le
    CHARACTER(LEN=*), PARAMETER :: substr = 'check_rerun_event'
    INTEGER                     :: status
    CHARACTER(LEN=STRLEN_ULONG) :: message

    status = 0
    IF (event_is_active(event)) THEN

       ! CHECK if restart stop date is reached
       CALL if_less(rerun_stop_date, current_date, le, status)
       CALL timer_message(status, substr)

       ! DEACTIVE EVENT
       IF (le) THEN
          CALL event_init(event, TRIG_NONE, status, message)
          CALL timer_message(status, substr, message)
       ENDIF
    ENDIF

  END SUBROUTINE check_rerun_event

  ! ---------------------------------------------------------------------------
#endif

  ! -------------------------------------------------------------------------

  SUBROUTINE write_date(date, text)

    ! write date in constant format to standard output
    ! input date can be different declared

    IMPLICIT NONE
    INTRINSIC :: PRESENT, TRIM

    TYPE (time_days), INTENT(in)           :: date
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: text

    ! LOCAL
    CHARACTER(len=STRLEN_ULONG) :: date_mess1, date_mess2
    INTEGER                     :: status

    IF (p_parallel_io) THEN
       CALL print_date_components(date,status,date_mess1)
       CALL timer_message(status, ' ')

       IF (PRESENT(text)) THEN
          date_mess2 = TRIM(text)//' '//TRIM(date_mess1)
       ELSE
          date_mess2 = TRIM(date_mess1)
       END IF
       CALL timer_message(0,' ',date_mess2)
    ENDIF

  END SUBROUTINE write_date

  ! --------------------------------------------------------------------------

  SUBROUTINE write_date_debug(date, text)

    ! write date in constant format to standard output
    ! input date can be differently declared

    USE messy_main_constants_mem, ONLY: iouerr

    IMPLICIT NONE
    INTRINSIC :: PRESENT, TRIM

    TYPE (time_days), INTENT(in)           :: date
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: text

    ! LOCAL
    CHARACTER(len=STRLEN_ULONG) :: date_mess1, date_mess2
    INTEGER                     :: status

    IF (p_parallel_io) THEN
       CALL print_date_components(date,status,date_mess1)
       CALL timer_message(status, ' ')

       IF (PRESENT(text)) THEN
          date_mess2 = TRIM(text)//' '//TRIM(date_mess1)
       ELSE
          date_mess2 = TRIM(date_mess1)
       END IF
       write(iouerr,*) 'WRITE_DATE: ', TRIM(date_mess2)
   ENDIF

  END SUBROUTINE write_date_debug

  ! --------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  ! PUBLIC HELPER ROUTINES
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  SUBROUTINE timer_set_rerun_event(status, counter, unit, adjustment, offset)

    IMPLICIT NONE

    ! I/O
    INTEGER,  INTENT(OUT)        :: status
    INTEGER,  INTENT(IN)         :: counter    ! No. of steps in given unit
    CHARACTER(LEN=*), INTENT(IN) :: unit       ! counter unit type
    CHARACTER(LEN=*), INTENT(IN) :: adjustment ! adjustment in side the unit
    ! offset to initial date in seconds
    REAL(dp), INTENT(IN)         :: offset
    !
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'timer_set_rerun_event'

    IO_RERUN_EV%counter    = counter
    IO_RERUN_EV%unit       = unit
    IO_RERUN_EV%adjustment = adjustment
    IO_RERUN_EV%offset     = offset

    status = 0
    RETURN

  END SUBROUTINE timer_set_rerun_event

  ! ---------------------------------------------------------------------------

  SUBROUTINE timer_get_rerun_event(status, counter, unit, adjustment, offset &
       , ev_adjust)

    IMPLICIT NONE

    ! I/O
    INTEGER,  INTENT(OUT)         :: status
    INTEGER,  INTENT(OUT)         :: counter    ! No. of steps in given unit
    CHARACTER(LEN=*), INTENT(OUT) :: unit       ! counter unit type
    CHARACTER(LEN=*), INTENT(OUT) :: adjustment ! adjustment in side the unit
    ! offset to initial date in seconds
    REAL(dp), INTENT(OUT)         :: offset
    CHARACTER(LEN=*), INTENT(OUT) :: ev_adjust

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'timer_get_rerun_event'

    status = 1

    counter    = event_count(RERUN_EVENT)
    unit       = event_unit(RERUN_EVENT)
    adjustment = event_adjust(RERUN_EVENT)
    offset     = event_offset(RERUN_EVENT)

    ev_adjust = 'next'

    status = 0

  END SUBROUTINE timer_get_rerun_event

  ! ---------------------------------------------------------------------------

  SUBROUTINE main_timer_read_nml_cpl(status, iou)

    ! MODULE ROUTINE (SMIL)
    !
    ! READ NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! Author: Astrid Kerkweg, UNI-MZ, Jun 2009

    USE messy_main_tools,   ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit
    ! local variable required, POINTERs must NOT be part of a NAMELIST
    REAL(dp) :: delta_time

    NAMELIST /CPL/ CAL_TYPE, MODEL_START, MODEL_STOP, IO_RERUN_EV &
                 , delta_time, lresume, NO_CYCLES, LABORT &
                 , NO_DAYS, NO_STEPS, nmonth, lprintnewdate

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_timer_read_nml_cpl'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    MODEL_START(:)= 0
    MODEL_STOP(:) = 0

    ! default:
    delta_time    = -999._dp

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! CHECK NAMELIST
    CALL read_nml_close(substr, iou, modstr)

    clock(1)%delta_time = delta_time

    status = 0  ! no ERROR

  END SUBROUTINE main_timer_read_nml_cpl

  ! ---------------------------------------------------------------------------

  SUBROUTINE timer_message(status,substr, message)

    IMPLICIT NONE
    INTRINSIC :: PRESENT, TRIM

    INTEGER, INTENT(IN)          :: status
    CHARACTER(LEN=*), INTENT(IN) :: substr
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: message

    ! LOCAL
    CHARACTER(LEN=STRLEN_ULONG)  :: errorstr
    CHARACTER(LEN=*), PARAMETER  :: m_text= ' '

    SELECT CASE (status)
    CASE (0)
       ! no error
       IF (PRESENT(message)) THEN
          IF (.NOT. TRIM(message) == TRIM(m_text)) &
               CALL info_bi(TRIM(message), substr)
       ENDIF
       RETURN
    CASE (1)
       ! WARNING
       IF (PRESENT(message)) CALL warning_bi(TRIM(message),substr)
       RETURN
    CASE(3400)
       IF (PRESENT(message)) THEN
          CALL error_bi(TRIM(message), substr)
       ELSE
          CALL error_bi('ERROR 3400' , substr)
       ENDIF
       RETURN
    CASE(3434)
       errorstr = 'manager not initialized'
    CASE(3435)
       errorstr = 'date not initialized'
    CASE(3436)
       errorstr = 'messy_main_timer manager_reinit_days//&
            &// Time manager was not initialized, no reinit possible.'
    CASE(3437)
       errorstr = 'messy_main_timer:manager_step negative steps invalid'
    CASE(3438)
       errorstr = 'new date before start date'
    CASE(3439)
       errorstr = 'new date do not fit to a time step'
    CASE(3440)
       errorstr = 'combination of offset/unit invalid'
    CASE(3441)
       errorstr = 'Start date larger/equal than stop date ....'
    CASE(3442)
       errorstr = 'Current date larger than stop date ....'
    CASE(3443)
       errorstr = 'Unknown date ....'
    CASE(3444)
       errorstr = 'Unknown calendar type ....'
    CASE(3445)
       errorstr = 'Time difference between start and resume date negative...'
    CASE(3446)
       errorstr = 'Start date in future'
    CASE(3447)
       errorstr = 'Set start date to resume date, force initial run'
    CASE(3450)
       errorstr = 'Unit months not unambiguously convertable to seconds'
    CASE(3451)
       errorstr = 'Unit years not unambiguously convertable to seconds'
    CASE(3452)
       errorstr = 'dt in COSMO RUNCTL /= delta_time in TIMER namelist'
    CASE(3453)
       errorstr = 'dtime in ICON run_nml /= delta_time in TIMER namelist'
    CASE DEFAULT
       errorstr = ' '
    END SELECT

    CALL error_bi(errorstr, substr)

  END SUBROUTINE timer_message

  ! ----------------------------------------------------------------------------

! ****************************************************************************
END MODULE messy_main_timer_bi
! ****************************************************************************
