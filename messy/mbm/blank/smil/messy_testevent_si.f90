! **********************************************************************
MODULE messy_testevent_si
! **********************************************************************
! This submodel is based on the CHANNEL BMIL
! Authors: Astrid Kerkweg, UNI-BN, BONN (eventtest added)
! **********************************************************************

  ! BMIL
  USE messy_main_blather_bi,    ONLY: error_bi
  USE messy_main_mpi_bi,        ONLY: p_parallel_io, p_io, p_bcast
  USE messy_main_timer_bi,      ONLY: p_bcast_event, timer_event_init
  ! SMCL
  USE messy_main_timer_event,   ONLY: time_event, io_time_event &
                                    , TIME_INC_HOURS,TRIG_FIRST
  USE messy_main_constants_mem, ONLY: STRLEN_VLONG
  USE messy_testevent

  IMPLICIT NONE
  INTRINSIC :: NULL
  PRIVATE

  ! MODULE VARIABLES

  INTEGER, PARAMETER :: NMAXEVENTS = 100

  TYPE t_event_timer
     CHARACTER(LEN=STRLEN_VLONG) :: name     = ''       ! CHANNEL NAME
     TYPE(io_time_event)         :: io_event = &
          io_time_event(1, TIME_INC_HOURS,TRIG_FIRST,0) ! DEFAULT
  END TYPE t_event_timer

  TYPE(t_event_timer), DIMENSION(NMAXEVENTS), SAVE :: TIMER_EVENT
  TYPE(time_event),    DIMENSION(NMAXEVENTS), SAVE :: EVENT

  CHARACTER(LEN=8) :: evaldate

  PUBLIC :: testevent_init_coupling
  PUBLIC :: testevent_global_end

CONTAINS

  ! --------------------------------------------------------------------
  SUBROUTINE testevent_init_coupling

    USE messy_main_tools,    ONLY: find_next_free_unit

    IMPLICIT NONE

    INTEGER :: i, iou, status
    CHARACTER(LEN=*), PARAMETER :: substr = 'testevent_init_coupling'

    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL testevent_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF

    CALL p_bcast(evaldate, p_io)
    DO i=1, NMAXEVENTS
       CALL p_bcast(TIMER_EVENT(i)%name, p_io)
       IF (TRIM(TIMER_EVENT(i)%name) == '') CYCLE
       CALL p_bcast_event(TIMER_EVENT(i)%io_event, p_io)
       write (*,*) 'NNL CPL ',i, TRIM(TIMER_EVENT(i)%name) &
            , TIMER_EVENT(i)%io_event%counter &
            , TIMER_EVENT(i)%io_event%unit &
            , TIMER_EVENT(i)%io_event%adjustment &
            , TIMER_EVENT(i)%io_event%offset 

       CALL timer_event_init(EVENT(i), TIMER_EVENT(i)%io_event &
            , TRIM(TIMER_EVENT(i)%name), TRIM(evaldate) )  
    END DO

  END SUBROUTINE testevent_init_coupling
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE testevent_global_end
    
    USE messy_main_timer_bi,    ONLY: event_state
    USE messy_main_timer,       ONLY: next_date, current_date, time_days &
                                    , print_date_components
    USE messy_main_tools,       ONLY: int2str

    IMPLICIT NONE

    INTEGER          :: i, status
    LOGICAL          :: ltrig
    TYPE(time_days)  :: date 
    CHARACTER(LEN=3)  :: str
    
    IF (TRIM(evaldate) == 'present') THEN
       date = current_date
    ELSE
       date = next_date
    END IF

    IF ( p_parallel_io) write (*,*) '********************************************************'
    CALL print_date_components(date, status, info='TESTEVENT: ')
    DO i=1, NMAXEVENTS
       IF (TRIM(TIMER_EVENT(i)%name) == '') CYCLE
       CALL int2str(str, i)
       write (*,*) '--------------------------------'
       write (*,*) 'TEV '//STR//' '//TRIM(TIMER_EVENT(i)%name)
       ltrig =  event_state(EVENT(i), date)

       IF (ltrig) THEN
          CALL print_date_components(date, status, info='TEV '//STR//': '//TRIM(TIMER_EVENT(i)%name)//':    ')
       END IF
       write (*,*) '--------------------------------'
    END DO
    write (*,*) '********************************************************'

  END SUBROUTINE testevent_global_end
  ! -------------------------------------------------------------------


  ! -------------------------------------------------------------------
  ! -------------------------------------------------------------------
  ! PRIVATE
  ! -------------------------------------------------------------------
  ! -------------------------------------------------------------------
  SUBROUTINE testevent_read_nml_cpl(status, iou)

    ! MODULE ROUTINE (SMIL)
    !
    ! READ NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! Author: Patrick Joeckel, MPICH, Dec 2004

    USE messy_main_tools,   ONLY: read_nml_open, read_nml_check, read_nml_close
    USE messy_testevent,    ONLY: modstr

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    NAMELIST /CPL/  evaldate, TIMER_EVENT

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'testevent_read_nml_cpl'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    ! -> OTHER DEFAULT VALUES ARE SET AT DECLARATION ABOVE

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE testevent_read_nml_cpl
  ! -------------------------------------------------------------------

! **********************************************************************
END MODULE messy_testevent_si
! **********************************************************************
