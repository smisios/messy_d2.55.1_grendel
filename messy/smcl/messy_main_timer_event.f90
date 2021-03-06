! ****************************************************************************
MODULE messy_main_timer_event
! ****************************************************************************

!#define LDEBUG

  !--------------------------------------------------------------------------
  ! EVENT MANAGER
  !--------------------------------------------------------------------------

  USE messy_main_constants_mem, ONLY: dp, OneDay, STRLEN_SHORT    &
                                    , STRLEN_MEDIUM, STRLEN_ULONG 
  USE messy_main_tools,         ONLY: int2str
  USE messy_main_timer

  IMPLICIT NONE
  PRIVATE
  SAVE

  !
  ! TIME_INC_*     predefined counter units
  !
  CHARACTER(len=*), PUBLIC, PARAMETER :: &
       TIME_INC_MILSECS = 'milsecs'  ,&!  milliseconds
       TIME_INC_SECONDS = 'seconds'  ,&!
       TIME_INC_MINUTES = 'minutes'  ,&!
#ifndef _LINUX64
       TIME_INC_HOURS   = 'hours'    ,&!
       TIME_INC_DAYS    = 'days'     ,&!
       TIME_INC_MONTHS  = 'months'   ,&!
       TIME_INC_YEARS   = 'years'      !
#else
! workaround for compiler bug on opteron
       TIME_INC_HOURS   = 'hours  '  ,&!
       TIME_INC_DAYS    = 'days   '  ,&!
       TIME_INC_MONTHS  = 'months '  ,&!
       TIME_INC_YEARS   = 'years  '    !
#endif

  INTEGER, PARAMETER, PUBLIC :: I_INC_MILSECS = 1
  INTEGER, PARAMETER, PUBLIC :: I_INC_SECONDS = 2
  INTEGER, PARAMETER, PUBLIC :: I_INC_MINUTES = 3
  INTEGER, PARAMETER, PUBLIC :: I_INC_HOURS   = 4
  INTEGER, PARAMETER, PUBLIC :: I_INC_DAYS    = 5
  INTEGER, PARAMETER, PUBLIC :: I_INC_MONTHS  = 6
  INTEGER, PARAMETER, PUBLIC :: I_INC_YEARS   = 7

  INTEGER,  PARAMETER :: CONVERT(5) = (/0,1,60,3600,84600/) 
  REAL(dp), PARAMETER :: CONVERTMS(5) = &
       (/0.001_dp,1._dp,60._dp,3600._dp,84600._dp/) 

  CHARACTER(len=*), PUBLIC, PARAMETER :: &
#ifndef _LINUX64
       TRIG_FIRST  = 'first'  ,&! trigger in first step of counter unit
       TRIG_LAST   = 'last'   ,&! trigger in last step of counter unit
       TRIG_EXACT  = 'exact'  ,&! trigger without adjustment inside the unit
       TRIG_NONE   = 'off'      ! event trigger non active
#else
! workaround for compiler bug on opteron
       TRIG_FIRST  = 'first'  ,&! trigger in first step of counter unit
       TRIG_LAST   = 'last '  ,&! trigger in last step of counter unit
       TRIG_EXACT  = 'exact'  ,&! trigger without adjustment in side the unit
       TRIG_NONE   = 'off  '    ! event trigger non active
#endif
  !
  INTEGER, PARAMETER, PUBLIC :: I_TRIG_FIRST = 1
  INTEGER, PARAMETER, PUBLIC :: I_TRIG_LAST  = 2
  INTEGER, PARAMETER, PUBLIC :: I_TRIG_EXACT = 3
  INTEGER, PARAMETER, PUBLIC :: I_TRIG_NONE  = 4


  CHARACTER(len=*), PARAMETER, PUBLIC :: &
  ! the adjustment of events for RERUN is dependent on the trigger step
  ! the trigger step can be the present date or the next date
       EV_TLEV_PRES    = 'present'        ,&! check event with present date
       EV_TLEV_NEXT    = 'next'           ,&! check event with next date
       EV_TLEV_PREV    = 'previous'       ,&! check event with previous date
       TIME_INC_STEPS  = 'steps'          ,&! special event interval unit
       TIME_INC_ALWAYS = 'always'           ! special event always used  

  INTEGER, PARAMETER, PUBLIC :: I_EV_TLEV_PRES = 1
  INTEGER, PARAMETER, PUBLIC :: I_EV_TLEV_NEXT = 2
  INTEGER, PARAMETER, PUBLIC :: I_EV_TLEV_PREV = 3

  ! structures
  TYPE, PUBLIC :: io_time_event            
     ! external given event properties
     ! - No. of steps in given unit
     INTEGER     :: counter    = 0           
     ! - counter unit type
     CHARACTER(len=STRLEN_SHORT) :: unit  = TIME_INC_SECONDS 
     ! - adjustment inside the unit
     CHARACTER(len=STRLEN_SHORT) :: adjustment = TRIG_EXACT  
     ! - offset to initial date in milli (!) seconds
     REAL(dp)     :: offset     = 0.0_dp
  END TYPE io_time_event

  TYPE, PUBLIC :: time_event          

     !   hold all event relevant informations
     PRIVATE
     LOGICAL    :: init   = .FALSE.    ! event state for access control
     LOGICAL    :: active = .FALSE.    ! event active true=on false=off
     INTEGER    :: count  = 1          ! increment between the action
!     INTEGER    :: offset = 0          ! offset (seconds)
     REAL(dp)   :: offset = 0.0_dp      ! offset (seconds)
     REAL(dp)   :: delta      = 2.0_dp ! adjustment interval in seconds
     REAL(dp)   :: half_delta = 1.0_dp
     ! descriptive label of the event
     CHARACTER(len=STRLEN_MEDIUM) :: label  = '' 
     ! tag to distinguish events with same label == meaning for different clocks
     INTEGER                      :: tag  = -99
     ! name of the basic units
     CHARACTER(len=STRLEN_SHORT)  :: unit   = TIME_INC_SECONDS 
     CHARACTER(len=STRLEN_SHORT)  :: adjust = TRIG_EXACT ! type of triggering

     TYPE (time_days)   ::  initial_date     ! initial data for event
     TYPE (time_days)   ::    cycle_date     ! without offset
     TYPE (time_days)   :: previous_trigger  ! previous trigger date
     TYPE (time_days)   ::  current_trigger  ! current trigger date
     TYPE (time_days)   ::     next_trigger  ! next trigger date
     
  END TYPE time_event

  ! INTERFACE SUBROUTINES AND FUNCTIONS

  PUBLIC :: event_init        !   initialise parts of an event
  INTERFACE event_init
    MODULE PROCEDURE event_init
                            ! (IO:event,I:name,I:count,I:unit,I:adj[,I:delta])
    MODULE PROCEDURE event_set_day    
                            ! (IO:event,I:time_days,I:lshift,O:status[,O:message])
    MODULE PROCEDURE event_toggle ! (IO:event,I:adjust,O:status,O:message)
  END INTERFACE

  PUBLIC :: event_reinit     ! reset initial trigger dates
                             ! (IO:event,I:time_days)

  PUBLIC  :: print_event_name    ! (I:event) returns event name 
  PUBLIC  :: event_initial_date  ! (I:event,O:initial_date)
  PUBLIC  :: event_cycle_date    ! (I:event,O:cycle_date)
  PUBLIC  :: event_previous_date ! (I:event,O:previous_trigger)
  PUBLIC  :: event_current_date  ! (I:event,O:current_trigger)
  PUBLIC  :: event_next_date     ! (I:event,O:next_trigger)

  PUBLIC  :: event_is_active
  PUBLIC  :: event_is_init
  PUBLIC  :: event_delta
  PUBLIC  :: event_half_delta
  PUBLIC  :: event_adjust
  PUBLIC  :: event_label
  PUBLIC  :: event_name 
  PUBLIC  :: event_tag
  PUBLIC  :: event_tagstr
  PUBLIC  :: event_unit
  PUBLIC  :: event_count
  PUBLIC  :: event_offset

  PUBLIC  :: convert_steps2unit
  PUBLIC  :: convert_unit2seconds

  PRIVATE :: get_next_trigger ! (I:event) find next trigger

CONTAINS

  SUBROUTINE event_init (event, name, count, unit, adjust, delta, offset &
       , status, message, tag)

    ! initialize elements of an event structure

    IMPLICIT NONE
    INTRINSIC :: ADJUSTL, INT, LEN, MIN, MOD, PRESENT, REAL, TRIM

    TYPE (time_event), INTENT(inout) :: event   ! new event
    CHARACTER(len=*),  INTENT(in)    :: name    ! name of event
    INTEGER,           INTENT(inout) :: count   ! no. of units between events
    CHARACTER(len=*),  INTENT(inout) :: unit    ! unit
    CHARACTER(len=*),  INTENT(in)    :: adjust  ! type of adjustment
    REAL(dp), OPTIONAL,INTENT(in)    :: delta   ! adjustment interval [seconds]
    REAL(dp), OPTIONAL, INTENT(in)   :: offset  ! offset relativ to initial date
    INTEGER, OPTIONAL, INTENT(OUT)   :: status    ! error/message status
    CHARACTER(LEN=STRLEN_ULONG), OPTIONAL, INTENT(OUT) :: message ! message text
    INTEGER, OPTIONAL, INTENT(IN)    :: tag

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='main_timer: event_init'
    INTEGER                     :: is, isec1, isec2, ioffset
    REAL(dp)                    :: smallest
    CHARACTER(LEN=STRLEN_ULONG) :: m_text

    IF (PRESENT(status)) status = 0
    IF (PRESENT(message)) message = ' '

    IF (event%init) THEN
       IF (PRESENT(status)) status = 1
       IF (PRESENT(message)) THEN
          write (message,*) substr,' Event was initialized, nothing changed'
          RETURN
       ELSE
          write (*,*) substr,' Event was initialized, nothing changed'
       ENDIF
    ELSE
       DO is=1,STRLEN_MEDIUM
          event%label(is:is) = ' '
       END DO
       is = MIN(LEN(TRIM(name)),STRLEN_MEDIUM)

       event%label(1:is) = name(1:is)

       IF (PRESENT(tag)) THEN
          event%tag = tag
       END IF

       ! evaluate offset 
       event%offset = 0.0_dp
       IF (PRESENT(offset)) event%offset = offset
       IF ( ABS(event%offset) > TINY(0._dp) ) THEN 
          IF (TRIM(adjust) /= TRIM(TRIG_EXACT)) THEN

             SELECT CASE(TRIM(ADJUSTL(unit)))
             CASE(TRIM(TIME_INC_MONTHS),TRIM(TIME_INC_YEARS))
                IF (PRESENT(status)) THEN
                   status = 3440
                ELSE
                   write (*,*) substr, 'combination of offset/unit invalid'
                END IF
             END SELECT

             ! correct unit, dependent on the offset
             ioffset = NINT(event%offset)
             IF (event%offset -REAL(ioffset,dp) < TINY(0._dp)) THEN
                ! NO MILLISECONDS
                IF (MOD(ioffset,CONVERT(I_INC_DAYS)) == 0) THEN 
                   ! it's fine
                ELSE IF (MOD(ioffset,CONVERT(I_INC_HOURS)) == 0) THEN
                   SELECT CASE(TRIM(unit))
                   CASE(TRIM(TIME_INC_DAYS))
                      unit = TIME_INC_HOURS;   count = count*24
                   END SELECT
                ELSE IF (MOD(ioffset,CONVERT(I_INC_MINUTES)) == 0) THEN
                   SELECT CASE(TRIM(unit))
                   CASE(TRIM(TIME_INC_DAYS))
                      unit = TIME_INC_MINUTES
                      count = count*24*60
                   CASE(TRIM(TIME_INC_HOURS))
                      unit = TIME_INC_MINUTES
                      count = count*60
                   END SELECT
                ELSE 
                   SELECT CASE(TRIM(unit))
                   CASE(TRIM(TIME_INC_DAYS))
                      unit = TIME_INC_SECONDS
                      count = count*24*3600
                   CASE(TRIM(TIME_INC_HOURS))
                      unit = TIME_INC_SECONDS
                      count = count*3600
                   CASE(TRIM(TIME_INC_MINUTES))
                      unit = TIME_INC_SECONDS
                      count = count*60
                   END SELECT
                END IF

             ELSE 

                SELECT CASE(TRIM(unit))
                CASE(TRIM(TIME_INC_DAYS))
                   count = count*24*3600*1000
                CASE(TRIM(TIME_INC_HOURS))
                   count = count*3600*1000
                CASE(TRIM(TIME_INC_MINUTES))
                   count = count*60*1000
                CASE(TRIM(TIME_INC_SECONDS))
                   count = count*1000
                END SELECT
                unit  = TIME_INC_MILSECS

             END IF
          END IF
       END IF
       event%count = count

       ! define the smallest delta needed for adjustment 
       SELECT CASE(TRIM(unit))
       CASE(TRIM(TIME_INC_MILSECS))
          smallest = 0.0015_dp
       CASE(TRIM(TIME_INC_SECONDS))
          smallest = 1.5_dp
       CASE(TRIM(TIME_INC_MINUTES))
          smallest = 60.0_dp
       CASE(TRIM(TIME_INC_HOURS))
          smallest = 3600.0_dp
       CASE(TRIM(TIME_INC_DAYS))
          smallest = 24.0_dp*3600.0_dp
       CASE(TRIM(TIME_INC_MONTHS))
          smallest = 28.0_dp*24.0_dp*3600.0_dp
       CASE(TRIM(TIME_INC_YEARS))
          smallest = 360.0_dp*28.0_dp*24.0_dp*3600.0_dp
       CASE default
          m_text = 'Counter unit unknown ::' // TRIM(unit)
          IF (PRESENT(status)) status = 3400
          IF (PRESENT(message)) THEN
             message = m_text
          ELSE
             write (*,*) substr, m_text
          ENDIF
          RETURN
       END SELECT
       IF (smallest > 10._dp) THEN
          smallest   = smallest * REAL(count,dp) - 1.0_dp
       ELSE
          smallest   = smallest * REAL(count,dp) - 0.00001_dp
       END IF
       event%unit = unit

       ! valuate epsilon interval and adjustment      
       SELECT CASE(TRIM(adjust))
       CASE (TRIM(TRIG_FIRST), TRIM(TRIG_LAST),TRIM(TRIG_EXACT),TRIM(TRIG_NONE))
          IF (PRESENT(delta)) THEN
             IF (0._dp < delta .AND. delta <= smallest) THEN
                event%delta = delta
             ELSE
                event%delta = smallest
                WRITE(m_text,*) 'Preset adjustment interval (a) ::',smallest
                IF (PRESENT(status)) status = 1
                IF  (PRESENT(message)) THEN
                   message = TRIM(m_text)
                ELSE
                   WRITE (*,*) substr, m_text
                ENDIF
             END IF

          ELSE
             event%delta = smallest
             WRITE(m_text,*) 'Preset adjustment interval (b) ::',smallest
             IF (PRESENT(message)) THEN
                message = m_text
             ELSE
                write (*,*) substr,m_text 
             ENDIF

          END IF

       CASE default
          m_text = 'Event adjustment not defined ::' // TRIM(adjust)

          IF (PRESENT(status)) status = 3400
          IF (PRESENT(message)) THEN
             message = m_text
          ELSE
             write (*,*) substr, m_text
          ENDIF
          RETURN

       END SELECT

       SELECT CASE(TRIM(adjust))
       CASE (TRIM(TRIG_NONE))
          event %adjust = TRIG_EXACT
          event %active = .FALSE.
       CASE default
          event %adjust = TRIM(adjust)
          event %active = .TRUE.
       END SELECT

       ! calculates half delta in milli seconds     
       IF (event%delta < 0.001_dp) THEN
          IF (PRESENT(status)) status = 3400
          IF (PRESENT(message)) THEN
             message = m_text
          ELSE
             write (*,*) substr, m_text
          ENDIF
          RETURN
       ENDIF
       IF ( ABS(REAL(INT(event%delta),dp)-event%delta) < TINY(1._dp)) THEN
          isec1 = INT(0.5_dp * event%delta)
          isec2 = INT(event%delta) - isec1
          event%half_delta = REAL(isec1,dp)
          IF (isec1 < isec2 ) event%half_delta = REAL(isec2,dp)
       ELSE
          isec1 = INT(0.5_dp * event%delta*1000._dp)
          isec2 = INT(event%delta*1000._dp) - isec1
          event%half_delta = REAL(isec1,dp)/1000._dp
          IF (isec1 < isec2 ) event%half_delta = REAL(isec2,dp)/1000._dp
       ENDIF

       ! set initial date to calendar start     
       CALL date_set (0,0._dp,event %initial_date)
       CALL date_set (0,0._dp,event %cycle_date)
       CALL date_set (0,0._dp,event %current_trigger)
       CALL date_set (0,0._dp,event %previous_trigger)
       CALL date_set (0,0._dp,event %next_trigger)

       event %init    = .TRUE.
    END IF

  END SUBROUTINE event_init

! -------------------------------------------------------------------

  SUBROUTINE event_toggle (event, adjust, status, message)  

    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! reset the adjustment type
    ! used for activate/deactive an event

    TYPE (time_event), INTENT(inout) :: event   ! new event
    CHARACTER(len=*),  INTENT(in)    :: adjust  ! type of adjustment
    INTEGER,                     INTENT(OUT) :: status
    CHARACTER(LEN=STRLEN_ULONG), INTENT(OUT) :: message

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='main_timer: event_toggle'

    status    = 0
    message = ' '

    SELECT CASE(TRIM(adjust))
    CASE (TRIM(TRIG_FIRST), TRIM(TRIG_LAST), TRIM(TRIG_EXACT))
       write (message,*) substr &
            ,' Reset event adjustment for ... '// TRIM(event_name(event))
       event%adjust = TRIM(adjust)
       event%active = .TRUE.

    CASE (TRIG_NONE)
       write (message,*) 'Deactivate event ... '// TRIM(event_name(event))
       event%active = .FALSE.

    CASE default
       write (message,*) 'Event adjustment not defined <' // TRIM(adjust) // '>'
       status    = 3400
       RETURN
    END SELECT

  END SUBROUTINE event_toggle

  ! ----------------------------------------------------------------------------

  SUBROUTINE event_set_day (event, date, lshift, status, message) 

    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! set first trigger date

    TYPE (time_event), INTENT(inout) :: event
    TYPE (time_days),  INTENT(in)    :: date
    LOGICAL,           INTENT(in)    :: lshift
    INTEGER,           INTENT(OUT)   :: status
    CHARACTER(LEN=STRLEN_ULONG), INTENT(OUT) ::  message

    ! LOCAL
    !CHARACTER(LEN=*), PARAMETER :: substr='main_timer: event_set_day'

    status = 0
    message = ' '

    IF (.NOT. event%init) THEN
       message = 'Initialize event <' // TRIM(event_name(event)) // '>first'
       status =3400
    ELSE
       IF (lshift) THEN
          CALL copy_date(event% current_trigger,event% previous_trigger,status)
          IF (status /= 0) RETURN
       ELSE
          CALL copy_date (date,event%  initial_date,status)
          IF (status /= 0) RETURN
          CALL copy_date (date,event%  cycle_date, status)
          IF (status /= 0) RETURN
          CALL copy_date (date,event%  previous_trigger,status)
          IF (status /= 0) RETURN
       END IF

       CALL copy_date(date, event% current_trigger,status)
       IF (status /= 0) RETURN
       CALL get_next_trigger(event, status, message)
       IF (status /= 0 ) RETURN

    END IF
    status = 0
    message = '  '

  END SUBROUTINE event_set_day

  ! ---------------------------------------------------------------------------

  SUBROUTINE event_reinit (event, date, status, message) 

    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! reset first trigger date without evaluation of next trigger

    TYPE (time_event), INTENT(inout) :: event
    TYPE (time_days),  INTENT(in)    :: date
    INTEGER,           INTENT(OUT)   :: status
    CHARACTER(LEN=STRLEN_ULONG), INTENT(OUT) ::  message

    !LOCAL
    !CHARACTER(LEN=*), PARAMETER :: substr='TIMER: event_reinit'
    LOGICAL :: learlier
   
    status = 0
    message = ' ' 

    CALL if_less(event%next_trigger,date,learlier, status)
    IF (status /= 0) RETURN

    IF (learlier) THEN
      message = 'New trigger date of <' // TRIM(event_name(event)) // &
           '> in the past. Reset not possible.'
      status = 3400
      RETURN
    END IF
    
    CALL copy_date(date, event% previous_trigger, status)
    IF ( status /= 0 ) RETURN
    CALL copy_date(date, event% current_trigger, status)
    IF ( status /= 0 ) RETURN

  END SUBROUTINE event_reinit

  ! ---------------------------------------------------------------------------

  FUNCTION print_event_name (ev) RESULT (str)  !*****************************

    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! return event name
    !-
    TYPE(time_event), INTENT(in)   :: ev
    CHARACTER(len=STRLEN_MEDIUM+4) :: str
    CHARACTER(len=2) :: tagstr
    
    tagstr= event_tagstr(ev)
    IF (TRIM(tagstr) == '') THEN
       str = TRIM(ev%label)
    ELSE
       str = TRIM(ev%label)//'_T'//TRIM(tagstr)
    END IF

  END FUNCTION print_event_name

  ! ---------------------------------------------------------------------------

  FUNCTION print_event_label (ev) RESULT (str)  !*****************************

    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! return event label
    !-
    TYPE(time_event), INTENT(in) :: ev
    CHARACTER(len=STRLEN_MEDIUM) :: str

    str = TRIM(ev%label)

  END FUNCTION print_event_label

  ! ---------------------------------------------------------------------------

  SUBROUTINE get_next_trigger (event, status, message)

    ! calculates next trigger from present trigger

    IMPLICIT NONE
    INTRINSIC :: TRIM

    TYPE (time_event), INTENT(inout) :: event
    INTEGER,           INTENT(OUT)   :: status
    CHARACTER(LEN=STRLEN_ULONG), INTENT(OUT) :: message
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='main_timer get_next_trigger'
    TYPE (time_days)  :: my_date
    LOGICAL           :: first_call, learly

    status    = 0
    message = ' '

    ! evaluate if called for the first time
    CALL if_equal(event%previous_trigger, event%current_trigger &
         ,first_call,status)
    IF (status /=0) RETURN

    ! calculates the end of the next time interval
#ifdef LDEBUG
    write (*,*) '** cycle date 01 ** ', first_call
    CALL print_date_components (event%cycle_date, status &
         , info='** cycle date 01 ** ')
#endif

    CALL copy_date(event%cycle_date, my_date,status)
    IF (status /=0) RETURN

#ifdef LDEBUG
    CALL print_date_components (my_date, status, info= '** my_date 01 ** ')
#endif

    IF ( (.NOT. first_call) .OR. (event%offset<=event%delta)) THEN
       CALL get_next_date (event%unit, event%count, my_date, status)
       IF (status /=0) RETURN
    ENDIF

#ifdef LDEBUG
    CALL print_date_components (my_date, status, info= '** my_date 02 ** ')
    CALL print_date_components (event% current_trigger, status &
         , info='** current trigger 01 ** ' )
    CALL print_date_components (event% previous_trigger, status &
         , info= '** previous trigger 01 ** ' )
#endif

    ! adjust the interval date
    CALL adjust_date(first_call, event%adjust, event%unit, event%delta &
         , my_date, status)
    IF (status /=0) RETURN

#ifdef LDEBUG
    CALL print_date_components (my_date, status, info='** my_date 03 ** ')
#endif

    ! store the adjusted interval end
    CALL copy_date(my_date, event%cycle_date, status)
    IF (status /=0) RETURN

#ifdef LDEBUG
    CALL print_date_components (event%cycle_date, status&
         , info='** cycle_date 02 ** ')
#endif

    ! perform the offset calculation
    CALL add_date (0, event%offset, my_date,status)
    IF (status /=0) RETURN

#ifdef LDEBUG
    CALL print_date_components (my_date, status, info='** my_date 05 ** ')
#endif

    CALL copy_date (my_date,event%next_trigger,status) 
    IF (status /=0) RETURN

#ifdef LDEBUG
    CALL print_date_components (event%next_trigger, status &
         , info='** next trigger 01 ** ')
#endif

    CALL if_less(event%current_trigger, event%next_trigger, learly,status)
    IF (status /=0) RETURN

    IF (.NOT. learly) THEN
       write (*,*) 'Event <'// TRIM(event_name(event))//'> [current date > next date]'
       CALL print_date_components (event% current_trigger, status &
            , info='****** CURRENT TRIGGER ********')
       IF (status /=0) RETURN
       CALL print_date_components (event%    next_trigger, status &
            , info='****** NEXT TRIGGER ********')
       STATUS = 3400
       write (message,*) substr,' can not correctly evaluate next_trigger'
       RETURN
    END IF

    message = '  '
    status = 0

  CONTAINS

    SUBROUTINE get_next_date (unit, count, my_date, status)

      IMPLICIT NONE
      INTRINSIC :: INT, MIN, MOD, REAL

      CHARACTER (len=*),  INTENT(in)    :: unit
      INTEGER,            INTENT(in)    :: count
      TYPE (time_days), INTENT(inout)   :: my_date 
      INTEGER,            INTENT(OUT)   :: status

      !LOCAL
      INTEGER        :: yr, mo, dy, hr, mn, se, iday
      INTEGER        :: ms
      REAL (kind=dp) :: zsec, rsec

      status = 0
      
      CALL timer_get_date(status, my_date, yr, mo, dy, hr, mn, se, ms)
      IF (status /=0) RETURN

      SELECT CASE(TRIM(unit))
      CASE(TRIM(TIME_INC_MILSECS), TRIM(TIME_INC_SECONDS) &
           , TRIM(TIME_INC_MINUTES), TRIM(TIME_INC_HOURS) , TRIM(TIME_INC_DAYS))

        SELECT CASE(TRIM(unit))       ! calculate increment 
        CASE(TRIM(TIME_INC_MILSECS))
           zsec = 0.001_dp
        CASE(TRIM(TIME_INC_SECONDS))
           zsec = 1._dp
        CASE(TRIM(TIME_INC_MINUTES))
           zsec = 60._dp
        CASE(TRIM(TIME_INC_HOURS))
           zsec = 3600._dp
        CASE(TRIM(TIME_INC_DAYS))
           zsec = 24._dp*3600._dp
        END SELECT

        ! convert seconds into days and seconds
        zsec = zsec * REAL(count,dp)
        iday = INT((zsec+0.00001_dp)/OneDay)
        rsec = zsec - REAL(iday,dp)*OneDay

        ! add  day and time to last adjustment trigger date
        CALL add_date(iday, rsec, my_date,status)
        IF (status /= 0) RETURN

      CASE(TIME_INC_MONTHS)
        mo = mo + count
        yr = yr + INT((REAL(mo,dp)-0.5_dp)/12._dp)
        mo = MOD(mo,12); IF (mo == 0) mo = 12

        SELECT CASE(CAL_TYPE)
        CASE (CAL_JULIAN)    
           dy = MIN(dy,JulianMonthLength(yr,mo))
        CASE (CAL_360D)    
           dy = MIN(dy,30)  ! Month length in  360 year
        END SELECT
        CALL timer_set_date(status, my_date, yr, mo, dy, hr, mn, se, ms)

      CASE(TIME_INC_YEARS)
        yr = yr + count
        SELECT CASE(CAL_TYPE)
        CASE (CAL_JULIAN)
           dy = MIN(dy,JulianMonthLength(yr,mo))
        CASE (CAL_360D)    
           dy = MIN(dy,30) ! Ly360
        END SELECT
        CALL timer_set_date(status, my_date, yr, mo, dy, hr, mn, se, ms)

      END SELECT

    END SUBROUTINE get_next_date

    ! ------------------------------------------------------

    SUBROUTINE adjust_date (first_call, adjust, unit, delta, my_date, status)

      IMPLICIT NONE

      LOGICAL,            INTENT(in)    :: first_call
      CHARACTER (len=*),  INTENT(in)    :: adjust, unit
      REAL (kind=dp),     INTENT(in)    :: delta
      TYPE (time_days), INTENT(inout)   :: my_date 
      INTEGER,            INTENT(OUT)   :: status

      ! LOCAL
      INTEGER  :: yr, mo, dy, hr, mn, se, iday
      INTEGER  :: ms
      REAL(dp) :: rsec

      STATUS = 0

#ifdef LDEBUG
      write (*,*) 'adjust ' ,adjust, 'unit ', unit, 'delta ', delta
      CALL print_date_components (my_date, status)
#endif

      SELECT CASE(TRIM(adjust))
      CASE(TRIM(TRIG_FIRST))

         CALL timer_get_date(status, my_date, yr, mo, dy, hr, mn, se, ms)
         IF (status /= 0) RETURN
#ifdef LDEBUG
         write (*,*) 'adjust_date 1a ', yr, mo, dy, hr, mn, se, ms
#endif
         ! adjust next trigger date to first delta environment in unit
         ! interval = <first second, first second + delta>

         SELECT CASE(TRIM(unit))
         CASE(TRIM(TIME_INC_SECONDS))
            ms = 0
         CASE(TRIM(TIME_INC_MINUTES))
            se = 0; ms = 0
         CASE(TRIM(TIME_INC_HOURS))
            mn = 0; se = 0; ms = 0
         CASE(TRIM(TIME_INC_DAYS))
            hr = 0; mn = 0; se = 0; ms = 0
         CASE(TRIM(TIME_INC_MONTHS))
            dy=1; hr = 0; mn = 0; se = 0; ms = 0
         CASE(TRIM(TIME_INC_YEARS))
            mo = 1; dy=1; hr = 0; mn = 0; se = 0; ms = 0
         END SELECT

#ifdef LDEBUG
         write (*,*) 'adjust_date 1b ', yr, mo, dy, hr, mn, se, ms
#endif

         CALL timer_set_date (status, my_date,yr, mo, dy, hr, mn, se, ms)

#ifdef LDEBUG
         CALL print_date_components (my_date, status)
#endif

      CASE(TRIM(TRIG_LAST))

#ifdef LDEBUG
         write (*,*) ' first_call', first_call
#endif
         IF (first_call) THEN ! skip back (one fraction - delta_adjustment)

            iday = 0
            rsec = delta
            CALL add_date (iday, rsec, my_date,status)
            IF (status /= 0) RETURN
            CALL timer_get_date   (status, my_date, yr, mo, dy, hr, mn, se, ms)
            IF (status /= 0) RETURN
#ifdef LDEBUG
            write (*,*) 'adjust_date 2 ', yr, mo, dy, hr, mn, se, ms
#endif

            SELECT CASE(TRIM(unit))
            CASE(TRIM(TIME_INC_MILSECS))
               rsec = -0.001_dp
               CALL add_date (iday, rsec, my_date,status)
#ifdef LDEBUG
               CALL print_date_components (my_date, status, info='adj milsecs')
#endif            
               CALL timer_get_date(status, my_date, yr, mo, dy, hr, mn, se, ms)
               IF (status /= 0) RETURN
#ifdef LDEBUG
               write (*,*) 'adjust_date milsecs ', yr, mo, dy, hr, mn, se, ms, rsec
#endif 
            CASE(TRIM(TIME_INC_SECONDS))
               rsec = -1.000_dp
               CALL add_date (iday, rsec, my_date,status)
#ifdef LDEBUG
               CALL print_date_components (my_date, status, info='adj seconds')
#endif            
               CALL timer_get_date(status, my_date, yr, mo, dy, hr, mn, se, ms)
               IF (status /= 0) RETURN
#ifdef LDEBUG
               write (*,*) 'adjust_date secs ', yr, mo, dy, hr, mn, se, ms, rsec
#endif 
            CASE(TRIM(TIME_INC_MINUTES))
               rsec =   -60._dp
               CALL add_date (iday, rsec, my_date,status)
               IF (status /= 0) RETURN
#ifdef LDEBUG
               CALL print_date_components (my_date, status, info='adj mins')
#endif            
               CALL timer_get_date(status, my_date, yr, mo, dy, hr, mn, se, ms)
               IF (status /= 0) RETURN
#ifdef LDEBUG
               write (*,*) 'adjust_date min ', yr, mo, dy, hr, mn, se, ms, rsec
#endif 
            CASE(TRIM(TIME_INC_HOURS))
               rsec = -3600._dp
               CALL add_date (iday, rsec, my_date,status)
               IF (status /= 0) RETURN
#ifdef LDEBUG
               CALL print_date_components (my_date, status, info='adj hours')
#endif
               CALL timer_get_date  (status, my_date, yr, mo, dy, hr, mn, se, ms)
               IF (status /= 0) RETURN
#ifdef LDEBUG
               write (*,*) 'adjust_date hr ', yr, mo, dy, hr, mn, se, ms, rsec, iday
#endif

            CASE(TRIM(TIME_INC_DAYS))
               iday =    -1
               rsec = 0._dp
               CALL add_date (iday, rsec, my_date,status)
               IF (status /= 0) RETURN
#ifdef LDEBUG
               CALL print_date_components (my_date, status, info='adj days')
#endif
               CALL timer_get_date(status, my_date, yr, mo, dy, hr, mn, se, ms)
               IF (status /= 0) RETURN
#ifdef LDEBUG
               write (*,*) 'adjust_date days ', yr, mo, dy, hr, mn, se, ms, rsec, iday
#endif

            CASE(TIME_INC_MONTHS)           ! skip back one month
               IF (mo == 1) THEN
                  mo = 12; yr = yr - 1
               ELSE
                  mo = mo - 1
               END IF
#ifdef LDEBUG
               write (*,*) 'adjust_date month ', yr, mo, dy, hr, mn, se, ms
#endif

            CASE(TIME_INC_YEARS)
               yr = yr - 1
#ifdef LDEBUG
               write (*,*) 'adjust_date year ', yr, mo, dy, hr, mn, se, ms
#endif
            END SELECT

         ELSE
            CALL timer_get_date (status, my_date, yr, mo, dy, hr, mn, se, ms)
            IF (status /= 0) RETURN
#ifdef LDEBUG
            write (*,*) 'adjust_date 9 ', yr, mo, dy, hr, mn, se, ms
#endif
         END IF

         ! adjust next trigger date to last environment in unit
         ! interval = <last second - delta, last second>
         SELECT CASE(TRIM(unit))
         CASE(TRIM(TIME_INC_SECONDS))
            ms = 999 
         CASE(TRIM(TIME_INC_MINUTES))
            se = 59; ms = 999 
         CASE(TRIM(TIME_INC_HOURS))
            mn = 59; se = 59; ms = 999
         CASE(TRIM(TIME_INC_DAYS))
            hr = 23; mn = 59; se = 59; ms = 999
         CASE(TRIM(TIME_INC_MONTHS))
            hr = 23; mn = 59; se = 59; ms = 999
            SELECT CASE(CAL_TYPE)
            CASE (CAL_JULIAN)
               dy = JulianMonthLength(yr,mo)  ! reset to the end of the month
            CASE (CAL_360D)
               dy = 30                        ! reset to the end of the month
            END SELECT
         CASE(TRIM(TIME_INC_YEARS))
            mo = 12; dy = 31; hr = 23; mn = 59; se = 59; ms = 999
         END SELECT

#ifdef LDEBUG
         write (*,*) 'adjust_date 10 ', yr, mo, dy, hr, mn, se, ms
#endif
         CALL timer_set_date(status, my_date, yr, mo, dy, hr, mn, se, ms)
#ifdef LDEBUG
         CALL print_date_components (my_date, status)
#endif             
      END SELECT

    END SUBROUTINE adjust_date

  END SUBROUTINE get_next_trigger

  ! --------------------------------------------------------------------------

  SUBROUTINE event_next_date (event, date)

    TYPE (time_event), INTENT (in)  :: event
    TYPE (time_days),  INTENT (out) :: date

    date = event %next_trigger

  END SUBROUTINE event_next_date

  ! ------------------------------------------------------------------------

  SUBROUTINE event_previous_date (event, date)

    TYPE (time_event), INTENT (in)  :: event
    TYPE (time_days),  INTENT (out) :: date

    date = event %previous_trigger

  END SUBROUTINE event_previous_date

  ! --------------------------------------------------------------------------

  SUBROUTINE event_current_date (event, date) 

    TYPE (time_event), INTENT (in)  :: event
    TYPE (time_days),  INTENT (out) :: date

    date = event %current_trigger

  END SUBROUTINE event_current_date

  ! --------------------------------------------------------------------------

  SUBROUTINE event_cycle_date (event, date) 

    TYPE (time_event), INTENT (in)  :: event
    TYPE (time_days),  INTENT (out) :: date

    date = event %cycle_date

  END SUBROUTINE event_cycle_date

  ! --------------------------------------------------------------------------

  SUBROUTINE event_initial_date (event, date) 

    TYPE (time_event), INTENT (in)  :: event
    TYPE (time_days),  INTENT (out) :: date

    date = event %initial_date

  END SUBROUTINE event_initial_date

  ! --------------------------------------------------------------------------

  LOGICAL FUNCTION event_is_active (event) 

    TYPE (time_event), INTENT (in)  :: event

    event_is_active = event %active

  END FUNCTION event_is_active

  ! --------------------------------------------------------------------------

  LOGICAL FUNCTION event_is_init (event) 

    TYPE (time_event), INTENT (in)  :: event

    event_is_init = event %init

  END FUNCTION event_is_init

  ! --------------------------------------------------------------------------

  REAL(dp) FUNCTION event_delta (event) 

    TYPE (time_event), INTENT (in)  :: event

    event_delta = event %delta

  END FUNCTION event_delta

  ! --------------------------------------------------------------------------

  REAL(dp) FUNCTION event_half_delta (event) 

    TYPE (time_event), INTENT (in)  :: event

    event_half_delta = event %half_delta

  END FUNCTION event_half_delta

  ! --------------------------------------------------------------------------

  CHARACTER(LEN=STRLEN_SHORT) FUNCTION event_adjust (event) 

    TYPE (time_event), INTENT (in)  :: event

    event_adjust = event %adjust

  END FUNCTION event_adjust

  ! --------------------------------------------------------------------------

  CHARACTER(LEN=STRLEN_MEDIUM+4) FUNCTION event_name (event) 

    TYPE (time_event), INTENT (in)  :: event

    IF (TRIM(event_tagstr(event)) == '') THEN
       event_name = TRIM(event %label)
    ELSE
       event_name = TRIM(event %label)//'_T'//TRIM(event_tagstr(event))
    END IF

  END FUNCTION event_name

  ! --------------------------------------------------------------------------

  INTEGER FUNCTION event_tag (event) 

    TYPE (time_event), INTENT (in)  :: event

    event_tag = event %tag

  END FUNCTION event_tag

  ! --------------------------------------------------------------------------

  CHARACTER(LEN=2) FUNCTION event_tagstr (event) 

    TYPE (time_event), INTENT (in)  :: event
 
    IF (event %tag >= 0 ) THEN
       CALL INT2STR(event_tagstr, event %tag,'0')
    ELSE
       event_tagstr = '  '
    END IF

  END FUNCTION event_tagstr

  ! --------------------------------------------------------------------------

  CHARACTER(LEN=STRLEN_MEDIUM) FUNCTION event_label (event) 

    TYPE (time_event), INTENT (in)  :: event

    event_label = event %label

  END FUNCTION event_label

  ! --------------------------------------------------------------------------

  CHARACTER(LEN=STRLEN_SHORT) FUNCTION event_unit (event) 

    TYPE (time_event), INTENT (in)  :: event

    event_unit = event %unit

  END FUNCTION event_unit

  ! ---------------------------------------------------------------------------

  INTEGER FUNCTION event_count (event) 

    TYPE (time_event), INTENT (in)  :: event

    event_count = event %count

  END FUNCTION event_count

  ! --------------------------------------------------------------------------

  REAL(dp) FUNCTION event_offset (event)

    TYPE (time_event), INTENT (in)  :: event

    event_offset = event %offset

  END FUNCTION event_offset

  ! --------------------------------------------------------------------------

  SUBROUTINE convert_steps2unit (steps, dt, unit, count, status, message) 

    IMPLICIT NONE
    INTRINSIC :: NINT, REAL

    ! convert steps into normal event units
    INTEGER,          INTENT(in)  :: steps  ! interval in steps
    REAL(dp),         INTENT(in)  :: dt     ! length of on step in seconds
    CHARACTER(len=*), INTENT(out) :: unit   ! new unit
    INTEGER,          INTENT(out) :: count  ! new counter
    INTEGER,          INTENT(OUT) :: status   ! error status
    CHARACTER(len=STRLEN_ULONG), INTENT(OUT) ::message  ! error message


    CHARACTER(LEN=*), PARAMETER:: substr = 'convert_steps2unit' 
    INTEGER :: isdt, imdt, ihdt, iddt, isteps
    REAL(dp):: rdt

    INTEGER, PARAMETER :: imax = 2000000000

    status = 0
    message = ''

    isteps = steps
    isdt   = NINT(dt)
    rdt    = REAL(isdt,dp)

    IF (ABS(rdt-dt) > TINY(1._dp)) THEN
       IF (isteps > imax ) THEN
          WRITE(message,*) 'step number too large: ',&
               ' maximum is ',imax/NINT(dt),' < ',isteps
          status = 3400
          RETURN
       END IF

       ! in this case only milliseconds are possible
       unit = TIME_INC_MILSECS
       count= NINT(dt*real(isteps,dp)/ CONVERTMS(I_INC_MILSECS))
    ELSE
       IF (isteps > imax/NINT(dt) ) THEN
          WRITE(message,*) 'step number too large: ',&
               ' maximum is ',imax/NINT(dt),' < ',isteps
          status = 3400
          RETURN
       END IF

       isdt = NINT(dt*real(isteps,dp))   
       imdt = isdt / CONVERT(I_INC_MINUTES)
       ihdt = isdt / CONVERT(I_INC_HOURS)
       iddt = isdt / CONVERT(I_INC_DAYS)

       count = 0

       IF (iddt * CONVERT(I_INC_DAYS) == isdt) THEN 
          ! multiple of seconds per day
          unit  = TIME_INC_DAYS;    count = iddt
       ELSE IF (ihdt * CONVERT(I_INC_HOURS) == isdt) THEN 
          ! seconds per hour
          unit  = TIME_INC_HOURS;   count = ihdt
       ELSE IF (imdt * CONVERT(I_INC_MINUTES) == isdt) THEN
          ! seconds per minute
          unit  = TIME_INC_MINUTES; count = imdt
       ELSE
          unit = TIME_INC_SECONDS; count = isdt
       END IF

    END IF

  END SUBROUTINE convert_steps2unit

  ! -----------------------------------------------------------------------

  SUBROUTINE convert_unit2seconds (status, seconds, TIMER) 

    USE messy_main_timer, ONLY: delta_time

    IMPLICIT NONE

    INTEGER, INTENT(OUT)            :: status
    REAL(dp), INTENT(OUT)           :: seconds
    TYPE(IO_TIME_EVENT), INTENT(IN) :: TIMER

    ! LOCAL
    !      INTEGER :: factor
    REAL(dp) :: factor

    SELECT CASE(TRIM(TIMER%unit))
    CASE(TRIM(TIME_INC_STEPS))
       !         factor = INT(delta_time)
       factor = delta_time
    CASE(TRIM(TIME_INC_MILSECS))
       factor = 0.001_dp
    CASE(TRIM(TIME_INC_SECONDS))
       !         factor = 1
       factor = 1._dp
    CASE(TRIM(TIME_INC_MINUTES))
       !         factor = 60
       factor = 60._dp
    CASE(TRIM(TIME_INC_HOURS))
       !         factor = 3600
       factor = 3600._dp
    CASE(TRIM(TIME_INC_DAYS))
       !         factor = 86400
       factor = 86400._dp
    CASE(TRIM(TIME_INC_MONTHS))
       status = 3450
       RETURN
    CASE(TRIM(TIME_INC_YEARS))
       status = 3451
       RETURN
    END SELECT

    !      seconds = factor * TIMER%counter
    seconds = factor * REAL(TIMER%counter,dp)

    status  = 0

  END SUBROUTINE convert_unit2seconds

  ! -----------------------------------------------------------------------

! ****************************************************************************
END MODULE messy_main_timer_event
! ****************************************************************************
