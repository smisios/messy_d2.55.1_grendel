! ************************************************************************
MODULE messy_main_timer
! ************************************************************************

  ! MESSy
  USE messy_main_constants_mem, ONLY: dp, OneDay, STRLEN_ULONG

  IMPLICIT NONE
  PRIVATE
  SAVE

  CHARACTER(LEN=*), PUBLIC, PARAMETER :: modstr='timer'
  CHARACTER(LEN=*), PUBLIC, PARAMETER :: modver='0.1'

  ! CALENDER TYPE (per default Julian)
  INTEGER, PARAMETER, PUBLIC :: CAL_JULIAN = 0
  INTEGER, PARAMETER, PUBLIC :: CAL_360D   = 1
  !
  INTEGER, PUBLIC :: CAL_TYPE = 0

  CHARACTER(LEN=8), PARAMETER, PUBLIC :: CALENDER(0:1) = &
       (/ 'standard', '360_day '/)

  CHARACTER(len=3), PARAMETER, PUBLIC :: CMONTHS(12) = &
       (/ 'Jan','Feb','Mar','Apr','May','Jun',&
          'Jul','Aug','Sep','Oct','Nov','Dec' /)

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! DATE MANAGEMENT
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  LOGICAL, PUBLIC :: LDEBUG2 = .FALSE.
!  LOGICAL, PUBLIC :: LDEBUG2 = .TRUE.
  LOGICAL, PUBLIC :: LDEBUG  = .FALSE.
!  LOGICAL, PUBLIC :: LDEBUG  = .TRUE.

  TYPE, PUBLIC :: time_days 
     !
     ! relative calendar date and time format
     !
     ! time_days [structure]
     !   day    [integer]  (day in calendar, -2147483648 ..... 2147483647
     !                      approx. +/-5.8 Mio. years)
     !   second [integer]  (seconds of day, 0,...,86399)
     !
     ! op_bk_20161013+
     !   milliseconds [integer]
     ! op_bk_20161013-
     !PRIVATE
     LOGICAL :: init   = .FALSE.
     INTEGER :: day    = 0
     INTEGER :: second = 0
     INTEGER :: millisecond = 0
  END TYPE time_days

  ! number of system clocks
  INTEGER, PUBLIC           :: nclock

  ! DEFINE SPECIFIC INFORMATION about SIMULATION STATE
  INTEGER,  PUBLIC, POINTER ::    INIT_STEP  => NULL()
  REAL(dp), PUBLIC, POINTER ::    delta_time => NULL() ! BM model time step
  REAL(dp), PUBLIC, POINTER :: time_step_len => NULL() ! leap frog BM 
  !                                                    ! model time step
  LOGICAL,PUBLIC    :: lresume = .FALSE.  ! .TRUE. during rerun step
  LOGICAL,PUBLIC    :: lbreak  = .FALSE.  ! .TRUE. at end of one time segment
  LOGICAL,PUBLIC    :: lstop   = .FALSE.  ! .TRUE. during the last time step
  LOGICAL,PUBLIC,POINTER :: lstart       => NULL() ! T for the first time step
  LOGICAL,PUBLIC,POINTER :: lfirst_cycle => NULL() ! T during 1st ts/rerun cycle
  LOGICAL,PUBLIC,POINTER :: l_rerun  => NULL()

  LOGICAL,PUBLIC    :: L_TRIGGER_RESTART = .FALSE.
  ! The following two logicals are needed to control the restart behaviour
  ! of coupled (MMD) simulations.
  ! regular/cycle induced break i.e. lbreak without L_TRIGGER_RESTART
  LOGICAL,PUBLIC    :: lcycbreak = .FALSE.  ! .TRUE. at end of one time segment
  ! LOGICAL INIDCATING IF TIMING IS DRIVEN BY OTHER MODEL (SERVER)
  LOGICAL, PUBLIC   :: lforcedtime = .FALSE.

  ! no_cycles: DUMMY needed for adjustment to ECHAM5 rerun structure
  INTEGER, PUBLIC :: NO_CYCLES = 9999

  INTEGER, PUBLIC :: NO_DAYS  = -9999
  INTEGER, PUBLIC :: NO_STEPS = -9999
  LOGICAL, PUBLIC :: LABORT = .FALSE.

  INTEGER, PUBLIC, POINTER :: YEAR              => NULL()
  INTEGER, PUBLIC, POINTER :: MONTH             => NULL()
  INTEGER, PUBLIC, POINTER :: DAY               => NULL()
  INTEGER, PUBLIC, POINTER :: HOUR              => NULL()
  INTEGER, PUBLIC, POINTER :: MINUTE            => NULL()
  INTEGER, PUBLIC, POINTER :: SECOND            => NULL()
  INTEGER, PUBLIC, POINTER :: MILLISECOND       => NULL()

  INTEGER, PUBLIC, POINTER :: YEAR_START        => NULL()
  INTEGER, PUBLIC, POINTER :: MONTH_START       => NULL()
  INTEGER, PUBLIC, POINTER :: DAY_START         => NULL()
  INTEGER, PUBLIC, POINTER :: HOUR_START        => NULL()
  INTEGER, PUBLIC, POINTER :: MINUTE_START      => NULL()
  INTEGER, PUBLIC, POINTER :: SECOND_START      => NULL()
  INTEGER, PUBLIC, POINTER :: MILLISECOND_START => NULL()

  INTEGER, PUBLIC, POINTER :: YEAR_NEXT         => NULL()
  INTEGER, PUBLIC, POINTER :: MONTH_NEXT        => NULL()
  INTEGER, PUBLIC, POINTER :: DAY_NEXT          => NULL()
  INTEGER, PUBLIC, POINTER :: HOUR_NEXT         => NULL()
  INTEGER, PUBLIC, POINTER :: MINUTE_NEXT       => NULL()
  INTEGER, PUBLIC, POINTER :: SECOND_NEXT       => NULL()
  INTEGER, PUBLIC, POINTER :: MILLISECOND_NEXT  => NULL()

  INTEGER, PUBLIC, POINTER :: current_time_step => NULL()
  
  ! day of year [day] 
  INTEGER,  PUBLIC, POINTER :: DAYOFYEAR        => NULL()     

  REAL(dp),         PUBLIC  :: JULIAN_DATE_START 

  TYPE(time_days), PUBLIC   :: start_date         ! start date
  TYPE(time_days), PUBLIC   :: stop_date          ! stop date
  TYPE(time_days), PUBLIC   :: resume_date        ! rerun date
  TYPE(time_days), PUBLIC   :: rerun_stop_date  

  ! DEFINE SPECIFIC DATE INFORMATION
  
  ! date at (time - delta_time)
  TYPE(time_days),PUBLIC, POINTER  ::  previous_date   => NULL()   
  ! date at (time)
  TYPE(time_days),PUBLIC, POINTER  ::   current_date   => NULL()   
  ! date at (time + delta_time)
  TYPE(time_days),PUBLIC, POINTER  ::      next_date   => NULL()  

  TYPE t_clock_settings
     INTEGER,  POINTER ::        INIT_STEP  => NULL()
     REAL(dp), POINTER ::        delta_time => NULL() 
     REAL(dp), POINTER ::     time_step_len => NULL() 
     INTEGER,  POINTER ::        YEAR       => NULL()
     INTEGER,  POINTER ::       MONTH       => NULL()
     INTEGER,  POINTER ::         DAY       => NULL()
     INTEGER,  POINTER ::        HOUR       => NULL()
     INTEGER,  POINTER ::      MINUTE       => NULL()
     INTEGER,  POINTER ::      SECOND       => NULL()
     INTEGER,  POINTER :: MILLISECOND       => NULL()

     INTEGER,  POINTER ::        YEAR_START => NULL()
     INTEGER,  POINTER ::       MONTH_START => NULL()
     INTEGER,  POINTER ::         DAY_START => NULL()
     INTEGER,  POINTER ::        HOUR_START => NULL()
     INTEGER,  POINTER ::      MINUTE_START => NULL()
     INTEGER,  POINTER ::      SECOND_START => NULL()
     INTEGER,  POINTER :: MILLISECOND_START => NULL()

     INTEGER,  POINTER ::        YEAR_NEXT  => NULL()
     INTEGER,  POINTER ::       MONTH_NEXT  => NULL()
     INTEGER,  POINTER ::         DAY_NEXT  => NULL()
     INTEGER,  POINTER ::        HOUR_NEXT  => NULL()
     INTEGER,  POINTER ::      MINUTE_NEXT  => NULL()
     INTEGER,  POINTER ::      SECOND_NEXT  => NULL()
     INTEGER,  POINTER :: MILLISECOND_NEXT  => NULL()

     INTEGER, POINTER :: current_time_step  => NULL()

     ! day of year [day] 
     INTEGER, POINTER :: DAYOFYEAR          => NULL()

     ! DEFINE SPECIFIC DATE INFORMATION
     ! date at (time - delta_time)
     TYPE(time_days), POINTER  ::  previous_date  => NULL()
     ! date at (time)
     TYPE(time_days), POINTER  ::   current_date  => NULL()  
     ! date at (time + delta_time)
     TYPE(time_days), POINTER  ::      next_date  => NULL()   

     LOGICAL, POINTER :: lstart       => NULL()   
     LOGICAL, POINTER :: lfirst_cycle => NULL() 
     LOGICAL, POINTER :: l_rerun      => NULL()     

  END type t_clock_settings

  TYPE(t_clock_settings), DIMENSION(:), POINTER, PUBLIC :: clock => NULL()

  INTEGER, PUBLIC :: current_clock = -99

  PUBLIC :: timer_alloc_clock
  PUBLIC :: timer_dealloc_clock
  PUBLIC :: main_timer_set_clock
  PUBLIC :: t_clock_settings

  PUBLIC :: date_set
  PUBLIC :: date_get
  !PRIVATE :: date_get_components
  !PRIVATE :: date_set_components
  INTERFACE add_date
     MODULE PROCEDURE add_date_int
     MODULE PROCEDURE add_date_real
  END INTERFACE
  PUBLIC :: add_date
  PUBLIC :: copy_date
  PUBLIC :: if_less
  PUBLIC :: if_equal
  PUBLIC :: is_init
  PUBLIC :: print_date
  PUBLIC :: print_date_components

  ! -----------------------------------------------------------------
  ! HELPER ROUTINES FOR DATE CONVERSIONS / TIME DISTANCE CALCULATIONS
  ! -----------------------------------------------------------------
  !
  PUBLIC :: MonthLength         ! function to calculate the length of a year
  !                             ! in days
  PUBLIC :: JulianMonthLength   ! function to calculate the length of a 
  !                             ! julian month in days
  PUBLIC :: YearLength          ! function to calculate the length of a year
  !                             ! in days
  PUBLIC :: JulianYearLength    ! function to calculate the length of a 
  !                             ! julian year in days
  PUBLIC :: YearDay             ! function to calculate the current number of a 
  !                             ! day in the current year 
  PUBLIC :: julian_day          ! function to calculate the Julian day
  !
  INTERFACE time_span_s
     MODULE PROCEDURE time_span_s_int
     MODULE PROCEDURE time_span_s_real
  END INTERFACE
  PUBLIC :: time_span_s         ! time difference [s] between two greg.dates
  INTERFACE time_span_d
     MODULE PROCEDURE time_span_d_s
     MODULE PROCEDURE time_span_d_ms
  END INTERFACE
  PUBLIC :: time_span_d         ! time difference [d] between two greg.dates
  !
  PUBLIC :: gregor2julian       ! convert gregorian date + time to julian date
  PUBLIC :: julian2gregor       ! convert julian date to gregorian date + time
  PUBLIC :: utc2lt              ! calculate local time from UTC time and
                                ! longitude in degrees
  PUBLIC :: eval_time_str       ! evaluate netcdf time string
  PUBLIC :: calc_sza            ! calculate solar zenith angle based on time
  !                               and longitude
  !
  ! --------------------------------------------------------------
  ! INTERFACE ROUTINES FOR TIMER EXTERNAL USE
  ! --------------------------------------------------------------
  ! +++ TRANSPORT NAMELIST INFORMATION IN BOTH DIRECTIONS
  INTERFACE timer_set_date
     MODULE PROCEDURE timer_set_date_str
     MODULE PROCEDURE timer_set_date_myd
     MODULE PROCEDURE timer_set_date_str_ds_int
     MODULE PROCEDURE timer_set_date_str_ds_real
     MODULE PROCEDURE timer_set_date_myd_ds_int
     MODULE PROCEDURE timer_set_date_myd_ds_real
  END INTERFACE
  PUBLIC :: timer_set_date

  INTERFACE timer_get_date
     MODULE PROCEDURE timer_get_date_str
     MODULE PROCEDURE timer_get_date_myd
  END INTERFACE
  PUBLIC :: timer_get_date

  PUBLIC :: timer_set_calendar
  PUBLIC :: timer_get_calendar
  PUBLIC :: timer_set_delta_time
  PUBLIC :: timer_get_delta_time
  PUBLIC :: timer_set_lstart
  PUBLIC :: timer_set_lfirst_cycle
  PUBLIC :: timer_set_lresume
  PUBLIC :: timer_get_lresume
  PUBLIC :: timer_get_l_rerun
  PUBLIC :: timer_set_no_cycles
  PUBLIC :: timer_get_no_cycles
  PUBLIC :: timer_set_labort
  PUBLIC :: timer_get_labort
  PUBLIC :: timer_get_init_step
  PUBLIC :: timer_get_next_date

  ! --- TRANSPORT NAMELIST INFORMATION IN BOTH DIRECTIONS

  INTERFACE timer_add_date
     MODULE PROCEDURE timer_add_date_int
     MODULE PROCEDURE timer_add_date_real
     MODULE PROCEDURE timer_add_seconds_to_date_real
  END INTERFACE
  PUBLIC :: timer_add_date
  PUBLIC :: timer_set_time_step_len
  PUBLIC :: timer_get_time_step_len

  INTERFACE date_set
     MODULE PROCEDURE date_set_int
     MODULE PROCEDURE date_set_real
  END INTERFACE

CONTAINS

  !--------------------------------------------------------------------------
  SUBROUTINE timer_alloc_clock(numclock)

    IMPLICIT NONE
    
    INTEGER, INTENT(IN), OPTIONAL :: numclock

    INTEGER :: nc, ic

    IF (PRESENT(numclock)) THEN
       nc = numclock
    ELSE
       nc = 1
    END IF
    
    nclock = nc

    ALLOCATE(clock(nc)) 

    DO ic = 1, nc

       ALLOCATE(clock(ic)%INIT_STEP);        clock(ic)%INIT_STEP     = 0
       ALLOCATE(clock(ic)%delta_time);       clock(ic)%delta_time    = -999._dp
       ALLOCATE(clock(ic)%time_step_len);    clock(ic)%time_step_len = 0._dp

       ALLOCATE(clock(ic)%YEAR             ); clock(ic)%YEAR        = 0
       ALLOCATE(clock(ic)%MONTH            ); clock(ic)%MONTH       = 0
       ALLOCATE(clock(ic)%DAY              ); clock(ic)%DAY         = 0
       ALLOCATE(clock(ic)%HOUR             ); clock(ic)%HOUR        = 0
       ALLOCATE(clock(ic)%MINUTE           ); clock(ic)%MINUTE      = 0
       ALLOCATE(clock(ic)%SECOND           ); clock(ic)%SECOND      = 0
       ALLOCATE(clock(ic)%MILLISECOND      ); clock(ic)%MILLISECOND = 0

       ALLOCATE(clock(ic)%YEAR_START       );  clock(ic)%YEAR_START        = 0
       ALLOCATE(clock(ic)%MONTH_START      );  clock(ic)%MONTH_START       = 0
       ALLOCATE(clock(ic)%DAY_START        );  clock(ic)%DAY_START         = 0
       ALLOCATE(clock(ic)%HOUR_START       );  clock(ic)%HOUR_START        = 0
       ALLOCATE(clock(ic)%MINUTE_START     );  clock(ic)%MINUTE_START      = 0
       ALLOCATE(clock(ic)%SECOND_START     );  clock(ic)%SECOND_START      = 0
       ALLOCATE(clock(ic)%MILLISECOND_START);  clock(ic)%MILLISECOND_START = 0

       ALLOCATE(clock(ic)%YEAR_NEXT        );  clock(ic)%YEAR_NEXT         = 0
       ALLOCATE(clock(ic)%MONTH_NEXT       );  clock(ic)%MONTH_NEXT        = 0
       ALLOCATE(clock(ic)%DAY_NEXT         );  clock(ic)%DAY_NEXT          = 0
       ALLOCATE(clock(ic)%HOUR_NEXT        );  clock(ic)%HOUR_NEXT         = 0
       ALLOCATE(clock(ic)%MINUTE_NEXT      );  clock(ic)%MINUTE_NEXT       = 0
       ALLOCATE(clock(ic)%SECOND_NEXT      );  clock(ic)%SECOND_NEXT       = 0
       ALLOCATE(clock(ic)%MILLISECOND_NEXT );  clock(ic)%MILLISECOND_NEXT  = 0

       ALLOCATE(clock(ic)%current_time_step); clock(ic)%current_time_step  = 0
       ALLOCATE(clock(ic)%DAYOFYEAR        )

       ALLOCATE(clock(ic)%previous_date    )
       ALLOCATE(clock(ic)%current_date     )
       ALLOCATE(clock(ic)%next_date        )

       ALLOCATE(clock(ic)%lstart           ) ; clock(ic)%lstart  = .TRUE.
       ALLOCATE(clock(ic)%lfirst_cycle     ) ; clock(ic)%lfirst_cycle =.TRUE.
       !ALLOCATE(clock(ic)%lresume          ) ; clock(ic)%lresume = .FALSE.
       !ALLOCATE(clock(ic)%lbreak           ) ; clock(ic)%lbreak  = .FALSE.
       !ALLOCATE(clock(ic)%lstop            ) ; clock(ic)%lstop  = .FALSE.

       ALLOCATE(clock(ic)%l_rerun          ) ; clock(ic)%l_rerun = .FALSE.
    END DO

  END SUBROUTINE timer_alloc_clock
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  SUBROUTINE timer_dealloc_clock

    IMPLICIT NONE
    
    INTEGER :: ic

    DO ic = 1, nclock

       DEALLOCATE(clock(ic)%INIT_STEP        ); NULLIFY(clock(ic)%INIT_STEP)   
       DEALLOCATE(clock(ic)%delta_time       ); NULLIFY(clock(ic)%delta_time)
       DEALLOCATE(clock(ic)%time_step_len    ); NULLIFY(clock(ic)%time_step_len)

       DEALLOCATE(clock(ic)%YEAR             )
       DEALLOCATE(clock(ic)%MONTH            )
       DEALLOCATE(clock(ic)%DAY              )
       DEALLOCATE(clock(ic)%HOUR             )
       DEALLOCATE(clock(ic)%MINUTE           )
       DEALLOCATE(clock(ic)%SECOND           )
       DEALLOCATE(clock(ic)%MILLISECOND      )

       DEALLOCATE(clock(ic)%YEAR_START       )
       DEALLOCATE(clock(ic)%MONTH_START      )
       DEALLOCATE(clock(ic)%DAY_START        )
       DEALLOCATE(clock(ic)%HOUR_START       )
       DEALLOCATE(clock(ic)%MINUTE_START     )
       DEALLOCATE(clock(ic)%SECOND_START     )
       DEALLOCATE(clock(ic)%MILLISECOND_START)

       DEALLOCATE(clock(ic)%YEAR_NEXT        )
       DEALLOCATE(clock(ic)%MONTH_NEXT       )
       DEALLOCATE(clock(ic)%DAY_NEXT         )
       DEALLOCATE(clock(ic)%HOUR_NEXT        )
       DEALLOCATE(clock(ic)%MINUTE_NEXT      )
       DEALLOCATE(clock(ic)%SECOND_NEXT      )
       DEALLOCATE(clock(ic)%MILLISECOND_NEXT )

       DEALLOCATE(clock(ic)%current_time_step)
       DEALLOCATE(clock(ic)%DAYOFYEAR        )

       DEALLOCATE(clock(ic)%previous_date    )
       DEALLOCATE(clock(ic)%current_date     )
       DEALLOCATE(clock(ic)%next_date        )

       DEALLOCATE(clock(ic)%lstart           )
       DEALLOCATE(clock(ic)%lfirst_cycle     )
       DEALLOCATE(clock(ic)%l_rerun          )

       NULLIFY(clock(ic)%YEAR             )
       NULLIFY(clock(ic)%MONTH            )
       NULLIFY(clock(ic)%DAY              )
       NULLIFY(clock(ic)%HOUR             )
       NULLIFY(clock(ic)%MINUTE           )
       NULLIFY(clock(ic)%SECOND           )
       NULLIFY(clock(ic)%MILLISECOND      )

       NULLIFY(clock(ic)%YEAR_START       )
       NULLIFY(clock(ic)%MONTH_START      )
       NULLIFY(clock(ic)%DAY_START        )
       NULLIFY(clock(ic)%HOUR_START       )
       NULLIFY(clock(ic)%MINUTE_START     )
       NULLIFY(clock(ic)%SECOND_START     )
       NULLIFY(clock(ic)%MILLISECOND_START)

       NULLIFY(clock(ic)%YEAR_NEXT        )
       NULLIFY(clock(ic)%MONTH_NEXT       )
       NULLIFY(clock(ic)%DAY_NEXT         )
       NULLIFY(clock(ic)%HOUR_NEXT        )
       NULLIFY(clock(ic)%MINUTE_NEXT      )
       NULLIFY(clock(ic)%SECOND_NEXT      )
       NULLIFY(clock(ic)%MILLISECOND_NEXT )

       NULLIFY(clock(ic)%current_time_step)
       NULLIFY(clock(ic)%DAYOFYEAR        )

       NULLIFY(clock(ic)%previous_date    )
       NULLIFY(clock(ic)%current_date     )
       NULLIFY(clock(ic)%next_date        )

       NULLIFY(clock(ic)%lstart           )
       NULLIFY(clock(ic)%lfirst_cycle     )
       !NULLIFY(clock(ic)%lresume          )
       !NULLIFY(clock(ic)%lbreak           )
       !NULLIFY(clock(ic)%lstop            )
       NULLIFY(clock(ic)%l_rerun          )
    END DO
    DEALLOCATE(clock) 
    NULLIFY(clock)

    NULLIFY(INIT_STEP        )  
    NULLIFY(delta_time       )
    NULLIFY(time_step_len    )
    
    NULLIFY(YEAR             )
    NULLIFY(MONTH            )
    NULLIFY(DAY              )
    NULLIFY(HOUR             )
    NULLIFY(MINUTE           )
    NULLIFY(SECOND           )
    NULLIFY(MILLISECOND      )
    
    NULLIFY(YEAR_START       )
    NULLIFY(MONTH_START      )
    NULLIFY(DAY_START        )
    NULLIFY(HOUR_START       )
    NULLIFY(MINUTE_START     )
    NULLIFY(SECOND_START     )
    NULLIFY(MILLISECOND_START)
    
    NULLIFY(YEAR_NEXT        )
    NULLIFY(MONTH_NEXT       )
    NULLIFY(DAY_NEXT         )
    NULLIFY(HOUR_NEXT        )
    NULLIFY(MINUTE_NEXT      )
    NULLIFY(SECOND_NEXT      )
    NULLIFY(MILLISECOND_NEXT )
    
    NULLIFY(current_time_step)
    NULLIFY(DAYOFYEAR        )
    
    NULLIFY(previous_date    )
    NULLIFY(current_date     )
    NULLIFY(next_date        )
    
    NULLIFY(lstart           )
    NULLIFY(lfirst_cycle     )
    !NULLIFY(lresume          )
    !NULLIFY(lbreak           )
    !NULLIFY(lstop            )
    NULLIFY(l_rerun          )
    
  END SUBROUTINE timer_dealloc_clock
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  SUBROUTINE date_set_int(day, second, time)
    
    INTEGER, INTENT(IN) :: day
    INTEGER, INTENT(IN) :: second
    TYPE(time_days), INTENT(OUT) :: time

    time %day    = day
    time %second = second
    time%init    =.TRUE.

  END SUBROUTINE date_set_int
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  SUBROUTINE date_set_real(day, second, time)
    
    INTEGER, INTENT(IN) :: day
    REAL(DP), INTENT(IN):: second
    TYPE(time_days), INTENT(OUT) :: time

    time %day         = day
    time %second      = FLOOR(second)
    time %millisecond = INT((second - AINT(second,dp))  * 1000._dp)
    time%init         =.TRUE.

  END SUBROUTINE date_set_real
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  SUBROUTINE date_set_components(nyr,nmo,ndy,nhr,nmin,nsec,nms,i_time)
    
    IMPLICIT NONE
    INTRINSIC :: AINT, INT

    INTEGER, INTENT(IN) :: nyr,nmo,ndy,nhr,nmin,nsec,nms
    TYPE(time_days), INTENT(OUT) :: i_time
    REAL(dp) :: julday, julfrac, zsec, julianday
    REAL(dp) :: zisec
    INTEGER  :: kday, imillisec, isecs
    
    SELECT CASE(CAL_TYPE)
    CASE (CAL_JULIAN)
       julianday = gregor2julian(nyr,nmo,ndy,nhr,nmin,nsec,nms)
       Julday = AINT(gregor2julian(nyr,nmo,ndy,nhr,nmin,nsec,nms),dp)
       Julfrac = julianday- Julday      

       IF (LDEBUG2) write (*,*) 'set', julday, julfrac, julianday, &
            ' : ',nyr,nmo,ndy,nhr,nmin,nsec,nms
       IF (julday  < 0.0_dp) THEN
          kday = INT(julday -0.00001_dp)
          IF (julfrac < -0.5_dp)       kday = kday - 1
       ELSE 
          kday = INT(julday+0.000001_dp)
          IF (julfrac > 0.0_dp) THEN
             IF (.NOT.(julfrac < 0.5_dp)) kday = kday + 1
          ELSE IF (julfrac < -0.5_dp) THEN
             kday = kday - 1
          ELSE IF (.NOT.(julfrac < 0.5_dp)) THEN
             kday = kday + 1
          END IF
       END IF
       
       IF  (julfrac < -0.5_dp) THEN
          zsec = julfrac + 1.5_dp
       ELSE IF (julfrac < 0.5_dp) THEN
          zsec = julfrac + 0.5_dp
       ELSE
          zsec = julfrac - 0.5_dp
       END IF

       zisec = AINT(zsec*Oneday,dp) 
       zsec  = zsec * OneDay        
       IF (LDEBUG2) write (*,*) 'set2', zsec, zisec, kday, julfrac
       imillisec = NINT((zsec - zisec + 0.00000001_dp)*1000._dp)

       IF (imillisec > 999) THEN
          isecs =  INT(REAL(imillisec,dp)/1000._dp)
          imillisec = imillisec - isecs * 1000
       ELSE
          isecs = 0
       END IF
       
       i_time %day    = kday
       i_time %second = INT(zisec)+isecs
       i_time %millisecond = imillisec
       IF (LDEBUG2) write (*,*) 'set3', zsec, zisec, imillisec, kday

    CASE (CAL_360D)

       i_time %second = nhr*3600 + nmo*60 + nsec
       i_time %millisecond = nms
       i_time %day    = 360*nyr + (nmo-1)*30 + (ndy-1)
       
    END SELECT
    
    i_time %init   = .TRUE.
    
  END SUBROUTINE date_set_components
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  SUBROUTINE date_get(time, day, second, millisecond, ierr)

    IMPLICIT NONE
    INTRINSIC :: PRESENT

    TYPE(time_days), INTENT(in)    :: time
    INTEGER, OPTIONAL, INTENT(OUT) :: day
    INTEGER, OPTIONAL, INTENT(OUT) :: second
    INTEGER, OPTIONAL, INTENT(OUT) :: millisecond
    INTEGER, OPTIONAL, INTENT(OUT) :: ierr

    IF (.NOT. time%init) THEN
       IF (PRESENT(ierr)) ierr = 3435
       RETURN
    END IF

    IF (PRESENT(day))      day  = time%day
    IF (PRESENT(second)) second = time%second
    IF (PRESENT(millisecond)) millisecond = time%millisecond

    IF (PRESENT(ierr))     ierr = 0

  END SUBROUTINE date_get
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  SUBROUTINE date_get_components(time, year, month, day, hour&
       , minute, second, millisecond, ierror)

    IMPLICIT NONE
    INTRINSIC :: AINT, INT, MOD, PRESENT, REAL

    TYPE (time_days)  ,INTENT(in)   :: time
    INTEGER           ,INTENT(out)  :: year, month, day, hour                   
    INTEGER ,OPTIONAL ,INTENT(out)  :: minute, second
    INTEGER, OPTIONAL, INTENT(OUT)  :: millisecond
    INTEGER, OPTIONAL, INTENT(out)  :: ierror

    REAL(dp)   :: juldate, julfrac, julday
    INTEGER    :: mn,se,rest, iday
    INTEGER    :: ms
    REAL(dp)   :: zsecs

    ! INIT
    year = 0
    month = 0
    day = 0
    hour = 0
    mn = 0
    se = 0
    ms = 0

    IF (.NOT.time%init) THEN
       IF (PRESENT(ierror)) ierror = 3435
       RETURN
    ENDIF

    SELECT CASE(CAL_TYPE)
    CASE (CAL_JULIAN)
      ! remapping of julian day (adjustment at 12UTC) 
      ! and day with adjustment at 00UTC
      !
      zsecs = (REAL(time %second,dp) + &
           REAL(time%millisecond,DP)/1000._dp)/Oneday

      iday  = time %day
      julfrac = zsecs - 0.5_dp
      
      IF ((time %day < 0).AND. (zsecs > 0.5_dp)) THEN
        iday     = iday + 1
        julfrac = zsecs - 1.5_dp
      ELSE IF ((time%day > 0) .AND. (zsecs < 0.5_dp)) THEN
        iday     = iday - 1
        julfrac  = zsecs + 0.5_dp
      END IF
      IF (iday < 0) THEN
        julday = AINT(REAL(iday,dp)-0.00001_dp)
      ELSE
        julday = AINT(REAL(iday,dp)+0.00001_dp)
      END IF

      IF (LDEBUG2) write (*,*) 'get1', time%day, time%second &
           , time%millisecond, julday, julfrac
      juldate = julday + julfrac
      IF (LDEBUG2) write (*,*) 'get2',juldate
      CALL julian2gregor(juldate,year, month, day, hour, mn, se, ms)
      IF (LDEBUG2) write (*,*) 'get', year, month, day, hour, mn, se, ms

   CASE (CAL_360D)
       IF (time%day < 0) THEN
          year = INT(time%day/360-1)
          rest  = MOD(time%day,360)
          month = 12-rest/30
          day = 30-(MOD(rest,30)+1)
       ELSE
          year  = time%day/360
          rest  = MOD(time%day, 360)
          month = rest/30+1
          day   = MOD(rest,30)+1
       ENDIF

       hour = time%second/3600 + time%millisecond/3600000
       mn = MOD(time%second,3600)/60 + MOD(time%millisecond,3600000)/60000
       se = MOD(time%second,60) + MOD(time%millisecond,60000)
       ms = MOD(time%millisecond,1000)
      
    END SELECT
   
    IF (PRESENT(minute)) minute = mn
    IF (PRESENT(second)) second = se
    IF (PRESENT(millisecond)) millisecond = ms

    IF (PRESENT(ierror)) ierror = 0

  END SUBROUTINE date_get_components
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  SUBROUTINE add_date_int (days, seconds, my_day, ierr)

    ! add a number of days and seconds to a given date and time
    ! adopted from mo_time_conversion (ECHAM5)

    IMPLICIT NONE
    INTRINSIC :: INT, REAL, PRESENT

    INTEGER,          INTENT(in)    :: days, seconds
    TYPE (time_days), INTENT(inout) :: my_day
    INTEGER, OPTIONAL, INTENT(out)  :: ierr

    !CHARACTER(LEN=*), PARAMETER :: substr='main_timer: add_date'
    INTEGER   :: idays, isecs

    IF (PRESENT(ierr)) ierr=0

    IF (.NOT. my_day%init) THEN
       IF (PRESENT(ierr)) ierr=3435
       RETURN
    ENDIF

    isecs = seconds + my_day%second
    IF (isecs < 0) THEN
       idays = INT((REAL(isecs,dp)-0.0001_dp)/OneDay)
    ELSE
       idays = INT((REAL(isecs,dp)+0.0001_dp)/OneDay)
    END IF
    isecs = isecs - idays*INT(OneDay)
    idays = my_day%day + days + idays

    IF (isecs < 0) THEN
       isecs = INT(OneDay) + isecs
       idays = idays - 1
    END IF

    my_day %day    = idays
    my_day %second = isecs

  END SUBROUTINE add_date_int
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  SUBROUTINE add_date_real (days, seconds, my_day, ierr)

    ! add a number of days and milliseconds to a given date and time
    ! adopted from mo_time_conversion (ECHAM5) + adaption to milliseconds

    IMPLICIT NONE
    INTRINSIC :: INT, REAL, PRESENT

    INTEGER,          INTENT(in)    :: days
    REAL(DP),         INTENT(IN)    :: seconds
    TYPE (time_days), INTENT(inout) :: my_day
    INTEGER, OPTIONAL, INTENT(out)  :: ierr

    !CHARACTER(LEN=*), PARAMETER :: substr='main_timer: add_date'
    INTEGER   :: idays, isecs, imillisecs

    IF (PRESENT(ierr)) ierr=0

    IF (.NOT. my_day%init) THEN
       IF (PRESENT(ierr)) ierr=3435
       RETURN
    ENDIF

    imillisecs = FLOOR((seconds - AINT(seconds,dp) + 0.00001_dp ) * 1000._dp)
    imillisecs = imillisecs + my_day%millisecond  

    IF (imillisecs < 0) THEN
       isecs = INT(REAL(imillisecs,dp)/1000._dp) - 1
    ELSE IF (imillisecs > 999) THEN
       isecs =  INT(REAL(imillisecs,dp)/1000._dp)
    ELSE
       isecs = 0
    END IF
    imillisecs = imillisecs - isecs * 1000

    isecs = INT(seconds) + my_day%second + isecs

    IF (isecs < 0) THEN
       idays = INT((REAL(isecs,dp)-0.0001_dp)/OneDay)
    ELSE
       idays = INT((REAL(isecs,dp)+0.0001_dp)/OneDay)
    END IF
    isecs = isecs - idays*INT(OneDay)
    idays = my_day%day + days + idays

    IF (isecs < 0) THEN
       isecs = INT(OneDay) + isecs
       idays = idays - 1
    END IF

    my_day %day    = idays
    my_day %second = isecs
    my_day %millisecond = imillisecs

  END SUBROUTINE add_date_real
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  SUBROUTINE copy_date(date1,date2, ierr)

    IMPLICIT NONE
    INTRINSIC :: PRESENT

    TYPE(time_days), INTENT(IN)  :: date1 
    TYPE(time_days), INTENT(OUT) :: date2
    INTEGER, OPTIONAL, INTENT(OUT) :: ierr

    IF (PRESENT(ierr)) ierr = 0

    IF (.NOT. date1%init) THEN
       IF (PRESENT(ierr)) ierr = 3435
       RETURN
    ENDIF

    date2%day   =date1%day
    date2%second=date1%second
    date2%millisecond=date1%millisecond
    date2%init  = .TRUE.

  END SUBROUTINE copy_date
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  SUBROUTINE print_date (day,ierr, mess)

    IMPLICIT NONE
    INTRINSIC :: PRESENT, TRIM

    ! print out the date/time information
    TYPE (time_days),           INTENT(in)  :: day
    INTEGER,                    INTENT(OUT) :: ierr
    CHARACTER(len=*), OPTIONAL, INTENT(out) :: mess

    ! LOCAL
    !CHARACTER(LEN=*), PARAMETER :: substr='main_timer print_date'
    CHARACTER(LEN=STRLEN_ULONG) :: message_text
    INTEGER                     :: iday, isec, imillisec

    CALL date_get(day,iday,isec,imillisec,ierr)
    IF (ierr /= 0) RETURN

    SELECT CASE(CAL_TYPE)
    CASE (CAL_JULIAN)      
       WRITE(message_text,'(a,i8,a,i8,a,i8)') &
            'modified Julian day (00 UT adjusted): ', &
            iday,' seconds: ', isec,' milliseconds: ', imillisec
    CASE (CAL_360D)
       WRITE(message_text,'(a,i8,a,i8,a,i8)') &
            '360 day year day (00 UT based): ', &
            iday,' seconds: ', isec, ' milliseconds: ', imillisec
    END SELECT

    IF (PRESENT(mess)) THEN
       mess = TRIM(message_text)
    ELSE
       WRITE(*,*) message_text
    END IF

    ierr = 0

  END SUBROUTINE print_date
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  SUBROUTINE print_date_components (day, ierr, mess, info)

    IMPLICIT NONE
    INTRINSIC :: PRESENT, TRIM

    TYPE (time_days),           INTENT(IN)   :: day
    INTEGER,                    INTENT(OUT)  :: ierr
    CHARACTER(len=STRLEN_ULONG), OPTIONAL, INTENT(OUT) :: mess
    CHARACTER(len=*),            OPTIONAL, INTENT(IN)  :: info

    ! LOCAL
    !CHARACTER(LEN=*), PARAMETER :: substr='main_timer print_date'
    CHARACTER(LEN=STRLEN_ULONG) :: message_text
    INTEGER                     :: yr,mo,dy,hr,mn,se, ms

    ierr = 0
    
    CALL date_get_components(day,yr,mo,dy,hr,mn,se,ms,ierr)
    IF (ierr /= 0) RETURN
    
    WRITE(message_text,'(i2,a,a3,a,i4,3(a,i2.2),a,i3.3)') &
         dy,'. ',CMONTHS(mo),' ',yr,' ',hr,':',mn,':',se,'.',ms

    IF (PRESENT(info)) THEN
       message_text = TRIM(info)//' '//TRIM(message_text)
    END IF

    IF (PRESENT(mess)) THEN
       mess = TRIM(message_text)
    ELSE
       write (*,*) TRIM(message_text)
    END IF
    
    ierr = 0

  END SUBROUTINE print_date_components
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  SUBROUTINE if_less (date1, date2, lless, ierr)

    IMPLICIT NONE
    INTRINSIC :: PRESENT

    TYPE (time_days), INTENT(IN)   :: date1
    TYPE (time_days), INTENT(IN)   :: date2
    LOGICAL, INTENT(OUT)           :: lless
    INTEGER, OPTIONAL, INTENT(OUT) :: ierr

    IF ((.NOT.date1%init) .OR. (.NOT.date2%init)) THEN
       IF (PRESENT(ierr)) ierr =3435
       RETURN
    ENDIF

    IF ( (date1%day == date2%day .AND. date1%second > date2%second) &
         .OR. (date1%day == date2%day .AND. date1%second == date2%second &
         .AND. date1%millisecond >= date2%millisecond) &
         .OR. (date1%day > date2%day) ) THEN
       lless = .FALSE.
    ELSE
       lless = .TRUE.
    END IF
    
    IF (PRESENT(ierr)) ierr = 0

  END SUBROUTINE if_less
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  SUBROUTINE if_equal (date1, date2, leq, ierr)

    IMPLICIT NONE
    INTRINSIC :: PRESENT

    TYPE (time_days), INTENT(IN)   :: date1
    TYPE (time_days), INTENT(IN)   :: date2
    LOGICAL, INTENT(OUT)           :: leq
    INTEGER, OPTIONAL, INTENT(OUT) :: ierr

    IF ((.NOT.date1%init) .OR. (.NOT.date2%init)) THEN
       IF (PRESENT(ierr)) ierr =3435
       RETURN
    ENDIF

    IF ( (date1%day == date2%day) .AND. (date1%second == date2%second) .AND. &
         (date1%millisecond == date2%millisecond)) THEN
       leq =.TRUE.
    ELSE
       leq = .FALSE.
    END IF
    
    IF (PRESENT(ierr)) ierr = 0

  END SUBROUTINE if_equal
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  SUBROUTINE is_init(date,linit)

    TYPE(time_days), INTENT(IN) :: date
    LOGICAL, INTENT(OUT)        :: linit

    linit = date%init

  END SUBROUTINE is_init
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  ! HELPER ROUTINES FOR DATE CONVERSIONS / TIME DISTANCE CALCULATIONS
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  INTEGER FUNCTION MonthLength(ky, km)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ky, km

    SELECT CASE(CAL_TYPE)
    CASE(CAL_JULIAN)
       MonthLength = JulianMonthLength(ky, km)
    CASE(CAL_360D)
       MonthLength = 30
    CASE DEFAULT
       MonthLength = 0
    END SELECT

  END FUNCTION MonthLength
  ! -------------------------------------------------------------------------
  
  ! -------------------------------------------------------------------------
  INTEGER FUNCTION JulianMonthLength(ky, km)
    !+
    !
    ! Get_JulianMonLen [function, integer]
    !    get the length of a months in a Julian year
    !    (
    !    year  [integer] input (Calendar year)
    !    month [integer] input (month of the year)
    !    )
    !
    !-
    IMPLICIT NONE
    INTRINSIC :: MOD

    INTEGER, INTENT(in) :: km, ky

    INTEGER :: idmax

    SELECT CASE(km)
    CASE(1,3,5,7,8,10,12);  idmax = 31
    CASE(4,6,9,11);                  idmax = 30
    CASE(2)
       IF ( (MOD(ky,4)==0 .AND. MOD(ky,100)/=0) .OR. MOD(ky,400)==0 ) THEN
          ! leap year found
          idmax = 29
       ELSE
          idmax = 28
       END IF

    CASE DEFAULT
       idmax = 0

    END SELECT
    JulianMonthLength = idmax

  END FUNCTION JulianMonthLength
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  INTEGER FUNCTION YearLength(yr)

    INTEGER, OPTIONAL, INTENT(IN) :: yr

    SELECT CASE(CAL_TYPE)
    CASE(CAL_JULIAN)
       IF (PRESENT(yr)) THEN
          YearLength= JulianYearLength(yr)
       ELSE
          YearLength= 365.2422_dp
       ENDIF
    CASE(CAL_360D)
       YearLength=360
    END SELECT

  END FUNCTION YearLength
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------  
  INTEGER FUNCTION JulianYearLength(yr)

    INTEGER, INTENT(IN) :: yr

    IF (yr == 1582) THEN
       JulianYearLength = 355
    ELSE IF ( (MOD(yr,4)==0 .AND. MOD(yr,100)/=0) .OR. MOD(yr,400)==0 ) THEN
       JulianYearLength = 366
    ELSE
       JulianYearLength = 365
    END IF

  END FUNCTION JulianYearLength
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  FUNCTION YearDay (date) RESULT (dayno)   

    ! returns the day of the year for a given date
    ! the seconds of the day are contained fractional

    TYPE(time_days) :: date    ! evaluate for this date
    REAL(dp)        :: dayno

    ! LOCAL
    INTEGER         :: yr, mo, dy, hr, mi, se, ms
    INTEGER         :: day, day01, day_diff
    REAL(DP)        :: nseconds
    REAL(DP)        :: frac
    INTEGER         :: status

    CALL timer_get_date(status,date,yr,mo,dy,hr,mi,se,ms)
    nseconds = 3600*hr+ 60* mi + se + ms/1000._dp

    SELECT CASE (CAL_TYPE)
    CASE (CAL_JULIAN)
       ! Julian date for current date
       day      = Julian_day(REAL(dy,dp),mo,yr)
       ! Julian date for 1st January (same year)
       day01    = Julian_day(REAL(1,dp),01,yr)
       day_diff = day-day01 +1 
       frac = (nseconds+0.000001_dp)/OneDay
       dayno = REAL(day_diff,dp) + frac
    CASE (CAL_360D)
       day   = 360*yr+(mo-1)*30+(dy-1)
       frac  = nseconds/OneDay
       dayno = REAL(day,dp) + frac
    END SELECT

  END FUNCTION YearDay
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  FUNCTION julian_day(DD, MM, YY)

    !
    ! [x] = the greatest integer that does not exceed x.
    !       For example, [-1.5]=-2. This is sometimes called the floor
    !       function (for example in C/C++). 
    ! INT(x) = [x] NOTE: some computer languages have a different definition.
    ! FIX(x) = the number x without its fraction. For example, FIX(-1.5)=-1 
    ! x\y    = FIX(x/y) 
    !
    ! 2.3.1 Gregorian Date to Julian Day Number
    ! For a Gregorian Date specified as day D (a real number),
    ! month M (an integer with January = 1), and year Y (any integer),
    ! the Julian Day Number JD can be calculated as follows: 
    ! IF M < 3 THEN 
    !       M = M + 12 
    !       Y=Y-1 
    ! END IF 
    ! JD = D + (153 * M - 457) \ 5 + 365 * Y + [Y / 4] - [Y / 100] + 
    !      [Y / 400] + 1721118.5
    !

    IMPLICIT NONE

    INTRINSIC :: FLOOR, REAL, INT

    ! I/O
    REAL(DP)             :: julian_day
    REAL(dp), INTENT(IN) :: DD
    INTEGER,  INTENT(IN) :: MM, YY

    ! LOCAL
    REAL(dp) :: D
    INTEGER  :: M, Y

    D = DD
    IF (MM < 3) THEN
       M = MM + 12
       Y = YY - 1
    ELSE
       M = MM
       Y = YY
    END IF

    julian_day = D &
         + INT(REAL((153 * M - 457), DP) / 5.0_DP) &
         + 365.0_DP * REAL(Y, DP) &
         + FLOOR(REAL(Y, DP) / 4.0_DP) &
         - FLOOR(REAl(Y, DP) / 100.0_DP) &
         + FLOOR(REAL(Y, DP) / 400.0_DP) &
         + 1721118.5_DP

  END FUNCTION julian_day
  ! -------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  SUBROUTINE time_span_s_real(dts               &
       , yy1, mo1, dy1, hr1, mi1, se1, ms1  &
       , yy2, mo2, dy2, hr2, mi2, se2, ms2 )

    IMPLICIT NONE

    ! I/O
    REAL(DP), INTENT(OUT) :: dts ! dime span [s]
    INTEGER,  INTENT(IN)  :: yy1, mo1, dy1, hr1, mi1, se1, ms1
    INTEGER,  INTENT(IN)  :: yy2, mo2, dy2, hr2, mi2, se2, ms2

    ! LOCAL
    REAL(dp) :: day_1, day_2

    day_1 = gregor2julian(yy1, mo1, dy1, hr1, mi1, se1, ms1)
    day_2 = gregor2julian(yy2, mo2, dy2, hr2, mi2, se2, ms2)

    dts = 86400.0_dp * (day_2 - day_1)

  END SUBROUTINE time_span_s_real
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  SUBROUTINE time_span_s_int(dts       &
       , yy1, mo1, dy1, hr1, mi1, se1  &
       , yy2, mo2, dy2, hr2, mi2, se2 )

    IMPLICIT NONE

    INTRINSIC :: NINT

    ! I/O
    INTEGER, INTENT(OUT) :: dts ! dime span [s]
    INTEGER, INTENT(IN)  :: yy1, mo1, dy1, hr1, mi1, se1
    INTEGER, INTENT(IN)  :: yy2, mo2, dy2, hr2, mi2, se2

    ! LOCAL
    REAL(dp) :: day_1, day_2

    day_1 = gregor2julian(yy1, mo1, dy1, hr1, mi1, se1)
    day_2 = gregor2julian(yy2, mo2, dy2, hr2, mi2, se2)

    dts = NINT(86400.0_dp * (day_2 - day_1))

  END SUBROUTINE time_span_s_int
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  SUBROUTINE time_span_d_s(dtd           &
       , yy1, mo1, dy1, hr1, mi1, se1  &
       , yy2, mo2, dy2, hr2, mi2, se2 )


    IMPLICIT NONE

    ! I/O
    REAL(dp), INTENT(OUT) :: dtd ! dime span [d]
    INTEGER,  INTENT(IN)  :: yy1, mo1, dy1, hr1, mi1, se1
    INTEGER,  INTENT(IN)  :: yy2, mo2, dy2, hr2, mi2, se2

    ! LOCAL
    REAL(dp) :: day_1, day_2

    day_1 = gregor2julian(yy1, mo1, dy1, hr1, mi1, se1)
    day_2 = gregor2julian(yy2, mo2, dy2, hr2, mi2, se2)

    dtd = (day_2 - day_1)

  END SUBROUTINE time_span_d_s
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  SUBROUTINE time_span_d_ms(dtd           &
       , yy1, mo1, dy1, hr1, mi1, se1, ms1  &
       , yy2, mo2, dy2, hr2, mi2, se2, ms2 )

    IMPLICIT NONE

    ! I/O
    REAL(dp), INTENT(OUT) :: dtd ! dime span [d]
    INTEGER,  INTENT(IN)  :: yy1, mo1, dy1, hr1, mi1, se1, ms1
    INTEGER,  INTENT(IN)  :: yy2, mo2, dy2, hr2, mi2, se2, ms2

    ! LOCAL
    REAL(dp) :: day_1, day_2

    day_1 = gregor2julian(yy1, mo1, dy1, hr1, mi1, se1, ms1)
    day_2 = gregor2julian(yy2, mo2, dy2, hr2, mi2, se2, ms2)

    dtd = (day_2 - day_1)

  END SUBROUTINE time_span_d_ms
  !--------------------------------------------------------------------------

  ! ------------------------------------------------------------------------- 
  FUNCTION gregor2julian(YY, MM, DD, hr, mi, se, ms) RESULT(julian_date)

    ! adapted from 'Numerical Recipes in Fortran'
    ! calculates Julian day and Julian date, output: Julian date
    ! [x]    = The greatest integer that does not exceed x, e.g., [-1.5]=-2.
    !          This is sometimes called the floor function, e.g., in C/C++.
    ! INT(x) = [x] NOTE: some computer languages have a different definition.
    !          HERE: INT = FIX
    ! FIX(x) = the number x without its fraction, e.g., FIX(-1.5)=-1.
    !          HERE: FIX = INT
    ! x\y    = FIX(x/y)
    !
    ! Gregorian Date to Julian Day Number
    ! The Julian Day Count is a uniform count of days from a remote epoch
    ! in the past:
    !  -4712     January 1,   12 hours Greenwich Mean Time 
    !                                  (Julian proleptic Calendar)
    ! = 4713 BCE January 1,   12 hours GMT                 
    !                                  (Julian proleptic Calendar)
    ! = 4714 BCE November 24, 12 hours GMT                 
    !                                  (Gregorian proleptic Calendar)).
    ! At this instant, the Julian Day Number is 0.
    ! see also: http://aa.usno.navy.mil/data/docs/JulianDate.html
    !
    ! For a Gregorian Date specified as
    !   day D   (a real number),
    !   month M (an integer with January = 1), and
    !   year Y  (any integer),
    ! the Julian Day Number JD can be calculated as follows:
    ! IF M < 3 THEN
    !   M = M + 12
    !   Y=Y-1
    ! END IF
    ! JD= D + (153*M - 457)\5 + 365*Y + [Y/4] - [Y/100] + [Y/400] + 1721118.5

    IMPLICIT NONE

    INTRINSIC :: FLOOR, REAL, INT

    ! I/O
    REAL(DP)             :: julian_day, julian_date
    INTEGER,  INTENT(IN) :: DD, MM, YY, hr, mi, se
    INTEGER, INTENT(IN), OPTIONAL :: ms

    ! LOCAL
    INTEGER  :: D, M, Y
    INTEGER :: millisec

    IF (PRESENT(ms)) THEN
       millisec = ms
    ELSE
       millisec = 0
    END IF

    D = DD
    IF (MM < 3) THEN
       M = MM + 12
       Y = YY - 1
    ELSE
       M = MM
       Y = YY
    END IF

    julian_day  = D &
         + INT(REAL((153 * M - 457), DP) / 5.0_DP) &
         + 365.0_DP * REAL(Y, DP) &
         + FLOOR(REAL(Y, DP) / 4.0_DP) &
         - FLOOR(REAL(Y, DP) / 100.0_DP) &
         + FLOOR(REAL(Y, DP) / 400.0_DP) &
         + 1721118.5_DP

    julian_date = REAL(julian_day, dp)                  &
         + REAL(hr, DP) / 24.0_DP                & ! hour day fraction
         + REAL(mi, DP) / (60.0_DP * 24.0_DP)    & ! minute day fraction
                                ! seconds day fraction
         + REAL(se, DP) / (60.0_DP * 60.0_DP * 24.0_DP) &
                                ! milliseconds day fraction
         + REAL(millisec, DP) / (1000.0_DP * 60.0_DP * 60.0_DP * 24.0_DP)

  END FUNCTION gregor2julian
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE julian2gregor(jdate, year, month, day, hrs, mins, secs, ms)

    ! Here julian is input as a Julian Date, and the routine outputs
    ! iyyy (year), mm (month), id (day), hh (hour), min (minute), and
    ! ss (second) of the Gregorian date
    !
    ! adapted from 'Numerical Recipes in Fortran'

    IMPLICIT NONE
    INTRINSIC :: INT, REAL

    ! I/O
    INTEGER, INTENT(OUT) :: year, month, day, hrs, mins, secs
    REAL(dp), INTENT(IN) :: jdate
    INTEGER, INTENT(OUT), OPTIONAL :: ms

    ! LOCAL
    INTEGER, PARAMETER :: IGREG = 2299161
    INTEGER            :: julian, ja,alpha,jb,jc,jd,je
    REAL(dp)           :: fract

    ! add 0.5 to julian date to make the decimals 0 at midnight, not noon
    fract  = jdate + 0.5_dp
    julian = INT(fract)

    ! treat fraction to calculate time
    fract = fract - REAL(julian,dp)  ! day fraction

    fract = fract * 24.0_dp ! hour fraction
    hrs   = INT(fract)

    fract = (fract - hrs) * 60.0_dp
    mins  = INT(fract)

    fract = (fract - mins) * 60.0_dp 
    secs  = INT(fract)

    IF (PRESENT(ms)) THEN
       fract = (fract - secs) * 1000.0_dp
       ms = INT(fract)
       fract = fract - REAL(INT(fract),dp)
       IF (fract >= 0.5_dp) THEN
          ms = ms + 1
       END IF
       IF (ms >= 1000) THEN
          ms = ms - 1000
          secs = secs + 1
       END IF
    ELSE
       fract = fract - REAL(INT(fract),dp)
       IF (fract >= 0.5_dp) THEN
          secs = secs + 1
       END IF
    END IF
    IF (secs >= 60) THEN
       mins = mins + 1
       secs = secs - 60
    END IF
    IF (mins >= 60) THEN
       hrs  = hrs + 1
       mins = mins - 60
    END IF

    ! treat integer julian day number for date
    IF (julian >= IGREG) THEN 
       ! correction because of cross-over to Gregorian Calendar in 1582
       alpha = int( ( (julian-1867216) - 0.25_dp) / 36524.25_dp )
       ja    = julian + 1 + alpha - int(0.25_dp*alpha)
    ELSE IF (julian < 0) THEN 
       ! make day number pos by adding int number of Julian centuries,
       ! subtract off at the end
       ja = julian + 36525 * ( 1 - julian/36525 )
    ELSE
       ja = julian
    END IF

    jb = ja + 1524
    jc = int(6680.0_dp + ( (jb-2439870) - 122.1_dp) / 365.25_dp)
    jd = 365 * jc + int(0.25_dp*jc)
    je = int( (jb-jd) / 30.6001_dp )

    day   = jb - jd - int(30.6001_dp*je)

    month = je - 1
    IF (month > 12) month = month - 12

    year = jc - 4715
    IF (month > 2)  year = year - 1
    IF (year <= 0)  year = year - 1
    if (julian < 0) year = year - 100 * ( 1 - julian / 36525)

  END SUBROUTINE julian2gregor
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  FUNCTION utc2lt(status, model_time, degree_lon)

    USE messy_main_constants_mem, ONLY: OneDay

    IMPLICIT NONE

    ! I/O
    REAL(dp)              :: utc2lt
    INTEGER,  INTENT(OUT) :: status
    REAL(dp), INTENT(IN)  :: model_time
    REAL(dp), INTENT(IN)  :: degree_lon
    ! LOCAL
    REAL(dp) :: dateback

    status = 0

    IF ( (degree_lon < 0.0_dp) .OR. (degree_lon >= 360.0_dp) ) THEN
       WRITE(*,*) 'ERROR in utc2lt: Longitude out of bounds (lon = ', &
            degree_lon, ')!'
       status = 1
       RETURN
    ELSE IF (degree_lon <= 180.0_dp) THEN
       dateback = 0._dp
    ELSE IF (degree_lon > 180.0_dp) THEN
       dateback = 1._dp
    ENDIF

    ! [s]       [s]           [°]    [°/day]                [s/day]
    utc2lt = model_time + (degree_lon/360.0_dp - dateback) * OneDay

  END FUNCTION utc2lt
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  ! ROUTINES FOR TIMER EXTERNAL USE
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE timer_set_date_str(status, strflag, yr, mo, dy, hr, mi, se, ms)

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    CHARACTER(LEN=*), INTENT(IN)  :: strflag
    INTEGER,          INTENT(IN)  :: yr, mo, dy, hr, mi, se
    INTEGER,          INTENT(IN), OPTIONAL :: ms

    INTEGER                       :: milliseconds

    IF (PRESENT(ms)) THEN
       milliseconds = ms
    ELSE
       milliseconds = 0
    END IF

    SELECT CASE(strflag)
    CASE('start')
       CALL date_set_components(yr, mo, dy, hr, mi, se, milliseconds &
            , start_date)

       ! SET DATE COMPONENTS
       YEAR_START   = yr
       MONTH_START  = mo
       DAY_START    = dy
       HOUR_START   = hr
       MINUTE_START = mi
       SECOND_START = se
       MILLISECOND_START = milliseconds

       ! SET JULIAN START DATE
       JULIAN_DATE_START = &
            gregor2julian( YEAR_START,MONTH_START, DAY_START &
            , HOUR_START,MINUTE_START,SECOND_START, MILLISECOND_START )

    CASE('previous')
       CALL date_set_components(yr, mo, dy, hr, mi, se, milliseconds &
            , previous_date)
    CASE('current')
       CALL date_set_components(yr, mo, dy, hr, mi, se, milliseconds &
            , current_date)
    CASE('next')
       CALL date_set_components(yr, mo, dy, hr, mi, se, milliseconds &
            , next_date)
    CASE('stop')
       CALL date_set_components(yr, mo, dy, hr, mi, se, milliseconds &
            , stop_date)
    CASE('resume')
       CALL date_set_components(yr, mo, dy, hr, mi, se, milliseconds &
            , resume_date)
    CASE DEFAULT
       status = 3443 ! unknown date 
       RETURN
    END SELECT

    status = 0

  END SUBROUTINE timer_set_date_str
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE timer_set_date_myd(status, my_date, yr, mo, dy, hr, mi, se, ms)

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    TYPE(time_days),  INTENT(OUT) :: my_date
    INTEGER,          INTENT(IN)  :: yr, mo, dy, hr, mi, se
    INTEGER,          INTENT(IN), OPTIONAL :: ms

    INTEGER                       :: milliseconds

    IF (PRESENT(ms)) THEN
       milliseconds = ms
    ELSE
       milliseconds = 0
    END IF

    CALL date_set_components(yr, mo, dy, hr, mi, se, milliseconds, my_date)
    status = 0

  END SUBROUTINE timer_set_date_myd
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE timer_set_date_str_ds_int(status, strflag, day, second)

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    CHARACTER(LEN=*), INTENT(IN)  :: strflag
    INTEGER,          INTENT(IN)  :: day, second

    ! LOCAL 
    INTEGER  ::  yr, mo, dy, hr, mi, se, ms

    ms = 0

    SELECT CASE(strflag)
    CASE('start')
       CALL date_set(day, second, start_date)
       CALL date_get_components(start_date, yr, mo, dy, hr, mi, se, ms, status)

       ! SET DATE COMPONENTS
       YEAR_START   = yr
       MONTH_START  = mo
       DAY_START    = dy
       HOUR_START   = hr
       MINUTE_START = mi
       SECOND_START = se
       MILLISECOND_START = ms

       ! SET JULIAN START DATE
       JULIAN_DATE_START = &
            gregor2julian( YEAR_START,MONTH_START, DAY_START &
            , HOUR_START,MINUTE_START,SECOND_START, MILLISECOND_START )

    CASE('previous')
       CALL date_set(day, second, previous_date)
    CASE('current')
       CALL date_set(day, second, current_date)
    CASE('next')
       CALL date_set(day, second, next_date)
    CASE('stop')
       CALL date_set(day, second, stop_date)
    CASE('resume')
       CALL date_set(day, second, resume_date)
    CASE DEFAULT
       status = 3443 ! unknown date 
       RETURN
    END SELECT

    status = 0
    
  END SUBROUTINE timer_set_date_str_ds_int
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE timer_set_date_str_ds_real(status, strflag, day, second)

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    CHARACTER(LEN=*), INTENT(IN)  :: strflag
    INTEGER,          INTENT(IN)  :: day
    REAL(DP),         INTENT(IN)  :: second

    ! LOCAL 
    INTEGER  ::  yr, mo, dy, hr, mi, se, ms

    se = FLOOR(second)
    ms = INT((second - AINT(second,dp) ) * 1000._dp)

    SELECT CASE(strflag)
    CASE('start')
       CALL date_set(day, second, start_date)
       CALL date_get_components(start_date, yr, mo, dy, hr, mi, se, ms, status)

       ! SET DATE COMPONENTS
       YEAR_START   = yr
       MONTH_START  = mo
       DAY_START    = dy
       HOUR_START   = hr
       MINUTE_START = mi
       SECOND_START = se
       MILLISECOND_START = ms

       ! SET JULIAN START DATE
       JULIAN_DATE_START = &
            gregor2julian( YEAR_START,MONTH_START, DAY_START &
            , HOUR_START,MINUTE_START,SECOND_START, MILLISECOND_START )

    CASE('previous')
       CALL date_set(day, second, previous_date)
    CASE('current')
       CALL date_set(day, second, current_date)
    CASE('next')
       CALL date_set(day, second, next_date)
    CASE('stop')
       CALL date_set(day, second, stop_date)
    CASE('resume')
       CALL date_set(day, second, resume_date)
    CASE DEFAULT
       status = 3443 ! unknown date 
       RETURN
    END SELECT

    status = 0
    
  END SUBROUTINE timer_set_date_str_ds_real
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE timer_set_date_myd_ds_int(status, my_date, day, second)

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    TYPE(time_days),  INTENT(OUT) :: my_date
    INTEGER,          INTENT(IN)  :: day, second

    CALL date_set(day, second, my_date)
    status = 0

  END SUBROUTINE timer_set_date_myd_ds_int
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE timer_set_date_myd_ds_real(status, my_date, day, second)

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    TYPE(time_days),  INTENT(OUT) :: my_date
    INTEGER,          INTENT(IN)  :: day
    REAL(DP),         INTENT(IN)  ::second

    CALL date_set(day, second, my_date)
    status = 0

  END SUBROUTINE timer_set_date_myd_ds_real
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE timer_get_date_str(status, strflag, yr, mo, dy, hr, mi, se, ms)

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    CHARACTER(LEN=*), INTENT(IN)  :: strflag
    INTEGER,          INTENT(OUT) :: yr, mo, dy, hr, mi, se
    INTEGER,          INTENT(OUT), OPTIONAL :: ms

    SELECT CASE(strflag)
    CASE('start')
       CALL date_get_components(start_date, yr, mo, dy, hr, mi, se, ms &
            , status)
    CASE('previous')
       CALL date_get_components(previous_date, yr, mo, dy, hr, mi, se, ms &
            , status)
    CASE('current')
       CALL date_get_components(current_date, yr, mo, dy, hr, mi, se, ms &
            , status)
    CASE('next')
       CALL date_get_components(next_date, yr, mo, dy, hr, mi, se, ms &
            , status)
    CASE('stop')
       CALL date_get_components(stop_date, yr, mo, dy, hr, mi, se, ms &
            , status)
    CASE('resume')
       CALL date_get_components(resume_date, yr, mo, dy, hr, mi, se, ms &
            , status)
    CASE('rerun_stop')
       CALL date_get_components(rerun_stop_date, yr, mo, dy, hr, mi, se, ms &
            , status)
    CASE DEFAULT
       status = 3443 ! unknown date 
       RETURN
    END SELECT
    
    !status = 0

  END SUBROUTINE timer_get_date_str
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE timer_get_date_myd(status, my_date, yr, mo, dy, hr, mi, se, ms)

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    TYPE(time_days),  INTENT(IN)  :: my_date
    INTEGER,          INTENT(OUT) :: yr, mo, dy, hr, mi, se
    INTEGER,          INTENT(OUT), OPTIONAL :: ms

    CALL date_get_components(my_date, yr, mo, dy, hr, mi, se, ms, status)

  END SUBROUTINE timer_get_date_myd
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE timer_set_calendar(status, strcal)

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    CHARACTER(LEN=*), INTENT(IN)  :: strcal    

    SELECT CASE(strcal)
       CASE('julian')
          CAL_TYPE = CAL_JULIAN
       CASE('days360')
          CAL_TYPE = CAL_360D
       CASE DEFAULT
          status = 3444 ! unknown calendar type
          RETURN
    END SELECT

    status = 0

  END SUBROUTINE timer_set_calendar
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE timer_get_calendar(status, ical)

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status
    INTEGER, INTENT(OUT) :: ical    

    ical = CAL_TYPE
    status = 0

  END SUBROUTINE timer_get_calendar
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE timer_set_delta_time(status, dt)

    IMPLICIT NONE

    ! I/O
    INTEGER,  INTENT(OUT) :: status
    REAL(dp), INTENT(IN)  :: dt

    delta_time = dt
    
    status = 0

  END SUBROUTINE timer_set_delta_time
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE timer_get_delta_time(status, dt, clock_id)

    IMPLICIT NONE

    ! I/O
    INTEGER,  INTENT(OUT) :: status
    REAL(dp), INTENT(OUT) :: dt
    INTEGER,  INTENT(IN), OPTIONAL  :: clock_id

    IF (.NOT. PRESENT(clock_id)) THEN
       dt = delta_time
    ELSE
       dt = clock(clock_id)%delta_time
    END IF
    
    status = 0

  END SUBROUTINE timer_get_delta_time
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE timer_add_date_int(status, add_seconds &
       , iyr, imo, idy, ihr, imi, ise               &
       , oyr, omo, ody, ohr, omi, ose)

    IMPLICIT NONE

    ! I/O
    INTEGER,  INTENT(OUT) :: status
    INTEGER,  INTENT(IN)  :: add_seconds ! seconds to add to date
    INTEGER,  INTENT(IN)  :: iyr, imo, idy, ihr, imi, ise
    INTEGER,  INTENT(OUT) :: oyr, omo, ody, ohr, omi, ose

    ! LOCAL
    TYPE(time_days)       :: my_date

    ! calculate date
    CALL date_set_components(iyr, imo, idy, ihr, imi, ise, 0, my_date)

    CALL add_date(0, add_seconds, my_date, status)
    IF (status /= 0) RETURN

    CALL date_get_components(my_date, oyr, omo, ody, ohr, omi, ose)

    status = 0

  END SUBROUTINE timer_add_date_int
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE timer_add_date_real(status, add_seconds, iyr &
       , imo, idy, ihr, imi, ise, ims       &
       , oyr, omo, ody, ohr, omi, ose, oms)

    IMPLICIT NONE

    ! I/O
    INTEGER,  INTENT(OUT) :: status
    REAL(DP), INTENT(IN)  :: add_seconds ! seconds to add to date
    INTEGER,  INTENT(IN)  :: iyr, imo, idy, ihr, imi, ise, ims
    INTEGER,  INTENT(OUT) :: oyr, omo, ody, ohr, omi, ose, oms

    ! LOCAL
    TYPE(time_days)       :: my_date

    ! calculate date
    CALL date_set_components(iyr, imo, idy, ihr, imi, ise, ims, my_date)

    CALL add_date(0, add_seconds, my_date, status)
    IF (status /= 0) RETURN

    CALL date_get_components(my_date, oyr, omo, ody, ohr, omi, ose, oms)

    status = 0

  END SUBROUTINE timer_add_date_real
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE timer_add_seconds_to_date_real(status, add_seconds, idate, odate)

    IMPLICIT NONE

    ! I/O
    INTEGER,  INTENT(OUT) :: status
    REAL(DP), INTENT(IN)  :: add_seconds ! seconds to add to date
    TYPE(time_days),  INTENT(IN)  :: idate
    TYPE(time_days),  INTENT(OUT) :: odate

    ! calculate date
    CALL copy_date(idate, odate, status)
    IF (status /= 0) RETURN
    
    CALL add_date(0, add_seconds, odate, status)
    IF (status /= 0) RETURN

    status = 0

  END SUBROUTINE timer_add_seconds_to_date_real
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE timer_set_lstart(llstart, clock_id)

    IMPLICIT NONE
    
    LOGICAL, INTENT(IN)           :: llstart
    INTEGER, INTENT(IN), OPTIONAL :: clock_id
    
    IF (PRESENT(clock_id)) THEN
       clock(clock_id)%lstart = llstart 
    ELSE
       lstart = llstart
    END IF

  END SUBROUTINE timer_set_lstart
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE timer_set_lfirst_cycle(llfirst_cycle, clock_id)

    IMPLICIT NONE
    
    LOGICAL, INTENT(IN)           :: llfirst_cycle
    INTEGER, INTENT(IN), OPTIONAL :: clock_id
    
    IF (PRESENT(clock_id)) THEN
       clock(clock_id)%lfirst_cycle = llfirst_cycle 
    ELSE
       lfirst_cycle = llfirst_cycle
    END IF

  END SUBROUTINE timer_set_lfirst_cycle
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE timer_set_lresume

    IMPLICIT NONE

    lresume = .TRUE.

  END SUBROUTINE timer_set_lresume
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE timer_get_lresume(lr)

    IMPLICIT NONE

    LOGICAL, INTENT(OUT) :: lr

    lr = lresume

  END SUBROUTINE timer_get_lresume
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  LOGICAL FUNCTION timer_get_l_rerun(clock_id)

    IMPLICIT NONE

    INTEGER, INTENT(IN), OPTIONAL :: clock_id
    
    IF (PRESENT(clock_id)) THEN
       timer_get_l_rerun = clock(clock_id)%l_rerun
    ELSE
       timer_get_l_rerun = clock(current_clock)%l_rerun
    END IF

  END FUNCTION timer_get_l_rerun
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE timer_set_time_step_len(l2tls, clock_id)

    IMPLICIT NONE

    LOGICAL, INTENT(IN) :: l2tls
    INTEGER, INTENT(IN), OPTIONAL :: clock_id

    IF (PRESENT(clock_id)) THEN
       IF (lstart) THEN
          clock(clock_id)%time_step_len = clock(clock_id)%delta_time
       ELSE
          IF (l2tls) THEN
             clock(clock_id)%time_step_len = clock(clock_id)%delta_time
          ELSE
             clock(clock_id)%time_step_len = 2.0_dp*clock(clock_id)%delta_time
          ENDIF
       END IF
    ELSE
       IF (lstart) THEN
          time_step_len = delta_time
       ELSE
          IF (l2tls) THEN
             time_step_len = delta_time
          ELSE
             time_step_len = 2.0_dp*delta_time
          ENDIF
       END IF
    END IF

  END SUBROUTINE timer_set_time_step_len
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE timer_get_time_step_len(dt, clock_id)

    IMPLICIT NONE

    REAL(dp), INTENT(OUT)          :: dt
    INTEGER,  INTENT(IN), OPTIONAL :: clock_id

    IF (PRESENT(clock_id)) THEN
       dt = clock(clock_id)%time_step_len
    ELSE
       dt = time_step_len   
    END IF

  END SUBROUTINE timer_get_time_step_len
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE timer_set_no_cycles(ncyc)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ncyc

    no_cycles = ncyc

  END SUBROUTINE timer_set_no_cycles
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE timer_get_no_cycles(ncyc)

    IMPLICIT NONE

    INTEGER, INTENT(out) :: ncyc
    
    ncyc = no_cycles

  END SUBROUTINE timer_get_no_cycles
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE timer_set_labort(la)

    IMPLICIT NONE

    LOGICAL, INTENT(IN) :: la
    
    labort = la

  END SUBROUTINE timer_set_labort
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE timer_get_labort(la)

    IMPLICIT NONE

    LOGICAL, INTENT(OUT) :: la
    
    la = labort

  END SUBROUTINE timer_get_labort
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE timer_get_init_step(clock_istep, clock_id)

    IMPLICIT NONE

    INTEGER, INTENT(OUT)           :: clock_istep
    INTEGER, INTENT(IN), OPTIONAL  :: clock_id
    
    IF (PRESENT(clock_id)) THEN
       clock_istep = clock(clock_id)%INIT_STEP
    ELSe
       clock_istep = clock(current_clock)%INIT_STEP
    END IF

  END SUBROUTINE timer_get_init_step
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  TYPE(time_days) FUNCTION timer_get_next_date(clock_id)

    IMPLICIT NONE

    INTEGER, INTENT(IN), OPTIONAL  :: clock_id
    
    IF (PRESENT(clock_id)) THEN
       timer_get_next_date = clock(clock_id)%next_date
    ELSE
       timer_get_next_date = clock(current_clock)%next_date
    END IF

  END FUNCTION timer_get_next_date
  ! -------------------------------------------------------------------------
  
  ! -------------------------------------------------------------------------
  SUBROUTINE eval_time_str(status, z_time_string, z_tuf, z_year, z_month, &
       z_day, z_hour, z_min, z_sec)

    ! crack common netcdf time string format into usable bits
    ! example: 'seconds since 2000-01-01 00:00:00'

    USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM
    USE messy_main_tools,         ONLY: strcrack, ucase
    IMPLICIT NONE

    INTRINSIC TRIM, ADJUSTL, LEN_TRIM

    ! I/O
    INTEGER, INTENT(OUT) :: status
    INTEGER, INTENT(OUT), OPTIONAL :: z_year, z_month, z_day, z_hour &
                                    , z_min, z_sec
    REAL(DP), INTENT(OUT) :: z_tuf
    CHARACTER(LEN=*),  INTENT(INOUT) :: z_time_string

    ! LOCAL
    INTEGER                      :: nosub
    CHARACTER(LEN=STRLEN_MEDIUM) :: tunit, helpvar, helpvar2 &
         , helpvar3, odate, otime, newtime
    CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(:), POINTER :: field => NULL()
    LOGICAL :: l_formtok   = .FALSE. ! time in time string format HH:MM:SS?
    LOGICAL :: l_uniformat = .TRUE.  ! change time string format to
                                     ! YYYY-MM-DD HH:MM:SS if different?

    status = 0

    ! crack time_unit string into field of nosub substrings
    ! <tunit> since <date> <time>
    ! 4: with time
    CALL strcrack(TRIM(z_time_string), " ", field, nosub)
    IF ((nosub /= 3) .AND. (nosub /= 4)) THEN
      WRITE(*,*) 'eval_time_str ERROR: '// &
        'time_string in incompatible format: '
      WRITE(*,*) '  ', z_time_string
      WRITE(*,*) '  (use YYYY-MM-DD or DD-MONTH-YYYY)'
      status = 1
      RETURN
    ENDIF

    ! determine time unit factor tuf <time>*tuf=time(s)
    tunit = TRIM(ADJUSTL(field(1)))
    CALL ucase(tunit)
    SELECT CASE (tunit)
    CASE ('SECONDS')
      z_tuf = 1.0_dp
    CASE ('MINUTES')
      z_tuf = 60.0_dp
    CASE ('HOURS')
      z_tuf = 3600.0_dp
    CASE ('DAYS')
      z_tuf = 86400.0_dp
    CASE DEFAULT
      WRITE(*,*) 'eval_time_str ERROR: '// &
        'time_string unit ', tunit, ' not recognized!'
      status = 1
      RETURN
    END SELECT

    ! needs to be here, 'field' later modified for time
    odate  = field(3)

    ! processing time
    IF (nosub == 4) THEN
      ! <units> <since> <date> <time>
      otime  = field(4)
      IF (ASSOCIATED(field)) THEN
         DEALLOCATE(field); NULLIFY(field)
      END IF
      CALL strcrack(otime, ":", field, nosub)
      IF (nosub == 3) THEN
        ! <hh>:<mm>:<ss>
        l_formtok = .TRUE.
        helpvar      = field(1)(1:2)
        IF (PRESENT(z_hour)) READ(helpvar, *) z_hour
        helpvar      = field(2)(1:2)
        IF (PRESENT(z_min))  READ(helpvar, *) z_min
        helpvar      = field(3)(1:2)
        IF (PRESENT(z_sec))  READ(helpvar, *) z_sec
        IF (l_uniformat) THEN
          newtime = TRIM(otime)
        ENDIF
      ELSEIF (nosub == 2) THEN
        ! <hh>:<mm>
        helpvar      = field(1)(1:2)
        IF (PRESENT(z_hour)) READ(helpvar, *) z_hour
        helpvar      = field(2)(1:2)
        IF (PRESENT(z_min))  READ(helpvar, *) z_min
        IF (l_uniformat) THEN
          newtime = TRIM(otime) // ':00'
        ENDIF
      ELSE
        WRITE(*,*) 'eval_time_str ERROR: '// &
          'TIME in time_string in incompatible format:'
        WRITE(*,*) '  ', z_time_string
        WRITE(*,*) '  (use hh:mm or hh:mm:ss)'
        status = 1
        RETURN
      ENDIF
    ELSEIF (nosub == 3) THEN
      ! <units> <since> <date>
      IF (PRESENT(z_hour)) z_hour = 0
      IF (PRESENT(z_min))  z_min  = 0
      IF (PRESENT(z_sec))  z_sec  = 0
      WRITE(*,*) 'eval_time_str: no time given, 00:00:00 assumed'
      IF (l_uniformat) THEN
        newtime = '00:00:00'
      ENDIF
    ENDIF

    IF (ASSOCIATED(field)) THEN
       DEALLOCATE(field); NULLIFY(field)
    END IF

    ! processing date
    CALL strcrack(odate, "-", field, nosub)
    IF (nosub /= 3) THEN
      WRITE(*,*) 'eval_time_str ERROR:'// &
        ' DATE in time_string in incompatible format:'
      WRITE(*,*) '  ', z_time_string
      WRITE(*,*) '  (use YYYY-MM-DD or DD-MONTH-YYYY)'
      status = 1
      RETURN
    ENDIF

    ! YYYY-MM-DD and DD-MONTH-YYYY accepted
    IF (LEN_TRIM(field(1)) == 2 .AND. LEN_TRIM(field(2)) == 3 .AND. &
         LEN_TRIM(field(3)) == 4) THEN
      ! DD-MONTH-YYYY (month is a string)
      WRITE(*,*) 'eval_time_str NOTE: DATE FORMAT DD-MONTH-YYYY'
      helpvar3     = field(3)(1:4)
      IF (PRESENT(z_year))  READ(helpvar3, *) z_year
      helpvar2     = TRIM(ADJUSTL(field(2)(1:3)))
      CALL ucase(helpvar2)
      SELECT CASE (TRIM(helpvar2))
      CASE ('JAN')
        helpvar2 = '01'
      CASE ('FEB')
        helpvar2 = '02'
      CASE ('MAR')
        helpvar2 = '03'
      CASE ('APR')
        helpvar2 = '04'
      CASE ('MAY')
        helpvar2 = '05'
      CASE ('JUN')
        helpvar2 = '06'
      CASE ('JUL')
        helpvar2 = '07'
      CASE ('AUG')
        helpvar2 = '08'
      CASE ('SEP')
        helpvar2 = '09'
      CASE ('OCT')
        helpvar2 = '10'
      CASE ('NOV')
        helpvar2 = '11'
      CASE ('DEC')
        helpvar2 = '12'
      CASE DEFAULT
        WRITE(*,*) 'eval_time_str ERROR:'// &
          ' month ', helpvar2, 'not recognized'// &
          ' (use three-letter acronym with DD-MONTH-YYYY format)'
        status = 1
        RETURN
      END SELECT

      IF (PRESENT(z_month)) READ(helpvar2, *) z_month
      helpvar     = field(1)(1:2)
      IF (PRESENT(z_day))   READ(helpvar, *) z_day

      IF (l_uniformat) THEN
        ! build new time string
        WRITE(*,*) 'eval_time_str NOTE: unifying time string to '// &
          'YYYY-MM-DD hh:mm:ss format'
        z_time_string = TRIM(tunit) // ' since ' // TRIM(helpvar3) // '-' // &
        TRIM(helpvar2) // '-' // TRIM(helpvar) // ' ' //  TRIM(newtime)
      ENDIF
    ELSE IF (LEN_TRIM(field(1)) == 4 .AND. LEN_TRIM(field(2)) == 2 .AND. &
      LEN_TRIM(field(3)) == 2) THEN
      ! YYYY-MM-DD
      !WRITE(*,*) 'eval_time_str NOTE: DATE FORMAT YYYY-MM-DD'
      helpvar      = field(1)(1:4)
      IF (PRESENT(z_year))  READ(helpvar, *) z_year
      helpvar2     = field(2)(1:2)
      IF (PRESENT(z_month))   READ(helpvar2, *) z_month
      helpvar3     = field(3)(1:2)
      IF (PRESENT(z_day))   READ(helpvar3, *) z_day

      IF (l_uniformat .AND. (.NOT. l_formtok)) THEN
        ! build new time string
        WRITE(*,*) 'eval_time_str NOTE: unifying time string to '// &
          'YYYY-MM-DD hh:mm:ss format'
        z_time_string = TRIM(tunit) // ' since ' // TRIM(helpvar) // '-' // &
          TRIM(helpvar2) // '-' // TRIM(helpvar3) // ' ' //  TRIM(newtime)
      ENDIF
    ELSE
      WRITE(*,*) 'eval_time_str ERROR:'// &
        ' DATE in time_string in incompatible format:'
      WRITE(*,*) '  ', z_time_string
      WRITE(*,*) '  (use YYYY-MM-DD or DD-MONTH-YYYY)'
      status = 1
      RETURN
    ENDIF

    IF (ASSOCIATED(field)) THEN
      DEALLOCATE(field)
      NULLIFY(field)
    ENDIF

  END SUBROUTINE eval_time_str
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE calc_sza(status, cossza, modeltime, t0str, deglon, deglat)

    ! calculate solar zenith angle from
    ! - model time (UTC)
    ! - time origin as string, e.g., 'seconds since 2000-01-01 00:00:00'
    !                             or 'minutes since 12-Feb-2007 01:12'
    ! - longitude in degrees
    ! - latitude in degrees

    USE messy_main_constants_mem, ONLY: OneDay, PI

    IMPLICIT NONE
    INTRINSIC :: SIN, COS

    INTEGER, INTENT(OUT) :: status
    REAL(DP), INTENT(OUT) :: cossza
    REAL(DP), INTENT(IN)  :: modeltime, deglon, deglat
    CHARACTER(LEN=*), INTENT(INOUT) :: t0str

    INTEGER  :: t0year, t0month, t0day, t0hour, t0min, t0sec ! start time vars
    REAL(DP) :: localtime, time0_jul, localtime_jul, firstjan_jul,&
                dayreal, sodecli, rad_lat
    REAL(DP), PARAMETER :: InitSpr = 80._DP !< first day of spring in NH [d]
    REAL(DP), PARAMETER :: cancer  = 23.441_DP*(PI/180._DP) !< Earth obliquity

    status = 0

    ! calculate local time from UTC
    localtime = utc2lt(status, modeltime, deglon)
    IF (status /= 0) THEN
      WRITE(*,*) 'ERROR in calc_sza: utc2lt'
      STOP
    ENDIF

    ! determine start time year, month, day, hour, minute, second
    ! time0_jul used as dummy argument, time factor output not needed here
    CALL eval_time_str(status, t0str, time0_jul, &
                       t0year, t0month, t0day, t0hour, t0min, t0sec)
    IF (status /= 0) THEN
      WRITE(*,*) 'ERROR in calc_sza: eval_time_str'
      STOP
    ENDIF

    ! time origin as Julian day
    time0_jul = gregor2julian(t0year, t0month, t0day, t0hour, t0min, t0sec)

    ! localtime as Julian day
    localtime_jul = time0_jul + localtime/OneDay

    ! start of current year as Julian day
    firstjan_jul = gregor2julian(t0year, 1, 1, 0, 0, 0)

    ! day as a real value (at start of spring dayreal = 0):
    ! = (julian date)-(1st January of current year)-(start of spring)
    dayreal = localtime_jul - firstjan_jul - InitSpr

    ! seasonal cycle:
    sodecli = cancer * SIN (2._DP * PI * dayreal / 365.25_DP)
    ! diurnal cycle of psi, the solar elevation angle
    rad_lat  = deglat * (PI/180.)
    cossza = SIN(rad_lat) * SIN(sodecli) &
      - COS(rad_lat) * COS(sodecli) * COS(2._DP * PI * dayreal)

  END SUBROUTINE calc_sza
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE main_timer_set_clock(clock_id)

    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: clock_id

    IF (current_clock == clock_id) RETURN

    INIT_STEP          => clock(clock_id)%INIT_STEP
    delta_time         => clock(clock_id)%delta_time
    time_step_len      => clock(clock_id)%time_step_len

    YEAR               => clock(clock_id)% YEAR             
    MONTH              => clock(clock_id)% MONTH            
    DAY                => clock(clock_id)% DAY              
    HOUR               => clock(clock_id)% HOUR             
    MINUTE             => clock(clock_id)% MINUTE           
    SECOND             => clock(clock_id)% SECOND           
    MILLISECOND        => clock(clock_id)% MILLISECOND      

    YEAR_START         => clock(clock_id)% YEAR_START       
    MONTH_START        => clock(clock_id)% MONTH_START      
    DAY_START          => clock(clock_id)% DAY_START        
    HOUR_START         => clock(clock_id)% HOUR_START       
    MINUTE_START       => clock(clock_id)% MINUTE_START     
    SECOND_START       => clock(clock_id)% SECOND_START     
    MILLISECOND_START  => clock(clock_id)% MILLISECOND_START

    YEAR_NEXT          => clock(clock_id)% YEAR_NEXT        
    MONTH_NEXT         => clock(clock_id)% MONTH_NEXT       
    DAY_NEXT           => clock(clock_id)% DAY_NEXT         
    HOUR_NEXT          => clock(clock_id)% HOUR_NEXT        
    MINUTE_NEXT        => clock(clock_id)% MINUTE_NEXT      
    SECOND_NEXT        => clock(clock_id)% SECOND_NEXT      
    MILLISECOND_NEXT   => clock(clock_id)% MILLISECOND_NEXT 

    current_time_step  => clock(clock_id)% current_time_step

    DAYOFYEAR          => clock(clock_id)%DAYOFYEAR
    previous_date      => clock(clock_id)%previous_date 
    current_date       => clock(clock_id)%current_date 
    next_date          => clock(clock_id)%next_date

    lstart             => clock(clock_id)%lstart
    lfirst_cycle       => clock(clock_id)%lfirst_cycle
    !lresume            => clock(clock_id)%lresume
    !lbreak             => clock(clock_id)%lbreak
    !lstop              => clock(clock_id)%lstop 

    l_rerun            => clock(clock_id)%l_rerun

    current_clock = clock_id

  END SUBROUTINE main_timer_set_clock
  ! -------------------------------------------------------------------------

! ************************************************************************
END MODULE messy_main_timer
! ************************************************************************
