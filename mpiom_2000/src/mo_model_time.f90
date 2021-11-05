MODULE mo_model_time
  USE mo_kind, ONLY: i8
  USE mo_io_config, ONLY: next_free_unit
  USE mo_mpi, ONLY: p_bcast
  USE mo_parallel, ONLY: p_pe, p_io

  IMPLICIT NONE 

  type calendar
    integer :: seconds_per_step
    integer :: days_per_year = 365
  end type

  TYPE time_desc
    SEQUENCE
    INTEGER :: year=0, month=0,  mday=0, &
               hour=0, minute=0, second=0
  END TYPE time_desc
  
  INTEGER, PRIVATE, PARAMETER :: basemonlen(12, 2) = &
       reshape((/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, &
       &          31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /), &
       &       (/12, 2/))
  
  !>
  !! Module variables.
  !!
  type(calendar), private, save :: cal

  !>
  !! Operator definitions.
  !!
  INTERFACE OPERATOR(+)
    MODULE PROCEDURE add_time
    MODULE PROCEDURE add_seconds
  END INTERFACE

  INTERFACE OPERATOR(-)
    MODULE PROCEDURE sub_time
    MODULE PROCEDURE sub_seconds
  END INTERFACE

  INTERFACE OPERATOR(.diff.)
    MODULE PROCEDURE time_diff
  END INTERFACE
  
  interface operator(<)
    module procedure time_lt
  end interface
  interface operator(<=)
    module procedure time_le
  end interface
  interface operator(==)
    module procedure time_eq
  end interface
  interface operator(>=)
    module procedure time_ge
  end interface
  interface operator(>)
    module procedure time_gt
  end interface
  interface operator(/=)
    module procedure time_ne
  end interface

CONTAINS

  !>
  !! Initialize model time settings.
  !!
  !! Sets the calendar information.
  !!
  subroutine model_time_init(seconds_per_step, days_per_year)

    integer, intent(in) :: seconds_per_step
    integer, intent(in) :: days_per_year

    cal = calendar(seconds_per_step, days_per_year)

  end subroutine model_time_init
  
  

  !> Convert model date and time to character string
  ELEMENTAL FUNCTION format_model_time(model_time) RESULT(cmodel_time)
    TYPE(time_desc), INTENT(IN) :: model_time
    CHARACTER(LEN=40) :: cmodel_time

    WRITE(cmodel_time,'(i5.4,"-",i2.2,"-",i2.2,1x,i2.2,":",i2.2,":",i2.2)') &
       model_time%year, model_time%month, model_time%mday, model_time%hour,  &
       model_time%minute, model_time%second
  END FUNCTION format_model_time

!------------------------------------------------------------------------------

  !>
  !! Check if year is a leap year.
  !!
  elemental function is_leap_year(year)
    integer, intent(in) :: year
    logical :: is_leap_year
    is_leap_year = ( modulo(year, 4) == 0 .and. modulo(year, 100) /= 0 ) .or. &
                   modulo(year, 400) == 0
  end function

  !>
  !! Determine month length (i.e. number of days) for specified month and year
  !!
  ELEMENTAL FUNCTION monlen(mon, year)
    INTEGER, INTENT(IN) :: mon
    INTEGER, INTENT(IN) :: year
    INTEGER :: monlen
    IF (cal%days_per_year .EQ. 365) THEN
       monlen = basemonlen(mon, 1)
    ELSE IF (cal%days_per_year .EQ. 360) THEN
       monlen = 30
    ELSE
       monlen = MERGE(basemonlen(mon, 2), basemonlen(mon, 1), is_leap_year(year))
    END IF
  END FUNCTION monlen

!------------------------------------------------------------------------------

  !>
  !! Calculate number of days between 1st day of monstart and
  !! last day of monend inclusive
  !!
  ELEMENTAL FUNCTION monlen_sum(monstart, monend, year)
    INTEGER, INTENT(IN) :: monstart, monend
    INTEGER, INTENT(IN) :: year
    INTEGER :: monlen_sum, i
    monlen_sum = 0
    DO i = monstart, monend
       monlen_sum = monlen_sum + monlen(i, year)
    END DO
  END FUNCTION monlen_sum

!------------------------------------------------------------------------------

  !> Calculate number of days between January 1st of yearstart and
  !> last December day (31/30) of yearend inclusive
  ELEMENTAL FUNCTION yearlen_sum(yearstart, yearend)
    INTEGER :: yearlen_sum
    INTEGER, INTENT(IN) :: yearstart, yearend
    INTEGER :: year
    yearlen_sum = 0
    DO year = yearstart, yearend
       yearlen_sum = yearlen_sum + monlen_sum(1, 12, year)
    END DO
  END FUNCTION yearlen_sum

!------------------------------------------------------------------------------

  !> Add one day
  SUBROUTINE inc_day(year, month, day)
    INTEGER, INTENT(INOUT) :: year, month, day

    IF (month .EQ. 12 .AND. day .EQ. monlen(12, year)) THEN
       year = year + 1
       month = 1
       day = 1
    ELSE IF (day .EQ. monlen(month, year)) THEN
       day = 1
       month = month + 1
    ELSE
       day = day + 1
    END IF
  END SUBROUTINE inc_day

!------------------------------------------------------------------------------

  !> Write file model_date.asc at the end of run with start date and end date 
  !> of the current run and start date of the next run
  SUBROUTINE write_model_date(time1, time2, time3)

    TYPE(time_desc), INTENT(IN) :: time1, time2, time3
    INTEGER :: io_out_date
    CHARACTER(LEN=1), PARAMETER :: del="-"

    io_out_date = next_free_unit()

    OPEN(io_out_date,file='model_date.asc',status='unknown',                  &
         form='formatted',action='write')

    WRITE(io_out_date,'(3(i8.4,a1,i2.2,a1,i2.2))') &
       time1%year, del, time1%month, del, time1%mday, &
       time2%year, del, time2%month, del, time2%mday, &
       time3%year, del, time3%month, del, time3%mday

    CLOSE(io_out_date)

  END SUBROUTINE write_model_date

!------------------------------------------------------------------------------
  
  !> Convert model time to seconds since
  !>   i)  begin of year startyear
  !>   ii) begin of current year if startyear not specified
  ELEMENTAL INTEGER(i8) FUNCTION model_seconds(req_time, startyear)

    TYPE(time_desc), INTENT(IN)   :: req_time
    INTEGER, INTENT(IN), OPTIONAL :: startyear

    IF (PRESENT(startyear)) THEN

      model_seconds = INT(yearlen_sum(startyear, req_time%year - 1) &
                    + monlen_sum(1,req_time%month - 1,req_time%year) &
                    + req_time%mday - 1,i8)   * 86400_i8 &
                    + INT(req_time%hour,i8)   * 3600_i8 &
                    + INT(req_time%minute,i8) * 60_i8 &
                    + INT(req_time%second,i8)
    ELSE
      model_seconds = INT(monlen_sum(1,req_time%month - 1,req_time%year) &
                    + req_time%mday - 1,i8 )  * 86400_i8 &
                    + INT(req_time%hour,i8)   * 3600_i8 &
                    + INT(req_time%minute,i8) * 60_i8 &
                    + INT(req_time%second,i8)
     ENDIF

  END FUNCTION model_seconds

!------------------------------------------------------------------------------
  
  !> Calculate number of seconds between two times
  ELEMENTAL INTEGER(i8) FUNCTION seconds_between_two_times(time1, time2)

    TYPE(time_desc), INTENT(IN)   :: time1, time2
    INTEGER(i8) :: time1_seconds, time2_seconds
    INTEGER :: ref_year

    ref_year=MIN(time1%year,time2%year)
    time1_seconds = model_seconds(time1, ref_year)
    time2_seconds = model_seconds(time2, ref_year)
    seconds_between_two_times = time1_seconds - time2_seconds

  END FUNCTION seconds_between_two_times

!------------------------------------------------------------------------------

  !> Add time interval add_sec (in seconds) to model time time1 and 
  !> determine new date and time time2
  FUNCTION add_seconds(time1, add_sec) RESULT(time2)

    TYPE(time_desc), INTENT(IN)  :: time1
    INTEGER(i8),     INTENT(IN)  :: add_sec
    TYPE(time_desc) :: time2, time

    time = time_desc(0,0,0,INT(add_sec/3600_i8),&
           INT(MODULO(add_sec/60_i8,60_i8)),INT(MODULO(add_sec,60_i8)))
           
    time2 = time1 + time

  END FUNCTION add_seconds

!------------------------------------------------------------------------------

  !> Substract time interval sub_sec (in seconds) from model time time1 and 
  !> determine new date and time time2
  FUNCTION sub_seconds(time1, sub_sec) RESULT(time2)

    TYPE(time_desc), INTENT(IN)  :: time1
    INTEGER(i8),     INTENT(IN)  :: sub_sec
    TYPE(time_desc) :: time2, time

    time = time_desc(0,0,0,INT(sub_sec/3600_i8),&
           INT(MODULO(sub_sec/60_i8,60_i8)),INT(MODULO(sub_sec,60_i8)))
           
    time2 = time1 - time

  END FUNCTION sub_seconds

!------------------------------------------------------------------------------

  !> Add time interval dtime to model time time1 and 
  !> determine new date and time time2
  FUNCTION add_time(time1, dtime) RESULT(time2)

    TYPE(time_desc), INTENT(IN)  :: time1, dtime
    TYPE(time_desc) :: time2
    INTEGER :: add_minute, add_hour, add_day, add_year, day_in_month

    time2%second = time1%second + dtime%second
    add_minute = time2%second / 60 + dtime%minute 
    time2%second = MODULO(time2%second,60)

    time2%minute = time1%minute + add_minute
    add_hour = time2%minute / 60 + dtime%hour
    time2%minute = MODULO(time2%minute,60)
    
    time2%hour = time1%hour + add_hour
    add_day = time2%hour / 24 + dtime%mday
    time2%hour = MODULO(time2%hour,24)

    time2%mday  = time1%mday + add_day
    time2%month = time1%month
    time2%year = time1%year

    DO
      day_in_month = monlen(time2%month, time2%year)
      IF (time2%mday .GT. day_in_month) THEN
        time2%mday = time2%mday - day_in_month
        time2%month = time2%month + 1
        IF (time2%month .GT. 12) THEN
          time2%month = 1
          time2%year = time2%year + 1
        ENDIF
      ELSE
        EXIT
      ENDIF
    ENDDO

    time2%month = time2%month + dtime%month
    add_year = (time2%month - 1) / 12 + dtime%year
    time2%month = MODULO(time2%month - 1,12) + 1

    time2%year  = time2%year + add_year

  END FUNCTION add_time

!------------------------------------------------------------------------------

  !> Substract time interval dtime from model time time1 and 
  !> determine new date and time time2
  ELEMENTAL FUNCTION sub_time(time1, dtime) RESULT(time2)

    TYPE(time_desc), INTENT(IN)  :: time1, dtime
    TYPE(time_desc) :: time2
    INTEGER, PARAMETER :: minref = 1

    time2%year = time1%year - dtime%year

    time2%month = time1%month - dtime%month
    CALL adjust_month_neg(time2, minref)

    time2%mday = time1%mday - dtime%mday
    CALL adjust_day_neg(time2, minref)

    time2%hour = time1%hour - dtime%hour
    CALL adjust_hour_neg(time2, minref)

    time2%minute = time1%minute - dtime%minute
    CALL adjust_minute_neg(time2, minref)

    time2%second = time1%second - dtime%second
    CALL adjust_second_neg(time2, minref)

  END FUNCTION sub_time

!------------------------------------------------------------------------------

  !> Calculate time difference dtime between time1 and time2
  !> (time1 > time2 currently)
  ELEMENTAL FUNCTION time_diff(time1, time2) RESULT(dtime)

    TYPE(time_desc), INTENT(IN)  :: time1, time2
    TYPE(time_desc) :: dtime
    INTEGER, PARAMETER :: minref = 0

    dtime%year = time1%year -time2%year

    dtime%month = time1%month - time2%month
    CALL adjust_month_neg(dtime, minref)

    dtime%mday = time1%mday - time2%mday
    CALL adjust_day_neg(dtime, minref)

    dtime%hour = time1%hour - time2%hour
    CALL adjust_hour_neg(dtime, minref)

    dtime%minute = time1%minute - time2%minute
    CALL adjust_minute_neg(dtime, minref)

    dtime%second = time1%second - time2%second
    CALL adjust_second_neg(dtime, minref)

  END FUNCTION time_diff

!------------------------------------------------------------------------------

  PURE SUBROUTINE adjust_month_neg(time, minref)
    TYPE(time_desc), INTENT(INOUT) :: time
    INTEGER, INTENT(IN) :: minref

    DO
      IF (time%month .LT. minref) THEN
        time%month = time%month + 12
        time%year = time%year -1
      ELSE
        EXIT
      ENDIF
    ENDDO
  END SUBROUTINE adjust_month_neg

  PURE SUBROUTINE adjust_day_neg(time, minref)
    TYPE(time_desc), INTENT(INOUT) :: time
    INTEGER, INTENT(IN) :: minref
    
    INTEGER :: day_in_month
    DO
      IF (time%mday .LT. minref) THEN
         time%month = time%month - 1
         CALL adjust_month_neg(time, minref)
         day_in_month = monlen(time%month, time%year)
         time%mday = time%mday + day_in_month
      ELSE
        EXIT
      ENDIF
      END DO
  END SUBROUTINE adjust_day_neg

  PURE SUBROUTINE adjust_hour_neg(time, minref)
    TYPE(time_desc), INTENT(INOUT) :: time
    INTEGER, INTENT(IN) :: minref
    
    DO
      IF (time%hour .LT. 0) THEN
         time%hour = time%hour + 24
         time%mday = time%mday -1
         CALL adjust_day_neg(time, minref)
      ELSE
        EXIT
      ENDIF
    END DO
  END SUBROUTINE adjust_hour_neg

  PURE SUBROUTINE adjust_minute_neg(time, minref)
    TYPE(time_desc), INTENT(INOUT) :: time
    INTEGER, INTENT(IN) :: minref
      
    DO
      IF (time%minute .LT. 0) THEN
         time%minute = time%minute + 60
         time%hour = time%hour -1
         CALL adjust_hour_neg(time, minref)
      ELSE
        EXIT
      ENDIF
    END DO
  END SUBROUTINE adjust_minute_neg

  PURE SUBROUTINE adjust_second_neg(time, minref)
    TYPE(time_desc), INTENT(INOUT) :: time
    INTEGER, INTENT(IN) :: minref
    
    DO
      IF (time%second .LT. 0) THEN
         time%second = time%second + 60
         time%minute = time%minute -1
         CALL adjust_minute_neg(time, minref)
      ELSE
        EXIT
      ENDIF
    END DO
  END SUBROUTINE adjust_second_neg

!------------------------------------------------------------------------------
  
  !>
  !! Compare two time descriptions.
  !!
  !! Returns a value < 0, if this is before that, 
  !! a value > 0 if this is after that, or 0 if this is equal to that. 
  !!
  elemental function time_cmp(this, that)
    type(time_desc), intent(in) :: this, that
    integer :: time_cmp
    time_cmp = this%year - that%year
    if(time_cmp /= 0) return
    time_cmp = this%month - that%month
    if(time_cmp /= 0) return
    time_cmp = this%mday - that%mday
    if(time_cmp /= 0) return
    time_cmp = this%hour - that%hour
    if(time_cmp /= 0) return
    time_cmp = this%minute - that%minute
    if(time_cmp /= 0) return
    time_cmp = this%second - that%second
  end function time_cmp
  
  !>
  !! Check if one time is before (less than) the other.
  !!
  elemental function time_lt(this, that)    
    type(time_desc), intent(in) :: this, that
    logical :: time_lt
    time_lt = time_cmp(this, that) < 0
  end function time_lt  
  
  !>
  !! Check if one time is before (less than) or equal to the other.
  !!
  elemental function time_le(this, that)    
    type(time_desc), intent(in) :: this, that
    logical :: time_le
    time_le = time_cmp(this, that) <= 0
  end function time_le  
  
  !>
  !! Check if two times are equal.
  !!
  elemental function time_eq(this, that)
    type(time_desc), intent(in) :: this, that
    logical :: time_eq
    time_eq = time_cmp(this, that) == 0
  end function time_eq  

  !>
  !! Check if one time is after (greater than) or equal to the other.
  !!
  elemental function time_ge(this, that)    
    type(time_desc), intent(in) :: this, that
    logical :: time_ge
    time_ge = time_cmp(this, that) >= 0
  end function time_ge  
  
  !>
  !! Check if one time is after (greater than) the other.
  !!
  elemental function time_gt(this, that)    
    type(time_desc), intent(in) :: this, that
    logical :: time_gt
    time_gt = time_cmp(this, that) > 0
  end function time_gt  
  
  !>
  !! Check if two times are not equal.
  !!
  elemental function time_ne(this, that)    
    type(time_desc), intent(in) :: this, that
    logical :: time_ne
    time_ne = time_cmp(this, that) /= 0
  end function time_ne  
  
!------------------------------------------------------------------------------

  !>
  !! Check if we are at end of year.
  !!
  function is_end_of_year(this)
     type(time_desc), intent(in) :: this
     logical :: is_end_of_year
     type(time_desc) :: next, line
     next = add_seconds( this, int(cal%seconds_per_step, i8) )
     line = time_desc(next%year, 1, 1, 0, 0, 0)
     is_end_of_year = this < line .and. line <= next
  end function is_end_of_year
  
  !>
  !! Get current step in year, beginning with 1.
  !!
  function get_step_of_year(this)
     type(time_desc), intent(in) :: this
     integer :: get_step_of_year
     get_step_of_year = int( model_seconds(this)/ &
                             int(cal%seconds_per_step, i8) ) + 1
  end function get_step_of_year
  
  !> 
  !! Check if we are at end of month.
  !!
  function is_end_of_month(this)
     type(time_desc), intent(in) :: this
     logical is_end_of_month
     type(time_desc) :: next, line
     next = add_seconds( this, int(cal%seconds_per_step, i8) )
     line = time_desc(next%year, next%month, 1, 0, 0, 0)
     is_end_of_month = this < line .and. line <= next
  end function is_end_of_month
  
  !> 
  !! Get current step of month, beginning with 1.
  !!
  function get_step_of_month(this)
     type(time_desc), intent(in) :: this
     integer :: get_step_of_month
     type(time_desc) :: line
     line = time_desc(this%year, this%month, 1, 0, 0, 0)
     get_step_of_month = int( seconds_between_two_times(this, line)/ &
                              int(cal%seconds_per_step, i8) ) + 1
  end function get_step_of_month
  
  !> 
  !! Check if we are at end of day.
  !!
  function is_end_of_day(this)
     type(time_desc), intent(in) :: this
     logical :: is_end_of_day
     integer :: steps_per_day
     steps_per_day = 86400/cal%seconds_per_step
     is_end_of_day = get_step_of_day(this) == steps_per_day
  end function is_end_of_day
  
  !>
  !! Get current step of day, beginning with 1.
  !!
  function get_step_of_day(this)
     type(time_desc), intent(in) :: this
     integer :: get_step_of_day
     type(time_desc) :: line
     line = time_desc(this%year, this%month, this%mday, 0, 0, 0)
     get_step_of_day = int( seconds_between_two_times(this, line)/ &
                            int(cal%seconds_per_step, i8) ) + 1
  end function get_step_of_day
  
  !>
  !! Check if calendar supports a certain part of a day.
  !!
  function is_partial_day_supported(parts)
     integer :: parts
     logical :: is_partial_day_supported
     is_partial_day_supported = modulo(86400, parts) == 0 .and. & 
                                modulo(86400/parts, cal%seconds_per_step) == 0
  end function is_partial_day_supported
  
  !> 
  !! Check if we are at end of some part of a day.
  !!
  function is_end_of_partial_day(this, parts)
     type(time_desc), intent(in) :: this
     integer, intent(in) :: parts
     logical is_end_of_partial_day
     integer :: steps_per_part
     steps_per_part = 86400/parts/cal%seconds_per_step
     is_end_of_partial_day = &
          get_step_of_partial_day(this, parts) ==  steps_per_part
  end function is_end_of_partial_day
  
  !>
  !! Get current step of day, beginning with 1.
  !!
  function get_step_of_partial_day(this, parts)
     type(time_desc), intent(in) :: this
     integer, intent(in) :: parts
     integer :: get_step_of_partial_day
     type(time_desc) :: line
     integer :: seconds_per_part
     integer :: seconds_of_day
     integer :: seconds_of_part
     seconds_per_part = 86400/parts
     line = time_desc(this%year, this%month, this%mday, 0, 0, 0)
     seconds_of_day = int( seconds_between_two_times(this, line) )
     seconds_of_part = modulo(seconds_of_day, seconds_per_part) 
     get_step_of_partial_day = seconds_of_part/cal%seconds_per_step + 1
  end function get_step_of_partial_day
  

!------------------------------------------------------------------------------
  
  !> Convert integer array time_arr to structure type(time_desc) time
  FUNCTION array2time_desc(time_arr) RESULT(time)

    INTEGER, INTENT(IN) :: time_arr(6)
    TYPE(time_desc)     :: time

    time%year   = time_arr(1)
    time%month  = time_arr(2)
    time%mday   = time_arr(3)
    time%hour   = time_arr(4)
    time%minute = time_arr(5)
    time%second = time_arr(6)

  END FUNCTION array2time_desc 

!------------------------------------------------------------------------------
  
  !> Convert structure type(time_desc) time to integer array time_arr
  FUNCTION time_desc2array(time) RESULT(time_arr)

    TYPE(time_desc), INTENT(IN)     :: time
    INTEGER :: time_arr(6)

    time_arr(1) = time%year
    time_arr(2) = time%month
    time_arr(3) = time%mday
    time_arr(4) = time%hour
    time_arr(5) = time%minute
    time_arr(6) = time%second

  END FUNCTION time_desc2array

!------------------------------------------------------------------------------

  SUBROUTINE broadcast_time_desc(time)

    TYPE(time_desc) :: time
    INTEGER         :: buf(6)

    IF (p_pe == p_io) THEN
      buf = time_desc2array(time)
    ENDIF

    CALL p_bcast(buf,p_io)

    time = array2time_desc(buf)

  END SUBROUTINE broadcast_time_desc

END MODULE mo_model_time
