MODULE mo_util_string
  !----------------------------------------------
  ! This module holds string conversion utilities
  !----------------------------------------------
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: tolower        ! Conversion   : 'ABCXYZ' -> 'abcxyz'
  PUBLIC :: toupper        ! Conversion   : 'abcxyz' -> 'ABCXYZ'
  PUBLIC :: char2          ! Conversion   : INTEGER  -> CHAR (LEN=2)
  PUBLIC :: separator      ! Format string: (/"-----...-----"/)
  PUBLIC :: translate
  PUBLIC :: capitalize
  public :: logical_to_string
  !-----------------
  ! module variables
  !-----------------
  CHARACTER(len=*) ,PARAMETER :: separator = '("'//REPEAT('-',79)//'")'
  !==============================================================================
CONTAINS
  !==============================================================================
  FUNCTION tolower (upper)
    !-----------------------------------
    ! Conversion: Uppercase -> Lowercase
    !-----------------------------------
    CHARACTER(LEN=*)              ,INTENT(in) :: upper
    CHARACTER(LEN=LEN_TRIM(upper))            :: tolower

    INTEGER            :: i
    INTEGER ,PARAMETER :: idel = ICHAR('a')-ICHAR('A')

    DO i=1,LEN_TRIM(upper)
       IF (ICHAR(upper(i:i)) >= ICHAR('A') .AND. &
            ICHAR(upper(i:i)) <= ICHAR('Z')) THEN
          tolower(i:i) = CHAR( ICHAR(upper(i:i)) + idel )
       ELSE
          tolower(i:i) = upper(i:i)
       END IF
    END DO

  END FUNCTION tolower
  !------------------------------------------------------------------------------
  FUNCTION toupper (lower)
    !-----------------------------------
    ! Conversion: Lowercase -> Uppercase
    !-----------------------------------
    CHARACTER(LEN=*)              ,INTENT(in) :: lower
    CHARACTER(LEN=LEN_TRIM(lower))            :: toupper

    INTEGER            :: i
    INTEGER ,PARAMETER :: idel = ICHAR('A')-ICHAR('a')

    DO i=1,LEN_TRIM(lower)
       IF (ICHAR(lower(i:i)) >= ICHAR('a') .AND. &
            ICHAR(lower(i:i)) <= ICHAR('z')) THEN
          toupper(i:i) = CHAR( ICHAR(lower(i:i)) + idel )
       ELSE
          toupper(i:i) = lower(i:i)
       END IF
    END DO

  END FUNCTION toupper
  !------------------------------------------------------------------------------
  FUNCTION char2 (i, zero)
    !----------------------------------------
    ! Conversion: INTEGER -> CHARACTER(LEN=2)
    !----------------------------------------
    CHARACTER(LEN=2)                       :: char2 ! result
    INTEGER          ,INTENT(in)           :: i     ! argument
    CHARACTER        ,INTENT(in) ,OPTIONAL :: zero  ! padding instead of '0'

    INTEGER ,PARAMETER :: i0 = ICHAR ('0')

    IF (i>99 .OR. i<0) THEN
       char2 = '**'
    ELSE
       char2(1:1) = CHAR(    i/10  + i0)
       char2(2:2) = CHAR(MOD(i,10) + i0)
    ENDIF

    IF(PRESENT(zero)) THEN
       IF(char2(1:1) == '0') char2(1:1) = zero
       IF(char2(2:2) == '0') char2(2:2) = zero
    ENDIF
  END FUNCTION char2
  !------------------------------------------------------------------------------

  FUNCTION translate(the_string, before, after)
    CHARACTER(*), INTENT(in) :: the_string
    CHARACTER(*), INTENT(in) :: before
    CHARACTER(LEN(before)), INTENT(in) :: after
    CHARACTER(LEN(the_string)) :: translate

    INTEGER :: offset, i, tr_index

    translate = the_string
    offset = 1
    DO WHILE(offset <= LEN(translate))
       i = SCAN(translate(offset:), before)
       IF(i == 0) THEN
          RETURN
       END IF
       offset = offset + i - 1
       ! This cannot fail!
       tr_index = INDEX(before, translate(offset:offset))
       translate(offset:offset) = after(tr_index:tr_index)
       offset = offset + 1
    END DO
  END FUNCTION translate

  FUNCTION capitalize(the_string)
    CHARACTER(*), INTENT(in) :: the_string
    CHARACTER(LEN(the_string)) :: capitalize
    capitalize = the_string
    if(len(the_string) > 0) then
        capitalize(1:1) = toupper(capitalize(1:1))
    end if
  END FUNCTION capitalize
  
  function logical_to_string(value)
    logical, intent(in) :: value
    character(5) :: logical_to_string
    logical_to_string = merge('true ', 'false', value)
  end function logical_to_string

END MODULE mo_util_string
