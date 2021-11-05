! **********************************************************************+
MODULE messy_main_tools_bi
! **********************************************************************+

  ! THIS MODULE PROVIDES ROUTINES IN A PSEUDO-PARALLEL ENVIRONMENT
  
  IMPLICIT NONE
  PUBLIC

  INTERFACE start_message_bi
    MODULE PROCEDURE start_message_si
  END INTERFACE

  INTERFACE end_message_bi
    MODULE PROCEDURE end_message_si
  END INTERFACE

CONTAINS

  ! --------------------------------------------------------------------
  SUBROUTINE start_message_si(modstr, str, substr)

    USE messy_main_mpi_bi,    ONLY: p_parallel_io
    
    IMPLICIT NONE
    INTRINSIC :: TRIM
    
    ! I/O
    CHARACTER(LEN=*), INTENT(IN) :: modstr, str, substr
    
    IF (p_parallel_io) THEN
       WRITE (*,'(75("*"))')
       WRITE(*,*) '*** START ',TRIM(modstr),': ',TRIM(str), &
            ' (',TRIM(substr),')'
    END IF
    
  END SUBROUTINE start_message_si
  ! --------------------------------------------------------------------
  
  ! --------------------------------------------------------------------
  SUBROUTINE end_message_si(modstr, str, substr)
    
    ! ECHAM5
    USE messy_main_mpi_bi,    ONLY: p_parallel_io
    
    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! I/O
    CHARACTER(LEN=*), INTENT(IN) :: modstr, str, substr
    
    IF (p_parallel_io) THEN
       WRITE(*,*) '*** END   ',TRIM(modstr),': ',TRIM(str), &
            ' (',TRIM(substr),')'
       WRITE (*,'(75("*"))')
    END IF
    
  END SUBROUTINE end_message_si
  ! --------------------------------------------------------------------

  !---------------------------------------------------------------------
  FUNCTION find_next_free_unit(istart,istop) RESULT(unit)

    USE messy_main_mpi_bi, ONLY: finish

    ! I/O
    INTEGER :: istart, istop, unit

    ! LOCAL
    LOGICAL        :: found, opened
    INTEGER        :: i
    CHARACTER(256) :: info

    found = .FALSE.
    DO i=istart,istop
       INQUIRE(unit=i,opened=opened)
       IF (.NOT.opened) THEN
          unit = i
          found = .TRUE.
          EXIT
       END IF
    END DO

    IF (.NOT. found) THEN
       WRITE(info,'(a,i2.2,a,i2.2,a)') &
         'No unit in range <',istart,':',istop,'> free.'
       CALL finish('find_next_free_unit',info)
    END IF

  END FUNCTION find_next_free_unit
!-------------------------------------------------------------------------

! um_ak_20080328+
!-------------------------------------------------------------------------
  SUBROUTINE error(string, substr)

    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: string, substr

    WRITE (*,'(/,78("*"),/)')
    WRITE (*,'("ERROR: substr = ",A)')  TRIM(substr)
    WRITE (*,'("ERROR: ",A)')           TRIM(string)
    WRITE (*,'(/,78("*"),/)')

    STOP

  END SUBROUTINE error
!-------------------------------------------------------------------------
! um_ak_20080328-

! **********************************************************************+
END MODULE messy_main_tools_bi
! **********************************************************************+
