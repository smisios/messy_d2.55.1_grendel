module mo_getenv
#ifdef HAVE_CONFIG_H
#include "config.h"
  IMPLICIT NONE
#ifndef HAVE_FORTRAN_ROUTINE_GET_ENVIRONMENT_VARIABLE
  interface
     subroutine getenv_wrapper(name, value, length, status, trim_name, &
       trailing_blanks)
       character(len=*), intent(in) :: name
       character(len=*), intent(out) :: value
       integer, intent(out) :: length, status
       logical, intent(in) :: trim_name
       integer, intent(in) :: trailing_blanks
     end subroutine getenv_wrapper
  end interface
  LOGICAL, PARAMETER :: USE_GETENV_WRAPPER = .TRUE.
contains
  subroutine get_environment_variable(name, value, length, status, trim_name)
    character(len=*), intent(in) :: name
    character(len=*), optional, intent(out) :: value
    integer, optional, intent(out) :: length, status
    logical, optional, intent(in) :: trim_name
    character(len=1) :: value_dummy
    integer :: length_dummy, status_dummy, trailing_blanks
    logical :: trim_name_dummy
    trim_name_dummy = .true.
    if (present(trim_name)) trim_name_dummy = trim_name
    trailing_blanks = len(name) - len_trim(name)
    if (present(value)) then
       length_dummy = len(value)
       call getenv_wrapper(name, value, length_dummy, status_dummy, &
            trim_name_dummy, trailing_blanks)
    else
       length_dummy = len(value_dummy)
       call getenv_wrapper(name, value_dummy, length_dummy, status_dummy, &
            trim_name_dummy, trailing_blanks)
    end if
    if (present(status)) then
       status = status_dummy
       if (present(value) .and. len(value) < length_dummy) status = -1
    end if
    if (present(length)) length = length_dummy
  end subroutine get_environment_variable
#else
  LOGICAL, PARAMETER :: USE_GETENV_WRAPPER = .FALSE.
#endif
#else
  IMPLICIT NONE
  LOGICAL, PARAMETER :: USE_GETENV_WRAPPER = .FALSE.
#endif
end module mo_getenv
