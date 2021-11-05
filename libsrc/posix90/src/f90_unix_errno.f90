!
! Interface to the facilities from ISO/IEC 9945-1:1990 sec. 2.4
!
!   This file is part of Posix90
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 51 Franklin Street, Fifth Floor,
!   Boston, MA 02110-1301, USA. 
!
!
module f90_unix_errno
  use f90_unix_tools
  implicit none

  include 'f90_unix_errno_const.inc'

contains

  character(len=80) function strerror(err, errno)
    integer, intent(in) :: err
    integer, intent(out), optional :: errno
    integer :: errc
    external c_strerror

    call c_strerror(strerror, err, errc)
    call fortranify_string(strerror)
    if(present(errno)) then
       errno = errc
    else
       if(errc/=0) then
          print *, "Unknown error code in strerror. stop."
          stop
       end if
    end if
    
  end function strerror


  subroutine perror(str, errc)
    character(len=*), intent(in) :: str
    integer, intent(in), optional :: errc
    integer err

    if(present(errc)) then
       err = errc
    else
       err = get_errno()
    end if
    
    call c_perror(str//C0, err)
  end subroutine perror
  

  integer function get_errno()
    integer, external :: c_get_errno

    get_errno = c_get_errno()
  end function get_errno


  subroutine set_errno(errc)
    integer, intent(in), optional :: errc
    external c_set_errno

    if(present(errc)) then
       call c_set_errno(errc)
    else
       call c_set_errno(0)
    end if
  end subroutine set_errno
  
  
end module f90_unix_errno
