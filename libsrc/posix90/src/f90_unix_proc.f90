!
! Interface to the facilities from ISO/IEC 9945-1:1990 sec. 3
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
module f90_unix_proc
  use f90_unix_errno
  use f90_unix_tools, only : C0
  implicit none
  include 'f90_unix_proc_const.inc'

contains

  subroutine abort(msg)
    character(len=*), intent(in) :: msg
    print '(A)', msg
    call c_abort()
  end subroutine abort

  subroutine alarm(seconds, routine, secleft)
    integer, intent(in) :: seconds
    interface
       subroutine routine()
       end subroutine routine
    end interface
    optional routine
    integer, optional, intent(out) :: secleft
  end subroutine alarm
  
  subroutine system(cmd, status, errno)
    character(len=*) :: cmd
    integer, intent(out), optional :: status, errno
    integer :: stat

    if(present(errno)) errno = 0
    call set_errno

    call c_system(cmd//C0, stat)
    if(present(status)) status = stat

    if(present(errno)) then
       errno = get_errno()
    else
       if(get_errno() /= 0) then
          call perror("f90_unix_proc::system")
          stop
       end if
    endif

  end subroutine system
  
end module f90_unix_proc

