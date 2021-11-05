!
! Interface to the facilities from ISO/IEC 9945-1:1990 sec. 5.2-5.5
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
module f90_unix_dir
  use f90_unix_errno
  use f90_unix_tools, only : C0
  implicit none
  include 'f90_unix_dir_const.inc'

contains

  subroutine chdir(path, errno)
    character(len=*), intent(in) :: path
    integer, intent(out), optional :: errno

    if(present(errno)) errno = 0
    call set_errno

    call c_chdir(path//C0)

    if(present(errno)) then
       errno = get_errno()
    else
       if(get_errno() /= 0) then
          call perror("f90_unix_dir::chdir")
          stop
       end if
    endif
  end subroutine chdir
  
  subroutine getcwd(path, lenpath, errno)
    character(len=*), intent(out) :: path
    integer, intent(out), optional :: lenpath, errno

    if(present(errno)) errno = 0
    call set_errno

    call c_getcwd(path)
    call fortranify_string(path, lenpath)

    if(present(errno)) then
       errno = get_errno()
    else
       if(get_errno() /= 0) then
          call perror("f90_unix_dir::getcwd")
          stop
       end if
    endif
  end subroutine getcwd
  
  subroutine link(existing, new, errno)
    character(len=*), intent(in) :: existing, new
    integer, intent(out), optional :: errno

    if(present(errno)) errno = 0
    call set_errno

    call c_link(existing//C0, new//C0)

    if(present(errno)) then
       errno = get_errno()
    else
       if(get_errno() /= 0) then
          call perror("f90_unix_dir::link")
          stop
       end if
    endif
  end subroutine link

  subroutine mkdir(path, mode, errno)
    character(len=*), intent(in) :: path
    integer(mode_kind), intent(in) :: mode
    integer, intent(out), optional :: errno

    if(present(errno)) errno = 0
    call set_errno

    call c_mkdir(path//C0, int(mode, mode_kind))

    if(present(errno)) then
       errno = get_errno()
    else
       if(get_errno() /= 0) then
          call perror("f90_unix_dir::mkdir")
          stop
       end if
    endif

  end subroutine mkdir
  
  subroutine mkfifo(path, mode, errno)
    character(len=*), intent(in) :: path
    integer(mode_kind), intent(in) :: mode
    integer, intent(out), optional :: errno

    if(present(errno)) errno = 0
    call set_errno

    call c_mkfifo(path//C0, int(mode, mode_kind))

    if(present(errno)) then
       errno = get_errno()
    else
       if(get_errno() /= 0) then
          call perror("f90_unix_dir::mkfifo")
          stop
       end if
    endif

  end subroutine mkfifo
  
  subroutine rmdir(path, errno)
    character(len=*), intent(in) :: path
    integer, intent(out), optional :: errno

    if(present(errno)) errno = 0
    call set_errno

    call c_rmdir(path//C0)

    if(present(errno)) then
       errno = get_errno()
    else
       if(get_errno() /= 0) then
          call perror("f90_unix_dir::rmdir")
          stop
       end if
    endif

  end subroutine rmdir
  
  subroutine unlink(path, errno)
    character(len=*), intent(in) :: path
    integer, intent(out), optional :: errno

    if(present(errno)) errno = 0
    call set_errno

    call c_unlink(path//C0)

    if(present(errno)) then
       errno = get_errno()
    else
       if(get_errno() /= 0) then
          call perror("f90_unix_dir::unlink")
          stop
       end if
    endif
  end subroutine unlink
  
end module f90_unix_dir

