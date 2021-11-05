!
! Interface to the facilities from ISO/IEC 9945-1:1990 sec. 5.1.2
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

module f90_unix_dirent
  use f90_unix_errno
  use f90_unix_tools, only : C0
  implicit none
  include 'f90_unix_dirent_const.inc'

  type DIR
     integer(dir_kind):: dir
  end type DIR
  
contains

  subroutine closedir(dirp, errno)
    type(dir), intent(inout) :: dirp
    integer, intent(out), optional :: errno

    if(present(errno)) errno = 0
    call set_errno

    call c_closedir(dirp%dir)

    if(present(errno)) then
       errno = get_errno()
    else
       if(get_errno() /= 0) then
          call perror("f90_unix_dirent::closedir")
          stop
       end if
    endif

  end subroutine closedir
  

  subroutine opendir(dirname, dirp, errno)
    character(len=*), intent(in) :: dirname
    type(dir), intent(inout) :: dirp
    integer, intent(out), optional :: errno

    if(present(errno)) errno = 0
    call set_errno

    call c_opendir(dirname//C0, dirp%dir)

    if(present(errno)) then
       errno = get_errno()
    else
       if(get_errno() /= 0) then
          call perror("f90_unix_dirent::opendir")
          stop
       end if
    endif

  end subroutine opendir
  

  subroutine readdir(dirp, name, lenname, errno)
    type(dir), intent(inout) :: dirp
    character(len=*), intent(out) :: name
    integer, intent(out) :: lenname
    integer, intent(out), optional :: errno

    if(present(errno)) errno = 0
    call set_errno

    call c_readdir(dirp%dir, name, lenname)

    if(present(errno)) then
       errno = get_errno()
    else
       if(get_errno() /= 0) then
          call perror("f90_unix_dirent::readdir")
          stop
       end if
    endif
  end subroutine readdir
  

  subroutine rewinddir(dirp, errno)
    type(dir), intent(inout) :: dirp
    integer, intent(out), optional :: errno

    if(present(errno)) errno = 0
    call set_errno

    call c_rewinddir(dirp%dir)

    if(present(errno)) then
       errno = get_errno()
    else
       if(get_errno() /= 0) then
          call perror("f90_unix_dirent::rewinddir")
          stop
       end if
    endif
  end subroutine rewinddir
  
end module f90_unix_dirent

