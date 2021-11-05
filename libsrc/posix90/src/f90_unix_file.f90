!
! Interface to the facilities from ISO/IEC 9945-1:1990 sec. 5.6
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
module f90_unix_file
  use f90_unix_errno
  use f90_unix_tools, only : C0
  use f90_unix_dir, only : mode_kind
  use f90_unix_env, only : uid_kind, gid_kind
  use f90_unix_time, only : time_kind
  implicit none
  include 'f90_unix_file_const.inc'
  type stat_t
     integer(dev_kind) ::     st_dev
     integer(ino_kind) ::     st_ino
     integer(mode_kind) ::    st_mode
     integer(nlink_kind) ::   st_nlink
     integer(uid_kind) ::     st_uid
     integer(gid_kind) ::     st_gid
     integer(dev_kind) ::     st_rdev
     integer(off_kind) ::     st_size
     integer(time_kind) ::    st_atime
     integer(time_kind) ::    st_mtime
     integer(time_kind) ::    st_ctime
  end type stat_t
  type utimbuf
     integer(time_kind) :: actime, modtime
  end type utimbuf
  
contains

  subroutine access(path, amode, errno)
    character(len=*), intent(in) :: path
    integer, intent(in) :: amode
    integer, intent(out) :: errno

    errno = 0
    call set_errno

    call c_access(path//C0, amode)

    errno = get_errno()
  end subroutine access

  subroutine chmod(path, mode, errno)
    character(len=*), intent(in) :: path
    integer(mode_kind), intent(in) :: mode
    integer, optional, intent(out) :: errno

    if(present(errno)) errno = 0
    call set_errno

    call c_chmod(path//C0, mode)

    if(present(errno)) then
       errno = get_errno()
    else
       if(get_errno() /= 0) then
          call perror("f90_unix_file::chmod")
          stop
       end if
    endif

  end subroutine chmod

  subroutine chown(path, owner, group, errno)
    character(len=*), intent(in) :: path
    integer(UID_KIND), intent(in) :: owner
    integer(GID_KIND), intent(in) :: group
    integer, optional, intent(out) :: errno

    if(present(errno)) errno = 0
    call set_errno

    call c_chown(path//C0, owner, group)

    if(present(errno)) then
       errno = get_errno()
    else
       if(get_errno() /= 0) then
          call perror("f90_unix_file::chown")
          stop
       end if
    endif

  end subroutine chown

  subroutine stat(path, buf, errno)
    character(len=*), intent(in) :: path
    type(stat_t), intent(out) :: buf
    integer, optional, intent(out) :: errno

    if(present(errno)) errno = 0
    call set_errno

    call c_stat(path//C0, buf%st_dev, buf%st_ino, buf%st_mode, buf%st_nlink,&
         & buf%st_uid, buf%st_gid, buf%st_rdev, buf%st_size, buf%st_atime,&
         & buf%st_mtime, buf%st_ctime)

    if(present(errno)) then
       errno = get_errno()
    else
       if(get_errno() /= 0) then
          call perror("f90_unix_file::stat")
          stop
       end if
    endif

  end subroutine stat

  
end module f90_unix_file
