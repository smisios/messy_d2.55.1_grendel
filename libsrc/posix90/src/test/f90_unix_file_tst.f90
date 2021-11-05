!
! Test module f90_unix_file
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
program f90_unix_file_tst
  use f90_unix_file
  use f90_unix_proc, only : system
  implicit none

  type(stat_t) :: buf
  integer :: errno

  call access("f90_unix_file_tst", X_OK, errno)
  call check_errno(errno, 'access("f90_unix_file_tst", X_OK, errno)')

  call system("touch tstfile")
  call chmod("tstfile", S_IXUSR, errno)
  call check_errno(errno, 'chmod("tstfile", S_IXUSR, errno)')

  call chown("tstfile", 100, 100, errno)
  call check_errno(errno, 'chown("tstfile", 100, 100, errno)',1)

  call stat("tstfile", buf, errno)
  call check_errno(errno, 'stat(tstfile, buf, errno)')
  print 1000, buf%st_dev, buf%st_ino, buf%st_mode, buf%st_nlink,&
         & buf%st_uid, buf%st_gid, buf%st_rdev, buf%st_size, buf%st_atime,&
         & buf%st_mtime, buf%st_ctime
1000 format("dev=", I4, " ino=", I8, " mode=", I6, " nlink=", I2, " uid=", I4,&
          &" gid=", I4, " rdev=" ,I8, " size=", I2, " atime", I10,&
          &" mtime", I10, " ctime", I10)
  call system("rm -f tstfile")
contains

  subroutine check_errno(errno, msg, flag)
    integer, intent(in) :: errno
    character(len=*), intent(in) :: msg
    integer,intent(in),optional:: flag

    if(errno/=0) then
       call perror(msg, errno)
       if(present(flag)) then
        print *, msg//' failed as it should'        
       else
          stop
       end if
    else
       if(present(flag)) then
          print *, msg//' succeeded, but should not. STOP'
          stop
       else
          print *, msg//' succeeded'
       end if       
    endif
  end subroutine check_errno

end program f90_unix_file_tst
