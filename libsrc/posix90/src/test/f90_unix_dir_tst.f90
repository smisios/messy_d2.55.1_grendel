!
! Test module f90_unix_dir
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
program f90_unix_dir_tst
  use f90_unix_dir
  use f90_unix_io
  use f90_unix_file, only :s_iwusr,s_irusr,s_ixusr
  implicit none

  character(len=128) :: name
  integer :: lenname
  integer :: errno
  type(FILE) :: fp

  call chdir("..", errno)
  call check_errno(errno, 'chdir("..")')

  call chdir("BAD", errno)
  call check_errno(errno, 'chdir("BAD")', 1)

  call chdir("test", errno)
  call check_errno(errno, 'chdir("test")')

  call getcwd(name, lenname, errno)
  call check_errno(errno, 'getcwd(name, lenname, errno)')
  print *, 'cwd is:>'//name(1:lenname)//'<'

  name = "testfile"
  fp = fopen(trim(name), "w"); call fclose(fp)
  call link(trim(name), trim(name)//'.lnk', errno)
  call check_errno(errno, "link(trim(name), trim(name)//'.lnk', errno)")

  call unlink(trim(name)//'.lnk', errno)
  call check_errno(errno, "unlink(trim(name)//'.lnk', errno)")
  call unlink(trim(name), errno)
  call check_errno(errno, "unlink(trim(name), errno)")

  name = "testdir"
  call mkdir(trim(name), s_iwusr+s_irusr+s_ixusr, errno)
  call check_errno(errno, "mkdir(trim(name), s_iwusr+s_irusr+s_ixusr, errno)")
  name = "testdir/testdir"
  call mkdir(trim(name), s_iwusr+s_irusr+s_ixusr, errno)
  call check_errno(errno, "mkdir(trim(name), s_iwusr+s_irusr+s_ixusr, errno)")

  name = "testdir"
  call rmdir(trim(name), errno)
  call check_errno(errno, "rmdir(trim(name), errno)",1)
  name = "testdir/testdir"
  call rmdir(trim(name), errno)
  call check_errno(errno, "rmdir(trim(name), errno)")
  name = "testdir"
  call rmdir(trim(name), errno)
  call check_errno(errno, "rmdir(trim(name), errno)")
  
  name = "testfifo"
  call mkfifo(trim(name), s_iwusr+s_irusr, errno)
  call check_errno(errno, "mkfifo(trim(name), s_iwusr+s_irusr, errno)")
  call unlink(trim(name), errno)
  call check_errno(errno, "unlink(trim(name), errno)")

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

end program f90_unix_dir_tst
