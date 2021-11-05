!
! Test module f90_unix_proc
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
program f90_unix_proc_tst
  use f90_unix_proc
  implicit none

  integer :: status, errno

  call system("ls", status, errno)
  call check_errno(errno, 'system("ls", status, errno)')

  call abort("Creating core dump");
contains

  subroutine handle_sig()
  end subroutine handle_sig
  
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

end program f90_unix_proc_tst
