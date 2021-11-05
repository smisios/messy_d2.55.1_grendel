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
program f90_unix_signal_tst
  use f90_unix_signal
  implicit none

  integer :: status, errno
  type(sigset_type) :: set
  type(sigaction_type) :: action, oldaction
  external sig_handler, c_sig_handler

  call sigfillset(set)
  call sigemptyset(set)
  call sigaddset(set, SIGINT, SIGQUIT)
  print *, sigismember(set, SIGINT)
  call sigdelset(set, SIGINT, SIGQUIT)
  print *, sigismember(set, SIGINT)

  call sigaddset(set, SIGINT, SIGQUIT)
  action = sigaction_compile(sig_handler, set)
  call sigaction(SIGUSR1, action, oldaction)
  print *, 'Installed signal ', sigusr1
  print *, 'Old action was: ', oldaction%sa_handler
  call sigaction(SIGUSR1, oldaction=oldaction)
  print *, 'queried signal ', sigusr1
  print *, 'Installed action is: ', oldaction%sa_handler
  call raise(SIGUSR1)

  ! This is how a signal can be set to SIG_IGN or SIG_DFL
  action = sigaction_compile(SIG_IGN)
  call sigaction(SIGINT, action)

  ! C signal handlers should be installed with SA_NOWRAPPER
  action = sigaction_compile(c_sig_handler, flags = SA_NOWRAPPER)
  call sigaction(SIGUSR2, action)
  call raise(SIGUSR2)

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

end program f90_unix_signal_tst

subroutine sig_handler(sig)
  integer, intent(in) :: sig
  print *, 'recieved signal ', sig
end subroutine sig_handler
  
