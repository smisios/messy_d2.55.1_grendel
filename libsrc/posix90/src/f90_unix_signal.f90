!
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
module f90_unix_signal
  use f90_unix_errno
  use f90_unix_tools, only : C0
  use f90_unix_env, only : NULL, pid_kind
  implicit none
  include 'f90_unix_signal_const.inc'

  type sigaction_type
     integer(funcp_kind) :: sa_handler
     type(sigset_type) :: sa_mask
     integer :: sa_flags 
  end type sigaction_type
  
  private manage_errno

  interface sigaction_compile
     module procedure sigaction_compile_handler
     module procedure sigaction_compile_integer
  end interface
  
contains

  type(sigaction_type) function sigaction_compile_handler(handler, mask, flags, errno)
    interface
       subroutine handler(sig)
         integer, intent(in) :: sig
       end subroutine handler       
    end interface
    type(sigset_type), intent(in), optional :: mask
    integer, intent(in), optional :: flags
    integer, intent(out), optional :: errno
    type(sigset_type) :: lmask
    integer :: lflags

    if(present(errno)) errno = 0
    lflags = 0 ; if (present(flags)) lflags = flags
    call sigemptyset(lmask) ; if(present(mask)) lmask = mask

    call c_sigaction_compile(sigaction_compile_handler%sa_handler, &
         &sigaction_compile_handler%sa_mask%sigset,&
         &sigaction_compile_handler%sa_flags, handler, lmask%sigset, lflags)
  end function sigaction_compile_handler
  
  type(sigaction_type) function sigaction_compile_integer(action_code, mask, flags, errno)
    integer, intent(in) :: action_code
    type(sigset_type), intent(in), optional :: mask
    integer, intent(in), optional :: flags
    integer, intent(out), optional :: errno
    type(sigset_type) :: lmask
    integer :: lflags

    if(present(errno)) errno = 0
    lflags = 0 ; if (present(flags)) lflags = flags
    call sigemptyset(lmask) ; if(present(mask)) lmask = mask

    if( (action_code/=SIG_IGN) .and. (action_code/=SIG_DFL)) go to 8000
    call c_sigaction_compile(sigaction_compile_integer%sa_handler, &
         &sigaction_compile_integer%sa_mask%sigset,&
         &sigaction_compile_integer%sa_flags, int(action_code, funcp_kind),&
         & lmask%sigset, lflags)
    go to 9000

8000 continue
    ! action_code is neither SIG_IGN nor SIG_DFL. This is illegal
    call set_errno(EINVAL)
    call manage_errno(errno, "f90_unix_signal::sigaction_compile_integer")

9000 continue

  end function sigaction_compile_integer
  
  subroutine sigaction(signum, action, oldaction, errno)
    integer :: signum
    type(sigaction_type), intent(in), optional :: action
    type(sigaction_type), intent(out), optional :: oldaction
    integer, intent(out), optional :: errno
    type(sigaction_type):: loldaction
    
    if(present(errno)) errno = 0
    call set_errno
   
    if(present(action)) then
       call c_sigaction(signum,&
            & action%sa_handler, action%sa_mask%sigset,action%sa_flags,&
            & loldaction%sa_handler,loldaction%sa_mask%sigset,loldaction%sa_flags)
    else
       call c_sigaction(signum,&
            & NULL, '', NULL,&
            & loldaction%sa_handler,loldaction%sa_mask%sigset,loldaction%sa_flags)
    end if
    
    if(present(oldaction)) oldaction = loldaction

    call manage_errno(errno, "f90_unix_signal::sigaction")
   
  end subroutine sigaction
  
  subroutine sigemptyset(set)
    type(sigset_type), intent(in):: set
    call c_sigemptyset(set%sigset)
  end subroutine sigemptyset
  
  subroutine sigfillset(set)
    type(sigset_type), intent(in):: set
    call c_sigfillset(set%sigset)
  end subroutine sigfillset

  subroutine sigaddset(set, sig1, sig2, sig3, sig4, sig5, sig6, sig7, sig8)
    type(sigset_type), intent(in):: set
    integer, intent(in) :: sig1
    integer, intent(in), optional :: sig2, sig3, sig4, sig5, sig6, sig7, sig8
    call c_sigaddset(set%sigset,sig1)
    if(present(sig2)) call c_sigaddset(set%sigset,sig2)
    if(present(sig3)) call c_sigaddset(set%sigset,sig3)
    if(present(sig4)) call c_sigaddset(set%sigset,sig4)
    if(present(sig5)) call c_sigaddset(set%sigset,sig5)
    if(present(sig6)) call c_sigaddset(set%sigset,sig6)
    if(present(sig7)) call c_sigaddset(set%sigset,sig7)
    if(present(sig8)) call c_sigaddset(set%sigset,sig8)
  end subroutine sigaddset
  
  subroutine sigdelset(set, sig1, sig2, sig3, sig4, sig5, sig6, sig7, sig8)
    type(sigset_type), intent(in):: set
    integer, intent(in) :: sig1
    integer, intent(in), optional :: sig2, sig3, sig4, sig5, sig6, sig7, sig8
    call c_sigdelset(set%sigset,sig1)
    if(present(sig2)) call c_sigdelset(set%sigset,sig2)
    if(present(sig3)) call c_sigdelset(set%sigset,sig3)
    if(present(sig4)) call c_sigdelset(set%sigset,sig4)
    if(present(sig5)) call c_sigdelset(set%sigset,sig5)
    if(present(sig6)) call c_sigdelset(set%sigset,sig6)
    if(present(sig7)) call c_sigdelset(set%sigset,sig7)
    if(present(sig8)) call c_sigdelset(set%sigset,sig8)
  end subroutine sigdelset
  
  logical function sigismember(set, sig)
    type(sigset_type), intent(in):: set
    integer, intent(in) :: sig
    logical, external :: c_sigismember
    sigismember = c_sigismember(set%sigset, sig)
  end function sigismember

  subroutine kill(pid, sig, errno)
    integer(pid_kind), intent(in) :: pid
    integer, intent(in) :: sig
    integer, intent(out), optional :: errno
    if(present(errno)) errno = 0
    call set_errno
    call c_kill(pid, sig)
    call manage_errno(errno, "f90_unix_signal::kill")
  end subroutine kill
  
  subroutine raise(sig,errno)
    integer, intent(in) :: sig
    integer, intent(out), optional :: errno
    if(present(errno)) errno = 0
    call set_errno
    call c_raise(sig)
    call manage_errno(errno, "f90_unix_signal::raise")
 end subroutine raise
  
  subroutine manage_errno(errno, loc)
    integer, intent(OUT), optional :: errno
    character(len=*), intent(in) :: loc
    
    if(present(errno)) then
       errno = get_errno()
    else
       if(get_errno() /= 0) then
          call perror(loc)
          stop
       end if
    endif
  end subroutine manage_errno
  
end module f90_unix_signal

