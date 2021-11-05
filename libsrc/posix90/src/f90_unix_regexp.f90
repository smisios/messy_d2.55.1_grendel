!
! Interface to the regexp facilities from ISO/IEC 9945-1:1990 
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
module f90_unix_regexp
  use f90_unix_errno
  use f90_unix_tools, only : C0, fortranify_string
  implicit none
  include 'f90_unix_regexp_const.inc'

  type regex_t
     integer(regex_kind) :: rp
  end type regex_t
  
  type regmatch_t
     integer(regoff_kind) :: rm_so, rm_eo
  end type regmatch_t
  
contains

  subroutine regcomp(preg, regex, cflags, errc)
    type(regex_t) :: preg
    character(len=*), intent(in) :: regex
    integer, intent(in) :: cflags
    integer, intent(out), optional :: errc
    integer :: e
    character(len=80) :: msg

    if(present(errc)) errc = 0

    call c_regcomp(preg%rp, regex//C0, cflags,e)
    if(e/=0) then
       if(present(errc)) then
          errc = e
       else
          call regerror(e, preg, msg)
          call fortranify_string(msg)
          print '(A)', trim(msg)
          stop
       end if
    end if    
  end subroutine regcomp
  
  subroutine regexec(preg, string, pmatch, eflags, errc)
    type(regex_t), intent(in) :: preg
    character(len=*), intent(in) :: string
    type(regmatch_t) :: pmatch(:)
    integer, intent(in) :: eflags
    integer, intent(out), optional :: errc
    integer :: e
    character(len=80) :: msg
    integer :: i

    if(present(errc)) errc = 0

    call c_regexec(preg%rp, string//C0, size(pmatch,1), pmatch(1), eflags, e)
    
    do i=1, size(pmatch,1)
        if(pmatch(i)%rm_so>=0) pmatch(i)%rm_so =  pmatch(i)%rm_so + 1
    end do
    
    if(e/=0) then
       if(present(errc)) then
          errc = e
       else
          call regerror(e, preg, msg)
          call fortranify_string(msg)
          print '(A)', trim(msg)
          stop
       end if
    end if    
  end subroutine regexec

  subroutine regerror(e, preg, msg)
    integer, intent(in) :: e
    type(regex_t), intent(in) :: preg
    character(len=*), intent(out) :: msg
   
     call c_regerror(e, preg%rp, msg)
     call fortranify_string(msg)
  end subroutine regerror
  
  subroutine regfree(preg)
    type(regex_t), intent(in) :: preg
    call c_regfree(preg%rp)
  end subroutine regfree
  
end module f90_unix_regexp

