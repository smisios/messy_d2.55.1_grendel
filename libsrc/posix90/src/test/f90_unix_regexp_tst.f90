!
! Test module f90_unix_regexp
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
program f90_unix_regexp_tst
  use f90_unix_regexp
  implicit none

  character(len=80), parameter :: re = "([a-z]*)([0-9]*)([a-z]*)"
  character(len=80), parameter :: string = "abcd123efg"
  type(regex_t) :: preg
  type(regmatch_t) :: pmatch(5)
  integer :: i,j1,j2,err

  call regcomp(preg, trim(re), REG_EXTENDED, err)
  call check_errc(err, preg, "regcomp "//trim(re))

  call regexec(preg, trim(string), pmatch, 0, err)
  call check_errc(err, preg, trim(re)//"~="//trim(string))

  if(err==0) then
     do i=1,size(pmatch,1)
        j1 = pmatch(i)%rm_so
        j2 = pmatch(i)%rm_eo
        if(j1>=0) print *, j1,j2,string(j1:j2)
     end do
  end if
  
  call regfree(preg)
contains

  subroutine check_errc(errc, preg, msg, flag)
    integer, intent(in) :: errc
    type(regex_t) :: preg
    character(len=*), intent(in) :: msg
    integer,intent(in),optional:: flag
    character(len=128) :: emsg

    if(errc/=0) then
       call regerror(errc, preg, emsg)
       if(present(flag)) then
        print *, msg//' failed as it should:'//trim(emsg)
       else
        print *, msg//' failed, but should not:'//trim(emsg)
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
  end subroutine check_errc

end program f90_unix_regexp_tst
