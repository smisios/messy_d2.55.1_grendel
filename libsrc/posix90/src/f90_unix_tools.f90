!
! utility routines and constants
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
module f90_unix_tools
  implicit none
  character(len=1), parameter :: C0=char(0) 
  character(len=1), parameter :: CR=char(13)
  character(len=1), parameter :: LF=char(10)
  character(len=2), parameter :: CRLF=CR//LF
contains
 subroutine fortranify_string(str, len)
   character(len=*), intent(inout) :: str
   integer, intent(out), optional :: len
   integer :: l

   l = index(str, char(0))
   if(present(len)) len = l-1
   if(l>0) str(l:)=' '
 end subroutine fortranify_string
 
end module f90_unix_tools
