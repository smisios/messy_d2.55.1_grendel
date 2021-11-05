!
! Interface to the facilities from ISO/IEC 9945-1:1990 sec.
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
module f90_unix_time
  use f90_unix_errno
  use f90_unix_tools, only : C0, fortranify_string
  implicit none
  include 'f90_unix_time_const.inc'

contains

  character(len=30) function ctime(time)
    integer(time_kind), intent(in) :: time

    call c_ctime(ctime, time)
    call fortranify_string(ctime)
  end function ctime
  
  
end module f90_unix_time

