!
! Test module f90_unix_dirent
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
program f90_unix_dirent_tst
  use f90_unix_dirent
  implicit none

  type(dir) :: dirp
  integer :: errno
  character(len=128) :: name
  integer :: lenname
  
  call opendir('.', dirp, errno)
  call rewinddir(dirp, errno)

  do 
     call readdir(dirp, name, lenname, errno)
     if(lenname<0) go to 1000
     print *, name(1:lenname)
  end do

1000 continue  
  call closedir(dirp, errno)

end program f90_unix_dirent_tst
