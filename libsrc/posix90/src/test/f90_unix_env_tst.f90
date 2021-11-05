!
! Test module f90_unix_env
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
program f90_unix_env_tst
  use f90_unix_env
  implicit none
  integer(LONG_KIND) :: value
  integer :: errno
  character(len=80) :: arg
  integer :: lenarg, ngroups
  integer, allocatable:: grouplist(:)
  integer(GID_KIND) :: gid
  integer(PID_KIND) :: pid
  integer(UID_KIND) :: uid
  type(tms) :: buffer
  type(utsname) :: name

  print *, "CLK_TCK:", clk_tck()

  print *, "ctermid:", ctermid()

  call getarg2(0, arg, lenarg, errno)
  print *, "getarg2((0, arg, lenarg, errno), arg=", arg(1:lenarg), &
       &" lenarg=", lenarg,  " errno=", errno

  print *, "getegid():",getegid()

  call getenv2("DISPLAY", arg, lenarg, errno)
  print *, "getenv2('DISPLAY', arg, lenarg, errno), arg=", arg(1:lenarg), &
       &" lenarg=", lenarg,  " errno=", errno

  print *, "geteuid():", geteuid()
  
  print *, "getgid():", getgid()

  call getgroups(ngroups=ngroups, errno = errno)
  print *, "getgroups(ngroups=ngroups, errno = errno); ngroups =", ngroups, errno
  allocate(grouplist(ngroups))
  call getgroups(grouplist, ngroups, errno)
  print *, "getgroups(grouplist, ngroups, errno); grouplist(:)=", grouplist(:)

  call gethostname(arg, lenarg, errno)
  print *, "gethostname(arg, lenarg, errno), arg=",trim(arg), " lenarg=", &
       &lenarg, "errno=", errno

  call getlogin(arg, lenarg, errno)
  print *, "getlogin(arg, lenarg, errno), arg=",trim(arg), " lenarg=", &
       &lenarg, "errno=", errno

  print *, "getpgrp():", getpgrp()

  print *, "getpid():", getpid()

  print *, "getppid():", getppid()

  print *, "getuid():", getuid()

  gid = -1
  call setgid(gid, errno)
  print *, 'setgid(gid, errno), gid=',gid, ' errno=', errno, &
       &trim(strerror(errno))
  call set_errno(0)

  gid = -1
  call setpgid(gid, gid, errno)
  print *, 'setpgid(gid, gid, errno), gid=',gid, ' errno=', errno, &
       &trim(strerror(errno))
  call set_errno(0)

  call setsid(errno)
  print *, 'setsid(errno), errno=', errno, &
       &trim(strerror(errno))
  call set_errno(0)

  uid = -1
  call setuid(uid, errno)
  print *, 'setuid(uid, errno), uid=',uid, ' errno=', errno, &
       &trim(strerror(errno))
  call set_errno(0)

  call sysconf(SC_VERSION, value, errno)
  print *, "sysconf(SC_VERSION, value, errno), value = ", &
       &value, " errno=", errno

  call sysconf(-1, value, errno)
  print *, "sysconf(-1, value, errno), value = ", value, " errno=", errno
  call set_errno (0);

  print *, "time(errno)= ", time(errno), " errno=", errno

  print *, "times(buffer,errno)= ", times(buffer, errno), " errno=", errno
  print *, "buffer%utime: ", buffer%utime
  print *, "buffer%stime: ", buffer%stime
  print *, "buffer%cutime: ", buffer%cutime
  print *, "buffer%cstime: ", buffer%cstime

  call uname(name, errno)
  print *, "uname(name, errno)", " errno=", errno
  print *, "name%sysname: ", trim(name%sysname)
  print *, "name%nodename: ", trim(name%nodename)
  print *, "name%release: ", trim(name%release) 
  print *, "name%version: ", trim(name%version)
  print *, "name%machine: ", trim(name%machine)

end program f90_unix_env_tst 
