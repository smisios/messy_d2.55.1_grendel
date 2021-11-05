!
! Interface to the facilities from ISO/IEC 9945-1:1990 sec. 4
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
module f90_unix_env
  use f90_unix_errno
  use f90_unix_time, only : time_kind
  implicit none
  include 'f90_unix_env_const.inc'

  type tms
     sequence
     integer(clock_kind):: utime, stime,cutime,cstime
  end type tms
  
  type utsname
     sequence
     character(len=80):: sysname, nodename, release, version, machine
  end type utsname
  

contains

  integer(kind=clock_kind) function clk_tck()
    integer(kind=clock_kind), external:: c_clk_tck
    clk_tck=c_clk_tck()
  end function clk_tck
  
  character(len=L_CTERMID) function ctermid(len)
    integer, intent(out), optional :: len
    character(len=L_CTERMID) :: s
    call c_ctermid(s)
    if(present(len)) len = len_trim(s)
    ctermid = s
  end function ctermid
  
  subroutine getarg2(k, arg, lenarg, errno)
    integer, intent(in) :: K
    character(len=*), intent(out), optional :: arg
    integer, intent(out), optional :: lenarg, errno

    if(present(errno) ) errno = 0
    call set_errno
    
    if( (k<0) .or. (k>iargc())) go to 8000

    call getarg(k,arg)

    if(present(lenarg)) lenarg = len_trim(arg)
    go to 9000

8000 continue
    if(present(errno)) then
       errno = EINVAL
    else
       call set_errno(EINVAL)
       call perror("f90_unix_env::getarg2");
       stop
    end if
    
9000 continue
  end subroutine getarg2
  
  integer(GID_KIND) function getegid()
    integer(GID_KIND), external :: c_getegid
    getegid = c_getegid()
  end function getegid
  
  subroutine getenv2(name, value, lenvalue, errno)
    character(len=*), intent(in) :: name
    character(len=*), intent(out), optional :: value
    integer, intent(out), optional :: lenvalue, errno

    if(present(errno) ) errno = 0
    call set_errno
    
    call getenv(name,value)

    if(present(lenvalue)) lenvalue = len_trim(value)
    go to 9000

    if(present(errno)) then
       errno = EINVAL
    else
       call set_errno(EINVAL)
       call perror("f90_unix_env::getenv2");
       stop
    end if
    
9000 continue
  end subroutine getenv2

  integer(uid_kind) function geteuid()
    integer(uid_kind),external::c_geteuid
    geteuid = c_geteuid()
  end function geteuid
  
  integer(gid_kind) function getgid()
    integer(gid_kind),external::c_getgid
    getgid= c_getgid()
  end function getgid
  
  subroutine getgroups(grouplist, ngroups, errno)
    integer(gid_kind), optional :: grouplist(:)
    integer, optional, intent(out) :: ngroups, errno
    integer :: ngr, err
    external c_getgroups

    if(present(errno)) errno = 0
    call set_errno

    if(present(grouplist)) then
       call c_getgroups(grouplist, size(grouplist,1), ngr, err)
    else
       call c_getgroups(grouplist, 0, ngr, err)
    end if
   
    if(err /=0) go to 8000
    if(present(ngroups)) ngroups = ngr
    go to 9000

8000 continue
    
    if(present(errno)) then
       errno = err
    else
       call perror("f90_unix_env::GetGroups:", err)
       stop
    end if
    
9000 continue
  end  subroutine getgroups

  subroutine gethostname(name, lenname, errno)
    character(len=*), optional, intent(out) :: name
    integer, optional, intent(out) :: lenname, errno
    
    call c_gethostname(name)
    call fortranify_string(name, lenname)

    if(present(errno)) then
       errno = get_errno()
    else
       if(get_errno()/=0) then
          call perror("f90_unix_env::Gethostname:" );
          stop
       end if
    end if
    
  end subroutine gethostname


  subroutine getlogin(name, lenname, errno)
    character(len=*), optional, intent(out) :: name
    integer, optional, intent(out) :: lenname, errno
    
    if(present(errno)) errno = 0
    call set_errno

    call c_getlogin(name)
    call fortranify_string(name, lenname)

    if(present(errno)) then
       errno = get_errno()
    else
       if(get_errno()/=0) then
          call perror("f90_unix_env::Getlogin:" );
          stop
       end if
    end if
    
  end subroutine getlogin

  integer(PID_KIND) function getpgrp()
    integer(PID_KIND),external:: c_getpgrp
    getpgrp = c_getpgrp()
  end function getpgrp
  
  integer(PID_KIND) function getpid()
    integer(PID_KIND),external:: c_getpid
    getpid = c_getpid()
  end function getpid
  
  integer(PID_KIND) function getppid()
    integer(PID_KIND),external:: c_getppid
    getppid = c_getppid()
  end function getppid
  
  integer(UID_KIND) function getuid()
    integer(UID_KIND),external:: c_getuid
    getuid = c_getuid()
  end function getuid
  
  subroutine setgid(gid, errno)
    integer(GID_KIND), intent(in) :: gid
    integer, intent(out), optional :: errno
    integer :: errc

    if(present(errno)) errno = 0
    call set_errno

    call c_setgid(gid)

    errc= get_errno()
    if(errc ==0) go to 9000

    if(present(errno)) then
       errno=errc
    else
       call perror("f90_unix_env::setgid")
       stop
    end if
  
9000 continue  
  end subroutine setgid
    

  subroutine setpgid(gid, pgid, errno)
    integer(GID_KIND), intent(in) :: gid, pgid
    integer, intent(out), optional :: errno
    integer :: errc

    if(present(errno)) errno = 0
    call set_errno

    call c_setpgid(gid,pgid)

    errc= get_errno()
    if(errc ==0) go to 9000

    if(present(errno)) then
       errno=errc
    else
       call perror("f90_unix_env::setpgid")
       stop
    end if
  
9000 continue  
  end subroutine setpgid
  

  subroutine setsid(errno)
    integer, intent(out), optional :: errno
    integer :: errc

    if(present(errno)) errno = 0
    call set_errno

    call c_setsid()

    errc= get_errno()
    if(errc ==0) go to 9000

    if(present(errno)) then
       errno=errc
    else
       call perror("f90_unix_env::setsid")
       stop
    end if
  
9000 continue  
  end subroutine setsid
  

  subroutine setuid(gid, errno)
    integer(UID_KIND), intent(in) :: gid
    integer, intent(out), optional :: errno
    integer :: errc

    if(present(errno)) errno = 0
    call set_errno

    call c_setuid(gid)

    errc= get_errno()
    if(errc ==0) go to 9000

    if(present(errno)) then
       errno=errc
    else
       call perror("f90_unix_env::setuid")
       stop
    end if
  
9000 continue  
  end subroutine setuid
  
  subroutine sysconf(name, val,errno)
    integer, intent(in) :: name
    integer(long_kind), intent(out) :: val
    integer, intent(out), optional :: errno
    integer :: err

    if(present(errno)) errno = 0
    call set_errno

    call c_sysconf(name,val,err)
    if(err/=0) go to 8000
    go to 9000

8000 continue
    if(present(errno)) then
       errno = err
    else
       call perror("f90_unix_env::sysconf");
       stop
    end if
    
9000 continue
  end subroutine sysconf
  

  integer(TIME_KIND) function time(errno)
    integer, optional, intent(out) :: errno
    integer(TIME_KIND), external :: c_time
 
    if(present(errno)) errno = 0
    call set_errno

    time = c_time()
    if(get_errno()/=0) go to 8000
    go to 9000

8000 continue
    if(present(errno)) then
       errno = get_errno()
    else
       call perror("f90_unix_env::time");
       stop
    end if
    
9000 continue
    
  end function time

  
  integer(CLOCK_KIND) function times(buffer, errno)
    type(tms) :: buffer
    integer, optional, intent(out) :: errno
    integer(clock_kind), external :: c_times

    if(present(errno)) errno = 0
    call set_errno

    times = c_times(buffer%utime, buffer%stime, buffer%cutime, buffer%cstime)

    if(get_errno()/=0) go to 8000
    go to 9000

8000 continue
    if(present(errno)) then
       errno = get_errno()
    else
       call perror("f90_unix_env::times");
       stop
    end if
    
9000 continue
    
  end function times
  
  subroutine uname(name, errno)
    type(utsname), intent(out) :: name
    integer, optional, intent(out) :: errno

    if(present(errno)) errno = 0
    call set_errno

    call c_uname(name%sysname, name%nodename, name%release, &
         &       name%version, name%machine)

    call fortranify_string(name%sysname)
    call fortranify_string(name%nodename)
    call fortranify_string(name%release)
    call fortranify_string(name%version)
    call fortranify_string(name%machine)

    if(get_errno()/=0) go to 8000
    go to 9000

8000 continue
    if(present(errno)) then
       errno = get_errno()
    else
       call perror("f90_unix_env::uname");
       stop
    end if
    
9000 continue
       
  end subroutine uname
  
end module f90_unix_env

