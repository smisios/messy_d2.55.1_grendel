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
module f90_unix_io
  use f90_unix_errno
  use f90_unix_env, only : sizet_kind, NULL
  use f90_unix_tools, only : C0
  implicit none
  
  include 'f90_unix_io_const.inc'
  
  type FILE
     integer(file_kind) :: fp
  end type FILE
  
  interface fread
     module procedure fread_str, fread_str_array
     module procedure fread_real, fread_real_array
     module procedure fread_double, fread_double_array
     module procedure fread_int, fread_int_array
  end interface
  
  interface fwrite
     module procedure fwrite_str, fwrite_str_array
     module procedure fwrite_real, fwrite_real_array
     module procedure fwrite_double, fwrite_double_array
     module procedure fwrite_int, fwrite_int_array
  end interface
  
  interface associated
     module procedure fassociated
  end interface
  
  private manage_errno
  integer, parameter :: double = kind(1.D0)

contains

  logical function feof(stream, errno)
    type(FILE), intent(in) :: stream
    integer, intent(out), optional :: errno
    logical, external :: c_feof

    if(present(errno)) errno = 0
    call set_errno

    feof = c_feof(stream)
    call manage_errno(errno,"f90_unix_io::feof")
  end function feof
  
  subroutine rewind(stream, errno)
    type(FILE), intent(inout) :: stream
    integer, intent(out), optional :: errno

    call fseek(stream, 0_long_kind, SEEK_SET, errno)
  end subroutine rewind
  
  subroutine fseek(stream, offset, whence, errno)
    type(FILE), intent(inout) :: stream
    integer(long_kind), intent(in) :: offset
    integer, intent(in), optional :: whence
    integer, intent(out), optional :: errno
    integer:: w

    if(present(errno)) errno = 0
    call set_errno

    w=SEEK_SET
    if(present(whence)) w = whence

    call c_fseek(stream, offset, w)
    call manage_errno(errno, "f90_unix_io::fseek")
  end subroutine fseek
  
  integer(long_kind) function ftell(stream, errno)
    type(FILE), intent(inout) :: stream
    integer, intent(out), optional :: errno
    integer(long_kind), external :: c_ftell

    if(present(errno)) errno = 0
    call set_errno

    ftell = c_ftell(stream)
    call manage_errno(errno, "f90_unix_io::ftell")
  end function ftell
  
  logical function fassociated(fp)
    type(FILE), intent(in) :: fp
    fassociated = (fp%fp /= NULL)
  end function fassociated
  
  type(FILE) function fopen(path, mode, errno)
    character(len=*), intent(in) :: path, mode
    integer, intent(out), optional :: errno
    type(FILE) :: fp

    if(present(errno)) errno = 0
    call set_errno

    call c_fopen(fp%fp, path//C0, mode//C0)
    fopen%fp = fp%fp

    call manage_errno(errno, &
         &"f90_unix_io::fopen('"//trim(path)//"','"//trim(mode)//"')")

  end function fopen

  subroutine fclose(fp, errno)
    type(FILE), intent(inout) :: fp
    integer, intent(out), optional :: errno

    if(present(errno)) errno = 0
    call set_errno

    call c_fclose(fp%fp)

    call manage_errno(errno, "f90_unix_io::fclose")
  end subroutine fclose
  
  type(FILE) function popen(command, mode, errno)
    character(len=*), intent(in) :: command, mode
    integer, intent(out), optional :: errno
    type(FILE) :: fp

    if(present(errno)) errno = 0
    call set_errno

    call c_popen(fp%fp,command//C0, mode//C0)
    popen = fp

    call manage_errno(errno, &
         &"f90_unix_io::popen('"//trim(command)//"','"//trim(mode)//"')")
  end function popen

  subroutine pclose(fp, errno)
    type(FILE), intent(inout) :: fp
    integer, intent(out), optional :: errno

    if(present(errno)) errno = 0
    call set_errno

    call c_pclose(fp%fp)

    call manage_errno(errno, "f90_unix_io::pclose")

  end subroutine pclose
  
  integer(sizet_kind) function fread_str(str, length, fp, errno)
    character(len=*), intent(out) :: str
    integer(sizet_kind), intent(in) :: length
    type(FILE), intent(inout) :: fp
    integer, intent(out), optional :: errno
    integer(sizet_kind), external :: c_fread_str

    if(present(errno)) errno = 0
    call set_errno

    if(length>len(str)) then
       call set_errno(EINVAL) ; go to 8000
    end if
    
    fread_str = c_fread_str(str, length, fp%fp)

8000 continue
    call manage_errno(errno, "f90_unix_io::fread_str")

  end function fread_str

  integer(sizet_kind) function fread_str_array(str, length, fp, errno)
    character(len=*), intent(out) :: str(:)
    integer(sizet_kind), intent(in) :: length
    type(FILE), intent(inout) :: fp
    integer, intent(out), optional :: errno
    integer, external :: c_fread_str_array

    if(present(errno)) errno = 0
    call set_errno

    if(length>len(str(1))) then
       call set_errno(EINVAL) ; go to 8000
    end if

    fread_str_array = c_fread_str_array(str, length, size(str,1), fp%fp) 

8000 continue
    call manage_errno(errno, "f90_unix_io::fread_str_array")
    
  end function fread_str_array

  integer(sizet_kind) function fread_real(r, fp, errno)
    real, intent(out) :: r
    type(FILE), intent(inout) :: fp
    integer, intent(out), optional :: errno
    integer(sizet_kind), external :: c_fread_real

    if(present(errno)) errno = 0
    call set_errno

    fread_real = c_fread_real(r, 1, fp%fp)

    call manage_errno(errno, "f90_unix_io::fread_real")
  end function fread_real

  integer(sizet_kind) function fread_real_array(r, fp, errno)
    real, intent(out) :: r(:)
    type(FILE), intent(inout) :: fp
    integer, intent(out), optional :: errno
    integer, external :: c_fread_real

    if(present(errno)) errno = 0
    call set_errno

    fread_real_array = c_fread_real(r, size(r,1), fp%fp) 

    call manage_errno(errno, "f90_unix_io::fread_real_array")
  end function fread_real_array
  
  integer(sizet_kind) function fread_double(d, fp, errno)
    real(double), intent(out) :: d
    type(FILE), intent(inout) :: fp
    integer, intent(out), optional :: errno
    integer(sizet_kind), external :: c_fread_double

    if(present(errno)) errno = 0
    call set_errno

    fread_double = c_fread_double(d, 1, fp%fp)

    call manage_errno(errno, "f90_unix_io::fread_double")
  end function fread_double

  integer(sizet_kind) function fread_double_array(d, fp, errno)
    real(double), intent(out) :: d(:)
    type(FILE), intent(inout) :: fp
    integer, intent(out), optional :: errno
    integer, external :: c_fread_double

    if(present(errno)) errno = 0
    call set_errno

    fread_double_array = c_fread_double(d, size(d,1), fp%fp) 

    call manage_errno(errno, "f90_unix_io::fread_double_array")
  end function fread_double_array
  
  integer(sizet_kind) function fread_int(i, fp, errno)
    integer, intent(out) :: i
    type(FILE), intent(inout) :: fp
    integer, intent(out), optional :: errno
    integer(sizet_kind), external :: c_fread_int

    if(present(errno)) errno = 0
    call set_errno

    fread_int = c_fread_int(i, 1, fp%fp)

    call manage_errno(errno, "f90_unix_io::fread_int")
  end function fread_int

  integer(sizet_kind) function fread_int_array(i, fp, errno)
    integer, intent(out) :: i(:)
    type(FILE), intent(inout) :: fp
    integer, intent(out), optional :: errno
    integer, external :: c_fread_int

    if(present(errno)) errno = 0
    call set_errno

    fread_int_array = c_fread_int(i, size(i,1), fp%fp) 

    call manage_errno(errno, "f90_unix_io::fread_int_array")
  end function fread_int_array
  
  integer(sizet_kind) function fwrite_str(str, length, fp, errno)
    character(len=*), intent(in) :: str
    integer(sizet_kind), intent(in) :: length
    type(FILE), intent(inout) :: fp
    integer, intent(out), optional :: errno
    integer, external :: c_fwrite_str
 
    if(present(errno)) errno = 0
    call set_errno

    if(length>len(str)) then
       call set_errno(EINVAL) ; go to 8000
    end if
    
    fwrite_str = c_fwrite_str(str, length, fp%fp)

8000 continue
    call manage_errno(errno,"f90_unix_io::fwrite_str")

    
  end function fwrite_str

  integer(sizet_kind) function fwrite_str_array(str, length, fp, errno)
    character(len=*), intent(in) :: str(:)
    integer(sizet_kind), intent(in) :: length
    type(FILE), intent(inout) :: fp
    integer, intent(out), optional :: errno
    integer, external :: c_fwrite_str_array

    if(present(errno)) errno = 0
    call set_errno

    if(length>len(str(1))) then
       call set_errno(EINVAL) ; go to 8000
    end if
    
    fwrite_str_array = c_fwrite_str_array(str, length, size(str,1), fp%fp) 

8000 continue
    call manage_errno(errno,"f90_unix_io::fwrite_str_array")

  end function fwrite_str_array
  
  integer(sizet_kind) function fwrite_real(r, fp, errno)
    real, intent(in) :: r
    type(FILE), intent(inout) :: fp
    integer, intent(out), optional :: errno
    integer(sizet_kind), external :: c_fwrite_real

    if(present(errno)) errno = 0
    call set_errno

    fwrite_real = c_fwrite_real(r, 1, fp%fp)

    call manage_errno(errno, "f90_unix_io::fwrite_real")
  end function fwrite_real

  integer(sizet_kind) function fwrite_real_array(r, fp, errno)
    real, intent(in) :: r(:)
    type(FILE), intent(inout) :: fp
    integer, intent(out), optional :: errno
    integer, external :: c_fwrite_real

    if(present(errno)) errno = 0
    call set_errno

    fwrite_real_array = c_fwrite_real(r, size(r,1), fp%fp) 

    call manage_errno(errno, "f90_unix_io::fwrite_real_array")
  end function fwrite_real_array
  
  integer(sizet_kind) function fwrite_double(r, fp, errno)
    real(double), intent(in) :: r
    type(FILE), intent(inout) :: fp
    integer, intent(out), optional :: errno
    integer(sizet_kind), external :: c_fwrite_double

    if(present(errno)) errno = 0
    call set_errno

    fwrite_double = c_fwrite_double(r, 1, fp%fp)

    call manage_errno(errno, "f90_unix_io::fwrite_double")
  end function fwrite_double

  integer(sizet_kind) function fwrite_double_array(r, fp, errno)
    real(double), intent(in) :: r(:)
    type(FILE), intent(inout) :: fp
    integer, intent(out), optional :: errno
    integer, external :: c_fwrite_double

    if(present(errno)) errno = 0
    call set_errno

    fwrite_double_array = c_fwrite_double(r, size(r,1), fp%fp) 

    call manage_errno(errno, "f90_unix_io::fwrite_double_array")
  end function fwrite_double_array
  
  integer(sizet_kind) function fwrite_int(i, fp, errno)
    integer, intent(in) :: i
    type(FILE), intent(inout) :: fp
    integer, intent(out), optional :: errno
    integer(sizet_kind), external :: c_fwrite_int

    if(present(errno)) errno = 0
    call set_errno

    fwrite_int = c_fwrite_int(i, 1, fp%fp)

    call manage_errno(errno, "f90_unix_io::fwrite_int")
  end function fwrite_int

  integer(sizet_kind) function fwrite_int_array(i, fp, errno)
    integer, intent(in) :: i(:)
    type(FILE), intent(inout) :: fp
    integer, intent(out), optional :: errno
    integer, external :: c_fwrite_int

    if(present(errno)) errno = 0
    call set_errno

    fwrite_int_array = c_fwrite_int(i, size(i,1), fp%fp) 

    call manage_errno(errno, "f90_unix_io::fwrite_int_array")
  end function fwrite_int_array
  
  subroutine fgets(str, strlen, fp, errno)
    character(len=*), intent(inout) :: str
    integer, intent(out) :: strlen
    type(FILE), intent(in) :: fp
    integer, intent(out), optional :: errno

    if(present(errno)) errno = 0
    call set_errno

    call c_fgets(str,fp%fp)
    call fortranify_string(str, strlen)
    if(present(errno)) then
       errno = get_errno()
    else
       if(get_errno() /= 0) then
          call perror(&
               &"f90_unix_io::fgets")
          stop
       end if
    endif
    
  end subroutine fgets
  
  subroutine fputs(str, fp, errno)
    character(len=*), intent(in) :: str
    type(FILE):: fp
    integer, intent(out), optional :: errno

    if(present(errno)) errno = 0
    call set_errno

    call c_fputs(str//C0, fp%fp)
    if(present(errno)) then
       errno = get_errno()
    else
       if(get_errno() /= 0) then
          call perror(&
               &"f90_unix_io::fputs")
          stop
       end if
    endif
  end subroutine fputs
  
  type(FILE) function stdin()
    call c_stdin(stdin%fp)
  end function stdin
  
  type(FILE) function stdout()
    call c_stdout(stdout%fp)
  end function stdout
  
  type(FILE) function stderr()
    call c_stderr(stderr%fp)
  end function stderr

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
  
end module f90_unix_io
