!
! Test module f90_unix_io
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
program f90_unix_io_tst
  use f90_unix_io
  use f90_unix_tools, only : CRLF
  implicit none
  
  type(FILE) :: fp
  character(len=128) :: line
  integer :: linelen, errno !, nrec
  integer(sizet_kind) :: nrec

  ! Arrays for fread/fwrite
  character(len=10) :: strio(2) = (/ "0123456789", "9876543210" /)
  real(4) :: rio(2) = (/ 3.145, 2.718 /)
  real(8) :: dio(2) = (/ 3.145D0, 2.718D0 /)
  integer :: iio(2) = (/ 1, 2 /)

  fp = popen("ls", "r")
  if(.not. associated(fp)) then
     print *, "popen failed."
  end if

1000 continue
  call fgets(line, linelen, fp)
  if(linelen>0) print *, line(1:linelen-1)
  if(linelen>0) go to 1000

  call pclose(fp)

  fp = fopen("bla", "r", errno)
  if(.not. associated(fp)) then
     print *, "fopen failed as it should."
     call perror('fopen("bla", "r")');
  end if
  call check_errno(errno, 'fopen("bla", "r", errno)', 1)
  call set_errno
  
  fp = fopen("testfile", "w", errno)
  call check_errno(errno, "fopen")
  call fputs("Teststring"//CRLF, fp, errno)
  call check_errno(errno, "fputs")
  call fclose(fp, errno)
  call check_errno(errno, "fclose")

  call fputs("stdio-test"//CRLF, stdout())
  call fputs("Enter string:", stdout())
  call fgets(line, linelen, stdin())
  call fputs('>'//line(1:linelen)//'<', stderr())

  print '(A)', "Fread/fwrite tests:"
  fp = fopen("testfile", "w", errno)
  call check_errno(errno, "fopen")
  nrec = fwrite("0123456789", 10_sizet_kind, fp, errno) ;  if(nrec/=1) go to 8000
  call check_errno(errno, "fwrite-str")
  nrec = fwrite(strio, 10_sizet_kind, fp, errno) ;  if(nrec/=2) go to 8000
  call check_errno(errno, "fwrite-str-array")
  nrec = fwrite(1.0,  fp, errno) ;  if(nrec/=1) go to 8000
  call check_errno(errno, "fwrite-real")
  nrec = fwrite(rio,  fp, errno) ;  if(nrec/=2) go to 8000
  call check_errno(errno, "fwrite-real-array")
  nrec = fwrite(2.0D0,  fp, errno) ;  if(nrec/=1) go to 8000
  call check_errno(errno, "fwrite-double")
  nrec = fwrite(dio,  fp, errno) ;  if(nrec/=2) go to 8000
  call check_errno(errno, "fwrite-double-array")
  nrec = fwrite(4711,  fp, errno) ;  if(nrec/=1) go to 8000
  call check_errno(errno, "fwrite-integer")
  nrec = fwrite(iio,  fp, errno) ;  if(nrec/=2) go to 8000
  call check_errno(errno, "fwrite-integer-array")

  call fclose(fp)
  fp = fopen("testfile", "r", errno)
  
  strio(1)=''
  nrec = fread(strio(1), 10_sizet_kind, fp, errno) ; if(nrec/=1) go to 8010
  if(strio(1)/="0123456789") go to 8020
  call check_errno(errno, "fread-str")

  nrec = fread(strio, 10_sizet_kind, fp, errno) ; if(nrec/=2) go to 8010
  if(strio(1)/="0123456789") go to 8020
  if(strio(2)/="9876543210") go to 8020
  call check_errno(errno, "fread-str")
  
  nrec = fread(rio(1), fp, errno) ; if(nrec/=1) go to 8010
  if(rio(1)/=1.0) go to 8020
  call check_errno(errno, "fread-real")

  nrec = fread(rio, fp, errno) ; if(nrec/=2) go to 8010
  if(rio(1)/=3.145) go to 8020
  if(rio(2)/=2.718) go to 8020
  call check_errno(errno, "fread-real-array")
  
  nrec = fread(dio(1), fp, errno) ; if(nrec/=1) go to 8010
  if(dio(1)/=2.0D0) go to 8020
  call check_errno(errno, "fread-double")

  nrec = fread(dio, fp, errno) ; if(nrec/=2) go to 8010
  if(dio(1)/=3.145D0) go to 8020
  if(dio(2)/=2.718D0) go to 8020
  call check_errno(errno, "fread-double-array")

  nrec = fread(iio(1), fp, errno) ; if(nrec/=1) go to 8010
  if(iio(1)/=4711) go to 8020
  call check_errno(errno, "fread-integer")

  nrec = fread(iio, fp, errno) ; if(nrec/=2) go to 8010
  if(iio(1)/=1) go to 8020
  if(iio(2)/=2) go to 8020
  call check_errno(errno, "fread-integer-array")
  
  call fseek(fp, 20_long_kind, SEEK_SET)
  nrec = fread(strio(1), 10_sizet_kind, fp, errno) ; if(nrec/=1) go to 8010
  if(strio(1)/="9876543210") go to 8030
  call check_errno(errno, "fseek")

  print *, feof(fp)
  go to 9000

8000 continue
  print '(A,I3)', 'fwrite failed.', nrec
  stop

8010 continue
  print '(A, I3)', "fread failed.", nrec
  stop

8020 continue
  print '(A)', "fread /= fwrite."
  stop

8030 continue
  print '(A)', "fseek failed."
  stop

9000 continue
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
end program f90_unix_io_tst
