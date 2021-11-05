PROGRAM regexp

#if defined(HAVE_POSIX90)
  USE messy_main_tools, ONLY: match_regexp

  INTEGER           :: status
  CHARACTER(LEN=80) :: pattern = ''
  CHARACTER(LEN=80) :: string  = ''

  WRITE(*,*) 'PRESS CTRL-C TO STOP.'
  WRITE(*,*)

  WRITE(*,*) 'STRING ?:'
  READ(*,*) string

  DO
     WRITE(*,*) 'PATTERN ?:'
     READ(*,*) pattern

     WRITE(*,*) match_regexp(status, pattern, string)
     WRITE(*,*) "status = ", status
  END DO
  
#else
  WRITE(*,*) 'ERROR: regexp.exe requires POSIX90 library.'
  WRITE(*,*) '       reconfigure with --enable-POSIX90 and recompile.'
#endif
  
END PROGRAM regexp
