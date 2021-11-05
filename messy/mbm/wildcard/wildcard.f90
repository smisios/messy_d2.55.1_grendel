PROGRAM wildcard

  USE messy_main_tools, ONLY: match_wild

  CHARACTER(LEN=80) :: pattern = ''
  CHARACTER(LEN=80) :: string  = ''

  WRITE(*,*) 'PRESS CTRL-C TO STOP.'
  WRITE(*,*)

  WRITE(*,*) 'STRING ?:'
  READ(*,*) string

  DO
     WRITE(*,*) 'PATTERN ?:'
     READ(*,*) pattern

     WRITE(*,*) match_wild(pattern,string)
  END DO

END PROGRAM wildcard
