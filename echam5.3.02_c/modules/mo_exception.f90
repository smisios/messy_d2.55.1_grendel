MODULE mo_exception

  USE mo_doctor, ONLY: nerr

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: message_text
  PUBLIC :: message, finish
  PUBLIC :: em_none, em_info, em_warn

  INTEGER, PARAMETER :: em_none = 0 
  INTEGER, PARAMETER :: em_info = 1
  INTEGER, PARAMETER :: em_warn = 2

  CHARACTER(512) :: message_text = ''

CONTAINS

  SUBROUTINE finish (name, text, exit_no)

#ifndef MESSY
    USE mo_mpi,           ONLY: p_abort, p_parallel 
#else
    USE mo_mpi,           ONLY: p_abort, p_parallel &
                              , p_pe     ! mz_pj_20041025
#endif

    CHARACTER(*) :: name
    CHARACTER(*), OPTIONAL :: text
    INTEGER, OPTIONAL :: exit_no
    INTEGER           :: iexit
#ifdef MESSY
    INTEGER           :: iou       ! mz_pj_20030707
    LOGICAL           :: opened    ! mz_pj_20030707
    LOGICAL           :: lex            ! mz_pj_20041025
    CHARACTER(LEN=7)  :: efile     ! END0001
#endif

    EXTERNAL util_exit

#ifndef MESSY
    WRITE (nerr,'(/,80("*"),/)')
#endif
    IF (PRESENT(exit_no)) THEN
       iexit = exit_no
    ELSE
       iexit = 1
    END IF

#ifndef MESSY
    IF (PRESENT(text)) THEN
      WRITE (nerr,'(1x,a,a,a)') TRIM(name), ': ', TRIM(text)
    ELSE
      WRITE (nerr,'(1x,a,a)') TRIM(name), ': '
    ENDIF

    WRITE (nerr,'(/,80("*"),/)')
#else

    IF (iexit /=0) WRITE (nerr,'(/,80("*"),/)')  ! mz_pj_20030707

    IF (PRESENT(text)) THEN
       WRITE (nerr,'(1x,a,i4,a,a,a,a)') '(p_pe=',p_pe,') ', &
            TRIM(name), ': ', TRIM(text)
    ELSE
       WRITE (nerr,'(1x,a,i4,a,a,a)') '(p_pe=',p_pe,') ', &
            TRIM(name), ': '
    ENDIF

    IF (iexit /=0) WRITE (nerr,'(/,80("*"),/)')  ! mz_pj_20030707

    ! mz_pj_20030707+
    ! WRITE FILE 'END' FOR BREAKING RERUN-CHAIN IN RUNSCRIPT xmessy
    ! NOTE:
    !   USE mo_filename,  ONLY: find_next_free_unit
    !   cannot be used, due to a circular dependency !
    DO iou=100,300
       INQUIRE(unit=iou,opened=opened)
       IF (.NOT.opened) EXIT
    END DO
    WRITE(efile,'(a3,i4.4)') 'END',p_pe
    !
    ! CHECK IF FILE EXISTS
    INQUIRE(FILE=efile, EXIST=lex)
    IF (lex) THEN  ! FILE EXISTS ...
       ! -> OPEN EXISTING FILE
       OPEN(iou, FILE=efile, STATUS='OLD', POSITION='APPEND')
    ELSE
       ! -> OPEN NEW FILE
       OPEN(iou, FILE=efile, STATUS='UNKNOWN', POSITION='APPEND')
    END IF
    ! OUTPUT MESSAGE
    IF (PRESENT(text)) THEN
       WRITE (iou,'(1x,a,i4,a,a,a,a)') '(p_pe=',p_pe,') ', &
            TRIM(name), ': ', TRIM(text)
    ELSE
       WRITE (iou,'(1x,a,i4,a,a,a)') '(p_pe=',p_pe,') ', &
            TRIM(name), ': '
    ENDIF
    ! CLOSE FILE
    CLOSE(iou)
    ! mz_pj_20030707-
#endif

    IF (p_parallel) THEN 
      CALL p_abort
    ELSE
       CALL util_exit(iexit)
    END IF

  END SUBROUTINE finish

  SUBROUTINE message (name, text, kout, klevel, p_io_c)

#ifdef DEBUG
    USE mo_mpi, ONLY: p_parallel, p_pe
#else
    USE mo_mpi, ONLY: p_parallel_io, p_io, p_pe
#endif

    CHARACTER (*) :: name, text
    INTEGER, INTENT(in), OPTIONAL :: kout
    INTEGER, INTENT(in), OPTIONAL :: klevel
    ! p_io_c is alternative PE for IO to p_io 
    INTEGER, INTENT(in), OPTIONAL :: p_io_c

    INTEGER :: iout
    INTEGER :: ilevel
    INTEGER :: p_io_up

    IF (PRESENT(kout)) THEN
      iout = kout
    ELSE
      iout = nerr
    END IF

    IF (PRESENT(klevel)) THEN
      ilevel = klevel
    ELSE
      ilevel = em_none
    END IF

    SELECT CASE (ilevel)
    CASE (em_none)
    CASE (em_info)
    CASE (em_warn)
    END SELECT

#ifdef DEBUG
    IF (p_parallel) THEN
       IF (name == '') THEN
          WRITE(iout,'(1x,a,i4," - ",a)') 'PE ',p_pe,TRIM(text)
       ELSE
          WRITE(iout,'(1x,a,i4," - ",a,": ",a)') 'PE ', p_pe, TRIM(name), TRIM(text)
       END IF
    ELSE
       IF (name == '') THEN
          WRITE(iout,'(1x,a)') TRIM(text)
       ELSE
          WRITE(iout,'(1x,a,": ",a)') TRIM(name), TRIM(text)
       END IF
    END IF
#else
    p_io_up = p_io
    IF (PRESENT(p_io_c)) p_io_up = p_io_c
    IF (p_pe==p_io_up) THEN
       IF (name == '') THEN
#if defined (_SX)
!mz_bs_20051129+
          WRITE(iout,'(1x,a132)') TRIM(text)
!mz_bs_20051129-
#else
          WRITE(iout,'(1x,a)') TRIM(text)
#endif
       ELSE
          WRITE(iout,'(1x,a,": ",a)') TRIM(name), TRIM(text)
       END IF
    END IF
#endif

  END SUBROUTINE message

END MODULE mo_exception
