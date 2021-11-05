MODULE messy_main_channel_error_bi

  IMPLICIT NONE

  PUBLIC :: channel_halt

CONTAINS

  ! -------------------------------------------------------------------
  SUBROUTINE channel_halt(substr, status)

#ifndef MBM_CHANNEL
    USE messy_main_blather_bi,     ONLY: error_bi
#endif
    USE messy_main_constants_mem,  ONLY: STRLEN_VLONG, iouerr
    USE messy_main_channel_error,  ONLY: channel_error_str

    IMPLICIT NONE

    ! I/O
    CHARACTER(LEN=*), INTENT(IN)  :: substr
    INTEGER,          INTENT(IN)  :: status

    ! LOCAL
    CHARACTER(LEN=STRLEN_VLONG)   :: errstr

    IF (status == 0) RETURN

    errstr = channel_error_str(status)

#ifndef MBM_CHANNEL
    CALL error_bi(errstr,substr)
#else
    WRITE(iouerr,*) errstr
    STOP
#endif
  END SUBROUTINE channel_halt
  ! -------------------------------------------------------------------

END MODULE messy_main_channel_error_bi
