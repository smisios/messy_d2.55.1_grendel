SUBROUTINE savehis

  ! Description:
  !
  ! This routine saves history files for future use.
  !
  ! Method:
  !
  ! The history files are put together with "tar".
  ! The name of the tar file is: exp.number"re"yymm
  !
  ! *savehis* is called from *posts2*
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, September 1991, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! U. Schlese, DKRZ, July 1999, modifications for ECHAM5
  ! I. Kirchner, MPI, March 2001, date and time control
  ! A. Rhodin, DWD/MPI, November 2001, changed for output stream interface
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_exception,     ONLY: finish
  USE mo_filename,      ONLY: yexp, standard_rerun_file
  USE mo_doctor,        ONLY: nout
  USE mo_time_control,  ONLY: current_date, str_date
  USE mo_memory_base,   ONLY: ostreams, nstreams

  IMPLICIT NONE

  INTEGER         :: ilenc, ilenf2, i, number
  CHARACTER (512) :: yocmnd
  CHARACTER (40)  :: yofile2

  INTEGER, EXTERNAL :: util_system

  yofile2 = yexp(1:5) // 're' // TRIM(str_date(4,current_date))
  ilenf2 = MIN(LEN_TRIM(yofile2),40)

#if defined(__SX__)
    WRITE (nout,*) " savehis: tar of history files failed! (f90 bug)"
    CALL finish('savehis','Run terminated.')
#else
    !
    ! write tar command, name of the tar file
    !
    yocmnd = 'tar cvf ' // yofile2(:ilenf2)
    !
    ! add names of the rerun files from the output buffers
    !
    number = 0
    DO i = 1, nstreams
       if (.not. ostreams(i)% lrerun)                               cycle 
       if (any(ostreams(1:i-1)% rest_suf == ostreams(i)% rest_suf)) cycle
       yocmnd =  TRIM(yocmnd) // ' ' // &
                 TRIM(standard_rerun_file) // ostreams(i)% rest_suf
       number = number + 1
    END DO

#endif

  ilenc = MIN(LEN_TRIM(yocmnd),512)
  IF (util_system(yocmnd(:ilenc)) /= 0) THEN
    WRITE (nout,*) ' ERROR: tar of history files failed! Job aborted!'
    CALL finish('savehis','Run terminated.')
  END IF

  WRITE (nout,*) number, ' history files saved as ' // yofile2(:ilenf2)

END SUBROUTINE savehis

