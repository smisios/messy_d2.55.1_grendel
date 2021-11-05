SUBROUTINE iorestart

  ! Description:
  !
  ! Reads netCDF history files for a resumed run.
  !
  ! Method:
  !
  ! *iorestart* positions data sets at the beginning of a rerun,
  ! writing data description records, and setting up necessary work
  ! files.
  !
  ! Information is written to the data description records of
  ! appropriate files, and work files are written if necessary.
  !
  !
  ! Authors:
  !
  ! L. Kornblueh, MPI, May 1999, f90 rewrite
  ! U. Schulzweida, MPI, May 1999, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

#ifndef MESSY
  USE mo_tracer,        ONLY: xtini ! mz_pj_20040329
#endif
  USE mo_doctor,        ONLY: nout
  USE mo_io,            ONLY: IO_read_streams
  USE mo_memory_gl,     ONLY: xt
  USE mo_memory_g1a,    ONLY: xtm1
  USE mo_control,       ONLY: ltimer
  USE mo_timer,         ONLY: timer_start, timer_stop, timer_netcdf
  USE mo_mpi,           ONLY: p_pe, p_io

  IMPLICIT NONE


  !  Executable statements 

  ! Restart from history files

  ! Read all restart files

  IF (ltimer) CALL timer_start(timer_netcdf)
  CALL IO_read_streams
  IF (ltimer) CALL timer_stop(timer_netcdf)
#ifdef MESSY
  ! mz_pj_20050419+
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  CALL messy_read_restart
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! mz_pj_20050419-
#endif

#ifndef MESSY
  CALL xtini (xt, xtm1)
#else
! mz_pj_20040422+
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  CALL messy_init_tracer
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! mz_pj_20040422-
#endif

  IF (p_pe == p_io) &
       WRITE (nout,'(/)')

END SUBROUTINE iorestart
