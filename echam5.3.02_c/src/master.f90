PROGRAM master

  !----------------------------------------------------------------------------
  !
  ! Copyright 2000-2005 by Max Planck Institute for Meteorology
#ifdef MESSY
  ! Modular Earth Submodel System (MESSy):
  ! Copyright 2000-2005 by MESSy Consortium, see http://www.messy-interface.org
#endif
  !
  ! This software can be used only under the conditions of the 
  ! "MPI-M Software Licence Agreement", which must be signed 
  ! by each user.
  !
  ! The "MPI-M Software Licence Agreement" and the information on
  ! the distribution of MPI-M models are given on:
  !
  ! http://www.mpimet.mpg.de/en/extra/models/distribution/index.php
  !
  ! The ECHAM5-webpage can be found at:
  !
  ! http://www.mpimet.mpg.de/en/extra/models/echam/echam5.php
  !
  !----------------------------------------------------------------------------
  !
  ! Call the control subroutine (*control*).
  !
  ! Externals:
  !
  ! *control*   called to control the run.
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, February 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! I. Kirchner, MPI, October 2000, date/time control
  ! S. Legutke, MPI M&D, Juli 2001, redirect stdout for coupling
  ! 
  ! for more details see file AUTHORS

  USE mo_kind,         ONLY: dp
  USE mo_doctor,       ONLY: nout, nin
#ifdef MESSY
  USE messy_main_constants_mem, ONLY: modver
  ! op_pj_20110216+
  USE messy_main_compilerinfo_mem, ONLY: compiler_version, compiler_call &
                                       , compiler_flags &
                                       , compiler_cppdefs, compiler_includes
  ! op_pj_20110216-
#endif
  USE mo_mpi,          ONLY: p_start, p_stop, p_pe, p_io
  USE mo_time_control, ONLY: lbreak, lstop, labort
  USE mo_exception,    ONLY: finish, message
  USE mo_util_string,  ONLY: separator
#if defined (__oasis)
  USE mo_couple,       ONLY: couple_quit
#endif

  IMPLICIT NONE

  !  External functions 
  REAL(dp), EXTERNAL :: util_walltime

  !  External subroutines 
  EXTERNAL control
#ifdef __XT3__
  EXTERNAL :: util_base_iobuf
#endif

  REAL(dp) :: zwtime

#ifdef __XT3__
  ! set buffer for stderr and stdout

  CALL util_base_iobuf
#endif

  ! Initialize wallclock timer

!$OMP PARALLEL
!$OMP MASTER
  zwtime = util_walltime()
!$OMP END MASTER
!$OMP END PARALLEL

  ! Start MPI

#ifndef MESSY
  CALL p_start
#else
  CALL p_start('EMAC')
#endif

  IF (p_pe == p_io) THEN
#if defined (__prism) 
!--  Redirect standard output file to atmout if coupled run.
 
        OPEN (UNIT=nout,FILE='atmout',STATUS='UNKNOWN',FORM ='FORMATTED')
        WRITE(nout,*)' '
        WRITE(nout,*)'Atmosphere standard output is assigned to file atmout.'
        WRITE(nout,*)' '
#endif

  ! Print version

        
     WRITE (nout,separator)

#ifndef MESSY
     WRITE (nout,'(/,a,/,a,/)')                                        &
          '  ECHAM         - Release  5.3.02 ',                        &
          '  Copyright by Max-Planck-Institute for Meteorology, 2005', &
          '  Read master.f90 and licence.pdf before using ECHAM5'

#else
     ! mz_rs_20030704+
     WRITE (nout,*) '  ECHAM         - Release 5.3.02'
     WRITE (nout,*) '  Copyright by Max-Planck-Institute for Meteorology, 2005'
     WRITE (nout,*) '  Read master.f90 and licence.pdf before using ECHAM5'
     WRITE (nout,*) '  MESSy         - Release '//modver
     WRITE (nout,*) '  Copyright by MESSy Consortium, 2005'
     WRITE (nout,*) '  Read www.messy-interface.org before using MESSy'
     ! op_pj_20110216+
     WRITE (nout,*) ''
     WRITE (nout,*) '  f95 compiler version : '//compiler_version
     WRITE (nout,*) '  f95 compiler call    : '//compiler_call
     WRITE (nout,*) '  f95 compiler flags   : '//compiler_flags
     WRITE (nout,*) '  f95 compiler pp.defs : '//compiler_cppdefs
     WRITE (nout,*) '  f95 compiler includes: '//compiler_includes
     WRITE (nout,*) ''
     ! op_pj_20110216-
     ! mz_rs_20030704-
#endif
     WRITE (nout,*)
     WRITE (nout,separator)
     WRITE (nout,*)
  END IF

  DO                               ! Loop over rerun cycles
     CALL control
     IF (lbreak .OR. lstop) EXIT
     CALL message('','Start next rerun cycle.')
  END DO

 CALL p_stop                      ! Stop MPI

  IF (lstop) THEN
#ifndef MESSY
    CALL message('','Experiment finished.')
    IF (labort) CALL finish('master','Run terminated.',1)
  ELSE
    CALL message('','Experiment stopped.')
#else
    ! mz_pj_20030707+
    CALL message('','Simulation finished.')          ! mz_pj_20030707
    CALL finish('master','Simulation finished.',0)   ! mz_pj_20030707
    ! Notes: 
    !    - simulation is finished (lstop), if final simulation date/step
    !      is reached
    !    - finish is needed in addition to message, in order to trigger
    !      the output of file 'END'
    !      (this breaks always the rerun chain, labort has no meaning)
    ! mz_pj_20030707-
  ELSE
    CALL message('','Simulation stopped.')
    ! mz_pj_20030707+
    ! Notes:
    !    - simulation is stopped (lbreak) and a rerun is started
    !      (continue rerun chain)
    !    - if labort=.true. in the namelist, the chain is broken in any case
    !      (test modus)
    IF (labort) CALL finish('master','Simulation terminated.')
    ! mz_pj_20030707-
#endif
  END IF

END PROGRAM master
