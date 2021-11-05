! ****************************************************************************
!                Time-stamp: <2014-02-13 12:05:31 sander>
!     Author: Rolf Sander, Max-Planck Institute, Mainz, Germany, 2009-2014
!*****************************************************************************

! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with this program; if not, get it from:
! http://www.gnu.org/copyleft/gpl.html

! ****************************************************************************

PROGRAM jvpp

  USE jvpp_mem,   ONLY: IO_NML, IO_LOG, DIR_JNL, inputdir
  USE jvpp_step1, ONLY: conv_sig_176
  USE jvpp_step2, ONLY: conv_176_eff
  USE jvpp_step3, ONLY: conv_eff_poly

  IMPLICIT NONE
  CHARACTER(LEN=80), PARAMETER :: NML_FILENAME = 'jvpp.nml'
  NAMELIST /CTRL/ inputdir

  ! read namelist file:
  OPEN(IO_NML, FILE=NML_FILENAME)
  READ(IO_NML, NML=CTRL)
  CLOSE(IO_NML)
  SELECT CASE(inputdir)
  CASE ("dat_m17", "dat_lit") ! ok
    PRINT *, 'inputdir = ' // inputdir
  CASE DEFAULT ! not ok:
    PRINT *, 'error while reading inputdir from namelist ' // NML_FILENAME
    STOP
  END SELECT

  OPEN(IO_LOG, FILE='jvpp_detail.log', STATUS='UNKNOWN')

  CALL conv_sig_176  ! step 1: sigma from literature --> 176 intervals
  CALL conv_176_eff  ! step 2: 176 intervals --> effective values
  CALL conv_eff_poly ! step 3: effective values --> polynomial fits

  CLOSE(IO_LOG)

  WRITE(*,*)
  WRITE(*,'(A)') 'To view results with ferret, type:'
  WRITE(*,'(A)') 'cd '//DIR_JNL//' ; ls -la *.jnl ; ferret'

END PROGRAM jvpp

! ****************************************************************************
