! Time-stamp: <2019-07-25 13:14:30 sander>

!> \image latex "../../caaba_mecca_logo_print.pdf" "" width=\textwidth
!> \mainpage Introduction
!> CAABA is a box model that uses MECCA chemistry, plus simplified
!> calculations for emission, deposition, and photolysis.
!> For more information, see \cite 2405
! cite doesn't work if label starts with number

!> \brief CAABA = Chemistry As A Boxmodel Application

!> \authors Rolf Sander, MPICH, Mainz, 2003-2015
!> \authors Hella Riede, MPICH, Mainz, 2007
!> \authors Sergey Gromov, MPICH, Mainz, 2009-2019

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

!*****************************************************************************

PROGRAM caaba

  ! This is the base model layer (BML), i.e. the main box model program

  USE caaba_module,          ONLY: caaba_init, caaba_result, caaba_physc, &
                                   caaba_finish
  USE caaba_mem,             ONLY: model_time, model_end, timesteplen, &
                                   timestep
  USE messy_main_control_cb, ONLY: messy_init, messy_result, messy_physc, &
                                   messy_finish

  IMPLICIT NONE
  CHARACTER(LEN=*), PARAMETER :: statusfile = 'status.log'

  ! write exit status 1 to statusfile:
  OPEN(10, FILE=statusfile, status='UNKNOWN')
  WRITE(10,'(A)') "1"
  CLOSE(10)

  ! initialization:
  CALL caaba_init
  CALL messy_init
  CALL caaba_result
  CALL messy_result

  ! time loop:
  DO WHILE (model_time < model_end)
     CALL caaba_physc
     CALL messy_physc
     timestep = timestep + 1
     model_time = model_time + timesteplen
     CALL caaba_result
     CALL messy_result
  ENDDO

  ! reset timestep so that all output files in *_finish will be written,
  ! independent of output_step_freq and output_sync_freq:
  timestep = 0
  ! final clean up:
  CALL messy_finish
  CALL caaba_finish

  ! program has finished successfully, write exit status 0 to statusfile:
  OPEN(10, FILE=statusfile, status='UNKNOWN')
  WRITE(10,'(A)') "0"
  CLOSE(10)

END PROGRAM caaba

!*****************************************************************************
