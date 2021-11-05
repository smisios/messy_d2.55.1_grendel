!*******************************************************************************
!                Time-stamp: <2006-08-11 15:37:58 akerkweg>
!*******************************************************************************

! Author:
! Astrid Kerkweg,    MPICH, 2006

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

!*****************************************************************************

PROGRAM sedi

  USE messy_main_constants_mem, ONLY: dp
  USE messy_sedi_column,        ONLY: sedi_column_initialize &
                                    , sedi_column_init_memory  & 
                                    , sedi_column_init_column  &
                                    , sedi_column_physc        &
                                    , sedi_column_result       &
                                    , sedi_column_free_memory  &
                                    , NT

  IMPLICIT NONE

  INTEGER :: t                  ! timeloop

  CALL sedi_column_initialize   ! read namelists
  CALL sedi_column_init_column  ! init meteorology of column

  CALL sedi_column_init_memory  ! initialize model_time, concentrations etc.
                                ! define aerosol profile
  DO t=1,NT                     ! IMITATE TIMELOOP
     CALL sedi_column_result(t)
     CALL sedi_column_physc
  ENDDO
  
  CALL sedi_column_result(t)    ! print results
  CALL sedi_column_free_memory  ! close output files, DEALLOCATE FIELDS
  
END PROGRAM sedi
