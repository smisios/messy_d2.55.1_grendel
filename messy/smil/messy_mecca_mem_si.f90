!*******************************************************************************
!                Time-stamp: <2009-11-12 12:22:14 akerkweg>
!*******************************************************************************

! Author:
! Rolf Sander, MPICH, 2003:   original code

! This module provides echam-specific variables for the mecca submodel

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

MODULE messy_mecca_mem_si

  ! MESSy
  USE messy_main_constants_mem, ONLY: dp

  IMPLICIT NONE

  ! PUBLIC is already default
  PRIVATE :: dp

  SAVE

  ! the INCLUDE file is produced via xmecca. It contains the idt_* declarations
  INCLUDE 'messy_mecca_idt_si.inc'

  INTEGER :: jk, jp

END MODULE messy_mecca_mem_si
