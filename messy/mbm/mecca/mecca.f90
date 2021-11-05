! Time-stamp: <2018-01-20 14:42:04 sander>

! Author:
! Rolf Sander, MPICH, Mainz, 2014-...

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

MODULE mecca_module

  USE messy_cmn_photol_mem, ONLY: IP_MAX, ip_O1D, ip_NO2
  USE messy_mecca_kpp, ONLY: DP, initialize, kpp_integrate, NSPEC, initialize_kpp_ctrl, &
    fill_temp, fill_cair, fill_press, fill_jx, &
    ind_H2O, ind_OH, ind_O3, ind_NO, ind_NO2, ind_O2, ind_N2, ind_CH4, ind_CO, ind_CO2
  IMPLICIT NONE

  ! declarations for all variables that are transported to the SMCL via fill
  ! subroutines (or via kpp_integrate for "C")
  REAL(DP) :: jx(IP_MAX) = 0._DP  ! J-values
  REAL(DP) :: c(NSPEC)            ! concentrations [mcl/cc]
  REAL(DP) :: cair                ! concentration of air [mcl/cc]
  REAL(DP) :: temp   = 293._DP    ! temperature [K]
  REAL(DP) :: press  = 101325._DP ! pressure [Pa]

  PRIVATE
  PUBLIC :: DP
  PUBLIC :: mecca_init     ! initialize
  PUBLIC :: update_jvalues ! calculate photolysis rate coefficients
  PUBLIC :: mecca_physc    ! calculate chemistry
  PUBLIC :: mecca_result   ! print results

CONTAINS

  !***************************************************************************

  SUBROUTINE mecca_init

    !USE messy_mecca, ONLY: mecca_read_nml_ctrl
    IMPLICIT NONE
    INTEGER :: status ! error status

    CALL initialize_kpp_ctrl(status, 999, 'mecca') ! read kpp ctrl namelist
    IF (status /= 0) STOP
    CALL initialize ! define tolerances rtol and atol

    ! calculate concentration of "air molecules" in [molecules/cm3]:
    cair = 6.02214129E23 * press / ( 1E6 * 8.3144621 *temp )
    c(:) = 0. ! default value unless explicitly initialized
    ! The numbers given here are mixing ratios [mol/mol]. Multiplication
    ! with cair converts to [molecules/cm3].
    c(ind_H2O)     =   1.E-04 * cair
    c(ind_O3)      =  10.E-09 * cair
    c(ind_NO)      =   1.E-09 * cair
    c(ind_O2)      = 210.E-03 * cair
    c(ind_N2)      = 780.E-03 * cair
    c(ind_CH4)     =  1.8E-06 * cair
    c(ind_CO)      =  70.E-09 * cair
    c(ind_CO2)     = 350.E-06 * cair

    WRITE (*,'(5A12)') 'time [s]', 'O3', 'NO', 'NO2', 'OH'

  END SUBROUTINE mecca_init

  !***************************************************************************

  SUBROUTINE update_jvalues

    jx(:) = 0.
    jx(ip_O1D) = 1E-5_DP
    jx(ip_NO2) = 1E-4_DP

  END SUBROUTINE update_jvalues

  !***************************************************************************

  SUBROUTINE mecca_physc(timesteplen)

    IMPLICIT NONE

    REAL, INTENT(IN)  :: timesteplen
    INTEGER, PARAMETER :: NBL = 1 ! N_block_length
    REAL(DP), DIMENSION(NBL,NSPEC) :: cbl
    INTEGER :: status

    ! transfer of data in mecca to kpp:
    CALL fill_temp (status, SPREAD(temp, 1,NBL))
    CALL fill_cair (status, SPREAD(cair, 1,NBL))
    CALL fill_press(status, SPREAD(press,1,NBL))
    CALL fill_jx   (status, SPREAD(jx,   1,NBL))

    c(:) = MAX(c(:),0._DP) ! set negative values to zero
    cbl = SPREAD(c,1,NBL) ! add one dummy dimension
    CALL kpp_integrate(REAL(timesteplen,DP),cbl)  ! main kpp call
    c = cbl(1,:)          ! remove the dummy dimension

  END SUBROUTINE mecca_physc

  !***************************************************************************

  SUBROUTINE mecca_result(time)
    IMPLICIT NONE
    REAL, INTENT(IN)  :: time

    WRITE (*,'(F12.0,4(1PE12.4))') &
      time, c(ind_O3)/cair, c(ind_NO)/cair, c(ind_NO2)/cair, c(ind_OH)/cair

  END SUBROUTINE mecca_result

  !***************************************************************************

END MODULE mecca_module

!*****************************************************************************

PROGRAM mecca

  ! This is the base model layer (BML), i.e. the main box model program

  USE mecca_module
  IMPLICIT NONE
  INTEGER :: timestep
  REAL :: time
  REAL, PARAMETER :: timesteplen = 1200. ! 20 min

  ! initialization:
  time = 0.
  CALL mecca_init
  CALL mecca_result(time) ! print some concentrations at current model time

  ! time loop:
  DO timestep = 1,10
    CALL update_jvalues
    CALL mecca_physc(timesteplen)
    time = time + timesteplen
    CALL mecca_result(time) ! print some concentrations at current model time
  END DO

END PROGRAM mecca

!*****************************************************************************
