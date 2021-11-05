!*****************************************************************************
!                Time-stamp: <2018-05-15 21:28:03 sander>
!*****************************************************************************

! dissoc.f90: a simple box model that calls the dissoc subroutines

! Author:
! Rolf Sander, MPICH, 2018-...

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

MODULE dissoc_column

  USE messy_main_constants_mem, ONLY: DP
  USE messy_cmn_photol_mem,     ONLY: IP_MAX, jname
  USE messy_main_tools,         ONLY: PTR_2D_ARRAY
  USE messy_dissoc

  IMPLICIT NONE
  PRIVATE

  ! PUBLIC dissoc INTERFACE ROUTINES
  PUBLIC :: dissoc_init
  PUBLIC :: dissoc_physc
  PUBLIC :: dissoc_result
  PUBLIC :: dissoc_free_memory

  INTEGER, PARAMETER :: nsza = 10 ! number of zenith angles
  INTEGER, PARAMETER :: nlev = 19 ! number of levels
  REAL, DIMENSION(nsza) :: sza = &
    (/ 0., 20., 40., 60., 70., 80., 85., 90., 92., 94. /)
  REAL(DP), DIMENSION(nlev) :: press_out = &
    (/10.00,     12.92,     16.68,     21.54,     27.83,  &
      35.94,     46.42,     59.95,     77.43,     100.00, &
      129.15,    166.81,    215.44,    278.26,    359.38, &
      464.16,    599.48,    774.26,    900.00/)  ! pressure [hPa]
  REAL(DP), DIMENSION(nlev) :: temp_out = &
    (/ 228.36,    224.99,    221.91,    219.11,    215.84, &
       211.72,    206.74,    200.57,    197.04,    196.14, &
       203.06,    211.95,    222.80,    235.88,    249.06, &
       262.34,    274.18,    285.72,    292.35/) ! temperature [K]

  REAL(DP), DIMENSION(IP_MAX)  :: photoarr

CONTAINS

  ! --------------------------------------------------------------------------

  SUBROUTINE dissoc_init

    IMPLICIT NONE

    ! LOCAL
    INTEGER :: status ! status flag
    INTEGER :: ip
    INTEGER :: month=7

    ! INTITIALIZE GLOBAL SWITCHES / PARAMETERS
    CALL dissoc_read_nml_ctrl(status, 99)
    IF (status /= 0) STOP

    ALLOCATE(dissoc_2d(IP_MAX))
    DO ip=1, ip_MAX
      ALLOCATE(dissoc_2d(ip)%ptr(nsza,nlev))
    END DO

    ALLOCATE(dtemp(jpslevall,jplats))
    ALLOCATE(dtempc(jpslev,jplats))
    ALLOCATE(alt(jpslevall,jplats))
    ALLOCATE(altc(jpslev,jplats))
    ALLOCATE(do2(jpslevall,jplats))
    ALLOCATE(do2c(jpslev,jplats))
    ALLOCATE(do3(jpslevall,jplats))
    ALLOCATE(do3c(jpslev,jplats))
    ALLOCATE(angdeg(jpschi))
    ALLOCATE(wavenm(jpwave))
    ALLOCATE(lats(jplats))
    ALLOCATE(pres(jpslevall))
    ALLOCATE(presc(jpslev))
    ALLOCATE(tabs_davg(jpslev, jpwave, jplats))
    ALLOCATE(tabs(jpslev, jpschi, jpwave, jplats))
    CALL iniphoto(0)
    CALL reado3(o3dat,month,nlats)
    CALL setp
    CALL settab

  END SUBROUTINE dissoc_init

  ! --------------------------------------------------------------------------

  SUBROUTINE dissoc_physc

    IMPLICIT NONE

    INTEGER :: i,j, ip
    REAL(DP) :: lat_calc = 17.5000

    !ju_jg_20181212
    ! Don't set the albedo here. It should be read in through the namelist dissoc.nml
    !!albedo  = 0.07

    ! calculate jvalues:
    DO i=1,nsza
      DO j=1,nlev
        CALL dissoc(temp_out(j), press_out(j), lat_calc, sza(i)*dtr, photoarr)
        DO ip=1, IP_MAX
          IF (jcalc(ip)) dissoc_2d(ip)%ptr(i,j) = photoarr(ip)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE dissoc_physc

  ! --------------------------------------------------------------------------

  SUBROUTINE dissoc_result

    USE mo_netcdf, ONLY: open_dissoc_nc_file,  &
      write_dissoc_nc_file, &
      close_dissoc_nc_file

    INTEGER :: ncid_dissoc ! netcdf id for dissoc.nc
    INTEGER :: ip

    CALL open_dissoc_nc_file(ncid_dissoc, nlev, REAL(press_out), nsza, sza)

    DO ip=1, ip_MAX
      IF (jcalc(ip)) CALL write_dissoc_nc_file(ncid_dissoc, 'J_'//TRIM(jname(ip)), &
        dissoc_2d(ip)%ptr(:,:))
    ENDDO

    CALL close_dissoc_nc_file(ncid_dissoc)

  END SUBROUTINE dissoc_result

  ! --------------------------------------------------------------------------

  SUBROUTINE dissoc_free_memory

    IMPLICIT NONE

    INTEGER :: ip

    DO ip=1, ip_MAX
      DEALLOCATE(dissoc_2d(ip)%ptr)
      NULLIFY(dissoc_2d(ip)%ptr)
    ENDDO
    DEALLOCATE(dissoc_2d)

    DEALLOCATE(dtemp)
    DEALLOCATE(dtempc)
    DEALLOCATE(alt)
    DEALLOCATE(altc)
    DEALLOCATE(do2)
    DEALLOCATE(do2c)
    DEALLOCATE(do3)
    DEALLOCATE(do3c)
    DEALLOCATE(angdeg)
    DEALLOCATE(wavenm)
    DEALLOCATE(lats)
    DEALLOCATE(pres)
    DEALLOCATE(presc)
    DEALLOCATE(tabs_davg)
    DEALLOCATE(tabs)
    ! DEALLOCATE variables from messy_dissoc.f90:
    DEALLOCATE(pres_inp)
    DEALLOCATE(alt_inp)
    DEALLOCATE(temp_inp)
    DEALLOCATE(o3_inp)

  END SUBROUTINE dissoc_free_memory

  ! --------------------------------------------------------------------------

END MODULE dissoc_column

!*****************************************************************************

PROGRAM dissoc

  USE dissoc_column, ONLY: dissoc_init, &
    dissoc_physc, dissoc_result, dissoc_free_memory

  IMPLICIT NONE

  PRINT *, 'Starting dissoc_init...'
  CALL dissoc_init ! read CTRL namelist, ALLOCATE, iniphoto, reado3, setp, settab
  PRINT *, 'Starting dissoc_physc...'
  CALL dissoc_physc        ! calculate J values
  PRINT *, 'Starting dissoc_result...'
  CALL dissoc_result       ! print results
  PRINT *, 'Starting dissoc_free_memory...'
  CALL dissoc_free_memory  ! DEALLOCATE

END PROGRAM dissoc

!*****************************************************************************
