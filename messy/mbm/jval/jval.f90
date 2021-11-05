!*****************************************************************************
!                Time-stamp: <2015-06-22 19:16:28 sander>
!*****************************************************************************

! Author:
! Rolf Sander, MPICH, 2004-...

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

MODULE jval_column

  USE messy_main_constants_mem, ONLY: DP
  USE messy_cmn_photol_mem, ONLY: IP_MAX, jname
  USE messy_jval ! ONLY: jval_2d, jname, ip_MAX, jval_read_nml_ctrl, 
                 !       aerosol_data, jvalues, lp

  IMPLICIT NONE
  PRIVATE

  ! PUBLIC jval INTERFACE ROUTINES
  PUBLIC :: jval_initialize
  PUBLIC :: jval_init_memory
  PUBLIC :: jval_physc
  PUBLIC :: jval_result
  PUBLIC :: jval_free_memory

  INTEGER, PARAMETER :: nsza = 10
  INTEGER, PARAMETER :: nlev = 19 ! number of levels
  REAL, DIMENSION(nsza) :: sza = &
    (/ 0., 10., 20., 30., 40., 50., 60., 70., 80., 90. /)

CONTAINS

  ! --------------------------------------------------------------------------

  SUBROUTINE jval_initialize

    IMPLICIT NONE

    ! LOCAL
    INTEGER :: status ! status flag

    ! INTITIALIZE GLOBAL SWITCHES / PARAMETERS
    CALL jval_read_nml_ctrl(status, 99)
    IF (status /= 0) STOP

    ! intialize aerosol data
    CALL aerosol_data ! aerosol optical data (Shettle, Fenn)

    lp(:)=.TRUE.

  END SUBROUTINE jval_initialize

  ! --------------------------------------------------------------------------

  SUBROUTINE jval_init_memory

    IMPLICIT NONE

    INTEGER :: jt

    ALLOCATE(jval_2d(ip_MAX))
    DO jt=1, ip_MAX
       ALLOCATE(jval_2d(jt)%ptr(nsza,nlev))
    END DO

    ALLOCATE(rh_o2_2d(nsza,nlev))
    ALLOCATE(rh_o3_2d(nsza,nlev))

    ALLOCATE(fhuv_2d(nsza,nlev))
    ALLOCATE(fhuvdna_2d(nsza,nlev))

  END SUBROUTINE jval_init_memory

  ! --------------------------------------------------------------------------

  SUBROUTINE jval_physc

    USE messy_main_constants_mem, ONLY: pi

    IMPLICIT NONE

    INTEGER :: i

    REAL, DIMENSION(nsza,nlev+1) :: v3
    REAL, DIMENSION(nsza,nlev)   :: press
    REAL, DIMENSION(nsza,nlev+1) :: relo3
    REAL, DIMENSION(nsza,nlev)   :: rhum
    REAL, DIMENSION(nsza,nlev)   :: temp
    REAL, DIMENSION(nsza)        :: albedo
    REAL, DIMENSION(nsza,nlev)   :: aclc
    REAL, DIMENSION(nsza)        :: slf
    REAL, DIMENSION(nsza,nlev)   :: clp
    LOGICAL                      :: lmidatm
    LOGICAL                      :: l_heating
    INTEGER                      :: status
    INTEGER                      :: pbllev   ! number of levels in pbl

    DO i=1,nsza
      ! global average values are extracted with ferret from messy and
      ! jval_diag streams using e.g.: "list rhum[i=@ave,j=@ave,l=1]"

      ! vertical ozone column [mcl/cm2]
      v3(i,:)    = (/ &
        3.366E+17, 1.437E+18, 4.085E+18, 5.428E+18, 6.157E+18, 6.583E+18, &
        6.860E+18, 7.070E+18, 7.227E+18, 7.343E+18, 7.436E+18, 7.523E+18, &
        7.605E+18, 7.678E+18, 7.740E+18, 7.788E+18, 7.822E+18, 7.844E+18, &
        7.857E+18, 7.862E+18 /)
      ! relative ozone, i.e. ozone mixing ratio [mol/mol]
      ! Note that although relo3 has the dimension 1:nlev+1, the value
      ! relo3(1) is not used at all here. Also, note that relo3 is
      ! _only_ used for the heating rates. For the calculation of the
      ! J-values, only v3 is used.
      relo3(i,:) = (/ &
        7.182E-06, 8.319E-06, 4.172E-06, 2.041E-06, 9.525E-07, 4.334E-07, &
        2.571E-07, 1.514E-07, 9.760E-08, 5.775E-08, 5.064E-08, 4.394E-08, &
        3.980E-08, 3.636E-08, 3.209E-08, 2.807E-08, 2.479E-08, 2.242E-08, &
        2.105E-08, 2.065E-08 /)
      ! pressure [Pa]
      press(i,:) = (/ &
        1000., 3000., 5040., 7339., 10248., 14053., 18935., 24966., 32107., &
        40212., 49027., 58204., 67317., 75897., 83472., 89631., 94099.,     &
        96838., 98169. /)
      ! relative humidity [%]
      rhum(i,:)  = (/ &
        0.23, 1.57, 3.52, 11.73, 24.55, 25.31, 27.45, 36.46, 44.52, 46.27,  &
        46.48, 49.18, 51.73, 57.95, 72.82, 80.71, 81.66, 77.65, 76.18 /)
      ! temperature [K]
      temp(i,:)  = (/ &
        230.6, 218.2, 211.7, 207.0, 205.6, 210.9, 218.1, 225.8, 235.7, 246.6, &
        256.4, 264.2, 270.6, 275.4, 278.2, 280.9, 283.2, 284.9, 285.7 /)
    ENDDO
    albedo(:)  = 0.07
    aclc(:,:)  = 0.            ! assume clear sky
    slf(:)     = 0.            ! 0 = sea
    ! clp = cloud liquid water path per layer [g/m^2]
    clp(:,:)   = 0.            ! assume clear sky
    lmidatm    = .FALSE.
    l_heating  = .TRUE.
    pbllev     = 5

    ! use r_sol [0,...1] in CTRL for solar cycle
    ! orbital parameter is set to 1.0 AU here (no orbital variation)
    ! no external solar cycle data provided here
    CALL jval_solar_time_control(status, 1.0_DP)
    IF (status /= 0) STOP

    ! calculate jvalues and uv-heating rates
    CALL jvalues(                           &
      v3,                                   &
      COS(sza*REAL(pi)/180.), press, relo3, &
      rhum, temp, albedo,                   &
      aclc, slf, clp,                       &
      lmidatm, l_heating, pbllev)

  END SUBROUTINE jval_physc

  ! --------------------------------------------------------------------------

  SUBROUTINE jval_result

    USE mo_netcdf, ONLY: open_jval_nc_file,  &
                         write_jval_nc_file, &
                         close_jval_nc_file

    INTEGER :: ncid_jval ! netcdf id for jval.nc
    INTEGER :: jt

    CALL open_jval_nc_file(ncid_jval, nlev, nsza, sza)

    DO jt=1, ip_MAX
       CALL write_jval_nc_file(ncid_jval, 'J_'//TRIM(jname(jt)), &
            jval_2d(jt)%ptr(:,:))
    ENDDO

    CALL write_jval_nc_file(ncid_jval, 'RH_O2', rh_O2_2d(:,:))
    CALL write_jval_nc_file(ncid_jval, 'RH_O3', rh_O3_2d(:,:))

    CALL close_jval_nc_file(ncid_jval)

  END SUBROUTINE jval_result

  ! --------------------------------------------------------------------------

  SUBROUTINE jval_free_memory

    IMPLICIT NONE

    INTEGER :: jt

    DO jt=1, ip_MAX
       DEALLOCATE(jval_2d(jt)%ptr)
       NULLIFY(jval_2d(jt)%ptr)
    ENDDO
    DEALLOCATE(jval_2d)
    NULLIFY(jval_2d)

    DEALLOCATE(rh_o2_2d)      ; NULLIFY(rh_o2_2d)       
    DEALLOCATE(rh_o3_2d)      ; NULLIFY(rh_o3_2d)       

    DEALLOCATE(fhuv_2d)       ; NULLIFY(fhuv_2d)
    DEALLOCATE(fhuvdna_2d)    ; NULLIFY(fhuvdna_2d)

  END SUBROUTINE jval_free_memory

  ! --------------------------------------------------------------------------

END MODULE jval_column

!*****************************************************************************

PROGRAM jval

  USE jval_column, ONLY: jval_init_memory, jval_initialize, &
                         jval_physc, jval_result, jval_free_memory

  IMPLICIT NONE

  CALL jval_initialize   ! read CTRL namelist, intialize aerosol
  CALL jval_init_memory
  CALL jval_physc        ! calculate J values
  CALL jval_result       ! print results
  CALL jval_free_memory

END PROGRAM jval

!*****************************************************************************
