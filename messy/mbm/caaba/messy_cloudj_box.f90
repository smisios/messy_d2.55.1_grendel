!*****************************************************************************
!                Time-stamp: <2018-05-16 14:19:10 sander>
!*****************************************************************************

! cloudj = calculation of J-values (photolysis rate coefficients)

! Authors:
! Rolf Sander, MPICH, 2016

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

MODULE messy_cloudj_box

  USE messy_main_constants_mem, ONLY: DP
  USE caaba_io,  ONLY: open_output_file, write_output_file, &
                       close_file
  USE messy_cmn_photol_mem      ! IP_MAX, ip_*, jname
  USE messy_cloudj

  IMPLICIT NONE
  PRIVATE

  INTEGER :: ncid_cloudj

  INTEGER, PARAMETER :: nsza =  1
  INTEGER :: nlev ! number of levels

  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: v3
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: jpress  ! given pressure levels
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: relo3
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: jrhum
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: jtemp
  REAL(DP), ALLOCATABLE, DIMENSION(:)   :: albedo
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: aclc
  REAL(DP), ALLOCATABLE, DIMENSION(:)   :: slf
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: clp

  PUBLIC :: cloudj_init   ! initialize J-values
  PUBLIC :: cloudj_physc  ! calculate J values
  PUBLIC :: cloudj_result
  PUBLIC :: cloudj_finish

CONTAINS

  ! --------------------------------------------------------------------------

  SUBROUTINE cloudj_init

    USE caaba_mem ! DEFAULT_*, photo_scenario, l_input_jval
    IMPLICIT NONE
    INTEGER :: status ! status flag
    INTEGER :: ip
    ! longjname must be 3 characters longer than jname:
    CHARACTER(LEN=15), DIMENSION(IP_MAX) :: longjname

    ! read cloudj ctrl namelist:
    CALL cloudj_read_nml_ctrl(status, 999)
    IF (status /= 0) STOP 1

    lp(:)=.TRUE.

    DO ip=1, ip_MAX
      longjname(ip) = 'J('//TRIM(jname(ip))//')'
    ENDDO
    CALL open_output_file(ncid_cloudj, 'caaba_cloudj', &
      (/ ('J_'//jname(ip), ip=1,IP_MAX) /), &
      (/ ('1/s',           ip=1,IP_MAX) /), &
      (/ (longjname(ip),   ip=1,IP_MAX) /))

    nlev = DEFAULT_NLEV
    ALLOCATE(cloudj_2d(ip_MAX))
    DO ip=1, ip_MAX
      ALLOCATE(cloudj_2d(ip)%ptr(nsza,nlev))
    ENDDO

    ALLOCATE(v3(nsza,nlev+1))
    ALLOCATE(jpress(nsza,nlev))
    ALLOCATE(relo3(nsza,nlev+1))
    ALLOCATE(jrhum(nsza,nlev))
    ALLOCATE(jtemp(nsza,nlev))
    ALLOCATE(albedo(nsza))
    ALLOCATE(aclc(nsza,nlev))
    ALLOCATE(slf(nsza))
    ALLOCATE(clp(nsza,nlev))

    v3(1,:)     = DEFAULT_V3     ! vertical ozone column [mcl/cm2]
    relo3(1,:)  = DEFAULT_RELO3  ! rel. ozone, i.e. O3 mixing ratio [mol/mol]
    jpress(1,:) = DEFAULT_JPRESS ! pressure [Pa]
    jrhum(1,:)  = DEFAULT_JRHUM  ! relative humidity [%]
    jtemp(1,:)  = DEFAULT_JTEMP  ! temperature [K]
    albedo(1)   = DEFAULT_ALBEDO ! albedo
    aclc(1,:)   = 0.             ! assume clear sky
    slf(1)      = 0.             ! 0 = sea
    clp(1,:)    = 0.             ! cloud liquid water path [g/m^2] (clear sky)

    CALL init_cloudj

  END SUBROUTINE cloudj_init

  ! --------------------------------------------------------------------------

  SUBROUTINE cloudj_physc

    USE caaba_mem,        ONLY: cossza
    USE messy_main_tools, ONLY: nn_index

    IMPLICIT NONE

    ! calculate J-values
    CALL cloudjvalues(                               &
      REAL(v3),                                      &
      REAL((/ cossza /)), REAL(jpress), REAL(relo3), &
      REAL(jrhum), REAL(jtemp), REAL(albedo),        &
      REAL(aclc), REAL(slf), REAL(clp))

  END SUBROUTINE cloudj_physc

  ! --------------------------------------------------------------------------

  SUBROUTINE cloudj_finish

    IMPLICIT NONE

    INTEGER :: ip

    DO ip = 1, IP_MAX
      DEALLOCATE(cloudj_2d(ip)%ptr)
      NULLIFY(cloudj_2d(ip)%ptr)
    ENDDO
    DEALLOCATE(cloudj_2d)    ; NULLIFY(cloudj_2d)

    DEALLOCATE(v3)
    DEALLOCATE(jpress)
    DEALLOCATE(relo3)
    DEALLOCATE(jrhum)
    DEALLOCATE(jtemp)
    DEALLOCATE(albedo)
    DEALLOCATE(aclc)
    DEALLOCATE(slf)
    DEALLOCATE(clp)

    CALL close_file(ncid_cloudj)

  END SUBROUTINE cloudj_finish

  ! --------------------------------------------------------------------------

  SUBROUTINE cloudj_result

    USE caaba_mem, ONLY: model_time, percent_done, photol_clev
    USE messy_main_constants_mem, ONLY: FLAGGED_BAD

    IMPLICIT NONE

    REAL(DP), DIMENSION(IP_MAX) :: cloudj_array
    INTEGER :: ip

    IF (percent_done > 0.) THEN
      ! The following loop is necessary because the component (after the
      ! "%") to the right of a part reference with nonzero rank cannot
      ! have the POINTER attribute. In other words, we cannot use
      ! "cloudj_2d(:)%ptr(1,photol_clev,1)"
      DO ip = 1, IP_MAX
        cloudj_array(ip) = cloudj_2d(ip)%ptr(1,photol_clev)
      ENDDO
      CALL write_output_file(ncid_cloudj, model_time, cloudj_array)
    ELSE
      CALL write_output_file(ncid_cloudj, model_time, &
        (/ (FLAGGED_BAD, ip=1,IP_MAX) /))
    ENDIF

  END SUBROUTINE cloudj_result

  ! --------------------------------------------------------------------------

END MODULE messy_cloudj_box

!*****************************************************************************
