!*****************************************************************************
!                Time-stamp: <2018-05-17 16:37:17 sander>
!*****************************************************************************

! dissoc_box = smil file for DISSOC photolysis submodel

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

MODULE messy_dissoc_box

  USE messy_main_constants_mem, ONLY: DP
  USE caaba_io,  ONLY: open_output_file, write_output_file, close_file
  USE messy_cmn_photol_mem ! IP_MAX, ip_*, jname
  USE messy_dissoc

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dissoc_init
  PUBLIC :: dissoc_physc
  PUBLIC :: dissoc_result
  PUBLIC :: dissoc_finish

  INTEGER :: ncid_dissoc
  INTEGER, PARAMETER :: nsza =  1
  INTEGER :: nlev ! number of levels
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: jtemp
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: jpress
  REAL(DP), DIMENSION(IP_MAX)  :: photoarr

CONTAINS

  ! --------------------------------------------------------------------------

  SUBROUTINE dissoc_init

    USE caaba_mem
    USE messy_main_constants_mem, ONLY: OneDay
    
    IMPLICIT NONE
    INTEGER :: ip, status
    INTEGER :: month
    ! longjname must be 3 characters longer than jname:
    CHARACTER(LEN=15), DIMENSION(IP_MAX) :: longjname

    ! read dissoc ctrl namelist:
    CALL dissoc_read_nml_ctrl(status, 999)
    IF (status /= 0) STOP 1

    jcalc(:)=.TRUE.

    DO ip=1, ip_MAX
      longjname(ip) = 'J('//TRIM(jname(ip))//')'
    ENDDO
    CALL open_output_file(ncid_dissoc, 'caaba_dissoc', &
      (/ ('J_'//jname(ip), ip=1,IP_MAX) /), &
      (/ ('1/s',           ip=1,IP_MAX) /), &
      (/ (longjname(ip),   ip=1,IP_MAX) /))

    nlev = DEFAULT_NLEV
    ALLOCATE(jtemp(nsza,nlev))
    ALLOCATE(jpress(nsza,nlev))
    jtemp(1,:)  = DEFAULT_JTEMP  ! temperature [K]
    jpress(1,:) = DEFAULT_JPRESS ! pressure [Pa]
    ALLOCATE(dissoc_2d(IP_MAX))
    DO ip = 1, IP_MAX
      ALLOCATE(dissoc_2d(ip)%ptr(nsza,nlev))
      dissoc_2d(ip)%ptr(:,:) = 0.
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
    ! approximate current month for O3 climatology:
    month = INT(model_time * 12. / (365.25*OneDay)) + 1
    CALL reado3(o3dat,month,nlats)
    CALL setp
    CALL settab

  END SUBROUTINE dissoc_init

  ! --------------------------------------------------------------------------

  SUBROUTINE dissoc_physc

    USE caaba_mem,        ONLY: cossza, press, photol_clev, DEFAULT_ALBEDO, degree_lat
    USE messy_main_tools, ONLY: nn_index

    IMPLICIT NONE
    INTEGER :: i, j, ip

    albedo = DEFAULT_ALBEDO
    ! calculate jvalues:
    DO i=1,nsza
      DO j=1,nlev
        ! note that SUBROUTINE dissoc needs p in hPa
        CALL dissoc(jtemp(i,j), jpress(i,j)/100., degree_lat, ACOS(cossza), photoarr)
        DO ip=1, IP_MAX
          IF (jcalc(ip)) dissoc_2d(ip)%ptr(i,j) = photoarr(ip)
        ENDDO
      ENDDO
    ENDDO

    ! calculate pressure level in jpress according to current pressure
    CALL nn_index(jpress(1,:), press, photol_clev)
    IF (press > 1.1 * jpress(1,nlev)) THEN
      WRITE (*,*) 'Warning dissoc_physc: pressure more than 10% above '// &
        'highest standard atmosphere pressure'
    ENDIF
    IF (press < 0.9 * jpress(1,1)) THEN
      WRITE (*,*) 'Warning dissoc_physc: pressure more than 10% below '// &
        'lowest standard atmosphere pressure'
    ENDIF

  END SUBROUTINE dissoc_physc

  ! --------------------------------------------------------------------------

  SUBROUTINE dissoc_finish

    IMPLICIT NONE
    INTEGER :: ip

    DEALLOCATE(jtemp)
    DEALLOCATE(jpress)

    DO ip = 1, IP_MAX
      DEALLOCATE(dissoc_2d(ip)%ptr)
      NULLIFY(dissoc_2d(ip)%ptr)
    ENDDO
    DEALLOCATE(dissoc_2d)    ; NULLIFY(dissoc_2d)

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

    CALL close_file(ncid_dissoc)

  END SUBROUTINE dissoc_finish

  ! --------------------------------------------------------------------------

  SUBROUTINE dissoc_result

    USE caaba_mem, ONLY: model_time, percent_done, photol_clev
    USE messy_main_constants_mem, ONLY: FLAGGED_BAD

    IMPLICIT NONE

    REAL(DP), DIMENSION(IP_MAX) :: dissoc_array
    INTEGER :: ip

    IF (percent_done > 0.) THEN
      ! The following loop is necessary because the component (after the
      ! "%") to the right of a part reference with nonzero rank cannot
      ! have the POINTER attribute. In other words, we cannot use
      ! "dissoc_2d(:)%ptr(1,photol_clev,1)"
      DO ip = 1, IP_MAX
        dissoc_array(ip) = dissoc_2d(ip)%ptr(1,photol_clev)
      ENDDO
      CALL write_output_file(ncid_dissoc, model_time, dissoc_array)
    ELSE
      ! photol_clev not defined in first call of dissoc_result (before dissoc_physc)
      CALL write_output_file(ncid_dissoc, model_time, &
        (/ (FLAGGED_BAD, ip=1,IP_MAX) /))
    ENDIF

  END SUBROUTINE dissoc_result

  ! --------------------------------------------------------------------------

END MODULE messy_dissoc_box

!*****************************************************************************
