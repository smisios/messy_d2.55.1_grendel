!*****************************************************************************
!                Time-stamp: <2016-09-27 19:04:29 sander>
!*****************************************************************************

MODULE cloudj_column

  USE messy_main_constants_mem, ONLY: DP
  USE messy_cmn_photol_mem, ONLY: IP_MAX, jname
  USE messy_cloudj ! init_cloudj, cloudj_2d, cloudj_read_nml_ctrl, cloudjvalues, lp


  IMPLICIT NONE
  PRIVATE

  ! PUBLIC cloudj INTERFACE ROUTINES
  PUBLIC :: cloudj_initialize
  PUBLIC :: cloudj_init_memory
  PUBLIC :: cloudj_physc
  PUBLIC :: cloudj_result
  PUBLIC :: cloudj_free_memory

  INTEGER, PARAMETER :: nsza = 10
  INTEGER, PARAMETER :: nlev = 19 ! number of levels
  REAL, DIMENSION(nsza) :: sza = &
    (/ 0., 10., 20., 30., 40., 50., 60., 70., 80., 90. /)

CONTAINS

  ! --------------------------------------------------------------------------

  SUBROUTINE cloudj_initialize

    IMPLICIT NONE

    ! LOCAL
    INTEGER :: status ! status flag

    ! INTITIALIZE GLOBAL SWITCHES / PARAMETERS
    CALL cloudj_read_nml_ctrl(status, 99)
    IF (status /= 0) STOP

    lp(:) = .TRUE.

  END SUBROUTINE cloudj_initialize

  ! --------------------------------------------------------------------------

  SUBROUTINE cloudj_init_memory

    IMPLICIT NONE

    INTEGER :: jt

    ALLOCATE(cloudj_2d(ip_MAX))
    DO jt=1, ip_MAX
      ALLOCATE(cloudj_2d(jt)%ptr(nsza,nlev))
    END DO

    CALL init_cloudj

  END SUBROUTINE cloudj_init_memory

  ! --------------------------------------------------------------------------

  SUBROUTINE cloudj_physc

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

    DO i=1,nsza
      ! global average values are extracted with ferret from messy and
      ! cloudj_diag streams using e.g.: "list rhum[i=@ave,j=@ave,l=1]"

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

    ! calculate cloudjvalues
    CALL cloudjvalues(                      &
      v3,                                   &
      COS(sza*REAL(pi)/180.), press, relo3, &
      rhum, temp, albedo,                   &
      aclc, slf, clp)

  END SUBROUTINE cloudj_physc

  ! --------------------------------------------------------------------------

  SUBROUTINE cloudj_result

    USE mo_netcdf, ONLY: open_cloudj_nc_file,  &
      write_cloudj_nc_file, &
      close_cloudj_nc_file

    INTEGER :: ncid_cloudj ! netcdf id for cloudj.nc
    INTEGER :: jt

    CALL open_cloudj_nc_file(ncid_cloudj, nlev, nsza, sza)

    DO jt=1, ip_MAX
      CALL write_cloudj_nc_file(ncid_cloudj, 'J_'//TRIM(jname(jt)), &
        cloudj_2d(jt)%ptr(:,:))
    ENDDO

    CALL close_cloudj_nc_file(ncid_cloudj)

  END SUBROUTINE cloudj_result

  ! --------------------------------------------------------------------------

  SUBROUTINE cloudj_free_memory

    IMPLICIT NONE

    INTEGER :: jt

    DO jt=1, ip_MAX
      DEALLOCATE(cloudj_2d(jt)%ptr)
      NULLIFY(cloudj_2d(jt)%ptr)
    ENDDO
    DEALLOCATE(cloudj_2d)
    NULLIFY(cloudj_2d)

  END SUBROUTINE cloudj_free_memory

  ! --------------------------------------------------------------------------

END MODULE cloudj_column

!*****************************************************************************

PROGRAM cloudj

  USE cloudj_column, ONLY: cloudj_init_memory, cloudj_initialize, &
    cloudj_physc, cloudj_result, cloudj_free_memory

  IMPLICIT NONE

  CALL cloudj_initialize   ! read CTRL namelist, intialize aerosol
  CALL cloudj_init_memory
  CALL cloudj_physc        ! calculate J values
  CALL cloudj_result       ! print results
  CALL cloudj_free_memory

END PROGRAM cloudj

!*****************************************************************************
