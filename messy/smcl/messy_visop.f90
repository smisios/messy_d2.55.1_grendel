! *************************************************************************
MODULE messy_visop
! *************************************************************************

  ! VISOP
  ! visible satellite image forward operator
  ! 
  ! 2014 Leonhard Scheck
  ! 2019 Bastian Kern: Split in two SUBROUTINES to use QNI from model

  USE messy_main_constants_mem, ONLY: dp

  IMPLICIT NONE
  PRIVATE

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'visop'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '0.6'

  PUBLIC :: OPTPROP_COLUMN
  PUBLIC :: OPTPROP_COLUMN_QNI

CONTAINS

  SUBROUTINE OPTPROP_COLUMN( nz, z, temp, rho, qc, qi,&
                             tau06w, tau06i, tau08w, tau08i, &
                             reffw, reffi, cwtop, cwbase, citop, cibase )

    INTEGER,                   INTENT(IN)    :: nz
    REAL(dp), DIMENSION(nz+1), INTENT(IN)    :: z
    REAL(dp), DIMENSION(nz),   INTENT(IN)    :: temp, rho, qc, qi
    REAL(dp),                  INTENT(INOUT) :: tau06w, tau06i, tau08w, tau08i, &
                                                reffw, reffi, cwtop, cwbase, citop, cibase

    ! minimum/maximum particle sizes [microns]
    REAL(dp), PARAMETER :: rw_min =  2.5_dp, rw_max = 25.0_dp
    REAL(dp), PARAMETER :: ri_min = 20.0_dp, ri_max = 60.0_dp

    ! droplet concentration, corresponds to qnc, which is set to 2e8 (hard coded) in ICON
    REAL(dp), PARAMETER :: ndroplet = 2.0e8_dp

    REAL(dp) :: LWC, IWC, Rw, Ri, B
    REAL(dp) :: betaw_by_LWC, betai_by_IWC, betaw, betai, dtauw, dtaui
    REAL(dp) :: pi, wavelength
    INTEGER  :: k, band06i, band06w, band08i, band08w

    pi = 4._dp*DATAN(1._dp)

    tau06w  = 0.0_dp
    tau06i  = 0.0_dp
    tau08w  = 0.0_dp
    tau08i  = 0.0_dp
    reffw = 0.0_dp
    reffi = 0.0_dp

    band06i = 0
    band06w = 0
    band08i = 0
    band08w = 0

    DO k = 1, nz
       ! compute liquid and frozen water content
       IWC = 1000._dp * qi(k) * rho(k) ! ddensity   ! g/m3
       LWC = 1000._dp * qc(k) * rho(k) ! density    ! g/m3

       ! compute effective particle sizes

       ! water: Martin, Johnson, Spice (1994)
       Rw = 1.0e+6_dp * ( 0.75_dp*LWC / (pi*0.67_dp* ndroplet *1.0e+6_dp) )**(1._dp/3._dp)
       Rw = max( rw_min, min( Rw, rw_max) )

       ! ice: Wyser (1998) parameterization
       IF( ( IWC > 0 ).and.( temp(k) < 273._dp ) ) THEN
          B = -2.0_dp + log10( IWC/50.0_dp ) * ((273.0_dp-temp(k))**1.5_dp) * 1e-3_dp
          Ri = (4.0_dp/(sqrt(3.0_dp)+4.0_dp)) * ( 377.4_dp + 203.3_dp*B + 37.91_dp*B**2 + 2.3696_dp*B**3 )
       ELSE
          Ri = 0._dp
       END IF
       Ri = max( ri_min, min( Ri, ri_max) )

       ! VIS008 CHANNEL ----------------------------------------------
       wavelength = 0.808_dp

       ! compute extinction coefficients
       CALL interp_betaw_by_LWC( wavelength, Rw, betaw_by_LWC, band08w )
       CALL interp_betai_by_IWC( wavelength, Ri, betai_by_IWC, band08i )
       betaw = LWC * betaw_by_LWC
       betai = IWC * betai_by_IWC

       ! add up optical depths, weighted sizes
       dtauw = betaw * (z(k)-z(k+1))
       dtaui = betai * (z(k)-z(k+1))

       tau08w = tau08w + dtauw
       tau08i = tau08i + dtaui

       ! VIS006 CHANNEL ----------------------------------------------
       wavelength = 0.638_dp

       ! compute extinction coefficients
       CALL interp_betaw_by_LWC( wavelength, Rw, betaw_by_LWC, band06w )
       CALL interp_betai_by_IWC( wavelength, Ri, betai_by_IWC, band06i )
       betaw = LWC * betaw_by_LWC
       betai = IWC * betai_by_IWC

       ! add up optical depths, weighted sizes
       dtauw = betaw * (z(k)-z(k+1))
       dtaui = betai * (z(k)-z(k+1))

       tau06w = tau06w + dtauw
       tau06i = tau06i + dtaui

       ! Use dtau for VIS006 channel as weight in effective radius calculation
       reffw = reffw + Rw*dtauw
       reffi = reffi + Ri*dtaui

       ! detect water cloud top and base
       IF( (tau06w > 0.1_dp ).and.(cwtop < 1e-6_dp) ) THEN
          cwtop  = z(k)
          cwbase = z(k)
       END IF
       IF( (dtauw > 0.01_dp).and.(cwtop > 1e-6_dp) ) cwbase = z(k)

       ! detect ice cloud top and base
       IF( (tau06i > 0.03_dp).and.(citop < 1e-6_dp) ) THEN
          citop  = z(k)
          cibase = z(k)
       END IF
       IF( (dtaui > 0.003_dp).and.(citop > 1e-6_dp) ) cibase = z(k)
    END DO

    ! normalize effective sizes
    IF( tau06w > 1e-5_dp ) THEN
       reffw = reffw / tau06w
    ELSE
       reffw =  rw_min
    END IF
    if( tau06i > 1e-5_dp ) THEN
       reffi = reffi / tau06i
    ELSE
       reffi = ri_min
    END IF

  END SUBROUTINE OPTPROP_COLUMN

  SUBROUTINE OPTPROP_COLUMN_QNI( nz, z, rho, qc, qi, qni, &
                             tau06w, tau06i, tau08w, tau08i, &
                             reffw, reffi, cwtop, cwbase, citop, cibase )

    INTEGER,                   INTENT(IN)    :: nz
    REAL(dp), DIMENSION(nz+1), INTENT(IN)    :: z
    REAL(dp), DIMENSION(nz),   INTENT(IN)    :: rho, qc, qi, qni
    REAL(dp),                  INTENT(INOUT) :: tau06w, tau06i, tau08w, tau08i, &
                                                reffw, reffi, cwtop, cwbase, citop, cibase

    ! minimum/maximum particle sizes [microns]
    REAL(dp), PARAMETER :: rw_min =  2.5_dp, rw_max = 25.0_dp
    REAL(dp), PARAMETER :: ri_min = 20.0_dp, ri_max = 60.0_dp

    ! droplet concentration, corresponds to qnc, which is set to 2e8 (hard coded) in ICON
    REAL(dp), PARAMETER :: ndroplet = 2.0e8_dp

    ! parameters for ice from Seifert & Beheng 2006
    REAL(dp), PARAMETER :: ageo    = 0.217_dp      ! m kg^-bgeo
    REAL(dp), PARAMETER :: bgeo    = 0.302115_dp   ! --
    REAL(dp), PARAMETER :: nugam   = 1._dp         ! --
    REAL(dp), PARAMETER :: mugam   = 1.0_dp/3.0_dp ! --
    REAL(dp), PARAMETER :: rho_ice = 0.9167e+3_dp  ! kg/m3 = 0.9167 g/cm3 


    REAL(dp) :: LWC, IWC, Rw, Ri, D_ge, xbar
    REAL(dp) :: betaw_by_LWC, betai_by_IWC, betaw, betai, dtauw, dtaui
    REAL(dp) :: pi, wavelength
    INTEGER  :: k, band06i, band06w, band08i, band08w

    pi = 4._dp*DATAN(1._dp)

    tau06w  = 0.0_dp
    tau06i  = 0.0_dp
    tau08w  = 0.0_dp
    tau08i  = 0.0_dp
    reffw = 0.0_dp
    reffi = 0.0_dp

    band06i = 0
    band06w = 0
    band08i = 0
    band08w = 0

    DO k = 1, nz
       ! compute liquid and frozen water content
       IWC = 1000._dp * qi(k) * rho(k) ! ddensity   ! g/m3
       LWC = 1000._dp * qc(k) * rho(k) ! density    ! g/m3

       ! compute effective particle sizes

       ! water: Martin, Johnson, Spice (1994)
       Rw = 1.0e+6_dp * ( 0.75_dp*LWC / (pi*0.67_dp* ndroplet *1.0e+6_dp) )**(1._dp/3._dp)
       Rw = max( rw_min, min( Rw, rw_max) )

       ! ice: compute particle size using number concentration from two moment scheme
       IF( qni(k) > 1.0_dp ) THEN
          xbar = qi(k) / qni(k)  ! SB06, Eq. 80
          D_ge = qi(k) * rho(k) &
               / ( sqrt(3._dp*sqrt(3._dp)*rho_ice*ageo/8._dp) &
               * mgam( (1._dp+bgeo)/2._dp, nugam, mugam, xbar, qni(k), rho(k) ) &
               + (sqrt(3._dp)/(4._dp*ageo)) * mgam( 1._dp-bgeo, nugam, mugam, xbar, qni(k), rho(k) ) )
       ELSE
          D_ge = 0.0_dp
       END IF
       Ri = 1e6_dp*(4._dp/(sqrt(3._dp)+4._dp)) * D_ge ! [micron]
       Ri = max( ri_min, min( Ri, ri_max) )

       ! VIS008 CHANNEL ----------------------------------------------
       wavelength = 0.808_dp

       ! compute extinction coefficients
       CALL interp_betaw_by_LWC( wavelength, Rw, betaw_by_LWC, band08w )
       CALL interp_betai_by_IWC( wavelength, Ri, betai_by_IWC, band08i )
       betaw = LWC * betaw_by_LWC
       betai = IWC * betai_by_IWC

       ! add up optical depths, weighted sizes
       dtauw = betaw * (z(k)-z(k+1))
       dtaui = betai * (z(k)-z(k+1))

       tau08w = tau08w + dtauw
       tau08i = tau08i + dtaui

       ! VIS006 CHANNEL ----------------------------------------------
       wavelength = 0.638_dp

       ! compute extinction coefficients
       CALL interp_betaw_by_LWC( wavelength, Rw, betaw_by_LWC, band06w )
       CALL interp_betai_by_IWC( wavelength, Ri, betai_by_IWC, band06i )
       betaw = LWC * betaw_by_LWC
       betai = IWC * betai_by_IWC

       ! add up optical depths, weighted sizes
       dtauw = betaw * (z(k)-z(k+1))
       dtaui = betai * (z(k)-z(k+1))

       tau06w = tau06w + dtauw
       tau06i = tau06i + dtaui

       ! Use dtau for VIS006 channel as weight in effective radius calculation
       reffw = reffw + Rw*dtauw
       reffi = reffi + Ri*dtaui

       ! detect water cloud top and base
       IF( (tau06w > 0.1_dp ).and.(cwtop < 1e-6_dp) ) THEN
          cwtop  = z(k)
          cwbase = z(k)
       END IF
       IF( (dtauw > 0.01_dp).and.(cwtop > 1e-6_dp) ) cwbase = z(k)

       ! detect ice cloud top and base
       IF( (tau06i > 0.03_dp).and.(citop < 1e-6_dp) ) THEN
          citop  = z(k)
          cibase = z(k)
       END IF
       IF( (dtaui > 0.003_dp).and.(citop > 1e-6_dp) ) cibase = z(k)
    END DO

    ! normalize effective sizes
    IF( tau06w > 1e-5_dp ) THEN
       reffw = reffw / tau06w
    ELSE
       reffw =  rw_min
    END IF
    if( tau06i > 1e-5_dp ) THEN
       reffi = reffi / tau06i
    ELSE
       reffi = ri_min
    END IF

  END SUBROUTINE OPTPROP_COLUMN_QNI

  SUBROUTINE interp_betaw_by_LWC( wavelength, Rw, betaw_by_LWC, band )

    ! Water extinction coefficient according to Hu & Stamnes 1992

    REAL(dp), INTENT(IN)             :: wavelength, Rw
    REAL(dp), INTENT(OUT)            :: betaw_by_LWC
    INTEGER, OPTIONAL, INTENT(INOUT) :: band

    INTEGER                              :: k, kband, krbin
    INTEGER, PARAMETER                   :: nbands = 25, nrbins = 3
    REAL(dp)                             :: fband, bl, bh
    REAL(dp), DIMENSION(4,nbands,nrbins) :: tab
    ! from <libradtran_data_path>/wc/wc.ext (only wavelengths <= 3.9 micron)
    !
    ! Parameterisation of extinction coefficient, b_ext = a * r**b  + c
    ! 
    ! lambda(micron)  a          b          c
    DATA tab / 3.900_dp, 6.40E+03_dp, -1.79E+00_dp,  7.03E+01_dp, & ! 2.5 micron < r < 12.5 micron
               3.690_dp, 5.29E+03_dp, -1.73E+00_dp,  7.34E+01_dp, &
               3.145_dp, 2.71E+03_dp, -1.27E+00_dp,  2.35E+01_dp, &
               2.618_dp, 4.56E+03_dp, -1.61E+00_dp,  5.74E+01_dp, &
               2.247_dp, 3.26E+03_dp, -1.46E+00_dp,  5.42E+01_dp, &
               1.855_dp, 2.15E+03_dp, -1.15E+00_dp,  1.42E+01_dp, &
               1.587_dp, 2.01E+03_dp, -1.11E+00_dp,  8.80E+00_dp, &
               1.393_dp, 1.98E+03_dp, -1.11E+00_dp,  7.61E+00_dp, &
               1.232_dp, 1.96E+03_dp, -1.11E+00_dp,  9.29E+00_dp, &
               1.142_dp, 1.94E+03_dp, -1.11E+00_dp,  1.01E+01_dp, &
               1.046_dp, 1.91E+03_dp, -1.10E+00_dp,  7.51E+00_dp, &
               0.929_dp, 1.87E+03_dp, -1.09E+00_dp,  8.41E+00_dp, &
               0.821_dp, 1.86E+03_dp, -1.09E+00_dp,  8.61E+00_dp, &
               0.766_dp, 1.84E+03_dp, -1.09E+00_dp,  8.81E+00_dp, &
               0.719_dp, 1.81E+03_dp, -1.08E+00_dp,  6.85E+00_dp, &
               0.664_dp, 1.79E+03_dp, -1.07E+00_dp,  5.98E+00_dp, &
               0.603_dp, 1.76E+03_dp, -1.06E+00_dp,  5.01E+00_dp, &
               0.544_dp, 1.75E+03_dp, -1.07E+00_dp,  5.95E+00_dp, &
               0.499_dp, 1.73E+03_dp, -1.06E+00_dp,  5.13E+00_dp, &
               0.459_dp, 1.72E+03_dp, -1.06E+00_dp,  4.99E+00_dp, &
               0.419_dp, 1.70E+03_dp, -1.05E+00_dp,  4.49E+00_dp, &
               0.379_dp, 1.68E+03_dp, -1.05E+00_dp,  4.26E+00_dp, &
               0.344_dp, 1.67E+03_dp, -1.04E+00_dp,  3.49E+00_dp, &
               0.314_dp, 1.67E+03_dp, -1.04E+00_dp,  3.83E+00_dp, &
               0.290_dp, 1.63E+03_dp, -1.03E+00_dp,  7.66E-01_dp, &
               3.900_dp, 2.24E+03_dp, -1.11E+00_dp,  3.32E+00_dp, & ! 12.5 micron < r < 30.0 micron
               3.690_dp, 2.17E+03_dp, -1.10E+00_dp,  3.01E+00_dp, &
               3.145_dp, 2.02E+03_dp, -1.08E+00_dp,  2.24E+00_dp, &
               2.618_dp, 2.05E+03_dp, -1.09E+00_dp,  2.66E+00_dp, &
               2.247_dp, 1.99E+03_dp, -1.08E+00_dp,  2.54E+00_dp, &
               1.855_dp, 1.91E+03_dp, -1.07E+00_dp,  1.96E+00_dp, &
               1.587_dp, 1.87E+03_dp, -1.06E+00_dp,  1.93E+00_dp, &
               1.393_dp, 1.83E+03_dp, -1.06E+00_dp,  1.63E+00_dp, &
               1.232_dp, 1.80E+03_dp, -1.05E+00_dp,  1.50E+00_dp, &
               1.142_dp, 1.78E+03_dp, -1.05E+00_dp,  1.32E+00_dp, &
               1.046_dp, 1.77E+03_dp, -1.05E+00_dp,  1.46E+00_dp, &
               0.929_dp, 1.74E+03_dp, -1.04E+00_dp,  1.19E+00_dp, &
               0.821_dp, 1.73E+03_dp, -1.04E+00_dp,  1.16E+00_dp, &
               0.766_dp, 1.71E+03_dp, -1.04E+00_dp,  1.01E+00_dp, &
               0.719_dp, 1.70E+03_dp, -1.04E+00_dp,  1.04E+00_dp, &
               0.664_dp, 1.69E+03_dp, -1.03E+00_dp,  9.89E-01_dp, &
               0.603_dp, 1.68E+03_dp, -1.03E+00_dp,  9.28E-01_dp, &
               0.544_dp, 1.67E+03_dp, -1.03E+00_dp,  8.73E-01_dp, &
               0.499_dp, 1.66E+03_dp, -1.03E+00_dp,  8.13E-01_dp, &
               0.459_dp, 1.65E+03_dp, -1.03E+00_dp,  7.23E-01_dp, &
               0.419_dp, 1.64E+03_dp, -1.02E+00_dp,  6.44E-01_dp, &
               0.379_dp, 1.64E+03_dp, -1.03E+00_dp,  8.33E-01_dp, &
               0.344_dp, 1.62E+03_dp, -1.02E+00_dp,  6.34E-01_dp, &
               0.314_dp, 1.61E+03_dp, -1.02E+00_dp,  5.44E-01_dp, &
               0.290_dp, 1.63E+03_dp, -1.03E+00_dp,  9.90E-01_dp, &
               3.900_dp, 1.20E+03_dp, -8.70E-01_dp, -8.51E+00_dp, & !  30.0 micron < r < 60.0 micron
               3.690_dp, 1.17E+03_dp, -8.64E-01_dp, -8.67E+00_dp, &
               3.145_dp, 1.12E+03_dp, -8.52E-01_dp, -8.99E+00_dp, &
               2.618_dp, 1.12E+03_dp, -8.52E-01_dp, -8.94E+00_dp, &
               2.247_dp, 1.09E+03_dp, -8.46E-01_dp, -9.08E+00_dp, &
               1.855_dp, 1.07E+03_dp, -8.40E-01_dp, -9.22E+00_dp, &
               1.587_dp, 1.05E+03_dp, -8.36E-01_dp, -9.31E+00_dp, &
               1.393_dp, 1.04E+03_dp, -8.32E-01_dp, -9.40E+00_dp, &
               1.232_dp, 1.03E+03_dp, -8.30E-01_dp, -9.44E+00_dp, &
               1.142_dp, 1.01E+03_dp, -8.26E-01_dp, -9.61E+00_dp, &
               1.046_dp, 1.01E+03_dp, -8.24E-01_dp, -9.65E+00_dp, &
               0.929_dp, 9.99E+02_dp, -8.22E-01_dp, -9.69E+00_dp, &
               0.821_dp, 9.91E+02_dp, -8.20E-01_dp, -9.74E+00_dp, &
               0.766_dp, 9.84E+02_dp, -8.18E-01_dp, -9.79E+00_dp, &
               0.719_dp, 9.78E+02_dp, -8.16E-01_dp, -9.89E+00_dp, &
               0.664_dp, 9.76E+02_dp, -8.16E-01_dp, -9.84E+00_dp, &
               0.603_dp, 9.70E+02_dp, -8.14E-01_dp, -9.91E+00_dp, &
               0.544_dp, 9.63E+02_dp, -8.12E-01_dp, -9.98E+00_dp, &
               0.499_dp, 9.62E+02_dp, -8.12E-01_dp, -9.95E+00_dp, &
               0.459_dp, 9.55E+02_dp, -8.10E-01_dp, -1.00E+01_dp, &
               0.419_dp, 9.54E+02_dp, -8.10E-01_dp, -9.99E+00_dp, &
               0.379_dp, 9.48E+02_dp, -8.08E-01_dp, -1.01E+01_dp, &
               0.344_dp, 9.42E+02_dp, -8.06E-01_dp, -1.02E+01_dp, &
               0.314_dp, 9.41E+02_dp, -8.06E-01_dp, -1.01E+01_dp, &
               0.290_dp, 9.40E+02_dp, -8.06E-01_dp, -1.01E+01_dp  /

    kband = 0
    IF (PRESENT(band)) THEN
       kband = band
    END IF

    IF (kband == 0) THEN
       ! find band index
       DO k = 1, nbands-1
          IF( (wavelength <= tab(1,k,1)).and.(wavelength >= tab(1,k+1,1)) ) THEN
             kband = k
             EXIT
          END IF
       END DO
       IF (PRESENT(band)) THEN
          band = kband
       END IF
    END IF

    krbin = 0
    IF      ((Rw >=  2.5_dp).and.(Rw <  12.5_dp)) THEN
       krbin = 1
    ELSE IF ((Rw >= 12.5_dp).and.(Rw <  30.0_dp)) THEN
       krbin = 2
    ELSE IF ((Rw >= 30.0_dp).and.(Rw <= 60.0_dp)) THEN
       krbin = 3
    END IF

    IF( krbin == 0 ) THEN
       krbin = 1
       print *, 'krbin error in interp_betaw_by_LWC, Rw=', RW
    END IF
    IF( kband == 0 ) THEN
       kband = 1
       print *, 'kband error in interp_betaw_by_LWC, wavelength=', wavelength
    END IF

    !betaw_by_LWC = (tab(2,kband,krbin) * Rw**tab(3,kband,krbin) + tab(4,kband,krbin))*1e-3_dp

    fband = (wavelength-tab(1,kband,1)) / (tab(1,kband+1,1)-tab(1,kband,1))
    bl = (tab(2,kband  ,krbin) * Rw**tab(3,kband  ,krbin) + tab(4,kband  ,krbin))*1e-3_dp
    bh = (tab(2,kband+1,krbin) * Rw**tab(3,kband+1,krbin) + tab(4,kband+1,krbin))*1e-3_dp

    betaw_by_LWC = bl*(1.0_dp-fband) + bh*fband

  END SUBROUTINE interp_betaw_by_LWC


  SUBROUTINE interp_betai_by_IWC( wavelength, Ri, betai_by_IWC, band )

    ! Ice extinction coefficient according to Fu 1996

    REAL(dp),          INTENT(IN)    :: wavelength, Ri
    REAL(dp),          INTENT(OUT)   :: betai_by_IWC
    INTEGER, OPTIONAL, INTENT(INOUT) :: band

    INTEGER                       :: k, kband
    INTEGER, PARAMETER            :: nbands = 25
    REAL(dp)                      :: fband, bl, bh
    REAL(dp), DIMENSION(4,nbands) :: tab
    !
    ! from <libradtran_data_path>/ic/fu96/fu96.ext
    !
    !# Q. Fu, An Accurate Parameterization of the Solar Radiative 
    !# Properties of Cirrus Clouds for Climate Models, Journal of 
    !# Climate 9, 2058-2082, 1996.
    !#
    !# Extinction coefficient b = IWC * (a0 + a1/D), eq. 3.9a
    !#
    !#        Band limit     a0             a1
    !#        (um)
    !#
    DATA tab /0.25_dp, 0.30_dp,  -0.236447e-03_dp,  0.253817e+01_dp, &
              0.30_dp, 0.33_dp,  -0.266955e-03_dp,  0.254179e+01_dp, &
              0.33_dp, 0.36_dp,  -0.293599e-03_dp,  0.254540e+01_dp, &
              0.36_dp, 0.40_dp,  -0.258858e-03_dp,  0.253815e+01_dp, &
              0.40_dp, 0.44_dp,  -0.106451e-03_dp,  0.252684e+01_dp, &
              0.44_dp, 0.48_dp,   0.129121e-03_dp,  0.250410e+01_dp, &
              0.48_dp, 0.52_dp,  -0.945458e-04_dp,  0.252061e+01_dp, &
              0.52_dp, 0.57_dp,  -0.303108e-04_dp,  0.251805e+01_dp, &
              0.57_dp, 0.64_dp,   0.982244e-04_dp,  0.250875e+01_dp, &
              0.64_dp, 0.69_dp,   0.161983e-03_dp,  0.250746e+01_dp, &
              0.69_dp, 0.75_dp,  -0.304991e-03_dp,  0.254412e+01_dp, &
              0.75_dp, 0.78_dp,   0.226539e-03_dp,  0.249909e+01_dp, &
              0.78_dp, 0.87_dp,   0.810443e-04_dp,  0.251619e+01_dp, &
              0.87_dp, 1.00_dp,   0.737638e-04_dp,  0.251051e+01_dp, &
              1.00_dp, 1.10_dp,  -0.614288e-03_dp,  0.256520e+01_dp, &
              1.10_dp, 1.19_dp,   0.413595e-03_dp,  0.248783e+01_dp, &
              1.19_dp, 1.41_dp,   0.651659e-04_dp,  0.251660e+01_dp, &
              1.41_dp, 1.53_dp,  -0.805155e-03_dp,  0.257600e+01_dp, &
              1.53_dp, 1.64_dp,   0.644675e-03_dp,  0.247060e+01_dp, &
              1.64_dp, 2.13_dp,  -0.837325e-04_dp,  0.252504e+01_dp, &
              2.13_dp, 2.38_dp,   0.489477e-03_dp,  0.248776e+01_dp, &
              2.38_dp, 2.91_dp,   0.234245e-03_dp,  0.248573e+01_dp, &
              2.91_dp, 3.42_dp,   0.297295e-03_dp,  0.248895e+01_dp, &
              3.42_dp, 4.00_dp,   0.187598e-03_dp,  0.251396e+01_dp, &
              4.00_dp, 4.99_dp,  -0.254823e-03_dp,  0.252909e+01_dp  /

    kband = 0
    IF (PRESENT(band)) THEN
       kband = band
    END IF

    IF (kband == 0) THEN
       ! find band index
       DO k = 1, nbands
          IF( (wavelength >= tab(1,k)).and.(wavelength <= tab(2,k)) ) THEN
             kband = k
             EXIT
          END IF
       END DO
       IF (PRESENT(band)) THEN
          band = kband
       END IF
    END IF

    IF( kband == 0 ) THEN
       kband = 1
       print *, 'kband error in interp_betai_by_IWC, wavelength=', wavelength
    END IF

    !betai_by_IWC = ( tab(3,kband) + tab(4,kband)/(2.0*Ri) ) * 1.299

    fband = (wavelength-tab(1,kband)) / (tab(2,kband)-tab(1,kband))
    bl = ( tab(3,kband  ) + tab(4,kband  )/(2._dp*Ri) ) * 1.299_dp
    bh = ( tab(3,kband+1) + tab(4,kband+1)/(2._dp*Ri) ) * 1.299_dp

    betai_by_IWC = bl*(1._dp-fband) + bh*fband

  END SUBROUTINE interp_betai_by_IWC

  FUNCTION gammln(xx)
    ! logarithm of gamma function from Numerical Recipes

    implicit none
    REAL(dp), intent(in) :: xx
    REAL(dp)             :: gammln

    integer :: j
    REAL(dp) :: ser, stp, tmp, x, y, cof(6)

    data cof, stp / 76.18009172947146_dp, -86.50532032941677_dp, 24.01409824083091_dp, &
         -1.231739572450155_dp, .1208650973866179E-2_dp, -.5395239384953E-5_dp, 2.5066282746310005_dp /

    x = xx
    y = x
    tmp = x + 5.5_dp
    tmp = (x + 0.5_dp)*log(tmp) - tmp
    ser = 1.000000000190015_dp
    do j = 1, 6
       y   = y + 1._dp
       ser = ser + cof(j)/y
    end do
    gammln = tmp + log(stp*ser/x)
    return

  END FUNCTION gammln

  FUNCTION mgam( n, nu, mu, xbar, qni, ddensity )
    ! moment definition from Seifert & Beheng 2006, Eq. 82

    implicit none
    REAL(dp), intent(in) :: n, nu, mu, xbar, qni, ddensity
    REAL(dp)             :: mgam

    mgam = exp( gammln((n+nu+1._dp)/mu) - gammln((nu+1._dp)/mu) ) &
         * ( max(xbar,1E-11_dp) * exp( gammln((nu+1._dp)/mu) - gammln((nu+2._dp)/mu)) )**n &
         * (qni*ddensity)
    return

  END FUNCTION mgam


! *************************************************************************
END MODULE messy_visop
! *************************************************************************
