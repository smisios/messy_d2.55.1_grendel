!*****************************************************************************
!                Time-stamp: <2020-06-19 13:08:23 b302010>
!*****************************************************************************

! SAPPHO = Simplified And Parameterized PHOtolysis rates

! Authors:
! Rolf Sander,    MPICH, Mainz, 2003-2007
! Hella Riede,    MPICH, Mainz, 2007

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

MODULE messy_sappho

  USE messy_main_constants_mem, ONLY: DP, PI
  USE messy_cmn_photol_mem      ! IP_MAX, ip_*, jname

  IMPLICIT NONE
  PRIVATE
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'sappho'
  PUBLIC :: jvalues

  ! photolysis rate coefficients:
  REAL(DP), PUBLIC, DIMENSION(IP_MAX), SAVE :: jx

  ! **************************************************************************

CONTAINS

  SUBROUTINE jvalues(cossza, photo_scenario, rad_lat, photon, efact)

    IMPLICIT NONE

    ! I/O
    REAL(DP),          INTENT(IN)  :: cossza
    CHARACTER(LEN=12), INTENT(IN)  :: photo_scenario
    REAL(DP),          INTENT(IN)  :: rad_lat
    REAL(DP),          INTENT(OUT) :: photon
    REAL(DP),          INTENT(IN)  :: efact

    ! The value of DUSK was chosen so that photon=5E-2 when cossza=0,
    ! i.e. at sunrise and sunset.
    REAL(DP), PARAMETER :: DUSK = 0.0721347_DP ! 5.E-2/LOG(2.)
    REAL(DP) :: DAY8090, FCT, dn, Lichtfaktor

    ! photon is approximately the positive part of cossza but
    !   - avoid sharp switching on and off
    !   - add some light at dawn before sunrise and at dusk after sunset
    photon = dusk*LOG(1._DP+EXP(cossza/dusk))
    ! dn is a day/night switch to set photolyses to about 0 at night
    fct=50._DP
    dn = EXP(fct*cossza)/(EXP(-fct*cossza)+EXP(fct*cossza))

    SELECT CASE (TRIM(photo_scenario))
    CASE ('CUMULUS')   
      CALL photo_cumulus
    CASE ('FF_ANTARCTIC','FF_ARCTIC')
      CALL photo_ff
    CASE ('LAB')
      CALL photo_lab
    CASE ('LAB_C15')
      jx(:) = 0. ! dark
    CASE ('','MBL','MOM','OOMPH','VOLCANO')
      CALL photo_mbl
    CASE DEFAULT
      PRINT *, 'ERROR, photo_scenario '//TRIM(photo_scenario)// &
        ' is not defined in '//TRIM(modstr)//'.'
      STOP 1
    END SELECT

    !-------------------------------------------------------------------------

  CONTAINS

    !-------------------------------------------------------------------------

    SUBROUTINE photo_ff
      ! quick and dirty conversion of J values by Landgraf et al.
      ! DAY8090 = ratio between first day of spring (JD80) and 1 April (JD90);
      DAY8090 = 0.674_DP;
      ! photon/COS(rad_lat) reaches about 1 on noon of first day of spring;
      ! PI/2. converts from 12h-mean to noon maximum;

      jx(ip_O1D)    = 6.4697E-04 / EXP( 2.8200    /( 1.1388E-01+photon)) ! J01
      jx(ip_O3P)    = 7.2132E-04 / EXP( 1.7075    /( 1.9411E-01+photon)) ! J02
      jx(ip_H2O2)   = 5.3105E-05 / EXP( 1.3588    /( 1.8467E-01+photon)) ! J03
      jx(ip_NO2)    = 3.9861E-02 / EXP( 7.4565E-01/( 1.2747E-01+photon)) ! J04
      jx(ip_NOO2)   = 2.2000E-02 / EXP( 1.0000E-02/( 1.0000E-03+photon)) ! J05
      jx(ip_NO2O)   = 1.8000E-01 / EXP( 1.0000E-02/( 1.0000E-03+photon)) ! J06
      jx(ip_N2O5)   = 2.1621E-04 / EXP( 1.2521    /( 1.8418E-01+photon)) ! J07
      jx(ip_HNO3)   = 4.7367E-06 / EXP( 2.1397    /( 1.9219E-01+photon)) ! J08
      jx(ip_CH3OOH) = 4.7946E-05 / EXP( 1.3607    /( 1.8532E-01+photon)) ! J09
      jx(ip_CHOH)   = 2.1913E-04 / EXP( 1.4250    /( 1.7789E-01+photon)) ! J10
      jx(ip_COH2)   = 2.7881E-04 / EXP( 1.1487    /( 1.6675E-01+photon)) ! J11
      jx(ip_HOCl)   = 1.8428E-03 / EXP( 1.2575    /( 1.8110E-01+photon)) ! J12
      jx(ip_Cl2O2)  = 1.8302E-02 / EXP( 9.9560E-01/( 1.5859E-01+photon)) ! J13
      jx(ip_ClNO2)  = 2.4861E-03 / EXP( 1.2157    /( 1.7917E-01+photon)) ! J14
      jx(ip_ClNO3)  = 2.1778E-04 / EXP( 9.4755E-01/( 1.5296E-01+photon)) ! J15
      jx(ip_Cl2)    = 1.3150E-02 / EXP( 9.1210E-01/( 1.4929E-01+photon)) ! J16
      jx(ip_HOBr)   = 4.0575E-03 / EXP( 8.9913E-01/( 1.5034E-01+photon)) ! J17
      jx(ip_BrNO2)  = 5.0826E-02 / EXP( 6.5344E-01/( 1.1219E-01+photon)) ! J18
      jx(ip_BrNO3)  = 6.0737E-03 / EXP( 7.5901E-01/( 1.2785E-01+photon)) ! J19
      jx(ip_Br2)    = 6.9038E-02 / EXP( 4.9563E-01/( 8.3580E-02+photon)) ! J20
      jx(ip_BrCl)   = 3.4235E-02 / EXP( 6.3421E-01/( 1.0965E-01+photon)) ! J21
      jx(ip_IO)     = 1.2E-01 *PI/2. *DAY8090* photon/COS(rad_lat)       ! J22
      jx(ip_HOI)    = 2.0E-03 *PI/2. *DAY8090* photon/COS(rad_lat)       ! J23
      jx(ip_INO3)   = 1.4E-03 *PI/2. *DAY8090* photon/COS(rad_lat)       ! J24
      jx(ip_CH3I)   = 9.1E-07 *PI/2. *DAY8090* photon/COS(rad_lat)       ! J25
      jx(ip_I2)     = 6.1E-02 *PI/2. *DAY8090* photon/COS(rad_lat)       ! J26
      jx(ip_BrO)    = 2.3319E-01 / EXP( 1.1094    /( 1.6736E-01+photon)) ! J28
      jx(ip_ICl)    = 9.3E-03 *PI/2. *DAY8090* photon/COS(rad_lat)       ! J29
      jx(ip_IBr)    = 2.7E-02 *PI/2. *DAY8090* photon/COS(rad_lat)       ! J30
      jx(ip_INO2)   = jx(ip_INO3) ! assumed by R.V.                      ! J31
      jx(ip_C3H7I)  = 2.8E-06 *PI/2. *DAY8090* photon/COS(rad_lat)       ! J32
      jx(ip_CH2ClI) = 2.9E-05 *PI/2. *DAY8090* photon/COS(rad_lat)       ! J33
      jx(ip_CH2I2)  = 2.1E-03 *PI/2. *DAY8090* photon/COS(rad_lat)       ! J34
      jx(ip_OClO)   = 3.2837E-01 / EXP( 6.9733E-01/( 1.1907E-01+photon)) ! J89
      jx(ip_HNO4)   = 6.0657E-05 / EXP( 1.5901    /( 1.8577E-01+photon)) ! J91
      jx(ip_HONO)   = 1.0218E-02 / EXP( 8.7261E-01/( 1.4850E-01+photon)) ! J92
      jx(ip_CH3Br)  = 0.                                                 ! J99
      jx(ip_CH3CHO) = jx(ip_CHOH)+jx(ip_COH2) ! mz_rs_20050911 assumed
    END SUBROUTINE photo_ff

    !-------------------------------------------------------------------------

    SUBROUTINE photo_lab
      ! email S. Bleicher from 2013-12-04: Die Lampen haben im neuen
      ! Zustand eine etwa vierfache Intensität im Vergleich mit Lampen
      ! die mehr als 600 h gebrannt haben. Mit dem Lichtfaktor skaliere
      ! ich alle j linear zum abgebrannten Zustand der Lampen. Es sind
      ! insgesamt sieben Lampen vorhanden deren Alter gemischt ist. Ich
      ! führe zwar ein Buch über das Alter der Lampen, der Lichtfaktor=2
      ! hat sich jedoch als Mittelwert etabliert.
      Lichtfaktor = 2
      jx(ip_CH3CHO)   = 1.98547E-6     * Lichtfaktor
      jx(ip_CH3COCH3) = 6.3392E-06     * Lichtfaktor
      jx(ip_CH3OOH)   = 1.300389E-05   * Lichtfaktor
      jx(ip_CHOH)     = 8.030251E-05   * Lichtfaktor
      jx(ip_COH2)     = 0.0001400505   * Lichtfaktor
      jx(ip_H2O2)     = 1.461451E-05   * Lichtfaktor
      jx(ip_HNO3)     = 1.232162E-06   * Lichtfaktor
      jx(ip_HNO4)     = 1.255586E-05   * Lichtfaktor
      jx(ip_HONO)     = 0.00108        * Lichtfaktor
      jx(ip_N2O5)     = 1.94451E-5     * Lichtfaktor
      jx(ip_NO2)      = 0.00522        * Lichtfaktor
      jx(ip_NO2O)     = 0.07025        * Lichtfaktor
      jx(ip_NOO2)     = 0.00838        * Lichtfaktor
      jx(ip_O1D)      = 2.01887E-5     * Lichtfaktor
      jx(ip_O3P)      = 1.59546E-4     * Lichtfaktor
      jx(ip_PAN)      = 9.793765E-07   * Lichtfaktor
      jx(ip_OClO)    = 0.038           * Lichtfaktor
      jx(ip_NO3NOO)  = 1.05727E-8      * Lichtfaktor
      jx(ip_HNO4)    = 1.68178E-6      * Lichtfaktor
      jx(ip_Cl2)     = 0.0012          * Lichtfaktor
      jx(ip_Cl2O2)   = 6.40443E-4      * Lichtfaktor
      jx(ip_HCl)     = 0.              * Lichtfaktor
      jx(ip_HOCl)    = 1.32871E-4      * Lichtfaktor
      !jx(ip_HOCl)   = 10.32871E-3     * Lichtfaktor
      jx(ip_ClNO2)   = 1.51193E-4      * Lichtfaktor
      jx(ip_ClNO3)   = 2.09834E-5      * Lichtfaktor
      jx(ip_ClONO2)  = 3.08242E-6      * Lichtfaktor
      jx(ip_Br2)     = 0.01294         * Lichtfaktor
      jx(ip_BrO)     = 0.01777         * Lichtfaktor
      jx(ip_HOBr)    = 0.00113         * Lichtfaktor
      !jx(ip_HOBr)   = 0.0613          * Lichtfaktor
      jx(ip_BrNO2)   = 0.00288         * Lichtfaktor
      jx(ip_BrNO3)   = 5.49697E-4      * Lichtfaktor
      jx(ip_CH3Br)   = 0.              * Lichtfaktor
      jx(ip_CH2Br2)  = 0.              * Lichtfaktor
      jx(ip_CHBr3)   = 2.89048E-7      * Lichtfaktor
      jx(ip_BrCl)    = 0.0054          * Lichtfaktor
      jx(ip_CH2ClBr) = 1.14484E-9      * Lichtfaktor
      jx(ip_CHCl2Br) = 1.44636E-7      * Lichtfaktor
      jx(ip_CHClBr2) = 1.52024E-7      * Lichtfaktor
      jx(ip_I2)      = 0.05595         * Lichtfaktor
      jx(ip_IO)      = 0.08634         * Lichtfaktor
      jx(ip_HOI)     = 0.0047          * Lichtfaktor
      jx(ip_INO2)    = 0.00166         * Lichtfaktor
      jx(ip_INO3)    = 0.02617         * Lichtfaktor
      jx(ip_CH2I2)   = 0.00321         * Lichtfaktor
      jx(ip_CH3I)    = 2.73648E-6      * Lichtfaktor
      jx(ip_C3H7I)   = 3.80498E-6      * Lichtfaktor
      jx(ip_CH2ClI)  = 5.79601E-5      * Lichtfaktor
      jx(ip_ICl)     = 0.00829         * Lichtfaktor
      jx(ip_IBr)     = 0.02314         * Lichtfaktor
      jx(ip_H2O2)    = 2.68726E-6      * Lichtfaktor
      jx(ip_HNO3)    = 2.66045E-7      * Lichtfaktor
      jx(ip_COH2)    = 1.60143E-5      * Lichtfaktor
    END SUBROUTINE photo_lab

    !-------------------------------------------------------------------------
      
    ! from Klaus-Peter Heue cumulus cloud DOAS simulation
    ! CARIBIC flight 353 above Guadeloupe
    SUBROUTINE photo_cumulus
      ! from TUV model simulation Klaus-Peter Heue
      ! enhancement factor, EF = 2 already included (cloud white box approx.,
      ! light path down to ground and back up); [Brasseur2002]: EF up to 5
      JX(ip_COH2)     = efact * 1.10E-04
      JX(ip_CHOH)     = efact * 1.20E-04
      JX(ip_CH3COCH3) = efact * 6.70E-06
      JX(ip_HONO)     = efact * 0.00424
      JX(ip_O1D)      = efact * 1.79E-04
      JX(ip_NO2)      = efact * 0.025
      !! from photo_lab
      JX(ip_CH3CHO)   = efact * 5.339559E-05
      !JX(ip_CH3COCH3) = efact * 6.3392E-06   ! = 95% of Heue, 275nm
      JX(ip_CH3OOH)   = efact * 1.300389E-05
      !JX(ip_CHOH)     = efact * 8.030251E-05 ! = 67% of Heue, 300nm max quantum yield
      !JX(ip_COH2)     = efact * 0.0001400505 ! = 127% of Heue, 335nm max q. yield
      JX(ip_H2O2)     = efact * 1.461451E-05
      JX(ip_HNO3)     = efact * 1.232162E-06
      JX(ip_HNO4)     = efact * 1.255586E-05
      !JX(ip_HONO)     = efact * 0.00411277   ! = 97% of Heue, 210nm, 355nm
      JX(ip_N2O5)     = efact * 4.650903E-05
      !JX(ip_NO2)      = efact * 0.01843288 ! = 73% of Heue
      JX(ip_NO2O)     = efact * 0.2858042
      JX(ip_NOO2)     = efact * 0.03579563
      !JX(ip_O1D)      = efact * 7.111431E-05 ! = 40% of Heue, 255nm
      JX(ip_O3P)      = efact * 0.0007458425
      JX(ip_PAN)      = efact * 9.793765E-07
      !!! no Cl chemistry recorded by DOAS, CH3Cl tau >= 3 months
      !! from photo_mbl; factor 2.5 since equiv. values only 40% of photo_lab
      !! comment: equivalent values at constant SZA = 34.46 degrees
      !!jx(ip_O1D)    = dn*4.4916E-04*EXP(-3.5807E+00 /(3.2382E-01+photon)) !1.987210E-05, 28% photo_lab, 11% Heue
      !!jx(ip_O3P)    = dn*4.3917E-04*EXP(-1.6665E-01 /(5.0375E-02+photon)) !3.630040E-04, 49% photo_lab
      !!jx(ip_H2O2)   = dn*1.8182E-05*EXP(-1.2565E+00 /(1.7904E-01+photon)) !5.198816E-06, 36% photo_lab
      !!jx(ip_NO2)    = dn*1.2516E-02*EXP(-4.5619E-01 /(6.9151E-02+photon)) !7.512433E-03, 41% photo_lab
      !!jx(ip_NOO2)   = dn*2.2959E-02*EXP(-1.0873E-01 /(2.1442E-02+photon)) !2.018999E-02, 56% photo_lab
      !!jx(ip_NO2O)   = dn*1.9176E-01*EXP(-1.2557E-01 /(1.9375E-02+photon)) !0.16524882  , 58% photo_lab
      !!jx(ip_N2O5)   = dn*1.1375E-04*EXP(-1.1092E+00 /(1.7004E-01+photon)) !3.729180E-05, 80% photo_lab
      !!jx(ip_HNO3)   = dn*3.4503E-06*EXP(-2.1412E+00 /(2.6143E-01+photon)) !4.803654E-07, 39% photo_lab
      !!jx(ip_CH3OOH) = dn*1.2858E-05*EXP(-1.1739E+00 /(1.7044E-01+photon)) !3.951753E-06, 30% photo_lab
      !!jx(ip_CHOH)   = dn*8.6933E-05*EXP(-1.3293E+00 /(1.6405E-01+photon)) !2.265801E-05, 28% photo_lab, 19% Heue
      !!jx(ip_COH2)   = dn*8.9896E-05*EXP(-8.4745E-01 /(1.2056E-01+photon)) !3.667108E-05, 26% photo_lab, 33% Heue
      !jx(ip_HOCl)   = dn*4.9596E-04*EXP(-8.4303E-01 /(1.3121E-01+photon))*2.5 !2.052947E-04*2.5
      !jx(ip_Cl2O2)  = dn*2.8605E-03*EXP(-8.8027E-01 /(1.4813E-01+photon))*2.5 !1.157202E-03*2.5
      !jx(ip_ClNO2)  = dn*8.0989E-04*EXP(-9.7319E-01 /(1.4656E-01+photon))*2.5 !2.973054E-04*2.5
      !jx(ip_ClNO3)  = dn*7.4610E-05*EXP(-7.2496E-01 /(1.2847E-01+photon))*2.5 !3.486835E-05*2.5
      !jx(ip_Cl2)    = dn*3.9398E-03*EXP(-6.3121E-01 /(1.0305E-01+photon))*2.5 !1.995024E-03*2.5
      !jx(ip_OClO)   = dn*1.0933E-01*EXP(-4.4797E-01 /(7.2470E-02+photon))*2.5 !6.635205E-02*2.5
      !!jx(ip_HNO4)   = dn*2.1532E-05*EXP(-1.9648E+00 /(2.1976E-01+photon)) !3.280883E-06, 26% photo_lab,
      !!jx(ip_HONO)   = dn*2.9165E-03*EXP(-5.1317E-01 /(7.4940E-02+photon)) !1.648514E-03, 40% photo_lab, 40% Heue
    END SUBROUTINE photo_cumulus

    !-------------------------------------------------------------------------

    SUBROUTINE photo_mbl
      ! J values from PAPER model by Landgraf et al.
      ! J value as a function of photon: J = A*exp(-B/(C+photon))
      ! photon=sin(PSI)=cos(THETA)
      ! THETA=zenith angle, PSI=90-THETA
      ! temp profile of atmos. is: data/prof.AFGL.midl.sum
      ! surface albedo :  0.00
      ! ozone column (Dobson) :300.00
      ! cloud cover: OFF
      ! paper model was started on: 00-01-18
      jx(ip_O1D)    = dn*4.4916E-04*EXP(-3.5807E+00 /(3.2382E-01+photon)) !J01
      jx(ip_O3P)    = dn*4.3917E-04*EXP(-1.6665E-01 /(5.0375E-02+photon)) !J02
      jx(ip_H2O2)   = dn*1.8182E-05*EXP(-1.2565E+00 /(1.7904E-01+photon)) !J03
      jx(ip_NO2)    = dn*1.2516E-02*EXP(-4.5619E-01 /(6.9151E-02+photon)) !J04
      jx(ip_NOO2)   = dn*2.2959E-02*EXP(-1.0873E-01 /(2.1442E-02+photon)) !J05
      jx(ip_NO2O)   = dn*1.9176E-01*EXP(-1.2557E-01 /(1.9375E-02+photon)) !J06
      jx(ip_N2O5)   = dn*1.1375E-04*EXP(-1.1092E+00 /(1.7004E-01+photon)) !J07
      jx(ip_HNO3)   = dn*3.4503E-06*EXP(-2.1412E+00 /(2.6143E-01+photon)) !J08
      jx(ip_CH3OOH) = dn*1.2858E-05*EXP(-1.1739E+00 /(1.7044E-01+photon)) !J09
      jx(ip_CHOH)   = dn*8.6933E-05*EXP(-1.3293E+00 /(1.6405E-01+photon)) !J10
      jx(ip_MGLYOX) = 2.*jx(ip_CHOH) ! based on code from Domenico
      jx(ip_COH2)   = dn*8.9896E-05*EXP(-8.4745E-01 /(1.2056E-01+photon)) !J11
      jx(ip_HOCl)   = dn*4.9596E-04*EXP(-8.4303E-01 /(1.3121E-01+photon)) !J12
      jx(ip_Cl2O2)  = dn*2.8605E-03*EXP(-8.8027E-01 /(1.4813E-01+photon)) !J13
      jx(ip_ClNO2)  = dn*8.0989E-04*EXP(-9.7319E-01 /(1.4656E-01+photon)) !J14
      jx(ip_ClNO3)  = dn*7.4610E-05*EXP(-7.2496E-01 /(1.2847E-01+photon)) !J15
      jx(ip_Cl2)    = dn*3.9398E-03*EXP(-6.3121E-01 /(1.0305E-01+photon)) !J16
      jx(ip_HOBr)   = dn*3.2114E-03*EXP(-4.6104E-01 /(8.9870E-02+photon)) !J17
      jx(ip_BrNO2)  = dn*7.6992E-03*EXP(-3.4144E-01 /(6.2243E-02+photon)) !J18
      jx(ip_BrNO3)  = dn*1.9670E-03*EXP(-5.2431E-01 /(9.9522E-02+photon)) !J19
      jx(ip_Br2)    = dn*3.6899E-02*EXP(-2.1097E-01 /(3.3855E-02+photon)) !J20
      jx(ip_BrCl)   = dn*1.3289E-02*EXP(-3.0734E-01 /(5.3134E-02+photon)) !J21
      jx(ip_IO)     = dn*3.6977E-01*EXP(-2.6523E-01 /(3.7541E-02+photon)) !J22
      jx(ip_HOI)    = dn*1.2681E-02*EXP(-3.7515E-01 /(6.1696E-02+photon)) !J23
      jx(ip_INO3)   = dn*5.7985E-03*EXP(-5.2214E-01 /(1.0423E-01+photon)) !J24
      jx(ip_CH3I)   = dn*2.6536E-05*EXP(-1.7937E+00 /(2.4260E-01+photon)) !J25
      jx(ip_I2)     = dn*1.6544E-01*EXP(-1.3475E-01 /(2.1859E-02+photon)) !J26
      jx(ip_BrO)    = dn*6.9639E-02*EXP(-7.4569E-01 /(1.0859E-01+photon)) !J28
      jx(ip_ICl)    = dn*2.6074E-02*EXP(-1.9215E-01 /(3.0177E-02+photon)) !J29
      jx(ip_IBr)    = dn*7.4918E-02*EXP(-1.5931E-01 /(2.3807E-02+photon)) !J30
      jx(ip_INO2)   = dn*5.5440E-03*EXP(-6.6060E-01 /(1.0034E-01+photon)) !J31
      jx(ip_C3H7I)  = dn*8.9493E-05*EXP(-1.8899E+00 /(2.5853E-01+photon)) !J32
      jx(ip_CH2ClI) = dn*4.9715E-04*EXP(-1.4267E+00 /(2.0610E-01+photon)) !J33
      jx(ip_CH2I2)  = dn*2.0184E-02*EXP(-1.0349E+00 /(1.4961E-01+photon)) !J34
      jx(ip_OClO)   = dn*1.0933E-01*EXP(-4.4797E-01 /(7.2470E-02+photon)) !J89
      jx(ip_HNO4)   = dn*2.1532E-05*EXP(-1.9648E+00 /(2.1976E-01+photon)) !J91
      jx(ip_HONO)   = dn*2.9165E-03*EXP(-5.1317E-01 /(7.4940E-02+photon)) !J92
      jx(ip_CH3Br)  = dn*7.3959E-10*EXP(-2.9326E+01 /(1.0132E-01+photon)) !J99
    END SUBROUTINE photo_mbl

    !-------------------------------------------------------------------------

  END SUBROUTINE jvalues

  ! **************************************************************************

END MODULE messy_sappho
