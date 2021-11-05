! ****************************************************************************
!                Time-stamp: <2017-02-01 10:19:47 joec_pa>
!     Author: Rolf Sander, Max-Planck Institute, Mainz, Germany, 2009-2014
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

! ****************************************************************************

MODULE jvpp_step2

  USE messy_main_constants_mem, ONLY: DP
  USE messy_cmn_photol_mem ! IP_MAX, ip_*, jname
  USE jvpp_mem, ONLY: IO_IN, IO_OUT, IO_TDP, IO_TMP, HLINE2, DIR_176, DIR_EFF, &
                      DIR_JNL, NWAV, NBIN, l_sig, l_s2t, l_s3t, l_phi, &
                      l_tc1, l_tc2, l_tc3, l_species, l_tdep, wave_nm, inputdir

  IMPLICIT NONE

  REAL(DP), PARAMETER :: GRAV  = 9.81_dp      ! acceleration due to gravity [m/s^2]
  REAL(DP), PARAMETER :: RELO2 = 0.2095_dp    ! O2 mixing ratio
  REAL(DP), PARAMETER :: AVOGA = 6.022E23_dp  ! Avogadro constant [1/mol]
  REAL(DP), PARAMETER :: AIR_M = 28.97E-3_dp  ! mol mass of air [kg/mol]
  REAL(DP), PARAMETER :: CONST = RELO2*AVOGA/(AIR_M*GRAV)*1.E-2_dp ! [1/(cm2*hPa)]
  REAL(DP), PARAMETER :: PRESS   = 1000.000_dp  ! pressure [hPa]
  REAL(DP), PARAMETER :: BOLTZ   = 1.381E-23_dp ! Boltzmann constant [J/K]
  REAL(DP), PARAMETER :: PLANCK  = 6.626E-34_dp ! Planck constant [J*s]
  REAL(DP), PARAMETER :: C_LIGHT = 2.9979E8_dp  ! speed of light [m/s]

  INTEGER :: N_tdep(IP_MAX) = 1 ! N_tdep = 2 or 3
  CHARACTER(LEN=80) :: filename(IP_MAX), filename_tdep(IP_MAX)

  INTEGER :: li ! initial interval number
  INTEGER :: lf ! final interval number
  INTEGER :: ji
  CHARACTER(LEN=6) :: suffix

  ! cross sections cs_xxx_tdep [cm^2] at temperatures T_xxx [K]
  REAL(DP) :: cs_xxx_tdep(NWAV,IP_MAX,3)
  REAL(DP) :: t_xxx(IP_MAX,3)

  ! miscellaneous, wavelength-dependent temperature coefficients:
  REAL(DP), DIMENSION(NWAV,IP_MAX,3) :: tc_xxx = 0._dp

  ! cross sections [cm^2/part.]:
  REAL(DP), DIMENSION(NWAV,IP_MAX) :: cs_xxx, cst_xxx

  ! Chebyshev coefficients A and B:
  REAL(DP) :: cheb_coeff_a(20,13), cheb_coeff_b(20,13)

  ! quantum yields phi [1]:
  REAL(DP) :: phi_xxx(IP_MAX, NWAV) = 1._dp

  ! extraterrestrial flux from flux.ATLAS2 [photons/(cm^2 s)]:
  REAL(DP), PARAMETER :: flx(NWAV) = (/ &
    2.41E+11_dp, 3.14E+11_dp, 3.56E+11_dp, 3.46E+11_dp, 4.22E+11_dp, 5.37E+11_dp, 6.45E+11_dp, &
    7.41E+11_dp, 7.36E+11_dp, 1.03E+12_dp, 1.22E+12_dp, 1.32E+12_dp, 1.63E+12_dp, 1.91E+12_dp, &
    2.34E+12_dp, 2.79E+12_dp, 5.03E+12_dp, 8.11E+12_dp, 9.27E+12_dp, 8.92E+12_dp, 1.16E+13_dp, &
    1.22E+13_dp, 1.68E+13_dp, 1.38E+13_dp, 1.46E+13_dp, 1.67E+13_dp, 1.40E+13_dp, 1.70E+13_dp, &
    1.50E+13_dp, 2.33E+13_dp, 2.08E+13_dp, 2.10E+13_dp, 2.03E+13_dp, 2.85E+13_dp, 5.21E+13_dp, &
    4.33E+13_dp, 1.12E+14_dp, 1.23E+14_dp, 1.15E+14_dp, 1.06E+14_dp, 7.56E+13_dp, 1.56E+14_dp, &
    2.21E+14_dp, 3.64E+14_dp, 2.79E+14_dp, 1.37E+14_dp, 7.07E+13_dp, 6.87E+13_dp, 7.43E+13_dp, &
    8.76E+13_dp, 9.69E+13_dp, 9.52E+13_dp, 9.17E+13_dp, 9.37E+13_dp, 9.84E+13_dp, 8.92E+13_dp, &
    9.04E+13_dp, 1.10E+14_dp, 1.13E+14_dp, 1.11E+14_dp, 1.12E+14_dp, 1.06E+14_dp, 1.04E+14_dp, &
    1.16E+14_dp, 1.19E+14_dp, 1.14E+14_dp, 1.26E+14_dp, 1.57E+14_dp, 2.01E+14_dp, 2.59E+14_dp, &
    4.43E+14_dp, 7.18E+14_dp, 7.91E+14_dp, 8.38E+14_dp, 8.20E+14_dp, 8.61E+14_dp, 9.16E+14_dp, &
    8.32E+14_dp, 1.04E+15_dp, 1.11E+15_dp, 1.00E+15_dp, 1.13E+15_dp, 9.15E+14_dp, 1.17E+15_dp, &
    9.78E+14_dp, 1.62E+15_dp, 1.75E+15_dp, 1.84E+15_dp, 1.87E+15_dp, 1.95E+15_dp, 1.81E+15_dp, &
    1.67E+15_dp, 1.98E+15_dp, 2.02E+15_dp, 2.18E+15_dp, 2.36E+15_dp, 2.31E+15_dp, 2.39E+15_dp, &
    2.38E+15_dp, 2.39E+15_dp, 2.44E+15_dp, 2.51E+15_dp, 2.30E+15_dp, 2.39E+15_dp, 2.48E+15_dp, &
    2.40E+15_dp, 2.46E+15_dp, 2.49E+15_dp, 2.32E+15_dp, 2.39E+15_dp, 2.42E+15_dp, 2.55E+15_dp, &
    2.51E+15_dp, 2.49E+15_dp, 2.55E+15_dp, 2.53E+15_dp, 2.54E+15_dp, 2.50E+15_dp, 2.57E+15_dp, &
    2.58E+15_dp, 2.67E+15_dp, 2.67E+15_dp, 2.70E+15_dp, 2.62E+15_dp, 2.69E+15_dp, 2.63E+15_dp, &
    2.68E+15_dp, 2.66E+15_dp, 2.59E+15_dp, 2.69E+15_dp, 2.61E+15_dp, 2.62E+15_dp, 2.62E+15_dp, &
    2.63E+15_dp, 2.60E+15_dp, 2.55E+15_dp, 2.48E+15_dp, 2.57E+15_dp, 2.61E+15_dp, 2.61E+15_dp, &
    2.62E+15_dp, 2.62E+15_dp, 2.57E+15_dp, 2.52E+15_dp, 2.60E+15_dp, 2.58E+15_dp, 2.52E+15_dp, &
    2.51E+15_dp, 2.48E+15_dp, 2.45E+15_dp, 2.48E+15_dp, 2.45E+15_dp, 2.44E+15_dp, 2.39E+15_dp, &
    2.40E+15_dp, 2.41E+15_dp, 2.40E+15_dp, 2.38E+15_dp, 2.34E+15_dp, 2.32E+15_dp, 2.30E+15_dp, &
    2.33E+15_dp, 2.34E+15_dp, 2.29E+15_dp, 2.29E+15_dp, 2.27E+15_dp, 2.27E+15_dp, 2.20E+15_dp, &
    2.22E+15_dp, 2.18E+15_dp, 2.20E+15_dp, 2.14E+15_dp, 2.14E+15_dp, 2.13E+15_dp, 2.09E+15_dp, &
    2.05E+15_dp  /)

  REAL(DP) :: flx_sum    ! integrated extraterrestrial flux [photons/(cm^2 s)]
  REAL(DP) :: T_ref      ! reference temperature [K]
  REAL(DP) :: v2min, v2max, v3min_DU, v3max_DU
  INTEGER :: maxv3, maxv2
  LOGICAL :: l_tdep_loop

CONTAINS

  ! --------------------------------------------------------------------------

  SUBROUTINE read_tc(jp, filetype)

    IMPLICIT NONE
    INTEGER,          INTENT(IN)  :: jp
    CHARACTER(LEN=*), INTENT(IN)  :: filetype
    INTEGER :: zN_tdep, jw

    zN_tdep = 1
    IF (filetype=='tc2') zN_tdep = 2
    IF (filetype=='tc3') zN_tdep = 3

    ! read wavelength-dependent temperature coefficients:
    OPEN(IO_IN, FILE=DIR_176//'/'//TRIM(jname(jp))//'.'//filetype//'_176', &
      STATUS='OLD')
    DO jw = 1,NWAV
      READ(IO_IN,*) tc_xxx(jw,jp,1:zN_tdep)
    ENDDO
    CLOSE(IO_IN)
    l_tdep(jp) = .TRUE.

  END SUBROUTINE read_tc

  ! --------------------------------------------------------------------------

  SUBROUTINE read_phi(jp)

    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: jp
    INTEGER :: jw

    OPEN(IO_IN, FILE=DIR_176//'/'//TRIM(jname(jp))//'.phi_176', STATUS='OLD')
    DO jw = 1,NWAV
      READ(IO_IN,*) phi_xxx(jp,jw)
    ENDDO
    CLOSE(IO_IN)

  END SUBROUTINE read_phi

  ! --------------------------------------------------------------------------

  SUBROUTINE initialize

    IMPLICIT NONE

    INTEGER :: jp, i, jw

    DO jp = 1, IP_MAX ! loop over all photolysis reactions

      ! check if tdep file exists:
      IF (l_s2t(jp)) filename_tdep(jp) = &
        DIR_176//'/'//TRIM(jname(jp))//'.s2t_176'
      IF (l_s3t(jp)) filename_tdep(jp) = &
        DIR_176//'/'//TRIM(jname(jp))//'.s3t_176'
      IF (l_s2t(jp).OR.l_s3t(jp)) THEN
        OPEN(IO_IN, FILE=TRIM(filename_tdep(jp)), STATUS='OLD')
        READ(IO_IN,*) N_tdep(jp)
        READ(IO_IN,*) (T_XXX(jp,i),i=1,N_tdep(jp))
        DO jw = 1,NWAV
          READ(IO_IN,*) cs_xxx_tdep(jw,jp,1:N_tdep(jp))
        ENDDO
        CLOSE(IO_IN)
        l_tdep(jp) = .TRUE.
      ENDIF

      ! check if file for one temperature exists:
      IF (l_sig(jp)) THEN
        filename(jp) = DIR_176//'/'//TRIM(jname(jp))//'.sig_176'
        OPEN(IO_IN, FILE=TRIM(filename(jp)), STATUS='OLD')
        DO jw = 1,NWAV
          READ(IO_IN,*) cs_xxx(jw, jp)
        ENDDO
        CLOSE(IO_IN)
      ENDIF

      l_species(jp) = l_sig(jp) .OR. l_s2t(jp) .OR. l_s3t(jp)
      ! N2O does not have an input file but is calculated analytically later:
      l_species(ip_n2o) = .TRUE.

      ! Tdep using temperature coefficients:
      IF (l_tc1(jp)) CALL read_tc(jp, 'tc1')
      IF (l_tc2(jp)) CALL read_tc(jp, 'tc2')
      IF (l_tc3(jp)) CALL read_tc(jp, 'tc3')

      ! quantum yields phi:
      IF (l_phi(jp)) CALL read_phi(jp)

    ENDDO

    ! Chebyshev coefficient A and B for Koppers & Murtagh parameterization:
    OPEN(IO_IN,file='cheb_coeff.txt',status='old')
    DO I=1,20
      READ(IO_IN,*) cheb_coeff_a(i,:)
    ENDDO
    DO I=1,20
      READ(IO_IN,*) cheb_coeff_b(i,:)
    ENDDO
    CLOSE(IO_IN)

  END SUBROUTINE initialize

  ! --------------------------------------------------------------------------

  SUBROUTINE calc_tdep

    ! calculation of temperature-dependent cross sections

    IMPLICIT NONE

    REAL(DP), DIMENSION(NWAV) :: c2, c3
    REAL(DP) :: chi, suma, sumb, T_clip, T_red, wave_nm_clip, xa, xb, factor, B
    INTEGER :: i, jp, jw
    ! for H2O2 (JPL2011):
    REAL(DP), PARAMETER :: A_H2O2(0:7) = (/                 &
         6.4761000E+04_dp, -9.2170972E+02_dp,    4.5356490E+00_dp,  &
         -4.4589016E-03_dp, -4.0351010E-05_dp,    1.6878206E-07_dp, &
         -2.6520140E-10_dp,  1.5534675E-13_dp/)
    REAL(DP), PARAMETER :: B_H2O2(0:4) = (/                 &
         6.8123000E+03_dp, -5.1351000E+01_dp,    1.1522000E-01_dp,  &
         -3.0493000E-05_dp, -1.0924000E-07_dp/)
    ! for CH3COCH3:
    REAL(DP), PARAMETER :: A_CH3COCH3(6) = (/ 13.1852_dp, -0.252309_dp, &
         0.00192398_dp, -7.30847e-6_dp, 1.3832e-8_dp, -1.04375e-11_dp /)
    REAL(DP), PARAMETER :: B_CH3COCH3(6) = (/ -0.0445238_dp, 0.000851385_dp, &
         -6.48766e-6_dp, 2.46273e-8_dp, -4.65788e-11_dp, 3.51253e-14_dp /)
    ! for CH3Cl:
    REAL(DP), PARAMETER :: A_CH3Cl(0:4) = (/                 &
         -299.80_dp, 5.1047_dp, -3.3630E-2_dp, 9.5805E-5_dp, -1.0135E-7_dp /)
    REAL(DP), PARAMETER :: B_CH3Cl(0:4) = (/                 &
         -7.1727_dp, 1.4837E-1_dp, -1.1463E-3_dp, 3.9188E-6_dp, -4.9994E-9_dp /)
    ! for CH3CCl3:
    REAL(DP), PARAMETER :: A_CH3CCl3(0:4) = (/                 &
         341.085191_dp, -7.273362_dp, 5.498387E-2_dp, -1.827552E-4_dp, &
         2.238640E-7_dp /)
    REAL(DP), PARAMETER :: B_CH3CCl3(0:4) = (/ & 
         -1.660090_dp, 3.079969E-2_dp, -2.106719E-4_dp, 6.264984E-7_dp, &
         -6.781342E-10_dp /)
    ! op_dc_20161027+
    ! for CHF2Cl (JPL2011)
    REAL(DP), PARAMETER :: A_CHF2Cl(0:3) = (/                  &
         -103.029_dp, 1.5038_dp, -8.2476e-3_dp, 1.4206e-5_dp/)
    REAL(DP), PARAMETER :: B_CHF2Cl(0:3) = (/                  &
         -1.3399e-1_dp, 2.7405e-3_dp, -1.8028e-5_dp, 3.8504e-8_dp/)
    ! for CH3CFCl2 (JPL2011)
    REAL(DP), PARAMETER :: A_CH3CFCl2(0:4) = (/                 &
         -682.913042_dp, 12.122290_dp, -8.187699e-2_dp, 2.437244e-4_dp, &
         -2.719103e-7_dp /)
    REAL(DP), PARAMETER :: B_CH3CFCl2(0:4) = (/                 &
         4.074747_dp, -8.053899e-2_dp, 5.946552e-4_dp, -1.945048e-6_dp, &
         2.380143e-9_dp /)
    ! for CF2ClCF2Cl (JPL2011;CFC-114)
    REAL(DP), PARAMETER :: A_CF2ClCF2Cl(0:4) = (/                 &
         -160.50_dp, 2.4807_dp, -1.5202e-2_dp, 3.8412e-5_dp, -3.4373e-8_dp /)
    REAL(DP), PARAMETER :: B_CF2ClCF2Cl(0:4) = (/                 &
         -1.5296_dp, 3.5248e-2_dp, -2.9951e-4_dp, 1.1129e-6_dp, -1.5259e-9_dp/)
    ! for CF2ClCFCl2 (JPL2011,CFC-113)
    REAL(DP), PARAMETER :: A_CF2ClCFCl2(0:4) = (/                 &
         -1087.9_dp, 20.004_dp, -1.3920e-1_dp, 4.2828e-4_dp, -4.9384e-7_dp /)
    REAL(DP), PARAMETER :: B_CF2ClCFCl2(0:4) = (/                 &
         12.493_dp, -2.3937e-1_dp, 1.7142e-3_dp, -5.4393e-6_dp, 6.4548e-9_dp/)
    ! for CF3CF2Cl (JPL2011,CFC-115)
    REAL(DP), PARAMETER :: A_CF3CF2Cl (0:3) = (/                   &
         5.8281_dp, -2.990e-2_dp, 1.3525e-3_dp, -2.6851e-6_dp /)
    ! for CH2Cl2 (JPL 2011)
    REAL(DP), PARAMETER :: A_CH2Cl2 (0:4) = (/                   &
         -1421.8_dp, 27.395_dp, -1.9807e-1_dp, 6.3468e-4_dp, -7.6298e-7_dp /)
    REAL(DP), PARAMETER :: B_CH2Cl2 (0:4) = (/                   &
         -3.1171_dp, 6.7874e-2_dp, -5.5000e-4_dp, 1.9649e-6_dp, -2.6101e-9_dp/)
    ! for CHCl3 (JPL 2011)
    REAL(DP), PARAMETER :: A_CHCl3 (0:4) = (/                   &
         269.80_dp, -6.0908_dp, 4.7830e-2_dp, -1.6427e-4_dp, 2.0682e-7_dp/)
    REAL(DP), PARAMETER :: B_CHCl3 (0:4) = (/                   &
         3.7973_dp, -7.0913e-2_dp, 4.9397e-4_dp, -1.5226e-6_dp, 1.7555e-9_dp/)
    ! op_dc_20161027-
    CHARACTER(LEN=80) :: filesuffix

    WRITE(filesuffix, '(A,I3.3,A,I1.1,A)') '_', NINT(T_ref), 'K_', ji, '.tmp'
    IF (l_tdep_loop) filesuffix = '_Tdep'//TRIM(filesuffix)

    ! first, set cst_xxx to cs_xxx:
    cst_xxx(:,:) = cs_xxx(:,:)

    ! Tdep via cs_xxx_tdep:
    DO jp = 1, IP_MAX ! loop over all photolysis reactions
      IF (N_tdep(jp)>=2) THEN
        c2(:) = (cs_xxx_tdep(:,jp,2)-cs_xxx_tdep(:,jp,1)) / &
          (T_xxx(jp,2)-T_xxx(jp,1))
        IF (N_tdep(jp)>=3) THEN
          c3(:) = ((cs_xxx_tdep(:,jp,3) - cs_xxx_tdep(:,jp,2)) / &
            (T_xxx(jp,3)-T_xxx(jp,2)) - c2(:)) / (T_xxx(jp,3)-T_xxx(jp,1))
        ELSE
          c3(:) = 0._dp
        ENDIF
        cst_xxx(:,jp) = MAX(0._dp, ((T_ref-T_xxx(jp,2))*c3(:) + &
          c2(:)) * (T_ref-T_xxx(jp,1)) + cs_xxx_tdep(:,jp,1))
      ENDIF
    ENDDO

    ! Tdep via wavelength-dependent temperature coefficients tc:
    ! 9) HNO3 (JPL2011, p. 282):
    IF (l_tc1(ip_hno3))     cst_xxx(:,ip_hno3)       = & ! tc1
      cs_xxx(:,ip_hno3)     * EXP(tc_xxx(:,ip_hno3,1)  * (T_ref-298._dp))
    ! 11) PAN (ln sigma = ln sigma(298K) + B(T-298), JPL2011, p. 310):
    IF (l_tc1(ip_pan))      cst_xxx(:,ip_pan)        = & ! tc1
      cs_xxx(:,ip_pan)      * EXP(tc_xxx(:,ip_pan,1)   * (T_ref-298._dp))
    ! 14) COH2 (HCHO):
    IF (l_tc1(ip_coh2))     cst_xxx(:,ip_coh2)       = & ! tc1
      cs_xxx(:,ip_coh2)     + tc_xxx(:,ip_coh2,1)      * (T_ref-298._dp)
    ! 15) CHOH (HCHO):
    IF (l_tc1(ip_choh))     cst_xxx(:,ip_choh)       = & ! tc1
      cs_xxx(:,ip_choh)     + tc_xxx(:,ip_choh,1)      * (T_ref-298._dp)
    ! currently, hardcoded CH3COCH3 is used
    ! ! 18) CH3COCH3 (JPL2011, p. 318):
    ! IF (l_tc3(ip_ch3coch3)) cst_xxx(:,ip_ch3coch3)   = & ! tc3
    !   cs_xxx(:,ip_ch3coch3) * (1. + tc_xxx(:,ip_ch3coch3,1)*T_ref + &
    !   tc_xxx(:,ip_ch3coch3,2)*T_ref**2 + tc_xxx(:,ip_ch3coch3,3)*T_ref**3)
    ! 23) ClNO3 (JPL2011, p. 389):
    IF (l_tc2(ip_clno3))    cst_xxx(:,ip_clno3)      = & ! tc2
      cs_xxx(:,ip_clno3)    * (1._dp + tc_xxx(:,ip_clno3,1)*(T_ref-296._dp) + &
      tc_xxx(:,ip_clno3,2)*(T_ref-296._dp)**2)
    ! 24) ClNO2 (equation (VII) of Ghosh et al.):
    IF (l_tc2(ip_clno2))    cst_xxx(:,ip_clno2)      = & ! tc2
      cs_xxx(:,ip_clno2)    * (1._dp + tc_xxx(:,ip_clno2,1)*(T_ref-296._dp) + &
      tc_xxx(:,ip_clno2,2)*(T_ref-296._dp)**2)
    ! 29) BrNO3 (JPL2011, p. 453):
    IF (l_tc2(ip_brno3))    cst_xxx(:,ip_brno3)      = & ! tc2
      cs_xxx(:,ip_brno3)    * (1._dp + tc_xxx(:,ip_brno3,1)*(T_ref-296._dp) + &
      tc_xxx(:,ip_brno3,2)*(T_ref-296._dp)**2)
    ! 41) C3H7I (sigma = sigma(298)*(1+a1(T-298)+a2(T-298)^2), JPL2011, p. 507):
    IF (l_tc2(ip_c3h7i))    cst_xxx(:,ip_c3h7i)      = & ! tc2
      cs_xxx(:,ip_c3h7i)    * (1._dp + tc_xxx(:,ip_c3h7i,1)*(T_ref-298._dp) + &
      tc_xxx(:,ip_c3h7i,2)*(T_ref-298._dp)**2)
    ! 63) CH2Br2 (JPL2011, p. 458):
    IF (l_tc1(ip_ch2br2))   cst_xxx(:,ip_ch2br2)     = & ! tc1
      cs_xxx(:,ip_ch2br2)   * EXP(tc_xxx(:,ip_ch2br2,1) * (T_ref-298._dp))
    ! 67) ClONO2 (JPL2011, p. 389):
    IF (l_tc2(ip_clono2))    cst_xxx(:,ip_clono2)    = & ! tc2
      cs_xxx(:,ip_clono2)    * (1._dp + tc_xxx(:,ip_clono2,1)*(T_ref-296._dp)+ &
      tc_xxx(:,ip_clono2,2)*(T_ref-296._dp)**2)
    !104) CH3NO3 (sigma = sigma(298K) * exp(B(T-298)), JPL2011, p. 309):
    ! (note that the tc1 file contains 1E3*B, not B)
    IF (l_tc1(ip_ch3no3))    cst_xxx(:,ip_ch3no3)    = & ! tc1
      cs_xxx(:,ip_ch3no3)    * EXP(tc_xxx(:,ip_ch3no3,1)*1E-3_dp*(T_ref-298._dp))

    ! Tdep via individual wavelength-dependent functions:

    ! 1) O2
    l_tdep(ip_o2) = .TRUE.
    ! Koppers & Murtagh parameterization for Schumann-Runge band from
    ! 179.37 nm <= lambda <= 201.01 nm (index 1-13):
    CALL sr_o2_km(CONST*PRESS,cst_xxx(1:13,ip_o2))
    ! spectrum from input file for larger wavelengths:
    cst_xxx(14:NWAV,ip_o2) = cs_xxx(14:NWAV,ip_o2)
    ! mz_rs_20130204+
    ! JPL-2011 only has data for >= 205 nm
    ! data from Koppers & Murtagh is only used for <= 201.01 nm (index 1-13)
    ! simple fix for gap in interpolation:
    cst_xxx(14,ip_o2) = cst_xxx(13,ip_o2) + &
      (cst_xxx(16,ip_o2)-cst_xxx(13,ip_o2)) * &
      (wave_nm(14)-wave_nm(13)) / (wave_nm(16)-wave_nm(13))
    cst_xxx(15,ip_o2) = cst_xxx(13,ip_o2) + &
      (cst_xxx(16,ip_o2)-cst_xxx(13,ip_o2)) * &
      (wave_nm(15)-wave_nm(13)) / (wave_nm(16)-wave_nm(13))
    ! DO jw = 1,NWAV
    !   PRINT *, jw, wave_nm(jw), cst_xxx(jw,ip_o2)
    ! ENDDO
    ! mz_rs_20130204-

    ! 4) H2O2 (JPL2011)
    l_tdep(ip_h2o2) = .TRUE.
    DO jw = 1,NWAV
      T_clip = MIN(400._dp,MAX(200._dp,T_ref))
      chi = 1._dp/(1.+EXP(-1265._dp/T_clip))
      IF (jw<=35) THEN
        cst_xxx(jw,ip_h2o2) = cs_xxx(jw,ip_h2o2)
      ENDIF
      IF ((jw>=36).AND.(jw<=76)) THEN ! 261.44 nm <= lambda <= 350 nm
        xa=0._dp
        DO i = 0,7
          xa = A_H2O2(i) * wave_nm(jw)**i + xa
        ENDDO
        xb=0._dp
        DO i = 0,4
          xb = B_H2O2(i) * wave_nm(jw)**i + xb
        ENDDO
        cst_xxx(jw,ip_h2o2) = ( chi*xa + (1._dp-chi)*xb ) * 1.E-21_dp
      ENDIF
      IF (jw>76) THEN
        cst_xxx(jw,ip_h2o2) = cs_xxx(jw,ip_h2o2)
      ENDIF
    ENDDO

    ! 6,7) NO3 (NO2O and NOO2)
    l_tdep(ip_no2o) = .TRUE.
    l_tdep(ip_noo2) = .TRUE.
    SELECT CASE(inputdir)
    CASE ("dat_m17")
      ! T-dep via *.s2t input files
    CASE ("dat_lit")
      factor = (1._dp - EXP(-1096.4_dp/T_ref) - 2._dp*EXP(-529.5_dp/T_ref)) / &
        (1._dp - EXP(-1096.4_dp/298.0_dp) - 2._dp*EXP(-529.5_dp/298.0_dp))
      cst_xxx(:,ip_no2o) = cs_xxx(:,ip_no2o) * factor
      cst_xxx(:,ip_noo2) = cs_xxx(:,ip_noo2) * factor
    CASE DEFAULT ! not ok:
      STOP "ERROR: inputdir must be dat_m17 or dat_lit"
    END SELECT

    ! 8) N2O5
    l_tdep(ip_n2o5) = .TRUE.
    SELECT CASE(inputdir)
    CASE ("dat_m17") ! JPL1997
      DO jw = 1,NWAV
        cst_xxx(jw,ip_n2o5) = 0._dp
        IF (jw<=42) THEN ! lambda <= 283.69 nm
          cst_xxx(jw,ip_n2o5) = cs_xxx(jw,ip_n2o5)
        ENDIF
        IF ((jw>=43).AND.(jw<=82)) THEN ! 287.77 nm <= lambda <= 380 nm
          cst_xxx(jw,ip_n2o5) = &
            EXP(2.735_dp+(4728.5_dp-17.127_dp*wave_nm(jw))/T_ref)*1.E-20_dp
        ENDIF
      ENDDO
    CASE ("dat_lit") ! JPL2011:
      DO jw = 1,NWAV
        IF ((jw>=36).AND.(jw<=88)) THEN ! 261.44 nm <= lambda <= 410 nm
          ! mz_rs_20130301+
          ! A = tc_xxx(jw,ip_n2o5,1)
          B = tc_xxx(jw,ip_n2o5,2)
          ! WRONG: cst_xxx(jw,ip_n2o5) = 10.**(A + 1000.*B/T_ref)
          ! Bug fix: Do not use A because it has been extrapolated
          ! towards zero at the edges of the wavelength range (jw=36 and
          ! jw=88) which causes BIG errors. Instead, use the value
          ! sigma_300K and multiply with a correction factor based on B:
          ! lg(sigma) = A + 1000 * B/T
          ! => sigma = sigma_300K * 10 ^ (1000 * B * (1/T-1/300K))
          cst_xxx(jw,ip_n2o5) = &
            cs_xxx(jw,ip_n2o5) * 10.**(1000._dp*B*(1._dp/T_ref-1._dp/300._dp))
          ! mz_rs_20130301-
        ELSE
          cst_xxx(jw,ip_n2o5) = cs_xxx(jw,ip_n2o5)
        ENDIF
      ENDDO
    CASE DEFAULT ! not ok:
      STOP "ERROR: inputdir must be dat_m17 or dat_lit"
    END SELECT
    ! activate to save data to a file:
    ! IF (.NOT.l_tdep_loop) THEN
    !   OPEN(IO_TMP, FILE='N2O5'//filesuffix, STATUS='UNKNOWN')
    !   DO jw = 1,NWAV
    !     WRITE(IO_TMP, '(F8.2,1PE16.8)') wave_nm(jw), cst_xxx(jw,ip_n2o5)
    !   ENDDO
    !   CLOSE(IO_TMP)
    ! ENDIF

    ! 18) CH3COCH3
    SELECT CASE(inputdir)
    CASE ("dat_m17")
      l_tdep(ip_ch3coch3) = .TRUE.
      T_clip = T_ref
      IF (T_clip < 235._dp) T_clip = 235._dp
      IF (T_clip > 298._dp) T_clip = 298._dp
      T_red = 298._dp - T_clip
      DO jw = 1,NWAV
        wave_nm_clip = wave_nm(jw)
        IF (wave_nm_clip < 215._dp) wave_nm_clip = 215._dp
        IF (wave_nm_clip > 330._dp) wave_nm_clip = 330._dp
        suma=0._dp
        sumb=0._dp
        DO i = 1,6
          suma = suma + a_ch3coch3(i)*wave_nm_clip**(i-1)
          sumb = sumb + b_ch3coch3(i)*wave_nm_clip**(i-1)
        ENDDO
        cst_xxx(jw,ip_ch3coch3) = &
          cs_xxx(jw,ip_ch3coch3) * (1._dp+suma*T_red+sumb*T_red*T_red)
      ENDDO
    CASE ("dat_lit")
      ! hardcoded is used here
    CASE DEFAULT ! not ok:
      STOP "ERROR: inputdir must be dat_m17 or dat_lit"
    END SELECT

    ! 33) CH3Cl
    ! pp. 395-396 of JPL2011
    l_tdep(ip_ch3cl) = .TRUE.
    T_clip = T_ref
    IF (T_clip < 210._dp) T_clip = 210._dp
    IF (T_clip > 300._dp) T_clip = 300._dp
    DO jw = 1,19 ! 176 nm <= lambda <= 216 nm
      suma = 0._dp
      sumb = 0._dp
      DO i = 0,4
        suma = suma + A_ch3cl(i) * wave_nm(jw)**i
        sumb = sumb + B_ch3cl(i) * wave_nm(jw)**i
      ENDDO
      cst_xxx(jw,ip_ch3cl) = 10._dp**(suma+(T_clip-273._dp)*sumb)
    ENDDO
    DO jw = 20,NWAV ! lambda > 216 nm
      cst_xxx(jw,ip_ch3cl) = cs_xxx(jw,ip_ch3cl)
    ENDDO

    ! 34) CH3CCl3
    ! pp. 396-397 of JPL2011
    l_tdep(ip_ch3ccl3) = .TRUE.
    DO jw = 1,2 ! lambda < 182 nm
      cst_xxx(jw,ip_ch3ccl3) = cs_xxx(jw,ip_ch3ccl3)
    ENDDO
    T_clip = T_ref
    IF (T_clip < 210._dp) T_clip = 210._dp
    IF (T_clip > 300._dp) T_clip = 300._dp
    DO jw = 3,29 ! 182 nm <= lambda <= 240 nm
      suma = 0._dp
      sumb = 0._dp
      DO i = 0,4
        suma = suma + A_ch3ccl3(i) * wave_nm(jw)**i
        sumb = sumb + B_ch3ccl3(i) * wave_nm(jw)**i
      ENDDO
      cst_xxx(jw,ip_ch3ccl3) = 10._dp**(suma+(T_clip-273._dp)*sumb)
    ENDDO
    DO jw = 30,NWAV ! lambda > 240 nm
      cst_xxx(jw,ip_ch3ccl3) = cs_xxx(jw,ip_ch3ccl3)
    ENDDO

    ! 35) CFCl3
    ! p. 188 of JPL1997
    l_tdep(ip_cfcl3) = .TRUE.
    DO jw = 1,NWAV
      cst_xxx(jw,ip_cfcl3) = 0._dp
      IF (jw<=36) THEN ! lambda <= 261.44 nm
        cst_xxx(jw,ip_cfcl3) = &
          cs_xxx(jw,ip_cfcl3)  * EXP(1.0E-4_dp*(wave_nm(jw)-184.9_dp)*(T_ref-298._dp))
      ENDIF
    ENDDO

    ! 36) CF2Cl2
    ! p. 188 of JPL1997
    l_tdep(ip_cf2cl2) = .TRUE.
    DO jw = 1,NWAV
      cst_xxx(jw,ip_cf2cl2) = 0._dp
      IF (jw<=36) THEN ! lambda <= 261.44 nm
        cst_xxx(jw,ip_cf2cl2) = &
          cs_xxx(jw,ip_cf2cl2) * EXP(4.1E-4_dp*(wave_nm(jw)-184.9_dp)*(T_ref-298._dp))
      ENDIF
    ENDDO

    ! 56) N2O
    l_tdep(ip_n2o) = .TRUE.
    DO jw = 1,NWAV
      cst_xxx(jw,ip_n2o) = 0._dp
      IF (jw<=29) THEN ! lambda <= 239.52 nm
        cst_xxx(jw,ip_n2o) = EXP(6.821023E+01_dp + &
          wave_nm(jw) * (-4.071805_dp + wave_nm(jw) * (4.301146E-02_dp + &
          wave_nm(jw) * (-1.777846E-04_dp + wave_nm(jw) * 2.520672E-07_dp))) + &
          (T_ref - 300._dp) * EXP(1.234014E+02_dp + wave_nm(jw) * (-2.116255_dp + &
          wave_nm(jw) * (1.111572E-02_dp - wave_nm(jw) * 1.881058E-05_dp))))
      ENDIF
    ENDDO

    ! 64) CHBr3
    l_tdep(ip_chbr3) = .TRUE.
    DO jw = 1,NWAV
      IF (wave_nm(jw)<= 290._dp) THEN
        cst_xxx(jw,ip_chbr3) = cs_xxx(jw,ip_chbr3)
      ELSE
        T_clip=T_ref
        IF (T_clip < 210._dp) T_clip = 210._dp
        IF (T_clip > 300._dp) T_clip = 300._dp
        ! p. 4G-24 of JPL2011:
        cst_xxx(jw,ip_chbr3) = EXP( (0.06183_dp - 0.000241_dp * wave_nm(jw)) * &
          (273._dp - T_clip) - (2.376_dp + 0.14757_dp * wave_nm(jw)))
      ENDIF
    ENDDO
    
    ! op_dc_20161027+
    ! 65) CHF2Cl (JPL2011)
    l_tdep(ip_CHF2Cl) = .TRUE.
    T_clip = T_ref
    IF (T_clip < 210._dp) T_clip = 210._dp
    IF (T_clip > 300._dp) T_clip = 300._dp
    DO jw = 3,18 ! 174 nm <= lambda <= 204 nm
      suma = 0._dp
      sumb = 0._dp
      DO i = 0,3
        suma = suma + A_CHF2Cl(i) * wave_nm(jw)**i
        sumb = sumb + B_CHF2Cl(i) * wave_nm(jw)**i
      ENDDO
      cst_xxx(jw,ip_CHF2Cl) = 10._dp**(suma+(T_clip-273._dp)*sumb)
    ENDDO
    DO jw = 19,26 ! lambda >= 206 nm
      cst_xxx(jw,ip_CHF2Cl) = cs_xxx(jw,ip_CHF2Cl)
    ENDDO

    ! 66) CH3CFCl2 (JPL2011)
    l_tdep(ip_CH3CFCl2) = .TRUE.
    T_clip = T_ref
    IF (T_clip < 210._dp) T_clip = 210._dp
    IF (T_clip > 300._dp) T_clip = 300._dp
    DO jw = 2,36 ! 172 nm <= lambda <= 240 nm
      suma = 0._dp
      sumb = 0._dp
      DO i = 0,4
        suma = suma + A_CH3CFCl2(i) * wave_nm(jw)**i
        sumb = sumb + B_CH3CFCl2(i) * wave_nm(jw)**i
      ENDDO
      cst_xxx(jw,ip_CH3CFCl2) = 10._dp**(suma+(T_clip-273._dp)*sumb)
    ENDDO

    ! 67) CF2ClCF2Cl (JPL2011)
    l_tdep(ip_CF2ClCF2Cl) = .TRUE.
    T_clip = T_ref
    IF (T_clip < 210._dp) T_clip = 210._dp
    IF (T_clip > 300._dp) T_clip = 300._dp
    DO jw = 1,25 ! 172 nm <= lambda <= 220 nm
      suma = 0._dp
      sumb = 0._dp
      DO i = 0,4
        suma = suma + A_CF2ClCF2Cl(i) * wave_nm(jw)**i
        sumb = sumb + B_CF2ClCF2Cl(i) * wave_nm(jw)**i
      ENDDO
      cst_xxx(jw,ip_CF2ClCF2Cl) = 10._dp**(suma+(T_clip-273._dp)*sumb)
    ENDDO
    DO jw = 26,33 ! lambda >= 222 nm
      cst_xxx(jw,ip_CF2ClCF2Cl) = cs_xxx(jw,ip_CF2ClCF2Cl)
    ENDDO

    ! 68) CF2ClCFCl2 (JPL2011)
    l_tdep(ip_CF2ClCFCl2) = .TRUE.    
    T_clip = T_ref
    IF (T_clip < 210._dp) T_clip = 210._dp
    IF (T_clip > 300._dp) T_clip = 300._dp
    DO jw = 3,26 ! 182 nm <= lambda <= 230 nm
      suma = 0._dp
      sumb = 0._dp
      DO i = 0,4
        suma = suma + A_CF2ClCFCl2(i) * wave_nm(jw)**i
        sumb = sumb + B_CF2ClCFCl2(i) * wave_nm(jw)**i
      ENDDO
      cst_xxx(jw,ip_CF2ClCFCl2) = 10._dp**(suma+(T_clip-273._dp)*sumb)
    ENDDO
    DO jw = 27,36 ! lambda >= 232 nm
      cst_xxx(jw,ip_CF2ClCFCl2) = cs_xxx(jw,ip_CF2ClCFCl2)
    ENDDO
    
    ! 69) CF3CF2Cl (JPL2011;CFC-115) no t-dependence
    !l_tdep(ip_CF3CF2Cl) = .TRUE.    
    !T_clip = T_ref
    !IF (T_clip < 210._dp) T_clip = 210._dp
    !IF (T_clip > 300._dp) T_clip = 300._dp
    !DO jw = 1,17 ! 172 nm <= lambda <= 204 nm
    !   suma = 0._dp
    !  sumb = 0._dp
    !  DO i = 0,3
    !    suma = suma + A_CF3CF2Cl(i) * wave_nm(jw)**i
    !  ENDDO
    !  cst_xxx(jw,ip_CF3CF2Cl) = 10._dp**(suma)
    !ENDDO
    !DO jw = 18,23 ! lambda >= 205 nm
    !  cst_xxx(jw,ip_CF3CF2Cl) = cs_xxx(jw,ip_CF3CF2Cl)
    !ENDDO
 
    !69) CH2Cl2 (JPL2011)
    l_tdep(ip_CH2Cl2) = .TRUE.    
    T_clip = T_ref
    IF (T_clip < 210._dp) T_clip = 210._dp
    IF (T_clip > 300._dp) T_clip = 300._dp
    DO jw = 1,23 ! 176 nm <= lambda <= 220 nm
      suma = 0._dp
      sumb = 0._dp
      DO i = 0,4
        suma = suma + A_CH2Cl2(i) * wave_nm(jw)**i
        sumb = sumb + B_CH2Cl2(i) * wave_nm(jw)**i
      ENDDO
      cst_xxx(jw,ip_CH2Cl2) = 10._dp**(suma+(T_clip-273._dp)*sumb)
    ENDDO
    DO jw = 24,41 ! lambda >= 222 nm
      cst_xxx(jw,ip_CH2Cl2) = EXP (-2.1337_dp - (0.08439_dp * wave_nm(jw)))
    ENDDO

    !70) CHCl3 (JPL2011)
    l_tdep(ip_CHCl3) = .TRUE.    
    T_clip = T_ref
    IF (T_clip < 210._dp) T_clip = 210._dp
    IF (T_clip > 300._dp) T_clip = 300._dp
    DO jw = 6,31 ! 190 nm <= lambda <= 240 nm
      suma = 0._dp
      sumb = 0._dp
      DO i = 0,4
        suma = suma + A_CHCl3(i) * wave_nm(jw)**i
        sumb = sumb + B_CHCl3(i) * wave_nm(jw)**i
      ENDDO
      cst_xxx(jw,ip_CHCl3) = 10._dp**(suma+(T_clip-273._dp)*sumb)
    ENDDO
    DO jw = 32,39 ! lambda >= 222 nm
      cst_xxx(jw,ip_CHCl3) = EXP (-1.2277_dp - (0.0844_dp * wave_nm(jw)))
    ENDDO
   ! op_dc_20161027-

  END SUBROUTINE calc_tdep

  ! --------------------------------------------------------------------------

  SUBROUTINE calc_phi

    ! calculation of quantum yields

    IMPLICIT NONE

    REAL(DP), PARAMETER :: A_O1D(20) = (/ &
      0.804_dp,   0.728_dp,    0.865_dp, 0.768_dp, 1.31_dp,   &
      2.29_dp,    5.48_dp,     10.33_dp, 16.21_dp, 25.07_dp,  &
      27.65_dp,   28.43_dp,    24.54_dp, 29.96_dp, 43.91_dp,  &
      37.14_dp,   46.54_dp,    9.18_dp,  0.067_dp, 0.045_dp /)
    REAL(DP), PARAMETER :: B_O1D(20) = (/ &
      9.84_dp,      -16.37_dp,    53.10_dp,    75.57_dp,  305.51_dp,  &
      588.24_dp,    906.80_dp,  1167.66_dp,  1349.77_dp, 1497.28_dp,  &
      1527.01_dp,   1560.50_dp, 1552.14_dp,  1616.80_dp, 1835.12_dp,  &
      1836.45_dp,   2092.44_dp, 1654.54_dp,   238.72_dp,  436.85_dp /)
    REAL(DP) :: T_clip, q1, q2
    INTEGER :: jw
    CHARACTER(LEN=80) :: filesuffix

    WRITE(filesuffix, '(A,I3.3,A,I1.1,A)') '_', NINT(T_ref), 'K_', ji, '.phi'
    IF (l_tdep_loop) filesuffix = '_Tdep'//TRIM(filesuffix)

    ! 2,3) O3
    SELECT CASE(inputdir)
    CASE ("dat_m17")
      T_clip=T_ref
      IF(T_clip<200._dp) T_clip=200._dp
      IF(T_clip>320._dp) T_clip=320._dp
      DO jw=1,38
        phi_xxx(ip_o1d,jw)=0.87_dp
      ENDDO
      DO jw=39,52
        phi_xxx(ip_o1d,jw)=1.98_dp - 301._dp/(wave_nm(jw))
      ENDDO
      DO jw=53,72
        phi_xxx(ip_o1d,jw)=0.06_dp + A_O1D(jw-52)*EXP(-B_O1D(jw-52)/T_clip)
      ENDDO
      DO jw=73,85
        phi_xxx(ip_o1d,jw)=0.06_dp
      ENDDO
      DO jw=86,NWAV
        phi_xxx(ip_o1d,jw)=0._dp
      ENDDO
    CASE ("dat_lit") ! JPL2011:
      DO jw=1,8     ! < 193 nm
        phi_xxx(ip_o1d,jw) = 0.48_dp ! same value as for 193 nm assumed
      ENDDO
      DO jw=9,23    ! 193-225 nm
        ! linear decrease from 0.9 at 225 nm to 0.48 at 193 nm:
        phi_xxx(ip_o1d,jw) = 1.3125E-2_dp * wave_nm(jw) - 2.053125_dp
      ENDDO
      DO jw=24,52   ! 225-305 nm
        phi_xxx(ip_o1d,jw) = 0.9_dp
      ENDDO
      q1 = EXP(     -0._dp/(0.695_dp*T_ref)) ! why not say directly that q1=1 ?
      q2 = EXP(-825.518_dp/(0.695_dp*T_ref))
      DO jw=53,71   ! 306-328 nm
        phi_xxx(ip_o1d,jw) = &
          (q1/(q1+q2)) * 0.8036_dp * &
          EXP(-((304.225_dp-wave_nm(jw))/5.576_dp)**4) + &
          (q2/(q1+q2)) * 8.9061_dp * (T_ref/300._dp)**2 * &
          EXP(-((314.957_dp-wave_nm(jw))/6.601_dp)**2) + &
          0.1192_dp * (T_ref/300._dp)**1.5_dp * &
          EXP(-((310.737_dp-wave_nm(jw))/2.187_dp)**2) + &
          0.0765_dp
      ENDDO
      DO jw=72,74   ! 329-340 nm
        phi_xxx(ip_o1d,jw) = 0.08_dp
      ENDDO
      DO jw=75,NWAV ! > 340 nm
        phi_xxx(ip_o1d,jw) = 0._dp
      ENDDO
    CASE DEFAULT ! not ok:
      STOP "ERROR: inputdir must be dat_m17 or dat_lit"
    END SELECT
    phi_xxx(ip_o3p,:) = 1._dp - phi_xxx(ip_o1d,:)
    ! activate to save quantum yields to a file:
    ! OPEN(IO_TMP, FILE='O1D'//filesuffix, STATUS='UNKNOWN')
    ! DO jw = 1,NWAV
    !   WRITE(IO_TMP, '(F8.2,1PE16.8)') wave_nm(jw), phi_xxx(ip_o1d,jw)
    ! ENDDO
    ! CLOSE(IO_TMP)

    ! 23,67) ClNO3 (ClNO3 and ClONO2)
    SELECT CASE(inputdir)
    CASE ("dat_m17")
      ! not calculated
    CASE ("dat_lit")
      ! JPL2011, p. 4F-24:
      DO jw=1,54 ! < 308 nm
        phi_xxx(ip_clno3,jw) = 0.6_dp
      ENDDO
      DO jw=55,78 ! 308-364 nm
        phi_xxx(ip_clno3,jw) = 7.143E-3_dp * wave_nm(jw) - 1.6_dp ! = lambda/140-1.6
      ENDDO
      DO jw=79,NWAV ! > 364 nm
        phi_xxx(ip_clno3,jw) = 1.0_dp
      ENDDO
      phi_xxx(ip_clono2,:) = 1._dp - phi_xxx(ip_clno3,:)
    CASE DEFAULT ! not ok:
      STOP "ERROR: inputdir must be dat_m17 or dat_lit"
    END SELECT

  END SUBROUTINE calc_phi

  ! --------------------------------------------------------------------------

  SUBROUTINE sr_o2_km(v2s, sro2)

    ! G. A. A. Koppers and D. P. Murtagh
    ! Model studies of the influence of O2 photodissociation
    ! parameterizations in the Schumann-Runge bands on ozone related
    ! photolysis in the upper atmosphere
    ! Ann. Geophys., 14, 68-79, 1996
    ! www.ann-geophys.net/14/68/1996/

    IMPLICIT NONE

    REAL(DP), INTENT(IN)  :: v2s ! slant O2 column
    REAL(DP), INTENT(OUT) :: sro2(13) ! effective O2 cross section

    REAL(DP) :: a, b, dl_o2
    INTEGER :: i

    ! calculation of coefficients A and B and effective O2 cross section:
    DO i=1,13
      IF (v2s>=EXP(38._dp)) THEN ! approx 3e16
        dl_o2 =MIN(56._dp,LOG(v2s))
        a = chebyshev(38._dp,56._dp,cheb_coeff_a(:,i),dl_o2)
        b = chebyshev(38._dp,56._dp,cheb_coeff_b(:,i),dl_o2)
        sro2(i) = EXP(a*(T_ref-220._dp) + b)
      ELSE
        sro2(i) = 0._dp
      ENDIF
    ENDDO

  CONTAINS

    REAL(DP) FUNCTION chebyshev(a,b,c,x)
      ! Chebyshev polynomial evaluation: C(1:M) is an array of Chebyshev
      ! coefficients. The Chebyshev ploynomial sum_{k=1}^{M}
      ! C_k*T_{k-1}(y)-c_1/2 is evaluated at a point
      ! y=(x-(a+b)/2)/(a-b)/2.
      ! See: Numerical recipes in fortran chapter 5.8 P.184 ff
      REAL(DP), INTENT(IN) :: a, b, c(:), x
      INTEGER :: i
      REAL(DP) :: d, dd, sv, y, y2
      IF ((x-a)*(x-b)>0._dp) THEN
        PRINT *, 'ERROR: x not in range in Chebyshev'
        STOP
      ENDIF
      d  = 0._dp
      dd = 0._dp
      y  = (2._dp*x-a-b)/(b-a) ! change of variable to range [-1,1]
      y2 = 2._dp*y
      DO i = SIZE(c),2,-1 ! Clenshaw recurrence
        sv = d
        d  = y2*d-dd+c(i)
        dd = sv
      ENDDO
      chebyshev = y*d-dd + 0.5_dp*c(1) ! last step is different
    END FUNCTION chebyshev

  END SUBROUTINE sr_o2_km

  ! --------------------------------------------------------------------------

  SUBROUTINE sig_eff

    IMPLICIT NONE

    REAL(DP) :: factor, factor2, tau, constant, dlv2, dlv2max, dlv2min, &
      v2max_DU, v2min_DU, v2_m

    INTEGER :: jp, k2, k3, jw, maxv2_or_1

    REAL(DP), ALLOCATABLE, DIMENSION(:) :: &
      V3,    & ! O3 column [part./cm^2]
      V3_DU, & ! O3 column [Dob. units]
      V2       ! O2 column [part./cm^2]
    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: rj_xxx
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: h_no2, h_o2, h_o3
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: SRO2, CRO2
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: &
      TAU_O3, & ! optical depth O3
      TAU_O2, & ! optical depth O2
      FINT      ! integrated actinic flux

    ! allocate arrays:
    ALLOCATE(V3(maxv3), V3_DU(maxv3), V2(maxv2))
    ALLOCATE(rj_xxx(maxv3,maxv2,IP_MAX))
    ALLOCATE(h_no2(maxv3,maxv2))
    ALLOCATE(h_o2(maxv3,maxv2))
    ALLOCATE(h_o3(maxv3,maxv2))
    ALLOCATE(SRO2(13,maxv2), CRO2(maxv2,NWAV))
    ALLOCATE(TAU_O3(maxv3,NWAV), TAU_O2(maxv2,NWAV), FINT(maxv2,maxv3))

    ! change the units ( Dobson Units - part./cm^2) for O3
    ! Boltzmann constant K = 1.38E-23 J/K
    ! normal conditions T0 = 273 K, P0 = 1000 mbar
    ! CONSTANT=K*T0/P0=3.767e-20 cm^3
    constant=3.767e-20_dp

    ! ozone column
    v2max_DU = 0._dp
    v2min_DU = 0._dp
    IF (ji==1) THEN
      v3max_DU = 300._dp
      v3min_DU = 0.00_dp
    ENDIF
    IF (ji==2) THEN
      v3max_DU  = 300._dp
      v3min_DU  = 0.5_dp
      v2max_DU  = 500._dp
      v2min_DU  = 10._dp
    ENDIF
    IF (ji==3) THEN
      v3max_DU = 300._dp
      v3min_DU = 0.5_dp
    ENDIF
    IF (ji>=4) THEN
      v3max_DU = 7000._dp
      v3min_DU = 0.5_dp
    ENDIF

    DO k3=1,maxv3
      v3_du(k3) = v3min_DU + (v3max_DU-v3min_DU)/REAL(maxv3-1,dp)*REAL(k3-1,dp)
    ENDDO

    ! special treatment for Tdep calculations with only 4 values:
    IF (maxv3==4) THEN
      v3_du(1) =   0.01_dp
      v3_du(2) =  70._dp
      v3_du(3) = 120._dp
      v3_du(4) = 180._dp
      v3min_DU = MINVAL(v3_du)
      v3max_DU = MAXVAL(v3_du)
    ENDIF

    DO k3=1,maxv3
      v3(k3)=v3_du(k3)*1._dp/constant*1.E-3_dp
    ENDDO

    ! O2 column:
    IF (ji==1) THEN
      dlv2max = 55.5_dp
      dlv2min = 44.0_dp
      DO k2=1,maxv2 ! ln(V2(K))-ln(V2(K-1)) = 0.05_dp
        dlv2   = dlv2min+(dlv2max-dlv2min)/REAL(maxv2-1,dp)*REAL(k2-1,dp)
        v2(k2) = EXP(dlv2)
      ENDDO
    ELSE
      DO k2=1,maxv2
        v2_m   = v2min_DU + (v2max_DU-v2min_DU)/REAL(maxv2-1,dp)*REAL(k2-1,dp)
        v2(k2) = v2_m * 1._dp/constant*1.E2_dp
      ENDDO
    ENDIF

    ! special treatment for Tdep calculations with only 2 values:
    IF (maxv2==2) THEN
      v2(1) =  20._dp * 1._dp/constant*1.E2_dp
      v2(2) = 200._dp * 1._dp/constant*1.E2_dp
    ENDIF

    ! O2 optical depth:
    IF (ji==1) THEN
      ! Schumann-Runge bands, Koppers and Murtagh parameterization
      DO k2 = 1,maxv2
        CALL sr_o2_km(V2(k2),SRO2(:,k2))
      ENDDO
    ENDIF
    DO k2 = 1,maxv2
      cro2(k2,:) = cst_xxx(:,ip_o2)
      DO jw = 1,13
        ! SRO2 is only defined if ji=1
        IF (ji==1) THEN
          cro2(k2,jw)=sro2(jw,k2)
        ELSE
          cro2(k2,jw)=0._dp ! dummy value
        ENDIF
      ENDDO
    ENDDO

    IF (lf<=13) THEN
      ! optical depths (approximation formula from Koppers and Murtagh)
      DO jw=li,lf
        tau_o2(1,jw)  = sro2(jw,1)*v2(1)
        DO k2=2,maxv2
          tau_o2(k2,jw) = tau_o2(k2-1,jw) - &
            (sro2(jw,k2-1)*v2(k2-1) - sro2(jw,k2)*v2(k2)) / &
            LOG(sro2(jw,k2-1)*v2(k2-1)/(sro2(jw,k2)*v2(k2))) * &
            LOG(v2(k2-1)/v2(k2))
        ENDDO
      ENDDO
    ELSE
      DO jw=li,lf
        DO k2=1,maxv2
          tau_o2(k2,jw) = v2(k2) * cst_xxx(jw,ip_o2)
        ENDDO
      ENDDO
    ENDIF

    ! o3 optical depth
    DO jw=li,lf
      DO k3=1,maxv3
        tau_o3(k3,jw) = v3(k3) * cst_xxx(jw,ip_o1d)
      ENDDO
    ENDDO

    rj_xxx(:,:,:) = 0._dp

    h_no2(:,:) = 0._dp
    h_o2(:,:)  = 0._dp
    h_o3(:,:)  = 0._dp

    DO k3=1,maxv3
      DO k2=1,maxv2
        fint(k2,k3)=0._dp
        DO jw=li,lf
          tau = tau_o3(k3,jw) + tau_o2(k2,jw) ! total optical depth
          factor = flx(jw) * EXP(-tau)
          DO jp = 1,ip_max
            IF ((jp==ip_o2).AND.(jw<=13)) cst_xxx(jw,jp) = sro2(jw,k2)
            rj_xxx(k3,k2,jp) = &
              rj_xxx(k3,k2,jp) + cst_xxx(jw,jp) * factor * phi_xxx(jp,jw)
          ENDDO
          factor2 = 2.E9_dp/(7._dp*BOLTZ) * PLANCK * C_LIGHT / wave_nm(jw) * factor
          h_o2(k3,k2)  = h_o2(k3,k2)  + factor2 * cro2(k2,jw)
          h_o3(k3,k2)  = h_o3(k3,k2)  + factor2 * cst_xxx(jw,ip_o1d) * 1.e-6_dp
          h_no2(k3,k2) = h_no2(k3,k2) + factor2 * cst_xxx(jw,ip_no2) * 1.e-9_dp
          fint(k2,k3)   = fint(k2,k3) + factor
        ENDDO
      ENDDO
    ENDDO

    flx_sum = SUM(flx(li:lf)) ! integrated extraterrestrial flux

    IF(ji<=2) THEN
      v2min = V2(1)
      v2max = V2(maxv2)
      maxv2_or_1 = maxv2
    ELSE
      maxv2_or_1 = 1
    ENDIF

    IF(ji/=1) THEN
      OPEN(IO_OUT,file=DIR_EFF//'/tau.ef'//TRIM(suffix),status='unknown')
      WRITE(IO_OUT,*)'# tau'
      DO k2=1,maxv2_or_1
        WRITE(IO_OUT,'(1P,6E12.4)')(-LOG(fint(k2,k3)/flx_sum), k3=1,maxv3)
      ENDDO
      CLOSE(IO_OUT)
    ENDIF

    OPEN(IO_OUT,file=DIR_EFF//'/v3_du.ef'//TRIM(suffix),status='unknown')
    WRITE(IO_OUT,'(A)') '# v3_du'
    DO k3=1,maxv3
      WRITE(IO_OUT,'(1P,E12.4)') v3_du(k3)
    ENDDO
    CLOSE(IO_OUT)

    IF (ji<=2) THEN
      OPEN(IO_OUT,file=DIR_EFF//'/v2.ef'//TRIM(suffix),status='unknown')
      WRITE(IO_OUT,'(A)') '# v2'
      DO k2=1,maxv2
        WRITE(IO_OUT,'(1P,E12.4)') v2(k2)
      ENDDO
      CLOSE(IO_OUT)
    ENDIF

    ! loop over all photolysis reactions:
    DO jp = 1, IP_MAX
      IF ((l_species(jp).AND..NOT.l_tdep_loop) .OR. &
        (l_tdep(jp).AND.l_tdep_loop)) THEN
        CALL write_xxx(TRIM(jname(jp)), rj_xxx(:,:,jp), fint, maxv2_or_1, maxv3)
      ENDIF
    ENDDO
    CALL write_xxx('h_O2',  h_o2,  fint, maxv2_or_1, maxv3)
    CALL write_xxx('h_O3',  h_o3,  fint, maxv2_or_1, maxv3)
    CALL write_xxx('h_NO2', h_no2, fint, maxv2_or_1, maxv3)

    DEALLOCATE(v3, v3_du, v2, sro2, cro2, tau_o3, tau_o2, fint, rj_xxx)
    DEALLOCATE(h_o2, h_o3, h_no2)

  END SUBROUTINE sig_eff

  !---------------------------------------------------------------------------

  SUBROUTINE write_xxx(pspecies, rj, fint, maxv2_or_1, maxv3)

    IMPLICIT NONE
    CHARACTER(LEN=*),     INTENT(IN) :: pspecies
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: rj, fint
    INTEGER,              INTENT(IN) :: maxv2_or_1, maxv3

    INTEGER :: k2, k3

    OPEN(IO_OUT,FILE=DIR_EFF//'/'//pspecies//'.ef'// &
      TRIM(suffix),status='UNKNOWN')
    WRITE(IO_OUT,'(A)') '# '//pspecies
    DO k2=1,maxv2_or_1
      DO k3=1,maxv3
        WRITE(IO_OUT,'(1P,E12.4)') rj(k3,k2)/fint(k2,k3)
      ENDDO
    ENDDO
    CLOSE(IO_OUT)

  END SUBROUTINE write_xxx

  !---------------------------------------------------------------------------

  SUBROUTINE print_info_header(title)

    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: title

    WRITE(*,'(A)')
    WRITE(*,'(A)') HLINE2
    WRITE(*,'(A)') title
    WRITE(*,'(A)') HLINE2
    WRITE(*,'(A)') '# T_ref   wave_i wave_f flx_sum'// &
      '  v2min   v2max  v3min v3max '
    WRITE(*,'(A)') HLINE2
    WRITE(*,'(A)') '                        photons'// &
      '   mcl     mcl               '
    WRITE(*,'(A)') '     K      nm     nm   -------'// &
      '   ---     ---    DU    DU   '
    WRITE(*,'(A)') '                         cm2*s '// &
      '   cm2     cm2               '
    WRITE(*,'(A)') HLINE2

  END SUBROUTINE print_info_header

  !---------------------------------------------------------------------------

  SUBROUTINE print_info(T_range)

    IMPLICIT NONE
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: T_range

    WRITE(*,'(I1)', ADVANCE='NO') ji
    IF (PRESENT(T_range)) THEN
      WRITE(*,'(A8)', ADVANCE='NO') T_range
    ELSE
      WRITE(*,'(F8.2)', ADVANCE='NO') T_ref
    ENDIF
    WRITE(*,'(2F7.2,1P,E8.1)', ADVANCE='NO') &
      wave_nm(li), wave_nm(lf), flx_sum
    IF(ji<=2)THEN
      WRITE(*,'(1P,2E8.1)', ADVANCE='NO') v2min, v2max
    ELSE
      WRITE(*,'(A16)', ADVANCE='NO') ''
    ENDIF
    WRITE(*,'(F5.2,F7.0)') v3min_DU, v3max_DU

  END SUBROUTINE print_info

  !---------------------------------------------------------------------------

  SUBROUTINE calc_interval

    IMPLICIT NONE
    INTEGER, PARAMETER :: &
      lini(NBIN) = (/  1, 14, 30, 44, 53, 61, 74,  91 /), & ! initial
      lfin(NBIN) = (/ 13, 29, 43, 52, 60, 73, 90, 142 /)    ! final

    ! note that out of the 176 wave lengths, only 1...142 are used
    li = lini(ji)
    lf = lfin(ji)

    CALL calc_tdep  ! temperature-dependent cross sections
    CALL calc_phi   ! quantum yields
    CALL sig_eff    ! effective values

  END SUBROUTINE calc_interval

  !---------------------------------------------------------------------------

  SUBROUTINE conv_176_eff

    USE jvpp_mem, ONLY: HLINE1
    IMPLICIT NONE
    INTEGER, PARAMETER :: NTEMP = 15 ! number of temperatures
    INTEGER :: jtemp, jp, jw
    CHARACTER(LEN=80) :: tdepfile

    ! ------------------------------------------------------------------------

    WRITE(*,'(A)') HLINE1
    WRITE(*,'(A)') '*** STEP 2'
    WRITE(*,'(A)') HLINE1

    CALL initialize

    ! ------------------------------------------------------------------------

    CALL print_info_header('Calculations at a fixed reference temperature:')
    l_tdep_loop = .FALSE.
    DO ji=1,NBIN
      WRITE(suffix,'(I1.1)') ji
      IF (ji==1) THEN ! Schumann-Runge bands:
        T_ref = 240._dp
        maxv3 = 127
        maxv2 = 231
      ELSE ! interval 2-8:
        T_ref = 250._dp
        maxv3 = 600
        maxv2 =  50
      ENDIF
      CALL calc_interval
      CALL print_info()
    ENDDO
    WRITE(*,'(A)') HLINE2

    ! ------------------------------------------------------------------------

    CALL print_info_header('Calculations at several temperatures:')
    l_tdep_loop = .TRUE.
    DO ji=1,NBIN
      DO jtemp=1,NTEMP
        WRITE(suffix,'(I1.1,A,I3.3,A)') ji, '_', 10*jtemp+170, 'K'
        T_ref = REAL(10*jtemp+170,dp) ! = 180, 190, 200, ... 320 K
        IF (ji==1) THEN ! Schumann-Runge bands:
          maxv3 =   4
          maxv2 = 231
        ELSE ! interval 2-8:
          maxv3 = 600
          maxv2 =   2
        ENDIF
        CALL calc_interval
      ENDDO
      CALL print_info('180-320')
    ENDDO
    WRITE(*,'(A)') HLINE2

    ! print temperature-dependent data to files for plotting:
    DO jtemp=1,NTEMP
      T_ref = REAL(10*jtemp+170,dp) ! = 180, 190, 200, ... 320 K
      CALL calc_tdep  ! temperature-dependent cross sections
      ! loop over all photolysis reactions:
      DO jp=1,IP_MAX
        WRITE(tdepfile,'(4A,I3.3,A)') &
          DIR_JNL, '/', TRIM(jname(jp)), '_', 10*jtemp+170, 'K.tdp'
        OPEN(IO_TDP, FILE=TRIM(tdepfile), STATUS='UNKNOWN')
        DO jw = 1,NWAV
          WRITE(IO_TDP,*) wave_nm(jw), cst_xxx(jw,jp)
        ENDDO
        CLOSE(IO_TDP)
      ENDDO
    ENDDO

    ! ------------------------------------------------------------------------

    WRITE(*,*)
    DO jp = 1, IP_MAX ! loop over all photolysis reactions
      WRITE(*,'(A)', ADVANCE='NO') jname(jp)//' '
      IF (l_s2t(jp).OR.l_s3t(jp)) THEN
        WRITE(*,'(A)', ADVANCE='NO')   '(Tdep from file) '
        WRITE(*,'(A)') TRIM(filename_tdep(jp))
      ELSEIF (l_sig(jp)) THEN
        IF (l_tdep(jp)) THEN
          WRITE(*,'(A)', ADVANCE='NO') '(other Tdep)     '
        ELSE
          WRITE(*,'(A)', ADVANCE='NO') '(no Tdep)        '
        ENDIF
        WRITE(*,'(A)') TRIM(filename(jp))
      ELSEIF (l_species(jp)) THEN
        IF (l_tdep(jp)) THEN
          WRITE(*,'(A)', ADVANCE='NO') '(other Tdep)     '
        ELSE
          WRITE(*,'(A)', ADVANCE='NO') '(no Tdep)        '
        ENDIF
        WRITE(*,'(A)') 'calculated spectrum'
      ELSE
        WRITE(*,'(A)') '---              ---'
      ENDIF
    ENDDO

    PRINT *

  END SUBROUTINE conv_176_eff

END MODULE jvpp_step2

! ****************************************************************************
