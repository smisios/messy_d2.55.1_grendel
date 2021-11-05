MODULE messy_main_tools_kinetics

  USE MESSY_MAIN_CONSTANTS_MEM,     ONLY: DP

  IMPLICIT NONE

  PUBLIC
  PRIVATE  :: DP

CONTAINS

  !=============================================================
  ! FUNCTIONS FOR KPP
  !=============================================================

  ELEMENTAL REAL(DP) FUNCTION K_SIV_H2O2 (K_298,TDEP,CHP,TEMP)
    ! SPECIAL RATE FUNCTION FOR S(IV) + H2O2

    REAL,     INTENT(IN) :: K_298 ! K AT T = 298.15K
    REAL,     INTENT(IN) :: TDEP  ! TEMPERATURE DEPENDENCE
    REAL(DP), INTENT(IN) :: CHP   ! C(H+)
    REAL(DP), INTENT(IN) :: TEMP  ! TEMPERATURE

    INTRINSIC :: REAL

    K_SIV_H2O2 = REAL(K_298,DP) * EXP(REAL(TDEP,DP)*(1._DP/TEMP-3.3540E-3_DP)) &
      * CHP / (CHP+0.1_DP)

  END FUNCTION K_SIV_H2O2

  !=================================================================

  ELEMENTAL REAL(dp) FUNCTION k_limited (k3rd,cHp)
    ! diffusion limitation caps 3rd order rate coefficients
    REAL(dp), INTENT(IN) :: k3rd  ! 3rd order rate coefficient
    REAL(dp), INTENT(IN) :: cHp   ! c(H+)
    REAL(dp), PARAMETER  :: DiffLimit = 1E10 ! diffusion limitation [M-1s-1]
    INTRINSIC :: EXP
    k_limited = 1._dp / ( 1._dp/k3rd + cHp/DiffLimit)
  END FUNCTION k_limited

  !=================================================================

  ELEMENTAL REAL(dp) FUNCTION k_3rd_iupac(temp,cair,k0_300K,n,kinf_300K,m,fc)
    ! IUPAC three body reaction formula (iupac.pole-ether.fr)
    INTRINSIC :: LOG10
    REAL(dp), INTENT(IN) :: temp      ! temperature [K]
    REAL(dp), INTENT(IN) :: cair      ! air concentration [molecules/cm3]
    REAL,     INTENT(IN) :: k0_300K   ! low pressure limit at 300 K
    REAL,     INTENT(IN) :: n         ! exponent for low pressure limit
    REAL,     INTENT(IN) :: kinf_300K ! high pressure limit at 300 K
    REAL,     INTENT(IN) :: m         ! exponent for high pressure limit
    REAL,     INTENT(IN) :: fc        ! broadening factor (e.g. 0.45 or 0.6...)
    REAL                 :: nu        ! N
    REAL                 :: zt_help, k0_T, kinf_T, k_ratio

    zt_help = 300._dp/temp
    k0_T    = k0_300K   * zt_help**(n) * cair ! k_0   at current T
    kinf_T  = kinf_300K * zt_help**(m)        ! k_inf at current T
    k_ratio = k0_T/kinf_T
    nu      = 0.75-1.27*LOG10(fc)
    k_3rd_iupac = k0_T/(1._dp+k_ratio)* &
      fc**(1._dp/(1._dp+(LOG10(k_ratio)/nu)**2))
  END FUNCTION k_3rd_iupac

  !=================================================================

  ELEMENTAL REAL(dp) FUNCTION alpha_AN(n,ro2type,bcarb,gcarb,abic,temp,cair)
  ! Alkyl nitrate yields dependent on T and P according to
  ! Arey. ref3202 and Teng, ref3189
    INTRINSIC :: LOG10
    INTEGER,  INTENT(IN) :: n         ! number of heavy atoms (C, O, N) except the O atom of beta-carbonyls
    INTEGER,  INTENT(IN) :: ro2type   ! 1, 2 or 3 for primary, secondary and tertiary RO2
    INTEGER,  INTENT(IN) :: bcarb     ! 1 for beta-carbonyl group, 0 for none
    INTEGER,  INTENT(IN) :: gcarb     ! 1 for gamma-carbonyl group, 0 for none
    INTEGER,  INTENT(IN) :: abic      ! 1 for bicyclic peroxy from aromatics, 0 for none
    REAL(dp), INTENT(IN) :: temp      ! temperature [K]
    REAL(dp), INTENT(IN) :: cair      ! air concentration [molecules/cm3]
    REAL(dp), PARAMETER  :: alpha=2.E-22, beta=1.0, Yinf_298K=0.43, F=0.41, m0=0., minf=8.0
    REAL(dp)             :: m           ! factor for primary, secondary and tertiary RO2
    REAL(dp)             :: bcf,gcf,abf ! reduction factors for the presence of beta-carbonyl and gamma-carbonyl group and a bicyclic peroxy from aromatics
    REAL(dp)             :: Y0_298K, Y0_298K_tp, Yinf_298K_t, zeta, k_ratio

    m = 1. ! According to Teng, ref3189
    IF (bcarb .EQ. 1) THEN
      bcf = 0.19 ! derived from Praske, ref3190: alpha_AN = 0.03 for the secondary HMKO2 relative to alpha_AN for 6C RO2 (0.16)
    ELSE IF (bcarb .EQ. 0) THEN
      bcf = 1.
    ELSE
      bcf = 1.
    ENDIF
    IF (gcarb .EQ. 1) THEN
      gcf = 0.44 ! derived from Praske, ref3190: alpha_AN = 0.07 for the primary HMKO2 relative to alpha_AN for 6C RO2 (0.16)
    ELSE IF (gcarb .EQ. 0) THEN
      gcf = 1.
    ELSE
      gcf = 1.
    ENDIF
    IF (abic .EQ. 1) THEN
      abf = 0.24 ! ratio of AN-yield for toluene (Elrod et al., ref3180), 5.5% at 200 torr, and this SAR for linear alkyl RO2 with 9 heavy atoms, 23.3%
    ELSE IF (abic .EQ. 0) THEN
      abf = 1.
    ENDIF
    Y0_298K     = alpha*EXP(beta*n)
    Y0_298K_tp  = Y0_298K * cair * (temp/298)**(-m0)
    Yinf_298K_t = Yinf_298K * (temp/298)**(-minf)
    zeta        = 1/(1+LOG10(Y0_298K_tp/Yinf_298K_t)**2)
    k_ratio     = (Y0_298K_tp/(1+Y0_298K_tp/Yinf_298K_t))*F**zeta
    alpha_AN    = k_ratio/(1+k_ratio) * m * bcf * gcf * abf
  END FUNCTION alpha_AN

  !=================================================================

  ELEMENTAL REAL(dp) FUNCTION k_N2_O(temp,temp_ion)
    REAL(dp), INTENT(IN) :: temp      ! temperature [K]
    REAL(dp), INTENT(IN) :: temp_ion  ! ion temperature [K]
    REAL                 :: temp_mean
    temp_mean = (temp_ion + temp)/2
    k_N2_O = 1.4E-10*(300./temp_mean)**0.44
  END FUNCTION k_N2_O

  !=================================================================

  ELEMENTAL REAL(dp) FUNCTION k_Op_O2(temp,temp_ion)
    REAL(dp), INTENT(IN) :: temp      ! temperature [K]
    REAL(dp), INTENT(IN) :: temp_ion  ! ion temperature [K]
    REAL(dp)             :: temp_mean
    temp_mean  = 0.667*temp_ion + 0.333*temp
    k_Op_O2  = 2.82E-11 - 7.74E-12*(temp_mean/300.) + &
      1.073E-12*(temp_mean/300.)**2 - 5.17E-14*(temp_mean/300.)**3 + &
      9.65E-16*(temp_mean/300.)**4
  END FUNCTION k_Op_O2

  !=================================================================

  ELEMENTAL REAL(dp) FUNCTION k_Op_N2(temp,temp_ion)
    REAL(dp), INTENT(IN) :: temp      ! temperature [K]
    REAL(dp), INTENT(IN) :: temp_ion  ! ion temperature [K]
    REAL(dp)             :: temp_mean
    temp_mean = 0.6363*temp_ion + 0.3637*temp
    k_Op_N2 = 1.533E-12 - 5.92E-13*(temp_mean/300.) + &
         8.6E-14*(temp_mean/300.)**2
  END FUNCTION k_Op_N2

  !=================================================================

  ELEMENTAL REAL(dp) FUNCTION k_RO2_HO2(temp,nC)
    real(dp), intent(in) :: temp      ! temperature [K]
    integer, intent(in)  :: nC        ! q-ty of C atoms
    k_RO2_HO2 = 2.91E-13*EXP(1300./temp)*(1.-EXP(-0.245*REAL(nC,dp))) ! ref1630
  END FUNCTION k_RO2_HO2

  !=================================================================

  ! This function has been introduced for sensitivity studies on OClO
  ! within the PRACE project CASiMIR (Chemistry of the Atmosphere Simulated
  ! with an Earth System Model for the Interpretation of Satellite based
  ! Remote sensing observations)
  ELEMENTAL REAL(dp) FUNCTION uef(temp)

    ! uncertainty estimate function
    ! (approximated by 3rd order polynomial) for G7603

    REAL(dp), INTENT(IN) :: temp      ! temperature [K]

    REAL(dp), PARAMETER :: a0 =  8.4096792e+000_dp
    REAL(dp), PARAMETER :: a1 = -6.8484593e-002_dp
    REAL(dp), PARAMETER :: a2 =  2.3044184e-004_dp
    REAL(dp), PARAMETER :: a3 = -2.7257885e-007_dp

    uef = (a0 + a1 * temp + a2 * temp**2 + a3 * temp**3)

  END FUNCTION uef

  !=================================================================

END MODULE messy_main_tools_kinetics
