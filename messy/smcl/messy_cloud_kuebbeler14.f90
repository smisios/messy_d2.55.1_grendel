MODULE MESSY_CLOUD_KUEBBELER14

  ! ---------------------------------------------------------------------------
  ! Author
  !   Mattia Righi (DLR-IPA, Oberpfaffenhofen)
  !   based on messy_cloud_lohmann10 structure by Holger Tost
  !
  ! Original code
  !   Ulrike Lohmann (ETH, Zurich)
  !   Sylvaine Ferrachat (ETH, Zurich)
  !   Miriam Kuebbeler (ETH, Zurich)
  !
  ! Source
  !   URL: https://svn.iac.ethz.ch/external/echam-hammoz/echam5-ham/echam-
  !        5.5.00_rc2/branches/ethz/branches/miriam
  !   Repository Root: https://svn.iac.ethz.ch/external/echam-hammoz
  !   Repository UUID: 5abd231c-38af-11dd-ab7a-a2415d6bb57b
  !   Revision: 1928
  ! ---------------------------------------------------------------------------


  USE messy_main_constants_mem, ONLY: dp, api => pi, g, ak => k_B, m_air, &
                                      M_H2O, STRLEN_MEDIUM
  USE messy_cloud_ori,          ONLY: crhoi, cqtmin
  USE messy_cloud_mem

  IMPLICIT NONE
  PRIVATE
  SAVE

  ! --- Default settings for namelist parameters (&CTRL_K14) ------------------
  ! ---------------------------------------------------------------------------

  ! Number of freezing aerosol particle types in cirrus (including hom. freez.)
  ! (only for cloud_param = 5 and icnc = 2)
  INTEGER, PUBLIC  :: nfrzaer = 4

  ! Number of external ice modes
  ! 0 = without pre-existing ice crystals
  ! 1 = with pre-existing ice crystals
  INTEGER, PUBLIC  :: nexti = 1

  ! INP properties
  INTEGER, PUBLIC :: ninp = 0
  INTEGER, PARAMETER, PUBLIC :: max_ninp = 8
  TYPE nml_set_inp
     ! INP name
     CHARACTER(LEN=STRLEN_MEDIUM) :: name = ''
     ! Critical supersaturation
     REAL(dp)                     :: Scrit = 1.
     ! Fraction of active INP
     REAL(dp)                     :: f_active = -1.
  END TYPE nml_set_inp
  TYPE(nml_set_inp), DIMENSION(max_ninp), PUBLIC :: inp_properties

  ! Minimum cloud droplet number concentration [m^(-3)]
  REAL(dp), PUBLIC :: cdncmin = 20.e6_dp

  ! Vertical velocity scheme (default: large-scale + TKE)
  INTEGER, PUBLIC :: vervel_scheme = 1

  ! Scaling factor for large-scale vertical velocity in cirrus
  REAL(dp), PUBLIC:: scale_v_ls_cirrus = 1._dp

  ! Scaling factor for tke vertical velocity in cirrus
  REAL(dp), PUBLIC :: scale_v_tke_cirrus = 0.7_dp

  ! Scaling factor for orographic vertical velocity in cirrus
  REAL(dp), PUBLIC:: scale_v_orogw_cirrus = 1._dp

  ! Scaling factor for orographic vertical velocity in cirrus

  ! Transform IC larger than 100 micron to snow
  LOGICAL, PUBLIC :: l_ic2snow  = .TRUE.

  ! --- Other parameters  -----------------------------------------------------
  ! ---------------------------------------------------------------------------

  ! Auto-conversion rate
  ! 1 = Beheng (Atmospheric Res., 1994) - ECHAM5 standard
  ! 2 = Khairoutdinov and Kogan (Mon. Weather Rev., 2000)
  INTEGER, PUBLIC :: nauto = 2

  ! Consider convective microphyisical scheme (added by Junhua Zhang)
  LOGICAL :: ncvmicro = .FALSE.

  ! Logical for thermophoresis in freezing of cloud water
  LOGICAL :: lthermo = .FALSE.


  ! --- Physical constants and limits -----------------------------------------
  ! ---------------------------------------------------------------------------

  ! Epsilons
  REAL(dp), PARAMETER :: zepsec = 1.0e-12_dp
  REAL(dp), PARAMETER :: zxsec = 1._dp - zepsec
  REAL(dp), PARAMETER :: zqsec = 1._dp - cqtmin
  REAL(dp), PARAMETER :: zeps = EPSILON(1.0_dp)

  ! Volume mean radius of ice crystals (to estimate CDNC from melting IC)
  REAL(dp), PARAMETER :: zcri = 10.E-6_dp  ! [m]

  ! Assumed mass of ice crystals with zcri radius [kg]
  REAL(dp), PARAMETER :: zmi = 4._dp/3._dp*zcri**3*api*crhoi

  ! Minimum and maximum ICNC
  REAL(dp), PARAMETER :: zicemin = 10._dp   ! [m-3]
  REAL(dp), PARAMETER :: zicemax = 1.e7_dp  ! [m-3]

  ! Some constants used in Levkov et al. (1992)
  REAL(dp), PARAMETER :: zsigmaw = 0.28_dp
  REAL(dp), PARAMETER :: zcdi = 0.6_dp
  REAL(dp), PARAMETER :: zmw0 = 4.19e-12_dp
  REAL(dp), PARAMETER :: zmi0 = 1.e-12_dp
  REAL(dp), PARAMETER :: zmi0_rcp = 1.e12_dp
  REAL(dp), PARAMETER :: zka = 0.024_dp

  ! Thermal conductivity of carbon (Seinfeld and Pandis, 1996; p. 481)
  REAL(dp), PARAMETER :: zkbc = 4.2_dp

  ! Thermal conductivity of clay (Seinfeld and Pandis, 1996; p. 481)
  REAL(dp), PARAMETER :: zkdu = 0.72_dp

  ! Boltzmann constant
  REAL(dp), PARAMETER :: zkb = 1.38e-23_dp ! [m2 kg s-2 K-1]

  ! Deposition coefficient alpha
  REAL(dp), PARAMETER :: zalpha = 0.5_dp

  ! Mass of a water molecule
  REAL(dp), PARAMETER :: zxmw = 2.992e-26_dp ! [kg]

  ! Correction factor applied to fall velocity
  REAL(dp), PARAMETER :: zfall = 3._dp

  ! Specific density of ice
  REAL(dp), PARAMETER :: zrhoice = 925._dp  ! [kg/m3]


  ! --- Public subroutines ----------------------------------------------------
  ! ---------------------------------------------------------------------------
  PUBLIC :: cloud_cdnc_icnc3
  PUBLIC :: cloud_read_nml_ctrl_k14


CONTAINS


  ! ===========================================================================


  ! --- SUBROUTINE cloud_cdnc_icnc3 -------------------------------------------
  ! ---------------------------------------------------------------------------
  !
  ! Description
  !   Computes large-scale water phase changes, precipitation, cloud cover, and
  !   vertical integrals of specific humidity, cloud liquid water content and
  !   cloud ice (diagnostics).
  !
  !   This routine computes the tendencies of the four prognostic variables
  !   (temperature t, specific humidity q, cloud liquid water xl, cloud ice xi)
  !   due to phase changes (condensation/deposition, evaporation/sublimation of
  !   rain/snow falling into the unsaturated part of the grid box, melting of
  !   snow, melting/freezing of cloud ice/cloud water, sedimentation of cloud
  !   ice, and precipitation formation in warm, cold and mixed phase clouds.
  !
  !   The precipitation at the surface (rain and snow) is used in later for
  !   computing the land surface hydrology in *surf*. The cloud parameters
  !   (cloud cover, cloud liquid water and cloud ice are used for the
  !   calculation of radiation at the next timestep.
  !
  !   In the current version the advective tendencies of skewness and variance
  !   are set to zero.
  !
  ! References
  !   Lohmann and Roeckner, 1996: Clim. Dyn. 557-572
  !   Levkov et al., 1992: Beitr. Phys. Atm. 35-58 (ice phase)
  !   Beheng, 1994: Atmos. Res. 193-206 (warm phase)
  !   Lenderink et al., 1998; KNMI-REPORT NO. 98-13 (condensation)
  !   Tompkins 2002, J. Atmos. Sci. (cloud cover)
  !
  ! Authors
  !   M. Esch (MPI, Hamburg, 1999)
  !   G. Lenderink (KNMI, de Bilt, 1998)
  !   U. Lohmann (MPI, Hamburg, 1995)
  !
  ! Modifications
  !   E. Roeckner (MPI, Hamburg,  2000)
  !   A. Tompkins (MPI, Hamburg,  2000)
  !   U. Schlese (MPI, Hamburg,  2003)
  !   U. Lohmann (Dalhousie University, 2002-2006)
  !   P. Stier (MPI, Hamburg, 2002-2006)
  !   J. Zhang (Dalhousie University, 2004)
  !   S. Ferrachat (ETH, Zurich, 2008)
  !
  ! Structure
  !
  !  0 - INITIALISATION
  !
  !  1 - TOP BOUNDARY CONDITIONS, AIR DENSITY
  !      1.1 - Set precipitation fluxes to zero
  !      1.2 - Air density
  !      1.3 - Supersaturations
  !      1.4 - Liquid phase
  !      1.5 - Ice phase
  !      1.6 - Cirrus clouds parametrization
  !
  !  2 - SET TO ZERO LOCAL TENDENCIES
  !      2.1 - Transfer cirrus IC with R >= 100 micron to snow
  !      2.2 - Transformation of multimodal ICNC and Rice to unimodal
  !
  !  3 - MODIFICATION OF INCOMING PRECIPITATION FLUXES
  !      3.1 - Melting of snow and ice
  !      3.2 - Submimation of snow (zsub) and ice (zxisub)
  !      3.3 - Evaporation of rain
  !
  !  4 - SEDIMENTATION OF CLOUD ICE FROM GRID-MEAN VALUES
  !
  !  5 - CONDENSATION/DEPOSITION AND EVAPORATION/SUBLIMATION
  !      5.1 - Turbulence: Skewness
  !      5.2 - Turbulence: Variance
  !      5.3 - Deposition/sublimation of cloud ice and condensation/evaporation
  !            of liquid water due to changes in water vapour and temperature
  !      5.4 - Account for cloud evaporation in clear air and check for
  !            supersaturation
  !      5.5 - Change of in-cloud water due to deposition/sublimation and
  !            condensation/evaporation
  !
  !  6 - FREEZING OF CLOUD WATER
  !      6.1 - Freezing of all cloud water for T < 238 K
  !      6.2 - Freezing of cloud water between 238 and 273 K
  !
  !  7 - CLOUD PHYSICS AND PRECIPITATION FLUXES AT THE SURFACE
  !      7.1 - Warm clouds
  !      7.2 - Cold clouds
  !      7.3 - Update precipitation fluxes
  !
  !  8 - UPDATE TENDENCIES OF T, Q, XL, XI AND FINAL CLOUD COVER
  !      8.1 - Cloud cover scheme tendencies
  !      8.2 - Tendencies of thermodynamic variables
  !      8.4 - Corrections: avoid negative cloud water and cloud ice
  !      8.4 - Cirrus-specific diagnostics to compare with in-situ data
  !
  !  9 - WET CHEMISTRY AND IN-CLOUD SCAVENGING
  !
  ! 10 - DIAGNOSTICS
  !      10.1 - Accumulated precipitation at the surface
  !      10.2 - Total cloud cover
  !      10.3 - Vertical integrals of humidity, cloud water and cloud ice

  ! ---------------------------------------------------------------------------
  SUBROUTINE cloud_cdnc_icnc3( &

       ! GENERAL INDICES
       kproma, kbdim,       & ! kproma and nproma indices
       ktdia, klev, klevp1, & ! level indices (ktdia = 1)
       ztmst,               & ! time-step [s]
       ktrac,               & ! number of cloud tracers (= 2, CDNC and ICNC)
       krow,                & ! jrow index
       !
       ! VARIABLES FOR OROGRAPHIC WAVES
       sqcst_2d,  & ! cos(lat) [1]
       w_gwd_kpr, & ! vertical vel. induced by orogr. waves [m/s]
       itest,     & ! gridpoints where orogr. waves are present
       l_z,       & ! half-wavelength seen by the incident flow [m]
       ampl_gwd,  & ! wave amplitude [m]
       !
       ! THERMODYNAMIC VARIABLES
       paphm1,      & ! pressure at half levels (t-1) [Pa]
       pvervel,     & ! large-scale vertical velocity [Pa/s]
       pvervel_p18, & ! sub-scale vertical velocity (Penner Ansatz) [cm/s]
       papm1,       & ! pressure at full levels (t-1) [Pa]
       papp1,       & ! pressure at full levels (t+1) [Pa]
       pacdnc,      & ! standard ECHAM5 cloud droplet number concentration [m-3]
       pqm1,        & ! specific humidity [kg/kg]
       ptm1,        & ! air temperature [K]
       ptvm1,       & ! virtual temperature [K]
       pxlm1,       & ! cloud water content [kg/kg] - grid-box value
       pxim1,       & ! coud ice content [kg/kg] - grid-box value
       pxtec,       & ! detrained convective cloud water or ice [kg/kg]
       pxvar,       & ! variance of total water amount qv+qi+ql [kg/kg]
       pxskew,      & ! skewness of total water amount qv+qi+ql [1]
       pqtec,       & ! convective detrained humidity [kg/kg]
       pbetaa,      & ! beta distribution minimum a
       pbetab,      & ! beta distribution maximum b
       pvdiffp,     & ! rate of change of q due to vdiff scheme [1/s]
       phmixtau,    & ! inverse mixing timescale for horiz. turbulence [1/s]
       pvmixtau,    & ! inverse mixing timescale for vertical turbulence [1/s]
       pgeo,        & ! geopotential [m2/s2]
       pbetass,     & !
       pxtm1,       & ! tracer for CDNC (ktrac=1) and ICNC (ktrac=2) [1/mol]
       ptkem1,      & ! turbulent kinetic energy [m2/s2]
       pcvcbot,     & ! bottom level of convection / convective cloud base [1]
       pxtecl,      & ! convective detrained water [kg/kg]
       pxteci,      & ! convective detrained ice [kg/kg]
       pxtecnl,     & ! convective detrained water numbers [1/kg]
       pxtecni,     & ! convective detrained ice numbers [1/kg]
       knvb,        & !
       !
       ! CLOUD VARIABLES (OUTPUT)
       paclc,   & ! large-scale cloud cover [1]
       paclcac, & ! large-scale cloud cover accumulated [1]
       prelhum, & ! relative humidity (0-1), diagnostic [1]
       paclcov, & ! total cloud cover [1]
       paprl,   & ! accumulated large-scale surface-level precipiatation rate
                  ! (rain+snow) [kg/m2/s]
       pqvi,    & ! vertically integrated water vapor [kg/m2]
       pxlvi,   & ! vertically integrated cloud water (LWP) [kg/m2]
       pxivi,   & ! vertically integrated cloud ice (IWP) [kg/m2]
       pssfl,   & ! large-scale surface-level snow rate [kg/m2/s]
       prsfl,   & ! large-scale surface-level rain rate [kg/m2/s]
       !
       ! TENDENCIES OF THERMODYNAMIC AND CLOUD VARIABLES
       pqte,          & ! tendency of specific humidity [kg/kg/s]
       ptte,          & ! tendency of temperature [K/s]
       pxlte,         & ! tendency of cloud water [kg/kg/s]
       pxite,         & ! tendency of cloud ice [kg/kg/s]
       pxtte,         & ! tendency of CDNC and ICNC (ktrac = 1 or 2) [1/mol/s]
       !
       ! INPUT FROM THE BASE MODEL
       paprs,         & ! accum. large-scale surface-level snow rate [kg/m2/s]
       status_string, & ! status string for look-up error
       cirrus_status, & ! status from the cirrus parametrization
       cirrus_string, & ! error message from the cirrus parametrization
       lcover,        & ! namelist switch for extra cloud cover calculation
       slm,           & ! land-sea mask [1]
       glac,          & ! fraction of land covered by glaciers [1]
       pcdncact,      & ! CDNC from the activation scheme [m-3]
       !
       ! CHANNEL OUTPUT
       plwc,          & ! large-scale cloud water content [kg/kg]
       piwc,          & ! large-scale cloud ice content [kg/kg]
       pfrain,        & ! large-scale rain flux [kg/m2/s]
       pfsnow,        & ! large-scale snow flux [kg/m2/s]
       pfrain_no,     & ! large-scale rain flux without cloud production of
                        ! new rain [kg/m2/s]
       pfsnow_no,     & ! large-scale snow flux without cloud production of
                        ! new snow [kg/m2/s]
       prate_r,       & ! large-scale rain formation inside cloud [kg/kg]
       prate_s,       & ! large-scale snow formation inside cloud [kg/kg]
       prevap,        & ! large-scale rain evaporation [kg/kg]
       pssubl,        & ! large-scale snow sublimation [kg/kg]
       pr_cover,      & ! large-scale precipitation cloud cover [1]
       pcond,         & ! condensate in cloud-covered part of gridbox [1]
       pimelt,        & ! large-scale frozen precipitation melting [kg/m2/s]
       pisedi         & ! large-scale ice sedimentation [kg/kg]
      )


  ! --- Use statements --------------------------------------------------------
  ! ---------------------------------------------------------------------------
  USE messy_cloud_ori,  ONLY : cqtmin, tmelt, cvtfall, crhosno, cn0s,    &
                               cthomi, csecfrl, ncctop, cvarmin,         &
                               cbeta_pq, cbeta_pq_max, nbetaq, cbetaqs,  &
                               rbetak, nbetax, tbetai0, tbetai1,         &
                               clmax, clmin, jbmin, jbmax, lonacc,       &
                               ccraut, crhoi, ccsaut,                    &
                               cbeta_cs, lookupoverflow

  USE messy_main_tools, ONLY : jptlucu1, jptlucu2, tlucua, tlucuaw, tlucub

  USE messy_main_constants_mem, ONLY: ceffmin, ceffmax, ccwmin, rd, rv,   &
                                      vtmpc1, vtmpc2, cpd => cp_air,      &
                                      tmelt, rhoh2o => rho_h2o

  ! -- Declarations -----------------------------------------------------------
  ! ---------------------------------------------------------------------------
  IMPLICIT NONE

  CHARACTER(LEN=32) :: status_string
  CHARACTER(LEN=70) :: cirrus_string

  LOGICAL :: lcover

  INTEGER :: krow, ktdia, kproma, kbdim, klev, klevp1, ktrac, jrow, jl, &
       jk, jkk, iqidx, ixidx, ntype, ntt, tbin, ibin, cirrus_status

  LOGICAL, DIMENSION(kbdim) :: lo_1d, lo2_1d, locc_1d, ll1_1d, ll2_1d, &
       ll3_1d, ll4_1d, ll5_1d, ll6_1d, ll7_1d, ll8_1d

  LOGICAL, DIMENSION(kbdim, klev) :: ll_look, ll_cv, ll_ice, ll1, ll1_2d, &
       ll2_2d, ll3_2d

  INTEGER, DIMENSION(kbdim) :: knvb, it1_1d, it_1d

  INTEGER, DIMENSION(kbdim, klev) :: itop, ibas, icl_minusbas, icl_minustop, &
       iclbas, itm1_look, itm1p1_look

  REAL(dp), PARAMETER :: zdummy = 1._dp

  REAL(dp) :: zesw, zf1, zrieff, zradl, zdisp, zdw0, ztmst, ztmst_rcp,        &
       zcons1, zcons2, zcons3, zcons4, zcons5, zdtime, zdtime_rcp, ztmstdt,   &
       zdt, g_rcp, zpirho, zpirho_rcp, zb2, zdv, zgtp, zexm1_1, zexp_1,       &
       zexm1_2, zexp_2, zrih, ztte, zomega, ztc, zvth, zfuchs, zfre, zre,     &
       zfracdusol, zfracduinsolai, zfracduinsolci, zfracbcsol, zfracbcinsol,  &
       zfrzcntdu, zfrzcntbc, zfrzcnt, znaimmdu, znaimmbc, zfrzimm,            &
       zfrzthermo, zqsm1, zqst1, zdqsdtime,         &
       zdeltatemp, zf2, zcap, zdiag

  REAL(dp), DIMENSION(kbdim) :: sqcst_2d, itest, l_z, pxlvi, pxivi, paclcov,  &
       paprl, pqvi, pssfl, paprs, prsfl, slm, glac, zclcpre, zcnd, zdep,      &
       zxievap, zxlevap, zimlt, zsmlt, zspr, zxlte, zxite, zxiflux, zsacl,    &
       zxlte2, zxite2, zlsdcp, zlvdcp, zximlt, ztp1tmp, zqp1tmp, zxisub,      &
       zxlb, zxib, zxibold_1d, zclcov, zclcaux, zqvi, zxlvi, zxivi, zbetaqt,  &
       zwide, zbetacl, zturbvar, zturbskew, zconvvar, zconvskew, zvartg,      &
       zmicroskew, zgenti, zgentl, zxvarte, zxskewte, pcvcbot, zrprn, &
       zsacln, zfrln, zcdnc_burden, zreffct, ztau1i, zicnc_burden, zscnc_1d,  &
       zdplanar_1d, zusnow_1d, zsacl2in_1d, zxsp2_1d, zxsp_1d, zstcrit_1d,    &
       zstokes_1d, zudrop_1d, zrey_1d, zrac1_1d, zrac2_1d, zrautself_1d,      &
       zxsp1_1d, zraut_1d, zsaut_1d, zsaci2_1d, zsacl2_1d, zris_1d,           &
       zcolleffi_1d, zc1_1d, zxlp1_1d, zreffl_1d, zdxlcor_1d, zdxicor_1d,     &
       zsecprod_1d, zcsacl_1d, zpretot_1d, zpredel_1d, zpresum_1d, zzdrr_1d,  &
       zzdrs_1d, zxidtstar_1d, zlc_1d, zxldtstar_1d, zxlm1evp_1d,             &
       zxim1evp_1d, zxldt_1d, zxidt_1d, zqsm1_1d, zqp1_1d, zdtdt_1d,          &
       zdqsat_1d, ztp1_1d, zdqsdt_1d, zqst1_1d, zqvdt_1d, zxip1_1d, zcorw_1d, &
       zesw_1d, zoversatw_1d, zqsp1tmpw_1d, zqsp1tmp_1d, zcor_1d, zrhtest_1d, &
       zoversat_1d, zlcdqsdt_1d, zclcstar_1d, zxrp1_1d, zauloc_1d,            &
       zqsp1tmpcirrus_1d, zqcon_1d, zrelhum_1d, iqidx_1d,                     &
       zbap1_1d, zgent_1d, ixidx_1d, zqtau_1d, zxilb_1d, zbbap1_1d, zbqp1_1d, &
       zqcdif_1d, zlucuawp1_1d, zlucuap1_1d, zes_1d, zlucub_1d, zlucuaw_1d,   &
       zlucua_1d, zicncp1_1d, zetaair_1d,  &
       zkair_1d, zdfarbcki_1d, zdfarduai_1d, zdfarduci_1d, zknbcki_1d,        &
       zknduai_1d, zknduci_1d, zftbcki_1d, zftduai_1d, zftduci_1d,            &
       zvervmax_1d, zrice_1d, zeta_1d, zdv_1d, ztmp1_1d, ztmp2_1d, ztmp3_1d,  &
       ztmp4_1d, ztmp5_1d

  REAL(dp), DIMENSION(kbdim, klev) :: w_gwd_kpr, ampl_gwd, pvervel,           &
       pvervel_p18, papm1, pqm1, papp1, ptm1,ptvm1, pxlm1, pxim1, pxtec,      &
       pqtec, pxvar, pxskew, pbetaa, pbetab, pvdiffp, phmixtau, pvmixtau,     &
       pgeo, pbetass, paclc, paclcac, pacdnc, prelhum, ptte, pqte, pxlte,     &
       pxite, zfrl, zrho, zrho_rcp, zmratepr, zmrateps, zfrain,   &
       zfsnow, zfevapr, zfsubls, zmlwc, zmiwc, zmsnowacl, pcdncact, ptkem1,   &
       pxtecl, pxteci, pxtecnl, pxtecni, zsprn, zcdnc, zqlnuc, zqlnuccvh,     &
       zqlnuccv, zqlnuc_bas, zcdnc_bas, zicnc, zicncq,                        &
       zascs, zrid_2d, zqsi_2d, zninucl, zqinucl, zri, znidetr, zdz_2d,       &
       zdp_2d, zdpg_2d, zaaa_2d, zviscos_2d, zqswp1_2d, zqsip1_2d, zqsw_2d,   &
       zastbstw, zastbsti, zmmean_2d, zxifallmc_2d, zalfased_2d, zbetased_2d, &
       zxifallnc_2d, zxised_2d, zicesed_2d, zesw_2d, zesi_2d, zsusatw_2d,     &
       zsusatw_evap, zsusatix_2d, zsusatix_2d_lim, zvervx_2d, zicesub,        &
       zqrho_2d, ztmp, ztmp1, ztmp2, wnucl

  REAL(dp), DIMENSION(kbdim,klevp1) :: paphm1, zgeoh

  REAL(dp), DIMENSION(kbdim,klev,ktrac) :: pxtm1, pxtte

  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev)         :: plwc,     piwc
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev), TARGET :: pfrain,   pfsnow
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev)         :: pfrain_no,pfsnow_no
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev), TARGET :: prevap,   pssubl
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev), TARGET :: prate_r,  prate_s
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev)         :: pr_cover
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev)         :: pimelt,   pisedi
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev)         :: pcond

  REAL(dp), POINTER, DIMENSION(:) :: zrfl, zsfl, zsub, zevp, zrpr

  ! --- Some variable definitions ---------------------------------------------
  ! ---------------------------------------------------------------------------
  ! zmratepr      rain formation rate in cloudy part of the grid box [kg/kg]
  ! zmrateps      ice  formation rate in cloudy part of the grid box [kg/kg]
  ! zfrain        rain flux before evaporation [kg/m2/s]
  ! zfsnow        snow flux before sublimation [kg/m2/s]
  ! zfevapr       evaporation of rain [kg/m2/s]
  ! zfsubls       sublimation of snow [kg/m2/s]
  ! zmlwc         in-cloud liquid water mmr before rain formation [kg/kg]
  ! zmiwc         in-cloud ice mmr before snow formation [kg/kg]
  ! zmsnowacl     accretion rate of snow with cloud droplets in cloudy
  !               part of the box [kg/kg]
  ! zcdnc         cloud droplet number concentration (in-cloud) [m-3]
  ! zicnc         ice crystal number concentration (in-cloud)  [m-3]
  ! zqlnuccvh     nucleated CDNC from convective detrainment
  ! zw_strat      stratiform updraft velocity, large-scale+TKE (>0.0) [m s-1]
  ! zw_conv       convective updraft velocity, large-scale+CAPE (>0.0) [m s-1]
  ! zninucl       number conc. of newly nucleated IC
  ! zqinucl       mass mixing ratio of newly nucleated IC
  ! zri           size of newly nucleated IC
  ! zreffct       cloud top effective radius
  ! ztau1i        ice cloud optical depth - visible wavelength
  ! znidetr       iC number from detrainment
  ! zsusatw_2d    supersaturation with respect to water
  ! zvervx_2d     updraft [cm/s]
  ! zcd2ic        amount of cdnc erroneously present at T < T_hom that should
  !               be credited to zicnc
  ! zic2cd        amount if icnc erroneously present at T > T_melt that should
  !               be credited to zcdnc
  ! zlsdcp        latent heat of sublimation divided by the specific heat at
  !               constant pressure
  ! zlvdcp        latent heat of vaporization divided by the specific heat at
  !               constant pressure
  !
  ! --- Variables for the cirrus scheme and the coupling to aerosol -----------
  ! ---------------------------------------------------------------------------
  ! ntype         total number of ice particle types
  ! znicex        ice crystal number concentration from the cirrus scheme [m-3]
  ! znicex_3d     ice crystal number concentration for each freezing mode [m-3]
  ! zri_3d        ice crystal radius for each freezing mode [m]
  ! znicex_tmp    temporary array for ICNC
  ! zri_tmp       temporary array for IC radius
  ! zri_crt       mean volume radius of IC of last tstep and detrained ones [m]
  ! zicncq1       fraction of zicncq weighted with the sum of IWC
  ! znidetr1      fraction of znidetr weighted with the sum of IWC
  ! zsicirrus_2d  highest freezing threshold reached within a nucleation event
  ! zapnx         aerosol number conc. available for each freezing mode [cm-3]
  ! zaprx         aerosol radius available for each freezing mode [cm]
  ! zapsx         stdandard deviation of aerosol size distribution available
  !               for each freezing mode [1]
  ! f_du_active   active fraction of dust particles [1]
  ! scrhet        critical saturation ratio for heterogeneous freezing [1]
  ! zap           sum of zapnx over the freezing modes [cm-3]
  ! zcaer         total aerosol number concentration per type [cm-3]
  ! zraer         radius of the smallest particle per type [cm]
  ! zsaut2_1d     variable for IC conversion to snow if R>100 micron
  ! zlwc_strat    liquid water content in the stratiform part [kg/kg]
  ! zlwc_detr     liquid water content in the detrained part [kg/kg]
  ! zlwc_tot      liquid water content in the strat. + detrained parts [kg/kg]
  ! zqinucl_3d    mass mixing ratio of IC from depositional growth [kg/kg]
  ! zvice         fictitious downdraft accounting for depositional growth of
  !               preexisting ice crystals

  REAL(dp), DIMENSION(kbdim) :: zsaut2_1d, zcd2ic, zic2cd, zcd2unphys, &
       zxifluxn

  REAL(dp), DIMENSION(kbdim, klev) :: znicex, znicex_tmp, zri_tmp, zri_crt,   &
       zicncq1, znidetr1, zsicirrus_2d, f_du_active, zap, zxtec, zlwc_strat,  &
       zlwc_detr, zlwc_tot, zqinucl_tmp, zicncq_vice, zri_vice, zvipr, zvit,  &
       zvit1, zvirhoa, zvipi, zvicis, zvia1, zvia2, zvia3, zvisi, zvivsc,     &
       zviflx, zvib1, zvib2, zvidl, zvifall, zvirey, zvivent, zvice,          &
       zrwetki_2d, zrwetai_2d, zrwetci_2d

  REAL(dp), DIMENSION(kbdim, klev, nfrzaer) :: zcaer, zraer

  REAL(dp), DIMENSION(kbdim, klev, nfrzaer+nexti) :: znicex_3d, zri_3d

  REAL(dp), DIMENSION(kbdim, nfrzaer, nfrzmod) :: zapnx, zaprx, zapsx

  REAL(dp), DIMENSION(kbdim, klev, nfrzaer+1) :: zqinucl_3d

  REAL(dp), DIMENSION(2, nfrzaer-1) :: scrhet

  ! For INP properties
  CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(ninp) :: inp_name
  REAL(dp), DIMENSION(ninp) :: inp_Scrit, inp_f_active

  ! --- Executable statements -------------------------------------------------
  ! ---------------------------------------------------------------------------

  jrow = krow

  ! === 0 - INITIALISATION ====================================================
  ! ===========================================================================

  lookupoverflow = .FALSE.

  ! Initialize channgel objects
  plwc(:,:)      = 0.0_dp
  piwc(:,:)      = 0.0_dp
  pfrain(:,:)    = 0.0_dp
  pfsnow(:,:)    = 0.0_dp
  pfrain_no(:,:) = 0.0_dp
  pfsnow_no(:,:) = 0.0_dp
  prate_r(:,:)   = 0.0_dp
  prate_s(:,:)   = 0.0_dp
  prevap(:,:)    = 0.0_dp
  pssubl(:,:)    = 0.0_dp
  pr_cover(:,:)  = 0.0_dp
  pimelt(:,:)    = 0.0_dp
  pisedi(:,:)    = 0.0_dp

  ICNC_burden_acc(:,jrow) = 0._dp
  icnc_burden(:,jrow)     = 0._dp
  cdnc_burden(:,jrow)     = 0._dp
  cdnc_burden_acc(:,jrow) = 0._dp
  icnc_acc(:,:,jrow)      = 0._dp
  cdnc_acc(:,:,jrow)      = 0._dp
  cdnc_insitu(:,:,jrow)   = 0._dp
  reffi(:,:,jrow)         = ceffmax ! = 150 um
  reffl(:,:,jrow)         = 40._dp  ! [um]
  cdnc(:,:,jrow)          = 0._dp
  icnc(:,:,jrow)          = 0._dp

  ! Diagnostics
  vervel_ls(:,:,jrow)          = 0._dp
  vervel_tke(:,:,jrow)         = 0._dp
  vervel_gw(:,:,jrow)          = 0._dp
  CIRRUS_IWC(:,:,:,jrow)       = 0._dp
  CIRRUS_Nice(:,:,:,jrow)      = 0._dp
  CIRRUS_Nice_ML(:,:,:,jrow)   = 0._dp
  CIRRUS_Rice(:,:,:,jrow)      = 0._dp
  CIRRUS_RHi_cloud(:,:,:,jrow) = 0._dp
  CIRRUS_RHi_clear(:,:,:,jrow) = 0._dp
  CIRRUS_vervel(:,:,:,jrow)    = 0._dp
  Nice_preex(:,:,jrow)     = 0._dp
  Nice_DUdep(:,:,jrow)     = 0._dp
  Nice_DUimm(:,:,jrow)     = 0._dp
  Nice_BC(:,:,jrow)        = 0._dp
  Nice_BCtag(:,:,jrow)     = 0._dp
  Nice_homog(:,:,jrow)     = 0._dp

  zmratepr(:,:)  = 0._dp
  zmrateps(:,:)  = 0._dp
  zfrain(:,:)    = 0._dp
  zfsnow(:,:)    = 0._dp
  zfevapr(:,:)   = 0._dp
  zfsubls(:,:)   = 0._dp
  zmlwc(:,:)     = 0._dp
  zmiwc(:,:)     = 0._dp
  zmsnowacl(:,:) = 0._dp

  zqinucl(:,:) = 0._dp
  zninucl(:,:) = 0._dp
  zreffct(:)   = 0._dp
  ztau1i(:)    = 0._dp

  zxifluxn(:)       = 0.0_dp
  zmmean_2d(:,:)    = zmi  ! [kg]
  zxifallmc_2d(:,:) = 0.0_dp
  zxifallnc_2d(:,:) = 0.0_dp
  zxised_2d(:,:)    = 0.0_dp
  zicesed_2d(:,:)   = 0.0_dp

  zcdnc_burden(1:kproma) = 0.0_dp
  zcdnc_bas(1:kproma,:)  = 0.0_dp
  zqlnuc(1:kproma,:)     = 0.0_dp
  zqlnuc_bas(1:kproma,:) = 0.0_dp
  zicnc_burden(1:kproma) = 0.0_dp
  zri(1:kproma,:)        = 1.e-6_dp

  zvipr(kbdim, klev)     = 0.0_dp
  zvit(kbdim, klev)      = 0.0_dp
  zvit1(kbdim, klev)     = 0.0_dp
  zvirhoa(kbdim, klev)   = 0.0_dp
  zvipi(kbdim, klev)     = 0.0_dp
  zvicis(kbdim, klev)    = 0.0_dp
  zvia1(kbdim, klev)     = 0.0_dp
  zvia2(kbdim, klev)     = 0.0_dp
  zvia3(kbdim, klev)     = 0.0_dp
  zvisi(kbdim, klev)     = 0.0_dp
  zvivsc(kbdim, klev)    = 0.0_dp
  zviflx(kbdim, klev)    = 0.0_dp
  zvib1(kbdim, klev)     = 0.0_dp
  zvib2(kbdim, klev)     = 0.0_dp
  zvidl(kbdim, klev)     = 0.0_dp
  zvifall(kbdim, klev)   = 0.0_dp
  zvirey(kbdim, klev)    = 0.0_dp
  zvivent(kbdim, klev)   = 0.0_dp
  zvice(kbdim, klev)     = 0.0_dp

  ! Computational constant
  zdisp = EXP(0.5_dp * zsigmaw**2)
  zdw0 = 10.e-6_dp * zdisp  ! modal diameter (10 um) * dispersion parameter
  zdtime = ztmst / 2._dp
  zdtime_rcp = 1._dp / zdtime
  ztmst_rcp = 1._dp / ztmst
  ztmstdt = ztmst * zdtime
  zdt = ztmst_rcp * zdtime
  g_rcp = 1._dp / g
  zcons1 = cpd * vtmpc2
  zcons2 = ztmst_rcp * g_rcp
  zexm1_1 = 2.47_dp - 1._dp
  zexp_1 = -1._dp / zexm1_1
  zexm1_2 = 4.7_dp - 1._dp
  zexp_2 = -1._dp / zexm1_2
  zpirho = api * rhoh2o
  zpirho_rcp = 1._dp / zpirho
  zcap = 2._dp / api
  zcons3 = 1._dp / (api * crhosno * cn0s *cvtfall**(1._dp/1.16_dp))**0.25_dp
  zcons4 = 1._dp / (api * crhosno * cn0s )**0.8125_dp
  zcons5 = 1._dp / (api * crhosno * cn0s )**0.875_dp

  zrwetki_2d(1:kproma,:) = mode(iaiti)%wetrad(1:kproma,:)
  zrwetai_2d(1:kproma,:) = mode(iacci)%wetrad(1:kproma,:)
  zrwetci_2d(1:kproma,:) = mode(icoai)%wetrad(1:kproma,:)

  ! Cloud utility variables
  CALL get_util_var(kproma, kbdim, ktdia, klev, klevp1, paphm1, pgeo, papm1, &
                    ptm1, zgeoh, zdp_2d, zdpg_2d, zdz_2d, zaaa_2d, zviscos_2d)


  ! === 1 - TOP BOUNDARY CONDITIONS, AIR DENSITY ==============================
  ! ===========================================================================

  ! --- 1.1 - Set precipitation fluxes to zero --------------------------------
  ! ---------------------------------------------------------------------------

  zclcpre(:) = 0._dp ! fraction of grid box covered by precip
  zxiflux(:) = 0._dp ! IC flux falling into the grid box from above

  ! --- 1.2 - Air density -----------------------------------------------------
  ! ---------------------------------------------------------------------------

  DO 122 jk = ktdia, klev  ! vertical loop (bottom to top)

     ! Air density definitions
     zrho(:,jk) = papm1(:,jk) / (rd * ptvm1(:,jk))  ! [kg/m3]
     zrho_rcp(:,jk) = 1._dp / zrho(:,jk)  ! 1/rho [m3/kg]
     zqrho_2d(:,jk) = 1.3_dp * zrho_rcp(:,jk)
     zfrl(:,jk) = 0.0_dp

     ! Added for convective microphysical scheme
     IF (.not.ncvmicro) THEN

        pxtec(:,jk)  = MAX(pxtec(:,jk), 0.0_dp)

        ll1_1d(:) = (pxtec(:,jk) > 0.0_dp)

        ztmp1_1d(:) = twc_conv(1:kproma,jk,jrow) + ztmstdt*pxtec(:,jk)
        twc_conv(1:kproma,jk,jrow) = &
             MERGE(ztmp1_1d(:), twc_conv(1:kproma,jk,jrow), ll1_1d(:))

        ztmp1_1d(:) = conv_time(1:kproma,jk,jrow) + zdtime
        conv_time(1:kproma,jk,jrow) = &
             MERGE(ztmp1_1d(:), conv_time(1:kproma,jk,jrow), ll1_1d(:))

     ELSE

        pxtecl(:,jk)  = MAX(pxtecl(:,jk),  0.0_dp)
        pxteci(:,jk)  = MAX(pxteci(:,jk),  0.0_dp)
        pxtecnl(:,jk) = MAX(pxtecnl(:,jk), 0.0_dp)
        pxtecni(:,jk) = MAX(pxtecni(:,jk), 0.0_dp)

        ll1_1d(:) = (pxtecl(:,jk) > 0._dp) .AND. (ptm1(:,jk)   < cthomi)

        ztmp1_1d(:)  = pxteci(:,jk) + pxtecl(:,jk)
        pxteci(:,jk) = MERGE(ztmp1_1d(:), pxteci(:,jk), ll1_1d(:))
        pxtecl(:,jk) = MERGE(0._dp, pxtecl(:,jk), ll1_1d(:))

        ztmp1_1d(:)   = pxtecni(:,jk) + pxtecnl(:,jk)
        pxtecni(:,jk) = MERGE(ztmp1_1d(:), pxtecni(:,jk), ll1_1d(:))
        pxtecnl(:,jk) = MERGE(0._dp, pxtecnl(:,jk), ll1_1d(:))

        zfrl(:,jk) = MERGE(pxtecl(:,jk), zfrl(:,jk), ll1_1d(:))

        ll1_1d(:) = (pxteci(:,jk) > 0.0_dp) .AND. (ptm1(:,jk)   > tmelt)

        IF (ANY(ll1_1d(1:kproma))) THEN
           ztmp1_1d(:)  = pxteci(:,jk) + pxtecl(:,jk)
           pxtecl(:,jk)  = MERGE(ztmp1_1d(:), pxtecl(:,jk), ll1_1d(:))
           pxteci(:,jk)  = MERGE(0._dp, pxteci(:,jk), ll1_1d(:))

           ztmp1_1d(:)   = pxtecni(:,jk) + pxtecnl(:,jk)
           pxtecnl(:,jk) = MERGE(ztmp1_1d(:), pxtecnl(:,jk), ll1_1d(:))
           pxtecni(:,jk) = MERGE(0.0_dp, pxtecni(:,jk), ll1_1d(:))
        ENDIF

        ll1_1d(:) = (pxtecl(:,jk) > 0.0_dp)
        pxtecnl(:,jk) = MERGE(pxtecnl(:,jk), 0._dp, ll1_1d(:))

        ll1_1d(:) = (pxteci(:,jk) > 0.0_dp)
        pxtecni(:,jk) = MERGE(pxtecni(:,jk), 0._dp, ll1_1d(:))

     ENDIF

     pqtec(:,jk)  = MAX(pqtec(:,jk), 0.0_dp)

     ! Cloud droplet number concentration [1/mol] --> [1/m3]
     ! NOTE: the tendency must be added for a correct start value
     ! (same as in Lohmann07 and in the original ECHAM code)
!!$  zcdnc(:,jk) = zrho(:,jk) * (pxtm1(:,jk,1) / m_air * 1000._dp)
     zcdnc(:,jk) = zrho(:,jk) * ( (pxtm1(:,jk,1) + pxtte(:,jk,1)*ztmst) &
          / m_air * 1000._dp )
     zcdnc(:,jk) = MAX(zcdnc(:,jk), cqtmin)  ! lower limit 1.e-12

     ! Ice crystal number concentration [1/mol] --> [1/m3]
     ! As above, the tendency must be added
!!$  zicnc(:,jk) = zrho(:,jk) * (pxtm1(:,jk,2) / m_air * 1000._dp)
     zicnc(:,jk) = zrho(:,jk) * ( (pxtm1(:,jk,2)  + pxtte(:,jk,2)*ztmst ) &
          / m_air * 1000._dp )
     zicnc(:,jk)  = MAX(zicnc(:,jk), cqtmin)  ! lower limit 1.e-12

     ! Correct for inconsistencies relative to phases and temperature
     ! Refined correction for crediting to ice *only* when cloud cover and
     ! LWC are non-zero (see https://redmine.hammoz.ethz.ch/issues/477)
     ll1_1d(:) = (ptm1(:,jk) <= cthomi) .AND. (zcdnc(:,jk) > cdncmin)
     ll2_1d(:) = (paclc(:,jk) >= zepsec) .AND. &
                 (pxlm1(:,jk) + ztmst * pxlte(:,jk) >= zepsec)

     zcd2ic(:) = MERGE(zcdnc(:,jk), 0._dp, ll1_1d(:) .AND. ll2_1d(:))
     zcd2unphys(:) = MERGE(zcdnc(:,jk), 0._dp, ll1_1d(:) .AND. .not. ll2_1d(:))

     ll1_1d(:) = (ptm1(:,jk) > tmelt) .AND. (zicnc(:,jk) > zicemin)
     zic2cd(:) = MERGE(zicnc(:,jk), 0._dp, ll1_1d(:))

     zcdnc(:,jk) = zcdnc(:,jk) - zcd2ic(:) - zcd2unphys(:) + zic2cd(:)
     zicnc(:,jk) = zicnc(:,jk) + zcd2ic(:)                 - zic2cd(:)

     zcdnc(:,jk) = MAX(zcdnc(:,jk), cqtmin)
     zicnc(:,jk) = MAX(zicnc(:,jk), cqtmin)

     zicncq(:,jk) = zicnc(:,jk)  ! initialize

122 ENDDO  ! vertical loop


  ! --- 1.3 - Supersaturations ------------------------------------------------
  ! ---------------------------------------------------------------------------

  ! Cloud bottom and cloud top
  CALL get_cloud_bounds(kproma, kbdim, ktdia, klev, paclc, itop, ibas, &
                        icl_minustop, icl_minusbas)

  ! Search current temperature in the look-up table
  ztmp(:,:)      = 1000._dp * ptm1(:,:)
  itm1_look(:,:) = NINT(ztmp(:,:))
  ll_look(:,:) = (itm1_look(:,:) < jptlucu1 .OR. itm1_look(:,:) > jptlucu2)
  IF (ANY(ll_look(1:kproma,ktdia:klev))) lookupoverflow = .TRUE.
  itm1_look(:,:) = MAX(MIN(itm1_look(:,:), jptlucu2), jptlucu1)
  itm1p1_look(:,:) = itm1_look(:,:) + 1
  itm1p1_look(:,:) = MAX(MIN(itm1p1_look(:,:), jptlucu2), jptlucu1)

  DO jk = klev, ktdia, -1  ! reverse vertical loop (top to bottom)
     DO jl = 1, kproma

        ! Water
        zesw = tlucuaw(itm1_look(jl,jk)) / papm1(jl,jk)
        zesw = MIN(zesw, 0.5_dp)
        zqsw_2d(jl,jk) = zesw / (1._dp - vtmpc1 * zesw)
        zsusatw_2d(jl,jk) = pqm1(jl,jk) / zqsw_2d(jl,jk) - 1._dp
        zsusatw_2d(jl,jk) = MAX(zsusatw_2d(jl,jk), 0._dp)

        ! For evaporation of rain in 3.3:
        zsusatw_evap(jl,jk) = pqm1(jl,jk) / zqsw_2d(jl,jk) - 1._dp
        zsusatw_evap(jl,jk) = MIN(zsusatw_evap(jl,jk), 0._dp)

        ! Saturation water vapour pressure
        zesw_2d(jl,jk) = zesw * papm1(jl,jk) * rv / rd

        ! For later use in 5
        zqswp1_2d(jl,jk) = tlucuaw(itm1p1_look(jl,jk)) / papm1(jl,jk)
        zqswp1_2d(jl,jk) = MIN(zqswp1_2d(jl,jk), 0.5_dp)
        zqswp1_2d(jl,jk) = &
             zqswp1_2d(jl,jk) / (1._dp-vtmpc1*zqswp1_2d(jl,jk))

        ! ICE
        zqsi_2d(jl,jk) = tlucua(itm1_look(jl,jk)) / papm1(jl,jk)
        zqsi_2d(jl,jk) = MIN(zqsi_2d(jl,jk), 0.5_dp)
        zesi_2d(jl,jk) = zqsi_2d(jl,jk) * papm1(jl,jk) * rv / rd
        zqsi_2d(jl,jk) = zqsi_2d(jl,jk) / (1._dp - vtmpc1 * zqsi_2d(jl,jk))
        sice(jl,jk,jrow) = pqm1(jl,jk) / zqsi_2d(jl,jk) - 1._dp
        sice(jl,jk,jrow) = MAX(sice(jl,jk,jrow), 0._dp)

        ! For sublimation of snow and ice in 3.2
        zicesub(jl,jk) = pqm1(jl,jk) / zqsi_2d(jl,jk) - 1._dp
        zicesub(jl,jk) = MIN(zicesub(jl,jk), 0._dp)

        ! For later use in 5
        zqsip1_2d(jl,jk) = tlucua(itm1p1_look(jl,jk)) / papm1(jl,jk)
        zqsip1_2d(jl,jk) = MIN(zqsip1_2d(jl,jk), 0.5_dp)
        zqsip1_2d(jl,jk) = &
             zqsip1_2d(jl,jk) / (1._dp - vtmpc1 * zqsip1_2d(jl,jk))

     ENDDO
  ENDDO  ! reverse vertical loop

  ! Store supersaturation with respect to water
  swat(1:kproma,:,jrow)  = zsusatw_2d(1:kproma,:)

  ! Some utility variables
  zastbstw(:,:) = & ! water
       alv * (alv / (rv * ptm1(:,:)) - 1.0_dp) / (zka * ptm1(:,:)) + &
       rv * ptm1(:,:) / (2.21_dp / papm1(:,:) * zesw_2d(:,:))

  ztmp1(:,:) = 4.1867e-3_dp * (5.69_dp + 0.017_dp * (ptm1(:,:) - tmelt))

  zastbsti(:,:) = & ! ice
       als * (als / (rv * ptm1(:,:)) - 1.0_dp) / (ztmp1(:,:) * ptm1(:,:)) + &
       rv * ptm1(:,:) / (2.21_dp / papm1(:,:) * zesi_2d(:,:))


  ! --- 1.4 - Liquid phase ----------------------------------------------------
  ! ---------------------------------------------------------------------------

  ! Conditions for the formation of new cloud droplets
  ll1_2d(:,:) = &
       (ibas(:,:) > 0)      .AND. &  ! cloud base level
       (ptm1(:,:) > cthomi) .AND. &  ! T > T_hom
       ((zcdnc(:,:) <= cdncmin) .OR. & ! CDNC >= CDNC_min (10 cm-3)
       (pxlte(:,:) > 0._dp)     .OR. & ! d(LWC)/dt > 0
       (paclc(:,:) > cloud_tm1(1:kproma,:,jrow)) .OR. & ! d(cover)/dt > 0
       (zsusatw_2d(:,:) >  zeps)) ! water supersaturation > 0 (ratio > 1)

  ! Convert the aerosol activation into number of newly formed cloud droplets
  ! First computes newly formed cloud droplets at cloud bases:
  ztmp(:,:) = pcdncact(:,:) - zcdnc(:,:)  ! CDNC(ARG) - CDNC(preexisting) [m-3]
  ztmp(:,:) = MAX(0._dp, ztmp(:,:))
  zqlnuc(:,:) = MERGE(ztmp(:,:), 0._dp, ll1_2d(:,:))

  ztmp(:,:) = zcdnc(:,:) + zqlnuc(:,:)
  zcdnc(:,:) = MERGE(ztmp(:,:), zcdnc(:,:), ll1_2d(:,:))
  qnuc(1:kproma,:,jrow) = &
       qnuc(1:kproma,:,jrow) + zdt * zqlnuc(1:kproma,:)  ! nucl. rate [m-3 s-1]

  ! Then computes newly formed cloud droplets above cloud base by assuming that
  ! the number of nucleated cloud droplets is constant above cloud base
  ! (adiabatic parcel theory)
  DO jk = ktdia, klev  ! vertical loop (bottom to top)
     DO jl = 1, kproma

        ! Sets the level index to either relevant cloud base or itself
        jkk = MAX(jk,icl_minusbas(jl,jk))

        ! In each cloud level, hold the number of newly formed droplets at base
        zqlnuc_bas(jl,jk) = zqlnuc(jl,jkk)
        zcdnc_bas(jl,jk) = zcdnc(jl,jkk)
     ENDDO
  ENDDO ! vertical loop

  ll1_2d(:,:) = (icl_minusbas(:,:) > 0) .AND. & ! all cloud levels above base
                (zqlnuc_bas(:,:)   > 0._dp)     ! newly cloud droplets at base

  zqlnuc(:,:) = MERGE(zqlnuc_bas(:,:), zqlnuc(:,:), ll1_2d(:,:))
  ztmp(1:kproma,:) = qnuc(1:kproma,:,jrow) + &
       zdt * (zcdnc_bas(1:kproma,:) - zcdnc(1:kproma,:))
  qnuc(1:kproma,:,jrow) = &
       MERGE(ztmp(1:kproma,:), qnuc(1:kproma,:,jrow), ll1_2d(:,:))
  zcdnc(:,:) = MERGE(zcdnc_bas(:,:), zcdnc(:,:), ll1_2d(:,:))

  ! Correct for inconsistencies relative to phases and temperature
  ! Refined correction for crediting to ice *only* when cloud cover and
  ! LWC are non-zero (see https://redmine.hammoz.ethz.ch/issues/477)

  DO jk = ktdia, klev  ! vertical loop

     ll1_1d(:) = (ptm1(:,jk) <= cthomi) .AND. (zcdnc(:,jk) > cdncmin)
     ll2_1d(:) = (paclc(:,jk) >= zepsec) .AND. &
                 (pxlm1(:,jk) + ztmst*pxlte(:,jk) >= zepsec)

     zcd2ic(:) = MERGE(zcdnc(:,jk), 0._dp, ll1_1d(:) .AND. ll2_1d(:))
     zcd2unphys(:) = MERGE(zcdnc(:,jk), 0._dp, ll1_1d(:) .AND. .not.ll2_1d(:))

     zcdnc(:,jk) = zcdnc(:,jk) - zcd2ic(:) - zcd2unphys(:)
     zicnc(:,jk) = zicnc(:,jk) + zcd2ic(:)

     zcdnc(:,jk) = MAX(zcdnc(:,jk),cqtmin)
     zicncq(:,jk) = zicnc(:,jk)

  ENDDO  ! vertical loop

  ! Number of cloud droplets for he detrained cloud water from conv. anvils
  DO jk = ktdia, klev
     iclbas(:,jk) = NINT(pcvcbot(:))
     ll1_2d(:,jk) = (jk == iclbas(:,jk))
  ENDDO

  ll2_2d(:,:) = (iclbas(:,:) > 0)

  IF (.not.ncvmicro) THEN

     zqlnuccvh(:,:) = &
          MERGE(cdncact_cv(1:kproma,:,jrow) - zcdnc(:,:), 0._dp, ll1_2d(:,:))

     zxtec(:,:) = pxtec(:,:)

     ll3_2d(:,:) = ll2_2d(:,:) .AND.           &
                   (zxtec(:,:) > 0._dp)  .AND. &
                   (ptm1(:,:)  > cthomi)

     ztmp1(:,:) = cdncact_cv(1:kproma,:,jrow) - zcdnc(:,:)

     DO jk = ktdia, klev
        DO jl = 1, kproma
           jkk = iclbas(jl,jk)
           jkk = MAX(1, jkk)
           ztmp2(jl,jk) = zqlnuccvh(jl,jkk)
        ENDDO
     ENDDO

     ztmp(:,:) = MIN(ztmp1(:,:),ztmp2(:,:))
     ztmp(:,:) = MAX(0._dp,ztmp(:,:))

  ELSE

     zxtec(:,:) = pxtecl(:,:)

     ll3_2d(:,:) = ll2_2d(:,:) .AND. (zxtec(:,:) > 0._dp)

     ztmp(:,:) = pxtecnl(:,:) - zcdnc(:,:)
     ztmp(:,:) = MAX(0._dp, ztmp(:,:))

  ENDIF
  zqlnuccv(:,:)  = MERGE(ztmp(:,:), 0._dp, ll3_2d(:,:))
  zlwc_detr(:,:) = MERGE(zxtec(:,:), 0._dp, ll3_2d(:,:))
  qnuc(1:kproma,:,jrow) = qnuc(1:kproma,:,jrow) + zdt * zqlnuccv(1:kproma,:)

  ll1_2d(:,:) = (paclc(:,:) >= zepsec) .AND. &
                (pxlm1(:,:) + ztmst * pxlte(:,:) >= zepsec)

  ztmp1(:,:) = MERGE(paclc(:,:), zdummy, ll1_2d(:,:))
  ztmp2(:,:) = (pxlm1(:,:) + ztmst * pxlte(:,:)) / ztmp1(:,:)

  zlwc_strat(:,:) = MERGE(ztmp2(:,:), 0._dp, ll1_2d(:,:))

  ll2_2d(:,:) = ll1_2d(:,:) .OR. ll3_2d(:,:)

  zlwc_tot(:,:) = zlwc_strat(:,:) + zlwc_detr(:,:)
  zlwc_tot(:,:) = MERGE(zlwc_tot(:,:), zdummy, ll2_2d(:,:))

  ! Weight stratiform and detrained contribs by their respective LWC:
  zcdnc(:,:) = (zcdnc(:,:) * zlwc_strat(:,:) + &
                zqlnuccv(:,:) * zlwc_detr(:,:)) / zlwc_tot(:,:)

  zcdnc(:,:) = MAX(zcdnc(:,:),cqtmin)


  ! --- 1.5 - Ice phase -------------------------------------------------------
  ! ---------------------------------------------------------------------------

  ! If cloud ice is detrained from convective clouds at temperatures below
  ! -35 deg C, the corrsponding crystal number concentration is otained with
  ! an empirical formulation of the mean volume radius

  DO jk = klev, ktdia, -1  ! reverse vertical loop (top to bottom)
     DO jl = 1, kproma

        ! Mean volume radius [um]
        ! Boudala et al. (Int. J. Climatol., 2002), Eq. (10b)
        ztc = ptm1(jl,jk) - tmelt
        ztc = MIN(0._dp, ztc)
        zrieff = 0.015_dp * ztc
        zrieff = EXP(zrieff)
        zrieff = 23.2_dp * zrieff
        zrieff = MAX(zrieff, 1.0_dp)

        ! Effective ice crystal radius (Moss et al., priv. comm., 1995) [um]
        zrih = 5113188.044_dp + 2809._dp * zrieff**3
        zrih = SQRT(zrih)
        zrih = -2261.236_dp + zrih
        zrid_2d(jl,jk) = 1.e-6_dp * zrih**(1._dp/3._dp)
        zrid_2d(jl,jk) = MAX(1.e-6_dp, zrid_2d(jl,jk)) ! [m]

     ENDDO
  ENDDO  ! reverse vertical loop

  znidetr(:,:) = 0._dp

  IF (ncvmicro) THEN

     ll_cv(:,:) = (pxteci(:,:) > 0._dp)
     ztmp(:,:) = pxtecni(:,:) - zicncq(:,:)
     ztmp(:,:) = MAX(ztmp(:,:), 0._dp)
     znidetr(:,:) = MERGE(ztmp(:,:), 0._dp, ll_cv(:,:))
     zicncq(:,:)  = zicncq(:,:) + znidetr(:,:)

  ELSE

     ll_cv(:,:) = (pxtec(:,:) > 0._dp) .AND. &
                  (ptm1(:,:)  < cthomi)

     ll1_2d(:,:) = (paclc(1:kproma,:) >= 0.01_dp)
     ztmp2(:,:)   = ztmst * pxtec(:,:) / MAX(paclc, 0.01_dp)
     ztmp2(:,:)   = MERGE(ztmp2(:,:), ztmst * pxtec(:,:), ll1_2d(:,:))

     ztmp1(:,:) = 0.75_dp * zrho(:,:) * ztmp2(:,:) / &
          (api * zrhoice * zrid_2d(:,:)**3)  ! Eq. (6) documentation
     ztmp1(:,:) = MAX(ztmp1(:,:), 0._dp)

     znidetr(:,:) = MERGE( ztmp1(:,:), 0._dp, ll_cv(:,:) )
     zicncq1(:,:) = 0._dp
     znidetr1(:,:)= 0._dp

     ll2_2d(:,:) = &
          ((pxim1(:,:) + ztmst * pxite(:,:) + ztmst * pxtec(:,:)) >= zepsec)

     ! Total IWC
     ztmp(:,:) = (pxim1(:,:) + ztmst * pxite(:,:) + ztmst * pxtec(:,:)) / &
          MAX(paclc(:,:), 0.01_dp)
     ztmp(:,:) = &
          MERGE(ztmp(:,:), &
                (pxim1(:,:) + ztmst * pxite(:,:) + &
                ztmst * pxtec(:,:)), ll1_2d(:,:))

     ! Stratiform IWC
     ztmp1(:,:) = &
          (pxim1(:,:) + ztmst * pxite(:,:)) / MAX(paclc(:,:), 0.01_dp)
     ztmp1(:,:)  = MERGE(ztmp1, (pxim1(:,:) + ztmst * pxite(:,:)), ll1_2d(:,:))

     zicncq1(:,:) = zicncq(:,:) * ztmp1(:,:) / MAX(ztmp(:,:), zepsec)
     zicncq1(:,:) = MERGE(zicncq1(:,:), cqtmin, ll2_2d(:,:))  ! IWC-weigthed

     znidetr1(:,:) = znidetr(:,:) * ztmp2(:,:) / MAX(ztmp(:,:), zepsec)
     znidetr1(:,:) = MERGE(znidetr1(:,:), 0._dp, ll2_2d(:,:))  ! IWC-weighted

     zicncq(:,:)   = MAX((zicncq1(:,:) + znidetr1(:,:)), cqtmin)
     znidetr(:,:)  = znidetr1(:,:)
     zicnc(:,:)    = zicncq1(:,:)

     ! Calculate mean-volume radius of all currently existing IC
     ztmp(:,:) = ztmp(:,:) * 0.75_dp * zrho(:,:)
     ll_ice(:,:)  = (ll2_2d(:,:) .AND. (zicncq(:,:) > zicemin))
     ztmp1(:,:)   = api * zrhoice * MAX(zicncq(:,:), zicemin)

     zri_crt(:,:) = MERGE((ztmp(:,:) / ztmp1(:,:)), 1.e-18_dp , ll_ice(:,:))
     zri_crt(:,:) = zri_crt(:,:)**(1./3.)
     zri_crt(:,:) = MAX(zri_crt(:,:), 1.e-6_dp)  ! [m]

  ENDIF

  ! Aerosol number concentration for homogeneous nucleation
  zascs = 0._dp
  zascs(1:kproma,:) = naersol(1:kproma,:,jrow)
  zascs(:,:) = MAX(zascs(:,:), 10.E6_dp*zrho(:,:))  ! limit to 10.e6 [1/kg]

  ! Calculate vertical velocity for the cirrus formation
  zsusatix_2d(1:kproma,:) = sice(1:kproma,:,jrow)  ! ice supersaturation

  ztmp1(:,:) = -100._dp * scale_v_ls_cirrus * g_rcp * pvervel(:,:) * &
       zrho_rcp(:,:) ! [cm/s]
  ztmp2(:,:) = 100._dp * scale_v_tke_cirrus * SQRT(ptkem1(:,:)) ! [cm/s]
  ztmp2(:,klev) = 0._dp

  ! Diagnostic output
  vervel_ls(1:kproma,:,jrow) = ztmp1(1:kproma,:)
  vervel_tke(1:kproma,:,jrow) = ztmp2(1:kproma,:)

  SELECT CASE(vervel_scheme)
  CASE(1)

     ! Large-scale + TKE
     zvervx_2d(:,:) = ztmp1(:,:) + ztmp2(:,:)

  CASE(2)

     CALL calc_vervel_orogw(kproma, klev, l_z, ptm1, papm1, pqm1, &
          zsusatix_2d, ampl_gwd, w_gwd_kpr, itest, wnucl)

     wnucl(:,:) = scale_v_orogw_cirrus * wnucl(:,:)

     ! Large scale + TKE or OROGW (depending on orogr. mask)
     ll2_1d(:) = (itest(:).eq.1)
     DO jk = 1, klev
        zvervx_2d(:,jk) = &
             ztmp1(:,jk) + MERGE(wnucl(:,jk), ztmp2(:,jk), ll2_1d(:))
     ENDDO

     ! Diagnostic output
     vervel_gw(1:kproma,:,jrow) = wnucl(1:kproma,:)

  CASE(3)

     CALL calc_vervel_orogw(kproma, klev, l_z, ptm1, papm1, pqm1, &
          zsusatix_2d, ampl_gwd, w_gwd_kpr, itest, wnucl)

     wnucl(:,:) = scale_v_orogw_cirrus * wnucl(:,:)

     ! Large-scale + TKE + OROGW
     zvervx_2d(:,:) = ztmp1(:,:) + ztmp2(:,:) + wnucl(:,:)

     ! Diagnostic output
     vervel_gw(1:kproma,:,jrow) = wnucl(1:kproma,:)

  CASE(4)

     ! Large-scale + Penner
     zvervx_2d(:,:) = ztmp1(:,:) + pvervel_p18(1:kproma,:)

     ! Diagnostic output for pvervel_p18 is saved in cloud_si

  CASE DEFAULT

     WRITE(*,*) "vervel_scheme ranges from 1 to 3 in messy_cloud_Kuebbeler14"
     STOP

  END SELECT


  ! --- 1.6 - Cirrus clouds parametrization -----------------------------------
  ! ---------------------------------------------------------------------------

  SELECT CASE(nicnc)
  CASE(1)  ! Scheme by Lohmann (J. Atmos. Sci., 2002)
     ll_ice(:,:) = (sice(1:kproma,:,jrow) > 0._dp) .AND. (ptm1(:,:) < cthomi)

     ztmp1(:,:) = 0.75_dp * zrho(:,:) * sice(1:kproma,:,jrow) * &
          zqsi_2d(:,:) / (api * zrhoice * zrid_2d(:,:)**3) - zicncq(:,:)
     ztmp2(:,:) = zascs(:,:) - zicncq(:,:)
     ztmp(:,:) = MIN(ztmp1(:,:), ztmp2(:,:))
     ztmp(:,:) = MAX(ztmp(:,:), 0._dp)

     zninucl(:,:) = MERGE(ztmp(:,:), 0._dp, ll_ice(:,:))
     zicncq(:,:)   = zicncq(:,:) + zninucl(:,:)

  CASE(2)  ! Scheme by Kaercher et al. (J. Geophys. Res., 2006)

     ! Supersaturation with respect to ice is allowed, thus the depositional
     ! growth equation needs to be solved

    ! Check monodisperse assumption
     IF (nfrzmod.NE.1) THEN
        WRITE(*,*) "NFRZMOD must be 1 in messy_cloud_Kuebbeler14 for nicnc=2"
        STOP
     ENDIF

     ! Total number of ice modes
     ntype  = nexti + nfrzaer  ! pre-existing ice + ice modes

     ! Initialize local arrays
     zapnx(:,:,:) = 0._dp
     zaprx(:,:,:) = 0._dp
     zapsx(:,:,:) = 0._dp
     scrhet(:,:)  = 0._dp
     zcaer(:,:,:) = 0.0_dp
     zraer(:,:,:) = 0.0_dp
     znicex(:,:) = 0._dp
     znicex_3d(:,:,:) = 0._dp
     zri_3d(:,:,:) = 1.0e-6_dp

     ! Default active fraction for dust deposition
     ! Moehler et al. (Atmos. Chem. Phys., 2006)
     ll1_2d(:,:) = ( ptm1(:,:) <= 238._dp )
     ll2_2d(:,:) = ( ptm1(:,:) <= 220._dp )

     ztmp1(:,:) = MERGE(0.5_dp, 0._dp, ll1_2d(:,:))
     ztmp1(:,:) = MERGE(2._dp, ztmp1(:,:), ll2_2d(:,:))
     ztmp2(:,:) = MERGE(1.2_dp, 0._dp, ll1_2d(:,:))
     ztmp2(:,:) = MERGE(1.1_dp, ztmp2(:,:), ll2_2d(:,:))

     ll3_2d(:,:) = (zsusatix_2d(:,:) .ge. 0._dp)
     ztmp(:,:) = (zsusatix_2d(:,:) + 1._dp) - ztmp2(:,:)
     ztmp(:,:) = MERGE(ztmp(:,:), 0._dp, ll3_2d(:,:))

     f_du_active(:,:) = EXP(ztmp1(:,:) * ztmp(:,:)) - 1._dp
     f_du_active(:,:) = MAX(f_du_active(:,:), 0.001_dp)
     f_du_active(:,:) = MIN(f_du_active(:,:), 1._dp)

     ! Read INP properties
     DO ntt = 1, ninp
        inp_name(ntt) = inp_properties(ntt)%name
        inp_Scrit(ntt) = inp_properties(ntt)%Scrit
        inp_f_active(ntt) = inp_properties(ntt)%f_active
     ENDDO

     ! Set INP properties for each freezing mode and call the parametrization
     DO jk = klev, ktdia, -1  ! reverse vertical loop

        DO ntt = 1, ninp

           ! Critical supersaturation
           scrhet(1,ntt) = inp_Scrit(ntt)
           scrhet(2,ntt) = inp_Scrit(ntt)

           SELECT CASE(TRIM(inp_name(ntt)))

           CASE('DUdep')

              ! Special case: Scrit is temperature dependent
              ! Moehler et al.(Atmos.Chem. Phys., 2006)
              scrhet(2,ntt) = 1.1_dp  ! T < 220 K

              ! Dust active fraction (use default if not provided)
              IF (inp_f_active(ntt).ge.0._dp) THEN
                 f_du_active(:,jk) = inp_f_active(ntt)
              ENDIF

              ! Number [cm-3]
              zapnx(:,ntt,1) = 1.e-6_dp * f_du_active(:,jk) * &
                   (nduinsolai(1:kproma,jk,jrow) + &
                    nduinsolci(1:kproma,jk,jrow) - &
                    zicncq(:,jk))
              zapnx(:,ntt,1) = MAX(zapnx(:,ntt,1), 0._dp)

              ! Radius [cm]
              zaprx(:,ntt,1) = 100._dp * mode(iacci)%wetrad(1:kproma,jk)
              zaprx(:,ntt,1) = MAX(zaprx(:,ntt,1), mode(iacci)%crdiv)

              ! Sigma
              zapsx(:,ntt,1) = 1.01_dp

           CASE('DUimm')

              ! Number [cm-3]
              zapnx(:,ntt,1) = 1.e-6_dp * inp_f_active(ntt) * &
                   (ndusol_cirrus(1:kproma,jk,jrow) - zicncq(:,jk))
              zapnx(:,ntt,1) = MAX(zapnx(:,ntt,1), 0._dp)

              ! Radius [cm]
              zaprx(:,ntt,1) = 100._dp * mode(iaccm)%wetrad(1:kproma,jk)
              zaprx(:,ntt,1) = MAX(zaprx(:,ntt,1), mode(iaccm)%crdiv)

              ! Sigma
              zapsx(:,ntt,1) = 1.01_dp

           CASE('BC')

              ! Number [cm-3]
              zapnx(:,ntt,1) = 1.e-6_dp * inp_f_active(ntt) * &
                   (nbcsol_cirrus(1:kproma,jk,jrow) + &
                    nbcinsol(1:kproma,jk,jrow) - &
                    zicncq(:,jk))
              zapnx(:,ntt,1) = MAX(zapnx(:,ntt,1), 0._dp)

              ! Radius [cm]
              zaprx(:,ntt,1) = 100._dp * mode(iaccm)%wetrad(1:kproma,jk)
              zaprx(:,ntt,1) = MAX(zaprx(:,ntt,1), mode(iaccm)%crdiv)

              ! Sigma
              zapsx(:,ntt,1) = 1.01_dp

           CASE('BCtag')

              ! Number [cm-3]
              zapnx(:,ntt,1) = 1.e-6_dp * inp_f_active(ntt) * &
                   (nbctagsol_cirrus(1:kproma,jk,jrow) + &
                    nbctaginsol(1:kproma,jk,jrow) - &
                    zicncq(:,jk))
              zapnx(:,ntt,1) = MAX(zapnx(:,ntt,1), 0._dp)

              ! Radius [cm]
              zaprx(:,ntt,1) = 100._dp * mode(iaccm)%wetrad(1:kproma,jk)
              zaprx(:,ntt,1) = MAX(zaprx(:,ntt,1), mode(iaccm)%crdiv)

              ! Sigma
              zapsx(:,ntt,1) = 1.01_dp

           CASE DEFAULT

              WRITE(*,*) "INP name",inp_name(ntt), &
                   "not implemented in messy_cloud_Kuebbeler14"
              STOP

           END SELECT

        ENDDO

        ! Homogeneous freezing always at the end
        zapnx(:,nfrzaer,1) = 1.e-6_dp * (zascs(:,jk) - zicncq(:,jk) )
        zapnx(:,nfrzaer,1) = MAX(zapnx(:,nfrzaer,1), 0._dp)
        zaprx(:,nfrzaer,1) = 100._dp * mode(iaccs)%wetrad(1:kproma,jk)
        zaprx(:,nfrzaer,1) = MAX(zaprx(:,nfrzaer,1), mode(iaccs)%crdiv)
        zapsx(:,nfrzaer,1) = 1.01_dp

        zap(:,jk) = SUM(zapnx(:,:,1), 2)

        ! Fictitious downdraft to account for preexisting ice crsystals
        DO  jl = 1, kproma

           IF (nexti.eq.1) THEN  ! with preexisting ice crystals

              ! Number and mean-volume radius of all currently existing ice
              ! crystals are transformed to mode 1 of znicex_3d(jl,jk, ntt=1)
              ! and zri_3d(jl,jk, ntt=1)
              IF (zicncq(jl,jk) .LE. cqtmin) THEN
                 znicex_3d(jl,jk,1) = 0._dp
                 ZRI_3d(jl,jk,1) = 1.e-6_dp  ! minimum radius 1 micron
              ELSE
                 znicex_3d(jl,jk,1) =  zicncq(jl,jk)
                 zri_3d(jl,jk,1) =  zri_crt(jl,jk)
              ENDIF

              ! Fictitious downdraft of preexisting ice crystals (zvice),
              ! based on calculations for depositional growth in subroutine
              ! XICE of Kaercher et al. (J. Geophys. Res., 2006)
              ! parametrization.
              ! Use CGS units here
              zicncq_vice(jl,jk) = zicncq(jl,jk) * 1.e-6_dp  ! [cm-3]
              zri_vice(jl,jk) = zri_crt(jl,jk) * 100.0_dp ! [cm]
              zvipr(jl,jk) = MIN(0.01_dp * papm1(jl,jk), 1013.25_dp) ! [hPa]
              zvit(jl,jk) = MAX(ptm1(jl,jk), 170.0_dp)  ! [K]
              zvit1(jl,jk) = 1.0_dp/zvit(jl,jk) ! [1/K]
              zvirhoa(jl,jk) = 0.35_dp * zvipr(jl,jk) * zvit1(jl,jk) ! dens.
              zvipi(jl,jk) = &
                   3.6e10_dp * EXP(-(6145._dp * zvit1(jl,jk))) ! sat. press.
              zvicis(jl,jk) = 7.24637701e+18_dp * zvipi(jl,jk) * &
                   zvit1(jl,jk)
              zvia1(jl,jk) = (0.601272523_dp * zvit1(jl,jk) - &
                   0.000342181855_dp)* zvit1(jl,jk) ! a1
              zvia2(jl,jk) = 1.0_dp / zvicis(jl,jk) ! a2
              zvia3(jl,jk) = &
                   1.49236645e-12_dp * zvit1(jl,jk) / zvipr(jl,jk) ! a3
              zvisi(jl,jk) = &
                   1._dp + zsusatix_2d(jl,jk) ! ice saturation ratio
              zvivsc(jl,jk) = 1.e-5_dp * &  ! air viscosity
                   (1.512_dp + 0.0052_dp * (zvit(jl,jk)-233.15_dp))
              zviflx(jl,jk) = 0.25_dp * 0.5_dp * &
                   SQRT(11713803.0_dp * zvit(jl,jk)) ! v(thermal) * alpha / 4
              zvib1(jl,jk) = zviflx(jl,jk) * 3.23e-23_dp * zvicis(jl,jk) * &
                   (zvisi(jl,jk) - 1._dp) ! b1
              zvib2(jl,jk) = zviflx(jl,jk) * 249.239822_dp * &
                   zvipr(jl,jk) * zvit1(jl,jk)**1.94_dp  ! b2
              IF (zicncq_vice(jl,jk) > 1.e-5_dp .AND. &
                   zri_vice(jl,jk) > 1e-4_dp) THEN
                 zvifall(jl,jk) = &
                      14._dp * zri_vice(jl,jk) * &
                      (1.3_dp/zvirhoa(jl,jk))**0.35_dp ! fall velocity of IC
                 zvirey(jl,jk) = &
                      zvirhoa(jl,jk) * 2._dp * zri_vice(jl,jk) * &
                      zvifall(jl,jk) / zvivsc(jl,jk) ! Reynolds number
                 zvivent(jl,jk) = &
                      1._dp + 0.229_dp * SQRT(zvirey(jl,jk)) ! ventil. effect
                 ! Freezing integral (Eq. 12)
                 zvidl(jl,jk) = 3.89051704e+23_dp * zicncq_vice(jl,jk) * &
                      zvivent(jl,jk) * zvib1(jl,jk) * &
                      zri_vice(jl,jk)**2._dp / &
                      (1._dp + zvib2(jl,jk) * zri_vice(jl,jk))
                 ! Fictitious downdraft (Eq. 13)
                 zvice(jl,jk) = &
                      (zvia2(jl,jk) + zvia3(jl,jk) * zvisi(jl,jk)) * &
                      zvidl(jl,jk) / (zvia1(jl,jk) * zvisi(jl,jk))
              ELSE
                 zvice(jl,jk) = 0.0_dp
              ENDIF
           ENDIF

           IF (nexti.eq.0) THEN  ! no preexisting ice crystals
              zvice(jl,jk) = 0._dp
           ENDIF

        ENDDO

        ! Limit supersaturation to 1 (i.e., saturation ratio to 2) to avoid
        ! instabilities in the parametrization (Kaercher, priv. comm., 2017)
        zsusatix_2d_lim(:,jk) = MIN(zsusatix_2d(:,jk), 1._dp)

        ! Call Kaercher parametrization main driver
        CALL ACIDRV(jk, kbdim, kproma, klev, nexti, nfrzaer, nfrzmod,    &
             zvervx_2d(:,jk), ztmst, papm1, ptm1, zsusatix_2d_lim(:,jk), &
             zapnx, zaprx, zapsx, scrhet, zcaer, zraer, znicex_3d,     &
             zri_3d, zvice(:,jk), cthomi, zsicirrus_2d(:,jk), &
             cirrus_status, cirrus_string)
        IF (cirrus_status /= 0) RETURN

     ENDDO ! reverse vertical loop

     ! Store output in diagnostic variables for each ice mode [cm-3]
     IF (nexti.eq.1) Nice_preex(1:kproma,:,jrow) = znicex_3d(:,:,1) * 1.e-6_dp
     DO ntt = 1, nfrzaer - 1
        SELECT CASE(TRIM(inp_name(ntt)))
        CASE('DUdep')
           Nice_DUdep(1:kproma,:,jrow) = znicex_3d(:,:,ntt+nexti) * 1.e-6_dp
        CASE('DUimm')
           Nice_DUimm(1:kproma,:,jrow) = znicex_3d(:,:,ntt+nexti) * 1.e-6_dp
        CASE('BC')
           Nice_BC(1:kproma,:,jrow) = znicex_3d(:,:,ntt+nexti) * 1.e-6_dp
        CASE('BCtag')
           Nice_BCtag(1:kproma,:,jrow) = znicex_3d(:,:,ntt+nexti) * 1.e-6_dp
        CASE DEFAULT
           WRITE(*,*) "INP name",inp_name(ntt), &
                "not implemented in messy_cloud_Kuebbeler14"
           STOP
        END SELECT
     ENDDO
     Nice_homog(1:kproma,:,jrow) = znicex_3d(:,:,nfrzaer+nexti) * 1.e-6_dp

     ! Calculate depositional growth of the multimodal IC

     ! Update ICNC
     zri_3d(:,:,:) = MAX(zri_3d(:,:,:), 1.e-6_dp)
     DO jk = klev, ktdia, -1

        IF (nexti.eq.0) THEN
           DO ntt = 1, ntype
              ztmp1(:,jk) = 1.e6_dp * zapnx(:,ntt,1)  ! [m-3]
              ztmp2(:,jk) = znicex_3d(:,jk,ntt)
              znicex_3d(:,jk,ntt) = MIN(ztmp1(:,jk), ztmp2(:,jk))
              znicex_3d(:,jk,ntt) = MAX(znicex_3d(:,jk,ntt), 0._dp)
           ENDDO
        ENDIF

        IF (nexti.eq.1) THEN

           ! Preexisting IC
           ztmp1(:,jk) = zicncq(:,jk)
           ztmp2(:,jk) = znicex_3d(:,jk,1)
           znicex_3d(:,jk,1) = MIN(ztmp1(:,jk), ztmp2(:,jk))
           znicex_3d(:,jk,1) = MAX(znicex_3d(:,jk,1), 0._dp)

           ! Other modes
           DO ntt = 1, nfrzaer
              ztmp1(:,jk) = 1.e6_dp * zapnx(:,ntt,1)  ! [m-3]
              ztmp2(:,jk) = znicex_3d(:,jk,ntt+1)
              znicex_3d(:,jk,ntt+1) = MIN(ztmp1(:,jk), ztmp2(:,jk))
              znicex_3d(:,jk,ntt+1) = MAX(znicex_3d(:,jk,ntt+1), 0._dp)
           ENDDO
        ENDIF
     ENDDO

     ! Calculate deposition rate taking the ventilation (zfre) into account
     IF (nexti.eq.1) THEN  ! preexisting IC
        DO ntt = 1, ntype  ! = nexti + nfrzaer

           ! Mass of a single IC [kg]
           ll_ice(:,:) = (znicex_3d(:,:,ntt) > cqtmin)
           ztmp(:,:) = 4._dp/3._dp * api * zri_3d(:,:,ntt)**3 * zrhoice
           ztmp(:,:) = MAX(ztmp(:,:), zmi)  ! minimum mass
           zmmean_2d(:,:) = MERGE(ztmp(:,:), zmmean_2d(:,:), ll_ice(:,:))

           ! Eq. (14) documentation
           ll1_2d(:,:) = (zmmean_2d(:,:) < 2.166E-9_dp)
           ll2_2d(:,:) = (zmmean_2d(:,:) >= 2.166E-9_dp) .AND. &
                (zmmean_2d(:,:) < 4.264E-8_dp)
           zalfased_2d(:,:) = MERGE(63292.4_dp, 8.78_dp, ll1_2d(:,:))
           zalfased_2d(:,:) = MERGE(329.75_dp, zalfased_2d(:,:), ll2_2d(:,:))
           zbetased_2d(:,:) = MERGE(0.5727_dp, 0.0954_dp, ll1_2d(:,:))
           zbetased_2d(:,:) = MERGE(0.3091_dp, zbetased_2d(:,:), ll2_2d(:,:))

           ! Fall velocity
           ! Spichtinger and Gierens (Atmos. Chem. Phys., 2009) - Eq. (18)
           zxifallmc_2d(:,:) = zfall * zalfased_2d(:,:) * &
                (zmmean_2d(:,:)**zbetased_2d(:,:)) * zaaa_2d(:,:)

           ztmp(:,:) = pqm1(:,:) - zqsi_2d(:,:)

           DO jk = ktdia, klev
              DO jl = 1,kproma
                 zdv = 2.21_dp / papm1(jl,jk)
                 zgtp = 1._dp/(zrho(jl,jk) * zastbsti(jl,jk))
                 zvth = SQRT(8._dp * zkb * ptm1(jl,jk) / (api * zxmw))
                 zb2 = 0.25_dp * zalpha * zvth / zdv
                 zfuchs = 1._dp / (1._dp + zb2 * zri_3d(jl,jk,ntt))
                 zre = 2._dp * zrho(jl,jk) * zri_3d(jl,jk,ntt) * &
                      zxifallmc_2d(jl,jk) / zviscos_2d(jl,jk)
                 zfre = 1._dp + 0.229_dp * SQRT(zre)
                 zfre = MERGE(zfre, 1._dp, ll_ice(jl,jk))
                 zqinucl(jl,jk) = 4._dp * api * zri_3d(jl,jk,ntt) * &
                      zsusatix_2d(jl,jk) * znicex_3d(jl,jk,ntt) * &
                      zfre * zgtp * zfuchs * zalpha * ztmst
                 zqinucl(jl,jk) = MIN(zqinucl(jl,jk), ztmp(jl,jk))
                 zqinucl(jl,jk) = MAX(zqinucl(jl,jk), -pxim1(jl,jk))
              ENDDO
           ENDDO
           zqinucl_3d(:,:,ntt)  = zqinucl(:,:)

        ENDDO
     ENDIF

     IF (nexti.eq.0) THEN  ! no preexisting IC

        ! Mass of a single IC [kg]
        ll_ice(:,:) = (pxim1(:,:) > 0._dp)  ! IWC > 0
        ztmp(:,:) = MAX(zicncq(:,:), zicemin)
        zicncq(:,:) = MERGE(ztmp(:,:), zicncq(:,:), ll_ice(:,:))
        ll1_2d(:,:) = (paclc(:,:) > 0.01_dp)
        ztmp1(:,:) = pxim1(:,:) / MAX(paclc, 0.01_dp)  ! in-cloud IWC
        ztmp1(:,:) = MERGE(ztmp1, pxim1(:,:) , ll1_2d )
        ztmp(:,:) = zrho(:,:) * ztmp1(:,:) / zicncq(:,:)  ! [kg]
        ztmp(:,:) = MAX(ztmp(:,:), zmi)
        zmmean_2d(:,:) = MERGE(ztmp(:,:), zmmean_2d(:,:), ll_ice(:,:))

        ! Eq. (14) documentation
        ll1_2d(:,:) = (zmmean_2d(:,:) < 2.166E-9_dp)
        ll2_2d(:,:) = (zmmean_2d(:,:) >= 2.166E-9_dp) .AND. &
             (zmmean_2d(:,:) < 4.264E-8_dp)
        zalfased_2d(:,:) = MERGE(63292.4_dp, 8.78_dp, ll1_2d(:,:))
        zalfased_2d(:,:) = MERGE(329.75_dp, zalfased_2d(:,:), ll2_2d(:,:))
        zbetased_2d(:,:) = MERGE(0.5727_dp, 0.0954_dp, ll1_2d(:,:))
        zbetased_2d(:,:) = MERGE(0.3091_dp, zbetased_2d(:,:), ll2_2d(:,:))

        ! Fall velocity
        ! Spichtinger and Gierens (Atmos. Chem. Phys., 2009) - Eq. (18)
        zxifallmc_2d(:,:) = zfall * zalfased_2d(:,:) * &
             (zmmean_2d(:,:)**zbetased_2d(:,:)) * zaaa_2d(:,:)

        ztmp(:,:) = pqm1(:,:) - zqsi_2d(:,:)

        DO jk = ktdia, klev
           DO jl = 1,kproma
              zdv = 2.21_dp / papm1(jl,jk)
              zgtp = 1._dp / (zrho(jl,jk) * zastbsti(jl,jk))
              zvth = SQRT(8._dp * zkb * ptm1(jl,jk) / (api*zxmw) )
              zb2 = 0.25_dp * zalpha * zvth / zdv
              zfuchs = 1._dp / (1._dp + zb2 * zrid_2d(jl,jk))
              zre = 2._dp * zrho(jl,jk) * zrid_2d(jl,jk) * &
                   zxifallmc_2d(jl,jk) / zviscos_2d(jl,jk)
              zfre = 1._dp + 0.229_dp * SQRT(zre)
              zfre = MERGE(zfre, 1._dp, ll_ice(jl,jk))
              zqinucl(jl,jk) = 4._dp * api * zrid_2d(jl,jk) * &
                   zsusatix_2d(jl,jk) * zicncq(jl,jk) * &
                   zfre * zgtp * zfuchs * zalpha * ztmst
              zqinucl(jl,jk) = MIN(zqinucl(jl,jk), ztmp(jl,jk))
              zqinucl(jl,jk) = MAX(zqinucl(jl,jk), -pxim1(jl,jk))
           ENDDO
        ENDDO
        zqinucl_3d(:,:,1)  = zqinucl(:,:)

        DO ntt = 1, nfrzaer

           ! Mass of a single IC [kg]
           ll_ice(:,:) = (znicex_3d(:,:,ntt) > cqtmin)
           ztmp(:,:) = 4_dp/3._dp *api * zri_3d(:,:,ntt)**3 * zrhoice
           ztmp(:,:) = MAX(ztmp(:,:), zmi)
           zmmean_2d(:,:) = MERGE(ztmp(:,:), zmmean_2d(:,:), ll_ice(:,:))
           ll1_2d(:,:) = (zmmean_2d(:,:) < 2.166E-9_dp)
           ll2_2d(:,:) = (zmmean_2d(:,:) >= 2.166E-9_dp) .AND. &
                (zmmean_2d(:,:) < 4.264E-8_dp)

           ! Eq. (14) documentation
           zalfased_2d(:,:) = MERGE(63292.4_dp, 8.78_dp, ll1_2d(:,:))
           zalfased_2d(:,:) = MERGE(329.75_dp, zalfased_2d(:,:), ll2_2d(:,:))
           zbetased_2d(:,:) = MERGE(0.5727_dp, 0.0954_dp, ll1_2d(:,:))
           zbetased_2d(:,:) = MERGE(0.3091_dp, zbetased_2d(:,:), ll2_2d(:,:))

           ! Fall velocity
           ! Spichtinger and Gierens (Atmos. Chem. Phys., 2009) - Eq. (18)
           zxifallmc_2d(:,:) = zfall * zalfased_2d(:,:) * &
                (zmmean_2d(:,:)**zbetased_2d(:,:)) * zaaa_2d(:,:)

           ztmp(:,:) = pqm1(:,:) - zqsi_2d(:,:)

           DO jk = ktdia, klev
              DO jl = 1,kproma
                 zdv = 2.21_dp / papm1(jl,jk)
                 zgtp = 1._dp / (zrho(jl,jk) * zastbsti(jl,jk))
                 zvth = SQRT(8._dp * zkb * ptm1(jl,jk) / (api * zxmw))
                 zb2 = 0.25_dp * zalpha * zvth / zdv
                 zfuchs = 1._dp / (1._dp + zb2 * zri_3d(jl,jk,ntt))
                 zre = 2._dp * zrho(jl,jk) * zri_3d(jl,jk,ntt) * &
                      zxifallmc_2d(jl,jk) / zviscos_2d(jl,jk)
                 zfre = 1._dp + 0.229_dp * SQRT(zre)
                 zfre = MERGE(zfre, 1._dp, ll_ice(jl,jk))
                 zqinucl(jl,jk) = 4._dp * api * zri_3d(jl,jk,ntt) * &
                      zsusatix_2d(jl,jk) * znicex_3d(jl,jk,ntt) * &
                      zfre * zgtp * zfuchs * zalpha * ztmst
                 zqinucl(jl,jk) = MIN(zqinucl(jl,jk), ztmp(jl,jk))
                 zqinucl(jl,jk) = MAX(zqinucl(jl,jk), -pxim1(jl,jk))
              ENDDO
           ENDDO
           zqinucl_3d(:,:,ntt+1)  = zqinucl(:,:)

        ENDDO
     ENDIF

  CASE DEFAULT

     WRITE(*,*) "nicnc can be either 1 or 2 in messy_cloud_Kuebbeler14"
     STOP

  END SELECT  ! ice-cirrus parametrization selection


  ! Cloud cover enhancement due to orographic waves
  IF (vervel_scheme == 2 .or. vervel_scheme == 3) THEN
     DO jk = ktdia, klev
        ll1_2d(:,jk) = ((itest(:).eq.1) .AND. &
                       (w_gwd_kpr(:,jk).gt.0.0_dp) .AND. &
                       (zqinucl(:,jk).gt.ccwmin) .AND.(nicnc > 1))
        ll2_2d(:,jk) = (zsusatix_2d(:,jk).gt.0._dp)
        ztmp1_1d(:) = MAX(sqcst_2d(:) * 111325._dp * 1.8_dp, 0.1_dp)

        ! Cloud cover increment due to orographic waves
        ztmp(:,jk) = MIN(2._dp * l_z(:) / ztmp1_1d(:), 1._dp)
        ztmp(:,jk) = MERGE(1._dp, ztmp(:,jk), ll2_2d(:,jk))
        paclc(:,jk) = MERGE(paclc(:,jk)+ztmp(:,jk), paclc(:,jk), ll1_2d(:,jk))
        paclc(:,jk) = MIN(paclc(:,jk), 1._dp)
     ENDDO
  ENDIF

  ! ===========================================================================
  ! ===========================================================================

  DO 831 jk = ktdia, klev


     ! === 2 - SET TO ZERO LOCAL TENDENCIES ===================================
     ! ========================================================================

     ! Pointer to the rain and snowflux through the bottom of each layer
     zrfl => pfrain(:,jk)
     zsfl => pfsnow(:,jk)
     IF (jk > 1) THEN
        zrfl(1:kproma) = pfrain(1:kproma, jk-1)
        zsfl(1:kproma) = pfsnow(1:kproma, jk-1)
     ENDIF

     zevp => prevap(:,jk)  ! pointer to rain evaporation
     zsub => pssubl(:,jk)  ! pointer to snow sublimation
     zrpr => prate_r(:,jk) ! pointer to rain production

     ! Initialisations
     zcnd(:)       = 0.0_dp
     zdep(:)       = 0.0_dp
     zrpr(:)       = 0.0_dp
     zspr(:)       = 0.0_dp
     zimlt(:)      = 0.0_dp
     zevp(:)       = 0.0_dp
     zsub(:)       = 0.0_dp
     zximlt(:)     = 0.0_dp
     zxisub(:)     = 0.0_dp
     zsmlt(:)      = 0.0_dp
     zsacl(:)      = 0.0_dp
     zgenti(:)     = 0.0_dp
     zgentl(:)     = 0.0_dp
     zxievap(:)    = 0.0_dp
     zxlevap(:)    = 0.0_dp
     zvartg(:)     = 0.0_dp
     zconvvar(:)   = 0.0_dp
     zconvskew(:)  = 0.0_dp
     zturbvar(:)   = 0.0_dp
     zturbskew(:)  = 0.0_dp
     zmicroskew(:) = 0.0_dp
     zxvarte(:)    = 0.0_dp
     zxskewte(:)   = 0.0_dp
     zzdrr_1d(:)   = 0.0_dp
     zzdrs_1d(:)   = 0.0_dp
     zrprn(:)      = 0.0_dp
     zfrln(:)      = 0.0_dp
     zsacln(:)     = 0.0_dp
     zsprn(:,jk)   = 0.0_dp

     ! Latent heats
     ztmp(:,jk) = MAX(pqm1(:,jk), 0._dp)
     ztmp(:,jk) = 1._dp / (cpd + zcons1 * ztmp(:,jk))
     zlvdcp(:) = alv * ztmp(:,jk)
     zlsdcp(:) = als * ztmp(:,jk)

     zsaut2_1d(:)  = 0.0_dp

     IF (jk > 1) THEN

        ! --- 2.1 - Transfer cirrus IC with R >= 100 micron to snow -----------
        ! ---------------------------------------------------------------------

        ! Levkov et al. (Beitr. Phys. Atmosph., 1992)
        IF (l_ic2snow) THEN

           IF (nexti.eq.1) THEN

              ! Use zqinucl: mmr of all currently existing IC
              DO ntt = 1, ntype
                 ll1_1d(:) = (zri_3d(:,jk,ntt) >= 1.0e-4_dp)  ! R >= 100 um

                 ! IC mass flux in this mode [kg/m2/s]
                 ! FIX-ME: add ice mass of all IC modes (?)
                 ztmp(:,jk) = MAX(zqinucl_3d(:,jk,ntt), 0._dp)
                 ztmp1_1d(:) = zdp_2d(:,jk) * ztmp(:,jk) / (g * ztmst)

                 zsaut2_1d(:) = &
                      MERGE(zsaut2_1d(:) + &
                            (ztmp1_1d(:) * g * ztmst / zdp_2d(:,jk)), &
                            zsaut2_1d(:), ll1_1d(:))

                 ztmp2_1d(:) = zsfl(:) + ztmp1_1d(:)

                 ! Mass
                 zsfl(:) = MERGE(ztmp2_1d(:), zsfl(:), ll1_1d(:)) ! Add to snow
                 ztmp2_1d(:) = zxiflux(:) - ztmp1_1d(:)
                 zxiflux(:) = &
                      MERGE(ztmp2_1d(:), zxiflux(:), ll1_1d(:)) ! Remove from IC

                 ! IC number flux in this mode
                 ztmp1_1d(:) = zdp_2d(:,jk) * znicex_3d(:,jk,ntt) * &
                      zrho_rcp(:,jk) / (g * ztmst)

                 ! Number
                 ztmp2_1d(:) = zxifluxn(:) - ztmp1_1d(:)
                 zxifluxn(:) = MERGE(ztmp2_1d(:), zxifluxn(:), ll1_1d(:))

                 ! Set IC mass mixing ratio to zero if R > 100 micron
                 zqinucl_3d(:,jk,ntt) = &
                      MERGE(0._dp, zqinucl_3d(:,jk,ntt), ll1_1d(:))
                 znicex_3d(:,jk,ntt) = &
                      MERGE(0._dp, znicex_3d(:,jk,ntt), ll1_1d(:))
                 zri_3d(:,jk,ntt) = &
                      MERGE(0.0_dp, zri_3d(:,jk,ntt), ll1_1d(:))
              ENDDO

           ENDIF

           IF (nexti.eq.0) THEN
              ll1_1d(:)     = ( zrid_2d(:,jk) >= 1.0e-4_dp )  ! R >= 100 um

              ! IC mass flux in this mode [kg/m2/s]
              ! FIX-ME: add ice mass of all IC modes (?)
              ztmp(:,jk) = MAX(zqinucl_3d(:,jk,1), 0._dp)
              ztmp1_1d(:) = zdp_2d(:,jk) * ztmp(:,jk) / (g * ztmst)

              zsaut2_1d(:) = &
                   MERGE(zsaut2_1d(:) + &
                         (ztmp1_1d(:) * g * ztmst / zdp_2d(:,jk)), &
                         zsaut2_1d(:), ll1_1d(:))

              ztmp2_1d(:) = zsfl(:) + ztmp1_1d(:)

              ! Mass
              zsfl(:) = MERGE(ztmp2_1d(:), zsfl(:), ll1_1d(:)) ! Add to snow
              ztmp2_1d(:) = zxiflux(:) - ztmp1_1d(:)
              zxiflux(:) = &
                   MERGE(ztmp2_1d(:), zxiflux(:), ll1_1d(:)) ! Remove from IC

              ! IC number flux in this mode
              ztmp1_1d(:) = zdp_2d(:,jk) * zicncq(:,jk) * &
                   zrho_rcp(:,jk) / (g * ztmst)

              ! Number
              ztmp2_1d(:) = zxifluxn(:) - ztmp1_1d(:)
              zxifluxn(:) = MERGE(ztmp2_1d(:), zxifluxn(:), ll1_1d(:))

              ! Set IC mass mixing ratio to zero if R > 100 micron
              zqinucl_3d(:,jk,1) = &
                   MERGE(0._dp, zqinucl_3d(:,jk,1), ll1_1d(:))
              zicncq(:,jk) = &
                   MERGE(0.0_dp, zicncq(:,jk), ll1_1d(:))
              zrid_2d(:,jk) = &
                   MERGE(0.0_dp, zrid_2d(:,jk), ll1_1d(:))

              DO ntt = 1, nfrzaer
                 ll1_1d(:) = (zri_3d(:,jk,ntt+1) >= 1.0e-4_dp)

                 ! IC mass flux in this mode [kg/m2/s]
                 ! FIX-ME: add ice mass of all IC modes (?)
                 ztmp(:,jk) = MAX(zqinucl_3d(:,jk,ntt+1), 0._dp)
                 ztmp1_1d(:) = zdp_2d(:,jk) * ztmp(:,jk) / (g * ztmst)

                 zsaut2_1d(:) = &
                      MERGE (zsaut2_1d(:) + &
                             (ztmp1_1d(:) * g * ztmst / zdp_2d(:,jk)), &
                             zsaut2_1d(:), ll1_1d(:))

                 ztmp2_1d(:) = zsfl(:) + ztmp1_1d(:)

                 ! Mass
                 zsfl(:) = MERGE(ztmp2_1d(:), zsfl(:), ll1_1d(:)) ! Add to snow
                 ztmp2_1d(:) = zxiflux(:) - ztmp1_1d(:)
                 zxiflux(:) = &
                      MERGE(ztmp2_1d(:), zxiflux(:), ll1_1d(:)) ! Remove from IC

                 ! IC number flux in this mode
                 ztmp1_1d(:) = zdp_2d(:,jk) * znicex_3d(:,jk,ntt+1) * &
                      zrho_rcp(:,jk) / (g * ztmst)

                 ! Number
                 ztmp2_1d(:) = zxifluxn(:) - ztmp1_1d(:)
                 zxifluxn(:) = MERGE(ztmp2_1d(:), zxifluxn(:), ll1_1d(:))

                 ! Set IC mass mixing ratio to zero if R > 100 micron
                 zqinucl_3d(:,jk,ntt+1) = &
                      MERGE(0.0_dp, zqinucl_3d(:,jk,ntt+1), ll1_1d(:))
                 znicex_3d(:,jk,ntt+1) = &
                      MERGE(0.0_dp, znicex_3d(:,jk,ntt+1), ll1_1d(:))
                 zri_3d(:,jk,ntt+1) = &
                      MERGE(0.0_dp, zri_3d(:,jk,ntt+1), ll1_1d(:))

              ENDDO
           ENDIF

        ENDIF ! l_ic2snow

        ! --- 2.2 - Transformation of multimodal ICNC and Rice to unimodal ----
        ! ---------------------------------------------------------------------
        znicex_tmp(:,jk)   = 0.0_dp
        zri_tmp(:,jk)      = 0.0_dp
        zqinucl_tmp(:,jk)  = 0.0_dp

        DO jl = 1, kproma

           ! Sum over all ntypes (pre-existing + freezing modes)
           DO ntt = 1, ntype
              znicex_tmp(jl,jk) = znicex_tmp(jl,jk) + znicex_3d(jl,jk,ntt)
              zri_tmp(jl,jk) = zri_tmp(jl,jk) + &
                   (znicex_3d(jl,jk,ntt) * zri_3d(jl,jk,ntt))
           ENDDO

           IF (nexti.eq.1) THEN
              DO ntt = 1, ntype
                 zqinucl_tmp(jl,jk) = &
                      zqinucl_tmp(jl,jk) + zqinucl_3d(jl,jk,ntt)
              ENDDO
           ENDIF

           IF (nexti.eq.0) THEN
              DO ntt = 1, nfrzaer + 1
                 zqinucl_tmp(jl,jk) = &
                      zqinucl_tmp(jl,jk) + zqinucl_3d(jl,jk,ntt)
              ENDDO
           ENDIF

           ! Number-weigthed radius
           IF (znicex_tmp(jl,jk).gt.0._dp) THEN
              zri_tmp(jl,jk) = zri_tmp(jl,jk) / znicex_tmp(jl,jk)
           ELSE
              zri_tmp(jl,jk) = 1.e-6_dp  ! 1 [um]
           ENDIF

           zri(jl,jk) = zri_tmp(jl,jk)  ! [m]
           znicex(jl,jk) = znicex_tmp(jl,jk)  ! [m-3]
           zqinucl(jl,jk) = zqinucl_tmp(jl,jk)  ! [kg/kg]
        ENDDO

        ! Update znicex, zri, zninucl, zquinucl
        zri(:,jk)=MAX(zri(:,jk), 1.e-6_dp)

        ll_ice(:,jk) = (zsusatix_2d(:,jk) > 0._dp) .AND. &
                       (ptm1(:,jk) < (cthomi-1._dp))

        IF (nexti.eq.1)  ztmp1(:,jk) = 1.e6_dp * zap(:,jk) + zicncq(:,jk)
        IF (nexti.eq.0)  ztmp1(:,jk) = 1.e6_dp * zap(:,jk)

        ztmp2(:,jk) = znicex(:,jk)
        ztmp(:,jk) = MIN(ztmp1(:,jk),ztmp2(:,jk))  ! min(N_IC, N_IN)
        ztmp(:,jk) = MAX(ztmp(:,jk), 0._dp)

        IF (nexti.eq.1) &
             zninucl(:,jk) = &
             MERGE(ztmp(:,jk) - zicncq(:,jk), 0._dp, ll_ice(:,jk))
        IF (nexti.eq.0) &
             zninucl(:,jk) = MERGE(ztmp(:,jk), 0._dp, ll_ice(:,jk))
        zninucl(:,jk) = Max(zninucl(:,jk), 0._dp)

        zqinucl(:,jk) = MERGE(zqinucl(:,jk), 0._dp, ll_ice(:,jk))
        zicncq(:,jk)  = zicncq(:,jk) + zninucl(:,jk)
        zicncq(:,jk)  = Max(zicncq(:,jk), cqtmin)

        ! Set CDNC and ICNC to minimum if nucleation is not strong enough
        ll1(:,jk) = (paclc(:,jk) >= zepsec ) .AND. (ptm1(:,jk) > cthomi)
        ztmp(:,jk) = MAX(zcdnc(:,jk), cdncmin)
        zcdnc(:,jk) = MERGE(ztmp(:,jk), zcdnc(:,jk), ll1(:,jk))

        ll1(:,jk) = (paclc(:,jk) >= zepsec ) .AND. (ptm1(:,jk) < tmelt)
        ztmp(:,jk) = MAX(zicncq(:,jk), zicemin)
        zicncq(:,jk) = MERGE(ztmp(:,jk), zicncq(:,jk), ll1(:,jk))

        ! === 3 - MODIFICATION OF INCOMING PRECIPITATION FLUXES ===============
        ! =====================================================================

        ! --- 3.1 - Melting of snow and ice -----------------------------------
        ! ---------------------------------------------------------------------

        ztmp1_1d(:) = ptm1(:,jk) - tmelt

        ll1_1d(:) = (ztmp1_1d(:) > 0._dp)

        ztmp1_1d(:) = MAX(0._dp, ztmp1_1d(:))
        ztmp1_1d(:) = zcons2 * ztmp1_1d(:) * zdp_2d(:,jk) / &
             (zlsdcp(:) - zlvdcp(:))

        ztmp2_1d(:) = zxsec * zsfl(:)
        ztmp2_1d(:) = MIN(ztmp2_1d(:), ztmp1_1d(:))
        pimelt(:,jk) = ztmp2_1d(:)  ! save in channel output

        zrfl(:)  = zrfl(:) + ztmp2_1d(:)
        zsfl(:)  = zsfl(:) - ztmp2_1d(:)

        ! Eq. (8) documentation
        zsmlt(:) = ztmst * g * ztmp2_1d(:) / zdp_2d(:,jk)

        ztmp2_1d(:) = zxsec * zxiflux(:)
        ztmp2_1d(:) = MIN(ztmp2_1d(:), ztmp1_1d(:))

        ll2_1d(:) = (zxiflux(:) > zepsec)

        ztmp3_1d(:) = zxifluxn(:) * ztmp2_1d(:) / MAX(zxiflux(:),zepsec)
        ztmp3_1d(:) = MERGE(ztmp3_1d(:), 0._dp, ll2_1d(:))

        zxiflux(:)  = zxiflux(:)  - ztmp2_1d(:)
        zxifluxn(:) = zxifluxn(:) - ztmp3_1d(:)

        ll2_1d(:) = (zxiflux(:) < zepsec)

        zxifluxn(:) = MERGE(0._dp, zxifluxn(:), ll2_1d(:))

        ! Eq. (9) documentation
        zximlt(:)   = ztmst*g*ztmp2_1d(:) / zdp_2d(:,jk)

        ztmp1_1d(:) = pxim1(:,jk) + ztmst * pxite(:,jk)
        ztmp1_1d(:) = MAX(0._dp, ztmp1_1d(:))
        zimlt(:) = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))


        ! If T > T_melt all IC are transfered to cloud droplets
        ztmp1_1d(:) = MERGE(zicncq(:,jk), 0._dp, ll1_1d(:))
        zicnc(:,jk) = MERGE(zicemin, zicnc(:,jk), ll1_1d(:))
        zcdnc(:,jk) = zcdnc(:,jk) + ztmp1_1d(:)
        qmel(1:kproma,jk,jrow) = qmel(1:kproma,jk,jrow) + zdt*ztmp1_1d(:)

        ! --- 3.2 - Submimation of snow (zsub) and ice (zxisub) following -----
        ! ---------------------------------------------------------------------

        ! Based on Lin et al., (J. Clim. Appl. Meteorol., 1983)
        ll1_1d(:) = (zclcpre(:) > 0._dp)  ! fraction of grid covered by precip.

        ztmp1_1d(:) = 1._dp / (2.43e-2_dp * rv) * zlsdcp(:)**2 / ptm1(:,jk)**2
        ztmp1_1d(:) = ztmp1_1d(:) + &
             1._dp / 0.211e-4_dp * zrho_rcp(:,jk) / zqsi_2d(:,jk)
        ztmp1_1d(:) = 3.e6_dp * 2._dp * api * zicesub(:,jk) * zrho_rcp(:,jk) / &
             ztmp1_1d(:)

        ztmp2_1d(:) = MERGE(zclcpre(:), 1._dp, ll1_1d(:))

        ! Snow
        ll2_1d(:) = (zsfl(:) > cqtmin) .AND. ll1_1d(:)

        ztmp3_1d(:) = zcons3 * (zsfl(:) / ztmp2_1d(:))**(0.25_dp / 1.16_dp)
        ztmp3_1d(:) = 0.78_dp * ztmp3_1d(:)**2 + &
             232.19_dp * zqrho_2d(:,jk)**0.25_dp * ztmp3_1d(:)**2.625_dp
        ztmp3_1d(:) = ztmp3_1d(:) * ztmp1_1d(:) * zdpg_2d(:,jk)
        ztmp4_1d(:) = -zxsec * zsfl(:) / ztmp2_1d(:)
        ztmp4_1d(:) = MAX(ztmp4_1d(:), ztmp3_1d(:))
        ztmp4_1d(:) = -ztmst * ztmp4_1d(:) / zdpg_2d(:,jk) * ztmp2_1d(:)
        ztmp5_1d(:) = zxsec * (zqsi_2d(:,jk)-pqm1(:,jk))
        ztmp5_1d(:) = MAX(ztmp5_1d(:), 0._dp)
        ztmp4_1d(:) = MIN(ztmp4_1d(:), ztmp5_1d(:))
        ztmp4_1d(:) = MAX(ztmp4_1d(:), 0._dp)
        zsub(:) = MERGE(ztmp4_1d(:), zsub(:), ll2_1d(:)) ! Eq. (10) docum.

        ! Ice
        ll2_1d(:) = (zxiflux(:) > cqtmin) .AND. ll1_1d(:)

        ztmp3_1d(:) = zcons3 * (zxiflux(:) / ztmp2_1d(:))**(0.25_dp / 1.16_dp)
        ztmp3_1d(:) = 0.78_dp * ztmp3_1d(:)**2 + 232.19_dp * &
             zqrho_2d(:,jk)**0.25_dp * ztmp3_1d(:)**2.625_dp
        ztmp3_1d(:) = ztmp3_1d(:) * ztmp1_1d(:) * zdpg_2d(:,jk)
        ztmp4_1d(:) = -zxsec * zxiflux(:) / ztmp2_1d(:)
        ztmp4_1d(:) = MAX(ztmp4_1d(:), ztmp3_1d(:))
        ztmp4_1d(:) = -ztmst * ztmp4_1d(:) / zdpg_2d(:,jk) * ztmp2_1d(:)
        ztmp5_1d(:) = zxsec * (zqsi_2d(:,jk) - pqm1(:,jk))
        ztmp5_1d(:) = MAX(ztmp5_1d(:), 0._dp)
        ztmp4_1d(:) = MIN(ztmp4_1d(:), ztmp5_1d(:))
        ztmp4_1d(:) = MAX(ztmp4_1d(:), 0._dp)
        zxisub(:) = MERGE(ztmp4_1d(:), zxisub(:), ll2_1d(:))

        ztmp5_1d(:) = zxisub(:) * zxifluxn(:) / MAX(zxiflux(:), cqtmin)
        ztmp5_1d(:) = zcons2 * ztmp5_1d(:) * zdp_2d(:,jk)
        ztmp5_1d(:) = MERGE(ztmp5_1d(:), 0._dp, ll2_1d(:))
        zxifluxn(:) = zxifluxn(:) - ztmp5_1d(:)
        zxiflux(:) = zxiflux(:) - zcons2 * zxisub(:) * zdp_2d(:,jk)

        ll3_1d(:) = ll2_1d(:) .AND. (zxiflux(:) < zepsec)
        zxifluxn(:) = MERGE(0._dp, zxifluxn(:), ll3_1d(:))

        ! --- 3.3 - Evaporation of rain ---------------------------------------
        ! ---------------------------------------------------------------------

        ! Rotstayn (Q. J. Roy. Meteor. Soc., 1997)
        ll2_1d(:) = (zrfl(:) > cqtmin) .AND. ll1_1d(:)

        ztmp3_1d(:) = 870._dp * zsusatw_evap(:,jk) * zdpg_2d(:,jk) * &
             (zrfl(:) / ztmp2_1d(:))**0.61_dp / &
             (SQRT(zrho(:,jk)) * zastbstw(:,jk))

        ztmp4_1d(:) = -zxsec * zrfl(:) / ztmp2_1d(:)
        ztmp4_1d(:) = MAX(ztmp4_1d(:), ztmp3_1d(:))
        ztmp4_1d(:) = -ztmst * ztmp4_1d(:) * ztmp2_1d(:) / zdpg_2d(:,jk)

        ztmp5_1d(:) = zxsec * (zqsw_2d(:,jk) - pqm1(:,jk))
        ztmp5_1d(:) = MAX(ztmp5_1d(:), 0._dp)

        ztmp4_1d(:) = MIN(ztmp4_1d(:), ztmp5_1d(:))
        ztmp4_1d(:) = MAX(ztmp4_1d(:), 0._dp)

        ! Eq. (13) documentation
        zevp(:) = MERGE(ztmp4_1d(:), zevp(:), ll2_1d(:))

        IF (lookupoverflow) THEN
           status_string = 'lookuperror: cdnc - cloud (1)'
           RETURN
        ENDIF

     ENDIF ! jk > 1


     ! === 4 - SEDIMENTATION OF CLOUD ICE FROM GRID-MEAN VALUES ===============
     ! ========================================================================

     ! Updating the tendency 'pxite' to include sedimentation. At jk=klev, the
     ! sedimentation sink is balanced by precipitation at the surface (through
     ! 'zzdrs', see 7.3). Finally: In-cloud cloud water/ice.

     zxip1_1d(:) = pxim1(:,jk) + ztmst * pxite(:,jk) - zimlt(:)
     zxip1_1d(:) = MAX(zxip1_1d(:), EPSILON(1._dp))

     zicncp1_1d(:) = zicnc(:,jk) * paclc(:,jk)
     zicncp1_1d(:) = MAX(zicncp1_1d(:), zicemin)

     zmmean_2d(:,jk) = zrho(:,jk) * zxip1_1d(:) / zicncp1_1d(:)
     zmmean_2d(:,jk) = MAX(zmmean_2d(:,jk), zmi)

     ll1_1d(:) = (zmmean_2d(:,jk) < 2.166E-9_dp )
     ll2_1d(:) = (.not. ll1_1d(:)) .AND. (zmmean_2d(:,jk) < 4.264E-8_dp )

     zalfased_2d(:,jk) = MERGE(63292.4_dp, 8.78_dp, ll1_1d(:))
     zalfased_2d(:,jk) = MERGE(329.75_dp, zalfased_2d(:,jk), ll2_1d(:))

     zbetased_2d(:,jk) = MERGE(0.5727_dp, 0.0954_dp, ll1_1d(:))
     zbetased_2d(:,jk) = MERGE(0.3091_dp, zbetased_2d(:,jk), ll2_1d(:))

     ! Fall velocity
     ! Spichtinger and Gierens (Atmos. Chem. Phys., 2009) - Eq. (18)
     zxifallmc_2d(:,jk) = zfall * zalfased_2d(:,jk) * &
          (zmmean_2d(:,jk)**zbetased_2d(:,jk)) * zaaa_2d(:,jk)

     ! Limit fall velocity to [0.001; 2] m/s
     zxifallmc_2d(:,jk) = MAX(0.001_dp,zxifallmc_2d(:,jk))
     zxifallmc_2d(:,jk) = MIN(2._dp,zxifallmc_2d(:,jk))

     zxifallnc_2d(:,jk) = zxifallmc_2d(:,jk)

     ztmp1_1d(:) = ztmst * g * zxifallmc_2d(:,jk) * zrho(:,jk) / zdp_2d(:,jk)
     ll1_1d(:) = (zxifallmc_2d(:,jk) > zeps)
     ztmp2_1d(:) = zxiflux(:) * zrho_rcp(:,jk) / MAX(zxifallmc_2d(:,jk), zeps)
     ztmp2_1d(:) = MERGE(ztmp2_1d(:), 0._dp, ll1_1d(:))
     ztmp3_1d(:) = g * ztmst * zxifallnc_2d(:,jk) * zrho(:,jk) / zdp_2d(:,jk)

     ll1_1d(:)   = (zxifallnc_2d(:,jk) > zeps)
     ztmp4_1d(:) = zxifluxn(:) / MAX(zxifallnc_2d(:,jk), zeps)
     ztmp4_1d(:) = MERGE(ztmp4_1d(:), 0._dp, ll1_1d(:))

     zxised_2d(:,jk) = zxip1_1d(:) * EXP(-ztmp1_1d(:)) + ztmp2_1d(:) * &
          (1._dp - EXP(-ztmp1_1d(:)))
     zicesed_2d(:,jk) = zicncp1_1d(:) * EXP(-ztmp3_1d(:)) + ztmp4_1d(:) * &
          (1._dp - EXP(-ztmp3_1d(:)))

     ll1_1d(:) = (paclc(:,jk) > 0.01_dp)

     ztmp1_1d(:) = zicesed_2d(:,jk) / MAX(paclc(:,jk), 0.01_dp)
     ztmp1_1d(:) = MERGE(ztmp1_1d(:), zicesed_2d(:,jk), ll1_1d(:))

     zicnc(:,jk) = ztmp1_1d(:) + znidetr(:,jk) + zninucl(:,jk)
     zicnc(:,jk) = MIN(zicnc(:,jk), zicemax)
     zicnc(:,jk) = MAX(zicnc(:,jk), zicemin)

     ! Correct for inconsistencies relative to phases and temperature
     ll1_1d(:) = (ptm1(:,jk) > tmelt) .AND. (zicnc(:,jk) > zicemin)
     zic2cd(:) = MERGE(zicnc(:,jk), 0._dp, ll1_1d(:))
     zicnc(:,jk) = zicnc(:,jk) - zic2cd(:)
     zcdnc(:,jk) = zcdnc(:,jk) + zic2cd(:)
     zicnc(:,jk) = MAX(zicnc(:,jk),zicemin)


     zxiflux(:) = zxiflux(:) + &
          zcons2 * (zxip1_1d(:) - zxised_2d(:,jk)) * zdp_2d(:,jk)
     zxifluxn(:) = zxifluxn(:) + &
          zcons2 * (zicncp1_1d(:) - zicesed_2d(:,jk)) * zdp_2d(:,jk) * &
          zrho_rcp(:,jk)

     pxite(:,jk) = ztmst_rcp * (zxised_2d(:,jk) - pxim1(:,jk))

     zmrateps(:,jk) = zmrateps(:,jk) + zxip1_1d(:) - zxised_2d(:,jk)
     pisedi(:,jk) = zxised_2d(:,jk)  ! save in channel output

     ! In-cloud water/ice calculated from respective grid-means, partial cloud
     ! cover, advective/diffusive tendencies, detrained cloud water/ice and ice
     ! sedimentation. In-cloud values are required for cloud microphysics.

     zclcaux(:) = paclc(:,jk)  ! Cloud cover
     locc_1d(:) = (zclcaux(:) > zeps)

     ! Parameters from Korolev and Mazin (J. Atmos. Sci., 2003)
     ! ztmp1_1d(:) = a2
     ! ztmp2_1d(:) = a0
     ! ztmp3_1d(:) = Ai

     ! Eq. 13-18a Pruppacher and Klett (1997)
     ztmp1_1d(:) = 1._dp / MAX(pqm1(:,jk), zeps) + zlsdcp(:) * &
          alv / (rv * ptm1(:,jk)**2)
     ztmp2_1d(:) = &
          g * (zlvdcp(:) * rd / rv / ptm1(:,jk) - 1._dp) / (rd * ptm1(:,jk))
     zkair_1d(:) = 4.1867e-3_dp * (5.69_dp + 0.017_dp * (ptm1(:,jk) - tmelt))

     zdv_1d(:) = 2.21_dp / papm1(:,jk)
     ztmp3_1d(:) = 1._dp / (crhoi * als**2 / (zkair_1d * rv * ptm1(:,jk)**2) +&
          crhoi * rv * ptm1(:,jk) / (zesi_2d(:,jk) * zdv_1d(:)))

     zrice_1d(:) = (0.75_dp * zxised_2d(:,jk) / &
          (api * zicnc(:,jk) * crhoi))**(1._dp / 3._dp)
     zeta_1d(:) = ztmp1_1d(:) / ztmp2_1d(:) * ztmp3_1d(:) * &
          4._dp * api * crhoi * zcap / zrho(:,jk)
     zvervmax_1d(:) = (zesw_2d(:,jk) - zesi_2d(:,jk)) / zesi_2d(:,jk) * &
          zicnc(:,jk) * zrice_1d(:) * zeta_1d(:)

     lo2_1d(:) = (ptm1(:,jk) < cthomi) .OR. &
          (ptm1(:,jk) < tmelt .AND. 0.01_dp * zvervx_2d(:,jk) < zvervmax_1d(:))

     ll1_1d(:) =  (ptm1(:,jk) < tmelt .AND. ptm1(:,jk) > cthomi .AND. &
          0.01_dp * zvervx_2d(:,jk) < zvervmax_1d(:) .AND. zclcaux(:) > 0._dp)

     IF(ncvmicro) THEN
        zxite(:) = pxteci(:,jk)
        zxlte(:) = pxtecl(:,jk)
     ELSE
        zxite(:) = MERGE(pxtec(:,jk), 0._dp      , lo2_1d(:))
        zxlte(:) = MERGE(0._dp      , pxtec(:,jk), lo2_1d(:))
     ENDIF

     ztmp1_1d(:) = pxim1(:,jk) / MAX(zclcaux(:),zeps)
     ztmp2_1d(:) = pxlm1(:,jk) / MAX(zclcaux(:),zeps)

     zxib(:) = MERGE(ztmp1_1d(:), 0._dp, locc_1d(:))
     zxlb(:) = MERGE(ztmp2_1d(:), 0._dp, locc_1d(:))

     zxim1evp_1d(:) = MERGE(0._dp, pxim1(:,jk), locc_1d(:))
     zxlm1evp_1d(:) = MERGE(0._dp, pxlm1(:,jk), locc_1d(:))

     zxite2(:) = MERGE(zxite(:), 0._dp   , lo2_1d(:))
     zxlte2(:) = MERGE(0._dp   , zxlte(:), lo2_1d(:))

     zxidt_1d(:) = ztmst * (pxite(:,jk) + zxite2(:))
     zxldt_1d(:) = ztmst * (pxlte(:,jk) + zxlte2(:)) + zximlt(:) + zimlt(:)

     ! Ice cloud
     ll1_1d(:) = (zxidt_1d(:) > 0._dp)
     zxidtstar_1d(:) = MERGE(zxidt_1d(:), 0._dp, ll1_1d(:))

     ztmp1_1d(:) = zxidt_1d(:) / MAX(zclcaux(:), zeps)
     ztmp1_1d(:) = MAX(ztmp1_1d(:), -zxib(:))
     ztmp1_1d(:) = MERGE(zxidt_1d(:), ztmp1_1d(:), ll1_1d(:))
     ztmp1_1d(:) = MERGE(ztmp1_1d(:), 0._dp, locc_1d(:))
     zxib(:) = zxib(:) + ztmp1_1d(:)

     ztmp1_1d(:) = - (ztmst_rcp * pxim1(:,jk) + zxite2(:))
     ztmp1_1d(:) = MAX(pxite(:,jk), ztmp1_1d(:))
     pxite(:,jk) = MERGE(pxite(:,jk), ztmp1_1d(:), ll1_1d(:))

     ll2_1d(:) = (.not. locc_1d(:)) .AND. (.not. ll1_1d(:))
     ztmp1_1d(:) = ztmst * ( pxite(:,jk) + zxite2(:) )
     ztmp1_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll2_1d(:))
     zxim1evp_1d(:) = zxim1evp_1d(:) + ztmp1_1d(:)

     ! Water cloud
     ll1_1d(:) = (zxldt_1d(:) > 0._dp)
     zxldtstar_1d(:) = MERGE(zxldt_1d(:), 0._dp, ll1_1d(:))

     ztmp1_1d(:) = zxldt_1d(:) / MAX(zclcaux(:), zeps)
     ztmp1_1d(:) = MAX(ztmp1_1d(:), -zxlb(:))
     ztmp1_1d(:) = MERGE(zxldt_1d(:), ztmp1_1d(:), ll1_1d(:))
     ztmp1_1d(:) = MERGE(ztmp1_1d(:), 0._dp, locc_1d(:))
     zxlb(:) = zxlb(:) + ztmp1_1d(:)

     ztmp1_1d(:) = -(pxlm1(:,jk)/ztmst+zxlte2(:))
     ztmp1_1d(:) = MAX(pxlte(:,jk), ztmp1_1d(:))
     pxlte(:,jk) = MERGE(pxlte(:,jk), ztmp1_1d(:), ll1_1d(:))

     ll2_1d(:) = (.not. locc_1d(:)) .AND. (.not. ll1_1d(:))
     ztmp1_1d(:) = ztmst * ( pxlte(:,jk) + zxlte2(:) )
     ztmp1_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll2_1d(:))
     zxlm1evp_1d(:) = zxlm1evp_1d(:) + ztmp1_1d(:)


     ! === 5 - CONDENSATION/DEPOSITION AND EVAPORATION/SUBLIMATION ============
     ! ========================================================================

     ! zlc =  L_{v/s} / c_p
     ! L latent heat of evaporation (liquid)  or sublimation (ice)

     zlc_1d(:)   = MERGE(zlsdcp(:), zlvdcp(:), lo2_1d(:))
     zqsm1_1d(:) = MERGE(zqsi_2d(:,jk)  , zqsw_2d(:,jk)  , lo2_1d(:))
     zqst1_1d(:) = MERGE(zqsip1_2d(:,jk), zqswp1_2d(:,jk), lo2_1d(:))

     zdqsdt_1d(:) = 1000._dp * (zqst1_1d(:) - zqsm1_1d(:))

     zxievap(:) = (1._dp - zclcaux(:)) * zxidtstar_1d(:) + zxim1evp_1d(:)
     zxlevap(:) = (1._dp - zclcaux(:)) * zxldtstar_1d(:) + zxlm1evp_1d(:)

     zqvdt_1d(:) = ztmst * pqte(:,jk) + zevp(:) + zsub(:) + zxievap(:) + &
          zxlevap(:) + zxisub(:)

     zdtdt_1d(:) = ztmst * ptte(:,jk) -  zlvdcp(:) * (zevp(:) + zxlevap(:)) - &
          (zlsdcp(:) - zlvdcp(:)) * (zsmlt(:) + zximlt(:) + zimlt(:))

     zdtdt_1d(:) = zdtdt_1d(:) - zlsdcp(:) * (zsub(:) + zxievap(:) + zxisub(:))

     zqp1_1d(:) = pqm1(:,jk) + zqvdt_1d(:)
     zqp1_1d(:) = MAX(zqp1_1d(:), 0._dp)

     ztp1_1d(:) = ptm1(:,jk) + zdtdt_1d(:)
     zdqsat_1d(:) = zdtdt_1d(:) + &
          zclcaux(:) * (ztmst * zlc_1d(:) * pqte(:,jk) + &
          zlvdcp(:) * (zevp(:) + zxlevap(:)) + zlsdcp(:) * &
          (zsub(:) + zxievap(:)+zxisub(:)) )

     zdqsat_1d(:) = zdqsat_1d(:) * zdqsdt_1d(:) / (1._dp + zclcaux(:) * &
          zlc_1d(:) * zdqsdt_1d(:))

     zxib(:) = MAX(zxib(:), 0._dp)
     zxlb(:) = MAX(zxlb(:), 0._dp)
     zxilb_1d(:) = zxib(:) + zxlb(:)

     ! Diagnostic output of relative humidity
     prelhum(:,jk) = pqm1(:,jk) / zqsm1_1d(:)
     prelhum(:,jk) = MAX(MIN(prelhum(:,jk), 1._dp), 0._dp)

     ! Extra cloud cover calculation depending on scheme in use
     ! lcover = T/F for Tomkins/Sundqvist
     IF (lcover .AND. jk >= ncctop) THEN

        ! Define variables needed for the cloud cover scheme
        !   zbetaqt = total water
        !   zbetass = saturation mixing ratio adjusted to match qv
        !   zwide   = current diagnosed distribution width
        zbetacl(:) = MAX(0.0_dp, pxlm1(:,jk)) + MAX(0._dp, pxim1(:,jk))
        zbetaqt(:) = MAX(cqtmin, pqm1(:,jk))  + zbetacl(:)
        zvartg(:)  = MAX(cqtmin, cvarmin * pqm1(:,jk))
        zwide(:)   = MAX(zvartg(:), pbetab(:,jk) - pbetaa(:,jk))

        ! --- 5.1 - Turbulence: Skewness - equation solved implicitly ---------
        ! ---------------------------------------------------------------------

        ! This solver only works if phmixtau has non-zero timescale
        zqtau_1d(:) = phmixtau(:,jk) + pvmixtau(:,jk)

        zbqp1_1d(:) = -zdtime * zqtau_1d(:)
        zbqp1_1d(:) = EXP(zbqp1_1d(:))
        zbqp1_1d(:) = cbeta_pq - (cbeta_pq - pxskew(:,jk)) * zbqp1_1d(:)
        zbqp1_1d(:) = MIN(zbqp1_1d(:), cbeta_pq_max)
        zbqp1_1d(:) = MAX(zbqp1_1d(:), cbeta_pq)
        zturbskew(:) = zdtime_rcp*(zbqp1_1d(:)-pxskew(:,jk))

        ! --- 5.2 - Turbulence: Variance - equation solved implicitly ---------
        ! ---------------------------------------------------------------------

        zbbap1_1d(:) = (cbeta_pq + pxskew(:,jk))**2 *   &
             (cbeta_pq+pxskew(:,jk) + 1._dp) / (cbeta_pq * pxskew(:,jk))
        zbbap1_1d(:) = zbbap1_1d(:) * pvdiffp(:,jk) / zwide(:)

        zbbap1_1d(:) = zbbap1_1d(:) / zqtau_1d(:) + zvartg(:) - &
             (zbbap1_1d(:) / zqtau_1d(:) + zvartg(:) - zwide(:)) * &
             EXP(-zdtime * zqtau_1d(:))
        zbbap1_1d(:) = MAX(zbbap1_1d(:), zvartg(:))
        ztmp1_1d(:) = 1._dp / cbeta_pq * zbetaqt(:) * (cbeta_pq + zbqp1_1d(:))
        zbbap1_1d(:) = MIN(zbbap1_1d(:), ztmp1_1d(:))

        zturbvar(:) = zdtime_rcp * (zbbap1_1d(:) - zwide(:))

        zbap1_1d(:) = zbetaqt(:) - &
             cbeta_pq * zbbap1_1d(:) / (cbeta_pq + zbqp1_1d(:))

        ! Translated into apparent xl, xi, q and heat sources first order
        ! effect only, effect of evaporation of cloud on qsat taken into
        ! account in thermodynamic budget but does not change the mixing term
        ! here since that would require iteration and is therefore neglected

        ! Calculate values after one timestep
        iqidx_1d(:) = (zbqp1_1d(:) - cbeta_pq) / rbetak + 1._dp
        iqidx_1d(:) = LOG(iqidx_1d(:))
        iqidx_1d(:) = (nbetaq / cbetaqs) * iqidx_1d(:) + 0.5_dp
        iqidx_1d(:) = INT(iqidx_1d(:))

        ztmp1_1d(:) = cbeta_pq * (pbetass(:,jk) - zbap1_1d(:)) / &
             ((zbetaqt(:) - zbap1_1d(:)) * (cbeta_pq + zbqp1_1d(:)))
        ztmp1_1d(:) = MAX(MIN(ztmp1_1d(:), 1._dp), 0._dp)
        ztmp1_1d(:) = nbetax * ztmp1_1d(:)
        ixidx_1d(:) = INT(ztmp1_1d(:))

        ll1_1d(:) = (ixidx_1d(:) == nbetax)

        ixidx_1d(:) = MERGE(nbetax - 1._dp, ixidx_1d(:), ll1_1d(:))

        DO jl=1,kproma
           iqidx = iqidx_1d(jl)
           ixidx = ixidx_1d(jl)

           ztmp2_1d(jl) = &
                (ztmp1_1d(jl) - ixidx_1d(jl)) * tbetai0(iqidx, ixidx+1) + &
                (ixidx_1d(jl) + 1._dp - ztmp1_1d(jl)) * tbetai0(iqidx, ixidx)
           ztmp3_1d(jl) = &
                (ztmp1_1d(jl) - ixidx_1d(jl)) * tbetai1(iqidx,ixidx+1) + &
                (ixidx_1d(jl) + 1._dp - ztmp1_1d(jl)) * tbetai1(iqidx, ixidx)
        ENDDO

        ztmp4_1d(:) = MERGE(1._dp, ztmp2_1d(:), ll1_1d(:))
        ztmp5_1d(:) = MERGE(1._dp, ztmp3_1d(:), ll1_1d(:))

        ztmp1_1d(:) = -zxilb_1d(:) * zclcaux(:)
        ztmp2_1d(:) = zqsec*zqp1_1d(:)

        zgent_1d(:) = pqm1(:,jk) - (zbetaqt(:) - zbap1_1d(:)) * ztmp5_1d(:) + &
             (pbetass(:,jk) - zbap1_1d(:)) * ztmp4_1d(:) - pbetass(:,jk)

        zgent_1d(:) = MAX(zgent_1d(:), ztmp1_1d(:))
        zgent_1d(:) = MIN(zgent_1d(:), ztmp2_1d(:))

        ztmp1_1d(:) = MAX(zepsec, zxilb_1d(:))
        ztmp1_1d(:) = zxib(:) / ztmp1_1d(:)
        ztmp1_1d(:) = MAX(MIN(ztmp1_1d(:), 1._dp), 0._dp)
        ztmp2_1d(:) = 1._dp - ztmp1_1d(:)

        ztmp3_1d(:) = zgent_1d(:) / MAX(zclcaux(:), zeps)

        zgenti(:) = zgent_1d(:) * ztmp1_1d(:)
        zgentl(:) = zgent_1d(:) * ztmp2_1d(:)

        ztmp4_1d(:) = zxib(:) + ztmp3_1d(:) * ztmp1_1d(:)
        ztmp4_1d(:) = MAX(ztmp4_1d(:), 0._dp)
        zxib(:) = MERGE(ztmp4_1d(:), zxib(:), locc_1d(:))

        ztmp4_1d(:) = zxlb(:) + ztmp3_1d(:) * ztmp2_1d(:)
        ztmp4_1d(:) = MAX(ztmp4_1d(:), 0.0_dp)
        zxlb(:) = MERGE(ztmp4_1d(:), zxlb(:), locc_1d(:))

        zxilb_1d(:) = zxib(:) + zxlb(:)

        ! --- 5.3 - Deposition/sublimation of cloud ice and condensation/
        ! --- evaporation of liquid water due to changes in water vapour and
        ! --- temperature (advection, convection, turbulent mixing,
        ! --- evaporation, of rain, sublimation and melting of snow)
        ! ---------------------------------------------------------------------

        ! Translate PDF laterally to calculate cloud after one timestep
        zqvdt_1d(:) = zqvdt_1d(:) - zgent_1d(:)
        zdtdt_1d(:) = &
             zdtdt_1d(:) + zlvdcp(:) * zgentl(:) + zlsdcp(:) * zgenti(:)

        zqp1_1d(:) = pqm1(:,jk) + zqvdt_1d(:)
        zqp1_1d(:) = MAX(zqp1_1d(:), 0._dp)

        ztp1_1d(:) = ptm1(:,jk) + zdtdt_1d(:)

        zdqsat_1d(:) = zdtdt_1d(:) + &
             zclcaux(:) * (ztmst * zlc_1d(:) * pqte(:,jk) + &
             zlvdcp(:) * (zevp(:) + zxlevap(:) - zgentl(:)) + &
             zlsdcp(:) * (zsub(:) + zxievap(:) + zxisub(:) - zgenti(:)))

        zdqsat_1d(:) = zdqsat_1d(:) * zdqsdt_1d(:) / (1._dp + zclcaux(:) * &
             zlc_1d(:) * zdqsdt_1d(:))

        ztmp1_1d(:) = (pbetass(:,jk) - zqvdt_1d(:) + zdqsat_1d(:) - &
             zbap1_1d(:)) / zbbap1_1d(:)
        ztmp1_1d(:) = MAX(MIN(ztmp1_1d(:), 1._dp), 0._dp)
        ztmp1_1d(:) = nbetax * ztmp1_1d(:)
        ixidx_1d(:) = INT(ztmp1_1d(:))

        ll1_1d(:) = (ixidx_1d(:) == nbetax)

        ixidx_1d(:) = MERGE(nbetax - 1._dp, ixidx_1d(:), ll1_1d(:))

        DO jl=1,kproma

           iqidx = iqidx_1d(jl)
           ixidx = ixidx_1d(jl)

           ztmp2_1d(jl) = &
                (ztmp1_1d(jl)-ixidx_1d(jl)) * tbetai0(iqidx,ixidx+1) + &
                (ixidx_1d(jl) + 1._dp - ztmp1_1d(jl)) * tbetai0(iqidx, ixidx)
           ztmp3_1d(jl) = &
                (ztmp1_1d(jl)-ixidx_1d(jl)) * tbetai1(iqidx, ixidx+1) + &
                (ixidx_1d(jl) + 1._dp - ztmp1_1d(jl)) * tbetai1(iqidx, ixidx)
        ENDDO

        ztmp4_1d(:) = MERGE(1._dp, ztmp2_1d(:), ll1_1d(:))
        ztmp5_1d(:) = MERGE(1._dp, ztmp3_1d(:), ll1_1d(:))

        zqcdif_1d(:) = (zbetaqt(:) - pbetaa(:,jk)) * (1._dp - ztmp5_1d(:)) + &
             (pbetaa(:,jk) + zqvdt_1d(:) - pbetass(:,jk) - zdqsat_1d(:)) * &
             (1._dp-ztmp4_1d(:))

        zqcdif_1d(:) = MAX(0._dp, zqcdif_1d(:))
        zqcdif_1d(:) = zqcdif_1d(:)-zbetacl(:)

        ztmp1_1d(:) = -zxilb_1d(:) * zclcaux(:)
        ztmp2_1d(:) = zqsec * zqp1_1d(:)

        zqcdif_1d(:) = MAX(zqcdif_1d(:), ztmp1_1d(:))
        zqcdif_1d(:) = MIN(zqcdif_1d(:), ztmp2_1d(:))

     ELSE ! lcover = F or jk < ncctop

        zqcdif_1d(:) = (zqvdt_1d(:) - zdqsat_1d(:)) * zclcaux(:)

        ztmp1_1d(:) = -zxilb_1d(:) * zclcaux(:)
        ztmp2_1d(:) = zqsec * zqp1_1d(:)

        zqcdif_1d(:) = MAX(zqcdif_1d(:), ztmp1_1d(:))
        zqcdif_1d(:) = MIN(zqcdif_1d(:), ztmp2_1d(:))  ! limit to qv

     ENDIF ! lcover

     ll1_1d(:) = (zqcdif_1d(:) < 0._dp)  ! cloud dissipation

     ztmp1_1d(:) = MAX(zepsec, zxilb_1d(:))
     ztmp1_1d(:) = zxib(:) / ztmp1_1d(:)
     ztmp1_1d(:) = MAX(MIN(ztmp1_1d(:), 1._dp), 0._dp)

     ztmp1_1d(:) = MERGE(ztmp1_1d(:), 1._dp, ll1_1d(:))
     ztmp2_1d(:) = zqcdif_1d(:)*(1.0_dp-ztmp1_1d(:))
     zcnd(:) = MERGE(ztmp2_1d(:), 0._dp, ll1_1d(:))

     IF (nicnc <= 1) THEN
        zdep(:) = zqcdif_1d(:) * ztmp1_1d(:)
     ELSE
        zdep(:) = zqinucl(:,jk) * ztmp1_1d(:)
     ENDIF

     ! No cloud dissipation and condensation
     ll2_1d(:) = (.not. ll1_1d(:)) .AND. &
          (.not. lo2_1d(:))
     zdep(:) = MERGE(0._dp, zdep(:), ll2_1d(:))

     ! Use standard condensation for empirical Lin & Leaitch approach and
     ! explicit condensation after Levkov et al. (1992) for the explicit
     ! activation schemes that allow for supersaturation:

     SELECT CASE(ncdnc)
     CASE(1) ! Standard condensation for Lin and Leaitch activation scheme

        zcnd(:) = MERGE(zqcdif_1d(:), zcnd(:), ll2_1d(:))

     CASE(2) ! Explicit condensation after Levkov et al. (1992) for ARG
             ! activation scheme that allow for supersaturation
        ztmp1_1d(:) = 0.5_dp * api * zdw0 * ztmst * &
             swat(1:kproma,jk,jrow) * zcdnc(:,jk) * zclcaux(:) * &
             zrho_rcp(:,jk) / zastbstw(:,jk)
        ztmp1_1d(:) = MIN(ztmp1_1d(:), zqcdif_1d(:))
        zcnd(:) = MERGE(ztmp1_1d(:), zcnd(:), ll2_1d(:))

     CASE DEFAULT

        WRITE(*,*) "ncdnc can be either 1 or 2 in messy_cloud_Kuebbeler14"
        STOP

     END SELECT

     ! --- 5.4 - Account for cloud evaporation in clear air and check for
     ! --- supersaturation ----------------------------------------------------
     ! ------------------------------------------------------------------------

     ztp1tmp(:) = ztp1_1d(:) + zlvdcp(:) * zcnd(:) + zlsdcp(:) * zdep(:)
     zqp1tmp(:) = zqp1_1d(:) - zcnd(:) - zdep(:)

     zxip1_1d(:) = &
          zxised_2d(:,jk) + ztmst*zxite(:) - zxievap(:) + zgenti(:) + zdep(:)
     zxip1_1d(:) = MAX(zxip1_1d(:), 0._dp)

     lo2_1d(:) = (ztp1tmp(:) < cthomi) .OR. &
          (ztp1tmp(:) < tmelt .AND. zxip1_1d(:) > csecfrl  &
          .AND. zsusatw_2d(:,jk) < zeps)

     ztmp1_1d(:) = 1000._dp * ztp1tmp(:)
     it_1d(:) = NINT(ztmp1_1d(:))

     ll_look(:,jk) = (it_1d(:) < jptlucu1 .OR. it_1d(:) > jptlucu2)
     IF (ANY(ll_look(1:kproma,jk))) lookupoverflow = .TRUE.

     it_1d(:) = MAX(MIN(it_1d(:),jptlucu2),jptlucu1)

     DO jl=1,kproma
        zlucua_1d(jl)  = tlucua(it_1d(jl))
        zlucuaw_1d(jl) = tlucuaw(it_1d(jl))

        zlucuap1_1d(jl)  = tlucua(it_1d(jl)+1)
        zlucuawp1_1d(jl) = tlucuaw(it_1d(jl)+1)

        zlucub_1d(jl) = tlucub(it_1d(jl))
     ENDDO

     zes_1d(:) = MERGE(zlucua_1d(:), zlucuaw_1d(:), lo2_1d(:))
     zes_1d(:) = zes_1d(:) / papp1(:,jk)
     zes_1d(:) = MIN(zes_1d(:), 0.5_dp)

     ll1_1d(:) = (zes_1d(:) < 0.4_dp)

     zcor_1d(:) = 1._dp/(1._dp-vtmpc1*zes_1d(:))

     zqsp1tmp_1d(:) = zes_1d(:)*zcor_1d(:)
     zoversat_1d(:) = zqsp1tmp_1d(:)*0.01_dp

     zrhtest_1d(:) = pqm1(:,jk)/zqsm1_1d(:)
     zrhtest_1d(:) = MIN(zrhtest_1d(:), 1._dp)
     zrhtest_1d(:) = zrhtest_1d(:)*zqsp1tmp_1d(:)

     ! For cirrus scheme
     zesw_1d(:) = zlucuaw_1d(:) / papp1(:,jk)
     zesw_1d(:) = MIN(zesw_1d(:), 0.5_dp)
     zcorw_1d(:) = 1._dp / (1._dp - vtmpc1 * zesw_1d(:))
     zqsp1tmpw_1d(:) = zesw_1d(:) * zcorw_1d(:)
     zoversatw_1d(:) = 0.01_dp * zqsp1tmpw_1d(:)

     zqst1_1d(:) = MERGE(zlucuap1_1d(:), zlucuawp1_1d(:), lo2_1d(:))
     zqst1_1d(:) = zqst1_1d(:) / papp1(:,jk)
     zqst1_1d(:) = MIN(zqst1_1d(:), 0.5_dp)
     zqst1_1d(:) = zqst1_1d(:) / (1._dp - vtmpc1 * zqst1_1d(:))

     zdqsdt_1d(:) = 1000._dp * (zqst1_1d(:) - zqsp1tmp_1d(:))

     zlc_1d(:) = MERGE(zlsdcp(:), zlvdcp(:), lo2_1d(:))

     ztmp1_1d(:) = zlc_1d(:) * zdqsdt_1d(:)
     ztmp2_1d(:) = zqsp1tmp_1d(:) * zcor_1d(:) * zlucub_1d(:)
     zlcdqsdt_1d(:) = MERGE(ztmp1_1d(:), ztmp2_1d(:), ll1_1d(:))

     zqcon_1d(:) = 1._dp / (1._dp + zlcdqsdt_1d(:))

     ztmp1_1d(:) = zqsp1tmpw_1d(:) + zoversatw_1d(:)
     ztmp1_1d(:) = MIN(ztmp1_1d(:), zqsp1tmp_1d(:) * zsicirrus_2d(:,jk))
     zqsp1tmpcirrus_1d(:) = MERGE(ztmp1_1d(:), 0._dp, lo2_1d(:))

     ! Reduce possible supersaturation above the homogeneous freezing threshold
     ll1_1d(:) = (nicnc <= 1)

     ! For the Kaercher et al. (2006) parametrization (nicnc = 2) which
     ! considers several competing freezing mechanisms, this means reducing the
     ! supersaturation above the highest supersaturation reached in the
     ! parametrization

     ! S_ice above ice supersaturation
     ll2_1d(:) = (zqp1tmp(:) > (zqsp1tmp_1d(:)  + zoversat_1d(:) ) )

     ! S_ice above highest S_ice in Kaercher parametrization
     ll3_1d(:) = (zqp1tmp(:) > zqsp1tmpcirrus_1d(:))

     ! Above homogeneous freezing threshold
     ll5_1d(:) = (ztp1tmp(:) >= cthomi)

     ! Ice mass above ice supersaturation
     ztmp1_1d(:) = (zqp1tmp(:) - zqsp1tmp_1d(:) - zoversat_1d(:) ) * zqcon_1d(:)

     ! Ice mass above highest S_ice in Kaercher parametrization
     ztmp3_1d(:) = (zqp1tmp(:) - zqsp1tmpcirrus_1d(:)) * zqcon_1d(:)

     ! Ice cloud cases

     ! nicnc = 1
     ll6_1d(:) = (lo2_1d(:) .AND. ll1_1d(:) .AND. ll2_1d(:)) .OR. &
          (lo2_1d(:) .AND. (.not. ll1_1d(:)) .AND. ll2_1d(:) .AND. ll5_1d(:))
     ztmp4_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll6_1d(:))

     ! nicnc = 2
     ll6_1d(:) = lo2_1d(:) .AND. (.not. ll1_1d(:)) .AND. ll3_1d(:)    &
          .AND. (.not. ll5_1d(:))
     ztmp4_1d(:) = MERGE(ztmp3_1d(:), ztmp4_1d(:), ll6_1d(:))

     zdep(:) = zdep(:) + ztmp4_1d(:)

     ! Water cloud cases
     ll6_1d(:) = (.not. lo2_1d(:)) .AND. ll2_1d(:)
     ztmp4_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll6_1d(:))
     zcnd(:) = zcnd(:) + ztmp4_1d(:)

     ! Final corrections
     ztmp5_1d(:) = zqp1_1d(:) - zrhtest_1d(:)
     ztmp5_1d(:) = MAX(ztmp5_1d(:), 0._dp)

     ll1_1d(:) = (zdep(:) > 0._dp)
     ll2_1d(:) = (zcnd(:) > 0._dp)
     ll3_1d(:) = (zqp1tmp(:) < zrhtest_1d(:))
     ll4_1d(:) = (zqsp1tmp_1d(:) <= zqsm1_1d(:))
     ll5_1d(:) = lo2_1d(:) .AND. ll1_1d(:) .AND. ll3_1d(:) .AND. ll4_1d(:)

     zdep(:) = MERGE(ztmp5_1d(:), zdep(:), ll5_1d(:))

     ll5_1d(:) = &
          (.not.lo2_1d(:)) .AND. ll2_1d(:) .AND. ll3_1d(:) .AND. ll4_1d(:)
     zcnd(:) = MERGE(ztmp5_1d(:), zcnd(:), ll5_1d(:))

     ! --- 5.5 - Change of in-cloud water due to deposition/sublimation and
     ! --- condensation/evaporation (input for cloud microphysics) ------------
     ! ------------------------------------------------------------------------

     zrelhum_1d(:) = zqp1tmp(:) / zqsp1tmp_1d(:)

     ztmp1_1d(:) = zdep(:) + zgenti(:)
     ztmp1_1d(:) = MAX(ztmp1_1d(:), 0._dp)

     ztmp2_1d(:) = zcnd(:) + zgentl(:)
     ztmp2_1d(:) = MAX(ztmp2_1d(:), 0._dp)

     ztmp3_1d(:) = zxib(:) + zdep(:) / MAX(zclcaux(:), zeps)
     ztmp3_1d(:) = MAX(ztmp3_1d(:), 0._dp)

     ztmp4_1d(:) = zxlb(:) + zcnd(:) / MAX(zclcaux(:), zeps)
     ztmp4_1d(:) = MAX(ztmp4_1d(:), 0._dp)

     zxib(:) = MERGE(ztmp3_1d(:), zxib(:), locc_1d(:))
     zxlb(:) = MERGE(ztmp4_1d(:), zxlb(:), locc_1d(:))

     ll1_1d(:) = (.not. locc_1d(:)) &
          .AND. ((ztmp1_1d(:) > 0._dp) .OR. (ztmp2_1d(:) > 0._dp))

     ztmp3_1d(:) = MAX(MIN(zrelhum_1d(:), 1._dp), 0.01_dp)
     zclcaux(:)  = MERGE(ztmp3_1d(:), zclcaux(:), ll1_1d(:))

     ztmp3_1d(:) = ztmp1_1d(:) / MAX(zclcaux(:), zeps)
     ztmp4_1d(:) = ztmp2_1d(:) / MAX(zclcaux(:), zeps)

     zxib(:) = MERGE(ztmp3_1d(:), zxib(:), ll1_1d(:))
     zxlb(:) = MERGE(ztmp4_1d(:), zxlb(:), ll1_1d(:))

     ! Redefine locc_1d which was no longer true:
     locc_1d(:) = (zclcaux(:) > 0._dp)

     ! Update cdnc:
     ll1_1d(:) = locc_1d(:) .AND. (zxlb(:) > cqtmin)

     ! No previous nucleation
     ll2_1d(:) = ll1_1d(:) .AND. (zcdnc(:,jk) <= cdncmin) .AND. &
          (ptm1(:,jk) > cthomi)

     ztmp1_1d(:) = pcdncact(:,jk) - zcdnc(:,jk)
     ztmp1_1d(:) = MAX(0._dp, ztmp1_1d(:))
     zqlnuc(:,jk) = MERGE(ztmp1_1d(:), zqlnuc(:,jk), ll2_1d(:))

     ztmp1_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll2_1d(:))
     zcdnc(:,jk) = zcdnc(:,jk) + ztmp1_1d(:)
     qnuc(1:kproma,jk,jrow) = qnuc(1:kproma,jk,jrow) + zdt*ztmp1_1d(:)

     ztmp1_1d(:) = MAX(zcdnc(:,jk), cdncmin)
     zcdnc(:,jk) = MERGE(ztmp1_1d(:), cqtmin, ll1_1d(:))

     ! Update icnc:
     ll1_1d(:) = locc_1d(:) .AND. (zxib(:) > cqtmin)
     ll2_1d(:) = ll1_1d(:)  .AND. (zicnc(:,jk) <= zicemin)

     IF (nicnc <= 1) THEN
        ztmp1_1d(:) = 0.75_dp / &
             (api * zrhoice) * zrho(:,jk) * zxib(:) / zrid_2d(:,jk)**3
     ELSE
        ztmp1_1d(:) = MIN(znicex(:,jk), (zap(:,jk)*1.e6_dp))
     ENDIF

     zicnc(:,jk) = MERGE(ztmp1_1d(:), zicnc(:,jk), ll2_1d(:))

     ztmp1_1d(:) = MAX(zicnc(:,jk), zicemin)
     zicnc(:,jk) = MERGE(ztmp1_1d(:), cqtmin, ll1_1d(:))

     ztp1tmp(:) = ztp1_1d(:) + zlvdcp(:) * zcnd(:) + zlsdcp(:) * zdep(:)
     zqp1tmp(:) = zqp1_1d(:) - zcnd(:) - zdep(:)


     ! === 6 - FREEZING OF CLOUD WATER ========================================
     ! ========================================================================

     ! --- 6.1 - Freezing of all cloud water for T < 238 K --------------------
     ! ------------------------------------------------------------------------

     ll1_1d(:) = (ztp1tmp(:) <= cthomi)

     ztmp1_1d(:) = zfrl(:,jk) + zxlb(:) * zclcaux(:)
     zfrl(:,jk)  = MERGE(ztmp1_1d(:), zfrl(:,jk), ll1_1d(:))

     ztmp1_1d(:) = zxib(:)+zxlb(:)
     zxib(:) = MERGE(ztmp1_1d(:), zxib(:), ll1_1d(:))
     zxlb(:) = MERGE(0._dp, zxlb(:), ll1_1d(:))

     ztmp1_1d(:) = zcdnc(:,jk) - cdncmin
     ztmp1_1d(:) = MAX(ztmp1_1d(:), 0._dp)
     ztmp2_1d(:) = qfre(1:kproma,jk,jrow) - zdt * ztmp1_1d(:)
     qfre(1:kproma,jk,jrow) = &
          MERGE(ztmp2_1d(:), qfre(1:kproma,jk,jrow), ll1_1d(:))

     ztmp2_1d(:) = zicnc(:,jk) + ztmp1_1d(:)
     zicnc(:,jk) = MERGE(ztmp2_1d(:), zicnc(:,jk), ll1_1d(:))
     zcdnc(:,jk) = MERGE(cqtmin, zcdnc(:,jk), ll1_1d(:))

     ! --- 6.2 - Freezing of cloud water between 238 and 273 K ----------------
     ! ------------------------------------------------------------------------

     lo_1d(:) = (zxlb(:) > cqtmin)       .AND. &
                (ztp1tmp(:)  < tmelt)    .AND. &
                (ztp1tmp(:)  > cthomi  ) .AND. &
                (zcdnc(:,jk) >= cdncmin) .AND. &
                locc_1d(:)

     ! Mixed-phase clouds parametrization by
     ! Lohmann and Diehl (J. Atmos. Sci., 2006)

     ! Cunningham factors in eaech size mode
     ztmp1_1d(:) = 1._dp + 1.26_dp * 6.6E-8_dp / (zrwetki_2d(:,jk) + zeps) * &
          (101325._dp / papp1(:,jk)) * (ztp1tmp(:) / tmelt)
     ztmp2_1d(:) = 1._dp + 1.26_dp * 6.6E-8_dp / (zrwetai_2d(:,jk) + zeps) * &
          (101325._dp/papp1(:,jk)) * (ztp1tmp(:) / tmelt)
     ztmp3_1d(:) = 1._dp + 1.26_dp * 6.6E-8_dp / (zrwetci_2d(:,jk) + zeps) * &
          (101325._dp/papp1(:,jk)) * (ztp1tmp(:) / tmelt)

     ! Air viscosity
     zetaair_1d(:) = 1.e-5_dp * (1.718_dp + 0.0049_dp * (ztp1tmp(:) - tmelt) -&
          1.2e-5_dp * (ztp1tmp(:) - tmelt) * (ztp1tmp(:) - tmelt))

     ! For use in thermophoresis below (Eq. 13-18a, Pruppacher and Klett, 1997)
     zkair_1d(:) = 4.1867e-3_dp * (5.69_dp + 0.017_dp * (ztp1tmp(:) - tmelt))

     ! Brownian diffusions in each size mode
     ll1_1d(:)= (zrwetki_2d(:,jk) < zeps)
     zdfarbcki_1d(:) = ak * ztp1tmp(:) * ztmp1_1d(:) / &
          (6._dp * api * zetaair_1d(:) * (zrwetki_2d(:,jk) + zeps))
     zdfarbcki_1d(:) = MERGE(0._dp, zdfarbcki_1d(:), ll1_1d(:))

     ll2_1d(:) = (zrwetai_2d(:,jk) < zeps)
     zdfarduai_1d(:) = ak * ztp1tmp(:) * ztmp2_1d(:) / &
          (6._dp * api * zetaair_1d(:) * (zrwetai_2d(:,jk) + zeps))
     zdfarduai_1d(:) = MERGE(0._dp, zdfarduai_1d(:), ll2_1d(:))

     ll3_1d(:) = (zrwetci_2d(:,jk) < zeps)
     zdfarduci_1d(:) = ak * ztp1tmp(:) * ztmp3_1d(:) / &
          (6._dp * api * zetaair_1d(:) * (zrwetci_2d(:,jk) + zeps))
     zdfarduci_1d(:) = MERGE(0._dp, zdfarduci_1d(:), ll3_1d(:))

     ! Thermophoresis
     IF (lthermo) THEN
        zknbcki_1d(:) = 7.37_dp * ztp1tmp(:) / &
             (2.88e5_dp * papp1(:,jk) * (zrwetki_2d(:,jk) + zeps))
        zknbcki_1d(:) = MERGE(1._dp, zknbcki_1d(:), ll1_1d(:))
        zftbcki_1d(:) = 0.4_dp * (1._dp + zknbcki_1d(:) * &
             (1.45_dp + 0.4_dp * EXP(-1._dp / zknbcki_1d(:)))) * &
             (zkair_1d(:) + 2.5_dp * zknbcki_1d(:) * zkbc) / &
             ((1._dp + 3._dp * zknbcki_1d(:)) * &
             (2._dp * zkair_1d(:) + zkbc * (5._dp * zknbcki_1d(:) + 1._dp)))
        zftbcki_1d(:) = MERGE(0._dp, zftbcki_1d(:), ll1_1d(:))

        zknduai_1d(:) = 7.37_dp * ztp1tmp(:) / &
             (2.88e5_dp * papp1(:,jk) * (zrwetai_2d(:,jk) + zeps))
        zknduai_1d(:) = MERGE(1._dp, zknduai_1d(:), ll2_1d(:))
        zftduai_1d(:) = 0.4_dp * (1._dp + zknduai_1d(:) * &
             (1.45_dp + 0.4_dp * EXP(-1._dp / zknduai_1d(:)))) * &
             (zkair_1d(:) + 2.5_dp * zknduai_1d(:) * zkdu) / &
             ((1._dp + 3._dp * zknduai_1d(:)) * &
             (2._dp * zkair_1d(:) + zkdu * (5._dp * zknduai_1d(:) + 1._dp)))
        zftduai_1d(:) = MERGE(0._dp, zftduai_1d(:), ll2_1d(:))

        zknduci_1d(:) = 7.37_dp * ztp1tmp(:) / &
             (2.88e5_dp * papp1(:,jk) * (zrwetci_2d(:,jk) + zeps))
        zknduci_1d(:) = MERGE(1._dp, zknduci_1d(:), ll3_1d(:))
        zftduci_1d(:) = 0.4_dp * (1._dp + zknduci_1d(:) * &
             (1.45_dp + 0.4_dp * EXP(-1._dp / zknduci_1d(:)))) * &
             (zkair_1d(:) + 2.5_dp * zknduci_1d(:) * zkdu) / &
             ( (1._dp+3._dp * zknduci_1d(:)) * &
             (2._dp * zkair_1d(:) + zkdu * (5._dp * zknduci_1d(:) + 1._dp)))
        zftduci_1d(:) = MERGE(0._dp, zftduci_1d(:), ll3_1d(:))
     ENDIF


     DO jl =1, kproma

        zfracdusol = &  ! Dust immersion freezing
             MIN(ndusol_strat(jl,jk,jrow) / (pcdncact(jl,jk) + zeps), 1._dp)
        zfracduinsolai = &  ! Dust (accumulation mode )contact freezing
             MIN(nduinsolai(jl,jk,jrow)  /(naerinsol(jl,jk,jrow) + zeps), 1._dp)
        zfracduinsolci = &  ! Dust (coarse mode) contact freezing
             MIN(nduinsolci(jl,jk,jrow)  /(naerinsol(jl,jk,jrow) + zeps), 1._dp)
        zfracbcsol = &  ! BC immersion freezing
             MIN(nbcsol_strat(jl,jk,jrow) / (pcdncact(jl,jk) + zeps), 1._dp)
        zfracbcinsol = &  ! BC contact freezing
             MIN(nbcinsol(jl,jk,jrow) + nbctaginsol(jl,jk,jrow) / &
             (naerinsol(jl,jk,jrow) + zeps), 1._dp)

        ! Droplet radius [m]
        zradl = (0.75_dp * zxlb(jl) * zrho(jl,jk) / &
             (api * rhoh2o * zcdnc(jl,jk)))**(1._dp/3._dp)

        zf1 = 4._dp*api * zradl * zcdnc(jl,jk) * zrho_rcp(jl,jk)

        ! Dust
        zfrzcntdu = &  ! Montmorillonite
             MIN(1._dp, &
             MAX(0._dp, -(0.1014_dp * (ztp1tmp(jl) - tmelt) + 0.3277_dp)))

        ! Black carbon
        zfrzcntbc = 0._dp ! disable BC contact freezing

        zfrzcnt = &
             zxlb(jl) / zcdnc(jl,jk) * zrho(jl,jk) * zf1 * (zfrzcntdu * &
             (zdfarduai_1d(jl) * zfracduinsolai +  &
             zdfarduci_1d(jl) * zfracduinsolci) + &
             zfrzcntbc * zdfarbcki_1d(jl) * &
             zfracbcinsol) * (zcdnc(jl,jk) + zicnc(jl,jk))

        zfrzcnt = &
             zxlb(jl) * (1._dp - EXP(-zfrzcnt / MAX(zxlb(jl), cqtmin) * ztmst))


        ! Add thermophoresis
        zfrzthermo = 0._dp
        IF (lthermo) THEN

           lo2_1d(jl) = (ztp1tmp(jl) < cthomi) .OR.                     &
                (ztp1tmp(jl) < tmelt .AND. zxip1_1d(jl) > csecfrl .AND. &
                zsusatw_2d(jl,jk) < zeps)

           zqsm1 = MERGE(tlucua(itm1_look(jl,jk)), &
                tlucuaw(itm1_look(jl,jk)), lo2_1d(jl)) / papm1(jl,jk)
           zqsm1 = MIN(zqsm1, 0.5_dp)
           zqsm1 = zqsm1 / (1._dp - vtmpc1 * zqsm1)

           it1_1d(jl) = NINT(ztp1tmp(jl) * 1000._dp)
           it1_1d(jl) = MAX(MIN(it1_1d(jl), jptlucu2), jptlucu1)
           IF (it1_1d(jl) < jptlucu1 .OR. &
                it1_1d(jl) > jptlucu2) lookupoverflow = .TRUE.

           zqst1 = &
                MERGE(tlucua(it1_1d(jl)), tlucuaw(it1_1d(jl)), lo2_1d(jl)) / &
                      papm1(jl,jk)
           zqst1 = MIN(zqst1, 0.5_dp)
           zqst1 = zqst1 / (1._dp - vtmpc1 * zqst1)
           zdqsdtime = (zqst1-zqsm1)
           zdeltatemp = MAX((zrho(jl,jk) * alv * (zdep(jl) + zdqsdtime)) / &
                (ztmst * zkair_1d(jl) * &
                4._dp * api * zcdnc(jl,jk) * zradl + zeps), 0._dp)
           zf2 = zkair_1d(jl) / papp1(jl,jk) * zdeltatemp
           zfrzthermo = zf1 * zf2 * &
                (zfrzcntbc * zftbcki_1d(jl) * zfracbcinsol + &
                zfrzcntdu * (zftduai_1d(jl) * zfracduinsolai + &
                zftduci_1d(jl) * zfracduinsolci)) * &
                (zcdnc(jl,jk) + zicnc(jl,jk)) * ztmst / &
                (zcdnc(jl,jk) * zrho(jl,jk))

           zfrzthermo = zxlb(jl) * (1._dp - EXP(-zfrzthermo))

        ENDIF

        znaimmdu = 32.3_dp * zfracdusol  ! montmorillonite
        ! znaimmdu = 6.15E-2_dp * zfracdusol  ! kaolinite
        znaimmbc  = 2.91E-3_dp * zfracbcsol  ! black carbon

        ! Vertical velocity (large-scale + 1.33 SQRT(TKE)) [Pa/s]
        zomega = &
             pvervel(jl,jk) - 1.33_dp * SQRT(ptkem1(jl,jk)) * zrho(jl,jk) * g

        ! Temperature tendency [K/s]
        ztte   = zomega / cpd * zrho_rcp(jl,jk)

        zfrzimm = -(znaimmdu + znaimmbc) * zrho(jl,jk) / rhoh2o * &
             EXP(tmelt - ztp1tmp(jl)) * MIN(ztte, 0._dp)

        zfrzimm = zxlb(jl) * &
             (1._dp - EXP(-zfrzimm * zxlb(jl) / zcdnc(jl,jk) * ztmst))

        ztmp1_1d(jl) = zfrzcnt + zfrzimm + zfrzthermo
        ztmp1_1d(jl) = MAX(0._dp, MIN(ztmp1_1d(jl), zxlb(jl)))

        ztmp2_1d(jl) = zcdnc(jl,jk) * ztmp1_1d(jl) / (zxlb(jl) + zeps)
        ztmp3_1d(jl) = ndusol_strat(jl,jk,jrow) + nbcsol_strat(jl,jk,jrow) + &
             nduinsolai(jl,jk,jrow) + nduinsolci(jl,jk,jrow) + &
             nbcinsol(jl,jk,jrow) + nbctaginsol(jl,jk,jrow) - zicnc(jl,jk)

        ztmp2_1d(jl) = MIN(ztmp2_1d(jl), ztmp3_1d(jl))
        ztmp2_1d(jl) = MAX(ztmp2_1d(jl), 0._dp)

     ENDDO  ! kproma

     zfrl(:,jk) = MERGE(ztmp1_1d(:), zfrl(:,jk) , lo_1d(:))

     ztmp1_1d(:) = MIN(ztmp2_1d(:), zcdnc(:,jk)-cdncmin)
     ztmp1_1d(:) = MAX(ztmp1_1d(:), 0._dp)
     zfrln(:) = MERGE(ztmp1_1d(:), zfrln(:), lo_1d(:))

     ztmp1_1d(:) = zcdnc(:,jk)-zfrln(:)
     ztmp1_1d(:) = MAX(ztmp1_1d(:), cqtmin)
     zcdnc(:,jk) = MERGE(ztmp1_1d(:), zcdnc(:,jk), lo_1d(:))

     ztmp1_1d(:) = zicnc(:,jk)+zfrln(:)
     ztmp1_1d(:) = MAX(ztmp1_1d(:), cqtmin)
     zicnc(:,jk) = MERGE(ztmp1_1d(:), zicnc(:,jk), lo_1d(:))

     ztmp1_1d(:) = zxlb(:)-zfrl(:,jk)
     zxlb(:) = MERGE(ztmp1_1d(:), zxlb(:), lo_1d(:))

     ztmp1_1d(:) = zxib(:)+zfrl(:,jk)
     zxib(:) = MERGE(ztmp1_1d(:), zxib(:), lo_1d(:))

     ztmp1_1d(:) = zfrl(:,jk) * zclcaux(:)
     zfrl(:,jk) = MERGE(ztmp1_1d(:), zfrl(:,jk), lo_1d(:))

     ! Bergeron-Findeisen process:
     ll1_1d(:) = lo_1d(:)   .AND. &  ! T_hom < T < 0
          locc_1d(:)        .AND. &  ! cloud cover > 0
          (zdep(:) > 0._dp) .AND. &  ! deposition rate > 0
          (zxlb(:) > 0._dp) .AND. &  ! in-cloud LWP > 0
          (0.01_dp * zvervx_2d(:,jk) < zvervmax_1d(:))  ! 0.01 < v_vert < v_max

     ztmp1_1d(:) = ztmst_rcp * zxlb(:) * zclcaux(:)
     ztmp2_1d(:) = pxlte(:,jk) - ztmp1_1d(:)
     pxlte(:,jk) = MERGE(ztmp2_1d(:), pxlte(:,jk), ll1_1d(:))

     ztmp2_1d(:) = pxite(:,jk) + ztmp1_1d(:)
     pxite(:,jk) = MERGE(ztmp2_1d(:), pxite(:,jk), ll1_1d(:))

     ztmp2_1d(:) = ptte(:,jk) + (zlsdcp(:) - zlvdcp(:)) * ztmp1_1d(:)
     ptte(:,jk) = MERGE(ztmp2_1d(:), ptte(:,jk), ll1_1d(:))

     zcdnc(:,jk) = MERGE(cqtmin, zcdnc(:,jk), ll1_1d(:))

     ztmp2_1d(:) = zxib(:) + zxlb(:)
     zxib(:) = MERGE(ztmp2_1d(:), zxib(:), ll1_1d(:))  ! water --> ice
     zxlb(:) = MERGE(0._dp, zxlb(:), ll1_1d(:))  ! water --> 0

     IF (lookupoverflow) THEN
        status_string = 'lookuperror: cdnc - cloud (1)'
        RETURN
     ENDIF


     ! === 7 - CLOUD PHYSICS AND PRECIPITATION FLUXES AT THE SURFACE ==========
     ! ========================================================================

     zclcstar_1d(:) = MIN(zclcaux(:), zclcpre(:)) ! cloud cover / precip. cover
     zauloc_1d(:) = 3. / 5000._dp * zdz_2d(:,jk)
     zauloc_1d(:) = MAX(MIN(zauloc_1d(:), clmax), clmin)

     ll1_1d(:) = (knvb(:) >= jbmin) .AND. &  ! Above lowest inversion level
                 (knvb(:) <= jbmax) .AND. &  ! Below highest inversion level
                 (pvervel(:,jk) > 0._dp)     ! v_ver > 0

     ll2_1d(:) = (jk == knvb(:)) .OR. (jk == knvb(:) + 1)

     ll3_1d(:) = ll1_1d(:) .AND. ll2_1d(:) .AND. lonacc

     zauloc_1d(:) = MERGE(0._dp, zauloc_1d(:), ll3_1d(:))

     zxlb(:) = MAX(zxlb(:), 1.e-20_dp)
     zxib(:) = MAX(zxib(:), 1.e-20_dp)

     ! Liquid water and snow content are stored before the reduction by
     ! outfalling rain (necessary for nucleation scavenging)
     plwc(:,jk) = zxlb(:)
     piwc(:,jk) = zxib(:)

     zmlwc(:,jk) = zxlb(:)
     zmiwc(:,jk) = zxib(:)

     ! Calculate the rain/snow water content [kg/kg] from the rain/snow flux
     ll1_1d(:) = (zclcpre(:) > zeps)  ! Fraction of gridbox covered by precip.
     ll2_1d(:) = ll1_1d(:) .AND. (zrfl(:) > cqtmin)
     ll3_1d(:) = ll1_1d(:) .AND. (zsfl(:) > cqtmin)

     ztmp1_1d(:) = (MAX(zrfl(:), cqtmin) / (12.45_dp * MAX(zclcpre(:), zeps) * &
          SQRT(zqrho_2d(:,jk))))**(8._dp/9._dp)
     ztmp2_1d(:) = (MAX(zsfl(:), cqtmin) / &
          (cvtfall * MAX(zclcpre(:), zeps)))**(1._dp/1.16_dp)

     zxrp1_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll2_1d(:))
     zxsp1_1d(:) = MERGE(ztmp2_1d(:), 0._dp, ll3_1d(:))

     ! --- 7.1 - Warm clouds --------------------------------------------------
     ! ------------------------------------------------------------------------

     ! Autoconversion of clouds droplets and collection of cloud droplets by
     ! falling rain

     ll1_1d(:) = locc_1d(:)         .AND. &  ! Cloud cover > 0
                 (zxlb(:) > cqtmin) .AND. &  ! LWP > 0
                 (zcdnc(:,jk) >= cdncmin)    ! CDNC >= CDNC_min

     SELECT CASE(nauto)

     CASE(1) ! Beheng (Atmospheric Res., 1994)

        ztmp1_1d(:) = ccraut * 1.2e27_dp * zrho_rcp(:,jk) * &
             (zcdnc(:,jk) * 1.e-6_dp)**(-3.3_dp) * &
             (zrho(:,jk) * 1.e-3_dp)**4.7_dp
        zraut_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))
        zraut_1d(:) = zxlb(:) * (1._dp - &
             (1._dp + ztmst * zexm1_2 * zraut_1d(:) * &
             zxlb(:)**zexm1_2)**zexp_2)
        zraut_1d(:) = MIN(zxlb(:), zraut_1d(:))

        ztmp1_1d(:) = 7.7e9_dp * zraut_1d(:) * zrho(:,jk)
        ztmp2_1d(:) = 1.289e10_dp * 1.e-6_dp * ztmst * &
             (zrho(:,jk)*zxlb(:))**2
        ztmp3_1d(:) = ztmp1_1d(:) + ztmp2_1d(:)
        ztmp3_1d(:) = MIN(ztmp3_1d(:), zcdnc(:,jk))
        zrautself_1d(:) = MERGE(ztmp3_1d(:), 0._dp, ll1_1d(:))

        ztmp1_1d(:) = zcdnc(:,jk) - zrautself_1d(:)
        ztmp1_1d(:) = MAX(ztmp1_1d(:), cqtmin)
        zcdnc(:,jk) = MERGE(ztmp1_1d(:), zcdnc(:,jk), ll1_1d(:))

        ztmp1_1d(:) = zxlb(:) - zraut_1d(:)
        zxlb(:)     = MERGE(ztmp1_1d(:), zxlb(:), ll1_1d(:))

        ! Accretion with rain from above
        zrac1_1d(:) = -6._dp * ztmst * zxrp1_1d(:)
        zrac1_1d(:) = EXP(zrac1_1d(:))
        zrac1_1d(:) = zxlb(:) * (1._dp-zrac1_1d(:))
        ztmp1_1d(:) = zxlb(:) - zrac1_1d(:)
        ztmp2_1d(:) = zxlb(:)
        zxlb(:) = MERGE(ztmp1_1d(:), zxlb(:), ll1_1d(:))  ! update LWC

        ! Accretion with newly formed rain inside the gridbox
        zrac2_1d(:) = &
             -6._dp * ztmst * zauloc_1d(:) * zrho(:,jk) * zraut_1d(:)
        zrac2_1d(:) = EXP(zrac2_1d(:))
        zrac2_1d(:) = zxlb(:) * (1._dp-zrac2_1d(:))
        ztmp1_1d(:) = zxlb(:)-zrac2_1d(:)
        zxlb(:)     = MERGE(ztmp1_1d(:), zxlb(:), ll1_1d(:))  ! update LWC

        ztmp1_1d(:) = zrpr(:) + zclcaux(:) * (zraut_1d(:) + zrac2_1d(:)) + &
             zclcstar_1d(:) * zrac1_1d(:)
        zrpr(:) = MERGE(ztmp1_1d(:), zrpr(:), ll1_1d(:))

        ztmp1_1d(:) = zraut_1d(:) + zrac1_1d(:) + zrac2_1d(:)
        zmratepr(:,jk) = MERGE(ztmp1_1d(:), zmratepr(:,jk), ll1_1d(:))

        ! Autoconversion also affects the number of cloud droplets (zrprn)
        ztmp1_1d(:) = (zrac1_1d(:) + zrac2_1d(:)) / (ztmp2_1d(:) + zeps)
        ztmp1_1d(:) = zcdnc(:,jk) * ztmp1_1d(:)
        ztmp1_1d(:) = MIN(ztmp1_1d(:), zcdnc(:,jk))
        ztmp2_1d(:) = zcdnc(:,jk)-ztmp1_1d(:)
        ztmp2_1d(:) = MAX(ztmp2_1d(:), cqtmin)
        zcdnc(:,jk) = MERGE(ztmp2_1d(:), zcdnc(:,jk), ll1_1d(:))
        ztmp2_1d(:) = zrautself_1d(:) + ztmp1_1d(:)
        zrprn(:) = MERGE(ztmp2_1d(:), zrprn(:), ll1_1d(:))

     CASE(2) ! Khairoutdinov and Kogan (Mon. Weather Rev., 2000)

        ! Eq. (32) documentation
        ztmp1_1d(:) = ccraut * 1350._dp * (1.e-6_dp * zcdnc(:,jk))**(-1.79_dp)

        ztmp1_1d(:) = zxlb(:) * (1._dp - &
             (1._dp + ztmst * zexm1_1 * ztmp1_1d(:) * &
                 zxlb(:)**zexm1_1)**zexp_1)

        ztmp1_1d(:) = MIN(zxlb(:), ztmp1_1d(:))
        zraut_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))

        ztmp1_1d(:) = zxlb(:) - zraut_1d(:)
        ztmp2_1d(:) = zxlb(:)
        zxlb(:) = MERGE(ztmp1_1d(:), zxlb(:), ll1_1d(:))

        ! Accretion with rain from above
        ztmp1_1d(:) = -3.7_dp * ztmst * zxrp1_1d(:)
        ztmp1_1d(:) = EXP(ztmp1_1d(:))
        ztmp1_1d(:) = zxlb(:) * (1._dp - ztmp1_1d(:))
        zrac1_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))
        zxlb(:) = zxlb(:) - zrac1_1d(:)  ! update LWC

        ! Accretion with newly formed rain inside the gridbox
        ztmp1_1d(:) = &
             -3.7_dp * ztmst * zauloc_1d(:) * zrho(:,jk) * zraut_1d(:)
        ztmp1_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))
        ztmp1_1d(:) = zxlb(:) * (1._dp - EXP(ztmp1_1d(:)))
        zrac2_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))
        zxlb(:) = zxlb(:) - zrac2_1d(:)  ! update LWC

        zrpr(:) = zrpr(:) + zclcaux(:) * (zraut_1d(:) + zrac2_1d(:)) + &
             zclcstar_1d(:) * zrac1_1d(:)

        ztmp1_1d(:) = zraut_1d(:) + zrac1_1d(:) + zrac2_1d(:)
        zmratepr(:,jk) = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))

        ! Autoconversion also affects the number of cloud droplets (zrprn)
        ztmp1_1d(:) = &
             (zraut_1d(:) + zrac1_1d(:) + zrac2_1d(:)) / (ztmp2_1d(:) + zeps)
        zrprn(:) = MERGE(ztmp1_1d(:), zrprn(:), ll1_1d(:))

        ll2_1d(:) = ll1_1d(:) .AND. (zxlb(:) > cqtmin)

        ztmp1_1d(:) = MERGE(cdncmin, 0._dp, ll2_1d(:))
        ztmp1_1d(:) = zcdnc(:,jk)-ztmp1_1d(:)
        ztmp2_1d(:) = zcdnc(:,jk)*zrprn(:)
        ztmp3_1d(:) = MIN(ztmp1_1d(:), ztmp2_1d(:))
        zrprn(:) = MERGE(ztmp3_1d(:), zrprn(:), ll1_1d(:))

        ztmp1_1d(:) = zcdnc(:,jk)-zrprn(:)
        ztmp1_1d(:) = MAX(ztmp1_1d(:), cqtmin)
        zcdnc(:,jk) = MERGE(ztmp1_1d(:), zcdnc(:,jk), ll1_1d(:))

     CASE DEFAULT

        WRITE(*,*) "nauto can be either 1 or 2 in messy_cloud_Kuebbeler14"
        STOP

     END SELECT

     ! --- 7.2 - Cold clouds --------------------------------------------------
     ! ------------------------------------------------------------------------

     ! Conversion of cloud ice to snow after Levkov et al. (1992)
     ! Aggregation of ice crystals to snow assuming plates (zsaut) and
     ! accretion of ice by falling snow (zsaci)
     ! Accretion of cloud droplets by falling snow. (zsacl)

     ll1_1d(:) = locc_1d(:) .AND. (zxib(:) > cqtmin)

     ! Effective radius of IC assuming plates, Eq. (59) documentation
     ztmp1_1d(:) = 0.5e4_dp * (1000._dp / 0.0376_dp * zxib(:) * zrho(:,jk) &
          / zicnc(:,jk))**0.302  ! [um]
     ztmp1_1d(:) = MIN(MAX(ztmp1_1d(:), ceffmin), ceffmax)

     ztmp2_1d(:) = 5113188._dp + 2809._dp * ztmp1_1d(:)**3
     ztmp2_1d(:) = SQRT(ztmp2_1d(:))
     ztmp2_1d(:) = -2261._dp + ztmp2_1d(:)

     ztmp3_1d(:) = 1.e-6_dp * ztmp2_1d(:)**(1._dp/3._dp)
     zris_1d(:) = MERGE(ztmp3_1d(:), 1._dp, ll1_1d(:))

     ! Temperature dependent collision efficiency (zcolleffi)
     ztmp1_1d(:) = 0.025_dp * (ztp1tmp(:) - tmelt)
     ztmp1_1d(:) = EXP(ztmp1_1d(:))
     zcolleffi_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))

     zc1_1d(:) = 17.5_dp / crhoi * zrho(:,jk) * zqrho_2d(:,jk)**0.33_dp

     ztmp1_1d(:) = -6._dp / zc1_1d(:) * LOG10(1.e4_dp*zris_1d(:))
     ztmp1_1d(:) = ccsaut / ztmp1_1d(:)

     zsaut_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))
     zsaut_1d(:) = &
          zxib(:) * (1._dp - 1._dp / (1._dp + zsaut_1d(:) * ztmst * zxib(:)))

     ztmp1_1d(:) = zxib(:) - zsaut_1d(:)
     zxibold_1d(:) = zxib(:)
     zxib(:) = MERGE(ztmp1_1d(:), zxib(:), ll1_1d(:))

     zsaci2_1d(:) = 0.0_dp
     zsacl2_1d(:) = 0.0_dp

     ! Aaccretion also affects the number of cloud droplets
     zsacl2in_1d(:)  = 0.0_dp

     ztmp1_1d(:) = zauloc_1d(:)*zrho(:,jk)*zsaut_1d(:)

     ! Snow formed inside the grid box
     zxsp2_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))
     zxsp_1d(:) = zxsp1_1d(:) + zxsp2_1d(:)

     ll2_1d(:) = ll1_1d(:)               .AND. &
                 (zxsp_1d(:)  >  cqtmin) .AND. &
                 (zxlb(:)     >  cqtmin) .AND. &
                 (zcdnc(:,jk) >= cdncmin)

     ! The riming of snow with droplets is calculated assuming planar snow
     ! flakes (Lohmann, J. Atmos. Sci., 2004). It depends on the droplet
     ! (zudrop) and snow flake (zusnow) fall velocity, the Stokes number
     ! (zstokes) and the Reynolds number (zrey)

     zscnc_1d(:) = 0.5_dp * zicnc(:,jk)
     IF (jk > 1) THEN
        zscnc_1d(:) = MAX(zicemin, MIN(zsprn(:,jk-1), zscnc_1d(:)))
     ELSE
        zscnc_1d(:) = MAX(zicemin, zscnc_1d(:))
     ENDIF

     ztmp1_1d(:) = (6._dp * zpirho_rcp * zrho(:,jk) * zxlb(:) / &
          zcdnc(:,jk) )**(1._dp/3._dp)
     ztmp1_1d(:) = MAX(ztmp1_1d(:), 1.e-6_dp) !SF zdw

     zudrop_1d(:) = 1.19e4_dp * 2500._dp * ztmp1_1d(:)**2 * &
          (1.3_dp * zrho_rcp(:,jk))**0.35_dp

     zdplanar_1d(:) = 1.e3_dp / 3.8e-4_dp * zxsp_1d(:) / zscnc_1d(:)
     zdplanar_1d(:) = SQRT(zdplanar_1d(:))
     zdplanar_1d(:) = 1.e-2_dp * zdplanar_1d(:)
     zdplanar_1d(:) = MAX(20.e-6_dp, zdplanar_1d(:))

     zusnow_1d(:) = 2.34_dp * (100._dp * zdplanar_1d(:))**0.3_dp * &
          (1.3_dp * zrho_rcp(:,jk))**0.35_dp

     zstokes_1d(:) = 2._dp * g_rcp * (zusnow_1d(:) - zudrop_1d(:)) * &
          zudrop_1d(:) / zdplanar_1d(:)
     zstokes_1d(:) = MAX(zstokes_1d(:), cqtmin)

     zrey_1d(:) = zrho(:,jk) * zdplanar_1d(:) * zusnow_1d(:) / zviscos_2d(:,jk)
     zrey_1d(:) = MAX(zrey_1d(:), cqtmin)

     ll3_1d(:) = (zrey_1d(:) <=  5._dp)
     ll4_1d(:) = (zrey_1d(:) >   5._dp) .AND. (zrey_1d(:) <  40._dp)
     ll5_1d(:) = (zrey_1d(:) >= 40._dp)

     ztmp1_1d(:) = 5.52_dp * zrey_1d(:)**(-1.12_dp)
     ztmp2_1d(:) = 1.53_dp * zrey_1d(:)**(-0.325_dp)
     zstcrit_1d(:) = 1._dp
     zstcrit_1d(:) = MERGE(ztmp1_1d(:), zstcrit_1d(:), ll3_1d(:))
     zstcrit_1d(:) = MERGE(ztmp2_1d(:), zstcrit_1d(:), ll4_1d(:))

     zcsacl_1d(:) = 0.2_dp * &
          (LOG10(zstokes_1d(:)) - LOG10(zstcrit_1d(:)) - 2.236_dp )**2
     zcsacl_1d(:) = MIN(zcsacl_1d(:), 1._dp - cqtmin)
     zcsacl_1d(:) = MAX(zcsacl_1d(:), 0._dp)
     zcsacl_1d(:) = SQRT(1._dp - zcsacl_1d(:))

     ll6_1d(:) = ll5_1d(:) .AND. (zstokes_1d(:) <= 0.06_dp)
     ll7_1d(:) = ll5_1d(:) .AND. (zstokes_1d(:) >  0.06_dp) .AND. &
          (zstokes_1d(:) <= 0.25_dp)
     ll8_1d(:) = ll5_1d(:) .AND. (zstokes_1d(:) >  0.25_dp) .AND. &
          (zstokes_1d(:) <= 1.00_dp)

     WHERE (ll6_1d(:))
        zcsacl_1d(:) = 1.034_dp * zstokes_1d(:)**1.085_dp
     ELSEWHERE (ll7_1d(:))
        zcsacl_1d(:) = 0.787_dp * zstokes_1d(:)**0.988_dp
     ELSEWHERE (ll8_1d(:))
        zcsacl_1d(:) = 0.7475_dp * LOG10(zstokes_1d(:))+0.65_dp
     ELSEWHERE (ll5_1d(:))
        zcsacl_1d(:) = (zstokes_1d(:)+1.1_dp)**2 / (zstokes_1d(:) + 1.6_dp)**2
     ENDWHERE

     zcsacl_1d(:) = MAX(MIN(zcsacl_1d(:), 1._dp), 0.01_dp)
     zcsacl_1d(:)  = MERGE(zcsacl_1d(:), 0._dp, ll2_1d(:))

     ztmp1_1d(:)  = zcons4 * zxsp_1d(:)**0.8125_dp
     ztmp1_1d(:)  = api * cn0s * 3.078_dp * ztmp1_1d(:) * zqrho_2d(:,jk)**0.5_dp
     zsaci2_1d(:) = MERGE(ztmp1_1d(:), zsaci2_1d(:), ll2_1d(:))

     ztmp1_1d(:)  = -ztmst * zsaci2_1d(:) * zcsacl_1d(:)
     ztmp1_1d(:)  = EXP(ztmp1_1d(:))
     ztmp1_1d(:)  = zxlb(:) * (1._dp-ztmp1_1d(:))

     zsacl2in_1d(:) = MERGE(ztmp1_1d(:), zsacl2in_1d(:), ll2_1d(:))

     ztmp2_1d(:) = zxlb(:)-ztmp1_1d(:)
     ztmp3_1d(:) = zxlb(:)
     zxlb(:) = MERGE(ztmp2_1d(:), zxlb(:), ll2_1d(:))

     ztmp2_1d(:) = zclcaux(:) * ztmp1_1d(:)
     zsacl2_1d(:)   = MERGE(ztmp2_1d(:), zsacl2_1d(:), ll2_1d(:))

     ll2_1d(:) = ll1_1d(:)             .AND. &
                 (zxsp_1d(:) > cqtmin) .AND. &
                 (zxib(:)    > cqtmin)

     ztmp1_1d(:) = zcons4 * zxsp_1d(:)**0.8125_dp
     ztmp1_1d(:) = api * cn0s * 3.078_dp * ztmp1_1d(:) * zqrho_2d(:,jk)**0.5_dp
     ztmp1_1d(:) = -ztmst * ztmp1_1d(:) * zcolleffi_1d(:)
     ztmp1_1d(:) = EXP(ztmp1_1d(:))
     ztmp1_1d(:) = zxib(:) * (1._dp-ztmp1_1d(:))

     zsaci2_1d(:) = MERGE(ztmp1_1d(:), zsaci2_1d(:), ll2_1d(:))

     zxib(:) = zxib(:)-zsaci2_1d(:)

     zsacl(:) = MERGE(zsacl2_1d(:), zsacl(:), ll1_1d(:))

     ztmp1_1d(:) = zspr(:) + zclcaux(:) * (zsaut_1d(:) + zsaci2_1d(:))
     zspr(:)  = MERGE(ztmp1_1d(:), zspr(:), ll1_1d(:))

     ztmp1_1d(:) = zsaut_1d(:) + zsaci2_1d(:)
     zmrateps(:,jk) = MERGE(ztmp1_1d(:), zmrateps(:,jk), ll1_1d(:))

     ll2_1d(:) = (zxlb(:) > cqtmin)

     ztmp1_1d(:) = zcdnc(:,jk) * zsacl2in_1d(:) / (ztmp3_1d(:) + zeps)
     ztmp1_1d(:) = MIN(ztmp1_1d(:), zcdnc(:,jk) - cdncmin)
     ztmp1_1d(:) = MAX(ztmp1_1d(:), 0._dp)

     ztmp2_1d(:) = zcdnc(:,jk) * zsacl2in_1d(:) / (ztmp3_1d(:) + zeps)
     ztmp2_1d(:) = MIN(ztmp2_1d(:), zcdnc(:,jk))

     ztmp3_1d(:)  = MERGE(ztmp1_1d(:), ztmp2_1d(:), ll2_1d(:))
     zsacln(:) = MERGE(ztmp3_1d(:), zsacln(:), ll1_1d(:))

     ztmp1_1d(:)     = zcdnc(:,jk)-zsacln(:)
     zcdnc(:,jk)     = MERGE(ztmp1_1d(:), zcdnc(:,jk), ll1_1d(:))
     zmsnowacl(:,jk) = MERGE(zsacl2in_1d(:), zmsnowacl(:,jk), ll1_1d(:))

     ! Secondary ice crystal production (zsecprod) after Levkov et al. (1992).
     ! Sink for snow, source for ice crystals.
     ! Accretion rate is size dependent (Lohmann, J. Atmos. Sci., 2004)

     zsecprod_1d(:) = 0._dp

     ll2_1d(:) = ll1_1d(:)               .AND. &
                 (zxsp_1d(:) > zepsec)   .AND. &
                 (zxlb(:)    > zepsec)   .AND. &
                 (ztp1tmp(:) > 265.2_dp) .AND. &
                 (ztp1tmp(:) < 270.2_dp)

     ztmp1_1d(:) = (6._dp * zpirho_rcp * zrho(:,jk) * zxlb(:) / &
          zcdnc(:,jk))**(1._dp / 3._dp)
     ztmp1_1d(:) = MAX(ztmp1_1d(:), 1.e-6_dp)

     zudrop_1d(:) = 1.19e4_dp * (50._dp * ztmp1_1d(:))**2 * &
          (1.3_dp * zrho_rcp(:,jk))**0.35_dp

     zstokes_1d(:) = 2._dp * g_rcp * (zusnow_1d(:) - zudrop_1d(:)) * &
          zudrop_1d(:) / zdplanar_1d(:)
     zstokes_1d(:) = MAX(zstokes_1d(:), cqtmin)

     ztmp1_1d(:) = &
          0.2_dp * (LOG10(zstokes_1d(:)) - LOG10(zstcrit_1d(:)) - 2.236_dp )**2
     ztmp1_1d(:) = MIN(ztmp1_1d(:), 1._dp-cqtmin)
     ztmp1_1d(:) = MAX(ztmp1_1d(:), 0._dp)
     ztmp1_1d(:) = SQRT(1._dp - ztmp1_1d(:))

     ll6_1d(:) = ll5_1d(:) .AND. (zstokes_1d(:) <= 0.06_dp)
     ll7_1d(:) = ll5_1d(:) .AND. (zstokes_1d(:) >  0.06_dp) .AND. &
          (zstokes_1d(:) <= 0.25_dp)
     ll8_1d(:) = ll5_1d(:) .AND. (zstokes_1d(:) >  0.25_dp) .AND. &
          (zstokes_1d(:) <= 1.00_dp)

     WHERE (ll6_1d(:))
        ztmp1_1d(:) = 1.034_dp * zstokes_1d(:)**1.085_dp
     ELSEWHERE (ll7_1d(:))
        ztmp1_1d(:) = 0.787_dp * zstokes_1d(:)**0.988_dp
     ELSEWHERE (ll8_1d(:))
        ztmp1_1d(:) = 0.7475_dp * LOG10(zstokes_1d(:))+0.65_dp
     ELSEWHERE (ll5_1d(:))
        ztmp1_1d(:) = (zstokes_1d(:)+1.1_dp)**2/(zstokes_1d(:)+1.6_dp)**2
     ENDWHERE

     ztmp1_1d(:)  = MAX(MIN(ztmp1_1d(:), 1._dp), 0.01_dp)
     zcsacl_1d(:) = MERGE(ztmp1_1d(:), zcsacl_1d(:), ll2_1d(:))

     ztmp1_1d(:) = zcons5 * zxsp_1d(:)**0.875_dp

     ztmp2_1d(:) = cn0s * 0.831_dp * api / zmw0 * zcsacl_1d(:) * zrho(:,jk) * &
          zxlb(:) * ztmp1_1d(:) * &
          (g * crhosno / (0.75_dp * zcdi * zrho(:,jk)) )**0.5_dp

     ztmp2_1d(:) = MAX(0.00285_dp * ztmp2_1d(:), 0._dp)

     ztmp2_1d(:) = ztmst * zmi0 * ztmp2_1d(:) * zrho_rcp(:,jk)
     ztmp3_1d(:) = zxsp_1d(:) * zrho_rcp(:,jk)
     ztmp2_1d(:) = MIN(ztmp3_1d(:), ztmp2_1d(:))
     ztmp2_1d(:) = MAX(ztmp2_1d(:), 0._dp)

     zsecprod_1d(:) = MERGE(ztmp2_1d(:), zsecprod_1d(:), ll2_1d(:))

     ztmp1_1d(:) = zxib(:) + zsecprod_1d(:)
     zxib(:) = MERGE(ztmp1_1d(:), zxib(:), ll2_1d(:))

     ztmp1_1d(:) = zspr(:) - zclcstar_1d(:) * zsecprod_1d(:)
     ztmp1_1d(:) = MAX(ztmp1_1d(:), 0._dp)
     zspr(:) = MERGE(ztmp1_1d(:), zspr(:), ll2_1d(:))

     ztmp1_1d(:) = zmrateps(:,jk) - zsecprod_1d(:)
     zmrateps(:,jk) = MERGE(ztmp1_1d(:), zmrateps(:,jk), ll2_1d(:))

     ! Store of snow production before it is transformed into a flux
     prate_s(:,jk) = zsaut_1d(:) + zsaci2_1d(:)


     ! Change the number of ice crystals due to the break-up of snow flakes
     ll1_1d(:) = locc_1d(:)               .AND. &
                 (zxib(:)     > zepsec  ) .AND. &
                 (zicnc(:,jk) >= zicemin)

     ztmp1_1d(:) = zxibold_1d(:)
     ztmp1_1d(:) = MAX(ztmp1_1d(:), 0._dp)

     ztmp2_1d(:) = &
          zicnc(:,jk) * (zsaci2_1d(:) + zsaut_1d(:)) / (ztmp1_1d(:) + zeps)

     ztmp3_1d(:) = 0.5_dp * ztmst * zc1_1d(:) * zicnc(:,jk) * zxib(:)

     ztmp4_1d(:) = zmi0_rcp * zrho(:,jk) * zsecprod_1d(:)

     ztmp5_1d(:) = ztmp2_1d(:) + ztmp3_1d(:) - ztmp4_1d(:)
     ztmp5_1d(:) = MIN(ztmp5_1d(:), zicnc(:,jk))
     zsprn(:,jk) = MERGE(ztmp5_1d(:), zsprn(:,jk), ll1_1d(:))

     ztmp1_1d(:) = zicnc(:,jk) - zsprn(:,jk)
     ztmp1_1d(:) = MAX(ztmp1_1d(:), cqtmin)
     zicnc(:,jk) = MERGE(ztmp1_1d(:), zicnc(:,jk), ll1_1d(:))


     ! --- 7.3 - Update precipitation fluxes ----------------------------------
     ! ------------------------------------------------------------------------

     ! In the lowest layer (klev), the sedimentation sink of cloud ice is
     ! balanced by precipitation at the surface (through 'zzdrs'). The fraction
     ! of precipitating clouds (zclcpre) is used for the calculation of
     ! evaporation/sublimation of rain/snow in the next layer

     zzdrr_1d(:)    = zcons2*zdp_2d(:,jk)*zrpr(:)
     zzdrs_1d(:)    = zcons2*zdp_2d(:,jk)*(zspr(:)+zsacl(:))

     IF (jk.EQ.klev) THEN
        zzdrs_1d(:) = zzdrs_1d(:) + zxiflux(:)
        ztmp1_1d(:) = zcons2 * zdp_2d(:,jk) / (zlsdcp(:) - zlvdcp(:)) * &
             MAX(0._dp, (ztp1tmp(:) - tmelt))
        ztmp2_1d(:) = MIN(zxsec * zzdrs_1d(:), ztmp1_1d(:))
        zzdrr_1d(:) = zzdrr_1d(:) + ztmp2_1d(:)
        zzdrs_1d(:) = zzdrs_1d(:) - ztmp2_1d(:)
        zsmlt(:) = zsmlt(:) + ztmp2_1d(:) / (zcons2 * zdp_2d(:,jk))
     ENDIF

     zpretot_1d(:) = zrfl(:) + zsfl(:)
     zpredel_1d(:) = zzdrr_1d(:) + zzdrs_1d(:)

     ll1_1d(:) = (zpretot_1d(:) > zpredel_1d(:))

     zclcpre(:) = MERGE(zclcpre(:), zclcaux(:), ll1_1d(:))

     zpresum_1d(:)  = zpretot_1d(:)+zpredel_1d(:)

     ll1_1d(:) = (zpresum_1d(:) > cqtmin)

     ztmp1_1d(:) = (zclcaux(:) * zpredel_1d(:) + &
                    zclcpre(:) * zpretot_1d(:)) / &
                    MAX(zpresum_1d(:), cqtmin)
     ztmp1_1d(:) = MIN(ztmp1_1d(:), 1.0_dp)
     ztmp1_1d(:) = MAX(ztmp1_1d(:), 0.0_dp)

     zclcpre(:) = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))

     ll1_1d(:) = (zclcpre(:) > zepsec)
     ztmp1_1d(:) = (zrfl(:) + zzdrr_1d(:)) / MAX(zclcpre(:), zepsec)
     ztmp2_1d(:) = (zsfl(:) + zzdrs_1d(:)) / MAX(zclcpre(:), zepsec)
     ztmp3_1d(:) = (zcons2 * zdp_2d(:,jk) * zevp(:)) / MAX(zclcpre(:), zepsec)
     ztmp4_1d(:) = (zcons2 * zdp_2d(:,jk) * zsub(:)) / MAX(zclcpre(:), zepsec)

     zfrain(:,jk)  = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))
     zfsnow(:,jk)  = MERGE(ztmp2_1d(:), 0._dp, ll1_1d(:))
     zfevapr(:,jk) = MERGE(ztmp3_1d(:), 0._dp, ll1_1d(:))
     zfsubls(:,jk) = MERGE(ztmp4_1d(:), 0._dp, ll1_1d(:))


     ! Rain and snow flux considering incoming rain, melting of snow, droplet
     ! evaporation/sublimation, but no new production of rain or snow in that
     ! layer.... (neccessary for impaction scavenging)
     pfrain_no(:,jk)   = zrfl(:) - zcons2 * zdp_2d(:,jk) * zevp(:)
     pfsnow_no(:,jk)   = zsfl(:) - zcons2 * zdp_2d(:,jk) * zsub(:)

     ! Precipitating cloud cover of this layer is used for the next lower layer
     ! to estimate the part of the cloud cover in which rain impacts
     pr_cover(:,jk) = zclcpre(:)

     zrfl(:) = zrfl(:) + zzdrr_1d(:) - zcons2 * zdp_2d(:,jk) * zevp(:)
     zsfl(:) = zsfl(:) + zzdrs_1d(:) - zcons2 * zdp_2d(:,jk) * zsub(:)


     ! === 8 - UPDATE TENDENCIES OF T, Q, XL, XI AND FINAL CLOUD COVER ========
     ! ========================================================================

     ! --- 8.1 - Cloud cover scheme tendencies --------------------------------
     ! ------------------------------------------------------------------------

     IF (lcover .AND. jk >= ncctop) THEN

        ! Source terms from convection: skewness
        IF (ncvmicro) THEN
           zconvskew(:) = cbeta_cs * &
              (pxtecl(:,jk) + pxteci(:,jk) + pqtec(:,jk)) / pbetass(:,jk)
        ELSE
           zconvskew(:) = cbeta_cs * &
                (pxtec(:,jk) +pqtec(:,jk)) / pbetass(:,jk)
        ENDIF

        ztmp1_1d(:) = zdtime_rcp * (cbeta_pq_max - pxskew(:,jk))
        zconvskew(:) = MIN(zconvskew(:), ztmp1_1d(:))

        ! Convective width now diagnosed, assuming 'a' unchanged:
        ll1_1d(:) = (pqm1(:,jk) >= pbetass(:,jk))

        ztmp1_1d(:) = pxskew(:,jk) + zdtime * zconvskew(:)
        ztmp2_1d(:) = zwide(:) * (cbeta_pq+ztmp1_1d(:)) &
             / (cbeta_pq+pxskew(:,jk))
        ztmp3_1d(:) = zdtime_rcp*(ztmp2_1d(:)-zwide(:))
        zconvvar(:) = MERGE(ztmp3_1d(:), 0._dp, ll1_1d(:))

        ! Simple linearized effect of microphysics on skewness
        ll1_1d(:) = &
             (pbetaa(:,jk) < pbetass(:,jk)).AND.(pbetab(:,jk) > pbetass(:,jk))

        ztmp1_1d(:) = ztmst * (zxlte(:) + zxite(:)) - zrpr(:) - zsacl(:) - &
             zspr(:) + zcnd(:) + zdep(:) + zgenti(:) + zgentl(:)

        ztmp1_1d(:) = - ztmp1_1d(:) / MAX(zepsec, zbetacl(:))
        ztmp1_1d(:) = MIN(1._dp, ztmp1_1d(:))
        ztmp1_1d(:) = MAX(0._dp, ztmp1_1d(:))
        ztmp1_1d(:) = (pbetass(:,jk) - pbetab(:,jk)) * ztmp1_1d(:)

        ll2_1d(:) = (ABS(zbetaqt(:) - pbetaa(:,jk)) < zeps)
        ztmp4_1d(:) = MERGE(zeps, (zbetaqt(:) - pbetaa(:,jk)), ll2_1d(:))

        ztmp2_1d(:) = (pbetab(:,jk) + ztmp1_1d(:) - pbetaa(:,jk)) * &
             cbeta_pq / ztmp4_1d(:) - cbeta_pq
        ztmp2_1d(:) = MAX(MIN(ztmp2_1d(:), cbeta_pq_max), cbeta_pq)

        ztmp3_1d(:) = zdtime_rcp * (ztmp2_1d(:) - pxskew(:,jk))
        ztmp3_1d(:) = MIN(0._dp, ztmp3_1d(:))

        zmicroskew(:) = MERGE(ztmp3_1d(:), zmicroskew(:), ll1_1d(:))

        ! New skewness and variance
        zxskewte(:) = zconvskew(:) + zmicroskew(:) + zturbskew(:)
        zxvarte(:) = zconvvar(:) + zturbvar(:)

        ztmp1_1d(:) = pxvar(:,jk)  + zdtime * zxvarte(:)
        ztmp2_1d(:) = pxskew(:,jk) + zdtime * zxskewte(:)
        pxskew(:,jk)   = MAX(MIN(ztmp2_1d(:), cbeta_pq_max), cbeta_pq)
        ztmp3_1d(:) = zbetaqt(:) * (1._dp + pxskew(:,jk) / cbeta_pq)
        pxvar(:,jk) = MAX(MIN(ztmp1_1d(:), ztmp3_1d(:)), zvartg(:))

     ENDIF

     ! --- 8.2 - Tendencies of thermodynamic variables ------------------------
     ! ------------------------------------------------------------------------

     ! The terms zxisub and zximlt do not appear in pxite because these
     ! processes have already been included in pxite via changes in cloud ice
     ! sedimentation (see 3.1, 3.2 and 4)

     pqte(:,jk) = pqte(:,jk) + &
          ztmst_rcp * (-zcnd(:) - zgentl(:) + zevp(:) + zxlevap(:) -&
           zdep(:) - zgenti(:) + zsub(:) + zxievap(:) + zxisub(:))

     ptte(:,jk) = ptte(:,jk) + &
          ztmst_rcp * (zlvdcp(:) * (zcnd(:) + zgentl(:) - zevp(:) - &
          zxlevap(:)) + zlsdcp(:) * (zdep(:) + zgenti(:) - zsub(:) - &
          zxievap(:) - zxisub(:)) + (zlsdcp(:) - zlvdcp(:)) * (-zsmlt(:) - &
          zimlt(:) - zximlt(:) + zfrl(:,jk) + zsacl(:)))

     ztmp1_1d(:) = pxlte(:,jk) + zxlte(:)

     ztmp2_1d(:) = zimlt(:) + zximlt(:) - zfrl(:,jk) - zrpr(:) - zsacl(:) + &
          zcnd(:)   + zgentl(:)  - zxlevap(:)

     zxlp1_1d(:) = pxlm1(:,jk) + ztmst * ztmp1_1d(:) + ztmp2_1d(:)

     pxlte(:,jk) = ztmp1_1d(:) + ztmst_rcp * ztmp2_1d(:)

     ztmp1_1d(:) = pxite(:,jk) + zxite(:)

     ztmp2_1d(:) = zfrl(:,jk) - zspr(:) + zdep(:) + zgenti(:) - zxievap(:) - &
          zsaut2_1d(:) ! zsaut2_1d extra term from IC to snow for R>100 um

     zxip1_1d(:) = pxim1(:,jk) + ztmst * ztmp1_1d(:)  + ztmp2_1d(:)
     pxite(:,jk) = ztmp1_1d(:) + ztmst_rcp * ztmp2_1d(:)

     ! Calculate new total tendency of CDNC
     ! NOTE: this is indeed the total tendency, since zcdnc is initialised
     ! above including the tendency, and here the difference to t-1
     ! is calculated
     pxtte(:,jk,1) = ztmst_rcp * (zcdnc(:,jk) * m_air / 1000._dp * &
          zrho_rcp(:,jk) - pxtm1(:,jk,1))  ! [1/mol/s]
     
     ! Update CDNC for radiation:
     pacdnc(:,jk) = zcdnc(:,jk)

     ! Calculate new total tendency of ICNC
     ! As above
     pxtte(:,jk,2) = ztmst_rcp * (zicnc(:,jk)* m_air/1000._dp * &
          zrho_rcp(:,jk) - pxtm1(:,jk,2))  ! [1/mol/s]
     
     ! Diagnostics
     qaut(1:kproma,jk,jrow) = qaut(1:kproma,jk,jrow) - zdt * zrprn(:)
     qfre(1:kproma,jk,jrow) = qfre(1:kproma,jk,jrow) - zdt * zfrln(:)
     qacc(1:kproma,jk,jrow) = qacc(1:kproma,jk,jrow) - zdt * zsacln(:)
     cloud_tm1(1:kproma,jk,jrow) = paclc(1:kproma,jk)

     ! Liquid cloud conditions
     ll1_1d(:) = (zxlb(:) > zeps) .AND. (zcdnc(:,jk) >= cdncmin)

     ! CDNC (in-cloud) [m-3]
     cdnc_acc(1:kproma,jk,jrow) = &
          MERGE(zcdnc(:,jk), cdnc_acc(1:kproma,jk,jrow), ll1_1d(:))

     ! CDNC burden (in-cloud) [m-2]
     ztmp1_1d(:) = zcdnc_burden(:) + zcdnc(:,jk) * zdz_2d(:,jk)
     zcdnc_burden(:) = MERGE(ztmp1_1d(:), zcdnc_burden(:), ll1_1d(:))
     CDNC_burden_acc(1:kproma,jrow) = zcdnc_burden(:)

     ! CDNC (grid-box) [m-3]
     cdnc(1:kproma,jk,jrow) = &
          MERGE(zcdnc(:,jk) * zclcaux(:), cdnc(1:kproma,jk,jrow), ll1_1d(:))

     ! CDNC burden (grid-box) [m-2]
     cdnc_burden(1:kproma,jrow) = &
          MERGE(cdnc_burden(1:kproma,jrow) + zcdnc(:,jk) * zdz_2d(:,jk) * &
          zclcaux(:), cdnc_burden(1:kproma,jrow), ll1_1d(:))

     ! Cloud droplet effective radius [um]
     ! (Peng and Lohmann, Geophys. Res. Lett., 2003)
     ztmp1_1d(:)  = 0.00045e-6_dp * zcdnc(:,jk) + 1.18_dp  ! Eq. (6)
     zreffl_1d(:) = 1.E6_dp * ztmp1_1d(:) * ((3._dp / (4._dp *api * rhoh2o)) * &
          zxlb(:) * zrho(:,jk) / zcdnc(:,jk))**(1._dp/3._dp)  ! Eq. (1)

     zreffl_1d(:) = MAX(4._dp, MIN(zreffl_1d(:), 40._dp))
     reffl(1:kproma,jk,jrow) = &
          MERGE(zreffl_1d(:), reffl(1:kproma,jk,jrow), ll1_1d(:))

     ! Ice cloud condition
     ll1_1d(:) = (zxib(:) > zeps)

     ! ICNC (in-cloud) [m-3]
     icnc_acc(1:kproma,jk,jrow) = &
          MERGE(zicnc(:,jk), icnc_acc(1:kproma,jk,jrow), ll1_1d(:))

     ! Ice crystal effective radius [um]
     ! Lohmann et al. (Atmos. Chem. Phys., 2007)
     ll2_1d(:) = (ztp1tmp(:) > cthomi)
     ztmp1_1d(:) = &  ! plate assumption - Eq. (59) documentation
          0.5e4_dp * (1000._dp / 0.0376_dp * MAX(zxib(:), zeps) * &
          zrho(:,jk) / zicnc(:,jk))**0.302_dp
     ztmp2_1d(:) = &  ! volume mean radius - Eq. (6) Lohmann et al. (2007)
          1.0e6_dp * (3._dp / (4._dp * api * zrhoice) * &
          MAX(zxib(:), zeps) * zrho(:,jk) / zicnc(:,jk))**0.333_dp
     ztmp2_1d(:) = &  ! Eq. (5) Lohmann et al. (2007)
          (1.61_dp * ztmp2_1d(:)**3 + 3.56e-4_dp * ztmp2_1d(:)**6)**0.333_dp

     ztmp3_1d(:) = MERGE(ztmp1_1d(:), ztmp2_1d(:), ll2_1d(:))
     ztmp3_1d(:) = MAX(ztmp3_1d(:), ceffmin)
     ztmp3_1d(:) = MIN(ztmp3_1d(:), ceffmax)

     reffi(1:kproma,jk,jrow) = &
          MERGE(ztmp3_1d(:), reffi(1:kproma,jk,jrow), ll1_1d(:))

     ! Ice cloud condition
     ll2_1d(:) = ll1_1d(:) .AND. (zicnc(:,jk) >= zicemin)

     ! ICNC burden (in-cloud) [m-2]
     ztmp1_1d(:) = zicnc_burden(:)+zicnc(:,jk)*zdz_2d(:,jk)
     zicnc_burden(:) = MERGE(ztmp1_1d(:), zicnc_burden(:), ll2_1d(:))
     ICNC_burden_acc(1:kproma,jrow) = zicnc_burden(:)

     ! ICNC (grid-box) [m-3]
     icnc(1:kproma,jk,jrow) = &
          MERGE(zicnc(:,jk) * zclcaux(:), icnc(1:kproma,jk,jrow), ll2_1d(:))

     ! ICNC burden (grid-box) [m-2]
     icnc_burden(1:kproma,jrow) = &
          MERGE(icnc_burden(1:kproma,jrow) + zicnc(:,jk) * zdz_2d(:,jk) * &
          zclcaux(:), icnc_burden(1:kproma,jrow), ll2_1d(:))

     ! --- 8.3 - Corrections: avoid negative cloud water and cloud ice --------
     ! ------------------------------------------------------------------------

     ! Liquid clouds
     ll1_1d(:) = (zxlp1_1d(:) < ccwmin)  ! LWC below minimum LWC for cover>0
     zdxlcor_1d(:) = -ztmst_rcp * zxlp1_1d(:)
     zdxlcor_1d(:) = MERGE(zdxlcor_1d(:), 0._dp, ll1_1d(:))
     pxlte(:,jk) = pxlte(:,jk) + zdxlcor_1d(:)  ! correction

     ! CDNC tendency [1/mol/s]
     ztmp1_1d(:) = pxtte(:,jk,1) - &
          ztmst_rcp * m_air / 1000._dp * zcdnc(:,jk) * zrho_rcp(:,jk)
     pxtte(:,jk,1) = MERGE(ztmp1_1d(:), pxtte(:,jk,1), ll1_1d(:))

     ! Ice clouds
     ll2_1d(:) = (zxip1_1d(:) < ccwmin)  ! IWC below minimum IWC for cover>0
     zdxicor_1d(:) = -ztmst_rcp * zxip1_1d(:)
     zdxicor_1d(:) = MERGE(zdxicor_1d(:), 0._dp, ll2_1d(:))
     pxite(:,jk)   = pxite(:,jk) + zdxicor_1d(:)  ! correction

     ! ICNC tendency [1/mol/s]
     ztmp1_1d(:) = pxtte(:,jk,2) - &
          ztmst_rcp * m_air / 1000._dp * zicnc(:,jk) * zrho_rcp(:,jk)
     pxtte(:,jk,2) = MERGE(ztmp1_1d(:), pxtte(:,jk,2), ll2_1d(:))

     ! Cloud cover
     paclc(:,jk)   = MERGE(0.0_dp, paclc(:,jk), ll1_1d(:) .AND. ll2_1d(:))
     paclcac(:,jk) = paclcac(:,jk) + zdtime*paclc(:,jk)

     ! Specific humidity tendency [kg/kg/s]
     pqte(:,jk) = pqte(:,jk) - zdxlcor_1d(:) - zdxicor_1d(:)

     ! Temperature tendency [K/s]
     ptte(:,jk) = ptte(:,jk) + &
          zlvdcp(:) * zdxlcor_1d(:) + zlsdcp(:) * zdxicor_1d(:)

     ! --- 8.4 - Specific diagnostics to compare with in-situ data ------------
     ! ------------------------------------------------------------------------

     ! CDNC for LWC>0.01 g/m3 (DACCIWA field campaign)
     ztmp(:,jk) = 1000._dp * zxlb(:) * zrho(:, jk)  ! [g/m3]
     ll1_1d(:) = ztmp(:, jk) .gt. 0.01_dp
     cdnc_insitu(1:kproma, jk, jrow) = &
          MERGE(cdnc_acc(1:kproma, jk, jrow), 0._dp, ll1_1d(:))

     ! p_ice - Murphy and Koop (Q. J. R. Meteorol. Soc., 2005)
     ztmp(:,jk) = 100._dp * exp(24.7219_dp - 6024.5282_dp / ztp1tmp(:) +   &
          1.0613868e-2_dp * ztp1tmp(:) - 1.3198825e-5_dp * ztp1tmp(:)**2 - &
          0.49382577_dp * log(ztp1tmp(:)))

     ! q_ice
     ztmp1(:,jk) = 0.622_dp * ztmp(:,jk)/(papp1(:,jk) - 0.378_dp * ztmp(:,jk))

     ! q / q_ice
     ztmp2(:,jk) = zqp1tmp(:) / ztmp1(:,jk)


     ! Conditions for cirrus regime
     ll1_1d(:) = (ztp1tmp(:) .lt. 246.) .AND. (papp1(:,jk) .gt. 10000.)
     ll2_1d(:) = ll1_1d .AND. (zxib(:) > zeps)  ! for consistency with ICNC_acc
     ll3_1d(:) = ll1_1d .AND. (zxib(:).gt.csecfrl)  ! in-cloud criterion
     ll4_1d(:) = ll1_1d .AND. (zxib(:).le.csecfrl)  ! clearsky criterion

     DO jl = 1, kproma

        ! Temperature [K]
        IF (ztp1tmp(jl).le.(MINVAL(CIRRUS_TEMP) - 0.5_dp).OR. &
            ztp1tmp(jl).gt.(MAXVAL(CIRRUS_TEMP) + 0.5_dp)) CYCLE
        tbin = MAXLOC(CIRRUS_TEMP, DIM=1, &
                      MASK=(ztp1tmp(jl).gt.CIRRUS_TEMP - 0.5_dp))

        ! IWC in-cloud [ppmv]
        zdiag = MERGE(zxib(jl) * m_air / m_h2o * 1.e6_dp, 0._dp, ll2_1d(jl))
        IF (zdiag.gt.MINVAL(CIRRUS_IBIN_IWC) .AND. &
            zdiag.le.MAXVAL(CIRRUS_IBIN_IWC)) THEN
           ibin = MAXLOC(CIRRUS_IBIN_IWC, DIM=1, &
                         MASK=(zdiag.gt.CIRRUS_IBIN_IWC))
           CIRRUS_IWC(jl,tbin,ibin,jrow) = CIRRUS_IWC(jl,tbin,ibin,jrow) + 1
        ENDIF

        ! Rice [um] - Volume mean radius - Eq. (6) Lohmann et al. (2007)
        zdiag = MERGE(1.0e6_dp * (3._dp/(4._dp * api * zrhoice) *             &
                      MAX(zxib(jl),zeps)*zrho(jl,jk)/zicnc(jl,jk))**0.333_dp, &
                      0._dp, ll2_1d(jl))
        IF (zdiag.gt.MINVAL(CIRRUS_IBIN_Rice) .AND. &
            zdiag.le.MAXVAL(CIRRUS_IBIN_Rice)) THEN
           ibin = MAXLOC(CIRRUS_IBIN_Rice, DIM=1, &
                         MASK=(zdiag.gt.CIRRUS_IBIN_Rice))
           CIRRUS_Rice(jl,tbin,ibin,jrow) = CIRRUS_Rice(jl,tbin,ibin,jrow) + 1
        ENDIF

        ! Filter for the measured range of the NIXE-CAPS instrument
        ! 1.5 um < R < 480 um (3 um < D < 960 um), Kraemer priv. comm. (2019)
        ll5_1d(jl) = zdiag.ge.1.5_dp .AND. zdiag.le.480._dp

        ! Filter for the measured range of the ML-CIRRUS instruments
        ! 1.5 um < R < 3200 um (3 um < D < 6400 um), Heller priv. comm. (2019)
        ll6_1d(jl) = zdiag.ge.1.5_dp .AND. zdiag.le.3200._dp

        ! Nice in-cloud [cm-3] - Kraemer climatology
        zdiag = MERGE(zicnc(jl,jk)*1.e-6_dp, 0._dp, ll2_1d(jl).AND.ll5_1d(jl))
        IF (zdiag.gt.MINVAL(CIRRUS_IBIN_Nice) .AND. &
            zdiag.le.MAXVAL(CIRRUS_IBIN_Nice)) THEN
           ibin = MAXLOC(CIRRUS_IBIN_Nice, DIM=1, &
                         MASK=(zdiag.gt.CIRRUS_IBIN_Nice))
           CIRRUS_Nice(jl,tbin,ibin,jrow) = &
                CIRRUS_Nice(jl,tbin,ibin,jrow) + 1
        ENDIF

        ! Nice in-cloud [cm-3] - ML-CIRRUS
        zdiag = MERGE(zicnc(jl,jk)*1.e-6_dp, 0._dp, ll2_1d(jl).AND.ll6_1d(jl))
        IF (zdiag.gt.MINVAL(CIRRUS_IBIN_Nice) .AND. &
            zdiag.le.MAXVAL(CIRRUS_IBIN_Nice)) THEN
           ibin = MAXLOC(CIRRUS_IBIN_Nice, DIM=1, &
                         MASK=(zdiag.gt.CIRRUS_IBIN_Nice))
           CIRRUS_Nice_ML(jl,tbin,ibin,jrow) = &
                CIRRUS_Nice_ML(jl,tbin,ibin,jrow) + 1
        ENDIF


        ! RHice (in-cloud) [%]
        zdiag = MERGE(100._dp * ztmp2(jl,jk), 0._dp, ll3_1d(jl))
        IF (zdiag.gt.MINVAL(CIRRUS_IBIN_RHi) .AND. &
            zdiag.le.MAXVAL(CIRRUS_IBIN_RHi)) THEN
           ibin = MAXLOC(CIRRUS_IBIN_RHi, DIM=1, &
                         MASK=(zdiag.gt.CIRRUS_IBIN_RHi))
           CIRRUS_RHi_cloud(jl,tbin,ibin,jrow) = &
                CIRRUS_RHi_cloud(jl,tbin,ibin,jrow) + 1
        ENDIF

        ! RHice (clear) [%]
        zdiag = MERGE(100._dp * ztmp2(jl,jk), 0._dp, ll4_1d(jl))
        IF (zdiag.gt.MINVAL(CIRRUS_IBIN_RHi) .AND. &
            zdiag.le.MAXVAL(CIRRUS_IBIN_RHi)) THEN
           ibin = MAXLOC(CIRRUS_IBIN_RHi, DIM=1, &
                         MASK=(zdiag.gt.CIRRUS_IBIN_RHi))
           CIRRUS_RHi_clear(jl,tbin,ibin,jrow) = &
                CIRRUS_RHi_clear(jl,tbin,ibin,jrow) + 1
        ENDIF

        ! Vertical velocity [cm/s]
        zdiag = zvervx_2d(jl,jk)
        IF (zdiag.gt.MINVAL(CIRRUS_IBIN_VEL) .AND. &
            zdiag.le.MAXVAL(CIRRUS_IBIN_VEL)) THEN
           ibin = MAXLOC(CIRRUS_IBIN_VEL, DIM=1, &
                         MASK=(zdiag.gt.CIRRUS_IBIN_VEL))
           CIRRUS_vervel(jl,tbin,ibin,jrow) = &
                CIRRUS_vervel(jl,tbin,ibin,jrow) + 1
        ENDIF

     ENDDO

831 ENDDO  ! Vertical loop


  ! === 9 - WET CHEMISTRY AND IN-CLOUD SCAVENGING =============================
  ! ===========================================================================

  ! Not needed in EMAC, since it is done in SCAV


  ! === 10 - DIAGNOSTICS ======================================================
  ! ===========================================================================

  ! --- 10.1 - Accumulated precipitation at the surface -----------------------
  ! ---------------------------------------------------------------------------

   prsfl(:) = zrfl(:)
   pssfl(:) = zsfl(:)
   paprl(:) = paprl(:) + zdtime * (prsfl(:) + pssfl(:))
   paprs(:) = paprs(:) + zdtime * pssfl(:)


   ! --- 10.2 - Total cloud cover ---------------------------------------------
   ! --------------------------------------------------------------------------

   zclcov(:) = 1.0_dp - paclc(:,1)

   DO 923 jk = 2,klev
      ztmp1_1d(:) = MAX(paclc(:,jk), paclc(:,jk-1))
      ztmp2_1d(:) = MIN(paclc(:,jk-1), zxsec)

      zclcov(:) = zclcov(:) * (1._dp - ztmp1_1d(:)) / (1._dp - ztmp2_1d(:))
923 ENDDO

   zclcov(:)  = 1.0_dp - zclcov(:)
   paclcov(:) = paclcov(:) + zdtime * zclcov(:)


   ! --- 10.3 - Vertical integrals of humidity, cloud water and cloud ice -----
   ! --------------------------------------------------------------------------

   zqvi(:)  = 0.0_dp
   zxlvi(:) = 0.0_dp
   zxivi(:) = 0.0_dp

   DO 933 jk = ktdia,klev
      zqvi(:)  = zqvi(:)  + pqm1(:,jk) *zdpg_2d(:,jk)
      zxlvi(:) = zxlvi(:) + pxlm1(:,jk)*zdpg_2d(:,jk)
      zxivi(:) = zxivi(:) + pxim1(:,jk)*zdpg_2d(:,jk)
933 ENDDO

   pqvi(:)  = pqvi(:)  + zdtime*zqvi(:)
   pxlvi(:) = pxlvi(:) + zdtime*zxlvi(:)
   pxivi(:) = pxivi(:) + zdtime*zxivi(:)

   RETURN

  END SUBROUTINE cloud_cdnc_icnc3


  ! ===========================================================================

  SUBROUTINE get_util_var(kproma, kbdim, ktdia, klev, klevp1, paphm1, pgeo, &
                          papm1, ptm1, pgeoh, pdp, pdpg, pdz, paaa, pviscos)

    ! Get several utility variables:
    !   geopotential at half levels (pgeoh)
    !   pressure- and height-differences (pdp, pdz)
    !   air density correction for computing ice crystal fall velocity (paaa)
    !   dynamic viscosity of air (pviscos)

    INTEGER, INTENT(IN) :: kproma, kbdim, ktdia, klev, klevp1

    REAL(dp), INTENT(IN)  :: paphm1(kbdim,klevp1), pgeo(kbdim,klev), &
         papm1(kbdim,klev), ptm1(kbdim,klev)
    REAL(dp), INTENT(OUT) :: pgeoh(kbdim,klevp1),  pdp(kbdim,klev), &
         pdpg(kbdim,klev), pdz(kbdim,klev), paaa(kbdim,klev), &
         pviscos(kbdim,klev)

    REAL(dp) :: g_rcp

    g_rcp = 1._dp / g

    ! Geopotential at half levels:
    pgeoh(:,ktdia+1:klev) = &
         0.5_dp * (pgeo(:,ktdia+1:klev) + pgeo(:,ktdia:klev-1))
    pgeoh(:,ktdia) = pgeo(:,ktdia) + (pgeo(:,ktdia)-pgeoh(:,ktdia+1))
    pgeoh(:,klevp1) = 0.0_dp

    ! Pressure differences [Pa]
    pdp(:,ktdia:klev) = paphm1(:,ktdia+1:klevp1) - paphm1(:,ktdia:klev)
    pdpg(:,ktdia:klev) = g_rcp * pdp(:,ktdia:klev)

    ! Height differences [m]
    pdz(:,ktdia:klev) = g_rcp * (pgeoh(:,ktdia:klev) - pgeoh(:,ktdia+1:klevp1))

    ! Air density correction
    ! Spichtinger and Gierens (Atmos. Chem. Phys., 2009) - Eq. (19)
    paaa(:,:) = ((papm1(:,:) / 30000._dp)**(-0.178_dp) ) * &
         ((ptm1(:,:) / 233.0_dp)**(-0.394_dp))

    ! Dynamic viscosity of air
    pviscos(:,:) = (1.512_dp + 0.0052_dp * (ptm1(:,:) - 233.15_dp)) * 1.e-5_dp

  END SUBROUTINE get_util_var

  ! ===========================================================================

  SUBROUTINE get_cloud_bounds(kproma, kbdim, ktdia, klev, paclc, ktop, kbas, &
                              kcl_minustop, kcl_minusbas)

    ! Flag top, base and intermediate levels for each cloud

    INTEGER, INTENT(IN)  :: kproma, kbdim, ktdia, klev
    INTEGER, INTENT(OUT) :: ktop(kbdim,klev), &         ! cloud top
                            kbas(kbdim,klev), &         ! cloud base
                            kcl_minustop(kbdim,klev), & ! all levs except top
                            kcl_minusbas(kbdim,klev)    ! all levs except base

    REAL(dp), INTENT(IN) :: paclc(kbdim,klev)  ! cloud cover

    INTEGER  :: jk, jl, jm, jnumb, jtop, jbas
    INTEGER  :: iindex(kbdim,klev),    &  ! level index
                 iclnb(kbdim),         &  ! number of clouds per column
                 iclbounds(2,klev/2+1)    ! boundary infos per cloud

    REAL(dp), DIMENSION(kbdim,klev) :: zaclcm, zaclcp

    LOGICAL, DIMENSION(kbdim, klev) :: ll, llm, llp, lltop, llbas

    ! Initialization
    zaclcm(:,:)       = 0._dp
    zaclcp(:,:)       = 0._dp
    iclnb(:)          = 0
    kcl_minustop(:,:) = 0
    kcl_minusbas(:,:) = 0

    ll(:,:)    = .FALSE.
    llm(:,:)   = .FALSE.
    llp(:,:)   = .FALSE.
    lltop(:,:) = .FALSE.
    llbas(:,:) = .FALSE.

    DO jk = ktdia, klev
       iindex(:,jk) = jk
    ENDDO

    ! Duplicate paclc at level-1 and level+1
    zaclcm(:,ktdia+1:klev) = paclc(:,ktdia:klev-1)
    zaclcp(:,ktdia:klev-1) = paclc(:,ktdia+1:klev)

    ! Set logical switches
    ll(:,:)  = (paclc(:,:)  >= zepsec)
    llm(:,:) = (zaclcm(:,:) <  zepsec)
    llp(:,:) = (zaclcp(:,:) <  zepsec)

    lltop(:,:) = (ll(:,:) .AND. llm(:,:)) ! true if cloud top detected
    llbas(:,:) = (ll(:,:) .AND. llp(:,:)) ! true if cloud base detected

    ! Set itop and ibas
    ktop(:,:) = MERGE(iindex(:,:), 0, lltop(:,:))
    kbas(:,:) = MERGE(iindex(:,:), 0, llbas(:,:))

    ! Reset the logical switches
    lltop(:,:) = .FALSE.
    llbas(:,:) = .FALSE.

    lltop(:,:) = (ktop(:,:) > 0)
    llbas(:,:) = (kbas(:,:) > 0)

    ! Counts the number of clouds per column
    iclnb(:) = COUNT(lltop(:,:),2)

    DO jl=1,kproma
       jnumb = iclnb(jl)

       ! Sets the bounds in a compact array
       iclbounds(:,:) = 0

       iclbounds(1,1:jnumb) = PACK(ktop(jl,:), lltop(jl,:))  ! cloud tops
       iclbounds(2,1:jnumb) = PACK(kbas(jl,:), llbas(jl,:))  ! cloud bases

       ! Flag cloud levels excepted their base (or top)
       DO jm = 1, jnumb
          jtop = iclbounds(1,jm)
          jbas = iclbounds(2,jm)
          kcl_minusbas(jl,jtop:jbas-1) = jbas
          kcl_minustop(jl,jtop+1:jbas) = jtop
       ENDDO
    ENDDO

  END SUBROUTINE get_cloud_bounds

  ! ===========================================================================
  SUBROUTINE calc_vervel_orogw(kproma, klev, l_z, tt, pp, qq, icesat, ampl, &
       wgwd, imask, vervel)

    USE messy_main_constants_mem, ONLY: g, cp_air

    INTEGER, INTENT(IN) :: kproma, klev
    REAL(dp), DIMENSION(kproma), INTENT(IN) :: l_z, imask
    REAL(dp), DIMENSION(kproma,klev), INTENT(IN) :: tt, pp, qq, icesat, ampl, &
         wgwd
    REAL(dp), DIMENSION(kproma,klev), INTENT(INOUT) :: vervel
    LOGICAL, DIMENSION(kproma) :: ll1, ll2
    INTEGER :: jl, jk
    REAL(dp) :: lambda, Scr0, Scr, rhival, i, tempgw, wsin, pice, rhi

    vervel = 0._dp

    ! Reduction of max vertical velocity (from boxmodel calculations)
    DO jk = klev, 1, -1  ! vertical loop
       DO jl = 1, kproma

          ! Wavelength
          lambda = 2._dp * l_z(jl)

          ! Critical supersaturation for homogeneous freezing
          ! Kaercher and Lohmann (J. Geophys. Res., 2002a) - Eq (4)
          Scr0 = (2.583_dp - tt(jl,jk) / 207.83_dp) * 100._dp

          ! Relative humidity w.r.t. ice
          rhival = icesat(jl,jk) * 100._dp + 100._dp
          i = 1._dp

          IF (rhival .gt. Scr0) THEN
             vervel(jl,jk) = 0.0_dp
          ELSE
             DO WHILE (rhival.le.Scr0 .and. i.lt.lambda)
                tempgw = MAX(100._dp, &
                             tt(jl,jk) - ((g/cp_air*sqrt(ampl(jl,jk))) * &
                             SIN(3._dp*api/2._dp + 2._dp*api*i/lambda) + &
                             (g/cp_air*SQRT(ampl(jl,jk)))))
                wsin = wgwd(jl,jk) * SIN(2._dp * api * i / lambda)

                ! Murphy and Koop (Q. J. R. Meteorol. Soc., 2005) - Eq. (7)
                pice = EXP(9.550426_dp - (5723.265_dp / tempgw) + &
                     (3.53068_dp * LOG(tempgw)) - (0.00728332_dp * tempgw))
                rhi = (100._dp * pp(jl,jk) * qq(jl,jk)) / &
                     (0.622_dp * pice)

                ! Critical supersaturation for homogeneous freezing
                ! Koop et al. (Nature, 2000)
                Scr = (2.418_dp - tempgw / 245.68_dp) * 100._dp
                rhival = rhi
                Scr0 = Scr
                vervel(jl,jk) = wsin ! [m/s]
                i = i + 1000._dp
             ENDDO
          ENDIF
       ENDDO

       ll1(:) = (l_z(:).lt.100._dp)
       vervel(:,jk) = MERGE(0._dp, vervel(:,jk), ll1(:))

       ll2(:) = (imask(:).eq.1) ! mask for orographic waves
       vervel(:,jk) = MERGE(vervel(:,jk), 0._dp, ll2(:))

       vervel(:,jk) = MAX(0.001_dp, vervel(:,jk)) ! [m/s]
       vervel(:,jk) = vervel(:,jk) * 100.0_dp ! [cm/s]

    ENDDO

  ENDSUBROUTINE calc_vervel_orogw

  ! ===========================================================================


  ! ===========================================================================
  ! === KAERCHER ET AL. (2006) PARAMETRIZATION FOR CIRRUS CLOUDS ==============
  ! ===========================================================================
  !
  ! Structure:
  !   - FUNCTION SCRHOM: calculate the critical ice saturation ratio for
  !                      homogeneous freezing
  !   - FUNCTION PISAT: calculate vapour pressure over ice [hPa]
  !   - FUNCTION XEERFC: calculate the product EXP(1/Y) * ERFC(1/SQRT(Y))
  !   - FUNCTION XERF: error function with rational approximation
  !   - FUNCTION TAUG: dimensionless growth time scale
  !   - SUBROUTINE ACIDRV
  !   - SUBROUTINE XFRZMSTR
  !   - SUBROUTINE XICE
  !   - SUBROUTINE XFRZHET
  !   - SUBROUTINE XICEHET
  !   - SUBROUTINE XFRZHOM
  !   - SUBROUTINE XICEHOM
  !   - SUBROUTINE NUCTYPE
  !
  ! ===========================================================================

  REAL(dp) FUNCTION SCRHOM(T)

    ! ***** CRITICAL ICE SATURATION RATIO FOR HOMOGENEOUS FREEZING
    ! ***** EVALUATED FOR AN INTERMEDIATE PARTICLE SIZE OF 0.25 MUM
    ! ***** TH.KOOP, B.P.LUO, A.TSIAS, TH.PETER, NATURE 406, 611-614, 2000

    IMPLICIT NONE

    ! function parameters

    REAL(dp), INTENT(in) :: T

    ! local variables

    REAL(dp) :: TEMP
    TEMP = MIN(240.0_dp, MAX(170.0_dp, T))
    SCRHOM = 2.418_dp - (TEMP/245.68_dp)

    RETURN

  END FUNCTION SCRHOM


  ! ===========================================================================

  REAL(dp) FUNCTION PISAT(T)

    ! ***** VAPOR PRESSURE OVER ICE IN MBAR (T IN K)
    ! ***** J.MARTI AND K.MAUERSBERGER, GRL 20(5), 363-366, 1993

    IMPLICIT NONE

    ! function parameters

    REAL(dp), INTENT(in) :: T

    ! parameters

    REAL(dp), PARAMETER :: A = 0.01_dp
    REAL(dp), PARAMETER :: B = 12.537_dp
    REAL(dp), PARAMETER :: C = -2663.5_dp

    PISAT = A * 10.0_dp**( B + C/T )

    RETURN

  END FUNCTION PISAT


  ! ===========================================================================

  REAL(dp) FUNCTION XEERFC(Y)

    ! ***** PRODUCT EXP(1/Y) * ERFC(1/SQRT(Y)) AND ASYMPTOTES

    IMPLICIT NONE

    ! function parameters

    REAL(dp), INTENT(in) :: Y

    ! parameters

    REAL(dp), PARAMETER :: SQPI = 1.7724539_dp

    ! local variables

    REAL(dp) :: Y1, SQY1, PPROD, POLY
    INTEGER :: K

    Y1   = 1.0_dp / Y
    SQY1 = SQRT(Y1)
    IF (Y.LE.0.2_dp) THEN
       PPROD = 1.0_dp
       POLY  = 1.0_dp
       DO K = 1, 4
          PPROD = PPROD * REAL(2*K-1,dp)
          POLY  = POLY  + PPROD * (-0.5_dp*Y)**K
       ENDDO
       XEERFC = POLY / SQY1 / SQPI
    ELSEIF (Y.GE.2.0_dp) THEN
       PPROD  = 1.0_dp
       POLY   = 1.0_dp
       DO K = 1, 4
          PPROD = PPROD * REAL(2*K+1,dp)
          POLY  = POLY  + REAL(2**K,dp) * Y1**K / PPROD
       ENDDO
       XEERFC = EXP(Y1) - 2.0_dp*SQY1/SQPI * POLY
    ELSE
       XEERFC = ( 1.0_dp - XERF(SQY1) ) * EXP(Y1)
    ENDIF

    RETURN

  END FUNCTION XEERFC


  ! ===========================================================================

  REAL(dp) FUNCTION XERF(X)

    ! ***** ERROR FUNCTION WITH RATIONAL APPROXIMATION
    ! ***** M.ABRAMOWITZ, I.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS
    ! ***** DOVER, 10th PRINTING 1972, p.299, # 7.1.26

    IMPLICIT NONE

    ! function parameters

    REAL(dp), INTENT(in) :: X

    ! parameters

    REAL(dp), PARAMETER :: P  = 0.3275911_dp
    REAL(dp), PARAMETER :: A1 = 0.254829592_dp
    REAL(dp), PARAMETER :: A2 =-0.284496736_dp
    REAL(dp), PARAMETER :: A3 = 1.421413741_dp
    REAL(dp), PARAMETER :: A4 =-1.453152027_dp
    REAL(dp), PARAMETER :: A5 = 1.061405429_dp

    ! local variables

    REAL(dp) :: T, ARG

    T    = 1.0_dp / ( 1.0_dp + P * ABS(X) )
    ARG  = MIN( X**2, 75.0_dp )
    XERF = 1.0_dp - (A1*T + A2*T**2 + A3*T**3 + A4*T**4 + A5*T**5) * EXP(-ARG)
    XERF = SIGN( XERF, X )

    RETURN

  END FUNCTION XERF


  ! ===========================================================================

  REAL(dp) FUNCTION TAUG(B,Y,X0)

    ! ***** DIMENSIONLESS GROWTH TIME SCALE
    ! ***** B.KÄRCHER AND S.SOLOMON, JGR 104(D22), 27441-27459, 1999

    IMPLICIT NONE

    REAL(dp), INTENT(in) :: B, Y, X0

    REAL(dp), PARAMETER :: SQ31 = 0.577350269_dp
    REAL(dp), PARAMETER :: SIX1 = 0.166666666_dp

    REAL(dp) :: X, F1, F10, F2, F20

    TAUG = 0.0_dp

    IF (Y.LE.X0) RETURN

    X    = MIN(Y, 0.999_dp )
    F1   = SIX1 * LOG( (1.0_dp+X +X**2.0_dp) / (1.0_dp-X)**2.0_dp )
    F10  = SIX1 * LOG( (1.0_dp+X0+X0**2._dp) / (1.0_dp-X0)**2.0_dp )
    F2   = SQ31 * ATAN( SQ31*(1.0_dp+2.0_dp*X) )
    F20  = SQ31 * ATAN( SQ31*(1.0_dp+2.0_dp*X0) )
    TAUG = (B+1.0_dp)*(F1-F10) + (B-1.0_dp)*(F2-F20)

    RETURN

  END FUNCTION TAUG


  ! ===========================================================================

  SUBROUTINE ACIDRV(JK,        &  ! current level
                    KBDIM,     &  ! kproma
                    KPROMA,    &  ! kproma
                    KLEV,      &  ! number of levels
                    NEXT,      &  ! =0 no pre-existing, =1 pre-existing IC
                    NAER,      &  ! number of freezing aerosol particles
                    NFRZMOD,   &  ! largest number of log-normal modes per type
                    ZVERVX,    &  ! constant updraft velocity [cm/s]
                    ZTMST,     &  ! time-step of the calling model [s]
                    PAPM1,     &  ! initial air pressure [hPa]
                    PTM1,      &  ! initial air temperature [K]
                    ZSUSATIX,  &  ! initial ice supersaturation
                    ZAPNX,     &  ! number density of aerosol modes [cm-3]
                    ZAPRX,     &  ! mean radius of aerosol modes [cm]
                    ZAPSIGX,   &  ! standard deviation of aerosol modes
                    SCRHET,    &  ! critical saturation ratio for het. freezing
                    CAER_2d,   &  ! total aerosol number density [cm-3]
                    RAER_2d,   &  ! radius of the smallest freezing aerosol [cm]
                    ZNICEX_3d, &  ! number density of IC [m-3]
                    ZRI_3d,    &  ! mean radius of IC [m]
                    ZVICE,     &  ! fictitious downdraft [cm/s]
                    CTHOMI,    &  ! homogeneous freezing temperature [K]
                    ZSICIRRUS, &  ! highest exceeded freezing threshold
                    STATUS,    &  ! return status
                    ERROR)        ! error message (if status /= 0)

    ! ***** DRIVER TO DEMONSTRATE USE OF CIRRUS PARAMETERIZATION

    IMPLICIT NONE

    ! SUBROUTINE parameters

    ! INPUT
    INTEGER ::  JL
    INTEGER ::  JK      ! CURRENT LEVEL
    INTEGER ::  KBDIM   ! SIZE OF ARRAYS
    INTEGER ::  KLEV    ! NUMBER OF LEVELS number of levels
    INTEGER ::  KPROMA  ! NUMBER OF CELLS IN ARRAY
    INTEGER ::  NEXT    ! TOTAL NUMBER OF ICE TYPES ORIGINATING FROM SOURCES
                        ! OTHER THAN AEROSOL PARTICLES (MIXED-PHASE CLOUD,
                        ! CONTRAIL CIRRUS, CONVECTIVE ICE)
    INTEGER ::  NAER    ! TOTAL NUMBER OF FREEZING AEROSOL PARTICLE
                        ! TYPES WITH UNLIKE FREEZING THRESHOLDS,
                        ! STARTING WITH THE MOST EFFICIENT IN AND
                        ! ENDING WITH HOMOGENEOUS NUCLEI
    INTEGER ::  NFRZMOD ! LARGEST NUMBER OF LOG-NORMAL MODES PER NAER-TYPE

    REAL(dp) :: ZVERVX(kbdim) ! CONSTANT UPDRAFT SPEED [CM/S]
    REAL(dp) :: ZTMST         ! TIME STEP of CALLING MODEL [S]
    REAL(dp) :: DTGM          ! TIME STEP of the PARAMETERIZATION AFTER WHICH
                              ! THE SOLUTION IS RETURNED [S]
    REAL(dp) :: PAPM1(kbdim, klev)  ! INITIAL AIR PRESSURE [hPa]
    REAL(dp) :: PTM1(kbdim, klev)   ! INITIAL AIR TEMPERATURE [K]
    REAL(dp) :: ZSUSATIX(kbdim)     ! INITIAL ICE SUPERSATURATION [1]
    REAL(dp) :: PAIR
    REAL(dp) :: TAIR
    REAL(dp) :: ZAPNX(kbdim,1:NAER,1:NFRZMOD)   ! NUMBER DENSITY OF AEROSOL
                                                ! MODE PER TYPE [CM-3]
    REAL(dp) :: ZAPRX(kbdim,1:NAER,1:NFRZMOD)   ! MEAN RADIUS OF AEROSOL MODE
                                                ! PER TYPE [CM]
    REAL(dp) :: ZAPSIGX(kbdim,1:NAER,1:NFRZMOD) ! STANDARD DEVIATION OF AEROSOL
                                                ! MODE PER TYPE [1]
    REAL(dp) :: SCRHET(2,1:NAER-1)              ! CRITICAL SATURATION RATIO FOR
                                                ! HETEROGENEOUS FREEZING [1]
    REAL(dp) :: ZNICEX_3d(kbdim,klev,1:(NEXT+NAER)) ! INITIAL NUMBER DENSITY
                                                    ! OF IC [M-3]
    REAL(dp) :: ZRI_3d(kbdim,klev,1:(NEXT+NAER))  ! INITIAL MEAN RADIUS OF
                                                  ! PREEXISTING IC [M]
    REAL(dp) :: CAER_2d(kbdim,klev,1:NAER)        ! TOTAL AEROSOL NUMBER
                                                  ! DENSITY PER TYPE [CM-3]
    REAL(dp) :: RAER_2d(kbdim,klev, 1:NAER)       ! RADIUS OF SMALLEST FREEZING
                                                  ! PARTICLE PER TYPE [CM]
    REAL(dp) :: ZVICE(kbdim)        ! FICTITIOUS DOWNDRAFT
    REAL(dp), INTENT(in) :: CTHOMI  ! TEMPERATURE OF HOMOGENEOUS FREEZING [K]

    ! OUTPUT
    INTEGER :: STATUS      ! STATUS
    CHARACTER(70) :: ERROR ! ERROR MESSAGE
    REAL(dp), INTENT(out):: ZSICIRRUS(kbdim) ! HIGHEST EXCEEDED FREEZING
                                             ! THRESHOLD [1]
    REAL(dp) :: SSTOP_CIRRUS ! HIGHEST FREEZING THRESHOLD OF NUCLEATION EVENT
    REAL(dp) :: CAER(1:NAER) ! TOTAL AEROSOL NUMBER DENSITY PER TYPE [CM-3]
    REAL(dp) :: RAER(1:NAER) ! RADIUS OF SMALLEST FREEZING PART. PER TYPE [CM±
    REAL(dp) :: CICE(1:(NEXT+NAER)) ! NUMBER DENS. OF IC PER AEROSOL TYPE [CM-3]
    REAL(dp) :: RICE(1:(NEXT+NAER)) ! MEAN RADIUS OF IC PER AEROSOL TYPE [CM]
    REAL(dp) :: SICE  ! FINAL ICE SATURATION RATIO [1]
    REAL(dp) :: V     ! CONSTANT UPDRAFT SPEED [CM/S]
    REAL(dp) :: P     ! INITIAL AIR PRESSURE [hPa]
    REAL(dp) :: T     ! INITIAL AIR TEMPERATURE [K]
    REAL(dp) :: S     ! INITIAL ICE SATURATION RATIO [-]
    REAL(dp) :: C(1:NAER, 1:NFRZMOD)
    REAL(dp) :: R(1:NAER, 1:NFRZMOD)
    REAL(dp) :: SIG(1:NAER, 1:NFRZMOD)


    ! local parameters
    REAL(dp) :: DZMAX  ! maximum vertical distance an air parcel is allowed
                       ! to rise per timestep [CM] (= 1000000, 10 km)
    REAL(dp) :: DTMAX  !
    REAL(dp) :: PIN    ! initial air pressure [hPa]
    REAL(dp) :: TIN    ! initial air temperature [K]
    REAL(dp) :: fad
    INTEGER ::  NTYPE  ! total number of ice particle types
    INTEGER ::  NT     ! index running from 1 to NTYPE

    !--------------------------------------------------------------------------
    ! HETEROGENEOUS ICE NUCLEATION PROPERTIES ARE DEFINED IN THE SCRHET INPUT
    ! ARGUMENT IN CONJUNCTION WITH SUBROUTINE NUCTYPE AND CONSISTENT WITH THE
    ! CHOICES FOR NAER AND NEXT
    !
    ! IN GCM APPLICATIONS, LIMIT DTGM SUCH THAT PARCEL DOES NOT RISE MORE THAN
    ! ~1 KM AND LIMIT V SUCH THAT IT IS LARGER THAN VICE ON INPUT IN XICE  AND
    ! IS LARGER THAN ~0.1 MM/S. THESE MEASURES HELP PREVENT PROBLEMS WITH
    ! INFINITE TIME LOOPS IN  XICE (see DZMAX above)
    !--------------------------------------------------------------------------

    DO JL = 1, kproma

       DTGM = ZTMST

       IF (zvervx(JL).le.0.01_dp   .OR. &
           zvervx(JL).le.zvice(JL) .OR. &
           zsusatix(JL).le.0.0_dp  .OR. &
           ptm1(JL,JK).ge.cthomi) THEN

          ZSICIRRUS(JL) = ZSUSATIX(JL) + 1  ! set actual S as minimum S
          ZSICIRRUS(JL) = MAX(ZSICIRRUS(JL), MINVAL(SCRHET))
          CYCLE
       ENDIF

       NTYPE = NEXT + NAER

       ! Pass EMAC variable to local copies
       RICE(:)  = ZRI_3d(JL,JK,:)
       CICE(:)  = ZNICEX_3d(JL,JK,:)
       S        = ZSUSATIX(JL)
       V        = ZVERVX(JL)
       P        = PAPM1(JL,JK)
       T        = MAX(PTM1(JL,JK), 150._dp)  ! Kaercher, priv. comm. (2019)
       CAER(:)  = CAER_2d(JL,JK,:)
       RAER(:)  = RAER_2d(JL,JK,:)
       C(:,:)   = ZAPNX(JL,:,:)
       R(:,:)   = ZAPRX(JL,:,:)
       SIG(:,:) = ZAPSIGX(JL,:,:)

       ! Convert MKS (EMAC) to CGS (parametrization) units
       P = P * 0.01_dp  ! [hPa]
       S = S + 1._dp    ! supersaturation as saturation ratio
       CICE(:) = CICE(:) * 1.E-6_dp  ! [cm-3]
       RICE(:) = RICE(:) * 100._dp   ! [cm]

       ! Set actual S as minimum, will be overwritten in XICE
       SSTOP_CIRRUS = S

       ! Limit DTGM such that air parcel does not rise more than 1 km per tstep
       DZMAX = 1.E5_dp
       DTMAX = DZMAX / V
       DTGM = MIN(DTGM, DTMAX)
       DTGM = FLOAT(INT(DTGM + 0.5_dp))
       PIN = P
       TIN = T

       CALL XFRZMSTR(JL, JK, NEXT, NAER, NFRZMOD, DTGM, P, T, V, S, C, R, &
                     SIG, SCRHET, CAER, RAER, CICE, RICE, SICE, PAIR, TAIR, &
                     SSTOP_CIRRUS, STATUS, ERROR)
       IF (STATUS /= 0) RETURN

       fad = (TAIR / TIN)**2.5_dp

       DO NT = 1, NTYPE
          CICE(NT)=CICE(NT) / fad
          IF (NT.GT.NEXT) THEN
             CAER(NT - NEXT) = CAER(NT - NEXT) / fad
          ENDIF
       ENDDO

       ! Convert back to MKS (EMAC) units
       ZNICEX_3d(JL,JK,:) = CICE(:) * 1.0E6_dp  ! [cm-3]
       ZRI_3d(JL,JK,:)    = RICE(:) * 0.01_dp   ! [m]  !MK: in m

       ZSICIRRUS(JL)   = SSTOP_CIRRUS
       ZSICIRRUS(JL)   = MAX(ZSICIRRUS(JL), MINVAL(SCRHET))

    ENDDO

    STATUS = 0
    RETURN

  END SUBROUTINE ACIDRV


  ! ===========================================================================

  SUBROUTINE XFRZMSTR(jl, jk, NEXT, NAER, NFRZMOD, DTGM, P, T, V, S, C, R, &
                      SIG, SCRHET, CAER, RAER, CICE, RICE, SICE, PAIR, TAIR, &
                      SSTOP_CIRRUS, STATUS, ERROR)

    ! -------------------------------------------------------------------------
    ! ***** BERND KAERCHER APR 05 2006
    ! ***** bernd.kaercher@dlr.de  http://www.pa.op.dlr.de/~pa1c/
    !
    ! ***** B. Kaercher, J. Hendricks, and U. Lohmann
    ! ***** Physically-based parameterization of cirrus cloud formation for use
    ! ***** in global atmospheric models
    ! ***** J. Geophys. Res., 111, D01205, doi:10.1029/2005JD006219, 2006.
    !
    ! ***** CHANGES IN ROUTINE XICE AFTER PUBLICATION:
    !       1/ NO FURTHER COOLING AFTER FINAL NUCLEATION EVENT
    !       2/ VENTILATION CORRECTION FOR SPHERICAL CRYSTALS CONSISTENT WITH
    !          ECHAM4 MICROPHYSICS
    !       3/ EXPLCIT TERMINATION OF INTEGRATION AFTER THE GLOBAL MODEL TIME
    !          STEP HAS BEEN REACHED OR A STEADY-STATE DEVELOPS IN SICE(TIME)
    !--------------------------------------------------------------------------

    IMPLICIT NONE

    ! SUBROUTINE parameters
    INTEGER :: jl, jk
    INTEGER :: NEXT
    INTEGER :: NAER
    INTEGER :: NFRZMOD

    REAL(dp) :: C(1:NAER,1:NFRZMOD)
    REAL(dp) :: R(1:NAER,1:NFRZMOD)
    REAL(dp) :: SIG(1:NAER,1:NFRZMOD)
    REAL(dp) :: SCRHET(2,1:NAER-1)
    REAL(dp) :: CAER(1:NAER)
    REAL(dp) :: RAER(1:NAER)
    REAL(dp) :: CICE(1:(NEXT+NAER))
    REAL(dp) :: RICE(1:(NEXT+NAER))
    REAL(dp) :: DTGM
    REAL(dp) :: P
    REAL(dp) :: T
    REAL(dp) :: V
    REAL(dp) :: S
    REAL(dp) :: SICE
    REAL(dp) :: PAIR
    REAL(dp) :: TAIR

    ! local parameters
    LOGICAL :: NOSIZE ! ignore for aerosol size effects on homogeneour freezing

    REAL(dp) :: CI(1:(NEXT+NAER))
    REAL(dp) :: RI(1:(NEXT+NAER))
    REAL(dp) :: CA(1:NFRZMOD)
    REAL(dp) :: RA(1:NFRZMOD)
    REAL(dp) :: SIGA(1:NFRZMOD)
    REAL(dp) :: DT
    REAL(dp) :: SI
    REAL(dp) :: S1
    REAL(dp) :: CVF
    REAL(dp) :: DTBACK
    REAL(dp) :: ADCHG
    REAL(dp) :: TIMER
    REAL(dp) :: PRESS
    REAL(dp) :: TEMP
    REAL(dp) :: TIME1
    REAL(dp) :: VICE
    REAL(dp) :: SINT
    REAL(dp) :: VEFF
    REAL(dp) :: CINT
    REAL(dp) :: RINT
    REAL(dp) :: RS
    REAL(dp) :: CASUM
    REAL(dp) :: SSTOP_CIRRUS

    INTEGER :: STATUS      ! STATUS
    CHARACTER(70) :: ERROR ! ERROR MESSAGE

    INTEGER :: NTYPE
    INTEGER :: NT
    INTEGER :: NM
    INTEGER :: INUC
    INTEGER :: NTT
    INTEGER :: NIN
    INTEGER :: STOP_MK

    REAL(dp) :: TIMINC
    COMMON /OUT/ TIMINC

    ! ***** INITIALIZATION
    NTYPE      = NEXT + NAER
    NOSIZE     = .FALSE.
    DT         = 0._dp
    SI         = S
    PAIR       = P
    TAIR       = T
    DO 10 NT = 1, NAER
       RAER(NT)  = 0._dp
       CAER(NT)  = 0._dp
       DO 10 NM  = 1, NFRZMOD
          CAER(NT) = CAER(NT) + C(NT,NM)
10  CONTINUE
    IF (V.LE.0._dp) THEN
       SICE      = S
       STATUS = 0
       RETURN
    ENDIF

    ! ***** CHECK WHETHER S EXCEEDS LOWEST IN THRESHOLD S1 ON INPUT
    ! ***** IF TRUE, STEP BACK IN TIME TO INCLUDE THE FIRST IN-TYPE


    IF (NAER.EQ.1) THEN
       S1 = SCRHOM(T)
    ELSE
       NT = 1
       CALL NUCTYPE(NT,NAER,INUC,CVF,STATUS,ERROR)
       IF (STATUS /= 0) RETURN
       IF (T .lt. 220._dp) THEN
          S1 = SCRHET(2, INUC)
       ELSE
          S1 = SCRHET(1, INUC)
       ENDIF
    ENDIF

    IF (SI.GE.S1) THEN
       DTBACK = 1.E5_dp * LOG(SI/(0.95_dp*S1)) / V
       SI         = 0.95_dp * S1
       DTGM       = DTGM + DTBACK
       TAIR       = T + 1.E-4_dp * V * DTBACK
       PAIR       = P * (TAIR/T)**3.5_dp
       ADCHG      = (TAIR/T)**2.5_dp
       DO NT      = 1, NTYPE
          CICE(NT)  = ADCHG * CICE(NT)
       ENDDO
       DO 20 NT   = 1, NAER
          CAER(NT)  = 0._dp
          DO 20 NM  = 1, NFRZMOD
             C(NT,NM) = ADCHG * C(NT,NM)
             CAER(NT) = CAER(NT) + C(NT,NM)
 20    CONTINUE
       GOTO 50
    ENDIF
50  CONTINUE

    TIMER       = DTGM
    TIMINC      = 0._dp

    ! ***** EVALUATE ICE SATURATION HISTORY
    DO 100 NT  = 1, NAER

       PRESS     = PAIR
       TEMP      = TAIR
       IF (NT.EQ.NAER) THEN
          S1       = 0._dp
       ELSE
          CALL NUCTYPE(NT,NAER,INUC,CVF,STATUS,ERROR)
          IF (STATUS /= 0) RETURN
          IF (T .lt. 220._dp) THEN
             S1 = SCRHET(2, INUC)
          ELSE
             S1 = SCRHET(1, INUC)
          ENDIF
       ENDIF
       DO NTT    = 1, NTYPE
          RI(NTT)  = RICE(NTT)
          CI(NTT)  = CICE(NTT)
       ENDDO
       TIME1     = -1._dp

       CALL XICE(NTYPE, DTGM, TIME1, TIMER, S1, V, PRESS, TEMP, RI, CI, SI, &
                 PAIR, TAIR, VICE, RICE, CICE, SICE, STOP_MK, SSTOP_CIRRUS, &
                 STATUS, ERROR)
       IF (STATUS /= 0) RETURN

       TIMER = TIMER - TIME1
       IF (STOP_MK == 1) THEN
          WRITE(*,*) 'SICE<1 in XICE; SICE=', SICE
          WRITE(*,*) 'Stop cirrus scheme and return'
          GOTO 200
       ENDIF
       IF (TIMER.LE.0._dp) GOTO 200
       TIMINC = TIMINC + TIME1
       ADCHG  = (TAIR / TEMP)**2.5_dp
       DO NM  = 1, NFRZMOD
          CA(NM)   = ADCHG * C(NT,NM)
          RA(NM)   = R(NT,NM)
          SIGA(NM) = SIG(NT,NM)
       ENDDO
       SINT = SICE
       VEFF = MAX(0._dp, (V-VICE))

       IF (NT.EQ.NAER) THEN
          CALL XFRZHOM(NOSIZE, NFRZMOD, DT, CA, RA, SIGA, PAIR, TAIR, VEFF, &
                       SINT, CINT, RINT, RS, STATUS, ERROR)
       ELSE
          CALL XFRZHET(NAER, NFRZMOD, INUC, CVF, DT, CA, RA, SIGA, SCRHET, &
                       PAIR, TAIR, VEFF, SINT, CINT, RINT, RS, T, STATUS, ERROR)
       ENDIF
       IF (STATUS /= 0) RETURN

       CASUM = 0._dp
       DO NM = 1, NFRZMOD
          CASUM    = CASUM + CA(NM)
       ENDDO
       CAER(NT) = MAX(0._dp,(CASUM-CINT))
       RAER(NT) = RS
       NIN = NEXT + NT
       RICE(NIN) = (RINT * CINT + RICE(NIN) * CICE(NIN)) / &
            (MAX(1.E-20_dp, (CINT + CICE(NIN))))
       CICE(NIN) = CINT + CICE(NIN)
       SI = SINT
100    CONTINUE
       PRESS = PAIR
       TEMP = TAIR
       DO NT = 1, NTYPE
          RI(NT) = RICE(NT)
          CI(NT) = CICE(NT)
       ENDDO
       TIME1 = TIMER
       CALL XICE(NTYPE, DTGM, TIME1, TIMER, S1, V, PRESS, TEMP, RI, CI, SI, &
                 PAIR, TAIR, VICE, RICE, CICE, SICE, STOP_MK, SSTOP_CIRRUS, &
                 STATUS, ERROR)
       IF (STATUS /= 0) RETURN

200 CONTINUE

    STATUS = 0
    RETURN

  END SUBROUTINE XFRZMSTR


  ! ===========================================================================

  SUBROUTINE XICE(NTYPE, DTGM, TIME1, TIMER, S1, V, P, T, R, C, S, &
                  PAIR, TAIR, VICE, RICE, CICE, SICE, STOP_MK, SSTOP_CIRRUS, &
                  STATUS, ERROR)

    ! ***** BERND KAERCHER  APR 05 2006
    ! ***** bernd.kaercher@dlr.de  http://www.op.dlr.de/~pa1c/

    IMPLICIT NONE

    ! SUBROUTINE parameters

    ! INPUT
    INTEGER :: NTYPE  ! TOTAL NUMBER OF ICE PARTICLE TYPES

    REAL(dp) :: R(1:NTYPE)  ! INITIAL MEAN  ICE CRYSTAL NUMBER RADIUS [CM]
    REAL(dp) :: C(1:NTYPE)  ! INITIAL TOTAL ICE CRYSTAL NUMBER DENSITY [CM-3]
    REAL(dp) :: DTGM        ! GLOBAL MODEL TIME STEP [s]
    REAL(dp) :: V           ! CONSTANT VERTICAL WIND SPEED > 0 [CM/S]
    REAL(dp) :: T           ! INITIAL AIR TEMPERATURE [K]
    REAL(dp) :: P           ! INITIAL AIR PRESSURE [hPa]
    REAL(dp) :: S           ! INITIAL ICE SATURATION RATIO [1]
    REAL(dp) :: TIME1       ! TIME AFTER WHICH SOLUTION IS RETURNED [S]
                            ! IF TIME1 > 0 ON INPUT AFTER THE LAST NUCLEATION
                            ! EVENT; SIMULTANEOUSLY, FURTHER COOLING IS
                            ! INHIBITED AND S1 IS UNUSED
    REAL(dp) :: TIMER       ! TIME AFTER WHICH SOLUTION IS RETURNED [S]
                            ! IF TIME1 < 0 ON INPUT
    REAL(dp) :: S1          ! IF TIME1 < 0 ON INPUT, THE SOLUTION IS RETURNED
                            ! WHEN SICE = S1, IN WHICH CASE TIME1 IS REPLACED
                            ! BY THE TIME NEEDED TO REACH S1 - THIS TIME IS
                            ! BOUNDED BY TIMER;
                            ! IF S1 = 0 ON INPUT, SCRHOM(T) IS USED FOR THE
                            ! CALCULATION

    ! OUTPUT
    REAL(dp) :: RICE(1:NTYPE) ! FINAL MEAN ICE CRYSTAL NUMBER RADIUS [CM]
    REAL(dp) :: CICE(1:NTYPE) ! FINAL MEAN ICE CRYSTAL NUMBER DENSITY [CM-3]
    REAL(dp) :: SICE          ! FINAL ICE SATURATION RATIO
    REAL(dp) :: TAIR          ! FINAL AIR TEMPERATURE
    REAL(dp) :: PAIR          ! FINAL AIR PRESSURE
    REAL(dp) :: VICE          ! EFFECTIVE VERTICAL VELOCITY CALCULATED
                              ! FROM DEPOSITIONAL LOSS TERM [CM/S]
    REAL(dp) :: SSTOP_CIRRUS
    INTEGER :: STATUS      ! STATUS
    CHARACTER(70) :: ERROR ! ERROR MESSAGE

    ! local parameters
    REAL(dp) :: DTIME
    REAL(dp) :: SSTOP
    REAL(dp) :: ALP4
    REAL(dp) :: T1
    REAL(dp) :: COOLR
    REAL(dp) :: TIME
    REAL(dp) :: SCHG
    REAL(dp) :: SCHGO
    REAL(dp) :: TEMP
    REAL(dp) :: TEMP1
    REAL(dp) :: PICE
    REAL(dp) :: ADCHG
    REAL(dp) :: PRESS
    REAL(dp) :: RHOA
    REAL(dp) :: VISC
    REAL(dp) :: FLUX
    REAL(dp) :: CISAT
    REAL(dp) :: A1, A2, A3
    REAL(dp) :: B1, B2
    REAL(dp) :: DLOSS
    REAL(dp) :: VFALL
    REAL(dp) :: REYN
    REAL(dp) :: VENT
    REAL(dp) :: DELP1
    REAL(dp) :: ETA
    REAL(dp) :: DSICE
    REAL(dp) :: YTIME
    REAL(dp) :: XTIME
    REAL(dp) :: STRY

    INTEGER :: IV
    INTEGER :: NT
    INTEGER :: M
    INTEGER :: M_MK
    INTEGER :: STOP_MK

    ! parameters
    REAL(dp), PARAMETER :: ALPHA = 0.5_dp
    REAL(dp), PARAMETER :: DTPLUS=1.01_dp
    REAL(dp), PARAMETER :: DTMINUS=0.1_dp
    REAL(dp), PARAMETER :: THOUBK = 7.24637701E+18_dp
    REAL(dp), PARAMETER :: FPIVOL = 3.89051704E+23_dp
    REAL(dp), PARAMETER :: FA1 = 0.601272523_dp
    REAL(dp), PARAMETER :: FA2 = 0.000342181855_dp
    REAL(dp), PARAMETER :: FA3 = 1.49236645E-12_dp
    REAL(dp), PARAMETER :: FC  = 9.80999976E-05_dp
    REAL(dp), PARAMETER :: FD = 249.239822_dp
    REAL(dp), PARAMETER :: FVTH = 11713803.0_dp
    REAL(dp), PARAMETER :: SVOL = 3.23E-23_dp
    REAL(dp), PARAMETER :: WVP1 = 3.6E+10_dp
    REAL(dp), PARAMETER :: WVP2 = 6145.0_dp
    REAL(dp), PARAMETER :: TEPS = 1.E-3_dp
    REAL(dp), PARAMETER :: SEPS = 1.E-6_dp
    REAL(dp), PARAMETER :: SERR = 1.E-3_dp

    REAL(dp) :: TIMINC
    COMMON /OUT/ TIMINC


    ! ***** CHECK INPUT
    IF (NTYPE.LT.1) THEN
       ERROR = 'K14: invalid NTYPE value in XICE'
       STATUS = 1
       RETURN
    ENDIF
    DO NT = 1, NTYPE
       IF (R(NT).LT.0._dp) THEN
          ERROR = 'K14: invalid R value in XICE'
          STATUS = 1
          RETURN
       ENDIF
       IF (C(NT).LT.0._dp) THEN
          ERROR = 'K14: invalid C value in XICE'
          STATUS = 1
          RETURN
       ENDIF
    ENDDO
    IF (S.LE.0._dp) THEN
       ERROR = 'K14: invalid S value in XICE'
       STATUS = 1
       RETURN
    ENDIF
    IF (TIME1.LT.0._dp) THEN
       IF (S1.NE.0._dp .AND. S.GT.S1) THEN
          ERROR = 'K14: S exceeds S1 in XICE'
          STATUS = 1
          RETURN
       ENDIF
       IF (V.LT.0._dp) THEN
          ERROR = 'K14: invalid V value in XICE'
          STATUS = 1
          RETURN
       ENDIF
    ENDIF

    ! ***** DO NOTHING
    IF (((TIME1.GE.0._dp).AND.(TIME1.LT.TEPS))) THEN
       TAIR  = T
       PAIR  = P
       VICE  = 0._dp
       DO NT = 1, NTYPE
          RICE(NT) = R(NT)
          CICE(NT) = C(NT)
       ENDDO
       SICE      = S
       STATUS = 0
       RETURN
    ENDIF

    ! ***** SCALE TIME STEP WITH UPDRAFT SPEED AND INHIBIT COOLING
    ! ***** AFTER LAST NUCLEATION EVENT (THEN TIME1 > 0 ON INPUT)
    IV = 1
    DTIME  = 0.0005_dp * (1.E5_dp/V)
    IF (TIME1.GE.0._dp) THEN
       IV = 0
       DTIME = MIN(DTIME,TIME1)
    ENDIF

    ! ***** SET CONSTANTS
    SSTOP = S1
    ALP4  = 0.25_dp * ALPHA
    T1    = 1._dp/ T
    COOLR = FC * V*REAL(IV)

    ! ***** SET INITIAL VALUES
    M = 0
    M_MK = 0
    STOP_MK = 0
    TIME = 0._dp
    SCHG = 100._dp
    DO NT = 1, NTYPE
       RICE(NT) = R(NT)
    ENDDO
    SICE  = S

    ! ***** TIME LOOP
100 CONTINUE

    M = M + 1
    IF (M.GT.9999) THEN
       ERROR = 'K14: infinite time loop in XICE (1)'
       STATUS = 1
       RETURN
    ENDIF
    SCHGO = SCHG

    ! ***** REPEAT WITH REDUCED TIME STEP
10  CONTINUE
    M_MK = M_MK + 1
    IF (M_MK.GT.9999) THEN
       IF (SICE.lt.1._dp) THEN
          STOP_MK = 1
          GOTO 999
       ELSE
          ERROR = 'K14: infinite time loop in XICE (2)'
          STATUS = 1
          RETURN
       ENDIF
    ENDIF

    ! ***** COMPUTE SATURATION RATIO INCREMENT
    TEMP  = T - (TIME+DTIME) * COOLR
    IF (S1.EQ.0._dp) SSTOP=SCRHOM(TEMP)
    TEMP1 = 1._dp/ TEMP
    PICE  = WVP1 * EXP(-(WVP2*TEMP1))
    ADCHG = ( TEMP * T1 )**2.5_dp
    PRESS = P * ( TEMP * T1 ) * ADCHG
    RHOA  = 0.35_dp * PRESS * TEMP1
    VISC  = 1.E-5_dp * ( 1.512_dp + 0.0052_dp*(TEMP-233.15_dp) )
    FLUX  = ALP4 * SQRT(FVTH*TEMP)
    CISAT = THOUBK * PICE * TEMP1
    A1    = ( FA1 * TEMP1 - FA2 ) * TEMP1
    A2    = 1._dp/ CISAT
    A3    = FA3 * TEMP1 / PRESS
    B1    = FLUX * SVOL * CISAT * ( SICE-1._dp )
    B2    = FLUX * FD * PRESS * TEMP1**1.94_dp

    ! ***** SUM UP DEPOSITION TERM
    DLOSS      = 0._dp
    DO NT      = 1, NTYPE
       IF (C(NT).GT.0._dp) THEN
          CICE(NT) = ADCHG * C(NT)
          VFALL    = 14._dp* RICE(NT) * (1.3_dp/RHOA)**0.35_dp
          REYN     = RHOA * 2._dp* RICE(NT) * VFALL / VISC
          VENT     = 1._dp+ 0.229_dp * SQRT(REYN)
          DELP1    = 1._dp+ B2 * RICE(NT)
          ETA      = 1._dp+ 2._dp* VENT * B1 * B2 * DTIME / DELP1**2._dp
          RICE(NT) = ( DELP1 * SQRT(ETA) - 1._dp ) / B2
          DLOSS    = DLOSS + FPIVOL * CICE(NT) * VENT   &
               * B1 * RICE(NT)**2._dp / ( 1._dp+ B2 * RICE(NT) )
       ENDIF
    ENDDO
    VICE = ( A2 + A3 * SICE ) * DLOSS / ( A1 * SICE )  ! Eq. (13)
    DSICE = A1 * SICE * ( V*REAL(IV) - VICE ) * DTIME

    ! ***** CHECK ACCURACY AND PROVIDE UPDATE
    SCHG     = ABS(DSICE) / SICE
    IF (SCHG.GE.SERR) THEN
       DTIME   = DTMINUS * DTIME
       GOTO 10
    ELSE
       IF (TIME1.GE.0._dp) THEN
          SICE   = SICE + DSICE
          TIME   = TIME + DTIME
          IF (TIME.EQ.TIME1) GOTO 999
          YTIME  = DTPLUS * DTIME
          DTIME  = YTIME
          XTIME  = TIME + YTIME
          IF (XTIME.GT.TIME1) DTIME=YTIME-XTIME+TIME1
       ELSE
          STRY   = SICE + DSICE
          IF (STRY.EQ.SSTOP) THEN
             SICE  = STRY
             TIME  = TIME + DTIME
             GOTO 999
          ELSEIF (STRY.LT.SSTOP) THEN
             SICE  = STRY
             TIME  = TIME + DTIME
             IF (TIME.EQ.TIMER) GOTO 999
             YTIME = DTPLUS * DTIME
             DTIME = YTIME
             XTIME = TIME + YTIME
             IF (XTIME.GT.TIMER) DTIME=YTIME-XTIME+TIMER
          ELSEIF (STRY.LE.(SSTOP+SERR)) THEN
             SICE  = SSTOP + SERR
             TIME  = TIME + DTIME
             GOTO 999
          ELSE
             DTIME = DTMINUS * DTIME
             GOTO 10
          ENDIF
       ENDIF
    ENDIF

    ! ***** LOCATE STEADY-STATE IN S(T) THAT OCCURS WHEN VICE >> V,
    ! ***** NO FURTHER NUCLEATION EVENTS ARE ALLOWED
    IF ((ABS(SCHG).LE.SEPS).AND.(ABS(SCHGO).LE.SEPS)) THEN
       TIME = DTGM + TEPS
       GOTO 999
    ENDIF

    ! ***** COUNT THE TOTAL TIME ELAPSED DURING INTEGRATION,
    ! ***** WHICH MUST NOT EXCEED DTGM IN THIS ROUTINE
    IF ((TIME+TIMINC).GE.(DTGM-TEPS)) GOTO 999
    GOTO 100

    ! ***** EXIT
999 CONTINUE

    TAIR = TEMP
    PAIR = PRESS
    SSTOP_CIRRUS = SSTOP

    IF (TIME1.LT.0._dp) TIME1 = TIME
    STATUS = 0
    RETURN

  END SUBROUTINE XICE


  ! ===========================================================================

  SUBROUTINE XFRZHET(NAER, NFRZMOD, IFRZ, CVF, DT, C, R, SIG, SCRHET, &
                     P, T, V, SI, CI, RI, RS, T_start, STATUS, ERROR)

    ! ***** BERND KAERCHER  NOV 08  2004
    ! ***** bernd.kaercher@dlr.de  http://www.op.dlr.de/~pa1c/

    IMPLICIT NONE

    ! SUBROUTINE parameters
    INTEGER :: NAER
    INTEGER :: NFRZMOD
    INTEGER :: IFRZ

    INTEGER :: STATUS      ! STATUS
    CHARACTER(70) :: ERROR ! ERROR MESSAGE

    REAL(dp) :: CVF
    REAL(dp) :: DT
    REAL(dp) :: C(1:NFRZMOD)
    REAL(dp) :: R(1:NFRZMOD)
    REAL(dp) :: SIG(1:NFRZMOD)
    REAL(dp) :: SCRHET(2,1:NAER-1)
    REAL(dp) :: P
    REAL(dp) :: T
    REAL(dp) :: V
    REAL(dp) :: SI
    REAL(dp) :: CI
    REAL(dp) :: RI
    REAL(dp) :: RS
    REAL(dp) :: T_start

    ! local parameters
    REAL(dp) :: A(1:3)
    REAL(dp) :: B(1:2)
    REAL(dp) :: CCR(1:NFRZMOD)
    REAL(dp) :: XMI
    REAL(dp) :: CTOT
    REAL(dp) :: COOLR
    REAL(dp) :: PW
    REAL(dp) :: TMIN
    REAL(dp) :: SCR
    REAL(dp) :: PWCR
    REAL(dp) :: SISVE
    REAL(dp) :: TEMP
    REAL(dp) :: PICE
    REAL(dp) :: T1
    REAL(dp) :: T2
    REAL(dp) :: F1
    REAL(dp) :: F2
    REAL(dp) :: DS
    REAL(dp) :: PCORR
    REAL(dp) :: VOLF
    REAL(dp) :: PCR
    REAL(dp) :: SUPSCR
    REAL(dp) :: CTAU
    REAL(dp) :: TAU
    REAL(dp) :: DIFFC
    REAL(dp) :: BKT
    REAL(dp) :: VTH
    REAL(dp) :: FLUX
    REAL(dp) :: CISAT
    REAL(dp) :: THETA
    REAL(dp) :: XMI0
    REAL(dp) :: XMIMAX
    REAL(dp) :: RIMAX
    REAL(dp) :: TGROW
    REAL(dp) :: RIHAT
    REAL(dp) :: ZF
    REAL(dp) :: XMISAT
    REAL(dp) :: XMFP
    REAL(dp) :: BETA
    REAL(dp) :: X0
    REAL(dp) :: X
    REAL(dp) :: DX
    REAL(dp) :: XMID
    REAL(dp) :: FMID
    REAL(dp) :: YK

    INTEGER :: N
    INTEGER :: K
    INTEGER :: J

    ! parameters
    REAL(dp), PARAMETER :: KMAX   = 120.0_dp
    REAL(dp), PARAMETER :: RHOICE = 0.925_dp
    REAL(dp), PARAMETER :: PI  = 3.1415927_dp
    REAL(dp), PARAMETER :: DSCRMIN = 0.05_dp
    REAL(dp), PARAMETER :: THIRD  = 0.3333333_dp
    REAL(dp), PARAMETER :: RGAS  = 8.3145E7_dp
    REAL(dp), PARAMETER :: ALPHA =  0.5_dp
    REAL(dp), PARAMETER :: BK     = 1.3807E-16_dp
    REAL(dp), PARAMETER :: CPAIR = 1.00467E7_dp
    REAL(dp), PARAMETER :: HEAT =2830.3E7_dp
    REAL(dp), PARAMETER :: AVOG   = 6.02213E23_dp
    REAL(dp), PARAMETER :: GRAV  = 981.0_dp
    REAL(dp), PARAMETER :: SVOL =3.23E-23_dp
    REAL(dp), PARAMETER :: XMW    = 2.992E-23_dp
    REAL(dp), PARAMETER :: XWW   = 18.016_dp
    REAL(dp), PARAMETER :: XWA  = 28.966_dp

    ! ***** CHECK INPUT
    IF (NFRZMOD.LT.1) THEN
       ERROR = 'K14: invalid NFRZMOD value in XFRZHET'
       STATUS = 1
       RETURN
    ENDIF
    IF ((CVF.LT.0._dp).OR.(CVF.GT.0.99_dp)) THEN
       ERROR = 'K14: invalid CVF value in XFRZHET'
       STATUS = 1
       RETURN
    ENDIF
    DO N = 1, NFRZMOD
       IF (R(N).LE.1.E-7_dp) THEN
          ERROR = 'K14: invalid R value in XFRZHET'
          STATUS = 1
          RETURN
       ENDIF
       IF (SIG(N).LE.1._dp) THEN
          ERROR = 'K14: invalid SIG value in XFRZHET'
          STATUS = 1
          RETURN
       ENDIF
       IF ((NFRZMOD.GT.1).AND.(SIG(N).LT.1.1_dp)) THEN
          ERROR = 'K14: increase SIGMA above 1.1 in XFRZHET'
          STATUS = 1
          RETURN
       ENDIF
    ENDDO
    IF (DT.LT.0._dp) THEN
       ERROR = 'K14: invalid DT value in XFRZHET'
       STATUS = 1
       RETURN
    ENDIF
    IF (SI.LE.0._dp) THEN
       ERROR = 'K14: invalid SI value in XFRZHET'
       STATUS = 1
       RETURN
    ENDIF
    IF (P.LE.0._dp.OR.P.GT.1300.25_dp) THEN
       WRITE(*,*) 'K14: invalid P value = ',P
       ERROR = 'K14: invalid P value in XFRZHET'
       STATUS = 1
       RETURN
    ENDIF
    ! The saturation ratio is limited to S=2 before calling the parametrization
    ! This should always limit temperature below 246 K even after step-back
    IF (T.LT.150._dp .OR. T.GT.246._dp) THEN
       WRITE(*,*) 'K14: invalid T value = ',T
       ERROR = 'K14: invalid T value in XFRZHET'
       STATUS = 1
       RETURN
    ENDIF

    ! ***** DO NOTHING CRITERIA
    CI    = 0._dp
    RI    = 0._dp
    XMI   = 0._dp
    CTOT  = 0._dp
    DO N  = 1, NFRZMOD
       CTOT = CTOT + C(N)
    ENDDO
    IF ((CTOT.LE.0.0_dp).OR.(V.LE.0.0_dp)) THEN
       STATUS = 0
       RETURN
    ENDIF
    COOLR = GRAV * V / CPAIR
    PW    = SI * PISAT(T)

    ! ***** CHECK IF FREEZING CONDITIONS ARE MET WITHIN DT > 0
    IF (DT.GT.0._dp) THEN
       TMIN = T - COOLR * DT
       TMIN = MAX(TMIN, 130.0_dp)

       IF (T_start .lt. 220._dp) THEN
          SCR = SCRHET(2,IFRZ)
       ELSE
          SCR = SCRHET(1,IFRZ)
       ENDIF

       PWCR = PW * ( TMIN / T )**3.5_dp
       IF ((PWCR/PISAT(TMIN)).LT.SCR) THEN
          STATUS = 0
          RETURN
       ENDIF
    ENDIF

    ! ***** FREEZING TEMPERATURE
    IF (T_start.lt.220._dp) THEN
       SCR = SCRHET(2, IFRZ)
    ELSE
       SCR = SCRHET(1, IFRZ)
    ENDIF

    IF ((SI.GE.SCR).OR.(DT.EQ.0._dp)) THEN
       SISVE = SCR
       TEMP  = T
       PICE  = PISAT(TEMP)
       PW    = SCR * PICE
       PWCR  = PW
       GOTO 10
    ELSE
       SISVE = SI
       T1    = 130._dp
       T2    = T
       K     = 0
11     CONTINUE
       K = K + 1

       IF (T_start .lt. 220._dp) THEN
          F1    = SCRHET(2,IFRZ) / (PW*(T1/T)**3.5_dp / PISAT(T1))
       ELSE
          F1    = SCRHET(1,IFRZ) / (PW*(T1/T)**3.5_dp / PISAT(T1))
       ENDIF

       IF (T_start .ge. 220._dp) THEN
          F2    = SCRHET(2,IFRZ) / (PW*(T2/T)**3.5_dp / PISAT(T2))
       ELSE
          F2    = SCRHET(1,IFRZ) / (PW*(T2/T)**3.5_dp / PISAT(T2))
       ENDIF

       DS    = LOG(F2/F1)
       TEMP  = T2 - LOG(F2) * (T2-T1) / DS
       IF (ABS(TEMP-T2).LT.1.E-3_dp) GOTO 10
       T1    = T2
       T2    = TEMP
       IF (K.GT.10) THEN
          ERROR = 'K14: undetermined het. freezing temp. in XFRZHET'
          STATUS = 1
          RETURN
       ENDIF
       GOTO 11
    ENDIF
    STATUS = 0
    RETURN
10  CONTINUE

    IF (T_start .lt. 220._dp) THEN
       SCR = SCRHET(2, IFRZ)
    ELSE
       SCR = SCRHET(1, IFRZ)
    ENDIF

    PICE    = PISAT(TEMP)
    PWCR    = PW * (TEMP/T)**3.5_dp
    PCORR   = PWCR / PW
    VOLF    = PI * RHOICE / 0.75_dp
    PCR     = P  * PCORR
    SCR     = MAX( SCR, 1.001_dp )
    SUPSCR  = SCR - 1.0_dp + DSCRMIN
    DO N    = 1, NFRZMOD
       CCR(N) = C(N) * PCORR
    ENDDO

    ! ***** TIMESCALE OF THE FREEZING EVENT
    CTAU   = TEMP * ( 0.004_dp*TEMP-2._dp ) + 304.4_dp
    TAU    = 1.0_dp / ( CTAU*COOLR )

    ! ***** ICE CRYSTAL PROPERTIES AFTER THE FREEZING EVENT
    DIFFC  = 4.0122E-3_dp * TEMP**1.94_dp / PCR
    BKT    = BK   * TEMP
    VTH    = SQRT( 8.0_dp*BKT / (PI*XMW) )
    FLUX   = 0.25_dp * ALPHA * VTH
    CISAT  = 1.E3_dp * PICE / BKT
    THETA  = HEAT * XWW / (RGAS*TEMP)
    A(1)   = ( THETA/CPAIR - XWA/RGAS ) * ( GRAV/TEMP )
    A(2)   = 1.0_dp   / CISAT
    A(3)   = 0.001_dp*XWW**2._dp*HEAT**2._dp / ( AVOG*XWA*PCR*TEMP*CPAIR )
    B(1)   = SVOL * FLUX * CISAT * SUPSCR
    B(2)   = FLUX / DIFFC

    CALL XICEHET(NFRZMOD, V, TEMP, TAU, SCR, A, B, CCR, R, SIG, CI, RIHAT, &
                 RS, YK, STATUS, ERROR)
    IF (STATUS /= 0) RETURN

    XMI0   = VOLF * CI  * RIHAT**3._dp  * CVF
    XMIMAX = XMI0 + XMW * CISAT     * ( SCR-1.0_dp )
    RIMAX  = ( XMIMAX / (VOLF*CI) )**THIRD
    TGROW  = 0.75_dp / ( PI*DIFFC*CI*RIMAX )

    ! ***** VAPOR RELAXATION: ICE CRYSTAL PROPERTIES AFTER DT > 0
    IF (DT.EQ.0._dp) THEN
       RI    = RIHAT
       XMI   = XMI0
       SI    = SCR
       GOTO 22
    ENDIF
    ZF     = DT    / TGROW
    XMISAT = XMW   * CISAT
    XMFP   = 3.0_dp    * DIFFC / VTH
    BETA   = XMFP  / ( 0.75_dp*ALPHA*RIMAX )
    X0     = RIHAT / RIMAX
    X      = X0
    DX     = 1._dp - X0
    J      = 0
21  CONTINUE
    J      = J  + 1
    IF (J.GT.40) THEN
       ERROR = 'K14: undetermined ice crystal size after DT in HFRZHET'
       STATUS = 1
       RETURN
    ENDIF

    DX     = DX * 0.5_dp
    XMID   = X  + DX
    FMID   = TAUG(BETA,XMID,X0) - ZF
    IF (FMID.LT.0._dp) X=XMID
    IF ( (ABS(DX/X).LT.1.E-4_dp) .OR. &
         (FMID.EQ.0._dp) .OR. (X.GT.0.9999999_dp) ) GOTO 20
    GOTO 21
20  CONTINUE
    RI     = X    * RIMAX
    XMI    = VOLF * CI * RI**3._dp
    SI     = SCR  - ( XMI-XMI0 ) / XMISAT
22  CONTINUE

    STATUS = 0
    RETURN

  END SUBROUTINE XFRZHET


  ! ===========================================================================

  SUBROUTINE XICEHET(NFRZMOD, V, T, TAU, SCR, A, B, C, R, SIG, CI, RIHAT, &
                     RS, YK, STATUS, ERROR)

    ! ***** ICE CRYSTAL CONCENTRATION AND SIZE AFTER FREEZING EVENT

    IMPLICIT NONE

    ! SUBROUTINE parameters
    INTEGER :: NFRZMOD
    INTEGER :: STATUS      ! STATUS
    CHARACTER(70) :: ERROR ! ERROR MESSAGE

    REAL(dp) :: A(1:3)
    REAL(dp) :: B(1:2)
    REAL(dp) :: C(1:NFRZMOD)
    REAL(dp) :: R(1:NFRZMOD)
    REAL(dp) :: SIG(1:NFRZMOD)
    REAL(dp) :: V
    REAL(dp) :: T
    REAL(dp) :: TAU
    REAL(dp) :: SCR
    REAL(dp) :: CI
    REAL(dp) :: RIHAT
    REAL(dp) :: RS
    REAL(dp) :: YK

    ! parameters
    INTEGER,  PARAMETER :: IBIN = 40.0_dp

    ! local parameters
    REAL(dp) :: SHAPE(1:IBIN)
    REAL(dp) :: R0(1:IBIN+1)
    REAL(dp) :: CRINT(1:IBIN+1)
    REAL(dp) :: CIINT(1:IBIN+1)
    REAL(dp) :: VOLF
    REAL(dp) :: PHI
    REAL(dp) :: RIMFC
    REAL(dp) :: TBBT
    REAL(dp) :: CTOT
    REAL(dp) :: DELTA
    REAL(dp) :: DELP1
    REAL(dp) :: SYK
    REAL(dp) :: REN
    REAL(dp) :: RIM
    REAL(dp) :: RIRAT
    REAL(dp) :: XMIHAT
    REAL(dp) :: DLOGR0
    REAL(dp) :: SUMSH
    REAL(dp) :: SIGL
    REAL(dp) :: SHAPFC1
    REAL(dp) :: SHAPFC2
    REAL(dp) :: ARG
    REAL(dp) :: RMEAN
    REAL(dp) :: CSH
    REAL(dp) :: RLOGRAT
    REAL(dp) :: SLOPER
    REAL(dp) :: SLOPE

    INTEGER  :: NBIN
    INTEGER  :: N
    INTEGER  :: I
    INTEGER  :: II

    ! parameters
    REAL(dp), PARAMETER :: RMIN   = 1.E-7_dp
    REAL(dp), PARAMETER :: RMAX   = 1.E-3_dp
    REAL(dp), PARAMETER :: VRAT   = 2.0_dp
    REAL(dp), PARAMETER :: PI     = 3.1415927_dp
    REAL(dp), PARAMETER :: TWOPI  = 6.2831853_dp
    REAL(dp), PARAMETER :: RHOICE = 0.925_dp
    REAL(dp), PARAMETER :: SVOL   = 3.23E-23_dp
    REAL(dp), PARAMETER :: THIRD  = 0.3333333_dp
    REAL(dp), PARAMETER :: XMW    = 2.992E-23_dp
    REAL(dp), PARAMETER :: SQPI   = 1.7724539_dp


    ! ***** CONSTANTS
    NBIN = 1 + INT( LOG( (RMAX/RMIN)**3 ) / LOG(VRAT) )

    IF (NBIN.GT.IBIN) THEN
       ERROR = 'K14: increase value for IBIN in XICEHET'
       STATUS = 1
       RETURN
    ENDIF

    VOLF  = PI * RHOICE    / 0.75_dp
    PHI   = V  * A(1)*SCR / ( A(2) + A(3)*SCR )
    RIMFC = 4.0_dp * PI   * B(1)/B(2)**2._dp / SVOL
    TBBT  = 2.0_dp * B(1) * B(2) * TAU
    CTOT  = 0.0_dp
    DO N  = 1, NFRZMOD
       CTOT = CTOT + C(N)
    ENDDO

    ! ***** MONODISPERSE APPROXIMATION (SINGLE MODE ONLY)
    IF ((NFRZMOD.EQ.1).AND.(SIG(NFRZMOD).LT.1.1_dp)) THEN
       RS     = R(NFRZMOD)
       DELTA  = B(2) * RS
       DELP1  = 1.0_dp   + DELTA
       YK     = TBBT / DELP1**2._dp
       SYK    = SQRT(YK)
       REN    = 3.0_dp * SYK / ( 2.0_dp + SQRT(1.0_dp+9.0_dp*YK/PI) )
       RIM    = RIMFC / DELP1 * ( DELTA**2._dp - 1.0_dp + &
            (1.0_dp+0.5_dp*YK* DELP1**2._dp) * REN/SYK )

       CI     = PHI / RIM
       RIRAT  = 1.0_dp + 0.5_dp * SYK * REN
       RIHAT  = ( RIRAT * DELP1 - 1.0_dp ) / B(2)
       XMIHAT = VOLF * CI * RIHAT**3._dp
       CI     = MIN( CI, CTOT )
       RIHAT  = ( XMIHAT / (VOLF*CI) )**THIRD
       STATUS = 0
       RETURN
    ENDIF

    ! ***** SIZE DISTRIBUTION PROPERTIES
    R0(NBIN+1)    = RMAX * VRAT**THIRD
    CIINT(NBIN+1) = 1.E-35_dp
    CRINT(NBIN+1) = 1.E-25_dp
    DLOGR0     = 2.0_dp**THIRD * (VRAT**THIRD-1.0_dp) / (VRAT+1.0_dp)**THIRD
    SUMSH      = 0.0_dp
    DO I       = 1, NBIN
       CIINT(I)  = 0._dp
       CRINT(I)  = 0._dp
       SHAPE(I)  = 0._dp
       R0(I)     = RMIN * VRAT**( THIRD*FLOAT(I-1) )
       DO N      = 1, NFRZMOD
          SIGL     = LOG( MAX( SIG(N), 1.1_dp ) )
          SHAPFC1  = 1.0_dp   / ( SQRT(TWOPI) * SIGL )
          SHAPFC2  = 0.5_dp / SIGL**2._dp
          ARG      = SHAPFC2  * ( LOG(R0(I)/R(N)) )**2._dp
          ARG      = MIN( ARG, 75.0_dp )
          SHAPE(I) = SHAPE(I) + DLOGR0 * C(N) * SHAPFC1 * EXP(-ARG)
       ENDDO
       SUMSH = SUMSH + SHAPE(I)
    ENDDO

    ! ***** ICE CRYSTAL PROPERTIES
    RMEAN = 0.0_dp
    DO I  = NBIN, 1, -1
       DELTA = B(2) * R0(I)
       DELP1 = 1.0_dp   + DELTA
       YK    = TBBT / DELP1**2._dp
       SYK   = SQRT(YK)
       REN   = 3.0_dp * SYK / ( 2.0_dp + SQRT(1.0_dp+9.0_dp*YK/PI) )
       RIM   = RIMFC / DELP1 * ( DELTA**2._dp - 1.0_dp + &
            (1.0_dp+0.5_dp*YK*DELP1**2._dp) * REN/SYK )
       CSH      = SHAPE(I) / SUMSH   * CTOT
       CRINT(I) = CRINT(I+1) + RIM   * CSH
       CIINT(I) = CIINT(I+1) +         CSH
       RMEAN    = RMEAN      + R0(I) * CSH
       IF (CRINT(I).GE.PHI) GOTO 10
    ENDDO
    RS      = R0(1)
    RMEAN   = RMEAN / CTOT
    CI      = CTOT  * PHI   / CRINT(1)
    DELP1   = 1.0_dp    + B(2)  * RMEAN
    YK      = TBBT  / DELP1**2._dp
    SYK     = SQRT(YK)
    REN     = 3.0_dp * SYK / ( 2.0_dp + SQRT(1.0_dp+9.0_dp*YK/PI) )
    RIRAT   = 1.0_dp + 0.5_dp * SYK * REN
    RIHAT   = ( RIRAT * DELP1 - 1.0_dp ) / B(2)
    XMIHAT  = VOLF    * CI * RIHAT**3._dp
    CI      = CTOT
    RIHAT   = ( XMIHAT / (VOLF*CI) )**THIRD

    STATUS = 0
    RETURN

10  CONTINUE
    II = MAX(I, 1)
    RLOGRAT = LOG( R0(II) / R0(II+1) )
    SLOPER  = LOG( CRINT(II) / CRINT(II+1) ) / RLOGRAT
    RS      = R0(II+1)    * ( PHI / CRINT(II+1) )**(1._dp/SLOPER)
    SLOPE   = LOG( CIINT(II) / CIINT(II+1) ) / RLOGRAT
    CI      = CIINT(II+1) * ( RS / R0(II+1) )**SLOPE
    RMEAN   = RMEAN / CTOT
    DELP1   = 1.0_dp   + B(2) * MAX( RS, RMEAN )
    YK      = TBBT / DELP1**2.0_dp
    SYK     = SQRT(YK)
    REN     = 3.0_dp * SYK / ( 2.0_dp + SQRT(1.0_dp+9.0_dp*YK/PI) )
    RIRAT   = 1.0_dp + 0.5_dp * SYK * REN
    RIHAT   = ( RIRAT * DELP1 - 1.0_dp ) / B(2)

    STATUS = 0
    RETURN

  END SUBROUTINE XICEHET


  ! ===========================================================================

  SUBROUTINE XFRZHOM(NOSIZE, NFRZMOD, DT, C, R, SIG, P, T, V, &
                     SI, CI, RI, RS, STATUS, ERROR)

    ! ***** BERND KAERCHER  DEC 09  2003
    ! ***** bernd.kaercher@dlr.de  http://www.op.dlr.de/~pa1c/

    IMPLICIT NONE

    ! SUBROUTINE parameters
    LOGICAL :: CHECK
    LOGICAL :: NOSIZE

    INTEGER :: NFRZMOD
    INTEGER :: STATUS      ! STATUS
    CHARACTER(70) :: ERROR ! ERROR MESSAGE

    REAL(dp) :: C(1:NFRZMOD)
    REAL(dp) :: R(1:NFRZMOD)
    REAL(dp) :: SIG(1:NFRZMOD)
    REAL(dp) :: P
    REAL(dp) :: T
    REAL(dp) :: V
    REAL(dp) :: SI
    REAL(dp) :: CI
    REAL(dp) :: RI
    REAL(dp) :: RS
    REAL(dp) :: DT

    ! local parameters
    REAL(dp) :: A(1:3)
    REAL(dp) :: B(1:2)
    REAL(dp) :: CCR(1:NFRZMOD)
    REAL(dp) :: XMI
    REAL(dp) :: CTOT
    REAL(dp) :: COOLR
    REAL(dp) :: PW
    REAL(dp) :: TMIN
    REAL(dp) :: SCR
    REAL(dp) :: PWCR
    REAL(dp) :: SISVE
    REAL(dp) :: TEMP
    REAL(dp) :: PICE
    REAL(dp) :: T1
    REAL(dp) :: T2
    REAL(dp) :: F1
    REAL(dp) :: F2
    REAL(dp) :: DS
    REAL(dp) :: PCORR
    REAL(dp) :: VOLF
    REAL(dp) :: PCR
    REAL(dp) :: CTAU
    REAL(dp) :: TAU
    REAL(dp) :: DIFFC
    REAL(dp) :: BKT
    REAL(dp) :: VTH
    REAL(dp) :: FLUX
    REAL(dp) :: CISAT
    REAL(dp) :: THETA
    REAL(dp) :: RIHAT
    REAL(dp) :: YK
    REAL(dp) :: XMI0
    REAL(dp) :: XMIMAX
    REAL(dp) :: XMISAT
    REAL(dp) :: RIMAX
    REAL(dp) :: TGROW
    REAL(dp) :: ZF
    REAL(dp) :: XMFP
    REAL(dp) :: BETA
    REAL(dp) :: X0
    REAL(dp) :: X
    REAL(dp) :: DX
    REAL(dp) :: XMID
    REAL(dp) :: FMID


    INTEGER :: N
    INTEGER :: K
    INTEGER :: J

    ! parameters
    REAL(dp), PARAMETER :: KMAX   = 120.0_dp
    REAL(dp), PARAMETER :: RHOICE = 0.925_dp
    REAL(dp), PARAMETER :: PI     = 3.1415927_dp
    REAL(dp), PARAMETER :: ALPHA  = 0.5_dp
    REAL(dp), PARAMETER :: THIRD  = 0.3333333_dp
    REAL(dp), PARAMETER :: RGAS   = 8.3145E7_dp
    REAL(dp), PARAMETER :: BK     = 1.3807E-16_dp
    REAL(dp), PARAMETER :: CPAIR  = 1.00467E7_dp
    REAL(dp), PARAMETER :: HEAT   = 2830.3E7_dp
    REAL(dp), PARAMETER :: AVOG   = 6.02213E23_dp
    REAL(dp), PARAMETER :: GRAV   = 981.0_dp
    REAL(dp), PARAMETER :: SVOL   = 3.23E-23_dp
    REAL(dp), PARAMETER :: XMW    = 2.992E-23_dp
    REAL(dp), PARAMETER :: XWW    = 18.016_dp
    REAL(dp), PARAMETER :: XWA    = 28.966_dp

    ! ***** CHECK INPUT
    IF (NFRZMOD.LT.1) THEN
       ERROR = 'K14: invalid NFRZMOD value in XFRZHOM'
       STATUS = 1
       RETURN
    ENDIF
    DO N = 1, NFRZMOD
       IF (R(N).LE.1.E-7_dp) THEN
          ERROR = 'K14: invalid R value in XFRZHOM'
          STATUS = 1
          RETURN
       ENDIF
       IF (SIG(N).LE.1._dp) THEN
          ERROR = 'K14: invalid SIG value in XFRZHOM'
          STATUS = 1
          RETURN
       ENDIF
       IF ((NFRZMOD.GT.1).AND.(SIG(N).LT.1.1_dp)) THEN
          ERROR = 'K14: increase SIGMA above 1.1 in XFRZHOM'
          STATUS = 1
          RETURN
       ENDIF
    ENDDO
    IF (DT.LT.0._dp) THEN
       ERROR = 'K14: invalid DT value in XFRZHOM'
       STATUS = 1
       RETURN
    ENDIF
    IF (SI.LE.0._dp) THEN
       ERROR = 'K14: invalid SI value in XFRZHOM'
       STATUS = 1
       RETURN
    ENDIF
    IF (P.LE.0._dp.OR.P.GT.1300.25_dp) THEN
       WRITE(*,*) 'K14: invalid P value = ',P
       ERROR = 'K14: invalid P value in XFRZHOM'
       STATUS = 1
       RETURN
    ENDIF
    ! The saturation ratio is limited to S=2 before calling the parametrization
    ! This should always limit temperature below 246 K even after step-back
    IF (T.LT.150._dp .OR. T.GT.246._dp) THEN
       WRITE(*,*) 'K14: invalid T value = ',T
       ERROR = 'K14: invalid T value in XFRZHET'
       STATUS = 1
       RETURN
    ENDIF

    ! ***** DO NOTHING CRITERIA
    CI    = 0.0_dp
    RI    = 0.0_dp
    XMI   = 0.0_dp
    CTOT  = 0.0_dp
    DO N  = 1, NFRZMOD
       CTOT = CTOT + C(N)
    ENDDO
    IF ((CTOT.LE.0._dp).OR.(V.LE.0._dp)) THEN
       STATUS = 0
       RETURN
    ENDIF
    COOLR = GRAV * V / CPAIR
    PW    = SI * PISAT(T)

    ! ***** CHECK IF FREEZING CONDITIONS ARE MET WITHIN DT > 0
    IF (DT.GT.0._dp) THEN
       TMIN = T - COOLR * DT
       TMIN = MAX( TMIN, 130.0_dp )
       SCR  = SCRHOM(TMIN)
       PWCR = PW * ( TMIN / T )**3.5_dp
       IF ((PWCR/PISAT(TMIN)).LT.SCR) THEN
          STATUS = 0
          RETURN
       ENDIF
    ENDIF

    ! ***** FREEZING TEMPERATURE
    SCR    = SCRHOM(T)
    IF ((SI.GE.SCR).OR.(DT.EQ.0._dp)) THEN
       SISVE = SCR
       TEMP  = T
       PICE  = PISAT(TEMP)
       PW    = SCR * PICE
       PWCR  = PW
       GOTO 10
    ELSE
       SISVE = SI
       T1    = 130.0_dp
       T2    = T
       K     = 0
11     CONTINUE
       K     = K + 1
       F1    = SCRHOM(T1) / (PW*(T1/T)**3.5_dp / PISAT(T1))
       F2    = SCRHOM(T2) / (PW*(T2/T)**3.5_dp / PISAT(T2))
       DS    = LOG(F2/F1)
       TEMP  = T2 - LOG(F2) * (T2-T1) / DS
       IF (ABS(TEMP-T2).LT.1.E-3_dp) GOTO 10
       T1    = T2
       T2    = TEMP
       IF (K.GT.10) THEN
          ERROR = 'K14: undetermined het. freezing temp. in XFRZHOM'
          STATUS = 1
          RETURN
       ENDIF
       GOTO 11
    ENDIF
    STATUS = 0
    RETURN
10  CONTINUE
    SCR     = SCRHOM(TEMP)
    PICE    = PISAT(TEMP)
    PWCR    = PW * (TEMP/T)**3.5_dp
    PCORR   = PWCR / PW
    VOLF    = PI * RHOICE / 0.75_dp
    PCR     = P  * PCORR
    DO N    = 1, NFRZMOD
       CCR(N) = C(N) * PCORR
    ENDDO

    ! ***** TIMESCALE OF THE FREEZING EVENT
    CTAU = TEMP * ( 0.004_dp*TEMP-2._dp ) + 304.4_dp
    TAU  = 1._dp / ( CTAU*COOLR )

    ! ***** ICE CRYSTAL PROPERTIES AFTER THE FREEZING EVENT
    DIFFC  = 4.0122E-3_dp * TEMP**1.94_dp / PCR
    BKT    = BK   * TEMP
    VTH    = SQRT( 8.0_dp*BKT / (PI*XMW) )
    FLUX   = 0.25_dp * ALPHA * VTH
    CISAT  = 1.E3_dp * PICE  / BKT
    THETA  = HEAT * XWW   / (RGAS*TEMP)
    A(1)   = ( THETA/CPAIR - XWA/RGAS ) * ( GRAV/TEMP )
    A(2)   = 1.0_dp   / CISAT
    A(3)   = 0.001_dp*XWW**2._dp*HEAT**2._dp / ( AVOG*XWA*PCR*TEMP*CPAIR )
    B(1)   = SVOL * FLUX * CISAT * ( SCR-1.0_dp )
    B(2)   = FLUX / DIFFC
    CALL XICEHOM(CHECK, NOSIZE, NFRZMOD, V, TEMP, TAU, SCR, A, B, CCR, R, &
                 SIG, CI, RIHAT, RS, YK, STATUS, ERROR)
    IF (STATUS /= 0) RETURN
    XMI0   = VOLF * CI  * RIHAT**3._dp
    XMIMAX = XMI0 + XMW * CISAT * ( SCR-1.0_dp )
    RIMAX  = ( XMIMAX / (VOLF*CI) )**THIRD
    TGROW  = 0.75_dp / ( PI*DIFFC*CI*RIMAX )

    ! ***** VAPOR RELAXATION: ICE CRYSTAL PROPERTIES AFTER DT > 0
    IF (DT.EQ.0._dp) THEN
       RI    = RIHAT
       XMI   = XMI0
       SI    = SCR
       GOTO 22
    ENDIF
    ZF     = DT    / TGROW
    XMISAT = XMW   * CISAT
    XMFP   = 3.0_dp    * DIFFC / VTH
    BETA   = XMFP  / ( 0.75_dp*ALPHA*RIMAX )
    X0     = RIHAT / RIMAX
    X      = X0
    DX     = 1._dp - X0
    J      = 0
21  CONTINUE
    J      = J  + 1
    IF (J.GT.30) THEN
       ERROR = 'K14: undetermined ice crystal size after DT in HFRZHOM'
       STATUS = 1
       RETURN
    ENDIF

    DX     = DX * 0.5_dp
    XMID   = X  + DX
    FMID   = TAUG(BETA,XMID,X0) - ZF
    IF (FMID.LT.0._dp) X=XMID
    IF ( (ABS(DX/X).LT.1.E-4_dp) .OR.      &
         (FMID.EQ.0._dp) .OR. (X.GT.0.9999999_dp) ) GOTO 20
    GOTO 21
20  CONTINUE
    RI     = X    * RIMAX
    XMI    = VOLF * CI * RI**3._dp
    SI     = SCR  - ( XMI-XMI0 ) / XMISAT
22  CONTINUE

    STATUS = 0
    RETURN

  END SUBROUTINE XFRZHOM


  ! ===========================================================================

  SUBROUTINE XICEHOM (CHECK, NOSIZE, NFRZMOD, V, T, TAU, SCR, A, B, C, R, &
                      SIG, CI, RIHAT, RS, YK, STATUS, ERROR)

    ! ***** ICE CRYSTAL CONCENTRATION AND SIZE AFTER FREEZING EVENT

    IMPLICIT NONE

    ! SUBROUTINE parameters
    LOGICAL :: CHECK
    LOGICAL :: NOSIZE

    INTEGER :: NFRZMOD
    INTEGER :: STATUS      ! STATUS
    CHARACTER(70) :: ERROR ! ERROR MESSAGE

    REAL(dp) :: V
    REAL(dp) :: T
    REAL(dp) :: TAU
    REAL(dp) :: SCR
    REAL(dp) :: A(1:3)
    REAL(dp) :: B(1:2)
    REAL(dp) :: C(1:NFRZMOD)
    REAL(dp) :: R(1:NFRZMOD)
    REAL(dp) :: SIG(1:NFRZMOD)
    REAL(dp) :: CI
    REAL(dp) :: RIHAT
    REAL(dp) :: RS
    REAL(dp) :: YK

    ! local parameters
    INTEGER, PARAMETER :: IBIN = 40
    REAL(dp) :: SHAPE(1:IBIN)
    REAL(dp) :: R0(1:IBIN+1)
    REAL(dp) :: CRINT(1:IBIN+1)
    REAL(dp) :: CIINT(1:IBIN+1)
    REAL(dp) :: VOLF
    REAL(dp) :: PHI
    REAL(dp) :: RIMFC
    REAL(dp) :: TBBT
    REAL(dp) :: CTOT
    REAL(dp) :: XMIHAT
    REAL(dp) :: DELTA
    REAL(dp) :: DELP1
    REAL(dp) :: SYK
    REAL(dp) :: REN
    REAL(dp) :: RIM
    REAL(dp) :: RIRAT
    REAL(dp) :: DLOGR0
    REAL(dp) :: SUMSH
    REAL(dp) :: SIGL
    REAL(dp) :: SHAPFC1
    REAL(dp) :: SHAPFC2
    REAL(dp) :: ARG
    REAL(dp) :: RMEAN
    REAL(dp) :: CSH
    REAL(dp) :: RLOGRAT
    REAL(dp) :: SLOPER
    REAL(dp) :: SLOPE

    INTEGER :: NBIN
    INTEGER :: N
    INTEGER :: I
    INTEGER :: II

    ! parameters
    REAL(dp), PARAMETER :: RMIN = 1.E-7_dp
    REAL(dp), PARAMETER :: RMAX = 1.E-3_dp
    REAL(dp), PARAMETER :: VRAT = 2.0_dp
    REAL(dp), PARAMETER :: PI   = 3.1415927_dp
    REAL(dp), PARAMETER :: TWOPI = 6.2831853_dp
    REAL(dp), PARAMETER :: RHOICE = 0.925_dp
    REAL(dp), PARAMETER :: SVOL = 3.23E-23_dp
    REAL(dp), PARAMETER :: THIRD = 0.3333333_dp
    REAL(dp), PARAMETER :: XMW = 2.992E-23_dp
    REAL(dp), PARAMETER :: SQPI   = 1.7724539_dp

    ! ***** CONSTANTS
    NBIN = 1 + INT( LOG( (RMAX/RMIN)**3._dp ) / LOG(VRAT) )
    IF (NBIN.GT.IBIN) THEN
       ERROR = 'K14: increate IBIN value in XICEHOM'
       STATUS = 1
       RETURN
    ENDIF
    VOLF  = PI * RHOICE    / 0.75_dp
    PHI   = V  * A(1)*SCR / ( A(2) + A(3)*SCR )
    RIMFC = 4.0_dp * PI   * B(1)/B(2)**2._dp / SVOL
    TBBT  = 2.0_dp * B(1) * B(2) * TAU
    CTOT  = 0.0_dp
    DO N  = 1, NFRZMOD
       CTOT = CTOT + C(N)
    ENDDO

    ! ***** NO SIZE EFFECTS
    IF (NOSIZE) THEN
       RS = 0.25E-4_dp
       YK = TAU  * ( B(1)/RS ) / ( 1.0_dp + B(2)*RS )
       CI = SVOL * ( B(2) / (TWOPI*B(1)) )**1.5_dp * PHI / SQRT(TAU)
       CI = MIN( CI, CTOT )
       XMIHAT = XMW * PI * PHI * TAU / 6.0_dp
       RIHAT  = ( XMIHAT / (VOLF*CI) )**THIRD
       STATUS = 0
       RETURN
    ENDIF

    ! ***** MONODISPERSE APPROXIMATION (SINGLE MODE ONLY)
    IF ((NFRZMOD.EQ.1).AND.(SIG(NFRZMOD).LT.1.1_dp)) THEN
       RS     = R(NFRZMOD)
       DELTA  = B(2) * RS
       DELP1  = 1.0_dp   + DELTA
       YK     = TBBT / DELP1**2._dp
       SYK    = SQRT(YK)
       REN    = 3.0_dp * SYK / ( 2.0_dp + SQRT(1.0_dp+9.0_dp*YK/PI) )
       RIM    = RIMFC / DELP1 * ( DELTA**2._dp - 1.0_dp + &
            (1.0_dp+0.5_dp*YK*DELP1**2._dp) * REN/SYK )
       CI     = PHI / RIM
       RIRAT  = 1.0_dp + 0.5_dp * SYK * REN
       RIHAT  = ( RIRAT * DELP1 - 1.0_dp ) / B(2)
       XMIHAT = VOLF * CI * RIHAT**3._dp
       CI     = MIN( CI, CTOT )
       RIHAT  = ( XMIHAT / (VOLF*CI) )**THIRD
       STATUS = 0
       RETURN
    ENDIF

    ! ***** SIZE DISTRIBUTION PROPERTIES
    R0(NBIN+1)    = RMAX * VRAT**THIRD
    CIINT(NBIN+1) = 1.E-35_dp
    CRINT(NBIN+1) = 1.E-25_dp
    DLOGR0 = 2.0_dp**THIRD * (VRAT**THIRD-1.0_dp) / (VRAT+1.0_dp)**THIRD
    SUMSH = 0.0_dp
    DO I = 1, NBIN
       CIINT(I)  = 0.0_dp
       CRINT(I)  = 0.0_dp
       SHAPE(I)  = 0.0_dp
       R0(I)     = RMIN * VRAT**( THIRD*FLOAT(I-1) )
       DO N  = 1, NFRZMOD
          SIGL     = LOG( MAX( SIG(N), 1.1_dp ) )
          SHAPFC1  = 1.0_dp   / ( SQRT(TWOPI) * SIGL )
          SHAPFC2  = 0.5_dp / SIGL**2._dp
          ARG      = SHAPFC2  * ( LOG(R0(I)/R(N)) )**2._dp
          ARG      = MIN( ARG, 75.0_dp )
          SHAPE(I) = SHAPE(I) + DLOGR0 * C(N) * SHAPFC1 * EXP(-ARG)
       ENDDO
       SUMSH = SUMSH + SHAPE(I)
    ENDDO

    ! ***** ICE CRYSTAL PROPERTIES
    RMEAN     = 0.0_dp
    DO I = NBIN, 1, -1
       DELTA    = B(2) * R0(I)
       DELP1    = 1.0_dp   + DELTA
       YK       = TBBT / DELP1**2._dp
       SYK      = SQRT(YK)
       REN      = 3.0_dp * SYK / ( 2.0_dp + SQRT(1.0_dp+9.0_dp*YK/PI) )
       RIM      = RIMFC / DELP1 * ( DELTA**2._dp - 1.0_dp + &
            (1.0_dp+0.5_dp*YK*DELP1**2._dp) * REN/SYK )
       CSH      = SHAPE(I) / SUMSH   * CTOT
       CRINT(I) = CRINT(I+1) + RIM   * CSH
       CIINT(I) = CIINT(I+1) +         CSH
       RMEAN    = RMEAN      + R0(I) * CSH
       IF (CRINT(I).GE.PHI) GOTO 10
    ENDDO
    RS     = R0(1)
    RMEAN  = RMEAN / CTOT
    CI     = CTOT  * PHI   / CRINT(1)
    DELP1  = 1.0_dp + B(2)  * RMEAN
    YK     = TBBT  / DELP1**2._dp
    SYK    = SQRT(YK)
    REN    = 3.0_dp * SYK / ( 2.0_dp + SQRT(1.0_dp+9.0_dp*YK/PI) )
    RIRAT  = 1.0_dp + 0.5_dp * SYK * REN
    RIHAT  = ( RIRAT * DELP1 - 1.0_dp ) / B(2)
    XMIHAT = VOLF    * CI * RIHAT**3._dp
    CI     = CTOT
    RIHAT  = ( XMIHAT / (VOLF*CI) )**THIRD
    STATUS = 0
    RETURN

10  CONTINUE

    II      = MAX( I, 1)
    RLOGRAT = LOG( R0(II) / R0(II+1) )
    SLOPER  = LOG( CRINT(II) / CRINT(II+1) ) / RLOGRAT
    RS      = R0(II+1)    * ( PHI / CRINT(II+1) )**(1.0_dp/SLOPER)
    SLOPE   = LOG( CIINT(II) / CIINT(II+1) ) / RLOGRAT
    CI      = CIINT(II+1) * ( RS / R0(II+1) )**SLOPE
    RMEAN   = RMEAN / CTOT
    DELP1   = 1.0_dp   + B(2) * MAX( RS, RMEAN )
    YK      = TBBT / DELP1**2._dp
    SYK     = SQRT(YK)
    REN     = 3.0_dp * SYK / ( 2.0_dp + SQRT(1.0_dp+9.0_dp*YK/PI) )
    RIRAT   = 1.0_dp + 0.5_dp * SYK * REN
    RIHAT   = ( RIRAT * DELP1 - 1.0_dp ) / B(2)

    ! ***** EXIT
    STATUS = 0
    RETURN

  END SUBROUTINE XICEHOM

  ! ===========================================================================

  SUBROUTINE NUCTYPE (NIN, NAER, INUC, CVF, STATUS, ERROR)

    IMPLICIT NONE

    ! SUBROUTINE parameters
    INTEGER :: NIN
    INTEGER :: NAER
    INTEGER :: INUC
    INTEGER :: STATUS      ! STATUS
    CHARACTER(70) :: ERROR ! ERROR MESSAGE
    REAL(dp) :: CVF

    IF (NIN.LT.NAER) THEN

       INUC = NIN
       CVF = 0.99_dp
       STATUS = 0
       RETURN

    ELSE

       ERROR = 'K14: NINC value not valid for het. IN in NUCTYPE'
       STATUS = 1
       RETURN

    ENDIF

  END SUBROUTINE NUCTYPE

  ! ===========================================================================

  SUBROUTINE cloud_read_nml_ctrl_k14(iou, status, modstr)

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

   ! I/O
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit
    INTEGER, INTENT(OUT) :: status ! error status
    CHARACTER(LEN=*), INTENT(IN) :: modstr

    NAMELIST /CTRL_K14/ cdncmin, nfrzaer, nexti, inp_properties, &
         vervel_scheme, scale_v_ls_cirrus, scale_v_tke_cirrus, &
         scale_v_orogw_cirrus, l_ic2snow, du_nfrac_thr

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'cloud_read_nml_ctrl_k14'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    ! INITIALIZE GLOBAL CONTROL VARIABLES
    ! -> DEFAULT VALUES ARE SET AT DECLARATION ABOVE

    CALL read_nml_open(lex, substr, iou, 'CTRL_K14', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL_K14, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL_K14', modstr)

    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE cloud_read_nml_ctrl_k14

END MODULE MESSY_CLOUD_KUEBBELER14
