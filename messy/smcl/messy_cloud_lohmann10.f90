MODULE MESSY_CLOUD_LOHMANN10

  ! Author: Holger Tost
  ! original code from:  Ulrike Lohmann, ETH Zurich
  !                      Axel Lauer
  !                      Sylvaine Ferrachat

  USE messy_main_constants_mem,   ONLY: dp, api => pi, g, ak => k_B, m_air 
  USE messy_cloud_ori,            ONLY: crhoi, cqtmin, crt
  USE messy_cloud_mem

  IMPLICIT NONE
  PRIVATE
  SAVE

  ! new &CTRL_L10 namelist parameter
  ! minimum cloud droplet number concentration [m^(-3)]
  REAL(dp), PUBLIC :: cdncmin = 40.e6_dp
  !# immersion nucleation is computed by Barahona and Nenes, 2009b 
  !# (F = by Lohmann and Diehl, 2006)
  LOGICAL, PUBLIC  :: limm_BN09 = .FALSE.

  INTEGER, PUBLIC :: nauto = 2
  INTEGER, PUBLIC :: ncvmicro = 0

  LOGICAL :: lomonodisp  = .FALSE.     ! switch to turn on simplified freezing with monodisperse
                                       ! freezing aerosol modes
                                       ! lomonodisp = .FALSE. (default, full physics)
                                       !            = .TRUE. monodisperse distrib (faster with 
                                       !                    no substancial degradation in physics)
  LOGICAL :: lthermo      = .FALSE.    ! switch to take into account thermophoresis
                                       ! in freezing of cloud water
                                       ! lthermo = .FALSE. (default, no thermophoresis)
                                       !         = .TRUE. (thermophoresis term)
!----------------
! Public entities
!----------------

  REAL(dp), PARAMETER :: zepsec = 1.0e-12_dp
  REAL(dp), PARAMETER :: zxsec  = 1._dp - zepsec
  REAL(dp), PARAMETER :: zqsec  = 1._dp - cqtmin 
  REAL(dp), PARAMETER :: zeps   = EPSILON(1.0_dp)
  REAL(dp), PARAMETER :: zcri   = 10.E-6_dp  ! to estimate the number of produced  
                                             ! cloud droplets from ice melting in  
                                             ! case of licnc=.FALSE. [m]=> 10 um
  REAL(dp), PARAMETER :: zmi = 4._dp/3._dp*zcri**3*api*crhoi ! assumed mass of ice crystals with 
                                                             ! corresponding volume mean radius zcri

  REAL(dp), PARAMETER :: zicemin  = 10._dp
  REAL(dp), PARAMETER :: zicemax  = 1.e7_dp
  REAL(dp), PARAMETER :: zsigmaw  = 0.28_dp 
  REAL(dp), PARAMETER :: zcdi     = 0.6_dp
  REAL(dp), PARAMETER :: zmw0     = 4.19e-12_dp
  REAL(dp), PARAMETER :: zmi0     = 1.e-12_dp
  REAL(dp), PARAMETER :: zmi0_rcp = 1.e12_dp
  REAL(dp), PARAMETER :: zka      = 0.024_dp
  REAL(dp), PARAMETER :: zkbc     = 4.2_dp     !thermal conductivity of carbon (Seinfeld&Pandis, p.481)
  REAL(dp), PARAMETER :: zkdu     = 0.72_dp    !thermal conductivity of clay (Seinfeld&Pandis, p.481)
  REAL(dp), PARAMETER :: zkb      = 1.38e-23_dp
  REAL(dp), PARAMETER :: zalpha   = 0.5_dp        ! deposition coefficent alpha
  REAL(dp), PARAMETER :: zxmw     = 2.992e-26_dp  ! mass of a h2o molecule [kg]
  REAL(dp), PARAMETER :: zfall    = 3._dp
  REAL(dp), PARAMETER :: zrhoice  = 925._dp

    
  PUBLIC :: cloud_cdnc_icnc2
  PUBLIC :: cloud_cdnc_icnc_cc
  PUBLIC :: cloud_read_nml_ctrl_l10

CONTAINS
!===============================================================================

! main microphysics cloud rountine

  SUBROUTINE cloud_cdnc_icnc2(kproma, kbdim, ktdia, klev, klevp1, ztmst &
!---Included for in-cloud scavenging (Philip Stier, 28/03/01):----------
                         , ktrac,  krow                                &
!---End Included for scavenging-----------------------------------------
!-----------------------------------------------------------------------
! - INPUT  2D .
                         , paphm1,   pvervel                           &
                         , papm1,    papp1,    pacdnc                  &
                         , pqm1,     ptm1,     ptvm1                   &
                         , pxlm1,    pxim1,    pxtec                   &
                         , pxvar,    pxskew,   pqtec                   &
                         , pbetaa,   pbetab                            &
                         , pvdiffp,  phmixtau, pvmixtau                &
                         , pgeo,     pbetass                           &
!---Included for in-cloud scavenging (Philip Stier, 28/03/01):----------
                         , pxtm1                                       &
!---End Included for scavenging-----------------------------------------
!--- Included for prognostic CDNC/IC scheme ----------------------------
                         , ptkem1,   pcvcbot,  pwcape                  &
                         , pxtecl,   pxteci                            &
                         , pxtecnl,  pxtecni                           &
!--- End included for CDNC/IC scheme -----------------------------------
! - INPUT  1D .
                         , knvb                                        &
! - OUTPUT 2D .
                         , paclc,    paclcac,  prelhum                 &
! - OUTPUT 3D (vertical velocity)
                         , w_sub, w_ave, w_grid                        &
! - INPUT 3D (for BN09 scheme)
                         , sigwBN,     ndropBN,     dsulfBN            &
                         , ndust_aiBN, ddust_aiBN                      &
                         , ndust_ciBN, ddust_ciBN                      &
                         , norgBN,     dorgBN                          &
                         , nsootBN,    dsootBN                         &
                         , nsolo_ksBN, dsolo_ksBN                      &
                         , nsolo_asBN, dsolo_asBN                      &
                         , nsolo_csBN, dsolo_csBN                      &
! - OUTPUT 3D (BN09 scheme)
                         , smaxice_cirrusBN, sc_ice_cirrusBN           &
                         , nlim_cirrusBN, nhet_cirrusBN, nice_cirrusBN &
                         , dice_cirrusBN, sigwpre_cirrusBN             &
                         , smaxice_immBN, sc_ice_immBN                 &
                         , nice_immBN, dice_immBN, sigwpre_immBN       &
! - OUTPUT 3D (ice crystals)
                         , newIC_cirri, newICR_cirri, newIC_cnt_therm  &
                         , newIC_imm,  newIC_mix,  newIC_mix_fin       &
! - INPUT/OUTPUT 1D .
                         , paclcov,  paprl,    pqvi                    &
                         , pxlvi,    pxivi                             &
! - OUTPUT 1D .
                         , pssfl,    prsfl                             &
! - INPUT/OUTPUT 2D .
                         , pqte,     ptte                              &
                         , pxlte,    pxite                             &
!---Included for in-cloud scavenging (Philip Stier, 28/03/01):----------
                         , pxtte                                       &
!---End Included for scavenging-----------------------------------------
! - INPUT/OUTPUT 1D .
                         , paprs                                       &
                         , status_string,      lcover                  &
                         , slm,      glac,     pcdncact                &
                         , plwc,     piwc                              &
                         , pfrain,   pfsnow                            &
                         , pfrain_no,pfsnow_no                         &
                         , prate_r,  prate_s                           &
                         , prevap,   pssubl                            &
                         , pr_cover, pcond                             &
                         , pimelt,   pisedi                            &
                         ) 
!
!     *Cloud* computes large-scale water phase changes, precipitation,
!             cloud cover, and vertical integrals of specific humidity,
!             cloud liquid water content and cloud ice (diagnostics).
!
!     Subject.
!     --------
!
!          This routine computes the tendencies of the four prognostic
!          variables (temperature t, specific humidity q, cloud liquid
!          water xl, cloud ice xi) due to phase changes (condensation/
!          deposition, evaporation/sublimation of rain/snow falling
!          into the unsaturated part of the grid box, melting of snow,
!          melting/freezing of cloud ice/cloud water, sedimentation of
!          cloud ice, and precipitation formation in warm, cold and
!          mixed phase clouds.
!          The precipitation at the surface (rain and snow) is used in
!          later for computing the land surface hydrology in *surf*.
!          The cloud parameters (cloud cover, cloud liquid water and
!          cloud ice are used for the calculation of radiation at the
!          next timestep.
!          Attention: 
!          In the current version the advective tendencies of skewness 
!          and variance are set to zero.
!
!     INTERFACE.
!     ----------
!
!     *Call cloud*
!
!     Input arguments.
!     ----- ----------
!  - 2D
!  paphm1   : pressure at half levels                              (n-1)
!  papm1    : pressure at full levels                              (n-1)
!  papp1    : pressure at full levels                              (n-1)
!  pacdnc   : cloud droplet number concentration (specified)
!  pqm1     : specific humidity                                    (n-1)
!  ptm1     : temperature                                          (n-1)
!  pxlm1    : cloud liquid water                                   (n-1)
!  pxim1    : cloud ice                                            (n-1)
!  pxtec    : detrained convective cloud liquid water or cloud ice (n)
!  pxvar    : distribution width (b-a)                             (n-1)
!  pxskew   : beta shape parameter "q"                             (n-1)
!  pbetaa   : the beta distribution minimum a                      (n-1)
!  pbetab   : the beta distribution maximum b                      (n-1)
!  pvdiffp  : the rate of change of q due to vdiff scheme          (n-1)
!  phmixtau : mixing timescale**-1 for horizontal turbulence       (n)
!  pvmixtau : mixing timescale**-1 for horizontal turbulence       (n)
!--- Included for prognostic CDNC/IC scheme ----------------------------
!  zcdnc    : in-cloud value, cloud droplet number concentration [m-3] 
!             = pxtm1(1,1,idt_cdnc)*density(air)
!  zicnc    : in-cloud value, ice crystal number concentration   [m-3] 
!             = pxtm1(1,1,idt_icnc)*density(air)
!--- End included for CDNC/IC scheme -----------------------------------
!
!  - 1D
!  knvb     :
!
!     Output arguments.
!     ------ ----------
!  - 1D
!  prsfl    : surface rain flux
!  pssfl    : surface snow flux
!
!     Input, Output arguments.
!     ------------ ----------
!  - 2D
!  paclc    : cloud cover  (now diagnosed in cover)
!  paclcac  : cloud cover, accumulated
!  paclcov  : total cloud cover
!  paprl    : total stratiform precipitation (rain+snow), accumulated
!  pqvi     : vertically integrated spec. humidity, accumulated
!  pxlvi    : vertically integrated cloud liquid water, accumulated
!  pxivi    : vertically integrated cloud ice, accumulated
!  ptte     : tendency of temperature
!  pqte     : tendency of specific humidity
!  pxlte    : tendency of cloud liquid water
!  pxite    : tendency of cloud ice
!--- Included for prognostic CDNC/IC scheme ----------------------------
!  pxtte    : tendency of cloud droplet number, in-cloud 
!  pxtte    : tendency of ice crystal number,in-cloud
!--- End included for CDNC/IC scheme -----------------------------------
!  - 1D
!  paprs    : Snowfall, accumulated
!
!     Externals.
!     ----------
!
!     Method.
!     -------
!     see References
!
!     References.
!     ----------
!
!     Lohmann and Roeckner, 1996: Clim. Dyn. 557-572
!     Levkov et al., 1992: Beitr. Phys. Atm. 35-58.          (ice phase)
!     Beheng, 1994: Atmos. Res. 193-206.                    (warm phase)
!     Lenderink et al., 1998; KNMI-REPORT NO. 98-13       (condensation)
!     Tompkins 2002, J. Atmos. Sci.                        (cloud cover)
!
!     Authors.
!     -------
!     M.Esch        MPI-Hamburg  1999
!     G.Lenderink   KNMI, de Bilt 1998
!     U.Lohmann     MPI-Hamburg  1995
!
!     Modifications.
!     --------------
!     E.Roeckner    MPI-Hamburg  2000
!     A.Tompkins    MPI-Hamburg  2000
!     U.Schlese     MPI-Hamburg  2003
!     U.Lohmann     Dalhousie University 2002-2006: Prognostic CDNC/IC scheme
!     P.Stier       MPI-Hamburg          2002-2006: Prognostic CDNC/IC scheme
!                                                   Scavenging parameters
!     J.Zhang       Dalhousie University      2004: Prognostic CDNC/IC scheme
!     S. Ferrachat  ETH Zuerich  2008: complete re-writing to allow vectorization 
!                                      cleanup of the code 
!
    USE messy_cloud_ori,  ONLY : cqtmin, tmelt, cvtfall, crhosno, cn0s    &
                               , cthomi, csecfrl, ncctop, cvarmin         &
                               , cbeta_pq, cbeta_pq_max, nbetaq, cbetaqs  &
                               , rbetak, nbetax, tbetai0, tbetai1         &
                               , clmax, clmin, jbmin, jbmax, lonacc       &
                               , ccraut, crhoi, ccsaut                    &
                               , cbeta_cs, LOOKUPOVERFLOW
    USE messy_main_tools, ONLY : jptlucu1, jptlucu2,                      &
                                 tlucua, tlucuaw, tlucub

    USE messy_main_constants_mem, ONLY: ceffmin, ceffmax, ccwmin, &
                                        rd, rv, vtmpc1, vtmpc2,   &
                                        cpd => cp_air,  &
                                        tmelt, rhoh2o => rho_h2o
   USE messy_cloud_ice_BN09, ONLY: ice_activate

  IMPLICIT NONE
 
  INTEGER :: krow, ktdia, kproma, kbdim, klev, klevp1, ktrac 

  REAL(dp):: paphm1(kbdim,klevp1), pvervel(kbdim,klev)  &
            ,papm1(kbdim,klev)   , pqm1(kbdim,klev)     &
            ,papp1(kbdim,klev)   , ptm1(kbdim,klev)     &
            ,ptvm1(kbdim,klev)   , pxlm1(kbdim,klev)    &
            ,pxim1(kbdim,klev)   , pxtec(kbdim,klev)    &
            ,pqtec(kbdim,klev)   , pxvar(kbdim,klev)    &
            ,pxskew(kbdim,klev)  , pbetaa(kbdim,klev)   &
            ,pbetab(kbdim,klev)  , pvdiffp(kbdim,klev)  &
            ,phmixtau(kbdim,klev), pvmixtau(kbdim,klev) &
            ,pgeo(kbdim,klev)    , pbetass(kbdim,klev)  &
            ,pxlvi(kbdim)        ,pxivi(kbdim)          & 
            ,paclc(kbdim,klev)   ,paclcac(kbdim,klev)   &
            ,pacdnc(kbdim,klev)  ,prelhum(kbdim,klev)   &
            ,paclcov(kbdim)      ,paprl(kbdim)          &
            ,pqvi(kbdim)         ,pssfl(kbdim)          &
            ,ptte(kbdim,klev)    ,pqte(kbdim,klev)      & 
            ,pxlte(kbdim,klev)   ,pxite(kbdim,klev)     &
            ,paprs(kbdim)        ,prsfl(kbdim)

!---Included for in-cloud scavenging (Philip Stier, 28/03/01):----------
  REAL(dp) ::    pxtm1(kbdim,klev,ktrac), pxtte(kbdim,klev,ktrac) 
!---End Included for scavenging-----------------------------------------
  CHARACTER(LEN=32) :: status_string
  LOGICAL  :: lcover
  REAL(dp) :: slm(kbdim), glac(kbdim)

  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev)         :: plwc,     piwc
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev), TARGET :: pfrain,   pfsnow
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev)         :: pfrain_no,pfsnow_no
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev), TARGET :: prevap,   pssubl
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev), TARGET :: prate_r,  prate_s 
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev)         :: pr_cover
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev)         :: pimelt,   pisedi
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev)         :: pcond 

  ! Local variables
  REAL(dp), POINTER, DIMENSION(:) :: zrfl, zsfl, zsub, zevp, zrpr
!
!   Temporary arrays
!

  REAL(dp)   :: zclcpre_2d(kbdim,klev)    ,zclcpre(kbdim)              &
           ,zcnd(kbdim)         ,zdep(kbdim)                           &
           ,                     zxievap(kbdim)     ,zxlevap(kbdim)    &
           ,zfrl(kbdim,klev)    ,zimlt(kbdim)       ,zsmlt(kbdim)      &
           ,                     zspr(kbdim)                           &
           ,zxlte(kbdim)        ,zxite(kbdim)       ,zxiflux(kbdim)    &
           ,zsacl(kbdim)        ,zxlte2(kbdim)      ,zxite2(kbdim)     &
           ,zlsdcp(kbdim)       ,zlvdcp(kbdim)      ,zximlt(kbdim)     &
           ,ztp1tmp(kbdim)      ,zqp1tmp(kbdim)     ,zxisub(kbdim)     &
           ,zxlb(kbdim)         ,zxib(kbdim)        ,zxibold_1d(kbdim) &
           ,zrho(kbdim,klev)    ,zclcov(kbdim)      ,zclcaux(kbdim)    &
           ,zrho_rcp(kbdim,klev)                                       &
           ,zqvi(kbdim)         ,zxlvi(kbdim)       ,zxivi(kbdim)      &
           ,zbetaqt(kbdim)      ,zwide(kbdim)                          &
           ,zbetacl(kbdim)      ,zturbvar(kbdim)    ,zturbskew(kbdim)  &
           ,zconvvar(kbdim)     ,zconvskew(kbdim)   ,zvartg(kbdim)     &
           ,zmicroskew(kbdim)   ,zgenti(kbdim)      ,zgentl(kbdim)     &
           ,zxvarte(kbdim)      ,zxskewte(kbdim)                       &
           ,zgeoh(kbdim,klevp1) 
!
  INTEGER   knvb(kbdim)

!--- Included for dust emissions (Philip Stier 10/01/02)-----------------
  INTEGER :: jrow
!--- End Included for dust emissions ------------------------------------

!---Included for in-cloud scavenging (Philip Stier, 28/03/01):----------
  REAL(dp)    :: zmratepr(kbdim,klev), & ! Rain formation rate in cloudy part
                                     ! of the grid box [kg/kg]
             zmrateps(kbdim,klev), & ! Ice  formation rate in cloudy part
                                     ! of the grid box  [kg/kg]
             zfrain(kbdim,klev),   & ! Rain flux before evaporation
                                     ! [kg/m2/s]
             zfsnow(kbdim,klev),   & ! Snow flux before sublimation
                                     ! [kg/m2/s]
             zfevapr(kbdim,klev),  & ! Evaporation of rain [kg/m2/s]
             zfsubls(kbdim,klev),  & ! Sublimation of snow [kg/m2/s]
             zmlwc(kbdim,klev),    & ! In-cloud liquid water mass mixing
                                     ! ratio before rain formation [kg/kg]
             zmiwc(kbdim,klev),    & ! In-cloud ice mass mixing ratio
                                     ! before snow formation [kg/kg]
             zmsnowacl(kbdim,klev)   ! Accretion rate of snow with cloud
                                     ! droplets in cloudy part of the 
                                     ! grid box  [kg/kg]
!---End Included for scavenging-----------------------------------------

!--- Included for prognostic CDNC/IC scheme ----------------------------

  REAL(dp) :: pcdncact(kbdim,klev),     ptkem1(kbdim,klev),             &
             pcvcbot(kbdim),           pwcape(kbdim),                  &
             pxtecl(kbdim,klev),       pxteci(kbdim,klev),             &
             pxtecnl(kbdim,klev),       pxtecni(kbdim,klev)

  INTEGER ::         &
             itop(kbdim,klev), ibas(kbdim,klev), &
             icl_minusbas(kbdim,klev), icl_minustop(kbdim,klev),    &
             iclbas(kbdim,klev), itm1_look(kbdim,klev), it1_1d(kbdim) 

  REAL(dp)    :: zrprn(kbdim),             zsprn(kbdim,klev),              &
             zsacln(kbdim),            zfrln(kbdim),                   &
             zcdnc(kbdim,klev),     & !
             zcdnc_burden(kbdim),   &
             zqlnuc(kbdim,klev),    & !
             zqlnuccvh(kbdim,klev), & ! Nucleated CDNC from convective detrainment
             zqlnuccv(kbdim,klev), & 
             zqlnuc_bas(kbdim,klev), &
             zcdnc_bas(kbdim,klev)

  !---------------------------------------------------------------
  !--- Variables used by BN09 ice scheme (in ice_activate subroutine)
  !--- or needed for the implementation of BN09
  !---------------------------------------------------------------
  REAL(dp) :: sigw(kbdim,klev),       & !width of the distribution of updraft velocity [m/s]    
              ndrop(kbdim,klev),      & !total cloud droplet number concentration [m-3]
              dsulf(kbdim,klev),      & !diameter [m] of sulfate
              ndust(kbdim,klev,2),    & !number concentration [m-3] of dust in 2 modes (ai, ci)
              ddust(kbdim,klev,2),    & !diameter [m] of dust in 2 modes (ai, ci)
              sigdust(2),             & !aerosol lognormal distribution of ai and ci  
              norg(kbdim,klev),       & !number concentration [m-3] of organics
              dorg(kbdim,klev),       & !diameter [m] of organics
              sigorg,                 & !aerosol lognormal distribution of ki  
              nsoot(kbdim,klev),      & !number concentration [m-3] of soot
              dsoot(kbdim,klev),      & !diameter [m] of soot
              sigsoot,                & !aerosol lognormal distribution of ki  
              nbio(kbdim,klev),       & !number concentration [m-3] of biological aerosols
              dbio(kbdim,klev),       & !diameter [m] of biological aerosols
              sigbio,                 & !aerosol lognormal distribution of ki  
              nsolo(kbdim,klev,3),    & !number concentration [m-3] of soluble organics in 3 modes (ks, as, cs)
              dsolo(kbdim,klev,3),    & !diameter [m] of soluble OC in 3 modes
              sigsolo(3),             & !aerosol lognormal distribution of ks, as, cs
              smaxice,                & !maximum supersaturation with respect to ice
              sc_ice,                 & !characteristic freezing point of the aerosol population 
              nlim,                   & !limiting IN concentration [m-3]
              nhet,                   & !ice crystal concentration by het freezing [m-3]
              dice,                   & !diameter of ice crystals [m]
              nice,                   & !nucleated ice crystal number concentration [m-3]
              sigwparc,               & !corrected vertical velocity accounting for the effect of preexisting ice crystals [m/s]
              rdryki(kbdim,klev),     & !dry radius of ki mode
              rdryai(kbdim,klev),     & !dry radius of ai mode
              rdryci(kbdim,klev),     & !dry radius of ci mode 
              rdryks(kbdim,klev),     & !dry radius of ks mode
              rdryas(kbdim,klev),     & !dry radius of as mode
              rdrycs(kbdim,klev),     & !dry radius of cs mode 
              sigmaki,                & !aerosol lognormal distribution of ki 
              sigmaai,                & !aerosol lognormal distribution of ai 
              sigmaci,                & !aerosol lognormal distribution of ci 
              sigmaks,                & !aerosol lognormal distribution of ks 
              sigmaas,                & !aerosol lognormal distribution of as 
              sigmacs                   !aerosol lognormal distribution of cs 
              
  !Input variables 
  REAL(dp) :: sigwBN(kbdim,klev), ndropBN(kbdim,klev), dsulfBN(kbdim,klev),                                   &
              ndust_aiBN(kbdim,klev), ddust_aiBN(kbdim,klev), ndust_ciBN(kbdim,klev), ddust_ciBN(kbdim,klev), &
              norgBN(kbdim,klev), dorgBN(kbdim,klev), nsootBN(kbdim,klev), dsootBN(kbdim,klev),               &
              nsolo_ksBN(kbdim,klev), dsolo_ksBN(kbdim,klev), nsolo_asBN(kbdim,klev), dsolo_asBN(kbdim,klev), &
              nsolo_csBN(kbdim,klev), dsolo_csBN(kbdim,klev)
  !Output variables 
  REAL(dp) :: smaxice_cirrusBN(kbdim,klev), sc_ice_cirrusBN(kbdim,klev), nlim_cirrusBN(kbdim,klev),  &
              nhet_cirrusBN(kbdim,klev), nice_cirrusBN(kbdim,klev), dice_cirrusBN(kbdim,klev),       &
              sigwpre_cirrusBN(kbdim,klev),                                                          &
              smaxice_immBN(kbdim,klev), sc_ice_immBN(kbdim,klev), nice_immBN(kbdim,klev),           &
              dice_immBN(kbdim,klev), sigwpre_immBN(kbdim,klev),                                     &     
              newIC_cirri(kbdim,klev), newICR_cirri(kbdim,klev), newIC_cnt_therm(kbdim,klev),        &
              newIC_imm(kbdim,klev), newIC_mix(kbdim,klev), newIC_mix_fin(kbdim,klev),               &
              w_sub(kbdim,klev), w_ave(kbdim,klev), w_grid(kbdim,klev)
  INTEGER  :: I, J
  !---------------------------------------------------------------

  !REAL(dp)   zw_strat(kbdim,klev) ,  & ! stratiform updraft velocity, large-scale+TKE (>0.0) [m s-1]
  !           zw_conv(kbdim,klev)      ! convective updraft velocity, large-scale+CAPE (>0.0) [m s-1]

 ! Temporary fields
  INTEGER  :: jl, jk, jkk, jmod,  &
!SFcray
              iqidx, ixidx
!SFcray

  REAL(dp) :: zesw,    &                  
              zf1, &        
              zrieff, &
              zradl,  &
              zdisp, zdw0, ztmst, ztmst_rcp, zcons1, zcons2, zcons3, zcons4, zcons5, zdtime, zdtime_rcp, ztmstdt, zdt, g_rcp,    &
              zpirho, zpirho_rcp, &
              zb2, &
              zdv,  &
              zgtp, & 
              zexm1_1, zexp_1, &
              zexm1_2, zexp_2, &
              zrih, &
              ztte, zomega, &
              ztc, zvth, zfuchs, zfre, zre, &
              zfracdusol, zfracduinsolai, zfracduinsolci, zfracbcsol, zfracbcinsol, &
              zfrzcntdu, zfrzcntbc, zfrzcnt, znaimmdu, znaimmbc, zfrzimm, zfrzthermo, &
              zqsm1, zqst1, zdqsdtime, zdeltatemp, zf2, zcap

  REAL(dp)    :: zicnc(kbdim,klev),     & !
             zicncq(kbdim,klev),    &
             zicnc_burden(kbdim),   & !
             zascs(kbdim,klev),     & !
             zrid_2d(kbdim,klev),   &
             zqsi_2d(kbdim,klev),   &
             zninucl(kbdim, klev),  & ! number conc. of newly nucleated IC
             zqinucl(kbdim, klev),  & ! mixing ratio of newly nucleated IC
             zri(kbdim, klev),      & ! size of newly nucleated IC
             zreffct(kbdim),        &  ! cloud top effective radius
             ztau1i(kbdim),         &  ! ice cloud optical depth - visible wavelength
             znidetr(kbdim, klev)     ! IC number from detrainment

 !SFnew arrays
   !SF 2D
 INTEGER    :: itm1p1_look(kbdim,klev), it_1d(kbdim)

 REAL(dp)   :: zdz_2d(kbdim,klev), zdp_2d(kbdim,klev), zdpg_2d(kbdim,klev), zaaa_2d(kbdim,klev), zviscos_2d(kbdim,klev), &
               zqswp1_2d(kbdim,klev), zqsip1_2d(kbdim,klev), zqsw_2d(kbdim,klev), &
               zastbstw(kbdim,klev), zastbsti(kbdim,klev), zmmean_2d(kbdim,klev), &
               zxifallmc_2d(kbdim,klev), zalfased_2d(kbdim,klev), zbetased_2d(kbdim,klev), &
               zxifallnc_2d(kbdim,klev), &
               zxised_2d(kbdim,klev),zicesed_2d(kbdim,klev), &
             zesw_2d(kbdim,klev),   & !
             zesi_2d(kbdim,klev),   & !
             zsusatw_2d(kbdim,klev),& ! supersat. with respect to water
             zsusatw_evap(kbdim,klev),& 
             zsusatix_2d(kbdim,klev),& 
             zvervx_2d(kbdim,klev),  &     ! updraft [cm/s]
             zicesub(kbdim,klev), zqrho_2d(kbdim,klev)
               
   !SF 1D
 REAL(dp)   :: zscnc_1d(kbdim), zdplanar_1d(kbdim), zusnow_1d(kbdim), zsacl2in_1d(kbdim), &
               zxsp2_1d(kbdim), zxsp_1d(kbdim), zstcrit_1d(kbdim),  &
               zstokes_1d(kbdim), zudrop_1d(kbdim), zrey_1d(kbdim), &
               zrac1_1d(kbdim), zrac2_1d(kbdim), zrautself_1d(kbdim), zxsp1_1d(kbdim), &
               zraut_1d(kbdim), zsaut_1d(kbdim), zsaci2_1d(kbdim), zsacl2_1d(kbdim), &
               zris_1d(kbdim), zcolleffi_1d(kbdim), zc1_1d(kbdim), zxlp1_1d(kbdim), &
               zreffl_1d(kbdim), zdxlcor_1d(kbdim), &
               zdxicor_1d(kbdim), zsecprod_1d(kbdim), &
               zcsacl_1d(kbdim), zpretot_1d(kbdim), zpredel_1d(kbdim), zpresum_1d(kbdim), &
               zzdrr_1d(kbdim), zzdrs_1d(kbdim), zxidtstar_1d(kbdim), &
               zlc_1d(kbdim), zxldtstar_1d(kbdim), zxlm1evp_1d(kbdim), zxim1evp_1d(kbdim), &
               zxldt_1d(kbdim), zxidt_1d(kbdim), zqsm1_1d(kbdim), zqp1_1d(kbdim), zdtdt_1d(kbdim), &
               zdqsat_1d(kbdim), ztp1_1d(kbdim), zdqsdt_1d(kbdim), zqst1_1d(kbdim), &
               zqvdt_1d(kbdim), &
               zxip1_1d(kbdim), zcorw_1d(kbdim), zesw_1d(kbdim), zoversatw_1d(kbdim), &
               zqsp1tmpw_1d(kbdim), zqsp1tmp_1d(kbdim), zcor_1d(kbdim), zrhtest_1d(kbdim), &
               zoversat_1d(kbdim), zlcdqsdt_1d(kbdim), zclcstar_1d(kbdim), zxrp1_1d(kbdim), &
               zauloc_1d(kbdim), zqsp1tmphet_1d(kbdim), zqcon_1d(kbdim), zrelhum_1d(kbdim), &
               iqidx_1d(kbdim), zbap1_1d(kbdim), zgent_1d(kbdim), ixidx_1d(kbdim), &
               zqtau_1d(kbdim), zxilb_1d(kbdim), zbbap1_1d(kbdim), zbqp1_1d(kbdim), &
               zqcdif_1d(kbdim), zlucuawp1_1d(kbdim), zlucuap1_1d(kbdim), zes_1d(kbdim), &
               zlucub_1d(kbdim), zlucuaw_1d(kbdim), zlucua_1d(kbdim), zicncp1_1d(kbdim), &
               zetaair_1d(kbdim), zkair_1d(kbdim), &
               zdfarbcki_1d(kbdim), zdfarduai_1d(kbdim), zdfarduci_1d(kbdim), &
               zknbcki_1d(kbdim), zknduai_1d(kbdim), zknduci_1d(kbdim), &
               zftbcki_1d(kbdim), zftduai_1d(kbdim), zftduci_1d(kbdim), &
               zvervmax_1d(kbdim), zrice_1d(kbdim), zeta_1d(kbdim), zdv_1d(kbdim)
               
  REAL(dp)  :: ztmp(kbdim,klev), ztmp1(kbdim,klev), ztmp2(kbdim,klev), &
               ztmp1_1d(kbdim), ztmp2_1d(kbdim), ztmp3_1d(kbdim), &
               ztmp4_1d(kbdim), ztmp5_1d(kbdim)

  LOGICAL :: ll_look(kbdim,klev), &
             ll_cv(kbdim,klev), &
             ll_ice(kbdim,klev), ll1(kbdim,klev), &
             ll1_2d(kbdim,klev), ll2_2d(kbdim,klev), ll3_2d(kbdim,klev), lo_1d(kbdim), &
             lo2_1d(kbdim), locc_1d(kbdim), &
             ll1_1d(kbdim), ll2_1d(kbdim), ll3_1d(kbdim), ll4_1d(kbdim), & 
             ll5_1d(kbdim), ll6_1d(kbdim), ll7_1d(kbdim), ll8_1d(kbdim) 
 !SF end new arrays

  !--- Included for cirrus scheme -------------------------------------------
  REAL(dp) :: znicex(kbdim, klev)     ! 
  REAL(dp) :: zapnx(kbdim,nfrzmod) ! aerosol number available for homogeneous freezing [1/cm3]
  REAL(dp) :: zaprx(kbdim,nfrzmod) ! radius of aerosols available for homogeneous freezing  [cm]
  REAL(dp) :: zapsigx(kbdim,nfrzmod)! std. dev. of aerosols available for homogeneous freezing 
  REAL(dp) :: zap(kbdim,klev)

  LOGICAL, PARAMETER :: nosize = .false. ! .true. ---> aerosol size effects
                                         ! are ignored for homogeneous
                                         ! freezing

    ! !!! lhet = .true. NOT YET IMPLEMENTED !!!

  LOGICAL, PARAMETER :: lhet = .false.   ! .true. ---> heterogeneous
                                         ! freezing of aerosol particles
                                         ! below 235K is considered
  !--- End included for cirrus scheme -------------------------------------------

  REAL(dp) :: zxifluxn(kbdim)

!--- End Included for CDNC -------------------------------------------
  REAL(dp) :: zrwetki_2d(kbdim,klev), zrwetai_2d(kbdim,klev), zrwetci_2d(kbdim,klev)
!--- Ulrike: for additional disgnostics ---------------------------------------
  REAL(dp) :: zrdry(kbdim,klev,nmod) ! dry radius for each mode
                                    !@@@ currently re-stored, check if avoidable
!--- end Ulrike: for additional diagnostics -----------------------------------
!
! Executable statements
!
  jrow = krow

!-- 0. Initialisations:

  lookupoverflow = .FALSE.

  ! initializing new values for channel objects
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
  reffi(:,:,jrow)         = ceffmax
  reffl(:,:,jrow)         = 40._dp
  cdnc(:,:,jrow)          = 0._dp
  icnc(:,:,jrow)          = 0._dp

  w_sub(:,:)            = 0._dp
  w_ave(:,:)            = 0._dp
  w_grid(:,:)           = 0._dp
  sigwBN(:,:)           = 0._dp
  ndropBN(:,:)          = 0._dp
  dsulfBN(:,:)          = 0._dp
  ndust_aiBN(:,:)       = 0._dp
  ddust_aiBN(:,:)       = 0._dp
  ndust_ciBN(:,:)       = 0._dp 
  ddust_ciBN(:,:)       = 0._dp
  norgBN(:,:)           = 0._dp
  dorgBN(:,:)           = 0._dp
  nsootBN(:,:)          = 0._dp
  dsootBN(:,:)          = 0._dp
  nsolo_ksBN(:,:)       = 0._dp
  dsolo_ksBN(:,:)       = 0._dp
  nsolo_asBN(:,:)       = 0._dp
  dsolo_asBN(:,:)       = 0._dp
  nsolo_csBN(:,:)       = 0._dp 
  dsolo_csBN(:,:)       = 0._dp
  smaxice_cirrusBN(:,:) = 0._dp
  sc_ice_cirrusBN(:,:)  = 0._dp
  nlim_cirrusBN(:,:)    = 0._dp
  nhet_cirrusBN(:,:)    = 0._dp
  nice_cirrusBN(:,:)    = 0._dp
  dice_cirrusBN(:,:)    = 0._dp
  sigwpre_cirrusBN(:,:) = 0._dp 
  smaxice_immBN(:,:)    = 0._dp
  sc_ice_immBN(:,:)     = 0._dp 
  nice_immBN(:,:)       = 0._dp
  dice_immBN(:,:)       = 0._dp
  sigwpre_immBN(:,:)    = 0._dp
  newIC_cirri(:,:)      = 0._dp
  newICR_cirri(:,:)     = 0._dp
  newIC_cnt_therm(:,:)  = 0._dp
  newIC_imm(:,:)        = 0._dp
  newIC_mix(:,:)        = 0._dp
  newIC_mix_fin(:,:)    = 0._dp

!---Included for in-cloud scavenging (Philip Stier, 28/03/01):----------
  zmratepr(:,:)  = 0._dp
  zmrateps(:,:)  = 0._dp
  zfrain(:,:)    = 0._dp
  zfsnow(:,:)    = 0._dp
  zfevapr(:,:)   = 0._dp
  zfsubls(:,:)   = 0._dp
  zmlwc(:,:)     = 0._dp
  zmiwc(:,:)     = 0._dp
  zmsnowacl(:,:) = 0._dp
!---End Included for scavenging-----------------------------------------

!--- Included for prognostic CDNC/IC scheme ----------------------------
  zqinucl(:,:) = 0._dp
  zninucl(:,:) = 0._dp
  zreffct(:)   = 0._dp
  ztau1i(:)    = 0._dp
!--- End Included for CDNC --------------------------------------------

!--- Included for prognostic CDNC/IC scheme ----------------------------

!cyjstart----------------------------
  zxifluxn(:)       = 0.0_dp
  zmmean_2d(:,:)    = zmi
  zxifallmc_2d(:,:) = 0.0_dp
  zxifallnc_2d(:,:) = 0.0_dp
  zxised_2d(:,:)    = 0.0_dp
  zicesed_2d(:,:)   = 0.0_dp
!cyjend----------------------------------------

  zcdnc_burden(1:kproma) = 0.0_dp
  zcdnc_bas(1:kproma,:)  = 0.0_dp 
  zqlnuc(1:kproma,:)     = 0.0_dp
  zqlnuc_bas(1:kproma,:) = 0.0_dp
  zicnc_burden(1:kproma) = 0.0_dp
  zri(1:kproma,:)        = 1.e-6_dp

!--- End included for CDNC/IC scheme -----------------------------------

!---Computational constants:
  zdisp      = EXP(0.5_dp * zsigmaw**2)
  zdw0       = 10.e-6_dp*zdisp       ! modal diameter times dispersion parameter
  zdtime     = ztmst / 2._dp
  zdtime_rcp = 1._dp / zdtime
  ztmst_rcp  = 1._dp / ztmst

  ztmstdt    = ztmst * zdtime
  zdt        = ztmst_rcp * zdtime 
  g_rcp      = 1._dp / g
  zcons1     = cpd*vtmpc2
  zcons2     = ztmst_rcp * g_rcp
  zexm1_1    = 2.47_dp-1.0_dp
  zexp_1     = -1._dp / zexm1_1
  zexm1_2    = 4.7_dp-1.0_dp
  zexp_2     = -1._dp / zexm1_2
  zpirho     = api*rhoh2o
  zpirho_rcp = 1._dp / zpirho
  zcap       = 2._dp / api

  zcons3 = 1._dp / ( api*crhosno*cn0s*cvtfall**(1._dp/1.16_dp) )**0.25_dp
  zcons4 = 1._dp / ( api*crhosno*cn0s )**0.8125_dp
  zcons5 = 1._dp / ( api*crhosno*cn0s )**0.875_dp

!SF store wet radius for freezing calculations (6.):

  zrwetki_2d(1:kproma,:) = mode(iaiti)%wetrad(1:kproma,:)
  zrwetai_2d(1:kproma,:) = mode(iacci)%wetrad(1:kproma,:)
  zrwetci_2d(1:kproma,:) = mode(icoai)%wetrad(1:kproma,:)

! Ulrike: for thermophoresis
  DO jmod = 1,nmod
    zrdry(1:kproma,:,jmod) = mode(jmod)%dryrad(1:kproma,:)
  ENDDO
! end Ulrike: for thermophoresis

!--- Get several utility variables:
  CALL get_util_var(kproma, kbdim, ktdia, klev, klevp1,           &
                    paphm1, pgeo, papm1, ptm1, zgeoh,             &
                    zdp_2d, zdpg_2d, zdz_2d, zaaa_2d, zviscos_2d)


!
!     ------------------------------------------------------------------
!
!       1.   Top boundary conditions, air density
!
!       1.1   Set to zero precipitation fluxes etc.
!
  zclcpre(:) = 0.0_dp   ! fraction of grid box covered by precip
  zxiflux(:) = 0.0_dp   ! flux of ice crystals falling into the grid box from above
!
!       1.2   Air density


  DO 122 jk = ktdia,klev

     zrho(:,jk)     = papm1(:,jk)/(rd*ptvm1(:,jk))
     zrho_rcp(:,jk) = 1._dp / zrho(:,jk)
     zqrho_2d(:,jk) = 1.3_dp * zrho_rcp(:,jk)
     zfrl(:,jk)     = 0.0_dp

!-----------------Added by Junhua Zhang for CONV Micro--------------------
     IF (ncvmicro == 0) THEN

        pxtec(:,jk)  = MAX(pxtec(:,jk),0.0_dp)

        ll1_1d(:) = (pxtec(:,jk) > 0.0_dp)

        ztmp1_1d(:)         = twc_conv(1:kproma,jk,jrow) + ztmstdt*pxtec(:,jk)
        twc_conv(1:kproma,jk,jrow) = MERGE(ztmp1_1d(:), twc_conv(1:kproma,jk,jrow), ll1_1d(:))

        ztmp1_1d(:)          = conv_time(1:kproma,jk,jrow) + zdtime
        conv_time(1:kproma,jk,jrow) = MERGE(ztmp1_1d(:), conv_time(1:kproma,jk,jrow), ll1_1d(:))
 
     ELSE !SF ncvmicro .ne. 0

           pxtecl(:,jk)  = MAX(pxtecl(:,jk),  0.0_dp)
           pxteci(:,jk)  = MAX(pxteci(:,jk),  0.0_dp)
           pxtecnl(:,jk) = MAX(pxtecnl(:,jk), 0.0_dp)
           pxtecni(:,jk) = MAX(pxtecni(:,jk), 0.0_dp)

           ll1_1d(:) = (pxtecl(:,jk) > 0._dp) .AND. &
                       (ptm1(:,jk)   < cthomi)

           ztmp1_1d(:)  = pxteci(:,jk) + pxtecl(:,jk)
           pxteci(:,jk) = MERGE(ztmp1_1d(:), pxteci(:,jk), ll1_1d(:))

           pxtecl(:,jk) = MERGE(0._dp, pxtecl(:,jk), ll1_1d(:))

           ztmp1_1d(:)   = pxtecni(:,jk) + pxtecnl(:,jk)
           pxtecni(:,jk) = MERGE(ztmp1_1d(:), pxtecni(:,jk), ll1_1d(:))

           pxtecnl(:,jk) = MERGE(0._dp, pxtecnl(:,jk), ll1_1d(:))

           zfrl(:,jk) = MERGE(pxtecl(:,jk), zfrl(:,jk), ll1_1d(:))

           ll1_1d(:) = (pxteci(:,jk) > 0.0_dp) .AND. &
                       (ptm1(:,jk)   > tmelt)

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

     ENDIF  !SF end ncvmicro

!------------------------------End----------------------------------------

     pqtec(:,jk)  = MAX(pqtec(:,jk),0.0_dp)

!--- Included for prognostic CDNC/IC scheme ----------------------------

     ! calculate cloud droplet number concentration 
!qqq+ include tendency for correct start value (see Lohman07)?
!!$     zcdnc(:,jk) = zrho(:,jk)*(pxtm1(:,jk,1)/m_air * 1000._dp)
     zcdnc(:,jk) = zrho(:,jk)*((pxtm1(:,jk,1)+pxtte(:,jk,1)*ztmst) &
          /m_air * 1000._dp)
!qqq-
     zcdnc(:,jk) = MAX(zcdnc(:,jk),cqtmin)

     ! calculate ice crystal number
!qqq+ see above
!!$     zicnc(:,jk)  = zrho(:,jk)*(pxtm1(:,jk,2)/m_air * 1000._dp)
     zicnc(:,jk)  = zrho(:,jk)*((pxtm1(:,jk,2)+pxtte(:,jk,2)*ztmst) &
          /m_air * 1000._dp)
!qqq-
     zicnc(:,jk)  = MAX(zicnc(:,jk),cqtmin)
     zicncq(:,jk) = zicnc(:,jk)

     ! END changed S.K.
     ! end corinna
122 END DO
!
!       1.3   Geopotential at half levels
!

  CALL get_cloud_bounds(kproma, kbdim, ktdia, klev, paclc, itop, ibas, icl_minustop, icl_minusbas)

  ztmp(:,:)      = 1000._dp*ptm1(:,:)
  itm1_look(:,:) = NINT(ztmp(:,:))

  ll_look(:,:) = (itm1_look(:,:)<jptlucu1 .OR. itm1_look(:,:)>jptlucu2)
  IF (ANY(ll_look(1:kproma,ktdia:klev))) lookupoverflow = .TRUE.


  itm1_look(:,:) = MAX(MIN(itm1_look(:,:),jptlucu2),jptlucu1)

! for later use in 5.:
  itm1p1_look(:,:) = itm1_look(:,:) + 1
  itm1p1_look(:,:) = MAX(MIN(itm1p1_look(:,:),jptlucu2),jptlucu1)

  DO jk=klev,ktdia,-1
     DO jl=1,kproma
!SF water:
        zesw              = tlucuaw(itm1_look(jl,jk))/papm1(jl,jk)
        zesw              = MIN(zesw,0.5_dp)
        zqsw_2d(jl,jk)    = zesw/(1._dp-vtmpc1*zesw)
        zsusatw_2d(jl,jk) = pqm1(jl,jk)/zqsw_2d(jl,jk)-1.0_dp
        zsusatw_2d(jl,jk) = MAX(zsusatw_2d(jl,jk),0.0_dp)
        !SF for evaporation of rain in 3.3:
        zsusatw_evap(jl,jk) = pqm1(jl,jk)/zqsw_2d(jl,jk)-1.0_dp
        zsusatw_evap(jl,jk) = MIN(zsusatw_evap(jl,jk), 0._dp)

        !--- Saturation water vapour pressure:
        zesw_2d(jl,jk) = zesw*papm1(jl,jk)*rv/rd

        !SF for later use in 5.:
        zqswp1_2d(jl,jk) = tlucuaw(itm1p1_look(jl,jk))/papm1(jl,jk)
        zqswp1_2d(jl,jk) = MIN(zqswp1_2d(jl,jk),0.5_dp)
        zqswp1_2d(jl,jk) = zqswp1_2d(jl,jk)/(1._dp-vtmpc1*zqswp1_2d(jl,jk))   
        !SFend for later use in 5.

!SF ice:

        zqsi_2d(jl,jk)   = tlucua(itm1_look(jl,jk))/papm1(jl,jk)
        zqsi_2d(jl,jk)   = MIN(zqsi_2d(jl,jk),0.5_dp)
        zesi_2d(jl,jk)   = zqsi_2d(jl,jk)*papm1(jl,jk)*rv/rd
        zqsi_2d(jl,jk)   = zqsi_2d(jl,jk)/(1._dp-vtmpc1*zqsi_2d(jl,jk))
        sice(jl,jk,jrow) = pqm1(jl,jk)/zqsi_2d(jl,jk)-1.0_dp
        sice(jl,jk,jrow) = MAX(sice(jl,jk,jrow), 0._dp)

        !SF for sublimation of snow and ice in 3.2:
        zicesub(jl,jk) = pqm1(jl,jk)/zqsi_2d(jl,jk)-1.0_dp
        zicesub(jl,jk) = MIN(zicesub(jl,jk), 0._dp) 

        !SF for later use in 5.:
        zqsip1_2d(jl,jk) = tlucua(itm1p1_look(jl,jk))/papm1(jl,jk)
        zqsip1_2d(jl,jk) = MIN(zqsip1_2d(jl,jk),0.5_dp)
        zqsip1_2d(jl,jk) = zqsip1_2d(jl,jk)/(1._dp-vtmpc1*zqsip1_2d(jl,jk))
        !SFend for later use in 5.
 
     END DO !jl
  END DO !jk

  !--- Store supersaturation with respect to water in stream:
  swat(1:kproma,:,jrow)  = zsusatw_2d(1:kproma,:)

!SF Set some utility variables

  zastbstw(:,:) = alv * (alv/(rv*ptm1(:,:)) - 1.0_dp) / (zka*ptm1(:,:))  & !SF zast for water
                + rv*ptm1(:,:)/(2.21_dp/papm1(:,:)*zesw_2d(:,:))           !SF zbst for water

  ztmp1(:,:) = 4.1867e-3_dp*(5.69_dp + 0.017_dp*(ptm1(:,:)-tmelt)) !SF zkair

  zastbsti(:,:) = als * (als/(rv*ptm1(:,:)) - 1.0_dp) / (ztmp1(:,:)*ptm1(:,:))  & !SF zast for ice
                + rv*ptm1(:,:)/(2.21_dp/papm1(:,:)*zesi_2d(:,:))                     !SF zbst for ice

!--- Aerosol activation (so far only the Lin and Leaitch scheme is tested)


!--- Convert the aerosol activation into the number of newly formed cloud droplets

  ll1_2d(:,:) = ( ibas(:,:) > 0 )                                .AND.  &  !cloud base level
                ( (zcdnc(:,:)      <= cdncmin)             .OR.         &
                  (pxlte(:,:)      >  0._dp)               .OR.         &
                  (paclc(:,:)      >  cloud_tm1(1:kproma,:,jrow)) .OR.         &
                  (zsusatw_2d(:,:) >  zeps)                     )

  !SF: first computes newly formed cloud droplets at cloud bases:
  ztmp(:,:)      = pcdncact(:,:) - zcdnc(:,:) 
  ztmp(:,:)      = MAX(0._dp, ztmp(:,:))
  zqlnuc(:,:)    = MERGE(ztmp(:,:), 0._dp, ll1_2d(:,:))

  ztmp(:,:)      = zcdnc(:,:) + zqlnuc(:,:)
  zcdnc(:,:)     = MERGE(ztmp(:,:), zcdnc(:,:), ll1_2d(:,:))
  qnuc(1:kproma,:,jrow) = qnuc(1:kproma,:,jrow) + zdt*zqlnuc(1:kproma,:)
 
  !   then computes newly formed cloud droplets above cloud base
  !    by assuming that the number of nucleated cloud droplets is constant above cloud base
  !    (adiabatic parcel theory)

  DO jk=ktdia,klev
     DO jl=1,kproma
        jkk = MAX(jk,icl_minusbas(jl,jk))  !sets the level index to either relevant cloud base or itself
        zqlnuc_bas(jl,jk) = zqlnuc(jl,jkk) !holds in each cloud level the number of
                                           !newly formed cloud droplets at the base

        zcdnc_bas(jl,jk) = zcdnc(jl,jkk)   !holds in each cloud level the number of
                                           !activated cloud droplets at the base
     ENDDO !end jl
  ENDDO !end jk

  ll1_2d(:,:) = (icl_minusbas(:,:) > 0)     .AND. & !all cloud levels minus base
                (zqlnuc_bas(:,:)   > 0._dp)         !newly cloud droplets at base

  zqlnuc(:,:)    = MERGE(zqlnuc_bas(:,:), zqlnuc(:,:), ll1_2d(:,:))
  ztmp(1:kproma,:)      = qnuc(1:kproma,:,jrow) + zdt*(zcdnc_bas(1:kproma,:) - zcdnc(1:kproma,:))
  qnuc(1:kproma,:,jrow) = MERGE(ztmp(1:kproma,:), qnuc(1:kproma,:,jrow), ll1_2d(:,:))
  zcdnc(:,:)     = MERGE(zcdnc_bas(:,:), zcdnc(:,:), ll1_2d(:,:))

!--- obtain a number of cloud droplets for the detrained cloud water from convective anvils
  
  DO jk=ktdia,klev
     iclbas(:,jk) = NINT(pcvcbot(:)) 
     ll1_2d(:,jk) = (jk == iclbas(:,jk))
  ENDDO

  ll2_2d(:,:) = (iclbas(:,:) > 0) 

  IF (ncvmicro == 0) THEN  !no CONV micro scheme

     zqlnuccvh(:,:) = MERGE(cdncact_cv(1:kproma,:,jrow)-zcdnc(:,:), 0._dp, ll1_2d(:,:))

     ll3_2d(:,:) = ll2_2d(:,:)           .AND. &
                   (pxtec(:,:) > 0._dp)  .AND. &
                   (ptm1(:,:)  > cthomi)

     ztmp1(:,:) = cdncact_cv(1:kproma,:,jrow)-zcdnc(:,:)

     DO jk=ktdia,klev
        DO jl=1,kproma
           jkk = iclbas(jl,jk)
           jkk = MAX(1, jkk) !SF prevents cases where iclbas=0
           ztmp2(jl,jk) = zqlnuccvh(jl,jkk) 
        ENDDO
     ENDDO

     ztmp(:,:) = MIN(ztmp1(:,:),ztmp2(:,:))
     ztmp(:,:) = MAX(0._dp,ztmp(:,:))

  ELSEIF (ncvmicro > 0) THEN  !CONV micro scheme ON

     ll3_2d(:,:) = ll2_2d(:,:)           .AND. &
                   (pxtecl(:,:) > 0._dp)

     ztmp(:,:) = pxtecnl(:,:)-zcdnc(:,:)
     ztmp(:,:) = MAX(0._dp, ztmp(:,:))

  ENDIF !end ncvmicro

  zqlnuccv(:,:)  = MERGE(ztmp(:,:), 0._dp, ll3_2d(:,:))
  qnuc(1:kproma,:,jrow) = qnuc(1:kproma,:,jrow) + zdt*zqlnuccv(1:kproma,:)
  zcdnc(:,:)     = zcdnc(:,:)     +     zqlnuccv(:,:)

!  Ice nucleation

  DO jk=klev,ktdia,-1
     DO jl=1,kproma

        ztc = ptm1(jl,jk)-tmelt
        ztc = MIN(0._dp,ztc)
!
!       test using Faisal's param
!
        zrieff = 0.015_dp*ztc
        zrieff = EXP(zrieff)
        zrieff = 23.2_dp*zrieff
        zrieff = MAX(zrieff,1.0_dp)
     
        zrih = 5113188.044_dp+2809._dp*zrieff**3
        zrih = SQRT(zrih)
        zrih = -2261.236_dp+zrih

        zrid_2d(jl,jk) = 1.e-6_dp*zrih**(1._dp/3._dp)
        zrid_2d(jl,jk) = MAX(1.e-6_dp,zrid_2d(jl,jk))

     END DO !jl
  END DO !jk

  znidetr(:,:)=0._dp  

  IF (ncvmicro > 0) THEN !SF conv micro. scheme ON

     ll_cv(:,:) = (pxteci(:,:)  > 0._dp) 

     ztmp(:,:) = pxtecni(:,:) - zicncq(:,:)
     ztmp(:,:) = MAX(ztmp(:,:), 0._dp) 


  ELSEIF (ncvmicro == 0) THEN !SF conv micro. scheme OFF

     ll_cv(:,:) = (pxtec(:,:) > 0._dp) .AND. &
                  (ptm1(:,:)  < cthomi)

     ztmp(:,:) = 0.75_dp*ztmst*zrho(:,:)*pxtec(:,:)/(api*zrhoice*zrid_2d(:,:)**3) - zicncq(:,:)
     ztmp(:,:) = MAX(ztmp(:,:), 0._dp)

  ENDIF

  znidetr(:,:) = MERGE(ztmp(:,:),0._dp,ll_cv(:,:))
  zicncq(:,:)  = zicncq(:,:) + znidetr(:,:)

  ! Provide soluble aerosol number number concentratios for the cirrus scheme:
  zascs(:,:) = 0._dp
    ! Calculate the sum by explicitly setting the involved modes
    ! (This is necessary for MADE3, to avoid including the mixed modes) 
     zascs(1:kproma,:) = mode(iaits)%aernum(1:kproma,:) + &
                         mode(iaccs)%aernum(1:kproma,:) + &
                         mode(icoas)%aernum(1:kproma,:)

  zascs(:,:) = MAX(zascs(:,:), 10.E6_dp*zrho(:,:))! change factor

!--- calculate the vertical velocity necessary for cirrus formation
!SF put this outside of the cirrus scheme choice, since this is also needed with nicnc=1

     ztmp(:,:)      = 70._dp*SQRT(ptkem1(:,:))
     ztmp(:,klev)   = 0.0_dp
     zvervx_2d(:,:) = -100._dp*g_rcp*pvervel(:,:)*zrho_rcp(:,:) + ztmp(:,:)

  !--------------- Variables for BN09 ------------------------------------------------------
  !Store variables
  w_sub(:,:)  = ztmp(:,:)                                   !turbulent contribution to vertical velocity [cm/s]
  w_ave(:,:)  = -100._dp*g_rcp*pvervel(:,:)*zrho_rcp(:,:)   !average (resolved) vertical velocity [cm/s]
  w_grid(:,:) = zvervx_2d(:,:)                              !vertical velocity in the grid cell (= w_ave + w_sub)

  !Prepare input variables to use BN09
  IF ( nicnc == 3 .OR. limm_BN09 ) THEN
  
     !Aerosol concentrations [m-3] and diameters [m]
     rdryki(1:kproma,:) = mode(iaiti)%dryrad(1:kproma,:)       !dry radius ki mode (iaiti=5)
     rdryai(1:kproma,:) = mode(iacci)%dryrad(1:kproma,:)       !dry radius ai mode (iacci=6)
     rdryci(1:kproma,:) = mode(icoai)%dryrad(1:kproma,:)       !dry radius ci mode (icoai=7)
     rdryks(1:kproma,:) = mode(iaits)%dryrad(1:kproma,:)       !dry radius ks mode (iaits=2)
     rdryas(1:kproma,:) = mode(iaccs)%dryrad(1:kproma,:)       !dry radius as mode (iaccs=3)
     rdrycs(1:kproma,:) = mode(icoas)%dryrad(1:kproma,:)       !dry radius cs mode (icoas=4)

     sigmaki = mode(iaiti)%sigma                               !standard deviation ki mode
     sigmaai = mode(iacci)%sigma                               !standard deviation ai mode
     sigmaci = mode(icoai)%sigma                               !standard deviation ci mode
     sigmaks = mode(iaits)%sigma                               !standard deviation ks mode
     sigmaas = mode(iaccs)%sigma                               !standard deviation as mode
     sigmacs = mode(icoas)%sigma                               !standard deviation cs mode

     !Turbulent contribution to vertical velocity (input of BN09)
     sigw(:,:) = ztmp(:,:)/100._dp            !from [cm/s] to [m/s]

     DO jk = ktdia,klev
        DO jl = 1,kproma
   
           !---- dust (CNT, PDA08, PDA13 spectra)
           DO I = 1,2                   !insoluble modes of DU cycle
              IF (I==1) THEN                              !insoluble DU (ai mode)
                 ndust(jl,jk,I) = MAX(nduinsolai(jl,jk,jrow), 0.0_dp)
                 ddust(jl,jk,I) = MAX(rdryai(jl,jk)*2, 1E-9_dp)
                 sigdust(I)     = sigmaai
              END IF
       
              IF (I==2) THEN                               !insoluble DU (ci mode)
                 ndust(jl,jk,I) = MAX(nduinsolci(jl,jk,jrow), 0.0_dp)
                 ddust(jl,jk,I) = MAX(rdryci(jl,jk)*2, 1E-8_dp)
                 sigdust(I)     = sigmaci
              END IF
           END DO
 
           !---- black carbon (CNT, PDA08, PDA13 spectra) 
           nsoot(jl,jk) = MAX(nbcinsol(jl,jk,jrow), 0.0_dp)         !insoluble BC (ki mode)
           dsoot(jl,jk) = MAX(rdryki(jl,jk)*2, 1E-9_dp)
           sigsoot      = sigmaki
       
           !---- insoluble organics (only for PDA08 spectrum)
           norg(jl,jk)  = MAX(nocinsol(jl,jk,jrow), 0.0_dp)        !insoluble OC (ki mode)
           dorg(jl,jk)  = MAX(rdryki(jl,jk)*2, 1E-9_dp)
           sigorg       = sigmaki     
       
           !---- soluble organics (only for PDA13 spectrum)
           DO J = 1,3                   !soluble modes of OC cycle  
             IF (J==1) THEN                              !soluble OC (ks mode)
                 nsolo(jl,jk,J) = MAX(nocsolks(jl,jk,jrow), 0.0_dp)
                 dsolo(jl,jk,J) = MAX(rdryks(jl,jk)*2, 1E-9_dp)
                 sigsolo(J)     = sigmaks
              END IF

              IF (J==2) THEN                             !soluble OC (as mode)
                 nsolo(jl,jk,J) = MAX(nocsolas(jl,jk,jrow), 0.0_dp)
                 dsolo(jl,jk,J) = MAX(rdryas(jl,jk)*2, 1E-9_dp)
                 sigsolo(J)     = sigmaas
              END IF

              IF (J==3) THEN                             !soluble OC (cs mode)
                 nsolo(jl,jk,J) = MAX(nocsolcs(jl,jk,jrow), 0.0_dp)
                 dsolo(jl,jk,J) = MAX(rdrycs(jl,jk)*2, 1E-8_dp)
                 sigsolo(J)     = sigmacs
              END IF
           END DO
       
           !---- biological aerosols (only for PDA13 spectrum)
           nbio(jl,jk) = 0.0_dp                            !insoluble BIO (ki mode)
           dbio(jl,jk) = MAX(rdryki(jl,jk)*2, 1E-9_dp)
           sigbio      = sigmaki     

           !---- droplets and sulfate
           ndrop(jl,jk) = MAX(pcdncact(jl,jk), 0.0_dp)      !total cloud droplet number conc
           dsulf(jl,jk) = MAX(rdryks(jl,jk)*2, 1E-9_dp)     !soluble SULF (ks mode) 

        END DO     !kproma
     END DO        !klev

     !INPUT variables saved in the new channel cloud_ice_BN09
     sigwBN(1:kproma,:)     = sigw(:,:) 
     ndropBN(1:kproma,:)    = ndrop(:,:)
     dsulfBN(1:kproma,:)    = dsulf(:,:)
     ndust_aiBN(1:kproma,:) = ndust(:,:,1)
     ddust_aiBN(1:kproma,:) = ddust(:,:,1)
     ndust_ciBN(1:kproma,:) = ndust(:,:,2) 
     ddust_ciBN(1:kproma,:) = ddust(:,:,2)
     norgBN(1:kproma,:)     = norg(:,:) 
     dorgBN(1:kproma,:)     = dorg(:,:) 
     nsootBN(1:kproma,:)    = nsoot(:,:) 
     dsootBN(1:kproma,:)    = dsoot(:,:) 
     nsolo_ksBN(1:kproma,:) = nsolo(:,:,1)
     dsolo_ksBN(1:kproma,:) = dsolo(:,:,1)
     nsolo_asBN(1:kproma,:) = nsolo(:,:,2)
     dsolo_asBN(1:kproma,:) = dsolo(:,:,2)
     nsolo_csBN(1:kproma,:) = nsolo(:,:,3) 
     dsolo_csBN(1:kproma,:) = dsolo(:,:,3)
  END IF 
  !--------------- END Variables for BN09------------------------------------------------------

  !----------------------------------------------------------------------------------------
  !--- Different CIRRUS PARAMETERIZATIONS
  !   
  !--- nicnc = 1
  !--- nicnc = 2  -> use Kaercher and Lohmann, 2002b
  !--- nicnc = 3  -> use Barahona and Nenes, 2009b
  !----------------------------------------------------------------------------------------

  IF ( nicnc == 1 ) THEN ! Use ICNC scheme after Lohmann, JAS, 2002
        
     ll_ice(:,:) = (sice(1:kproma,:,jrow) > 0._dp) .AND. &
                   (ptm1(:,:)      < cthomi)

     ztmp1(:,:) = 0.75_dp*zrho(:,:)*sice(1:kproma,:,jrow)*&
       zqsi_2d(:,:)/(api*zrhoice*zrid_2d(:,:)**3)-zicncq(:,:)
     ztmp2(:,:) = zascs(:,:)-zicncq(:,:)
     ztmp(:,:)  = MIN(ztmp1(:,:),ztmp2(:,:))
     ztmp(:,:)  = MAX(ztmp(:,:),0._dp)

     zninucl(:,:) = MERGE(ztmp(:,:), 0._dp, ll_ice(:,:))
     zicncq(:,:)   = zicncq(:,:) + zninucl(:,:)  

  ELSEIF ( nicnc > 1 ) THEN !--- Use Bernd's cirrus scheme 

     zsusatix_2d(:,:) = sice(1:kproma,:,jrow)
     znicex(:,:)      = 0._dp

     !-------- Cirrus Formation Parameterization by Kaercher and Lohmann--------------
     IF ( nicnc == 2 ) THEN
        !Kaercher and Lohmann, JGR, 2002b; Lohmann et al., JGR, 2004 (heterog)
        !This implies that supersaturation with respect to ice is allowed, thus the depositional
        !growth equation needs to be solved

        !The freezing rate is limited by the number of aerosols available in each mode. 
        !For homogeneous freezing, these are the soluble aerosols except for the nucleation mode

     IF (lhet) THEN

         DO jk = klev, ktdia, -1
             zapnx(:,1)   = 1.e-6_dp*(ndusol_strat(1:kproma,jk,jrow) - zicncq(:,jk))   ![1/cm3]
             zapnx(:,1)   = MAX(zapnx(:,1),1.e-6_dp)

             ! Soluble dust is in the mixed mode for MADE3 (--> iaccm)
             ! For GMXe iaccm=iaccs by definition (see cloud_init_coupling)
             zaprx(:,1)   = 100._dp*mode(iaccm)%wetrad(1:kproma,jk) ![cm]
             zaprx(:,1)   = MAX(zaprx(:,1),mode(iaccs)%crdiv)

             zapsigx(:,1) = 1.0_dp

             zap(:,jk) = zapnx(:,1)

!----------- the freezing parameterization returns the number of newly formed ice crystals
!----------- (znicex) and their size (zri)

             CALL xfrzmstr(lhet, nosize, ztmst,                      &
                           klev, kbdim, kproma, jk, nfrzmod,         &
                           zsusatix_2d(:,jk), zvervx_2d(:,jk), zapnx,                  &
                           zaprx, zapsigx, ptm1, tmelt, zeps, papm1, &
                           cthomi, zri, znicex)

         ENDDO !jk

     ELSE

         DO jk = klev, ktdia, -1
  
            IF (.NOT. lomonodisp) THEN !SF standard 3-modes freezing
      
                zapnx(1:kproma,1)   = 1.e-6_dp                     &
                             * ( mode(iaits)%aernum(1:kproma,jk) &
                               - zicncq(:,jk) )![1/cm3]
                zapnx(:,1)   = MAX(zapnx(:,1),1.e-6_dp)

                zapnx(1:kproma,2)   = 1.e-6_dp                     &
                             * ( mode(iaccs)%aernum(1:kproma,jk) &
                               - zicncq(:,jk) )![1/cm3]
                zapnx(:,2)   = MAX(zapnx(:,2),1.e-6_dp)

                zapnx(1:kproma,3)   = 1.e-6_dp                     &
                             * ( mode(icoas)%aernum(1:kproma,jk) &
                               - zicncq(:,jk) )![1/cm3]
                zapnx(:,3)   = MAX(zapnx(:,3),1.e-6_dp)

                zaprx(1:kproma,1)   = 100._dp*mode(iaits)%wetrad(1:kproma,jk) ![cm]
                zaprx(:,1)   = MAX(zaprx(:,1),mode(iaits)%crdiv)

                zaprx(1:kproma,2)   = 100._dp*mode(iaccs)%wetrad(1:kproma,jk) ![cm]
                zaprx(:,2)   = MAX(zaprx(:,2),mode(iaccs)%crdiv)

                zaprx(1:kproma,3)   = 100._dp*mode(icoas)%wetrad(1:kproma,jk) ![cm]
                zaprx(:,3)   = MAX(zaprx(:,3),mode(icoas)%crdiv)


                zapsigx(:,1) = mode(iaits)%sigma
                zapsigx(:,2) = mode(iaccs)%sigma
                zapsigx(:,3) = mode(icoas)%sigma

                zap(:,jk) = zapnx(:,1) + zapnx(:,2) + zapnx(:,3)

            ELSE !SF monodisperse aerosol distribution

                zapnx(:,1)   = 1.e-6_dp * ( zascs(:,jk) - zicncq(:,jk) ) ![1/cm3]
                zapnx(:,1)   = MAX(zapnx(:,1),1.e-6_dp)

                zaprx(1:kproma,1)   = 100._dp*mode(iaccs)%wetrad(1:kproma,jk) ![cm]
                zaprx(:,1)   = MAX(zaprx(:,1),mode(iaccs)%crdiv)

                zapsigx(:,1) = 1._dp

                zap(:,jk) = zapnx(:,1) 

            ENDIF !SF end lomonodisp

!----------- the freezing parameterization returns the number of newly formed ice crystals
!----------- (znicex) and their size (zri)

            CALL xfrzmstr(lhet, nosize, ztmst,                      &
                          klev, kbdim, kproma, jk, nfrzmod,         &
                          zsusatix_2d(:,jk), zvervx_2d(:,jk), zapnx,                  &
                          zaprx, zapsigx, ptm1, tmelt, zeps, papm1, &
                          cthomi, zri, znicex)

         ENDDO !jk

     ENDIF !lhet

     END IF  !nicnc=2 
     !-------- End Cirrus Formation Parameterization by Kaercher and Lohmann----------

     !-------- Cirrus Formation Parameterization by Barahona and Nenes----------------
     IF ( nicnc == 3 ) THEN
 
        !New ice crystals formed in the cirrus regime using the parameterization 
        !by Barahona and Nenes 2009b, when T <= 238 
        !=> homog and heterog competition.
         
        DO jk = ktdia,klev
           DO jl = 1,kproma

              !initialize output
              smaxice = 0.0_dp
              sc_ice  = 0.0_dp
              nlim    = 0.0_dp
              nhet    = 0.0_dp
              nice    = 0.0_dp
              dice    = 0.0_dp
              sigwparc= 0.0_dp
 
              zap(jl,jk) = ( (sum(ndust(jl,jk,:))+norg(jl,jk)+nsoot(jl,jk)+nbio(jl,jk)+sum(nsolo(jl,jk,:))+   &
                             ndrop(jl,jk)) - zicncq(jl,jk) ) * 1.E-6_dp                            ![zap]= [1/cm3] 

              !---------------------------------------
              !--- Call the BN09 scheme
              !---------------------------------------
              IF (ptm1(jl,jk) <= cthomi) THEN       !SB: IF T<=238, cirrus regime

                 CALL ice_activate (ptm1(jl,jk), paphm1(jl,jk), sigw(jl,jk), zicncq(jl,jk), & 
                                    ndrop(jl,jk),   dsulf(jl,jk),                           & 
                                    ndust(jl,jk,:), ddust(jl,jk,:), sigdust(:),             &
                                    norg(jl,jk),    dorg(jl,jk),    sigorg,                 &
                                    nsoot(jl,jk),   dsoot(jl,jk),   sigsoot,                &
                                    nbio(jl,jk),    dbio(jl,jk),    sigbio,                 &
                                    nsolo(jl,jk,:), dsolo(jl,jk,:), sigsolo(:),             &
                                    smaxice, nlim, sc_ice, nhet, nice, dice, sigwparc)  

              END IF  !condition for cirrus 

             !Variables saved in the new channel cloud_ice_BN09
             !Output:
             smaxice_cirrusBN(jl,jk)  = smaxice
             sc_ice_cirrusBN(jl,jk)   = sc_ice
             nlim_cirrusBN(jl,jk)     = nlim
             nhet_cirrusBN(jl,jk)     = nhet
             nice_cirrusBN(jl,jk)     = nice
             dice_cirrusBN(jl,jk)     = dice
             sigwpre_cirrusBN(jl,jk)  = sigwparc


             !The output variables of BN09 are assigned to the "old" variables znicex and zri
             znicex(jl,jk) = nice         !new ice crystals [1/m^3]
             zri(jl,jk)    = dice/2._dp   !ice crystal radius [m]

           END DO     !kproma
        END DO        !klev

     END IF  !nicnc=3
     !-------- End Cirrus Formation Parameterization by BN09b ------------------------

!--- Update ICNC

     zri(:,:)=MAX(zri(:,:), 1.e-6_dp)
     
     ll_ice(:,:) = (zsusatix_2d(:,:) > 0._dp) .AND. &
                   (ptm1(:,:)        < cthomi)

     ! the subtraction (-zicncq) is already done in the previous definitions
     !of zapnx (nicnc=2) or zap (nicnc=3) 
     !!ztmp1(:,:) = 1.e6_dp*zap(:,:) - zicncq(:,:)
     ztmp1(:,:) = 1.e6_dp*zap(:,:)                     !SB: [ztmp1]=[1/m3], [zap]=[1/cm3] 

     ztmp2(:,:) = znicex(:,:)
     ztmp(:,:)  = MIN(ztmp1(:,:),ztmp2(:,:)) 
     ztmp(:,:)  = MAX(ztmp(:,:), 0._dp)

     ! store variables in cirrus regime 
     !independentely on the scheme (nicnc) used.
     newIC_cirri(:,:)  = ztmp(:,:)     ![1/m3]     !newly formed IC in cirrus clouds (T<=238 K)
     newICR_cirri(:,:) = zri(:,:)      ![m]        !their radius in cirrus clouds (T<=238 K)

     zninucl(:,:) = MERGE(ztmp(:,:), 0._dp,ll_ice(:,:))
     zicncq(:,:)  = zicncq(:,:) + zninucl(:,:)

!---calculate the deposition rate taking ventilation (zfre) into account

     ll_ice(:,:) = (pxim1(:,:) > 0._dp)

     ztmp(:,:)   = MAX(zicncq(:,:), zicemin)
     zicncq(:,:) = MERGE(ztmp(:,:), zicncq(:,:), ll_ice(:,:)) 

     ztmp(:,:) = zrho(:,:)*pxim1(:,:)/zicncq(:,:)
     ztmp(:,:) = MAX(ztmp(:,:), zmi)
     zmmean_2d(:,:) = MERGE(ztmp(:,:), zmmean_2d(:,:), ll_ice(:,:))

     ll1_2d(:,:) = (zmmean_2d(:,:) < 2.166E-9_dp )
     ll2_2d(:,:) = (zmmean_2d(:,:) >= 2.166E-9_dp  ) .AND. (zmmean_2d(:,:) < 4.264E-8_dp )

     zalfased_2d(:,:) = MERGE(63292.4_dp, 8.78_dp, ll1_2d(:,:))
     zalfased_2d(:,:) = MERGE(329.75_dp, zalfased_2d(:,:), ll2_2d(:,:))

     zbetased_2d(:,:) = MERGE(0.5727_dp, 0.0954_dp, ll1_2d(:,:))
     zbetased_2d(:,:) = MERGE(0.3091_dp, zbetased_2d(:,:), ll2_2d(:,:))

     zxifallmc_2d(:,:) = zfall*zalfased_2d(:,:) &
                        *(zmmean_2d(:,:)**zbetased_2d(:,:))*zaaa_2d(:,:) !Fallgeschwindigkeit Masse

     ztmp(:,:)       = pqm1(:,:) - zqsi_2d(:,:)

     DO jk = ktdia, klev
        DO jl = 1,kproma
           zdv     = 2.21_dp/papm1(jl,jk)
           zgtp    = 1._dp/(zrho(jl,jk)*zastbsti(jl,jk))
           zvth    = SQRT( 8._dp*zkb*ptm1(jl,jk) / (api*zxmw) )
           zb2     = 0.25_dp * zalpha * zvth   / zdv
           zfuchs  = 1._dp/(1._dp+zb2*zri(jl,jk))
           zre     = 2._dp*zrho(jl,jk)*zri(jl,jk)*zxifallmc_2d(jl,jk)/zviscos_2d(jl,jk)
           zfre    = 1._dp + 0.229_dp*SQRT(zre)
           zfre    = MERGE(zfre, 1._dp, ll_ice(jl,jk))

           zqinucl(jl,jk)  = 4._dp*api*zri(jl,jk)*zsusatix_2d(jl,jk)*zicncq(jl,jk) &
                           *zfre*zgtp*zfuchs*zalpha*ztmst
           zqinucl(jl,jk) = MIN(zqinucl(jl,jk),ztmp(jl,jk))
           zqinucl(jl,jk) = MAX(zqinucl(jl,jk),-pxim1(jl,jk))

        ENDDO !jl
     ENDDO !jk
 
  ENDIF   !which nucleation scheme (nicnc)

  !corinna: set zcdnc and zicnc to minium now if nucleation is not strong enough

  ll1(:,:) = ( paclc(:,:) >= zepsec ) .AND. &
             ( ptm1(:,:)  >  cthomi )

  ztmp(:,:)  = MAX(zcdnc(:,:),cdncmin)
  zcdnc(:,:) = MERGE(ztmp(:,:),zcdnc(:,:),ll1(:,:))

  ll1(:,:) = ( paclc(:,:) >= zepsec ) .AND. &
             ( ptm1(:,:)  <  tmelt  )

  ztmp(:,:)   = MAX(zicncq(:,:), zicemin)
  zicncq(:,:) = MERGE(ztmp(:,:), zicncq(:,:), ll1(:,:))

!--- End included for CDNC/IC scheme -----------------------------------

  DO 831 jk=ktdia,klev
!
!     ------------------------------------------------------------------
!
!       2.    Set to zero local tendencies (increments)
    ! mz_ht_20120209+
! pointer on the rain and snowflux through the bottom of each layer
    zrfl => pfrain(:,jk)
    zsfl => pfsnow(:,jk)
    if (jk > 1 ) then
      zrfl(1:kproma) = pfrain(1:kproma, jk-1)
      zsfl(1:kproma) = pfsnow(1:kproma, jk-1)
    ENDIF
! pointer on the evaporation of rain
    zevp => prevap(:,jk)
! pointer on the sublimation of snow
    zsub => pssubl(:,jk)
! pointer on the rain production
    zrpr => prate_r(:,jk)
    ! mz_ht_20120209-

!
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
!--- Included for prognostic CDNC/IC scheme ----------------------------
!    additional arrays for transformation rate changes of CDNC and ICNC
     zrprn(:)      = 0.0_dp
     zfrln(:)      = 0.0_dp
     zsacln(:)     = 0.0_dp
     zsprn(:,jk)   = 0.0_dp
  
     ztmp(:,jk)    = MAX(pqm1(:,jk),0.0_dp)
     ztmp(:,jk)    = 1._dp/(cpd+zcons1*ztmp(:,jk))
     zlvdcp(:)     = alv*ztmp(:,jk)
     zlsdcp(:)     = als*ztmp(:,jk)
!
!     ------------------------------------------------------------------
!
!       3.   Modification of incoming precipitation fluxes by
!            melting, sublimation and evaporation
!
     IF (jk > 1) THEN
!
!       3.1   Melting of snow and ice
!
        ztmp1_1d(:) = ptm1(:,jk) - tmelt
       
        ll1_1d(:) = (ztmp1_1d(:) > 0._dp)
 
        ztmp1_1d(:) = MAX(0._dp, ztmp1_1d(:))   !SF ztdif
        ztmp1_1d(:) = zcons2*ztmp1_1d(:)*zdp_2d(:,jk) / ( zlsdcp(:)-zlvdcp(:) ) !SF zcons*ztdif
        
        ztmp2_1d(:) = zxsec*zsfl(:)
        ztmp2_1d(:) = MIN(ztmp2_1d(:), ztmp1_1d(:)) !SFzximelt
        pimelt(:,jk) = ztmp2_1d(:)

        zrfl(:)  = zrfl(:) + ztmp2_1d(:)
        zsfl(:)  = zsfl(:) - ztmp2_1d(:)


        zsmlt(:) = ztmst*g*ztmp2_1d(:) / zdp_2d(:,jk)

        ztmp2_1d(:) = zxsec*zxiflux(:)
        ztmp2_1d(:) = MIN(ztmp2_1d(:), ztmp1_1d(:))

        ll2_1d(:) = (zxiflux(:) > zepsec)

        ztmp3_1d(:) = zxifluxn(:)*ztmp2_1d(:)/MAX(zxiflux(:),zepsec)
        ztmp3_1d(:) = MERGE(ztmp3_1d(:), 0._dp, ll2_1d(:)) !SF zxinmelt

        zxiflux(:)  = zxiflux(:)  - ztmp2_1d(:)
        zxifluxn(:) = zxifluxn(:) - ztmp3_1d(:)

        ll2_1d(:) = (zxiflux(:) < zepsec)

        zxifluxn(:) = MERGE(0._dp, zxifluxn(:), ll2_1d(:))
 
        zximlt(:)   = ztmst*g*ztmp2_1d(:) / zdp_2d(:,jk)

        ztmp1_1d(:) = pxim1(:,jk) + ztmst*pxite(:,jk)
        ztmp1_1d(:) = MAX(0._dp, ztmp1_1d(:))
        zimlt(:)    = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))

!--- Included for prognostic CDNC/IC scheme (Philip Stier, 31/03/2004) -
!    If T > tmelt melt all ice crystals and transfer to cloud droplets 

        ztmp1_1d(:) = MERGE(zicncq(:,jk), 0._dp, ll1_1d(:))
        zicnc(:,jk) = MERGE(zicemin, zicnc(:,jk), ll1_1d(:)) 

        zcdnc(:,jk) = zcdnc(:,jk) + ztmp1_1d(:)
        qmel(1:kproma,jk,jrow) = qmel(1:kproma,jk,jrow) + zdt*ztmp1_1d(:)

!--- End included for CDNC/IC scheme -----------------------------------
!
!       3.2   Sublimation of snow (zsub) and ice (zxisub) following (Lin et al., 1983)
!
        ll1_1d(:) = (zclcpre(:) > 0._dp)

        ztmp1_1d(:) = 1._dp/(2.43e-2_dp*rv)*zlsdcp(:)**2/ptm1(:,jk)**2
        ztmp1_1d(:) = ztmp1_1d(:) &
                    + 1._dp/0.211e-4_dp*zrho_rcp(:,jk)/zqsi_2d(:,jk)
        ztmp1_1d(:) = 3.e6_dp*2._dp*api*zicesub(:,jk)*zrho_rcp(:,jk)/ztmp1_1d(:) !zcoeff

        ztmp2_1d(:) = MERGE(zclcpre(:), 1._dp, ll1_1d(:))

!SF     Snow:  
        ll2_1d(:) = (zsfl(:) > cqtmin) .AND. &
                    ll1_1d(:)

        ztmp3_1d(:) = zcons3*(zsfl(:) / ztmp2_1d(:))**(0.25_dp/1.16_dp)  !zclambs
        ztmp3_1d(:) = 0.78_dp  *ztmp3_1d(:)**2                                & 
                    + 232.19_dp*zqrho_2d(:,jk)**0.25_dp * ztmp3_1d(:)**2.625_dp     !zcfac4c
        ztmp3_1d(:) = ztmp3_1d(:) * ztmp1_1d(:) * zdpg_2d(:,jk)                     

        ztmp4_1d(:) = -zxsec * zsfl(:) / ztmp2_1d(:)
        ztmp4_1d(:) = MAX(ztmp4_1d(:), ztmp3_1d(:))                              !zzeps
        ztmp4_1d(:) = -ztmst*ztmp4_1d(:) / zdpg_2d(:,jk) * ztmp2_1d(:)

        ztmp5_1d(:) = zxsec*(zqsi_2d(:,jk)-pqm1(:,jk))
        ztmp5_1d(:) = MAX(ztmp5_1d(:),0._dp)

        ztmp4_1d(:) = MIN(ztmp4_1d(:),ztmp5_1d(:))
        ztmp4_1d(:) = MAX(ztmp4_1d(:),0._dp)
        zsub(:)     = MERGE(ztmp4_1d(:), zsub(:), ll2_1d(:))

!SF     Ice:
        ll2_1d(:) = (zxiflux(:) > cqtmin) .AND. &
                    ll1_1d(:)

        ztmp3_1d(:) = zcons3*(zxiflux(:) / ztmp2_1d(:))**(0.25_dp/1.16_dp) !zclambs
        ztmp3_1d(:) = 0.78_dp  *ztmp3_1d(:)**2                                &
                    + 232.19_dp*zqrho_2d(:,jk)**0.25_dp * ztmp3_1d(:)**2.625_dp     !zcfac4c
        ztmp3_1d(:) = ztmp3_1d(:) * ztmp1_1d(:) * zdpg_2d(:,jk)

        ztmp4_1d(:) = -zxsec * zxiflux(:) / ztmp2_1d(:)
        ztmp4_1d(:) = MAX(ztmp4_1d(:), ztmp3_1d(:))                              !zzeps
        ztmp4_1d(:) = -ztmst*ztmp4_1d(:) / zdpg_2d(:,jk) * ztmp2_1d(:)

        ztmp5_1d(:) = zxsec*(zqsi_2d(:,jk)-pqm1(:,jk))
        ztmp5_1d(:) = MAX(ztmp5_1d(:),0._dp)

        ztmp4_1d(:) = MIN(ztmp4_1d(:),ztmp5_1d(:))
        ztmp4_1d(:) = MAX(ztmp4_1d(:),0._dp)
        zxisub(:)   = MERGE(ztmp4_1d(:), zxisub(:), ll2_1d(:))

        ztmp5_1d(:) = zxisub(:) * zxifluxn(:) / MAX(zxiflux(:),cqtmin)  !SF zsubin
        ztmp5_1d(:) = zcons2 * ztmp5_1d(:) * zdp_2d(:,jk)
        ztmp5_1d(:) = MERGE(ztmp5_1d(:), 0._dp, ll2_1d(:))
        zxifluxn(:) = zxifluxn(:) - ztmp5_1d(:)

        zxiflux(:) = zxiflux(:) - zcons2*zxisub(:)*zdp_2d(:,jk)

        ll3_1d(:) = ll2_1d(:) .AND. (zxiflux(:) < zepsec)

        zxifluxn(:) = MERGE(0._dp, zxifluxn(:), ll3_1d(:))

!
!       3.3   Evaporation of rain (zevp following Rotstayn, 1997)
!
        ll2_1d(:) = (zrfl(:) > cqtmin) .AND. &
                    ll1_1d(:)

        ztmp3_1d(:) = 870._dp * zsusatw_evap(:,jk) * zdpg_2d(:,jk)      &
                              * (zrfl(:)/ztmp2_1d(:))**0.61_dp     &
                              / (SQRT(zrho(:,jk))*zastbstw(:,jk))                 

        ztmp4_1d(:) = -zxsec*zrfl(:)/ztmp2_1d(:)
        ztmp4_1d(:) = MAX(ztmp4_1d(:), ztmp3_1d(:))
        ztmp4_1d(:) = -ztmst*ztmp4_1d(:)*ztmp2_1d(:)/zdpg_2d(:,jk)

        ztmp5_1d(:) = zxsec*(zqsw_2d(:,jk)-pqm1(:,jk))
        ztmp5_1d(:) = MAX(ztmp5_1d(:), 0._dp)

        ztmp4_1d(:) = MIN(ztmp4_1d(:),ztmp5_1d(:))
        ztmp4_1d(:) = MAX(ztmp4_1d(:),0._dp)
        zevp(:)     = MERGE(ztmp4_1d(:), zevp(:), ll2_1d(:))

        IF (lookupoverflow) THEN
          status_string = 'lookuperror: cdnc - cloud (1)'
          RETURN
        ENDIF
      ENDIF !SF end jk > 1
!
!     --- End included for CDNC/IC scheme -----------------------------------
!
!     ------------------------------------------------------------------
!       4.    Sedimentation of cloud ice from grid-mean values (zicesed).
!             Updating the tendency 'pxite' to include sedimentation.
!             At jk=klev, the sedimentation sink is balanced by
!             precipitation at the surface (through 'zzdrs', see 7.3).
!             Finally: In-cloud cloud water/ice.
!
        zxip1_1d(:) = pxim1(:,jk) + ztmst*pxite(:,jk)-zimlt(:)
        zxip1_1d(:) = MAX(zxip1_1d(:), EPSILON(1._dp))

        zicncp1_1d(:) = zicnc(:,jk)*paclc(:,jk)
        zicncp1_1d(:) = MAX(zicncp1_1d(:),zicemin)

        zmmean_2d(:,jk) = zrho(:,jk)*zxip1_1d(:)/zicncp1_1d(:)
        zmmean_2d(:,jk) = MAX(zmmean_2d(:,jk), zmi)

        ll1_1d(:) = (zmmean_2d(:,jk) < 2.166E-9_dp )
        ll2_1d(:) = (.NOT. ll1_1d(:)) .AND. (zmmean_2d(:,jk) < 4.264E-8_dp )

        zalfased_2d(:,jk) = MERGE(63292.4_dp, 8.78_dp, ll1_1d(:))
        zalfased_2d(:,jk) = MERGE(329.75_dp, zalfased_2d(:,jk), ll2_1d(:))

        zbetased_2d(:,jk) = MERGE(0.5727_dp, 0.0954_dp, ll1_1d(:))
        zbetased_2d(:,jk) = MERGE(0.3091_dp, zbetased_2d(:,jk), ll2_1d(:))

        zxifallmc_2d(:,jk) = zfall*zalfased_2d(:,jk) & 
                           *(zmmean_2d(:,jk)**zbetased_2d(:,jk))*zaaa_2d(:,jk) !Fallgeschwindigkeit Masse

        !limit fall velocity to 1 cm/s - 2 m/s:
        zxifallmc_2d(:,jk) = MAX(0.001_dp,zxifallmc_2d(:,jk))
        zxifallmc_2d(:,jk) = MIN(2._dp,zxifallmc_2d(:,jk))

        zxifallnc_2d(:,jk) = zxifallmc_2d(:,jk)                                !Fallgeschwindigkeit Anzahl

        ztmp1_1d(:) = ztmst*g*zxifallmc_2d(:,jk)*zrho(:,jk)/zdp_2d(:,jk) !SF zal1

        ll1_1d(:) = (zxifallmc_2d(:,jk) > zeps)

        ztmp2_1d(:) = zxiflux(:)*zrho_rcp(:,jk)/MAX(zxifallmc_2d(:,jk), zeps)  
        ztmp2_1d(:) = MERGE(ztmp2_1d(:), 0._dp, ll1_1d(:)) !SF zal2
        
        ztmp3_1d(:) = g*ztmst*zxifallnc_2d(:,jk)*zrho(:,jk)/zdp_2d(:,jk) !SF zal3

        ll1_1d(:)   = (zxifallnc_2d(:,jk) > zeps)
        ztmp4_1d(:) = zxifluxn(:)/MAX(zxifallnc_2d(:,jk), zeps)
        ztmp4_1d(:) = MERGE(ztmp4_1d(:), 0._dp, ll1_1d(:)) !SF zal4

!        -----------change steffi-------------------
        zxised_2d(:,jk)  =  zxip1_1d(:)  *EXP(-ztmp1_1d(:))+ztmp2_1d(:)*(1._dp-EXP(-ztmp1_1d(:)))
        zicesed_2d(:,jk) =  zicncp1_1d(:)*EXP(-ztmp3_1d(:))+ztmp4_1d(:)*(1._dp-EXP(-ztmp3_1d(:)))
!        --------end change steffi--------------------

        ll1_1d(:) = (paclc(:,jk) > 0.01_dp)

        ztmp1_1d(:) = zicesed_2d(:,jk) / MAX(paclc(:,jk), 0.01_dp)
        ztmp1_1d(:) = MERGE(ztmp1_1d(:), zicesed_2d(:,jk), ll1_1d(:))

        zicnc(:,jk) = ztmp1_1d(:) + znidetr(:,jk) + zninucl(:,jk)
        zicnc(:,jk) = MIN(zicnc(:,jk), zicemax)
        zicnc(:,jk) = MAX(zicnc(:,jk), zicemin)

        zxiflux(:)  = zxiflux(:)  + zcons2*(zxip1_1d(:)   - zxised_2d(:,jk) )*zdp_2d(:,jk)
        zxifluxn(:) = zxifluxn(:) + zcons2*(zicncp1_1d(:) - zicesed_2d(:,jk))*zdp_2d(:,jk)*zrho_rcp(:,jk)

        pxite(:,jk) = ztmst_rcp*(zxised_2d(:,jk)-pxim1(:,jk))

!---Included for in-cloud scavenging (Philip Stier, 16/04/04):----------
        zmrateps(:,jk) = zmrateps(:,jk) + zxip1_1d(:) - zxised_2d(:,jk)
!---End Included for scavenging-----------------------------------------
        pisedi(:,jk) = zxised_2d(:,jk)

!
!             In-cloud water/ice calculated from respective grid-means,
!             partial cloud cover, advective/diffusive tendencies,
!             detrained cloud water/ice and ice sedimentation.
!             In-cloud values are required for cloud microphysics.
!
!       THIS PART NEEDS TO BE COMMENTED BY ERICH/ADRIAN
!
        zclcaux(:) = paclc(:,jk)
        locc_1d(:) = (zclcaux(:) > zeps)
        
        ztmp1_1d(:) = 1._dp/MAX(pqm1(:,jk),zeps) + zlsdcp(:)*alv/(rv*ptm1(:,jk)**2) !SF protected against divide by zero 
        ztmp2_1d(:) = g*(zlvdcp(:)*rd/rv/ptm1(:,jk)-1._dp)/(rd*ptm1(:,jk))
        zkair_1d(:) = 4.1867e-3_dp*(5.69_dp + 0.017_dp*(ptm1(:,jk)-tmelt)) ! eq. 13-18a P&K !SF replaced ztp1tmp 
                                                                           ! by ptm1
        zdv_1d(:)   = 2.21_dp/papm1(:,jk)
        ztmp3_1d(:) = 1._dp/(crhoi*als**2/(zkair_1d*rv*ptm1(:,jk)**2) & 
                    + crhoi*rv*ptm1(:,jk)/(zesi_2d(:,jk)*zdv_1d(:)))

        zrice_1d(:)    = (0.75_dp*zxised_2d(:,jk)/(api*zicnc(:,jk)*crhoi))**(1._dp/3._dp)
        zeta_1d(:)     = ztmp1_1d(:)/ztmp2_1d(:)*ztmp3_1d(:)*4._dp*api*crhoi*zcap/zrho(:,jk)
        zvervmax_1d(:) = (zesw_2d(:,jk)-zesi_2d(:,jk))/zesi_2d(:,jk)*zicnc(:,jk) &
                                                      *zrice_1d(:)*zeta_1d(:)

        lo2_1d(:)  = (ptm1(:,jk) < cthomi) .OR.                   &
                      (ptm1(:,jk) < tmelt .AND. 0.01_dp*zvervx_2d(:,jk) < zvervmax_1d(:))


        ll1_1d(:)  =  (ptm1(:,jk) < tmelt .AND. ptm1(:,jk) > cthomi .AND. & 
             0.01_dp*zvervx_2d(:,jk) < zvervmax_1d(:) .AND. zclcaux(:) > 0._dp)

        IF(ncvmicro>0) THEN
            zxite(:) = pxteci(:,jk)
            zxlte(:) = pxtecl(:,jk)
        ELSE
            zxite(:) = MERGE(pxtec(:,jk), 0._dp      , lo2_1d(:))
            zxlte(:) = MERGE(0._dp      , pxtec(:,jk), lo2_1d(:))
        ENDIF

        ztmp1_1d(:) = pxim1(:,jk)/MAX(zclcaux(:),zeps) !SF to protect from a division by 0
        ztmp2_1d(:) = pxlm1(:,jk)/MAX(zclcaux(:),zeps) 
    
        zxib(:) = MERGE(ztmp1_1d(:), 0._dp, locc_1d(:))
        zxlb(:) = MERGE(ztmp2_1d(:), 0._dp, locc_1d(:))

        zxim1evp_1d(:) = MERGE(0._dp, pxim1(:,jk), locc_1d(:))
        zxlm1evp_1d(:) = MERGE(0._dp, pxlm1(:,jk), locc_1d(:)) 

        zxite2(:) = MERGE(zxite(:), 0._dp   , lo2_1d(:)) !SF temporary var needed if ncvmicro>0
        zxlte2(:) = MERGE(0._dp   , zxlte(:), lo2_1d(:)) !SF temporary var needed if ncvmicro>0

        zxidt_1d(:) = ztmst*(pxite(:,jk)+zxite2(:))
        zxldt_1d(:) = ztmst*(pxlte(:,jk)+zxlte2(:)) + zximlt(:) + zimlt(:) 

        !SF ice cloud:

        ll1_1d(:) = (zxidt_1d(:) > 0._dp)
        zxidtstar_1d(:) = MERGE(zxidt_1d(:), 0._dp, ll1_1d(:))
       
        ztmp1_1d(:) = zxidt_1d(:)/MAX(zclcaux(:), zeps)
        ztmp1_1d(:) = MAX(ztmp1_1d(:), -zxib(:))
        ztmp1_1d(:) = MERGE(zxidt_1d(:), ztmp1_1d(:), ll1_1d(:))
        ztmp1_1d(:) = MERGE(ztmp1_1d(:), 0._dp, locc_1d(:))
        zxib(:) = zxib(:) + ztmp1_1d(:)

        ztmp1_1d(:) = -( ztmst_rcp*pxim1(:,jk) + zxite2(:) )
        ztmp1_1d(:) = MAX(pxite(:,jk), ztmp1_1d(:))
        pxite(:,jk) = MERGE(pxite(:,jk), ztmp1_1d(:), ll1_1d(:))

        ll2_1d(:) = (.NOT. locc_1d(:)) .AND. (.NOT. ll1_1d(:))
        ztmp1_1d(:) = ztmst * ( pxite(:,jk) + zxite2(:) )
        ztmp1_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll2_1d(:))
        zxim1evp_1d(:) = zxim1evp_1d(:) + ztmp1_1d(:)

        !SF water cloud:

        ll1_1d(:) = (zxldt_1d(:) > 0._dp)
        zxldtstar_1d(:) = MERGE(zxldt_1d(:), 0._dp, ll1_1d(:))

        ztmp1_1d(:) = zxldt_1d(:)/MAX(zclcaux(:), zeps)
        ztmp1_1d(:) = MAX(ztmp1_1d(:), -zxlb(:))
        ztmp1_1d(:) = MERGE(zxldt_1d(:), ztmp1_1d(:), ll1_1d(:))
        ztmp1_1d(:) = MERGE(ztmp1_1d(:), 0._dp, locc_1d(:))
        zxlb(:)     = zxlb(:) + ztmp1_1d(:)

        ztmp1_1d(:) = -(pxlm1(:,jk)/ztmst+zxlte2(:))
        ztmp1_1d(:) = MAX(pxlte(:,jk), ztmp1_1d(:))
        pxlte(:,jk) = MERGE(pxlte(:,jk), ztmp1_1d(:), ll1_1d(:))

        ll2_1d(:)      = (.NOT. locc_1d(:)) .AND. (.NOT. ll1_1d(:))
        ztmp1_1d(:)    = ztmst * ( pxlte(:,jk) + zxlte2(:) )
        ztmp1_1d(:)    = MERGE(ztmp1_1d(:), 0._dp, ll2_1d(:))
        zxlm1evp_1d(:) = zxlm1evp_1d(:) + ztmp1_1d(:)
!
!     ------------------------------------------------------------------
!       5.    Condensation/deposition and evaporation/sublimation
!
!             zlc       =  L_{v/s} / c_p
!
        zlc_1d(:)   = MERGE(zlsdcp(:),zlvdcp(:),lo2_1d(:))
        zqsm1_1d(:) = MERGE(zqsi_2d(:,jk)  , zqsw_2d(:,jk)  , lo2_1d(:))
        zqst1_1d(:) = MERGE(zqsip1_2d(:,jk), zqswp1_2d(:,jk), lo2_1d(:))

        zdqsdt_1d(:) = 1000._dp*(zqst1_1d(:)-zqsm1_1d(:))

        zxievap(:) = (1.0_dp-zclcaux(:)) * zxidtstar_1d(:) + zxim1evp_1d(:)
        zxlevap(:) = (1.0_dp-zclcaux(:)) * zxldtstar_1d(:) + zxlm1evp_1d(:)

        zqvdt_1d(:) = ztmst*pqte(:,jk) + zevp(:)    + zsub(:)   &
                    + zxievap(:)       + zxlevap(:) + zxisub(:)

        zdtdt_1d(:) = ztmst*ptte(:,jk)                                       &
                    -  zlvdcp(:)           *(zevp(:)  + zxlevap(:))          &
                    - (zlsdcp(:)-zlvdcp(:))*(zsmlt(:) + zximlt(:)+zimlt(:))

        zdtdt_1d(:) = zdtdt_1d(:) - zlsdcp(:)*(zsub(:)+zxievap(:)+zxisub(:))

        zqp1_1d(:)  = pqm1(:,jk)+zqvdt_1d(:)
        zqp1_1d(:)  = MAX(zqp1_1d(:),0.0_dp)
!SF to put in lcover=false?
        ztp1_1d(:)   = ptm1(:,jk)+zdtdt_1d(:)
        zdqsat_1d(:) = zdtdt_1d(:)                                             &
                     + zclcaux(:)*( ztmst*zlc_1d(:)*pqte(:,jk)                 &
                                  + zlvdcp(:)*(zevp(:)+zxlevap(:))             &
                                  + zlsdcp(:)*(zsub(:)+zxievap(:)+zxisub(:)) )   !zdtdtstar 

        zdqsat_1d(:) = zdqsat_1d(:)*zdqsdt_1d(:)/(1._dp+zclcaux(:)*zlc_1d(:)*zdqsdt_1d(:))
!SF end to put in lcover=false

        zxib(:)     = MAX(zxib(:),0.0_dp)
        zxlb(:)     = MAX(zxlb(:),0.0_dp)
        zxilb_1d(:) = zxib(:)+zxlb(:)
!
!       Diagnostics: relative humidity
!
        prelhum(:,jk) = pqm1(:,jk)/zqsm1_1d(:)
        prelhum(:,jk) = MAX(MIN(prelhum(:,jk),1._dp),0._dp)

        IF (jk >= ncctop) THEN
!
!       define variables needed for cover scheme
!
!       zbetaqt = total water
!       zbetass = saturation mixing ratio adjusted to match qv
!       zwide   = current diagnosed distribution width
!
           zbetacl(:) = MAX(0.0_dp, pxlm1(:,jk)) + MAX(0.0_dp, pxim1(:,jk))
           zbetaqt(:) = MAX(cqtmin, pqm1(:,jk))  + zbetacl(:)
           zvartg(:)  = MAX(cqtmin,cvarmin*pqm1(:,jk))
           zwide(:)   = MAX(zvartg(:), pbetab(:,jk)-pbetaa(:,jk))
!
!       5.1 Turbulence: Skewness - equation solved implicitly
!           This solver only works if phmixtau has non-zero timescale
!
           zqtau_1d(:) = phmixtau(:,jk)+pvmixtau(:,jk) 

           zbqp1_1d(:) = -zdtime * zqtau_1d(:)
           zbqp1_1d(:) = EXP(zbqp1_1d(:))
           zbqp1_1d(:) = cbeta_pq - (cbeta_pq-pxskew(:,jk))*zbqp1_1d(:)
           zbqp1_1d(:) = MIN(zbqp1_1d(:), cbeta_pq_max)
           zbqp1_1d(:) = MAX(zbqp1_1d(:), cbeta_pq)

           zturbskew(:) = zdtime_rcp*(zbqp1_1d(:)-pxskew(:,jk))
!
!       5.2 Turbulence: variance - equation solved implicitly
!
           zbbap1_1d(:) = (cbeta_pq+pxskew(:,jk))**2    &
                        * (cbeta_pq+pxskew(:,jk)+1._dp) &
                        / (cbeta_pq*pxskew(:,jk)      )       !SF zeta

           zbbap1_1d(:) = zbbap1_1d(:)*pvdiffp(:,jk)/zwide(:)  !SF zprod

           zbbap1_1d(:) = zbbap1_1d(:)/zqtau_1d(:)                          &
                        + zvartg(:)                                         &
                        - (zbbap1_1d(:)/zqtau_1d(:) + zvartg(:) - zwide(:)) &
                          * EXP(-zdtime*zqtau_1d(:))

           zbbap1_1d(:) = MAX(zbbap1_1d(:),zvartg(:))
          
           ztmp1_1d(:) = 1._dp / cbeta_pq * zbetaqt(:) * (cbeta_pq + zbqp1_1d(:))


           zbbap1_1d(:) = MIN(zbbap1_1d(:), ztmp1_1d(:))

           zturbvar(:) = zdtime_rcp*(zbbap1_1d(:)-zwide(:))

           zbap1_1d(:) = zbetaqt(:)-cbeta_pq*zbbap1_1d(:)/(cbeta_pq+zbqp1_1d(:))
!
           IF (lcover) THEN ! mz_jd_20161011
!          translated into apparent xl,xi,q and heat sources
!          first order effect only, effect of evaporation of
!          cloud on qsat taken into account in thermodynamic budget
!          but does not change the mixing term here since that
!          would require iteration and is therefore neglected
!
!          calculate values after one timestep
!
           iqidx_1d(:) = (zbqp1_1d(:)-cbeta_pq)/rbetak+1._dp
           iqidx_1d(:) = LOG(iqidx_1d(:))
           iqidx_1d(:) = (nbetaq/cbetaqs)*iqidx_1d(:) + 0.5_dp
           iqidx_1d(:) = INT(iqidx_1d(:))

           ztmp1_1d(:) =   cbeta_pq*(pbetass(:,jk) - zbap1_1d(:))                   &
                       / ( (zbetaqt(:)    - zbap1_1d(:))*(cbeta_pq + zbqp1_1d(:)) )
           ztmp1_1d(:) = MAX(MIN(ztmp1_1d(:), 1.0_dp), 0.0_dp)
           ztmp1_1d(:) = nbetax*ztmp1_1d(:)   !ztt
           ixidx_1d(:) = INT(ztmp1_1d(:))

           ll1_1d(:) = (ixidx_1d(:) == nbetax)

           DO jl=1,kproma

              iqidx = iqidx_1d(jl)
              ixidx = ixidx_1d(jl)

              ztmp2_1d(jl) = (ztmp1_1d(jl)-ixidx_1d(jl)      )*tbetai0(iqidx,ixidx+1)       &
                           + (ixidx_1d(jl)+1._dp-ztmp1_1d(jl))*tbetai0(iqidx,ixidx)
              ztmp3_1d(jl) = (ztmp1_1d(jl)-ixidx_1d(jl)      )*tbetai1(iqidx,ixidx+1)       &
                           + (ixidx_1d(jl)+1._dp-ztmp1_1d(jl))*tbetai1(iqidx,ixidx)
           ENDDO

           ztmp4_1d(:) = MERGE(1._dp, ztmp2_1d(:), ll1_1d(:))  !SF zbetai0
           ztmp5_1d(:) = MERGE(1._dp, ztmp3_1d(:), ll1_1d(:))  !SF zbetai1

           ztmp1_1d(:) = -zxilb_1d(:)*zclcaux(:) 
           ztmp2_1d(:) = zqsec*zqp1_1d(:)

           zgent_1d(:) = pqm1(:,jk)                                &
                       - (zbetaqt(:)    - zbap1_1d(:))*ztmp5_1d(:) &
                       + (pbetass(:,jk) - zbap1_1d(:))*ztmp4_1d(:) &
                       -  pbetass(:,jk)

           zgent_1d(:) = MAX(zgent_1d(:), ztmp1_1d(:))
           zgent_1d(:) = MIN(zgent_1d(:),ztmp2_1d(:))   !limit to qv

           ztmp1_1d(:) = MAX(zepsec,zxilb_1d(:))
           ztmp1_1d(:) = zxib(:)/ztmp1_1d(:)
           ztmp1_1d(:) = MAX(MIN(ztmp1_1d(:),1.0_dp),0.0_dp) !SFzifrac
           ztmp2_1d(:) = 1._dp - ztmp1_1d(:)

           ztmp3_1d(:) = zgent_1d(:)/MAX(zclcaux(:), zeps)

           zgenti(:) = zgent_1d(:)*ztmp1_1d(:)
           zgentl(:) = zgent_1d(:)*ztmp2_1d(:)

           ztmp4_1d(:) = zxib(:) + ztmp3_1d(:)*ztmp1_1d(:)
           ztmp4_1d(:) = MAX(ztmp4_1d(:),0.0_dp)
           zxib(:)     = MERGE(ztmp4_1d(:), zxib(:), locc_1d(:))

           ztmp4_1d(:) = zxlb(:) + ztmp3_1d(:)*ztmp2_1d(:)
           ztmp4_1d(:) = MAX(ztmp4_1d(:), 0.0_dp)
           zxlb(:)     = MERGE(ztmp4_1d(:), zxlb(:), locc_1d(:))

           zxilb_1d(:) = zxib(:)+zxlb(:)
!
!       5.3 Deposition/sublimation of cloud ice and condensation/
!           evaporation of liquid water due to changes in water vapour
!           and temperature (advection, convection, turbulent mixing,
!           evaporation of rain, sublimation and melting of snow).
!           Translate PDF laterally to calculate cloud
!           after one timestep
!
           zqvdt_1d(:) = zqvdt_1d(:)-zgent_1d(:)
           zdtdt_1d(:) = zdtdt_1d(:)+zlvdcp(:)*zgentl(:)+zlsdcp(:)*zgenti(:)

           zqp1_1d(:) = pqm1(:,jk)+zqvdt_1d(:)
           zqp1_1d(:) = MAX(zqp1_1d(:),0.0_dp)

           ztp1_1d(:) = ptm1(:,jk)+zdtdt_1d(:)

           zdqsat_1d(:) = zdtdt_1d(:)                                                       &
                        + zclcaux(:)*( ztmst*zlc_1d(:)*pqte(:,jk)                           &
                                     + zlvdcp(:)*(zevp(:)+zxlevap(:)-zgentl(:))             &
                                     + zlsdcp(:)*(zsub(:)+zxievap(:)+zxisub(:)-zgenti(:)) )   !zdtdtstar

           zdqsat_1d(:) = zdqsat_1d(:)*zdqsdt_1d(:)/(1._dp+zclcaux(:)*zlc_1d(:)*zdqsdt_1d(:))

           ztmp1_1d(:) = (pbetass(:,jk)-zqvdt_1d(:)+zdqsat_1d(:)-zbap1_1d(:))/zbbap1_1d(:)
           ztmp1_1d(:) = MAX(MIN(ztmp1_1d(:), 1.0_dp), 0.0_dp)
           ztmp1_1d(:) = nbetax*ztmp1_1d(:)  !ztt
           ixidx_1d(:) = INT(ztmp1_1d(:))

           ll1_1d(:) = (ixidx_1d(:) == nbetax)

           DO jl=1,kproma

              iqidx = iqidx_1d(jl)
              ixidx = ixidx_1d(jl)

              ztmp2_1d(jl) = (ztmp1_1d(jl)-ixidx_1d(jl)      )*tbetai0(iqidx,ixidx+1)       &
                           + (ixidx_1d(jl)+1._dp-ztmp1_1d(jl))*tbetai0(iqidx,ixidx)
              ztmp3_1d(jl) = (ztmp1_1d(jl)-ixidx_1d(jl)      )*tbetai1(iqidx,ixidx+1)       &
                           + (ixidx_1d(jl)+1._dp-ztmp1_1d(jl))*tbetai1(iqidx,ixidx)
           ENDDO

           ztmp4_1d(:) = MERGE(1._dp, ztmp2_1d(:), ll1_1d(:))  !SF zbetai0
           ztmp5_1d(:) = MERGE(1._dp, ztmp3_1d(:), ll1_1d(:))  !SF zbetai1

           zqcdif_1d(:) = (zbetaqt(:)-pbetaa(:,jk))                            *(1._dp-ztmp5_1d(:)) &
                        + (pbetaa(:,jk)+zqvdt_1d(:)-pbetass(:,jk)-zdqsat_1d(:))*(1._dp-ztmp4_1d(:))

           zqcdif_1d(:) = MAX(0.0_dp, zqcdif_1d(:))
           zqcdif_1d(:) = zqcdif_1d(:)-zbetacl(:)

           ztmp1_1d(:) = -zxilb_1d(:)*zclcaux(:)
           ztmp2_1d(:) = zqsec*zqp1_1d(:)

            zqcdif_1d(:) = MAX(zqcdif_1d(:), ztmp1_1d(:))
            zqcdif_1d(:) = MIN(zqcdif_1d(:), ztmp2_1d(:))  ! limit to qv
          END IF !lcover  
        END IF ! ncctop  

        IF((.NOT.lcover) .OR. jk < ncctop) THEN
          zqcdif_1d(:) = (zqvdt_1d(:)-zdqsat_1d(:))*zclcaux(:)

          ztmp1_1d(:) = -zxilb_1d(:)*zclcaux(:)
          ztmp2_1d(:) = zqsec*zqp1_1d(:)

          zqcdif_1d(:) = MAX(zqcdif_1d(:), ztmp1_1d(:))
          zqcdif_1d(:) = MIN(zqcdif_1d(:), ztmp2_1d(:))  ! limit to qv

        END IF !lcover

        ll1_1d(:) = (zqcdif_1d(:) < 0._dp)     !SF cloud dissipation
        
        ztmp1_1d(:) = MAX(zepsec, zxilb_1d(:))
        ztmp1_1d(:) = zxib(:) / ztmp1_1d(:)
        ztmp1_1d(:) = MAX(MIN(ztmp1_1d(:), 1.0_dp), 0.0_dp) !zifrac

        ztmp1_1d(:) = MERGE(ztmp1_1d(:), 1._dp, ll1_1d(:))

        ztmp2_1d(:) = zqcdif_1d(:)*(1.0_dp-ztmp1_1d(:))
        zcnd(:) = MERGE(ztmp2_1d(:), 0._dp, ll1_1d(:))

        IF (nicnc <= 1) THEN
           zdep(:) = zqcdif_1d(:)*ztmp1_1d(:)
        ELSE
           zdep(:) = zqinucl(:,jk)*ztmp1_1d(:)
        ENDIF

        ll2_1d(:) = (.NOT. ll1_1d(:)) .AND. &  !SF no cloud dissipation
                    (.NOT. lo2_1d(:))          !SF condensation

        zdep(:) = MERGE(0._dp, zdep(:), ll2_1d(:))

!--- Included/changed for prognostic CDNC/IC scheme --------------------
!    Use standard condensation for empirical Lin & Leaitch approach and 
!    explicit condensation after Levkov et al. 1992 for the explicit
!    activation schemes that allow for supersaturation:

        IF (ncdnc == 1) THEN
           zcnd(:) = MERGE(zqcdif_1d(:), zcnd(:), ll2_1d(:))
        ELSE IF (ncdnc > 1) THEN
           ztmp1_1d(:) = 0.5_dp*api*zdw0*ztmst                  &
                       * swat(1:kproma,jk,jrow)*zcdnc(:,jk)*zclcaux(:) &
                       * zrho_rcp(:,jk)/zastbstw(:,jk)
           ztmp1_1d(:) = MIN(ztmp1_1d(:), zqcdif_1d(:))
           zcnd(:) = MERGE(ztmp1_1d(:), zcnd(:), ll2_1d(:))
        ENDIF

!--- End included for CDNC/IC scheme -----------------------------------
!
!
!       5.4 Accounting for cloud evaporation in clear air and
!           checking for supersaturation
!
        ztp1tmp(:) = ztp1_1d(:)+zlvdcp(:)*zcnd(:)+zlsdcp(:)*zdep(:)
        zqp1tmp(:) = zqp1_1d(:)-zcnd(:)-zdep(:)

        zxip1_1d(:) = zxised_2d(:,jk) + ztmst*zxite(:) - zxievap(:) + zgenti(:) + zdep(:)
        zxip1_1d(:) = MAX(zxip1_1d(:), 0._dp)

        lo2_1d(:)   = (ztp1tmp(:) < cthomi) .OR.                        &
                       (ztp1tmp(:) < tmelt .AND. zxip1_1d(:) > csecfrl  &
                       .AND. zsusatw_2d(:,jk) < zeps)

        ztmp1_1d(:) = 1000._dp*ztp1tmp(:)
        it_1d(:)    = NINT(ztmp1_1d(:))

        ll_look(:,jk)     = (it_1d(:)<jptlucu1 .OR. it_1d(:)>jptlucu2)
        IF (ANY(ll_look(1:kproma,jk))) lookupoverflow = .TRUE.

        it_1d(:)    = MAX(MIN(it_1d(:),jptlucu2),jptlucu1)

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

        ll1_1d(:) = (zes_1d(:) < 0.4_dp)  !SF LO

        zcor_1d(:) = 1._dp/(1._dp-vtmpc1*zes_1d(:))

        zqsp1tmp_1d(:) = zes_1d(:)*zcor_1d(:)
        zoversat_1d(:) = zqsp1tmp_1d(:)*0.01_dp
       
        zrhtest_1d(:) = pqm1(:,jk)/zqsm1_1d(:)
        zrhtest_1d(:) = MIN(zrhtest_1d(:), 1._dp)
        zrhtest_1d(:) = zrhtest_1d(:)*zqsp1tmp_1d(:)

!changed for cirrus scheme, Ulrike Lohmann, 17.12.2006
        zesw_1d(:) = zlucuaw_1d(:)/papp1(:,jk)
        zesw_1d(:) = MIN(zesw_1d(:), 0.5_dp)

        zcorw_1d(:) = 1._dp/(1._dp-vtmpc1*zesw_1d(:))

        zqsp1tmpw_1d(:) = zesw_1d(:)*zcorw_1d(:)
        zoversatw_1d(:) = 0.01_dp*zqsp1tmpw_1d(:)
!end changed for cirrus scheme

        zqst1_1d(:) = MERGE(zlucuap1_1d(:),zlucuawp1_1d(:),lo2_1d(:))
        zqst1_1d(:) = zqst1_1d(:)/papp1(:,jk)
        zqst1_1d(:) = MIN(zqst1_1d(:),0.5_dp)
        zqst1_1d(:) = zqst1_1d(:)/(1._dp-vtmpc1*zqst1_1d(:))

        zdqsdt_1d(:) = 1000._dp*(zqst1_1d(:)-zqsp1tmp_1d(:))

        zlc_1d(:) = MERGE(zlsdcp(:),zlvdcp(:),lo2_1d(:))

        ztmp1_1d(:) = zlc_1d(:)*zdqsdt_1d(:)
        ztmp2_1d(:) = zqsp1tmp_1d(:)*zcor_1d(:)*zlucub_1d(:)
        zlcdqsdt_1d(:) = MERGE(ztmp1_1d(:), ztmp2_1d(:), ll1_1d(:))

        zqcon_1d(:) = 1._dp/(1._dp+zlcdqsdt_1d(:))

        ztmp1_1d(:) = zqsp1tmpw_1d(:)+zoversatw_1d(:)
        ztmp1_1d(:) = MIN(ztmp1_1d(:),zqsp1tmp_1d(:)*1.3_dp)
        zqsp1tmphet_1d(:) = MERGE(ztmp1_1d(:), 0._dp, lo2_1d(:))

        ll1_1d(:) = (nicnc <= 1) 
! ll1_1d = false --> for cirrus scheme, Ulrike Lohmann, 17.12.2006

        ll2_1d(:) = (zqp1tmp(:) > (zqsp1tmp_1d(:)  + zoversat_1d(:) ) )
        ll3_1d(:) = (zqp1tmp(:) > (zqsp1tmpw_1d(:) + zoversatw_1d(:)) )
        ll4_1d(:) = (zqp1tmp(:) > zqsp1tmphet_1d(:))
        ll5_1d(:) = (ztp1tmp(:) >= cthomi)

        ztmp1_1d(:) = (zqp1tmp(:) - zqsp1tmp_1d(:)  - zoversat_1d(:) )*zqcon_1d(:)
        ztmp2_1d(:) = (zqp1tmp(:) - zqsp1tmpw_1d(:) - zoversatw_1d(:))*zqcon_1d(:)
        ztmp3_1d(:) = (zqp1tmp(:) - zqsp1tmphet_1d(:)                )*zqcon_1d(:)

        !SF ice cloud cases:         
        ll6_1d(:) = (lo2_1d(:) .AND. ll1_1d(:)         .AND. ll2_1d(:))                 &
                    .OR.                                                                &
                    (lo2_1d(:) .AND. (.NOT. ll1_1d(:)) .AND. ll2_1d(:) .AND. ll5_1d(:))
        ztmp4_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll6_1d(:))

        IF ( nicnc == 2 ) THEN 
           ll6_1d(:) = lo2_1d(:) .AND. ll3_1d(:) .AND. (.NOT. ll5_1d(:))          &    !SB: T<238 (cirrus regime)
                       .AND. (.NOT. lhet)               !SB: when lhet=false, since (.NOT.lhet)=TRUE, i.e. homog. freezing
           ztmp4_1d(:) = MERGE(ztmp2_1d(:), ztmp4_1d(:), ll6_1d(:))
    
           ll6_1d(:) = lo2_1d(:) .AND. ll4_1d(:) .AND. (.NOT. ll5_1d(:))          &    !SB: T<238 (cirrus regime)
                       .AND. lhet                       !SB: when lhet=true, i.e. heterog. freezing  
           ztmp4_1d(:) = MERGE(ztmp3_1d(:), ztmp4_1d(:), ll6_1d(:)) 
        END IF
    
        IF ( nicnc == 3 ) THEN
           !adding IF (nicnc == 3) => (.NOT. ll1_1d(:)) can be removed as well, besides the variable lhet
           !moreover, the condition ll3 for Si above water can be neglected, as in Kuebbeler.
           ll6_1d(:) = lo2_1d(:) .AND. ll4_1d(:) .AND. (.NOT. ll5_1d(:))        
           ztmp4_1d(:) = MERGE(ztmp3_1d(:), ztmp4_1d(:), ll6_1d(:)) 
        END IF

        zdep(:) = zdep(:) + ztmp4_1d(:)
       
        !SF water cloud cases:
        ll6_1d(:) = (.NOT. lo2_1d(:)) .AND. ll2_1d(:)
        ztmp4_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll6_1d(:))
        zcnd(:)     = zcnd(:) + ztmp4_1d(:)

        !SF final corrections for zdep and zcnd:
        ztmp5_1d(:) = zqp1_1d(:)-zrhtest_1d(:)
        ztmp5_1d(:) = MAX(ztmp5_1d(:), 0._dp)

        ll1_1d(:) = (zdep(:)        > 0._dp        )
        ll2_1d(:) = (zcnd(:)        > 0._dp        )
        ll3_1d(:) = (zqp1tmp(:)     < zrhtest_1d(:))
        ll4_1d(:) = (zqsp1tmp_1d(:) <= zqsm1_1d(:) )

        ll5_1d(:) = lo2_1d(:) .AND. ll1_1d(:) .AND. ll3_1d(:) .AND. ll4_1d(:)
        zdep(:) = MERGE(ztmp5_1d(:), zdep(:), ll5_1d(:)) 

        ll5_1d(:) = (.NOT. lo2_1d(:)) .AND. ll2_1d(:) .AND. ll3_1d(:) .AND. ll4_1d(:)
        zcnd(:) = MERGE(ztmp5_1d(:), zcnd(:), ll5_1d(:))
!
!       5.5 Change of in-cloud water due to deposition/sublimation and
!           condensation/evaporation (input for cloud microphysics)
!
        zrelhum_1d(:) = zqp1tmp(:)/zqsp1tmp_1d(:)

        ztmp1_1d(:) = zdep(:) + zgenti(:)  
        ztmp1_1d(:) = MAX(ztmp1_1d(:), 0._dp) !SF deposition zdepos
       
        ztmp2_1d(:) = zcnd(:) + zgentl(:)
        ztmp2_1d(:) = MAX(ztmp2_1d(:), 0._dp) !SF condensation zcond 

        ztmp3_1d(:) = zxib(:)+zdep(:)/MAX(zclcaux(:), zeps) 
        ztmp3_1d(:) = MAX(ztmp3_1d(:), 0._dp)

        ztmp4_1d(:) = zxlb(:)+zcnd(:)/MAX(zclcaux(:), zeps)
        ztmp4_1d(:) = MAX(ztmp4_1d(:), 0._dp)

        zxib(:) = MERGE(ztmp3_1d(:), zxib(:), locc_1d(:))
        zxlb(:) = MERGE(ztmp4_1d(:), zxlb(:), locc_1d(:))

        ll1_1d(:) = (.NOT. locc_1d(:))                                       &
                  .AND. ( (ztmp1_1d(:) > 0._dp) .OR. (ztmp2_1d(:) > 0._dp) )        

        ztmp3_1d(:) = MAX(MIN(zrelhum_1d(:), 1.0_dp), 0.01_dp)
        zclcaux(:)  = MERGE(ztmp3_1d(:), zclcaux(:), ll1_1d(:))

        ztmp3_1d(:) = ztmp1_1d(:) / MAX(zclcaux(:), zeps)
        ztmp4_1d(:) = ztmp2_1d(:) / MAX(zclcaux(:), zeps)
      
        zxib(:) = MERGE(ztmp3_1d(:), zxib(:), ll1_1d(:))
        zxlb(:) = MERGE(ztmp4_1d(:), zxlb(:), ll1_1d(:))

        !SF re-define locc_1d which was no longer true:
        locc_1d(:) = (zclcaux(:) > zeps)

        !SF update cdnc:
        ll1_1d(:) = locc_1d(:) .AND. (zxlb(:) > cqtmin)
        ll2_1d(:) = ll1_1d(:) .AND. (zcdnc(:,jk) <= cdncmin) !if there was no previous nucleation

        ztmp1_1d(:) = pcdncact(:,jk) - zcdnc(:,jk)
        ztmp1_1d(:) = MAX(0._dp, ztmp1_1d(:))
        zqlnuc(:,jk)    = MERGE(ztmp1_1d(:), zqlnuc(:,jk), ll2_1d(:))

        ztmp1_1d(:)      = MERGE(ztmp1_1d(:), 0._dp, ll2_1d(:))
        zcdnc(:,jk)      = zcdnc(:,jk) + ztmp1_1d(:)
        qnuc(1:kproma,jk,jrow)  = qnuc(1:kproma,jk,jrow) + zdt*ztmp1_1d(:)

        ztmp1_1d(:) = MAX(zcdnc(:,jk), cdncmin)
        zcdnc(:,jk) = MERGE(ztmp1_1d(:), cqtmin, ll1_1d(:)) 

        !SF update icnc:
        ll1_1d(:) = locc_1d(:) .AND. (zxib(:) > cqtmin) 
        ll2_1d(:) = ll1_1d(:)  .AND. (zicnc(:,jk) <= zicemin)

        IF (nicnc <= 1) THEN
           ztmp1_1d(:) = 0.75_dp/(api*zrhoice)*zrho(:,jk)*zxib(:)/zrid_2d(:,jk)**3
        ELSE
           ztmp1_1d(:) = MIN(znicex(:,jk),(zap(:,jk)*1.e6_dp))
        ENDIF

        zicnc(:,jk) = MERGE(ztmp1_1d(:), zicnc(:,jk), ll2_1d(:))

        ztmp1_1d(:) = MAX(zicnc(:,jk), zicemin)
        zicnc(:,jk) = MERGE(ztmp1_1d(:), cqtmin, ll1_1d(:)) 

        ztp1tmp(:) = ztp1_1d(:) + zlvdcp(:)*zcnd(:) + zlsdcp(:)*zdep(:)
!
!     ------------------------------------------------------------------
!       6.    Freezing of cloud water
!
!       6.1   Freezing of all cloud water for T < 238 K
!
        ll1_1d(:) = (ztp1tmp(:) <= cthomi)

        ztmp1_1d(:) = zfrl(:,jk) + zxlb(:)*zclcaux(:)
        zfrl(:,jk)  = MERGE(ztmp1_1d(:), zfrl(:,jk), ll1_1d(:))

        ztmp1_1d(:) = zxib(:)+zxlb(:)
        zxib(:)     = MERGE(ztmp1_1d(:), zxib(:), ll1_1d(:))

        zxlb(:)     = MERGE(0._dp, zxlb(:), ll1_1d(:))

!--- Included for prognostic CDNC/IC scheme ----------------------------
        ztmp1_1d(:)     = zcdnc(:,jk)-cdncmin
        ztmp1_1d(:)     = MAX(ztmp1_1d(:), 0._dp)
        ztmp2_1d(:)     = qfre(1:kproma,jk,jrow) - zdt*ztmp1_1d(:)
        qfre(1:kproma,jk,jrow) = MERGE(ztmp2_1d(:), qfre(1:kproma,jk,jrow), ll1_1d(:))

        ztmp2_1d(:) = zicnc(:,jk) + ztmp1_1d(:)
        zicnc(:,jk) = MERGE(ztmp2_1d(:), zicnc(:,jk), ll1_1d(:))

        zcdnc(:,jk) = MERGE(cqtmin, zcdnc(:,jk), ll1_1d(:)) 

!--- End included for CDNC/IC scheme -----------------------------------
!
!       6.2   Freezing of cloud water between 238 and 273 K
!
        lo_1d(:) =     (zxlb(:)     > cqtmin  ) &
                 .AND. (ztp1tmp(:)  < tmelt   ) &
                 .AND. (ztp1tmp(:)  > cthomi  ) & 
                 .AND. (zcdnc(:,jk) >= cdncmin) &
                 .AND. locc_1d(:)

!---Changed for prognostic CDNC/IC scheme ------------------------------
!   (Replaced pacdnc by zcdnc which is set above)

     !*************************************************************
     ! Compute the heterogenous freezing using:
     !
     ! - one spectrum included in Barahona & Nenes, 2009b (if limm_BN09 = TRUE) 
     !   for deposition and condensation/immersion nucleation 
     !   OR 
     !   the parameterization by Lohmann & Diehl, 2006 (if limm_BN09 = FALSE)
     !   (Corinna included for contact/immersion freezing by dust and soot)
     !
     ! - the parameterization by Lohmann & Diehl, 2006
     !   for contact nucleation
     !
     ! - the parameterization for thermopheresis (if lthermo = TRUE)
     !
     ! Notes: numbers of equations refer to Lohmann & Diehl, 2006 
     !*************************************************************

         !---Some computations needed for contact nucleation:
         !SB: Cunningham correction factors for ki, ai, ci modes
         !SB:  Cc    = 1 + 1.26 * lambda / radius_mode * (Po/P) * (T/To)           
         ztmp1_1d(:) = 1._dp + 1.26_dp*6.6E-8_dp / (zrwetki_2d(:,jk)+zeps) * (101325._dp/papp1(:,jk))  & !SF ccbcki
                                                 * (ztp1tmp(:)/tmelt)
         ztmp2_1d(:) = 1._dp + 1.26_dp*6.6E-8_dp / (zrwetai_2d(:,jk)+zeps) * (101325._dp/papp1(:,jk))  & !SF ccduai
                                                 * (ztp1tmp(:)/tmelt)
         ztmp3_1d(:) = 1._dp + 1.26_dp*6.6E-8_dp / (zrwetci_2d(:,jk)+zeps) * (101325._dp/papp1(:,jk))  & !SF ccduci
                                                 * (ztp1tmp(:)/tmelt)

         !SB: viscosity of air (eta)
         !SB: eta      = 10^(-5) * ( 1.718+0.0049*(T-To) - 1.2^10(-5)*(T-To)^2 )
         zetaair_1d(:) = 1.e-5_dp  &
                       * ( 1.718_dp + 0.0049_dp*(ztp1tmp(:)-tmelt)                      &
                                    - 1.2e-5_dp*(ztp1tmp(:)-tmelt)*(ztp1tmp(:)-tmelt) )

!Ulrike: for thermophoresis
         zkair_1d(:)    = 4.1867e-3_dp*(5.69_dp + 0.017_dp*(ztp1tmp(:)-tmelt)) ! eq. 13-18a P&K
!end thermophoresis

         !SB: Brownian aerosol diffusivity (eq. 2) for ki-BC, ai-DU, ci-DU modes
         ll1_1d(:)   = (zrwetki_2d(:,jk) < zeps)
         zdfarbcki_1d(:) = ak * ztp1tmp(:) * ztmp1_1d(:) &
                         / ( 6._dp*api*zetaair_1d(:)*(zrwetki_2d(:,jk)+zeps) )
         zdfarbcki_1d(:) = MERGE(0._dp, zdfarbcki_1d(:), ll1_1d(:))

         ll2_1d(:)   = (zrwetai_2d(:,jk) < zeps)
         zdfarduai_1d(:) = ak * ztp1tmp(:) * ztmp2_1d(:) &
                         / ( 6._dp*api*zetaair_1d(:)*(zrwetai_2d(:,jk)+zeps) )
         zdfarduai_1d(:) = MERGE(0._dp, zdfarduai_1d(:), ll2_1d(:))

         ll3_1d(:)   = (zrwetci_2d(:,jk) < zeps)
         zdfarduci_1d(:) = ak * ztp1tmp(:) * ztmp3_1d(:) &
                         / ( 6._dp*api*zetaair_1d(:)*(zrwetci_2d(:,jk)+zeps) )
         zdfarduci_1d(:) = MERGE(0._dp, zdfarduci_1d(:), ll3_1d(:))

!Ulrike: for thermophoresis
         IF (lthermo) THEN
             zknbcki_1d(:) = 7.37_dp*ztp1tmp(:)/(2.88e5_dp*papp1(:,jk)*(zrwetki_2d(:,jk)+zeps))
             zknbcki_1d(:) = MERGE(1._dp, zknbcki_1d(:), ll1_1d(:))
 
             zftbcki_1d(:) = 0.4_dp &
                           * ( 1._dp + zknbcki_1d(:)*(1.45_dp + 0.4_dp*EXP(-1._dp/zknbcki_1d(:)))) &
                           * (zkair_1d(:)+2.5_dp*zknbcki_1d(:)*zkbc)                                  &
                           / ( (1._dp+3._dp*zknbcki_1d(:))                                         &
                             * (2._dp*zkair_1d(:)+zkbc*(5._dp*zknbcki_1d(:)+1._dp)))
 
             zftbcki_1d(:) = MERGE(0._dp, zftbcki_1d(:), ll1_1d(:))

             zknduai_1d(:) = 7.37_dp*ztp1tmp(:)/(2.88e5_dp*papp1(:,jk)*(zrwetai_2d(:,jk)+zeps))
             zknduai_1d(:) = MERGE(1._dp, zknduai_1d(:), ll2_1d(:))
 
             zftduai_1d(:) = 0.4_dp &
                          * ( 1._dp + zknduai_1d(:)*(1.45_dp + 0.4_dp*EXP(-1._dp/zknduai_1d(:)))) &
                          * (zkair_1d(:)+2.5_dp*zknduai_1d(:)*zkdu)                                  &
                          / ( (1._dp+3._dp*zknduai_1d(:))                                         &
                            * (2._dp*zkair_1d(:)+zkdu*(5._dp*zknduai_1d(:)+1._dp)))
 
             zftduai_1d(:) = MERGE(0._dp, zftduai_1d(:), ll2_1d(:))

             zknduci_1d(:) = 7.37_dp*ztp1tmp(:)/(2.88e5_dp*papp1(:,jk)*(zrwetci_2d(:,jk)+zeps))
             zknduci_1d(:) = MERGE(1._dp, zknduci_1d(:), ll3_1d(:))
 
             zftduci_1d(:) = 0.4_dp &
                          * ( 1._dp + zknduci_1d(:)*(1.45_dp + 0.4_dp*EXP(-1._dp/zknduci_1d(:)))) &
                          * (zkair_1d(:)+2.5_dp*zknduci_1d(:)*zkdu)                                  &
                          / ( (1._dp+3._dp*zknduci_1d(:))                                         &
                            * (2._dp*zkair_1d(:)+zkdu*(5._dp*zknduci_1d(:)+1._dp)))
 
             zftduci_1d(:) = MERGE(0._dp, zftduci_1d(:), ll3_1d(:))
         ENDIF
!end thermophoresis

         DO jl=1,kproma

            !SB: Some computations needed for contact nucleation AND thermophoresis:
            zfracdusol     = MIN(ndusol_strat(jl,jk,jrow)/(pcdncact(jl,jk)      +zeps), 1._dp)
            zfracduinsolai = MIN(nduinsolai(jl,jk,jrow)  /(naerinsol(jl,jk,jrow)+zeps), 1._dp)
            zfracduinsolci = MIN(nduinsolci(jl,jk,jrow)  /(naerinsol(jl,jk,jrow)+zeps), 1._dp)
            zfracbcsol     = MIN(nbcsol_strat(jl,jk,jrow)/(pcdncact(jl,jk)      +zeps), 1._dp)
            zfracbcinsol   = MIN(nbcinsol(jl,jk,jrow)    /(naerinsol(jl,jk,jrow)+zeps), 1._dp)

            zradl = (0.75_dp*zxlb(jl) * zrho(jl,jk)           &
                  / (api*rhoh2o*zcdnc(jl,jk)))**(1._dp/3._dp)

            zf1 = 4._dp*api*zradl*zcdnc(jl,jk)*zrho_rcp(jl,jk)

            !SB: temperature dependence of individual species
            zfrzcntdu = MIN(1._dp,MAX(0._dp,-(0.1014_dp*(ztp1tmp(jl)-tmelt)+0.3277_dp)))  ! montmorillonite
            !zfrzcntdu = MIN(1._dp,MAX(0._dp,-(0.1007_dp*(ztp1tmp(jl)-tmelt)+0.6935_dp)))  ! kaolinite
            !zfrzcntbc = MIN(1._dp,MAX(0._dp,-(0.0614_dp*(ztp1tmp(jl)-tmelt)+0.5730_dp)))
            zfrzcntbc = 0._dp ! disable BC contact freezing

            !-------- Contact nucleation by dust and soot -----------------------------
            !SB: Contact freezing (similar to eq(1))
            zfrzcnt = zxlb(jl) / zcdnc(jl,jk) * zrho(jl,jk) * zf1                                &
                    * ( zfrzcntdu * (zdfarduai_1d(jl)*zfracduinsolai+zdfarduci_1d(jl)*zfracduinsolci)   &
                      + zfrzcntbc *  zdfarbcki_1d(jl)*zfracbcinsol                                 ) &
                    * ( zcdnc(jl,jk)+zicnc(jl,jk) )

            zfrzcnt        = zxlb(jl)*(1._dp-EXP(-zfrzcnt/MAX(zxlb(jl), cqtmin)*ztmst))
            !-------- End Contact nucleation ------------------------------------------

!--- ulrike: add thermophoresis (UL, 8.1.09)
            zfrzthermo=0._dp

            IF (lthermo) THEN

               lo2_1d(jl) = (ztp1tmp(jl) < cthomi) .OR.                      &
                            (ztp1tmp(jl) < tmelt .AND. zxip1_1d(jl) > csecfrl &
                            .AND. zsusatw_2d(jl,jk) < zeps)

               zqsm1      = MERGE(tlucua(itm1_look(jl,jk)),tlucuaw(itm1_look(jl,jk)),lo2_1d(jl)) &
                          / papm1(jl,jk)
               zqsm1      = MIN(zqsm1,0.5_dp)
               zqsm1      = zqsm1/(1._dp-vtmpc1*zqsm1)

               it1_1d(jl)  = NINT(ztp1tmp(jl)*1000._dp)
               it1_1d(jl)  = MAX(MIN(it1_1d(jl),jptlucu2),jptlucu1)
               IF (it1_1d(jl) < jptlucu1 .OR. it1_1d(jl) > jptlucu2) lookupoverflow = .TRUE.

               zqst1       = MERGE(tlucua(it1_1d(jl)),tlucuaw(it1_1d(jl)),lo2_1d(jl))/papm1(jl,jk)
               zqst1       = MIN(zqst1,0.5_dp)
               zqst1       = zqst1/(1._dp-vtmpc1*zqst1)
               zdqsdtime   = (zqst1-zqsm1)
 
               zdeltatemp = MAX( (zrho(jl,jk)*alv*(zdep(jl)+zdqsdtime))                &
                                 / (ztmst*zkair_1d(jl)*4._dp*api*zcdnc(jl,jk)*zradl+zeps), 0._dp )

               zf2        = zkair_1d(jl) / papp1(jl,jk) * zdeltatemp

               zfrzthermo = zf1*zf2                                      &
                          * ( zfrzcntbc* zftbcki_1d(jl)*zfracbcinsol     &
                            + zfrzcntdu*(zftduai_1d(jl)*zfracduinsolai   &
                                        +zftduci_1d(jl)*zfracduinsolci)) &
                          * ( zcdnc(jl,jk)+zicnc(jl,jk) )*ztmst          &
                          / ( zcdnc(jl,jk)*zrho(jl,jk) )

               zfrzthermo = zxlb(jl)*(1._dp-EXP(-zfrzthermo))

            ENDIF
!--- Ulrike: end thermophoresis

            !Only for diagnostic, store the contributions ...
            !... by contact and thermophoresis nucl.
            ztmp1_1d(jl) = MAX( 0.0_dp, (zfrzcnt + zfrzthermo) )
            ztmp4_1d(jl) = zcdnc(jl,jk)*ztmp1_1d(jl)/(zxlb(jl)+zeps)  !from [kg/kg] to [1/m3] 
            newIC_cnt_therm(jl,jk) = ztmp4_1d(jl)                !newly formed IC in via cnt & thermo nucl. [1/m3]

          IF (limm_BN09) THEN

            !-------- Immersion nucleation by BN09b at T>238.15 K-------------------
            !initialize output of BN09: 
            smaxice = 0.0_dp
            sc_ice  = 0.0_dp
            nlim    = 0.0_dp
            nhet    = 0.0_dp
            nice    = 0.0_dp
            dice    = 0.0_dp
            sigwparc= 0.0_dp

            !---------------------------------------
            !--- Call to the BN09 scheme
            !---------------------------------------
            IF (ptm1(jl,jk) > cthomi) THEN      !IF T > 238, mixed-phase regime
 
               CALL ice_activate (ptm1(jl,jk), paphm1(jl,jk), sigw(jl,jk), zicncq(jl,jk),     &  
                                  ndrop(jl,jk),   dsulf(jl,jk),                               &  
                                  ndust(jl,jk,:), ddust(jl,jk,:), sigdust(:),                 &  
                                  norg(jl,jk),    dorg(jl,jk),    sigorg,                     &  
                                  nsoot(jl,jk),   dsoot(jl,jk),   sigsoot,                    &  
                                  nbio(jl,jk),    dbio(jl,jk),    sigbio,                     &  
                                  nsolo(jl,jk,:), dsolo(jl,jk,:), sigsolo(:),                 &  
                                  smaxice, nlim, sc_ice, nhet, nice, dice, sigwparc)

            END IF   !condition T>238

            !Variables saved in the new channel cloud_ice_BN09
            smaxice_immBN(jl,jk)  = smaxice
            sc_ice_immBN(jl,jk)   = sc_ice
            nice_immBN(jl,jk)     = nice
            dice_immBN(jl,jk)     = dice
            sigwpre_immBN(jl,jk)  = sigwparc 
            !-------- End Immersion nucleation by BN09b -----------------------------

            !---- Only for diagnostic, store the contributions ...
            !... by immersion nucl. via BN09  
            newIC_imm(jl,jk) = nice                     !newly formed IC in via immersion nucl. [1/m3]
            !... TOTAL NEW ice crystals
            ztmp2_1d(jl) = ztmp4_1d(jl) + nice
            newIC_mix(jl,jk) = ztmp2_1d(jl)             !newly formed IC in mixed clouds [1/m3]
            !---- End

            ztmp1_1d(jl) = zfrzcnt + zfrzthermo                          !SB: [kg/kg]
            ztmp1_1d(jl) = MAX(0.0_dp,MIN(ztmp1_1d(jl),zxlb(jl)))        !SF zfrl surrogate
            ztmp4_1d(jl) = zcdnc(jl,jk)*ztmp1_1d(jl)/(zxlb(jl)+zeps)     !SB: from [kg/kg] to [1/m3]  
            ztmp2_1d(jl) = ztmp4_1d(jl) + nice  

            ztmp3_1d(jl) = (sum(ndust(jl,jk,:))+norg(jl,jk)+nsoot(jl,jk)+nbio(jl,jk)+ndrop(jl,jk)+ &  
                            sum(nsolo(jl,jk,:))) - zicnc(jl,jk)                                         !mz_sb_20170710

          ELSE   !limm_BN09=False

            !-------- Immersion nucleation by Lohmann & Diehl (2006) ---------------
            znaimmdu  = 32.3_dp*zfracdusol     ! montmorillonite     !SB: Table 2 from Diehl and Wurzler (2004) 
            !znaimmdu  = 6.15E-2_dp*zfracdusol  ! kaolinite 

            znaimmbc  = 2.91E-3_dp*zfracbcsol                        !SB: Table 2 from Diehl and Wurzler (2004)

            zomega = pvervel(jl,jk) - 1.33_dp*SQRT(ptkem1(jl,jk))*zrho(jl,jk)*g
            ztte   = zomega / cpd *zrho_rcp(jl,jk)

            zfrzimm = -(znaimmdu+znaimmbc)*zrho(jl,jk)/rhoh2o*EXP(tmelt-ztp1tmp(jl))*MIN(ztte,0._dp) 
            zfrzimm = zxlb(jl)*(1._dp-EXP(-zfrzimm*zxlb(jl)/zcdnc(jl,jk)*ztmst))        !SB: [kg/kg]
            !-------- End Immersion nucleation by Lohmann & Diehl (2006) ------------

            ztmp1_1d(jl) = zfrzcnt + zfrzimm + zfrzthermo                               !SB: [kg/kg]

            !Only for diagnostic, store the contributions ...
            !... by immersion nucl. via LD06  
            ztmp2_1d(jl) = MAX(0.0_dp, zfrzimm)
            ztmp3_1d(jl) = zcdnc(jl,jk)*ztmp2_1d(jl)/(zxlb(jl)+zeps)   !from [kg/kg] to [1/m3] 
            newIC_imm(jl,jk) = ztmp3_1d(jl)                        !newly formed IC in via immersion nucl. [1/m3]
            !... TOTAL NEW ice crystals
            ztmp2_1d(jl) = MAX(0.0_dp,ztmp1_1d(jl))  
            ztmp3_1d(jl) = zcdnc(jl,jk)*ztmp2_1d(jl)/(zxlb(jl)+zeps)   !from [kg/kg] to [1/m3] 
            newIC_mix(jl,jk) = ztmp3_1d(jl)                        !newly formed IC in mixed clouds [1/m3]

            ztmp1_1d(jl) = MAX(0.0_dp,MIN(ztmp1_1d(jl),zxlb(jl))) !SF zfrl surrogate

            ztmp2_1d(jl) = zcdnc(jl,jk)*ztmp1_1d(jl)/(zxlb(jl)+zeps)

            ztmp3_1d(jl) = ndusol_strat(jl,jk,jrow) + nbcsol_strat(jl,jk,jrow) + nduinsolai(jl,jk,jrow) &
                         + nduinsolci(jl,jk,jrow) + nbcinsol(jl,jk,jrow) - zicnc(jl,jk)

          END IF   !limm_BN09 

            ztmp2_1d(jl) = MIN(ztmp2_1d(jl), ztmp3_1d(jl))
            ztmp2_1d(jl) = MAX(ztmp2_1d(jl), 0._dp)  !SF zfrln surrogate

            ! Store the final values of NEW IC formed in the mixed-phase regime 
            newIC_mix_fin(jl,jk) = ztmp2_1d(jl)
            !

         ENDDO !SF end loop jl

         zfrl(:,jk) = MERGE(ztmp1_1d(:), zfrl(:,jk) , lo_1d(:))

         ztmp1_1d(:) = MIN(ztmp2_1d(:), zcdnc(:,jk)-cdncmin)
         ztmp1_1d(:) = MAX(ztmp1_1d(:), 0._dp)
         zfrln(:)    = MERGE(ztmp1_1d(:), zfrln(:), lo_1d(:))

         ztmp1_1d(:) = zcdnc(:,jk)-zfrln(:)
         ztmp1_1d(:) = MAX(ztmp1_1d(:), cqtmin)
         zcdnc(:,jk) = MERGE(ztmp1_1d(:), zcdnc(:,jk), lo_1d(:))

         ztmp1_1d(:) = zicnc(:,jk)+zfrln(:)
         ztmp1_1d(:) = MAX(ztmp1_1d(:), cqtmin)
         zicnc(:,jk) = MERGE(ztmp1_1d(:), zicnc(:,jk), lo_1d(:))

!--- End included for CDNC/IC scheme -----------------------------------

         ztmp1_1d(:) = zxlb(:)-zfrl(:,jk)
         zxlb(:)     = MERGE(ztmp1_1d(:), zxlb(:), lo_1d(:))

         ztmp1_1d(:) = zxib(:)+zfrl(:,jk)
         zxib(:)     = MERGE(ztmp1_1d(:), zxib(:), lo_1d(:))

         ztmp1_1d(:) = zfrl(:,jk)*zclcaux(:)
         zfrl(:,jk)     = MERGE(ztmp1_1d(:), zfrl(:,jk), lo_1d(:))

!--------corinna: Bergeron-Findeisen-Process:
         ll1_1d(:) =     lo_1d(:)                  &
                   .AND. locc_1d(:)                &
                   .AND. (zdep(:) > 0._dp)         &
                   .AND. (zxlb(:) > 0._dp)         &
                   .AND. (0.01_dp*zvervx_2d(:,jk) < zvervmax_1d(:))      

         ztmp1_1d(:) = ztmst_rcp*zxlb(:)*zclcaux(:) !SF zzevp
      
         ztmp2_1d(:) = pxlte(:,jk)-ztmp1_1d(:)
         pxlte(:,jk) = MERGE(ztmp2_1d(:), pxlte(:,jk), ll1_1d(:))

         ztmp2_1d(:) = pxite(:,jk)+ztmp1_1d(:)
         pxite(:,jk) = MERGE(ztmp2_1d(:), pxite(:,jk), ll1_1d(:))

         ztmp2_1d(:) = ptte(:,jk)+(zlsdcp(:)-zlvdcp(:))*ztmp1_1d(:)
         ptte(:,jk)  = MERGE(ztmp2_1d(:), ptte(:,jk), ll1_1d(:))
          
         zcdnc(:,jk) = MERGE(cqtmin, zcdnc(:,jk), ll1_1d(:))

         ztmp2_1d(:) = zxib(:)+zxlb(:)
         zxib(:)     = MERGE(ztmp2_1d(:), zxib(:), ll1_1d(:))

         zxlb(:) = MERGE(0._dp, zxlb(:), ll1_1d(:))

!-------End corinna: Bergeron-Findeisen-Process

        IF (lookupoverflow) THEN
          status_string = 'lookuperror: cdnc - cloud (1)'
          RETURN
        ENDIF
!
!     ------------------------------------------------------------------
!       7.  Cloud physics and precipitation fluxes at the surface
!
!ham_ps: cdir circumvents bug in sxf90 compiler

        zclcstar_1d(:) = MIN(zclcaux(:), zclcpre(:))
        zauloc_1d(:)   = 3./5000._dp*zdz_2d(:,jk)
        zauloc_1d(:)   = MAX(MIN(zauloc_1d(:), clmax), clmin)

        ll1_1d(:) = (knvb(:) >= jbmin) .AND. &
                    (knvb(:) <= jbmax) .AND. &
                    (pvervel(:,jk) > 0._dp)

        ll2_1d(:) = (jk == knvb(:)  ) .OR. &
                    (jk == knvb(:)+1)

        ll3_1d(:) = ll1_1d(:) .AND. ll2_1d(:) .AND. lonacc

        zauloc_1d(:) = MERGE(0._dp, zauloc_1d(:), ll3_1d(:))

        zxlb(:) = MAX(zxlb(:),1.e-20_dp)
        zxib(:) = MAX(zxib(:), 1.e-20_dp)

! liquid water and snow content are stored 
! before the reduction by outfalling rain
! (necessary for nucleation scavenging)
        plwc(:,jk) = zxlb(:)
        piwc(:,jk) = zxib(:)      

!---Included for in-cloud scavenging (Philip Stier, 25/11/03):----------
        zmlwc(:,jk) = zxlb(:)
        zmiwc(:,jk) = zxib(:)

!---End Included for scavenging-----------------------------------------
!
!---  Calculate the rain and snow water content in kg/kg from the rain and snow flux
!
       ll1_1d(:) = (zclcpre(:) > zeps )
       ll2_1d(:) = ll1_1d(:) .AND. (zrfl(:) > cqtmin)
       ll3_1d(:) = ll1_1d(:) .AND. (zsfl(:) > cqtmin)

       ztmp1_1d(:) = ( MAX(zrfl(:), cqtmin)/(12.45_dp*MAX(zclcpre(:),zeps)*SQRT(zqrho_2d(:,jk))) )**(8._dp/9._dp)
       ztmp2_1d(:) = ( MAX(zsfl(:), cqtmin)/(cvtfall *MAX(zclcpre(:),zeps)                     ) )**(1._dp/1.16_dp)

       zxrp1_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll2_1d(:))
       zxsp1_1d(:) = MERGE(ztmp2_1d(:), 0._dp, ll3_1d(:))

!
!       7.1   Warm clouds: Coalescence processes after Beheng (1994):
!             Autoconversion of cloud droplets and collection of cloud
!             droplets by falling rain. Accretion of cloud droplets by
!             falling snow (zsacl) is calculated under 7.2
!

        ll1_1d(:) = locc_1d(:)                .AND. &
                    (zxlb(:) > cqtmin)        .AND. &
                    (zcdnc(:,jk) >= cdncmin)

        IF (nauto == 2) THEN
!          Autoconversion rate from Khairoutdinov and Kogan, 2000

           ztmp1_1d(:) = ccraut*1350._dp*(1.e-6_dp*zcdnc(:,jk))**(-1.79_dp)
     
           ztmp1_1d(:) = zxlb(:) * (  1._dp &
                                   - (1._dp + ztmst*zexm1_1*ztmp1_1d(:)*zxlb(:)**zexm1_1)**zexp_1)

           ztmp1_1d(:) = MIN(zxlb(:), ztmp1_1d(:))
           zraut_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))
           
           ztmp1_1d(:) = zxlb(:) - zraut_1d(:)
           ztmp2_1d(:) = zxlb(:) !SF keeps zxlb for later use
           zxlb(:)     = MERGE(ztmp1_1d(:), zxlb(:), ll1_1d(:))
 
!--- zrac1 is formed by accretion with rain from above
!--- zrac2 is formed by accretion with newly formed rain inside the grid box

           ztmp1_1d(:) = -3.7_dp*ztmst*zxrp1_1d(:)
           ztmp1_1d(:) = EXP(ztmp1_1d(:))
           ztmp1_1d(:) = zxlb(:)*(1._dp-ztmp1_1d(:))
           zrac1_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))

           zxlb(:) = zxlb(:) - zrac1_1d(:)

           ztmp1_1d(:) = -3.7_dp*ztmst*zauloc_1d(:)*zrho(:,jk)*zraut_1d(:)
           ztmp1_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))
           ztmp1_1d(:) = zxlb(:)*(1._dp-EXP(ztmp1_1d(:)))
           zrac2_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))

           zxlb(:) = zxlb(:) - zrac2_1d(:)

           zrpr(:) = zrpr(:) + zclcaux(:)     * (zraut_1d(:)+zrac2_1d(:)) &
                             + zclcstar_1d(:) *  zrac1_1d(:)

!---Included for in-cloud scavenging (Philip Stier, 26/11/03):----------
           ztmp1_1d(:)    = zraut_1d(:)+zrac1_1d(:)+zrac2_1d(:)
           zmratepr(:,jk) = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))
!---End Included for scavenging-----------------------------------------

!--- Included for prognostic CDNC/IC scheme ----------------------------
!--- Autoconversion also changes the number of cloud droplets (zrprn)

           ztmp1_1d(:) = (zraut_1d(:)+zrac1_1d(:)+zrac2_1d(:))/(ztmp2_1d(:)+zeps) !SF ztmp2_1d=zxlb+zraut+zrac1+zrac2
           zrprn(:)    = MERGE(ztmp1_1d(:), zrprn(:), ll1_1d(:))

           ll2_1d(:) = ll1_1d(:)          .AND. &
                       (zxlb(:) > cqtmin)

           ztmp1_1d(:) = MERGE(cdncmin, 0._dp, ll2_1d(:))
           ztmp1_1d(:) = zcdnc(:,jk)-ztmp1_1d(:)
           ztmp2_1d(:) = zcdnc(:,jk)*zrprn(:)

           ztmp3_1d(:) = MIN(ztmp1_1d(:), ztmp2_1d(:))
           zrprn(:)    = MERGE(ztmp3_1d(:), zrprn(:), ll1_1d(:))

           ztmp1_1d(:) = zcdnc(:,jk)-zrprn(:)
           ztmp1_1d(:) = MAX(ztmp1_1d(:), cqtmin)
           zcdnc(:,jk) = MERGE(ztmp1_1d(:), zcdnc(:,jk), ll1_1d(:))

!--- End included for CDNC/IC scheme ------------------------------------
!--- End included alternative autoconversion parameterisation ----------
!--- Changed for alternative autoconversion parameterisation -----------

        ELSE   !SF( nauto ==1)
!          Beheng (1994) - ECHAM 5 standard

!---Changed for prognostic CDNC/IC scheme ------------------------------
!   (Replaced pacdnc by zcdnc which is set above)
           ztmp1_1d(:) = ccraut*1.2e27_dp * zrho_rcp(:,jk)                    & 
                                          * (zcdnc(:,jk)*1.e-6_dp)**(-3.3_dp) &
                                          * (zrho(:,jk)*1.e-3_dp)**4.7_dp
           zraut_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))
!--- End changed for prognostic CDNC/IC scheme -------------------------

           zraut_1d(:) = zxlb(:) * (  1._dp &
                                   - (1._dp + ztmst*zexm1_2*zraut_1d(:)*zxlb(:)**zexm1_2)**zexp_2)

           zraut_1d(:) = MIN(zxlb(:), zraut_1d(:))
  
!--- Included for prognostic CDNC/IC scheme ----------------------------
           ztmp1_1d(:) = 7.7e9_dp * zraut_1d(:) * zrho(:,jk)                      !SF zrautn
           ztmp2_1d(:) = 1.289e10_dp * 1.e-6_dp * ztmst * (zrho(:,jk)*zxlb(:))**2 !SF zself (1.e-6 comes
                                                                                  ! from a unit bug fix)
           ztmp3_1d(:) = ztmp1_1d(:) + ztmp2_1d(:)
           ztmp3_1d(:) = MIN(ztmp3_1d(:),zcdnc(:,jk))
           zrautself_1d(:) = MERGE(ztmp3_1d(:), 0._dp, ll1_1d(:))

           ztmp1_1d(:) = zcdnc(:,jk)-zrautself_1d(:)
           ztmp1_1d(:) = MAX(ztmp1_1d(:), cqtmin)
           zcdnc(:,jk) = MERGE(ztmp1_1d(:), zcdnc(:,jk), ll1_1d(:))

!--- End included for CDNC/IC scheme -----------------------------------

           ztmp1_1d(:) = zxlb(:) - zraut_1d(:)
           zxlb(:)     = MERGE(ztmp1_1d(:), zxlb(:), ll1_1d(:))
!
!--- zrac1 is formed by accretion with rain from above
!--- zrac2 is formed by accretion with newly formed rain inside the grid box
!
           zrac1_1d(:) = -6._dp*ztmst*zxrp1_1d(:)
           zrac1_1d(:) = EXP(zrac1_1d(:))
           zrac1_1d(:) = zxlb(:)*(1._dp-zrac1_1d(:))

           ztmp1_1d(:) = zxlb(:) - zrac1_1d(:)
           ztmp2_1d(:) = zxlb(:)  !SF keeps zxlb for later use
           zxlb(:)     = MERGE(ztmp1_1d(:), zxlb(:), ll1_1d(:))
 
           zrac2_1d(:) = -6._dp*ztmst*zauloc_1d(:)*zrho(:,jk)*zraut_1d(:)
           zrac2_1d(:) = EXP(zrac2_1d(:))
           zrac2_1d(:) = zxlb(:)*(1._dp-zrac2_1d(:))

           ztmp1_1d(:) = zxlb(:)-zrac2_1d(:)
           zxlb(:)     = MERGE(ztmp1_1d(:), zxlb(:), ll1_1d(:))

           ztmp1_1d(:) = zrpr(:) + zclcaux(:)     * (zraut_1d(:)+zrac2_1d(:)) & 
                                 + zclcstar_1d(:) * zrac1_1d(:)
           zrpr(:)     = MERGE(ztmp1_1d(:), zrpr(:), ll1_1d(:))

!---Included for in-cloud scavenging (Philip Stier, 26/11/03):----------
           ztmp1_1d(:)    = zraut_1d(:) + zrac1_1d(:) + zrac2_1d(:)
           zmratepr(:,jk) = MERGE(ztmp1_1d(:), zmratepr(:,jk), ll1_1d(:)) 
!---End Included for scavenging-----------------------------------------

!--- Included for prognostic CDNC/IC scheme ----------------------------
!--- Autoconversion also changes the number of cloud droplets (zrprn)
           ztmp1_1d(:) = (zrac1_1d(:)+zrac2_1d(:))/(ztmp2_1d(:)+zeps) !SF ztmp2_1d=zxlb+zrac1+zrac2
           ztmp1_1d(:) = zcdnc(:,jk)*ztmp1_1d(:)
           ztmp1_1d(:) = MIN(ztmp1_1d(:), zcdnc(:,jk)) !SF zraccn

           ztmp2_1d(:) = zcdnc(:,jk)-ztmp1_1d(:)
           ztmp2_1d(:) = MAX(ztmp2_1d(:), cqtmin)
           zcdnc(:,jk) = MERGE(ztmp2_1d(:), zcdnc(:,jk), ll1_1d(:))

           ztmp2_1d(:) = zrautself_1d(:) + ztmp1_1d(:)
           zrprn(:)    = MERGE(ztmp2_1d(:), zrprn(:), ll1_1d(:))
!--- End included for CDNC/IC scheme -----------------------------------
        ENDIF 

!       7.2  Cold clouds:
!            Conversion of cloud ice to snow after Levkov et al. 1992:
!            Aggregation of ice crystals to snow assuming plates (zsaut) and accretion of ice
!            by falling snow. (zsaci)
!            Accrection of cloud droplets by falling snow. (zsacl)
!            Effective radius of ice crystals assuming plates (Lohmann, ACPD, 2007)

        ll1_1d(:) = locc_1d(:) .AND. (zxib(:) > cqtmin)

        ztmp1_1d(:) = 0.5e4_dp * ( 1000._dp/0.0376_dp*zxib(:)*zrho(:,jk)/zicnc(:,jk) )**0.302 !plate !SF zrieff
        ztmp1_1d(:) = MIN(MAX(ztmp1_1d(:), ceffmin), ceffmax) !SF zrieff

        ztmp2_1d(:) = 5113188._dp+2809._dp*ztmp1_1d(:)**3
        ztmp2_1d(:) = SQRT(ztmp2_1d(:))
        ztmp2_1d(:) = -2261._dp + ztmp2_1d(:)  !SF zrih
        
        ztmp3_1d(:) = 1.e-6_dp*ztmp2_1d(:)**(1._dp/3._dp)
        zris_1d(:)  = MERGE(ztmp3_1d(:), 1._dp, ll1_1d(:))   !SF 1. could be whatever, just not 0.
       
!--- temperature dependent collision efficiency (zcolleffi)

        ztmp1_1d(:)    = 0.025_dp*(ztp1tmp(:)-tmelt)
        ztmp1_1d(:)    = EXP(ztmp1_1d(:))
        zcolleffi_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:)) 

        zc1_1d(:) = 17.5_dp / crhoi * zrho(:,jk) * zqrho_2d(:,jk)**0.33_dp

        ztmp1_1d(:) = -6._dp / zc1_1d(:) * LOG10(1.e4_dp*zris_1d(:)) !SF zdt2
        ztmp1_1d(:) = ccsaut / ztmp1_1d(:)

        zsaut_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))
        zsaut_1d(:) = zxib(:) * (1._dp - 1._dp/(1._dp+zsaut_1d(:)*ztmst*zxib(:)))

        ztmp1_1d(:) = zxib(:) - zsaut_1d(:)
        zxibold_1d(:) = zxib(:) !SF store zxib for later use at the end of 7.2
        zxib(:)     = MERGE(ztmp1_1d(:), zxib(:), ll1_1d(:))

        zsaci2_1d(:) = 0.0_dp
        zsacl2_1d(:) = 0.0_dp

!--- Included for prognostic CDNC/IC scheme ----------------------------
!    accretion also affects the number of cloud droplets
        zsacl2in_1d(:)  = 0.0_dp      ! 
!--- End included for CDNC/IC scheme -----------------------------------

        ztmp1_1d(:) = zauloc_1d(:)*zrho(:,jk)*zsaut_1d(:)
        zxsp2_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))  !snow that is formed inside the grid box

        zxsp_1d(:) = zxsp1_1d(:) + zxsp2_1d(:)

        ll2_1d(:) = ll1_1d(:)               .AND. &
                    (zxsp_1d(:)  >  cqtmin) .AND. &
                    (zxlb(:)     >  cqtmin) .AND. &
                    (zcdnc(:,jk) >= cdncmin)

!--- The riming of snow with droplets is calculated assuming planar snow flakes (Lohmann, JAS, 2004)
!--- It depends on the droplet (zudrop) and snow flake (zusnow) fall velocity, the Stokes
!--- number (zstokes) and the Reynolds number (zrey)

        zscnc_1d(:) = 0.5_dp*zicnc(:,jk)  !SF note: not put in the ll2_1d() condition as in orig code
        IF (jk > 1) THEN 
           zscnc_1d(:) = MAX(zicemin, MIN(zsprn(:,jk-1), zscnc_1d(:)))
        ELSE
           zscnc_1d(:) = MAX(zicemin, zscnc_1d(:)) !SF added protection for non-zero 
        ENDIF

        ztmp1_1d(:) = ( 6._dp*zpirho_rcp*zrho(:,jk)*zxlb(:)/zcdnc(:,jk) )**(1._dp/3._dp)
        ztmp1_1d(:) = MAX(ztmp1_1d(:), 1.e-6_dp) !SF zdw

        zudrop_1d(:) = 1.19e4_dp*2500._dp*ztmp1_1d(:)**2*(1.3_dp*zrho_rcp(:,jk))**0.35_dp

        zdplanar_1d(:) = 1.e3_dp/3.8e-4_dp*zxsp_1d(:)/zscnc_1d(:)
        zdplanar_1d(:) = SQRT(zdplanar_1d(:))
        zdplanar_1d(:) = 1.e-2_dp*zdplanar_1d(:)
        zdplanar_1d(:) = MAX(20.e-6_dp, zdplanar_1d(:)) !SF note: not put in the ll2_1d() condition as in orig code (useless)
     
        zusnow_1d(:) = 2.34_dp * (100._dp*zdplanar_1d(:))**0.3_dp &
                               * (1.3_dp*zrho_rcp(:,jk))**0.35_dp

        zstokes_1d(:) = 2._dp*g_rcp*(zusnow_1d(:)-zudrop_1d(:))*zudrop_1d(:)/zdplanar_1d(:)
        zstokes_1d(:) = MAX(zstokes_1d(:), cqtmin)

        zrey_1d(:) = zrho(:,jk)*zdplanar_1d(:)*zusnow_1d(:)/zviscos_2d(:,jk)
        zrey_1d(:) = MAX(zrey_1d(:),cqtmin)

        ll3_1d(:) = (zrey_1d(:) <=  5._dp)
        ll4_1d(:) = (zrey_1d(:) >   5._dp) .AND. (zrey_1d(:) <  40._dp)
        ll5_1d(:) = (zrey_1d(:) >= 40._dp)
     
        ztmp1_1d(:)   = 5.52_dp*zrey_1d(:)**(-1.12_dp)
        ztmp2_1d(:)   = 1.53_dp*zrey_1d(:)**(-0.325_dp)
        zstcrit_1d(:) = 1._dp
        zstcrit_1d(:) = MERGE(ztmp1_1d(:), zstcrit_1d(:), ll3_1d(:))
        zstcrit_1d(:) = MERGE(ztmp2_1d(:), zstcrit_1d(:), ll4_1d(:))

        zcsacl_1d(:) = 0.2_dp * ( LOG10(zstokes_1d(:)) - LOG10(zstcrit_1d(:)) - 2.236_dp )**2
        zcsacl_1d(:) = MIN(zcsacl_1d(:), 1._dp-cqtmin)
        zcsacl_1d(:) = MAX(zcsacl_1d(:), 0._dp)
        zcsacl_1d(:) = SQRT(1._dp - zcsacl_1d(:))

        ll6_1d(:) = ll5_1d(:) .AND. (zstokes_1d(:) <= 0.06_dp)
        ll7_1d(:) = ll5_1d(:) .AND. (zstokes_1d(:) >  0.06_dp) .AND. (zstokes_1d(:) <= 0.25_dp)
        ll8_1d(:) = ll5_1d(:) .AND. (zstokes_1d(:) >  0.25_dp) .AND. (zstokes_1d(:) <= 1.00_dp)

        WHERE (ll6_1d(:))
              zcsacl_1d(:) = 1.034_dp*zstokes_1d(:)**1.085_dp
        ELSEWHERE (ll7_1d(:))
              zcsacl_1d(:) = 0.787_dp*zstokes_1d(:)**0.988_dp
        ELSEWHERE (ll8_1d(:))
              zcsacl_1d(:) = 0.7475_dp*LOG10(zstokes_1d(:))+0.65_dp
        ELSEWHERE (ll5_1d(:))
              zcsacl_1d(:) = (zstokes_1d(:)+1.1_dp)**2/(zstokes_1d(:)+1.6_dp)**2
        ENDWHERE

        zcsacl_1d(:) = MAX(MIN(zcsacl_1d(:), 1._dp), 0.01_dp)
        zcsacl_1d(:)  = MERGE(zcsacl_1d(:), 0._dp, ll2_1d(:))

        ztmp1_1d(:)  = zcons4*zxsp_1d(:)**0.8125_dp !SF zlamsm
        ztmp1_1d(:)  = api*cn0s*3.078_dp*ztmp1_1d(:)*zqrho_2d(:,jk)**0.5_dp
        zsaci2_1d(:) = MERGE(ztmp1_1d(:), zsaci2_1d(:), ll2_1d(:)) 

        ztmp1_1d(:)  = -ztmst*zsaci2_1d(:)*zcsacl_1d(:)
        ztmp1_1d(:)  = EXP(ztmp1_1d(:))
        ztmp1_1d(:)  = zxlb(:)*(1._dp-ztmp1_1d(:))

!--- Included for prognostic CDNC/IC scheme ----------------------------
        zsacl2in_1d(:) = MERGE(ztmp1_1d(:), zsacl2in_1d(:), ll2_1d(:))
!--- End included for CDNC/IC scheme -----------------------------------

        ztmp2_1d(:) = zxlb(:)-ztmp1_1d(:)
        ztmp3_1d(:) = zxlb(:)  !SF keeps zxlb for later use
        zxlb(:)     = MERGE(ztmp2_1d(:), zxlb(:), ll2_1d(:))

        ztmp2_1d(:) = zclcaux(:)*ztmp1_1d(:)
        zsacl2_1d(:)   = MERGE(ztmp2_1d(:), zsacl2_1d(:), ll2_1d(:))

!SF end ll2_1d condition:
!         (ll1_1d(:) .AND. (zxsp_1d(:) > cqtmin) .AND. (zxlb(:) > cqtmin) .AND. (zcdnc(:,jk) > cdncmin))     

        ll2_1d(:) = ll1_1d(:)             .AND. &
                    (zxsp_1d(:) > cqtmin) .AND. &
                    (zxib(:)    > cqtmin)

        ztmp1_1d(:) = zcons4*zxsp_1d(:)**0.8125_dp
        ztmp1_1d(:) = api*cn0s*3.078_dp*ztmp1_1d(:)*zqrho_2d(:,jk)**0.5_dp
        ztmp1_1d(:) = -ztmst*ztmp1_1d(:)*zcolleffi_1d(:)
        ztmp1_1d(:) = EXP(ztmp1_1d(:))
        ztmp1_1d(:) = zxib(:) * (1._dp-ztmp1_1d(:))

        zsaci2_1d(:) = MERGE(ztmp1_1d(:), zsaci2_1d(:), ll2_1d(:))
        
        zxib(:) = zxib(:)-zsaci2_1d(:)

!SF end ll2_1d condition:
!     (ll1_1d(:) .AND. (zxsp_1d(:) > cqtmin) .AND. (zxib(:) > cqtmin))

        zsacl(:) = MERGE(zsacl2_1d(:), zsacl(:), ll1_1d(:))

        ztmp1_1d(:) = zspr(:) + zclcaux(:)*( zsaut_1d(:)+zsaci2_1d(:) )
        zspr(:)     = MERGE(ztmp1_1d(:), zspr(:), ll1_1d(:))

!---Included for in-cloud scavenging (Philip Stier, 25/11/03):----------
        ztmp1_1d(:)    = zsaut_1d(:)+zsaci2_1d(:)
        zmrateps(:,jk) = MERGE(ztmp1_1d(:), zmrateps(:,jk), ll1_1d(:))
!---End Included for scavenging-----------------------------------------

!--- Included for prognostic CDNC/IC scheme ----------------------------
        ll2_1d(:) = (zxlb(:) > cqtmin) 

        ztmp1_1d(:) = zcdnc(:,jk)*zsacl2in_1d(:)/(ztmp3_1d(:)+zeps) !SF ztmp3_1d =zxlb+zsacl2in
        ztmp1_1d(:) = MIN(ztmp1_1d(:), zcdnc(:,jk)-cdncmin)
        ztmp1_1d(:) = MAX(ztmp1_1d(:), 0._dp) 

        ztmp2_1d(:) = zcdnc(:,jk)*zsacl2in_1d(:)/(ztmp3_1d(:)+zeps) !SF ztmp3_1d = zxlb+zsacl2in
        ztmp2_1d(:) = MIN(ztmp2_1d(:), zcdnc(:,jk))

        ztmp3_1d(:)  = MERGE(ztmp1_1d(:), ztmp2_1d(:), ll2_1d(:))
        zsacln(:) = MERGE(ztmp3_1d(:), zsacln(:), ll1_1d(:))

!SF end ll2_1d condition (zxlb > cqtmin)

        ztmp1_1d(:)     = zcdnc(:,jk)-zsacln(:)
        zcdnc(:,jk)     = MERGE(ztmp1_1d(:), zcdnc(:,jk), ll1_1d(:))
        zmsnowacl(:,jk) = MERGE(zsacl2in_1d(:), zmsnowacl(:,jk), ll1_1d(:))

!       secondary ice crystal production (zsecprod) after Levkov et al. 1992
!       sink for snow, source for ice crystals
!uls    included size depedent accretion rate (Lohmann, JAS, 2004)

        zsecprod_1d(:) = 0._dp

        ll2_1d(:) = ll1_1d(:)               .AND. &
                    (zxsp_1d(:) > zepsec)   .AND. &
                    (zxlb(:)    > zepsec)   .AND. &
                    (ztp1tmp(:) > 265.2_dp) .AND. &
                    (ztp1tmp(:) < 270.2_dp)


        ztmp1_1d(:) = ( 6._dp*zpirho_rcp*zrho(:,jk)*zxlb(:)/zcdnc(:,jk) )**(1._dp/3._dp)
        ztmp1_1d(:) = MAX(ztmp1_1d(:), 1.e-6_dp) !SF zdw

        zudrop_1d(:) = 1.19e4_dp*(50._dp*ztmp1_1d(:))**2*(1.3_dp*zrho_rcp(:,jk))**0.35_dp
 
        zstokes_1d(:) = 2._dp*g_rcp*(zusnow_1d(:)-zudrop_1d(:))*zudrop_1d(:)/zdplanar_1d(:)
        zstokes_1d(:) = MAX(zstokes_1d(:), cqtmin)

        ztmp1_1d(:) = 0.2_dp * ( LOG10(zstokes_1d(:)) - LOG10(zstcrit_1d(:)) - 2.236_dp )**2
        ztmp1_1d(:) = MIN(ztmp1_1d(:), 1._dp-cqtmin)
        ztmp1_1d(:) = MAX(ztmp1_1d(:), 0._dp)
        ztmp1_1d(:) = SQRT(1._dp - ztmp1_1d(:))

        ll6_1d(:) = ll5_1d(:) .AND. (zstokes_1d(:) <= 0.06_dp)
        ll7_1d(:) = ll5_1d(:) .AND. (zstokes_1d(:) >  0.06_dp) .AND. (zstokes_1d(:) <= 0.25_dp)
        ll8_1d(:) = ll5_1d(:) .AND. (zstokes_1d(:) >  0.25_dp) .AND. (zstokes_1d(:) <= 1.00_dp)

        WHERE (ll6_1d(:))
              ztmp1_1d(:) = 1.034_dp*zstokes_1d(:)**1.085_dp
        ELSEWHERE (ll7_1d(:))
              ztmp1_1d(:) = 0.787_dp*zstokes_1d(:)**0.988_dp
        ELSEWHERE (ll8_1d(:))
              ztmp1_1d(:) = 0.7475_dp*LOG10(zstokes_1d(:))+0.65_dp
        ELSEWHERE (ll5_1d(:))
              ztmp1_1d(:) = (zstokes_1d(:)+1.1_dp)**2/(zstokes_1d(:)+1.6_dp)**2
        ENDWHERE

        ztmp1_1d(:)  = MAX(MIN(ztmp1_1d(:), 1._dp), 0.01_dp)
        zcsacl_1d(:) = MERGE(ztmp1_1d(:), zcsacl_1d(:), ll2_1d(:))

        ztmp1_1d(:) = zcons5*zxsp_1d(:)**0.875_dp !SF zlams2

        ztmp2_1d(:) = cn0s * 0.831_dp * api / zmw0                      &
                    * zcsacl_1d(:) * zrho(:,jk) * zxlb(:) * ztmp1_1d(:) &
                    * ( g*crhosno / (0.75_dp*zcdi*zrho(:,jk)) )**0.5_dp   !SF zj
    
        ztmp2_1d(:) = MAX(0.00285_dp*ztmp2_1d(:), 0._dp)  !SF zpn

        ztmp2_1d(:) = ztmst*zmi0*ztmp2_1d(:)*zrho_rcp(:,jk)
        ztmp3_1d(:) = zxsp_1d(:)*zrho_rcp(:,jk)
        ztmp2_1d(:) = MIN(ztmp3_1d(:), ztmp2_1d(:))
        ztmp2_1d(:) = MAX(ztmp2_1d(:), 0._dp)

        zsecprod_1d(:) = MERGE(ztmp2_1d(:), zsecprod_1d(:), ll2_1d(:))

        ztmp1_1d(:) = zxib(:)+zsecprod_1d(:)
        zxib(:)     = MERGE(ztmp1_1d(:), zxib(:), ll2_1d(:))

        ztmp1_1d(:) = zspr(:)-zclcstar_1d(:)*zsecprod_1d(:)
        ztmp1_1d(:) = MAX(ztmp1_1d(:), 0._dp)
        zspr(:)     = MERGE(ztmp1_1d(:), zspr(:), ll2_1d(:))
      
        ztmp1_1d(:)    = zmrateps(:,jk)-zsecprod_1d(:)
        zmrateps(:,jk) = MERGE(ztmp1_1d(:), zmrateps(:,jk), ll2_1d(:))

        ! storing of snow production before it is transformed into a flux
        prate_s(:,jk) = zsaut_1d(:) + zsaci2_1d(:)

!SF end ll2_1d condition (secondary ice production)

!SF end ll1_1d condition (locc and zxib(jl) > cqtmin) 
!
!--- Also change the number of ice crystals due to the break-up of snow flakes
!
        ll1_1d(:) = locc_1d(:)               .AND. &
                    (zxib(:)     > zepsec  ) .AND. &
                    (zicnc(:,jk) >= zicemin)
       
        ztmp1_1d(:) = zxibold_1d(:)
        ztmp1_1d(:) = MAX(ztmp1_1d(:), 0._dp) !SF zxibold

        ztmp2_1d(:) = zicnc(:,jk) * (zsaci2_1d(:)+zsaut_1d(:)) / (ztmp1_1d(:)+zeps) !SF zsprn1

        ztmp3_1d(:) = 0.5_dp * ztmst * zc1_1d(:) * zicnc(:,jk) * zxib(:) !SF zself

        ztmp4_1d(:) = zmi0_rcp * zrho(:,jk) * zsecprod_1d(:) !SF zsecprodn

        ztmp5_1d(:) = ztmp2_1d(:)+ztmp3_1d(:)-ztmp4_1d(:) 
        ztmp5_1d(:) = MIN(ztmp5_1d(:), zicnc(:,jk))    !SF zsprnself
        zsprn(:,jk) = MERGE(ztmp5_1d(:), zsprn(:,jk), ll1_1d(:))

        ztmp1_1d(:) = zicnc(:,jk)-zsprn(:,jk)
        ztmp1_1d(:) = MAX(ztmp1_1d(:), cqtmin)
        zicnc(:,jk) = MERGE(ztmp1_1d(:), zicnc(:,jk), ll1_1d(:))
!--- End included for CDNC/IC scheme -----------------------------------

!       7.3 Updating precipitation fluxes. In the lowest layer (klev),
!           the sedimentation sink of cloud ice is balanced
!           by precipitation at the surface (through 'zzdrs').
!           Fraction of precipitating clouds (zclcpre) used for the
!           calculation of evaporation/sublimation of rain/snow in
!           the next layer
!
        zzdrr_1d(:)    = zcons2*zdp_2d(:,jk)*zrpr(:)
        zzdrs_1d(:)    = zcons2*zdp_2d(:,jk)*(zspr(:)+zsacl(:))

        IF (jk .EQ. klev) THEN
           zzdrs_1d(:) = zzdrs_1d(:)+zxiflux(:)
           ztmp1_1d(:) = zcons2*zdp_2d(:,jk)/(zlsdcp(:)-zlvdcp(:)) &
                               *MAX(0._dp, (ztp1tmp(:)-tmelt))       !SF zcons
           ztmp2_1d(:) = MIN(zxsec*zzdrs_1d(:), ztmp1_1d(:))         !SF zsnmlt
           zzdrr_1d(:) = zzdrr_1d(:)+ztmp2_1d(:)
           zzdrs_1d(:) = zzdrs_1d(:)-ztmp2_1d(:)
           zsmlt(:)    = zsmlt(:)+ztmp2_1d(:)/(zcons2*zdp_2d(:,jk))
        END IF

        zpretot_1d(:)  = zrfl(:)+zsfl(:)
        zpredel_1d(:)  = zzdrr_1d(:)+zzdrs_1d(:)


        ll1_1d(:)      = (zpretot_1d(:) > zpredel_1d(:))

        zclcpre(:)     = MERGE(zclcpre(:), zclcaux(:), ll1_1d(:))

        zpresum_1d(:)  = zpretot_1d(:)+zpredel_1d(:)
        
        ll1_1d(:) = (zpresum_1d(:) > cqtmin)

        ztmp1_1d(:) = (zclcaux(:)*zpredel_1d(:) + zclcpre(:)*zpretot_1d(:)) / MAX(zpresum_1d(:), cqtmin)
        ztmp1_1d(:) = MIN(ztmp1_1d(:), 1.0_dp)
        ztmp1_1d(:) = MAX(ztmp1_1d(:), 0.0_dp)

        zclcpre(:) = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))

!---Included for in-cloud scavenging (Philip Stier, 28/03/01):----------
!   Corrected by Junhua Zhang, Philip Stier (01/2004)

        ll1_1d(:) = (zclcpre(:) > zepsec)

        ztmp1_1d(:) = (zrfl(:)+zzdrr_1d(:))/MAX(zclcpre(:), zepsec)
        ztmp2_1d(:) = (zsfl(:)+zzdrs_1d(:))/MAX(zclcpre(:), zepsec)
        ztmp3_1d(:) = (zcons2*zdp_2d(:,jk)*zevp(:))/MAX(zclcpre(:), zepsec)
        ztmp4_1d(:) = (zcons2*zdp_2d(:,jk)*zsub(:))/MAX(zclcpre(:), zepsec)

        zfrain(:,jk)  = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))
        zfsnow(:,jk)  = MERGE(ztmp2_1d(:), 0._dp, ll1_1d(:))
        zfevapr(:,jk) = MERGE(ztmp3_1d(:), 0._dp, ll1_1d(:))
        zfsubls(:,jk) = MERGE(ztmp4_1d(:), 0._dp, ll1_1d(:))

!---End Included for scavenging-----------------------------------------
        ! rain and snow flux considering incoming rain, melting of snow, 
        ! droplet evaporation / sublimation , but no new production of rain or snow 
        ! in that layer....
        ! (neccessary for impaction scavenging)
        pfrain_no(:,jk)   = zrfl(:) - zcons2*zdp_2d(:,jk)*zevp(:)  
        pfsnow_no(:,jk)   = zsfl(:) - zcons2*zdp_2d(:,jk)*zsub(:)
        ! precipitating cloud cover of this layer is used for the next lower layer 
        ! to estimate the part of the cloud cover in which rain impacts
        pr_cover(:,jk) = zclcpre(:)

        zrfl(:)       = zrfl(:)+zzdrr_1d(:)-zcons2*zdp_2d(:,jk)*zevp(:)
        zsfl(:)       = zsfl(:)+zzdrs_1d(:)-zcons2*zdp_2d(:,jk)*zsub(:)

!
!     ------------------------------------------------------------------
!       8.    Updating tendencies of t, q, xl, xi and final cloud cover
!
!
!       8.10   Cloud cover scheme tendencies
!
        IF (jk >= ncctop) THEN
!
!          Source terms from convection
!          Skewness:
!
!-------------------Added by Junhua Zhang for CONV Micro-----------------------
          IF (ncvmicro > 0) THEN
             zconvskew(:)   = cbeta_cs * (pxtecl(:,jk)+pxteci(:,jk)+pqtec(:,jk)) &
                                       / pbetass(:,jk)
          ELSE
             zconvskew(:)   = cbeta_cs * (pxtec(:,jk)+pqtec(:,jk)) &
                                       / pbetass(:,jk)
          ENDIF
!-------------------------------------end-------------------------------------

           ztmp1_1d(:)  = zdtime_rcp*(cbeta_pq_max-pxskew(:,jk))
           zconvskew(:) = MIN(zconvskew(:), ztmp1_1d(:))
!
!          Convective width now diagnosed, assuming 'a' unchanged:
!
           ll1_1d(:) = (pqm1(:,jk) >= pbetass(:,jk))

           ztmp1_1d(:) = pxskew(:,jk) + zdtime*zconvskew(:) !SF zskewp1
           ztmp2_1d(:) = zwide(:) * (cbeta_pq+ztmp1_1d(:))   &
                                  / (cbeta_pq+pxskew(:,jk))  !SF zbbap1
           ztmp3_1d(:) = zdtime_rcp*(ztmp2_1d(:)-zwide(:))

           zconvvar(:) = MERGE(ztmp3_1d(:), 0._dp, ll1_1d(:))

!
!       8.11 Simple linearized effect of microphysics on skewness
!
           ll1_1d(:) = (pbetaa(:,jk) < pbetass(:,jk)) .AND. &
                       (pbetab(:,jk) > pbetass(:,jk))

           ztmp1_1d(:) = ztmst*(zxlte(:)+zxite(:))                &
                       - zrpr(:)-zsacl(:)-zspr(:)+zcnd(:)+zdep(:) &
                       + zgenti(:)+zgentl(:)                       

           ztmp1_1d(:) = - ztmp1_1d(:) / MAX(zepsec,zbetacl(:))
           ztmp1_1d(:) = MIN(1._dp, ztmp1_1d(:))
           ztmp1_1d(:) = MAX(0._dp, ztmp1_1d(:))
           ztmp1_1d(:) = (pbetass(:,jk)-pbetab(:,jk)) * ztmp1_1d(:)  !SF zmdelb

           ztmp2_1d(:) = (pbetab(:,jk)+ztmp1_1d(:)-pbetaa(:,jk))  &
                       * cbeta_pq/(zbetaqt(:)-pbetaa(:,jk))       &
                       - cbeta_pq
           ztmp2_1d(:) = MAX(MIN(ztmp2_1d(:), cbeta_pq_max), cbeta_pq)

           ztmp3_1d(:) = zdtime_rcp*(ztmp2_1d(:)-pxskew(:,jk))
           ztmp3_1d(:) = MIN(0._dp, ztmp3_1d(:))

           zmicroskew(:) = MERGE(ztmp3_1d(:), zmicroskew(:), ll1_1d(:))

!
!       8.2   New skewness and variance
!
           zxskewte(:)    = zconvskew(:)                 &
                          + zmicroskew(:) + zturbskew(:)

           zxvarte(:)     = zconvvar(:) + zturbvar(:)
!
           ztmp1_1d(:)    = pxvar(:,jk)  + zdtime*zxvarte(:) !SF zvarp1
           ztmp2_1d(:)    = pxskew(:,jk) + zdtime*zxskewte(:) !SF zskew1
!
           pxskew(:,jk)   = MAX(MIN(ztmp2_1d(:), cbeta_pq_max), cbeta_pq)

           ztmp3_1d(:)    = zbetaqt(:)*(1._dp+pxskew(:,jk)/cbeta_pq) !SF zvarmx

           pxvar(:,jk)    = MAX(MIN(ztmp1_1d(:), ztmp3_1d(:)), zvartg(:))
!

        ENDIF !jk >= ncctop
!
!       8.3   Tendencies of thermodynamic variables
!             Attn: The terms zxisub and zximlt do not appear in
!                   pxite because these processes have already been
!                   included in pxite via changes in cloud ice
!                   sedimentation (see 3.1, 3.2 and 4)
!
        pqte(:,jk)  = pqte(:,jk)                                     &
                    + ztmst_rcp * ( -zcnd(:)-zgentl(:)+zevp(:)+zxlevap(:)    &
                                  -  zdep(:)-zgenti(:)+zsub(:)+zxievap(:)    &
                                  +  zxisub(:) )

        ptte(:,jk)  = ptte(:,jk)                                                                     &
                    + ztmst_rcp * ( zlvdcp(:) * (zcnd(:)+zgentl(:)-zevp(:)-zxlevap(:))               &
                                  + zlsdcp(:) * (zdep(:)+zgenti(:)-zsub(:)-zxievap(:)                &
                                                -zxisub(:))+(zlsdcp(:)-zlvdcp(:))                    &
                                              * (-zsmlt(:)-zimlt(:)-zximlt(:)+zfrl(:,jk)+zsacl(:)) )

        ztmp1_1d(:) = pxlte(:,jk) + zxlte(:)

        ztmp2_1d(:) = zimlt(:) + zximlt(:) - zfrl(:,jk) - zrpr(:)       &
                    - zsacl(:) + zcnd(:)   + zgentl(:)  - zxlevap(:)

        zxlp1_1d(:) = pxlm1(:,jk) + ztmst * ztmp1_1d(:) + ztmp2_1d(:)

        pxlte(:,jk) = ztmp1_1d(:) + ztmst_rcp * ztmp2_1d(:)

        ztmp1_1d(:) = pxite(:,jk) + zxite(:)
   
        ztmp2_1d(:) = zfrl(:,jk) - zspr(:) + zdep(:) + zgenti(:) - zxievap(:)

        zxip1_1d(:) = pxim1(:,jk) + ztmst * ztmp1_1d(:)  + ztmp2_1d(:) 

        pxite(:,jk) = ztmp1_1d(:) + ztmst_rcp * ztmp2_1d(:)

!
!--- Included for prognostic CDNC/IC scheme ----------------------------

        !--- Calculate new total tendency of CDNC:
!qqq+ this is indeed the total tendency, since zcdnc is initialised 
!     above including the tendency, and here the difference to t-1
!     is calculated ...
        pxtte(:,jk,1) = ztmst_rcp * (zcdnc(:,jk)* m_air/1000._dp &
                                           * zrho_rcp(:,jk) - pxtm1(:,jk,1))
!!$        pxtte(:,jk,1) = pxtte(:,jk,1) + &
!!$                        ztmst_rcp * (zcdnc(:,jk)* m_air/1000._dp &
!!$                                           * zrho_rcp(:,jk) - pxtm1(:,jk,1))
!qqq+-   
        !--- Update CDNC for radiation:
        pacdnc(:,jk)=zcdnc(:,jk)

        !--- Calculate new total tendency of ICNC:
!qqq+ see above
        pxtte(:,jk,2) = ztmst_rcp * (zicnc(:,jk)* m_air/1000._dp &
                                           * zrho_rcp(:,jk) - pxtm1(:,jk,2))
!!$        pxtte(:,jk,2) = pxtte(:,jk,2) + &
!!$                        ztmst_rcp * (zicnc(:,jk)* m_air/1000._dp &
!!$                                           * zrho_rcp(:,jk) - pxtm1(:,jk,2))
!qqq-
        !--- Diagnostics:
        qaut(1:kproma,jk,jrow)      = qaut(1:kproma,jk,jrow) - zdt*zrprn(:)
        qfre(1:kproma,jk,jrow)      = qfre(1:kproma,jk,jrow) - zdt*zfrln(:)
        qacc(1:kproma,jk,jrow)      = qacc(1:kproma,jk,jrow) - zdt*zsacln(:)
        cloud_tm1(1:kproma,jk,jrow) = paclc(1:kproma,jk)

        ll1_1d(:) = (zxlb(:)     >  zeps)    .AND. &
                    (zcdnc(:,jk) >= cdncmin)

        ztmp1_1d(:)         = cdnc_acc(1:kproma,jk,jrow) + zdtime*zcdnc(:,jk)
        cdnc_acc(1:kproma,jk,jrow) = MERGE(zcdnc(:,jk), cdnc_acc(1:kproma,jk,jrow), ll1_1d(:))

        !--- In-cloud CDNC burden:
        ztmp1_1d(:)     = zcdnc_burden(:)+zcdnc(:,jk)*zdz_2d(:,jk)
        zcdnc_burden(:) = MERGE(ztmp1_1d(:), zcdnc_burden(:), ll1_1d(:))
        CDNC_burden_acc(1:kproma,jrow) = zcdnc_burden(:)

        !---- CDNC and burden averaged over cloudy and cloud-free periods
        ztmp1_1d(:)     = cdnc(1:kproma,jk,jrow) + zdtime*zcdnc(:,jk)*zclcaux(:)
        cdnc(1:kproma,jk,jrow) = MERGE(zcdnc(:,jk)*zclcaux(:), cdnc(1:kproma,jk,jrow), ll1_1d(:))

        ztmp1_1d(:) = cdnc_burden(1:kproma,jrow) + zdtime*zcdnc(:,jk)*zdz_2d(:,jk)*zclcaux(:)
        cdnc_burden(1:kproma,jrow) = MERGE(cdnc_burden(1:kproma,jrow) + &
          zcdnc(:,jk)*zdz_2d(:,jk)*zclcaux(:), cdnc_burden(1:kproma,jrow), ll1_1d(:))

        !--- In-cloud effective radius [um]:
        ztmp1_1d(:)  = 0.00045e-6_dp*zcdnc(:,jk) + 1.18_dp !SF zkap
        zreffl_1d(:) = 1.E6_dp*ztmp1_1d(:)                                     &   
                     * ( (3._dp/(4._dp*api*rhoh2o)) * zxlb(:)                  &
                         * zrho(:,jk) / zcdnc(:,jk)           )**(1._dp/3._dp)   !SF zreffl

        zreffl_1d(:) = MAX(4._dp,MIN(zreffl_1d(:),40._dp))
        reffl(1:kproma,jk,jrow) =  MERGE(zreffl_1d(:), reffl(1:kproma,jk,jrow), ll1_1d(:))

        ll2_1d(:) = ll1_1d(:)                  .AND. &
                    (itop(:,jk)  == jk )       .AND. &
                    (ztp1tmp(:)   >  tmelt)    .AND. &
                    (zreffct(:)   <  4._dp)    .AND. &
                    (zreffl_1d(:) >= 4._dp)


        ll1_1d(:) = (zxib(:) > zeps)

        ztmp1_1d(:)         = icnc_acc(1:kproma,jk,jrow) + zdtime*zicnc(:,jk)
        icnc_acc(1:kproma,jk,jrow) = MERGE(zicnc(:,jk), icnc_acc(1:kproma,jk,jrow), ll1_1d(:))

        !--- Cirrus radii [um]:
        ll2_1d(:) = (ztp1tmp(:) > cthomi)

        ztmp1_1d(:) = 0.5e4_dp*( 1000._dp/0.0376_dp*MAX(zxib(:),zeps)*zrho(:,jk)/zicnc(:,jk)         )**0.302_dp !zrieff plate
        ztmp2_1d(:) = 1.0e6_dp*( 3._dp/(4._dp*api*zrhoice)*MAX(zxib(:),zeps)*zrho(:,jk)/zicnc(:,jk)) **0.333_dp
        ztmp2_1d(:) = (1.61_dp*ztmp2_1d(:)**3 + 3.56e-4_dp*ztmp2_1d(:)**6)**0.333_dp   !zrieff 2nd option

        ztmp3_1d(:) = MERGE(ztmp1_1d(:), ztmp2_1d(:), ll2_1d(:))
        ztmp3_1d(:) = MAX(ztmp3_1d(:), ceffmin)
        ztmp3_1d(:) = MIN(ztmp3_1d(:), ceffmax)   !zrieff

        reffi(1:kproma,jk,jrow) = MERGE(ztmp3_1d(:), reffi(1:kproma,jk,jrow), ll1_1d(:))

        !SF tovs diagnostics:
        ll2_1d(:) = ll1_1d(:)       .AND. &
                    .NOT. ll2_1d(:)

        ztmp1_1d(:) = 1000._dp*zxib(:)*zclcaux(:)*zdpg_2d(:,jk)  ! iwp in g/m2
        ztmp2_1d(:) = ztau1i(:) + 1.9787_dp*ztmp1_1d(:)*ztmp3_1d(:)**(-1.0365_dp)
        ztau1i(:)   = MERGE(ztmp2_1d(:), ztau1i(:), ll2_1d(:))
        
        ll3_1d(:) = ll2_1d(:)            .AND. &
                    (ztau1i(:) > 0.7_dp) .AND. &
                    (ztau1i(:) < 3.8_dp)

       !SF end tovs diagnostics 

        ll2_1d(:) = ll1_1d(:)                 .AND. &
                    (zicnc(:,jk) >= zicemin)

        !--- In-cloud ICNC burden:
        ztmp1_1d(:)     = zicnc_burden(:)+zicnc(:,jk)*zdz_2d(:,jk)
        zicnc_burden(:) = MERGE(ztmp1_1d(:), zicnc_burden(:), ll2_1d(:)) 
        ! mz_ht_20120404+
        ICNC_burden_acc(1:kproma,jrow) = zicnc_burden(:)
        ! mz_ht_20120404-


        !---- ICNC and burden averaged over cloudy and cloud-free periods
        ztmp1_1d(:)     = icnc(1:kproma,jk,jrow) + zdtime*zicnc(:,jk)*zclcaux(:)
        icnc(1:kproma,jk,jrow) = MERGE(zicnc(:,jk)*zclcaux(:), icnc(1:kproma,jk,jrow), ll2_1d(:))

        ztmp1_1d(:)         = icnc_burden(1:kproma,jrow) + zdtime*zicnc(:,jk)*zdz_2d(:,jk)*zclcaux(:)
        icnc_burden(1:kproma,jrow) =  MERGE(icnc_burden(1:kproma,jrow) + &
             zicnc(:,jk)*zdz_2d(:,jk)*zclcaux(:), icnc_burden(1:kproma,jrow), ll2_1d(:))  

!--- End included for CDNC/IC scheme -----------------------------------

!       8.4   Corrections: Avoid negative cloud water/ice
!
        ztmp1_1d(:)   = zxlp1_1d(:)      !SF zxlold

        ll1_1d(:)      = (zxlp1_1d(:) < ccwmin)

        zdxlcor_1d(:) = -ztmst_rcp * zxlp1_1d(:)
        zdxlcor_1d(:) = MERGE(zdxlcor_1d(:), 0._dp, ll1_1d(:))
        pxlte(:,jk)   = pxlte(:,jk) + zdxlcor_1d(:)

        ztmp1_1d(:)          = pxtte(:,jk,1) - ztmst_rcp*m_air*zcdnc(:,jk)*zrho_rcp(:,jk)*1.e-3_dp
        pxtte(:,jk,1) = MERGE(ztmp1_1d(:), pxtte(:,jk,1), ll1_1d(:))

        ll2_1d(:)     = (zxip1_1d(:) < ccwmin)

        zdxicor_1d(:) = -ztmst_rcp * zxip1_1d(:) 
        zdxicor_1d(:) = MERGE(zdxicor_1d(:), 0._dp, ll2_1d(:)) 
        pxite(:,jk)   = pxite(:,jk) + zdxicor_1d(:)

        ztmp1_1d(:)          = pxtte(:,jk,2) - ztmst_rcp*m_air*zicnc(:,jk)*zrho_rcp(:,jk)*1.e-3_dp
        pxtte(:,jk,2) = MERGE(ztmp1_1d(:), pxtte(:,jk,2), ll2_1d(:))

        paclc(:,jk)   = MERGE(0.0_dp,paclc(:,jk),ll1_1d(:) .AND. ll2_1d(:))
        paclcac(:,jk) = paclcac(:,jk) + zdtime*paclc(:,jk)

        pqte(:,jk)    = pqte(:,jk) - zdxlcor_1d(:) - zdxicor_1d(:)
        ptte(:,jk)    = ptte(:,jk) + zlvdcp(:)*zdxlcor_1d(:)                &
                                   + zlsdcp(:)*zdxicor_1d(:)

!SF for scavenging: create a 2d array of clcpre:
        zclcpre_2d(:,jk) = zclcpre(:)

831 END DO    ! Vertical loop
!

!     ------------------------------------------------------------------
!
!       10.    Diagnostics
!
!       10.1   Accumulated precipitation at the surface
!
   prsfl(:) = zrfl(:)
   pssfl(:) = zsfl(:)
   paprl(:) = paprl(:) + zdtime*(prsfl(:)+pssfl(:))
   paprs(:) = paprs(:) + zdtime*pssfl(:)

!
!       10.2   Total cloud cover
!
    zclcov(:) = 1.0_dp-paclc(:,1)

    DO 923 jk = 2,klev
       ztmp1_1d(:) = MAX(paclc(:,jk), paclc(:,jk-1))
       ztmp2_1d(:) = MIN(paclc(:,jk-1), zxsec)

       zclcov(:) = zclcov(:)*(1._dp - ztmp1_1d(:))  &
                            /(1._dp - ztmp2_1d(:))
923 END DO

    zclcov(:)  = 1.0_dp-zclcov(:)
    paclcov(:) = paclcov(:) + zdtime*zclcov(:)
!
!       10.3   Vertical integrals of humidity, cloud water and cloud ice
!
    zqvi(:)  = 0.0_dp
    zxlvi(:) = 0.0_dp
    zxivi(:) = 0.0_dp
!
    DO 933 jk = ktdia,klev
       zqvi(:)  = zqvi(:)  + pqm1(:,jk) *zdpg_2d(:,jk)
       zxlvi(:) = zxlvi(:) + pxlm1(:,jk)*zdpg_2d(:,jk) 
       zxivi(:) = zxivi(:) + pxim1(:,jk)*zdpg_2d(:,jk)
933 END DO
!
    pqvi(:)  = pqvi(:)  + zdtime*zqvi(:)
    pxlvi(:) = pxlvi(:) + zdtime*zxlvi(:)
    pxivi(:) = pxivi(:) + zdtime*zxivi(:)

!
  RETURN
END SUBROUTINE cloud_cdnc_icnc2

!===============================================================================
! main microphysics cloud rountine

  SUBROUTINE cloud_cdnc_icnc_cc(kproma, kbdim, ktdia, klev, klevp1, ztmst &
!---Included for in-cloud scavenging (Philip Stier, 28/03/01):----------
                         , ktrac,  krow                                &
!---End Included for scavenging-----------------------------------------
!-----------------------------------------------------------------------
! - INPUT  2D .
                         , pum1,     pvm1                              &
                         , paphm1,   pvervel                           &
                         , papm1,    papp1,    pacdnc                  &
                         , pqm1,     ptm1,     ptvm1                   &
                         , pxlm1,    pxim1,    pxtec                   &
                         , pxvar,    pxskew,   pqtec                   &
                         , pbetaa,   pbetab                            &
                         , pvdiffp,  phmixtau, pvmixtau                &
                         , pgeo,     pbetass                           &
!---Included for in-cloud scavenging (Philip Stier, 28/03/01):----------
                         , pxtm1                                       &
!---End Included for scavenging-----------------------------------------
!--- Included for prognostic CDNC/IC scheme ----------------------------
                         , ptkem1,   pcvcbot,  pwcape                  &
                         , pxtecl,   pxteci                            &
                         , pxtecnl,  pxtecni                           &
!--- End included for CDNC/IC scheme -----------------------------------
! - INPUT  1D .
                         , knvb                                        &
! - OUTPUT 2D .
                         , paclc,    paclcac,  prelhum                 &
! - INPUT/OUTPUT 1D .
                         , paclcov,  paprl,    pqvi                    &
                         , pxlvi,    pxivi                             &
! - OUTPUT 1D .
                         , pssfl,    prsfl                             &
! - INPUT/OUTPUT 2D .
                         , pqte,     ptte                              &
                         , pxlte,    pxite                             &
!---Included for in-cloud scavenging (Philip Stier, 28/03/01):----------
                         , pxtte                                       &
!---End Included for scavenging-----------------------------------------
! - INPUT/OUTPUT 1D .
                         , paprs                                       &
                         ! mz_ht_20120209+
                         , status_string,      lcover                  &
                         , slm,      glac,     pcdncact                &
                         , plwc,     piwc                              &
                         , pfrain,   pfsnow                            &
                         , pfrain_no,pfsnow_no                         &
                         , prate_r,  prate_s                           &
                         , prevap,   pssubl                            &
                         , pr_cover, pcond                             &
                         , pimelt,   pisedi                            &
                         ) 
!
!     *Cloud* computes large-scale water phase changes, precipitation,
!             cloud cover, and vertical integrals of specific humidity,
!             cloud liquid water content and cloud ice (diagnostics).
!
!     Subject.
!     --------
!
!          This rotuine computes the tendencies of the four prognostic
!          variables (temperature t, specific humidity q, cloud liquid
!          water xl, cloud ice xi) due to phase changes (condensation/
!          deposition, evaporation/sublimation of rain/snow falling
!          into the unsaturated part of the grid box, melting of snow,
!          melting/freezing of cloud ice/cloud water, sedimentation of
!          cloud ice, and precipitation formation in warm, cold and
!          mixed phase clouds.
!          The precipitation at the surface (rain and snow) is used in
!          later for computing the land surface hydrology in *surf*.
!          The cloud parameters (cloud cover, cloud liquid water and
!          cloud ice are used for the calculation of radiation at the
!          next timestep.
!          Attention: 
!          In the current version the advective tendencies of skewness 
!          and variance are set to zero.
!
!     INTERFACE.
!     ----------
!
!     *Call cloud*
!
!     Input arguments.
!     ----- ----------
!  - 2D
!  paphm1   : pressure at half levels                              (n-1)
!  papm1    : pressure at full levels                              (n-1)
!  papp1    : pressure at full levels                              (n-1)
!  pacdnc   : cloud droplet number concentration (specified)
!  pqm1     : specific humidity                                    (n-1)
!  ptm1     : temperature                                          (n-1)
!  pxlm1    : cloud liquid water                                   (n-1)
!  pxim1    : cloud ice                                            (n-1)
!  pxtec    : detrained convective cloud liquid water or cloud ice (n)
!  pxvar    : distribution width (b-a)                             (n-1)
!  pxskew   : beta shape parameter "q"                             (n-1)
!  pbetaa   : the beta distribution minimum a                      (n-1)
!  pbetab   : the beta distribution maximum b                      (n-1)
!  pvdiffp  : the rate of change of q due to vdiff scheme          (n-1)
!  phmixtau : mixing timescale**-1 for horizontal turbulence       (n)
!  pvmixtau : mixing timescale**-1 for horizontal turbulence       (n)
!--- Included for prognostic CDNC/IC scheme ----------------------------
!  zcdnc    : in-cloud value, cloud droplet number concentration [m-3] 
!             = pxtm1(1,1,idt_cdnc)*density(air)
!  zicnc    : in-cloud value, ice crystal number concentration   [m-3] 
!             = pxtm1(1,1,idt_icnc)*density(air)
!--- End included for CDNC/IC scheme -----------------------------------
!
!  - 1D
!  knvb     :
!
!     Output arguments.
!     ------ ----------
!  - 1D
!  prsfl    : surface rain flux
!  pssfl    : surface snow flux
!
!     Input, Output arguments.
!     ------------ ----------
!  - 2D
!  paclc    : cloud cover  (now diagnosed in cover)
!  paclcac  : cloud cover, accumulated
!  paclcov  : total cloud cover
!  paprl    : total stratiform precipitation (rain+snow), accumulated
!  pqvi     : vertically integrated spec. humidity, accumulated
!  pxlvi    : vertically integrated cloud liquid water, accumulated
!  pxivi    : vertically integrated cloud ice, accumulated
!  ptte     : tendency of temperature
!  pqte     : tendency of specific humidity
!  pxlte    : tendency of cloud liquid water
!  pxite    : tendency of cloud ice
!--- Included for prognostic CDNC/IC scheme ----------------------------
!  pxtte    : tendency of cloud droplet number, in-cloud 
!  pxtte    : tendency of ice crystal number,in-cloud
!--- End included for CDNC/IC scheme -----------------------------------
!  - 1D
!  paprs    : Snowfall, accumulated
!
!     Externals.
!     ----------
!
!     Method.
!     -------
!     see References
!
!     References.
!     ----------
!
!     Lohmann and Roeckner, 1996: Clim. Dyn. 557-572
!     Levkov et al., 1992: Beitr. Phys. Atm. 35-58.          (ice phase)
!     Beheng, 1994: Atmos. Res. 193-206.                    (warm phase)
!     Lenderink et al., 1998; KNMI-REPORT NO. 98-13       (condensation)
!     Tompkins 2002, J. Atmos. Sci.                        (cloud cover)
!
!     Authors.
!     -------
!     M.Esch        MPI-Hamburg  1999
!     G.Lenderink   KNMI, de Bilt 1998
!     U.Lohmann     MPI-Hamburg  1995
!
!     Modifications.
!     --------------
!     E.Roeckner    MPI-Hamburg  2000
!     A.Tompkins    MPI-Hamburg  2000
!     U.Schlese     MPI-Hamburg  2003
!     U.Lohmann     Dalhousie University 2002-2006: Prognostic CDNC/IC scheme
!     P.Stier       MPI-Hamburg          2002-2006: Prognostic CDNC/IC scheme
!                                                   Scavenging parameters
!     J.Zhang       Dalhousie University      2004: Prognostic CDNC/IC scheme
!     S. Ferrachat  ETH Zuerich  2008: complete re-writing to allow vectorization 
!                                      cleanup of the code 
!
    USE messy_cloud_ori,  ONLY : cqtmin, tmelt, cvtfall, crhosno, cn0s    &
                               , cthomi, csecfrl, ncctop, cvarmin         &
                               , cbeta_pq, cbeta_pq_max, nbetaq, cbetaqs  &
                               , rbetak, nbetax, tbetai0, tbetai1         &
                               , clmax, clmin, jbmin, jbmax, lonacc       &
                               , ccraut, crhoi, ccsaut                    &
                               , cbeta_cs, LOOKUPOVERFLOW
    USE messy_main_tools, ONLY : jptlucu1, jptlucu2,                      &
                                 tlucua, tlucuaw, tlucub

    USE messy_main_constants_mem, ONLY: ceffmin, ceffmax, ccwmin, &
                                        rd, rv, vtmpc1, vtmpc2,   &
                                        cpd => cp_air,  &
                                        tmelt, rhoh2o => rho_h2o

  IMPLICIT NONE
 
  INTEGER :: krow, ktdia, kproma, kbdim, klev, klevp1, ktrac 

  REAL(dp):: paphm1(kbdim,klevp1), pvervel(kbdim,klev)  &
            ,papm1(kbdim,klev)   , pqm1(kbdim,klev)     &
            ,papp1(kbdim,klev)   , ptm1(kbdim,klev)     &
            ,ptvm1(kbdim,klev)   , pxlm1(kbdim,klev)    &
            ,pxim1(kbdim,klev)   , pxtec(kbdim,klev)    &
            ,pqtec(kbdim,klev)   , pxvar(kbdim,klev)    &
            ,pxskew(kbdim,klev)  , pbetaa(kbdim,klev)   &
            ,pbetab(kbdim,klev)  , pvdiffp(kbdim,klev)  &
            ,phmixtau(kbdim,klev), pvmixtau(kbdim,klev) &
            ,pgeo(kbdim,klev)    , pbetass(kbdim,klev)  &
            ,pxlvi(kbdim)        ,pxivi(kbdim)          & 
            ,paclc(kbdim,klev)   ,paclcac(kbdim,klev)   &
            ,pacdnc(kbdim,klev)  ,prelhum(kbdim,klev)   &
            ,paclcov(kbdim)      ,paprl(kbdim)          &
            ,pqvi(kbdim)         ,pssfl(kbdim)          &
            ,ptte(kbdim,klev)    ,pqte(kbdim,klev)      & 
            ,pxlte(kbdim,klev)   ,pxite(kbdim,klev)     &
            ,paprs(kbdim)        ,prsfl(kbdim)

!---Included for in-cloud scavenging (Philip Stier, 28/03/01):----------
  REAL(dp) ::    pxtm1(kbdim,klev,ktrac), pxtte(kbdim,klev,ktrac) 
  REAL(dp) ::    pum1(kbdim,klev), pvm1(kbdim,klev)
!---End Included for scavenging-----------------------------------------
  CHARACTER(LEN=32) :: status_string
  LOGICAL  :: lcover
  REAL(dp) :: slm(kbdim), glac(kbdim)

  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev)         :: plwc,     piwc
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev), TARGET :: pfrain,   pfsnow
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev)         :: pfrain_no,pfsnow_no
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev), TARGET :: prevap,   pssubl
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev), TARGET :: prate_r,  prate_s 
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev)         :: pr_cover
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev)         :: pimelt,   pisedi
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev)         :: pcond 

  ! Local variables
  REAL(dp), POINTER, DIMENSION(:) :: zrfl, zsfl, zsub, zevp, zrpr

!
!   Temporary arrays
!

  REAL(dp)   :: zclcpre_2d(kbdim,klev)    ,zclcpre(kbdim)              &
           ,zcnd(kbdim)         ,zdep(kbdim)                           &
           ,                     zxievap(kbdim)     ,zxlevap(kbdim)    &
           ,zfrl(kbdim,klev)    ,zimlt(kbdim)       ,zsmlt(kbdim)      &
           ,                     zspr(kbdim)                           &
           ,zxlte(kbdim)        ,zxite(kbdim)       ,zxiflux(kbdim)    &
           ,zsacl(kbdim)        ,zxlte2(kbdim)      ,zxite2(kbdim)     &
           ,zlsdcp(kbdim)       ,zlvdcp(kbdim)      ,zximlt(kbdim)     &
           ,ztp1tmp(kbdim)      ,zqp1tmp(kbdim)     ,zxisub(kbdim)     &
           ,zxlb(kbdim)         ,zxib(kbdim)        ,zxibold_1d(kbdim) &
           ,zrho(kbdim,klev)    ,zclcov(kbdim)      ,zclcaux(kbdim)    &
           ,zrho_rcp(kbdim,klev)                                       &
           ,zqvi(kbdim)         ,zxlvi(kbdim)       ,zxivi(kbdim)      &
           ,zbetaqt(kbdim)      ,zwide(kbdim)                          &
           ,zbetacl(kbdim)      ,zturbvar(kbdim)    ,zturbskew(kbdim)  &
           ,zconvvar(kbdim)     ,zconvskew(kbdim)   ,zvartg(kbdim)     &
           ,zmicroskew(kbdim)   ,zgenti(kbdim)      ,zgentl(kbdim)     &
           ,zxvarte(kbdim)      ,zxskewte(kbdim)                       &
           ,zgeoh(kbdim,klevp1) 
!
  INTEGER   knvb(kbdim)

!--- Included for dust emissions (Philip Stier 10/01/02)-----------------
  INTEGER :: jrow
!--- End Included for dust emissions ------------------------------------

!---Included for in-cloud scavenging (Philip Stier, 28/03/01):----------
  REAL(dp)    :: zmratepr(kbdim,klev), & ! Rain formation rate in cloudy part
                                     ! of the grid box [kg/kg]
             zmrateps(kbdim,klev), & ! Ice  formation rate in cloudy part
                                     ! of the grid box  [kg/kg]
             zfrain(kbdim,klev),   & ! Rain flux before evaporation
                                     ! [kg/m2/s]
             zfsnow(kbdim,klev),   & ! Snow flux before sublimation
                                     ! [kg/m2/s]
             zfevapr(kbdim,klev),  & ! Evaporation of rain [kg/m2/s]
             zfsubls(kbdim,klev),  & ! Sublimation of snow [kg/m2/s]
             zmlwc(kbdim,klev),    & ! In-cloud liquid water mass mixing
                                     ! ratio before rain formation [kg/kg]
             zmiwc(kbdim,klev),    & ! In-cloud ice mass mixing ratio
                                     ! before snow formation [kg/kg]
             zmsnowacl(kbdim,klev)   ! Accretion rate of snow with cloud
                                     ! droplets in cloudy part of the 
                                     ! grid box  [kg/kg]
!---End Included for scavenging-----------------------------------------

!--- Included for prognostic CDNC/IC scheme ----------------------------

  REAL(dp) :: pcdncact(kbdim,klev),     ptkem1(kbdim,klev),             &
             pcvcbot(kbdim),           pwcape(kbdim),                  &
             pxtecl(kbdim,klev),       pxteci(kbdim,klev),             &
             pxtecnl(kbdim,klev),       pxtecni(kbdim,klev)

  INTEGER ::         &
             itop(kbdim,klev), ibas(kbdim,klev), &
             icl_minusbas(kbdim,klev), icl_minustop(kbdim,klev),    & !SFtemp: iclminustop: useless for now
             iclbas(kbdim,klev), itm1_look(kbdim,klev), it1_1d(kbdim) 

  REAL(dp)    :: zrprn(kbdim),             zsprn(kbdim,klev),              &
             zsacln(kbdim),            zfrln(kbdim),                   &
             zcdnc(kbdim,klev),     & !
             zcdnc_burden(kbdim),   &
             zqlnuc(kbdim,klev),    & !
             zqlnuccvh(kbdim,klev), & ! Nucleated CDNC from convective detrainment
             zqlnuccv(kbdim,klev), & 
             zqlnuc_bas(kbdim,klev), &
             zcdnc_bas(kbdim,klev)

!  REAL(dp)   zw_strat(kbdim,klev)!,  & ! stratiform updraft velocity, large-scale+TKE (>0.0) [m s-1]
!             zw_conv(kbdim,klev)      ! convective updraft velocity, large-scale+CAPE (>0.0) [m s-1]

 ! Temporary fields
  INTEGER  :: jl, jk, jkk, jmod, &
!SFcray
              iqidx, ixidx
!SFcray

  REAL(dp) :: zesw,    &                  
              zf1, &        
              zrieff, &
              zradl,  &
              zdisp, zdw0, ztmst, ztmst_rcp, zcons1, zcons2, zcons3, zcons4, zcons5, zdtime, zdtime_rcp, ztmstdt, zdt, g_rcp,    &
              zpirho, zpirho_rcp, &
              zb2, &
              zdv,  &
              zgtp, & 
              zexm1_1, zexp_1, &
              zexm1_2, zexp_2, &
              zrih, &
              ztte, zomega, &
              ztc, zvth, zfuchs, zfre, zre, &
              zfracdusol, zfracduinsolai, zfracduinsolci, zfracbcsol, zfracbcinsol, &
              zfrzcntdu, zfrzcntbc, zfrzcnt, znaimmdu, znaimmbc, zfrzimm, zfrzthermo, &
              zqsm1, zqst1, zdqsdtime, zdeltatemp, zf2, zcap

  REAL(dp)    :: zicnc(kbdim,klev),     & !
             zicncq(kbdim,klev),    &
             zicnc_burden(kbdim),   & !
             zascs(kbdim,klev),     & !
             zrid_2d(kbdim,klev),   &
             zqsi_2d(kbdim,klev),   &
             zninucl(kbdim, klev),  & ! number conc. of newly nucleated IC
             zqinucl(kbdim, klev),  & ! mixing ratio of newly nucleated IC
             zri(kbdim, klev),      & ! size of newly nucleated IC
             zreffct(kbdim),        &  ! cloud top effective radius
             ztau1i(kbdim),         &  ! ice cloud optical depth - visible wavelength
             znidetr(kbdim, klev)     ! IC number from detrainment

 !SFnew arrays
   !SF 2D
 INTEGER    :: itm1p1_look(kbdim,klev), it_1d(kbdim)

 REAL(dp)   :: zdz_2d(kbdim,klev), zdp_2d(kbdim,klev), zdpg_2d(kbdim,klev), zaaa_2d(kbdim,klev), zviscos_2d(kbdim,klev), &
               zqswp1_2d(kbdim,klev), zqsip1_2d(kbdim,klev), zqsw_2d(kbdim,klev), &
               zastbstw(kbdim,klev), zastbsti(kbdim,klev), zmmean_2d(kbdim,klev), &
               zxifallmc_2d(kbdim,klev), zalfased_2d(kbdim,klev), zbetased_2d(kbdim,klev), &
               zxifallnc_2d(kbdim,klev), &
               zxised_2d(kbdim,klev),zicesed_2d(kbdim,klev), &
             zesw_2d(kbdim,klev),   & !
             zesi_2d(kbdim,klev),   & !
             zsusatw_2d(kbdim,klev),& ! supersat. with respect to water
             zsusatw_evap(kbdim,klev),& 
             zsusatix_2d(kbdim,klev),& 
             zvervx_2d(kbdim,klev),  &     ! updraft [cm/s]
             zicesub(kbdim,klev), zqrho_2d(kbdim,klev)
               
   !SF 1D
 REAL(dp)   :: zscnc_1d(kbdim), zdplanar_1d(kbdim), zusnow_1d(kbdim), zsacl2in_1d(kbdim), &
               zxsp2_1d(kbdim), zxsp_1d(kbdim), zstcrit_1d(kbdim),  &
               zstokes_1d(kbdim), zudrop_1d(kbdim), zrey_1d(kbdim), &
               zrac1_1d(kbdim), zrac2_1d(kbdim), zrautself_1d(kbdim), zxsp1_1d(kbdim), &
               zraut_1d(kbdim), zsaut_1d(kbdim), zsaci2_1d(kbdim), zsacl2_1d(kbdim), &
               zris_1d(kbdim), zcolleffi_1d(kbdim), zc1_1d(kbdim), zxlp1_1d(kbdim), &
               zreffl_1d(kbdim), zdxlcor_1d(kbdim), &
               zdxicor_1d(kbdim), zsecprod_1d(kbdim), &
               zcsacl_1d(kbdim), zpretot_1d(kbdim), zpredel_1d(kbdim), zpresum_1d(kbdim), &
               zzdrr_1d(kbdim), zzdrs_1d(kbdim), zxidtstar_1d(kbdim), &
               zlc_1d(kbdim), zxldtstar_1d(kbdim), zxlm1evp_1d(kbdim), zxim1evp_1d(kbdim), &
               zxldt_1d(kbdim), zxidt_1d(kbdim), zqsm1_1d(kbdim), zqp1_1d(kbdim), zdtdt_1d(kbdim), &
               zdqsat_1d(kbdim), ztp1_1d(kbdim), zdqsdt_1d(kbdim), zqst1_1d(kbdim), &
               zqvdt_1d(kbdim), &
               zxip1_1d(kbdim), zcorw_1d(kbdim), zesw_1d(kbdim), zoversatw_1d(kbdim), &
               zqsp1tmpw_1d(kbdim), zqsp1tmp_1d(kbdim), zcor_1d(kbdim), zrhtest_1d(kbdim), &
               zoversat_1d(kbdim), zlcdqsdt_1d(kbdim), zclcstar_1d(kbdim), zxrp1_1d(kbdim), &
               zauloc_1d(kbdim), zqsp1tmphet_1d(kbdim), zqcon_1d(kbdim), zrelhum_1d(kbdim), &
               iqidx_1d(kbdim), zbap1_1d(kbdim), zgent_1d(kbdim), ixidx_1d(kbdim), &
               zqtau_1d(kbdim), zxilb_1d(kbdim), zbbap1_1d(kbdim), zbqp1_1d(kbdim), &
               zqcdif_1d(kbdim), zlucuawp1_1d(kbdim), zlucuap1_1d(kbdim), zes_1d(kbdim), &
               zlucub_1d(kbdim), zlucuaw_1d(kbdim), zlucua_1d(kbdim), zicncp1_1d(kbdim), &
               zetaair_1d(kbdim), zkair_1d(kbdim), &
               zdfarbcki_1d(kbdim), zdfarduai_1d(kbdim), zdfarduci_1d(kbdim), &
               zknbcki_1d(kbdim), zknduai_1d(kbdim), zknduci_1d(kbdim), &
               zftbcki_1d(kbdim), zftduai_1d(kbdim), zftduci_1d(kbdim), &
               zvervmax_1d(kbdim), zrice_1d(kbdim), zeta_1d(kbdim), zdv_1d(kbdim)
               
  REAL(dp)  :: ztmp(kbdim,klev), ztmp1(kbdim,klev), ztmp2(kbdim,klev), &
               ztmp1_1d(kbdim), ztmp2_1d(kbdim), ztmp3_1d(kbdim), &
               ztmp4_1d(kbdim), ztmp5_1d(kbdim)

  LOGICAL :: ll_look(kbdim,klev), &
             ll_cv(kbdim,klev), &
             ll_ice(kbdim,klev), ll1(kbdim,klev), &
             ll1_2d(kbdim,klev), ll2_2d(kbdim,klev), ll3_2d(kbdim,klev), lo_1d(kbdim), &
             lo2_1d(kbdim), locc_1d(kbdim), &
             ll1_1d(kbdim), ll2_1d(kbdim), ll3_1d(kbdim), ll4_1d(kbdim), & 
             ll5_1d(kbdim), ll6_1d(kbdim), ll7_1d(kbdim), ll8_1d(kbdim) 
 !SF end new arrays

  !--- Included for cirrus scheme -------------------------------------------
  REAL(dp) :: znicex(kbdim, klev)     ! 
  REAL(dp) :: zapnx(kbdim,nfrzmod) ! aerosol number available for homogeneous freezing [1/cm3]
  REAL(dp) :: zaprx(kbdim,nfrzmod) ! radius of aerosols available for homogeneous freezing  [cm]
  REAL(dp) :: zapsigx(kbdim,nfrzmod)! std. dev. of aerosols available for homogeneous freezing 
  REAL(dp) :: zap(kbdim,klev)

  LOGICAL, PARAMETER :: nosize = .false. ! .true. ---> aerosol size effects
                                         ! are ignored for homogeneous
                                         ! freezing

    ! !!! lhet = .true. NOT YET IMPLEMENTED !!!

  LOGICAL, PARAMETER :: lhet = .false.   ! .true. ---> heterogeneous
                                         ! freezing of aerosol particles
                                         ! below 235K is considered
  !--- End included for cirrus scheme -------------------------------------------

  REAL(dp) :: zxifluxn(kbdim)

!--- End Included for CDNC -------------------------------------------
  REAL(dp) :: zrwetki_2d(kbdim,klev), zrwetai_2d(kbdim,klev), zrwetci_2d(kbdim,klev)
!--- Ulrike: for additional disgnostics ---------------------------------------
  REAL(dp) :: zrdry(kbdim,klev,nmod) ! dry radius for each mode
                                    !@@@ currently re-stored, check if avoidable
!--- end Ulrike: for additional diagnostics -----------------------------------

  REAL(dp) :: coolr, scr, pice, pw, znicex0, zri0, vice
  REAL(dp) :: zri1(kbdim, klev), zri2(kbdim, klev), zri3(kbdim, klev), &
              znicex1(kbdim, klev), znicex2(kbdim, klev), znicex3(kbdim, klev), &
              zvervxmod_2d(kbdim, klev), znicex_kor1(kbdim, klev), znicex_kor2(kbdim, klev), &
              zmult_1d(kbdim), zqinucl_new(kbdim,klev), zicncqi(kbdim,klev)
  !CCMod
  LOGICAL  :: ll_bcc(kbdim)
  REAL(dp) :: zumf, zwidthini, zheightini, zztune2, zmi_cc
  REAL(dp) :: zmult3(kbdim), zmult7(kbdim),                &
              zxidt_cc(kbdim), zxidtsed_cc(kbdim), zcorc2(kbdim), zcorb(kbdim),          &
              zxiflux_cc(kbdim), zxifluxn_cc(kbdim), zximlt_cc(kbdim), zimlt_cc(kbdim),  &
              zxievap_cc(kbdim), zqccdif_1d(kbdim), zqcodif_1d(kbdim), zxip1_cc(kbdim),  &
              zdepcc(kbdim), zdepco(kbdim), zxibco(kbdim), zxibcc(kbdim), zcorclc(kbdim),& 
              zcolleffi_cc(kbdim), zc1_cc(kbdim), zsaut_cc(kbdim), zxibold_cc(kbdim),    &
              zspr_cc(kbdim), zxisub_cc(kbdim), zicncp1_cc(kbdim), zccxite(kbdim)
  REAL(dp) :: zfkme(kbdim, klev), zfh2o(kbdim, klev),                                    &
              zflaegb(kbdim, klev), zcctau(kbdim, klev),                                 &
              zcccov_new(kbdim, klev), zccvol_new(kbdim, klev), zcclen_new(kbdim, klev), &
              zccicnc_new(kbdim, klev), zcciwc_new(kbdim, klev),                         &
              zcccov_ges(kbdim, klev), zccvol_ges(kbdim, klev),                          &
              zcclen_ges(kbdim, klev), zcclen_ges2(kbdim, klev),                         &
              zcccov_ori(kbdim, klev), zcccov_all(kbdim, klev), zcccov_tau(kbdim, klev), &
              zccvol_ges_ori(kbdim, klev), zccvol_new_ori(kbdim, klev),                  &
              zccvol_ges2(kbdim, klev), zccvol_all(kbdim, klev),                         &
              zcccov_1dt(kbdim, klev), zcccov_2dt(kbdim, klev), zcccov_3dt(kbdim, klev), &
              zcccov_4dt(kbdim, klev), zcccov_5dt(kbdim, klev), zcccov(kbdim, klev),     &
              zccvol_1dt(kbdim, klev), zccvol_2dt(kbdim, klev), zccvol_3dt(kbdim, klev), &
              zccvol_4dt(kbdim, klev), zccvol_5dt(kbdim, klev), zccvol(kbdim, klev),     &
              zcclen_1dt(kbdim, klev), zcclen_2dt(kbdim, klev), zcclen_3dt(kbdim, klev), &
              zcclen_4dt(kbdim, klev), zcclen_5dt(kbdim, klev), zcclen(kbdim, klev),     &
              zcclen_ou(kbdim, klev),                                                    &
              zccicnc(kbdim, klev), zcciwc(kbdim, klev),                                 &
              zconcovinitial(kbdim, klev), zdepco_form(kbdim, klev),                     &
              zcovini(kbdim, klev), zcovini_k(kbdim, klev), zcorcov(kbdim, klev),        &
              zxised_cc(kbdim, klev), zicesed_cc(kbdim, klev),                           &
              zxisedf_cc(kbdim, klev), zicesedf_cc(kbdim, klev),                         &
              zb_stern(kbdim, klev), za_stern(kbdim, klev), za_stern2(kbdim, klev),      & 
              zfre_2d(kbdim, klev), zcap_2d(kbdim, klev), zsprn_cc(kbdim, klev),         &
              zhelp(kbdim, klev), zcccor(kbdim, klev),                                   &
              zcccor_vol(kbdim, klev), zxim1evp_cc(kbdim, klev),                         &
              zlmit(kbdim, klev), zmult2(kbdim, klev), zdeltu(kbdim, klev),              &
              zmult(kbdim, klev), zcor44(kbdim, klev), zfrac1(kbdim, klev)
!
! Executable statements
!
  jrow = krow

!-- 0. Initialisations:

  lookupoverflow = .FALSE.

  ! initializing new values for channel objects
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
  reffi(:,:,jrow)         = ceffmax
  reffl(:,:,jrow)         = 40._dp
  cdnc(:,:,jrow)          = 0._dp
  icnc(:,:,jrow)          = 0._dp

!---Included for in-cloud scavenging (Philip Stier, 28/03/01):----------
  zmratepr(:,:)  = 0._dp
  zmrateps(:,:)  = 0._dp
  zfrain(:,:)    = 0._dp
  zfsnow(:,:)    = 0._dp
  zfevapr(:,:)   = 0._dp
  zfsubls(:,:)   = 0._dp
  zmlwc(:,:)     = 0._dp
  zmiwc(:,:)     = 0._dp
  zmsnowacl(:,:) = 0._dp
!---End Included for scavenging-----------------------------------------

!--- Included for prognostic CDNC/IC scheme ----------------------------
  zqinucl(:,:) = 0._dp
  zninucl(:,:) = 0._dp
  zreffct(:)   = 0._dp
  ztau1i(:)    = 0._dp
!--- End Included for CDNC --------------------------------------------

!--- Included for prognostic CDNC/IC scheme ----------------------------

!cyjstart----------------------------
  zxifluxn(:)       = 0.0_dp
  zmmean_2d(:,:)    = zmi
  zxifallmc_2d(:,:) = 0.0_dp
  zxifallnc_2d(:,:) = 0.0_dp
  zxised_2d(:,:)    = 0.0_dp
  zicesed_2d(:,:)   = 0.0_dp
!cyjend----------------------------------------

  zcdnc_burden(1:kproma) = 0.0_dp
  zcdnc_bas(1:kproma,:)  = 0.0_dp 
  zqlnuc(1:kproma,:)     = 0.0_dp
  zqlnuc_bas(1:kproma,:) = 0.0_dp
  zicnc_burden(1:kproma) = 0.0_dp
  zri(1:kproma,:)        = 1.e-6_dp

!--- End included for CDNC/IC scheme -----------------------------------

!---Computational constants:
  zdisp      = EXP(0.5_dp * zsigmaw**2)
  zdw0       = 10.e-6_dp*zdisp       ! modal diameter times dispersion parameter

  zdtime     = ztmst / 2._dp
  zdtime_rcp = 1._dp / zdtime
  ztmst_rcp  = 1._dp / ztmst
  ztmstdt    = ztmst * zdtime
  zdt        = ztmst_rcp * zdtime 
  g_rcp      = 1._dp / g
  zcons1     = cpd*vtmpc2
  zcons2     = ztmst_rcp * g_rcp
  zexm1_1    = 2.47_dp-1.0_dp
  zexp_1     = -1._dp / zexm1_1
  zexm1_2    = 4.7_dp-1.0_dp
  zexp_2     = -1._dp / zexm1_2
  zpirho     = api*rhoh2o
  zpirho_rcp = 1._dp / zpirho
  zcap       = 2._dp / api

  zcons3 = 1._dp / ( api*crhosno*cn0s*cvtfall**(1._dp/1.16_dp) )**0.25_dp
  zcons4 = 1._dp / ( api*crhosno*cn0s )**0.8125_dp
  zcons5 = 1._dp / ( api*crhosno*cn0s )**0.875_dp

  !Initialisierung fuer Zirren

  zvervxmod_2d(:,:) = 0._dp
  zri1(:,:) = 0._dp
  znicex1(:,:) = 0._dp
  zri2(:,:) = 0._dp
  znicex2(:,:) = 0._dp
  zri3(:,:) = 0._dp
  znicex3(:,:) = 0._dp
  znicex_kor1(:,:) = 0._dp
  znicex_kor2(:,:) = 0._dp
  zqinucl_new(:,:) = 0._dp

!--- CCMod -------------------------------------------------------------
!Initialisierung fuer Kondensstreifenzirren

  zxised_cc(:,:) = 0._dp
  zxisedf_cc(:,:) = 0._dp
  zicesed_cc(:,:) = 0._dp
  zicesedf_cc(:,:) = 0._dp
  zxiflux_cc(:) = 0.0_dp
  zxifluxn_cc(:) = 0.0_dp
  zxim1evp_cc(:,:) = 0._dp
  za_stern(:,:) = 0._dp
  za_stern2(:,:) = 0._dp
  zb_stern(:,:) = 0._dp
  zccvol_ges2(:,:) = 0._dp

!SF store wet radius for freezing calculations (6.):

  zrwetki_2d(1:kproma,:) = mode(iaiti)%wetrad(1:kproma,:)
  zrwetai_2d(1:kproma,:) = mode(iacci)%wetrad(1:kproma,:)
  zrwetci_2d(1:kproma,:) = mode(icoai)%wetrad(1:kproma,:)

! Ulrike: for thermophoresis
  DO jmod = 1,nmod
    zrdry(1:kproma,:,jmod) = mode(jmod)%dryrad(1:kproma,:)
  ENDDO
! end Ulrike: for thermophoresis

!--- Get several utility variables:
  CALL get_util_var(kproma, kbdim, ktdia, klev, klevp1,           &
                    paphm1, pgeo, papm1, ptm1, zgeoh,             &
                    zdp_2d, zdpg_2d, zdz_2d, zaaa_2d, zviscos_2d)


!
!     ------------------------------------------------------------------
!
!       1.   Top boundary conditions, air density
!
!       1.1   Set to zero precipitation fluxes etc.
!
  zclcpre(:) = 0.0_dp   ! fraction of grid box covered by precip
  zxiflux(:)      = 0.0_dp   ! flux of ice crystals falling into the grid box from above
!
!       1.2   Air density


  DO 122 jk = ktdia,klev

     zrho(:,jk)     = papm1(:,jk)/(rd*ptvm1(:,jk))
     zrho_rcp(:,jk) = 1._dp / zrho(:,jk)
     zqrho_2d(:,jk) = 1.3_dp * zrho_rcp(:,jk)
     zfrl(:,jk)     = 0.0_dp

!-----------------Added by Junhua Zhang for CONV Micro--------------------
     IF (ncvmicro == 0) THEN

        pxtec(:,jk)  = MAX(pxtec(:,jk),0.0_dp)

        ll1_1d(:) = (pxtec(:,jk) > 0.0_dp)

        ztmp1_1d(:)         = twc_conv(1:kproma,jk,jrow) + ztmstdt*pxtec(:,jk)
        twc_conv(1:kproma,jk,jrow) = MERGE(ztmp1_1d(:), twc_conv(1:kproma,jk,jrow), ll1_1d(:))

        ztmp1_1d(:)          = conv_time(1:kproma,jk,jrow) + zdtime
        conv_time(1:kproma,jk,jrow) = MERGE(ztmp1_1d(:), conv_time(1:kproma,jk,jrow), ll1_1d(:))
 
     ELSE !SF ncvmicro .ne. 0

           pxtecl(:,jk)  = MAX(pxtecl(:,jk),  0.0_dp)
           pxteci(:,jk)  = MAX(pxteci(:,jk),  0.0_dp)
           pxtecnl(:,jk) = MAX(pxtecnl(:,jk), 0.0_dp)
           pxtecni(:,jk) = MAX(pxtecni(:,jk), 0.0_dp)

           ll1_1d(:) = (pxtecl(:,jk) > 0._dp) .AND. &
                       (ptm1(:,jk)   < cthomi)

           ztmp1_1d(:)  = pxteci(:,jk) + pxtecl(:,jk)
           pxteci(:,jk) = MERGE(ztmp1_1d(:), pxteci(:,jk), ll1_1d(:))

           pxtecl(:,jk) = MERGE(0._dp, pxtecl(:,jk), ll1_1d(:))

           ztmp1_1d(:)   = pxtecni(:,jk) + pxtecnl(:,jk)
           pxtecni(:,jk) = MERGE(ztmp1_1d(:), pxtecni(:,jk), ll1_1d(:))

           pxtecnl(:,jk) = MERGE(0._dp, pxtecnl(:,jk), ll1_1d(:))

           zfrl(:,jk) = MERGE(pxtecl(:,jk), zfrl(:,jk), ll1_1d(:))

           ll1_1d(:) = (pxteci(:,jk) > 0.0_dp) .AND. &
                       (ptm1(:,jk)   > tmelt)

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

     ENDIF  !SF end ncvmicro

!------------------------------End----------------------------------------

     pqtec(:,jk)  = MAX(pqtec(:,jk),0.0_dp)

!--- Included for prognostic CDNC/IC scheme ----------------------------

     ! calculate cloud droplet number concentration
!qqq+ include tendency for correct start value (see Lohman07)?
     zcdnc(:,jk) = zrho(:,jk)*(pxtm1(:,jk,1)/m_air * 1000._dp)
     zcdnc(:,jk) = zrho(:,jk)*((pxtm1(:,jk,1)+pxtte(:,jk,1)*ztmst) &
          /m_air * 1000._dp)
!qqq-
     zcdnc(:,jk) = MAX(zcdnc(:,jk),cqtmin)

     ! calculate ice crystal number
!qqq+ see above
!!$     zicnc(:,jk)  = zrho(:,jk)*(pxtm1(:,jk,2)/m_air * 1000._dp)
     zicnc(:,jk)  = zrho(:,jk)*((pxtm1(:,jk,2)+pxtte(:,jk,2)*ztmst) &
          /m_air * 1000._dp)
!qqq-
     zicnc(:,jk)  = MAX(zicnc(:,jk),cqtmin)
     zicncq(:,jk) = zicnc(:,jk)

     ! END changed S.K.
     ! end corinna
122 END DO
!
!       1.3   Geopotential at half levels  !SF moved in mo_cloud_utils
!
!SF end moved section

  CALL get_cloud_bounds(kproma, kbdim, ktdia, klev, paclc, itop, ibas, icl_minustop, icl_minusbas)

  ztmp(:,:)      = 1000._dp*ptm1(:,:)
  itm1_look(:,:) = NINT(ztmp(:,:))

  ll_look(:,:) = (itm1_look(:,:)<jptlucu1 .OR. itm1_look(:,:)>jptlucu2)
  IF (ANY(ll_look(1:kproma,ktdia:klev))) lookupoverflow = .TRUE.


  itm1_look(:,:) = MAX(MIN(itm1_look(:,:),jptlucu2),jptlucu1)

!SF for later use in 5.:
  itm1p1_look(:,:) = itm1_look(:,:) + 1
  itm1p1_look(:,:) = MAX(MIN(itm1p1_look(:,:),jptlucu2),jptlucu1)
!SFend for later use in 5.

  DO jk=klev,ktdia,-1
     DO jl=1,kproma
!SF water:
        zesw              = tlucuaw(itm1_look(jl,jk))/papm1(jl,jk)
        zesw              = MIN(zesw,0.5_dp)
        zqsw_2d(jl,jk)    = zesw/(1._dp-vtmpc1*zesw)
        zsusatw_2d(jl,jk) = pqm1(jl,jk)/zqsw_2d(jl,jk)-1.0_dp
        zsusatw_2d(jl,jk) = MAX(zsusatw_2d(jl,jk),0.0_dp)
        !SF for evaporation of rain in 3.3:
        zsusatw_evap(jl,jk) = pqm1(jl,jk)/zqsw_2d(jl,jk)-1.0_dp
        zsusatw_evap(jl,jk) = MIN(zsusatw_evap(jl,jk), 0._dp)

        !--- Saturation water vapour pressure:
        zesw_2d(jl,jk) = zesw*papm1(jl,jk)*rv/rd

        !SF for later use in 5.:
        zqswp1_2d(jl,jk) = tlucuaw(itm1p1_look(jl,jk))/papm1(jl,jk)
        zqswp1_2d(jl,jk) = MIN(zqswp1_2d(jl,jk),0.5_dp)
        zqswp1_2d(jl,jk) = zqswp1_2d(jl,jk)/(1._dp-vtmpc1*zqswp1_2d(jl,jk))   
        !SFend for later use in 5.

!SF ice:

        zqsi_2d(jl,jk)   = tlucua(itm1_look(jl,jk))/papm1(jl,jk)
        zqsi_2d(jl,jk)   = MIN(zqsi_2d(jl,jk),0.5_dp)
        zesi_2d(jl,jk)   = zqsi_2d(jl,jk)*papm1(jl,jk)*rv/rd
        zqsi_2d(jl,jk)   = zqsi_2d(jl,jk)/(1._dp-vtmpc1*zqsi_2d(jl,jk))
        sice(jl,jk,jrow) = pqm1(jl,jk)/zqsi_2d(jl,jk)-1.0_dp
        sice(jl,jk,jrow) = MAX(sice(jl,jk,jrow), 0._dp)

        !SF for sublimation of snow and ice in 3.2:
        zicesub(jl,jk) = pqm1(jl,jk)/zqsi_2d(jl,jk)-1.0_dp
        zicesub(jl,jk) = MIN(zicesub(jl,jk), 0._dp) 

        !SF for later use in 5.:
        zqsip1_2d(jl,jk) = tlucua(itm1p1_look(jl,jk))/papm1(jl,jk)
        zqsip1_2d(jl,jk) = MIN(zqsip1_2d(jl,jk),0.5_dp)
        zqsip1_2d(jl,jk) = zqsip1_2d(jl,jk)/(1._dp-vtmpc1*zqsip1_2d(jl,jk))
        !SFend for later use in 5.
 
     END DO !jl
  END DO !jk

  !--- Store supersaturation with respect to water in stream:
  swat(1:kproma,:,jrow)  = zsusatw_2d(1:kproma,:)

!SF Set some utility variables

  zastbstw(:,:) = alv * (alv/(rv*ptm1(:,:)) - 1.0_dp) / (zka*ptm1(:,:))  & !SF zast for water
                + rv*ptm1(:,:)/(2.21_dp/papm1(:,:)*zesw_2d(:,:))           !SF zbst for water

  ztmp1(:,:) = 4.1867e-3_dp*(5.69_dp + 0.017_dp*(ptm1(:,:)-tmelt)) !SF zkair

  zastbsti(:,:) = als * (als/(rv*ptm1(:,:)) - 1.0_dp) / (ztmp1(:,:)*ptm1(:,:))  & !SF zast for ice
                + rv*ptm1(:,:)/(2.21_dp/papm1(:,:)*zesi_2d(:,:))                     !SF zbst for ice

!--- Aerosol activation (so far only the Lin and Leaitch scheme is tested)


!--- Convert the aerosol activation into the number of newly formed cloud droplets

  ll1_2d(:,:) = ( ibas(:,:) > 0 )                                .AND.  &  !cloud base level
                ( (zcdnc(:,:)      <= cdncmin)             .OR.         &
                  (pxlte(:,:)      >  0._dp)               .OR.         &
                  (paclc(:,:)      >  cloud_tm1(1:kproma,:,jrow)) .OR.         &
                  (zsusatw_2d(:,:) >  zeps)                     )

  !SF: first computes newly formed cloud droplets at cloud bases:
  ztmp(:,:)      = pcdncact(:,:) - zcdnc(:,:) 
  ztmp(:,:)      = MAX(0._dp, ztmp(:,:))
  zqlnuc(:,:)    = MERGE(ztmp(:,:), 0._dp, ll1_2d(:,:))

  ztmp(:,:)      = zcdnc(:,:) + zqlnuc(:,:)
  zcdnc(:,:)     = MERGE(ztmp(:,:), zcdnc(:,:), ll1_2d(:,:))
  qnuc(1:kproma,:,jrow) = qnuc(1:kproma,:,jrow) + zdt*zqlnuc(1:kproma,:)
 
  !SF: then computes newly formed cloud droplets above cloud base
  !    by assuming that the number of nucleated cloud droplets is constant above cloud base
  !    (adiabatic parcel theory)

  DO jk=ktdia,klev
     DO jl=1,kproma
        jkk = MAX(jk,icl_minusbas(jl,jk))  !sets the level index to either relevant cloud base or itself
        zqlnuc_bas(jl,jk) = zqlnuc(jl,jkk) !holds in each cloud level the number of
                                           !newly formed cloud droplets at the base

        zcdnc_bas(jl,jk) = zcdnc(jl,jkk)   !holds in each cloud level the number of
                                           !activated cloud droplets at the base
     ENDDO !end jl
  ENDDO !end jk

  ll1_2d(:,:) = (icl_minusbas(:,:) > 0)     .AND. & !all cloud levels minus base
                (zqlnuc_bas(:,:)   > 0._dp)         !newly cloud droplets at base

  zqlnuc(:,:)    = MERGE(zqlnuc_bas(:,:), zqlnuc(:,:), ll1_2d(:,:))
  ztmp(1:kproma,:)      = qnuc(1:kproma,:,jrow) + zdt*(zcdnc_bas(1:kproma,:) - zcdnc(1:kproma,:))
  qnuc(1:kproma,:,jrow) = MERGE(ztmp(1:kproma,:), qnuc(1:kproma,:,jrow), ll1_2d(:,:))
  zcdnc(:,:)     = MERGE(zcdnc_bas(:,:), zcdnc(:,:), ll1_2d(:,:))

!--- obtain a number of cloud droplets for the detrained cloud water from convective anvils
  
  DO jk=ktdia,klev
     iclbas(:,jk) = NINT(pcvcbot(:)) 
     ll1_2d(:,jk) = (jk == iclbas(:,jk))
  ENDDO

  ll2_2d(:,:) = (iclbas(:,:) > 0) 

  IF (ncvmicro == 0) THEN  !no CONV micro scheme

     zqlnuccvh(:,:) = MERGE(cdncact_cv(1:kproma,:,jrow)-zcdnc(:,:), 0._dp, ll1_2d(:,:))

     ll3_2d(:,:) = ll2_2d(:,:)           .AND. &
                   (pxtec(:,:) > 0._dp)  .AND. &
                   (ptm1(:,:)  > cthomi)

     ztmp1(:,:) = cdncact_cv(1:kproma,:,jrow)-zcdnc(:,:)

     DO jk=ktdia,klev
        DO jl=1,kproma
           jkk = iclbas(jl,jk)
           jkk = MAX(1, jkk) !SF prevents cases where iclbas=0
           ztmp2(jl,jk) = zqlnuccvh(jl,jkk) 
        ENDDO
     ENDDO

     ztmp(:,:) = MIN(ztmp1(:,:),ztmp2(:,:))
     ztmp(:,:) = MAX(0._dp,ztmp(:,:))

  ELSEIF (ncvmicro > 0) THEN  !CONV micro scheme ON

     ll3_2d(:,:) = ll2_2d(:,:)           .AND. &
                   (pxtecl(:,:) > 0._dp)

     ztmp(:,:) = pxtecnl(:,:)-zcdnc(:,:)
     ztmp(:,:) = MAX(0._dp, ztmp(:,:))

  ENDIF !end ncvmicro

  zqlnuccv(:,:)  = MERGE(ztmp(:,:), 0._dp, ll3_2d(:,:))
  qnuc(1:kproma,:,jrow) = qnuc(1:kproma,:,jrow) + zdt*zqlnuccv(1:kproma,:)
  zcdnc(:,:)     = zcdnc(:,:)     +     zqlnuccv(:,:)

!  Ice nucleation

  DO jk=klev,ktdia,-1
     DO jl=1,kproma

        ztc = ptm1(jl,jk)-tmelt
        ztc = MIN(0._dp,ztc)
!
!       test using Faisal's param
!
        zrieff = 0.015_dp*ztc
        zrieff = EXP(zrieff)
        zrieff = 23.2_dp*zrieff
        zrieff = MAX(zrieff,1.0_dp)
     
        zrih = 5113188.044_dp+2809._dp*zrieff**3
        zrih = SQRT(zrih)
        zrih = -2261.236_dp+zrih

        zrid_2d(jl,jk) = 1.e-6_dp*zrih**(1._dp/3._dp)
        zrid_2d(jl,jk) = MAX(1.e-6_dp,zrid_2d(jl,jk))

     END DO !jl
  END DO !jk

  znidetr(:,:)=0._dp  

  IF (ncvmicro > 0) THEN !SF conv micro. scheme ON

     ll_cv(:,:) = (pxteci(:,:)  > 0._dp) 

     ztmp(:,:) = pxtecni(:,:) - zicncq(:,:)
     ztmp(:,:) = MAX(ztmp(:,:), 0._dp) 


  ELSEIF (ncvmicro == 0) THEN !SF conv micro. scheme OFF

     ll_cv(:,:) = (pxtec(:,:) > 0._dp) .AND. &
                  (ptm1(:,:)  < cthomi)

     ztmp(:,:) = 0.75_dp*ztmst*zrho(:,:)*pxtec(:,:)/(api*zrhoice*zrid_2d(:,:)**3) - zicncq(:,:)
     ztmp(:,:) = MAX(ztmp(:,:), 0._dp)

  ENDIF

  znidetr(:,:) = MERGE(ztmp(:,:),0._dp,ll_cv(:,:))
  zicncq(:,:)  = zicncq(:,:) + znidetr(:,:)

  ! Provide soluble aerosol number number concentratios for the cirrus scheme:
  zascs(:,:) = 0._dp
  ! Calculate the sum by explicitly using the involved modes
  ! (This is necessary for MADE3, to avoid including the mixed modes, 
  !  but it is also correct for GMXe)
  zascs(1:kproma,:) = mode(iaits)%aernum(1:kproma,:) + &
       mode(iaccs)%aernum(1:kproma,:) + &
       mode(icoas)%aernum(1:kproma,:)

  zascs(:,:) = MAX(zascs(:,:), 10.E6_dp*zrho(:,:))! change factor

!--- calculate the vertical velocity necessary for cirrus formation
!SF put this outside of the cirrus scheme choice, since this is also needed with nicnc=1

     ztmp(:,:)      = 70._dp*SQRT(ptkem1(:,:))
     ztmp(:,klev)   = 0.0_dp
     zvervx_2d(:,:) = -100._dp*g_rcp*pvervel(:,:)*zrho_rcp(:,:) + ztmp(:,:)

  IF ( nicnc == 1 ) THEN ! Use ICNC scheme after Lohmann, JAS, 2002
        
     ll_ice(:,:) = (sice(1:kproma,:,jrow) > 0._dp) .AND. &
                   (ptm1(:,:)      < cthomi)

     ztmp1(:,:) = 0.75_dp*zrho(:,:)*sice(1:kproma,:,jrow)*&
       zqsi_2d(:,:)/(api*zrhoice*zrid_2d(:,:)**3)-zicncq(:,:)
     ztmp2(:,:) = zascs(:,:)-zicncq(:,:)
     ztmp(:,:)  = MIN(ztmp1(:,:),ztmp2(:,:))
     ztmp(:,:)  = MAX(ztmp(:,:),0._dp)

     zninucl(:,:) = MERGE(ztmp(:,:), 0._dp, ll_ice(:,:))
     zicncq(:,:)   = zicncq(:,:) + zninucl(:,:)  

  ELSEIF ( nicnc > 1 ) THEN !--- Use Bernd's cirrus scheme 
!
!--- Kaercher and Lohmann, JGR, 2002b; Lohmann et al., JGR, 2004)
!--- This implies that supersaturation with respect to ice is allowed, thus the depositional
!--- growth equation needs to be solved
!


     zsusatix_2d(:,:) = sice(1:kproma,:,jrow)
     znicex(:,:)      = 0._dp

!--- the freezing rate is limited by the number of aerosols available in each mode. 
!    For homogeneous freezing, these are the soluble aerosols except for the nucleation mode

     IF (lhet) THEN

         DO jk = klev, ktdia, -1
             zapnx(:,1)   = 1.e-6_dp*(ndusol_strat(1:kproma,jk,jrow) - zicncq(:,jk))   ![1/cm3]
             zapnx(:,1)   = MAX(zapnx(:,1),1.e-6_dp)

             zaprx(:,1)   = 100._dp*mode(iaccm)%wetrad(1:kproma,jk) ![cm]
             zaprx(:,1)   = MAX(zaprx(:,1),mode(iaccs)%crdiv)

             zapsigx(:,1) = 1.0_dp

             zap(:,jk) = zapnx(:,1)

!----------- the freezing parameterization returns the number of newly formed ice crystals
!----------- (znicex) and their size (zri)

             CALL xfrzmstr(lhet, nosize, ztmst,                      &
                           klev, kbdim, kproma, jk, nfrzmod,         &
                           zsusatix_2d(:,jk), zvervx_2d(:,jk), zapnx,                  &
                           zaprx, zapsigx, ptm1, tmelt, zeps, papm1, &
                           cthomi, zri, znicex)

         ENDDO !jk

     ELSE

         DO jk = klev, ktdia, -1
  
            IF (.NOT. lomonodisp) THEN !SF standard 3-modes freezing
      
               zapnx(1:kproma,1)   = 1.e-6_dp                     &
                             * ( mode(iaits)%aernum(1:kproma,jk) &
                               - zicncq(:,jk) )![1/cm3]
                zapnx(:,1)   = MAX(zapnx(:,1),1.e-6_dp)

                zapnx(1:kproma,2)   = 1.e-6_dp                     &
                             * ( mode(iaccs)%aernum(1:kproma,jk) &
                               - zicncq(:,jk) )![1/cm3]
                zapnx(:,2)   = MAX(zapnx(:,2),1.e-6_dp)

                zapnx(1:kproma,3)   = 1.e-6_dp                     &
                             * ( mode(icoas)%aernum(1:kproma,jk) &
                               - zicncq(:,jk) )![1/cm3]
                zapnx(:,3)   = MAX(zapnx(:,3),1.e-6_dp)

                zaprx(1:kproma,1)   = 100._dp*mode(iaits)%wetrad(1:kproma,jk) ![cm]
                zaprx(:,1)   = MAX(zaprx(:,1),mode(iaits)%crdiv)

                zaprx(1:kproma,2)   = 100._dp*mode(iaccs)%wetrad(1:kproma,jk) ![cm]

                zaprx(:,2)   = MAX(zaprx(:,2),mode(iaccs)%crdiv)

                zaprx(1:kproma,3)   = 100._dp*mode(icoas)%wetrad(1:kproma,jk) ![cm]

                zaprx(:,3)   = MAX(zaprx(:,3),mode(icoas)%crdiv)

                zapsigx(:,1) = mode(iaits)%sigma
                zapsigx(:,2) = mode(iaccs)%sigma
                zapsigx(:,3) = mode(icoas)%sigma

                zap(:,jk) = zapnx(:,1) + zapnx(:,2) + zapnx(:,3)

            ELSE !SF monodisperse aerosol distribution

                zapnx(:,1)   = 1.e-6_dp * ( zascs(:,jk) - zicncq(:,jk) ) ![1/cm3]
                zapnx(:,1)   = MAX(zapnx(:,1),1.e-6_dp)

                zaprx(1:kproma,1)   = 100._dp*mode(iaccs)%wetrad(1:kproma,jk) ![cm]

                zaprx(:,1)   = MAX(zaprx(:,1),mode(iaccs)%crdiv)

                zapsigx(:,1) = 1._dp

                zap(:,jk) = zapnx(:,1) 

            ENDIF !SF end lomonodisp

!----------- the freezing parameterization returns the number of newly formed ice crystals
!----------- (znicex) and their size (zri)

              !alter Wolkenteil:

              DO jl = 1, kproma

                IF (ptm1(jl,jk) < cthomi .AND. paclc(jl,jk) .ge. pxtm1(jl,jk,4) .AND. &
                     pxtm1(jl,jk,4) .gt. zeps  .AND.           &
                     zvervx_2d(jl,jk) .gt. 0._dp) THEN

                  IF(zicnc(jl,jk) .gt. zicemin) THEN

                    !Bewoelkung gestiegen oder gleich, T<-35C
                    !modification of the vertical updraft because of
                    !pre-existing ice particles

                    !Volumenradius Kugel
                    zri0 = (0.75_dp * zrho(jl,jk) * (pxim1(jl,jk)+ztmst*pxite(jl,jk))  &
                             /(api*zrhoice*zicnc(jl,jk)*paclc(jl,jk)))**(1._dp/3._dp)
                    zri0 = MAX(zri0, 0._dp)

                    !critical ice saturation ratio for homogenous freezing
                    scr = SCRHOM(ptm1(jl,jk))

                    CALL XVICE(ztmst, zvervx_2d(jl,jk), papm1(jl,jk)*0.01_dp,  &
                                ptm1(jl,jk), zri0*1.e2_dp, zicnc(jl,jk)*1.e-6_dp, &
                                scr, vice)
                    !vice = 0._dp

                    !modifizierte Vertikalgeschwindigkeit
                    zvervxmod_2d(jl,jk) = zvervx_2d(jl,jk) - vice

                  ELSE

                    zvervxmod_2d(jl,jk) = zvervx_2d(jl,jk)

                  END IF

                END IF

              END DO

              WHERE (zvervxmod_2d(:,jk) >= 1.e-2_dp .AND.  &
                      zvervxmod_2d(:,jk) <= zvervx_2d(:,jk))

                !w groesser 0.1 mm/s
                !Sekundaernukleation mit Ueberpruefung, ob
                !Vertikalgeschwindigkeit ausreicht von
                !Saettigung zur kritischen Feuchte zu gelangen

                zsusatix_2d(:,jk) = 1._dp

              ELSEWHERE

                zsusatix_2d(:,jk) = 0._dp

              END WHERE

              CALL XFRZMSTR(lhet, nosize, ztmst, klev, kbdim, kproma, jk, nfrzmod, &
                      zsusatix_2d(:,jk), zvervxmod_2d(:,jk), zapnx,                &
                      zaprx, zapsigx, ptm1, tmelt, zeps, papm1, &
                      cthomi, zri, znicex)

              zri1(:,jk) = zri(:,jk)
              znicex1(:,jk) = znicex(:,jk)

              !neuer Wolkenteil:

              DO jl = 1, kproma

                IF ( ptm1(jl,jk) < cthomi .AND. (paclc (jl,jk) >= zeps .AND.    &
                      paclc(jl,jk) > pxtm1(jl,jk,4)) ) THEN
                      !Bedeckungsgrad steigt

                  IF (zvervx_2d(jl,jk) >= 1.e-2_dp) THEN  !groesser 0.1 mm/s

                    !critical ice saturation ratio for homogenous freezing
                    scr = SCRHOM(ptm1(jl,jk))
                    !vapor pressure over ice in mbar
                    pice = PISAT(ptm1(jl,jk))
                    pw = scr * pice
                    coolr = 981.0_dp * zvervx_2d(jl,jk) / 1.00467E7_dp

                    !nucleation
                    zri0 = 0.0_dp
                    znicex0 = 0.0_dp
                    CALL XFRZHOM(nosize, nfrzmod, ztmst, zapnx(jl,:), zaprx(jl,:), &
                                  zapsigx(jl,:), papm1(jl,jk)*0.01_dp, ptm1(jl,jk), &
                                  zvervx_2d(jl,jk), scr, znicex0, zri0, coolr, pw, &
                                  scr, ptm1(jl,jk), pice, pw )
                    zri2(jl,jk) = zri0 * 1.0E-2_dp         !Umrechnung in m
                    znicex2(jl,jk) = MAX(znicex0 * 1.0E6_dp,0._dp) !Umrechnung in /m3

                  END IF

                END IF

              END DO

              !Anteil neue Wolke
              ztmp1(:,jk) = (paclc(:,jk)-pxtm1(:,jk,4))/MAX(paclc(:,jk),zeps)
              ztmp1(:,jk) = MERGE(ztmp1(:,jk),0._dp,paclc(:,jk).gt.zeps)

              !neuer mittlerer in-cloud Wert
              znicex(:,jk) = (1._dp - ztmp1(:,jk)) * znicex1(:,jk)  &
                                 + ztmp1(:,jk) * znicex2(:,jk)
              zri(:,jk) = (1._dp - ztmp1(:,jk)) * zri1(:,jk) &
                                 + ztmp1(:,jk) * zri2(:,jk)

         ENDDO !jk

     ENDIF !lhet

!--- Update ICNC

     zri(:,:)=MAX(zri(:,:), 1.e-6_dp)
     
     ll_ice(:,:) = (ptm1(:,:)        < cthomi)

     ztmp1(:,:) = 1.e6_dp*zap(:,:)                     !SB: [ztmp1]=[1/m3], [zap]=[1/cm3] 

     ztmp2(:,:) = znicex(:,:)
     ztmp(:,:)  = MIN(ztmp1(:,:),ztmp2(:,:)) 
     ztmp(:,:)  = MAX(ztmp(:,:), 0._dp)

     zninucl(:,:) = MERGE(ztmp(:,:), 0._dp,ll_ice(:,:))
     zicncq(:,:)  = zicncq(:,:) + zninucl(:,:)

!---calculate the deposition rate taking ventilation (zfre) into account

     ll_ice(:,:) = (pxim1(:,:) > 0._dp)

     ztmp(:,:)   = MAX(pxtm1(:,:,2), zicemin)
     zicncqi(:,:) = MERGE(ztmp(:,:), zicncq(:,:), ll_ice(:,:))

     ztmp(:,:) = zrho(:,:)*pxim1(:,:)/zicncqi(:,:)

     ztmp(:,:) = zrho(:,:)*pxim1(:,:)/zicncq(:,:)
     ztmp(:,:) = MAX(ztmp(:,:), zmi)
     zmmean_2d(:,:) = MERGE(ztmp(:,:), zmmean_2d(:,:), ll_ice(:,:))

     ll1_2d(:,:) = (zmmean_2d(:,:) < 2.166E-9_dp )
     ll2_2d(:,:) = (zmmean_2d(:,:) >= 2.166E-9_dp  ) .AND. (zmmean_2d(:,:) < 4.264E-8_dp )

     zalfased_2d(:,:) = MERGE(63292.4_dp, 8.78_dp, ll1_2d(:,:))
     zalfased_2d(:,:) = MERGE(329.75_dp, zalfased_2d(:,:), ll2_2d(:,:))

     zbetased_2d(:,:) = MERGE(0.5727_dp, 0.0954_dp, ll1_2d(:,:))
     zbetased_2d(:,:) = MERGE(0.3091_dp, zbetased_2d(:,:), ll2_2d(:,:))

     zxifallmc_2d(:,:) = zfall*zalfased_2d(:,:) &
                        *(zmmean_2d(:,:)**zbetased_2d(:,:))*zaaa_2d(:,:) !Fallgeschwindigkeit Masse

     ztmp(:,:)       = pqm1(:,:) - zqsi_2d(:,:)
     zsusatix_2d(:,:) = pqm1(:,:)+pqte(:,:) - zqsi_2d(:,:)

     DO jk = ktdia, klev
        DO jl = 1,kproma
           zdv     = 2.21_dp/papm1(jl,jk)
           zgtp    = 1._dp/(zrho(jl,jk)*zastbsti(jl,jk))
           zvth    = SQRT( 8._dp*zkb*ptm1(jl,jk) / (api*zxmw) )
           zb2     = 0.25_dp * zalpha * zvth   / zdv
           zfuchs  = 1._dp/(1._dp+zb2*zri(jl,jk))
           zre     = 2._dp*zrho(jl,jk)*zri(jl,jk)*zxifallmc_2d(jl,jk)/zviscos_2d(jl,jk)
           zfre    = 1._dp + 0.229_dp*SQRT(zre)
           zfre    = MERGE(zfre, 1._dp, ll_ice(jl,jk))

           zqinucl(jl,jk)  = 4._dp*api*zri(jl,jk)*zsusatix_2d(jl,jk)*zicncqi(jl,jk) &
                           *zfre*zgtp*zfuchs*zalpha*ztmst
           zqinucl(jl,jk) = MAX(zqinucl(jl,jk),0._dp)

        ENDDO !jl
     ENDDO !jk
 
     ztmp(:,:)   = MAX(znicex2(:,:), zicemin)
     zicncqi(:,:) = MERGE(ztmp(:,:), znicex2(:,:), ll_ice(:,:))

     !Radius 5 um
     ztmp(:,:) = 4._dp/3._dp * 5.e-6_dp ** 3._dp * api * crhoi
     ztmp(:,:) = MAX(ztmp(:,:), zmi)
     zmmean_2d(:,:) = MERGE(ztmp(:,:), zmmean_2d(:,:), ll_ice(:,:))

     ll1_2d(:,:) = (zmmean_2d(:,:) < 2.166E-9_dp )
     ll2_2d(:,:) = (zmmean_2d(:,:) >= 2.166E-9_dp  ) .AND. (zmmean_2d(:,:) < 4.264E-8_dp )

     zalfased_2d(:,:) = MERGE(63292.4_dp, 8.78_dp, ll1_2d(:,:))
     zalfased_2d(:,:) = MERGE(329.75_dp, zalfased_2d(:,:), ll2_2d(:,:))

     zbetased_2d(:,:) = MERGE(0.5727_dp, 0.0954_dp, ll1_2d(:,:))
     zbetased_2d(:,:) = MERGE(0.3091_dp, zbetased_2d(:,:), ll2_2d(:,:))

     zxifallmc_2d(:,:) = zfall*zalfased_2d(:,:) &
                        *(zmmean_2d(:,:)**zbetased_2d(:,:))*zaaa_2d(:,:) !Fallgeschwindigkeit Masse

     DO jk = ktdia, klev
        DO jl = 1,kproma
           zdv     = 2.21_dp/papm1(jl,jk)
           zgtp    = 1._dp/(zrho(jl,jk)*zastbsti(jl,jk))
           zvth    = SQRT( 8._dp*zkb*ptm1(jl,jk) / (api*zxmw) )
           zb2     = 0.25_dp * zalpha * zvth   / zdv
           zfuchs  = 1._dp/(1._dp+zb2*zri(jl,jk))
           zre     = 2._dp*zrho(jl,jk)*zri(jl,jk)*zxifallmc_2d(jl,jk)/zviscos_2d(jl,jk)
           zfre    = 1._dp + 0.229_dp*SQRT(zre)
           zfre    = MERGE(zfre, 1._dp, ll_ice(jl,jk))

           zqinucl_new(jl,jk)  = 4._dp*api*5.e-6_dp*zsusatix_2d(jl,jk)*zicncqi(jl,jk) &
                           *zfre*zgtp*zfuchs*zalpha*ztmst
           zqinucl_new(jl,jk) = MAX(zqinucl_new(jl,jk),0._dp)

        ENDDO !jl
     ENDDO !jk

  ENDIF   !which nucleation scheme (nicnc)

  !corinna: set zcdnc and zicnc to minium now if nucleation is not strong enough

  ll1(:,:) = ( paclc(:,:) >= zepsec ) .AND. &
             ( ptm1(:,:)  >  cthomi )

  ztmp(:,:)  = MAX(zcdnc(:,:),cdncmin)
  zcdnc(:,:) = MERGE(ztmp(:,:),zcdnc(:,:),ll1(:,:))

  ll1(:,:) = ( paclc(:,:) >= zepsec ) .AND. &
             ( ptm1(:,:)  <  tmelt  )

  ztmp(:,:)   = MAX(zicncq(:,:), zicemin)
  zicncq(:,:) = MERGE(ztmp(:,:), zicncq(:,:), ll1(:,:))

!--- End included for CDNC/IC scheme -----------------------------------

!--- CCMod -------------------------------------------------------------

  !-------------------------
  !flight inventory
  !-------------------------

  zfkme(:,:) = 0._dp
  zfh2o(:,:) = 0._dp

  zfkme(:,:) = ccfkme_invent(:,:,jrow) * ztmst !km
  zfh2o(:,:) = ccfh2o_invent(:,:,jrow) * ztmst            !kg/kg 

  !Umrechnung flight dist in Flugverkehrsdichte km/(s*m2)
  !und Wasserdampfausstoss in kg/m3
  ccfkme(:,:,jrow) = 1.e-3_dp * ccfkme_invent(:,:,jrow)
  ccfh2o(:,:,jrow) = ccfh2o_invent(:,:,jrow)*zrho(:,:)


  !-------------------------
  !tracer
  !-------------------------

  zcccov_1dt(:,:) = MAX(pxtm1(:,:,5), 0._dp)
  zcccov_2dt(:,:) = MAX(pxtm1(:,:,6), 0._dp)
  zcccov_3dt(:,:) = MAX(pxtm1(:,:,7), 0._dp)
  zcccov_4dt(:,:) = MAX(pxtm1(:,:,8), 0._dp)
  zcccov_5dt(:,:) = MAX(pxtm1(:,:,9), 0._dp)
  zcccov(:,:)     = MAX(pxtm1(:,:,10), 0._dp)

  zccvol_1dt(:,:) = MAX(pxtm1(:,:,11), 0._dp)
  zccvol_2dt(:,:) = MAX(pxtm1(:,:,12), 0._dp)
  zccvol_3dt(:,:) = MAX(pxtm1(:,:,13), 0._dp)
  zccvol_4dt(:,:) = MAX(pxtm1(:,:,14), 0._dp)
  zccvol_5dt(:,:) = MAX(pxtm1(:,:,15), 0._dp)
  zccvol(:,:)     = MAX(pxtm1(:,:,16), 0._dp)

  zcclen_1dt(:,:) = MAX(pxtm1(:,:,17), 0._dp)
  zcclen_2dt(:,:) = MAX(pxtm1(:,:,18), 0._dp)
  zcclen_3dt(:,:) = MAX(pxtm1(:,:,19), 0._dp)
  zcclen_4dt(:,:) = MAX(pxtm1(:,:,20), 0._dp)
  zcclen_5dt(:,:) = MAX(pxtm1(:,:,21), 0._dp)
  zcclen(:,:)     = MAX(pxtm1(:,:,22), 0._dp)

  zcclen_ou(:,:)  = MAX(pxtm1(:,:,23), 0._dp)

  !Abgleich gescherter Bedeckungsgrad
  zcccov_1dt(:,:) = MAX(zccvol_1dt(:,:),zcccov_1dt(:,:))
  zcccov_2dt(:,:) = MAX(zccvol_2dt(:,:),zcccov_2dt(:,:))
  zcccov_3dt(:,:) = MAX(zccvol_3dt(:,:),zcccov_3dt(:,:))
  zcccov_4dt(:,:) = MAX(zccvol_4dt(:,:),zcccov_4dt(:,:))
  zcccov_5dt(:,:) = MAX(zccvol_5dt(:,:),zcccov_5dt(:,:))
  zcccov(:,:)     = MAX(zccvol(:,:),zcccov(:,:))
 
  zcccov_ges(:,:) = zcccov_1dt(:,:) + zcccov_2dt(:,:) + zcccov_3dt(:,:) &
                    + zcccov_4dt(:,:) + zcccov_5dt(:,:) + zcccov(:,:)

  zccvol_ges(:,:) = zccvol_1dt(:,:) + zccvol_2dt(:,:) + zccvol_3dt(:,:) &
                    + zccvol_4dt(:,:) + zccvol_5dt(:,:) + zccvol(:,:)

  zcclen_ges(:,:) = zcclen_1dt(:,:) + zcclen_2dt(:,:) + zcclen_3dt(:,:) &
                    + zcclen_4dt(:,:) + zcclen_5dt(:,:) + zcclen(:,:)

  zccicnc(:,:)    = MAX(zrho(:,:)*(pxtm1(:,:,24)/m_air  * 1000._dp), 0._dp)
  zccicnc(:,:)    = zccicnc(:,:) / MAX(zccvol_ges(:,:), zeps)
  zcciwc(:,:)     = MAX(pxtm1(:,:,25), 0._dp)


  zumf = 1.e-30_dp
  WHERE(zcccov_ges(:,:).le.zumf .OR. zccvol_ges(:,:).le.zumf .OR. zcclen_ges(:,:).le.zumf .OR. &
      zcciwc(:,:).le.zumf .OR. zccicnc(:,:).le.zumf) 
    zcccov_ges(:,:) = 0._dp
    zccvol_ges(:,:) = 0._dp
    zcclen_ges(:,:) = 0._dp
    zcciwc(:,:)     = 0._dp
    zccicnc(:,:)    = 0._dp
  END WHERE  

  !Begrenzung der Contrailcirrus-Bedeckung potentiellen Bedeckungsgrad
  !     Korrektur nach Zeitschritt:
  !     Der oben berechnete Bedeckungsgrad fuer Contrailcirrus darf
  !     den potentiellen Bedeckungsgrad nicht uebersteigen.
  !     Wenn korrigiert werden muss, muss folgender Faktor zcor44
  !     fuer die Korrektur der Eiswassergehalte gemerkt werden.
  
  !Bedeckungsgrade werden zu gleichem Anteil reduziert
  !berechneter Bedeckungsgrad darf nicht kleiner 0 sein
  where ((zcccov_ges(:,:) - B_cc(1:kproma,:,jrow)) .gt. 0._dp)

    zcor44(:,:) = (zcccov_ges(:,:)-B_cc(1:kproma,:,jrow))/MAX(zcccov_ges(:,:),zeps)
    zmult(:,:) = 1._dp - zcor44(:,:)
  
    zccvol_ges(:,:) = MAX(zccvol_ges(:,:) * zmult(:,:), 0._dp)
    zccvol_1dt(:,:) = MAX(zccvol_1dt(:,:) * sqrt(zmult(:,:)), 0._dp)
    zccvol_2dt(:,:) = MAX(zccvol_2dt(:,:) * sqrt(zmult(:,:)), 0._dp)
    zccvol_3dt(:,:) = MAX(zccvol_3dt(:,:) * sqrt(zmult(:,:)), 0._dp)
    zccvol_4dt(:,:) = MAX(zccvol_4dt(:,:) * sqrt(zmult(:,:)), 0._dp)
    zccvol_5dt(:,:) = MAX(zccvol_5dt(:,:) * sqrt(zmult(:,:)), 0._dp)
    zccvol(:,:)     = MAX(zccvol(:,:) * sqrt(zmult(:,:)), 0._dp)
  
    zcccov_ges(:,:) = MAX(zcccov_ges(:,:) * zmult(:,:), 0._dp)
    zcccov_1dt(:,:) = MAX(zcccov_1dt(:,:) * sqrt(zmult(:,:)), 0._dp)
    zcccov_2dt(:,:) = MAX(zcccov_2dt(:,:) * sqrt(zmult(:,:)), 0._dp)
    zcccov_3dt(:,:) = MAX(zcccov_3dt(:,:) * sqrt(zmult(:,:)), 0._dp)
    zcccov_4dt(:,:) = MAX(zcccov_4dt(:,:) * sqrt(zmult(:,:)), 0._dp)
    zcccov_5dt(:,:) = MAX(zcccov_5dt(:,:) * sqrt(zmult(:,:)), 0._dp)
    zcccov(:,:)     = MAX(zcccov(:,:) * sqrt(zmult(:,:)), 0._dp)
  
    zcclen_ges(:,:) = MAX(zcclen_ges(:,:) * sqrt(zmult(:,:)), 0._dp)
    zcclen_1dt(:,:) = MAX(zcclen_1dt(:,:) * sqrt(zmult(:,:)), 0._dp)
    zcclen_2dt(:,:) = MAX(zcclen_2dt(:,:) * sqrt(zmult(:,:)), 0._dp)
    zcclen_3dt(:,:) = MAX(zcclen_3dt(:,:) * sqrt(zmult(:,:)), 0._dp)
    zcclen_4dt(:,:) = MAX(zcclen_4dt(:,:) * sqrt(zmult(:,:)), 0._dp)
    zcclen_5dt(:,:) = MAX(zcclen_5dt(:,:) * sqrt(zmult(:,:)), 0._dp)
    zcclen(:,:)     = MAX(zcclen(:,:) * sqrt(zmult(:,:)), 0._dp)
  
    zcclen_ou(:,:)  = MAX(zcclen_ou(:,:) * sqrt(zmult(:,:)), 0._dp)
  
    zcor44(:,:) = zcor44(:,:) * zcciwc(:,:)
    zcciwc(:,:) = MAX(zcciwc(:,:) * zmult(:,:), 0._dp)
  elsewhere
    zcor44(:,:) = 0._dp
    zmult(:,:) = 1._dp
  end where


  !-------------------------------
  !Contrail formation
  !-------------------------------
  
  !initial contrail width 200m
  zwidthini = 200._dp
  !initial contrail height 200m
  zheightini = 200._dp
  !Area GB in m**2
  DO jk = 1, klev
    zflaegb(:,jk) = gboxarea_2d(:,jrow)
  END DO
  !Maximale initiale Bedeckung zconcovinitial
  where (B_co(1:kproma,:,jrow) .gt. 0._dp)
    zconcovinitial(:,:) = MIN(zfkme(:,:)*zwidthini/zflaegb(:,:), 1._dp)
  elsewhere
    zconcovinitial(:,:) = 0._dp
  end where

  !b_co kann nur kleiner/gleich als pot B_co sein
  !und Contrailcirren existieren teilweise schon in B_co
  where (B_co(1:kproma,:,jrow) .gt. 0._dp .AND. B_cc(1:kproma,:,jrow) .gt. 0._dp)
    zfrac1(:,:) = B_co(1:kproma,:,jrow) * (1._dp-zccvol_ges(:,:) &
                                /MAX(B_cc(1:kproma,:,jrow),zeps))
    zcccov_new(:,:) = zconcovinitial(:,:) * zfrac1(:,:)
  elsewhere
    zfrac1(:,:) = 0._dp
    zcccov_new(:,:) = 0._dp
  end where
  
  !initial contrail volume
  zccvol_new(:,:) = zcccov_new(:,:) * MIN(zheightini/zdz_2d(:,:), 1._dp)

  !initial contrail length
  zcclen_new(:,:) = zcccov_new(:,:) / zwidthini
  
  !initial contrail ice crystal number concentration
  where(zccvol_new(:,:).gt.0._dp)
    zccicnc_new(:,:) = 150.e6_dp
  elsewhere
    zccicnc_new(:,:) = 0._dp
  end where

  !initial contrail ice water content
  where(zccvol_new(:,:).gt.0._dp)
    zcciwc_new(:,:) = zfrac1(:,:) * zfh2o(:,:)
  elsewhere
    zcciwc_new(:,:) = 0._dp
  end where
  zcciwc(:,:) = zcciwc(:,:) + (B_co(1:kproma,:,jrow) - zfrac1(:,:)) * zfh2o(:,:)

  !deposition of water when ice supersaturation in initial contrails
  zdepco_form(:,:) = 0._dp
  where(zccvol_new(:,:).gt.0._dp)
    ztmp2(:,:)=crt/(0.8_dp*(2.349_dp - ptm1(:,:)/259._dp))
    where(paclc(:,:).gt.0._dp)
      ztmp1(:,:)=0.5_dp*(crt-ztmp2(:,:))*zqsi_2d(:,:)
    elsewhere
      ztmp1(:,:)=0.5_dp*(pqm1(:,:)-ztmp2(:,:)*zqsi_2d(:,:))
    end where
    zdepco_form(:,:)=MAX(ztmp1(:,:),0._dp)*zccvol_new(:,:)
    zcciwc_new(:,:) = zcciwc_new(:,:) + zdepco_form(:,:)
  end where

  !Beschraenkung des gescherten Bedeckungsgrads neuer Kondensstreifen, 
  !falls sich diese unter- oder oberhalb der alten Kondensstreifen bilden
  
  where (B_cc(1:kproma,:,jrow) .gt. zccvol_ges(:,:))
    ztmp1(:,:) = (B_cc(1:kproma,:,jrow)-zcccov_ges(:,:))/(B_cc(1:kproma,:,jrow)-zccvol_ges(:,:))
    ztmp1(:,:) = MIN(MAX(ztmp1(:,:),0._dp),1._dp)
    zcccov_new(:,:) = ztmp1(:,:) * zcccov_new(:,:)
    zcclen_new(:,:) = ztmp1(:,:) * zcclen_new(:,:)
    !Laenge der KS die sich ober oder unterhalb bilden
    zcclen_ou(:,:) = zcclen_ou(:,:)+(1._dp-ztmp1(:,:)) * zcclen_new(:,:)
  end where
  

  !------------------
  !Spreading
  !------------------

  zcccov_ori(:,:) = zcccov_ges(:,:)

  !siehe Burkhardt und Kaercher - Tuning for spreading
  zztune2 = 0.72_dp * 0.75_dp
  !delta uv [m/s] (at time step - not accumulated)
  DO jk=ktdia,klev
    zdeltu(:,jk) = (pum1(:,jk)-pum1(:,jk+1))**2 + (pvm1(:,jk)-pvm1(:,jk+1))**2
  END DO
  zdeltu(:,:) = sqrt(MAX(zdeltu(:,:), 1.e-5_dp))

  where (B_cc(1:kproma,:,jrow) > 0._dp)

    zcccov_new(:,:) = zcccov_new(:,:)  &
                   + zdeltu(:,:)*zheightini/zdz_2d(:,:)*zztune2*ztmst*0.75_dp &
                   *(B_cc(1:kproma,:,jrow)-(zcccov_new(:,:)+zcccov_ges(:,:)))*zcclen_new(:,:)
    zcccov_ges(:,:) = zcccov_ges(:,:)  &
                   + zdeltu(:,:)*zztune2*ztmst &
                   *(B_cc(1:kproma,:,jrow)-(zcccov_new(:,:)+zcccov_ges(:,:)))*zcclen_ges(:,:)
  
  end where
  
  
  !Begrenzung auf potentiellen Bedeckungsgrad
  zcovini(:,:) = zcccov_new(:,:)+zcccov_ges(:,:)
  zcovini_k(:,:) = MIN(zcovini(:,:), B_cc(1:kproma,:,jrow))
  
  where (zcovini(:,:) .gt. zeps)
    zcorcov(:,:) = zcovini_k(:,:) / zcovini(:,:)
  elsewhere
    zcorcov(:,:) = 1._dp
  end where
  
  
  zcccov_new(:,:) = MAX(zcccov_new(:,:)*zcorcov(:,:), 0._dp)
  zcclen_new(:,:) = MAX(zcclen_new(:,:)*sqrt(zcorcov(:,:)), 0._dp)
  zcccov_ges(:,:) = MAX(zcccov_ges(:,:)*zcorcov(:,:), 0._dp)
  zccvol_ges(:,:) = MAX(zccvol_ges(:,:)*zcorcov(:,:), 0._dp)
  zccvol_new(:,:) = MAX(zccvol_new(:,:)*zcorcov(:,:), 0._dp)
  zcclen_ges2(:,:)= MAX(zcclen_ges(:,:)*sqrt(zcorcov(:,:)), 0._dp)
  zcciwc_new(:,:) = MAX(zcciwc_new(:,:)*zcorcov(:,:), 0._dp)
  zcciwc(:,:)     = MAX(zcciwc(:,:)*zcorcov(:,:), 0._dp)
  zcor44(:,:) = zcor44(:,:)+(1._dp-zcorcov(:,:))*(zcciwc_new(:,:)+zcciwc(:,:))
  
  !Verteilung der Raten auf die verschiedenen contrail Klassen
  !gescherter Bedeckungsgrad
  where (zcccov_ges(:,:) .gt. 0._dp)
    zmult2(:,:) = zcccov_ges(:,:) - zcccov_ori(:,:)
  elsewhere
    zmult2(:,:) = 0._dp
  end where
  where (zcclen_ges(:,:) .gt. zeps)
    zlmit(:,:) = 1._dp/zcclen_ges(:,:)
    zcccov_1dt(:,:) = zcccov_1dt(:,:)+zmult2(:,:)*zcclen_1dt(:,:)*zlmit(:,:)
    zcccov_2dt(:,:) = zcccov_2dt(:,:)+zmult2(:,:)*zcclen_2dt(:,:)*zlmit(:,:)
    zcccov_3dt(:,:) = zcccov_3dt(:,:)+zmult2(:,:)*zcclen_3dt(:,:)*zlmit(:,:)
    zcccov_4dt(:,:) = zcccov_4dt(:,:)+zmult2(:,:)*zcclen_4dt(:,:)*zlmit(:,:)
    zcccov_5dt(:,:) = zcccov_5dt(:,:)+zmult2(:,:)*zcclen_5dt(:,:)*zlmit(:,:)
    zcccov(:,:)     = zcccov(:,:)    +zmult2(:,:)*zcclen(:,:)    *zlmit(:,:)
  elsewhere
    zcccov_1dt(:,:) = 0._dp
    zcccov_2dt(:,:) = 0._dp
    zcccov_3dt(:,:) = 0._dp
    zcccov_4dt(:,:) = 0._dp
    zcccov_5dt(:,:) = 0._dp
    zcccov(:,:)     = 0._dp
  end where
  
  zccvol_ges_ori(:,:) = zccvol_ges(:,:)
  zccvol_new_ori(:,:) = zccvol_new(:,:)
  zcccov_ori(:,:) = zcccov_ges(:,:)

  !------------------
  !Turbulent diffusion
  !------------------

  where (zccvol_new(:,:) > 0._dp)
  
    ztmp1(:,:) = zdeltu(:,:)/zdz_2d(:,:)
    ztmp2(:,:) = ztmst*0.5_dp
    ztmp(:,:) = 2._dp*api*sqrt(0.0075_dp*ztmp1(:,:)**2*ztmp2(:,:)**4 &
                  + 444.9_dp*ztmp1(:,:)**2*ztmp2(:,:)**3 + 1334.7_dp*ztmp2(:,:) + 19.8e6_dp)
    ztmp1(:,:) = 1._dp - MAX(0._dp,MIN(1._dp,(0.75_dp*zccvol_new(:,:)+zccvol_ges(:,:))/(0.75_dp*zcccov_new(:,:)+zcccov_ges(:,:))))
    ztmp2(:,:) = 2._dp*api*sqrt(19.8e6_dp) 
    zccvol_new(:,:) = zccvol_new(:,:) + 0.5_dp*ztmp1(:,:)*(ztmp(:,:)-ztmp2(:,:))/ztmp2(:,:)*zccvol_new(:,:) 
    zccvol_new(:,:) = MIN(zccvol_new(:,:),zcccov_new(:,:))
  
  endwhere
  
  where (zccvol_ges(:,:) > 0._dp)
  
    ztmp1(:,:) = zdeltu(:,:)/zdz_2d(:,:)
    ztmp2(:,:) = ztmst*0.5_dp
    ztmp(:,:) = 2._dp*api*sqrt(0.0075_dp*ztmp1(:,:)**2*ztmp2(:,:)**4 &
                  + 444.9_dp*ztmp1(:,:)**2*ztmp2(:,:)**3 + 1334.7_dp*ztmp2(:,:) + 19.8e6_dp)
    ztmp1(:,:) = 1._dp - MAX(0._dp,MIN(1._dp,(0.75_dp*zccvol_new(:,:)+zccvol_ges(:,:))/(0.75_dp*zcccov_new(:,:)+zcccov_ges(:,:))))
    zdeltu(:,:) = zdeltu(:,:) / zdz_2d(:,:)
  
    ztmp2(:,:) = 2._dp*api*sqrt(0.0075_dp*zdeltu(:,:)**2*(2._dp*ztmst)**4 &
                  + 444.9_dp*zdeltu(:,:)**2*(2._dp*ztmst)**3 + 1334.7_dp*(2._dp*ztmst) + 19.8e6_dp)
    zccvol_1dt(:,:) = zccvol_1dt(:,:) + ztmp1(:,:)*(ztmp2(:,:)-ztmp(:,:))/ztmp(:,:)*zccvol_1dt(:,:)
    zccvol_1dt(:,:) = MIN(zccvol_1dt(:,:),zcccov_1dt(:,:))
  
    ztmp(:,:) = 2._dp*api*sqrt(0.0075_dp*zdeltu(:,:)**2*(3._dp*ztmst)**4 &
                  + 444.9_dp*zdeltu(:,:)**2*(3._dp*ztmst)**3 + 1334.7_dp*(3._dp*ztmst) + 19.8e6_dp)
    zccvol_2dt(:,:) = zccvol_2dt(:,:) + ztmp1(:,:)*(ztmp(:,:)-ztmp2(:,:))/ztmp2(:,:)*zccvol_2dt(:,:)
    zccvol_2dt(:,:) = MIN(zccvol_2dt(:,:),zcccov_2dt(:,:))
  
    ztmp2(:,:) = 2._dp*api*sqrt(0.0075_dp*zdeltu(:,:)**2*(4._dp*ztmst)**4 &
                  + 444.9_dp*zdeltu(:,:)**2*(4._dp*ztmst)**3 + 1334.7_dp*(4._dp*ztmst) + 19.8e6_dp)
    zccvol_3dt(:,:) = zccvol_3dt(:,:) + ztmp1(:,:)*(ztmp2(:,:)-ztmp(:,:))/ztmp(:,:)*zccvol_3dt(:,:)
    zccvol_3dt(:,:) = MIN(zccvol_3dt(:,:),zcccov_3dt(:,:))
  
    ztmp(:,:) = 2._dp*api*sqrt(0.0075_dp*zdeltu(:,:)**2*(5._dp*ztmst)**4 &
                  + 444.9_dp*zdeltu(:,:)**2*(5._dp*ztmst)**3 + 1334.7_dp*(5._dp*ztmst) + 19.8e6_dp)
    zccvol_4dt(:,:) = zccvol_4dt(:,:) + ztmp1(:,:)*(ztmp(:,:)-ztmp2(:,:))/ztmp2(:,:)*zccvol_4dt(:,:)
    zccvol_4dt(:,:) = MIN(zccvol_4dt(:,:),zcccov_4dt(:,:))
  
    ztmp2(:,:) = 2._dp*api*sqrt(0.0075_dp*zdeltu(:,:)**2*(6._dp*ztmst)**4 &
                  + 444.9_dp*zdeltu(:,:)**2*(6._dp*ztmst)**3 + 1334.7_dp*(6._dp*ztmst) + 19.8e6_dp)
    zccvol_5dt(:,:) = zccvol_5dt(:,:) + ztmp1(:,:)*(ztmp2(:,:)-ztmp(:,:))/ztmp(:,:)*zccvol_5dt(:,:)
    zccvol_5dt(:,:) = MIN(zccvol_5dt(:,:),zcccov_5dt(:,:))
  
    zccvol(:,:) = zccvol(:,:) + ztmp1(:,:)*0.2_dp*zccvol(:,:)
    zccvol(:,:) = MIN(zccvol(:,:),zcccov(:,:))
  
  endwhere
  
  zccvol_ges(:,:) = zccvol_1dt(:,:) + zccvol_2dt(:,:) + zccvol_3dt(:,:) + zccvol_4dt(:,:)  &
                    + zccvol_5dt(:,:) + zccvol(:,:)
  
  WHERE(zccvol_new(:,:).gt.1.e-10_dp)
    zccicnc_new(:,:) = zccvol_new_ori(:,:)/zccvol_new(:,:) * zccicnc_new(:,:)
  ELSEWHERE
    zccicnc_new(:,:) = 0._dp
  END WHERE
  WHERE(zccvol_ges(:,:).gt.1.e-10_dp)
    zccicnc(:,:) = MAX(0._dp,MIN(1._dp,zccvol_ges_ori(:,:)/zccvol_ges(:,:))) * zccicnc(:,:)
  ELSEWHERE
    zccicnc(:,:) = 0._dp
  END WHERE
  
  zccvol_ges_ori(:,:) = zccvol_ges(:,:)
  

!--- End CCMod -------------------------------------------------------------

  DO 831 jk=ktdia,klev
!
!     ------------------------------------------------------------------
!
!       2.    Set to zero local tendencies (increments)
! pointer on the rain and snowflux through the bottom of each layer
    zrfl => pfrain(:,jk)
    zsfl => pfsnow(:,jk)
    if (jk > 1 ) then
      zrfl(1:kproma) = pfrain(1:kproma, jk-1)
      zsfl(1:kproma) = pfsnow(1:kproma, jk-1)
    ENDIF
! pointer on the evaporation of rain
    zevp => prevap(:,jk)
! pointer on the sublimation of snow
    zsub => pssubl(:,jk)
! pointer on the rain production
    zrpr => prate_r(:,jk)
!
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
!--- Included for prognostic CDNC/IC scheme ----------------------------
!    additional arrays for transformation rate changes of CDNC and ICNC
     zrprn(:)      = 0.0_dp
     zfrln(:)      = 0.0_dp
     zsacln(:)     = 0.0_dp
     zsprn(:,jk)   = 0.0_dp
  
!--- CCMod -------------------------------------------------------------

     zximlt_cc(:) = 0._dp
     zimlt_cc(:)  = 0._dp
     zxisub_cc(:) = 0._dp
     zxievap_cc(:) = 0._dp
     zdepcc(:)    = 0._dp
     zdepco(:)    = 0._dp
     zspr_cc(:)   = 0.0_dp
     zsprn_cc(:,jk)   = 0.0_dp
     zcorclc(:) = 0._dp

!--- End CCMod -------------------------------------------------------------

     ztmp(:,jk)    = MAX(pqm1(:,jk),0.0_dp)
     ztmp(:,jk)    = 1._dp/(cpd+zcons1*ztmp(:,jk))
     zlvdcp(:)     = alv*ztmp(:,jk)
     zlsdcp(:)     = als*ztmp(:,jk)
!
!     ------------------------------------------------------------------
!
!       3.   Modification of incoming precipitation fluxes by
!            melting, sublimation and evaporation
!
     IF (jk > 1) THEN
!
!       3.1   Melting of snow and ice
!
        ztmp1_1d(:) = ptm1(:,jk) - tmelt
       
        ll1_1d(:) = (ztmp1_1d(:) > 0._dp)
 
        ztmp1_1d(:) = MAX(0._dp, ztmp1_1d(:))   !SF ztdif
        ztmp1_1d(:) = zcons2*ztmp1_1d(:)*zdp_2d(:,jk) / ( zlsdcp(:)-zlvdcp(:) ) !SF zcons*ztdif
        
        ztmp2_1d(:) = zxsec*zsfl(:)
        ztmp2_1d(:) = MIN(ztmp2_1d(:), ztmp1_1d(:)) !SFzximelt
        pimelt(:,jk) = ztmp2_1d(:)

        zrfl(:)  = zrfl(:) + ztmp2_1d(:)
        zsfl(:)  = zsfl(:) - ztmp2_1d(:)


        zsmlt(:) = ztmst*g*ztmp2_1d(:) / zdp_2d(:,jk)

        ztmp2_1d(:) = zxsec*zxiflux(:)
        ztmp2_1d(:) = MIN(ztmp2_1d(:), ztmp1_1d(:))

        ll2_1d(:) = (zxiflux(:) > zepsec)

        ztmp3_1d(:) = zxifluxn(:)*ztmp2_1d(:)/MAX(zxiflux(:),zepsec)
        ztmp3_1d(:) = MERGE(ztmp3_1d(:), 0._dp, ll2_1d(:)) !SF zxinmelt

        zxiflux(:)  = zxiflux(:)  - ztmp2_1d(:)
        zxifluxn(:) = zxifluxn(:) - ztmp3_1d(:)

        ll2_1d(:) = (zxiflux(:) < zepsec)

        zxifluxn(:) = MERGE(0._dp, zxifluxn(:), ll2_1d(:))
 
        zximlt(:)   = ztmst*g*ztmp2_1d(:) / zdp_2d(:,jk)

!--- CCMod -------------------------------------------------------------

        !melting of sedimentating contrail ice

        ztmp2_1d(:) = zxsec*zxiflux_cc(:)
        ztmp2_1d(:) = MIN(ztmp2_1d(:), ztmp1_1d(:))

        ll2_1d(:) = (zxiflux_cc(:) > zepsec)

        ztmp3_1d(:) = zxifluxn_cc(:)*ztmp2_1d(:)/MAX(zxiflux_cc(:),zepsec)
        ztmp3_1d(:) = MERGE(ztmp3_1d(:), 0._dp, ll2_1d(:)) !SF zxinmelt

        zxiflux_cc(:)  = zxiflux_cc(:)  - ztmp2_1d(:)
        zxifluxn_cc(:) = zxifluxn_cc(:) - ztmp3_1d(:)

        ll2_1d(:) = (zxiflux_cc(:) < zeps)

        zxifluxn_cc(:) = MERGE(0._dp, zxifluxn_cc(:), ll2_1d(:))

        zximlt_cc(:)   = ztmst*g*ztmp2_1d(:) / zdp_2d(:,jk)

!--- End CCMod -------------------------------------------------------------

        ztmp1_1d(:) = pxim1(:,jk) + ztmst*pxite(:,jk)
        ztmp1_1d(:) = MAX(0._dp, ztmp1_1d(:))
        zimlt(:)    = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))

!--- Included for prognostic CDNC/IC scheme (Philip Stier, 31/03/2004) -
!    If T > tmelt melt all ice crystals and transfer to cloud droplets 

        ztmp1_1d(:) = MERGE(zicncq(:,jk), 0._dp, ll1_1d(:))
        zicnc(:,jk) = MERGE(zicemin, zicnc(:,jk), ll1_1d(:)) 

!--- CCMod -------------------------------------------------------------

        !melting of contrail ice
        zimlt_cc(:) = MERGE(zcciwc(:,jk), 0._dp, ll1_1d(:))

        ztmp2_1d(:) = MERGE(zccicnc(:,jk), 0._dp, ll1_1d(:))
        zccicnc(:,jk) = MERGE(0._dp, zccicnc(:,jk), ll1_1d(:)) 
        zccvol_ges(:,jk) = MERGE(0._dp, zccvol_ges(:,jk), ll1_1d(:))
        zcccov_ges(:,jk) = MERGE(0._dp, zcccov_ges(:,jk), ll1_1d(:))

!--- End CCMod -------------------------------------------------------------

!--- CCMod -------------------------------------------------------------

        zcdnc(:,jk) = zcdnc(:,jk) + ztmp1_1d(:) + ztmp2_1d(:)
        qmel(1:kproma,jk,jrow) = qmel(1:kproma,jk,jrow) + zdt*(ztmp1_1d(:)+ztmp2_1d(:))

!--- End included for CDNC/IC scheme -----------------------------------
!
!       3.2   Sublimation of snow (zsub) and ice (zxisub) following (Lin et al., 1983)
!
        ll1_1d(:) = (zclcpre(:) > 0._dp)

        ztmp1_1d(:) = 1._dp/(2.43e-2_dp*rv)*zlsdcp(:)**2/ptm1(:,jk)**2
        ztmp1_1d(:) = ztmp1_1d(:) &
                    + 1._dp/0.211e-4_dp*zrho_rcp(:,jk)/zqsi_2d(:,jk)
        ztmp1_1d(:) = 3.e6_dp*2._dp*api*zicesub(:,jk)*zrho_rcp(:,jk)/ztmp1_1d(:) !zcoeff

        ztmp2_1d(:) = MERGE(zclcpre(:), 1._dp, ll1_1d(:))

!SF     Snow:  
        ll2_1d(:) = (zsfl(:) > cqtmin) .AND. &
                    ll1_1d(:)

        ztmp3_1d(:) = zcons3*(zsfl(:) / ztmp2_1d(:))**(0.25_dp/1.16_dp)  !zclambs
        ztmp3_1d(:) = 0.78_dp  *ztmp3_1d(:)**2                                & 
                    + 232.19_dp*zqrho_2d(:,jk)**0.25_dp * ztmp3_1d(:)**2.625_dp     !zcfac4c
        ztmp3_1d(:) = ztmp3_1d(:) * ztmp1_1d(:) * zdpg_2d(:,jk)                     

        ztmp4_1d(:) = -zxsec * zsfl(:) / ztmp2_1d(:)
        ztmp4_1d(:) = MAX(ztmp4_1d(:), ztmp3_1d(:))                              !zzeps
        ztmp4_1d(:) = -ztmst*ztmp4_1d(:) / zdpg_2d(:,jk) * ztmp2_1d(:)

        ztmp5_1d(:) = zxsec*(zqsi_2d(:,jk)-pqm1(:,jk))
        ztmp5_1d(:) = MAX(ztmp5_1d(:),0._dp)

        ztmp4_1d(:) = MIN(ztmp4_1d(:),ztmp5_1d(:))
        ztmp4_1d(:) = MAX(ztmp4_1d(:),0._dp)
        zsub(:)     = MERGE(ztmp4_1d(:), zsub(:), ll2_1d(:))

!SF     Ice:
        ll2_1d(:) = (zxiflux(:) > cqtmin) .AND. &
                    ll1_1d(:)

        ztmp3_1d(:) = zcons3*(zxiflux(:) / ztmp2_1d(:))**(0.25_dp/1.16_dp) !zclambs
        ztmp3_1d(:) = 0.78_dp  *ztmp3_1d(:)**2                                &
                    + 232.19_dp*zqrho_2d(:,jk)**0.25_dp * ztmp3_1d(:)**2.625_dp     !zcfac4c
        ztmp3_1d(:) = ztmp3_1d(:) * ztmp1_1d(:) * zdpg_2d(:,jk)

        ztmp4_1d(:) = -zxsec * zxiflux(:) / ztmp2_1d(:)
        ztmp4_1d(:) = MAX(ztmp4_1d(:), ztmp3_1d(:))                              !zzeps
        ztmp4_1d(:) = -ztmst*ztmp4_1d(:) / zdpg_2d(:,jk) * ztmp2_1d(:)

        ztmp5_1d(:) = zxsec*(zqsi_2d(:,jk)-pqm1(:,jk))
        ztmp5_1d(:) = MAX(ztmp5_1d(:),0._dp)

        ztmp4_1d(:) = MIN(ztmp4_1d(:),ztmp5_1d(:))
        ztmp4_1d(:) = MAX(ztmp4_1d(:),0._dp)
        zxisub(:)   = MERGE(ztmp4_1d(:), zxisub(:), ll2_1d(:))

        ztmp5_1d(:) = zxisub(:) * zxifluxn(:) / MAX(zxiflux(:),cqtmin)  !SF zsubin
        ztmp5_1d(:) = zcons2 * ztmp5_1d(:) * zdp_2d(:,jk)
        ztmp5_1d(:) = MERGE(ztmp5_1d(:), 0._dp, ll2_1d(:))
        zxifluxn(:) = zxifluxn(:) - ztmp5_1d(:)

        zxiflux(:) = zxiflux(:) - zcons2*zxisub(:)*zdp_2d(:,jk)

        ll3_1d(:) = ll2_1d(:) .AND. (zxiflux(:) < zepsec)

        zxifluxn(:) = MERGE(0._dp, zxifluxn(:), ll3_1d(:))

!--- CCMod -------------------------------------------------------------

        !Sublimation of sed. Contrail Ice:
        ll2_1d(:) = (zxiflux_cc(:) > cqtmin) .AND. ll1_1d(:)

        ztmp3_1d(:) = zcons3*(zxiflux_cc(:) / ztmp2_1d(:))**(0.25_dp/1.16_dp) !zclambs
        ztmp3_1d(:) = 0.78_dp  *ztmp3_1d(:)**2                                &
                    + 232.19_dp*zqrho_2d(:,jk)**0.25_dp * ztmp3_1d(:)**2.625_dp     !zcfac4c
        ztmp3_1d(:) = ztmp3_1d(:) * ztmp1_1d(:) * zdpg_2d(:,jk)

        ztmp4_1d(:) = -zxsec * zxiflux_cc(:) / ztmp2_1d(:)
        ztmp4_1d(:) = MAX(ztmp4_1d(:), ztmp3_1d(:))                              !zzeps
        ztmp4_1d(:) = -ztmst*ztmp4_1d(:) / zdpg_2d(:,jk) * ztmp2_1d(:)

        ztmp5_1d(:) = zxsec*(zqsi_2d(:,jk)-pqm1(:,jk))
        ztmp5_1d(:) = MAX(ztmp5_1d(:),0._dp)

        ztmp4_1d(:) = MIN(ztmp4_1d(:),ztmp5_1d(:))
        ztmp4_1d(:) = MAX(ztmp4_1d(:),0._dp)
        zxisub_cc(:)   = MERGE(ztmp4_1d(:), zxisub_cc(:), ll2_1d(:))

        ztmp5_1d(:) = zxisub_cc(:) * zxifluxn_cc(:) / MAX(zxiflux_cc(:),cqtmin)  !SF zsubin
        ztmp5_1d(:) = zcons2 * ztmp5_1d(:) * zdp_2d(:,jk)
        ztmp5_1d(:) = MERGE(ztmp5_1d(:), 0._dp, ll2_1d(:))
        zxifluxn_cc(:) = zxifluxn_cc(:) - ztmp5_1d(:)

        zxiflux_cc(:) = zxiflux_cc(:) - zcons2*zxisub_cc(:)*zdp_2d(:,jk)

        ll3_1d(:) = ll2_1d(:) .AND. (zxiflux_cc(:) < zepsec)

        zxifluxn_cc(:) = MERGE(0._dp, zxifluxn_cc(:), ll3_1d(:))

!--- End CCMod -------------------------------------------------------------

!
!       3.3   Evaporation of rain (zevp following Rotstayn, 1997)
!
        ll2_1d(:) = (zrfl(:) > cqtmin) .AND. &
                    ll1_1d(:)

        ztmp3_1d(:) = 870._dp * zsusatw_evap(:,jk) * zdpg_2d(:,jk)      &
                              * (zrfl(:)/ztmp2_1d(:))**0.61_dp     &
                              / (SQRT(zrho(:,jk))*zastbstw(:,jk))                 

        ztmp4_1d(:) = -zxsec*zrfl(:)/ztmp2_1d(:)
        ztmp4_1d(:) = MAX(ztmp4_1d(:), ztmp3_1d(:))
        ztmp4_1d(:) = -ztmst*ztmp4_1d(:)*ztmp2_1d(:)/zdpg_2d(:,jk)

        ztmp5_1d(:) = zxsec*(zqsw_2d(:,jk)-pqm1(:,jk))
        ztmp5_1d(:) = MAX(ztmp5_1d(:), 0._dp)

        ztmp4_1d(:) = MIN(ztmp4_1d(:),ztmp5_1d(:))
        ztmp4_1d(:) = MAX(ztmp4_1d(:),0._dp)
        zevp(:)     = MERGE(ztmp4_1d(:), zevp(:), ll2_1d(:))

        IF (lookupoverflow) THEN
          status_string = 'lookuperror: cdnc - cloud (1)'
          RETURN
        ENDIF
      ENDIF !SF end jk > 1
!
!     --- End included for CDNC/IC scheme -----------------------------------
!
!     ------------------------------------------------------------------
!       4.    Sedimentation of cloud ice from grid-mean values (zicesed).
!             Updating the tendency 'pxite' to include sedimentation.
!             At jk=klev, the sedimentation sink is balanced by
!             precipitation at the surface (through 'zzdrs', see 7.3).
!             Finally: In-cloud cloud water/ice.
!
        zxip1_1d(:) = pxim1(:,jk) + ztmst*pxite(:,jk)-zimlt(:)
        zxip1_1d(:) = MAX(zxip1_1d(:), EPSILON(1._dp))

        zicncp1_1d(:) = zicnc(:,jk)*paclc(:,jk)
        zicncp1_1d(:) = MAX(zicncp1_1d(:),zicemin)

        zmmean_2d(:,jk) = zrho(:,jk)*zxip1_1d(:)/zicncp1_1d(:)
        zmmean_2d(:,jk) = MAX(zmmean_2d(:,jk), zmi)

        ll1_1d(:) = (zmmean_2d(:,jk) < 2.166E-9_dp )
        ll2_1d(:) = (.NOT. ll1_1d(:)) .AND. (zmmean_2d(:,jk) < 4.264E-8_dp )

        zalfased_2d(:,jk) = MERGE(63292.4_dp, 8.78_dp, ll1_1d(:))
        zalfased_2d(:,jk) = MERGE(329.75_dp, zalfased_2d(:,jk), ll2_1d(:))

        zbetased_2d(:,jk) = MERGE(0.5727_dp, 0.0954_dp, ll1_1d(:))
        zbetased_2d(:,jk) = MERGE(0.3091_dp, zbetased_2d(:,jk), ll2_1d(:))

        zxifallmc_2d(:,jk) = zfall*zalfased_2d(:,jk) & 
                           *(zmmean_2d(:,jk)**zbetased_2d(:,jk))*zaaa_2d(:,jk) !Fallgeschwindigkeit Masse

        !limit fall velocity to 1 cm/s - 2 m/s:
        zxifallmc_2d(:,jk) = MAX(0.001_dp,zxifallmc_2d(:,jk))
        zxifallmc_2d(:,jk) = MIN(2._dp,zxifallmc_2d(:,jk))

        zxifallnc_2d(:,jk) = zxifallmc_2d(:,jk)                                !Fallgeschwindigkeit Anzahl

        ztmp1_1d(:) = ztmst*g*zxifallmc_2d(:,jk)*zrho(:,jk)/zdp_2d(:,jk) !SF zal1

        ll1_1d(:) = (zxifallmc_2d(:,jk) > zeps)

        ztmp2_1d(:) = zxiflux(:)*zrho_rcp(:,jk)/MAX(zxifallmc_2d(:,jk), zeps)  
        ztmp2_1d(:) = MERGE(ztmp2_1d(:), 0._dp, ll1_1d(:)) !SF zal2
        
        ztmp3_1d(:) = g*ztmst*zxifallnc_2d(:,jk)*zrho(:,jk)/zdp_2d(:,jk) !SF zal3

        ll1_1d(:)   = (zxifallnc_2d(:,jk) > zeps)
        ztmp4_1d(:) = zxifluxn(:)/MAX(zxifallnc_2d(:,jk), zeps)
        ztmp4_1d(:) = MERGE(ztmp4_1d(:), 0._dp, ll1_1d(:)) !SF zal4

!        -----------change steffi-------------------
        zxised_2d(:,jk)  =  zxip1_1d(:)  *EXP(-ztmp1_1d(:))+ztmp2_1d(:)*(1._dp-EXP(-ztmp1_1d(:)))
        zicesed_2d(:,jk) =  zicncp1_1d(:)*EXP(-ztmp3_1d(:))+ztmp4_1d(:)*(1._dp-EXP(-ztmp3_1d(:)))
!        --------end change steffi--------------------

        ll1_1d(:) = (paclc(:,jk) > 0.01_dp)

        ztmp1_1d(:) = zicesed_2d(:,jk) / MAX(paclc(:,jk), 0.01_dp)
        ztmp1_1d(:) = MERGE(ztmp1_1d(:), zicesed_2d(:,jk), ll1_1d(:))

        zicnc(:,jk) = ztmp1_1d(:) + znidetr(:,jk) + zninucl(:,jk)
        zicnc(:,jk) = MIN(zicnc(:,jk), zicemax)
        zicnc(:,jk) = MAX(zicnc(:,jk), zicemin)

        zxiflux(:)  = zxiflux(:)  + zcons2*(zxip1_1d(:)   - zxised_2d(:,jk) )*zdp_2d(:,jk)
        zxifluxn(:) = zxifluxn(:) + zcons2*(zicncp1_1d(:) - zicesed_2d(:,jk))*zdp_2d(:,jk)*zrho_rcp(:,jk)

        pxite(:,jk) = ztmst_rcp*(zxised_2d(:,jk)-pxim1(:,jk))

!--- CCMod -------------------------------------------------------------

        !Sedimentation von Contrail-Eis

        IF(jk>1) THEN

        zcciwc(:,jk) = zcciwc(:,jk)-zimlt_cc(:)
        zcciwc(:,jk) = MAX(zcciwc(:,jk), 0._dp)

        zicncp1_cc(:) = zccicnc(:,jk)*zccvol_ges(:,jk)
        zicncp1_cc(:) = MAX(zicncp1_cc(:),zeps)

        !mass of particles with radius of 0.1mu
        zmi_cc = 4._dp/3._dp*1.e-7_dp**3*api*crhoi
        zmmean_2d(:,jk) = zrho(:,jk)*zcciwc(:,jk)/zicncp1_cc(:)
        zmmean_2d(:,jk) = MAX(zmmean_2d(:,jk), zmi_cc)

        ll1_1d(:) = (zmmean_2d(:,jk) < 2.166E-9_dp )
        ll2_1d(:) = (.NOT. ll1_1d(:)) .AND. (zmmean_2d(:,jk) < 4.264E-8_dp )

        zalfased_2d(:,jk) = MERGE(63292.4_dp, 8.78_dp, ll1_1d(:))
        zalfased_2d(:,jk) = MERGE(329.75_dp, zalfased_2d(:,jk), ll2_1d(:))

        zbetased_2d(:,jk) = MERGE(0.5727_dp, 0.0954_dp, ll1_1d(:))
        zbetased_2d(:,jk) = MERGE(0.3091_dp, zbetased_2d(:,jk), ll2_1d(:))

        zxifallmc_2d(:,jk) = zfall*zalfased_2d(:,jk) &
                           *(zmmean_2d(:,jk)**zbetased_2d(:,jk))*zaaa_2d(:,jk) !Fallgeschwindigkeit Masse 
   
        !limit fall velocity to 1 cm/s - 2 m/s:
        zxifallmc_2d(:,jk) = MAX(0.001_dp,zxifallmc_2d(:,jk))
        zxifallmc_2d(:,jk) = MIN(2._dp,zxifallmc_2d(:,jk))

        zxifallnc_2d(:,jk) = zxifallmc_2d(:,jk)                                !Fallgeschwindigkeit Anzahl

        ztmp1_1d(:) = ztmst*g*zxifallmc_2d(:,jk)*zrho(:,jk)/zdp_2d(:,jk) !SF zal1

        ll1_1d(:) = (zxifallmc_2d(:,jk) > zeps)

        ztmp2_1d(:) = zxiflux_cc(:)*zrho_rcp(:,jk)/MAX(zxifallmc_2d(:,jk), zeps)
        ztmp2_1d(:) = MERGE(ztmp2_1d(:), 0._dp, ll1_1d(:)) !SF zal2

        ztmp3_1d(:) = g*ztmst*zxifallnc_2d(:,jk)*zrho(:,jk)/zdp_2d(:,jk) !SF zal3

        ll1_1d(:)   = (zxifallnc_2d(:,jk) > zeps)
        ztmp4_1d(:) = zxifluxn_cc(:)/MAX(zxifallnc_2d(:,jk), zeps)
        ztmp4_1d(:) = MERGE(ztmp4_1d(:), 0._dp, ll1_1d(:)) !SF zal4

        zxised_cc(:,jk)  =  zcciwc(:,jk)  *EXP(-ztmp1_1d(:))
        zicesed_cc(:,jk) =  zicncp1_cc(:)*EXP(-ztmp3_1d(:))

        ztmp1_1d(:) = ztmst*g*zxifallmc_2d(:,jk-1)*zrho(:,jk)/zdp_2d(:,jk) !SF zal1
        ll1_1d(:) = (zxifallmc_2d(:,jk-1) > zeps)
        ztmp2_1d(:) = zxiflux_cc(:)*zrho_rcp(:,jk)/MAX(zxifallmc_2d(:,jk-1), zeps)
        ztmp2_1d(:) = MERGE(ztmp2_1d(:), 0._dp, ll1_1d(:)) !SF zal2

        ztmp3_1d(:) = g*ztmst*zxifallnc_2d(:,jk-1)*zrho(:,jk)/zdp_2d(:,jk) !SF zal3
        ll1_1d(:)   = (zxifallnc_2d(:,jk-1) > zeps)
        ztmp4_1d(:) = zxifluxn_cc(:)/MAX(zxifallnc_2d(:,jk-1), zeps)
        ztmp4_1d(:) = MERGE(ztmp4_1d(:), 0._dp, ll1_1d(:)) !SF zal4

        zxisedf_cc(:,jk)  =  ztmp2_1d(:)*(1._dp-EXP(-ztmp1_1d(:)))
        zicesedf_cc(:,jk) = ztmp4_1d(:)*(1._dp-EXP(-ztmp3_1d(:)))

        zxiflux_cc(:)  = zxiflux_cc(:) + zcons2*(zcciwc(:,jk) &
                           - (zxised_cc(:,jk)+zxisedf_cc(:,jk)))*zdp_2d(:,jk)
        zxifluxn_cc(:) = zxifluxn_cc(:) + zcons2*(zicncp1_cc(:) &
                           - (zicesed_cc(:,jk)+zicesedf_cc(:,jk)))*zdp_2d(:,jk)*zrho_rcp(:,jk)

        zccxite(:) = ztmst_rcp*((zxised_cc(:,jk)+zxisedf_cc(:,jk))-pxtm1(:,jk,25))

        END IF !(jk > 1)

!--- End CCMod -------------------------------------------------------------

!---Included for in-cloud scavenging (Philip Stier, 16/04/04):----------
        zmrateps(:,jk) = zmrateps(:,jk) + zxip1_1d(:) - zxised_2d(:,jk)
!---End Included for scavenging-----------------------------------------
        pisedi(:,jk) = zxised_2d(:,jk)

!
!             In-cloud water/ice calculated from respective grid-means,
!             partial cloud cover, advective/diffusive tendencies,
!             detrained cloud water/ice and ice sedimentation.
!             In-cloud values are required for cloud microphysics.
!
!       THIS PART NEEDS TO BE COMMENTED BY ERICH/ADRIAN
!
        zclcaux(:) = paclc(:,jk)
        locc_1d(:) = (zclcaux(:) > zeps)

        ztmp1_1d(:) = 1._dp/MAX(pqm1(:,jk),zeps) + zlsdcp(:)*alv/(rv*ptm1(:,jk)**2) !SF protected against divide by zero 
        ztmp2_1d(:) = g*(zlvdcp(:)*rd/rv/ptm1(:,jk)-1._dp)/(rd*ptm1(:,jk))
        zkair_1d(:) = 4.1867e-3_dp*(5.69_dp + 0.017_dp*(ptm1(:,jk)-tmelt)) ! eq. 13-18a P&K !SF replaced ztp1tmp 
                                                                           ! by ptm1
        zdv_1d(:)   = 2.21_dp/papm1(:,jk)
        ztmp3_1d(:) = 1._dp/(crhoi*als**2/(zkair_1d*rv*ptm1(:,jk)**2) & 
                    + crhoi*rv*ptm1(:,jk)/(zesi_2d(:,jk)*zdv_1d(:)))

        zrice_1d(:)    = (0.75_dp*zxised_2d(:,jk)/(api*zicnc(:,jk)*crhoi))**(1._dp/3._dp)
        zeta_1d(:)     = ztmp1_1d(:)/ztmp2_1d(:)*ztmp3_1d(:)*4._dp*api*crhoi*zcap/zrho(:,jk)
        zvervmax_1d(:) = (zesw_2d(:,jk)-zesi_2d(:,jk))/zesi_2d(:,jk)*zicnc(:,jk) &
                                                      *zrice_1d(:)*zeta_1d(:)

        lo2_1d(:)  = (ptm1(:,jk) < cthomi) .OR.                   &
                      (ptm1(:,jk) < tmelt .AND. 0.01_dp*zvervx_2d(:,jk) < zvervmax_1d(:))


        ll1_1d(:)  =  (ptm1(:,jk) < tmelt .AND. ptm1(:,jk) > cthomi .AND. & 
             0.01_dp*zvervx_2d(:,jk) < zvervmax_1d(:) .AND. zclcaux(:) > 0._dp)

        IF(ncvmicro>0) THEN
            zxite(:) = pxteci(:,jk)
            zxlte(:) = pxtecl(:,jk)
        ELSE
            zxite(:) = MERGE(pxtec(:,jk), 0._dp      , lo2_1d(:))
            zxlte(:) = MERGE(0._dp      , pxtec(:,jk), lo2_1d(:))
        ENDIF

        ztmp1_1d(:) = pxim1(:,jk)/MAX(zclcaux(:),zeps) !SF to protect from a division by 0
        ztmp2_1d(:) = pxlm1(:,jk)/MAX(zclcaux(:),zeps) 
    
        zxib(:) = MERGE(ztmp1_1d(:), 0._dp, locc_1d(:))
        zxlb(:) = MERGE(ztmp2_1d(:), 0._dp, locc_1d(:))

        zxim1evp_1d(:) = MERGE(0._dp, pxim1(:,jk), locc_1d(:))
        zxlm1evp_1d(:) = MERGE(0._dp, pxlm1(:,jk), locc_1d(:)) 

        zxite2(:) = MERGE(zxite(:), 0._dp   , lo2_1d(:)) !SF temporary var needed if ncvmicro>0
        zxlte2(:) = MERGE(0._dp   , zxlte(:), lo2_1d(:)) !SF temporary var needed if ncvmicro>0

        zxidt_1d(:) = ztmst*(pxite(:,jk)+zxite2(:))
        zxldt_1d(:) = ztmst*(pxlte(:,jk)+zxlte2(:)) + zximlt(:) + zimlt(:) 

        !SF ice cloud:

        ll1_1d(:) = (zxidt_1d(:) > 0._dp)
        zxidtstar_1d(:) = MERGE(zxidt_1d(:), 0._dp, ll1_1d(:))
       
        ztmp1_1d(:) = zxidt_1d(:)/MAX(zclcaux(:), zeps)
        ztmp1_1d(:) = MAX(ztmp1_1d(:), -zxib(:))
        ztmp1_1d(:) = MERGE(zxidt_1d(:), ztmp1_1d(:), ll1_1d(:))
        ztmp1_1d(:) = MERGE(ztmp1_1d(:), 0._dp, locc_1d(:))
        zxib(:) = zxib(:) + ztmp1_1d(:)


!--- CCMod -------------------------------------------------------------
        !Sedimentation&Transport:(1-b_ci-b_cc) verdampft; Detrainment:(1-b_ci) verdampft
        ztmp1_1d(:) = zccvol_ges(:,jk)/MAX(zclcaux(:), zeps)*ztmst*MAX(pxite(:,jk),0._dp)
        zxib(:) = MERGE(zxib(:) + ztmp1_1d(:), 0._dp, zclcaux(:).gt.zeps)
        zxidtstar_1d(:) = (1._dp-zclcaux(:))*zxidtstar_1d(:)+zccvol_ges(:,jk)*MAX(pxite(:,jk),0._dp)
!--- End CCMod -------------------------------------------------------------

        ztmp1_1d(:) = -( ztmst_rcp*pxim1(:,jk) + zxite2(:) )
        ztmp1_1d(:) = MAX(pxite(:,jk), ztmp1_1d(:))
        pxite(:,jk) = MERGE(pxite(:,jk), ztmp1_1d(:), ll1_1d(:))

        ll2_1d(:) = (.NOT. locc_1d(:)) .AND. (.NOT. ll1_1d(:))
        ztmp1_1d(:) = ztmst * ( pxite(:,jk) + zxite2(:) )
        ztmp1_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll2_1d(:))
        zxim1evp_1d(:) = zxim1evp_1d(:) + ztmp1_1d(:)

        !SF water cloud:

        ll1_1d(:) = (zxldt_1d(:) > 0._dp)
        zxldtstar_1d(:) = MERGE(zxldt_1d(:), 0._dp, ll1_1d(:))

        ztmp1_1d(:) = zxldt_1d(:)/MAX(zclcaux(:), zeps)
        ztmp1_1d(:) = MAX(ztmp1_1d(:), -zxlb(:))
        ztmp1_1d(:) = MERGE(zxldt_1d(:), ztmp1_1d(:), ll1_1d(:))
        ztmp1_1d(:) = MERGE(ztmp1_1d(:), 0._dp, locc_1d(:))
        zxlb(:)     = zxlb(:) + ztmp1_1d(:)

        ztmp1_1d(:) = -(pxlm1(:,jk)/ztmst+zxlte2(:))
        ztmp1_1d(:) = MAX(pxlte(:,jk), ztmp1_1d(:))
        pxlte(:,jk) = MERGE(pxlte(:,jk), ztmp1_1d(:), ll1_1d(:))

        ll2_1d(:)      = (.NOT. locc_1d(:)) .AND. (.NOT. ll1_1d(:))
        ztmp1_1d(:)    = ztmst * ( pxlte(:,jk) + zxlte2(:) )
        ztmp1_1d(:)    = MERGE(ztmp1_1d(:), 0._dp, ll2_1d(:))
        zxlm1evp_1d(:) = zxlm1evp_1d(:) + ztmp1_1d(:)

!--- CCMod -------------------------------------------------------------

        IF(jk>1) THEN

        ll_bcc(:) = zccvol_ges(:,jk) .gt. 1.e-10_dp 


        where(ll_bcc(:))
        !neues Volumen durch Sedimentation innerhalb der Gitterbox
          ztmp1_1d(:) = zdz_2d(:,jk)*zccvol_ges(:,jk)/MAX(zcccov_ges(:,jk),zeps)
          ztmp1_1d(:) = MIN(MAX(ztmp1_1d(:),0._dp),zdz_2d(:,jk)) !h
          ztmp2_1d(:) = MIN(zxifallnc_2d(:,jk)*ztmst, zdz_2d(:,jk)-ztmp1_1d(:))
          ztmp2_1d(:) = MAX(ztmp2_1d(:),0._dp) !h*
          ztmp4_1d(:) = 0.5_dp*ztmp2_1d(:)*zcccov_ges(:,jk) &
                         *(2._dp - ztmp2_1d(:)/MAX(zdz_2d(:,jk)-ztmp1_1d(:),zeps)) !A1
          ztmp4_1d(:) = MERGE(ztmp4_1d(:),0._dp,ztmp1_1d(:).lt.zdz_2d(:,jk))
          zb_stern(:,jk) = zcccov_ges(:,jk)*ztmp2_1d(:)/MAX(zdz_2d(:,jk)-ztmp1_1d(:),zeps)
          zb_stern(:,jk) = MERGE(zb_stern(:,jk),zcccov_ges(:,jk),ztmp1_1d(:).lt.zdz_2d(:,jk)) !b*
          za_stern(:,jk) = zxifallnc_2d(:,jk)*ztmst*zcccov_ges(:,jk) - ztmp4_1d(:)
          za_stern(:,jk) = MAX(za_stern(:,jk),0._dp) !A*
          zccvol_ges(:,jk) = zccvol_ges(:,jk) + MAX(ztmp4_1d(:)/zdz_2d(:,jk),0._dp)
          zccvol_ges(:,jk) = MIN(zccvol_ges(:,jk),zcccov_ges(:,jk))

        !mehr Volumenzuwachs von KS, die sich ober oder unterhalb gebildet haben
          ztmp3_1d(:) = 0._dp
          WHERE(zcclen_ges(:,jk).gt.zeps .AND. zcclen_ou(:,jk).gt.0._dp)
            ztmp2_1d(:) = 1._dp-MIN(1._dp,MAX(zccvol_ges(:,jk)/zcccov_ges(:,jk),0._dp))
            ztmp1_1d(:) = ztmp4_1d(:)/zdz_2d(:,jk) * ztmp2_1d(:) 
            ztmp3_1d(:) = zcclen_ou(:,jk)/zcclen_ges(:,jk) * ztmp1_1d(:)
            zcclen_ou(:,jk) = zcclen_ou(:,jk) * ztmp2_1d(:)
            zccvol_ges(:,jk) = zccvol_ges(:,jk) + ztmp3_1d(:)
            zccvol_ges(:,jk) = MIN(zccvol_ges(:,jk),zcccov_ges(:,jk))
          ENDWHERE

        !Test, ob Zeitskala fuer Partikel groesser
          ztmp2_1d(:) = zicesed_cc(:,jk)/zccvol_ges(:,jk) !n
          ztmp3_1d(:) = (3._dp*zxised_cc(:,jk)*zrho(:,jk)/(4._dp*api*zrhoice*MAX(zicesed_cc(:,jk),zeps)))**0.333 !r
          ztmp4_1d(:) = 1.e-4_dp * 0.211_dp * (ptm1(:,jk)/273.15_dp)**1.94_dp * 1013.15_dp/(0.01_dp*papm1(:,jk)) !D_w
          !Koeffizienten fuer Zeitskala
          zfre_2d(:,jk) = 2._dp*zrho(:,jk)*ztmp3_1d(:)*zxifallnc_2d(:,jk)/zviscos_2d(:,jk)
          zfre_2d(:,jk) = 1._dp + 0.229_dp * sqrt(zfre_2d(:,jk))
          zcap_2d(:,jk) = 1.1_dp
          ztmp1_1d(:) = 1._dp / MAX(4._dp*api*zfre_2d(:,jk)*zcap_2d(:,jk)*ztmp2_1d(:)*ztmp3_1d(:)*ztmp4_1d(:),zeps)
          ztmp1_1d(:) = MERGE(zdtime/ztmp1_1d(:),1._dp,ztmp1_1d(:).gt.zdtime)
          ztmp1_1d(:) = MAX(MIN(ztmp1_1d(:),1._dp),0._dp)
          zccvol_ges2(:,jk) = ztmp1_1d(:) * zccvol_ges(:,jk)

        ENDWHERE

        !neuer Bedeckungsgrad wenn Sedimentation von oben
        ztmp1_1d(:) = 0._dp
        ztmp2_1d(:) = 0._dp
        ztmp3_1d(:) = 0._dp

        !Test, ob Zeitskala fuer Partikel groesser
        ztmp5_1d(:) = MIN(za_stern(:,jk-1),zb_stern(:,jk-1)*zdz_2d(:,jk))
        where(ztmp5_1d(:).gt.zeps)
          ztmp2_1d(:) = zicesedf_cc(:,jk)/ztmp5_1d(:) !n
          ztmp3_1d(:) = (3._dp*zxisedf_cc(:,jk)*zrho(:,jk)/(4._dp*api*zrhoice*MAX(zicesedf_cc(:,jk),zeps)))**0.333 !r
          ztmp4_1d(:) = 1.e-4_dp * 0.211_dp * (ptm1(:,jk)/273.15_dp)**1.94_dp * 1013.15_dp/(0.01_dp*papm1(:,jk)) !D_w
          !Koeffizienten fuer Zeitskala
          zfre_2d(:,jk) = 2._dp*zrho(:,jk)*ztmp3_1d(:)*zxifallnc_2d(:,jk-1)/zviscos_2d(:,jk)
          zfre_2d(:,jk) = 1._dp + 0.229_dp * sqrt(zfre_2d(:,jk))
          zcap_2d(:,jk) = 1.1_dp
          ztmp1_1d(:) = 1._dp / MAX(4._dp*api*zfre_2d(:,jk)*zcap_2d(:,jk)*ztmp2_1d(:)*ztmp3_1d(:)*ztmp4_1d(:),zeps)
          ztmp1_1d(:) = MERGE(zdtime/ztmp1_1d(:),1._dp,ztmp1_1d(:).gt.zdtime)
          ztmp1_1d(:) = MAX(MIN(ztmp1_1d(:),1._dp),0._dp)
          za_stern2(:,jk) = ztmp1_1d(:) * ztmp5_1d(:)

          where(zb_stern(:,jk-1).gt.zccvol_ges(:,jk) .AND. zb_stern(:,jk-1).gt.1.e-10_dp)

            ztmp3_1d(:) = B_cc(:,jk,jrow)/B_cc(:,jk-1,jrow)
            ztmp3_1d(:) = MAX(MIN(ztmp3_1d(:), 1._dp),0._dp)

            ztmp1_1d(:) = MAX(MIN(zccvol_ges(:,jk)/zb_stern(:,jk-1),1._dp),0._dp)
            ztmp1_1d(:) = ztmp5_1d(:)*(1._dp - ztmp1_1d(:))
            zccvol_ges(:,jk) = zccvol_ges(:,jk) + ztmp3_1d(:)*ztmp1_1d(:)/zdz_2d(:,jk)

            zccvol_ges2(:,jk) = zccvol_ges2(:,jk) + ztmp3_1d(:)*za_stern2(:,jk)/zdz_2d(:,jk)

            zcccov_ges(:,jk) = MAX(zcccov_ges(:,jk),zb_stern(:,jk-1)*ztmp3_1d(:))

            zccvol_ges(:,jk)  = MIN(zccvol_ges(:,jk),  zcccov_ges(:,jk))
            zccvol_ges2(:,jk) = MIN(zccvol_ges2(:,jk), zccvol_ges(:,jk))

            ztmp1_1d(:) = MAX(zcclen_ges(:,jk-1) - zcclen_ges(:,jk), 0._dp)
            zcclen_ges2(:,jk) = zcclen_ges(:,jk) + ztmp3_1d(:)* ztmp1_1d(:)
 
          end where

        end where

        ll_bcc(:) = zcccov_ges(:,jk) .gt. 1.e-10_dp

        !Begrenzung auf potentiellen Bedeckungsgrad
        zhelp(:,jk) = 0._dp
        where (B_cc(:,jk,jrow) > 0._dp)
          ztmp1_1d(:) = zcccov_ges(:,jk)+zcccov_new(:,jk)-B_cc(:,jk,jrow)
          where(ztmp1_1d(:).gt.0._dp)
            zhelp(:,jk) = ztmp1_1d(:)/MAX(zcccov_ges(:,jk),zeps)
            zcccor(:,jk) = zcccov_ges(:,jk)*zhelp(:,jk)
            zcccor_vol(:,jk) = zccvol_ges(:,jk)*zhelp(:,jk)
          elsewhere
            zcccor(:,jk) = 0._dp
            zcccor_vol(:,jk) = 0._dp
          end where
        elsewhere
          zcccor(:,jk) = zcccov_ges(:,jk)
          zcccor_vol(:,jk) = zccvol_ges(:,jk)
        end where

        zcccov_ges(:,jk) = zcccov_ges(:,jk) - zcccor(:,jk)
        zccvol_ges(:,jk) = zccvol_ges(:,jk) - zcccor_vol(:,jk)

        !Verteilung der Raten auf die verschiedenen contrail Klassen
        !gescherter Bedeckungsgrad
        where (zcccov_ges(:,jk) .gt. 0._dp)
          zmult2(:,jk) = zcccov_ges(:,jk) - zcccov_ori(:,jk)
        elsewhere
          zmult2(:,jk) = 0._dp
        end where
        where (zcclen_ges(:,jk) .gt. zeps)
          zlmit(:,jk) = 1._dp/zcclen_ges(:,jk)
          zcccov_1dt(:,jk) = zcccov_1dt(:,jk)+zmult2(:,jk)*zcclen_1dt(:,jk)*zlmit(:,jk)
          zcccov_2dt(:,jk) = zcccov_2dt(:,jk)+zmult2(:,jk)*zcclen_2dt(:,jk)*zlmit(:,jk)
          zcccov_3dt(:,jk) = zcccov_3dt(:,jk)+zmult2(:,jk)*zcclen_3dt(:,jk)*zlmit(:,jk)
          zcccov_4dt(:,jk) = zcccov_4dt(:,jk)+zmult2(:,jk)*zcclen_4dt(:,jk)*zlmit(:,jk)
          zcccov_5dt(:,jk) = zcccov_5dt(:,jk)+zmult2(:,jk)*zcclen_5dt(:,jk)*zlmit(:,jk)
          zcccov(:,jk)     = zcccov(:,jk)    +zmult2(:,jk)*zcclen(:,jk)    *zlmit(:,jk)
        elsewhere
          zcccov_1dt(:,jk) = 0._dp
          zcccov_2dt(:,jk) = 0._dp
          zcccov_3dt(:,jk) = 0._dp
          zcccov_4dt(:,jk) = 0._dp
          zcccov_5dt(:,jk) = 0._dp
          zcccov(:,jk)     = zcccov(:,jk)
        end where
        !Volumenbedeckungsgrad
        where (zcccov_ges(:,jk) .gt. 0._dp)
          zmult2(:,jk) = zccvol_ges(:,jk) - zccvol_ges_ori(:,jk)
        elsewhere
          zmult2(:,jk) = 0._dp
        end where
        where (zcclen_ges(:,jk) .gt. zeps)
          zlmit(:,jk) = 1._dp/zcclen_ges(:,jk)
          zccvol_1dt(:,jk) = zccvol_1dt(:,jk)+zmult2(:,jk)*zcclen_1dt(:,jk)*zlmit(:,jk)
          zccvol_2dt(:,jk) = zccvol_2dt(:,jk)+zmult2(:,jk)*zcclen_2dt(:,jk)*zlmit(:,jk)
          zccvol_3dt(:,jk) = zccvol_3dt(:,jk)+zmult2(:,jk)*zcclen_3dt(:,jk)*zlmit(:,jk)
          zccvol_4dt(:,jk) = zccvol_4dt(:,jk)+zmult2(:,jk)*zcclen_4dt(:,jk)*zlmit(:,jk)
          zccvol_5dt(:,jk) = zccvol_5dt(:,jk)+zmult2(:,jk)*zcclen_5dt(:,jk)*zlmit(:,jk)
          zccvol(:,jk)     = zccvol(:,jk)    +zmult2(:,jk)*zcclen(:,jk)    *zlmit(:,jk)
        elsewhere
          zccvol_1dt(:,jk) = 0._dp
          zccvol_2dt(:,jk) = 0._dp
          zccvol_3dt(:,jk) = 0._dp
          zccvol_4dt(:,jk) = 0._dp
          zccvol_5dt(:,jk) = 0._dp
          zccvol(:,jk)     = zccvol_ges(:,jk)
        end where

        !Merken bis Aufteilung am Ende
        zcccov_ori(:,jk) = zcccov_ges(:,jk)
        zccvol_ges_ori(:,jk) = zccvol_ges(:,jk)


        !Update In-cloud Wert
        !Eispartikelanzahldichte und Eismasse

        WHERE(zccvol_ges(:,jk) > 1.e-10_dp)
 
          !Tendenz im ganzen Zeitschritt nach Sedimentation
          zxidt_cc(:) = ztmst * zccxite(:)

          zxidtsed_cc(:) = 0._dp

          !Sedimentation in Kondensstreifenfreien Teil
          ll1_1d(:) = B_cc(:,jk-1,jrow).gt.B_cc(:,jk,jrow)
          WHERE (ll1_1d(:))
            ztmp1_1d(:) = B_cc(:,jk,jrow)/MAX(B_cc(:,jk-1,jrow), zeps)
            ztmp1_1d(:) = MAX(MIN(ztmp1_1d(:),1._dp),0._dp)
            zxidtsed_cc(:) = zxisedf_cc(:,jk) * (1._dp - ztmp1_1d(:))
            zicesedf_cc(:,jk) = zicesedf_cc(:,jk) * ztmp1_1d(:)
            zxisedf_cc(:,jk) = zxisedf_cc(:,jk) * ztmp1_1d(:)
          ELSEWHERE
            zxidtsed_cc(:) = 0._dp
          END WHERE

          zicesed_cc(:,jk) = zicesed_cc(:,jk) + zicesedf_cc(:,jk)
          zxised_cc(:,jk) = zxised_cc(:,jk) + zxisedf_cc(:,jk)

          !Eispartikelanzahl in-cloud
          zccicnc(:,jk) = zicesed_cc(:,jk) / zccvol_ges(:,jk)
          !Eismasse in-cloud
          zxibcc(:) = zxised_cc(:,jk)/zccvol_ges(:,jk)

        ELSEWHERE

          !alles verdampft
          zxidtsed_cc(:) = zxised_cc(:,jk) + zxisedf_cc(:,jk)
          zxised_cc(:,jk) = 0._dp
          zicesed_cc(:,jk) = 0._dp
          zccicnc(:,jk) = 0._dp
          zxibcc(:) = 0._dp

        END WHERE

        !Verlust nach Sedimentation kann nicht mehr sein als 
        !Eis vorhanden
        ztmp1_1d(:) = MAX(zccxite(:), -(ztmst_rcp*pxtm1(:,jk,25)))
        zccxite(:) = MERGE(zccxite(:), ztmp1_1d(:), zccxite(:).gt.0._dp)

        END IF


        !--------------------------------------------------------
        !Abgleich

        zcorc2(:) = 0._dp
        
        !Wolke ohne Partikel
        ll1_1d(:) = (zccvol_ges(:,jk)*zccicnc(:,jk) > 1.e-10_dp)
        where (.not. ll1_1d(:))
          where(paclc(:,jk) .gt. 1.e-10_dp .AND. pxim1(:,jk) .gt. 0._dp)
            zcorc2(:) = zxibcc(:)*zccvol_ges(:,jk)
            zxib(:) = zxib(:) + zcorc2(:)/paclc(:,jk)
          elsewhere
            zxim1evp_cc(:,jk) = zxim1evp_cc(:,jk)+zxibcc(:)*zccvol_ges(:,jk)
          end where
        
          zcccov_ges(:,jk) = 0._dp
          zccvol_ges(:,jk) = 0._dp
          zcciwc(:,jk) = 0._dp
          zccicnc(:,jk) = 0._dp
          zcclen_ges2(:,jk) = 0._dp
        end where

!--- End CCMod -------------------------------------------------------------

!
!     ------------------------------------------------------------------
!       5.    Condensation/deposition and evaporation/sublimation
!
!             zlc       =  L_{v/s} / c_p
!
        zlc_1d(:)   = MERGE(zlsdcp(:),zlvdcp(:),lo2_1d(:))
        zqsm1_1d(:) = MERGE(zqsi_2d(:,jk)  , zqsw_2d(:,jk)  , lo2_1d(:))
        zqst1_1d(:) = MERGE(zqsip1_2d(:,jk), zqswp1_2d(:,jk), lo2_1d(:))

        zdqsdt_1d(:) = 1000._dp*(zqst1_1d(:)-zqsm1_1d(:))

!--- CCMod -------------------------------------------------------------

        zxievap(:) = zxidtstar_1d(:) + zxim1evp_1d(:)

        !Zuwachs durch Transport verdampft nicht, nur bei 
        !Sedimentation
        zxievap_cc(:) = zxim1evp_cc(:,jk) + zxidtsed_cc(:) &
                              + zcor44(:,jk)

!--- End CCMod -------------------------------------------------------------

        zxlevap(:) = (1.0_dp-zclcaux(:)) * zxldtstar_1d(:) + zxlm1evp_1d(:)

        zqvdt_1d(:) = ztmst*pqte(:,jk) + zevp(:)    + zsub(:)   &
                    + zxievap(:)       + zxlevap(:) + zxisub(:)

!--- CCMod -------------------------------------------------------------

        zqvdt_1d(:) = zqvdt_1d(:) + zxievap_cc(:) + zxisub_cc(:)              &
                        !emittierter Wasserdampf ohne Kondensstreifen
                        + (1._dp - B_co(:,jk,jrow))*zfh2o(:,jk)

!--- End CCMod -------------------------------------------------------------

        zdtdt_1d(:) = ztmst*ptte(:,jk)                                       &
                    -  zlvdcp(:)           *(zevp(:)  + zxlevap(:))          &
                    - (zlsdcp(:)-zlvdcp(:))*(zsmlt(:) + zximlt(:)+zimlt(:))

        zdtdt_1d(:) = zdtdt_1d(:) - zlsdcp(:)*(zsub(:)+zxievap(:)+zxisub(:))

!--- CCMod -------------------------------------------------------------

        zdtdt_1d(:) = zdtdt_1d(:) - (zlsdcp(:)-zlvdcp(:))*(zximlt_cc(:)    &
                                                                + zimlt_cc(:))
        zdtdt_1d(:) = zdtdt_1d(:) - zlsdcp(:)*(zxievap_cc(:)+zxisub_cc(:))

!--- End CCMod -------------------------------------------------------------

        zqp1_1d(:)  = pqm1(:,jk)+zqvdt_1d(:)
        zqp1_1d(:)  = MAX(zqp1_1d(:),0.0_dp)
!SF to put in lcover=false?
        ztp1_1d(:)   = ptm1(:,jk)+zdtdt_1d(:)

!--- CCMod -------------------------------------------------------------

        zccvol_all(:,jk) = 0.75_dp * zccvol_new(:,jk) + zccvol_ges(:,jk)
        zdqsat_1d(:) = zdtdt_1d(:)                                             &
                       + (zclcaux(:)+zccvol_all(:,jk))*( ztmst*zlc_1d(:)*pqte(:,jk)     &
                                + zlvdcp(:)*(zevp(:)+zxlevap(:))             &
                                + zlsdcp(:)*(zsub(:)+zxievap(:)+zxisub(:)  &  !zdtdtstar
                                                    +zxievap_cc(:)+zxisub_cc(:)) )
        zdqsat_1d(:) = zdqsat_1d(:)*zdqsdt_1d(:)/(1._dp+(zclcaux(:) &
                                + zccvol_all(:,jk))*zlc_1d(:)*zdqsdt_1d(:))

!--- End CCMod -------------------------------------------------------------

!SF end to put in lcover=false

        zxib(:)     = MAX(zxib(:),0.0_dp)
        zxlb(:)     = MAX(zxlb(:),0.0_dp)
        zxilb_1d(:) = zxib(:)+zxlb(:)
!
!       Diagnostics: relative humidity
!
        prelhum(:,jk) = pqm1(:,jk)/zqsm1_1d(:)
        prelhum(:,jk) = MAX(MIN(prelhum(:,jk),1._dp),0._dp)

        IF (lcover .AND. jk >= ncctop) THEN
!
!       define variables needed for cover scheme
!
!       zbetaqt = total water
!       zbetass = saturation mixing ratio adjusted to match qv
!       zwide   = current diagnosed distribution width
!
           zbetacl(:) = MAX(0.0_dp, pxlm1(:,jk)) + MAX(0.0_dp, pxim1(:,jk))
           zbetaqt(:) = MAX(cqtmin, pqm1(:,jk))  + zbetacl(:)
           zvartg(:)  = MAX(cqtmin,cvarmin*pqm1(:,jk))
           zwide(:)   = MAX(zvartg(:), pbetab(:,jk)-pbetaa(:,jk))
!
!       5.1 Turbulence: Skewness - equation solved implicitly
!           This solver only works if phmixtau has non-zero timescale
!
           zqtau_1d(:) = phmixtau(:,jk)+pvmixtau(:,jk) 

           zbqp1_1d(:) = -zdtime * zqtau_1d(:)
           zbqp1_1d(:) = EXP(zbqp1_1d(:))
           zbqp1_1d(:) = cbeta_pq - (cbeta_pq-pxskew(:,jk))*zbqp1_1d(:)
           zbqp1_1d(:) = MIN(zbqp1_1d(:), cbeta_pq_max)
           zbqp1_1d(:) = MAX(zbqp1_1d(:), cbeta_pq)

           zturbskew(:) = zdtime_rcp*(zbqp1_1d(:)-pxskew(:,jk))
!
!       5.2 Turbulence: variance - equation solved implicitly
!
           zbbap1_1d(:) = (cbeta_pq+pxskew(:,jk))**2    &
                        * (cbeta_pq+pxskew(:,jk)+1._dp) &
                        / (cbeta_pq*pxskew(:,jk)      )       !SF zeta

           zbbap1_1d(:) = zbbap1_1d(:)*pvdiffp(:,jk)/zwide(:)  !SF zprod

           zbbap1_1d(:) = zbbap1_1d(:)/zqtau_1d(:)                          &
                        + zvartg(:)                                         &
                        - (zbbap1_1d(:)/zqtau_1d(:) + zvartg(:) - zwide(:)) &
                          * EXP(-zdtime*zqtau_1d(:))

           zbbap1_1d(:) = MAX(zbbap1_1d(:),zvartg(:))
          
           ztmp1_1d(:) = 1._dp / cbeta_pq * zbetaqt(:) * (cbeta_pq + zbqp1_1d(:))


           zbbap1_1d(:) = MIN(zbbap1_1d(:), ztmp1_1d(:))

           zturbvar(:) = zdtime_rcp*(zbbap1_1d(:)-zwide(:))

           zbap1_1d(:) = zbetaqt(:)-cbeta_pq*zbbap1_1d(:)/(cbeta_pq+zbqp1_1d(:))
!
!          translated into apparent xl,xi,q and heat sources
!          first order effect only, effect of evaporation of
!          cloud on qsat taken into account in thermodynamic budget
!          but does not change the mixing term here since that
!          would require iteration and is therefore neglected
!
!          calculate values after one timestep
!
           iqidx_1d(:) = (zbqp1_1d(:)-cbeta_pq)/rbetak+1._dp
           iqidx_1d(:) = LOG(iqidx_1d(:))
           iqidx_1d(:) = (nbetaq/cbetaqs)*iqidx_1d(:) + 0.5_dp
           iqidx_1d(:) = INT(iqidx_1d(:))

           ztmp1_1d(:) =   cbeta_pq*(pbetass(:,jk) - zbap1_1d(:))                   &
                       / ( (zbetaqt(:)    - zbap1_1d(:))*(cbeta_pq + zbqp1_1d(:)) )
           ztmp1_1d(:) = MAX(MIN(ztmp1_1d(:), 1.0_dp), 0.0_dp)
           ztmp1_1d(:) = nbetax*ztmp1_1d(:)   !ztt
           ixidx_1d(:) = INT(ztmp1_1d(:))

           ll1_1d(:) = (ixidx_1d(:) == nbetax)

           DO jl=1,kproma

              iqidx = iqidx_1d(jl)
              ixidx = ixidx_1d(jl)

              ztmp2_1d(jl) = (ztmp1_1d(jl)-ixidx_1d(jl)      )*tbetai0(iqidx,ixidx+1)       &
                           + (ixidx_1d(jl)+1._dp-ztmp1_1d(jl))*tbetai0(iqidx,ixidx)
              ztmp3_1d(jl) = (ztmp1_1d(jl)-ixidx_1d(jl)      )*tbetai1(iqidx,ixidx+1)       &
                           + (ixidx_1d(jl)+1._dp-ztmp1_1d(jl))*tbetai1(iqidx,ixidx)
           ENDDO

           ztmp4_1d(:) = MERGE(1._dp, ztmp2_1d(:), ll1_1d(:))  !SF zbetai0
           ztmp5_1d(:) = MERGE(1._dp, ztmp3_1d(:), ll1_1d(:))  !SF zbetai1

           ztmp1_1d(:) = -zxilb_1d(:)*zclcaux(:) 
           ztmp2_1d(:) = zqsec*zqp1_1d(:)

           zgent_1d(:) = pqm1(:,jk)                                &
                       - (zbetaqt(:)    - zbap1_1d(:))*ztmp5_1d(:) &
                       + (pbetass(:,jk) - zbap1_1d(:))*ztmp4_1d(:) &
                       -  pbetass(:,jk)

           zgent_1d(:) = MAX(zgent_1d(:), ztmp1_1d(:))
           zgent_1d(:) = MIN(zgent_1d(:),ztmp2_1d(:))   !limit to qv

           ztmp1_1d(:) = MAX(zepsec,zxilb_1d(:))
           ztmp1_1d(:) = zxib(:)/ztmp1_1d(:)
           ztmp1_1d(:) = MAX(MIN(ztmp1_1d(:),1.0_dp),0.0_dp) !SFzifrac
           ztmp2_1d(:) = 1._dp - ztmp1_1d(:)

           ztmp3_1d(:) = zgent_1d(:)/MAX(zclcaux(:), zeps)

           zgenti(:) = zgent_1d(:)*ztmp1_1d(:)
           zgentl(:) = zgent_1d(:)*ztmp2_1d(:)

           ztmp4_1d(:) = zxib(:) + ztmp3_1d(:)*ztmp1_1d(:)
           ztmp4_1d(:) = MAX(ztmp4_1d(:),0.0_dp)
           zxib(:)     = MERGE(ztmp4_1d(:), zxib(:), locc_1d(:))

           ztmp4_1d(:) = zxlb(:) + ztmp3_1d(:)*ztmp2_1d(:)
           ztmp4_1d(:) = MAX(ztmp4_1d(:), 0.0_dp)
           zxlb(:)     = MERGE(ztmp4_1d(:), zxlb(:), locc_1d(:))

           zxilb_1d(:) = zxib(:)+zxlb(:)
!
!       5.3 Deposition/sublimation of cloud ice and condensation/
!           evaporation of liquid water due to changes in water vapour
!           and temperature (advection, convection, turbulent mixing,
!           evaporation of rain, sublimation and melting of snow).
!           Translate PDF laterally to calculate cloud
!           after one timestep
!
           zqvdt_1d(:) = zqvdt_1d(:)-zgent_1d(:)
           zdtdt_1d(:) = zdtdt_1d(:)+zlvdcp(:)*zgentl(:)+zlsdcp(:)*zgenti(:)

           zqp1_1d(:) = pqm1(:,jk)+zqvdt_1d(:)
           zqp1_1d(:) = MAX(zqp1_1d(:),0.0_dp)

           ztp1_1d(:) = ptm1(:,jk)+zdtdt_1d(:)

           zdqsat_1d(:) = zdtdt_1d(:)                                                       &
                        + zclcaux(:)*( ztmst*zlc_1d(:)*pqte(:,jk)                           &
                                     + zlvdcp(:)*(zevp(:)+zxlevap(:)-zgentl(:))             &
                                     + zlsdcp(:)*(zsub(:)+zxievap(:)+zxisub(:)-zgenti(:)) )   !zdtdtstar

           zdqsat_1d(:) = zdqsat_1d(:)*zdqsdt_1d(:)/(1._dp+zclcaux(:)*zlc_1d(:)*zdqsdt_1d(:))

           ztmp1_1d(:) = (pbetass(:,jk)-zqvdt_1d(:)+zdqsat_1d(:)-zbap1_1d(:))/zbbap1_1d(:)
           ztmp1_1d(:) = MAX(MIN(ztmp1_1d(:), 1.0_dp), 0.0_dp)
           ztmp1_1d(:) = nbetax*ztmp1_1d(:)  !ztt
           ixidx_1d(:) = INT(ztmp1_1d(:))

           ll1_1d(:) = (ixidx_1d(:) == nbetax)

           DO jl=1,kproma

              iqidx = iqidx_1d(jl)
              ixidx = ixidx_1d(jl)

              ztmp2_1d(jl) = (ztmp1_1d(jl)-ixidx_1d(jl)      )*tbetai0(iqidx,ixidx+1)       &
                           + (ixidx_1d(jl)+1._dp-ztmp1_1d(jl))*tbetai0(iqidx,ixidx)
              ztmp3_1d(jl) = (ztmp1_1d(jl)-ixidx_1d(jl)      )*tbetai1(iqidx,ixidx+1)       &
                           + (ixidx_1d(jl)+1._dp-ztmp1_1d(jl))*tbetai1(iqidx,ixidx)
           ENDDO

           ztmp4_1d(:) = MERGE(1._dp, ztmp2_1d(:), ll1_1d(:))  !SF zbetai0
           ztmp5_1d(:) = MERGE(1._dp, ztmp3_1d(:), ll1_1d(:))  !SF zbetai1

           zqcdif_1d(:) = (zbetaqt(:)-pbetaa(:,jk))                            *(1._dp-ztmp5_1d(:)) &
                        + (pbetaa(:,jk)+zqvdt_1d(:)-pbetass(:,jk)-zdqsat_1d(:))*(1._dp-ztmp4_1d(:))

           zqcdif_1d(:) = MAX(0.0_dp, zqcdif_1d(:))
           zqcdif_1d(:) = zqcdif_1d(:)-zbetacl(:)

           ztmp1_1d(:) = -zxilb_1d(:)*zclcaux(:)
           ztmp2_1d(:) = zqsec*zqp1_1d(:)

           zqcdif_1d(:) = MAX(zqcdif_1d(:), ztmp1_1d(:))
           zqcdif_1d(:) = MIN(zqcdif_1d(:), ztmp2_1d(:))  ! limit to qv

        ELSE !lcover=.false. or jk < ncctop

           zqcdif_1d(:) = (zqvdt_1d(:)-zdqsat_1d(:))*zclcaux(:)

           ztmp1_1d(:) = -zxilb_1d(:)*zclcaux(:)
           ztmp2_1d(:) = zqsec*zqp1_1d(:)

           zqcdif_1d(:) = MAX(zqcdif_1d(:), ztmp1_1d(:))
           zqcdif_1d(:) = MIN(zqcdif_1d(:), ztmp2_1d(:))  ! limit to qv

!--- CCMod -------------------------------------------------------------

           zccvol_all(:,jk) = zccvol_new(:,jk) + zccvol_ges(:,jk)
 
           zqccdif_1d(:) = (zqvdt_1d(:)-zdqsat_1d(:))*zccvol_ges2(:,jk)
           
           ztmp1_1d(:) = -zxibcc(:)*zccvol_ges(:,jk)
           ztmp2_1d(:) = zqsec*zqp1_1d(:)*zccvol_ges(:,jk)/MAX(zclcaux(:) &
                                                        +zccvol_all(:,jk),zeps)

           zqccdif_1d(:) = MAX(zqccdif_1d(:), ztmp1_1d(:))
           zqccdif_1d(:) = MIN(zqccdif_1d(:), ztmp2_1d(:))  ! limit to qv

           zqcodif_1d(:) = (zqvdt_1d(:)-zdqsat_1d(:))*0.75_dp*zccvol_new(:,jk)
           
           ztmp1_1d(:) = -zcciwc_new(:,jk)
           ztmp2_1d(:) = zqsec*zqp1_1d(:)*zccvol_new(:,jk)/MAX(zclcaux(:) &
                                                         +zccvol_all(:,jk),zeps)

           zqcodif_1d(:) = MAX(zqcodif_1d(:), ztmp1_1d(:))
           zqcodif_1d(:) = MIN(zqcodif_1d(:), ztmp2_1d(:))  ! limit to qv

!--- End CCMod -------------------------------------------------------------

        END IF !lcover

        ll1_1d(:) = (zqcdif_1d(:) < 0._dp)     !SF cloud dissipation
        
        ztmp1_1d(:) = MAX(zepsec, zxilb_1d(:))
        ztmp1_1d(:) = zxib(:) / ztmp1_1d(:)
        ztmp1_1d(:) = MAX(MIN(ztmp1_1d(:), 1.0_dp), 0.0_dp) !zifrac

        ztmp1_1d(:) = MERGE(ztmp1_1d(:), 1._dp, ll1_1d(:))

        ztmp2_1d(:) = zqcdif_1d(:)*(1.0_dp-ztmp1_1d(:))
        zcnd(:) = MERGE(ztmp2_1d(:), 0._dp, ll1_1d(:))

        IF (nicnc <= 1) THEN
           zdep(:) = zqcdif_1d(:)*ztmp1_1d(:)
        ELSE

!--- CCMod -------------------------------------------------------------

!           zdep(:) = zqinucl(:,jk)*ztmp1_1d(:)

          ! for saturation adjustment
          zdep(:) = zqcdif_1d(:)*ztmp1_1d(:)
          zdepcc(:) = zqccdif_1d(:)
          zdepco(:) = zqcodif_1d(:)

!--- End CCMod -------------------------------------------------------------

        ENDIF

        ll2_1d(:) = (.NOT. ll1_1d(:)) .AND. &  !SF no cloud dissipation
                    (.NOT. lo2_1d(:))          !SF condensation

        zdep(:) = MERGE(0._dp, zdep(:), ll2_1d(:))

!--- CCMod -------------------------------------------------------------

        zdepcc(:) = MERGE(0._dp, zdepcc(:), ll2_1d(:))
        zdepco(:) = MERGE(0._dp, zdepco(:), ll2_1d(:))

!--- End CCMod -------------------------------------------------------------

!--- Included/changed for prognostic CDNC/IC scheme --------------------
!    Use standard condensation for empirical Lin & Leaitch approach and 
!    explicit condensation after Levkov et al. 1992 for the explicit
!    activation schemes that allow for supersaturation:

        IF (ncdnc == 1) THEN
           zcnd(:) = MERGE(zqcdif_1d(:), zcnd(:), ll2_1d(:))
        ELSE IF (ncdnc > 1) THEN
           ztmp1_1d(:) = 0.5_dp*api*zdw0*ztmst                  &
                       * swat(1:kproma,jk,jrow)*zcdnc(:,jk)*zclcaux(:) &
                       * zrho_rcp(:,jk)/zastbstw(:,jk)
           ztmp1_1d(:) = MIN(ztmp1_1d(:), zqcdif_1d(:))
           zcnd(:) = MERGE(ztmp1_1d(:), zcnd(:), ll2_1d(:))
        ENDIF

!--- End included for CDNC/IC scheme -----------------------------------
!
!
!       5.4 Accounting for cloud evaporation in clear air and
!           checking for supersaturation
!
!--- CCMod -------------------------------------------------------------

        ztp1tmp(:) = ztp1_1d(:)+zlvdcp(:)*zcnd(:)+zlsdcp(:)*(zdep(:) &
                                                     +zdepcc(:)+zdepco(:))
        zqp1tmp(:) = zqp1_1d(:)-zcnd(:)-zdep(:)-zdepcc(:)-zdepco(:)

!--- End CCMod -------------------------------------------------------------

        zxip1_1d(:) = zxised_2d(:,jk) + ztmst*zxite(:) - zxievap(:) + zgenti(:) + zdep(:)
        zxip1_1d(:) = MAX(zxip1_1d(:), 0._dp)

!--- CCMod -------------------------------------------------------------

        zxip1_cc(:) = zxised_cc(:,jk) - zxievap_cc(:) + zdepcc(:) + zdepco(:)
        zxip1_cc(:) = MAX(zxip1_cc(:), 0._dp)

        lo2_1d(:)   = (ztp1tmp(:) < cthomi) .OR.                            &
                        (ztp1tmp(:) < tmelt .AND. ((zxip1_1d(:) > csecfrl   &
                                              .OR. zxip1_cc(:) > csecfrl))  &
                       .AND. zsusatw_2d(:,jk) < zeps)        

!--- End CCMod -------------------------------------------------------------

        ztmp1_1d(:) = 1000._dp*ztp1tmp(:)
        it_1d(:)    = NINT(ztmp1_1d(:))

        ll_look(:,jk)     = (it_1d(:)<jptlucu1 .OR. it_1d(:)>jptlucu2)
        IF (ANY(ll_look(1:kproma,jk))) lookupoverflow = .TRUE.

        it_1d(:)    = MAX(MIN(it_1d(:),jptlucu2),jptlucu1)

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

        ll1_1d(:) = (zes_1d(:) < 0.4_dp)  !SF LO

        zcor_1d(:) = 1._dp/(1._dp-vtmpc1*zes_1d(:))

        zqsp1tmp_1d(:) = zes_1d(:)*zcor_1d(:)
        zoversat_1d(:) = zqsp1tmp_1d(:)*0.01_dp
       
        zrhtest_1d(:) = pqm1(:,jk)/zqsm1_1d(:)
        zrhtest_1d(:) = MIN(zrhtest_1d(:), 1._dp)
        zrhtest_1d(:) = zrhtest_1d(:)*zqsp1tmp_1d(:)

!changed for cirrus scheme, Ulrike Lohmann, 17.12.2006
        zesw_1d(:) = zlucuaw_1d(:)/papp1(:,jk)
        zesw_1d(:) = MIN(zesw_1d(:), 0.5_dp)

        zcorw_1d(:) = 1._dp/(1._dp-vtmpc1*zesw_1d(:))

        zqsp1tmpw_1d(:) = zesw_1d(:)*zcorw_1d(:)
        zoversatw_1d(:) = 0.01_dp*zqsp1tmpw_1d(:)
!end changed for cirrus scheme

        zqst1_1d(:) = MERGE(zlucuap1_1d(:),zlucuawp1_1d(:),lo2_1d(:))
        zqst1_1d(:) = zqst1_1d(:)/papp1(:,jk)
        zqst1_1d(:) = MIN(zqst1_1d(:),0.5_dp)
        zqst1_1d(:) = zqst1_1d(:)/(1._dp-vtmpc1*zqst1_1d(:))

        zdqsdt_1d(:) = 1000._dp*(zqst1_1d(:)-zqsp1tmp_1d(:))

        zlc_1d(:) = MERGE(zlsdcp(:),zlvdcp(:),lo2_1d(:))

        ztmp1_1d(:) = zlc_1d(:)*zdqsdt_1d(:)
        ztmp2_1d(:) = zqsp1tmp_1d(:)*zcor_1d(:)*zlucub_1d(:)
        zlcdqsdt_1d(:) = MERGE(ztmp1_1d(:), ztmp2_1d(:), ll1_1d(:))

        zqcon_1d(:) = 1._dp/(1._dp+zlcdqsdt_1d(:))

        ztmp1_1d(:) = zqsp1tmpw_1d(:)+zoversatw_1d(:)
        ztmp1_1d(:) = MIN(ztmp1_1d(:),zqsp1tmp_1d(:)*1.3_dp)
        zqsp1tmphet_1d(:) = MERGE(ztmp1_1d(:), 0._dp, lo2_1d(:))

        ll1_1d(:) = (nicnc <= 1) 
! ll1_1d = false --> for cirrus scheme, Ulrike Lohmann, 17.12.2006

        ll2_1d(:) = (zqp1tmp(:) > (zqsp1tmp_1d(:)  + zoversat_1d(:) ) )
        ll3_1d(:) = (zqp1tmp(:) > (zqsp1tmpw_1d(:) + zoversatw_1d(:)) )
        ll4_1d(:) = (zqp1tmp(:) > zqsp1tmphet_1d(:))
        ll5_1d(:) = (ztp1tmp(:) >= cthomi)

        ztmp1_1d(:) = (zqp1tmp(:) - zqsp1tmp_1d(:)  - zoversat_1d(:) )*zqcon_1d(:)
        ztmp2_1d(:) = (zqp1tmp(:) - zqsp1tmpw_1d(:) - zoversatw_1d(:))*zqcon_1d(:)
        ztmp3_1d(:) = (zqp1tmp(:) - zqsp1tmphet_1d(:)                )*zqcon_1d(:)

        !SF ice cloud cases:         
        ll6_1d(:) = (lo2_1d(:) .AND. ll1_1d(:)         .AND. ll2_1d(:))                 &
                    .OR.                                                                &
                    (lo2_1d(:) .AND. (.NOT. ll1_1d(:)) .AND. ll2_1d(:) .AND. ll5_1d(:))
        ztmp4_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll6_1d(:))

        ll6_1d(:) = lo2_1d(:) .AND. (.NOT. ll1_1d(:)) .AND. ll3_1d(:)    &
                              .AND. (.NOT. ll5_1d(:)) .AND. (.NOT. lhet)
        ztmp4_1d(:) = MERGE(ztmp2_1d(:), ztmp4_1d(:), ll6_1d(:))

        ll6_1d(:) = lo2_1d(:) .AND. (.NOT. ll1_1d(:)) .AND. ll4_1d(:)    &
                              .AND. (.NOT. ll5_1d(:)) .AND. lhet
        ztmp4_1d(:) = MERGE(ztmp3_1d(:), ztmp4_1d(:), ll6_1d(:))

        ll6_1d(:) = (lo2_1d(:) .AND. ll2_1d(:))

        ztmp4_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll6_1d(:))

        zdep(:) = zdep(:) + ztmp4_1d(:)
       
        !SF water cloud cases:
        ll6_1d(:) = (.NOT. lo2_1d(:)) .AND. ll2_1d(:)
        ztmp4_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll6_1d(:))
        zcnd(:)     = zcnd(:) + ztmp4_1d(:)

        !SF final corrections for zdep and zcnd:
        ztmp5_1d(:) = zqp1_1d(:)-zrhtest_1d(:)
        ztmp5_1d(:) = MAX(ztmp5_1d(:), 0._dp)

!--- CCMod -------------------------------------------------------------

!        ll1_1d(:) = (zdep(:)        > 0._dp        )

        ll1_1d(:) = (zdep(:) > 0._dp .OR. zdepcc(:) > 0._dp &
                                       .OR. zdepco(:) > 0._dp)

!--- End CCMod -------------------------------------------------------------

        ll2_1d(:) = (zcnd(:)        > 0._dp        )
        ll3_1d(:) = (zqp1tmp(:)     < zrhtest_1d(:))
        ll4_1d(:) = (zqsp1tmp_1d(:) <= zqsm1_1d(:) )

        ll5_1d(:) = lo2_1d(:) .AND. ll1_1d(:) .AND. ll3_1d(:) .AND. ll4_1d(:)

!--- CCMod -------------------------------------------------------------

!        zdep(:) = MERGE(ztmp5_1d(:), zdep(:), ll5_1d(:)) 

        ztmp1_1d(:) = zccvol_ges(:,jk)/MAX(zclcaux(:)+zccvol_all(:,jk),zeps)
        ztmp2_1d(:) = zccvol_new(:,jk)/MAX(zclcaux(:)+zccvol_all(:,jk),zeps)
        zdep(:) = MERGE((1._dp-ztmp1_1d(:)-ztmp2_1d(:))*ztmp5_1d(:), zdep(:), &
                                                                      ll5_1d(:))
        zdepcc(:) = MERGE(ztmp1_1d(:)*ztmp5_1d(:), zdepcc(:), ll5_1d(:))
        zdepco(:) = MERGE(ztmp2_1d(:)*ztmp5_1d(:), zdepco(:), ll5_1d(:))

!--- End CCMod -------------------------------------------------------------

        ll5_1d(:) = (.NOT. lo2_1d(:)) .AND. ll2_1d(:) .AND. ll3_1d(:) .AND. ll4_1d(:)
        zcnd(:) = MERGE(ztmp5_1d(:), zcnd(:), ll5_1d(:))
!
!       5.5 Change of in-cloud water due to deposition/sublimation and
!           condensation/evaporation (input for cloud microphysics)
!
        zrelhum_1d(:) = zqp1tmp(:)/zqsp1tmp_1d(:)

        ztmp1_1d(:) = zdep(:) + zgenti(:)  
        ztmp1_1d(:) = MAX(ztmp1_1d(:), 0._dp) !SF deposition zdepos
       
        ztmp2_1d(:) = zcnd(:) + zgentl(:)
        ztmp2_1d(:) = MAX(ztmp2_1d(:), 0._dp) !SF condensation zcond 

        ztmp3_1d(:) = zxib(:)+zdep(:)/MAX(zclcaux(:), zeps) 
        ztmp3_1d(:) = MAX(ztmp3_1d(:), 0._dp)

        ztmp4_1d(:) = zxlb(:)+zcnd(:)/MAX(zclcaux(:), zeps)
        ztmp4_1d(:) = MAX(ztmp4_1d(:), 0._dp)

        zxib(:) = MERGE(ztmp3_1d(:), zxib(:), locc_1d(:))
        zxlb(:) = MERGE(ztmp4_1d(:), zxlb(:), locc_1d(:))

        ll1_1d(:) = (.NOT. locc_1d(:))                                       &
                  .AND. ( (ztmp1_1d(:) > 0._dp) .OR. (ztmp2_1d(:) > 0._dp) )        

        ztmp3_1d(:) = MAX(MIN(zrelhum_1d(:), 1.0_dp), 0.01_dp)
        zclcaux(:)  = MERGE(ztmp3_1d(:), zclcaux(:), ll1_1d(:))

        ztmp3_1d(:) = ztmp1_1d(:) / MAX(zclcaux(:), zeps)
        ztmp4_1d(:) = ztmp2_1d(:) / MAX(zclcaux(:), zeps)
      
        zxib(:) = MERGE(ztmp3_1d(:), zxib(:), ll1_1d(:))
        zxlb(:) = MERGE(ztmp4_1d(:), zxlb(:), ll1_1d(:))

!--- CCMod -------------------------------------------------------------

        !Update In-cloud-Werte
        ztmp1_1d(:) = MERGE(zdepcc(:)/zccvol_ges(:,jk),0._dp, & 
                                          zccvol_ges(:,jk).gt.1.e-10_dp)
        ztmp5_1d(:) = zxibcc(:)+ztmp1_1d(:)
        ztmp5_1d(:) = MAX(ztmp5_1d(:), 0._dp)
        ztmp2_1d(:) = MAX(zxibcc(:)-(B_co(:,jk,jrow) - zfrac1(:,jk))*zfh2o(:,jk)*zccvol_ges(:,jk),0._dp)
        zxibcc(:) = MERGE(ztmp5_1d(:), zxibcc(:), zccvol_ges(:,jk).gt.zeps)
        ztmp1_1d(:) = MERGE(zdepco(:)/zccvol_new(:,jk),0._dp, &
                                          zccvol_new(:,jk).gt.1.e-10_dp)
        ztmp5_1d(:) = MERGE(zcciwc_new(:,jk)/zccvol_new(:,jk)+ztmp1_1d(:),0._dp, &
                                          zccvol_new(:,jk).gt.1.e-10_dp)
        ztmp5_1d(:) = MAX(ztmp5_1d(:), 0._dp)
        zxibco(:) = MERGE(ztmp5_1d(:), 0._dp, zccvol_new(:,jk).gt.zeps)
        ll3_1d(:) = zdepcc(:) .lt. 0._dp .OR. zdepco(:) .lt. 0._dp
        zdepcc(:) = zdepcc(:) + zdepco(:) + zdepco_form(:,jk)


!--- End CCMod -------------------------------------------------------------

        !SF re-define locc_1d which was no longer true:
        locc_1d(:) = (zclcaux(:) > zeps)

        !SF update cdnc:
        ll1_1d(:) = locc_1d(:) .AND. (zxlb(:) > cqtmin)
        ll2_1d(:) = ll1_1d(:) .AND. (zcdnc(:,jk) <= cdncmin) !if there was no previous nucleation

        ztmp1_1d(:) = pcdncact(:,jk) - zcdnc(:,jk)
        ztmp1_1d(:) = MAX(0._dp, ztmp1_1d(:))
        zqlnuc(:,jk)    = MERGE(ztmp1_1d(:), zqlnuc(:,jk), ll2_1d(:))

        ztmp1_1d(:)      = MERGE(ztmp1_1d(:), 0._dp, ll2_1d(:))
        zcdnc(:,jk)      = zcdnc(:,jk) + ztmp1_1d(:)
        qnuc(1:kproma,jk,jrow)  = qnuc(1:kproma,jk,jrow) + zdt*ztmp1_1d(:)

        ztmp1_1d(:) = MAX(zcdnc(:,jk), cdncmin)
        zcdnc(:,jk) = MERGE(ztmp1_1d(:), cqtmin, ll1_1d(:)) 

        !SF update icnc:
        ll1_1d(:) = locc_1d(:) .AND. (zxib(:) > cqtmin) 
        ll2_1d(:) = ll1_1d(:)  .AND. (zicnc(:,jk) <= zicemin)

        IF (nicnc <= 1) THEN
           ztmp1_1d(:) = 0.75_dp/(api*zrhoice)*zrho(:,jk)*zxib(:)/zrid_2d(:,jk)**3
        ELSE
           !ztmp1_1d(:) = MIN(znicex(:,jk),zap(:,jk))             !mz_sb_20170710: [znicex]=[1/m3] while [zap]=[1/cm3]
           ztmp1_1d(:) = MIN(znicex(:,jk),(zap(:,jk)*1.e6_dp))    !mz_sb_20170710
        ENDIF

        zicnc(:,jk) = MERGE(ztmp1_1d(:), zicnc(:,jk), ll2_1d(:))

        ztmp1_1d(:) = MAX(zicnc(:,jk), zicemin)
        zicnc(:,jk) = MERGE(ztmp1_1d(:), cqtmin, ll1_1d(:)) 

        ztp1tmp(:) = ztp1_1d(:) + zlvdcp(:)*zcnd(:) + zlsdcp(:)*zdep(:)

!--- CCMod -------------------------------------------------------------

        zcccov_all(:,jk) = zcccov_ges(:,jk) + zcccov_new(:,jk)
        WHERE(zcccov_all(:,jk).gt.1.e-10_dp .AND. &
                            (zclcaux(:)+zcccov_all(:,jk)) .gt. 1._dp)
          ztmp1_1d(:) = (1._dp - zclcaux(:))/zcccov_all(:,jk) 
          zcorclc(:) = (1._dp-ztmp1_1d(:))*(zccvol_ges(:,jk)*zxibcc(:) &
                                              +zccvol_new(:,jk)*zxibco(:))
          zxib(:) = zxib(:)+zcorclc(:)/MAX(zclcaux(:),zeps)
          zccvol_all(:,jk) = ztmp1_1d(:)*zccvol_all(:,jk)
          zccvol_ges(:,jk) = ztmp1_1d(:)*zccvol_ges(:,jk)
          zccvol_new(:,jk) = ztmp1_1d(:)*zccvol_new(:,jk)
          zcccov_all(:,jk) = ztmp1_1d(:)*zcccov_all(:,jk)
          zcccov_ges(:,jk) = ztmp1_1d(:)*zcccov_ges(:,jk)
          zcccov_new(:,jk) = ztmp1_1d(:)*zcccov_new(:,jk) 
        END WHERE

        !zu kleine Kondensstreifen verschwinden
        !Abgleich 
        !Grenze liegt bei 1.e-20_dp kg/kg im gittermittel
        ztmp1_1d(:) = 1.e-20_dp
        ztmp2_1d(:) = zxibcc(:)*zccvol_ges(:,jk)
        ll1_1d(:) = (ztmp2_1d(:) .ge. ztmp1_1d(:))
        where (.not. ll1_1d(:))
          where(zclcaux(:) .gt. 1.e-10_dp .AND. pxim1(:,jk) .gt. 0._dp)
            zcorc2(:) =  zcorc2(:) + ztmp2_1d(:)
            zxib(:) = zxib(:) + ztmp2_1d(:)/zclcaux(:)
          elsewhere
            zxievap_cc(:) = zxievap_cc(:)+ztmp2_1d(:)
          end where

          zxibcc(:) = 0._dp
          zcccov_ges(:,jk) = 0._dp
          zccvol_ges(:,jk) = 0._dp
          zccicnc(:,jk) = 0._dp
          zcclen_ges2(:,jk) = 0._dp
        end where
        ztmp2_1d(:) = zxibco(:)*zccvol_new(:,jk)
        ll1_1d(:) = (ztmp2_1d(:) .ge. ztmp1_1d(:))
        where (.not. ll1_1d(:))
          where(zclcaux(:) .gt. 1.e-10_dp .AND. pxim1(:,jk) .gt. 0._dp)
            zcorc2(:) =  zcorc2(:) + ztmp2_1d(:)
            zxib(:) = zxib(:) + ztmp2_1d(:)/zclcaux(:)
          elsewhere
            zxievap_cc(:) = zxievap_cc(:)+ztmp2_1d(:)
          end where

          zxibco(:) = 0._dp
          zcccov_new(:,jk) = 0._dp
          zccvol_new(:,jk) = 0._dp
          zccicnc_new(:,jk) = 0._dp
          zcclen_new(:,jk) = 0._dp
        end where

!--- End CCMod -------------------------------------------------------------

        !Bildung von Eispartikeln wenn Vertikalgeschwindigkeit zu klein
        !Anteil neue Wolke
        znicex(:,jk) = 0._dp
        zmult_1d(:) = (paclc(:,jk) - pxtm1(:,jk,4)) / MAX(paclc(:,jk),zeps)
        zmult_1d(:) = MIN(MAX(zmult_1d(:), 0._dp), 1._dp)
        zmult_1d(:) = MERGE(zmult_1d(:), 0._dp, paclc(:,jk).gt.zeps)
        WHERE (ptm1(:,jk) < cthomi .AND. (paclc(:,jk) >= zeps .AND.    &
                  paclc(:,jk) > pxtm1(:,jk,4)) .AND. &
                  zvervx_2d(:,jk) < 1.e-2_dp) !kleiner 0.1 mm/s

          !Nukleation ist noetig, aber w ist zu klein
          !Berechnung durch IWC

          !Radius aus Fit von Wang und Sassen 2002
          ztmp1_1d(:) = MAX(ptm1(:,jk), 195._dp)
          ztmp1_1d(:) = MIN(ztmp1_1d(:), 240._dp)
          zri3(:,jk) = (1.2_dp * ztmp1_1d(:) - 221.6_dp) * 0.37_dp
          zri3(:,jk) = zri3(:,jk) * 1.e-6_dp
          !Anteil an Wolkeneis
          ztmp2_1d(:) = zdep(:) / paclc(:,jk)
          znicex3(:,jk) = 0.75_dp * zrho(:,jk) * ztmp2_1d(:)    &
                             / (api * zrhoice * zri3(:,jk)**3._dp)
          znicex3(:,jk) = MAX(znicex3(:,jk), 0._dp)
          znicex(:,jk) = zmult_1d(:)*znicex3(:,jk)
        END WHERE

        !Kontrolle Eispartikelanzahl
        !Wurden zu wenig Eispartikel gebildet?
        !alte Wolke:

        WHERE(ptm1(:,jk).lt.cthomi .AND. paclc(:,jk).gt.zeps &
                   .AND. pxtm1(:,jk,4).gt.zeps .AND. zdep(:) .gt. 0._dp)

          !Abgleich Depositionswachstum und deponierbares Wasser
          ztmp(:,jk) = (1._dp-zmult_1d(:))*zdep(:)
          ztmp1_1d(:) = ztmp(:,jk) - zqinucl(:,jk)

          WHERE(ztmp1_1d(:) .gt. zeps)
            ztmp2_1d(:) = ztmp(:,jk)/MAX(znicex1(:,jk)+zicncq(:,jk),zeps)
            ztmp4_1d(:) = MAX(zqinucl(:,jk)/zicncq(:,jk), &
                               4._dp/3._dp * 10.e-6_dp**3._dp * api * zrhoice / zrho(:,jk))
            WHERE( ztmp2_1d(:) .gt. 2.5_dp*ztmp4_1d(:))
              znicex_kor1(:,jk) = 2._dp/5._dp * ztmp1_1d(:)/MAX(ztmp4_1d(:),zeps) &
                                         - znicex1(:,jk)
              znicex_kor1(:,jk) = MAX(znicex_kor1(:,jk), 0._dp)
              znicex(:,jk) = (1._dp-zmult_1d(:))*znicex_kor1(:,jk) + znicex(:,jk)
            END WHERE
          END WHERE

        END WHERE


        !neue Wolke:

        WHERE(ptm1(:,jk).lt.cthomi .AND. paclc(:,jk).gt.zeps &
                   .AND. zdep(:) .gt. 0._dp)

          !Abgleich Depositionswachstum und deponierbares Wasser
          ztmp(:,jk) = zmult_1d(:)*zdep(:)
          ztmp1_1d(:) = ztmp(:,jk) - zqinucl_new(:,jk)

          WHERE(ztmp1_1d(:) .gt. zeps)
            ztmp5_1d(:) = MAX(znicex2(:,jk), znicex3(:,jk))
            ztmp2_1d(:) = ztmp(:,jk)/MAX(ztmp5_1d(:),zeps)
            ztmp4_1d(:) = MAX(zqinucl_new(:,jk)/MAX(ztmp5_1d(:),zeps), &
                               4._dp/3._dp * 10.e-6_dp**3._dp * api * zrhoice / zrho(:,jk))
            WHERE( ztmp2_1d(:) .gt. 2.5_dp*ztmp4_1d(:))
              znicex_kor2(:,jk) = 2._dp/5._dp * ztmp(:,jk)/MAX(ztmp4_1d(:),zeps) &
                                         - MAX(znicex2(:,jk),znicex3(:,jk))
              znicex(:,jk) = zmult_1d(:)*znicex_kor2(:,jk) + znicex(:,jk)
            END WHERE
          END WHERE

        END WHERE


        !--- Update ICNC

        zri(:,jk)=MAX(zri(:,jk), 1.e-6_dp)

        ll_ice(:,jk) = (ptm1(:,jk) < cthomi)

        ztmp1(:,jk) = 1.e6_dp*zap(:,jk) - zicnc(:,jk)
        ztmp2(:,jk) = znicex(:,jk)
        ztmp(:,jk)  = MIN(ztmp1(:,jk),ztmp2(:,jk))
        ztmp(:,jk)  = MAX(ztmp(:,jk), 0._dp)

        zninucl(:,jk) = MERGE(zninucl(:,jk)+ztmp(:,jk), 0._dp,ll_ice(:,jk))

        zicnc(:,jk) = zicnc(:,jk) + znicex(:,jk)
!
!     ------------------------------------------------------------------
!       6.    Freezing of cloud water
!
!       6.1   Freezing of all cloud water for T < 238 K
!
        ll1_1d(:) = (ztp1tmp(:) <= cthomi)

        ztmp1_1d(:) = zfrl(:,jk) + zxlb(:)*zclcaux(:)
        zfrl(:,jk)  = MERGE(ztmp1_1d(:), zfrl(:,jk), ll1_1d(:))

        ztmp1_1d(:) = zxib(:)+zxlb(:)
        zxib(:)     = MERGE(ztmp1_1d(:), zxib(:), ll1_1d(:))

        zxlb(:)     = MERGE(0._dp, zxlb(:), ll1_1d(:))

!--- Included for prognostic CDNC/IC scheme ----------------------------
        ztmp1_1d(:)     = zcdnc(:,jk)-cdncmin
        ztmp1_1d(:)     = MAX(ztmp1_1d(:), 0._dp)
        ztmp2_1d(:)     = qfre(1:kproma,jk,jrow) - zdt*ztmp1_1d(:)
        qfre(1:kproma,jk,jrow) = MERGE(ztmp2_1d(:), qfre(1:kproma,jk,jrow), ll1_1d(:))

        ztmp2_1d(:) = zicnc(:,jk) + ztmp1_1d(:)
        zicnc(:,jk) = MERGE(ztmp2_1d(:), zicnc(:,jk), ll1_1d(:))

        zcdnc(:,jk) = MERGE(cqtmin, zcdnc(:,jk), ll1_1d(:)) 

!--- End included for CDNC/IC scheme -----------------------------------
!
!       6.2   Freezing of cloud water between 238 and 273 K
!
        lo_1d(:) =     (zxlb(:)     > cqtmin  ) &
                 .AND. (ztp1tmp(:)  < tmelt   ) &
                 .AND. (ztp1tmp(:)  > cthomi  ) & 
                 .AND. (zcdnc(:,jk) >= cdncmin) &
                 .AND. locc_1d(:)

!---Changed for prognostic CDNC/IC scheme ------------------------------
!   (Replaced pacdnc by zcdnc which is set above)

!--------- new freezing parameterisations (Lohmann & Diehl, 2006)
           ! corinna: included for contact/immersion freezing by dust and soot

         ztmp1_1d(:) = 1._dp + 1.26_dp*6.6E-8_dp / (zrwetki_2d(:,jk)+zeps) * (101325._dp/papp1(:,jk))  & !SF ccbcki
                                                 * (ztp1tmp(:)/tmelt)
         ztmp2_1d(:) = 1._dp + 1.26_dp*6.6E-8_dp / (zrwetai_2d(:,jk)+zeps) * (101325._dp/papp1(:,jk))  & !SF ccduai
                                                 * (ztp1tmp(:)/tmelt)
         ztmp3_1d(:) = 1._dp + 1.26_dp*6.6E-8_dp / (zrwetci_2d(:,jk)+zeps) * (101325._dp/papp1(:,jk))  & !SF ccduci
                                                 * (ztp1tmp(:)/tmelt)

         zetaair_1d(:) = 1.e-5_dp  &
                       * ( 1.718_dp + 0.0049_dp*(ztp1tmp(:)-tmelt)                      &
                                    - 1.2e-5_dp*(ztp1tmp(:)-tmelt)*(ztp1tmp(:)-tmelt) )

!Ulrike: for thermophoresis
         zkair_1d(:)    = 4.1867e-3_dp*(5.69_dp + 0.017_dp*(ztp1tmp(:)-tmelt)) ! eq. 13-18a P&K
!end thermophoresis

         ll1_1d(:)   = (zrwetki_2d(:,jk) < zeps)
         zdfarbcki_1d(:) = ak * ztp1tmp(:) * ztmp1_1d(:) &
                         / ( 6._dp*api*zetaair_1d(:)*(zrwetki_2d(:,jk)+zeps) )
         zdfarbcki_1d(:) = MERGE(0._dp, zdfarbcki_1d(:), ll1_1d(:))

         ll2_1d(:)   = (zrwetai_2d(:,jk) < zeps)
         zdfarduai_1d(:) = ak * ztp1tmp(:) * ztmp2_1d(:) &
                         / ( 6._dp*api*zetaair_1d(:)*(zrwetai_2d(:,jk)+zeps) )
         zdfarduai_1d(:) = MERGE(0._dp, zdfarduai_1d(:), ll2_1d(:))

         ll3_1d(:)   = (zrwetci_2d(:,jk) < zeps)
         zdfarduci_1d(:) = ak * ztp1tmp(:) * ztmp3_1d(:) &
                         / ( 6._dp*api*zetaair_1d(:)*(zrwetci_2d(:,jk)+zeps) )
         zdfarduci_1d(:) = MERGE(0._dp, zdfarduci_1d(:), ll3_1d(:))

!Ulrike: for thermophoresis
         IF (lthermo) THEN
             zknbcki_1d(:) = 7.37_dp*ztp1tmp(:)/(2.88e5_dp*papp1(:,jk)*(zrwetki_2d(:,jk)+zeps))
             zknbcki_1d(:) = MERGE(1._dp, zknbcki_1d(:), ll1_1d(:))
 
             zftbcki_1d(:) = 0.4_dp &
                           * ( 1._dp + zknbcki_1d(:)*(1.45_dp + 0.4_dp*EXP(-1._dp/zknbcki_1d(:)))) &
                           * (zkair_1d(:)+2.5_dp*zknbcki_1d(:)*zkbc)                                  &
                           / ( (1._dp+3._dp*zknbcki_1d(:))                                         &
                             * (2._dp*zkair_1d(:)+zkbc*(5._dp*zknbcki_1d(:)+1._dp)))
 
             zftbcki_1d(:) = MERGE(0._dp, zftbcki_1d(:), ll1_1d(:))

             zknduai_1d(:) = 7.37_dp*ztp1tmp(:)/(2.88e5_dp*papp1(:,jk)*(zrwetai_2d(:,jk)+zeps))
             zknduai_1d(:) = MERGE(1._dp, zknduai_1d(:), ll2_1d(:))
 
             zftduai_1d(:) = 0.4_dp &
                          * ( 1._dp + zknduai_1d(:)*(1.45_dp + 0.4_dp*EXP(-1._dp/zknduai_1d(:)))) &
                          * (zkair_1d(:)+2.5_dp*zknduai_1d(:)*zkdu)                                  &
                          / ( (1._dp+3._dp*zknduai_1d(:))                                         &
                            * (2._dp*zkair_1d(:)+zkdu*(5._dp*zknduai_1d(:)+1._dp)))
 
             zftduai_1d(:) = MERGE(0._dp, zftduai_1d(:), ll2_1d(:))

             zknduci_1d(:) = 7.37_dp*ztp1tmp(:)/(2.88e5_dp*papp1(:,jk)*(zrwetci_2d(:,jk)+zeps))
             zknduci_1d(:) = MERGE(1._dp, zknduci_1d(:), ll3_1d(:))
 
             zftduci_1d(:) = 0.4_dp &
                          * ( 1._dp + zknduci_1d(:)*(1.45_dp + 0.4_dp*EXP(-1._dp/zknduci_1d(:)))) &
                          * (zkair_1d(:)+2.5_dp*zknduci_1d(:)*zkdu)                                  &
                          / ( (1._dp+3._dp*zknduci_1d(:))                                         &
                            * (2._dp*zkair_1d(:)+zkdu*(5._dp*zknduci_1d(:)+1._dp)))
 
             zftduci_1d(:) = MERGE(0._dp, zftduci_1d(:), ll3_1d(:))
         ENDIF
!end thermophoresis

         DO jl=1,kproma

            zfracdusol     = MIN(ndusol_strat(jl,jk,jrow)/(pcdncact(jl,jk)      +zeps), 1._dp)
            zfracduinsolai = MIN(nduinsolai(jl,jk,jrow)  /(naerinsol(jl,jk,jrow)+zeps), 1._dp)
            zfracduinsolci = MIN(nduinsolci(jl,jk,jrow)  /(naerinsol(jl,jk,jrow)+zeps), 1._dp)
            zfracbcsol     = MIN(nbcsol_strat(jl,jk,jrow)/(pcdncact(jl,jk)      +zeps), 1._dp)
            zfracbcinsol   = MIN(nbcinsol(jl,jk,jrow)    /(naerinsol(jl,jk,jrow)+zeps), 1._dp)

            zradl = (0.75_dp*zxlb(jl) * zrho(jl,jk)           &
                  / (api*rhoh2o*zcdnc(jl,jk)))**(1._dp/3._dp)

            zf1 = 4._dp*api*zradl*zcdnc(jl,jk)*zrho_rcp(jl,jk)

            zfrzcntdu = MIN(1._dp,MAX(0._dp,-(0.1014_dp*(ztp1tmp(jl)-tmelt)+0.3277_dp)))  ! montmorillonite
            !zfrzcntdu = MIN(1._dp,MAX(0._dp,-(0.1007_dp*(ztp1tmp(jl)-tmelt)+0.6935_dp)))  ! kaolinite
         
!            zfrzcntbc = MIN(1._dp,MAX(0._dp,-(0.0614_dp*(ztp1tmp(jl)-tmelt)+0.5730_dp)))
            zfrzcntbc = 0._dp ! disable BC contact freezing

            zfrzcnt = zxlb(jl) / zcdnc(jl,jk) * zrho(jl,jk) * zf1                                &
                    * ( zfrzcntdu * (zdfarduai_1d(jl)*zfracduinsolai+zdfarduci_1d(jl)*zfracduinsolci)   &
                      + zfrzcntbc *  zdfarbcki_1d(jl)*zfracbcinsol                                 ) &
                    * ( zcdnc(jl,jk)+zicnc(jl,jk) )

            zfrzcnt        = zxlb(jl)*(1._dp-EXP(-zfrzcnt/MAX(zxlb(jl), cqtmin)*ztmst))

!--- ulrike: add thermophoresis (UL, 8.1.09)
            zfrzthermo=0._dp

            IF (lthermo) THEN

               lo2_1d(jl) = (ztp1tmp(jl) < cthomi) .OR.                      &
                            (ztp1tmp(jl) < tmelt .AND. zxip1_1d(jl) > csecfrl &
                            .AND. zsusatw_2d(jl,jk) < zeps)

               zqsm1      = MERGE(tlucua(itm1_look(jl,jk)),tlucuaw(itm1_look(jl,jk)),lo2_1d(jl)) &
                          / papm1(jl,jk)
               zqsm1      = MIN(zqsm1,0.5_dp)
               zqsm1      = zqsm1/(1._dp-vtmpc1*zqsm1)

               it1_1d(jl)  = NINT(ztp1tmp(jl)*1000._dp)
               it1_1d(jl)  = MAX(MIN(it1_1d(jl),jptlucu2),jptlucu1)
               IF (it1_1d(jl) < jptlucu1 .OR. it1_1d(jl) > jptlucu2) lookupoverflow = .TRUE.

               zqst1       = MERGE(tlucua(it1_1d(jl)),tlucuaw(it1_1d(jl)),lo2_1d(jl))/papm1(jl,jk)
               zqst1       = MIN(zqst1,0.5_dp)
               zqst1       = zqst1/(1._dp-vtmpc1*zqst1)
               zdqsdtime   = (zqst1-zqsm1)
 
               zdeltatemp = MAX( (zrho(jl,jk)*alv*(zdep(jl)+zdqsdtime))                &
                                 / (ztmst*zkair_1d(jl)*4._dp*api*zcdnc(jl,jk)*zradl+zeps), 0._dp )

               zf2        = zkair_1d(jl) / papp1(jl,jk) * zdeltatemp

               zfrzthermo = zf1*zf2                                      &
                          * ( zfrzcntbc* zftbcki_1d(jl)*zfracbcinsol     &
                            + zfrzcntdu*(zftduai_1d(jl)*zfracduinsolai   &
                                        +zftduci_1d(jl)*zfracduinsolci)) &
                          * ( zcdnc(jl,jk)+zicnc(jl,jk) )*ztmst          &
                          / ( zcdnc(jl,jk)*zrho(jl,jk) )

               zfrzthermo = zxlb(jl)*(1._dp-EXP(-zfrzthermo))

            ENDIF
!--- Ulrike: end thermophoresis

            znaimmdu  = 32.3_dp*zfracdusol     ! montmorillonite 
            !znaimmdu  = 6.15E-2_dp*zfracdusol  ! kaolinite 

            znaimmbc  = 2.91E-3_dp*zfracbcsol

            zomega = pvervel(jl,jk) - 1.33_dp*SQRT(ptkem1(jl,jk))*zrho(jl,jk)*g
            ztte   = zomega / cpd *zrho_rcp(jl,jk)

            zfrzimm = -(znaimmdu+znaimmbc)*zrho(jl,jk)/rhoh2o*EXP(tmelt-ztp1tmp(jl))*MIN(ztte,0._dp) 
            zfrzimm = zxlb(jl)*(1._dp-EXP(-zfrzimm*zxlb(jl)/zcdnc(jl,jk)*ztmst))

            ztmp1_1d(jl) = zfrzcnt + zfrzimm + zfrzthermo
            ztmp1_1d(jl) = MAX(0.0_dp,MIN(ztmp1_1d(jl),zxlb(jl))) !SF zfrl surrogate

            ztmp2_1d(jl) = zcdnc(jl,jk)*ztmp1_1d(jl)/(zxlb(jl)+zeps)

            ztmp3_1d(jl) = ndusol_strat(jl,jk,jrow) + nbcsol_strat(jl,jk,jrow) + nduinsolai(jl,jk,jrow) &
                         + nduinsolci(jl,jk,jrow) + nbcinsol(jl,jk,jrow) - zicnc(jl,jk)

            ztmp2_1d(jl) = MIN(ztmp2_1d(jl), ztmp3_1d(jl))
            ztmp2_1d(jl) = MAX(ztmp2_1d(jl), 0._dp)  !SF zfrln surrogate

         ENDDO !SF end loop jl

         zfrl(:,jk) = MERGE(ztmp1_1d(:), zfrl(:,jk) , lo_1d(:))

         ztmp1_1d(:) = MIN(ztmp2_1d(:), zcdnc(:,jk)-cdncmin)
         ztmp1_1d(:) = MAX(ztmp1_1d(:), 0._dp)
         zfrln(:)    = MERGE(ztmp1_1d(:), zfrln(:), lo_1d(:))

         ztmp1_1d(:) = zcdnc(:,jk)-zfrln(:)
         ztmp1_1d(:) = MAX(ztmp1_1d(:), cqtmin)
         zcdnc(:,jk) = MERGE(ztmp1_1d(:), zcdnc(:,jk), lo_1d(:))

         ztmp1_1d(:) = zicnc(:,jk)+zfrln(:)
         ztmp1_1d(:) = MAX(ztmp1_1d(:), cqtmin)
         zicnc(:,jk) = MERGE(ztmp1_1d(:), zicnc(:,jk), lo_1d(:))

!--- End included for CDNC/IC scheme -----------------------------------

         ztmp1_1d(:) = zxlb(:)-zfrl(:,jk)
         zxlb(:)     = MERGE(ztmp1_1d(:), zxlb(:), lo_1d(:))

         ztmp1_1d(:) = zxib(:)+zfrl(:,jk)
         zxib(:)     = MERGE(ztmp1_1d(:), zxib(:), lo_1d(:))

         ztmp1_1d(:) = zfrl(:,jk)*zclcaux(:)
         zfrl(:,jk)     = MERGE(ztmp1_1d(:), zfrl(:,jk), lo_1d(:))

!--------corinna: Bergeron-Findeisen-Process:
         ll1_1d(:) =     lo_1d(:)                  &
                   .AND. locc_1d(:)                &
                   .AND. (zdep(:) > 0._dp)         &
                   .AND. (zxlb(:) > 0._dp)         &
                   .AND. (0.01_dp*zvervx_2d(:,jk) < zvervmax_1d(:))      

         ztmp1_1d(:) = ztmst_rcp*zxlb(:)*zclcaux(:) !SF zzevp
      
         ztmp2_1d(:) = pxlte(:,jk)-ztmp1_1d(:)
         pxlte(:,jk) = MERGE(ztmp2_1d(:), pxlte(:,jk), ll1_1d(:))

         ztmp2_1d(:) = pxite(:,jk)+ztmp1_1d(:)
         pxite(:,jk) = MERGE(ztmp2_1d(:), pxite(:,jk), ll1_1d(:))

         ztmp2_1d(:) = ptte(:,jk)+(zlsdcp(:)-zlvdcp(:))*ztmp1_1d(:)
         ptte(:,jk)  = MERGE(ztmp2_1d(:), ptte(:,jk), ll1_1d(:))
          
         zcdnc(:,jk) = MERGE(cqtmin, zcdnc(:,jk), ll1_1d(:))

         ztmp2_1d(:) = zxib(:)+zxlb(:)
         zxib(:)     = MERGE(ztmp2_1d(:), zxib(:), ll1_1d(:))

         zxlb(:) = MERGE(0._dp, zxlb(:), ll1_1d(:))

!-------End corinna: Bergeron-Findeisen-Process

!!$        IF (lookupoverflow) CALL lookuperror ('cloud (2)    ') ! mz_ht_20120321
        IF (lookupoverflow) THEN
          status_string = 'lookuperror: cdnc - cloud (1)'
          RETURN
        ENDIF
!
!     ------------------------------------------------------------------
!       7.  Cloud physics and precipitation fluxes at the surface
!
!ham_ps: cdir circumvents bug in sxf90 compiler

        zclcstar_1d(:) = MIN(zclcaux(:), zclcpre(:))
        zauloc_1d(:)   = 3./5000._dp*zdz_2d(:,jk)
        zauloc_1d(:)   = MAX(MIN(zauloc_1d(:), clmax), clmin)

        ll1_1d(:) = (knvb(:) >= jbmin) .AND. &
                    (knvb(:) <= jbmax) .AND. &
                    (pvervel(:,jk) > 0._dp)

        ll2_1d(:) = (jk == knvb(:)  ) .OR. &
                    (jk == knvb(:)+1)

        ll3_1d(:) = ll1_1d(:) .AND. ll2_1d(:) .AND. lonacc

        zauloc_1d(:) = MERGE(0._dp, zauloc_1d(:), ll3_1d(:))

        zxlb(:) = MAX(zxlb(:),1.e-20_dp)
        zxib(:) = MAX(zxib(:), 1.e-20_dp)

! liquid water and snow content are stored 
! before the reduction by outfalling rain
! (necessary for nucleation scavenging)
        plwc(:,jk) = zxlb(:)
        piwc(:,jk) = zxib(:)      

!---Included for in-cloud scavenging (Philip Stier, 25/11/03):----------
        zmlwc(:,jk) = zxlb(:)
        zmiwc(:,jk) = zxib(:)

!---End Included for scavenging-----------------------------------------
!
!---  Calculate the rain and snow water content in kg/kg from the rain and snow flux
!
       ll1_1d(:) = (zclcpre(:) > zeps )
       ll2_1d(:) = ll1_1d(:) .AND. (zrfl(:) > cqtmin)
       ll3_1d(:) = ll1_1d(:) .AND. (zsfl(:) > cqtmin)

       ztmp1_1d(:) = ( MAX(zrfl(:), cqtmin)/(12.45_dp*MAX(zclcpre(:),zeps)*SQRT(zqrho_2d(:,jk))) )**(8._dp/9._dp)
       ztmp2_1d(:) = ( MAX(zsfl(:), cqtmin)/(cvtfall *MAX(zclcpre(:),zeps)                     ) )**(1._dp/1.16_dp)

       zxrp1_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll2_1d(:))
       zxsp1_1d(:) = MERGE(ztmp2_1d(:), 0._dp, ll3_1d(:))

!
!       7.1   Warm clouds: Coalescence processes after Beheng (1994):
!             Autoconversion of cloud droplets and collection of cloud
!             droplets by falling rain. Accretion of cloud droplets by
!             falling snow (zsacl) is calculated under 7.2
!

        ll1_1d(:) = locc_1d(:)                .AND. &
                    (zxlb(:) > cqtmin)        .AND. &
                    (zcdnc(:,jk) >= cdncmin)

        IF (nauto == 2) THEN
!          Autoconversion rate from Khairoutdinov and Kogan, 2000

           ztmp1_1d(:) = ccraut*1350._dp*(1.e-6_dp*zcdnc(:,jk))**(-1.79_dp)
     
           ztmp1_1d(:) = zxlb(:) * (  1._dp &
                                   - (1._dp + ztmst*zexm1_1*ztmp1_1d(:)*zxlb(:)**zexm1_1)**zexp_1)

           ztmp1_1d(:) = MIN(zxlb(:), ztmp1_1d(:))
           zraut_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))
           
           ztmp1_1d(:) = zxlb(:) - zraut_1d(:)
           ztmp2_1d(:) = zxlb(:) !SF keeps zxlb for later use
           zxlb(:)     = MERGE(ztmp1_1d(:), zxlb(:), ll1_1d(:))
 
!--- zrac1 is formed by accretion with rain from above
!--- zrac2 is formed by accretion with newly formed rain inside the grid box

           ztmp1_1d(:) = -3.7_dp*ztmst*zxrp1_1d(:)
           ztmp1_1d(:) = EXP(ztmp1_1d(:))
           ztmp1_1d(:) = zxlb(:)*(1._dp-ztmp1_1d(:))
           zrac1_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))

           zxlb(:) = zxlb(:) - zrac1_1d(:)

           ztmp1_1d(:) = -3.7_dp*ztmst*zauloc_1d(:)*zrho(:,jk)*zraut_1d(:)
           ztmp1_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))
           ztmp1_1d(:) = zxlb(:)*(1._dp-EXP(ztmp1_1d(:)))
           zrac2_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))

           zxlb(:) = zxlb(:) - zrac2_1d(:)

           zrpr(:) = zrpr(:) + zclcaux(:)     * (zraut_1d(:)+zrac2_1d(:)) &
                             + zclcstar_1d(:) *  zrac1_1d(:)

!---Included for in-cloud scavenging (Philip Stier, 26/11/03):----------
           ztmp1_1d(:)    = zraut_1d(:)+zrac1_1d(:)+zrac2_1d(:)
           zmratepr(:,jk) = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))
!---End Included for scavenging-----------------------------------------

!--- Included for prognostic CDNC/IC scheme ----------------------------
!--- Autoconversion also changes the number of cloud droplets (zrprn)

           ztmp1_1d(:) = (zraut_1d(:)+zrac1_1d(:)+zrac2_1d(:))/(ztmp2_1d(:)+zeps) !SF ztmp2_1d=zxlb+zraut+zrac1+zrac2
           zrprn(:)    = MERGE(ztmp1_1d(:), zrprn(:), ll1_1d(:))

           ll2_1d(:) = ll1_1d(:)          .AND. &
                       (zxlb(:) > cqtmin)

           ztmp1_1d(:) = MERGE(cdncmin, 0._dp, ll2_1d(:))
           ztmp1_1d(:) = zcdnc(:,jk)-ztmp1_1d(:)
           ztmp2_1d(:) = zcdnc(:,jk)*zrprn(:)

           ztmp3_1d(:) = MIN(ztmp1_1d(:), ztmp2_1d(:))
           zrprn(:)    = MERGE(ztmp3_1d(:), zrprn(:), ll1_1d(:))

           ztmp1_1d(:) = zcdnc(:,jk)-zrprn(:)
           ztmp1_1d(:) = MAX(ztmp1_1d(:), cqtmin)
           zcdnc(:,jk) = MERGE(ztmp1_1d(:), zcdnc(:,jk), ll1_1d(:))

!--- End included for CDNC/IC scheme ------------------------------------
!--- End included alternative autoconversion parameterisation ----------
!--- Changed for alternative autoconversion parameterisation -----------

        ELSE   !SF( nauto ==1)
!          Beheng (1994) - ECHAM 5 standard

!---Changed for prognostic CDNC/IC scheme ------------------------------
!   (Replaced pacdnc by zcdnc which is set above)
           ztmp1_1d(:) = ccraut*1.2e27_dp * zrho_rcp(:,jk)                    & 
                                          * (zcdnc(:,jk)*1.e-6_dp)**(-3.3_dp) &
                                          * (zrho(:,jk)*1.e-3_dp)**4.7_dp
           zraut_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))
!--- End changed for prognostic CDNC/IC scheme -------------------------

           zraut_1d(:) = zxlb(:) * (  1._dp &
                                   - (1._dp + ztmst*zexm1_2*zraut_1d(:)*zxlb(:)**zexm1_2)**zexp_2)

           zraut_1d(:) = MIN(zxlb(:), zraut_1d(:))
  
!--- Included for prognostic CDNC/IC scheme ----------------------------
           ztmp1_1d(:) = 7.7e9_dp * zraut_1d(:) * zrho(:,jk)                      !SF zrautn
           ztmp2_1d(:) = 1.289e10_dp * 1.e-6_dp * ztmst * (zrho(:,jk)*zxlb(:))**2 !SF zself (1.e-6 comes
                                                                                  ! from a unit bug fix)
           ztmp3_1d(:) = ztmp1_1d(:) + ztmp2_1d(:)
           ztmp3_1d(:) = MIN(ztmp3_1d(:),zcdnc(:,jk))
           zrautself_1d(:) = MERGE(ztmp3_1d(:), 0._dp, ll1_1d(:))

           ztmp1_1d(:) = zcdnc(:,jk)-zrautself_1d(:)
           ztmp1_1d(:) = MAX(ztmp1_1d(:), cqtmin)
           zcdnc(:,jk) = MERGE(ztmp1_1d(:), zcdnc(:,jk), ll1_1d(:))

!--- End included for CDNC/IC scheme -----------------------------------

           ztmp1_1d(:) = zxlb(:) - zraut_1d(:)
           zxlb(:)     = MERGE(ztmp1_1d(:), zxlb(:), ll1_1d(:))
!
!--- zrac1 is formed by accretion with rain from above
!--- zrac2 is formed by accretion with newly formed rain inside the grid box
!
           zrac1_1d(:) = -6._dp*ztmst*zxrp1_1d(:)
           zrac1_1d(:) = EXP(zrac1_1d(:))
           zrac1_1d(:) = zxlb(:)*(1._dp-zrac1_1d(:))

           ztmp1_1d(:) = zxlb(:) - zrac1_1d(:)
           ztmp2_1d(:) = zxlb(:)  !SF keeps zxlb for later use
           zxlb(:)     = MERGE(ztmp1_1d(:), zxlb(:), ll1_1d(:))
 
           zrac2_1d(:) = -6._dp*ztmst*zauloc_1d(:)*zrho(:,jk)*zraut_1d(:)
           zrac2_1d(:) = EXP(zrac2_1d(:))
           zrac2_1d(:) = zxlb(:)*(1._dp-zrac2_1d(:))

           ztmp1_1d(:) = zxlb(:)-zrac2_1d(:)
           zxlb(:)     = MERGE(ztmp1_1d(:), zxlb(:), ll1_1d(:))

           ztmp1_1d(:) = zrpr(:) + zclcaux(:)     * (zraut_1d(:)+zrac2_1d(:)) & 
                                 + zclcstar_1d(:) * zrac1_1d(:)
           zrpr(:)     = MERGE(ztmp1_1d(:), zrpr(:), ll1_1d(:))

!---Included for in-cloud scavenging (Philip Stier, 26/11/03):----------
           ztmp1_1d(:)    = zraut_1d(:) + zrac1_1d(:) + zrac2_1d(:)
           zmratepr(:,jk) = MERGE(ztmp1_1d(:), zmratepr(:,jk), ll1_1d(:)) 
!---End Included for scavenging-----------------------------------------

!--- Included for prognostic CDNC/IC scheme ----------------------------
!--- Autoconversion also changes the number of cloud droplets (zrprn)
           ztmp1_1d(:) = (zrac1_1d(:)+zrac2_1d(:))/(ztmp2_1d(:)+zeps) !SF ztmp2_1d=zxlb+zrac1+zrac2
           ztmp1_1d(:) = zcdnc(:,jk)*ztmp1_1d(:)
           ztmp1_1d(:) = MIN(ztmp1_1d(:), zcdnc(:,jk)) !SF zraccn

           ztmp2_1d(:) = zcdnc(:,jk)-ztmp1_1d(:)
           ztmp2_1d(:) = MAX(ztmp2_1d(:), cqtmin)
           zcdnc(:,jk) = MERGE(ztmp2_1d(:), zcdnc(:,jk), ll1_1d(:))

           ztmp2_1d(:) = zrautself_1d(:) + ztmp1_1d(:)
           zrprn(:)    = MERGE(ztmp2_1d(:), zrprn(:), ll1_1d(:))
!--- End included for CDNC/IC scheme -----------------------------------
        ENDIF 

!       7.2  Cold clouds:
!            Conversion of cloud ice to snow after Levkov et al. 1992:
!            Aggregation of ice crystals to snow assuming plates (zsaut) and accretion of ice
!            by falling snow. (zsaci)
!            Accrection of cloud droplets by falling snow. (zsacl)
!            Effective radius of ice crystals assuming plates (Lohmann, ACPD, 2007)

        ll1_1d(:) = locc_1d(:) .AND. (zxib(:) > cqtmin)

        ztmp1_1d(:) = 0.5e4_dp * ( 1000._dp/0.0376_dp*zxib(:)*zrho(:,jk)/zicnc(:,jk) )**0.302 !plate !SF zrieff
        ztmp1_1d(:) = MIN(MAX(ztmp1_1d(:), ceffmin), ceffmax) !SF zrieff

        ztmp2_1d(:) = 5113188._dp+2809._dp*ztmp1_1d(:)**3
        ztmp2_1d(:) = SQRT(ztmp2_1d(:))
        ztmp2_1d(:) = -2261._dp + ztmp2_1d(:)  !SF zrih
        
        ztmp3_1d(:) = 1.e-6_dp*ztmp2_1d(:)**(1._dp/3._dp)
        zris_1d(:)  = MERGE(ztmp3_1d(:), 1._dp, ll1_1d(:))   !SF 1. could be whatever, just not 0.
       
!--- temperature dependent collision efficiency (zcolleffi)

        ztmp1_1d(:)    = 0.025_dp*(ztp1tmp(:)-tmelt)
        ztmp1_1d(:)    = EXP(ztmp1_1d(:))
        zcolleffi_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:)) 

        zc1_1d(:) = 17.5_dp / crhoi * zrho(:,jk) * zqrho_2d(:,jk)**0.33_dp

        ztmp1_1d(:) = -6._dp / zc1_1d(:) * LOG10(1.e4_dp*zris_1d(:)) !SF zdt2
        ztmp1_1d(:) = ccsaut / ztmp1_1d(:)

        zsaut_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))
        zsaut_1d(:) = zxib(:) * (1._dp - 1._dp/(1._dp+zsaut_1d(:)*ztmst*zxib(:)))

        ztmp1_1d(:) = zxib(:) - zsaut_1d(:)
        zxibold_1d(:) = zxib(:) !SF store zxib for later use at the end of 7.2
        zxib(:)     = MERGE(ztmp1_1d(:), zxib(:), ll1_1d(:))

        zsaci2_1d(:) = 0.0_dp
        zsacl2_1d(:) = 0.0_dp

!--- Included for prognostic CDNC/IC scheme ----------------------------
!    accretion also affects the number of cloud droplets
        zsacl2in_1d(:)  = 0.0_dp      ! 
!--- End included for CDNC/IC scheme -----------------------------------

        ztmp1_1d(:) = zauloc_1d(:)*zrho(:,jk)*zsaut_1d(:)
        zxsp2_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))  !snow that is formed inside the grid box

        zxsp_1d(:) = zxsp1_1d(:) + zxsp2_1d(:)

!--- CCMod -------------------------------------------------------------

        ll_bcc(:) = zcccov_ges(:,jk) .gt. 1.e-10_dp 

        !Precipitation for contrails

        !Aggretion
        ll1_1d(:) = ll_bcc(:) .AND. (zxibcc(:) > cqtmin)

        ztmp1_1d(:) = 0.5e4_dp * ( 1000._dp/0.0376_dp*zxibcc(:)*zrho(:,jk)/zccicnc(:,jk) )**0.302 !plate !SF zrieff
        ztmp1_1d(:) = MIN(MAX(ztmp1_1d(:), ceffmin), ceffmax) !SF zrieff

        ztmp2_1d(:) = 5113188._dp+2809._dp*ztmp1_1d(:)**3
        ztmp2_1d(:) = SQRT(ztmp2_1d(:))
        ztmp2_1d(:) = -2261._dp + ztmp2_1d(:)  !SF zrih

        ztmp3_1d(:) = 1.e-6_dp*ztmp2_1d(:)**(1._dp/3._dp)
        zris_1d(:)  = MERGE(ztmp3_1d(:), 1._dp, ll1_1d(:))   !SF 1. could be whatever, just not 0.

!--- temperature dependent collision efficiency (zcolleffi)

        ztmp1_1d(:)    = 0.025_dp*(ztp1tmp(:)-tmelt)
        ztmp1_1d(:)    = EXP(ztmp1_1d(:))
        zcolleffi_cc(:) = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))

        zc1_cc(:) = 17.5_dp / crhoi * zrho(:,jk) * zqrho_2d(:,jk)**0.33_dp

        ztmp1_1d(:) = -6._dp / zc1_cc(:) * LOG10(1.e4_dp*zris_1d(:)) !SF zdt2
        ztmp1_1d(:) = ccsaut / ztmp1_1d(:)

        zsaut_cc(:) = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))
        zsaut_cc(:) = zxibcc(:) * (1._dp - 1._dp/(1._dp+zsaut_cc(:)*ztmst*zxibcc(:)))

        ztmp1_1d(:) = zxibcc(:) - zsaut_cc(:)
        zxibold_cc(:) = zxibcc(:) !SF store zxib for later use at the end of 7.2
        zxibcc(:)     = MERGE(ztmp1_1d(:), zxibcc(:), ll1_1d(:))

        ztmp1_1d(:) = zauloc_1d(:)*zrho(:,jk)*zsaut_cc(:)
        zxsp2_1d(:) = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))  !snow that is formed inside the grid box

        zxsp_1d(:) = zxsp_1d(:) + zxsp2_1d(:)

!--- End CCMod -------------------------------------------------------------

        ll2_1d(:) = ll1_1d(:)               .AND. &
                    (zxsp_1d(:)  >  cqtmin) .AND. &
                    (zxlb(:)     >  cqtmin) .AND. &
                    (zcdnc(:,jk) >= cdncmin)

!--- The riming of snow with droplets is calculated assuming planar snow flakes (Lohmann, JAS, 2004)
!--- It depends on the droplet (zudrop) and snow flake (zusnow) fall velocity, the Stokes
!--- number (zstokes) and the Reynolds number (zrey)

        zscnc_1d(:) = 0.5_dp*zicnc(:,jk)  !SF note: not put in the ll2_1d() condition as in orig code
        IF (jk > 1) THEN 
           zscnc_1d(:) = MAX(zicemin, MIN(zsprn(:,jk-1), zscnc_1d(:)))
        ELSE
           zscnc_1d(:) = MAX(zicemin, zscnc_1d(:)) !SF added protection for non-zero 
        ENDIF

        ztmp1_1d(:) = ( 6._dp*zpirho_rcp*zrho(:,jk)*zxlb(:)/zcdnc(:,jk) )**(1._dp/3._dp)
        ztmp1_1d(:) = MAX(ztmp1_1d(:), 1.e-6_dp) !SF zdw

        zudrop_1d(:) = 1.19e4_dp*2500._dp*ztmp1_1d(:)**2*(1.3_dp*zrho_rcp(:,jk))**0.35_dp

        zdplanar_1d(:) = 1.e3_dp/3.8e-4_dp*zxsp_1d(:)/zscnc_1d(:)
        zdplanar_1d(:) = SQRT(zdplanar_1d(:))
        zdplanar_1d(:) = 1.e-2_dp*zdplanar_1d(:)
        zdplanar_1d(:) = MAX(20.e-6_dp, zdplanar_1d(:)) !SF note: not put in the ll2_1d() condition as in orig code (useless)
     
        zusnow_1d(:) = 2.34_dp * (100._dp*zdplanar_1d(:))**0.3_dp &
                               * (1.3_dp*zrho_rcp(:,jk))**0.35_dp

        zstokes_1d(:) = 2._dp*g_rcp*(zusnow_1d(:)-zudrop_1d(:))*zudrop_1d(:)/zdplanar_1d(:)
        zstokes_1d(:) = MAX(zstokes_1d(:), cqtmin)

        zrey_1d(:) = zrho(:,jk)*zdplanar_1d(:)*zusnow_1d(:)/zviscos_2d(:,jk)
        zrey_1d(:) = MAX(zrey_1d(:),cqtmin)

        ll3_1d(:) = (zrey_1d(:) <=  5._dp)
        ll4_1d(:) = (zrey_1d(:) >   5._dp) .AND. (zrey_1d(:) <  40._dp)
        ll5_1d(:) = (zrey_1d(:) >= 40._dp)
     
        ztmp1_1d(:)   = 5.52_dp*zrey_1d(:)**(-1.12_dp)
        ztmp2_1d(:)   = 1.53_dp*zrey_1d(:)**(-0.325_dp)
        zstcrit_1d(:) = 1._dp
        zstcrit_1d(:) = MERGE(ztmp1_1d(:), zstcrit_1d(:), ll3_1d(:))
        zstcrit_1d(:) = MERGE(ztmp2_1d(:), zstcrit_1d(:), ll4_1d(:))

        zcsacl_1d(:) = 0.2_dp * ( LOG10(zstokes_1d(:)) - LOG10(zstcrit_1d(:)) - 2.236_dp )**2
        zcsacl_1d(:) = MIN(zcsacl_1d(:), 1._dp-cqtmin)
        zcsacl_1d(:) = MAX(zcsacl_1d(:), 0._dp)
        zcsacl_1d(:) = SQRT(1._dp - zcsacl_1d(:))

        ll6_1d(:) = ll5_1d(:) .AND. (zstokes_1d(:) <= 0.06_dp)
        ll7_1d(:) = ll5_1d(:) .AND. (zstokes_1d(:) >  0.06_dp) .AND. (zstokes_1d(:) <= 0.25_dp)
        ll8_1d(:) = ll5_1d(:) .AND. (zstokes_1d(:) >  0.25_dp) .AND. (zstokes_1d(:) <= 1.00_dp)

        WHERE (ll6_1d(:))
              zcsacl_1d(:) = 1.034_dp*zstokes_1d(:)**1.085_dp
        ELSEWHERE (ll7_1d(:))
              zcsacl_1d(:) = 0.787_dp*zstokes_1d(:)**0.988_dp
        ELSEWHERE (ll8_1d(:))
              zcsacl_1d(:) = 0.7475_dp*LOG10(zstokes_1d(:))+0.65_dp
        ELSEWHERE (ll5_1d(:))
              zcsacl_1d(:) = (zstokes_1d(:)+1.1_dp)**2/(zstokes_1d(:)+1.6_dp)**2
        ENDWHERE

        zcsacl_1d(:) = MAX(MIN(zcsacl_1d(:), 1._dp), 0.01_dp)
        zcsacl_1d(:)  = MERGE(zcsacl_1d(:), 0._dp, ll2_1d(:))

        ztmp1_1d(:)  = zcons4*zxsp_1d(:)**0.8125_dp !SF zlamsm
        ztmp1_1d(:)  = api*cn0s*3.078_dp*ztmp1_1d(:)*zqrho_2d(:,jk)**0.5_dp
        zsaci2_1d(:) = MERGE(ztmp1_1d(:), zsaci2_1d(:), ll2_1d(:)) 

        ztmp1_1d(:)  = -ztmst*zsaci2_1d(:)*zcsacl_1d(:)
        ztmp1_1d(:)  = EXP(ztmp1_1d(:))
        ztmp1_1d(:)  = zxlb(:)*(1._dp-ztmp1_1d(:))

!--- Included for prognostic CDNC/IC scheme ----------------------------
        zsacl2in_1d(:) = MERGE(ztmp1_1d(:), zsacl2in_1d(:), ll2_1d(:))
!--- End included for CDNC/IC scheme -----------------------------------

        ztmp2_1d(:) = zxlb(:)-ztmp1_1d(:)
        ztmp3_1d(:) = zxlb(:)  !SF keeps zxlb for later use
        zxlb(:)     = MERGE(ztmp2_1d(:), zxlb(:), ll2_1d(:))

        ztmp2_1d(:) = zclcaux(:)*ztmp1_1d(:)
        zsacl2_1d(:)   = MERGE(ztmp2_1d(:), zsacl2_1d(:), ll2_1d(:))

!SF end ll2_1d condition:

        ll2_1d(:) = ll1_1d(:)             .AND. &
                    (zxsp_1d(:) > cqtmin) .AND. &
                    (zxib(:)    > cqtmin)

        ztmp1_1d(:) = zcons4*zxsp_1d(:)**0.8125_dp
        ztmp1_1d(:) = api*cn0s*3.078_dp*ztmp1_1d(:)*zqrho_2d(:,jk)**0.5_dp
        ztmp1_1d(:) = -ztmst*ztmp1_1d(:)*zcolleffi_1d(:)
        ztmp1_1d(:) = EXP(ztmp1_1d(:))
        ztmp1_1d(:) = zxib(:) * (1._dp-ztmp1_1d(:))

        zsaci2_1d(:) = MERGE(ztmp1_1d(:), zsaci2_1d(:), ll2_1d(:))
        
        zxib(:) = zxib(:)-zsaci2_1d(:)

!SF end ll2_1d condition:

        zsacl(:) = MERGE(zsacl2_1d(:), zsacl(:), ll1_1d(:))

        ztmp1_1d(:) = zspr(:) + zclcaux(:)*( zsaut_1d(:)+zsaci2_1d(:) )
        zspr(:)     = MERGE(ztmp1_1d(:), zspr(:), ll1_1d(:))

!---Included for in-cloud scavenging (Philip Stier, 25/11/03):----------
        ztmp1_1d(:)    = zsaut_1d(:)+zsaci2_1d(:)
        zmrateps(:,jk) = MERGE(ztmp1_1d(:), zmrateps(:,jk), ll1_1d(:))
!---End Included for scavenging-----------------------------------------

!--- Included for prognostic CDNC/IC scheme ----------------------------
        ll2_1d(:) = (zxlb(:) > cqtmin) 

        ztmp1_1d(:) = zcdnc(:,jk)*zsacl2in_1d(:)/(ztmp3_1d(:)+zeps) !SF ztmp3_1d =zxlb+zsacl2in
        ztmp1_1d(:) = MIN(ztmp1_1d(:), zcdnc(:,jk)-cdncmin)
        ztmp1_1d(:) = MAX(ztmp1_1d(:), 0._dp) 

        ztmp2_1d(:) = zcdnc(:,jk)*zsacl2in_1d(:)/(ztmp3_1d(:)+zeps) !SF ztmp3_1d = zxlb+zsacl2in
        ztmp2_1d(:) = MIN(ztmp2_1d(:), zcdnc(:,jk))

        ztmp3_1d(:)  = MERGE(ztmp1_1d(:), ztmp2_1d(:), ll2_1d(:))
        zsacln(:) = MERGE(ztmp3_1d(:), zsacln(:), ll1_1d(:))

!SF end ll2_1d condition (zxlb > cqtmin)

        ztmp1_1d(:)     = zcdnc(:,jk)-zsacln(:)
        zcdnc(:,jk)     = MERGE(ztmp1_1d(:), zcdnc(:,jk), ll1_1d(:))
        zmsnowacl(:,jk) = MERGE(zsacl2in_1d(:), zmsnowacl(:,jk), ll1_1d(:))

!       secondary ice crystal production (zsecprod) after Levkov et al. 1992
!       sink for snow, source for ice crystals
!uls    included size depedent accretion rate (Lohmann, JAS, 2004)

        zsecprod_1d(:) = 0._dp

        ll2_1d(:) = ll1_1d(:)               .AND. &
                    (zxsp_1d(:) > zepsec)   .AND. &
                    (zxlb(:)    > zepsec)   .AND. &
                    (ztp1tmp(:) > 265.2_dp) .AND. &
                    (ztp1tmp(:) < 270.2_dp)


        ztmp1_1d(:) = ( 6._dp*zpirho_rcp*zrho(:,jk)*zxlb(:)/zcdnc(:,jk) )**(1._dp/3._dp)
        ztmp1_1d(:) = MAX(ztmp1_1d(:), 1.e-6_dp) !SF zdw

        zudrop_1d(:) = 1.19e4_dp*(50._dp*ztmp1_1d(:))**2*(1.3_dp*zrho_rcp(:,jk))**0.35_dp
 
        zstokes_1d(:) = 2._dp*g_rcp*(zusnow_1d(:)-zudrop_1d(:))*zudrop_1d(:)/zdplanar_1d(:)
        zstokes_1d(:) = MAX(zstokes_1d(:), cqtmin)

        ztmp1_1d(:) = 0.2_dp * ( LOG10(zstokes_1d(:)) - LOG10(zstcrit_1d(:)) - 2.236_dp )**2
        ztmp1_1d(:) = MIN(ztmp1_1d(:), 1._dp-cqtmin)
        ztmp1_1d(:) = MAX(ztmp1_1d(:), 0._dp)
        ztmp1_1d(:) = SQRT(1._dp - ztmp1_1d(:))

        ll6_1d(:) = ll5_1d(:) .AND. (zstokes_1d(:) <= 0.06_dp)
        ll7_1d(:) = ll5_1d(:) .AND. (zstokes_1d(:) >  0.06_dp) .AND. (zstokes_1d(:) <= 0.25_dp)
        ll8_1d(:) = ll5_1d(:) .AND. (zstokes_1d(:) >  0.25_dp) .AND. (zstokes_1d(:) <= 1.00_dp)

        WHERE (ll6_1d(:))
              ztmp1_1d(:) = 1.034_dp*zstokes_1d(:)**1.085_dp
        ELSEWHERE (ll7_1d(:))
              ztmp1_1d(:) = 0.787_dp*zstokes_1d(:)**0.988_dp
        ELSEWHERE (ll8_1d(:))
              ztmp1_1d(:) = 0.7475_dp*LOG10(zstokes_1d(:))+0.65_dp
        ELSEWHERE (ll5_1d(:))
              ztmp1_1d(:) = (zstokes_1d(:)+1.1_dp)**2/(zstokes_1d(:)+1.6_dp)**2
        ENDWHERE

        ztmp1_1d(:)  = MAX(MIN(ztmp1_1d(:), 1._dp), 0.01_dp)
        zcsacl_1d(:) = MERGE(ztmp1_1d(:), zcsacl_1d(:), ll2_1d(:))

        ztmp1_1d(:) = zcons5*zxsp_1d(:)**0.875_dp !SF zlams2

        ztmp2_1d(:) = cn0s * 0.831_dp * api / zmw0                      &
                    * zcsacl_1d(:) * zrho(:,jk) * zxlb(:) * ztmp1_1d(:) &
                    * ( g*crhosno / (0.75_dp*zcdi*zrho(:,jk)) )**0.5_dp   !SF zj
    
        ztmp2_1d(:) = MAX(0.00285_dp*ztmp2_1d(:), 0._dp)  !SF zpn

        ztmp2_1d(:) = ztmst*zmi0*ztmp2_1d(:)*zrho_rcp(:,jk)
        ztmp3_1d(:) = zxsp_1d(:)*zrho_rcp(:,jk)
        ztmp2_1d(:) = MIN(ztmp3_1d(:), ztmp2_1d(:))
        ztmp2_1d(:) = MAX(ztmp2_1d(:), 0._dp)

        zsecprod_1d(:) = MERGE(ztmp2_1d(:), zsecprod_1d(:), ll2_1d(:))

        ztmp1_1d(:) = zxib(:)+zsecprod_1d(:)
        zxib(:)     = MERGE(ztmp1_1d(:), zxib(:), ll2_1d(:))

        ztmp1_1d(:) = zspr(:)-zclcstar_1d(:)*zsecprod_1d(:)
        ztmp1_1d(:) = MAX(ztmp1_1d(:), 0._dp)
        zspr(:)     = MERGE(ztmp1_1d(:), zspr(:), ll2_1d(:))
      
        ztmp1_1d(:)    = zmrateps(:,jk)-zsecprod_1d(:)
        zmrateps(:,jk) = MERGE(ztmp1_1d(:), zmrateps(:,jk), ll2_1d(:))

        ! storing of snow production before it is transformed into a flux
        prate_s(:,jk) = zsaut_1d(:) + zsaci2_1d(:)

!SF end ll2_1d condition (secondary ice production)

!SF end ll1_1d condition (locc and zxib(jl) > cqtmin) 
!
!--- Also change the number of ice crystals due to the break-up of snow flakes
!
        ll1_1d(:) = locc_1d(:)               .AND. &
                    (zxib(:)     > zepsec  ) .AND. &
                    (zicnc(:,jk) >= zicemin)
       
        ztmp1_1d(:) = zxibold_1d(:)
        ztmp1_1d(:) = MAX(ztmp1_1d(:), 0._dp) !SF zxibold

        ztmp2_1d(:) = zicnc(:,jk) * (zsaci2_1d(:)+zsaut_1d(:)) / (ztmp1_1d(:)+zeps) !SF zsprn1

        ztmp3_1d(:) = 0.5_dp * ztmst * zc1_1d(:) * zicnc(:,jk) * zxib(:) !SF zself

        ztmp4_1d(:) = zmi0_rcp * zrho(:,jk) * zsecprod_1d(:) !SF zsecprodn

        ztmp5_1d(:) = ztmp2_1d(:)+ztmp3_1d(:)-ztmp4_1d(:) 
        ztmp5_1d(:) = MIN(ztmp5_1d(:), zicnc(:,jk))    !SF zsprnself
        zsprn(:,jk) = MERGE(ztmp5_1d(:), zsprn(:,jk), ll1_1d(:))

        ztmp1_1d(:) = zicnc(:,jk)-zsprn(:,jk)
        ztmp1_1d(:) = MAX(ztmp1_1d(:), cqtmin)
        zicnc(:,jk) = MERGE(ztmp1_1d(:), zicnc(:,jk), ll1_1d(:))

!--- CCMod -------------------------------------------------------------

        !Precipitation for contrails

        !Accretion

        ll1_1d(:) = ll_bcc(:) .AND. (zxibcc(:) > cqtmin)
        ll2_1d(:) = ll1_1d(:) .AND. (zxsp_1d(:) > cqtmin) 

        zsaci2_1d(:) = 0.0_dp

        ztmp1_1d(:) = zcons4*zxsp_1d(:)**0.8125_dp
        ztmp1_1d(:) = api*cn0s*3.078_dp*ztmp1_1d(:)*zqrho_2d(:,jk)**0.5_dp
        ztmp1_1d(:) = -ztmst*ztmp1_1d(:)*zcolleffi_cc(:)
        ztmp1_1d(:) = EXP(ztmp1_1d(:))
        ztmp1_1d(:) = zxibcc(:) * (1._dp-ztmp1_1d(:))

        zsaci2_1d(:) = MERGE(ztmp1_1d(:), zsaci2_1d(:), ll2_1d(:))
        
        zxibcc(:) = zxibcc(:)-zsaci2_1d(:)

        ztmp1_1d(:) = zspr_cc(:) + zccvol_ges(:,jk)*(zsaut_cc(:)+zsaci2_1d(:))
        zspr_cc(:)  = MERGE(ztmp1_1d(:), zspr_cc(:), ll1_1d(:))

!--- Also change the number of ice crystals due to the break-up of snow flakes
!
        ll1_1d(:) = ll_bcc(:)               .AND. &
                    (zxibcc(:)     > zepsec  ) .AND. &
                    (zccicnc(:,jk) >= zicemin)
       
        ztmp1_1d(:) = zxibold_cc(:)
        ztmp1_1d(:) = MAX(ztmp1_1d(:), 0._dp) !SF zxibold

        ztmp2_1d(:) = zccicnc(:,jk) * (zsaci2_1d(:)+zsaut_cc(:)) / (ztmp1_1d(:)+zeps) !SF zsprn1

        ztmp3_1d(:) = 0.5_dp * ztmst * zc1_cc(:) * zccicnc(:,jk) * zxibcc(:) !SF zself

        ztmp5_1d(:) = ztmp2_1d(:)+ztmp3_1d(:) 
        ztmp5_1d(:) = MIN(ztmp5_1d(:), zccicnc(:,jk))    !SF zsprnself
        zsprn_cc(:,jk) = MERGE(ztmp5_1d(:), zsprn_cc(:,jk), ll1_1d(:))

        ztmp1_1d(:) = zccicnc(:,jk)-zsprn_cc(:,jk)
        ztmp1_1d(:) = MAX(ztmp1_1d(:), cqtmin)
        zccicnc(:,jk) = MERGE(ztmp1_1d(:), zccicnc(:,jk), ll1_1d(:))

!--- End CCMod -------------------------------------------------------------

!--- End included for CDNC/IC scheme -----------------------------------

!       7.3 Updating precipitation fluxes. In the lowest layer (klev),
!           the sedimentation sink of cloud ice is balanced
!           by precipitation at the surface (through 'zzdrs').
!           Fraction of precipitating clouds (zclcpre) used for the
!           calculation of evaporation/sublimation of rain/snow in
!           the next layer
!
        zzdrr_1d(:)    = zcons2*zdp_2d(:,jk)*zrpr(:)
!--- CCMod -------------------------------------------------------------

        zzdrs_1d(:)    = zcons2*zdp_2d(:,jk)*(zspr(:)+zspr_cc(:)+zsacl(:))

!--- End CCMod -------------------------------------------------------------

        IF (jk .EQ. klev) THEN
!--- CCMod -------------------------------------------------------------

           zzdrs_1d(:) = zzdrs_1d(:)+zxiflux(:)+zxiflux_cc(:)

!--- End CCMod -------------------------------------------------------------

           ztmp1_1d(:) = zcons2*zdp_2d(:,jk)/(zlsdcp(:)-zlvdcp(:)) &
                               *MAX(0._dp, (ztp1tmp(:)-tmelt))       !SF zcons
           ztmp2_1d(:) = MIN(zxsec*zzdrs_1d(:), ztmp1_1d(:))         !SF zsnmlt
           zzdrr_1d(:) = zzdrr_1d(:)+ztmp2_1d(:)
           zzdrs_1d(:) = zzdrs_1d(:)-ztmp2_1d(:)
           zsmlt(:)    = zsmlt(:)+ztmp2_1d(:)/(zcons2*zdp_2d(:,jk))
        END IF

        zpretot_1d(:)  = zrfl(:)+zsfl(:)
        zpredel_1d(:)  = zzdrr_1d(:)+zzdrs_1d(:)


        ll1_1d(:)      = (zpretot_1d(:) > zpredel_1d(:))

!--- CCMod -------------------------------------------------------------

        zccvol_all(:,jk) = zccvol_new(:,jk) + zccvol_ges2(:,jk)
        zclcpre(:)     = MERGE(zclcpre(:), (zclcaux(:)+zccvol_all(:,jk)), ll1_1d(:))

!--- End CCMod -------------------------------------------------------------

        zpresum_1d(:)  = zpretot_1d(:)+zpredel_1d(:)
        
        ll1_1d(:) = (zpresum_1d(:) > cqtmin)

!--- CCMod -------------------------------------------------------------

        ztmp1_1d(:) = ((zclcaux(:)+zccvol_all(:,jk))*zpredel_1d(:)  &
                   + zclcpre(:)*zpretot_1d(:)) / MAX(zpresum_1d(:), cqtmin)

!--- End CCMod -------------------------------------------------------------

        ztmp1_1d(:) = MIN(ztmp1_1d(:), 1.0_dp)
        ztmp1_1d(:) = MAX(ztmp1_1d(:), 0.0_dp)

        zclcpre(:) = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))

!---Included for in-cloud scavenging (Philip Stier, 28/03/01):----------
!   Corrected by Junhua Zhang, Philip Stier (01/2004)

        ll1_1d(:) = (zclcpre(:) > zepsec)

        ztmp1_1d(:) = (zrfl(:)+zzdrr_1d(:))/MAX(zclcpre(:), zepsec)
        ztmp2_1d(:) = (zsfl(:)+zzdrs_1d(:))/MAX(zclcpre(:), zepsec)
        ztmp3_1d(:) = (zcons2*zdp_2d(:,jk)*zevp(:))/MAX(zclcpre(:), zepsec)
        ztmp4_1d(:) = (zcons2*zdp_2d(:,jk)*zsub(:))/MAX(zclcpre(:), zepsec)

        zfrain(:,jk)  = MERGE(ztmp1_1d(:), 0._dp, ll1_1d(:))
        zfsnow(:,jk)  = MERGE(ztmp2_1d(:), 0._dp, ll1_1d(:))
        zfevapr(:,jk) = MERGE(ztmp3_1d(:), 0._dp, ll1_1d(:))
        zfsubls(:,jk) = MERGE(ztmp4_1d(:), 0._dp, ll1_1d(:))

!---End Included for scavenging-----------------------------------------
        ! mz_ht_20120403+
        ! rain and snow flux considering incoming rain, melting of snow, 
        ! droplet evaporation / sublimation , but no new production of rain or snow 
        ! in that layer....
        ! (neccessary for impaction scavenging)
        pfrain_no(:,jk)   = zrfl(:) - zcons2*zdp_2d(:,jk)*zevp(:)  
        pfsnow_no(:,jk)   = zsfl(:) - zcons2*zdp_2d(:,jk)*zsub(:)
        ! precipitating cloud cover of this layer is used for the next lower layer 
        ! to estimate the part of the cloud cover in which rain impacts
        pr_cover(:,jk) = zclcpre(:)
        ! mz_ht_20120403-

        zrfl(:)       = zrfl(:)+zzdrr_1d(:)-zcons2*zdp_2d(:,jk)*zevp(:)
        zsfl(:)       = zsfl(:)+zzdrs_1d(:)-zcons2*zdp_2d(:,jk)*zsub(:)

!
!     ------------------------------------------------------------------
!       8.    Updating tendencies of t, q, xl, xi and final cloud cover
!
!
!       8.10   Cloud cover scheme tendencies
!
        IF (lcover .AND. jk >= ncctop) THEN
!
!          Source terms from convection
!          Skewness:
!
!-------------------Added by Junhua Zhang for CONV Micro-----------------------
          IF (ncvmicro > 0) THEN
             zconvskew(:)   = cbeta_cs * (pxtecl(:,jk)+pxteci(:,jk)+pqtec(:,jk)) &
                                       / pbetass(:,jk)
          ELSE
             zconvskew(:)   = cbeta_cs * (pxtec(:,jk)+pqtec(:,jk)) &
                                       / pbetass(:,jk)
          ENDIF
!-------------------------------------end-------------------------------------

           ztmp1_1d(:)  = zdtime_rcp*(cbeta_pq_max-pxskew(:,jk))
           zconvskew(:) = MIN(zconvskew(:), ztmp1_1d(:))
!
!          Convective width now diagnosed, assuming 'a' unchanged:
!
           ll1_1d(:) = (pqm1(:,jk) >= pbetass(:,jk))

           ztmp1_1d(:) = pxskew(:,jk) + zdtime*zconvskew(:) !SF zskewp1
           ztmp2_1d(:) = zwide(:) * (cbeta_pq+ztmp1_1d(:))   &
                                  / (cbeta_pq+pxskew(:,jk))  !SF zbbap1
           ztmp3_1d(:) = zdtime_rcp*(ztmp2_1d(:)-zwide(:))

           zconvvar(:) = MERGE(ztmp3_1d(:), 0._dp, ll1_1d(:))

!
!       8.11 Simple linearized effect of microphysics on skewness
!
           ll1_1d(:) = (pbetaa(:,jk) < pbetass(:,jk)) .AND. &
                       (pbetab(:,jk) > pbetass(:,jk))

           ztmp1_1d(:) = ztmst*(zxlte(:)+zxite(:))                &
                       - zrpr(:)-zsacl(:)-zspr(:)+zcnd(:)+zdep(:) &
                       + zgenti(:)+zgentl(:)                       

           ztmp1_1d(:) = - ztmp1_1d(:) / MAX(zepsec,zbetacl(:))
           ztmp1_1d(:) = MIN(1._dp, ztmp1_1d(:))
           ztmp1_1d(:) = MAX(0._dp, ztmp1_1d(:))
           ztmp1_1d(:) = (pbetass(:,jk)-pbetab(:,jk)) * ztmp1_1d(:)  !SF zmdelb

           ztmp2_1d(:) = (pbetab(:,jk)+ztmp1_1d(:)-pbetaa(:,jk))  &
                       * cbeta_pq/(zbetaqt(:)-pbetaa(:,jk))       &
                       - cbeta_pq
           ztmp2_1d(:) = MAX(MIN(ztmp2_1d(:), cbeta_pq_max), cbeta_pq)

           ztmp3_1d(:) = zdtime_rcp*(ztmp2_1d(:)-pxskew(:,jk))
           ztmp3_1d(:) = MIN(0._dp, ztmp3_1d(:))

           zmicroskew(:) = MERGE(ztmp3_1d(:), zmicroskew(:), ll1_1d(:))

!
!       8.2   New skewness and variance
!
           zxskewte(:)    = zconvskew(:)                 &
                          + zmicroskew(:) + zturbskew(:)

           zxvarte(:)     = zconvvar(:) + zturbvar(:)
!
           ztmp1_1d(:)    = pxvar(:,jk)  + zdtime*zxvarte(:) !SF zvarp1
           ztmp2_1d(:)    = pxskew(:,jk) + zdtime*zxskewte(:) !SF zskew1
!
           pxskew(:,jk)   = MAX(MIN(ztmp2_1d(:), cbeta_pq_max), cbeta_pq)

           ztmp3_1d(:)    = zbetaqt(:)*(1._dp+pxskew(:,jk)/cbeta_pq) !SF zvarmx

           pxvar(:,jk)    = MAX(MIN(ztmp1_1d(:), ztmp3_1d(:)), zvartg(:))
!
        ENDIF !lcover and jk >= ncctop
!
!       8.3   Tendencies of thermodynamic variables
!             Attn: The terms zxisub and zximlt do not appear in
!                   pxite because these processes have already been
!                   included in pxite via changes in cloud ice
!                   sedimentation (see 3.1, 3.2 and 4)
!
!--- CCMod -------------------------------------------------------------

        !zu kleine Kondensstreifen verschwinden
        !Abgleich
        !Eis ohne Wolke oder Wolke ohne Partikel oder kein Eis
        ztmp1_1d(:) = zccicnc_new(:,jk)*zccvol_new(:,jk)+zccicnc(:,jk)*zccvol_ges(:,jk)
        ztmp2_1d(:) = zxibco(:)*zccvol_new(:,jk)+zxibcc(:)*zccvol_ges(:,jk)
        ll1_1d(:) = (ztmp1_1d(:) > 1.e-10_dp .AND. ztmp2_1d(:).gt.1.e-20_dp)
        where (.not. ll1_1d(:))
          where(paclc(:,jk) .gt. 1.e-10_dp .AND. pxim1(:,jk) .gt. 0._dp)
            zcorc2(:) = zcorc2(:) + ztmp2_1d(:)
            zxib(:) = zxib(:) + ztmp2_1d(:)/zclcaux(:)
          elsewhere
            zxievap_cc(:) = zxievap_cc(:) + ztmp2_1d(:)
          end where

          zcccov_ges(:,jk) = 0._dp
          zccvol_ges(:,jk) = 0._dp
          zcciwc(:,jk) = 0._dp
          zxibcc(:) = 0._dp
          zccicnc(:,jk) = 0._dp
          zcclen_ges2(:,jk) = 0._dp
          zcccov_new(:,jk) = 0._dp
          zccvol_new(:,jk) = 0._dp
          zcciwc_new(:,jk) = 0._dp
          zxibco(:) = 0._dp
          zccicnc_new(:,jk) = 0._dp
          zcclen_new(:,jk) = 0._dp
        end where

!--- End CCMod -------------------------------------------------------------

        pqte(:,jk)  = pqte(:,jk)                                     &
                    + ztmst_rcp * ( -zcnd(:)-zgentl(:)+zevp(:)+zxlevap(:)    &
                                  -  zdep(:)-zgenti(:)+zsub(:)+zxievap(:)    &
                                  +  zxisub(:) )

!--- CCMod -------------------------------------------------------------

        pqte(:,jk)  = pqte(:,jk)                                   &
                      + ztmst_rcp * (- zdepcc(:)+ zxisub_cc(:)+ zxievap_cc(:)  &
                                  + (1._dp - B_co(:,jk,jrow))*zfh2o(:,jk) )

!--- End CCMod -------------------------------------------------------------

        ptte(:,jk)  = ptte(:,jk)                                                                     &
                    + ztmst_rcp * ( zlvdcp(:) * (zcnd(:)+zgentl(:)-zevp(:)-zxlevap(:))               &
                                  + zlsdcp(:) * (zdep(:)+zgenti(:)-zsub(:)-zxievap(:)                &
                                                -zxisub(:))+(zlsdcp(:)-zlvdcp(:))                    &
                                              * (-zsmlt(:)-zimlt(:)-zximlt(:)+zfrl(:,jk)+zsacl(:)) )

!--- CCMod -------------------------------------------------------------

        ptte(:,jk)  = ptte(:,jk)                    &
                    + ztmst_rcp * (zlsdcp(:) * (zdepcc(:)-zxisub_cc(:)-zxievap_cc(:)) &
                    +(zlsdcp(:)-zlvdcp(:)) * (-zimlt_cc(:)-zximlt_cc(:)))

!--- End CCMod -------------------------------------------------------------

        ztmp1_1d(:) = pxlte(:,jk) + zxlte(:)

        ztmp2_1d(:) = zimlt(:) + zximlt(:) - zfrl(:,jk) - zrpr(:)       &
                    - zsacl(:) + zcnd(:)   + zgentl(:)  - zxlevap(:)

!--- CCMod -------------------------------------------------------------

        ztmp2_1d(:) = ztmp2_1d(:) + zimlt_cc(:) + zximlt_cc(:)

!--- End CCMod -------------------------------------------------------------

        zxlp1_1d(:) = pxlm1(:,jk) + ztmst * ztmp1_1d(:) + ztmp2_1d(:)

        pxlte(:,jk) = ztmp1_1d(:) + ztmst_rcp * ztmp2_1d(:)

        ztmp1_1d(:) = pxite(:,jk) + zxite(:)
   
        ztmp2_1d(:) = zfrl(:,jk) - zspr(:) + zdep(:) + zgenti(:) - zxievap(:)

!--- CCMod -------------------------------------------------------------

        !Anteil FLugzeugwasser wenn bges=0
        zcorb(:) = MERGE(zcciwc_new(:,jk)*zqcon_1d(:), 0._dp, &
                                     (zcccov_ges(:,jk)+zccvol_new(:,jk)).le.0._dp)
        ztmp2_1d(:) = ztmp2_1d(:) + zcorc2(:) + zcorb(:) + zcorclc(:)

!--- End CCMod -------------------------------------------------------------

        zxip1_1d(:) = pxim1(:,jk) + ztmst * ztmp1_1d(:)  + ztmp2_1d(:) 

        pxite(:,jk) = ztmp1_1d(:) + ztmst_rcp * ztmp2_1d(:)

!
!--- Included for prognostic CDNC/IC scheme ----------------------------

        !--- Calculate new total tendency of CDNC:
!qqq+ this is indeed the total tendency, since zcdnc is initialised 
!     above including the tendency, and here the difference to t-1
!     is calculated ...
        pxtte(:,jk,1) = ztmst_rcp * (zcdnc(:,jk)* m_air/1000._dp &
                                           * zrho_rcp(:,jk) - pxtm1(:,jk,1))
!!$        pxtte(:,jk,1) = pxtte(:,jk,1) + &
!!$             ztmst_rcp * (zcdnc(:,jk)* m_air/1000._dp &
!!$                                           * zrho_rcp(:,jk) - pxtm1(:,jk,1))
!qqq-   
        !--- Update CDNC for radiation:
        pacdnc(:,jk)=zcdnc(:,jk)

        !--- Calculate new total tendency of ICNC:
!qqq+ see above
        pxtte(:,jk,2) = ztmst_rcp * (zicnc(:,jk)* m_air/1000._dp &
                                           * zrho_rcp(:,jk) - pxtm1(:,jk,2))
!!$        pxtte(:,jk,2) = pxtte(:,jk,2) + &
!!$             ztmst_rcp * (zicnc(:,jk)* m_air/1000._dp &
!!$                                           * zrho_rcp(:,jk) - pxtm1(:,jk,2))
!qqq-
        !--- Diagnostics:
        qaut(1:kproma,jk,jrow)      = qaut(1:kproma,jk,jrow) - zdt*zrprn(:)
        qfre(1:kproma,jk,jrow)      = qfre(1:kproma,jk,jrow) - zdt*zfrln(:)
        qacc(1:kproma,jk,jrow)      = qacc(1:kproma,jk,jrow) - zdt*zsacln(:)
        cloud_tm1(1:kproma,jk,jrow) = paclc(1:kproma,jk)

        ll1_1d(:) = (zxlb(:)     >  zeps)    .AND. &
                    (zcdnc(:,jk) >= cdncmin)

        ztmp1_1d(:)         = cdnc_acc(1:kproma,jk,jrow) + zdtime*zcdnc(:,jk)

        cdnc_acc(1:kproma,jk,jrow) = MERGE(zcdnc(:,jk), cdnc_acc(1:kproma,jk,jrow), ll1_1d(:))

        !--- In-cloud CDNC burden:
        ztmp1_1d(:)     = zcdnc_burden(:)+zcdnc(:,jk)*zdz_2d(:,jk)
        zcdnc_burden(:) = MERGE(ztmp1_1d(:), zcdnc_burden(:), ll1_1d(:))
        CDNC_burden_acc(1:kproma,jrow) = zcdnc_burden(:)

        !---- CDNC and burden averaged over cloudy and cloud-free periods
        ztmp1_1d(:)     = cdnc(1:kproma,jk,jrow) + zdtime*zcdnc(:,jk)*zclcaux(:)

        cdnc(1:kproma,jk,jrow) = MERGE(zcdnc(:,jk)*zclcaux(:), cdnc(1:kproma,jk,jrow), ll1_1d(:))

        ztmp1_1d(:) = cdnc_burden(1:kproma,jrow) + zdtime*zcdnc(:,jk)*zdz_2d(:,jk)*zclcaux(:)
        cdnc_burden(1:kproma,jrow) = MERGE(cdnc_burden(1:kproma,jrow) + &
          zcdnc(:,jk)*zdz_2d(:,jk)*zclcaux(:), cdnc_burden(1:kproma,jrow), ll1_1d(:))

        !--- In-cloud effective radius [um]:
        ztmp1_1d(:)  = 0.00045e-6_dp*zcdnc(:,jk) + 1.18_dp !SF zkap
        zreffl_1d(:) = 1.E6_dp*ztmp1_1d(:)                                     &   
                     * ( (3._dp/(4._dp*api*rhoh2o)) * zxlb(:)                  &
                         * zrho(:,jk) / zcdnc(:,jk)           )**(1._dp/3._dp)   !SF zreffl

        zreffl_1d(:) = MAX(4._dp,MIN(zreffl_1d(:),40._dp))
        reffl(1:kproma,jk,jrow) =  MERGE(zreffl_1d(:), reffl(1:kproma,jk,jrow), ll1_1d(:))

        ll2_1d(:) = ll1_1d(:)                  .AND. &
                    (itop(:,jk)  == jk )       .AND. &
                    (ztp1tmp(:)   >  tmelt)    .AND. &
                    (zreffct(:)   <  4._dp)    .AND. &
                    (zreffl_1d(:) >= 4._dp)

        ll1_1d(:) = (zxib(:) > zeps)

        ztmp1_1d(:)         = icnc_acc(1:kproma,jk,jrow) + zdtime*zicnc(:,jk)
        icnc_acc(1:kproma,jk,jrow) = MERGE(zicnc(:,jk), icnc_acc(1:kproma,jk,jrow), ll1_1d(:))

        !--- Cirrus radii [um]:
        ll2_1d(:) = (ztp1tmp(:) > cthomi)

        ztmp1_1d(:) = 0.5e4_dp*( 1000._dp/0.0376_dp*MAX(zxib(:),zeps)*zrho(:,jk)/zicnc(:,jk)         )**0.302_dp !zrieff plate
        ztmp2_1d(:) = 1.0e6_dp*( 3._dp/(4._dp*api*zrhoice)*MAX(zxib(:),zeps)*zrho(:,jk)/zicnc(:,jk)) **0.333_dp

        ztmp2_1d(:) = (1.61_dp*ztmp2_1d(:)**3 + 3.56e-4_dp*ztmp2_1d(:)**6)**0.333_dp   !zrieff 2nd option

        ztmp3_1d(:) = MERGE(ztmp1_1d(:), ztmp2_1d(:), ll2_1d(:))
        ztmp3_1d(:) = MAX(ztmp3_1d(:), ceffmin)
        ztmp3_1d(:) = MIN(ztmp3_1d(:), ceffmax)   !zrieff

        reffi(1:kproma,jk,jrow) = MERGE(ztmp3_1d(:), reffi(1:kproma,jk,jrow), ll1_1d(:))

        !SF tovs diagnostics:
        ll2_1d(:) = ll1_1d(:)       .AND. &
                    .NOT. ll2_1d(:)

        ztmp1_1d(:) = 1000._dp*zxib(:)*zclcaux(:)*zdpg_2d(:,jk)  ! iwp in g/m2
        ztmp2_1d(:) = ztau1i(:) + 1.9787_dp*ztmp1_1d(:)*ztmp3_1d(:)**(-1.0365_dp)
        ztau1i(:)   = MERGE(ztmp2_1d(:), ztau1i(:), ll2_1d(:))
        
        ll3_1d(:) = ll2_1d(:)            .AND. &
                    (ztau1i(:) > 0.7_dp) .AND. &
                    (ztau1i(:) < 3.8_dp)

        !SF end tovs diagnostics 

        ll2_1d(:) = ll1_1d(:)                 .AND. &
                    (zicnc(:,jk) >= zicemin)

        !--- In-cloud ICNC burden:
        ztmp1_1d(:)     = zicnc_burden(:)+zicnc(:,jk)*zdz_2d(:,jk)
        zicnc_burden(:) = MERGE(ztmp1_1d(:), zicnc_burden(:), ll2_1d(:)) 
        ICNC_burden_acc(1:kproma,jrow) = zicnc_burden(:)

        !---- ICNC and burden averaged over cloudy and cloud-free periods
        ztmp1_1d(:)     = icnc(1:kproma,jk,jrow) + zdtime*zicnc(:,jk)*zclcaux(:)
        icnc(1:kproma,jk,jrow) = MERGE(zicnc(:,jk)*zclcaux(:), icnc(1:kproma,jk,jrow), ll2_1d(:))

        ztmp1_1d(:)         = icnc_burden(1:kproma,jrow) + zdtime*zicnc(:,jk)*zdz_2d(:,jk)*zclcaux(:)
        icnc_burden(1:kproma,jrow) =  MERGE(icnc_burden(1:kproma,jrow) + &
          zicnc(:,jk)*zdz_2d(:,jk)*zclcaux(:), icnc_burden(1:kproma,jrow), ll2_1d(:))  

!--- CCMod -------------------------------------------------------------        

        !Initiales Kondensstreifen-Eis
        zcciwc_new(:,jk) = zcciwc_new(:,jk)*zqcon_1d(:)

        where ((zcciwc_new(:,jk)).le.0._dp)
          zcccov_new(:,jk) = 0._dp
          zcclen_new(:,jk) = 0._dp
        end where

        !Contrail-Cirrus-Tendencies

        zccvol_all(:,jk) = zccvol_new(:,jk) + zccvol_ges(:,jk)
        ll_bcc(:) = (zccvol_all(:,jk) .gt. 1.e-5_dp)

        !Eiswassergehalt
        zccxite(:) = zccxite(:) + (zcciwc_new(:,jk)-zdepco_form(:,jk) &
                         +zdepcc(:)-zcorc2(:)-zxievap_cc(:) &
                         -zspr_cc(:)-zcorb(:)-zcorclc(:))*ztmst_rcp

!qqq tendency accumulation?
        pxtte(:,jk,25) = zccxite(:) 
        ztmp1_1d(:) = MAX(pxtm1(:,jk,25) + ztmst*pxtte(:,jk,25), 0._dp)    
        ztmp4_1d(:) = MERGE(ztmp1_1d(:)/zccvol_all(:,jk),0._dp,ll_bcc(:))
        cciwc(1:kproma,jk,jrow) = zrho(:,jk)*ztmp4_1d(:)*zccvol_all(:,jk)

        WHERE(ztmp1_1d(:).le.0._dp)
          zxibcc(:) = 0._dp
          zxibco(:) = 0._dp
          zcccov_ges(:,jk) = 0._dp
          zccvol_ges(:,jk) = 0._dp
          zcccov_new(:,jk) = 0._dp
          zccvol_new(:,jk) = 0._dp
          zccicnc(:,jk) = 0._dp
          zccicnc_new(:,jk) = 0._dp
        ENDWHERE

        !Update
        !Aufteilen der Bedeckungsgradaenderung
        !gescherter Bedeckungsgrad
        where (zcccov_ori(:,jk) .gt. 0._dp)
          zmult3(:) = zcccov(:,jk) / zcccov_ori(:,jk)
        elsewhere
          zmult3(:) = 0._dp
        end where
        zcccov_1dt(:,jk) = zcccov_1dt(:,jk)*zmult3(:)
        zcccov_2dt(:,jk) = zcccov_2dt(:,jk)*zmult3(:)
        zcccov_3dt(:,jk) = zcccov_3dt(:,jk)*zmult3(:)
        zcccov_4dt(:,jk) = zcccov_4dt(:,jk)*zmult3(:)
        zcccov_5dt(:,jk) = zcccov_5dt(:,jk)*zmult3(:)
        zcccov(:,jk)     = zcccov(:,jk)    *zmult3(:)
        !Volumenbedeckungsgrad
        where (zccvol_ges_ori(:,jk) .gt. 0._dp)
          zmult3(:) = zccvol_ges(:,jk) / zccvol_ges_ori(:,jk)
        elsewhere
          zmult3(:) = 0._dp
        end where
        zccvol_1dt(:,jk) = zccvol_1dt(:,jk)*zmult3(:)
        zccvol_2dt(:,jk) = zccvol_2dt(:,jk)*zmult3(:)
        zccvol_3dt(:,jk) = zccvol_3dt(:,jk)*zmult3(:)
        zccvol_4dt(:,jk) = zccvol_4dt(:,jk)*zmult3(:)
        zccvol_5dt(:,jk) = zccvol_5dt(:,jk)*zmult3(:)
        zccvol(:,jk)     = zccvol(:,jk)    *zmult3(:)
        !Aufteilen der Laengenaenderungen
        where (zcclen_ges(:,jk) .gt. 0._dp)
          zmult7(:) = zcclen_ges2(:,jk)/zcclen_ges(:,jk)
        elsewhere
          zmult7(:) = 0._dp
        end where
        zcclen_1dt(:,jk) = zcclen_1dt(:,jk)*zmult7(:)
        zcclen_2dt(:,jk) = zcclen_2dt(:,jk)*zmult7(:)
        zcclen_3dt(:,jk) = zcclen_3dt(:,jk)*zmult7(:)
        zcclen_4dt(:,jk) = zcclen_4dt(:,jk)*zmult7(:)
        zcclen_5dt(:,jk) = zcclen_5dt(:,jk)*zmult7(:)
        zcclen(:,jk)     = zcclen(:,jk)    *zmult7(:)

        ! tracer tendencies
        !coverage
!qqq+
!qqq: all tracer tendencies must be accumulated?, i.e.
!qqq:    pxtte = pxtte + ... ???
!qqq: otherwise previous tendencies (e.g. advection) get lost        
        pxtte(:,jk,5)  = ztmst_rcp * (zcccov_new(:,jk) - pxtm1(:,jk,5))
        pxtte(:,jk,6)  = ztmst_rcp * (zcccov_1dt(:,jk) - pxtm1(:,jk,6))
        pxtte(:,jk,7)  = ztmst_rcp * (zcccov_2dt(:,jk) - pxtm1(:,jk,7))
        pxtte(:,jk,8)  = ztmst_rcp * (zcccov_3dt(:,jk) - pxtm1(:,jk,8))
        pxtte(:,jk,9)  = ztmst_rcp * (zcccov_4dt(:,jk) - pxtm1(:,jk,9))
        pxtte(:,jk,10) = ztmst_rcp * (zcccov_5dt(:,jk) + zcccov(:,jk) - pxtm1(:,jk,10))
        !volume
        pxtte(:,jk,11) = ztmst_rcp * (zccvol_new(:,jk) - pxtm1(:,jk,11))
        pxtte(:,jk,12) = ztmst_rcp * (zccvol_1dt(:,jk) - pxtm1(:,jk,12))
        pxtte(:,jk,13) = ztmst_rcp * (zccvol_2dt(:,jk) - pxtm1(:,jk,13))
        pxtte(:,jk,14) = ztmst_rcp * (zccvol_3dt(:,jk) - pxtm1(:,jk,14))
        pxtte(:,jk,15) = ztmst_rcp * (zccvol_4dt(:,jk) - pxtm1(:,jk,15))
        pxtte(:,jk,16) = ztmst_rcp * (zccvol_5dt(:,jk) + zccvol(:,jk) - pxtm1(:,jk,16))
        !length
        pxtte(:,jk,17) = ztmst_rcp * (zcclen_new(:,jk) - pxtm1(:,jk,17))
        pxtte(:,jk,18) = ztmst_rcp * (zcclen_1dt(:,jk) - pxtm1(:,jk,18))
        pxtte(:,jk,19) = ztmst_rcp * (zcclen_2dt(:,jk) - pxtm1(:,jk,19))
        pxtte(:,jk,20) = ztmst_rcp * (zcclen_3dt(:,jk) - pxtm1(:,jk,20))
        pxtte(:,jk,21) = ztmst_rcp * (zcclen_4dt(:,jk) - pxtm1(:,jk,21))
        pxtte(:,jk,22) = ztmst_rcp * (zcclen_5dt(:,jk) + zcclen(:,jk) - pxtm1(:,jk,22))
        pxtte(:,jk,23) = ztmst_rcp * (zcclen_ou(:,jk) - pxtm1(:,jk,23))
!qqq-

        zcccov_all(:,jk) = zcccov_ges(:,jk) + zcccov_new(:,jk)
        zccvol_all(:,jk) = zccvol_ges(:,jk) + zccvol_new(:,jk)

        cccov(1:kproma,jk,jrow) = zcccov_all(:,jk)
        ccvol(1:kproma,jk,jrow) = zccvol_all(:,jk)
        cctau(1:kproma,jk,jrow) = zcccov_new(:,jk)

        ll_bcc(:) = (zccvol_all(:,jk) .gt. 1.e-5_dp)


        !ice crystal number concentration
        ztmp1_1d(:) = zccicnc(:,jk)     * (zccvol_ges(:,jk)-zconcovinitial(:,jk)*zccvol_ges(:,jk))
        ztmp2_1d(:) = zccicnc_new(:,jk) * (zccvol_new(:,jk)+zconcovinitial(:,jk)*zccvol_ges(:,jk))
        ztmp1_1d(:) = ztmp1_1d(:)+ztmp2_1d(:)
!qqq+
        pxtte(:,jk,24) = ztmst_rcp * (ztmp1_1d(:)* m_air/1000._dp &
                                           * zrho_rcp(:,jk) - pxtm1(:,jk,24))
!!$        pxtte(:,jk,24) = pxtte(:,jk,24) + &
!!$             ztmst_rcp * (ztmp1_1d(:)* m_air/1000._dp &
!!$                                           * zrho_rcp(:,jk) - pxtm1(:,jk,24))
!qqq+
        ztmp1_1d(:) = MERGE(ztmp1_1d(:)/zccvol_all(:,jk),0._dp,ll_bcc(:))
        ccicnc(1:kproma,jk,jrow) = ztmp1_1d(:)

        !Volumenradius
        where(ztmp1_1d(:).gt.zeps)
          ztmp3_1d(:) = 1.e6_dp*(3._dp*ztmp4_1d(:)*zrho(:,jk)/(4._dp*api*zrhoice*ztmp1_1d(:)))**0.333
          ztmp3_1d(:) = MAX(ztmp3_1d(:), 0.01_dp)
          ztmp3_1d(:) = MIN(ztmp3_1d(:), 1.e6_dp)
        elsewhere
          ztmp3_1d(:) = 0._dp
        endwhere

        !--- Effective Cirrus radii [um]:
        ll2_1d(:) = (ztp1tmp(:) > cthomi)
        where(ztmp1_1d(:).gt.zeps)
          ztmp2_1d(:) = 0.5e4_dp*( 1000._dp/0.0376_dp*ztmp4_1d(:)*zrho(:,jk)/ztmp1_1d(:))**0.302_dp !zrieff plate
          ztmp5_1d(:) = 1.0e6_dp*( 3._dp/(4._dp*api*zrhoice) *ztmp4_1d(:)*zrho(:,jk)/ztmp1_1d(:))**0.333_dp
          ztmp5_1d(:) = (1.61_dp*ztmp5_1d(:)**3 + 3.56e-4_dp*ztmp5_1d(:)**6)**0.333_dp   !zrieff 2nd option
          ztmp3_1d(:) = MERGE(ztmp2_1d(:), ztmp5_1d(:), ll2_1d(:))
          ztmp3_1d(:) = MAX(ztmp3_1d(:), ceffmin)
          ztmp3_1d(:) = MIN(ztmp3_1d(:), ceffmax)   !zrieff
        elsewhere
          ztmp3_1d(:) = 0._dp
        end where
 
        ccreffi(1:kproma,jk,jrow)= ztmp3_1d(:)

        !optische Dicke
        ztmp1_1d(:) = 1000._dp*ztmp4_1d(:)*zccvol_all(:,jk)/MAX(zeps,zcccov_all(:,jk))*zdpg_2d(:,jk)  !iwp in g/m2
        ztmp2_1d(:) = 1.9787_dp*ztmp1_1d(:)*ztmp3_1d(:)**(-1.0365_dp)
        zcctau(:,jk) = ztmp2_1d(:)
 
!--- End CCMod ---------------------------------------------------------     

!--- End included for CDNC/IC scheme -----------------------------------

!       8.4   Corrections: Avoid negative cloud water/ice
!
        ztmp1_1d(:)   = zxlp1_1d(:)      !SF zxlold

        ll1_1d(:)      = (zxlp1_1d(:) < ccwmin)

        zdxlcor_1d(:) = -ztmst_rcp * zxlp1_1d(:)
        zdxlcor_1d(:) = MERGE(zdxlcor_1d(:), 0._dp, ll1_1d(:))
        pxlte(:,jk)   = pxlte(:,jk) + zdxlcor_1d(:)

        ztmp1_1d(:)          = pxtte(:,jk,1) - ztmst_rcp*m_air*zcdnc(:,jk)*zrho_rcp(:,jk)*1.e-3_dp
        pxtte(:,jk,1) = MERGE(ztmp1_1d(:), pxtte(:,jk,1), ll1_1d(:))

        ll2_1d(:)     = (zxip1_1d(:) < ccwmin)

        zdxicor_1d(:) = -ztmst_rcp * zxip1_1d(:) 
        zdxicor_1d(:) = MERGE(zdxicor_1d(:), 0._dp, ll2_1d(:)) 
        pxite(:,jk)   = pxite(:,jk) + zdxicor_1d(:)

        ztmp1_1d(:)          = pxtte(:,jk,2) - ztmst_rcp*m_air*zicnc(:,jk)*zrho_rcp(:,jk)*1.e-3_dp
        pxtte(:,jk,2) = MERGE(ztmp1_1d(:), pxtte(:,jk,2), ll2_1d(:))

        paclc(:,jk)   = MERGE(0.0_dp,paclc(:,jk),ll1_1d(:) .AND. ll2_1d(:))
        paclcac(:,jk) = paclcac(:,jk) + zdtime*paclc(:,jk)

        pqte(:,jk)    = pqte(:,jk) - zdxlcor_1d(:) - zdxicor_1d(:)
        ptte(:,jk)    = ptte(:,jk) + zlvdcp(:)*zdxlcor_1d(:)                &
                                   + zlsdcp(:)*zdxicor_1d(:)

!SF for scavenging: create a 2d array of clcpre:
        zclcpre_2d(:,jk) = zclcpre(:)

831 END DO    ! Vertical loop
!
!     ------------------------------------------------------------------
!
!       10.    Diagnostics
!
!       10.1   Accumulated precipitation at the surface
!
   prsfl(:) = zrfl(:)
   pssfl(:) = zsfl(:)
   paprl(:) = paprl(:) + zdtime*(prsfl(:)+pssfl(:))
   paprs(:) = paprs(:) + zdtime*pssfl(:)

!
!       10.2   Total cloud cover
!
    zclcov(:) = 1.0_dp-paclc(:,1)

    DO 923 jk = 2,klev
       ztmp1_1d(:) = MAX(paclc(:,jk), paclc(:,jk-1))
       ztmp2_1d(:) = MIN(paclc(:,jk-1), zxsec)

       zclcov(:) = zclcov(:)*(1._dp - ztmp1_1d(:))  &
                            /(1._dp - ztmp2_1d(:))
923 END DO

    zclcov(:)  = 1.0_dp-zclcov(:)
    paclcov(:) = paclcov(:) + zdtime*zclcov(:)

!--- CCMod ---------------------------------------------------------     

    !Total contrail cirrus cloud coverage spreaded
    zclcov(:) = 1.0_dp-zcccov_all(:,1)
    DO jk = 2,klev
      ztmp1_1d(:) = MAX(zcccov_all(:,jk), zcccov_all(:,jk-1))
      ztmp2_1d(:) = MIN(zcccov_all(:,jk-1), zxsec)

      zclcov(:) = zclcov(:)*(1._dp - ztmp1_1d(:))  &
                             /(1._dp - ztmp2_1d(:))
    END DO
    zclcov(:)  = 1.0_dp-zclcov(:)
    cccov_tot(1:kproma,jrow) = zclcov(:)

    !Total contrail cirrus cloud coverage spreaded tau>0.05
    zcccov_tau(:,:) = MERGE(zcccov_all(:,:), 0._dp, zcctau(:,:)>0.05)
    zclcov(:) = 1.0_dp-zcccov_tau(:,1)
    DO jk = 2,klev
      ztmp1_1d(:) = MAX(zcccov_tau(:,jk), zcccov_tau(:,jk-1))
      ztmp2_1d(:) = MIN(zcccov_tau(:,jk-1), zxsec)

      zclcov(:) = zclcov(:)*(1._dp - ztmp1_1d(:))  &
                             /(1._dp - ztmp2_1d(:))
    END DO
    zclcov(:)  = 1.0_dp-zclcov(:)
    cccov_tot2(1:kproma,jrow) = zclcov(:)

!--- End CCMod ---------------------------------------------------------     
!
!       10.3   Vertical integrals of humidity, cloud water and cloud ice
!
    zqvi(:)  = 0.0_dp
    zxlvi(:) = 0.0_dp
    zxivi(:) = 0.0_dp
!
    DO 933 jk = ktdia,klev
       zqvi(:)  = zqvi(:)  + pqm1(:,jk) *zdpg_2d(:,jk)
       zxlvi(:) = zxlvi(:) + pxlm1(:,jk)*zdpg_2d(:,jk) 
       zxivi(:) = zxivi(:) + pxim1(:,jk)*zdpg_2d(:,jk)
933 END DO
!
    pqvi(:)  = pqvi(:)  + zdtime*zqvi(:)
    pxlvi(:) = pxlvi(:) + zdtime*zxlvi(:)
    pxivi(:) = pxivi(:) + zdtime*zxivi(:)

!qqq: here tendecies are again overwritten, and a m1-value as well ???
    !tracer for cloud fraction
    pxtm1(:,:,4) = pxtm1(:,:,3)
    pxtte(:,:,4) = pxtte(:,:,3)
    pxtm1(:,:,3) = paclc(:,:)
    pxtte(:,:,3) = 0.0_dp
!
  RETURN
END SUBROUTINE cloud_cdnc_icnc_cc


!===============================================================================

! SUBROUTINES originally from mo_cloud_utils

 SUBROUTINE get_util_var(kproma, kbdim, ktdia, klev, klevp1, &
                         paphm1, pgeo, papm1, ptm1,          &
                         pgeoh, pdp, pdpg, pdz, paaa, pviscos)

 !Get several utility variables:
 !    geopotential at half levels (pgeoh)
 !    pressure- and height-differences (pdp, pdz)
 !    air density correction for computing ice crystal fall velocity (paaa)
 !    dynamic viscosity of air (pviscos)
 !
   INTEGER, INTENT(IN) :: kproma, kbdim, ktdia, klev, klevp1

   REAL(dp), INTENT(IN)  :: paphm1(kbdim,klevp1), pgeo(kbdim,klev), papm1(kbdim,klev), ptm1(kbdim,klev)
   REAL(dp), INTENT(OUT) :: pgeoh(kbdim,klevp1),  pdp(kbdim,klev), pdpg(kbdim,klev), pdz(kbdim,klev), &
                            paaa(kbdim,klev), pviscos(kbdim,klev)

   REAL(dp) :: g_rcp

   g_rcp = 1._dp / g

   !Geopotential at half levels:
   pgeoh(:,ktdia+1:klev) = 0.5_dp * (pgeo(:,ktdia+1:klev) + pgeo(:,ktdia:klev-1))
   pgeoh(:,ktdia)        = pgeo(:,ktdia) + (pgeo(:,ktdia)-pgeoh(:,ktdia+1))
   pgeoh(:,klevp1)       = 0.0_dp

   !Pressure differences:
   pdp(:,ktdia:klev)           = paphm1(:,ktdia+1:klevp1) - paphm1(:,ktdia:klev)
   pdpg(:,ktdia:klev)          = g_rcp * pdp(:,ktdia:klev)

   !Height differences:
   pdz(:,ktdia:klev)           = g_rcp * (pgeoh(:,ktdia:klev) - pgeoh(:,ktdia+1:klevp1))

   !Air density correction:
   paaa(:,:) = ( (papm1(:,:)/30000._dp)**(-0.178_dp) ) * ( (ptm1(:,:)/233.0_dp)**(-0.394_dp) )

   !Dynamic viscosity of air:
   pviscos(:,:) = (1.512_dp + 0.0052_dp*(ptm1(:,:)-233.15_dp))*1.e-5_dp

 END SUBROUTINE get_util_var

!-----------------------------------------------------------------------------------

 SUBROUTINE get_cloud_bounds(kproma, kbdim, ktdia, klev, paclc, ktop, kbas, kcl_minustop, kcl_minusbas)

 ! *get_cloud_bounds* flags the top, base and intermediate levels for each cloud  
 !
   INTEGER, INTENT(IN)  :: kproma, kbdim, ktdia, klev
   INTEGER, INTENT(OUT) :: ktop(kbdim,klev), &      !flag for cloud tops 
                           kbas(kbdim,klev), &      !flag for cloud bases
                           kcl_minustop(kbdim,klev), & !flag for all cloud levels excepted their top
                          kcl_minusbas(kbdim,klev) !flag for all cloud levels excepted their base

   REAL(dp), INTENT(IN) :: paclc(kbdim,klev)    !cloud cover

   !Local variables:

   INTEGER  :: jk, jl, jm, jnumb, jtop, jbas
   INTEGER  :: iindex(kbdim,klev), &  !index of levels
               iclnb(kbdim),       &  !number of clouds per column
               iclbounds(2,klev/2+1)  !bounds infos per cloud

   REAL(dp) :: zaclcm(kbdim,klev), zaclcp(kbdim,klev)

   LOGICAL  :: ll(kbdim,klev), llm(kbdim,klev), llp(kbdim,klev), lltop(kbdim,klev), llbas(kbdim,klev)

   !Initialization:
   zaclcm(:,:)       = 0._dp
   zaclcp(:,:)       = 0._dp
   iclnb(:)          = 0
   kcl_minustop(:,:) = 0
   kcl_minusbas(:,:) = 0

   ll(:,:)    = .false.
   llm(:,:)   = .false.
   llp(:,:)   = .false.
   lltop(:,:) = .false.
   llbas(:,:) = .false.
 
   DO jk=ktdia,klev
      iindex(:,jk)=jk
   ENDDO 

   !Duplicates paclc at level-1 and level+1
   zaclcm(:,ktdia+1:klev) = paclc(:,ktdia:klev-1)
   zaclcp(:,ktdia:klev-1) = paclc(:,ktdia+1:klev)

   !Sets logical switches
   ll(:,:)  = (paclc(:,:)  >= zepsec)   
   llm(:,:) = (zaclcm(:,:) <  zepsec)  
   llp(:,:) = (zaclcp(:,:) <  zepsec)

   lltop(:,:) = (ll(:,:) .AND. llm(:,:)) !true if cloud top detected
   llbas(:,:) = (ll(:,:) .AND. llp(:,:)) !true if cloud base detected

   !Sets itop and ibas
   ktop(:,:) = MERGE(iindex(:,:),0,lltop(:,:))
   kbas(:,:) = MERGE(iindex(:,:),0,llbas(:,:))

   !Resets the logical switches
   lltop(:,:) = .false.
   llbas(:,:) = .false.

   lltop(:,:) = (ktop(:,:) > 0)
   llbas(:,:) = (kbas(:,:) > 0)
 
   !Counts the number of clouds per column
   iclnb(:) = COUNT(lltop(:,:),2)

   DO jl=1,kproma
      jnumb = iclnb(jl)

      !Sets the bounds in a compact array
      iclbounds(:,:) = 0

      iclbounds(1,1:jnumb) = PACK(ktop(jl,:),lltop(jl,:))   !cloud tops
      iclbounds(2,1:jnumb) = PACK(kbas(jl,:),llbas(jl,:))   !cloud bases

      !Flags cloud levels excepted their base (or top -this later is not yet needed-)
      DO jm=1,jnumb
         jtop=iclbounds(1,jm)
         jbas=iclbounds(2,jm)
         kcl_minusbas(jl,jtop:jbas-1)=jbas
         kcl_minustop(jl,jtop+1:jbas)=jtop
      ENDDO
   ENDDO

 END SUBROUTINE get_cloud_bounds


!===============================================================================

! SUBROUTINES originally from mo_cirrus

! ------------------------------------------------------------------------------

  REAL(dp) FUNCTION SCRHOM(T)

     ! ***** CRITICAL ICE SATURATION RATIO FOR HOMOGENEOUS FREEZING
     ! ***** EVALUATED FOR AN INTERMEDIATE PARTICLE SIZE OF 0.25 MUM
     ! ***** TH.KOOP, B.P.LUO, A.TSIAS, TH.PETER, NATURE 406, 611-614, 2000

     IMPLICIT NONE

     ! function parameters

     REAL(dp), INTENT(in) :: T

     ! local variables

     REAL(dp) :: TEMP

     TEMP = MIN(240.0_dp,MAX(170.0_dp,T))
     SCRHOM = 2.418_dp - (TEMP/245.68_dp)

     RETURN

  END FUNCTION SCRHOM

! ------------------------------------------------------------------------------

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

! ------------------------------------------------------------------------------

  REAL(dp) FUNCTION XEERFC(Y)

     ! ***** PRODUCT EXP(1/Y) * ERFC(1/SQRT(Y)) AND ASYMPTOTES

     IMPLICIT NONE

     ! function parameters

     REAL(dp), INTENT(in) :: Y

     ! parameters

     REAL(dp), PARAMETER :: SQPI = 1.7724539_dp

     ! local variables

     REAL(dp) :: Y1,SQY1,PPROD,POLY
     INTEGER :: K

!     ! functions
!
!     REAL(dp) :: XERF ! error function with rational approximation

     Y1   = 1.0_dp / Y
     SQY1 = SQRT(Y1)
     IF (Y.LE.0.2_dp) THEN
        PPROD = 1.0_dp
        POLY  = 1.0_dp
        DO K = 1, 4
           PPROD = PPROD * REAL(2*K-1,dp)
           POLY  = POLY  + PPROD * (-0.5_dp*Y)**K
        END DO
        XEERFC = POLY / SQY1 / SQPI
     ELSEIF (Y.GE.2.0_dp) THEN
        PPROD  = 1.0_dp
        POLY   = 1.0_dp
        DO K = 1, 4
           PPROD = PPROD * REAL(2*K+1,dp)
           POLY  = POLY  + REAL(2**K,dp) * Y1**K / PPROD
        END DO
        XEERFC = EXP(Y1) - 2.0_dp*SQY1/SQPI * POLY
     ELSE
        XEERFC = ( 1.0_dp - XERF(SQY1) ) * EXP(Y1)
     END IF

     RETURN

  END FUNCTION XEERFC

! ------------------------------------------------------------------------------

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

     REAL(dp) :: T,ARG

     T    = 1.0_dp / ( 1.0_dp + P * ABS(X) )
     ARG  = MIN( X**2, 75.0_dp )
     XERF = 1.0_dp - ( A1*T + A2*T**2 + A3*T**3 + A4*T**4 + A5*T**5 ) &
            * EXP( -ARG )
     XERF = SIGN( XERF, X )

     RETURN

     END FUNCTION XERF

! ------------------------------------------------------------------------------

 REAL(dp) FUNCTION SCRHET(T, NFRZ)

     IMPLICIT NONE

     ! function parameters

     REAL(dp), INTENT(in) :: T
     INTEGER,  INTENT(in) :: NFRZ

     ! ***** CRITICAL ICE SATURATION RATIO FOR HETEROGENEOUS FREEZING

     IF (NFRZ.EQ.0) THEN

        ! ***** CONSTANT THRESHOLD (FREEZING OF 1 PTCL PER SEC, RADIUS 0.25 MUM)
        SCRHET = 1.3_dp

     ELSE
        STOP 'VALUE IFRZ NOT IMPLEMENTED.'
     END IF

     RETURN

  END FUNCTION SCRHET

! ------------------------------------------------------------------------------

  REAL(dp) FUNCTION TAUG(B,Y,X0)

     ! ***** DIMENSIONLESS GROWTH TIME SCALE
     ! ***** B.K??RCHER AND S.SOLOMON, JGR 104(D22), 27441-27459, 1999

     IMPLICIT NONE

     REAL(dp), INTENT(in) :: B, Y, X0

     REAL(dp), PARAMETER :: SQ31 = 0.577350269_dp
     REAL(dp), PARAMETER :: SIX1 = 0.166666666_dp

     REAL(dp) :: X, F1, F10, F2, F20

     TAUG = 0.0_dp

     IF (Y.LE.X0) RETURN

     X    = MIN( Y, 0.999_dp )
     F1   = SIX1 * LOG( (1.0_dp+X +X**2)  / (1.0_dp-X)**2 )
     F10  = SIX1 * LOG( (1.0_dp+X0+X0**2) / (1.0_dp-X0)**2 )
     F2   = SQ31 * ATAN( SQ31*(1.0_dp+2.0_dp*X) )
     F20  = SQ31 * ATAN( SQ31*(1.0_dp+2.0_dp*X0) )
     TAUG = (B+1.0_dp)*(F1-F10) + (B-1.0_dp)*(F2-F20)

     RETURN

  END FUNCTION TAUG

! ------------------------------------------------------------------------------

  SUBROUTINE XFRZMSTR(LHET,NOSIZE,ZTMST,                   &
                      KLEV,kbdim,kproma,JK, nfrzmod,       &
                      ZSUSATIX,ZVERVX,ZAPNX,               &
                      ZAPRX,ZAPSIGX,PTM1,TMELT,ZMIN,PAPM1, &
                      ZTHOMI,ZRI,ZNICEX)

    ! This is the calling loop for XFRZHOM/XFRZHET extracted from CLOUD

    IMPLICIT NONE

    ! SUBROUTINE parameters

    LOGICAL,  INTENT(in)  :: NOSIZE  ! .true. ---> aerosol size effects are
                                     !             ignored for homogeneous
                                     !             freezing
    LOGICAL,  INTENT(in)  :: LHET    ! .true. ---> heterogeneous freezing of
                                     !             aerosol particles below 235K
                                     !             is considered
                                     !             aerosol particles below 235K
                                     !             is considered

    INTEGER,  INTENT(in)  :: KLEV    ! number of levels
    INTEGER,  INTENT(in)  :: KBDIM   ! size of arrays
    INTEGER,  INTENT(in)  :: KPROMA  ! number of cells in arrays
    INTEGER,  INTENT(in)  :: JK      ! current level
    INTEGER,  INTENT(in)  :: NFRZMOD  ! number of aerosol modes
    REAL(dp), INTENT(in)  :: ZSUSATIX(kbdim)  ! ice supersaturation pressure
    REAL(dp), INTENT(in)  :: ZVERVX(kbdim)    ! updraft [cm/s]
    REAL(dp), INTENT(in)  :: ZAPNX(kbdim,NFRZMOD)  ! aerosol number conc. [1/cm3]
    REAL(dp), INTENT(in)  :: ZAPRX(kbdim,NFRZMOD)  ! aerosol radius [cm]
    REAL(dp), INTENT(in)  :: ZAPSIGX(kbdim,NFRZMOD)! aerosol standard deviation
                                                  ! (log-normal distribution)
    REAL(dp), INTENT(in)  :: PTM1(kbdim,KLEV) ! temperature [K]
    REAL(dp), INTENT(in)  :: PAPM1(kbdim,KLEV)! pressure [Pa]
    REAL(dp), INTENT(in)  :: ZTMST            ! time step [s]
    REAL(dp), INTENT(in)  :: TMELT            ! melting temperature
    REAL(dp), INTENT(in)  :: ZMIN             ! "epsilon"
    REAL(dp), INTENT(in)  :: ZTHOMI           ! temperature of hom. freezing

    REAL(dp) ZNICEX(kbdim,klev)    ! number of ice nuclei
    REAL(dp) ZRI(kbdim,KLEV)  ! 

    ! local variables

    REAL(dp) :: ZAPN(NFRZMOD)    ! aerosol number concentration [1/cm3]
    REAL(dp) :: ZAPR(NFRZMOD)    ! aerosol radius [cm]
    REAL(dp) :: ZAPSIG(NFRZMOD)  ! aerosol standard deviation (log-normal distr.)

    INTEGER  :: IXLIST(kbdim)
    INTEGER  :: IX2LIST(kbdim)
    INTEGER  :: KLIST(kbdim)
    REAL(dp) :: PWX(kbdim), COOLRX(kbdim)

    INTEGER  :: IXC,JL,IX,I,IX0,IXC2,K
    REAL(dp) :: ZSUSATI,ZVERV,ZSATI,COOLR,T,DT,SI,TMIN,SCR
    REAL(dp) :: PW,PWCR,SL,PISATTMIN,ZNICE,ZRI1,PICE,ZPHPA,TEMP,CTOT

    ! parameters

    INTEGER,  PARAMETER :: KMAX   = 120

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

!    ! functions
!
!    REAL(dp) :: PISAT  ! vapor pressure over ice [hPa]
!    REAL(dp) :: SCRHOM ! critical ice saturation ratio for homogeneous freezinh

    IXC=0

    DO 361 JL = 1, KPROMA
       ZSUSATI = ZSUSATIX(JL)
       ZVERV   = ZVERVX(JL)
       ZAPN(:) = ZAPNX(JL,:)
       zri(jl,jk)=0._dp
       znicex(jl,jk)=0._dp

       IF (ZSUSATI.GT.0.0_dp.AND.PTM1(JL,JK).LT.ZTHOMI.AND.ZVERV.GT.ZMIN) THEN
          ZSATI=ZSUSATI+1.0_dp
          ZNICE=0.0_dp
          ZRI1=0.0_dp
          IF (PTM1(JL,JK).LT.ZTHOMI) THEN
             CTOT=0.0_dp
             DO I=1,NFRZMOD
                CTOT=CTOT+ZAPN(I)
             END DO
             IF ( CTOT > 0.0_dp .AND. ZVERV > 0.0_dp ) THEN

                !*****CHECK IF FREEZING CONDITIONS ARE MET WITHIN DT

                COOLR  = GRAV * ZVERV / CPAIR
                T = PTM1(JL,JK)
                DT = ZTMST
                SI = ZSATI
                TMIN   = T    - COOLR * DT
                TMIN   = MAX( TMIN, 170.0_dp )
                SCR    = SCRHOM(TMIN)
                PW     = SI * PISAT(T)
                PWCR   = PW * ( TMIN / T )**3.5_dp
                PISATTMIN=PISAT(TMIN)
                IF ((PWCR/PISATTMIN).GE. SCR) THEN
                   COOLRX(JL)=COOLR
                   PWX(JL)=PW
                   IXC=IXC+1
                   IXLIST(IXC)=JL
                END IF
             END IF
          END IF
       END IF
 361 END DO ! JL-loop

    ! *****FREEZING TEMPERATURE logic

    IX0=IXC
    IXC=0
    IXC2=0
    DO IX=1,IX0
       JL=IXLIST(IX)
       ZSUSATI=ZSUSATIX(JL)
       ZSATI=ZSUSATI+1.0_dp
       T = PTM1(JL,JK)
       SI=ZSATI
       SCR    = SCRHOM(T)
       IF (SI.GE.SCR) THEN
          ! Need not to figure out a K
          IXC=IXC+1
          IXLIST(IXC)=JL
       ELSE
          ! Need to figure out a K
          IXC2=IXC2+1
          IX2LIST(IXC2)=JL
       END IF
    END DO ! IX-loop

    DO K=1,KMAX+1
       IX0=IXC2
       IXC2=0
       IF ( IX0 > 0 ) THEN
          DO IX=1,IX0
             JL=IX2LIST(IX)
             T = PTM1(JL,JK)
             PW = PWX(JL)
             SL    = ( T - 170.0_dp ) / REAL(KMAX,dp)
             TEMP = T  - SL * REAL(K-1,dp)
             SCR  = SCRHOM(TEMP)
             PICE = PISAT(TEMP)
             PWCR = PW * (TEMP/T)**3.5_dp
             IF ((PWCR/PICE).GE.SCR) THEN
                ! The K condition has been reached.
                ! Keep the K in mind and add this JL point to the list
                ! of points to be passed to XFRZHOM
                KLIST(JL)=K
                IXC=IXC+1
                IXLIST(IXC)=JL
             ELSE
                ! Keep for further K iteration
                IXC2=IXC2+1
                IX2LIST(IXC2)=JL
             END IF
          END DO ! IX-loop
       END IF
    END DO ! K-loop

    DO IX=1,IXC
       JL=IXLIST(IX)
       ZSUSATI=ZSUSATIX(JL)
       ZVERV=ZVERVX(JL)
       ZAPN(:)=ZAPNX(JL,:)
       ZAPR(:)=ZAPRX(JL,:)
       ZAPSIG(:)=ZAPSIGX(JL,:)
       COOLR=COOLRX(JL)
       PW=PWX(JL)

       ZPHPA = PAPM1(JL,JK)*0.01_dp
       ZSATI = ZSUSATI+1.0_dp
       ZNICE = 0.0_dp
       ZRI1  = 0.0_dp

       T = PTM1(JL,JK)
       SI = ZSATI

       ! ***** FREEZING TEMPERATURE

       SCR    = SCRHOM(T)
       IF (SI.GE.SCR) THEN
          TEMP  = T
          PICE  = PISAT(TEMP)
          PW    = SCR * PICE
          PWCR  = PW
       ELSE
          SL   = ( T - 170.0_dp ) / REAL(KMAX,dp)
          K    = KLIST(JL)
          TEMP = T  - SL * REAL(K-1,dp)
          SCR  = SCRHOM(TEMP)
          PICE = PISAT(TEMP)
          PWCR = PW * (TEMP/T)**3.5_dp
       END IF

       IF (LHET) THEN
          CALL XFRZHET (NOSIZE,nfrzmod,ZTMST,ZAPN,ZAPR,ZAPSIG,ZPHPA,    &
                        PTM1(JL,JK),ZVERV,ZSATI,ZNICE,ZRI1, COOLR, PW, &
                        SCR, TEMP, PICE, PWCR)
       ELSE
          CALL XFRZHOM (NOSIZE,nfrzmod,ZTMST,ZAPN,ZAPR,ZAPSIG,ZPHPA,    &
                        PTM1(JL,JK),ZVERV,ZSATI,ZNICE,ZRI1, COOLR, PW, &
                        SCR, TEMP, PICE, PWCR)
       END IF

       ZRI(JL,JK) = ZRI1 *1.0E-2_dp
       ZNICEX(JL,JK) = ZNICE*1.0E6_dp

    END DO ! IX-loop

    RETURN

  END SUBROUTINE XFRZMSTR

! ------------------------------------------------------------------------------

  SUBROUTINE XFRZHET (NOSIZE,nfrzmod,DT,C,R,SIG,P,T,V,SI,CI,RI,COOLR, PW, &
                      SCR, TEMP, PICE, PWCR)

     ! ***** HETEROGENEOUS FREEZING (ADIABATIC ASCENT)
     !
     ! ***** BERND K??RCHER  APR 14  2003
     ! ***** bernd.kaercher@dlr.de  http://www.op.dlr.de/~pa1c/
     !
     ! ***** Ulrike Lohmann      Dalhousie Univ   Apr 14 03
     !       Ulrike.Lohmann@Dal.Ca
     ! ***** Johannes Hendricks  DLR-IPA          Apr 14 03
     !       Johannes.Hendricks@dlr.de
     !
     ! ***** References
     !       K??rcher, B. and U. Lohmann
     !       A parameterization of cirrus cloud formation: Heterogeneous
     !       freezing.
     !       J. Geophys. Res. 108, doi:10.1029/2002JD003220, in press, 2003.

     IMPLICIT NONE

     ! SUBROUTINE parameters

     LOGICAL,  INTENT(in)  :: NOSIZE     ! .true. ---> aerosol size effects are
                                         !             ignored for homogeneous
                                         !             freezing
     INTEGER,  INTENT(in)  :: NFRZMOD     ! number of aerosol modes
     REAL(dp), INTENT(in)  :: DT         ! time step [s]
     REAL(dp), INTENT(in)  :: C(NFRZMOD)  ! aerosol number concentration [1/cm3]
     REAL(dp), INTENT(in)  :: R(NFRZMOD)  ! aerosol radius [cm]
     REAL(dp), INTENT(in)  :: SIG(NFRZMOD)! aerosol standard deviation
                                         ! (log-normal distribution)
     REAL(dp), INTENT(in)  :: P          ! pressure [hPa]
     REAL(dp), INTENT(in)  :: T          ! temperature [K]
     REAL(dp), INTENT(in)  :: V          ! updraft [cm/s]

     REAL(dp), INTENT(inout) :: SI       ! saturation pressure (ice)

     REAL(dp), INTENT(out) :: COOLR      ! cooling rate
     REAL(dp), INTENT(out) :: RI         ! 
     REAL(dp), INTENT(out) :: PICE       ! 
     REAL(dp), INTENT(out) :: TEMP       ! 
     REAL(dp), INTENT(out) :: CI         ! ice nuclei
     REAL(dp), INTENT(out) :: SCR        ! 
     REAL(dp), INTENT(out) :: PW         ! 
     REAL(dp), INTENT(out) :: PWCR       ! 

     ! local variables

     INTEGER, SAVE :: NFRZ = 0

     INTEGER :: IFRZ,N,K

     REAL(dp) :: A(3)
     REAL(dp) :: B(2)
     REAL(dp) :: CCR(nfrzmod)

     REAL(dp) :: TMIN,SISVE,SL,VOLF,PCR,SUPSCR,DLNJDT,TAU,DIFFC
     REAL(dp) :: BKT,VTH,CISAT,THETA,XMI0,XMIMAX,RIMAX,XMISAT,CTAU
     REAL(dp) :: TGROW,ZF,XMFP,BETA,RIHAT,X0,Z,XMI,YK,CVF,CTOT,DSCRMIN
     REAL(dp) :: RS,X

     ! parameters

     INTEGER,  PARAMETER :: KMAX   = 120

     REAL(dp), PARAMETER :: RHOICE = 0.925_dp
     REAL(dp), PARAMETER :: PI     = 3.1415927_dp
     REAL(dp), PARAMETER :: SQPI   = 1.7724539_dp
     REAL(dp), PARAMETER :: THIRD  = 0.3333333_dp
     REAL(dp), PARAMETER :: RGAS   = 8.3145E7_dp
     REAL(dp), PARAMETER :: ALPHA  = 0.5_dp
     REAL(dp), PARAMETER :: BK     = 1.3807E-16_dp
     REAL(dp), PARAMETER :: CPAIR  = 1.00467E7_dp
     REAL(dp), PARAMETER :: HEAT   = 2830.3E7_dp
     REAL(dp), PARAMETER :: AVOG   = 6.02213E23_dp
     REAL(dp), PARAMETER :: GRAV   = 981.0_dp
     REAL(dp), PARAMETER :: SVOL   = 3.23E-23_dp
     REAL(dp), PARAMETER :: XMW    = 2.992E-23_dp
     REAL(dp), PARAMETER :: XWW    = 18.016_dp
     REAL(dp), PARAMETER :: XWA    = 28.966_dp

!     ! functions
!
!     REAL(dp) :: SCRHET ! critical ice saturation ratio for hom. freezing
!     REAL(dp) :: TAUG   ! dimensionless growth time scale
!     REAL(dp) :: PISAT  ! vapor pressure over ice [hPa]

     IFRZ    = 0
     CVF     = 0.99_dp
     CTAU    = 50.0_dp
     DSCRMIN = 0.05_dp

     ! ***** DO NOTHING CRITERIA

     CI    = 0.0_dp
     RI    = 0.0_dp
     XMI   = 0.0_dp
     CTOT  = 0.0_dp
     DO N  = 1, NFRZMOD
        CTOT = CTOT + C(N)
     END DO
     IF ((CTOT.LE.0.0_dp).OR.(V.LE.0.0_dp)) RETURN
     NFRZ  = IFRZ

     ! ***** CHECK IF FREEZING CONDITIONS ARE MET WITHIN DT

     COOLR  = GRAV * V     / CPAIR
     TMIN   = T    - COOLR * DT
     TMIN   = MAX( TMIN, 170.0_dp )
     SCR    = SCRHET(TMIN,NFRZ)
     PW     = SI * PISAT(T)
     PWCR   = PW * ( TMIN / T )**3.5_dp
     IF ((PWCR/PISAT(TMIN)).LT.SCR) RETURN

     ! ***** FREEZING TEMPERATURE

     SCR    = SCRHET(T,NFRZ)
     IF (SI.GE.SCR) THEN
        SISVE = SCR
        TEMP  = T
        PICE  = PISAT(TEMP)
        PW    = SCR * PICE
        PWCR  = PW
        GOTO 10
     ELSE
        SISVE = SI
        SL    = ( T - 170.0_dp ) / REAL(KMAX,dp)
        DO K = 1, KMAX+1
           TEMP = T  - SL * REAL(K-1,dp)
           SCR  = SCRHET(TEMP,NFRZ)
           PICE = PISAT(TEMP)
           PWCR = PW * (TEMP/T)**3.5_dp
           IF ((PWCR/PICE).GE.SCR) GOTO 10
        END DO ! K-loop
     END IF

     RETURN

 10  CONTINUE

     VOLF    = PI * RHOICE / 0.75_dp
     PCR     = P  * PWCR   / PW
     SCR     = MAX( SCR, 1.001_dp )
     SUPSCR  = SCR - 1.0_dp + DSCRMIN
     DO N    = 1, NFRZMOD
        CCR(N) = C(N) * PWCR / PW
     END DO ! N-loop

     ! ***** TIMESCALE OF THE FREEZING EVENT (FROM HOMOGENEOUS RATE)

     DLNJDT = ABS( 4.37_dp - 0.03_dp*TEMP )
     TAU    = 1.0_dp / ( CTAU*DLNJDT*COOLR )

     ! ***** ICE CRYSTAL PROPERTIES AFTER THE FREEZING EVENT

     DIFFC  = 4.0122E-3_dp * TEMP**1.94_dp / PCR
     BKT    = BK   * TEMP
     VTH    = SQRT( 8.0_dp*BKT / (PI*XMW) )
     CISAT  = 1.0E3_dp * PICE / BKT
     THETA  = HEAT * XWW / (RGAS*TEMP)
     A(1)   = ( THETA/CPAIR - XWA/RGAS ) * ( GRAV/TEMP )
     A(2)   = 1.0_dp   / CISAT
     A(3)   = 0.001_dp*XWW**2*HEAT**2 / ( AVOG*XWA*PCR*TEMP*CPAIR )
     B(1)   = SVOL * 0.25_dp  * ALPHA * VTH * CISAT * SUPSCR
     B(2)   = 0.25_dp * ALPHA * VTH   / DIFFC

     CALL XICEHET (NFRZMOD,V,TEMP,TAU,SCR,A,B,CCR,R,SIG,CI,RIHAT,RS,YK)

     XMI0   = VOLF * CI  * RIHAT**3  * CVF
     XMIMAX = XMI0 + XMW * CISAT * ( SCR - 1.0_dp )
     RIMAX  = ( XMIMAX / (VOLF*CI) )**THIRD

     ! ***** VAPOR RELAXATION: ICE CRYSTAL PROPERTIES AFTER DT

     XMISAT = XMW   * CISAT
     TGROW  = 0.75_dp  / ( PI*DIFFC*CI*RIMAX )
     ZF     = TGROW / DT
     XMFP   = 3.0_dp    * DIFFC / VTH
     BETA   = XMFP  / ( 0.75_dp*ALPHA*RIMAX )
     X0     = RIHAT / RIMAX

     X = 1.0_dp
     DO WHILE (X >= X0)
        Z     = ZF * TAUG(BETA,X,X0)
        IF (Z.LE.1.0_dp) EXIT
        X = X - 0.01_dp
     END DO

     RI     = X    * RIMAX
     XMI    = VOLF * CI * RI**3
     SI     = SCR  - ( XMI - XMI0 )/XMISAT

     RETURN

  END SUBROUTINE XFRZHET

! ------------------------------------------------------------------------------

  SUBROUTINE XICEHET (NFRZMOD,V,T,TAU,SCR,A,B,C,R,SIG,CI,RIHAT,RS,YK)

    ! ***** ICE CRYSTAL CONCENTRATION AND SIZE AFTER FREEZING EVENT

    IMPLICIT NONE

    ! subroutine parameters

    INTEGER,  INTENT(in)  :: NFRZMOD      ! number of aerosol modes
    REAL(dp), INTENT(in)  :: V           ! updraft [cm/s]
    REAL(dp), INTENT(in)  :: T           ! temperature [K]
    REAL(dp), INTENT(in)  :: TAU         ! time scale of the freezing event
    REAL(dp), INTENT(in)  :: SCR         ! 
    REAL(dp), INTENT(in)  :: A(3)        ! 
    REAL(dp), INTENT(in)  :: B(2)        ! 
    REAL(dp), INTENT(in)  :: C(nfrzmod)   ! 
    REAL(dp), INTENT(in)  :: R(nfrzmod)   ! aerosol radius [cm]
    REAL(dp), INTENT(in)  :: SIG(nfrzmod) ! aerosol standard deviation
                                         ! (log-normal distribution)
    REAL(dp), INTENT(out) :: CI          ! ice nuclei
    REAL(dp), INTENT(out) :: RIHAT       ! 
    REAL(dp), INTENT(out) :: RS          ! 
    REAL(dp), INTENT(out) :: YK          ! 

    ! parameters

    INTEGER,  PARAMETER :: IBIN   = 80

    REAL(dp), PARAMETER :: RMIN   = 1.0E-7_dp
    REAL(dp), PARAMETER :: RMAX   = 1.0E-3_dp
    REAL(dp), PARAMETER :: VRAT   = 1.5_dp
    REAL(dp), PARAMETER :: PI     = 3.1415927_dp
    REAL(dp), PARAMETER :: SQPI   = 1.7724539_dp
    REAL(dp), PARAMETER :: TWOPI  = 6.2831853_dp
    REAL(dp), PARAMETER :: SVOL   = 3.23E-23_dp
    REAL(dp), PARAMETER :: THIRD  = 0.3333333_dp
    REAL(dp), PARAMETER :: XMW    = 2.992E-23_dp
    REAL(dp), PARAMETER :: RHOICE = 0.925_dp

    ! local variables

    REAL(dp) :: SHAPE(IBIN),   R0(IBIN+1)
    REAL(dp) :: CRINT(IBIN+1), CIINT(IBIN+1)
    REAL(dp) :: VOLF,PHI,RIMFC,TBBT,CTOT,DELTA,DELP1,SYK,EERFC,RIM,ARG,RMEAN
    REAL(dp) :: SLOPE,SLOPER,RLOGRAT,CSH,SHAPFC1,SHAPFC2,SIGL,RIRAT,SUMSH
    REAL(dp) :: XMIHAT,DLOGR0

    INTEGER :: NBIN,N,I,II

!    ! functions
!
!    REAL(dp) :: XEERFC

    ! ***** CONSTANTS

    NBIN  = 1 + INT( LOG( (RMAX/RMIN)**3 ) / LOG(VRAT) )
    VOLF  = PI * RHOICE    / 0.75_dp
    PHI   = V  * A(1)*SCR / ( A(2) + A(3)*SCR )
    RIMFC = 4.0_dp * PI   * B(1)/B(2)**2 / SVOL
    TBBT  = 2.0_dp * B(1) * B(2) * TAU
    CTOT  = 0.0_dp
    DO N  = 1, NFRZMOD
       CTOT = CTOT + C(N)
    END DO

    ! ***** MONODISPERSE APPROXIMATION (SINGLE MODE ONLY)

    IF ((NFRZMOD.EQ.1).AND.(SIG(NFRZMOD).LT.1.1)) THEN
       RS     = R(NFRZMOD)
       DELTA  = B(2) * RS
       DELP1  = 1.0_dp   + DELTA
       YK     = TBBT / DELP1**2
       SYK    = SQRT(YK)
       EERFC  = XEERFC(YK)
       RIM    = RIMFC / DELP1 * ( DELTA**2 - 1.0_dp &
                                  + (1.0_dp+0.5_dp*YK*DELP1**2)*SQPI*EERFC/SYK )
       CI     = PHI / RIM
       RIRAT  = 1.0_dp + 0.5_dp * SQPI * SYK * EERFC
       RIHAT  = ( RIRAT * DELP1 - 1.0_dp ) / B(2)
       XMIHAT = VOLF * CI * RIHAT**3
       CI     = MIN( CI, CTOT )
       RIHAT  = ( XMIHAT / (VOLF*CI) )**THIRD
       RETURN
    END IF

    ! ***** SIZE DISTRIBUTION PROPERTIES

    R0(NBIN+1)    = RMAX * VRAT**THIRD
    CIINT(NBIN+1) = 1.0E-35_dp
    CRINT(NBIN+1) = 1.0E-25_dp
    DLOGR0     = 2.0_dp**THIRD * (VRAT**THIRD-1.0_dp) / (VRAT+1.0_dp)**THIRD
    SUMSH      = 0.0_dp
    DO I       = 1, NBIN
       CIINT(I)  = 0.0_dp
       CRINT(I)  = 0.0_dp
       SHAPE(I)  = 0.0_dp
       R0(I)     = RMIN * VRAT**( THIRD*REAL(I-1,dp) )
       DO N      = 1, NFRZMOD
          SIGL     = LOG( MAX( SIG(N), 1.1_dp ) )
          SHAPFC1  = 1.0_dp   / ( SQRT(TWOPI) * SIGL )
          SHAPFC2  = 0.5_dp / SIGL**2
          ARG      = SHAPFC2  * ( LOG(R0(I)/R(N)) )**2
          ARG      = MIN( ARG, 75.0_dp )
          SHAPE(I) = SHAPE(I) + DLOGR0 * C(N) * SHAPFC1 * EXP( -ARG )
       END DO ! N-loop
       SUMSH     = SUMSH    + SHAPE(I)
    END DO ! I-loop

    ! ***** ICE CRYSTAL PROPERTIES

    RMEAN     = 0.0_dp
    DO I  = NBIN, 1, -1
       DELTA    = B(2) * R0(I)
       DELP1    = 1.0_dp   + DELTA
       YK       = TBBT / DELP1**2
       SYK      = SQRT(YK)
       EERFC    = XEERFC(YK)
       RIM      = RIMFC / DELP1 * ( DELTA**2 - 1.0_dp &
                                  + (1.0_dp+0.5_dp*YK*DELP1**2)*SQPI*EERFC/SYK )
       CSH      = SHAPE(I) / SUMSH   * CTOT
       CRINT(I) = CRINT(I+1) + RIM   * CSH
       CIINT(I) = CIINT(I+1) +         CSH
       RMEAN    = RMEAN      + R0(I) * CSH
       IF (CRINT(I).GE.PHI) GOTO 10
    END DO ! I-loop
    RS        = R0(1)
    RMEAN     = RMEAN / CTOT
    CI        = CTOT  * PHI   / CRINT(1)
    DELP1     = 1.0_dp    + B(2)  * RMEAN
    YK        = TBBT  / DELP1**2
    SYK       = SQRT(YK)
    EERFC     = XEERFC(YK)
    RIRAT     = 1.0_dp + 0.5_dp * SQPI * SYK * EERFC
    RIHAT     = ( RIRAT * DELP1 - 1.0_dp ) / B(2)
    XMIHAT    = VOLF    * CI * RIHAT**3
    CI        = CTOT
    RIHAT     = ( XMIHAT / (VOLF*CI) )**THIRD
    RETURN
 10 CONTINUE
    II        = MAX( I, 1 )
    RLOGRAT   = LOG( R0(II) / R0(II+1) )
    SLOPER    = LOG( CRINT(II) / CRINT(II+1) ) / RLOGRAT
    RS        = R0(II+1)    * ( PHI / CRINT(II+1) )**(1.0_dp/SLOPER)
    SLOPE     = LOG( CIINT(II) / CIINT(II+1) ) / RLOGRAT
    CI        = CIINT(II+1) * ( RS / R0(II+1) )**SLOPE
    RMEAN     = RMEAN / CTOT
    DELP1     = 1.0_dp   + B(2) * MAX( RS, RMEAN )
    YK        = TBBT / DELP1**2
    SYK       = SQRT(YK)
    EERFC     = XEERFC(YK)
    RIRAT     = 1.0_dp + 0.5_dp * SQPI * SYK * EERFC
    RIHAT     = ( RIRAT * DELP1 - 1.0_dp ) / B(2)

    ! ***** EXIT
    RETURN

  END SUBROUTINE XICEHET

! ------------------------------------------------------------------------------

  SUBROUTINE XFRZHOM (NOSIZE,nfrzmod,DT,C,R,SIG,P,T,V,SI,CI,RI,  &
                      COOLR, PW, SCR, TEMP, PICE, PWCR)

     ! ***** HOMOGENEOUS FREEZING OF SUPERCOOLED AEROSOL (ADIABATIC ASCENT)
     !
     ! ***** BERND K??RCHER  APR 30  2002
     ! ***** bernd.kaercher@dlr.de  http://www.op.dlr.de/~pa1c/
     ! ***** Ulrike Lohmann  Dalhousie Univ  Apr 02  Ulrike.Lohmann@Dal.Ca
     !
     ! ***** References
     ! K??rcher, B. and U. Lohmann
     ! A parameterization of cirrus cloud formation:
     ! Homogeneous freezing including effects of aerosol size.
     ! J. Geophys. Res. 107 (D ), in press, 2002.

     IMPLICIT NONE

     ! subroutine parameters
 
     LOGICAL,  INTENT(in)  :: NOSIZE      ! take into account aerosol size?
     INTEGER,  INTENT(in)  :: NFRZMOD      ! number of aerosol modes
     REAL(dp), INTENT(in)  :: DT          ! time step [s]
     REAL(dp), INTENT(in)  :: C(NFRZMOD)   !
     REAL(dp), INTENT(in)  :: R(NFRZMOD)   ! aerosol radius [cm]
     REAL(dp), INTENT(in)  :: SIG(NFRZMOD) ! aerosol standard deviation
                                          ! (log-normal size distribution)
     REAL(dp), INTENT(in)  :: P           ! pressure [hPa]
     REAL(dp), INTENT(in)  :: T           ! temperature [K]
     REAL(dp), INTENT(in)  :: V           ! updraft [cm/s]
     REAL(dp), INTENT(in)  :: COOLR       ! 
     REAL(dp), INTENT(in)  :: PW          ! 
     REAL(dp), INTENT(in)  :: SCR         ! 
     REAL(dp), INTENT(in)  :: TEMP        ! 
     REAL(dp), INTENT(in)  :: PICE        ! 
     REAL(dp), INTENT(in)  :: PWCR        ! 

     REAL(dp), INTENT(out) :: SI          ! saturation pressure (ice)
     REAL(dp), INTENT(out) :: CI          ! ice nuclei
     REAL(dp), INTENT(out) :: RI          ! 

     ! parameters

     INTEGER,  PARAMETER :: KMAX   = 120

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

     ! local variables

     REAL(dp) :: A(3)        ! 
     REAL(dp) :: B(2)        ! 
     REAL(dp) :: CCR(NFRZMOD) ! 

     REAL(dp) :: XMI,VOLF,PCR,CTAU,DLNJDT,TAU,DIFFC,BKT,VTH,CISAT,THETA
     REAL(dp) :: YK,RS,RIHAT,XMI0,XMIMAX,RIMAX,XMISAT,TGROW,ZF,XMFP
     REAL(dp) :: BETA,X0,X,Z

     INTEGER  :: N

!     ! functions
!
!     REAL(dp) :: PISAT  ! vapor pressure over ice [hPa]
!     REAL(dp) :: TAUG   ! dimensionless growth time scale

     ! ***** DO NOTHING CRITERIA

     CI    = 0.0_dp
     RI    = 0.0_dp
     XMI   = 0.0_dp

     ! ***** FREEZING TEMPERATURE

!!$CC      SCR    = SCRHOM(T)
!!$CC      IF (SI.GE.SCR) THEN
!!$CC       SISVE = SCR
!!$CC       TEMP  = T
!!$CC       PICE  = PISAT(TEMP)
!!$CC       PW    = SCR * PICE
!!$CC       PWCR  = PW
!!$CC       GOTO 10
!!$CC      ELSE
!!$CC       SISVE = SI
!!$CC       SL    = ( T - 170. ) / REAL(KMAX,dp)
!!$CC       DO  K = 1, KMAX+1
!!$CC        TEMP = T  - SL * REAL(K-1,dp)
!!$CC        SCR  = SCRHOM(TEMP)
!!$CC        PICE = PISAT(TEMP)
!!$CC        PWCR = PW * (TEMP/T)**3.5
!!$CC        IF ((PWCR/PICE).GE.SCR) GOTO 10
!!$CC       ENDDO
!!$CC      ENDIF
!!$CC      RETURN
!!$CC 10   CONTINUE

     VOLF    = PI * RHOICE / 0.75_dp
     PCR     = P  * PWCR   / PW
     DO N    = 1, NFRZMOD
        CCR(N) = C(N) * PWCR / PW
     END DO

     ! ***** TIMESCALE OF THE FREEZING EVENT

     IF (NOSIZE) THEN
        CTAU  = MAX( (2260.0_dp-10.0_dp*TEMP), 100.0_dp )
     ELSE
        CTAU  = 50.0_dp
     END IF
     DLNJDT = ABS( 4.37_dp - 0.03_dp*TEMP )
     TAU    = 1.0_dp / ( CTAU*DLNJDT*COOLR )

     ! ***** ICE CRYSTAL PROPERTIES AFTER THE FREEZING EVENT

     DIFFC  = 4.0122E-3_dp * TEMP**1.94_dp / PCR
     BKT    = BK   * TEMP
     VTH    = SQRT( 8.0_dp*BKT / (PI*XMW) )
     CISAT  = 1.E3_dp * PICE / BKT
     THETA  = HEAT * XWW / (RGAS*TEMP)
     A(1)   = ( THETA/CPAIR - XWA/RGAS ) * ( GRAV/TEMP )
     A(2)   = 1.0_dp   / CISAT
     A(3)   = 0.001_dp*XWW**2*HEAT**2 / ( AVOG*XWA*PCR*TEMP*CPAIR )
     B(1)   = SVOL * 0.25_dp  * ALPHA * VTH * CISAT * ( SCR - 1.0_dp )
     B(2)   = 0.25_dp * ALPHA * VTH   / DIFFC

     CALL XICEHOM (NOSIZE,NFRZMOD,V,TEMP,TAU,SCR,A,B,CCR,R,SIG,CI,RIHAT,RS,YK)

     XMI0   = VOLF * CI  * RIHAT**3
     XMIMAX = XMI0 + XMW * CISAT * ( SCR - 1.0_dp )
     RIMAX  = ( XMIMAX / (VOLF*CI) )**THIRD

     ! ***** VAPOR RELAXATION: ICE CRYSTAL PROPERTIES AFTER DT

     XMISAT = XMW   * CISAT
     TGROW  = 0.75_dp  / ( PI*DIFFC*CI*RIMAX )
     ZF     = TGROW / DT
     XMFP   = 3.0_dp    * DIFFC / VTH
     BETA   = XMFP  / ( 0.75_dp*ALPHA*RIMAX )
     X0     = RIHAT / RIMAX
!     DO X   = 1.0_dp, X0, -0.01_dp
     X = 1.0_dp
     DO WHILE (X >= X0)
        Z = ZF * TAUG(BETA,X,X0)
        IF (Z.LE.1.0_dp) EXIT
        X = X - 0.01_dp
     END DO

     RI     = X    * RIMAX
     XMI    = VOLF * CI * RI**3
     SI     = SCR  - ( XMI - XMI0 )/XMISAT

     RETURN

  END SUBROUTINE XFRZHOM

! ------------------------------------------------------------------------------

  SUBROUTINE XICEHOM (NOSIZE,NFRZMOD,V,T,TAU,SCR,A,B,C,R,SIG,CI,RIHAT,RS,YK)

     ! ***** ICE CRYSTAL CONCENTRATION AND SIZE AFTER FREEZING EVENT

     IMPLICIT NONE

     ! subroutine parameters

     LOGICAL,  INTENT(in)  :: NOSIZE      ! take into account aerosol size?
     INTEGER,  INTENT(in)  :: NFRZMOD      ! number of aerosol modes
     REAL(dp), INTENT(in)  :: V           ! updraft [cm/s]
     REAL(dp), INTENT(in)  :: T           ! temperature [K]
     REAL(dp), INTENT(in)  :: TAU         ! 
     REAL(dp), INTENT(in)  :: SCR         ! 
     REAL(dp), INTENT(in)  :: A(3)        ! 
     REAL(dp), INTENT(in)  :: B(2)        ! 
     REAL(dp), INTENT(in)  :: C(NFRZMOD)
     REAL(dp), INTENT(in)  :: R(NFRZMOD)   ! aerosol radius [cm]
     REAL(dp), INTENT(in)  :: SIG(NFRZMOD) ! aerosol standard deviation
                                          ! (log-normal size distribution)
     REAL(dp), INTENT(out) :: CI          ! ice nuclei
     REAL(dp), INTENT(out) :: RIHAT       ! 
     REAL(dp), INTENT(out) :: RS          ! 
     REAL(dp), INTENT(out) :: YK          ! 

     ! parameters

     INTEGER,  PARAMETER :: IBIN   = 80
     REAL(dp), PARAMETER :: RMIN   = 1.0E-7_dp
     REAL(dp), PARAMETER :: RMAX   = 1.0E-3_dp
     REAL(dp), PARAMETER :: VRAT   = 1.5_dp
     REAL(dp), PARAMETER :: PI     = 3.1415927_dp
     REAL(dp), PARAMETER :: SQPI   = 1.7724539_dp
     REAL(dp), PARAMETER :: TWOPI  = 6.2831853_dp
     REAL(dp), PARAMETER :: SVOL   = 3.23E-23_dp
     REAL(dp), PARAMETER :: THIRD  = 0.3333333_dp
     REAL(dp), PARAMETER :: XMW    = 2.992E-23_dp
     REAL(dp), PARAMETER :: RHOICE = 0.925_dp

     ! local variables

     REAL(dp) :: SHAPE(IBIN), R0(IBIN+1)
     REAL(dp) :: CRINT(IBIN+1), CIINT(IBIN+1)
     REAL(dp) :: RIMX(IBIN+1), CSHX(IBIN+1)

     REAL(dp) :: VOLF,PHI,RIMFC,TBBT,CTOT,XMIHAT,DELTA,DELP1,SYK,RIM,RIRAT
     REAL(dp) :: DLOGR0,SUMSH,SIGL,SHAPFC1,SHAPFC2,ARG,RMEAN,CSH,RLOGRAT
     REAL(dp) :: SLOPER,SLOPE,EERFC

     INTEGER  :: NBIN,N,I,II

!     ! functions
!
!     REAL(dp) :: XEERFC

     ! ***** CONSTANTS

     NBIN  = 1 + INT( LOG( (RMAX/RMIN)**3 ) / LOG(VRAT) )
     VOLF  = PI * RHOICE    / 0.75_dp
     PHI   = V  * A(1)*SCR / ( A(2) + A(3)*SCR )
     RIMFC = 4.0_dp * PI   * B(1)/B(2)**2 / SVOL
     TBBT  = 2.0_dp * B(1) * B(2) * TAU
     CTOT  = 0.0_dp
     DO N  = 1, NFRZMOD
        CTOT = CTOT + C(N)
     END DO

     ! ***** NO SIZE EFFECTS

     IF (NOSIZE) THEN
        RS     = 0.25E-4_dp
        YK     = TAU  * ( B(1)/RS ) / ( 1.0_dp + B(2)*RS )
        CI     = SVOL * ( B(2) / (TWOPI*B(1)) )**1.5_dp * PHI / SQRT(TAU)
        CI     = MIN( CI, CTOT )
        XMIHAT = XMW * PI * PHI * TAU / 6.0_dp
        RIHAT  = ( XMIHAT / (VOLF*CI) )**THIRD
        RETURN
     END IF

     ! ***** MONODISPERSE AEROSOL (SINGLE MODE ONLY)

     IF ((NFRZMOD.EQ.1).AND.(SIG(NFRZMOD).LT.1.1)) THEN
        RS     = R(NFRZMOD)
        DELTA  = B(2) * RS
        DELP1  = 1.0_dp   + DELTA
        YK     = TBBT / DELP1**2
        SYK    = SQRT(YK)
        EERFC  = XEERFC(YK)
        RIM    = RIMFC / DELP1 * ( DELTA**2 - 1.0_dp   &
                                  + (1.0_dp+0.5_dp*YK*DELP1**2)*SQPI*EERFC/SYK )
        CI     = PHI / RIM
        RIRAT  = 1.0_dp + 0.5_dp * SQPI * SYK * EERFC
        RIHAT  = ( RIRAT * DELP1 - 1.0_dp ) / B(2)
        XMIHAT = VOLF * CI * RIHAT**3
        CI     = MIN( CI, CTOT )
        RIHAT  = ( XMIHAT / (VOLF*CI) )**THIRD
        RETURN
     END IF

     ! ***** SIZE DISTRIBUTION PROPERTIES

     R0(NBIN+1)    = RMAX * VRAT**THIRD
     CIINT(NBIN+1) = 1.E-35_dp
     CRINT(NBIN+1) = 1.E-25_dp
     DLOGR0     = 2.0_dp**THIRD * (VRAT**THIRD-1.0_dp) / (VRAT+1.0_dp)**THIRD
     SUMSH      = 0.0_dp
     DO I       = 1, NBIN
        CIINT(I)  = 0.0_dp
        CRINT(I)  = 0.0_dp
        SHAPE(I)  = 0.0_dp
        R0(I)     = RMIN * VRAT**( THIRD*REAL(I-1,dp) )
        DO N      = 1, NFRZMOD
           SIGL     = LOG( MAX( SIG(N), 1.1_dp ) )
           SHAPFC1  = 1.0_dp   / ( SQRT(TWOPI) * SIGL )
           SHAPFC2  = 0.5_dp / SIGL**2
           ARG      = SHAPFC2  * ( LOG(R0(I)/R(N)) )**2
           ARG      = MIN( ARG, 75.0_dp )
           SHAPE(I) = SHAPE(I) + DLOGR0 * C(N) * SHAPFC1 * EXP( -ARG )
        END DO
        SUMSH     = SUMSH    + SHAPE(I)
     END DO

     ! ***** ICE CRYSTAL PROPERTIES

     RMEAN     = 0.0_dp
     DO     I  = NBIN, 1, -1
        DELTA    = B(2) * R0(I)
        DELP1    = 1.0_dp   + DELTA
        YK       = TBBT / DELP1**2
        SYK      = SQRT(YK)
        EERFC    = XEERFC(YK)
        RIMX(I)  = RIMFC / DELP1 * ( DELTA**2 - 1.0_dp  &
                                  + (1.0_dp+0.5_dp*YK*DELP1**2)*SQPI*EERFC/SYK )
        CSHX(I)  = SHAPE(I) / SUMSH   * CTOT
     END DO
     DO     I  = NBIN, 1, -1
        RIM = RIMX(I)
        CSH = CSHX(I)
        CRINT(I) = CRINT(I+1) + RIM   * CSH
        CIINT(I) = CIINT(I+1) +         CSH
        RMEAN    = RMEAN      + R0(I) * CSH
        IF (CRINT(I).GE.PHI) GOTO 10
     END DO
     RS        = R0(1)
     RMEAN     = RMEAN / CTOT
     CI        = CTOT  * PHI   / CRINT(1)
     DELP1     = 1.0_dp    + B(2)  * RMEAN
     YK        = TBBT  / DELP1**2
     SYK       = SQRT(YK)
     EERFC     = XEERFC(YK)
     RIRAT     = 1.0_dp + 0.5_dp * SQPI * SYK * EERFC
     RIHAT     = ( RIRAT * DELP1 - 1.0_dp ) / B(2)
     XMIHAT    = VOLF    * CI * RIHAT**3
     CI        = CTOT
     RIHAT     = ( XMIHAT / (VOLF*CI) )**THIRD

     RETURN

 10  CONTINUE

     II        = MAX( I, 1 )
     RLOGRAT   = LOG( R0(II) / R0(II+1) )
     SLOPER    = LOG( CRINT(II) / CRINT(II+1) ) / RLOGRAT
     RS        = R0(II+1)    * ( PHI / CRINT(II+1) )**(1.0_dp/SLOPER)
     SLOPE     = LOG( CIINT(II) / CIINT(II+1) ) / RLOGRAT
     CI        = CIINT(II+1) * ( RS / R0(II+1) )**SLOPE
     RMEAN     = RMEAN / CTOT
     DELP1     = 1.0_dp   + B(2) * MAX( RS, RMEAN )
     YK        = TBBT / DELP1**2
     SYK       = SQRT(YK)
     EERFC     = XEERFC(YK)
     RIRAT     = 1.0_dp + 0.5_dp * SQPI * SYK * EERFC
     RIHAT     = ( RIRAT * DELP1 - 1.0_dp ) / B(2)

     ! **** EXIT
     RETURN

  END SUBROUTINE XICEHOM

! ------------------------------------------------------------------------------

  SUBROUTINE XVICE(TIME, V, P, T, R, C, S, VICE)

    !adjust updraft speed because of pre-existed ice

    !*****       V [CM/S]        CONSTANT VERTICAL WIND SPEED > 0
    !*****       T [K], P [MB]   INITIAL AIR TEMPERATURE, PRESSURE
    !*****       R [CM]       INITIAL MEAN  ICE CRYSTAL NUMBER RADIUS
    !*****       C [CM-3]     INITIAL TOTAL ICE CRYSTAL NUMBER DENSITY
    !*****       S               INITIAL ICE SATURATION RATIO

    !*****       VICE [CM/S]     EFFECTIVE VERTICAL VELOCITY CALCULATED
    !*****                       FROM DEPOSITIONAL LOSS TERM

    !*****       RICE(*) [CM]    FINAL MEAN ICE CRYSTAL NUMBER RADIUS
    !*****       CICE(*) [CM-3]  FINAL MEAN ICE CRYSTAL NUMBER DENSITY
    !*****       SICE            FINAL ICE SATURATION RATIO

    IMPLICIT NONE

    ! subroutine parameters

    REAL(dp), INTENT(in)  :: TIME        ! 
    REAL(dp), INTENT(in)  :: V           ! updraft [cm/s]
    REAL(dp), INTENT(in)  :: P           ! pressure
    REAL(dp), INTENT(in)  :: T           ! temperature [K]
    REAL(dp), INTENT(in)  :: R           ! initial mean ice crystal number radius
    REAL(dp), INTENT(in)  :: C           ! initial total ice crystal number density
    REAL(dp), INTENT(in)  :: S           ! initial super saturation ratio

    REAL(dp), INTENT(out) :: VICE        ! effective vertical velocity

    ! parameters


    REAL(dp), PARAMETER :: ALPHA = 0.5_dp
    REAL(dp), PARAMETER :: THOUBK = 7.24637701E+18_dp
    REAL(dp), PARAMETER :: FPIVOL = 3.89051704E+23_dp
    REAL(dp), PARAMETER :: FA1 = 0.601272523_dp
    REAL(dp), PARAMETER :: FA2 = 0.000342181855_dp
    REAL(dp), PARAMETER :: FA3 = 1.49236645E-12_dp
    REAL(dp), PARAMETER :: FC  = 9.80999976E-05_dp
    REAL(dp), PARAMETER :: FD = 249.239822_dp
    REAL(dp), PARAMETER :: FVTH = 11713803._dp
    REAL(dp), PARAMETER :: SVOL = 3.23E-23_dp
    REAL(dp), PARAMETER :: WVP1 = 3.6E+10_dp
    REAL(dp), PARAMETER :: WVP2 = 6145._dp

    REAL(dp) :: ALP4, T1, COOLR, TEMP, TEMP1, ADCHG, PRESS, FLUX, CISAT
    REAL(dp) :: PICE, RICE, CICE, SICE
    REAL(dp) :: A1, A2, A3, B1, B2, DLOSS, DELP1, ETA

    !*****SET CONSTANTS
      ALP4  = 0.25 * ALPHA
      T1    = 1./ T
      COOLR = FC * V

    !*****SET INITIAL VALUES
      RICE = R
      SICE = S

    !*****COMPUTE SATURATION RATIO INCREMENT

      TEMP  = T - TIME * COOLR
      TEMP1 = 1./ TEMP
      PICE  = WVP1 * EXP(-(WVP2*TEMP1))
      ADCHG = ( TEMP * T1 )**2.5
      PRESS = P * ( TEMP * T1 )**3.5
      FLUX  = ALP4 * SQRT(FVTH*TEMP)
      CISAT = THOUBK * PICE * TEMP1
      A1    = ( FA1 * TEMP1 - FA2 ) * TEMP1
      A2    = 1./ CISAT
      A3    = FA3 * TEMP1 / PRESS
      B1    = FLUX * SVOL * CISAT * ( SICE-1. )
      B2    = FLUX * FD * PRESS * TEMP1**1.94

    !*****SUM UP DEPOSITION TERM
      DLOSS     = 0.
      IF (C.GT.0.) THEN
       CICE = ADCHG * C
       DELP1    = 1.+ B2 * RICE
       ETA      = 1.+ 2.* B1 * B2 * TIME / DELP1**2
       RICE = ( DELP1 * SQRT(ETA) - 1. ) / B2
       DLOSS    = DLOSS + FPIVOL * CICE                  &
                    * B1 * RICE**2 / ( 1.+ B2 * RICE)
      ENDIF

      VICE      = ( A2 + A3 * SICE ) * DLOSS / ( A1 * SICE )

      RETURN

  END SUBROUTINE XVICE
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
 SUBROUTINE cloud_read_nml_ctrl_l10(iou, status, modstr)

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

   ! I/O
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit
    INTEGER, INTENT(OUT) :: status ! error status
    CHARACTER(LEN=*), INTENT(IN) :: modstr

    NAMELIST /CTRL_L10/ cdncmin, limm_BN09

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'cloud_read_nml_ctrl_l10'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    ! INITIALIZE GLOBAL CONTROL VARIABLES
    ! -> DEFAULT VALUES ARE SET AT DECLARATION ABOVE

    CALL read_nml_open(lex, substr, iou, 'CTRL_L10', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL_L10, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL_L10', modstr)

    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE cloud_read_nml_ctrl_l10
! ------------------------------------------------------------------------------

END MODULE MESSY_CLOUD_LOHMANN10
