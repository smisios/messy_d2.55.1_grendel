MODULE MESSY_CLOUD_Lohmann07

!       Author:  Holger Tost
!       last modified: 26.06.2007

  USE messy_main_constants_mem,      ONLY: dp, ceffmin, ceffmax, ccwmin, &
                                           rd, rv, vtmpc1, vtmpc2,       &
                                           cpd => cp_air, rhoh2o => rho_h2o
  USE messy_cloud_mem

  IMPLICIT NONE
  SAVE
  
  INTEGER, PUBLIC :: idt_ncs 
  INTEGER, PUBLIC :: idt_nas
  INTEGER, PUBLIC :: idt_nks
  INTEGER, PUBLIC :: idt_nci 
  INTEGER, PUBLIC :: idt_nai
  INTEGER, PUBLIC :: idt_nki

  INTEGER, PUBLIC :: idt_moccs
  INTEGER, PUBLIC :: idt_mocas
  INTEGER, PUBLIC :: idt_mocks
  INTEGER, PUBLIC :: idt_mocki
  INTEGER, PUBLIC :: idt_mbccs
  INTEGER, PUBLIC :: idt_mbcas
  INTEGER, PUBLIC :: idt_mbcks
  INTEGER, PUBLIC :: idt_mbcki
  INTEGER, PUBLIC :: idt_ms4ks
  INTEGER, PUBLIC :: idt_ms4as
  INTEGER, PUBLIC :: idt_ms4cs
  INTEGER, PUBLIC :: idt_mduas
  INTEGER, PUBLIC :: idt_mducs
  INTEGER, PUBLIC :: idt_mssas
  INTEGER, PUBLIC :: idt_msscs

  INTEGER, PUBLIC :: nauto = 2
  INTEGER, PUBLIC :: ncvmicro = 0

  PRIVATE   

  PUBLIC :: cloud_cdnc, cloud_cdnc_icnc

  REAL(dp),PARAMETER :: cpbot  = 50000.0_dp  ! max. pressure level for 

! Warning: lohmann07 is not yet coupled to other aerosol submodels!
  REAL(dp),PARAMETER :: &
    crdiv(4)=(/ 0.0005E-4_dp, 0.005E-4_dp, 0.05E-4_dp, 0.5E-4_dp /)

CONTAINS


!========================================================

  SUBROUTINE cloud_cdnc (    kproma, kbdim, ktdia, klev, klevp1, ztmst &
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
                         ! mz_ht_20070629+
                         , status_string,      lcover                  &
                         , slm,      glac,     pcdncact                &
                         , zna                                         &
                         , plwc,     piwc                              &
                         , pfrain,   pfsnow                            &
                         , pfrain_no,pfsnow_no                         &
                         , prate_r,  prate_s                           &
                         , prevap,   pssubl                            &
! mz_ak_20051221 zcond added
                         , pr_cover, pcond                             &
                         , pimelt,   pisedi                            &
                         ! mz_ht_20070629-
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
!
!

    USE messy_cloud_ori,  ONLY : cqtmin, tmelt, cvtfall, crhosno, cn0s &
                            , cthomi, csecfrl, ncctop, cvarmin         &
                            , cbeta_pq, cbeta_pq_max, nbetaq, cbetaqs  &
                            , rbetak, nbetax, tbetai0, tbetai1, cauloc &
                            , clmax, clmin, jbmin, jbmax, lonacc       &
                            , ccraut, crhoi, ccsaut  &
                            , ccsacl, cbeta_cs, LOOKUPOVERFLOW

    USE messy_main_constants_mem, ONLY: API => PI, g, M_air
    USE messy_main_tools,         ONLY: jptlucu1, jptlucu2, &
                                        tlucua, tlucuaw, tlucub
    

  IMPLICIT NONE
!
  INTEGER :: kbdim, klevp1, klev, kproma, ktdia, krow, ktrac, jl, jk, it

  CHARACTER(LEN=32) :: status_string
  LOGICAL  :: lcover
  REAL(dp) :: slm(kbdim), glac(kbdim)

  REAL(dp)::paphm1(kbdim,klevp1),pvervel(kbdim,klev)                   &
           ,papm1(kbdim,klev)   ,pqm1(kbdim,klev)   ,papp1(kbdim,klev) &
           ,ptm1(kbdim,klev)    ,ptvm1(kbdim,klev)  ,pxlm1(kbdim,klev) &
           ,pxim1(kbdim,klev)   ,pxtec(kbdim,klev)  ,pqtec(kbdim,klev) &
           ,pxvar(kbdim,klev)   ,pxskew(kbdim,klev)                    &
           ,pbetaa(kbdim,klev)  ,pbetab(kbdim,klev)                    &
           ,pvdiffp(kbdim,klev) ,phmixtau(kbdim,klev)                  &
           ,pvmixtau(kbdim,klev),pgeo(kbdim,klev)   ,pbetass(kbdim,klev)
  REAL(dp)::pxlvi(kbdim)        ,pxivi(kbdim)
  REAL(dp)::paclc(kbdim,klev)   ,paclcac(kbdim,klev)
  REAL(dp)::pacdnc(kbdim,klev)  ,prelhum(kbdim,klev)
  REAL(dp)::paclcov(kbdim)      ,paprl(kbdim)                          &
           ,pqvi(kbdim)         ,pssfl(kbdim)
  REAL(dp)::ptte(kbdim,klev)    ,pqte(kbdim,klev)
  REAL(dp)::pxlte(kbdim,klev)   ,pxite(kbdim,klev)
  REAL(dp)::paprs(kbdim)        ,prsfl(kbdim)

  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev)         :: plwc,     piwc
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev), TARGET :: pfrain,   pfsnow
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev)         :: pfrain_no,pfsnow_no
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev), TARGET :: prevap,   pssubl
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev), TARGET :: prate_r,  prate_s 
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev)         :: pr_cover
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev)         :: pimelt,   pisedi
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev)         :: pcond 

  REAL(dp), POINTER, DIMENSION(:) :: zrfl, zsfl, zsub, zevp, zrpr
!---Included for in-cloud scavenging (Philip Stier, 28/03/01):----------
  REAL(dp) ::    pxtm1(kbdim,klev,ktrac), pxtte(kbdim,klev,ktrac) 
!---End Included for scavenging-----------------------------------------
!
!   Temporary arrays
!
  REAL(dp):: zclcpre(kbdim)                                             &
            ,zcnd(kbdim)         ,zdep(kbdim)        ,zdp(kbdim)        &
            ,                     zxievap(kbdim)     ,zxlevap(kbdim)    &
            ,zfrl(kbdim)         ,zimlt(kbdim)       ,zsmlt(kbdim)      &
            ,                     zspr(kbdim)                           &
            ,zxlte(kbdim)        ,zxite(kbdim)       ,zxiflux(kbdim)    &
            ,zsacl(kbdim)        ,zdz(kbdim)                            &
            ,zlsdcp(kbdim)       ,zlvdcp(kbdim)      ,zximlt(kbdim)     &
            ,ztp1tmp(kbdim)      ,zqp1tmp(kbdim)     ,zxisub(kbdim)     &
            ,zxlb(kbdim)         ,zxib(kbdim)                           &
            ,zrho(kbdim,klev)    ,zclcov(kbdim)      ,zclcaux(kbdim)    &
            ,zqvi(kbdim)         ,zxlvi(kbdim)       ,zxivi(kbdim)      &
            ,zbetaqt(kbdim)      ,zwide(kbdim)                          &
            ,zbetacl(kbdim)      ,zturbvar(kbdim)    ,zturbskew(kbdim)  &
            ,zconvvar(kbdim)     ,zconvskew(kbdim)   ,zvartg(kbdim)     &
            ,zmicroskew(kbdim)   ,zgenti(kbdim)      ,zgentl(kbdim)     &
            ,zxvarte(kbdim)      ,zxskewte(kbdim)                       &
            ,zgeoh(kbdim,klevp1)
!
  REAL(dp)::zbqp1,zbbap1,zbap1,ztt,zgent,zqcdif,zqp1,ztp1,             &
            zqp1b,zskew,zbetai0,zbetai1,zskewp1,zvarp1,zifrac,zvarmx
  REAL(dp)::zdtime,zauloc
  INTEGER:: iqidx,ixidx, jb
  LOGICAL   lo,lo1,lo2,locc
  INTEGER   knvb(kbdim)

!--- Included for dust emissions (Philip Stier 10/01/02)-----------------
  INTEGER :: jrow
!--- End Included for dust emissions ------------------------------------

!--- Included for prognostic CDNC/IC scheme ----------------------------

  LOGICAL, PARAMETER :: losacl=.TRUE.! Size-dependent acretion rate
                                      ! (Lohmann, JAS, 2004)

  LOGICAL, PARAMETER :: lodetr=.TRUE. ! Convective detrainment activation

  INTEGER :: itop(kbdim,klev),         ibas(kbdim,klev),               &
             iclbas (kbdim)  

  REAL(dp):: zeps, zkap, zreffl, zqlnuccv

  REAL(dp):: zrprn(kbdim),          &
             zsacln(kbdim),            zfrln(kbdim),                   &
             zcdnc(kbdim,klev),     & !
             !zesw_2d(kbdim,klev),   & !
             zsusatw_2d(kbdim,klev),& ! supersat. with respect to water
             zcdnc_burden(kbdim),   &
             zna(kbdim,klev),       &
             zqlnuc(kbdim,klev),    & !
             zqlnuccvh(kbdim,klev)    ! Nucleated CDNC from convective
                                      ! detrainment []
  ! Temporary fields
  INTEGER  :: it1
  REAL(dp) :: zesw,  zvervc, zcdncnew, zcdnmin, zxsec, zcons, ztdif,  &
              zesi , zcfac4c,zsubi, zzepr, zxip1, zqsm1, zqst1, zrelhum, &
              zcond, zes, zlc, zf1, zraut, zrautself, zrautn, &
              zrieff, zhelp, zstokes, zstcrit, zcsacl1, zsecprod, zmdelb, &
              zepsec, zxilb, zmqp1, zxlp1, zxrp1, zxsp1, zqrho, zqsec, ztmst,&
              zcons1, zcons2, zradl, zsigmaw, zdisp, zdw0, zcsacl,  zrsnow, &
              zmw0, zmi0,   zqsw , zrcp, zsnmlt, zsnowmin, zcdi, zximelt, zicncmelt, &
              zclcstar, zdpg, zqsi, zb1, zsusati, zb2, zcoeff, zclambs,&
              zzeps, zesat, zsusatw, zdv, zast, zbst, zxifall, zal1, zal2, zxised, &
              zxim1evp, zxlm1evp,zxidt, zxldt, zxidtstar, zxldtstar, zdqsdt, &
              zqvdt, zdtdt, zdtdtstar, zdqsat, zqtau, zpp, zqq, zeta, &
              zprod, zaa, zgtp, zcor, zqsp1tmp, zoversat, zrhtest, zqcon, zdepcor, zcndcor, &
              zlcdqsdt, zzevp, zdepos, znfrl, zexm1, zexp, zrac1, zrac2, &
              zxlbold, zraccn, zrih, zris, zcolleffi, zc1, zsaut, zsaci1, zsaci2, &
              zsacl1,zsacl2, zsacl1in, zsacl2in, zscnc, zdw,  zusnow, zudrop, zviscos,&
              zrey, zlamsm, zxsp2, zself, zdt2, zdplanar, zlams2, zj, zpn, zzdrr, zzdrs, &
              zpretot, zpredel, zpresum, zxlold, zdxlcor, zxiold, zdxicor, zcdnold, &
              zdnlcor, zcdnp

  REAL(dp) ::    pcdncact(kbdim,klev),     ptkem1(kbdim,klev)
  REAL(dp) ::    pcvcbot(kbdim),           pwcape(kbdim)

  REAL(dp)            :: zmi                ! assumed mass of ice crystals with 
                                            ! corresponding volume mean radius zcri 
  REAL(dp), PARAMETER :: zcri=10.E-6_dp     ! to estimate the number of produced  
                                            ! cloud droplets from ice melting in  
                                            ! case of licnc=.FALSE. [m]=> 10 um
!--- End Included for CDNC -------------------------------------------

!--- Included for dust emissions (Philip Stier 10/01/02)-----------------
   jrow = krow
!--- End Included for dust emissions ------------------------------------

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
!
! Executable statements
!
  lookupoverflow = .FALSE.
!
!   Security parameters
!
  zepsec = 1.0e-12_dp
  zxsec  = 1.0_dp-zepsec
  zqsec  = 1.0_dp-cqtmin
!
!   Computational constants
!
! commented since imported via parameter list  
  zdtime  = ztmst / 2._dp
  zcons1 = cpd*vtmpc2
  zcons2 = 1._dp/(ztmst*g)
!--- Included for prognostic CDNC/IC scheme ----------------------------
  zeps=EPSILON(1.0_dp)

  zmi=4._dp/3._dp*zcri**3._dp*api*crhoi

  zcdnmin=40.e6_dp      !
  zradl=10.e-6_dp      !
  zsigmaw=0.28_dp
  zdisp=EXP(zsigmaw**2/2._dp)
  zdw0=10.e-6_dp*zdisp ! modal diameter times dispersion parameter
  zcdnc_burden(1:kproma) = 0.0_dp

  cdnc_burden(1:kproma,jrow) = 0.0_dp
  cdnc_acc(1:kproma,:,jrow)  = 0._dp
  cdnc(1:kproma,:,jrow)      = 0._dp
  qnuc(1:kproma,:,jrow)      = 0._dp

  zcsacl = 0.01_dp

  zqlnuc(1:kproma,:) = 0.0_dp
  zsnowmin=1._dp
  zrsnow=1.e-3_dp
  zcdi=0.6_dp          !
  zmw0=4.19e-12_dp     !
  zmi0=1.e-12_dp       !
!--- End included for CDNC/IC scheme -----------------------------------
!
!     ------------------------------------------------------------------
!
!       1.   Top boundary conditions, air density and geopotential
!            height at half levels
!
!       1.1   Set to zero precipitation fluxes etc.
!
  DO 111 jl = 1,kproma
     zclcpre(jl)   = 0.0_dp
     zxiflux(jl)   = 0.0_dp
111 END DO
!
!       1.2   Air density
!
  DO 122 jk = ktdia,klev
     DO 121 jl = 1,kproma
        zrho(jl,jk)   = papm1(jl,jk)/(rd*ptvm1(jl,jk))
        pxtec(jl,jk)  = MAX(pxtec(jl,jk),0.0_dp)
        pqtec(jl,jk)  = MAX(pqtec(jl,jk),0.0_dp)
!--- Included for prognostic CDNC/IC scheme ----------------------------

        ! calculate cloud droplet number concentration 
        zcdnc(jl,jk)=MAX( ( pxtm1(jl,jk,idt_cdnc) + &
                            pxtte(jl,jk,idt_cdnc) * ztmst ) &
                            / M_air *1000._dp * zrho(jl,jk) , cqtmin)

!--- End included for CDNC/IC scheme -----------------------------------

121  END DO
122 END DO
!
!       1.3   Geopotential at half levels
!
  DO 132 jk = 2,klev
     DO 131 jl = 1,kproma
        zgeoh(jl,jk)   = 0.5_dp*(pgeo(jl,jk)+pgeo(jl,jk-1))
131  END DO
132 END DO
   DO 133 jl = 1,kproma
      zgeoh(jl,1)      = pgeo(jl,1)+(pgeo(jl,1)-zgeoh(jl,2))
      zgeoh(jl,klevp1) = 0.0_dp
133 END DO

!--- Included for prognostic CDNC/IC scheme ----------------------------
  !
  !       1.4  CLOUD TOPS
  !
  DO 230 jl=1,kproma
     IF (paclc(jl,1) .GE. zepsec) THEN
        itop(jl,1)=1
     ELSE
        itop(jl,1)=0
     ENDIF
230 END DO
  !
  DO 232 jk=ktdia+1,klev
     DO 231 jl=1,kproma
        IF (paclc(jl,jk) .GE. zepsec .AND. itop(jl,jk-1) .EQ. 0) THEN
           itop(jl,jk)=jk
        ELSE IF (paclc(jl,jk) .GE. zepsec) THEN
           itop(jl,jk)=itop(jl,jk-1)
        ELSE
           itop(jl,jk)=0
        ENDIF
231  END DO
232 END DO
  !
  !         1.5 CLOUD BASES
  !
  DO 233 jl=1,kproma
     IF (paclc(jl,klev) .GE. zepsec) THEN
        ibas(jl,klev)=klev
     ELSE
        ibas(jl,klev)=0
     ENDIF
233 END DO
!
  DO 235 jk=klev-1,ktdia,-1
     DO 234 jl=1,kproma
        IF (paclc(jl,jk) .GE. zepsec .AND. ibas(jl,jk+1) .EQ. 0) THEN
           ibas(jl,jk)=jk
        ELSE IF (paclc(jl,jk) .GE. zepsec) THEN
           ibas(jl,jk)=ibas(jl,jk+1)
        ELSE
           ibas(jl,jk)=0
        ENDIF
234  END DO
235 END DO

  DO jk=klev,ktdia,-1
     DO jl=1,kproma
        !
        it       = NINT(ptm1(jl,jk)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zesw     = tlucuaw(it)/papm1(jl,jk)
        zesw     = MIN(zesw,0.5_dp)
        zqsw     = zesw/(1._dp-vtmpc1*zesw)
        zsusatw_2d(jl,jk)=MAX(pqm1(jl,jk)/zqsw-1.0_dp,0.0_dp)
        !--- Store supersaturation with respect to water in stream:
        swat(jl,jk,jrow)=zsusatw_2d(jl,jk)
        !--- Saturation water vapour pressure:
        !zesw_2d(jl,jk)=zesw*papm1(jl,jk)*rv/rd

     END DO !jl
  END DO !jk

  DO jk=klev,ktdia,-1
     DO jl=1,kproma

        IF (ibas(jl,jk) .EQ. jk) THEN
           IF (zcdnc(jl,jk).LE.zcdnmin .OR. pxlte(jl,jk).GT.0._dp &
               .OR. paclc(jl,jk).GT.cloud_tm1(jl,jk,jrow) &
               .OR. zsusatw_2d(jl,jk).GT.zeps) THEN
             zqlnuc(jl,jk)=MAX(0._dp,pcdncact(jl,jk))
             zcdnc(jl,jk)=zcdnc(jl,jk)+zqlnuc(jl,jk)
             qnuc(jl,jk,jrow)=zqlnuc(jl,jk)/ztmst
           ENDIF
        ENDIF
        !
        IF (jk < klev) THEN
           IF (itop(jl,jk) .NE. ibas(jl,jk) .AND. ibas(jl,jk) .GT. jk &
                .AND. zqlnuc(jl,jk+1)>0._dp) THEN
              zqlnuc(jl,jk)=zqlnuc(jl,jk+1)
              qnuc(jl,jk,jrow)=(zcdnc(jl,jk+1)-zcdnc(jl,jk))/ztmst
              zcdnc(jl,jk)=zcdnc(jl,jk+1)
           ENDIF
        ENDIF

        IF (lodetr) THEN

           iclbas(jl)=NINT(pcvcbot(jl))

           zqlnuccvh(jl,jk)=0._dp 
           IF (jk.EQ.iclbas(jl)) THEN
              zvervc= (-pvervel(jl,jk)/(zrho(jl,jk)*g) + pwcape(jl) &
                   + 1.33_dp*SQRT(ptkem1(jl,jk)))
              zvervc=MAX(zvervc,1._dp)
              zcdncnew=zna(jl,jk)*zvervc/(zvervc+2.3E-10_dp*zna(jl,jk))
              zcdncnew=1.e5_dp*(1.0E-6_dp*zcdncnew)**1.27_dp
              zqlnuccvh(jl,jk)=MAX(0._dp,MIN(zna(jl,jk),zcdncnew))
           ENDIF

        END IF

     END DO
  END DO

  IF (lodetr) THEN

     DO jk=klev,ktdia,-1
        DO jl=1,kproma
           IF (pxtec(jl,jk).GT.0._dp.AND.iclbas(jl).GT.0.AND. &
                ptm1(jl,jk).GT.cthomi) THEN
              zqlnuccv=MAX(0._dp,MIN(zna(jl,jk),zqlnuccvh(jl,iclbas(jl))))
              qnuc(jl,jk,jrow)=zqlnuccv/ztmst
              zcdnc(jl,jk)=zcdnc(jl,jk)+zqlnuccv
           ENDIF
        END DO
     END DO

  END IF

  !corinna: set zcdnc to minium now if nucleation is not strong enough
  DO jk = ktdia, klev
     DO jl = 1, kproma
        IF (paclc(jl,jk).GT.0._dp .AND. ptm1(jl,jk).GT.cthomi) zcdnc(jl,jk)=MAX(zcdnc(jl,jk),zcdnmin)
     ENDDO !jl
  ENDDO ! jk

!--- End included for CDNC/IC scheme -----------------------------------

  DO 831 jk=ktdia,klev
!
!     ------------------------------------------------------------------
!
!       2.    Set to zero local tendencies (increments)
!
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

    DO 201 jl = 1,kproma
        zcnd(jl)       = 0.0_dp
        zdep(jl)       = 0.0_dp
        zfrl(jl)       = 0.0_dp
        zspr(jl)       = 0.0_dp
        zimlt(jl)      = 0.0_dp
        zximlt(jl)     = 0.0_dp
        zxisub(jl)     = 0.0_dp
        zsmlt(jl)      = 0.0_dp
        zsacl(jl)      = 0.0_dp
        zgenti(jl)     = 0.0_dp
        zgentl(jl)     = 0.0_dp
        zxievap(jl)    = 0.0_dp
        zxlevap(jl)    = 0.0_dp
        zvartg(jl)     = 0.0_dp
        zconvvar(jl)   = 0.0_dp
        zconvskew(jl)  = 0.0_dp
        zturbvar(jl)   = 0.0_dp
        zturbskew(jl)  = 0.0_dp
        zmicroskew(jl) = 0.0_dp
        zxvarte(jl)    = 0.0_dp 
        zxskewte(jl)   = 0.0_dp   
        zdp(jl)        = paphm1(jl,jk+1)-paphm1(jl,jk)
        zdz(jl)        = (zgeoh(jl,jk)-zgeoh(jl,jk+1))/g
        zrcp           = 1._dp/(cpd+zcons1*MAX(pqm1(jl,jk),0.0_dp))
        zlvdcp(jl)     = alv*zrcp
        zlsdcp(jl)     = als*zrcp
!--- Included for prognostic CDNC/IC scheme ----------------------------
!    additional arrays for transformation rate changes of CDNC and ICNC
        zrprn(jl)=0._dp
        zfrln(jl)=0._dp
        zsacln(jl)=0._dp
!--- End included for CDNC/IC scheme -----------------------------------
201  END DO
!
!     ------------------------------------------------------------------
!
!       3.   Modification of incoming precipitation fluxes by
!            melting, sublimation and evaporation
!
     IF (jk .GT. 1) THEN
!
!DIR$ CONCURRENT
        DO 331 jl = 1,kproma
!
!       3.1   Melting of snow and ice
!
           zcons     = zcons2*zdp(jl)/(zlsdcp(jl)-zlvdcp(jl))
           ztdif     = MAX(0.0_dp,ptm1(jl,jk)-tmelt)
           zsnmlt    = MIN(zxsec*zsfl(jl),zcons*ztdif)
           zrfl(jl)  = zrfl(jl)+zsnmlt
           zsfl(jl)  = zsfl(jl)-zsnmlt
           zsmlt(jl) = zsnmlt/(zcons2*zdp(jl))
           zximelt   = MIN(zxsec*zxiflux(jl),zcons*ztdif)
           pimelt(jl,jk) = zximelt ! mz_ht_20070611
           zxiflux(jl)=zxiflux(jl)-zximelt
           zximlt(jl) =zximelt/(zcons2*zdp(jl))
           IF (ztdif.GT.0.0_dp) THEN
            zimlt(jl) = MAX(0.0_dp,pxim1(jl,jk)+pxite(jl,jk)*ztmst)
!--- Included for prognostic CDNC/IC scheme (Philip Stier, 31/03/2004) -
!    If T > tmelt melt all ice crystals and transfer to cloud droplets
            IF (paclc(jl,jk) .GE. zepsec) THEN
               zicncmelt=(zimlt(jl)*zrho(jl,jk)/zmi)/paclc(jl,jk)
            ELSE
               zicncmelt=(zimlt(jl)*zrho(jl,jk)/zmi)
            END IF
            zcdnc(jl,jk)=zcdnc(jl,jk)+zicncmelt
            qmel(jl,jk,jrow)=zicncmelt/ztmst
!--- End included for CDNC/IC scheme -----------------------------------
           ELSE
            zimlt(jl) = 0.0_dp
           END IF
!
!       3.2   Sublimation of snow and ice (Lin et al., 1983)
!
           IF (zclcpre(jl) .GT. 0.0_dp) THEN
              zclcstar = zclcpre(jl)
              zdpg     = zdp(jl)/g
              zqrho    = 1.3_dp/zrho(jl,jk)
              it       = NINT(ptm1(jl,jk)*1000._dp)
              IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
              it = MAX(MIN(it,jptlucu2),jptlucu1)
              zesi     = tlucua(it)/papm1(jl,jk)
              zesi     = MIN(zesi,0.5_dp)
              zqsi     = zesi/(1._dp-vtmpc1*zesi)
              zsusati  = MIN(pqm1(jl,jk)/zqsi-1.0_dp,0.0_dp)
              zb1      = zlsdcp(jl)**2/(2.43e-2_dp*rv*(ptm1(jl,jk)**2))
              zb2      = 1._dp/(zrho(jl,jk)*zqsi*0.211e-4_dp)
              zcoeff   = 3.e6_dp*2._dp*api*zsusati/(zrho(jl,jk)*(zb1+zb2))
!
              IF (zsfl(jl) .GT. cqtmin) THEN
               zxsp1    = (zsfl(jl)/(zclcpre(jl)*cvtfall))**(1._dp/1.16_dp)
               zclambs  = (zxsp1/(api*crhosno*cn0s))**0.25_dp
               zcfac4c  = 0.78_dp*zclambs**2+232.19_dp*zqrho**0.25_dp           &
                                                    *zclambs**2.625_dp
               zzeps    = MAX(-zxsec*zsfl(jl)/zclcpre(jl),             &
                                                  zcoeff*zcfac4c*zdpg)
               zsub(jl) = -zzeps/zdpg*ztmst*zclcstar
               zsub(jl) = MIN(zsub(jl),                                &
                                    MAX(zxsec*(zqsi-pqm1(jl,jk)),0.0_dp))
               zsub(jl) = MAX(zsub(jl),0.0_dp)
              END IF
!
              IF (zxiflux(jl) .GT. cqtmin) THEN
               zxsp1    = (zxiflux(jl)/(zclcpre(jl)*cvtfall))**(1._dp/1.16_dp)
               zclambs  = (zxsp1/(api*crhosno*cn0s))**0.25_dp
               zcfac4c  = 0.78_dp*zclambs**2+232.19_dp*zqrho**0.25_dp           &
                                                    *zclambs**2.625_dp
               zzeps    = MAX(-zxsec*zxiflux(jl)/zclcpre(jl),          &
                                                  zcoeff*zcfac4c*zdpg)
               zsubi    = -zzeps/zdpg*ztmst*zclcstar
               zsubi    = MIN(zsubi,MAX(zxsec*(zqsi-pqm1(jl,jk)),0.0_dp))
               zsubi    = MAX(zsubi,0.0_dp)
               zxiflux(jl)=zxiflux(jl)-zsubi*zcons2*zdp(jl)
               zxisub(jl) =zsubi
              END IF
           END IF
!
!       3.3   Evaporation of rain (Rotstayn, 1997)
!
           IF (zclcpre(jl) .GT. 0.0_dp .AND. zrfl(jl) .GT. cqtmin) THEN
              zclcstar = zclcpre(jl)
              zdpg     = zdp(jl)/g
              zqrho    = 1.3_dp/zrho(jl,jk)
              zxrp1    = (zrfl(jl)/(zclcpre(jl)*12.45_dp*SQRT(zqrho)))    &
                                                       **(8._dp/9._dp)
              it       = NINT(ptm1(jl,jk)*1000._dp)
              IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
              it = MAX(MIN(it,jptlucu2),jptlucu1)
              zesw     = tlucuaw(it)/papm1(jl,jk)
              zesat    = zesw*papm1(jl,jk)*rv/rd
              zesw     = MIN(zesw,0.5_dp)
              zqsw     = zesw/(1._dp-vtmpc1*zesw)
              zsusatw  = MIN(pqm1(jl,jk)/zqsw-1.0_dp,0.0_dp)
              zdv      = 2.21_dp/papm1(jl,jk)
              zast     = alv*(alv/(rv*ptm1(jl,jk))-1.0_dp)/               &
                                            (0.024_dp*ptm1(jl,jk))
              zbst     = rv*ptm1(jl,jk)/(zdv*zesat)
              zzepr    = 870._dp*zsusatw*(zrfl(jl)/zclcpre(jl))**0.61_dp     &
                                     /(SQRT(zrho(jl,jk))*(zast+zbst))
              zzepr    = MAX(-zxsec*zrfl(jl)/zclcpre(jl),zzepr*zdpg)
              zevp(jl) = -zzepr/zdpg*ztmst*zclcstar
              zevp(jl) = MIN(zevp(jl),MAX(zxsec*(zqsw-pqm1(jl,jk)),0.0_dp))
              zevp(jl) = MAX(zevp(jl),0.0_dp)
           END IF
331     END DO
!
!        IF (lookupoverflow) CALL lookuperror ('cloud (1)    ')
!
         IF (lookupoverflow) THEN
           status_string = 'lookuperror: cdnc - cloud (1)'
           RETURN
         ENDIF
     END IF
!
!DIR$ CONCURRENT
     DO 610 jl=1,kproma
!
!     ------------------------------------------------------------------
!       4.    Sedimentation of cloud ice from grid-mean values.
!             Updating the tendency 'pxite' to include sedimentation.
!             At jk=klev, the sedimentation sink is balanced by
!             precipitation at the surface (through 'zzdrs', see 7.3).
!             Finally: In-cloud cloud water/ice.
!
        zxip1         = pxim1(jl,jk)+pxite(jl,jk)*ztmst-zimlt(jl)
        zxip1         = MAX(zxip1,EPSILON(1._dp))
        zxifall       = cvtfall*(zrho(jl,jk)*zxip1)**0.16_dp
        zal1          = zxifall*g*zrho(jl,jk)*ztmst/zdp(jl)
        zal2          = zxiflux(jl)/(zrho(jl,jk)*zxifall)
        zxised        = zxip1*EXP(-zal1)+zal2*(1._dp-EXP(-zal1))
        pisedi(jl,jk) = zxised! mz_ht_20071114
        zxiflux(jl)   = zxiflux(jl)+(zxip1-zxised)*zcons2*zdp(jl)
        pxite(jl,jk)  = (zxised-pxim1(jl,jk))/ztmst
!
!             In-cloud water/ice calculated from respective grid-means,
!             partial cloud cover, advective/diffusive tendencies,
!             detrained cloud water/ice and ice sedimentation.
!             In-cloud values are required for cloud microphysics.
!
        zclcaux(jl) = paclc(jl,jk)
        locc        = zclcaux(jl) .GT. 0.0_dp
        lo2         = (ptm1(jl,jk) .LT. cthomi) .OR.                   &
                      (ptm1(jl,jk) .LT. tmelt .AND. zxised .GT. csecfrl &
                      .AND. zsusatw_2d(jl,jk) .LT. 1._dp)
        IF (lo2) THEN                 !     ice cloud
           zxite(jl)          = pxtec(jl,jk)
           zxlte(jl)          = 0.0_dp
           IF (locc) THEN
              zxib(jl)        = pxim1(jl,jk)/zclcaux(jl)
              zxlb(jl)        = pxlm1(jl,jk)/zclcaux(jl)
              zxim1evp        = 0.0_dp
              zxlm1evp        = 0.0_dp
              zxidt           = (pxite(jl,jk)+zxite(jl))*ztmst
              zxldt           =  pxlte(jl,jk)*ztmst+zximlt(jl)+zimlt(jl)
              IF (zxidt .GT. 0.0_dp) THEN
                 zxidtstar    = zxidt
                 zxib(jl)     = zxib(jl)+zxidt
              ELSE
                 zxidtstar    = 0.0_dp
                 zxib(jl)     = zxib(jl)+MAX(zxidt/zclcaux(jl),        &
                                                       -zxib(jl))
                 pxite(jl,jk) = MAX(pxite(jl,jk),                      &
                                 -(pxim1(jl,jk)/ztmst+zxite(jl)))
              END IF
              IF (zxldt .GT. 0.0_dp) THEN
                 zxldtstar    = zxldt
                 zxlb(jl)     = zxlb(jl)+zxldt
              ELSE
                 zxldtstar    = 0.0_dp
                 zxlb(jl)     = zxlb(jl)+MAX(zxldt/zclcaux(jl),        &
                                                       -zxlb(jl))
                 pxlte(jl,jk) = MAX(pxlte(jl,jk),-pxlm1(jl,jk)/ztmst)
              END IF
           ELSE                      !    cloud cover = 0.
              zxib(jl)        = 0.0_dp
              zxlb(jl)        = 0.0_dp
              zxidt           = (pxite(jl,jk)+zxite(jl))*ztmst
              zxldt           =  pxlte(jl,jk)*ztmst+zximlt(jl)+zimlt(jl)
              IF (zxidt .GT. 0.0_dp) THEN
                 zxidtstar    = zxidt
                 zxim1evp     = pxim1(jl,jk)
              ELSE
                 zxidtstar    = 0.0_dp
                 pxite(jl,jk) = MAX(pxite(jl,jk),                      &
                                  -(pxim1(jl,jk)/ztmst+zxite(jl)))
                 zxim1evp     = pxim1(jl,jk)+(pxite(jl,jk)+zxite(jl))  &
                                                                 *ztmst
              END IF
              IF (zxldt .GT. 0.0_dp) THEN
                 zxldtstar    = zxldt
                 zxlm1evp     = pxlm1(jl,jk)
              ELSE
                 zxldtstar    = 0.0_dp
                 pxlte(jl,jk) = MAX(pxlte(jl,jk),-pxlm1(jl,jk)/ztmst)
                 zxlm1evp     = pxlm1(jl,jk)+pxlte(jl,jk)*ztmst
              END IF
           END IF
        ELSE                           !    water cloud
           zxlte(jl)          = pxtec(jl,jk)
           zxite(jl)          = 0.0_dp
           IF (locc) THEN
              zxlb(jl)        = pxlm1(jl,jk)/zclcaux(jl)
              zxib(jl)        = pxim1(jl,jk)/zclcaux(jl)
              zxlm1evp        = 0.0_dp
              zxim1evp        = 0.0_dp
              zxldt           = (pxlte(jl,jk)+zxlte(jl))*ztmst         &
                                             +zximlt(jl)+zimlt(jl)
              zxidt           =  pxite(jl,jk)*ztmst
              IF (zxldt .GT. 0.0_dp) THEN
                 zxldtstar    = zxldt
                 zxlb(jl)     = zxlb(jl)+zxldt
              ELSE
                 zxldtstar    = 0.0_dp
                 zxlb(jl)     = zxlb(jl)+MAX(zxldt/zclcaux(jl),        &
                                                          -zxlb(jl))
                 pxlte(jl,jk) = MAX(pxlte(jl,jk),                      &
                                 -(pxlm1(jl,jk)/ztmst+zxlte(jl)))
              END IF
              IF (zxidt .GT. 0.0_dp) THEN
                 zxidtstar    = zxidt
                 zxib(jl)     = zxib(jl)+zxidt
              ELSE
                 zxidtstar    = 0.0_dp
                 zxib(jl)     = zxib(jl)+MAX(zxidt/zclcaux(jl),        &
                                                          -zxib(jl))
                 pxite(jl,jk) = MAX(pxite(jl,jk),-pxim1(jl,jk)/ztmst)
              END IF
           ELSE                          !    cloud cover = 0.
              zxlb(jl)        = 0.0_dp
              zxib(jl)        = 0.0_dp
              zxldt           = (pxlte(jl,jk)+zxlte(jl))*ztmst         &
                                             +zximlt(jl)+zimlt(jl)
              zxidt           =  pxite(jl,jk)*ztmst
              IF (zxldt .GT. 0.0_dp) THEN
                 zxldtstar    = zxldt
                 zxlm1evp     = pxlm1(jl,jk)
              ELSE
                 zxldtstar    = 0.0_dp
                 pxlte(jl,jk) = MAX(pxlte(jl,jk),                      &
                                  -(pxlm1(jl,jk)/ztmst+zxlte(jl)))
                 zxlm1evp     = pxlm1(jl,jk)+(pxlte(jl,jk)+zxlte(jl))  &
                                                                *ztmst
              END IF
              IF (zxidt .GT. 0.0_dp) THEN
                 zxidtstar    = zxidt
                 zxim1evp     = pxim1(jl,jk)
              ELSE
                 zxidtstar    = 0.0_dp
                 pxite(jl,jk) = MAX(pxite(jl,jk),-pxim1(jl,jk)/ztmst)
                 zxim1evp     = pxim1(jl,jk)+pxite(jl,jk)*ztmst
              END IF
           END IF
        END IF
!
!     ------------------------------------------------------------------
!       5.    Condensation/deposition and evaporation/sublimation
!
!             zlc       =  L_{v/s} / c_p
!             zlcdqsdt  = L dq_sat / c_p dT
!             zdqsdt    = dq_sat / dT
!
        zrcp        = 1._dp/(cpd+zcons1*MAX(pqm1(jl,jk),0.0_dp))
        zlvdcp(jl)  = alv*zrcp
        zlsdcp(jl)  = als*zrcp
        zlc         = MERGE(zlsdcp(jl),zlvdcp(jl),lo2)
        it          = NINT(ptm1(jl,jk)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zqsm1       = MERGE(tlucua(it),tlucuaw(it),lo2)/papm1(jl,jk)
        zqsm1       = MIN(zqsm1,0.5_dp)
        zqsm1       = zqsm1/(1._dp-vtmpc1*zqsm1)
        it1         = it+1
        it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
        zqst1       = MERGE(tlucua(it1),tlucuaw(it1),lo2)/papm1(jl,jk)
        zqst1       = MIN(zqst1,0.5_dp)
        zqst1       = zqst1/(1._dp-vtmpc1*zqst1)
        zdqsdt      = (zqst1-zqsm1)*1000._dp
        zlcdqsdt    = zlc*zdqsdt
        zxievap(jl) = (1.0_dp-zclcaux(jl))*zxidtstar+zxim1evp
        zxlevap(jl) = (1.0_dp-zclcaux(jl))*zxldtstar+zxlm1evp
        zqvdt       = pqte(jl,jk)*ztmst+zevp(jl)+zsub(jl)              &
                             +zxievap(jl)+zxlevap(jl)+zxisub(jl)
        zdtdt       = ptte(jl,jk)*ztmst-zlvdcp(jl)*(zevp(jl)           &
                       +zxlevap(jl)) -(zlsdcp(jl)-zlvdcp(jl))          &
                       *(zsmlt(jl)+zximlt(jl)+zimlt(jl))       
        zqp1        = MAX(pqm1(jl,jk)+zqvdt,0.0_dp)
        ztp1        = ptm1(jl,jk)+zdtdt
        zdtdtstar   = zdtdt+zclcaux(jl)*(zlc*pqte(jl,jk)*ztmst         &
                           +zlvdcp(jl)*(zevp(jl)+zxlevap(jl))          &
                         +zlsdcp(jl)*(zsub(jl)+zxievap(jl)+zxisub(jl)))
        zdqsat      = zdtdtstar*zdqsdt/(1._dp+zclcaux(jl)*zlcdqsdt)
        zxib(jl)    = MAX(zxib(jl),0.0_dp)
        zxlb(jl)    = MAX(zxlb(jl),0.0_dp)
        zxilb       = zxib(jl)+zxlb(jl)
!
!       Diagnostics: relative humidity
!
        zrelhum        = pqm1(jl,jk)/zqsm1
        zrelhum        = MAX(MIN(zrelhum,1._dp),0._dp)
        prelhum(jl,jk) = zrelhum

        IF (jk.GE.ncctop) THEN
!
!       define variables needed for cover scheme
!
!       zbetaqt = total water
!       zbetass = saturation mixing ratio adjusted to match qv
!       zwide   = current diagnosed distribution width
!
           zbetacl(jl) = MAX(0.0_dp,pxlm1(jl,jk))+MAX(0.0_dp,pxim1(jl,jk))
           zbetaqt(jl) = MAX(cqtmin,pqm1(jl,jk))+zbetacl(jl)
           zvartg(jl)  = MAX(cqtmin,cvarmin*pqm1(jl,jk))
           zwide(jl)   = MAX(zvartg(jl),pbetab(jl,jk)-pbetaa(jl,jk))
           zskew       = MAX(MIN(pxskew(jl,jk),cbeta_pq_max),cbeta_pq)
           iqidx       = INT((nbetaq/cbetaqs)*                         &
                                 LOG((zskew-cbeta_pq)/rbetak+1._dp)+0.5_dp)
!
!
!       5.1 Turbulence: Skewness - equation solved implicitly
!           This solver only works if phmixtau has non-zero timescale
!
           zqtau         = phmixtau(jl,jk)+pvmixtau(jl,jk)
           zbqp1         = cbeta_pq-(cbeta_pq-pxskew(jl,jk))           &
                                                    *EXP(-zqtau*zdtime)
           zbqp1         = MAX(MIN(zbqp1,cbeta_pq_max),cbeta_pq)
           zturbskew(jl) = (zbqp1-pxskew(jl,jk))/zdtime
!
!       5.2 Turbulence: variance - equation solved implicitly
!
           zpp          = cbeta_pq
           zqq          = pxskew(jl,jk)
           zeta         = (zpp+zqq)**2*(zpp+zqq+1._dp)/(zpp*zqq)
           zprod        = zeta*pvdiffp(jl,jk)/zwide(jl)
           zbbap1       = zprod/zqtau+zvartg(jl)-(zprod/zqtau          &
                            +zvartg(jl)-zwide(jl))*EXP(-zqtau*zdtime)
           zbbap1       = MAX(zbbap1,zvartg(jl))
           zbbap1       = MIN(zbbap1,zbetaqt(jl)*(cbeta_pq+zbqp1)      &
                                                            /cbeta_pq)
           zturbvar(jl) = (zbbap1-zwide(jl))/zdtime
           zbap1        = zbetaqt(jl)-zbbap1*cbeta_pq/(cbeta_pq+zbqp1)
!
           IF (lcover) THEN
!              translated into apparent xl,xi,q and heat sources
!              first order effect only, effect of evaporation of
!              cloud on qsat taken into account in thermodynamic budget
!              but does not change the mixing term here since that
!              would require iteration and is therefore neglected
!
!              calculate values after one timestep
!
           iqidx   = INT((nbetaq/cbetaqs)*                             &
                     LOG((zbqp1-cbeta_pq)/rbetak+1._dp)+0.5_dp)
           ztt     = (pbetass(jl,jk)-zbap1)*cbeta_pq/                     &
                     ((zbetaqt(jl)-zbap1)*(cbeta_pq+zbqp1))
           ztt     = nbetax*MAX(MIN(ztt,1.0_dp),0.0_dp)
           ixidx   = INT(ztt)
           IF (ixidx == nbetax) THEN
              zbetai0 = 1.0_dp
              zbetai1 = 1.0_dp
           ELSE
              zbetai0 = (ztt-ixidx)*tbetai0(iqidx,ixidx+1)             &
                             +(ixidx+1._dp-ztt)*tbetai0(iqidx,ixidx)
              zbetai1 = (ztt-ixidx)*tbetai1(iqidx,ixidx+1)             &
                             +(ixidx+1._dp-ztt)*tbetai1(iqidx,ixidx)
           ENDIF
           zqp1b      = (zbetaqt(jl)-zbap1)*zbetai1 -                  &
                        (pbetass(jl,jk)-zbap1)*zbetai0 + pbetass(jl,jk)
           zgent      = MAX(pqm1(jl,jk)-zqp1b,-zxilb*zclcaux(jl))
           zgent      = MIN(zgent,zqsec*zqp1)              ! limit to qv
           zifrac     = zxib(jl)/MAX(zepsec,zxilb)
           zifrac     = MAX(MIN(zifrac,1.0_dp),0.0_dp)
           zgenti(jl) = zgent*zifrac
           zgentl(jl) = zgent*(1.0_dp-zifrac)
           IF (locc) THEN
              zxib(jl) = MAX(zxib(jl)+zgenti(jl)/zclcaux(jl),0.0_dp)
              zxlb(jl) = MAX(zxlb(jl)+zgentl(jl)/zclcaux(jl),0.0_dp)
           END IF
           zxilb       = zxib(jl)+zxlb(jl)
!
!       5.3 Deposition/sublimation of cloud ice and condensation/
!           evaporation of liquid water due to changes in water vapour
!           and temperature (advection, convection, turbulent mixing,
!           evaporation of rain, sublimation and melting of snow).
!           Translate PDF laterally to calculate cloud
!           after one timestep
!
           zqvdt     = zqvdt-zgent
           zdtdt     = zdtdt+zlvdcp(jl)*zgentl(jl)+zlsdcp(jl)*zgenti(jl)
           zqp1      = MAX(pqm1(jl,jk)+zqvdt,0.0_dp)
           ztp1      = ptm1(jl,jk)+zdtdt
           zdtdtstar = zdtdt+zclcaux(jl)*(zlc*pqte(jl,jk)*ztmst        &
              +zlvdcp(jl)*(zevp(jl)+zxlevap(jl)-zgentl(jl))            &
              +zlsdcp(jl)*(zsub(jl)+zxievap(jl)+zxisub(jl)-zgenti(jl)))
           zdqsat    = zdtdtstar*zdqsdt/(1._dp+zclcaux(jl)*zlcdqsdt)
           ztt       = (pbetass(jl,jk)-zqvdt+zdqsat-zbap1)/zbbap1
           ztt       = nbetax*MAX(MIN(ztt,1.0_dp),0.0_dp)
           ixidx     = INT(ztt)
           IF (ixidx == nbetax) THEN
              zbetai0 = 1.0_dp
              zbetai1 = 1.0_dp
           ELSE
              zbetai0 = (ztt-ixidx)*tbetai0(iqidx,ixidx+1)             &
                       +(ixidx+1._dp-ztt)*tbetai0(iqidx,ixidx)
              zbetai1 = (ztt-ixidx)*tbetai1(iqidx,ixidx+1)             &
                       +(ixidx+1._dp-ztt)*tbetai1(iqidx,ixidx)
           ENDIF
           zaa        = pbetaa(jl,jk)
           zqcdif     = (zbetaqt(jl)-zaa)*(1._dp-zbetai1)               &
                         +(zaa+zqvdt-pbetass(jl,jk)-zdqsat)*(1._dp-zbetai0)
           zqcdif     = MAX(0.0_dp,zqcdif)-zbetacl(jl)
           zqcdif     = MAX(zqcdif,-zxilb*zclcaux(jl))
           zqcdif     = MIN(zqcdif,zqsec*zqp1)             ! limit to qv
!
           IF (zqcdif .LT. 0.0_dp) THEN                 ! cloud dissipation
              zifrac   = zxib(jl)/MAX(zepsec,zxilb)
              zifrac   = MAX(MIN(zifrac,1.0_dp),0.0_dp)
              zdep(jl) = zqcdif*zifrac
              zcnd(jl) = zqcdif*(1.0_dp-zifrac)
           ELSE                                      ! cloud generation

              IF (lo2) THEN                          ! deposition
                 zdep(jl) = zqcdif
                 zcnd(jl) = 0.0_dp
              ELSE                                   ! condensation
!--- Included/changed for prognostic CDNC/IC scheme --------------------
!    Use standard condensation for empirical Lin & Leaitch approach and 
!    explicit condensation after Levkov et al. 1992 for the explicit
!    activation schemes that allow for supersaturation:
                 IF (ncdnc==1) THEN
                    zcnd(jl) = zqcdif
                 ELSE IF (ncdnc>1) THEN
                    zdv      = 2.21_dp/papm1(jl,jk)
                    it       = NINT(ptm1(jl,jk)*1000._dp)
                    IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
                    it = MAX(MIN(it,jptlucu2),jptlucu1)
                    zesw     = tlucuaw(it)/papm1(jl,jk)
                    zesat    = zesw*papm1(jl,jk)*rv/rd
                    zast     = alv*(alv/(rv*ptm1(jl,jk))-1.0_dp)/               &
                               (0.024_dp*ptm1(jl,jk))
                    zbst     = rv*ptm1(jl,jk)/(zdv*zesat)
                    zgtp=1._dp/(zrho(jl,jk)*(zast+zbst))
                    zcond=0.5_dp*api*zdw0*swat(jl,jk,jrow)*zcdnc(jl,jk)  &
                          *zgtp*ztmst*zclcaux(jl)
                    zcnd(jl)= MIN(zcond,zqcdif)
                 END IF
!--- End included for CDNC/IC scheme -----------------------------------
                 zdep(jl) = 0.0_dp
              END IF
           END IF
          END IF !lcover
        END IF !ncctop
!
        IF ((.NOT.lcover) .OR. jk < ncctop) THEN
           zqcdif         = (zqvdt-zdqsat)*zclcaux(jl)
           zqcdif         = MAX(zqcdif,-zxilb*zclcaux(jl))
           zqcdif         = MIN(zqcdif,zqsec*zqp1)
           IF (zqcdif .LT. 0.0_dp) THEN                 ! cloud dissipation
              zifrac      = zxib(jl)/MAX(zepsec,zxilb)
              zifrac      = MAX(MIN(zifrac,1.0_dp),0.0_dp)
              zdep(jl)    = zqcdif*zifrac
              zcnd(jl)    = zqcdif*(1.0_dp-zifrac)
           ELSE                                      ! cloud generation
              IF (lo2) THEN                          ! deposition
                 zdep(jl) = zqcdif
                 zcnd(jl) = 0.0_dp
              ELSE                                   ! condensation
!--- Included/changed for prognostic CDNC/IC scheme --------------------
!    Use standard condensation for empirical Lin & Leaitch approach and 
!    explicit condensation after Levkov et al. 1992 for the explicit
!    activation schemes that allow for supersaturation:
                 IF (ncdnc==1) THEN
                    zcnd(jl) = zqcdif
                 ELSE IF (ncdnc>1) THEN
                    zdv      = 2.21_dp/papm1(jl,jk)
                    it       = NINT(ptm1(jl,jk)*1000._dp)
                    IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
                    it = MAX(MIN(it,jptlucu2),jptlucu1)
                    zesw     = tlucuaw(it)/papm1(jl,jk)
                    zesat    = zesw*papm1(jl,jk)*rv/rd
                    zast     = alv*(alv/(rv*ptm1(jl,jk))-1.0_dp)/               &
                               (0.024_dp*ptm1(jl,jk))
                    zbst     = rv*ptm1(jl,jk)/(zdv*zesat)
                    zgtp=1._dp/(zrho(jl,jk)*(zast+zbst))
                    zcond=0.5_dp*api*zdw0*swat(jl,jk,jrow)*zcdnc(jl,jk)  &
                          *zgtp*ztmst*zclcaux(jl)
                    zcnd(jl)= MIN(zcond,zqcdif)
                 END IF
!--- End included for CDNC/IC scheme -----------------------------------
                 zdep(jl) = 0.0_dp
              END IF
           END IF
        END IF !lcover
!
!       5.4 Accounting for cloud evaporation in clear air and
!           checking for supersaturation
!
        ztp1tmp(jl) = ztp1+zlvdcp(jl)*zcnd(jl)+zlsdcp(jl)*zdep(jl)
        zqp1tmp(jl) = zqp1-zcnd(jl)-zdep(jl)
        zxip1       = MAX(zxised+zxite(jl)*ztmst-zxievap(jl)           &
                               +zgenti(jl)+zdep(jl),0.0_dp)   
        lo2         = (ztp1tmp(jl) .LT. cthomi) .OR.                   &
                      (ztp1tmp(jl) .LT. tmelt .AND. zxip1 .GT. csecfrl &
                      .AND. zsusatw_2d(jl,jk) .LT. 1._dp)
        it          = NINT(ztp1tmp(jl)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zes         = MERGE(tlucua(it),tlucuaw(it),lo2)/papp1(jl,jk)
        zes         = MIN(zes,0.5_dp)
        LO          = zes<0.4_dp
        zcor        = 1._dp/(1._dp-vtmpc1*zes)
        zqsp1tmp    = zes*zcor
        zoversat    = zqsp1tmp*0.01_dp
        zrhtest     = MIN(pqm1(jl,jk)/zqsm1,1._dp)*zqsp1tmp
        it1         = it+1
        it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
        zqst1       = MERGE(tlucua(it1),tlucuaw(it1),lo2)/papp1(jl,jk)
        zqst1       = MIN(zqst1,0.5_dp)
        zqst1       = zqst1/(1._dp-vtmpc1*zqst1)
        zdqsdt      = (zqst1-zqsp1tmp)*1000._dp
        zlc         = MERGE(zlsdcp(jl),zlvdcp(jl),lo2)
        zlcdqsdt    = MERGE(zlc*zdqsdt,zqsp1tmp*zcor*tlucub(it),LO)
        zqcon       = 1._dp/(1._dp+zlcdqsdt)
!
        IF (lo2) THEN                                    ! ice cloud
           IF (zqp1tmp(jl) .GT. zqsp1tmp+zoversat) THEN
              zdepcor     = (zqp1tmp(jl)-zqsp1tmp-zoversat)*zqcon
              zdep(jl)    = zdep(jl)+zdepcor
           END IF
           IF (zdep(jl) .GT. 0.0_dp .AND. zqp1tmp(jl) .LT. zrhtest        &
                                 .AND. zqsp1tmp .LE. zqsm1) THEN
              zdep(jl)    = zqp1-zrhtest
              zdep(jl)    = MAX(zdep(jl),0.0_dp)
           END IF
        ELSE                                             ! water cloud
           IF (zqp1tmp(jl) .GT. zqsp1tmp+zoversat) THEN
              zcndcor     = (zqp1tmp(jl)-zqsp1tmp-zoversat)*zqcon
              zcnd(jl)    = zcnd(jl)+zcndcor
           END IF
           IF (zcnd(jl) .GT. 0.0_dp .AND. zqp1tmp(jl) .LT. zrhtest        &
                                 .AND. zqsp1tmp .LE. zqsm1) THEN
              zcnd(jl)    = zqp1-zrhtest
              zcnd(jl)    = MAX(zcnd(jl),0.0_dp)
           END IF
        END IF
!
!       5.5 Change of in-cloud water due to deposition/sublimation and
!           condensation/evaporation (input for cloud microphysics)
!
        zrelhum=zqp1tmp(jl)/zqsp1tmp
        zdepos =MAX(zdep(jl)+zgenti(jl),0.0_dp)
        zcond  =MAX(zcnd(jl)+zgentl(jl),0.0_dp)
        pcond(jl,jk) = zcond
        IF (locc) THEN
           zxib(jl) = MAX(zxib(jl)+zdep(jl)/zclcaux(jl),0.0_dp)
           zxlb(jl) = MAX(zxlb(jl)+zcnd(jl)/zclcaux(jl),0.0_dp)
        ELSEIF (zdepos>0.0_dp .OR. zcond>0.0_dp) THEN
           zclcaux(jl)=MAX(MIN(zrelhum,1.0_dp),0.01_dp)
           zxib(jl)   = zdepos/zclcaux(jl)
           zxlb(jl)   = zcond /zclcaux(jl)
        END IF
!--- Included for prognostic CDNC/IC scheme ----------------------------
        IF (zclcaux(jl).GT.0._dp .AND. zxlb(jl).GT. cqtmin) THEN
           IF (zqlnuc(jl,jk) .LE. 0._dp) THEN    !if there was no previous nucleation
              zqlnuc(jl,jk)=MAX(0._dp,pcdncact(jl,jk))
              zcdnc(jl,jk)=zcdnc(jl,jk)+zqlnuc(jl,jk)
              qnuc(jl,jk,jrow)=zqlnuc(jl,jk)/ztmst
           ENDIF
           zcdnc(jl,jk)=MAX(zcdnc(jl,jk),zcdnmin)
        ELSE
           zcdnc(jl,jk)=cqtmin
        ENDIF
!--- End included for CDNC/IC scheme -----------------------------------
        ztp1tmp(jl) = ztp1+zlvdcp(jl)*zcnd(jl)+zlsdcp(jl)*zdep(jl)
!
!     ------------------------------------------------------------------
!       6.    Freezing of cloud water
!
!       6.1   Freezing of cloud water for T < 238 K
!
        IF (ztp1tmp(jl) .LE. cthomi) THEN
           zfrl(jl)  = zxlb(jl)*zclcaux(jl)
           zxib(jl)  = zxib(jl)+zxlb(jl)
           zxlb(jl)  = 0.0_dp
!--- Included for prognostic CDNC/IC scheme ----------------------------
           znfrl=zcdnc(jl,jk)
           zcdnc(jl,jk)=cqtmin
           qfre(jl,jk,jrow)=-znfrl/ztmst
!--- End included for CDNC/IC scheme -----------------------------------
        END IF
!
!       6.2   Freezing of cloud water between 238 and 273 K
!
!---Changed for prognostic CDNC/IC scheme ------------------------------

        lo           = zxlb(jl) .GT. cqtmin .AND. ztp1tmp(jl) .LT. tmelt  &
                       .AND. ztp1tmp(jl) .GT. cthomi .AND. locc
        IF (lo) THEN
!   (Replaced pacdnc by zcdnc which is set above)
         IF (zcdnc(jl,jk) .GE. zcdnmin) THEN
           zfrl(jl)  = 100._dp*(EXP(0.66_dp*(tmelt-ztp1tmp(jl)))-1._dp)         &
                                  *zrho(jl,jk)/(rhoh2o*zcdnc(jl,jk))
           zfrl(jl)  = zxlb(jl)*(1._dp-1._dp/(1._dp+zfrl(jl)*ztmst*zxlb(jl)))
           zradl     = (0.75_dp*zxlb(jl)*zrho(jl,jk)                      &
                               /(api*rhoh2o*zcdnc(jl,jk)))**(1._dp/3._dp)
           zf1       = 4._dp*api*zradl*zcdnc(jl,jk)*2.e5_dp                  &
                                 *(tmelt-3._dp-ztp1tmp(jl))/zrho(jl,jk)
!--- End changed for prognostic CDNC/IC scheme -------------------------
           zf1       = MAX(0.0_dp,zf1)
!---Changed for prognostic CDNC/IC scheme ------------------------------
           zfrl(jl)  = zfrl(jl) + zxlb(jl)*(1._dp-                        &
                exp(-ztmst*1.4e-20_dp*zf1/zxlb(jl)))
!--- End changed for prognostic CDNC/IC scheme -------------------------
           zfrl(jl)  = MAX(0.0_dp,MIN(zfrl(jl),zxlb(jl)))
!--- Included for prognostic CDNC/IC scheme ----------------------------
           !
           ! freezing of cloud droplets
           !
           IF (zxlb(jl)-zfrl(jl) .GT. cqtmin) THEN
              zfrln(jl)=MIN(zcdnc(jl,jk)*zfrl(jl)/(zxlb(jl)+zeps),zcdnc(jl,jk)-zcdnmin)
           ELSE
              zfrln(jl)=MAX(0.0_dp,MIN(zcdnc(jl,jk)*zfrl(jl)/(zxlb(jl)+zeps),zcdnc(jl,jk)-cqtmin))
           END IF
           zcdnc(jl,jk)=MAX(zcdnc(jl,jk)-zfrln(jl),cqtmin)
!--- End included for CDNC/IC scheme -----------------------------------
           zxlb(jl)  = zxlb(jl)-zfrl(jl)
           zxib(jl)  = zxib(jl)+zfrl(jl)
           zfrl(jl)  = zfrl(jl)*zclcaux(jl)
!-------corinna: Bergeron-Findeisen-Process
           IF (locc .AND. zdep(jl)>0._dp .AND. zxlb(jl)>0._dp .AND. zxib(jl)>csecfrl &
                .AND. zsusatw_2d(jl,jk) .LT. 1._dp) THEN
              zzevp=zxlb(jl)/ztmst*zclcaux(jl)
              pxlte(jl,jk)=pxlte(jl,jk)-zzevp
              pqte(jl,jk) =pqte(jl,jk)+zzevp
              ptte(jl,jk) =ptte(jl,jk)-zlvdcp(jl)*zzevp
              zcdnc(jl,jk)=cqtmin
              zxlb(jl)=0._dp
           END IF
!-------end corinna: Bergeron-Findeisen-Process
        END IF
!--- Included for prognostic CDNC/IC scheme ----------------------------
     ENDIF
!--- End included for CDNC/IC scheme ----------------------------------- 
!
610  END DO
!
!     IF (lookupoverflow) CALL lookuperror ('cloud (2)    ')
!
     IF (lookupoverflow) THEN
       status_string = 'lookuperror: cdnc - cloud (2)'
       RETURN
     ENDIF
!     ------------------------------------------------------------------
!       7.  Cloud physics and precipitation fluxes at the surface
!
!ham_ps: cdir circumvents bug in sxf90 compiler
!CDIR NOMOVEDIV
     DO 701 jl = 1,kproma
        locc     = zclcaux(jl) .GT. 0.0_dp
        zclcstar = MIN(zclcaux(jl),zclcpre(jl))
        zauloc   = cauloc*zdz(jl)/5000._dp
        zauloc   = MAX(MIN(zauloc,clmax),clmin)
!
        jb=knvb(jl)
        lo=(jb.GE.jbmin .AND. jb.LE.jbmax .AND. pvervel(jl,jk).GT.0._dp)
        lo1=(jk.EQ.jb .OR. jk.EQ.jb+1)
        IF(lo .AND. lo1 .AND. lonacc) zauloc= 0.0_dp
!
        zqrho    = 1.3_dp/zrho(jl,jk)
        zxlb(jl) = MAX(zxlb(jl),1.e-20_dp)
        zxib(jl) = MAX(zxib(jl),1.e-20_dp)

! liquid water and snow content are stored 
! before the reduction by outfalling rain
! (necessary for nucleation scavenging)
        plwc(jl,jk) = zxlb(jl)
        piwc(jl,jk) = zxib(jl)      

        IF (zclcpre(jl) .GT. 0.0_dp .AND. zrfl(jl) .GT. cqtmin) THEN
           zxrp1 = (zrfl(jl)/(zclcpre(jl)*12.45_dp*SQRT(zqrho)))**(8._dp/9._dp)
         ELSE    
           zxrp1 = 0.0_dp
         END IF

        IF (zclcpre(jl) .GT. 0.0_dp .AND. zsfl(jl) .GT. cqtmin) THEN
           zxsp1 = (zsfl(jl)/(zclcpre(jl)*cvtfall))**(1._dp/1.16_dp)
        ELSE
           zxsp1 = 0.0_dp
        END IF
!
!       7.1   Warm clouds: Coalescence processes after Beheng (1994):
!             Autoconversion of cloud droplets and collection of cloud
!             droplets by falling rain. Accretion of cloud droplets by
!             falling snow (zsacl) is calculated under 7.2
!
!---Changed for prognostic CDNC/IC scheme ------------------------------
        IF (locc) THEN
           IF (zxlb(jl) > cqtmin .AND. zcdnc(jl,jk) .GE. zcdnmin) THEN
!--- End changed for prognostic CDNC/IC scheme -------------------------
!--- Included alternative autoconversion parameterisation --------------

           IF (nauto==2) THEN

!           Autoconversion rate from Khairoutdinov and Kogan, 2000

              zraut    = ccraut*1350._dp*(zcdnc(jl,jk)*1.e-6_dp)**(-1.79_dp)
              zexm1    = 2.47_dp-1.0_dp
              zexp     = -1._dp/zexm1
              zraut    = zxlb(jl)*(1._dp-(1._dp+zraut*ztmst*zexm1*zxlb(jl)      &
                                                       **zexm1)**zexp)
              zraut    = MIN(zxlb(jl),zraut)
              zxlb(jl) = zxlb(jl)-zraut
              zrac1    = 6._dp*zxrp1*ztmst
              zrac1    = zxlb(jl)*(1._dp-EXP(-zrac1))
              zxlb(jl) = zxlb(jl)-zrac1
              zrac2    = 6._dp*zauloc*zrho(jl,jk)*zraut*ztmst
              zrac2    = zxlb(jl)*(1._dp-EXP(-zrac2))
              zxlb(jl) = zxlb(jl)-zrac2
              zrpr(jl) = zrpr(jl)+zclcaux(jl)*(zraut+zrac2)+zclcstar*zrac1

!--- Included for prognostic CDNC/IC scheme ----------------------------
              zxlbold=zxlb(jl)+zrac1+zrac2+zraut
              zrprn(jl)=(zraut+zrac1+zrac2)/(zxlbold+zeps)
              IF (zxlb(jl) .GT. cqtmin) THEN
                 zrprn(jl)=MIN(zcdnc(jl,jk)*zrprn(jl),zcdnc(jl,jk)-zcdnmin)
              ELSE
                 zrprn(jl)=MIN(zcdnc(jl,jk)*zrprn(jl),zcdnc(jl,jk))
              END IF
              zcdnc(jl,jk)=MAX(zcdnc(jl,jk)-zrprn(jl),cqtmin)
!--- End included for CDNC/IC scheme ------------------------------------
!--- End included alternative autoconversion parameterisation ----------
!--- Changed for alternative autoconversion parameterisation -----------

           ELSE IF (nauto==1) THEN

!---Changed for prognostic CDNC/IC scheme ------------------------------
!   (Replaced pacdnc by zcdnc which is set above)
              zraut    = ccraut*1.2e27_dp/zrho(jl,jk)*(zcdnc(jl,jk)*1.e-6_dp)  &
                         **(-3.3_dp)*(zrho(jl,jk)*1.e-3_dp)**4.7_dp
!--- End changed for prognostic CDNC/IC scheme -------------------------

              zexm1    = 4.7_dp-1.0_dp
              zexp     = -1._dp/zexm1
              zraut    = zxlb(jl)*(1._dp-(1._dp+zraut*ztmst*zexm1*zxlb(jl)      &
                                                       **zexm1)**zexp)
!--- Included for prognostic CDNC/IC scheme ----------------------------
              zrautn=zraut*7.7e9_dp*zrho(jl,jk)
              zself=1.289e10_dp*(zrho(jl,jk)*zxlb(jl))**2*ztmst
              zrautself=MIN(zrautn+zself,zcdnc(jl,jk))
              zcdnc(jl,jk)=MAX(zcdnc(jl,jk)-zrautself,cqtmin)
!--- End included for CDNC/IC scheme -----------------------------------
              zxlb(jl) = zxlb(jl)-zraut
              zrac1    = 6._dp*zxrp1*ztmst
              zrac1    = zxlb(jl)*(1._dp-EXP(-zrac1))
              zxlb(jl) = zxlb(jl)-zrac1
              zrac2    = 6._dp*zauloc*zrho(jl,jk)*zraut*ztmst
              zrac2    = zxlb(jl)*(1._dp-EXP(-zrac2))
              zxlb(jl) = zxlb(jl)-zrac2
              zrpr(jl) = zrpr(jl)+zclcaux(jl)*(zraut+zrac2)+zclcstar*zrac1
!--- Included for prognostic CDNC/IC scheme ----------------------------
              zxlbold=zxlb(jl)+zrac1+zrac2
              zraccn=(zrac1+zrac2)/(zxlbold+zeps)
              zraccn=MIN(zcdnc(jl,jk)*zraccn,zcdnc(jl,jk))
              zcdnc(jl,jk)=MAX(zcdnc(jl,jk)-zraccn,cqtmin)
              zrprn(jl)=zrautself+zraccn
!--- End included for CDNC/IC scheme -----------------------------------
           END IF
           !--- End changed for alternative autoconversion parameterisation -------
        ENDIF
        IF (zxib(jl) > cqtmin) THEN

!       7.2  Cold clouds:
!            Conversion of cloud ice to snow after Levkov et al. 1992:
!            Aggregation of ice crystals to snow and accretion of ice
!            by falling snow.
!            Accrection of cloud droplets by falling snow.
!            Effective radius of ice crystals after Moss (1995)
!
           zrieff    = 83.8_dp*(zxib(jl)*zrho(jl,jk)*1000._dp)**0.216_dp
           zrieff    = MIN(MAX(zrieff,ceffmin),ceffmax)
           zrih      = -2261._dp+SQRT(5113188._dp+2809._dp*zrieff*zrieff*zrieff)
           zris       = 1.e-6_dp*zrih**(1._dp/3._dp)
           zcolleffi = EXP(0.025_dp*(ztp1tmp(jl)-tmelt))
           zc1       = 17.5_dp*zrho(jl,jk)/crhoi*zqrho**0.33_dp
           zdt2      = -6._dp/zc1*LOG10(zris*1.e4_dp)
           zsaut     = ccsaut/zdt2
           zsaut     = zxib(jl)*(1._dp-1._dp/(1._dp+zsaut*ztmst*zxib(jl)))
           zxib(jl)  = zxib(jl)-zsaut
           zsaci1    = 0.0_dp
           zsaci2    = 0.0_dp
           zsacl1    = 0.0_dp
           zsacl2    = 0.0_dp
!--- Included for prognostic CDNC/IC scheme ----------------------------
           zsacl1in  = 0.0_dp      ! in-cloud values saved
           zsacl2in  = 0.0_dp      ! 
!--- End included for CDNC/IC scheme -----------------------------------
           IF (zxsp1 .GT. cqtmin .AND. zxlb(jl).GT.cqtmin .AND. zcdnc(jl,jk).GE. zcdnmin) THEN
!--- Included for prognostic CDNC/IC scheme ----------------------------
              IF (losacl) THEN
                 !
                 !uls - included size depedent accretion rate (Lohmann, JAS, 2004)
                 !
                 zscnc=0.75_dp*zrho(jl,jk)*zsaut/(api*zrsnow**3*crhosno)
                 zscnc=MAX(zsnowmin,zscnc)
                 zdw=MAX((6._dp*zrho(jl,jk)*zxlb(jl)/(api*rhoh2o*zcdnc(jl,jk))) &
                      **(1._dp/3._dp),1.e-6_dp)
                 zdplanar=MAX(20.e-6_dp,SQRT(zxsp1*1.e3_dp/(zscnc*3.8e-4_dp))*1.e-2_dp)
                 zusnow = 2.34_dp*(zdplanar*100._dp)**0.3_dp*(1.3_dp/zrho(jl,jk))**0.35_dp
                 zudrop = 1.19e4_dp*(0.5_dp*zdw*100._dp)**2*(1.3_dp/zrho(jl,jk))**0.35_dp
                 zstokes=MAX(2._dp*(zusnow-zudrop)*zudrop/(zdplanar*g),cqtmin)
                 zviscos=(1.512_dp + 0.0052_dp*(ptm1(jl,jk)-233.15_dp))*1.e-5_dp
                 zrey=MAX(zrho(jl,jk)*zdplanar*zusnow/zviscos,cqtmin)
                 !
                 zstcrit=1._dp
                 IF (zrey.LE.5._dp) zstcrit=5.52_dp*zrey**(-1.12_dp)
                 IF (zrey.GT.5._dp.AND.zrey.LT.40._dp) zstcrit=1.53_dp*zrey**(-0.325_dp)

                 zhelp=MAX(MIN(0.2_dp*(LOG10(zstokes)-LOG10(zstcrit) &
                      -2.236_dp)**2,1._dp-cqtmin),0._dp)
                 zcsacl1=SQRT(1._dp-zhelp)

                 IF (zrey.GE.40._dp) THEN
                    IF (zstokes .LE. 0.06_dp) THEN
                       zcsacl1=1.034_dp*zstokes**1.085_dp
                    ELSE IF (zstokes.GT.0.06_dp .AND. zstokes .LE. 0.25_dp) THEN
                       zcsacl1=0.787_dp*zstokes**0.988_dp
                    ELSE IF (zstokes.GT.0.25_dp .AND. zstokes .LE. 1._dp) THEN
                       zcsacl1=0.7475_dp*LOG10(zstokes)+0.65_dp
                    ELSE
                       zcsacl1=(zstokes+1.1_dp)**2/(zstokes+1.6_dp)**2
                    ENDIF
                 ENDIF
                 zcsacl=MAX(0.01_dp,MIN(zcsacl1,1._dp))
              ELSE
                 zcsacl=ccsacl
              ENDIF
!--- End included for CDNC/IC scheme -----------------------------------
              zlamsm    = (zxsp1/(api*crhosno*cn0s))**0.8125_dp
              zsaci1    = api*cn0s*3.078_dp*zlamsm*zqrho**0.5_dp
              zsacl1    = zxlb(jl)*(1._dp-EXP(-zsaci1*zcsacl*ztmst))
!--- Included for prognostic CDNC/IC scheme ----------------------------
              zsacl1in  = zsacl1
!--- End included for CDNC/IC scheme -----------------------------------
              zxlb(jl)  = zxlb(jl)-zsacl1
              zsacl1    = zclcstar*zsacl1
              zsaci1    = zsaci1*zcolleffi*ztmst
              zsaci1    = zxib(jl)*(1._dp-EXP(-zsaci1))
              zxib(jl)  = zxib(jl)-zsaci1
           END IF
           zxsp2        = zauloc*zrho(jl,jk)*zsaut
           IF (zxsp2 .GT. cqtmin .AND. zxlb(jl).GT.cqtmin .AND. zcdnc(jl,jk).GE. zcdnmin) THEN
              IF (losacl) THEN
                 !
                 !uls - included size depedent accretion rate (Lohmann, JAS, 2004)
                 !
                 zscnc=0.75_dp*zrho(jl,jk)*zsaci1/(api*zrsnow**3*crhosno)
                 zscnc=MAX(zsnowmin,zscnc)
                 zdw=MAX((6._dp*zrho(jl,jk)*zxlb(jl)/(api*rhoh2o*zcdnc(jl,jk))) &
                      **(1._dp/3._dp),1.e-6_dp)
                 zdplanar=MAX(20.e-6_dp,SQRT(zxsp2*1.e3_dp/(zscnc*3.8e-4_dp))*1.e-2_dp)
                 zusnow = 2.34_dp*(zdplanar*100._dp)**0.3_dp*(1.3_dp/zrho(jl,jk))**0.35_dp
                 zudrop = 1.19e4_dp*(0.5_dp*zdw*100._dp)**2*(1.3_dp/zrho(jl,jk))**0.35_dp
                 zstokes=MAX(2._dp*(zusnow-zudrop)*zudrop/(zdplanar*g),cqtmin)
                 zviscos=(1.512_dp + 0.0052_dp*(ptm1(jl,jk)-233.15_dp))*1.e-5_dp
                 zrey=MAX(zrho(jl,jk)*zdplanar*zusnow/zviscos,cqtmin)
                 !
                 zstcrit=1._dp
                 IF (zrey.LE.5._dp) zstcrit=5.52_dp*zrey**(-1.12_dp)
                 IF (zrey.GT.5._dp.AND.zrey.LT.40._dp) zstcrit=1.53_dp*zrey**(-0.325_dp)

                 zhelp=MAX(MIN(0.2_dp*(LOG10(zstokes)-LOG10(zstcrit) &
                      -2.236_dp)**2,1._dp-cqtmin),0._dp)
                 zcsacl1=SQRT(1._dp-zhelp)

                 IF (zrey.GE.40._dp) THEN
                    IF (zstokes .LE. 0.06_dp) THEN
                       zcsacl1=1.034_dp*zstokes**1.085_dp
                    ELSE IF (zstokes.GT.0.06_dp .AND. zstokes .LE. 0.25_dp) THEN
                       zcsacl1=0.787_dp*zstokes**0.988_dp
                    ELSE IF (zstokes.GT.0.25_dp .AND. zstokes .LE. 1._dp) THEN
                       zcsacl1=0.7475_dp*LOG10(zstokes)+0.65_dp
                    ELSE
                       zcsacl1=(zstokes+1.1_dp)**2/(zstokes+1.6_dp)**2
                    ENDIF
                 ENDIF
                 zcsacl=MAX(0.01_dp,MIN(zcsacl1,1._dp))
              ELSE
                 zcsacl=ccsacl
              ENDIF
              zlamsm    = (zxsp2/(api*crhosno*cn0s))**0.8125_dp
              zsaci2    = api*cn0s*3.078_dp*zlamsm*zqrho**0.5_dp
              zsacl2    = zxlb(jl)*(1._dp-EXP(-zsaci2*zcsacl*ztmst))
!--- Included for prognostic CDNC/IC scheme ----------------------------
              zsacl2in  = zsacl2
!--- End included for CDNC/IC scheme -----------------------------------
              zxlb(jl)  = zxlb(jl)-zsacl2
              zsacl2    = zclcaux(jl)*zsacl2
              zsaci2    = zsaci2*zcolleffi*ztmst
              zsaci2    = zxib(jl)*(1._dp-EXP(-zsaci2))
              zxib(jl)  = zxib(jl)-zsaci2
           END IF
           zsacl(jl)    = zsacl1+zsacl2
           zspr(jl)     = zspr(jl)+zclcaux(jl)*(zsaut+zsaci2)          &
                                  +zclcstar*zsaci1
!--- Included for prognostic CDNC/IC scheme ----------------------------
           zxlbold=zxlb(jl)+zsacl1in+zsacl2in
           IF (zxlb(jl) .GT. cqtmin) THEN
              zsacln(jl)=MAX( MIN(zcdnc(jl,jk)*(zsacl1in+zsacl2in)/(zxlbold+zeps),&
                   zcdnc(jl,jk)-zcdnmin ),0._dp )
           ELSE
              zsacln(jl)=MIN(zcdnc(jl,jk)*(zsacl1in+zsacl2in)/(zxlbold+zeps),&
                   zcdnc(jl,jk))
           END IF
           zcdnc(jl,jk)=MAX(zcdnc(jl,jk)-zsacln(jl),cqtmin)
!
!          secondary ice crystal production after Levkov et al. 1992
!          sink for snow, source for ice crystals
!
           zsecprod=0._dp
           IF ((zxsp1+zxsp2).GT.zepsec.AND.zxlb(jl).GT.zepsec.AND.   &
                ztp1tmp(jl).GT.265.2_dp.AND.ztp1tmp(jl).LT.270.2_dp) THEN
              zlams2=((zxsp1+zxsp2)/(api*crhosno*cn0s))**0.875_dp
              zj=zcsacl*api*zrho(jl,jk)*zxlb(jl)*cn0s*0.831_dp*zlams2*         &
                 (g*crhosno/(0.75_dp*zrho(jl,jk)*zcdi))**0.5_dp/zmw0            
              zpn=MAX(0.00285_dp*zj,0._dp)
              zsecprod=ztmst*zpn*zmi0/zrho(jl,jk)
              zsecprod=MAX(0._dp,MIN((zxsp1+zxsp2)/zrho(jl,jk),zsecprod))
              zxib(jl)=zxib(jl)+zsecprod
           END IF
!--- End included for CDNC/IC scheme -----------------------------------
           !
           zspr(jl)=MAX(zspr(jl)-zclcstar*zsecprod,0._dp)           
! storing of snow production before it is transformed into a flux
           prate_s(jl,jk) = zspr(jl) + zsacl(jl)         
        END IF
!--- Included for prognostic CDNC/IC scheme ----------------------------
        END IF
!--- End included for CDNC/IC scheme -----------------------------------
!
!       7.3 Updating precipitation fluxes. In the lowest layer (klev),
!           the sedimentation sink of cloud ice is balanced
!           by precipitation at the surface (through 'zzdrs').
!           Fraction of precipitating clouds (zclcpre) used for the
!           calculation of evaporation/sublimation of rain/snow in
!           the next layer
!
        zzdrr          = zcons2*zdp(jl)*zrpr(jl)
        zzdrs          = zcons2*zdp(jl)*(zspr(jl)+zsacl(jl))
        IF (jk .EQ. klev) THEN
           zzdrs       = zzdrs+zxiflux(jl)
           zcons       = zcons2*zdp(jl)/(zlsdcp(jl)-zlvdcp(jl))
           zsnmlt      = MIN(zxsec*zzdrs,zcons                         &
                                         *MAX(0._dp,(ztp1tmp(jl)-tmelt)))
           zzdrr       = zzdrr+zsnmlt
           zzdrs       = zzdrs-zsnmlt
           zsmlt(jl)   = zsmlt(jl)+zsnmlt/(zcons2*zdp(jl))
        END IF
        zpretot        = zrfl(jl)+zsfl(jl)
        zpredel        = zzdrr+zzdrs
        lo=(zpretot .GT. zpredel)
        zclcpre(jl)    = MERGE(zclcpre(jl),zclcaux(jl),lo)
        zpresum        = zpretot+zpredel
        IF (zpresum .GT. cqtmin) THEN
           zclcpre(jl) = MAX(zclcpre(jl),(zclcaux(jl)*zpredel          &
                                         +zclcpre(jl)*zpretot)/zpresum)
           zclcpre(jl) = MIN(zclcpre(jl),1.0_dp)
           zclcpre(jl) = MAX(zclcpre(jl),0.0_dp)
        ELSE
           zclcpre(jl) = 0.0_dp
        END IF

! rain and snow flux considering incoming rain, melting of snow, 
! droplet evaporation / sublimation , but no new production of rain or snow 
! in that layer....
! (neccessary for impaction scavenging)
        pfrain_no(jl,jk)   = zrfl(jl) - zcons2*zdp(jl)*zevp(jl)  
        pfsnow_no(jl,jk)   = zsfl(jl) - zcons2*zdp(jl)*zsub(jl)
! precipitating cloud cover of this layer is used for the next lower layer 
! to estimate the part of the cloud cover in which rain impacts
        pr_cover(jl,jk) = zclcpre(jl)

        zrfl(jl)       = zrfl(jl)+zzdrr-zcons2*zdp(jl)*zevp(jl)
        zsfl(jl)       = zsfl(jl)+zzdrs-zcons2*zdp(jl)*zsub(jl)
!
! rain and snow flux out of the bottom of this layer
        pfrain(jl,jk) = zrfl(jl)
        pfsnow(jl,jk) = zsfl(jl)
701  END DO
!
!     ------------------------------------------------------------------
!       8.    Updating tendencies of t, q, xl, xi and final cloud cover
!
     DO 811 jl = 1,kproma
!
!       8.10   Cloud cover scheme tendencies
!
        IF (jk.GE.ncctop) THEN
           locc            = zclcaux(jl) .GT. 0.0_dp
!
!          Source terms from convection
!          Skewness:
!
           zconvskew(jl)   = cbeta_cs * (pxtec(jl,jk)+pqtec(jl,jk))    &
                                  /pbetass(jl,jk)
           zconvskew(jl)   = MIN(zconvskew(jl),                        &
                                (cbeta_pq_max-pxskew(jl,jk))/zdtime)
!
!          Convective width now diagnosed, assuming 'a' unchanged:
!
           IF (pqm1(jl,jk) >= pbetass(jl,jk)) THEN
              zskewp1      = pxskew(jl,jk)+zconvskew(jl)*zdtime
              zbbap1       = zwide(jl)*(cbeta_pq+zskewp1)/             &
                                       (cbeta_pq+pxskew(jl,jk))
              zconvvar(jl) = (zbbap1-zwide(jl))/zdtime
           ELSE
              zconvvar(jl) = 0.0_dp
           ENDIF
!
!       8.11 Simple linearized effect of microphysics on skewness
!
           IF (pbetaa(jl,jk) < pbetass(jl,jk) .AND.                       &
               pbetab(jl,jk) > pbetass(jl,jk)) THEN
              zmdelb = (zxlte(jl)+zxite(jl))*ztmst                     &
                       -zrpr(jl)-zsacl(jl)-zspr(jl)+zcnd(jl)+zdep(jl)  &
                       +zgenti(jl)+zgentl(jl)
              zmdelb = MAX(0.0_dp,MIN(1.0_dp,-zmdelb/MAX(zepsec,zbetacl(jl))))
              zmdelb = (pbetass(jl,jk)-pbetab(jl,jk))*zmdelb
              zmqp1  = (pbetab(jl,jk)+zmdelb-pbetaa(jl,jk))            &
                        *cbeta_pq/(zbetaqt(jl)-pbetaa(jl,jk))          &
                                                          - cbeta_pq
              zmqp1  = MAX(MIN(zmqp1,cbeta_pq_max),cbeta_pq)
              zmicroskew(jl) = MIN(0.0_dp,(zmqp1-pxskew(jl,jk))/zdtime)
           ENDIF
!
!       8.2   New skewness and variance
!
           zxskewte(jl)    = zconvskew(jl)                             &
                             +zmicroskew(jl)+zturbskew(jl)
           zxvarte(jl)     = zconvvar(jl)+zturbvar(jl)
!
           zvarp1          = pxvar(jl,jk)+zxvarte(jl)*zdtime
           zskewp1         = pxskew(jl,jk)+zxskewte(jl)*zdtime
!
           pxskew(jl,jk)   = MAX(MIN(zskewp1,cbeta_pq_max),cbeta_pq)
           zvarmx          = zbetaqt(jl)*(1._dp+pxskew(jl,jk)/cbeta_pq)
           pxvar(jl,jk)    = MAX(MIN(zvarp1,zvarmx),zvartg(jl))
!
        ENDIF !jk >= ncctop ! mz_jd_20161011
!
!       8.3   Tendencies of thermodynamic variables
!             Attn: The terms zxisub and zximlt do not appear in
!                   pxite because these processes have already been
!                   included in pxite via changes in cloud ice
!                   sedimentation (see 3.1, 3.2 and 4)
!
        pqte(jl,jk)  = pqte(jl,jk)                                     &
                        +(-zcnd(jl)-zgentl(jl)+zevp(jl)+zxlevap(jl)    &
                          -zdep(jl)-zgenti(jl)+zsub(jl)+zxievap(jl)    &
                                          +zxisub(jl))/ztmst
        ptte(jl,jk)  = ptte(jl,jk)+(zlvdcp(jl)                         &
                        *(zcnd(jl)+zgentl(jl)-zevp(jl)-zxlevap(jl))    &
                                  +zlsdcp(jl)                          &
                        *(zdep(jl)+zgenti(jl)-zsub(jl)-zxievap(jl)     &
                        -zxisub(jl))+(zlsdcp(jl)-zlvdcp(jl))           &
                        *(-zsmlt(jl)-zimlt(jl)-zximlt(jl)+zfrl(jl)     &
                                           +zsacl(jl)))/ztmst
        pxlte(jl,jk) = pxlte(jl,jk)+zxlte(jl)                          &
                        +(zimlt(jl)+zximlt(jl)-zfrl(jl)-zrpr(jl)       &
                        -zsacl(jl)+zcnd(jl)+zgentl(jl)-zxlevap(jl))    &
                                                       /ztmst
        pxite(jl,jk) = pxite(jl,jk)+zxite(jl)+(zfrl(jl)-zspr(jl)       &
                          +zdep(jl)+zgenti(jl)-zxievap(jl))/ztmst
        ztp1         = ptm1(jl,jk)+ptte(jl,jk)*ztmst
        zqp1         = pqm1(jl,jk)+pqte(jl,jk)*ztmst
        zxlp1        = pxlm1(jl,jk)+pxlte(jl,jk)*ztmst
        zxip1        = pxim1(jl,jk)+pxite(jl,jk)*ztmst
!
!--- Included for prognostic CDNC/IC scheme ----------------------------

        !--- Calculate new total tendency of CDNC:

        pxtte(jl,jk,idt_cdnc) = pxtte(jl,jk,idt_cdnc) + ( zcdnc(jl,jk) /  &
                                ( zrho(jl,jk) * 1000._dp / M_air) -       &
                                ( pxtm1(jl,jk,idt_cdnc) +                 &
                                  pxtte(jl,jk,idt_cdnc) * ztmst ) ) / ztmst
        
        !--- Update CDNC for radiation:

        pacdnc(jl,jk)=zcdnc(jl,jk)

        !--- Diagnostics:

        qaut(jl,jk,jrow) = -zrprn(jl)  / ztmst
        qfre(jl,jk,jrow) = -zfrln(jl)  / ztmst
        qacc(jl,jk,jrow) = -zsacln(jl) / ztmst

        cloud_tm1(jl,jk,jrow) = paclc(jl,jk)

        IF (zxlb(jl)>zeps.AND.zcdnc(jl,jk).GE.zcdnmin)THEN
           cdnc_acc(jl,jk,jrow)   = zcdnc(jl,jk)
           zcdnc_burden(jl) = zcdnc_burden(jl)+zcdnc(jl,jk)*zdz(jl)

           !---- CDNC and burden averaged over cloudy and cloud-free periods

           cdnc(jl,jk,jrow) = zcdnc(jl,jk)*zclcaux(jl)
           cdnc_burden(jl,jrow) = cdnc_burden(jl,jrow) + &
                                  zcdnc(jl,jk) * zdz(jl) * zclcaux(jl)

           !--- In-cloud effective radius [um]:

          IF (slm(jl).GT.0._dp .AND.(.NOT. glac(jl).GT.0._dp)) THEN 
              zkap=1.143_dp
           ELSE
              zkap=1.07_dp
           ENDIF

           zreffl=1.E6_dp*zkap*((3._dp/(4._dp*api*rhoh2o))*zxlb(jl)*zrho(jl,jk)/zcdnc(jl,jk))**(1._dp/3._dp)

           reffl(jl,jk,jrow) = zreffl
        END IF

!--- End included for CDNC/IC scheme -----------------------------------

!       8.4   Corrections: Avoid negative cloud water/ice
!
        zxlold = zxlp1
        lo             = (zxlp1 .LT. ccwmin)
        zxlp1          = MERGE(0.0_dp,zxlp1,lo)
        zdxlcor        = (zxlp1-zxlold)/ztmst
        zxiold         = zxip1
        lo1            = (zxip1 .LT. ccwmin)
        zxip1          = MERGE(0.0_dp,zxip1,lo1)
        zdxicor        = (zxip1-zxiold)/ztmst
        paclc(jl,jk)   = MERGE(0.0_dp,paclc(jl,jk),lo.AND.lo1)
        paclcac(jl,jk) = paclcac(jl,jk)+paclc(jl,jk)*zdtime
        pxlte(jl,jk)   = pxlte(jl,jk)+zdxlcor
        pxite(jl,jk)   = pxite(jl,jk)+zdxicor
        pqte(jl,jk)    = pqte(jl,jk)-zdxlcor-zdxicor
        ptte(jl,jk)    = ptte(jl,jk)+zlvdcp(jl)*zdxlcor                &
                                    +zlsdcp(jl)*zdxicor
!--- Included for prognostic CDNC/IC scheme ----------------------------
        zcdnold        = zcdnc(jl,jk)
        zcdnp          = MERGE(0.0_dp,zcdnc(jl,jk),lo)
        zdnlcor        = (zcdnp-zcdnold)/ztmst
        pxtte(jl,jk,idt_cdnc) = pxtte(jl,jk,idt_cdnc) + &
                                zdnlcor /( zrho(jl,jk) * 1000._dp / M_air)
!--- End included for CDNC/IC scheme -----------------------------------
!
811  END DO
831 END DO    ! Vertical loop
!
!--- Included for prognostic CDNC/IC scheme ----------------------------
  DO jl=1, kproma
     IF (zcdnc_burden(jl)>zepsec) THEN 
       cdnc_burden_acc(jl,jrow) = zcdnc_burden(jl)
     END IF
  END DO 
!--- End included for CDNC/IC scheme -----------------------------------
!
!     ------------------------------------------------------------------
!
!       10.    Diagnostics
!
!       10.1   Accumulated precipitation at the surface
!
  DO 911 jl    = 1,kproma
     prsfl(jl) = zrfl(jl)
     pssfl(jl) = zsfl(jl)
     paprl(jl) = paprl(jl)+zdtime*(prsfl(jl)+pssfl(jl))
     paprs(jl) = paprs(jl)+zdtime*pssfl(jl)
911 END DO
!
!       10.2   Total cloud cover
!
    DO 921 jl    = 1,kproma
      zclcov(jl) = 1.0_dp-paclc(jl,1)
921 END DO
    DO 923 jk      = 2,klev
      DO 922 jl    = 1,kproma
        zclcov(jl) = zclcov(jl)*(1._dp-MAX(paclc(jl,jk),paclc(jl,jk-1)))  &
                               /(1._dp-MIN(paclc(jl,jk-1),zxsec))
922   END DO
923 END DO
    DO 924 jl     = 1,kproma
      zclcov(jl)  = 1.0_dp-zclcov(jl)
      paclcov(jl) = paclcov(jl)+zdtime*zclcov(jl)
924 END DO
!
!       10.3   Vertical integrals of humidity, cloud water and cloud ice
!
    DO 931 jl   = 1,kproma
      zqvi(jl)  = 0.0_dp
      zxlvi(jl) = 0.0_dp
      zxivi(jl) = 0.0_dp
931 END DO
!
    DO 933 jk     = ktdia,klev
      DO 932 jl   = 1,kproma
        zdpg      = (paphm1(jl,jk+1)-paphm1(jl,jk))/g
        zqvi(jl)  = zqvi(jl)+pqm1(jl,jk)*zdpg
        zxlvi(jl) = zxlvi(jl)+pxlm1(jl,jk)*zdpg
        zxivi(jl) = zxivi(jl)+pxim1(jl,jk)*zdpg
932   END DO
933 END DO
!
    DO 934 jl   = 1,kproma
      pqvi(jl)  = pqvi(jl)+zdtime*zqvi(jl)
      pxlvi(jl) = pxlvi(jl)+zdtime*zxlvi(jl)
      pxivi(jl) = pxivi(jl)+zdtime*zxivi(jl)
934 END DO
!
  RETURN

END SUBROUTINE cloud_cdnc

!========================================================
SUBROUTINE cloud_cdnc_icnc(kproma, kbdim, ktdia, klev, klevp1, ztmst   &
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
                          ! mz_ht_20070629+
                         , status_string,      lcover                  &
                         , slm,      glac,     pcdncact                &
                         , zna                                         &
                         , plwc,     piwc                              &
                         , pfrain,   pfsnow                            &
                         , pfrain_no,pfsnow_no                         &
                         , prate_r,  prate_s                           &
                         , prevap,   pssubl                            &
                         , pr_cover, pcond                             &
                         , pimelt,   pisedi                            &
                         , sigma,    wetradius                         &
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
!
!+    
  USE messy_cloud_ori,  ONLY : cqtmin, tmelt, cvtfall, crhosno, cn0s    &
                             , cthomi, csecfrl, ncctop, cvarmin         &
                             , cbeta_pq, cbeta_pq_max, nbetaq, cbetaqs  &
                             , rbetak, nbetax, tbetai0, tbetai1, cauloc &
                             , clmax, clmin, jbmin, jbmax, lonacc       &
                             , ccraut, crhoi, ccsaut  &
                             , ccsacl, cbeta_cs, LOOKUPOVERFLOW

  USE messy_main_constants_mem, ONLY: API => PI, g, M_air, ak => k_B
  USE messy_main_tools,         ONLY: jptlucu1, jptlucu2, &
                                      tlucua, tlucuaw, tlucub
    

  IMPLICIT NONE
!
  INTEGER :: kbdim, klevp1, klev, kproma, ktdia, krow, ktrac, jl, jk, it

  CHARACTER(LEN=32) :: status_string
  LOGICAL  :: lcover
  REAL(dp) :: slm(kbdim), glac(kbdim)

!
  REAL(dp)::    paphm1(kbdim,klevp1),pvervel(kbdim,klev)                   &
           ,papm1(kbdim,klev)   ,pqm1(kbdim,klev)   ,papp1(kbdim,klev) &
           ,ptm1(kbdim,klev)    ,ptvm1(kbdim,klev)  ,pxlm1(kbdim,klev) &
           ,pxim1(kbdim,klev)   ,pxtec(kbdim,klev)  ,pqtec(kbdim,klev) &
           ,pxvar(kbdim,klev)   ,pxskew(kbdim,klev)                    &
           ,pbetaa(kbdim,klev)  ,pbetab(kbdim,klev)                    &
           ,pvdiffp(kbdim,klev) ,phmixtau(kbdim,klev)                  &
           ,pvmixtau(kbdim,klev),pgeo(kbdim,klev)   ,pbetass(kbdim,klev)
  REAL(dp)::    pxlvi(kbdim)        ,pxivi(kbdim)
  REAL(dp)::    paclc(kbdim,klev)   ,paclcac(kbdim,klev)
  REAL(dp)::    pacdnc(kbdim,klev)  ,prelhum(kbdim,klev)
  REAL(dp)::    paclcov(kbdim)      ,paprl(kbdim)                          &
           ,pqvi(kbdim)         ,pssfl(kbdim)
  REAL(dp)::    ptte(kbdim,klev)    ,pqte(kbdim,klev)
  REAL(dp)::    pxlte(kbdim,klev)   ,pxite(kbdim,klev)
  REAL(dp)::    paprs(kbdim)        ,prsfl(kbdim)

!---Included for in-cloud scavenging (Philip Stier, 28/03/01):----------
  REAL(dp) ::    pxtm1(kbdim,klev,ktrac), pxtte(kbdim,klev,ktrac) 
!---End Included for scavenging-----------------------------------------
!
!   Temporary arrays
!
  REAL(dp)   :: zclcpre(kbdim)                                             &
           ,zcnd(kbdim)         ,zdep(kbdim)        ,zdp(kbdim)        &
           ,                     zxievap(kbdim)     ,zxlevap(kbdim)    &
           ,zfrl(kbdim)         ,zimlt(kbdim)       ,zsmlt(kbdim)      &
           ,                     zspr(kbdim)                           &
           ,zxlte(kbdim)        ,zxite(kbdim)       ,zxiflux(kbdim)    &
           ,zsacl(kbdim)        ,zdz(kbdim)                            &
           ,zlsdcp(kbdim)       ,zlvdcp(kbdim)      ,zximlt(kbdim)     &
           ,ztp1tmp(kbdim)      ,zqp1tmp(kbdim)     ,zxisub(kbdim)     &
           ,zxlb(kbdim)         ,zxib(kbdim)                           &
           ,zrho(kbdim,klev)    ,zclcov(kbdim)      ,zclcaux(kbdim)    &
           ,zqvi(kbdim)         ,zxlvi(kbdim)       ,zxivi(kbdim)      &
           ,zbetaqt(kbdim)      ,zwide(kbdim)                          &
           ,zbetacl(kbdim)      ,zturbvar(kbdim)    ,zturbskew(kbdim)  &
           ,zconvvar(kbdim)     ,zconvskew(kbdim)   ,zvartg(kbdim)     &
           ,zmicroskew(kbdim)   ,zgenti(kbdim)      ,zgentl(kbdim)     &
           ,zxvarte(kbdim)      ,zxskewte(kbdim)                       &
           ,zgeoh(kbdim,klevp1)
!
  REAL(dp)   :: zbqp1,zbbap1,zbap1,ztt,zgent,zqcdif,zqp1,ztp1,             &
            zqp1b,zskew,zbetai0,zbetai1,zskewp1,zvarp1,zifrac,zvarmx
  REAL(dp)   :: zdtime,zauloc
  INTEGER:: iqidx,ixidx, jb
  LOGICAL   lo,lo1,lo2,locc
  INTEGER   knvb(kbdim)

!--- Included for dust emissions (Philip Stier 10/01/02)-----------------
  INTEGER :: jrow
!--- End Included for dust emissions ------------------------------------

!---Included for in-cloud scavenging (Philip Stier, 28/03/01):----------
  REAL(dp)    :: zmratepr(kbdim,klev), & ! Rain formation rate in cloudy part
                                     ! of the grid box [kg/kg]
             zmsnowacl(kbdim,klev)   ! Accretion rate of snow with cloud
                                     ! droplets in cloudy part of the 
                                     ! grid box  [kg/kg]
!---End Included for scavenging-----------------------------------------

!--- Included for prognostic CDNC/IC scheme ----------------------------

  LOGICAL, PARAMETER :: losacl=.TRUE.  ! Size-dependent acretion rate
                                       ! (Lohmann, JAS, 2004)

  REAL(dp)    :: pcdncact(kbdim,klev),     ptkem1(kbdim,klev),             &
             pcvcbot(kbdim),           pwcape(kbdim),                  &
             pxtecl(kbdim,klev),       pxteci(kbdim,klev)

  INTEGER :: itop(kbdim,klev),         ibas(kbdim,klev),               &
             iclbas (kbdim)  

  REAL(dp)    :: zeps, zkap, zreffl, zicemin, zqlnuccv,    &
             znicenew

  REAL(dp)    :: zrprn(kbdim),             zsprn(kbdim,klev),              &
             zsacln(kbdim),            zfrln(kbdim),                   &
             zcdnc(kbdim,klev),     & !
             !zesw_2d(kbdim,klev),   & !
             zsusatw_2d(kbdim,klev),& ! supersat. with respect to water
             zcdnc_burden(kbdim),   &
             zna(kbdim,klev),       &
             zqlnuc(kbdim,klev),    & !
             zqlnuccvh(kbdim,klev)    ! Nucleated CDNC from convective detrainment 

 ! Temporary fields
  INTEGER  :: it1
  REAL(dp) :: zesw,  zvervc, zcdncnew, zcdnmin, zxsec, zcons, ztdif,  &
              zesi , zcfac4c,zsubi, zzepr, zxip1, zqsm1, zqst1, zrelhum, &
              zcond, zes, zlc, zf1, zraut, zrautself, zrautn, &
              zrieff, zhelp, zhelp1, zhelp2, zhelp3, zstokes, zstcrit, zcsacl1, zsecprod, zmdelb, &
              zepsec, zxilb, zmqp1, zxlp1, zxrp1, zxsp1, zqrho, zqsec, ztmst,&
              zcons1, zcons2, zradl, zsigmaw, zdisp, zdw0, zcsacl,  zrsnow, &
              zmw0, zmi0,   zqsw , zrcp, zsnmlt, zsnowmin, zcdi, zximelt, zicncmelt, &
              zclcstar, zdpg, zqsi, zb1, zsusati, zb2, zcoeff, zclambs,&
              zzeps, zesat, zsusatw, zdv, zast, zbst, zxifall, zal1, zal2,  &
              zxim1evp, zxlm1evp,zxidt, zxldt, zxidtstar, zxldtstar, zdqsdt, &
              zqvdt, zdtdt, zdtdtstar, zdqsat, zqtau, zpp, zqq, zeta, &
              zprod, zaa, zgtp, zcor, zqsp1tmp, zoversat, zrhtest, zqcon, zdepcor, zcndcor, &
              zlcdqsdt, zzevp, zdepos, znfrl, zexm1, zexp, zrac1, zrac2, &
              zxlbold, zraccn, zrih, zris, zcolleffi, zc1, zsaut, zsaci1, zsaci2, &
              zsacl1,zsacl2, zsacl1in, zsacl2in, zscnc, zdw,  zusnow, zudrop, zviscos,&
              zrey, zlamsm, zxsp2, zself, zdt2, zdplanar, zlams2, zj, zpn, zzdrr, zzdrs, &
              zpretot, zpredel, zpresum, zxlold, zdxlcor, zxiold, zdxicor, zcdnold, &
              zdnlcor, zcdnp, zrid

  REAL(dp)    :: zicnc(kbdim,klev),     & !
             zicnc_burden(kbdim),   & !
             zascs(kbdim,klev),     & !
             zninucl(kbdim, klev),  & ! number conc. of newly nucleated IC
             zqinucl(kbdim, klev),  & ! mixing ratio of newly nucleated IC
             zri(kbdim, klev),      & ! size of newly nucleated IC
             zndusol(kbdim,klev), znduinsolai(kbdim,klev),  &
             znduinsolci(kbdim,klev),  znbcsol(kbdim,klev), &
             znbcinsol(kbdim,klev),                         &
             zreffct(kbdim)

  !--- Included for cirrus scheme -------------------------------------------
  REAL(dp) :: znicex(kbdim)     ! 
  REAL(dp) :: zsusatix(kbdim)   ! ice supersaturation pressure
  REAL(dp) :: zvervx(kbdim)     ! updraft [cm/s]
  REAL(dp) :: zapnx(kbdim,nfrzmod) ! aerosol number available for homogeneous freezing [1/cm3]
  REAL(dp) :: zaprx(kbdim,nfrzmod) ! radius of aerosols available for homogeneous freezing  [cm]
  REAL(dp) :: zapsigx(kbdim,nfrzmod)! std. dev. of aerosols available for homogeneous freezing 

  LOGICAL, PARAMETER :: nosize = .false. ! .true. ---> aerosol size effects
                                         ! are ignored for homogeneous
                                         ! freezing

    ! !!! lhet = .true. NOT YET IMPLEMENTED !!!

  LOGICAL, PARAMETER :: lhet = .FALSE.   ! .true. ---> heterogeneous
                                         ! freezing of aerosol particles
                                         ! below 235K is considered
  !--- End included for cirrus scheme -------------------------------------------

  !cyjs-----------------------------------------
  REAL(dp) :: zxifluxn(klev),zicesed(klev),zaaa
  REAL(dp) :: zsigma_m2, zmm(klev),zmmean(klev), zalfased, zbetased
  REAL(dp) :: zxifallmc(klev), zxifallnc(klev), zxised(klev)
  !cyje-----------------------------------------

  REAL(dp)            :: zmi             ! assumed mass of ice crystals with 
                                     ! corresponding volume mean radius zcri 
  REAL(dp), PARAMETER :: zcri=10.E-6_dp     ! to estimate the number of produced  
                                     ! cloud droplets from ice melting in  
                                     ! case of licnc=.FALSE. [m]=> 10 um

!--- End Included for CDNC -------------------------------------------
  REAL(dp) :: zdfarduci, zfrzcntdu, zfrzcnt, zccbcki, zccduai
  REAL(dp) :: zdfarbcki, zdfarduai, zccduci, zfrzcntbc, zetaair
  REAL(dp) :: znaimmdu, znaimmbc, zfrzimm
  REAL(dp) :: zxibold, zicnold, zicnp, zdnicor
  REAL(dp) :: zqsp1tmpw, zoversatw, znaerinsol, zfracdusol
  REAL(dp) :: zfracduinsolci, zfracduinsolai, zfracbcsol, zfracbcinsol
  REAL(dp) :: zverv, zvth, zfuchs, zfre, zre, zcorw, ztc, znicenew1
  REAL(dp) :: zrdryki, zrdryai, zrdryci
  REAL(dp) :: zsprn1, zsecprodn, zsprnself
  REAL(dp) :: zfall, zrhoice, zka, zkb, zalpha, zxmw
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev)         :: plwc,     piwc
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev), TARGET :: pfrain,   pfsnow
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev)         :: pfrain_no,pfsnow_no
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev), TARGET :: prevap,   pssubl
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev), TARGET :: prate_r,  prate_s 
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev)         :: pr_cover
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev)         :: pimelt,   pisedi
  REAL(dp), INTENT(INOUT), DIMENSION(kbdim, klev)         :: pcond 

  REAL(dp), POINTER, DIMENSION(:) :: zrfl, zsfl, zsub, zevp, zrpr

  REAL(dp), INTENT(IN), DIMENSION(7)            :: sigma
  REAL(dp), INTENT(IN), DIMENSION(kbdim,klev,7) :: wetradius
!--- Included for dust emissions (Philip Stier 10/01/02)-----------------
  jrow = krow

!--- End Included for dust emissions ------------------------------------

!---Included for in-cloud scavenging (Philip Stier, 28/03/01):----------
  zmratepr(:,:) = 0._dp
  zmsnowacl(:,:)= 0._dp
!---End Included for scavenging-----------------------------------------
!--- Included for prognostic CDNC/IC scheme ----------------------------
  zqinucl(:,:) = 0._dp
  zninucl(:,:) = 0._dp
  zreffct(:)   = 0._dp
!--- End Included for CDNC -------------------------------------------
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

! Executable statements
!
  lookupoverflow = .FALSE.
!
!   Security parameters
!
  zepsec = 1.0e-12_dp
  zxsec  = 1.0_dp-zepsec
  zqsec  = 1.0_dp-cqtmin
!
!   Computational constants
!
  zdtime  = ztmst / 2._dp
  zcons1 = cpd*vtmpc2
  zcons2 = 1._dp/(ztmst*g)
!--- Included for prognostic CDNC/IC scheme ----------------------------
  zeps=EPSILON(1.0_dp)*10._dp

  zmi=4._dp/3._dp*zcri**3._dp*api*crhoi

  zsigma_m2=(LOG(3.25_dp))**2
!
  zmm(:)=zmi
  zmmean(:)=zmi
  zxifallmc(:)=0._dp
  zxifallnc(:)=0._dp
  zxised(:)=0._dp
  zxifluxn(:)=0._dp
  zicesed(:)=0._dp
  zfall=8._dp
  zrhoice=925._dp
!

  zcdnmin=40.e6_dp      !
  zsigmaw=0.28_dp
  zdisp=EXP(zsigmaw**2/2._dp)
  zdw0=10.e-6_dp*zdisp ! modal diameter times dispersion parameter
  zcdnc_burden(1:kproma)=0.0_dp
  zicemin=10._dp
  zcsacl = 0.01_dp

  zqlnuc(1:kproma,:) = 0.0_dp
  zicnc_burden(1:kproma)=0.0_dp
  zri(1:kproma,:) = 1.e-6_dp
  zsnowmin=1._dp
  zrsnow=1.e-3_dp
  zcdi=0.6_dp         !
  zmw0=4.19e-12_dp     !
  zmi0=1.e-12_dp       !
  zka=0.024_dp
  zkb=1.38e-23_dp
  zalpha=0.5_dp        ! deposition coefficent alpha
  zxmw=2.992e-26_dp    ! mass of a h2o molecule [kg]

  cdnc_burden(1:kproma,jrow) = 0.0_dp
  cdnc_acc(1:kproma,:,jrow)  = 0._dp
  cdnc(1:kproma,:,jrow)      = 0._dp
  qnuc(1:kproma,:,jrow)      = 0._dp
!--- End included for CDNC/IC scheme -----------------------------------
!
!     ------------------------------------------------------------------
!
!       1.   Top boundary conditions, air density and geopotential
!            height at half levels
!
!       1.1   Set to zero precipitation fluxes etc.
!
  DO 111 jl = 1,kproma
     zclcpre(jl)   = 0.0_dp
     zxiflux(jl)   = 0.0_dp
111 END DO
!
!       1.2   Air density
!
  DO 122 jk = ktdia,klev
     DO 121 jl = 1,kproma
        zrho(jl,jk)   = papm1(jl,jk)/(rd*ptvm1(jl,jk))
!-----------------Added by Junhua Zhang for CONV Micro--------------------
        IF(ncvmicro>0) THEN
         pxtecl(jl,jk)  = MAX(pxtecl(jl,jk),0.0_dp)
         pxteci(jl,jk)  = MAX(pxteci(jl,jk),0.0_dp)
        ELSE
         pxtec(jl,jk)  = MAX(pxtec(jl,jk),0.0_dp)
        ENDIF
!------------------------------End----------------------------------------
        pqtec(jl,jk)  = MAX(pqtec(jl,jk),0.0_dp)
!--- Included for prognostic CDNC/IC scheme ----------------------------

        ! calculate cloud droplet number concentration 
        zcdnc(jl,jk)=MAX( ( pxtm1(jl,jk,idt_cdnc) + &
                            pxtte(jl,jk,idt_cdnc) * ztmst ) &
                            / M_air *1000._dp * zrho(jl,jk) , cqtmin)

        ! calculate ice crystal number
        zicnc(jl,jk)=MAX( ( pxtm1(jl,jk,idt_icnc) + &
                            pxtte(jl,jk,idt_icnc) * ztmst ) / &
                            M_air * 1000._dp * zrho(jl,jk) , cqtmin)
        ! Provide soluble aerosol number number concentratios for the cirrus scheme:

        zascs(jl,jk)=MAX( (pxtm1(jl,jk,idt_nks) + pxtm1(jl,jk,idt_nas) + &
                           pxtm1(jl,jk,idt_ncs) +                        &
                          (pxtte(jl,jk,idt_nks) + pxtte(jl,jk,idt_nas) + &
                           pxtte(jl,jk,idt_ncs) ) * ztmst ) / M_air    * &
                           1000._dp * zrho(jl,jk), 10.E6_dp)
        ! corinna: provide total aerosol, dust and BC number concentrations
        ! insoluble modes for contact freezing
        znduinsolai(jl,jk)= (pxtm1(jl,jk,idt_nai) + pxtte(jl,jk,idt_nai) * &
                             ztmst) / M_air * 1000._dp * zrho(jl,jk)
        znduinsolci(jl,jk)= (pxtm1(jl,jk,idt_nci) + pxtte(jl,jk,idt_nci) * &
                             ztmst) / M_air * 1000._dp * zrho(jl,jk)
        znbcinsol(jl,jk)= MAX((pxtm1(jl,jk,idt_mbcki) + pxtte(jl,jk,idt_mbcki)*&
                            ztmst) /                                          &
                           (pxtm1(jl,jk,idt_mbcki) + pxtte(jl,jk,idt_mbcki) * &
                            ztmst +                                           &
                            pxtm1(jl,jk,idt_mocki) + pxtte(jl,jk,idt_mocki) * &
                            ztmst + zeps),0._dp)**(2._dp/3._dp)             * &
                           (pxtm1(jl,jk,idt_nki)   + pxtte(jl,jk,idt_nki)   * &
                           ztmst) / M_air * 1000._dp * zrho(jl,jk) 
        ! soluble modes for immersion freezing: hydrophilic BC
        ! changed S.K.
        zhelp1 = ( (pxtm1(jl,jk,idt_mbcks) + pxtte(jl,jk,idt_mbcks) * &
                    ztmst) /                                          &
                   (pxtm1(jl,jk,idt_mbcks) + pxtte(jl,jk,idt_mbcks) * &
                    ztmst +                                           &
                    pxtm1(jl,jk,idt_mocks) + pxtte(jl,jk,idt_mocks) * &
                    ztmst +                                           &
                    pxtm1(jl,jk,idt_ms4ks) + pxtte(jl,jk,idt_ms4ks) * &
                    ztmst + zeps) )
        
        zhelp2 = ( (pxtm1(jl,jk,idt_mbcas) + pxtte(jl,jk,idt_mbcas) * &
                    ztmst) /                                          &
                   (pxtm1(jl,jk,idt_mbcas) + pxtte(jl,jk,idt_mbcas) * &
                    ztmst +                                           &
                    pxtm1(jl,jk,idt_mocas) + pxtte(jl,jk,idt_mocas) * &
                    ztmst +                                           &
                    pxtm1(jl,jk,idt_ms4as) + pxtte(jl,jk,idt_ms4as) * &
                    ztmst +                                           &
                    pxtm1(jl,jk,idt_mssas) + pxtte(jl,jk,idt_mssas) * &
                    ztmst +                                           &
                    pxtm1(jl,jk,idt_mduas) + pxtte(jl,jk,idt_mduas) * &
                    ztmst + zeps) )

        zhelp3 = ( (pxtm1(jl,jk,idt_mbccs) + pxtte(jl,jk,idt_mbccs) * &
                    ztmst) /                                          &
                   (pxtm1(jl,jk,idt_mbccs) + pxtte(jl,jk,idt_mbccs) * &
                    ztmst +                                           &
                    pxtm1(jl,jk,idt_moccs) + pxtte(jl,jk,idt_moccs) * &
                    ztmst +                                           &
                    pxtm1(jl,jk,idt_ms4cs) + pxtte(jl,jk,idt_ms4cs) * &
                    ztmst +                                           &
                    pxtm1(jl,jk,idt_msscs) + pxtte(jl,jk,idt_msscs) * &
                    ztmst +                                           &
                    pxtm1(jl,jk,idt_mducs) + pxtte(jl,jk,idt_mducs) * &
                    ztmst + zeps) )

        IF (zhelp1 .LT. zeps) zhelp1 = 0.0_dp
        IF (zhelp2 .LT. zeps) zhelp2 = 0.0_dp
        IF (zhelp3 .LT. zeps) zhelp3 = 0.0_dp

        znbcsol(jl,jk)=( zhelp1**(2._dp/3._dp) * &
                         pxtm1(jl,jk,idt_nks) + pxtte(jl,jk,idt_nks) * ztmst + &
                         zhelp2**(2._dp/3._dp) * &
                         pxtm1(jl,jk,idt_nas) + pxtte(jl,jk,idt_nas) * ztmst + &
                         zhelp3**(2._dp/3._dp) * &
                         pxtm1(jl,jk,idt_ncs) + pxtte(jl,jk,idt_ncs) * ztmst)/ &
                         M_air * 1000._dp * zrho(jl,jk) 
        
        zhelp1 = ( (pxtm1(jl,jk,idt_mduas) + pxtte(jl,jk,idt_mduas) * &
                    ztmst) /                                          &
                   (pxtm1(jl,jk,idt_mbcas) + pxtte(jl,jk,idt_mbcas) * &
                    ztmst +                                           &
                    pxtm1(jl,jk,idt_mocas) + pxtte(jl,jk,idt_mocas) * &
                    ztmst +                                           &
                    pxtm1(jl,jk,idt_ms4as) + pxtte(jl,jk,idt_ms4as) * &
                    ztmst +                                           &
                    pxtm1(jl,jk,idt_mssas) + pxtte(jl,jk,idt_mssas) * &
                    ztmst +                                           &
                    pxtm1(jl,jk,idt_mduas) + pxtte(jl,jk,idt_mduas) * &
                    ztmst + zeps) )

        zhelp2 = ( (pxtm1(jl,jk,idt_mducs) + pxtte(jl,jk,idt_mducs) * &
                    ztmst) /                                          &
                   (pxtm1(jl,jk,idt_mbccs) + pxtte(jl,jk,idt_mbccs) * &
                    ztmst +                                           &
                    pxtm1(jl,jk,idt_moccs) + pxtte(jl,jk,idt_moccs) * &
                    ztmst +                                           &
                    pxtm1(jl,jk,idt_ms4cs) + pxtte(jl,jk,idt_ms4cs) * &
                    ztmst +                                           &
                    pxtm1(jl,jk,idt_msscs) + pxtte(jl,jk,idt_msscs) * &
                    ztmst +                                           &
                    pxtm1(jl,jk,idt_mducs) + pxtte(jl,jk,idt_mducs) * &
                    ztmst + zeps) )

        IF (zhelp1 .LT. zeps) zhelp1 = 0.0_dp
        IF (zhelp2 .LT. zeps) zhelp2 = 0.0_dp

        zndusol(jl,jk)=( zhelp1**(2._dp/3._dp) *                               &
                         pxtm1(jl,jk,idt_nas) + pxtte(jl,jk,idt_nas) * ztmst + &
                         zhelp2**(2._dp/3._dp) *                               &
                         pxtm1(jl,jk,idt_ncs) + pxtte(jl,jk,idt_ncs) * ztmst)/ &
                         M_air * 1000._dp * zrho(jl,jk)
121  END DO
122 END DO
!
!       1.3   Geopotential at half levels
!
  DO 132 jk = 2,klev
     DO 131 jl = 1,kproma
        zgeoh(jl,jk)   = 0.5_dp*(pgeo(jl,jk)+pgeo(jl,jk-1))
131  END DO
132 END DO
   DO 133 jl = 1,kproma
      zgeoh(jl,1)      = pgeo(jl,1)+(pgeo(jl,1)-zgeoh(jl,2))
      zgeoh(jl,klevp1) = 0.0_dp
133 END DO

!--- Included for prognostic CDNC/IC scheme ----------------------------
  !
  !       1.4  CLOUD TOPS
  !
  DO 230 jl=1,kproma
     IF (paclc(jl,1) .GE. zepsec) THEN
        itop(jl,1)=1
     ELSE
        itop(jl,1)=0
     ENDIF
230 END DO
  !
  DO 232 jk=ktdia+1,klev
     DO 231 jl=1,kproma
        IF (paclc(jl,jk) .GE. zepsec .AND. itop(jl,jk-1) .EQ. 0) THEN
           itop(jl,jk)=jk
        ELSE IF (paclc(jl,jk) .GE. zepsec) THEN
           itop(jl,jk)=itop(jl,jk-1)
        ELSE
           itop(jl,jk)=0
        ENDIF
231  END DO
232 END DO
  !
  !         1.5 CLOUD BASES
  !
  DO 233 jl=1,kproma
     IF (paclc(jl,klev) .GE. zepsec) THEN
        ibas(jl,klev)=klev
     ELSE
        ibas(jl,klev)=0
     ENDIF
233 END DO
!
  DO 235 jk=klev-1,ktdia,-1
     DO 234 jl=1,kproma
        IF (paclc(jl,jk) .GE. zepsec .AND. ibas(jl,jk+1) .EQ. 0) THEN
           ibas(jl,jk)=jk
        ELSE IF (paclc(jl,jk) .GE. zepsec) THEN
           ibas(jl,jk)=ibas(jl,jk+1)
        ELSE
           ibas(jl,jk)=0
        ENDIF
234  END DO
235 END DO

  DO jk=klev,ktdia,-1
     DO jl=1,kproma
        !
        it       = NINT(ptm1(jl,jk)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zesw     = tlucuaw(it)/papm1(jl,jk)
        zesw     = MIN(zesw,0.5_dp)
        zqsw     = zesw/(1._dp-vtmpc1*zesw)
        zsusatw_2d(jl,jk)=MAX(pqm1(jl,jk)/zqsw-1.0_dp,0.0_dp)
        !--- Store supersaturation with respect to water in stream:
        swat(jl,jk,jrow)=zsusatw_2d(jl,jk)
        !--- Saturation water vapour pressure:
        !zesw_2d(jl,jk)=zesw*papm1(jl,jk)*rv/rd

     END DO !jl
  END DO !jk

  DO jk=klev,ktdia,-1
     DO jl=1,kproma
        IF (ibas(jl,jk) .EQ. jk) THEN
           IF (zcdnc(jl,jk).LE.zcdnmin .OR. pxlte(jl,jk).GT.0._dp &
               .OR. paclc(jl,jk).GT.cloud_tm1(jl,jk,jrow) &
               .OR. zsusatw_2d(jl,jk).GT.zeps) THEN
             zqlnuc(jl,jk)=MAX(0._dp,pcdncact(jl,jk))
             zcdnc(jl,jk)=zcdnc(jl,jk)+zqlnuc(jl,jk)
             qnuc(jl,jk,jrow)=zqlnuc(jl,jk)/ztmst
           ENDIF
        ENDIF
        !
        IF (jk < klev) THEN
           IF (itop(jl,jk) .NE. ibas(jl,jk) .AND. ibas(jl,jk) .GT. jk &
                .AND. zqlnuc(jl,jk+1)>0._dp) THEN
              zqlnuc(jl,jk)=zqlnuc(jl,jk+1)
              qnuc(jl,jk,jrow)=(zcdnc(jl,jk+1)-zcdnc(jl,jk))/ztmst
              zcdnc(jl,jk)=zcdnc(jl,jk+1)
           ENDIF
        ENDIF

        iclbas(jl)=NINT(pcvcbot(jl))
        zqlnuccvh(jl,jk)=0._dp
        IF (jk.EQ.iclbas(jl)) THEN
           zvervc= (-pvervel(jl,jk)/(zrho(jl,jk)*g) + pwcape(jl) &
                + 1.33_dp*SQRT(ptkem1(jl,jk)))
           zvervc=MAX(zvervc,1._dp)
           zcdncnew=zna(jl,jk)*zvervc/(zvervc+2.3E-10_dp*zna(jl,jk))
           zcdncnew=1.e5_dp*(1.0E-6_dp*zcdncnew)**1.27_dp
           zqlnuccvh(jl,jk)=MAX(0._dp,MIN(zna(jl,jk),zcdncnew))
        ENDIF

     END DO
  END DO

     DO jk=klev,ktdia,-1
        DO jl=1,kproma
           IF (ncvmicro>0) THEN
              IF (pxtecl(jl,jk).GT.0._dp.AND.iclbas(jl).GT.0._dp) THEN
                 zqlnuccv=MAX(0._dp,MIN(zna(jl,jk),zqlnuccvh(jl,iclbas(jl))))
                 qnuc(jl,jk,jrow)=zqlnuccv/ztmst
                 zcdnc(jl,jk)=zcdnc(jl,jk)+zqlnuccv
              ENDIF
           ELSE
              IF (pxtec(jl,jk).GT.0._dp.AND.iclbas(jl).GT.0 .AND. &
                   ptm1(jl,jk).GT.cthomi) THEN
                 zqlnuccv=MAX(0._dp,MIN(zna(jl,jk),zqlnuccvh(jl,iclbas(jl))))
                 qnuc(jl,jk,jrow)=zqlnuccv/ztmst
                 zcdnc(jl,jk)=zcdnc(jl,jk)+zqlnuccv
              ENDIF
           ENDIF
        END DO
     END DO
  
   DO jk=klev,ktdia,-1
      DO jl=1,kproma
           !     Ice nucleation
           !
        it= NINT(ptm1(jl,jk)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zesi     = tlucua(it)/papm1(jl,jk)
        zesi     = MIN(zesi,0.5_dp)
        zqsi     = zesi/(1._dp-vtmpc1*zesi)
        zsusati  = MAX(pqm1(jl,jk)/zqsi-1.0_dp,0.0_dp)
        sice(jl,jk,jrow)=zsusati
        ztc=MIN(0._dp,ptm1(jl,jk)-tmelt)

!       test using Faisal's param
!
        zrieff=23.2_dp*EXP(0.015_dp*ztc)
        zrieff=MAX(zrieff,1.0_dp)
        zrih=-2261.236_dp+SQRT(5113188.044_dp+2809._dp*zrieff*zrieff*zrieff)
        zrid=1.e-6_dp*zrih**(1._dp/3._dp)
        zrid=MAX(1.e-6_dp,zrid)

        znicenew=0._dp
        IF (ncvmicro > 0) THEN
           IF (pxteci(jl,jk) .GT. 0._dp) THEN
              znicenew=0.75_dp*zrho(jl,jk)*pxteci(jl,jk)*ztmst/(api*zrid**3*zrhoice)
              zicnc(jl,jk)=zicnc(jl,jk)+znicenew
           ENDIF
        ELSE
           IF (pxtec(jl,jk) .GT. 0._dp .AND. ptm1(jl,jk) .LT. cthomi) THEN
              znicenew=0.75_dp*zrho(jl,jk)*pxtec(jl,jk)*ztmst/(api*zrid**3*zrhoice)
              zicnc(jl,jk)=zicnc(jl,jk)+znicenew
           ENDIF
        ENDIF
        IF ( nicnc .EQ. 1 ) THEN 
           ! Use ICNC scheme after Lohmann, JAS, 2002
           !
           !
           IF (zsusati .GT. 0._dp .AND. ptm1(jl,jk) .LT. cthomi) THEN
              znicenew1=0.75_dp*zrho(jl,jk)*zsusati*zqsi/(api*zrid**3*zrhoice)
              zninucl(jl,jk)=MAX(MIN(znicenew1-zicnc(jl,jk),zascs(jl,jk)-zicnc(jl,jk)),0._dp)
              zicnc(jl,jk)=zicnc(jl,jk)+zninucl(jl,jk)
           ENDIF
        ENDIF
        !
     END DO ! jl
  END DO ! jk

  IF (nicnc > 1) THEN

    DO jk = klev, ktdia, -1
       DO jl = 1, kproma
     ! use Bernd's cirrus scheme (Kaercher and Lohmann, JGR, 2002b; Lohmann et al., JGR, 2004)
          it= NINT(ptm1(jl,jk)*1000._dp)
          IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
          it = MAX(MIN(it,jptlucu2),jptlucu1)
          zesi     = tlucua(it)/papm1(jl,jk)
          zesi     = MIN(zesi,0.5_dp)
          zqsi     = zesi/(1._dp-vtmpc1*zesi)
          zsusati  = MAX(pqm1(jl,jk)/zqsi-1.0_dp,0.0_dp)

          IF (jk.NE.klev) THEN
             zvervx(jl) =100.0_dp*(-pvervel(jl,jk)/(zrho(jl,jk)*g) &
             +0.7_dp*SQRT(ptkem1(jl,jk)))
          ELSE
             zvervx(jl) =100.0_dp*(-pvervel(jl,jk))/(zrho(jl,jk)*g)
          END IF

          zdp(jl) = paphm1(jl,jk+1)-paphm1(jl,jk)
          zsusatix(jl)  = zsusati
          znicex(jl)    = zicemin

          zaprx(jl,1)   = MAX(wetradius(jl,jk,iaits)*100.0_dp,crdiv(2))![cm]
          zaprx(jl,2)   = MAX(wetradius(jl,jk,iaccs)*100.0_dp,crdiv(3))![cm]
          zaprx(jl,3)   = MAX(wetradius(jl,jk,icoas)*100.0_dp,crdiv(4))![cm]

          zapsigx(jl,1) = sigma(2)
          zapsigx(jl,2) = sigma(3)
          zapsigx(jl,3) = sigma(4)

          zapnx(jl,1)   = MAX((pxtm1(jl,jk,idt_nks)+pxtte(jl,jk,idt_nks) &
               *ztmst+pxtm1(jl,jk,idt_nki)+pxtte(jl,jk,idt_nki) &
               *ztmst)/M_air * 1000._dp * zrho(jl,jk) * 1.0e-6_dp,5._dp)![1/cm3]
          zapnx(jl,2)   = MAX((pxtm1(jl,jk,idt_nas)+pxtte(jl,jk,idt_nas) &
               *ztmst+pxtm1(jl,jk,idt_nai)+pxtte(jl,jk,idt_nai) &
               *ztmst)/M_air * 1000._dp * zrho(jl,jk) * 1.0e-6_dp,5._dp)![1/cm3]
          zapnx(jl,3)   = MAX((pxtm1(jl,jk,idt_ncs)+pxtte(jl,jk,idt_ncs) &
               *ztmst+pxtm1(jl,jk,idt_nci)+pxtte(jl,jk,idt_nci) &
               *ztmst)/M_air * 1000._dp * zrho(jl,jk) * 1.0e-6_dp,0.01_dp)![1/cm3]

       END DO ! jl-loop

       CALL xfrzmstr(lhet,nosize,ztmst,                   &
                     klev,kbdim,kproma,jk,nfrzmod,        &
                     zsusatix,zvervx,zapnx,               &
                     zaprx,zapsigx,ptm1,tmelt,zeps,papm1, &
                     cthomi,zri,znicex)

! Update ICNC

       DO jl = 1, kproma
          zri(jl,jk)=MAX(zri(jl,jk),1.e-6_dp)
          it= NINT(ptm1(jl,jk)*1000._dp)
          IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
          it = MAX(MIN(it,jptlucu2),jptlucu1)
          zqsi     = tlucua(it)/papm1(jl,jk)
          zesi     = zqsi*papm1(jl,jk)*rv/rd
          zqsi     = MIN(zqsi,0.5_dp)
          zqsi     = zqsi/(1._dp-vtmpc1*zqsi)
          zsusati=zsusatix(jl)
          zverv=zvervx(jl)
          IF (zsusati.GT.0.0_dp.AND.ptm1(jl,jk).LT.cthomi.AND.zsusati.GT.0.0_dp) THEN
             zninucl(jl,jk) = MAX(0.0_dp,MIN((zapnx(jl,1)+zapnx(jl,2)+  &
                    zapnx(jl,3))*1.e6_dp-zicnc(jl,jk),znicex(jl)-zicnc(jl,jk)))
             zicnc(jl,jk)=zicnc(jl,jk)+zninucl(jl,jk)
             zdv=2.21_dp/papm1(jl,jk)
             zast=als*(als/(rv*ptm1(jl,jk))-1._dp)/(zka*ptm1(jl,jk))
             zbst=rv*ptm1(jl,jk)/(zdv*zesi)
             zgtp=1._dp/(zrho(jl,jk)*(zast+zbst))
             zvth    = SQRT( 8._dp*zkb*ptm1(jl,jk) / (api*zxmw) )
             zb2= 0.25_dp * zalpha * zvth   / zdv
             zfuchs=1._dp/(1._dp+zb2*zri(jl,jk))
             zfre=1._dp
             IF (pxim1(jl,jk) > 0._dp) THEN
                zviscos=(1.512_dp + 0.0052_dp*(ptm1(jl,jk)-233.15_dp))*1.e-5_dp
                zxifall = cvtfall*(zrho(jl,jk)*pxim1(jl,jk))**0.16_dp
                zre=zrho(jl,jk)*2._dp*zri(jl,jk)*zxifall/zviscos
                zfre=1._dp+ 0.229_dp*SQRT(zre)
             ENDIF

             zdepos=MAX(-pxim1(jl,jk),MIN(4._dp*api*zri(jl,jk)*zsusati*zicnc(jl,jk) &
                  *zfre*zgtp*zfuchs*ztmst,pqm1(jl,jk)-zqsi))
             zqinucl(jl,jk)=zdepos
          ENDIF
       ENDDO !jl
    ENDDO ! jk
  ENDIF ! which nucleation scheme (nicnc)

  !corinna: set zcdnc and zicnc to minium now if nucleation is not strong enough
  DO jk = ktdia, klev
     DO jl = 1, kproma
        IF (paclc(jl,jk).GT.0._dp .AND. ptm1(jl,jk).GT.cthomi) zcdnc(jl,jk)=MAX(zcdnc(jl,jk),zcdnmin)
        IF (paclc(jl,jk).GT.0._dp .AND. ptm1(jl,jk).LT.tmelt .AND. nicnc > 0) &
      zicnc(jl,jk)=MAX(zicnc(jl,jk),zicemin)
     ENDDO !jl
  ENDDO ! jk

!--- End included for CDNC/IC scheme -----------------------------------

  DO 831 jk=ktdia,klev
!
!     ------------------------------------------------------------------
!
!       2.    Set to zero local tendencies (increments)
!
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

     DO 201 jl = 1,kproma
        zcnd(jl)       = 0.0_dp
        zdep(jl)       = 0.0_dp
        zfrl(jl)       = 0.0_dp
        zspr(jl)       = 0.0_dp
        zimlt(jl)      = 0.0_dp
        zximlt(jl)     = 0.0_dp
        zxisub(jl)     = 0.0_dp
        zsmlt(jl)      = 0.0_dp
        zsacl(jl)      = 0.0_dp
        zgenti(jl)     = 0.0_dp
        zgentl(jl)     = 0.0_dp
        zxievap(jl)    = 0.0_dp
        zxlevap(jl)    = 0.0_dp
        zvartg(jl)     = 0.0_dp
        zconvvar(jl)   = 0.0_dp
        zconvskew(jl)  = 0.0_dp
        zturbvar(jl)   = 0.0_dp
        zturbskew(jl)  = 0.0_dp
        zmicroskew(jl) = 0.0_dp
        zxvarte(jl)    = 0.0_dp
        zxskewte(jl)   = 0.0_dp  
        zdp(jl)        = paphm1(jl,jk+1)-paphm1(jl,jk)
        zdz(jl)        = (zgeoh(jl,jk)-zgeoh(jl,jk+1))/g
        zrcp           = 1._dp/(cpd+zcons1*MAX(pqm1(jl,jk),0.0_dp))
        zlvdcp(jl)     = alv*zrcp
        zlsdcp(jl)     = als*zrcp
!--- Included for prognostic CDNC/IC scheme ----------------------------
!    additional arrays for transformation rate changes of CDNC and ICNC
        zrprn(jl)=0._dp
        zfrln(jl)=0._dp
        zsacln(jl)=0._dp
        zsprn(jl,jk)=0._dp
!--- End included for CDNC/IC scheme -----------------------------------
201  END DO
!
!     ------------------------------------------------------------------
!
!       3.   Modification of incoming precipitation fluxes by
!            melting, sublimation and evaporation
!
     IF (jk .GT. 1) THEN
!
!DIR$ CONCURRENT
        DO 331 jl = 1,kproma
!
!       3.1   Melting of snow and ice
!
           zcons     = zcons2*zdp(jl)/(zlsdcp(jl)-zlvdcp(jl))
           ztdif     = MAX(0.0_dp,ptm1(jl,jk)-tmelt)
           zsnmlt    = MIN(zxsec*zsfl(jl),zcons*ztdif)
           zrfl(jl)  = zrfl(jl)+zsnmlt
           zsfl(jl)  = zsfl(jl)-zsnmlt
           zsmlt(jl) = zsnmlt/(zcons2*zdp(jl))
           zximelt   = MIN(zxsec*zxiflux(jl),zcons*ztdif)
           pimelt(jl,jk) = zximelt ! mz_ht_20071114
           zxiflux(jl)=zxiflux(jl)-zximelt
           zximlt(jl) =zximelt/(zcons2*zdp(jl))
           IF (ztdif.GT.0.0_dp) THEN
            zimlt(jl) = MAX(0.0_dp,pxim1(jl,jk)+pxite(jl,jk)*ztmst)
!--- Included for prognostic CDNC/IC scheme (Philip Stier, 31/03/2004) -
!    If T > tmelt melt all ice crystals and transfer to cloud droplets
            IF (nicnc > 0) THEN
               zicncmelt=zicnc(jl,jk)
               zicnc(jl,jk)=cqtmin
            ELSE
               IF (paclc(jl,jk) .GE. zepsec) THEN
                  zicncmelt=(zimlt(jl)*zrho(jl,jk)/zmi)/paclc(jl,jk)
               ELSE
                  zicncmelt=(zimlt(jl)*zrho(jl,jk)/zmi)
               END IF
            ENDIF
            zcdnc(jl,jk)=zcdnc(jl,jk)+zicncmelt
            qmel(jl,jk,jrow)=zicncmelt/ztmst
!--- End included for CDNC/IC scheme -----------------------------------
           ELSE
            zimlt(jl) = 0.0_dp
           END IF
!
!       3.2   Sublimation of snow and ice (Lin et al., 1983)
!
           IF (zclcpre(jl) .GT. 0.0_dp) THEN
              zclcstar = zclcpre(jl)
              zdpg     = zdp(jl)/g
              zqrho    = 1.3_dp/zrho(jl,jk)
              it       = NINT(ptm1(jl,jk)*1000._dp)
              IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
              it = MAX(MIN(it,jptlucu2),jptlucu1)
              zesi     = tlucua(it)/papm1(jl,jk)
              zesi     = MIN(zesi,0.5_dp)
              zqsi     = zesi/(1._dp-vtmpc1*zesi)
              zsusati  = MIN(pqm1(jl,jk)/zqsi-1.0_dp,0.0_dp)
              zb1      = zlsdcp(jl)**2/(2.43e-2_dp*rv*(ptm1(jl,jk)**2))
              zb2      = 1._dp/(zrho(jl,jk)*zqsi*0.211e-4_dp)
              zcoeff   = 3.e6_dp*2._dp*api*zsusati/(zrho(jl,jk)*(zb1+zb2))
!
              IF (zsfl(jl) .GT. cqtmin) THEN
               zxsp1    = (zsfl(jl)/(zclcpre(jl)*cvtfall))**(1._dp/1.16_dp)
               zclambs  = (zxsp1/(api*crhosno*cn0s))**0.25_dp
               zcfac4c  = 0.78_dp*zclambs**2+232.19_dp*zqrho**0.25_dp         &
                                                    *zclambs**2.625_dp
               zzeps    = MAX(-zxsec*zsfl(jl)/zclcpre(jl),             &
                                                  zcoeff*zcfac4c*zdpg)
               zsub(jl) = -zzeps/zdpg*ztmst*zclcstar
               zsub(jl) = MIN(zsub(jl),                                &
                                    MAX(zxsec*(zqsi-pqm1(jl,jk)),0.0_dp))
               zsub(jl) = MAX(zsub(jl),0.0_dp)
              END IF
!
              IF (zxiflux(jl) .GT. cqtmin) THEN
               zxsp1    = (zxiflux(jl)/(zclcpre(jl)*cvtfall))**(1._dp/1.16_dp)
               zclambs  = (zxsp1/(api*crhosno*cn0s))**0.25_dp
               zcfac4c  = 0.78_dp*zclambs**2+232.19_dp*zqrho**0.25_dp           &
                                                    *zclambs**2.625_dp
               zzeps    = MAX(-zxsec*zxiflux(jl)/zclcpre(jl),          &
                                                  zcoeff*zcfac4c*zdpg)
               zsubi    = -zzeps/zdpg*ztmst*zclcstar
               zsubi    = MIN(zsubi,MAX(zxsec*(zqsi-pqm1(jl,jk)),0.0_dp))
               zsubi    = MAX(zsubi,0.0_dp)
               zxiflux(jl)=zxiflux(jl)-zsubi*zcons2*zdp(jl)
               zxisub(jl) =zsubi
              END IF
           END IF
!
!       3.3   Evaporation of rain (Rotstayn, 1997)
!
           IF (zclcpre(jl) .GT. 0.0_dp .AND. zrfl(jl) .GT. cqtmin) THEN
              zclcstar = zclcpre(jl)
              zdpg     = zdp(jl)/g
              zqrho    = 1.3_dp/zrho(jl,jk)
              it       = NINT(ptm1(jl,jk)*1000._dp)
              IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
              it = MAX(MIN(it,jptlucu2),jptlucu1)
              zesw     = tlucuaw(it)/papm1(jl,jk)
              zesat    = zesw*papm1(jl,jk)*rv/rd
              zesw     = MIN(zesw,0.5_dp)
              zqsw     = zesw/(1._dp-vtmpc1*zesw)
              zsusatw  = MIN(pqm1(jl,jk)/zqsw-1.0_dp,0.0_dp)
              zdv      = 2.21_dp/papm1(jl,jk)
              zast     = alv*(alv/(rv*ptm1(jl,jk))-1.0_dp)/               &
                                            (0.024_dp*ptm1(jl,jk))
              zbst     = rv*ptm1(jl,jk)/(zdv*zesat)
              zzepr    = 870._dp*zsusatw*(zrfl(jl)/zclcpre(jl))**0.61_dp     &
                                     /(SQRT(zrho(jl,jk))*(zast+zbst))
              zzepr    = MAX(-zxsec*zrfl(jl)/zclcpre(jl),zzepr*zdpg)
              zevp(jl) = -zzepr/zdpg*ztmst*zclcstar
              zevp(jl) = MIN(zevp(jl),MAX(zxsec*(zqsw-pqm1(jl,jk)),0.0_dp))
              zevp(jl) = MAX(zevp(jl),0.0_dp)
           END IF
331     END DO
!
!        IF (lookupoverflow) CALL lookuperror ('cloud (1)    ')
         IF (lookupoverflow) THEN
           status_string = 'lookuperror: cdnc - icnc - cloud (1)'
           RETURN
         ENDIF
!
     END IF
!
!DIR$ CONCURRENT
     DO 610 jl=1,kproma
!
!     ------------------------------------------------------------------
!       4.    Sedimentation of cloud ice from grid-mean values.
!             Updating the tendency 'pxite' to include sedimentation.
!             At jk=klev, the sedimentation sink is balanced by
!             precipitation at the surface (through 'zzdrs', see 7.3).
!             Finally: In-cloud cloud water/ice.
!
        zxip1         = pxim1(jl,jk)+pxite(jl,jk)*ztmst-zimlt(jl)
        zxip1         = MAX(zxip1,EPSILON(1._dp))
        zxifall       = cvtfall*(zrho(jl,jk)*zxip1)**0.16_dp
        zicnc(jl,jk)  = MAX(zicemin,zicnc(jl,jk))
        zmmean(jk)=MAX((zrho(jl,jk)*zxip1)/zicnc(jl,jk)*1.e3_dp,zmi*1000._dp)
        zmm(jk)=MAX((zrho(jl,jk)*zxip1)/zicnc(jl,jk)/EXP(0.5_dp*(zsigma_m2))*1.e3_dp,zmi*1000._dp)
        IF (zmm(jk) .LT. 2.166E-6_dp) THEN 
           zalfased=121107.08_dp
           zbetased=0.5727_dp
        ELSE
           zalfased=3898.65_dp
           zbetased=0.3091_dp
        ENDIF
        zaaa=( (papm1(jl,jk)/30000._dp)**(-0.178_dp) )* ( (ptm1(jl,jk)/233.0_dp)**(-0.394_dp) )
        zxifallmc(jk)=zfall*zalfased*(zmm(jk)**zbetased)* &
             EXP(0.5_dp*(zbetased**2+2._dp*zbetased)*zsigma_m2)*zaaa/100._dp
        zxifallnc(jk)=zfall*zxifallmc(jk)/EXP(zbetased*zsigma_m2)
!S.K. changed to
        IF (jk .GT. 1 ) THEN
           IF (zmmean(jk-1) .GT. 0._dp) THEN
              zxifluxn(jk)=MAX(0.0_dp,zxiflux(jl)/(zmmean(jk-1)/1000.0_dp*EXP(zbetased*zsigma_m2)))
           ELSE
              zxifluxn(jk)=0._dp
           ENDIF
        ELSE
           zxifluxn(jk)=0._dp
        ENDIF

        zxifallmc(jk) = zxifall
        zal1          = zxifallmc(jk)*g*zrho(jl,jk)*ztmst/zdp(jl)
        zal2          = 0._dp
        IF (zxifallmc(jk)>0._dp) zal2 = zxiflux(jl)/(zrho(jl,jk)*zxifallmc(jk))
        zxised(jk)    =  zxip1*EXP(-zal1)+zal2*(1._dp-EXP(-zal1))
!S.K. changed to :
        IF(zxised(jk).LE.zxip1.AND.zxip1.GT.0._dp) THEN
           zicesed(jk)=MAX(zicemin,zicnc(jl,jk)*(zxised(jk)/zxip1))  
        ELSE 
           zicesed(jk)=0._dp
        ENDIF
!END S.K. 
        zicnc(jl,jk) = zicesed(jk)
        pisedi(jl,jk) = zxised(jk)
        zxiflux(jl)  = zxiflux(jl)+(zxip1-zxised(jk))*zcons2*zdp(jl)
        pxite(jl,jk) = (zxised(jk)-pxim1(jl,jk))/ztmst

!
!             In-cloud water/ice calculated from respective grid-means,
!             partial cloud cover, advective/diffusive tendencies,
!             detrained cloud water/ice and ice sedimentation.
!             In-cloud values are required for cloud microphysics.
!
        zclcaux(jl) = paclc(jl,jk)
        locc        = zclcaux(jl) .GT. 0.0_dp
        lo2         = (ptm1(jl,jk) .LT. cthomi) .OR.                   &
                      (ptm1(jl,jk) .LT. tmelt .AND. zxised(jk) .GT. csecfrl &
                      .AND. zsusatw_2d(jl,jk) .LT. 1._dp)
        IF (lo2) THEN                 !     ice cloud
!---------------Added by Junhua Zhang for CONV Micro----------------------
           IF(ncvmicro>0) THEN
            zxite(jl)          = pxteci(jl,jk)
            zxlte(jl)          = pxtecl(jl,jk)
           ELSE
            zxite(jl)          = pxtec(jl,jk)
            zxlte(jl)          = 0.0_dp
           ENDIF
!----------------------------------end-------------------------------------
           IF (locc) THEN
              zxib(jl)        = pxim1(jl,jk)/zclcaux(jl)
              zxlb(jl)        = pxlm1(jl,jk)/zclcaux(jl)
              zxim1evp        = 0.0_dp
              zxlm1evp        = 0.0_dp
              zxidt           = (pxite(jl,jk)+zxite(jl))*ztmst
              zxldt           =  pxlte(jl,jk)*ztmst+zximlt(jl)+zimlt(jl)
              IF (zxidt .GT. 0.0_dp) THEN
                 zxidtstar    = zxidt
                 zxib(jl)     = zxib(jl)+zxidt
              ELSE
                 zxidtstar    = 0.0_dp
                 zxib(jl)     = zxib(jl)+MAX(zxidt/zclcaux(jl),        &
                                                       -zxib(jl))
                 pxite(jl,jk) = MAX(pxite(jl,jk),                      &
                                 -(pxim1(jl,jk)/ztmst+zxite(jl)))
              END IF
              IF (zxldt .GT. 0.0_dp) THEN
                 zxldtstar    = zxldt
                 zxlb(jl)     = zxlb(jl)+zxldt
              ELSE
                 zxldtstar    = 0.0_dp
                 zxlb(jl)     = zxlb(jl)+MAX(zxldt/zclcaux(jl),        &
                                                       -zxlb(jl))
                 pxlte(jl,jk) = MAX(pxlte(jl,jk),-pxlm1(jl,jk)/ztmst)
              END IF
           ELSE                      !    cloud cover = 0._dp
              zxib(jl)        = 0.0_dp
              zxlb(jl)        = 0.0_dp
              zxidt           = (pxite(jl,jk)+zxite(jl))*ztmst
              zxldt           =  pxlte(jl,jk)*ztmst+zximlt(jl)+zimlt(jl)
              IF (zxidt .GT. 0.0_dp) THEN
                 zxidtstar    = zxidt
                 zxim1evp     = pxim1(jl,jk)
              ELSE
                 zxidtstar    = 0.0_dp
                 pxite(jl,jk) = MAX(pxite(jl,jk),                      &
                                  -(pxim1(jl,jk)/ztmst+zxite(jl)))
                 zxim1evp     = pxim1(jl,jk)+(pxite(jl,jk)+zxite(jl))  &
                                                                 *ztmst
              END IF
              IF (zxldt .GT. 0.0_dp) THEN
                 zxldtstar    = zxldt
                 zxlm1evp     = pxlm1(jl,jk)
              ELSE
                 zxldtstar    = 0.0_dp
                 pxlte(jl,jk) = MAX(pxlte(jl,jk),-pxlm1(jl,jk)/ztmst)
                 zxlm1evp     = pxlm1(jl,jk)+pxlte(jl,jk)*ztmst
              END IF
           END IF
        ELSE                           !    water cloud
!---------------Added by Junhua Zhang for CONV Micro----------------------
           IF(ncvmicro>0) THEN
            zxlte(jl)          = pxtecl(jl,jk)
            zxite(jl)          = pxteci(jl,jk)
           ELSE
            zxlte(jl)          = pxtec(jl,jk)
            zxite(jl)          = 0.0_dp
           ENDIF
!----------------------------------end-------------------------------------
           IF (locc) THEN
              zxlb(jl)        = pxlm1(jl,jk)/zclcaux(jl)
              zxib(jl)        = pxim1(jl,jk)/zclcaux(jl)
              zxlm1evp        = 0.0_dp
              zxim1evp        = 0.0_dp
              zxldt           = (pxlte(jl,jk)+zxlte(jl))*ztmst         &
                                             +zximlt(jl)+zimlt(jl)
              zxidt           =  pxite(jl,jk)*ztmst
              IF (zxldt .GT. 0.0_dp) THEN
                 zxldtstar    = zxldt
                 zxlb(jl)     = zxlb(jl)+zxldt
              ELSE
                 zxldtstar    = 0.0_dp
                 zxlb(jl)     = zxlb(jl)+MAX(zxldt/zclcaux(jl),        &
                                                          -zxlb(jl))
                 pxlte(jl,jk) = MAX(pxlte(jl,jk),                      &
                                 -(pxlm1(jl,jk)/ztmst+zxlte(jl)))
              END IF
              IF (zxidt .GT. 0.0_dp) THEN
                 zxidtstar    = zxidt
                 zxib(jl)     = zxib(jl)+zxidt
              ELSE
                 zxidtstar    = 0.0_dp
                 zxib(jl)     = zxib(jl)+MAX(zxidt/zclcaux(jl),        &
                                                          -zxib(jl))
                 pxite(jl,jk) = MAX(pxite(jl,jk),-pxim1(jl,jk)/ztmst)
              END IF
           ELSE                          !    cloud cover = 0._dp
              zxlb(jl)        = 0.0_dp
              zxib(jl)        = 0.0_dp
              zxldt           = (pxlte(jl,jk)+zxlte(jl))*ztmst         &
                                             +zximlt(jl)+zimlt(jl)
              zxidt           =  pxite(jl,jk)*ztmst
              IF (zxldt .GT. 0.0_dp) THEN
                 zxldtstar    = zxldt
                 zxlm1evp     = pxlm1(jl,jk)
              ELSE
                 zxldtstar    = 0.0_dp
                 pxlte(jl,jk) = MAX(pxlte(jl,jk),                      &
                                  -(pxlm1(jl,jk)/ztmst+zxlte(jl)))
                 zxlm1evp     = pxlm1(jl,jk)+(pxlte(jl,jk)+zxlte(jl))  &
                                                                *ztmst
              END IF
              IF (zxidt .GT. 0.0_dp) THEN
                 zxidtstar    = zxidt
                 zxim1evp     = pxim1(jl,jk)
              ELSE
                 zxidtstar    = 0.0_dp
                 pxite(jl,jk) = MAX(pxite(jl,jk),-pxim1(jl,jk)/ztmst)
                 zxim1evp     = pxim1(jl,jk)+pxite(jl,jk)*ztmst
              END IF
           END IF
        END IF
!
!     ------------------------------------------------------------------
!       5.    Condensation/deposition and evaporation/sublimation
!
!             zlc       =  L_{v/s} / c_p
!             zlcdqsdt  = L dq_sat / c_p dT
!             zdqsdt    = dq_sat / dT
!
        zrcp        = 1._dp/(cpd+zcons1*MAX(pqm1(jl,jk),0.0_dp))
        zlvdcp(jl)  = alv*zrcp
        zlsdcp(jl)  = als*zrcp
        zlc         = MERGE(zlsdcp(jl),zlvdcp(jl),lo2)
        it          = NINT(ptm1(jl,jk)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zqsm1       = MERGE(tlucua(it),tlucuaw(it),lo2)/papm1(jl,jk)
        zqsm1       = MIN(zqsm1,0.5_dp)
        zqsm1       = zqsm1/(1._dp-vtmpc1*zqsm1)
        it1         = it+1
        it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
        zqst1       = MERGE(tlucua(it1),tlucuaw(it1),lo2)/papm1(jl,jk)
        zqst1       = MIN(zqst1,0.5_dp)
        zqst1       = zqst1/(1._dp-vtmpc1*zqst1)
        zdqsdt      = (zqst1-zqsm1)*1000._dp
        zlcdqsdt    = zlc*zdqsdt
        zxievap(jl) = (1.0_dp-zclcaux(jl))*zxidtstar+zxim1evp
        zxlevap(jl) = (1.0_dp-zclcaux(jl))*zxldtstar+zxlm1evp
        zqvdt       = pqte(jl,jk)*ztmst+zevp(jl)+zsub(jl)              &
                             +zxievap(jl)+zxlevap(jl)+zxisub(jl)
        zdtdt       = ptte(jl,jk)*ztmst-zlvdcp(jl)*(zevp(jl)           &
                       +zxlevap(jl)) -(zlsdcp(jl)-zlvdcp(jl))          &
                       *(zsmlt(jl)+zximlt(jl)+zimlt(jl))
        zdtdt    = zdtdt-zlsdcp(jl)*(zsub(jl)+zxievap(jl)+zxisub(jl))

        zqp1        = MAX(pqm1(jl,jk)+zqvdt,0.0_dp)
        ztp1        = ptm1(jl,jk)+zdtdt
        zdtdtstar   = zdtdt+zclcaux(jl)*(zlc*pqte(jl,jk)*ztmst         &
                           +zlvdcp(jl)*(zevp(jl)+zxlevap(jl))          &
                         +zlsdcp(jl)*(zsub(jl)+zxievap(jl)+zxisub(jl)))
        zdqsat      = zdtdtstar*zdqsdt/(1._dp+zclcaux(jl)*zlcdqsdt)
        zxib(jl)    = MAX(zxib(jl),0.0_dp)
        zxlb(jl)    = MAX(zxlb(jl),0.0_dp)
        zxilb       = zxib(jl)+zxlb(jl)
!
!       Diagnostics: relative humidity
!
        zrelhum        = pqm1(jl,jk)/zqsm1
        zrelhum        = MAX(MIN(zrelhum,1._dp),0._dp)
        prelhum(jl,jk) = zrelhum
!
        IF (jk.GE.ncctop) THEN
!
!       define variables needed for cover scheme
!
!       zbetaqt = total water
!       zbetass = saturation mixing ratio adjusted to match qv
!       zwide   = current diagnosed distribution width
!
           zbetacl(jl) = MAX(0.0_dp,pxlm1(jl,jk))+MAX(0.0_dp,pxim1(jl,jk))
           zbetaqt(jl) = MAX(cqtmin,pqm1(jl,jk))+zbetacl(jl)
           zvartg(jl)  = MAX(cqtmin,cvarmin*pqm1(jl,jk))
           zwide(jl)   = MAX(zvartg(jl),pbetab(jl,jk)-pbetaa(jl,jk))
           zskew       = MAX(MIN(pxskew(jl,jk),cbeta_pq_max),cbeta_pq)
           iqidx       = INT((nbetaq/cbetaqs)*                         &
                                 LOG((zskew-cbeta_pq)/rbetak+1._dp)+0.5_dp)
!
!
!       5.1 Turbulence: Skewness - equation solved implicitly
!           This solver only works if phmixtau has non-zero timescale
!
           zqtau         = phmixtau(jl,jk)+pvmixtau(jl,jk)
           zbqp1         = cbeta_pq-(cbeta_pq-pxskew(jl,jk))           &
                                                    *EXP(-zqtau*zdtime)
           zbqp1         = MAX(MIN(zbqp1,cbeta_pq_max),cbeta_pq)
           zturbskew(jl) = (zbqp1-pxskew(jl,jk))/zdtime
!
!       5.2 Turbulence: variance - equation solved implicitly
!
           zpp          = cbeta_pq
           zqq          = pxskew(jl,jk)
           zeta         = (zpp+zqq)**2*(zpp+zqq+1._dp)/(zpp*zqq)
           zprod        = zeta*pvdiffp(jl,jk)/zwide(jl)
           zbbap1       = zprod/zqtau+zvartg(jl)-(zprod/zqtau          &
                            +zvartg(jl)-zwide(jl))*EXP(-zqtau*zdtime)
           zbbap1       = MAX(zbbap1,zvartg(jl))
           zbbap1       = MIN(zbbap1,zbetaqt(jl)*(cbeta_pq+zbqp1)      &
                                                            /cbeta_pq)
           zturbvar(jl) = (zbbap1-zwide(jl))/zdtime
           zbap1        = zbetaqt(jl)-zbbap1*cbeta_pq/(cbeta_pq+zbqp1)
!
           IF (lcover) THEN ! mz_jd_20161011
!              translated into apparent xl,xi,q and heat sources
!              first order effect only, effect of evaporation of
!              cloud on qsat taken into account in thermodynamic budget
!              but does not change the mixing term here since that
!              would require iteration and is therefore neglected
!
!              calculate values after one timestep
!
           iqidx   = INT((nbetaq/cbetaqs)*                             &
                     LOG((zbqp1-cbeta_pq)/rbetak+1._dp)+0.5_dp)
           ztt     = (pbetass(jl,jk)-zbap1)*cbeta_pq/                     &
                     ((zbetaqt(jl)-zbap1)*(cbeta_pq+zbqp1))
           ztt     = nbetax*MAX(MIN(ztt,1.0_dp),0.0_dp)
           ixidx   = INT(ztt)
           IF (ixidx == nbetax) THEN
              zbetai0 = 1.0_dp
              zbetai1 = 1.0_dp
           ELSE
              zbetai0 = (ztt-ixidx)*tbetai0(iqidx,ixidx+1)             &
                             +(ixidx+1._dp-ztt)*tbetai0(iqidx,ixidx)
              zbetai1 = (ztt-ixidx)*tbetai1(iqidx,ixidx+1)             &
                             +(ixidx+1._dp-ztt)*tbetai1(iqidx,ixidx)
           ENDIF
           zqp1b      = (zbetaqt(jl)-zbap1)*zbetai1 -                  &
                        (pbetass(jl,jk)-zbap1)*zbetai0 + pbetass(jl,jk)
           zgent      = MAX(pqm1(jl,jk)-zqp1b,-zxilb*zclcaux(jl))
           zgent      = MIN(zgent,zqsec*zqp1)              ! limit to qv
           zifrac     = zxib(jl)/MAX(zepsec,zxilb)
           zifrac     = MAX(MIN(zifrac,1.0_dp),0.0_dp)
           zgenti(jl) = zgent*zifrac
           zgentl(jl) = zgent*(1.0_dp-zifrac)
           IF (locc) THEN
              zxib(jl) = MAX(zxib(jl)+zgenti(jl)/zclcaux(jl),0.0_dp)
              zxlb(jl) = MAX(zxlb(jl)+zgentl(jl)/zclcaux(jl),0.0_dp)
           END IF
           zxilb       = zxib(jl)+zxlb(jl)
!
!       5.3 Deposition/sublimation of cloud ice and condensation/
!           evaporation of liquid water due to changes in water vapour
!           and temperature (advection, convection, turbulent mixing,
!           evaporation of rain, sublimation and melting of snow).
!           Translate PDF laterally to calculate cloud
!           after one timestep
!
           zqvdt     = zqvdt-zgent
           zdtdt     = zdtdt+zlvdcp(jl)*zgentl(jl)+zlsdcp(jl)*zgenti(jl)
           zqp1      = MAX(pqm1(jl,jk)+zqvdt,0.0_dp)
           ztp1      = ptm1(jl,jk)+zdtdt
           zdtdtstar = zdtdt+zclcaux(jl)*(zlc*pqte(jl,jk)*ztmst        &
              +zlvdcp(jl)*(zevp(jl)+zxlevap(jl)-zgentl(jl))            &
              +zlsdcp(jl)*(zsub(jl)+zxievap(jl)+zxisub(jl)-zgenti(jl)))
           zdqsat    = zdtdtstar*zdqsdt/(1._dp+zclcaux(jl)*zlcdqsdt)
           ztt       = (pbetass(jl,jk)-zqvdt+zdqsat-zbap1)/zbbap1
           ztt       = nbetax*MAX(MIN(ztt,1.0_dp),0.0_dp)
           ixidx     = INT(ztt)
           IF (ixidx == nbetax) THEN
              zbetai0 = 1.0_dp
              zbetai1 = 1.0_dp
           ELSE
              zbetai0 = (ztt-ixidx)*tbetai0(iqidx,ixidx+1)             &
                       +(ixidx+1._dp-ztt)*tbetai0(iqidx,ixidx)
              zbetai1 = (ztt-ixidx)*tbetai1(iqidx,ixidx+1)             &
                       +(ixidx+1._dp-ztt)*tbetai1(iqidx,ixidx)
           ENDIF
           zaa        = pbetaa(jl,jk)
           zqcdif     = (zbetaqt(jl)-zaa)*(1._dp-zbetai1)               &
                         +(zaa+zqvdt-pbetass(jl,jk)-zdqsat)*(1._dp-zbetai0)
           zqcdif     = MAX(0.0_dp,zqcdif)-zbetacl(jl)
           zqcdif     = MAX(zqcdif,-zxilb*zclcaux(jl))
           zqcdif     = MIN(zqcdif,zqsec*zqp1)             ! limit to qv
!
           IF (zqcdif .LT. 0.0_dp) THEN                 ! cloud dissipation
              zifrac   = zxib(jl)/MAX(zepsec,zxilb)
              zifrac   = MAX(MIN(zifrac,1.0_dp),0.0_dp)
              zdep(jl) = zqcdif*zifrac
              zcnd(jl) = zqcdif*(1.0_dp-zifrac)
           ELSE                                      ! cloud generation

              IF (lo2) THEN                          ! deposition
!
! changed for cirrus scheme, UL, 17.12.2006
!
                 IF (nicnc .LE. 1) THEN
                    zdep(jl) = zqcdif
                 ELSE IF (nicnc>1) THEN
                    zdep(jl) = zqinucl(jl,jk)
!                    zdep(jl) = MIN(zqinucl(jl,jk),zqcdif)
                 ENDIF
! end changed for cirrus scheme
                 zcnd(jl) = 0.0_dp
              ELSE                                   ! condensation
!--- Included/changed for prognostic CDNC/IC scheme --------------------
!    Use standard condensation for empirical Lin & Leaitch approach and 
!    explicit condensation after Levkov et al. 1992 for the explicit
!    activation schemes that allow for supersaturation:
                 IF (ncdnc==1) THEN
                    zcnd(jl) = zqcdif
                 ELSE IF (ncdnc>1) THEN
                    zdv      = 2.21_dp/papm1(jl,jk)
                    it       = NINT(ptm1(jl,jk)*1000._dp)
                    IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
                    it = MAX(MIN(it,jptlucu2),jptlucu1)
                    zesw     = tlucuaw(it)/papm1(jl,jk)
                    zesat    = zesw*papm1(jl,jk)*rv/rd
                    zast     = alv*(alv/(rv*ptm1(jl,jk))-1.0_dp)/               &
                               (0.024_dp*ptm1(jl,jk))
                    zbst     = rv*ptm1(jl,jk)/(zdv*zesat)
                    zgtp=1._dp/(zrho(jl,jk)*(zast+zbst))
                    zcond=0.5_dp*api*zdw0*swat(jl,jk,jrow)*zcdnc(jl,jk)  &
                          *zgtp*ztmst*zclcaux(jl)
                    zcnd(jl)= MIN(zcond,zqcdif)
                 END IF
!--- End included for CDNC/IC scheme -----------------------------------
                 zdep(jl) = 0.0_dp
              END IF
           END IF
          END IF !lcover
        END IF !ncctop
!
        IF ((.NOT.lcover) .OR. jk < ncctop) THEN
           zqcdif         = (zqvdt-zdqsat)*zclcaux(jl)
           zqcdif         = MAX(zqcdif,-zxilb*zclcaux(jl))
           zqcdif         = MIN(zqcdif,zqsec*zqp1)
           IF (zqcdif .LT. 0.0_dp) THEN                 ! cloud dissipation
              zifrac      = zxib(jl)/MAX(zepsec,zxilb)
              zifrac      = MAX(MIN(zifrac,1.0_dp),0.0_dp)
              zdep(jl)    = zqcdif*zifrac
              zcnd(jl)    = zqcdif*(1.0_dp-zifrac)
           ELSE                                      ! cloud generation
              IF (lo2) THEN                          ! deposition
!
! changed for cirrus scheme, UL, 17.12.2006
!
                 IF (nicnc .LE. 1) THEN
                    zdep(jl) = zqcdif
                 ELSE IF (nicnc>1) THEN
                    zdep(jl) = zqinucl(jl,jk)
                 ENDIF
! end changed for cirrus scheme
                 zcnd(jl) = 0.0_dp
              ELSE                                   ! condensation
!--- Included/changed for prognostic CDNC/IC scheme --------------------
!    Use standard condensation for empirical Lin & Leaitch approach and 
!    explicit condensation after Levkov et al. 1992 for the explicit
!    activation schemes that allow for supersaturation:
                 IF (ncdnc==1) THEN
                    zcnd(jl) = zqcdif
                 ELSE IF (ncdnc>1) THEN
                    zdv      = 2.21_dp/papm1(jl,jk)
                    it       = NINT(ptm1(jl,jk)*1000._dp)
                    IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
                    it = MAX(MIN(it,jptlucu2),jptlucu1)
                    zesw     = tlucuaw(it)/papm1(jl,jk)
                    zesat    = zesw*papm1(jl,jk)*rv/rd
                    zast     = alv*(alv/(rv*ptm1(jl,jk))-1.0_dp)/               &
                               (0.024_dp*ptm1(jl,jk))
                    zbst     = rv*ptm1(jl,jk)/(zdv*zesat)
                    zgtp=1._dp/(zrho(jl,jk)*(zast+zbst))
                    zcond=0.5_dp*api*zdw0*swat(jl,jk,jrow)*zcdnc(jl,jk)  &
                          *zgtp*ztmst*zclcaux(jl)
                    zcnd(jl)= MIN(zcond,zqcdif)
                 END IF
!--- End included for CDNC/IC scheme -----------------------------------
                 zdep(jl) = 0.0_dp
              END IF
           END IF
        END IF !lcover
!
!       5.4 Accounting for cloud evaporation in clear air and
!           checking for supersaturation
!
        ztp1tmp(jl) = ztp1+zlvdcp(jl)*zcnd(jl)+zlsdcp(jl)*zdep(jl)
        zqp1tmp(jl) = zqp1-zcnd(jl)-zdep(jl)
        zxip1       = MAX(zxised(jk)+zxite(jl)*ztmst-zxievap(jl)           &
                               +zgenti(jl)+zdep(jl),0.0_dp)
        lo2         = (ztp1tmp(jl) .LT. cthomi) .OR.                   &
                      (ztp1tmp(jl) .LT. tmelt .AND. zxip1 .GT. csecfrl &
                      .AND. zsusatw_2d(jl,jk) .LT. 1._dp)
        it          = NINT(ztp1tmp(jl)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zes         = MERGE(tlucua(it),tlucuaw(it),lo2)/papp1(jl,jk)
        zes         = MIN(zes,0.5_dp)
        LO          = zes<0.4_dp
        zcor        = 1._dp/(1._dp-vtmpc1*zes)
        zqsp1tmp    = zes*zcor
        zoversat    = zqsp1tmp*0.01_dp
        zrhtest     = MIN(pqm1(jl,jk)/zqsm1,1._dp)*zqsp1tmp
!
!changed for cirrus scheme, Ulrike Lohmann, 17.12.2006
!
        zesw        = tlucuaw(it)/papp1(jl,jk)
        zesw         = MIN(zesw,0.5_dp)
        zcorw        = 1._dp/(1._dp-vtmpc1*zesw)
        zqsp1tmpw    = zesw*zcorw
        zoversatw    = zqsp1tmpw*0.01_dp
!end changed for cirrus scheme
        it1         = it+1
        it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
        zqst1       = MERGE(tlucua(it1),tlucuaw(it1),lo2)/papp1(jl,jk)
        zqst1       = MIN(zqst1,0.5_dp)
        zqst1       = zqst1/(1._dp-vtmpc1*zqst1)
        zdqsdt      = (zqst1-zqsp1tmp)*1000._dp
        zlc         = MERGE(zlsdcp(jl),zlvdcp(jl),lo2)
        zlcdqsdt    = MERGE(zlc*zdqsdt,zqsp1tmp*zcor*tlucub(it),LO)
        zqcon       = 1._dp/(1._dp+zlcdqsdt)
!
        IF (lo2) THEN                                    ! ice cloud
           IF (nicnc .LE. 1) THEN
              IF (zqp1tmp(jl) .GT. zqsp1tmp+zoversat) THEN
                 zdepcor     = (zqp1tmp(jl)-zqsp1tmp-zoversat)*zqcon
                 zdep(jl)    = zdep(jl)+zdepcor
              END IF
!
!changed for cirrus scheme, Ulrike Lohmann, 17.12.2006
!
           ELSE IF (nicnc > 1) THEN
              IF (zqp1tmp(jl) .GT. zqsp1tmp+zoversat.AND.ztp1tmp(jl).GE.cthomi) THEN
                 zdepcor     = (zqp1tmp(jl)-zqsp1tmp-zoversat)*zqcon
                 zdep(jl)    = zdep(jl)+zdepcor
              END IF
              IF (zqp1tmp(jl) .GT. zqsp1tmpw+zoversatw.AND.ztp1tmp(jl).LT.cthomi) THEN
                 zdepcor     = (zqp1tmp(jl)-zqsp1tmpw-zoversatw)*zqcon
                 zdep(jl)    = zdep(jl)+zdepcor
              END IF
           ENDIF
           IF (zdep(jl) .GT. 0.0_dp .AND. zqp1tmp(jl) .LT. zrhtest        &
                .AND. zqsp1tmp .LE. zqsm1) THEN
              zdep(jl)    = zqp1-zrhtest
              zdep(jl)    = MAX(zdep(jl),0.0_dp)
           END IF
!end changed for cirrus scheme
        ELSE                                             ! water cloud
           IF (zqp1tmp(jl) .GT. zqsp1tmp+zoversat) THEN
              zcndcor     = (zqp1tmp(jl)-zqsp1tmp-zoversat)*zqcon
              zcnd(jl)    = zcnd(jl)+zcndcor
           END IF
           IF (zcnd(jl) .GT. 0.0_dp .AND. zqp1tmp(jl) .LT. zrhtest        &
                                 .AND. zqsp1tmp .LE. zqsm1) THEN
              zcnd(jl)    = zqp1-zrhtest
              zcnd(jl)    = MAX(zcnd(jl),0.0_dp)
           END IF
        END IF
!
!       5.5 Change of in-cloud water due to deposition/sublimation and
!           condensation/evaporation (input for cloud microphysics)
!
        zrelhum=zqp1tmp(jl)/zqsp1tmp
        zdepos =MAX(zdep(jl)+zgenti(jl),0.0_dp)
        
        zcond  =MAX(zcnd(jl)+zgentl(jl),0.0_dp)
        IF (locc) THEN
           zxib(jl) = MAX(zxib(jl)+zdep(jl)/zclcaux(jl),0.0_dp)
           zxlb(jl) = MAX(zxlb(jl)+zcnd(jl)/zclcaux(jl),0.0_dp)
        ELSEIF (zdepos>0.0_dp .OR. zcond>0.0_dp) THEN
           zclcaux(jl)=MAX(MIN(zrelhum,1.0_dp),0.01_dp)
           zxib(jl)   = zdepos/zclcaux(jl)
           zxlb(jl)   = zcond /zclcaux(jl)
        END IF
        IF (zclcaux(jl).GT.0._dp .AND. zxlb(jl).GT. cqtmin) THEN
!if there was no previous nucleation
           IF (zqlnuc(jl,jk) .LE. 0._dp) THEN    
              zqlnuc(jl,jk)=MAX(0._dp,pcdncact(jl,jk))
              zcdnc(jl,jk)=zcdnc(jl,jk)+zqlnuc(jl,jk)
              qnuc(jl,jk,jrow)=zqlnuc(jl,jk)/ztmst
           ENDIF
           zcdnc(jl,jk)=MAX(zcdnc(jl,jk),zcdnmin)
        ELSE
           zcdnc(jl,jk)=cqtmin
        ENDIF
        IF (zclcaux(jl).GT.0._dp .AND. zxib(jl).GT. cqtmin .AND. nicnc > 0) THEN
           ztc=MIN(0._dp,ptm1(jl,jk)-tmelt)
           zrieff=23.2_dp*EXP(0.015_dp*ztc)
           zrieff=MAX(zrieff,1._dp)
           zrih=-2261.236_dp+SQRT(5113188.044_dp+2809._dp*zrieff*zrieff*zrieff)
           zrid=1.e-6_dp*zrih**(1._dp/3._dp)
           zrid=MAX(1.e-6_dp,zrid)
           IF (zninucl(jl,jk) .LE. 0._dp) zicnc(jl,jk)=zicnc(jl,jk)+0.75_dp*zrho(jl,jk)*zxib(jl)/(api*zrid**3*zrhoice)
           zicnc(jl,jk)=MAX(zicnc(jl,jk),zicemin)
        ELSE
           zicnc(jl,jk)=cqtmin
        ENDIF
        ztp1tmp(jl) = ztp1+zlvdcp(jl)*zcnd(jl)+zlsdcp(jl)*zdep(jl)
!
!     ------------------------------------------------------------------
!       6.    Freezing of cloud water
!
!       6.1   Freezing of cloud water for T < 238 K
!
        IF (ztp1tmp(jl) .LE. cthomi) THEN
           zfrl(jl)  = zxlb(jl)*zclcaux(jl)
           zxib(jl)  = zxib(jl)+zxlb(jl)
           zxlb(jl)  = 0.0_dp
!--- Included for prognostic CDNC/IC scheme ----------------------------
           znfrl=zcdnc(jl,jk)
           zcdnc(jl,jk)=cqtmin
           qfre(jl,jk,jrow)=-znfrl/ztmst

           IF (nicnc .LE. 1) zicnc(jl,jk)=zicnc(jl,jk)+znfrl
!--- End included for CDNC/IC scheme -----------------------------------
     ENDIF

!       6.2   Freezing of cloud water between 238 and 273 K
!
        lo           = zxlb(jl) .GT. cqtmin .AND. ztp1tmp(jl) .LT. tmelt  &
                       .AND. ztp1tmp(jl) .GT. cthomi .AND. locc
        IF (lo) THEN
!---Changed for prognostic CDNC/IC scheme ------------------------------
!   (Replaced pacdnc by zcdnc which is set above)
!
           IF (zcdnc(jl,jk) .GE. zcdnmin) THEN
!--------- new freezing parameterisations (Lohmann & Diehl, 2006)
           ! corinna: included for contact/immersion freezing by dust and soot
           IF  (nicnc > 0) THEN  
              znaerinsol   = (pxtm1(jl,jk,idt_nai) +  pxtm1(jl,jk,idt_nci)  + &
                              pxtm1(jl,jk,idt_nki) + (pxtte(jl,jk,idt_nai)  + &
                              pxtte(jl,jk,idt_nci) +  pxtte(jl,jk,idt_nki)) * &
                              ztmst) / M_air * 1000._dp * zrho(jl,jk)

              zfracdusol   = zndusol(jl,jk)  /(zna(jl,jk)+zeps)
              zfracduinsolai = MIN(znduinsolai(jl,jk)/(znaerinsol+zeps),1._dp)
              zfracduinsolci = MIN(znduinsolci(jl,jk)/(znaerinsol+zeps),1._dp)
              zfracbcsol   = znbcsol(jl,jk)  /(zna(jl,jk)+zeps)
              zfracbcinsol = MIN(znbcinsol(jl,jk)/(znaerinsol+zeps),1._dp)
              zradl     = (0.75_dp*zxlb(jl)*zrho(jl,jk)                      &
                   /(api*rhoh2o*zcdnc(jl,jk)))**(1._dp/3._dp)
              zf1       = 4._dp*api*zradl*zcdnc(jl,jk)/zrho(jl,jk)          
              zetaair      = 1.e-5_dp*(1.718_dp+0.0049_dp*(ztp1tmp(jl)-tmelt)  &
                            -1.2e-5_dp*(ztp1tmp(jl)-tmelt)*(ztp1tmp(jl)-tmelt))
              zrdryki      = wetradius(jl,jk,iaiti)
              zrdryai      = wetradius(jl,jk,iacci)
              zrdryci      = wetradius(jl,jk,icoai)
              zccbcki      = 1._dp+1.26_dp*6.6E-8_dp/(zrdryki+zeps) *(101325._dp/papp1(jl,jk))   &
                                                *(ztp1tmp(jl)/tmelt)
              zccduai      = 1._dp+1.26_dp*6.6E-8_dp/(zrdryai+zeps) *(101325._dp/papp1(jl,jk))   &
                                                *(ztp1tmp(jl)/tmelt)
              zccduci      = 1._dp+1.26_dp*6.6E-8_dp/(zrdryci+zeps) *(101325._dp/papp1(jl,jk))   &
                                                *(ztp1tmp(jl)/tmelt)
              IF (zrdryki==0._dp) THEN 
                  zdfarbcki=0._dp
              ELSE
                  zdfarbcki    = ak*ztp1tmp(jl)*zccbcki/(6._dp*api*zetaair*(zrdryki+zeps))
              ENDIF
              IF (zrdryai==0._dp) THEN
                  zdfarduai=0._dp
              ELSE
                  zdfarduai    = ak*ztp1tmp(jl)*zccduai/(6._dp*api*zetaair*(zrdryai+zeps))
              ENDIF
              IF (zrdryci==0._dp) THEN
                  zdfarduci=0._dp
              ELSE
                  zdfarduci    = ak*ztp1tmp(jl)*zccduci/(6._dp*api*zetaair*(zrdryci+zeps))
              ENDIF
              zfrzcntdu = MIN(1._dp,MAX(0._dp,-(0.1014_dp*(ztp1tmp(jl)-tmelt)+0.3277_dp)))  ! montmorillonite
              !zfrzcntdu = MIN(1._dp,MAX(0._dp,-(0.1007_dp*(ztp1tmp(jl)-tmelt)+0.6935_dp)))  ! kaolinite 
              zfrzcntbc = MIN(1._dp,MAX(0._dp,-(0.0614_dp*(ztp1tmp(jl)-tmelt)+0.5730_dp)))
              zfrzcnt   = zxlb(jl)/zcdnc(jl,jk)*zrho(jl,jk)*zf1*                     &
                    (zfrzcntdu*(zdfarduai*zfracduinsolai+zdfarduci*zfracduinsolci)   &
                    +zfrzcntbc*zdfarbcki*zfracbcinsol)*(zcdnc(jl,jk)+zicnc(jl,jk))
              zfrzcnt   = zxlb(jl)*(1._dp-EXP(-zfrzcnt/zxlb(jl)*ztmst))
              znaimmdu  = 32.3_dp*zfracdusol     ! montmorillonite 
              !znaimmdu  = 6.15E-2_dp*zfracdusol     ! kaolinite 
              znaimmbc  = 2.91E-3_dp*zfracbcsol
              zfrzimm   = -(znaimmdu+znaimmbc)*EXP(tmelt-ztp1tmp(jl))*MIN(ptte(jl,jk),0._dp) 
              zfrzimm   = zxlb(jl)*(1._dp-EXP(-zfrzimm*zxlb(jl)/zcdnc(jl,jk)*ztmst))
              zfrl(jl)  = zfrzcnt + zfrzimm
              zfrl(jl)  = MAX(0.0_dp,MIN(zfrl(jl),zxlb(jl)))
              zfrln(jl) = MAX(0.0_dp,MIN(zcdnc(jl,jk)*zfrl(jl)/(zxlb(jl)+zeps), &
                zndusol(jl,jk)+znbcsol(jl,jk)+znduinsolai(jl,jk)            &
                +znduinsolci(jl,jk)+znbcinsol(jl,jk)-zicnc(jl,jk)))
              IF (zxlb(jl)-zfrl(jl) .GT. cqtmin) THEN
                 zfrln(jl)=MAX(0.0_dp,MIN(zfrln(jl),zcdnc(jl,jk)-zcdnmin))
              ELSE
                 zfrln(jl)=MAX(0.0_dp,MIN(zfrln(jl),zcdnc(jl,jk)-cqtmin))
              END IF
           ELSE
              zfrl(jl)  = 100._dp*(EXP(0.66_dp*(tmelt-ztp1tmp(jl)))-1._dp)         &
                   *zrho(jl,jk)/(rhoh2o*zcdnc(jl,jk))
              zfrl(jl)  = zxlb(jl)*(1._dp-EXP(-zfrl(jl)*ztmst*zxlb(jl)))
              zradl     = (0.75_dp*zxlb(jl)*zrho(jl,jk)                      &
                   /(api*rhoh2o*zcdnc(jl,jk)))**(1._dp/3._dp)
              zf1       = 4._dp*api*zradl*zcdnc(jl,jk)*2.e5_dp                  &
                   *(tmelt-3._dp-ztp1tmp(jl))/zrho(jl,jk)
              zf1       = MAX(0.0_dp,zf1)
              zfrl(jl)  = zfrl(jl) + zxlb(jl)*(1.-                        &
                   EXP(-ztmst*1.4e-20_dp*zf1/zxlb(jl)))
              zfrl(jl)  = MAX(0.0_dp,MIN(zfrl(jl),zxlb(jl)))
              !
              ! freezing of cloud droplets
              !
              IF (zxlb(jl)-zfrl(jl) .GT. cqtmin) THEN
                 zfrln(jl)=MAX(0.0_dp,MIN(zcdnc(jl,jk)*zfrl(jl)/(zxlb(jl)+zeps),zcdnc(jl,jk)-zcdnmin))
              ELSE
                 zfrln(jl)=MAX(0.0_dp,MIN(zcdnc(jl,jk)*zfrl(jl)/(zxlb(jl)+zeps),zcdnc(jl,jk)-cqtmin))
              END IF
           ENDIF
           zcdnc(jl,jk)=MAX(zcdnc(jl,jk)-zfrln(jl),cqtmin)
           IF (nicnc > 0) zicnc(jl,jk)=MAX(zicnc(jl,jk)+zfrln(jl),cqtmin)
!--- End included for CDNC/IC scheme -----------------------------------
           zxlb(jl)  = zxlb(jl)-zfrl(jl)
           zxib(jl)  = zxib(jl)+zfrl(jl)
           zfrl(jl)  = zfrl(jl)*zclcaux(jl)
!-------corinna: Bergeron-Findeisen-Process
           IF (locc .AND. zdep(jl)>0._dp .AND. zxlb(jl)>0._dp .AND. zxib(jl)>csecfrl &
                .AND. zsusatw_2d(jl,jk) .LT. 1._dp) THEN
              zzevp=zxlb(jl)/ztmst*zclcaux(jl)
              pxlte(jl,jk)=pxlte(jl,jk)-zzevp
              pxite(jl,jk)=pxite(jl,jk)+zzevp
              ptte(jl,jk) =ptte(jl,jk)+(zlsdcp(jl)-zlvdcp(jl))*zzevp
              zcdnc(jl,jk)=cqtmin
              zxlb(jl)=0._dp
           END IF
!-------End corinna: Bergeron-Findeisen-Process
        END IF
     ENDIF
!
610  END DO
!
!     IF (lookupoverflow) CALL lookuperror ('cloud (2)    ')
     IF (lookupoverflow) THEN
       status_string = 'lookuperror: icnc/cdnc - cloud (2)'
       RETURN
     ENDIF
!
!     ------------------------------------------------------------------
!       7.  Cloud physics and precipitation fluxes at the surface
!
!ham_ps: cdir circumvents bug in sxf90 compiler
!CDIR NOMOVEDIV
     DO 701 jl = 1,kproma
        locc     = zclcaux(jl) .GT. 0.0_dp
        zclcstar = MIN(zclcaux(jl),zclcpre(jl))
        zauloc   = cauloc*zdz(jl)/5000._dp
        zauloc   = MAX(MIN(zauloc,clmax),clmin)
!
        jb=knvb(jl)
        lo=(jb.GE.jbmin .AND. jb.LE.jbmax .AND. pvervel(jl,jk).GT.0._dp)
        lo1=(jk.EQ.jb .OR. jk.EQ.jb+1)
        IF(lo .AND. lo1 .AND. lonacc) zauloc= 0.0_dp
!
        zqrho    = 1.3_dp/zrho(jl,jk)
        zxlb(jl) = MAX(zxlb(jl),1.e-20_dp)
        zxib(jl) = MAX(zxib(jl),1.e-20_dp)
! liquid water and snow content are stored 
! before the reduction by outfalling rain
! (necessary for nucleation scavenging)
        plwc(jl,jk) = zxlb(jl)
        piwc(jl,jk) = zxib(jl)      
        IF (zclcpre(jl) .GT. 0.0_dp .AND. zrfl(jl) .GT. cqtmin) THEN
           zxrp1 = (zrfl(jl)/(zclcpre(jl)*12.45_dp*SQRT(zqrho)))**(8._dp/9._dp)
        ELSE
           zxrp1 = 0.0_dp
        ENDIF

        IF (zclcpre(jl) .GT. 0.0_dp .AND. zsfl(jl) .GT. cqtmin) THEN
           zxsp1 = (zsfl(jl)/(zclcpre(jl)*cvtfall))**(1._dp/1.16_dp)
        ELSE
           zxsp1 = 0.0_dp
        END IF
!
!       7.1   Warm clouds: Coalescence processes after Beheng (1994):
!             Autoconversion of cloud droplets and collection of cloud
!             droplets by falling rain. Accretion of cloud droplets by
!             falling snow (zsacl) is calculated under 7.2
!
        IF (locc) THEN
           IF (zxlb(jl) > cqtmin .AND. zcdnc(jl,jk) .GE. zcdnmin) THEN

!--- Included alternative autoconversion parameterisation --------------

           IF (nauto==2) THEN

!           Autoconversion rate from Khairoutdinov and Kogan, 2000

              zraut    = ccraut*1350._dp*(zcdnc(jl,jk)*1.e-6_dp)**(-1.79_dp)
              zexm1    = 2.47_dp-1.0_dp
              zexp     = -1._dp/zexm1
              zraut    = zxlb(jl)*(1._dp-(1._dp+zraut*ztmst*zexm1*zxlb(jl)      &
                                                       **zexm1)**zexp)
              zraut    = MIN(zxlb(jl),zraut)
              zxlb(jl) = zxlb(jl)-zraut
              zrac1    = 6._dp*zxrp1*ztmst
              zrac1    = zxlb(jl)*(1._dp-EXP(-zrac1))
              zxlb(jl) = zxlb(jl)-zrac1
              zrac2    = 6._dp*zauloc*zrho(jl,jk)*zraut*ztmst
              zrac2    = zxlb(jl)*(1._dp-EXP(-zrac2))
              zxlb(jl) = zxlb(jl)-zrac2
              zrpr(jl) = zrpr(jl)+zclcaux(jl)*(zraut+zrac2)+zclcstar*zrac1

!---Included for in-cloud scavenging (Philip Stier, 26/11/03):----------
              zmratepr(jl,jk)=zraut+zrac1+zrac2
!---End Included for scavenging-----------------------------------------
!--- Included for prognostic CDNC/IC scheme ----------------------------
              zxlbold=zxlb(jl)+zrac1+zrac2+zraut
              zrprn(jl)=(zraut+zrac1+zrac2)/(zxlbold+zeps)
              IF (zxlb(jl) .GT. cqtmin) THEN
                 zrprn(jl)=MIN(zcdnc(jl,jk)*zrprn(jl),zcdnc(jl,jk)-zcdnmin)
              ELSE
                 zrprn(jl)=MIN(zcdnc(jl,jk)*zrprn(jl),zcdnc(jl,jk))
              END IF
              zcdnc(jl,jk)=MAX(zcdnc(jl,jk)-zrprn(jl),cqtmin)
!--- End included for CDNC/IC scheme ------------------------------------
!--- End included alternative autoconversion parameterisation ----------
!--- Changed for alternative autoconversion parameterisation -----------

           ELSE IF (nauto==1) THEN

!---Changed for prognostic CDNC/IC scheme ------------------------------
!   (Replaced pacdnc by zcdnc which is set above)
              zraut    = ccraut*1.2e27_dp/zrho(jl,jk)*(zcdnc(jl,jk)*1.e-6_dp)  &
                         **(-3.3_dp)*(zrho(jl,jk)*1.e-3_dp)**4.7_dp
!--- End changed for prognostic CDNC/IC scheme -------------------------

              zexm1    = 4.7_dp-1.0_dp
              zexp     = -1._dp/zexm1
              zraut    = zxlb(jl)*(1._dp-(1._dp+zraut*ztmst*zexm1*zxlb(jl)      &
                                                       **zexm1)**zexp)
!--- Included for prognostic CDNC/IC scheme ----------------------------
              zrautn=zraut*7.7e9_dp*zrho(jl,jk)
              zself=1.289e10_dp*(zrho(jl,jk)*zxlb(jl))**2*ztmst
              zrautself=MIN(zrautn+zself,zcdnc(jl,jk))
              zcdnc(jl,jk)=MAX(zcdnc(jl,jk)-zrautself,cqtmin)
!--- End included for CDNC/IC scheme -----------------------------------
              zxlb(jl) = zxlb(jl)-zraut
              zrac1    = 6._dp*zxrp1*ztmst
              zrac1    = zxlb(jl)*(1._dp-EXP(-zrac1))
              zxlb(jl) = zxlb(jl)-zrac1
              zrac2    = 6._dp*zauloc*zrho(jl,jk)*zraut*ztmst
              zrac2    = zxlb(jl)*(1._dp-EXP(-zrac2))
              zxlb(jl) = zxlb(jl)-zrac2
              zrpr(jl) = zrpr(jl)+zclcaux(jl)*(zraut+zrac2)+zclcstar*zrac1
!---Included for in-cloud scavenging (Philip Stier, 26/11/03):----------
              zmratepr(jl,jk)=zraut+zrac1+zrac2
!---End Included for scavenging-----------------------------------------
!--- Included for prognostic CDNC/IC scheme ----------------------------
              zxlbold=zxlb(jl)+zrac1+zrac2
              zraccn=(zrac1+zrac2)/(zxlbold+zeps)
              zraccn=MIN(zcdnc(jl,jk)*zraccn,zcdnc(jl,jk))
              zcdnc(jl,jk)=MAX(zcdnc(jl,jk)-zraccn,cqtmin)
              zrprn(jl)=zrautself+zraccn
!--- End included for CDNC/IC scheme -----------------------------------
           END IF
           !--- End changed for alternative autoconversion parameterisation -------
        ENDIF
        IF (zxib(jl) > cqtmin) THEN

!       7.2  Cold clouds:
!            Conversion of cloud ice to snow after Levkov et al. 1992:
!            Aggregation of ice crystals to snow and accretion of ice
!            by falling snow.
!            Accrection of cloud droplets by falling snow.
!            Effective radius of ice crystals after Moss (1995)
!
           zrieff    = 83.8_dp*(zxib(jl)*zrho(jl,jk)*1000._dp)**0.216_dp
           zrieff    = MIN(MAX(zrieff,ceffmin),ceffmax)
           zrih      = -2261._dp+SQRT(5113188._dp+2809._dp*zrieff*zrieff*zrieff)
           zris       = 1.e-6_dp*zrih**(1._dp/3._dp)

           zcolleffi = EXP(0.025_dp*(ztp1tmp(jl)-tmelt))
           zc1       = 17.5_dp*zrho(jl,jk)/crhoi*zqrho**0.33_dp
           zdt2      = -6._dp/zc1*LOG10(zris*1.e4_dp)
           zsaut     = ccsaut/zdt2
           zsaut     = zxib(jl)*(1._dp-1._dp/(1._dp+zsaut*ztmst*zxib(jl)))
           zxib(jl)  = zxib(jl)-zsaut
           zsaci1    = 0.0_dp
           zsaci2    = 0.0_dp
           zsacl1    = 0.0_dp
           zsacl2    = 0.0_dp
!--- Included for prognostic CDNC/IC scheme ----------------------------
           zsacl1in  = 0.0_dp      ! in-cloud values saved
           zsacl2in  = 0.0_dp      ! 
!--- End included for CDNC/IC scheme -----------------------------------
           IF (zxsp1 .GT. cqtmin .AND. zxlb(jl).GT.cqtmin .AND. zcdnc(jl,jk).GE. zcdnmin) THEN
!--- Included for prognostic CDNC/IC scheme ----------------------------
              IF (losacl) THEN
                 !
                 ! included size depedent accretion rate (Lohmann, JAS, 2004)
                 !
                 IF (nicnc > 0) THEN
                    zscnc=zicnc(jl,jk)
                    IF (jk > 1) zscnc=MAX(zicemin,MIN(zsprn(jl,jk-1),zicnc(jl,jk)))
                 ELSE
                    zscnc=0.75_dp*zrho(jl,jk)*zsaut/(api*zrsnow**3*crhosno)
                    zscnc=MAX(zsnowmin,zscnc)
                 ENDIF
                 zdw=MAX((6._dp*zrho(jl,jk)*zxlb(jl)/(api*rhoh2o*zcdnc(jl,jk))) &
                      **(1._dp/3._dp),1.e-6_dp)
                 zdplanar=MAX(20.e-6_dp,SQRT(zxsp1*1.e3_dp/(zscnc*3.8e-4_dp))*1.e-2_dp)
                 zusnow = 2.34_dp*(zdplanar*100._dp)**0.3_dp*(1.3_dp/zrho(jl,jk))**0.35_dp
                 zudrop = 1.19e4_dp*(0.5_dp*zdw*100._dp)**2*(1.3_dp/zrho(jl,jk))**0.35_dp
                 zstokes=MAX(2._dp*(zusnow-zudrop)*zudrop/(zdplanar*g),cqtmin)
                 zviscos=(1.512_dp + 0.0052_dp*(ptm1(jl,jk)-233.15_dp))*1.e-5_dp
                 zrey=MAX(zrho(jl,jk)*zdplanar*zusnow/zviscos,cqtmin)
                 !
                 zstcrit=1._dp
                 IF (zrey.LE.5._dp) zstcrit=5.52_dp*zrey**(-1.12_dp)
                 IF (zrey.GT.5._dp.AND.zrey.LT.40._dp) zstcrit=1.53_dp*zrey**(-0.325_dp)

                 zhelp=MAX(MIN(0.2_dp*(LOG10(zstokes)-LOG10(zstcrit) &
                      -2.236_dp)**2,1._dp-cqtmin),0._dp)
                 zcsacl1=SQRT(1._dp-zhelp)

                 IF (zrey.GE.40._dp) THEN
                    IF (zstokes .LE. 0.06_dp) THEN
                       zcsacl1=1.034_dp*zstokes**1.085_dp
                    ELSE IF (zstokes.GT.0.06_dp .AND. zstokes .LE. 0.25_dp) THEN
                       zcsacl1=0.787_dp*zstokes**0.988_dp
                    ELSE IF (zstokes.GT.0.25_dp .AND. zstokes .LE. 1._dp) THEN
                       zcsacl1=0.7475_dp*LOG10(zstokes)+0.65_dp
                    ELSE
                       zcsacl1=(zstokes+1.1_dp)**2/(zstokes+1.6_dp)**2
                    ENDIF
                 ENDIF
                 zcsacl=MAX(0.01_dp,MIN(zcsacl1,1._dp))
              ELSE
                 zcsacl=ccsacl
              ENDIF
!--- End included for CDNC/IC scheme -----------------------------------
              zlamsm    = (zxsp1/(api*crhosno*cn0s))**0.8125_dp
              zsaci1    = api*cn0s*3.078_dp*zlamsm*zqrho**0.5_dp
              zsacl1    = zxlb(jl)*(1._dp-EXP(-zsaci1*zcsacl*ztmst))
!--- Included for prognostic CDNC/IC scheme ----------------------------
              zsacl1in  = zsacl1
!--- End included for CDNC/IC scheme -----------------------------------
              zxlb(jl)  = zxlb(jl)-zsacl1
              zsacl1    = zclcstar*zsacl1
              zsaci1    = zsaci1*zcolleffi*ztmst
              zsaci1    = zxib(jl)*(1._dp-EXP(-zsaci1))
              zxib(jl)  = zxib(jl)-zsaci1
           END IF
           zxsp2        = zauloc*zrho(jl,jk)*zsaut
           IF (zxsp2 .GT. cqtmin .AND. zxlb(jl).GT.cqtmin .AND. zcdnc(jl,jk).GE. zcdnmin) THEN
              IF (losacl) THEN
                 !
                 !uls - included size depedent accretion rate (Lohmann, JAS, 2004)
                 !
                 IF (nicnc > 0) THEN
                    zscnc=zicnc(jl,jk)
                    IF (jk > 1) zscnc=MAX(zicemin,MIN(zsprn(jl,jk-1),zicnc(jl,jk)))
                 ELSE
                    zscnc=0.75_dp*zrho(jl,jk)*zsaci1/(api*zrsnow**3*crhosno)
                    zscnc=MAX(zsnowmin,zscnc)
                 ENDIF
                 zdw=MAX((6._dp*zrho(jl,jk)*zxlb(jl)/(api*rhoh2o*zcdnc(jl,jk))) &
                      **(1._dp/3._dp),1.e-6_dp)
                 zdplanar=MAX(20.e-6_dp,SQRT(zxsp2*1.e3_dp/(zscnc*3.8e-4_dp))*1.e-2_dp)
                 zusnow = 2.34_dp*(zdplanar*100._dp)**0.3_dp*(1.3_dp/zrho(jl,jk))**0.35_dp
                 zudrop = 1.19e4_dp*(0.5_dp*zdw*100._dp)**2*(1.3_dp/zrho(jl,jk))**0.35_dp
                 zstokes=MAX(2._dp*(zusnow-zudrop)*zudrop/(zdplanar*g),cqtmin)
                 zviscos=(1.512_dp + 0.0052_dp*(ptm1(jl,jk)-233.15_dp))*1.e-5_dp
                 zrey=MAX(zrho(jl,jk)*zdplanar*zusnow/zviscos,cqtmin)
                 !
                 zstcrit=1._dp
                 IF (zrey.LE.5._dp) zstcrit=5.52_dp*zrey**(-1.12_dp)
                 IF (zrey.GT.5._dp .AND. zrey .LT. 40._dp) zstcrit=1.53_dp*zrey**(-0.325_dp)

                 zhelp=MAX(MIN(0.2_dp*(LOG10(zstokes)-LOG10(zstcrit) &
                      -2.236_dp)**2,1._dp-cqtmin),0._dp)
                 zcsacl1=SQRT(1._dp-zhelp)

                 IF (zrey.GE.40._dp) THEN
                    IF (zstokes .LE. 0.06_dp) THEN
                       zcsacl1=1.034_dp*zstokes**1.085_dp
                    ELSE IF (zstokes.GT.0.06_dp .AND. zstokes .LE. 0.25_dp) THEN
                       zcsacl1=0.787_dp*zstokes**0.988_dp
                    ELSE IF (zstokes.GT.0.25_dp .AND. zstokes .LE. 1._dp) THEN
                       zcsacl1=0.7475_dp*LOG10(zstokes)+0.65_dp
                    ELSE
                       zcsacl1=(zstokes+1.1_dp)**2/(zstokes+1.6_dp)**2
                    ENDIF
                 ENDIF
                 zcsacl=MAX(0.01_dp,MIN(zcsacl1,1._dp))
              ELSE
                 zcsacl=ccsacl
              ENDIF
              zlamsm    = (zxsp2/(api*crhosno*cn0s))**0.8125_dp
              zsaci2    = api*cn0s*3.078_dp*zlamsm*zqrho**0.5_dp
              zsacl2    = zxlb(jl)*(1._dp-EXP(-zsaci2*zcsacl*ztmst))
!--- Included for prognostic CDNC/IC scheme ----------------------------
              zsacl2in  = zsacl2
!--- End included for CDNC/IC scheme -----------------------------------
              zxlb(jl)  = zxlb(jl)-zsacl2
              zsacl2    = zclcaux(jl)*zsacl2
              zsaci2    = zsaci2*zcolleffi*ztmst
              zsaci2    = zxib(jl)*(1._dp-EXP(-zsaci2))
              zxib(jl)  = zxib(jl)-zsaci2
           END IF
           zsacl(jl)    = zsacl1+zsacl2
           zspr(jl)     = zspr(jl)+zclcaux(jl)*(zsaut+zsaci2)          &
                                  +zclcstar*zsaci1
!--- Included for prognostic CDNC/IC scheme ----------------------------
           zxlbold=zxlb(jl)+zsacl1in+zsacl2in
           IF (zxlb(jl) .GT. cqtmin) THEN
              zsacln(jl)=MAX( MIN(zcdnc(jl,jk)*(zsacl1in+zsacl2in)/(zxlbold+zeps),&
                   zcdnc(jl,jk)-zcdnmin ),0._dp )
           ELSE
              zsacln(jl)=MIN(zcdnc(jl,jk)*(zsacl1in+zsacl2in)/(zxlbold+zeps),&
                   zcdnc(jl,jk))
           END IF
           zcdnc(jl,jk)=MAX(zcdnc(jl,jk)-zsacln(jl),cqtmin)
!
!          secondary ice crystal production after Levkov et al. 1992
!          sink for snow, source for ice crystals
!
           zsecprod=0._dp
           IF ((zxsp1+zxsp2).GT.zepsec.AND.zxlb(jl).GT.zepsec.AND.   &
                ztp1tmp(jl).GT.265.2_dp.AND.ztp1tmp(jl).LT.270.2_dp) THEN
              IF (losacl .AND. nicnc > 0) THEN
                 !
                 !uls - included size depedent accretion rate (Lohmann, JAS, 2004)
                 !
                 zscnc=zicnc(jl,jk)
                 IF (jk > 1) zscnc=MAX(zicemin,MIN(zsprn(jl,jk-1),zicnc(jl,jk)))
                 zdw=MAX((6._dp*zrho(jl,jk)*zxlb(jl)/(api*rhoh2o*zcdnc(jl,jk))) &
                      **(1._dp/3._dp),1.e-6_dp)
                 zdplanar=MAX(20.e-6_dp,SQRT((zxsp1+zxsp2)*1.e3_dp/(zscnc*3.8e-4_dp))*1.e-2_dp)
                 zusnow = 2.34_dp*(zdplanar*100._dp)**0.3_dp*(1.3_dp/zrho(jl,jk))**0.35_dp
                 zudrop = 1.19e4_dp*(0.5_dp*zdw*100._dp)**2*(1.3_dp/zrho(jl,jk))**0.35_dp
                 zstokes=MAX(2._dp*(zusnow-zudrop)*zudrop/(zdplanar*g),cqtmin)
                 zviscos=(1.512_dp + 0.0052_dp*(ptm1(jl,jk)-233.15_dp))*1.e-5_dp
                 zrey=MAX(zrho(jl,jk)*zdplanar*zusnow/zviscos,cqtmin)
                 !
                 zstcrit=1._dp
                 IF (zrey.LE.5._dp) zstcrit=5.52_dp*zrey**(-1.12_dp)
                 IF (zrey.GT.5._dp.AND.zrey.LT.40._dp) zstcrit=1.53_dp*zrey**(-0.325_dp)

                 zhelp=MAX(MIN(0.2_dp*(LOG10(zstokes)-LOG10(zstcrit) &
                      -2.236_dp)**2,1._dp-cqtmin),0._dp)
                 zcsacl1=SQRT(1._dp-zhelp)

                 IF (zrey.GE.40._dp) THEN
                    IF (zstokes .LE. 0.06_dp) THEN
                       zcsacl1=1.034_dp*zstokes**1.085_dp
                    ELSE IF (zstokes.GT.0.06_dp .AND. zstokes .LE. 0.25_dp) THEN
                       zcsacl1=0.787_dp*zstokes**0.988_dp
                    ELSE IF (zstokes.GT.0.25_dp .AND. zstokes .LE. 1._dp) THEN
                       zcsacl1=0.7475_dp*LOG10(zstokes)+0.65_dp
                    ELSE
                       zcsacl1=(zstokes+1.1_dp)**2/(zstokes+1.6_dp)**2
                    ENDIF
                 ENDIF
                 zcsacl=MAX(0.01_dp,MIN(zcsacl1,1._dp))
              ELSE
                 zcsacl=ccsacl
              ENDIF
              zlams2=((zxsp1+zxsp2)/(api*crhosno*cn0s))**0.875_dp
              zj=zcsacl*api*zrho(jl,jk)*zxlb(jl)*cn0s*0.831_dp*zlams2*         &
                 (g*crhosno/(0.75_dp*zrho(jl,jk)*zcdi))**0.5_dp/zmw0
              zpn=MAX(0.00285_dp*zj,0._dp)
              zsecprod=ztmst*zpn*zmi0/zrho(jl,jk)
              zsecprod=MAX(0._dp,MIN((zxsp1+zxsp2)/zrho(jl,jk),zsecprod))
              zxib(jl)=zxib(jl)+zsecprod
           END IF
           !
           zspr(jl)=MAX(zspr(jl)-zclcstar*zsecprod,0._dp)           

           IF (zxib(jl).GT.zepsec .AND. zicnc(jl,jk).GE.zicemin .AND. nicnc > 0) THEN
              zxibold=MAX(zxib(jl)+zsaut+zsaci1+zsaci2-zsecprod,0._dp)
              zsprn1=zicnc(jl,jk)*(zsaci1+zsaci2+zsaut)/(zxibold+zeps)
              zself=zc1*0.5_dp*zicnc(jl,jk)*ztmst*zxib(jl)
              zsecprodn=zrho(jl,jk)/zmi0*zsecprod
              zsprnself=MIN(zsprn1+zself-zsecprodn,zicnc(jl,jk))
              zicnc(jl,jk)=MAX(zicnc(jl,jk)-zsprnself,cqtmin)
              zsprn(jl,jk)=zsprnself
           END IF
!--- End included for CDNC/IC scheme -----------------------------------
!---Included for in-cloud scavenging (Philip Stier, 25/11/03):----------
           prate_s(jl,jk) = prate_s(jl,jk) + zsaut+zsaci1+zsaci2-zsecprod 
!---End Included for scavenging-----------------------------------------

        END IF
        END IF
!
!       7.3 Updating precipitation fluxes. In the lowest layer (klev),
!           the sedimentation sink of cloud ice is balanced
!           by precipitation at the surface (through 'zzdrs').
!           Fraction of precipitating clouds (zclcpre) used for the
!           calculation of evaporation/sublimation of rain/snow in
!           the next layer
!
        zzdrr          = zcons2*zdp(jl)*zrpr(jl)
        zzdrs          = zcons2*zdp(jl)*(zspr(jl)+zsacl(jl))
        IF (jk .EQ. klev) THEN
           zzdrs       = zzdrs+zxiflux(jl)
           zcons       = zcons2*zdp(jl)/(zlsdcp(jl)-zlvdcp(jl))
           zsnmlt      = MIN(zxsec*zzdrs,zcons                         &
                                         *MAX(0._dp,(ztp1tmp(jl)-tmelt)))
           zzdrr       = zzdrr+zsnmlt
           zzdrs       = zzdrs-zsnmlt
           zsmlt(jl)   = zsmlt(jl)+zsnmlt/(zcons2*zdp(jl))
        END IF
        zpretot        = zrfl(jl)+zsfl(jl)
        zpredel        = zzdrr+zzdrs
        lo=(zpretot .GT. zpredel)
        zclcpre(jl)    = MERGE(zclcpre(jl),zclcaux(jl),lo)
        zpresum        = zpretot+zpredel
        IF (zpresum .GT. cqtmin) THEN
           zclcpre(jl) = MAX(zclcpre(jl),(zclcaux(jl)*zpredel          &
                                         +zclcpre(jl)*zpretot)/zpresum)
           zclcpre(jl) = MIN(zclcpre(jl),1.0_dp)
           zclcpre(jl) = MAX(zclcpre(jl),0.0_dp)
        ELSE
           zclcpre(jl) = 0.0_dp
        END IF
! rain and snow flux considering incoming rain, melting of snow, 
! droplet evaporation / sublimation , but no new production of rain or snow 
! in that layer....
! (neccessary for impaction scavenging)
        pfrain_no(jl,jk)   = zrfl(jl) - zcons2*zdp(jl)*zevp(jl)  
        pfsnow_no(jl,jk)   = zsfl(jl) - zcons2*zdp(jl)*zsub(jl)
! precipitating cloud cover of this layer is used for the next lower layer 
! to estimate the part of the cloud cover in which rain impacts
        pr_cover(jl,jk) = zclcpre(jl)
        ! mz_ht_20071113-

        zsfl(jl)       = zsfl(jl)+zzdrs-zcons2*zdp(jl)*zsub(jl)
! rain and snow flux out of the bottom of this layer
        pfrain(jl,jk) = zrfl(jl)
        pfsnow(jl,jk) = zsfl(jl)
!
701  END DO
!
!     ------------------------------------------------------------------
!       8.    Updating tendencies of t, q, xl, xi and final cloud cover
!
     DO 811 jl = 1,kproma
!
!       8.10   Cloud cover scheme tendencies
!
        IF (jk.GE.ncctop) THEN
           locc            = zclcaux(jl) .GT. 0.0_dp
!
!          Source terms from convection
!          Skewness:
!
!-------------------Added by Junhua Zhang for CONV Micro-----------------------
          IF(ncvmicro>0) THEN
           zconvskew(jl)   = cbeta_cs * (pxtecl(jl,jk)+pxteci(jl,jk)+pqtec(jl,jk)) &
                                  /pbetass(jl,jk)
                  ELSE
           zconvskew(jl)   = cbeta_cs * (pxtec(jl,jk)+pqtec(jl,jk))    &
                                  /pbetass(jl,jk)
          ENDIF
!-------------------------------------end-------------------------------------
           zconvskew(jl)   = MIN(zconvskew(jl),                        &
                                (cbeta_pq_max-pxskew(jl,jk))/zdtime)
!
!          Convective width now diagnosed, assuming 'a' unchanged:
!
           IF (pqm1(jl,jk) >= pbetass(jl,jk)) THEN
              zskewp1      = pxskew(jl,jk)+zconvskew(jl)*zdtime
              zbbap1       = zwide(jl)*(cbeta_pq+zskewp1)/             &
                                       (cbeta_pq+pxskew(jl,jk))
              zconvvar(jl) = (zbbap1-zwide(jl))/zdtime
           ELSE
              zconvvar(jl) = 0.0_dp
           ENDIF
!
!       8.11 Simple linearized effect of microphysics on skewness
!
           IF (pbetaa(jl,jk) < pbetass(jl,jk) .AND.                       &
               pbetab(jl,jk) > pbetass(jl,jk)) THEN
              zmdelb = (zxlte(jl)+zxite(jl))*ztmst                     &
                       -zrpr(jl)-zsacl(jl)-zspr(jl)+zcnd(jl)+zdep(jl)  &
                       +zgenti(jl)+zgentl(jl)
              zmdelb = MAX(0.0_dp,MIN(1.0_dp,-zmdelb/MAX(zepsec,zbetacl(jl))))
              zmdelb = (pbetass(jl,jk)-pbetab(jl,jk))*zmdelb
              zmqp1  = (pbetab(jl,jk)+zmdelb-pbetaa(jl,jk))            &
                        *cbeta_pq/(zbetaqt(jl)-pbetaa(jl,jk))          &
                                                          - cbeta_pq
              zmqp1  = MAX(MIN(zmqp1,cbeta_pq_max),cbeta_pq)
              zmicroskew(jl) = MIN(0.0_dp,(zmqp1-pxskew(jl,jk))/zdtime)
           ENDIF
!
!       8.2   New skewness and variance
!
           zxskewte(jl)    = zconvskew(jl)                             &
                             +zmicroskew(jl)+zturbskew(jl)
           zxvarte(jl)     = zconvvar(jl)+zturbvar(jl)
!
           zvarp1          = pxvar(jl,jk)+zxvarte(jl)*zdtime
           zskewp1         = pxskew(jl,jk)+zxskewte(jl)*zdtime
!
           pxskew(jl,jk)   = MAX(MIN(zskewp1,cbeta_pq_max),cbeta_pq)
           zvarmx          = zbetaqt(jl)*(1._dp+pxskew(jl,jk)/cbeta_pq)
           pxvar(jl,jk)    = MAX(MIN(zvarp1,zvarmx),zvartg(jl))
!
        ENDIF !jk >= ncctop
!
!       8.3   Tendencies of thermodynamic variables
!             Attn: The terms zxisub and zximlt do not appear in
!                   pxite because these processes have already been
!                   included in pxite via changes in cloud ice
!                   sedimentation (see 3.1, 3.2 and 4)
!
        pqte(jl,jk)  = pqte(jl,jk)                                     &
                        +(-zcnd(jl)-zgentl(jl)+zevp(jl)+zxlevap(jl)    &
                          -zdep(jl)-zgenti(jl)+zsub(jl)+zxievap(jl)    &
                                          +zxisub(jl))/ztmst
        ptte(jl,jk)  = ptte(jl,jk)+(zlvdcp(jl)                         &
                        *(zcnd(jl)+zgentl(jl)-zevp(jl)-zxlevap(jl))    &
                                  +zlsdcp(jl)                          &
                        *(zdep(jl)+zgenti(jl)-zsub(jl)-zxievap(jl)     &
                        -zxisub(jl))+(zlsdcp(jl)-zlvdcp(jl))           &
                        *(-zsmlt(jl)-zimlt(jl)-zximlt(jl)+zfrl(jl)     &
                                           +zsacl(jl)))/ztmst
        pxlte(jl,jk) = pxlte(jl,jk)+zxlte(jl)                          &
                        +(zimlt(jl)+zximlt(jl)-zfrl(jl)-zrpr(jl)       &
                        -zsacl(jl)+zcnd(jl)+zgentl(jl)-zxlevap(jl))    &
                                                       /ztmst
        pxite(jl,jk) = pxite(jl,jk)+zxite(jl)+(zfrl(jl)-zspr(jl)       &
                          +zdep(jl)+zgenti(jl)-zxievap(jl))/ztmst
        ztp1         = ptm1(jl,jk)+ptte(jl,jk)*ztmst
        zqp1         = pqm1(jl,jk)+pqte(jl,jk)*ztmst
        zxlp1        = pxlm1(jl,jk)+pxlte(jl,jk)*ztmst
        zxip1        = pxim1(jl,jk)+pxite(jl,jk)*ztmst
!
!--- Included for prognostic CDNC/IC scheme ----------------------------

        !--- Calculate new total tendency of CDNC:

        pxtte(jl,jk,idt_cdnc) = pxtte(jl,jk,idt_cdnc) + ( zcdnc(jl,jk) /  &
                                ( zrho(jl,jk) * 1000._dp / M_air) -       &
                                ( pxtm1(jl,jk,idt_cdnc) +                 &
                                  pxtte(jl,jk,idt_cdnc) * ztmst ) ) / ztmst
        !--- Update CDNC for radiation:

        pacdnc(jl,jk)=zcdnc(jl,jk)

        !--- Calculate new total tendency of ICNC:

        IF (nicnc > 0) &
          pxtte(jl,jk,idt_icnc) = pxtte(jl,jk,idt_icnc) + ( zicnc(jl,jk) /  &
                                ( zrho(jl,jk) * 1000._dp / M_air) -       &
                                ( pxtm1(jl,jk,idt_icnc) +                 &
                                  pxtte(jl,jk,idt_icnc) * ztmst ) ) / ztmst
        !--- Diagnostics:

        qaut(jl,jk,jrow) = -zrprn(jl)  / ztmst
        qfre(jl,jk,jrow) = -zfrln(jl)  / ztmst
        qacc(jl,jk,jrow) = -zsacln(jl) / ztmst
        cloud_tm1(jl,jk,jrow) = paclc(jl,jk)

        IF (zxlb(jl)>zeps.AND.zcdnc(jl,jk).GE.zcdnmin)THEN
           cdnc_acc(jl,jk,jrow)   = zcdnc(jl,jk)

           !--- In-cloud CDNC burden:

           zcdnc_burden(jl) = zcdnc_burden(jl)+zcdnc(jl,jk)*zdz(jl)

           !---- CDNC and burden averaged over cloudy and cloud-free periods

           cdnc(jl,jk,jrow) = zcdnc(jl,jk)*zclcaux(jl)
           cdnc_burden(jl,jrow) = cdnc_burden(jl,jrow) + &
                                  zcdnc(jl,jk) * zdz(jl) * zclcaux(jl)

           !--- In-cloud effective radius [um]:

           zkap=0.00045_dp*zcdnc(jl,jk)*1.e-6_dp+1.18_dp
           zreffl=1.E6_dp*zkap*((3._dp/(4._dp*api*rhoh2o))*zxlb(jl)*zrho(jl,jk)/zcdnc(jl,jk))**(1._dp/3._dp)

           reffl(jl,jk,jrow) = zreffl

        END IF

        IF (zxib(jl)>zeps)THEN
           icnc_acc(jl,jk,jrow) =zicnc(jl,jk)
           !--- Cirrus radii [um]:
           reffi(jl,jk,jrow)=zri(jl,jk)

           IF (nicnc > 0) THEN
              IF (zicnc(jl,jk) .GE. zicemin) THEN
                 !--- In-cloud ICNC burden:
                 zicnc_burden(jl)=zicnc_burden(jl)+zicnc(jl,jk)*zdz(jl)
                 !---- ICNC and burden averaged over cloudy and cloud-free periods
                 icnc(jl,jk,jrow) = zicnc(jl,jk)*zclcaux(jl)
                 icnc_burden(jl,jrow)=icnc_burden(jl,jrow) + &
                   zicnc(jl,jk)*zdz(jl)*zclcaux(jl)
              ENDIF
           ELSE
              zicnc_burden(jl)=0._dp
              icnc(jl,jk,jrow) =0._dp
              icnc_burden(jl,jrow)=0._dp
           ENDIF
        END IF

!--- End included for CDNC/IC scheme -----------------------------------

!       8.4   Corrections: Avoid negative cloud water/ice
!
        zxlold = zxlp1
        lo             = (zxlp1 .LT. ccwmin)
        zxlp1          = MERGE(0.0_dp,zxlp1,lo)
        zdxlcor        = (zxlp1-zxlold)/ztmst
        zcdnold        = zcdnc(jl,jk)
        zcdnp          = MERGE(0.0_dp,zcdnc(jl,jk),lo)
        zdnlcor        = (zcdnp-zcdnold)/ztmst

        zxiold         = zxip1
        lo1            = (zxip1 .LT. ccwmin)
        zxip1          = MERGE(0.0_dp,zxip1,lo1)
        zdxicor        = (zxip1-zxiold)/ztmst
        IF (nicnc > 0) THEN
           zicnold             = zicnc(jl,jk)
           zicnp               = MERGE(0.0_dp,zicnc(jl,jk),lo1)
           zdnicor             = (zicnp-zicnold)/ztmst
        ENDIF
        paclc(jl,jk)   = MERGE(0.0_dp,paclc(jl,jk),lo.AND.lo1)
        paclcac(jl,jk) = paclcac(jl,jk)+paclc(jl,jk)*zdtime
        pxlte(jl,jk)   = pxlte(jl,jk)+zdxlcor
        pxite(jl,jk)   = pxite(jl,jk)+zdxicor
        pqte(jl,jk)    = pqte(jl,jk)-zdxlcor-zdxicor
        ptte(jl,jk)    = ptte(jl,jk)+zlvdcp(jl)*zdxlcor                &
                                    +zlsdcp(jl)*zdxicor

        pxtte(jl,jk,idt_cdnc) = pxtte(jl,jk,idt_cdnc) + &
                                zdnlcor /( zrho(jl,jk) * 1000._dp / M_air)
        IF (nicnc > 0) pxtte(jl,jk,idt_icnc) = pxtte(jl,jk,idt_icnc) + &
                                zdnicor /( zrho(jl,jk) * 1000._dp / M_air)
!
811  END DO
831 END DO    ! Vertical loop
!
!--- Included for prognostic CDNC/IC scheme ----------------------------
  DO jl=1, kproma
     IF (zcdnc_burden(jl)>zepsec) THEN 
        cdnc_burden_acc(jl,jrow) = zcdnc_burden(jl)
     END IF
     IF (zicnc_burden(jl)>zepsec .AND. nicnc > 0) THEN 
        icnc_burden_acc(jl,jrow) = zicnc_burden(jl)
     END IF
  END DO 
!--- End included for CDNC/IC scheme -----------------------------------
!

!     ------------------------------------------------------------------
!
!       10.    Diagnostics
!
!       10.1   Accumulated precipitation at the surface
!
  DO 911 jl    = 1,kproma
     prsfl(jl) = zrfl(jl)
     pssfl(jl) = zsfl(jl)
     paprl(jl) = paprl(jl)+zdtime*(prsfl(jl)+pssfl(jl))
     paprs(jl) = paprs(jl)+zdtime*pssfl(jl)
911 END DO
!
!       10.2   Total cloud cover
!
    DO 921 jl    = 1,kproma
      zclcov(jl) = 1.0_dp-paclc(jl,1)
921 END DO
    DO 923 jk      = 2,klev
      DO 922 jl    = 1,kproma
        zclcov(jl) = zclcov(jl)*(1._dp-MAX(paclc(jl,jk),paclc(jl,jk-1)))  &
                               /(1._dp-MIN(paclc(jl,jk-1),zxsec))
922   END DO
923 END DO
    DO 924 jl     = 1,kproma
      zclcov(jl)  = 1.0_dp-zclcov(jl)
      paclcov(jl) = paclcov(jl)+zdtime*zclcov(jl)
924 END DO
!
!       10.3   Vertical integrals of humidity, cloud water and cloud ice
!
    DO 931 jl   = 1,kproma
      zqvi(jl)  = 0.0_dp
      zxlvi(jl) = 0.0_dp
      zxivi(jl) = 0.0_dp
931 END DO
!
    DO 933 jk     = ktdia,klev
      DO 932 jl   = 1,kproma
        zdpg      = (paphm1(jl,jk+1)-paphm1(jl,jk))/g
        zqvi(jl)  = zqvi(jl)+pqm1(jl,jk)*zdpg
        zxlvi(jl) = zxlvi(jl)+pxlm1(jl,jk)*zdpg
        zxivi(jl) = zxivi(jl)+pxim1(jl,jk)*zdpg
932   END DO
933 END DO
!
    DO 934 jl   = 1,kproma
      pqvi(jl)  = pqvi(jl)+zdtime*zqvi(jl)
      pxlvi(jl) = pxlvi(jl)+zdtime*zxlvi(jl)
      pxivi(jl) = pxivi(jl)+zdtime*zxivi(jl)
934 END DO
!
  RETURN
END SUBROUTINE cloud_cdnc_icnc

!===============================================================================


!=======================================================
! EXTRA FUNCTIONS AND SUBROUTINES FROM MO_CIRRUS

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
     Y1   = 1.0_dp / Y
     SQY1 = SQRT(Y1)
     IF (Y.LE.0.2_dp) THEN
        PPROD = 1.0_dp
        POLY  = 1.0_dp
        DO K = 1, 4
           PPROD = PPROD * FLOAT(2*K-1)
           POLY  = POLY  + PPROD * (-0.5_dp*Y)**K
        END DO
        XEERFC = POLY / SQY1 / SQPI
     ELSEIF (Y.GE.2.0_dp) THEN
        PPROD  = 1.0_dp
        POLY   = 1.0_dp
        DO K = 1, 4
           PPROD = PPROD * FLOAT(2*K+1)
           POLY  = POLY  + FLOAT(2**K) * Y1**K / PPROD
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
     ! ***** B.KÄRCHER AND S.SOLOMON, JGR 104(D22), 27441-27459, 1999

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

    REAL(dp), INTENT(out) :: ZNICEX(kbdim)    ! number of ice nuclei
    REAL(dp), INTENT(out) :: ZRI(kbdim,KLEV)  ! 

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
       znicex(jl)=0._dp

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
             SL    = ( T - 170.0_dp ) / FLOAT(KMAX)
             TEMP = T  - SL * FLOAT(K-1)
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
          SL   = ( T - 170.0_dp ) / FLOAT(KMAX)
          K    = KLIST(JL)
          TEMP = T  - SL * FLOAT(K-1)
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
       ZNICEX(JL) = ZNICE*1.0E6_dp

    END DO ! IX-loop

    RETURN

  END SUBROUTINE XFRZMSTR

! ------------------------------------------------------------------------------

  SUBROUTINE XFRZHET (NOSIZE,nfrzmod,DT,C,R,SIG,P,T,V,SI,CI,RI,COOLR, PW, &
                      SCR, TEMP, PICE, PWCR)

     ! ***** HETEROGENEOUS FREEZING (ADIABATIC ASCENT)
     !
     ! ***** BERND KÄRCHER  APR 14  2003
     ! ***** bernd.kaercher@dlr.de  http://www.op.dlr.de/~pa1c/
     !
     ! ***** Ulrike Lohmann      Dalhousie Univ   Apr 14 03
     !       Ulrike.Lohmann@Dal.Ca
     ! ***** Johannes Hendricks  DLR-IPA          Apr 14 03
     !       Johannes.Hendricks@dlr.de
     !
     ! ***** References
     !       Kärcher, B. and U. Lohmann
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
        SL    = ( T - 170.0_dp ) / FLOAT(KMAX)
        DO K = 1, KMAX+1
           TEMP = T  - SL * FLOAT(K-1)
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
       R0(I)     = RMIN * VRAT**( THIRD*FLOAT(I-1) )
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
     ! ***** BERND KÄRCHER  APR 30  2002
     ! ***** bernd.kaercher@dlr.de  http://www.op.dlr.de/~pa1c/
     ! ***** Ulrike Lohmann  Dalhousie Univ  Apr 02  Ulrike.Lohmann@Dal.Ca
     !
     ! ***** References
     ! Kärcher, B. and U. Lohmann
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
        R0(I)     = RMIN * VRAT**( THIRD*FLOAT(I-1) )
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
 

END MODULE MESSY_CLOUD_Lohmann07


