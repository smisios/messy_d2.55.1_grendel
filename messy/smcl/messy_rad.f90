! **********************************************************************
MODULE messy_rad
! **********************************************************************

  USE messy_main_constants_mem, ONLY: DP, cemiss
  USE messy_rad_albedo,         ONLY: albedos_rad   &
                                    , su_albedo_rad &
                                    , calbmns, calbmxs, calbmni
  USE messy_rad_short_cmn,      ONLY: nsw, rsun_scale, rad_sw_initialize
  USE messy_rad_long,           ONLY: jpband, rad_lw_initialize
  USE messy_rad_long,           ONLY: &
         absa_1, absb_1, fracrefa_1, fracrefb_1, forref_1, selfref_1, &
         absa_2, absb_2, fracrefa_2, fracrefb_2, forref_2, selfref_2, &
         refparam_2, &
         absa_3, absb_3, fracrefa_3, fracrefb_3, forref_3, selfref_3, &
         absn2oa_3, &
         absn2ob_3, etaref_3, h2oref_3, n2oref_3, co2ref_3, strrat_3, &
         absa_4, absb_4, fracrefa_4, fracrefb_4, selfref_4, strrat1_4, &
         strrat2_4, &
         absa_5, absb_5, ccl4_5, fracrefa_5, fracrefb_5, selfref_5, &
         strrat1_5, strrat2_5, &
         absa_6, absco2_6, cfc11adj_6, cfc12_6, fracrefa_6, selfref_6, &
         absa_7, absb_7, absco2_7, fracrefa_7, fracrefb_7, selfref_7, &
         strrat_7, &
         absa_8, absb_8, fracrefa_8, fracrefb_8, selfref_8, absco2a_8, &
         absco2b_8, &
         absn2oa_8, absn2ob_8, cfc12_8, cfc22adj_8, h2oref_8, n2oref_8, &
         o3ref_8, &
         absa_9, absb_9, fracrefa_9, fracrefb_9, selfref_9, absn2o_9, &
         ch4ref_9, &
         etaref_9, h2oref_9, n2oref_9, strrat_9, &
         absa_10, absb_10, fracrefa_10, fracrefb_10, &
         absa_11, absb_11, fracrefa_11, selfref_11, fracrefb_11, &
         absa_12, fracrefa_12, selfref_12, strrat_12, &
         absa_13, fracrefa_13, selfref_13, strrat_13, absa_14, absb_14, &
         fracrefa_14, fracrefb_14, selfref_14, &
         absa_15, fracrefa_15, selfref_15, strrat_15, &
         absa_16, fracrefa_16, selfref_16, strrat_16,        &
         NGC, NGS, NGM, NGN, NGB, WT,                        &
         corr1, corr2,                                       &
         PREF, PREFLOG, TREF,                                &
         NG, NSPA, NSPB, WAVENUM1, WAVENUM2, DELWAVE,        &
         TOTPLNK, TOTPLK16

  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC :: DP
  PUBLIC :: albedos_rad
  PUBLIC :: su_albedo_rad
  PUBLIC :: calbmns, calbmxs, calbmni
  PUBLIC :: nsw, rsun_scale, rad_sw_initialize
  PUBLIC :: jpband, rad_lw_initialize
  PUBLIC :: &
         absa_1, absb_1, fracrefa_1, fracrefb_1, forref_1, selfref_1, &
         absa_2, absb_2, fracrefa_2, fracrefb_2, forref_2, selfref_2, &
         refparam_2, &
         absa_3, absb_3, fracrefa_3, fracrefb_3, forref_3, selfref_3, &
         absn2oa_3, &
         absn2ob_3, etaref_3, h2oref_3, n2oref_3, co2ref_3, strrat_3, &
         absa_4, absb_4, fracrefa_4, fracrefb_4, selfref_4, strrat1_4, &
         strrat2_4, &
         absa_5, absb_5, ccl4_5, fracrefa_5, fracrefb_5, selfref_5, &
         strrat1_5, strrat2_5, &
         absa_6, absco2_6, cfc11adj_6, cfc12_6, fracrefa_6, selfref_6, &
         absa_7, absb_7, absco2_7, fracrefa_7, fracrefb_7, selfref_7, &
         strrat_7, &
         absa_8, absb_8, fracrefa_8, fracrefb_8, selfref_8, absco2a_8, &
         absco2b_8, &
         absn2oa_8, absn2ob_8, cfc12_8, cfc22adj_8, h2oref_8, n2oref_8, &
         o3ref_8, &
         absa_9, absb_9, fracrefa_9, fracrefb_9, selfref_9, absn2o_9, &
         ch4ref_9, &
         etaref_9, h2oref_9, n2oref_9, strrat_9, &
         absa_10, absb_10, fracrefa_10, fracrefb_10, &
         absa_11, absb_11, fracrefa_11, selfref_11, fracrefb_11, &
         absa_12, fracrefa_12, selfref_12, strrat_12, &
         absa_13, fracrefa_13, selfref_13, strrat_13, absa_14, absb_14, &
         fracrefa_14, fracrefb_14, selfref_14, &
         absa_15, fracrefa_15, selfref_15, strrat_15, &
         absa_16, fracrefa_16, selfref_16, strrat_16,        &
         NGC, NGS, NGM, NGN, NGB, WT,                        &
         corr1, corr2,                                       &
         PREF, PREFLOG, TREF,                                &
         NG, NSPA, NSPB, WAVENUM1, WAVENUM2, DELWAVE,        &
         TOTPLNK, TOTPLK16

  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modstr = 'rad'
  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modver = '2.2'

  TYPE t_rad_work
       REAL(DP),POINTER                         :: ptr0 => NULL()
       REAL(DP),POINTER,   DIMENSION(:,:)       :: ptr2 => NULL()
       REAL(DP),POINTER,   DIMENSION(:,:,:)     :: ptr3 => NULL()
       REAL(DP),POINTER,   DIMENSION(:,:,:,:)   :: ptr4 => NULL()
  END TYPE t_rad_work
  PUBLIC :: t_rad_work

  INTEGER, PARAMETER, PUBLIC :: NRADCALL=50    ! max. no. of radiation calls
  LOGICAL, DIMENSION(NRADCALL),PUBLIC :: l_switch = .FALSE.

  PUBLIC :: rad_rad_smcl
  PUBLIC :: rad_radheat_smcl
  !PRIVATE :: vrev1,vrev2

  CONTAINS

    ! =======================================================================
    SUBROUTINE rad_rad_smcl(                 &
         ! Input
         & kproma,kbdim,klev,                & ! dimensions
         & i_sw,                             & ! which SW scheme
         & ndaylen,                          & !
         & psct, pmu0, palb,                 & ! solar irradiation, zenith angle, surface albedo
         & ppf, pph,                         & ! p at full levels, half levels
         & ptf, pth,                         & ! T at full levels, half levels
         & pq,pqs,pclfr,                     & ! spec. hum., saturation spec. hum., cld fraction
         & ptclc,icldlyr,                    & ! total cld cover, clear-sky/cloudy idx
         & po3, pco2, pch4, pn2o,            & ! O3, CO2, CH4, N2O
         & pcfcs,                            & ! CFC species (vol.mix.ratio)
         & paot_lw, paot_sw,                 & ! aeosol opt. prop.
         & pomega_sw, pgamma_sw,             & ! aeosol opt. prop.
         & pcldtau_lw,pcldtau_sw,            & ! cloud optical prop.
         & pcld_gamma,pcld_omega,            & ! cloud optical prop.
         ! Output
         ! LW   , SW   , LWclear, SWclear
         & pflt , pfls , pfltc  , pflsc,     & ! net flux/transmissivity profiles
         & pflnir, pflnirc,                  & ! net NIR transmissivity profiles
         & pflsw1, pflsw1c,                  & ! net UVVIS transmissivity profiles
         ! LW(up), SW(up), LWclear(up), SWclear(up)
         & pflut , ptrus, pflutc, ptrusc,  & ! upward flux/transmissivity profiles
         ! NIR(up), NIRclear(up)
         & ptruni, ptrunic &
         &)


      USE messy_main_constants_mem,  ONLY: g, avo => N_A  &
                                         , MH, MN, MO, MC, M_O3, M_H2O &
                                         , M_air, api=>pi

      USE messy_rad_long,            ONLY: rad_lon_RRTM_RRTM_140GP &
                                         , jpinpx & ! number of gases in WKL
                                         , jpxsec   ! number of tracer cross
      !                                             !  sections
      USE messy_rad_short_v1,           ONLY: rad_sw_SW_v1 => rad_sw_SW
      USE messy_rad_short_v2,           ONLY: rad_sw_SW_v2 => rad_sw_SW

      IMPLICIT NONE
      INTRINSIC :: EPSILON, INT, MAX

      ! Input
      ! vector length
      INTEGER,INTENT(in)                           :: kproma
      ! first dimension of 2D arrays
      INTEGER,INTENT(in)                           :: kbdim
      ! number of levels
      INTEGER,INTENT(in)                           :: klev
      ! which SW scheme
      INTEGER,INTENT(in)                           :: i_sw
      INTEGER,INTENT(in)                           :: ndaylen
      ! for rf-diagn.
      !
      ! solar irradiation in W/m2
      REAL(dp),INTENT(in)                          :: psct
      ! mu0 for solar zenith angle
      REAL(dp),INTENT(in), DIMENSION(kbdim)        :: pmu0
      ! surface albedo
      REAL(dp),INTENT(in), DIMENSION(kbdim)        :: palb
      ! full level pressure in Pa
      REAL(dp),INTENT(in), DIMENSION(kbdim,klev)   :: ppf
      ! half level pressure in Pa
      REAL(dp),INTENT(in), DIMENSION(kbdim,klev+1) :: pph
      ! full level temperature in K
      REAL(dp),INTENT(in), DIMENSION(kbdim,klev)   :: ptf
      ! half level temperature in K
      REAL(dp),INTENT(in), DIMENSION(kbdim,klev+1) :: pth
      ! specific humidity in g/g
      REAL(dp),INTENT(in), DIMENSION(kbdim,klev)   :: pq
      ! satur. specific humidity in g/g
      REAL(dp),INTENT(in), DIMENSION(kbdim,klev)   :: pqs
      ! fractional cloud cover in m2/m2
      REAL(dp),INTENT(in), DIMENSION(kbdim,klev)   :: pclfr
      ! index for clear or cloudy layers:
      REAL(dp),INTENT(in), DIMENSION(kbdim,klev)   :: icldlyr
      ! total cloud cover in m2/m2
      REAL(dp),INTENT(in), DIMENSION(kbdim)        :: ptclc
      ! o3 mass mixing ratio
      REAL(dp),INTENT(in), DIMENSION(kbdim,klev)   :: po3
      ! co2 mass mixing ratio:
      REAL(dp),INTENT(in), DIMENSION(kbdim,klev)   :: pco2
      ! ch4 mass mixing ratio
      REAL(dp),INTENT(in), DIMENSION(kbdim,klev)   :: pch4
      ! n2o mass mixing ratio
      REAL(dp),INTENT(in), DIMENSION(kbdim,klev)   :: pn2o
      ! cfc volume mixing ratio
      REAL(dp),INTENT(in), DIMENSION(kbdim,klev,2) :: pcfcs
      ! aerosol optical thickness lw
      REAL(dp),INTENT(in), DIMENSION(kbdim,klev,jpband):: paot_lw
      ! aerosol optical thickness sw
      REAL(dp),INTENT(in), DIMENSION(kbdim,klev,nsw)   :: paot_sw
      ! single scattering albedo (aerosol)
      REAL(dp),INTENT(in), DIMENSION(kbdim,klev,nsw)   :: pomega_sw
      ! assymetry factor (aerosol)
      REAL(dp),INTENT(in), DIMENSION(kbdim,klev,nsw)   :: pgamma_sw
      ! cloud emissivity
      REAL(dp),INTENT(in), DIMENSION(kbdim,klev,jpband):: pcldtau_lw
      ! extinction tau (cloud)
      REAL(dp),INTENT(in), DIMENSION(kbdim,klev,nsw)   :: pcldtau_sw
      ! single scattering albedo (cloud)
      REAL(dp),INTENT(in), DIMENSION(kbdim,klev,nsw)   :: pcld_gamma
      ! assymetry factor (cloud)
      REAL(dp),INTENT(in), DIMENSION(kbdim,klev,nsw)   :: pcld_omega

      ! OUT
      ! net flux profile, LW all sky
      REAL(dp),INTENT(out),DIMENSION(kbdim,klev+1) :: pflt
      ! net transmissivity profile, SW all sky
      REAL(dp),INTENT(out),DIMENSION(kbdim,klev+1) :: pfls
      ! net flux profile, LW clear sky
      REAL(dp),INTENT(out),DIMENSION(kbdim,klev+1) :: pfltc
      ! net transmissivity profile, SW clear sky
      REAL(dp),INTENT(out),DIMENSION(kbdim,klev+1) :: pflsc
      ! fb_mk_20111115+
      ! net transmissivity profile, NIR total sky:
      REAL(dp),INTENT(out),DIMENSION(kbdim,klev+1):: pflnir
      ! net transmissivity profile, NIR clear sky:
      REAL(dp),INTENT(out),DIMENSION(kbdim,klev+1):: pflnirc
      ! net transmissivity profile, UVVIS total sky:
      REAL(dp),INTENT(out),DIMENSION(kbdim,klev+1):: pflsw1
      ! net transmissivity profile, UVVIS clear sky:
      REAL(dp),INTENT(out),DIMENSION(kbdim,klev+1):: pflsw1c
      ! fb_mk_20111115-
      ! op_mk_20180103+
      ! upward flux profile, LW all sky
      REAL(dp),INTENT(out),DIMENSION(kbdim,klev+1) :: pflut
      ! upward transmissivity profile, SW all sky
      REAL(dp),INTENT(out),DIMENSION(kbdim,klev+1) :: ptrus
      ! upward flux profile, LW clear sky
      REAL(dp),INTENT(out),DIMENSION(kbdim,klev+1) :: pflutc
      ! upward transmissivity profile, SW clear sky
      REAL(dp),INTENT(out),DIMENSION(kbdim,klev+1) :: ptrusc
      ! upward transmissivity profile, NIR total sky:
      REAL(dp),INTENT(out),DIMENSION(kbdim,klev+1) :: ptruni
      ! upward transmissivity profile, NIR clear sky:
      REAL(dp),INTENT(out),DIMENSION(kbdim,klev+1) :: ptrunic
      ! op_mk_20180103-

      !=======================================
      !        local
      !======================================
      INTEGER::jp
      ! -------------- RRTM input -------------------------------
      ! full level pressure [mb]
      REAL(dp), DIMENSION(kbdim,klev)       :: zpave
      ! half level pressure [mb]
      REAL(dp), DIMENSION(kbdim,klev+1)     :: zpmb
      ! full level temperature [K]
      REAL(dp), DIMENSION(kbdim,klev)       :: ztave
      ! half level temperature [K]
      REAL(dp), DIMENSION(kbdim,klev+1)     :: ztl
      ! pressure thickness in Pa
      REAL(dp), DIMENSION(kbdim,klev)       :: zdp
      ! surface temp.
      REAL(dp), DIMENSION(kbdim)            :: ztbound
      ! surface emissivity in each band []:
      REAL(dp), DIMENSION(kbdim,jpband)     :: zsemiss
      ! secure cloud fraction
      REAL(dp), DIMENSION(kbdim,klev)       :: zclfr
      ! number of molecules/cm2 of dry air and water vapor [#/cm2]
      REAL(dp), DIMENSION(kbdim,klev)       :: zcoldry
      ! index for clear or cloudy layers
      INTEGER , DIMENSION(kbdim,klev)       :: zcldidx
      ! aerosol optical prop. lw
      REAL(dp), DIMENSION(kbdim,klev,jpband):: zaot_lw
      ! cloud optical prop. lw
      REAL(dp), DIMENSION(kbdim,klev,jpband):: zcldtau_lw

      ! number of molecules/cm2 of N species [#/cm2],
      ! N=JPINPX
      !  1: H2O
      !  2: CO2
      !  3: O3
      !  4: N2O
      !  5: - empty -
      !  6: CH4
      !  7... : - empty -
      REAL(dp), DIMENSION(kbdim,jpinpx,klev):: zwkl

      ! number of molecules/cm2 of N species in [1e20/cm2],
      ! N=JPXSEC
      !  1: - empty -
      !  2: CFC11
      !  3: CFC12
      !  4... : - empty -
      REAL(dp), DIMENSION(kbdim,jpxsec,klev):: zwx

      ! for LW, RRTM
!!$   ! molecular weight of moist air
!!$   REAL(dp), DIMENSION(kbdim,klev)       :: amm
      REAL(dp), PARAMETER                   :: amco2 = MC+2._dp*MO   &
                                             , amch4 = MC+4._dp*MH   &
                                             , amn2o = 2._dp*MN+MO
      ! to convert fluxes to W/m2
      REAL(dp)                              :: zfluxfac
      ! clear sky fraction of total column
      REAL(dp),   DIMENSION(kbdim)          :: ztclear
      ! layer cloud fract. with respect to cloud fract. of total column
      REAL(dp),   DIMENSION(kbdim,klev)     :: zcldfrac

      ! --------   RRTM output -----------------------------
      ! upward flux, total sky
      REAL(dp), DIMENSION(kbdim,klev+1)     :: ztotuflux
      ! upward flux, clear sky
      REAL(dp), DIMENSION(kbdim,klev+1)     :: ztotufluc
      ! downward flux, total sky
      REAL(dp), DIMENSION(kbdim,klev+1)     :: ztotdflux
      ! downward flux, clear sky
      REAL(dp), DIMENSION(kbdim,klev+1)     :: ztotdfluc
      ! surface emissivity
      REAL(dp), DIMENSION(kbdim)            :: zsemit

      !------------- SW F&B ---------------------------------------------------
      INTEGER                          :: jb,jk,jl
      ! -----------  SW F&B input ---------------------------------------------
      INTEGER                          :: kidia   ! index of first longitude
      INTEGER                          :: kfdia   ! index of last longitude
      REAL(dp)                         :: zsct    ! solar irradiation in W/m2
      REAL(dp)                         :: zcardi  ! co2 mass mixing ratio
      REAL(dp), DIMENSION(kbdim)       :: zpsol   ! half level pressure in Pa
      REAL(dp), DIMENSION(kbdim,nsw)   :: zalbd   ! diffuse surface albedo
      REAL(dp), DIMENSION(kbdim,nsw)   :: zalbp   ! parallel surface albedo
      REAL(dp), DIMENSION(kbdim,klev)  :: zwv     ! specific humidity in g/g
      REAL(dp), DIMENSION(kbdim,klev)  :: zqs     ! satur. spec. hum. in g/g
      REAL(dp), DIMENSION(kbdim)       :: zrmu0   ! mu0 for solar zenith angle
      ! cld asymmetry factor (vrev2, as from cloudopt)
      REAL(dp), DIMENSION(kbdim,nsw,klev) :: zcg, zcg_o
      ! cld single scattering albedo
      REAL(dp), DIMENSION(kbdim,nsw,klev) :: zomega, zomega_o
      REAL(dp), DIMENSION(kbdim,nsw,klev) :: ztau, ztau_o ! cld extincion
      REAL(dp), DIMENSION(kbdim,klev,nsw) :: zaot_sw
      ! aerosol asymmetry factor
      REAL(dp), DIMENSION(kbdim,klev,nsw) :: zgamma_sw
      ! aerosol single scattering albedo
      REAL(dp), DIMENSION(kbdim,klev,nsw) :: zomega_sw
      ! ozone mass mixing ratio *
      REAL(dp), DIMENSION(kbdim,klev)     :: zoz

      ! --------------SW output------------------
      REAL(dp), DIMENSION(kbdim,klev+1) :: zfsdwn  ! U ! downward flux total sky
      REAL(dp), DIMENSION(kbdim,klev+1) :: zfsup   ! U ! upward flux total sky
      REAL(dp), DIMENSION(kbdim,klev+1) :: zfcdwn  ! downward flux clear sky
      REAL(dp), DIMENSION(kbdim,klev+1) :: zfcup   ! upward flux clear sky

      REAL(dp), DIMENSION(kbdim,klev+1) :: zfnird   ! NIR flux total sky, down
      REAL(dp), DIMENSION(kbdim,klev+1) :: zfniru   ! NIR flux total sky, up
      REAL(dp), DIMENSION(kbdim,klev+1) :: zfcnird  ! NIR flux clear sky, down
      REAL(dp), DIMENSION(kbdim,klev+1) :: zfcniru  ! NIR flux clear sky, up
      REAL(dp), DIMENSION(kbdim,klev+1) :: zfsw1d   ! UVVIS flux total sky, down
      REAL(dp), DIMENSION(kbdim,klev+1) :: zfsw1u   ! UVVIS flux total sky, up
      REAL(dp), DIMENSION(kbdim,klev+1) :: zfcsw1d  ! UVVIS flux clear sky, down
      REAL(dp), DIMENSION(kbdim,klev+1) :: zfcsw1u  ! UVVIS flux clear sky, up

      ! ====================================
      ! 0. Initialisation for kproma < kbdim
      ! ====================================
      pflt(:,:) = 0.0_dp
      pfls(:,:) = 0.0_dp
      pfltc(:,:) = 0.0_dp
      pflsc(:,:) = 0.0_dp
      pflnir(:,:) = 0.0_dp
      pflnirc(:,:) = 0.0_dp
      pflsw1(:,:) = 0.0_dp
      pflsw1c(:,:) = 0.0_dp

      pflut(:,:) = 0.0_dp
      pflutc(:,:) = 0.0_dp
      ptrus(:,:) = 0.0_dp
      ptrusc(:,:) = 0.0_dp
      ptruni(:,:) = 0.0_dp
      ptrunic(:,:) = 0.0_dp

      ! ========================
      !1. General preprocessing
      ! ========================
      ! pressure
      zpave(1:kproma,:) =vrev2(ppf(1:kproma,:))/100._dp
      zpmb(1:kproma,:)  =vrev2(pph(1:kproma,:))/100._dp

      ! Pressure thickness in Pa
      zdp(1:kproma,:)=vrev2(pph(1:kproma,2:klev+1)-pph(1:kproma,1:klev))

      ! temperature
      ztave(1:kproma,:)=vrev2(ptf(1:kproma,:))
      ztl(1:kproma,:)  =vrev2(pth(1:kproma,:))

      ! secure cloud fraction
      zclfr(1:kproma,:)=vrev2(pclfr(1:kproma,:))
      zclfr(1:kproma,:)=MAX(zclfr(1:kproma,:),EPSILON(1.0_dp))

      ! ========================
      ! 2. LW computation
      ! =======================
      !
      ! 2.2 LW, RRTM
      ! ------------
      ! 2.2.1 preprocessing

      ! surface temperature
      ztbound(1:kproma)=pth(1:kproma,klev+1) !!

      ! surface emissivity
      zsemiss(1:kproma,1:5)= cemiss
      zsemiss(1:kproma,6:8)= cemiss ! atmospheric window !
      zsemiss(1:kproma,9:16)=cemiss
      !
      DO jp=1,jpband
         zaot_lw(1:kproma,:,jp) = vrev2(paot_lw(1:kproma,:,jp))
         zcldtau_lw(1:kproma,:,jp) = vrev2(pcldtau_lw(1:kproma,:,jp))
      END DO

      zcldidx(1:kproma,:)=INT(vrev2(icldlyr(1:kproma,1:klev)))

      ! number of molecules of dry air (incl. water vapour) and JPINPX species
      zwkl(1:kproma,:,:)=0._dp
! bug-fix+
!!$   zwkl(1:kproma,1,:)=vrev2(pq(1:kproma,:))*M_air/M_H2O
      ! conversion of spec. humidity into vmr
      zwkl(1:kproma,1,:)=vrev2(pq(1:kproma,:) / &
           (1.0_dp-pq(1:kproma,:)))*M_air/M_H2O
! bug-fix-
      zwkl(1:kproma,2,:)=vrev2(pco2(1:kproma,:))*M_air/amco2
      zwkl(1:kproma,3,:)=vrev2(po3(1:kproma,:))*M_air/M_O3
      zwkl(1:kproma,4,:)=vrev2(pn2o(1:kproma,:))*M_air/amn2o
      zwkl(1:kproma,6,:)=vrev2(pch4(1:kproma,:))*M_air/amch4

! bug-fix+
!!$   amm(1:kproma,:)=(1.0_dp-zwkl(1:kproma,1,:))*M_air+zwkl(1:kproma,1,:)*M_H2O
!!$   zcoldry(1:kproma,:)= (zdp(1:kproma,:)/100.0_dp)    &
!!$        *10.0_dp*avo/g/amm(1:kproma,:) &
!!$        /(1.0_dp+zwkl(1:kproma,1,:))
      ! dry air molecules per cm^2
      zcoldry(1:kproma,:) = zdp(1:kproma,:) / g * avo / M_air &
           * (1.0_dp - vrev2(pq(1:kproma,:)) ) / 10.0_dp
! bug-fix-

      DO jp=1,6
         zwkl(1:kproma,jp,:)=zcoldry(1:kproma,:)*zwkl(1:kproma,jp,:)
      END DO

      ! number of molecules of JPXSEC species, pcfcs is in vol.mixing ratio
      zwx(1:kproma,:,:)=0._dp
      zwx(1:kproma,2,:)=zcoldry(1:kproma,:)*vrev2(pcfcs(1:kproma,:,1))*1.e-20_dp
      zwx(1:kproma,3,:)=zcoldry(1:kproma,:)*vrev2(pcfcs(1:kproma,:,2))*1.e-20_dp

      ! 2.2.2 Call RRTM ----------------------
      !
      CALL rad_lon_RRTM_RRTM_140GP(kproma, kbdim, klev,& !--- input -----
           & zpave, ztave, ztl,                        &
           & ztbound, zsemiss,                         &
           & zclfr,zcldidx,                            &
           & zcoldry, zwkl, zwx,                       &
           & zcldtau_lw, zaot_lw,                      &
           & ztotuflux, ztotufluc,                     & !--- output ----
           & ztotdflux, ztotdfluc,                     &
           & zsemit )

      ! 2.2.3 Postprocessing -------------------
      !
      ! convert W/cm2/sr to W/m2
      zfluxfac  =2._dp*api*10000._dp
      pflt(1:kproma,:) =zfluxfac*vrev2(ztotdflux(1:kproma,:) &
           - ztotuflux(1:kproma,:)) ! net downward flux, total sky
      pfltc(1:kproma,:)=zfluxfac*vrev2(ztotdfluc(1:kproma,:) &
           - ztotufluc(1:kproma,:)) ! net downward flux, clear sky
      ! upward LW flux, total sky
      pflut(1:kproma,:) =zfluxfac*vrev2(-ztotuflux(1:kproma,:))
      ! upward LW flux, clear sky
      pflutc(1:kproma,:)=zfluxfac*vrev2(-ztotufluc(1:kproma,:))

      ! ========================
      ! 3. SW computation
      ! =======================
      !
      !  SW, F&B(3&1), 3 bands near IR, 1 band vis
      !  ---------------------------------------------

      ! 3.1.1 preprocess

      ! longitude loop start and end
      kidia=1
      kfdia=kproma

      ! solar irradiation at Top of atmosphere
      zsct=psct

      ! co2 mass mixing ration, uniformly mixed
      zcardi=pco2(1,1)

      ! surface albedo
      DO jb=1,nsw
         zalbd(1:kproma,jb)=palb(1:kproma)
         zalbp(1:kproma,jb)=palb(1:kproma)
      END DO

      ! water vapor specific humidity
      zwv(1:kproma,:)=pq(1:kproma,:)

      ! saturation specific humidity
      zqs(1:kproma,:)=pqs(1:kproma,:)

      ! zenith angle
      zrmu0(1:kproma)=pmu0(1:kproma)

      ! Pressure thickness in Pa
      zdp(1:kproma,:)=pph(1:kproma,2:klev+1)-pph(1:kproma,1:klev)

      ! surface pressure in Pa
      zpsol(1:kproma)=pph(1:kproma,klev+1)

      ! cloud cover
      ztclear(1:kproma)=1.0_dp-ptclc(1:kproma)
      DO jl=1,kproma
         IF (ptclc(jl)<EPSILON(1._dp)) THEN
            zcldfrac(jl,:)=0._dp
         ELSE
            zcldfrac(jl,:)=zclfr(jl,:)/ptclc(jl)
         END IF
      END DO

      ! ozon mass mixing ratio * pressure thickness of layer
      zoz(1:kproma,:)=vrev2(po3(1:kproma,:)*zdp(1:kproma,:))*46.6968_dp/g

      ! aerosol optical properties
      DO jb=1,nsw
         zaot_sw(1:kproma,:,jb)   = vrev2(paot_sw(1:kproma,:,jb))
         zomega_sw(1:kproma,:,jb) = vrev2(pomega_sw(1:kproma,:,jb))
         zgamma_sw(1:kproma,:,jb) = vrev2(pgamma_sw(1:kproma,:,jb))
      END DO

      ! cloud optical prop.
      ! resort index klev und nsw for sw-Input -> ztau(1:kproma,1:nsw,1:klev)
      DO jk=1,klev
         DO jb=1,nsw
            ztau_o(1:kproma,jb,jk)=pcldtau_sw(1:kproma,jk,jb)
            zomega_o(1:kproma,jb,jk)=pcld_omega(1:kproma,jk,jb)
            zcg_o(1:kproma,jb,jk)=pcld_gamma(1:kproma,jk,jb)
         END DO
      END DO

      DO jb=1,nsw
         ztau(1:kproma,jb,:) = vrev2(ztau_o(1:kproma,jb,:))
         zomega(1:kproma,jb,:) = vrev2(zomega_o(1:kproma,jb,:))
         zcg(1:kproma,jb,:) = vrev2(zcg_o(1:kproma,jb,:))
      END DO

      SELECT CASE(i_sw)
      CASE(1)
         CALL rad_sw_SW_v1 ( ndaylen                                     &
           !----input---------
           &, kidia, kfdia, kbdim, klev                                  &
           &, zgamma_sw, zomega_sw, zaot_sw                              &
           &, zsct, pco2, zpsol, zalbd, zalbp, zwv, zqs                  &
           &, zrmu0, zcg, ztclear, zcldfrac, zdp, zomega, zoz, zpmb      &
           &, ztau, ztave                                                &
           !----output--------
           &, zfsdwn, zfsup &
           &, zfcdwn, zfcup &
           &, zfnird, zfniru, zfcnird, zfcniru & ! NIR, up, down
           &, zfsw1d, zfsw1u, zfcsw1d, zfcsw1u & ! UVVIS, up, down
           ! TOA and surface diagnostic variables are obsolete.
           )
      CASE(2)
         CALL rad_sw_SW_v2 ( ndaylen                                     &
           !----input---------
           &, kidia, kfdia, kbdim, klev                                  &
           &, zgamma_sw, zomega_sw, zaot_sw                              &
           &, zsct, pco2, zpsol, zalbd, zalbp, zwv, zqs                  &
           &, zrmu0, zcg, ztclear, zcldfrac, zdp, zomega, zoz, zpmb      &
           &, ztau, ztave                                                &
           !----output--------
           &, zfsdwn, zfsup &
           &, zfcdwn, zfcup &
           &, zfnird, zfniru, zfcnird, zfcniru & !  NIR, up, down
           &, zfsw1d, zfsw1u, zfcsw1d, zfcsw1u & !  UVVIS, up, down
           ! TOA and surface diagnostic variables are obsolete.
           )
      END SELECT

      ! 3.1.3 Postprocessing

      ! net SW transmissivity, total sky, top to surface
      pfls(1:kproma,:) =vrev2(zfsdwn(1:kproma,:)-zfsup(1:kproma,:))
      ! net SW transmissivity, clear sky, top to surface
      pflsc(1:kproma,:)=vrev2(zfcdwn(1:kproma,:)-zfcup(1:kproma,:))

      ! net NIR transmissivity, total sky, top to surface
      pflnir(1:kproma,:) =vrev2(zfnird(1:kproma,:) -zfniru(1:kproma,:))
      ! net NIR transmissivity, clear sky, top to surface
      pflnirc(1:kproma,:)=vrev2(zfcnird(1:kproma,:)-zfcniru(1:kproma,:))
      ! net UVVIS transmissivity, total sky, top to surface
      pflsw1(1:kproma,:) =vrev2(zfsw1d(1:kproma,:) -zfsw1u(1:kproma,:))
      ! net UVVIS transmissivity, clear sky, top to surface
      pflsw1c(1:kproma,:)=vrev2(zfcsw1d(1:kproma,:)-zfcsw1u(1:kproma,:))

      ! upward transmissivity, total sky, top to surface
      ptrus(1:kproma,:) =vrev2(-zfsup(1:kproma,:))
      ! upward transmissivity, clear sky, top to surface
      ptrusc(1:kproma,:)=vrev2(-zfcup(1:kproma,:))
      ! upward NIR transmissivity, total sky, top to surface
      ptruni(1:kproma,:) =vrev2(-zfniru(1:kproma,:))
      ! upward NIR transmissivity, clear sky, top to surface
      ptrunic(1:kproma,:)=vrev2(-zfcniru(1:kproma,:))

    END SUBROUTINE rad_rad_smcl
    ! =======================================================================

    ! =======================================================================
    SUBROUTINE rad_radheat_smcl (kproma, kbdim, klev,  klevp1      &
         , pi0                                        &
         , ptm1,         pqm1                         &
         , ptrsof,       ptrsol                       &
         , ptrsw1,       ptrnir                       &
         , ptrs1f,       ptrnif                       &
         , ptrfll,       ptrflw,      ptrfli          &
         , psofll,       psoflw,      psofli          &
         , ptrfllac,     ptrflwac,    ptrfliac        &
         , psofllac,     psoflwac,    psofliac        &
         , psrad0u,      psradsu                      &
         , ptradsu                                    &
         , ptslm1,       ptsi,        ptsw            &
         , palbedo                                    &
         , palsol,       palsow,      palsoi          &
         , paphm1,       papm1                        &
         , ptslnew,      ptte                         &
         , pfrl,         pfrw,        pfri            &
         , pflxs,        pflxt                        &
         , pflxsf,       pflxtf                       &
         , pheatsw,      pheatlw                      &
         , pflxs1,       pflxni                       &
         , pheats1,      pheatni                      &
         , pheatswc,     pheatlwc                     &
         , ledith                                     &
         )

      !
      !**** *RADHEAT* - COMPUTES TEMPERATURE CHANGES DUE TO RADIATION.
      !
      !
      !     SUBJECT.
      !     --------
      !
      !          THIS ROUTINE COMPUTES THE TENDENCIES OF THE ATMOSPHERE'S
      !     TEMPERATURE DUE TO THE EFFECTS OF LONG WAVE AND SHORT WAVE
      !     RADIATION. THE COMPUTATION IS DONE ON THE T-1 TIME LEVEL USING
      !     VALUES OF ATMOSPHERIC TRANSMISIVITIES AND EMISSIVITIES THAT HAVE
      !     BEEN STORED AT THE LAST FULL RADIATION TIME STEP. THE SURFACE
      !     SOLAR FLUX LATER TO BE USED IN THE SOIL PROCESS CALCULATIONS IS
      !     ALSO STORED.
      !
      !**   INTERFACE.
      !     ----------
      !
      !          *RADHEAT* IS CALLED FROM *PHYSC*.
      !
      !     INPUT ARGUMENTS.
      !     ----- ---------
      !
      !
      !     OUTPUT ARGUMENTS.
      !     ------ ---------
      !
      !
      !     METHOD.
      !     -------
      !
      !     PRODUCT OF SOLAR
      !     INFLUX BY TRANSMISSIVITIES LEADS TO SOLAR FLUXES. THEN THE
      !     TEMPERATURES ARE INTERPOLATED/EXTRAPOLATED TO THE LAYER BOUNDARIES
      !     (AT THE BOTTOM ONE TAKES THE SURFACE TEMPERATURE) AND A PRODUCT BY
      !     EMISSIVITIES OF SIGMA*T**4 GIVES THERMAL FLUXES. THE TWO FLUXES
      !     ARE ADDED AND DIVERGENCES COMPUTED TO GIVE HEATING RATES.
      !
      !     EXTERNALS.
      !     ----------
      !
      !          *SOLANG*.
      !
      !     AUTHOR.
      !     ------
      !
      !     U. SCHLESE    DKRZ-HAMBURG    JUNE 1995

      !     Modifications
      !     U. Schlese, December 1999:  version for coupling
      !     U. Schlese, July 2000, *solang* removed, weighted surface fluxes
      !     I. Kirchner, May 2002, tendency diagnose bugfix for surface fluxes

      USE messy_main_constants_mem, ONLY : stbo, cpd => cp_air &
           , vtmpc2, g

      IMPLICIT NONE

      INTEGER,INTENT(IN)  :: kproma,kbdim,klev,klevp1
      REAL(dp),INTENT(IN) :: pi0(kbdim),                                    &
           ptm1(kbdim,klev),      pqm1(kbdim,klev),                         &
           ptrsol(kbdim,klevp1),  ptrsof(kbdim,klevp1),                     &
           ptrsw1(kbdim,klevp1),  ptrnir(kbdim,klevp1),                     &
           ptrs1f(kbdim,klevp1),  ptrnif(kbdim,klevp1),                     &
           ptslm1(kbdim),         ptsi(kbdim),            ptsw(kbdim),      &
           paphm1(kbdim,klevp1),  papm1(kbdim,klev),                        &
           ptslnew(kbdim),                                                  &
           pfrl(kbdim),           pfrw(kbdim),            pfri(kbdim),      &
           palbedo(kbdim),                                                  &
           palsol(kbdim),         palsow(kbdim),          palsoi(kbdim)

      REAL(dp),INTENT(OUT) :: ptte(kbdim,klev),                             &
           psrad0u(kbdim),        psradsu(kbdim),                           &
           ptradsu(kbdim),                                                  &
           psofllac(kbdim),       psoflwac(kbdim),        psofliac(kbdim),  &
           ptrfllac(kbdim),       ptrflwac(kbdim),        ptrfliac(kbdim)

      REAL(dp),INTENT(OUT) :: ptrfll(kbdim), ptrflw(kbdim), ptrfli(kbdim),  &
           psofll(kbdim), psoflw(kbdim), psofli(kbdim)

      ! pemter and pemtef are equal to pflxt and pflxtf and therefore
      ! obsolete, as they do not contain the emissivities.
      ! pflxt and pflxtf are declared as INTENT(IN OUT).
      REAL(dp), INTENT(IN OUT) :: pflxt(kbdim,klevp1), pflxtf(kbdim,klevp1)
      REAL(dp), INTENT(OUT) :: pflxs(kbdim,klevp1)
      REAL(dp), INTENT(OUT) :: pflxsf(kbdim,klevp1)
      !
      !  pflxs1  - UV-Vis fluxes
      !  pflxni  - NIR - fluxes
      REAL(dp), DIMENSION(kbdim,klevp1), INTENT(OUT) :: pflxs1, pflxni

      !
      !  pheatlw  - long-wave heating rates in K/s
      !  pheatsw  - short-wave heating rates in K/s
      !  pheats1  - UV-Vis heating rates in K/s
      !  pheatni  - NIR heating rates in K/s
      REAL(dp),INTENT(out), DIMENSION(kbdim,klev) :: pheatlw
      REAL(dp),INTENT(out), DIMENSION(kbdim,klev) :: pheatsw
      REAL(dp),INTENT(out), DIMENSION(kbdim,klev) :: pheats1
      REAL(dp),INTENT(out), DIMENSION(kbdim,klev) :: pheatni

      !  pheatlwc  - long-wave heating rates clear sky in K/s
      !  pheatswc  - short-wave heating rates clear sky in K/s
      REAL(dp),INTENT(out), DIMENSION(kbdim,klev) :: pheatlwc
      REAL(dp),INTENT(out), DIMENSION(kbdim,klev) :: pheatswc
      LOGICAL, INTENT(in)                         :: ledith

      !local
      INTEGER :: jk, jl
      REAL(dp) :: zcons3, zdtdt, zffact, zflbot, zfltop,  &
           zsr0u, zsrsu, ztrdown, ztrsu, ztrdownf

      REAL(dp) :: zti(kbdim,klevp1),   ztsnew(kbdim),       zteffl4(kbdim)
      !
      ! ----------------------------------------------------------------------
      !
      !*     1.   COMPUTATIONAL CONSTANTS.
      !           ------------- ----------
      !
      zcons3=g/cpd
      !
      ! INIT
      ptte(:,:) = 0.0_dp
      !     ------------------------------------------------------------------
      !
      !*         3.     TEMPERATURES AT LAYERS' BOUDARIES.
      !                 ------------ -- ------- ----------
      !
      !
      !*         3.1     INTERPOLATION PROPER.
      DO  jk=2,klev
         DO jl=1,kproma
            zti(jl,jk)=(ptm1(jl,jk-1)*papm1(jl,jk-1)                 &
                 *(papm1(jl,jk)-paphm1(jl,jk))                  &
                 +ptm1(jl,jk)*papm1(jl,jk)                      &
                 *(paphm1(jl,jk)-papm1(jl,jk-1)))               &
                 /(paphm1(jl,jk)*(papm1(jl,jk)-papm1(jl,jk-1)))
         END DO
      END DO
      !
      !*        3.2     SURFACE AND TOP OF ATMOSPHERE TEMPERATURE.
      !
      DO  jl=1,kproma
         !  -  fractional surface coverage:
         zti(jl,klevp1)=(pfrl(jl)*ptslm1(jl)**4                      &
              +pfri(jl)*ptsi(jl)**4                        &
              +pfrw(jl)*ptsw(jl)**4)**0.25_dp
         zteffl4(jl)=ptslm1(jl)**3*(4.0_dp*ptslnew(jl)-3.0_dp*ptslm1(jl))
         ztsnew(jl)=(pfrl(jl)*zteffl4(jl)                            &
              +pfri(jl)*ptsi(jl)**4                            &
              +pfrw(jl)*ptsw(jl)**4)**0.25_dp
         zti(jl,1)=ptm1(jl,1)-papm1(jl,1)*(ptm1(jl,1)-zti(jl,2))     &
              /(papm1(jl,1)-paphm1(jl,2))
      END DO

      !     ------------------------------------------------------------------
      !
      !*         4.    UPDATE FLUXES AND COMPUTE HEATING RATES.
      !

      DO jl=1,kproma
         !
         ! When ptrsol and ptrsof have been corrected with transmissivity
         ! from FUBRAD, they cannot be used here to calculate the heating
         ! rates. The separated transmissivities (UVVIS and NIR) will be
         ! used instead.
         !
         pflxs(jl,1)=pi0(jl)*(ptrsw1(jl,1)+ptrnir(jl,1))
         pflxsf(jl,1)=pi0(jl)*(ptrs1f(jl,1)+ptrnif(jl,1))

         pflxs1(jl,1)=pi0(jl)*ptrsw1(jl,1)
         pflxni(jl,1)=pi0(jl)*ptrnir(jl,1)

      END DO
      !
      !
      !     4.2  Fluxes and heating rates except for lowest layer
      !
      DO  jk=1,klev-1
         DO  jl=1,kproma
            zfltop=pflxs(jl,jk)+pflxt(jl,jk)

            pflxs(jl,jk+1)=pi0(jl)*(ptrsw1(jl,jk+1)+ptrnir(jl,jk+1))
            pflxsf(jl,jk+1)=pi0(jl)*(ptrs1f(jl,jk+1)+ptrnif(jl,jk+1))

            zflbot=pflxs(jl,jk+1)+pflxt(jl,jk+1)

            ! LW-"blending" from HAMMONIA
            IF (ledith) THEN
               IF ( (10._dp > papm1(jl,jk+1)) .AND. &
                    (papm1(jl,jk+1) > 2._dp) ) THEN
                  zflbot=pflxs(jl,jk+1)+pflxt(jl,jk+1)* &
                       (log(100000._dp/2._dp)-log(100000._dp/papm1(jl,jk+1)))/ &
                       (log(100000._dp/2._dp)-log(100000._dp/10._dp))
                  zfltop=pflxs(jl,jk)+pflxt(jl,jk)* &
                       (log(100000._dp/2._dp)-log(100000._dp/papm1(jl,jk+1)))/ &
                       (log(100000._dp/2._dp)-log(100000._dp/10._dp))
               END IF
            ENDIF

            zdtdt=-zcons3*(zflbot-zfltop)/((paphm1(jl,jk+1)            &
                 -paphm1(jl,jk))*(1.0_dp+vtmpc2*pqm1(jl,jk)))

            IF (.NOT. ledith) THEN
               ptte(jl,jk)=zdtdt
            ELSE

               IF (papm1(jl,jk) >= 1.0_dp) ptte(jl,jk)=zdtdt

            ENDIF
            !
            ! Update the UV-Vis and NIR fluxes
            !
            pflxs1(jl,jk+1) = pi0(jl)*ptrsw1(jl,jk+1)
            pflxni(jl,jk+1) = pi0(jl)*ptrnir(jl,jk+1)
            !
            ! Save the short-wave and long-wave heating rates
            !
            zffact = -zcons3 /((paphm1(jl,jk+1) - paphm1(jl,jk)) &
                 * (1.0_dp + vtmpc2*pqm1(jl,jk)) )
            pheatlw(jl,jk) = zffact * (pflxt(jl,jk+1) - pflxt(jl,jk))
            IF (.NOT. ledith) THEN
               pheatsw(jl,jk) = zffact * (pflxs(jl,jk+1) - pflxs(jl,jk))
            ELSE

               IF (papm1(jl,jk) >= 1.0_dp) THEN
                  pheatsw(jl,jk) = zffact * (pflxs(jl,jk+1) - pflxs(jl,jk))
               END IF

            END IF
            pheats1(jl,jk) = zffact * (pflxs1(jl,jk+1) - pflxs1(jl,jk))
            pheatni(jl,jk) = zffact * (pflxni(jl,jk+1) - pflxni(jl,jk))

            pheatlwc(jl,jk) = zffact * (pflxtf(jl,jk+1) - pflxtf(jl,jk))
            pheatswc(jl,jk) = zffact * (pflxsf(jl,jk+1) - pflxsf(jl,jk))

         END DO
      END DO
      !
      !
      !     4.3  Lowest layer
      !
      DO jl=1,kproma
         ! Calculation of ztrsu moved from 5.0 to here, to use the
         ! original value of pflxt.
         ztrsu = -cemiss*stbo*ztsnew(jl)**4+(cemiss-1.0_dp)               &
              *(stbo*zti(jl,klevp1)**4+pflxt(jl,klevp1)/cemiss)
         ptradsu(jl)= ztrsu
         ztrdown=pflxt(jl,klevp1)+cemiss*stbo*zti(jl,klevp1)**4
         pflxt(jl,klevp1)=ztrdown-cemiss*stbo*ztsnew(jl)**4
         ztrdownf=pflxtf(jl,klevp1)+cemiss*stbo*zti(jl,klevp1)**4
         pflxtf(jl,klevp1)=ztrdownf-cemiss*stbo*ztsnew(jl)**4
         ptrfll(jl)=ztrdown-cemiss*stbo*zteffl4(jl)
         ptrflw(jl)=ztrdown-cemiss*stbo*ptsw(jl)**4
         ptrfli(jl)=ztrdown-cemiss*stbo*ptsi(jl)**4
         pflxs(jl,klevp1)=pi0(jl)*(ptrsw1(jl,klevp1)+ptrnir(jl,klevp1))
         pflxsf(jl,klevp1)=pi0(jl)*(ptrs1f(jl,klevp1)+ptrnif(jl,klevp1))

         psofll(jl)=(1.0_dp-palsol(jl))*pflxs(jl,klevp1)/(1.0_dp-palbedo(jl))
         psoflw(jl)=(1.0_dp-palsow(jl))*pflxs(jl,klevp1)/(1.0_dp-palbedo(jl))
         psofli(jl)=(1.0_dp-palsoi(jl))*pflxs(jl,klevp1)/(1.0_dp-palbedo(jl))
         zfltop=pflxs(jl,klev)+pflxt(jl,klev)
         zflbot=pflxs(jl,klevp1)+pflxt(jl,klevp1)
         zdtdt=-zcons3*(zflbot-zfltop)/((paphm1(jl,klevp1)              &
              -paphm1(jl,klev))*(1.0_dp+vtmpc2*pqm1(jl,klev)))

         ptte(jl,klev)=zdtdt

         pflxs1(jl,klevp1) = pi0(jl)*ptrsw1(jl,klevp1)
         pflxni(jl,klevp1) = pi0(jl)*ptrnir(jl,klevp1)
         ! Save the short-wave and long-wave heating rates
         zffact = -zcons3 /((paphm1(jl,klevp1) - paphm1(jl,klev)) &
              * (1.0_dp + vtmpc2*pqm1(jl,klev)))
         pheatlw(jl,klev) = zffact * (pflxt(jl,klevp1) - pflxt(jl,klev))
         pheatsw(jl,klev) = zffact * (pflxs(jl,klevp1) - pflxs(jl,klev))
         pheats1(jl,klev) = zffact * (pflxs1(jl,klevp1) - pflxs1(jl,klev))
         pheatni(jl,klev) = zffact * (pflxni(jl,klevp1) - pflxni(jl,klev))

         pheatlwc(jl,klev) = zffact * (pflxtf(jl,klevp1) - pflxtf(jl,klev))
         pheatswc(jl,klev) = zffact * (pflxsf(jl,klevp1) - pflxsf(jl,klev))

      END DO

      !     ------------------------------------------------------------------
      !
      !*         5.     Diagnostics of additional fluxes
      !
      !NOTE: With new memory channel interface, laccu (requiring
      !       'user'-accumulation) is obsolete !mz_pj_20050525
      DO jl = 1, kproma
         zsr0u = pflxs(jl,1) - pi0(jl)
         psrad0u(jl) = zsr0u
         zsrsu = -pflxs(jl,klevp1)*(1.0_dp/(1.0_dp-palbedo(jl))-1.0_dp)
         psradsu(jl) = zsrsu
      END DO

      !        Area weighted surface fluxes over land/water/ice

      DO jl=1,kproma
         psofllac(jl)=psofll(jl)*pfrl(jl)
         psoflwac(jl)=psoflw(jl)*pfrw(jl)
         psofliac(jl)=psofli(jl)*pfri(jl)
         ptrfllac(jl)=ptrfll(jl)*pfrl(jl)
         ptrflwac(jl)=ptrflw(jl)*pfrw(jl)
         ptrfliac(jl)=ptrfli(jl)*pfri(jl)
      END DO

      RETURN

    END SUBROUTINE rad_radheat_smcl
    ! =======================================================================

    ! =======================================================================
    ! PRIVATE SUBROUTINES AND FUNCTIONS
    ! =======================================================================

    ! =======================================================================
    FUNCTION vrev1(p)

      IMPLICIT NONE
      INTRINSIC :: SIZE

      ! inverts array p(klev)

      REAL(dp), INTENT(in), DIMENSION(:) :: p
      REAL(dp), DIMENSION(SIZE(p))       :: vrev1

      INTEGER :: klev, jk

      klev=SIZE(p)

      DO jk=1,klev
         vrev1(jk)=p(klev+1-jk)
      END DO

    END FUNCTION vrev1
    ! =======================================================================

    ! =======================================================================
    FUNCTION vrev2(p)

      IMPLICIT NONE
      INTRINSIC :: SIZE

      ! inverts array p(kproma,klev) in the second index

      REAL(dp), INTENT(in), DIMENSION(:,:)     :: p
      REAL(dp), DIMENSION(SIZE(p,1),SIZE(p,2)) :: vrev2

      INTEGER :: klev, jk

      klev=SIZE(p,2)

      DO jk=1,klev
         vrev2(:,jk)=p(:,klev+1-jk)
      END DO

    END FUNCTION vrev2
    ! =======================================================================

! **********************************************************************
END MODULE messy_rad
! **********************************************************************
