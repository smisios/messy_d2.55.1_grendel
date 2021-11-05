#if defined (NAG)
#define ARGCHECK 1
#endif
!OCL NOALIAS

SUBROUTINE physc(krow, kglat)
!
! Description:
!
! Controls the calls to the various physical subroutines.
!
!
!  *physc* is called from *gpc*.
!
!  Externals:
!
!  *geopot*    computes full level geopotentials.
!  *pres*      computes half level pressure.
!  *presf*     computes full level pressure.
!  *radiation* controls radiation computations.
!  *radheat*   computes radiation tendencies.
!  *vdiff*     computes vertical exchange by turbulence.
!  *ssodrag*   computes gravity wave drag.
!  *cucall*    controls mass-flux scheme.
!  *cloud*     computes large scale water phase changes and cloud cover.
!  *surf*      computes soiltemperature and moisture.
!  *lake*      computes water temperature and ice of lakes
!  *ml_ocean*  computes mixed layer ocean
!  *licetemp*  computes ice temperature of lakes
!  *sicetemp*  computes sea ice skin temperature.
!
!
!  Authors:
!
!  M. Jarraud, ECMWF, January 1982, original source
!
!     Modifications.
!     --------------
!  M.A. Giorgetta, MPI-Hamburg, May 2000, modified for ECHAM5
!  A.Tompkins,     MPI-Hamburg, June  2000, cloud cover scheme
!  U. Schlese,     MPI-Hamburg, July  2000, calling sequence changed
!  I. Kirchner,    MPI, July 2000, tendency diagnostics revision
!  L. Kornblueh,   MPI, August 2001, changes for different advection
!                       schemes, adapt polefilter call 
!  U. Schulzweida, MPI, May 2002, blocking (nproma)
!  I. Kirchner,    MPI, August 2002, nmi revision/extension
!  U. Schlese, M. Esch, MPI, September 2002, mixed layer ocean
!  L. Kornblueh,   MPI, August 2004, new tropopause calculation
!
!
USE mo_kind,              ONLY: dp
USE mo_exception,         ONLY: message
USE mo_memory_g1a,        ONLY: xlm1, xim1, tm1, qm1, alpsm1, xtm1
USE mo_memory_g2a,        ONLY: vm1, um1
USE mo_memory_g3a
USE mo_memory_g3b
#ifdef OBSOLETE
USE mo_control,           ONLY: ltdiag, lcouple, lmidatm, lhd,   &
                                nlev, nlevp1, ldiagamip, ltimer
#else
USE mo_control,           ONLY: ltdiag, lcouple, lmidatm,   &
                                nlev, nlevp1, ltimer
#endif
USE mo_hyb,               ONLY: delb, nlevm1
#ifndef MESSY
USE mo_param_switches,    ONLY: lcond, lconv, lsurf, lgwdrag, lrad
USE mo_constants,         ONLY: cpd, vtmpc1, vtmpc2, g, tmelt
#else
!!$USE mo_param_switches,    ONLY: lgwdrag ! op_pj_20160618
USE mo_constants,         ONLY: cpd, vtmpc1, vtmpc2, g, tmelt          &
                              , rhoh2o ! mz_lg_20030902
#endif
#ifndef MESSY
USE mo_radiation,         ONLY: solc
USE mo_cloud,             ONLY: ctaus, ctaul, ctauk
USE mo_hydrology,         ONLY: hydrology_collect
#endif
USE mo_scan_buffer,       ONLY: vo, vol, vom, qte, xlte, xite, tte, xtte, &
                                alnpr, alpha, alpste, vervel
#ifndef MESSY
USE mo_physc1,            ONLY: cdisse
#endif
#ifdef OBSOLETE
USE mo_physc2,            ONLY: cqsncr, cwlmax
#endif
#ifndef MESSY
USE mo_tracer,            ONLY: ntrac, xtdriver1, xtdriver2, xtdiagn
#else
!mz_rs_20040326+
USE messy_main_tracer_mem_bi, ONLY: ntrac => ntrac_gp
!mz_rs_20040326-
#endif
!
#ifndef MESSY
USE mo_midatm,            ONLY: gwspectrum
#endif
#ifndef MESSY
USE mo_ssortns,           ONLY: ssodrag
#endif
!
USE mo_decomposition,     ONLY: ldc => local_decomposition
USE mo_diag_tendency,     ONLY: pdiga
#ifdef OBSOLETE
USE mo_column,            ONLY: get_col_pol, lcotra
#endif
#ifndef MESSY
USE mo_geoloc,            ONLY: amu0_x, rdayl_x, sqcst_2d, &
                                philat_2d, budw_2d, twomu_2d
!
USE mo_time_control,      ONLY: lstart, lresume, delta_time, l_trigrad,&
                                lfirst_day, time_step_len, l_getocean
#else
USE mo_geoloc,            ONLY: sqcst_2d, philat_2d
!
USE mo_time_control,      ONLY: lstart, lresume, delta_time,&
                                lfirst_day, time_step_len
#endif
!
USE mo_advection
USE mo_spitfire,          ONLY: pole_filter
#ifdef MESSY
! mz_rs_20020828+
USE messy_main_data_bi, ONLY: geopot_3d, geopoti_3d, &
!!$                             cvs, cvw & ! ub_ak_20190208
                                prc, prl,                        &
!!$                             icecov,  &                  ! mz_pj_20030714 ! ub_ak_20190208+
                                ilab, qtec, qflux, &        ! mz_ht_20040405
                                vmixtau, vdiffp, &          ! mz_ht_20041026
! ub_ak_20190208+
!!$                             cvsc,  seacov, &            ! mz_bs_20050520
!!$                             landcov, &                  ! op_pj_20130410
! ub_ak_20190208-
                                tslnew, &                   ! mz_pj_20050522
                                srfl_2d, &                  ! mz_pj_20071130
                                aphm1, aphp1, apm1, app1,              &
                                ! ub_ak_20190208+
!!$                             geom1, loglac,                         &
!!$                             loland, rsfc, rsfl, ssfc, ssfl,        &
                                !geom1, rsfc, rsfl, ssfc, ssfl,        &
                                ! ub_ak_20190208-
! mz_ap_20161103+
!!$                                ! fields for exchange with ocean:
!!$                                awhea_2d, aicon_2d, &  ! mz_ap_20070820
!!$                                aiqre_2d,           &  ! mz_ap_20070820
!!$                                awfre_2d, aifre_2d, &  ! mz_ap_20070820
!!$                                aiust_2d,aivst_2d,  &  ! mz_ap_20070820
!!$                                awust_2d, awvst_2d, &  ! mz_ap_20070820
!!$                                apmebco_2d, &          ! mz_ap_20090519
!!$                                evapi, cvsi, &         ! fb_mk_20100215
!!$                                evapw, &               ! mz_ho_20160412
!!$                                ahfsw, & !!$ ahfsi,  &       ! fb_mk_20100514
!!$                                radflxw_2d,         &  ! op_pj_20130408
!!$                                cvsi, & ! fb_mk_20100215 ! ub_ak_20190208
! mz_ap_20161103-
                                ! op_sb_20120927+
                                rsfl_2d, rsfc_2d, ssfl_2d, ssfc_2d !, & 
! ub_ak_20190208+
!!$                             wlmx_2d, & !!$ evapl_2d, evapot_2d,       &
!!$                             loland_2d, loglac_2d !,               & 
! ub_ak_20190208-
                                ! op_sb_20121119-
!!$                             ahflw, ustri, vstri, ustrw, vstrw ! op_pj_20160614
                              
! mz_rs_20020828-
#ifdef MESSYTENDENCY
   USE messy_main_tendency_bi, ONLY: mtend_set_sqcst_scal
   USE messy_main_mpi_bi,      ONLY: p_pe
#endif
#endif
USE mo_timer,             ONLY: timer_start, timer_stop, timer_radiation, timer_cloud
!
USE mo_nmi,               ONLY: nmi_phase, NMI_ACCU, NMI_USE_AVG, &
     dh_t, dh_m, dh_l, buf_t, buf_m, buf_l, lnmi_run, lnmi_cloud
#ifdef OBSOLETE
USE mo_diag_amip2,        ONLY: collect_amip2_diag 
#endif
#ifndef MESSY
USE mo_diag_radiation,    ONLY: diag_rad1, diag_rad2
USE mo_tropopause,        ONLY: WMO_tropopause
#endif
#ifdef MESSY
!mz_ht_20040415+
USE messy_main_switch,    ONLY: USE_CONVECT, USE_CLOUD
!!#D crm +
USE messy_main_switch,    ONLY: USE_CRM
!!#D crm -
!mz_ht_20040415-
#endif

IMPLICIT NONE
!
INTEGER :: krow, kglat
! 
! ub_ak_20190208+
!!$#ifndef OBSOLETE
!!$! op_pj_20160619: this has been temporarily moved here to
!!$!                 get rid of mo_physc2.f90, iniphy and finaly physc ...
!!$!  *inverse of equivalent water height when snow
!!$!  is considered to cover completely the ground
!!$!  in the box.
!!$REAL(dp), PARAMETER :: cqsncr = 0.95_dp   
!!$REAL(dp), PARAMETER :: cwlmax = 2.E-4_dp  !  *maximum moisture content of
!!$                                          !   the skin reservoir
!!$#endif
! ub_ak_20190208-
!  Local array bounds
INTEGER :: nglon, nproma, nbdim
!
!  Local scalars:
REAL(dp) :: zcst, zrcst, ztwodt, zprat,  zcdnc,  zn1,  zn2, ztsurf         &
       , zepsec, zsigfac, zsigh, zsn_mm
INTEGER :: jlev, jl, kfdia, kidia, ktdia, jk, nexp
INTEGER :: jglat
LOGICAL :: loconv, locond
! 
!  Local arrays:
#ifdef MESSY
REAL(dp) ::  zalpha_hp(ldc%nproma,nlev) ! mz_lg_20030127
REAL(dp) ::  zdadc(ldc%nproma), zdpsdt(ldc%nproma), zgeo(ldc%nproma)          &
        ,ztvm1(ldc%nproma,nlev)                                               &
#else
REAL(dp) ::  zdadc(ldc%nproma), zdpsdt(ldc%nproma), zgeo(ldc%nproma)          &
        ,ztvm1(ldc%nproma,nlev), ztslnew(ldc%nproma)                     &
#endif
! zbetaa: qt distribution minimum in beta
! zbetab: qt distribution maximum in beta
! zvdiffp:  dq/dt from vdiff scheme needed
! zhmixtau:  timescale of mixing for horizontal eddies
! zvmixtau:  timescale of mixing for vertical turbulence
        ,zbetaa(ldc%nproma,nlev)                                        &
        ,zbetab(ldc%nproma,nlev)
#ifndef MESSY
 REAL(dp) :: zvdiffp(ldc%nproma,nlev)                                   &
        ,zhmixtau(ldc%nproma,nlev)                                      &
        ,zvmixtau(ldc%nproma,nlev)                                      &
        ,zi0(ldc%nproma)
REAL(dp) ::  zqtec(ldc%nproma,nlev), zbetass(ldc%nproma,nlev)
#else
REAL(dp) ::  zbetass(ldc%nproma,nlev)
#endif
REAL(dp) ::  zqtold(ldc%nproma)

!
#ifndef MESSY
LOGICAL :: lonorth(ldc%nproma)
#endif
!
!  Surface fluxes over land/water/ice
!
#ifndef MESSY
REAL(dp) :: zhfsl(ldc%nproma),  zhfsw(ldc%nproma),  zhfsi(ldc%nproma),  &
        zhfll(ldc%nproma),  zhflw(ldc%nproma),  zhfli(ldc%nproma),  &
        zevapl(ldc%nproma), zevapw(ldc%nproma), zevapi(ldc%nproma), &
        ztrfll(ldc%nproma), ztrflw(ldc%nproma), ztrfli(ldc%nproma), &
        zsofll(ldc%nproma), zsoflw(ldc%nproma), zsofli(ldc%nproma)
#else
! fb_mk_20100215+
REAL(dp) :: zhfsl(ldc%nproma), zhfll(ldc%nproma)
! sensible/latent heat flux over water in W m-2
!!$REAL(dp), DIMENSION(:), POINTER :: zhfsw, zhflw ! mz_ap_20161103
! sensible/latent heat flux over ice in W m-2
!!$REAL(dp), DIMENSION(:), POINTER :: zhfsi, zhfli ! op_pj_20160614
! evaporation over ice
!!$REAL(dp), DIMENSION(:), POINTER :: zevapi ! mz_ap_20161103
! fb_mk_20100215-
!!$REAL(dp), DIMENSION(:), POINTER :: zevapl ! op_sb_20120927
!!$REAL(dp), DIMENSION(:), POINTER :: zevapw ! mz_ho_20160412 ! mz_ap_20161103
#endif

#ifndef MESSY
REAL(dp) ::  zfrl(ldc%nproma),  zfrw(ldc%nproma),   zfri(ldc%nproma)          &
        ,zcvsi(ldc%nproma), zcvsc(ldc%nproma),  zwlmx(ldc%nproma)         &
        ,zcvs(ldc%nproma),  zcvw(ldc%nproma)
#else
! ub_ak_20190208+
!!$! mz_bs_20050114+
!!$REAL(dp), DIMENSION(:), POINTER :: zfrl     ! op_pj_20130410
!!$REAL(dp), DIMENSION(:), POINTER :: zcvsc
!!$REAL(dp), DIMENSION(:), POINTER :: zfrw
!!$! mz_bs_20050114-
!!$REAL(dp), DIMENSION(:), POINTER :: zfri
!!$REAL(dp), DIMENSION(:), POINTER :: zcvs
!!$REAL(dp), DIMENSION(:), POINTER :: zcvw
!!$REAL(dp), DIMENSION(:), POINTER :: zcvsi ! fb_mk_20100215
!!$REAL(dp), DIMENSION(:), POINTER :: zwlmx ! op_sb_20121205
! ub_ak_20190208-
#endif
! 
!  Local arrays for the HD-model and glacier calving model
!

#ifndef MESSY     
REAL(dp) :: zros_hd(ldc%nproma), zdrain_hd(ldc%nproma)
REAL(dp) :: zalac(ldc%nproma)
#endif

#ifndef MESSY
INTEGER :: ilab(ldc%nproma,nlev), itype(ldc%nproma)
!!$#else
!!$INTEGER :: itype(ldc%nproma)
#endif
INTEGER :: invb(ldc%nproma)
!
INTEGER :: itrpwmo(ldc%nproma), itrpwmop1(ldc%nproma)
!
!
!    Arrays internal to physics
!
#ifndef MESSY
REAL(dp) :: qhfla(ldc%nproma), evapot(ldc%nproma), zprecip(ldc%nproma)
#else
REAL(dp) :: zprecip(ldc%nproma)
! mz_ht_20041026+
REAL(dp), DIMENSION(:,:), POINTER :: zilab    ! movable to cloud or convect?
REAL(dp), DIMENSION(:,:), POINTER :: zqtec    ! obsolete?
REAL(dp), DIMENSION(:,:), POINTER :: zvmixtau ! for vdiff 
REAL(dp), DIMENSION(:,:), POINTER :: zvdiffp  ! for vdiff
REAL(dp), DIMENSION(:),   POINTER :: qhfla    ! for vdiff
REAL(dp), DIMENSION(:),   POINTER :: ztslnew  ! for vdiff ! mz_pj_20050522
! mz_ht_20041026-
!!$REAL(dp), DIMENSION(:), POINTER :: evapot     ! for vdiff ! op_sb_20121205
#endif
!
REAL(dp) :: zdtime
!
! Local arrays
!
#ifndef MESSY
REAL(dp) :: srfl(ldc%nproma)
LOGICAL  :: loland(ldc%nproma), loglac(ldc%nproma)
REAL(dp) :: geom1(ldc%nproma,nlev)
REAL(dp) :: aphm1(ldc%nproma,nlevp1), apm1(ldc%nproma,nlev) 
REAL(dp) :: aphp1(ldc%nproma,nlevp1), app1(ldc%nproma,nlev) 
REAL(dp) :: rsfc(ldc%nproma), ssfc(ldc%nproma)
REAL(dp) :: rsfl(ldc%nproma), ssfl(ldc%nproma)
#else
REAL(dp), DIMENSION(:,:), POINTER :: geom1 => NULL()
REAL(dp), DIMENSION(:),   POINTER :: srfl  => NULL()
REAL(dp), DIMENSION(:),   POINTER :: rsfc  => NULL()
REAL(dp), DIMENSION(:),   POINTER :: ssfc  => NULL()
REAL(dp), DIMENSION(:),   POINTER :: ssfl  => NULL()
REAL(dp), DIMENSION(:),   POINTER :: rsfl  => NULL()
#endif
! 
!  External subroutines
#ifndef MESSY
EXTERNAL geopot, pres, presf, radiation, vdiff, surf, cloud, &
         sicetemp, cucall, radheat, collect, &
         lake, ml_ocean, licetemp ! op_pj_20130407
#else
! mz_jb_20040511+
EXTERNAL geopot, pres, presf, &
         messy_convec, messy_physc, messy_vdiff, &
         messy_mixlo, &        ! fb_mk_20110129
         messy_radiation, messy_radheat, messy_gwdrag ! op_pj_20130407
! mz_jb_20040511-
#endif
! 
!  Intrinsic functions
INTRINSIC EXP, &
  ABS, EPSILON, MAX, MIN, NINT, SQRT, TANH   ! mz_jb_20040512
!
!  Local array bounds
  jglat = kglat           ! global continuous latitude index
  nglon = ldc% nglon      ! local number of longitudes

  nbdim = ldc% nproma

  IF ( krow == ldc% ngpblks ) THEN
    nproma = ldc% npromz
  ELSE
    nproma = ldc% nproma
  END IF

#ifdef MESSY
! mz_ht_20040405+
  zilab    => ilab(:,:,krow)
  zqtec    => qtec(:,:,krow)
  zvmixtau => vmixtau(:,:,krow)
  zvdiffp  => vdiffp(:,:,krow)
  qhfla    => qflux(:,krow)
  zilab(:,:)    = 0.0_dp
  zqtec(:,:)    = 0.0_dp
  qhfla(:)      = 0.0_dp
  zvdiffp(:,:)  = 0.0_dp
  zvmixtau(:,:) = 0.0_dp
! mz_ht_20040405-
! ub_ak_20190208+
!!$ mz_pj_20041026+
!!$  zfri => icecov(:,krow)
!!$  zcvs => cvs(:,krow)
!!$  zcvw => cvw(:,krow)
!!$ mz_pj_20041026-
!!$! mz_bs_20050114+
!!$  zfrw => seacov(:,krow)
!!$  zcvsc => cvsc(:,krow)
!!$! mz_bs_20050114-
!!$  zfrl => landcov(:,krow)   ! op_pj_20130410
! ub_ak_20190208-
  ztslnew => tslnew(:,krow) ! mz_pj_20050522
  ! op_sb_20120927+ for surface
  rsfl  => rsfl_2d(:,krow)
  rsfc  => rsfc_2d(:,krow)
  ssfl  => ssfl_2d(:,krow)
  ssfc  => ssfc_2d(:,krow)
!!$  zwlmx  => wlmx_2d(:,krow)   ! ub_ak_20190208
!!$  zevapl => evapl_2d(:,krow)  ! op_pj_20160615
!!$  evapot => evapot_2d(:,krow) ! op_pj_20160615
! ub_ak_20190208+
!!$ loland    => loland_2d(:,krow)
!!$ loglac    => loglac_2d(:,krow)
!!$  ! op_sb_20120927-
!!$  ! fb_mk_20100212+
! ub_ak_20190208-
! mz_ap_20161103+
!!$  zhfsw  => ahfsw(:,krow)
!!$  zhflw  => ahflw(:,krow)
!!$!!!$  zhfsi  => ahfsi(:,krow) ! op_pj_20160614
!!$!!!$  zhfli  => ahfli(:,krow) ! op_pj_20160614
!!$  zevapi => evapi(:,krow)
!!$  zevapw => evapw(:,krow) ! mz_ho_20160412
! mz_ap_20161103-
!!$  zcvsi  => cvsi(:,krow) ! ub_ak_20190208
! mz_ap_20161103+
!!$  zhfsw(:)  = 0.0_dp
!!$  zhflw(:)  = 0.0_dp
! mz_ap_20161103-
!!$  zhfsi(:)  = 0.0_dp ! op_pj_20160614
!!$  zhfli(:)  = 0.0_dp ! op_pj_20160614
! mz_ap_20161103+
!!$  zevapi(:) = 0.0_dp
!!$  zevapw(:) = 0.0_dp      ! mz_ho_20160412
! mz_ap_20161103-
!!$  zcvsi(:)  = 0.0_dp ! ub_ak_20190208
  ! fb_mk_20100212-
  ! op_sb_20120927+
  rsfl(:) = 0.0_dp
  rsfc(:) = 0.0_dp
  ssfl(:) = 0.0_dp
  ssfc(:) = 0.0_dp
!!$  zwlmx(:) = 0.0_dp ! ub_ak_20190208
!!$  zevapl(:)= 0.0_dp
!!$  evapot(:) = 0.0_dp ! op_pj_20160615
  ! op_sb_20120927-

#endif

!*    COMPUTATIONAL CONSTANTS.
!     ------------- ----------
!
  zdtime = delta_time
  zepsec=1.E-12_dp
  zsigfac=0.15_dp
!
!     ----------------------------------------------------------------
!
!*        2.    ALLOCATE STORAGE 
!               -------- ------- 
!
200 CONTINUE
!
#ifdef MESSY
  ALLOCATE (aphm1(ldc%nproma,nlevp1))
  ALLOCATE (apm1(ldc%nproma,nlev))
  ALLOCATE (aphp1(ldc%nproma,nlevp1))
  ALLOCATE (app1(ldc%nproma,nlev))
  srfl => srfl_2d(:,krow)
  ! mz_pj_20050616+
  aphm1(:,:) = 0.0_dp
  apm1(:,:)  = 0.0_dp
  aphp1(:,:) = 0.0_dp
  app1(:,:)  = 0.0_dp
! ub_ak_20190208+
!!$ loland(:)  = .FALSE.
!!$ loglac(:)  = .FALSE.
! ub_ak_20190208-
  rsfc(:)    = 0.0_dp
  ssfc(:)    = 0.0_dp
  rsfl(:)    = 0.0_dp
  ssfl(:)    = 0.0_dp
  ! mz_pj_20050616-
#endif
!
!     ------------------------------------------------------------
!
!*        3.    COMPUTE SOME FIELDS NEEDED BY THE PHYSICAL ROUTINES.
!               ------- ---- ------ ------ -- --- -------- ---------
!
300 CONTINUE
!
!*        3.1   COMPUTE VIRTUAL TEMPERATURE AT T-DT AND SET *ZGEO* TO 0.
!
310 CONTINUE
!
  ztvm1(1:nproma,:) = tm1(1:nproma,:,krow)*(1._dp+vtmpc1*qm1(1:nproma,:,krow) &
                      -(xlm1(1:nproma,:,krow)+xim1(1:nproma,:,krow)))
!
  zgeo(1:nproma)=0._dp
!
!*        3.2   COMPUTE (PHI-PHIS) AT T-DT USING LN(P) AT T.
!
320 CONTINUE
!
#ifdef MESSY
  ! mz_lg_20030127+ 
  ! added the calculation of the half pressure
  ! geopotential height, which is needed to calculate the thickness
  ! the layers by calling the routine geopot and setting the
  ! parameter alpha to zero.
  zalpha_hp(:,:) = 0.0_dp
  geom1 => geopoti_3d(:,:,krow)
  CALL geopot(geom1,ztvm1,alnpr(:,:,krow),zalpha_hp(:,:),zgeo,nbdim,nproma)
  geom1 => geopot_3d(:,:,krow)
  ! mz_lg_20030127-
#endif
  CALL geopot(geom1,ztvm1,alnpr(:,:,krow),alpha(:,:,krow),zgeo,nbdim,nproma)
!
!
!*        3.3    COMPUTE PRESSURE AT FULL AND HALF LEVELS AT T-DT.
!
330 CONTINUE
!
  aphm1(1:nproma,nlevp1)=EXP(alpsm1(1:nproma,krow))
!
  CALL pres(aphm1,nbdim,aphm1(1,nlevp1),nproma)
!
  CALL presf(apm1,nbdim,aphm1,nproma)
!
!*        3.4   COMPUTE REAL WINDS AND WIND TENDENCIES.
!
340 CONTINUE
!
  DO jlev = 1, nlev
!DIR$ CONCURRENT
    DO jl = 1, nproma
      zrcst=1._dp/sqcst_2d(jl,krow)
      um1(jl,jlev,krow) = zrcst*um1(jl,jlev,krow)
      vm1(jl,jlev,krow) = zrcst*vm1(jl,jlev,krow)
      vol(jl,jlev,krow) = zrcst*vol(jl,jlev,krow)
      vom(jl,jlev,krow) = zrcst*vom(jl,jlev,krow)
    END DO
  END DO
#ifdef MESSYTENDENCY
  ! From here on horiz. wind vector not scaled with cos latitude
  call mtend_set_sqcst_scal(.false.)
#endif

!
! Horizontal wind shear for horizontal mixing of variance
!
  DO jlev=1,nlev
     DO jl=1,nproma
#ifndef MESSY       
        zhmixtau(jl,jlev)=ctauk*ABS(vo(jl,jlev,krow)) 
        zhmixtau(jl,jlev)=MIN(ctaus,MAX(ctaul,zhmixtau(jl,jlev)))
#endif        
        zvmixtau(jl,jlev)=0._dp
     END DO
  END DO
!
! ---------------------------!!!!!-----------------------------------
!   no moisture computations at the first day to prevent instability
!   at the beginning caused by unbalanced initial fields  !!!!
!
#ifndef MESSY
!     IF(lfirst_day) THEN
!       locond=.false.
!       loconv=.false.
!     ELSE
       locond=lcond
       loconv=lconv
!     END IF
#else
!mz_ht_20040422+
       IF (lstart) THEN
          locond=.false.
          loconv=.false.
!!#D crm +
          IF (USE_CRM) loconv = .TRUE.
!!#D crm -
       ELSE
          locond=USE_CLOUD
          loconv=USE_CONVECT
!!#D crm +
          loconv=USE_CONVECT .OR. USE_CRM
!!#D crm -
       END IF
!mz_ht_20040422-
#endif

     IF (lnmi_run) THEN
       ! control cloud parameterisation in nmi initialization mode
       locond=lnmi_cloud
       loconv=lnmi_cloud
       SELECT CASE(nmi_phase)
       CASE(NMI_ACCU)     ! store tendencies
!DIR$ CONCURRENT
         buf_t(:,:) = tte(:,:,krow)
!DIR$ CONCURRENT
         buf_m(:,:) = vom(:,:,krow)
!DIR$ CONCURRENT
         buf_l(:,:) = vol(:,:,krow)
       CASE(NMI_USE_AVG)  ! store tendencies
!DIR$ CONCURRENT
         buf_t(:,:) = tte(:,:,krow)
!DIR$ CONCURRENT
         buf_m(:,:) = vom(:,:,krow)
!DIR$ CONCURRENT
         buf_l(:,:) = vol(:,:,krow)
       CASE default
       END SELECT
     END IF
!
! ------------------------------------------------------------------
!
!*        3.5   ESTIMATE ADIABATIC CONVERSION OF POTENTIAL ENERGY.
!
350 CONTINUE
!
  zdadc(1:nproma)=0._dp
!
  DO 351 jl=1,nproma
     zdpsdt(jl)=aphm1(jl,nlevp1)*alpste(jl,krow)
     zdadc(jl)=zdadc(jl)+geospm(jl,krow)*zdpsdt(jl)
351 END DO
!
  DO 353 jlev=1,nlev
     DO 352 jl=1,nproma
        zdadc(jl)=zdadc(jl)+(1._dp+vtmpc2*qm1(jl,jlev,krow))*cpd*   &
                (tte(jl,jlev,krow)*(aphm1(jl,jlev+1)-aphm1(jl,jlev)) &
                        +tm1(jl,jlev,krow)*delb(jlev)*zdpsdt(jl))
352  END DO
353 END DO
!
#ifndef MESSY
!
!*        3.6   COMPUTE LOGICAL MASK FOR LAND AND GLACIER.
!
360 CONTINUE
!
  DO 365 jl=1,nproma
     loland(jl)=slm(jl,krow).GT.0._dp
     loglac(jl)=loland(jl).AND.glac(jl,krow).GT.0._dp
     lonorth(jl)=philat_2d(jl,krow).GT.0 ! true in northern hemisphere
365 END DO
!
!       3.7 Weighting factors for fractional surface coverage
!           Accumulate ice portion for diagnostics
!
!DIR$ CONCURRENT
   DO jl=1,nproma
      zfrl(jl)=slm(jl,krow)
      zfrw(jl)=(1._dp-slm(jl,krow))*(1._dp-seaice(jl,krow))
      zfri(jl)=1._dp-zfrl(jl)-zfrw(jl) 
!!$#ifndef MESSY
      friac(jl,krow)=friac(jl,krow)+zdtime*zfri(jl)
!!$#endif
   END DO
  IF(lcouple) THEN
    DO jl=1,nproma
       IF(slf(jl,krow).GT.1.0_dp-zepsec) THEN
         tsi(jl,krow)=tmelt
         tsw(jl,krow)=tmelt
       END IF
    END DO
  ENDIF
!
!      3.8  Skin reservoir, wet skin fraction and snow cover
!           (bare land, canopy, lake ice)
!
! op_pj_20160619: qqq this block should be moved to SURFACE?
   DO jl=1,nproma
     IF (.NOT.loglac(jl)) THEN
        zwlmx(jl)=cwlmax*(1._dp+vlt(jl,krow))
        zcvw(jl)=MIN(wl(jl,krow)/zwlmx(jl),1.0_dp)
        zsn_mm=1000._dp*sn(jl,krow)
        zsigh=SQRT(zsn_mm/(zsn_mm+zepsec+zsigfac*orostd(jl,krow)))
        zcvs(jl)=cqsncr*TANH(zsn_mm/10._dp)*zsigh
        zcvsc(jl)=MIN(1._dp,snc(jl,krow)/(zwlmx(jl)-cwlmax+EPSILON(1._dp)))
        IF (zcvs(jl).LT.EPSILON(1._dp) .AND. zcvsc(jl).GE.EPSILON(1._dp)) THEN
           zcvs(jl)=zcvsc(jl)
        END IF
     ELSE
        zwlmx(jl)=0._dp
        zcvw(jl)=0._dp
        zcvs(jl)=1._dp
        zcvsc(jl)=0._dp
     END IF
     IF (.NOT. loland(jl)) THEN
        zcvsi(jl)=TANH(sni(jl,krow)*100._dp)
     ELSE
        zcvsi(jl)=0._dp
     END IF
   END DO
#endif
!
370 CONTINUE
!
!
!*         3.8    SET LOOP  VALUES FOR PHYSICAL PARAMETERIZATIONS
!
  kidia=1
  kfdia=nproma
  ktdia=1
!
!         3.9   DETERMINE TROPOPAUSE HEIGHT AND MASS BUDGETS
!
!
#ifndef MESSY
  CALL WMO_tropopause (nproma, nbdim, nlev,                 &
                       tm1(:,:,krow),  apm1, tropo(:,krow), &
                       itrpwmo, itrpwmop1)
#endif
!
!
#ifndef MESSY
!*    3.12 INITIALISATION OF CLOUD DROPLET NUMBER CONCENTRATION 
!          (1/M**3) USED IN RADLSW AND CLOUD 
!
  IF (lstart) THEN
     DO 4103 jk=ktdia,nlev        
!DIR$ CONCURRENT
        DO 4102 jl=kidia,kfdia       
           nexp=2
           zprat=(MIN(8._dp,80000._dp/apm1(jl,jk)))**nexp
           IF (loland(jl).AND.(.NOT.loglac(jl))) THEN
              zn1= 50._dp
              zn2=220._dp
           ELSE 
              zn1= 50._dp
              zn2= 80._dp
           ENDIF
           IF (apm1(jl,jk).LT.80000._dp) THEN
              zcdnc=1.e6_dp*(zn1+(zn2-zn1)*(EXP(1._dp-zprat)))
           ELSE
              zcdnc=zn2*1.e6_dp
           ENDIF
           acdnc(jl,jk,krow)=zcdnc
           acdncm(jl,jk,krow)=zcdnc
4102    END DO
4103 END DO
  ENDIF
!
#ifndef MESSY
     DO 4104 jl=kidia,kfdia       
           itype(jl)=NINT(rtype(jl,krow))
4104 END DO
#endif
!
     DO 4106 jk=ktdia,nlev
        DO 4105 jl=kidia,kfdia
           zqtec(jl,jk)=0.0_dp
4105    END DO
4106 END DO
!
!*        3.13   DIAGNOSE CURRENT CLOUD COVER
!
  IF(locond) THEN

     IF (ltimer) CALL timer_start(timer_cloud)

#ifdef ARGCHECK
     CALL cover( nproma, nbdim, ktdia, nlev, nlevp1                   &
               , itype,             zfrw                              &
               , invb,              rintop(:,krow)                    &
               , aphm1,             apm1                              &
               , qm1(:,:,krow),     tm1(:,:,krow)                     &
               , xlm1(:,:,krow),    xim1(:,:,krow)                    &
               , vervel(:,:,krow)                                     &
               , xvar(:,:,krow),    xskew(:,:,krow)                   &
               , aclc(:,:,krow)                                       &
               , zbetaa,            zbetab                            &
               , zbetass                                              &
                )
#else
     CALL cover( nproma, nbdim, ktdia, nlev, nlevp1                   &
               , itype,             zfrw                              &
               , invb,              rintop(1,krow)                    &
               , aphm1,             apm1                              &
               , qm1(1,1,krow),     tm1(1,1,krow)                     &
               , xlm1(1,1,krow),    xim1(1,1,krow)                    &
               , vervel(1,1,krow)                                     &
               , xvar(1,1,krow),    xskew(1,1,krow)                   &
               , aclc(1,1,krow)                                       &
               , zbetaa,            zbetab                            &
               , zbetass                                              &
                )
#endif
     IF (ltimer) CALL timer_stop(timer_cloud)

  ENDIF
#else
  ! mz_ht_20041102+
  !  separate routine for cloud droplet numbers and cover from 
  !  MESSY cloud called in MESSY_radiation
  CALL messy_radiation
  ! mz_ht_20041102-
#endif
!
!*        4.    RADIATION PARAMETERISATION.
!               --------- -----------------
!
400 CONTINUE
!
!
410 CONTINUE
!
#ifndef MESSY
  IF (lrad) THEN

  IF (l_trigrad) THEN

    IF (ltimer) CALL timer_start(timer_radiation)

#ifdef ARGCHECK
    CALL radiation(nproma,nbdim,nlev,nlevp1,krow,jglat,         &
                   twomu_2d(:,krow),budw_2d(:,krow),            &
                   aphm1,apm1,tm1(:,:,krow),                    &
                   qm1(:,:,krow),                               &
                   acdnc(:,:,krow),aclc(:,:,krow),              &
                   xlm1(:,:,krow),xim1(:,:,krow),               &
                   loland,loglac,                               &
                   forest(:,krow),seaice(:,krow),               &
                   zcvs,sni(:,krow),                            &
                   zcvsc, vlt(:,krow),                          &
                   zfrl,zfri,zfrw,                              &
                   tslm1(:,krow),tsi(:,krow),tsw(:,krow),       &
                   alb(:,krow),ao3(:,:,krow),                   &
                   aclcv(:,krow),                               &
                   emter(:,:,krow),trsol(:,:,krow),             &
                   emtef(:,:,krow),trsof(:,:,krow),             &
                   emtef0(:,:,krow),trsof0(:,:,krow),           &
                   alsol(:,krow),alsoi(:,krow),                 &
                   alsow(:,krow),albedo(:,krow),                &
                   so4nat(:,:,krow),so4all(:,:,krow),           &
                   diag_rad1,diag_rad2)
#else
    CALL radiation(nproma,nbdim,nlev,nlevp1,krow,jglat,         &
                   twomu_2d(1,krow),budw_2d(1,krow),            &
                   aphm1,apm1,tm1(1,1,krow),                    &
                   qm1(1,1,krow),                               &
                   acdnc(1,1,krow),aclc(1,1,krow),              &
                   xlm1(1,1,krow),xim1(1,1,krow),               &
                   loland,loglac,                               &
                   forest(1,krow),seaice(1,krow),               &
                   zcvs,sni(1,krow),                            &
                   zcvsc, vlt(1,krow),                          &
                   zfrl,zfri,zfrw,                              &
                   tslm1(1,krow),tsi(1,krow),tsw(1,krow),       &
                   alb(1,krow),ao3(1,1,krow),                   &
                   aclcv(1,krow),                               &
                   emter(1,1,krow),trsol(1,1,krow),             &
                   emtef(1,1,krow),trsof(1,1,krow),             &
                   emtef0(1,1,krow),trsof0(1,1,krow),           &
                   alsol(1,krow),alsoi(1,krow),                 &
                   alsow(1,krow),albedo(1,krow),                &
                   so4nat(1,1,krow),so4all(1,1,krow),           &
                   diag_rad1,diag_rad2)
#endif

    IF (ltimer) CALL timer_stop(timer_radiation)

  END IF

  ELSE
     ! lrad=.FALSE.
     ! --> no radiative effect of the atmosphere or the surface
     emter(:,:,krow)=0._dp
     trsol(:,:,krow)=1._dp
     emtef(:,:,krow)=0._dp
     trsof(:,:,krow)=1._dp
     emtef0(:,:,krow)=0._dp
     trsof0(:,:,krow)=1._dp
     ! --> no radiative heating or cooling of 
     !     (1) the surface, as computed in vdiff, or
     !     (2) the atmosphere, as computed in radheat
     END IF
#endif
! ... of #ifndef MESSY ! op_pj_20130407

!
! ----------------------------------------------------------------------
!
!       Update solar incidence *zi0* and compute solar surface flux *srfl*
!     
!
#ifndef MESSY
  IF (lrad) THEN
     DO  jl=1,nproma
        zi0(jl)=cdisse*solc*amu0_x(jl,krow)*rdayl_x(jl,krow)
        srfl(jl)=zi0(jl)*trsol(jl,nlevp1,krow)
     END DO
  ELSE
     DO  jl=1,nproma
        zi0(jl)=0._dp
        srfl(jl)=0._dp
     END DO
  END IF
#endif
! ... of #ifndef MESSY ! op_pj_20130407
!
!     ------------------------------------------------------------
!
!*              VERTICAL EXCHANGE OF U,V,T,Q BY TURBULENCE.
!               -------- -------- -- - - - - -- -----------
!
!
!      COMPUTE PRESSURE AT FULL AND HALF LEVELS AT T+DT.
!
  ztwodt=time_step_len
  DO 522 jl=1,nproma
     aphp1(jl,nlevp1)=EXP(alpsm1(jl,krow)+ztwodt*alpste(jl,krow))
522 END DO
!
  CALL pres(aphp1,nbdim,aphp1(1,nlevp1),nproma)
!
  CALL presf(app1,nbdim,aphp1,nproma)
!
  IF (ltdiag) THEN
! prepare next fields for VDIFF
#ifndef MESSY
    pdiga(1:nglon,:, 3,krow) = pdiga(1:nglon,:, 3,krow) - vom(1:nproma,:,krow)
    pdiga(1:nglon,:, 8,krow) = pdiga(1:nglon,:, 8,krow) - vol(1:nproma,:,krow)
    pdiga(1:nglon,:,16,krow) = pdiga(1:nglon,:,16,krow) - tte(1:nproma,:,krow)
#else
    ! mz_jb_20040830+
    pdiga(1:nproma,:, 3,krow) = pdiga(1:nproma,:, 3,krow)- vom(1:nproma,:,krow)
    pdiga(1:nproma,:, 8,krow) = pdiga(1:nproma,:, 8,krow)- vol(1:nproma,:,krow)
    pdiga(1:nproma,:,16,krow) = pdiga(1:nproma,:,16,krow)- tte(1:nproma,:,krow)
    ! mz_jb_20040830-
#endif
  ENDIF
!!
! mz_ho_20160412+
#ifdef MESSY
  CALL messy_vdiff
#else
! mz_ho_20160412-
#ifdef ARGCHECK
  CALL vdiff (nproma, nbdim, ktdia, nlev, nlevm1, nlevp1, ntrac       &
            , krow                                                    &
!
            , xtm1(:,:,:,krow)                                        &
!
            , qm1(:,:,krow),     tm1(:,:,krow),     um1(:,:,krow)     &
            , vm1(:,:,krow),     xlm1(:,:,krow),    xim1(:,:,krow)    &
            , xvar(:,:,krow)                                          &
!
            , ahfl(:,krow),      ahfs(:,krow),      az0(:,krow)       &
            , dew2(:,krow),      evap(:,krow),      forest(:,krow)    &
            , temp2(:,krow),     t2max(:,krow)                        &
            , t2min(:,krow),     wind10w(:,krow),   vdis(:,krow)      &
            , u10(:,krow),       v10(:,krow),       ustr(:,krow)      &
            , vstr(:,krow),      wimax(:,krow),     wind10(:,krow)    &
            , wsmx(:,krow),      vlt(:,krow),       grndcapc(:,krow)  &
            , grndhflx(:,krow),  vgrat(:,krow)                        &
            , tsl(:,krow),       tsw(:,krow),       tsi(:,krow)       &
            , ocu(:,krow),       ocv(:,krow)                          &
            , az0l(:,krow),      az0w(:,krow),      az0i(:,krow)      &
            , zhfsl,             zhfsw,             zhfsi             &
            , zhfll,             zhflw,             zhfli             &
            , zevapl,            zevapw,            zevapi            &
            , ahfslac(:,krow),   ahfswac(:,krow),   ahfsiac(:,krow)   &
            , ahfllac(:,krow),   ahflwac(:,krow),   ahfliac(:,krow)   &
            , evaplac(:,krow),   evapwac(:,krow),   evapiac(:,krow)   &
            , ustrl(:,krow),     ustrw(:,krow),     ustri(:,krow)     &
            , vstrl(:,krow),     vstrw(:,krow),     vstri(:,krow)     &
            , sn(:,krow),        snc(:,krow),       tslm1(:,krow)     &
            , ws(:,krow),        albedo(:,krow),    alsol(:,krow)     &
!
            , tke(:,:,krow),     tkem1(:,:,krow),   tkem(:,:,krow)    &
            , aclc(:,:,krow),    emter(:,:,krow)                      &
!
            , aphm1,             apm1,              geom1             &
            , ztvm1                                                   &
!!$#ifndef MESSY
            , zvdiffp,           zvmixtau          &
!!$#endif
!
            , zcvs,              zcvw,              srfl              &
            , qhfla,             evapot                               &
            , ztslnew,           zwlmx                                &
            , zfrl,              zfrw,              zfri              &
            , loland,            loglac                               &
!
!!$#ifndef MESSY
            , xtte(:,:,:,krow)                                        &
!!$#endif
!
            , vol(:,:,krow),     vom(:,:,krow),     qte(:,:,krow)     &
            , tte(:,:,krow),     xlte(:,:,krow),    xite(:,:,krow))
#else
  CALL vdiff (nproma, nbdim, ktdia, nlev, nlevm1, nlevp1, ntrac       &
            , krow                                                    &
!
            , xtm1(:,:,:,krow)                                        &
!
            , qm1(1,1,krow),     tm1(1,1,krow),     um1(1,1,krow)     &
            , vm1(1,1,krow),     xlm1(1,1,krow),    xim1(1,1,krow)    &
            , xvar(1,1,krow)                                          &
!
            , ahfl(1,krow),      ahfs(1,krow),      az0(1,krow)       &
            , dew2(1,krow),      evap(1,krow),      forest(1,krow)    &
            , temp2(1,krow),     t2max(1,krow)                        &
            , t2min(1,krow),     wind10w(1,krow),   vdis(1,krow)      &
            , u10(1,krow),       v10(1,krow),       ustr(1,krow)      &
            , vstr(1,krow),      wimax(1,krow),     wind10(1,krow)    &
            , wsmx(1,krow),      vlt(1,krow),       grndcapc(1,krow)  &
            , grndhflx(1,krow),  vgrat(1,krow)                        &
            , tsl(1,krow),       tsw(1,krow),       tsi(1,krow)       &
            , ocu(1,krow),       ocv(1,krow)                          &
            , az0l(1,krow),      az0w(1,krow),      az0i(1,krow)      &
            , zhfsl,             zhfsw,             zhfsi             &
            , zhfll,             zhflw,             zhfli             &
            , zevapl,            zevapw,            zevapi            &
            , ahfslac(1,krow),   ahfswac(1,krow),   ahfsiac(1,krow)   &
            , ahfllac(1,krow),   ahflwac(1,krow),   ahfliac(1,krow)   &
            , evaplac(1,krow),   evapwac(1,krow),   evapiac(1,krow)   &
            , ustrl(1,krow),     ustrw(1,krow),     ustri(1,krow)     &
            , vstrl(1,krow),     vstrw(1,krow),     vstri(1,krow)     &
            , sn(1,krow),        snc(1,krow),       tslm1(1,krow)     &
            , ws(1,krow),        albedo(1,krow),    alsol(1,krow)     &
!
            , tke(1,1,krow),     tkem1(1,1,krow),   tkem(1,1,krow)    &
            , aclc(1,1,krow),    emter(1,1,krow)                      &
!
            , aphm1,             apm1,              geom1             &
            , ztvm1                                                   &
!!$#ifndef MESSY
            , zvdiffp,           zvmixtau          &
!!$#endif
!
            , zcvs,              zcvw,              srfl              &
            , qhfla,             evapot                               &
            , ztslnew,           zwlmx                                &
            , zfrl,              zfrw,              zfri              &
            , loland,            loglac                               &
!
!!$#ifndef MESSY
!!$            , xtte(1,1,1,krow)                                        &
            , xtte(1,1,:,krow)                                        &
!!$#endif
!
            , vol(1,1,krow),     vom(1,1,krow),     qte(1,1,krow)     &
            , tte(1,1,krow),     xlte(1,1,krow),    xite(1,1,krow))
#endif
! mz_ho_20160412+
#endif
! mz_ho_20160412-
!
!
  IF (ltdiag) THEN
! store VDIFF increment
#ifndef MESSY
    pdiga(1:nglon,:, 3,krow) = pdiga(1:nglon,:, 3,krow) + vom(1:nproma,:,krow)
    pdiga(1:nglon,:, 8,krow) = pdiga(1:nglon,:, 8,krow) + vol(1:nproma,:,krow)
    pdiga(1:nglon,:,16,krow) = pdiga(1:nglon,:,16,krow) + tte(1:nproma,:,krow)

! prepare next fields for RADHEAT
    pdiga(1:nglon,:,15,krow) = pdiga(1:nglon,:,15,krow) - tte(1:nproma,:,krow)
#else
    ! mz_jb_20040830+
    pdiga(1:nproma,:, 3,krow) = pdiga(1:nproma,:, 3,krow)+ vom(1:nproma,:,krow)
    pdiga(1:nproma,:, 8,krow) = pdiga(1:nproma,:, 8,krow)+ vol(1:nproma,:,krow)
    pdiga(1:nproma,:,16,krow) = pdiga(1:nproma,:,16,krow)+ tte(1:nproma,:,krow)
    ! prepare next fields for RADHEAT
    pdiga(1:nproma,:,15,krow) = pdiga(1:nproma,:,15,krow)- tte(1:nproma,:,krow)
    ! mz_jb_20040830-
#endif
  ENDIF
!
#ifndef MESSY
  IF (ntrac.GT.0) THEN
     CALL xtdriver1 (nproma, nbdim, nlev, nlevp1, ntrac, aphp1(:,:),   &
                     app1(:,:), tte(:,:,krow), tm1(:,:,krow),          &
                     xtm1(:,:,:,krow), xtte(:,:,:,krow))
  ENDIF
#endif
!
!------------------------------------------------------------------------------
!
!*            ADD RADIATION TENDENCIES EVERY TIME STEP.
!
420 CONTINUE
!
#ifndef MESSY
  IF (lrad) THEN
#ifdef ARGCHECK
  CALL radheat (nproma, nbdim, nlev, nlevp1,                           &
          krow,                                                        &
          zi0,                                                         &
          tm1(:,:,krow)  ,   qm1(:,:,krow),                            &
          trsof(:,:,krow),   trsol(:,:,krow),                          &
          emtef(:,:,krow),   emter(:,:,krow),                          &
          emtef0(:,:,krow),  trsof0(:,:,krow),                         &
          ztrfll,            ztrflw,           ztrfli,                 &
          zsofll,            zsoflw,           zsofli,                 &
          trfllac(:,krow),   trflwac(:,krow),  trfliac(:,krow),        &
          sofllac(:,krow),   soflwac(:,krow),  sofliac(:,krow),        &
          srad0(:,krow),     srads(:,krow),                            &
          sradl(:,krow),     srafl(:,krow),                            &
          srad0u(:,krow),    sradsu(:,krow),                           &
          sraf0(:,krow),     srafs(:,krow),                            &
          srad0d(:,krow),                                              &
          trad0(:,krow),     trads(:,krow),                            &
          tradl(:,krow),     trafl(:,krow),                            &
          traf0(:,krow),     trafs(:,krow),                            &
          tradsu(:,krow),                                              &
          tslm1(:,krow),     tsi(:,krow),      tsw(:,krow),            &
          albedo(:,krow),                                              &
          alsol(:,krow),     alsow(:,krow),    alsoi(:,krow),          &
          aphm1,             apm1,                                     &
          ztslnew,           tte(:,:,krow),                            &
          zfrl,              zfrw,             zfri)   
#else
 CALL radheat (nproma, nbdim, nlev, nlevp1,                           &
          krow,                                                        &
          zi0,                                                         &
          tm1(1,1,krow)  ,   qm1(1,1,krow),                            &
          trsof(1,1,krow),   trsol(1,1,krow),                          &
          emtef(1,1,krow),   emter(1,1,krow),                          &
          emtef0(1,1,krow),  trsof0(1,1,krow),                         &
          ztrfll,            ztrflw,           ztrfli,                 &
          zsofll,            zsoflw,           zsofli,                 &
          trfllac(1,krow),   trflwac(1,krow),  trfliac(1,krow),        &
          sofllac(1,krow),   soflwac(1,krow),  sofliac(1,krow),        &
          srad0(1,krow),     srads(1,krow),                            &
          sradl(1,krow),     srafl(1,krow),                            &
          srad0u(1,krow),    sradsu(1,krow),                           &
          sraf0(1,krow),     srafs(1,krow),                            &
          srad0d(1,krow),                                              &
          trad0(1,krow),     trads(1,krow),                            &
          tradl(1,krow),     trafl(1,krow),                            &
          traf0(1,krow),     trafs(1,krow),                            &
          tradsu(1,krow),                                              &
          tslm1(1,krow),     tsi(1,krow),      tsw(1,krow),            &
          albedo(1,krow),                                              &
          alsol(1,krow),     alsow(1,krow),    alsoi(1,krow),          &
          aphm1,             apm1,                                     &
          ztslnew,           tte(1,1,krow),                            &
          zfrl,              zfrw,             zfri)
#endif
  ELSE
     ! lrad=.FALSE.
     ! --> no radiative effect at the surface
     ztrfll(:)=0._dp
     ztrflw(:)=0._dp
     ztrfli(:)=0._dp
     zsofll(:)=0._dp
     zsoflw(:)=0._dp
     zsofli(:)=0._dp
     ! --> lake, ml_ocean, licetemp and sicetemp get zero fluxes 
  END IF
#endif
! ... of #ifndef MESSY ! op_pj_20130407
!
#ifdef MESSY
! mz_bs_20050228+
   CALL messy_radheat
! mz_bs_20050228-
#endif
  IF (ltdiag) THEN
! store RADHEAT increment
#ifndef MESSY
    pdiga(1:nglon,:,15,krow) = pdiga(1:nglon,:,15,krow) + tte(1:nproma,:,krow)

! prepare next fields for GWDRAG
    pdiga(1:nglon,:, 4,krow) = pdiga(1:nglon,:, 4,krow) - vom(1:nproma,:,krow)
    pdiga(1:nglon,:, 9,krow) = pdiga(1:nglon,:, 9,krow) - vol(1:nproma,:,krow)
    pdiga(1:nglon,:,17,krow) = pdiga(1:nglon,:,17,krow) - tte(1:nproma,:,krow)
#else
    ! mz_jb_20040830+
    pdiga(1:nproma,:,15,krow) = pdiga(1:nproma,:,15,krow)+ tte(1:nproma,:,krow)
    ! prepare next fields for GWDRAG
    pdiga(1:nproma,:, 4,krow) = pdiga(1:nproma,:, 4,krow)- vom(1:nproma,:,krow)
    pdiga(1:nproma,:, 9,krow) = pdiga(1:nproma,:, 9,krow)- vol(1:nproma,:,krow)
    pdiga(1:nproma,:,17,krow) = pdiga(1:nproma,:,17,krow)- tte(1:nproma,:,krow)
    ! mz_jb_20040830-
#endif
  ENDIF
!
!
!     ------------------------------------------------------------
!
!        ***  GRAVITY WAVE DRAG PARAMETERISATION  ***
!
  IF (lmidatm) THEN
    IF (lresume) THEN
       aprflux(1:nproma,krow) = 0._dp
       aprfluxm(1:nproma,krow) = 0._dp
    ENDIF

#ifndef MESSY
    CALL gwspectrum ( krow, nproma,  nbdim,  nlev,                             &
                  aphm1,            apm1,                                      &
                  tm1(:,:,krow),    um1(:,:,krow),   vm1(:,:,krow),            &
                  aprflux(:,krow),                                             &
                  tte(:,:,krow),    vol(:,:,krow),   vom(:,:,krow) )
#endif
  END IF

#ifdef MESSY
  CALL messy_gwdrag
#endif

! op_re_20160608+
#ifndef MESSY
  IF (lgwdrag) THEN

    CALL ssodrag ( nproma,  nbdim,  nlev,                                      &
                  aphm1,            apm1,            geom1,                    &
                  tm1(:,:,krow),    um1(:,:,krow),   vm1(:,:,krow),            &
                  oromea(:,krow),   orostd(:,krow),  orosig(:,krow),           &
                  orogam(:,krow),                                              &
                  orothe(:,krow),   oropic(:,krow),  oroval(:,krow),           &
                  ustrgw(:,krow),   vstrgw(:,krow),  vdisgw(:,krow)            &
!!$#ifndef MESSYTENDENCY
                  , tte(:,:,krow),    vol(:,:,krow),   vom(:,:,krow)           &
!!$#else
                  )
  ! NOTE: The tendencies must not be passed to the subroutine, since with
  ! MESSYTENDENCY they are not longer modified directly (but rather indirectly
  ! via mtend_add_l internally).  Passing the tendency in addition causes the
  ! compiler to copy back the result of the calculation operating on the local
  ! INTENT(inout) variable (e.g. ptte), which then overwrites the indirectly
  ! modified tendencies ... weird!  
  ! CALL ssodrag ( ..., tte, ) COMPILER: copy tte into ptte 
  ! ssodrag: no direct modification of ptte, but mtend_add_l to tte 
  ! COMPILER: copy ppte back to tte (this overwrites the modified tte of
  !                                  the step before)
!!$#endif

  END IF
#endif
! op_re_20160608-

  IF (ltdiag) THEN
! store GWDRAG increment
#ifndef MESSY
    pdiga(1:nglon,:, 4,krow) = pdiga(1:nglon,:, 4,krow) + vom(1:nproma,:,krow)
    pdiga(1:nglon,:, 9,krow) = pdiga(1:nglon,:, 9,krow) + vol(1:nproma,:,krow)
    pdiga(1:nglon,:,17,krow) = pdiga(1:nglon,:,17,krow) + tte(1:nproma,:,krow)

! prepare next fields for CUCALL
    pdiga(1:nglon,:, 5,krow) = pdiga(1:nglon,:, 5,krow) - vom(1:nproma,:,krow)
    pdiga(1:nglon,:,10,krow) = pdiga(1:nglon,:,10,krow) - vol(1:nproma,:,krow)
    pdiga(1:nglon,:,18,krow) = pdiga(1:nglon,:,18,krow) - tte(1:nproma,:,krow)
#else
    ! mz_jb_20040830+
    pdiga(1:nproma,:, 4,krow) = pdiga(1:nproma,:, 4,krow)+ vom(1:nproma,:,krow)
    pdiga(1:nproma,:, 9,krow) = pdiga(1:nproma,:, 9,krow)+ vol(1:nproma,:,krow)
    pdiga(1:nproma,:,17,krow) = pdiga(1:nproma,:,17,krow)+ tte(1:nproma,:,krow)
    ! prepare next fields for CUCALL
    ! this is done in MESSY_CONVECT_E5 ! mz_ht_20050330
    ! mz_jb_20040830-
#endif
  ENDIF
!
!     ------------------------------------------------------------------
!
!*        6.    CONVECTION PARAMETERISATION.
!               ---------- -----------------
!
600 CONTINUE
!
!
!*        6.3    COMPUTE *T* AND *Q* TENDENCIES BY MOIST CONVECTION.
!*                AND SHALLOW CONVECTION.
!
630 CONTINUE
!
!
#ifndef MESSY
   itype(1:nproma)=0  !!!!!
#endif
!
!
!
!*         6.3.1   INITIALIZE ARRAYS FOR CONVECTIVE PRECIPITATION
!*                 AND COPY ARRAYS FOR CONVECTIVE CLOUD PARAMETERS
!*                 -----------------------------------------------
!
631 CONTINUE
!
!
   xtec(1:nproma,:,krow)=0._dp  !!!!! 
!
   IF (.NOT. USE_CRM) THEN
     DO 632 jl=1,nproma
       rsfc(jl)=0._dp
       ssfc(jl)=0._dp
632  END DO
   END IF
!
!
!*         6.3.3   CALL SUBROUTINE CUCALL FOR CUMULUS PARAMETERIZATION
!                  ---------------------------------------------------
!
633 CONTINUE
!
!
#ifndef MESSY
  IF (loconv) THEN
!
#ifdef ARGCHECK
     CALL cucall(nproma, nbdim, nlev, nlevp1, nlevm1, ilab,            &
                ntrac,                                                 &
                xtm1(:,:,:,krow), xtte(:,:,:,krow),                    &
                tm1(:,:,krow),    qm1(:,:,krow),    um1(:,:,krow),     &
                vm1(:,:,krow),    xlm1(:,:,krow),   xim1(:,:,krow),    &
                tte(:,:,krow),    qte(:,:,krow),    vom(:,:,krow),     &
                vol(:,:,krow),    xlte(:,:,krow),   xite(:,:,krow),    &
                vervel(:,:,krow), xtec(:,:,krow),                      &
                zqtec,            qhfla,                               &
                app1,             aphp1,            geom1,             &
                rsfc,             ssfc,                                &
                aprc(:,krow),     aprs(:,krow),                        &
                itype,            loland,                              &
                topmax(:,krow)  )
#else
     CALL cucall(nproma, nbdim, nlev, nlevp1, nlevm1, ilab,            &
                ntrac,                                                 &
!!$                xtm1(1,1,1,krow), xtte(1,1,1,krow),                    &
                xtm1(1,1,:,krow), xtte(1,1,:,krow),                    &
                tm1(1,1,krow),    qm1(1,1,krow),    um1(1,1,krow),     &
                vm1(1,1,krow),    xlm1(1,1,krow),   xim1(1,1,krow),    &
                tte(1,1,krow),    qte(1,1,krow),    vom(1,1,krow),     &
                vol(1,1,krow),    xlte(1,1,krow),   xite(1,1,krow),    &
                vervel(1,1,krow), xtec(1,1,krow),                      &
                zqtec,            qhfla,                               &
                app1,             aphp1,            geom1,             &
                rsfc,             ssfc,                                &
                aprc(1,krow),     aprs(1,krow),                        &
                itype,            loland,                              &
                topmax(1,krow)  )
#endif
!
     DO jl=kidia,kfdia       
       rtype(jl,krow)=REAL(itype(jl),dp)
     END DO
!
  ELSE
!       NECESSARY COMPUTATIONS IF MASSFLUX IS BY-PASSED
!
     ilab(1:nproma,1:nlev)=0
!
  ENDIF
#else
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if (loconv.or.locond) then
    CALL messy_convec
  ELSE
!       NECESSARY COMPUTATIONS IF MASSFLUX IS BY-PASSED

    if (.not.loconv) zilab(1:nproma,1:nlev)=0.0_dp ! mz_ht_20040405
  ENDIF
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif
!
#ifndef MESSY
  IF (ltdiag) THEN
! store CUCALL increment (massflux)
    pdiga(1:nglon,:, 5,krow) = pdiga(1:nglon,:, 5,krow) + vom(1:nproma,:,krow)
    pdiga(1:nglon,:,10,krow) = pdiga(1:nglon,:,10,krow) + vol(1:nproma,:,krow)
    pdiga(1:nglon,:,18,krow) = pdiga(1:nglon,:,18,krow) + tte(1:nproma,:,krow)

! prepare next fields for COND
    pdiga(1:nglon,:,19,krow) = pdiga(1:nglon,:,19,krow) - tte(1:nproma,:,krow)
  ENDIF
#endif
!
!     ------------------------------------------------------------
!
!*       7.    LARGE SCALE CONDENSATION.
!              ----- ----- -------------
!
700 CONTINUE
!
!
  IF(locond) THEN
!
!     POLE FILTER FOR TENDENCIES
!
     IF (iadvec == spitfire .and. .not. ldc%col_1d)                    & 
       CALL pole_filter(tte(:,:,krow), qte(:,:,krow), krow)
#ifdef OBSOLETE
     IF (lcotra) CALL get_col_pol (tte(:,:,krow), qte(:,:,krow), krow)
#endif

#ifdef MESSY
! mz_ht_20041102+
! separate routine from MESSY cloud called at the end of MESSY convec
#else
     IF (ltimer) CALL timer_start(timer_cloud)
#ifdef ARGCHECK
     CALL cloud( nproma, nbdim, ktdia, nlev, nlevp1                   &
               , aphm1,             vervel(:,:,krow)                  &
               , apm1,              app1,             acdnc(:,:,krow) &
               , qm1(:,:,krow),     tm1(:,:,krow),    ztvm1           &
               , xlm1(:,:,krow),    xim1(:,:,krow),   xtec(:,:,krow)  &
               , xvar(:,:,krow),    xskew(:,:,krow),  zqtec           &
               , zbetaa,            zbetab                            &
               , zvdiffp,           zhmixtau,         zvmixtau        &
               , geom1,             zbetass                           &
               , invb                                                 &
               , aclc(:,:,krow),    aclcac(:,:,krow)                  &
               , relhum(:,:,krow)                                     &
               , aclcov(:,krow),    aprl(:,krow),     qvi(:,krow)     &
               , xlvi(:,krow),      xivi(:,krow)                      &
               , ssfl,              rsfl                              &
               , qte(:,:,krow),     tte(:,:,krow)                     &
               , xlte(:,:,krow),    xite(:,:,krow)                    &
               , aprs(:,krow)     )
#else
     CALL cloud( nproma, nbdim, ktdia, nlev, nlevp1                   &
               , aphm1,             vervel(1,1,krow)                  &
               , apm1,              app1,             acdnc(1,1,krow) &
               , qm1(1,1,krow),     tm1(1,1,krow),    ztvm1           &
               , xlm1(1,1,krow),    xim1(1,1,krow),   xtec(1,1,krow)  &
               , xvar(1,1,krow),    xskew(1,1,krow),  zqtec           &
               , zbetaa,            zbetab                            &
               , zvdiffp,           zhmixtau,         zvmixtau        &
               , geom1,             zbetass                           &
               , invb                                                 &
               , aclc(1,1,krow),    aclcac(1,1,krow)                  &
               , relhum(1,1,krow)                                     &
               , aclcov(1,krow),    aprl(1,krow),     qvi(1,krow)     &
               , xlvi(1,krow),      xivi(1,krow)                      &
               , ssfl,              rsfl                              &
               , qte(1,1,krow),     tte(1,1,krow)                     &
               , xlte(1,1,krow),    xite(1,1,krow)                    &
               , aprs(1,krow)     )
#endif
     IF (ltimer) CALL timer_stop(timer_cloud)
#endif
  ELSE
!
!              NECESSARY COMPUTATIONS IF *CLOUD* IS BY-PASSED.
!
! um_hr_20190301+
!!#D crm +
#ifdef MESSY
     IF (.NOT. USE_CRM) THEN
#endif
!!#D crm -
! um_hr_20190301-
    ssfl(1:nproma) = 0._dp
    rsfl(1:nproma) = 0._dp
!
    aclc(1:nproma,:,krow) = 0._dp
! um_hr_20190301+
!!#D crm +
#ifdef MESSY
    END IF
#endif
!!#D crm -
! um_hr_20190301-

!
  ENDIF
!
! store COND increment
#ifndef MESSY
  IF (ltdiag) &
       pdiga(1:nglon,:, 19,krow) = pdiga(1:nglon,:, 19,krow) + tte(1:nproma,:,krow)
#else
  ! for MESSy this is done in messy_cloud_e5 ! mz_ht_20050330
#endif
 
  IF (lmidatm) THEN
    DO jl=1,nproma
       aprflux(jl,krow)=rsfl(jl)+ssfl(jl)+rsfc(jl)+ssfc(jl)
    END DO
  END IF
  DO jl=1,nproma
    zprecip(jl)=rsfl(jl)+ssfl(jl)+rsfc(jl)+ssfc(jl)
  END DO
!
!     ----------------------------------------------------------------
!
!*        8.    Computation of new surface values over land points
!
!
#ifndef MESSY
! op_sb_20121205
  IF (lsurf) THEN   

#ifdef ARGCHECK
    CALL surf ( nproma,          nbdim,             nlev, nlevp1      &
            , tsl(:,krow),       tslm(:,krow),      tslm1(:,krow)     &
            , ws(:,krow),        wl(:,krow),        wsmx(:,krow)      &
            , sn(:,krow),        snmel(:,krow),     gld(:,krow)       &
            , snc(:,krow),       u10(:,krow),       v10(:,krow)       &
            , runoff(:,krow),    rogl(:,krow),      drain(:,krow)     &
            , apmegl(:,krow),    snacl(:,krow),     orostd(:,krow)    &
            , rgcgn(:,krow),     sodif(:,krow),     slm(:,krow)       &
            , grndcapc(:,krow),  grndhflx(:,krow),  grndflux(:,krow)  &
! 
            , tsoil(:,:,krow),   grndd(:,:,krow),   grndc(:,:,krow)   &
            , tm1(:,:,krow),     qm1(:,:,krow),     tte(:,:,krow)     &
            , aphm1                                                   &
! 
            , zcvs,              zcvw,              zwlmx             &
            , zevapl,            evapot                               &
            , rsfl,              rsfc                                 &
            , ssfl,              ssfc                                 &
            , zros_hd,           zdrain_hd                            &
            , zalac                                                   &
            , loland,            loglac     )
#else
    CALL surf ( nproma,          nbdim,             nlev, nlevp1      &
            , tsl(1,krow),       tslm(1,krow),      tslm1(1,krow)     &
            , ws(1,krow),        wl(1,krow),        wsmx(1,krow)      &
            , sn(1,krow),        snmel(1,krow),     gld(1,krow)       &
            , snc(1,krow),       u10(1,krow),       v10(1,krow)       &
            , runoff(1,krow),    rogl(1,krow),      drain(1,krow)     &
            , apmegl(1,krow),    snacl(1,krow),     orostd(1,krow)    &
            , rgcgn(1,krow),     sodif(1,krow),     slm(1,krow)       &
            , grndcapc(1,krow),  grndhflx(1,krow),  grndflux(1,krow)  &
! 
            , tsoil(1,1,krow),   grndd(1,1,krow),   grndc(1,1,krow)   &
            , tm1(1,1,krow),     qm1(1,1,krow),     tte(1,1,krow)     &
            , aphm1                                                   &
! 
            , zcvs,              zcvw,              zwlmx             &
            , zevapl,            evapot                               &
            , rsfl,              rsfc                                 &
            , ssfl,              ssfc                                 &
            , zros_hd,           zdrain_hd                            &
            , zalac                                                   &
            , loland,            loglac     )
#endif
!
!
!*     9.  Lake physics
!
!
#ifdef ARGCHECK
    CALL lake ( nproma                                                 &
            , seaice(:,krow),    siced(:,krow),    alake(:,krow)      &
            , tsi(:,krow),       tsw(:,krow)                          &
            , zhflw,             zhfsw,            fluxres(:,krow)    &
            , ztrflw,            zsoflw                               &
            , zevapi,            sni(:,krow),      zcvsi              &
            , ahfres(:,krow),    zfri           )
#else
    CALL lake ( nproma                                                 &
            , seaice(1,krow),    siced(1,krow),    alake(1,krow)      &
            , tsi(1,krow),       tsw(1,krow)                          &
            , zhflw,             zhfsw,            fluxres(1,krow)    &
            , ztrflw,            zsoflw                               &
            , zevapi,            sni(1,krow),      zcvsi              &
            , ahfres(1,krow),    zfri           )
#endif

#ifdef ARGCHECK
    CALL licetemp ( nproma                                             &
               , siced(:,krow),     sni(:,krow),      alake(:,krow)   &
               , tsi(:,krow),       ztrfli,           zsofli          &
               , ahfice(:,krow),    fluxres(:,krow)                   &
               , ahfcon(:,krow),    ahfres(:,krow),   zevapi          &
               , ssfl,              ssfc                              &
               , zhfsi,             zhfli,            zcvsi           &
               , zfri                            )
#else
    CALL licetemp ( nproma                                             &
               , siced(1,krow),     sni(1,krow),      alake(1,krow)   &
               , tsi(1,krow),       ztrfli,           zsofli          &
               , ahfice(1,krow),    fluxres(1,krow)                   &
               , ahfcon(1,krow),    ahfres(1,krow),   zevapi          &
               , ssfl,              ssfc                              &
               , zhfsi,             zhfli,            zcvsi           &
               , zfri                            )
#endif 
!
!
!*   10.1  Compute mixed layer ocean physics
!
!
    IF (lmlo) THEN
#ifdef ARGCHECK
      CALL ml_ocean ( nproma                                          &
            , slm(:,krow)                                             &
            , lonorth                                                 &
            , seaice(:,krow),    siced(:,krow),    alake(:,krow)      &
            , tsi(:,krow),       tsw(:,krow)                          &
            , zhflw,             zhfsw,            fluxres(:,krow)    &
            , ztrflw,            zsoflw                               &
            , amlcorr(:,krow),   amlcorac(:,krow), amlheatac(:,krow)  &
            , zevapi,            sni(:,krow),      zcvsi              &
            , ahfres(:,krow),    zfri           )
#else
      CALL ml_ocean ( nproma                                          &
            , slm(1,krow)                                             &
            , lonorth                                                 &
            , seaice(1,krow),    siced(1,krow),    alake(1,krow)      &
            , tsi(1,krow),       tsw(1,krow)                          &
            , zhflw,             zhfsw,            fluxres(1,krow)    &
            , ztrflw,            zsoflw                               &
            , amlcorr(1,krow),   amlcorac(1,krow), amlheatac(1,krow)  &
            , zevapi,            sni(1,krow),      zcvsi              &
            , ahfres(1,krow),    zfri           )
#endif
    END IF 

!*   10.2  Compute sea ice skin temperature

#ifdef ARGCHECK
    CALL sicetemp ( nproma                                            &
               , siced(:,krow),     sni(:,krow),      alake(:,krow)   &
               , slf(:,krow)                                          &
               , tsi(:,krow),       ztrfli,           zsofli          &
               , ahfice(:,krow),    fluxres(:,krow),  qres(:,krow)    &
               , ahfcon(:,krow),    ahfres(:,krow)                    &
               , zhfsi,             zhfli                             &
               , zfri                            )
#else
    CALL sicetemp ( nproma                                            &
               , siced(1,krow),     sni(1,krow),      alake(1,krow)   &
               , slf(1,krow)                                          &
               , tsi(1,krow),       ztrfli,           zsofli          &
               , ahfice(1,krow),    fluxres(1,krow),  qres(1,krow)    &
               , ahfcon(1,krow),    ahfres(1,krow)                    &
               , zhfsi,             zhfli                             &
               , zfri                            )
#endif
  END IF
#else
! op_sb_20121008+
! moved from above
    ! fb_mk_20100407+
    ! MESSy interface to the mixed layer ocean and surf ...
    CALL messy_mixlo
    ! fb_mk_20100407-
! op_sb_20121008-
#endif
! end of #ifndef MESSY
!
!
!    compute surface temperature for diagnostics
!
#ifndef MESSY
!DIR$ CONCURRENT
    DO jl=1,nproma
      ztsurf= zfrl(jl)*tslm1(jl,krow)              &
             +zfri(jl)*tsi(jl,krow)                &
             +zfrw(jl)*tsw(jl,krow)
      tsurf(jl,krow)=tsurf(jl,krow)+zdtime*ztsurf
    END DO
#else
! moved to SURFACE
#endif

    IF (lnmi_run) THEN
      SELECT CASE(nmi_phase)
      CASE(NMI_ACCU)    ! accumulate tendencies
!DIR$ CONCURRENT
        dh_t(:,:,krow) = dh_t(:,:,krow) + tte(:,:,krow) - buf_t(:,:)
!DIR$ CONCURRENT
        dh_m(:,:,krow) = dh_m(:,:,krow) + vom(:,:,krow) - buf_m(:,:)
!DIR$ CONCURRENT
        dh_l(:,:,krow) = dh_l(:,:,krow) + vol(:,:,krow) - buf_l(:,:)

      CASE(NMI_USE_AVG) ! prepare/reset tendencies
!DIR$ CONCURRENT
        tte(:,:,krow) = buf_t(:,:) + dh_t(:,:,krow)
!DIR$ CONCURRENT
        vom(:,:,krow) = buf_m(:,:) + dh_m(:,:,krow)
!DIR$ CONCURRENT
        vol(:,:,krow) = buf_l(:,:) + dh_l(:,:,krow)

      END SELECT
    END IF
!
!      p-e budget correction for coupling
!
      DO jl=1,nproma
        zqtold(jl)=qtnew(jl,krow)
        qtnew(jl,krow)=0._dp
      END DO
      DO jk=ktdia,nlev
!DIR$ CONCURRENT
        DO jl=1,nproma
          qtnew(jl,krow)=qtnew(jl,krow)+(qm1(jl,jk,krow)              &
                         +xlm1(jl,jk,krow)+xim1(jl,jk,krow))          &
                         *(aphm1(jl,jk+1)-aphm1(jl,jk))/g
        END DO
      END DO
!     
!     Accumulate p-e correction for standard diagnostics
!
!DIR$ CONCURRENT
      DO jl=1,nproma
        apmeb(jl,krow)=apmeb(jl,krow)-(qtnew(jl,krow)-zqtold(jl))     &
                       -(zprecip(jl)+qhfla(jl))*zdtime
      END DO
!
#ifndef MESSY
!     Vertical integral of anthropogenic sulfur burden 
!
      IF(lso4) THEN
        DO jk=ktdia,nlev
          DO jl=1,nproma
            abso4(jl,krow)=abso4(jl,krow)                               &
                    +(so4all(jl,jk,krow)-so4nat(jl,jk,krow))            &
                    *(aphm1(jl,jk+1)-aphm1(jl,jk))/g*zdtime
          END DO
        END DO
      ENDIF
!
#endif
!
#ifndef MESSY
!mz_ap_20070830+
 IF (lcouple) THEN
!
      IF (l_getocean) apmebco(:,krow) = 0._dp ! set to zero after coupling

!     Accumulate p-e correction for coupling
!DIR$ CONCURRENT
      DO jl=1,nproma
        apmebco(jl,krow)=apmebco(jl,krow)-(qtnew(jl,krow)-zqtold(jl)) &
                       -(zprecip(jl)+qhfla(jl))*zdtime
      END DO
!
!
!       collect data needed as input for the ocean model
!
#ifdef ARGCHECK
   CALL collect ( nproma                                               &
        , zhflw,          zhfsw,         ahfice(:,krow)                &
        , ztrflw,         zsoflw                                       &
        , qres(:,krow),   zevapw,        zevapi                        &
        , ustrw(:,krow),  vstrw(:,krow), ustri(:,krow), vstri(:,krow)  &
        , alake(:,krow),  slf(:,krow),   seaice(:,krow)                &
        , wind10w(:,krow)                                              &
        , awhea(:,krow),  awsol(:,krow), awfre(:,krow), awust(:,krow)  &
        , awvst(:,krow),  aicon(:,krow), aiqre(:,krow), aifre(:,krow)  &
        , aiust(:,krow),  aivst(:,krow), awsta(:,krow)                 &
        , rsfc,           ssfc,          rsfl,          ssfl         )
#else
   CALL collect ( nproma                                               &
        , zhflw,          zhfsw,         ahfice(1,krow)                &
        , ztrflw,         zsoflw                                       &
        , qres(1,krow),   zevapw,        zevapi                        &
        , ustrw(1,krow),  vstrw(1,krow), ustri(1,krow), vstri(1,krow)  &
        , alake(1,krow),  slf(1,krow),   seaice(1,krow)                &
        , wind10w(1,krow)                                              &
        , awhea(1,krow),  awsol(1,krow), awfre(1,krow), awust(1,krow)  &
        , awvst(1,krow),  aicon(1,krow), aiqre(1,krow), aifre(1,krow)  &
        , aiust(1,krow),  aivst(1,krow), awsta(1,krow)                 &
        , rsfc,           ssfc,          rsfl,          ssfl         )
#endif
 END IF
!
!     collect data needed as input for the HD-model 
!     and glacier calving model
!    (includes diagnostics of discharge and calving on atmospheric grid)
!
 IF (lhd) THEN
     CALL hydrology_collect ( nproma                                   &
                            , aros(:,krow),    adrain(:,krow)          &
                            , apmecal(:,krow)                          &
                            , disch(:,krow),   runtoc(:,krow)          &
                            , zros_hd,         zdrain_hd               &
                            , zalac )
 END IF
!mz_ap_20070920: end of #ifdef MESSY
#endif

#ifdef MESSY
! mz_ap_20161103+
!!$!mz_ap_20070820+
!!$!!!$ IF (lcouple) THEN ! mz_pj_20071129
!!$    ! save some istantaneous data for later use in the ocean exchange
!!$    ! (a2o_global_end, via a2o.nml)
!!$    ! latent + sensible heat+ SW + LW radiation 
!!$    awhea_2d(1:nproma,krow) = ( zhflw(1:nproma) + zhfsw(1:nproma)  &
!!$         + radflxw_2d(1:nproma,krow) )
!!$    aicon_2d(1:nproma,krow) = ahfice(1:nproma,krow)
!!$    aiqre_2d(1:nproma,krow) = qres(1:nproma,krow) 
!!$    awfre_2d(1:nproma,krow) = (  (rsfl(1:nproma)+rsfc(1:nproma)) + &
!!$         (ssfl(1:nproma)+ssfc(1:nproma)+zevapw(1:nproma))*         &
!!$         (1._dp-seaice(1:nproma,krow))) / &
!!$         rhoh2o
!!$    aifre_2d(1:nproma,krow) = ((ssfl(1:nproma)+ssfc(1:nproma) + &
!!$         zevapi(1:nproma)) * (seaice(1:nproma,krow))) / rhoh2o
!!$    aiust_2d(1:nproma,krow) = ustri(1:nproma,krow)
!!$    aivst_2d(1:nproma,krow) = vstri(1:nproma,krow)
!!$    awust_2d(1:nproma,krow) = ustrw(1:nproma,krow)
!!$    awvst_2d(1:nproma,krow) = vstrw(1:nproma,krow)
!!$    
!!$    ! save some istantaneous data for later use in the
!!$    ! hydrological discharge model (hd_global_end)
!!$    ! Note: aros_2d, adrain_2d, apmecal_2d are channel objects
!!$    !       if ECHAM5 (i.e. in messy_main_data_e5.f90) and
!!$    !       calculated in submodel SURFACE
!!$    apmebco_2d(1:nproma,krow)   =                    &
!!$         -(qtnew(1:nproma,krow)-zqtold(1:nproma))    &
!!$         -(zprecip(1:nproma)+qhfla(1:nproma))*zdtime
!!$!!!$ END IF ! mz_pj_20071129
!!$!mz_ap_20070820-
! mz_ap_20161103-

! mz_rs_20040326+
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CALL messy_physc
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!mz_rs_20040326-
#endif
!
! allow submodels to calculate some diagnostics
!

#ifndef MESSY
  IF (ntrac>0) THEN
    CALL xtdriver2 (nproma, nbdim, nlev, nlevp1, ntrac, aphp1(:,:),    &
     app1(:,:), tte(:,:,krow), tm1(:,:,krow), xtm1(:,:,:,krow), xtte(:,:,:,krow))
    CALL xtdiagn
  END IF
#endif
!
! daily block statistics - AMIP2 global diagnostics
!  
#ifdef OBSOLETE
  IF(ldiagamip) THEN
    CALL collect_amip2_diag(nproma,nbdim,nlev,nlevp1,krow                    &
       ,tm1(:,:,krow),um1(:,:,krow),vm1(:,:,krow),aphm1(:,:),geospm(:,krow)  &
       ,ustr(:,krow),ustrgw(:,krow),ustrm(:,krow),ustrgwm(:,krow)       &
       ,tslm1(:,krow),seaicem(:,krow),loland)
  END IF
#endif
#ifdef MESSY
  ! mz_lg_20030116+ save the convective and large-scale rainfall
  prc(1:nproma,krow)=rsfc(1:nproma)*zdtime/rhoh2o
  prl(1:nproma,krow)=rsfl(1:nproma)*zdtime/rhoh2o
  ! mz_lg_20030116-
#endif
!
!      
!*        9.    RESTORE WINDS AND WIND TENDENCIES.
!
900 CONTINUE
!
!
!
!*        9.2   RESTORE WINDS AND WIND TENDENCIES.
!
920 CONTINUE
!
  DO jlev = 1, nlev
!DIR$ CONCURRENT
    DO jl = 1, nproma
      zcst = sqcst_2d(jl,krow)
      um1(jl,jlev,krow) = zcst*um1(jl,jlev,krow)
      vm1(jl,jlev,krow) = zcst*vm1(jl,jlev,krow)
      vol(jl,jlev,krow) = zcst*vol(jl,jlev,krow)
      vom(jl,jlev,krow) = zcst*vom(jl,jlev,krow)
    END DO
  END DO
#ifdef MESSYTENDENCY
  ! from here horiz. wind vector scaled again with cos latitude
  call mtend_set_sqcst_scal(.true.)
#endif
!
!
!     ------------------------------------------------------------
!
!*       10.    RELEASE SPACE.
!               ------- ------
!
1200 CONTINUE

#ifdef MESSY
  DEALLOCATE (aphm1)
  DEALLOCATE (apm1)
  DEALLOCATE (aphp1)
  DEALLOCATE (app1)
  NULLIFY(geom1)
  NULLIFY(ssfc)
  NULLIFY(ssfl)
  NULLIFY(rsfc)
  NULLIFY(rsfl)
  NULLIFY(srfl)
#endif
!
!     ------------------------------------------------------------
!
  RETURN
END SUBROUTINE physc
