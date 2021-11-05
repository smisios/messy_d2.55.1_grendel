SUBROUTINE vdiff ( kproma, kbdim, ktdia, klev, klevm1, klevp1, ktrac   &
         , krow                                                        &
!-----------------------------------------------------------------------
! - 3D from mo_memory_g1a
         , pxtm1                                                       &
! - 2D from mo_memory_g1a
!!$#ifndef MESSY
         , pqm1,       ptm1,       pum1                                &
         , pvm1,       pxlm1,      pxim1                               &
!!$#else
!!$! mz_ht_20041209+        
!!$
!!$         , pqmm1,      ptmm1,      pumm1                               &
!!$         , pvmm1,      pxlmm1,     pximm1                              &
!!$! mz_ht_20041209-
!!$#endif
         , pxvar                                                       &
! - 1D
         , pahfl,      pahfs,      paz0                                &
         , pdew2,      pevap,      pforest                             &
         , ptemp2,     pt2max                                          &
         , pt2min,     pwind10w,   pvdis                               &
         , pu10,       pv10,       pustr                               &
         , pvstr,      pwimax,     pwind10                             &
         , pwsmx,      pvlt,       pgrndcapc                           &
         , pgrndhflx,  pvgrat                                          &
         , ptsl,       ptsw,       ptsi                                &
         , pocu,       pocv                                            &
         , paz0l,      paz0w,      paz0i                               &
         , pahfsl,     pahfsw,     pahfsi                              &
         , pahfll,     pahflw,     pahfli                              &
         , pevapl,     pevapw,     pevapi                              &
         , pahfslac,   pahfswac,   pahfsiac                            &
         , pahfllac,   pahflwac,   pahfliac                            &
         , pevaplac,   pevapwac,   pevapiac                            &
         , pustrl,     pustrw,     pustri                              &
         , pvstrl,     pvstrw,     pvstri                              &
         , psn,        psnc,       ptslm1                              &
         , pws,        palbedo,    palsol                              &
! - 2D from mo_memory_g3b
         , ptke,       ptkem1,     ptkem                               &
         , paclc,      pemter                                          &
! - 2D within physics only
         , paphm1,     papm1,      pgeom1                              &
         , ptvm1                                                       &
!!$#ifndef MESSY    
         , pvdiffp,    pvmixtau                                        &
!!$#endif 
! - 1D within physics only.
         , pcvs,       pcvw,       psrfl                               &
         , pqhfla,     pevapot                                         &
         , ptslnew,    pwlmx                                           &
         , pfrl,       pfrw,       pfri                                &
         , lpland,     lpglac                                          &
! - Tendencies
! - 3D
!!$#ifndef MESSY
         , pxtte                                                       &
!!$#endif
! - 2D
         , pvol,       pvom,       pqte                                &
         , ptte,       pxlte,      pxite)
!
!**** *vdiff* - does the vertical exchange of u, v, t, q, xl, xi and xt
!               by turbulence.
!
!
!     Subject.
!
!       This routine computes the physical tendencies of the seven
!   prognostic variables u, v, t, q, xl, xi and xt due to the vertical
!   exchange by turbulent (= non-moist convective) processes.
!   These tendencies are obtained as the difference between
!   the result of an implicit time-step starting from values at t-1
!   and these t-1 (???) values.
!   all the diagnostic computations (exchange coefficients, ...) are
!   done from the t-1 values. As a by-product the roughness length
!   over sea is updated accordingly to the *charnock formula. heat and
!   moisture surface fluxes and their derivatives against ts and ws,
!   later to be used for soil processes treatment, are also
!   computed as well as a stability value to be used as a diagnostic
!   of the depth of the well mixed layer in convective computations.
!
!**   Interface.
!     ----------
!
!          *vdiff* is called from *physc*.
!
!      Arguments.
!      ----------
!
!  - 3d from mo_memory_g1a
!
!  pxtm1    : tracer variables (t-dt)
!
!  - 2d from mo_memory_g1a
!
!  pqm1     : humidity (t-dt)
!  ptm1     : temperature (t-dt)
!  pum1     : zonal wind (t-dt)
!  pvm1     : meridional wind (t-dt)
!  pxlm1    : cloud water (t-dt)
!  pxim1    : cloud ice (t-dt)
!  pxvar    : distribution width (b-a) (t-dt)

! - 2d from mo_memory_g3
!
!  paclc    : cloud cover
!
!  ptke     : turbulent kinetic energy at t+dt (unfiltered)
!  ptkem    :            "             at t    (unfiltered)
!  ptkem1   :            "             at t-dt   (filtered)
!
!  - 1d from mo_memory_g3
!
!  ptsl     : surface temperature over land
!  ptsw     :             "       over water
!  ptsi     :             "       over ice
!
!  pocu     : ocean u-velocity
!  pocv     : ocean v-velocity
!
!  pahfs    : surface sensible heat flux (accumulated)
!  pahfsl   :             "              over land
!  pahfsw   :             "              over water
!  pahfsi   :             "              over ice
!
!  pahfl    : surface latent heat flux   (accumulated)
!  pahfll   :             "              over land
!  pahflw   :             "              over water
!  pahfli   :             "              over ice
!
!  pevap    : surface evaporation (accumulated)
!  pevapl   :             "        over land
!  pevapw   :             "        over water
!  pevapi   :             "        over ice
!
!  paz0     : roughness length for momentum
!  paz0l    :      "            over land
!  paz0w    :      "            over water
!  paz0i    :      "            over ice
!
!  pustr    : u-stress (accumulated)
!  pustrl   :     "     over land
!  pustrw   :     "     over sea
!  pustri   :     "     over ice
!
!  pvstr    : v-stress (accumulated)
!  pvstrl   :     "     over land
!  pvstrw   :     "     over water
!  pvstri   :     "     over ice
!
!  pdew2    : dew point temperature at 2 meter
!  peforest : forest coverage
!  psn      : snow depth
!  psnc     : snow depth on canopy
!  ptemp2   : temperature at 2 meter
!  ptsm1    : surface temperature (t-dt)
!  pt2max   : maximum temp. at 2 m between output intervals
!  pt2min   : minimun temp. at 2 m between output intervals
!  pwind10w : 10m wind over water
!  pu10     : u-wind at 10 meter
!  pv10     : v-wind at 10 meter
!  pwind10  : wind speed at 10 meter (accumulated)
!  pwimax   : maximum windspeed at 10 m. between output intervals
!  pvdis    : boundary layer dissipation (accumulated)
!  pws      : surface soil wetness
!  pwsmx    : field capacity of soil
!  pwlmx    : skin reservoir
!  pvlt     : leaf area index
!  pvgrat   : vegetation ratio
!
! - 2d within physics only
!
!  paphm1   : half level pressure (t-dt)
!  papm1    : full level pressure (t-dt)
!  ptvm1    : virtual temperature at t-dt
!  pvdiffp  : rate of change of qv due to vdiff routine for cover
!  pvmixtau : vdiff mixing timescale for variance and skewness
!
! - 1d within physics only
!
!  pgeom1   : geopotential above surface (t-dt)
!  psrfl    : net solar radiative flux at the surface
!  pqhfla   : moisture flux at the surface
!  pevapot  : potential evaporation
!  pcvs     : fractional snow cover (defined in *physc*)
!  pcvw     : wet skin fraction
!  ktropo   : tropopause index
!  lpland   : land-sea flag
!
!        Tendencies
!
!  - 3d
!
!  pxtte    : tendencies of tracer variables
!
!  - 2d
!  pvol     : tendency of meridional wind
!  pvom     : tendency of zonal wind
!  pqte     : tendency of humidity
!  ptte     : tendency of temperature
!  pxlte    : tendency of cloud water
!  pxite    : tendency of cloud ice
!
!
!     Method.
!     -------
!
!        First an auxialiary variable cp(q)t+gz is created on which
!   the vertical diffusion process will work like on u,v and q. then
!   along the vertical and at the surface, exchange coefficients (with
!   the dimension of a pressure thickness) are computed for momentum
!   and for heat (sensible plus latent). the letters m and h are used
!   to distinguish them. the diffusioncoefficents depend on the
!   turbulent kinetic energy (tke) calculated by an additional
!   prognostic equation, which considers advection of tke.
!        In the second part of the routine the implicit linear
!   systems for u,v first and t,q second are solved by a *gaussian
!   elimination back-substitution method. for t and q the lower
!   boundary condition depends on the surface state.
!   for tke the lower boundary condition depends on the square of
!   the frictional velocity.
!   over land, two different regimes of evaporation prevail:
!   a stomatal resistance dependent one over the vegetated part
!   and a soil relative humidity dependent one over the
!   bare soil part of the grid mesh.
!   potential evaporation takes place over the sea, the snow
!   covered part and the liquid water covered part of the
!   grid mesh as well as in case of dew deposition.
!        Finally one returns to the variable temperature to compute
!   its tendency and the later is modified by the dissipation's effect
!   (one assumes no storage in the turbulent kinetic energy range) and
!   the effect of moisture diffusion on cp. z0 is updated and the
!   surface fluxes of t and q and their derivatives are prepared and
!   stored like the difference between the implicitely obtained
!   cp(q)t+gz and cp(q)t at the surface.
!
!
!     Reference.

!
!          See vertical diffusion's part of the model's documentation
!     for details about the mathematics of this routine.
!
!     Authors.
!
!     u. schlese     dkrz-hamburg  feb-93
!       modified     e. roeckner  - 1994
!
!     j.-p. schulz   mpi - 1997 : implementation of implicit
!                                 coupling between land surface
!                                 and atmosphere.
!     m. esch, mpi, june 1999, echam5-modifications
!
!
!     based  on  original ecmwf version by j.f. geleyn  - 1982
!                              modified by c.b. blondin - 1986
!                                          h. feichter  - 1991
!                                          s. brinkop   - 1992
!                                          m. claussen  - 1993
USE mo_kind,             ONLY: dp
!!$#ifndef MESSY
USE mo_geoloc,           ONLY: coriol_2d
!!$#else
!!$USE mo_geoloc,           ONLY: coriol_2d,philon_2d,philat_2d
!!$#endif
USE mo_param_switches,   ONLY: lvdiff
USE mo_physc2,           ONLY: clam, ckap, cb, cc, cchar, cvdifts,     &
                               cfreec, cgam, cz0ice, csncri
USE mo_constants,        ONLY: vtmpc1, cpd, rd, g, vtmpc2, tmelt, alv, &
                               als, api, rhoh2o, stbo, c3les, c3ies,   &
                               c4les, c4ies, c2es
USE mo_vegetation,       ONLY: cva, cvb, cvc, cvbc, cvk, cvkc, cvabc,  &
                               cvrad
!!$#ifndef MESSY
USE mo_tracer,           ONLY: trlist
!!$#else
!!$! mz_rs_20040329+
!!$  USE messy_main_tracer_mem_bi, ONLY: ntrac_gp, ti_gp, ON, I_Vdiff
!!$! mz_rs_20040329-
!!$#ifdef MESSYTENDENCY
!!$  USE messy_main_tendency_bi,   ONLY: mtend_get_handle,               &
!!$                                      mtend_get_start_l,              &
!!$                                      mtend_add_l,                    &
!!$                                      mtend_id_t,                     &
!!$                                      mtend_id_q,                     &
!!$                                      mtend_id_xl,                    &
!!$                                      mtend_id_xi,                    &
!!$                                      mtend_id_u,                     &
!!$                                      mtend_id_v,                     &
!!$                                      mtend_id_tracer
!!$#endif
!!$#endif
USE mo_convect_tables,   ONLY: tlucua, jptlucu1, jptlucu2,             &
                               lookuperror, lookupoverflow
!!$#ifndef MESSY
USE mo_radiation,        ONLY: cemiss
!!$#else
!!$USE messy_main_constants_mem, ONLY: cemiss
!!$#endif
USE mo_exception,        ONLY: finish
USE mo_time_control,     ONLY: delta_time, lstart, time_step_len
USE mo_semi_impl,        ONLY: eps
!
#ifdef __ibm__
USE mo_specfun,          ONLY: merge
#endif
!!$#ifdef MESSY
!!$! mz_rs_20020911+
!!$USE messy_main_data_bi, ONLY: pxtems, &
!!$                              cfml,  cfmw,  cfmi, cfncl, &
!!$                              cfncw, cfnci, cdnl, cdnw, cdni, &
!!$                              riw,   rii,   ril,  tvir, tvl,  &
!!$                              tvw,   tvi,   srfl, fws,        &
!!$                              rco_leaf, rh_2m, wind10_2d,     &
!!$                              wind10w_2d,                & ! mz_ap_20070913
!!$                              rinum_3d,tpot_3d,          & ! mz_mt_20031009
!!$                              zust_2d ,zsenkf_2d ,       & ! mz_mt_20031119
!!$                              zlatkf_2d , zcdh_2d,       & ! mz_mt_20031119
!!$                              vmixtau, vdiffp,           & ! mz_ht_20041025
!!$                              s_heatflux,                & ! mz_ht_20050222
!!$                              l_heatflux,                & ! mz_ht_20050222
!!$                              qsl, qsw, qsi, phum, chl,  & ! op_re_20140904
!!$                              cfhl, cfhw, cfhi,          & ! op_re_20140904
!!$                              cfh, ebsh,                 & ! op_re_20140904
!!$                              qslnew, wet_tmp, cair        ! op_re_20140904
!!$
!!$! mz_rs_20020911-
!!$! mz_lg_20040603+, mz_ht_20040603+ added
!!$USE messy_main_data_bi, ONLY: pxtte => xtte
!!$! mz_lg_20040603-, mz_ht_20040603-
!!$#endif
!
IMPLICIT NONE
!
! Variables declared
!
  INTEGER :: it, itop, itopp1, jk, jl, jt, klev, klevm1, klevp1
  INTEGER :: kproma, kbdim, ktdia, ktrac, it1, krow
  REAL(dp):: z0h, z1dgam, z2geomf, zabcs, zalf, zalh2, zalo, zaloh
  REAL(dp):: zaph2m, zb, zbet, zblend, zbuoy, zc
  REAL(dp):: zcbn, zcbs, zcbu, zcdn2m, zcdnr, zcfm2m, zchar, zchneu
  REAL(dp):: zchsnow, zchland, zchmean, zepot, zsnfac, zwlfac
  REAL(dp):: zcoefi, zcoefl, zcoefw, zcons, zcons11, zcons12
  REAL(dp):: zcons13, zcons14, zcons15, zcons16, zcons17, zcons18
  REAL(dp):: zcons2, zcons23, zcons25, zcons29, zcons3, zcons30, zcons5
  REAL(dp):: zcons6, zcons8, zcons9, zconvs, zcor, zcpd, zcpt, zcvm3
  REAL(dp):: zcvm4, zda1, zdew2i, zdew2l, zdew2w
  REAL(dp):: zdisc, zdisci, zdiscl, zdiscw, zdisl
  REAL(dp):: zdisqi, zdisql, zdisqw, zdisx, zdisxt, zdivv, zdivv1
  REAL(dp):: zdqdt, zdqtot, zds, zdtdt, zdthv, zdudt, zdus1
  REAL(dp):: zdus2, zdvdt, zdximdt, zdxlmdt, zdxtdt, zdz
  REAL(dp):: zepcor, zepdu2, zepevap, zephum, zeps, zepsec, zepshr
  REAL(dp):: zepsr, zepz0o, zepzzo, zes, zfac, zfox, zfrac
  REAL(dp):: zfreec, zfux, zgam, zghabl, zh1, zh2, zh2m, zhexp, zhtq
  REAL(dp):: zgtl, zgtw, zgti, zgtsum, zgql, zgqw, zgqi, zgqsum
  REAL(dp):: zhuv, zkap, zkappa, zktest, zlam, zln1, zln2, zm1, zm2
  REAL(dp):: zm4, zmix, zmonob, zmult1, zmult2, zmult3, zmult4
  REAL(dp):: zmult5, zplmax, zplmin, zq2m, zqddif, zqdp, zqklevi
  REAL(dp):: zqklevl, zqklevw, zqlwi1, zqlwi2, zqmitte, zqnlev
  REAL(dp):: zqs1, zqs2, zqsmit, zqst1, zqtmit, zqvhfl, zqwevap, zrat, zrd
  REAL(dp):: zsdep2, zsh, zshear, zshn, zsm, zsmn, zsoil, zspeedi
  REAL(dp):: zrdrv, zred, zrh2m, zri, zrsi, zrvrd, zsdep1
  REAL(dp):: zspeedl, zspeedw, zsrfl, zsrfld, zstabf, zt2i, zt2l
  REAL(dp):: zt2w, zteff4, zteldif, ztemitte, ztest, ztkemin
  REAL(dp):: ztkesq, ztkevi, ztkevl, ztkevw, ztklevi, ztklevl
  REAL(dp):: ztklevw, ztmit, ztmst, ztnlev, ztpfac1, ztpfac2
  REAL(dp):: ztpfac3, ztpfac4, zqnew, ztrfll, zu10i, zu10l, zu10w
  REAL(dp):: zucf, zust, zustarm, zustf, zusti, zustl, zustw, zusus1
  REAL(dp):: zv10i, zv10l, zv10w, zva, zvabc, zvb, zvbc, zvc
  REAL(dp):: zvirmitte, zvk, zvkc, zvklt, zvrad, zvxmklt, zvxpklt
  REAL(dp):: zwcrit, zwpwp, zwslev, zwstf, zwstop, zz2geo
  REAL(dp):: zzb, zzcpts, zzqs, zztvm, zzzlam
  REAL(dp):: zcptlcorr, zwslim, ztrdown, zlai, zdtime
  REAL(dp):: zqsurf, zrhodz
  REAL(dp):: ztvlan, ztvsea, ztvice, ztvh
  REAL(dp):: zdisv, zfav

!     ------------------------------------------------------------------
!
  REAL(dp):: pxtm1(kbdim,klev,ktrac)
!
  REAL(dp)::                                                           &
       pqm1(kbdim,klev),     ptm1(kbdim,klev),     pum1(kbdim,klev)    &
      ,pvm1(kbdim,klev),     pxlm1(kbdim,klev),    pxim1(kbdim,klev)   &
      ,pxvar(kbdim,klev)
!!$#ifdef MESSY
!!$! mz_ht_20041209+
!!$! these are the real t-dt values
!!$  REAL(dp) ::                                                             &
!!$       pqmm1(kbdim,klev),    ptmm1(kbdim,klev),    pumm1(kbdim,klev)      &
!!$      ,pvmm1(kbdim,klev),    pxlmm1(kbdim,klev),   pximm1(kbdim,klev)
!!$! mz_ht_20041209-
!!$#endif
!!$#ifdef MESSYTENDENCY
!!$   REAL(dp),DIMENSION(kbdim,klev)          :: zdtdt_2d
!!$   REAL(dp),DIMENSION(kbdim,klev)          :: zdqdt_2d
!!$   REAL(dp),DIMENSION(kbdim,klev)          :: zdudt_2d
!!$   REAL(dp),DIMENSION(kbdim,klev)          :: zdvdt_2d
!!$   REAL(dp),DIMENSION(kbdim,klev)          :: zdxlmdt_2d
!!$   REAL(dp),DIMENSION(kbdim,klev)          :: zdximdt_2d
!!$   REAL(dp),DIMENSION(kbdim,klev,ntrac_gp) :: zdxtdt_3d
!!$   LOGICAL,SAVE                            :: lfirst=.true.
!!$   INTEGER,SAVE                            :: my_handle
!!$#endif
  REAL(dp)::                                                           &
       pahfl(kbdim),         pahfs(kbdim),         paz0(kbdim)         &
      ,pdew2(kbdim),         pevap(kbdim),         pforest(kbdim)      &
      ,ptemp2(kbdim),        pt2max(kbdim)                             &
      ,pt2min(kbdim),        pwind10w(kbdim),      pvdis(kbdim)        &
      ,pu10(kbdim),          pv10(kbdim),          pustr(kbdim)        &
      ,pvstr(kbdim),         pwimax(kbdim),        pwind10(kbdim)      &
      ,pwsmx(kbdim),         pvlt(kbdim),          pgrndcapc(kbdim)    &
      ,pgrndhflx(kbdim),     pvgrat(kbdim)                             &
      ,ptsl(kbdim),          ptsw(kbdim),          ptsi(kbdim)         &
      ,pocu(kbdim),          pocv(kbdim)                               &
      ,paz0l(kbdim),         paz0w(kbdim),         paz0i(kbdim)        &
      ,pahfsl(kbdim),        pahfsw(kbdim),        pahfsi(kbdim)       &
      ,pahfll(kbdim),        pahflw(kbdim),        pahfli(kbdim)       &
      ,pevapl(kbdim),        pevapw(kbdim),        pevapi(kbdim)       &
      ,pahfslac(kbdim),      pahfswac(kbdim),      pahfsiac(kbdim)     &
      ,pahfllac(kbdim),      pahflwac(kbdim),      pahfliac(kbdim)     &
      ,pevaplac(kbdim),      pevapwac(kbdim),      pevapiac(kbdim)     &
      ,pustrl(kbdim),        pustrw(kbdim),        pustri(kbdim)       &
      ,pvstrl(kbdim),        pvstrw(kbdim),        pvstri(kbdim)       &
      ,psn(kbdim),           psnc(kbdim),          ptslm1(kbdim)       &
      ,pws(kbdim),           palbedo(kbdim),       palsol(kbdim)
  REAL(dp)::                                                           &
       ptke(kbdim,klev),     ptkem1(kbdim,klev),   ptkem(kbdim,klev)   &
      ,paclc(kbdim,klev),    pemter(kbdim,klevp1)
  REAL(dp)::                                                           &
       paphm1(kbdim,klevp1), papm1(kbdim,klev),    pgeom1(kbdim,klev)  &
!!$#ifndef MESSY
      ,ptvm1(kbdim,klev),    pvdiffp(kbdim,klev),  pvmixtau(kbdim,klev)
!!$#else
!!$      ,ptvm1(kbdim,klev)
!!$#endif
  REAL(dp)::                                                           &
       pcvs(kbdim),          pcvw(kbdim),          psrfl(kbdim)        &
      ,pqhfla(kbdim),        pevapot(kbdim)                            &
      ,ptslnew(kbdim),       pwlmx(kbdim)
  LOGICAL ::                lpland(kbdim),        lpglac(kbdim)
  REAL(dp)::                                                           &
       pfrl(kbdim),          pfrw(kbdim),          pfri(kbdim)
!
!!$#ifndef MESSY
  REAL(dp)::                                                           &
       pxtte(kbdim,klev,ktrac)
!!$#endif
!
  REAL(dp)::                                                           &
       pvol(kbdim,klev),     pvom(kbdim,klev),     pqte(kbdim,klev)    &
      ,ptte(kbdim,klev),     pxlte(kbdim,klev),    pxite(kbdim,klev)
!
!     local variables
!
  LOGICAL :: lo,lo1
!
  REAL(dp):: zxtvn(kbdim,klev,ktrac) ! op_ck_20031001
!!$#ifndef MESSY
  REAL(dp)::                                                           &
          zxtdif(kbdim,klev,ktrac), zxtems(kbdim,ktrac)
  REAL(dp)::                                                           &
          zcfm(kbdim,klev),    zdis(kbdim,klev)                        &
         ,zcfh(kbdim,klev),    zcptgz(kbdim,klev),   zebsm(kbdim,klev) &
         ,zudif(kbdim,klev),   zvdif(kbdim,klev)                       &
         ,zwet(kbdim),         zdqsl(kbdim),         zsrfll(kbdim)     &
         ,ztcoe(kbdim),        zchnw(kbdim),         zcfnchw(kbdim)
  REAL(dp)::                                                           &
          zcr(kbdim),          zhum(kbdim)                             &
         ,zcsat(kbdim),        zcair(kbdim),         ztdif(kbdim,klev) &
         ,zqdif(kbdim,klev),   zebsh(kbdim,klev),    zvidis(kbdim)     &
         ,z1mxtm1(kbdim)
  REAL(dp)::                                                           &
          zhdyn(kbdim),        zteta1(kbdim,klev)                      &
         ,zlteta1(kbdim,klev)                                          &
         ,ztvir1(kbdim,klev),  zhh(kbdim,klevm1),    zqss(kbdim,klev)  &
         ,zxldif(kbdim,klev),  zxidif(kbdim,klev),   zedif(kbdim,klev) &
         ,ztkevn(kbdim,klev),  zx(kbdim,klev),       zhsoil(kbdim)     &
         ,zvardif(kbdim,klev)
  REAL(dp)::                                                           &
          zqssm(kbdim,klevm1), ztmitte(kbdim,klevm1)                   &
         ,zqmit(kbdim,klevm1), ztvirmit(kbdim,klevm1)                  &
         ,zfaxen(kbdim,klevm1),zfaxe(kbdim,klev)                       &
         ,ztemit(kbdim,klevm1),zccover(kbdim,klevm1)                   &
         ,zcdum(kbdim,klev),   zlwcmit(kbdim,klevm1)                   &
         ,zcfv(kbdim,klev),    zebsv(kbdim,klev)
!
  REAL(dp)::                                                           &
          zqsl(kbdim),         zqsw(kbdim),          zqsi(kbdim)       &
         ,ztesl(kbdim),        ztesw(kbdim),         ztesi(kbdim)      &
         ,ztvl(kbdim),         ztvw(kbdim),          ztvi(kbdim)       &
         ,zril(kbdim),         zriw(kbdim),          zrii(kbdim)       &
         ,zucfl(kbdim),        zucfw(kbdim),         zucfi(kbdim)      &
         ,zscfl(kbdim),        zscfw(kbdim),         zscfi(kbdim)      &
         ,zcfncl(kbdim),       zcfncw(kbdim),        zcfnci(kbdim)     &
         ,zcfml(kbdim),        zcfmw(kbdim),         zcfmi(kbdim)      &
         ,zcfhl(kbdim),        zcfhw(kbdim),         zcfhi(kbdim)      &
         ,zcdnl(kbdim),        zcdnw(kbdim),         zcdni(kbdim)      &
         ,zwstl(kbdim),        zwstw(kbdim),         zwsti(kbdim)      &
         ,zchl(kbdim),         zchw(kbdim),          zchi(kbdim)       &
         ,zbnl(kbdim),         zbnw(kbdim),          zbni(kbdim)       &
         ,zbml(kbdim),         zbmw(kbdim),          zbmi(kbdim)       &
         ,zbhl(kbdim),         zbhw(kbdim),          zbhi(kbdim)       &
         ,zustarl(kbdim),      zustarw(kbdim),       zustari(kbdim)    &
         ,zcptl(kbdim),        zcptw(kbdim),         zcpti(kbdim)      &
         ,zqhfll(kbdim),       zqhflw(kbdim),        zqhfli(kbdim)     &
         ,zthfll(kbdim),       zthflw(kbdim),        zthfli(kbdim)     &
         ,zchnl(kbdim),        zcfnchl(kbdim),       zbhnl(kbdim)      &
         ,zucfhl(kbdim),       zdu2(kbdim),          zdu2oc(kbdim)

  REAL(dp)::                                                           &
          zqshear(kbdim,klev), zrho(kbdim,klevp1)                      &
         ,zqflux(kbdim,klevp1),zvarpr(kbdim,klevp1)
!!$#else
!!$  REAL(dp)::                                                           &
!!$          zxtdif(kbdim,klev,ktrac)
!!$  REAL(dp)::                                                           &
!!$          zcfm(kbdim,klev),    zdis(kbdim,klev)                        &
!!$!op_re_20140904+
!!$! !!$         ,zcfh(kbdim,klev),    zcptgz(kbdim,klev),   zebsm(kbdim,klev) &
!!$         ,zcptgz(kbdim,klev),   zebsm(kbdim,klev) &
!!$!op_re_20140904-
!!$         ,zudif(kbdim,klev),   zvdif(kbdim,klev)                       &
!!$         ,zwet(kbdim),         zdqsl(kbdim)                            &
!!$         ,ztcoe(kbdim),        zchnw(kbdim),         zcfnchw(kbdim)
!!$!op_re_20140904+                                                       &
!!$! !!$  REAL(dp)::                                                           &  
!!$! !!$          zcr(kbdim),          zhum(kbdim)                             &
!!$! !!$         ,zcsat(kbdim),        zcair(kbdim),         ztdif(kbdim,klev) &
!!$! !!$         ,zqdif(kbdim,klev),   zebsh(kbdim,klev),    zvidis(kbdim)     &
!!$  REAL(dp):: &
!!$          zcr(kbdim)                                                   &
!!$         ,zcsat(kbdim),                              ztdif(kbdim,klev) &
!!$         ,zqdif(kbdim,klev),                         zvidis(kbdim)     &
!!$!op_re_20140904-
!!$         ,z1mxtm1(kbdim)
!!$  REAL(dp)::                                                           &
!!$          zhdyn(kbdim)                                                 &
!!$         ,zlteta1(kbdim,klev)                                          &
!!$         ,ztvir1(kbdim,klev),  zhh(kbdim,klevm1),    zqss(kbdim,klev)  &
!!$         ,zxldif(kbdim,klev),  zxidif(kbdim,klev),   zedif(kbdim,klev) &
!!$         ,ztkevn(kbdim,klev),  zx(kbdim,klev),       zhsoil(kbdim)     &
!!$         ,zvardif(kbdim,klev)
!!$  REAL(dp)::                                                           &
!!$          zqssm(kbdim,klevm1), ztmitte(kbdim,klevm1)                   &
!!$         ,zqmit(kbdim,klevm1), ztvirmit(kbdim,klevm1)                  &
!!$         ,zfaxen(kbdim,klevm1),zfaxe(kbdim,klev)                       &
!!$         ,ztemit(kbdim,klevm1),zccover(kbdim,klevm1)                   &
!!$         ,zcdum(kbdim,klev),   zlwcmit(kbdim,klevm1)                   &
!!$         ,zcfv(kbdim,klev),    zebsv(kbdim,klev)
!!$!
!!$!op_re_20140904+
!!$! !!$  REAL(dp)::                                                           &
!!$! !!$          zqsl(kbdim),         zqsw(kbdim),          zqsi(kbdim)       &
!!$!op_re_20140904-
!!$  REAL(dp)::                                                           &
!!$          ztesl(kbdim),        ztesw(kbdim),         ztesi(kbdim)      &
!!$         ,zucfl(kbdim),        zucfw(kbdim),         zucfi(kbdim)      &
!!$         ,zscfl(kbdim),        zscfw(kbdim),         zscfi(kbdim)      &
!!$!op_re_20140904+
!!$! !!$         ,zcfhl(kbdim),        zcfhw(kbdim),         zcfhi(kbdim)      &
!!$!op_re_20140904-
!!$         ,zwstl(kbdim),        zwstw(kbdim),         zwsti(kbdim)      &
!!$!op_re_20140904+
!!$! !!$         ,zchl(kbdim),         zchw(kbdim),          zchi(kbdim)       &
!!$         ,                     zchw(kbdim),          zchi(kbdim)       &
!!$!op_re_20140904-
!!$         ,zbnl(kbdim),         zbnw(kbdim),          zbni(kbdim)       &
!!$         ,zbml(kbdim),         zbmw(kbdim),          zbmi(kbdim)       &
!!$         ,zbhl(kbdim),         zbhw(kbdim),          zbhi(kbdim)       &
!!$         ,zustarl(kbdim),      zustarw(kbdim),       zustari(kbdim)    &
!!$         ,zcptl(kbdim),        zcptw(kbdim),         zcpti(kbdim)      &
!!$         ,zqhfll(kbdim),       zqhflw(kbdim),        zqhfli(kbdim)     &
!!$         ,zthfll(kbdim),       zthflw(kbdim),        zthfli(kbdim)     &
!!$         ,zchnl(kbdim),        zcfnchl(kbdim),       zbhnl(kbdim)      &
!!$         ,zucfhl(kbdim),       zdu2(kbdim),          zdu2oc(kbdim)
!!$
!!$  REAL(dp)::                                                           &
!!$          zqshear(kbdim,klev), zrho(kbdim,klevp1)                      &
!!$         ,zqflux(kbdim,klevp1),zvarpr(kbdim,klevp1)
!!$
!!$  ! CHANGES:
!!$  REAL(dp), POINTER, DIMENSION(:,:) :: zxtems
!!$  REAL(dp), POINTER, DIMENSION(:) :: zcfml  
!!$  REAL(dp), POINTER, DIMENSION(:) :: zcfmw  
!!$  REAL(dp), POINTER, DIMENSION(:) :: zcfmi  
!!$  REAL(dp), POINTER, DIMENSION(:) :: zcfncl 
!!$  REAL(dp), POINTER, DIMENSION(:) :: zcfncw 
!!$  REAL(dp), POINTER, DIMENSION(:) :: zcfnci 
!!$  REAL(dp), POINTER, DIMENSION(:) :: zcdnl  
!!$  REAL(dp), POINTER, DIMENSION(:) :: zcdnw  
!!$  REAL(dp), POINTER, DIMENSION(:) :: zcdni  
!!$  REAL(dp), POINTER, DIMENSION(:) :: zril   
!!$  REAL(dp), POINTER, DIMENSION(:) :: zriw   
!!$  REAL(dp), POINTER, DIMENSION(:) :: zrii   
!!$  REAL(dp), POINTER, DIMENSION(:) :: ztvl   
!!$  REAL(dp), POINTER, DIMENSION(:) :: ztvw   
!!$  REAL(dp), POINTER, DIMENSION(:) :: ztvi   
!!$  REAL(dp), POINTER, DIMENSION(:) :: zsrfll 
!!$
!!$  ! op_re_20140904+
!!$  ! Variables converted to pointers for H2OISO
!!$  REAL(dp), POINTER, DIMENSION(:) :: zqsl
!!$  REAL(dp), POINTER, DIMENSION(:) :: zqsi
!!$  REAL(dp), POINTER, DIMENSION(:) :: zqsw
!!$  REAL(dp), POINTER, DIMENSION(:) :: zhum
!!$  REAL(dp), POINTER, DIMENSION(:) :: zchl
!!$  REAL(dp), POINTER, DIMENSION(:) :: zcfhl
!!$  REAL(dp), POINTER, DIMENSION(:) :: zcfhw
!!$  REAL(dp), POINTER, DIMENSION(:) :: zcfhi
!!$  REAL(dp), POINTER, DIMENSION(:) :: zqslnew
!!$  REAL(dp), POINTER, DIMENSION(:) :: zcair
!!$  !re Additional temp result for H2OISO
!!$  REAL(dp), POINTER, DIMENSION(:) :: zwet_tmp
!!$
!!$  ! Variables converted to pointers for H2OISO
!!$  REAL(dp), POINTER, DIMENSION(:,:) :: zcfh
!!$  REAL(dp), POINTER, DIMENSION(:,:) :: zebsh
!!$  ! op_re_20140904-
!!$
!!$  REAL(dp), POINTER, DIMENSION(:,:) :: zteta1   
!!$  REAL(dp), POINTER, DIMENSION(:,:) :: pvdiffp  
!!$  REAL(dp), POINTER, DIMENSION(:,:) :: pvmixtau 
!!$#endif
!
  INTEGER :: ihpbl(kbdim),  ihpblc(kbdim),   ihpbld(kbdim)
!!$#ifdef MESSY
!!$  REAL    :: zcdm(kbdim)                      ! mz_mt_20031120
!!$  REAL    :: dens        ! density [kg m-3]   ! mz_mt_20031120
!!$  REAL    :: pahfsm1(kbdim)                   ! mz_mt_20031120
!!$#endif
!
!     ------------------------------------------------------------------
!
!     THE FOLLOWING VARIABLES ARE NEEDED FOR THE SUBROUTINES SURFTEMP
!
  REAL(dp)::    zetnl(kbdim)    ! richtmyer-morton-coefficients
  REAL(dp)::    zetnw(kbdim)    !             "
  REAL(dp)::    zetni(kbdim)    !             "
  REAL(dp)::    zftnl(kbdim)    ! for dry static energy.
  REAL(dp)::    zftnw(kbdim)    !
  REAL(dp)::    zftni(kbdim)    !
  REAL(dp)::    zeqnl(kbdim)    ! richtmyer-morton-coefficients
  REAL(dp)::    zeqnw(kbdim)    !
  REAL(dp)::    zeqni(kbdim)    !
  REAL(dp)::    zfqnl(kbdim)    ! for moisture.
  REAL(dp)::    zfqnw(kbdim)    !
  REAL(dp)::    zfqni(kbdim)    !
  REAL(dp)::    zcpq(kbdim)     ! specific heat of air as used in *vdiff*
  REAL(dp)::    znetr(kbdim)    ! surface net radiation (old)
  REAL(dp)::    zcdrag(kbdim)   ! drag coefficient for heat and moisture
  REAL(dp)::    zcptlnew(kbdim) ! new surface dry static energy
!!$#ifndef MESSY
  REAL(dp)::    zqslnew(kbdim)  ! new surf. sat. spec. humidity
!!$#endif

!!$#ifdef MESSY
!!$  REAL(dp)::    t2l_1d(kbdim)
!!$  REAL(dp)::    t2w_1d(kbdim)
!!$  REAL(dp)::    t2i_1d(kbdim)
!!$#endif

!  Executable statements 

!!$#ifdef MESSYTENDENCY
!!$  IF(lfirst) THEN
!!$     my_handle=mtend_get_handle('vdiff')
!!$     lfirst=.FALSE.
!!$  ENDIF
!!$  zdtdt_2d(:,:)    = 0.0_dp
!!$  zdqdt_2d(:,:)    = 0.0_dp
!!$  zdudt_2d(:,:)    = 0.0_dp
!!$  zdvdt_2d(:,:)    = 0.0_dp
!!$  zdxlmdt_2d(:,:)  = 0.0_dp
!!$  zdximdt_2d(:,:)  = 0.0_dp
!!$  zdxtdt_3d(:,:,:) = 0.0_dp
!!$#endif

  lookupoverflow = .FALSE.
!!$#ifdef MESSY
!!$  zxtems   => pxtems(:,1,:,krow) ! mz_pj_20040401
!!$  ! mz_ht_20041026+
!!$  ! 2D
!!$  zcfml    => cfml(:,krow)
!!$  zcfmw    => cfmw(:,krow)
!!$  zcfmi    => cfmi(:,krow)
!!$  zcfncl   => cfncl(:,krow)
!!$  zcfncw   => cfncw(:,krow)
!!$  zcfnci   => cfnci(:,krow)
!!$  zcdnl    => cdnl(:,krow)
!!$  zcdnw    => cdnw(:,krow)
!!$  zcdni    => cdni(:,krow)
!!$  zril     => ril(:,krow)
!!$  zriw     => riw(:,krow)
!!$  zrii     => rii(:,krow)
!!$  ztvl     => tvl(:,krow)
!!$  ztvw     => tvw(:,krow)
!!$  ztvi     => tvi(:,krow)
!!$  zsrfll   => srfl(:,krow)
!!$  ! op_re_20140904+
!!$  !Variables converted to pointers for H2OISO
!!$  zqsl     => qsl(1:kproma,krow)
!!$  zqsw     => qsw(1:kproma,krow)
!!$  zqsi     => qsi(1:kproma,krow)
!!$  zhum     => phum(1:kproma,krow)
!!$  zchl     => chl(1:kproma,krow)
!!$  zcfhl    => cfhl(1:kproma,krow)
!!$  zcfhw    => cfhw(1:kproma,krow)
!!$  zcfhi    => cfhi(1:kproma,krow)
!!$  zqslnew  => qslnew(1:kproma,krow)
!!$  zcair    => cair(1:kproma,krow)
!!$  !re Additional pointer to store temp result for H2OISO
!!$  zwet_tmp => wet_tmp(1:kproma,krow)
!!$  ! op_re_20140904-
!!$  ! 3D
!!$  zteta1   => tpot_3d(:,:,krow)
!!$  pvmixtau => vmixtau(:,:,krow) 
!!$  pvdiffp  => vdiffp(:,:,krow)  
!!$  ! mz_ht_20041026-
!!$  ! mz_ht_20050222+
!!$  s_heatflux(1:kproma,krow) = 0._dp
!!$  ! mz_ht_20050222-
!!$
!!$  ! op_re_20140904+
!!$  !Variables converted to pointers for H2OISO
!!$  zcfh     => cfh(1:kproma,:,krow)
!!$  zebsh    => ebsh(1:kproma,:,krow)
!!$  ! op_re_20140904-
!!$#endif

!
!     ------------------------------------------------------------------
!
!*    PHYSICAL CONSTANTS.
!
!          *ZLAM* IS THE ASYMPTOTIC MIXING LENGTH FOR MOMENTUM EXCHANGE,
!     *ZKAP* IS THE VON KARMAN CONSTANT, *ZB*, *ZC* AND *ZD* ARE SOME
!     CONSTANTS FOR THE FORMULAE ABOUT STABILITY DEPENDENCY RESPECTIVELY
!     NEAR THE NEUTRAL CASE, IN THE UNSTABLE CASE AND IN THE STABLE
!     CASE AND *ZCHAR* IS THE CONSTANT OF THE *CHARNOCK FORMULA.
!     *ZQWSSAT* AND *ZQSNCR* ARE THE INVERSES OF CRITICAL VALUES FOR
!     SOIL WATER AND SNOW DEPTH THAT ARE USED IN THE COMPUTATION OF THE
!     EVAPOTRANSPIRATION'S EFFICIENCY.
!
  zlam=clam
  zkap=ckap
  zb=cb
  zc=cc
  zchar=cchar
  zva=cva
  zvb=cvb
  zvc=cvc
  zvbc=cvbc
  zvk=cvk
  zvkc=cvkc
  zvabc=cvabc
  zvrad=cvrad
  zrvrd=vtmpc1+1._dp
  zrdrv=1._dp/zrvrd
  zcpd=cpd
  zrd=rd
  zkappa=zrd/zcpd
!
!*      PARAMETERS FOR BOUNDARY LAYER DIAGNOSTICS
!
  zhuv=10._dp*g
  zhtq=2._dp*g
  zephum=5.e-2_dp
!
!*    SECURITY PARAMETERS.
!
  zepdu2=1.0_dp
  zepshr=1.e-5_dp
  zepzzo=1.5e-05_dp
  zepz0o=2._dp
  zepcor=5.e-05_dp
  ztkemin=1.e-10_dp
  zepsr=1.e-10_dp
  zepevap=1.e-10_dp
  zepsec=1.e-2_dp
!
!*    COMPUTATIONAL CONSTANTS.
!
  zdtime = delta_time
  ztmst  = time_step_len
  ztpfac1=cvdifts
  ztpfac2=1._dp/ztpfac1
  ztpfac3=1._dp-ztpfac2
  ztpfac4=1._dp+ztpfac3
  zzzlam=1._dp
  zcons2=0.5_dp*zkap/g
  zcons3=zlam
  zcons5=3._dp*zb*zc*g**2
  zcons6=1._dp/3._dp
  zcons8=2._dp*zb
  zcons9=3._dp*zb
  zcons11=3._dp*zb*zc
  zcons12=ztpfac1*ztmst*g/rd
  zcons13=1._dp/ztmst
  zcons14=zchar*rd/(g**2*ztmst)
  zcons15=1._dp/(g*ztmst)
  zcons16=cpd*vtmpc2
  zcons18=ztpfac1*ztmst*g**2
  zcons17=1._dp/zkap**2
  zcons25=zcons2/zcons3
  zcons29=ztpfac1*ztmst
  zcons30=1._dp/(ztpfac1*ztmst*g)
  zplmax=0.75_dp
  zplmin=0.35_dp
  zwslim=zplmin
  zblend=100._dp
  zchneu=.3_dp
  zfreec=cfreec
  zgam=cgam
  z1dgam=1._dp/zgam
  zh1= 2.22_dp
  zh2= 0.22_dp
  zm1= 1.24_dp
  zm2= 2.37_dp
  zm4= 3.69_dp
  zshn=zh1*zh2*SQRT(2._dp)
  zsmn=zshn*zm1*zm2/zm4
  zda1=1._dp/zsmn**3
  zustf=1._dp/zsmn**2
  zwstf=0.2_dp
!
  itop=1
  itopp1=itop+1

!!$#ifdef MESSY
!!$! mz_ht_20041209+
!!$! now the actual value of temperature, humidity, wind, cloud water and
!!$! cloud ice are calculated and used all through the routine
!!$! ATTENTION:  the code reads still as if the minus 1 values (t-dt)
!!$! are used as in the original ECHAM5 code! However, with this addition
!!$! it is consistent with the leapfrog operator splitting integration scheme
!!$
!!$#ifndef MESSYTENDENCY
!!$  do jk=1,klev
!!$    ptm1(1:kproma,jk)  = ptmm1(1:kproma,jk)  + ptte(1:kproma,jk)  * ztmst
!!$    pqm1(1:kproma,jk)  = pqmm1(1:kproma,jk)  + pqte(1:kproma,jk)  * ztmst
!!$    pum1(1:kproma,jk)  = pumm1(1:kproma,jk)  + pvom(1:kproma,jk)  * ztmst
!!$    pvm1(1:kproma,jk)  = pvmm1(1:kproma,jk)  + pvol(1:kproma,jk)  * ztmst
!!$    pxlm1(1:kproma,jk) = pxlmm1(1:kproma,jk) + pxlte(1:kproma,jk) * ztmst
!!$    pxim1(1:kproma,jk) = pximm1(1:kproma,jk) + pxite(1:kproma,jk) * ztmst
!!$  enddo
!!$! mz_ht_20041209-
!!$#else
!!$  ! get compute start values 
!!$  call mtend_get_start_l (mtend_id_t,  v0 = ptm1)
!!$  call mtend_get_start_l (mtend_id_q,  v0 = pqm1)
!!$  call mtend_get_start_l (mtend_id_xl, v0 = pxlm1)  
!!$  call mtend_get_start_l (mtend_id_xi, v0 = pxim1)
!!$  call mtend_get_start_l (mtend_id_u,  v0 = pum1)  
!!$  call mtend_get_start_l (mtend_id_v,  v0 = pvm1) 
!!$#endif
!!$#endif
!
!     ------------------------------------------------------------------
!
!*         2.     NEW THERMODYNAMIC VARIABLE AND BOUNDARY CONDITIONS.
!
!*         2.1     REPLACE T BY CP(Q)*T+GZ IN THE ATMOSPHERE.
!
210 CONTINUE
  DO 212 jk=ktdia,klev
     DO 211 jl=1,kproma
        zx(jl,jk)=pxlm1(jl,jk)+pxim1(jl,jk)   ! total cloud water
        zcptgz(jl,jk)=pgeom1(jl,jk)+ptm1(jl,jk)                        &
                                   *cpd*(1._dp+vtmpc2*pqm1(jl,jk))
        zteta1(jl,jk)=ptm1(jl,jk)*(100000._dp/papm1(jl,jk))**zkappa
        ztvir1(jl,jk)=zteta1(jl,jk)*(1._dp+vtmpc1*pqm1(jl,jk)-zx(jl,jk))
        lo=ptm1(jl,jk).GE.tmelt
        zfaxe(jl,jk)=MERGE(alv,als,lo)
        zbet=zfaxe(jl,jk)/zcpd
        zusus1=zbet*zteta1(jl,jk)/ptm1(jl,jk)*zx(jl,jk)
        zlteta1(jl,jk)=zteta1(jl,jk)-zusus1
        it = NINT(ptm1(jl,jk)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zes=tlucua(it)/papm1(jl,jk)
        zes=MIN(zes,0.5_dp)
        zqss(jl,jk)=zes/(1._dp-vtmpc1*zes)
211  END DO
212 END DO

!!$#ifdef MESSY
!!$  tpot_3d(1:kproma,:,krow) = zteta1(1:kproma,:)  !mz_mt_20031009
!!$  ! mz_ak_20060816+
!!$  IF (lookupoverflow) THEN 
!!$     do jk=1,ktdia,klev
!!$        do jl=1,kproma
!!$           if ( NINT(ptm1(jl,jk)*1000.) <jptlucu1 .OR. &
!!$                NINT(ptm1(jl,jk)*1000.) >jptlucu2)     &
!!$                print*, ' vdiff(1) ptm1jk: ',ptm1(jl,jk)  &
!!$                ,philon_2d(jl,krow),philat_2d(jl,krow)
!!$        enddo
!!$     enddo
!!$  ENDIF
!!$  ! mz_ak_20060816-
!!$#endif
  IF (lookupoverflow) CALL lookuperror ('vdiff (1)   ')

  DO 214 jk=ktdia,klevm1
     DO 213 jl=1,kproma
        zhh(jl,jk)=(pgeom1(jl,jk)-pgeom1(jl,jk+1))
        zsdep1=(paphm1(jl,jk)-paphm1(jl,jk+1))                         &
              /(paphm1(jl,jk)-paphm1(jl,jk+2))
        zsdep2=(paphm1(jl,jk+1)-paphm1(jl,jk+2))                       &
              /(paphm1(jl,jk)  -paphm1(jl,jk+2))
        zqssm(jl,jk)=zsdep1*zqss(jl,jk)+zsdep2*zqss(jl,jk+1)
        ztmitte(jl,jk)=zsdep1*ptm1(jl,jk)+zsdep2*ptm1(jl,jk+1)
        ztvirmit(jl,jk)=zsdep1*ztvir1(jl,jk)+zsdep2*ztvir1(jl,jk+1)
        zfaxen(jl,jk)=zsdep1*zfaxe(jl,jk)+zsdep2*zfaxe(jl,jk+1)
        zlwcmit(jl,jk)=zsdep1*zx(jl,jk)+zsdep2*zx(jl,jk+1)
        zqmit(jl,jk)=zsdep1*pqm1(jl,jk)+zsdep2*pqm1(jl,jk+1)
        ztemit(jl,jk)=zsdep1*zteta1(jl,jk)+zsdep2*zteta1(jl,jk+1)
        zccover(jl,jk)=paclc(jl,jk)*zsdep1+paclc(jl,jk+1)*zsdep2
213  END DO
214 END DO

!
!*      2.2   surface humidity and virtual temperature
!                  for land, water and ice
!
  DO 221 jl=1,kproma
!
!    land ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
     it = NINT(ptslm1(jl)*1000._dp)
     IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
     it = MAX(MIN(it,jptlucu2),jptlucu1)
     zes=tlucua(it)/paphm1(jl,klevp1)
     zqsl(jl)=zes/(1._dp-vtmpc1*zes)
     it1=it+1
     it1= MAX(MIN(it1,jptlucu2),jptlucu1)
     zqst1=tlucua(it1)/paphm1(jl,klevp1)
     zqst1=zqst1/(1._dp-vtmpc1*zqst1)
     zdqsl(jl)=(zqst1-zqsl(jl))*1000._dp
     pws(jl)=MIN(pws(jl),pwsmx(jl))
     zwstop=MIN(0.1_dp,pwsmx(jl))
     zwslev=pwsmx(jl)-zwstop
     IF(pws(jl).GT.zwslev.AND.pws(jl).GT.zwslim*pwsmx(jl)) THEN
        zhum(jl)=0.5_dp*(1._dp-COS((pws(jl)-zwslev)*api/zwstop))
     ELSE
        zhum(jl)=0._dp
     END IF
     zhsoil(jl)=pcvs(jl)+(1._dp-pcvs(jl))                              &
                                  *(pcvw(jl)+(1._dp-pcvw(jl))*zhum(jl))
     lo=pqm1(jl,klev).GT.zqsl(jl)
     zhsoil(jl)=MERGE(1._dp,zhsoil(jl),lo)
     ztesl(jl)=ptslm1(jl)*(1.e5_dp/paphm1(jl,klevp1))**zkappa
     ztvl(jl)=ztesl(jl)*(1._dp+vtmpc1*zhsoil(jl)*zqsl(jl))
!
!    water   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      it = NINT(ptsw(jl)*1000._dp)
      IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
      it = MAX(MIN(it,jptlucu2),jptlucu1)
      zes=tlucua(it)/paphm1(jl,klevp1)
      zqsw(jl)=zes/(1._dp-vtmpc1*zes)
      zcptw(jl)=ptsw(jl)*cpd*(1._dp+vtmpc2*zqsw(jl))
      ztesw(jl)=ptsw(jl)*(1.e5_dp/paphm1(jl,klevp1))**zkappa
      ztvw(jl)=ztesw(jl)*(1._dp+vtmpc1*zqsw(jl))
!
!    ice   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      it = NINT(ptsi(jl)*1000._dp)
      IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
      it = MAX(MIN(it,jptlucu2),jptlucu1)
      zes=tlucua(it)/paphm1(jl,klevp1)
      zqsi(jl)=zes/(1._dp-vtmpc1*zes)
      zcpti(jl)=ptsi(jl)*cpd*(1._dp+vtmpc2*zqsi(jl))
      ztesi(jl)=ptsi(jl)*(1.e5_dp/paphm1(jl,klevp1))**zkappa
      ztvi(jl)=ztesi(jl)*(1._dp+vtmpc1*zqsi(jl))
!
221 END DO

!!$#ifdef MESSY
!!$  ! mz_ak_20060816+
!!$  IF (lookupoverflow) THEN 
!!$     do jl=1,kproma
!!$        if ( NINT(ptslm1(jl)*1000.) <jptlucu1 .OR. &
!!$             NINT(ptslm1(jl)*1000.) >jptlucu2)     &
!!$             print*, 'vdiff(2) ptslm1: ',ptslm1(jl) &
!!$                ,philon_2d(jl,krow),philat_2d(jl,krow)
!!$        if ( NINT(ptsw(jl)*1000.) <jptlucu1 .OR. &
!!$             NINT(ptsw(jl)*1000.) >jptlucu2)     &
!!$             print*, 'vdiff(2) ptsw: ',ptsw(jl) &
!!$                ,philon_2d(jl,krow),philat_2d(jl,krow)
!!$        if ( NINT(ptsi(jl)*1000.) <jptlucu1 .OR. &
!!$             NINT(ptsi(jl)*1000.) >jptlucu2)     &
!!$             print*, 'vdiff(2) ptsi: ',ptsi(jl) &
!!$                ,philon_2d(jl,krow),philat_2d(jl,krow)
!!$     enddo
!!$  ENDIF
!!$  ! mz_ak_20060816-
!!$#endif 
    IF (lookupoverflow) CALL lookuperror ('vdiff (2)   ')

!
!     surface net solar radiation flux over land
!
     DO jl=1,kproma
        zsrfld=psrfl(jl)/(1._dp-palbedo(jl))
        zsrfll(jl)=zsrfld*(1._dp-palsol(jl))
     END DO
!
!*         2.3     DEFINITION OF THE STOMATAL RESISTANCE
!
  DO 231 jl=1,kproma
     zwcrit=zplmax*pwsmx(jl)
     zwpwp=zplmin*pwsmx(jl)
     zqwevap=1._dp/(zwcrit-zwpwp)
     zsoil=MAX(0._dp,MIN(1._dp,(pws(jl)-zwpwp)*zqwevap))
     zsrfl=MAX(zepsr,zsrfll(jl)*zvrad)
     zabcs=(zva+zvbc)/(zvc*zsrfl)
     zlai=pvlt(jl)
     zvklt=zvk*zlai
     zvxpklt=EXP(zvklt)
     zvxmklt=EXP(-zvklt)
     zln1=LOG((zabcs*zvxpklt+1._dp)/(zabcs+1._dp))
     zln2=LOG((zabcs+zvxmklt)/(zabcs+1._dp))
     zrsi=(zvb*zln1/zvabc-zln2)/zvkc
     zwet(jl)=1._dp/(zrsi*zsoil+zepevap)
     lo=pqm1(jl,klev).GT.zqsl(jl)
     zwet(jl)=MERGE(0._dp,zwet(jl),lo)
!!$#ifdef MESSY
!!$     ! mz_lg_20020129+
!!$     fws(jl,krow) = zsoil
!!$     ! calculation of leaf stomatal resistance using the echam
!!$     ! equation applying an LAI of 1, included by Laurens Ganzeveld, 18/10/01 
!!$     rco_leaf(jl,krow)=1./((zvb*(LOG((zabcs*EXP(zvk)+1.)/(zabcs+1.))) &
!!$               /zvabc-(LOG((zabcs+EXP(-zvk))/(zabcs+1.))))/zvkc)
!!$     ! mz_lg_20020129-
!!$     zwet_tmp(jl)=zwet(jl) !store temporary result for H2OISO op_re_20140904
!!$#endif
231 END DO
!
  IF (lvdiff) THEN
!
!     ------------------------------------------------------------------
!
!*         3.     COMPUTATION OF THE EXCHANGE COEFFICIENTS.
!
!        THE SURFACE LAYER IS NOW COMPUTED BEFORE THE OTHER LEVELS
!
!        3.1       COMPUTATION OF BASIC QUANTITIES: WIND SHEAR,
!                  RICHARDSON NUMBER,SQUARED MIXING LENGTHS, UNSTABLE
!                  AND STABLE CASE COMMON FACTORS AND NEUTRAL CASE
!                  COMMON PART OF THE DRAG COEFFICIENTS.
!
     DO 311 jl=1,kproma
!
!     land   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      zdu2(jl)=MAX(zepdu2,pum1(jl,klev)**2+pvm1(jl,klev)**2)
      zqmitte=(pqm1(jl,klev)+zqsl(jl)*zhsoil(jl))/2._dp
      zqtmit=zx(jl,klev)*0.5_dp+zqmitte
      ztmit=(ptm1(jl,klev)+ptslm1(jl))/2._dp
      zqsmit=(zqss(jl,klev)+zqsl(jl))/2._dp
      ztemitte=(zteta1(jl,klev)+ztesl(jl))/2._dp
      zvirmitte=(ztvir1(jl,klev)+ztvl(jl))/2._dp
      zqlwi1=pqm1(jl,klev)+zx(jl,klev)
      zqlwi2=zqsl(jl)*zhsoil(jl)
      zfux=zfaxe(jl,klev)/(zcpd*ztmit)
      zfox=zfaxe(jl,klev)/(zrd*ztmit)
      zmult1=1._dp+vtmpc1*zqtmit
      zmult2=zfux*zmult1-zrvrd
      zmult3=zrdrv*zfox*zqsmit/(1._dp+zrdrv*zfox*zfux*zqsmit)
      zmult5=zmult1-zmult2*zmult3
      zmult4=zfux*zmult5-1._dp
      zdus1=paclc(jl,klev)*zmult5+(1._dp-paclc(jl,klev))*zmult1
      zdus2=paclc(jl,klev)*zmult4+(1._dp-paclc(jl,klev))*vtmpc1
      zteldif=zlteta1(jl,klev)-ztesl(jl)
      zqddif=zqlwi1-zqlwi2
      zbuoy=zdus1*zteldif+zdus2*ztemitte*zqddif
      zril(jl)=pgeom1(jl,klev)*zbuoy/(zvirmitte*zdu2(jl))
      z0h=MIN(1._dp,paz0l(jl))
      IF(pcvs(jl).GT.0._dp) THEN
        zchsnow=(LOG(zblend/paz0i(jl)))**2
        zchland=(LOG(zblend/z0h))**2
        zchsnow=pcvs(jl)/zchsnow
        zchland=(1._dp-pcvs(jl))/zchland
        zchmean=1._dp/SQRT(zchsnow+zchland)
        z0h=zblend*EXP(-zchmean)
      END IF
      zcdnl(jl)=(zkap/LOG(1._dp+pgeom1(jl,klev)/(g*paz0l(jl))))**2
      zchnl(jl)=(zkap/LOG(1._dp+pgeom1(jl,klev)/(g*z0h)))**2
      zucfl(jl)=1._dp/(1._dp+zcons11*zcdnl(jl)*SQRT(ABS(zril(jl))      &
                             *(1._dp +pgeom1(jl,klev)/(g*paz0l(jl)))))
      zucfhl(jl)=1._dp/(1._dp+zcons11*zchnl(jl)*SQRT(ABS(zril(jl))     &
                             *(1._dp +pgeom1(jl,klev)/(g*z0h))))
      zscfl(jl)=SQRT(1._dp+ABS(zril(jl)))
      zcons=zcons12*paphm1(jl,klevp1)/                                 &
            (ptm1(jl,klev)*(1._dp+vtmpc1*pqm1(jl,klev)-zx(jl,klev)))
      zcfncl(jl)=zcons*SQRT(zdu2(jl))*zcdnl(jl)
      zcfnchl(jl)=zcons*SQRT(zdu2(jl))*zchnl(jl)
      zdthv=MAX(0._dp,(ztvl(jl)-ztvir1(jl,klev)))
      zwstl(jl)=zdthv*SQRT(zdu2(jl))/zvirmitte
!
!    water   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! correction for water and ice points
!
      zdu2oc(jl)=MAX(zepdu2,(pum1(jl,klev)-pocu(jl))**2                &
                     +(pvm1(jl,klev)-pocv(jl))**2)
!
      zqmitte=(pqm1(jl,klev)+zqsw(jl))/2._dp
      zqtmit=zx(jl,klev)*0.5_dp+zqmitte
      ztmit=(ptm1(jl,klev)+ptsw(jl))/2._dp
      zqsmit=(zqss(jl,klev)+zqsw(jl))/2._dp
      ztemitte=(zteta1(jl,klev)+ztesw(jl))/2._dp
      zvirmitte=(ztvir1(jl,klev)+ztvw(jl))/2._dp
      zqlwi1=pqm1(jl,klev)+zx(jl,klev)
      zqlwi2=zqsw(jl)
      zfux=zfaxe(jl,klev)/(zcpd*ztmit)
      zfox=zfaxe(jl,klev)/(zrd*ztmit)
      zmult1=1._dp+vtmpc1*zqtmit
      zmult2=zfux*zmult1-zrvrd
      zmult3=zrdrv*zfox*zqsmit/(1._dp+zrdrv*zfox*zfux*zqsmit)
      zmult5=zmult1-zmult2*zmult3
      zmult4=zfux*zmult5-1._dp
      zdus1=paclc(jl,klev)*zmult5+(1._dp-paclc(jl,klev))*zmult1
      zdus2=paclc(jl,klev)*zmult4+(1._dp-paclc(jl,klev))*vtmpc1
      zteldif=zlteta1(jl,klev)-ztesw(jl)
      zqddif=zqlwi1-zqlwi2
      zbuoy=zdus1*zteldif+zdus2*ztemitte*zqddif
      zriw(jl)=pgeom1(jl,klev)*zbuoy/(zvirmitte*zdu2oc(jl))
      z0h=paz0w(jl)*EXP(2._dp-86.276_dp*paz0w(jl)**0.375_dp)
      zalo=LOG(1._dp+pgeom1(jl,klev)/(g*paz0w(jl)))
      zaloh=LOG(1._dp+pgeom1(jl,klev)/(g*z0h))
      zcdnw(jl)=(zkap/zalo)**2
      zchnw(jl)=zkap**2/(zalo*zaloh)
      zucfw(jl)=1._dp/(1._dp+zcons11*zcdnw(jl)*SQRT(ABS(zriw(jl))      &
                              *(1._dp+pgeom1(jl,klev)/(g*paz0w(jl)))))
      zscfw(jl)=SQRT(1._dp+ABS(zriw(jl)))
      zcfncw(jl)=zcons*SQRT(zdu2oc(jl))*zcdnw(jl)
      zcfnchw(jl)=zcons*SQRT(zdu2oc(jl))*zchnw(jl)
      zdthv=MAX(0._dp,(ztvw(jl)-ztvir1(jl,klev)))
      zwstw(jl)=zdthv*SQRT(zdu2oc(jl))/zvirmitte
      zcr(jl)=(zfreec/(zchnw(jl)*SQRT(zdu2oc(jl))))*ABS(zbuoy)**zcons6
!
!     ice   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      zqmitte=(pqm1(jl,klev)+zqsi(jl))/2._dp
      zqtmit=zx(jl,klev)*0.5_dp+zqmitte
      ztmit=(ptm1(jl,klev)+ptsi(jl))/2._dp
      zqsmit=(zqss(jl,klev)+zqsi(jl))/2._dp
      ztemitte=(zteta1(jl,klev)+ztesi(jl))/2._dp
      zvirmitte=(ztvir1(jl,klev)+ztvi(jl))/2._dp
      zqlwi1=pqm1(jl,klev)+zx(jl,klev)
      zqlwi2=zqsi(jl)
      zfux=zfaxe(jl,klev)/(zcpd*ztmit)
      zfox=zfaxe(jl,klev)/(zrd*ztmit)
      zmult1=1._dp+vtmpc1*zqtmit
      zmult2=zfux*zmult1-zrvrd
      zmult3=zrdrv*zfox*zqsmit/(1._dp+zrdrv*zfox*zfux*zqsmit)
      zmult5=zmult1-zmult2*zmult3
      zmult4=zfux*zmult5-1._dp
      zdus1=paclc(jl,klev)*zmult5+(1._dp-paclc(jl,klev))*zmult1
      zdus2=paclc(jl,klev)*zmult4+(1._dp-paclc(jl,klev))*vtmpc1
      zteldif=zlteta1(jl,klev)-ztesi(jl)
      zqddif=zqlwi1-zqlwi2
      zbuoy=zdus1*zteldif+zdus2*ztemitte*zqddif
      zrii(jl)=pgeom1(jl,klev)*zbuoy/(zvirmitte*zdu2oc(jl))
      zalo=LOG(1._dp+pgeom1(jl,klev)/(g*paz0i(jl)))
      zcdni(jl)=(zkap/zalo)**2
      zucfi(jl)=1._dp/(1._dp+zcons11*zcdni(jl)*SQRT(ABS(zrii(jl))      &
                  *(1._dp+pgeom1(jl,klev)/(g*paz0i(jl)))))
      zscfi(jl)=SQRT(1._dp+ABS(zrii(jl)))
      zcfnci(jl)=zcons*SQRT(zdu2oc(jl))*zcdni(jl)
      zdthv=MAX(0._dp,(ztvi(jl)-ztvir1(jl,klev)))
      zwsti(jl)=zdthv*SQRT(zdu2oc(jl))/zvirmitte
311  END DO

!!$#ifdef MESSY
!!$    ! mz_rs_20020910+
!!$    tvir(1:kproma,krow) = ztvir1(1:kproma,klev)
!!$    ! mz_rs_20020910-
!!$#endif

!
!     3.2  DIMENSIONLESS HEAT TRANSFER COEFFICIENTS MULTIPLIED
!          BY PRESSURE THICKNESSES FOR MOMENTUM AND HEAT EXCHANGE
!
     DO 321 jl=1,kproma
!
!     land   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      IF(zril(jl).GT.0._dp) THEN
         zcfml(jl)=zcfncl(jl)/(1._dp+zcons8*zril(jl)/zscfl(jl))
         zcfhl(jl)=zcfnchl(jl)/(1._dp+zcons8*zril(jl)*zscfl(jl))
         zchl(jl)=zcfhl(jl)/zcfnchl(jl)*zchnl(jl)
      ELSE
         zcfml(jl)=zcfncl(jl)*(1._dp-zcons8*zril(jl)*zucfl(jl))
         zcfhl(jl)=zcfnchl(jl)*(1._dp-zcons9*zril(jl)*zucfhl(jl))
         zchl(jl)=zcfhl(jl)/zcfnchl(jl)*zchnl(jl)
      END IF
!
!     water   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      IF(zriw(jl).GT.0._dp) THEN
         zcfmw(jl)=zcfncw(jl)/(1._dp+zcons8*zriw(jl)/zscfw(jl))
         zcfhw(jl)=zcfnchw(jl)/(1._dp+zcons8*zriw(jl)*zscfw(jl))
         zchw(jl)=zcfhw(jl)/zcfnchw(jl)*zchnw(jl)
      ELSE
         zcfmw(jl)=zcfncw(jl)*(1._dp-zcons8*zriw(jl)*zucfw(jl))
         zcfhw(jl)=zcfnchw(jl)*(1._dp+zcr(jl)**zgam)**z1dgam
         zchw(jl)=zcfhw(jl)/zcfnchw(jl)*zchnw(jl)
      END IF
!
!     ice   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      IF(zrii(jl).GT.0._dp) THEN
         zcfmi(jl)=zcfnci(jl)/(1._dp+zcons8*zrii(jl)/zscfi(jl))
         zcfhi(jl)=zcfnci(jl)/(1._dp+zcons8*zrii(jl)*zscfi(jl))
         zchi(jl)=zcfhi(jl)/zcfnci(jl)*zcdni(jl)
      ELSE
         zcfmi(jl)=zcfnci(jl)*(1._dp-zcons8*zrii(jl)*zucfi(jl))
         zcfhi(jl)=zcfnci(jl)*(1._dp-zcons9*zrii(jl)*zucfi(jl))
         zchi(jl)=zcfhi(jl)/zcfnci(jl)*zcdni(jl)
      END IF
!
!     aggregated exchange coefficient for momentum
!
      zcfm(jl,klev)=pfrl(jl)*zcfml(jl)+pfrw(jl)*zcfmw(jl)              &
                   +pfri(jl)*zcfmi(jl)
      zcdum(jl,klev)=zcfm(jl,klev)
!
!     interpolation functions for diagnostics
!
!     land   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      zbnl(jl)=zkap/SQRT(zcdnl(jl))
      zbhnl(jl)=zkap/SQRT(zchnl(jl))
      zbml(jl)=MAX(zepsec,SQRT(zcfml(jl)*zcdnl(jl)*zcons17/zcfncl(jl)))
      zbhl(jl)=MAX(zepsec,zchl(jl)/zbml(jl)*zcons17)
      zbml(jl)=1._dp/zbml(jl)
      zbhl(jl)=1._dp/zbhl(jl)
!
!     water   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      zbnw(jl)=zkap/SQRT(zcdnw(jl))
      zbmw(jl)=MAX(zepsec,SQRT(zcfmw(jl)*zcdnw(jl)*zcons17/zcfncw(jl)))
      zbhw(jl)=MAX(zepsec,zchw(jl)/zbmw(jl)*zcons17)
      zbmw(jl)=1._dp/zbmw(jl)
      zbhw(jl)=1._dp/zbhw(jl)
!
!     ice   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      zbni(jl)=zkap/SQRT(zcdni(jl))
      zbmi(jl)=MAX(zepsec,SQRT(zcfmi(jl)*zcdni(jl)*zcons17/zcfnci(jl)))
      zbhi(jl)=MAX(zepsec,zchi(jl)/zbmi(jl)*zcons17)
      zbmi(jl)=1._dp/zbmi(jl)
      zbhi(jl)=1._dp/zbhi(jl)
!!$#ifdef MESSY
!!$      !mz_mt_20031120+
!!$      dens= paphm1(jl,klevp1)/zrd         &
!!$           /ptm1(jl,klev)*(1.+vtmpc1*pqm1(jl,klev))  ! density [kg m-3]
!!$      zcdh_2d(jl,krow)=zcfhl(jl)/(zcons12*SQRT(zdu2(jl))*zrd*dens) 
!!$      !mz_mt_20031120-
!!$#endif
321  END DO
!
!    INITIALIZE SURFACE EMISSION FOR TRACERS
!
!!$#ifndef MESSY
     IF(ktrac.GT.0) THEN
!!$#endif
        DO 3231 jt=1,ktrac
           DO 3230 jl=1,kproma
              zxtems(jl,jt)=0._dp
3230       END DO
3231    END DO
!
!
!     SURFACE EMISSIONS AND DRY DEPOSITION
!
!!        IF(lemis) THEN
!!           DO 324 jl=1,kproma
!!              z1mxtm1(jl)=papm1(jl,klev)                               &
!!                   /(ptm1(jl,klev)*rd*(1._dp+vtmpc1*pqm1(jl,klev)))
!!324        END DO
!!!
!!           CALL xtemiss (kproma,   klev,     irow,   cvdifts,  zdtime, &
!!                         pxtm1,  zxtems,   z1mxtm1,                    &
!!                         lpland, pforest,  psn)
!!        END IF

!!$#ifndef MESSY
     CALL call_chem_bcond(kproma, kbdim, klev, pxtte)
!!$#else
!!$      ! mz_lg_20020115+ 
!!$      ! The new parameter parameter zxtems contains the emission
!!$      ! flux, corrected for dry deposition resulting in the definition of the
!!$      ! net surface exchange flux. zxtems is used within vdiff in the
!!$      ! loop of vertical transport calculations as the lower boundary
!!$      ! condition.
!!$      ! mz_lg_20020115- 
!!$
!!$         !mz_pj_20040326+
!!$! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$         CALL messy_vdiff
!!$! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$         !mz_pj_20040326-
!!$#endif

!!$#ifndef MESSY
     END IF
!!$#endif
!
!     3.3 equivalent evapotranspiration efficiency coeff. over land
!
     DO 330 jl=1,kproma
        IF(lpland(jl).AND..NOT.lpglac(jl)) THEN
          zepot=zcons30*zcfhl(jl)*(zqsl(jl)-pqm1(jl,klev))
          IF(pcvs(jl).GT.0._dp) THEN
            zsnfac=pcvs(jl)*zepot*zdtime/(rhoh2o*(psn(jl)+psnc(jl)))
            IF(zsnfac.GT.1._dp) pcvs(jl)=pcvs(jl)/zsnfac
          END IF
          IF(pcvw(jl).GT.0._dp) THEN
            zwlfac=(1._dp-pcvs(jl))*zepot*zdtime/(rhoh2o*pwlmx(jl))
            IF(zwlfac.GT.1._dp) pcvw(jl)=pcvw(jl)/zwlfac
          END IF
        END IF
 330 CONTINUE
     DO 331 jl=1,kproma
        IF (pws(jl).GT.zwslim*pwsmx(jl)) THEN
           zwet(jl)=pcvs(jl)+(1._dp-pcvs(jl))*(pcvw(jl)                &
                            +(1._dp-pcvw(jl))/                         &
                             (1._dp+zchl(jl)*SQRT(zdu2(jl))*zwet(jl)))
        ELSE
           zwet(jl)=pcvs(jl)+(1._dp-pcvs(jl))*pcvw(jl)
        END IF
        lo=zhum(jl).LE.pqm1(jl,klev)/zqsl(jl) .OR. zhum(jl).LT.zepevap
        zcsat(jl)=pcvs(jl)+(1._dp-pcvs(jl))                            &
                *(pcvw(jl)+(1._dp-pcvw(jl))*MERGE(0._dp,zhum(jl),lo))
        zcair(jl)=pcvs(jl)+(1._dp-pcvs(jl))                            &
                *(pcvw(jl)+(1._dp-pcvw(jl))*MERGE(0._dp,1._dp,lo))
        lo=pqm1(jl,klev).GT.zqsl(jl)
        zcsat(jl)=MERGE(1._dp,zcsat(jl),lo)
        zcair(jl)=MERGE(1._dp,zcair(jl),lo)
        zcsat(jl)=pvgrat(jl)*zwet(jl)+(1._dp-pvgrat(jl))*zcsat(jl)
        zcair(jl)=pvgrat(jl)*zwet(jl)+(1._dp-pvgrat(jl))*zcair(jl)
        zcptl(jl)=ptslm1(jl)*cpd*(1._dp+vtmpc2*(zcsat(jl)*zqsl(jl)     &
                                +(1._dp-zcair(jl))*pqm1(jl,klev)))
        zcpq(jl)=zcptl(jl)/ptslm1(jl)
331  END DO
!
!*       3.4       COMPUTATION OF THE PBL EXTENSION.
!
     DO 341 jl=1,kproma
!
!    land   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      lo=paz0l(jl).GT.zepz0o
      zcdn2m=MERGE((zkap/LOG(1._dp+pgeom1(jl,klev)/(g*zepz0o)))**2,    &
                      zcdnl(jl),lo)
      zcdnr=zcdn2m/zcdnl(jl)
      zcfm2m=MERGE(zcfncl(jl)*zcdnr*(1._dp-zcons8*zril(jl)             &
                    /(1._dp+zcons11*zcdn2m*SQRT(ABS(zril(jl))          &
                    *(1._dp+pgeom1(jl,klev)/(g*zepz0o))))),            &
                     zcfml(jl)*zcdnr,lo.AND.zril(jl).LT.0._dp)
      zustl=zcfm2m*SQRT(zdu2(jl))
      zustarl(jl)=SQRT(zustl*ptm1(jl,klev)                             &
                     *(1._dp+vtmpc1*pqm1(jl,klev)-zx(jl,klev))         &
                     /(zcons12*paphm1(jl,klevp1)))
!
!     water   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      lo=paz0w(jl).GT.zepz0o
      zcdn2m=MERGE((zkap/LOG(1._dp+pgeom1(jl,klev)/(g*zepz0o)))**2,    &
                    zcdnw(jl),lo)
      zcdnr=zcdn2m/zcdnw(jl)
      zcfm2m=MERGE(zcfncw(jl)*zcdnr*(1._dp-zcons8*zriw(jl))            &
                  /(1._dp+zcons11*zcdn2m*SQRT(ABS(zriw(jl))            &
                  *(1._dp+pgeom1(jl,klev)/(g*zepz0o)))),               &
                   zcfmw(jl)*zcdnr,lo.AND.zriw(jl).LT.0._dp)
      zustw=zcfm2m*SQRT(zdu2oc(jl))
      zustarw(jl)=SQRT(zustw*ptm1(jl,klev)                             &
                      *(1._dp+vtmpc1*pqm1(jl,klev)-zx(jl,klev))        &
                      /(zcons12*paphm1(jl,klevp1)))
!
!    ice   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      lo=paz0i(jl).GT.zepz0o
      zcdn2m=MERGE((zkap/LOG(1._dp+pgeom1(jl,klev)/(g*zepz0o)))**2,    &
                    zcdni(jl),lo)
      zcdnr=zcdn2m/zcdni(jl)
      zcfm2m=MERGE(zcfnci(jl)*zcdnr*(1._dp-zcons8*zrii(jl))            &
                  /(1._dp+zcons11*zcdn2m*SQRT(ABS(zrii(jl))            &
                  *(1._dp+pgeom1(jl,klev)/(g*zepz0o)))),               &
                   zcfmi(jl)*zcdnr,lo.AND.zrii(jl).LT.0._dp)
      zusti=zcfm2m*SQRT(zdu2oc(jl))
      zustari(jl)=SQRT(zusti*ptm1(jl,klev)                             &
                      *(1._dp+vtmpc1*pqm1(jl,klev)-zx(jl,klev))        &
                      /(zcons12*paphm1(jl,klevp1)))
!
      zust=pfrl(jl)*zustl+pfrw(jl)*zustw+pfri(jl)*zusti
!!$#ifdef MESSY
!!$      ! mz_mt_20031119+
!!$      zust_2d(jl,krow) =zustarl(jl)*pfrl(jl) &
!!$                       +zustarw(jl)*pfrw(jl) &
!!$                       +zustari(jl)*pfri(jl) 
!!$      ! mz_mt_20031119-
!!$#endif
!
      zustarm=SQRT(zust*ptm1(jl,klev)*                                 &
                  (1._dp+vtmpc1*pqm1(jl,klev)-zx(jl,klev))             &
                  /(zcons12*paphm1(jl,klevp1)))
      zcor=MAX(ABS(coriol_2d(jl,krow)),zepcor)
      zhdyn(jl)=MIN(pgeom1(jl,1)/g,zchneu*zustarm/zcor)
!
     ihpblc(jl)=klev
     ihpbld(jl)=klev
341  END DO
!
     DO 343 jk=klevm1,1,-1
        DO 342 jl=1,kproma
           zds=zcptgz(jl,jk)-zcptgz(jl,klev)
           zdz=pgeom1(jl,jk)/g-zhdyn(jl)
           ihpblc(jl)=MERGE(jk,ihpblc(jl),                             &
                            ihpblc(jl).EQ.klev.AND.zds.GT.0._dp)
           ihpbld(jl)=MERGE(jk,ihpbld(jl),                             &
                            ihpbld(jl).EQ.klev.AND.zdz.GE.0._dp)
342     END DO
343  END DO
!
!      CONVECTIVE VELOCITY SCALE, MONIN-OBUKHOV LENGTH AND
!      TKE BOUNDARY CONDITION (MAILHOT/BENOIT, 1982)
!
     DO 344 jl=1,kproma
        ihpbl(jl)=MIN(ihpblc(jl),ihpbld(jl))
        zghabl=MIN(50000._dp,pgeom1(jl,ihpbl(jl)))
!
!       land   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        IF(zwstl(jl).GT.zepsr) THEN
           zconvs=(zwstl(jl)*zchl(jl)*zghabl)**zcons6
           zmonob=(zustarl(jl)**3)/(zkap*g*zwstl(jl)*zchl(jl))
           zstabf=(pgeom1(jl,klev)/(g*zmonob))**(zcons6*2._dp)
           zstabf=MIN(zustf*3._dp,zstabf)
        ELSE
           zconvs=0._dp
           zstabf=0._dp
        END IF
        ztkevl=(zustf+zstabf)*(zustarl(jl)**2)+zwstf*(zconvs**2)
!
!       water   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        IF(zwstw(jl).GT.zepsr) THEN
           zconvs=(zwstw(jl)*zchw(jl)*zghabl)**zcons6
           zmonob=(zustarw(jl)**3)/(zkap*g*zwstw(jl)*zchw(jl))
           zstabf=(pgeom1(jl,klev)/(g*zmonob))**(zcons6*2._dp)
           zstabf=MIN(zustf*3._dp,zstabf)
        ELSE
           zconvs=0._dp
           zstabf=0._dp
        END IF
        ztkevw=(zustf+zstabf)*(zustarw(jl)**2)+zwstf*(zconvs**2)
!
!       ice   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        IF(zwsti(jl).GT.zepsr) THEN
           zconvs=(zwsti(jl)*zchi(jl)*zghabl)**zcons6
           zmonob=(zustari(jl)**3)/(zkap*g*zwsti(jl)*zchi(jl))
           zstabf=(pgeom1(jl,klev)/(g*zmonob))**(zcons6*2._dp)
           zstabf=MIN(zustf*3._dp,zstabf)
        ELSE
           zconvs=0._dp
           zstabf=0._dp
        END IF
        ztkevi=(zustf+zstabf)*(zustari(jl)**2)+zwstf*(zconvs**2)
        ztkevn(jl,klev)=pfrl(jl)*ztkevl+pfrw(jl)*ztkevw+pfri(jl)*ztkevi
        ztkevn(jl,klev)=MAX(ztkemin,ztkevn(jl,klev))
344  END DO
!
     IF(lstart) THEN
        DO 345 jl=1,kproma
           ptkem1(jl,klev)=ztkevn(jl,klev)
           ptkem(jl,klev)=ztkevn(jl,klev)
345     END DO
!!$#ifdef MESSY
!!$        pahfsm1(:)  = 0.0 ! mz_mt_20031120
!!$#endif
     END IF
!
!     ==================================================================
!
!*       3.5   Vertical loop: Computation of basic quantities:
!              wind shear, buoyancy, Ri-number, mixing length
!
     DO 372 jk=ktdia,klevm1
        DO 361 jl=1,kproma
           zqtmit=zlwcmit(jl,jk)+zqmit(jl,jk)
           zfux=zfaxen(jl,jk)/(zcpd*ztmitte(jl,jk))
           zfox=zfaxen(jl,jk)/(zrd*ztmitte(jl,jk))
           zmult1=1._dp+vtmpc1*zqtmit
           zmult2=zfux*zmult1-zrvrd
           zmult3=zrdrv*zfox*zqssm(jl,jk)                              &
                  /(1._dp+zrdrv*zfux*zfox*zqssm(jl,jk))
           zmult5=zmult1-zmult2*zmult3
           zmult4=zfux*zmult5-1._dp
           zdus1=zccover(jl,jk)*zmult5+(1._dp-zccover(jl,jk))*zmult1
           zdus2=zccover(jl,jk)*zmult4+(1._dp-zccover(jl,jk))*vtmpc1
           zteldif=(zlteta1(jl,jk)-zlteta1(jl,jk+1))/zhh(jl,jk)*g
           zdqtot=(pqm1(jl,jk)+zx(jl,jk))-(pqm1(jl,jk+1)+zx(jl,jk+1))
           zqddif=zdqtot/zhh(jl,jk)*g
           zqshear(jl,jk)=zqddif !store for variance production
           zbuoy=(zteldif*zdus1+ztemit(jl,jk)*zdus2*zqddif)            &
                 *g/ztvirmit(jl,jk)
           zdivv=(pum1(jl,jk)-pum1(jl,jk+1))**2
           zdivv1=(pvm1(jl,jk)-pvm1(jl,jk+1))**2
           zshear=(zdivv+zdivv1)*(g/zhh(jl,jk))**2
           zri=zbuoy/MAX(zshear,zepshr)
!!$#ifdef MESSY
!!$           rinum_3d(jl,jk,krow) = zri
!!$#endif
!
!      ASYMPTOTIC MIXING LENGTH FOR MOMENTUM AND
!      HEAT (ZLAM) ABOVE THE PBL AS A FUNCTION OF HEIGHT
!      ACCORDING TO HOLTSLAG AND BOVILLE (1992), J. CLIMATE.
!
           zhexp=EXP(1._dp-pgeom1(jl,jk)/pgeom1(jl,ihpbl(jl)))
           zlam=zzzlam+(zcons3-zzzlam)*zhexp
           IF(jk.GE.ihpbl(jl)) THEN
              zcons23=zcons25
           ELSE
              zcons23=zcons2/zlam
           END IF
!
!     MIXING LENGTH (BLACKADAR) + STABILITY DEPENDENT FUNCTION
!
           z2geomf=pgeom1(jl,jk)+pgeom1(jl,jk+1)
           zz2geo=zcons2*z2geomf
           zmix=zz2geo/(1._dp+zcons23*z2geomf)
!
!      STABILITY FUNCTIONS (LOUIS, 1979)
!
           IF(zri.LT.0._dp) THEN
              zalh2=zmix*zmix
              zucf=1._dp/                                              &
                   (1._dp+zcons5*zalh2*SQRT(ABS(zri)*(((pgeom1(jl,jk)  &
                       /pgeom1(jl,jk+1))**zcons6-1._dp)/(pgeom1(jl,jk) &
                       -pgeom1(jl,jk+1)))**3/(pgeom1(jl,jk+1))))
              zsh=zshn*(1._dp-zcons9*zri*zucf)*zmix
              zsm=zsmn*(1._dp-zcons8*zri*zucf)*zmix
           ELSE
              zsh=zshn/(1._dp+zcons8*zri*SQRT(1._dp+zri))*zmix
              zsm=zsmn/(1._dp+zcons8*zri/SQRT(1._dp+zri))*zmix
           END IF
!
!       Dimensionless coefficients multiplied by pressure
!            thicknesses for momentum and heat exchange
!
           zzb=zshear*zsm-zbuoy*zsh
           zdisl=zda1*zmix/ztmst
           zktest=1._dp+(zzb*ztmst+SQRT(ptkem1(jl,jk))*2._dp)/zdisl
           IF (zktest.LE.1._dp) THEN
              ztkevn(jl,jk)=ztkemin
           ELSE
              ztkevn(jl,jk)=MAX(ztkemin,(zdisl*(SQRT(zktest)-1._dp))**2)
           END IF
           IF(lstart) THEN
              ptkem1(jl,jk)=ztkevn(jl,jk)
              ptkem(jl,jk)=ztkevn(jl,jk)
           END IF
           ztkesq=SQRT(MAX(ztkemin,ptkem1(jl,jk)))
           zztvm=(ptvm1(jl,jk)+ptvm1(jl,jk+1))*0.5_dp
           zalf=paphm1(jl,jk+1)/(zztvm*zhh(jl,jk)*zrd)
           zcfm(jl,jk)=zsm*ztkesq*zcons18*zalf
           zcfh(jl,jk)=zsh*ztkesq*zcons18*zalf
           zcfv(jl,jk)=0.5_dp*zcfh(jl,jk)
           zcdum(jl,jk)=zcfm(jl,jk)/ztkesq*SQRT(ztkevn(jl,jk))
361     END DO
372  END DO
!
!     ==================================================================
!
!*       3.8        DIFFUSION IMPLICIT COMPUTATIONS FOR TKE
!
     DO 380 jk=ktdia,klev
        DO 381 jl=1,kproma
           zedif(jl,jk)=ztpfac2*ztkevn(jl,jk)
381     END DO
380  END DO
!
     DO 385 jl=1,kproma
        ztcoe(jl)=(zcdum(jl,itop)+zcdum(jl,itopp1))*0.5_dp
        zqdp=1._dp/(papm1(jl,itopp1)-papm1(jl,itop))
        zdisc=1._dp/(1._dp+(zcdum(jl,itop)+zcdum(jl,itopp1))           &
                                                         *0.5_dp*zqdp)
        zebsm(jl,itop)=zdisc*(zcdum(jl,itop)+zcdum(jl,itopp1))         &
                                                         *0.5_dp*zqdp
        zedif(jl,itop)=zdisc*zedif(jl,itop)
385  END DO
!
     DO 386 jk=itopp1,klev-2
        DO 387 jl=1,kproma
           zqdp=1._dp/(papm1(jl,jk+1)-papm1(jl,jk))
           zfac=ztcoe(jl)*zqdp
           ztcoe(jl)=(zcdum(jl,jk+1)+zcdum(jl,jk))*0.5_dp
           zdisc=1._dp/(1._dp+zfac*(1.-zebsm(jl,jk-1))                 &
                     +(zcdum(jl,jk+1)+zcdum(jl,jk))*0.5_dp*zqdp)
           zebsm(jl,jk)=zdisc*(zcdum(jl,jk+1)+zcdum(jl,jk))*0.5_dp*zqdp
           zedif(jl,jk)=zdisc*(zedif(jl,jk)+zfac*zedif(jl,jk-1))
387     END DO
386  END DO
!
     DO 390 jl=1,kproma
        zqdp=1._dp/(papm1(jl,klev)-papm1(jl,klevm1))
        zfac=ztcoe(jl)*zqdp
        ztcoe(jl)=(zcdum(jl,klev)+zcdum(jl,klevm1))*0.5_dp
        zdisc=1._dp/(1._dp+zfac*(1._dp-zebsm(jl,klev-2))               &
                  +(zcdum(jl,klev)+zcdum(jl,klevm1))*0.5_dp*zqdp)
        zedif(jl,klevm1)=zdisc*((zcdum(jl,klev)+zcdum(jl,klevm1))      &
                         *0.5_dp*zqdp*zedif(jl,klev)+zedif(jl,klevm1)  &
                         +zfac*zedif(jl,klev-2))
390  END DO
!
     DO 392 jk=klev-2,itop,-1
        DO 393 jl=1,kproma
           zedif(jl,jk)=zedif(jl,jk)+zebsm(jl,jk)*zedif(jl,jk+1)
393     END DO
392  END DO
!
!*    TIME INTEGRATION OF TURBULENT KINETIC ENERGY AND CHECK
!
     DO 394 jk=itop,klev
        ztest=0._dp
        DO 395 jl=1,kproma
           ptke(jl,jk)=zedif(jl,jk)+ztpfac3*ztkevn(jl,jk)
           ztest=ztest+MERGE(1._dp,0._dp,ptke(jl,jk)<0._dp)
395     END DO
        IF(ztest.NE.0._dp) CALL finish('vdiff','TKE IS NEGATIVE')
394  END DO
!
!*    TIME FILTER FOR TURBULENT KINETIC ENERGY
!
     IF(.NOT.lstart) THEN
       zeps=eps
     ELSE
       zeps=0._dp
     END IF
     DO 397 jk=ktdia,klev
       DO 396 jl=1,kproma
         ptkem1(jl,jk)=ptkem(jl,jk)                                    &
                   +zeps*(ptkem1(jl,jk)-2._dp*ptkem(jl,jk)+ptke(jl,jk))
         ptkem(jl,jk)=ptke(jl,jk)
396     END DO
397  END DO
!
!     ------------------------------------------------------------------
!
!*         4.     DIFFUSION IMPLICIT COMPUTATIONS FOR MOMENTUM.
!
!*         4.1     SETTING OF RIGHT HAND SIDES.
!
     DO 412 jk=itop,klevm1
        DO 411 jl=1,kproma
           zudif(jl,jk)=ztpfac2*pum1(jl,jk)
           zvdif(jl,jk)=ztpfac2*pvm1(jl,jk)
411     END DO
412  END DO

        DO 413 jl=1,kproma
           zudif(jl,klev)=ztpfac2*(pum1(jl,klev)-                      &
                         pocu(jl)*(1._dp-pfrl(jl)))
           zvdif(jl,klev)=ztpfac2*(pvm1(jl,klev)-                      &
                         pocv(jl)*(1._dp-pfrl(jl)))
413     ENDDO
!
!*         4.2     TOP LAYER ELIMINATION.
!
     DO 421 jl=1,kproma
        zqdp=1._dp/(paphm1(jl,itopp1)-paphm1(jl,itop))
        zdisc=1._dp/(1._dp+zcfm(jl,itop)*zqdp)
        zebsm(jl,itop)=zdisc*(zcfm(jl,itop)*zqdp)
        zudif(jl,itop)=zdisc*zudif(jl,itop)
        zvdif(jl,itop)=zdisc*zvdif(jl,itop)
421  END DO
!
!*         4.3     ELIMINATION FOR MIDDLE LAYERS.
!
     DO 432 jk=itopp1,klevm1
        DO 431 jl=1,kproma
           zqdp=1._dp/(paphm1(jl,jk+1)-paphm1(jl,jk))
           zfac=zcfm(jl,jk-1)*zqdp
           zdisc=1._dp/(1._dp+zfac*(1._dp-zebsm(jl,jk-1))              &
                                                   +zcfm(jl,jk)*zqdp)
           zebsm(jl,jk)=zdisc*(zcfm(jl,jk)*zqdp)
           zudif(jl,jk)=zdisc*(zudif(jl,jk)+zfac*zudif(jl,jk-1))
           zvdif(jl,jk)=zdisc*(zvdif(jl,jk)+zfac*zvdif(jl,jk-1))
431     END DO
432  END DO
!
!*         4.4     BOTTOM LAYER ELIMINATION.
!
     DO 441 jl=1,kproma
        zqdp=1._dp/(paphm1(jl,klevp1)-paphm1(jl,klev))
        zfac=zcfm(jl,klevm1)*zqdp
        zdisc=1._dp/(1._dp+zfac*(1._dp-zebsm(jl,klevm1))               &
                                                  +zcfm(jl,klev)*zqdp)
        zudif(jl,klev)=zdisc*(zudif(jl,klev)+zfac*zudif(jl,klevm1))
        zvdif(jl,klev)=zdisc*(zvdif(jl,klev)+zfac*zvdif(jl,klevm1))
441  END DO
!
!*         4.5     BACK-SUBSTITUTION.
!
     DO 452 jk=klevm1,itop,-1
        DO 451 jl=1,kproma
           zudif(jl,jk)=zudif(jl,jk)+zebsm(jl,jk)*zudif(jl,jk+1)
           zvdif(jl,jk)=zvdif(jl,jk)+zebsm(jl,jk)*zvdif(jl,jk+1)
451     END DO
452  END DO
!
!*         4.6     INCREMENTATION OF U AND V TENDENCIES AND STORAGE OF
!*                 THE DISSIPATION.
!
     DO 461 jl=1,kproma
        zvidis(jl)=0._dp
461  END DO
!
     DO 471 jk=itop,klev
        DO 462 jl=1,kproma
           zdudt=(zudif(jl,jk)-ztpfac2*pum1(jl,jk))*zcons13
!!$#ifndef MESSYTENDENCY
           pvom(jl,jk)=pvom(jl,jk)+zdudt
!!$#else
!!$           zdudt_2d(jl,jk)=zdudt
!!$#endif
           zdvdt=(zvdif(jl,jk)-ztpfac2*pvm1(jl,jk))*zcons13
!!$#ifndef MESSYTENDENCY
           pvol(jl,jk)=pvol(jl,jk)+zdvdt
!!$#else
!!$           zdvdt_2d(jl,jk)=zdvdt
!!$#endif
           zdis(jl,jk)=0.5_dp*((ztpfac2*pum1(jl,jk)-zudif(jl,jk))      &
                              *(ztpfac4*pum1(jl,jk)+zudif(jl,jk))      &
                              +(ztpfac2*pvm1(jl,jk)-zvdif(jl,jk))      &
                              *(ztpfac4*pvm1(jl,jk)+zvdif(jl,jk)))
           zvidis(jl)=zvidis(jl)+                                      &
                      zdis(jl,jk)*(paphm1(jl,jk+1)-paphm1(jl,jk))
462     END DO
471  END DO
!!$#ifdef MESSYTENDENCY
!!$     CALL mtend_add_l (my_handle, mtend_id_u, px = zdudt_2d)
!!$     CALL mtend_add_l (my_handle, mtend_id_v, px = zdvdt_2d)
!!$#endif
!
!*         4.8     UPDATING OF Z0 FOR OPEN SEA.
!
     DO  481 jl=1,kproma
      paz0w(jl)=MAX(zepzzo,                                            &
                    zcons14*zcfmw(jl)                                  &
                   *SQRT(zudif(jl,klev)**2+zvdif(jl,klev)**2)          &
                   *ptm1(jl,klev)                                      &
                   *(1._dp+vtmpc1*pqm1(jl,klev)-zx(jl,klev))           &
                   /paphm1(jl,klevp1))
      paz0i(jl)=cz0ice
      paz0(jl)=pfrl(jl)*paz0l(jl)+pfrw(jl)*paz0w(jl)                   &
              +pfri(jl)*paz0i(jl)
!
!     windstress
!
!      land   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      pustrl(jl)=zcons15*zcfml(jl)*zudif(jl,klev)
      pvstrl(jl)=zcons15*zcfml(jl)*zvdif(jl,klev)
!
!      water   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      pustrw(jl)=zcons15*zcfmw(jl)*zudif(jl,klev)
      pvstrw(jl)=zcons15*zcfmw(jl)*zvdif(jl,klev)
!
!      ice   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      pustri(jl)=zcons15*zcfmi(jl)*zudif(jl,klev)
      pvstri(jl)=zcons15*zcfmi(jl)*zvdif(jl,klev)
!
!      average   +++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      pustr(jl)=pustr(jl)+(pfrl(jl)*pustrl(jl)                         &
                          +pfrw(jl)*pustrw(jl)                         &
                          +pfri(jl)*pustri(jl))*zdtime
      pvstr(jl)=pvstr(jl)+(pfrl(jl)*pvstrl(jl)                         &
                          +pfrw(jl)*pvstrw(jl)                         &
                          +pfri(jl)*pvstri(jl))*zdtime
      pvdis(jl)=pvdis(jl)+zdtime*zcons15*zvidis(jl)
481  END DO
!
!     ------------------------------------------------------------------
!
!*         5.     DIFFUSION IMPLICIT COMPUTATIONS FOR HEAT (S.+L.).
!
500  CONTINUE
     DO 502 jk=1,klev
        DO 501 jl=1,kproma
           ztdif(jl,jk)=0._dp
           zqdif(jl,jk)=0._dp
           zxldif(jl,jk)=0._dp
           zxidif(jl,jk)=0._dp
           zvardif(jl,jk)=0._dp
501     END DO
502  END DO
!
!!$#ifndef MESSY
     DO 505 jt=1,trlist% ntrac
        IF (trlist% ti(jt)% nvdiff /= 1) CYCLE
        DO 504 jk=1,klev
           DO 503 jl=1,kproma
              zxtdif(jl,jk,jt)=0._dp
503        END DO
504     END DO
505  END DO
!!$#else
!!$       ! mz_rs_20040329+
!!$       DO jt = 1, ntrac_gp
!!$         IF (ti_gp(jt)%tp%meta%cask_i(I_Vdiff) == ON) THEN
!!$
!!$           DO jk=1,klev
!!$             DO jl=1,kproma
!!$               zxtdif(jl,jk,jt)=0._dp
!!$             END DO
!!$           END DO
!!$
!!$         END IF
!!$       END DO
!!$       ! mz_rs_20040329-
!!$#endif
!
!*         5.1     SETTING OF RIGHT HAND SIDES.
!
510  CONTINUE
     DO 512 jk=itop,klev
        DO 511 jl=1,kproma
           ztdif(jl,jk)=ztpfac2*zcptgz(jl,jk)
           zqdif(jl,jk)=ztpfac2*pqm1(jl,jk)
           zxldif(jl,jk)=ztpfac2*pxlm1(jl,jk)
           zxidif(jl,jk)=ztpfac2*pxim1(jl,jk)
           zvardif(jl,jk)=ztpfac2*pxvar(jl,jk)
511     END DO
512  END DO
!
!!$#ifndef MESSY
     DO 518 jt=1,trlist% ntrac
        IF (trlist% ti(jt)% nvdiff /= 1) CYCLE
        DO 516 jk=itop,klev
           DO 514 jl=1,kproma
               ! op_ck_20031001+
               zxtvn(jl,jk,jt)=pxtm1(jl,jk,jt)+pxtte(jl,jk,jt)*ztmst
               !!$ zxtdif(jl,jk,jt)=ztpfac2*pxtm1(jl,jk,jt)
               zxtdif(jl,jk,jt)=ztpfac2*zxtvn(jl,jk,jt)
               ! op_ck_20031001-
514        END DO
516     END DO
518  END DO
!!$#else
!!$#ifdef MESSYTENDENCY
!!$      call  mtend_get_start_l(mtend_id_tracer, v0t = zxtvn)
!!$#endif
!!$       ! mz_rs_20040329+
!!$       DO jt = 1, ntrac_gp
!!$         IF (ti_gp(jt)%tp%meta%cask_i(I_Vdiff) == ON) THEN
!!$
!!$           DO jk=itop,klev
!!$             DO jl=1,kproma
!!$               ! op_ck_20031001+
!!$#ifndef MESSYTENDENCY
!!$               zxtvn(jl,jk,jt)=pxtm1(jl,jk,jt)+pxtte(jl,jk,jt)*ztmst
!!$#endif
!!$               !!$ zxtdif(jl,jk,jt)=ztpfac2*pxtm1(jl,jk,jt)
!!$               zxtdif(jl,jk,jt)=ztpfac2*zxtvn(jl,jk,jt)
!!$               ! op_ck_20031001-
!!$             END DO
!!$           END DO
!!$         END IF
!!$       END DO
!!$       ! mz_rs_20040329-
!!$#endif
!
!*         5.2     TOP LAYER ELIMINATION.
!
520  CONTINUE
!
     DO 521 jl=1,kproma
        zqdp=1._dp/(paphm1(jl,itopp1)-paphm1(jl,itop))
        zdisc=1._dp/(1._dp+zcfh(jl,itop)*zqdp)
        zebsh(jl,itop)=zdisc*(zcfh(jl,itop)*zqdp)
        zdisv=1._dp/(1._dp+zcfv(jl,itop)*zqdp)
        zebsv(jl,itop)=zdisv*(zcfv(jl,itop)*zqdp)
        ztdif(jl,itop)=zdisc*ztdif(jl,itop)
        zqdif(jl,itop)=zdisc*zqdif(jl,itop)
        zxldif(jl,itop)=zdisc*zxldif(jl,itop)
        zxidif(jl,itop)=zdisc*zxidif(jl,itop)
        zvardif(jl,itop)=zdisv*zvardif(jl,itop)
521  END DO
!
!!$#ifndef MESSY
     DO 528 jt=1,trlist% ntrac
        IF (trlist% ti(jt)% nvdiff /= 1) CYCLE
        DO 526 jl=1,kproma
           zqdp=1._dp/(paphm1(jl,itopp1)-paphm1(jl,itop))
           zdisc=1._dp/(1._dp+zcfh(jl,itop)*zqdp)
           zxtdif(jl,itop,jt)=zdisc*zxtdif(jl,itop,jt)
526     END DO
528  END DO
!!$#else
!!$      ! mz_rs_20040329+
!!$      DO jt = 1, ntrac_gp
!!$        IF (ti_gp(jt)%tp%meta%cask_i(I_Vdiff) == ON) THEN
!!$          DO jl=1,kproma
!!$            zqdp=1._dp/(paphm1(jl,itopp1)-paphm1(jl,itop))
!!$            zdisc=1._dp/(1._dp+zcfh(jl,itop)*zqdp)
!!$            zxtdif(jl,itop,jt)=zdisc*zxtdif(jl,itop,jt)
!!$          END DO
!!$        END IF
!!$      END DO
!!$      ! mz_rs_20040329-
!!$#endif
!
!*         5.3     ELIMINATION FOR MIDDLE LAYERS.
!
530  CONTINUE
!
     DO 532 jk=itopp1,klevm1
        DO 531 jl=1,kproma
           zqdp=1._dp/(paphm1(jl,jk+1)-paphm1(jl,jk))
           zfac=zcfh(jl,jk-1)*zqdp
           zdisc=1._dp/(1._dp+zfac*(1._dp-zebsh(jl,jk-1))              &
                                                    +zcfh(jl,jk)*zqdp)
           zebsh(jl,jk)=zdisc*(zcfh(jl,jk)*zqdp)
           zfav=zcfv(jl,jk-1)*zqdp
           zdisv=1._dp/(1._dp+zfav*(1.-zebsv(jl,jk-1))                 &
                                                    +zcfv(jl,jk)*zqdp)
           zebsv(jl,jk)=zdisv*(zcfv(jl,jk)*zqdp)
           ztdif(jl,jk)=zdisc*(ztdif(jl,jk)+zfac*ztdif(jl,jk-1))
           zqdif(jl,jk)=zdisc*(zqdif(jl,jk)+zfac*zqdif(jl,jk-1))
           zxldif(jl,jk)=zdisc*(zxldif(jl,jk)+zfac*zxldif(jl,jk-1))
           zxidif(jl,jk)=zdisc*(zxidif(jl,jk)+zfac*zxidif(jl,jk-1))
           zvardif(jl,jk)=zdisv*(zvardif(jl,jk)+zfav*zvardif(jl,jk-1))
531     END DO
532  END DO
!
!!$#ifndef MESSY
     DO 538 jt=1,trlist% ntrac
        IF (trlist% ti(jt)% nvdiff /= 1) CYCLE
        DO 536 jk=itopp1,klevm1
           DO 534 jl=1,kproma
              zqdp=1._dp/(paphm1(jl,jk+1)-paphm1(jl,jk))
              zfac=zcfh(jl,jk-1)*zqdp
              zdisc=1._dp/(1._dp+zfac*(1._dp-zebsh(jl,jk-1))           &
                                                    +zcfh(jl,jk)*zqdp)
              zxtdif(jl,jk,jt)=zdisc*        &
                              (zxtdif(jl,jk,jt)+zfac*zxtdif(jl,jk-1,jt))
534        END DO
536     END DO
538  END DO
!!$#else
!!$       ! mz_rs_20040329+
!!$       DO jt = 1, ntrac_gp
!!$         IF (ti_gp(jt)%tp%meta%cask_i(I_Vdiff) == ON) THEN
!!$           DO jk=itopp1,klevm1
!!$             DO jl=1,kproma
!!$               zqdp=1._dp/(paphm1(jl,jk+1)-paphm1(jl,jk))
!!$               zfac=zcfh(jl,jk-1)*zqdp
!!$               zdisc=1._dp/(1._dp+zfac*(1._dp-zebsh(jl,jk-1))+zcfh(jl,jk)*zqdp)
!!$               zxtdif(jl,jk,jt)=zdisc*        &
!!$                 (zxtdif(jl,jk,jt)+zfac*zxtdif(jl,jk-1,jt))
!!$             END DO
!!$           END DO
!!$         END IF
!!$       END DO
!!$       ! mz_rs_20040329-
!!$#endif
!
!*         5.4     BOTTOM LAYER ELIMINATION.
!
540  CONTINUE
!
     DO 541 jl=1,kproma
        zqdp=1._dp/(paphm1(jl,klevp1)-paphm1(jl,klev))
        zfac=zcfh(jl,klevm1)*zqdp
        zdisx=1._dp/(1._dp+zfac*(1._dp-zebsh(jl,klevm1)))
        zfav=zcfv(jl,klevm1)*zqdp
        zdisv=1._dp/(1._dp+zfav*(1._dp-zebsv(jl,klevm1)))
        zxldif(jl,klev)=zdisx*(zxldif(jl,klev)+zfac*zxldif(jl,klevm1))
        zxidif(jl,klev)=zdisx*(zxidif(jl,klev)+zfac*zxidif(jl,klevm1))
        zvardif(jl,klev)=zdisv*(zvardif(jl,klev)                       &
                                             +zfav*zvardif(jl,klevm1))
!
!*  CALCULATION OF THE EN AND FN COEFFICIENTS OF THE RICHTMYER-
!*  MORTON-SCHEME CONCERNING THE EQUATION:
!
!*  XN = EN * XS + FN
!
!*  WITH XN = S_ATM  OR  XN = QATM : ATM. VALUE OF S OR Q
!*  AND  XS = SSURF  OR  XS = QSAT : SURFACE VALUE OF S OR SAT. SPEC.
!*                                   HUM. AT THE SURFACE
!
!       land   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        zdiscl=1._dp/(1._dp+zfac*(1._dp-zebsh(jl,klevm1))              &
                                         +zcfhl(jl)*zqdp)
        zdisql=1._dp/(1._dp+zfac*(1._dp-zebsh(jl,klevm1))              &
                                         +zcair(jl)*zcfhl(jl)*zqdp)
        zetnl(jl)=zdiscl*zcfhl(jl)*zqdp
        zftnl(jl)=zdiscl*(ztdif(jl,klev)+zfac*ztdif(jl,klevm1))*ztpfac1
        zeqnl(jl)=zdisql*zcsat(jl)*zcfhl(jl)*zqdp
        zfqnl(jl)=zdisql*(zqdif(jl,klev)+zfac*zqdif(jl,klevm1))*ztpfac1
!
!       water   +++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        zdiscw=1._dp/(1._dp+zfac*(1._dp-zebsh(jl,klevm1))              &
                                         +zcfhw(jl)*zqdp)
        zdisqw=zdiscw
        zetnw(jl)=zdiscw*zcfhw(jl)*zqdp
        zftnw(jl)=zdiscw*(ztdif(jl,klev)+zfac*ztdif(jl,klevm1))*ztpfac1
        zeqnw(jl)=zdisqw*zcfhw(jl)*zqdp
        zfqnw(jl)=zdisqw*(zqdif(jl,klev)+zfac*zqdif(jl,klevm1))*ztpfac1
!
!       ice   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        zdisci=1._dp/(1._dp+zfac*(1._dp-zebsh(jl,klevm1))              &
                                         +zcfhi(jl)*zqdp)
        zdisqi=zdisci
        zetni(jl)=zdisci*zcfhi(jl)*zqdp
        zftni(jl)=zdisci*(ztdif(jl,klev)+zfac*ztdif(jl,klevm1))*ztpfac1
        zeqni(jl)=zdisqi*zcfhi(jl)*zqdp
        zfqni(jl)=zdisqi*(zqdif(jl,klev)+zfac*zqdif(jl,klevm1))*ztpfac1
541  END DO
!
!!$#ifndef MESSY
     DO 544 jt=1,trlist% ntrac
        IF (trlist% ti(jt)% nvdiff /= 1) CYCLE
        DO 543 jl=1,kproma
           zqdp=1._dp/(paphm1(jl,klevp1)-paphm1(jl,klev))
           zfac=zcfh(jl,klevm1)*zqdp
           zdisxt=1._dp/(1._dp+zfac*(1._dp-zebsh(jl,klevm1)))
           zxtdif(jl,klev,jt)=zdisxt*(zxtdif(jl,klev,jt)+              &
                                    ztmst*g*zqdp*zxtems(jl,jt)         &
                                   +zfac*zxtdif(jl,klevm1,jt))
543     END DO
544  END DO
!!$#else
!!$      ! mz_rs_20040329+
!!$      DO jt = 1, ntrac_gp
!!$        IF (ti_gp(jt)%tp%meta%cask_i(I_Vdiff) == ON) THEN
!!$          DO jl=1,kproma
!!$            zqdp=1._dp/(paphm1(jl,klevp1)-paphm1(jl,klev))
!!$            zfac=zcfh(jl,klevm1)*zqdp
!!$            zdisxt=1._dp/(1._dp+zfac*(1._dp-zebsh(jl,klevm1)))
!!$            ! mz_lg_20030702+
!!$            ! The term CVDIFTS (1/ZTPFAC2) is normally 
!!$            ! introduced in the calculation of the term ZXTEMS.
!!$            ! But since in this version of vdiff.f90,
!!$            ! the inclusion of the this term has been changed, it is not 
!!$            ! included yet in ZXTEMS 
!!$        !!$   zxtdif(jl,klev,jt)=zdisxt*(zxtdif(jl,klev,jt)+              &
!!$        !!$                            ztmst*g*zqdp*zxtems(jl,jt)         &
!!$        !!$                           +zfac*zxtdif(jl,klevm1,jt))
!!$            zxtdif(jl,klev,jt)=zdisxt*(zxtdif(jl,klev,jt)+          &
!!$              ztmst*g*zqdp*zxtems(jl,jt)+  &
!!$              zfac*zxtdif(jl,klevm1,jt)  )
!!$            ! mz_lg_20030702-
!!$          END DO
!!$        END IF
!!$      END DO
!!$      ! mz_rs_20040329-
!!$#endif
!
     DO 545 jl = 1,kproma
        zteff4=pfrl(jl)*ptslm1(jl)**4                                  &
              +pfri(jl)*ptsi(jl)**4                                    &
              +pfrw(jl)*ptsw(jl)**4
        ztrdown=pemter(jl,klevp1)+cemiss*stbo*zteff4
        ztrfll=ztrdown-cemiss*stbo*zteff4
        znetr(jl) =zsrfll(jl)+ztrfll
        zcdrag(jl)=zcons30*zcfhl(jl)    ! rho*ch*|v|
545  END DO
!
!    Land surface temperature (cp*Ts=zcptlnew) and
!    saturation specific humidity (qsat=zqslnew)
!
    IF (.NOT.lstart) THEN
       CALL surftemp(kproma, zcons29,              &
        cemiss, stbo, zcpq, zcons16, alv, als, & ! constants
        zftnl, zetnl, zfqnl, zeqnl, &! coefficients from the elimination
        zcptl, zqsl, zdqsl,         &! old values at the surface
        znetr, pgrndhflx,           &! other fluxes
        zcdrag, zcair, zcsat, pcvs, pgrndcapc, & ! diff. coefficients,
                                              ! cair and csat for evap
                                              ! and soil heat capacity
        lpland,                             & ! logical land mask
        zcptlnew, zqslnew)                    ! output
    ELSE
       DO 546 jl = 1,kproma
         zcptlnew(jl)=zcptl(jl)
         zqslnew(jl)=zqsl(jl)
546    END DO
    END IF
!
!   Land surface temperature and 'zcptlnew' correction for snowmelt
!
     DO 547 jl = 1,kproma
       zcpt=ztpfac2*zcptlnew(jl)+ztpfac3*zcptl(jl)
       IF (lpland(jl)) THEN
         ptsl(jl)=zcpt/zcpq(jl)
       ELSE
         ptsl(jl)=tmelt
       END IF
       IF (psn(jl).GT.csncri.OR.lpglac(jl)) THEN
          zcptlcorr=MIN(ptsl(jl),tmelt)*zcpq(jl)
          zcptlnew(jl)=(zcptlcorr-ztpfac3*zcptl(jl))/ztpfac2
       ENDIF
!
!   New land surface temperature for sensible heat flux and *radheat*
!
       ptslnew(jl)=zcptlnew(jl)/zcpq(jl)
547  END DO
!
!*  CALCULATION OF SKLEV AND QKLEV USING THE NEW SURFACE VALUES
!*  ZSNEW AND ZQSNEW WHICH WERE CALCULATED IN SUBROUTINE SURFTEMP
!
     DO 548 jl = 1,kproma
!
!    land   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        ztklevl=zetnl(jl)*zcptlnew(jl)+zftnl(jl)
        zqklevl=zeqnl(jl)*zqslnew(jl)+zfqnl(jl)
!
!    water   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        ztklevw=zetnw(jl)*zcptw(jl)+zftnw(jl)
        zqklevw=zeqnw(jl)*zqsw(jl)+zfqnw(jl)
!
!    ice   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        ztklevi=zetni(jl)*zcpti(jl)+zftni(jl)
        zqklevi=zeqni(jl)*zqsi(jl)+zfqni(jl)
!
!    Grid-mean dry static energy and specific humidity at the
!    'blending height' (here: lowest model level 'klev').
!
        zgtl=pfrl(jl)*zcfhl(jl)
        zgtw=pfrw(jl)*zcfhw(jl)
        zgti=pfri(jl)*zcfhi(jl)
        zgtsum=(zgtl+zgtw+zgti)/ztpfac2
        IF (pfrl(jl).LT.1._dp) THEN
           zgql=zgtl*zcair(jl)
        ELSE
           zgql=zgtl
        ENDIF
        zgqw=zgtw
        zgqi=zgti
        zgqsum=(zgql+zgqw+zgqi)/ztpfac2
        ztdif(jl,klev)=(zgtl*ztklevl+zgtw*ztklevw+zgti*ztklevi)/zgtsum
        zqdif(jl,klev)=(zgql*zqklevl+zgqw*zqklevw+zgqi*zqklevi)/zgqsum
548  END DO
!
!*         5.5     BACK-SUBSTITUTION.
!
550  CONTINUE
!
     DO 552 jk=klevm1,itop,-1
        DO 551 jl=1,kproma
           ztdif(jl,jk)=ztdif(jl,jk)+zebsh(jl,jk)*ztdif(jl,jk+1)
           zqdif(jl,jk)=zqdif(jl,jk)+zebsh(jl,jk)*zqdif(jl,jk+1)
           zxldif(jl,jk)=zxldif(jl,jk)+zebsh(jl,jk)*zxldif(jl,jk+1)
           zxidif(jl,jk)=zxidif(jl,jk)+zebsh(jl,jk)*zxidif(jl,jk+1)
           zvardif(jl,jk)=zvardif(jl,jk)+zebsv(jl,jk)*zvardif(jl,jk+1)
551     END DO
552  END DO
!
!!$#ifndef MESSY
     DO 558 jt=1,trlist% ntrac
        IF (trlist% ti(jt)% nvdiff /= 1) CYCLE
        DO 556 jk=klevm1,itop,-1
           DO 554 jl=1,kproma
              zxtdif(jl,jk,jt)=zxtdif(jl,jk,jt)+                       &
                               zebsh(jl,jk)*zxtdif(jl,jk+1,jt)
554        END DO
556     END DO
558  END DO
!!$#else
!!$       ! mz_rs_20040329+
!!$       DO jt = 1, ntrac_gp
!!$         IF (ti_gp(jt)%tp%meta%cask_i(I_Vdiff) == ON) THEN
!!$           DO jk=klevm1,itop,-1
!!$             DO jl=1,kproma
!!$               zxtdif(jl,jk,jt)=zxtdif(jl,jk,jt)+                       &
!!$                 zebsh(jl,jk)*zxtdif(jl,jk+1,jt)
!!$             END DO
!!$           END DO
!!$         END IF
!!$       END DO
!!$       ! mz_rs_20040329-
!!$#endif
!
!*         5.6     INCREMENTATION OF T AND Q TENDENCIES.
!
560  CONTINUE
!
     DO 571 jk=itop,klev
        DO 561 jl=1,kproma
           zqdif(jl,jk)=zqdif(jl,jk)+ztpfac3*pqm1(jl,jk)
           zdqdt=(zqdif(jl,jk)-pqm1(jl,jk))*zcons13
!!$#ifndef MESSYTENDENCY
           pqte(jl,jk)=pqte(jl,jk)+zdqdt
!!$#else
!!$           zdqdt_2d(jl,jk)=zdqdt
!!$#endif
           ztdif(jl,jk)=ztdif(jl,jk)+ztpfac3*zcptgz(jl,jk)
           zdtdt=((ztdif(jl,jk)+zdis(jl,jk)-pgeom1(jl,jk))             &
                 /(cpd*(1._dp+vtmpc2*zqdif(jl,jk)))-ptm1(jl,jk))*zcons13
!!$#ifndef MESSYTENDENCY
           ptte(jl,jk)=ptte(jl,jk)+zdtdt
!!$#else
!!$           zdtdt_2d(jl,jk)=zdtdt
!!$#endif
           zxldif(jl,jk)=zxldif(jl,jk)+ztpfac3*pxlm1(jl,jk)
           zxidif(jl,jk)=zxidif(jl,jk)+ztpfac3*pxim1(jl,jk)
           zdxlmdt=(zxldif(jl,jk)-pxlm1(jl,jk))*zcons13
           zdximdt=(zxidif(jl,jk)-pxim1(jl,jk))*zcons13
!!$#ifndef MESSYTENDENCY
           pxlte(jl,jk)=pxlte(jl,jk)+zdxlmdt
           pxite(jl,jk)=pxite(jl,jk)+zdximdt
!!$#else
!!$           zdxlmdt_2d(jl,jk)=zdxlmdt
!!$           zdximdt_2d(jl,jk)=zdximdt
!!$#endif
           pxvar(jl,jk)=zvardif(jl,jk)+ztpfac3*pxvar(jl,jk)
           pvdiffp(jl,jk)=zdqdt+zdxlmdt+zdximdt !store for production
561     END DO
571  END DO
!!$#ifdef MESSYTENDENCY
!!$     CALL mtend_add_l (my_handle, mtend_id_t, px = zdtdt_2d)
!!$     CALL mtend_add_l (my_handle, mtend_id_q, px = zdqdt_2d)
!!$     CALL mtend_add_l (my_handle, mtend_id_xl, px = zdxlmdt_2d)
!!$     CALL mtend_add_l (my_handle, mtend_id_xi, px = zdximdt_2d)
!!$#endif
!
!!$#ifndef MESSY
     IF (trlist% anyvdiff /= 0) THEN
        DO 577 jt=1,trlist% ntrac
           IF (trlist% ti(jt)% nvdiff /= 1) CYCLE
           DO 575 jk=itop,klev
              DO 573 jl=1,kproma
               ! op_ck_20031001+
           !!$      zxtdif(jl,jk,jt)=zxtdif(jl,jk,jt)+                    &
           !!$                       ztpfac3*pxtm1(jl,jk,jt)
           !!$      zdxtdt=(zxtdif(jl,jk,jt)-pxtm1(jl,jk,jt))*zcons13
               zxtdif(jl,jk,jt)=zxtdif(jl,jk,jt) +                    &
                 ztpfac3*zxtvn(jl,jk,jt)
               zdxtdt=(zxtdif(jl,jk,jt)-zxtvn(jl,jk,jt))*zcons13
               ! op_ck_20031001-
                 pxtte(jl,jk,jt)=pxtte(jl,jk,jt)+zdxtdt
573           END DO
575        END DO
577     END DO
     END IF
!!$#else       
!!$       ! mz_rs_20040329+
!!$       DO jt = 1, ntrac_gp
!!$         IF (ti_gp(jt)%tp%meta%cask_i(I_Vdiff) == ON) THEN
!!$           DO jk=itop,klev
!!$             DO jl=1,kproma
!!$               ! op_ck_20031001+
!!$           !!$      zxtdif(jl,jk,jt)=zxtdif(jl,jk,jt)+                    &
!!$           !!$                       ztpfac3*pxtm1(jl,jk,jt)
!!$           !!$      zdxtdt=(zxtdif(jl,jk,jt)-pxtm1(jl,jk,jt))*zcons13
!!$               zxtdif(jl,jk,jt)=zxtdif(jl,jk,jt) +                    &
!!$                 ztpfac3*zxtvn(jl,jk,jt) 
!!$               zdxtdt=(zxtdif(jl,jk,jt)-zxtvn(jl,jk,jt))*zcons13 
!!$               ! op_ck_20031001-
!!$               ! mz_ht_20030702+
!!$               ! discard negative values qqq
!!$               !zdxtdt = &
!!$               !  MAX(-pxtm1(jl,jk,jt)/time_step_len-pxtte(jl,jk,jt), zdxtdt)
!!$               IF (zxtvn(jl,jk,jt) > 0.0_DP) THEN
!!$                  zdxtdt = & 
!!$                      MAX(-zxtvn(jl,jk,jt)/time_step_len, zdxtdt)
!!$               ELSE
!!$                  zdxtdt = & 
!!$                      MAX(0.0_DP, zdxtdt)
!!$               END IF
!!$               ! mz_ht_20030702-
!!$#ifndef MESSYTENDENCY
!!$               pxtte(jl,jk,jt)=pxtte(jl,jk,jt)+zdxtdt
!!$#else
!!$               zdxtdt_3d(jl,jk,jt)=zdxtdt
!!$#endif
!!$             END DO
!!$           END DO
!!$         END IF
!!$       END DO
!!$       ! mz_rs_20040329-
!!$#ifdef MESSYTENDENCY
!!$       IF (ntrac_gp > 0) &
!!$            CALL mtend_add_l(my_handle, mtend_id_tracer, pxt=zdxtdt_3d)
!!$#endif
!!$#endif
!
! back out moisture flux
!
    DO jl=1,kproma
       zqflux(jl,itop)=0._dp
       zvarpr(jl,itop)=0._dp
       ztvlan=ptslm1(jl)*(1._dp+vtmpc1*zhsoil(jl)*zqsl(jl))
       ztvsea=ptsw(jl)*(1._dp+vtmpc1*zqsw(jl))
       ztvice=ptsi(jl)*(1._dp+vtmpc1*zqsi(jl))
       ztvh=pfrl(jl)*ztvlan+pfrw(jl)*ztvsea+pfri(jl)*ztvice
       zrho(jl,klevp1)=paphm1(jl,klevp1)/(rd*ztvh)        !surf density
       zqsurf=pfrl(jl)*zqsl(jl)*zcsat(jl)+pfrw(jl)*zqsw(jl)            &
              +pfri(jl)*zqsi(jl)
       zdqtot=(pqm1(jl,klev)+zx(jl,klev))-zqsurf
       zqshear(jl,klev)=zdqtot*g/pgeom1(jl,klev) !qshear in lev 19
    ENDDO !jl
    DO jk=itop+1,klevp1
       DO jl=1,kproma
          IF (jk<klevp1) THEN
             ztvh=(ptm1(jl,jk)*(1._dp+vtmpc1*pqm1(jl,jk)-zx(jl,jk))    &
                  +ptm1(jl,jk-1)*(1._dp+vtmpc1*pqm1(jl,jk-1)           &
                   -zx(jl,jk-1)))/2._dp
             zrho(jl,jk)=paphm1(jl,jk)/(rd*ztvh)
          ENDIF
          zrhodz=-(paphm1(jl,jk)-paphm1(jl,jk-1))/g
          zqflux(jl,jk)=zrhodz*pvdiffp(jl,jk-1)+zqflux(jl,jk-1)
          zvarpr(jl,jk)=zqshear(jl,jk-1)*zqflux(jl,jk)/zrho(jl,jk)
       ENDDO !jl
    ENDDO !jk
    DO jk=itop,klev
       DO jl=1,kproma
          pvdiffp(jl,jk)=(zvarpr(jl,jk)+zvarpr(jl,jk+1))/2._dp
          zhexp=EXP(1._dp-pgeom1(jl,jk)/pgeom1(jl,ihpbl(jl)))
          zlam=zzzlam+(zcons3-zzzlam)*zhexp
          IF(jk.GE.ihpbl(jl)) THEN
             zcons23=zcons25
          ELSE
             zcons23=zcons2/zlam
          END IF
          z2geomf=2._dp*pgeom1(jl,jk)
          zz2geo=zcons2*z2geomf
          zmix=zz2geo/(1._dp+zcons23*z2geomf)
          IF(jk.EQ.1) THEN
             ztkesq=SQRT(MAX(ztkemin,ptkem1(jl,1)))
          ELSE
             ztkesq=SQRT(MAX(ztkemin,0.5_dp*(ptkem1(jl,jk-1)           &
                                           +ptkem1(jl,jk))))
          END IF
          pvmixtau(jl,jk)=ztkesq/(zmix*zda1)
       ENDDO !jl
    ENDDO !jk
!
!*         5.8     Surface fluxes of heat and moisture
!
580  CONTINUE
!
     DO 581 jl=1,kproma
!
      zcoefl=zcons15*zcfhl(jl)
      zcoefw=zcons15*zcfhw(jl)
      zcoefi=zcons15*zcfhi(jl)
!
!*    Moisture fluxes
!
        zqnlev=zqdif(jl,klev)-ztpfac3*pqm1(jl,klev)
!
!       land   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        zzqs=ztpfac2*zqslnew(jl)
        zqhfll(jl)=zcoefl*(zcair(jl)*zqnlev-zcsat(jl)*zzqs)
!
!       water   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        zqhflw(jl)=zcoefw*(zqnlev-ztpfac2*zqsw(jl))

!       ice   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        zqhfli(jl)=zcoefi*(zqnlev-ztpfac2*zqsi(jl))

!       Area mean moisture flux (=evaporation)

        pqhfla(jl)=pfrl(jl)*zqhfll(jl)+pfrw(jl)*zqhflw(jl)             &
                  +pfri(jl)*zqhfli(jl)
!
!*    Sensible heat fluxes

        ztnlev=ztdif(jl,klev)-ztpfac3*zcptgz(jl,klev)
!
!       land   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        zqnew=pqm1(jl,klev)+zcair(jl)*(zqnlev*ztpfac1-pqm1(jl,klev))
        zzcpts=ztpfac2*ptslnew(jl)*cpd*(1._dp+vtmpc2*zqnew)
        zthfll(jl)=zcoefl*(ztnlev-zzcpts)
!
!       water   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        zthflw(jl)=zcoefw*(ztnlev-ztpfac2*zcptw(jl))
        zthflw(jl)=zthflw(jl)-ptsw(jl)*zcons16*zqhflw(jl)
!
!       ice   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        zthfli(jl)=zcoefi*(ztnlev-ztpfac2*zcpti(jl))
        zthfli(jl)=zthfli(jl)-ptsi(jl)*zcons16*zqhfli(jl)
!
!    Accumulated sensible heat flux and evaporation
!

        pevapl(jl)=zqhfll(jl)
        pevapw(jl)=zqhflw(jl)
        pevapi(jl)=zqhfli(jl)

        pahfsl(jl)=zthfll(jl)
        pahfsw(jl)=zthflw(jl)
        pahfsi(jl)=zthfli(jl)

!!$#ifdef MESSY
!!$        ! mz_mt_20031120+
!!$        ! SAVE pahfs of last timestep
!!$        pahfsm1(jl)=pahfs(jl)
!!$        ! mz_mt_20031120-
!!$        ! mz_ht_20050222+
!!$        s_heatflux(jl,krow) = pfrl(jl) * zthfll(jl)                    &
!!$                           + pfrw(jl) * zthflw(jl)                     &
!!$                           + pfri(jl) * zthfli(jl)
!!$        ! mz_ht_20050222-
!!$#endif
        pahfs(jl)=pahfs(jl)+(pfrl(jl)*zthfll(jl)                       &
                           +pfrw(jl)*zthflw(jl)                        &
                           +pfri(jl)*zthfli(jl))*zdtime
        pevap(jl)=pevap(jl)+(pfrl(jl)*zqhfll(jl)                       &
                           +pfrw(jl)*zqhflw(jl)                        &
                           +pfri(jl)*zqhfli(jl))*zdtime
!
!*    Potential evaporation over land (for snow and skin reservoir)
!
        pevapot(jl)=zcoefl*(zqnlev-zzqs)
!
!!$#ifdef MESSY
!!$        !mz_mt_20031120+
!!$        ! CALCULATION OF SURFACE KINEMATIC HEAT FLUX [mk/s]
!!$        dens= paphm1(jl,klevp1)/zrd/ptm1(jl,klev)    &
!!$            * (1.+vtmpc1*pqm1(jl,klev))                  !density [kg m-3]
!!$
!!$!        zsenkf_2d(jl,krow)=(pahfs(jl)-pahfsm1(jl))/(zdtime*dens*zcpd)
!!$        zsenkf_2d(jl,krow)= s_heatflux(jl,krow) /(dens*zcpd)
!!$        
!!$        ! CALCULATION OF SURFACE KINEMATIC MOISTURE FLUX [m/s]
!!$        zlatkf_2d(jl,krow)=((alv*pqhfla(jl))+als*pcvs(jl)*pevapot(jl))  &
!!$             /(dens*alv)
!!$        !mz_mt_20031120-
!!$#endif
!     Latent heat fluxes
!
!     land   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        zqvhfl=zqhfll(jl)-pcvs(jl)*pevapot(jl)
        pahfll(jl)=alv*zqvhfl+als*pcvs(jl)*pevapot(jl)
!
!     water
!
        pahflw(jl)=alv*zqhflw(jl)
!
!     ice   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        pahfli(jl)=als*zqhfli(jl)
!
!     Accumulated area mean latent heat flux
!
!!$#ifdef MESSY
!!$        ! mz_ht_20070421+
!!$        l_heatflux(jl,krow) = (pfrl(jl)*pahfll(jl)                     &
!!$                  +pfrw(jl)*pahflw(jl)                                 &
!!$                  +pfri(jl)*pahfli(jl))
!!$        ! mz_ht_20070421-
!!$#endif
        pahfl(jl)=pahfl(jl)                                            &
                  +(pfrl(jl)*pahfll(jl)                                &
                  +pfrw(jl)*pahflw(jl)                                 &
                  +pfri(jl)*pahfli(jl))*zdtime
!
!     Accumulated area weighted heat fluxes and evaporations
!                   for land/water/ice
!
       pahfslac(jl)=pahfslac(jl)+pfrl(jl)*pahfsl(jl)*zdtime
       pahfswac(jl)=pahfswac(jl)+pfrw(jl)*pahfsw(jl)*zdtime
       pahfsiac(jl)=pahfsiac(jl)+pfri(jl)*pahfsi(jl)*zdtime
!
       pahfllac(jl)=pahfllac(jl)+pfrl(jl)*pahfll(jl)*zdtime
       pahflwac(jl)=pahflwac(jl)+pfrw(jl)*pahflw(jl)*zdtime
       pahfliac(jl)=pahfliac(jl)+pfri(jl)*pahfli(jl)*zdtime
!
       pevaplac(jl)=pevaplac(jl)+pfrl(jl)*pevapl(jl)*zdtime
       pevapwac(jl)=pevapwac(jl)+pfrw(jl)*pevapw(jl)*zdtime
       pevapiac(jl)=pevapiac(jl)+pfri(jl)*pevapi(jl)*zdtime
581  END DO
!
!     Compute new t2m, t2m_max t2m_min
!
     DO 597 jl=1,kproma
!
!     land   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        lo1=zril(jl).GT.0._dp
        zrat=zhtq/pgeom1(jl,klev)
        zcbn=LOG(1._dp+(EXP(zbhnl(jl))-1._dp)*zrat)
        zcbs=-(zbhnl(jl)-zbhl(jl))*zrat
        zcbu=-LOG(1._dp+(EXP(zbhnl(jl)-zbhl(jl))-1._dp)*zrat)
        zred=(zcbn+MERGE(zcbs,zcbu,lo1))/zbhl(jl)
        zh2m=zcptl(jl)+zred*(zcptgz(jl,klev)-zcptl(jl))
        zt2l=(zh2m-zhtq)/(cpd*(1._dp+vtmpc2*pqm1(jl,klev)))
!
!     water   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        lo1=zriw(jl).GT.0._dp
        zcbn=LOG(1._dp+(EXP(zbnw(jl))-1._dp)*zrat)
        zcbs=-(zbnw(jl)-zbhw(jl))*zrat
        zcbu=-LOG(1._dp+(EXP(zbnw(jl)-zbhw(jl))-1._dp)*zrat)
        zred=(zcbn+MERGE(zcbs,zcbu,lo1))/zbhw(jl)
        zh2m=zcptw(jl)+zred*(zcptgz(jl,klev)-zcptw(jl))
        zt2w=(zh2m-zhtq)/(cpd*(1._dp+vtmpc2*pqm1(jl,klev)))
!
!     ice   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        lo1=zrii(jl).GT.0._dp
        zcbn=LOG(1._dp+(EXP(zbni(jl))-1._dp)*zrat)
        zcbs=-(zbni(jl)-zbhi(jl))*zrat
        zcbu=-LOG(1._dp+(EXP(zbni(jl)-zbhi(jl))-1._dp)*zrat)
        zred=(zcbn+MERGE(zcbs,zcbu,lo1))/zbhi(jl)
        zh2m=zcpti(jl)+zred*(zcptgz(jl,klev)-zcpti(jl))
        zt2i=(zh2m-zhtq)/(cpd*(1._dp+vtmpc2*pqm1(jl,klev)))
        ptemp2(jl)=pfrl(jl)*zt2l+pfrw(jl)*zt2w+pfri(jl)*zt2i
        pt2max(jl)=MAX(pfrl(jl)*zt2l+pfrw(jl)*zt2w                     &
                   +pfri(jl)*zt2i,pt2max(jl))
        pt2min(jl)=MIN(pfrl(jl)*zt2l+pfrw(jl)*zt2w                     &
                   +pfri(jl)*zt2i,pt2min(jl))
!!$#ifdef MESSY
!!$        t2l_1d(jl)=zt2l
!!$        t2w_1d(jl)=zt2w
!!$        t2i_1d(jl)=zt2i
!!$#endif
!
!           5.96   2M DEW POINT
!
        it = NINT(ptm1(jl,klev)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zqs1=tlucua(it)/papm1(jl,klev)
        zqs1=zqs1/(1._dp-vtmpc1*zqs1)
        zrh2m=MAX(zephum,pqm1(jl,klev)/zqs1)
!!$#ifdef MESSY
!!$        ! mz_rs_20030128+
!!$        rh_2m(jl,krow) = MIN(1.,zrh2m) ! set to a maximum of 1
!!$        ! mz_rs_20030128-
!!$#endif
!
!       land   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        lo=zt2l.GT.tmelt
        zcvm3=MERGE(c3les,c3ies,lo)
        zcvm4=MERGE(c4les,c4ies,lo)
        zaph2m=paphm1(jl,klevp1)*(1._dp-zhtq/(rd*zt2l                  &
                           *(1._dp+vtmpc1*pqm1(jl,klev)-zx(jl,klev))))
        it = NINT(zt2l*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zqs2=tlucua(it)/zaph2m
        zqs2=zqs2/(1._dp-vtmpc1*zqs2)
        zq2m=zrh2m*zqs2
        zfrac=LOG(zaph2m*zq2m/(c2es*(1._dp+vtmpc1*zq2m)))/zcvm3
        zdew2l=MIN(zt2l,(tmelt-zfrac*zcvm4)/(1._dp-zfrac))
!
!     water   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        lo=zt2w.GT.tmelt
        zcvm3=MERGE(c3les,c3ies,lo)
        zcvm4=MERGE(c4les,c4ies,lo)
        zaph2m=paphm1(jl,klevp1)*(1._dp-zhtq/(rd*zt2w                  &
                           *(1._dp+vtmpc1*pqm1(jl,klev)-zx(jl,klev))))
        it = NINT(zt2w*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zqs2=tlucua(it)/zaph2m
        zqs2=zqs2/(1._dp-vtmpc1*zqs2)
        zq2m=zrh2m*zqs2
        zfrac=LOG(zaph2m*zq2m/(c2es*(1._dp+vtmpc1*zq2m)))/zcvm3
        zdew2w=MIN(zt2w,(tmelt-zfrac*zcvm4)/(1._dp-zfrac))
!
!     ice   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        lo=zt2i.GT.tmelt
        zcvm3=MERGE(c3les,c3ies,lo)
        zcvm4=MERGE(c4les,c4ies,lo)
        zaph2m=paphm1(jl,klevp1)*(1._dp-zhtq/(rd*zt2i                  &
                           *(1._dp+vtmpc1*pqm1(jl,klev)-zx(jl,klev))))
        it = NINT(zt2i*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zqs2=tlucua(it)/zaph2m
        zqs2=zqs2/(1._dp-vtmpc1*zqs2)
        zq2m=zrh2m*zqs2
        zfrac=LOG(zaph2m*zq2m/(c2es*(1._dp+vtmpc1*zq2m)))/zcvm3
        zdew2i=MIN(zt2i,(tmelt-zfrac*zcvm4)/(1._dp-zfrac))
        pdew2(jl)=pfrl(jl)*zdew2l+pfrw(jl)*zdew2w+pfri(jl)*zdew2i

!
!*          5.97   10M WIND COMPONENTS, MAX 10M WINDSPEED
!
!    land   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        lo1=zril(jl).GT.0._dp
        zrat=zhuv/pgeom1(jl,klev)
        zcbn=LOG(1._dp+(EXP(zbnl(jl))-1._dp)*zrat)
        zcbs=-(zbnl(jl)-zbml(jl))*zrat
        zcbu=-LOG(1._dp+(EXP(zbnl(jl)-zbml(jl))-1._dp)*zrat)
        zred=(zcbn+MERGE(zcbs,zcbu,lo1))/zbml(jl)
        zu10l=zred*pum1(jl,klev)
        zv10l=zred*pvm1(jl,klev)
        zspeedl=SQRT(zu10l**2+zv10l**2)
!
!     water   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        lo1=zriw(jl).GT.0._dp
        zcbn=LOG(1._dp+(EXP(zbnw(jl))-1._dp)*zrat)
        zcbs=-(zbnw(jl)-zbmw(jl))*zrat
        zcbu=-LOG(1._dp+(EXP(zbnw(jl)-zbmw(jl))-1._dp)*zrat)
        zred=(zcbn+MERGE(zcbs,zcbu,lo1))/zbmw(jl)
        zu10w=zred*pum1(jl,klev)
        zv10w=zred*pvm1(jl,klev)
        zspeedw=SQRT(zu10w**2+zv10w**2)
        pwind10w(jl)=zred*SQRT((pum1(jl,klev)-pocu(jl))**2             &
                              +(pvm1(jl,klev)-pocv(jl))**2)
!
!     ice   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        lo1=zrii(jl).GT.0._dp
        zcbn=LOG(1._dp+(EXP(zbni(jl))-1._dp)*zrat)
        zcbs=-(zbni(jl)-zbmi(jl))*zrat
        zcbu=-LOG(1._dp+(EXP(zbni(jl)-zbmi(jl))-1._dp)*zrat)
        zred=(zcbn+MERGE(zcbs,zcbu,lo1))/zbmi(jl)
        zu10i=zred*pum1(jl,klev)
        zv10i=zred*pvm1(jl,klev)
        zspeedi=SQRT(zu10i**2+zv10i**2)
        pu10(jl)=pfrl(jl)*zu10l+pfrw(jl)*zu10w+pfri(jl)*zu10i
        pv10(jl)=pfrl(jl)*zv10l+pfrw(jl)*zv10w+pfri(jl)*zv10i
        pwimax(jl)=MAX(pwimax(jl),                                     &
             (pfrl(jl)*zspeedl+pfrw(jl)*zspeedw+pfri(jl)*zspeedi))
        pwind10(jl)=pwind10(jl)+zdtime*                                &
             (pfrl(jl)*zspeedl+pfrw(jl)*zspeedw+pfri(jl)*zspeedi)
!!$#ifdef MESSY
!!$        ! mz_rs_20030708+
!!$        ! Save the non-accumulated wind speed at 10m in physc stream.
!!$        ! It is calculated seperately for land, water and ice  
!!$        ! and weighted according to the respective surface fractions.
!!$        wind10_2d(jl,krow)=pfrl(jl)*zspeedl+pfrw(jl)*zspeedw+pfri(jl)*zspeedi
!!$        wind10w_2d(jl,krow)=pwind10w(jl) ! mz_ap_20070913
!!$        ! mz_rs_20030708-
!!$#endif
597  END DO

!!$#ifdef MESSY
!!$  ! mz_ak_20060816+
!!$  IF (lookupoverflow) THEN 
!!$     do jl=1,kproma
!!$        if ( NINT(ptm1(jl,klev)*1000.) <jptlucu1 .OR. &
!!$             NINT(ptm1(jl,klev)*1000.) >jptlucu2)     &
!!$             print*, 'vdiff(3) ptm1: ',ptm1(jl,klev) &
!!$                ,philon_2d(jl,krow),philat_2d(jl,krow)
!!$        if ( NINT(t2l_1d(jl)*1000.) <jptlucu1 .OR. &
!!$             NINT(t2l_1d(jl)*1000.) >jptlucu2)     &
!!$             print*, 'vdiff(3) t2l: ',t2l_1d(jl) &
!!$                ,philon_2d(jl,krow),philat_2d(jl,krow)
!!$        if ( NINT(t2w_1d(jl)*1000.) <jptlucu1 .OR. &
!!$             NINT(t2w_1d(jl)*1000.) >jptlucu2)     &
!!$             print*, 'vdiff(3) t2w: ',t2w_1d(jl)  &
!!$                ,philon_2d(jl,krow),philat_2d(jl,krow)
!!$        if ( NINT(t2i_1d(jl)*1000.) <jptlucu1 .OR. &
!!$             NINT(t2i_1d(jl)*1000.) >jptlucu2)     &
!!$             print*, 'vdiff(3) t2i: ',t2i_1d(jl) &
!!$                ,philon_2d(jl,krow),philat_2d(jl,krow)
!!$     enddo
!!$  ENDIF
!!$  ! mz_ak_20060816-
!!$#endif 
     IF (lookupoverflow) CALL lookuperror ('vdiff (3)   ')
!
!     ------------------------------------------------------------------
!
!*         6.     NECESSARY COMPUTATIONS IF SUBROUTINE IS BY-PASSED.
!
600  CONTINUE
!
  ELSE
     DO  601 jl=1,kproma
        pevapot(jl)=0._dp
        pqhfla(jl)=0._dp
601  END DO
  END IF
!
!     ------------------------------------------------------------------
!
  RETURN
END SUBROUTINE vdiff
