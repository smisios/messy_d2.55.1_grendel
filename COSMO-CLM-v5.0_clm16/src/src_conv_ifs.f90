MODULE src_conv_ifs

!-------------------------------------------------------------------------------
!
! Description:
!   The module "src_conv_ifs" provides the interface from the COSMO model to 
!   the ECMWF IFS moist convection scheme (which has to be included as the 
!   library "libconvifs.a").
!
! Warning:
!   The namelist parameters lconf_avg has a different effect here than
!   in the Tiedtke scheme. Only lconf_avg=.FALSE. has been tested.
!
! How2Use:
!   (1) compile the library libconvifs.a
!HJP IMK, 2009-11-18 Begin Change
!   related to the following line: include a blank between "/" and "*/ otherwise compilation on NEC failed
!   (2) compile COSMO with linked library and include *.mod/ *.o of the library
!HJP IMK, 2009-11-18 End   Change
!   (3) ltype_conv=4 in namelist PHYCTL
!
! Publications:
!   Brockhaus, P.; Bechtold, P.; Fuhrer, O.; Luethi, D. & Schaer, C. 
!    (in preparation): The ECMWF IFS convection scheme applied to a 
!    limited-area model; QJRMS
!   Bechtold, P.; Koehler, M.; Jung, T.; Reyes, F.D.; Leutbecher, M.; 
!    Rodwell, M.J.; Vitart, F. & Balsamo, G. (2008): Advances in simulating 
!    atmospheric variability with the ECMWF model: from synoptic to decadal 
!    time-scales; QJRMS, 134, 1337-1351
!
! In-depth Documentation:
!   http.//www.ecmwf.int/newsevents/training/
!   lecture_notes/pdf_files/PARAM/MoistConvection.pdf
!
!------------------------------------------------------------------------------
!
! Notes for updates with future IFS cycles:
!
!  1) Input parameters of CUMASTRN that are not (yet) used in IFS Cy33r1
!      may have to be added to accomodate for future cycles. They are 
!      currently set to -123.
!  2) The only changes in the original convection library from ECMWF are as
!      follows (they are marked with !!**$$):
!        a) rprcon in sucumf.F90 reset from 1.4E-3 to 1.5E-3
!        b) allowing negative surface geopotential heights in cubasen.F90
!
!------------------------------------------------------------------------------
!
! Current Code Owner: ETH Zurich, Daniel Luethi
!  phone:  +41 44 632773
!  fax:    +41 44 6321311
!  email:  daniel.luethi@env.ethz.ch
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.0        2009/07/10 Peter Brockhaus/Peter Bechtold
!  Initial release (ECMWF IFS version Cy33r1)
!
!            2009/12/23 Hans-Juergen Panitz, IMK/TRO
!  Use p0hl (reference pressure at half levels) for full consistency with
!  new reference atmosphere implementation (Guenther)
!  Initialize PTENT with zero instead of ttens. Otherwise results
!  of a restart run won't be consistent with a comparable continuous run
!
!            2010/06/30 Hans-Juergen Panitz, IMK/TRO
!  Comment the warnings related to lconf_avg to avoid its
!  printing by every CPU and at every time step
!  The warning  are now printed in organize_physics.f90 after reading the
!  Namelist-Group PHYCTL
!  The warning related to llhn is now printed only once by processor no. 0
! V4_25_clm1        2012/10/26 Anne Roches
!  Replaced qx-variables by using them from the tracer module
! V4_25_clm1        2012/10/26 B. Rockel, HZG
!  include #ifdef NUDGING in relationship to variable llhn
!  correct the name of src_tracer
!
!==============================================================================

  USE data_parameters,    ONLY : &
       ireals,     &
       iintegers
     ! both only used by organize_conv_ifs (libconvifs.a uses those in PARKIND1)

  USE parallel_utilities, ONLY :  &
       global_values

  USE data_parallel,      ONLY :  &
       num_compute,    &
       icomm_cart,     &
       my_cart_id,     &
       my_cart_pos, &
       isubpos,     &
       imp_integers,   &
       imp_reals

IMPLICIT NONE

  INTEGER (KIND=iintegers) ::  &
    izerror

  CHARACTER (LEN=80)       ::  &
    yzerrmsg

  LOGICAL, PRIVATE ::       &
    lfirsti = .TRUE.   ! switch for initialization

CONTAINS

!********************************************************************************
!
! The interface routine from COSMO/CLM to the ECMWF IFS convection routines
!
!********************************************************************************

SUBROUTINE organize_conv_ifs

  USE data_fields     , ONLY :   &
       u,         & ! zonal wind speed                              ( m/s ) (:,:,:,:)
       v,         & ! meridional wind speed                         ( m/s ) (:,:,:,:)
       w,         & ! vertical wind speed (defined on half levels)  ( m/s ) (:,:,:,:)
       t,         & ! temperature                                   (  k  ) (:,:,:,:)
       shfl_s,    & ! sensible heat flux (surface, downward)        (W/m2 ) (:,:)
       lhfl_s,    & ! latent heat flux (surface, downward)          (W/m2 ) (:,:)
       p0,        & ! reference pressure at full levels             ( Pa  ) (:,:,:)
!
!HJP 2009/12/23
       p0hl,      & ! reference pressure at half levels             ( Pa  ) (:,:,:)
!HJP 2009/12/23
!
       rho0,      & ! reference density at the full model levels    (kg/m3) (:,:,:)
       hhl,       & ! geometrical height of half levels             (  m  ) (:,:,:)
       dp0,       & ! reference pressure thickness of layer         ( Pa  ) (:,:,:)
       pp,        & ! deviation from the reference pressure         ( Pa  ) (:,:,:,:)
       ps,        & ! surface pressure                              ( Pa  ) (:,:,:)
       pptens,    & ! pressure tendency without sound-wave terms    (Pa/s ) (:,:,:)
       ttens,     & ! temperature tendency without sound-wave terms ( K/s ) (:,:,:)
       utens,     & ! u-tendency without sound-wave terms           (m/s2 ) (:,:,:)
       vtens,     & ! v-tendency without sound-wave terms           (m/s2 ) (:,:,:)
       tt_conv,   & ! temperature tendency due to convection        ( K/s ) (:,:,:)
       qvt_conv,  & ! specific humidity tendency due to convection  ( 1/s ) (:,:,:)
       qct_conv,  & ! cloud liquid water tendency due to convection ( 1/s ) (:,:,:)
       qit_conv,  & ! qi-tendency tendency due to convection        ( 1/s ) (:,:,:)
       ut_conv,   & ! u-tendency due to convection                  (m/s2 ) (:,:,:)
       vt_conv,   & ! v-tendency due to convection                  (m/s2 ) (:,:,:)
       prr_con,   & ! precipitation rate of rain, convective        (kg/m2s)(:,:)
       prs_con,   & ! precipitation rate of snow, convective        (kg/m2s)(:,:)
       prne_con,  & ! precipitation rate, no evaporat., convective  (kg/m2s)(:,:)
       clc_con,   & ! cloud cover due to convection                 (  1  ) (:,:,:)
       top_con,   & ! level index of convective cloud top (real number)     (:,:)
       bas_con,   & ! level index of convective cloud base (real number)    (:,:)
       mflx_con,  & ! cloud base massflux                           (kg/m2s)(:,:)
       cape_con,  & ! convective available energy                   (J/kg ) (:,:)
       clw_con,   & ! convective cloud liquid water                 (kg/kg?)(:,:,:)

       crlat,     &
       llandmask    ! landpoint mask (:,:,:)

  USE data_runcontrol , ONLY :   &
       ! timestep variables:
       nstart,       & ! first time step of the forecast
       ntstep,       & ! actual time step
       nold,         & ! corresponds to ntstep-1
       nnow,         & ! corresponds to ntstep
       nnew,         & ! corresponds to ntstep + 1
       nincconv,     & ! time step increment for running the convection scheme

       ! physics:
       itype_gscp,   & ! type of grid-scale precipitation physics

       ! logicals:
       l2tls,        & ! forecast with 2-TL integration scheme
       lconf_avg       ! average convective forcings

  USE data_modelconfig, ONLY :   &
       ie,     & ! number of grid points in zonal direction (for the active processor)
       je,     & ! number of grid points in meridional direction
       ke,     & ! number of vertical levels
       istart, & ! start index (zonal) for the forecast of w, t, qd, qw and pp
       iend,   & ! end index (zonal) for the forecast of w, t, qd, qw and pp
       jstart, & ! start index (meridional) for the forecast of w, t, qd, qw and pp
       jend,   & ! end index (meridional) for the forecast of w, t, qd, qw and pp
       dlon,   & ! grid point distance in zonal direction (in degrees)
       dlat,   & ! grid point distance in meridional direction (in degrees)
       dt,     & ! long time-step
       dt2       ! 2 * dt

  USE data_constants  , ONLY :   &
       pi,      & ! circle constant 3.14159...
       r_earth, & ! mean radius of the earth
       g,       & ! acceleration due to gravity
       lh_v,    & ! latent heat of vapourization (J/kg)
       lh_f,    & ! latent heat of fusion (J/kg)
       cpdr       ! 1 / cp_d (1/(J/kgK))

  USE pp_utilities,             ONLY :  &
       calomega   ! routine that computes vertical velocity in pressure coordinates

#ifdef NUDGING
  USE data_lheat_nudge, ONLY :  &
       llhn       ! on/off switch for latent heat nudging, for warning message only
#endif

  ! global_values() has to be applied to some ECMWF IFS variables in libconvifs.a/yoecumf:
  USE yoecumf, ONLY :  &
       RTAU, NJKT1, NJKT2, NJKT3, NJKT4, NJKT5

  USE src_tracer,      ONLY : trcr_get, trcr_errorstr 

  USE environment,      ONLY : model_abort


  IMPLICIT NONE

  INTEGER (KIND=iintegers)    :: &
       KIDIA,       &    ! start point (zonal direction)
       KFDIA,       &    ! end point (zonal direction)
       KLON,        &    ! number of grid points per j-slice
       KTDIA,       &    ! start of vertical loop
       KLEV,        &    ! number of vertical levels
       KTRAC,       &    ! number of chemical tracers
       KSTEP,       &    ! NOT USED (initialized with -1234_iintegers)
       KSTART,      &    ! NOT USED (initialized with -1234_iintegers)
       KTYPE (ie),  &    ! type of convection (0=none,1=penetrative,2=shallow,3=midlevel)
       KCBOT (ie),  &    ! cloud base level
       KCTOP (ie),  &    ! cloud top level
       KBOTSC (ie), &    ! cloud base level for shallow cumulus
       NSMAX             ! approximate equivalent spectral resolution

  REAL (KIND=ireals) ::    &
       PTSPHY,             &  ! time step for convection
       PTEN (ie,ke),       &  ! environment temperature (K)
       PQEN (ie,ke),       &  ! environment specific humidity (kg/kg)
       PUEN (ie,ke),       &  ! environment u-velocity (m/s)
       PVEN (ie,ke),       &  ! environment v-velocity (m/s)
       PLITOT (ie,ke),     &  ! grid mean liquid water and ice content
       PVERVEL (ie,ke),    &  ! vertical velocity (Pa/s)
       PQHFL (ie,ke+1),    &  ! moisture flux w/o snow evaporation   (kg/m2, LH/L_v)
       PAHFS (ie,ke+1),    &  ! sensible heat flux (W/m2)
       PSSTRU (ie),        &  ! NOT USED
       PSSTRV (ie),        &  ! NOT USED
       PAP (ie,ke),        &  ! provisional pressure on full levels (Pa)
       PAPH (ie,ke+1),     &  ! provisional pressure on half levels (Pa)
       PGEO (ie,ke),       &  ! provisional geopotential on full levels (m2/s2)
       PGEOH (ie,ke+1),    &  ! provisional geopotential on half levels (m2/s2)
       PTENT (ie,ke),      &  ! temperature tendency (K/s)
       PTENQ (ie,ke),      &  ! moisture tendency (kg/(kg*s))
       PTENU (ie,ke),      &  ! u-wind tendency (m/s2)
       PTENV (ie,ke),      &  ! v-wind tendency (m/s2)
       PTENL (ie,ke),      &  ! liquid water tendency (kg/(kg*s))
       PTENI (ie,ke),      &  ! ice condensate tendency (kg/(kg*s))
       PTENC (ie,ke),      &  ! NOT USED (chemical tracer tendency, 1/s)
       PTU (ie,ke),        &  ! temperature in updrafts (K)
       PQU (ie,ke),        &  ! specific humidity in updrafts (kg/kg)
       PLU (ie,ke),        &  ! liquid water + ice content in updrafts (kg/kg)
       PMFLXR (ie,ke+1),   &  ! convective rain flux (m/s !!)
       PMFLXS (ie,ke+1),   &  ! convective snow flux (m/s !!)
       PRAIN (ie),         &  ! total precip. flux without evaporation (kg/(m2*s))
       PMFU (ie,ke),       &  ! massflux updrafts (kg/(m2*s))
       PMFD (ie,ke),       &  ! massflux downdrafts (kg/(m2*s))
       PMFUDE_RATE (ie,ke),&  ! updraft detrainment rate (kg/(m3*s))
       PMFDDE_RATE (ie,ke),&  ! downdraft detrainment rate (kg/(m3*s))
       PCAPE (ie),         &  ! convective available potential energy (J/kg)
       PCEN (ie,ke,1),     &  ! provisional tracer concentration (kg/kg), last dim = KTRAC, but >0
       PREF (ke),          &  ! mean pressure on model levels (used for NJKT* in SUCUMF)
       ZDX                    ! model resolution (m)

  LOGICAL  :: &
       LDLAND (ie), &   ! land-sea mask for a j-slice
       LDCUM (ie),  &   ! grid points with convection
       LDSC (ie)        ! grid points with boundary layer clouds

  REAL (KIND=ireals) ::     &
       zdt, &
       zcent, zside, zedge, &  ! weights for horizontal averaging
       omega (ie,je,ke),    &  ! vertical velocity (Pa/s)
       lhflx (ie,ke+1),     &  ! latent heat flux (W/m2, downward)
       shflx (ie,ke+1),     &  ! sensible heat flux (W/m2, downward)
       ! tcorr (ie,ke),       &  ! temperature correction due to qi+qc->PLITOT
       ppres (ie,je,ke),    &  ! aux. pressure variable
       ppress (je,ke),      &  ! aux. pressure variable
       ztop, zbas,          &
       zdtbot, zdttop

  INTEGER (KIND=iintegers) ::  &
       i, j, k, nx, km1, &
       mtop, mbas


  REAL (KIND=ireals), POINTER :: &
    qv     (:,:,:)   => NULL() , &     ! QV at nx
    qv_tens(:,:,:)   => NULL()         ! QV tendency

  CHARACTER (LEN=25)       :: yzroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  izerror  = 0
  yzerrmsg = '   '
  yzroutine= 'organize_conv_ifs'


  ! initialization of auxiliary variables
  zcent = 0.2500_ireals   ! centre weight in a nine point stencil
  zside = 0.1250_ireals   ! weight for side points
  zedge = 0.0625_ireals   ! weight for edge points

  ! initialize input parameters of CUMASTRN that are NOT USED:
  KSTEP      = -123_iintegers
  KSTART     = -123_iintegers
  PSSTRU (:) = -123.0_ireals
  PSSTRV (:) = -123.0_ireals
!
!HJP 2010-06-30 Begin
! Comment the subsequent warning
! it now appearsin organize_phyiscs.f90
!
  ! warnings
! IF (lconf_avg) THEN
!   PRINT *, '*** WARNING: Namelist parameter lconf_avg has a different effect ***'
!   PRINT *, '***          on the IFS scheme than with the Tiedtke scheme.     ***'
!   PRINT *, '***          Only lconf_avg=.FALSE. has yet been tested.         ***'
! ENDIF
!
! Subsequent warning is, hopefully printed only once by Processor No. 0
#ifdef NUDGING
  IF (ntstep == nstart) THEN
   IF (my_cart_id == 0 ) then
    IF (llhn) THEN
     PRINT *, '*** WARNING: Tiedtke-IFS convection scheme may  ***'
     PRINT *, '***          conflict with latent heat nudging! ***'
    ENDIF
   ENDIF
  ENDIF
#endif
!
!HJP 2010-06-30 End

  ! Select timelevel and timestep of the computation
  IF ( l2tls ) THEN
     ! Runge-Kutta
     nx  = nnow
     zdt = dt
  ELSE
     ! Leapfrog
     nx  = nold
     zdt = dt2
  ENDIF


  ! Retrieve the required microphysics tracers
  CALL trcr_get(izerror, 'QV', ptr_tens = qv_tens, ptr_tlev = nx, ptr = qv)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF


  ! initialize integer input parameters:
  KTRAC = 0_iintegers  ! tracers switched off, overpowers LMFTRAC=.TRUE. in SUCUMF
  KTDIA = 1_iintegers  ! start of vertical loop (used only in SUCUMF)
  KLEV  = ke
  KLON  = ie
  KIDIA = istart
  KFDIA = iend

  ! initialize real input parameter:
  PTSPHY = zdt*nincconv
  PCEN (:,:,:) = 0.0_ireals  ! chemical tracer concentration if(KTRAC.ne.0)

  ! mean pressure of each model level PREF (used for NJKT* in SUCUMF)
  DO k = 1, ke
    DO  j = 1, je
      ppres (:,j,k) = p0(:,j,k) + pp(:,j,k,nx)
      ppress (j,k) = sum(ppres(:,j,k))/real(ie)
    ENDDO
    PREF (k) = sum(ppress(:,k))/real(je)
  ENDDO

  ! run initialization routines:
  IF (lfirsti) THEN
    CALL SUCST(54,20020211,0,0)
    ! calculate approximate model resolution (m):
    !  (WARNING: assumes a small variance of resolution over the domain)
    ZDX=2*pi*r_earth*sqrt((dlon/360.0_ireals)*(dlat/360_ireals)* &
      crlat(INT(je/2.0_ireals),1))
    ! make sure ZDX (and thereby RTAU) is identical for all processors:
    IF (num_compute > 1) THEN
      ZDX = ZDX / num_compute
      CALL global_values(ZDX,1,'SUM',imp_reals,icomm_cart,-1,yzerrmsg,izerror)
    ENDIF
    NSMAX=INT(3.14159*r_earth/ZDX) ! cartesian -> spectral resolution
    CALL SUCUMF(NSMAX,KLEV,PREF)
    ! make sure some parameters are identical for all processors:
    IF (num_compute > 1) THEN
      CALL global_values(NJKT1,1,'MIN',imp_integers,icomm_cart,-1,yzerrmsg,izerror)
      CALL global_values(NJKT2,1,'MAX',imp_integers,icomm_cart,-1,yzerrmsg,izerror)
      CALL global_values(NJKT3,1,'MIN',imp_integers,icomm_cart,-1,yzerrmsg,izerror)
      CALL global_values(NJKT4,1,'MIN',imp_integers,icomm_cart,-1,yzerrmsg,izerror)
      CALL global_values(NJKT5,1,'MIN',imp_integers,icomm_cart,-1,yzerrmsg,izerror)
    ENDIF
    CALL SU_YOETHF
    CALL SUPHLI
    CALL SUVDF
    CALL SUVDFS
    CALL SUCLDP
    lfirsti=.FALSE.
  ENDIF

  ! calculate vertical velocity omega (Pa/s):
  omega (:,:,:) = 0.0_ireals
  CALL calomega (omega(:,:,:), pp(:,:,:,nnew),          &
                 pp(:,:,:,nnow), pptens(:,:,:), w(:,:,:,nx),     &
                 rho0(:,:,:), ie, je, ke, dt, g )

  !!!!!!!!!!!!!!!!!!!!!
  ! South-north loop: !
  !!!!!!!!!!!!!!!!!!!!!

  DO j = jstart, jend

    LDLAND(:)=llandmask(:,j)

    ! latent and sensible heat flux forcing (surface fluxes only)
    lhflx (:,ke+1) = lhfl_s(:,j)
    shflx (:,ke+1) = shfl_s(:,j)
    lhflx (:,1:ke) = 0.0_ireals
    shflx (:,1:ke) = 0.0_ireals
    PQHFL (:,:) = lhflx(:,:)/lh_v  ! convention in ECMWF IFS: w*q*
    PAHFS (:,:) = shflx(:,:)       ! convention in ECMWF IFS: c_p w*T*

    ! pressure and geopotential on full and half levels:
    DO k = 1, ke
      km1 = MAX ( 1, k-1 )
      PAP   (:,k) = p0(:,j,k) + pp(:,j,k,nx)
!
!HJP 2009/12/23
!     PAPH  (:,k) = p0(:,j,k) - 0.5_ireals*dp0(:,j,k) &
      PAPH  (:,k) = p0hl(:,j,k)                       &
                     + 0.5_ireals*(pp(:,j,k,nx) + pp(:,j,km1,nx))
!HJP 2009/12/23
!
      PGEO  (:,k) = 0.5_ireals*g*( hhl(:,j,k) + hhl(:,j,k+1) )
      PGEOH (:,k) = g*hhl(:,j,k) 
    ENDDO
    PAPH  (:,ke+1) = ps(:,j,nx)      ! surface pressure
    PGEOH (:,ke+1) = g*hhl(:,j,ke+1) ! surface geopotential

    ! input tendencies
!
!HJP 2009/12/23
!   PTENT (:,:) = ttens  (:,j,:)
    PTENT (:,:) = 0.0_ireals
!HJP 2009/12/23
!
    PTENQ (:,:) = qv_tens(:,j,:)
    PTENU (:,:) = utens  (:,j,:)  ! not necessarily needed as input
    PTENV (:,:) = vtens  (:,j,:)  ! not necessarily needed as input
    PTENL (:,:) = 0.0_ireals      ! recommendation by Peter Bechtold
    PTENI (:,:) = 0.0_ireals      ! recommendation by Peter Bechtold
    PTENC (:,:) = 0.0_ireals

    ! transfer state variables:
    PLITOT (:,:) = 0.0_ireals     ! recommendation by Peter Bechtold

    IF (.NOT.lconf_avg) THEN

      PTEN    (:,:) = t   (:,j,:,nx)
      PQEN    (:,:) = qv  (:,j,:)
      PUEN    (:,:) = u   (:,j,:,nx)  ! interpolation to unstaggered grid omitted
      PVEN    (:,:) = v   (:,j,:,nx)  ! interpolation to unstaggered grid omitted
      PVERVEL (:,:) = omega(:,j,:)

    ELSE

      DO k = 1, ke
        DO i = istart, iend
          PTEN (i,k) = zcent * t (i,j,k,nx)                    &
               + zside*( t (i-1,j  ,k,nx) + t (i+1,j  ,k,nx)   &
               +         t (i  ,j-1,k,nx) + t (i  ,j+1,k,nx) ) &
               + zedge*( t (i-1,j-1,k,nx) + t (i+1,j-1,k,nx)   &
               +         t (i-1,j+1,k,nx) + t (i+1,j+1,k,nx) )
          PQEN(i,k) = zcent * qv(i,j,k)                        &
               + zside*( qv(i-1,j  ,k) + qv(i+1,j  ,k)         &
               +         qv(i  ,j-1,k) + qv(i  ,j+1,k) )       &
               + zedge*( qv(i-1,j-1,k) + qv(i+1,j-1,k)         &
               +         qv(i-1,j+1,k) + qv(i+1,j+1,k) )
          PUEN (i,k) = zcent * u (i,j,k,nx)                    &
               + zside*( u (i-1,j  ,k,nx) + u (i+1,j  ,k,nx)   &
               +         u (i  ,j-1,k,nx) + u (i  ,j+1,k,nx) ) &
               + zedge*( u (i-1,j-1,k,nx) + u (i+1,j-1,k,nx)   &
               +         u (i-1,j+1,k,nx) + u (i+1,j+1,k,nx) )
          PVEN (i,k) = zcent * v (i,j,k,nx)                    &
               + zside*( v (i-1,j  ,k,nx) + v (i+1,j  ,k,nx)   &
               +         v (i  ,j-1,k,nx) + v (i  ,j+1,k,nx) ) &
               + zedge*( v (i-1,j-1,k,nx) + v (i+1,j-1,k,nx)   &
               +         v (i-1,j+1,k,nx) + v (i+1,j+1,k,nx) )
          PVERVEL (i,k) = zcent * omega (i,j,k)                  &
               + zside*( omega (i-1,j  ,k) + omega (i+1,j  ,k)   &
               +         omega (i  ,j-1,k) + omega (i  ,j+1,k) ) &
               + zedge*( omega (i-1,j-1,k) + omega (i+1,j-1,k)   &
               +         omega (i-1,j+1,k) + omega (i+1,j+1,k) )
        ENDDO
      ENDDO

    ENDIF

    ! variables that are not initialized:
    ! (but that are initialized in CUMASTRN or in subroutines thereof)
    ! LDCUM,KTYPE,KCBOT,KCTOP,KBOTSC,LDSC,PCAPE,
    ! PTU,PQU,PMFLXR,PMFLXS,PRAIN,PMFU,PMFD,PMFUDE_RATE,PMFDDE_RATE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    CALL CUMASTRN &
     & (  KIDIA,       KFDIA,       KLON,     KTDIA,    KLEV,&
     &    KSTEP,       KSTART,      LDLAND,   PTSPHY,&
     &    PTEN,        PQEN,        PUEN,     PVEN,     PLITOT,&
     &    PVERVEL,     PQHFL,       PAHFS,&
     &    PSSTRU,      PSSTRV,&
     &    PAP,         PAPH,        PGEO,     PGEOH,&
     &    PTENT,       PTENQ,       PTENU,    PTENV,&
     &    PTENL,       PTENI, &
     &    LDCUM,       KTYPE,       KCBOT,    KCTOP,&
     &    KBOTSC,      LDSC,&
     &    PTU,         PQU,         PLU,&
     &    PMFLXR,      PMFLXS,      PRAIN,&
     &    PMFU,        PMFD,&
     &    PMFUDE_RATE, PMFDDE_RATE, PCAPE,&
     &    KTRAC,       PCEN,        PTENC )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    tt_conv  (:,j,:) = 0.0_ireals
    qvt_conv (:,j,:) = 0.0_ireals
    qct_conv (:,j,:) = 0.0_ireals
    qit_conv (:,j,:) = 0.0_ireals
    ut_conv  (:,j,:) = 0.0_ireals
    vt_conv  (:,j,:) = 0.0_ireals

    ! note that the tendencies in the IFS scheme are incrementally updated, 
    !  whereas in COSMO the convective tendencies are added to
    !  the tendencies in the dynamics routines
    DO  i = istart, iend
      IF( LDCUM(i) ) THEN
        tt_conv  (i,j,:) = PTENT(i,:) - ttens (i,j,:)
        qvt_conv (i,j,:) = PTENQ(i,:) - qv_tens (i,j,:)
        qct_conv (i,j,:) = PTENL(i,:)   ! PTENL was initialized with 0
        qit_conv (i,j,:) = PTENI(i,:)   ! PTENI was initialized with 0
        ut_conv  (i,j,:) = PTENU(i,:) - utens(i,j,:)
        vt_conv  (i,j,:) = PTENV(i,:) - vtens(i,j,:)
      ENDIF
    ENDDO

    ! precipitation fluxes/rates:
    prr_con  (:,j) = PMFLXR (:,ke+1) * 1000.0_ireals ! m/s -> kg/m2s
    prs_con  (:,j) = PMFLXS (:,ke+1) * 1000.0_ireals ! m/s -> kg/m2s
    prne_con (:,j) = PRAIN  (:) ! already in kg/m2s

    ! convective cloud cover following Ritter, as in COSMO Tiedtke scheme:
    clc_con (:,j,:) = 0.0_ireals
    DO  i = istart, iend
      IF( LDCUM(i) .AND. KCTOP(i)>0 .AND. KCTOP(i)<(KLEV-1)) THEN
        mtop = KCTOP(i)
        mbas = KCBOT(i)
        zbas = 0.5_ireals*( hhl(i,j,mbas) + hhl(i,j,mbas+1) )
        ztop = 0.5_ireals*( hhl(i,j,mtop) + hhl(i,j,mtop+1) )
        DO  k = mtop, mbas-1
          clc_con(i,j,k) = 0.35_ireals*(ztop-zbas)/5000.0_ireals 
          IF ( k == mtop ) THEN
            zdtbot = t(i,j,k+1,nx) - t(i,j,k  ,nx)
            zdttop = t(i,j,k  ,nx) - t(i,j,k-1,nx)
            IF ( zdtbot > 0.0_ireals .AND. zdttop <= 0.0_ireals ) THEN
              clc_con(i,j,k) = 2.0_ireals*clc_con(i,j,k)  ! cloud top inversion => anvil
            ENDIF
          ENDIF
          clc_con(i,j,k) = MIN ( 1.0_ireals, MAX(0.05_ireals, clc_con(i,j,k)) )
        ENDDO
      ENDIF
    ENDDO

    ! other diagnostics:
    DO i = istart, iend
      IF(LDCUM(i)) THEN
        top_con  (i,j)   = REAL(KCTOP(i))
        bas_con  (i,j)   = REAL(KCBOT(i))
        mflx_con (i,j)   = PMFU(i,KCBOT(i)) ! updraft massflux at cloud base
      ELSE
        top_con  (i,j)   = 0.0_ireals
        bas_con  (i,j)   = 0.0_ireals
        mflx_con (i,j)   = 0.0_ireals
      ENDIF
    ENDDO
    cape_con (:,j)   = 0.0_ireals   ! use ldiagnos=.TRUE. & CAPE_ML etc. instead!
    clw_con  (:,j,:) = PLU(:,:)     ! cloud liquid water, not really used

    ! diagnostic convective fields in module 'data_fields' that are not assigned:
    ! tke_con, qcvg_con, vmax_10m, w0avg, qit_conv/qitens, qrt_conv, qst_conv, qitens

  ENDDO  ! End of loop south-north

END SUBROUTINE organize_conv_ifs

END MODULE src_conv_ifs
