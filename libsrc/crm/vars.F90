MODULE vars

USE grid

implicit none
!--------------------------------------------------------------------
! prognostic variables:

REAL(dp), POINTER, DIMENSION(:,:,:) :: u         => NULL()   ! x-wind
REAL(dp), POINTER, DIMENSION(:,:,:) :: v         => NULL()   ! y-wind
REAL(dp), POINTER, DIMENSION(:,:,:) :: w         => NULL()   ! z-wind
REAL(dp), POINTER, DIMENSION(:,:,:) :: t         => NULL()   ! moist static energy
REAL(dp), POINTER, DIMENSION(:,:,:) :: tke       => NULL()   ! sub-grid scale TKE

!--------------------------------------------------------------------
! diagnostic variables:

REAL(dp), POINTER, DIMENSION(:,:,:) :: p         => NULL()   ! pressure
REAL(dp), POINTER, DIMENSION(:,:,:) :: tabs      => NULL()   ! temperature
REAL(dp), POINTER, DIMENSION(:,:,:) :: qv        => NULL()   ! water vapour
REAL(dp), POINTER, DIMENSION(:,:,:) :: qcl       => NULL()   ! liquid water (condensate)
REAL(dp), POINTER, DIMENSION(:,:,:) :: qpl       => NULL()   ! liquid water (precipitation)
REAL(dp), POINTER, DIMENSION(:,:,:) :: qci       => NULL()   ! ice water (condensate)
REAL(dp), POINTER, DIMENSION(:,:,:) :: qpi       => NULL()   ! ice water (precipitation)
REAL(dp), POINTER, DIMENSION(:,:,:) :: tk        => NULL()   ! SGS eddy viscosity
REAL(dp), POINTER, DIMENSION(:,:,:) :: tkh       => NULL()   ! SGS eddy conductivity

!--------------------------------------------------------------------
! time-tendencies for prognostic variables

REAL(dp), POINTER, DIMENSION(:,:,:,:) :: dudt       => NULL()   
REAL(dp), POINTER, DIMENSION(:,:,:,:) :: dvdt       => NULL()   
REAL(dp), POINTER, DIMENSION(:,:,:,:) :: dwdt       => NULL()   

!----------------------------------------------------------------
! Temporary storage array:

REAL(dp), POINTER, DIMENSION(:,:,:) :: misc         => NULL()   

!------------------------------------------------------------------
! fluxes at the top and bottom of the domain:

REAL(dp), POINTER, DIMENSION(:,:) :: fluxbu       => NULL()
REAL(dp), POINTER, DIMENSION(:,:) :: fluxbv       => NULL()
REAL(dp), POINTER, DIMENSION(:,:) :: fluxbt       => NULL()
REAL(dp), POINTER, DIMENSION(:,:) :: fluxbq       => NULL()
REAL(dp), POINTER, DIMENSION(:,:) :: fluxtu       => NULL()
REAL(dp), POINTER, DIMENSION(:,:) :: fluxtv       => NULL()
REAL(dp), POINTER, DIMENSION(:,:) :: fluxtt       => NULL()
REAL(dp), POINTER, DIMENSION(:,:) :: fluxtq       => NULL()
REAL(dp), POINTER, DIMENSION(:,:) :: fzero        => NULL()
REAL(dp), POINTER, DIMENSION(:,:) :: precsfc      => NULL()  ! surface precip. rate
REAL(dp), POINTER, DIMENSION(:,:) :: precssfc     => NULL()  ! surface ice precip rate

               
!-----------------------------------------------------------------
! profiles 

REAL(dp), POINTER, DIMENSION(:) :: t0       => NULL()
REAL(dp), POINTER, DIMENSION(:) :: q0       => NULL()
REAL(dp), POINTER, DIMENSION(:) :: qv0      => NULL()
REAL(dp), POINTER, DIMENSION(:) :: tabs0    => NULL()
REAL(dp), POINTER, DIMENSION(:) :: tl0      => NULL()
REAL(dp), POINTER, DIMENSION(:) :: tv0      => NULL()
REAL(dp), POINTER, DIMENSION(:) :: u0       => NULL()
REAL(dp), POINTER, DIMENSION(:) :: v0       => NULL()
REAL(dp), POINTER, DIMENSION(:) :: tg0      => NULL()
REAL(dp), POINTER, DIMENSION(:) :: qg0      => NULL()
REAL(dp), POINTER, DIMENSION(:) :: ug0      => NULL()
REAL(dp), POINTER, DIMENSION(:) :: vg0      => NULL()
REAL(dp), POINTER, DIMENSION(:) :: p0       => NULL()
REAL(dp), POINTER, DIMENSION(:) :: tke0     => NULL()
REAL(dp), POINTER, DIMENSION(:) :: t01      => NULL()
REAL(dp), POINTER, DIMENSION(:) :: q01      => NULL()
REAL(dp), POINTER, DIMENSION(:) :: qp0      => NULL()
REAL(dp), POINTER, DIMENSION(:) :: qn0      => NULL()

!----------------------------------------------------------------
! "observed" (read from snd file) surface characteristics 

REAL(dp) :: sstobs, lhobs, shobs

!----------------------------------------------------------------
!  Domain top stuff:

REAL(dp) ::   gamt0      ! gradient of t() at the top,K/m
REAL(dp) ::   gamq0      ! gradient of q() at the top,g/g/m
REAL(dp) ::   tau_min    ! minimum damping time-scale (at the top)
REAL(dp) ::   tau_max    ! maximum damping time-scale (base of damping layer)
REAL(dp) ::   damp_depth ! damping depth as a fraction of the domain height
!-----------------------------------------------------------------
! reference vertical profiles:

REAL(dp), POINTER, DIMENSION(:) :: prespot    => NULL() ! (1000./pres)**R/cp
REAL(dp), POINTER, DIMENSION(:) :: rho        => NULL() ! air density at pressure levels, kg/m3
REAL(dp), POINTER, DIMENSION(:) :: rhow       => NULL() ! air density at vertical velocity levels, kg/m3
REAL(dp), POINTER, DIMENSION(:) :: bet        => NULL() ! ggr/tv0
REAL(dp), POINTER, DIMENSION(:) :: gamaz      => NULL() ! ggr/cp*z
REAL(dp), POINTER, DIMENSION(:) :: wsub       => NULL() ! large-scale subsidence velocity, m/s
REAL(dp), POINTER, DIMENSION(:) :: qtend      => NULL() ! large-scale tendency for total water
REAL(dp), POINTER, DIMENSION(:) :: ttend      => NULL() ! large-scale tendency for temperature
REAL(dp), POINTER, DIMENSION(:) :: utend      => NULL() ! large-scale tendency for u
REAL(dp), POINTER, DIMENSION(:) :: vtend      => NULL() ! large-scale tendency for v

!---------------------------------------------------------------------
! Large-scale and surface forcing:

integer nlsf    ! number of large-scale forcing profiles
integer nrfc    ! number of radiative forcing profiles
integer nsfc    ! number of surface forcing profiles
integer nsnd    ! number of observed soundings
integer nzlsf   ! number of large-scale forcing profiles
integer nzrfc   ! number of radiative forcing profiles
integer nzsnd   ! number of observed soundings

!!!!!!! VARIABLES NOT USED IN CAM OR SPCAM !!!!!!!
!REAL(dp), POINTER, DIMENSION(:,:) :: dqls      => NULL() ! large-scale tendency for total water
!REAL(dp), POINTER, DIMENSION(:,:) :: dtls      => NULL() ! large-scale tendency for temperature
!REAL(dp), POINTER, DIMENSION(:,:) :: ugls      => NULL() ! large-scale wind in x-direction
!REAL(dp), POINTER, DIMENSION(:,:) :: vgls      => NULL() ! large-scale wind in y-direction
!REAL(dp), POINTER, DIMENSION(:,:) :: wgls      => NULL() ! large-scale subsidence velocity (m/s)
!REAL(dp), POINTER, DIMENSION(:)   :: pres0ls   => NULL() ! surface pressure (mb)
!REAL(dp), POINTER, DIMENSION(:,:) :: zls       => NULL() ! height
!REAL(dp), POINTER, DIMENSION(:,:) :: pls       => NULL() ! pressure
!REAL(dp), POINTER, DIMENSION(:)   :: dayls     => NULL() ! large-scale forcing arrays time (days)
!REAL(dp), POINTER, DIMENSION(:,:) :: dtrfc     => NULL() ! radiative tendency for pot. temp.
!REAL(dp), POINTER, DIMENSION(:)   :: dayrfc    => NULL() ! radiative forcing arrays time (days)
!REAL(dp), POINTER, DIMENSION(:,:) :: prfc      => NULL() ! pressure/height
!REAL(dp), POINTER, DIMENSION(:)   :: sstsfc    => NULL() ! SSTs
!REAL(dp), POINTER, DIMENSION(:)   :: shsfc     => NULL() ! sensible heat flux (W/m2)
!REAL(dp), POINTER, DIMENSION(:,:) :: lhsfc     => NULL() ! latent heat flux (W/m2)
!REAL(dp), POINTER, DIMENSION(:)   :: tausfc    => NULL() ! surface drag (m2/s2)
!REAL(dp), POINTER, DIMENSION(:)   :: daysfc    => NULL() ! surface forcing arrays time (days)
!REAL(dp), POINTER, DIMENSION(:,:) :: usnd      => NULL() ! observed zonal wind
!REAL(dp), POINTER, DIMENSION(:,:) :: vsnd      => NULL() ! observed meridional wind
!REAL(dp), POINTER, DIMENSION(:,:) :: tsnd      => NULL() ! observed absolute temperature
!REAL(dp), POINTER, DIMENSION(:,:) :: qsnd      => NULL() ! observed moisture
!REAL(dp), POINTER, DIMENSION(:,:) :: zsnd      => NULL() ! height
!REAL(dp), POINTER, DIMENSION(:,:) :: psnd      => NULL() ! pressure
!REAL(dp), POINTER, DIMENSION(:)   :: daysnd    => NULL() ! number of sounding samples

!---------------------------------------------------------------------
!  Horizontally varying stuff (as a function of xy)
!

REAL(dp), POINTER, DIMENSION(:,:)   :: sstxy      => NULL() ! surface temperature xy-distribution
REAL(dp), POINTER, DIMENSION(:)     :: fcory      => NULL() ! coriolis parameter xy-distribution
REAL(dp), POINTER, DIMENSION(:)     :: fcorzy     => NULL() ! z-coriolis parameter xy-distribution
REAL(dp), POINTER, DIMENSION(:,:)   :: latitude   => NULL() ! latitude (degrees)
REAL(dp), POINTER, DIMENSION(:,:)   :: longitude  => NULL() ! longitude (degrees)
REAL(dp), POINTER, DIMENSION(:,:)   :: prec_xy    => NULL() ! mean precip. rate for output
REAL(dp), POINTER, DIMENSION(:,:)   :: shf_xy     => NULL() ! sensible heat flux for output
REAL(dp), POINTER, DIMENSION(:,:)   :: lhf_xy     => NULL() ! latent heat flux for output
REAL(dp), POINTER, DIMENSION(:,:)   :: lwns_xy    => NULL() ! mean net lw at surface
REAL(dp), POINTER, DIMENSION(:,:)   :: swns_xy    => NULL() ! mean net sw at surface
REAL(dp), POINTER, DIMENSION(:,:)   :: lwnsc_xy   => NULL() ! clear-sky mean net lw at surface
REAL(dp), POINTER, DIMENSION(:,:)   :: swnsc_xy   => NULL() ! clear-sky mean net sw at surface
REAL(dp), POINTER, DIMENSION(:,:)   :: lwnt_xy    => NULL() ! mean net lw at top of atmosphere
REAL(dp), POINTER, DIMENSION(:,:)   :: swnt_xy    => NULL() ! mean net sw at top of atmosphere
REAL(dp), POINTER, DIMENSION(:,:)   :: lwntc_xy   => NULL() ! clear-sky mean net lw at TOA
REAL(dp), POINTER, DIMENSION(:,:)   :: swntc_xy   => NULL() ! clear-sky mean net sw at TOA
REAL(dp), POINTER, DIMENSION(:,:)   :: solin_xy   => NULL() ! solar TOA insolation
REAL(dp), POINTER, DIMENSION(:,:)   :: pw_xy      => NULL() ! precipitable water
REAL(dp), POINTER, DIMENSION(:,:)   :: cw_xy      => NULL() ! cloud water path
REAL(dp), POINTER, DIMENSION(:,:)   :: iw_xy      => NULL() ! ice water path
REAL(dp), POINTER, DIMENSION(:,:)   :: u200_xy    => NULL() ! u-wind at 200 mb
REAL(dp), POINTER, DIMENSION(:,:)   :: usfc_xy    => NULL() ! u-wind at surface
REAL(dp), POINTER, DIMENSION(:,:)   :: v200_xy    => NULL() ! v-wind at 200 mb
REAL(dp), POINTER, DIMENSION(:,:)   :: vsfc_xy    => NULL() ! v-wind at surface
REAL(dp), POINTER, DIMENSION(:,:)   :: w500_xy    => NULL() ! w at 500 mb
REAL(dp), POINTER, DIMENSION(:,:)   :: qocean_xy  => NULL() ! ocean cooling in W/m2

!----------------------------------------------------------------------
!       Vertical profiles of quantities sampled for statitistics purposes:


!!$REAL     ::  twle(nz), twsb(nz), tkewle(nz), &
!!$             tkewsb(nz), precflux(nz), &
!!$             uwle(nz), uwsb(nz), vwle(nz), vwsb(nz), &
!!$             tkeleadv(nz), tkelepress(nz), tkelediss(nz), tkelediff(nz), &
!!$             tkesbbuoy(nz), tkesbshear(nz),tkesbdiss(nz), tkesbdiff(nz), &
!!$             tkelebuoy(nz), radlwup(nz), radlwdn(nz), radswup(nz), radswdn(nz), &
!!$             radqrlw(nz), radqrsw(nz), w_max, s_acld, s_acldcold, s_ar, p_conv, p_strat,&
!!$             s_acldl, s_acldm, s_acldh, s_acldisccp, &
!!$             s_acldlisccp, s_acldmisccp, s_acldhisccp, &
!!$             s_flns,s_flnt,s_flnsc,s_flntc,s_flds,s_fsns, &
!!$             s_fsnt,s_fsnsc,s_fsntc,s_fsds,s_solin, & 
!!$             t2leadv(nz),t2legrad(nz),t2lediff(nz),t2leprec(nz),t2lediss(nz), &
!!$             q2leadv(nz),q2legrad(nz),q2lediff(nz),q2leprec(nz),q2lediss(nz), &
!!$             s2leadv(nz),s2legrad(nz),s2lediff(nz),s2lediss(nz), &
!!$             twleadv(nz),twlediff(nz),twlepres(nz),twlebuoy(nz),twleprec(nz), &
!!$             qwleadv(nz),qwlediff(nz),qwlepres(nz),qwlebuoy(nz),qwleprec(nz), &
!!$             swleadv(nz),swlediff(nz),swlepres(nz),swlebuoy(nz), &
!!$             momleadv(nz,3),momlepress(nz,3),momlebuoy(nz,3), &
!!$             momlediff(nz,3),tadv(nz),tdiff(nz),tlat(nz), tlatqi(nz),qifall(nz), qpfall(nz)

REAL(dp) , DIMENSION(:), POINTER     ::  &
             twle, twsb, tkewle, &
             tkewsb, precflux, &
             uwle, uwsb, vwle, vwsb, &
             tkeleadv, tkelepress, tkelediss, tkelediff, &
             tkesbbuoy, tkesbshear,tkesbdiss, tkesbdiff, &
             tkelebuoy, radlwup, radlwdn, radswup, radswdn, &
             radqrlw, radqrsw

REAL(dp) :: &
             w_max, s_acld, s_acldcold, s_ar, p_conv, p_strat,&
             s_acldl, s_acldm, s_acldh, s_acldisccp, &
             s_acldlisccp, s_acldmisccp, s_acldhisccp, &
             s_flns,s_flnt,s_flnsc,s_flntc,s_flds,s_fsns, &
             s_fsnt,s_fsnsc,s_fsntc,s_fsds,s_solin
REAL(dp) , DIMENSION(:), POINTER     ::  &
             t2leadv,t2legrad,t2lediff,t2leprec,t2lediss, &
             q2leadv,q2legrad,q2lediff,q2leprec,q2lediss, &
             s2leadv,s2legrad,s2lediff,s2lediss, &
             twleadv,twlediff,twlepres,twlebuoy,twleprec, &
             qwleadv,qwlediff,qwlepres,qwlebuoy,qwleprec, &
             swleadv,swlediff,swlepres,swlebuoy
REAL(dp) , DIMENSION(:,:), POINTER   ::  &
             momleadv,momlepress,momlebuoy, &
             momlediff
REAL(dp) , DIMENSION(:), POINTER     ::  &
             tadv,tdiff,tlat, tlatqi,qifall, qpfall


! register functions:

real, external :: esatw_crm,esati_crm,dtesatw_crm,dtesati_crm
real, external :: qsatw_crm,qsati_crm,dtqsatw_crm,dtqsati_crm
integer, external :: lenstr

! energy conservation diagnostics:
 
real(dp) total_water_before, total_water_after
real(dp) total_water_evap, total_water_prec, total_water_ls

!===========================================================================
! UW ADDITIONS

! conditional average statistics, subsumes cloud_factor, core_factor, coredn_factor
integer :: ncondavg, icondavg_cld, icondavg_cor, icondavg_cordn, &
           icondavg_satdn, icondavg_satup, icondavg_env

!!!!!! FOR USE IN MORRISON2005 MICROPHYSICS  !!!!!!   
!real, allocatable :: condavg_factor(:,:)             ! replaces cloud_factor, core_factor
!real, allocatable :: condavg_mask(:,:,:,:)           ! indicator array for various conditional averages
!character(LEN=8),  allocatable :: condavgname(:)     ! array of short names
!character(LEN=25), allocatable :: condavglongname(:) ! array of long names

REAL(dp), POINTER, DIMENSION(:)   :: qlsvadv  => NULL() ! Large-scale vertical advection tendency for total water
REAL(dp), POINTER, DIMENSION(:)   :: tlsvadv  => NULL() ! Large-scale vertical advection tendency for temperature
REAL(dp), POINTER, DIMENSION(:)   :: ulsvadv  => NULL() ! Large-scale vertical advection tendency for zonal wind velocity
REAL(dp), POINTER, DIMENSION(:)   :: vlsvadv  => NULL() ! Large-scale vertical advection tendency for meridional wind velocity

REAL(dp), POINTER, DIMENSION(:)   :: qnudge   => NULL() ! Nudging of horiz.-averaged total water profile
REAL(dp), POINTER, DIMENSION(:)   :: tnudge   => NULL() ! Nudging of horiz.-averaged temperature profile
REAL(dp), POINTER, DIMENSION(:)   :: unudge   => NULL() ! Nudging of horiz.-averaged zonal wind velocity
REAL(dp), POINTER, DIMENSION(:)   :: vnudge   => NULL() ! Nudging of horiz.-averaged meridional wind velocity

REAL(dp), POINTER, DIMENSION(:)   :: qstor    => NULL() ! Storage of horiz.-averaged total water profile
REAL(dp), POINTER, DIMENSION(:)   :: tstor    => NULL() ! Storage of horiz.-averaged temperature profile
REAL(dp), POINTER, DIMENSION(:)   :: ustor    => NULL() ! Storage of horiz.-averaged zonal wind profile
REAL(dp), POINTER, DIMENSION(:)   :: vstor    => NULL() ! Storage of horiz.-averaged meridional wind profile

REAL(dp), POINTER, DIMENSION(:)   :: utendcor => NULL() ! coriolis acceleration of zonal wind velocity
REAL(dp), POINTER, DIMENSION(:)   :: vtendcor => NULL() ! coriolis acceleration of meridional wind velocity

! 850 mbar horizontal winds
REAL(dp), POINTER, DIMENSION(:,:)   :: u850_xy  => NULL() ! zonal wind velocity at 850 mb
REAL(dp), POINTER, DIMENSION(:,:)   :: v850_xy  => NULL() ! meridional wind velocity at 850 mb

! Surface pressure
REAL(dp), POINTER, DIMENSION(:,:)   :: psfc_xy  => NULL() ! pressure (mb) at lowest grid point

! Saturated water vapor path, useful for computing column relative humidity
REAL(dp), POINTER, DIMENSION(:,:)   :: swvp_xy  => NULL() ! saturated water vapor path (wrt water)

! Cloud and echo top heights, and cloud top temperature (instantaneous)
REAL(dp), POINTER, DIMENSION(:,:)   :: cloudtopheight =>NULL()
REAL(dp), POINTER, DIMENSION(:,:)   :: echotopheight  =>NULL()
REAL(dp), POINTER, DIMENSION(:,:)   :: cloudtoptemp   =>NULL()

! END UW ADDITIONS
!===========================================================================
! Initial bubble parameters. Activated when perturb_type = 2
  REAL(dp) bubble_x0  
  REAL(dp) bubble_y0 
  REAL(dp) bubble_z0 
  REAL(dp) bubble_radius_hor 
  REAL(dp) bubble_radius_ver 
  REAL(dp) bubble_dtemp 
  REAL(dp) bubble_dq 


END MODULE vars

! ====================================================================

  SUBROUTINE crm_allocate_vars

    use vars

    IMPLICIT NONE
    
    ALLOCATE(z(nz)); z=0.
    ALLOCATE(pres(nzm)); pres=0.
    ALLOCATE(zi(nz)); zi=0.
    ALLOCATE(presi(nz)); presi=0.
    ALLOCATE(adz(nzm)); adz=0.
    ALLOCATE(adzw(nz)); adzw=0.
    ALLOCATE(grdf_x(nzm)); grdf_x=0.
    ALLOCATE(grdf_y(nzm)); grdf_y=0.
    ALLOCATE(grdf_z(nzm)); grdf_z=0.

    ALLOCATE(u(dimx1_u:dimx2_u,dimy1_u:dimy2_u,nzm)); u=0.
    ALLOCATE(v(dimx1_v:dimx2_v,dimy1_v:dimy2_v,nzm)); v=0.
    ALLOCATE(w(dimx1_w:dimx2_w,dimy1_w:dimy2_w,nz )); w=0.
    ALLOCATE(t(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm)); t=0.
    ALLOCATE(tke(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm)); tke=0.

    ALLOCATE(p(0:nx,(1-YES3D):ny,nzm)); p=0.
    ALLOCATE(tabs(nx,ny,nzm)); tabs=0.
    ALLOCATE(qv(nx,ny,nzm)); qv=0.
    ALLOCATE(qcl(nx,ny,nzm)); qcl=0.
    ALLOCATE(qpl(nx,ny,nzm)); qpl=0.
    ALLOCATE(qci(nx,ny,nzm)); qci=0.
    ALLOCATE(qpi(nx,ny,nzm)); qpi=0.
    ALLOCATE(tk(0:nxp1,(1-YES3D):nyp1,nzm)); tk   = 0.
    ALLOCATE(tkh(0:nxp1,(1-YES3D):nyp1,nzm)); tkh = 0.

    ALLOCATE(dudt(nxp1,ny,nzm,3)); dudt=0.
    ALLOCATE(dvdt(nx,nyp1,nzm,3)); dvdt=0.
    ALLOCATE(dwdt(nx,ny,nz,3));    dwdt=0.

    ALLOCATE(misc(nx,ny,nz));      misc=0.

    ALLOCATE(fluxbu(nx,ny));fluxbu= 0.
    ALLOCATE(fluxbv(nx,ny));fluxbv= 0.
    ALLOCATE(fluxbt(nx,ny)); fluxbt= 0.
    ALLOCATE(fluxbq(nx,ny)); fluxbq= 0.
    ALLOCATE(fluxtu(nx,ny)); fluxtu= 0.
    ALLOCATE(fluxtv(nx,ny)); fluxtv= 0.
    ALLOCATE(fluxtt(nx,ny)); fluxtt= 0.
    ALLOCATE(fluxtq(nx,ny)); fluxtq= 0.
    ALLOCATE(fzero(nx,ny)) ; fzero= 0.
    ALLOCATE(precsfc(nx,ny));precsfc= 0.
    ALLOCATE(precssfc(nx,ny));precssfc= 0.

    ALLOCATE(t0(nzm));    t0= 0.
    ALLOCATE(q0(nzm));    q0= 0.
    ALLOCATE(qv0(nzm));   qv0= 0.
    ALLOCATE(tabs0(nzm)); tabs0= 0.
    ALLOCATE(tl0(nzm));   tl0= 0.
    ALLOCATE(tv0(nzm));   tv0= 0.
    ALLOCATE(u0(nzm));    u0= 0.
    ALLOCATE(v0(nzm));    v0= 0.
    ALLOCATE(tg0(nzm));   tg0= 0.
    ALLOCATE(qg0(nzm));   qg0= 0.
    ALLOCATE(ug0(nzm));   ug0= 0.
    ALLOCATE(vg0(nzm));   vg0= 0.
    ALLOCATE(p0(nzm));    p0= 0.
    ALLOCATE(tke0(nzm));  tke0= 0.
    ALLOCATE(t01(nzm));   t01= 0.
    ALLOCATE(q01(nzm));   q01= 0.
    ALLOCATE(qp0(nzm));   qp0= 0.
    ALLOCATE(qn0(nzm));   qn0= 0.

    ALLOCATE(prespot(nzm)); prespot= 0.
    ALLOCATE(rho(nzm));     rho= 0.
    ALLOCATE(rhow(nz));     rhow= 0.
    ALLOCATE(bet(nzm));     bet= 0.
    ALLOCATE(gamaz(nzm));   gamaz= 0.
    ALLOCATE(wsub(nz));     wsub= 0.
    ALLOCATE(qtend(nzm));   qtend= 0.
    ALLOCATE(ttend(nzm));   ttend= 0.
    ALLOCATE(utend(nzm));   utend= 0.
    ALLOCATE(vtend(nzm));   vtend= 0.

    ALLOCATE(sstxy(0:nx,(1-YES3D):ny)); sstxy = 0.
    ALLOCATE(fcory(ny)); fcory = 0.
    ALLOCATE(fcorzy(ny)); fcorzy = 0.
    ALLOCATE(latitude(nx,ny)); latitude = 0.
    ALLOCATE(longitude(nx,ny)); longitude = 0.
    ALLOCATE(prec_xy(nx,ny)); prec_xy = 0.
    ALLOCATE(shf_xy(nx,ny)); shf_xy = 0.
    ALLOCATE(lhf_xy(nx,ny)); lhf_xy = 0.
    ALLOCATE(lwns_xy(nx,ny)); lwns_xy = 0.
    ALLOCATE(swns_xy(nx,ny)); swns_xy = 0.
    ALLOCATE(lwnsc_xy(nx,ny)); lwnsc_xy = 0.
    ALLOCATE(swnsc_xy(nx,ny)); swnsc_xy = 0.
    ALLOCATE(lwnt_xy(nx,ny)); lwnt_xy = 0.
    ALLOCATE(swnt_xy(nx,ny)); swnt_xy = 0.
    ALLOCATE(lwntc_xy(nx,ny)); lwntc_xy = 0.
    ALLOCATE(swntc_xy(nx,ny)); swntc_xy = 0.
    ALLOCATE(solin_xy(nx,ny)); solin_xy = 0.
    ALLOCATE(pw_xy(nx,ny)); pw_xy = 0.
    ALLOCATE(cw_xy(nx,ny)); cw_xy = 0.
    ALLOCATE(iw_xy(nx,ny)); iw_xy = 0.
    ALLOCATE(u200_xy(nx,ny)); u200_xy = 0.
    ALLOCATE(usfc_xy(nx,ny)); usfc_xy = 0.
    ALLOCATE(v200_xy(nx,ny)); v200_xy = 0.
    ALLOCATE(vsfc_xy(nx,ny)); vsfc_xy = 0.
    ALLOCATE(w500_xy(nx,ny)); w500_xy = 0.
    ALLOCATE(qocean_xy(nx,ny)); qocean_xy = 0.

    ALLOCATE(qlsvadv(nzm));qlsvadv= 0.
    ALLOCATE(tlsvadv(nzm));tlsvadv= 0.
    ALLOCATE(ulsvadv(nzm));ulsvadv= 0.
    ALLOCATE(vlsvadv(nzm));vlsvadv= 0.

    ALLOCATE(qnudge(nzm));qnudge= 0.
    ALLOCATE(tnudge(nzm));tnudge= 0.
    ALLOCATE(unudge(nzm));unudge= 0.
    ALLOCATE(vnudge(nzm));vnudge= 0.

    ALLOCATE(qstor(nzm));qstor= 0.
    ALLOCATE(tstor(nzm));tstor= 0.
    ALLOCATE(ustor(nzm));ustor= 0.
    ALLOCATE(vstor(nzm));vstor= 0.

    ALLOCATE(utendcor(nzm)); utendcor = 0.
    ALLOCATE(vtendcor(nzm)); vtendcor = 0.

    ALLOCATE(u850_xy(nx,ny));u850_xy= 0.
    ALLOCATE(v850_xy(nx,ny));v850_xy= 0.
    ALLOCATE(psfc_xy(nx,ny));psfc_xy= 0.
    ALLOCATE(swvp_xy(nx,ny));swvp_xy= 0.

    ALLOCATE(cloudtopheight(nx,ny)); cloudtopheight = 0.
    ALLOCATE(echotopheight(nx,ny));  echotopheight  = 0.
    ALLOCATE(cloudtoptemp(nx,ny));   cloudtoptemp   = 0.

    ALLOCATE(twle(nz));twle = 0.
    ALLOCATE(twsb(nz));twsb  = 0.
    ALLOCATE(tkewle(nz));tkewle  = 0.
    ALLOCATE(tkewsb(nz));tkewsb  = 0.
    ALLOCATE(precflux(nz));precflux = 0.
    ALLOCATE(uwle(nz));uwle      = 0.
    ALLOCATE(uwsb(nz));uwsb = 0.
    ALLOCATE(vwle(nz));vwle = 0.
    ALLOCATE(vwsb(nz));vwsb = 0.
    ALLOCATE(tkeleadv(nz));tkeleadv = 0.
    ALLOCATE(tkelepress(nz));tkelepress = 0.
    ALLOCATE(tkelediss(nz));tkelediss = 0.
    ALLOCATE(tkelediff(nz));tkelediff = 0.
    ALLOCATE(tkesbbuoy(nz));tkesbbuoy = 0.
    ALLOCATE(tkesbshear(nz));tkesbshear = 0.
    ALLOCATE(tkesbdiss(nz));tkesbdiss = 0.
    ALLOCATE(tkesbdiff(nz));tkesbdiff = 0.
    ALLOCATE(tkelebuoy(nz));tkelebuoy = 0.
    ALLOCATE(radlwup(nz));radlwup = 0.
    ALLOCATE(radlwdn(nz));radlwdn = 0.
    ALLOCATE(radswup(nz));radswup = 0.
    ALLOCATE(radswdn(nz));radswdn  = 0.
    ALLOCATE(radqrlw(nz));radqrlw = 0.
    ALLOCATE(radqrsw(nz));radqrsw = 0.
    ALLOCATE(t2leadv(nz));t2leadv = 0.
    ALLOCATE(t2legrad(nz));t2legrad = 0.
    ALLOCATE(t2lediff(nz));t2lediff = 0.
    ALLOCATE(t2leprec(nz));t2leprec = 0.
    ALLOCATE(t2lediss(nz));t2lediss = 0.
    ALLOCATE(q2leadv(nz));q2leadv = 0.
    ALLOCATE(q2legrad(nz));q2legrad = 0.
    ALLOCATE(q2lediff(nz));q2lediff = 0.
    ALLOCATE(q2leprec(nz));q2leprec = 0.
    ALLOCATE(q2lediss(nz));q2lediss = 0.
    ALLOCATE(s2leadv(nz));s2leadv = 0.
    ALLOCATE(s2legrad(nz));s2legrad = 0.
    ALLOCATE(s2lediff(nz));s2lediff = 0.
    ALLOCATE(s2lediss(nz));s2lediss = 0.
    ALLOCATE(twleadv(nz));twleadv = 0.
    ALLOCATE(twlediff(nz));twlediff = 0.
    ALLOCATE(twlepres(nz));twlepres = 0.
    ALLOCATE(twlebuoy(nz));twlebuoy = 0.
    ALLOCATE(twleprec(nz));twleprec = 0.
    ALLOCATE(qwleadv(nz));qwleadv = 0.
    ALLOCATE(qwlediff(nz));qwlediff = 0.
    ALLOCATE(qwlepres(nz));qwlepres = 0.
    ALLOCATE(qwlebuoy(nz));qwlebuoy = 0.
    ALLOCATE(qwleprec(nz));qwleprec = 0.
    ALLOCATE(swleadv(nz));swleadv = 0.
    ALLOCATE(swlediff(nz));swlediff = 0.
    ALLOCATE(swlepres(nz));swlepres = 0.
    ALLOCATE(swlebuoy(nz));swlebuoy = 0.
    ALLOCATE(momleadv(nz,3));momleadv = 0.
    ALLOCATE(momlepress(nz,3));momlepress = 0.
    ALLOCATE(momlebuoy(nz,3));momlebuoy = 0.
    ALLOCATE(momlediff(nz,3));momlediff = 0.
    ALLOCATE(tadv(nz));tadv = 0.
    ALLOCATE(tdiff(nz));tdiff = 0.
    ALLOCATE(tlat(nz));tlat = 0.
    ALLOCATE(tlatqi(nz));tlatqi = 0.
    ALLOCATE(qifall(nz));qifall = 0.
    ALLOCATE(qpfall(nz));qpfall = 0.


  END SUBROUTINE crm_allocate_vars

! ====================================================================


  SUBROUTINE crm_deallocate_vars

    use vars 

    IMPLICIT NONE
    
    DEALLOCATE(z)
    DEALLOCATE(pres)
    DEALLOCATE(zi)
    DEALLOCATE(presi)
    DEALLOCATE(adz)
    DEALLOCATE(adzw)
    DEALLOCATE(grdf_x)
    DEALLOCATE(grdf_y)
    DEALLOCATE(grdf_z)

    DEALLOCATE(u)
    DEALLOCATE(v)
    DEALLOCATE(w)
    DEALLOCATE(t)
    DEALLOCATE(tke)

    DEALLOCATE(p)
    DEALLOCATE(tabs)
    DEALLOCATE(qv)
    DEALLOCATE(qcl)
    DEALLOCATE(qpl)
    DEALLOCATE(qci)
    DEALLOCATE(qpi)
    DEALLOCATE(tk)
    DEALLOCATE(tkh)

    DEALLOCATE(dudt)
    DEALLOCATE(dvdt)
    DEALLOCATE(dwdt)

    DEALLOCATE(misc)

    DEALLOCATE(fluxbu)
    DEALLOCATE(fluxbv)
    DEALLOCATE(fluxbt)
    DEALLOCATE(fluxbq)
    DEALLOCATE(fluxtu)
    DEALLOCATE(fluxtv)
    DEALLOCATE(fluxtt)
    DEALLOCATE(fluxtq)
    DEALLOCATE(fzero)
    DEALLOCATE(precsfc)
    DEALLOCATE(precssfc)

    DEALLOCATE(t0)
    DEALLOCATE(q0)
    DEALLOCATE(qv0)
    DEALLOCATE(tabs0)
    DEALLOCATE(tl0)
    DEALLOCATE(tv0)
    DEALLOCATE(u0)
    DEALLOCATE(v0)
    DEALLOCATE(tg0)
    DEALLOCATE(qg0)
    DEALLOCATE(ug0)
    DEALLOCATE(vg0)
    DEALLOCATE(p0)
    DEALLOCATE(tke0)
    DEALLOCATE(t01)
    DEALLOCATE(q01)
    DEALLOCATE(qp0)
    DEALLOCATE(qn0)

    DEALLOCATE(prespot)
    DEALLOCATE(rho)
    DEALLOCATE(rhow)
    DEALLOCATE(bet)
    DEALLOCATE(gamaz)
    DEALLOCATE(wsub)
    DEALLOCATE(qtend)
    DEALLOCATE(ttend)
    DEALLOCATE(utend)
    DEALLOCATE(vtend)

    DEALLOCATE(sstxy)
    DEALLOCATE(fcory)
    DEALLOCATE(fcorzy)
    DEALLOCATE(latitude)
    DEALLOCATE(longitude)
    DEALLOCATE(prec_xy)
    DEALLOCATE(shf_xy)
    DEALLOCATE(lhf_xy)
    DEALLOCATE(lwns_xy)
    DEALLOCATE(swns_xy)
    DEALLOCATE(lwnsc_xy)
    DEALLOCATE(swnsc_xy)
    DEALLOCATE(lwnt_xy)
    DEALLOCATE(swnt_xy)
    DEALLOCATE(lwntc_xy)
    DEALLOCATE(swntc_xy)
    DEALLOCATE(solin_xy)
    DEALLOCATE(pw_xy)
    DEALLOCATE(cw_xy)
    DEALLOCATE(iw_xy)
    DEALLOCATE(u200_xy)
    DEALLOCATE(usfc_xy)
    DEALLOCATE(v200_xy)
    DEALLOCATE(vsfc_xy)
    DEALLOCATE(w500_xy)
    DEALLOCATE(qocean_xy)

    DEALLOCATE(qlsvadv)
    DEALLOCATE(tlsvadv)
    DEALLOCATE(ulsvadv)
    DEALLOCATE(vlsvadv)

    DEALLOCATE(qnudge)
    DEALLOCATE(tnudge)
    DEALLOCATE(unudge)
    DEALLOCATE(vnudge)

    DEALLOCATE(qstor)
    DEALLOCATE(tstor)
    DEALLOCATE(ustor)
    DEALLOCATE(vstor)

    DEALLOCATE(utendcor)
    DEALLOCATE(vtendcor)

    DEALLOCATE(u850_xy)
    DEALLOCATE(v850_xy)
    DEALLOCATE(psfc_xy)
    DEALLOCATE(swvp_xy)

    DEALLOCATE(cloudtopheight)
    DEALLOCATE(echotopheight)
    DEALLOCATE(cloudtoptemp)

    DEALLOCATE(twle)
    DEALLOCATE(twsb) 
    DEALLOCATE(tkewle) 
    DEALLOCATE(tkewsb) 
    DEALLOCATE(precflux)
    DEALLOCATE(uwle)     
    DEALLOCATE(uwsb)
    DEALLOCATE(vwle)
    DEALLOCATE(vwsb)
    DEALLOCATE(tkeleadv)
    DEALLOCATE(tkelepress)
    DEALLOCATE(tkelediss)
    DEALLOCATE(tkelediff)
    DEALLOCATE(tkesbbuoy)
    DEALLOCATE(tkesbshear)
    DEALLOCATE(tkesbdiss)
    DEALLOCATE(tkesbdiff)
    DEALLOCATE(tkelebuoy)
    DEALLOCATE(radlwup)
    DEALLOCATE(radlwdn)
    DEALLOCATE(radswup)
    DEALLOCATE(radswdn) 
    DEALLOCATE(radqrlw)
    DEALLOCATE(radqrsw)
    DEALLOCATE(t2leadv)
    DEALLOCATE(t2legrad)
    DEALLOCATE(t2lediff)
    DEALLOCATE(t2leprec)
    DEALLOCATE(t2lediss)
    DEALLOCATE(q2leadv)
    DEALLOCATE(q2legrad)
    DEALLOCATE(q2lediff)
    DEALLOCATE(q2leprec)
    DEALLOCATE(q2lediss)
    DEALLOCATE(s2leadv)
    DEALLOCATE(s2legrad)
    DEALLOCATE(s2lediff)
    DEALLOCATE(s2lediss)
    DEALLOCATE(twleadv)
    DEALLOCATE(twlediff)
    DEALLOCATE(twlepres)
    DEALLOCATE(twlebuoy)
    DEALLOCATE(twleprec)
    DEALLOCATE(qwleadv)
    DEALLOCATE(qwlediff)
    DEALLOCATE(qwlepres)
    DEALLOCATE(qwlebuoy)
    DEALLOCATE(qwleprec)
    DEALLOCATE(swleadv)
    DEALLOCATE(swlediff)
    DEALLOCATE(swlepres)
    DEALLOCATE(swlebuoy)
    DEALLOCATE(momleadv)
    DEALLOCATE(momlepress)
    DEALLOCATE(momlebuoy)
    DEALLOCATE(momlediff)
    DEALLOCATE(tadv)
    DEALLOCATE(tdiff)
    DEALLOCATE(tlat)
    DEALLOCATE(tlatqi)
    DEALLOCATE(qifall)
    DEALLOCATE(qpfall)


  END SUBROUTINE crm_deallocate_vars

! ====================================================================
