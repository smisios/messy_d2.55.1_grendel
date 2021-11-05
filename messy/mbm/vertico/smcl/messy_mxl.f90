
MODULE messy_mxl

  !-------------------------------------------------------------------------------------------------------
  !  MXL: MiXed Layer model for the diurnal dynamics of the convective boundary layer  
  !
  !  AUTHORS:  Ruud Janssen, Andrea Pozzer, MPIC, Sept 2013
  !-------------------------------------------------------------------------------------------------------


  ! MESSy
  USE messy_mxl_mem
  USE messy_main_constants_mem,   ONLY: dp, pi, TWOPI, PI_2, OneDay, DTR, &
                                        c_vKar, cp_air, rho_air, rho_H2O, &
                                        alv, g, rd, rv, cpv,              &
                                        stbo, solc, OneDay, vtmpc1 
  USE messy_main_grid_def_mem_bi, ONLY : nlon, nlat, nlev 
  USE messy_main_timer,           ONLY : current_date, start_date, stop_date

  IMPLICIT NONE 

  PRIVATE

  PUBLIC :: modstr
  PUBLIC :: modver

  PUBLIC :: mxl_read_nml_ctrl
  PUBLIC :: mxl_algorithm, mxl_radiation, mxl_surfacelayer, mxl_momentumflux, mxl_landsurface
  PUBLIC :: emis_simple

  REAL(dp), PARAMETER :: g2kg = 1.e-3_dp      ! convert g/kg to g/g

CONTAINS

  ! ---------------------------------------------------------------------------
SUBROUTINE mxl_algorithm(hbl, hsl, mpress_3d, mpressi_3d, wthetasmax, wqsmax, beta, omega, &
                         thetam, thetate, dtheta, gammatheta, lgamma, gammath2, hcrit, advthetamax, &
                         qm, qte,    dq,     gammaq,     advqmax, rh, &
                         um,     du,     gammau, uws, &
                         vm,     dv,     gammav, vws, &      
                         we, ws, thetav, dthetav, wthetas, wqs, wthetavs, wthetave,wqe, betaq, SH, LE, ueff, &
                         f_wthetas, starttime_wths,stoptime_wths,&
                         f_wqs, starttime_wqs, stoptime_wqs, &
                         starttime_adv, stoptime_adv)

  USE messy_main_timer,      ONLY: delta_time, lstart 
  USE messy_main_blather_bi, ONLY: error_bi
   
  implicit none

  ! mixed layer variables & parameters
  REAL(dp), DIMENSION(:,:),   INTENT(INOUT) :: hbl
  REAL(dp), DIMENSION(:,:),   INTENT(out)   :: hsl
  REAL(dp), DIMENSION(:,:,:), INTENT(IN)    :: mpress_3d
  REAL(dp), DIMENSION(:,:,:), INTENT(IN)    :: mpressi_3d
  REAL(dp), DIMENSION(:,:)  , INTENT(IN)    :: wthetasmax  
  REAL(dp), DIMENSION(:,:)  , INTENT(IN)    :: wqsmax  
  REAL(dp)                  , INTENT(IN)    :: beta
  REAL(dp)                  , INTENT(IN)    :: omega
  REAL(dp), DIMENSION(:,:,:), INTENT(INOUT) :: thetam
  REAL(dp), DIMENSION(:,:,:), INTENT(OUT)   :: thetate
  REAL(dp), DIMENSION(:,:)  , INTENT(INOUT) :: dtheta
  REAL(dp)                  , INTENT(INOUT) :: gammatheta
  REAL(dp)                  , INTENT(IN)    :: gammath2
  REAL(dp), DIMENSION(:,:)  , INTENT(IN)    :: hcrit
  LOGICAL                   , INTENT(IN)    :: lgamma
  REAL(dp), DIMENSION(nlon,nlev)  , INTENT(IN)    :: advthetamax
  REAL(dp), DIMENSION(:,:,:), INTENT(INOUT) :: qm
  REAL(dp), DIMENSION(:,:,:), INTENT(OUT)   :: qte
  REAL(dp), DIMENSION(:,:)  , INTENT(INOUT) :: dq
  REAL(dp)                  , INTENT(IN)    :: gammaq
  REAL(dp), DIMENSION(nlon,nlev)  , INTENT(IN)    :: advqmax
  REAL(dp), DIMENSION(:,:,:), INTENT(OUT)   :: rh
  REAL(dp), DIMENSION(:,:,:), INTENT(INOUT) :: um
  REAL(dp), DIMENSION(:,:)  , INTENT(INOUT) :: du
  REAL(dp)                  , INTENT(IN)    :: gammau
  REAL(dp), DIMENSION(:,:)  , INTENT(IN)    :: uws
  REAL(dp), DIMENSION(:,:,:), INTENT(INOUT) :: vm
  REAL(dp), DIMENSION(:,:)  , INTENT(INOUT) :: dv
  REAL(dp)                  , INTENT(IN)    :: gammav
  REAL(dp), DIMENSION(:,:)  , INTENT(IN)    :: vws
  REAL(dp), DIMENSION(:,:)  , INTENT(out)   :: ueff              ! effective wind speed (m s-1)

  REAL(dp), DIMENSION(:,:)  , INTENT(OUT) :: we         ! entrainment velocity (m s-1)
  REAL(dp), DIMENSION(:,:)  , INTENT(OUT) :: ws         ! subsidence velocity (m s-1)
  REAL(dp), DIMENSION(:,:,:), INTENT(OUT) :: thetav     ! mixed layer virtual potential temp. (K)
  REAL(dp), DIMENSION(:,:)  , INTENT(OUT) :: dthetav    ! virtual potential temperature jump (K)
  REAL(dp), DIMENSION(:,:)  , INTENT(OUT) :: wthetas    ! surface heat flux (K m s-1)
  REAL(dp), DIMENSION(:,:)  , INTENT(OUT) :: wqs        ! surface moisture flux (g kg-1 m s-1) 
  REAL(dp), DIMENSION(:,:)  , INTENT(OUT) :: wthetavs   ! surface buoyancy flux (K m s-1)
  REAL(dp), DIMENSION(:,:)  , INTENT(OUT) :: wqe        ! entrainment moisture flux (g kg-1 m s-1)
  REAL(dp), DIMENSION(:,:)  , INTENT(OUT) :: wthetave   ! entrainment buoyancy flux (K m s-1)
  REAL(dp), DIMENSION(:,:)  , INTENT(OUT) :: betaq      ! ratio of entrainment/surface moisture flux (-) 
  REAL(dp), DIMENSION(:,:)  , INTENT(INOUT) :: SH         ! sensible heat flux (W m-2)
  REAL(dp), DIMENSION(:,:)  , INTENT(INOUT) :: LE         ! latent heat flux (W m-2)

  REAL(dp), DIMENSION(nlon,nlat)      :: uwe  ! entrainment momentum flux x-direction (m2 s-2) 
  REAL(dp), DIMENSION(nlon,nlat)      :: vwe  ! entrainment momentum flux y-direction (m2 s-2)
  REAL(dp), DIMENSION(nlon,nlat)      :: hbl_p1
  REAL(dp), DIMENSION(nlon,nlat,nlev) :: thetam_p1
  REAL(dp), DIMENSION(nlon,nlat)      :: dtheta_p1
  REAL(dp), DIMENSION(nlon,nlat,nlev) :: qm_p1
  REAL(dp), DIMENSION(nlon,nlat)      :: dq_p1
  REAL(dp), DIMENSION(nlon,nlev)      :: advq
  REAL(dp), DIMENSION(nlon,nlev)      :: advtheta
  REAL(dp), DIMENSION(nlon,nlat,nlev) :: um_p1
  REAL(dp), DIMENSION(nlon,nlat)      :: du_p1
  REAL(dp), DIMENSION(nlon,nlat,nlev) :: vm_p1
  REAL(dp), DIMENSION(nlon,nlat)      :: dv_p1
  REAL(dp), DIMENSION(nlon,nlat,nlev) :: esat
  REAL(dp), DIMENSION(nlon,nlat,nlev) :: e

  INTEGER  :: i, j
  CHARACTER(LEN=STRLEN_MEDIUM)  :: f_wthetas
  CHARACTER(LEN=STRLEN_MEDIUM)  :: f_wqs 
  REAL(dp), INTENT(in) :: starttime_wths
  REAL(dp), INTENT(in) :: stoptime_wths
  REAL(dp), INTENT(in) :: starttime_wqs 
  REAL(dp), INTENT(in) :: stoptime_wqs 
  REAL(dp), INTENT(in) :: starttime_adv 
  REAL(dp), INTENT(in) :: stoptime_adv
      
  REAL(dp)                           :: fc   ! Coriolis parameter [s-1], see Stull 1988, pp. 78
  REAL(dp), DIMENSION(nlon,nlat)     :: wstar             ! free convection scaling velocity (m s-1)

  hsl = 0.1_dp * hbl

  thetav(1:nlon,1:nlat,nlev) = thetam(1:nlon,1:nlat,nlev) * (1. + vtmpc1 * qm(1:nlon,1:nlat,nlev) * g2kg) 
  wstar(1:nlon,1:nlat)       = ( g / thetav(1:nlon,1:nlat,nlev) * hbl(1:nlon,1:nlat) * &
                                wthetavs(1:nlon,1:nlat) )  ** (1._dp / 3._dp)  
  ueff(1:nlon,1:nlat)        = sqrt(um(1:nlon,1:nlat,nlev)**2._dp + vm(1:nlon,1:nlat,nlev)**2._dp &
                               + wstar(1:nlon,1:nlat)**2._dp)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! calculate surface heat fluxes (sensible and latent)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  select case(f_wthetas)
    case('NOFLUX') ! no heat flux -> no BL growth -> box model with hbl = hbl_ic
      wthetas = 0.
    case('CONST') ! constant heat flux with start & end time
      if ((current_date%second .le. start_date%second+starttime_wths) .or. &
          (current_date%second .ge. stop_date%second+stoptime_wths)) then
        wthetas = 0. 
      else
        wthetas = wthetasmax 
      end if
    case('SINE') ! sinusoidal heat flux
      if ((current_date%second .le. start_date%second + starttime_wths) .or. &
          (current_date%second .ge. stop_date%second + stoptime_wths)) then
        wthetas = 0.
      else
        wthetas = wthetasmax * sin(pi * (current_date%second - start_date%second - starttime_wths)/ &
                                        (stop_date%second + stoptime_wths - start_date%second - starttime_wths)) 
      end if
    case('COSINE') ! cosine heat flux
      if ((current_date%second .le. start_date%second + starttime_wths) .or. &
          (current_date%second .ge. stop_date%second + stoptime_wths)) then
        wthetas = 0.
      else
        wthetas = (wthetasmax/2._dp)*(1._dp - cos(2._dp*pi*(current_date%second - start_date%second &
                - starttime_wths)/(stop_date%second + stoptime_wths - start_date%second - starttime_wths))) 
      end if
    case('INTERACT') ! heat flux calculated interactively
      wthetas(1:nlon,1:nlat)= SH(1:nlon,1:nlat) / (rho_air * cp_air)
    case('IMPORT') ! import heat flux from measurements
      write(*,*) 'import heat flux from measurements: not yet available'   
      stop
    case('SEA') ! heat fluxes over sea surface
      !wthetas = wthetas_sea
    case default 
      write(*,*) 'Flag for the function of wthetas is invalid', f_wthetas
      stop 
  end select

  select case(f_wqs)
    case('NOFLUX') ! no moisture flux -> no BL growth -> box model with hbl = hbl0
      wqs = 0._dp
    case('CONST') ! constant moisture flux with start & end time
      if ((current_date%second .le. start_date%second+starttime_wqs) .or. &
          (current_date%second .ge. stop_date%second+stoptime_wqs)) then
        wqs = 0._dp 
      else
        wqs = wqsmax 
      end if
    case('SINE') ! sinusoidal moisture flux
      if ((current_date%second .le. start_date%second + starttime_wqs) .or. &
          (current_date%second .ge. stop_date%second + stoptime_wqs)) then
        wqs = 0._dp
      else
        wqs = wqsmax * sin(pi * (current_date%second - start_date%second - starttime_wqs)/ &
                                        (stop_date%second + stoptime_wqs - start_date%second - starttime_wqs)) 
      end if
    case('COSINE') ! cosine moisture flux
      if ((current_date%second .le. start_date%second + starttime_wqs) .or. &
          (current_date%second .ge. stop_date%second + stoptime_wqs)) then
        wqs = 0._dp
      else
        wqs = (wqsmax/2._dp)*(1._dp - cos(2._dp*pi*(current_date%second - start_date%second &
                - starttime_wqs)/(stop_date%second + stoptime_wqs - start_date%second - starttime_wqs))) 
      end if
    case('INTERACT') ! moisture flux calculated interactively
        wqs(1:nlon,1:nlat)    = LE(1:nlon,1:nlat) / (rho_air * alv * g2kg)
    case('IMPORT') ! import moisture flux from measurements
      write(*,*) 'import moisture flux from measurements: not yet available'   
      stop
    case('SEA') ! heat fluxes over sea surface
      !wqs = wqs_sea
    case default 
      write(*,*) 'Flag for the function of wqs is invalid', f_wqs
      stop 
  end select

  wthetavs(1:nlon,1:nlat)    = wthetas(1:nlon, 1:nlat) + vtmpc1 * thetam(1:nlon,1:nlat,nlev) &
                                  *wqs(1:nlon,1:nlat) * g2kg  

  ! convert kinematic to dynamic heat fluxes
  SH                         = wthetas * (rho_air * cp_air)
  LE                         = wqs *  (rho_air * alv * g2kg)
  
  ! introduce advection of heat and moisture
  if ((current_date%second .le. start_date%second + starttime_adv) .or. &
          (current_date%second .ge. stop_date%second + stoptime_adv)) then
    advq(1:nlon,1:nlev)     = 0._dp
    advtheta(1:nlon,1:nlev) = 0._dp
  else 
    advq(1:nlon,1:nlev)     = advqmax
    advtheta(1:nlon,1:nlev) = advthetamax
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! calculate boundary layer dynamics
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! switch for applying different lapse rate
  if (lgamma .eqv. .true.) then
    do i = 1,nlon
     do j = 1,nlat
        if (hbl(i,j)>hcrit(i,j)) then
          gammatheta = gammath2
        endif
      enddo
    enddo
  endif  

  dthetav(1:nlon,1:nlat)     = dtheta(1:nlon,1:nlat) + vtmpc1*(qm(1:nlon,1:nlat,nlev)*dtheta(1:nlon,1:nlat) &
                                   + thetam(1:nlon,1:nlat,nlev)*dq(1:nlon,1:nlat) + dtheta*dq) * g2kg 

  we(1:nlon,1:nlat)          = beta * wthetavs(1:nlon,1:nlat)/dthetav(1:nlon,1:nlat)   

  ws(1:nlon,1:nlat)          = -omega * hbl(1:nlon,1:nlat)

  hbl_p1(1:nlon,1:nlat)      = hbl(1:nlon,1:nlat) + (we(1:nlon,1:nlat)+ws(1:nlon,1:nlat)) * delta_time

  dtheta_p1(1:nlon,1:nlat)   = dtheta(1:nlon,1:nlat) + ( gammatheta * we(1:nlon,1:nlat) &
                               - ((wthetas(1:nlon,1:nlat) + we(1:nlon,1:nlat)*dtheta(1:nlon,1:nlat)) &
                               /hbl(1:nlon,1:nlat) + advtheta(1:nlon,1:nlat)) ) * delta_time

  thetam_p1(1:nlon,1:nlat,nlev) = thetam(1:nlon,1:nlat,nlev) + ( (wthetas(1:nlon,1:nlat) &
                               + we(1:nlon,1:nlat)*dtheta(1:nlon,1:nlat))/hbl(1:nlon,1:nlat) & 
                               + advtheta(1:nlon,1:nlat) ) * delta_time 

  wqe(1:nlon,1:nlat)         = -we(1:nlon,1:nlat) * dq(1:nlon,1:nlat) 

  dq_p1(1:nlon,1:nlat)       = dq(1:nlon,1:nlat) + (gammaq * we(1:nlon,1:nlat) &
                               - (wqs(1:nlon,1:nlat) - wqe(1:nlon,1:nlat))/hbl(1:nlon,1:nlat) ) * delta_time

  qm_p1(1:nlon,1:nlat,nlev)  = qm(1:nlon,1:nlat,nlev) + ((wqs(1:nlon,1:nlat) &
                               - wqe(1:nlon,1:nlat))/hbl(1:nlon,1:nlat) + advq(1:nlon,1:nlat)) * delta_time
  
  uwe(1:nlon,1:nlat)         = -we(1:nlon,1:nlat) * du(1:nlon,1:nlat)

  vwe(1:nlon,1:nlat)         = -we(1:nlon,1:nlat) * dv(1:nlon,1:nlat)

  fc                         = TWOPI/(OneDay)*2._dp*sin(DTR*lat) 

  du_p1(1:nlon,1:nlat)       = du(1:nlon,1:nlat) + ((gammau*we(1:nlon,1:nlat)) - (uws(1:nlon,1:nlat) &
                               - uwe(1:nlon,1:nlat))/hbl(1:nlon,1:nlat) + fc*dv(1:nlon,1:nlat)) &
                               * delta_time

  dv_p1(1:nlon,1:nlat)       = dv(1:nlon,1:nlat) + ((gammav*we(1:nlon,1:nlat)) - (vws(1:nlon,1:nlat) &
                               - vwe(1:nlon,1:nlat))/hbl(1:nlon,1:nlat) - fc*du(1:nlon,1:nlat)) &
                               * delta_time
  
  um_p1(1:nlon,1:nlat,nlev)  = um(1:nlon,1:nlat,nlev) + (-fc*dv(1:nlon,1:nlat) &
                               + (uws(1:nlon,1:nlat)-uwe(1:nlon,1:nlat))/hbl(1:nlon,1:nlat)) &
                               * delta_time

  vm_p1(1:nlon,1:nlat,nlev)  = vm(1:nlon,1:nlat,nlev) + ( fc*du(1:nlon,1:nlat) &
                               + (vws(1:nlon,1:nlat)-vwe(1:nlon,1:nlat))/hbl(1:nlon,1:nlat)) &
                               * delta_time
  
  wthetave(1:nlon,1:nlat)    = -beta * wthetavs(1:nlon,1:nlat)
  betaq                      = -wqe/wqs


  ! temperature, moisture and pressure in the free troposphere (to be used in chemical reaction rates)
  ! theta_FT = theta_BL + dtheta: convert to T
  thetam_p1(1:nlon,1:nlat,nlev-1) = thetam_p1(1:nlon,1:nlat,nlev) + dtheta_p1(1:nlon,1:nlat) 
  ! q_FT = q_BL + dq         
  qm_p1(1:nlon,1:nlat,nlev-1)= qm_p1(1:nlon,1:nlat,nlev) + dq_p1(1:nlon,1:nlat)         

  um_p1(1:nlon,1:nlat,nlev-1)= um_p1(1:nlon,1:nlat,nlev) + du_p1(1:nlon,1:nlat)         
  vm_p1(1:nlon,1:nlat,nlev-1)= vm_p1(1:nlon,1:nlat,nlev) + dv_p1(1:nlon,1:nlat)         

  thetate(1:nlon,1:nlat,:)   = (thetam_p1(1:nlon,1:nlat,:)-thetam(1:nlon,1:nlat,:))/delta_time ! temp. tendency

  qte(1:nlon,1:nlat,:)       = (qm_p1(1:nlon,1:nlat,:)-qm(1:nlon,1:nlat,:))/delta_time ! sp. moisture tendency

  hbl                        = hbl_p1  
  dtheta                     = dtheta_p1
  thetam                     = thetam_p1
  dq                         = dq_p1
  qm                         = qm_p1
  du                         = du_p1
  um                         = um_p1
  dv                         = dv_p1
  vm                         = vm_p1

  ! saturation vapor pressure   
  esat(1:nlon,1:nlat,:)      = 0.611_dp/g2kg * exp(17.2694_dp * (thetam(1:nlon,1:nlat,:) - 273.16_dp) &
                                 / (thetam(1:nlon,1:nlat,:) - 35.86_dp))! Teten's formula
  ! vapor pressure at BL midpoint  
  e(1:nlon,1:nlat,nlev)      = qm(1:nlon,1:nlat,nlev)   * g2kg * (mpress_3d(1:nlon,1:nlat,nlev)) / (rd/rv) 
  !vapor pressure at BL-FT interface
  e(1:nlon,1:nlat,nlev-1)    = qm(1:nlon,1:nlat,nlev-1) * g2kg * (mpressi_3d(1:nlon,1:nlat,nlev)) / (rd/rv) 
  rh(1:nlon,1:nlat,nlev)     = 100.*(e(1:nlon,1:nlat,nlev)   / esat(1:nlon,1:nlat,nlev))
  rh(1:nlon,1:nlat,nlev-1)   = 100.*(e(1:nlon,1:nlat,nlev-1) / esat(1:nlon,1:nlat,nlev-1))
END SUBROUTINE mxl_algorithm

! -------------------------------------------------------------------------------------------------------

subroutine mxl_surfacelayer(pressi_3d, thetav, thetam, qm, um, vm, hsl, wthetas, wqs, ueff, &
                            rs, &
                            z0m, z0h, &
                            thetasurf, thetavsurf, qsurf,  uws, vws, T2m, q2m, u2m, v2m, &
                              u10m, v10m, rh2m, Rib, e2m, esat2m, ustar, Cm, Ch)

  USE messy_main_timer,         ONLY: lstart
  
  IMPLICIT NONE
  
  REAL(dp), DIMENSION(nlon,nlat,nlev+1), INTENT(in)   :: pressi_3d
  REAL(dp), DIMENSION(nlon,nlat,nlev),   INTENT(in)   :: thetav 
  REAL(dp), DIMENSION(nlon,nlat,nlev),   INTENT(in)   :: thetam
  REAL(dp), DIMENSION(nlon,nlat,nlev),   INTENT(in)   :: qm
  REAL(dp), DIMENSION(nlon,nlat,nlev),   INTENT(in)   :: um
  REAL(dp), DIMENSION(nlon,nlat,nlev),   INTENT(in)   :: vm    
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(in)   :: hsl     
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(in)   :: rs     
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(in)   :: z0h     
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(in)   :: z0m
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(in)   :: wthetas
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(in)   :: wqs

  REAL(dp), DIMENSION(nlon,nlat),        INTENT(out)  :: uws
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(out)  :: vws
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(out)  :: T2m     
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(out)  :: q2m    
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(out)  :: u2m    
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(out)  :: v2m    
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(out)  :: u10m
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(out)  :: v10m        
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(out)  :: rh2m    
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(out)  :: Rib     
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(out)  :: e2m         ! 2 meter vapor pressure (Pa)
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(out)  :: esat2m      ! 2 meter saturation vapor pressure (Pa)
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(out)  :: ustar       ! friction velocity (m s-1)
  REAL(dp), DIMENSION(nlon,nlat),        intent(inout):: Ch          ! drag coefficient for heat (-)
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(in)   :: ueff        ! effective wind velocity (m s-1)
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(inout):: thetasurf   ! surface potential temperature (K)
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(out)  :: thetavsurf ! surface virtual potential temp. (K) 

  REAL(dp), DIMENSION(nlon,nlat)        :: L           ! Obukhov length (m)
  REAL(dp), DIMENSION(nlon,nlat)        :: L0
  REAL(dp), DIMENSION(nlon,nlat)        :: Cq          ! drag coefficient for moisture (-)
  REAL(dp), DIMENSION(nlon,nlat)        :: Cm          ! drag coefficient for momentum (-)
  REAL(dp), DIMENSION(nlon,nlat)        :: esatsurf    ! surface saturation vapor pressure (Pa) 
  REAL(dp), DIMENSION(nlon,nlat)        :: qsatsurf    ! surface saturation specific humidity (g kg-1)
  REAL(dp), DIMENSION(nlon,nlat)        :: qsurf       ! surface specific humidity (g kg-1)

  REAL(dp)   :: Lbegin
  REAL(dp)   :: Lend
  REAL(dp)   :: fx
  REAL(dp)   :: fxdif
  INTEGER    :: iter, i, j

  thetasurf(1:nlon,1:nlat) = thetam(1:nlon,1:nlat,nlev) + wthetas(1:nlon,1:nlat) &
                               / (Ch(1:nlon,1:nlat) * ueff(1:nlon,1:nlat))

  esatsurf(1:nlon,1:nlat)  = 0.611_dp/g2kg * exp(17.2694_dp * (thetasurf(1:nlon,1:nlat) - 273.16_dp) &
                               / (thetasurf(1:nlon,1:nlat) - 35.86_dp))! Teten's formula
  qsatsurf(1:nlon,1:nlat)  = rd/rv * esatsurf(1:nlon,1:nlat) / pressi_3d(1:nlon,1:nlat,nlev+1)

  Cq(1:nlon,1:nlat)        = (1._dp + Ch(1:nlon,1:nlat) * ueff(1:nlon,1:nlat) * rs(1:nlon,1:nlat)) ** (-1._dp)

  qsurf(1:nlon,1:nlat)     = (1._dp - Cq(1:nlon,1:nlat)) * qm(1:nlon,1:nlat,nlev) + Cq(1:nlon,1:nlat) &
                               * qsatsurf(1:nlon,1:nlat) / g2kg 

  thetavsurf(1:nlon,1:nlat)= thetasurf(1:nlon,1:nlat) * (1._dp + vtmpc1 * qsurf(1:nlon,1:nlat) * g2kg)
  
  Rib(1:nlon,1:nlat)       = min(0.2_dp, g/thetav(1:nlon,1:nlat,nlev) * hsl(1:nlon,1:nlat) &
                  * (thetav(1:nlon,1:nlat,nlev) - thetavsurf(1:nlon,1:nlat)) / (ueff(1:nlon,1:nlat) ** 2._dp))
  L(1:nlon,1:nlat)         = sign(0.01_dp, Rib(1:nlon,1:nlat))
  L0(1:nlon,1:nlat)        = sign(0.1_dp, Rib(1:nlon,1:nlat))

  iter = 0
  do i = 1,nlon
    do j = 1,nlat
      do while(.true.) ! iteratively calculate Obukhov lenght
        iter      = iter + 1
        L0(i,j)   = L(i,j)
        fx        = Rib(i,j) - hsl(i,j) / L(i,j) * (log(hsl(i,j) / z0h(i,j))  - psih(hsl(i,j) / L(i,j)) &
                    + psih(z0h(i,j) / L(i,j))) / (log(hsl(i,j) / z0m(i,j)) - psim(hsl(i,j) / L(i,j)) &
                    + psim(z0m(i,j) / L(i,j))) ** 2._dp
        Lbegin    = L(i,j) - 0.001_dp*L(i,j)
        Lend      = L(i,j) + 0.001_dp*L(i,j)
        fxdif     = ( (- hsl(i,j) / Lbegin * (log(hsl(i,j) / z0h(i,j)) - psih(hsl(i,j) / Lbegin) &
                    + psih(z0h(i,j) / Lbegin)) / (log(hsl(i,j) / z0m(i,j)) - psim(hsl(i,j) / Lbegin) &
                    + psim(z0m(i,j) / Lbegin)) ** 2._dp) - (-hsl(i,j) / Lend * (log(hsl(i,j) / z0h(i,j)) &
                    - psih(hsl(i,j) / Lend) + psih(z0h(i,j) / Lend)) / (log(hsl(i,j) / z0m(i,j)) &
                    - psim(hsl(i,j) / Lend) + psim(z0m(i,j) / Lend)) ** 2._dp) ) / (Lbegin - Lend)
        L(i,j)    = L(i,j) - fx / fxdif
        L(i,j)    = sign(min(abs(L(i,j)),1e6_dp),L(i,j)) ! capping L
      
        if(abs((L(i,j) - L0(i,j))/L(i,j)) < 1e-4_dp) exit  
        if(abs((L(i,j) - L0(i,j))) < 1e-3_dp) exit        
      enddo
    enddo
  enddo

  Cm(1:nlon,1:nlat)     = c_vKar ** 2._dp / (log(hsl(1:nlon,1:nlat) / z0m(1:nlon,1:nlat)) &
                          - psim(hsl(1:nlon,1:nlat) / L(1:nlon,1:nlat)) + psim(z0m(1:nlon,1:nlat) &
                          / L(1:nlon,1:nlat))) ** 2._dp
  Ch(1:nlon,1:nlat)     = c_vKar ** 2._dp / (log(hsl(1:nlon,1:nlat) / z0m(1:nlon,1:nlat)) &
                          - psim(hsl(1:nlon,1:nlat) / L(1:nlon,1:nlat)) + psim(z0m(1:nlon,1:nlat) &
                          / L(1:nlon,1:nlat))) / (log(hsl(1:nlon,1:nlat) / z0h(1:nlon,1:nlat)) &
                          - psih(hsl(1:nlon,1:nlat) / L(1:nlon,1:nlat)) &
                          + psih(z0h(1:nlon,1:nlat) / L(1:nlon,1:nlat)))

  ustar(1:nlon,1:nlat)  = sqrt(Cm(1:nlon,1:nlat)) * ueff(1:nlon,1:nlat)

  uws(1:nlon,1:nlat)    = - Cm(1:nlon,1:nlat) * ueff(1:nlon,1:nlat) * um(1:nlon,1:nlat,nlev)
  vws(1:nlon,1:nlat)    = - Cm(1:nlon,1:nlat) * ueff(1:nlon,1:nlat) * vm(1:nlon,1:nlat,nlev)

  T2m(1:nlon,1:nlat)    = thetasurf(1:nlon,1:nlat) - wthetas(1:nlon,1:nlat) / ustar(1:nlon,1:nlat) & 
                          / c_vKar * (log(2._dp / z0h(1:nlon,1:nlat)) - psih(2._dp &
                         / L(1:nlon,1:nlat)) + psih(z0h(1:nlon,1:nlat) / L(1:nlon,1:nlat)))
  q2m(1:nlon,1:nlat)    = qsurf(1:nlon,1:nlat)     - wqs(1:nlon,1:nlat)     / ustar(1:nlon,1:nlat) &
                          / c_vKar * (log(2._dp / z0h(1:nlon,1:nlat)) - psih(2._dp &
                         / L(1:nlon,1:nlat)) + psih(z0h(1:nlon,1:nlat) / L(1:nlon,1:nlat)))

  u2m(1:nlon,1:nlat)    =                          - uws(1:nlon,1:nlat)     / ustar(1:nlon,1:nlat) &
                          / c_vKar * (log(2._dp / z0m(1:nlon,1:nlat)) - psim(2._dp &
                         / L(1:nlon,1:nlat)) + psim(z0m(1:nlon,1:nlat) / L(1:nlon,1:nlat)))
  v2m(1:nlon,1:nlat)    =                          - vws(1:nlon,1:nlat)     / ustar(1:nlon,1:nlat) &
                          / c_vKar * (log(2._dp / z0m(1:nlon,1:nlat)) - psim(2._dp &
                         / L(1:nlon,1:nlat)) + psim(z0m(1:nlon,1:nlat) / L(1:nlon,1:nlat)))
  u10m(1:nlon,1:nlat)    =                          - uws(1:nlon,1:nlat)     / ustar(1:nlon,1:nlat) &
                          / c_vKar * (log(10._dp / z0m(1:nlon,1:nlat)) - psim(10._dp &
                         / L(1:nlon,1:nlat)) + psim(z0m(1:nlon,1:nlat) / L(1:nlon,1:nlat)))

  v10m(1:nlon,1:nlat)    =                          - vws(1:nlon,1:nlat)     / ustar(1:nlon,1:nlat) &
                          / c_vKar * (log(10._dp / z0m(1:nlon,1:nlat)) - psim(10._dp &
                         / L(1:nlon,1:nlat)) + psim(z0m(1:nlon,1:nlat) / L(1:nlon,1:nlat)))
  esat2m(1:nlon,1:nlat) = 0.611_dp/g2kg * exp(17.2694_dp * (T2m(1:nlon,1:nlat) - 273.16_dp) &
                          / (T2m(1:nlon,1:nlat) - 35.86_dp)) ! Teten's formula
  e2m(1:nlon,1:nlat)    = q2m(1:nlon,1:nlat) * g2kg * (pressi_3d(1:nlon,1:nlat,nlev+1)) / (rd/rv)  
  rh2m(1:nlon,1:nlat)   = e2m(1:nlon,1:nlat) / esat2m(1:nlon,1:nlat)
end subroutine mxl_surfacelayer  
!------------------------------------------------------------------------------------------------------  
subroutine mxl_momentumflux(l_ustconst, um, vm, z0m, uws, vws, ustar)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! calculate surface momentum fluxes
  !!! 2 options: constant momentum fluxes, 
  !!!            or time-dependent uws, vws and ustar (z0 constant)  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

  IMPLICIT NONE

  INTEGER                                              :: i,j
  LOGICAL, intent(in)                                  :: l_ustconst
  REAL(dp), DIMENSION(nlon,nlat,nlev),   INTENT(in)    :: um
  REAL(dp), DIMENSION(nlon,nlat,nlev),   INTENT(in)    :: vm
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(in)    :: z0m
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(inout) :: uws
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(inout) :: vws
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(out)   :: ustar
  REAL(dp), DIMENSION(nlon,nlat)                       :: zp
  REAL(dp), DIMENSION(nlon,nlat)                       :: alpha

  zp  = 10. ! Height for the calculation of the logarithmic equation for u* (m)

  if (l_ustconst) then ! constant
        ustar = ((uws**2._dp) + (vws**2._dp))**0.25_dp
  else ! time-dependent
        ustar(1:nlon,1:nlat) = c_vKar*sqrt(um(1:nlon,1:nlat,nlev)**2._dp + vm(1:nlon,1:nlat,nlev)**2._dp) &
                               /log(zp(1:nlon,1:nlat)/z0m(1:nlon,1:nlat))
  
        do i = 1,nlon
          do j = 1,nlat
            if (um(i,j,nlev) .ne. 0.) then
              alpha(i,j) = atan(vm(i,j,nlev)/um(i,j,nlev))
            else
              if (vm(i,j,nlev) .ne. 0.) then
                alpha(i,j) = PI_2
              endif
            endif
          enddo
        enddo
        uws = -ustar**2._dp * cos(alpha)
        vws = -ustar**2._dp * sin(alpha)
  endif

end subroutine mxl_momentumflux
!------------------------------------------------------------------------------------------------------
subroutine mxl_landsurface(l_surfacelayer, l_radiation, &
                           rsveg,ra,rssoil,SH,LE,GR,Tsurf, &
                           pressi_3d,thetam,qm,ueff, &
                           Swin,Rn, &
                           e2m,esat2m,T2m,ustar, &
                           LAI,wwilt,w2,w1,wfc,wsat,CLa,CLb,CLc,C1sat,C2ref,gD,rsmin, &
                             rssoilmin,Lambda,cveg,Wl,CGsat,Tsoil1,Tsoil2)

  USE messy_main_timer,      ONLY: delta_time
  
  IMPLICIT NONE
  
  LOGICAL,                               INTENT(in)   :: l_surfacelayer
  LOGICAL,                               INTENT(in)   :: l_radiation
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(out)  :: SH
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(out)  :: LE
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(out)  :: GR
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(out)  :: ra
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(out)  :: rsveg
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(out)  :: rssoil
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(out)  :: Tsurf

  REAL(dp), DIMENSION(nlon,nlat,nlev),   INTENT(in)   :: thetam
  REAL(dp), DIMENSION(nlon,nlat,nlev),   INTENT(in)   :: qm
  REAL(dp), DIMENSION(nlon,nlat,nlev+1), INTENT(in)   :: pressi_3d
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(in)   :: Swin
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(in)   :: Rn
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(in)   :: esat2m
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(in)   :: e2m
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(in)   :: T2m
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(in)   :: LAI
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(in)   :: ustar
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(in)   :: ueff
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(in)   :: wwilt
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(in)   :: w2
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(in)   :: wfc
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(in)   :: wsat
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(in)   :: CLa
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(in)   :: CLb
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(in)   :: CLc
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(in)   :: C1sat
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(in)   :: C2ref
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(in)   :: gD
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(in)   :: rsmin
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(in)   :: rssoilmin
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(in)   :: Lambda
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(in)   :: cveg
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(in)   :: CGsat

  REAL(dp), DIMENSION(nlon,nlat),        INTENT(inout)   :: Tsoil1
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(inout)   :: Tsoil2
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(inout)   :: Wl
  REAL(dp), DIMENSION(nlon,nlat),        INTENT(inout)   :: w1

  REAL(dp), DIMENSION(nlon,nlat) :: desatdT
  REAL(dp), DIMENSION(nlon,nlat) :: dqsatdT
  REAL(dp), DIMENSION(nlon,nlat) :: efinal
  REAL(dp), DIMENSION(nlon,nlat) :: Wlmx
  REAL(dp), DIMENSION(nlon,nlat) :: Wmax 
  REAL(dp), DIMENSION(nlon,nlat) :: cliq
  REAL(dp), DIMENSION(nlon,nlat) :: LEveg
  REAL(dp), DIMENSION(nlon,nlat) :: LEliq
  REAL(dp), DIMENSION(nlon,nlat) :: LEsoil
  REAL(dp), DIMENSION(nlon,nlat) :: CG
  REAL(dp), DIMENSION(nlon,nlat) :: C1
  REAL(dp), DIMENSION(nlon,nlat) :: C2
  REAL(dp), DIMENSION(nlon,nlat) :: dTsoil1_dt    
  REAL(dp), DIMENSION(nlon,nlat) :: w1eq
  REAL(dp), DIMENSION(nlon,nlat) :: dw1_dt    
  REAL(dp), DIMENSION(nlon,nlat) :: esat    
  REAL(dp), DIMENSION(nlon,nlat) :: qsat 
  REAL(dp), DIMENSION(nlon,nlat) :: f1    
  REAL(dp), DIMENSION(nlon,nlat) :: f2    
  REAL(dp), DIMENSION(nlon,nlat) :: f3    
  REAL(dp), DIMENSION(nlon,nlat) :: f4    
  REAL(dp), DIMENSION(nlon,nlat) :: esatsurf
  REAL(dp), DIMENSION(nlon,nlat) :: qsatsurf
  REAL(dp), DIMENSION(nlon,nlat) :: dWl_dt
  integer :: i,j

  Wmax =  2.0e-4_dp

  ra(1:nlon,1:nlat)       = ueff(1:nlon,1:nlat) / (ustar(1:nlon,1:nlat) ** 2._dp)
  esat(1:nlon,1:nlat)     = 0.611e3_dp * exp(17.2694_dp * (thetam(1:nlon,1:nlat,nlev) - 273.16_dp) / &
                                                          (thetam(1:nlon,1:nlat,nlev) - 35.86_dp))
  qsat(1:nlon,1:nlat)     = (rd/rv) * esat(1:nlon,1:nlat) / (pressi_3d(1:nlon,1:nlat,nlev+1))
  desatdT(1:nlon,1:nlat)  = esat(1:nlon,1:nlat) * (17.2694_dp / (thetam(1:nlon,1:nlat,nlev) - 35.86_dp) & 
                            - 17.2694_dp * (thetam(1:nlon,1:nlat,nlev) - 273.16_dp) / &
                            (thetam(1:nlon,1:nlat,nlev) - 35.86_dp)**2._dp)
  dqsatdT(1:nlon,1:nlat)  = (rd/rv) * desatdT(1:nlon,1:nlat) / pressi_3d(1:nlon,1:nlat,nlev+1)
  efinal(1:nlon,1:nlat)   = qm(1:nlon,1:nlat,nlev) * g2kg * pressi_3d(1:nlon,1:nlat,nlev+1) / (rd/rv) 

  if (l_radiation) then
    f1(1:nlon,1:nlat) = 1.0_dp / ((0.004_dp * Swin(1:nlon,1:nlat) + 0.05_dp) / &
                       (0.81_dp * (0.004_dp * Swin(1:nlon,1:nlat) + 1._dp)))
  else
    f1(1:nlon,1:nlat) = 1.0_dp
  endif

  do i = 1,nlon
    do j = 1,nlat
      if (w2(i,j) .gt. wwilt(i,j)) then
        f2(i,j)     = (wfc(i,j) - wwilt(i,j)) / (w2(i,j) - wwilt(i,j))
      else 
        f2(i,j)     = 1.e8_dp
      endif
    enddo
  enddo
  f2     = min(1.e8_dp, f2)

  if (l_surfacelayer) then
    f3(1:nlon,1:nlat) = 1.0_dp / exp( - gD * (esat2m(1:nlon,1:nlat) - e2m(1:nlon,1:nlat)) )
    f4(1:nlon,1:nlat) = 1.0_dp / (1.0_dp - 0.0016_dp * (298.0_dp - T2m(1:nlon,1:nlat)) ** 2._dp)
  else ! use BL-averaged values instead of 2 m values   
    f3(1:nlon,1:nlat) = 1.0_dp / exp( - gD * (esat(1:nlon,1:nlat) - efinal(1:nlon,1:nlat)) )
    f4(1:nlon,1:nlat) = 1.0_dp / (1.0_dp - 0.0016_dp * (298.0_dp - thetam(1:nlon,1:nlat,nlev)) ** 2._dp)
  endif

  rsveg(1:nlon,1:nlat)   = (rsmin(1:nlon,1:nlat) / LAI(1:nlon,1:nlat)) * f1(1:nlon,1:nlat) * & 
                                  f2(1:nlon,1:nlat) * f3(1:nlon,1:nlat) * f4(1:nlon,1:nlat)   

  ! recompute f2 using w1 instead of w2
  f2     = 1.e8_dp
  do i = 1,nlon
    do j = 1,nlat
      if (w1(i,j) .gt. wwilt(i,j)) then
        f2(i,j)   = (wfc(i,j) - wwilt(i,j)) / (w1(i,j) - wwilt(i,j))
      endif       
    enddo
  enddo
  f2     = min(1.e8_dp, f2)

  rssoil(1:nlon,1:nlat) = rssoilmin(1:nlon,1:nlat) * f2(1:nlon,1:nlat)

  Wlmx(1:nlon,1:nlat)   = LAI(1:nlon,1:nlat) * Wmax(1:nlon,1:nlat)
  cliq(1:nlon,1:nlat)   = min(1.0_dp, Wl(1:nlon,1:nlat) / Wlmx(1:nlon,1:nlat))

  Tsurf(1:nlon,1:nlat)  = (Rn + rho_air * cp_air / ra(1:nlon,1:nlat) * thetam(1:nlon,1:nlat,nlev) &
      + cveg(1:nlon,1:nlat) * (1.0_dp-cliq(1:nlon,1:nlat)) * rho_air * alv / (ra(1:nlon,1:nlat) &
      + rsveg(1:nlon,1:nlat)) * (dqsatdT(1:nlon,1:nlat) * thetam(1:nlon,1:nlat,nlev) - qsat(1:nlon,1:nlat) &
      + qm(1:nlon,1:nlat,nlev) * g2kg) + (1.0_dp - cveg(1:nlon,1:nlat)) * rho_air * alv / (ra(1:nlon,1:nlat) &
      + rssoil(1:nlon,1:nlat)) * (dqsatdT(1:nlon,1:nlat) * thetam(1:nlon,1:nlat,nlev) - qsat(1:nlon,1:nlat) &
      + qm(1:nlon,1:nlat,nlev) * g2kg) &
      + cveg(1:nlon,1:nlat) * cliq(1:nlon,1:nlat) * rho_air * alv / ra(1:nlon,1:nlat) &
       * (dqsatdT(1:nlon,1:nlat) * thetam(1:nlon,1:nlat,nlev) &
      - qsat(1:nlon,1:nlat) + qm(1:nlon,1:nlat,nlev) * g2kg) &
      + Lambda * Tsoil1(1:nlon,1:nlat)) * (rho_air * cp_air / ra(1:nlon,1:nlat) &
      + cveg(1:nlon,1:nlat) * (1._dp - cliq(1:nlon,1:nlat)) * rho_air * alv / (ra(1:nlon,1:nlat) &
      + rsveg(1:nlon,1:nlat)) * dqsatdT(1:nlon,1:nlat) + (1._dp - cveg(1:nlon,1:nlat)) * rho_air * &
         alv / (ra(1:nlon,1:nlat) &
      + rssoil(1:nlon,1:nlat)) * dqsatdT(1:nlon,1:nlat) &
      + cveg(1:nlon,1:nlat) * cliq(1:nlon,1:nlat) * rho_air * alv / ra(1:nlon,1:nlat) * &
        dqsatdT(1:nlon,1:nlat) + Lambda) ** (-1._dp)

  esatsurf(1:nlon,1:nlat) = 0.611e3_dp * exp(17.2694_dp * (Tsurf(1:nlon,1:nlat) - 273.16_dp) &
                            / (Tsurf(1:nlon,1:nlat) - 35.86_dp))
  qsatsurf(1:nlon,1:nlat) = (rd/rv) * esatsurf(1:nlon,1:nlat) / pressi_3d(1:nlon,1:nlat,nlev+1)

  LEveg(1:nlon,1:nlat)  = (1.0_dp - cliq(1:nlon,1:nlat)) * cveg(1:nlon,1:nlat)  * rho_air * alv &
                          / (ra(1:nlon,1:nlat) + rsveg(1:nlon,1:nlat)) * (dqsatdT(1:nlon,1:nlat) * &
                          (Tsurf(1:nlon,1:nlat) - thetam(1:nlon,1:nlat,nlev)) + qsat(1:nlon,1:nlat) &
                          - qm(1:nlon,1:nlat,nlev) * g2kg)
  LEliq(1:nlon,1:nlat)  = cliq(1:nlon,1:nlat)  * cveg(1:nlon,1:nlat)  * rho_air * alv /  ra(1:nlon,1:nlat) &
                          * (dqsatdT(1:nlon,1:nlat) * (Tsurf(1:nlon,1:nlat) - thetam(1:nlon,1:nlat,nlev)) &
                          + qsat(1:nlon,1:nlat) - qm(1:nlon,1:nlat,nlev) * g2kg)
  LEsoil(1:nlon,1:nlat) = (1.0_dp - cveg(1:nlon,1:nlat)) * rho_air * alv / (ra(1:nlon,1:nlat) &
                          + rssoil(1:nlon,1:nlat)) * (dqsatdT(1:nlon,1:nlat) * (Tsurf(1:nlon,1:nlat) &
                        - thetam(1:nlon,1:nlat,nlev)) + qsat(1:nlon,1:nlat) - qm(1:nlon,1:nlat,nlev) * g2kg)

  dWl_dt(1:nlon,1:nlat) = - LEliq(1:nlon,1:nlat) / (rho_H2O * alv)
  Wl(1:nlon,1:nlat)     = Wl(1:nlon,1:nlat) + dWl_dt(1:nlon,1:nlat) * delta_time

  LE(1:nlon,1:nlat)     = LEsoil(1:nlon,1:nlat) + LEveg(1:nlon,1:nlat) + LEliq(1:nlon,1:nlat)
  SH(1:nlon,1:nlat)     = rho_air * cp_air / ra(1:nlon,1:nlat) * &
                              (Tsurf(1:nlon,1:nlat) - thetam(1:nlon,1:nlat,nlev))
  GR(1:nlon,1:nlat)     = Lambda * (Tsurf(1:nlon,1:nlat) - Tsoil1(1:nlon,1:nlat))


  CG(1:nlon,1:nlat)     = CGsat * (wsat / w2(1:nlon,1:nlat)) ** (CLb / (2.0_dp * log(10.0_dp)))

  dTsoil1_dt(1:nlon,1:nlat) = CG(1:nlon,1:nlat) * GR(1:nlon,1:nlat) - TWOPI / OneDay &
                             * (Tsoil1(1:nlon,1:nlat) - Tsoil2(1:nlon,1:nlat))
  Tsoil1(1:nlon,1:nlat) = Tsoil1(1:nlon,1:nlat) + dTsoil1_dt(1:nlon,1:nlat) * delta_time

  C1(1:nlon,1:nlat)     = C1sat * (wsat / w1(1:nlon,1:nlat)) ** (CLb / 2.0_dp + 1.0_dp)
  C2(1:nlon,1:nlat)     = C2ref * (w2(1:nlon,1:nlat) / (wsat - w2(1:nlon,1:nlat)) )
  w1eq(1:nlon,1:nlat)   = w2(1:nlon,1:nlat) - wsat * CLa * ( (w2(1:nlon,1:nlat) / wsat)**CLc &
                          * (1.0_dp - (w2(1:nlon,1:nlat) / wsat) ** (8.0_dp * CLc)) )

  dw1_dt(1:nlon,1:nlat) = - C1 / (rho_H2O * 0.1_dp) * LEsoil(1:nlon,1:nlat) / alv &
                          - C2 / OneDay * (w1(1:nlon,1:nlat) - w1eq(1:nlon,1:nlat))

  w1(1:nlon,1:nlat)     = w1(1:nlon,1:nlat) + dw1_dt(1:nlon,1:nlat) * delta_time

end subroutine mxl_landsurface
!------------------------------------------------------------------------------------------------------
subroutine mxl_radiation(shsl,ssurfpress,sthetam, &
                         sCc,salbedo,sTsurf, &
                         Swin,Swout,Lwin,Lwout,Rn,costh)

  use messy_main_timer,   only : DAYOFYEAR
  
  IMPLICIT NONE
  
  REAL(dp), DIMENSION(nlon,nlat),        intent(in) :: shsl
  REAL(dp), DIMENSION(nlon,nlat,nlev+1), intent(in) :: ssurfpress
  REAL(dp), DIMENSION(nlon,nlat,nlev),   intent(in) :: sthetam
  REAL(dp), DIMENSION(nlon,nlat),        intent(in) :: sCc
  REAL(dp), DIMENSION(nlon,nlat),        intent(in) :: salbedo
  REAL(dp), DIMENSION(nlon,nlat),        intent(in) :: sTsurf
  REAL(dp) :: DOYreal
  REAL(dp) :: timeofday
  REAL(dp), dimension(nlon,nlat) :: Ta
  REAL(dp), dimension(nlon,nlat) :: Tr
  REAL(dp), dimension(nlon,nlat), intent(out) :: Swin
  REAL(dp), dimension(nlon,nlat), intent(out) :: Swout
  REAL(dp), dimension(nlon,nlat), intent(out) :: Lwin
  REAL(dp), dimension(nlon,nlat), intent(out) :: Lwout
  REAL(dp), dimension(nlon,nlat), intent(out) :: Rn
  REAL(dp),                       intent(out) :: costh  ! cosine of solar zenith angel
  REAL(dp), parameter :: etaa = 0.8      ! emissivity of the atmosphere
  REAL(dp), parameter :: etas = 1.0      ! emissivity of the surface           

  timeofday = current_date%second/3600._dp
  DOYreal   = real(DAYOFYEAR)
  costh     = cos(getsoz(lat, lon, DOYreal, timeofday)) 

  Ta(1:nlon,1:nlat) = sthetam(1:nlon,1:nlat,nlev)*((((ssurfpress(1:nlon,1:nlat,nlev+1)) &
                  - shsl(1:nlon,1:nlat)*rho_air*g) / (ssurfpress(1:nlon,1:nlat,nlev+1)))**(rd/cpv))
  Tr    = (0.6_dp + 0.2_dp * max(0.0_dp,costh)) * (1._dp - 0.4_dp * sCc)
 
  Swin(1:nlon,1:nlat)  = solc * Tr(1:nlon,1:nlat) * max(0.0_dp,costh)
  Swout(1:nlon,1:nlat) = -salbedo(1:nlon,1:nlat) * Swin
  Lwin(1:nlon,1:nlat)  =  etaa * stbo * (Ta(1:nlon,1:nlat) ** 4_dp)
  Lwout(1:nlon,1:nlat) = -etas * stbo * (sTsurf(1:nlon,1:nlat) ** 4_dp)

  Rn  = Swin + Swout + Lwin + Lwout

end subroutine mxl_radiation
!-----------------------------------------------------------------------------
SUBROUTINE emis_simple(f_emisfunc,maxflux, starttime_emis, stoptime_emis, actflux)

  IMPLICIT NONE
  
  CHARACTER(LEN=STRLEN_MEDIUM),   INTENT(in)  :: f_emisfunc
!  CHARACTER(LEN=STRLEN_MEDIUM),   INTENT(in)  :: tracer
  REAL(dp),                       INTENT(in)  :: maxflux
  REAL(dp),                       INTENT(in)  :: starttime_emis
  REAL(dp),                       INTENT(in)  :: stoptime_emis
  REAL(dp),                       INTENT(out) :: actflux
  
  SELECT CASE(f_emisfunc)
    CASE('NOEMIS') ! no emission flux
      actflux = 0.
    CASE('CONST') ! constant emission flux with start & end time
      IF ((current_date%second .le. start_date%second+starttime_emis) .or. &
          (current_date%second .ge. stop_date%second+stoptime_emis)) then
        actflux = 0. 
      ELSE
        actflux = maxflux 
      END IF
    CASE('SINE') ! sinusoidal emission flux
      IF ((current_date%second .le. start_date%second + starttime_emis) .or. &
          (current_date%second .ge. stop_date%second + stoptime_emis)) then
        actflux = 0.
      ELSE
        actflux = maxflux * sin(pi * (current_date%second - start_date%second - starttime_emis)/ &
                           (stop_date%second + stoptime_emis - start_date%second - starttime_emis)) 
      END IF
    CASE('COSINE') ! cosine emission flux
      IF ((current_date%second .le. start_date%second + starttime_emis) .or. &
          (current_date%second .ge. stop_date%second + stoptime_emis)) then
        actflux = 0.
      ELSE
        actflux = (maxflux/2._dp)*(1._dp - cos(2._dp*pi*(current_date%second - start_date%second &
                - starttime_emis)/(stop_date%second + stoptime_emis - start_date%second - starttime_emis))) 
      END IF
    CASE default 
      write(*,*) 'Flag for the emission function is invalid ', f_emisfunc
      stop 
  END SELECT
  
END SUBROUTINE emis_simple
!-------------------------------------------------------------------------------

  REAL FUNCTION getsoz(latd, lond, DoY, hourofday)
    ! Function to calculate solar zenith angle for -65< lat <65 degrees north
    ! and lon between -60 and + 60

    IMPLICIT NONE

    REAL(dp), intent(in) :: DoY
    REAL(dp), intent(in) :: latd,lond
    REAL(dp), intent(in) :: hourofday 
    REAL(dp) :: houra
    REAL(dp) :: phir, deday, delta, lonr, latr
    REAL(dp) :: hutc ! time of day in hours UTC
    INTEGER  :: status
    REAL(dp) :: dateback

    status = 0

    IF ( (lond < 0.0_dp) .OR. (lond >= 360.0_dp) ) THEN
       WRITE(*,*) 'ERROR in lt2utc: Longitude out of bounds (lon = ', lond, ')!'
       status = 1
       RETURN
    ELSE IF (lond <= 180.0_dp) THEN
       dateback = 0._dp
    ELSE IF (lond > 180.0_dp) THEN
       dateback = 1._dp
    ENDIF

    ! [h]       [h]           [deg]    [°/day]          [hour/day]
    hutc  = hourofday - nint((lond/360.0_dp - dateback) * (OneDay/3600._dp)) ! convert local time to UTC
    lonr  = lond*DTR
    latr  = latd*DTR
    phir  = 23.45_dp * DTR
    deday = 4.88_dp + TWOPI/365_dp * DoY
    delta = asin(sin(phir)*sin(deday))
    houra = lonr - pi + hutc * (TWOPI/24._dp)
    getsoz = acos(sin(delta)*sin(latr) + cos(delta)*cos(latr)*cos(houra))
  END FUNCTION getsoz

!----------------------------------------------------------------------------
  elemental real(dp) function psim(zeta)
    
    implicit none
    
    REAL(dp), intent(in) :: zeta
    REAL(dp) :: x

    if(zeta <= 0_dp) then
      x    = (1._dp - 16._dp * zeta) ** (0.25_dp)
      psim = pi / 2._dp - 2._dp * atan(x) + log( (1._dp+x) ** 2._dp * (1._dp + x ** 2._dp) / 8._dp)
    else
      psim = -2._dp/3._dp * (zeta - 5._dp/0.35_dp) * exp(-0.35_dp * zeta) &
              - zeta - (10._dp/3._dp) / 0.35_dp
    endif
    return
  end function
!----------------------------------------------------------------------------
  elemental real(dp) function psih(zeta)
    
    implicit none
    
    REAL(dp), intent(in) :: zeta
    REAL(dp) :: x

    if(zeta <= 0_dp) then
      x    = (1._dp - 16._dp * zeta) ** (0.25_dp)
      psih = 2._dp * log( (1._dp + x ** 2._dp) / 2._dp )
    else
      psih = -2._dp/3._dp * (zeta - 5._dp/0.35_dp) * exp(-0.35_dp * zeta) &
             - (1._dp + (2._dp/3._dp) * zeta) ** (1.5_dp) - (10._dp/3._dp) / 0.35_dp + 1._dp
    endif
    return
  end function
!----------------------------------------------------------------------------

  SUBROUTINE mxl_read_nml_ctrl(status, iou)
 
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /CTRL/ lat, lon, l_verbose, l_chem_ft

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='mxl_read_nml_ctrl'
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    ! INITIALIZE NAMELIST VARIABLES
    l_verbose = .false.

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist
    
    CALL read_nml_close(substr, iou, modstr)

    status = 0 ! NO ERROR
    
  END SUBROUTINE mxl_read_nml_ctrl


END MODULE
