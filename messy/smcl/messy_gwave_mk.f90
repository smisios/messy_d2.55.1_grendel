MODULE messy_gwave_mk
!r                        
!r This module calculates gravity wave deposition and eddy diffusion 
!r coefficient, based on the Parameterisation of Medvedev & Klaassen 95.
!r 
!r The code was kindly sent to me by Gary Klassen, I converted to
!r Fortran 90 with some refactoring. I also added functionality in
!r this module to set atmospheric parameters below the lower boundary and
!r interpolation wrappers. I also added latitudinal variation of spectral
!r density and eddy diffusion calculation.
!r
!r The gravity wave calculation is made on a finer grid than CMAT2. This
!r is required by the parameterization. Also the gwave grid may extend below
!r the lower boundary of CMAT2. The specification of this grid is set
!r by variables plevs, lev_step, and p0_gw below. Note, the low boundary
!r of the gwave grid is where the gwave spectrum is defined, and is 
!r typically set at 165mb.
!r
!r See spectrum definition below after comment 'GRAVITY WAVE SPECTRUM'.
!r 


  USE messy_main_constants_mem, ONLY: g & 
       , prcn=>dp, pi &
       , rd  ! Gas constant

  IMPLICIT NONE

  PRIVATE

  REAL(prcn), PUBLIC :: gwave_int_fac = 1.

  !v
  !v Note, the parameterization is sensative to winds on the lower
  !v boundary. If you set anisotropic=0 winds at u(0) are assumed
  !v 0. Anisotropic also activates the projection of winds onto local wind
  !v wind vector, basically setting NS drag to zero, only including
  !v waves that travel EW.
  !v
  !v anistropic = 1/0
  !v

  INTEGER               :: anisotropic = 0

  !v
  !v These variables define the grid used for gwave calculation. Note
  !v it has to have at least 0.25 scale height steps for the parameterization
  !v to work.

  !v Lower boundary pressure of Gwave Grid . Note, this
  !v should correspond to the spectrum definition height
  REAL(prcn), PUBLIC            :: p0_gw = 165.*1.e2

  !v
  !v GRAVITY WAVE SPECTRUM at p0_gw
  !v
  INTEGER, PARAMETER, PUBLIC :: nhar   = 15   ! Number of harmonics in spectrum
  INTEGER, PUBLIC            :: levgws = 1    ! Level where spectrum defined gwave grid
  REAL(prcn) :: xkzl  = 2.*pi/19.e+3 ! lower kz corresponds to lambda_z=19 km
  REAL(prcn) :: xkzu  = 2.*pi/0.9e+3 ! upper kz corresponds to lambda_z=900 m
  REAL(prcn) :: mstar = 0.006        ! Modified DeSaubies spectrum 0.006
  REAL(prcn) :: amp   = 50.          ! Av. ampl of Modified DeSaubies spectrum
  !REAL(prcn) :: amp_ratio = 0.2      ! Ratio of meridional to zonal amp
  REAL(prcn) :: gwkh  = 2.*pi/3.0e+5 ! 300 km wave length hoz wavenumber

  !v Linear damping term factor. This is used for calculating beta dissipation
  !v passed the turbopause. mjh
  REAL(prcn) :: visc_fac = 1.e-2

  !+
  ! common variables for MK95 gravity wave parameterization
  !-
  REAL(prcn)   ::             &
         rg,                  & ! gas constant/acceleration of gravity
         sqrt2,               & ! sqrt(2)
         sqrtp,               & ! sqrt(pi)
         sqrt2p,              & ! sqrt(2pi)
         sqrtkh,              & ! sqrt(2)*horizontal wavenumber (gwkh)
         dxkz,                & ! delta kz of the spectrum at the lower bndry
         zkz(nhar),           & ! kz's of the spectrum at the lower bndry
         tndmax,              & ! max value of gw tendency
         dtempmax,            & ! max value of gw heat flux divergence
         eddymax                ! max value of eddy diffusion. mjh

  ! squared amplitudes of the l.b. spectrum
  REAL(prcn), ALLOCATABLE, DIMENSION(:,:,:), PUBLIC  :: usql
  ! Latitudinal strength variation factor
  REAL(prcn), ALLOCATABLE, DIMENSION(:), PUBLIC      :: fac


  PUBLIC :: mkgwintr
  PUBLIC :: mkgwinti
  PUBLIC :: splineon
  PUBLIC :: mk_read_nml_ctrl

  
CONTAINS

  !-----------------------------------------------------------------------
  ! Interface for MK95 gravity wave drag parameterization.
  ! Isotropic GW source spectrum.        A. Medvedev   Oct 12/1999
  !
  ! Rewritten a little to use 1d arrays instead of 2d CCM arrays. 
  ! mjh Jun 2004.
  !-----------------------------------------------------------------------

   SUBROUTINE mkgwintr (uu, vv, tf, th, pf, ph, dpf, bvfr, drag_u,           &
        drag_v, temp_gw, lat, fac, scp, mk_eddy, svisc, plev)


     IMPLICIT NONE

     INTEGER, INTENT(IN) :: plev

     real(prcn), INTENT(INOUT) ::     &
      bvfr(plev),                  & ! "full-level" Brunt-Vaisalla frequency
      dpf(plev),                   & ! "full-level" delta p (=-dpm)
      pf(plev),                    & ! "full-level" pressure
      ph(plev),                    & ! "half-level" pressure
      tf(plev),                    & ! "full-level" temperature
      th(plev),                    & ! "half-level" temperature
      uu(plev),                    & ! "full-level" background zonal wind
      vv(plev),                    & ! "full-level" background merid wind
      fac(:),                      & ! Latitudinal variation in strength. mjh
      scp(plev),                   & ! Specific heat capacity
      svisc(plev)                    ! Coefficient of viscosity

     ! Local
     real(prcn) ::     &
      ax(plev),                    & ! gwd drag in the vertical column
      omeg(nhar,plev),             & ! intrinsic frequency
      ttpgw(plev),                 & ! gw temp tendency in the vert column
      usq(nhar,plev),              & ! U' amplitudes squared
      ed(plev),                    & ! Eddy mjh
      om(nhar)!,                    & ! observed freq at the source level

     ! Output variables. 
     REAL(prcn), INTENT(OUT) :: &
          drag_u(plev), drag_v(plev), temp_gw(plev), mk_eddy(plev)
     
     ! Local workspace and input variables
     integer  :: j,k         ! loop indices

     ! mjh latitude grid point
     integer :: lat

     ! Filtering variables. mjh
     real(prcn)   :: ubm, yv=1., xv=1.
     real(prcn)   :: factor=0.

     !-----------------------------------------------------------------------
     ! Determine the source layer wind and unit vectors, then project winds.
     !-----------------------------------------------------------------------

     ! Anisotropic spectrum
     IF(anisotropic == 1) THEN

        k   = levgws 
        ubm = sqrt (uu(k)**2 + vv(k)**2)
        IF(ubm == 0) ubm = 1.

        xv  = ABS(uu(k) / ubm)
        yv  = ABS(vv(k) / ubm)
        
        ! Project the local wind at midpoints onto the source wind.
        do k = levgws, plev+1-levgws
           ubm   = uu(k) * xv + vv(k) * yv
           uu(k) = ubm
           vv(k) = 0.  ! NS drag not calculated below for anisotropic
        end do

     ELSE
        
        xv = 1.
        yv = 1.
     ENDIF

     
     ! Initialize gravity wave drag tendencies to zero
     
     do k=1,plev
        drag_u(k)  = 0.
        drag_v(k)  = 0.
        temp_gw(k) = 0.
        mk_eddy(k) = 0.  ! mjh
     end do

     ! mjh
     factor = gwave_int_fac*fac(lat)

     !-------------------------------------------------------------------------
     ! Driver for mkgwd.
     !-------------------------------------------------------------------------
     
     do j=1,nhar                    ! Initialize squared amplitudes of
        ! ! the spectrum at the bottom, levgws,
        ! ! (where wave source is assumed)
        usq(j,levgws)=usql(j,lat,1) !
     end do                         !

     !+
     ! Eastward waves (c>0), ............ x-direction
     !-
     omeg(:,:)=0. ! mz_ab_20100901
     do j=1,nhar
        om(j) = bvfr(levgws) * gwkh / zkz(j)
        omeg(j,levgws) = om(j) - gwkh * uu(levgws)
     end do

     call mkgwd(uu, pf, ph, tf, th, bvfr, dpf, usq, omeg,                     &
          ax, ttpgw, scp, ed, svisc, plev, plev)

     Do k=1,plev                   ! accumulate GW tendency
        drag_u(k)  = drag_u(k)  + ax(k)
        temp_gw(k) = temp_gw(k) + ttpgw(k)
        mk_eddy(k) = mk_eddy(k) + ed(k)  ! mjh
     end do

     !+
     ! Westward waves (c<0), ............ x-direction
     !-
     do j=1,nhar
        omeg(j,levgws) = -om(j) - gwkh * uu(levgws)
     end do
     
     call mkgwd(uu, pf, ph, tf, th, bvfr, dpf, usq, omeg, &
          ax, ttpgw, scp, ed, svisc, plev, plev)

     do k=1,plev                   ! accumulate GW tendency
        drag_u(k)  = drag_u(k)  + ax(k)
        temp_gw(k) = temp_gw(k) + ttpgw(k)        
        mk_eddy(k) = mk_eddy(k) + ed(k)  ! mjh

        ! mjh
        drag_u(k)  = drag_u(k)*factor*xv
        temp_gw(k) = temp_gw(k)*factor*xv
        mk_eddy(k) = mk_eddy(k)*factor*xv

        if (abs(drag_u(k)) > tndmax) then
           drag_u(k) = sign(tndmax,drag_u(k))
        endif

        if (abs(temp_gw(k)) > dtempmax) then
           temp_gw(k) = sign(dtempmax,temp_gw(k))
        endif

        ! mjh
        if (abs(mk_eddy(k)) > eddymax) then
           mk_eddy(k) = sign(eddymax,mk_eddy(k))
        endif

     end do

     do j=1,nhar                    ! Initialize squared amplitudes of
        ! ! the spectrum at the bottom, levgws,
        ! ! (where wave source is assumed)
        usq(j,levgws)=usql(j,lat,2) !
     end do   

     ! Only do E/W for anisotropic spectrum mjh
     IF(anisotropic == 1) RETURN
     ! RETURN
     !+
     ! Northward waves (c>0), ............ y-direction
     !-
     do j=1,nhar
        omeg(j,levgws) = om(j) - gwkh * vv(levgws)
     end do
     
     call mkgwd(vv, pf, ph, tf, th, bvfr, dpf, usq, omeg,                     &
          ax, ttpgw, scp, ed, svisc, plev, plev)
     
     do k=1,plev                   ! accumulate GW tendency
        drag_v(k)  = drag_v(k)  + ax(k)
        temp_gw(k) = temp_gw(k) + ttpgw(k)
        mk_eddy(k) = mk_eddy(k) + ed(k)  ! mjh
     end do
     
     !+
     ! Southward waves (c<0), ............ y-direction
     !-
     do j=1,nhar
        omeg(j,levgws) = -om(j) - gwkh * uu(levgws)
     end do
     
     call mkgwd(vv, pf, ph, tf, th, bvfr, dpf, usq, omeg, &
          ax, ttpgw, scp, mk_eddy, svisc, plev, plev)
     

     do k=1,plev                   ! accumulate GW tendency
        
        drag_v(k)  = drag_v(k)  + ax(k)
        temp_gw(k) = temp_gw(k) + ttpgw(k)
        mk_eddy(k) = mk_eddy(k) + ed(k)  ! mjh

        ! mjh
        drag_v(k)  = drag_v(k)*factor*yv
        temp_gw(k) = temp_gw(k)*factor*yv
        mk_eddy(k) = mk_eddy(k)*factor*yv

        if (abs(drag_v(k)) > tndmax) then
           drag_v(k) = sign(tndmax,drag_v(k))
        endif
        
        if (abs(temp_gw(k)) > dtempmax) then
           temp_gw(k) = sign(dtempmax,temp_gw(k))
        endif

        ! mjh
        if (abs(mk_eddy(k)) > eddymax) then
           mk_eddy(k) = sign(eddymax,mk_eddy(k))
        endif


     end do

     return
   end subroutine mkgwintr
   
   
   
   subroutine mkgwinti (latitude, lat_dim) 
     !-----------------------------------------------------------------------   
     ! Time independent initialization for MK95 gravity wave parameterization.
     ! Isotropic GW source spectrum.        A. Medvedev   Oct 12/1999
     !-----------------------------------------------------------------------
     implicit none
     !-----------------------------------------------------------------------
     
     ! Input variables     
     REAL(prcn), DIMENSION(:), INTENT(IN)  :: latitude ! latitude (deg north)
     INTEGER, INTENT(IN)                   :: lat_dim

     ! Local variables 
     integer :: k,m
     real(prcn)   ::  ddkz(nhar) !,          & ! variable dkz to convert PSD into amplitudes

     !-----------------------------------------------------------------------

     !+
     ! Set MKGWD constants
     !-
     rg      = rd/g             ! gas constant/acceleration of gravity
     sqrt2   = sqrt(2.)
     sqrtp   = sqrt(pi)
     sqrt2p  = sqrt(2.*pi)
     sqrtkh  = sqrt2*gwkh

     !+
     ! Prescribe "universal" wave spectrum at the lower boundary 
     ! (appr. 165mb)
     !-
     dxkz = (log(xkzu) - log(xkzl)) / REAL(nhar-1,prcn)
     dxkz = exp(dxkz)
     
     zkz(1) = xkzu
     do k = 1, nhar-1
        zkz(k+1) = zkz(k) / dxkz
        ddkz(k)  = zkz(k) - zkz(k+1)
     end do
     ddkz(nhar) = zkz(nhar)*(dxkz - 1.)/dxkz
     
     ! latitudinal dep factor. mjh
     
     do m = 1, lat_dim

        ! Sinusoid Peak at mid latitudes
        fac(m) = 0.5*(1 + COS((ABS(latitude(m))-45.)/90.*pi))

        ! Vary spectral energy density with latitude
        do k = 1, nhar
           
           ! Zonal
           usql(k,m,1) = ddkz(k)*amp*fac(m)*                                  &
                zkz(k)/mstar /(1.+(zkz(k)/mstar)**4)
           
           ! Meridional
           !usql(k,m,2) = ddkz(k)*amp*amp_ratio*                              &
           !     zkz(k)/mstar /(1.+(zkz(k)/mstar)**4)
           usql(k,m,2) = usql(k,m,1) 

        end do

        ! Latitudinal variation in efficiency. Medvedev shouldn't need
        ! this, so set to 1.
        fac(m) = 1.

     enddo

     !+
     ! Other constants
     !-
     tndmax   = 500. / 86400. ! Wave drag shouldn't exceed 500 m/s/day
     dtempmax = 80./ 86400.   ! gw heat tendency shouldn't exceed 80 K/day
     eddymax  = 400.          ! Maximum allowed eddy diffusion coefficient. mjh
     
     ! mjh
     print *, ' '
     print *, 'Initialising Medvedev and Klassen gravity wave drag  ...'
     IF(anisotropic == 0) print *, 'mkgwinti:  Isotropic source spectrum'
     IF(anisotropic == 1) print *, 'mkgwinti:  Anisotropic source spectrum'
     print *, 'mkgwinti:  Defined at pressure = ',p0_gw,' mb'
     print *, 'mkgwinti:  nhar  = ', nhar
     print *, 'mkgwinti:  mstar = ', mstar, ' ampl = ', amp
     print *, 'mkgwinti:  Lz_min = ', 2*pi/xkzl, ' Lz_max = ', 2*pi/xkzu
     print *, 'mkgwinti:  kh = ', 2*pi/gwkh
     print *, 'mkgwinti:  Modified DeSaubius spectrum,  s = 1'
     

   end subroutine mkgwinti
   


   
   subroutine mkgwd(u, pf, ph, tf, th, bvfr, dp, usq, &
        omeg, ax, ttpgw, scp, ed, svisc, plev, lupp)

     ! Temperature tendency added.   A. Medvedev  (Oct 12/1999)
     
     !-----------------------------------------------------------------------
     implicit none
     

     ! Input variables

     INTEGER, INTENT(IN) :: plev, lupp

     real(prcn)::                 &
      bvfr(plev),                  & ! "full-level" Brunt-Vaisalla frequency
      dp(plev),                    & ! "full-level" delta p (=-dpm)
      omeg(nhar,plev),             & ! intrinsic frequency
      pf(plev),                    & ! "full-level" pressure
      ph(plev),                    & ! "half-level" pressure
      tf(plev),                    & ! "full-level" temperature
      th(plev),                    & ! "half-level" temperature
      usq(nhar,plev),              & ! U' amplitudes squared
      u(plev),                     & ! "full-level" background wind
      scp(plev),                   & ! Specific heat capacity. mjh
      ed(plev),                    & ! Local eddy
      svisc(plev)                    ! Coeff of turbulent+molecular viscosity

     ! Output variables

     real(prcn)::                 &
      ax(plev),                    & ! gwd drag in the vertical column
      ttpgw(plev)                    ! gw temp tend in the vert column

     ! Local workspace

     integer ::                    &
          k,l,lm,lp,nh,            & ! loop indices
          l1p

     real(prcn) :: a1,a11,         & 
          alpha(nhar),             &  ! parameter alpha for the spectrum
          alpha1

     real(prcn)::                  &
          bet,                     &
          axx(plev),               &
          beta(nhar,plev),         & ! nonlinear induced damping rate
          bvkh,                    &
          drag(nhar),              & ! GW drag
          dz, dz2,                 &
          gwkhu,                   & ! K_h dot U
          rms(plev),               & ! RMS^2(=variance) of the entire spect
          sigma(nhar),             & ! observed frequencies
          ttnd(plev),              & !
          wt(plev),                & ! workspace for gw temp tendency
          wt_k(nhar),              & !
          vv, vv2,                 &
          v2(nhar),                & ! velocity variance (=sigma(m)^2)
          x, xx,                   &
          zkzr(nhar,plev),         & ! vertical WN
          eddy_1(nhar)               ! Eddy stuff mjh
         

     ! machine zero. Actually a bit more than the real one, to catch
     ! zero*zero where zero=1.e<big number> mjh
     REAL(prcn), PARAMETER      :: zero = 1.e-30

     !-----------------------------------------------------------------------
     !+
     ! Initialize gravity wave drag tendency to zero
     !-
     l1p = levgws + 1
!     do l=l1p, lupp 
     do l=1, lupp ! mz_ab_20100830

        ax(l)    = zero
        axx(l)   = zero
        wt(l)    = zero
        ed(l)    = zero ! mjh
        ttnd(l)  = zero
        ttpgw(l) = zero
        rms(l)   = zero
     end do
     !+
     ! Compute parameters at the lower boundary (ie source level=levgws)
     !-
     gwkhu = gwkh * u(levgws)
     vv2   = zero

     do k=1,nhar

        beta(k,levgws) = zero
        sigma(k)       = omeg(k,levgws) + gwkhu
        zkzr(k,levgws) = gwkh * bvfr(levgws) / omeg(k,levgws)
        vv2            = vv2 + usq(k,levgws)
        v2(k)          = vv2
        vv             = sqrt(vv2)
        alpha(k)       = omeg(k,levgws) / sqrt2 / (gwkh*vv)
     end do
     !+
     ! Step-by-step integration upward ..................................
     !-
     do l=l1p, plev
        lm    = l - 1
        gwkhu = gwkh*u(l)
        xx    = tf(l)*pf(lm) / (tf(lm)*pf(l))

        dz    = -rg * dp(lm) * th(l) / ph(l)

        bvkh  = bvfr(l) * gwkh

        do nh=1,nhar

           !+
           ! Exclude harmonics which have been filtered out by critical levels
           !-
           if(usq(nh,lm) <= zero) then

              usq(nh,l)  = zero
              drag(nh)   = zero
              wt_k(nh)   = zero
              eddy_1(nh) = zero   ! mjh
           else
              !+
              ! Calculate only for harmonics with nonzero amplitude, 
              ! i.e. which haven't met critical levels yet
              !-
              vv         = sqrt(v2(nh))
              omeg(nh,l) = sigma(nh) - gwkhu
              alpha1     = omeg(nh,l) / (sqrtkh*vv)
              a1         = alpha1 * alpha(nh) 
              a11        = alpha1 * alpha1 

              if((a1 <= zero).OR. (a11 <= 0.5d0)) then  ! catch FPE mjh
                 usq(nh,l)  = zero
                 drag(nh)   = zero
                 wt_k(nh)   = zero
                 eddy_1(nh) = zero   ! mjh
              else
                 alpha(nh)  = alpha1
                 beta(nh,l) = sqrt2p * exp(-a11) * bvfr(l) / vv
                 zkzr(nh,l) = bvkh / omeg(nh,l)

                 ! Add in a damping term due to viscosity. mjh
                 beta(nh,l) = beta(nh,l) + visc_fac*ABS(zkzr(nh,l)*           &
                                                        svisc(l)*dz/ph(l))

                 x          = xx * zkzr(nh,l) / zkzr(nh,lm)
                 bet        = 0.5 * (beta(nh,l) + beta(nh,lm))
                 usq(nh,l)  = x * exp(-bet*dz) * usq(nh,lm)
                 !+
                 drag(nh)   = bet * usq(nh,l) * gwkh / zkzr(nh,l)
                 !+
                 ! Sensible heat flux: w'_kT'_k
                 !-
                 wt_k(nh)   = -0.5 * bvfr(l) * drag(nh) / zkzr(nh,l)          &
                      * tf(l) / g


                 ! Integrate eddy diffusion coeffient. Equation from thermal 
                 ! effects paper. (2*m*bvfr)-1*drag mjh
                 eddy_1(nh)  =  drag(nh)/(2.*zkzr(nh,l)*bvfr(l))

              endif
           endif
        end do
        
        vv2 = zero
        do nh=1,nhar
           axx(l) = axx(l) + drag(nh)
           wt (l) = wt (l) + wt_k(nh)
           ed(l)  = ed(l) + eddy_1(nh)  ! mjh
           vv2    = vv2 + usq(nh,l)
           v2(nh) = vv2
        end do
        rms(l)=v2(nhar)
        
     end do          ! ....end of the altitude loop.................

     do l=2, plev-1 
        lm = l - 1
        lp = l + 1
        !+
        ! (2 del z) filter on ax
        !-
        ax(l) = (axx(lm) + 2.*axx(l) + axx(lp)) * 0.25
        dz2   = -rg * (dp(lm) * th(l) /ph(l) &
             +dp(l) * th(lp) / ph(lp) )


        !+
        ! calculate temperature tendency
        !-
        ttnd(l) = wt(l) / rg / tf(l) - &
             (wt(lp) - wt(lm) ) / dz2
        ttnd(l) = ttnd(l) + 2.*abs(wt(l)) * g / tf(l) / scp(l)
        
     end do
     ax(1)=axx(1)
     ax(plev)=axx(plev)
     !+
     ! (2 del z) filter on ttpgw
     !-
     do l=2, plev-1
        ttpgw(l)=(ttnd(l-1)+2.*ttnd(l)+ttnd(l+1))*0.25
     end do
     ttpgw(1)=ttnd(1)
     ttpgw(plev)=ttnd(plev)
     
     ! 2 del z filter on eddies too. mjh
     do l=2, plev-1
        ed(l)=(ed(l-1)+2.*ed(l)+ed(l+1))*0.25
     end do

     return
   end subroutine mkgwd


    
   ! ****************************************************************
   !
   ! Wrapper around spline of parameters from cmat to Gwave grid. 
   !

   SUBROUTINE SplineOn(cmat_field, interp_field, field_flag, zero_winds,       &
       placemat, placegwave, ht_dim, lev, status)

     USE messy_main_tools, ONLY: Spline1d, Splint1d

     IMPLICIT NONE

    INTEGER,    INTENT(IN)  :: ht_dim
    INTEGER,    INTENT(IN)  :: lev
    INTEGER,    INTENT(IN)  :: field_flag
    INTEGER,    INTENT(IN)  :: zero_winds
    REAL(prcn), INTENT(IN)  :: cmat_field(ht_dim) 
    REAL(prcn), INTENT(IN)  :: placemat(ht_dim)
    REAL(prcn), INTENT(IN)  :: placegwave(lev)

    REAL(prcn), INTENT(OUT) :: interp_field(lev)
    INTEGER,    INTENT(OUT) :: status


    ! arrays for the interpolation
    REAL(prcn) :: splineno(ht_dim),splineno11(ht_dim),splineno1(lev)

    INTEGER    :: k

    ! first zonal winds...(bottom up)
    DO k=1,ht_dim
       splineno(k)=cmat_field(ht_dim+1-k)
    ENDDO
    
    CALL Spline1d(placemat,splineno,ht_dim,0._prcn,0._prcn,splineno11,.TRUE.)
    
    DO k=1,lev
       CALL Splint1d(placemat,splineno,splineno11,ht_dim,placegwave(k),       &
            splineno1(k),status)
    ENDDO
    
    ! We make code handle any altitude range. mjh
    DO k=1,lev

       IF(placemat(ht_dim) >= placegwave(k)) THEN
          interp_field(lev+1-k) = splineno1(k)
       ELSE

          ! No gradient or zero below lower b
          interp_field(lev+1-k) = interp_field(lev+1-k +1)
          
          ! It's a velocity field and zero_winds is on
          IF(zero_winds == 1 .AND. field_flag == 1) interp_field(lev+1-k) = 0.

          ! It's a height field
          IF(field_flag == 2) interp_field(lev+1-k) = interp_field(lev+1-k +1) &
               - (interp_field(lev+1-k +2) - interp_field(lev+1-k +1))

       ENDIF
    ENDDO

    
  END SUBROUTINE SplineOn
 ! ****************************************************************



  ! ****************************************************************
  SUBROUTINE mk_read_nml_ctrl(status, iou, modstr)

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    NAMELIST /CTRL_MK/ gwave_int_fac

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER       :: substr='mk_read_nml_ctrl'
    LOGICAL                           :: lex          ! file exists ?
    INTEGER                           :: fstat        ! file status
    CHARACTER(LEN=*), INTENT(IN) :: modstr

    ! INITIALIZE
    status = 1 ! ERROR

    CALL read_nml_open(lex, substr, iou, 'CTRL_MK', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL_MK, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL_MK', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST
    WRITE(*,*) ' Gwave intermittency fac  : ', gwave_int_fac

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE mk_read_nml_ctrl
  ! ****************************************************************


END MODULE messy_gwave_mk

