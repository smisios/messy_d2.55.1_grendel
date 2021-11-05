! **********************************************************************

! SUBMODEL CORE LAYER (SMCL) ROUTINES FOR MESSy SUBMODEL RELAX
!
! THIS SUBMODEL IS USED TO RELAX VARIABLES TO A EQUILIBRIUM VALUE. 
! IT CONTAINS A) NEWTONIAN COOLING, ...
! I.E. RELAXES TEMPERATURE AGAINST A SPECIFIED BACKGROUND TEMPERATURE
! AND B) RAYLEIGH FRICTION,
! I.E. RELAXES HORIZONTAL WINDS AGAINST ZERO
!
!
! Author : Hella Garrny,     DLR-IPA, June 2016, original code
!          Roland Walz,      DLR-IPA, 2018, equilibrium temperature
!                                         , idealized heating for planetary wave generation
!                                     2019, climate change-like tropical upper-tropospheric heating
!                                     2019, EH sponge for vertical grids other
!                                           than L90MA
!                                     2020, tropical upper-tropospheric heating with log-p 
!                                           instead of linear p inside gaussian heating function
!          Matthias Nuetzel, DLR-IPA, 2018, (localized) idealized heating for monsoon
!
! References:
! Held & Suarez (1994), A proposal for the intercomparison of the dynamical cores of atmospheric general 
!   circulation models, BAMS, Vol. 75, p. 1825-1830
! Polvani & Kushner (2002), Tropospheric response to stratospheric perturbations in a relatively simple 
!   general circulation model, GRL, Vol. 29, No. 7, p. 10.1029/2001GL014284
! Butler, Thompson & Heikes (2010), The Steady-State Atmospheric Circulation Response to Climate Change-like
!   Thermal Forcings in a Simple General Circulation Model, Journal of Climate, Vol. 23, No.13, p. 3474-3496,
!   https://doi.org/10.1175/2010JCLI3228.1
! Lindgren, Sheshadri & Plumb (2018), Sudden stratospheric warming formation in an idealized general
!   circulation model using three types of tropospheric forcing,
!   Journal of Geophysical Research: Atmospheres, 123, https://doi.org/10.1029/2018JD028537
! Schubert, W. H., and M. T. Masarik, 2006: Potential vorticity aspects of the MJO.
!   Dynam. Atmos. Oceans, 42, 127-151. Doi: 10.1016/j.dynatmoce.2006.02.003.
!   https://doi.org/10.1016/j.dynatmoce.2006.02.003
!
!
!
! **********************************************************************

! **********************************************************************
MODULE messy_relax
! **********************************************************************

  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: DP

  ! ----------- <

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'relax'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '1.1'
 
  ! PUBLIC SUBROUTINES (to be called from messy_relax_e5.f90)
  ! ### add your own public subroutines here
  PUBLIC :: relax_rayfr_smcl
  PUBLIC :: relax_newco_smcl
  ! functions for TEQU, KAPPA AND KDAMP
  PUBLIC :: relax_kdamphs
  PUBLIC :: relax_kdamppk
  PUBLIC :: relax_kdampeh
  PUBLIC :: relax_tequhs
  PUBLIC :: relax_kappahs
  PUBLIC :: relax_tequpk
  PUBLIC :: relax_intpol_p
  
  ! op_rw_20190212
  ! function for climate change-like tropical upper-tropospheric heating
  PUBLIC :: relax_tteh_cc_tropics_smcl
 
  ! op_rw_20181102
  ! function for planetary wave generation by idealized heating
  PUBLIC :: relax_tteh_waves_smcl 

  ! op_mn_20180620
  ! function for monsoon generation by idealized heating
  PUBLIC :: relax_tteh_mons_smcl


  ! PRIVATE SUBROUTINES
  ! ### add your own private subroutines here

CONTAINS

  ! =========================================================================
  ! ### add your own public subroutines here
  ! =========================================================================

  ! =========================================================================
  SUBROUTINE relax_rayfr_smcl(u,v,kdamp,ute,vte)

    IMPLICIT NONE
    INTRINSIC :: SIZE

    ! I/O
    REAL(DP), DIMENSION(:,:), INTENT(OUT) :: ute,vte  ! wind tendencies
    REAL(DP), DIMENSION(:,:), INTENT(IN)  :: u,v      ! winds
    REAL(DP), DIMENSION(:,:), INTENT(IN)  :: kdamp    ! damping coefficient

    
    ! LOCAL
    INTEGER :: jp, jk, lkproma, nlev

    lkproma = SIZE(u,1)
    nlev   = SIZE(u,2)

    ! INIT
    ute(:,:) = 0.0_dp
    vte(:,:) = 0.0_dp

    !calculate tendency
    DO jp = 1, lkproma
        DO jk = 1, nlev
            ute(jp,jk) = -kdamp(jp,jk)*u(jp,jk)
            vte(jp,jk) = -kdamp(jp,jk)*v(jp,jk)
        END DO
    END DO

  END SUBROUTINE relax_rayfr_smcl
  ! =========================================================================

  ! =========================================================================
  SUBROUTINE relax_newco_smcl(temp, tequ, kappa, tte1)

    IMPLICIT NONE
    INTRINSIC :: SIZE

    ! I/O
    REAL(DP), DIMENSION(:,:), INTENT(OUT) :: tte1  ! perturbation temperature
    REAL(DP), DIMENSION(:,:), INTENT(IN)  :: temp  ! temperature
    REAL(DP), DIMENSION(:,:), INTENT(IN)  :: tequ  ! equilibrium temperature
    REAL(DP), DIMENSION(:,:), INTENT(IN)  :: kappa ! inverse relaxation time scale


    ! LOCAL
    INTEGER :: jp, jk, lkproma, nlev

    lkproma = SIZE(temp,1)
    nlev   = SIZE(temp,2)

    ! INIT
    tte1(:,:) = 0.0_dp

    !calculate tendency
    DO jp = 1, lkproma
        DO jk = 1, nlev
            tte1(jp,jk) = -kappa(jp,jk)*(temp(jp,jk)-tequ(jp,jk))
        END DO
    END DO

  END SUBROUTINE relax_newco_smcl
  ! =========================================================================

  ! =========================================================================
  SUBROUTINE relax_kdamphs(kdamp1,p,p_surf,kmaxHS,sig0)

    ! calculate relaxation damping coefficients for drag close to surface (p/ps <= 0.7) according
    ! to formula by Held & Suarez, 1994, BAMS

    IMPLICIT NONE
    INTRINSIC :: SIZE, max

    ! I/O
    REAL(DP), DIMENSION(:,:), INTENT(OUT) :: kdamp1  ! perturbation temperature
    REAL(DP), DIMENSION(:,:), INTENT(IN)  :: p       ! pressure
    REAL(DP), DIMENSION(:),   INTENT(IN)  :: p_surf  ! surface pressure
    REAL(DP),                 INTENT(IN)  :: kmaxHS  ! maximum wind damping at surface
    REAL(DP),                 INTENT(IN)  :: sig0    ! sigma level at which surface wind damping stops


    ! LOCAL
    INTEGER  :: jp, jk, lkproma, nlev

    lkproma = SIZE(kdamp1,1)
    nlev    = SIZE(kdamp1,2)

    ! INIT
    kdamp1(:,:) = 0.0_dp

    !calculate damp coeff as fct(p)
    DO jp = 1, lkproma
        DO jk = 1, nlev
            kdamp1(jp,jk) = kmaxHS*max(0.0_dp,(p(jp,jk)/p_surf(jp)-sig0)/(1.0_dp-sig0))
        END DO
    END DO


  END SUBROUTINE relax_kdamphs
  ! =========================================================================

  ! =========================================================================
  SUBROUTINE relax_kdamppk(kdamp1,p,kmaxPK,psp)

    ! calculate relaxation damping coefficients for drag at layer at model top according
    ! to formula by Polvani & Kushner, GRL, 2002

    IMPLICIT NONE
    INTRINSIC :: SIZE

    ! I/O
    REAL(DP), DIMENSION(:,:), INTENT(OUT) :: kdamp1   ! perturbation temperature
    REAL(DP), DIMENSION(:,:), INTENT(IN)  :: p        ! pressure
    REAL(DP),                 INTENT(IN)  :: kmaxPK   ! approximate wind damping at model top
    REAL(DP),                 INTENT(IN)  :: psp      ! sponge layer (damping above this layer)


    ! LOCAL
    INTEGER  :: jp, jk, lkproma, nlev

    lkproma = SIZE(kdamp1,1)
    nlev    = SIZE(kdamp1,2)

    ! INIT
    kdamp1(:,:) = 0.0_dp

    !calculate damp coeff as fct(p)
    DO jp = 1, lkproma
        DO jk = 1, nlev
            IF (p(jp,jk) .le. psp) THEN
                kdamp1(jp,jk) = kmaxPK*(1.0_dp - p(jp,jk)/psp)**2
            ENDIF
        END DO
    END DO


  END SUBROUTINE relax_kdamppk
  ! =========================================================================

  ! =========================================================================
  SUBROUTINE relax_kdampeh(kdamp1,pressure,spdrag,enfac)

    ! calculate relaxation damping coefficients for drag at layer at model top
    ! according
    ! to echam-like sponge

    IMPLICIT NONE
    INTRINSIC :: SIZE, log

    ! I/O
    REAL(DP), DIMENSION(:,:), INTENT(OUT) :: kdamp1     ! damping coefficient
    REAL(DP), DIMENSION(:,:), INTENT(IN)  :: pressure   ! pressure on model levels
    REAL(DP),                 INTENT(IN)  :: spdrag     ! damping prefactor
    REAL(DP),                 INTENT(IN)  :: enfac      ! enhancement factor (enspodi)

    ! LOCAL
    INTEGER             :: jk, jk90        ! loop variables
    INTEGER             :: nlevs           ! total number of pressure levels of grid
    INTEGER             :: nlevs_sp        ! number of sponge-layer levels of grid (if grid is not L90MA)
    INTEGER, PARAMETER  :: nlevsL90sp = 10 ! number of sponge-layer levels of L90MA-grid
    REAL(DP)            :: psp             ! current sponge-layer pressure in loop of kdamp1 calculation

    ! sponge-layer pressure levels of L90MA-grid (pressures smaller than 50 Pa)
    REAL(DP), DIMENSION(nlevsL90sp), PARAMETER  :: psp_L90 &
                                                   = (/0.9945911_dp, 3.182691_dp, 5.808412_dp, 8.959276_dp &
                                                    , 12.74032_dp,  17.17444_dp, 22.27369_dp, 28.13782_dp  &
                                                    , 34.88157_dp,  42.63688_dp /)

    ! damping coefficient of L90MA-grid
    REAL(DP), DIMENSION(nlevsL90sp)             :: kdamp_L90 

    ! First of all compute kdamp_L90
    DO jk = 1, nlevsL90sp
        kdamp_L90(jk) = spdrag*(enfac)**(nlevsL90sp-jk)
        !print *, '#*# kdamp_L90(', jk, ') = ', kdamp_L90(jk)
    END DO   

    ! INIT
    kdamp1(:,:) = 0.0_dp

    nlevs    = SIZE(kdamp1,2)   ! total number of pressure levels

    IF (nlevs .EQ. 90) THEN
        ! L90MA grid:
        DO jk = 1, nlevsL90sp
            kdamp1(:,jk) = kdamp_L90(jk)
        END DO
    ELSE
        ! Other vertical grids: calculate damp coeff from interpolation
        
        ! Count number of sponge-layer levels
        nlevs_sp = 0
        jk       = 1

        DO WHILE ( (pressure(1,jk) .LE. 50.0_dp) .AND. (jk .LE. nlevs) )
            jk   = jk + 1
        END DO

        nlevs_sp = jk - 1

        ! In the while-loop above it is justified to use pressures
        ! at any value for lkproma (here, the first one is used)
        ! since pressures of the grids L90MA and L47MA
        ! are independent of the surface pressure for
        ! levels 1 to 54 and for levels 1 to 20, respectively.
        ! L90MA: nlevs_sp = 10
        ! L47MA: nlevs_sp =  5


        IF (nlevs_sp .GT. 0) THEN

            ! Compute kdamp1
            DO jk = 1, nlevs_sp

                ! Assign sponge-layer pressures
                psp = pressure(1,jk)

                ! Get nearest neighbors of L90MA pressure levels
                jk90 = 1

                DO WHILE ( (psp .GT. psp_L90(jk90)) .AND. (jk90 .LE. nlevsL90sp) )
                    jk90 = jk90 + 1
                END DO

                ! Case of psp smaller than smallest psp_L90 pressure
                IF ( jk90 == 1 ) THEN
                    kdamp1(:,jk) = kdamp_L90(1)
            
                ! Case of psp larger than largest psp_L90 pressure
                ELSE IF ( jk90 .GT. nlevsL90sp ) THEN
                    kdamp1(:,jk) = kdamp_L90(nlevsL90sp)
            
                ! Other cases: psp_L90(jk90 - 1) < psp < psp_L90(jk90)
                ELSE
                    kdamp1(:,jk) = kdamp_L90(jk90) + ( kdamp_L90(jk90 - 1) - kdamp_L90(jk90) ) &
                                                 / log(  psp_L90(jk90 - 1) /   psp_L90(jk90) ) &
                                                 * log(  psp               /   psp_L90(jk90) )
            
                END IF

                !print *, '#*# kdamp1(1,', jk, ') = ', kdamp1(1,jk)
 
            END DO

        END IF

    END IF


  END SUBROUTINE relax_kdampeh
  ! =========================================================================


  ! =========================================================================
  SUBROUTINE relax_tequhs(tequ1,lat,p,hfac,p0,T0,T1,Ty,Tz,eps_abs)

    ! calculate equilibrium temperature and inverse relaxation time scale according 
    ! to formula by Held & Suarez, 1994, BAMS
    ! which originally is hemispherically symmetric.
    ! But you can choose a winter hemisphere via hfac and an asymmetry via eps_abs as well.

    USE messy_main_constants_mem, ONLY: pi

    IMPLICIT NONE
    INTRINSIC :: SIZE, sin, cos, log, max, sign

    ! I/O
    REAL(DP), DIMENSION(:,:), INTENT(OUT) :: tequ1   ! perturbation temperature
    REAL(DP), DIMENSION(:,:), INTENT(IN)  :: p       ! pressure
    REAL(DP), DIMENSION(:),   INTENT(IN)  :: lat     ! latitude
    REAL(DP),                 INTENT(IN)  :: hfac    ! hemispheric factor
    REAL(DP),                 INTENT(IN)  :: p0      ! reference pressure
    REAL(DP),                 INTENT(IN)  :: T0      ! minimum temperature
    REAL(DP),                 INTENT(IN)  :: T1      ! maximum temperature
    REAL(DP),                 INTENT(IN)  :: Ty      ! meridional temperature gradient in troposphere
    REAL(DP),                 INTENT(IN)  :: Tz      ! vertical temperature gradient in troposphere
    REAL(DP),                 INTENT(IN)  :: eps_abs ! absolute value of asymmetry factor in troposphere


    ! LOCAL
    REAL(DP), PARAMETER     ::  k =  0.2857_dp

    INTEGER  :: jp, jk, lkproma, nlev
    REAL(DP) :: sinlat, coslat, eps

    lkproma = SIZE(tequ1,1)
    nlev    = SIZE(tequ1,2)

    IF (hfac==0.0_dp) THEN
        eps = 0.0_dp
    ELSE
        eps = sign(eps_abs,hfac)
    END IF

    !calculate tendency
    DO jp = 1, lkproma
        sinlat = sin(lat(jp)/180.0_dp*pi)
        coslat = cos(lat(jp)/180.0_dp*pi)
        DO jk = 1, nlev
            tequ1(jp,jk) = max(T0, (T1 - Ty*sinlat**2 - &
                 Tz*log(p(jp,jk)/p0)*coslat**2 - &
                 eps*sinlat)*(p(jp,jk)/p0)**k)
        END DO
    END DO


  END SUBROUTINE relax_tequhs
  ! =========================================================================

  ! =========================================================================
  SUBROUTINE relax_kappahs(kappa1,lat,p,p_surf,ta,ts,sigb)

    ! calculate equilibrium temperature and inverse relaxation time scale according 
    ! to formula by Held & Suarez, 1994, BAMS

    USE messy_main_constants_mem, ONLY: pi

    IMPLICIT NONE
    INTRINSIC :: SIZE, cos, max

    ! I/O
    REAL(DP), DIMENSION(:,:), INTENT(OUT) :: kappa1 ! inverse relaxaion time scale
    REAL(DP), DIMENSION(:),   INTENT(IN)  :: lat    ! latitude
    REAL(DP), DIMENSION(:,:), INTENT(IN)  :: p      ! pressure on model level
    REAL(DP), DIMENSION(:),   INTENT(IN)  :: p_surf ! surface pressure
    REAL(DP),                 INTENT(IN)  :: ta     ! relaxation time outside of trop. troposphere
    REAL(DP),                 INTENT(IN)  :: ts     ! relaxation time at surface of trop. troposphere
    REAL(DP),                 INTENT(IN)  :: sigb   ! sigma level below which the short relax time is used



    ! LOCAL
    INTEGER  :: jp, jk, lkproma, nlev
    REAL(DP) :: coslat, ka, ks

    lkproma = SIZE(kappa1,1)
    nlev   = SIZE(kappa1,2)

    ! convert to inverse relaxation time in 1/day
    ka = 1.0_dp/(ta*24.0_dp*60.0_dp*60.0_dp)
    ks = 1.0_dp/(ts*24.0_dp*60.0_dp*60.0_dp)

    !calculate tendency
    DO jp = 1, lkproma
        coslat = cos(lat(jp)/180._dp*pi)
        DO jk = 1, nlev
            kappa1(jp,jk)= ka + (ks-ka)*max(0.0_dp,(p(jp,jk)/p_surf(jp)-sigb)/(1.0_dp - sigb))*coslat**4
        END DO
    END DO


  END SUBROUTINE relax_kappahs
  ! =========================================================================

  ! =========================================================================
  SUBROUTINE relax_tequpk(tequ1,lat,p,gammapk,hfac,p0,T1,Ty,Tz,eps_abs,l0_abs,dl,pT_SH,pT_WH,no_polar_vortex)

    ! calculate equilibrium temperature according to
    ! modified formula by Polvani & Kushner, GRL, 2002
    ! The value of several parameters is set
    ! via the namelist 

    USE messy_main_constants_mem, ONLY: pi

    IMPLICIT NONE
    INTRINSIC :: SIZE, sin, cos, log, tanh, max, sign

    ! I/O
    REAL(DP), DIMENSION(:,:), INTENT(OUT) :: tequ1     ! perturbation temperature
    REAL(DP), DIMENSION(:),   INTENT(IN)  :: lat       ! latitude
    REAL(DP), DIMENSION(:,:), INTENT(IN)  :: p         ! pressure on model level
    REAL(DP),                 INTENT(IN)  :: gammapk   ! gamma value in PK02 (in K/km) 
    REAL(DP),                 INTENT(IN)  :: hfac      ! hemispheric factor
    REAL(DP),                 INTENT(IN)  :: p0        ! surface pressure
    REAL(DP),                 INTENT(IN)  :: T1        ! max. temp. for troposphere
    REAL(DP),                 INTENT(IN)  :: Ty        ! meridional temp. gradient in troposphere
    REAL(DP),                 INTENT(IN)  :: Tz        ! vertical temp. gradient in troposphere
    REAL(DP),                 INTENT(IN)  :: eps_abs   ! shift of max. temp. to summer hemisphere
    REAL(DP),                 INTENT(IN)  :: l0_abs    ! transition latitude where polar vortex starts
    REAL(DP),                 INTENT(IN)  :: dl        ! rapidity of transition
    REAL(DP),                 INTENT(IN)  :: pT_SH     ! tropopause pressure in summer hemisphere
    REAL(DP),                 INTENT(IN)  :: pT_WH     ! tropopause pressure in winter hemisphere
    LOGICAL,                  INTENT(IN)  :: no_polar_vortex ! if .TRUE., sets W_function to zero


    ! Constants
    REAL(DP),               PARAMETER ::  R = 287.0_dp  ! specific gas constant for dry air
    REAL(DP),               PARAMETER ::  g = 9.81_dp   ! gravitational acceleration
    REAL(DP),               PARAMETER ::  k = 0.2857_dp ! R / c_p
    REAL(DP),               PARAMETER :: ph = 0.1_dp    ! top boundary pressure
        
    ! Values of the US standard atmosphere that are used as background in the stratosphere
    REAL(DP), DIMENSION(9), PARAMETER :: pUS = (/101325.0_dp, 22632.1_dp, 5474.89_dp, 868.019_dp, 110.906_dp &
                                                            , 66.9389_dp, 3.95642_dp, 0.39814_dp, 0.0999_dp/)
    REAL(DP), DIMENSION(9), PARAMETER :: TUS = (/288.15_dp,    216.65_dp,  216.65_dp,  228.65_dp,  270.65_dp &
                                                            ,  270.65_dp, 214.65_dp, 187.65_dp, 187.65_dp /)

    INTEGER  :: jp, jk, lkproma, nlev, iUS
    REAL(DP) :: lat_current, sinlat, coslat, TT, Tint, TPV, eps, l0, p_current, p_tropopause, W_function

    lkproma = SIZE(tequ1,1)             ! number of latitudes
    nlev    = SIZE(tequ1,2)             ! number of levels

    
    ! Set eps and l0 (in the namelist, only eps_abs and l0_abs are set).
    eps = sign(1.0_dp,hfac) * eps_abs
    l0  = sign(1.0_dp,hfac) * l0_abs
    
    
    
    ! Temperature from TUS at pT_SH
    iUS = 0

    DO WHILE (pT_SH .LT. pUS(iUS+1) .AND. iUS .LT. 9)
        
        iUS = iUS + 1
        
    ENDDO    
    ! now the temperature TT = TUS(pT_SH) lies between TUS(iUS) and TUS(iUS+1)

    TT =  TUS(iUS) + (TUS(iUS+1) - TUS(iUS))/log(pUS(iUS+1)/pUS(iUS))*log(pT_SH/pUS(iUS))

    ! calculate tendency
    ! 1. loop over latitudes
    DO jp = 1, lkproma
    
        ! often needed values
        lat_current = lat(jp)
        sinlat = sin(lat_current/180._dp*pi)
        coslat = cos(lat_current/180._dp*pi)

        ! variable tropopause pressure
        p_tropopause = ( ( pT_WH - pT_SH ) / 2.0_dp ) * ( 1.0_dp + sign(1.0_dp,hfac) * tanh( (lat_current - l0)/dl ) ) + pT_SH;

        ! weight function for stratospheric temperature
        W_function = 0.5_dp * ( 1.0_dp + sign(1.0_dp,hfac) * tanh( (lat_current - l0)/dl ) )

        IF ( no_polar_vortex ) THEN
            W_function = 0.0_dp
        END IF

        ! 2. loop over levels
        DO jk = 1, nlev

            ! current pressure
            p_current = p(jp,jk)

            IF (p_current .LT. p_tropopause) THEN
            ! 1. stratospheric temperature
            ! -------------------------------------------------------------------------

                iUS = 0

                DO WHILE (p_current .LT. pUS(iUS+1) .AND. iUS .LT. 9)
                    iUS = iUS + 1
                ENDDO

                Tint = TUS(iUS) + (TUS(iUS+1) - TUS(iUS))/log(pUS(iUS+1)/pUS(iUS))*log(p_current/pUS(iUS))
            ! -------------------------------------------------------------------------

                IF (p_current .le. ph) THEN

                    ! if pressure is lower than some boundary pressure
                    TPV = TT*(ph/p_tropopause)**(-R*gammapk*(-1.0_dp*1e-3)/g)

                ELSE

                    TPV = TT*(p_current/p_tropopause)**(-R*gammapk*(-1.0_dp*1e-3)/g)

                END IF
            ! -------------------------------------------------------------------------

                tequ1(jp,jk) = ( 1.0_dp - W_function )*Tint + W_function*TPV

            ELSE
            ! 2. tropospheric temperature

                tequ1(jp,jk) = max(TT, (T1 - Ty*sinlat**2 - &
                     Tz*log(p_current/p0)*coslat**2 - &
                     eps*sinlat)*(p_current/p0)**k)

            ENDIF

        END DO

    END DO

  END SUBROUTINE relax_tequpk
  ! =========================================================================

  ! =========================================================================
  SUBROUTINE relax_intpol_p(tequin,p,hyam,hybm,tequout)

    ! interpolate tequ profile from import to actual pressure profile 
    ! of current time step

    USE messy_main_tools, ONLY: iso2ind, ind2val

    IMPLICIT NONE
    INTRINSIC :: SIZE

    ! I/O
    REAL(DP), DIMENSION(:,:), INTENT(IN)  :: tequin   ! temperature in
    REAL(DP), DIMENSION(:,:), INTENT(IN)  :: p        ! pressure on model level
    REAL(DP), DIMENSION(:),   INTENT(IN)  :: hyam, hybm ! hybrid level coefficients
    REAL(DP), DIMENSION(:,:), INTENT(OUT) :: tequout  ! temperature out

    ! LOCAL
    REAL(DP), PARAMETER     :: p0 = 101325._dp

    INTEGER  :: jp, jk, lkproma, nlev, k
    REAL(DP) :: f
    REAL(DP), DIMENSION(:), ALLOCATABLE   :: pin

    lkproma = SIZE(tequin,1)
    nlev    = SIZE(tequin,2)

    ALLOCATE(pin(nlev))

    pin(:) = hyam(:) + hybm(:)*p0
    DO jp = 1, lkproma
        DO jk = 1, nlev
            CALL iso2ind(pin/p0,p(jp,jk)/p0, k, f)
            CALL ind2val(tequout(jp,jk), tequin(jp,:), k, f)
        ENDDO
    ENDDO 

    DEALLOCATE(pin)

  END SUBROUTINE relax_intpol_p



  ! +op_rw_20190212
  ! =========================================================================
  SUBROUTINE relax_tteh_cc_tropics_smcl(pressure, p_surf, lat, q0_cct, lat0, sigma_lat, z0, sigma_z, Butler_heat, tteh_cc_tropics1)

    ! Calculate a temperature tendency due to idealized heating as in first line of Table 1 of
    ! Butler, Thompson & Heikes (2010), The Steady-State Atmospheric Circulation Response to Climate Change-like
    !   Thermal Forcings in a Simple General Circulation Model, Journal of Climate, Vol. 23, No.13, p. 3474-3496,
    !   https://doi.org/10.1175/2010JCLI3228.1
    ! if Butler_heat is set to True
    ! Otherwise, use a different shape with log-p instead of linear pressure inside gaussian heating function

    USE messy_main_constants_mem, ONLY: pi

    IMPLICIT NONE
    INTRINSIC :: SIZE, exp, log

    ! I/O
    REAL(DP), DIMENSION(:,:), INTENT(IN)  :: pressure         ! pressure on model levels
    REAL(DP), DIMENSION(:),   INTENT(IN)  :: p_surf           ! surface pressure
    REAL(DP), DIMENSION(:),   INTENT(IN)  :: lat              ! latitudes
    REAL(DP),                 INTENT(IN)  :: q0_cct           ! heating amplitude                 [K/day]
    REAL(DP),                 INTENT(IN)  :: lat0             ! latitudinal center of heating     [degree]
    REAL(DP),                 INTENT(IN)  :: sigma_lat        ! latitudinal half width of heating [rad]
    REAL(DP),                 INTENT(IN)  :: z0               ! sigma level center of heating     [1]
    REAL(DP),                 INTENT(IN)  :: sigma_z          ! vertical half width of heating    [1]
    LOGICAL,                  INTENT(IN)  :: Butler_heat      ! if .FALSE., use exp(-log(p)**2) decay instead of Butler heating
    REAL(DP), DIMENSION(:,:), INTENT(OUT) :: tteh_cc_tropics1 ! temperature tendency due to heating for wave generation

    ! Parameters
    REAL(DP), PARAMETER                   :: daytosec = 86400.0_dp

    ! Local
    INTEGER  :: jp, jk, lkproma, nlev
    REAL(DP) :: q0_cct_Ks                                     ! heating amplitude in K/s

    ! Init
    tteh_cc_tropics1(:,:) = 0.0_dp

    ! number of current grid points
    lkproma = SIZE(tteh_cc_tropics1,1)                        ! number of latitudes*longitudes
    nlev    = SIZE(tteh_cc_tropics1,2)                        ! number of levels

    ! heating amplitude in Kelvin per second
    q0_cct_Ks = q0_cct / daytosec


    IF ( Butler_heat ) THEN

        DO jp = 1, lkproma
            DO jk = 1, nlev

                tteh_cc_tropics1(jp,jk) = q0_cct_Ks * exp( -0.5_dp * ( ( (lat(jp) - lat0)*pi/180.0_dp ) / sigma_lat )**2 &
                                                       -0.5_dp * ( (pressure(jp,jk)/p_surf(jp) - z0)/ sigma_z   )**2 )

            END DO
        END DO

    ELSE

        DO jp = 1, lkproma
            DO jk = 1, nlev

                tteh_cc_tropics1(jp,jk) = q0_cct_Ks * exp( -0.5_dp * ( ( (lat(jp) - lat0)*pi/180.0_dp ) / sigma_lat )**2 &
                                                       -0.5_dp * ( ( log(pressure(jp,jk)/p_surf(jp)) - log(z0) )/ sigma_z   )**2 )

            END DO
        END DO

    END IF


  END SUBROUTINE relax_tteh_cc_tropics_smcl
  ! =========================================================================
  ! -op_rw_20190212



  ! +op_rw_20181102
  ! =========================================================================
  SUBROUTINE relax_tteh_waves_smcl(pressure, lat, lon, q0, m_WN, phi0, sigma_phi, p_bot, p_top, tteh_waves1)

    ! Calculate a temperature tendency due to idealized heating as in equation (1) of
    ! Lindgren, Sheshadri & Plumb (2018), Sudden stratospheric warming formation in an idealized general
    ! circulation model using three types of tropospheric forcing,
    ! Journal of Geophysical Research: Atmospheres, 123, https://doi.org/10.1029/2018JD028537

    USE messy_main_constants_mem, ONLY: pi

    IMPLICIT NONE
    INTRINSIC :: SIZE, sin, exp, log

    ! I/O
    REAL(DP), DIMENSION(:,:), INTENT(IN)  :: pressure    ! pressure on model levels
    REAL(DP), DIMENSION(:),   INTENT(IN)  :: lat         ! latitudes
    REAL(DP), DIMENSION(:),   INTENT(IN)  :: lon         ! longitudes
    REAL(DP),                 INTENT(IN)  :: q0          ! heating amplitude             [K/day]
    REAL(DP),                 INTENT(IN)  :: m_WN        ! longitudinal wave number      [1]
    REAL(DP),                 INTENT(IN)  :: phi0        ! latitudinal center of heating [degree]
    REAL(DP),                 INTENT(IN)  :: sigma_phi   ! latitudinal decay rate        [rad]
    REAL(DP),                 INTENT(IN)  :: p_bot       ! bottom pressure boundary      [Pa]
    REAL(DP),                 INTENT(IN)  :: p_top       ! top pressure boundary         [Pa]
    REAL(DP), DIMENSION(:,:), INTENT(OUT) :: tteh_waves1 ! temperature tendency due to heating for wave generation

    ! Parameters
    REAL(DP), PARAMETER                   :: daytosec = 86400.0_dp

    ! Local
    INTEGER  :: jp, jk, lkproma, nlev
    REAL(DP) :: q0_Ks                        ! heating amplitude in K/s
    REAL(DP) :: press                        ! current pressure

    ! Init
    tteh_waves1(:,:) = 0.0_dp

    ! number of current grid points
    lkproma = SIZE(tteh_waves1,1)            ! number of latitudes*longitudes
    nlev    = SIZE(tteh_waves1,2)            ! number of levels

    ! heating amplitude in Kelvin per second
    q0_KS = q0 / daytosec


    DO jp = 1, lkproma
        DO jk = 1, nlev

            press = pressure(jp,jk)

            IF ( (press .LE. p_bot) .AND. (press .GE. p_top) ) THEN

                tteh_waves1(jp,jk) = q0_Ks * sin( m_WN * lon(jp) * pi / 180.0_dp ) &
                                           * exp( -0.5_dp * ( (lat(jp)-phi0) / (sigma_phi*180.0_dp/pi) )**2 ) &
                                           * sin( pi * log(press/p_bot) / log(p_top/p_bot) )

            END IF
        END DO
    END DO


  END SUBROUTINE relax_tteh_waves_smcl
  ! =========================================================================
  ! -op_rw_20181102



  
  
  ! +op_mn_20180620
  ! CURRENTLY UNDER DEVELOPMENT
  ! =========================================================================
  SUBROUTINE relax_tteh_mons_smcl(pres, lat, lon, delt, tstep, p_bot, p_top, lath, declat &
                                    , lonh, declon, offh, hperiod, amph, spinup, tteh)
 

 
    ! Calculate localized temperature tendency for generation of monsoon anticyclones.
    ! Vertical heating profile based on Lindgren et al. (2018); horizontal heating
    ! based on Schubert and Masarik (2006).
    
    ! Schubert, W. H., and M. T. Masarik, 2006: Potential vorticity aspects of the MJO.
    ! Dynam. Atmos. Oceans, 42, 127-151. Doi: 10.1016/j.dynatmoce.2006.02.003.
    ! https://doi.org/10.1016/j.dynatmoce.2006.02.003
    
    ! Lindgren, Sheshadri & Plumb (2018), Sudden stratospheric warming formation in an idealized general
    ! circulation model using three types of tropospheric forcing,
    ! Journal of Geophysical Research: Atmospheres, 123, https://doi.org/10.1029/2018JD028537


    USE messy_main_constants_mem, ONLY: pi



    IMPLICIT NONE
    INTRINSIC :: SIZE, abs, sin, cos, exp

    ! I/O
    REAL(DP),DIMENSION(:,:),  INTENT(OUT)   :: tteh       ! temperature tendency due to monsoonal heating

    REAL(DP),                 INTENT(IN)    :: spinup     ! length of spin up phase to gradually increase the forcing [day]
    REAL(DP),                 INTENT(IN)    :: amph       ! amplitude of heating modulated by frequency [K/day]
    REAL(DP),                 INTENT(IN)    :: hperiod    ! period of the heating signal [day]
    REAL(DP),                 INTENT(IN)    :: offh       ! constant heating [K/day]

    REAL(DP),                 INTENT(IN)    :: declon     ! decay of heating in longitude direction [deg]
    REAL(DP),                 INTENT(IN)    :: lonh       ! longitude center for the heating [deg]
    REAL(DP),                 INTENT(IN)    :: declat     ! decay of heating in latitude direction [deg]    
    REAL(DP),                 INTENT(IN)    :: lath       ! latitude center for the heating [deg]
    REAL(DP),                 INTENT(IN)    :: p_top      ! top pressure boundary         [Pa]
    REAL(DP),                 INTENT(IN)    :: p_bot      ! bottom pressure boundary      [Pa]

    REAL(DP),                 INTENT(IN)    :: tstep      ! timestep of the integration
    REAL(DP),                 INTENT(IN)    :: delt       ! length of time step
    REAL(DP), DIMENSION(:),   INTENT(IN)    :: lon        ! longitude
    REAL(DP), DIMENSION(:),   INTENT(IN)    :: lat        ! latitude
    REAL(DP), DIMENSION(:,:),   INTENT(IN)  :: pres       ! pressure

    ! Parameters
    REAL(DP), PARAMETER     ::  daytosec =  86400.0_dp  ! days to seconds

        ! LOCAL
    INTEGER :: jp, jk, lkproma, nlev
    REAL(DP):: supscale         ! spin up scaling
    REAL(DP):: perscale         ! periodic scaling
    REAL(DP):: londec           ! longitudinal decay factor
    REAL(DP):: latdec           ! latitudinal decay factor
    REAL(DP):: presdec          ! pressure decay factor
   
    REAL(DP):: amphKs           ! heating in K/s    
    REAL(DP):: offhKs           ! heating in K/s
    REAL(DP):: spinups          ! spin up time in s    
    REAL(DP):: hperiods         ! period of heating in s


    lkproma = SIZE(lat,1) !!!!!!!!!!!!!!!!!!!! CHECK if there is a lat with kproma numbering  !SIZE(temp,1)
    nlev    = SIZE(pres,2) !SIZE(temp,2)

    ! INIT
    tteh(:,:) = 0.0_dp

    ! heating rates (K/day) to (K/s)
    amphKs = amph/daytosec
    offhKs = offh/daytosec
    
    ! time periods (days) to (sec)
    spinups = spinup*daytosec
    hperiods = hperiod*daytosec


    ! spin up scaling (linear)
    IF( (delt*tstep) .LE. spinups ) THEN
           supscale = 1.0_dp - 1.0_dp*((spinups-(delt*tstep))/spinups)
      ELSE
           supscale = 1.0_dp
    END IF

    ! periodic scaling
    IF (hperiods .GT. 0.0) THEN
       perscale = offhKs + amphKs*sin(2.0_dp*pi*(delt*tstep)/hperiods)
       ELSE
       perscale = offhKs
    END IF 

    ! loop through horizontal dims
    DO jp = 1, lkproma


       ! longitudinal decay factor as in Schubert and Masarik (2006)
       
       ! deal with lonh + declon > 360 or lonh-declon < 0

       IF ( abs(lon(jp)-lonh) .LE. declon) THEN

               londec = 0.5_dp*(1.0_dp+cos(pi*(lon(jp)-lonh)/declon))

       ELSE IF ( lonh + declon .GE. 360.0_dp .AND. lon(jp)+360.0_dp .LE. lonh + declon) THEN

               londec = 0.5_dp*(1.0_dp+cos(pi*(lon(jp)+360.0_dp-lonh)/declon))
                
       ELSE IF (lonh - declon .LT. 0.0_dp .AND. lon(jp)-360.0_dp .GE. lonh - declon) THEN
                
               londec = 0.5_dp*(1.0_dp+cos(pi*(lon(jp)-360.0_dp-lonh)/declon))

       ELSE
               londec = 0.0_dp
       END IF
       
       

      ! latitudinal decay factor as in Schubert and Masarik (2006)
      latdec = exp(-((lat(jp)- lath)/declat)**2)
       

       
      ! loop through vertical
      DO jk = 1, nlev

         ! vertical/pressure decay factor as in Lindgren et al. (2018)
         IF ( (pres(jp,jk) .LE. p_bot) .AND. (pres(jp,jk) .GE. p_top) ) THEN

                   presdec = sin( pi * log(pres(jp,jk)/p_bot) / log(p_top/p_bot) )
            ELSE
                   presdec = 0.0_dp
         END IF
 

        ! Final calculation of the temperature tendency
        tteh(jp,jk) = supscale*perscale*presdec*londec*latdec


      END DO
    END DO

  
  END SUBROUTINE relax_tteh_mons_smcl
  ! =========================================================================
! -op_mn_20180620 
  
  
  
  ! =========================================================================
  ! ### add your own private subroutines here
  ! =========================================================================

! **********************************************************************
END MODULE messy_relax
! **********************************************************************

