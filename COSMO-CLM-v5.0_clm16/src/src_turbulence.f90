!+ Source module  "turbulence"
!------------------------------------------------------------------------------

MODULE src_turbulence

!------------------------------------------------------------------------------
!
! Description:
!
!   The module "turbulence" calculates the diffusion coefficients for turbulent
!   vertial transport of momentum (TKVM) and heat (TKVH). At the surface, the
!   dimensionless transfer coefficients for momentum and heat are calculated.
!
!   The module contains the subroutines partura and prankolmo_rk (for the
!   diffusion coefficients) and parturs (for the transfer coefficients).
!   These routines are called from the basic driving routine of the model.
!   The subroutine prankolmo_rk specifies the diffusion coefficients following
!   a TKE-depending approach by Prandtl and Kolmogorov.
!
! Current Code Owner: DWD, Matthias Raschendorfer
!  phone:  +49  69  8062 2708
!  fax:    +49  69  8062 3721
!  email:  Matthias.Raschendorfer@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        1998/03/11 Guenter Doms     
!  Initial release
! 1.4        1998/05/22 Guenther Doms
!  Some modifications to use the routines of this module in the
!  two time level integration scheme.
!  A new subroutine "vertical diffusion" is included to calculate
!  the diffusional time tendencies explicitly. These tendencies
!  are added to the total tendency fields.
! 1.9        1998/09/16 Guenther Doms
!  The diffusion coefficients are now limited according to the 
!  stability criterion for the explicit numerical solution.
! 1.20       1999/01/07 Guether Doms
!  Renaming of some global variables
! 1.30       1999/06/24 Matthias Raschendofer
!  Introduction of the field 'tke' in sub. 'partura' and 'parturs'
!  as an output variable.
! 1.39       2000/05/03 Ulrich Schaettler
!  Eliminated lphys: if this routine is called, lphys is .TRUE. in any case.
! 2.2        2000/08/18 Matthias Raschendorfer
!  Correction of an error in the tke diagnosis in routine 'parturs'.
! 2.8        2001/07/06 Ulrich Schaettler
!  Now use actual surface fluxes instead of summarized
! 2.18       2002/07/16 Guenther Doms
!  Calculation of roughness lenght for scalar quantities over sea
!  follows surface parameterization similar to ECMWF formulation (as in GME).
!  The Charnock constant for z0 is increased from 0.0123 to 0.0150.
!  Changes by Almut Gassmann
!  New routine vertical_diffusion_impl for 2 time level scheme;
!  Limited fluxes in explicit version
! 2.19       2002/10/24 Ulrich Schaettler
!  Eliminated subroutine vertical_diffusion_impl for 2 time level scheme
!  (is now solved in a different way)
! 3.7        2004/02/18 Ulrich Schaettler
!  Eliminated parameter itype_gscp.
! 3.13       2004/12/03 Jochen Foerstner
!  Introduction of 3D turbulence
! 3.14       2005/01/25 Jochen Foerstner
!  Modifications for the 3D turbulence scheme
! 3.18       2006/03/03 Jochen Foerstner
!  Modifications for the 3D turbulence scheme
! 3.21       2006/12/04 Dmitrii Mironov
!  Changes to use the FLake model
! V3_23        2007/03/30 Matthias Raschendorfer
!  'akt' is now used form MODULE data_turbulence
! V4_3         2008/02/25 Matthias Raschendorfer
!  Calculation of the 3D diagnostic field 'edr' for the eddy dissipotion rate.
! V4_8         2009/02/16 Ulrich Blahak
!  Bug fix in SR prankolmo_rk: the resolution in meter was given with a fixed
!  value here: now it is computed depending on dlat, dlon
! V4_10        2009/09/11 Jan-Peter Schulz
!  Modifications for new sea-ice model
! V4_12        2010/05/11 Ulrich Schaettler
!  Renamed t0 to t0_melt because of conflicting names
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_17        2011/02/24 Ulrich Blahak
!  itype_turb=7/8: Changed constant zcs in the LES turbulence scheme
!    from 0.25 to 0.15 and set horiz. diff. coeff. equal to 3.0 times 
!    the vertical diff. coeffs.
!  Eliminated my_peri_neigh + other adaptions for periodic BCs
! V4_18        2011/05/26 Ulrich Schaettler
!  Changed the code owner
! V4_20        2011/08/31 Matthias Raschendorfer
!  Eliminated unnecessary MAX-construct for index computations
! V4_23        2012/05/10 Ulrich Schaettler
!  Use field sqrtg_r_w from new module grid_metrics_utilities
!  Some editorial changes
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer
!  Replaced qx-variables by using them from the tracer module
! V4_27        2013/03/19 Astrid Kerkweg, Ulrich Schaettler
!  Introduced MESSy interface
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler, Ulrich Blahak
!  Unification of MESSy interfaces and COSMO Tracer structure
!  Bugfixes in subroutine prankolmo_rk:
!  - zlen has to be a 3D field instead of a scalar
!    because it is carried from one loop to a subsequent loop.
!  - zsgesh and zsgesv are already the horizontal and vertical
!    components of the squared deformation. Therefore remove
!    the square operator from calculations of ztpm and ztpm_v.
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:

USE data_parameters , ONLY :   &
    ireals,     & ! KIND-type parameters for real variables
    iintegers     ! KIND-type parameter for standard integer variables

!------------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &

! 2. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------

    ie,           & ! number of grid points in zonal direction
    je,           & ! number of grid points in meridional direction
    ke,           & ! number of grid points in vertical direction
    ke1,          & ! ke + 1

! 3. start- and end-indices for the computations in the horizontal layers
! -----------------------------------------------------------------------
!    These variables give the start- and the end-indices of the 
!    forecast for the prognostic variables in a horizontal layer.
!    Note, that the indices for the wind-speeds u and v differ from 
!    the other ones because of the use of the staggered Arakawa-B-grid.
!    
    ie_tot,       & ! number of grid points in zonal direction
    je_tot,       & ! number of grid points in meridional direction
    istart,       & ! start index for the forecast of w, t, qv, qc and pp
    iend,         & ! end index for the forecast of w, t, qv, qc and pp
    istartu,      & ! start index for the forecast of u
    iendu,        & ! end index for the forecast of u
    istartv,      & ! start index for the forecast of v
    iendv,        & ! end index for the forecast of v
    istartpar,    & ! start index for computations in the parallel program
    iendpar,      & ! end index for computations in the parallel program
    jstart,       & ! start index for the forecast of w, t, qv, qc and pp
    jend,         & ! end index for the forecast of w, t, qv, qc and pp
    jstartu,      & ! start index for the forecast of u
    jendu,        & ! start index for the forecast of u
    jstartv,      & ! start index for the forecast of v
    jendv,        & ! end index for the forecast of v
    jstartpar,    & ! start index for computations in the parallel program
    jendpar,      & ! end index for computations in the parallel program

! 4. constants for the horizontal rotated grid and related variables
! ------------------------------------------------------------------
    dlon,         & ! grid point distance in zonal direction (in degrees)
    dlat,         & ! grid point distance in meridional direction (in degrees)
    startlat_tot, & ! transformed latitude of the lower left grid point
                    ! of the total domain (in degrees, N>0)
    eddlon,       & ! 1 / dlon
    edadlat,      & ! 1 / (radius of the earth * dlat)
    degrad,       & ! factor for transforming degree to rad

! 5. variables for the time discretization and related variables
! --------------------------------------------------------------
    dt,           & !

! 8. Organizational variables to handle the COSMO humidity tracers
! ----------------------------------------------------------------
    idt_qv,  idt_qc

! end of data_modelconfig
!------------------------------------------------------------------------------

USE data_constants  , ONLY :   &

! 2. physical constants and related variables
! -------------------------------------------
    t0_melt,      & ! 273.15 K
    r_d,          & ! gas constant for dry air
    r_v,          & ! gas constant for water vapor
    rvd_m_o,      & ! r_v/r_d - 1
    rdocp,        & ! r_d / cp_d
    lh_v,         & ! latent heat of vapourization
    lh_s,         & ! latent heat of sublimation
    cp_d,         & ! specific heat of dry air at constant pressure
    cpdr ,        & ! 1.0/cp_d
    g,            & ! acceleration due to gravity
    gq,           & ! g*g
    gr,           & ! 1./g
    r_earth         ! mean radius of the earth

! end of data_constants

!------------------------------------------------------------------------------

USE data_turbulence     , ONLY :   &

! 1. constants for parametrizations
! ---------------------------------
    akt,           &! von Karman-constant
    d_mom           ! = 16.6 (dissipation constant)

!------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------
    rho0,           & ! reference density at the full model levels    (kg/m3)
    hhl ,           & ! height of model half levels                   ( m )
    p0  ,           & ! reference pressure at main levels             (pa )
    dp0 ,           & ! pressure thickness of layer                   (pa )

! 2. external parameter fields                                        (unit)
! ----------------------------
    gz0,            & ! roughness length * g
    fr_land,        & ! fraction of land in a grid element
    depth_lk ,      & ! lake depth                                    ( m  )
    rmy,            & ! Davis-parameter for boundary relaxation         --

! 3. prognostic variables                                             (unit)
! -----------------------
    u          ,    & ! zonal wind speed                              ( m/s )
    v          ,    & ! meridional wind speed                         ( m/s )
    w          ,    & ! vertical velocity                             ( m/s )
    t          ,    & ! temperature                                   (  k  )
    pp         ,    & ! deviation from the reference pressure         ( pa  )
    tke        ,    & ! SQRT(2 * turbulent kinetik energy)            ( m/s )
                      ! (defined on half levels)

! 4. tendency fields for the prognostic variables                     (unit )
! -----------------------------------------------
!    timely deviation  by diabatic and adiabatic processes
!    without sound-wave terms
    tketens,        & ! tke-tendency                                  (m2/s3)

! 5. fields for surface values and soil model variables               (unit )
! -----------------------------------------------------
    ps        ,     & ! surface pressure                              ( pa  )
    h_ice     ,     & ! ice thickness                                 (  m  )
    t_s       ,     & ! temperature of the ground surface             (  k  )
    t_g       ,     & ! weighted surface temperature                  (  k  )
    qv_s      ,     & ! specific water vapor content on the surface   (kg/kg)

! 6. fields that are computed in the parametrization and dynamics     (unit )
! ---------------------------------------------------------------
    rho        ,    & ! density of moist air (at timelevel nnow)      (kg/m^3)
    a1t, a2t   ,    & ! implicit weight of vertical diffusion
    tkvm, tkhm ,    & ! turbulent diffusion coefficients for momentum (m^2/s)
    tkvh, tkhh ,    & ! turbulent diffusion coefficients for heat     (m^2/s)
                      ! and moisture
    edr        ,    & ! eddy dissipation rate of TKE                  (m^2/s^3)
    tcm        ,    & ! turbulent diffusion coefficients for momentum   --
    tch        ,    & ! turbulent diffusion coefficients for heat       --
    qrs        ,    & ! precipitation water (water loading)           (kg/kg)

! 7. external parameter fields                                        (unit )
! ----------------------------
    crlat      ,    & ! cosine of transformed latitude
    acrlat            ! 1 / ( crlat * radius of the earth )           ( 1/m )

! end of data_fields

!------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 1. start and end of the forecast
! --------------------------------
    ntstep,       & ! actual time step
                    ! indices for permutation of three time levels
    nold,         & ! corresponds to ntstep - 1
    nnow,         & ! corresponds to ntstep
    nnew,         & ! corresponds to ntstep + 1
    ntke,         & ! time index of SQRT(2*TKE)

! 3. controlling the physics
! --------------------------
    lphys,        & ! forecast with physical parametrizations
    lconv,        & ! forecast with convection
    lgsp,         & ! forecast with grid-scale precipitation
    lseaice,      & ! forecast with sea ice model
    llake,        & ! forecst with lake model FLake
    itype_gscp,   & ! type of grid-scale precipitation physics
    itype_turb,   & ! type of turbulent diffusion parametrization
    l3dturb,      & ! 3D-turbulence (additional horizontal diffusion)
    lprog_tke,    & ! prognostic treatment of TKE (for itype_turb=5/7)

! 5. additional control variables
! -------------------------------
    l2tls           ! forecast with 2-TL integration scheme

! end of data_runcontrol 
!------------------------------------------------------------------------------

USE data_flake, ONLY : &
    ! flake_parameters
    h_Ice_min_flk      ! Minimum ice thickness [m]

!------------------------------------------------------------------------------

USE data_parallel,      ONLY :  &
    my_cart_neigh, & ! neighbors of this subdomain in the cartesian grid
    my_cart_id

!------------------------------------------------------------------------------

USE grid_metrics_utilities, ONLY: &
    sqrtg_r_w         ! reciprocal square root of G at w points       ( 1/m )

!------------------------------------------------------------------------------

USE environment,        ONLY :  &
    model_abort

!------------------------------------------------------------------------------

USE src_tracer,         ONLY: trcr_get, trcr_errorstr

!==============================================================================

IMPLICIT NONE

LOGICAL                  ::  &
    l_ls_ice               ! Logical switch to discriminate between ice-covered
                           ! and open-water grid boxes.
                           ! Notice that this is done differently for lake
                           ! and sea grid boxes if the lake model is used.

!==============================================================================

CONTAINS

!==============================================================================
!+ Module procedure partura in "src_turbulence" for computing the coefficients
!+ for vertical diffusion
!------------------------------------------------------------------------------

SUBROUTINE partura             

!------------------------------------------------------------------------------
!
! Description:
!   This module procedure calculates the coefficients TKVM and TKVH for 
!   vertical turbulent diffusion of momentum and heat, respectively.
!   TKVM and TKVH are defined on half levels.                              
!
! Method:
!   The closure is second order on level 2 in the Mellor/Yamada notation.
!   The impact of subgrid condensation/evaporation in not taken into account
!   in this version.
!
! Subroutine history:
! S-Version Date       Comment
! -------   ----       -------
! 1.0       26.06.96   Original code based on EM/DM.    Guenter Doms
! 1.1       21.04.98   Selection of time levels depending on the
!                      time integration scheme used. G.D.
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!
!------------------------------------------------------------------------------
!
! Declarations:
!
!------------------------------------------------------------------------------

! Modules used: These are declared in the module declaration section
! -------------

! Subroutine arguments: None
! --------------------
!
! Local parameters:
! ----------------
  REAL    (KIND=ireals   ), PARAMETER ::  &
    ! basic constants of the parameterization scheme
    zric   = 0.380_ireals,    & !
    zalhn  = 1.00_ireals ,    & !   
    zdze   = 500.0_ireals,    & !     
    zdvbmin= 0.01_ireals ,    & !
    ! constants for the numerical scheme 
    zalf   = 0.900_ireals,    & !
    zemalf = 0.100_ireals,    & !  
    ztmmin = 0.010_ireals,    & ! min value as part of neutral value for TKVM
    zthmin = 0.007_ireals,    & ! min value as part of neutral value for TKVH
    ztkmmin= 1.000_ireals,    & ! mimimum absolut value TKVM
    ztkhmin= 0.400_ireals       ! minimum alsolut value TKVH

! Local scalars:
! -------------
  INTEGER (KIND=iintegers) ::  &
    k     ,            & ! loop index in vertical direction              
    i     ,            & ! loop index in x-direction              
    j     ,            & ! loop index in y-direction              
    jzloc, izloc,      & ! local start indices for some computations
    izerror,           & !
    ms    ,            & !
    nx    ,            & ! time-level for computation
    izg    ,           & !
    jzg    ,           & !
    mzg    ,           & !
    nzp    ,           & !
    nztgpk ,           & !
    nzpa

  REAL    (KIND=ireals   ) ::  &
    zumo  ,            & ! 
    zvmo  ,            & ! 
    zumu  ,            & ! 
    zvmu  ,            & ! 
    zdpo  ,            & !
    zdpu  ,            & !  
    zpio  ,  z1o, z2o, & !
    zpiu  ,  z1u, z2u, & !
    zdpn2 ,            & !
    zthvo ,            & !
    zthvu ,            & !
    zthn  ,            & !
    zdza  ,            & !
    zldi  ,            & !
    zmzb  ,            & !
    zdzq  ,            & !
    zrf   ,            & !
    zgam  ,            & !
    zgs   ,            & !
    zsm   ,            & !
    zalh  ,            & !
    ztkvhom,           & !
    ztkvmom              !

  CHARACTER (LEN=80)         :: yzerrmsg
  CHARACTER (LEN=25)         :: yzroutine

! Local (automatic) arrays:
! -------------------------
  REAL    (KIND=ireals   ) ::  &
    zvbq    (ie,je,2:ke),& !  
    ztpm    (ie,je,2:ke),& ! 
    ztph    (ie,je,2:ke),& !
    zrinum  (ie,je,2:ke),& !  Richardson number                            
    zaa     (2)         ,& !
    zab1    (2)         ,& !
    zab2    (2)         ,& !
    zac1    (2)         ,& !
    zac2    (2)         ,& !
    zac3    (2)            !

! Tracer pointers
!----------------
  REAL (KIND=ireals), POINTER :: &
    qv  (:,:,:)   => NULL()   ! QV at nx

!------------ End of header ---------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine partura             
!------------------------------------------------------------------------------
 
yzroutine = 'partura'

! Initialize parameterization constants for stable and unstable cases

!               (stable) (unstable)
DATA zaa      / 3.7_ireals     , 4.025_ireals   /
DATA zab1     / 2.5648_ireals  , 3.337_ireals   /
DATA zab2     / 1.1388_ireals  , 0.688_ireals   /
DATA zac1     / 0.8333_ireals  , 1.285_ireals   /
DATA zac2     / 0.2805_ireals  , 0.2305_ireals  /
DATA zac3     / 0.1122_ireals  ,-0.1023_ireals  /

! Select time-level for TKVM, TKVH computations

IF (l2tls) THEN
  nx = nnow
ELSE
  nx = nold
ENDIF

ntke = 1
 
! retrieve the required microphysics tracers (at timelevel nx)
CALL trcr_get(izerror, idt_qv, ptr_tlev = nx, ptr = qv)
IF (izerror /= 0) THEN
  yzerrmsg = trcr_errorstr(izerror)
  CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
ENDIF

! Start of calculation for each layer from top to bottom  
! ------------------------------------------------------
!FPP$ SELECT CONCUR
 
loop_over_levels: DO k = 2, ke
 
  !----------------------------------------------------------------------------
  ! Section 1: Averaging the horizontal velocities (C-grid)
  !----------------------------------------------------------------------------

  ! In the interior of the subdomain and the northern and eastern boundary of 
  ! the total domain

  ! What do we need for periodic boundary conditions ???

  IF (my_cart_neigh(1) == -1) THEN
    izloc = 2
  ELSE
    izloc = istartpar
  ENDIF
  IF (my_cart_neigh(4) == -1) THEN
    jzloc = 2
  ELSE
    jzloc = jstartpar
  ENDIF

  DO   j   = jzloc , jendpar
    DO i   = izloc , iendpar
      zumo  = 0.5*(  u(i,j  ,k-1,nx) + u(i-1,j  ,k-1,nx)  )
      zvmo  = 0.5*(  v(i,j  ,k-1,nx) + v(i  ,j-1,k-1,nx)  )
      zumu  = 0.5*(  u(i,j  ,k  ,nx) + u(i-1,j  ,k  ,nx)  )
      zvmu  = 0.5*(  v(i,j  ,k  ,nx) + v(i  ,j-1,k  ,nx)  )
      zvbq(i,j,k) = (zumu - zumo)**2 + (zvmu - zvmo)**2  
    ENDDO
  ENDDO                   
 
  ! Special treatment at the southern boundary of the total domain
  IF (my_cart_neigh(4) == -1) THEN
    DO i = izloc , iendpar
      zvbq(i,1,k) = zvbq(i,2,k)
    ENDDO       
  ENDIF
 
  ! Special treatment at the western boundary of the total domain
  IF (my_cart_neigh(1) == -1) THEN
    DO j = 1 , jendpar
      zvbq(1,j,k) = zvbq(2,j,k)
    ENDDO       
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 2: Calculation of the forcing functions           
  !----------------------------------------------------------------------------

  DO   j = jstartpar , jendpar
    DO i = istartpar , iendpar
       zdpu      = dp0(i,j,k  )
       zdpo      = dp0(i,j,k-1)
       zdpn2     = zdpu + zdpo
       zpiu      =(1.0E-5_ireals * (p0(i,j,k  ) + pp(i,j,k  ,nx)))**rdocp
       zpio      =(1.0E-5_ireals * (p0(i,j,k-1) + pp(i,j,k-1,nx)))**rdocp
       zthvu     = t(i,j,k  ,nx)*(1.0+rvd_m_o*qv(i,j,k  ))/zpiu
       zthvo     = t(i,j,k-1,nx)*(1.0+rvd_m_o*qv(i,j,k-1))/zpio
       zthn      = (zdpu*zthvo + zdpo*zthvu)/zdpn2
       zdza      = 0.5*( hhl(i,j,k+1) - hhl(i,j,k-1) )
       ztpm  (i,j,k) = zvbq(i,j,k) / (zdza**2)
       ztpm  (i,j,k) = MAX( ztpm(i,j,k) , (zdvbmin/zdza)**2 )
       ztph  (i,j,k) = g*  ((zthvu - zthvo)/zdza )/zthn
       zrinum(i,j,k) = ztph(i,j,k)/ztpm(i,j,k)
    ENDDO
  ENDDO        

  !----------------------------------------------------------------------------
  ! Section 3: Calculation of the diffusion coefficients for time level nx
  !----------------------------------------------------------------------------

  DO   j = jstartpar , jendpar
    DO i = istartpar , iendpar
 
      zmzb = hhl(i,j,k) - hhl(i,j,ke1)
      zdzq = (akt*zmzb/(1.+zmzb/zdze))**2

          ! very stable case (ri>ric)
      IF (zrinum(i,j,k) >= zric) THEN
          tkvm(i,j,k) = ztmmin*zdzq*SQRT( ztpm(i,j,k) )
          tkvh(i,j,k) = zthmin*zdzq*SQRT( ztpm(i,j,k) )

          tke(i,j,k,ntke)=0.0
          edr(i,j,k)=0.0
      ELSE
 
          ! stable and unstable cases 
          ms = 1
          IF(zrinum(i,j,k) <= 0.0) THEN
              ms = 2
          ENDIF
 
          ! Flux-Richardson number as function
          ! of the Gradient-Richardson number
          zrf  = zac1(ms) * ( zrinum(i,j,k) + zac2(ms) -  &
                  SQRT(zrinum(i,j,k)**2 - zac3(ms)*zrinum(i,j,k) +   &
                  zac2(ms)**2) )
          zrf  = MIN ( zrf, 0.99999_ireals )
 
          ! Gam(rf)
          zgam = zrf/(1.0-zrf)
 
          ! SM(Gam(rf))**1.5 and alh(Gam(rf))
          zgs  = (1.0-zab1(ms)*zgam)/(1.0-zab2(ms)*zgam)
          zsm  = SQRT(((1.-zaa(ms)*zgam)/zgs)**3)
          zalh = zalhn*zgs
 
          ! the diffusion coefficients
          ztkvmom     = zdzq*zsm*SQRT(ztpm(i,j,k) - zalh*ztph(i,j,k))
          ztkvhom     = ztkvmom*zalh
          tkvm(i,j,k) = zalf*ztkvmom + zemalf*tkvm(i,j,k) 
          tkvh(i,j,k) = zalf*ztkvhom + zemalf*tkvh(i,j,k) 

          zldi = d_mom*SQRT(zdzq)
          edr(i,j,k) = ( zldi*(ztkvmom*ztpm(i,j,k) &
                       - ztkvhom*ztph(i,j,k)) )
          tke(i,j,k,ntke) = edr(i,j,k)**(1./3.)
          edr(i,j,k) = edr(i,j,k)/zldi
      ENDIF

      ! Set lower limit for diffusion coefficients
      tkvm (i,j,k) = MAX ( ztkmmin, tkvm(i,j,k) )
      tkvh (i,j,k) = MAX ( ztkhmin, tkvh(i,j,k) )

    ENDDO       
  ENDDO       

ENDDO loop_over_levels

!------------------------------------------------------------------------------
! End of module procedure partura
!------------------------------------------------------------------------------

END SUBROUTINE partura

!------------------------------------------------------------------------------
!+ Module procedure parturs in "src_turbulence" for computing the transfer
!+ coefficient
!------------------------------------------------------------------------------

SUBROUTINE  parturs

!------------------------------------------------------------------------------
!
! Description:
!   This module procedure calculates the transfer coefficients for momentum
!   (TCM) and heat and moisture (TCH) at the surface as well as the roughness
!   length over sea (GZ0).
!
! Method:
!   Dyer-Businger equations as modified by Louis.
!
!------------------------------------------------------------------------------

! Subroutine arguments: None
! --------------------
!
! Local parameters:
! ----------------
  REAL    (KIND=ireals   ), PARAMETER ::  &
    ! basic constants of the parameterization scheme and meteorological and
    ! numerical parameters
    zah     = 5.300_ireals,    & !
    zgz0hh  = 0.980_ireals,    & ! upper limit for roughness length for heat
    zalphaf = 1.000_ireals,    & !
    zalpha0 = 0.0150_ireals,   & ! Charnock constant for momentum
    zalphah = 0.60_ireals  ,   & !
    zviscos = 1.5E-5       ,   & ! kinematic viscosity constant (m**2/s)
    zbeta10 = 0.0420_ireals,   & !
    z10     = 10.0_ireals  ,   & !
    zvmin   = 0.01_ireals  ,   & ! Minimum wind velocity
    ztmmin  = 0.140_ireals ,   & ! Minimum value for TCM
    zthmin  = 0.010_ireals ,   & ! Minimum value for TCH
    zed3    = 1.0_ireals/3.0_ireals,  & !
    zd3     = 2.0_ireals/3.0_ireals,  & !
    zctke   = 0.516_ireals

! Local scalars:
! -------------
  INTEGER (KIND=iintegers) ::  &
    i     ,            & ! loop index in x-direction              
    j     ,            & ! loop index in y-direction              
    izloc , jzloc,     & ! start-indices for local computations
    izerror,           & !
    nx    ,            & ! time-level for computation
    izg   ,            & !
    jzg   ,            & !
    mzg   ,            & !
    nzpa                 !

  REAL    (KIND=ireals   ) ::  &
    zum   ,            & ! 
    zvm   ,            & ! 
    ztvg  ,            & ! 
    ztvke ,            & ! 
    zxi   ,            & ! 
    zxih  ,            & ! 
    zy    ,            & !
    zgz0d ,            & !
    zgz0dd,            & !
    zustar,            & !
    zdz                  !

  CHARACTER (LEN=80)            :: yzerrmsg
  CHARACTER (LEN=25)            :: yzroutine

! Local (automatic) arrays:
! -------------------------
  REAL    (KIND=ireals   ) ::  &
    zvpb    (ie,je),   & !  
    zx      (ie,je),   & ! 
    ztcm    (ie,je),   & ! 
    ztch    (ie,je),   & ! 
    zdfip   (ie,je),   & ! 
    zris    (ie,je),   & !
    zgz0m   (ie,je),   & ! Roughness lenght for momentum
    zgz0h   (ie,je)      ! Roughness lenght for scalars (heat, moisture)

! Tracer pointers
!----------------
  REAL (KIND=ireals), POINTER :: &
    qv  (:,:,:) => NULL()         ! QV at nx

!------------ End of header ---------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine parturs             
!------------------------------------------------------------------------------

  yzroutine = 'parturs'

  ! Select time-level for TCM, TCH computations
  IF (l2tls) THEN
    nx = nnow
  ELSE
    nx = nold
  ENDIF

  ! retrieve the required microphysics tracers (at timelevel nx)
  CALL trcr_get(izerror, idt_qv, ptr_tlev = nx, ptr = qv)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

!------------------------------------------------------------------------------
! Section 1: Averaging the horizontal velocities (C-grid) at lowest layer
!------------------------------------------------------------------------------

  ! In the interior of the subdomain and the northern and eastern boundary of 
  ! the total domain

  IF (my_cart_neigh(1) == -1) THEN
    izloc = 2
  ELSE
    izloc = istartpar
  ENDIF
  IF (my_cart_neigh(4) == -1) THEN
    jzloc = 2
  ELSE
    jzloc = jstartpar
  ENDIF

  DO   j   = jzloc , jendpar
    DO i   = izloc , iendpar
      zum  = 0.5*(  u(i,j  ,ke,nx) + u(i-1,j  ,ke,nx) )
      zvm  = 0.5*(  v(i,j  ,ke,nx) + v(i  ,j-1,ke,nx) )
      zvpb(i,j) =  MAX( SQRT(zum**2 + zvm**2), zvmin )
    ENDDO
  ENDDO                   
 
  ! Special treatment at the southern boundary of the total domain
  IF (my_cart_neigh(4) == -1) THEN
    DO i = izloc , iendpar
      zvpb(i,1) = zvpb(i,2)
    ENDDO       
  ENDIF
 
  ! Special treatment at the western boundary of the total domain
  IF (my_cart_neigh(1) == -1) THEN
    DO j = 1 , jendpar
      zvpb(1,j) = zvpb(2,j)
    ENDDO       
  ENDIF

!------------------------------------------------------------------------------
! Section 2: Calculation of the forcing functions and transfer coefficients    
!------------------------------------------------------------------------------

  DO j = jstartpar , jendpar
    DO i = istartpar , iendpar
 
      ztvg       = t_g(i,j,   nx)*(1.0 + rvd_m_o*qv_s(i,j,   nx))
      ztvke      = t  (i,j,ke,nx)*(1.0 + rvd_m_o*qv  (i,j,ke))
      zdfip(i,j) = 0.5*g*( hhl(i,j,ke) - hhl(i,j,ke1) ) 
      zx   (i,j) = ( ztvke - ztvg + cpdr*zdfip(i,j) )*zdfip(i,j)/t_g(i,j,nx)
      zris (i,j) = zx(i,j)/zvpb(i,j)**2
      zx   (i,j) = ABS(zx(i,j))

      ! Set the logical switch l_ls_ice to distinguish between ice covered
      ! and open water sea or lake grid points.

      IF (fr_land(i,j) < 0.5) THEN
        ! Water point.
        IF (.NOT. lseaice) THEN
          ! Sea ice model is not used.
          ! Ice surface if SST is less than the salt water freezing temperature.
          l_ls_ice = t_g(i,j,nx) < t0_melt - 1.7_ireals
        ELSE
          ! Sea ice model is used.
          ! Ice surface if ice is present.
          l_ls_ice = h_ice(i,j,nnow) > 0.0_ireals
        ENDIF
        IF (llake) THEN
          ! Lake model is used.
          ! Ice surface if this is a lake point AND ice is present.
          IF ((depth_lk(i,j) > 0.0_ireals) .AND. (h_ice(i,j,nnow) >= h_Ice_min_flk)) &
          l_ls_ice = .TRUE.
        ENDIF
      ENDIF

      ! Distinguish land and sea points for z0 calculations
      ! ---------------------------------------------------
      IF ( fr_land(i,j) < 0.5 .AND. gz0(i,j) <= 0.0 ) THEN

        ! Use ice surface roughness or open-water surface roughness according
        ! to l_ls_ice

        ! set z0 to 0.001m over sea ice      
        IF ( l_ls_ice ) THEN
          gz0(i,j) = 0.001*g
  
        ! set z0 over water for first iteration
        ELSE
          zgz0d = zalpha0*(zvpb(i,j)/(1.0/zbeta10       &
                  + LOG( zdfip(i,j)/(g*z10) )/akt))**2
          IF ( zris(i,j) < 0.0 ) THEN
             zgz0dd = zalphaf**3 / ( zah*SQRT(zdfip(i,j)) ) &
                      * ( zalpha0*zx(i,j) )**1.5
          ELSE
             zgz0dd = 0.0
          ENDIF
        ! gz0(i,j) = MAX(zgz0d , zgz0dd)
          gz0(i,j) = zgz0d + zgz0dd
        ENDIF
      ENDIF

      ! Limit the roughness length for momentum to a value smaller than
      ! half of the lower layer thickness

      zgz0m(i,j) = MIN( gz0(i,j), zdfip(i,j) - g )

      ! Derive roughness length for heat over open sea using the friction
      ! velocity ustar derived from z0 for momentum (Charnock formula)
      ! Over land and sea ice, z0h = z0m.

!_cdm replaced by calculations including the FLake model
!     IF (  fr_land(i,j) < 0.5 .AND.  t_g(i,j,nx) >=  T0-1.7 ) THEN
!       zustar = SQRT ( zgz0m(i,j)/zalpha0 )
!       zgz0h(i,j) = g*zviscos*zalphah / MAX( 1.0E-8_ireals, zustar )
!     ELSE
!       zgz0h(i,j) = zgz0m(i,j)
!     ENDIF

      IF ( fr_land(i,j) < 0.5 ) THEN
        ! Open-water or ice surface

        ! Use ice surface roughness or open-water surface roughness according
        ! to l_ls_ice
        IF ( l_ls_ice ) THEN
          ! Ice surface
          zgz0h(i,j) = zgz0m(i,j)
        ELSE
          ! Open-water surface
          zustar = SQRT ( zgz0m(i,j)/zalpha0 )
          zgz0h(i,j) = g*zviscos*zalphah / MAX( 1.0E-8_ireals, zustar )
        ENDIF

      ELSE
        ! Land surface
        zgz0h(i,j) = zgz0m(i,j)
      ENDIF

      ! Limit the roughness length for momentum to a fixed value

      zgz0h(i,j)  = MIN ( zgz0h(i,j), zgz0hh)

      zxi         = zdfip(i,j)/ zgz0m(i,j)
      zxih        = zdfip(i,j)/ zgz0h(i,j)
      zy          = (akt/LOG(zxi))**2
   
      ! Stable case for land and water 
      ! ------------------------------
      IF ( zris(i,j) >= 0.0 ) THEN
  
        ztcm(i,j) = zy*zvpb(i,j)*MAX ( ztmmin, 1.0/  &
                   (1.0 + 10.0*zris(i,j)/SQRT(1.0 + 5.0*zris(i,j)) ) ) 
        ztch(i,j) = akt**2/(LOG(zxi)*LOG(zxih))*zvpb(i,j)*    &
                     MAX ( zthmin, 1.0/(1.0 + 15.0*zris(i,j)*  &
                     SQRT( 1.0 + 5.0*zris(i,j) ) ) ) 
   
        ! new z0 value over sea
        IF ( fr_land(i,j) < 0.5 ) THEN

          ! Use ice surface roughness or open-water surface roughness 
          ! according to l_ls_ice
          IF ( l_ls_ice ) THEN
            ! set z0 to 0.001m over sea ice
            gz0(i,j) = 0.001*g
          ELSE
            ! calculate z0 over water by Charnock formula
            gz0 (i,j) = zalpha0*ztcm(i,j)*zvpb(i,j)
          ENDIF
        ENDIF
   
      ! Unstable case
      ELSE
   
        ! Land
        IF ( fr_land(i,j) >= 0.5 ) THEN
          ztcm(i,j) = zy*zvpb(i,j)*(1.0 - 10.0*zris(i,j)/   &
                      (1.0 + 75.0*zy*(zxi**zed3-1.0)**1.5*SQRT(-zris(i,j)) ))
          ztch(i,j) = akt**2/(LOG(zxi)*LOG(zxih))*zvpb(i,j)*   &
                      (1.0-15.0*zris(i,j)/(1.0 + 75.0*SQRT(zy)*akt/LOG(zxih)* &
                      (zxih**zed3-1.0)**1.5*SQRT(-zris(i,j)) ))
   
        ! Water
        ELSE
          ztcm(i,j) = zy*zvpb(i,j)*(1.0 - 10.0*zris(i,j)/   &
                      (1.0 + 75.0*zy*SQRT(-zris(i,j)*zxi) ))
          ztch(i,j) = akt**2/(LOG(zxi)*LOG(zxih))*zvpb(i,j)*  &
                      (1.0-15.0*zris(i,j)/(1.0 + 75.0*SQRT(zy)*akt/LOG(zxih)* &
                      SQRT(-zris(i,j)*zxih) ))

          ! Use ice surface roughness or open-water surface roughness 
          ! according to l_ls_ice
          IF ( l_ls_ice ) THEN
            ! set z0 to 0.001m over sea ice
            gz0(i,j)= 0.001*g
          ELSE
            ! calculate z0 over water by Charnock formula
            gz0(i,j) = zalpha0 * ( ztcm(i,j)*zvpb(i,j) + zalphaf**2 &
                       *( ztch(i,j)*zx(i,j) )**zd3 )
          ENDIF
        ENDIF

      ENDIF

      ! Store TCM, TCH, GZ0 and the skin layer values of TKVH, TKVM and TKE

      tcm  (i,j) = ztcm(i,j)/zvpb(i,j)
      tch  (i,j) = ztch(i,j)/zvpb(i,j)
      gz0  (i,j) = MAX( gz0 (i,j), 1.0E-10_ireals )

    ENDDO
  ENDDO

  SELECT CASE( itype_turb )
  CASE( 1:2,7:8 )

    ! Store the skin layer values of TKVH, TKVM and TKE

    ntke = 1
    DO j = jstartpar , jendpar
      DO i = istartpar , iendpar
        IF (itype_turb <= 2) THEN
           tke(i,j,ke1,ntke) = SQRT(2.0*ztcm(i,j)*zvpb(i,j))/zctke
           edr(i,j,ke1) = tke(i,j,ke1,ntke)**3/(d_mom*akt*gz0(i,j)/g)
        END IF

        zdz=gz0(i,j)/g*LOG(zdfip(i,j)/gz0(i,j)+1.0)
        tkvm(i,j,ke1)=ztcm(i,j)*zdz
        tkvh(i,j,ke1)=ztch(i,j)*zdz
      ENDDO
    ENDDO
  END SELECT

!------------------------------------------------------------------------------
! End of module procedure parturs
!------------------------------------------------------------------------------

END SUBROUTINE parturs


!------------------------------------------------------------------------------
!+ Module procedure prankolmo_rk in "src_turbulence" for computing the 
!+ coefficients for vertical diffusion
!------------------------------------------------------------------------------

SUBROUTINE prankolmo_rk( lstfnct )

!------------------------------------------------------------------------------
!
! Description:
!   This module procedure calculates the coefficients TKVM and TKVH for 
!   vertical turbulent diffusion of momentum and heat, respectively following
!   an approach by Prandtl and Kolmogorov.
!   TKVM and TKVH are defined on half levels.
!
! Method:
!   The diffusion coefficient are assumed to be depending on the turbulent
!   kinetic energy following an approach by Prandtl and Kolmogorov.
!   The stability functions are specified either by use of Smagorinski
!   stability functions (100m-scale) or by those of the modified 
!   MELLOR & YAMADA-scheme (partura). The impact of subgrid
!   condensation/evaporation is not taken into account in this version.
!
! Subroutine history:
! S-Version Date       Comment
! -------   ----       -------
! 2.16      2003/01/29 Gerd Vogel
!                      Code created 
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!
!------------------------------------------------------------------------------
!
! Declarations:
!
!------------------------------------------------------------------------------

! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments: None
! --------------------
LOGICAL, INTENT (IN) ::      &
  lstfnct                 ! calculate stability functions - flag

!
! Local scalars:
! -------------
INTEGER (KIND=iintegers) ::  &
  i, j, k,              & ! loop indices
  izstata,              & ! error status at allocation
  izstatd,              & ! error status at deallocation
  izerror,              & !
  ms                      !

REAL    (KIND=ireals   ) ::  &
  zs11, zs22, zs33,     & !
  zs12, zs21,           & !
  zs13, zs31,           & !
  zs23, zs32,           & !
  zshear, zbuoy, zdiss, & !
  zpi, zdza, ztempr,    & !
  zth_hl, zdthdz,       & !
  zlhvocp, zlhsocp,     & !
  zqc_hl, zqis_hl,      & !
  zlhv_qv, zlhs_qv,     & !
  za_buoy, zdqxdz,      & !
  zlenq

! constants used in both parameteterization schemes
REAL    (KIND=ireals   ) ::  &
  zphim,                & !
  za,                   & !
  zdvbmin,              & !
  zric,                 & !
  z1d3, z2d3,           & !
  zceps, zceps1d3,      & !
  ztkmmin,              & ! mimimum absolut value TKVM
  ztkhmin                 ! minimum alsolut value TKVH

! constants for diffusion parameterization (MELLOR & YAMADA)
REAL    (KIND=ireals   ) ::  &
  zmzb  ,               & !
  zdze  ,               & !
  zrf   ,               & !
  zgam  ,               & !
  zsh   ,               & !
  zsm   ,               & !
  zalhn ,               & !
  zcepsm2d3,            & !
  ztmmin,               & !
  zthmin,               & !
  zalf, zemalf,         & !
  ztkvmom,ztkvhom         !

! constants for diffusion  parameterization (Smagorinski)
REAL    (KIND=ireals   ) ::  &
  zd,                   & ! Achtung: zd von aussen setzen
  zfm, zfm1, zfh,       & !
  zkh,                  & !
  zkhdzd, zfdzd,        & !
  zceps2d3,             & !
  zd4d3,                & !
  zprn, zb, zc,         & !
  zvw, zg, zh, zr,      & !
  zcs, zcs2, zcs4d3       !

CHARACTER (LEN=80)            :: yzerrmsg
CHARACTER (LEN=25)            :: yzroutine

! Local (automatic) arrays:
! -------------------------
REAL    (KIND=ireals   ) ::  &
  zlen    (ie,je,2:ke), & ! turbulent length scale
  zsgesh  (ie,je,2:ke), & !
  zsgesv  (ie,je,2:ke), & !
  ztpm    (ie,je,2:ke), & ! squared 3D-deformation
  ztpm_v  (ie,je,2:ke), & ! squared vertical 1D-deformation
  ztph    (ie,je,2:ke), & ! squared Brunt-Vaisaelae-frequency
  zrinum  (ie,je,2:ke), & ! 3D-Richardson number: N^2/S^2
  zrinum_v(ie,je,2:ke), & ! 1D-Richardson number: N^2/S_V^2
  zaa     (2),          & !
  zab1    (2),          & !
  zab2    (2),          & !
  zac1    (2),          & !
  zac2    (2),          & !
  zac3    (2)             !

! Local allocatable arrays:
! -------------------------
REAL(KIND=ireals), ALLOCATABLE ::  &
  zuvw_i  (:,:,:),     & ! 
  zuvw_j  (:,:,:),     & !
  ztheta  (:,:,:),     & !
  ztheta_v(:,:,:),     & !
  ztheta_l(:,:,:),     & !
  ztheta_i(:,:,:),     & !
  zqx_sum (:,:,:)        !

! Tracer pointers
!----------------
REAL (KIND=ireals), POINTER :: &
  qv  (:,:,:)=> NULL(),         & ! QV at tlev=nnow
  qc  (:,:,:)=> NULL()            ! QC at tlev=nnow

!------------ End of header ---------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine prankolmo_rk
!------------------------------------------------------------------------------

  yzroutine = 'prankolmo_rk'
  
  zdvbmin  = 0.01_ireals

  ALLOCATE( zuvw_i(ie,je,ke), zuvw_j(ie,je,ke), STAT=izstata )

  ! retrieve the required microphysics tracers (at nnow)
  CALL trcr_get(izerror, idt_qv, ptr_tlev = nnow, ptr = qv)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qc, ptr_tlev = nnow, ptr = qc)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

  IF ( l3dturb ) THEN ! 3D-turbulence
    
    DO k = 2, ke
      DO j = jstart, jend
        DO i = istart, iend
          zs11 = ( u(i,j,k  ,nnow) - u(i-1,j,k  ,nnow)   +       &
                   u(i,j,k-1,nnow) - u(i-1,j,k-1,nnow) )         &
                   * 0.5_ireals * eddlon * acrlat(j,1)
          zs22 = ( v(i,j,k  ,nnow) - v(i,j-1,k  ,nnow)   +       &
                   v(i,j,k-1,nnow) - v(i,j-1,k-1,nnow) )         &
                   * 0.5_ireals * edadlat
          zs33 = ( w(i,j,k+1,nnow) - w(i,j,k-1,nnow) )           &
                   / ( hhl(i,j,k-1) - hhl(i,j,k+1) )
          zsgesh(i,j,k) = 2.0_ireals * ( zs11**2 + zs22**2 )
          zsgesv(i,j,k) = 2.0_ireals * zs33**2
        END DO
      END DO
    END DO

    DO k = 2, ke
      DO j = jstart-1, jend
        DO i = istart-1, iend
          zuvw_i(i,j,k) = ( u(i,j  ,k  ,nnow) + u(i-1,j  ,k  ,nnow) +     &
                            u(i,j+1,k  ,nnow) + u(i-1,j+1,k  ,nnow) +     &
                            u(i,j  ,k-1,nnow) + u(i-1,j  ,k-1,nnow) +     &
                            u(i,j+1,k-1,nnow) + u(i-1,j+1,k-1,nnow) )     &
                            * 0.125_ireals
          zuvw_j(i,j,k) = ( v(i  ,j,k  ,nnow) + v(i  ,j-1,k  ,nnow) +     &
                            v(i+1,j,k  ,nnow) + v(i+1,j-1,k  ,nnow) +     &
                            v(i  ,j,k-1,nnow) + v(i  ,j-1,k-1,nnow) +     &
                            v(i+1,j,k-1,nnow) + v(i+1,j-1,k-1,nnow) )     &
                            * 0.125_ireals
        END DO
      END DO
    END DO
    DO k = 2, ke
      DO j = jstart, jend
        DO i = istart, iend
          zs12 = ( zuvw_i(i,j,k) - zuvw_i(i,j-1,k) ) * edadlat
          zs21 = ( zuvw_j(i,j,k) - zuvw_j(i-1,j,k) ) * acrlat(j,1) * eddlon
          zsgesh(i,j,k) = zsgesh(i,j,k) + ( zs12 + zs21 )**2
        END DO
      END DO
    END DO

    DO k = 1, ke
      DO j = jstart, jend
        DO i = istart-1, iend
          zuvw_i(i,j,k) = 0.5_ireals * ( u(i,j,k,nnow) + u(i-1,j,k,nnow) )
          zuvw_j(i,j,k) = 0.5_ireals * ( w(i,j,k,nnow) + w(i+1,j,k,nnow) )
        END DO
      END DO
    END DO
    DO k = 2, ke
      DO j = jstart, jend
        DO i = istart, iend
          zs13 = - 2.0_ireals * ( zuvw_i(i,j,k) - zuvw_i(i,j,k-1) )       &
                              / ( hhl(i,j,k-1)-hhl(i,j,k+1) )
          zs31 = ( zuvw_j(i,j,k) - zuvw_j(i-1,j,k) ) * acrlat(j,1) * eddlon
          zsgesv(i,j,k) = zsgesv(i,j,k) + ( zs13 + zs31 )**2
        END DO
      END DO
    END DO

    DO k = 1, ke
      DO j = jstart-1, jend
        DO i = istart, iend
          zuvw_i(i,j,k) = 0.5_ireals * ( v(i,j,k,nnow) + v(i,j-1,k,nnow) )
          zuvw_j(i,j,k) = 0.5_ireals * ( w(i,j,k,nnow) + w(i,j+1,k,nnow) )
        END DO
      END DO
    END DO
    DO k = 2, ke
      DO j = jstart, jend
        DO i = istart, iend
          zs23 = - 2.0_ireals * ( zuvw_i(i,j,k) - zuvw_i(i,j,k-1) )       &
                              / ( hhl(i,j,k-1)-hhl(i,j,k+1) )
          zs32 = ( zuvw_j(i,j,k) - zuvw_j(i,j-1,k) ) * edadlat
          zsgesv(i,j,k) = zsgesv(i,j,k) + ( zs23 + zs32 )**2

          ! mechanical production of TKE
!!$ UB          ztpm(i,j,k) = ( zsgesh(i,j,k) + zsgesv(i,j,k) )**2
!!$ UB          ztpm_v(i,j,k) = zsgesv(i,j,k)**2
          ztpm  (i,j,k) = zsgesh(i,j,k) + zsgesv(i,j,k)
          ztpm_v(i,j,k) = zsgesv(i,j,k)
        END DO
      END DO
    END DO

  ELSE ! 1D-turbulence
    
    DO k = 1, ke
      DO j = jstart, jend
        DO i = istart, iend
          zuvw_i(i,j,k) = 0.5_ireals * ( u(i,j,k,nnow) + u(i-1,j,k,nnow) )
          zuvw_j(i,j,k) = 0.5_ireals * ( v(i,j,k,nnow) + v(i,j-1,k,nnow) )
        END DO
      END DO
    END DO
    DO k = 2, ke
      DO j = jstart, jend
        DO i = istart, iend
          zs13 = - 2.0_ireals * ( zuvw_i(i,j,k) - zuvw_i(i,j,k-1) )       &
                              / ( hhl(i,j,k-1)-hhl(i,j,k+1) )
          zs23 = - 2.0_ireals * ( zuvw_j(i,j,k) - zuvw_j(i,j,k-1) )       &
                              / ( hhl(i,j,k-1)-hhl(i,j,k+1) )
          zsgesv(i,j,k) = zs13**2 + zs23**2

          ! mechanical production of TKE
!!$ UB          ztpm(i,j,k) = zsgesv(i,j,k)**2
          ztpm  (i,j,k) = zsgesv(i,j,k)
          ztpm_v(i,j,k) = ztpm  (i,j,k)
        
        END DO
      END DO
    END DO

  END IF
  
  DEALLOCATE( zuvw_i, zuvw_j, STAT=izstatd )

  
  ALLOCATE( ztheta(ie,je,ke), ztheta_v(ie,je,ke), STAT=izstata )
  
  DO k = 1, ke
    DO j = jstart, jend
      DO i = istart, iend
        zpi = ( 1.0E-5_ireals * ( p0(i,j,k) + pp(i,j,k,nnow) ) )**rdocp
        ! potential temperature
        ztheta(i,j,k)   = t(i,j,k,nnow) / zpi
        ! virtual potential temperature
        ztheta_v(i,j,k) = ztheta(i,j,k) * ( 1.0_ireals +                 &
                          rvd_m_o*qv(i,j,k) - qc(i,j,k) - qrs(i,j,k) )
      END DO
    END DO
  END DO

  
  SELECT CASE( itype_turb )
    
  CASE( 5, 7 )

    DO k = 2, ke
      DO j = jstart, jend
        DO i = istart, iend
          
          zdza   = 0.5_ireals * ( hhl(i,j,k+1) - hhl(i,j,k-1) )
          zth_hl = 0.5_ireals * ( ztheta_v(i,j,k) + ztheta_v(i,j,k-1) )
          zdthdz = ztheta_v(i,j,k) - ztheta_v(i,j,k-1)
          
          ! buoyant production of TKE
          ztph(i,j,k) = g / zth_hl * zdthdz / zdza

          ! set minimum values for mechanical production of TKE
          ztpm(i,j,k) = MAX( ztpm(i,j,k), ( zdvbmin/zdza )**2 )
          ztpm_v(i,j,k) = MAX( ztpm_v(i,j,k), ( zdvbmin/zdza )**2 )

          ! Richardson number
          zrinum(i,j,k)   = ztph(i,j,k) / ztpm(i,j,k)
          zrinum_v(i,j,k) = ztph(i,j,k) / ztpm_v(i,j,k)

        END DO
      END DO
    END DO
    
  CASE( 6, 8 )
    
    ALLOCATE( ztheta_l(ie,je,ke), zqx_sum(ie,je,ke), STAT=izstata )
! JF:     IF ( itype_gscp >= 3 ) ALLOCATE( ztheta_i(ie,je,ke), STAT=izstata )

    zlhvocp = lh_v / cp_d
! JF:     zlhsocp = lh_s / cp_d
  
    DO k = 1, ke
      DO j = jstart, jend
        DO i = istart, iend
          
! JF:           ztempr = 1.0_ireals / t(i,j,k,nnow)
          
          ! equivalent potential temperature
! JF:           ztheta_l(i,j,k) = ztheta(i,j,k)  &
! JF:                           * EXP( zlhvocp * qv(i,j,k) * ztempr )
          ! liquid water potential temperature
! JF:           ztheta_l(i,j,k) = ztheta(i,j,k)  &
! JF:                           * EXP( - zlhvocp * qc(i,j,k) * ztempr )

          ! Exner function
          zpi = ( 1.0E-5_ireals * ( p0(i,j,k) + pp(i,j,k,nnow) ) )**rdocp

          ! approx. equivalent potential temperature
!          ztheta_l(i,j,k) = ztheta(i,j,k) + zlhvocp * qv(i,j,k) / zpi

          ! approx. liquid water potential temperature
          ztheta_l(i,j,k) = ztheta(i,j,k) - zlhvocp * qc(i,j,k) / zpi

          ! sum of liquid and solid substances
          zqx_sum(i,j,k) = qc(i,j,k) + qrs(i,j,k)
          
        END DO
      END DO
    END DO

! JF:     IF ( itype_gscp >= 3 ) THEN
! JF:       DO k = 1, ke
! JF:         DO j = jstart, jend
! JF:           DO i = istart, iend
! JF:             
! JF:             ztempr = 1.0_ireals / t(i,j,k,nnow)
! JF:             
! JF:             ! ice potential temperature
! JF:             ztheta_i(i,j,k) = ztheta(i,j,k)  &
! JF:                             * EXP( - zlhsocp * qi(i,j,k,nnow) * ztempr )
! JF:             
! JF:           END DO
! JF:         END DO
! JF:       END DO
! JF:     END IF

    DO k = 2, ke
      DO j = jstart, jend
        DO i = istart, iend
          
          zdza   = 0.5_ireals * ( hhl(i,j,k+1) - hhl(i,j,k-1) )
          zth_hl = 0.5_ireals * ( ztheta(i,j,k) + ztheta(i,j,k-1) )

          zqc_hl = 0.5_ireals * ( qc(i,j,k) + qc(i,j,k-1) )
          ! water cloud
          IF ( zqc_hl > 0.0_ireals ) THEN

            zlhv_qv = lh_v * qv(i,j,k)
            za_buoy = 1.0_ireals + zlhv_qv / ( r_d * t(i,j,k,nnow) )
            za_buoy = za_buoy / zth_hl /                               &
                    ( 1.0_ireals + zlhvocp * zlhv_qv /                 &
                    ( r_v * t(i,j,k,nnow)**2 ) )

            zdthdz = ztheta_l(i,j,k) - ztheta_l(i,j,k-1)
            zdqxdz = zqx_sum (i,j,k) - zqx_sum (i,j,k-1)

            ! buoyant production of TKE
            ztph(i,j,k) = g * ( za_buoy * zdthdz - zdqxdz ) / zdza
!            ztph(i,j,k) = g * ( zdthdz / zth_hl - zdqxdz ) / zdza
            
          ELSE

! JF:             IF ( itype_gscp >= 3 ) THEN
! JF:               zqis_hl = 0.5_ireals * ( qi(i,j,k,nnow)   + qs(i,j,k,nnow)      &
! JF:                                      + qi(i,j,k-1,nnow) + qs(i,j,k-1,nnow) )
! JF:             ELSE
! JF:               zqis_hl = 0.0_ireals
! JF:             END IF
! JF:             
! JF:             ! ice cloud
! JF:             IF ( zqis_hl > 0.0_ireals ) THEN
! JF: 
! JF:               zlhs_qv = lh_s * qv(i,j,k)
! JF:               za_buoy = 1.0_ireals + zlhs_qv / ( r_d * t(i,j,k,nnow) )
! JF:               za_buoy = za_buoy / zth_hl /                               &
! JF:                       ( 1.0_ireals + zlhsocp * zlhs_qv /                 &
! JF:                       ( r_v * t(i,j,k,nnow)**2 ) )
! JF: 
! JF:               zdthdz = ztheta_i(i,j,k) - ztheta_i(i,j,k-1)
! JF:               zdqxdz = zqx_sum (i,j,k) - zqx_sum (i,j,k-1)
! JF: 
! JF:               ! buoyant production of TKE
! JF:               ztph(i,j,k) = g * ( za_buoy * zdthdz - zdqxdz ) / zdza
! JF: !              ztph(i,j,k) = g * ( zdthdz / zth_hl - zdqxdz ) / zdza
! JF:             
! JF:             ELSE
              
              zdthdz = ztheta_v(i,j,k) - ztheta_v(i,j,k-1)
          
              ! buoyant production of TKE
              ztph(i,j,k) = g / zth_hl * zdthdz / zdza

! JF:             END IF

          END IF
          
          ! set minimum values for mechanical production of TKE
          ztpm(i,j,k) = MAX( ztpm(i,j,k), ( zdvbmin/zdza )**2 )
          ztpm_v(i,j,k) = MAX( ztpm_v(i,j,k), ( zdvbmin/zdza )**2 )

          ! Richardson number
          zrinum(i,j,k)   = ztph(i,j,k) / ztpm(i,j,k)
          zrinum_v(i,j,k) = ztph(i,j,k) / ztpm_v(i,j,k)

        END DO
      END DO
    END DO

    DEALLOCATE( ztheta_l, zqx_sum, STAT=izstatd )
! JF:     IF ( itype_gscp >= 3 ) DEALLOCATE( ztheta_i, STAT=izstatd )

  END SELECT

  DEALLOCATE( ztheta, ztheta_v, STAT=izstatd )


  IF ( lprog_tke ) THEN
    ntke = nnow
  ELSE
    ntke = 1
  END IF


  IF ( lstfnct ) THEN
    
    z1d3 = 1.0_ireals/3.0_ireals
    z2d3 = 2.0_ireals/3.0_ireals

    SELECT CASE( itype_turb )

    CASE( 5:6 )
      !
      ! diffusion approach using the stability functions of the
      ! modified MELLOR & YAMADA scheme (Mueller...)
      !

      ! basic constants of the parameterization scheme
      zric    = 0.380_ireals      !
      zalhn   = 1.00_ireals       !
      zceps   = 0.060241_ireals   ! = 1 / 16.6
      zcepsm2d3 = zceps**(-z2d3)  !
      zceps1d3 = zceps**z1d3      !
      zdze    = 500.0_ireals      !     
      ztmmin  = 0.010_ireals      ! min value as part of neutral value for TKVM
      zthmin  = 0.007_ireals      ! min value as part of neutral value for TKVH
      ztkmmin = 1.000_ireals      ! mimimum absolut value TKVM
      ztkhmin = 0.400_ireals      ! minimum alsolut value TKVH

      ! zalf    = 0.9_ireals        !
      ! zemalf  = 0.1_ireals        !
      zalf    = 1.0_ireals        !
      zemalf  = 0.0_ireals        !
      zaa(1)  = 3.7_ireals        !
      zaa(2)  = 4.025_ireals      !
      zab1(1) = 2.5648_ireals     !
      zab1(2) = 3.337_ireals      !
      zab2(1) = 1.1388_ireals     !
      zab2(2) = 0.688_ireals      !
      zac1(1) = 0.8333_ireals     !
      zac1(2) = 1.285_ireals      !
      zac2(1) = 0.2805_ireals     !
      zac2(2) = 0.2305_ireals     !
      zac3(1) = 0.1122_ireals     !
      zac3(2) =-0.1023_ireals     !

      ! Start of calculation for each layer from top to bottom  
      ! ------------------------------------------------------

      DO k = 2, ke

        DO j = jstart, jend
          DO i = istart, iend

            zmzb = hhl(i,j,k) - hhl(i,j,ke1)
            zlen(i,j,k) = akt*zmzb / ( 1.0_ireals + akt*zmzb/zdze )
            zlenq = zlen(i,j,k) * zlen(i,j,k)

            ! very stable case (ri>ric)
            IF ( zrinum_v(i,j,k) >= zric ) THEN

              tkvm(i,j,k) = ztmmin * zlenq * SQRT( ztpm_v(i,j,k) )
              tkvh(i,j,k) = zthmin * zlenq * SQRT( ztpm_v(i,j,k) )
              IF ( .NOT.lprog_tke .OR. ntstep == 0 ) THEN
                tke(i,j,k,ntke) = 0.0_ireals
              ENDIF

            ELSE

              IF ( zrinum_v(i,j,k) > 0.0_ireals ) THEN
                ! stable case
                ms = 1
              ELSE
                ! unstable case
                ms = 2
              ENDIF

              ! Flux-Richardson number as function
              ! of the Gradient-Richardson number
              zrf  = zac1(ms) * ( zrinum_v(i,j,k) + zac2(ms) -                &
                     SQRT( zrinum_v(i,j,k)**2 - zac3(ms)*zrinum_v(i,j,k)      &
                                            + zac2(ms)**2) )
              zrf  = MIN ( zrf, 0.99999_ireals )

              ! Gam(rf)
              zgam = zrf / (1.0_ireals - zrf)

              ! Sm(Gam(rf)) and Sh(Gam(rf))
              zsh = ( 1.0_ireals - zab1(ms)*zgam )                            &
                  / ( 1.0_ireals - zab2(ms)*zgam )
              zsm = ( 1.0_ireals - zaa(ms)*zgam ) / zsh

              !
              ! turbulent kinetic energy and diffusion coefficients
              !
              IF ( .NOT.lprog_tke .OR. ntstep == 0 ) THEN
                tke(i,j,k,ntke) = ( 1.0_ireals - zalhn*zsh*zrinum_v(i,j,k) )  &
                                * zcepsm2d3 * zsm * zlenq * ztpm_v(i,j,k)
                tke(i,j,k,ntke) = MAX( 0.0_ireals, tke(i,j,k,ntke) )
              ENDIF
              
              zphim   = zceps1d3 * zsm
              ztkvmom = zphim * zlen(i,j,k) * SQRT( tke(i,j,k,ntke) )
              ztkvhom = zalhn * zsh * ztkvmom

              tkvm(i,j,k) = zalf*ztkvmom + zemalf*tkvm(i,j,k) 
              tkvh(i,j,k) = zalf*ztkvhom + zemalf*tkvh(i,j,k) 

            ENDIF

            ! Set lower limit for diffusion coefficients
            tkvm(i,j,k) = MAX( ztkmmin, tkvm(i,j,k) )
            tkvh(i,j,k) = MAX( ztkhmin, tkvh(i,j,k) )

          ENDDO
        ENDDO

      END DO

    CASE( 7:8 )
      !
      ! diffusion approach using stability functions of Smagorinski
      !

      ! Bug Fix by Uli Blakah
      ! from fixed value to real model value of horizontal grid spacing
      ! (96.5 m was from original LLM-setup):
      ! zd = 96.5_ireals ! Achtung: von aussen setzen!!!!
      zd  =  SQRT( dlat * dlon * r_earth * r_earth * degrad * degrad * &
           COS( (startlat_tot+(INT(je_tot/2)+1)*dlat)*degrad) )

      zric     = 0.25_ireals
      zprn     = 0.7_ireals
      zb       = 40.0_ireals
      zc       = 16.0_ireals
      zvw      = 0.5_ireals
      zg       = 1.2_ireals
      zh       = 0.0_ireals
      zr       = 4.0_ireals
      ! UB>> Better results with zcs=0.15, together with tkhx = 3 * tkvx
      !      zcs      = 0.25_ireals
      zcs      = 0.15_ireals
      zceps    = 0.93_ireals 
      za       = 1.0 / zprn    
      ztkmmin  = 0.10_ireals  ! mimimum absolut value TKVM
      ztkhmin  = 0.10_ireals  ! minimum absolut value TKVH
      zd4d3    = zd**(4.0_ireals/3.0_ireals)
      zcs2     = zcs * zcs
      zcs4d3   = zcs**(4.0_ireals/3.0_ireals)
      zceps2d3 = zceps**z2d3
      zceps1d3 = zceps**z1d3

      ! Start of calculation for each layer from top to bottom  
      ! ------------------------------------------------------

      DO k = 2, ke

        DO j = jstart, jend
          DO i = istart, iend

            zkh = 0.5_ireals*( hhl(i,j,k-1)-hhl(i,j,k+1) )
            zkhdzd  = zkh / zd
            IF (zkhdzd > 1.0_ireals) THEN
              zkhdzd  = zd / zkh
            ENDIF
            IF (zkhdzd == 1.0_ireals) THEN
              zfdzd = 1.0_ireals
            ELSE
              zfdzd = COSH( 0.3849_ireals * LOG(zkhdzd) )
            ENDIF
            zlenq = zfdzd * zfdzd * zd4d3 * zkh**z2d3
            zlen(i,j,k)  = SQRT( zlenq )

            IF (zrinum(i,j,k) < 0.0_ireals) THEN
              ! unstable case
              zfm = ( 1.0_ireals - zc*zrinum(i,j,k) )**zvw
              zfh = za * ( 1.0_ireals - zb * zrinum(i,j,k) )**zvw
            ELSE IF (zrinum(i,j,k) >= zric) THEN
              zfm = 0.0_ireals
              zfh = 0.0_ireals
            ELSE
              ! stable case
              zfm1 = ( 1.0_ireals - zrinum(i,j,k)/zric )**zr
              zfm  = zfm1 * ( 1.0_ireals - zh * zrinum(i,j,k) )
              zfh  = za * zfm1 * ( 1.0_ireals - zg * zrinum(i,j,k) )
            ENDIF

            !
            ! turbulent kinetic energy and diffusion coefficients
            !
            IF (zrinum(i,j,k) >= zric) THEN
              
              zphim = 0.0_ireals
              IF ( ntstep == 0 ) tke(i,j,k,ntke) = 0.0_ireals
              tkvm(i,j,k) = 0.0_ireals
              tkvh(i,j,k) = 0.0_ireals
              
            ELSE
              
              zphim  = ( zfm**z2d3 ) * zcs4d3 * zceps1d3 /                    &
                       ( ( 1.0 - zfh/zfm * zrinum(i,j,k))**z1d3 )
              IF ( .NOT.lprog_tke .OR. ntstep == 0) THEN
                tke(i,j,k,ntke) = ( zfm - zfh*zrinum(i,j,k) )**z2d3           &
                                * ztpm(i,j,k) * zlenq * zcs4d3 / zceps2d3
                tke(i,j,k,ntke) = MAX( 0.0_ireals, tke(i,j,k,ntke) )
              ENDIF
              tkvm(i,j,k) = zphim * zlen(i,j,k) * SQRT( tke(i,j,k,ntke) )
              tkvh(i,j,k) = tkvm(i,j,k) * zfh/zfm
              
            ENDIF

            ! Set lower limit for diffusion coefficients
            tkvm(i,j,k) = MAX( ztkmmin, tkvm(i,j,k) )
            tkvh(i,j,k) = MAX( ztkhmin, tkvh(i,j,k) )

          ENDDO
        ENDDO

      END DO

    END SELECT

  
    ! Special treatment at the southern boundary of the total domain
    IF (my_cart_neigh(4) == -1) THEN
      DO k = 2, ke
        DO j = 1, jstart-1
          DO i = istartpar, iendpar
            tke (i,j,k,ntke) = tke (i,jstart,k,ntke)
            tkvm(i,j,k)      = tkvm(i,jstart,k)
            tkvh(i,j,k)      = tkvh(i,jstart,k)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    ! Special treatment at the northern boundary of the total domain
    IF (my_cart_neigh(2) == -1) THEN
      DO k = 2, ke
        DO j = jend+1, je
          DO i = istartpar, iendpar
            tke (i,j,k,ntke) = tke (i,jend,k,ntke)
            tkvm(i,j,k)      = tkvm(i,jend,k)
            tkvh(i,j,k)      = tkvh(i,jend,k)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    ! Special treatment at the western boundary of the total domain
    IF (my_cart_neigh(1) == -1) THEN
      DO k = 2, ke
        DO j = 1, jendpar
          DO i = 1, istart-1
            tke (i,j,k,ntke) = tke (istart,j,k,ntke)
            tkvm(i,j,k)      = tkvm(istart,j,k)
            tkvh(i,j,k)      = tkvh(istart,j,k)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    ! Special treatment at the eastern boundary of the total domain
    IF (my_cart_neigh(3) == -1) THEN
      DO k = 2, ke
        DO j = 1, jendpar
          DO i = iend+1, ie
            tke (i,j,k,ntke) = tke (iend,j,k,ntke)
            tkvm(i,j,k)      = tkvm(iend,j,k)
            tkvh(i,j,k)      = tkvh(iend,j,k)
          ENDDO
        ENDDO
      ENDDO
    ENDIF

  END IF
  
! JF:   ! Bottom boundary condition for TKE: v. Neumann BC
! JF:   DO j = jstartpar, jendpar
! JF:     DO i = istartpar, iendpar
! JF:       tke(i,j,ke1,nnew) = tke(i,j,ke,nnew)
! JF:     ENDDO
! JF:   ENDDO
      
  
  IF ( lprog_tke ) THEN

    tketens(:,:,1) = 0.0_ireals
    
    !
    ! Add tke-tendencies due to (vertical) shear, buoyancy and dissipation
    !
    DO k = 2, ke
      DO j = jstart, jend
        DO i = istart, iend
          
          zshear =   tkvm(i,j,k) * zsgesv(i,j,k)
          zbuoy  = - tkvh(i,j,k) * zrinum(i,j,k) * ztpm(i,j,k)
          zdiss  = - zceps / zlen(i,j,k) * tke(i,j,k,ntke)**1.5_ireals

          tketens(i,j,k) = zshear + zbuoy + zdiss
          
        ENDDO
      ENDDO
    ENDDO
    
    IF ( l3dturb ) THEN

      ! Compute horizontal diffusion coefficients
      CALL horizontal_diffcoeffs
      
      ! Add tke-tendency due to horizontal shear
      
      DO k = 2, ke
        DO j = jstart, jend
          DO i = istart, iend
            
            zshear = tkhm(i,j,k) * zsgesh(i,j,k)

            tketens(i,j,k) = tketens(i,j,k) + zshear
            
          ENDDO
        ENDDO
      ENDDO

    END IF
    
  ENDIF
  
!------------------------------------------------------------------------------
! End of module procedure prankolmo_rk
!------------------------------------------------------------------------------

END SUBROUTINE prankolmo_rk

!==============================================================================

SUBROUTINE horizontal_diffcoeffs

  ! Local scalars:
  ! -------------
  INTEGER (KIND=iintegers) ::  &
       i, j, k              ! loop indices

  REAL    (KIND=ireals   ) ::  &
       zdelt1, zdelt2,    & !
       za_dhodz, zdhodz

  ! Local (automatic) arrays:
  ! -------------------------
  REAL    (KIND=ireals   ) ::  &
       zdelt   (je)


  ! JF:   za_dhodz = 0.1_ireals
  ! JF:   
  ! JF:   DO j = jstart-1, jend+1
  ! JF:     zdelt1 = r_earth*crlat(j,1)*degrad*dlon
  ! JF:     zdelt2 = r_earth*degrad*dlat
  ! JF:     zdelt(j) = SQRT( zdelt1*zdelt1 + zdelt2*zdelt2 )
  ! JF:   ENDDO

  !
  ! Compute horizontal diffusion coefficients
  !

  SELECT CASE (itype_turb)

  CASE (7:8)

    DO k = 2, ke
      DO j = jstart-1, jend+1
        DO  i = istart-1, iend+1

          ! UB>>
          ! Better results with tkhx = 3 * tkvx, 
          ! together WITH zcs=0.15 in CASE of itype_turb=7/8
          tkhm(i,j,k) = 3.0 * tkvm(i,j,k)
          tkhh(i,j,k) = 3.0 * tkvh(i,j,k)

        ENDDO
      ENDDO
    END DO

  CASE default

    DO k = 2, ke
      DO j = jstart-1, jend+1
        DO  i = istart-1, iend+1

          ! JF:         zdhodz = za_dhodz * zdelt(j) * sqrtg_r_w(i,j,k)
          ! JF:
          ! JF:         tkhm(i,j,k) = zdhodz * tkvm(i,j,k)
          ! JF:         tkhh(i,j,k) = zdhodz * tkvh(i,j,k)

          ! This has to be further tuned in the future!
          tkhm(i,j,k) = 1.0 * tkvm(i,j,k)
          tkhh(i,j,k) = 1.0 * tkvh(i,j,k)

        ENDDO
      ENDDO
    END DO


  END SELECT


END SUBROUTINE horizontal_diffcoeffs

END MODULE src_turbulence
