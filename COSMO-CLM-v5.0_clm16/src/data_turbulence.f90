!+ Data module for variables of the turbulence parameterization
!------------------------------------------------------------------------------

MODULE data_turbulence

!------------------------------------------------------------------------------
!
! Description:
!  This module contains variables that are used in the turbulence
!  parameterizations. With these variables a tuning of the scheme is
!  possible.
!
! Current Code Owner: DWD, Matthias Raschendorfer
!  phone:  +49  69  8062 2708
!  fax:    +49  69  8062 3721
!  email:  Matthias.Raschendorfer@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 3.21       2006/12/04 Ulrich Schaettler
!  Initial release
! V3_23        2007/03/30 Matthias Raschendorfer
!  Importing 'rat_lam' from data_soil
!  and 'clc_diag', 'q_crit', 'akt' from data_constants.
!  Introduction of some parameters from turb_param.incf.
!  Initialization of all parameters with default values.
! V4_10        2009/09/11 Matthias Raschendorfer
!  Introduction of 'a_hshr' and 'a_stab'.
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_15        2010/11/19 Oliver Fuhrer
!  Reduced maximal value of securi to 0.5 because of possible numerical
!  instabilities otherwise
! V4_18        2011/05/26 Ulrich Schaettler
!  Changed the code owner
! V4_20        2011/08/31 Matthias Raschendorfer
!  Introduction of INTEGER parameter 'it_end'
!  Introduction of additional declarations (for turb_param)  (Uli Schaettler)
! V4_27        2013/03/19 Ulrich Schaettler
!  Modified default values of some tuning constants to reflect settings of COSMO-EU
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
USE data_parameters, ONLY :   &
    ireals,    & ! KIND-type parameter for real variables
    iintegers    ! KIND-type parameter for standard integer variables

!==============================================================================

IMPLICIT NONE

!==============================================================================


! Variables for tuning the turbulence parameterizations
! -----------------------------------------------------

! Attention:
! The given initializations are default settings of the boundary layer 
! parameters. Some of these initial parameter values may be changed afterwards
! by model input NAMELISTs!

!==============================================================================

! 1. Parameters describing physical properties of the lower boundary 
!    of the atmosphere:
!------------------------------------------

REAL (KIND=ireals) ::          &
  rlam_mom   =  0.0_ireals,    & ! scaling factor of the laminar boundary layer for momentum
  rlam_heat  =  1.0_ireals,    & ! scaling factor of the laminar boundary layer for heat

  rat_lam    =  1.0_ireals,    & ! ratio of laminar scaling factors for vapour and heat
  rat_can    =  1.0_ireals,    & ! ratio of canopy height over z0m
  rat_sea    = 20.0_ireals,    & ! ratio of laminar scaling factors for heat over sea and land

  z0m_dia    =  0.2_ireals,    & ! roughness length of a typical synoptic station [m]

  alpha0     =  0.0123_ireals, & ! Charnock-parameter
  alpha1     =  0.0000_ireals    ! parameter scaling the molek. roughness of water waves


! 2. Parameters that should be external parameter fields being not yet 
!    available:
!------------------------------------------

REAL (KIND=ireals) ::          &
  c_lnd      = 2.0_ireals,     & ! surface area density of the roughness elements over land
  c_sea      = 1.5_ireals,     & ! surface area density of the waves over sea
  c_soil     = 1.0_ireals,     & ! surface area density of the (evaporative) soil surface
  e_surf     = 1.0_ireals        ! exponent to get the effective surface area


! 3. Parameters that should be dynamical fields being not yet available:
!------------------------------------------

REAL (KIND=ireals) ::          &
  zt_ice     = -1.7_ireals,    & !freezing temperature of sea ice
  z0_ice     =  0.001_ireals     !roughness length of sea ice


! 4. Parameters for modelling turbulent diffusion:
!------------------------------------------

REAL (KIND=ireals) ::         &
  tur_len    = 500.0_ireals,  & ! assymtotic maximal turbulent length scale [m]
  pat_len    = 500.0_ireals,  & ! lenth scale of subscale surface patterns over land [m]
                                ! (should be dependent on location)
  len_min    =  1.0E-6_ireals,& ! minimal turbulent length scale [m]

  vel_min    =  0.01_ireals,  & ! minimal velocity scale [m/s]

  akt        =  0.4_ireals,   & ! von Karman-constant

  ! Length scale factors for pressure destruction of turbulent
  a_heat     =  0.74_ireals,  & ! scalar (heat) transport
  a_mom      =  0.92_ireals,  & ! momentum transport

  ! Length scale factors for dissipation of
  d_heat     =  10.1_ireals,  & ! scalar (temperature) variance
  d_mom      =  16.6_ireals,  & ! momentum variance

  ! Length scale factors for turbulent transport (vertical diffusion)
  c_diff     =  0.20_ireals,  & ! of TKE

  ! Length scale factor for separate horizontal shear production
  a_hshr     =  0.20_ireals,  & ! of TKE

  ! Length scale factor for the stability correction
  a_stab     =  0.00_ireals,  & ! no stability correction so far

  ! Dimensionless parameters used in the sub grid scale condensation scheme
  ! (statistical cloud scheme):
  clc_diag   =  0.5_ireals,   & !cloud cover at saturation
  q_crit     =  4.0_ireals,   & !critical value for normalized over-saturation
  c_scld     =  1.00_ireals,  & !factor for liquid water flux density in sub grid scale clouds

  ! Minimal diffusion coefficients in [m^2/s] for vertical
  tkhmin     =  0.4_ireals,   & ! scalar (heat) transport
  tkmmin     =  0.4_ireals      ! momentum transport


! 5. Numerical parameters:
!-------------------------

REAL (KIND=ireals) ::         &
  epsi       =  1.0E-6_ireals,& ! relative limit of accuracy for comparison of numbers
  tkesmot    =  0.15_ireals,  & ! time smoothing factor for TKE and diffusion coefficients
  wichfakt   =  0.00_ireals,  & ! vertical smoothing factor for explicit diffusion tendencies
  securi     =  0.50_ireals     ! security factor for maximal diffusion coefficients


! 6. Deduced constants for turbdiff and turbtran schemes
!-------------------------------------------------------

REAL (KIND=ireals), TARGET :: &
     c_tke,tet_g,c_g,rim, &
     d_0,d_1,d_2,d_3,d_4,d_5,d_6, &
     a_3,a_5,a_6,b_1,b_2, &
     l_scal,                      &
     l_hori                         ! horizontal grid spacing

INTEGER (KIND=iintegers) ::   &
  it_end     =  1               ! number of initialization iterations (>=0)

!==============================================================================

END MODULE data_turbulence
