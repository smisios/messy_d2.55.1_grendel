! **********************************************************************
!
! SUBMODEL CORE LAYER (SMCL) ROUTINES FOR MESSy SUBMODEL MLOCEAN_PLASIM
!
! Author : Markus Kunze, FUB-IFM, January-February  2011
!
! References:
!
! * P. Joeckel, R. Sander, A. Kerkweg, H. Tost, and J. Lelieveld,
!   Technical Note: The Modular Earth Submodel System (MESSy) - a new
!   approach towards Earth System Modeling,
!   Atmos. Chem. Phys., 5, 433-444, 2005.
!   http://www.atmos-chem-phys.net/5/433 
!
! **********************************************************************

! **********************************************************************
MODULE messy_mlocean_plasim
  ! **********************************************************************

  ! ----------- >
  
  USE messy_main_constants_mem, ONLY : dp  &
       , stbo       & ! Stephan-Boltzmann constant [W/m2/K4]
       , tmelt      & ! melting temperature of ice/snow
       , alf        & ! latent heat for fusion in J/kg (L_f)
       , rho_H2O    & ! density of H2O [kg/m3]
       , csw        & ! specific heat of sea waterJ/K/kg
       , ctfreez      ! temperature at which sea starts freezing/melting

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: dp

  ! ----------- <
  !
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'mlocean_plasim' ! submodel name
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '0.1'            ! submodel version
  !
  INTEGER, PARAMETER :: nlev_oce = 1         ! Number of Layers
  !
  ! CTRL-NAMELIST PARAMETERS  public for broadcast in SMIL
  INTEGER,  PUBLIC :: ndiag            = 480 ! diagnostics each ndiag timesteps
  INTEGER,  PUBLIC :: nout             = 32  ! afterburner output each nout timesteps
  INTEGER,  PUBLIC :: nocean           = 1   ! compute ocean yes/no
  INTEGER,  PUBLIC :: newsurf          = 0   ! update surface arrays at restart
  INTEGER,  PUBLIC :: nice             = 1   ! compute ice yes/no (1/0)
  INTEGER,  PUBLIC :: nsnow            = 1   ! allow snow on ice yes/no (1/0)
  INTEGER,  PUBLIC :: ntspd            = 32  ! ocean timesteps per day
  INTEGER,  PUBLIC :: nperpetual_ocean = 0   ! perpetual climate conditions
  INTEGER,  PUBLIC :: nprint           = 0   ! print debug information
  INTEGER,  PUBLIC :: nprhor           = 0   ! gp to print debug information
  !
  REAL(dp), PUBLIC :: dlayer(nlev_oce) = 50._dp   ! layer depth (m)
  REAL(dp), PUBLIC :: taunc            =  0._dp   ! newtonian cooling timescale (d)
  REAL(dp), PUBLIC :: vdiffk           = 1.E-4_dp ! vertikal diffusion coeff. [m**2/s]
  REAL(dp), PUBLIC :: xmind            = 0.1_dp   ! minimal ice thickness (m)
  !
  REAL(dp), PARAMETER :: zalpha    = 2.1656_dp ! thermal conductivity of ice in J m^-1 K^-1
  REAL(dp), PARAMETER :: zalphas   = 0.31_dp   ! thermal conductivity of snow in J m^-1 K^-1
  REAL(dp), PARAMETER :: zrho_sn   = 330._dp   ! density of snow in kg m^-3
  REAL(dp), PARAMETER :: zrhoice   = 910._dp   ! density of ice in kg m^-3
  REAL(dp), PARAMETER :: zrho_sea  = 1025._dp  ! density of sea water in kg m^-3

  REAL(dp), PARAMETER :: zcpice    = 2106._dp  ! specific heat of ice Ws kg^-1 K^-1
  REAL(dp), PARAMETER :: zdice     = 0.10_dp   ! minimum thickness of ice slab in meter
  ! -
  REAL(dp), PARAMETER :: CRHOS     = 1030._dp  ! Density of sea water (kg/m**3)
  REAL(dp), PARAMETER :: CRHOI     = 920._dp   ! Density of ice (kg/m**3)
  REAL(dp), PARAMETER :: CRHOF     = 1003.8_dp ! Density of fresh water AT S=5 (kg/m**3)
  REAL(dp), PARAMETER :: CRHOSN    = 330._dp   ! Density of snow
  REAL(dp), PARAMETER :: CPS       = 4180._dp  ! Specific heat of sea water (J/(kg*K))
  REAL(dp), PARAMETER :: CPI       = 2070._dp  ! Specific heat of ice  (J/(kg*K))
  REAL(dp), PARAMETER :: CPSN      = 2090._dp  ! Specific heat of snow (J/(kg*K))
  REAL(dp), PARAMETER :: CKAPI     = 2.03_dp   ! heat conductivity in ice (W/(m*K))
  REAL(dp), PARAMETER :: CKAPSN    = 0.31_dp   ! heat conductivity in snow (W/(m*K))
  REAL(dp), PARAMETER :: TFREEZE   = 271.25_dp ! Freezing point (K)
  REAL(dp), PARAMETER :: CLFI      = 3.28E5_dp ! heat of fusion of ice (J/kg)
  REAL(dp), PARAMETER :: CLFSN     = 3.32E5_dp ! heat of fusion of snow (J/kg)
  !
  ! PUBLIC SUBROUTINES (to be called from messy_mlocean_e5.f90)
  !
  PUBLIC :: mlocean_plasim
  PUBLIC :: mlocean_plasim_ice
  PUBLIC :: mlocean_plasim_sst

CONTAINS

  SUBROUTINE mlocean_plasim()
    !  ---------------------------------------------------------------------
    RETURN
  END SUBROUTINE mlocean_plasim
  ! =========================================================================
  SUBROUTINE mlocean_plasim_ice ()
    !  ---------------------------------------------------------------------
    RETURN
  END SUBROUTINE mlocean_plasim_ice
  ! =========================================================================
  SUBROUTINE mlocean_plasim_sst ()  
    !  ---------------------------------------------------------------------
    RETURN
  END SUBROUTINE mlocean_plasim_sst
  ! **********************************************************************
END MODULE messy_mlocean_plasim
! **********************************************************************
