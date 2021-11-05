!>
!! Contains parameters needed in several other modules.
!!
!! Values are initialised in program mpiom
!!
!! @ingroup common
!!
MODULE mo_commoau1

  USE mo_kind, ONLY: wp
  IMPLICIT NONE

  INTEGER :: isnflg

  !> melting enthalpy of ice = 320*10**6 ws/m**3
  REAL(wp), PARAMETER :: entmel = 320.e6_wp
  !> ocean ice/water interface temperature [degrees C]
  REAL(wp), PARAMETER :: tfreeze = -1.9_wp
  !> stefan boltzmann constant =5.67 10**(-8)w/m**2/k**4
  REAL(wp), PARAMETER :: stebol = 5.67e-8_wp
  ! hth00: arbitrary constant for discrimination between thin and thick ice
  ! set to 1.5 m,   not in use
  ! REAL(wp) :: hth00 = 0.5_wp

  REAL(wp) :: h0
  REAL(wp) :: armin
  REAL(wp) :: armax
  REAL(wp) :: hmin
  REAL(wp), PARAMETER :: cc = 4.2e6_wp
  REAL(wp) :: sice
  REAL(wp) :: sicthmin
  REAL(wp) :: hsntoice !< Maximal ice thickness for which snow to ice conversion is done
  !> Enthalpy of vaporization for water at 0°C [J/kg]
  REAL(wp), PARAMETER :: vapl = 2.5e6_wp
  !< Enthalpy of sublimation for water at 0°C [J/kg]
  !! (H_sublimation - H_vaporization = H_fusion ~= 334 kJ/kg)
  REAL(wp), PARAMETER :: subl = 2.834e6_wp

  !> mo_planetary_constants::rhoref_ice * lat. heat of fusion at ice surface [J/m³]
  !! Tabulated value would be rhoice * enthalpy of fusion =
  !! 910 kg/m³ * 334 kJ/kg =  303940 kJ/m³ ~= 3.04E8
  REAL(wp), PARAMETER :: clo = 3.02e8_wp
  !> mo_planetary_constants::rhoref_ice * lat. heat of fusion at ice bottom [J/m³]
  !! This is the reduced value for ice bottom (see Semtner76)
  !! This is now used for both surface and bottom of ice!
  REAL(wp), PARAMETER :: clb = 2.70e8_wp
  !> melting temperature of water ice [K]
  REAL(wp), PARAMETER :: tmelt = 273.16_wp
  REAL(wp), PARAMETER :: con = 2.1656_wp
  REAL(wp), PARAMETER :: consn = 0.31_wp

END MODULE mo_commoau1
