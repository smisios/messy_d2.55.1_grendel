!> module mo_planetary_constants
!> provides the model with parameters specific to the planet
!> for which the simulation runs
!> @author Thomas Jahns <jahns@dkrz.de>
!> @date 2009-04-29
MODULE mo_planetary_constants
  USE mo_kind, ONLY: wp
  USE mo_constants, ONLY: api
  IMPLICIT NONE
  !
  !>  GRAVITATIONAL ACCELERATION (CONSTANT 9.81 M/S**2)
  !     TO ACCOUNT FOR THE EXISTENCE OF MINI BLACK HOLES
  REAL(wp), PARAMETER :: G = 9.81_wp
  !
  !  ROCD   :  WIND STRESS DRAG COEFFICIENT
  !REAL(wp), PARAMETER :: ROCD = 1.2*1.5E-3
  !
  !>  RADIUS OF THE ABOVE MENTIONED PLANET
  REAL(wp), PARAMETER :: radius = 6371.E+3_wp
  !
  !> time per rotation in seconds, sidereal day
  REAL(wp), PARAMETER :: trot = 86164._wp
  !
  !>  ANGULAR VELOCITY OF OUR NICE BLUE PLANET
#ifdef MESSY
  REAL(wp), PARAMETER :: OMEGA = 2._wp * api/trot
#else
  REAL(wp), PARAMETER :: OMEGA = 2._wp * REAL(api, wp)/trot
#endif
  !>  HEAT CAPACITY PER CUBICMETER
  REAL(wp), PARAMETER :: ROCP = 4.E06_wp

  !>  OCEAN REFERENCE DENSITY [kg m-3]
  REAL(wp), PARAMETER :: rhoref_water = 1025.0_wp
  REAL(wp), PARAMETER :: inv_rhoref_water = 1._wp/rhoref_water

  !>  SEA LEVEL REFERENCE PRESSURE [Pa]
  REAL(wp), PARAMETER :: slpref = 101300._wp

  !>  SEA ICE REFERENCE DENSITY [kg m-3]
  REAL(wp), PARAMETER :: rhoref_ice = 910.0_wp

  !>  SNOW REFERENCE DENSITY [kg m-3]
  REAL(wp), PARAMETER :: rhoref_snow = 330.0_wp

  !>  AIR REFERENCE DENSITY [kg m-3]
  REAL(wp), PARAMETER :: rhoref_air = 1.3_wp

  !>  RATIO REFERENCE DENSITY ICE TO REFERENCE DENSITY WATER
  REAL(wp), PARAMETER :: rhoicwa = rhoref_ice/rhoref_water
  !>  RATIO REFERENCE DENSITY SNOW TO REFERENCE DENSITY WATER
  REAL(wp), PARAMETER :: rhosnwa = rhoref_snow/rhoref_water
  !>  RATIO REFERENCE DENSITY SNOW TO REFERENCE DENSITY ICE
  REAL(wp), PARAMETER :: rhosnic = rhoref_snow/rhoref_ice

  !> heat conductivity of ice = 2 w/m/k   cuwe not in use
  REAL(wp), PARAMETER :: sichec = 2._wp

END MODULE mo_planetary_constants
