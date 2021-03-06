MODULE mo_constants

!- Description:
!
!  This module contains basic constants and derived constants
!
!- Author:
!
!  M. Giorgetta, MPI, April 1999
!  I. Kirchner, MPI, December 2000, time control
!  L. Kornblueh, MPI, January 2001, cleanup

  USE mo_kind, ONLY: dp

! mz_pj_20040405+
#ifdef MESSY
  USE messy_main_constants_mem, ONLY: ak => k_B,      amd => M_air, &
                                      amw => M_H2O,   api => pi,    &
                                      argas => R_gas, avo => N_A,   &
                                      g, cpd=>cp_air, a=>radius_earth , &
                                      rd, rv, stbo, MC, MO, MN, M_O3, MH, &
                                      rhoh2o => rho_H2O
#endif
! mz_pj_20040405-
  IMPLICIT NONE

! Universal constants
#ifndef MESSY
  REAL(dp), PARAMETER :: api   = 3.14159265358979323846_dp ! pi
! REAL(dp), PARAMETER :: ar    = 8.314e3_dp       ! universal gas constant in J/K/kmol
  REAL(dp), PARAMETER :: argas = 8.314409_dp      ! universal gas constant in J/K/mol  
  REAL(dp), PARAMETER :: avo   = 6.022045e23_dp   ! Avogadro constant in 1/mol
  REAL(dp), PARAMETER :: ak    = 1.380662e-23_dp  ! Boltzmann constant in J/K
  REAL(dp), PARAMETER :: stbo  = 5.67E-8_dp  ! Stephan-Boltzmann constant [W/m2/K4]

! Molar weights in g/mol
  REAL(dp),PARAMETER :: amco2 = 44.011_dp   ! molecular weight of carbon dioxide
  REAL(dp),PARAMETER :: amch4 = 16.043_dp   ! molecular weight of methane
  REAL(dp),PARAMETER :: amo3  = 47.9982_dp  ! molecular weight of ozone
  REAL(dp),PARAMETER :: amn2o = 44.013_dp   ! molecular weight of N2O
#else
  REAL(dp),PARAMETER :: amco2 = MC+2._dp*MO ! molecular weight of carbon dioxide
  REAL(dp),PARAMETER :: amch4 = MC+4._dp*MH ! molecular weight of methane
  REAL(dp),PARAMETER :: amo3  = M_O3        ! molecular weight of ozone
  REAL(dp),PARAMETER :: amn2o = 2._dp*MN+MO ! molecular weight of N2O
#endif
!!$  REAL(dp),PARAMETER :: amc11 = 137.3686_dp ! molecular weight of CFC11
!!$  REAL(dp),PARAMETER :: amc12 = 120.9140_dp ! molecular weight of CFC12
#ifndef MESSY
  REAL(dp),PARAMETER :: amw   = 18.0154_dp  ! molecular weight of water vapor
  REAL(dp),PARAMETER :: amd   = 28.970_dp   ! molecular weight of dry air
#endif

! Dry air and water vapour thermodynamic constants
#ifndef MESSY
  REAL(dp),PARAMETER :: cpd   = 1005.46_dp  ! specific heat of dry air at
                                            ! constant pressure in J/K/kg
#endif
  REAL(dp),PARAMETER :: cpv   = 1869.46_dp  ! specific heat of water vapour at 
                                            ! constant pressure in J/K/kg
#ifndef MESSY
  REAL(dp),PARAMETER :: rd    = 287.05_dp   ! gas constant for dry air in J/K/kg
  REAL(dp),PARAMETER :: rv    = 461.51_dp   ! gas constant for water vapour 
                                            ! in J/K/kg
#endif
  REAL(dp),PARAMETER :: rcpd  = 1.0_dp/cpd  ! auxiliary constant in K*kg/J
  REAL(dp)           :: vtmpc1= rv/rd-1.0_dp! dimensionless auxiliary constant
  REAL(dp)           :: vtmpc2= cpv/cpd-1.0_dp ! dimensionless auxiliary constant

! H2O related constants, liquid density, phase change constants
#ifndef MESSY
  REAL(dp), PARAMETER :: rhoh2o= 1000.0_dp        ! density of liquid water in kg/m3
#endif
  REAL(dp), PARAMETER :: alv   = 2.5008e6_dp      ! latent heat for vaporisation in J/kg
  REAL(dp), PARAMETER :: als   = 2.8345e6_dp      ! latent heat for sublimation in J/kg
  REAL(dp), PARAMETER :: alf   = als-alv          ! latent heat for fusion in J/kg
  REAL(dp), PARAMETER :: clw   = 4186.84_dp       ! specific heat for liquid waterJ/K/kg
  REAL(dp), PARAMETER :: tmelt = 273.15_dp        ! melting temperature of ice/snow

! Earth and earth orbit parameters
#ifndef MESSY
  REAL(dp),PARAMETER :: a     = 6371000.0_dp! radius of the earth in m
#endif
  REAL(dp),PARAMETER :: omega = 0.7292E-4_dp! solid rotation velocity of the
                                            ! earth in 1/s
#ifndef MESSY
  REAL(dp),PARAMETER :: g     = 9.80665_dp  ! gravity acceleration in m/s2
#endif

! Constants used for computation of saturation mixing ratio
! over liquid water (*c_les*) or ice(*c_ies*)
  REAL(dp), PARAMETER :: c1es  = 610.78_dp           ! 
  REAL(dp), PARAMETER :: c2es  = c1es*rd/rv          ! 
  REAL(dp), PARAMETER :: c3les = 17.269_dp           ! 
  REAL(dp), PARAMETER :: c3ies = 21.875_dp           ! 
  REAL(dp), PARAMETER :: c4les = 35.86_dp            ! 
  REAL(dp), PARAMETER :: c4ies = 7.66_dp             ! 
  REAL(dp), PARAMETER :: c5les = c3les*(tmelt-c4les) ! 
  REAL(dp), PARAMETER :: c5ies = c3ies*(tmelt-c4ies) ! 
  REAL(dp), PARAMETER :: c5alvcp = c5les*alv/cpd     ! 
  REAL(dp), PARAMETER :: c5alscp = c5ies*als/cpd     ! 
  REAL(dp), PARAMETER :: alvdcp  = alv/cpd           ! 
  REAL(dp), PARAMETER :: alsdcp  = als/cpd           ! 

! IK
! for later use, ice and snow constants
!  REAL(dp),parameter :: ice_dmin = 0.10_dp     ! minimal ice thicknessss
!  REAL(dp),parameter :: ice_rho  = 910.0_dp    ! ice density
!  REAL(dp),parameter :: ice_cp   = 2106.0_dp   ! specific ice heat capacity ??
!  REAL(dp),parameter :: ice_alpha= 2.1656_dp   ! ??
!  REAL(dp),parameter :: ice_rici = 2.09e+06_dp ! volumetric heat capacity of
!                                               ! ice [j/m**3/k]
!  REAL(dp),parameter :: ice_difiz= 12.e-07_dp  ! temperature diffusivity of ice
!                                               !  [m**2/s]
!  REAL(dp),parameter :: sn_cond  = 0.22_dp     ! snow thermal conductivity
!                                               ! [j/s/m/k]
!  REAL(dp),parameter :: sn_dens  = 300.0_dp    ! snow density   [kg/m**3]
!  REAL(dp),parameter :: sn_capa  = 634500.0_dp ! snow vol. heat capacity
!                                               !  [j/m**3/k]
  
END MODULE mo_constants
