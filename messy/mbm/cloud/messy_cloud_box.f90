PROGRAM MESSY_CLOUD_BOX

  USE MESSY_MAIN_CONSTANTS_MEM,   ONLY: dp
  USE MESSY_MAIN_TOOLS,           ONLY: init_convect_tables
  USE MESSY_CLOUD_DROPLET


  IMPLICIT NONE

  INTEGER  :: kproma, klev
  INTEGER, PARAMETER :: nmod  = 1
!  INTEGER, PARAMETER :: nmod  = 4

!  INTEGER, PARAMETER :: nspec = 1
!  INTEGER, PARAMETER :: nspec = 2
  INTEGER, PARAMETER :: nspec = 3

  INTEGER  :: sss = 1 ! op_pj_20160815 sup_sat_scheme
  REAL(dp) :: kappa(nspec) ! kappa parameter of species

  REAL(DP) :: TEMP, PRESS, velo
  REAL(dp) :: sigma(nmod), aer_rad(nmod), aer_density(nspec)
  REAL(dp) :: mass(nmod,nspec)
  REAL(dp) :: number(nmod)

  REAL(dp) :: S_max, nfrac(nmod), aer_crit(nmod), s_crit(nmod)
  REAL(dp) :: cmr_to_mmr(nmod), kfrac(nmod), mfrac(nmod)
  REAL(dp) :: N_activ

  REAL(dp) :: nu(nspec)           ! ion number of a species
  REAL(dp) :: eps(nspec)          ! mass fraction of water soluble substance
  REAL(dp) :: Phi(nspec)          ! osmotic coefficient of species
  REAL(dp) :: M_ap(nspec)         ! molecular weight of the aerosol species

  INTEGER  :: ji

  temp  = 283.15_dp
  press = 80000._dp
  velo  = 5._dp

  sss = 1                         ! op_pj_20160815
  kappa = (/ 0.7, 1.2, 0.001 /)   ! mz_ht_20161010

  nu(1)   = 3._dp   ! (NH_4)_2SO_4 decomposes into three ions
  eps(:)  = 1._dp   ! completely soluble
  Phi(:)  = 0.7_dp
  M_ap(1) = 132._dp 

  nu(2)   = 2._dp   ! NaCl decomposes into two ions
  M_ap(2) = 58.44_dp 

  nu(3)   = 1._dp   ! BC does not decompose into ions
  M_ap(3) = 12._dp 
  eps(3)  = 1.e-30_dp

  sigma(1) = 2.5_dp
!!$  sigma(1)   = 1.5_dp
!!$  sigma(2)   = 1.5_dp
!!$  sigma(3)   = 1.5_dp
!!$  sigma(4)   = 2.5_dp

  cmr_to_mmr(1) = EXP(3.5*(LOG(sigma(1)))**2)
!!$  cmr_to_mmr(2) = EXP(3.5*(LOG(sigma(2)))**2)
!!$  cmr_to_mmr(3) = EXP(3.5*(LOG(sigma(3)))**2)
!!$  cmr_to_mmr(4) = EXP(3.5*(LOG(sigma(4)))**2)!

  aer_rad(1)  = 0.01_dp * 1.e-6_dp
!!$  aer_rad(1)  = 0.001_dp * 1.e-6_dp
!!$  aer_rad(2)  = 0.01_dp * 1.e-6_dp
!!$  aer_rad(3)  = 0.1_dp * 1.e-6_dp
!!$  aer_rad(4)  = 1.0_dp * 1.e-6_dp

  number(1) = 200._dp * 1.e6_dp ! 200 pro cm^3 
!!$  number(1) = 1000._dp * 1.e6_dp ! 1000 pro cm^3 
!!$  number(2) = 500._dp * 1.e6_dp  ! 500 pro cm^3 
!!$  number(3) = 200._dp * 1.e6_dp  ! 300 pro cm^3 
!!$  number(4) = 50._dp * 1.e6_dp  ! 200 pro cm^3 


  mass(1,1) = 1.e-11_dp
  mass(1,2) = 2.e-11_dp
  mass(1,3) = 1.e-12_dp

!!$  mass(2,1) = 5.e-10_dp
!!$  mass(2,2) = 8.e-10_dp

!!$  mass(3,1) = 1.e-10_dp
!!$  mass(3,2) = 2.e-10_dp

!!$  mass(4,1) = 1.e-9_dp
!!$  mass(4,2) = 2.e-9_dp

  aer_density(1) = 1770._dp   ![kg/m^3]  ! (NH4)2SO4
  aer_density(2) = 2156._dp   ![kg/m^3]  ! NaCl
  aer_density(3) = 2000._dp   ![kg/m^3]  ! BC

  CALL init_convect_tables

  do ji=1,10000
    aer_rad = 0.0001_dp * real(ji,dp) * 1.e-6_dp ! converted from µm to m
!    velo = 0.01_dp * real(ji,dp) 
    CALL cloud_droplet_ARG(temp, press, velo,                              &
                           aer_rad, aer_density, sigma, number, mass,      &
! op_pj_20160815+
!!$                        nu, Phi, eps, M_ap,                             &
!!$                        s_crit, nfrac, S_max, aer_crit, N_activ,        &
                        nu, Phi, eps, kappa, M_ap,                         &
                        s_crit, sss, nfrac, S_max, aer_crit, N_activ,        &
! op_pj_20160815-
                           kproma, klev, nmod, nspec,                      &
                           kfrac, mfrac, cmr_to_mmr)
  print*, aer_rad, nfrac, s_max, aer_crit, N_activ
!  print*, velo, nfrac, s_max, aer_crit, N_activ
  enddo




END PROGRAM MESSY_CLOUD_BOX
