MODULE MESSY_GMXE_AERCHEM_INP_KPP

! Author: Holger Tost, Mainz, Sept 2019
!
! This module contains the parameters required to be used by KPP
! 
! To avoid circular dependencies these variables are defined here, 
! but are filled with values in the module MESSY_SCAV_LIQ

  USE MESSY_MAIN_CONSTANTS_MEM,    ONLY: dp
  USE MESSY_MAIN_TOOLS_KINETICS

  IMPLICIT NONE
  !!!PRIVATE
  SAVE
!  PUBLIC :: K_ARR, K_SIV_H2O2
  
  PRIVATE                                   :: DP

  INTEGER, PUBLIC, PARAMETER                :: IP_H2O2 = 1
  INTEGER, PUBLIC, PARAMETER                :: IP_HNO3 = 2
  INTEGER, PUBLIC, PARAMETER                :: IP_O3   = 3
  INTEGER, PUBLIC, PARAMETER                :: IP_MAX  = 3

  REAL(DP), PUBLIC, PARAMETER :: &
                      F_EQUIL_HO2   = 1.E5_DP, F_EQUIL_HONO   = 1.E5_DP, &
                      F_EQUIL_HNO4  = 1.E5_DP, F_EQUIL_HCL    = 1.E2_DP, &
                      F_EQUIL_CL2M  = 1.E2_DP, F_EQUIL_HBR    = 1.E0_DP, &
                      F_EQUIL_BRCL2 = 1.E0_DP, F_EQUIL_BR2CL  = 1.E0_DP, &
                      F_EQUIL_HOCL  = 1.E0_DP, F_EQUIL_BR2M   = 1.E0_DP, &
                      F_EQUIL_HOBR  = 1.E0_DP, F_EQUIL_BR2CLM = 1.E0_DP, &
                      F_EQUIL_ICL   = 1.E2_DP, F_EQUIL_IBR    = 1.E2_DP, &
                      F_EQUIL_ICLBR = 1.E2_DP
  REAL(DP), PUBLIC, PARAMETER ::  &
                      F_EQUIL_H2O   = 1.0_DP,   F_EQUIL_NH3   = 1.E4_DP, &
                      F_EQUIL_HNO3  = 1.E2_DP,  F_EQUIL_CO2   = 1.E1_DP, &
                      F_EQUIL_HCOOH = 1.E1_DP,  F_EQUIL_SO2   = 1.E2_DP, &
                      F_EQUIL_HSO3M = 1.E2_DP,  F_EQUIL_HSO4M = 1.E2_DP, &
                      F_EQUIL_H2SO4 = 1.E1_DP,  F_EQUIL_CH3CO2H = 5.E-1_DP

  REAL(dp), PARAMETER :: thres_lwc = 5.e-8_dp ! minimum threshold lwc for aerchem calculations
CONTAINS


  ELEMENTAL REAL(DP) FUNCTION K_olig_glyx(c_H2O, jx_h2o2)
    ! SPECIAL RATE FUNCTION FOR oligomerisation of glyx in aerosol water
    ! according to Ervens et al., 4 s^-1
    ! scaling with photolysis of H2O2 in liquid phase (-> factor 2.33
    ! upscaling by 1e3 due to low yield
    ! scald with liquid water content -> more efficient at the lowest lwc
    USE MESSY_MAIN_CONSTANTS_MEM,      ONLY: AVO => N_A
    REAL(dp), INTENT(IN) :: c_h2o
    REAL(dp), INTENT(IN) :: jx_h2o2

    K_OLIG_GLYX = 4._dp * JX_H2O2 * 2.33_dp * 1.e3_dp * &
                  thres_lwc/ (C_H2O / (55.5_DP * AVO * 1.e-6_dp))


  END FUNCTION K_OLIG_GLYX

END MODULE MESSY_GMXE_AERCHEM_INP_KPP
