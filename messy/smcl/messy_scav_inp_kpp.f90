MODULE MESSY_SCAV_INP_KPP

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
!!$  REAL(DP), PUBLIC, POINTER, DIMENSION(:)   :: JX
!!$
!!$  REAL(dp), PUBLIC                          :: TEMP, CV_l
!!$  REAL(dp), PUBLIC, POINTER, DIMENSION(:)   :: k_exf(:), k_exb(:)
!!$  REAL(DP), PUBLIC                          :: K_EXF_N2O5, K_EXF_CLNO3, K_EXF_BRNO3
!!$  REAL(DP), PUBLIC                          :: K_EXF_Hg, K_EXF_RGM
!!$
  INTEGER, PUBLIC, PARAMETER                :: IP_H2O2 = 1
  INTEGER, PUBLIC, PARAMETER                :: IP_HNO3 = 2
  INTEGER, PUBLIC, PARAMETER                :: IP_O3   = 3
  INTEGER, PUBLIC, PARAMETER                :: IP_MAX   = 3

!!$  ! VARIABLES WITH VARIABLE VECTOR LENGTH
!!$  REAL(DP), PUBLIC, POINTER, DIMENSION(:,:) :: VK_EXF, VK_EXB
!!$  REAL(DP), PUBLIC, POINTER, DIMENSION(:)   :: VCV_L, VK_EXF_N2O5,          &
!!$                                               VK_EXF_CLNO3, VK_EXF_BRNO3,  &
!!$                                               VK_EXF_Hg, VK_EXF_RGM
!!$  REAL(DP), PUBLIC, POINTER, DIMENSION(:)   :: VTEMP, VPRESS
!!$  REAL(DP), PUBLIC, POINTER, DIMENSION(:,:) :: JX_VEC

  REAL(DP), PUBLIC, PARAMETER :: &
                      F_EQUIL_HO2   = 1.E5_DP, F_EQUIL_HONO   = 1.E5_DP, &
                      F_EQUIL_HNO4  = 1.E5_DP, F_EQUIL_HCL    = 1.E2_DP, &
                      F_EQUIL_CL2M  = 1.E2_DP, F_EQUIL_HBR    = 1.E7_DP, &
                      F_EQUIL_BRCL2 = 1.E5_DP, F_EQUIL_BR2CL  = 1.E5_DP, &
                      F_EQUIL_HOCL  = 1.E2_DP, F_EQUIL_BR2M   = 1.E2_DP, &
                      F_EQUIL_HOBR  = 1.E2_DP, F_EQUIL_BR2CLM = 1.E5_DP, &
                      F_EQUIL_ICL   = 1.E2_DP, F_EQUIL_IBR    = 1.E2_DP, &
                      F_EQUIL_ICLBR = 1.E2_DP
  REAL(DP), PUBLIC, PARAMETER ::  &
                      F_EQUIL_H2O   = 1.0_DP,   F_EQUIL_NH3   = 1.E4_DP, &
                      F_EQUIL_HNO3  = 1.E2_DP,  F_EQUIL_CO2   = 1.E1_DP, &
                      F_EQUIL_HCOOH = 1.E1_DP,  F_EQUIL_SO2   = 1.E2_DP, &
                      F_EQUIL_HSO3M = 1.E2_DP,  F_EQUIL_HSO4M = 1.E2_DP, &
                      F_EQUIL_H2SO4 = 1.E1_DP,  F_EQUIL_CH3CO2H = 5.E-1_DP

CONTAINS


END MODULE MESSY_SCAV_INP_KPP
