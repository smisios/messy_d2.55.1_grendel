MODULE MESSY_GMXE_AERCHEM_LIQ


  USE MESSY_GMXE_AERCHEM_KPP,    ONLY: lspec => nspec, DP

  USE messy_main_tools,          ONLY: PTR_3D_ARRAY, PTR_2D_ARRAY

  IMPLICIT NONE
  REAL(DP), PUBLIC               :: TEMP, CV_L
  REAL(DP), PUBLIC               :: K_EXF(0:LSPEC), K_EXB(0:LSPEC)
  REAL(DP), PUBLIC               :: K_EXF_N2O5, K_EXF_CLNO3, K_EXF_BRNO3
  REAL(DP), PUBLIC               :: K_EXF_Hg, K_EXF_RGM
  REAL(dp), PUBLIC               :: HS0(0:LSPEC), DHT(0:LSPEC)
  REAL(dp), PUBLIC               :: alpha0(0:LSPEC), alpha_T(0:LSPEC)

  INTEGER, PUBLIC, PARAMETER     :: IP_H2O2 = 1
  INTEGER, PUBLIC, PARAMETER     :: IP_NO3  = 2
  INTEGER, PUBLIC, PARAMETER     :: IP_O3   = 3

! VARIABLES WITH VARIABLE VECTOR LENGTH
  REAL(DP), ALLOCATABLE, DIMENSION(:,:), PUBLIC ::       &
            VK_EXF, VK_EXB
  REAL(DP), ALLOCATABLE, DIMENSION(:), PUBLIC ::         &
            VCV_L, VK_EXF_N2O5, VK_EXF_CLNO3, VK_EXF_BRNO3, &
            VK_EXF_Hg, VK_EXF_RGM
! PHOTOLYSIS_VALUES
  REAL(DP), PUBLIC, DIMENSION(:,:,:), POINTER   :: JVAL_H2O2
  REAL(DP), PUBLIC, DIMENSION(:,:), POINTER     :: JVAL_H2O2_2d
  REAL(DP), PUBLIC, DIMENSION(:,:,:), POINTER   :: JVAL_NO3
  REAL(DP), PUBLIC, DIMENSION(:,:), POINTER     :: JVAL_NO3_2d
  REAL(DP), PUBLIC, DIMENSION(:,:,:), POINTER   :: JVAL_O3
  REAL(DP), PUBLIC, DIMENSION(:,:), POINTER     :: JVAL_O3_2d
  REAL(DP), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: JX_VEC
 ! TEMPERATURE , PRESSURE, DENSITY
  REAL(DP), PUBLIC, ALLOCATABLE, DIMENSION(:) :: VTEMP, VPRESS, VRHO, VRDRAD
  REAL(DP), PUBLIC, ALLOCATABLE, DIMENSION(:) :: VLIQ_LWC
  REAL(DP), ALLOCATABLE, DIMENSION(:,:)       :: ALPHA_V, VMEAN_V
  REAL(dp), PARAMETER :: scal_fac=1e-2_dp
!  REAL(DP), PUBLIC, PARAMETER :: &
!                      F_EQUIL_HO2   = 1.E5_DP, F_EQUIL_HONO   = 1.E5_DP, &
!                      F_EQUIL_HNO4  = 1.E5_DP, F_EQUIL_HCL    = 1.E2_DP, &
!                      F_EQUIL_CL2M  = 1.E2_DP, F_EQUIL_HBR    = 1.E7_DP, &
!                      F_EQUIL_BRCL2 = 1.E5_DP, F_EQUIL_BR2CL  = 1.E5_DP, &
!                      F_EQUIL_HOCL  = 1.E2_DP, F_EQUIL_BR2M   = 1.E2_DP, &
!                      F_EQUIL_HOBR  = 1.E2_DP, F_EQUIL_BR2CLM = 1.E5_DP, &
!                      F_EQUIL_ICL   = 1.E2_DP, F_EQUIL_IBR    = 1.E2_DP, &
!                      F_EQUIL_ICLBR = 1.E2_DP
!  REAL(DP), PUBLIC, PARAMETER ::  &
!!                      F_EQUIL_H2O   = 1.0_DP,   F_EQUIL_NH3   = 1.E4_DP, &
!                      F_EQUIL_H2O   = 1.0e-3_DP,   F_EQUIL_NH3   = 1.E4_DP, &
!                      F_EQUIL_HNO3  = 1.E2_DP,  F_EQUIL_CO2   = 1.E1_DP, &
!                      F_EQUIL_HCOOH = 1.E1_DP,  F_EQUIL_SO2   = 1.E2_DP, &
!                      F_EQUIL_HSO3M = 1.E2_DP,  F_EQUIL_HSO4M = 1.E2_DP, &
!                      F_EQUIL_H2SO4 = 1.E-3_DP,  F_EQUIL_CH3CO2H = 5.E-1_DP
!             F_EQUIL_H2SO4 = 1.E1_DP,  F_EQUIL_CH3CO2H = 5.E-1_DP
 REAL(DP), PUBLIC, PARAMETER :: &
                      F_EQUIL_HO2   = scal_fac * 1.E5_DP, F_EQUIL_HONO   = scal_fac * 1.E5_DP, &
                      F_EQUIL_HNO4  = scal_fac * 1.E5_DP, F_EQUIL_HCL    = scal_fac * 1.E2_DP, &
                      F_EQUIL_CL2M  = scal_fac * 1.E2_DP, F_EQUIL_HBR    = scal_fac * 1.E7_DP, &
                      F_EQUIL_BRCL2 = scal_fac * 1.E5_DP, F_EQUIL_BR2CL  = scal_fac * 1.E5_DP, &
                      F_EQUIL_HOCL  = scal_fac * 1.E2_DP, F_EQUIL_BR2M   = scal_fac * 1.E2_DP, &
                      F_EQUIL_HOBR  = scal_fac * 1.E2_DP, F_EQUIL_BR2CLM = scal_fac * 1.E5_DP, &
                      F_EQUIL_ICL   = scal_fac * 1.E2_DP, F_EQUIL_IBR    = scal_fac * 1.E2_DP, &
                      F_EQUIL_ICLBR = scal_fac * 1.E2_DP
  REAL(DP), PUBLIC, PARAMETER ::  &
!                      F_EQUIL_H2O   = scal_fac * 1.0_DP,   F_EQUIL_NH3   = scal_fac * 1.E4_DP, &
                      F_EQUIL_H2O   = scal_fac * 1.0e-3_DP,   F_EQUIL_NH3   = scal_fac * 1.E4_DP, &
                      F_EQUIL_HNO3  = scal_fac * 1.E2_DP,  F_EQUIL_CO2   = scal_fac * 1.E1_DP, &
                      F_EQUIL_HCOOH = scal_fac * 1.E1_DP,  F_EQUIL_SO2   = scal_fac * 1.E2_DP, &
                      F_EQUIL_HSO3M = scal_fac * 1.E2_DP,  F_EQUIL_HSO4M = scal_fac * 1.E2_DP, &
                      F_EQUIL_H2SO4 = scal_fac * 1.E-3_DP,  F_EQUIL_CH3CO2H = scal_fac * 5.E-1_DP
!             F_EQUIL_H2SO4 = scal_fac * 1.E1_DP,  F_EQUIL_CH3CO2H = scal_fac * 5.E-1_DP
  REAL(dp), PUBLIC, DIMENSION(:,:), ALLOCATABLE :: L_SPEC
  INTEGER :: nspec_gas, nspec_aer

  INTEGER, PARAMETER :: NMAX_SPEC=999
 
  INTEGER :: NVEC

  TYPE IDX_L_2D_PTR_ARRAY
    INTEGER, DIMENSION(:,:), POINTER ::  GAS_SPEC => NULL()
    INTEGER, DIMENSION(:,:), POINTER ::  LIQ_SPEC => NULL()
    REAL(DP),DIMENSION(:,:), POINTER ::  GAS_ATTR => NULL()
    REAL(DP),DIMENSION(:,:), POINTER ::  LIQ_ATTR => NULL()
  END TYPE IDX_L_2D_PTR_ARRAY
  TYPE(IDX_L_2D_PTR_ARRAY), SAVE :: KPP_L_IDX 
  INTEGER :: LSPEC_GAS, LSPEC_LIQ
! FOR GAS_SPEC
  INTEGER, PARAMETER  :: GAS_IDX  = 1   ! index of gas phase compounds 
                                        !     within the species array
  INTEGER, PARAMETER  :: GAS2TRAC = 2   ! index of tracer for gas phase
                                        !     compound of gas species
  INTEGER, PARAMETER  :: GAS_G2P  = 3   ! index of corresponding species
                                        !     in gas phase and liquid / ice
                                        !     phase
  INTEGER, PARAMETER  :: GAS_MAX  = 3
! FOR GAS_ATTR
  INTEGER, PARAMETER  :: GAS_MW   = 1   ! molar weight / 100.
!  INTEGER, PARAMETER  :: GAS_HS0  = 2   ! henry coefficient for gas species
!  INTEGER, PARAMETER  :: GAS_DHT  = 3   ! henry coefficient temperature 
                                        ! dependence for gas species
!  INTEGER, PARAMETER  :: GAS_ALPHA= 4   ! alpha coefficient for gas species
  INTEGER, PARAMETER  :: GATT_MAX = 1
! FOR LIQ_SPEC  
  INTEGER, PARAMETER  :: LIQ_IDX  = 1   ! index of liquid phase compounds
                                        !     within the species array
  INTEGER, PARAMETER  :: LIQ2TRAC = 2   ! index of tracer for liquid phase
                                        !     compound of aerosol species
  INTEGER, PARAMETER  :: LIQ_MAX  = 2
! FOR LIQ_ATTR
  INTEGER, PARAMETER  :: LIQ_MW   = 1   ! index of molar weight
  INTEGER, PARAMETER  :: LIQ_DENS = 2   ! index of liquid density
  INTEGER, PARAMETER  :: LATT_MAX = 2

  INTEGER   :: NMAXKPP = 10000   ! MAXIMUM NUMBER OF KPP CALCULATED TIMESTEPS

  INTEGER   :: nmod_aerchem = 1

  
  TYPE PROD_AERCHEM
    REAL(dp), DIMENSION(:,:,:),        POINTER :: ptr_3d
    REAL(dp), DIMENSION(:,:),          POINTER :: ptr_2d
  END type PROD_AERCHEM
  TYPE DIAGAERCHEM
     TYPE(PROD_AERCHEM), DIMENSION(:), POINTER :: PROD
     INTEGER, DIMENSION(:),            POINTER :: PROD_IDX
     CHARACTER(LEN=2), DIMENSION(:),   POINTER :: PROD_CHAR
     REAL(dp),DIMENSION(:,:,:),        POINTER :: PH
     REAL(dp),DIMENSION(:,:),          POINTER :: PH_2D
     REAL(dp),DIMENSION(:,:,:),        POINTER :: LWC
     REAL(dp),DIMENSION(:,:),          POINTER :: LWC_2D
  END TYPE DIAGAERCHEM
  INTEGER                                   :: nprod = 0
  TYPE (DIAGAERCHEM), DIMENSION(:), POINTER :: diag_aerchem

  INTRINSIC :: EXP, TINY
CONTAINS

!------------------------------------------------------------------------------

   SUBROUTINE TRANSFER_COEFFICIENT_LIQ

    USE MESSY_MAIN_CONSTANTS_MEM, ONLY: R_GAS
    USE MESSY_GMXE_AERCHEM_KPP

    IMPLICIT NONE

    INTEGER  :: JT, K, JL, IDX

    REAL(DP) :: ZDGAIR
    REAL(DP) :: ZRDRAD(NVEC)
    REAL(DP) :: ZTX(NVEC), KMT(NVEC)
    REAL(DP) :: ZHETT
    REAL(DP) :: ZLAMBDA(NVEC)
    REAL(DP) :: HENRY(LSPEC,NVEC)
 
   !CALCULATES THE TRANSFER COEFFICIENTS IN AND OUT OF THE WATER DROPLET
    VK_EXF(:,:)     = 0.
    VK_EXB(:,:)     = 0.
    VK_EXF_N2O5(:)  = 0.
    VK_EXF_CLNO3(:) = 0.
    VK_EXF_BRNO3(:) = 0.
    VK_EXF_Hg(:)    = 0.
    VK_EXF_RGM(:)   = 0.
    
    KMT(:) = 0._DP

    
    DO JL=1,NVEC
      ZLAMBDA(JL) = 2.28E-5_DP * VTEMP(jl) / VPRESS(jl)
      ZRDRAD(JL)  = VRDRAD(jl)
          ! RECALCULATION RADIUS FROM MM TO M
    ENDDO
    CALL AERCHEM_ALPHA_V
    CALL AERCHEM_VMEAN_V

    ZTX(1:NVEC) = 1._DP/VTEMP(1:NVEC) - 1._DP/298._DP
    HENRY(:,:) = 0.

    DO JT=1,LSPEC_GAS
      IDX = kpp_L_idx%gas_spec(jt,gas_idx)
      ! CALCULATE HENRY CONSTANT
      IF (HS0(IDX) > 0._DP ) THEN
        DO JL=1,NVEC
          HENRY(IDX,JL) = HS0(IDX) * EXP( DHT(IDX) * ZTX(JL) )
      ! UNIT: MOL/(L*ATM) --> MOL(AQ)/M3(AQ) / MOL(G)/M3(G) (I.E. DIMENSIONLESS)
      ! FCT=1.E3*8.3145*T/P_0=0.082*T
      ! PLUS CONVERSION TO INVERSE HENRY CONSTANT:  H_(H,INV)^CC = 1/K_H^CC
      ! I.E. MOL(AQ)/M3(AQ) / MOL(G)/M3(AIR) ->  MOL(G)/M3(AIR) / MOL(AQ)/M3(AQ)
          HENRY(IDX,JL) = 1._DP / (HENRY(IDX,JL) * 0.082 * VTEMP(JL))
        ENDDO

        ! HT     ALTERNATIVE FOLLOWING SCHWARTZ           
        DO JL=1,NVEC
          KMT(JL) = VMEAN_V(JL,IDX) /  (ZRDRAD(JL) * (ZRDRAD(JL) / &
                    ZLAMBDA(JL) + 4._DP / (3._DP * ALPHA_V(JL,IDX) )))
        ENDDO
         
        !  SETTING OF EXCHANGE FUNCTION            
        DO JL=1,NVEC
          VK_EXF(JL,IDX) = KMT(JL) * VLIQ_LWC(JL) * 1.E-3_DP 
          VK_EXB(JL,IDX) = KMT(JL) * HENRY(IDX,JL)
        ENDDO

      ENDIF
    ENDDO

    ! THIS USED TO BE FHET:
    DO JL = 1,NVEC
      ZHETT = L_SPEC(IND_H2O_L,JL)
      IF (IND_CLM_L/=0) ZHETT = ZHETT + 5.E2 * L_SPEC(IND_CLM_L,JL)
      IF (IND_BRM_L/=0) ZHETT = ZHETT + 3.E5 * L_SPEC(IND_BRM_L,JL)
      VK_EXF_N2O5(JL)  = VK_EXF(JL,IND_N2O5)  / ZHETT
      VK_EXF_CLNO3(JL) = VK_EXF(JL,IND_CLNO3) / ZHETT
      VK_EXF_BRNO3(JL) = VK_EXF(JL,IND_BRNO3) / ZHETT
      VK_EXF_Hg(JL)    = VK_EXF(JL,IND_Hg)    / ZHETT
      VK_EXF_RGM(JL)   = VK_EXF(JL,IND_HgO)   / ZHETT
    ENDDO
        ! END OF OLD FHET
    VK_EXF(1:NVEC,0) = 0._DP
    VK_EXB(1:NVEC,0) = 0._DP


  END SUBROUTINE TRANSFER_COEFFICIENT_LIQ
  
!-------------------------------------------------------------------------------

  SUBROUTINE EQUIL_DISSOCIATION_COEFF

    USE MESSY_MAIN_CONSTANTS_MEM,    ONLY: AVO => N_A

    REAL(DP) :: FAC
    INTEGER  :: JL
    
    VCV_L(:)      = 0.0_DP
    DO JL=1,NVEC
      IF (VLIQ_LWC(JL).GT.TINY(0._DP)) THEN
        VCV_L(JL) = 1.E6_DP / VLIQ_LWC(JL) / AVO
      ELSE
        VCV_L(JL) = 0.E0
      ENDIF
    ENDDO

  END SUBROUTINE EQUIL_DISSOCIATION_COEFF
!--------------------------------------------------------------------------
  SUBROUTINE LIQ_PHYSC_INIT(LON_LWC, LON_RDRAD,&
                            LON_TEMP, LON_PRESS)

! INITIALIZATION OF THE BOX VALUES FOR TEMPERATURE, PRESSURE, LIQUID WATER, 
! CHEMICAL CONCENTRATION OF WATER, DROPLET DIAMETER
    USE MESSY_MAIN_CONSTANTS_MEM,      ONLY: AVO => N_A, M_AIR, R_GAS
    USE MESSY_GMXE_AERCHEM_KPP,        ONLY: IND_H2O_L

    INTRINSIC :: REAL
    REAL(DP), INTENT(IN)  :: LON_LWC(NVEC), LON_RDRAD(NVEC),&
                             LON_TEMP(NVEC), LON_PRESS(NVEC)

    INTEGER :: JL

    DO JL=1,NVEC
      VLIQ_LWC(JL)     = LON_LWC(JL)
      IF (VLIQ_LWC(JL)*REAL(IND_H2O_L,DP) > 0._DP) &
        L_SPEC(IND_H2O_L,JL) = 55.5_DP * AVO * VLIQ_LWC(JL) * 1.E-6_DP
      VRDRAD(JL)   = LON_RDRAD(JL)
      VTEMP(JL)    = LON_TEMP(JL)
      VPRESS(JL)   = LON_PRESS(JL)
      VRHO(JL)     = VPRESS(JL) * M_AIR / (R_GAS * VTEMP(JL) * 1.E3_DP)
    ENDDO

  END SUBROUTINE LIQ_PHYSC_INIT

!------------------------------------------------------------------------------
  SUBROUTINE ALLOC_CHEM_FIELDS(LPROMA)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: LPROMA

    NVEC  = LPROMA

    ALLOCATE(VTEMP(NVEC))
    ALLOCATE(VPRESS(NVEC))
    ALLOCATE(VRHO(NVEC))
    ALLOCATE(VLIQ_LWC(NVEC))
    ALLOCATE(VRDRAD(NVEC))
    ALLOCATE(JX_VEC(NVEC,3))
    ALLOCATE(VCV_L(NVEC))
   
    ALLOCATE(VMEAN_V(NVEC,0:LSPEC))
    ALLOCATE(ALPHA_V(NVEC,0:LSPEC))

    ALLOCATE(VK_EXF(NVEC,0:LSPEC))
    ALLOCATE(VK_EXB(NVEC,0:LSPEC))
    ALLOCATE(VK_EXF_N2O5(NVEC))
    ALLOCATE(VK_EXF_CLNO3(NVEC))
    ALLOCATE(VK_EXF_BRNO3(NVEC))
    ALLOCATE(VK_EXF_Hg(NVEC))
    ALLOCATE(VK_EXF_RGM(NVEC))

  END SUBROUTINE ALLOC_CHEM_FIELDS
!-------------------------------------------------------------------
  SUBROUTINE DEALLOC_CHEM_FIELDS

    IMPLICIT NONE

    DEALLOCATE(VTEMP)
    DEALLOCATE(VPRESS)
    DEALLOCATE(VRHO)
    DEALLOCATE(VLIQ_LWC)
    DEALLOCATE(VRDRAD)
    DEALLOCATE(JX_VEC)
    DEALLOCATE(VCV_L)
    
    DEALLOCATE(VMEAN_V)
    DEALLOCATE(ALPHA_V)

    DEALLOCATE(VK_EXF)
    DEALLOCATE(VK_EXB)
    DEALLOCATE(VK_EXF_N2O5)
    DEALLOCATE(VK_EXF_CLNO3)
    DEALLOCATE(VK_EXF_BRNO3)
    DEALLOCATE(VK_EXF_Hg)
    DEALLOCATE(VK_EXF_RGM)

  END SUBROUTINE DEALLOC_CHEM_FIELDS

!------------------------------------------------------------------------------

  SUBROUTINE  AERCHEM_ALPHA_V

    USE messy_main_constants_mem,   ONLY: T0

    INTEGER  :: JT, JL, IDX
    REAL(dp) :: alpha_T0, alpha_Tdep

    DO jt=1,lspec_gas
      IDX = KPP_L_IDX%GAS_SPEC(JT,GAS_IDX)
      alpha_T0    = alpha0(IDX)
      alpha_Tdep  = alpha_T(IDX)
      DO jl=1,nvec
        alpha_v(jl,idx) = 1._dp / (1._dp+(1._dp/alpha_T0 - 1._dp) * &
                        EXP(alpha_Tdep * ((1._dp/T0)-(1._dp/vtemp(jl)))))
      ENDDO
    ENDDO

  END SUBROUTINE AERCHEM_ALPHA_V
!------------------------------------------------------------
  SUBROUTINE AERCHEM_VMEAN_V

    INTEGER  :: JL, JT, IDX
    REAL(dp) :: MW

    VMEAN_V(:,:) = -999.999_DP ! DUMMY VALUE

    DO JT=1,lspec_gas
      IDX = KPP_L_IDX%GAS_SPEC(JT,GAS_IDX)
      MW  = KPP_L_IDX%GAS_ATTR(JT,GAS_MW)*1.e3_dp
      DO JL=1,NVEC
        VMEAN_V(JL,IDX) = VMEAN_F(VTEMP(JL),MW)
      ENDDO
    ENDDO
  END SUBROUTINE AERCHEM_VMEAN_V
!---------------------------------------------------  
  ELEMENTAL REAL(DP) FUNCTION VMEAN_F(ZT,MOLWEIGHT)

      ! MEAN MOLECULAR SPEED FROM MAXWELL-BOLTZMANN DISTRIBUTION:
      ! VMEAN=SQRT(8*R_GAS*T/(M*PI))      (M IN KG/MOL)
      ! VMEAN IN M/S
      ! SQRT(8*R_GAS/PI)=4.60138

      REAL(DP), INTENT(IN) :: MOLWEIGHT
      REAL(DP), INTENT(IN) :: ZT

      INTRINSIC SQRT

      VMEAN_F = SQRT(ZT/MOLWEIGHT)*4.60138_dp

    END FUNCTION VMEAN_F

!!!FUNCTIONS FOR KPP
!=============================================================


!-----------------------------------------------------------------------------
  SUBROUTINE SET_PHOTO(IDX,JK)

    INTEGER, INTENT(IN) :: IDX(NVEC),JK
    INTEGER :: JL

    INTRINSIC :: ASSOCIATED

    JX_VEC(1:NVEC,:) = 0._DP

    IF (ASSOCIATED(JVAL_H2O2)) THEN
      DO JL=1,NVEC
        JX_VEC(JL,IP_H2O2) = JVAL_H2O2_2d(IDX(JL),JK)
      END DO
    ENDIF
    IF (ASSOCIATED(JVAL_NO3)) THEN
      DO JL=1,NVEC
        JX_VEC(JL,IP_NO3) = JVAL_NO3_2d(IDX(JL),JK)
      END DO
    ENDIF
    IF (ASSOCIATED(JVAL_O3)) THEN
      DO JL=1,NVEC
        JX_VEC(JL,IP_O3) = JVAL_O3_2d(IDX(JL),JK)
      END DO
    ENDIF

  END SUBROUTINE SET_PHOTO
!-----------------------------------------------------------------------------
END MODULE MESSY_GMXE_AERCHEM_LIQ
