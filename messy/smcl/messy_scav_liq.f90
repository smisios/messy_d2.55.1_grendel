MODULE MESSY_SCAV_LIQ


!  AUTHOR: HOLGER TOST, MAINZ, JULY 2004
!
! - IN THIS MODULE THE REACTION RATES WHICH ARE NOT DEFINED DIRECTLY IN THE EQUATION FILE
!   ARE CALCULATED
! - THE FUNCTIONS APPEARING IN THE EQUATION FILE ARE DEFINED IN THIS MODULE
! - THE CALCULATION OF THE TRANSFER COEFFICIENTS VIA 2 PARAMETRIZATIONS IS INCLUDED HERE:
!   - PARAMETRIZATION FOLLOWING ROELOFS PHD THESIS WITH DIFFUSION COEFFICIENTS FROM FULLER, 86
!   - PARAMETRIZATION FOLLOWING SCHWARTZ WITH DIFFUSION COEFFICIENTS CALCULATED FROM BOLTZMANN
!     MOLECULAR SPEED AND FREE PATH LENGTH
! - TESTFACS SHIFT THE DISSOCIATION EQUILIBRIA TOWARDS A NUMERICAL STABLE, BUT FAST INTEGRATING CODE
! - THE ALLOCATE SCAV_VALUE - ROUTINE WHICH IS CALLED DURING INITIALIZATION IS THE DATA ROUTINE FOR 
!   THE CALCULATION OF THE TRANSFER COEFFICIENT
! - THE SUBROUTINE SET_EASY_SCAV IS ALSO ONLY CALLED DURING INITIALIZATION AND DEFINES THE SCAV_VALUES 
!   FOR THE EASY_SCAV_APPROACH (CONSTANT SCAVENGING COEFFICIENT)
!
! THIS MODULE IS USED WITHIN KPP AS ALL REACTION RATES ARE DEFINED HERE !!!!!

! THIS CONTAINS THE PROCESSOR DIRECTIVE IF A VECTORIZED SOLVER IS USED
!#include "messy_scav_l_vec.inc"


  USE MESSY_SCAV_L_KPP,             ONLY: DP, LSPEC=>NSPEC
  USE MESSY_SCAV_MEM,               ONLY: LSPEC_GAS
  USE MESSY_SCAV_INP_KPP
  USE MESSY_MAIN_TOOLS_KINETICS
  IMPLICIT NONE

  PRIVATE 

  SAVE

!  PRIVATE :: JX, ip_h2o2, ip_hno3, ip_o3

  PUBLIC :: LIQ_PHYSC_INIT
  PUBLIC :: EQUIL_DISSOCIATION_COEFF
  PUBLIC :: TRANSFER_COEFFICIENT_LIQ
  PUBLIC :: K_ARR
!  PUBLIC :: K_SIV_H2O2 
  PUBLIC :: ALLOC_SCAV_VALUES_L
  PUBLIC :: SET_EASY_SCAV_COEFFICIENTS, DEFINE_KS_FOR_EFFECTIVE_HENRY 
  PUBLIC :: SET_PHOTO
  PUBLIC :: ALLOC_CHEM_FIELDS, DEALLOC_CHEM_FIELDS
  PUBLIC :: VMEAN_F


  REAL(DP), PUBLIC, POINTER, DIMENSION(:)   :: JX

  REAL(dp), PUBLIC                          :: TEMP, CV_l
  REAL(dp), PUBLIC, POINTER, DIMENSION(:)   :: k_exf, k_exb
  REAL(DP), PUBLIC                          :: K_EXF_N2O5, K_EXF_CLNO3, K_EXF_BRNO3
  REAL(DP), PUBLIC                          :: K_EXF_Hg, K_EXF_RGM

  ! VARIABLES WITH VARIABLE VECTOR LENGTH
  REAL(DP), PUBLIC, POINTER, DIMENSION(:,:) :: VK_EXF, VK_EXB
  REAL(DP), PUBLIC, POINTER, DIMENSION(:)   :: VCV_L, VK_EXF_N2O5,          &
                                               VK_EXF_CLNO3, VK_EXF_BRNO3,  &
                                               VK_EXF_Hg, VK_EXF_RGM
  REAL(DP), PUBLIC, POINTER, DIMENSION(:)   :: VTEMP, VPRESS
  REAL(DP), PUBLIC, POINTER, DIMENSION(:,:) :: JX_VEC

! PHOTOLYSIS_VALUES
  REAL(DP), PUBLIC, DIMENSION(:,:,:), POINTER   :: JVAL_H2O2
  REAL(DP), PUBLIC, DIMENSION(:,:), POINTER     :: JVAL_H2O2_2d
  REAL(DP), PUBLIC, DIMENSION(:,:,:), POINTER   :: JVAL_HNO3
  REAL(DP), PUBLIC, DIMENSION(:,:), POINTER     :: JVAL_HNO3_2d
  REAL(DP), PUBLIC, DIMENSION(:,:,:), POINTER   :: JVAL_O3
  REAL(DP), PUBLIC, DIMENSION(:,:), POINTER     :: JVAL_O3_2d

 ! TEMPERATURE , PRESSURE, DENSITY
  REAL(DP), PUBLIC, ALLOCATABLE, DIMENSION(:) :: VRHO, VRDRAD
  REAL(DP), PUBLIC, ALLOCATABLE, DIMENSION(:) :: VLIQ_LWC, VLIQ_LWC_EFF
  REAL(DP), ALLOCATABLE, DIMENSION(:,:)       :: ALPHA_V, VMEAN_V
 
  ! SWITCH FOR PARAMETRISATION OF TRANSFER COEFFICIENT FOLLOWING
  ! SCHWARTZ (ON) OR ROELOFS(OFF)
  LOGICAL,  PUBLIC ::  USE_SCHWARTZ

  INTEGER :: NVEC, SET

  INTRINSIC :: EXP, TINY
 

CONTAINS
!==============================================================

!-------------------------------------------------------------------------------
     ELEMENTAL REAL(DP) FUNCTION K_ARR (K_298,TDEP,TEMP)
       ! ARRHENIUS FUNCTION

       REAL, INTENT(IN) :: K_298 ! K AT T = 298.15K
       REAL, INTENT(IN) :: TDEP  ! TEMPERATURE DEPENDENCE
       REAL(DP), INTENT(IN) :: TEMP  ! TEMPERATURE
       
       INTRINSIC :: REAL
       
       K_ARR = REAL(K_298,DP) * EXP(REAL(TDEP,DP)*(1._DP/TEMP-3.3540E-3_DP)) ! 1/298.15=3.3540E-3

     END FUNCTION K_ARR

!-------------------------------------------------------------------------------

  SUBROUTINE ALLOC_CHEM_FIELDS(LPROMA)


!#ifdef SCAV_VECTOR
!    USE MESSY_SCAV_L_KPP, ONLY : C, RCONST, VAR, RAD, FIX, F,  &
!                                 NVAR,      NRAD,     NFIX,    &
!                                 NVARST,    NRADST,   NFIXST,  &
!                                 NVECT,     NREACT
!
!#endif
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: LPROMA

!#ifdef SCAV_VECTOR
!    NVECT = LPROMA
!#endif    
    NVEC  = LPROMA

    ALLOCATE(VTEMP(NVEC))
    ALLOCATE(VPRESS(NVEC))
    ALLOCATE(VRHO(NVEC))
    ALLOCATE(VLIQ_LWC(NVEC))
    ALLOCATE(VLIQ_LWC_EFF(NVEC))
    ALLOCATE(VRDRAD(NVEC))
    ALLOCATE(JX_VEC(NVEC,IP_MAX))
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

!#ifdef SCAV_VECTOR
!    ALLOCATE (C(NVEC,LSPEC))
!    ALLOCATE (RCONST(NVEC,NREACT))
!    
!    VAR => C(:,1+NVARST:NVARST+NVAR)
!    RAD => C(:,1+NRADST:NRADST+NRAD)
!    FIX => C(:,1+NFIXST:NFIXST+NFIX)
!    F   => C(:,1+NFIXST:NFIXST+NFIX)
!#endif
  END SUBROUTINE ALLOC_CHEM_FIELDS
!-------------------------------------------------------------------
  SUBROUTINE DEALLOC_CHEM_FIELDS

!#ifdef SCAV_VECTOR
!    USE MESSY_SCAV_L_KPP, ONLY : C, RCONST
!#endif
    IMPLICIT NONE

    DEALLOCATE(VTEMP)
    DEALLOCATE(VPRESS)
    DEALLOCATE(VRHO)
    DEALLOCATE(VLIQ_LWC)
    DEALLOCATE(VLIQ_LWC_EFF)
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

!#ifdef SCAV_VECTOR
!    DEALLOCATE(C)
!    DEALLOCATE(RCONST)
!#endif
  END SUBROUTINE DEALLOC_CHEM_FIELDS
!--------------------------------------------------------------------------------------------
  
  SUBROUTINE SET_PHOTO(IDX,JK,JROW)

    INTEGER, INTENT(IN) :: IDX(NVEC),JK,JROW
    INTEGER :: JL

    INTRINSIC :: ASSOCIATED

    JX_VEC(1:NVEC,:) = 0._DP

    IF (ASSOCIATED(JVAL_H2O2)) THEN
      DO JL=1,NVEC
        JX_VEC(JL,IP_H2O2) = JVAL_H2O2_2d(IDX(JL),JK)
      END DO
    ENDIF
    IF (ASSOCIATED(JVAL_HNO3)) THEN
      DO JL=1,NVEC
        JX_VEC(JL,IP_HNO3) = JVAL_HNO3_2d(IDX(JL),JK)
      END DO
    ENDIF
    IF (ASSOCIATED(JVAL_O3)) THEN
      DO JL=1,NVEC
        JX_VEC(JL,IP_O3) = JVAL_O3_2d(IDX(JL),JK)
      END DO
    ENDIF

  END SUBROUTINE SET_PHOTO

!-------------------------------------------------------------------------------

    SUBROUTINE TRANSFER_COEFFICIENT_LIQ(PROC)

    USE MESSY_SCAV_MEM,           ONLY: DHT, HS0, L_SPEC, KPP_L_IDX, GAS_IDX &
                                      , loverwrite_alpha
    USE MESSY_SCAV_AER,           ONLY: LAMBDA, AIRVISC, REY_RAIN
    USE MESSY_SCAV_L_KPP


    IMPLICIT NONE
    INTEGER, INTENT(IN) :: PROC     ! SWITCH FOR DESTINCTION BETWEEN CLOUD AND RAIN
                                    ! 1 = CLOUD, 2 = RAIN
    INTEGER  :: JT, K, JL, IDX
    REAL(DP) :: ZDGAIR, ZRU, ZNRE, ZNSC
    REAL(DP) :: ZRDRAD(NVEC), ZNSH(NVEC)
    REAL(DP) :: ZTX(NVEC)
    REAL(DP) :: ZHETT
    REAL(DP) :: VISC(NVEC), ZLAMBDA(NVEC)
    REAL(DP) :: REY_R
    REAL(DP) :: KC, P1, P3, FLUX
    REAL(DP) :: FL1(NVEC), FL2(NVEC)
    REAL(DP) :: P2(NVEC)
    REAL(DP) :: HENRY(LSPEC,NVEC)
    INTEGER, PARAMETER  :: IS = 6     ! ITERATION STEP FOR KMT CALCULATION
    REAL(DP) :: ZKMT_I(IS,LSPEC,NVEC), KMT(NVEC)
    REAL(DP) :: RADFIELD(IS+1), DEL_NUM 
    REAL(DP) :: NUMBER(IS,NVEC)
!BS-27022007+
    REAL(DP) :: NUMBER_HELP(NVEC)
!BS-27022007-
    USE_SCHWARTZ = .TRUE.
 
   !CALCULATES THE TRANSFER COEFFICIENTS IN AND OUT OF THE WATER DROPLET
    VK_EXF(:,:)     = 0.
    VK_EXB(:,:)     = 0.
    VK_EXF_N2O5(:)  = 0.
    VK_EXF_CLNO3(:) = 0.
    VK_EXF_BRNO3(:) = 0.
    VK_EXF_Hg(:)    = 0.
    VK_EXF_RGM(:)   = 0.
    
    KMT(:) = 0._DP

    SELECT CASE(PROC)

    CASE(1)

      IF (.NOT.USE_SCHWARTZ) THEN   
        ZDGAIR=0.133_DP
        DO JL=1,NVEC
          ZRDRAD(JL) = VRDRAD(JL)
          ZRU        = 9.58 * (1._DP - EXP (- (ZRDRAD(JL)/0.885)**1.147 ) )
          ZNRE       = 20.* ZRDRAD(JL) * ZRU / ZDGAIR
          ZNSC       = 1._DP
          ZNSH(JL)   = 1._DP + 0.3 * (ZNRE**0.5) * (ZNSC**(1._DP/3._DP))
        ENDDO
      ELSE
!!$        ZLAMBDA(1:NVEC) = 2.28E-5_DP * VTEMP(1:NVEC) / VPRESS(1:NVEC)
!!$          ! RECALCULATION RADIUS FROM MM TO M
!!$        ZRDRAD(1:NVEC)  = 1.E-3_DP * VRDRAD(1:NVEC)
        DO JL=1,NVEC
          ZLAMBDA(JL) = LAMBDA(VTEMP(JL),VPRESS(JL))
          ! RECALCULATION RADIUS FROM MM TO M
          ZRDRAD(JL)  = 1.E-3_DP * VRDRAD(JL)
        ENDDO
        IF (loverwrite_alpha) THEN
           CALL SCAV_ALPHA_V
        ELSE
           CALL SCAV_ALPHA_V_GEN
        END IF
        CALL SCAV_VMEAN_V
      ENDIF

      ZTX(1:NVEC) = 1._DP/VTEMP(1:NVEC) - 1._DP/298._DP
      HENRY(:,:) = 0.

      DO JT=1,LSPEC_GAS
        IDX = KPP_L_IDX(SET)%gas_spec(jt,gas_idx)

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

          ! HT   CALCULATION OF THE PARAMETER CONTROLLING THE DISSOLUTION IN THE
          !      WATER PHASE FOLLOWING ROELOFS

          IF (.NOT.USE_SCHWARTZ) THEN
            DO JL=1,NVEC
              KMT(JL) = VMEAN_V(JL,IDX) * ZLAMBDA(JL) * 1.E4_DP/ &
                        ( (VRDRAD(JL) *0.1)**2) * ZNSH(JL)
            ENDDO
          ELSE
            ! HT     ALTERNATIVE FOLLOWING SCHWARTZ           
            DO JL=1,NVEC
              KMT(JL) = VMEAN_V(JL,IDX) /  (ZRDRAD(JL) * (ZRDRAD(JL) / &
                ZLAMBDA(JL) + 4._DP / (3._DP * ALPHA_V(JL,IDX) )))
            ENDDO
          ENDIF
         
        !  SETTING OF EXCHANGE FUNCTION            
          
!!$          VK_EXF(1:NVEC,JT) = KMT(1:NVEC) * VLIQ_LWC_EFF(1:NVEC) * 1.E-3_DP
!!$          VK_EXB(1:NVEC,JT) = KMT(1:NVEC) * HENRY(JT,1:NVEC) 
          DO JL=1,NVEC
            VK_EXF(JL,IDX) = KMT(JL) * VLIQ_LWC_EFF(JL) * 1.E-3_DP
            VK_EXB(JL,IDX) = KMT(JL) * HENRY(IDX,JL)
          ENDDO

        ENDIF
      ENDDO

    CASE(2)
      !   RAIN TRANSFER COEFFICIENT
      !   ACCORDING TO FROESLING(1938), BIRD (1960)
      !   KC[CM/S] -> KMT[1/S] => KC * 3/R  FROM GEOMETRY (SURFACE TO VOLUME)
      !   KC = DG/DIAMETER * (2 + 0.6 *(SQRT(REY_R)) &
      !                           * ( AIRVISC / (RHO * DG) )**1/3 )
      !   DG = LAMBDA * VMEAN / 3. 
      ! MEAN FREE PATH IN M
      DO JL = 1,NVEC
        ZLAMBDA(JL) = LAMBDA(VTEMP(JL), VPRESS(JL) )
        ! RECALCULATION RADIUS FROM MM TO M
        ZRDRAD(JL)  = 1.E-3_DP * VRDRAD(JL)
        ZTX(JL)     = 1._DP / VTEMP(JL) - 1._DP/298._DP
        VISC(JL)    = AIRVISC(VTEMP(JL))
      ENDDO
      CALL SCAV_VMEAN_V
      HENRY(:,:) = 0.
        
      ! CALCULATE VALUES INDEPENDENT OF THE SPECIES

      ! BUILD RADIUS FIELD IN MM
      RADFIELD(IS+1) = 0._DP

      RADFIELD(1) = 5.00E-3_DP
      RADFIELD(2) = 2.00E-3_DP
      RADFIELD(3) = 1.00E-3_DP
      RADFIELD(4) = 0.50E-3_DP
      RADFIELD(5) = 0.20E-3_DP
      RADFIELD(6) = 0.10E-3_DP

      DO JL=1,NVEC
        FLUX        = ZRDRAD(JL) * 1E3
        FLUX        = (FLUX / 0.3659_DP) ** (1._DP/0.21_DP)

        FL1(JL)     = 2.8E-5 * FLUX**0.324_DP
        FL2(JL)     = - 98.5_DP * FLUX**(-0.522_DP)
      ENDDO

      DO JT=1,LSPEC_GAS
        IDX = KPP_L_IDX(SET)%gas_spec(jt,gas_idx)
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
        ENDIF
      ENDDO
      
      DO K = 1,IS
        DO JL=1,NVEC
          REY_R     = REY_RAIN(RADFIELD(K), VRHO(JL), VISC(JL) )
          ! RADFIELD IN M
          ! NUMBER IN CM^-3 CM^-1 OF ALL DROPLETS LARGER THAN RADFIELD(K)
        ! WARNING: DEL_NUM IS NOT THE ABSOLUTE NUMBER BUT DN/DR
        ! BEST
          DEL_NUM   = FL1(JL) * (2._DP * RADFIELD(K) * 10._DP)**(-1.75_DP) *&
                      EXP( FL2(JL) * (RADFIELD(K) * 10._DP * 2._DP)**2.25)
          NUMBER(K,JL) = DEL_NUM * (RADFIELD(K) - RADFIELD(K+1))
          P2(JL) = SQRT(REY_R)
        ENDDO
        ! NOW INDIVIDUAL FOR ALL SPECIES
        DO JT = 1,LSPEC_GAS
          IDX = KPP_L_IDX(SET)%gas_spec(jt,gas_idx)
          IF (HS0(IDX) > 0._dp) THEN        
            DO JL = 1,NVEC
              if ((3._DP * VISC(JL) / (VRHO(JL) *                  &
                   VMEAN_V(JL,IDX) * ZLAMBDA(JL)) ) < tiny(0._dp) ) &
                   print*,jt, VISC(JL),VRHO(JL), VMEAN_V(JL,IDX), ZLAMBDA(JL)
              P1 = (3._DP * VISC(JL) / (VRHO(JL) * &
                   VMEAN_V(JL,IDX) * ZLAMBDA(JL)) )**(1._DP/3._DP)             
              P3 = 2._DP + 0.6_DP * P2(JL) * P1
              KC = VMEAN_V(JL,IDX) * ZLAMBDA(JL) / (6._DP * RADFIELD(K)) * P3
              ZKMT_I(K,IDX,JL) = 3._DP * KC / RADFIELD(K)
            ENDDO
          ENDIF
        END DO
      END DO
      !  SETTING OF EXCHANGE FUNCTION            
        
      DO JT = 1,LSPEC_GAS
        IDX = KPP_L_IDX(SET)%gas_spec(jt,gas_idx)
        IF (HS0(IDX) > 0._dp) THEN   
          KMT(:) = 0._DP
          DO K = 1,IS
!BS-27022007+
            DO JL=1,NVEC
              KMT(JL) = KMT(JL) + &
                           ZKMT_I(K,IDX,JL) * NUMBER(K,JL)
            ENDDO
          ENDDO
          NUMBER_HELP(:)=0._DP
          DO K = 1,IS
            DO JL = 1,NVEC
              NUMBER_HELP(JL)=NUMBER_HELP(JL)+NUMBER(K,JL)
            ENDDO
          ENDDO
          DO JL = 1,NVEC
            KMT(JL) = KMT(JL) / NUMBER_HELP(JL)
          ENDDO
!BS-27022007-
          VK_EXF(1:NVEC,IDX) = KMT(1:NVEC) * VLIQ_LWC_EFF(1:NVEC) * 1.E-3_DP
          VK_EXB(1:NVEC,IDX) = KMT(1:NVEC) * HENRY(IDX,1:NVEC)           
        ENDIF
      ENDDO

    END SELECT

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

!-------------------------------------------------------------------------

  SUBROUTINE LIQ_PHYSC_INIT(LON_LWC, LON_LWC2, LON_RDRAD,&
                            LON_TEMP, LON_PRESS, JS)

! INITIALIZATION OF THE BOX VALUES FOR TEMPERATURE, PRESSURE, LIQUID WATER, 
! CHEMICAL CONCENTRATION OF WATER, DROPLET DIAMETER
    USE MESSY_MAIN_CONSTANTS_MEM,    ONLY: AVO => N_A, M_AIR, R_GAS
    USE MESSY_SCAV_MEM,              ONLY: L_SPEC
    USE MESSY_SCAV_L_KPP,            ONLY: IND_H2O_L

INTRINSIC :: REAL
    REAL(DP), INTENT(IN)  :: LON_LWC(NVEC) ,LON_LWC2(NVEC), LON_RDRAD(NVEC),&
                             LON_TEMP(NVEC), LON_PRESS(NVEC)
    INTEGER,  INTENT(IN)  :: JS
    INTEGER :: JL

    SET = JS
    DO JL=1,NVEC
      VLIQ_LWC_EFF(JL) = LON_LWC(JL)
      VLIQ_LWC(JL)     = LON_LWC2(JL)
      IF (VLIQ_LWC(JL)*REAL(IND_H2O_L,DP) > 0._DP) &
        L_SPEC(IND_H2O_L,JL) = 55.5_DP * AVO * VLIQ_LWC(JL) * 1.E-6_DP
      VRDRAD(JL)   = LON_RDRAD(JL)
      VTEMP(JL)    = LON_TEMP(JL)
      VPRESS(JL)   = LON_PRESS(JL)
      VRHO(JL)     = VPRESS(JL) * M_AIR / (R_GAS * VTEMP(JL) * 1.E3_DP)
    ENDDO

  END SUBROUTINE LIQ_PHYSC_INIT
!--------------------------------------------------------------------------


!======================================================================
  SUBROUTINE ALLOC_SCAV_VALUES_L

    USE MESSY_SCAV_L_KPP
    USE MESSY_SCAV_MEM,     ONLY: HS0, DHT


    HS0(IND_OH)      = 3.0E1_DP          !RS_MOCCA 
    HS0(IND_HO2)     = 3.9E3_DP          !RS_MOCCA 
    HS0(IND_H2O2)    = 1.E+05_DP         !RS_MOCCA
    HS0(IND_O2)      = 1.3E-3_DP         !RS_LIDE95 
    HS0(IND_O3)      = 1.2E-02_DP
    
    HS0(IND_NH3)     = 58._DP            !RS_MOCCA
    HS0(IND_NO)      = 1.9E-03_DP        !RS_LIDE
    HS0(IND_NO2)     = 6.4E-03_DP        !RS_MOCCA 
    HS0(IND_NO3)     = 2._DP             !RS_MOCCA
    HS0(IND_HONO)    = 4.9E+01_DP        !RS_MOCCA
    HS0(IND_HNO3)    = 2.5E6_DP/1.5E1_DP !RS_MOCCA_EXAKT
    HS0(IND_N2O5)    = 1.40_DP
    HS0(IND_HNO4)    = 1.2E4_DP
    
    HS0(IND_CH4)     = 1.4E-3_DP         !RS_LIDE
    HS0(IND_HCHO)    = 7.0E+03_DP        !RS_MOCCA
    HS0(IND_HCOOH)   = 3.7E+03_DP        !RS_MOCCA
    HS0(IND_CH3OH)   = 2.20E2_DP
    HS0(IND_CH3OOH)  = 3.0E+02_DP        !RS_MOCCA CH3OOH
    HS0(IND_CH3O2)   = 6._DP             !RS_JACOB
    HS0(IND_CO)      = 9.9E-4_DP         !RS_LIDE
    HS0(IND_CO2)     = 3.1E-02_DP        !RS_MOCCA

    HS0(IND_C2H6)    = 1.9E-3_DP         !RS_LIDE
    HS0(IND_C2H4)    = 4.7E-3_DP         !RS_SEINFELD
    HS0(IND_CH3CO2H) = 5.50E+3_DP
    HS0(IND_CH3CHO)  = 1.14E1_DP
    HS0(IND_C2H5OOH)   = 3.4E2_DP
    HS0(IND_C2H5O2)    = 6._DP
    HS0(IND_CH3CO3H) = 8.4E2_DP          !RS_OSULLIVAN
    HS0(IND_CH3CO3)  = 0.1_DP            !RS
    HS0(IND_PAN)     = 4.1_DP            !RS_KAMES
    
    HS0(IND_C3H6)    = 4.8E-3_DP
    HS0(IND_C3H8)    = 1.5E-3_DP         !RS_LIDE   
    HS0(IND_CH3COCH3)= 3.52E1_DP     
    HS0(IND_MGLYOX)  = 3.2E4_DP          !RS_ZHOU
    HS0(IND_HYPERACET)  = 1.3E2_DP          !RS_SNIDER
    HS0(IND_CH3COCH2O2) = 6._DP
    HS0(IND_ACETOL)  = 1.3E2_DP          !RS_SNIDER
    HS0(IND_LHOC3H6OOH) = 3.4E2_DP
    HS0(IND_LHOC3H6O2)  = 6_DP
!   HS0(IND_NOA)     = 1.0E3_DP          !RS_KAMES
    HS0(IND_MPAN)    = 1._DP             !RS_KAMES
    HS0(IND_IC3H7OOH)   = 3.4E2_DP
    HS0(IND_IC3H7O2)    = 6._DP
    HS0(IND_IC3H7NO3)  = 6.2E-1_DP         !RS_HAUFF

    HS0(IND_NC4H10)   = 1.1E-3_DP
    HS0(IND_LC4H9OOH) = 3.4E2_DP
    HS0(IND_LC4H9O2)  = 6._DP
    HS0(IND_MEK)     = 2.0E1_DP          !RS_STAUDINGER
    HS0(IND_LMEKO2)  = 6._DP
    HS0(IND_LMEKOOH) = 3.4E2_DP
    HS0(IND_MVK)     = 2.1E1_DP          !RS_ALLEN
    HS0(IND_MVKO2)   = 6._DP
    HS0(IND_MVKOOH)  = 3.4E2_DP
    HS0(IND_BIACET)  = 7.4E1_DP          !RS_BETTERTON
    HS0(IND_LC4H9NO3) = 1._DP             !RS_KAMES

    HS0(IND_C5H8)        = 2.8E-2_DP     !RS_KARL
!    HS0(IND_NISOPOOH)    = 6.E3_DP       !RS_SHEPSON
!    HS0(IND_NISOPO2)     = 6.E1_DP
!    HS0(IND_ISOPBO2)     = 6._DP
!    HS0(IND_NISOOH)      = 3.4E2_DP
!    HS0(IND_ISOPBOOH)    = 3.4E2_DP
!    HS0(IND_ISOPBNO3)    = 4.5E-1_DP  
!    HS0(IND_ISOPACDOOH)  = 3.4E2_DP
!    HS0(IND_ISOPACDNO3)  = 4.5E-1_DP
!    HS0(IND_ISOPACDO2)   = 6._DP
!    HS0(IND_NC4CHO)      = 5._DP
!    HS0(IND_NISO3)       = 6._DP
!    HS0(IND_C5PAN18)     = 5._DP         ! PAN LIKE

    HS0(IND_HCL)     = 2._DP/1.7E0_DP  !RS_MOCCA_EXAKT
    HS0(IND_HOCL)    = 6.7E2_DP        !RS_MOCCA
    HS0(IND_CLNO2)   = 4.6E-2_DP
    HS0(IND_CLNO3)   = 1.E11_DP
    HS0(IND_CL2)     = 9.1E-2_DP       !RS_MOCCA
    HS0(IND_CL)      = 2.0E-1_DP
!    HS0(IND_CLO)     = 
!    HS0(IND_OCLO)    = 
    HS0(IND_HBR)     = 1.3E0_DP        !RS_MOCCA
    HS0(IND_HOBR)    = 9.3E1_DP        !RS_MOCCA
    HS0(IND_BRNO2)   = 3.0E-1_DP
    HS0(IND_BRNO3)   = 1.E11_DP
    HS0(IND_BR2)     = 7.6E-1_DP       !RS_MOCCA
    HS0(IND_BRCL)    = 9.4E-1_DP       !RS_MOCCA
    HS0(IND_BR)      = 1.2_DP
!    HS0(IND_BRO)     = 

    HS0(IND_HI)      = 2.2E4_DP
    HS0(IND_HOI)     = 4.5E2_DP        !RS_MOCCA 
    HS0(IND_I2)      = 3._DP           !RS_MOCCA
    HS0(IND_ICL)     = 1.1E2_DP        !RS_MOCCA 
    HS0(IND_IBR)     = 2.4E1_DP        !RS_MOCCA
    HS0(IND_I)       = 8.0E-2_DP
    HS0(IND_IO)      = 4.5E2_DP        !RS_MOCCA
!    HS0(IND_INO2)    = 
!    HS0(IND_I2O2)    = 
!    HS0(IND_INO3)    = 
!    HS0(IND_CH3I)    = 
!    HS0(IND_CH2I2)   = 
!    HS0(IND_CH2CLI)  = 
!    HS0(IND_C3H7I)   = 
!    HS0(IND_OIO)     = 
!    HS0(IND_HIO3)    = 

    HS0(IND_SO2)     = 1.2E0_DP          !RS_MOCCA
    HS0(IND_H2SO4)   = 1.0E11_DP
    HS0(IND_DMS)     = 4.8E-1_DP       !RS_DEBRYUN
    HS0(IND_DMSO)    = 5.E4_DP         !RS_MOCCA
    HS0(IND_CH3SO3H) = 1.E8_DP
!    HS0(IND_CH3SO2)  = 
!    HS0(IND_CH3SO3)  = 

    HS0(IND_Hg)       = 0.13_dp   ! (see gas.tex of MECCA)
    HS0(IND_HgO)      = 2.4E7_dp  ! (see gas.tex of MECCA)
    HS0(IND_HgCl)     = 2.4E7_dp
    HS0(IND_HgCl2)    = 2.4E7_dp
    HS0(IND_HgBr)     = 2.4E7_dp
    HS0(IND_HgBr2)    = 2.4E7_dp
    HS0(IND_ClHgBr)   = 2.4E7_dp
    HS0(IND_BrHgOBr)  = 2.4E7_dp
    HS0(IND_ClHgOBr)  = 2.4E7_dp

! TEMPERATURE DEPEDENCE DHT OF HENRY COEFFICIENTS
    DHT(IND_O2)      = 1500._DP
    DHT(IND_O3)      = 2560._DP !RS_MOCCA

    DHT(IND_OH)      = 4300._DP !RS_MOCCA
    DHT(IND_HO2)     = 5900._DP !RS_MOCCA 
    DHT(IND_H2O2)    = 6338._DP !RS_MOCCA

    DHT(IND_NH3)     = 4085._DP !RS_MOCCA
    DHT(IND_NO)      = 1480._DP !RS_LIDE
    DHT(IND_NO2)     = 2500._DP !RS_MOCCA
    DHT(IND_NO3)     = 2000._DP !RS_MOCCA
    DHT(IND_HONO)    = 4780._DP !RS_MOCCA
    DHT(IND_HNO3)    = 8694._DP !RS_MOCCA_EXAKT
    DHT(IND_N2O5)    = 0._DP
    DHT(IND_HNO4)    = 6900._DP

    DHT(IND_CH4)     = 1600._DP
    DHT(IND_HCHO)    = 6425._DP !RS_MOCCA
    DHT(IND_HCOOH)   = 5700._DP !RS_MOCCA
    DHT(IND_CH3OH)   = 5390._DP
    DHT(IND_CH3OOH)  = 5322._DP !RS_MOCCA CH3OOH
    DHT(IND_CH3O2)   = 5600._DP !RS_JACOB
    DHT(IND_CO)      = 1600._DP
    DHT(IND_CO2)     = 2423._DP !RS_MOCCA

    DHT(IND_C2H6)    = 2300._DP
    DHT(IND_C2H4)    = 0._DP
    DHT(IND_CH3CO2H) = 5894._DP
    DHT(IND_CH3CHO)  = 6254._DP
    DHT(IND_C2H5OOH)   = 5322._DP
    DHT(IND_C2H5O2)    = 87._DP
    DHT(IND_CH3CO3H) = 5300._DP
    DHT(IND_CH3CO3)  = 0._DP
    DHT(IND_PAN)     = 0._DP

    DHT(IND_C3H6)    = 0._DP
    DHT(IND_C3H8)    = 2700._DP
    DHT(IND_CH3COCH3)= 3800._DP
    DHT(IND_MGLYOX)  = 0._DP
    DHT(IND_HYPERACET)  = 7500._DP
    DHT(IND_CH3COCH2O2) = 0._DP
    DHT(IND_ACETOL)  = 7500._DP
    DHT(IND_LHOC3H6OOH) = 5322._DP
    DHT(IND_LHOC3H6O2)  = 0._DP
!    DHT(IND_NOA)     = 0._DP
    DHT(IND_MPAN)    = 0._DP
    DHT(IND_IC3H7OOH)   = 5322._DP
    DHT(IND_IC3H7O2)    = 0._DP
    DHT(IND_IC3H7NO3)  = 0._DP

    DHT(IND_NC4H10)   = 0._DP
    DHT(IND_LC4H9OOH) = 5322._DP
    DHT(IND_LC4H9O2)  = 0._DP
    DHT(IND_MEK)     = 5000._DP
    DHT(IND_LMEKO2)   = 0._DP
    DHT(IND_LMEKOOH)  = 5322._DP
    DHT(IND_MVK)     = 7800._DP
    DHT(IND_MVKO2)   = 0._DP
    DHT(IND_MVKOOH)  = 5322._DP
    DHT(IND_BIACET)  = 5700._DP
    DHT(IND_LC4H9NO3) = 5800._DP

    DHT(IND_C5H8)        = 0._DP
!    DHT(IND_NISOPOOH)    = 5322._DP 
!    DHT(IND_NISOPO2)     = 0._DP
!    DHT(IND_ISOPBO2)     = 0._DP
!    DHT(IND_NISOOH)      = 5322._DP
!    DHT(IND_ISOPBOOH)    = 5322._DP
!    DHT(IND_ISOPBNO3)    = 0._DP
!    DHT(IND_ISOPACDOOH)  = 5322._DP 
!    DHT(IND_ISOPACDNO3)  = 0._DP
!    DHT(IND_ISOPACDO2)   = 0._DP
!    DHT(IND_NC4CHO)      = 0._DP
!    DHT(IND_NISO3)       = 0._DP
!    DHT(IND_C5PAN18)     = 0._DP

    DHT(IND_HCL)     = 9001._DP  !RS_MOCCA_EXAKT
    DHT(IND_HOCL)    = 5862._DP  !RS_MOCCA
    DHT(IND_CLNO2)   = 0._DP
    DHT(IND_CLNO3)   = 0._DP
    DHT(IND_CL2)     = 2500._DP  !RS_MOCCA
    DHT(IND_CL)      = 0._DP
!    DHT(IND_CLO)     =
!    DHT(IND_OCLO)    =

    DHT(IND_HBR)     = 10239._DP !RS_MOCCA
    DHT(IND_HOBR)    = 5862._DP  !RS_MOCCA 
    DHT(IND_BRNO2)   = 0._DP
    DHT(IND_BRNO3)   = 0._DP
    DHT(IND_BR2)     = 4094._DP  !RS_MOCCA
    DHT(IND_BRCL)    = 5600._DP  !RS_MOCCA
    DHT(IND_BR)      = 0._DP
!    DHT(IND_BRO)     =

    DHT(IND_HI)      = 9800._DP
    DHT(IND_HOI)     = 5862._DP  !RS_MOCCA 
    DHT(IND_I2)      = 4431._DP  !RS_MOCCA
    DHT(IND_ICL)     = 5600._DP  !RS_MOCCA 
    DHT(IND_IBR)     = 5600._DP  !RS_MOCCA
    DHT(IND_I)       = 0._DP
    DHT(IND_IO)      = 5862._DP  !RS_MOCCA
!    DHT(IND_INO2)    =
!    DHT(IND_I2O2)    =
!    DHT(IND_INO3)    =
!    DHT(IND_CH3I)    =
!    DHT(IND_CH2I2)   =
!    DHT(IND_CH2CLI)  =
!    DHT(IND_C3H7I)   =
!    DHT(IND_OIO)     =
!    DHT(IND_HIO3)    =

    DHT(IND_SO2)     = 3120._DP !RS_MOCCA
    DHT(IND_H2SO4)   = 0._DP
    DHT(IND_DMS)     = 3100._DP
    DHT(IND_DMSO)    = 6425._DP  !RS_MOCCA
    DHT(IND_CH3SO3H) = 0._DP
!    DHT(IND_CH3SO2)  =
!    DHT(IND_CH3SO3)  =

  END SUBROUTINE ALLOC_SCAV_VALUES_L

!-------------------------------------------------------------------------------
  SUBROUTINE SCAV_ALPHA_V

      ! CALCULATION OF ACCOMMODATION COEFFICIENTS NEEDED FOR K_EX CALCULATIONS
      USE MESSY_MAIN_CONSTANTS_MEM, ONLY: CAL2J, R_GAS, T0
      USE MESSY_SCAV_L_KPP

      REAL(DP) :: ZTCORR, ZTEMP
      INTEGER  :: JL
      
      INTRINSIC EXP

      ! INIT
      ! STANDARD VALUE
      ALPHA_V(:,:)      = 0.1_DP
      DO JL=1,NVEC
        ZTEMP= VTEMP(JL)

        ZTCORR=1./ZTEMP-1./T0
      ! FOR SOME SPECIES A TEMPERATURE DEPENDENCE IS CALCULATED/ESTIMATED
      ! USING:
      ! ALPHA = 1./(1.+1./(1./(1./ALPHA(T0)-1.)*EXP((-DELTAH/R_GAS)*ZTCORR)))
      ! EXPERIMENTALLY DETERMINED VALUES (FROM MOCCA)

        ALPHA_V(JL,IND_O2)      = &
          1.0/(1.0_DP+1.0/(1.0/(1.0/1.0E-2 -1.0_DP)*EXP(2000.*ZTCORR)))
        ALPHA_V(JL,IND_O3)      = 2.0E-03_DP

        ALPHA_V(JL,IND_OH)      = 1.0E-2_DP
        ALPHA_V(JL,IND_HO2)     = 2.0E-01_DP
        ALPHA_V(JL,IND_H2O2)    = &
          1./(EXP(-26.E3/(R_GAS*ZTEMP)+107.8456/R_GAS) + 1.0_DP)

        ALPHA_V(JL,IND_NH3)     = 6.0E-02_DP
        ALPHA_V(JL,IND_NO2)     = 1.5E-03_DP
        ALPHA_V(JL,IND_NO3)     = 4.0E-02_DP
        ALPHA_V(JL,IND_HONO)    = 4.0E-02_DP
        ALPHA_V(JL,IND_HNO3)    = 5.0E-01_DP
        ALPHA_V(JL,IND_N2O5)    = 1.0E-01_DP
        ALPHA_V(JL,IND_HNO4)    = 0.1_DP

        ALPHA_V(JL,IND_HCHO)    = 4.0E-02_DP
        ALPHA_V(JL,IND_HCOOH)   =                                 &
          1.0_DP/(EXP(-7.9E3*CAL2J/(R_GAS*ZTEMP)+34.9*CAL2J/R_GAS)+1.0_DP)
        ALPHA_V(JL,IND_CH3OOH)  =                                 &
          1.0_DP/(EXP(-6.5E3*CAL2J/(R_GAS*ZTEMP)+32.5*CAL2J/R_GAS) + 1.0_DP)
        ALPHA_V(JL,IND_CO2)     =                                 &
          1.0_DP/(1.0_DP+1.0_DP/(1.0_DP/(1.0_DP/1.0E-2 -1.0_DP) * &
          EXP(2000.*ZTCORR)))
        ALPHA_V(JL,IND_CH3O2)   =                                 &
          1.0_DP/(1.0_DP+1.0_DP/(1.0_DP/(1.0_DP/1.0E-2 -1.0_DP) * &
          EXP(2000.*ZTCORR)))

        ALPHA_V(JL,IND_CH3CO2H) = 1.9E-02_DP
        ALPHA_V(JL,IND_CH3CHO)  = 3.0E-02_DP

        ALPHA_V(JL,IND_CH3COCH3)= 1.9E-02_DP

        ALPHA_V(JL,IND_HCL)     = &
          1.0_DP/(EXP(-3.072E3/ZTEMP + 1.283E1_DP)+1.0_DP)
        ALPHA_V(JL,IND_HOCL)    = 5.0E-01_DP
        !T=290: 0.096, T=270: 0.190
        !OLD VALUE: ALPHA(IND_CLNO3)   = 1.0E-01_DP
        ALPHA_V(JL,IND_CLNO3)   = 0.108_DP !MZ_RS_20040308 NEW REF1647
        ALPHA_V(JL,IND_CL2)     = &
          1.0_DP/(EXP(-1.3E4*CAL2J/(R_GAS*ZTEMP)+50.*CAL2J/R_GAS)+1.0_DP)

        ALPHA_V(JL,IND_HBR)     = &
          1.0_DP/(EXP(-3.94E3/ZTEMP + 1.664E1_DP) + 1.0_DP)
        !T=290K: 0.017, T=270K: 0.130
        ALPHA_V(JL,IND_HOBR)    = 5.0E-01_DP
        !OLD VALUE: ALPHA(IND_BRNO3)   = 8.0E-01_DP
        ALPHA_V(JL,IND_BRNO3)   = 0.063_DP !MZ_RS_20040308 NEW REF1647
        ALPHA_V(JL,IND_BR2)     = &
          1.0_DP/(EXP(-1.3E4*CAL2J/(R_GAS*ZTEMP)+50.*CAL2J/R_GAS)+1.0_DP)
        ALPHA_V(JL,IND_BRCL)    = 0.33_DP  ! #840 ALPHA(IND_CL2)

        ALPHA_V(JL,IND_HI)      = &
          1.0_DP/(EXP(-4.13E3/ZTEMP + 1.715E1_DP)+1.0_DP)
        ALPHA_V(JL,IND_HOI)     = ALPHA_V(JL,IND_HOBR)
        ALPHA_V(JL,IND_INO3)    =                                 &
          1.0_DP/(1.0_DP+1.0_DP/(1.0_DP/(1.0_DP/1.0E-1 -1.0_DP) * &
          EXP(2000.*ZTCORR)))
        ALPHA_V(JL,IND_I2)      =                                 &
          1.0_DP/(1.0_DP+1.0_DP/(1.0_DP/(1.0_DP/1.0E-2 -1.0_DP) * &
          EXP(2000.*ZTCORR)))
        ALPHA_V(JL,IND_IO)      =                                 &
          1.0_DP/(1.0_DP+1.0_DP/(1.0_DP/(1.0_DP/5.0E-1 -1.0_DP) * &
          EXP(2000.*ZTCORR)))
        ALPHA_V(JL,IND_I2O2)    =                                 &
          1.0_DP/(1.0_DP+1.0_DP/(1.0_DP/(1.0_DP/1.0E-1 -1.0_DP) * &
          EXP(2000.*ZTCORR)))
        ALPHA_V(JL,IND_ICL)     =                                 &
          1.0_DP/(1.0_DP+1.0_DP/(1.0_DP/(1.0_DP/1.0E-2 -1.0_DP) * &
          EXP(2000.*ZTCORR)))
        ALPHA_V(JL,IND_IBR)     =                                 &
          1.0_DP/(1.0_DP+1.0_DP/(1.0_DP/(1.0_DP/1.0E-2 -1.0_DP) * &
          EXP(2000.*ZTCORR)))
        ALPHA_V(JL,IND_INO2)    =                                 &
          1.0_DP/(1.0_DP+1.0_DP/(1.0_DP/(1.0_DP/1.0E-1 -1.0_DP) * &
          EXP(2000.*ZTCORR)))
        ALPHA_V(JL,IND_OIO)     = 0.01_DP
        ALPHA_V(JL,IND_HIO3)    = 0.01_DP

        ALPHA_V(JL,IND_SO2)     = 1.1E-01_DP
        ALPHA_V(JL,IND_H2SO4)   = 0.65_DP
        ALPHA_V(JL,IND_DMSO)    = &
          1.0_DP/(EXP(-5.12E3*CAL2J/(R_GAS*ZTEMP)+23.1*CAL2J/R_GAS)+1.0_DP)
        ALPHA_V(JL,IND_CH3SO3H) = &
          1.0_DP/(EXP(-3.50E3*CAL2J/(R_GAS*ZTEMP)+16.7*CAL2J/R_GAS)+1.0_DP)

        ALPHA_V(JL,IND_Hg)      = 0.001_DP
        ALPHA_V(JL,IND_HgO)     = 0.1_DP

      END DO

    END SUBROUTINE SCAV_ALPHA_V
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  SUBROUTINE SCAV_ALPHA_V_GEN

      ! CALCULATION OF ACCOMMODATION COEFFICIENTS NEEDED FOR K_EX CALCULATIONS
      USE MESSY_MAIN_CONSTANTS_MEM, ONLY: CAL2J, R_GAS, T0
      USE MESSY_SCAV_MEM,           ONLY: luse_empiric_alpha, alpha0, alpha_T
      USE MESSY_SCAV_L_KPP

      REAL(DP) :: ZTCORR, ZTEMP
      INTEGER  :: JL, IX
      
      INTRINSIC EXP

      ! INIT
      ! STANDARD VALUE
      ALPHA_V(:,:)      = 0.1_DP

      DO JL=1,NVEC
        ZTEMP= VTEMP(JL)

        ZTCORR=1./ZTEMP-1./T0

        ! FOR SOME SPECIES A TEMPERATURE DEPENDENCE IS CALCULATED/ESTIMATED
        ! USING:
        ! ALPHA = 1./(1.+1./(1./(1./ALPHA(T0)-1.)*EXP((-DELTAH/R_GAS)*ZTCORR)))
        DO ix = 1, LSPEC
           IF (alpha0(ix) > 0._dp) THEN
              IF (alpha_T(ix) > 0._dp) THEN
                 ALPHA_V(JL,ix) =&
                      1.0_dp/(1.0_DP+1.0/(1.0/(1.0/alpha0(ix) -1.0_DP)&
                         *EXP(alpha_T(ix)*ZTCORR)))
              ELSE
                 ALPHA_V(JL,ix)  =alpha0(ix)
              END IF
           END IF
        END DO


        IF (LUSE_EMPIRIC_ALPHA) THEN
        ! EXPERIMENTALLY DETERMINED VALUES (FROM MOCCA)
           ALPHA_V(JL,IND_H2O2)    = &
                1./(EXP(-26.E3/(R_GAS*ZTEMP)+107.8456/R_GAS) + 1.0_DP)
           ALPHA_V(JL,IND_HCOOH)   =                                 &
                1.0_DP/(EXP(-7.9E3*CAL2J/(R_GAS*ZTEMP)+34.9*CAL2J/R_GAS)+1.0_DP)
           ALPHA_V(JL,IND_CH3OOH)  =                                 &
                1.0_DP/(EXP(-6.5E3*CAL2J/(R_GAS*ZTEMP)+32.5*CAL2J/R_GAS) + 1.0_DP)
           ALPHA_V(JL,IND_CO2)     =                                 &
                1.0_DP/(1.0_DP+1.0_DP/(1.0_DP/(1.0_DP/1.0E-2 -1.0_DP) * &
                EXP(2000.*ZTCORR)))
           ALPHA_V(JL,IND_CH3O2)   =                                 &
                1.0_DP/(1.0_DP+1.0_DP/(1.0_DP/(1.0_DP/1.0E-2 -1.0_DP) * &
                EXP(2000.*ZTCORR)))

           ALPHA_V(JL,IND_HCL)     = &
                1.0_DP/(EXP(-3.072E3/ZTEMP + 1.283E1_DP)+1.0_DP)

           ALPHA_V(JL,IND_CL2)     = &
                1.0_DP/(EXP(-1.3E4*CAL2J/(R_GAS*ZTEMP)+50.*CAL2J/R_GAS)+1.0_DP)

           ALPHA_V(JL,IND_HBR)     = &
                1.0_DP/(EXP(-3.94E3/ZTEMP + 1.664E1_DP) + 1.0_DP)
           
           ALPHA_V(JL,IND_BR2)     = &
                1.0_DP/(EXP(-1.3E4*CAL2J/(R_GAS*ZTEMP)+50.*CAL2J/R_GAS)+1.0_DP)

           ALPHA_V(JL,IND_HI)      = &
                1.0_DP/(EXP(-4.13E3/ZTEMP + 1.715E1_DP)+1.0_DP)
           ALPHA_V(JL,IND_HOI)     = ALPHA_V(JL,IND_HOBR)
           ALPHA_V(JL,IND_DMSO)    = &
                1.0_DP/(EXP(-5.12E3*CAL2J/(R_GAS*ZTEMP)+23.1*CAL2J/R_GAS)+1.0_DP)
           ALPHA_V(JL,IND_CH3SO3H) = &
                1.0_DP/(EXP(-3.50E3*CAL2J/(R_GAS*ZTEMP)+16.7*CAL2J/R_GAS)+1.0_DP)
           
        END IF
      END DO

    END SUBROUTINE SCAV_ALPHA_V_GEN
!-------------------------------------------------------------------------------
    SUBROUTINE SCAV_VMEAN_V

      USE MESSY_SCAV_MEM,     ONLY: LSPEC_GAS, GAS_IDX, GAS_MW, KPP_L_IDX
      INTEGER  :: JL, JT, IDX
      REAL(dp) :: MW

      VMEAN_V(:,:) = -999.999_DP ! DUMMY VALUE

      DO JT=1,lspec_gas
        IDX = KPP_L_IDX(SET)%GAS_SPEC(JT,GAS_IDX)
        MW  = KPP_L_IDX(SET)%GAS_ATTR(JT,GAS_MW)
        DO JL=1,NVEC
          VMEAN_V(JL,IDX) = VMEAN_F(VTEMP(JL),MW)
        ENDDO
      ENDDO
    END SUBROUTINE SCAV_VMEAN_V

    ! ------------------------------------------------------------------------

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

!===============================================================================

      SUBROUTINE SET_EASY_SCAV_COEFFICIENTS


!    SETTING OF ESTIMATED SCAVENGING COEFFICIENTS FOR FAST, EASY, FIRST GUESS 
!    OF SCAVENGING (ESTIMATES ARE DUE TO HS0, BUT MY OWN OPINION):
!    VERY SOLUBLE SPECIES GET HIGH VALUES (ALMOST 90%)
!    MODERATELY SOLUBLE SPECIES GET VALUES ABOUT 50%
!    HARDLY SOLUBLE SPECIES GET VALUES AROUND 1 - 10%
!    INSOLUBLE SPECIES GET A VALUE OF 0
!    IT IS POSSIBLE THAT SOME VALUES (E.G. FOR SO2 MUST BE ADJUSTED)

        USE MESSY_SCAV_MEM
        USE MESSY_SCAV_L_KPP

        SCAV_VALUE(:)=0.

        SCAV_VALUE(IND_HNO3)    = 0.95_DP
        SCAV_VALUE(IND_H2O2)    = 0.97_DP
        SCAV_VALUE(IND_CH3OOH)  = 0.30_DP
        SCAV_VALUE(IND_HCHO)    = 0.70_DP
        SCAV_VALUE(IND_HCOOH)   = 0.85_DP
        SCAV_VALUE(IND_CH3CO2H) = 0.85_DP
        SCAV_VALUE(IND_O3)      = 0.02_DP
        SCAV_VALUE(IND_SO2)     = 0.25_DP
       !HS0 IS SMALL, BUT TAKING HETEROGENEOUS REACTION INTO ACCOUNT
        SCAV_VALUE(IND_N2O5)    = 0.90_DP 
        SCAV_VALUE(IND_PAN)     = 0.30_DP
        SCAV_VALUE(IND_OH)      = 0.25_DP
        SCAV_VALUE(IND_HONO)    = 0.40_DP
       !HS0 IS SMALL, AND CO2 IS A FIXED SPECIES IN OUR MODEL SETUP
        SCAV_VALUE(IND_CO2)     = 0.00_DP
        SCAV_VALUE(IND_NH3)     = 0.40_DP
        SCAV_VALUE(IND_HO2)     = 0.85_DP
        SCAV_VALUE(IND_NO3)     = 0.25_DP
        SCAV_VALUE(IND_NO2)     = 0.02_DP
        SCAV_VALUE(IND_HNO4)    = 0.95_DP
        SCAV_VALUE(IND_CH3OH)   = 0.33_DP
        SCAV_VALUE(IND_CH3CHO)  = 0.25_DP
        SCAV_VALUE(IND_CH3O2)   = 0.30_DP
        SCAV_VALUE(IND_C2H5O2)    = 0.30_DP
        SCAV_VALUE(IND_H2SO4)   = 0.9999_DP
        SCAV_VALUE(IND_NO)      = 0.02_DP
        SCAV_VALUE(IND_CH3COCH3)= 0.25_DP

        SCAV_VALUE(IND_Hg)      = 0.001_DP
        SCAV_VALUE(IND_HgO)     = 0.9999_DP
        SCAV_VALUE(IND_HgCl )   = 0.9999_DP
        SCAV_VALUE(IND_HgCl2)   = 0.9999_DP
        SCAV_VALUE(IND_HgBr)    = 0.9999_DP
        SCAV_VALUE(IND_HgBr2)   = 0.9999_DP
        SCAV_VALUE(IND_ClHgBr)  = 0.9999_DP
        SCAV_VALUE(IND_BrHgOBr) = 0.9999_DP
        SCAV_VALUE(IND_ClHgOBr) = 0.9999_DP

        RETURN

      END SUBROUTINE SET_EASY_SCAV_COEFFICIENTS

!----------------------------------------------------------------------------

      SUBROUTINE DEFINE_KS_FOR_EFFECTIVE_HENRY


        USE MESSY_SCAV_MEM
        USE MESSY_SCAV_L_KPP
        K1(:)   = 0._DP
        K2(:)   = 0._DP
        K1_T(:) = 0._DP
        K2_T(:) = 0._DP

        K1(IND_HNO3)    = 15._DP
        K1(IND_HO2)     = 1.6E-5_DP
        K1(IND_H2O2)    = 2.2E-12_DP
        K1(IND_NH3)     = 1.75E-5_DP
        K1(IND_HONO)    = 5.1E-4_DP
        K1(IND_HNO4)    = 1.E-5_DP 
        K1(IND_HCOOH)   = 1.8E-4_DP
        K1(IND_CH3CO2H) = 1.75E-5_DP
        K1(IND_HCL)     = 1.7E6_DP
        K1(IND_HOCL)    = 3.2E-8_DP
        K1(IND_HBR)     = 1.0E9_DP
        K1(IND_HOBR)    = 2.3E-9_DP
        K1(IND_SO2)     = 1.7E-2_DP
        K1(IND_H2SO4)   = 1.0E3_DP

        K1_T(IND_HNO3)    =  8700._DP
        K1_T(IND_H2O2)    = -3730._DP
        K1_T(IND_NH3)     = -450._DP
        K1_T(IND_HONO)    = -1260._DP
        K1_T(IND_CH3CO2H) = -46._DP
        K1_T(IND_HCL)     =  6896._DP
        K1_T(IND_HOBR)    = -3091._DP
        K1_T(IND_SO2)     =  2090._DP

        K2(IND_SO2)     = 6.0E-8_DP
        K2(IND_H2SO4)   = 1.2E-2_DP

        K2_T(IND_SO2)     =  1120._DP
        K2_T(IND_H2SO4)   =  2720._DP


      END SUBROUTINE DEFINE_KS_FOR_EFFECTIVE_HENRY
!-------------------------------------------------------------------------------
    END MODULE MESSY_SCAV_LIQ
