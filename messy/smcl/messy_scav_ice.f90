MODULE MESSY_SCAV_ICE


!  AUTHOR: HOLGER TOST, MAINZ, SEP 2007

  USE MESSY_SCAV_I_KPP,           ONLY: DP, ISPEC=>NSPEC
  USE MESSY_SCAV_MEM,             ONLY: ISPEC_GAS
  USE MESSY_SCAV_INP_KPP

  IMPLICIT NONE

  PRIVATE 

  SAVE

  PUBLIC :: ALLOC_CHEM_FIELDS, DEALLOC_CHEM_FIELDS
  PUBLIC :: ALLOC_SCAV_VALUES_I
  PUBLIC :: ICE_PHYSC_INIT
  PUBLIC :: TRANSFER_COEFFICIENT_ICE
  PUBLIC :: SET_EASY_SCAV_COEFFICIENTS_I
  PUBLIC :: EQUIL_DISSOCIATION_COEFF

  REAL(DP), PUBLIC               :: TEMP, CV_I
  REAL(DP), PUBLIC               :: K_EXF(0:ISPEC), K_EXB(0:ISPEC)
  
! VARIABLES WITH VARIABLE VECTOR LENGTH
  REAL(DP), ALLOCATABLE, DIMENSION(:,:), PUBLIC ::       &
            VK_EXF, VK_EXB
  REAL(DP), ALLOCATABLE, DIMENSION(:), PUBLIC ::         &
            VCV_I, VK_EXF_N2O5, VK_EXF_CLNO3, VK_EXF_BRNO3

  ! PHOTOLYSIS_VALUES
  REAL(DP), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: JX_VEC

 ! TEMPERATURE , PRESSURE, DENSITY
  REAL(DP), PUBLIC, ALLOCATABLE, DIMENSION(:) :: VTEMP, VPRESS, VRHO, VICERAD
  REAL(DP), PUBLIC, ALLOCATABLE, DIMENSION(:) :: VIWC, VIWC_EFF, V_SAD
  REAL(DP), ALLOCATABLE, DIMENSION(:,:)       :: ALPHA_V, VMEAN_V

  INTEGER :: NVEC, SET

CONTAINS

!===============================================================================
  SUBROUTINE ALLOC_CHEM_FIELDS(LPROMA)

!!$#ifdef SCAV_VECTOR
!!$    USE MESSY_SCAV_I_KPP, ONLY : C, RCONST, VAR, RAD, FIX, F,  &
!!$                                 NVAR,      NRAD,     NFIX,    &
!!$                                 NVARST,    NRADST,   NFIXST,  &
!!$                                 NVECT,     NREACT
!!$#endif
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: LPROMA

!!$#ifdef SCAV_VECTOR
!!$    NVECT = LPROMA
!!$#endif    
    NVEC  = LPROMA

    ALLOCATE(VTEMP(NVEC))
    ALLOCATE(VPRESS(NVEC))
    ALLOCATE(VRHO(NVEC))
    ALLOCATE(VIWC(NVEC))
    ALLOCATE(VIWC_EFF(NVEC))
    ALLOCATE(VICERAD(NVEC))
    ALLOCATE(VCV_I(NVEC))
    ALLOCATE(V_SAD(NVEC))
   
    ALLOCATE(VMEAN_V(NVEC,0:ISPEC))
    ALLOCATE(ALPHA_V(NVEC,0:ISPEC))

    ALLOCATE(VK_EXF(NVEC,0:ISPEC))
    ALLOCATE(VK_EXB(NVEC,0:ISPEC))

    ALLOCATE(JX_VEC(NVEC,1:IP_MAX))

!!$#ifdef SCAV_VECTOR
!!$    ALLOCATE (C(NVEC,ISPEC))
!!$    ALLOCATE (RCONST(NVEC,NREACT))
!!$    
!!$    VAR => C(:,1+NVARST:NVARST+NVAR)
!!$    RAD => C(:,1+NRADST:NRADST+NRAD)
!!$    FIX => C(:,1+NFIXST:NFIXST+NFIX)
!!$    F   => C(:,1+NFIXST:NFIXST+NFIX)
!!$#endif
  END SUBROUTINE ALLOC_CHEM_FIELDS
!-------------------------------------------------------------------
  SUBROUTINE DEALLOC_CHEM_FIELDS
    
!#ifdef SCAV_VECTOR
!    USE MESSY_SCAV_I_KPP, ONLY : C, RCONST
!#endif
    IMPLICIT NONE

    DEALLOCATE(VTEMP)
    DEALLOCATE(VPRESS)
    DEALLOCATE(VRHO)
    DEALLOCATE(VIWC)
    DEALLOCATE(VIWC_EFF)
    DEALLOCATE(VICERAD)
    DEALLOCATE(VCV_I)
    DEALLOCATE(V_SAD)

    DEALLOCATE(VMEAN_V)
    DEALLOCATE(ALPHA_V)

    DEALLOCATE(VK_EXF)
    DEALLOCATE(VK_EXB)

    DEALLOCATE(JX_VEC)

!#ifdef SCAV_VECTOR
!    DEALLOCATE(C)
!    DEALLOCATE(RCONST)
!#endif
  END SUBROUTINE DEALLOC_CHEM_FIELDS

!===============================================================================

  SUBROUTINE ICE_PHYSC_INIT(LON_IWC, LON_IWC2, LON_ICERAD,&
                            LON_TEMP, LON_PRESS, JS)

! INITIALIZATION OF THE BOX VALUES FOR TEMPERATURE, PRESSURE, LIQUID WATER, 
! CHEMICAL CONCENTRATION OF WATER, DROPLET DIAMETER
    USE MESSY_MAIN_CONSTANTS_MEM,    ONLY: AVO => N_A, M_AIR, R_GAS
    USE MESSY_SCAV_MEM,              ONLY: I_SPEC
    USE MESSY_SCAV_I_KPP,            ONLY: IND_H2O_I

    INTRINSIC :: REAL
    REAL(DP), INTENT(IN)  :: LON_IWC(NVEC) ,LON_IWC2(NVEC), LON_ICERAD(NVEC),&
                             LON_TEMP(NVEC), LON_PRESS(NVEC)
    INTEGER,  INTENT(IN)  :: JS
    INTEGER :: JL

    SET = JS
    DO JL=1,NVEC
      VIWC_EFF(JL) = LON_IWC(JL)
      VIWC(JL)     = LON_IWC2(JL)
      IF (VIWC(JL)*REAL(IND_H2O_I,DP) >  0._DP) &
        I_SPEC(IND_H2O_I,JL) = 55.5_DP * AVO * VIWC(JL) * 1.E-6_DP
      VICERAD(JL)  = LON_ICERAD(JL)
      VTEMP(JL)    = LON_TEMP(JL)
      VPRESS(JL)   = LON_PRESS(JL)
      VRHO(JL)     = VPRESS(JL) * M_AIR / (R_GAS * VTEMP(JL) * 1.E3_DP)
      V_SAD(JL)    = 2.E-4_DP * (VIWC_EFF(JL) * 1.E3_DP)**0.9_DP
    ENDDO

  END SUBROUTINE ICE_PHYSC_INIT

!===============================================================================

  SUBROUTINE TRANSFER_COEFFICIENT_ICE(PROC)

    USE MESSY_MAIN_CONSTANTS_MEM, ONLY: R_GAS
    USE MESSY_SCAV_MEM,           ONLY: K_ICE, H_ADS, GAS_IDX, KPP_I_IDX
    USE MESSY_SCAV_AER,           ONLY: LAMBDA, AIRVISC, REY_RAIN
    USE MESSY_SCAV_I_KPP

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: PROC     ! SWITCH FOR DESTINCTION BETWEEN CLOUD AND RAIN
                                    ! 1 = CLOUD, 2 = RAIN
    INTEGER  :: JT, K, JL, IDX

    REAL(DP) :: ZRDRAD(NVEC)
    REAL(DP) :: ZTX(NVEC)
    REAL(DP) :: VISC(NVEC), ZLAMBDA(NVEC)
    REAL(DP) :: REY_R
    REAL(DP) :: KC, P1, P3, FLUX
    REAL(DP) :: FL1(NVEC), FL2(NVEC)
    REAL(DP) :: P2(NVEC)
    REAL(DP) :: K_EQ(ISPEC,NVEC)
    INTEGER, PARAMETER  :: IS = 6     ! ITERATION STEP FOR KMT CALCULATION
    REAL(DP) :: ZKMT_I(IS,ISPEC,NVEC), KMT(NVEC)
    REAL(DP) :: RADFIELD(IS+1), DEL_NUM 
    REAL(DP) :: NUMBER(IS,NVEC)
    REAL(DP) :: NUMBER_HELP(NVEC)

    VK_EXF(:,:)     = 0.
    VK_EXB(:,:)     = 0.
    
    KMT(:) = 0._DP

    SELECT CASE(PROC)

    CASE(1)

      DO JL=1,NVEC
        ZLAMBDA(JL) = LAMBDA(VTEMP(JL),VPRESS(JL))
          ! RECALCULATION RADIUS FROM MM TO M
        ZRDRAD(JL)  = 1.E-3_DP * VICERAD(JL)
      ENDDO
      !for the time being keep ICE part as hard-coded as it is
      CALL SCAV_ALPHA_V
      CALL SCAV_VMEAN_V

      ZTX(1:NVEC) = 1._DP/VTEMP(1:NVEC) - 1._DP/298._DP
      K_EQ(:,:) = 0.

      DO JT=1,ISPEC_GAS
        IDX = KPP_I_IDX(SET)%gas_spec(jt,gas_idx)
        ! CALCULATE equilibrium constant
        IF (K_ICE(IDX) > 0._DP ) THEN
          DO JL=1,NVEC
            K_EQ(IDX,JL) = K_ICE(IDX) * EXP( H_ADS(IDX) * ZTX(JL) )    
          ENDDO
        ! calculate diffusion limitation FOLLOWING SCHWARTZ           
          DO JL=1,NVEC
            KMT(JL) = VMEAN_V(JL,IDX) /  (ZRDRAD(JL) * (ZRDRAD(JL) / &
              ZLAMBDA(JL) + 4._DP / (3._DP * ALPHA_V(JL,IDX) )))
          ENDDO
         
        !  SETTING OF EXCHANGE FUNCTION            
        !  Attention: for backward function the uptake isotherm is included as 
        !             a funtion in the equation file, since it must be updated 
        !             every integration substep !!!
          DO JL=1,NVEC
            VK_EXF(JL,IDX) = KMT(JL) * VIWC_EFF(JL) * 1.E-3_DP
            VK_EXB(JL,IDX) = VK_EXF(JL,IDX) /                     &
                           ( K_EQ(IDX,JL) * R_GAS * VTEMP(JL) )   
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
        ZRDRAD(JL)  = 1.E-3_DP * VICERAD(JL)
        ZTX(JL)     = 1._DP / VTEMP(JL) - 1._DP/298._DP
        VISC(JL)    = AIRVISC(VTEMP(JL))
      ENDDO
      CALL SCAV_VMEAN_V
      K_EQ(:,:) = 0.
        
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
      
      DO JT=1,ISPEC_GAS
        IDX = KPP_I_IDX(SET)%gas_spec(jt,gas_idx)
        ! CALCULATE equilibrium constant
        IF (K_ICE(IDX) > 0._DP ) THEN
          DO JL=1,NVEC
            K_EQ(IDX,JL) = K_ICE(IDX) * EXP( H_ADS(IDX) * ZTX(JL) )    
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
        DO JT = 1,ISPEC_GAS
          IDX = KPP_I_IDX(SET)%gas_spec(jt,gas_idx)
          IF (K_ICE(IDX) > 0._DP ) THEN
            DO JL = 1,NVEC
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
        
      DO JT = 1,ISPEC_GAS
        IDX = KPP_I_IDX(SET)%gas_spec(jt,gas_idx)
        IF (K_ICE(IDX) > 0._DP ) THEN
          KMT(:) = 0._DP
          DO K = 1,IS
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

            VK_EXF(JL,IDX) = KMT(JL) * VIWC_EFF(JL) * 1.E-3_DP
            VK_EXB(JL,IDX) = VK_EXF(JL,IDX) /                     &
                           ( K_EQ(IDX,JL) * R_GAS * VTEMP(JL) )   
          ENDDO
        ENDIF
      ENDDO

    END SELECT

    VK_EXF(1:NVEC,0) = 0._DP
    VK_EXB(1:NVEC,0) = 0._DP


  END SUBROUTINE TRANSFER_COEFFICIENT_ICE

!===============================================================================

  SUBROUTINE ALLOC_SCAV_VALUES_I

    USE MESSY_SCAV_MEM,            ONLY: K_ICE, H_ADS, HI
    USE MESSY_SCAV_I_KPP

    IMPLICIT NONE

!------------------------------------------------------------------------------

!   FOR ICE PHASE
    HI(:)    = 0._DP
    K_ICE(:) = 0._DP
    H_ADS(:) = 0._DP

! SIMPLIFIED FOR EQUIVALENT HENRY EQUILIBRIUM OVER ICE

    HI(IND_H2O2)      = 5.E-3_DP

    HI(IND_HNO3)      = 5.E-3_DP



! UPTAKE RATES FOR ICE UPTAKE

!!$        K_ICE(IND_OH)      = 
!!$        K_ICE(IND_HO2)     = 
!!$        K_ICE(IND_H2O2)    = 
!!$        K_ICE(IND_O2)      = 
!!$        
!!$        K_ICE(IND_NH3)     = 
!!$        K_ICE(IND_NO)      = 
!!$        K_ICE(IND_NO2)     = 
!!$        K_ICE(IND_NO3)     = 
!!$        K_ICE(IND_HONO)    = 
    K_ICE(IND_HNO3)    = 6.77E+4_DP
!!$        K_ICE(IND_N2O5)    = 
!!$        K_ICE(IND_HNO4)    = 
!!$
!!$        K_ICE(IND_CH4)     = 
!!$        K_ICE(IND_HCHO)    = 
!!$        K_ICE(IND_HCOOH)   = 
!!$        K_ICE(IND_CH3OH)   = 
!!$        K_ICE(IND_CH3OOH)  = 
!!$        K_ICE(IND_CH3O2)   = 
!!$        K_ICE(IND_CO)      = 
!!$        K_ICE(IND_CO2)     = 
!!$
!!$        K_ICE(IND_C2H6)    = 
!!$        K_ICE(IND_C2H4)    = 
!!$        K_ICE(IND_CH3CO2H) = 
!!$        K_ICE(IND_CH3CHO)  = 
!!$        K_ICE(IND_C2H5OOH)   = 
!!$        K_ICE(IND_C2H5O2)    = 
!!$        K_ICE(IND_CH3CO3H) = 
!!$        K_ICE(IND_PA)      = 
!!$        K_ICE(IND_PAN)     = 
!!$
!!$        K_ICE(IND_C3H6)    = 
!!$        K_ICE(IND_C3H8)    = 
!!$        K_ICE(IND_CH3COCH3)= 
!!$        K_ICE(IND_MGLYOX)     = 
!!$        K_ICE(IND_HYPERACET)  = 
!!$        K_ICE(IND_CH3COCH2O2) = 
!!$        K_ICE(IND_ACETOL)  = 
!!$        K_ICE(IND_LHOC3H6OOH) = 
!!$        K_ICE(IND_LHOC3H6O2)  = 
!!$        K_ICE(IND_NOA)     =
!!$        K_ICE(IND_MPAN)    = 
!!$        K_ICE(IND_IC3H7OOH)   = 
!!$        K_ICE(IND_IC3H7O2)    = 
!!$        K_ICE(IND_IC3H7NO3)  = 
!!$
!!$        K_ICE(IND_NC4H10)   = 
!!$        K_ICE(IND_LC4H9OOH) = 
!!$        K_ICE(IND_LC4H9O2)  = 
!!$        K_ICE(IND_MEK)     = 
!!$        K_ICE(IND_LMEKO2)  = 
!!$        K_ICE(IND_LMEKOOH) = 
!!$        K_ICE(IND_MVK)     = 
!!$        K_ICE(IND_MVKO2)   = 
!!$        K_ICE(IND_MVKOOH)  = 
!!$        K_ICE(IND_BIACET)  = 
!!$        K_ICE(IND_LC4H9NO3)    = 
!!$
!!$        K_ICE(IND_C5H8)        = 
!!$        K_ICE(IND_NISOPOOH)    = 
!!$        K_ICE(IND_NISOPO2)     = 
!!$        K_ICE(IND_ISOPBO2)     = 
!!$        K_ICE(IND_NISOOH)      = 
!!$        K_ICE(IND_ISOPBOOH)    = 
!!$        K_ICE(IND_ISOPBNO3)    = 
!!$        K_ICE(IND_ISOPACDOOH)  = 
!!$        K_ICE(IND_ISOPACDNO3)  = 
!!$        K_ICE(IND_ISOPACDO2)   = 
!!$        K_ICE(IND_NC4CHO)      = 
!!$        K_ICE(IND_NISO3)       = 
!!$        K_ICE(IND_C5PAN18)     = 
!!$
!!$        K_ICE(IND_HCL)     = 
!!$        K_ICE(IND_HOCL)    = 
!!$        K_ICE(IND_CLNO2)   = 
!!$        K_ICE(IND_CLNO3)   = 
!!$        K_ICE(IND_CL2)     = 
!!$        K_ICE(IND_CL)      = 
!!$        K_ICE(IND_CLO)     =
!!$        K_ICE(IND_OCLO)    =
!!$
!!$        K_ICE(IND_HBR)     = 
!!$        K_ICE(IND_HOBR)    = 
!!$        K_ICE(IND_BRNO2)   = 
!!$        K_ICE(IND_BRNO3)   = 
!!$        K_ICE(IND_BR2)     = 
!!$        K_ICE(IND_BRCL)    = 
!!$        K_ICE(IND_BR)      = 
!!$        K_ICE(IND_BRO)     =
!!$
!!$        K_ICE(IND_HI)      = 
!!$        K_ICE(IND_HOI)     = 
!!$        K_ICE(IND_I2)      = 
!!$        K_ICE(IND_ICL)     = 
!!$        K_ICE(IND_IBR)     = 
!!$        K_ICE(IND_I)       = 
!!$        K_ICE(IND_IO)      = 
!!$        K_ICE(IND_INO2)    =
!!$        K_ICE(IND_I2O2)    =
!!$        K_ICE(IND_INO3)    =
!!$        K_ICE(IND_CH3I)    =
!!$        K_ICE(IND_CH2I2)   =
!!$        K_ICE(IND_CH2CLI)  =
!!$        K_ICE(IND_C3H7I)   =
!!$        K_ICE(IND_OIO)     =
!!$        K_ICE(IND_HIO3)    =
!!$
!!$        K_ICE(IND_SO2)     = 
!!$        K_ICE(IND_H2SO4)   = 
!!$        K_ICE(IND_DMS)     = 
!!$        K_ICE(IND_DMSO)    = 
!!$        K_ICE(IND_CH3SO3H) = 
!!$        K_ICE(IND_CH3SO2)  =
!!$        K_ICE(IND_CH3SO3)  =

! TEMPERATURE DEPEDENCE H_ADS OF UPTAKE COEFFICIENTS

!!$        H_ADS(IND_O2)      = 
!!$        H_ADS(IND_O3)      = 
!!$
!!$        H_ADS(IND_OH)      = 
!!$        H_ADS(IND_HO2)     = 
!!$        H_ADS(IND_H2O2)    = 
!!$
!!$        H_ADS(IND_NH3)     = 
!!$        H_ADS(IND_NO)      = 
!!$        H_ADS(IND_NO2)     = 
!!$        H_ADS(IND_NO3)     = 
!!$        H_ADS(IND_HONO)    = 
    H_ADS(IND_HNO3)    = 44000._DP
!!$        H_ADS(IND_N2O5)    = 
!!$        H_ADS(IND_HNO4)    = 
!!$
!!$        H_ADS(IND_CH4)     = 
!!$        H_ADS(IND_HCHO)    = 
!!$        H_ADS(IND_HCOOH)   = 
!!$        H_ADS(IND_CH3OH)   = 
!!$        H_ADS(IND_CH3OOH)  = 
!!$        H_ADS(IND_CH3O2)   = 
!!$        H_ADS(IND_CO)      = 
!!$        H_ADS(IND_CO2)     = 
!!$
!!$        H_ADS(IND_C2H6)    = 
!!$        H_ADS(IND_C2H4)    = 
!!$        H_ADS(IND_CH3CO2H) = 
!!$        H_ADS(IND_CH3CHO)  = 
!!$        H_ADS(IND_C2H5OOH)   = 
!!$        H_ADS(IND_C2H5O2)    = 
!!$        H_ADS(IND_CH3CO3H) = 
!!$        H_ADS(IND_PA)      = 
!!$        H_ADS(IND_PAN)     = 
!!$
!!$        H_ADS(IND_C3H6)    = 
!!$        H_ADS(IND_C3H8)    = 
!!$        H_ADS(IND_CH3COCH3)= 
!!$        H_ADS(IND_MGLYOX)     = 
!!$        H_ADS(IND_HYPERACET)  = 
!!$        H_ADS(IND_CH3COCH2O2) = 
!!$        H_ADS(IND_ACETOL)  = 
!!$        H_ADS(IND_LHOC3H6OOH) = 
!!$        H_ADS(IND_LHOC3H6O2)  = 
!!$        H_ADS(IND_NOA)     =
!!$        H_ADS(IND_MPAN)    = 
!!$        H_ADS(IND_IC3H7OOH)   = 
!!$        H_ADS(IND_IC3H7O2)    = 
!!$        H_ADS(IND_IC3H7NO3)  = 
!!$
!!$        H_ADS(IND_NC4H10)   = 
!!$        H_ADS(IND_LC4H9OOH) = 
!!$        H_ADS(IND_LC4H9O2)  = 
!!$        H_ADS(IND_MEK)     = 
!!$        H_ADS(IND_LMEKO2)  = 
!!$        H_ADS(IND_LMEKOOH) = 
!!$        H_ADS(IND_MVK)     = 
!!$        H_ADS(IND_MVKO2)   = 
!!$        H_ADS(IND_MVKOOH)  = 
!!$        H_ADS(IND_BIACET)  = 
!!$        H_ADS(IND_LC4H9NO3)    = 
!!$
!!$        H_ADS(IND_C5H8)        = 
!!$        H_ADS(IND_NISOPOOH)    = 
!!$        H_ADS(IND_NISOPO2)     = 
!!$        H_ADS(IND_ISOPBO2)     = 
!!$        H_ADS(IND_NISOOH)      = 
!!$        H_ADS(IND_ISOPBOOH)    = 
!!$        H_ADS(IND_ISOPBNO3)    = 
!!$        H_ADS(IND_ISOPACDOOH)  = 
!!$        H_ADS(IND_ISOPACDNO3)  = 
!!$        H_ADS(IND_ISOPACDO2)   = 
!!$        H_ADS(IND_NC4CHO)      = 
!!$        H_ADS(IND_NISO3)       = 
!!$        H_ADS(IND_C5PAN18)     = 
!!$
!!$        H_ADS(IND_HCL)     = 
!!$        H_ADS(IND_HOCL)    = 
!!$        H_ADS(IND_CLNO2)   = 
!!$        H_ADS(IND_CLNO3)   = 
!!$        H_ADS(IND_CL2)     = 
!!$        H_ADS(IND_CL)      = 
!!$        H_ADS(IND_CLO)     =
!!$        H_ADS(IND_OCLO)    =
!!$
!!$        H_ADS(IND_HBR)     = 
!!$        H_ADS(IND_HOBR)    = 
!!$        H_ADS(IND_BRNO2)   = 
!!$        H_ADS(IND_BRNO3)   = 
!!$        H_ADS(IND_BR2)     = 
!!$        H_ADS(IND_BRCL)    = 
!!$        H_ADS(IND_BR)      = 
!!$        H_ADS(IND_BRO)     =
!!$
!!$        H_ADS(IND_HI)      = 
!!$        H_ADS(IND_HOI)     = 
!!$        H_ADS(IND_I2)      = 
!!$        H_ADS(IND_ICL)     = 
!!$        H_ADS(IND_IBR)     = 
!!$        H_ADS(IND_I)       = 
!!$        H_ADS(IND_IO)      = 
!!$        H_ADS(IND_INO2)    =
!!$        H_ADS(IND_I2O2)    =
!!$        H_ADS(IND_INO3)    =
!!$        H_ADS(IND_CH3I)    =
!!$        H_ADS(IND_CH2I2)   =
!!$        H_ADS(IND_CH2CLI)  =
!!$        H_ADS(IND_C3H7I)   =
!!$        H_ADS(IND_OIO)     =
!!$        H_ADS(IND_HIO3)    =
!!$
!!$        H_ADS(IND_SO2)     = 
!!$        H_ADS(IND_H2SO4)   = 
!!$        H_ADS(IND_DMS)     = 
!!$        H_ADS(IND_DMSO)    = 
!!$        H_ADS(IND_CH3SO3H) = 
!!$        H_ADS(IND_CH3SO2)  =
!!$        H_ADS(IND_CH3SO3)  = 
    RETURN

  END SUBROUTINE ALLOC_SCAV_VALUES_I

!===============================================================================
  SUBROUTINE SCAV_ALPHA_V

      ! CALCULATION OF ACCOMMODATION COEFFICIENTS NEEDED FOR K_EX CALCULATIONS
      USE MESSY_MAIN_CONSTANTS_MEM, ONLY: CAL2J, R_GAS, T0
      USE MESSY_SCAV_I_KPP

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

      END DO

    END SUBROUTINE SCAV_ALPHA_V
!-------------------------------------------------------------------------------
    SUBROUTINE SCAV_VMEAN_V

      USE MESSY_SCAV_MEM,           ONLY: GAS_MW, GAS_IDX, KPP_I_IDX
      USE MESSY_SCAV_LIQ,           ONLY: VMEAN_F
      INTEGER  :: JL, JT, IDX
      REAL(dp) :: MW

      VMEAN_V(:,:) = -999.999_DP ! DUMMY VALUE

      DO JT=1,ispec_gas
        IDX = KPP_I_IDX(SET)%GAS_SPEC(JT,GAS_IDX)
        MW  = KPP_I_IDX(SET)%GAS_ATTR(JT,GAS_MW)
        DO JL=1,NVEC
          VMEAN_V(JL,IDX) = VMEAN_F(VTEMP(JL),MW)
        ENDDO
      ENDDO
    END SUBROUTINE SCAV_VMEAN_V

!===============================================================================
    SUBROUTINE SET_EASY_SCAV_COEFFICIENTS_I


!    SETTING OF ESTIMATED SCAVENGING COEFFICIENTS FOR FAST, EASY, FIRST GUESS 
!    OF SCAVENGING (ESTIMATES ARE DUE TO HS0, BUT MY OWN OPINION):
!    VERY SOLUBLE SPECIES GET HIGH VALUES (ALMOST 90%)
!    MODERATELY SOLUBLE SPECIES GET VALUES ABOUT 50%
!    HARDLY SOLUBLE SPECIES GET VALUES AROUND 1 - 10%
!    INSOLUBLE SPECIES GET A VALUE OF 0
!    IT IS POSSIBLE THAT SOME VALUES (E.G. FOR SO2 MUST BE ADJUSTED)

        USE MESSY_SCAV_MEM,           ONLY: SCAV_VALUE_I
        USE MESSY_SCAV_I_KPP

        SCAV_VALUE_I(:)=0.

        SCAV_VALUE_I(IND_HNO3)    = 0.2_DP
!!$        SCAV_VALUE_I(IND_H2O2)    = 0.97_DP
!!$        SCAV_VALUE_I(IND_CH3OOH)  = 0.30_DP
!!$        SCAV_VALUE_I(IND_HCHO)    = 0.70_DP
!!$        SCAV_VALUE_I(IND_HCOOH)   = 0.85_DP
!!$        SCAV_VALUE_I(IND_CH3CO2H) = 0.85_DP
!!$        SCAV_VALUE_I(IND_O3)      = 0.02_DP
!!$        SCAV_VALUE_I(IND_SO2)     = 0.25_DP
!!$       !HS0 IS SMALL, BUT TAKING HETEROGENEOUS REACTION INTO ACCOUNT
!!$        SCAV_VALUE_I(IND_N2O5)    = 0.90_DP 
!!$        SCAV_VALUE_I(IND_PAN)     = 0.30_DP
!!$        SCAV_VALUE_I(IND_OH)      = 0.25_DP
!!$        SCAV_VALUE_I(IND_HONO)    = 0.40_DP
!!$       !HS0 IS SMALL, AND CO2 IS A FIXED SPECIES IN OUR MODEL SETUP
!!$        SCAV_VALUE_I(IND_CO2)     = 0.00_DP
!!$        SCAV_VALUE_I(IND_NH3)     = 0.40_DP
!!$        SCAV_VALUE_I(IND_HO2)     = 0.85_DP
!!$        SCAV_VALUE_I(IND_NO3)     = 0.25_DP
!!$        SCAV_VALUE_I(IND_NO2)     = 0.02_DP
!!$        SCAV_VALUE_I(IND_HNO4)    = 0.95_DP
!!$        SCAV_VALUE_I(IND_CH3OH)   = 0.33_DP
!!$        SCAV_VALUE_I(IND_CH3CHO)  = 0.25_DP
!!$        SCAV_VALUE_I(IND_CH3O2)   = 0.30_DP
!!$        SCAV_VALUE_I(IND_C2H5O2)    = 0.30_DP
!!$        SCAV_VALUE_I(IND_H2SO4)   = 0.9999_DP
!!$        SCAV_VALUE_I(IND_NO)      = 0.02_DP
!!$        SCAV_VALUE_I(IND_CH3COCH3)= 0.25_DP

        RETURN

      END SUBROUTINE SET_EASY_SCAV_COEFFICIENTS_I
!----------------------------------------------------------------------------
      SUBROUTINE EQUIL_DISSOCIATION_COEFF

      USE MESSY_MAIN_CONSTANTS_MEM,    ONLY: AVO => N_A

      INTEGER  :: JL
      
      VCV_I(:)      = 0.0_DP
      DO JL=1,NVEC
         IF (VIWC(JL).GT.TINY(0._DP)) THEN
            VCV_I(JL) = 1.E6_DP / VIWC(JL) / AVO
         ELSE
            VCV_I(JL) = 0.E0
         ENDIF
      ENDDO
      
    END SUBROUTINE EQUIL_DISSOCIATION_COEFF
!----------------------------------------------------------------------------

END MODULE MESSY_SCAV_ICE
