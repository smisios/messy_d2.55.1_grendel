MODULE MESSY_SCAV

!  AUTHOR:   HOLGER TOST, MPI CHEMIE, MAINZ
!            LAST MODIFIED 10.01.2005

  USE MESSY_SCAV_L_KPP,         ONLY: DP
  USE MESSY_SCAV_MEM,           ONLY: LSCAV,                                  &
                                      LSCAV_CV, LSCAV_LS,                     &
                                      LSCAV_GAS, LSCAV_AER,                   &
                                      LSCAV_NUC, LSCAV_IMP,                   &
                                      LSCAV_I, LSCAV_L,                       &
                                      LSCAV_EASY, ISCAV_EASY,                 &
                                      L_SCAV_EASY, I_SCAV_EASY,               &
                                      CPL_AEROSOL, I_EVAP, FRAC_RESNUM,       &
                                      ISCAV_RATE_HITMP,                       &
                                      ISCAV_RATE_LOTMP_HET,                   &
                                      ISCAV_RATE_LOTMP_HOM,                   &
                                      LOVERWRITE_HENRY, LOVERWRITE_ALPHA,     &
                                      LUSE_EMPIRIC_ALPHA
  USE MESSY_SCAV_AER,           ONLY: COEFF_PARA
  USE MESSY_SCAV_LIQ,           ONLY: USE_SCHWARTZ


  IMPLICIT NONE

  SAVE
  PRIVATE

  INTRINSIC :: MAX, SQRT, NULL
  CHARACTER(LEN=*), PUBLIC, PARAMETER :: MODSTR='scav'
  CHARACTER(LEN=*), PUBLIC, PARAMETER :: MODVER='2.4.0'

  ! MZ_HT_20020115 AND OUTPUT PARAMETERS

! APPROX. AND USED CONVECTIVE LWC
  REAL(DP) ,PUBLIC, POINTER :: LWC_CV(:,:,:)              => NULL() 
    
 ! GLOBAL NAMELIST PARAMETERS ('COUPLING' PARAMETERS)
   
  REAL(DP) ,PUBLIC, POINTER     :: UPDRAFT(:,:,:)          => NULL()
  REAL(DP) ,PUBLIC, POINTER     :: TRACER_CONV(:,:,:)      => NULL()
  REAL(DP) ,PUBLIC, POINTER     :: TRAC_FIELD(:,:,:,:)     => NULL()
  REAL(DP) ,PUBLIC, POINTER     :: KBOT(:,:)               => NULL()
  REAL(DP) ,PUBLIC, POINTER     :: KRAINTOP(:,:)           => NULL()
  REAL(DP) ,PUBLIC, POINTER     :: COL                     => NULL()

! STRING FOR GAS SPECIES OUT FLUXES
  CHARACTER(LEN=3000), PUBLIC :: OUT_STRING                  
! STRING FOR AEROSOL OUT FLUXES
  CHARACTER(LEN=3000), PUBLIC :: OUT_STRING_AER              

! CHANNEL OBJECT FOR CONVECTIVE PRECIPITATION
  REAL (DP), PUBLIC, POINTER :: PRECFLX_CV(:,:,:)   => NULL()   
!                    CONVECTIVE SNOWFALL
  REAL (DP), PUBLIC, POINTER :: SNOWFLX_CV(:,:,:)   => NULL()   
!                    CONVECTIVE CLOUD COVER
  REAL (DP), PUBLIC, POINTER :: CONV_COVER(:,:,:)   => NULL()   
!                    CONVECTIVE CLOUD BASE LEVEL
  REAL (DP), PUBLIC, POINTER :: KCONBOT(:,:)        => NULL()   
!                    FRESHLY FORMED PRECIPITATION
  REAL (DP), PUBLIC, POINTER :: PCVDPREC(:,:,:)     => NULL()   
!                    FRESHLY FORMED SNOW
  REAL (DP), PUBLIC, POINTER :: PCVDSNOW(:,:,:)     => NULL() 
!                    CONVECTIVE LWC
  REAL (DP), PUBLIC, POINTER :: PCVLWC(:,:,:)       => NULL() 
!                    CONVECTIVE RAIN FORMATION
  REAL (DP), PUBLIC, POINTER :: PCVRFORM(:,:,:)     => NULL()  
!                    CONVECTIVE IWC
  REAL (DP), PUBLIC, POINTER :: PCVIWC(:,:,:)       => NULL() 
!                    CONVECTIVE RAIN FORMATION
  REAL (DP), PUBLIC, POINTER :: PCVSFORM(:,:,:)     => NULL() 

!                FOR LARGE SCALE CLOUD COVER
  REAL (DP), PUBLIC, POINTER :: PCLCOVER(:,:,:)     => NULL()   
!                    LARGE SCALE PRECIPITATING COVER
  REAL (DP), PUBLIC, POINTER :: PRCOVER(:,:,:)      => NULL()   
!                    PRECIPITATION FORMATION RATE     [KG KG-1]
  REAL (DP), PUBLIC, POINTER :: PMRATEP(:,:,:)      => NULL()   
!                    PRECIPITATION FLUX               [KG M-2 S-1]
  REAL (DP), PUBLIC, POINTER :: PFPREC(:,:,:)       => NULL()   
!                    SNOW/ICE FORMATION RATE          [KG KG-1]
  REAL (DP), PUBLIC, POINTER :: PMRATESI(:,:,:)     => NULL()   
!                    SNOW/ICE PRECIPITATION FLUX      [KG M-2 S-1]
  REAL (DP), PUBLIC, POINTER :: PFSI(:,:,:)         => NULL()   
!                    LIQUID WATER CONTENT             [KG KG-1]
  REAL (DP), PUBLIC, POINTER :: PMLWC(:,:,:)        => NULL()   
!                    ICE/SNOW WATER CONTENT           [KG KG-1]
  REAL (DP), PUBLIC, POINTER :: PMIWC(:,:,:)        => NULL()   
!                    MELTING OF FROZEN PRECIP         [KG M-2 S-1]
  REAL (DP), PUBLIC, POINTER :: PIMELT(:,:,:)       => NULL()   
!                    ICE SEDIMENTATION FORMATION      [KG KG-1]
  REAL (DP), PUBLIC, POINTER :: PISEDI(:,:,:)       => NULL()   
  

! CHANNEL OBJECT FOR GRID MASS
  REAL (DP), PUBLIC, POINTER :: GRMASS(:,:,:)       => NULL()   
!                FOR GRID VOLUME
  REAL (DP), PUBLIC, POINTER :: GRVOL(:,:,:)        => NULL()   
!                FOR GRID BOX PRESSURE
  REAL (DP), PUBLIC, POINTER :: PRESS_3D(:,:,:)     => NULL()   
!                FOR INTERFACE PRESSURE
  REAL (DP), PUBLIC, POINTER :: PRESSI_3D(:,:,:)    => NULL()   

! CHANNEL OBJECT FOR ATTILA NUMBER OF BOXES PER GRID CELL
  REAL (DP), PUBLIC, POINTER :: NCB(:,:,:)          => NULL()   

  PUBLIC :: INIT_SCAV, INIT_SCAV_SPEC_STRUCT
  PUBLIC :: SCAV_EASY_LIQ, SCAV_EASY_ICE
  PUBLIC :: CALC_LWC_BC, CALC_LWC_BC_CV, CALC_LWC_IC_CV
  PUBLIC :: SCAV_AQCHEM_KPP, SCAV_ICECHEM_KPP
  PUBLIC :: CALC_LWC_BC_N

CONTAINS

!-------------------------------------------------------------------------------

  SUBROUTINE INIT_SCAV(IOU, FSTAT)

    ! SCAV MODULE ROUTINE (CORE)
    !
    ! READ SCAV NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! AUTHOR: PATRICK JOECKEL, MPICH, FEB 2002
    ! MODIFIED: LAURENS GANZEVELD, MPICH, 14-11-2002

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN)  :: IOU   ! LOGICAL I/O UNIT

    NAMELIST /CTRL/     LSCAV_LS,   LSCAV_CV,               &
                        LSCAV_NUC,  LSCAV_IMP,              & 
                        LSCAV_L,    LSCAV_I,                &
                        LSCAV_GAS,  USE_SCHWARTZ,           &
                        LSCAV_EASY, L_SCAV_EASY,            &
                        ISCAV_EASY, I_SCAV_EASY,            &
                        LSCAV_AER,  COEFF_PARA, CPL_AEROSOL,&
                        I_EVAP, FRAC_RESNUM,                &
                        ISCAV_RATE_HITMP,                   &
                        ISCAV_RATE_LOTMP_HET,               &
                        ISCAV_RATE_LOTMP_HOM,               &
                        LOVERWRITE_HENRY, LOVERWRITE_ALPHA, &
                        LUSE_EMPIRIC_ALPHA


    ! LOCAL
    LOGICAL                           :: LEX          ! FILE EXISTS ?
    INTEGER                           :: FSTAT        ! FILE STATUS

    ! INITIALIZE GLOBAL CONTROL VARIABLES

    LSCAV        = .FALSE.
    LSCAV_LS     = .FALSE.
    LSCAV_CV     = .FALSE.
    LSCAV_NUC    = .FALSE.
    LSCAV_IMP    = .FALSE.
    LSCAV_L      = .FALSE.
    LSCAV_I      = .FALSE.
    LSCAV_GAS    = .FALSE.
    LSCAV_AER    = .FALSE.
    USE_SCHWARTZ = .FALSE.
    LSCAV_EASY   = .FALSE.
    ISCAV_EASY   = .FALSE.
    L_SCAV_EASY  = 0
    I_SCAV_EASY  = 0
    COEFF_PARA   = 0
    CPL_AEROSOL  = 0
    I_EVAP       = 1
    FRAC_RESNUM  = 0.6_dp
    ISCAV_RATE_HITMP     = 0.1_dp                   
    ISCAV_RATE_LOTMP_HET = 0.1_dp
    ISCAV_RATE_LOTMP_HOM = 0.05_dp

    ! INPUT NAMELIST
    WRITE(*,*) '*******************************************************'
    WRITE(*,*) 'START SCAV MODULE INITIALISATION (INIT_SCAV)'
    WRITE(*,*) '*******************************************************'
    ! CHECK IF FILE EXISTS, YES: SWITCH SCAV ON
    !                       NO : KEEP SCAV SWITCHED OFF
    INQUIRE(FILE=MODSTR//'.nml', EXIST=LEX)
    IF (.NOT.LEX) THEN
       WRITE(*,*) 'WARNING *** FILE '//MODSTR//'.nml'//'  NOT FOUND !'
       WRITE(*,*) ' SCAV SWITCHED OFF !'
       WRITE(*,*) '******************************************************'
       RETURN
    END IF
    ! SET GLOBAL SWITCH
    LSCAV = .TRUE.

    ! READ NAMELIST
    OPEN(IOU,FILE=MODSTR//'.nml')
    WRITE(*,*) 'READING NAMELIST FROM '//MODSTR//'.nml', &
         ' (UNIT ',IOU,') ...'
    READ(IOU, NML=CTRL, IOSTAT=FSTAT)
    IF (FSTAT /= 0) THEN
       WRITE(*,*) 'ERROR *** READ ERROR IN NAMELIST ', &
            MODSTR//'.nml'//' !'
       WRITE(*,*) '******************************************************'
       RETURN
     END IF
    CLOSE(IOU)
    
    WRITE(*,*) '******************************************************'
    WRITE(*,*) 'END SCAV MODULE INITIALISATION (INIT_SCAV)'
    WRITE(*,*) '******************************************************'

  END SUBROUTINE INIT_SCAV

!==============================================================================
  
  SUBROUTINE INIT_SCAV_SPEC_STRUCT

    USE messy_scav_l_kpp,       ONLY: lspec => nspec
    USE messy_scav_i_kpp,       ONLY: ispec => nspec
    USE messy_main_tools,       ONLY: strcrack
    USE messy_scav_mem
    USE messy_scav_l_kpp,       ONLY: spc_l_names => spc_names    
    USE messy_scav_i_kpp,       ONLY: spc_i_names => spc_names

    IMPLICIT NONE

    INTEGER                        :: jt
    INTEGER                        :: count_g, count_l, dummy
    CHARACTER(LEN=26), POINTER, DIMENSION(:)     :: strname => NULL()
    
    ! liquid phase

    count_g = 0
    count_l = 0
    LSPEC_GAS = 0
    LSPEC_LIQ = 0
    DO jt=1,lspec
      if (associated(strname)) DEALLOCATE (strname)
      NULLIFY(strname)     
      call strcrack(spc_l_names(jt),'_', strname, dummy)
      IF (TRIM(strname(1)) == 'Prod') CYCLE
      IF (dummy==1) THEN 
        LSPEC_GAS = LSPEC_GAS + 1
      ELSE
        LSPEC_LIQ = LSPEC_LIQ + 1
      ENDIF
    ENDDO

    ALLOCATE(kpp0_l_idx%gas_spec(lspec_gas, gas_max))
    ALLOCATE(kpp0_l_idx%liq_spec(lspec_liq, liq_max))

    DO jt=1,lspec
      dummy=0
      if (associated(strname)) DEALLOCATE (strname)
      NULLIFY(strname)     
      call strcrack(spc_l_names(jt),'_', strname, dummy)
      IF (TRIM(strname(1)) == 'Prod') CYCLE
      IF (dummy==1) THEN 
        count_g = count_g + 1
        kpp0_l_idx%gas_spec(count_g,gas_idx) = jt
      ELSE
        count_l = count_l + 1
        kpp0_l_idx%liq_spec(count_l,liq_idx) = jt
        
      ENDIF
    ENDDO


    ! ice phase
    count_g = 0
    count_l = 0
    ISPEC_GAS = 0
    ISPEC_ICE = 0
    DO jt=1,ispec
      dummy=0
      if (associated(strname)) DEALLOCATE (strname)
      NULLIFY(strname)     
      call strcrack(spc_i_names(jt),'_', strname, dummy)
      IF (TRIM(strname(1)) == 'Prod') CYCLE
      IF (dummy==1) THEN 
        ISPEC_GAS = ISPEC_GAS + 1
      ELSE
        ISPEC_ICE = ISPEC_ICE + 1
      ENDIF
    ENDDO

    ALLOCATE(kpp0_i_idx%gas_spec(ispec_gas, gas_max))
    ALLOCATE(kpp0_i_idx%ice_spec(ispec_ice, ice_max))

    DO jt=1,ispec
      dummy=0
      if (associated(strname)) DEALLOCATE (strname)
      NULLIFY(strname)     
      call strcrack(spc_i_names(jt),'_', strname, dummy)
      IF (TRIM(strname(1)) == 'Prod') CYCLE
      IF (dummy==1) THEN 
        count_g = count_g + 1
        kpp0_i_idx%gas_spec(count_g,gas_idx) = jt
      ELSE
        count_l = count_l + 1
        kpp0_i_idx%ice_spec(count_l,ice_idx) = jt
        
      ENDIF
    ENDDO

  END SUBROUTINE INIT_SCAV_SPEC_STRUCT

!==============================================================================

!=============================== LIQUID EASY SCAVENGING =======================

  SUBROUTINE SCAV_EASY_LIQ(LWC, TEMP, LPROMA, JS)

    IMPLICIT NONE

    INTEGER,  INTENT(IN)    :: LPROMA, JS
    REAL(DP), INTENT(IN)    :: LWC(LPROMA), TEMP(LPROMA)

    SELECT CASE (L_SCAV_EASY)

    CASE(1)
      CALL SCAV_EASY_L1(LWC,TEMP,LPROMA,JS)
    CASE(2)
      CALL SCAV_EASY_L2(LWC,TEMP,LPROMA,JS)
    CASE(3)
      CALL SCAV_EASY_L3(LWC,TEMP,LPROMA,JS)
    CASE DEFAULT
      RETURN
    END SELECT
      
    RETURN
    
  END SUBROUTINE SCAV_EASY_LIQ
!----------------------------------------------------
  SUBROUTINE SCAV_EASY_L1(LWC,TEMP,LPROMA,JS)

! CALCULATES THE "EASY" SCAVENGING IN CLOUDS AND PRECIPITATION
! ACCORDING TO FIXED COEFFICIENTS

    USE MESSY_SCAV_MEM,         ONLY: SCAV_VALUE, HS0, &
                                      L_SPEC, LSPEC_GAS, KPP_L_IDX, GAS_IDX, &
                                      GAS_G2P

    INTEGER,  INTENT(IN)    :: LPROMA, JS
    REAL(DP), INTENT(IN)    :: LWC(LPROMA), TEMP(LPROMA)

    INTEGER                 :: JT, JL, IDX1, IDX2
      
    DO JT=1,LSPEC_GAS
       IDX1 = KPP_L_IDX(js)%gas_spec(jt,gas_idx)
       IDX2 = KPP_L_IDX(js)%gas_spec(jt,gas_g2p)
       IF (HS0(IDX1) > 0._DP ) THEN
          DO JL=1,LPROMA
             L_SPEC(IDX2,JL) = L_SPEC(IDX2,JL) +&
                SCAV_VALUE(IDX1) * L_SPEC(IDX1,JL)
             L_SPEC(IDX1,JL) = L_SPEC(IDX1,JL) * (1._DP-SCAV_VALUE(IDX1))
          ENDDO
       ENDIF
    ENDDO

    
    RETURN
  END SUBROUTINE SCAV_EASY_L1
!----------------------------------------------------
  SUBROUTINE SCAV_EASY_L2(LWC, TEMP, LPROMA, JS)

! CALCULATES THE "EASY" SCAVENGING IN PRECIPITATION
    USE MESSY_SCAV_MEM,         ONLY: HS0, DHT, L_SPEC, &
                                      LSPEC_GAS, KPP_L_IDX, GAS_IDX,     &
                                      GAS_G2P
    USE MESSY_SCAV_L_KPP,       ONLY: NSPEC

    INTEGER,  INTENT(IN)    :: LPROMA, JS
    REAL(DP), INTENT(IN)    :: LWC(LPROMA), TEMP(LPROMA)

    INTEGER                 :: JT, JL, IDX1, IDX2
    REAL(DP) :: HENRY(NSPEC, LPROMA)
    REAL(DP) :: ZTX(LPROMA), C_SUM
     


    ZTX(1:LPROMA) = 1._DP / TEMP(1:LPROMA) - 1._DP/298._DP
    HENRY = 0.
    DO JT=1,LSPEC_GAS
       IDX1 = KPP_L_IDX(js)%gas_spec(jt,gas_idx)
       IDX2 = KPP_L_IDX(js)%gas_spec(jt,gas_g2p)
      ! CALCULATE HENRY CONSTANT ! SEE SCAV_LIQ.F90 FOR COMMENTS
       IF (HS0(IDX1) > 0._DP ) THEN
          DO JL=1,LPROMA
             HENRY(IDX1,JL) = HS0(IDX1) * EXP(DHT(IDX1) * ZTX(JL))
             HENRY(IDX1,JL) = 1._DP / (HENRY(IDX1,JL) * 0.082 * TEMP(JL))
             HENRY(IDX1,JL) = HENRY(IDX1,JL)/(LWC(JL) * 1.E-3_DP)
 
             C_SUM = L_SPEC(IDX1,JL) + L_SPEC(IDX2,JL)
             L_SPEC(IDX2,JL) = C_SUM / (HENRY(IDX1,JL) + 1._DP)
             L_SPEC(IDX1,JL) = L_SPEC(IDX2,JL) * HENRY(IDX1,JL)
          ENDDO
       ENDIF
    ENDDO

    
    RETURN
  END SUBROUTINE SCAV_EASY_L2
!--------------------------------------------------------------
  SUBROUTINE SCAV_EASY_L3(LWC, TEMP, LPROMA, JS)

! CALCULATES THE "EASY" SCAVENGING IN PRECIPITATION
    USE MESSY_SCAV_MEM,         ONLY: HS0, DHT,  &
                                      L_SPEC, K1, K1_T, K2, K2_T, &
                                      KPP_l_IDX, GAS_IDX, GAS_G2P,&
                                      LSPEC_GAS
    USE MESSY_SCAV_L_KPP,       ONLY: NSPEC, IND_NH3

    INTEGER,  INTENT(IN)    :: LPROMA, JS
    REAL(DP), INTENT(IN)    :: LWC(LPROMA), TEMP(LPROMA)

    INTEGER                 :: JT, JL, IDX1, IDX2
    REAL(DP) :: HENRY(NSPEC,LPROMA), HENRY1(NSPEC,LPROMA), HENRY2(NSPEC)
    REAL(DP) :: ZTX(LPROMA), C_SUM
    REAL(DP) :: KS1, KS2
    REAL(DP), PARAMETER :: HPLUS = 1.E-5_DP   ! PH = 5 => 10^-5 MOL H+ / L

    ZTX(1:LPROMA) = 1._DP / TEMP(1:LPROMA) - 1._DP/298._DP
    HENRY = 0.

    DO JT=1,LSPEC_GAS
       IDX1 = KPP_L_IDX(js)%gas_spec(jt,gas_idx)
       IDX2 = KPP_L_IDX(js)%gas_spec(jt,gas_g2p)
       ! CALCULATE HENRY CONSTANT ! SEE SCAV_LIQ.F90 FOR COMMENTS
      IF (HS0(IDX1) > 0._DP ) THEN
        DO JL=1,LPROMA
          HENRY1(IDX1,JL) = HS0(IDX1) * EXP(DHT(IDX1)*ZTX(JL))
          KS1           = K1(IDX1)  * EXP(K1_T(IDX1)*ZTX(JL))
          KS2           = K2(IDX1)  * EXP(K2_T(IDX1)*ZTX(JL))
       ! EFFECTIVE HENRY'S LAW ACCORDING TO YIN, ACP, 2001
          HENRY2(IDX1) = HENRY1(IDX1,JL) * &
           ( 1._DP + KS1/HPLUS + (KS1*KS2)/(HPLUS*HPLUS) )
         
          HENRY(IDX1,JL) = 1._DP / (HENRY2(IDX1) * 0.082 * TEMP(JL))
          HENRY(IDX1,JL) = HENRY(IDX1,JL)/(LWC(JL)*1.E-3_DP)
        ENDDO
      ENDIF
    ENDDO
    DO JL=1,LPROMA
      KS1 = K1(IND_NH3) *  EXP(K1_T(IND_NH3)*ZTX(JL))
      HENRY1(IND_NH3,JL) = HS0(IND_NH3) * EXP(DHT(IND_NH3)*ZTX(JL))
      HENRY2(IND_NH3)    = HENRY1(IND_NH3,JL) * (1._DP + KS1 * HPLUS/1.E-14_DP)
      HENRY(IND_NH3,JL)  = 1._DP / (HENRY2(IND_NH3) * 0.082 * TEMP(JL))
      HENRY(IND_NH3,JL)  = HENRY(IND_NH3,JL)/(LWC(JL)*1.E-3_DP)
    ENDDO
    DO JT=1,LSPEC_GAS
       IDX1 = KPP_L_IDX(js)%gas_spec(jt,gas_idx)
       IDX2 = KPP_L_IDX(js)%gas_spec(jt,gas_g2p)
       IF (HS0(IDX1) > 0._DP ) THEN
          DO JL=1,LPROMA
             C_SUM = L_SPEC(IDX1,JL) + L_SPEC(IDX2,JL)
             L_SPEC(IDX2,JL) = C_SUM / (HENRY(IDX1,JL) + 1._DP)
             L_SPEC(IDX1,JL) = L_SPEC(IDX2,JL) * HENRY(IDX1,JL)
          ENDDO
       ENDIF
    ENDDO


    
    RETURN
  END SUBROUTINE SCAV_EASY_L3

!------------------------------------------------------------------------------


!=============================== ICE PHASE EASY SCAVENGING ====================
  SUBROUTINE SCAV_EASY_ICE(IWC, TEMP, PRESS,              &
                           phi_i, mju_i, iwc_i, iwc_T_i,  &
                           LPROMA, JS)

    IMPLICIT NONE

    INTEGER,  INTENT(IN)    :: LPROMA, JS
    REAL(DP), INTENT(IN)    :: IWC(LPROMA), TEMP(LPROMA)
    REAL(DP), INTENT(IN)    :: PRESS(LPROMA)
    REAL(DP), INTENT(INOUT) :: phi_i(lproma), mju_i(lproma)
    REAL(DP), INTENT(INOUT) :: iwc_i(lproma), iwc_T_i(lproma)

    SELECT CASE (I_SCAV_EASY)

    CASE(0)
      RETURN
    CASE(1)
      CALL SCAV_EASY_I1(IWC, TEMP, LPROMA, JS)
    CASE(2)
      CALL SCAV_EASY_I2(IWC, TEMP, LPROMA, JS)
    CASE(3)
       ! for output of partitioning parameters
      CALL SCAV_EASY_I3(IWC, TEMP, PRESS, phi_i, mju_i, &
                        iwc_i, iwc_T_i, LPROMA, js)
    CASE(4)
      CALL SCAV_EASY_I4(IWC, TEMP, PRESS, phi_i, mju_i, &
                        iwc_i, iwc_T_i, LPROMA, js)

    CASE DEFAULT
      RETURN
    END SELECT
      
    RETURN
    
  END SUBROUTINE SCAV_EASY_ICE
!----------------------------------------------------
  SUBROUTINE SCAV_EASY_I1(IWC,TEMP,LPROMA,JS)

! CALCULATES THE "EASY" SCAVENGING IN CLOUDS AND PRECIPITATION
! ACCORDING TO FIXED COEFFICIENTS

    USE MESSY_SCAV_MEM,         ONLY: I_SPEC, ISPEC_GAS,    &
                                      KPP_I_IDX, GAS_IDX,   &
                                      SCAV_VALUE_I, gas_g2p

    INTEGER,  INTENT(IN)    :: LPROMA, JS
    REAL(DP), INTENT(IN)    :: IWC(LPROMA), TEMP(LPROMA)

    INTEGER                 :: JT, JL, IDX, IDX2
      
    
    DO jt=1,ispec_gas
       IDX2 = KPP_I_IDX(js)%gas_spec(jt,gas_idx)
       IDX  = KPP_I_IDX(js)%gas_spec(jt,gas_g2p)
       IF (SCAV_VALUE_I(IDX2) > 0._DP ) THEN
          DO JL=1,LPROMA
             I_SPEC(IDX,JL) = I_SPEC(IDX,JL) +                     &
                              SCAV_VALUE_I(IDX2) * I_SPEC(IDX2,JL)
             I_SPEC(IDX2,JL) = I_SPEC(IDX2,JL) * (1._DP-SCAV_VALUE_I(IDX2))
          ENDDO
       END IF
    END DO

    RETURN

  END SUBROUTINE SCAV_EASY_I1
!----------------------------------------------------
  SUBROUTINE SCAV_EASY_I2(IWC,TEMP,LPROMA,JS)
    

! CALCULATES THE "EASY" SCAVENGING BY ICE/SNOW
    USE MESSY_SCAV_MEM,         ONLY: HS0, DHT, HI, &
                                      I_SPEC, ISPEC_GAS,             &
                                      KPP_I_IDX, GAS_IDX, GAS_I2L,   &
                                      GAS_G2P

    INTEGER,  INTENT(IN)    :: LPROMA, JS
    REAL(DP), INTENT(IN)    :: IWC(LPROMA), TEMP(LPROMA)

    INTEGER                 :: JT, JL, IDX1, IDX2, IDX3
    REAL(DP) :: HENRY(ISPEC_GAS, LPROMA)
    REAL(DP) :: ZTX(LPROMA), C_SUM
    REAL(DP) :: HSI(ISPEC_GAS)


    ZTX(1:LPROMA) = 1._DP / TEMP(1:LPROMA) - 1._DP/298._DP
    HENRY = 0.
    DO JT=1,ISPEC_GAS
       IDX1 = KPP_I_IDX(JS)%GAS_SPEC(JT,GAS_I2L)
       IDX2 = KPP_I_IDX(JS)%GAS_SPEC(JT,GAS_IDX)
       IDX3 = KPP_I_IDX(JS)%GAS_SPEC(JT,GAS_G2P)

       HSI(JT) = HS0(IDX1) * HI(IDX2) !* 1e3
       ! CALCULATE HENRY CONSTANT ! SEE SCAV_LIQ.F90 FOR COMMENTS
       IF (HSI(JT) > 0._DP ) THEN
          DO JL=1,LPROMA
             HENRY(JT,JL) = HSI(JT) * EXP(DHT(IDX1) * ZTX(JL))
             HENRY(JT,JL) = 1._DP / (HENRY(JT,JL) * 0.082 * TEMP(JL))
             HENRY(JT,JL) = HENRY(JT,JL) / (IWC(JL) * 1.E-3_DP)
          ENDDO
          DO JL=1,LPROMA
             C_SUM = I_SPEC(IDX2,JL) + I_SPEC(IDX3,JL)
             I_SPEC(IDX3,JL)   = C_SUM / (HENRY(JT,JL) + 1._DP)
             I_SPEC(IDX2,JL) = I_SPEC(IDX3,JL) * HENRY(JT,JL)
          ENDDO
       ENDIF
    ENDDO
    
    RETURN

  END SUBROUTINE SCAV_EASY_I2
  
!------------------------------------------------------------------------------
  SUBROUTINE SCAV_EASY_I3(IWC, TEMP, PRESS, phi_i, mju_i, &
                          iwc_i, iwc_T_i, LPROMA, JS)

    USE MESSY_MAIN_CONSTANTS_MEM, ONLY: ATM2PA, R_GAS, AVO => N_A, MN, MH, MO
    USE MESSY_SCAV_MEM,           ONLY: I_SPEC, ISPEC_GAS, KPP_I_IDX, GAS_IDX,   &
                                        H_ADS, K_ICE, GAS_G2P
    USE MESSY_SCAV_I_KPP,         ONLY: ISPEC => NSPEC, str_field_kpp_i => spc_names

    INTRINSIC               :: TRIM
     
    INTEGER,  INTENT(IN)    :: LPROMA, JS
    REAL(DP), INTENT(IN)    :: IWC(LPROMA), TEMP(LPROMA)
    REAL(DP), INTENT(IN)    :: PRESS(LPROMA)
    REAL(DP), INTENT(INOUT) :: phi_i(lproma), mju_i(lproma)   
    REAL(DP), INTENT(INOUT) :: iwc_i(lproma), iwc_T_i(lproma)
    
    INTEGER  :: JT, JL, IDX, ITER, IDX2
    REAL(DP) :: K_EQ(ISPEC_GAS,LPROMA)
    REAL(DP) :: P_GAS(ISPEC_GAS,LPROMA), P_TOT(ISPEC_GAS,LPROMA)
    REAL(DP) :: P_OLDICE(ISPEC_GAS,LPROMA), P_ICE(ISPEC_GAS,LPROMA)
    REAL(DP) :: ZTX(LPROMA), SAD(LPROMA), MC(LPROMA)
    REAL(DP) :: THETA, P_DELICE, ospec(ISPEC,lPROMA)
    REAL(DP) :: P_DIAG

    REAL(DP), PARAMETER :: TORR2PA = ATM2PA/760._DP
    REAL(DP), PARAMETER :: PA2TORR = 1._DP/TORR2PA
    REAL(DP), PARAMETER :: ALPHA_T = 0.27_DP !(SEE TABAZADEH, 99)
    REAL(DP), PARAMETER :: SIGMA_T = 1.E15_DP
    REAL(DP), PARAMETER :: BETA    = 9.656490576E18_DP

    REAL(DP), PARAMETER :: M_HNO3  = (MN + 3._dp *MO + MH)*0.001_dp
                                   ! 0.06301284_dp ![kg/mol] HNO3 molar weight

    ospec(1:ispec,1:lproma) = i_spec(1:ispec,1:lproma)
    ztx(1:lproma) = 1._dp/r_gas * (1._dp / temp(1:lproma) - 1._dp/228._dp)
    mc(1:lproma)  = (avo/1.e6_dp) * press(1:lproma) / (r_gas * temp(1:lproma))

!   determine equilibrium constants for all species and boxes  
!   convert concentrations to partial pressures in torr
    DO jt=1,ispec_gas
      idx = kpp_i_idx(js)%gas_spec(jt,gas_idx)
      idx2 = kpp_i_idx(js)%gas_spec(jt,gas_g2p)
      DO jl=1,lproma
        k_eq(jt,jl)  = SQRT(228._dp/temp(jl)) * k_ice(idx) * &
                       EXP(h_ads(idx) * ztx(jl))
        p_gas(jt,jl) = MAX(i_spec(idx,jl) * press(jl) / mc(jl) * pa2torr, 0._dp)
      ENDDO
      DO jl=1,lproma
         p_gas(jt,jl) = p_gas(jt,jl) + &
              i_spec(idx2,jl) * press(jl) / mc(jl) * pa2torr
      ENDDO
      DO jl=1,lproma
        p_tot(jt,jl) = p_gas(jt,jl)
        p_ice(jt,jl) = 0._dp
      ENDDO
    ENDDO

!   determine surface area of the ice
!   eqn. 2 of Kuhlmann & Lawrence, 2006; iwc in [kg/m3], sad in [g/m3]
    DO jl=1,lproma
      sad(jl) = 2.e-4_dp * (iwc(jl) * 1.e3_dp)**0.9_dp
    END DO
  
    DO jt=1,ispec_gas
      DO jl=1,lproma
        DO iter=1,20
          theta = alpha_t * SQRT( k_eq(jt,jl) * p_gas(jt,jl) ) / &
                  (1._dp  + SQRT( k_eq(jt,jl) * p_gas(jt,jl) ) )
          p_oldice(jt,jl) = p_ice(jt,jl)
          p_ice(jt,jl)    = MIN(theta * sigma_t * sad(jl) * &
                                temp(jl) / beta, p_tot(jt,jl) )
          IF (p_ice(jt,jl) < 0._dp) THEN
            p_ice(jt,jl) = 0._dp
          ENDIF
          p_delice     = p_ice(jt,jl) - p_oldice(jt,jl)
          p_gas(jt,jl) = p_gas(jt,jl) - p_delice
          IF (p_gas(jt,jl) <= 0._dp) THEN
            p_gas(jt,jl) = 0._dp
            EXIT
          ENDIF
        ENDDO
!        if (jt==1) print*, "in iteration",jt,jl,&
!          p_gas(jt,jl), p_ice(jt,jl), p_oldice(jt,jl), p_tot(jt,jl)
      ENDDO
    ENDDO

    DO jt=1,ispec_gas
       idx  = kpp_i_idx(js)%gas_spec(jt,gas_idx)
       idx2 = kpp_i_idx(js)%gas_spec(jt,gas_g2p)
       DO jl=1,lproma
          i_spec(idx,jl) = p_gas(jt,jl) * torr2pa * mc(jl) / press(jl)
          i_spec(idx2,jl) = p_ice(jt,jl) * torr2pa * mc(jl) / press(jl)
       ENDDO
    ENDDO

    ! for output of HNO3 partitioning parameters
    DO jl=1,lproma
      iwc_i(jl) = IWC(jl)
      iwc_T_i(jl) = 1.E-6_dp*EXP((temp(jl)/15._dp)-13.33_dp) 
    ENDDO
    DO jt=1,ispec_gas
      idx = kpp_i_idx(js)%gas_spec(jt,gas_idx)  !gas_idx: gas phase part
      IF (TRIM(str_field_kpp_i(idx)) == 'HNO3') THEN
        DO jl=1,lproma
           !  catch pathological case
          p_diag=p_ice(jt,jl)+p_gas(jt,jl)
          IF (p_diag <= 0._dp) THEN
            phi_i(jl) = -1._dp * p_ice(jt,jl)
          ELSE  
            phi_i(jl) = p_ice(jt,jl)/p_diag  !fraction of HNO3 molecules in the ice phase
          ENDIF
          mju_i(jl) = 1.e6_dp * M_HNO3 / avo * &
                      p_ice(jt,jl) * torr2pa * mc(jl) / press(jl) / &
                      iwc(jl)   ![HNO3_ice]/[H2O_ice]
        ENDDO
      ENDIF
    ENDDO

  END SUBROUTINE SCAV_EASY_I3
!------------------------------------------------------------------------------
! trapping of HNO3, according to Kaercher et al. (2009)
! * possible for HNO3 only (hardcoded)
  SUBROUTINE SCAV_EASY_I4(IWC, TEMP, PRESS, phi_i, mju_i, &
                          iwc_i, iwc_T_i, LPROMA, JS)

    USE MESSY_MAIN_CONSTANTS_MEM, ONLY: PI, K_B
    USE MESSY_SCAV_MEM,           ONLY: I_SPEC, ISPEC_GAS, KPP_I_IDX, GAS_IDX
    USE MESSY_SCAV_I_KPP,         ONLY: ISPEC => NSPEC, str_field_kpp_i => spc_names

    INTRINSIC               :: TRIM

    INTEGER,  INTENT(IN)    :: LPROMA, JS
    REAL(DP), INTENT(IN)    :: IWC(LPROMA), TEMP(LPROMA)
    REAL(DP), INTENT(IN)    :: PRESS(LPROMA)
    REAL(DP), INTENT(INOUT) :: phi_i(lproma), mju_i(lproma)
    REAL(DP), INTENT(INOUT) :: iwc_i(lproma), iwc_T_i(lproma)

    INTEGER                 :: JT, JL, IDX, ITER
    REAL(DP) :: P_GAS(LPROMA), P_TOT(LPROMA), P_ICE(LPROMA), P_OLDICE(LPROMA)
    REAL(DP) :: P_DELICE
    REAL(DP) :: e_ice(lproma), n_ice(lproma), s_ice(lproma)
    REAL(DP) :: Dw(lproma), D(lproma), uw(lproma), a(lproma), lw(lproma)
    REAL(DP) :: beta_w(lproma), a_dot(lproma), u(lproma), v(lproma)
    REAL(DP) :: nw(lproma), theta(lproma), temp1(lproma), temp2(lproma)
    REAL(DP) :: n_minus(lproma), lambda(lproma), kappa(lproma), beta(lproma)
    REAL(DP) :: aa_dot  !auxiliary for avoiding division by zero

!   constants
    REAL(DP), PARAMETER :: u_C12 = 1.660538782e-27_dp  ![kg] atomic mass unit
    REAL(DP), PARAMETER :: mw = 18.01528_dp*u_C12      ![kg] molecular weight of H2O
    REAL(DP), PARAMETER :: m = 63.01284_dp*u_C12       ![kg] molecular weight of HNO3
    REAL(DP), PARAMETER :: p0 = 101325._dp             ![Pa] reference pressure
    REAL(DP), PARAMETER :: T0 = 273.15_dp              ![K] reference temperature
    REAL(DP), PARAMETER :: vw = 3.E-29_dp              ![m^3] volume of H2O molecules in ice
    REAL(DP), PARAMETER :: alpha_w = 0.5_dp            !deposition coefficient of H2O on ice
    REAL(DP), PARAMETER :: alpha = 0.3_dp              !HNO3 deposition coefficient on ice
    REAL(DP), PARAMETER :: asp = 0.8_dp                !aspect ratio of hexagonal ice crystals
    REAL(DP), PARAMETER :: n_max = 1.8E23_dp           ![#/m^3] max. trapped number density
    REAL(DP), PARAMETER :: sigma=3.7e-15_dp            !surface area per adsorption site (mystic number)

!   parameters depending on 3d fields (temp, press) 
    DO jl=1,lproma
      e_ice(jl) = EXP(28.868_dp-(6132.9_dp/temp(jl)))  ![Pa] saturation vapour pressure over ice 
      !check also messy_main_tools.f90: tlucua (saturation mixing ratio, temp-Abh. in centi K)
      n_ice(jl) = e_ice(jl)/(k_B*temp(jl))              ![1/m^3] ice particle number density
      !iwc(jl) = 1.E-6*EXP((temp(jl)/15.)-13.33)       ![kg/m^3] ice water content used by BK09
!     iwc in cloud is in [kg/kg], but here it is in [kg/m^3], see messy_scav_e5.f90:scav_cv_new:zriwc & zrhoa
      s_ice(jl) = iwc(jl)/(mw*n_ice(jl))               ![] net (cirrus particle lifetime) ice supersaturation
      Dw(jl) = 2.11e-5_dp*((temp(jl)/T0)**1.94_dp)* &
               (p0/press(jl))                          ![m^2/s]  diffusivity of water molecules in air
      D(jl) = Dw(jl) * SQRT(mw/m)                      ![m^2/s]  diffusivity of HNO3 in air
      uw(jl) = SQRT((8._dp*k_B*temp(jl))/(mw*pi))      ![m/s] mean thermal speed
      a(jl) = 1.E-6_dp*(1.2_dp*temp(jl)-221.6_dp)* & 
              (((9._dp*SQRT(3._dp)) / & 
              (32_dp*pi*asp))**(1._dp/3._dp))          ![m] ice particle radius
      lw(jl) = 4._dp*Dw(jl)/uw(jl)
      beta_w(jl) = 1._dp/(a(jl)+(lw(jl)/alpha_w))
      a_dot(jl) = vw*Dw(jl)*beta_w(jl)*e_ice(jl)* &
                  s_ice(jl) / (k_B*temp(jl))            ![m/s] ice particle growth rate
      u(jl) = 18.32_dp*SQRT(temp(jl))                  ![m/s] mean thermal speed of HNO3
      v(jl)=(alpha*u(jl)*0.25_dp)/ &
            (n_max*sigma*7.5e-11_dp* &                 !7.5E-11 ... mystic number
            EXP(4585._dp/temp(jl)))                    ![m/s] escape speed
      lambda(jl) = 4._dp*v(jl)/(alpha*u(jl))
      beta(jl) = 1._dp/(1._dp+(4._dp*D(jl)/ &
             (alpha*u(jl)*a(jl))))
      nw(jl) = iwc(jl)/mw                              ![#/m^3] number density of water molecules in ice particles 
      aa_dot=a(jl)*a_dot(jl)                           !auxiliary var to catch pathological case  
      IF (aa_dot .GT. 0._dp) THEN
        kappa(jl) = D(jl)/(a(jl)*a_dot(jl))            !ratio of gas phase diffusion speed (D/a) and ice particle growth rate
      ELSE
        kappa(jl)=-1._dp*beta(jl)  !flag invalid kappa: makes temp2=0
      ENDIF

    ENDDO

!   initial HNO3 number densities: i_spec [#/cm^3], p_* [#/m^3] 
    DO jt=1,ispec_gas
      idx = kpp_i_idx(js)%gas_spec(jt,gas_idx)
      IF (TRIM(str_field_kpp_i(idx)) == 'HNO3') THEN
        DO jl=1,lproma
          !n_inf = NA*p_HNO3(ip)*(1.-(0.378*xv))/(Rd*T*Mair) ... [#/m^3] ambient HNO3 number density
          p_gas(jl) = MAX(1.E6_dp*i_spec(idx,jl), 0._dp)  !this is n_inf of BK09 
          p_tot(jl) = p_gas(jl)
          p_ice(jl) = 0._dp
        ENDDO
      ENDIF
    ENDDO

!   parameters depending also on HNO3 partial pressure -> part of iteration
    DO jl=1,lproma 
!     Within 195K < T < 240K and 0.1e-8mb < p_tot < 50.e-8mb (range of fig. 5 in BK09) no more
!     than 7 iterations are needed to get p_delice/p_ice < 1.e-9 (tested offline: trapping_SCAV.pro).     
      DO iter=1,15
        theta(jl) = p_gas(jl)/n_max
        temp1(jl) = (1._dp+beta(jl)*kappa(jl)*(lambda(jl)+theta(jl))) 
        temp2(jl) = (4._dp*beta(jl)*theta(jl)*(1._dp+beta(jl)*kappa(jl)))/ &
                    ((temp1(jl))**2._dp)  !the flag kappa=-beta means kappa=inf, which translates to temp2=0
        IF (temp2(jl) .LT. 1.e5_dp) THEN  !arbitrary treshold from trap_field_Kaercher.doc; trapping_SCAV.pro: MIN(temp2)=0.00095
          !n_minus(jl) = (1._dp+beta(jl)*kappa(jl))/temp1(jl)  !trap_field_Kaercher.doc
          n_minus(jl) = 1._dp/(lambda(jl)+theta(jl))   !limit of above eqn for kappa -> inf
        ELSE  !catch invalid root within eqn with biggest possible value (=1)
          n_minus(jl) = temp1(jl)/(2._dp*beta(jl)*theta(jl)) * (1._dp-SQRT(1._dp-MIN(temp2(jl),1._dp)))  !eqn. 9a in Kaercher et al.,2009
        ENDIF
        mju_i(jl) = vw * p_gas(jl) * n_minus(jl) !HNO3_i/H2O_i
        phi_i(jl) = (nw(jl)*vw*n_minus(jl))/(1._dp + nw(jl)*vw*n_minus(jl)) !fraction of HNO3 molecules in the ice phase
        p_oldice(jl) = p_ice(jl)
        p_ice(jl) = MIN(p_tot(jl)*phi_i(jl), p_tot(jl))
        IF (p_ice(jl) < 0._dp) THEN
          p_ice(jl) = 0._dp
        ENDIF
        p_delice = p_ice(jl) - p_oldice(jl)
        p_gas(jl) = p_gas(jl) - p_delice
        IF (p_gas(jl) <= 0._dp) THEN
          p_gas(jl) = 0._dp
          EXIT
        ENDIF      
      ENDDO 
!      print*, "in iteration",jl,&
!          p_gas(jl), p_ice(jl), p_oldice(jl), p_tot(jl)
    ENDDO

!   conversion back to [#/cm3]
    DO jt=1,ispec_gas
      idx = kpp_i_idx(js)%gas_spec(jt,gas_idx)  !gas_idx: gas phase part
      IF (TRIM(str_field_kpp_i(idx)) == 'HNO3') THEN
        DO jl=1,lproma
          i_spec(idx,jl) = p_gas(jl)*1.E-6_dp   ![#/cm^3] 
        ENDDO
      ENDIF
    ENDDO
    DO jt = 1,ispec
      IF (TRIM(str_field_kpp_i(jt)) == 'HNO3_i') THEN
        DO jl=1,lproma
          i_spec(jt,jl) = p_ice(jl)*1.E-6_dp   ![#/cm^3]
        ENDDO
      ENDIF
    ENDDO
    DO jl=1,lproma
      iwc_i(jl) = IWC(jl)
      iwc_T_i(jl) = 1.E-6_dp*EXP((temp(jl)/15._dp)-13.33_dp)
    ENDDO

  END SUBROUTINE SCAV_EASY_I4
!op_kg_20110126-
!------------------------------------------------------------------------------
!================================ END EASY SCAVENGING =========================
  
!=============================== KPP LIQUID SCAVENGING ========================

  SUBROUTINE SCAV_AQCHEM_KPP (PTEMP, PPRESS, AQRLWC_EFF, AQRLWC,   &
                              AQRDRAD, TIME_STEP_LEN, PROC, LPROMA, &
                              LIDX, JK, JROW, JS, NSTEPS, RSTEPS)

    ! NEW SUBROUTINE IN WHICH THE AQUEOUS CHEMISTRY IS BEING
    ! CALCULATED PRODUCED BY KPP 
    USE MESSY_MAIN_CONSTANTS_MEM, ONLY: R_GAS, N_A

    USE MESSY_SCAV_MEM,           ONLY: L_SPEC
    USE MESSY_SCAV_LIQ
! KPP VARIABLES
    USE MESSY_SCAV_L_KPP
! KPP Input variables

! THIS CONTAINS THE PREPROCESSOR DIRECTIVE IF A VECTORIZED SOLVER IS USED
!#include "messy_scav_l_vec.inc"
!#DEFINE SCAV_VECTOR
    IMPLICIT NONE
    SAVE
    
    INTEGER,  INTENT(IN) :: PROC                ! 1 = CLOUD, 2 = RAIN
    INTEGER,  INTENT(IN) :: LPROMA, JK, JROW, LIDX(LPROMA), JS
    REAL(DP), INTENT(IN) :: TIME_STEP_LEN
    REAL(DP), INTENT(IN) :: PTEMP(LPROMA), AQRLWC(LPROMA), AQRLWC_EFF(LPROMA), &
                            AQRDRAD(LPROMA), PPRESS(LPROMA)

    REAL(dp), INTENT(INOUT) :: NSTEPS(lproma), RSTEPS(lproma)
    INTEGER  :: JL, status

    REAL(dp) :: cair(lproma)
    REAL(dp) :: conc(lproma,1:nspec)

    INTEGER  :: kpp_stat(lproma), kpp_steps(lproma), kpp_rsteps(lproma)

    
    CALL ALLOC_CHEM_FIELDS(LPROMA)
! SETTING UP PHYSICAL PARAMETERS FOR EACH BOX
    CALL LIQ_PHYSC_INIT(AQRLWC_EFF, AQRLWC , AQRDRAD, PTEMP, PPRESS, JS)

    CALL SET_PHOTO(LIDX, JK, JROW)
! CALCULATES THE EQUILIBRIUM COEFFICIENTS / TESTFACS
    CALL EQUIL_DISSOCIATION_COEFF
! CALCULATES THE TRANSFER COEFFICIENTS
    CALL TRANSFER_COEFFICIENT_LIQ(PROC)

    conc(:,:) = 0._dp
    cair(:) = (N_A/1.E6_dp) * ppress(:) / (R_gas*ptemp(:))
    IF (.NOT.ivec) THEN       ! solver is not created with kp4 and uses a non-vector format
       DO jl=1,lproma

          CALL fill_jx (status, jx_vec(jl:jl,:))
          IF (status /= 0) STOP "fill_jx array size"
          
          CALL fill_temp (status,ptemp(jl:jl))
          IF (status /= 0) STOP "fill_temp array size"
          
          CALL fill_press (status, ppress(jl:jl))
          IF (status /= 0) STOP "fill_press array size"
          
          CALL fill_cair (status, cair(jl:jl))
          IF (status /= 0) STOP "fill_cair array size"
          
          CALL fill_lwc (status, vliq_lwc(jl:jl))
          IF (status /= 0) STOP "fill_lwc array size"
          
          CALL fill_cv_l (status, vcv_l(jl:jl))
          IF (status /= 0) STOP "fill_cv_l"
          
          CALL fill_k_exf(status, vk_exf(jl:jl,1:nspec))
          IF (status /= 0) STOP "fill_k_exf array size"
          CALL fill_k_exb(status, vk_exb(jl:jl,1:nspec))
          IF (status /= 0) STOP "fill_k_exb array size"
          CALL fill_k_exf_N2O5(status, vk_exf_N2O5(jl:jl))
          IF (status /= 0) STOP "fill_k_exf_N2O5 array size"
          CALL fill_k_exf_ClNO3(status, vk_exf_ClNO3(jl:jl))
          IF (status /= 0) STOP "fill_k_exf_ClNO3 array size"
          CALL fill_k_exf_BrNO3(status, vk_exf_BrNO3(jl:jl))
          IF (status /= 0) STOP "fill_k_exf_BrNO3 array size"

          conc(1,1:NSPEC) = L_SPEC(1:NSPEC,JL)
    ! KPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPP
          CALL kpp_integrate(TIME_STEP_LEN, conc, ierrf=kpp_stat, xNacc=kpp_steps, xNrej=kpp_rsteps)
    ! KPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPP
          conc(:,:) = MAX(conc(:,:), 0.0_DP) ! force positive results
          L_SPEC(1:NSPEC,JL) = conc(1,1:NSPEC)

          nsteps(jl) = kpp_steps(1)
          rsteps(jl) = kpp_rsteps(1)
          if (kpp_stat(1) /= 1) print*, "Error in SCAV_L_KPP: ", kpp_stat(1)
       END DO

    ELSE  ! solver uses a vectorised format

       CALL fill_jx (status, jx_vec)
       IF (status /= 0) STOP "fill_jx array size"

       CALL fill_temp (status,ptemp(1:LPROMA))
       IF (status /= 0) STOP "fill_temp array size"
       
       CALL fill_press (status, ppress)
       IF (status /= 0) STOP "fill_press array size"
       
       CALL fill_cair (status, cair)
       IF (status /= 0) STOP "fill_cair array size"
       
       CALL fill_lwc (status, vliq_LWC)
       IF (status /= 0) STOP "fill_lwc array size"
       
       CALL fill_cv_l (status, VCV_L)
       IF (status /= 0) STOP "fill_cv_l"
       
       CALL fill_k_exf(status, vk_exf(:,1:NSPEC))
       IF (status /= 0) STOP "fill_k_exf array size"
       CALL fill_k_exb(status, vk_exb(:,1:NSPEC))
       IF (status /= 0) STOP "fill_k_exb array size"
       CALL fill_k_exf_N2O5(status, vk_exf_N2O5)
       IF (status /= 0) STOP "fill_k_exf_N2O5 array size"
       CALL fill_k_exf_ClNO3(status, vk_exf_ClNO3)
       IF (status /= 0) STOP "fill_k_exf_ClNO3 array size"
       CALL fill_k_exf_BrNO3(status, vk_exf_BrNO3)
       IF (status /= 0) STOP "fill_k_exf_BrNO3 array size"

       DO jl=1,lproma
          conc(jl,1:NSPEC) = L_SPEC(1:NSPEC,JL)
       ENDDO

    ! KPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPP
       !        CALL kpp_integrate(TIME_STEP_LEN, conc)
       CALL kpp_integrate(TIME_STEP_LEN, conc, ierrf=kpp_stat, xNacc=kpp_steps, xNrej=kpp_rsteps)
    ! KPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPP
       conc(:,:) = MAX(conc(:,:), 0.0_DP) ! force positive results
       DO jl=1,lproma
          L_SPEC(1:NSPEC,JL) = conc(jl,1:NSPEC)
          if (kpp_stat(jl) /= 1) print*, "Error in SCAV_L_KPP: ", kpp_stat(jl),&
               "in Box: ", jl 
       ENDDO
       nsteps(:) = kpp_steps
       rsteps(:) = kpp_rsteps

    END IF

    CALL DEALLOC_CHEM_FIELDS

    RETURN

  END SUBROUTINE SCAV_AQCHEM_KPP
!==============================================================================
!================================= KPP ICE SCAVENGING =========================
  
  SUBROUTINE SCAV_ICECHEM_KPP (PTEMP, PPRESS, IWC_EFF, IWC,          &
                               ICERAD, TIME_STEP_LEN, PROC, LPROMA,  &
                               LIDX, JK, JROW, JS)

    ! NEW SUBROUTINE IN WHICH THE ICE PHASE CHEMISTRY IS BEING
    ! CALCULATED PRODUCED BY KPP 
    USE MESSY_MAIN_CONSTANTS_MEM, ONLY: R_GAS, N_A

    USE MESSY_SCAV_MEM,           ONLY: I_SPEC !, L_LG
    USE MESSY_SCAV_ICE
    USE MESSY_SCAV_LIQ,           ONLY: set_photo
! KPP VARIABLES
    USE MESSY_SCAV_I_KPP
    USE MESSY_SCAV_INP_KPP

! THIS CONTAINS THE PREPROCESSOR DIRECTIVE IF A VECTORIZED SOLVER IS USED
!#include "messy_scav_i_vec.inc"
!#DEFINE SCAV_VECTOR
    IMPLICIT NONE
    SAVE
    
    INTEGER,  INTENT(IN) :: PROC                ! 1 = CLOUD, 2 = RAIN
    INTEGER,  INTENT(IN) :: LPROMA, JK, JROW, LIDX(LPROMA), JS
    REAL(DP), INTENT(IN) :: TIME_STEP_LEN
    REAL(DP), INTENT(IN) :: PTEMP(LPROMA), IWC(LPROMA), IWC_EFF(LPROMA), &
                            ICERAD(LPROMA), PPRESS(LPROMA)
   ! INTEGER  :: I, JL, IDX2, status
    INTEGER  :: JL, status

    REAL(dp) :: cair(lproma)
    REAL(dp) :: conc(lproma,1:nspec)

    INTEGER  :: kpp_stat(lproma), kpp_steps(lproma), kpp_rsteps(lproma)

    !INTRINSIC :: REAL, SUM, TINY
    
    CALL ALLOC_CHEM_FIELDS(LPROMA)
! SETTING UP PHYSICAL PARAMETERS FOR EACH BOX
    CALL ICE_PHYSC_INIT(IWC_EFF, IWC , ICERAD, PTEMP, PPRESS, JS)

    CALL SET_PHOTO(LIDX, JK, JROW)
! CALCULATES THE EQUILIBRIUM COEFFICIENTS / TESTFACS
    CALL EQUIL_DISSOCIATION_COEFF
! CALCULATES THE TRANSFER COEFFICIENTS
    CALL TRANSFER_COEFFICIENT_ICE(PROC)

    CALL fill_jx (status, jx_vec)
    IF (status /= 0) STOP "fill_jx array size"

    CALL fill_temp (status,ptemp(1:LPROMA))
    IF (status /= 0) STOP "fill_temp array size"

    CALL fill_press (status, ppress)
    IF (status /= 0) STOP "fill_press array size"

    cair(:) = (N_A/1.E6_dp) * ppress(:) / (R_gas*ptemp(:))
    CALL fill_cair (status, cair)
    IF (status /= 0) STOP "fill_cair array size"

    CALL fill_lwc (status, viwc)
    IF (status /= 0) STOP "fill_lwc array size"

    CALL fill_cv_l (status, VCV_I)
    IF (status /= 0) STOP "fill_cv_l"
    
    CALL fill_k_exf(status, vk_exf(:,1:NSPEC))
    IF (status /= 0) STOP "fill_k_exf array size"
    CALL fill_k_exb(status, vk_exb(:,1:NSPEC))
    IF (status /= 0) STOP "fill_k_exb array size"
    CALL fill_k_exf_N2O5(status, vk_exf_N2O5)
    IF (status /= 0) STOP "fill_k_exf_N2O5 array size"
    CALL fill_k_exf_ClNO3(status, vk_exf_ClNO3)
    IF (status /= 0) STOP "fill_k_exf_ClNO3 array size"
    CALL fill_k_exf_BrNO3(status, vk_exf_BrNO3)
    IF (status /= 0) STOP "fill_k_exf_BrNO3 array size"

    DO jl=1,lproma
      conc(jl,1:NSPEC) = I_SPEC(1:NSPEC,JL)
    ENDDO 

    ! KPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPP
    CALL kpp_integrate(TIME_STEP_LEN, conc, ierrf=kpp_stat, xNacc=kpp_steps, xNrej=kpp_rsteps)
    ! KPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPP
    conc(:,:) = MAX(conc(:,:), 0.0_DP) ! force positive results
    DO jl=1,lproma
       I_SPEC(1:NSPEC,JL) = conc(jl,1:NSPEC)
       if (kpp_stat(jl) /= 1) print*, "Error in SCAV_I_KPP: ", kpp_stat(jl),&
            "in Box: ", jl        
    ENDDO
!    nsteps(:) = kpp_steps
!    rsteps(:) = kpp_rsteps

!!$#ifndef SCAV_VECTOR
!!$
!!$    DO JL=1,LPROMA
!!$!      PRINT*, "IN_KPP", L_LG
!!$      IF (L_LG) THEN
!!$        IDX2 = LIDX(JL)
!!$        IF (NINT(NCB(IDX2,JK,JROW)) < 1) THEN
!!$!          PRINT*, "SKIPPING LG_KPP", L_SPEC(IND_HNO3,JL), &
!!$!            IDX2, JK, JL, JROW, NCB(IDX2,JK,JROW)
!!$          CYCLE
!!$!        ELSE
!!$!          PRINT*, "KPP_LG NOT SKIPPED, BEFORE", L_SPEC(IND_HNO3,JL), &
!!$!            L_SPEC(IND_NO3M_L,JL),  L_SPEC(IND_HNO3_L,JL) 
!!$        ENDIF
!!$      ENDIF
!!$      C(1:NSPEC) = I_SPEC(1:NSPEC,JL)
!!$      K_EXF(1:NSPEC) = VK_EXF(JL,1:NSPEC)
!!$      K_EXB(1:NSPEC) = VK_EXB(JL,1:NSPEC)
!!$      CV_I           = VCV_I(JL)
!!$      TEMP           = VTEMP(JL)    
!!$
!!$      TIME = 0.
!!$      CALL UPDATE_RCONST
!!$      DO I=1,NSUBSTEPS
!!$        IF (LOGSTEPS) THEN
!!$          DT = TIME_STEP_LEN * SUBSTEP(I)
!!$        ELSE
!!$          DT=TIME_STEP_LEN/REAL(NSUBSTEPS,DP)
!!$          INFO(:) = 0
!!$          INFO(5) = NMAXKPP
!!$        ENDIF
!!$        ! FOR ROS2, THE DIFFERENCE BETWEEN THE TWO PARAMETERS _MUST_ BE DT:
!!$        CALL INTEGRATE(TIME, TIME+DT)
!!$      END DO
!!$    ENDDO
!!$#endif

  END SUBROUTINE SCAV_ICECHEM_KPP

!==============================================================================
! THIS SUBROUTINE CALCULATES THE LWC AND THE AVERAGE DROPLET DIAMETER 
! FOR ONE BOX FROM THE PRECIPITATION FLUX

  SUBROUTINE CALC_LWC_BC(FPREC,FSNOW, ZOLDCOV, &
                  ZRLWC, ZWASHFRAC, ZRDRAD, NLEV, JK, KCLTOP)

    IMPLICIT NONE
  
!   VARIABLES FROM THE INPUT
    REAL(DP), INTENT(IN)    :: FPREC, FSNOW, ZOLDCOV
    REAL(DP), INTENT(INOUT) :: ZRLWC, ZWASHFRAC, ZRDRAD
    INTEGER  :: NLEV, JK, KCLTOP
!   LOCAL VARIABLES
    REAL(DP) :: ZINFLUX,  ZRLWCX, ZRFLX, ZFWAT

    ZINFLUX=MAX(FPREC,0._DP)
    ZRFLX = 0._DP
    IF (ZINFLUX.LT.1.E-15_DP.AND. FSNOW.LT.1.E-15_DP) THEN                         ! ZINFLUX; [KG M-2 S-1]
       KCLTOP=NLEV
       ZRLWCX=0.
       ZRLWC=0.
    ENDIF
  
    ZWASHFRAC=ZOLDCOV
    IF (ZINFLUX.GT.1.E-15_DP .AND. &
         ZWASHFRAC.GT.0.01_DP) THEN

                ! -- CALCULATE RAIN PARAMETERS (KUMAR, 1985)
       ZRFLX=(ZINFLUX/ZOLDCOV)*3600. ! [MM HR-1], ([KG M-2 S-1]/[-])*3600
       ZRLWCX=72.*ZRFLX**0.88            ! [MG M-3] 

                ! ==============================================================
                ! MZ_HT_20030318+, THE FOLLOWING STATEMENTS CAN BE REMOVED 
                !     WHENEVER THE INPUT PARAMETERS PMLWC AND PMSIC HAVE BEEN 
                !     IMPLEMENTED PROPERLY. PMLWC REALLY REFLECTS THE LIQUID
                !     PHASE WHEREAS THE SNOW/ICE FRACTION IS REPRESENTED BY THE
                !     PARAMETER PMSIC
       ZFWAT = 1._DP     ! AS THE INPUT VARIABLE IS RAIN AND NOT SNOW IT CAN BE ASSUMED THAT THE RAIN IS LIQUID WATER
       ZRLWC=ZRLWCX*1E-6*ZFWAT  ! [L M-3]
      
    ENDIF
    
    IF (ZRFLX > 0._DP) THEN
       
       ZRDRAD=MAX(1E-20_DP,(0.3659*ZRFLX**0.21)) ! [MM], MASON, 1971
                                     ! SEE PAPER G-J. ROELOFS, PP 20,996, JGR, 1995
    ENDIF
                
  END SUBROUTINE CALC_LWC_BC

!==============================================================================
! this subroutine calculates the lwc and the average droplet diameter 
! for one box from the precipitation flux

  SUBROUTINE calc_lwc_bc_n(fprec, fsnow, cov, zrlwc, zrdrad, col_cov)

    IMPLICIT NONE
  
!   variables from the input
    REAL(dp), INTENT(in)    :: fprec, fsnow, cov
    REAL(dp), INTENT(inout) :: zrlwc, zrdrad, col_cov

!   local variables
    REAL(dp) :: zrlwcx, zrflx

    IF ( (fprec < 1.e-12_dp).AND.(fsnow < 1.e-12_dp) ) col_cov = 0._dp

    zrflx = 0._dp
    IF ( fprec < 1.e-12_dp) THEN
      zrlwcx = 0._dp
      zrlwc  = 0._dp
      zrdrad = 0._dp
    ELSE
      col_cov = MAX(col_cov,cov)
      IF (col_cov > 1.e-7_dp ) THEN
        ! -- calculate rain parameters (kumar, 1985)
        ! [mm hr-1], ([kg m-2 s-1]/[-])*3600
        zrflx=(fprec/col_cov)*3600._dp 
        ! [mg m-3] 
        zrlwcx=72._dp*zrflx**0.88_dp         
        ! [l m-3]
        zrlwc=zrlwcx*1.e-6_dp  
      ENDIF
    ENDIF
    
    IF (zrflx > 0._dp) THEN
       
       zrdrad=MAX(1e-20_dp,(0.3659_dp*zrflx**0.21_dp)) 
       ! radius is in [mm]
       ! calculation accoring to Mason, 1971
       ! see paper g-j. roelofs, pp 20,996, jgr, 1995
    ENDIF
                
  END SUBROUTINE calc_lwc_bc_n
  
!==============================================================================
! THIS SUBROUTINE CALCULATES THE LWC AND THE AVERAGE DROPLET DIAMETER
! FOR ONE BOX FROM THE PRECIPITATION FLUX

  SUBROUTINE CALC_LWC_BC_CV (ZRAINCV, ZCVCOVER, &
                  ZRLWC, ZWASHFRAC, ZRDRAD, JK, KCONBOT)

    IMPLICIT NONE
  
!   VARIABLES FROM THE INPUT
    REAL(DP), INTENT(IN) :: ZRAINCV, ZCVCOVER, KCONBOT
    REAL(DP), INTENT(INOUT) :: ZRLWC, ZWASHFRAC, ZRDRAD 
    INTEGER  :: JK
!   LOCAL VARIABLES
    REAL(DP) :: ZRLWCX, ZRFLX, ZFWAT

    ZRLWC=0._DP
    ZRFLX=0._DP
    ZRDRAD=0._DP
    ZWASHFRAC=ZCVCOVER

    IF (ZRAINCV > 1.0E-10_DP) THEN
     
       ZRFLX=ZRAINCV/ZWASHFRAC*3600.
       ZRLWCX=72.*ZRFLX**0.88

       ZFWAT = 1._DP     ! AS THE INPUT VARIABLE IS RAIN AND NOT SNOW IT CAN BE ASSUMED THAT THE RAIN IS LIQUID WATER
       ZRLWC=ZRLWCX*1E-6*ZFWAT

! ... REMOVED FOR VECTORIZATION (SX-6)
!       IF (ZWASHFRAC.LT. 1.E-30_DP) PRINT*, "WARNING, WASHFRAC", &
!         ZWASHFRAC, ZRAINCV, JK , ZCVCOVER, ZRLWC
    ENDIF
   
    
    IF (ZRFLX > 0._DP) THEN     
       ZRDRAD=MAX(1E-20_DP,(0.3659*ZRFLX**0.21)) ! [MM], MASON, 1971
                                     ! SEE PAPER G-J. ROELOFS, PP 20,996, JGR, 1995
    ENDIF
         

  END SUBROUTINE CALC_LWC_BC_CV

!===============================================================================
! THIS SUBROUTINE CALCULATES THE LWC AND THE AVERAGE DROPLET DIAMETER 
! FOR ONE BOX FROM THE PRECIPITATION FLUX

  SUBROUTINE CALC_LWC_IC_CV (RAINCVNEW, ZCVCOVER, &
                  ZRLWC, ZWASHFRAC, ZRDRAD, JK, KCONBOT)

    IMPLICIT NONE
  
!   VARIABLES FROM THE INPUT
    REAL(DP), INTENT(IN) :: RAINCVNEW, ZCVCOVER, KCONBOT
    REAL(DP), INTENT(INOUT) :: ZRLWC, ZWASHFRAC, ZRDRAD 
    INTEGER  :: JK
!   LOCAL VARIABLES
    REAL(DP) :: ZRLWCX, ZRFLX, ZFWAT

    ZRLWC=0._DP
    ZRFLX=0._DP
    ZRDRAD=0._DP
    ZWASHFRAC=ZCVCOVER
    IF (KCONBOT.GT.0..AND. RAINCVNEW.GT.1.0E-15_DP) THEN
     

       ZRFLX=RAINCVNEW/ZWASHFRAC*3600.
       ZRLWCX=72.*ZRFLX**0.88
       
       ZFWAT = 1._DP     
! AS THE INPUT VARIABLE IS RAIN AND NOT SNOW IT CAN BE ASSUMED 
! THAT THE RAIN IS LIQUID WATER
       ZRLWC=ZRLWCX*1E-6*ZFWAT

! ... REMOVED FOR VECTORIZATION (SX-6)
!       IF (ZWASHFRAC.LT. 1.E-30_DP) PRINT*, "WARNING, WASHFRAC IC", &
!         ZWASHFRAC, RAINCVNEW, JK , ZCVCOVER, ZRLWC

    ENDIF
   
    
    IF (ZRLWC.GT.1.E-7_DP) THEN     
       ZRDRAD=0.02_DP                       ! 20?M = 0.02MM
    ENDIF
         

  END SUBROUTINE CALC_LWC_IC_CV

  !=============================================================================


END MODULE MESSY_SCAV
