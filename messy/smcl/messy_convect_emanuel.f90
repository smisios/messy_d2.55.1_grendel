MODULE MESSY_CONVECT_EMANUEL

  USE messy_main_constants_mem,    ONLY: dp, rowl=>rho_H2O, cpd=>cp_air, &
                                         cpv, rd, rv, g, Tmelt

  IMPLICIT NONE
  PRIVATE

  INTEGER, PUBLIC :: ntrans
  PUBLIC :: convect_emanuel
  PUBLIC :: dp
!   ***                    SPECIFY PARAMETERS                        ***
!
!   *** ELCRIT IS THE AUTOCONVERSION THERSHOLD WATER CONTENT (gm/gm) ***
!   ***  TLCRIT IS CRITICAL TEMPERATURE BELOW WHICH THE AUTO-        ***
!   ***       CONVERSION THRESHOLD IS ASSUMED TO BE ZERO             ***
!   ***     (THE AUTOCONVERSION THRESHOLD VARIES LINEARLY            ***
!   ***               BETWEEN 0 C AND TLCRIT)                        ***
!   ***   ENTP IS THE COEFFICIENT OF MIXING IN THE ENTRAINMENT       ***
!   ***                       FORMULATION                            ***
!   ***  SIGD IS THE FRACTIONAL AREA COVERED BY UNSATURATED DNDRAFT  ***
!   ***  SIGS IS THE FRACTION OF PRECIPITATION FALLING OUTSIDE       ***
!   ***                        OF CLOUD                              ***
!   ***        OMTRAIN IS THE ASSUMED FALL SPEED (P/s) OF RAIN       ***
!   ***     OMTSNOW IS THE ASSUMED FALL SPEED (P/s) OF SNOW          ***
!   ***  COEFFR IS A COEFFICIENT GOVERNING THE RATE OF EVAPORATION   ***
!   ***                          OF RAIN                             ***
!   ***  COEFFS IS A COEFFICIENT GOVERNING THE RATE OF EVAPORATION   ***
!   ***                          OF SNOW                             ***
!   ***     CU IS THE COEFFICIENT GOVERNING CONVECTIVE MOMENTUM      ***
!   ***                         TRANSPORT                            ***
!   ***    DTMAX IS THE MAXIMUM NEGATIVE TEMPERATURE PERTURBATION    ***
!   ***        A LIFTED PARCEL IS ALLOWED TO HAVE BELOW ITS LFC      ***
!   ***    ALPHA AND DAMP ARE PARAMETERS THAT CONTROL THE RATE OF    ***
!   ***                 APPROACH TO QUASI-EQUILIBRIUM                ***
!   ***   (THEIR STANDARD VALUES ARE  0.20 AND 0.1, RESPECTIVELY)    ***
!   ***                   (DAMP MUST BE LESS THAN 1)                 ***
!
    ! original value
!  REAL(dp), PARAMETER ::      ELCRIT  =  0.0011_dp
!  REAL(dp), PARAMETER ::      ELCRIT  =  0.003_dp ! T63L87
!  REAL(dp), PARAMETER ::      ELCRIT  =  0.006_dp ! T63L87
  REAL(dp), PARAMETER ::      ELCRIT  =  0.0015_dp ! T63L87
  REAL(dp), PARAMETER ::      TLCRIT  = -55.0_dp
  REAL(dp), PARAMETER ::      ENTP    =  1.5_dp
  REAL(dp), PARAMETER ::      SIGD    =  0.05_dp
  REAL(dp), PARAMETER ::      SIGS    =  0.12_dp
  REAL(dp), PARAMETER ::      OMTRAIN =  50.0_dp
  REAL(dp), PARAMETER ::      OMTSNOW =  5.5_dp
  REAL(dp), PARAMETER ::      COEFFR  =  1.0_dp
  REAL(dp), PARAMETER ::      COEFFS  =  0.8_dp
  REAL(dp), PARAMETER ::      CU      =  0.7_dp
  REAL(dp), PARAMETER ::      BETA    =  10.0_dp
  REAL(dp), PARAMETER ::      DTMAX   =  0.9_dp
!!$  REAL(dp), PARAMETER ::      ALPHA   =  0.2_dp
!!$  REAL(dp), PARAMETER ::      DAMP    =  0.1_dp
  REAL(dp), PARAMETER ::      ALPHA   =  0.05_dp
  REAL(dp), PARAMETER ::      DAMP    =  0.2_dp

!   ***        ASSIGN VALUES OF THERMODYNAMIC CONSTANTS,        ***
!   ***            GRAVITY, AND LIQUID WATER DENSITY.           ***
!   ***             THESE SHOULD BE CONSISTENT WITH             ***
!   ***              THOSE USED IN CALLING PROGRAM              ***
!   ***     NOTE: THESE ARE ALSO SPECIFIED IN SUBROUTINE TLIFT  ***
!
! from messy_main_constants_mem:
!      G=9.8  
!      ROWL=1000.0
!      CPD=1005.7
!      CPV=1870.0
!      RV=461.5
!      RD=287.04
  REAL(dp), PARAMETER ::    CL     = 2500.0_dp
  REAL(dp), PARAMETER ::    LV0    = 2.501E6_dp
!
  REAL(dp), PARAMETER ::    CPVMCL = CL - CPV 
  REAL(dp), PARAMETER ::    EPS    = RD / RV
  REAL(dp), PARAMETER ::    EPSI   = 1._dp / EPS
  REAL(dp), PARAMETER ::    GINV   = 1._dp / G
CONTAINS
!----------------------------------------------------------------------------

!***************************************************************************
!*****                       SUBROUTINE CONVECT                        *****
!*****                          VERSION 4.3c                           *****
!*****                          20 May, 2002                           *****
!*****                          Kerry Emanuel                          *****
!***************************************************************************
  SUBROUTINE CONVECT_EMANUEL(NA,                             &
          T,   Q,    QS,     U,    V,      TRA,    P,    PH, &
          ND,  NL,   NTRA,   DELT, IFLAG,  FT,     FQ,   FU, &
          FV,  FTRA, PRECIP, WD,   TPRIME, QPRIME, CBMF,     &
          MINORIG,   ICB,    INB,  PUMF,   PDMF,   CV_PREC,  &
          CV_FORM,   PUER,   PUDR, PDER,   PDDR,   CAPE,     &
          LWC, IWC,  RFORM,  SFORM )

!-----------------------------------------------------------------------------
!    *** On input:      ***
!
!     T:   Array of absolute temperature (K) of dimension ND, with first
!           index corresponding to lowest model level. Note that this array
!           will be altered by the subroutine if dry convective adjustment
!           occurs and if IPBL is not equal to 0.
!
!     Q:   Array of specific humidity (gm/gm) of dimension ND, with first
!            index corresponding to lowest model level. Must be defined
!            at same grid levels as T. Note that this array will be altered
!            if dry convective adjustment occurs and if IPBL is not equal to 0.
!
!     QS:  Array of saturation specific humidity of dimension ND, with first
!            index corresponding to lowest model level. Must be defined
!            at same grid levels as T. Note that this array will be altered
!            if dry convective adjustment occurs and if IPBL is not equal to 0.
!
!     U:   Array of zonal wind velocity (m/s) of dimension ND, witth first
!            index corresponding with the lowest model level. Defined at
!            same levels as T. Note that this array will be altered if
!            dry convective adjustment occurs and if IPBL is not equal to 0.
!
!     V:   Same as U but for meridional velocity.
!
!     TRA: Array of passive tracer mixing ratio, of dimensions (ND,NTRA),
!            where NTRA is the number of different tracers. If no
!            convective tracer transport is needed, define a dummy
!            input array of dimension (ND,1). Tracers are defined at
!            same vertical levels as T. Note that this array will be altered
!            if dry convective adjustment occurs and if IPBL is not equal to 0.
!
!     P:   Array of pressure (mb) of dimension ND, with first
!            index corresponding to lowest model level. Must be defined
!            at same grid levels as T.
!
!     PH:  Array of pressure (mb) of dimension ND+1, with first index
!            corresponding to lowest level. These pressures are defined at
!            levels intermediate between those of P, T, Q and QS. The first
!            value of PH should be greater than (i.e. at a lower level than)
!            the first value of the array P.
!
!     ND:  The dimension of the arrays T,Q,QS,P,PH,FT and FQ
!
!     NL:  The maximum number of levels to which convection can
!            penetrate, plus 1.
!            NL MUST be less than or equal to ND-1.
!
!     NTRA:The number of different tracers. If no tracer transport
!            is needed, set this equal to 1. (On most compilers, setting
!            NTRA to 0 will bypass tracer calculation, saving some CPU.)  
!
!     DELT: The model time step (sec) between calls to CONVECT
!
!----------------------------------------------------------------------------
!    ***   On Output:         ***
!
!     IFLAG: An output integer whose value denotes the following:
!
!                VALUE                        INTERPRETATION
!                -----                        --------------
!                  0               No moist convection; atmosphere is not
!                                  unstable, or surface temperature is less
!                                  than 250 K or surface specific humidity
!                                  is non-positive.
!                  1               Moist convection occurs.
!                  2               No moist convection: lifted condensation
!                                  level is above the 200 mb level.
!                  3               No moist convection: cloud base is higher
!                                  then the level NL-1.
!                  4               Moist convection occurs, but a CFL condition
!                                  on the subsidence warming is violated. This
!                                  does not cause the scheme to terminate.
!
!     FT:   Array of temperature tendency (K/s) of dimension ND, defined at same
!             grid levels as T, Q, QS and P.
!
!     FQ:   Array of specific humidity tendencies ((gm/gm)/s) of dimension ND,
!             defined at same grid levels as T, Q, QS and P.
!
!     FU:   Array of forcing of zonal velocity (m/s^2) of dimension ND,
!             defined at same grid levels as T.
!
!     FV:   Same as FU, but for forcing of meridional velocity.
!
!     FTRA: Array of forcing of tracer content, in tracer mixing ratio per
!             second, defined at same levels as T. Dimensioned (ND,NTRA).
!
!     PRECIP: Scalar convective precipitation rate (mm/day).
!
!     WD:    A convective downdraft velocity scale. For use in surface
!             flux parameterizations. See convect.ps file for details.
!
!     TPRIME: A convective downdraft temperature perturbation scale (K).
!              For use in surface flux parameterizations. See convect.ps
!              file for details.
!
!     QPRIME: A convective downdraft specific humidity
!              perturbation scale (gm/gm).
!              For use in surface flux parameterizations. See convect.ps
!              file for details.
!
!     CBMF:   The cloud base mass flux ((kg/m**2)/s). THIS SCALAR VALUE MUST
!              BE STORED BY THE CALLING PROGRAM AND RETURNED TO CONVECT AT
!              ITS NEXT CALL. That is, the value of CBMF must be "remembered"
!              by the calling program between calls to CONVECT.
!
!------------------------------------------------------------------------------
!
!    ***  THE PARAMETER NA SHOULD IN GENERAL BE GREATER THAN   ***
!    ***                OR EQUAL TO  ND + 1                    ***
!
    ! mz_ht_20070214+
      !PARAMETER (NA=70)           
      INTEGER, INTENT(IN) :: NA
      INTEGER, INTENT(IN) :: ND
      INTEGER, INTENT(IN) :: NL
      INTEGER, INTENT(IN) :: NTRA
    ! mz_ht_20070214-
!
      INTEGER  :: NENT(NA)
      REAL(dp) :: T(ND),Q(ND),QS(ND),U(ND),V(ND),TRA(ND,NTRA),P(ND),PH(ND)
      REAL(dp) :: FT(ND),FQ(ND),FU(ND),FV(ND) !FTRA(ND,NTRA)
      ! op_mm_20140327+
      ! compiler workaround (G95 0.92,0.93) 
       REAL(dp) ::FTRA(ND,0:NTRA)
      ! op_mm_20140327-

      REAL(dp) :: UENT(NA,NA),VENT(NA,NA),TRAENT(NA,NA,NTRA),TRATM(NA)
      REAL(dp) :: UP(NA),VP(NA),TRAP(NA,NTRA)
      REAL(dp) :: M(NA),MP(NA),MENT(NA,NA),QENT(NA,NA),ELIJ(NA,NA)
      REAL(dp) :: SIJ(NA,NA),TVP(NA),TV(NA),WATER(NA)
      REAL(dp) :: QP(NA),EP(NA),TH(NA),WT(NA),EVAP(NA),CLW(NA)
      REAL(dp) :: SIGP(NA),TP(NA),TOLD(NA),CPN(NA)
      REAL(dp) :: LV(NA),LVCP(NA),H(NA),HP(NA),GZ(NA),HM(NA)


      ! mz_ht_20070214+
      ! additional variable declaration required because of IMPLICIT NONE
      INTEGER  :: I, J, K
      INTEGER  :: INB1, IHMIN,  IPBL, IFLAG 
      INTEGER  :: JC, JN, JTT, NK

      REAL(dp) :: DAMPS, DELT0, DTMA, DTMIN, DTPBL, DPHINV, DELTI
      REAL(dp) :: CBMFOLD
      REAL(dp) :: DBO, DBOSUM
      REAL(dp) :: TVPPLCL, TVAPLCL
      REAL(dp) :: frac, DEFRAC, BY, BYP
      REAL(dp) :: CAPE, CAPEM, AHMIN, ALV, ALVNEW
      REAL(dp) :: ELACRIT, EPMAX, TCA, CHI, PLCL
      REAL(dp) :: FUOLD, FVOLD, UAV, VAV, FQOLD, FTOLD
      REAL(dp) :: QNEW, TNEW, TVX, TVY, TG, TC, X
      REAL(dp) :: AHM, AHMAX, RM, UM, VM, THBAR
      REAL(dp) :: RDCP, CPINV
      REAL(dp) :: SMIN, SMID, SJMAX, SJMIN, SUM
      REAL(dp) :: A2, B6, C6, BF2, AM, AMDE
      REAL(dp) :: RH, QTI, ANUM, DENOM, DEI, ALTEM, ALT
      REAL(dp) :: CWAT, STEMP, QP1, SCRIT, QSUM, COEFF
      REAL(dp) :: DELP, DELM, BSUM, WDTRAIN, AWAT, AFAC
      REAL(dp) :: REVAP, DHDP, FAC, RAT, QSTM, AMP1, AD
      REAL(dp) :: ENTS, DPINV, SIGT, QSM, ASIJ
      REAL(dp) :: FTRAOLD, TRAAV

      INTEGER  :: ICB, INB, MINORIG
      REAL(dp) :: CBMF, QPRIME, TPRIME, WD, PRECIP, DELT
      REAL(dp) :: PUMF(ND), PDMF(ND)           ! storing of mass fluxes
      REAL(dp) :: PUER(ND), PUDR(ND)           ! storing updraft entrainment
                                               ! and detrainment of mass fluxes
      REAL(dp) :: PDER(ND), PDDR(ND)           ! storing downdraft entrainment 
                                               ! and detrainment of mass fluxes
      REAL(dp) :: CV_PREC(ND)                  ! storing of cv precipitation
      REAL(dp) :: CV_FORM(ND)                  ! storing of cv precip. formation
      REAL(DP) :: LWC(ND), IWC(ND), RFORM(ND), SFORM(ND)

      ! mz_ht_20070214-

!
! -----------------------------------------------------------------------
!
!   ***                     Specify Switches                         ***
!
!   ***   IPBL: Set to zero to bypass dry adiabatic adjustment       ***
!   ***    Any other value results in dry adiabatic adjustment       ***
!   ***     (Zero value recommended for use in models with           ***
!   ***                   boundary layer schemes)                    ***
!
!   ***   MINORIG: Lowest level from which convection may originate  ***
!   ***     (Should be first model level at which T is defined       ***
!   ***      for models using bulk PBL schemes; otherwise, it should ***
!   ***      be the first model level at which T is defined above    ***
!   ***                      the surface layer)                      ***
!
        IPBL=0
!        MINORIG=1   is now set in the driving routine ! mz_ht_20070215
!
!------------------------------------------------------------------------------
!
!
        ICB = 0
        INB = 0
        DELTI=1.0_dp/DELT
!
!           ***  INITIALIZE OUTPUT ARRAYS AND PARAMETERS  ***
!
        DO 5 I=1,ND
         FT(I)=0.0_dp
         FQ(I)=0.0_dp
         FU(I)=0.0_dp
         FV(I)=0.0_dp
         DO 4 J=1,NTRA
          FTRA(I,J)=0.0_dp
    4    CONTINUE
    5   CONTINUE
        DO 7 I=1,NL+1
         RDCP=(RD*(1._dp-Q(I))+Q(I)*RV)/(CPD*(1._dp-Q(I))+Q(I)*CPV)
         TH(I)=T(I)*(1000.0_dp/P(I))**RDCP
    7   CONTINUE
        PRECIP=0.0_dp
        WD=0.0_dp
        TPRIME=0.0_dp
        QPRIME=0.0_dp
        IFLAG=0
!
        IF(IPBL.NE.0)THEN
!
!     ***            PERFORM DRY ADIABATIC ADJUSTMENT            ***
!
        JC=0
        DO 30 I=NL-1,1,-1
         JN=0
          SUM=TH(I)*(1._dp+Q(I)*EPSI-Q(I))
         DO 10 J=I+1,NL
          SUM=SUM+TH(J)*(1._dp+Q(J)*EPSI-Q(J))
          THBAR=SUM/FLOAT(J+1-I)
          IF((TH(J)*(1._dp+Q(J)*EPSI-Q(J))).LT.THBAR) JN=J
   10    CONTINUE
         IF(I.EQ.1)JN=MAX(JN,2)
         IF(JN.EQ.0)GOTO 30
   12    CONTINUE
         AHM=0.0_dp
         RM=0.0_dp
         UM=0.0_dp
         VM=0.0_dp
         DO K=1,NTRA
          TRATM(K)=0.0_dp
         END DO
         DO 15 J=I,JN
          AHM=AHM+(CPD*(1._dp-Q(J))+Q(J)*CPV)*T(J)*(PH(J)-PH(J+1))
          RM=RM+Q(J)*(PH(J)-PH(J+1))
          UM=UM+U(J)*(PH(J)-PH(J+1))
          VM=VM+V(J)*(PH(J)-PH(J+1))
          DO K=1,NTRA
           TRATM(K)=TRATM(K)+TRA(J,K)*(PH(J)-PH(J+1))
          END DO
   15    CONTINUE
         DPHINV=1._dp/(PH(I)-PH(JN+1))
         RM=RM*DPHINV
         UM=UM*DPHINV
         VM=VM*DPHINV
         DO K=1,NTRA
          TRATM(K)=TRATM(K)*DPHINV
         END DO
         A2=0.0_dp
         DO 20 J=I,JN
          Q(J)=RM
          U(J)=UM
          V(J)=VM
          DO K=1,NTRA
           TRA(J,K)=TRATM(K)
          END DO
          RDCP=(RD*(1._dp-Q(J))+Q(J)*RV)/(CPD*(1._dp-Q(J))+Q(J)*CPV)  
          X=(0.001_dp*P(J))**RDCP
          TOLD(J)=T(J)
          T(J)=X
          A2=A2+(CPD*(1._dp-Q(J))+Q(J)*CPV)*X*(PH(J)-PH(J+1))
   20    CONTINUE
         DO 25 J=I,JN
          TH(J)=AHM/A2
          T(J)=T(J)*TH(J)
          TC=TOLD(J)-TMELT
          ALV=LV0-CPVMCL*TC
          QS(J)=QS(J)+QS(J)*(1._dp+QS(J)*(EPSI-1._dp))*ALV*(T(J)-   &
                TOLD(J))/(RV*TOLD(J)*TOLD(J))
   25    CONTINUE
         IF((TH(JN+1)*(1._dp+Q(JN+1)*EPSI-Q(JN+1))).LT.          &
            (TH(JN)*(1._dp+Q(JN)*EPSI-Q(JN))))THEN
          JN=JN+1
          GOTO 12
         END IF
         IF(I.EQ.1)JC=JN 
   30   CONTINUE
!
!   ***   Remove any supersaturation that results from adjustment ***
!
      IF(JC.GT.1)THEN
       DO 38 J=1,JC
          IF(QS(J).LT.Q(J))THEN 
           ALV=LV0-CPVMCL*(T(J)-TMELT)  
           TNEW=T(J)+ALV*(Q(J)-QS(J))/(CPD*(1._dp-Q(J))+                  &
                CL*Q(J)+QS(J)*(CPV-CL+ALV*ALV/(RV*T(J)*T(J))))
           ALVNEW=LV0-CPVMCL*(TNEW-TMELT)
           QNEW=(ALV*Q(J)-(TNEW-T(J))*(CPD*(1._dp-Q(J))+CL*Q(J)))/ALVNEW
           PRECIP=PRECIP+24._dp*3600._dp*1.0E5_dp*(PH(J)-PH(J+1))*        &
                  (Q(J)-QNEW)/(G*DELT*ROWL)
           T(J)=TNEW
           Q(J)=QNEW
           QS(J)=QNEW
          END IF     
   38  CONTINUE  
      END IF
!
      END IF
!
!  *** CALCULATE ARRAYS OF GEOPOTENTIAL, HEAT CAPACITY AND STATIC ENERGY
!  
        GZ(1)=0.0_dp
        CPN(1)=CPD*(1._dp-Q(1))+Q(1)*CPV
        H(1)=T(1)*CPN(1)
        LV(1)=LV0-CPVMCL*(T(1)-TMELT)
        HM(1)=LV(1)*Q(1)
        TV(1)=T(1)*(1._dp+Q(1)*EPSI-Q(1))
        AHMIN=1.0E12_dp
        IHMIN=NL
        DO 40 I=2,NL+1
          TVX=T(I)*(1._dp+Q(I)*EPSI-Q(I))
          TVY=T(I-1)*(1._dp+Q(I-1)*EPSI-Q(I-1))
          GZ(I)=GZ(I-1)+0.5_dp*RD*(TVX+TVY)*(P(I-1)-P(I))/PH(I)
          CPN(I)=CPD*(1._dp-Q(I))+CPV*Q(I)
          H(I)=T(I)*CPN(I)+GZ(I)
          LV(I)=LV0-CPVMCL*(T(I)-TMELT)
          HM(I)=(CPD*(1._dp-Q(I))+CL*Q(I))*(T(I)-T(1))+LV(I)*Q(I)+GZ(I)
          TV(I)=T(I)*(1._dp+Q(I)*EPSI-Q(I))
!
!  ***  Find level of minimum moist static energy    ***
!
          IF(I.GE.MINORIG.AND.HM(I).LT.AHMIN.AND.HM(I).LT.HM(I-1))THEN
           AHMIN=HM(I)
           IHMIN=I
          END IF
   40   CONTINUE
        IHMIN=MIN(IHMIN, NL-1)
!
!  ***     Find that model level below the level of minimum moist       ***
!  ***  static energy that has the maximum value of moist static energy ***
!
        AHMAX=0.0_dp
        DO 42 I=MINORIG,IHMIN
         IF(HM(I).GT.AHMAX)THEN
          NK=I
          AHMAX=HM(I)
         END IF
   42   CONTINUE
!
!  ***  CHECK WHETHER PARCEL LEVEL TEMPERATURE AND SPECIFIC HUMIDITY   ***
!  ***                          ARE REASONABLE                         ***
!  ***      Skip convection if HM increases monotonically upward       ***
!
        IF(T(NK).LT.250.0_dp.OR.Q(NK).LE.0.0_dp.OR.IHMIN.EQ.(NL-1))THEN
         IFLAG=0
         CBMF=0.0_dp
!         print*, "return1"
         RETURN
        END IF

!
!   ***  CALCULATE LIFTED CONDENSATION LEVEL OF AIR AT PARCEL ORIGIN LEVEL ***
!   ***       (WITHIN 0.2% OF FORMULA OF BOLTON, MON. WEA. REV.,1980)      ***
!
        RH=Q(NK)/QS(NK)
        CHI=T(NK)/(1669.0_dp-122.0_dp*RH-T(NK))
        PLCL=P(NK)*(RH**CHI)
        IF(PLCL.LT.200.0_dp.OR.PLCL.GE.2000.0_dp)THEN
         IFLAG=2
         CBMF=0.0_dp
!         print*, "return2"
         RETURN
        END IF
!
!   ***  CALCULATE FIRST LEVEL ABOVE LCL (=ICB)  ***
!
        ICB=NL-1
        DO 50 I=NK+1,NL
         IF(P(I).LT.PLCL)THEN
          ICB=MIN(ICB,I)
         END IF
   50   CONTINUE
        IF(ICB.GE.(NL-1))THEN
         IFLAG=3
         CBMF=0.0_dp
!         print*, "return3"
         RETURN
        END IF
!
!   *** FIND TEMPERATURE UP THROUGH ICB AND TEST FOR INSTABILITY           ***
!
!   *** SUBROUTINE TLIFT CALCULATES PART OF THE LIFTED PARCEL VIRTUAL      ***
!   ***  TEMPERATURE, THE ACTUAL TEMPERATURE AND THE ADIABATIC             ***
!   ***                   LIQUID WATER CONTENT                             ***
!
        CALL TLIFT(P,T,Q,QS,GZ,ICB,NK,TVP,TP,CLW,ND,NL,1)
        DO 54 I=NK,ICB
         TVP(I)=TVP(I)-TP(I)*Q(NK)
   54   CONTINUE
!
!   ***  If there was no convection at last time step and parcel    ***
!   ***       is stable at ICB then skip rest of calculation        ***
!
        IF(CBMF.EQ.0.0_dp.AND.TVP(ICB).LE.(TV(ICB)-DTMAX))THEN
         IFLAG=0
!         print*, "return4"
         RETURN
        END IF

!
!   ***  IF THIS POINT IS REACHED, MOIST CONVECTIVE ADJUSTMENT IS NECESSARY ***
!
        IF(IFLAG.NE.4)IFLAG=1
!
!   ***  FIND THE REST OF THE LIFTED PARCEL TEMPERATURES          ***
!
        CALL TLIFT(P,T,Q,QS,GZ,ICB,NK,TVP,TP,CLW,ND,NL,2)
!
!   ***  SET THE PRECIPITATION EFFICIENCIES AND THE FRACTION OF   ***
!   ***          PRECIPITATION FALLING OUTSIDE OF CLOUD           ***
!   ***      THESE MAY BE FUNCTIONS OF TP(I), P(I) AND CLW(I)     ***
!                 
        DO 57 I=1,NK
         EP(I)=0.0_dp
         SIGP(I)=SIGS
   57   CONTINUE
        DO 60 I=NK+1,NL
         TCA=TP(I)-TMELT
         IF(TCA.GE.0.0_dp)THEN
          ELACRIT=ELCRIT
         ELSE
          ELACRIT=ELCRIT*(1.0_dp-TCA/TLCRIT)
         END IF
         ELACRIT=MAX(ELACRIT,0.0_dp)
         EPMAX=0.999_dp
         EP(I)=EPMAX*(1.0_dp-ELACRIT/MAX(CLW(I),1.0E-8_dp))
         EP(I)=MAX(EP(I),0.0_dp)
         EP(I)=MIN(EP(I),EPMAX)
         SIGP(I)=SIGS
   60   CONTINUE
!
!   ***       CALCULATE VIRTUAL TEMPERATURE AND LIFTED PARCEL     ***
!   ***                    VIRTUAL TEMPERATURE                    ***
!
        DO 64 I=ICB+1,NL
         TVP(I)=TVP(I)-TP(I)*Q(NK)
   64   CONTINUE
        TVP(NL+1)=TVP(NL)-(GZ(NL+1)-GZ(NL))/CPD
!
!   ***        NOW INITIALIZE VARIOUS ARRAYS USED IN THE COMPUTATIONS       ***
!
        DO 70 I=1,NL+1
         HP(I)=H(I)
         NENT(I)=0
         WATER(I)=0.0_dp
         EVAP(I)=0.0_dp
         WT(I)=OMTSNOW
         MP(I)=0.0_dp
         M(I)=0.0_dp
         LVCP(I)=LV(I)/CPN(I)
         DO 70 J=1,NL+1
          QENT(I,J)=Q(J)
          ELIJ(I,J)=0.0_dp
          MENT(I,J)=0.0_dp
          SIJ(I,J)=0.0_dp
          UENT(I,J)=U(J)
          VENT(I,J)=V(J)
          DO 70 K=1,NTRA
           TRAENT(I,J,K)=TRA(J,K)
   70   CONTINUE
        QP(1)=Q(1)
        UP(1)=U(1)
        VP(1)=V(1)
        DO 71 I=1,NTRA
         TRAP(1,I)=TRA(1,I)
   71      CONTINUE
        DO 72 I=2,NL+1
         QP(I)=Q(I-1)
         UP(I)=U(I-1)
         VP(I)=V(I-1)
         DO 72 J=1,NTRA
          TRAP(I,J)=TRA(I-1,J)
   72         CONTINUE
!
!  ***  FIND THE FIRST MODEL LEVEL (INB1) ABOVE THE PARCEL'S      ***
!  ***          HIGHEST LEVEL OF NEUTRAL BUOYANCY                 ***
!  ***     AND THE HIGHEST LEVEL OF POSITIVE CAPE (INB)           ***
!
        CAPE=0.0_dp
        CAPEM=0.0_dp
        INB=ICB+1
        INB1=INB
        BYP=0.0_dp
        DO 82 I=ICB+1,NL-1
         BY=(TVP(I)-TV(I))*(PH(I)-PH(I+1))/P(I)
         CAPE=CAPE+BY
         IF(BY.GE.0.0_dp)INB1=I+1
         IF(CAPE.GT.0.0_dp)THEN
          INB=I+1
          BYP=(TVP(I+1)-TV(I+1))*(PH(I+1)-PH(I+2))/P(I+1)
          CAPEM=CAPE
         END IF
   82       CONTINUE
        INB=MAX(INB,INB1)
        CAPE=CAPEM+BYP
        DEFRAC=CAPEM-CAPE
        DEFRAC=MAX(DEFRAC,0.001_dp)
        FRAC=-CAPE/DEFRAC
        FRAC=MIN(FRAC,1.0_dp)
        FRAC=MAX(FRAC,0.0_dp)

!
!   ***   CALCULATE LIQUID WATER STATIC ENERGY OF LIFTED PARCEL   ***
!
        DO 95 I=ICB,INB
         HP(I)=H(NK)+(LV(I)+(CPD-CPV)*T(I))*EP(I)*CLW(I)
   95   CONTINUE                  
!
!   ***  CALCULATE CLOUD BASE MASS FLUX AND RATES OF MIXING, M(I),  ***
!   ***                   AT EACH MODEL LEVEL                       ***
!
        DBOSUM=0.0_dp
!   
!   ***     INTERPOLATE DIFFERENCE BETWEEN LIFTED PARCEL AND      ***
!   ***  ENVIRONMENTAL TEMPERATURES TO LIFTED CONDENSATION LEVEL  ***
!      
        TVPPLCL=TVP(ICB-1)-RD*TVP(ICB-1)*(P(ICB-1)-PLCL)/(CPN(ICB-1)*P(ICB-1))
        TVAPLCL=TV(ICB)+(TVP(ICB)-TVP(ICB+1))*(PLCL-P(ICB))/(P(ICB)-P(ICB+1))
        DTPBL=0.0_dp
        DO 96 I=NK,ICB-1
         DTPBL=DTPBL+(TVP(I)-TV(I))*(PH(I)-PH(I+1))
   96   CONTINUE
        DTPBL=DTPBL/(PH(NK)-PH(ICB))
        DTMIN=TVPPLCL-TVAPLCL+DTMAX+DTPBL
        DTMA=DTMIN
!
!   ***  ADJUST CLOUD BASE MASS FLUX   ***
!
      CBMFOLD=CBMF
      DELT0=300.0_dp
      DAMPS=DAMP*DELT/DELT0 
      CBMF=(1._dp-DAMPS)*CBMF+0.1_dp*ALPHA*DTMA 
      CBMF=MAX(CBMF,0.0_dp)
!
!   *** If cloud base mass flux is zero, skip rest of calculation  ***
!
      IF(CBMF.EQ.0.0_dp.AND.CBMFOLD.EQ.0.0_dp)THEN
!        print*, "return5"
       RETURN
      END IF

!
!   ***   CALCULATE RATES OF MIXING,  M(I)   ***
!
      M(ICB)=0.0_dp
      DO 103 I=ICB+1,INB
       K=MIN(I,INB1)
       DBO=ABS(TV(K)-TVP(K)) + ENTP*0.02_dp*(PH(K)-PH(K+1))
       DBOSUM=DBOSUM+DBO
       M(I)=CBMF*DBO
  103 CONTINUE
      DO 110 I=ICB+1,INB
       M(I)=M(I)/DBOSUM  
  110 CONTINUE     
!
!   ***  CALCULATE ENTRAINED AIR MASS FLUX (MENT), TOTAL WATER MIXING  ***
!   ***     RATIO (QENT), TOTAL CONDENSED WATER (ELIJ), AND MIXING     ***
!   ***                        FRACTION (SIJ)                          ***
!
        DO 170 I=ICB+1,INB
         QTI=Q(NK)-EP(I)*CLW(I)
         ! mz_ht_20070919+
         IF (T(J) > TMELT) THEN
           LWC(J)   = CLW(J)
           RFORM(j) = EP(j) * CLW(j)
         ELSE
           IWC(J)   = CLW(J)
           SFORM(j) = EP(j) * CLW(j)
         ENDIF
         ! mz_ht_20070919-
         DO 160 J=ICB,INB
          BF2=1._dp+LV(J)*LV(J)*QS(J)/(RV*T(J)*T(J)*CPD)
          ANUM=H(J)-HP(I)+(CPV-CPD)*T(J)*(QTI-Q(J))
          DENOM=H(I)-HP(I)+(CPD-CPV)*(Q(I)-QTI)*T(J)
          DEI=DENOM
          IF(ABS(DEI).LT.0.01_dp)DEI=0.01_dp
          SIJ(I,J)=ANUM/DEI
          SIJ(I,I)=1.0_dp
          ALTEM=SIJ(I,J)*Q(I)+(1._dp-SIJ(I,J))*QTI-QS(J)
          ALTEM=ALTEM/BF2
          CWAT=CLW(J)*(1._dp-EP(J))
          STEMP=SIJ(I,J)
          IF((STEMP.LT.0.0_dp.OR.STEMP.GT.1.0_dp.OR.        &
              ALTEM.GT.CWAT).AND.J.GT.I)THEN
           ANUM=ANUM-LV(J)*(QTI-QS(J)-CWAT*BF2)
           DENOM=DENOM+LV(J)*(Q(I)-QTI)
           IF(ABS(DENOM).LT.0.01_dp)DENOM=0.01_dp
           SIJ(I,J)=ANUM/DENOM
           ALTEM=SIJ(I,J)*Q(I)+(1._dp-SIJ(I,J))*QTI-QS(J)
           ALTEM=ALTEM-(BF2-1._dp)*CWAT
          END IF
          IF(SIJ(I,J).GT.0.0_dp.AND.SIJ(I,J).LT.0.9_dp)THEN
           QENT(I,J)=SIJ(I,J)*Q(I)+(1._dp-SIJ(I,J))*QTI
           UENT(I,J)=SIJ(I,J)*U(I)+(1._dp-SIJ(I,J))*U(NK)
           VENT(I,J)=SIJ(I,J)*V(I)+(1._dp-SIJ(I,J))*V(NK)
           DO K=1,NTRA
            TRAENT(I,J,K)=SIJ(I,J)*TRA(I,K)+(1._dp-SIJ(I,J))*TRA(NK,K)
           END DO
           ELIJ(I,J)=ALTEM
           ELIJ(I,J)=MAX(0.0_dp,ELIJ(I,J))
           MENT(I,J)=M(I)/(1._dp-SIJ(I,J))
           NENT(I)=NENT(I)+1
          END IF
          SIJ(I,J)=MAX(0.0_dp,SIJ(I,J))
          SIJ(I,J)=MIN(1.0_dp,SIJ(I,J))
  160    CONTINUE
!
!   ***   IF NO AIR CAN ENTRAIN AT LEVEL I ASSUME THAT UPDRAFT DETRAINS  ***
!   ***   AT THAT LEVEL AND CALCULATE DETRAINED AIR FLUX AND PROPERTIES  ***
!
         IF(NENT(I).EQ.0)THEN
          MENT(I,I)=M(I)
          QENT(I,I)=Q(NK)-EP(I)*CLW(I)
          UENT(I,I)=U(NK)
          VENT(I,I)=V(NK)
          DO J=1,NTRA
           TRAENT(I,I,J)=TRA(NK,J)
          END DO
          ELIJ(I,I)=CLW(I)
          SIJ(I,I)=1.0_dp
         END IF 
  170   CONTINUE
        SIJ(INB,INB)=1.0_dp
!
!   ***  NORMALIZE ENTRAINED AIR MASS FLUXES TO REPRESENT EQUAL  ***
!   ***              PROBABILITIES OF MIXING                     ***
!
        DO 200 I=ICB+1,INB
        IF(NENT(I).NE.0)THEN
         QP1=Q(NK)-EP(I)*CLW(I)
         ANUM=H(I)-HP(I)-LV(I)*(QP1-QS(I))
         DENOM=H(I)-HP(I)+LV(I)*(Q(I)-QP1)
         IF(ABS(DENOM).LT.0.01_dp)DENOM=0.01_dp
         SCRIT=ANUM/DENOM
         ALT=QP1-QS(I)+SCRIT*(Q(I)-QP1)
         IF(ALT.LT.0.0_dp)SCRIT=1.0_dp
         SCRIT=MAX(SCRIT,0.0_dp)
         ASIJ=0.0_dp
         SMIN=1.0_dp
         DO 175 J=ICB,INB
          IF(SIJ(I,J).GT.0.0_dp.AND.SIJ(I,J).LT.0.9_dp)THEN
           IF(J.GT.I)THEN
            SMID=MIN(SIJ(I,J),SCRIT)
            SJMAX=SMID
            SJMIN=SMID
            IF(SMID.LT.SMIN.AND.SIJ(I,J+1).LT.SMID)THEN
             SMIN=SMID
             SJMAX=MIN(SIJ(I,J+1),SIJ(I,J),SCRIT)
             SJMIN=MAX(SIJ(I,J-1),SIJ(I,J))
             SJMIN=MIN(SJMIN,SCRIT)
            END IF
           ELSE
            SJMAX=MAX(SIJ(I,J+1),SCRIT)
            SMID=MAX(SIJ(I,J),SCRIT)
            SJMIN=0.0_dp
            IF(J.GT.1)SJMIN=SIJ(I,J-1)
            SJMIN=MAX(SJMIN,SCRIT)
           END IF
           DELP=ABS(SJMAX-SMID)
           DELM=ABS(SJMIN-SMID)
           ASIJ=ASIJ+(DELP+DELM)*(PH(J)-PH(J+1))
           MENT(I,J)=MENT(I,J)*(DELP+DELM)*(PH(J)-PH(J+1))
          END IF
  175    CONTINUE
         ASIJ=MAX(1.0E-21_dp,ASIJ)
         ASIJ=1.0_dp/ASIJ
         DO 180 J=ICB,INB
          MENT(I,J)=MENT(I,J)*ASIJ
  180    CONTINUE
         BSUM=0.0_dp
         DO 190 J=ICB,INB
          BSUM=BSUM+MENT(I,J)
  190    CONTINUE
         IF(BSUM.LT.1.0E-18_dp)THEN
          NENT(I)=0
          MENT(I,I)=M(I)
          QENT(I,I)=Q(NK)-EP(I)*CLW(I)
          UENT(I,I)=U(NK)
          VENT(I,I)=V(NK)
          DO J=1,NTRA
           TRAENT(I,I,J)=TRA(NK,J)
          END DO
          ELIJ(I,I)=CLW(I)
          SIJ(I,I)=1.0_dp
         END IF
        END IF
  200   CONTINUE
!
!   ***  CHECK WHETHER EP(INB)=0, IF SO, SKIP PRECIPITATING    ***
!   ***             DOWNDRAFT CALCULATION                      ***
!
        IF(EP(INB).LT.0.0001_dp)GOTO 405
!
!   ***  INTEGRATE LIQUID WATER EQUATION TO FIND CONDENSED WATER   ***
!   ***                AND CONDENSED WATER FLUX                    ***
!
        JTT=2
!
!    ***                    BEGIN DOWNDRAFT LOOP                    ***
!
        DO 400 I=INB,1,-1
!
!    ***              CALCULATE DETRAINED PRECIPITATION             ***
!
        WDTRAIN=G*EP(I)*M(I)*CLW(I)
        IF(I.GT.1)THEN
         DO 320 J=1,I-1
         AWAT=ELIJ(J,I)-(1.-EP(I))*CLW(I)
         AWAT=MAX(0.0_dp,AWAT)
  320    WDTRAIN=WDTRAIN+G*AWAT*MENT(J,I)
        END IF
!
!    ***    FIND RAIN WATER AND EVAPORATION USING PROVISIONAL   ***
!    ***              ESTIMATES OF QP(I)AND QP(I-1)             ***
!     
!
!  ***  Value of terminal velocity and coefficient of evaporation for snow   ***
! 
        COEFF=COEFFS
        WT(I)=OMTSNOW
!      
!  ***  Value of terminal velocity and coefficient of evaporation for rain   ***
!
        IF(T(I).GT.273.0_dp)THEN
         COEFF=COEFFR
         WT(I)=OMTRAIN
        END IF
        QSM=0.5_dp*(Q(I)+QP(I+1))
        AFAC=COEFF*PH(I)*(QS(I)-QSM)/(1.0E4_dp+2.0E3_dp*PH(I)*QS(I))
        AFAC=MAX(AFAC,0.0_dp)
        SIGT=SIGP(I)
        SIGT=MAX(0.0_dp,SIGT)
        SIGT=MIN(1.0_dp,SIGT)
        B6=100._dp*(PH(I)-PH(I+1))*SIGT*AFAC/WT(I)
        C6=(WATER(I+1)*WT(I+1)+WDTRAIN/SIGD)/WT(I)
        REVAP=0.5_dp*(-B6+SQRT(B6*B6+4._dp*C6))
        EVAP(I)=SIGT*AFAC*REVAP
        WATER(I)=REVAP*REVAP
!
!    ***  CALCULATE PRECIPITATING DOWNDRAFT MASS FLUX UNDER     ***
!    ***              HYDROSTATIC APPROXIMATION                 ***
!   
        IF(I.EQ.1)GOTO 360
        DHDP=(H(I)-H(I-1))/(P(I-1)-P(I))
        DHDP=MAX(DHDP,10.0_dp)
        MP(I)=100._dp*GINV*LV(I)*SIGD*EVAP(I)/DHDP
        MP(I)=MAX(MP(I),0.0_dp)
!
!   ***   ADD SMALL AMOUNT OF INERTIA TO DOWNDRAFT              ***
!
        FAC=20.0_dp/(PH(I-1)-PH(I))
        MP(I)=(FAC*MP(I+1)+MP(I))/(1._dp+FAC)
!   
!    ***      FORCE MP TO DECREASE LINEARLY TO ZERO                 ***
!    ***      BETWEEN ABOUT 950 MB AND THE SURFACE                  ***
!
          IF(P(I).GT.(0.949_dp*P(1)))THEN
           JTT=MAX(JTT,I)
           MP(I)=MP(JTT)*(P(1)-P(I))/(P(1)-P(JTT))
          END IF              
  360   CONTINUE
!
!    ***       FIND MIXING RATIO OF PRECIPITATING DOWNDRAFT     ***
!
        IF(I.EQ.INB)GOTO 400
        IF(I.EQ.1)THEN
         QSTM=QS(1)
        ELSE
         QSTM=QS(I-1)
        END IF
        IF(MP(I).GT.MP(I+1))THEN
          RAT=MP(I+1)/MP(I)
          QP(I)=QP(I+1)*RAT+Q(I)*(1.0_dp-RAT)+100._dp*GINV*             &
                SIGD*(PH(I)-PH(I+1))*(EVAP(I)/MP(I))
          UP(I)=UP(I+1)*RAT+U(I)*(1._dp-RAT)
          VP(I)=VP(I+1)*RAT+V(I)*(1._dp-RAT)
          DO J=1,NTRA
           TRAP(I,J)=TRAP(I+1,J)*RAT+TRAP(I,J)*(1._dp-RAT)
          END DO
         ELSE
          IF(MP(I+1).GT.0.0_dp)THEN
            QP(I)=(GZ(I+1)-GZ(I)+QP(I+1)*(LV(I+1)+T(I+1)*(              &
                  CL-CPD))+CPD*(T(I+1)-T(I)))/(LV(I)+T(I)*(CL-CPD))
            UP(I)=UP(I+1)
            VP(I)=VP(I+1)
            DO J=1,NTRA
             TRAP(I,J)=TRAP(I+1,J)
            END DO
          END IF
        END IF
        QP(I)=MIN(QP(I),QSTM)
        QP(I)=MAX(QP(I),0.0_dp)
  400   CONTINUE
!
!   ***  CALCULATE SURFACE PRECIPITATION IN MM/DAY     ***
!
        DO I=INB,1,-1
          CV_PREC(I) = WT(I)*SIGD*WATER(I)*3600._dp*24000._dp/(ROWL*G)
          CV_FORM(I) = MAX((CV_PREC(I) - CV_PREC(i+1)),0._dp)
        END DO
        PRECIP=PRECIP+WT(1)*SIGD*WATER(1)*3600._dp*24000._dp/(ROWL*G)
        
!
  405   CONTINUE
!
!   ***  CALCULATE DOWNDRAFT VELOCITY SCALE AND SURFACE TEMPERATURE AND  ***
!   ***                    WATER VAPOR FLUCTUATIONS                      ***
!
      WD=BETA*ABS(MP(ICB))*0.01_dp*RD*T(ICB)/(SIGD*P(ICB))
      QPRIME=0.5_dp*(QP(1)-Q(1))
      TPRIME=LV0*QPRIME/CPD
!
!   ***  CALCULATE TENDENCIES OF LOWEST LEVEL POTENTIAL TEMPERATURE  ***
!   ***                      AND MIXING RATIO                        ***
!
        DPINV=0.01_dp/(PH(1)-PH(2))
        AM=0.0_dp
        IF(NK.EQ.1)THEN
         DO 410 K=2,INB
  410    AM=AM+M(K)
        END IF
        IF((2._dp*G*DPINV*AM).GE.DELTI)IFLAG=4
        FT(1)=FT(1)+G*DPINV*AM*(T(2)-T(1)+(GZ(2)-GZ(1))/CPN(1))
        FT(1)=FT(1)-LVCP(1)*SIGD*EVAP(1)
        FT(1)=FT(1)+SIGD*WT(2)*(CL-CPD)*WATER(2)*(T(2)-T(1))*DPINV/CPN(1)
        FQ(1)=FQ(1)+G*MP(2)*(QP(2)-Q(1))*DPINV+SIGD*EVAP(1)
        FQ(1)=FQ(1)+G*AM*(Q(2)-Q(1))*DPINV
        FU(1)=FU(1)+G*DPINV*(MP(2)*(UP(2)-U(1))+AM*(U(2)-U(1)))
        FV(1)=FV(1)+G*DPINV*(MP(2)*(VP(2)-V(1))+AM*(V(2)-V(1)))
        DO J=1,NTRA
         FTRA(1,J)=FTRA(1,J)+G*DPINV*(MP(2)*(TRAP(2,J)-TRA(1,J))+         &
                   AM*(TRA(2,J)-TRA(1,J)))
        END DO
        AMDE=0.0_dp
        DO 415 J=2,INB
         FQ(1)=FQ(1)+G*DPINV*MENT(J,1)*(QENT(J,1)-Q(1))
         FU(1)=FU(1)+G*DPINV*MENT(J,1)*(UENT(J,1)-U(1))
         FV(1)=FV(1)+G*DPINV*MENT(J,1)*(VENT(J,1)-V(1))
         DO K=1,NTRA
          FTRA(1,K)=FTRA(1,K)+G*DPINV*MENT(J,1)*(TRAENT(J,1,K)-TRA(1,K))
         END DO
  415   CONTINUE
!
!   ***  CALCULATE TENDENCIES OF POTENTIAL TEMPERATURE AND MIXING RATIO  ***
!   ***               AT LEVELS ABOVE THE LOWEST LEVEL                   ***
!
!   ***  FIRST FIND THE NET SATURATED UPDRAFT AND DOWNDRAFT MASS FLUXES  ***
!   ***                      THROUGH EACH LEVEL                          ***
!
        PUMF(1) = AM
        DO 500 I=2,INB
        DPINV=0.01_dp/(PH(I)-PH(I+1))
        CPINV=1.0_dp/CPN(I)
        AMP1=0.0_dp
        AD=0.0_dp
        IF(I.GE.NK)THEN
         DO 440 K=I+1,INB+1
  440    AMP1=AMP1+M(K)
        END IF
        DO 450 K=1,I
        DO 450 J=I+1,INB+1
         AMP1=AMP1+MENT(K,J)
  450   CONTINUE
        ! mz_ht_20070220+
        PUMF(I) = AMP1

        DO K=I,INB+1
          PUER(I) = PUER(I) + MENT(I,K)
        ENDDO
        DO K=1,I
          PUDR(I) = PUDR(I) + MENT(K,I)
        ENDDO
        ! mz_ht_20070220-
        IF((2._dp*G*DPINV*AMP1).GE.DELTI)IFLAG=4
        DO 470 K=1,I-1
        DO 470 J=I,INB
         AD=AD+MENT(J,K)
  470   CONTINUE
        ! mz_ht_20070220+         
        PDMF(I) = AD
        DO K=1,I
          PDER(I) = PDER(I) + MENT(I,K)
        ENDDO
        DO K=I,INB+1
          PUDR(I) = PUDR(I) + MENT(K,I)
        ENDDO
        ! mz_ht_20070220-
        FT(I)=FT(I)+G*DPINV*(AMP1*(T(I+1)-T(I)+(GZ(I+1)-GZ(I))*     &
              CPINV)-AD*(T(I)-T(I-1)+(GZ(I)-GZ(I-1))*CPINV))        &
              -SIGD*LVCP(I)*EVAP(I)
        FT(I)=FT(I)+G*DPINV*MENT(I,I)*(HP(I)-H(I)+                  &
              T(I)*(CPV-CPD)*(Q(I)-QENT(I,I)))*CPINV
        FT(I)=FT(I)+SIGD*WT(I+1)*(CL-CPD)*WATER(I+1)*               &
              (T(I+1)-T(I))*DPINV*CPINV
        FQ(I)=FQ(I)+G*DPINV*(AMP1*(Q(I+1)-Q(I))-AD*(Q(I)-Q(I-1)))
        FU(I)=FU(I)+G*DPINV*(AMP1*(U(I+1)-U(I))-AD*(U(I)-U(I-1)))
        FV(I)=FV(I)+G*DPINV*(AMP1*(V(I+1)-V(I))-AD*(V(I)-V(I-1)))
        DO K=1,NTRA
         FTRA(I,K)=FTRA(I,K)+G*DPINV*(AMP1*(TRA(I+1,K)-             &
                   TRA(I,K))-AD*(TRA(I,K)-TRA(I-1,K)))
        END DO
        DO 480 K=1,I-1
         AWAT=ELIJ(K,I)-(1._dp-EP(I))*CLW(I)
         AWAT=MAX(AWAT,0.0_dp)
         FQ(I)=FQ(I)+G*DPINV*MENT(K,I)*(QENT(K,I)-AWAT-Q(I))
         FU(I)=FU(I)+G*DPINV*MENT(K,I)*(UENT(K,I)-U(I))
         FV(I)=FV(I)+G*DPINV*MENT(K,I)*(VENT(K,I)-V(I))
         DO J=1,NTRA
          FTRA(I,J)=FTRA(I,J)+G*DPINV*MENT(K,I)*(TRAENT(K,I,J)-TRA(I,J))
         END DO
  480   CONTINUE
        DO 490 K=I,INB
         FQ(I)=FQ(I)+G*DPINV*MENT(K,I)*(QENT(K,I)-Q(I))
         FU(I)=FU(I)+G*DPINV*MENT(K,I)*(UENT(K,I)-U(I))
         FV(I)=FV(I)+G*DPINV*MENT(K,I)*(VENT(K,I)-V(I))
         DO J=1,NTRA
          FTRA(I,J)=FTRA(I,J)+G*DPINV*MENT(K,I)*(TRAENT(K,I,J)-TRA(I,J))
         END DO
  490   CONTINUE
        FQ(I)=FQ(I)+SIGD*EVAP(I)+G*(MP(I+1)*                                &
              (QP(I+1)-Q(I))-MP(I)*(QP(I)-Q(I-1)))*DPINV
        FU(I)=FU(I)+G*(MP(I+1)*(UP(I+1)-U(I))-MP(I)*(UP(I)-U(I-1)))*DPINV
        FV(I)=FV(I)+G*(MP(I+1)*(VP(I+1)-V(I))-MP(I)*(VP(I)-V(I-1)))*DPINV
        DO J=1,NTRA
         FTRA(I,J)=FTRA(I,J)+G*DPINV*(MP(I+1)*(TRAP(I+1,J)-TRA(I,J))-       &
                   MP(I)*(TRAP(I,J)-TRA(I-1,J)))
        END DO
  500   CONTINUE
!
!   *** Adjust tendencies at top of convection layer to reflect  ***
!   ***       actual position of the level zero CAPE             ***
!
        FQOLD=FQ(INB)
        FQ(INB)=FQ(INB)*(1._dp-FRAC)
        FQ(INB-1)=FQ(INB-1)+FRAC*FQOLD*((PH(INB)-PH(INB+1))/                &
                  (PH(INB-1)-PH(INB)))*LV(INB)/LV(INB-1)
        FTOLD=FT(INB)
        FT(INB)=FT(INB)*(1._dp-FRAC)
        FT(INB-1)=FT(INB-1)+FRAC*FTOLD*((PH(INB)-PH(INB+1))/                &
                  (PH(INB-1)-PH(INB)))*CPN(INB)/CPN(INB-1)
        FUOLD=FU(INB)
        FU(INB)=FU(INB)*(1._dp-FRAC)
        FU(INB-1)=FU(INB-1)+FRAC*FUOLD*((PH(INB)-PH(INB+1))/                &
                  (PH(INB-1)-PH(INB)))
        FVOLD=FV(INB)
        FV(INB)=FV(INB)*(1._dp-FRAC)
        FV(INB-1)=FV(INB-1)+FRAC*FVOLD*((PH(INB)-PH(INB+1))/                &
                  (PH(INB-1)-PH(INB)))
        DO K=1,NTRA
         FTRAOLD=FTRA(INB,K)
         FTRA(INB,K)=FTRA(INB,K)*(1._dp-FRAC)
         FTRA(INB-1,K)=FTRA(INB-1,K)+FRAC*FTRAOLD*(PH(INB)-PH(INB+1))/      &
                      (PH(INB-1)-PH(INB))
        END DO
!
!   ***   Very slightly adjust tendencies to force exact   ***
!   ***     enthalpy, momentum and tracer conservation     ***
!
        ENTS=0.0_dp
        UAV=0.0_dp
        VAV=0.0_dp
        DO 680 I=1,INB
         ENTS=ENTS+(CPN(I)*FT(I)+LV(I)*FQ(I))*(PH(I)-PH(I+1))      
         UAV=UAV+FU(I)*(PH(I)-PH(I+1))
         VAV=VAV+FV(I)*(PH(I)-PH(I+1))
  680      CONTINUE
        ENTS=ENTS/(PH(1)-PH(INB+1))
        UAV=UAV/(PH(1)-PH(INB+1))
        VAV=VAV/(PH(1)-PH(INB+1))
        DO 640 I=1,INB
         FT(I)=FT(I)-ENTS/CPN(I)
         FU(I)=(1._dp-CU)*(FU(I)-UAV)
         FV(I)=(1._dp-CU)*(FV(I)-VAV)
  640      CONTINUE
        DO 700 K=1,NTRA
         TRAAV=0.0_dp
         DO 690 I=1,INB
          TRAAV=TRAAV+FTRA(I,K)*(PH(I)-PH(I+1))
  690    CONTINUE
         TRAAV=TRAAV/(PH(1)-PH(INB+1))
         DO 695 I=1,INB
          FTRA(I,K)=FTRA(I,K)-TRAAV
  695    CONTINUE
  700      CONTINUE
!
!   ***           RETURN           ***
!
!          print*, "final return"
        RETURN
!
      END SUBROUTINE CONVECT_EMANUEL
!
! ---------------------------------------------------------------------------
!
      SUBROUTINE TLIFT(P,T,Q,QS,GZ,ICB,NK,TVP,TPK,CLW,ND,NL,KK)
        ! mz_ht_20070214+
        ! explicit variable declaration
        INTEGER  :: ND, NL, NK, KK, ICB
        INTEGER  :: I, J, NST, NSB
        REAL(dp) :: RG, ES, DENOM, TC, AHG, AH0, S
        REAL(dp) :: TG, QG, ALV
        REAL(dp) :: CPP, CPINV
        ! mz_ht_20070214-


        REAL(dp) :: GZ(ND),TPK(ND),CLW(ND),P(ND)
        REAL(dp) :: T(ND),Q(ND),QS(ND),TVP(ND)


!
!!$!   ***   ASSIGN VALUES OF THERMODYNAMIC CONSTANTS     ***


!!$!   ***    already done in module header   ***  ! mz_ht_20070214

!!$        CPD=1005.7
!!$        CPV=1870.0
!!$        CL=2500.0
!!$        RV=461.5
!!$        RD=287.04
!!$        LV0=2.501E6
!!$!
!!$        CPVMCL=CL-CPV
!!$        EPS=RD/RV
!!$        EPSI=1./EPS
!
!   ***  CALCULATE CERTAIN PARCEL QUANTITIES, INCLUDING STATIC ENERGY   ***
!
        AH0=(CPD*(1._dp-Q(NK))+CL*Q(NK))*T(NK)+Q(NK)*(LV0-CPVMCL*(        &
            T(NK)-TMELT))+GZ(NK)
        CPP=CPD*(1._dp-Q(NK))+Q(NK)*CPV
        CPINV=1._dp/CPP
!
        IF(KK.EQ.1)THEN
!
!   ***   CALCULATE LIFTED PARCEL QUANTITIES BELOW CLOUD BASE   ***
!
        DO 50 I=1,ICB-1
         CLW(I)=0.0_dp
   50   CONTINUE
        DO 100 I=NK,ICB-1
         TPK(I)=T(NK)-(GZ(I)-GZ(NK))*CPINV
         TVP(I)=TPK(I)*(1._dp+Q(NK)*EPSI)
  100   CONTINUE
        END IF
!
!    ***  FIND LIFTED PARCEL QUANTITIES ABOVE CLOUD BASE    ***
!
        NST=ICB
        NSB=ICB
        IF(KK.EQ.2)THEN  
         NST=NL
         NSB=ICB+1
        END IF
        DO 300 I=NSB,NST
         TG=T(I)
         QG=QS(I)
         ALV=LV0-CPVMCL*(T(I)-TMELT)
         DO 200 J=1,2
          S=CPD+ALV*ALV*QG/(RV*T(I)*T(I))
          S=1._dp/S
          AHG=CPD*TG+(CL-CPD)*Q(NK)*T(I)+ALV*QG+GZ(I)
          TG=TG+S*(AH0-AHG)
          TG=MAX(TG,35.0_dp)
          TC=TG-TMELT
          DENOM=243.5_dp+TC
          IF(TC.GE.0.0_dp)THEN  
           ES=6.112_dp*EXP(17.67_dp*TC/DENOM)
          ELSE  
           ES=EXP(23.33086_dp-6111.72784_dp/TG+0.15215_dp*LOG(TG))
          END IF  
          QG=EPS*ES/(P(I)-ES*(1._dp-EPS))
  200    CONTINUE
         TPK(I)=(AH0-(CL-CPD)*Q(NK)*T(I)-GZ(I)-ALV*QG)/CPD
         CLW(I)=Q(NK)-QG
         CLW(I)=MAX(0.0_dp,CLW(I))
         RG=QG/(1._dp-Q(NK))
         TVP(I)=TPK(I)*(1._dp+RG*EPSI)
  300   CONTINUE
        RETURN
      END SUBROUTINE TLIFT



!----------------------------------------------------------------------------
END MODULE MESSY_CONVECT_EMANUEL
