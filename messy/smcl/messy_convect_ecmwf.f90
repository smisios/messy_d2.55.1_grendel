MODULE MESSY_CONVECT_ECMWF

!    This module contains the source routines for the Convection 
!    scheme from ECMWF (Version 29r1)

!    Author of this module: 
!    little adjustements due to implementation into the MESSy system
!    H. Tost,    MPI - Chemie, Mainz Feb. 2005

!    Original code as described in the routines

  USE MESSY_MAIN_CONSTANTS_MEM,      ONLY: dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC ::  CUININ, CUBASEN, CUASCN, CUDLFSN, CUDDRAFN, CUFLXN, &
             CUDTDQN, CUDUDV, CUCTRACER

CONTAINS


!==============================================================================

SUBROUTINE CUININ &
 & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
 & PTEN,     PQEN,     PQSEN,    PUEN,     PVEN,&
 & PVERVEL,  PGEO,     PAPH,     PAP,&
 & KLWMIN,   KLAB,&
 & PTENH,    PQENH,    PQSENH,   PGEOH,&
 & PTU,      PQU,      PTD,      PQD,&
 & PUU,      PVU,      PUD,      PVD,&
 & PLU  )  

!          M.TIEDTKE         E.C.M.W.F.     12/89

!          PURPOSE
!          -------

!          THIS ROUTINE INTERPOLATES LARGE-SCALE FIELDS OF T,Q ETC.
!          TO HALF LEVELS (I.E. GRID FOR MASSFLUX SCHEME),
!          DETERMINES LEVEL OF MAXIMUM VERTICAL VELOCITY
!          AND INITIALIZES VALUES FOR UPDRAFTS AND DOWNDRAFTS

!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.

!          METHOD.
!          --------
!          FOR EXTRAPOLATION TO HALF LEVELS SEE TIEDTKE(1989)

!     PARAMETER     DESCRIPTION                                   UNITS 
!     ---------     -----------                                   ----- 
!     INPUT PARAMETERS (INTEGER): 

!    *KIDIA*        START POINT 
!    *KFDIA*        END POINT 
!    *KLON*         NUMBER OF GRID POINTS PER PACKET 
!    *KTDIA*        START OF THE VERTICAL LOOP 
!    *KLEV*         NUMBER OF LEVELS 

!    INPUT PARAMETERS (REAL): 

!    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K 
!    *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1)  KG/KG 
!    *PQSEN*        ENVIRONMENT SPEC. SATURATION HUMIDITY (T+1)   KG/KG 
!    *PUEN*         PROVISIONAL ENVIRONMENT U-VELOCITY (T+1)       M/S
!    *PVEN*         PROVISIONAL ENVIRONMENT V-VELOCITY (T+1)       M/S
!    *PVERVEL*      VERTICAL VELOCITY                             PA/S
!    *PGEO*         GEOPOTENTIAL                                  M2/S2
!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS             PA
!    *PAP*          PROVISIONAL PRESSURE ON FULL LEVELS             PA

!    OUTPUT PARAMETERS (INTEGER):

!    *KLWMIN*       LEVEL OF MAXIMUM VERTICAL VELOCITY 
!    *KLAB*         FLAG KLAB=1 FOR SUBCLOUD LEVELS
!                        KLAB=2 FOR CONDENSATION LEVEL

!    OUTPUT PARAMETERS (REAL):

!    *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS         K
!    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS    KG/KG
!    *PQSENH*       ENV. SPEC. SATURATION HUMIDITY (T+1)
!                   ON HALF LEVELS                              KG/KG
!    *PTU*          TEMPERATURE IN UPDRAFTS                       K
!    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                  KG/KG
!    *PTD*          TEMPERATURE IN DOWNDRAFTS                     K
!    *PQU*          SPEC. HUMIDITY IN DOWNDRAFTS                KG/KG
!    *PUU*          U-VELOCITY IN UPDRAFTS                       M/S
!    *PVU*          V-VELOCITY IN UPDRAFTS                       M/S
!    *PUD*          U-VELOCITY IN DOWNDRAFTS                     M/S
!    *PVD*          V-VELOCITY IN DOWNDRAFTS                     M/S
!    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS            KG/KG

!          EXTERNALS
!          ---------
!          *CUADJTQ* TO SPECIFY QS AT HALF LEVELS

!          MODIFICATIONS
!          -------------
!             92-09-21 : Update to Cy44      J.-J. MORCRETTE
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!----------------------------------------------------------------------

USE MESSY_CONVECT_ECMWF_PARAM   , ONLY : RCPD,              &
                                         LPHYLIN,           &
                                         LHOOK, DR_HOOK          !from YOMHOOK

IMPLICIT NONE

INTEGER    ,INTENT(IN)    :: KLON 
INTEGER    ,INTENT(IN)    :: KLEV 
INTEGER    ,INTENT(IN)    :: KIDIA 
INTEGER    ,INTENT(IN)    :: KFDIA 
INTEGER                   :: KTDIA ! Argument NOT used
REAL(dp)   ,INTENT(IN)    :: PTEN(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PQEN(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PQSEN(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PUEN(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PVEN(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PVERVEL(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PGEO(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) 
REAL(dp)                  :: PAP(KLON,KLEV) ! Argument NOT used
INTEGER    ,INTENT(OUT)   :: KLWMIN(KLON) 
INTEGER    ,INTENT(OUT)   :: KLAB(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PTENH(KLON,KLEV) 
REAL(dp)   ,INTENT(OUT)   :: PQENH(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PQSENH(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PGEOH(KLON,KLEV+1) 
REAL(dp)   ,INTENT(OUT)   :: PTU(KLON,KLEV) 
REAL(dp)   ,INTENT(OUT)   :: PQU(KLON,KLEV) 
REAL(dp)   ,INTENT(OUT)   :: PTD(KLON,KLEV) 
REAL(dp)   ,INTENT(OUT)   :: PQD(KLON,KLEV) 
REAL(dp)   ,INTENT(OUT)   :: PUU(KLON,KLEV) 
REAL(dp)   ,INTENT(OUT)   :: PVU(KLON,KLEV) 
REAL(dp)   ,INTENT(OUT)   :: PUD(KLON,KLEV) 
REAL(dp)   ,INTENT(OUT)   :: PVD(KLON,KLEV) 
REAL(dp)   ,INTENT(OUT)   :: PLU(KLON,KLEV) 
REAL(dp)                  ::     ZWMAX(KLON)
REAL(dp)                  ::     ZPH(KLON)
LOGICAL                   ::  LLFLAG(KLON)

INTEGER                   :: icall, IK, JK, JL

REAL(dp)                  :: ZALFA, ZZS
REAL(dp)                  :: ZHOOK_HANDLE

INTRINSIC :: MAX, MIN, LOG

!----------------------------------------------------------------------

!*    1.           SPECIFY LARGE SCALE PARAMETERS AT HALF LEVELS
!*                 ADJUST TEMPERATURE FIELDS IF STATICLY UNSTABLE
!*                 FIND LEVEL OF MAXIMUM VERTICAL VELOCITY
!                  ----------------------------------------------

IF (LHOOK) CALL DR_HOOK('CUININ',0,ZHOOK_HANDLE)
ZALFA=LOG(2.0_dp)
DO JK=2,KLEV
  DO JL=KIDIA,KFDIA
    PTENH(JL,JK)=(MAX(RCPD*PTEN(JL,JK-1)+PGEO(JL,JK-1),&
     & RCPD*PTEN(JL,JK)+PGEO(JL,JK))-PGEOH(JL,JK))/RCPD  
    PQENH(JL,JK)=PQEN(JL,JK-1)
    PQSENH(JL,JK)=PQSEN(JL,JK-1)
    ZPH(JL)=PAPH(JL,JK)
    LLFLAG(JL)=.TRUE.
  ENDDO

!orig   IF(JK.GE.KLEV-1) GO TO 130
  IF(JK >= KLEV-1) CYCLE
  IK=JK
  if(lphylin)then
    icall=0
    CALL CUADJTQS &
     & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
     & IK,&
     & ZPH,      PTENH,    PQSENH,   LLFLAG,   ICALL)  
  else
    ICALL=3
    CALL CUADJTQ &
     & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
     & IK,&
     & ZPH,      PTENH,    PQSENH,   LLFLAG,   ICALL)  
  ENDIF

  DO JL=KIDIA,KFDIA
    PQENH(JL,JK)=MIN(PQEN(JL,JK-1),PQSEN(JL,JK-1))&
     & +(PQSENH(JL,JK)-PQSEN(JL,JK-1))  
    PQENH(JL,JK)=MAX(PQENH(JL,JK),0.0_dp)
  ENDDO
!orig  130   continue
ENDDO

DO JL=KIDIA,KFDIA
  PTENH(JL,KLEV)=(RCPD*PTEN(JL,KLEV)+PGEO(JL,KLEV)-PGEOH(JL,KLEV))/RCPD
  PQENH(JL,KLEV)=PQEN(JL,KLEV)
  PTENH(JL,1)=PTEN(JL,1)
  PQENH(JL,1)=PQEN(JL,1)
  KLWMIN(JL)=KLEV
  ZWMAX(JL)=0.0_dp
ENDDO

DO JK=KLEV-1,2,-1
  DO JL=KIDIA,KFDIA
    ZZS=MAX(RCPD*PTENH(JL,JK)+PGEOH(JL,JK),&
     & RCPD*PTENH(JL,JK+1)+PGEOH(JL,JK+1))  
    PTENH(JL,JK)=(ZZS-PGEOH(JL,JK))/RCPD
  ENDDO
ENDDO

DO JK=KLEV,3,-1
!DIR$ IVDEP
!OCL NOVREC
  DO JL=KIDIA,KFDIA
    IF(PVERVEL(JL,JK) < ZWMAX(JL)) THEN
      ZWMAX(JL)=PVERVEL(JL,JK)
      KLWMIN(JL)=JK
    ENDIF
  ENDDO
ENDDO

!-----------------------------------------------------------------------

!*    2.0          INITIALIZE VALUES FOR UPDRAFTS AND DOWNDRAFTS
!*                 ---------------------------------------------

DO JK=1,KLEV
  IK=JK-1
  IF(JK == 1) IK=1
  DO JL=KIDIA,KFDIA
    PTU(JL,JK)=PTENH(JL,JK)
    PTD(JL,JK)=PTENH(JL,JK)
    PQU(JL,JK)=PQENH(JL,JK)
    PQD(JL,JK)=PQENH(JL,JK)
    PLU(JL,JK)=0.0_dp
    PUU(JL,JK)=PUEN(JL,IK)
    PUD(JL,JK)=PUEN(JL,IK)
    PVU(JL,JK)=PVEN(JL,IK)
    PVD(JL,JK)=PVEN(JL,IK)
    KLAB(JL,JK)=0
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('CUININ',1,ZHOOK_HANDLE)
END SUBROUTINE CUININ

!==============================================================================

SUBROUTINE CUADJTQ &
 & (KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
 & KK,&
 & PSP,      PT,       PQ,       LDFLAG,   KCALL)  

!          M.TIEDTKE         E.C.M.W.F.     12/89

!          MODIFICATIONS
!          -------------
!          D.SALMOND         CRAY(UK))      12/8/91
!          J.J. MORCRETTE    ECMWF          92-09-18   Update to Cy44
!          J.F. MAHFOUF      ECMWF          96-06-11   Smoothing option
!          D.SALMOND & M.HAMRUD ECMWF       99-06-04   Optimisation
!          J.HAGUE                          03-01-13   MASS Vector Functions
!          J.HAGUE                          03-07-07   More MASS V.F.
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!          PURPOSE.
!          --------
!          TO PRODUCE T,Q AND L VALUES FOR CLOUD ASCENT

!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM SUBROUTINES:
!              *COND*     (T AND Q AT CONDENSATION LEVEL)
!              *CUBASE*   (T AND Q AT CONDENSATION LEVEL)
!              *CUASC*    (T AND Q AT CLOUD LEVELS)
!              *CUINI*    (ENVIRONMENTAL T AND QS VALUES AT HALF LEVELS)
!              *CUSTRAT*  (T AND Q AT CONDENSATION LEVEL)
!          INPUT ARE UNADJUSTED T AND Q VALUES,
!          IT RETURNS ADJUSTED VALUES OF T AND Q

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KTDIA*        START OF THE VERTICAL LOOP
!    *KLEV*         NUMBER OF LEVELS
!    *KK*           LEVEL
!    *KCALL*        DEFINES CALCULATION AS
!                      KCALL=0  ENV. T AND QS IN*CUINI*
!                      KCALL=1  CONDENSATION IN UPDRAFTS  (E.G. CUBASE, CUASC)
!                      KCALL=2  EVAPORATION IN DOWNDRAFTS (E.G. CUDLFS,CUDDRAF)

!     INPUT PARAMETERS (LOGICAL):

!    *LDLAND*       LAND-SEA MASK (.TRUE. FOR LAND POINTS)

!     INPUT PARAMETERS (REAL):

!    *PSP*          PRESSURE                                        PA

!     UPDATED PARAMETERS (REAL):

!    *PT*           TEMPERATURE                                     K
!    *PQ*           SPECIFIC HUMIDITY                             KG/KG

!          EXTERNALS   
!          ---------
!          3 LOOKUP TABLES ( TLUCUA, TLUCUB, TLUCUC )
!          FOR CONDENSATION CALCULATIONS.
!          THE TABLES ARE INITIALISED IN *SUPHEC*.

!----------------------------------------------------------------------

USE MESSY_CONVECT_ECMWF_PARAM   , ONLY : RETV, RLVTT, RLSTT, RTT,                 &
                                         LPHYLIN, RLPTRC, RLPAL1, RLPAL2,         & 
                                         R2ES, R3LES, R3IES, R4LES,               &
                                         R4IES, R5LES, R5IES, R5ALVCP, R5ALSCP,   &
                                         RALVDCP, RALSDCP, RTWAT, RTICE, RTICECU, &
                                         RTWAT_RTICE_R, RTWAT_RTICECU_R,          &
                                         N_VMASS,                                 &
                                         FOEALFA, FOEALFCU, FOEDEM, FOELDCPM,     &!functions
                                         FOELDCPMCU, FOEEWM, FOEEWMCU, FOEDEMCU,  &!functions
                                         VDIV, VEXP,                              &!subroutines (?)
                                         LHOOK, DR_HOOK                            !from YOMHOOK


IMPLICIT NONE

INTEGER,INTENT(IN)    :: KLON 
INTEGER,INTENT(IN)    :: KLEV 
INTEGER,INTENT(IN)    :: KIDIA 
INTEGER,INTENT(IN)    :: KFDIA 
INTEGER               :: KTDIA ! Argument NOT used
INTEGER,INTENT(IN)    :: KK 
REAL(dp)   ,INTENT(IN)    :: PSP(KLON) 
REAL(dp)   ,INTENT(INOUT) :: PT(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PQ(KLON,KLEV) 
LOGICAL           ,INTENT(IN)    :: LDFLAG(KLON) 
INTEGER,INTENT(IN)    :: KCALL 
INTEGER :: JL, JLEN

REAL(dp) :: ZTMP0(KFDIA-KIDIA+N_VMASS)
REAL(dp) :: ZTMP1(KFDIA-KIDIA+N_VMASS)
REAL(dp) :: ZTMP2(KFDIA-KIDIA+N_VMASS)
REAL(dp) :: ZTMP3(KFDIA-KIDIA+N_VMASS)
REAL(dp) :: ZTMP4(KFDIA-KIDIA+N_VMASS)
REAL(dp) :: ZTMP5(KFDIA-KIDIA+N_VMASS)
REAL(dp) :: ZTMP6(KFDIA-KIDIA+N_VMASS)

REAL(dp) :: Z1S, Z2S, ZCOND,ZCOND1, ZCOR, ZFOEEWI, ZFOEEWL,&
 & ZOEALFA, ZQMAX, ZQSAT, ZTARG, ZQP  


!     STATEMENT FUNCTIONS
!REAL(dp) :: FOEALFAJ,FOEDEMJ,FOELDCPMJ,FOEEWMJ

REAL(dp) :: MINJ, MAXJ, X, Y
REAL(dp) :: ZHOOK_HANDLE

INTRINSIC :: ABS, EXP, MAX, MIN, MOD, TANH

MINJ(X,Y) = Y - 0.5_dp*(ABS(X-Y)-(X-Y))
MAXJ(X,Y) = Y + 0.5_dp*(ABS(X-Y)+(X-Y))

!----------------------------------------------------------------------

!     1.           DEFINE CONSTANTS
!                  ----------------

IF (LHOOK) CALL DR_HOOK('CUADJTQ',0,ZHOOK_HANDLE)
IF(N_VMASS >  0) THEN
  JLEN=KFDIA-KIDIA+N_VMASS-MOD(KFDIA-KIDIA,N_VMASS)
  IF(KFDIA-KIDIA+1 /= JLEN) THEN 
    IF(KCALL == 3 .OR. KCALL == 5) THEN
      ZTMP1(KFDIA-KIDIA+2:JLEN)=1.0_dp
      ZTMP2(KFDIA-KIDIA+2:JLEN)=1.0_dp
      ZTMP3(KFDIA-KIDIA+2:JLEN)=1.0_dp
      ZTMP4(KFDIA-KIDIA+2:JLEN)=1.0_dp
    ENDIF
  ENDIF
ENDIF

ZQMAX=0.5_dp

!*********************************************
IF (.NOT.LPHYLIN) THEN
!*********************************************                 

!     2.           CALCULATE CONDENSATION AND ADJUST T AND Q ACCORDINGLY
!                  -----------------------------------------------------

  IF (KCALL == 1 ) THEN

!DIR$ IVDEP
!OCL NOVREC
    DO JL=KIDIA,KFDIA
      IF(LDFLAG(JL)) THEN
        ZQP    =1.0_dp/PSP(JL)
        ZQSAT=FOEEWMCU(PT(JL,KK))*ZQP    
        ZQSAT=MIN(0.5_dp,ZQSAT)
        ZCOR=1.0_dp-RETV*ZQSAT
        ZCOND=(PQ(JL,KK)*ZCOR**2-ZQSAT*ZCOR)/(ZCOR**2+ZQSAT*FOEDEMCU(PT(JL,KK)))
        ZCOND=MAX(ZCOND,0.0_dp)
        PT(JL,KK)=PT(JL,KK)+FOELDCPMCU(PT(JL,KK))*ZCOND
        PQ(JL,KK)=PQ(JL,KK)-ZCOND
        ZQSAT=FOEEWMCU(PT(JL,KK))*ZQP    
        ZQSAT=MINJ(0.5_dp,ZQSAT)
        ZCOR=1.0_dp-RETV*ZQSAT
        ZCOND1=(PQ(JL,KK)*ZCOR**2-ZQSAT*ZCOR)/(ZCOR**2+ZQSAT*FOEDEMCU(PT(JL,KK)))
        IF(ZCOND ==  0.0_dp)ZCOND1=0.0_dp
        PT(JL,KK)=PT(JL,KK)+FOELDCPMCU(PT(JL,KK))*ZCOND1
        PQ(JL,KK)=PQ(JL,KK)-ZCOND1
      ENDIF
    ENDDO

  ENDIF

  IF(KCALL == 2) THEN

!DIR$ IVDEP
!OCL NOVREC
    DO JL=KIDIA,KFDIA
      IF(LDFLAG(JL)) THEN
        ZQP    =1.0_dp/PSP(JL)
        ZQSAT=FOEEWMCU(PT(JL,KK))*ZQP    
        ZQSAT=MIN(0.5_dp,ZQSAT)
        ZCOR=1.0_dp/(1.0_dp-RETV  *ZQSAT)
        ZQSAT=ZQSAT*ZCOR
        ZCOND=(PQ(JL,KK)-ZQSAT)/(1.0_dp+ZQSAT*ZCOR*FOEDEMCU(PT(JL,KK)))
        ZCOND=MIN(ZCOND,0.0_dp)
        PT(JL,KK)=PT(JL,KK)+FOELDCPMCU(PT(JL,KK))*ZCOND
        PQ(JL,KK)=PQ(JL,KK)-ZCOND
        ZQSAT=FOEEWMCU(PT(JL,KK))*ZQP    
        ZQSAT=MIN(0.5_dp,ZQSAT)
        ZCOR=1.0_dp/(1.0_dp-RETV  *ZQSAT)
        ZQSAT=ZQSAT*ZCOR
        ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.0_dp+ZQSAT*ZCOR*FOEDEMCU(PT(JL,KK)))
        IF(ZCOND == 0.0_dp)ZCOND1=MIN(ZCOND1,0.0_dp)
        PT(JL,KK)=PT(JL,KK)+FOELDCPMCU(PT(JL,KK))*ZCOND1
        PQ(JL,KK)=PQ(JL,KK)-ZCOND1
      ENDIF
    ENDDO

  ENDIF

  IF(KCALL == 0) THEN

!DIR$ IVDEP
!OCL NOVREC
    DO JL=KIDIA,KFDIA
      ZQP    =1.0_dp/PSP(JL)
      ZQSAT=FOEEWM(PT(JL,KK))*ZQP    
      ZQSAT=MIN(0.5_dp,ZQSAT)
      ZCOR=1.0_dp/(1.0_dp-RETV  *ZQSAT)
      ZQSAT=ZQSAT*ZCOR
      ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.0_dp+ZQSAT*ZCOR*FOEDEM(PT(JL,KK)))
      PT(JL,KK)=PT(JL,KK)+FOELDCPM(PT(JL,KK))*ZCOND1
      PQ(JL,KK)=PQ(JL,KK)-ZCOND1
      ZQSAT=FOEEWM(PT(JL,KK))*ZQP    
      ZQSAT=MIN(0.5_dp,ZQSAT)
      ZCOR=1.0_dp/(1.0_dp-RETV  *ZQSAT)
      ZQSAT=ZQSAT*ZCOR
      ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.0_dp+ZQSAT*ZCOR*FOEDEM(PT(JL,KK)))
      PT(JL,KK)=PT(JL,KK)+FOELDCPM(PT(JL,KK))*ZCOND1
      PQ(JL,KK)=PQ(JL,KK)-ZCOND1
    ENDDO

  ENDIF

  IF(KCALL == 4 )THEN

!DIR$ IVDEP
!OCL NOVREC
    DO JL=KIDIA,KFDIA
      IF(LDFLAG(JL)) THEN
        ZQP    =1.0_dp/PSP(JL)
        ZQSAT=FOEEWM(PT(JL,KK))*ZQP    
        ZQSAT=MIN(0.5_dp,ZQSAT)
        ZCOR=1.0_dp/(1.0_dp-RETV  *ZQSAT)
        ZQSAT=ZQSAT*ZCOR
        ZCOND=(PQ(JL,KK)-ZQSAT)/(1.0_dp+ZQSAT*ZCOR*FOEDEM(PT(JL,KK)))
        PT(JL,KK)=PT(JL,KK)+FOELDCPM(PT(JL,KK))*ZCOND
        PQ(JL,KK)=PQ(JL,KK)-ZCOND
        ZQSAT=FOEEWM(PT(JL,KK))*ZQP    
        ZQSAT=MIN(0.5_dp,ZQSAT)
        ZCOR=1.0_dp/(1.0_dp-RETV  *ZQSAT)
        ZQSAT=ZQSAT*ZCOR
        ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.0_dp+ZQSAT*ZCOR*FOEDEM(PT(JL,KK)))
        PT(JL,KK)=PT(JL,KK)+FOELDCPM(PT(JL,KK))*ZCOND1
        PQ(JL,KK)=PQ(JL,KK)-ZCOND1
      ENDIF
    ENDDO
  ENDIF 

  IF(KCALL == 5) THEN  ! Same as 4 but with LDFLAG all true

!DIR$ IVDEP
!OCL NOVREC
    IF(N_VMASS <= 0)  THEN ! Not using Vector MASS
      DO JL=KIDIA,KFDIA
        ZQP    =1.0_dp/PSP(JL)
        ZQSAT=FOEEWM(PT(JL,KK))*ZQP    
        ZQSAT=MIN(0.5_dp,ZQSAT)
        ZCOR=1.0_dp/(1.0_dp-RETV  *ZQSAT)
        ZQSAT=ZQSAT*ZCOR
        ZCOND=(PQ(JL,KK)-ZQSAT)/(1.0_dp+ZQSAT*ZCOR*FOEDEM(PT(JL,KK)))
        PT(JL,KK)=PT(JL,KK)+FOELDCPM(PT(JL,KK))*ZCOND
        PQ(JL,KK)=PQ(JL,KK)-ZCOND
        ZQSAT=FOEEWM(PT(JL,KK))*ZQP    
        ZQSAT=MIN(0.5_dp,ZQSAT)
        ZCOR=1.0_dp/(1.0_dp-RETV  *ZQSAT)
        ZQSAT=ZQSAT*ZCOR
        ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.0_dp+ZQSAT*ZCOR*FOEDEM(PT(JL,KK)))
        PT(JL,KK)=PT(JL,KK)+FOELDCPM(PT(JL,KK))*ZCOND1
        PQ(JL,KK)=PQ(JL,KK)-ZCOND1
      ENDDO
    ELSE ! Using Vector VMASS
      DO JL=KIDIA,KFDIA
        ZTMP1(JL-KIDIA+1)=R3LES*(PT(JL,KK)-RTT)
        ZTMP2(JL-KIDIA+1)=R3IES*(PT(JL,KK)-RTT)
        ZTMP3(JL-KIDIA+1)=PT(JL,KK)-R4LES
        ZTMP4(JL-KIDIA+1)=PT(JL,KK)-R4IES
      ENDDO
      CALL VDIV(ZTMP5,ZTMP1,ZTMP3,JLEN)
      CALL VDIV(ZTMP6,ZTMP2,ZTMP4,JLEN)
      CALL VEXP(ZTMP1,ZTMP5,JLEN)
      CALL VEXP(ZTMP2,ZTMP6,JLEN)
      DO JL=KIDIA,KFDIA
        ZQP    =1.0_dp/PSP(JL)
        ZQSAT=R2ES*(FOEALFA(PT(JL,KK))*ZTMP1(JL-KIDIA+1)+&
         & (1.0_dp-FOEALFA(PT(JL,KK)))*ZTMP2(JL-KIDIA+1))*ZQP  
        ZQSAT=MINJ(0.5_dp,ZQSAT)
        ZCOR=1.0_dp-RETV*ZQSAT
        ZCOND=(PQ(JL,KK)*ZCOR**2-ZQSAT*ZCOR)/(ZCOR**2+ZQSAT*FOEDEM(PT(JL,KK)))
        PT(JL,KK)=PT(JL,KK)+FOELDCPM(PT(JL,KK))*ZCOND
        PQ(JL,KK)=PQ(JL,KK)-ZCOND
        ZTMP0(JL-KIDIA+1)=ZQP
        ZTMP1(JL-KIDIA+1)=R3LES*(PT(JL,KK)-RTT)
        ZTMP2(JL-KIDIA+1)=R3IES*(PT(JL,KK)-RTT)
        ZTMP3(JL-KIDIA+1)=PT(JL,KK)-R4LES
        ZTMP4(JL-KIDIA+1)=PT(JL,KK)-R4IES
      ENDDO
      CALL VDIV(ZTMP5,ZTMP1,ZTMP3,JLEN)
      CALL VDIV(ZTMP6,ZTMP2,ZTMP4,JLEN)
      CALL VEXP(ZTMP1,ZTMP5,JLEN)
      CALL VEXP(ZTMP2,ZTMP6,JLEN)
      DO JL=KIDIA,KFDIA
        ZQP  = ZTMP0(JL-KIDIA+1)
        ZQSAT=R2ES*(FOEALFA(PT(JL,KK))*ZTMP1(JL-KIDIA+1)+&
         & (1.0_dp-FOEALFA(PT(JL,KK)))*ZTMP2(JL-KIDIA+1))*ZQP  
        ZQSAT=MINJ(0.5_dp,ZQSAT)
        ZCOR=1.0_dp-RETV*ZQSAT
        ZCOND1=(PQ(JL,KK)*ZCOR**2-ZQSAT*ZCOR)/(ZCOR**2+ZQSAT*FOEDEM(PT(JL,KK)))
        PT(JL,KK)=PT(JL,KK)+FOELDCPM(PT(JL,KK))*ZCOND1
        PQ(JL,KK)=PQ(JL,KK)-ZCOND1
      ENDDO
    ENDIF
  ENDIF 

  IF(KCALL == 3) THEN 
!DIR$ IVDEP 
!OCL NOVREC 
    IF(N_VMASS <=  0)  THEN ! Not using Vector MASS
      DO JL=KIDIA,KFDIA
        ZQP    =1.0_dp/PSP(JL)
        ZQSAT=FOEEWMCU(PT(JL,KK))*ZQP
        ZQSAT=MIN(0.5_dp,ZQSAT)
        ZCOR=1.0_dp/(1.0_dp-RETV  *ZQSAT)
        ZQSAT=ZQSAT*ZCOR
        ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.0_dp+ZQSAT*ZCOR*FOEDEMCU(PT(JL,KK)))
        PT(JL,KK)=PT(JL,KK)+FOELDCPMCU(PT(JL,KK))*ZCOND1
        PQ(JL,KK)=PQ(JL,KK)-ZCOND1
        ZQSAT=FOEEWMCU(PT(JL,KK))*ZQP
        ZQSAT=MIN(0.5_dp,ZQSAT)
        ZCOR=1.0_dp/(1.0_dp-RETV  *ZQSAT)
        ZQSAT=ZQSAT*ZCOR
        ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.0_dp+ZQSAT*ZCOR*FOEDEMCU(PT(JL,KK)))
        PT(JL,KK)=PT(JL,KK)+FOELDCPMCU(PT(JL,KK))*ZCOND1
        PQ(JL,KK)=PQ(JL,KK)-ZCOND1
      ENDDO
    ELSE
      DO JL=KIDIA,KFDIA 
        ZTMP1(JL-KIDIA+1)=R3LES*(PT(JL,KK)-RTT) 
        ZTMP2(JL-KIDIA+1)=R3IES*(PT(JL,KK)-RTT) 
        ZTMP3(JL-KIDIA+1)=PT(JL,KK)-R4LES 
        ZTMP4(JL-KIDIA+1)=PT(JL,KK)-R4IES 
      ENDDO 
      CALL VDIV(ZTMP5,ZTMP1,ZTMP3,JLEN)
      CALL VDIV(ZTMP6,ZTMP2,ZTMP4,JLEN)
      CALL VEXP(ZTMP1,ZTMP5,JLEN)
      CALL VEXP(ZTMP2,ZTMP6,JLEN)
      DO JL=KIDIA,KFDIA
        ZQP    =1.0_dp/PSP(JL)
        ZQSAT=R2ES*(FOEALFCU(PT(JL,KK))*ZTMP1(JL-KIDIA+1)+&
         & (1.0_dp-FOEALFCU(PT(JL,KK)))*ZTMP2(JL-KIDIA+1))*ZQP  
        ZQSAT=MINJ(0.5_dp,ZQSAT)
        ZCOR=1.0_dp-RETV*ZQSAT
        ZCOND1=(PQ(JL,KK)*ZCOR**2-ZQSAT*ZCOR)/(ZCOR**2+ZQSAT*FOEDEMCU(PT(JL,KK)))
        PT(JL,KK)=PT(JL,KK)+FOELDCPMCU(PT(JL,KK))*ZCOND1
        PQ(JL,KK)=PQ(JL,KK)-ZCOND1
        ZTMP0(JL-KIDIA+1)=ZQP
        ZTMP1(JL-KIDIA+1)=R3LES*(PT(JL,KK)-RTT)
        ZTMP2(JL-KIDIA+1)=R3IES*(PT(JL,KK)-RTT)
        ZTMP3(JL-KIDIA+1)=PT(JL,KK)-R4LES
        ZTMP4(JL-KIDIA+1)=PT(JL,KK)-R4IES
      ENDDO
      CALL VDIV(ZTMP5,ZTMP1,ZTMP3,JLEN)
      CALL VDIV(ZTMP6,ZTMP2,ZTMP4,JLEN)
      CALL VEXP(ZTMP1,ZTMP5,JLEN)
      CALL VEXP(ZTMP2,ZTMP6,JLEN)
      DO JL=KIDIA,KFDIA
        ZQP  = ZTMP0(JL-KIDIA+1)
        ZQSAT=R2ES*(FOEALFCU(PT(JL,KK))*ZTMP1(JL-KIDIA+1)+&
         & (1.0_dp-FOEALFCU(PT(JL,KK)))*ZTMP2(JL-KIDIA+1))*ZQP  
        ZQSAT=MINJ(0.5_dp,ZQSAT)
        ZCOR=1.0_dp-RETV*ZQSAT
        ZCOND1=(PQ(JL,KK)*ZCOR**2-ZQSAT*ZCOR)/(ZCOR**2+ZQSAT*FOEDEMCU(PT(JL,KK)))
        PT(JL,KK)=PT(JL,KK)+FOELDCPMCU(PT(JL,KK))*ZCOND1
        PQ(JL,KK)=PQ(JL,KK)-ZCOND1
      ENDDO
    ENDIF

  ENDIF
!*********************************************
ELSE
!*********************************************                 

!     2.           CALCULATE CONDENSATION AND ADJUST T AND Q ACCORDINGLY
!                  -----------------------------------------------------

  IF (KCALL == 1 ) THEN

!DIR$ IVDEP
!OCL NOVREC
    DO JL=KIDIA,KFDIA
      IF(LDFLAG(JL)) THEN
        ZQP    =1.0_dp/PSP(JL)
        ZTARG=PT(JL,KK)
        ZOEALFA=0.5_dp*(TANH(RLPAL1*(ZTARG-RLPTRC))+1.0_dp)
        ZFOEEWL=R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
        ZFOEEWI=R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
        ZQSAT=ZQP    *(ZOEALFA*ZFOEEWL+(1.0_dp-ZOEALFA)*ZFOEEWI)
        Z1S=TANH(RLPAL2*(ZQSAT-ZQMAX))
        ZQSAT=0.5_dp*((1.0_dp-Z1S)*ZQSAT+(1.0_dp+Z1S)*ZQMAX)

        ZCOR=1.0_dp/(1.0_dp-RETV  *ZQSAT)
        ZQSAT=ZQSAT*ZCOR

        Z2S=    ZOEALFA *R5ALVCP*(1.0_dp/(ZTARG-R4LES)**2)+&
         & (1.0_dp-ZOEALFA)*R5ALSCP*(1.0_dp/(ZTARG-R4IES)**2)  
        ZCOND=(PQ(JL,KK)-ZQSAT)/(1.0_dp+ZQSAT*ZCOR*Z2S)

        ZCOND=MAX(ZCOND,0.0_dp)

        IF(ZCOND /= 0.0_dp) THEN

          PT(JL,KK)=PT(JL,KK)+&
           & (ZOEALFA*RALVDCP+(1.0_dp-ZOEALFA)*RALSDCP)*ZCOND  
          PQ(JL,KK)=PQ(JL,KK)-ZCOND
          ZTARG=PT(JL,KK)
          ZOEALFA=0.5_dp*(TANH(RLPAL1*(ZTARG-RLPTRC))+1.0_dp)
          ZFOEEWL=R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
          ZFOEEWI=R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
          ZQSAT=ZQP    *(ZOEALFA*ZFOEEWL+(1.0_dp-ZOEALFA)*ZFOEEWI)
          Z1S=TANH(RLPAL2*(ZQSAT-ZQMAX))
          ZQSAT=0.5_dp*((1.0_dp-Z1S)*ZQSAT+(1.0_dp+Z1S)*ZQMAX)
  
          ZCOR=1.0_dp/(1.0_dp-RETV  *ZQSAT)
          ZQSAT=ZQSAT*ZCOR
  
          Z2S=    ZOEALFA *R5ALVCP*(1.0_dp/(ZTARG-R4LES)**2)+&
           & (1.0_dp-ZOEALFA)*R5ALSCP*(1.0_dp/(ZTARG-R4IES)**2)  
          ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.0_dp+ZQSAT*ZCOR*Z2S)
  
          PT(JL,KK)=PT(JL,KK)+(ZOEALFA*RALVDCP+(1.0_dp-ZOEALFA)*RALSDCP)*ZCOND1
  
          PQ(JL,KK)=PQ(JL,KK)-ZCOND1
        ENDIF
      ENDIF
    ENDDO

  ENDIF

  IF(KCALL == 2) THEN

!DIR$ IVDEP
!OCL NOVREC
    DO JL=KIDIA,KFDIA
      IF(LDFLAG(JL)) THEN
        ZQP    =1.0_dp/PSP(JL)

        ZTARG=PT(JL,KK)
        ZOEALFA=0.5_dp*(TANH(RLPAL1*(ZTARG-RLPTRC))+1.0_dp)
        ZFOEEWL=R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
        ZFOEEWI=R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
        ZQSAT=ZQP    *(ZOEALFA*ZFOEEWL+(1.0_dp-ZOEALFA)*ZFOEEWI)
        Z1S=TANH(RLPAL2*(ZQSAT-ZQMAX))
        ZQSAT=0.5_dp*((1.0_dp-Z1S)*ZQSAT+(1.0_dp+Z1S)*ZQMAX)

        ZCOR=1.0_dp/(1.0_dp-RETV  *ZQSAT)
        ZQSAT=ZQSAT*ZCOR

        Z2S=    ZOEALFA *R5ALVCP*(1.0_dp/(ZTARG-R4LES)**2)+&
         & (1.0_dp-ZOEALFA)*R5ALSCP*(1.0_dp/(ZTARG-R4IES)**2)  
        ZCOND=(PQ(JL,KK)-ZQSAT)/(1.0_dp+ZQSAT*ZCOR*Z2S)

        ZCOND=MIN(ZCOND,0.0_dp)

        IF(ZCOND /= 0.0_dp) THEN

          PT(JL,KK)=PT(JL,KK)+&
           & (ZOEALFA*RALVDCP+(1.0_dp-ZOEALFA)*RALSDCP)*ZCOND  
          PQ(JL,KK)=PQ(JL,KK)-ZCOND
          ZTARG=PT(JL,KK)
          ZOEALFA=0.5_dp*(TANH(RLPAL1*(ZTARG-RLPTRC))+1.0_dp)
          ZFOEEWL=R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
          ZFOEEWI=R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
          ZQSAT=ZQP    *(ZOEALFA*ZFOEEWL+(1.0_dp-ZOEALFA)*ZFOEEWI)
          Z1S=TANH(RLPAL2*(ZQSAT-ZQMAX))
          ZQSAT=0.5_dp*((1.0_dp-Z1S)*ZQSAT+(1.0_dp+Z1S)*ZQMAX)
  
          ZCOR=1.0_dp/(1.0_dp-RETV  *ZQSAT)
          ZQSAT=ZQSAT*ZCOR
  
          Z2S=    ZOEALFA *R5ALVCP*(1.0_dp/(ZTARG-R4LES)**2)+&
           & (1.0_dp-ZOEALFA)*R5ALSCP*(1.0_dp/(ZTARG-R4IES)**2)  
          ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.0_dp+ZQSAT*ZCOR*Z2S)
  
          PT(JL,KK)=PT(JL,KK)+(ZOEALFA*RALVDCP+(1.0_dp-ZOEALFA)*RALSDCP)*ZCOND1
  
          PQ(JL,KK)=PQ(JL,KK)-ZCOND1
        ENDIF
      ENDIF
    ENDDO

  ENDIF

  IF(KCALL == 0) THEN

!DIR$ IVDEP
!OCL NOVREC
    DO JL=KIDIA,KFDIA
      ZQP    =1.0_dp/PSP(JL)

      ZTARG=PT(JL,KK)
      ZOEALFA=0.5_dp*(TANH(RLPAL1*(ZTARG-RLPTRC))+1.0_dp)
      ZFOEEWL=R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
      ZFOEEWI=R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
      ZQSAT=ZQP    *(ZOEALFA*ZFOEEWL+(1.0_dp-ZOEALFA)*ZFOEEWI)
      Z1S=TANH(RLPAL2*(ZQSAT-ZQMAX))
      ZQSAT=0.5_dp*((1.0_dp-Z1S)*ZQSAT+(1.0_dp+Z1S)*ZQMAX)

      ZCOR=1.0_dp/(1.0_dp-RETV  *ZQSAT)
      ZQSAT=ZQSAT*ZCOR

      Z2S=    ZOEALFA *R5ALVCP*(1.0_dp/(ZTARG-R4LES)**2)+&
       & (1.0_dp-ZOEALFA)*R5ALSCP*(1.0_dp/(ZTARG-R4IES)**2)  
      ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.0_dp+ZQSAT*ZCOR*Z2S)

      PT(JL,KK)=PT(JL,KK)+(ZOEALFA*RALVDCP+(1.0_dp-ZOEALFA)*RALSDCP)*ZCOND1

      PQ(JL,KK)=PQ(JL,KK)-ZCOND1

      ZTARG=PT(JL,KK)
      ZOEALFA=0.5_dp*(TANH(RLPAL1*(ZTARG-RLPTRC))+1.0_dp)
      ZFOEEWL=R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
      ZFOEEWI=R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
      ZQSAT=ZQP    *(ZOEALFA*ZFOEEWL+(1.0_dp-ZOEALFA)*ZFOEEWI)
      Z1S=TANH(RLPAL2*(ZQSAT-ZQMAX))
      ZQSAT=0.5_dp*((1.0_dp-Z1S)*ZQSAT+(1.0_dp+Z1S)*ZQMAX)

      ZCOR=1.0_dp/(1.0_dp-RETV  *ZQSAT)
      ZQSAT=ZQSAT*ZCOR

      Z2S=    ZOEALFA *R5ALVCP*(1.0_dp/(ZTARG-R4LES)**2)+&
       & (1.0_dp-ZOEALFA)*R5ALSCP*(1.0_dp/(ZTARG-R4IES)**2)  
      ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.0_dp+ZQSAT*ZCOR*Z2S)

      PT(JL,KK)=PT(JL,KK)+(ZOEALFA*RALVDCP+(1.0_dp-ZOEALFA)*RALSDCP)*ZCOND1

      PQ(JL,KK)=PQ(JL,KK)-ZCOND1
    ENDDO

  ENDIF

  IF(KCALL == 4) THEN

!DIR$ IVDEP
!OCL NOVREC
    DO JL=KIDIA,KFDIA
      IF(LDFLAG(JL)) THEN
        ZQP    =1.0_dp/PSP(JL)

        ZTARG=PT(JL,KK)
        ZOEALFA=0.5_dp*(TANH(RLPAL1*(ZTARG-RLPTRC))+1.0_dp)
        ZFOEEWL=R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
        ZFOEEWI=R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
        ZQSAT=ZQP    *(ZOEALFA*ZFOEEWL+(1.0_dp-ZOEALFA)*ZFOEEWI)
        Z1S=TANH(RLPAL2*(ZQSAT-ZQMAX))
        ZQSAT=0.5_dp*((1.0_dp-Z1S)*ZQSAT+(1.0_dp+Z1S)*ZQMAX)
        
        ZCOR=1.0_dp/(1.0_dp-RETV  *ZQSAT)
        ZQSAT=ZQSAT*ZCOR
        
        Z2S=    ZOEALFA *R5ALVCP*(1.0_dp/(ZTARG-R4LES)**2)+&
         & (1.0_dp-ZOEALFA)*R5ALSCP*(1.0_dp/(ZTARG-R4IES)**2)  
        ZCOND=(PQ(JL,KK)-ZQSAT)/(1.0_dp+ZQSAT*ZCOR*Z2S)
        
        PT(JL,KK)=PT(JL,KK)+(ZOEALFA*RALVDCP+(1.0_dp-ZOEALFA)*RALSDCP)*ZCOND
        
        PQ(JL,KK)=PQ(JL,KK)-ZCOND
        
        ZTARG=PT(JL,KK)
        ZOEALFA=0.5_dp*(TANH(RLPAL1*(ZTARG-RLPTRC))+1.0_dp)
        ZFOEEWL=R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
        ZFOEEWI=R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
        ZQSAT=ZQP    *(ZOEALFA*ZFOEEWL+(1.0_dp-ZOEALFA)*ZFOEEWI)
        Z1S=TANH(RLPAL2*(ZQSAT-ZQMAX))
        ZQSAT=0.5_dp*((1.0_dp-Z1S)*ZQSAT+(1.0_dp+Z1S)*ZQMAX)
        
        ZQSAT=MIN(ZQMAX,ZQSAT)
        ZCOR=1.0_dp/(1.0_dp-RETV  *ZQSAT)
        ZQSAT=ZQSAT*ZCOR
        
        Z2S=    ZOEALFA *R5ALVCP*(1.0_dp/(ZTARG-R4LES)**2)+&
         & (1.0_dp-ZOEALFA)*R5ALSCP*(1.0_dp/(ZTARG-R4IES)**2)  
        ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.0_dp+ZQSAT*ZCOR*Z2S)
        
        PT(JL,KK)=PT(JL,KK)+(ZOEALFA*RALVDCP+(1.0_dp-ZOEALFA)*RALSDCP)*ZCOND1

        PQ(JL,KK)=PQ(JL,KK)-ZCOND1
      ENDIF
    ENDDO

  ENDIF

!*********************************************
ENDIF
!*********************************************                 

IF (LHOOK) CALL DR_HOOK('CUADJTQ',1,ZHOOK_HANDLE)
END SUBROUTINE CUADJTQ

!==============================================================================

SUBROUTINE CUADJTQS &
 & (KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
 & KK,&
 & PSP,      PT,       PQ,       LDFLAG,   KCALL)  

!**   *CUADJTQS* - SIMPLIFIED VERSION OF MOIST ADJUSTMENT

!     J.F. MAHFOUF      ECMWF         

!     PURPOSE.
!     --------
!     TO PRODUCE T,Q AND L VALUES FOR CLOUD ASCENT

!     INTERFACE
!     ---------
!     THIS ROUTINE IS CALLED FROM SUBROUTINES:

!       *COND*       
!       *CUBMADJ*    
!       *CUBMD*      
!       *CONDAD*     
!       *CUBMADJAD*  
!       *CUBMDAD*    

!     INPUT ARE UNADJUSTED T AND Q VALUES,
!     IT RETURNS ADJUSTED VALUES OF T AND Q

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KTDIA*        START OF THE VERTICAL LOOP
!    *KLEV*         NUMBER OF LEVELS
!    *KK*           LEVEL
!    *KCALL*        DEFINES CALCULATION AS
!                      KCALL=0  ENV. T AND QS IN*CUINI*
!                      KCALL=1  CONDENSATION IN UPDRAFTS  (E.G. CUBASE, CUASC)
!                      KCALL=2  EVAPORATION IN DOWNDRAFTS (E.G. CUDLFS,CUDDRAF)

!     INPUT PARAMETERS (LOGICAL):

!    *LDLAND*       LAND-SEA MASK (.TRUE. FOR LAND POINTS)

!     INPUT PARAMETERS (REAL):

!    *PSP*          PRESSURE                                        PA

!     UPDATED PARAMETERS (REAL):

!    *PT*           TEMPERATURE                                     K
!    *PQ*           SPECIFIC HUMIDITY                             KG/KG

!          MODIFICATIONS
!          -------------
!          D.SALMOND & M.HAMRUD ECMWF       99-06-04   Optimisation
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!----------------------------------------------------------------------

USE MESSY_CONVECT_ECMWF_PARAM   , ONLY : RETV, RLVTT, RLSTT, RTT,                 &
                                         R2ES, R3LES, R3IES, R4LES,               &
                                         R4IES, R5LES, R5IES, R5ALVCP, R5ALSCP,   &
                                         RALVDCP, RALSDCP, RTWAT, RTICE, RTICECU, &
                                         RTWAT_RTICE_R, RTWAT_RTICECU_R,          &  
                                         LHOOK, DR_HOOK                            !from YOMHOOK

IMPLICIT NONE

INTEGER,INTENT(IN)    :: KLON 
INTEGER,INTENT(IN)    :: KLEV 
INTEGER,INTENT(IN)    :: KIDIA 
INTEGER,INTENT(IN)    :: KFDIA 
INTEGER               :: KTDIA ! Argument NOT used
INTEGER,INTENT(IN)    :: KK 
REAL(dp)   ,INTENT(IN)    :: PSP(KLON) 
REAL(dp)   ,INTENT(INOUT) :: PT(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PQ(KLON,KLEV) 
LOGICAL           ,INTENT(IN)    :: LDFLAG(KLON) 
INTEGER,INTENT(IN)    :: KCALL 
REAL(dp) ::     Z3ES(KLON),             Z4ES(KLON),&
 & Z5ALCP(KLON),           ZALDCP(KLON)  

INTEGER :: JL

REAL(dp) :: ZQMAX, ZQP, ZCOND, ZCOND1, ZTARG, ZCOR, ZQSAT, ZFOEEW, Z2S
REAL(dp) :: ZHOOK_HANDLE

INTRINSIC :: EXP, MAX, MIN
!----------------------------------------------------------------------

!     1.           DEFINE CONSTANTS
!                  ----------------

IF (LHOOK) CALL DR_HOOK('CUADJTQS',0,ZHOOK_HANDLE)
ZQMAX=0.5_dp

!     2.           CALCULATE CONDENSATION AND ADJUST T AND Q ACCORDINGLY
!                  -----------------------------------------------------

!*    ICE-WATER THERMODYNAMICAL FUNCTIONS

DO JL=KIDIA,KFDIA
  IF (PT(JL,KK) > RTT) THEN
    Z3ES(JL)=R3LES
    Z4ES(JL)=R4LES
    Z5ALCP(JL)=R5ALVCP
    ZALDCP(JL)=RALVDCP
  ELSE
    Z3ES(JL)=R3IES
    Z4ES(JL)=R4IES
    Z5ALCP(JL)=R5ALSCP
    ZALDCP(JL)=RALSDCP
  ENDIF
ENDDO

IF (KCALL == 1 ) THEN

!DIR$ IVDEP
!OCL NOVREC
  DO JL=KIDIA,KFDIA
    IF(LDFLAG(JL)) THEN
      ZQP    =1.0_dp/PSP(JL)
      ZTARG    =PT(JL,KK)
      ZFOEEW    =R2ES*EXP(Z3ES(JL)*(ZTARG    -RTT)/(ZTARG    -Z4ES(JL)))
      ZQSAT    =ZQP    *ZFOEEW    
      IF (ZQSAT     > ZQMAX) THEN
        ZQSAT    =ZQMAX
      ENDIF
      ZCOR    =1.0_dp/(1.0_dp-RETV*ZQSAT    )
      ZQSAT    =ZQSAT    *ZCOR    
      Z2S    =Z5ALCP(JL)/(ZTARG    -Z4ES(JL))**2
      ZCOND    =(PQ(JL,KK)-ZQSAT    )/(1.0_dp+ZQSAT    *ZCOR    *Z2S    )
      ZCOND    =MAX(ZCOND    ,0.0_dp)
!     IF(ZCOND /= 0._dp) THEN
      PT(JL,KK)=PT(JL,KK)+ZALDCP(JL)*ZCOND    
      PQ(JL,KK)=PQ(JL,KK)-ZCOND    
      ZTARG    =PT(JL,KK)
      ZFOEEW    =R2ES*EXP(Z3ES(JL)*(ZTARG    -RTT)/(ZTARG    -Z4ES(JL)))
      ZQSAT    =ZQP    *ZFOEEW    
      IF (ZQSAT     > ZQMAX) THEN
        ZQSAT    =ZQMAX
      ENDIF
      ZCOR    =1.0_dp/(1.0_dp-RETV*ZQSAT    )
      ZQSAT    =ZQSAT    *ZCOR    
      Z2S    =Z5ALCP(JL)/(ZTARG    -Z4ES(JL))**2
      ZCOND1    =(PQ(JL,KK)-ZQSAT    )/(1.0_dp+ZQSAT    *ZCOR    *Z2S    )
      IF(ZCOND ==  0.0_dp)ZCOND1=0.0_dp
      PT(JL,KK)=PT(JL,KK)+ZALDCP(JL)*ZCOND1    
      PQ(JL,KK)=PQ(JL,KK)-ZCOND1    
!     ENDIF
    ENDIF
  ENDDO

ENDIF

IF(KCALL == 2) THEN

!DIR$ IVDEP
!OCL NOVREC
  DO JL=KIDIA,KFDIA
    IF(LDFLAG(JL)) THEN
      ZQP    =1.0_dp/PSP(JL)
      ZTARG    =PT(JL,KK)
      ZFOEEW    =R2ES*EXP(Z3ES(JL)*(ZTARG    -RTT)/(ZTARG    -Z4ES(JL)))
      ZQSAT    =ZQP    *ZFOEEW    
      IF (ZQSAT     > ZQMAX) THEN
        ZQSAT    =ZQMAX
      ENDIF
      ZCOR    =1.0_dp/(1.0_dp-RETV  *ZQSAT    )
      ZQSAT    =ZQSAT    *ZCOR    
      Z2S    =Z5ALCP(JL)/(ZTARG    -Z4ES(JL))**2
      ZCOND    =(PQ(JL,KK)-ZQSAT    )/(1.0_dp+ZQSAT    *ZCOR    *Z2S    )
      ZCOND    =MIN(ZCOND    ,0.0_dp)
!     IF(ZCOND /= 0._dp) THEN
      PT(JL,KK)=PT(JL,KK)+ZALDCP(JL)*ZCOND    
      PQ(JL,KK)=PQ(JL,KK)-ZCOND    
      ZTARG    =PT(JL,KK)
      ZFOEEW    =R2ES*EXP(Z3ES(JL)*(ZTARG    -RTT)/(ZTARG    -Z4ES(JL)))
      ZQSAT    =ZQP    *ZFOEEW    
      IF (ZQSAT     > ZQMAX) THEN
        ZQSAT    =ZQMAX
      ENDIF
      ZCOR    =1.0_dp/(1.0_dp-RETV  *ZQSAT    )
      ZQSAT    =ZQSAT    *ZCOR    
      Z2S    =Z5ALCP(JL)/(ZTARG    -Z4ES(JL))**2
      ZCOND1    =(PQ(JL,KK)-ZQSAT    )/(1.0_dp+ZQSAT    *ZCOR    *Z2S    )
      IF(ZCOND ==  0.0_dp)ZCOND1=0.0_dp
      PT(JL,KK)=PT(JL,KK)+ZALDCP(JL)*ZCOND1    
      PQ(JL,KK)=PQ(JL,KK)-ZCOND1    
!     ENDIF
    ENDIF
  ENDDO

ENDIF

IF(KCALL == 0) THEN

!DIR$ IVDEP
!OCL NOVREC
  DO JL=KIDIA,KFDIA
    ZQP    =1.0_dp/PSP(JL)
    ZTARG    =PT(JL,KK)
    ZFOEEW    =R2ES*EXP(Z3ES(JL)*(ZTARG    -RTT)/(ZTARG    -Z4ES(JL)))
    ZQSAT    =ZQP    *ZFOEEW    
    IF (ZQSAT     > ZQMAX) THEN
      ZQSAT    =ZQMAX
    ENDIF
    ZCOR    =1.0_dp/(1.0_dp-RETV  *ZQSAT    )
    ZQSAT    =ZQSAT    *ZCOR    
    Z2S    =Z5ALCP(JL)/(ZTARG    -Z4ES(JL))**2
    ZCOND1    =(PQ(JL,KK)-ZQSAT    )/(1.0_dp+ZQSAT    *ZCOR    *Z2S    )
    PT(JL,KK)=PT(JL,KK)+ZALDCP(JL)*ZCOND1    
    PQ(JL,KK)=PQ(JL,KK)-ZCOND1    
    ZTARG    =PT(JL,KK)
    ZFOEEW    =R2ES*EXP(Z3ES(JL)*(ZTARG    -RTT)/(ZTARG    -Z4ES(JL)))
    ZQSAT    =ZQP    *ZFOEEW    
    IF (ZQSAT     > ZQMAX) THEN
      ZQSAT    =ZQMAX
    ENDIF
    ZCOR    =1.0_dp/(1.0_dp-RETV  *ZQSAT    )
    ZQSAT    =ZQSAT    *ZCOR    
    Z2S    =Z5ALCP(JL)/(ZTARG    -Z4ES(JL))**2
    ZCOND1    =(PQ(JL,KK)-ZQSAT    )/(1.0_dp+ZQSAT    *ZCOR    *Z2S    )
    PT(JL,KK)=PT(JL,KK)+ZALDCP(JL)*ZCOND1    
    PQ(JL,KK)=PQ(JL,KK)-ZCOND1    
  ENDDO

ENDIF

IF(KCALL == 4) THEN

!DIR$ IVDEP
!OCL NOVREC
  DO JL=KIDIA,KFDIA
    ZQP    =1.0_dp/PSP(JL)
    ZTARG    =PT(JL,KK)
    ZFOEEW    =R2ES*EXP(Z3ES(JL)*(ZTARG    -RTT)/(ZTARG    -Z4ES(JL)))
    ZQSAT    =ZQP    *ZFOEEW    
    IF (ZQSAT     > ZQMAX) THEN
      ZQSAT    =ZQMAX
    ENDIF
    ZCOR    =1.0_dp/(1.0_dp-RETV  *ZQSAT    )
    ZQSAT    =ZQSAT    *ZCOR    
    Z2S    =Z5ALCP(JL)/(ZTARG    -Z4ES(JL))**2
    ZCOND    =(PQ(JL,KK)-ZQSAT    )/(1.0_dp+ZQSAT    *ZCOR    *Z2S    )
    PT(JL,KK)=PT(JL,KK)+ZALDCP(JL)*ZCOND    
    PQ(JL,KK)=PQ(JL,KK)-ZCOND    
    ZTARG    =PT(JL,KK)
    ZFOEEW    =R2ES*EXP(Z3ES(JL)*(ZTARG    -RTT)/(ZTARG    -Z4ES(JL)))
    ZQSAT    =ZQP    *ZFOEEW    
    IF (ZQSAT     > ZQMAX) THEN
      ZQSAT    =ZQMAX
    ENDIF
    ZCOR    =1.0_dp/(1.0_dp-RETV  *ZQSAT    )
    ZQSAT    =ZQSAT    *ZCOR    
    Z2S    =Z5ALCP(JL)/(ZTARG    -Z4ES(JL))**2
    ZCOND1    =(PQ(JL,KK)-ZQSAT    )/(1.0_dp+ZQSAT    *ZCOR    *Z2S    )
    PT(JL,KK)=PT(JL,KK)+ZALDCP(JL)*ZCOND1    
    PQ(JL,KK)=PQ(JL,KK)-ZCOND1    
  ENDDO

ENDIF

IF (LHOOK) CALL DR_HOOK('CUADJTQS',1,ZHOOK_HANDLE)
END SUBROUTINE CUADJTQS

!==============================================================================

!OPTIONS XOPT(HSFUN)
SUBROUTINE CUASCN &
 & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
 & PTSPHY,&
 & PTENH,    PQENH,    PUEN,     PVEN,&
 & PTEN,     PQEN,     PQSEN,    PLITOT,&
 & PGEO,     PGEOH,    PAP,      PAPH,&
 & PTENQ,    PVERVEL,  pwubase,  KLWMIN,&
 & LDLAND,   LDCUM,    KTYPE,    KLAB,&
 & PTU,      PQU,      PLU,      PUU,      PVU,&
 & PMFU,     PMFUB,    PENTR,    PLGLAC,&
 & PMFUS,    PMFUQ,    PMFUL,    PLUDE,    PDMFUP,&
 & PDMFEN,&
 & KCBOT,    KCTOP,    KCTOP0,   KDPL,  &
 & PMFUDE_RATE,        PMFUEN_RATE,     &
 & lwc,      iwc,      rform,    sform   )  

!          THIS ROUTINE DOES THE CALCULATIONS FOR CLOUD ASCENTS
!          FOR CUMULUS PARAMETERIZATION

!          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89

!          PURPOSE.
!          --------
!          TO PRODUCE CLOUD ASCENTS FOR CU-PARAMETRIZATION
!          (VERTICAL PROFILES OF T,Q,L,U AND V AND CORRESPONDING
!           FLUXES AS WELL AS PRECIPITATION RATES)

!          INTERFACE
!          ---------

!          THIS ROUTINE IS CALLED FROM *CUMASTR*.

!          METHOD.
!          --------
!          LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD BASE
!          AND THEN CALCULATE MOIST ASCENT FOR
!          ENTRAINING/DETRAINING PLUME.
!          ENTRAINMENT AND DETRAINMENT RATES DIFFER FOR
!          SHALLOW AND DEEP CUMULUS CONVECTION.
!          IN CASE THERE IS NO PENETRATIVE OR SHALLOW CONVECTION
!          CHECK FOR POSSIBILITY OF MID LEVEL CONVECTION
!          (CLOUD BASE VALUES CALCULATED IN *CUBASMC*)

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KTDIA*        START OF THE VERTICAL LOOP
!    *KLEV*         NUMBER OF LEVELS
!    *KLWMIN*       LEVEL OF MAXIMUM VERTICAL VELOCITY
!    *KTYPE*        TYPE OF CONVECTION
!                       1 = PENETRATIVE CONVECTION
!                       2 = SHALLOW CONVECTION
!                       3 = MIDLEVEL CONVECTION
!    *KCBOT*        CLOUD BASE LEVEL
!    *KDPL*         DEPARTURE LEVEL FOR CONVECTION

!    INPUT PARAMETERS (REAL):

!    *PTSPHY*       TIME STEP FOR THE PHYSICS                      S
!    *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS          K
!    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS     KG/KG
!    *PUEN*         PROVISIONAL ENVIRONMENT U-VELOCITY (T+1)      M/S
!    *PVEN*         PROVISIONAL ENVIRONMENT V-VELOCITY (T+1)      M/S
!    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)      K
!    *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1) KG/KG
!    *PQSEN*        ENVIRONMENT SPEC. SATURATION HUMIDITY (T+1)  KG/KG
!    *PGEO*         GEOPOTENTIAL                                 M2/S2
!    *PLITOT*       GRID MEAN LIQUID WATER+ICE CONTENT           KG/KG
!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                  M2/S2
!    *PAP*          PROVISIONAL PRESSURE ON FULL LEVELS           PA
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS           PA
!    *PTENQ*        MOISTURE TENDENCY                            KG/(KG S)
!    *PVERVEL*      VERTICAL VELOCITY                            PA/S

!    INPUT PARAMETERS (LOGICAL):

!    *LDLAND*       LAND SEA MASK (.TRUE. FOR LAND)
!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS 

!    UPDATED PARAMETERS (INTEGER):

!    *KLAB*         FLAG KLAB=1 FOR SUBCLOUD LEVELS
!                        KLAB=2 FOR CLOUD LEVELS

!    UPDATED PARAMETERS (REAL):

!    *PTU*          TEMPERATURE IN UPDRAFTS                        K
!    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                   KG/KG
!    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS             KG/KG
!    *PUU*          U-VELOCITY IN UPDRAFTS                        M/S
!    *PVU*          V-VELOCITY IN UPDRAFTS                        M/S

!    OUTPUT PARAMETERS (INTEGER):

!    *KCTOP*        CLOUD TOP LEVEL
!    *KCTOP0*       FIRST GUESS OF CLOUD TOP LEVEL 

!    OUTPUT PARAMETERS (REAL):

!    *PMFU*         MASSFLUX IN UPDRAFTS                         KG/(M2*S)
!    *PMFUB*        MASSFLUX IN UPDRAFTS AT CLOUD BASE           KG/(M2*S)
!    *PENTR*        FRACTIONAL MASS ENTRAINMENT RATE              1/M
!    *PMFUS*        FLUX OF DRY STATIC ENERGY IN UPDRAFTS         J/(M2*S)
!    *PMFUQ*        FLUX OF SPEC. HUMIDITY IN UPDRAFTS           KG/(M2*S)
!    *PMFUL*        FLUX OF LIQUID WATER IN UPDRAFTS             KG/(M2*S)
!    *PLUDE*        DETRAINED LIQUID WATER                       KG/(M2*S)
!    *PLGLAC*       FROZEN CLOUD WATER CONTENT                   KG/KG
!    *PDMFUP*       FLUX DIFFERENCE OF PRECIP. IN UPDRAFTS       KG/(M2*S)
!    *PMFUDE_RATE*  UPDRAFT DETRAINMENT RATE                     KG/(M2*S)
!    *PMFUEN_RATE*  UPDRAFT ENTRAINMENT RATE                     KG/(M2*S)

!          EXTERNALS
!          ---------
!          *CUADJTQ* ADJUST T AND Q DUE TO CONDENSATION IN ASCENT
!          *CUENTR*  CALCULATE ENTRAINMENT/DETRAINMENT RATES
!          *CUBASMC* CALCULATE CLOUD BASE VALUES FOR MIDLEVEL CONVECTION

!          REFERENCE
!          ---------
!          (TIEDTKE,1989)

!          MODIFICATIONS
!          -------------
!             92-09-21 : Update to Cy44      J.-J. MORCRETTE
!             99-06-14 : Optimisation        D.SALMOND
!             01-05-22 : Modified flux limiter M.CULLEN
!             02-08-14 : Allow for departure level =/ KLEV  P.BECHTOLD
!             03-08-28 : Clean-up detrainment rates         P.BECHTOLD
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!----------------------------------------------------------------------

USE MESSY_CONVECT_ECMWF_PARAM   , ONLY : RG, RCPD, RETV, RLVTT, RLSTT, RTT,       &
                                         LPHYLIN, RLPTRC, RLPAL1,                 &
                                         R2ES, R3LES, R3IES, R4LES,               &
                                         R4IES, R5LES, R5IES, R5ALVCP, R5ALSCP,   &
                                         RALVDCP, RALSDCP, RALFDCP, RTWAT, RTBER, &
                                         RTBERCU, RTICE, RTICECU,                 &
                                         RTWAT_RTICECU_R, RTWAT_RTICE_R,          &
                                         RMFCMIN, RPRCON, LMFDUDV, RMFCFL,        &  ! from YOECUMF
                                         RMFSHCFL,                                &  ! from YOECUMF
                                         RLMIN,                                   &  ! from YOECLDP
                                         FOEALFCU,                                &  ! functions
                                         LHOOK, DR_HOOK                              ! from yomhook

IMPLICIT NONE

INTEGER,INTENT(IN)    :: KLON 
INTEGER,INTENT(IN)    :: KLEV 
INTEGER,INTENT(IN)    :: KIDIA 
INTEGER,INTENT(IN)    :: KFDIA 
INTEGER               :: KTDIA ! Argument NOT used
REAL(dp)   ,INTENT(IN)    :: PTSPHY 
REAL(dp)   ,INTENT(INOUT) :: PTENH(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PQENH(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PUEN(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PVEN(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PTEN(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PQEN(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PQSEN(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PLITOT(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PGEO(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PGEOH(KLON,KLEV+1) 
REAL(dp)   ,INTENT(IN)    :: PAP(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) 
REAL(dp)   ,INTENT(IN)    :: PTENQ(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PVERVEL(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PWUBASE(KLON) 
INTEGER,INTENT(IN)    :: KLWMIN(KLON) 
LOGICAL           ,INTENT(IN)    :: LDLAND(KLON) 
LOGICAL           ,INTENT(INOUT) :: LDCUM(KLON) 
INTEGER,INTENT(INOUT) :: KTYPE(KLON) 
INTEGER,INTENT(INOUT) :: KLAB(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PTU(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PQU(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PLU(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PUU(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PVU(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PMFU(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PMFUB(KLON) 
REAL(dp)   ,INTENT(INOUT) :: PENTR(KLON) 
REAL(dp)   ,INTENT(OUT)   :: PLGLAC(KLON,KLEV) 
REAL(dp)   ,INTENT(OUT)   :: PMFUS(KLON,KLEV) 
REAL(dp)   ,INTENT(OUT)   :: PMFUQ(KLON,KLEV) 
REAL(dp)   ,INTENT(OUT)   :: PMFUL(KLON,KLEV) 
REAL(dp)   ,INTENT(OUT)   :: PLUDE(KLON,KLEV) 
REAL(dp)   ,INTENT(OUT)   :: PDMFUP(KLON,KLEV) 
REAL(dp)   ,INTENT(OUT)   :: PDMFEN(KLON,KLEV) 
INTEGER,INTENT(INOUT) :: KCBOT(KLON) 
INTEGER,INTENT(OUT)   :: KCTOP(KLON) 
INTEGER,INTENT(INOUT) :: KCTOP0(KLON)
INTEGER,INTENT(IN)    :: KDPL(KLON) 
REAL(dp)   ,INTENT(OUT)   :: PMFUDE_RATE(KLON,KLEV), PMFUEN_RATE(KLON,KLEV) 
! mz_ht_20070622+
REAL(dp)   ,INTENT(INOUT) :: lwc(KLON,KLEV)
REAL(dp)   ,INTENT(INOUT) :: iwc(KLON,KLEV)
REAL(dp)   ,INTENT(INOUT) :: rform(KLON,KLEV)
REAL(dp)   ,INTENT(INOUT) :: sform(KLON,KLEV)
! mz_ht_20070622-
REAL(dp) ::     ZDMFEN(KLON),           ZDMFDE(KLON),&
 & ZMFUU(KLON),            ZMFUV(KLON),&
 & ZPBASE(KLON),           ZQOLD(KLON),&
 & ZLRAIN(KLON,KLEV),      ZKINE(KLON,KLEV),&
 & ZBUO(KLON,KLEV),        ZLUOLD(KLON),&
 & ZPRECIP(KLON)  
REAL(dp) ::     ZDLAND(KLON)
REAL(dp) ::     ZPH(KLON)
LOGICAL ::  LLFLAG(KLON), LLFLAGUV(KLON), LLO1(KLON), LLO3

INTEGER :: ICALL, IK, IS, JK, JL, IKB

REAL(dp) :: Z_CLCRIT, Z_CLDMAX, Z_CPRC2, Z_CWDRAG, Z_CWIFRAC, ZALFAW,&
 & ZBC, ZBE, ZBUOC, ZC, ZCBF, ZCONS2, ZD, ZDFI, &
 & ZDKBUO, ZDKEN, ZDMFDU, ZDMFEU, ZDNOPRC, ZDPHI, &
 & ZDT, ZFAC, ZFACBUO, ZINT, ZKEDKE, ZLCRIT, &
 & ZLEEN, ZLNEW, ZMFMAX, ZMFTEST, ZMFULK, ZMFUN, &
 & ZMFUQK, ZMFUSK, ZOEALFA, ZOEALFAP, ZPRCDGW, &
 & ZPRCON, ZQEEN, ZQUDE, ZRNEW, ZROLD, ZSCDE, &
 & ZSEEN, ZTGLACE, ZVI, ZVV, ZVW, ZWU, ZZ, ZZCO  

REAL(dp) ::  ZCHANGE,ZXS,ZXE
REAL(dp) :: ZHOOK_HANDLE

INTRINSIC :: EXP, MAX, MIN, SQRT, TANH
!----------------------------------------------------------------------

!*    1.           SPECIFY PARAMETERS
!                  ------------------

IF (LHOOK) CALL DR_HOOK('CUASCN',0,ZHOOK_HANDLE)
ZCONS2=RMFCFL/(RG*PTSPHY)
ZTGLACE=RTT  -13._dp
ZFACBUO=0.5_dp/(1.0_dp+0.5_dp)
ZPRCDGW=RPRCON/RG
Z_CLCRIT=0.5E-3_dp
Z_CLDMAX=5.E-3_dp
Z_CWIFRAC=0.5_dp
Z_CPRC2=0.5_dp
Z_CWDRAG=(3._dp/8._dp)*0.506_dp/0.2_dp/RG

!----------------------------------------------------------------------

!     2.           SET DEFAULT VALUES
!                  ------------------

LLO3=.FALSE.
DO JL=KIDIA,KFDIA
  ZMFUU(JL)=0.0_dp
  ZMFUV(JL)=0.0_dp
  ZLUOLD(JL)=0.0_dp
  IF(.NOT.LDCUM(JL)) THEN
    KCBOT(JL)=0
    PMFUB(JL)=0.0_dp
    PQU(JL,KLEV)=0.0_dp
    KTYPE(JL)=0
  ENDIF
ENDDO

! initalize various quantities
! note that liquid water and kinetic energy at cloud base is 
! preserved from cubase

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    IF (JK /= KCBOT(JL)) then 
      PLU(JL,JK)=0.0_dp
      ZKINE(JL,JK)=0.0_dp
    ENDIF
    PMFU(JL,JK)=0.0_dp
    PMFUS(JL,JK)=0.0_dp
    PMFUQ(JL,JK)=0.0_dp
    PMFUL(JL,JK)=0.0_dp
    PLUDE(JL,JK)=0.0_dp
    PLGLAC(JL,JK)=0.0_dp
    PDMFUP(JL,JK)=0.0_dp
    ZLRAIN(JL,JK)=0.0_dp
    ZBUO(JL,JK)=0.0_dp
    IF(.NOT.LDCUM(JL).OR.KTYPE(JL) == 3) KLAB(JL,JK)=0
    IF(.NOT.LDCUM(JL).AND.PAPH(JL,JK) < 4.E4_dp) KCTOP0(JL)=JK
    PDMFEN(JL,JK)=0.0_dp
    PMFUDE_RATE(JL,JK)=0.0_dp 
    PMFUEN_RATE(JL,JK)=0.0_dp
  ENDDO
ENDDO
!DIR$ IVDEP
!OCL NOVREC
DO JL=KIDIA,KFDIA
  IF(LDLAND(JL)) THEN
    ZDLAND(JL)=3.0E4_dp
    IF(LDCUM(JL)) THEN
      ZDPHI=PGEOH(JL,KCTOP0(JL))-PGEOH(JL,KCBOT(JL))
      IF(PTU(JL,KCTOP0(JL)) >= ZTGLACE) ZDLAND(JL)=ZDPHI
      ZDLAND(JL)=MAX(3.0E4_dp,ZDLAND(JL))
      ZDLAND(JL)=MIN(5.0E4_dp,ZDLAND(JL))
    ENDIF
  ENDIF
  IF(KTYPE(JL) == 3) LDCUM(JL)=.FALSE.
ENDDO

!----------------------------------------------------------------------

!     3.0          INITIALIZE VALUES AT cloud base LEVEL
!                  -------------------------------------

DO JL=KIDIA,KFDIA
  KCTOP(JL)=kcbot(jl)
  IF(LDCUM(JL)) THEN
    IKB=KCBOT(JL)
    ZKINE(JL,ikb)=0.5*PWUBASE(JL)**2
    PMFU(JL,IKB)=PMFUB(JL)
    PMFUS(JL,IKB)=PMFUB(JL)*(RCPD*PTU(JL,IKB)+PGEOH(JL,IKB))
    PMFUQ(JL,IKB)=PMFUB(JL)*PQU(JL,IKB)
    PMFUL(JL,IKB)=PMFUB(JL)*PLU(JL,IKB)
    IF(LMFDUDV) THEN
  ! IKB=KLEV
  ! ikb=kcbot(jl)
      ikb=kdpl(jl)
      ZMFUU(JL)=PMFUB(JL)*PUU(JL,ikb) !*UPG PB
      ZMFUV(JL)=PMFUB(JL)*PVU(JL,ikb)
    ENDIF
  endif
ENDDO

!----------------------------------------------------------------------

!     4.           DO ASCENT: SUBCLOUD LAYER (KLAB=1) ,CLOUDS (KLAB=2)
!                  BY DOING FIRST DRY-ADIABATIC ASCENT AND THEN
!                  BY ADJUSTING T,Q AND L ACCORDINGLY IN *CUADJTQ*,
!                  THEN CHECK FOR BUOYANCY AND SET FLAGS ACCORDINGLY
!                  -------------------------------------------------

DO JK=KLEV-1,3,-1

!                  SPECIFY CLOUD BASE VALUES FOR MIDLEVEL CONVECTION
!                  IN *CUBASMC* IN CASE THERE IS NOT ALREADY CONVECTION
!                  ----------------------------------------------------

  IK=JK
  CALL CUBASMCN &
   & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
   & IK,&
   & PTEN,     PQEN,     PQSEN,    PUEN,     PVEN,&
   & PVERVEL,  PGEO,     PGEOH,    LDCUM,    KTYPE,    KLAB,&
   & KCBOT,    PMFU,     PMFUB,    PENTR,    ZLRAIN,&
   & PTU,      PQU,      PLU,      PUU,      PVU,&
   & PMFUS,    PMFUQ,    PMFUL,    PDMFUP,   ZMFUU,    ZMFUV)  

  IS=0
  DO JL=KIDIA,KFDIA
    ZPRECIP(JL)=0.0_dp
    LLO1(JL)=.FALSE.
    IS=IS+KLAB(JL,JK+1)
    IF(KLAB(JL,JK+1) == 0) KLAB(JL,JK)=0
    IF((LDCUM(JL).AND.KLAB(JL,JK+1) == 2).OR.&
       & (KTYPE(JL) == 3 .and. KLAB(JL,JK+1) == 1)) THEN  
      LLFLAG(JL)=.TRUE.
    ELSE
      LLFLAG(JL)=.FALSE.
    ENDIF
    IF(KLAB(JL,JK+1) > 0) THEN
      LLFLAGUV(JL)=.TRUE.
    ELSE
      LLFLAGUV(JL)=.FALSE.
    ENDIF
    ZPH(JL)=PAPH(JL,JK)
    IF(KTYPE(JL) == 3.AND.JK == KCBOT(JL)) THEN
      ZMFMAX=(PAPH(JL,JK)-PAPH(JL,JK-1))*ZCONS2
      IF(PMFUB(JL) > ZMFMAX) THEN
        ZFAC=ZMFMAX/PMFUB(JL)
        PMFU(JL,JK+1)=PMFU(JL,JK+1)*ZFAC
        PMFUS(JL,JK+1)=PMFUS(JL,JK+1)*ZFAC
        PMFUQ(JL,JK+1)=PMFUQ(JL,JK+1)*ZFAC
        ZMFUU(JL)=ZMFUU(JL)*ZFAC
        ZMFUV(JL)=ZMFUV(JL)*ZFAC
        PMFUB(JL)=ZMFMAX
      ENDIF
    ENDIF
  ENDDO

  IF(IS > 0) LLO3=.TRUE.

!*                  SPECIFY ENTRAINMENT RATES IN *CUENTR*
!                   -------------------------------------

  IK=JK
  CALL CUENTR &
   & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
   & IK,       KLWMIN,   KTYPE,    KCBOT,    KCTOP0,&
   & LDCUM,    LLO3,&
   & PTENH,    PQENH,    PTENQ,    PAPH,     PAP,&
   & PMFU,     PENTR,&
   & ZPBASE,   ZDMFEN,   ZDMFDE )  

!                  DO ADIABATIC ASCENT FOR ENTRAINING/DETRAINING PLUME
!                  ---------------------------------------------------

  IF(LLO3) THEN

    DO JL=KIDIA,KFDIA
      IF(LLFLAG(JL)) THEN
        ZDMFDE(JL)=MIN(ZDMFDE(JL),0.75_dp*PMFU(JL,JK+1))
        IF(JK < KCBOT(JL)) THEN
          ZMFMAX=(PAPH(JL,JK)-PAPH(JL,JK-1))*ZCONS2
          if(ktype(jl)==2.and.ptsphy>1800._dp) ZMFMAX=ZMFMAX*RMFSHCFL
          ZXS=MAX(PMFU(JL,JK+1)-ZMFMAX,0.0_dp)
          ZMFTEST=PMFU(JL,JK+1)+ZDMFEN(JL)-ZDMFDE(JL)
          ZCHANGE=MAX(ZMFTEST-ZMFMAX,0.0_dp)
          ZXE=MAX(ZCHANGE-ZXS,0.0_dp)
          ZDMFEN(JL)=ZDMFEN(JL)-ZXE
          ZCHANGE=ZCHANGE-ZXE
          ZDMFDE(JL)=ZDMFDE(JL)+ZCHANGE
        ENDIF

        PDMFEN(JL,JK) = ZDMFEN(JL)-ZDMFDE(JL)

        PMFU(JL,JK)=PMFU(JL,JK+1)+ZDMFEN(JL)-ZDMFDE(JL)
        ZQEEN=PQENH(JL,JK+1)*ZDMFEN(JL)
        ZSEEN=(RCPD*PTENH(JL,JK+1)+PGEOH(JL,JK+1))*ZDMFEN(JL)
        IF(PLITOT(JL,JK) > RLMIN) THEN
          ZLEEN=PLITOT(JL,JK)*ZDMFEN(JL)
        ELSE
          ZLEEN=0.0_dp
        ENDIF
        ZSCDE=(RCPD*PTU(JL,JK+1)+PGEOH(JL,JK+1))*ZDMFDE(JL)
        ZQUDE=PQU(JL,JK+1)*ZDMFDE(JL)
        PLUDE(JL,JK)=PLU(JL,JK+1)*ZDMFDE(JL)
        ZMFUSK=PMFUS(JL,JK+1)+ZSEEN-ZSCDE
        ZMFUQK=PMFUQ(JL,JK+1)+ZQEEN-ZQUDE
        ZMFULK=PMFUL(JL,JK+1)+ZLEEN-PLUDE(JL,JK)
        PLU(JL,JK)=ZMFULK*(1.0_dp/MAX(RMFCMIN,PMFU(JL,JK)))
        PQU(JL,JK)=ZMFUQK*(1.0_dp/MAX(RMFCMIN,PMFU(JL,JK)))
        PTU(JL,JK)=(ZMFUSK*(1.0_dp/MAX(RMFCMIN,PMFU(JL,JK)))-&
         & PGEOH(JL,JK))/RCPD  
        PTU(JL,JK)=MAX(100._dp,PTU(JL,JK))
        PTU(JL,JK)=MIN(400._dp,PTU(JL,JK))
        ZQOLD(JL)=PQU(JL,JK)
        ZLRAIN(JL,JK)=ZLRAIN(JL,JK+1)*(PMFU(JL,JK+1)-ZDMFDE(JL))*&
         & (1.0_dp/MAX(RMFCMIN,PMFU(JL,JK)))  
        ZLUOLD(JL)=PLU(JL,JK)
      ELSE
        ZQOLD(JL)=0.0_dp
      ENDIF
        ! reset to environmental values if below departure level
      IF ( JK > KDPL(JL) ) THEN
        PTU(JL,JK)=PTENH(JL,JK)      
        PQU(JL,JK)=PQENH(JL,JK)      
        PLU(JL,JK)=0.0_dp
        ZLUOLD(JL)=PLU(JL,JK)
      ENDIF
    ENDDO

!                  DO CORRECTIONS FOR MOIST ASCENT
!                  BY ADJUSTING T,Q AND L IN *CUADJTQ*
!                  -----------------------------------

    IK=JK
    ICALL=1
    CALL CUADJTQ &
     & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
     & IK,&
     & ZPH,      PTU,      PQU,      LLFLAG,  ICALL )  

    IF (LPHYLIN) THEN

!DIR$ IVDEP
!OCL NOVREC
      DO JL=KIDIA,KFDIA
        IF(LLFLAG(JL).AND.PQU(JL,JK) /= ZQOLD(JL)) THEN
          ZOEALFA   = 0.5_dp*(TANH(RLPAL1*(PTU(JL,JK  )-RLPTRC))+1.0_dp)
          ZOEALFAP  = 0.5_dp*(TANH(RLPAL1*(PTU(JL,JK+1)-RLPTRC))+1.0_dp)
          PLGLAC(JL,JK)=PLU(JL,JK)*((1.0_dp-ZOEALFA)-(1.0_dp-ZOEALFAP))
          PTU(JL,JK)=PTU(JL,JK)+RALFDCP*PLGLAC(JL,JK)
        ENDIF
      ENDDO

    ELSE

!DIR$ IVDEP
!OCL NOVREC
      DO JL=KIDIA,KFDIA
        IF(LLFLAG(JL).AND.PQU(JL,JK) /= ZQOLD(JL)) THEN
          PLGLAC(JL,JK)=PLU(JL,JK)*((1.0_dp-FOEALFCU(PTU(JL,JK)))-&
           & (1.0_dp-FOEALFCU(PTU(JL,JK+1))))  
          PTU(JL,JK)=PTU(JL,JK)+RALFDCP*PLGLAC(JL,JK)
        ENDIF
      ENDDO

    ENDIF

    DO JL=KIDIA,KFDIA
      IF(LLFLAG(JL).AND.PQU(JL,JK) /= ZQOLD(JL)) THEN
        KLAB(JL,JK)=2
        PLU(JL,JK)=PLU(JL,JK)+ZQOLD(JL)-PQU(JL,JK)
        ! mz_ht_20070325+
        lwc(jl,jk) = plu(jl,jk) 
        iwc(jl,jk) = plglac(jl,jk) 
        ! mz_ht_20070325-
        ZBC=PTU(JL,JK)*(1.0_dp+RETV*PQU(JL,JK)-PLU(JL,JK+1)-ZLRAIN(JL,JK+1))
        ZBE=PTENH(JL,JK)*(1.0_dp+RETV*PQENH(JL,JK))
        ZBUO(JL,JK)=ZBC-ZBE

! set flags in case of midlevel convection

        IF(KTYPE(JL) == 3 .AND. KLAB(JL,JK+1)== 1) THEN
          IF(ZBUO(JL,JK) > -0.5_dp) THEN
            LDCUM(JL)=.TRUE.
            KCTOP(JL)=JK
            ZKINE(JL,JK)=0.5_dp
          ELSE
            KLAB(JL,JK)=0
            PMFU(JL,JK)=0.0_dp
            PLUDE(JL,JK)=0.0_dp
            PLU(JL,JK)=0.0_dp
          ENDIF
        ENDIF

        IF(KLAB(JL,JK+1) == 2) THEN

          IF(ZBUO(JL,JK) < 0.0_dp.AND.KLAB(JL,JK+1) == 2) THEN
            PTENH(JL,JK)=0.5_dp*(PTEN(JL,JK)+PTEN(JL,JK-1))
            PQENH(JL,JK)=0.5_dp*(PQEN(JL,JK)+PQEN(JL,JK-1))
            ZBUO(JL,JK)=ZBC-PTENH(JL,JK)*(1.0_dp+RETV*PQENH(JL,JK))
          ENDIF
          ZBUOC=(ZBUO(JL,JK)/(PTENH(JL,JK)*(1.0_dp+RETV*PQENH(JL,JK)))&
           & +ZBUO(JL,JK+1)/(PTENH(JL,JK+1)*(1.0_dp+RETV*&
           & PQENH(JL,JK+1))))*0.5_dp  
          ZDKBUO=(PGEOH(JL,JK)-PGEOH(JL,JK+1))*ZFACBUO*ZBUOC

! either use entrainment rate or if zero
! use detrainmnet rate as a subsitute for 
! mixing and "pressure" gradient term in upper
! troposphere

          IF(ZDMFEN(JL) > 0.0_dp)THEN
            ZDKEN=MIN(1.0_dp,(1._dp + RG*Z_CWDRAG)*&
             & ZDMFEN(JL)/MAX(RMFCMIN,PMFU(JL,JK+1)))  
          ELSE
            ZDKEN=MIN(1.0_dp,(1._dp + RG*Z_CWDRAG)*&
             & ZDMFDE(JL)/MAX(RMFCMIN,PMFU(JL,JK+1)))  
          ENDIF
          
          ZKINE(JL,JK)=(ZKINE(JL,JK+1)*(1._dp-ZDKEN)+ZDKBUO)/(1._dp+ZDKEN)
          IF(ZBUO(JL,JK) < 0.0_dp.AND.KLAB(JL,JK+1) == 2) THEN
            ZKEDKE=ZKINE(JL,JK)/MAX(1.E-10_dp,ZKINE(JL,JK+1))
            ZKEDKE=MAX(0.0_dp,MIN(1.0_dp,ZKEDKE))
            ZMFUN=SQRT(ZKEDKE)*PMFU(JL,JK+1)
            ZDMFDE(JL)=MAX(ZDMFDE(JL),PMFU(JL,JK+1)-ZMFUN)
            PLUDE(JL,JK)=PLU(JL,JK+1)*ZDMFDE(JL)
            PMFU(JL,JK)=PMFU(JL,JK+1)+ZDMFEN(JL)-ZDMFDE(JL)
          ENDIF
           ! Erase values if below departure level
          IF ( JK > KDPL(JL) ) THEN
            PMFU(JL,JK)=PMFU(JL,JK+1)
            ZKINE(JL,JK)=0.5_dp
          ENDIF
          IF(ZKINE(JL,JK) > 0.0_dp.AND.PMFU(JL,JK) > 0.0_dp) THEN
            KCTOP(JL)=JK
            LLO1(JL)=.TRUE.
          ELSE
            KLAB(JL,JK)=0
            PMFU(JL,JK)=0.0_dp
            ZKINE(JL,JK)=0.0_dp
            ZDMFDE(JL)=PMFU(JL,JK+1)
            PLUDE(JL,JK)=PLU(JL,JK+1)*ZDMFDE(JL)
          ENDIF
          
! store detrainment rates for updraught

          IF ( PMFU(JL,JK+1) > 0.0_dp ) THEN
            PMFUDE_RATE(JL,JK) = ZDMFDE(JL)
            PMFUEN_RATE(JL,JK) = ZDMFEN(JL)
          ENDIF
          
        ENDIF
      ENDIF
    ENDDO

!              CALCULATE PRECIPITATION RATE BY
!              ANALYTIC INTEGRATION OF EQUATION FOR L

    DO JL=KIDIA,KFDIA
      IF(LLO1(JL)) THEN
        IF(LDLAND(JL)) THEN
          ZDNOPRC=ZDLAND(JL)
        ELSE
          ZDNOPRC=1.5E4_dp
        ENDIF
        ZDNOPRC=0.0_dp
        IF(ZPBASE(JL)-PAPH(JL,JK) > ZDNOPRC .AND. PLU(JL,JK) > 1.E-10_dp) THEN
          ZWU=MIN(10._dp,SQRT(2.0_dp*MAX(0.1_dp,ZKINE(JL,JK+1))))
          ZPRCON=ZPRCDGW/(0.75_dp*ZWU)

!           PARAMETERS FOR BERGERON-FINDEISEN PROCESS (T < -5C)

          ZDT=MIN(RTBERCU-RTICECU,MAX(RTBER-PTU(JL,JK),0.0_dp))
          ZCBF=1._dp+Z_CPRC2*SQRT(ZDT)
          ZZCO=ZPRCON*ZCBF
          ZLCRIT=Z_CLCRIT/ZCBF

          ZDFI=PGEOH(JL,JK)-PGEOH(JL,JK+1)
          ZC=(PLU(JL,JK)-ZLUOLD(JL))
          ZD=ZZCO*(1.0_dp-EXP(-(PLU(JL,JK)/ZLCRIT)**2))*ZDFI
          ZINT=EXP(-ZD)
          ZLNEW=ZLUOLD(JL)*ZINT+ZC/ZD*(1.0_dp-ZINT)
          ZLNEW=MAX(0.0_dp,MIN(PLU(JL,JK),ZLNEW))
          ! mz_ht_20070325+
          IF (PTU(jl,jk) > RTT) THEN
            rform(jl,jk) = MAX(0._dp,(plu(jl,jk)-zlnew))
          ELSE
            sform(jl,jk) = MAX(0._dp,(plu(jl,jk)-zlnew))
          ENDIF
          ! mz_ht_20070325-
          ZLNEW=MIN(Z_CLDMAX,ZLNEW)
          ZPRECIP(JL)=MAX(0.0_dp,ZLUOLD(JL)+ZC-ZLNEW)
          PDMFUP(JL,JK)=ZPRECIP(JL)*PMFU(JL,JK)
          ZLRAIN(JL,JK)=ZLRAIN(JL,JK)+ZPRECIP(JL)
          PLU(JL,JK)=ZLNEW
        ENDIF
      ENDIF
    ENDDO

    IF (LPHYLIN) THEN

      DO JL=KIDIA,KFDIA
        IF(LLO1(JL)) THEN
          IF(ZLRAIN(JL,JK) > 0.0_dp) THEN
            ZVW=21.18_dp*ZLRAIN(JL,JK)**0.2_dp
            ZVI=Z_CWIFRAC*ZVW
            ZALFAW=0.5_dp*(TANH(RLPAL1*(PTU(JL,JK)-RLPTRC))+1.0_dp)
            ZVV=ZALFAW*ZVW+(1.0_dp-ZALFAW)*ZVI
            ZROLD=ZLRAIN(JL,JK)-ZPRECIP(JL)
            ZC=ZPRECIP(JL)
            ZWU=MIN(10._dp,SQRT(2.0_dp*MAX(0.1_dp,ZKINE(JL,JK))))
            ZD=ZVV/ZWU
            ZINT=EXP(-ZD)
            ZRNEW=ZROLD*ZINT+ZC/ZD*(1.0_dp-ZINT)
            ZRNEW=MAX(0.0_dp,MIN(ZLRAIN(JL,JK),ZRNEW))
            ZLRAIN(JL,JK)=ZRNEW
          ENDIF
        ENDIF
      ENDDO

    ELSE

      DO JL=KIDIA,KFDIA
        IF(LLO1(JL)) THEN
          IF(ZLRAIN(JL,JK) > 0.0_dp) THEN
            ZVW=21.18_dp*ZLRAIN(JL,JK)**0.2_dp
            ZVI=Z_CWIFRAC*ZVW
            ZALFAW=FOEALFCU(PTU(JL,JK))
            ZVV=ZALFAW*ZVW+(1.0_dp-ZALFAW)*ZVI
            ZROLD=ZLRAIN(JL,JK)-ZPRECIP(JL)
            ZC=ZPRECIP(JL)
            ZWU=MIN(10._dp,SQRT(2.0_dp*MAX(0.1_dp,ZKINE(JL,JK))))
            ZD=ZVV/ZWU
            ZINT=EXP(-ZD)
            ZRNEW=ZROLD*ZINT+ZC/ZD*(1.0_dp-ZINT)
            ZRNEW=MAX(0.0_dp,MIN(ZLRAIN(JL,JK),ZRNEW))
            ZLRAIN(JL,JK)=ZRNEW
          ENDIF
        ENDIF
      ENDDO

    ENDIF

    DO JL=KIDIA,KFDIA
      IF(LLFLAG(JL)) THEN
        PMFUL(JL,JK)=PLU(JL,JK)*PMFU(JL,JK)
        PMFUS(JL,JK)=(RCPD*PTU(JL,JK)+PGEOH(JL,JK))*PMFU(JL,JK)
        PMFUQ(JL,JK)=PQU(JL,JK)*PMFU(JL,JK)
      ENDIF
    ENDDO
    IF(LMFDUDV) THEN
      DO JL=KIDIA,KFDIA
        IF(LLFLAGUV(JL)) THEN
          IF(KTYPE(JL) == 1.OR.KTYPE(JL) == 3) THEN
            IF(ZDMFEN(JL) == 0.0_dp) THEN
              ZZ=3._dp
            ELSE
              ZZ=2.0_dp
            ENDIF
          ELSE
            IF(ZDMFEN(JL) == 0.0_dp) THEN
              ZZ=1.0_dp
            ELSE
              ZZ=0.0_dp
            ENDIF
          ENDIF
          ZDMFEU=ZDMFEN(JL)+ZZ*ZDMFDE(JL)
          ZDMFDU=ZDMFDE(JL)+ZZ*ZDMFDE(JL)
          ZMFUU(JL)=ZMFUU(JL)+ZDMFEU*PUEN(JL,JK)-ZDMFDU*PUU(JL,JK+1)
          ZMFUV(JL)=ZMFUV(JL)+ZDMFEU*PVEN(JL,JK)-ZDMFDU*PVU(JL,JK+1)
          IF(PMFU(JL,JK) > 0.001_dp) THEN
            PUU(JL,JK)=ZMFUU(JL)*(1.0_dp/PMFU(JL,JK))
            PVU(JL,JK)=ZMFUV(JL)*(1.0_dp/PMFU(JL,JK))
          ENDIF
        ENDIF
           ! Erase values if below departure level
        IF ( JK > KDPL(JL) ) THEN
          PUU(JL,JK)=PUEN(JL,JK-1)
          PVU(JL,JK)=PVEN(JL,JK-1)
        ENDIF
      ENDDO
    ENDIF

  ENDIF
ENDDO

!----------------------------------------------------------------------

!     5.           FINAL CALCULATIONS 
!                  ------------------

DO JL=KIDIA,KFDIA
  IF(KCTOP(JL) == -1) LDCUM(JL)=.FALSE.
  KCBOT(JL)=MAX(KCBOT(JL),KCTOP(JL))
ENDDO

IF (LHOOK) CALL DR_HOOK('CUASCN',1,ZHOOK_HANDLE)
END SUBROUTINE CUASCN
!==============================================================================

SUBROUTINE CUBASEN &
 & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
 & PTENH,    PQENH,    PGEOH,    PAPH,&
 & PQHFL,    PAHFS,    PSSTRU,   PSSTRV,   PWN,&
 & PTEN,     PQEN,     PGEO,&
 & PUEN,     PVEN,&
 & PTU,      PQU,      PLU,      PUU,      PVU ,  PWUBASE,&
 & KLAB,     LDCUM,    LDSC,     KCBOT,    KBOTSC,&
 & KCTOP,    KDPL,     PCAPE  )  

!          THIS ROUTINE CALCULATES CLOUD BASE FIELDS
!          CLOUD BASE HEIGHT AND CLOUD TOP HEIGHT

!          A. Pier Siebesma   KNMI ********      
!          modified C Jakob (ECMWF) (01/2001) 
!          modified P Bechtold (ECMWF) (08/2002) 
!          (include cycling over levels to find unstable departure/base level+
!           mixed layer properties +w Trigger)

!          PURPOSE.
!          --------
!          TO PRODUCE CLOUD BASE AND CLOUD TOP VALUES FOR CU-PARAMETRIZATION

!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!          INPUT ARE ENVIRONM. VALUES OF T,Q,P,PHI AT HALF LEVELS.
!          IT RETURNS CLOUD FIELDS VALUES AND FLAGS AS FOLLOWS;
!                 KLAB=0 FOR STABLE LAYERS
!                 KLAB=1 FOR SUBCLOUD LEVELS
!                 KLAB=2 FOR CLOUD LEVELS LEVEL

!          METHOD.
!          --------
!          LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD TOP
!          (ENTRAINING PLUME, WITH ENTRAINMENT PROPORTIONAL TO (1/Z))

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KTDIA*        START OF THE VERTICAL LOOP
!    *KLEV*         NUMBER OF LEVELS

!    INPUT PARAMETERS (REAL):

! not used at the moment because we want to use linear intepolation
! for fields on the half levels.

!    *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS           K
!    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS      KG/KG

!    *PQHFL*        MOISTURE FLUX (EXCEPT FROM SNOW EVAP.)        KG/(SM2)
!    *PAHFS*        SENSIBLE HEAT FLUX                            W/M2
!    *PSSTRU*       KINEMATIC surface U-MOMENTUM FLUX             (M/S)^2
!    *PSSTRV*       KINEMATIC surface V-MOMENTUM FLUX             (M/S)^2
!    *PWN*          NORMALIZED LARGE-SCALE VERTICAL VELOCITY      (M/S)
!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS             PA
!    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
!    *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1)  KG/KG
!    *PGEO*         GEOPOTENTIAL                                  M2/S2
!    *PUEN*         PROVISIONAL ENVIRONMENT U-VELOCITY (T+1)       M/S
!    *PVEN*         PROVISIONAL ENVIRONMENT V-VELOCITY (T+1)       M/S
!    *PQHFL*        MOISTURE FLUX (EXCEPT FROM SNOW EVAP.)        KG/(SM2)
!    *PAHFS*        SENSIBLE HEAT FLUX                            W/M2

!    UPDATED PARAMETERS (REAL):

!    *PTU*          TEMPERATURE IN UPDRAFTS                         K
!    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                    KG/KG
!    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS              KG/KG
!    *PUU*          U-VELOCITY IN UPDRAFTS                         M/S
!    *PVU*          V-VELOCITY IN UPDRAFTS                         M/S

!    UPDATED PARAMETERS (INTEGER):

!    *KLAB*         FLAG KLAB=1 FOR SUBCLOUD LEVELS
!                        KLAB=2 FOR CLOUD LEVELS

!    OUTPUT PARAMETERS (LOGICAL):

!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS 
!    *LDSC*         FLAG: .TRUE. IF BL-CLOUDS EXIST

!    OUTPUT PARAMETERS (INTEGER):

!    *KCBOT*       CLOUD BASE LEVEL !    
!    *KCTOP*       CLOUD TOP LEVEL = HEIGHEST HALF LEVEL 
!                  WITH A NON-ZERO CLOUD UPDRAFT.
!    *KBOTSC*      CLOUD BASE LEVEL OF BL-CLOUDS
!    *KDPL*        DEPARTURE LEVEL
!    *PCAPE*       PSEUDOADIABATIQUE max CAPE (J/KG)

!          EXTERNALS
!          ---------
!          *CUADJTQ* FOR ADJUSTING T AND Q DUE TO CONDENSATION IN ASCENT

!          MODIFICATIONS
!          -------------
!             92-09-21 : Update to Cy44      J.-J. MORCRETTE
!             02-11-02 : Use fixed last possible departure level and 
!                        last updraft computation level for bit-reproducibility
!                                            D.Salmond &  J. Hague
!             03-07-03 : Tuning for p690     J. Hague
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!----------------------------------------------------------------------

USE MESSY_CONVECT_ECMWF_PARAM   , ONLY : RCPD, RETV, RD, RG, RLVTT, RLSTT, RTT,      & 
                                         R2ES, R3LES, R3IES, R4LES,                  &
                                         R4IES, R5LES, R5IES, R5ALVCP, R5ALSCP,      &
                                         RALVDCP, RALSDCP, RALFDCP, RTWAT, RTICE,    &
                                         RTICECU, RTWAT_RTICECU_R, RTWAT_RTICE_R,    &
                                         LMFDUDV, ENTRPEN, RDEPTHS, NJKT1, NJKT2,    & !from YOECUMF
                                         RLMIN,                                      & !from YOECLDP
                                         RKAP,                                       & !from YOEVDF
                                         LHOOK, DR_HOOK,                             & !from yomhook
                                         FOEALFCU, FOEEWM, FOEALFA

IMPLICIT NONE

INTEGER,INTENT(IN)    :: KLON 
INTEGER,INTENT(IN)    :: KLEV 
INTEGER,INTENT(IN)    :: KIDIA 
INTEGER,INTENT(IN)    :: KFDIA 
INTEGER               :: KTDIA ! Argument NOT used
REAL(dp)   ,INTENT(IN)    :: PTENH(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PQENH(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PGEOH(KLON,KLEV+1) 
REAL(dp)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) 
REAL(dp)   ,INTENT(IN)    :: PQHFL(KLON,KLEV+1) 
REAL(dp)   ,INTENT(IN)    :: PAHFS(KLON,KLEV+1) 
REAL(dp)                  :: PSSTRU(KLON) ! Argument NOT used
REAL(dp)                  :: PSSTRV(KLON) ! Argument NOT used
REAL(dp)                  :: PWN(KLON,KLEV) ! Argument NOT used
REAL(dp)   ,INTENT(IN)    :: PTEN(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PQEN(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PGEO(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PUEN(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PVEN(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PTU(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PQU(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PLU(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PUU(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PVU(KLON,KLEV) 
REAL(dp)   ,INTENT(OUT)   :: PWUBASE(KLON) 
INTEGER,INTENT(INOUT) :: KLAB(KLON,KLEV) 
LOGICAL           ,INTENT(INOUT) :: LDCUM(KLON) 
LOGICAL           ,INTENT(OUT)   :: LDSC(KLON) 
INTEGER,INTENT(INOUT) :: KCBOT(KLON) 
INTEGER,INTENT(OUT)   :: KBOTSC(KLON) 
INTEGER,INTENT(OUT)   :: KCTOP(KLON) 
INTEGER,INTENT(OUT)   :: KDPL(KLON) 
REAL(dp)   ,INTENT(OUT)   :: PCAPE(KLON) 
INTEGER ::  ICTOP(KLON),            ICBOT(KLON),&
 & IBOTSC(KLON),           ILAB(KLON,KLEV),&
 & IDPL(KLON)  

!             LOCAL STORAGE
!             ----- -------

LOGICAL ::         LL_LDBASE(KLON),&
 & LLGO_ON(KLON),&
 & LLDEEP(KLON),    LLDCUM(KLON), &!*UPG PB
 & LLDSC(KLON),     LLFIRST(KLON)  
LOGICAL ::     LLRESET,        LLRESETJL(KLON)

INTEGER :: ICALL, IK, IKB, IS, JK, JL, JKK, JKT1, JKT2, JKT, JKB !*UPG PB

REAL(dp)    :: ZS(KLON,KLEV),&
 & ZSENH(KLON,KLEV+1),&
 & ZQENH(KLON,KLEV+1),&
 & ZSUH (KLON,KLEV),&
 & ZWU2H(KLON,KLEV),&
 & ZBUOH(KLON,KLEV)  
REAL(dp) :: ZQOLD(KLON),ZPH(KLON)
REAL(dp) :: ZMIX(KLON)
REAL(dp) :: ZDZ(KLON),zcbase(klon)

REAL(dp) ::    ZLU(KLON,KLEV),   ZQU(KLON,KLEV),&
 & ZTU(KLON,KLEV), &
 & ZUU(KLON,KLEV),   ZVU(KLON,KLEV)  

REAL(dp) :: ZCAPE(KLON,KLEV) ! local for CAPE at every departure level

REAL(dp) :: ZBUOF, ZZ, ZC2, ZEPSADD
REAL(dp) :: ZRHO      ! DENSITY AT SURFACE (KG/M^3) 
REAL(dp) :: ZKHVFL    ! SURFACE BUOYANCY FLUX (K M/S)
REAL(dp) :: ZWS       ! SIGMA_W            AT LOWEST MODEL HALFLEVEL (M/S)
REAL(dp) :: ZQEXC     ! HUMIDITY    EXCESS AT LOWEST MODEL HALFLEVEL (#)
REAL(dp) :: ZTEXC     ! TEMPERATURE EXCESS AT LOWEST MODEL HALFLEVEL (K)
REAL(dp) :: ZEPS      ! FRACTIONAL ENTRAINMENT RATE   [M^-1]
REAL(dp) :: ZTVENH    ! ENVIRONMENT VIRTUAL TEMPERATURE AT HALF LEVELS (K)  
REAL(dp) :: ZTVUH     ! UPDRAFT VIRTUAL TEMPERATURE AT HALF LEVELS     (K)
REAL(dp) :: ZLGLAC    ! UPDRAFT LIQUID WATER FROZEN IN ONE LAYER
REAL(dp) :: zqsu, ZCOR, zdq, zalfaw, zfacw, zfaci, zfac,&
 & zesdp, zdqsdt, zdtdp, zdp,zpdifftop,zpdiffbot,ZSF,ZQF,zaw,zbw  
REAL(dp) :: ZTVEN1, ZTVEN2, ZTVU1, ZTVU2 ! pseudoadiabatique T_v
REAL(dp) :: ZDTVTRIG(KLON) ! virtual temperatures
REAL(dp) :: ZWORK1, ZWORK2! work arrays for T and w perturbations
REAL(dp) :: ZRCPD, ZRG, ZTMP
REAL(dp) :: ZTINY
REAL(dp) :: ZHOOK_HANDLE

INTRINSIC :: EXP, MAXVAL, MIN, SQRT, TINY
!----------------------------------------------------------------------
!     0.           INITIALIZE CONSTANTS AND FIELDS
!                  -------------------------------
!----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CUBASEN',0,ZHOOK_HANDLE)
ZC2    = 0.55_dp
ZAW    = 1.0_dp
ZBW    = 1.0_dp
ZEPSADD= 1.E-4_dp

DO JL=KIDIA,KFDIA
  PWUBASE(JL)=0.0_dp
  LLGO_ON(JL)=.TRUE.
  LLFIRST(JL)=.TRUE.
  KDPL(JL)=KLEV
ENDDO

! Set last possible departure level and last updraft computation level
! NOT Bit-reproducible
!* change to operations
DO JK=KLEV+1,2,-1
  IF((PAPH(KIDIA,KLEV+1)-PAPH(KIDIA,JK)) < 350.E2_dp ) JKT1=JK
  IF(PAPH(KIDIA,JK) > 60.E2_dp ) JKT2=JK
END DO

!JKT1=NJKT1
!JKT2=NJKT2
!* change to operations

ZRG=1.0_dp/RG
ZRCPD=1.0_dp/RCPD
! The function TINY is standard Fortran90, but is not available on all systems
!ZTINY=TINY(ZRG)
ZTINY=1.0E-35_dp

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZTU(JL,JK) = PTU(JL,JK)
    ZQU(JL,JK) = PQU(JL,JK)
    ZLU(JL,JK) = PLU(JL,JK)
    ZUU(JL,JK) = PUU(JL,JK)
    ZVU(JL,JK) = PVU(JL,JK)
    ILAB(JL,JK)= KLAB(JL,JK)
    ZCAPE(JL,JK)= 0.0_dp
  ENDDO
ENDDO

!----------------------------------------------------------------------
!       -----------------------------------------------------------
!       1.1  PREPARE FIELDS ON HALF LEVELS BY LINEAR INTERPOLATION
!             OF SPECIFIC HUMIDITY AND STATIC ENERGY
!       -----------------------------------------------------------

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZWU2H(JL,JK)=0.0_dp
    ZS   (JL,JK) = RCPD*PTEN(JL,JK) + PGEO(JL,JK)
    ZQENH(JL,JK) = PQENH(JL,JK)
    ZSENH(JL,JK) = RCPD*PTENH(JL,JK)+PGEOH(JL,JK)
  ENDDO
ENDDO

DO JKK=KLEV,JKT1,-1 ! Big external loop for level testing:
                    ! find first departure level that produces deepest cloud top
                    ! or take surface level for shallow convection and Sc
   !
   !        ---------------------------------------------------------
   !        1.2    INITIALISE FIELDS AT DEPARTURE HALF MODEL LEVEL
   !        ---------------------------------------------------------
   !
  IS=0
  DO JL=KIDIA,KFDIA
    IF (LLGO_ON(JL)) THEN
      IS=IS+1
      IDPL(JL)    =JKK      ! departure level
      ICBOT  (JL) =JKK      ! cloud base level for convection, (-1 if not found)
      IBOTSC (JL) =KLEV-1   ! sc    base level for sc-clouds , (-1 if not found)
      ICTOP(JL)   =KLEV-1   ! cloud top for convection (-1 if not found)
      LLDCUM(JL)  =.FALSE.  ! on exit: true if cloudbase=found
      LLDSC (JL)  =.FALSE.  ! on exit: true if cloudbase=found
      LL_LDBASE(JL)   =.FALSE. ! on exit: true if cloudbase=found
      ZDTVTRIG(JL) =0.0_dp
      ZUU(JL,JKK) =PUEN(JL,JKK)*(PAPH(JL,JKK+1)-PAPH(JL,JKK))
      ZVU(JL,JKK) =PVEN(JL,JKK)*(PAPH(JL,JKK+1)-PAPH(JL,JKK))
    ENDIF 
  ENDDO

  IF(IS /= 0) THEN

    IF(JKK == KLEV) THEN

      DO JL=KIDIA,KFDIA
        IF (LLGO_ON(JL)) THEN
          ZRHO  = PAPH(JL,JKK+1)/(RD*(PTEN(JL,JKK)*(1._dp+RETV*PQEN(JL,JKK))))
          ZKHVFL= (PAHFS(JL,JKK+1)/RCPD+RETV*PTEN(JL,JKK)*PQHFL(JL,JKK+1))/ZRHO
        ! ZUST  = MAX(SQRT(PSSTRU(JL)**2 + PSSTRV(JL)**2),REPUST)
        ! ZWS=ZUST**3._dp- 1.5_dp*RKAP*ZKHVFL*PGEOH(JL,KLEV)/PTEN(JL,KLEV)
          ZWS=0.001_dp- 1.5_dp*RKAP*ZKHVFL*PGEOH(JL,KLEV)/PTEN(JL,KLEV)
          IF(ZWS >= ZTINY) THEN
            ZWS=1.2_dp*ZWS**.3333_dp
            ILAB(JL,JKK)= 1
            ZTEXC     = MAX(-1.5_dp*PAHFS(JL,JKK+1)/(ZRHO*ZWS*RCPD),0.0_dp)
            ZQEXC     = MAX(-1.5_dp*PQHFL(JL,JKK+1)/(ZRHO*ZWS),0.0_dp)
            ZQU (JL,JKK) = ZQENH(JL,JKK) + ZQEXC
            ZSUH (JL,JKK) = ZSENH(JL,JKK) + RCPD*ZTEXC
            ZTU (JL,JKK) = (ZSENH(JL,JKK)-PGEOH(JL,JKK))/RCPD + ZTEXC
            ZWU2H(JL,JKK) = ZWS**2
        !
        !  determine buoyancy at lowest half level
        !
            ZTVENH            = (1.0_dp+RETV*ZQENH(JL,JKK)) &
             & *(ZSENH(JL,JKK)-PGEOH(JL,JKK))/RCPD  
            ZTVUH             = (1.0_dp+RETV*ZQU(JL,JKK))*ZTU(JL,JKK)
            ZBUOH(JL,JKK) = (ZTVUH-ZTVENH)*RG/ZTVENH
          ELSE
            LLGO_ON(JL)=.FALSE.      ! non-convective point
          ENDIF
        ENDIF
      ENDDO
   
    ELSE

      DO JL=KIDIA,KFDIA
        IF (LLGO_ON(JL)) THEN
          ZRHO  = PAPH(JL,JKK+1)/(RD*(PTEN(JL,JKK)*(1._dp+RETV*PQEN(JL,JKK))))
          ILAB(JL,JKK)= 1
          ZTEXC=.2_dp
          ZQEXC=1.E-4_dp
          ZQU (JL,JKK) = ZQENH(JL,JKK) + ZQEXC
          ZSUH (JL,JKK) = ZSENH(JL,JKK) + RCPD*ZTEXC
          ZTU (JL,JKK) = (ZSENH(JL,JKK)-PGEOH(JL,JKK))*ZRCPD + ZTEXC
         ! construct mixed layer for parcels emanating in lowest 60 hPa
          IF (PAPH(JL,KLEV+1)-PAPH(JL,JKK-1)<60.E2_dp) THEN
            ZQU(JL,JKK) =0.0_dp
            ZSUH(JL,JKK)=0.0_dp
            ZWORK1      =0.0_dp
            DO JK=JKK+1,JKK-1,-1
              if( ZWORK1 < 50.E2_dp ) then
                ZWORK2=PAPH(JL,JK)-PAPH(JL,JK-1)
                ZWORK1      =ZWORK1+ZWORK2
                ZQU(JL,JKK) =ZQU(JL,JKK) +ZQENH(JL,JK)*ZWORK2
                ZSUH(JL,JKK)=ZSUH(JL,JKK)+ZSENH(JL,JK)*ZWORK2
              endif
            ENDDO
            ZQU(JL,JKK) =ZQU(JL,JKK) /ZWORK1+ZQEXC
            ZSUH(JL,JKK)=ZSUH(JL,JKK)/ZWORK1+RCPD*ZTEXC
            ZTU(JL,JKK) =(ZSUH(JL,JKK)-PGEOH(JL,JKK))/RCPD+ZTEXC
          ENDIF
          ZWU2H(JL,JKK) = 1.0_dp
      !
      !  determine buoyancy at lowest half level
      !
          ZTVENH            = (1.0_dp+RETV*ZQENH(JL,JKK)) &
           & *(ZSENH(JL,JKK)-PGEOH(JL,JKK))/RCPD  
          ZTVUH             = (1.0_dp+RETV*ZQU(JL,JKK))*ZTU(JL,JKK)
          ZBUOH(JL,JKK) = (ZTVUH-ZTVENH)*RG/ZTVENH
        ENDIF
      ENDDO
   
    ENDIF

  ENDIF
   
   !----------------------------------------------------------------------
   !     2.0          DO ASCENT IN SUBCLOUD AND LAYER,
   !                  CHECK FOR EXISTENCE OF CONDENSATION LEVEL,
   !                  ADJUST T,Q AND L ACCORDINGLY IN *CUADJTQ*,
   !                  CHECK FOR BUOYANCY AND SET FLAGS
   !                  -------------------------------------
   !       ------------------------------------------------------------
   !        1.2  DO THE VERTICAL ASCENT UNTIL VELOCITY BECOMES NEGATIVE
   !       ------------------------------------------------------------
  DO JK=JKK-1,JKT2,-1
    IS=0

    if(jkk==klev) then ! 1/z mixing for shallow

      DO JL=KIDIA,KFDIA
        IF (LLGO_ON(JL)) THEN
          IS         = IS+1
          ZDZ(JL)        = (PGEOH(JL,JK) - PGEOH(JL,JK+1))*ZRG
          if ((PGEO(JL,JK)*ZRG + ZDZ(jl)).lt.1e4*tiny(0.0_dp)) &
            print*,(PGEO(JL,JK)*ZRG + ZDZ(jl)),PGEO(JL,JK)*ZRG, ZDZ(jl),PGEO(JL,JK),jl,jk 
          ZEPS       = ZC2/(PGEO(JL,JK)*ZRG + ZDZ(jl)) + ZEPSADD
          ZMIX(JL)       = 0.5_dp*ZDZ(JL)*ZEPS
          ZQF = (PQENH(JL,JK+1) + PQENH(JL,JK))*0.5_dp
          ZSF = (ZSENH(JL,JK+1) + ZSENH(JL,JK))*0.5_dp
          ZTMP = 1.0_dp/(1.0_dp+ZMIX(JL))
          ZQU(JL,JK)= (ZQU(JL,JK+1)*(1.0_dp-ZMIX(JL))&
           & +2.0_dp*ZMIX(jl)*ZQF) * ZTMP  
          ZSUH (JL,JK)= (ZSUH(JL,JK+1)*(1.0_dp-ZMIX(JL))&
           & +2.0_dp*ZMIX(jl)*ZSF) * ZTMP  
          ZQOLD(JL)  = ZQU(JL,JK)
          ZTU (JL,JK) = (ZSUH(JL,JK)-PGEOH(JL,JK))*ZRCPD
          ZPH  (JL)    = PAPH(JL,JK)
        ENDIF
      ENDDO
     
    else

      DO JL=KIDIA,KFDIA
        IF (LLGO_ON(JL)) THEN
          IS         = IS+1
          ZDZ(JL)        = (PGEOH(JL,JK) - PGEOH(JL,JK+1))*ZRG
          ZQF = (PQENH(JL,JK+1) + PQENH(JL,JK))*0.5_dp
          ZSF = (ZSENH(JL,JK+1) + ZSENH(JL,JK))*0.5_dp
        ! ZMIX(JL)=5.E-2_dp  
          ZMIX(JL)=ENTRPEN*(PGEOH(JL,JK)-PGEOH(JL,JK+1))/RG
          ZQU(JL,JK)= ZQU(JL,JK+1)*(1.0_dp-ZMIX(JL))+ ZQF*ZMIX(JL)
          ZSUH(JL,JK)= ZSUH(JL,JK+1)*(1.0_dp-ZMIX(JL))+ ZSF*ZMIX(JL)
          ZQOLD(JL)  = ZQU(JL,JK)
          ZTU (JL,JK) = (ZSUH(JL,JK)-PGEOH(JL,JK))*ZRCPD
          ZPH  (JL)    = PAPH(JL,JK)
        ENDIF
      ENDDO
  
    endif
       
    IF (IS == 0) EXIT
     
    IK=JK
    ICALL=1
     
    CALL CUADJTQ &
     & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
     & IK,&
     & ZPH,      ZTU,      ZQU,      LLGO_ON,   ICALL)  
   
!DIR$ IVDEP
!OCL NOVREC
   
    DO JL=KIDIA,KFDIA
      IF(LLGO_ON(JL)) THEN
   
   ! add condensation to water
   
        ZDQ=MAX(ZQOLD(JL)-ZQU(JL,JK),0.0_dp)
        ZLU(JL,JK)=ZLU(JL,JK+1)+ZDQ

   ! freezing
   
        ZLGLAC=ZDQ*((1.0_dp-FOEALFCU(ZTU(JL,JK)))-&
         & (1.0_dp-FOEALFCU(ZTU(JL,JK+1))))  
   
   ! pseudo-microphysics
   
        if(jkk==klev) then  ! no precip for shallow
          ZLU(JL,JK)=MIN(ZLU(JL,JK),5.E-3_dp)
   !* chose a more pseudo-adiabatic formulation as original overestimates
   !* water loading efect and therefore strongly underestimates cloud thickness
        else 
          ZLU(JL,JK)=0.5_dp*ZLU(JL,JK) 
        endif
   
   ! update dry static energy after condensation + freezing
   
        ZSUH(JL,JK)    = RCPD*(ZTU(JL,JK)+RALFDCP*ZLGLAC)+PGEOH(JL,JK)
         
   ! Buoyancy on half and full levels
            
        ZTVUH           = (1.0_dp+RETV*ZQU(JL,JK)-ZLU(JL,JK))*ZTU(JL,JK)&
         & +RALFDCP*ZLGLAC  
        ZTVENH          = (1.0_dp+RETV*ZQENH(JL,JK)) &
         & *(ZSENH(JL,JK)-PGEOH(JL,JK))*ZRCPD  
        ZBUOH(JL,JK)   = (ZTVUH-ZTVENH)*RG/ZTVENH
        ZBUOF          = (ZBUOH(JL,JK) + ZBUOH(JL,JK+1))*0.5_dp
   
   ! solve kinetic energy equation
   
        ZTMP=1.0_dp/(1.0_dp+2.0_dp*ZBW*ZMIX(jl))
        ZWU2H(JL,JK) = (ZWU2H(JL,JK+1)*(1.0_dp-2.0_dp*ZBW*ZMIX(jl))&
         & +2.0_dp*ZAW*ZBUOF*ZDZ(jl)) * ZTMP  
   
   ! compute pseudoadiabatique CAPE for diagnostics
   
        ZTVU2 = ZTU(JL,JK)  *(1.0_dp+RETV*ZQU(JL,JK))
        ZTVEN2= PTENH(JL,JK)*(1.0_dp+RETV*PQENH(JL,JK))
        IF (JK == JKK-1) THEN
          ZTVU1  = ZTVU2
          ZTVEN1 = ZTVEN2
        ENDIF
        ZBUOF = (ZTVU2+ZTVU1-ZTVEN1-ZTVEN2)/ZTVEN2
        ZBUOF = ZBUOF*ZDZ(JL)*RG
        ZCAPE(JL,JKK)  = ZCAPE(JL,JKK) + MAX(0.0_dp,ZBUOF)
        ZTVU1=ZTVU2
        ZTVEN1=ZTVEN2
   
   ! first layer with liquid water - find exact cloud base
   
        IF(ZLU(JL,JK) >0.0_dp.AND.ILAB(JL,JK+1)==1) THEN
           
          IK=JK+1
          ZQSU=FOEEWM(ZTU(JL,IK))/PAPH(JL,IK)
          ZQSU=MIN(0.5_dp,ZQSU)
          ZCOR=1.0_dp/(1.0_dp-RETV  *ZQSU)
          ZQSU=ZQSU*ZCOR
          ZDQ=MIN(0._dp,ZQU(JL,IK)-ZQSU)
          ZALFAW=FOEALFA(ZTU(JL,IK))
          ZFACW=R5LES/((ZTU(JL,IK)-R4LES)**2)
          ZFACI=R5IES/((ZTU(JL,IK)-R4IES)**2)
          ZFAC=ZALFAW*ZFACW+(1._dp-ZALFAW)*ZFACI
          ZESDP=FOEEWM(ZTU(JL,IK))/PAPH(JL,IK)
          ZCOR=1.0_dp/(1.0_dp-RETV*ZESDP)
          ZDQSDT=ZFAC*ZCOR*ZQSU
          ZDTDP=RD*ZTU(JL,IK)/(RCPD*PAPH(JL,IK))
          ZDP=ZDQ/(ZDQSDT*ZDTDP)
          ZCBASE(JL)=PAPH(JL,IK)+ZDP
           
   ! chose nearest half level as cloud base
   
          ZPDIFFTOP=ZCBASE(JL)-PAPH(JL,JK)
          ZPDIFFBOT=PAPH(JL,JK+1)-ZCBASE(JL)
           
          IF(ZPDIFFTOP > ZPDIFFBOT.AND.ZWU2H(JL,JK+1)>0.0_dp) THEN
            JKB=MIN(KLEV-1,JK+1)
            ILAB(JL,JKB)=2 !*UPG
            ILAB(JL,JK)=2
            LL_LDBASE(JL) =.TRUE.
            LLDSC(JL)   =.TRUE.
            IBOTSC(JL) =JKB
            ICBOT(JL)  =JKB
            ZLU(JL,JK+1) = RLMIN
          ELSEIF(ZPDIFFTOP <= ZPDIFFBOT.AND.ZWU2H(JL,JK)>0.0_dp) THEN
            ILAB(JL,JK)=2
            LL_LDBASE(JL) =.TRUE.
            LLDSC(JL)   =.TRUE.
            IBOTSC(JL) =JK
            ICBOT(JL)  =JK
          ENDIF
          JKB=ICBOT(JL)
         ! vertical velocity trigger for deep conv
         ! was only used for test
          !IF( JKK<KLEV .AND. ZWU2H(JL,JKB)>0._dp ) THEN 
          !  ZWORK1=0.5_dp*(PWN(JL,JKB)+&
          !     PWN(JL,JKK)*PGEOH(JL,JK)/PGEOH(JL,JKK)) ! mean w in PBL
          !  ZTVENDPL  =PTENH(JL,JKK)*(1._dp+RETV*PQENH(JL,JKK))
          !  ZTVENBASE =PTENH(JL,JKB)*(1._dp+RETV*PQENH(JL,JKB))
          !  ZTVUBASE  =ZTU(JL,JKB)*(1._dp+RETV*ZQU(JL,JKB))
          !  ZWORK2    =SIGN(1._dp,ZWORK1)
          !  ZWORK1    =6._dp*ZWORK2*ABS(ZWORK1)**.3333_dp
          !  ZWORK2    =ZWORK1*(PGEOH(JL,JK)-PGEOH(JL,JKK))&
          !             &/(ZTVENDPL+ZTVENBASE)
          !  ZDTVTRIG(JL)=ZTVUBASE-ZTVENBASE+ZWORK1
          !  ZWU2H(JL,JK) =1._dp+0.5_dp*ZWORK2  ! erase previous values
          !  ZWU2H(JL,JKB)=1._dp+0.5_dp*ZWORK2
          !ENDIF
   
        ENDIF
   
   ! decide on presence of convection, cloud base and cloud top based on
   ! kinetic energy
   
        IF (ZWU2H(JL,JK) < 0.0_dp) THEN
          LLGO_ON(JL) = .FALSE.             
          IF (ZLU(JL,JK+1)>0.0_dp) THEN
            ICTOP(JL)   = JK
            LLDCUM(JL)   = .TRUE.
          ELSE
            LLDCUM(JL)   = .FALSE.
          ENDIF
        ELSE
          IF (ZLU(JL,JK)>0.0_dp) then
            ILAB(JL,JK) = 2
          ELSE
            ILAB(JL,JK) = 1
          ENDIF
        ENDIF
      ENDIF
    ENDDO
   
    IF(LMFDUDV.AND.JKK==KLEV) THEN
      DO JL=KIDIA,KFDIA
        IF(.NOT.LL_LDBASE(JL).AND.LLGO_ON(JL)) THEN
          ZUU(JL,JKK)=ZUU(JL,JKK)+PUEN(JL,JK)*(PAPH(JL,JK+1)-PAPH(JL,JK))
          ZVU(JL,JKK)=ZVU(JL,JKK)+PVEN(JL,JK)*(PAPH(JL,JK+1)-PAPH(JL,JK))
        ENDIF
      ENDDO
    ENDIF
   
!     IF (IS == 0) EXIT
  ENDDO
   
  IF( JKK==KLEV) THEN
      ! set values for departure level for PBL clouds = first model level
    DO JL=KIDIA,KFDIA
      LDSC(JL)  = LLDSC(JL)
      IF(LDSC(JL)) THEN
        KBOTSC(JL)= IBOTSC(JL)
      ELSE
        KBOTSC(JL)=-1
      ENDIF
    
      LLGO_ON(JL) = .FALSE.
      JKT=ICTOP(JL)
      JKB=ICBOT(JL)
      LLDEEP(JL)=PAPH(JL,JKB)-PAPH(JL,JKT)>RDEPTHS
      IF(LLDEEP(JL)) LLDCUM(JL)=.FALSE. ! no deep allowed for KLEV
      lldeep(jl)=.false. ! for deep convection start only at level KLEV-1
                            ! and form mixed layer, so go on
      ! test further for deep convective columns as not yet found
      IF ( LLDEEP(JL) ) LLFIRST(JL)=.FALSE.
      LLGO_ON(JL) = .NOT.LLDEEP(JL)
      IF(LLDCUM(JL)) THEN
        KCBOT(JL)= ICBOT(JL)
        KCTOP(JL)= ICTOP(JL)
        KDPL(JL)  = IDPL(JL)
        LDCUM(JL) = LLDCUM(JL)
        PWUBASE(JL)=SQRT(MAX(ZWU2H(JL,JKB),0.0_dp))
      ELSE
        KCTOP(JL)=-1
        KCBOT(JL)=-1
        KDPL(JL) =KLEV-1
        LDCUM(JL)=.FALSE.
        PWUBASE(JL)=0.0_dp
      ENDIF
    ENDDO
    DO JK=KLEV,1,-1
      DO JL=KIDIA,KFDIA
        JKT=ICTOP(JL)
        IF ( JK>=JKT ) THEN
          KLAB(JL,JK)=ILAB(JL,JK)
          PTU(JL,JK)=ZTU(JL,JK)
          PQU(JL,JK)=ZQU(JL,JK)
          PLU(JL,JK)=ZLU(JL,JK)
        ENDIF
      ENDDO
    ENDDO
    IF(LMFDUDV) THEN
      DO JL=KIDIA,KFDIA
        IF(LDCUM(JL)) THEN
          IKB=KCBOT(JL)
          ZZ=1.0_dp/(PAPH(JL,JKK+1)-PAPH(JL,IKB))
          PUU(JL,JKK)=ZUU(JL,JKK)*ZZ
          PVU(JL,JKK)=ZVU(JL,JKK)*ZZ
!        ELSE
!          PUU(JL,JKK)=PUEN(JL,JKK-1)
!          PVU(JL,JKK)=PVEN(JL,JKK-1)
        ENDIF
      ENDDO
    ENDIF
  ENDIF
   
  IF( JKK < KLEV ) THEN
    LLRESET=.FALSE.
    DO JL=KIDIA,KFDIA
      IF ( .NOT.LLDEEP(JL) ) THEN
        JKT=ICTOP(JL)
        JKB=ICBOT(JL)
           ! test on cloud thickness and buoyancy
        LLDEEP(JL)=PAPH(JL,JKB)-PAPH(JL,JKT)>=RDEPTHS
       ! LLDEEP(JL)=PAPH(JL,JKB)-PAPH(JL,JKT)>=RDEPTHS &
       !   &.AND. ZDTVTRIG(JL)>0._dp
      ENDIF
      LLRESETJL(JL)=LLDEEP(JL).AND.LLFIRST(JL)
      LLRESET=LLRESET.OR.LLRESETJL(JL)
    ENDDO


    IF(LLRESET) THEN
      DO JK=KLEV,1,-1
        DO JL=KIDIA,KFDIA
         ! keep first departure level that produces deep cloud
!          IF ( LLDEEP(JL) .AND. LLFIRST(JL) ) THEN 
          IF ( LLRESETJL(JL) ) THEN 
            JKT=ICTOP(JL)
            JKB=IDPL(JL)
            IF ( JK<=JKB .AND. JK>=JKT ) THEN
              KLAB(JL,JK)=ILAB(JL,JK)
              PTU(JL,JK)=ZTU(JL,JK)
              PQU(JL,JK)=ZQU(JL,JK)
              PLU(JL,JK)=ZLU(JL,JK)
            ELSE 
              KLAB(JL,JK)=1
              PTU(JL,JK)=PTENH(JL,JK)
              PQU(JL,JK)=PQENH(JL,JK)
              PLU(JL,JK)=0.0_dp
            ENDIF
            IF ( JK<JKT ) KLAB(JL,JK)=0
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    DO JL=KIDIA,KFDIA
      IF ( LLDEEP(JL) .AND. LLFIRST(JL) ) THEN
        KDPL(JL)  = IDPL(JL)
        KCTOP(JL) = ICTOP(JL)
        KCBOT(JL) = ICBOT(JL)
        LDCUM(JL) = LLDCUM(JL)
        LDSC(JL)  = .FALSE.
        KBOTSC(JL)= -1
        JKB=KCBOT(JL)
        PWUBASE(JL)=SQRT(MAX(ZWU2H(JL,JKB),0.0_dp))
        IF(LMFDUDV) THEN
              ! IKB=KCBOT(JL)
              ! ZZ=1._dp/(PAPH(JL,JKK+1)-PAPH(JL,IKB))
              ! PUU(JL,JKK)=ZUU(JL,JKK)*ZZ
              ! PVU(JL,JKK)=ZVU(JL,JKK)*ZZ
          puu(jl,jkk)=puen(jl,jkk-1)
          pvu(jl,jkk)=pven(jl,jkk-1)
        ELSE
          PUU(JL,JKK)=PUEN(JL,JKK-1)
          PVU(JL,JKK)=PVEN(JL,JKK-1)
        ENDIF
        LLFIRST(JL)=.FALSE.
      ENDIF
      LLGO_ON(JL) = .NOT.LLDEEP(JL)
    ENDDO
  ENDIF

ENDDO ! end of big loop for search of departure level     

      ! chose maximum CAPE value
DO JL=KIDIA,KFDIA
  PCAPE(JL) = MAXVAL(ZCAPE(JL,:))
ENDDO

IF (LHOOK) CALL DR_HOOK('CUBASEN',1,ZHOOK_HANDLE)
END SUBROUTINE CUBASEN
!==============================================================================

SUBROUTINE CUBASMCN &
 & (KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
 & KK,&
 & PTEN,     PQEN,     PQSEN,    PUEN,     PVEN,&
 & PVERVEL,  PGEO,     PGEOH,    LDCUM,    KTYPE,    KLAB,&
 & KCBOT,    PMFU,     PMFUB,    PENTR,    PLRAIN,&
 & PTU,      PQU,      PLU,      PUU,      PVU,&
 & PMFUS,    PMFUQ,    PMFUL,    PDMFUP,   PMFUU,    PMFUV)  

!          M.TIEDTKE         E.C.M.W.F.     12/89

!          PURPOSE.
!          --------
!          THIS ROUTINE CALCULATES CLOUD BASE VALUES
!          FOR MIDLEVEL CONVECTION

!          INTERFACE
!          ---------

!          THIS ROUTINE IS CALLED FROM *CUASC*.
!          INPUT ARE ENVIRONMENTAL VALUES T,Q ETC
!          IT RETURNS CLOUDBASE VALUES FOR MIDLEVEL CONVECTION

!          METHOD.
!          --------
!          S. TIEDTKE (1989)

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KTDIA*        START OF THE VERTICAL LOOP
!    *KLEV*         NUMBER OF LEVELS
!    *KK*           ACTUAL LEVEL

!    INPUT PARAMETERS (REAL):

!    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
!    *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1)  KG/KG
!    *PQSEN*        ENVIRONMENT SPEC. SATURATION HUMIDITY (T+1)   KG/KG
!    *PUEN*         PROVISIONAL ENVIRONMENT U-VELOCITY (T+1)       M/S
!    *PVEN*         PROVISIONAL ENVIRONMENT V-VELOCITY (T+1)       M/S
!    *PVERVEL*      VERTICAL VELOCITY                             PA/S
!    *PGEO*         GEOPOTENTIAL                                  M2/S2
!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2
!    *PLRAIN*       RAIN WATER CONTENT IN UPDRAFTS                KG/KG

!    INPUT PARAMETERS (LOGICAL):

!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS 

!    UPDATED PARAMETERS (INTEGER):

!    *KTYPE*        TYPE OF CONVECTION
!                       1 = PENETRATIVE CONVECTION
!                       2 = SHALLOW CONVECTION
!                       3 = MIDLEVEL CONVECTION
!    *KLAB*         FLAG KLAB=1 FOR SUBCLOUD LEVELS
!                        KLAB=2 FOR CLOUD LEVELS
!    *KCBOT*        CLOUD BASE LEVEL

!    OUTPUT PARAMETERS (REAL):

!    *PMFU*         MASSFLUX IN UPDRAFTS                          KG/(M2*S)
!    *PMFUB*        MASSFLUX IN UPDRAFTS AT CLOUD BASE            KG/(M2*S)
!    *PENTR*        FRACTIONAL MASS ENTRAINMENT RATE               1/M
!    *PTU*          TEMPERATURE IN UPDRAFTS                         K
!    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                    KG/KG
!    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS              KG/KG
!    *PUU*          U-VELOCITY IN UPDRAFTS                         M/S
!    *PVU*          V-VELOCITY IN UPDRAFTS                         M/S
!    *PMFUS*        FLUX OF DRY STATIC ENERGY IN UPDRAFTS          J/(M2*S)
!    *PMFUQ*        FLUX OF SPEC. HUMIDITY IN UPDRAFTS            KG/(M2*S)
!    *PMFUL*        FLUX OF LIQUID WATER IN UPDRAFTS              KG/(M2*S)
!    *PDMFUP*       FLUX DIFFERENCE OF PRECIP. IN UPDRAFTS        KG/(M2*S)
!    *PMFUU*        FLUX OF U-MOMENTUM IN UPDRAFTS          (M/S)*KG/(M2*S)
!    *PMFUV*        FLUX OF V-MOMENTUM IN UPDRAFTS          (M/S)*KG/(M2*S)

!          EXTERNALS
!          ---------
!          NONE

!----------------------------------------------------------------------

USE MESSY_CONVECT_ECMWF_PARAM   , ONLY : RG, RCPD,                       &
                                         ENTRMID, RMFCMAX, RMFCMIN,      & !from YOECUMF
                                         LMFMID, LMFDUDV,                & !from YOECUMF
                                         LHOOK, DR_HOOK                    !from YOMHOOK

IMPLICIT NONE

INTEGER,INTENT(IN)    :: KLON 
INTEGER,INTENT(IN)    :: KLEV 
INTEGER,INTENT(IN)    :: KIDIA 
INTEGER,INTENT(IN)    :: KFDIA 
INTEGER               :: KTDIA ! Argument NOT used
INTEGER,INTENT(IN)    :: KK 
REAL(dp)   ,INTENT(IN)    :: PTEN(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PQEN(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PQSEN(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PUEN(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PVEN(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PVERVEL(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PGEO(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PGEOH(KLON,KLEV+1) 
LOGICAL           ,INTENT(IN)    :: LDCUM(KLON) 
INTEGER,INTENT(INOUT)   :: KTYPE(KLON) 
INTEGER,INTENT(INOUT)   :: KLAB(KLON,KLEV) 
INTEGER,INTENT(INOUT)   :: KCBOT(KLON) 
REAL(dp)   ,INTENT(INOUT)   :: PMFU(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT)   :: PMFUB(KLON) 
REAL(dp)   ,INTENT(INOUT)   :: PENTR(KLON) 
REAL(dp)   ,INTENT(INOUT)   :: PLRAIN(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT)   :: PTU(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT)   :: PQU(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT)   :: PLU(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT)   :: PUU(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT)   :: PVU(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT)   :: PMFUS(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT)   :: PMFUQ(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT)   :: PMFUL(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT)   :: PDMFUP(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT)   :: PMFUU(KLON) 
REAL(dp)   ,INTENT(INOUT)   :: PMFUV(KLON) 
INTEGER :: JL

REAL(dp) :: ZZZMB
REAL(dp) :: ZHOOK_HANDLE

INTRINSIC :: MAX, MIN
!----------------------------------------------------------------------

!*    1.           CALCULATE ENTRAINMENT AND DETRAINMENT RATES
!                  -------------------------------------------

!DIR$ IVDEP
!OCL NOVREC
IF (LHOOK) CALL DR_HOOK('CUBASMCN',0,ZHOOK_HANDLE)
DO JL=KIDIA,KFDIA
  IF(LMFMID.AND.PGEO(JL,KK) > 5000.0_dp.AND.PGEO(JL,KK)<1.E5_dp) THEN
    IF(.NOT.LDCUM(JL).AND.KLAB(JL,KK+1) == 0 &
       & .AND.PQEN(JL,KK) > 0.80_dp*PQSEN(JL,KK)) THEN  
      PTU(JL,KK+1)=(RCPD*PTEN(JL,KK)+PGEO(JL,KK)-PGEOH(JL,KK+1))/RCPD
      PQU(JL,KK+1)=PQEN(JL,KK)
      PLU(JL,KK+1)=0.0_dp
      ZZZMB=MAX(RMFCMIN,-PVERVEL(JL,KK)/RG)
      ZZZMB=MIN(ZZZMB,RMFCMAX)
      PMFUB(JL)=ZZZMB
      PMFU(JL,KK+1)=PMFUB(JL)
      PMFUS(JL,KK+1)=PMFUB(JL)*(RCPD*PTU(JL,KK+1)+PGEOH(JL,KK+1))
      PMFUQ(JL,KK+1)=PMFUB(JL)*PQU(JL,KK+1)
      PMFUL(JL,KK+1)=0.0_dp
      PDMFUP(JL,KK+1)=0.0_dp
      PLRAIN(JL,KK+1)=0.0_dp
      KCBOT(JL)=KK
      KLAB(JL,KK+1)=1
      KTYPE(JL)=3
      PENTR(JL)=ENTRMID
      IF(LMFDUDV) THEN
        PUU(JL,KK+1)=PUEN(JL,KK)
        PVU(JL,KK+1)=PVEN(JL,KK)
        PMFUU(JL)=PMFUB(JL)*PUU(JL,KK+1)
        PMFUV(JL)=PMFUB(JL)*PVU(JL,KK+1)
      ENDIF
    ENDIF
  ENDIF
ENDDO

IF (LHOOK) CALL DR_HOOK('CUBASMCN',1,ZHOOK_HANDLE)
END SUBROUTINE CUBASMCN

!==============================================================================

SUBROUTINE CUBIDIAG &
 & ( KIDIA, KFDIA, KLON, KLEV,&
 & KCTOP, LD_LCUMASK,&
 & PA,    PB,   PR,   PU )  

!          P. Bechtold         E.C.M.W.F.     07/03

!          PURPOSE.
!          --------
!          SOLVES BIDIAGONAL SYSTEM
!          FOR IMPLICIT SOLUTION OF ADVECTION EQUATION

!          INTERFACE
!          ---------

!          THIS ROUTINE IS CALLED FROM *CUDUDV* AND CUDTDQ.
!          IT RETURNS UPDATED VALUE OF QUANTITY

!          METHOD.
!          --------
!          NUMERICAL RECIPES (Cambridge Press)
!          DERIVED FROM TRIDIAGONAL ALGORIHM WITH C=0.
!          (ONLY ONE FORWARD SUBSTITUTION NECESSARY)
!          M  x  U  = R
!          ELEMENTS OF MATRIX M ARE STORED IN VECTORS A, B, C

!          (  B(kctop-1)    C(kctop-1)    0          0        )
!          (  A(kctop)      B(kctop)      C(kctop)   0        )
!          (  0             A(jk)         B(jk)      C(jk)    )
!          (  0             0             A(klev)    B(klev)  )

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----

!    INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEV*         NUMBER OF LEVELS
!    *KCTOP*        CLOUD TOP LEVELS

!    INPUT PARAMETERS (REAL):

!    *PA, PB*       VECTORS CONTAINING DIAGONAL ELEMENTS 
!    *PR*           RHS VECTOR CONTAINING "CONSTANTS"

!    OUTPUT PARAMETERS (REAL):

!    *PU*            SOLUTION VECTOR = UPDATED VALUE OF QUANTITY

!          EXTERNALS
!          ---------
!          NONE

!----------------------------------------------------------------------

  USE MESSY_CONVECT_ECMWF_PARAM,  ONLY: LHOOK, DR_HOOK                    !from YOMHOOK

IMPLICIT NONE

!     DUMMY INTEGER 
INTEGER,INTENT(IN)    :: KIDIA 
INTEGER,INTENT(IN)    :: KFDIA 
INTEGER,INTENT(IN)    :: KLON 
INTEGER,INTENT(IN)    :: KLEV 
INTEGER,INTENT(IN)    :: KCTOP(KLON) 
LOGICAL           ,INTENT(IN)    :: LD_LCUMASK(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PA(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PB(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PR(KLON,KLEV) 
REAL(dp)   ,INTENT(OUT)   :: PU(KLON,KLEV) 
!     DUMMY REALS
!     DUMMY LOGICALS
!     LOCALS
INTEGER :: JK, JL
REAL(dp) :: ZBET
REAL(dp) :: ZHOOK_HANDLE

!----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CUBIDIAG',0,ZHOOK_HANDLE)

PU(KIDIA:KFDIA,:)=0._dp

! Forward Substitution

DO JK = 2, KLEV
  DO JL = KIDIA,KFDIA
    IF ( LD_LCUMASK(JL,JK) ) THEN
      IF ( JK==KCTOP(JL)-1 ) THEN
        ZBET      =1.0_dp/(PB(JL,JK)+1.E-20_dp)
        PU(JL,JK) = PR(JL,JK) * ZBET
      ELSEIF ( JK>KCTOP(JL)-1 ) THEN
        ZBET      = 1.0_dp/(PB(JL,JK) + 1.E-20_dp)
        PU(JL,JK) =(PR(JL,JK)-PA(JL,JK)*PU(JL,JK-1))*ZBET
      ENDIF
    ENDIF
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('CUBIDIAG',1,ZHOOK_HANDLE)
END SUBROUTINE CUBIDIAG

!==============================================================================
SUBROUTINE CUBIDIAGAD &
 & ( KIDIA, KFDIA, KLON, KLEV,&
 & KCTOP, LD_LCUMASK,&
 & PA5,    PB5,   PR5,   PU5,&
 & PA ,    PB ,   PR ,   PU )  

!          P. Bechtold         E.C.M.W.F.     07/03

!          PURPOSE.
!          --------
!          SOLVES BIDIAGONAL SYSTEM
!          FOR IMPLICIT SOLUTION OF ADVECTION EQUATION
!          ADJOINT VERSION

!          INTERFACE
!          ---------

!          THIS ROUTINE IS CALLED FROM *CUDUDV* AND CUDTDQ.
!          IT RETURNS UPDATED VALUE OF QUANTITY

!          METHOD.
!          --------
!          NUMERICAL RECIPES (Cambridge Press)
!          DERIVED FROM TRIDIAGONAL ALGORIHM WITH C=0.
!          (ONLY ONE FORWARD SUBSTITUTION NECESSARY)
!          M  x  U  = R
!          ELEMENTS OF MATRIX M ARE STORED IN VECTORS A, B, C

!          (  B(kctop-1)    C(kctop-1)    0          0        )
!          (  A(kctop)      B(kctop)      C(kctop)   0        )
!          (  0             A(jk)         B(jk)      C(jk)    )
!          (  0             0             A(klev)    B(klev)  )

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----

!    INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEV*         NUMBER OF LEVELS
!    *KCTOP*        CLOUD TOP LEVELS

!    INPUT PARAMETERS (REAL):

!    Trajectory 
!    *PA5, PB5*     VECTORS CONTAINING DIAGONAL ELEMENTS 
!    *PR5*          RHS VECTOR CONTAINING "CONSTANTS"

!    Adjoint 
!    *PA, PB*       VECTORS CONTAINING DIAGONAL ELEMENTS 
!    *PR*           RHS VECTOR CONTAINING "CONSTANTS"

!    OUTPUT PARAMETERS (REAL):

!    Trajectory 
!    *PU5*           SOLUTION VECTOR = UPDATED VALUE OF QUANTITY

!    Adjoint
!    *PU*            SOLUTION VECTOR = UPDATED VALUE OF QUANTITY

!          EXTERNALS
!          ---------
!          NONE

!     Modifications
!     -------------
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!       P. LOPEZ, ECMWF, October 2003: Adjoint version.

!----------------------------------------------------------------------

  USE MESSY_CONVECT_ECMWF_PARAM,  ONLY: LHOOK, DR_HOOK                    !from YOMHOOK

IMPLICIT NONE

!     DUMMY INTEGER 
INTEGER,INTENT(IN)    :: KIDIA 
INTEGER,INTENT(IN)    :: KFDIA 
INTEGER,INTENT(IN)    :: KLON 
INTEGER,INTENT(IN)    :: KLEV 
INTEGER,INTENT(IN)    :: KCTOP(KLON) 
LOGICAL           ,INTENT(IN)    :: LD_LCUMASK(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PA5(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PB5(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PR5(KLON,KLEV) 
REAL(dp)   ,INTENT(OUT)   :: PU5(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PA(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PB(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PR(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PU(KLON,KLEV) 
!     DUMMY REALS
!     DUMMY LOGICALS
!     LOCALS
INTEGER :: JK, JL
REAL(dp) :: ZBET5(KLON,KLEV)
REAL(dp) :: ZBET
REAL(dp) :: ZHOOK_HANDLE

!----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CUBIDIAGAD',0,ZHOOK_HANDLE)

PU5(:,:)=0._dp

! Forward Substitution

DO JK = 2, KLEV
  DO JL = KIDIA,KFDIA
    IF ( LD_LCUMASK(JL,JK) ) THEN
      IF ( JK==KCTOP(JL)-1 ) THEN
        ZBET5(JL,JK) =1.0_dp/(PB5(JL,JK)+1.E-20_dp)
        PU5(JL,JK) = PR5(JL,JK) * ZBET5(JL,JK)
      ELSEIF ( JK>KCTOP(JL)-1 ) THEN
        ZBET5(JL,JK) = 1.0_dp/(PB5(JL,JK) + 1.E-20_dp)
        PU5(JL,JK) =(PR5(JL,JK)-PA5(JL,JK)*PU5(JL,JK-1))*ZBET5(JL,JK)
      ENDIF
    ENDIF
  ENDDO
ENDDO

! --------------------------------------------------------

!                 ADJOINT CALCULATIONS

! --------------------------------------------------------

! Forward Substitution

DO JK = KLEV,2,-1
  DO JL = KIDIA,KFDIA
    IF ( LD_LCUMASK(JL,JK) ) THEN
      IF ( JK==KCTOP(JL)-1 ) THEN
        ZBET=0.0_dp

        PR(JL,JK) = PR(JL,JK) + ZBET5(JL,JK)*PU(JL,JK)
        ZBET      = ZBET      + PR5(JL,JK)  *PU(JL,JK)
        PU(JL,JK) =0.0_dp

        PB(JL,JK) = PB(JL,JK) - ZBET5(JL,JK)*ZBET5(JL,JK)*ZBET
        ZBET =0.0_dp
   
      ELSEIF ( JK>KCTOP(JL)-1 ) THEN
        ZBET=0.0_dp

        PR(JL,JK)  = PR(JL,JK) + ZBET5(JL,JK)*PU(JL,JK)
        PA(JL,JK)  = PA(JL,JK) - PU5(JL,JK-1)*ZBET5(JL,JK)*PU(JL,JK)
        PU(JL,JK-1)= PU(JL,JK-1) - PA5(JL,JK)*ZBET5(JL,JK)*PU(JL,JK)
        ZBET = ZBET + (PR5(JL,JK)-PA5(JL,JK)*PU5(JL,JK-1))*PU(JL,JK)
        PU(JL,JK) =0.0_dp
   
        PB(JL,JK) = PB(JL,JK) - ZBET5(JL,JK)*ZBET5(JL,JK)*ZBET
        ZBET =0.0_dp
      ENDIF
    ENDIF
  ENDDO
ENDDO

PU(:,:)=0._dp

IF (LHOOK) CALL DR_HOOK('CUBIDIAGAD',1,ZHOOK_HANDLE)
END SUBROUTINE CUBIDIAGAD
!==============================================================================

SUBROUTINE CUBIDIAGTL &
 & ( KIDIA, KFDIA, KLON, KLEV,&
 & KCTOP, LD_LCUMASK,&
 & PA5,    PB5,   PR5,   PU5,&
 & PA ,    PB ,   PR ,   PU )  

!          P. Bechtold         E.C.M.W.F.     07/03

!          PURPOSE.
!          --------
!          SOLVES BIDIAGONAL SYSTEM
!          FOR IMPLICIT SOLUTION OF ADVECTION EQUATION
!          TANGENT LINEAR VERSION

!          INTERFACE
!          ---------

!          THIS ROUTINE IS CALLED FROM *CUDUDV* AND CUDTDQ.
!          IT RETURNS UPDATED VALUE OF QUANTITY

!          METHOD.
!          --------
!          NUMERICAL RECIPES (Cambridge Press)
!          DERIVED FROM TRIDIAGONAL ALGORIHM WITH C=0.
!          (ONLY ONE FORWARD SUBSTITUTION NECESSARY)
!          M  x  U  = R
!          ELEMENTS OF MATRIX M ARE STORED IN VECTORS A, B, C

!          (  B(kctop-1)    C(kctop-1)    0          0        )
!          (  A(kctop)      B(kctop)      C(kctop)   0        )
!          (  0             A(jk)         B(jk)      C(jk)    )
!          (  0             0             A(klev)    B(klev)  )

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----

!    INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEV*         NUMBER OF LEVELS
!    *KCTOP*        CLOUD TOP LEVELS

!    INPUT PARAMETERS (REAL):

!    Trajectory 
!    *PA5, PB5*     VECTORS CONTAINING DIAGONAL ELEMENTS 
!    *PR5*          RHS VECTOR CONTAINING "CONSTANTS"

!    Tangent-linear perturbations
!    *PA, PB*       VECTORS CONTAINING DIAGONAL ELEMENTS 
!    *PR*           RHS VECTOR CONTAINING "CONSTANTS"

!    OUTPUT PARAMETERS (REAL):

!    Trajectory 
!    *PU5*           SOLUTION VECTOR = UPDATED VALUE OF QUANTITY

!    Tangent-linear perturbations
!    *PU*            SOLUTION VECTOR = UPDATED VALUE OF QUANTITY

!          EXTERNALS
!          ---------
!          NONE

!     Modifications
!     -------------
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!       P. LOPEZ, ECMWF, October 2003: Tangent-linear version.

!----------------------------------------------------------------------

  USE MESSY_CONVECT_ECMWF_PARAM,  ONLY: LHOOK, DR_HOOK                    !from YOMHOOK

IMPLICIT NONE

!     DUMMY INTEGER 
INTEGER,INTENT(IN)    :: KIDIA 
INTEGER,INTENT(IN)    :: KFDIA 
INTEGER,INTENT(IN)    :: KLON 
INTEGER,INTENT(IN)    :: KLEV 
INTEGER,INTENT(IN)    :: KCTOP(KLON) 
LOGICAL           ,INTENT(IN)    :: LD_LCUMASK(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PA5(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PB5(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PR5(KLON,KLEV) 
REAL(dp)   ,INTENT(OUT)   :: PU5(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PA(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PB(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PR(KLON,KLEV) 
REAL(dp)   ,INTENT(OUT)   :: PU(KLON,KLEV) 
!     DUMMY REALS
!     DUMMY LOGICALS
!     LOCALS
INTEGER :: JK, JL
REAL(dp) :: ZBET5
REAL(dp) :: ZBET
REAL(dp) :: ZHOOK_HANDLE

!----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CUBIDIAGTL',0,ZHOOK_HANDLE)

PU5(:,:)=0._dp
PU (:,:)=0._dp

! Forward Substitution

DO JK = 2, KLEV
  DO JL = KIDIA,KFDIA
    IF ( LD_LCUMASK(JL,JK) ) THEN
      IF ( JK==KCTOP(JL)-1 ) THEN
        ZBET5     =1.0_dp/(PB5(JL,JK)+1.E-20_dp)
        ZBET      =-PB(JL,JK)*ZBET5*ZBET5
        PU5(JL,JK) = PR5(JL,JK) * ZBET5
        PU (JL,JK) = PR (JL,JK) * ZBET5 + PR5(JL,JK) * ZBET
      ELSEIF ( JK>KCTOP(JL)-1 ) THEN
        ZBET5     = 1.0_dp/(PB5(JL,JK) + 1.E-20_dp)
        ZBET      =-PB(JL,JK)*ZBET5*ZBET5
        PU5(JL,JK) =(PR5(JL,JK)-PA5(JL,JK)*PU5(JL,JK-1))*ZBET5
        PU (JL,JK) =(PR (JL,JK) &
         & - PA (JL,JK)*PU5(JL,JK-1)-PA5(JL,JK)*PU (JL,JK-1))*ZBET5 &
         & + (PR5(JL,JK)-PA5(JL,JK)*PU5(JL,JK-1))*ZBET  
      ENDIF
    ENDIF
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('CUBIDIAGTL',1,ZHOOK_HANDLE)
END SUBROUTINE CUBIDIAGTL

!==============================================================================

SUBROUTINE CUCTRACER &
 & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,  KTRAC,&
 & KTYPE,    KCTOP,    KCBOT,    KDPL,     KDTOP,&
 & LDCUM,    LDDRAF,   PTSPHY,   PAPH,&
 & PMFU,     PMFD,     PUDRATE,  PDDRATE,&
 & PCEN,     PTENC  )  

!**** *CUCTRACER* - COMPUTE CONVECTIVE TRANSPORT OF CHEM. TRACERS
!                   IMPORTANT: ROUTINE IS FOR POSITIVE DEFINIT QUANTITIES

!          P.BECHTOLD        E.C.M.W.F.              11/02/2004

!**   INTERFACE.
!     ----------

!          *CUTRACER* IS CALLED FROM *CUMASTR*

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KTDIA*        START OF THE VERTICAL LOOP
!    *KLEV*         NUMBER OF LEVELS
!    *KTRAC*        NUMBER OF CHEMICAL TRACERS

!    *KTYPE*        CONVECTION TYPE (DEEP SHALLOW MID-LEVEL)
!    *KCTOP*        CLOUD TOP  LEVEL
!    *KCBOT*        CLOUD BASE LEVEL
!    *KDPL*         DEPARTURE LEVEL
!    *KDTOP*        DOWNDRAFT TOP LEVEL

!    INPUT PARAMETERS (LOGICAL):

!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS 
!    *LDDRAF*       FLAG: .TRUE. IF DOWNDRAFTS EXIST

!    INPUT PARAMETERS (REAL):

!    *PTSPHY*       PHYSICS TIME-STEP                              S
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS            PA
!    *PCEN*         PROVISIONAL ENVIRONMENT TRACER CONCENTRATION
!    *PMFU*         MASSFLUX UPDRAFTS                             KG/(M2*S)
!    *PMFD*         MASSFLUX DOWNDRAFTS                           KG/(M2*S)
!    *PUDRATE       UPDRAFT DETRAINMENT                           KG/(M2*S)
!    *PDDRATE       DOWNDRAFT DETRAINMENT                         KG/(M2*S)

!    UPDATED PARAMETERS (REAL):

!    *PTENC*        UPDATED TENDENCY OF CHEM. TRACERS              1/S

!          METHOD
!          -------
!     EXPLICIT UPSTREAM AND IMPLICIT SOLUTION OF VERTICAL ADVECTION
!     DEPENDING ON VALUE OF RMFSOLCT: 0=EXPLICIT 0-1 SEMI-IMPLICIT >=1 IMPLICIT

!     FOR EXPLICIT SOLUTION: ONLY ONE SINGLE ITERATION
!     FOR IMPLICIT SOLUTION: FIRST IMPLICIT SOLVER, THEN EXPLICIT SOLVER
!                            TO CORRECT TENDENCIES BELOW CLOUD BASE


!------------------------------------------------------------------------------------
!     COMMENTS FOR OFFLINE USERS IN CHEMICAL TRANSPORT MODELS
!     (i.e. reading mass fluxes and detrainment rates from ECMWF archive:
!      ------------------------------------------------------------------
!     KCTOP IS FIRST LEVEL FROM TOP WHERE PMFU>0      
!     KDTOP IS FIRST LEVEL FROM TOP WHERE PMFD<0      
!     KCBOT IS NOT NEEDED FOR EXPLICIT SOLUTION
!     ATTENTION: ON ARCHIVE DETRAINMENT RATES HAVE UNITS KG/(M3*S), SO FOR USE
!                IN CURRENT ROUTINE YOU HAVE TO MULTIPLY ARCHIVED VALUES BY DZ
!     LDCUM  IS TRUE IF CONVECTION EXISTS, i.e. IF PMFU>0 IN COLUMN OR IF
!                       KCTOP>0 AND KCTOP<KLEV
!     LDDRAF IS TRUE IF DOWNDRAUGHTS EXIST IF PMFD<0 IN COLUMN OR IF
!                       KDTOP>0 AND KDTOP<KLEV
!     IF MASSFLUX SATISFIES CFL CRITERIUM M<=DP/Dt IT IS SUFFICIENT TO 
!     ONLY CONSIDER EXPLICIT SOLUTION (RMFSOLCT=0), IN THIS CASE
!     YOU CAN IGNORE THE IMPLICIT PART 7.0 OF CURRENT ROUTINE
!------------------------------------------------------------------------------------

!          EXTERNALS
!          ---------
!          CUBIDIAG

!          MODIFICATIONS
!          -------------
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!----------------------------------------------------------------------

USE MESSY_CONVECT_ECMWF_PARAM   , ONLY : RG,                         &
                                         RMFSOLCT, RMFCMIN,          & !from YOECUMF
                                         LHOOK, DR_HOOK                !from YOMHOOK

IMPLICIT NONE

INTEGER,INTENT(IN)    :: KLON
INTEGER,INTENT(IN)    :: KLEV
INTEGER,INTENT(IN)    :: KTRAC 
INTEGER,INTENT(IN)    :: KIDIA
INTEGER,INTENT(IN)    :: KFDIA 
INTEGER               :: KTDIA ! Argument NOT used
INTEGER,INTENT(IN)    :: KTYPE(KLON) 
INTEGER,INTENT(IN)    :: KCTOP(KLON) 
INTEGER,INTENT(IN)    :: KCBOT(KLON) 
INTEGER,INTENT(IN)    :: KDPL(KLON) 
INTEGER,INTENT(IN)    :: KDTOP(KLON) 
LOGICAL           ,INTENT(IN)    :: LDCUM(KLON) 
LOGICAL           ,INTENT(IN)    :: LDDRAF(KLON) 
REAL(dp)   ,INTENT(IN)    :: PTSPHY 
REAL(dp)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) 
REAL(dp)   ,INTENT(IN)    :: PMFU(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PMFD(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PUDRATE(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PDDRATE(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PCEN(KLON,KLEV,KTRAC) 
REAL(dp)   ,INTENT(INOUT) :: PTENC(KLON,KLEV,KTRAC) 
INTEGER :: IK, IKB, JK, JL, JN, JITER, JIT 

REAL(dp) :: ZZP, ZMFA, ZIMP, ZERATE, ZPOSI, ZTSPHY

!     ALLOCATABLE ARAYS
REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: ZCEN, ZCU, ZCD, ZTENC, ZMFC
REAL(dp), DIMENSION(:,:),   ALLOCATABLE :: ZDP, ZB,  ZR1
LOGICAL, DIMENSION(:,:),  ALLOCATABLE :: LLCUMASK, LLCUMBAS
REAL(dp) :: ZHOOK_HANDLE

!----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CUCTRACER',0,ZHOOK_HANDLE)
ZIMP=1.0_dp-RMFSOLCT
ZTSPHY=1._dp/PTSPHY

ALLOCATE(ZCEN(KLON,KLEV,KTRAC)) !Half-level environmental values
ALLOCATE(ZCU(KLON,KLEV,KTRAC))  !Updraft values
ALLOCATE(ZCD(KLON,KLEV,KTRAC))  !Downdraft values
ALLOCATE(ZTENC(KLON,KLEV,KTRAC))!Tendency
ALLOCATE(ZMFC(KLON,KLEV,KTRAC)) !Fluxes
ALLOCATE(ZDP(KLON,KLEV))        !Pressure difference
ALLOCATE(LLCUMASK(KLON,KLEV))   !Mask for convection

! Initialize Cumulus mask + some setups

DO JK=2,KLEV
   DO JL=KIDIA,KFDIA
     LLCUMASK(JL,JK)=.FALSE.
     IF(LDCUM(JL)) THEN
       ZDP(JL,JK)=RG/(PAPH(JL,JK+1)-PAPH(JL,JK))
       IF(JK>=KCTOP(JL)-1) THEN
          LLCUMASK(JL,JK)=.TRUE.
       ENDIF
     ENDIF
   ENDDO
ENDDO
!----------------------------------------------------------------------

DO JN=1,KTRAC 

!*    1.0          DEFINE TRACERS AT HALF LEVELS
!                  -----------------------------

  DO JK=2,KLEV
    IK=JK-1
    DO JL=KIDIA,KFDIA
      ZCEN(JL,JK,JN)=PCEN(JL,JK,JN)
      ZCD(JL,JK,JN) =PCEN(JL,IK,JN)
      ZCU(JL,JK,JN) =PCEN(JL,IK,JN)
      ZMFC(JL,JK,JN)=0.0_dp
      ZTENC(JL,JK,JN)=0.0_dp
    ENDDO
  ENDDO
  DO JL=KIDIA,KFDIA
      ZCU(JL,KLEV,JN) =PCEN(JL,KLEV,JN)
  ENDDO

!*    2.0          COMPUTE UPDRAFT VALUES
!                  ----------------------

  DO JK=KLEV-1,3,-1
    IK=JK+1
    DO JL=KIDIA,KFDIA
      IF ( LLCUMASK(JL,JK) ) THEN
        ZERATE=PMFU(JL,JK)-PMFU(JL,IK)+PUDRATE(JL,JK)
        IF(ZERATE<0._dp) ZERATE=0._dp 
        ZMFA=1.0_dp/MAX(RMFCMIN,PMFU(JL,JK))
        IF (JK >=KCTOP(JL)-1 )  THEN
          ZCU(JL,JK,JN)=( PMFU(JL,IK)*ZCU(JL,IK,JN)+ZERATE*PCEN(JL,JK,JN) &
           & -PUDRATE(JL,JK)*ZCU(JL,IK,JN) )*ZMFA   
! if you have a source term dc/dt=dcdt write 
!             ZCU(JL,JK,JN)=( PMFU(JL,IK)*ZCU(JL,IK,JN)+ZERATE*PCEN(JL,JK,JN) &
!                           -PUDRATE(JL,JK)*ZCU(JL,IK,JN) )*ZMFA 
!                           +dcdt(jl,ik,jn)*ptsphy
        ENDIF
      ENDIF
    ENDDO
  ENDDO


!*    3.0          COMPUTE DOWNDRAFT VALUES
!                  ------------------------

  DO JK=3,KLEV
    IK=JK-1
    DO JL=KIDIA,KFDIA
        IF ( LDDRAF(JL).AND.JK==KDTOP(JL) ) THEN
         !Nota: in order to avoid final negative Tracer values at LFS the allowed value of ZCD
         !      depends on the jump in mass flux at the LFS
            !ZCD(JL,JK,JN)=0.5_dp*ZCU(JL,JK,JN)+0.5_dp*PCEN(JL,IK,JN)
             ZCD(JL,JK,JN)=0.1_dp*ZCU(JL,JK,JN)+0.9_dp*PCEN(JL,IK,JN)
        ELSEIF ( LDDRAF(JL).AND.JK>KDTOP(JL) ) THEN
             ZERATE=-PMFD(JL,JK)+PMFD(JL,IK)+PDDRATE(JL,JK)
             IF(ZERATE<0._dp) ZERATE=0._dp 
             ZMFA=1._dp/MIN(-RMFCMIN,PMFD(JL,JK))
             ZCD(JL,JK,JN)=( PMFD(JL,IK)*ZCD(JL,IK,JN)-ZERATE*PCEN(JL,IK,JN) &
                            +PDDRATE(JL,JK)*ZCD(JL,IK,JN) )*ZMFA 
! if you have a source term dc/dt=dcdt write 
!             ZCD(JL,JK,JN)=( PMFD(JL,IK)*ZCD(JL,IK,JN)-ZERATE*PCEN(JL,IK,JN) &
!                            +PDDRATE(JL,JK)*ZCD(JL,IK,JN) &
!                           +dcdt(jl,ik,jn)*ptsphy
        ENDIF
    ENDDO
  ENDDO

! In order to avoid negative Tracer at KLEV adjust ZCD
  JK=KLEV
  IK=JK-1
  DO JL=KIDIA,KFDIA
    IF (LDDRAF(JL)) THEN
     ZPOSI=-ZDP(JL,JK)*(PMFU(JL,JK)*ZCU(JL,JK,JN)+PMFD(JL,JK)*ZCD(JL,JK,JN)&
                       &-(PMFU(JL,JK)+PMFD(JL,JK))*PCEN(JL,IK,JN) )
     IF( PCEN(JL,JK,JN)+ZPOSI*PTSPHY<0.0_dp ) THEN
        ZMFA=1._dp/MIN(-RMFCMIN,PMFD(JL,JK))
        ZCD(JL,JK,JN)=( (PMFU(JL,JK)+PMFD(JL,JK))*PCEN(JL,IK,JN)-PMFU(JL,JK)*ZCU(JL,JK,JN)&
                    &+PCEN(JL,JK,JN)/(PTSPHY*ZDP(JL,JK)) )*ZMFA
     ENDIF
    ENDIF
  ENDDO

ENDDO

JITER=1
!IF(RMFSOLCT>0.0_dp) JITER=2
!----------------------------------------------------------------------
 
DO JIT=1,JITER

  IF(JIT==2) ZIMP=1.0_dp

  DO JN=1,KTRAC

!*    4.0          COMPUTE FLUXES
!                  --------------

    DO JK=2,KLEV
      IK=JK-1
      DO JL=KIDIA,KFDIA
        IF(LLCUMASK(JL,JK)) THEN
          ZMFA=PMFU(JL,JK)+PMFD(JL,JK)
          ZMFC(JL,JK,JN)=PMFU(JL,JK)*ZCU(JL,JK,JN)+PMFD(JL,JK)*ZCD(JL,JK,JN)&
           & -ZIMP*ZMFA*ZCEN(JL,IK,JN)   
    !     IF(JK > KCBOT(JL)) THEN
    !       IKB=KCBOT(JL)
    !       ZZP=((PAPH(JL,KLEV+1)-PAPH(JL,JK))/(PAPH(JL,KLEV+1)-PAPH(JL,IKB)))
    !       IF(KTYPE(JL) == 3) THEN
    !          ZZP=ZZP*ZZP
    !       ENDIF
    !       ZMFC(JL,JK,JN)=ZMFC(JL,IKB,JN)*ZZP
    !     ENDIF
        ENDIF
      ENDDO
    ENDDO

!*    5.0          COMPUTE TENDENCIES = RHS
!                  ------------------------

    DO JK=2,KLEV-1
      IK=JK+1
      DO JL=KIDIA,KFDIA
        IF(LLCUMASK(JL,JK)) THEN
          ZTENC(JL,JK,JN)=ZDP(JL,JK)*(ZMFC(JL,IK,JN)-ZMFC(JL,JK,JN))
        ENDIF
      ENDDO
    ENDDO

    JK=KLEV
       DO JL=KIDIA,KFDIA
         IF(LDCUM(JL)) THEN
            ZTENC(JL,JK,JN)=-ZDP(JL,JK)*ZMFC(JL,JK,JN)
         ENDIF
       ENDDO

  ENDDO

  IF ( ZIMP==1.0_dp ) THEN


!*    6.0          UPDATE TENDENCIES
!                  -----------------

    DO JN=1,KTRAC
      DO JK=2,KLEV
        DO JL=KIDIA,KFDIA
          IF(LLCUMASK(JL,JK)) THEN
            PTENC(JL,JK,JN)=PTENC(JL,JK,JN)+ZTENC(JL,JK,JN)
          ENDIF
        ENDDO
      ENDDO
    ENDDO

  ELSE

!---------------------------------------------------------------------------

!*    7.0          IMPLICIT SOLUTION
!                  -----------------
 
   ! Fill bi-diagonal Matrix vectors A=k-1, B=k;
   ! reuse ZMFC=A and ZB=B;
   ! ZTENC corresponds to the RHS ("constants") of the equation
   ! The solution is in ZR1 
 
    ALLOCATE(ZB(KLON,KLEV))
    ALLOCATE(ZR1(KLON,KLEV))   
    ALLOCATE(LLCUMBAS(KLON,KLEV))
    LLCUMBAS(:,:)=.FALSE.
    ZB(:,:)=1._dp
 
    DO JN=1,KTRAC

   ! Fill vectors A, B and RHS
 
      DO JK=2,KLEV
        IK=JK+1
        DO JL=KIDIA,KFDIA
        ! LLCUMBAS(JL,JK)=LLCUMASK(JL,JK).AND.JK<=KCBOT(JL)
          LLCUMBAS(JL,JK)=LLCUMASK(JL,JK)
          IF(LLCUMBAS(JL,JK)) THEN
            ZZP=RMFSOLCT*ZDP(JL,JK)*PTSPHY
            ZMFC(JL,JK,JN)=-ZZP*(PMFU(JL,JK)+PMFD(JL,JK))
            ZTENC(JL,JK,JN) = ZTENC(JL,JK,JN)*PTSPHY+PCEN(JL,JK,JN)
            IF(JK<KLEV) THEN
              ZB(JL,JK)=1.0_dp+ZZP*(PMFU(JL,IK)+PMFD(JL,IK))
            ELSE
              ZB(JL,JK)=1.0_dp
            ENDIF
          ENDIF
        ENDDO
      ENDDO
 
      CALL CUBIDIAG&
       & ( KIDIA, KFDIA, KLON, KLEV,&
       & KCTOP, LLCUMBAS,&
       & ZMFC(:,:,JN),  ZB,   ZTENC(:,:,JN),   ZR1 )  
 
  ! Compute tendencies
  ! If second iteration was used, the values below cloud base would be
  ! replaced by the explicit solution
 
      DO JK=2,KLEV
        DO JL=KIDIA,KFDIA
          IF(LLCUMBAS(JL,JK)) THEN
            PTENC(JL,JK,JN)=PTENC(JL,JK,JN)+(ZR1(JL,JK)-PCEN(JL,JK,JN))*ZTSPHY
            !ZCEN(JL,JK,JN)=ZR1(JL,JK)
          ENDIF
        ENDDO
      ENDDO

    ENDDO
 
    DEALLOCATE(LLCUMBAS)
    DEALLOCATE(ZB)
    DEALLOCATE(ZR1)

  ENDIF
!---------------------------------------------------------------------------

ENDDO

DEALLOCATE(LLCUMASK)  
DEALLOCATE(ZDP)
DEALLOCATE(ZMFC)
DEALLOCATE(ZTENC)
DEALLOCATE(ZCD)
DEALLOCATE(ZCU)
DEALLOCATE(ZCEN)

IF (LHOOK) CALL DR_HOOK('CUCTRACER',1,ZHOOK_HANDLE)
END SUBROUTINE CUCTRACER

!==============================================================================

SUBROUTINE CUCTRACERAD &
 & ( KIDIA,    KFDIA,    KLON,    KTDIA,    KLEV,  KTRAC,&
 & KTYPE,    KCTOP,    KCBOT,   KDPL,     KDTOP,&
 & LDCUM,    LDDRAF,   PTSPHY,&
 & PAPH5,    PMFU5,    PMFD5,&
 & PUDRATE5, PDDRATE5, PCEN5,&
 & PTENC5,&
 & PAPH,     PMFU,     PMFD,&
 & PUDRATE,  PDDRATE,  PCEN,&
 & PTENC  )  

!**** *CUCTRACERAD* - COMPUTE CONVECTIVE TRANSPORT OF CHEM. TRACERS
!                     IMPORTANT: ROUTINE IS FOR POSITIVE DEFINIT QUANTITIES
!                     ADJOINT VERSION

!          P.BECHTOLD        E.C.M.W.F.              11/02/2004

!**   INTERFACE.
!     ----------

!          *CUCTRACERAD* IS CALLED FROM *CUMASTRN*

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KTDIA*        START OF THE VERTICAL LOOP
!    *KLEV*         NUMBER OF LEVELS
!    *KTRAC*        NUMBER OF CHEMICAL TRACERS

!    *KTYPE*        CONVECTION TYPE (DEEP SHALLOW MID-LEVEL)
!    *KCTOP*        CLOUD TOP  LEVEL
!    *KCBOT*        CLOUD BASE LEVEL
!    *KDPL*         DEPARTURE LEVEL
!    *KDTOP*        DOWNDRAFT TOP LEVEL

!    INPUT PARAMETERS (LOGICAL):

!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS 
!    *LDDRAF*       FLAG: .TRUE. IF DOWNDRAFTS EXIST

!    INPUT PARAMETERS (REAL):

!    *PTSPHY*       PHYSICS TIME-STEP                              S
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS            PA

!    Trajectory arrays

!    *PCEN5*        PROVISIONAL ENVIRONMENT TRACER CONCENTRATION
!    *PMFU5*        MASSFLUX UPDRAFTS                             KG/(M2*S)
!    *PMFD5*        MASSFLUX DOWNDRAFTS                           KG/(M2*S)
!    *PUDRATE5      UPDRAFT DETRAINMENT                           KG/(M2*S)
!    *PDDRATE5      DOWNDRAFT DETRAINMENT                         KG/(M2*S)

!    Adjoint arrays

!    *PCEN*         PROVISIONAL ENVIRONMENT TRACER CONCENTRATION
!    *PMFU*         MASSFLUX UPDRAFTS                             KG/(M2*S)
!    *PMFD*         MASSFLUX DOWNDRAFTS                           KG/(M2*S)
!    *PUDRATE       UPDRAFT DETRAINMENT                           KG/(M2*S)
!    *PDDRATE       DOWNDRAFT DETRAINMENT                         KG/(M2*S)

!    UPDATED PARAMETERS (REAL):

!    Trajectory arrays

!    *PTENC5*       UPDATED TENDENCY OF CHEM. TRACERS              1/S

!    Adjoint arrays

!    *PTENC*        UPDATED TENDENCY OF CHEM. TRACERS              1/S

!          METHOD
!          -------
!     EXPLICIT CENTRED-UPSTREAM AND IMPLICIT SOLUTION OF VERTICAL ADVECTION
!     DEPENDING ON VALUE OF RMFSOLCT: 0=EXPLICIT 1=IMPLICIT 0-1 SEMI-IMPLICIT

!     FOR EXPLICIT SOLUTION: ONLY ONE SINGLE ITERATION
!     FOR IMPLICIT SOLUTION: FIRST IMPLICIT SOLVER, THEN EXPLICIT SOLVER
!                            TO CORRECT TENDENCIES BELOW CLOUD BASE

!          EXTERNALS
!          ---------
!          CUBIDIAG
!          CUBIDIAGAD

!          MODIFICATIONS
!          -------------
!          P.LOPEZ,   E.C.M.W.F.  11/02/2004 : Adjoint version
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!----------------------------------------------------------------------

USE MESSY_CONVECT_ECMWF_PARAM   , ONLY : RG,                   &
                                         RMFSOLCT, RMFCMIN,    &  !from YOECUMF
                                         LHOOK, DR_HOOK           !from YOMHOOK

IMPLICIT NONE

INTEGER,INTENT(IN)    :: KLON
INTEGER,INTENT(IN)    :: KLEV
INTEGER,INTENT(IN)    :: KTRAC 
INTEGER,INTENT(IN)    :: KIDIA
INTEGER,INTENT(IN)    :: KFDIA 
INTEGER               :: KTDIA ! Argument NOT used
INTEGER,INTENT(IN)    :: KTYPE(KLON) 
INTEGER,INTENT(IN)    :: KCTOP(KLON) 
INTEGER,INTENT(IN)    :: KCBOT(KLON) 
INTEGER,INTENT(IN)    :: KDPL(KLON) 
INTEGER,INTENT(IN)    :: KDTOP(KLON) 
LOGICAL           ,INTENT(IN)    :: LDCUM(KLON) 
LOGICAL           ,INTENT(IN)    :: LDDRAF(KLON) 
REAL(dp)   ,INTENT(IN)    :: PTSPHY 
REAL(dp)   ,INTENT(IN)    :: PAPH5(KLON,KLEV+1) 
REAL(dp)   ,INTENT(IN)    :: PMFU5(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PMFD5(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PUDRATE5(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PDDRATE5(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PCEN5(KLON,KLEV,KTRAC) 
REAL(dp)   ,INTENT(INOUT) :: PTENC5(KLON,KLEV,KTRAC) 
REAL(dp)   ,INTENT(INOUT) :: PAPH(KLON,KLEV+1) 
REAL(dp)   ,INTENT(INOUT) :: PMFU(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PMFD(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PUDRATE(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PDDRATE(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PCEN(KLON,KLEV,KTRAC) 
REAL(dp)   ,INTENT(INOUT) :: PTENC(KLON,KLEV,KTRAC) 
! Trajectory arrays

! Tangent-linear perturbation arrays

INTEGER :: IK, IKB, JK, JL, JN, JITER, JIT 

REAL(dp) :: ZMFB5, ZIMP, ZTSPHY
REAL(dp) :: ZZP , ZMFA , ZMFB , ZERATE , ZFACT, ZFACT1, ZPOSI5, ZMFA5

!     ALLOCATABLE ARAYS
REAL(dp), DIMENSION(:,:),   ALLOCATABLE :: ZMFC15, ZDP5, ZR15, ZZP15, ZZP25, ZZP35, ZCD15
REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: ZCEN5, ZCU5, ZCD5, ZTENC5
REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: ZMFA15, ZMFA25
REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: ZERATE15, ZERATE25, ZERATE35, ZERATE45
REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: ZMFA35, ZB5, ZMFC25
REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: ZMFC5, ZMFC45
REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: ZTENC15

REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: ZCEN, ZCU, ZCD, ZTENC, ZMFC
REAL(dp), DIMENSION(:,:),   ALLOCATABLE :: ZDP, ZB,  ZR1
LOGICAL, DIMENSION(:,:),  ALLOCATABLE :: LLCUMASK, LLCUMBAS
REAL(dp) :: ZHOOK_HANDLE

!----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CUCTRACERAD',0,ZHOOK_HANDLE)
ZIMP   = 1.0_dp - RMFSOLCT
ZTSPHY = 1.0_dp / PTSPHY

! Trajectory arrays

ALLOCATE(ZCEN5(KLON,KLEV,KTRAC))     !Half-level environmental values
ALLOCATE(ZCU5(KLON,KLEV,KTRAC))      !Updraft values
ALLOCATE(ZCD5(KLON,KLEV,KTRAC))      !Downdraft values
ALLOCATE(ZCD15(KLON,KTRAC))          !Downdraft values at KLEV
ALLOCATE(ZTENC5(KLON,KLEV,KTRAC))    !Tendency
ALLOCATE(ZMFA15(KLON,KLEV,KTRAC))    !Fluxes  
ALLOCATE(ZMFA25(KLON,KLEV,KTRAC))    !Fluxes  
ALLOCATE(ZMFA35(KLON,KLEV,KTRAC,2))  !Fluxes
ALLOCATE(ZERATE15(KLON,KLEV,KTRAC))  !Entrainment rates
ALLOCATE(ZERATE25(KLON,KLEV,KTRAC))  !Entrainment rates
ALLOCATE(ZERATE35(KLON,KLEV,KTRAC))  !Entrainment rates
ALLOCATE(ZERATE45(KLON,KLEV,KTRAC))  !Entrainment rates
ALLOCATE(ZMFC15(KLON,KLEV))          !Fluxes
ALLOCATE(ZMFC25(KLON,KLEV,KTRAC,2))  !Fluxes
ALLOCATE(ZMFC5(KLON,KLEV,KTRAC,2))   !Fluxes
ALLOCATE(ZMFC45(KLON,KLEV,KTRAC,2))  !Fluxes
ALLOCATE(ZDP5(KLON,KLEV))            !Pressure difference
ALLOCATE(ZZP35(KLON,KLEV))           !Pressure difference
ALLOCATE(LLCUMASK(KLON,KLEV))        !Mask for convection

IF (RMFSOLCT > 0._dp) THEN
  ALLOCATE(ZB5(KLON,KLEV,KTRAC,2))
  ALLOCATE(ZR15(KLON,KLEV))   
  ALLOCATE(ZTENC15(KLON,KLEV,KTRAC,2))
  ALLOCATE(LLCUMBAS(KLON,KLEV))
  LLCUMBAS = .FALSE.
  ZB5 = 1._dp
ENDIF


! Initialize Cumulus mask + some setups

DO JK=2,KLEV
   DO JL=KIDIA,KFDIA
     LLCUMASK(JL,JK)=.FALSE.
     IF(LDCUM(JL)) THEN
       ZDP5(JL,JK)=RG/(PAPH5(JL,JK+1)-PAPH5(JL,JK))
       IF(JK>=KCTOP(JL)-1) THEN
          LLCUMASK(JL,JK)=.TRUE.
       ENDIF
     ENDIF
   ENDDO 
ENDDO 
!----------------------------------------------------------------------

DO JN=1,KTRAC 

!*    1.0          DEFINE TRACERS AT HALF LEVELS
!                  -----------------------------

   DO JK=2,KLEV
      IK=JK-1
      DO JL=KIDIA,KFDIA
         ZCEN5(JL,JK,JN)=PCEN5(JL,JK,JN)
         ZCU5(JL,JK,JN) =PCEN5(JL,IK,JN)
         ZCD5(JL,JK,JN) =PCEN5(JL,IK,JN)
         ZTENC5(JL,JK,JN)=0._dp
      ENDDO
   ENDDO
   DO JL=KIDIA,KFDIA
         ZCU5(JL,KLEV,JN)=PCEN5(JL,KLEV,JN)
   ENDDO


!*    2.0          COMPUTE UPDRAFT VALUES
!                  ----------------------

    DO JK=KLEV-1,3,-1
       IK=JK+1
       DO JL=KIDIA,KFDIA
        IF ( LLCUMASK(JL,JK) ) THEN
          ZERATE15(JL,JK,JN)=PMFU5(JL,JK)-PMFU5(JL,IK)+PUDRATE5(JL,JK)
          ZERATE35(JL,JK,JN)=ZERATE15(JL,JK,JN)
          IF(ZERATE15(JL,JK,JN)<0._dp) ZERATE15(JL,JK,JN)=0._dp
          ZMFA15(JL,JK,JN)=1._dp/MAX(RMFCMIN,PMFU5(JL,JK))
          IF (JK >=KCTOP(JL)-1 )  THEN
               ZCU5(JL,JK,JN)=( PMFU5(JL,IK)*ZCU5(JL,IK,JN) &
                           & + ZERATE15(JL,JK,JN)*PCEN5(JL,JK,JN) &
                           & - PUDRATE5(JL,JK)*ZCU5(JL,IK,JN) )*ZMFA15(JL,JK,JN)
          ENDIF
        ENDIF
       ENDDO
    ENDDO


!*    3.0          COMPUTE DOWNDRAFT VALUES
!                  ------------------------

    DO JK=3,KLEV
      IK=JK-1
      DO JL=KIDIA,KFDIA
        IF ( LDDRAF(JL).AND.JK==KDTOP(JL)) THEN
        !Nota: in order to avoid final negative Tracer values at the LFS the allowed value of ZCD
        !      depends on the jump in mass flux at the LFS
        ! ZCD5(JL,JK,JN)=.5_dp*(ZCU5(JL,JK,JN)+PCEN5(JL,IK,JN))
          ZCD5(JL,JK,JN)=.1_dp*ZCU5(JL,JK,JN)+0.9_dp*PCEN5(JL,IK,JN)
        ELSEIF ( LDDRAF(JL).AND.JK>KDTOP(JL) ) THEN
          ZERATE25(JL,JK,JN)=-PMFD5(JL,JK)+PMFD5(JL,IK)+PDDRATE5(JL,JK)
          ZERATE45(JL,JK,JN)=ZERATE25(JL,JK,JN)
          IF(ZERATE25(JL,JK,JN)<0._dp) ZERATE25(JL,JK,JN)=0._dp !!Modif
          ZMFA25(JL,JK,JN)=1._dp/MIN(-RMFCMIN,PMFD5(JL,JK))
          ZCD5(JL,JK,JN)=( PMFD5(JL,IK)*ZCD5(JL,IK,JN) &
                      & - ZERATE25(JL,JK,JN)*PCEN5(JL,IK,JN) &
                      & + PDDRATE5(JL,JK)*ZCD5(JL,IK,JN) )*ZMFA25(JL,JK,JN)
        ENDIF
      ENDDO
    ENDDO
! In order to avoid negative Tracer at KLEV adjust ZCD
  JK=KLEV
  IK=JK-1
  DO JL=KIDIA,KFDIA
    IF( LDDRAF(JL) ) THEN
     ZPOSI5=-ZDP5(JL,JK)*(PMFU5(JL,JK)*ZCU5(JL,JK,JN)+PMFD5(JL,JK)*ZCD5(JL,JK,JN)&
         & -(PMFU5(JL,JK)+PMFD5(JL,JK))*PCEN5(JL,IK,JN))
     ZCD15(JL,JN)=ZCD5(JL,JK,JN)
     IF( PCEN5(JL,JK,JN)+ZPOSI5*PTSPHY<0.0_dp ) THEN
        ZCD5(JL,JK,JN)=( (PMFU5(JL,JK)+PMFD5(JL,JK))*PCEN5(JL,IK,JN)-PMFU5(JL,JK)*ZCU5(JL,JK,JN)&
                    &+PCEN5(JL,JK,JN)/(PTSPHY*ZDP5(JL,JK)) )*ZMFA25(JL,JK,JN)
     ENDIF
    ENDIF
  ENDDO


ENDDO

ZMFC5(:,:,:,:)=0._dp

JITER=1
!IF(RMFSOLCT>0._dp) JITER=2
!----------------------------------------------------------------------
 
DO JIT=1,JITER

  IF(JIT==2) ZIMP=1._dp
  

  DO JN=1,KTRAC

!*    3.0          COMPUTE FLUXES
!                  --------------

    DO JK=2,KLEV
      IK=JK-1
      DO JL=KIDIA,KFDIA
        IF(LLCUMASK(JL,JK)) THEN

          ZMFA35(JL,JK,JN,JIT)=PMFU5(JL,JK)+PMFD5(JL,JK)

          ZMFC15(JL,JK)=PMFU5(JL,JK)*ZCU5(JL,JK,JN) &
                &      +PMFD5(JL,JK)*ZCD5(JL,JK,JN) &
                &      -ZIMP*ZMFA35(JL,JK,JN,JIT)*ZCEN5(JL,IK,JN) 

   ! Linear decrease of Flux between Cloud base and surface
        ! IF(JK > KCBOT(JL)) THEN
        !    IKB=KCBOT(JL)
        !    ZZP15(JL,JK)=(PAPH5(JL,KLEV+1)-PAPH5(JL,JK))/(PAPH5(JL,KLEV+1)-PAPH5(JL,IKB))          
        !    IF(KTYPE(JL) == 3) THEN
        !      ZZP25(JL,JK)=ZZP15(JL,JK)*ZZP15(JL,JK)
        !    ELSE
        !      ZZP25(JL,JK)=ZZP15(JL,JK)
        !    ENDIF
        !    ZMFC15(JL,JK)=ZMFC5(JL,IKB,JN,JIT)*ZZP25(JL,JK)
        ! ENDIF
          ZMFC5(JL,JK,JN,JIT)=ZMFC15(JL,JK)
        ENDIF
      ENDDO
    ENDDO


!*    4.0          COMPUTE TENDENCIES = RHS
!                  ------------------------


    DO JK=2,KLEV-1
       IK=JK+1
       DO JL=KIDIA,KFDIA
         IF(LLCUMASK(JL,JK)) THEN
            ZTENC5(JL,JK,JN)=ZDP5(JL,JK)*(ZMFC5(JL,IK,JN,JIT)-ZMFC5(JL,JK,JN,JIT))
         ENDIF
       ENDDO
    ENDDO

    JK=KLEV
    DO JL=KIDIA,KFDIA
      IF(LDCUM(JL)) THEN
         ZTENC5(JL,JK,JN)=-ZDP5(JL,JK)*ZMFC5(JL,JK,JN,JIT)
      ENDIF
    ENDDO


  ENDDO


 IF ( ZIMP==1._dp ) THEN


!*    4.0          UPDATE TENDENCIES
!                  -----------------

   DO JN=1,KTRAC
     DO JK=2,KLEV
        DO JL=KIDIA,KFDIA
          IF(LLCUMASK(JL,JK)) THEN
             PTENC5(JL,JK,JN)=PTENC5(JL,JK,JN)+ZTENC5(JL,JK,JN)
          ENDIF
        ENDDO
     ENDDO
   ENDDO

 ELSE

!---------------------------------------------------------------------------

!*    5.0          IMPLICIT SOLUTION
!                  -----------------
 
   ! Fill bi-diagonal Matrix vectors A=k-1, B=k;
   ! reuse ZMFC=A and ZB=B=1;
   ! ZDUDT corresponds to the RHS ("constants") of the equation
   ! The solution is in ZR1 and ZR2
 
  DO JN=1,KTRAC

   ! Fill vectors A, B and RHS
 
    DO JK=2,KLEV
       IK=JK+1
       DO JL=KIDIA,KFDIA
       ! LLCUMBAS(JL,JK)=LLCUMASK(JL,JK).AND.JK<=KCBOT(JL)
         LLCUMBAS(JL,JK)=LLCUMASK(JL,JK)
         IF(LLCUMBAS(JL,JK)) THEN
           ZZP35(JL,JK)=RMFSOLCT*ZDP5(JL,JK)*PTSPHY
           ZMFC45(JL,JK,JN,JIT)=-ZZP35(JL,JK)*(PMFU5(JL,JK)+PMFD5(JL,JK))
           ZTENC5(JL,JK,JN) = ZTENC5(JL,JK,JN)*PTSPHY+PCEN5(JL,JK,JN)
           IF(JK<KLEV) THEN        
             ZB5(JL,JK,JN,JIT)=1._dp+ZZP35(JL,JK)*(PMFU5(JL,IK)+PMFD5(JL,IK))
           ELSE
             ZB5(JL,JK,JN,JIT)=1._dp         
           ENDIF             
         ENDIF
       ENDDO
    ENDDO
 
! Store trajectory tendency array for current iteration

    ZTENC15(:,:,JN,JIT)=ZTENC5(:,:,JN)

    CALL CUBIDIAG &
       &( KIDIA, KFDIA, KLON, KLEV &
       &, KCTOP, LLCUMBAS &
       &, ZMFC45(:,:,JN,JIT),  ZB5(:,:,JN,JIT),   ZTENC5(:,:,JN),   ZR15 )
 
   ! Compute tendencies
   ! If second iteration was used, the values below cloud base would be
   ! replaced by the explicit solution
 
    DO JK=2,KLEV
       DO JL=KIDIA,KFDIA
         IF(LLCUMBAS(JL,JK)) THEN
           PTENC5(JL,JK,JN)=PTENC5(JL,JK,JN)+(ZR15(JL,JK)-PCEN5(JL,JK,JN))*ZTSPHY
!          ZCEN5(JL,JK,JN)=ZR15(JL,JK)  
         ENDIF
       ENDDO
    ENDDO

  ENDDO
 
 ENDIF
!----------------------------------------------------------------------

ENDDO

!----------------------------------------------------------------------
!
!                     ADJOINT CALCULATIONS
!
!----------------------------------------------------------------------

! Allocate adjoint arrays

ALLOCATE(ZCEN(KLON,KLEV,KTRAC))  !Half-level environmental values
ALLOCATE(ZCU(KLON,KLEV,KTRAC))   !Updraft values
ALLOCATE(ZCD(KLON,KLEV,KTRAC))   !Downdraft values
ALLOCATE(ZTENC(KLON,KLEV,KTRAC)) !Tendency
ALLOCATE(ZMFC(KLON,KLEV,KTRAC))  !Fluxes
ALLOCATE(ZDP(KLON,KLEV))         !Pressure difference

! INITIALIZATION OF LOCAL ADJOINT VARIABLES
! -----------------------------------------

ZCEN(:,:,:)=0._dp
ZCU(:,:,:)=0._dp
ZCD(:,:,:)=0._dp
ZTENC(:,:,:)=0._dp
ZMFC(:,:,:)=0._dp
ZDP(:,:)=0._dp

DO JN=1,KTRAC
   DO JK=2,KLEV
      DO JL=KIDIA,KFDIA
         PTENC(JL,JK,JN)=PTENC(JL,JK,JN)+PCEN(JL,JK,JN)*PTSPHY
      ENDDO
   ENDDO
ENDDO

DO JIT=JITER,1,-1

  ZIMP=MAX(0._dp,1._dp-RMFSOLCT)
  IF(JIT==2) ZIMP=1._dp

  IF ( ZIMP==1._dp ) THEN
    
!*    4.0          UPDATE TENDENCIES
!                  -----------------

   DO JN=1,KTRAC
     DO JK=KLEV,2,-1
        DO JL=KIDIA,KFDIA
          IF(LLCUMASK(JL,JK)) THEN
             ZTENC(JL,JK,JN) = ZTENC(JL,JK,JN) + PTENC(JL,JK,JN)
          ENDIF
        ENDDO
     ENDDO
   ENDDO

  ELSE

!---------------------------------------------------------------------------

!*    5.0          IMPLICIT SOLUTION
!                  -----------------

   ! Fill bi-diagonal Matrix vectors A=k-1, B=k;
   ! reuse ZMFC=A and ZB=B;
   ! ZDUDT corresponds to the RHS ("constants") of the equation
   ! The solution is in ZR1 

    ALLOCATE(ZB(KLON,KLEV))
    ALLOCATE(ZR1(KLON,KLEV))   
 
   DO JN=1,KTRAC

     ZR1(:,:)=0._dp
     ZB(:,:)=0._dp

     ! Compute tendencies
     ! If second iteration was used, the values below cloud base would be
     ! replaced by the explicit solution 

     DO JK=KLEV,2,-1
       DO JL=KIDIA,KFDIA
         IF(LLCUMBAS(JL,JK)) THEN
           ZR1(JL,JK)=ZR1(JL,JK)+PTENC(JL,JK,JN)*ZTSPHY 
           PCEN(JL,JK,JN)=PCEN(JL,JK,JN)-PTENC(JL,JK,JN)*ZTSPHY 
         ENDIF
       ENDDO
    ENDDO
    
    CALL CUBIDIAGAD &
       &( KIDIA, KFDIA, KLON, KLEV &
       &, KCTOP, LLCUMBAS &
       &, ZMFC45(:,:,JN,JIT), ZB5(:,:,JN,JIT), ZTENC15(:,:,JN,JIT), ZR15 &
       &, ZMFC (:,:,JN), ZB , ZTENC (:,:,JN),  ZR1  )
       
    DO JK=KLEV,2,-1
       IK=JK+1
       DO JL=KIDIA,KFDIA
         IF(LLCUMBAS(JL,JK)) THEN
           IF(JK<KLEV) THEN      
             ZZP=0._dp

             PMFU(JL,IK)=PMFU(JL,IK)+ZZP35(JL,JK)*ZB(JL,JK)
             PMFD(JL,IK)=PMFD(JL,IK)+ZZP35(JL,JK)*ZB(JL,JK)
             ZZP=ZZP+(PMFU5(JL,IK)+PMFD5(JL,IK))*ZB(JL,JK)
             ZB(JL,JK)=0._dp
           ELSE
             ZB(JL,JK)=0._dp         
           ENDIF             

           PCEN (JL,JK,JN)=PCEN (JL,JK,JN)+ZTENC(JL,JK,JN)
           ZTENC(JL,JK,JN)=ZTENC(JL,JK,JN)*PTSPHY

           ZZP = ZZP - (PMFU5(JL,JK)+PMFD5(JL,JK))*ZMFC(JL,JK,JN)
           PMFU(JL,JK) = PMFU(JL,JK) - ZZP35(JL,JK)*ZMFC(JL,JK,JN) 
           PMFD(JL,JK) = PMFD(JL,JK) - ZZP35(JL,JK)*ZMFC(JL,JK,JN) 
           ZMFC(JL,JK,JN)=0._dp

           ZDP(JL,JK) = ZDP(JL,JK) + RMFSOLCT*PTSPHY*ZZP
           ZZP=0._dp 
         ENDIF
       ENDDO
    ENDDO

   ENDDO

   DEALLOCATE(ZB)
   DEALLOCATE(ZR1)

  ENDIF

!*    4.0          COMPUTE TENDENCIES = RHS
!                  ------------------------

  DO JN=1,KTRAC
  
    JK=KLEV
    DO JL=KIDIA,KFDIA
      IF(LDCUM(JL)) THEN
         ZMFC(JL,JK,JN)=ZMFC(JL,JK,JN)-ZDP5(JL,JK)*ZTENC(JL,JK,JN)
         ZDP(JL,JK)=ZDP(JL,JK)-ZMFC5(JL,JK,JN,JIT)*ZTENC(JL,JK,JN)
         ZTENC(JL,JK,JN)=0._dp
      ENDIF
    ENDDO

    DO JK=KLEV-1,2,-1
       IK=JK+1
       DO JL=KIDIA,KFDIA
         IF(LLCUMASK(JL,JK)) THEN
           ZMFC(JL,IK,JN)=ZMFC(JL,IK,JN)+ZDP5(JL,JK)*ZTENC(JL,JK,JN)
           ZMFC(JL,JK,JN)=ZMFC(JL,JK,JN)-ZDP5(JL,JK)*ZTENC(JL,JK,JN)
           ZDP(JL,JK)=ZDP(JL,JK)+(ZMFC5(JL,IK,JN,JIT)-ZMFC5(JL,JK,JN,JIT)) &
                   & *ZTENC(JL,JK,JN)
           ZTENC(JL,JK,JN)=0._dp 
         ENDIF
       ENDDO
    ENDDO
    
    
!*    3.0          COMPUTE FLUXES
!                  --------------

    DO JK=KLEV,2,-1
      IK=JK-1
      DO JL=KIDIA,KFDIA
        IF(LLCUMASK(JL,JK)) THEN

          ZMFA=0._dp

       ! Linear decrease of Flux between Cloud base and surface
       !  IF(JK > KCBOT(JL)) THEN
       !     IKB=KCBOT(JL)
       !
       !     ZZP=0._dp
       !
       !     ZZP=ZZP+ZMFC5(JL,IKB,JN,JIT)*ZMFC(JL,JK,JN)
       !     ZMFC(JL,IKB,JN)=ZMFC(JL,IKB,JN)+ZZP25(JL,JK)*ZMFC(JL,JK,JN)
       !     ZMFC(JL,JK,JN)=0._dp
       !
       !     IF(KTYPE(JL) == 3) THEN
       !       ZZP=2._dp*ZZP15(JL,JK)*ZZP
       !     ENDIF
       !
       !     ZFACT=1._dp/(PAPH5(JL,KLEV+1)-PAPH5(JL,IKB))
       !     PAPH(JL,KLEV+1) = PAPH(JL,KLEV+1) + (1._dp-ZZP15(JL,JK))*ZFACT*ZZP
       !     PAPH(JL,JK    ) = PAPH(JL,JK    ) - ZFACT*ZZP
       !     PAPH(JL,IKB   ) = PAPH(JL,IKB   ) + ZZP15(JL,JK)*ZFACT*ZZP
       !     ZZP=0._dp
       !  ENDIF

          PMFU(JL,JK)   = PMFU(JL,JK)   + ZCU5(JL,JK,JN)* ZMFC(JL,JK,JN)
          ZCU(JL,JK,JN) = ZCU(JL,JK,JN) + PMFU5(JL,JK)  * ZMFC(JL,JK,JN)
          PMFD(JL,JK)   = PMFD(JL,JK)   + ZCD5(JL,JK,JN)* ZMFC(JL,JK,JN)
          ZCD(JL,JK,JN) = ZCD(JL,JK,JN) + PMFD5(JL,JK)  * ZMFC(JL,JK,JN)
          ZCEN(JL,IK,JN)= ZCEN(JL,IK,JN)- ZMFA35(JL,JK,JN,JIT)*ZIMP*ZMFC(JL,JK,JN)
          ZMFA          = ZMFA          - ZCEN5(JL,IK,JN)*ZIMP*ZMFC(JL,JK,JN) 
          ZMFC(JL,JK,JN)=0._dp
     
          PMFU(JL,JK)=PMFU(JL,JK)+ZMFA
          PMFD(JL,JK)=PMFD(JL,JK)+ZMFA
          ZMFA=0._dp
        ENDIF
      ENDDO
    ENDDO

  ENDDO

! LOOP OVER ITERATIONS
ENDDO

ZMFC(:,:,:)=0._dp

!*    1.0          DEFINE TRACERS AT HALF LEVELS
!                  -----------------------------

DO JN=1,KTRAC


!*    3.0          COMPUTE DOWNDRAFT VALUES
!                  ------------------------

! In order to avoid negative Tracer at KLEV adjust ZCD
  JK=KLEV
  IK=JK-1
  DO JL=KIDIA,KFDIA
    IF( LDDRAF(JL) ) THEN
     ZPOSI5=-ZDP5(JL,JK)*(PMFU5(JL,JK)*ZCU5(JL,JK,JN)+PMFD5(JL,JK)*ZCD15(JL,JN)&
         & -(PMFU5(JL,JK)+PMFD5(JL,JK))*PCEN5(JL,IK,JN))
     IF( PCEN5(JL,JK,JN)+ZPOSI5*PTSPHY<0.0_dp ) THEN
        ZMFA5=1._dp/MIN(-RMFCMIN,PMFD5(JL,JK))
        IF (PMFD5(JL,JK) < -RMFCMIN) THEN
           ZMFA=-PMFD(JL,JK)/(PMFD5(JL,JK)*PMFD5(JL,JK))
        ELSE
           ZMFA=0._dp
        ENDIF
        ZFACT=(PMFU5(JL,JK)+PMFD5(JL,JK))*PCEN5(JL,IK,JN)-PMFU5(JL,JK)*ZCU5(JL,JK,JN)&
                    &+PCEN5(JL,JK,JN)/(PTSPHY*ZDP5(JL,JK)) 

        PMFD(JL,JK)=PMFD(JL,JK)-ZFACT/(PMFD5(JL,JK)*PMFD5(JL,JK))+ZCD(JL,JK,JN)
        PMFU(JL,JK)=PMFU(JL,JK)+ZMFA5*PCEN5(JL,IK,JN)*ZCD(JL,JK,JN)
        PMFD(JL,JK)=PMFD(JL,JK)+ZMFA5*PCEN5(JL,IK,JN)*ZCD(JL,JK,JN)
        PCEN(JL,IK,JN)=PCEN(JL,IK,JN)+ZMFA5*(PMFU5(JL,JK)+PMFD5(JL,JK))*ZCD(JL,JK,JN)
        PMFU(JL,JK)=PMFU(JL,JK)-ZMFA5*ZCU5(JL,JK,JN)*ZCD(JL,JK,JN)
        ZCU(JL,JK,JN)=ZCU(JL,JK,JN)-ZMFA5*PMFU5(JL,JK)*ZCD(JL,JK,JN)
        PCEN(JL,JK,JN)=PCEN(JL,JK,JN)+ZMFA5/(PTSPHY*ZDP5(JL,JK))*ZCD(JL,JK,JN)
        ZDP(JL,JK)=ZDP(JL,JK)-ZMFA5*PCEN5(JL,JK,JN)/(PTSPHY*ZDP5(JL,JK)*ZDP5(JL,JK))*ZCD(JL,JK,JN)
        ZCD(JL,JK,JN)=0._dp

     ENDIF
    ENDIF
  ENDDO

    DO JK=KLEV,3,-1
      IK=JK-1
      DO JL=KIDIA,KFDIA

        ZMFA=0._dp
        ZERATE=0._dp

        IF ( LDDRAF(JL).AND.JK==KDTOP(JL) ) THEN
        ! ZCU(JL,JK,JN)  = ZCU(JL,JK,JN) +.5_dp*ZCD(JL,JK,JN)
        ! ZCEN(JL,IK,JN) = ZCEN(JL,IK,JN)+.5_dp*ZCD(JL,JK,JN)
          ZCU(JL,JK,JN)  = ZCU(JL,JK,JN) +.1_dp*ZCD(JL,JK,JN)
          ZCEN(JL,IK,JN) = ZCEN(JL,IK,JN)+.9_dp*ZCD(JL,JK,JN)
          ZCD(JL,JK,JN)=0._dp

        ELSEIF ( LDDRAF(JL).AND.JK>KDTOP(JL) ) THEN
          ZFACT=ZMFA25(JL,JK,JN)*ZCD(JL,JK,JN)
          ZMFA = ZMFA + ( PMFD5(JL,IK)*ZCD5(JL,IK,JN) &
             &          - ZERATE25(JL,JK,JN)*PCEN5(JL,IK,JN) &
             &          + PDDRATE5(JL,JK)*ZCD5(JL,IK,JN) )*ZCD(JL,JK,JN)
          PMFD(JL,IK)   = PMFD(JL,IK)   + ZFACT * ZCD5(JL,IK,JN) 
          ZCD(JL,IK,JN) = ZCD(JL,IK,JN) + ZFACT * (PMFD5(JL,IK)+PDDRATE5(JL,JK))
          ZERATE        = ZERATE        - ZFACT * PCEN5(JL,IK,JN)
          PCEN(JL,IK,JN)= PCEN(JL,IK,JN)- ZFACT * ZERATE25(JL,JK,JN)        
          PDDRATE(JL,JK)= PDDRATE(JL,JK)+ ZFACT * ZCD5(JL,IK,JN)
          ZCD(JL,JK,JN)=0._dp
             
          IF (PMFD5(JL,JK) < -RMFCMIN) THEN
            PMFD(JL,JK)=PMFD(JL,JK)-ZMFA/(PMFD5(JL,JK)*PMFD5(JL,JK))
          ENDIF
          ZMFA=0._dp
             
          IF (ZERATE45(JL,JK,JN)<0._dp) ZERATE=0._dp    

          PMFD(JL,JK)=PMFD(JL,JK)-ZERATE
          PMFD(JL,IK)=PMFD(JL,IK)+ZERATE
          PDDRATE(JL,JK)=PDDRATE(JL,JK)+ZERATE
          ZERATE=0._dp
        ENDIF
      ENDDO
    ENDDO

  
!*    2.0          COMPUTE UPDRAFT VALUES
!                  ----------------------

    DO JK=3,KLEV-1
       IK=JK+1
       DO JL=KIDIA,KFDIA

        ZMFA=0._dp
        ZERATE=0._dp

        IF ( LLCUMASK(JL,JK) ) THEN

          IF ( JK>=KCTOP(JL)-1 )  THEN
               ZFACT=ZMFA15(JL,JK,JN)*ZCU(JL,JK,JN)
               ZMFA = ZMFA + ( PMFU5(JL,IK)*ZCU5(JL,IK,JN) &
                  &          + ZERATE15(JL,JK,JN)*PCEN5(JL,JK,JN) &
                  &          - PUDRATE5(JL,JK)*ZCU5(JL,IK,JN) )*ZCU(JL,JK,JN)
               PMFU(JL,IK)   = PMFU(JL,IK)   + ZFACT * ZCU5(JL,IK,JN) 
               ZCU(JL,IK,JN) = ZCU(JL,IK,JN) + ZFACT * (PMFU5(JL,IK)-PUDRATE5(JL,JK))
               ZERATE        = ZERATE        + ZFACT * PCEN5(JL,JK,JN)
               PCEN(JL,JK,JN)= PCEN(JL,JK,JN)+ ZFACT * ZERATE15(JL,JK,JN)
               PUDRATE(JL,JK)= PUDRATE(JL,JK)- ZFACT * ZCU5(JL,IK,JN)
               ZCU(JL,JK,JN)=0._dp
          ENDIF

          IF (PMFU5(JL,JK) > RMFCMIN) THEN
            PMFU(JL,JK)=PMFU(JL,JK)-ZMFA/(PMFU5(JL,JK)*PMFU5(JL,JK))
          ENDIF
          ZMFA=0._dp
          
          IF (ZERATE35(JL,JK,JN)<0._dp) ZERATE=0._dp

          PMFU(JL,JK)=PMFU(JL,JK)+ZERATE
          PMFU(JL,IK)=PMFU(JL,IK)-ZERATE
          PUDRATE(JL,JK)=PUDRATE(JL,JK)+ZERATE
          ZERATE=0._dp

        ENDIF
       ENDDO
    ENDDO
  
!*    1.0          DEFINE TRACERS AT HALF LEVELS
!                  -----------------------------

   DO JL=KIDIA,KFDIA
         PCEN(JL,KLEV,JN)=PCEN(JL,KLEV,JN)+ZCU(JL,KLEV,JN)
         ZCU(JL,KLEV,JN)=0._dp
   ENDDO
   DO JK=KLEV,2,-1
      IK=JK-1
      DO JL=KIDIA,KFDIA
         PCEN (JL,IK,JN)=PCEN (JL,IK,JN)+ZCD (JL,JK,JN)
         ZCD (JL,JK,JN)=0._dp
       
         PCEN(JL,IK,JN)=PCEN(JL,IK,JN)+ZCU(JL,JK,JN)
         ZCU(JL,JK,JN)=0._dp

         PCEN(JL,JK,JN)=PCEN(JL,JK,JN)+ZCEN(JL,JK,JN)
         ZCEN(JL,JK,JN)=0._dp
      ENDDO
   ENDDO

ENDDO
  
!----------------------------------------------------------------------

! Initialize Cumulus mask + some setups

DO JK=KLEV,2,-1
   DO JL=KIDIA,KFDIA
      IF(LDCUM(JL)) THEN
         ZFACT=ZDP(JL,JK)*ZDP5(JL,JK)/(PAPH5(JL,JK+1)-PAPH5(JL,JK))
         PAPH(JL,JK+1)=PAPH(JL,JK+1)-ZFACT
         PAPH(JL,JK  )=PAPH(JL,JK  )+ZFACT
         ZDP(JL,JK)=0._dp
      ENDIF
   ENDDO
ENDDO

DEALLOCATE(ZDP)
DEALLOCATE(ZMFC)
DEALLOCATE(ZTENC)
DEALLOCATE(ZCD)
DEALLOCATE(ZCU)
DEALLOCATE(ZCEN)

IF (RMFSOLCT > 0._dp) THEN
  DEALLOCATE(LLCUMBAS)
  DEALLOCATE(ZTENC15)
  DEALLOCATE(ZR15)
  DEALLOCATE(ZB5)
ENDIF

DEALLOCATE(LLCUMASK)  
DEALLOCATE(ZZP35)
DEALLOCATE(ZDP5)
DEALLOCATE(ZMFC45)
DEALLOCATE(ZMFC5)
DEALLOCATE(ZMFC25)
DEALLOCATE(ZMFC15)
DEALLOCATE(ZERATE45)
DEALLOCATE(ZERATE35)
DEALLOCATE(ZERATE25)
DEALLOCATE(ZERATE15)
DEALLOCATE(ZMFA35)
DEALLOCATE(ZMFA25)
DEALLOCATE(ZMFA15)
DEALLOCATE(ZTENC5)
DEALLOCATE(ZCD15)
DEALLOCATE(ZCD5)
DEALLOCATE(ZCU5)
DEALLOCATE(ZCEN5)

IF (LHOOK) CALL DR_HOOK('CUCTRACERAD',1,ZHOOK_HANDLE)
END SUBROUTINE CUCTRACERAD

!==============================================================================

SUBROUTINE CUCTRACERTL &
 & ( KIDIA,    KFDIA,    KLON,    KTDIA,    KLEV,  KTRAC,&
 & KTYPE,    KCTOP,    KCBOT,   KDPL,     KDTOP,&
 & LDCUM,    LDDRAF,   PTSPHY,&
 & PAPH5,    PMFU5,    PMFD5,&
 & PUDRATE5, PDDRATE5, PCEN5,&
 & PTENC5,&
 & PAPH,     PMFU,     PMFD,&
 & PUDRATE,  PDDRATE,  PCEN,&
 & PTENC  )  

!**** *CUCTRACERTL* - COMPUTE CONVECTIVE TRANSPORT OF CHEM. TRACERS
!                     IMPORTANT: ROUTINE IS FOR POSITIVE DEFINIT QUANTITIES
!                     TANGENT-LINEAR VERSION

!          P.BECHTOLD        E.C.M.W.F.              11/02/2004

!**   INTERFACE.
!     ----------

!          *CUTRACERTL* IS CALLED FROM *CUMASTRNTL*

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KTDIA*        START OF THE VERTICAL LOOP
!    *KLEV*         NUMBER OF LEVELS
!    *KTRAC*        NUMBER OF CHEMICAL TRACERS

!    *KTYPE*        CONVECTION TYPE (DEEP SHALLOW MID-LEVEL)
!    *KCTOP*        CLOUD TOP  LEVEL
!    *KCBOT*        CLOUD BASE LEVEL
!    *KDPL*         DEPARTURE LEVEL
!    *KDTOP*        DOWNDRAFT TOP LEVEL

!    INPUT PARAMETERS (LOGICAL):

!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS 
!    *LDDRAF*       FLAG: .TRUE. IF DOWNDRAFTS EXIST

!    INPUT PARAMETERS (REAL):

!    *PTSPHY*       PHYSICS TIME-STEP                              S
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS            PA

!    Trajectory arrays

!    *PCEN5*        PROVISIONAL ENVIRONMENT TRACER CONCENTRATION
!    *PMFU5*        MASSFLUX UPDRAFTS                             KG/(M2*S)
!    *PMFD5*        MASSFLUX DOWNDRAFTS                           KG/(M2*S)
!    *PUDRATE5      UPDRAFT DETRAINMENT                           KG/(M2*S)
!    *PDDRATE5      DOWNDRAFT DETRAINMENT                         KG/(M2*S)

!    TL perturbation arrays

!    *PCEN*         PROVISIONAL ENVIRONMENT TRACER CONCENTRATION
!    *PMFU*         MASSFLUX UPDRAFTS                             KG/(M2*S)
!    *PMFD*         MASSFLUX DOWNDRAFTS                           KG/(M2*S)
!    *PUDRATE       UPDRAFT DETRAINMENT                           KG/(M2*S)
!    *PDDRATE       DOWNDRAFT DETRAINMENT                         KG/(M2*S)

!    UPDATED PARAMETERS (REAL):

!    Trajectory arrays

!    *PTENC5*       UPDATED TENDENCY OF CHEM. TRACERS              1/S

!    TL perturbation arrays

!    *PTENC*        UPDATED TENDENCY OF CHEM. TRACERS              1/S

!          METHOD
!          -------
!     EXPLICIT CENTRED-UPSTREAM AND IMPLICIT SOLUTION OF VERTICAL ADVECTION
!     DEPENDING ON VALUE OF RMFSOLCT: 0=EXPLICIT 1=IMPLICIT 0-1 SEMI-IMPLICIT

!     FOR EXPLICIT SOLUTION: ONLY ONE SINGLE ITERATION
!     FOR IMPLICIT SOLUTION: FIRST IMPLICIT SOLVER, THEN EXPLICIT SOLVER
!                            TO CORRECT TENDENCIES BELOW CLOUD BASE

!          EXTERNALS
!          ---------
!          CUBIDIAGTL

!          MODIFICATIONS
!          -------------
!          P.LOPEZ,   E.C.M.W.F.  11/02/2004 : Tangent-linear version
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!----------------------------------------------------------------------

USE MESSY_CONVECT_ECMWF_PARAM   , ONLY : RG,                   &
                                         RMFSOLCT, RMFCMIN,    &  !from YOECUMF
                                         LHOOK, DR_HOOK           !from YOMHOOK

IMPLICIT NONE

INTEGER,INTENT(IN)    :: KLON
INTEGER,INTENT(IN)    :: KLEV
INTEGER,INTENT(IN)    :: KTRAC 
INTEGER,INTENT(IN)    :: KIDIA
INTEGER,INTENT(IN)    :: KFDIA 
INTEGER               :: KTDIA ! Argument NOT used
INTEGER,INTENT(IN)    :: KTYPE(KLON) 
INTEGER,INTENT(IN)    :: KCTOP(KLON) 
INTEGER,INTENT(IN)    :: KCBOT(KLON) 
INTEGER,INTENT(IN)    :: KDPL(KLON) 
INTEGER,INTENT(IN)    :: KDTOP(KLON) 
LOGICAL           ,INTENT(IN)    :: LDCUM(KLON) 
LOGICAL           ,INTENT(IN)    :: LDDRAF(KLON) 
REAL(dp)   ,INTENT(IN)    :: PTSPHY 
REAL(dp)   ,INTENT(IN)    :: PAPH5(KLON,KLEV+1) 
REAL(dp)   ,INTENT(IN)    :: PMFU5(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PMFD5(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PUDRATE5(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PDDRATE5(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PCEN5(KLON,KLEV,KTRAC) 
REAL(dp)   ,INTENT(INOUT) :: PTENC5(KLON,KLEV,KTRAC) 
REAL(dp)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) 
REAL(dp)   ,INTENT(IN)    :: PMFU(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PMFD(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PUDRATE(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PDDRATE(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PCEN(KLON,KLEV,KTRAC) 
REAL(dp)   ,INTENT(INOUT) :: PTENC(KLON,KLEV,KTRAC) 
! Trajectory arrays

! Tangent-linear perturbation arrays

INTEGER :: IK, IKB, JK, JL, JN, JITER, JIT 

REAL(dp) :: ZZP5, ZMFA5, ZMFB5, ZERATE5, ZPOSI5, ZIMP, ZTSPHY
REAL(dp) :: ZZP , ZMFA , ZMFB , ZERATE , ZPOSI, ZFACT1, ZFACT

!     ALLOCATABLE ARAYS
REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: ZCEN5, ZCU5, ZCD5, ZTENC5, ZMFC5
REAL(dp), DIMENSION(:,:),   ALLOCATABLE :: ZDP5, ZB5,  ZR15
REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: ZCEN, ZCU, ZCD, ZTENC, ZMFC
REAL(dp), DIMENSION(:,:),   ALLOCATABLE :: ZDP, ZB,  ZR1
LOGICAL, DIMENSION(:,:),  ALLOCATABLE :: LLCUMASK, LLCUMBAS
REAL(dp) :: ZHOOK_HANDLE

!----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CUCTRACERTL',0,ZHOOK_HANDLE)
ZIMP   = 1.0_dp - RMFSOLCT
ZTSPHY = 1.0_dp / PTSPHY
! Trajectory arrays

ALLOCATE(ZCEN5(KLON,KLEV,KTRAC)) !Half-level environmental values
ALLOCATE(ZCU5(KLON,KLEV,KTRAC))  !Updraft values
ALLOCATE(ZCD5(KLON,KLEV,KTRAC))  !Downdraft values
ALLOCATE(ZTENC5(KLON,KLEV,KTRAC))!Tendency
ALLOCATE(ZMFC5(KLON,KLEV,KTRAC)) !Fluxes
ALLOCATE(ZDP5(KLON,KLEV))        !Pressure difference
ALLOCATE(LLCUMASK(KLON,KLEV))    !Mask for convection

! Tangent-linear perturbation arrays

ALLOCATE(ZCEN(KLON,KLEV,KTRAC)) !Half-level environmental values
ALLOCATE(ZCU(KLON,KLEV,KTRAC))  !Updraft values
ALLOCATE(ZCD(KLON,KLEV,KTRAC))  !Downdraft values
ALLOCATE(ZTENC(KLON,KLEV,KTRAC))!Tendency
ALLOCATE(ZMFC(KLON,KLEV,KTRAC)) !Fluxes
ALLOCATE(ZDP(KLON,KLEV))        !Pressure difference

! Initialize Cumulus mask + some setups

DO JK=2,KLEV
   DO JL=KIDIA,KFDIA
      LLCUMASK(JL,JK)=.FALSE.
      IF( LDCUM(JL) ) THEN
          ZDP5(JL,JK)=RG/(PAPH5(JL,JK+1)-PAPH5(JL,JK))
          ZDP (JL,JK)=-(PAPH(JL,JK+1)-PAPH(JL,JK))*ZDP5(JL,JK) &
                   & / (PAPH5(JL,JK+1)-PAPH5(JL,JK))
          IF(JK>=KCTOP(JL)-1) THEN
             LLCUMASK(JL,JK)=.TRUE.
          ENDIF 
      ENDIF 
   ENDDO 
ENDDO 
!----------------------------------------------------------------------

DO JN=1,KTRAC 

!*    1.0          DEFINE TRACERS AT HALF LEVELS
!                  -----------------------------

   DO JK=2,KLEV
      IK=JK-1
      DO JL=KIDIA,KFDIA
         ZCEN5(JL,JK,JN)=PCEN5(JL,JK,JN)
         ZCEN (JL,JK,JN)=PCEN (JL,JK,JN)
         ZCU5(JL,JK,JN) =PCEN5(JL,IK,JN)
         ZCU (JL,JK,JN) =PCEN (JL,IK,JN)
         ZCD5(JL,JK,JN) =PCEN5(JL,IK,JN)
         ZTENC5(JL,JK,JN)=0._dp
         ZTENC(JL,JK,JN) =0._dp
      ENDDO
   ENDDO
   DO JL=KIDIA,KFDIA
         ZCU5(JL,KLEV,JN) =PCEN5(JL,KLEV,JN)
         ZCU (JL,KLEV,JN) =PCEN (JL,KLEV,JN)
   ENDDO


!*    2.0          COMPUTE UPDRAFT VALUES
!                  ----------------------

    DO JK=KLEV-1,3,-1
       IK=MIN(KLEV,JK+1)
       DO JL=KIDIA,KFDIA
        IF ( LLCUMASK(JL,JK) ) THEN
          ZERATE5=PMFU5(JL,JK)-PMFU5(JL,IK)+PUDRATE5(JL,JK)
          ZERATE =PMFU (JL,JK)-PMFU (JL,IK)+PUDRATE (JL,JK)
          IF(ZERATE5<0._dp) THEN
            ZERATE5=0._dp
            ZERATE =0._dp
          ENDIF
          ZMFA5=1._dp/MAX(RMFCMIN,PMFU5(JL,JK))
          IF (PMFU5(JL,JK) > RMFCMIN) THEN
            ZMFA=-PMFU(JL,JK)/(PMFU5(JL,JK)*PMFU5(JL,JK))
          ELSE
            ZMFA=0._dp
          ENDIF

          IF ( JK>=KCTOP(JL)-1 )  THEN
             ZCU5(JL,JK,JN)=( PMFU5(JL,IK)*ZCU5(JL,IK,JN)+ZERATE5*PCEN5(JL,JK,JN) &
                   &       - PUDRATE5(JL,JK)*ZCU5(JL,IK,JN) )*ZMFA5

             ZCU(JL,JK,JN)= ZMFA*( PMFU5(JL,IK)*ZCU5(JL,IK,JN)+ZERATE5*PCEN5(JL,JK,JN) &
                   &      - PUDRATE5(JL,JK)*ZCU5(JL,IK,JN) ) + ZMFA5 & 
                   &      * ( PMFU(JL,IK)*ZCU5(JL,IK,JN) + PMFU5(JL,IK)*ZCU(JL,IK,JN) &
                   &      + ZERATE*PCEN5(JL,JK,JN) + ZERATE5*PCEN(JL,JK,JN) &
                   &      - PUDRATE(JL,JK)*ZCU5(JL,IK,JN) - PUDRATE5(JL,JK)*ZCU(JL,IK,JN)) 
          ENDIF
        ENDIF
       ENDDO
    ENDDO


!*    3.0          COMPUTE DOWNDRAFT VALUES
!                  ------------------------

    DO JK=3,KLEV
       IK=JK-1
       DO JL=KIDIA,KFDIA
         IF ( LDDRAF(JL).AND.JK==KDTOP(JL) ) THEN
        !Nota: in order to avoid final negative Tracer values at the LFS the allowed value of ZCD
        !      depends on the jump in mass flux at the LFS
          !ZCD5(JL,JK,JN)=.5_dp*(ZCU5(JL,JK,JN)+PCEN5(JL,IK,JN))
          !ZCD(JL,JK,JN) =.5_dp*(ZCU(JL,JK,JN)+ZCEN(JL,IK,JN))
           ZCD5(JL,JK,JN)=.1_dp*ZCU5(JL,JK,JN)+.9_dp*PCEN5(JL,IK,JN)
           ZCD(JL,JK,JN) =.1_dp*ZCU(JL,JK,JN)+.9_dp*PCEN(JL,IK,JN)
         ELSEIF ( LDDRAF(JL).AND.JK>KDTOP(JL) ) THEN
           ZERATE5=-PMFD5(JL,JK)+PMFD5(JL,IK)+PDDRATE5(JL,JK)  
           ZERATE=-PMFD(JL,JK)+PMFD(JL,IK)+PDDRATE(JL,JK)
           IF(ZERATE5<0._dp) THEN
             ZERATE5=0._dp
             ZERATE =0._dp
           ENDIF
           ZMFA5=1._dp/MIN(-RMFCMIN,PMFD5(JL,JK))
           IF (PMFD5(JL,JK) < -RMFCMIN) THEN
             ZMFA=-PMFD(JL,JK)/(PMFD5(JL,JK)*PMFD5(JL,JK))
           ELSE
             ZMFA=0._dp
           ENDIF
           ZCD5(JL,JK,JN)=( PMFD5(JL,IK)*ZCD5(JL,IK,JN)-ZERATE5*PCEN5(JL,IK,JN) &
                         + PDDRATE5(JL,JK)*ZCD5(JL,IK,JN) )*ZMFA5 
           ZCD (JL,JK,JN)= ZMFA*( PMFD5(JL,IK)*ZCD5(JL,IK,JN)-ZERATE5*PCEN5(JL,IK,JN) &
                  &      + PDDRATE5(JL,JK)*ZCD5(JL,IK,JN) ) + ZMFA5 & 
                  &      * ( PMFD (JL,IK)*ZCD5(JL,IK,JN) + PMFD5(JL,IK)*ZCD (JL,IK,JN) &
                  &      - ZERATE *PCEN5(JL,IK,JN) - ZERATE5*PCEN (JL,IK,JN) &
                  &      + PDDRATE (JL,JK)*ZCD5(JL,IK,JN) + PDDRATE5(JL,JK)*ZCD (JL,IK,JN))
         ENDIF
       ENDDO
    ENDDO

! In order to avoid negative Tracer at KLEV adjust ZCD
  JK=KLEV
  IK=JK-1
  DO JL=KIDIA,KFDIA
    IF( LDDRAF(JL) ) THEN
     ZPOSI5=-ZDP5(JL,JK)*(PMFU5(JL,JK)*ZCU5(JL,JK,JN)+PMFD5(JL,JK)*ZCD5(JL,JK,JN)&
                    & -(PMFU5(JL,JK)+PMFD5(JL,JK))*PCEN5(JL,IK,JN))
     IF( PCEN5(JL,JK,JN)+ZPOSI5*PTSPHY<0.0_dp ) THEN
        ZMFA5=1._dp/MIN(-RMFCMIN,PMFD5(JL,JK))
        IF (PMFD5(JL,JK) < -RMFCMIN) THEN
           ZMFA=-PMFD(JL,JK)/(PMFD5(JL,JK)*PMFD5(JL,JK))
        ELSE
           ZMFA=0._dp
        ENDIF
        ZFACT=(PMFU5(JL,JK)+PMFD5(JL,JK))*PCEN5(JL,IK,JN)-PMFU5(JL,JK)*ZCU5(JL,JK,JN)&
                    &+PCEN5(JL,JK,JN)/(PTSPHY*ZDP5(JL,JK)) 
        ZCD5(JL,JK,JN)=ZFACT*ZMFA5
        ZCD(JL,JK,JN)=ZMFA*ZFACT+ZMFA5*(&
          &(PMFU(JL,JK)+PMFD(JL,JK))*PCEN5(JL,IK,JN)+(PMFU5(JL,JK)+PMFD5(JL,JK))*PCEN(JL,IK,JN)&
          &-PMFU(JL,JK)*ZCU5(JL,JK,JN)-PMFU5(JL,JK)*ZCU(JL,JK,JN)&
          &+PCEN(JL,JK,JN)/(PTSPHY*ZDP5(JL,JK))-PCEN5(JL,JK,JN)*ZDP(JL,JK) &
          &/(PTSPHY*ZDP5(JL,JK)*ZDP5(JL,JK))&
          &)
     ENDIF
    ENDIF
  ENDDO

ENDDO              

ZMFC5(:,:,:)=0._dp


JITER=1
!IF(RMFSOLCT>0._dp) JITER=2
!----------------------------------------------------------------------
 
DO JIT=1,JITER

   IF(JIT==2) ZIMP=1._dp

  DO JN=1,KTRAC

!*    3.0          COMPUTE FLUXES
!                  --------------

    DO JK=2,KLEV
      IK=JK-1
      DO JL=KIDIA,KFDIA
        IF(LLCUMASK(JL,JK)) THEN
          ZMFA5=PMFU5(JL,JK)+PMFD5(JL,JK)
          ZMFA =PMFU (JL,JK)+PMFD (JL,JK)

          ZMFC5(JL,JK,JN)=PMFU5(JL,JK)*ZCU5(JL,JK,JN)+PMFD5(JL,JK)*ZCD5(JL,JK,JN)&
                       & -ZIMP*ZMFA5*ZCEN5(JL,IK,JN) 
          ZMFC(JL,JK,JN)=PMFU (JL,JK)*ZCU5(JL,JK,JN)+PMFU5(JL,JK)*ZCU (JL,JK,JN) &
                      & +PMFD (JL,JK)*ZCD5(JL,JK,JN)+PMFD5(JL,JK)*ZCD (JL,JK,JN) &
                      & -ZIMP*(ZMFA5*ZCEN(JL,IK,JN)+ZMFA*ZCEN5(JL,IK,JN)) 

   ! Linear decrease of Flux between Cloud base and surface
        ! IF(JK > KCBOT(JL)) THEN
        !    IKB=KCBOT(JL)
        !    ZFACT=1._dp/(PAPH5(JL,KLEV+1)-PAPH5(JL,IKB))
        !    ZZP5=(PAPH5(JL,KLEV+1)-PAPH5(JL,JK))*ZFACT
        !    ZZP =((PAPH (JL,KLEV+1)-PAPH (JL,JK)) &
        !      & -  ZZP5*(PAPH (JL,KLEV+1)-PAPH (JL,IKB)))*ZFACT 
        !    IF(KTYPE(JL) == 3) THEN
        !      ZZP =2._dp*ZZP*ZZP5
        !      ZZP5=ZZP5*ZZP5 
        !    ENDIF
        !    ZMFC5(JL,JK,JN)=ZMFC5(JL,IKB,JN)*ZZP5
        !    ZMFC (JL,JK,JN)=ZMFC5(JL,IKB,JN)*ZZP +ZMFC (JL,IKB,JN)*ZZP5
        ! ENDIF
        ENDIF
      ENDDO
    ENDDO


!*    4.0          COMPUTE TENDENCIES = RHS
!                  ------------------------

       
    DO JK=2,KLEV-1
       IK=JK+1
       DO JL=KIDIA,KFDIA
         IF(LLCUMASK(JL,JK)) THEN
            ZTENC5(JL,JK,JN)=ZDP5(JL,JK)*(ZMFC5(JL,IK,JN)-ZMFC5(JL,JK,JN))
            ZTENC (JL,JK,JN)=ZDP (JL,JK)*(ZMFC5(JL,IK,JN)-ZMFC5(JL,JK,JN)) &
                          & +ZDP5(JL,JK)*(ZMFC (JL,IK,JN)-ZMFC (JL,JK,JN))
         ENDIF
       ENDDO
    ENDDO
    
    JK=KLEV
       DO JL=KIDIA,KFDIA
         IF(LDCUM(JL)) THEN
            ZTENC5(JL,JK,JN)=-ZDP5(JL,JK)*ZMFC5(JL,JK,JN)
            ZTENC (JL,JK,JN)=-ZDP5(JL,JK)*ZMFC (JL,JK,JN) &
                           & -ZDP (JL,JK)*ZMFC5(JL,JK,JN)
         ENDIF
       ENDDO

  ENDDO


 IF ( ZIMP==1._dp ) THEN


!*    4.0          UPDATE TENDENCIES
!                  -----------------

   DO JN=1,KTRAC
     DO JK=2,KLEV
        DO JL=KIDIA,KFDIA
          IF(LLCUMASK(JL,JK)) THEN
             PTENC5(JL,JK,JN)=PTENC5(JL,JK,JN)+ZTENC5(JL,JK,JN)
             PTENC (JL,JK,JN)=PTENC (JL,JK,JN)+ZTENC (JL,JK,JN)
          ENDIF
        ENDDO
     ENDDO
   ENDDO

 ELSE

!---------------------------------------------------------------------------

!*    5.0          IMPLICIT SOLUTION
!                  -----------------
 
   ! Fill bi-diagonal Matrix vectors A=k-1, B=k;
   ! reuse ZMFC=A and ZB=B;
   ! ZTENC corresponds to the RHS ("constants") of the equation
   ! The solution is in ZR1
 
    ALLOCATE(ZB5(KLON,KLEV))
    ALLOCATE(ZR15(KLON,KLEV))   
    ALLOCATE(LLCUMBAS(KLON,KLEV))
    LLCUMBAS(:,:)=.FALSE.
    ZB5(:,:)=1.0_dp

    ALLOCATE(ZB(KLON,KLEV))
    ALLOCATE(ZR1(KLON,KLEV))   
    ZB(:,:)=0.0_dp
 
  DO JN=1,KTRAC

   ! Fill vectors A, B and RHS
 
    DO JK=2,KLEV
       IK=JK+1
       DO JL=KIDIA,KFDIA
       ! LLCUMBAS(JL,JK)=LLCUMASK(JL,JK).AND.JK<=KCBOT(JL)
         LLCUMBAS(JL,JK)=LLCUMASK(JL,JK)
         IF(LLCUMBAS(JL,JK)) THEN
           ZZP5=RMFSOLCT*ZDP5(JL,JK)*PTSPHY
           ZZP =RMFSOLCT*ZDP (JL,JK)*PTSPHY

           ZMFC5(JL,JK,JN)=-ZZP5*(PMFU5(JL,JK)+PMFD5(JL,JK))
           ZMFC (JL,JK,JN)=-ZZP *(PMFU5(JL,JK)+PMFD5(JL,JK)) &
                        &  -ZZP5*(PMFU (JL,JK)+PMFD (JL,JK)) 

           ZTENC5(JL,JK,JN) = ZTENC5(JL,JK,JN)*PTSPHY+PCEN5(JL,JK,JN)
           ZTENC (JL,JK,JN) = ZTENC (JL,JK,JN)*PTSPHY+PCEN (JL,JK,JN)

           IF(JK<KLEV) THEN
             ZB5(JL,JK)=1._dp+ZZP5*(PMFU5(JL,IK)+PMFD5(JL,IK))
             ZB(JL,JK)=ZZP5*(PMFU (JL,IK)+PMFD (JL,IK)) &
                    & +ZZP *(PMFU5(JL,IK)+PMFD5(JL,IK))
           ELSE
             ZB5(JL,JK)=1._dp
             ZB(JL,JK)=0._dp
           ENDIF
         ENDIF
       ENDDO
    ENDDO
 
    CALL CUBIDIAGTL&
       &( KIDIA, KFDIA, KLON, KLEV &
       &, KCTOP, LLCUMBAS &
       &, ZMFC5(:,:,JN),  ZB5,   ZTENC5(:,:,JN),   ZR15 &
       &, ZMFC (:,:,JN),  ZB ,   ZTENC (:,:,JN),   ZR1  )
 
    ! Compute tendencies
    ! If second iteration was used, the values below cloud base would be
    ! replaced by the explicit solution
 
    DO JK=2,KLEV
       DO JL=KIDIA,KFDIA
         IF(LLCUMBAS(JL,JK)) THEN
           PTENC5(JL,JK,JN)=PTENC5(JL,JK,JN)+(ZR15(JL,JK)-PCEN5(JL,JK,JN))*ZTSPHY
           PTENC (JL,JK,JN)=PTENC (JL,JK,JN)+(ZR1 (JL,JK)-PCEN (JL,JK,JN))*ZTSPHY
           ! ZCEN5(JL,JK,JN)=ZR15(JL,JK)
           ! ZCEN (JL,JK,JN)=ZR1 (JL,JK)
         ENDIF
       ENDDO
    ENDDO

  ENDDO
 
 
  DEALLOCATE(ZB5)
  DEALLOCATE(ZR15)
  DEALLOCATE(LLCUMBAS)
  DEALLOCATE(ZB)
  DEALLOCATE(ZR1)

 ENDIF
!---------------------------------------------------------------------------


ENDDO


DEALLOCATE(LLCUMASK)  

DEALLOCATE(ZDP5)
DEALLOCATE(ZMFC5)
DEALLOCATE(ZTENC5)
DEALLOCATE(ZCD5)
DEALLOCATE(ZCU5)
DEALLOCATE(ZCEN5)

DEALLOCATE(ZDP)
DEALLOCATE(ZMFC)
DEALLOCATE(ZTENC)
DEALLOCATE(ZCD)
DEALLOCATE(ZCU)
DEALLOCATE(ZCEN)

IF (LHOOK) CALL DR_HOOK('CUCTRACERTL',1,ZHOOK_HANDLE)
END SUBROUTINE CUCTRACERTL

!==============================================================================

SUBROUTINE CUDDRAFN &
 & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
 & LDDRAF,&
 & PTENH,    PQENH,    PUEN,     PVEN,&
 & PGEO,     PGEOH,    PAPH,     PRFL,&
 & PTD,      PQD,      PUD,      PVD,      PMFU,&
 & PMFD,     PMFDS,    PMFDQ,    PDMFDP,&
 & PDMFDE,   PMFDDE_RATE,  PMFDEN_RATE )  

!          THIS ROUTINE CALCULATES CUMULUS DOWNDRAFT DESCENT

!          M.TIEDTKE         E.C.M.W.F.    12/86 MODIF. 12/89

!          PURPOSE.
!          --------
!          TO PRODUCE THE VERTICAL PROFILES FOR CUMULUS DOWNDRAFTS
!          (I.E. T,Q,U AND V AND FLUXES)

!          INTERFACE
!          ---------

!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!          INPUT IS T,Q,P,PHI,U,V AT HALF LEVELS.
!          IT RETURNS FLUXES OF S,Q AND EVAPORATION RATE
!          AND U,V AT LEVELS WHERE DOWNDRAFT OCCURS

!          METHOD.
!          --------
!          CALCULATE MOIST DESCENT FOR ENTRAINING/DETRAINING PLUME BY
!          A) MOVING AIR DRY-ADIABATICALLY TO NEXT LEVEL BELOW AND
!          B) CORRECTING FOR EVAPORATION TO OBTAIN SATURATED STATE.

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KTDIA*        START OF THE VERTICAL LOOP
!    *KLEV*         NUMBER OF LEVELS

!    INPUT PARAMETERS (LOGICAL):

!    *LDDRAF*       .TRUE. IF DOWNDRAFTS EXIST

!    INPUT PARAMETERS (REAL):

!    *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS          K
!    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS     KG/KG
!    *PUEN*         PROVISIONAL ENVIRONMENT U-VELOCITY (T+1)      M/S
!    *PVEN*         PROVISIONAL ENVIRONMENT V-VELOCITY (T+1)      M/S
!    *PGEO*         GEOPOTENTIAL                                  M2/S2
!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                  M2/S2
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS           PA
!    *PMFU*         MASSFLUX UPDRAFTS                           KG/(M2*S)

!    UPDATED PARAMETERS (REAL):

!    *PRFL*         PRECIPITATION RATE                           KG/(M2*S)

!    OUTPUT PARAMETERS (REAL):

!    *PTD*          TEMPERATURE IN DOWNDRAFTS                      K
!    *PQD*          SPEC. HUMIDITY IN DOWNDRAFTS                 KG/KG
!    *PUD*          U-VELOCITY IN DOWNDRAFTS                      M/S
!    *PVD*          V-VELOCITY IN DOWNDRAFTS                      M/S
!    *PMFD*         MASSFLUX IN DOWNDRAFTS                       KG/(M2*S)
!    *PMFDS*        FLUX OF DRY STATIC ENERGY IN DOWNDRAFTS       J/(M2*S)
!    *PMFDQ*        FLUX OF SPEC. HUMIDITY IN DOWNDRAFTS         KG/(M2*S)
!    *PDMFDP*       FLUX DIFFERENCE OF PRECIP. IN DOWNDRAFTS     KG/(M2*S)
!    *PMFDDE_RATE*  DOWNDRAFT DETRAINMENT RATE                   KG/(M2*S)
!    *PMFDEN_RATE*  DOWNDRAFT ENTRAINMENT RATE                   KG/(M2*S)

!          EXTERNALS
!          ---------
!          *CUADJTQ* FOR ADJUSTING T AND Q DUE TO EVAPORATION IN
!          SATURATED DESCENT

!          REFERENCE
!          ---------
!          (TIEDTKE,1989)

!          MODIFICATIONS
!          -------------
!             92-09-21 : Update to Cy44      J.-J. MORCRETTE
!             03-08-28 : Clean-up detrainment rates   P. BECHTOLD
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!----------------------------------------------------------------------

USE MESSY_CONVECT_ECMWF_PARAM   , ONLY : RG, RD, RCPD, RETV,                 &
                                         ENTRDD, RMFCMIN, LMFDUDV, NJKT3,    &  !from YOECUMF
                                         LHOOK, DR_HOOK                         !from YOMHOOK

IMPLICIT NONE

INTEGER,INTENT(IN)    :: KLON 
INTEGER,INTENT(IN)    :: KLEV 
INTEGER,INTENT(IN)    :: KIDIA 
INTEGER,INTENT(IN)    :: KFDIA 
INTEGER               :: KTDIA ! Argument NOT used
LOGICAL           ,INTENT(IN)    :: LDDRAF(KLON) 
REAL(dp)   ,INTENT(IN)    :: PTENH(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PQENH(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PUEN(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PVEN(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PGEO(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PGEOH(KLON,KLEV+1) 
REAL(dp)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) 
REAL(dp)   ,INTENT(INOUT) :: PRFL(KLON) 
REAL(dp)   ,INTENT(INOUT) :: PTD(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PQD(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PUD(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PVD(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PMFU(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PMFD(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PMFDS(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PMFDQ(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PDMFDP(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PDMFDE(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PMFDDE_RATE(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PMFDEN_RATE(KLON,KLEV) 
REAL(dp) ::     ZDMFEN(KLON),           ZDMFDE(KLON),&
 & ZCOND(KLON),            ZOENTR(KLON),&
 & ZBUOY(KLON)  
REAL(dp) ::     ZPH(KLON)
LOGICAL ::  LLO2(KLON)
INTEGER :: ICALL, IK, IS, ITOPDE, JK, JL

REAL(dp) :: ZBUO, ZBUOYZ, ZDMFDP, ZDZ, ZENTR, ZMFDQK,&
 & ZMFDSK, ZMFDUK, ZMFDVK, ZQDDE, ZQEEN, ZRAIN, &
 & ZSDDE, ZSEEN, ZZENTR  
REAL(dp) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('CUDDRAFN',0,ZHOOK_HANDLE)
ITOPDE=NJKT3
!----------------------------------------------------------------------
!* change to operations
! some preliminar to compute top of downdraft detrainment layer
ITOPDE=KLEV-2
do jk=klev,2,-1
   if (paph(1,klev)-paph(1,jk) <60.E2 ) itopde=jk
enddo
!* change to operations
njkt3=itopde


!     1.           CALCULATE MOIST DESCENT FOR CUMULUS DOWNDRAFT BY
!                     (A) CALCULATING ENTRAINMENT/DETRAINMENT RATES, 
!                         INCLUDING ORGANIZED ENTRAINMENT DEPENDENT ON
!                         NEGATIVE BUOYANCY AND ASSUMING
!                         LINEAR DECREASE OF MASSFLUX IN PBL
!                     (B) DOING MOIST DESCENT - EVAPORATIVE COOLING
!                         AND MOISTENING IS CALCULATED IN *CUADJTQ*
!                     (C) CHECKING FOR NEGATIVE BUOYANCY AND
!                         SPECIFYING FINAL T,Q,U,V AND DOWNWARD FLUXES
!                    -------------------------------------------------

DO JL=KIDIA,KFDIA
  ZOENTR(JL)=0.0_dp
  ZBUOY(JL)=0.0_dp
  ZDMFEN(JL)=0.0_dp
  ZDMFDE(JL)=0.0_dp
  PDMFDE(JL,:)=0.0_dp
  PMFDDE_RATE(JL,:)=0.0_dp
  PMFDEN_RATE(JL,:)=0.0_dp
ENDDO

DO JK=3,KLEV
  IS=0
  DO JL=KIDIA,KFDIA
    ZPH(JL)=PAPH(JL,JK)
    LLO2(JL)=LDDRAF(JL).AND.PMFD(JL,JK-1) < 0.0_dp
    IF(LLO2(JL)) THEN
      IS=IS+1
    ENDIF
  ENDDO
  IF(IS == 0) CYCLE

  DO JL=KIDIA,KFDIA
    IF(LLO2(JL)) THEN
      ZENTR=ENTRDD*PMFD(JL,JK-1)*RD*PTENH(JL,JK-1)/&
       & (RG*PAPH(JL,JK-1))*(PAPH(JL,JK)-PAPH(JL,JK-1))  
      ZDMFEN(JL)=ZENTR
      ZDMFDE(JL)=ZENTR
    ENDIF
  ENDDO

! IF(KLEV == 60.OR.KLEV == 40) THEN
!   ITOPDE=KLEV-6
! ELSE
!   ITOPDE=KLEV-2
! ENDIF

  IF(JK > ITOPDE) THEN
    DO JL=KIDIA,KFDIA
      IF(LLO2(JL)) THEN
        ZDMFEN(JL)=0.0_dp
        ZDMFDE(JL)=PMFD(JL,ITOPDE)*&
         & (PAPH(JL,JK)-PAPH(JL,JK-1))/&
         & (PAPH(JL,KLEV+1)-PAPH(JL,ITOPDE))  
      ENDIF
    ENDDO
  ENDIF

  IF(JK <= ITOPDE) THEN
    DO JL=KIDIA,KFDIA
      IF(LLO2(JL)) THEN
        ZDZ=-(PGEOH(JL,JK-1)-PGEOH(JL,JK))/RG
        ZZENTR=ZOENTR(JL)*ZDZ*PMFD(JL,JK-1)
        ZDMFEN(JL)=ZDMFEN(JL)+ZZENTR
        ZDMFEN(JL)=MAX(ZDMFEN(JL),0.3_dp*PMFD(JL,JK-1))
        ZDMFEN(JL)=MAX(ZDMFEN(JL),-0.75_dp*PMFU(JL,JK)-&
         & (PMFD(JL,JK-1)-ZDMFDE(JL)))  
        ZDMFEN(JL)=MIN(ZDMFEN(JL),0.0_dp)
      ENDIF

      PDMFDE(JL,JK)=ZDMFEN(JL)-ZDMFDE(JL)

    ENDDO
  ENDIF
  DO JL=KIDIA,KFDIA
    IF(LLO2(JL)) THEN
      PMFD(JL,JK)=PMFD(JL,JK-1)+ZDMFEN(JL)-ZDMFDE(JL)
      ZSEEN=(RCPD*PTENH(JL,JK-1)+PGEOH(JL,JK-1))*ZDMFEN(JL)
      ZQEEN=PQENH(JL,JK-1)*ZDMFEN(JL)
      ZSDDE=(RCPD*PTD(JL,JK-1)+PGEOH(JL,JK-1))*ZDMFDE(JL)
      ZQDDE=PQD(JL,JK-1)*ZDMFDE(JL)
      ZMFDSK=PMFDS(JL,JK-1)+ZSEEN-ZSDDE
      ZMFDQK=PMFDQ(JL,JK-1)+ZQEEN-ZQDDE
      PQD(JL,JK)=ZMFDQK*(1.0_dp/MIN(-RMFCMIN,PMFD(JL,JK)))
      PTD(JL,JK)=(ZMFDSK*(1.0_dp/MIN(-RMFCMIN,PMFD(JL,JK)))-PGEOH(JL,JK))/RCPD
      PTD(JL,JK)=MIN(400._dp,PTD(JL,JK))
      PTD(JL,JK)=MAX(100._dp,PTD(JL,JK))
      ZCOND(JL)=PQD(JL,JK)
    ENDIF
  ENDDO

  IK=JK
  ICALL=2
  CALL CUADJTQ &
   & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
   & IK,&
   & ZPH,      PTD,      PQD,      LLO2,     ICALL )  

  DO JL=KIDIA,KFDIA
    IF(LLO2(JL)) THEN
      ZCOND(JL)=ZCOND(JL)-PQD(JL,JK)
      ZBUO=PTD(JL,JK)*(1.0_dp+RETV  *PQD(JL,JK))-&
       & PTENH(JL,JK)*(1.0_dp+RETV  *PQENH(JL,JK))  
      IF(PRFL(JL) > 0.0_dp.AND.PMFU(JL,JK) > 0.0_dp) THEN
        ZRAIN=PRFL(JL)/PMFU(JL,JK)
        ZBUO=ZBUO-PTD(JL,JK)*ZRAIN
      ENDIF
      IF(ZBUO >= 0.0_dp.OR.PRFL(JL) <= (PMFD(JL,JK)*ZCOND(JL))) THEN
        PMFD(JL,JK)=0.0_dp
        ZBUO=0.0_dp
      ENDIF
      PMFDS(JL,JK)=(RCPD*PTD(JL,JK)+PGEOH(JL,JK))*PMFD(JL,JK)
      PMFDQ(JL,JK)=PQD(JL,JK)*PMFD(JL,JK)
      ZDMFDP=-PMFD(JL,JK)*ZCOND(JL)
      PDMFDP(JL,JK-1)=ZDMFDP
      PRFL(JL)=PRFL(JL)+ZDMFDP

! COMPUTE ORGANIZED ENTRAINMENT FOR USE AT NEXT LEVEL

      ZBUOYZ=ZBUO/PTENH(JL,JK)
      ZBUOYZ=MIN(ZBUOYZ,0.0_dp)
      ZDZ=-(PGEO(JL,JK-1)-PGEO(JL,JK))
      ZBUOY(JL)=ZBUOY(JL)+ZBUOYZ*ZDZ
      ZOENTR(JL)=RG*ZBUOYZ*0.5_dp/(1.0_dp+ZBUOY(JL))

! STORE DOWNDRAUGHT DETRAINMENT RATES

      PMFDDE_RATE(JL,JK)=-ZDMFDE(JL)
      PMFDEN_RATE(JL,JK)=-ZDMFEN(JL)

    ENDIF
  ENDDO
  
  IF(LMFDUDV) THEN
    DO JL=KIDIA,KFDIA
      IF(LLO2(JL).AND.PMFD(JL,JK) < 0.0_dp) THEN
        ZMFDUK=PMFD(JL,JK-1)*PUD(JL,JK-1)+&
         & ZDMFEN(JL)*PUEN(JL,JK-1)-ZDMFDE(JL)*PUD(JL,JK-1)  
        ZMFDVK=PMFD(JL,JK-1)*PVD(JL,JK-1)+&
         & ZDMFEN(JL)*PVEN(JL,JK-1)-ZDMFDE(JL)*PVD(JL,JK-1)  
        PUD(JL,JK)=ZMFDUK*(1.0_dp/MIN(-RMFCMIN,PMFD(JL,JK)))
        PVD(JL,JK)=ZMFDVK*(1.0_dp/MIN(-RMFCMIN,PMFD(JL,JK)))
      ENDIF
    ENDDO
  ENDIF

ENDDO

IF (LHOOK) CALL DR_HOOK('CUDDRAFN',1,ZHOOK_HANDLE)
END SUBROUTINE CUDDRAFN

!==============================================================================

SUBROUTINE CUDLFSN &
 & (KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
 & KCBOT,    KCTOP,    LDLAND,   LDCUM,&
 & PTENH,    PQENH,    PUEN,     PVEN,&
 & PTEN,     PQSEN,    PGEO,&
 & PGEOH,    PAPH,     PTU,      PQU,      PLU,&
 & PUU,      PVU,      PMFUB,    PRFL,&
 & PTD,      PQD,      PUD,      PVD,&
 & PMFD,     PMFDS,    PMFDQ,    PDMFDP,&
 & KDTOP,    LDDRAF)  

!          THIS ROUTINE CALCULATES LEVEL OF FREE SINKING FOR
!          CUMULUS DOWNDRAFTS AND SPECIFIES T,Q,U AND V VALUES

!          M.TIEDTKE         E.C.M.W.F.    12/86 MODIF. 12/89

!          PURPOSE.
!          --------
!          TO PRODUCE LFS-VALUES FOR CUMULUS DOWNDRAFTS
!          FOR MASSFLUX CUMULUS PARAMETERIZATION

!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!          INPUT ARE ENVIRONMENTAL VALUES OF T,Q,U,V,P,PHI
!          AND UPDRAFT VALUES T,Q,U AND V AND ALSO
!          CLOUD BASE MASSFLUX AND CU-PRECIPITATION RATE.
!          IT RETURNS T,Q,U AND V VALUES AND MASSFLUX AT LFS.

!          METHOD.

!          CHECK FOR NEGATIVE BUOYANCY OF AIR OF EQUAL PARTS OF
!          MOIST ENVIRONMENTAL AIR AND CLOUD AIR.

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KTDIA*        START OF THE VERTICAL LOOP
!    *KLEV*         NUMBER OF LEVELS
!    *KCBOT*        CLOUD BASE LEVEL
!    *KCTOP*        CLOUD TOP LEVEL

!    INPUT PARAMETERS (LOGICAL):

!    *LDLAND*       LAND SEA MASK (.TRUE. FOR LAND)
!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS

!    INPUT PARAMETERS (REAL):

!    *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS          K
!    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS     KG/KG
!    *PUEN*         PROVISIONAL ENVIRONMENT U-VELOCITY (T+1)      M/S
!    *PVEN*         PROVISIONAL ENVIRONMENT V-VELOCITY (T+1)      M/S
!    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
!    *PQSEN*        ENVIRONMENT SPEC. SATURATION HUMIDITY (T+1)   KG/KG
!    *PGEO*         GEOPOTENTIAL                                  M2/S2
!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                  M2/S2
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS           PA
!    *PTU*          TEMPERATURE IN UPDRAFTS                        K
!    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                   KG/KG
!    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS             KG/KG
!    *PUU*          U-VELOCITY IN UPDRAFTS                        M/S
!    *PVU*          V-VELOCITY IN UPDRAFTS                        M/S
!    *PMFUB*        MASSFLUX IN UPDRAFTS AT CLOUD BASE           KG/(M2*S)

!    UPDATED PARAMETERS (REAL):

!    *PRFL*         PRECIPITATION RATE                           KG/(M2*S)

!    OUTPUT PARAMETERS (REAL):

!    *PTD*          TEMPERATURE IN DOWNDRAFTS                      K
!    *PQD*          SPEC. HUMIDITY IN DOWNDRAFTS                 KG/KG
!    *PUD*          U-VELOCITY IN DOWNDRAFTS                      M/S
!    *PVD*          V-VELOCITY IN DOWNDRAFTS                      M/S
!    *PMFD*         MASSFLUX IN DOWNDRAFTS                       KG/(M2*S)
!    *PMFDS*        FLUX OF DRY STATIC ENERGY IN DOWNDRAFTS       J/(M2*S)
!    *PMFDQ*        FLUX OF SPEC. HUMIDITY IN DOWNDRAFTS         KG/(M2*S)
!    *PDMFDP*       FLUX DIFFERENCE OF PRECIP. IN DOWNDRAFTS     KG/(M2*S)

!    OUTPUT PARAMETERS (INTEGER):

!    *KDTOP*        TOP LEVEL OF DOWNDRAFTS

!    OUTPUT PARAMETERS (LOGICAL):

!    *LDDRAF*       .TRUE. IF DOWNDRAFTS EXIST

!          EXTERNALS
!          ---------
!          *CUADJTQ* FOR CALCULATING WET BULB T AND Q AT LFS

!          MODIFICATIONS
!          -------------
!             92-09-21 : Update to Cy44      J.-J. MORCRETTE
!             99-06-04 : Optimisation        D.SALMOND 
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!----------------------------------------------------------------------

USE MESSY_CONVECT_ECMWF_PARAM   , ONLY : RCPD, RETV, RLVTT, RLSTT, RTT,            &
                                         LPHYLIN, RLPTRC, RLPAL1,                  &
                                         R2ES, R3LES, R3IES, R4LES,                &
                                         R4IES, R5LES, R5IES, R5ALVCP, R5ALSCP,    &
                                         RALVDCP, RALSDCP, RTWAT, RTICE, RTICECU,  &
                                         RTWAT_RTICE_R, RTWAT_RTICECU_R,           &
                                         RMFDEPS, LMFDD, LMFDUDV,                  & !from YOECUMF
                                         LHOOK, DR_HOOK,                           & !from YOMHOOK
                                         FOELHMCU                                    !functions

IMPLICIT NONE

INTEGER,INTENT(IN)    :: KLON 
INTEGER,INTENT(IN)    :: KLEV 
INTEGER,INTENT(IN)    :: KIDIA 
INTEGER,INTENT(IN)    :: KFDIA 
INTEGER               :: KTDIA ! Argument NOT used
INTEGER               :: KCBOT(KLON) ! Argument NOT used
INTEGER               :: KCTOP(KLON) ! Argument NOT used
LOGICAL                          :: LDLAND(KLON) ! Argument NOT used
LOGICAL                          :: LDCUM(KLON) ! Argument NOT used
REAL(dp)   ,INTENT(IN)    :: PTENH(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PQENH(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PUEN(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PVEN(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PTEN(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PQSEN(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PGEO(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PGEOH(KLON,KLEV+1) 
REAL(dp)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) 
REAL(dp)   ,INTENT(IN)    :: PTU(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PQU(KLON,KLEV) 
REAL(dp)                  :: PLU(KLON,KLEV) ! Argument NOT used
REAL(dp)   ,INTENT(IN)    :: PUU(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PVU(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PMFUB(KLON) 
REAL(dp)   ,INTENT(INOUT) :: PRFL(KLON) 
REAL(dp)   ,INTENT(INOUT)   :: PTD(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT)   :: PQD(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT)   :: PUD(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT)   :: PVD(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PMFD(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT)   :: PMFDS(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT)   :: PMFDQ(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT)   :: PDMFDP(KLON,KLEV) 
INTEGER,INTENT(INOUT)   :: KDTOP(KLON) 
LOGICAL           ,INTENT(INOUT)   :: LDDRAF(KLON) 
INTEGER ::            IKHSMIN(KLON)
REAL(dp) ::     ZTENWB(KLON,KLEV),      ZQENWB(KLON,KLEV),&
 & ZCOND(KLON),            ZPH(KLON),&
 & ZHSMIN(KLON)  
LOGICAL ::  LLO2(KLON)

INTEGER :: ICALL, IK, IKE, IS, JK, JL

REAL(dp) :: ZBUO, ZFOEEWI, ZFOEEWL, ZHSK, ZMFTOP, ZOEALFA,&
 & ZOELHM, ZQTEST, ZTARG, ZTTEST  
REAL(dp) :: ZHOOK_HANDLE

!----------------------------------------------------------------------

!     1.           SET DEFAULT VALUES FOR DOWNDRAFTS
!                  ---------------------------------

IF (LHOOK) CALL DR_HOOK('CUDLFSN',0,ZHOOK_HANDLE)
DO JL=KIDIA,KFDIA
  LDDRAF(JL)=.FALSE.
  KDTOP(JL)=KLEV+1
  IKHSMIN(JL)=KLEV+1
  ZHSMIN(JL)=1.E8_dp
ENDDO

!orig IF(.NOT.LMFDD) GO TO 300
IF(LMFDD) THEN

!----------------------------------------------------------------------

!     2.           DETERMINE LEVEL OF FREE SINKING:
!                  DOWNDRAFTS SHALL START AT MODEL LEVEL OF MINIMUM
!                  OF SATURATION MOIST STATIC ENERGY OR BELOW
!                  RESPECTIVELY

!                  FOR EVERY POINT AND PROCEED AS FOLLOWS:

!                    (1) DETERMINE LEVEL OF MINIMUM OF HS
!                    (2) DETERMINE WET BULB ENVIRONMENTAL T AND Q
!                    (3) DO MIXING WITH CUMULUS CLOUD AIR
!                    (4) CHECK FOR NEGATIVE BUOYANCY
!                    (5) IF BUOYANCY>0 REPEAT (2) TO (4) FOR NEXT
!                        LEVEL BELOW

!                  THE ASSUMPTION IS THAT AIR OF DOWNDRAFTS IS MIXTURE
!                  OF 50% CLOUD AIR + 50% ENVIRONMENTAL AIR AT WET BULB
!                  TEMPERATURE (I.E. WHICH BECAME SATURATED DUE TO
!                  EVAPORATION OF RAIN AND CLOUD WATER)
!                  ----------------------------------------------------

  DO JK=3,KLEV-2

    IF (LPHYLIN) THEN

      DO JL=KIDIA,KFDIA
        ZTARG=PTEN(JL,JK)
        ZOEALFA=0.5_dp*(TANH(RLPAL1*(ZTARG-RLPTRC))+1.0_dp)
        ZFOEEWL=R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
        ZFOEEWI=R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
        ZOELHM =ZOEALFA*ZFOEEWL+(1.0_dp-ZOEALFA)*ZFOEEWI
        ZHSK=RCPD*PTEN(JL,JK)+PGEO(JL,JK)+ZOELHM             *PQSEN(JL,JK)
        IF(ZHSK < ZHSMIN(JL)) THEN
          ZHSMIN(JL)=ZHSK
          IKHSMIN(JL)=JK
        ENDIF
      ENDDO

    ELSE

      DO JL=KIDIA,KFDIA
        ZHSK=RCPD*PTEN(JL,JK)+PGEO(JL,JK)+FOELHMCU(PTEN(JL,JK))*PQSEN(JL,JK)
        IF(ZHSK < ZHSMIN(JL)) THEN
          ZHSMIN(JL)=ZHSK
          IKHSMIN(JL)=JK
        ENDIF
      ENDDO

    ENDIF

  ENDDO
  IKE=KLEV-3
  DO JK=3,IKE

!     2.1          CALCULATE WET-BULB TEMPERATURE AND MOISTURE
!                  FOR ENVIRONMENTAL AIR IN *CUADJTQ*
!                  -------------------------------------------

    IS=0
    DO JL=KIDIA,KFDIA
      ZTENWB(JL,JK)=PTENH(JL,JK)
      ZQENWB(JL,JK)=PQENH(JL,JK)
      ZPH(JL)=PAPH(JL,JK)
      LLO2(JL)=LDCUM(JL).AND.PRFL(JL) > 0.0_dp.AND..NOT.LDDRAF(JL).AND.&
       & (JK < KCBOT(JL).AND.JK > KCTOP(JL)).AND.&
       & JK >= IKHSMIN(JL)  
      IF(LLO2(JL))THEN
        IS=IS+1
      ENDIF
    ENDDO
!orig   IF(IS.EQ.0) GO TO 290
    IF(IS == 0) CYCLE

    IK=JK
    ICALL=2
    CALL CUADJTQ &
     & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
     & IK,&
     & ZPH,      ZTENWB,   ZQENWB,   LLO2,     ICALL)  

!     2.2          DO MIXING OF CUMULUS AND ENVIRONMENTAL AIR
!                  AND CHECK FOR NEGATIVE BUOYANCY.
!                  THEN SET VALUES FOR DOWNDRAFT AT LFS.
!                  ----------------------------------------

!DIR$ IVDEP
!OCL NOVREC
    DO JL=KIDIA,KFDIA
      IF(LLO2(JL)) THEN
        ZTTEST=0.5_dp*(PTU(JL,JK)+ZTENWB(JL,JK))
        ZQTEST=0.5_dp*(PQU(JL,JK)+ZQENWB(JL,JK))
        ZBUO=ZTTEST*(1.0_dp+RETV  *ZQTEST)-&
         & PTENH(JL,JK)*(1.0_dp+RETV  *PQENH(JL,JK))  
        ZCOND(JL)=PQENH(JL,JK)-ZQENWB(JL,JK)
        ZMFTOP=-RMFDEPS*PMFUB(JL)
        IF(ZBUO < 0.0_dp.AND.PRFL(JL) > 10._dp*ZMFTOP*ZCOND(JL)) THEN
          KDTOP(JL)=JK
          LDDRAF(JL)=.TRUE.
          PTD(JL,JK)=ZTTEST
          PQD(JL,JK)=ZQTEST
          PMFD(JL,JK)=ZMFTOP
          PMFDS(JL,JK)=PMFD(JL,JK)*(RCPD*PTD(JL,JK)+PGEOH(JL,JK))
          PMFDQ(JL,JK)=PMFD(JL,JK)*PQD(JL,JK)
          PDMFDP(JL,JK-1)=-0.5_dp*PMFD(JL,JK)*ZCOND(JL)
          PRFL(JL)=PRFL(JL)+PDMFDP(JL,JK-1)
        ENDIF
      ENDIF
    ENDDO
    IF(LMFDUDV) THEN
      DO JL=KIDIA,KFDIA
        IF(PMFD(JL,JK) < 0.0_dp) THEN
          PUD(JL,JK)=0.5_dp*(PUU(JL,JK)+PUEN(JL,JK-1))
          PVD(JL,JK)=0.5_dp*(PVU(JL,JK)+PVEN(JL,JK-1))
        ENDIF
      ENDDO
    ENDIF

! 290   continue
  ENDDO

!300  CONTINUE
ENDIF

IF (LHOOK) CALL DR_HOOK('CUDLFSN',1,ZHOOK_HANDLE)
END SUBROUTINE CUDLFSN

!==============================================================================

SUBROUTINE CUDTDQN &
 & (  KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
 & KTOPM2,   KTYPE,    LDLAND,   LDCUM,    PTSPHY,&
 & PAPH,     PTEN,     PLGLAC,   PLUDE,&
 & PMFUS,    PMFDS,    PMFUQ,    PMFDQ,&
 & PMFUL,    PDMFUP,   PDMFDP,   PDPMEL,&
 & PTENT,    PTENQ,    PENTH )  

!**** *CUDTDQ* - UPDATES T AND Q TENDENCIES, PRECIPITATION RATES
!                DOES GLOBAL DIAGNOSTICS

!          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89

!**   INTERFACE.
!     ----------

!          *CUDTDQ* IS CALLED FROM *CUMASTR*

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KTDIA*        START OF THE VERTICAL LOOP
!    *KLEV*         NUMBER OF LEVELS

!    INPUT PARAMETERS (LOGICAL): 

!    *LDLAND*       LAND SEA MASK (.TRUE. FOR LAND)
!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS 

!    INPUT PARAMETERS (REAL):

!    *PTSPHY*       TIME STEP FOR THE PHYSICS                       S
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS            PA
!    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
!    *PLGLAC*       FLUX OF FROZEN CLOUDWATER IN UPDRAFTS         KG/(M2*S) 
!    *PLUDE*        DETRAINED LIQUID WATER                        KG/(M3*S)
!    *PMFUS*        FLUX OF DRY STATIC ENERGY IN UPDRAFTS          J/(M2*S)
!    *PMFDS*        FLUX OF DRY STATIC ENERGY IN DOWNDRAFTS        J/(M2*S)
!    *PMFUQ*        FLUX OF SPEC. HUMIDITY IN UPDRAFTS            KG/(M2*S)
!    *PMFDQ*        FLUX OF SPEC. HUMIDITY IN DOWNDRAFTS          KG/(M2*S)
!    *PMFUL*        FLUX OF LIQUID WATER IN UPDRAFTS              KG/(M2*S)
!    *PDMFUP*       FLUX DIFFERENCE OF PRECIP. IN UPDRAFTS        KG/(M2*S)
!    *PDMFDP*       FLUX DIFFERENCE OF PRECIP. IN DOWNDRAFTS      KG/(M2*S)
!    *PDPMEL*       CHANGE IN PRECIP.-FLUXES DUE TO MELTING       KG/(M2*S)

!    UPDATED PARAMETERS (REAL):

!    *PTENT*        TEMPERATURE TENDENCY                           K/S
!    *PTENQ*        MOISTURE TENDENCY                             KG/(KG S)

!    OUTPUT PARAMETERS (REAL):

!    *PENTH*        INCREMENT OF DRY STATIC ENERGY                 J/(KG*S) 

!----------------------------------------------------------------------

!               MODIFICATIONS
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!       96-09-20       : changed so can be used with diagnost
!                        cloud scheme      D.GREGORY
!       99-06-04       : Optimisation      D.SALMOND   
!       03-08-28       : Clean up LINPHYS  P.BECHTOLD

!----------------------------------------------------------------------

USE MESSY_CONVECT_ECMWF_PARAM   , ONLY : RG, RCPD, RLVTT, RLSTT, RLMLT, RTT,   &
                                         LPHYLIN, RLPTRC, RLPAL1,              &
                                         R2ES, R3LES, R3IES, R4LES,            &
                                         R4IES, R5LES, R5IES, R5ALVCP, R5ALSCP,&
                                         RALVDCP, RALSDCP, RTWAT, RTICE,       &
                                         RTICECU, RTWAT_RTICE_R,               &
                                         RTWAT_RTICECU_R,                      &
                                         LHOOK, DR_HOOK,                       &
                                         LEPCLD, FOELHMCU    !functions
IMPLICIT NONE

INTEGER,INTENT(IN)    :: KLON 
INTEGER,INTENT(IN)    :: KLEV 
INTEGER,INTENT(IN)    :: KIDIA 
INTEGER,INTENT(IN)    :: KFDIA 
INTEGER               :: KTDIA ! Argument NOT used
INTEGER,INTENT(IN)    :: KTOPM2 
INTEGER,INTENT(IN)    :: KTYPE(KLON) 
LOGICAL                          :: LDLAND(KLON) ! Argument NOT used
LOGICAL           ,INTENT(INOUT) :: LDCUM(KLON) 
REAL(dp)                  :: PTSPHY ! Argument NOT used
REAL(dp)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) 
REAL(dp)   ,INTENT(IN)    :: PTEN(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PLGLAC(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PLUDE(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PMFUS(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PMFDS(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PMFUQ(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PMFDQ(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PMFUL(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PDMFUP(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PDMFDP(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PDPMEL(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PTENT(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PTENQ(KLON,KLEV) 
REAL(dp)   ,INTENT(OUT)   :: PENTH(KLON,KLEV) 
INTEGER :: JK, JL

REAL(dp) :: ZALV, ZDQDT, ZDTDT, ZOEALFA, ZTARG
REAL(dp) :: ZHOOK_HANDLE

!----------------------------------------------------------------------

!*    1.0          INCREMENTATION OF T AND Q TENDENCIES
!                  ------------------------------------

IF (LHOOK) CALL DR_HOOK('CUDTDQN',0,ZHOOK_HANDLE)
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    PENTH(JL,JK)=0.0_dp
  ENDDO
ENDDO

!         MASS-FLUX APPROACH SWITCHED ON FOR DEEP CONVECTION ONLY
!         IN THE TANGENT-LINEAR AND ADJOINT VERSIONS

do jl=kidia,kfdia
  if (ktype(jl) /= 1.and.lphylin) ldcum(jl)=.false.
enddo

! zero detrained liquid water if diagnostic cloud scheme to be used

! this means that detrained liquid water will be evaporated in the 
! cloud environment and not fed directly into a cloud liquid water 
! variable

if (.not.lepcld.or.lphylin) then
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      plude(jl,jk)=0.0_dp
    ENDDO
  ENDDO
ENDIF

DO JK=KTOPM2,KLEV

  IF(JK < KLEV) THEN
    DO JL=KIDIA,KFDIA
      IF(LDCUM(JL)) THEN
        IF (LPHYLIN) THEN
          ZTARG=PTEN(JL,JK)
          ZOEALFA=0.5_dp*(TANH(RLPAL1*(ZTARG-RLPTRC))+1.0_dp)
          ZALV=ZOEALFA*RLVTT+(1.0_dp-ZOEALFA)*RLSTT
        ELSE
          ZALV=FOELHMCU(PTEN(JL,JK))
        ENDIF
        ZDTDT=(RG/(PAPH(JL,JK+1)-PAPH(JL,JK)))/RCPD*&
         & (PMFUS(JL,JK+1)-PMFUS(JL,JK)+&
         & PMFDS(JL,JK+1)-PMFDS(JL,JK)&
         & +RLMLT*PLGLAC(JL,JK)&
         & -RLMLT*PDPMEL(JL,JK)&
         & -ZALV*(PMFUL(JL,JK+1)-PMFUL(JL,JK)-&
         & PLUDE(JL,JK)-&
         & (PDMFUP(JL,JK)+PDMFDP(JL,JK))))  

! increment of dry static energy
        PENTH(JL,JK)= ZDTDT*RCPD

        PTENT(JL,JK)=PTENT(JL,JK)+ZDTDT
        ZDQDT=(RG/(PAPH(JL,JK+1)-PAPH(JL,JK)))*&
         & (PMFUQ(JL,JK+1)-PMFUQ(JL,JK)+&
         & PMFDQ(JL,JK+1)-PMFDQ(JL,JK)+&
         & PMFUL(JL,JK+1)-PMFUL(JL,JK)-&
         & PLUDE(JL,JK)-&
         & (PDMFUP(JL,JK)+PDMFDP(JL,JK)))  
        PTENQ(JL,JK)=PTENQ(JL,JK)+ZDQDT
      ENDIF
    ENDDO

  ELSE
    DO JL=KIDIA,KFDIA
      IF(LDCUM(JL)) THEN
        IF (LPHYLIN) THEN
          ZTARG=PTEN(JL,JK)
          ZOEALFA=0.5_dp*(TANH(RLPAL1*(ZTARG-RLPTRC))+1.0_dp)
          ZALV=ZOEALFA*RLVTT+(1.0_dp-ZOEALFA)*RLSTT
        ELSE
          ZALV=FOELHMCU(PTEN(JL,JK))
        ENDIF
        ZDTDT=-(RG/(PAPH(JL,JK+1)-PAPH(JL,JK)))/RCPD*&
         & (PMFUS(JL,JK)+PMFDS(JL,JK)+RLMLT*PDPMEL(JL,JK)-ZALV*&
         & (PMFUL(JL,JK)+PDMFUP(JL,JK)+PDMFDP(JL,JK)))  
        PTENT(JL,JK)=PTENT(JL,JK)+ZDTDT

! increment of dry static energy
        PENTH(JL,JK)= ZDTDT*RCPD

        ZDQDT=-(RG/(PAPH(JL,JK+1)-PAPH(JL,JK)))*&
         & (PMFUQ(JL,JK)+PMFDQ(JL,JK)+&
         & (PMFUL(JL,JK)+PDMFUP(JL,JK)+PDMFDP(JL,JK)))  
        PTENQ(JL,JK)=PTENQ(JL,JK)+ZDQDT
      ENDIF
    ENDDO
  ENDIF

ENDDO

IF (LHOOK) CALL DR_HOOK('CUDTDQN',1,ZHOOK_HANDLE)
END SUBROUTINE CUDTDQN

!==============================================================================

SUBROUTINE CUDUDV &
 & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
 & KTOPM2,   KTYPE,    KCBOT,    KCTOP,    LDCUM,    PTSPHY,&
 & PAPH,     PUEN,     PVEN,     PMFU,     PMFD,&
 & PUU,      PUD,      PVU,      PVD,&
 & PTENU,    PTENV  )  

!**** *CUDUDV* - UPDATES U AND V TENDENCIES,
!                DOES GLOBAL DIAGNOSTIC OF DISSIPATION

!          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
!          P.BECHTOLD        E.C.M.W.F.    11/02/05 IMPLICIT SOLVER

!**   INTERFACE.
!     ----------

!          *CUDUDV* IS CALLED FROM *CUMASTR*

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KTDIA*        START OF THE VERTICAL LOOP
!    *KLEV*         NUMBER OF LEVELS
!    *KTYPE*        TYPE OF CONVECTION
!                       1 = PENETRATIVE CONVECTION
!                       2 = SHALLOW CONVECTION
!                       3 = MIDLEVEL CONVECTION
!    *KCBOT*        CLOUD BASE LEVEL
!    *KCTOP*        CLOUD TOP LEVEL

!    INPUT PARAMETERS (LOGICAL):

!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS 

!    INPUT PARAMETERS (REAL):

!    *PTSPHY*       TIME STEP FOR THE PHYSICS                      S
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS            PA
!    *PUEN*         PROVISIONAL ENVIRONMENT U-VELOCITY (T+1)       M/S
!    *PVEN*         PROVISIONAL ENVIRONMENT V-VELOCITY (T+1)       M/S
!    *PMFU*         MASSFLUX UPDRAFTS                             KG/(M2*S)
!    *PMFD*         MASSFLUX DOWNDRAFTS                           KG/(M2*S)
!    *PUU*          U-VELOCITY IN UPDRAFTS                         M/S
!    *PUD*          U-VELOCITY IN DOWNDRAFTS                       M/S
!    *PVU*          V-VELOCITY IN UPDRAFTS                         M/S
!    *PVD*          V-VELOCITY IN DOWNDRAFTS                       M/S

!    UPDATED PARAMETERS (REAL):

!    *PTENU*        TENDENCY OF U-COMP. OF WIND                    M/S2
!    *PTENV*        TENDENCY OF V-COMP. OF WIND                    M/S2

!            METHOD
!            -------
!       EXPLICIT UPSTREAM AND IMPLICIT SOLUTION OF VERTICAL ADVECTION
!       DEPENDING ON VALUE OF RMFSOLUV: 
!       0=EXPLICIT 0-1 SEMI-IMPLICIT >=1 IMPLICIT 
!
!       FOR EXPLICIT SOLUTION: ONLY ONE SINGLE ITERATION
!       FOR IMPLICIT SOLUTION: FIRST IMPLICIT SOLVER, THEN EXPLICIT SOLVER
!                              TO CORRECT TENDENCIES BELOW CLOUD BASE
!
!            EXTERNALS
!            ---------
!            CUBIDIAG

!          MODIFICATIONS
!          -------------
!             92-09-21 : Update to Cy44      J.-J. MORCRETTE
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!----------------------------------------------------------------------

USE MESSY_CONVECT_ECMWF_PARAM   , ONLY : RG,                &
                                         RMFSOLUV,          &  !from YOECUMF
                                         LHOOK, DR_HOOK        !from YOMHOOK


IMPLICIT NONE

INTEGER,INTENT(IN)    :: KLON 
INTEGER,INTENT(IN)    :: KLEV 
INTEGER,INTENT(IN)    :: KIDIA 
INTEGER,INTENT(IN)    :: KFDIA 
INTEGER               :: KTDIA ! Argument NOT used
INTEGER,INTENT(IN)    :: KTOPM2 
INTEGER,INTENT(IN)    :: KTYPE(KLON) 
INTEGER,INTENT(IN)    :: KCBOT(KLON) 
INTEGER,INTENT(IN)    :: KCTOP(KLON) 
LOGICAL           ,INTENT(IN)    :: LDCUM(KLON) 
REAL(dp)   ,INTENT(IN)    :: PTSPHY
REAL(dp)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) 
REAL(dp)   ,INTENT(IN)    :: PUEN(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PVEN(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PMFU(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PMFD(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PUU(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PUD(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PVU(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PVD(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PTENU(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PTENV(KLON,KLEV) 
REAL(dp) :: ZUEN(KLON,KLEV),     ZVEN(KLON,KLEV),&
 & ZMFUU(KLON,KLEV),    ZMFDU(KLON,KLEV),&
 & ZMFUV(KLON,KLEV),    ZMFDV(KLON,KLEV)

INTEGER :: IK, IKB, JK, JL, JIT, JITER

REAL(dp) :: ZZP, ZIMP, ZTSPHY
!       ALLOCATABLE ARAYS
REAL(dp),   DIMENSION(:,:), ALLOCATABLE :: ZDUDT, ZDVDT, ZDP
REAL(dp),   DIMENSION(:,:), ALLOCATABLE :: ZB,  ZR1,  ZR2
LOGICAL,DIMENSION(:,:),   ALLOCATABLE :: LLCUMBAS
REAL(dp) :: ZHOOK_HANDLE

!----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CUDUDV',0,ZHOOK_HANDLE)
ZIMP   = 1.0_dp - RMFSOLUV
ZTSPHY = 1.0_dp / PTSPHY

ALLOCATE(ZDUDT(KLON,KLEV))
ALLOCATE(ZDVDT(KLON,KLEV))
ALLOCATE(ZDP(KLON,KLEV))

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    IF(LDCUM(JL)) THEN
      ZUEN(JL,JK)=PUEN(JL,JK)
      ZVEN(JL,JK)=PVEN(JL,JK)
      ZDP(JL,JK)=RG/(PAPH(JL,JK+1)-PAPH(JL,JK))
    ENDIF
  ENDDO
ENDDO

JITER=1
!IF(RMFSOLUV>0.0_dp) JITER=2
!----------------------------------------------------------------------

DO JIT=1,JITER

  IF(JIT==2) ZIMP=1._dp

!*    1.0          CALCULATE FLUXES AND UPDATE U AND V TENDENCIES
!                  ----------------------------------------------


  DO JK=KTOPM2,KLEV
    IK=JK-1
    DO JL=KIDIA,KFDIA
      IF(LDCUM(JL)) THEN
        ZMFUU(JL,JK)=PMFU(JL,JK)*(PUU(JL,JK)-ZIMP*ZUEN(JL,IK))
        ZMFUV(JL,JK)=PMFU(JL,JK)*(PVU(JL,JK)-ZIMP*ZVEN(JL,IK))
        ZMFDU(JL,JK)=PMFD(JL,JK)*(PUD(JL,JK)-ZIMP*ZUEN(JL,IK))
        ZMFDV(JL,JK)=PMFD(JL,JK)*(PVD(JL,JK)-ZIMP*ZVEN(JL,IK))
      ENDIF
    ENDDO
  ENDDO
  
! linear fluxes below cloud
  if(ABS(rmfsoluv).lt.tiny(0.0_dp)) then
    DO JK=KTOPM2,KLEV
!DIR$ IVDEP
!OCL NOVREC
      DO JL=KIDIA,KFDIA
        IF(LDCUM(JL).AND.JK > KCBOT(JL)) THEN
          IKB=KCBOT(JL)
          ZZP=((PAPH(JL,KLEV+1)-PAPH(JL,JK))/(PAPH(JL,KLEV+1)-PAPH(JL,IKB)))
          IF(KTYPE(JL) == 3) THEN
            ZZP=ZZP*ZZP
          ENDIF
          ZMFUU(JL,JK)=ZMFUU(JL,IKB)*ZZP
          ZMFUV(JL,JK)=ZMFUV(JL,IKB)*ZZP
          ZMFDU(JL,JK)=ZMFDU(JL,IKB)*ZZP
          ZMFDV(JL,JK)=ZMFDV(JL,IKB)*ZZP
        ENDIF
      ENDDO
    ENDDO
  endif

!*    1.2          COMPUTE TENDENCIES
!                  ------------------

  DO JK=KTOPM2,KLEV
  
    IF(JK < KLEV) THEN
      IK=JK+1
      DO JL=KIDIA,KFDIA
        IF(LDCUM(JL)) THEN
          ZDUDT(JL,JK)=ZDP(JL,JK)*&
           & (ZMFUU(JL,IK)-ZMFUU(JL,JK)+ZMFDU(JL,IK)-ZMFDU(JL,JK))  
          ZDVDT(JL,JK)=ZDP(JL,JK)*&
           & (ZMFUV(JL,IK)-ZMFUV(JL,JK)+ZMFDV(JL,IK)-ZMFDV(JL,JK))  
        ENDIF
      ENDDO
  
    ELSE
      DO JL=KIDIA,KFDIA
        IF(LDCUM(JL)) THEN
          ZDUDT(JL,JK)=-ZDP(JL,JK)*(ZMFUU(JL,JK)+ZMFDU(JL,JK))
          ZDVDT(JL,JK)=-ZDP(JL,JK)*(ZMFUV(JL,JK)+ZMFDV(JL,JK))
        ENDIF
      ENDDO
    ENDIF
  
  ENDDO

  IF ( ZIMP==1.0_dp ) THEN

!*    1.3          UPDATE TENDENCIES
!                  -----------------

    DO JK=KTOPM2,KLEV
      DO JL=KIDIA,KFDIA
        IF(LDCUM(JL)) THEN
          PTENU(JL,JK)=PTENU(JL,JK)+ZDUDT(JL,JK)
          PTENV(JL,JK)=PTENV(JL,JK)+ZDVDT(JL,JK)
        ENDIF
      ENDDO
    ENDDO

  ELSE
!----------------------------------------------------------------------
      
!*      1.6          IMPLICIT SOLUTION
!                    -----------------
   
     ! Fill bi-diagonal Matrix vectors A=k-1, B=k;
     ! reuse ZMFUU=A and ZB=B; 
     ! ZDUDT and ZDVDT correspond to the RHS ("constants") of the equation
     ! The solution is in ZR1 and ZR2
     
      ALLOCATE(ZB(KLON,KLEV))
      ALLOCATE(ZR1(KLON,KLEV))
      ALLOCATE(ZR2(KLON,KLEV))
      ALLOCATE(LLCUMBAS(KLON,KLEV))
      LLCUMBAS(:,:)=.FALSE.
      ZB(:,:)=1.0_dp
      ZMFUU(:,:)=0.0_dp
     
     ! Fill vectors A, B and RHS 
     
      DO JK=KTOPM2,KLEV
         IK=JK+1
         DO JL=KIDIA,KFDIA
          ! LLCUMBAS(JL,JK)=LDCUM(JL).AND.JK>=KCTOP(JL)-1.AND.JK<=KCBOT(JL)
           LLCUMBAS(JL,JK)=LDCUM(JL).AND.JK>=KCTOP(JL)-1 
           IF(LLCUMBAS(JL,JK)) THEN
             ZZP=RMFSOLUV*ZDP(JL,JK)*PTSPHY
             ZMFUU(JL,JK)=-ZZP*(PMFU(JL,JK)+PMFD(JL,JK))
             ZDUDT(JL,JK) = ZDUDT(JL,JK)*PTSPHY+ZUEN(JL,JK)
             ZDVDT(JL,JK) = ZDVDT(JL,JK)*PTSPHY+ZVEN(JL,JK)
             IF(JK<KLEV) THEN
               ZB(JL,JK)=1._dp+ZZP*(PMFU(JL,IK)+PMFD(JL,IK))
             ELSE
               ZB(JL,JK)=1._dp
             ENDIF
           ENDIF
         ENDDO
      ENDDO
     
      CALL CUBIDIAG&
         &( KIDIA, KFDIA, KLON, KLEV &
         &, KCTOP, LLCUMBAS &
         &, ZMFUU,    ZB,    ZDUDT,   ZR1 )
     
      CALL CUBIDIAG&
         &( KIDIA, KFDIA, KLON, KLEV &
         &, KCTOP, LLCUMBAS &
         &, ZMFUU,    ZB,    ZDVDT,   ZR2 )
     
      ! Compute tendencies
      ! If second iteration was used, the values below cloud base would be
      ! replaced by the explicit solution
     
      DO JK=KTOPM2,KLEV
         DO JL=KIDIA,KFDIA
           IF(LLCUMBAS(JL,JK)) THEN
             PTENU(JL,JK)=PTENU(JL,JK)+(ZR1(JL,JK)-ZUEN(JL,JK))*ZTSPHY
             PTENV(JL,JK)=PTENV(JL,JK)+(ZR2(JL,JK)-ZVEN(JL,JK))*ZTSPHY
             ! ZUEN(JL,JK)=ZR1(JL,JK)
             ! ZVEN(JL,JK)=ZR2(JL,JK)
           ENDIF
         ENDDO
      ENDDO
     
      DEALLOCATE(LLCUMBAS)
      DEALLOCATE(ZR2)
      DEALLOCATE(ZR1)
      DEALLOCATE(ZB)
     
   ENDIF
!----------------------------------------------------------------------
ENDDO

DEALLOCATE(ZDP)
DEALLOCATE(ZDVDT)
DEALLOCATE(ZDUDT)

IF (LHOOK) CALL DR_HOOK('CUDUDV',1,ZHOOK_HANDLE)
END SUBROUTINE CUDUDV

!==============================================================================

SUBROUTINE CUENTR &
 & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
 & KK,       KLWMIN,   KTYPE,    KCBOT,    KCTOP0,&
 & LDCUM,    LDWORK,&
 & PTENH,    PQENH,    PTENQ,    PAPH,     PAP,&
 & PMFU,     PENTR,&
 & PCBASE,   PDMFEN,   PDMFDE )  

!          M.TIEDTKE         E.C.M.W.F.     12/89

!          PURPOSE.
!          --------
!          THIS ROUTINE CALCULATES ENTRAINMENT/DETRAINMENT RATES
!          FOR UPDRAFTS IN CUMULUS PARAMETERIZATION

!          INTERFACE
!          ---------

!          THIS ROUTINE IS CALLED FROM *CUASC*.
!          INPUT ARE ENVIRONMENTAL VALUES T,Q ETC
!          AND UPDRAFT VALUES T,Q ETC
!          IT RETURNS ENTRAINMENT/DETRAINMENT RATES

!          METHOD.
!          --------
!          M. TIEDTKE (1989)

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KTDIA*        START OF THE VERTICAL LOOP
!    *KLEV*         NUMBER OF LEVELS
!    *KK*           CURRENT LEVEL
!    *KLWMIN*       LEVEL OF MAXIMUM VERTICAL VELOCITY
!    *KTYPE*        TYPE OF CONVECTION
!                       1 = PENETRATIVE CONVECTION
!                       2 = SHALLOW CONVECTION
!                       3 = MIDLEVEL CONVECTION
!    *KCBOT*        CLOUD BASE LEVEL
!    *KCTOP0*       FIRST GUESS OF CLOUD TOP LEVEL

!    INPUT PARAMETERS (LOGICAL):

!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS

!    INPUT PARAMETERS (REAL):

!    *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS         K
!    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS    KG/KG
!    *PTENQ*        MOISTURE TENDENCY                           KG/(KG S)
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS          PA
!    *PAP*          PROVISIONAL PRESSURE ON FULL LEVELS          PA
!    *PMFU*         MASSFLUX IN UPDRAFTS                        KG/(M2*S)
!    *PENTR*        FRACTIONAL MASS ENTRAINMENT RATE             1/M

!    OUTPUT PARAMETERS (REAL):

!    *PCBASE*       PRESSURE AT CLOUD BASE                       PA
!    *PDMFEN*       ENTRAINMENT RATE                            KG/(M2*S)
!    *PDMFDE*       DETRAINMENT RATE                            KG/(M2*S)

!          EXTERNALS
!          ---------
!          NONE

!----------------------------------------------------------------------

USE MESSY_CONVECT_ECMWF_PARAM   , ONLY : RG, RD,        &
                                         RDEPTHS,       &    !from YOECUMF
                                         LHOOK, DR_HOOK      !from YOMHOOK

IMPLICIT NONE

INTEGER,INTENT(IN)    :: KLON 
INTEGER,INTENT(IN)    :: KLEV 
INTEGER,INTENT(IN)    :: KIDIA 
INTEGER,INTENT(IN)    :: KFDIA 
INTEGER               :: KTDIA ! Argument NOT used
INTEGER,INTENT(IN)    :: KK 
INTEGER,INTENT(IN)    :: KLWMIN(KLON) 
INTEGER,INTENT(IN)    :: KTYPE(KLON) 
INTEGER,INTENT(IN)    :: KCBOT(KLON) 
INTEGER,INTENT(IN)    :: KCTOP0(KLON) 
LOGICAL           ,INTENT(IN)    :: LDCUM(KLON) 
LOGICAL           ,INTENT(IN)    :: LDWORK 
REAL(dp)   ,INTENT(IN)    :: PTENH(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PQENH(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PTENQ(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) 
REAL(dp)   ,INTENT(IN)    :: PAP(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PMFU(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PENTR(KLON) 
REAL(dp)   ,INTENT(OUT)   :: PCBASE(KLON) 
REAL(dp)   ,INTENT(OUT)   :: PDMFEN(KLON) 
REAL(dp)   ,INTENT(OUT)   :: PDMFDE(KLON) 
REAL(dp) ::     ZRHO(KLON),             ZPTOP(KLON)

LOGICAL ::  LLO1,LLO2

INTEGER :: IKLWMIN, JL, IK1,IK2

REAL(dp) :: ZDPRHO, ZENTR, ZENTR2, ZPMID
REAL(dp) :: ZHOOK_HANDLE

!----------------------------------------------------------------------

!*    1.           CALCULATE ENTRAINMENT AND DETRAINMENT RATES
!                  -------------------------------------------

IF (LHOOK) CALL DR_HOOK('CUENTR',0,ZHOOK_HANDLE)
IF(LDWORK) THEN

  DO JL=KIDIA,KFDIA
    PDMFEN(JL)=0.0_dp
    PDMFDE(JL)=0.0_dp
    ZRHO(JL)=PAPH(JL,KK+1)/(RD*PTENH(JL,KK+1))
    IK1=max(1,KCBOT(JL))
    IK2=max(1,KCTOP0(JL))
    PCBASE(JL)=PAPH(JL,IK1)
    ZPTOP(JL)=PAPH(JL,IK2)
  ENDDO

!*    1.1          SPECIFY ENTRAINMENT RATES FOR SHALLOW CLOUDS
!                  --------------------------------------------

!*    1.2          SPECIFY ENTRAINMENT RATES FOR DEEP CLOUDS
!                  -----------------------------------------

  DO JL=KIDIA,KFDIA
    IF(LDCUM(JL)) THEN
      ZDPRHO=(PAPH(JL,KK+1)-PAPH(JL,KK))/(RG*ZRHO(JL))
      ZENTR=PENTR(JL)*PMFU(JL,KK+1)*ZDPRHO
      LLO1=KK < KCBOT(JL)
      IF(LLO1) PDMFDE(JL)=ZENTR
      ZPMID=0.5_dp*(PCBASE(JL)+ZPTOP(JL))
      LLO2=LLO1.AND.KTYPE(JL) == 2.AND.&
       & (PCBASE(JL)-PAPH(JL,KK) < RDEPTHS.OR.&
       & PAPH(JL,KK) > ZPMID)  
      IF(LLO2) PDMFEN(JL)=ZENTR
      IKLWMIN=MAX(KLWMIN(JL),KCTOP0(JL)+2)
      LLO2=LLO1.AND.(KTYPE(JL) == 1.OR.KTYPE(JL) == 3).AND.&
       & (KK >= IKLWMIN.OR.PAP(JL,KK) > ZPMID)  
      IF(LLO2) PDMFEN(JL)=ZENTR
      LLO1=PDMFEN(JL) > 0.0_dp.AND.(KTYPE(JL) == 1.OR.KTYPE(JL) == 2)
      IF(LLO1) THEN
        ZENTR2=&
         & (1.0_dp+3._dp*(1.0_dp-MIN(1.0_dp,(PCBASE(JL)-PAP(JL,KK))/ &
         & 1.5E4_dp)))  
        ZENTR=ZENTR*ZENTR2
        PDMFEN(JL)=PDMFEN(JL)*ZENTR2
        PDMFDE(JL)=PDMFDE(JL)*ZENTR2
      ENDIF
      IF(LLO2.AND.PQENH(JL,KK+1) > 1.E-5_dp)&
       & PDMFEN(JL)=ZENTR+MAX(PTENQ(JL,KK),0.0_dp)/PQENH(JL,KK+1)*&
       & ZRHO(JL)*ZDPRHO  
    ENDIF
  ENDDO

ENDIF

IF (LHOOK) CALL DR_HOOK('CUENTR',1,ZHOOK_HANDLE)
END SUBROUTINE CUENTR

!==============================================================================

!OPTIONS XOPT(HSFUN)
SUBROUTINE CUFLXN &
 & (  KIDIA,    KFDIA,    KLON,     KTDIA,   KLEV,&
 & PTSPHY,&
 & PTEN,     PQEN,     PQSEN,    PTENH,    PQENH,&
 & PAPH,     PAP,      PGEOH,    LDLAND,   LDCUM,&
 & KCBOT,    KCTOP,    KDTOP,    KTOPM2,&
 & KTYPE,    LDDRAF,&
 & PMFU,     PMFD,     PMFUS,    PMFDS,&
 & PMFUQ,    PMFDQ,    PMFUL,    PLUDE,&
 & PDMFUP,   PDMFDP,   PDPMEL,   PLGLAC,&
 & PMFLXR,   PMFLXS,   PRAIN,    PMFDDE_RATE )  

!          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89

!          PURPOSE
!          -------

!          THIS ROUTINE DOES THE FINAL CALCULATION OF CONVECTIVE
!          FLUXES IN THE CLOUD LAYER AND IN THE SUBCLOUD LAYER

!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KTDIA*        START OF THE VERTICAL LOOP
!    *KLEV*         NUMBER OF LEVELS
!    *KCBOT*        CLOUD BASE LEVEL
!    *KCTOP*        CLOUD TOP LEVEL
!    *KDTOP*        TOP LEVEL OF DOWNDRAFTS

!    INPUT PARAMETERS (LOGICAL):

!    *LDLAND*       LAND SEA MASK (.TRUE. FOR LAND)
!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS

!    INPUT PARAMETERS (REAL):

!    *PTSPHY*       TIME STEP FOR THE PHYSICS                       S
!    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
!    *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1)  KG/KG
!    *PQSEN*        ENVIRONMENT SPEC. SATURATION HUMIDITY (T+1)   KG/KG
!    *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS           K
!    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS      KG/KG
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS            PA
!    *PAP*          PROVISIONAL PRESSURE ON FULL LEVELS            PA
!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2

!    UPDATED PARAMETERS (INTEGER):

!    *KTYPE*        SET TO ZERO IF LDCUM=.FALSE.

!    UPDATED PARAMETERS (LOGICAL):

!    *LDDRAF*       SET TO .FALSE. IF LDCUM=.FALSE. OR KDTOP<KCTOP   

!    UPDATED PARAMETERS (REAL):

!    *PMFU*         MASSFLUX IN UPDRAFTS                          KG/(M2*S)
!    *PMFD*         MASSFLUX IN DOWNDRAFTS                        KG/(M2*S)
!    *PMFUS*        FLUX OF DRY STATIC ENERGY IN UPDRAFTS          J/(M2*S)
!    *PMFDS*        FLUX OF DRY STATIC ENERGY IN DOWNDRAFTS        J/(M2*S)
!    *PMFUQ*        FLUX OF SPEC. HUMIDITY IN UPDRAFTS            KG/(M2*S)
!    *PMFDQ*        FLUX OF SPEC. HUMIDITY IN DOWNDRAFTS          KG/(M2*S)
!    *PMFUL*        FLUX OF LIQUID WATER IN UPDRAFTS              KG/(M2*S)
!    *PLUDE*        DETRAINED LIQUID WATER                        KG/(M3*S)
!    *PDMFUP*       FLUX DIFFERENCE OF PRECIP. IN UPDRAFTS        KG/(M2*S)
!    *PDMFDP*       FLUX DIFFERENCE OF PRECIP. IN DOWNDRAFTS      KG/(M2*S)
!    *PMFDDE_RATE*  DOWNDRAFT DETRAINMENT RATE                    KG/(M2*S)

!    OUTPUT PARAMETERS (REAL):

!    *PDPMEL*       CHANGE IN PRECIP.-FLUXES DUE TO MELTING       KG/(M2*S)
!    *PLGLAC*       FLUX OF FROZEN CLOUD WATER IN UPDRAFTS        KG/(M2*S)
!    *PMFLXR*       CONVECTIVE RAIN FLUX                          KG/(M2*S)
!    *PMFLXS*       CONVECTIVE SNOW FLUX                          KG/(M2*S)
!    *PRAIN*        TOTAL PRECIP. PRODUCED IN CONV. UPDRAFTS      KG/(M2*S)
!                   (NO EVAPORATION IN DOWNDRAFTS)

!          EXTERNALS
!          ---------
!          NONE

!          MODIFICATIONS
!          -------------
!             99-06-14 : Optimisation        D.SALMOND
!             03-08-28 : Clean up LINPHYS    P.BECHTOLD
!                        Bugs in Evapor.
!             05-02-11 : Extend DDMflux to   P.BECHTOLD
!                        surface if zero
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!----------------------------------------------------------------------

USE MESSY_CONVECT_ECMWF_PARAM   ,ONLY : RG, RCPD, RLVTT, RLSTT, RLMLT, RTT,      &
                                        LPHYLIN, RLPTRC, RLPAL1,                 &
                                        R2ES, R3LES, R3IES, R4LES,               &
                                        R4IES, R5LES, R5IES, R5ALVCP, R5ALSCP,   &
                                        RALVDCP, RALSDCP, RTWAT, RTICE, RTICECU, &
                                        RTWAT_RTICE_R, RTWAT_RTICECU_R,          &
                                        RCUCOV, RCPECONS, RTAUMEL, RHEBC,        &!from YOECUMF
                                        LHOOK, DR_HOOK,                          &!from YOMHOOK
                                        FOELHMCU, FOEEWMCU, FOEALFCU              !functions
IMPLICIT NONE

INTEGER,INTENT(IN)    :: KLON 
INTEGER,INTENT(IN)    :: KLEV 
INTEGER,INTENT(IN)    :: KIDIA 
INTEGER,INTENT(IN)    :: KFDIA 
INTEGER               :: KTDIA ! Argument NOT used
REAL(dp)   ,INTENT(IN)    :: PTSPHY 
REAL(dp)   ,INTENT(IN)    :: PTEN(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PQEN(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PQSEN(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PTENH(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PQENH(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) 
REAL(dp)   ,INTENT(IN)    :: PAP(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PGEOH(KLON,KLEV+1) 
LOGICAL                          :: LDLAND(KLON) ! Argument NOT used
LOGICAL           ,INTENT(IN)    :: LDCUM(KLON) 
INTEGER,INTENT(IN)    :: KCBOT(KLON) 
INTEGER,INTENT(IN)    :: KCTOP(KLON) 
INTEGER,INTENT(IN)    :: KDTOP(KLON) 
INTEGER,INTENT(OUT)   :: KTOPM2 
INTEGER,INTENT(INOUT) :: KTYPE(KLON) 
LOGICAL    ,INTENT(INOUT) :: LDDRAF(KLON) 
REAL(dp)   ,INTENT(INOUT) :: PMFU(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PMFD(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PMFUS(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PMFDS(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PMFUQ(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PMFDQ(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PMFUL(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PLUDE(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PDMFUP(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PDMFDP(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PDPMEL(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PLGLAC(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PMFLXR(KLON,KLEV+1) 
REAL(dp)   ,INTENT(INOUT) :: PMFLXS(KLON,KLEV+1) 
REAL(dp)   ,INTENT(INOUT) :: PRAIN(KLON) 
REAL(dp)   ,INTENT(INOUT) :: PMFDDE_RATE(KLON,KLEV)

INTEGER :: IK, IKB, JK, JL
INTEGER :: IDBAS(KLON)
LOGICAL :: LLDDRAF
REAL(dp) :: ZALFAW, ZCONS1, ZCONS1A, ZCONS2, &
 & ZDENOM, ZDRFL, ZDRFL1, ZFAC, ZFOEEWI, &
 & ZFOEEWL, ZOEALFA, ZOEEWM, ZOELHM, ZPDR, ZPDS, &
 & ZRFL, ZRFLN, ZRMIN, ZRNEW, ZSNMLT, ZTARG, &
 & ZTMST, ZZP  
REAL(dp) :: ZMFS(KLON), ZMFMAX
REAL(dp) :: ZHOOK_HANDLE

!--------------------------------------------------------------------
!*             SPECIFY CONSTANTS

IF (LHOOK) CALL DR_HOOK('CUFLXN',0,ZHOOK_HANDLE)
ZTMST=PTSPHY
ZCONS1A=RCPD/(RLMLT*RG*RTAUMEL)
ZCONS2=1.0_dp/(RG*ZTMST)

!*    1.0          DETERMINE FINAL CONVECTIVE FLUXES
!                  ---------------------------------

DO JL=KIDIA,KFDIA
  PRAIN(JL)=0.0_dp
  IF(.NOT.LDCUM(JL).OR.KDTOP(JL) < KCTOP(JL)) LDDRAF(JL)=.FALSE.
  IF(.NOT.LDCUM(JL)) KTYPE(JL)=0
  IDBAS(JL)=KLEV
ENDDO
! TO GET IDENTICAL RESULTS FOR DIFFERENT NPROMA FORCE KTOPM2 TO 2
KTOPM2=2
DO JK=KTOPM2,KLEV
!DIR$ IVDEP
!OCL NOVREC
  IKB=MIN(JK+1,KLEV)
  DO JL=KIDIA,KFDIA
    PMFLXR(JL,JK)=0.0_dp
    PMFLXS(JL,JK)=0.0_dp
    PDPMEL(JL,JK)=0.0_dp
    IF(LDCUM(JL).AND.JK >= KCTOP(JL)-1) THEN
      PMFUS(JL,JK)=PMFUS(JL,JK)-PMFU(JL,JK)*(RCPD*PTENH(JL,JK)+PGEOH(JL,JK))
      PMFUQ(JL,JK)=PMFUQ(JL,JK)-PMFU(JL,JK)*PQENH(JL,JK)
      PLGLAC(JL,JK)=PMFU(JL,JK)*PLGLAC(JL,JK)
      LLDDRAF=LDDRAF(JL).AND.JK >= KDTOP(JL)
      IF(LLDDRAF) THEN
        PMFDS(JL,JK)=PMFDS(JL,JK)-PMFD(JL,JK)*&
         & (RCPD*PTENH(JL,JK)+PGEOH(JL,JK))  
        PMFDQ(JL,JK)=PMFDQ(JL,JK)-PMFD(JL,JK)*PQENH(JL,JK)
      ELSE
        PMFD(JL,JK)=0.0_dp
        PMFDS(JL,JK)=0.0_dp
        PMFDQ(JL,JK)=0.0_dp
        PDMFDP(JL,JK-1)=0.0_dp
      ENDIF
      IF(LLDDRAF.AND.PMFD(JL,JK)<0._dp .AND. PMFD(JL,IKB)==0._dp) THEN
        IDBAS(JL)=JK
      ENDIF
    ELSE
      PMFU(JL,JK)=0.0_dp
      PMFD(JL,JK)=0.0_dp
      PMFUS(JL,JK)=0.0_dp
      PMFDS(JL,JK)=0.0_dp
      PMFUQ(JL,JK)=0.0_dp
      PMFDQ(JL,JK)=0.0_dp
      PMFUL(JL,JK)=0.0_dp
      PLGLAC(JL,JK)=0.0_dp
      PDMFUP(JL,JK-1)=0.0_dp
      PDMFDP(JL,JK-1)=0.0_dp
      PLUDE(JL,JK-1)=0.0_dp
    ENDIF
  ENDDO
ENDDO
PMFLXR(:,KLEV+1)=0.0_dp
PMFLXS(:,KLEV+1)=0.0_dp

!*    1.5          SCALE FLUXES BELOW CLOUD BASE
!                  LINEAR DCREASE
!                  -----------------------------

!DIR$ IVDEP
!OCL NOVREC
DO JL=KIDIA,KFDIA
  IF(LDCUM(JL)) THEN
    IKB=KCBOT(JL)
    IK=IKB+1
    ZZP=((PAPH(JL,KLEV+1)-PAPH(JL,IK))/(PAPH(JL,KLEV+1)-PAPH(JL,IKB)))
    IF(KTYPE(JL) == 3) ZZP=ZZP*ZZP
    PMFU(JL,IK)=PMFU(JL,IKB)*ZZP
    IF (LPHYLIN) THEN
      ZTARG=PTENH(JL,IKB)
      ZOEALFA=0.5_dp*(TANH(RLPAL1*(ZTARG-RLPTRC))+1.0_dp)
      ZOELHM=ZOEALFA*RLVTT+(1.0_dp-ZOEALFA)*RLSTT
      PMFUS(JL,IK)=(PMFUS(JL,IKB)-ZOELHM*PMFUL(JL,IKB))*ZZP
    ELSE
      PMFUS(JL,IK)=(PMFUS(JL,IKB)-FOELHMCU(PTENH(JL,IKB))*PMFUL(JL,IKB))*ZZP
    ENDIF
    PMFUQ(JL,IK)=(PMFUQ(JL,IKB)+PMFUL(JL,IKB))*ZZP
    PMFUL(JL,IK)=0.0_dp
  ENDIF
ENDDO
DO JK=KTOPM2,KLEV
!DIR$ IVDEP
!OCL NOVREC
  DO JL=KIDIA,KFDIA
    IF(LDCUM(JL).AND.JK > KCBOT(JL)+1) THEN
      IKB=KCBOT(JL)+1
      ZZP=((PAPH(JL,KLEV+1)-PAPH(JL,JK))/(PAPH(JL,KLEV+1)-PAPH(JL,IKB)))
      PMFU(JL,JK)=PMFU(JL,IKB)*ZZP
      PMFUS(JL,JK)=PMFUS(JL,IKB)*ZZP
      PMFUQ(JL,JK)=PMFUQ(JL,IKB)*ZZP
      PMFUL(JL,JK)=0.0_dp
    ENDIF
    IK=IDBAS(JL)
    IF(LDDRAF(JL).AND.JK>IK.AND.IK<KLEV) THEN 
      ZZP=((PAPH(JL,KLEV+1)-PAPH(JL,JK))/(PAPH(JL,KLEV+1)-PAPH(JL,IK)))
    ! IF(KTYPE(JL) == 3) &
      ZZP=ZZP*ZZP
      PMFD(JL,JK)=PMFD(JL,IK)*ZZP
      PMFDS(JL,JK)=PMFDS(JL,IK)*ZZP
      PMFDQ(JL,JK)=PMFDQ(JL,IK)*ZZP
      PMFDDE_RATE(JL,JK)=-(PMFD(JL,JK-1)-PMFD(JL,JK))
    ENDIF
  ENDDO
ENDDO


!*    2.            CALCULATE RAIN/SNOW FALL RATES
!*                  CALCULATE MELTING OF SNOW
!*                  CALCULATE EVAPORATION OF PRECIP
!                   -------------------------------

DO JK=KTOPM2,KLEV
  DO JL=KIDIA,KFDIA
    IF(LDCUM(JL).AND.JK >= KCTOP(JL)-1) THEN
      PRAIN(JL)=PRAIN(JL)+PDMFUP(JL,JK)
      IF(PMFLXS(JL,JK) > 0.0_dp.AND.PTEN(JL,JK) > RTT) THEN
        ZCONS1=ZCONS1A*(1.0_dp+0.5_dp*(PTEN(JL,JK)-RTT))
        ZFAC=ZCONS1*(PAPH(JL,JK+1)-PAPH(JL,JK))
        ZSNMLT=MIN(PMFLXS(JL,JK),ZFAC*(PTEN(JL,JK)-RTT))
        PDPMEL(JL,JK)=ZSNMLT
        IF (LPHYLIN) THEN
          ZTARG=PTEN(JL,JK)-ZSNMLT/ZFAC
          ZOEALFA=0.5_dp*(TANH(RLPAL1*(ZTARG-RLPTRC))+1.0_dp)
          ZFOEEWL=R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
          ZFOEEWI=R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
          ZOEEWM=ZOEALFA*ZFOEEWL + (1.0_dp-ZOEALFA)*ZFOEEWI
          PQSEN(JL,JK)=ZOEEWM/PAP(JL,JK)
        ELSE
          PQSEN(JL,JK)=FOEEWMCU(PTEN(JL,JK)-ZSNMLT/ZFAC)/PAP(JL,JK)
        ENDIF
      ENDIF
      IF (LPHYLIN) THEN
        ZALFAW=0.5_dp*(TANH(RLPAL1*(PTEN(JL,JK)-RLPTRC))+1.0_dp)
      ELSE
        ZALFAW=FOEALFCU(PTEN(JL,JK))
      ENDIF
      PMFLXR(JL,JK+1)=PMFLXR(JL,JK)+ZALFAW*&
       & (PDMFUP(JL,JK)+PDMFDP(JL,JK))+PDPMEL(JL,JK)  
      PMFLXS(JL,JK+1)=PMFLXS(JL,JK)+(1.0_dp-ZALFAW)*&
       & (PDMFUP(JL,JK)+PDMFDP(JL,JK))-PDPMEL(JL,JK)  
      IF(PMFLXR(JL,JK+1)+PMFLXS(JL,JK+1) < 0.0_dp) THEN
        PDMFDP(JL,JK)=-(PMFLXR(JL,JK)+PMFLXS(JL,JK)+PDMFUP(JL,JK))
        PMFLXR(JL,JK+1)=0.0_dp
        PMFLXS(JL,JK+1)=0.0_dp
        PDPMEL(JL,JK)  =0.0_dp
      ELSEIF(PMFLXR(JL,JK+1) < 0.0_dp ) THEN
        PMFLXS(JL,JK+1)=PMFLXS(JL,JK+1)+PMFLXR(JL,JK+1)
        PMFLXR(JL,JK+1)=0.0_dp
      ELSEIF(PMFLXS(JL,JK+1) < 0.0_dp ) THEN
        PMFLXR(JL,JK+1)=PMFLXR(JL,JK+1)+PMFLXS(JL,JK+1)
        PMFLXS(JL,JK+1)=0.0_dp
      ENDIF
    ENDIF
  ENDDO
ENDDO

! Reminder for conservation: 
!    pdmfup(jl,jk)+pdmfdp(jl,jk)=pmflxr(jl,jk+1)+pmflxs(jl,jk+1)-pmflxr(jl,jk)-pmflxs(jl,jk)

DO JK=KTOPM2,KLEV
  DO JL=KIDIA,KFDIA
    IF(LDCUM(JL).AND.JK >= KCBOT(JL)) then
      ZRFL=PMFLXR(JL,JK)+PMFLXS(JL,JK)
      IF(ZRFL > 1.E-20_dp) THEN
        ZDRFL1=RCPECONS*MAX(0.0_dp,PQSEN(JL,JK)-PQEN(JL,JK))*RCUCOV*&
         & (SQRT(PAPH(JL,JK)/PAPH(JL,KLEV+1))/5.09E-3_dp*&
         & ZRFL/RCUCOV)**0.5777_dp*&
         & (PAPH(JL,JK+1)-PAPH(JL,JK))  
        ZRNEW=ZRFL-ZDRFL1
        ZRMIN=ZRFL-RCUCOV*MAX(0.0_dp,RHEBC*PQSEN(JL,JK)-PQEN(JL,JK))&
         & *ZCONS2*(PAPH(JL,JK+1)-PAPH(JL,JK))  
        ZRNEW=MAX(ZRNEW,ZRMIN)
        ZRFLN=MAX(ZRNEW,0.0_dp)
        ZDRFL=MIN(0.0_dp,ZRFLN-ZRFL)
        IF (LPHYLIN) THEN
          ZALFAW=0.5_dp*(TANH(RLPAL1*(PTEN(JL,JK)-RLPTRC))+1.0_dp)
        ELSE
          ZALFAW=FOEALFCU(PTEN(JL,JK))
        ENDIF
        ZPDR=ZALFAW*PDMFDP(JL,JK)
        ZPDS=(1.0_dp-ZALFAW)*PDMFDP(JL,JK)
        ZDENOM=1.0_dp/MAX(1.E-20_dp,PMFLXR(JL,JK)+PMFLXS(JL,JK))
        PMFLXR(JL,JK+1)=PMFLXR(JL,JK)+ZPDR &
         & +PDPMEL(JL,JK)+ZDRFL*PMFLXR(JL,JK)*ZDENOM  
        PMFLXS(JL,JK+1)=PMFLXS(JL,JK)+ZPDS &
         & -PDPMEL(JL,JK)+ZDRFL*PMFLXS(JL,JK)*ZDENOM  
        PDMFUP(JL,JK)=PDMFUP(JL,JK)+ZDRFL
        IF(PMFLXR(JL,JK+1)+PMFLXS(JL,JK+1) < 0.0_dp) THEN
          PDMFUP(JL,JK)=PDMFUP(JL,JK)-(PMFLXR(JL,JK+1)+PMFLXS(JL,JK+1))
          PMFLXR(JL,JK+1)=0.0_dp
          PMFLXS(JL,JK+1)=0.0_dp
          PDPMEL(JL,JK)  =0.0_dp
        ELSEIF(PMFLXR(JL,JK+1) < 0.0_dp ) THEN
          PMFLXS(JL,JK+1)=PMFLXS(JL,JK+1)+PMFLXR(JL,JK+1)
          PMFLXR(JL,JK+1)=0.0_dp
        ELSEIF(PMFLXS(JL,JK+1) < 0.0_dp ) THEN
          PMFLXR(JL,JK+1)=PMFLXR(JL,JK+1)+PMFLXS(JL,JK+1)
          PMFLXS(JL,JK+1)=0.0_dp
        ENDIF
      ELSE
        PMFLXR(JL,JK+1)=0.0_dp
        PMFLXS(JL,JK+1)=0.0_dp
        PDMFDP(JL,JK)  =0.0_dp
        PDPMEL(JL,JK)  =0.0_dp
      ENDIF
    ENDIF
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('CUFLXN',1,ZHOOK_HANDLE)
END SUBROUTINE CUFLXN

!==============================================================================

!------------------------------------------------------------------------------
END MODULE MESSY_CONVECT_ECMWF
