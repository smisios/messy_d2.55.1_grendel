SUBROUTINE SATUR ( KIDIA , KFDIA , KLON  , KTDIA , KLEV,&
 & PAPRSF, PT    , PQSAT , KFLAG)  

!***

! **   *SATUR* -  COMPUTES SPECIFIC HUMIDITY AT SATURATION

!       J.F. MAHFOUF       E.C.M.W.F.     15/05/96

!       Modified J. HAGUE          13/01/03 MASS Vector Functions       

!       PURPOSE.
!       --------

!       SPECIFIC HUMIDITY AT SATURATION IS USED BY THE
!       DIAGNOSTIC CLOUD SCHEME TO COMPUTE RELATIVE HUMIDITY
!       AND LIQUID WATER CONTENT  

!       INTERFACE
!       ---------

!       THIS ROUTINE IS CALLED FROM *CALLPAR*.

!       PARAMETER     DESCRIPTION                                 UNITS
!       ---------     -----------                                 -----
!       INPUT PARAMETERS (INTEGER):

!      *KIDIA*        START POINT
!      *KFDIA*        END POINT
!      *KLON*         NUMBER OF GRID POINTS PER PACKET
!      *KTDIA*        START OF THE VERTICAL LOOP
!      *KLEV*         NUMBER OF LEVELS

!       INPUT PARAMETERS (REAL):

!      *PAPRSF*        PRESSURE ON FULL LEVELS                      PA
!      *PT*            TEMPERATURE AT T-DT                          K

!       INPUT PARAMETERS (INTEGER):

!      *KFLAG*         FLAG TO DETECT CALL FROM

!                      CONVECTION  KFLAG=1
!                      OTHER       KFLAG=2

!       OUTPUT PARAMETER (REAL):

!      *PQSAT*         SATURATION SPECIFIC HUMIDITY                 KG/KG

!-------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK_IFS   ,ONLY : LHOOK,   DR_HOOK

USE YOMCST   , ONLY : RETV     ,RLVTT    ,RLSTT    ,RTT
USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
 & R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
 & RALVDCP  ,RALSDCP  ,RTWAT    ,RTICE    ,RTICECU  ,&
 & RTWAT_RTICE_R      ,RTWAT_RTICECU_R  
USE YOEPHLI  , ONLY : LPHYLIN  
USE YOMJFH   , ONLY : N_VMASS

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTDIA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRSF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQSAT(KLON,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLAG 
INTEGER(KIND=JPIM) :: JK, JL, JLEN

REAL(KIND=JPRB) :: ZCOR, ZEW, ZFOEEW, ZQMAX, ZQS, ZTARG
REAL(KIND=JPRB) :: ZALFA, ZFOEEWL, ZFOEEWI
REAL(KIND=JPRB) :: Z_EXPARG1(KIDIA:KFDIA)
REAL(KIND=JPRB) :: Z_EXPARG2(KIDIA:KFDIA)
REAL(KIND=JPRB) :: Z_EXPOUT1(KIDIA:KFDIA)
REAL(KIND=JPRB) :: Z_EXPOUT2(KIDIA:KFDIA)
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!DIR$ VFUNCTION EXPHF

#include "fcttre.h"

!----------------------------------------------------------------------

!*    1.           DEFINE CONSTANTS
!                  ----------------

IF (LHOOK) CALL DR_HOOK('SATUR',0,ZHOOK_HANDLE)
ZQMAX=0.5_JPRB

!     *
!----------------------------------------------------------------------

!     *    2.           CALCULATE SATURATION SPECIFIC HUMIDITY
!                       --------------------------------------

IF (LPHYLIN) THEN
  DO JK=KTDIA,KLEV
    DO JL=KIDIA, KFDIA
      ZTARG = PT(JL,JK)
      ZALFA = FOEALFA(ZTARG)

      ZFOEEWL = R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
      ZFOEEWI = R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
      ZFOEEW = ZALFA*ZFOEEWL+(1.0_JPRB-ZALFA)*ZFOEEWI

      ZQS    = ZFOEEW/PAPRSF(JL,JK)
      IF (ZQS > ZQMAX) THEN
        ZQS=ZQMAX
      ENDIF
      ZCOR = 1.0_JPRB/(1.0_JPRB-RETV*ZQS)
      PQSAT(JL,JK)=ZQS*ZCOR
    ENDDO
  ENDDO
ELSE

  IF(N_VMASS <= 0) THEN ! Not using Vector MASS

    DO JK=KTDIA,KLEV
      DO JL=KIDIA, KFDIA
        IF(KFLAG == 1) THEN
          ZEW  = FOEEWMCU(PT(JL,JK))
        ELSE
          ZEW  = FOEEWM(PT(JL,JK))
        ENDIF
        ZQS  = ZEW/PAPRSF(JL,JK)
        ZQS  = MIN(ZQMAX,ZQS)
        ZCOR = 1.0_JPRB/(1.0_JPRB-RETV*ZQS)
        PQSAT(JL,JK)=ZQS*ZCOR
      ENDDO
    ENDDO

  ELSE ! Using Vector MASS

    JLEN=KFDIA-KIDIA+1

    DO JK=KTDIA,KLEV
      DO JL=KIDIA, KFDIA
        Z_EXPARG1(JL)=FOELES_V(PT(JL,JK))
        Z_EXPARG2(JL)=FOEIES_V(PT(JL,JK))
      ENDDO
      CALL VEXP(Z_EXPOUT1,Z_EXPARG1,JLEN)
      CALL VEXP(Z_EXPOUT2,Z_EXPARG2,JLEN)
      DO JL=KIDIA, KFDIA
        IF(KFLAG == 1) THEN
          ZEW  = FOEEWMCU_V( PT(JL,JK),Z_EXPOUT1(JL),Z_EXPOUT2(JL) )
        ELSE
          ZEW  = FOEEWM_V  ( PT(JL,JK),Z_EXPOUT1(JL),Z_EXPOUT2(JL) )
        ENDIF
!       ZQS  = ZEW/PAPRSF(JL,JK)
!       ZQS  = MIN(ZQMAX,ZQS)
!!      ZCOR = _ONE_/(_ONE_-RETV*ZQS)
!       PQSAT(JL,JK)=ZQS/(_ONE_-RETV*ZQS)
        ZQS  = MIN(ZQMAX*PAPRSF(JL,JK),ZEW)
        PQSAT(JL,JK)=ZQS/(PAPRSF(JL,JK)-RETV*ZQS)
      ENDDO
    ENDDO
  
  ENDIF

ENDIF

IF (LHOOK) CALL DR_HOOK('SATUR',1,ZHOOK_HANDLE)
END SUBROUTINE SATUR
