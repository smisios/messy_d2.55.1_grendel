SUBROUTINE CUANCAPE2(KIDIA,  KFDIA,  KLON,  KTDIA,  KLEV,&
                   &PAP,    PAPH,   PT,    PQ,     PCAPE)

!***** CUANCAPE2 - COMPUTE APPROXIMATE CAPE USING THETAE AND THETAES

!     E. HOLM + P. BECHTOLD     E.C.M.W.F.     13/10/2005

!     PURPOSE 
!     -------
!                 ESTIMATE CAPE FIRST FOR A MIXED-LAYER PARCEL, THEN
!                 LOOP OVER SUBSEQUENT DEPARTURE LAYERS IN LOWEST 350 hPa
!                 Theta_e =Theta*exp[L q_v/(C_p T)] 
!                         = T*(P0/P)**(R_d/C_p) * exp[L q_v/(C_p T)]
!                 -> THIS WILL BE THE UPDRAUGHT PARCEL (CONSERVING ITS
!                 PROPERTIES)  (no entrainment)
!                 CAPE    = Int ( g dTheta_v/Theta_v dz ) = 
!                   aprox = Int ( g (Theta_e_up-Theta_e_sat)/Theta_e_sat ) dz
!                 WITH THIS FORMULATION THE ACTUAL CAPE IS OVERESTIMATED  BY
!                 ROUGHLY 20%. DEEP CONVECTION CAN BE CONSIDERED FOR CAPE
!                 VALUES ABOVE 200-500 J/KG            


!     PARAMETER     DESCRIPTION                              UNITS
!     ---------     -----------                              -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KTDIA*        START OF THE VERTICAL LOOP
!    *KLEV*         NUMBER OF LEVELS

!    INPUT PARAMETERS (REAL):

!    *PAP*          PRESSURE ON FULL LEVELS                    PA
!    *PAPH*         PRESSURE ON HALF LEVELS                    PA
!    *PT*           TEMPERATURE ON FULL LEVELS                 K   
!    *PQ*           SPECIFIC HUMIDITY ON FULL LEVELS          KG/KG

!    OUTPUT PARAMETERS (REAL):

!    *CAPE*         CONVECTIVE AVAILABLE POT. ENERGY          J/KG

!-------------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK_IFS   ,ONLY : LHOOK,   DR_HOOK
USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
 & R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
 & RALVDCP  ,RALSDCP  ,RTWAT    ,RTICE    ,RTICECU  ,&
 & RTWAT_RTICE_R      ,RTWAT_RTICECU_R

USE YOMCST   , ONLY : RETV     ,RLVTT    ,RLSTT    ,RTT,&
                     &RD       ,RKAPPA   ,RATM
USE YOECUMF  , ONLY : NJKT1, NJKT2

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM)               :: KTDIA ! Argument NOT used
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCAPE(KLON)

INTEGER(KIND=JPIM) :: JL, JK, JKK

REAL(KIND=JPRB), DIMENSION(KLON) :: ZPMIX, ZTMIX, ZTHMIX, ZQMIX, ZTHETEU
REAL(KIND=JPRB)                  :: ZCAPE(KLON,KLEV)
REAL(KIND=JPRB) :: ZDP, ZTHETES, ZQS, ZDZ, ZTEMP, ZRPAP

REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "fcttre.h"
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CUANCAPE2',0,ZHOOK_HANDLE)

PCAPE(:)=0.0_JPRB

DO JKK=KLEV-1,NJKT1,-1

  DO JL=KIDIA,KFDIA
    ZCAPE(JL,JKK)=0.0_JPRB
    IF (PAPH(JL,KLEV+1)-PAPH(JL,JKK-1)<60.E2_JPRB) THEN
      ZTMIX(JL)=0.0_JPRB
      ZTHMIX(JL)=0.0_JPRB
      ZQMIX(JL)=0.0_JPRB
      ZPMIX(JL)=0.0_JPRB
      DO JK=JKK+1,JKK-1,-1
        IF(ZPMIX(JL)<30.E2_JPRB) THEN
          ZDP=PAPH(JL,JK+1)-PAPH(JL,JK)
          ZPMIX(JL)=ZPMIX(JL)+ZDP
          ZTHMIX(JL)=ZTHMIX(JL)+PT(JL,JK)*ZDP*(RATM/PAP(JL,JK))**RKAPPA
          ZQMIX(JL)=ZQMIX(JL)+PQ(JL,JK)*ZDP
        ENDIF
      ENDDO
      ZDP=1.0_JPRB/ZPMIX(JL)
      ZQMIX(JL)=ZQMIX(JL)*ZDP
      ZPMIX(JL)=PAPH(JL,JKK+2)-0.5_JPRB*ZPMIX(JL)
      ZTHMIX(JL)=ZTHMIX(JL)*ZDP
      ZTMIX(JL)=ZTHMIX(JL)*(ZPMIX(JL)/RATM)**RKAPPA
    ELSE
      ZQMIX(JL)=PQ(JL,JKK)
      ZPMIX(JL)=PAP(JL,JKK)
      ZTMIX(JL)=PT(JL,JKK)
      ZTHMIX(JL)=PT(JL,JKK)*(RATM/ZPMIX(JL))**RKAPPA
    ENDIF
    ZTHETEU(JL)=ZTHMIX(JL)*&
               &EXP( FOELDCP(ZTMIX(JL))*ZQMIX(JL)/ZTMIX(JL) )
  ENDDO

  DO JK=JKK,NJKT2,-1
     DO JL=KIDIA,KFDIA
       IF(PAP(JL,JK)>80.E2_JPRB.AND. &
             & (PAPH(JL,KLEV+1)-PAPH(JL,JKK))<350.E2_JPRB) THEN
          ZRPAP=1.0_JPRB/PAP(JL,JK)
          ZQS = FOEEWM(PT(JL,JK))*ZRPAP
          ZQS = MAX(1.E-8_JPRB,ZQS)
          ZQS = ZQS/(1.0_JPRB-RETV*ZQS) ! small correction
          ZTHETES=PT(JL,JK)*(RATM*ZRPAP)**RKAPPA*&
                 &EXP( FOELDCP(PT(JL,JK))*ZQS/PT(JL,JK) ) 
          ZTEMP=ZTHETEU(JL)/ZTHETES-1.0_JPRB
          IF(ZTEMP > 0.0_JPRB)THEN
            ZDZ=(PAPH(JL,JK+1)-PAPH(JL,JK))*ZRPAP*RD*PT(JL,JK)*&
             &(1.0_JPRB+RETV*PQ(JL,JK))
            ZCAPE(JL,JKK)=ZCAPE(JL,JKK)+ZTEMP*ZDZ
          ENDIF
       ENDIF
     ENDDO
   ENDDO
   
ENDDO

      ! chose maximum CAPE value
DO JL=KIDIA,KFDIA
  PCAPE(JL)=MAXVAL(ZCAPE(JL,NJKT1:KLEV-1))
ENDDO

!-------------------------------------------------------------------------------
  
IF (LHOOK) CALL DR_HOOK('CUANCAPE2',1,ZHOOK_HANDLE)
END SUBROUTINE CUANCAPE2
