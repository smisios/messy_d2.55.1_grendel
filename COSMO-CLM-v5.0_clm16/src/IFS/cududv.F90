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

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK_IFS   ,ONLY : LHOOK,   DR_HOOK

USE YOMCST   , ONLY : RG
USE YOECUMF  , ONLY : RMFSOLUV

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM)               :: KTDIA ! Argument NOT used
INTEGER(KIND=JPIM),INTENT(IN)    :: KTOPM2 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTYPE(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCBOT(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCTOP(KLON) 
LOGICAL           ,INTENT(IN)    :: LDCUM(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUEN(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVEN(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFD(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUD(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVD(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENV(KLON,KLEV) 
REAL(KIND=JPRB) :: ZUEN(KLON,KLEV),     ZVEN(KLON,KLEV),&
 & ZMFUU(KLON,KLEV),    ZMFDU(KLON,KLEV),&
 & ZMFUV(KLON,KLEV),    ZMFDV(KLON,KLEV)

INTEGER(KIND=JPIM) :: IK, IKB, JK, JL

REAL(KIND=JPRB) :: ZZP, ZIMP, ZTSPHY
!       ALLOCATABLE ARAYS
REAL(KIND=JPRB),   DIMENSION(:,:), ALLOCATABLE :: ZDUDT, ZDVDT, ZDP
REAL(KIND=JPRB),   DIMENSION(:,:), ALLOCATABLE :: ZB,  ZR1,  ZR2
LOGICAL,DIMENSION(:,:),   ALLOCATABLE :: LLCUMBAS
REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "cubidiag.intfb.h"
!----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CUDUDV',0,ZHOOK_HANDLE)
ZIMP=1.0_JPRB-RMFSOLUV
ZTSPHY=1.0_JPRB/PTSPHY

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

!----------------------------------------------------------------------


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
  IF(RMFSOLUV==0.0_JPRB) THEN
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
  ENDIF

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

  IF ( RMFSOLUV==0.0_JPRB ) THEN

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
      ZB(:,:)=1.0_JPRB
      ZMFUU(:,:)=0.0_JPRB
     
     ! Fill vectors A, B and RHS 
     
      DO JK=KTOPM2,KLEV
         IK=JK+1
         DO JL=KIDIA,KFDIA
           LLCUMBAS(JL,JK)=LDCUM(JL).AND.JK>=KCTOP(JL)-1 
           IF(LLCUMBAS(JL,JK)) THEN
             ZZP=RMFSOLUV*ZDP(JL,JK)*PTSPHY
             ZMFUU(JL,JK)=-ZZP*(PMFU(JL,JK)+PMFD(JL,JK))
             ZDUDT(JL,JK) = ZDUDT(JL,JK)*PTSPHY+ZUEN(JL,JK)
             ZDVDT(JL,JK) = ZDVDT(JL,JK)*PTSPHY+ZVEN(JL,JK)
           ! ZDUDT(JL,JK) = (PTENU(JL,JK)+ZDUDT(JL,JK))*PTSPHY+ZUEN(JL,JK)
           ! ZDVDT(JL,JK) = (PTENV(JL,JK)+ZDVDT(JL,JK))*PTSPHY+ZVEN(JL,JK)
             IF(JK<KLEV) THEN
               ZB(JL,JK)=1.0_JPRB+ZZP*(PMFU(JL,IK)+PMFD(JL,IK))
             ELSE
               ZB(JL,JK)=1.0_JPRB
             ENDIF
           ENDIF
         ENDDO
      ENDDO
     
      CALL CUBIDIAG&
         &( KIDIA, KFDIA, KLON, KLEV, &
         &  KCTOP, LLCUMBAS, &
         &  ZMFUU,    ZB,    ZDUDT,   ZR1 )
     
      CALL CUBIDIAG&
         &( KIDIA, KFDIA, KLON, KLEV, &
         &  KCTOP, LLCUMBAS, &
         &  ZMFUU,    ZB,    ZDVDT,   ZR2 )
     
     ! Compute tendencies
     
      DO JK=KTOPM2,KLEV
         DO JL=KIDIA,KFDIA
           IF(LLCUMBAS(JL,JK)) THEN
             PTENU(JL,JK)=PTENU(JL,JK)+(ZR1(JL,JK)-ZUEN(JL,JK))*ZTSPHY
             PTENV(JL,JK)=PTENV(JL,JK)+(ZR2(JL,JK)-ZVEN(JL,JK))*ZTSPHY
           ! PTENU(JL,JK)=(ZR1(JL,JK)-ZUEN(JL,JK))*ZTSPHY
           ! PTENV(JL,JK)=(ZR2(JL,JK)-ZVEN(JL,JK))*ZTSPHY
           ENDIF
         ENDDO
      ENDDO
     
      DEALLOCATE(LLCUMBAS)
      DEALLOCATE(ZR2)
      DEALLOCATE(ZR1)
      DEALLOCATE(ZB)
     
   ENDIF
!----------------------------------------------------------------------

DEALLOCATE(ZDP)
DEALLOCATE(ZDVDT)
DEALLOCATE(ZDUDT)

IF (LHOOK) CALL DR_HOOK('CUDUDV',1,ZHOOK_HANDLE)
END SUBROUTINE CUDUDV


