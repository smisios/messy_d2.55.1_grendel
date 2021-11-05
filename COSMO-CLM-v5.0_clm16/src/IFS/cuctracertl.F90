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

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK_IFS   ,ONLY : LHOOK,   DR_HOOK

USE YOMCST   , ONLY : RG
USE YOECUMF  , ONLY : RMFSOLCT, RMFCMIN

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KTRAC 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM)               :: KTDIA ! Argument NOT used
INTEGER(KIND=JPIM),INTENT(IN)    :: KTYPE(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCTOP(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCBOT(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDPL(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDTOP(KLON) 
LOGICAL           ,INTENT(IN)    :: LDCUM(KLON) 
LOGICAL           ,INTENT(IN)    :: LDDRAF(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH5(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFU5(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFD5(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUDRATE5(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDDRATE5(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCEN5(KLON,KLEV,KTRAC) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENC5(KLON,KLEV,KTRAC) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFD(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUDRATE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDDRATE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCEN(KLON,KLEV,KTRAC) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENC(KLON,KLEV,KTRAC) 
! Trajectory arrays

! Tangent-linear perturbation arrays

INTEGER(KIND=JPIM) :: IK, IKB, JK, JL, JN, JITER, JIT 

REAL(KIND=JPRB) :: ZZP5, ZMFA5, ZERATE5, ZPOSI5, ZIMP, ZTSPHY
REAL(KIND=JPRB) :: ZZP , ZMFA , ZERATE , ZFACT

!     ALLOCATABLE ARAYS
REAL(KIND=JPRB), DIMENSION(:,:,:), ALLOCATABLE :: ZCEN5, ZCU5, ZCD5, ZTENC5, ZMFC5
REAL(KIND=JPRB), DIMENSION(:,:),   ALLOCATABLE :: ZDP5, ZB5,  ZR15
REAL(KIND=JPRB), DIMENSION(:,:,:), ALLOCATABLE :: ZCEN, ZCU, ZCD, ZTENC, ZMFC
REAL(KIND=JPRB), DIMENSION(:,:),   ALLOCATABLE :: ZDP, ZB,  ZR1
LOGICAL, DIMENSION(:,:),  ALLOCATABLE :: LLCUMASK, LLCUMBAS
REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "cubidiagtl.intfb.h"
!----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CUCTRACERTL',0,ZHOOK_HANDLE)
ZIMP=1.0_JPRB-RMFSOLCT
ZTSPHY=1._JPRB/PTSPHY

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
         ZTENC5(JL,JK,JN)=0._JPRB
         ZTENC(JL,JK,JN) =0._JPRB
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
        ! IF(ZERATE5<0._JPRB) THEN
        !   ZERATE5=0._JPRB
        !   ZERATE =0._JPRB
        ! ENDIF
          ZMFA5=1._JPRB/MAX(RMFCMIN,PMFU5(JL,JK))
          IF (PMFU5(JL,JK) > RMFCMIN) THEN
            ZMFA=-PMFU(JL,JK)/(PMFU5(JL,JK)*PMFU5(JL,JK))
          ELSE
            ZMFA=0._JPRB
          ENDIF

          IF ( JK>=KCTOP(JL) )  THEN
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
          !ZCD5(JL,JK,JN)=.5_JPRB*(ZCU5(JL,JK,JN)+PCEN5(JL,IK,JN))
          !ZCD(JL,JK,JN) =.5_JPRB*(ZCU(JL,JK,JN)+ZCEN(JL,IK,JN))
           ZCD5(JL,JK,JN)=.1_JPRB*ZCU5(JL,JK,JN)+.9_JPRB*PCEN5(JL,IK,JN)
           ZCD(JL,JK,JN) =.1_JPRB*ZCU(JL,JK,JN)+.9_JPRB*PCEN(JL,IK,JN)
         ELSEIF ( LDDRAF(JL).AND.JK>KDTOP(JL) ) THEN
           ZERATE5=-PMFD5(JL,JK)+PMFD5(JL,IK)+PDDRATE5(JL,JK)  
           ZERATE=-PMFD(JL,JK)+PMFD(JL,IK)+PDDRATE(JL,JK)
         ! IF(ZERATE5<0._JPRB) THEN
         !   ZERATE5=0._JPRB
         !   ZERATE =0._JPRB
         ! ENDIF
           ZMFA5=1._JPRB/MIN(-RMFCMIN,PMFD5(JL,JK))
           IF (PMFD5(JL,JK) < -RMFCMIN) THEN
             ZMFA=-PMFD(JL,JK)/(PMFD5(JL,JK)*PMFD5(JL,JK))
           ELSE
             ZMFA=0._JPRB
           ENDIF
           ZCD5(JL,JK,JN)=( PMFD5(JL,IK)*ZCD5(JL,IK,JN)-ZERATE5*PCEN5(JL,IK,JN) &
                       & + PDDRATE5(JL,JK)*ZCD5(JL,IK,JN) )*ZMFA5 
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
     IF( PCEN5(JL,JK,JN)+ZPOSI5*PTSPHY<0.0_JPRB ) THEN
        ZMFA5=1._JPRB/MIN(-RMFCMIN,PMFD5(JL,JK))
        IF (PMFD5(JL,JK) < -RMFCMIN) THEN
           ZMFA=-PMFD(JL,JK)/(PMFD5(JL,JK)*PMFD5(JL,JK))
        ELSE
           ZMFA=0._JPRB
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

ZMFC5(:,:,:)=0._JPRB


JITER=1
!IF(RMFSOLCT>0._JPRB) JITER=2
!----------------------------------------------------------------------
 
DO JIT=1,JITER

   IF(JIT==2) ZIMP=1._JPRB

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
        !    ZFACT=1._JPRB/(PAPH5(JL,KLEV+1)-PAPH5(JL,IKB))
        !    ZZP5=(PAPH5(JL,KLEV+1)-PAPH5(JL,JK))*ZFACT
        !    ZZP =((PAPH (JL,KLEV+1)-PAPH (JL,JK)) &
        !      & -  ZZP5*(PAPH (JL,KLEV+1)-PAPH (JL,IKB)))*ZFACT 
        !    IF(KTYPE(JL) == 3) THEN
        !      ZZP =2._JPRB*ZZP*ZZP5
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


 IF ( ZIMP==1._JPRB ) THEN


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
    ZB5(:,:)=1._JPRB

    ALLOCATE(ZB(KLON,KLEV))
    ALLOCATE(ZR1(KLON,KLEV))   
    ZB(:,:)=0._JPRB
 
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
             ZB5(JL,JK)=1._JPRB+ZZP5*(PMFU5(JL,IK)+PMFD5(JL,IK))
             ZB(JL,JK)=ZZP5*(PMFU (JL,IK)+PMFD (JL,IK)) &
                    & +ZZP *(PMFU5(JL,IK)+PMFD5(JL,IK))
           ELSE
             ZB5(JL,JK)=1._JPRB
             ZB(JL,JK)=0._JPRB
           ENDIF
         ENDIF
       ENDDO
    ENDDO
 
    CALL CUBIDIAGTL&
       &( KIDIA, KFDIA, KLON, KLEV, &
       &  KCTOP, LLCUMBAS, &
       &  ZMFC5(:,:,JN),  ZB5,   ZTENC5(:,:,JN),   ZR15, &
       &  ZMFC (:,:,JN),  ZB ,   ZTENC (:,:,JN),   ZR1  )
 
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
