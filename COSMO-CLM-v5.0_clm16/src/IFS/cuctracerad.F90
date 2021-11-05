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
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAPH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PMFU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PMFD(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUDRATE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDDRATE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCEN(KLON,KLEV,KTRAC) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENC(KLON,KLEV,KTRAC) 
! Trajectory arrays

! Tangent-linear perturbation arrays

INTEGER(KIND=JPIM) :: IK, JK, JL, JN, JIT

REAL(KIND=JPRB) :: ZIMP, ZTSPHY
REAL(KIND=JPRB) :: ZZP , ZMFA , ZERATE , ZFACT, ZPOSI5, ZMFA5

!     ALLOCATABLE ARAYS
REAL(KIND=JPRB), DIMENSION(:,:),   ALLOCATABLE :: ZMFC15, ZDP5, ZR15, ZZP35, ZCD15
REAL(KIND=JPRB), DIMENSION(:,:,:), ALLOCATABLE :: ZCEN5, ZCU5, ZCD5, ZTENC5
REAL(KIND=JPRB), DIMENSION(:,:,:), ALLOCATABLE :: ZMFA15, ZMFA25
REAL(KIND=JPRB), DIMENSION(:,:,:), ALLOCATABLE :: ZERATE15, ZERATE25, ZERATE35, ZERATE45
REAL(KIND=JPRB), DIMENSION(:,:,:,:), ALLOCATABLE :: ZMFA35, ZB5, ZMFC25
REAL(KIND=JPRB), DIMENSION(:,:,:,:), ALLOCATABLE :: ZMFC5, ZMFC45
REAL(KIND=JPRB), DIMENSION(:,:,:,:), ALLOCATABLE :: ZTENC15

REAL(KIND=JPRB), DIMENSION(:,:,:), ALLOCATABLE :: ZCEN, ZCU, ZCD, ZTENC, ZMFC
REAL(KIND=JPRB), DIMENSION(:,:),   ALLOCATABLE :: ZDP, ZB,  ZR1
LOGICAL, DIMENSION(:,:),  ALLOCATABLE :: LLCUMASK, LLCUMBAS
REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "cubidiag.intfb.h"
#include "cubidiagad.intfb.h"
!----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CUCTRACERAD',0,ZHOOK_HANDLE)
ZIMP=1.0_JPRB-RMFSOLCT
ZTSPHY=1.0_JPRB/PTSPHY
JIT=1


! Trajectory arrays

ALLOCATE(ZCEN5(KLON,KLEV,KTRAC))     !Half-level environmental values
ALLOCATE(ZCU5(KLON,KLEV,KTRAC))      !Updraft values
ALLOCATE(ZCD5(KLON,KLEV,KTRAC))      !Downdraft values
ALLOCATE(ZCD15(KLON,KTRAC))          !Downdraft values at KLEV
ALLOCATE(ZTENC5(KLON,KLEV,KTRAC))    !Tendency
ALLOCATE(ZMFA15(KLON,KLEV,KTRAC))    !Fluxes  
ALLOCATE(ZMFA25(KLON,KLEV,KTRAC))    !Fluxes  
ALLOCATE(ZMFA35(KLON,KLEV,KTRAC,1))  !Fluxes
ALLOCATE(ZERATE15(KLON,KLEV,KTRAC))  !Entrainment rates
ALLOCATE(ZERATE25(KLON,KLEV,KTRAC))  !Entrainment rates
ALLOCATE(ZERATE35(KLON,KLEV,KTRAC))  !Entrainment rates
ALLOCATE(ZERATE45(KLON,KLEV,KTRAC))  !Entrainment rates
ALLOCATE(ZMFC15(KLON,KLEV))          !Fluxes
ALLOCATE(ZMFC25(KLON,KLEV,KTRAC,1))  !Fluxes
ALLOCATE(ZMFC5(KLON,KLEV,KTRAC,1))   !Fluxes
ALLOCATE(ZMFC45(KLON,KLEV,KTRAC,1))  !Fluxes
ALLOCATE(ZDP5(KLON,KLEV))            !Pressure difference
ALLOCATE(ZZP35(KLON,KLEV))           !Pressure difference
ALLOCATE(LLCUMASK(KLON,KLEV))        !Mask for convection

IF (RMFSOLCT > 0._JPRB) THEN
  ALLOCATE(ZB5(KLON,KLEV,KTRAC,1))
  ALLOCATE(ZR15(KLON,KLEV))   
  ALLOCATE(ZTENC15(KLON,KLEV,KTRAC,1))
  ALLOCATE(LLCUMBAS(KLON,KLEV))
  LLCUMBAS=.FALSE.
  ZB5=1.0_JPRB
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
         ZTENC5(JL,JK,JN)=0._JPRB
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
        ! IF(ZERATE15(JL,JK,JN)<0._JPRB) ZERATE15(JL,JK,JN)=0._JPRB
          ZMFA15(JL,JK,JN)=1._JPRB/MAX(RMFCMIN,PMFU5(JL,JK))
          IF (JK >=KCTOP(JL) )  THEN
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
        ! ZCD5(JL,JK,JN)=.5_JPRB*(ZCU5(JL,JK,JN)+PCEN5(JL,IK,JN))
          ZCD5(JL,JK,JN)=.1_JPRB*ZCU5(JL,JK,JN)+0.9_JPRB*PCEN5(JL,IK,JN)
        ELSEIF ( LDDRAF(JL).AND.JK>KDTOP(JL) ) THEN
          ZERATE25(JL,JK,JN)=-PMFD5(JL,JK)+PMFD5(JL,IK)+PDDRATE5(JL,JK)
          ZERATE45(JL,JK,JN)=ZERATE25(JL,JK,JN)
        ! IF(ZERATE25(JL,JK,JN)<0._JPRB) ZERATE25(JL,JK,JN)=0._JPRB !!Modif
          ZMFA25(JL,JK,JN)=1._JPRB/MIN(-RMFCMIN,PMFD5(JL,JK))
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
     IF( PCEN5(JL,JK,JN)+ZPOSI5*PTSPHY<0.0_JPRB ) THEN
        ZCD5(JL,JK,JN)=( (PMFU5(JL,JK)+PMFD5(JL,JK))*PCEN5(JL,IK,JN)-PMFU5(JL,JK)*ZCU5(JL,JK,JN)&
                    &+PCEN5(JL,JK,JN)/(PTSPHY*ZDP5(JL,JK)) )*ZMFA25(JL,JK,JN)
     ENDIF
    ENDIF
  ENDDO


ENDDO

ZMFC5(:,:,:,:)=0._JPRB

!----------------------------------------------------------------------
 

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


 IF ( RMFSOLCT==0.0_JPRB ) THEN


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
         ! ZTENC5(JL,JK,JN) = (PTENC5(JL,JK)+ZTENC5(JL,JK,JN))*PTSPHY+PCEN5(JL,JK,JN)
           IF(JK<KLEV) THEN
             ZB5(JL,JK,JN,JIT)=1.0_JPRB+ZZP35(JL,JK)*(PMFU5(JL,IK)+PMFD5(JL,IK))
           ELSE
             ZB5(JL,JK,JN,JIT)=1.0_JPRB
           ENDIF
         ENDIF
       ENDDO
    ENDDO
 
! Store trajectory tendency array for current iteration

    ZTENC15(:,:,JN,JIT)=ZTENC5(:,:,JN)

    CALL CUBIDIAG &
       &( KIDIA, KFDIA, KLON, KLEV, &
       &  KCTOP, LLCUMBAS, &
       &  ZMFC45(:,:,JN,JIT),  ZB5(:,:,JN,JIT),   ZTENC5(:,:,JN),   ZR15 )
 
     ! Compute tendencies
 
    DO JK=2,KLEV
       DO JL=KIDIA,KFDIA
         IF(LLCUMBAS(JL,JK)) THEN
           PTENC5(JL,JK,JN)=PTENC5(JL,JK,JN)+(ZR15(JL,JK)-PCEN5(JL,JK,JN))*ZTSPHY
         ! PTENC5(JL,JK,JN)=(ZR15(JL,JK)-PCEN5(JL,JK,JN))*ZTSPHY
         ENDIF
       ENDDO
    ENDDO

  ENDDO
 
 ENDIF
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

ZCEN(:,:,:)=0._JPRB
ZCU(:,:,:)=0._JPRB
ZCD(:,:,:)=0._JPRB
ZTENC(:,:,:)=0._JPRB
ZMFC(:,:,:)=0._JPRB
ZDP(:,:)=0._JPRB

!LOP DO JN=1,KTRAC
!LOP    DO JK=2,KLEV
!LOP       DO JL=KIDIA,KFDIA
!LOP          PTENC(JL,JK,JN)=PTENC(JL,JK,JN)+PCEN(JL,JK,JN)*PTSPHY
!LOP       ENDDO
!LOP    ENDDO
!LOP ENDDO

  IF ( RMFSOLCT==0.0_JPRB ) THEN
    
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

     ZR1(:,:)=0._JPRB
     ZB(:,:)=0._JPRB

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
       &( KIDIA, KFDIA, KLON, KLEV, &
       &  KCTOP, LLCUMBAS, &
       &  ZMFC45(:,:,JN,JIT), ZB5(:,:,JN,JIT), ZTENC15(:,:,JN,JIT), ZR15, &
       &  ZMFC (:,:,JN), ZB , ZTENC (:,:,JN),  ZR1  )
       
     DO JK=KLEV,2,-1
       IK=JK+1
       DO JL=KIDIA,KFDIA
         IF(LLCUMBAS(JL,JK)) THEN
           ZZP=0._JPRB
           IF(JK<KLEV) THEN
!LOP Moved    ZZP=0._JPRB

             PMFU(JL,IK)=PMFU(JL,IK)+ZZP35(JL,JK)*ZB(JL,JK)
             PMFD(JL,IK)=PMFD(JL,IK)+ZZP35(JL,JK)*ZB(JL,JK)
             ZZP=ZZP+(PMFU5(JL,IK)+PMFD5(JL,IK))*ZB(JL,JK)
             ZB(JL,JK)=0._JPRB
           ELSE
             ZB(JL,JK)=0._JPRB
           ENDIF

           PCEN (JL,JK,JN)=PCEN (JL,JK,JN)+ZTENC(JL,JK,JN)
           ZTENC(JL,JK,JN)=ZTENC(JL,JK,JN)*PTSPHY

           ZZP = ZZP - (PMFU5(JL,JK)+PMFD5(JL,JK))*ZMFC(JL,JK,JN)
           PMFU(JL,JK) = PMFU(JL,JK) - ZZP35(JL,JK)*ZMFC(JL,JK,JN) 
           PMFD(JL,JK) = PMFD(JL,JK) - ZZP35(JL,JK)*ZMFC(JL,JK,JN) 
           ZMFC(JL,JK,JN)=0._JPRB

           ZDP(JL,JK) = ZDP(JL,JK) + RMFSOLCT*PTSPHY*ZZP
           ZZP=0._JPRB 
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
         ZTENC(JL,JK,JN)=0._JPRB
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
           ZTENC(JL,JK,JN)=0._JPRB 
         ENDIF
       ENDDO
    ENDDO
    
    
!*    3.0          COMPUTE FLUXES
!                  --------------

    DO JK=KLEV,2,-1
      IK=JK-1
      DO JL=KIDIA,KFDIA
        IF(LLCUMASK(JL,JK)) THEN

          ZMFA=0._JPRB

          PMFU(JL,JK)   = PMFU(JL,JK)   + ZCU5(JL,JK,JN)* ZMFC(JL,JK,JN)
          ZCU(JL,JK,JN) = ZCU(JL,JK,JN) + PMFU5(JL,JK)  * ZMFC(JL,JK,JN)
          PMFD(JL,JK)   = PMFD(JL,JK)   + ZCD5(JL,JK,JN)* ZMFC(JL,JK,JN)
          ZCD(JL,JK,JN) = ZCD(JL,JK,JN) + PMFD5(JL,JK)  * ZMFC(JL,JK,JN)
          ZCEN(JL,IK,JN)= ZCEN(JL,IK,JN)- ZMFA35(JL,JK,JN,JIT)*ZIMP*ZMFC(JL,JK,JN)
          ZMFA          = ZMFA          - ZCEN5(JL,IK,JN)*ZIMP*ZMFC(JL,JK,JN) 
          ZMFC(JL,JK,JN)=0._JPRB

          PMFU(JL,JK)=PMFU(JL,JK)+ZMFA
          PMFD(JL,JK)=PMFD(JL,JK)+ZMFA
          ZMFA=0._JPRB
        ENDIF
      ENDDO
    ENDDO

  ENDDO


ZMFC(:,:,:)=0._JPRB

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
     IF( PCEN5(JL,JK,JN)+ZPOSI5*PTSPHY<0.0_JPRB ) THEN
        ZMFA5=1._JPRB/MIN(-RMFCMIN,PMFD5(JL,JK))
        IF (PMFD5(JL,JK) < -RMFCMIN) THEN
           ZMFA=-PMFD(JL,JK)/(PMFD5(JL,JK)*PMFD5(JL,JK))
        ELSE
           ZMFA=0._JPRB
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
        ZCD(JL,JK,JN)=0._JPRB

     ENDIF
    ENDIF
  ENDDO

    DO JK=KLEV,3,-1
      IK=JK-1
      DO JL=KIDIA,KFDIA

        ZMFA=0._JPRB
        ZERATE=0._JPRB

        IF ( LDDRAF(JL).AND.JK==KDTOP(JL) ) THEN
        ! ZCU(JL,JK,JN)  = ZCU(JL,JK,JN) +.5_JPRB*ZCD(JL,JK,JN)
        ! ZCEN(JL,IK,JN) = ZCEN(JL,IK,JN)+.5_JPRB*ZCD(JL,JK,JN)
          ZCU(JL,JK,JN)  = ZCU(JL,JK,JN) +.1_JPRB*ZCD(JL,JK,JN)
          ZCEN(JL,IK,JN) = ZCEN(JL,IK,JN)+.9_JPRB*ZCD(JL,JK,JN)
          ZCD(JL,JK,JN)=0._JPRB

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
          ZCD(JL,JK,JN)=0._JPRB

          IF (PMFD5(JL,JK) < -RMFCMIN) THEN
            PMFD(JL,JK)=PMFD(JL,JK)-ZMFA/(PMFD5(JL,JK)*PMFD5(JL,JK))
          ENDIF
          ZMFA=0._JPRB

        ! IF (ZERATE45(JL,JK,JN)<0._JPRB) ZERATE=0._JPRB

          PMFD(JL,JK)=PMFD(JL,JK)-ZERATE
          PMFD(JL,IK)=PMFD(JL,IK)+ZERATE
          PDDRATE(JL,JK)=PDDRATE(JL,JK)+ZERATE
          ZERATE=0._JPRB
        ENDIF
      ENDDO
    ENDDO

  
!*    2.0          COMPUTE UPDRAFT VALUES
!                  ----------------------

    DO JK=3,KLEV-1
       IK=JK+1
       DO JL=KIDIA,KFDIA

        ZMFA=0._JPRB
        ZERATE=0._JPRB

        IF ( LLCUMASK(JL,JK) ) THEN

          IF ( JK>=KCTOP(JL) )  THEN
               ZFACT=ZMFA15(JL,JK,JN)*ZCU(JL,JK,JN)
               ZMFA = ZMFA + ( PMFU5(JL,IK)*ZCU5(JL,IK,JN) &
                  &          + ZERATE15(JL,JK,JN)*PCEN5(JL,JK,JN) &
                  &          - PUDRATE5(JL,JK)*ZCU5(JL,IK,JN) )*ZCU(JL,JK,JN)
               PMFU(JL,IK)   = PMFU(JL,IK)   + ZFACT * ZCU5(JL,IK,JN) 
               ZCU(JL,IK,JN) = ZCU(JL,IK,JN) + ZFACT * (PMFU5(JL,IK)-PUDRATE5(JL,JK))
               ZERATE        = ZERATE        + ZFACT * PCEN5(JL,JK,JN)
               PCEN(JL,JK,JN)= PCEN(JL,JK,JN)+ ZFACT * ZERATE15(JL,JK,JN)
               PUDRATE(JL,JK)= PUDRATE(JL,JK)- ZFACT * ZCU5(JL,IK,JN)
               ZCU(JL,JK,JN)=0._JPRB
          ENDIF

          IF (PMFU5(JL,JK) > RMFCMIN) THEN
            PMFU(JL,JK)=PMFU(JL,JK)-ZMFA/(PMFU5(JL,JK)*PMFU5(JL,JK))
          ENDIF
          ZMFA=0._JPRB
          
        ! IF (ZERATE35(JL,JK,JN)<0._JPRB) ZERATE=0._JPRB

          PMFU(JL,JK)=PMFU(JL,JK)+ZERATE
          PMFU(JL,IK)=PMFU(JL,IK)-ZERATE
          PUDRATE(JL,JK)=PUDRATE(JL,JK)+ZERATE
          ZERATE=0._JPRB

        ENDIF
       ENDDO
    ENDDO
  
!*    1.0          DEFINE TRACERS AT HALF LEVELS
!                  -----------------------------

   DO JL=KIDIA,KFDIA
         PCEN(JL,KLEV,JN)=PCEN(JL,KLEV,JN)+ZCU(JL,KLEV,JN)
         ZCU(JL,KLEV,JN)=0._JPRB
   ENDDO
   DO JK=KLEV,2,-1
      IK=JK-1
      DO JL=KIDIA,KFDIA
         PCEN (JL,IK,JN)=PCEN (JL,IK,JN)+ZCD (JL,JK,JN)
         ZCD (JL,JK,JN)=0._JPRB
       
         PCEN(JL,IK,JN)=PCEN(JL,IK,JN)+ZCU(JL,JK,JN)
         ZCU(JL,JK,JN)=0._JPRB

         PCEN(JL,JK,JN)=PCEN(JL,JK,JN)+ZCEN(JL,JK,JN)
         ZCEN(JL,JK,JN)=0._JPRB
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
         ZDP(JL,JK)=0._JPRB
      ENDIF
   ENDDO
ENDDO

DEALLOCATE(ZDP)
DEALLOCATE(ZMFC)
DEALLOCATE(ZTENC)
DEALLOCATE(ZCD)
DEALLOCATE(ZCU)
DEALLOCATE(ZCEN)

IF (RMFSOLCT > 0.0_JPRB) THEN
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
