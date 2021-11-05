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

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK_IFS   ,ONLY : LHOOK,   DR_HOOK

IMPLICIT NONE

!     DUMMY INTEGER 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCTOP(KLON) 
LOGICAL           ,INTENT(IN)    :: LD_LCUMASK(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PA5(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB5(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PR5(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PU5(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PA(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PR(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PU(KLON,KLEV) 
!     DUMMY REALS
!     DUMMY LOGICALS
!     LOCALS
INTEGER(KIND=JPIM) :: JK, JL
REAL(KIND=JPRB) :: ZBET5(KLON,KLEV)
REAL(KIND=JPRB) :: ZBET
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CUBIDIAGAD',0,ZHOOK_HANDLE)

PU5(:,:)=0._JPRB

! Forward Substitution

DO JK = 2, KLEV
  DO JL = KIDIA,KFDIA
    IF ( LD_LCUMASK(JL,JK) ) THEN
      IF ( JK==KCTOP(JL)-1 ) THEN
        ZBET5(JL,JK) =1.0_JPRB/(PB5(JL,JK)+1.E-35_JPRB)
        PU5(JL,JK) = PR5(JL,JK) * ZBET5(JL,JK)
      ELSEIF ( JK>KCTOP(JL)-1 ) THEN
        ZBET5(JL,JK) = 1.0_JPRB/(PB5(JL,JK) + 1.E-35_JPRB)
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
        ZBET=0.0_JPRB

        PR(JL,JK) = PR(JL,JK) + ZBET5(JL,JK)*PU(JL,JK)
        ZBET      = ZBET      + PR5(JL,JK)  *PU(JL,JK)
        PU(JL,JK) =0.0_JPRB

        PB(JL,JK) = PB(JL,JK) - ZBET5(JL,JK)*ZBET5(JL,JK)*ZBET
        ZBET =0.0_JPRB
   
      ELSEIF ( JK>KCTOP(JL)-1 ) THEN
        ZBET=0.0_JPRB

        PR(JL,JK)  = PR(JL,JK) + ZBET5(JL,JK)*PU(JL,JK)
        PA(JL,JK)  = PA(JL,JK) - PU5(JL,JK-1)*ZBET5(JL,JK)*PU(JL,JK)
        PU(JL,JK-1)= PU(JL,JK-1) - PA5(JL,JK)*ZBET5(JL,JK)*PU(JL,JK)
        ZBET = ZBET + (PR5(JL,JK)-PA5(JL,JK)*PU5(JL,JK-1))*PU(JL,JK)
        PU(JL,JK) =0.0_JPRB
   
        PB(JL,JK) = PB(JL,JK) - ZBET5(JL,JK)*ZBET5(JL,JK)*ZBET
        ZBET =0.0_JPRB
      ENDIF
    ENDIF
  ENDDO
ENDDO

PU(:,:)=0._JPRB

IF (LHOOK) CALL DR_HOOK('CUBIDIAGAD',1,ZHOOK_HANDLE)
END SUBROUTINE CUBIDIAGAD
