INTERFACE
SUBROUTINE CUBIDIAGAD&
 & ( KIDIA, KFDIA, KLON, KLEV,&
 & KCTOP, LD_LCUMASK,&
 & PA5, PB5, PR5, PU5,&
 & PA , PB , PR , PU ) 
USE PARKIND1 ,ONLY : JPIM ,JPRB
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KFDIA
INTEGER(KIND=JPIM) :: KLON
INTEGER(KIND=JPIM),INTENT(IN) :: KLEV
INTEGER(KIND=JPIM),INTENT(IN) :: KCTOP(KLON)
LOGICAL ,INTENT(IN) :: LD_LCUMASK(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PA5(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PB5(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PR5(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(OUT)   :: PU5(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PA(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PB(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PR(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(OUT)   :: PU(KLON,KLEV)
END SUBROUTINE CUBIDIAGAD
END INTERFACE
