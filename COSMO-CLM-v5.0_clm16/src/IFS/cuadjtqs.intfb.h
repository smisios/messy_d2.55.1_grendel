INTERFACE
SUBROUTINE CUADJTQS&
 & (KIDIA, KFDIA, KLON, KTDIA, KLEV,&
 & KK,&
 & PSP, PT, PQ, LDFLAG, KCALL) 
USE PARKIND1 ,ONLY : JPIM ,JPRB
INTEGER(KIND=JPIM),INTENT(IN) :: KLON
INTEGER(KIND=JPIM),INTENT(IN) :: KLEV
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KFDIA
INTEGER(KIND=JPIM) :: KTDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KK
REAL(KIND=JPRB) ,INTENT(IN) :: PSP(KLON)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PT(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PQ(KLON,KLEV)
LOGICAL ,INTENT(IN) :: LDFLAG(KLON)
INTEGER(KIND=JPIM),INTENT(IN) :: KCALL
END SUBROUTINE CUADJTQS
END INTERFACE
