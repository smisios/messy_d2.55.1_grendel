INTERFACE
SUBROUTINE CUSTRAT&
 & ( KIDIA, KFDIA, KLON, KTDIA, KLEV,&
 & LDCUM, PTSPHY,&
 & PAP, PAPH, PGEO,&
 & PTEN, PQEN, PQSAT, PENTH,&
 & PTENT, PTENQ ) 
USE PARKIND1 ,ONLY : JPIM ,JPRB
INTEGER(KIND=JPIM),INTENT(IN) :: KLON
INTEGER(KIND=JPIM),INTENT(IN) :: KLEV
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KFDIA
INTEGER(KIND=JPIM) :: KTDIA
LOGICAL ,INTENT(IN) :: LDCUM(KLON)
REAL(KIND=JPRB) ,INTENT(IN) :: PTSPHY
REAL(KIND=JPRB) ,INTENT(IN) :: PAP(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PAPH(KLON,KLEV+1)
REAL(KIND=JPRB) ,INTENT(IN) :: PGEO(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PTEN(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PQEN(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PQSAT(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(OUT) :: PENTH(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PTENT(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PTENQ(KLON,KLEV)
END SUBROUTINE CUSTRAT
END INTERFACE
