INTERFACE
SUBROUTINE CUININ&
 & ( KIDIA, KFDIA, KLON, KTDIA, KLEV,&
 & PTEN, PQEN, PQSEN, PUEN, PVEN,&
 & PVERVEL, PGEO, PAPH, PAP,&
 & KLWMIN, KLAB,&
 & PTENH, PQENH, PQSENH, PGEOH,&
 & PTU, PQU, PTD, PQD,&
 & PUU, PVU, PUD, PVD,&
 & PLU ) 
USE PARKIND1 ,ONLY : JPIM ,JPRB
INTEGER(KIND=JPIM),INTENT(IN) :: KLON
INTEGER(KIND=JPIM),INTENT(IN) :: KLEV
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KFDIA
INTEGER(KIND=JPIM) :: KTDIA
REAL(KIND=JPRB) ,INTENT(IN) :: PTEN(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PQEN(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PQSEN(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PUEN(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PVEN(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PVERVEL(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PGEO(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PAPH(KLON,KLEV+1)
REAL(KIND=JPRB) :: PAP(KLON,KLEV)
INTEGER(KIND=JPIM),INTENT(OUT) :: KLWMIN(KLON)
INTEGER(KIND=JPIM),INTENT(OUT) :: KLAB(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PTENH(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(OUT) :: PQENH(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PQSENH(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PGEOH(KLON,KLEV+1)
REAL(KIND=JPRB) ,INTENT(OUT) :: PTU(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(OUT) :: PQU(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(OUT) :: PTD(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(OUT) :: PQD(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(OUT) :: PUU(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(OUT) :: PVU(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(OUT) :: PUD(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(OUT) :: PVD(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(OUT) :: PLU(KLON,KLEV)
END SUBROUTINE CUININ
END INTERFACE
