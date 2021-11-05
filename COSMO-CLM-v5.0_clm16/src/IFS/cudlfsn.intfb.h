INTERFACE
SUBROUTINE CUDLFSN&
 & (KIDIA, KFDIA, KLON, KTDIA, KLEV,&
 & KCBOT, KCTOP, LDLAND, LDCUM,&
 & PTENH, PQENH, PUEN, PVEN,&
 & PTEN, PQSEN, PGEO,&
 & PGEOH, PAPH, PTU, PQU, PLU,&
 & PUU, PVU, PMFUB, PRFL,&
 & PTD, PQD,&
 & PMFD, PMFDS, PMFDQ, PDMFDP,&
 & KDTOP, LDDRAF) 
USE PARKIND1 ,ONLY : JPIM ,JPRB
INTEGER(KIND=JPIM),INTENT(IN) :: KLON
INTEGER(KIND=JPIM),INTENT(IN) :: KLEV
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KFDIA
INTEGER(KIND=JPIM) :: KTDIA
INTEGER(KIND=JPIM) :: KCBOT(KLON)
INTEGER(KIND=JPIM) :: KCTOP(KLON)
LOGICAL :: LDLAND(KLON)
LOGICAL :: LDCUM(KLON)
REAL(KIND=JPRB) ,INTENT(IN) :: PTENH(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PQENH(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PUEN(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PVEN(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PTEN(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PQSEN(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PGEO(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PGEOH(KLON,KLEV+1)
REAL(KIND=JPRB) ,INTENT(IN) :: PAPH(KLON,KLEV+1)
REAL(KIND=JPRB) ,INTENT(IN) :: PTU(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PQU(KLON,KLEV)
REAL(KIND=JPRB) :: PLU(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PUU(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PVU(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PMFUB(KLON)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PRFL(KLON)
REAL(KIND=JPRB) ,INTENT(OUT) :: PTD(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(OUT) :: PQD(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PMFD(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(OUT) :: PMFDS(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(OUT) :: PMFDQ(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(OUT) :: PDMFDP(KLON,KLEV)
INTEGER(KIND=JPIM),INTENT(OUT) :: KDTOP(KLON)
LOGICAL ,INTENT(OUT) :: LDDRAF(KLON)
END SUBROUTINE CUDLFSN
END INTERFACE
