INTERFACE
SUBROUTINE CUMASTRN&
 & ( KIDIA, KFDIA, KLON, KTDIA, KLEV,&
 & KSTEP, KSTART, LDLAND, PTSPHY,&
 & PTEN, PQEN, PUEN, PVEN, PLITOT,&
 & PVERVEL, PQHFL, PAHFS,&
 & psstru, psstrv,&
 & PAP, PAPH, PGEO, PGEOH,&
 & PTENT, PTENQ, PTENU, PTENV,&
 &,PTENL, PTENI &
 & LDCUM, KTYPE, KCBOT, KCTOP,&
 & KBOTSC, LDSC,&
 & PTU, PQU, PLU, &
 & PMFLXR, PMFLXS, PRAIN,&
 & PMFU, PMFD,&
 & PMFUDE_RATE, PMFDDE_RATE, PCAPE,&
 & KTRAC, PCEN, PTENC ) 
USE PARKIND1 ,ONLY : JPIM ,JPRB
INTEGER(KIND=JPIM),INTENT(IN) :: KLON
INTEGER(KIND=JPIM),INTENT(IN) :: KLEV
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KTRAC
INTEGER(KIND=JPIM) :: KTDIA
INTEGER(KIND=JPIM) :: KSTEP
INTEGER(KIND=JPIM) :: KSTART
LOGICAL ,INTENT(IN) :: LDLAND(KLON)
REAL(KIND=JPRB) ,INTENT(IN) :: PTSPHY
REAL(KIND=JPRB) ,INTENT(INOUT) :: PTEN(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PQEN(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PUEN(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PVEN(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PLITOT(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PVERVEL(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PQHFL(KLON,KLEV+1)
REAL(KIND=JPRB) ,INTENT(IN) :: PAHFS(KLON,KLEV+1)
REAL(KIND=JPRB) :: PSSTRU(KLON)
REAL(KIND=JPRB) :: PSSTRV(KLON)
REAL(KIND=JPRB) ,INTENT(IN) :: PAP(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PAPH(KLON,KLEV+1)
REAL(KIND=JPRB) ,INTENT(IN) :: PGEO(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PGEOH(KLON,KLEV+1)
REAL(KIND=JPRB) ,INTENT(IN) :: PCEN(KLON,KLEV,KTRAC)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PTENT(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PTENQ(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(OUT)   :: PTENL(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(OUT)   :: PTENI(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PTENU(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PTENV(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PTENC(KLON,KLEV,KTRAC)
LOGICAL ,INTENT(INOUT) :: LDCUM(KLON)
INTEGER(KIND=JPIM),INTENT(INOUT) :: KTYPE(KLON)
INTEGER(KIND=JPIM),INTENT(INOUT) :: KCBOT(KLON)
INTEGER(KIND=JPIM),INTENT(INOUT) :: KCTOP(KLON)
INTEGER(KIND=JPIM),INTENT(OUT) :: KBOTSC(KLON)
LOGICAL ,INTENT(OUT) :: LDSC(KLON)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PTU(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PQU(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PLU(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PMFLXR(KLON,KLEV+1)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PMFLXS(KLON,KLEV+1)
REAL(KIND=JPRB) ,INTENT(OUT) :: PRAIN(KLON)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PMFU(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PMFD(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PMFUDE_RATE(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PMFDDE_RATE(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(OUT) :: PCAPE(KLON)
END SUBROUTINE CUMASTRN
END INTERFACE
