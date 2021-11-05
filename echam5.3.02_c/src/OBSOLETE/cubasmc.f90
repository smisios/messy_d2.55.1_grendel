SUBROUTINE cubasmc(  kproma, kbdim, klev, kk, klab,                    &
           pten,     pqen,     pqsen,    puen,     pven,               &
           ktrac,                                                      &
           paphp1,   zpmfun,   & ! op_ck_20031001
           pxten,    pxtu,     pmfuxt,                                 &
           pverv,    pgeo,     pgeoh,    ldcum,    ktype,              &
           pmfu,     pmfub,    pentr,    kcbot,                        &
           ptu,      pqu,      plu,      puu,      pvu,                &
           pmfus,    pmfuq,    pmful,    pdmfup,   pmfuu,              &
           pcpen,                                                      &
           pmfuv)
!
!          M.TIEDTKE         E.C.M.W.F.     12/89
!
!          PURPOSE.
!          --------
!          THIS ROUTINE CALCULATES CLOUD BASE VALUES
!          FOR MIDLEVEL CONVECTION
!
!          INTERFACE
!          ---------
!
!          THIS ROUTINE IS CALLED FROM *CUASC*.
!          INPUT ARE ENVIRONMENTAL VALUES T,Q ETC
!          IT RETURNS CLOUDBASE VALUES FOR MIDLEVEL CONVECTION
!
!          METHOD.
!          --------
!          S. TIEDTKE (1989)
!
!          EXTERNALS
!          ---------
!          NONE
!
USE mo_kind,          ONLY : dp
USE mo_constants,     ONLY : g
USE mo_cumulus_flux,  ONLY : lmfdudv, entrmid, cmfcmin, cmfcmax
!
IMPLICIT NONE
!
INTEGER, INTENT (IN) :: kbdim, klev, ktrac, kproma, kk
REAL(dp) :: pten(kbdim,klev),        pqen(kbdim,klev),                 &
            puen(kbdim,klev),        pven(kbdim,klev),                 &
            pqsen(kbdim,klev),       pverv(kbdim,klev),                &
            pgeo(kbdim,klev),        pgeoh(kbdim,klev)
!
REAL(dp) :: ptu(kbdim,klev),         pqu(kbdim,klev),                  &
            puu(kbdim,klev),         pvu(kbdim,klev),                  &
            plu(kbdim,klev),         pmfu(kbdim,klev),                 &
            pmfub(kbdim),            pentr(kbdim),                     &
            pmfus(kbdim,klev),       pmfuq(kbdim,klev),                &
            pmful(kbdim,klev),       pdmfup(kbdim,klev),               &
            pmfuu(kbdim),            pmfuv(kbdim)
REAL(dp) :: zpmfun(kbdim,klev)   ! op_ck_20031001
REAL(dp) :: paphp1(kbdim,klev+1) ! op_ck_20031001
REAL(dp) :: zzp                  ! op_ck_20031001
INTEGER  :: ikb, jk              ! op_ck_20031001
REAL(dp) :: pcpen(kbdim,klev)
INTEGER  :: ktype(kbdim),            kcbot(kbdim),                     &
            klab(kbdim,klev)
LOGICAL  :: ldcum(kbdim)
!
REAL(dp) :: pxten(kbdim,klev,ktrac), pxtu(kbdim,klev,ktrac),           &
            pmfuxt(kbdim,klev,ktrac)
LOGICAL  :: llo3(kbdim)
INTEGER  :: jl, jt
REAL(dp) :: zzzmb
!
!
!----------------------------------------------------------------------
!
!*    1.           CALCULATE ENTRAINMENT AND DETRAINMENT RATES
!                  -------------------------------------------
!
100 CONTINUE
!DIR$ IVDEP
!OCL NOVREC
  DO 150 jl=1,kproma
     llo3(jl)=.FALSE.
     IF(.NOT.ldcum(jl).AND.klab(jl,kk+1).EQ.0                          &
              .AND.pqen(jl,kk).GT.0.90_dp*pqsen(jl,kk)) THEN
        llo3(jl)=.TRUE.
        ptu(jl,kk+1)=(pcpen(jl,kk)*pten(jl,kk)                         &
                        +pgeo(jl,kk)-pgeoh(jl,kk+1))/pcpen(jl,kk)
        pqu(jl,kk+1)=pqen(jl,kk)
        plu(jl,kk+1)=0._dp
        zzzmb=MAX(cmfcmin,-pverv(jl,kk)/g)
        zzzmb=MIN(zzzmb,cmfcmax)
        pmfub(jl)=zzzmb
        pmfu(jl,kk+1)=pmfub(jl)
        zpmfun(jl,kk+1)=pmfub(jl) ! op_ck_20031001
        pmfus(jl,kk+1)=pmfub(jl)*(pcpen(jl,kk)*ptu(jl,kk+1)            &
                                        +pgeoh(jl,kk+1))
        pmfuq(jl,kk+1)=pmfub(jl)*pqu(jl,kk+1)
        pmful(jl,kk+1)=0._dp
        pdmfup(jl,kk+1)=0._dp
        kcbot(jl)=kk
        klab(jl,kk+1)=1
        ktype(jl)=3
        pentr(jl)=entrmid
        IF(lmfdudv) THEN
           puu(jl,kk+1)=puen(jl,kk)
           pvu(jl,kk+1)=pven(jl,kk)
           pmfuu(jl)=pmfub(jl)*puu(jl,kk+1)
           pmfuv(jl)=pmfub(jl)*pvu(jl,kk+1)
        END IF
     END IF
150 END DO
! op_ck_20031001+
! !DIR$ IVDEP
! !OCL NOVREC
!  DO 1504 jt=1,ktrac
!     DO 1502 jl=1,kproma
!        IF (llo3(jl)) THEN
!           pxtu(jl,kk+1,jt)=pxten(jl,kk,jt)
!           pmfuxt(jl,kk+1,jt)=pmfub(jl)*pxtu(jl,kk+1,jt)
!        ENDIF
!1502 END DO
!1504 END DO

  DO 1500 jk=1,klev
!DIR$ IVDEP
!OCL NOVREC
     DO 1501 jl=1,kproma
        IF (llo3(jl)) THEN
          ikb=kcbot(jl)+1
          IF (jk.gt.ikb) THEN
            ZZP=(paphp1(jl,klev+1)-paphp1(jl,jk))                   &
               /(paphp1(jl,klev+1)-paphp1(jl,ikb))
            zpmfun(jl,jk)=ZPMFUN(jl,ikb)*zzp*zzp
          ENDIF
        ENDIF
1501    ENDDO
1500 ENDDO
!
  DO 1503 jt=1,ktrac
     DO 1504 jk=kk+1,klev
        DO 1505 jl=1,kproma
           IF (llo3(jl)) THEN
             pxtu(jl,jk,jt)=pxten(jl,jk,jt)
             pmfuxt(jl,jk,jt)=pxtu(jl,jk,jt)*zpmfun(jl,jk)
           ENDIF   
1505       ENDDO
1504    ENDDO
1503 ENDDO
! op_ck_20031001-
!
!
  RETURN
END SUBROUTINE cubasmc
