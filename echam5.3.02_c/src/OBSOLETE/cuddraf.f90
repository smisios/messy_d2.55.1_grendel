SUBROUTINE cuddraf(  kproma, kbdim, klev, klevp1,                      &
           ptenh,    pqenh,    puen,     pven,                         &
           ktrac,                                                      &
           pxtenh,   pxtd,     pmfdxt,                                 &
           pxten,    &                             ! op_ck_20031001
           pgeoh,    paphp1,   prfl,                                   &
           ptd,      pqd,      pud,      pvd,                          &
           pmfd,     pmfds,    pmfdq,    pdmfdp,                       &
           pcpcu,                                                      &
           lddraf)
!
!          THIS ROUTINE CALCULATES CUMULUS DOWNDRAFT DESCENT
!
!          M.TIEDTKE         E.C.M.W.F.    12/86 MODIF. 12/89
!
!          PURPOSE.
!          --------
!          TO PRODUCE THE VERTICAL PROFILES FOR CUMULUS DOWNDRAFTS
!          (I.E. T,Q,U AND V AND FLUXES)
!
!          INTERFACE
!          ---------
!
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!          INPUT IS T,Q,P,PHI,U,V AT HALF LEVELS.
!          IT RETURNS FLUXES OF S,Q AND EVAPORATION RATE
!          AND U,V AT LEVELS WHERE DOWNDRAFT OCCURS
!
!          METHOD.
!          --------
!          CALCULATE MOIST DESCENT FOR ENTRAINING/DETRAINING PLUME BY
!          A) MOVING AIR DRY-ADIABATICALLY TO NEXT LEVEL BELOW AND
!          B) CORRECTING FOR EVAPORATION TO OBTAIN SATURATED STATE.
!
!          EXTERNALS
!          ---------
!          *CUADJTQ* FOR ADJUSTING T AND Q DUE TO EVAPORATION IN
!          SATURATED DESCENT
!
!          REFERENCE
!          ---------
!          (TIEDTKE,1989)
!
USE mo_kind,         ONLY : dp
USE mo_constants,    ONLY : g, rd, vtmpc1
USE mo_cumulus_flux, ONLY : lmfdudv, cmfcmin, entrdd
#ifdef __ibm__
USE mo_specfun,      ONLY : merge
#endif
!
IMPLICIT NONE
!
INTEGER, INTENT (IN) :: kbdim, klev, ktrac, kproma, klevp1
!
REAL(dp) :: ptenh(kbdim,klev),       pqenh(kbdim,klev),                &
            puen(kbdim,klev),        pven(kbdim,klev),                 &
            pgeoh(kbdim,klev),       paphp1(kbdim,klevp1)
!
REAL(dp) :: ptd(kbdim,klev),         pqd(kbdim,klev),                  &
            pud(kbdim,klev),         pvd(kbdim,klev),                  &
            pmfd(kbdim,klev),        pmfds(kbdim,klev),                &
            pmfdq(kbdim,klev),       pdmfdp(kbdim,klev),               &
            prfl(kbdim)
REAL(dp) :: pcpcu(kbdim,klev)
LOGICAL  :: lddraf(kbdim)
!
REAL(dp) :: zdmfen(kbdim),           zdmfde(kbdim),                    &
            zcond(kbdim)
REAL(dp) :: zph(kbdim)
REAL(dp) :: pxten(kbdim,klev,ktrac) ! op_ck_20031001
LOGICAL  :: llo2(kbdim)
REAL(dp) :: pxtenh(kbdim,klev,ktrac),pxtd(kbdim,klev,ktrac),           &
            pmfdxt(kbdim,klev,ktrac)
LOGICAL  :: llo1
INTEGER  :: jk, is, jl, itopde, jt, ik, icall
REAL(dp) :: zentr, zseen, zqeen, zsdde, zqdde, zmfdsk, zmfdqk, zxteen  &
          , zxtdde, zmfdxtk, zbuo, zdmfdp, zmfduk, zmfdvk
!
!
!----------------------------------------------------------------------
!
!     1.           CALCULATE MOIST DESCENT FOR CUMULUS DOWNDRAFT BY
!                     (A) CALCULATING ENTRAINMENT RATES, ASSUMING
!                         LINEAR DECREASE OF MASSFLUX IN PBL
!                     (B) DOING MOIST DESCENT - EVAPORATIVE COOLING
!                         AND MOISTENING IS CALCULATED IN *CUADJTQ*
!                     (C) CHECKING FOR NEGATIVE BUOYANCY AND
!                         SPECIFYING FINAL T,Q,U,V AND DOWNWARD FLUXES
!                    -------------------------------------------------
!
100 CONTINUE
  DO 180 jk=3,klev
     is=0
     ! mz_ht_20040106+
     zdmfen(:) = 0._dp
     zdmfde(:) = 0._dp
     ! mz_ht_20040106-
     DO 110 jl=1,kproma
        zph(jl)=paphp1(jl,jk)
        llo2(jl)=lddraf(jl).AND.pmfd(jl,jk-1).LT.0._dp
        is=is+MERGE(1,0,llo2(jl))
110  END DO
     IF(is.EQ.0) go to 180
     DO 122 jl=1,kproma
        IF(llo2(jl)) THEN
           zentr=entrdd*pmfd(jl,jk-1)*rd*ptenh(jl,jk-1)/               &
                          (g*paphp1(jl,jk-1))*                         &
                          (paphp1(jl,jk)-paphp1(jl,jk-1))
           zdmfen(jl)=zentr
           zdmfde(jl)=zentr
        END IF
122  END DO
     itopde=klev-2
     IF(jk.GT.itopde) THEN
        DO 124 jl=1,kproma
           IF(llo2(jl)) THEN
              zdmfen(jl)=0._dp
              zdmfde(jl)=pmfd(jl,itopde)*                              &
                             (paphp1(jl,jk)-paphp1(jl,jk-1))/          &
                             (paphp1(jl,klevp1)-paphp1(jl,itopde))
           END IF
124     END DO
     END IF
!
     DO 126 jl=1,kproma
        IF(llo2(jl)) THEN
           pmfd(jl,jk)=pmfd(jl,jk-1)+zdmfen(jl)-zdmfde(jl)
           zseen=(pcpcu(jl,jk-1)*ptenh(jl,jk-1)+pgeoh(jl,jk-1))        &
                                                      *zdmfen(jl)
           zqeen=pqenh(jl,jk-1)*zdmfen(jl)
           zsdde=(pcpcu(jl,jk-1)*ptd(jl,jk-1)+pgeoh(jl,jk-1))*zdmfde(jl)
           zqdde=pqd(jl,jk-1)*zdmfde(jl)
           zmfdsk=pmfds(jl,jk-1)+zseen-zsdde
           zmfdqk=pmfdq(jl,jk-1)+zqeen-zqdde
           pqd(jl,jk)=zmfdqk*(1._dp/MIN(-cmfcmin,pmfd(jl,jk)))
           ptd(jl,jk)=(zmfdsk*(1._dp/MIN(-cmfcmin,pmfd(jl,jk)))-       &
                                   pgeoh(jl,jk))/pcpcu(jl,jk)
           ptd(jl,jk)=MIN(400._dp,ptd(jl,jk))
           ptd(jl,jk)=MAX(100._dp,ptd(jl,jk))
           zcond(jl)=pqd(jl,jk)
        END IF
126  END DO
!
     DO 1264 jt=1,ktrac
        DO 1262 jl=1,kproma
           IF(llo2(jl)) THEN
! op_ck_20031001+
!             zxteen=pxtenh(jl,jk-1,jt)*zdmfen(jl)
!             zxtdde=pxtd(jl,jk-1,jt)*zdmfde(jl)
!             zmfdxtk=pmfdxt(jl,jk-1,jt)+zxteen-zxtdde
!             pxtd(jl,jk,jt)=zmfdxtk*(1./MIN(-cmfcmin,pmfd(jl,jk)))
              zxteen=pxten(jl,jk-1,jt)*zdmfen(jl)
              zmfdxtk=pmfdxt(jl,jk-1,jt)+zxteen
              pmfdxt(jl,jk,jt)=zmfdxtk/(1.+zdmfde(JL) /                 &
                              (min(pmfd(jl,jk),-cmfcmin)))
              pxtd(jl,jk-1,jt)=pmfdxt(jl,jk,jt)       /                 &
                              (min(pmfd(JL,JK),-cmfcmin))
! op_ck_20031001-
           ENDIF
1262    END DO
1264 END DO
!
!
     ik=jk
     icall=2
     CALL cuadjtq(kproma,   kbdim,    klev,     ik,                    &
                  zph,      ptd,      pqd,      llo2,     icall)
!
     DO 150 jl=1,kproma
        IF(llo2(jl)) THEN
           zcond(jl)=zcond(jl)-pqd(jl,jk)
           zbuo=ptd(jl,jk)*(1._dp+vtmpc1*pqd(jl,jk))-                  &
                       ptenh(jl,jk)*(1._dp+vtmpc1*pqenh(jl,jk))
           llo1=zbuo.LT.0._dp.AND.                                     &
                             (prfl(jl)-pmfd(jl,jk)*zcond(jl).GT.0._dp)
           pmfd(jl,jk)=MERGE(pmfd(jl,jk),0._dp,llo1)
           pmfds(jl,jk)=(pcpcu(jl,jk)*ptd(jl,jk)+pgeoh(jl,jk))         &
                                                    *pmfd(jl,jk)
           pmfdq(jl,jk)=pqd(jl,jk)*pmfd(jl,jk)
           zdmfdp=-pmfd(jl,jk)*zcond(jl)
           pdmfdp(jl,jk-1)=zdmfdp
           prfl(jl)=prfl(jl)+zdmfdp
        END IF
150  END DO
!
     DO 1504 jt=1,ktrac
        DO 1502 jl=1,kproma
           IF(llo2(jl)) THEN
!             pmfdxt(jl,jk,jt)=pxtd(jl,jk,jt)*pmfd(jl,jk)   ! op_ck_20031001
              pmfdxt(jl,jk,jt)=pmfd(jl,jk)*pxtd(jl,jk-1,jt) ! op_ck_20031001
           ENDIF
1502    END DO
1504 END DO
!
     IF(lmfdudv) THEN
        DO 160 jl=1,kproma
           IF(llo2(jl).AND.pmfd(jl,jk).LT.0._dp) THEN
              zmfduk=pmfd(jl,jk-1)*pud(jl,jk-1)+                       &
                              zdmfen(jl)*puen(jl,jk-1)-                &
                              zdmfde(jl)*pud(jl,jk-1)
              zmfdvk=pmfd(jl,jk-1)*pvd(jl,jk-1)+                       &
                              zdmfen(jl)*pven(jl,jk-1)-                &
                              zdmfde(jl)*pvd(jl,jk-1)
              pud(jl,jk)=zmfduk*(1._dp/MIN(-cmfcmin,pmfd(jl,jk)))
              pvd(jl,jk)=zmfdvk*(1._dp/MIN(-cmfcmin,pmfd(jl,jk)))
           END IF
160     END DO
     END IF
!
180 END DO
!
  RETURN
END SUBROUTINE cuddraf
