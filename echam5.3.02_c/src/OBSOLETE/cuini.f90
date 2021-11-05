SUBROUTINE cuini(kproma, kbdim, klev, klevp1, klevm1,                  &
           pten,     pqen,     pqsen,    pxen,     puen,     pven,     &
           ptven,    ktrac,                                            &
           zpmfun,   &                        ! op_ck_20031001
           pxten,    pxtenh,   pxtu,     pxtd,     pmfuxt,   pmfdxt,   &
           pverv,    pgeo,     paphp1,   pgeoh,                        &
           ptenh,    pqenh,    pqsenh,   pxenh,    klwmin,             &
           ptu,      pqu,      ptd,      pqd,                          &
           puu,      pvu,      pud,      pvd,                          &
           pmfu,     pmfd,     pmfus,    pmfds,                        &
           pmfuq,    pmfdq,    pdmfup,   pdmfdp,                       &
           pcpen,    pcpcu,                                            &
           pdpmel,   plu,      plude,    pqude,    klab)
!
!          M.TIEDTKE         E.C.M.W.F.     12/89
!
!          PURPOSE
!          -------
!
!          THIS ROUTINE INTERPOLATES LARGE-SCALE FIELDS OF T,Q ETC.
!          TO HALF LEVELS (I.E. GRID FOR MASSFLUX SCHEME),
!          DETERMINES LEVEL OF MAXIMUM VERTICAL VELOCITY
!          AND INITIALIZES VALUES FOR UPDRAFTS AND DOWNDRAFTS
!
!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!
!          METHOD.
!          --------
!          FOR EXTRAPOLATION TO HALF LEVELS SEE TIEDTKE(1989)
!
!          EXTERNALS
!          ---------
!          *CUADJTQ* TO SPECIFY QS AT HALF LEVELS
!
USE mo_kind,      ONLY: dp
USE mo_constants, ONLY: rd, cpd

IMPLICIT NONE

INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, klevm1, ktrac

REAL(dp):: pten(kbdim,klev),          pqen(kbdim,klev),                &
           puen(kbdim,klev),          pven(kbdim,klev),                &
           pqsen(kbdim,klev),         pverv(kbdim,klev),               &
           pgeo(kbdim,klev),          pgeoh(kbdim,klev),               &
           paphp1(kbdim,klevp1),      ptenh(kbdim,klev),               &
           pxenh(kbdim,klev),         pxen(kbdim,klev),                &
           ptven(kbdim,klev),                                          &
           pqenh(kbdim,klev),         pqsenh(kbdim,klev)
REAL(dp):: pcpen(kbdim,klev),         pcpcu(kbdim,klev)
!
REAL(dp):: ptu(kbdim,klev),           pqu(kbdim,klev),                 &
           ptd(kbdim,klev),           pqd(kbdim,klev),                 &
           puu(kbdim,klev),           pud(kbdim,klev),                 &
           pvu(kbdim,klev),           pvd(kbdim,klev),                 &
           pmfu(kbdim,klev),          pmfd(kbdim,klev),                &
           pmfus(kbdim,klev),         pmfds(kbdim,klev),               &
           pmfuq(kbdim,klev),         pmfdq(kbdim,klev),               &
           pdmfup(kbdim,klev),        pdmfdp(kbdim,klev),              &
           plu(kbdim,klev),           plude(kbdim,klev),               &
           pqude(kbdim,klev)
REAL(dp):: pdpmel(kbdim,klev)
INTEGER :: klab(kbdim,klev),          klwmin(kbdim)
!
REAL(dp):: zwmax(kbdim)
REAL(dp):: zpmfun(kbdim,klev) ! op_ck_20031001
REAL(dp):: zph(kbdim)
LOGICAL :: loflag(kbdim)
REAL(dp):: pxten(kbdim,klev,ktrac),   pxtenh(kbdim,klev,ktrac),        &
           pxtu(kbdim,klev,ktrac),    pxtd(kbdim,klev,ktrac),          &
           pmfuxt(kbdim,klev,ktrac),  pmfdxt(kbdim,klev,ktrac)

INTEGER :: jk, jl, jt, ik, icall
REAL(dp):: zarg, zcpm, zzs
!
!  INTRINSIC FUNCTIONS
INTRINSIC MAX, MIN
!
!----------------------------------------------------------------------
!
!*    1.           SPECIFY LARGE SCALE PARAMETERS AT HALF LEVELS
!*                 ADJUST TEMPERATURE FIELDS IF STATICLY UNSTABLE
!*                 FIND LEVEL OF MAXIMUM VERTICAL VELOCITY
!                  ----------------------------------------------
!
100 CONTINUE
    DO 101 jk=1,klev
       DO 102 jl=1,kproma
!          pcpen(jl,jk)=cpd*(1.+vtmpc2*MAX(pqen(jl,jk),0.0_dp))
          pcpen(jl,jk)=cpd
102    END DO
101 END DO
  DO 105 jl=1,kproma
     zarg=paphp1(jl,klevp1)/paphp1(jl,klev)
     pgeoh(jl,klev)=rd*ptven(jl,klev)*LOG(zarg)
105 END DO
  DO 107 jk=klevm1,2,-1
     DO 106 jl=1,kproma
        zarg=paphp1(jl,jk+1)/paphp1(jl,jk)
        pgeoh(jl,jk)=pgeoh(jl,jk+1)+rd*ptven(jl,jk)*LOG(zarg)
106  END DO
107 END DO
  DO 130 jk=2,klev
     DO 110 jl=1,kproma
        zcpm=(pcpen(jl,jk)+pcpen(jl,jk-1))*0.5_dp
        ptenh(jl,jk)=(MAX(pcpen(jl,jk-1)*pten(jl,jk-1)+pgeo(jl,jk-1),  &
                   pcpen(jl,jk)*pten(jl,jk)+pgeo(jl,jk))-pgeoh(jl,jk)) &
                    /zcpm
        pqsenh(jl,jk)=pqsen(jl,jk-1)
        zph(jl)=paphp1(jl,jk)
        loflag(jl)=.TRUE.
110  END DO
!
     DO 1104 jt=1,ktrac
        DO 1102 jl=1,kproma
           pxtenh(jl,jk,jt)=(pxten(jl,jk,jt)+pxten(jl,jk-1,jt))*0.5_dp
1102    END DO
1104 END DO
!
!
     ik=jk
     icall=0
     CALL cuadjtq(kproma, kbdim, klev, ik,                             &
                  zph,      ptenh,    pqsenh,   loflag,   icall)
!
     DO 120 jl=1,kproma
        pxenh(jl,jk)=(pxen(jl,jk)+pxen(jl,jk-1))*0.5_dp
        pqenh(jl,jk)=MIN(pqen(jl,jk-1),pqsen(jl,jk-1))                 &
                          +(pqsenh(jl,jk)-pqsen(jl,jk-1))
        pqenh(jl,jk)=MAX(pqenh(jl,jk),0._dp)
!        pcpcu(jl,jk)=cpd*(1.+vtmpc2*pqenh(jl,jk))
        pcpcu(jl,jk)=cpd
120  END DO
130 END DO
!
  DO 140 jl=1,kproma
     ptenh(jl,klev)=(pcpen(jl,klev)*pten(jl,klev)+pgeo(jl,klev)-       &
                           pgeoh(jl,klev))/pcpen(jl,klev)
     pxenh(jl,klev)=pxen(jl,klev)
     pqenh(jl,klev)=pqen(jl,klev)
     pcpcu(jl,1)=pcpen(jl,1)
     ptenh(jl,1)=pten(jl,1)
     pxenh(jl,1)=pxen(jl,1)
     pqenh(jl,1)=pqen(jl,1)
     pgeoh(jl,1)=pgeo(jl,1)
     klwmin(jl)=klev
     zwmax(jl)=0._dp
140 END DO
!
  DO 1404 jt=1,ktrac
     DO 1402 jl=1,kproma
        pxtenh(jl,klev,jt)=pxten(jl,klev,jt)
        pxtenh(jl,1,jt)=pxten(jl,1,jt)
1402 END DO
1404 END DO
!
!
  DO 160 jk=klevm1,2,-1
     DO 150 jl=1,kproma
        zzs=MAX(pcpcu(jl,jk)*ptenh(jl,jk)+pgeoh(jl,jk),                &
                      pcpcu(jl,jk+1)*ptenh(jl,jk+1)+pgeoh(jl,jk+1))
        ptenh(jl,jk)=(zzs-pgeoh(jl,jk))/pcpcu(jl,jk)
150  END DO
160 END DO
!
  DO 190 jk=klev,3,-1
!DIR$ IVDEP
!OCL NOVREC
     DO 180 jl=1,kproma
        IF(pverv(jl,jk).LT.zwmax(jl)) THEN
           zwmax(jl)=pverv(jl,jk)
           klwmin(jl)=jk
        END IF
180  END DO
190 END DO
!
!
!-----------------------------------------------------------------------
!*    2.0          INITIALIZE VALUES FOR UPDRAFTS AND DOWNDRAFTS
!*                 ---------------------------------------------
!
200 CONTINUE
  DO 230 jk=1,klev
     ik=jk-1
     IF(jk.EQ.1) ik=1
     DO 220 jl=1,kproma
        ptu(jl,jk)=ptenh(jl,jk)
        ptd(jl,jk)=ptenh(jl,jk)
        pqu(jl,jk)=pqenh(jl,jk)
        pqd(jl,jk)=pqenh(jl,jk)
        plu(jl,jk)=0._dp
        puu(jl,jk)=puen(jl,ik)
        pud(jl,jk)=puen(jl,ik)
        pvu(jl,jk)=pven(jl,ik)
        pvd(jl,jk)=pven(jl,ik)
        pmfu(jl,jk)=0._dp
        pmfd(jl,jk)=0._dp
        zpmfun(jl,jk)=0._dp ! op_ck_20031001
        pmfus(jl,jk)=0._dp
        pmfds(jl,jk)=0._dp
        pmfuq(jl,jk)=0._dp
        pmfdq(jl,jk)=0._dp
        pdmfup(jl,jk)=0._dp
        pdmfdp(jl,jk)=0._dp
        pdpmel(jl,jk)=0._dp
        plude(jl,jk)=0._dp
        pqude(jl,jk)=0._dp
        klab(jl,jk)=0
220  END DO
!
     DO 2204 jt=1,ktrac
        DO 2202 jl=1,kproma
           pxtu(jl,jk,jt)=pxtenh(jl,jk,jt)
           pxtd(jl,jk,jt)=pxtenh(jl,jk,jt)
           pmfuxt(jl,jk,jt)=0._dp
           pmfdxt(jl,jk,jt)=0._dp
2202    END DO
2204 END DO
!
230 END DO
!
  RETURN
END SUBROUTINE cuini
