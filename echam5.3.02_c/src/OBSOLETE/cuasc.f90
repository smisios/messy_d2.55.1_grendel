SUBROUTINE cuasc(    kproma, kbdim, klev, klevp1, klevm1,               &
           ptenh,    pqenh,    puen,     pven,                          &
           ktrac,                                                       &
           pxtenh,   pxten,    pxtu,     pmfuxt,                        &
           zpmfun,   & ! op_ck_20031001
           pten,     pqen,     pqsen,                                   &
           pgeo,     pgeoh,    paphp1,                                  & 
           pqte,     pverv,    klwmin,                                  &
           ldcum,    ldland,   ktype,    klab,                          &
           ptu,      pqu,      plu,      puu,      pvu,                 &
           pmfu,     pmfub,    pentr,                                   &
           pmfus,    pmfuq,                                             &
           pmful,    plude,    pqude,    pdmfup,                        &
           khmin,    phhatt,   phcbase,  pqsenh,                        &
           pcpen,    pcpcu,                                             &
           kcbot,    kctop,    kctop0)
!
!          THIS ROUTINE DOES THE CALCULATIONS FOR CLOUD ASCENTS
!          FOR CUMULUS PARAMETERIZATION
!
!          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
!
!          PURPOSE.
!          --------
!          TO PRODUCE CLOUD ASCENTS FOR CU-PARAMETRIZATION
!          (VERTICAL PROFILES OF T,Q,L,U AND V AND CORRESPONDING
!           FLUXES AS WELL AS PRECIPITATION RATES)
!
!          INTERFACE
!          ---------
!
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!
!          METHOD.
!          --------
!          LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD BASE
!          AND THEN CALCULATE MOIST ASCENT FOR
!          ENTRAINING/DETRAINING PLUME.
!          ENTRAINMENT AND DETRAINMENT RATES DIFFER FOR
!          SHALLOW AND DEEP CUMULUS CONVECTION.
!          IN CASE THERE IS NO PENETRATIVE OR SHALLOW CONVECTION
!          CHECK FOR POSSIBILITY OF MID LEVEL CONVECTION
!          (CLOUD BASE VALUES CALCULATED IN *CUBASMC*)
!
!          EXTERNALS
!          ---------
!          *CUADJTQ* ADJUST T AND Q DUE TO CONDENSATION IN ASCENT
!          *CUENTR*  CALCULATE ENTRAINMENT/DETRAINMENT RATES
!          *CUBASMC* CALCULATE CLOUD BASE VALUES FOR MIDLEVEL CONVECTION
!          REFERENCE
!          ---------
!          (TIEDTKE,1989)
!
USE mo_kind,         ONLY : dp
USE mo_control,      ONLY : nn
USE mo_constants,    ONLY : g, tmelt, vtmpc1, rv, rd, alv, als
USE mo_cumulus_flux, ONLY : lmfdudv, lmfmid, nmctop, cmfcmin, cprcon   &
                          , cmfctop, centrmax
USE mo_time_control, ONLY : time_step_len
#ifdef __ibm__
USE mo_specfun,      ONLY : merge
#endif

IMPLICIT NONE

INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, klevm1, ktrac
INTEGER :: jl, jk, jt, ik, icall, ikb, ikt

REAL(dp) :: ptenh(kbdim,klev),       pqenh(kbdim,klev),                &
            puen(kbdim,klev),        pven(kbdim,klev),                 &
            pten(kbdim,klev),        pqen(kbdim,klev),                 &
            pgeo(kbdim,klev),        pgeoh(kbdim,klev),                &
            paphp1(kbdim,klevp1),                                      &
            pqsen(kbdim,klev),       pqte(kbdim,klev),                 &
            pverv(kbdim,klev)
!
REAL(dp) :: ptu(kbdim,klev),         pqu(kbdim,klev),                  &
            puu(kbdim,klev),         pvu(kbdim,klev),                  &
            pmfu(kbdim,klev),                                          &
            pmfub(kbdim),            pentr(kbdim),                     &
            pmfus(kbdim,klev),       pmfuq(kbdim,klev),                &
            plu(kbdim,klev),         plude(kbdim,klev),                &
            pqude(kbdim,klev),                                         &
            pmful(kbdim,klev),       pdmfup(kbdim,klev)
REAL(dp) :: zpmfun(kbdim,klev) ! op_ck_20031001
REAL(dp) :: zzp ! op_ck_20031001
REAL(dp) :: pcpen(kbdim,klev),       pcpcu(kbdim,klev)
INTEGER  :: klwmin(kbdim),           ktype(kbdim),                     &
            klab(kbdim,klev),        kcbot(kbdim),                     &
            kctop(kbdim),            kctop0(kbdim)
INTEGER  :: khmin(kbdim)
REAL(dp) :: phhatt(kbdim,klev)
REAL(dp) :: phcbase(kbdim)
REAL(dp) :: pqsenh(kbdim,klev)
LOGICAL  :: ldcum(kbdim),            ldland(kbdim)
!
REAL(dp) :: zdmfen(kbdim),           zdmfde(kbdim),                    &
            zmfuu(kbdim),            zmfuv(kbdim),                     &
            zpbase(kbdim),           zqold(kbdim)
REAL(dp) :: zph(kbdim)
REAL(dp) :: zodetr(kbdim,klev)
REAL(dp) :: zoentr(kbdim,klev)
REAL(dp) :: zbuoy(kbdim)
LOGICAL  :: loflag(kbdim)
REAL(dp) :: pxtenh(kbdim,klev,ktrac),pxten(kbdim,klev,ktrac),          &
            pxtu(kbdim,klev,ktrac),  pmfuxt(kbdim,klev,ktrac)
REAL(dp) :: zcons2, ztglace, zmfmax, zfac, zmftest, zqeen, zseen       &
          , zqude, zmfusk, zmfuqk, zmfulk, zxteen, zxtude, zmfuxtk     &
          , zbuo, zdnoprc, zprcon, zlnew, zz, zdmfeu, zdmfdu, zzdmf    &
          , zdz, zdrodz, zdprho, zalvs, zmse, znevn, zodmax, zga, zdt  &
          , zscod, zqcod, zbuoyz, zscde, zdlev
!
!
!      INTRINSIC FUNCTIONS
INTRINSIC MAX, MIN, LOG
!
!
!----------------------------------------------------------------------
!
!*    1.           SPECIFY PARAMETERS
!                  ------------------
!
100 CONTINUE
  zcons2=1._dp/(g*time_step_len)
  ztglace=tmelt-13._dp
  zqold(1:kproma) = 0.0_dp
  IF(klev == 11) THEN
    IF(nn == 21) THEN
      zdlev=1.5E4_dp
    ELSE IF(nn == 31) THEN
      zdlev=2.0E4_dp
    ELSE
      zdlev=3.0E4_dp
    ENDIF
  ELSE
   zdlev=3.0E4_dp
  ENDIF
!
!
!----------------------------------------------------------------------
!
!     2.           SET DEFAULT VALUES
!                  ------------------
!
200 CONTINUE
  DO 210 jl=1,kproma
     zmfuu(jl)=0._dp
     zmfuv(jl)=0._dp
     IF(.NOT.ldcum(jl)) ktype(jl)=0
210 END DO
  DO 230 jk=1,klev
     DO 220 jl=1,kproma
        plu(jl,jk)=0._dp
        pmfu(jl,jk)=0._dp
        zpmfun(jl,jk)=0._dp ! op_ck_20031001
        pmfus(jl,jk)=0._dp
        pmfuq(jl,jk)=0._dp
        pmful(jl,jk)=0._dp
        plude(jl,jk)=0._dp
        pqude(jl,jk)=0._dp
        pdmfup(jl,jk)=0._dp
        IF(.NOT.ldcum(jl).OR.ktype(jl).EQ.3) klab(jl,jk)=0
        IF(.NOT.ldcum(jl).AND.paphp1(jl,jk).LT.4.e4_dp) kctop0(jl)=jk
        IF(jk.LT.kcbot(jl)) klab(jl,jk)=0
220  END DO
     DO 2204 jt=1,ktrac
        DO 2202 jl=1,kproma
           pmfuxt(jl,jk,jt)=0._dp
2202    END DO
2204 END DO
!
230 END DO
  DO jk=1,klev
     DO jl=1,kproma
        zoentr(jl,jk)=0._dp
        zodetr(jl,jk)=0._dp
     ENDDO
  ENDDO
!
!
!----------------------------------------------------------------------
!
!     3.0          INITIALIZE VALUES AT LIFTING LEVEL
!                  ----------------------------------
!
300 CONTINUE
  DO 310 jl=1,kproma
     kctop(jl)=klevm1
     IF(.NOT.ldcum(jl)) THEN
        kcbot(jl)=klevm1
        pmfub(jl)=0._dp
        pqu(jl,klev)=0._dp
     END IF
     pmfu(jl,klev)=pmfub(jl)
     pmfus(jl,klev)=pmfub(jl)*(pcpcu(jl,klev)*ptu(jl,klev)             &
                                       +pgeoh(jl,klev))
     pmfuq(jl,klev)=pmfub(jl)*pqu(jl,klev)
     IF(lmfdudv) THEN
        zmfuu(jl)=pmfub(jl)*puu(jl,klev)
        zmfuv(jl)=pmfub(jl)*pvu(jl,klev)
     END IF
310 END DO
!
! op_ck_20031001+
  DO 3101 jk=1,klev
!DIR$ IVDEP
     DO 3102 jl=1,kproma
!
! LINEARISATION OF MASS FLUX FROM CLOUD BASE TO THE SURFACE:
! NEW VARIABLE FOR MASS FLUX (TRACER ONLY AT THIS STAGE): ZPMFUN
!
        IF (ktype(jl).ne.3) then
          IF (jk.eq.kcbot(jl)) zpmfun(jl,jk)=pmfub(jl)
          IF (jk.gt.kcbot(jl)) then
            ikb=kcbot(jl)
            zzp=(paphp1(jl,klevp1)-paphp1(jl,jk))/                &
                (paphp1(jl,klevp1)-paphp1(jl,ikb))
           zpmfun(jl,jk)=zpmfun(jl,ikb)*zzp
          END IF
        END IF
3102   END DO
3101 END DO
!
  DO 311 jl=1,kproma
     pmfu(jl,klev)=pmfub(jl)
     pmfus(jl,klev)=pmfub(jl)*(pcpcu(jl,klev)*ptu(jl,klev)          &
                                       +pgeoh(jl,klev))
     pmfuq(jl,klev)=pmfub(jl)*pqu(jl,klev)
     IF(lmfdudv) THEN
        zmfuu(jl)=pmfub(jl)*puu(jl,klev)
        zmfuv(jl)=pmfub(jl)*pvu(jl,klev)
     END IF
311 END DO
!
! op_ck_20031001-
  DO 3112 jt=1,ktrac
     DO 3110 jl=1,kproma
        IF(.NOT.ldcum(jl)) THEN
           pxtu(jl,klev,jt)=0._dp
        ENDIF
!       pmfuxt(jl,klev,jt)=pmfub(jl)*pxtu(jl,klev,jt)       ! op_ck_20031001
        pmfuxt(jl,klev,jt)=zpmfun(jl,klev)*pxtu(jl,klev,jt) ! op_ck_20031001
3110 END DO
3112 END DO
!
  DO 320 jl=1,kproma
     ldcum(jl)=.FALSE.
320 END DO
!
!
!
!----------------------------------------------------------------------
!
!     3.5          FIND ORGANIZED ENTRAINMENT AT CLOUD BASE
!                  ----------------------------------------
!
350 CONTINUE
  DO jl=1,kproma
     IF(ktype(jl).EQ.1) THEN
        ikb=kcbot(jl)
        zbuoy(jl)=g*(ptu(jl,ikb)-ptenh(jl,ikb))/ptenh(jl,ikb) +        &
                          g*vtmpc1*(pqu(jl,ikb)-pqenh(jl,ikb))
        IF(zbuoy(jl).GT.0._dp) THEN
           zdz=(pgeo(jl,ikb-1)-pgeo(jl,ikb))/g
           zdrodz=-LOG(pten(jl,ikb-1)/pten(jl,ikb))/zdz                &
                     -g/(rd*ptenh(jl,ikb)*(1._dp+vtmpc1*pqenh(jl,ikb)))
! nb zoentr is here a fractional value
           zoentr(jl,ikb-1)=zbuoy(jl)*0.5_dp/(1._dp+zbuoy(jl)*zdz)     &
                                                              + zdrodz
           zoentr(jl,ikb-1)=MIN(zoentr(jl,ikb-1),centrmax)
           zoentr(jl,ikb-1)=MAX(zoentr(jl,ikb-1),0._dp)
        ENDIF
     ENDIF
  ENDDO
!
!
!----------------------------------------------------------------------
!
!     4.           DO ASCENT: SUBCLOUD LAYER (KLAB=1) ,CLOUDS (KLAB=2)
!                  BY DOING FIRST DRY-ADIABATIC ASCENT AND THEN
!                  BY ADJUSTING T,Q AND L ACCORDINGLY IN *CUADJTQ*,
!                  THEN CHECK FOR BUOYANCY AND SET FLAGS ACCORDINGLY
!                  -------------------------------------------------
!
400 CONTINUE
  DO 480 jk=klevm1,2,-1
!
!                  SPECIFY CLOUD BASE VALUES FOR MIDLEVEL CONVECTION
!                  IN *CUBASMC* IN CASE THERE IS NOT ALREADY CONVECTION
!                  ----------------------------------------------------
!
     ik=jk
     IF(lmfmid.AND.ik.LT.klevm1.AND.ik.GT.nmctop) THEN
        CALL cubasmc(kproma, kbdim, klev, ik, klab,                     &
                     pten,     pqen,     pqsen,    puen,     pven,      &
                     ktrac,                                             &
                     paphp1,   zpmfun,   & ! op_ck_20031001
                     pxten,    pxtu,     pmfuxt,                        &
                     pverv,    pgeo,     pgeoh,    ldcum,    ktype,     &
                     pmfu,     pmfub,    pentr,    kcbot,               &
                     ptu,      pqu,      plu,      puu,      pvu,       &
                     pmfus,    pmfuq,    pmful,    pdmfup,   zmfuu,     &
                     pcpen,                                             &
                     zmfuv)
     ENDIF
!
     DO 410 jl=1,kproma
        IF(klab(jl,jk+1).EQ.0) klab(jl,jk)=0
        loflag(jl)=klab(jl,jk+1).GT.0
        zph(jl)=paphp1(jl,jk)
        IF(ktype(jl).EQ.3.AND.jk.EQ.kcbot(jl)) THEN
           zmfmax=(paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2
           IF(pmfub(jl).GT.zmfmax) THEN
              zfac=zmfmax/pmfub(jl)
              pmfu(jl,jk+1)=pmfu(jl,jk+1)*zfac
              pmfus(jl,jk+1)=pmfus(jl,jk+1)*zfac
              pmfuq(jl,jk+1)=pmfuq(jl,jk+1)*zfac
              zmfuu(jl)=zmfuu(jl)*zfac
              zmfuv(jl)=zmfuv(jl)*zfac
           END IF
        END IF
410  END DO
     DO 4102 jt=1,ktrac
        DO 4101 jl=1,kproma
           IF(ktype(jl).EQ.3.AND.jk.EQ.kcbot(jl)) THEN
              zmfmax=(paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2
              IF(pmfub(jl).GT.zmfmax) THEN
                 zfac=zmfmax/pmfub(jl)
                 pmfuxt(jl,jk+1,jt)=pmfuxt(jl,jk+1,jt)*zfac
              END IF
           END IF
4101    END DO
4102 END DO
!
! RESET PMFUB IF NECESSARY
!
     DO 4103 jl=1,kproma
        IF(ktype(jl).EQ.3.AND.jk.EQ.kcbot(jl)) THEN
           zmfmax=(paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2
           pmfub(jl)=MIN(pmfub(jl),zmfmax)
        END IF
4103 END DO
!
!
!*                 SPECIFY TURBULENT ENTRAINMENT AND DETRAINMENTS
!                  RATES PLUS ORGANIZED DETRAINMENT RATES IN *CUENTR*
!                   -------------------------------------
!
     ik=jk
     CALL cuentr(    kproma, kbdim, klev, klevp1, ik,                  &
           ptenh,    pqenh,    pqte,     paphp1,                       &
           klwmin,   ldcum,    ktype,    kcbot,    kctop0,             &
           zpbase,   pmfu,     pentr,    zodetr,                       &
           khmin,    pgeoh,                                            &
           zdmfen,   zdmfde)
!
!
!
!                  DO ADIABATIC ASCENT FOR ENTRAINING/DETRAINING PLUME
!                  THE CLOUD ENSEMBLE ENTRAINS ENVIRONMENTAL VALUES
!                  IN TURBULENT DETRAINMENT CLOUD ENSEMBLE VALUES
!                  ARE DETRAINED
!                  IN ORGANIZED DETRAINMENT THE DRY STATIC ENERGY AND
!                  MOISTURE THAT ARE NEUTRAL COMPARED TO THE
!                  ENVIRONMENTAL AIR ARE DETRAINED
!                  ---------------------------------------------------
!
     DO 420 jl=1,kproma
        IF(loflag(jl)) THEN
           IF(jk.LT.kcbot(jl)) THEN
              zmftest=pmfu(jl,jk+1)+zdmfen(jl)-zdmfde(jl)
              zmfmax=MIN(zmftest,(paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2)
              zdmfen(jl)=MAX(zdmfen(jl)-MAX(zmftest-zmfmax,0._dp),0._dp)
           END IF
           zdmfde(jl)=MIN(zdmfde(jl),0.75_dp*pmfu(jl,jk+1))
           pmfu(jl,jk)=pmfu(jl,jk+1)+zdmfen(jl)-zdmfde(jl)
           IF (ktype(jl).EQ.1 .AND. jk.LT.kcbot(jl)) THEN
              zdprho=(pgeoh(jl,jk)-pgeoh(jl,jk+1))/g
              zoentr(jl,jk)=zoentr(jl,jk)*zdprho*pmfu(jl,jk+1)
              zmftest=pmfu(jl,jk)+zoentr(jl,jk)-zodetr(jl,jk)
              zmfmax=MIN(zmftest,(paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2)
              zoentr(jl,jk)=MAX(zoentr(jl,jk)                          &
                                      -MAX(zmftest-zmfmax,0._dp),0._dp)
           ELSE
              zoentr(jl,jk)=0._dp
           ENDIF
           IF(ktype(jl).EQ.1.AND.jk.LT.kcbot(jl).AND.jk.LE.khmin(jl))  &
                                                                   THEN
!          limit organized detrainment to not allowing for too
!          deep clouds
              zalvs=MERGE(alv,als,ptu(jl,jk+1)>tmelt)
              zmse=pcpcu(jl,jk+1)*ptu(jl,jk+1)+zalvs*pqu(jl,jk+1)      &
                                                 +pgeoh(jl,jk+1)
              ikt=kctop0(jl)
              znevn=(pgeoh(jl,ikt)-pgeoh(jl,jk+1))                     &
                                              *(zmse-phhatt(jl,jk+1))/g
              IF(znevn.LE.0._dp) znevn=1.
              zdprho=(pgeoh(jl,jk)-pgeoh(jl,jk+1))/g
              zodmax=((phcbase(jl)-zmse)/znevn)*zdprho*pmfu(jl,jk+1)
              zodmax=MAX(zodmax,0._dp)
              zodetr(jl,jk)=MIN(zodetr(jl,jk),zodmax)
           ENDIF
           zodetr(jl,jk)=MIN(zodetr(jl,jk),0.75_dp*pmfu(jl,jk))
           pmfu(jl,jk)=pmfu(jl,jk)+zoentr(jl,jk)-zodetr(jl,jk)
           zqeen=pqenh(jl,jk+1)*zdmfen(jl)
           zqeen=zqeen+pqenh(jl,jk+1)*zoentr(jl,jk)
           zseen=(pcpcu(jl,jk+1)*ptenh(jl,jk+1)+pgeoh(jl,jk+1))        &
                                             *zdmfen(jl)
           zseen=zseen+(pcpcu(jl,jk+1)*ptenh(jl,jk+1)+pgeoh(jl,jk+1))  &
                                               *zoentr(jl,jk)
           zscde=(pcpcu(jl,jk+1)*ptu(jl,jk+1)+pgeoh(jl,jk+1))*zdmfde(jl)
! find moist static energy that give nonbuoyant air
           zalvs=MERGE(alv,als,ptenh(jl,jk+1)>tmelt)
           zga=zalvs*pqsenh(jl,jk+1)/(rv*(ptenh(jl,jk+1)**2))
           zdt=(plu(jl,jk+1)-vtmpc1*(pqsenh(jl,jk+1)-pqenh(jl,jk+1)))/ &
                          (1._dp/ptenh(jl,jk+1) + vtmpc1*zga)
           zscod=pcpcu(jl,jk+1)*ptenh(jl,jk+1)+pgeoh(jl,jk+1)          &
                                              +pcpcu(jl,jk+1)*zdt
           zscod=MAX(zscod,0._dp)
           zscde=zscde+zodetr(jl,jk)*zscod
           zqude=pqu(jl,jk+1)*zdmfde(jl)
           zqcod=pqsenh(jl,jk+1)+zga*zdt
           zqcod=MAX(zqcod,0._dp)
           zqude=zqude+zodetr(jl,jk)*zqcod
           pqude(jl,jk)=zqude
           plude(jl,jk)=plu(jl,jk+1)*zdmfde(jl)
           plude(jl,jk)=plude(jl,jk)+plu(jl,jk+1)*zodetr(jl,jk)
           zmfusk=pmfus(jl,jk+1)+zseen-zscde
           zmfuqk=pmfuq(jl,jk+1)+zqeen-zqude
           zmfulk=pmful(jl,jk+1)    -plude(jl,jk)
           plu(jl,jk)=zmfulk*(1._dp/MAX(cmfcmin,pmfu(jl,jk)))
           pqu(jl,jk)=zmfuqk*(1._dp/MAX(cmfcmin,pmfu(jl,jk)))
           ptu(jl,jk)=(zmfusk*(1._dp/MAX(cmfcmin,pmfu(jl,jk)))-        &
                                pgeoh(jl,jk))/pcpcu(jl,jk)
           ptu(jl,jk)=MAX(100._dp,ptu(jl,jk))
           ptu(jl,jk)=MIN(400._dp,ptu(jl,jk))
           zqold(jl)=pqu(jl,jk)
        END IF
420  END DO
!
!
     DO 4204 jt=1,ktrac
        DO 4202 jl=1,kproma
           IF(loflag(jl)) THEN
! op_ck_20031001+
!              zxteen=pxtenh(jl,jk+1,jt)*(zdmfen(jl)+zoentr(jl,jk))
!              zxtude=pxtu(jl,jk+1,jt)*(zdmfde(jl)+zodetr(jl,jk))
!              zmfuxtk=pmfuxt(jl,jk+1,jt)+zxteen-zxtude
!              pxtu(jl,jk,jt)=zmfuxtk*(1._dp/MAX(cmfcmin,pmfu(jl,jk)))
!          
             IF (jk.ge.kcbot(jl)) then
               IF (ktype(jl).ne.3) then
                 pmfuxt(jl,jk,jt)=(pxtu(jl,jk+1,jt)*zpmfun(jl,jk+1)    &
                   + (zpmfun(jl,jk)-zpmfun(jl,jk+1))*pxten(jl,jk,jt))
                 pxtu(jl,jk,jt)=pmfuxt(jl,jk,jt)                       &
                               /max(cmfcmin,zpmfun(jl,jk))
               ELSE
                 zpmfun(jl,jk)=pmfu(jl,jk)
                 pmfuxt(jl,jk,jt)=pmfuxt(jl,jk+1,jt)
                 pxtu(jl,jk,jt)=pmfuxt(jl,jk,jt)                       &
                               /max(cmfcmin,pmfu(jl,jk))
               ENDIF
             ELSE
               zpmfun(jl,jk)=pmfu(jl,jk)
               zxteen=pxten(jl,jk,jt)*(zdmfen(jl)+zoentr(jl,jk))
               zmfuxtk=pmfuxt(jl,jk+1,jt)+zxteen
               pmfuxt(jl,jk,jt)=zmfuxtk/(1.+(zdmfde(jl)+zodetr(jl,jk)) &
                               /max(cmfcmin,pmfu(jl,jk)))
               pxtu(jl,jk,jt)=pmfuxt(jl,jk,jt)/max(cmfcmin,pmfu(jl,jk))
             ENDIF
! op_ck_20031001-
           ENDIF
4202    END DO
4204 END DO
!
!
!                  DO CORRECTIONS FOR MOIST ASCENT
!                  BY ADJUSTING T,Q AND L IN *CUADJTQ*
!                  -----------------------------------
!
     ik=jk
     icall=1
     CALL cuadjtq(kproma, kbdim, klev, ik,                             &
                  zph,      ptu,      pqu,      loflag,   icall)
!
!DIR$ IVDEP
!OCL NOVREC
     DO 440 jl=1,kproma
       IF(loflag(jl)) THEN
         IF (pqu(jl,jk).LT.zqold(jl)) THEN
           klab(jl,jk)=2
           plu(jl,jk)=plu(jl,jk)+zqold(jl)-pqu(jl,jk)
           zbuo=ptu(jl,jk)*(1._dp+vtmpc1*pqu(jl,jk)-plu(jl,jk))-       &
                    ptenh(jl,jk)*(1._dp+vtmpc1*pqenh(jl,jk))
           IF(klab(jl,jk+1).EQ.1) zbuo=zbuo+0.5_dp
           IF(zbuo.GT.0._dp.AND.pmfu(jl,jk).GE.0.01_dp*pmfub(jl).AND.  &
                       jk.GE.kctop0(jl)) THEN
             kctop(jl)=jk
             ldcum(jl)=.TRUE.
             zdnoprc=MERGE(zdlev,1.5e4_dp,ldland(jl))
             zprcon=MERGE(0._dp,cprcon,                                &
                                   zpbase(jl)-paphp1(jl,jk).LT.zdnoprc)
             zlnew=plu(jl,jk)/                                         &
                           (1._dp+zprcon*(pgeoh(jl,jk)-pgeoh(jl,jk+1)))
             pdmfup(jl,jk)=MAX(0._dp,(plu(jl,jk)-zlnew)*pmfu(jl,jk))
             plu(jl,jk)=zlnew
           ELSE
             zpmfun(jl,jk)=0._dp ! op_ck_20031001
             klab(jl,jk)=0
             pmfu(jl,jk)=0._dp
           END IF
         END IF
       END IF
440  END DO
     DO 455 jl=1,kproma
        IF(loflag(jl)) THEN
           pmful(jl,jk)=plu(jl,jk)*pmfu(jl,jk)
           pmfus(jl,jk)=(pcpcu(jl,jk)*ptu(jl,jk)                       &
                                    +pgeoh(jl,jk))*pmfu(jl,jk)
           pmfuq(jl,jk)=pqu(jl,jk)*pmfu(jl,jk)
        END IF
455  END DO
     DO 4554 jt=1,ktrac
        DO 4552 jl=1,kproma
           IF(loflag(jl)) THEN
!             pmfuxt(jl,jk,jt)=pxtu(jl,jk,jt)*pmfu(jl,jk)   ! op_ck_20031001
              pmfuxt(jl,jk,jt)=pxtu(jl,jk,jt)*zpmfun(jl,jk) ! op_ck_20031001
           ENDIF
4552    END DO
4554 END DO
!
     IF(lmfdudv) THEN
        DO jl=1,kproma
           zdmfen(jl)=zdmfen(jl)+zoentr(jl,jk)
           zdmfde(jl)=zdmfde(jl)+zodetr(jl,jk)
        ENDDO
        DO 460 jl=1,kproma
           IF(loflag(jl)) THEN
              IF(ktype(jl).EQ.1.OR.ktype(jl).EQ.3) THEN
                 zz=MERGE(3._dp,2._dp,zdmfen(jl).EQ.0._dp)
              ELSE
                 zz=MERGE(1._dp,0._dp,zdmfen(jl).EQ.0._dp)
              END IF
              zdmfeu=zdmfen(jl)+zz*zdmfde(jl)
              zdmfdu=zdmfde(jl)+zz*zdmfde(jl)
              zdmfdu=MIN(zdmfdu,0.75_dp*pmfu(jl,jk+1))
              zmfuu(jl)=zmfuu(jl)+                                     &
                             zdmfeu*puen(jl,jk)-zdmfdu*puu(jl,jk+1)
              zmfuv(jl)=zmfuv(jl)+                                     &
                             zdmfeu*pven(jl,jk)-zdmfdu*pvu(jl,jk+1)
              IF(pmfu(jl,jk).GT.0._dp) THEN
                 puu(jl,jk)=zmfuu(jl)*(1._dp/pmfu(jl,jk))
                 pvu(jl,jk)=zmfuv(jl)*(1._dp/pmfu(jl,jk))
              END IF
           END IF
460     END DO
     END IF
!
!
!
!                  COMPUTE ORGANIZED ENTRAINMENT
!                  FOR USE AT NEXT LEVEL
!                  ------------------------------
!
     DO jl=1,kproma
        IF(loflag(jl).AND.ktype(jl).EQ.1) THEN
           zbuoyz=g*(ptu(jl,jk)-ptenh(jl,jk))/ptenh(jl,jk) +           &
                        g*vtmpc1*(pqu(jl,jk)-pqenh(jl,jk))-g*plu(jl,jk)
           zbuoyz=MAX(zbuoyz,0.0_dp)
           zdz=(pgeo(jl,jk-1)-pgeo(jl,jk))/g
           zdrodz=-LOG(pten(jl,jk-1)/pten(jl,jk))/zdz                  &
                       -g/(rd*ptenh(jl,jk)*(1._dp+vtmpc1*pqenh(jl,jk)))
           zbuoy(jl)=zbuoy(jl)+zbuoyz*zdz
           zoentr(jl,jk-1)=zbuoyz*0.5_dp/(1._dp+zbuoy(jl)) + zdrodz
           zoentr(jl,jk-1)=MIN(zoentr(jl,jk-1),centrmax)
           zoentr(jl,jk-1)=MAX(zoentr(jl,jk-1),0._dp)
!
        ENDIF
     ENDDO
!
!
480 END DO
!
!
!----------------------------------------------------------------------
!
!     5.           DETERMINE CONVECTIVE FLUXES ABOVE NON-BUOYANCY LEVEL
!                  ----------------------------------------------------
!                  (NOTE: CLOUD VARIABLES LIKE T,Q AND L ARE NOT
!                         AFFECTED BY DETRAINMENT AND ARE ALREADY KNOWN
!                         FROM PREVIOUS CALCULATIONS ABOVE)
!
500 CONTINUE
  DO 510 jl=1,kproma
     IF(kctop(jl).EQ.klevm1) ldcum(jl)=.FALSE.
     kcbot(jl)=MAX(kcbot(jl),kctop(jl))
510 END DO
!DIR$ IVDEP
  DO 530 jl=1,kproma
     IF(ldcum(jl)) THEN
        jk=kctop(jl)-1
        zzdmf=cmfctop
        zdmfde(jl)=(1._dp-zzdmf)*pmfu(jl,jk+1)
        plude(jl,jk)=zdmfde(jl)*plu(jl,jk+1)
        pqude(jl,jk)=zdmfde(jl)*pqu(jl,jk+1)
        pmfu(jl,jk)=pmfu(jl,jk+1)-zdmfde(jl)
        pdmfup(jl,jk)=0._dp
        pmfus(jl,jk)=(pcpcu(jl,jk)*ptu(jl,jk)+pgeoh(jl,jk))*pmfu(jl,jk)
        pmfuq(jl,jk)=pqu(jl,jk)*pmfu(jl,jk)
        pmful(jl,jk)=plu(jl,jk)*pmfu(jl,jk)
        IF(jk.GE.2) THEN
           plude(jl,jk-1)=pmful(jl,jk)
           pqude(jl,jk-1)=pmfuq(jl,jk)
        ELSE
           plude(jl,jk)=plude(jl,jk)+pmful(jl,jk)
           pqude(jl,jk)=pqude(jl,jk)+pmfuq(jl,jk)
        END IF
        zpmfun(jl,jk)=pmfu(jl,jk) ! op_ck_20031001
     END IF
530 END DO
  DO 5312 jt=1,ktrac
     DO 5310 jl=1,kproma
        IF(ldcum(jl)) THEN
           jk=kctop(jl)-1
!          pmfuxt(jl,jk,jt)=pxtu(jl,jk,jt)*pmfu(jl,jk)   ! op_ck_20031001
           pmfuxt(jl,jk,jt)=zpmfun(jl,jk)*pxtu(jl,jk,jt) ! op_ck_20031001
        ENDIF
5310 END DO
5312 END DO
!
  IF(lmfdudv) THEN
!DIR$      IVDEP
     DO 540 jl=1,kproma
        IF(ldcum(jl)) THEN
           jk=kctop(jl)-1
           puu(jl,jk)=puu(jl,jk+1)
           pvu(jl,jk)=pvu(jl,jk+1)
        END IF
540  END DO
  END IF
!
  RETURN
END SUBROUTINE cuasc
