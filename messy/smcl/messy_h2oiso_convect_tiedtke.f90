MODULE MESSY_H2OISO_CONVECT_TIEDTKE

!
  USE messy_main_constants_mem,     ONLY: g, R_gas, M_air, T0, wp, &
       rd, rv, cpd=>cp_air, cpv, vtmpc1, vtmpc2, tmelt
  USE messy_convect_tiedtke_param
  USE messy_convect_mem

  USE messy_main_tools,             ONLY: jptlucu1, jptlucu2, tlucua, tlucub, tlucuc   
  USE messy_convect_tiedtke,        ONLY: lookupoverflow
!---wiso-code

  USE messy_main_tools_wiso,      ONLY: tnat, cwisomin

!---wiso-code-end



!!$#ifdef __ibm__
!!$  USE mo_specfun, ONLY: merge
!!$#endif

 
  IMPLICIT NONE
  INTRINSIC :: MAX, MIN, SQRT, log, ABS, NINT, SUM, TAN, MERGE, TINY, REAL

  SAVE

! H2O related constants, liquid density, phase change constants
  REAL(dp),PARAMETER :: alv   = 2.5008e6_dp ! latent heat for vaporisation in J/kg
  REAL(dp),PARAMETER :: als   = 2.8345e6_dp ! latent heat for sublimation in J/kg
  REAL(dp),PARAMETER :: alf   = als - alv   

  LOGICAL :: lpos_def                                  ! switch for brinkop update for positive definite tracers

! Pointer for evapcu evaporation after Kuo
!  REAL(dp), ALLOCATABLE, DIMENSION(:):: cevapcu

!  LOGICAL :: lookupoverflow


  PUBLIC:: h2oiso_cuini


!======================================================================
CONTAINS
!======================================================================




SUBROUTINE h2oiso_cuini(kproma, kbdim, klev, klevp1, klevm1,                  &
           pten,     pqen,     pqsen,    pxen,     puen,     pven,     &
           ptven, &!re,    ktrac,                                            &
!re           pxten,    pxtenh,   pxtu,     pxtd,     pmfuxt,   pmfdxt,   &
           pverv,    pgeo,     paphp1,   pgeoh,                        &
           ptenh,    pqenh,    pqsenh,   pxenh,    klwmin,             &
           ptu,      pqu,      ptd,      pqd,                          &
           puu,      pvu,      pud,      pvd,                          &
           pmfu,     pmfd,     pmfus,    pmfds,                        &
           pmfuq,    pmfdq,    pdmfup,   pdmfdp,                       &
           pcpen,    pcpcu,                                            &
           pdpmel,   plu,      plude,    pqude,    klab,               &
!---wiso-code
           kwiso,                                                      &
           pwisoqen, pwisoxen,                                         &
           pwisoqenh,pwisoxenh,                                        &
           pwisoqu,  pwisoqd,                                          &
           pwisomfuq,pwisomfdq,pwisodmfup,pwisodmfdp,                  &
           pwisolu,  pwisolude,pwisoqude                               &
           )         
!---wiso-code-end
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


IMPLICIT NONE

INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, klevm1!re, ktrac

!---wiso-code

INTEGER, INTENT (IN) :: kwiso

!---wiso-code-end

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
REAL(dp):: zph(kbdim)
LOGICAL :: loflag(kbdim)
!re REAL(dp):: pxten(kbdim,klev,ktrac),   pxtenh(kbdim,klev,ktrac),        &
!re           pxtu(kbdim,klev,ktrac),    pxtd(kbdim,klev,ktrac),          &
!re           pmfuxt(kbdim,klev,ktrac),  pmfdxt(kbdim,klev,ktrac)

!---wiso-code

REAL(dp):: pwisoqen(kbdim,klev,kwiso) ,                                &
           pwisoxenh(kbdim,klev,kwiso), pwisoxen(kbdim,klev,kwiso),    &
           pwisoqenh(kbdim,klev,kwiso)
!
REAL(dp):: pwisoqu(kbdim,klev,kwiso),                                  &
           pwisoqd(kbdim,klev,kwiso),                                  &
           pwisomfuq(kbdim,klev,kwiso), pwisomfdq(kbdim,klev,kwiso),   &
           pwisodmfup(kbdim,klev,kwiso),pwisodmfdp(kbdim,klev,kwiso),  &
           pwisolu(kbdim,klev,kwiso),   pwisolude(kbdim,klev,kwiso),   &
           pwisoqude(kbdim,klev,kwiso)

!---wiso-code-end

INTEGER :: jk, jl, jt, ik, icall
REAL(dp):: zarg, zcpm, zzs
!
!  INTRINSIC FUNCTIONS
INTRINSIC MAX, MIN

!---wiso-code

  REAL(dp)            :: zdelta

!---wiso-code-end

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
!re     DO 1104 jt=1,ktrac
!re        DO 1102 jl=1,kproma
!re           pxtenh(jl,jk,jt)=(pxten(jl,jk,jt)+pxten(jl,jk-1,jt))*0.5_dp
!re1102    END DO
!re1104 END DO
!
!
     ik=jk
     icall=0
     CALL h2oiso_cuadjtq(kproma, kbdim, klev, ik,                             &
                  zph,      ptenh,    pqsenh,   loflag,   icall)
     IF (lookupoverflow) RETURN
!
     DO 120 jl=1,kproma
        pxenh(jl,jk)=(pxen(jl,jk)+pxen(jl,jk-1))*0.5_dp
        pqenh(jl,jk)=MIN(pqen(jl,jk-1),pqsen(jl,jk-1))                 &
                          +(pqsenh(jl,jk)-pqsen(jl,jk-1))
        pqenh(jl,jk)=MAX(pqenh(jl,jk),0._dp)
!        pcpcu(jl,jk)=cpd*(1.+vtmpc2*pqenh(jl,jk))
        pcpcu(jl,jk)=cpd
120  END DO

!---wiso-code

! Half level values of water isotopes in updrafts.
! Arith. Weighted Mean of the relations (Water isotope/normal Water) for
! and arith. Mean of the absolute values for cloud water

! do not use negative "normal water" values for delta calculation unless
! water on both levels is negative

     DO jt=1,kwiso
       DO jl=1,kproma
         pwisoxenh(jl,jk,jt)=(pwisoxen(jl,jk,jt)+pwisoxen(jl,jk-1,jt))*0.5_dp
         IF ((abs(pqen(jl,jk-1)).lt.cwisomin).and.(abs(pqen(jl,jk)).lt.cwisomin)) THEN
           zdelta=tnat(jt)    
         ELSEIF ((pqen(jl,jk-1).le.0).and.(pqen(jl,jk).gt.0)) THEN
           zdelta=MIN(pwisoqen(jl,jk,jt)/pqen(jl,jk),1.0_dp)!remin
         ELSEIF ((pqen(jl,jk-1).gt.0).and.(pqen(jl,jk).le.0)) THEN
           zdelta=MIN(pwisoqen(jl,jk-1,jt)/pqen(jl,jk-1),1.0_dp)!remin
         ELSE
           zdelta=MIN((pwisoqen(jl,jk,jt)+pwisoqen(jl,jk-1,jt))/(pqen(jl,jk)+pqen(jl,jk-1)),1.0_dp)!remin
         ENDIF
IF(jt.eq.1)zdelta=1.0_dp!rehc
         pwisoqenh(jl,jk,jt)=zdelta*pqenh(jl,jk)
       END DO
     END DO

!---wiso-code-end

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
!---wiso-code

     DO jt=1,kwiso
       DO jl=1,kproma
         pwisoxenh(jl,klev,jt)=pwisoxen(jl,klev,jt)
         pwisoqenh(jl,klev,jt)=pwisoqen(jl,klev,jt)
         pwisoxenh(jl,1,jt)=pwisoxen(jl,1,jt)
         pwisoqenh(jl,1,jt)=pwisoqen(jl,1,jt)   
       END DO
     END DO

!---wiso-code-end

!re  DO 1404 jt=1,ktrac
!re     DO 1402 jl=1,kproma
!re        pxtenh(jl,klev,jt)=pxten(jl,klev,jt)
!re        pxtenh(jl,1,jt)=pxten(jl,1,jt)
!re1402 END DO
!re1404 END DO
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
!---wiso-code

     DO jt=1,kwiso
       DO jl=1,kproma
         pwisoqu(jl,jk,jt)=pwisoqenh(jl,jk,jt)
         pwisoqd(jl,jk,jt)=pwisoqenh(jl,jk,jt)
         pwisolu(jl,jk,jt)=0._dp
         pwisomfuq(jl,jk,jt)=0._dp
         pwisomfdq(jl,jk,jt)=0._dp
         pwisodmfup(jl,jk,jt)=0._dp
         pwisodmfdp(jl,jk,jt)=0._dp
         pwisolude(jl,jk,jt)=0._dp
         pwisoqude(jl,jk,jt)=0._dp
       END DO
     END DO

!---wiso-code-end

!re     DO 2204 jt=1,ktrac
!re        DO 2202 jl=1,kproma
!re           pxtu(jl,jk,jt)=pxtenh(jl,jk,jt)
!re           pxtd(jl,jk,jt)=pxtenh(jl,jk,jt)
!re           pmfuxt(jl,jk,jt)=0._dp
!re           pmfdxt(jl,jk,jt)=0._dp
!re2202    END DO
!re2204 END DO
!
230 END DO
!
  RETURN
END SUBROUTINE h2oiso_cuini


!------------------------------------------------------------------------------
!==============================================================================
!------------------------------------------------------------------------------


SUBROUTINE h2oiso_cubase(   kproma, kbdim, klev, klevp1, klevm1,              &
           ptenh,    pqenh,    pgeoh,    paph,                         &
           ptu,      pqu,      plu,                                    &
           puen,     pven,     puu,      pvu,                          &
           pcpcu,                                                      &
           ldcum,    kcbot,    klab,                                   &
!---wiso-code
           kwiso,                                                      &
           pwisoqu,  pwisolu                                           &
           )
!---wiso-code-end
!
!          THIS ROUTINE CALCULATES CLOUD BASE VALUES (T AND Q)
!          FOR CUMULUS PARAMETERIZATION
!
!          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
!
!          PURPOSE.
!          --------
!          TO PRODUCE CLOUD BASE VALUES FOR CU-PARAMETRIZATION
!
!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!          INPUT ARE ENVIRONM. VALUES OF T,Q,P,PHI AT HALF LEVELS.
!          IT RETURNS CLOUD BASE VALUES AND FLAGS AS FOLLOWS;
!                 KLAB=1 FOR SUBCLOUD LEVELS
!                 KLAB=2 FOR CONDENSATION LEVEL
!
!          METHOD.
!          --------
!          LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD BASE
!          (NON ENTRAINING PLUME,I.E.CONSTANT MASSFLUX)
!
!          EXTERNALS
!          ---------
!          *CUADJTQ* FOR ADJUSTING T AND Q DUE TO CONDENSATION IN ASCENT
!

IMPLICIT NONE

INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, klevm1

!---wiso-code

INTEGER, INTENT (IN) :: kwiso

!---wiso-code-end

REAL(dp):: ptenh(kbdim,klev),       pqenh(kbdim,klev),                 &
           pgeoh(kbdim,klev),       paph(kbdim,klevp1)
!
REAL(dp):: ptu(kbdim,klev),         pqu(kbdim,klev),                   &
           plu(kbdim,klev)
REAL(dp):: puen(kbdim,klev),        pven(kbdim,klev),                  &
           puu(kbdim,klev),         pvu(kbdim,klev)
REAL(dp):: pcpcu(kbdim,klev)
INTEGER :: klab(kbdim,klev),        kcbot(kbdim)
LOGICAL :: ldcum(kbdim)
!
!---wiso-code

REAL(dp):: pwisoqu(kbdim,klev,kwiso),                                  &
           pwisolu(kbdim,klev,kwiso)

!---wiso-code-end

REAL(dp):: zqold(kbdim)
REAL(dp):: zph(kbdim)
LOGICAL :: loflag(kbdim)

INTEGER :: jl, jk, is, ik, ikb, icall
REAL(dp):: zbuo, zz

!---wiso-code

REAL(dp):: zqcond(kbdim,klev),zwisoqcond(kbdim,klev,kwiso),zqu_tmp(kbdim)
INTEGER :: jt

!---wiso-code-end

!
!
!----------------------------------------------------------------------
!
!     1.           INITIALIZE VALUES AT LIFTING LEVEL
!                  ----------------------------------
!
100 CONTINUE
  DO 110 jl=1,kproma
     klab(jl,klev)=1
     kcbot(jl)=klevm1
     ldcum(jl)=.FALSE.
     puu(jl,klev)=puen(jl,klev)*(paph(jl,klevp1)-paph(jl,klev))
     pvu(jl,klev)=pven(jl,klev)*(paph(jl,klevp1)-paph(jl,klev))
110 END DO
!
!
!----------------------------------------------------------------------
!
!     2.0          DO ASCENT IN SUBCLOUD LAYER,
!                  CHECK FOR EXISTENCE OF CONDENSATION LEVEL,
!                  ADJUST T,Q AND L ACCORDINGLY IN *CUADJTQ*,
!                  CHECK FOR BUOYANCY AND SET FLAGS
!                  -------------------------------------
!

!---wiso-code

  zqcond(:,:) = 0._dp                ! dummy condensate of normal water at cloud base 
  zwisoqcond(:,:,:) = 0._dp          ! dummy condensate of water isotope at cloud base

!---wiso-code-end

200 CONTINUE
  DO 290 jk=klevm1,2,-1
     is=0
     DO 210 jl=1,kproma
        is=is+MERGE(1,0,klab(jl,jk+1).EQ.1)
        loflag(jl)=klab(jl,jk+1).EQ.1
        zph(jl)=paph(jl,jk)
210  END DO
!     IF(is.EQ.0) go to 290 ! mz_ht_20071217
     DO 220 jl=1,kproma
        IF(loflag(jl)) THEN
           pqu(jl,jk)=pqu(jl,jk+1)
           ptu(jl,jk)=(pcpcu(jl,jk+1)*ptu(jl,jk+1)+pgeoh(jl,jk+1)      &
                           -pgeoh(jl,jk))/pcpcu(jl,jk)
           zbuo=ptu(jl,jk)*(1._dp+vtmpc1*pqu(jl,jk))-ptenh(jl,jk)      &
                                   *(1._dp+vtmpc1*pqenh(jl,jk))+0.5_dp
           IF(zbuo.GT.0._dp) klab(jl,jk)=1
           zqold(jl)=pqu(jl,jk)
        END IF
220  END DO
!
!---wiso-code

     DO jt=1,kwiso
       DO jl=1,kproma
         IF(loflag(jl)) THEN
           pwisoqu(jl,jk,jt)=pwisoqu(jl,jk+1,jt)
         END IF
       END DO
     END DO

!---wiso-code-end

     ik=jk
     icall=1
     zqu_tmp(:)=pqu(:,ik)                      !--- wiso-code: Store old value of pqu
     CALL h2oiso_cuadjtq(kproma, kbdim, klev, ik,                             &
                  zph,      ptu,      pqu,      loflag,   icall)
     IF (lookupoverflow) RETURN
!
!DIR$ IVDEP
!OCL NOVREC
     DO 240 jl=1,kproma
        zqcond(jl,jk)=zqu_tmp(jl)-pqu(jl,jk)   !--- wiso-code: Storing the condensate of normal water
        IF(loflag(jl).AND.pqu(jl,jk).LT.zqold(jl)) THEN
           klab(jl,jk)=2
           plu(jl,jk)=plu(jl,jk)+zqold(jl)-pqu(jl,jk)
           zbuo=ptu(jl,jk)*(1._dp+vtmpc1*pqu(jl,jk)-plu(jl,jk))-       &
                        ptenh(jl,jk)*(1._dp+vtmpc1*pqenh(jl,jk))+0.5_dp
           IF(zbuo.GT.0.) THEN
              kcbot(jl)=jk
              ldcum(jl)=.TRUE.
           END IF
        END IF
240  END DO

!---wiso-code

! Calculate the Condensate of water isotope (including fractionation)
! -> Subroutine cuadjwisoq returns new values of pwisoqu
    ik = jk
    CALL h2oiso_cuadjwisoq(kproma,   kbdim,      klev,      kwiso,    ik,     &
                    ptu,      pqu,        pwisoqu,   pwisolu,          &
                    zqcond,   zwisoqcond, loflag)

!---wiso-code-end
!
!             CALCULATE AVERAGES OF U AND V FOR SUBCLOUD ARA,.
!             THE VALUES WILL BE USED TO DEFINE CLOUD BASE VALUES.
!
     IF(lmfdudv) THEN
        DO 250 jl=1,kproma
           IF(jk.GE.kcbot(jl)) THEN
              puu(jl,klev)=puu(jl,klev)+                               &
                             puen(jl,jk)*(paph(jl,jk+1)-paph(jl,jk))
              pvu(jl,klev)=pvu(jl,klev)+                               &
                             pven(jl,jk)*(paph(jl,jk+1)-paph(jl,jk))
           END IF
250     END DO
     END IF
!
290 END DO
!
!
  IF(lmfdudv) THEN
     DO 310 jl=1,kproma
        IF(ldcum(jl)) THEN
           ikb=kcbot(jl)
           zz=1._dp/(paph(jl,klevp1)-paph(jl,ikb))
           puu(jl,klev)=puu(jl,klev)*zz
           pvu(jl,klev)=pvu(jl,klev)*zz
        ELSE
           puu(jl,klev)=puen(jl,klevm1)
           pvu(jl,klev)=pven(jl,klevm1)
        END IF
310  END DO
  END IF
!
  RETURN
END SUBROUTINE h2oiso_cubase

!======================================================================



SUBROUTINE h2oiso_cuasc(    kproma, kbdim, klev, klevp1, klevm1,              &
           ptenh,    pqenh,    puen,     pven,                         &
!re           ktrac,                                                      &
!re           pxtenh,   pxten,    pxtu,     pmfuxt,                       &
           pten,     pqen,     pqsen,                                  &
           pgeo,     pgeoh,    paphp1,                                 &
           pqte,     pverv,    klwmin,                                 &
           ldcum,    ldland,   ktype,    klab,                         &
           ptu,      pqu,      plu,      puu,      pvu,                &
           pmfu,     pmfub,    pentr,                                  &
           pmfus,    pmfuq,                                            &
           pmful,    plude,    pqude,    pdmfup,                       &
           khmin,    phhatt,   phcbase,  pqsenh,                       &
           pcpen,    pcpcu,                                            &
           kcbot,    kctop,    kctop0,                                 &
!---wiso-code
           kwiso,                                                      &
           pwisoqenh,                                                  &
           pwisoqen,                                                   &
           pwisoqu,  pwisolu,                                          &
           pwisomfuq,                                                  &
           pwisomful,pwisolude,pwisoqude,pwisodmfup                    &
           , nn, time_step_len)
!---wiso-code-end
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
!!$USE mo_constants,    ONLY : g, tmelt, vtmpc1, rv, rd, alv, als
!!$USE mo_cumulus_flux, ONLY : lmfdudv, lmfmid, nmctop, cmfcmin, cprcon   &
!!$                          , cmfctop, centrmax
!!$USE mo_time_control, ONLY : time_step_len

!---wiso-code

  USE messy_h2oiso,     ONLY : wiso_frac_liq_ice, fractcal
  USE messy_main_tools_wiso, ONLY: wiso_frac_liq_ice, tnat, cthomi &
       , cwisomin, cwisosec

!---wiso-code-end

IMPLICIT NONE

INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, klevm1!re, ktrac

REAL(DP),INTENT(IN):: time_step_len
INTEGER, INTENT (IN) :: nn
!---wiso-code

INTEGER, INTENT (IN) :: kwiso

!---wiso-code-end

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
!re REAL(dp) :: pxtenh(kbdim,klev,ktrac),pxten(kbdim,klev,ktrac),          &
!re            pxtu(kbdim,klev,ktrac),  pmfuxt(kbdim,klev,ktrac)
!
!---wiso-code

REAL(dp)::  pwisoqenh(kbdim,klev,kwiso),                               &
            pwisoqen(kbdim,klev,kwiso)

REAL(dp) :: pwisoqu(kbdim,klev,kwiso),                                 &
            pwisomfuq(kbdim,klev,kwiso),                               &
            pwisolu(kbdim,klev,kwiso),pwisolude(kbdim,klev,kwiso),     &
            pwisoqude(kbdim,klev,kwiso),                               &
            pwisomful(kbdim,klev,kwiso),pwisodmfup(kbdim,klev,kwiso)

!---wiso-code-end

REAL(dp) :: zcons2, ztglace, zmfmax, zfac, zmftest, zqeen, zseen       &
          , zqude, zmfusk, zmfuqk, zmfulk, zxteen, zxtude, zmfuxtk     &
          , zbuo, zdnoprc, zprcon, zlnew, zz, zdmfeu, zdmfdu, zzdmf    &
          , zdz, zdrodz, zdprho, zalvs, zmse, znevn, zodmax, zga, zdt  &
          , zscod, zqcod, zbuoyz, zscde, zdlev
!
!---wiso-code

REAL(dp) :: zprec_tmp(kbdim,klev),  zqcod_tmp(kbdim,klev)

REAL(dp) :: zwisofracliq(kbdim,kwiso), zwisofrac(kbdim,kwiso), zwisofracice(kbdim,kwiso)

REAL(dp) :: zwisoqeen, zwisoqcod,                                      &
            zwisoqude, zwisomfuqk, zwisomfulk

REAL(dp) :: zqliq, zqice,                                              &
            zql,  zqv,   zwisoql, zwisoqv,                             &
            zquot,zquot2,zqvo,    zdelta

LOGICAL  :: lo
            
!---wiso-code-end

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

!---wiso-code

     DO jt=1,kwiso
       DO jl=1,kproma
         pwisolu(jl,jk,jt)=0._dp
         pwisomfuq(jl,jk,jt)=0._dp
         pwisomful(jl,jk,jt)=0._dp
         pwisolude(jl,jk,jt)=0._dp
         pwisoqude(jl,jk,jt)=0._dp
         pwisodmfup(jl,jk,jt)=0._dp
       END DO
     END DO

!---wiso-code-end

!re     DO 2204 jt=1,ktrac
!re        DO 2202 jl=1,kproma
!re           pmfuxt(jl,jk,jt)=0._dp
!re2202    END DO
!re2204 END DO
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
!---wiso-code

  DO jt=1,kwiso
    DO jl=1,kproma
      IF(.NOT.ldcum(jl)) THEN
        pwisoqu(jl,klev,jt)=0._dp
      ENDIF
     pwisomfuq(jl,klev,jt)=pmfub(jl)*pwisoqu(jl,klev,jt)
    END DO
  END DO

!---wiso-code-end

!re  DO 3112 jt=1,ktrac
!re     DO 3110 jl=1,kproma
!re        IF(.NOT.ldcum(jl)) THEN
!re           pxtu(jl,klev,jt)=0._dp
!re        ENDIF
!re        pmfuxt(jl,klev,jt)=pmfub(jl)*pxtu(jl,klev,jt)
!re3110 END DO
!re3112 END DO
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
        CALL h2oiso_cubasmc(kproma, kbdim, klev, ik, klab,                    &
                     pten,     pqen,     pqsen,    puen,     pven,     &
!re                     ktrac,                                            &
!re                     pxten,    pxtu,     pmfuxt,                       &
                     pverv,    pgeo,     pgeoh,    ldcum,    ktype,    &
                     pmfu,     pmfub,    pentr,    kcbot,              &
                     ptu,      pqu,      plu,      puu,      pvu,      &
                     pmfus,    pmfuq,    pmful,    pdmfup,   zmfuu,    &
                     pcpen,                                            &
                     zmfuv,                                            &
!---wiso-code
                     kwiso,                                            &
                     pwisoqen,                                         &
                     pwisoqu,  pwisolu,                                &
                     pwisomfuq,pwisomful,pwisodmfup                    &
                     )
!---wiso-code-end
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

!---wiso-code

     DO jt=1,kwiso
       DO jl=1,kproma
         IF(ktype(jl).EQ.3.AND.jk.EQ.kcbot(jl)) THEN
           zmfmax=(paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2
           IF(pmfub(jl).GT.zmfmax) THEN
             zfac=zmfmax/pmfub(jl)
             pwisomfuq(jl,jk+1,jt)=pwisomfuq(jl,jk+1,jt)*zfac
           END IF
         END IF
       END DO
     END DO

!---wiso-code-end

!re     DO 4102 jt=1,ktrac
!re        DO 4101 jl=1,kproma
!re           IF(ktype(jl).EQ.3.AND.jk.EQ.kcbot(jl)) THEN
!re              zmfmax=(paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2
!re              IF(pmfub(jl).GT.zmfmax) THEN
!re                 zfac=zmfmax/pmfub(jl)
!re                 pmfuxt(jl,jk+1,jt)=pmfuxt(jl,jk+1,jt)*zfac
!re              END IF
!re           END IF
!re4101    END DO
!re4102 END DO
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
     CALL h2oiso_cuentr(    kproma, kbdim, klev, klevp1, ik,                  &
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
              IF(znevn.LE.0._dp) znevn=1._dp
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
           zqcod_tmp(jl,jk)=zqcod              !---wiso-code: store value of zqcod
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
!---wiso-code

     DO jt=1,kwiso
       DO jl=1,kproma
         IF(loflag(jl)) THEN
!re           zwisoqeen=pwisoqenh(jl,jk+1,jt)*(zdmfen(jl)+zoentr(jl,jk))
           zwisoqeen=pwisoqenh(jl,jk+1,jt)*zdmfen(jl)!re
           zwisoqeen=zwisoqeen+pwisoqenh(jl,jk+1,jt)*zoentr(jl,jk)!re

           zwisoqude=pwisoqu(jl,jk+1,jt)*zdmfde(jl)
           zdelta=tnat(jt)
           IF ((pwisoqenh(jl,jk+1,jt).GT.cwisomin).AND.&
                (pqenh(jl,jk+1).GT.cwisomin)) &
                zdelta=MIN(pwisoqenh(jl,jk+1,jt)/pqenh(jl,jk+1),1.0_dp)!remin
           IF (ABS(1._dp - zdelta).LT.cwisosec) zdelta = 1._dp
IF(jt.eq.1)zdelta=1.0_dp!rehc
           zwisoqcod=zqcod_tmp(jl,jk)*zdelta
           zwisoqude=zwisoqude+zodetr(jl,jk)*zwisoqcod
           pwisoqude(jl,jk,jt)=zwisoqude
           pwisolude(jl,jk,jt)=pwisolu(jl,jk+1,jt)*(zdmfde(jl)+zodetr(jl,jk))
           zwisomfuqk=pwisomfuq(jl,jk+1,jt)+zwisoqeen-zwisoqude
           zwisomfulk=pwisomful(jl,jk+1,jt)-pwisolude(jl,jk,jt)
           pwisolu(jl,jk,jt)=zwisomfulk*(1._dp/MAX(cmfcmin,pmfu(jl,jk)))
           pwisoqu(jl,jk,jt)=zwisomfuqk*(1._dp/MAX(cmfcmin,pmfu(jl,jk)))
         ENDIF
       END DO
     END DO

!---wiso-code-end

!re     DO 4204 jt=1,ktrac
!re        DO 4202 jl=1,kproma
!re           IF(loflag(jl)) THEN
!re              zxteen=pxtenh(jl,jk+1,jt)*(zdmfen(jl)+zoentr(jl,jk))
!re              zxtude=pxtu(jl,jk+1,jt)*(zdmfde(jl)+zodetr(jl,jk))
!re              zmfuxtk=pmfuxt(jl,jk+1,jt)+zxteen-zxtude
!re              pxtu(jl,jk,jt)=zmfuxtk*(1._dp/MAX(cmfcmin,pmfu(jl,jk)))
!re           ENDIF
!re4202    END DO
!re4204 END DO
!
!
!                  DO CORRECTIONS FOR MOIST ASCENT
!                  BY ADJUSTING T,Q AND L IN *CUADJTQ*
!                  -----------------------------------
!
     ik=jk
     icall=1
     CALL h2oiso_cuadjtq(kproma, kbdim, klev, ik,                             &
                  zph,      ptu,      pqu,      loflag,   icall)
!
!---wiso-code

! calculate fractionation factors for liquid-vapour and solid-vapour phase change

!reold     CALL fractcal(kproma,kbdim,kwiso,zwisofrac,zwisofracice,ptu(:,ik),tmelt)
     CALL wiso_frac_liq_ice(kproma,kbdim,kwiso,ptu(:,ik),zwisofracliq,zwisofracice)

!    adjusting isotope vapour *pwisoqu* depending on fract.
     DO jt=1,kwiso
       DO jl=1,kproma
         IF (loflag(jl)) THEN
! divide into ice and liquid
           IF (ptu(jl,jk).gt.tmelt) THEN
             zqliq=1._dp
             zqice=0._dp
           ELSEIF (ptu(jl,jk).lt.cthomi) THEN
             zqice=1._dp
             zqliq=0._dp
           ELSE
             zqice=(tmelt-ptu(jl,jk))/(tmelt-cthomi)
             zqliq=1._dp-zqice
           ENDIF
! fractionation of ice and liquid
           zql=zqliq*(plu(jl,jk)+zqold(jl)-pqu(jl,jk))
           zqv=zqold(jl)+zqliq*(pqu(jl,jk)-zqold(jl))
           zwisoql=zqliq*pwisolu(jl,jk,jt)
           zwisoqv=pwisoqu(jl,jk,jt)
!reold           zquot=zqv+zwisofrac(jl,jt)*zql
           zquot=zqv+zwisofracliq(jl,jt)*zql
           zqvo=zqv
           IF (zquot.gt.cwisomin.and.zqvo.gt.cwisomin) THEN
!reold             zdelta=(zwisofrac(jl,jt)*zql*zwisoqv-zwisoql*zqv)/zquot 
             zdelta=(zwisofracliq(jl,jt)*zql*zwisoqv-zwisoql*zqv)/zquot 
             zqv=pqu(jl,jk)
             zwisoqv=zwisoqv-zdelta
             zquot2=zqv/zqvo
             lo=abs(1._dp-zquot2).lt.cwisosec
             zquot2=MERGE(1._dp,zquot2,lo)
if(zquot2.gt.0.0_dp)then!rehc
             zdelta=zdelta+zwisoqv*(1._dp-zquot2**zwisofracice(jl,jt))
endif!rehc
             pwisoqu(jl,jk,jt)=pwisoqu(jl,jk,jt)-zdelta
             IF(loflag(jl).AND.pqu(jl,jk).LT.zqold(jl)) THEN  ! add condensate to liquid as it is done
               pwisolu(jl,jk,jt)=pwisolu(jl,jk,jt)+zdelta     ! for normal water in next loop 440
             END IF
           ENDIF
         ENDIF
       END DO
     END DO

!---wiso-code-end

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
!---wiso-code
! Store fraction of precipitation formed from cloud water
              IF(plu(jl,jk).gt.cwisomin) THEN
                zprec_tmp(jl,jk)=(plu(jl,jk)-zlnew)/plu(jl,jk)
              ENDIF
              IF (abs(1._dp-zprec_tmp(jl,jk)).lt.cwisosec) zprec_tmp(jl,jk)=1._dp
!---wiso-code-end
             plu(jl,jk)=zlnew
           ELSE
             klab(jl,jk)=0
             pmfu(jl,jk)=0._dp
!---wiso-code
! No precipitation formed from cloud water
             zprec_tmp(jl,jk)=0._dp
!---wiso-code-end
           END IF
         END IF
       END IF
440  END DO

!---wiso-code

     DO jt=1,kwiso
       DO jl=1,kproma
        IF(loflag(jl).AND.pqu(jl,jk).LT.zqold(jl)) THEN
! Calculate water isotope precipitation and new mass fluxes for the water isotopes in updrafts
           pwisodmfup(jl,jk,jt)=MAX(0._dp,pwisolu(jl,jk,jt)*pmfu(jl,jk)*zprec_tmp(jl,jk))
           pwisolu(jl,jk,jt)=(1._dp-zprec_tmp(jl,jk))*pwisolu(jl,jk,jt)
         END IF
       END DO
     END DO

!---wiso-code-end

     DO 455 jl=1,kproma
        IF(loflag(jl)) THEN
           pmful(jl,jk)=plu(jl,jk)*pmfu(jl,jk)
           pmfus(jl,jk)=(pcpcu(jl,jk)*ptu(jl,jk)                       &
                                    +pgeoh(jl,jk))*pmfu(jl,jk)
           pmfuq(jl,jk)=pqu(jl,jk)*pmfu(jl,jk)
        END IF
455  END DO

!---wiso-code

     DO jt=1,kwiso
       DO jl=1,kproma
         IF(loflag(jl)) THEN
           pwisomful(jl,jk,jt)=pwisolu(jl,jk,jt)*pmfu(jl,jk)
           pwisomfuq(jl,jk,jt)=pwisoqu(jl,jk,jt)*pmfu(jl,jk)
         END IF
       END DO
     END DO

!---wiso-code-end

!re     DO 4554 jt=1,ktrac
!re        DO 4552 jl=1,kproma
!re           IF(loflag(jl)) THEN
!re              pmfuxt(jl,jk,jt)=pxtu(jl,jk,jt)*pmfu(jl,jk)
!re           ENDIF
!re4552    END DO
!re4554 END DO
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
     END IF
530 END DO

!---wiso-code

! Calculate water isotope fluxes above the non-buoyancy level
     DO jt=1,kwiso
!DIR$ IVDEP
       DO jl=1,kproma
         IF (ldcum(jl)) THEN
           jk=kctop(jl)-1
           pwisolude(jl,jk,jt)=zdmfde(jl)*pwisolu(jl,jk+1,jt)
           pwisoqude(jl,jk,jt)=zdmfde(jl)*pwisoqu(jl,jk+1,jt)
           pwisodmfup(jl,jk,jt)=0._dp
           pwisomfuq(jl,jk,jt)=pwisoqu(jl,jk,jt)*pmfu(jl,jk)
           pwisomful(jl,jk,jt)=pwisolu(jl,jk,jt)*pmfu(jl,jk)
           IF(jk.GE.2) THEN
             pwisolude(jl,jk-1,jt)=pwisomful(jl,jk,jt)
             pwisoqude(jl,jk-1,jt)=pwisomfuq(jl,jk,jt)
           ELSE
             pwisolude(jl,jk,jt)=pwisolude(jl,jk,jt)+pwisomful(jl,jk,jt)
             pwisoqude(jl,jk,jt)=pwisoqude(jl,jk,jt)+pwisomfuq(jl,jk,jt)
           END IF
         END IF
       END DO
     END DO

!---wiso-code-end

!re  DO 5312 jt=1,ktrac
!re     DO 5310 jl=1,kproma
!re        IF(ldcum(jl)) THEN
!re           jk=kctop(jl)-1
!re           pmfuxt(jl,jk,jt)=pxtu(jl,jk,jt)*pmfu(jl,jk)
!re        ENDIF
!re5310 END DO
!re5312 END DO
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
END SUBROUTINE h2oiso_cuasc


!======================================================================


SUBROUTINE h2oiso_cuasct(   kproma, kbdim, klev, klevp1, klevm1,              &
           ptenh,    pqenh,    puen,     pven,                         &
!re           ktrac,                                                      &
!re           pxtenh,   pxten,    pxtu,     pmfuxt,                       &
           pten,     pqen,     pqsen,                                  &
           pgeo,     pgeoh,    paphp1,                                 &
           pqte,     pverv,    klwmin,                                 &
           ldcum,    ldland,   ktype,    klab,                         &
           ptu,      pqu,      plu,      puu,      pvu,                &
           pmfu,     pmfub,    pentr,                                  &
           pmfus,    pmfuq,                                            &
           pmful,    plude,    pqude,    pdmfup,                       &
           pcpen,    pcpcu,                                            &
           kcbot,    kctop,    kctop0,                                 &
!---wiso-code
           kwiso,                                                      &
           pwisoqenh,                                                  &
           pwisoqen,                                                   &
           pwisoqu,  pwisolu,                                          &
           pwisomfuq,                                                  &
           pwisomful,pwisolude,pwisoqude,pwisodmfup                    &
           , nn, time_step_len)
!---wiso-code-end
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
!          THIS ROUTINE IS CALLED FROM *CUMASTRT*.
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
!
!          REFERENCE
!          ---------
!          (TIEDTKE,1989)
!
!reUSE mo_kind,         ONLY : dp
!reUSE mo_control,      ONLY : nn
!reUSE mo_constants,    ONLY : g, tmelt, vtmpc1
!reUSE mo_cumulus_flux, ONLY : lmfdudv, lmfmid, nmctop, cmfcmin, cprcon   &
!re                          , cmfctop
!reUSE mo_time_control, ONLY : time_step_len

!---wiso-code

  USE messy_h2oiso,     ONLY : fractcal
  USE messy_main_tools_wiso, ONLY: wiso_frac_liq_ice, tnat, cthomi, cwisomin, cwisosec

!---wiso-code-end

IMPLICIT NONE

INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, klevm1!, ktrac

REAL(DP),INTENT(IN):: time_step_len
INTEGER, INTENT (IN) :: nn
!---wiso-code

INTEGER, INTENT (IN) :: kwiso

!---wiso-code-end

INTEGER :: jl, jk, jt, ik, icall 
!
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
REAL(dp) :: pcpen(kbdim,klev),       pcpcu(kbdim,klev)
INTEGER  :: klwmin(kbdim),           ktype(kbdim),                     &
            klab(kbdim,klev),        kcbot(kbdim),                     &
            kctop(kbdim),            kctop0(kbdim)
LOGICAL  :: ldcum(kbdim),            ldland(kbdim)
!
REAL(dp) :: zdmfen(kbdim),           zdmfde(kbdim),                    &
            zmfuu(kbdim),            zmfuv(kbdim),                     &
            zpbase(kbdim),           zqold(kbdim)
REAL(dp) :: zph(kbdim)
LOGICAL  :: loflag(kbdim)
!re REAL(dp) :: pxtenh(kbdim,klev,ktrac),pxten(kbdim,klev,ktrac),          &
!re            pxtu(kbdim,klev,ktrac),  pmfuxt(kbdim,klev,ktrac)
!
!---wiso-code

REAL(dp)::  pwisoqenh(kbdim,klev,kwiso),                               &
            pwisoqen(kbdim,klev,kwiso)

REAL(dp) :: pwisoqu(kbdim,klev,kwiso),                                 &
            pwisomfuq(kbdim,klev,kwiso),                               &
            pwisolu(kbdim,klev,kwiso),pwisolude(kbdim,klev,kwiso),     &
            pwisoqude(kbdim,klev,kwiso),                               &
            pwisomful(kbdim,klev,kwiso),pwisodmfup(kbdim,klev,kwiso)

!---wiso-code-end

REAL(dp) :: zcons2, ztglace, zmfmax, zfac, zmftest, zqeen, zseen       &
          , zscde, zqude, zmfusk, zmfuqk, zmfulk, zxteen, zxtude       &
          , zmfuxtk, zbuo, zdnoprc, zprcon, zlnew, zz, zdmfeu, zdmfdu  &
          , zzdmf, zdlev
!
!---wiso-code

REAL(dp) :: zprec_tmp(kbdim,klev)

REAL(dp) :: zwisofracliq(kbdim,kwiso), zwisofrac(kbdim,kwiso), zwisofracice(kbdim,kwiso)

REAL(dp) :: zwisoqeen,                                                 &
            zwisoqude, zwisomfuqk, zwisomfulk

REAL(dp) :: zqliq, zqice,                                              &
            zql,  zqv,   zwisoql, zwisoqv,                             &
            zquot,zquot2,zqvo,    zdelta

LOGICAL  :: lo
            
!---wiso-code-end

!      INTRINSIC FUNCTIONS
INTRINSIC MAX, MIN
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

!---wiso-code

     DO jt=1,kwiso
       DO jl=1,kproma
         pwisolu(jl,jk,jt)=0._dp
         pwisomfuq(jl,jk,jt)=0._dp
         pwisomful(jl,jk,jt)=0._dp
         pwisolude(jl,jk,jt)=0._dp
         pwisoqude(jl,jk,jt)=0._dp
         pwisodmfup(jl,jk,jt)=0._dp
       END DO
     END DO

!---wiso-code-end

!re     DO 2204 jt=1,ktrac
!re        DO 2202 jl=1,kproma
!re           pmfuxt(jl,jk,jt)=0._dp
!re2202    END DO
!re2204 END DO
!
230 END DO
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
!---wiso-code

  DO jt=1,kwiso
    DO jl=1,kproma
      IF(.NOT.ldcum(jl)) THEN
        pwisoqu(jl,klev,jt)=0._dp
      ENDIF
     pwisomfuq(jl,klev,jt)=pmfub(jl)*pwisoqu(jl,klev,jt)
    END DO
  END DO

!---wiso-code-end

!re  DO 3112 jt=1,ktrac
!re     DO 3110 jl=1,kproma
!re        IF(.NOT.ldcum(jl)) THEN
!re           pxtu(jl,klev,jt)=0._dp
!re        ENDIF
!re        pmfuxt(jl,klev,jt)=pmfub(jl)*pxtu(jl,klev,jt)
!re3110 END DO
!re3112 END DO
!
  DO 320 jl=1,kproma
     ldcum(jl)=.FALSE.
320 END DO
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
        CALL h2oiso_cubasmc(kproma,   kbdim,    klev,     ik,       klab,     &
                     pten,     pqen,     pqsen,    puen,     pven,     &
!re                     ktrac,                                            &
!re                     pxten,    pxtu,     pmfuxt,                       &
                     pverv,    pgeo,     pgeoh,    ldcum,    ktype,    &
                     pmfu,     pmfub,    pentr,    kcbot,              &
                     ptu,      pqu,      plu,      puu,      pvu,      &
                     pmfus,    pmfuq,    pmful,    pdmfup,   zmfuu,    &
                     pcpen,                                            &
                     zmfuv,                                            &
!---wiso-code
                     kwiso,                                            &
                     pwisoqen,                                         &
                     pwisoqu,  pwisolu,                                &
                     pwisomfuq,pwisomful,pwisodmfup                    &
                     )
!---wiso-code-end
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

!---wiso-code

     DO jt=1,kwiso
       DO jl=1,kproma
         IF(ktype(jl).EQ.3.AND.jk.EQ.kcbot(jl)) THEN
           zmfmax=(paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2
           IF(pmfub(jl).GT.zmfmax) THEN
             zfac=zmfmax/pmfub(jl)
             pwisomfuq(jl,jk+1,jt)=pwisomfuq(jl,jk+1,jt)*zfac
           END IF
         END IF
       END DO
     END DO

!---wiso-code-end

!re     DO 4102 jt=1,ktrac
!re        DO 4101 jl=1,kproma
!re           IF(ktype(jl).EQ.3.AND.jk.EQ.kcbot(jl)) THEN
!re              zmfmax=(paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2
!re              IF(pmfub(jl).GT.zmfmax) THEN
!re                 zfac=zmfmax/pmfub(jl)
!re                 pmfuxt(jl,jk+1,jt)=pmfuxt(jl,jk+1,jt)*zfac
!re              END IF
!re           END IF
!re4101    END DO
!re4102 END DO
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
!*                 SPECIFY ENTRAINMENT RATES IN *CUENTRT*
!                  --------------------------------------
!
     ik=jk
     CALL h2oiso_cuentrt(kproma,  kbdim,    klev,     klevp1,   ik,           &
                 ptenh,    pqenh,    pqte,     paphp1,                 &
                 klwmin,   ldcum,    ktype,    kcbot,    kctop0,       &
                 zpbase,   pmfu,     pentr,                            &
                 zdmfen,   zdmfde)
!
!
!
!                  DO ADIABATIC ASCENT FOR ENTRAINING/DETRAINING PLUME
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
           zqeen=pqenh(jl,jk+1)*zdmfen(jl)
           zseen=(pcpcu(jl,jk+1)*ptenh(jl,jk+1)+pgeoh(jl,jk+1))        &
                                               *zdmfen(jl)
           zscde=(pcpcu(jl,jk+1)*ptu(jl,jk+1)+pgeoh(jl,jk+1))          &
                                               *zdmfde(jl)
           zqude=pqu(jl,jk+1)*zdmfde(jl)
           pqude(jl,jk)=zqude
           plude(jl,jk)=plu(jl,jk+1)*zdmfde(jl)
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
!---wiso-code

     DO jt=1,kwiso
       DO jl=1,kproma
         IF(loflag(jl)) THEN
           zwisoqeen=pwisoqenh(jl,jk+1,jt)*zdmfen(jl)
           zwisoqude=pwisoqu(jl,jk+1,jt)*zdmfde(jl)
           pwisoqude(jl,jk,jt)=zwisoqude
           pwisolude(jl,jk,jt)=pwisolu(jl,jk+1,jt)*zdmfde(jl)
           zwisomfuqk=pwisomfuq(jl,jk+1,jt)+zwisoqeen-zwisoqude
           zwisomfulk=pwisomful(jl,jk+1,jt)-pwisolude(jl,jk,jt)
           pwisolu(jl,jk,jt)=zwisomfulk*(1._dp/MAX(cmfcmin,pmfu(jl,jk)))
           pwisoqu(jl,jk,jt)=zwisomfuqk*(1._dp/MAX(cmfcmin,pmfu(jl,jk)))
         ENDIF
       END DO
     END DO

!---wiso-code-end

!re     DO 4204 jt=1,ktrac
!re        DO 4202 jl=1,kproma
!re           IF(loflag(jl)) THEN
!re              zxteen=pxtenh(jl,jk+1,jt)*zdmfen(jl)
!re              zxtude=pxtu(jl,jk+1,jt)*zdmfde(jl)
!re              zmfuxtk=pmfuxt(jl,jk+1,jt)+zxteen-zxtude
!re              pxtu(jl,jk,jt)=zmfuxtk*(1._dp/MAX(cmfcmin,pmfu(jl,jk)))
!re           ENDIF
!re4202    END DO
!re4204 END DO
!
!
!                  DO CORRECTIONS FOR MOIST ASCENT
!                  BY ADJUSTING T,Q AND L IN *CUADJTQ*
!                  -----------------------------------
!
     ik=jk
     icall=1
     CALL h2oiso_cuadjtq(kproma,   kbdim,    klev,     ik,                    &
                  zph,      ptu,      pqu,      loflag,   icall)
!

!---wiso-code

! calculate fractionation factors for liquid-vapour and solid-vapour phase change

!reold     CALL fractcal(kproma,kbdim,kwiso,zwisofrac,zwisofracice,ptu(:,ik),tmelt)
     CALL wiso_frac_liq_ice(kproma,kbdim,kwiso,ptu(:,ik),zwisofracliq,zwisofracice)

!    adjusting isotope vapour *pwisoqu* depending on fract.
     DO jt=1,kwiso
       DO jl=1,kproma
         IF (loflag(jl)) THEN
! divide into ice and liquid
           IF (ptu(jl,jk).gt.tmelt) THEN
             zqliq=1._dp
             zqice=0._dp
           ELSEIF (ptu(jl,jk).lt.cthomi) THEN
             zqice=1._dp
             zqliq=0._dp
           ELSE
             zqice=(tmelt-ptu(jl,jk))/(tmelt-cthomi)
             zqliq=1._dp-zqice
           ENDIF
! fractionation of ice and liquid
           zql=zqliq*(plu(jl,jk)+zqold(jl)-pqu(jl,jk))
           zqv=zqold(jl)+zqliq*(pqu(jl,jk)-zqold(jl))
           zwisoql=zqliq*pwisolu(jl,jk,jt)
           zwisoqv=pwisoqu(jl,jk,jt)
!reold           zquot=zqv+zwisofrac(jl,jt)*zql
           zquot=zqv+zwisofracliq(jl,jt)*zql
           zqvo=zqv
           IF (zquot.gt.cwisomin.and.zqvo.gt.cwisomin) THEN
!reold             zdelta=(zwisofrac(jl,jt)*zql*zwisoqv-zwisoql*zqv)/zquot 
             zdelta=(zwisofracliq(jl,jt)*zql*zwisoqv-zwisoql*zqv)/zquot 
             zqv=pqu(jl,jk)
             zwisoqv=zwisoqv-zdelta
             zquot2=zqv/zqvo
             lo=abs(1._dp-zquot2).lt.cwisosec
             zquot2=MERGE(1._dp,zquot2,lo)

if(zquot2.gt.0.0_dp)then!rehc
             zdelta=zdelta+zwisoqv*(1._dp-zquot2**zwisofracice(jl,jt))
endif!rehc
             pwisoqu(jl,jk,jt)=pwisoqu(jl,jk,jt)-zdelta
             IF(loflag(jl).AND.pqu(jl,jk).LT.zqold(jl)) THEN  ! add condensate to liquid as it is done
               pwisolu(jl,jk,jt)=pwisolu(jl,jk,jt)+zdelta     ! for normal water in next loop 440
             END IF
           ENDIF
         ENDIF
       END DO
     END DO

!---wiso-code-end

!DIR$ IVDEP
!OCL NOVREC
     DO 440 jl=1,kproma
        IF(loflag(jl).AND.pqu(jl,jk).LT.zqold(jl)) THEN
           klab(jl,jk)=2
           plu(jl,jk)=plu(jl,jk)+zqold(jl)-pqu(jl,jk)
           zbuo=ptu(jl,jk)*(1._dp+vtmpc1*pqu(jl,jk)-plu(jl,jk))-       &
                    ptenh(jl,jk)*(1._dp+vtmpc1*pqenh(jl,jk))
           IF(klab(jl,jk+1).EQ.1) zbuo=zbuo+0.5_dp
           IF(zbuo.GT.0._dp.AND.pmfu(jl,jk).GE.0.1_dp*pmfub(jl)) THEN
              kctop(jl)=jk
              ldcum(jl)=.TRUE.
              zdnoprc=MERGE(zdlev,1.5e4_dp,ldland(jl))
              zprcon=MERGE(0._dp,cprcon,                               &
                               zpbase(jl)-paphp1(jl,jk).LT.zdnoprc)
              zlnew=plu(jl,jk)/                                        &
                           (1._dp+zprcon*(pgeoh(jl,jk)-pgeoh(jl,jk+1)))
              pdmfup(jl,jk)=MAX(0._dp,(plu(jl,jk)-zlnew)*pmfu(jl,jk))
!---wiso-code
! Store fraction of precipitation formed from cloud water
              IF(plu(jl,jk).gt.cwisomin) THEN
                zprec_tmp(jl,jk)=(plu(jl,jk)-zlnew)/plu(jl,jk)
              ENDIF
              IF (abs(1._dp-zprec_tmp(jl,jk)).lt.cwisosec) zprec_tmp(jl,jk)=1._dp
!---wiso-code-end
              plu(jl,jk)=zlnew
           ELSE
              klab(jl,jk)=0
              pmfu(jl,jk)=0._dp
!---wiso-code
! No precipitation formed from cloud water
              zprec_tmp(jl,jk)=0._dp
!---wiso-code-end
           END IF
        END IF
440  END DO

!---wiso-code

     DO jt=1,kwiso
       DO jl=1,kproma
        IF(loflag(jl).AND.pqu(jl,jk).LT.zqold(jl)) THEN
! Calculate water isotope precipitation and new mass fluxes for the water isotopes in updrafts
           pwisodmfup(jl,jk,jt)=MAX(0._dp,pwisolu(jl,jk,jt)*pmfu(jl,jk)*zprec_tmp(jl,jk))
           pwisolu(jl,jk,jt)=(1._dp-zprec_tmp(jl,jk))*pwisolu(jl,jk,jt)
         END IF
       END DO
     END DO

!---wiso-code-end

     DO 455 jl=1,kproma
        IF(loflag(jl)) THEN
           pmful(jl,jk)=plu(jl,jk)*pmfu(jl,jk)
           pmfus(jl,jk)=(pcpcu(jl,jk)*ptu(jl,jk)+pgeoh(jl,jk))         &
                                 *pmfu(jl,jk)
           pmfuq(jl,jk)=pqu(jl,jk)*pmfu(jl,jk)
        END IF
455  END DO

!---wiso-code

     DO jt=1,kwiso
       DO jl=1,kproma
         IF(loflag(jl)) THEN
           pwisomful(jl,jk,jt)=pwisolu(jl,jk,jt)*pmfu(jl,jk)
           pwisomfuq(jl,jk,jt)=pwisoqu(jl,jk,jt)*pmfu(jl,jk)
         END IF
       END DO
     END DO

!---wiso-code-end

!re     DO 4554 jt=1,ktrac
!re        DO 4552 jl=1,kproma
!re           IF(loflag(jl)) THEN
!re              pmfuxt(jl,jk,jt)=pxtu(jl,jk,jt)*pmfu(jl,jk)
!re           ENDIF
!re4552    END DO
!re4554 END DO
!
     IF(lmfdudv) THEN
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
        plude(jl,jk-1)=pmful(jl,jk)
        pqude(jl,jk-1)=pmfuq(jl,jk)
     END IF
530 END DO

!---wiso-code

! Calculate water isotope fluxes above the non-buoyancy level
     DO jt=1,kwiso
!DIR$ IVDEP
       DO jl=1,kproma
         IF (ldcum(jl)) THEN
           jk=kctop(jl)-1
           pwisolude(jl,jk,jt)=zdmfde(jl)*pwisolu(jl,jk+1,jt)
           pwisoqude(jl,jk,jt)=zdmfde(jl)*pwisoqu(jl,jk+1,jt)
           pwisodmfup(jl,jk,jt)=0._dp
           pwisomfuq(jl,jk,jt)=pwisoqu(jl,jk,jt)*pmfu(jl,jk)
           pwisomful(jl,jk,jt)=pwisolu(jl,jk,jt)*pmfu(jl,jk)
           pwisolude(jl,jk-1,jt)=pwisomful(jl,jk,jt)
           pwisoqude(jl,jk-1,jt)=pwisomfuq(jl,jk,jt)
         END IF
       END DO
     END DO

!---wiso-code-end

!re  DO 5312 jt=1,ktrac
!re     DO 5310 jl=1,kproma
!re        IF(ldcum(jl)) THEN
!re           jk=kctop(jl)-1
!re           pmfuxt(jl,jk,jt)=pxtu(jl,jk,jt)*pmfu(jl,jk)
!re        ENDIF
!re5310 END DO
!re5312 END DO
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
END SUBROUTINE h2oiso_cuasct



!======================================================================



SUBROUTINE h2oiso_cudlfs(   kproma, kbdim, klev, klevp1,                      &
           ptenh,    pqenh,    puen,     pven,                         &
!re           ktrac,                                                      &
!re           pxtenh,   pxtu,     pxtd,     pmfdxt,                       &
           pgeoh,    paphp1,                                           &
           ptu,      pqu,      puu,      pvu,                          &
           ldcum,    kcbot,    kctop,    pmfub,    prfl,               &
           ptd,      pqd,      pud,      pvd,                          &
           pmfd,     pmfds,    pmfdq,    pdmfdp,                       &
           pcpcu,                                                      &
           kdtop,    lddraf,                                           &
!---wiso-code
           kwiso,                                                      &
           pwisoqenh,                                                  &
           pwisoqu,                                                    &
           pwisoqd,                                                    &
           pwisomfdq,pwisodmfdp,                                       &
           pdmfup,   pwisodmfup                                        &
           )
!---wiso-code-end

!          THIS ROUTINE CALCULATES LEVEL OF FREE SINKING FOR
!          CUMULUS DOWNDRAFTS AND SPECIFIES T,Q,U AND V VALUES
!
!          M.TIEDTKE         E.C.M.W.F.    12/86 MODIF. 12/89
!
!          PURPOSE.
!          --------
!          TO PRODUCE LFS-VALUES FOR CUMULUS DOWNDRAFTS
!          FOR MASSFLUX CUMULUS PARAMETERIZATION
!
!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!          INPUT ARE ENVIRONMENTAL VALUES OF T,Q,U,V,P,PHI
!          AND UPDRAFT VALUES T,Q,U AND V AND ALSO
!          CLOUD BASE MASSFLUX AND CU-PRECIPITATION RATE.
!          IT RETURNS T,Q,U AND V VALUES AND MASSFLUX AT LFS.
!
!          METHOD.
!          --------
!          CHECK FOR NEGATIVE BUOYANCY OF AIR OF EQUAL PARTS OF
!          MOIST ENVIRONMENTAL AIR AND CLOUD AIR.
!
!          EXTERNALS
!          ---------
!          *CUADJTQ* FOR CALCULATING WET BULB T AND Q AT LFS
!
!!$ USE mo_kind,          ONLY : dp
!!$ USE mo_constants,     ONLY : vtmpc1
!!$ USE mo_cumulus_flux,  ONLY : lmfdudv, lmfdd, cmfdeps
!
!---wiso-code

!!$USE mo_constants,     ONLY : tmelt
!reoldUSE messy_h2oiso,          ONLY : talph1, talph2, talph3, cwisomin
USE messy_main_tools_wiso,     ONLY : talphal1, talphal2, talphal3, cwisomin
!---wiso-code-end

IMPLICIT NONE
!
INTEGER, INTENT (IN) :: kbdim, klev, kproma, klevp1!re, ktrac

!---wiso-code

INTEGER, INTENT (IN) :: kwiso

!---wiso-code-end

REAL(dp) :: ptenh(kbdim,klev),       pqenh(kbdim,klev),                &
            puen(kbdim,klev),        pven(kbdim,klev),                 &
            pgeoh(kbdim,klev),       paphp1(kbdim,klevp1),             &
            ptu(kbdim,klev),         pqu(kbdim,klev),                  &
            puu(kbdim,klev),         pvu(kbdim,klev),                  &
            pmfub(kbdim),            prfl(kbdim)
!
REAL(dp) :: ptd(kbdim,klev),         pqd(kbdim,klev),                  &
            pud(kbdim,klev),         pvd(kbdim,klev),                  &
            pmfd(kbdim,klev),        pmfds(kbdim,klev),                &
            pmfdq(kbdim,klev),       pdmfdp(kbdim,klev)
REAL(dp) :: pcpcu(kbdim,klev)
INTEGER  :: kcbot(kbdim),            kctop(kbdim),                     &
            kdtop(kbdim)
LOGICAL  :: ldcum(kbdim),            lddraf(kbdim)
!
!---wiso-code

REAL(dp) :: pwisoqenh(kbdim,klev,kwiso),                               &
            pwisoqu(kbdim,klev,kwiso)

REAL(dp) :: pwisoqd(kbdim,klev,kwiso),                                 &
            pwisomfdq(kbdim,klev,kwiso),pwisodmfdp(kbdim,klev,kwiso)

REAL(dp) :: pdmfup(kbdim,klev),      pwisodmfup(kbdim,klev,kwiso)

!---wiso-code-end

REAL(dp) :: ztenwb(kbdim,klev),      zqenwb(kbdim,klev),               &
            zcond(kbdim)
REAL(dp) :: zph(kbdim)
LOGICAL  :: llo2(kbdim)
!re REAL(dp) :: pxtenh(kbdim,klev,ktrac),pxtu(kbdim,klev,ktrac),           &
!re            pxtd(kbdim,klev,ktrac),  pmfdxt(kbdim,klev,ktrac)
LOGICAL  :: llo3(kbdim)
INTEGER  :: jl, jk, ke, is, ik, icall, jt
REAL(DP) :: zttest, zqtest, zbuo, zmftop
!
!---wiso-code

REAL(dp) :: zrain(kbdim),            zwisorain(kbdim,kwiso)

REAL(DP) :: zqv, zql, zwisoqv, zwisoql,                                &
            zt, zwisofracliq, zquot, zdelta!reold, zwisofrac

!---wiso-code-end

!
!----------------------------------------------------------------------
!
!     1.           SET DEFAULT VALUES FOR DOWNDRAFTS
!                  ---------------------------------
!
100 CONTINUE
  DO 110 jl=1,kproma
     lddraf(jl)=.FALSE.
     kdtop(jl)=klevp1
110 END DO
!
!---wiso-code

  zrain(:) = 0._dp
  zwisorain(:,:) = 0._dp

!---wiso-code-end

  IF(.NOT.lmfdd) go to 300
!
!
!----------------------------------------------------------------------
!
!     2.           DETERMINE LEVEL OF FREE SINKING BY
!                  DOING A SCAN FROM TOP TO BASE OF CUMULUS CLOUDS
!
!                  FOR EVERY POINT AND PROCEED AS FOLLOWS:
!
!                    (1) DETEMINE WET BULB ENVIRONMENTAL T AND Q
!                    (2) DO MIXING WITH CUMULUS CLOUD AIR
!                    (3) CHECK FOR NEGATIVE BUOYANCY
!
!                  THE ASSUMPTION IS THAT AIR OF DOWNDRAFTS IS MIXTURE
!                  OF 50% CLOUD AIR + 50% ENVIRONMENTAL AIR AT WET BULB
!                  TEMPERATURE (I.E. WHICH BECAME SATURATED DUE TO
!                  EVAPORATION OF RAIN AND CLOUD WATER)
!                  ----------------------------------------------------
!
200 CONTINUE
!
  ke=klev-3
  DO 290 jk=3,ke
!
!
!     2.1          CALCULATE WET-BULB TEMPERATURE AND MOISTURE
!                  FOR ENVIRONMENTAL AIR IN *CUADJTQ*
!                  -------------------------------------------
!
!---wiso-code

! Calculate Precipitation on this level
     DO jl=1,kproma
       zrain(jl)=zrain(jl)+pdmfup(jl,jk)
     END DO
     DO jt=1,kwiso
       DO jl=1,kproma
         zwisorain(jl,jt)=zwisorain(jl,jt)+pwisodmfup(jl,jk,jt)
       END DO
     END DO
     
!---wiso-code-end

210  CONTINUE
     is=0
     DO 212 jl=1,kproma
        ztenwb(jl,jk)=ptenh(jl,jk)
        zqenwb(jl,jk)=pqenh(jl,jk)
        zph(jl)=paphp1(jl,jk)
        llo2(jl)=ldcum(jl).AND.prfl(jl).GT.0._dp.AND..NOT.lddraf(jl)   &
                          .AND.(jk.LT.kcbot(jl).AND.jk.GT.kctop(jl))
        is=is+MERGE(1,0,llo2(jl))
212  END DO
!re     IF(is.EQ.0) go to 290 ! mz_ht_20071217
!
     ik=jk
     icall=2
     CALL h2oiso_cuadjtq(kproma, kbdim, klev, ik,                             &
                  zph,      ztenwb,   zqenwb,   llo2,     icall)
!
!
!     2.2          DO MIXING OF CUMULUS AND ENVIRONMENTAL AIR
!                  AND CHECK FOR NEGATIVE BUOYANCY.
!                  THEN SET VALUES FOR DOWNDRAFT AT LFS.
!                  ----------------------------------------
!
220  CONTINUE
!DIR$ IVDEP
!OCL NOVREC
     DO 222 jl=1,kproma
        llo3(jl)=.FALSE.
        IF(llo2(jl)) THEN
           zttest=0.5_dp*(ptu(jl,jk)+ztenwb(jl,jk))
           zqtest=0.5_dp*(pqu(jl,jk)+zqenwb(jl,jk))
           zbuo=zttest*(1._dp+vtmpc1*zqtest)-                          &
                         ptenh(jl,jk)*(1._dp+vtmpc1*pqenh(jl,jk))
           zcond(jl)=pqenh(jl,jk)-zqenwb(jl,jk)
           zmftop=-cmfdeps*pmfub(jl)
           IF(zbuo.LT.0._dp.AND.prfl(jl).GT.10._dp*zmftop*zcond(jl))   &
                                                                  THEN
              llo3(jl)=.TRUE.
              kdtop(jl)=jk
              lddraf(jl)=.TRUE.
              ptd(jl,jk)=zttest
              pqd(jl,jk)=zqtest
              pmfd(jl,jk)=zmftop
              pmfds(jl,jk)=pmfd(jl,jk)*(pcpcu(jl,jk)*ptd(jl,jk)        &
                                               +pgeoh(jl,jk))
              pmfdq(jl,jk)=pmfd(jl,jk)*pqd(jl,jk)
              pdmfdp(jl,jk-1)=-0.5_dp*pmfd(jl,jk)*zcond(jl)
              prfl(jl)=prfl(jl)+pdmfdp(jl,jk-1)
           END IF
        END IF
222  END DO
!
!---wiso-code

     DO jl=1,kproma
       IF(llo3(jl)) THEN
         zrain(jl)=zrain(jl)+pdmfdp(jl,jk-1)
       ENDIF
     END DO

     DO jt=1,kwiso
       DO jl=1,kproma
         IF(llo3(jl)) THEN
          zwisoql=zwisorain(jl,jt)
          zwisoqv=pwisoqu(jl,jk,jt)+pwisoqenh(jl,jk,jt)
          zql=zrain(jl)
          zqv=pqd(jl,jk)
          zt=ptd(jl,jk)

! fractionation over water
!reold          zwisofrac=exp(talph1(jt)/(zt**2)+talph2(jt)/zt+talph3(jt))
          zwisofracliq=exp(talphal1(jt)/(zt**2)+talphal2(jt)/zt+talphal3(jt))

!reold          zquot=0.5_dp*(pmfd(jl,jk)*zqv-zwisofrac*zql)
          zquot=0.5_dp*(pmfd(jl,jk)*zqv-zwisofracliq*zql)
          IF (abs(zquot).lt.cwisomin) GOTO 2220

!reold          zdelta=(zwisoql*zqv-0.5_dp*zwisofrac*zql*zwisoqv)/zquot
          zdelta=(zwisoql*zqv-0.5_dp*zwisofracliq*zql*zwisoqv)/zquot
          pwisoqd(jl,jk,jt)=0.5_dp*(zwisoqv-zdelta)
          pwisomfdq(jl,jk,jt)=pmfd(jl,jk)*pwisoqd(jl,jk,jt)
          pwisodmfdp(jl,jk-1,jt)=-0.5_dp*pmfd(jl,jk)*zdelta

        ENDIF
 2220   CONTINUE
 
      END DO
    END DO
!---wiso-code-end

!re     DO 2224 jt=1,ktrac
!re        DO 2222 jl=1,kproma
!re           IF(llo3(jl)) THEN
!re              pxtd(jl,jk,jt)=0.5_dp*(pxtu(jl,jk,jt)+pxtenh(jl,jk,jt))
!re              pmfdxt(jl,jk,jt)=pmfd(jl,jk)*pxtd(jl,jk,jt)
!re           ENDIF
!re2222    END DO
!re2224 END DO
!
     IF(lmfdudv) THEN
        DO 224 jl=1,kproma
           IF(pmfd(jl,jk).LT.0._dp) THEN
              pud(jl,jk)=0.5_dp*(puu(jl,jk)+puen(jl,jk-1))
              pvd(jl,jk)=0.5_dp*(pvu(jl,jk)+pven(jl,jk-1))
           END IF
224     END DO
     END IF
!
290 END DO
!
300 CONTINUE
  RETURN
END SUBROUTINE h2oiso_cudlfs




!======================================================================





SUBROUTINE h2oiso_cuddraf(  kproma, kbdim, klev, klevp1,                      &
           ptenh,    pqenh,    puen,     pven,                         &
!re           ktrac,                                                      &
!re           pxtenh,   pxtd,     pmfdxt,                                 &
           pgeoh,    paphp1,   prfl,                                   &
           ptd,      pqd,      pud,      pvd,                          &
           pmfd,     pmfds,    pmfdq,    pdmfdp,                       &
           pcpcu,                                                      &
           lddraf,                                                     &
!---wiso-code
           kwiso,                                                      &
           pwisoqenh,                                                  &
           pwisoqd,                                                    &
           pwisomfdq,pwisodmfdp,                                       &
           pdmfup,   pwisodmfup                                        &
           )
!---wiso-code-end
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
!!$USE mo_kind,         ONLY : dp
!!$USE mo_constants,    ONLY : g, rd, vtmpc1
!!$USE mo_cumulus_flux, ONLY : lmfdudv, cmfcmin, entrdd
!
!---wiso-code

!!$USE mo_constants,     ONLY : tmelt
!reoldUSE messy_h2oiso,          ONLY : talph1, talph2, talph3, cwisomin
USE messy_main_tools_wiso,   ONLY : talphal1, talphal2, talphal3, cwisomin

!---wiso-code-end

IMPLICIT NONE
!
INTEGER, INTENT (IN) :: kbdim, klev, kproma, klevp1!re, ktrac
!

!---wiso-code

INTEGER, INTENT (IN) :: kwiso

!---wiso-code-end

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
!---wiso-code

REAL(dp) :: pwisoqenh(kbdim,klev,kwiso)                                 

REAL(dp) :: pwisoqd(kbdim,klev,kwiso),                                 &
            pwisomfdq(kbdim,klev,kwiso),pwisodmfdp(kbdim,klev,kwiso)

REAL(dp) :: pdmfup(kbdim,klev),      pwisodmfup(kbdim,klev,kwiso)

!---wiso-code-end

REAL(dp) :: zdmfen(kbdim),           zdmfde(kbdim),                    &
            zcond(kbdim)
REAL(dp) :: zph(kbdim)
LOGICAL  :: llo2(kbdim)
!re REAL(dp) :: pxtenh(kbdim,klev,ktrac),pxtd(kbdim,klev,ktrac),           &
!re            pmfdxt(kbdim,klev,ktrac)
LOGICAL  :: llo1
INTEGER  :: jk, is, jl, itopde, jt, ik, icall
REAL(dp) :: zentr, zseen, zqeen, zsdde, zqdde, zmfdsk, zmfdqk, zxteen  &
          , zxtdde, zmfdxtk, zbuo, zdmfdp, zmfduk, zmfdvk
!
!---wiso-code

REAL(dp) :: zrain(kbdim),            zwisorain(kbdim,kwiso)

REAL(dp) :: zwisoqeen,                                                 &
            zwisoqdde, zwisomfdqk
            
REAL(dp) :: zqv, zql, zwisoqv, zwisoql,                                &
            zt, zwisofracliq, zquot, zdelta!reold, zwisofrac

!---wiso-code-end

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
!---wiso-code

  zrain(:) = 0._dp
  zwisorain(:,:) = 0._dp

!---wiso-code-end

100 CONTINUE
  DO 180 jk=3,klev

!---wiso-code

     DO jl=1,kproma
       zrain(jl)=zrain(jl)+pdmfup(jl,jk)
     END DO
     DO jt=1,kwiso
       DO jl=1,kproma
         zwisorain(jl,jt)=zwisorain(jl,jt)+pwisodmfup(jl,jk,jt)
       END DO
     END DO

!---wiso-code-end

     is=0
     DO 110 jl=1,kproma
        zph(jl)=paphp1(jl,jk)
        llo2(jl)=lddraf(jl).AND.pmfd(jl,jk-1).LT.0._dp
        is=is+MERGE(1,0,llo2(jl))
110  END DO
!re     IF(is.EQ.0) go to 180 ! mz_ht_20071217
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
!---wiso-code

     DO jt=1,kwiso
       DO jl=1,kproma
         IF (llo2(jl)) THEN
           zwisoqeen = pwisoqenh(jl,jk-1,jt)*zdmfen(jl)
           zwisoqdde = pwisoqd(jl,jk-1,jt)*zdmfde(jl)
           zwisomfdqk = pwisomfdq(jl,jk-1,jt)+zwisoqeen-zwisoqdde
           pwisoqd(jl,jk,jt) = zwisomfdqk*(1._dp/MIN(-cmfcmin,pmfd(jl,jk)))
         END IF
       END DO
     END DO

!---wiso-code-end

!re     DO 1264 jt=1,ktrac
!re        DO 1262 jl=1,kproma
!re           IF(llo2(jl)) THEN
!re              zxteen=pxtenh(jl,jk-1,jt)*zdmfen(jl)
!re              zxtdde=pxtd(jl,jk-1,jt)*zdmfde(jl)
!re              zmfdxtk=pmfdxt(jl,jk-1,jt)+zxteen-zxtdde
!re              pxtd(jl,jk,jt)=zmfdxtk*(1._dp/MIN(-cmfcmin,pmfd(jl,jk)))
!re           ENDIF
!re1262    END DO
!re1264 END DO
!
!
     ik=jk
     icall=2
     CALL h2oiso_cuadjtq(kproma,   kbdim,    klev,     ik,                    &
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
!---wiso-code

     DO jt=1,kwiso
       DO jl=1,kproma
         IF(llo2(jl)) THEN
           zwisoqv=pwisoqd(jl,jk,jt)
           zwisoql=zwisorain(jl,jt)    
           zqv=pqd(jl,jk)
           zql=zrain(jl)+pdmfdp(jl,jk-1)
           zt=ptenh(jl,jk)
! fractionation over water
!reold           zwisofrac=exp(talph1(jt)/(zt**2)+talph2(jt)/zt+talph3(jt))
!reold           zquot=zwisofrac*zql-pmfd(jl,jk)*zqv
           zwisofracliq=exp(talphal1(jt)/(zt**2)+talphal2(jt)/zt+talphal3(jt))
           zquot=zwisofracliq*zql-pmfd(jl,jk)*zqv
           IF (abs(zquot).lt.cwisomin) GOTO 1500
!reold           zdelta=(zwisofrac*zql*zwisoqv-zwisoql*zqv)/zquot 
           zdelta=(zwisofracliq*zql*zwisoqv-zwisoql*zqv)/zquot 
           pwisoqd(jl,jk,jt)=pwisoqd(jl,jk,jt)-zdelta
           pwisomfdq(jl,jk,jt)=pwisoqd(jl,jk,jt)*pmfd(jl,jk)  
           pwisodmfdp(jl,jk-1,jt)=-pmfd(jl,jk)*zdelta
         ENDIF
1500     CONTINUE
       END DO
     END DO

!---wiso-code-end

!re     DO 1504 jt=1,ktrac
!re        DO 1502 jl=1,kproma
!re           IF(llo2(jl)) THEN
!re              pmfdxt(jl,jk,jt)=pxtd(jl,jk,jt)*pmfd(jl,jk)
!re           ENDIF
!re1502    END DO
!re1504 END DO
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
END SUBROUTINE h2oiso_cuddraf




!======================================================================




SUBROUTINE h2oiso_cuflx(    kproma, kbdim, klev, klevp1,                      &
           pqen,     pqsen,    ptenh,    pqenh,                        &
!re           ktrac,                                                      &
!re           pxtenh,   pmfuxt,   pmfdxt,                                 &
           paphp1,   pgeoh,                                            &
           kcbot,    kctop,    kdtop,                                  &
           ktype,    lddraf,   ldcum,                                  &
           pmfu,     pmfd,     pmfus,    pmfds,                        &
           pmfuq,    pmfdq,    pmful,                                  &
           pdmfup,   pdmfdp,   prfl,     prain,                        &
           pcpcu,                                                      &
           pten,     psfl,     pdpmel,   ktopm2,                       &
!---wiso-code
           kwiso,                                                      &
           pwisoqenh,                                                  &
           pwisomfuq,pwisomfdq,pwisomful,                              &
           pwisodmfdp,                                                 &
           pmelt_tmp,pprec_frac,prain_tmp,                             &
           prfl_tmp, psfl_tmp                                          &
           , time_step_len)
!---wiso-code-end
!
!          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
!
!          PURPOSE
!          -------
!
!          THIS ROUTINE DOES THE FINAL CALCULATION OF CONVECTIVE
!          FLUXES IN THE CLOUD LAYER AND IN THE SUBCLOUD LAYER
!
!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!
!          EXTERNALS
!          ---------
!          NONE
!
!!$USE mo_kind,         ONLY: dp
!!$USE mo_constants,    ONLY: g, alf, cpd, tmelt, vtmpc2
!!$USE mo_physc2,       ONLY: cevapcu
USE messy_convect_tiedtke,       ONLY: cevapcu
!!$USE mo_time_control, ONLY: time_step_len

!---wiso-code
USE messy_main_tools_wiso,         ONLY: cwisosec
!---wiso-code-end

!
IMPLICIT NONE
!
INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1!re, ktrac
REAL(DP), INTENT(IN) :: time_step_len
!---wiso-code

INTEGER, INTENT (IN) :: kwiso

!---wiso-code-end

INTEGER, INTENT (OUT):: ktopm2
!
REAL(dp):: pqen(kbdim,klev),        pqsen(kbdim,klev),                 &
           ptenh(kbdim,klev),       pqenh(kbdim,klev),                 &
           paphp1(kbdim,klevp1),    pgeoh(kbdim,klev)
!
REAL(dp):: pmfu(kbdim,klev),        pmfd(kbdim,klev),                  &
           pmfus(kbdim,klev),       pmfds(kbdim,klev),                 &
           pmfuq(kbdim,klev),       pmfdq(kbdim,klev),                 &
           pdmfup(kbdim,klev),      pdmfdp(kbdim,klev),                &
           pmful(kbdim,klev),                                          &
           prfl(kbdim),             prain(kbdim)
REAL(dp):: pten(kbdim,klev),        pdpmel(kbdim,klev),                &
           psfl(kbdim)
REAL(dp):: pcpcu(kbdim,klev)
INTEGER :: kcbot(kbdim),            kctop(kbdim),                      &
           kdtop(kbdim),            ktype(kbdim)
LOGICAL :: lddraf(kbdim),           ldcum(kbdim)
!reREAL(dp):: pxtenh(kbdim,klev,ktrac),                                   &
!re           pmfuxt(kbdim,klev,ktrac),pmfdxt(kbdim,klev,ktrac)
!
!---wiso-code

REAL(dp):: pwisoqenh(kbdim,klev,kwiso)

REAL(dp):: pwisomfuq(kbdim,klev,kwiso),pwisomfdq(kbdim,klev,kwiso),    &
           pwisodmfdp(kbdim,klev,kwiso),                               &
           pwisomful(kbdim,klev,kwiso)

REAL(dp):: pmelt_tmp(kbdim,klev),pprec_frac(kbdim,klev),               &
           prain_tmp(kbdim,klev),                                      &
           prfl_tmp(kbdim), psfl_tmp(kbdim)

!---wiso-code-end

INTEGER :: jl, jk, jt, ikb
REAL(dp):: zcons1, zcons2, zcucov, ztmelp2, zzp, zfac, zsnmlt          &
         , zrfl, zrnew, zrmin, zrfln, zdrfl, zrsum, zdpevap
REAL(dp):: zpsubcl(kbdim)
!
!---wiso-code

REAL(dp):: zrfl_frac(kbdim)

!---wiso-code-end

!*             SPECIFY CONSTANTS
!

  zcons1=cpd/(alf*g*time_step_len)
  zcons2=1._dp/(g*time_step_len)
  zcucov=0.05_dp
  ztmelp2=tmelt+2._dp
!
!
!
!*    1.0          DETERMINE FINAL CONVECTIVE FLUXES
!                  ---------------------------------
!
100 CONTINUE
!  itop=klev
  DO 110 jl=1,kproma
!     itop=MIN(itop,kctop(jl))
     IF(.NOT.ldcum(jl).OR.kdtop(jl).LT.kctop(jl)) lddraf(jl)=.FALSE.
     IF(.NOT.ldcum(jl)) ktype(jl)=0
110 END DO
  ktopm2=1
  DO 120 jk=ktopm2,klev
!DIR$ IVDEP
!OCL NOVREC
     DO 115 jl=1,kproma
        IF(ldcum(jl).AND.jk.GE.kctop(jl)-1) THEN
           pmfus(jl,jk)=pmfus(jl,jk)-pmfu(jl,jk)*                      &
                               (pcpcu(jl,jk)*ptenh(jl,jk)+pgeoh(jl,jk))
           pmfuq(jl,jk)=pmfuq(jl,jk)-pmfu(jl,jk)*pqenh(jl,jk)
           IF(lddraf(jl).AND.jk.GE.kdtop(jl)) THEN
              pmfds(jl,jk)=pmfds(jl,jk)-pmfd(jl,jk)*                   &
                               (pcpcu(jl,jk)*ptenh(jl,jk)+pgeoh(jl,jk))
              pmfdq(jl,jk)=pmfdq(jl,jk)-pmfd(jl,jk)*pqenh(jl,jk)
           ELSE
              pmfd(jl,jk)=0._dp
              pmfds(jl,jk)=0._dp
              pmfdq(jl,jk)=0._dp
              pdmfdp(jl,jk-1)=0._dp
           END IF
        END IF
115  END DO
!

!---wiso-code

! Calculate final fluxes of the water isotopes (updrafts and downdrafts)
     DO jt=1,kwiso
       DO jl=1,kproma
         IF(ldcum(jl).AND.jk.GE.kctop(jl)-1) THEN
           pwisomfuq(jl,jk,jt)=pwisomfuq(jl,jk,jt)-pmfu(jl,jk)*pwisoqenh(jl,jk,jt)
           IF(lddraf(jl).AND.jk.GE.kdtop(jl)) THEN
             pwisomfdq(jl,jk,jt)=pwisomfdq(jl,jk,jt)-pmfd(jl,jk)*pwisoqenh(jl,jk,jt)
           ELSE
             pwisomfdq(jl,jk,jt)=0._dp
             pwisodmfdp(jl,jk-1,jt)=0._dp
           END IF
         END IF
       END DO
     END DO

!---wiso-code-end

!re     DO 1154 jt=1,ktrac
!re        DO 1152 jl=1,kproma
!re           IF(ldcum(jl).AND.jk.GE.kctop(jl)-1) THEN
!re              pmfuxt(jl,jk,jt)=pmfuxt(jl,jk,jt)                        &
!re                                     -pmfu(jl,jk)*pxtenh(jl,jk,jt)
!re              IF(lddraf(jl).AND.jk.GE.kdtop(jl)) THEN
!re                 pmfdxt(jl,jk,jt)=pmfdxt(jl,jk,jt)                     &
!re                                     -pmfd(jl,jk)*pxtenh(jl,jk,jt)
!re              ELSE
!re                 pmfdxt(jl,jk,jt)=0._dp
!re              ENDIF
!re           ELSE
!re              pmfuxt(jl,jk,jt)=0._dp
!re              pmfdxt(jl,jk,jt)=0._dp
!re           ENDIF
!re1152    END DO
!re1154 END DO
!
120 END DO
  DO 130 jk=ktopm2,klev
!DIR$ IVDEP
!OCL NOVREC
     DO 125 jl=1,kproma
        IF(ldcum(jl).AND.jk.GT.kcbot(jl)) THEN
           ikb=kcbot(jl)
           zzp=((paphp1(jl,klevp1)-paphp1(jl,jk))/                     &
                         (paphp1(jl,klevp1)-paphp1(jl,ikb)))
           zzp=MERGE(zzp**2,zzp,ktype(jl).EQ.3)
           pmfu(jl,jk)=pmfu(jl,ikb)*zzp
           pmfus(jl,jk)=pmfus(jl,ikb)*zzp
           pmfuq(jl,jk)=pmfuq(jl,ikb)*zzp
           pmful(jl,jk)=pmful(jl,ikb)*zzp
        END IF
125  END DO
!

!---wiso-code

! Correction of fluxes below the cloud base
     DO jt=1,kwiso
       DO jl=1,kproma
         IF(ldcum(jl).AND.jk.GT.kcbot(jl)) THEN
           ikb=kcbot(jl)
           zzp=((paphp1(jl,klevp1)-paphp1(jl,jk))/                     &
                         (paphp1(jl,klevp1)-paphp1(jl,ikb)))
           zzp=MERGE(zzp**2,zzp,ktype(jl).EQ.3)
           pwisomfuq(jl,jk,jt)=pwisomfuq(jl,ikb,jt)*zzp
           pwisomful(jl,jk,jt)=pwisomful(jl,ikb,jt)*zzp
         END IF
       END DO
     END DO

!---wiso-code-end

!re     DO 1254 jt=1,ktrac
!re!DIR$ IVDEP
!re!OCL NOVREC
!re        DO 1252 jl=1,kproma
!re           IF(ldcum(jl).AND.jk.GT.kcbot(jl)) THEN
!re              ikb=kcbot(jl)
!re              zzp=(paphp1(jl,klevp1)-paphp1(jl,jk))/                   &
!re                            (paphp1(jl,klevp1)-paphp1(jl,ikb))
!re              zzp=MERGE(zzp**2,zzp,ktype(jl).EQ.3)
!re              pmfuxt(jl,jk,jt)=pmfuxt(jl,ikb,jt)*zzp
!re           ENDIF
!re1252    END DO
!re1254 END DO
!re!

130 END DO
!
!
!*    2.            CALCULATE RAIN/SNOW FALL RATES
!*                  CALCULATE MELTING OF SNOW
!*                  CALCULATE EVAPORATION OF PRECIP
!                   -------------------------------
!
200 CONTINUE
  DO 210 jl=1,kproma
     prfl(jl)=0._dp
     psfl(jl)=0._dp
     prain(jl)=0._dp
210 END DO

!---wiso-code

! Set dummy field to zero
  DO jk=1,klev 
    DO jl=1,kproma
      pmelt_tmp(jl,jk)=0._dp
    END DO
  END DO
  
!---wiso-code-end

  DO 220 jk=ktopm2,klev
     DO 215 jl=1,kproma
        IF(ldcum(jl)) THEN
           prain(jl)=prain(jl)+pdmfup(jl,jk)
           IF(pten(jl,jk).GT.tmelt) THEN
              prfl(jl)=prfl(jl)+pdmfup(jl,jk)+pdmfdp(jl,jk)
              IF(psfl(jl).GT.0._dp.AND.pten(jl,jk).GT.ztmelp2) THEN
                 zfac=zcons1*(1._dp+vtmpc2*pqen(jl,jk))                &
                             *(paphp1(jl,jk+1)-paphp1(jl,jk))
                 zsnmlt=MIN(psfl(jl),zfac*(pten(jl,jk)-ztmelp2))
                 pdpmel(jl,jk)=zsnmlt
!---wiso-code
!                Store amount of melted snow
                 pmelt_tmp(jl,jk)=zsnmlt/psfl(jl)
                 pmelt_tmp(jl,jk)=MERGE(1._dp,pmelt_tmp(jl,jk),abs(1._dp-pmelt_tmp(jl,jk)).lt.cwisosec)
!---wiso-code-end
                 psfl(jl)=psfl(jl)-zsnmlt
                 prfl(jl)=prfl(jl)+zsnmlt
              END IF
           ELSE
              psfl(jl)=psfl(jl)+pdmfup(jl,jk)+pdmfdp(jl,jk)
           END IF
        END IF
215  END DO


220 END DO
  DO 230 jl=1,kproma
     prfl(jl)=MAX(prfl(jl),0._dp)
     psfl(jl)=MAX(psfl(jl),0._dp)
     zpsubcl(jl)=prfl(jl)+psfl(jl)
!---wiso-code
! store values of rainfall and snowfall
     prfl_tmp(jl) = prfl(jl)
     psfl_tmp(jl) = psfl(jl) 
!    Store fraction rain/(total precipitation)
     zrfl_frac(jl)=prfl(jl)/max(1.e-20_dp,prfl(jl)+psfl(jl))  
!---wiso-code-end
230 END DO

!---wiso-code

! Set some temporary fields to zero
  DO jk=1,klev
    DO jl=1,kproma
      pprec_frac(jl,jk)=0._dp
      prain_tmp(jl,jk)=0._dp
    END DO
  END DO

!---wiso-code-end

  DO 240 jk=ktopm2,klev
     DO 235 jl=1,kproma
       IF(ldcum(jl).AND.jk.GE.kcbot(jl).AND.zpsubcl(jl).GT.1.e-20_dp)  &
                                                                  THEN
           zrfl=zpsubcl(jl)
           zrnew=(MAX(0._dp,SQRT(zrfl/zcucov)-                         &
                        cevapcu(jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))*   &
                        MAX(0._dp,pqsen(jl,jk)-pqen(jl,jk))))**2*zcucov
           zrmin=zrfl-zcucov*MAX(0._dp,0.8_dp*pqsen(jl,jk)-pqen(jl,jk))&
                        *zcons2*(paphp1(jl,jk+1)-paphp1(jl,jk))
           zrnew=MAX(zrnew,zrmin)
           zrfln=MAX(zrnew,0._dp)
!---wiso-code
! Store (new amount of precipitation)/(old amount of prec.) in *pprec_frac*
! Calculate new amount of rain only and store it in *prain_tmp*
! (prain_tmp is needed to calculate isotope equilibrium of rain and
!  surrounding vapour for each vertical layer in CUXTEQU)
           pprec_frac(jl,jk)=zrfln/zrfl
           prain_tmp(jl,jk)=zrfln*zrfl_frac(jl)  
!---wiso-code-end
           zdrfl=MIN(0._dp,zrfln-zrfl)
           pdmfup(jl,jk)=pdmfup(jl,jk)+zdrfl
           zpsubcl(jl)=zrfln
       END IF
235  END DO
240 END DO
  DO 250 jl=1,kproma
     zrsum=prfl(jl)+psfl(jl)
     zdpevap=zpsubcl(jl)-zrsum
     prfl(jl)=prfl(jl)+zdpevap*prfl(jl)*(1._dp/MAX(1.e-20_dp,zrsum))
     psfl(jl)=psfl(jl)+zdpevap*psfl(jl)*(1._dp/MAX(1.e-20_dp,zrsum))
250 END DO
!
  RETURN
END SUBROUTINE h2oiso_cuflx



!======================================================================

SUBROUTINE h2oiso_cudtdq(kproma, kbdim, klev, klevp1, ktopm2, ldcum, &!re, ktrac,   &
                  paphp1,   pten,     ptte,     pqte,                  &
                  pxtec, &!re,    pxtte,    pmfuxt,   pmfdxt,                &
                  pmfus,    pmfds,    pmfuq,    pmfdq,                 &
                  pmful,    pdmfup,   pdmfdp,   plude,                 &
                  pdpmel,   prfl,     psfl,                            &
                  pcpen,    pqtec,    pqude,                           &
                  prsfc,    pssfc,    paprc,    paprs, delta_time)
!
!
!**** *CUDTDQ* - UPDATES T AND Q TENDENCIES, PRECIPITATION RATES
!                DOES GLOBAL DIAGNOSTICS
!
!          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
!
!**   INTERFACE.
!     ----------
!
!          *CUDTDQ* IS CALLED FROM *CUMASTR*
!
!!$USE mo_kind,         ONLY : dp
!!$USE mo_constants,    ONLY : alv, als, alf, tmelt, g
!!$USE mo_tracer,       ONLY : trlist
!!$USE mo_time_control, ONLY : delta_time
!
IMPLICIT NONE
!
INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, ktopm2!re, ktrac
!
REAL(DP),INTENT(IN):: delta_time
LOGICAL  llo1
!
REAL(dp) :: ptte(kbdim,klev),        pqte(kbdim,klev),                 &
            pten(kbdim,klev),        paphp1(kbdim,klevp1),             &
            paprc(kbdim),            paprs(kbdim),                     &
            prsfc(kbdim),            pssfc(kbdim)
REAL(dp) :: pmfus(kbdim,klev),       pmfds(kbdim,klev),                &
            pmfuq(kbdim,klev),       pmfdq(kbdim,klev),                &
            pmful(kbdim,klev),       plude(kbdim,klev),                &
            pdmfup(kbdim,klev),      pdmfdp(kbdim,klev),               &
            pqtec(kbdim,klev),       pqude(kbdim,klev),                &
            pxtec(kbdim,klev),       prfl(kbdim)
REAL(dp) :: pdpmel(kbdim,klev),      psfl(kbdim)
REAL(dp) :: pcpen(kbdim,klev)
LOGICAL  :: ldcum(kbdim)
!
REAL(dp) :: zmelt(kbdim)
REAL(dp) :: zsheat(kbdim)
!reREAL(dp) :: pxtte(kbdim,klev,ktrac), pmfuxt(kbdim,klev,ktrac),         &
!re            pmfdxt(kbdim,klev,ktrac)
!
REAL(dp) :: zrcpm ! reciprocal value of specific heat of moist air
INTEGER  :: jl, jk, jt
REAL(dp) :: zdiagt, zalv, zdtdt, zdqdt, zdxtdt
!
!----------------------------------------------------------------------
!
!*    1.0          SPECIFY PARAMETERS
!                  ------------------
!
100 CONTINUE
  zdiagt=delta_time
!
!
!----------------------------------------------------------------------
!
!*    2.0          INCREMENTATION OF T AND Q TENDENCIES
!                  ------------------------------------
!
200 CONTINUE
  DO 210 jl=1,kproma
     zmelt(jl)=0._dp
     zsheat(jl)=0._dp
210 END DO
!
  DO 250 jk=ktopm2,klev
!
     IF(jk.LT.klev) THEN
        DO 220 jl=1,kproma
           IF(ldcum(jl)) THEN
              llo1=(pten(jl,jk)-tmelt).GT.0._dp
              zalv=MERGE(alv,als,llo1)
              zrcpm=1._dp/pcpen(jl,jk)
              zdtdt=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*zrcpm*         &
                                  (pmfus(jl,jk+1)-pmfus(jl,jk)+        &
                                   pmfds(jl,jk+1)-pmfds(jl,jk)-        &
                                   alf*pdpmel(jl,jk)-                  &
                                   zalv*(pmful(jl,jk+1)-pmful(jl,jk)-  &
                                   plude(jl,jk)-                       &
                                  (pdmfup(jl,jk)+pdmfdp(jl,jk))))
              ptte(jl,jk)=ptte(jl,jk)+zdtdt
              zdqdt=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*               &
                                  (pmfuq(jl,jk+1)-pmfuq(jl,jk)+        &
                                   pmfdq(jl,jk+1)-pmfdq(jl,jk)+        &
                                   pmful(jl,jk+1)-pmful(jl,jk)-        &
                                   plude(jl,jk)-                       &
                                  (pdmfup(jl,jk)+pdmfdp(jl,jk)))
              pqte(jl,jk)=pqte(jl,jk)+zdqdt
              pxtec(jl,jk)=(g/(paphp1(jl,jk+1)-                        &
                            paphp1(jl,jk)))*plude(jl,jk)
              pqtec(jl,jk)=(g/(paphp1(jl,jk+1)-                        &
                            paphp1(jl,jk)))*pqude(jl,jk)
              zsheat(jl)=zsheat(jl)+zalv*(pdmfup(jl,jk)+pdmfdp(jl,jk))
              zmelt(jl)=zmelt(jl)+pdpmel(jl,jk)
           END IF
220     END DO
!
!re        IF (trlist% anyconv /= 0) THEN
!re            DO 2204 jt=1,ktrac
!re               IF (trlist% ti(jt)% nconv == 1) THEN
!re                 DO 2202 jl=1,kproma
!re                    IF(ldcum(jl)) THEN
!re                      zdxtdt=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))        &
!re                                  *(pmfuxt(jl,jk+1,jt)-pmfuxt(jl,jk,jt) &
!re                                   +pmfdxt(jl,jk+1,jt)-pmfdxt(jl,jk,jt))
!re                      pxtte(jl,jk,jt)=pxtte(jl,jk,jt)+zdxtdt
!re                    ENDIF
!re 2202            END DO
!re               ENDIF
!re 2204       END DO
!re         ENDIF
!
!
     ELSE
        DO 230 jl=1,kproma
           IF(ldcum(jl)) THEN
              llo1=(pten(jl,jk)-tmelt).GT.0._dp
              zalv=MERGE(alv,als,llo1)
              zrcpm=1._dp/pcpen(jl,jk)
              zdtdt=-(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*zrcpm*        &
                     (pmfus(jl,jk)+pmfds(jl,jk)+alf*pdpmel(jl,jk)-     &
                                 zalv*(pmful(jl,jk)+pdmfup(jl,jk)      &
                                +pdmfdp(jl,jk)+plude(jl,jk)))
              ptte(jl,jk)=ptte(jl,jk)+zdtdt
              zdqdt=-(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*              &
                        (pmfuq(jl,jk)+pmfdq(jl,jk)+plude(jl,jk)+       &
                        (pmful(jl,jk)+pdmfup(jl,jk)+pdmfdp(jl,jk)))
              pqte(jl,jk)=pqte(jl,jk)+zdqdt
              pxtec(jl,jk)=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))         &
                           *plude(jl,jk)
              pqtec(jl,jk)=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))         &
                           *pqude(jl,jk)
              zsheat(jl)=zsheat(jl)+zalv*(pdmfup(jl,jk)+pdmfdp(jl,jk))
              zmelt(jl)=zmelt(jl)+pdpmel(jl,jk)
           END IF
230     END DO


!
!re         IF (trlist% anyconv /= 0) THEN
!re           DO 2304 jt=1,ktrac
!re               IF (trlist% ti(jt)% nconv == 1) THEN
!re                 DO 2302 jl=1,kproma
!re                    IF(ldcum(jl)) THEN
!re                       zdxtdt=-(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))      &
!re                              *(pmfuxt(jl,jk,jt)+pmfdxt(jl,jk,jt))
!re                       pxtte(jl,jk,jt)=pxtte(jl,jk,jt)+zdxtdt
 !re                   ENDIF
!re 2302            END DO
!re               END IF
!re 2304       END DO
!re         ENDIF
!
     END IF
!
250 END DO
!
!
!---------------------------------------------------------------------
!
!      3.          UPDATE SURFACE FIELDS
!                  ---------------------
!
300 CONTINUE
  DO 310 jl=1,kproma
     prsfc(jl)=prfl(jl)
     pssfc(jl)=psfl(jl)
     paprc(jl)=paprc(jl)+zdiagt*(prfl(jl)+psfl(jl))
     paprs(jl)=paprs(jl)+zdiagt*psfl(jl)
310 END DO
!
  RETURN
END SUBROUTINE h2oiso_cudtdq


!======================================================================



SUBROUTINE h2oiso_cuwisodq(kproma, kbdim, klev, klevp1, ktopm2, ldcum, kwiso, &
                  paphp1,   pwisoqte,                                  &
                  pwisoxtec,                                           &
                  pwisomfuq,pwisomfdq,                                 &
                  pwisomful,pwisodmfup,pwisodmfdp,pwisolude,           &
                  pwisoqtec,pwisoqude, delta_time)
!
! DESCRIPTION:
!
! UPDATES WATER ISOTOPE TENDENCIES, PRECIPITATION RATES DOES GLOBAL DIAGNOSTICS
!
! METHOD:
!
! *CUWISODQ* IS CALLED FROM *CUMASTR*
!
! AUTHORS:
!
! G.HOFFMANN, MPI MET, HAMBURG, 1992 
! ADAPTED TO F90: M. WERNER, MPI BGC, JENA, 2004 
! ADAPTED TO ECHAM5: M. WERNER, AWI, BREMERHAVEN, 2009
!

!!$USE mo_kind,         ONLY : dp
!!$USE mo_constants,    ONLY : g
!!$USE mo_time_control, ONLY : delta_time

IMPLICIT NONE

INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, ktopm2, kwiso
REAL(dp), INTENT (IN) :: delta_time

REAL(dp) :: pwisoqte(kbdim,klev,kwiso),                                &
            paphp1(kbdim,klevp1)

REAL(dp) :: pwisomfuq(kbdim,klev,kwiso), pwisomfdq(kbdim,klev,kwiso),  &
            pwisomful(kbdim,klev,kwiso), pwisolude(kbdim,klev,kwiso),  &
            pwisodmfup(kbdim,klev,kwiso),pwisodmfdp(kbdim,klev,kwiso), &
            pwisoqtec(kbdim,klev,kwiso), pwisoqude(kbdim,klev,kwiso),  &
            pwisoxtec(kbdim,klev,kwiso)

LOGICAL  :: ldcum(kbdim)
!
REAL(dp) :: zwisodqdt

INTEGER  :: jl, jk, jt

!  Executable statements 

  DO jk=ktopm2,klev
!
    IF(jk.LT.klev) THEN
      DO jt=1,kwiso
        DO jl=1,kproma
          IF(ldcum(jl)) THEN
            zwisodqdt=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*             &
                      (pwisomfuq(jl,jk+1,jt)-pwisomfuq(jl,jk,jt)+      &
                       pwisomfdq(jl,jk+1,jt)-pwisomfdq(jl,jk,jt)+      &
                       pwisomful(jl,jk+1,jt)-pwisomful(jl,jk,jt)-      &
                       pwisolude(jl,jk,jt)-                            &
                      (pwisodmfup(jl,jk,jt)+pwisodmfdp(jl,jk,jt)))
            pwisoqte(jl,jk,jt)=pwisoqte(jl,jk,jt)+zwisodqdt
            pwisoxtec(jl,jk,jt)=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))    &
                                *pwisolude(jl,jk,jt)
            pwisoqtec(jl,jk,jt)=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))    &
                                *pwisoqude(jl,jk,jt)
          ENDIF
        END DO
      END DO
    ELSE
      DO jt=1,kwiso
        DO jl=1,kproma
          IF(ldcum(jl)) THEN
            zwisodqdt=-(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*            &
                       (pwisomfuq(jl,jk,jt)+pwisomfdq(jl,jk,jt)+       &
                        pwisolude(jl,jk,jt)+                           &
                       (pwisomful(jl,jk,jt)+pwisodmfup(jl,jk,jt)+      &
                        pwisodmfdp(jl,jk,jt)))
            pwisoqte(jl,jk,jt)=pwisoqte(jl,jk,jt)+zwisodqdt
            pwisoxtec(jl,jk,jt)=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))    &
                                *pwisolude(jl,jk,jt)
            pwisoqtec(jl,jk,jt)=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))    &
                                *pwisoqude(jl,jk,jt)
          ENDIF
        END DO
      END DO
    END IF
!
  END DO


  RETURN
END SUBROUTINE h2oiso_cuwisodq




!======================================================================







SUBROUTINE h2oiso_cuwisoequ(kproma, kbdim, klev, klevp1, ktopm2, ldcum, kwiso, &
                     kctop,          kcbot,       paphp1,              &
                     ptp1,           ptte,                             & 
                     pqp1,           pqte,                             &
                     pwisoqp1,       pwisoqte,                         &
                     ptu,            pqu,         pwisoqu,             &
                     pten,                                             &
                     prfl_tmp,       psfl_tmp,                         &
                     pwisorfl,       pwisosfl,                         &
                     pwisodmfup,     pwisodmfdp,                       &
                     pmelt_tmp,                                        &
                     pprec_frac,     prain_tmp,                        &
                     pwisorsfc,      pwisossfc,                        &
                     pwisoaprc,      pwisoaprs,                        &
                     time_step_len,  delta_time)
                     
!
! PURPOSE:
!
! GET ISOTOPE CONCENTRATION IN RAIN INTO EQUILIBRIUM WITH 
! THE CONCENTRATION IN VAPOUR.
!
! INTERFACE:
!
! THIS SUBROUTINE IS CALLED FROM
! *CUMASTR*
!
! AUTHORS:
!
! G.HOFFMANN, MPI MET, HAMBURG, 1992 
! ADAPTED TO F90: M. WERNER, MPI BGC, JENA, 2004 
! ADAPTED TO ECHAM5: M. WERNER, AWI, BREMERHAVEN, 2009
!


!!$USE mo_kind,           ONLY: dp
!!$USE mo_constants,      ONLY: tmelt, g, vtmpc1
!!$USE mo_convect_tables, ONLY: tlucua,                                   &
!!$                             tlucub,                                   &
!!$                             jptlucu1, jptlucu2,                       &
!!$                             lookuperror, lookupoverflow
!!$USE mo_time_control,   ONLY: time_step_len, delta_time
USE messy_main_tools_wiso, ONLY: &!reoldtalph1, talph2, talph3,        &
                             talphal1, talphal2, talphal3,             &
                             thumwiso1, thumwiso2, twisoeqcu,          &
                             tdifrel,   tdifexp,                       &
                             tnat, cthomi,                &
                             cwisomin,  cwisosec, i_HH18O, i_HDO, i_HHO

IMPLICIT NONE
!  
INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, ktopm2, kwiso
REAL(dp),INTENT (IN) :: time_step_len, delta_time

REAL(dp):: paphp1(kbdim,klevp1),                                       &
           ptp1(kbdim,klev),           ptte(kbdim,klev),               &             
           pqp1(kbdim,klev),           pqte(kbdim,klev),               &
           pwisoqp1(kbdim,klev,kwiso), pwisoqte(kbdim,klev,kwiso),     &
           ptu(kbdim,klev),            pqu(kbdim,klev),                &
           pwisoqu(kbdim,klev,kwiso),                                  &
           pten(kbdim,klev),                                           &
           prfl_tmp(kbdim),            psfl_tmp(kbdim),                &
           pwisorfl(kbdim,kwiso),      pwisosfl(kbdim,kwiso),          &
           pwisodmfup(kbdim,klev,kwiso), pwisodmfdp(kbdim,klev,kwiso), &
           pmelt_tmp(kbdim,klev),                                      &
           pprec_frac(kbdim,klev),    prain_tmp(kbdim,klev),           &
           pwisorsfc(kbdim,kwiso),    pwisossfc(kbdim,kwiso),          &
           pwisoaprc(kbdim,kwiso),    pwisoaprs(kbdim,kwiso)

INTEGER :: kctop(kbdim), kcbot(kbdim)

LOGICAL  :: ldcum(kbdim)
!
REAL(dp) :: zrelq(kbdim,klev),                                         &
            zwisopsubcl(kbdim,kwiso)

REAL(dp) :: zdiagt,  zwisofracliq, &!reold,  zwisofrac,                  &
            zes,      zcor,       zqsat,                               &
            zwisorfl, zwisorfln,  zwisodrfl,                           &
            zt,       zqv,        zwisoqv,    zql,      zwisoql,       &
            zhumass,  zqliq,      zqice,    zquot,  zdelt, &
            zpromill, zdmaxo18,   zdmaxdeu,   zdabsmax,                &
            zqob,     zwisoqob,   zdeltaob,                            &
            zqun,     zwisoqun,   zdeltaun,                            &
            zqtop,    zwisoqtop,  zdeltatop,                           &
            zqmiddle, zwisoqmiddle,zdeltamiddle,                       &
            zqbottom, zwisoqbottom,zdeltabottom,                       &
            zqsum,    zwisoqsum,                                       &
            zdeltasum,                                                 &
            zwisosum_tmp,                                              &
            zwisorsum,zwisodpevap

INTEGER  :: jl, jk, jt, it, jkoben, jtf1, jtf2

INTEGER  :: kmix(kbdim), knull(kbdim) 

LOGICAL  :: lomix1, lomix2, lomix(kbdim), lonot

REAL(dp):: qwisorfl(kproma,kwiso), qwisosfl(kproma,kwiso)!retest

! SPECIFY PARAMETERS

  zdiagt=delta_time

  pwisorfl(:,:) = 0._dp
  pwisosfl(:,:) = 0._dp
  zwisopsubcl(:,:) = 0._dp

! calculation of the relative humidity below the cloud base

!re#ifndef NOLOOKUP
  lookupoverflow = .FALSE.
!re#endif

  DO jk=1,klev
    DO jl=1,kproma

      ptp1(jl,jk)=ptp1(jl,jk)+ptte(jl,jk)*time_step_len
      pqp1(jl,jk)=pqp1(jl,jk)+pqte(jl,jk)*time_step_len
       
      IF (jk.ge.kcbot(jl)) THEN

        it = NINT(ptp1(jl,jk)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zes=tlucua(it)/paphp1(jl,jk)
        zes=MIN(0.5_dp,zes)
        zcor=1._dp/(1._dp-vtmpc1*zes)
        zqsat=zes*zcor
        IF (zqsat.lt.cwisomin) THEN
          zrelq(jl,jk)=0._dp
        ELSE IF (zqsat .LT. pqp1(jl,jk)) THEN
          zrelq(jl,jk)=1._dp
        ELSE
          zrelq(jl,jk)=pqp1(jl,jk)/zqsat
          IF ((1._dp-zrelq(jl,jk)).LT.cwisosec) zrelq(jl,jk)=1._dp
        ENDIF

      ELSE
        zrelq(jl,jk)=1._dp
      ENDIF
      
    END DO
  END DO

  IF (lookupoverflow) RETURN !re CALL lookuperror ('cuwisoequ')

! getting into equilibrium rain and vapour in and below the cloud

  DO jt=1,kwiso

! first part: calculate isotope precipitation like normal precipitation (done in *cuflx*)
 
    DO jk=ktopm2,klev
      DO jl=1,kproma
        IF (ldcum(jl)) THEN
          IF(pten(jl,jk).gt.tmelt) THEN
            pwisorfl(jl,jt)=pwisorfl(jl,jt)+pwisodmfup(jl,jk,jt)+pwisodmfdp(jl,jk,jt)
!           corrected rainfall for amount of melted water isotopes            
            pwisorfl(jl,jt)=pwisorfl(jl,jt)+pmelt_tmp(jl,jk)*pwisosfl(jl,jt)
            pwisosfl(jl,jt)=pwisosfl(jl,jt)-pmelt_tmp(jl,jk)*pwisosfl(jl,jt)
          ELSE
            pwisosfl(jl,jt)=pwisosfl(jl,jt)+pwisodmfup(jl,jk,jt)+pwisodmfdp(jl,jk,jt)
          END IF
        ENDIF
      END DO
    END DO

! check for negative rainfall and snowfall values
    DO jl=1,kproma
      if(prfl_tmp(jl).le.0.0_dp) pwisorfl(jl,jt)=0.0_dp!re _dp
      if(psfl_tmp(jl).le.0.0_dp) pwisosfl(jl,jt)=0.0_dp!re _dp
      zwisopsubcl(jl,jt)=pwisorfl(jl,jt)+pwisosfl(jl,jt)
    END DO
            
    DO jk=ktopm2,klev
      DO jl=1,kproma
        IF (ldcum(jl).AND.jk.ge.kcbot(jl).AND.ABS(zwisopsubcl(jl,jt)).gt.1.e-20_dp) THEN 

          zwisorfl=zwisopsubcl(jl,jt)
          zwisorfln=zwisorfl*pprec_frac(jl,jk)
          zwisodrfl=MIN(0._dp,zwisorfln-zwisorfl)
          pwisoqte(jl,jk,jt)=pwisoqte(jl,jk,jt)-zwisodrfl*g/(paphp1(jl,jk+1)-paphp1(jl,jk))
          zwisosum_tmp=pwisorfl(jl,jt)+pwisosfl(jl,jt)
          IF (ABS(zwisosum_tmp).LT.1.e-20_dp) THEN
            IF (zwisosum_tmp.GE.0._dp) THEN
              zwisosum_tmp=+1.e-20_dp
            ELSE 
              zwisosum_tmp=-1.e-20_dp
            ENDIF
          ENDIF
         pwisorfl(jl,jt)=zwisorfln*pwisorfl(jl,jt)/zwisosum_tmp          !retest  
       pwisosfl(jl,jt)=zwisorfln*pwisosfl(jl,jt)/zwisosum_tmp          !retest   

 !         qwisorfl(jl,jt)=zwisorfln*pwisorfl(jl,jt)/zwisosum_tmp !retest
 !         qwisosfl(jl,jt)=zwisorfln*pwisosfl(jl,jt)/zwisosum_tmp !retest


! second part: bring isotopes in precipitation in equilibrium with the surrounding water vapour  
! (only isotopes in rain water go into equilibrium)    
          IF (jk.lt.kcbot(jl).and.jk.ge.kctop(jl)) THEN ! if-condition = cloud layer
            zt=ptu(jl,jk)
            zqv=pqu(jl,jk)
            zwisoqv=pwisoqu(jl,jk,jt)
          ELSE                                          ! subcloud layer
            zt=ptp1(jl,jk)
            zqv=pqp1(jl,jk)
            zwisoqv=pwisoqp1(jl,jk,jt)+pwisoqte(jl,jk,jt)*time_step_len
          ENDIF
               
            zql=prain_tmp(jl,jk)/((paphp1(jl,jk+1)-paphp1(jl,jk))/g/time_step_len)     
            zwisoql=pwisorfl(jl,jt)/((paphp1(jl,jk+1)-paphp1(jl,jk))/g/time_step_len)
              
!         calculate fractionation coefficient for liquid-vapour phase change
!reold          zwisofrac=exp(talph1(jt)/(zt**2)+talph2(jt)/zt+talph3(jt))
          zwisofracliq=exp(talphal1(jt)/(zt**2)+talphal2(jt)/zt+talphal3(jt))

          IF (zqv.gt.cwisomin.and.zwisoqv.gt.cwisomin) THEN
      
! effective fractionation  below the cloud base
              zhumass=thumwiso1+thumwiso2*zrelq(jl,jk)
              !reold zwisofrac=zwisofrac*zhumass/(zwisofrac*(zhumass-1._dp)*(tdifrel(jt)**tdifexp)+1._dp)
                zwisofracliq=zwisofracliq*zhumass/(zwisofracliq*(zhumass-1._dp)*(tdifrel(jt)**tdifexp)+1._dp)

              IF (zt.GT.tmelt) THEN   ! liquid and ice fraction
                zqliq=1._dp
              ELSEIF (zt.LT.cthomi) THEN
                zqliq=0._dp
              ELSE
                zqice=(tmelt-zt)/(tmelt-cthomi)
                zqliq=1._dp-zqice
              ENDIF
      
!reold              zquot=zqv+zwisofrac*twisoeqcu*zqliq*zql
              zquot=zqv+zwisofracliq*twisoeqcu*zqliq*zql

            zdelt=0._dp
!reold              IF (zquot.gt.cwisomin) zdelt=twisoeqcu*zqliq*(zwisofrac*zql*zwisoqv-zwisoql*zqv)/zquot
              IF (zquot.gt.cwisomin) zdelt=twisoeqcu*zqliq*(zwisofracliq*zql*zwisoqv-zwisoql*zqv)/zquot

       pwisorfl(jl,jt)=pwisorfl(jl,jt)+zdelt*((paphp1(jl,jk+1)-paphp1(jl,jk))/g/time_step_len)!retest       
       pwisoqte(jl,jk,jt)=pwisoqte(jl,jk,jt)-zdelt/time_step_len!retest       

! for negative humidity values: set vapour in simple equilibrium to convective rain      
            ELSE
              
              zdelt=tnat(jt)
                    IF (zql.gt.cwisomin) zdelt=zwisoql/zql
                    
!reold         zwisoqv=(1._dp/zwisofrac)*zdelt*zqv
                zwisoqv=(1._dp/zwisofracliq)*zdelt*zqv
                pwisoqte(jl,jk,jt)=(zwisoqv-pwisoqp1(jl,jk,jt))/time_step_len

            ENDIF
            
         zwisopsubcl(jl,jt)=pwisorfl(jl,jt)+pwisosfl(jl,jt)          !retest  
 !         zwisopsubcl(jl,jt)=qwisorfl(jl,jt)+qwisosfl(jl,jt)!retest 

        ENDIF
        
      END DO
    END DO
  
  END DO

  DO jt=1,kwiso
    DO jl=1,kproma
      zwisorsum=pwisorfl(jl,jt)+pwisosfl(jl,jt)
      zwisodpevap=zwisopsubcl(jl,jt)-zwisorsum
      pwisorfl(jl,jt)=pwisorfl(jl,jt)+zwisodpevap*pwisorfl(jl,jt)*(1._dp/MAX(1.e-20_dp,zwisorsum))
      pwisosfl(jl,jt)=pwisosfl(jl,jt)+zwisodpevap*pwisosfl(jl,jt)*(1._dp/MAX(1.e-20_dp,zwisorsum))
    END DO
  END DO

! do additional vertical mixing of water isotopes
! if a certain vertical gradient of delta values is exceeded
! do mixing only if flag *lonot* is set
  lonot=.TRUE.
  IF (lonot) THEN
      
! set some constants 
    zpromill=1000._dp
!    zdmaxo18=1000._dp
!    zdmaxdeu=8000._dp  
    zdmaxo18=150._dp
    zdmaxdeu=1200._dp  
    zdabsmax=1000._dp
    jkoben=1
    kmix(:) = 0
    knull(:) = 0
    jtf1 = 0
    jtf2 = 0
      
! jtf1=tracerno. of h2-18o, jtf2=tracerno. of hdo 
    DO jt=1,kwiso
      IF (jt == i_HH18O) jtf1=jt   
      IF (jt == i_HDO) jtf2=jt
    END DO

! loop over all levels, isotopes and longitudes
    DO jk=jkoben+1,klev-1
      lomix(:) = .FALSE.
      DO jt=1,kwiso
!CDIR$ IVDEP
        DO jl=1,kproma
        
! check, if delta-o18-gradient between two levels is gt. zdmaxo18
          IF (jtf1.gt.0) THEN
            zwisoqob=pwisoqp1(jl,jk-1,jtf1)+pwisoqte(jl,jk-1,jtf1)*time_step_len
            zqob=pqp1(jl,jk-1)
            zdeltaob=0._dp
            IF ((zwisoqob.gt.cwisomin).and.(zqob.gt.cwisomin)) zdeltaob=(zwisoqob/(zqob*tnat(jtf1))-1._dp)*zpromill
            zwisoqun=pwisoqp1(jl,jk,jtf1)+pwisoqte(jl,jk,jtf1)*time_step_len
            zqun=pqp1(jl,jk)
            zdeltaun=0._dp
            IF ((zwisoqun.gt.cwisomin).and.(zqun.gt.cwisomin)) zdeltaun=(zwisoqun/(zqun*tnat(jtf1))-1._dp)*zpromill
          ELSE
            zdeltaob=0._dp
            zdeltaun=0._dp
          ENDIF
          lomix1=abs(zdeltaob-zdeltaun).gt.zdmaxo18

! check, if delta-hdo-gradient between two levels is gt. zdmaxdeu
          IF (jtf2.gt.0) THEN
            zwisoqob=pwisoqp1(jl,jk-1,jtf2)+pwisoqte(jl,jk-1,jtf2)*time_step_len
            zqob=pqp1(jl,jk-1)
            zdeltaob=0._dp
            IF ((zwisoqob.gt.cwisomin).and.(zqob.gt.cwisomin)) zdeltaob=(zwisoqob/(zqob*tnat(jtf2))-1._dp)*zpromill
            zwisoqun=pwisoqp1(jl,jk,jtf2)+pwisoqte(jl,jk,jtf2)*time_step_len
            zqun=pqp1(jl,jk)
            zdeltaun=0._dp
            IF ((zwisoqun.gt.cwisomin).and.(zqun.gt.cwisomin)) zdeltaun=(zwisoqun/(zqun*tnat(jtf2))-1._dp)*zpromill
          ELSE
            zdeltaob=0._dp
            zdeltaun=0._dp
          ENDIF
          lomix2=abs(zdeltaob-zdeltaun).gt.zdmaxdeu
          
! check exactly once if one of the two gradients is too big
! if yes: mix three levels (mix only water isotopes, not normal water)    

          IF ((lomix1.or.lomix2).and.(jt == i_HHO)) lomix(jl)=.TRUE.

          IF (lomix(jl).and. (jt /= i_HHO)) THEN
           
! get old delta values of three boxes (lev-1,lev,lev+1)
! (neglect any box with negative water values...)           
            zwisoqtop=pwisoqp1(jl,jk-1,jt)+pwisoqte(jl,jk-1,jt)*time_step_len
            zqtop=pqp1(jl,jk-1)
            zdeltatop=0._dp
            IF ((zwisoqtop.gt.cwisomin).and.(zqtop.gt.cwisomin))               & 
               zdeltatop=(zwisoqtop/(zqtop*tnat(jt))-1._dp)*zpromill
          
            zwisoqmiddle=pwisoqp1(jl,jk,jt)+pwisoqte(jl,jk,jt)*time_step_len
            zqmiddle=pqp1(jl,jk)
            zdeltamiddle=0._dp
            IF ((zwisoqmiddle.gt.cwisomin).and.(zqmiddle.gt.cwisomin))         &
               zdeltamiddle=(zwisoqmiddle/(zqmiddle*tnat(jt))-1._dp)*zpromill

            zwisoqbottom=pwisoqp1(jl,jk+1,jt)+pwisoqte(jl,jk+1,jt)*time_step_len
            zqbottom=pqp1(jl,jk+1)
            zdeltabottom=0._dp
            IF ((zwisoqbottom.gt.cwisomin).and.(zqbottom.gt.cwisomin))         &
               zdeltabottom=(zwisoqbottom/(zqbottom*tnat(jt))-1._dp)*zpromill

! Calculate mean delta value of the three boxes (exclude boxes with too high delta-values)
            zqsum=0._dp
            zwisoqsum=0._dp
              
            IF (abs(zdeltatop).lt.zdabsmax) THEN 
              zqsum=zqsum+zqtop
              zwisoqsum=zwisoqsum+zwisoqtop
            ENDIF

            IF (abs(zdeltamiddle).lt.zdabsmax) THEN               
              zqsum=zqsum+zqmiddle
              zwisoqsum=zwisoqsum+zwisoqmiddle
            ENDIF
         
            IF (abs(zdeltabottom).lt.zdabsmax) THEN
              zqsum=zqsum+zqbottom
              zwisoqsum=zwisoqsum+zwisoqbottom
            ENDIF
        
            zdeltasum=0._dp
            IF ((zwisoqsum.gt.cwisomin).and.(zqsum.gt.cwisomin)) zdeltasum=(zwisoqsum/(zqsum*tnat(jt))-1._dp)*zpromill

! Calculate new isotope values of the three boxes
            zwisoqtop=(zdeltasum/zpromill+1)*zqtop*tnat(jt)
            zwisoqmiddle=(zdeltasum/zpromill+1)*zqmiddle*tnat(jt)
            zwisoqbottom=(zdeltasum/zpromill+1)*zqbottom*tnat(jt)

! Calculate new isotope tendencies
            pwisoqte(jl,jk-1,jt)=(zwisoqtop-pwisoqp1(jl,jk-1,jt))/time_step_len
            pwisoqte(jl,jk,jt)=(zwisoqmiddle-pwisoqp1(jl,jk,jt))/time_step_len
            pwisoqte(jl,jk+1,jt)=(zwisoqbottom-pwisoqp1(jl,jk+1,jt))/time_step_len

          ENDIF
          
        END DO
      END DO
    END DO
      
  ENDIF

! UPDATE SURFACE FIELDS
  DO jt=1,kwiso
    DO jl=1,kproma
      pwisorsfc(jl,jt)=pwisorfl(jl,jt)
      pwisossfc(jl,jt)=pwisosfl(jl,jt)
      pwisoaprc(jl,jt)=pwisoaprc(jl,jt)+zdiagt*(pwisorfl(jl,jt)+pwisosfl(jl,jt))
      pwisoaprs(jl,jt)=pwisoaprs(jl,jt)+zdiagt*pwisosfl(jl,jt)
    END DO
  END DO

  RETURN

END SUBROUTINE h2oiso_cuwisoequ












!======================================================================




SUBROUTINE h2oiso_cuentr(   kproma, kbdim, klev, klevp1, kk,                  &
           ptenh,    pqenh,    pqte,     paphp1,                       &
           klwmin,   ldcum,    ktype,    kcbot,    kctop0,             &
           ppbase,   pmfu,     pentr,    podetr,                       &
           khmin,    pgeoh,                                            &
           pdmfen,   pdmfde)
!
!          M.TIEDTKE         E.C.M.W.F.     12/89
!
!          PURPOSE.
!          --------
!          THIS ROUTINE CALCULATES ENTRAINMENT/DETRAINMENT RATES
!          FOR UPDRAFTS IN CUMULUS PARAMETERIZATION
!
!          INTERFACE
!          ---------
!
!          THIS ROUTINE IS CALLED FROM *CUASC*.
!          INPUT ARE ENVIRONMENTAL VALUES T,Q ETC
!          AND UPDRAFT VALUES T,Q ETC
!          IT RETURNS ENTRAINMENT/DETRAINMENT RATES
!
!          METHOD.
!          --------
!          S. TIEDTKE (1989)
!
!          EXTERNALS
!          ---------
!          NONE
!
!!$USE mo_kind,           ONLY: dp
!!$USE mo_constants,      ONLY: g, rd, vtmpc1
!!$USE mo_cumulus_flux,   ONLY: centrmax, cmfcmin
!
IMPLICIT NONE
!
INTEGER, INTENT (IN) :: kbdim, klev, klevp1, kproma, kk
!
REAL(dp) :: ptenh(kbdim,klev),       pqenh(kbdim,klev),                &
            paphp1(kbdim,klevp1),                                      &
            pmfu(kbdim,klev),        pqte(kbdim,klev),                 &
            pentr(kbdim),            ppbase(kbdim)
REAL(dp) :: podetr(kbdim,klev)
REAL(dp) :: pgeoh (kbdim,klev)
INTEGER  :: khmin (kbdim)
INTEGER  :: klwmin(kbdim),           ktype(kbdim),                     &
            kcbot(kbdim),            kctop0(kbdim)
LOGICAL  :: ldcum(kbdim)
!
REAL(dp) :: pdmfen(kbdim),           pdmfde(kbdim)
!
LOGICAL  :: llo1,llo2
!
INTEGER  :: jl, ikb, ikt, ikh, iklwmin
REAL(dp) :: zrg, zrrho, zdprho, zpmid, zentr, zentest, zzmzk, ztmzk    &
          , zorgde, zarg
!
!----------------------------------------------------------------------
!
!*    1.           CALCULATE ENTRAINMENT AND DETRAINMENT RATES
!                  -------------------------------------------
!
!re100 CONTINUE
!
!
!*    1.1          SPECIFY ENTRAINMENT RATES FOR SHALLOW CLOUDS
!                  --------------------------------------------
!
!re110 CONTINUE
!
!
!*    1.2          SPECIFY ENTRAINMENT RATES FOR DEEP CLOUDS
!                  -----------------------------------------
!
!re120 CONTINUE
  zrg=1._dp/g
  DO 125 jl=1,kproma
     ppbase(jl)=paphp1(jl,kcbot(jl))
     zrrho=(rd*ptenh(jl,kk+1)*(1._dp+vtmpc1*pqenh(jl,kk+1)))           &
                    /paphp1(jl,kk+1)
     zdprho=(paphp1(jl,kk+1)-paphp1(jl,kk))*zrg
     zpmid=0.5_dp*(ppbase(jl)+paphp1(jl,kctop0(jl)))
     zentr=pentr(jl)*pmfu(jl,kk+1)*zdprho*zrrho
     llo1=kk.LT.kcbot(jl).AND.ldcum(jl)
     pdmfde(jl)=MERGE(zentr,0._dp,llo1)
     llo2=llo1.AND.ktype(jl).EQ.2.AND.                                 &
                   (ppbase(jl)-paphp1(jl,kk).LT.0.2e5_dp.OR.           &
                    paphp1(jl,kk).GT.zpmid)
     pdmfen(jl)=MERGE(zentr,0._dp,llo2)
     iklwmin=MAX(klwmin(jl),kctop0(jl)+2)
     llo2=llo1 .AND. ktype(jl).EQ.3 .AND. kk .GE. iklwmin
     IF(llo2) pdmfen(jl)=zentr
     IF(llo2 .AND. pqenh(jl,kk+1).GT.1.E-5_dp) THEN
        pmfu(jl,kk+1) = MAX(pmfu(jl,kk+1),cmfcmin)
        zentest = MAX(pqte(jl,kk),0._dp)/pqenh(jl,kk+1)
        zentest = MIN(centrmax,zentest/(pmfu(jl,kk+1)*zrrho))
        pdmfen(jl) = zentr+zentest*pmfu(jl,kk+1)*zrrho*zdprho
     ENDIF
     llo2=llo1 .AND. ktype(jl).EQ.1 .AND.                              &
                   (kk.GE.iklwmin.OR.paphp1(jl,kk).GT.zpmid)
     IF(llo2) pdmfen(jl)=zentr
!
!    organized detrainment, detrainment starts at khmin
!
     llo2=llo1 .AND. ktype(jl).EQ.1
!re     ikb=kcbot(jl)
     podetr(jl,kk)=0._dp
     IF(llo2.AND.kk.LE.khmin(jl).AND.kk.GE.kctop0(jl)) THEN
        ikt=kctop0(jl)
        ikh=khmin(jl)
        IF(ikh.GT.ikt) THEN
           zzmzk  =-(pgeoh(jl,ikh)-pgeoh(jl,kk))*zrg
           ztmzk  =-(pgeoh(jl,ikh)-pgeoh(jl,ikt))*zrg
           zarg  =3.1415_dp*(zzmzk/ztmzk)*0.5_dp
           zorgde=TAN(zarg)*3.1415_dp*0.5_dp/ztmzk
           zdprho=(paphp1(jl,kk+1)-paphp1(jl,kk))*(zrg*zrrho)
           podetr(jl,kk)=MIN(zorgde,centrmax)*pmfu(jl,kk+1)*zdprho
        ENDIF
     ENDIF
125 END DO
!
  RETURN
END SUBROUTINE h2oiso_cuentr


!======================================================================

SUBROUTINE h2oiso_cuentrt(  kproma, kbdim, klev, klevp1, kk,                  &
           ptenh,    pqenh,    pqte,     paphp1,                       &
           klwmin,   ldcum,    ktype,    kcbot,    kctop0,             &
           ppbase,   pmfu,     pentr,                                  &
           pdmfen,   pdmfde)
!
!          M.TIEDTKE         E.C.M.W.F.     12/89
!
!          PURPOSE.
!          --------
!          THIS ROUTINE CALCULATES ENTRAINMENT/DETRAINMENT RATES
!          FOR UPDRAFTS IN CUMULUS PARAMETERIZATION
!
!          INTERFACE
!          ---------
!
!          THIS ROUTINE IS CALLED FROM *CUASC*.
!          INPUT ARE ENVIRONMENTAL VALUES T,Q ETC
!          AND UPDRAFT VALUES T,Q ETC
!          IT RETURNS ENTRAINMENT/DETRAINMENT RATES
!
!          METHOD.
!          --------
!          S. TIEDTKE (1989)
!
!          EXTERNALS
!          ---------
!          NONE
!
!reUSE mo_kind,           ONLY: dp
!reUSE mo_constants,      ONLY: g, rd, vtmpc1
!reUSE mo_cumulus_flux,   ONLY: centrmax, cmfcmin
!
!reIMPLICIT NONE
!
INTEGER, INTENT (IN) :: kbdim, klev, klevp1, kproma, kk
!
REAL(dp) :: ptenh(kbdim,klev),       pqenh(kbdim,klev),                &
            paphp1(kbdim,klevp1),                                      &
            pmfu(kbdim,klev),        pqte(kbdim,klev),                 &
            pentr(kbdim),            ppbase(kbdim)
INTEGER  :: klwmin(kbdim),           ktype(kbdim),                     &
            kcbot(kbdim),            kctop0(kbdim)
LOGICAL  :: ldcum(kbdim)
!
REAL(dp) :: pdmfen(kbdim),           pdmfde(kbdim)
!
LOGICAL  :: llo1,llo2
!
INTEGER  :: jl, iklwmin
REAL(dp) :: zrg, zrrho, zdprho, zpmid, zentr, zentest
          
!
!----------------------------------------------------------------------
!
!*    1.           CALCULATE ENTRAINMENT AND DETRAINMENT RATES
!                  -------------------------------------------
!
!re100 CONTINUE
!
!
!*    1.1          SPECIFY ENTRAINMENT RATES FOR SHALLOW CLOUDS
!                  --------------------------------------------
!
!re110 CONTINUE
!
!
!*    1.2          SPECIFY ENTRAINMENT RATES FOR DEEP CLOUDS
!                  -----------------------------------------
!
!re120 CONTINUE
  zrg=1._dp/g
  DO 125 jl=1,kproma
     ppbase(jl)=paphp1(jl,kcbot(jl))
     zrrho=(rd*ptenh(jl,kk+1)*(1._dp+vtmpc1*pqenh(jl,kk+1)))           &
                    /paphp1(jl,kk+1)
     zdprho=(paphp1(jl,kk+1)-paphp1(jl,kk))*zrg
     zpmid=0.5_dp*(ppbase(jl)+paphp1(jl,kctop0(jl)))
     zentr=pentr(jl)*pmfu(jl,kk+1)*zdprho*zrrho
     llo1=kk.LT.kcbot(jl).AND.ldcum(jl)
     pdmfde(jl)=MERGE(zentr,0._dp,llo1)
     llo2=llo1.AND.ktype(jl).EQ.2.AND.                                 &
                   (ppbase(jl)-paphp1(jl,kk).LT.0.2e5_dp.OR.           &
                    paphp1(jl,kk).GT.zpmid)
     pdmfen(jl)=MERGE(zentr,0._dp,llo2)
     iklwmin=MAX(klwmin(jl),kctop0(jl)+2)
     llo2=llo1.AND.(ktype(jl).EQ.1 .OR. ktype(jl).EQ.3) .AND.          &
                   (kk.GE.iklwmin.OR.paphp1(jl,kk).GT.zpmid)
     IF(llo2) pdmfen(jl)=zentr
     IF(llo2 .AND. pqenh(jl,kk+1).GT.1.E-5_dp) THEN
        pmfu(jl,kk+1) = MAX(pmfu(jl,kk+1),cmfcmin)
        zentest = MAX(pqte(jl,kk),0._dp)/pqenh(jl,kk+1)
        zentest = MIN(centrmax,zentest/(pmfu(jl,kk+1)*zrrho))
        pdmfen(jl) = zentr+zentest*pmfu(jl,kk+1)*zrrho*zdprho
     ENDIF
125 END DO
!
  RETURN
END SUBROUTINE h2oiso_cuentrt



!======================================================================


SUBROUTINE h2oiso_cubasmc(  kproma, kbdim, klev, kk, klab,                    &
           pten,     pqen,     pqsen,    puen,     pven,               &
!re           ktrac,                                                      &
!re           pxten,    pxtu,     pmfuxt,                                 &
           pverv,    pgeo,     pgeoh,    ldcum,    ktype,              &
           pmfu,     pmfub,    pentr,    kcbot,                        &
           ptu,      pqu,      plu,      puu,      pvu,                &
           pmfus,    pmfuq,    pmful,    pdmfup,   pmfuu,              &
           pcpen,                                                      &
           pmfuv,                                                      &
!---wiso-code
           kwiso,                                                      &
           pwisoqen,                                                   &
           pwisoqu,  pwisolu,                                          &
           pwisomfuq,pwisomful,pwisodmfup                               &
           )
!---wiso-code-end
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
!!$USE mo_kind,          ONLY : dp
!!$USE mo_constants,     ONLY : g
!!$USE mo_cumulus_flux,  ONLY : lmfdudv, entrmid, cmfcmin, cmfcmax
!
IMPLICIT NONE
!
INTEGER, INTENT (IN) :: kbdim, klev, kproma, kk!re, ktrac

!---wiso-code

INTEGER, INTENT (IN) :: kwiso

!---wiso-code-end

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
REAL(dp) :: pcpen(kbdim,klev)
INTEGER  :: ktype(kbdim),            kcbot(kbdim),                     &
            klab(kbdim,klev)
LOGICAL  :: ldcum(kbdim)
!
!reREAL(dp) :: pxten(kbdim,klev,ktrac), pxtu(kbdim,klev,ktrac),           &
!re            pmfuxt(kbdim,klev,ktrac)

!---wiso-code

REAL(dp)::  pwisoqen(kbdim,klev,kwiso)

REAL(dp) :: pwisoqu(kbdim,klev,kwiso),                                 &
            pwisolu(kbdim,klev,kwiso),                                 &
            pwisomfuq(kbdim,klev,kwiso),                               &
            pwisomful(kbdim,klev,kwiso),pwisodmfup(kbdim,klev,kwiso)

!---wiso-code-end

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

!---wiso-code

!DIR$ IVDEP
!OCL NOVREC
  DO jt=1,kwiso
    DO jl=1,kproma
      IF (llo3(jl)) THEN
        pwisoqu(jl,kk+1,jt)=pwisoqen(jl,kk,jt)
        pwisolu(jl,kk+1,jt)=0._dp
        pwisomfuq(jl,kk+1,jt)=pmfub(jl)*pwisoqu(jl,kk+1,jt)
        pwisomful(jl,kk+1,jt)=0._dp
        pwisodmfup(jl,kk+1,jt)=0._dp
      ENDIF
    END DO
  END DO

!---wiso-code-end

!DIR$ IVDEP
!OCL NOVREC
!re  DO 1504 jt=1,ktrac
!re     DO 1502 jl=1,kproma
!re        IF (llo3(jl)) THEN
!re           pxtu(jl,kk+1,jt)=pxten(jl,kk,jt)
!re           pmfuxt(jl,kk+1,jt)=pmfub(jl)*pxtu(jl,kk+1,jt)
!re        ENDIF
!re1502 END DO
!re1504 END DO
!
!
  RETURN
END SUBROUTINE h2oiso_cubasmc




!======================================================================


SUBROUTINE h2oiso_cuadjtq(  kproma, kbdim, klev, kk,             &
           pp,       pt,       pq,       ldflag,   kcall)
!
!          M.TIEDTKE         E.C.M.W.F.     12/89
!          D.SALMOND         CRAY(UK))      12/8/91
!
!          PURPOSE.
!          --------
!          TO PRODUCE T,Q AND L VALUES FOR CLOUD ASCENT
!
!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM SUBROUTINES:
!              *CUBASE*   (T AND Q AT CONDENSTION LEVEL)
!              *CUASC*    (T AND Q AT CLOUD LEVELS)
!              *CUINI*    (ENVIRONMENTAL T AND QS VALUES AT HALF LEVELS)
!          INPUT ARE UNADJUSTED T AND Q VALUES,
!          IT RETURNS ADJUSTED VALUES OF T AND Q
!          NOTE: INPUT PARAMETER KCALL DEFINES CALCULATION AS
!               KCALL=0    ENV. T AND QS IN*CUINI*
!               KCALL=1  CONDENSATION IN UPDRAFTS  (E.G. CUBASE, CUASC)
!               KCALL=2  EVAPORATION IN DOWNDRAFTS (E.G. CUDLFS,CUDDRAF)
!
!          EXTERNALS
!          ---------
!          3 LOOKUP TABLES ( TLUCUA, TLUCUB, TLUCUC )
!          FOR CONDENSATION CALCULATIONS.
!          THE TABLES ARE INITIALISED IN *SETPHYS*.
!

!!$  USE mo_kind,           ONLY: dp
!!$  USE mo_constants,      ONLY: vtmpc1
!!$  USE mo_convect_tables, ONLY: tlucua,   & ! table a
!!$                               tlucub,   & ! table b
!!$                               tlucuc,   & ! table c
!!$                               jptlucu1, jptlucu2, &
!!$                               lookuperror, lookupoverflow

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: kcall, kk, klev, kproma, kbdim

  !  Array arguments with intent(In):
  REAL(dp), INTENT (IN) :: pp(kbdim)
  LOGICAL, INTENT (IN) :: ldflag(kbdim)

  !  Array arguments with intent(InOut):
  REAL(dp), INTENT (INOUT) :: pq(kbdim,klev), pt(kbdim,klev)

  !  Local scalars: 
  REAL(dp):: zcond1, zqst1, zdqsdt, zqsat, zes, zcor, zlcdqsdt
  INTEGER :: isum, jl, it, it1
  LOGICAL :: LO

  !  Local arrays: 
  REAL(dp):: zcond(kbdim)


  !  Executable statements 

  lookupoverflow = .FALSE.

  zcond = 0.0_dp
!
!----------------------------------------------------------------------
!
!     2.           CALCULATE CONDENSATION AND ADJUST T AND Q ACCORDINGLY
!                  -----------------------------------------------------
!
200 CONTINUE
  IF (kcall.EQ.1 ) THEN

     isum=0
!DIR$ IVDEP
!OCL NOVREC
     DO 210 jl=1,kproma
        IF(ldflag(jl)) THEN
           it = NINT(pt(jl,kk)*1000._dp)
           IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
           it = MAX(MIN(it,jptlucu2),jptlucu1)
           zes=tlucua(it)/pp(jl)
           zes=MIN(0.5_dp,zes)
           LO=zes<0.4_dp
           zcor=1._dp/(1._dp-vtmpc1*zes)
           zqsat=zes*zcor
           it1=it+1
           it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
           zqst1=tlucua(it1)/pp(jl)
           zqst1=MIN(0.5_dp,zqst1)
           zqst1=zqst1/(1._dp-vtmpc1*zqst1)
           zdqsdt=(zqst1-zqsat)*1000._dp
           zlcdqsdt=MERGE(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),LO)
           zcond(jl)=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
           zcond(jl)=MAX(zcond(jl),0._dp)
           pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond(jl)
           pq(jl,kk)=pq(jl,kk)-zcond(jl)
           IF(ABS(zcond(jl)) > 0._dp) isum=isum+1
        END IF
210  END DO

     IF (lookupoverflow) RETURN! CALL lookuperror ('h2oiso_cuadjtq (1) ')

     IF(isum.EQ.0) go to 230

!DIR$ IVDEP
!OCL NOVREC
     DO 220 jl=1,kproma
        IF(ldflag(jl)) THEN
           IF(ABS(zcond(jl)) > 0._dp) THEN
           it = NINT(pt(jl,kk)*1000._dp)
           IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
           it = MAX(MIN(it,jptlucu2),jptlucu1)
           zes=tlucua(it)/pp(jl)
           zes=MIN(0.5_dp,zes)
           LO=zes<0.4_dp
           zcor=1._dp/(1._dp-vtmpc1*zes)
           zqsat=zes*zcor
           it1=it+1
           it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
           zqst1=tlucua(it1)/pp(jl)
           zqst1=MIN(0.5_dp,zqst1)
           zqst1=zqst1/(1._dp-vtmpc1*zqst1)
           zdqsdt=(zqst1-zqsat)*1000._dp
           zlcdqsdt=MERGE(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),LO)
           zcond1=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
           pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond1
           pq(jl,kk)=pq(jl,kk)-zcond1
           END IF
        END IF
220  END DO

     IF (lookupoverflow) RETURN! CALL lookuperror ('h2oiso_cuadjtq (2) ')

230  CONTINUE

  END IF

  IF(kcall.EQ.2) THEN

     isum=0
!DIR$ IVDEP
!OCL NOVREC
     DO 310 jl=1,kproma
        IF(ldflag(jl)) THEN
           it = NINT(pt(jl,kk)*1000._dp)
           IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
           it = MAX(MIN(it,jptlucu2),jptlucu1)
           zes=tlucua(it)/pp(jl)
           zes=MIN(0.5_dp,zes)
           LO=zes<0.4_dp
           zcor=1._dp/(1._dp-vtmpc1*zes)
           zqsat=zes*zcor
           it1=it+1
           it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
           zqst1=tlucua(it1)/pp(jl)
           zqst1=MIN(0.5_dp,zqst1)
           zqst1=zqst1/(1._dp-vtmpc1*zqst1)
           zdqsdt=(zqst1-zqsat)*1000._dp
           zlcdqsdt=MERGE(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),LO)
           zcond(jl)=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
           zcond(jl)=MIN(zcond(jl),0._dp)
           pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond(jl)
           pq(jl,kk)=pq(jl,kk)-zcond(jl)
           IF(ABS(zcond(jl)) > 0._dp) isum=isum+1
        END IF
310  END DO

     IF (lookupoverflow) RETURN! CALL lookuperror ('h2oiso_cuadjtq (3) ')

     IF(isum.EQ.0) go to 330

!DIR$ IVDEP
!OCL NOVREC
     DO 320 jl=1,kproma
        IF(ldflag(jl) .AND. ABS(zcond(jl)) > 0._dp) THEN
           it = NINT(pt(jl,kk)*1000._dp)
           IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
           it = MAX(MIN(it,jptlucu2),jptlucu1)
           zes=tlucua(it)/pp(jl)
           zes=MIN(0.5_dp,zes)
           LO=zes<0.4_dp
           zcor=1._dp/(1._dp-vtmpc1*zes)
           zqsat=zes*zcor
           it1=it+1
           it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
           zqst1=tlucua(it1)/pp(jl)
           zqst1=MIN(0.5_dp,zqst1)
           zqst1=zqst1/(1._dp-vtmpc1*zqst1)
           zdqsdt=(zqst1-zqsat)*1000._dp
           zlcdqsdt=MERGE(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),LO)
           zcond1=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
           pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond1
           pq(jl,kk)=pq(jl,kk)-zcond1
        END IF
320  END DO

     IF (lookupoverflow) RETURN! CALL lookuperror ('h2oiso_cuadjtq (4) ')

330  CONTINUE

  END IF

  IF(kcall.EQ.0) THEN

     isum=0
!DIR$ IVDEP
!OCL NOVREC
     DO 410 jl=1,kproma
        it = NINT(pt(jl,kk)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zes=tlucua(it)/pp(jl)
        zes=MIN(0.5_dp,zes)
        LO=zes<0.4_dp
        zcor=1._dp/(1._dp-vtmpc1*zes)
        zqsat=zes*zcor
        it1=it+1
        it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
        zqst1=tlucua(it1)/pp(jl)
        zqst1=MIN(0.5_dp,zqst1)
        zqst1=zqst1/(1._dp-vtmpc1*zqst1)
        zdqsdt=(zqst1-zqsat)*1000._dp
        zlcdqsdt=MERGE(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),LO)
        zcond(jl)=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
        pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond(jl)
        pq(jl,kk)=pq(jl,kk)-zcond(jl)
        IF(ABS(zcond(jl)) > 0._dp) isum=isum+1
410  END DO

     IF (lookupoverflow) RETURN! CALL lookuperror ('h2oiso_cuadjtq (5) ')

     IF(isum.EQ.0) go to 430

!DIR$ IVDEP
!OCL NOVREC
     DO 420 jl=1,kproma
        it = NINT(pt(jl,kk)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zes=tlucua(it)/pp(jl)
        zes=MIN(0.5_dp,zes)
        LO=zes<0.4_dp
        zcor=1._dp/(1._dp-vtmpc1*zes)
        zqsat=zes*zcor
        it1=it+1
        it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
        zqst1=tlucua(it1)/pp(jl)
        zqst1=MIN(0.5_dp,zqst1)
        zqst1=zqst1/(1._dp-vtmpc1*zqst1)
        zdqsdt=(zqst1-zqsat)*1000._dp
        zlcdqsdt=MERGE(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),LO)
        zcond1=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
        pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond1
        pq(jl,kk)=pq(jl,kk)-zcond1
420  END DO

     IF (lookupoverflow) RETURN! CALL lookuperror ('h2oiso_cuadjtq (6) ')

430  CONTINUE

  END IF

  IF(kcall.EQ.4) THEN

!DIR$ IVDEP
!OCL NOVREC
     DO 510 jl=1,kproma
        it = NINT(pt(jl,kk)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zes=tlucua(it)/pp(jl)
        zes=MIN(0.5_dp,zes)
        LO=zes<0.4_dp
        zcor=1._dp/(1._dp-vtmpc1*zes)
        zqsat=zes*zcor
        it1=it+1
        it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
        zqst1=tlucua(it1)/pp(jl)
        zqst1=MIN(0.5_dp,zqst1)
        zqst1=zqst1/(1._dp-vtmpc1*zqst1)
        zdqsdt=(zqst1-zqsat)*1000._dp
        zlcdqsdt=MERGE(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),LO)
        zcond(jl)=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
        pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond(jl)
        pq(jl,kk)=pq(jl,kk)-zcond(jl)
510  END DO

     IF (lookupoverflow) RETURN! CALL lookuperror ('h2oiso_cuadjtq (7) ')

!DIR$ IVDEP
!OCL NOVREC
     DO 520 jl=1,kproma
        it = NINT(pt(jl,kk)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zes=tlucua(it)/pp(jl)
        zes=MIN(0.5_dp,zes)
        LO=zes<0.4_dp
        zcor=1._dp/(1._dp-vtmpc1*zes)
        zqsat=zes*zcor
        it1=it+1
        it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
        zqst1=tlucua(it1)/pp(jl)
        zqst1=MIN(0.5_dp,zqst1)
        zqst1=zqst1/(1._dp-vtmpc1*zqst1)
        zdqsdt=(zqst1-zqsat)*1000._dp
        zlcdqsdt=MERGE(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),LO)
        zcond1=(pq(jl,kk)-zqsat)/(1._dp+zlcdqsdt)
        pt(jl,kk)=pt(jl,kk)+tlucuc(it)*zcond1
        pq(jl,kk)=pq(jl,kk)-zcond1
520  END DO

     IF (lookupoverflow) RETURN! CALL lookuperror ('h2oiso_cuadjtq (8) ')

  END IF

  RETURN
END SUBROUTINE h2oiso_cuadjtq




!======================================================================




SUBROUTINE h2oiso_cuadjwisoq(kproma,   kbdim,      klev,      kwiso,    kk, &
                      ptu,      pqu,        pwisoqu,   pwisolu,      &
                      pqcond,   pwisoqcond, ldflag)
!
!      G.HOFFMANN          MPI MET, HAMBURG      1992 
!      M. WERNER,          MPI BGC, JENA,        2004
!      M. WERNER           AWI, BREMERHAVEN      2009
!
!      PURPOSE
!      -------
!      CALCULATES THE FRACTIONATION FACTOR (LIQUID/VAPOUR SOLID/VAPOUR)
!      AND PUTS THE LIQUID (OR SOLID) ISOTOPE CONDENSATE IN (EFFECTIVE) 
!      EQUILIBRIUM WITH THE SURROUNDING VAPOUR
!
!      INTERFACE
!      ---------
!      THIS SUBROUTINE IS CALLED FROM
!        *CUBASE* (CONDENSATION OF CLOUD LIQUID ISOTOPE WATER)

!      INPUT T, 
!            TWO WATER VALUES (E.G.,OLD VAPOUR AND CONDENSATE)
!            TWO ISOTOPE VALUES (E.G.,OLD ISOTOPE VAPOUR AND LIQUID)
!      OUTPUT ONE OR TWO ISOTOPE VALUES (E.G.,NEW ISOTOPE CONDENSATE)
!
!      EXTERNALS
!      ---------
!      NONE

!!$  USE mo_kind,           ONLY: dp
!!$  USE mo_constants,      ONLY: tmelt
  USE messy_h2oiso,      ONLY: fractcal
  USE messy_main_tools_wiso, ONLY: wiso_frac_liq_ice, cwisomin, cwisosec, cthomi
  
  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT(IN)      :: kproma, kbdim, klev, kwiso, kk
  
  !  Array arguments with intent(In):
  REAL(dp), INTENT(IN)     :: ptu(kbdim,klev), pqu(kbdim,klev), pqcond(kbdim,klev)
  LOGICAL, INTENT(IN)      :: ldflag(kbdim)

  !  Array arguments with intent(InOut):
  REAL(dp), INTENT(INOUT)  :: pwisoqu(kbdim,klev,kwiso), pwisolu(kbdim,klev,kwiso)
  REAL(dp), INTENT(INOUT)  :: pwisoqcond(kbdim,klev,kwiso)
  
  !  Local arrays: 
  REAL(dp) :: zwisofrac(kbdim,kwiso), zwisofracliq(kbdim,kwiso), zwisofracice(kbdim,kwiso)

  !  Local scalars: 
  REAL(dp) :: zqliq,  zqice
  REAL(dp) :: zqvapo, zqvapl,    zqconl, zqconi 
  REAL(dp) :: zwisoqvapo
  REAL(dp) :: zdenom, zevrel,    zquot
  INTEGER  :: jl, jt
  LOGICAL  :: LO

!  Executable statements 

! calculate fractionation factors for liq.-vap. and sol.-vap.
!reold  CALL fractcal(kproma,kbdim,kwiso,zwisofrac,zwisofracice,ptu(:,kk),tmelt)
  CALL wiso_frac_liq_ice(kproma,kbdim,kwiso,ptu(:,kk),zwisofracliq,zwisofracice)

! adjusting isotope vapour and liquid depending on fract.
  DO jt=1,kwiso
    DO jl=1,kproma

! divide into ice and liquid
      IF (ptu(jl,kk).gt.tmelt) THEN
        zqliq=1._dp
        zqice=0._dp
      ELSEIF (ptu(jl,kk).lt.cthomi) THEN
        zqice=1._dp
        zqliq=0._dp
      ELSE
        zqice=(tmelt-ptu(jl,kk))/(tmelt-cthomi)
        zqliq=1._dp-zqice
      ENDIF
     
! fractionation of ice and liquid
      IF (ldflag(jl).and.pqcond(jl,kk).gt.cwisomin) THEN
      
        zqvapo=pqu(jl,kk)+pqcond(jl,kk)
        zqconl=zqliq*pqcond(jl,kk)
        zqvapl=zqvapo-zqconl

!reold        zdenom=zqvapl+zwisofrac(jl,jt)*zqconl
        zdenom=zqvapl+zwisofracliq(jl,jt)*zqconl
        lo=abs(zdenom).lt.cwisomin
        IF (lo) THEN
          pwisoqcond(jl,kk,jt)=0._dp
        ELSE
!reold          pwisoqcond(jl,kk,jt)=zwisofrac(jl,jt)*zqconl*pwisoqu(jl,kk,jt)/zdenom
          pwisoqcond(jl,kk,jt)=zwisofracliq(jl,jt)*zqconl*pwisoqu(jl,kk,jt)/zdenom
        ENDIF
 
        zqconi=zqice*pqcond(jl,kk)
        zwisoqvapo=pwisoqu(jl,kk,jt)-pwisoqcond(jl,kk,jt)
          
        IF (zqvapl.gt.cwisomin) THEN
          zevrel=zqconi/zqvapl
          IF (abs(1.-zevrel).lt.cwisosec) zevrel=1._dp
        ELSE
          zevrel=0._dp
        ENDIF
           
        zquot=1._dp-zevrel

if(zquot.gt.0.0_dp)then!rehc
        zdenom=(1._dp-zquot**zwisofracice(jl,jt))
else!rehc
        zdenom=1._dp
endif!rehc

        pwisoqcond(jl,kk,jt)=pwisoqcond(jl,kk,jt)+zwisoqvapo*zdenom

        pwisoqu(jl,kk,jt)=pwisoqu(jl,kk,jt)-pwisoqcond(jl,kk,jt)
        pwisolu(jl,kk,jt)=pwisolu(jl,kk,jt)+pwisoqcond(jl,kk,jt)

      ELSE
      
        zqvapo=pqu(jl,kk)+pqcond(jl,kk)
        zwisoqvapo=pwisoqu(jl,kk,jt)
        
        zevrel=1._dp
        IF ((ABS(pqcond(jl,kk)).GT.cwisomin).AND.(ABS(zqvapo).GT.cwisomin)) zevrel=pqcond(jl,kk)/zqvapo
        IF (ABS(1._dp-zevrel).LT.cwisosec) zevrel=1._dp

        pwisoqu(jl,kk,jt)=pwisoqu(jl,kk,jt)-zwisoqvapo*zevrel
        
      ENDIF
    
    END DO
  END DO

  RETURN

END SUBROUTINE h2oiso_cuadjwisoq












!======================================================================

!======================================================================
END MODULE MESSY_H2OISO_CONVECT_TIEDTKE
!======================================================================
