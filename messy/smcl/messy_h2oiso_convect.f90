
! **********************************************************************
!
! SUBMODEL CORE LAYER (SMCL) ROUTINES FOR MESSy SUBMODEL H2OISO
!
! THIS SUBMODEL IS TO SIMULATE THE CONVECTION OF WATER ISOTOPES
! 
! Author : Roland Eichinger, DLR-IPA, oneday  oneyear
!
! References:
!
!
! **********************************************************************

! **********************************************************************
MODULE messy_h2oiso_convect
! **********************************************************************



   USE messy_main_constants_mem,    ONLY: DP, tmelt
   USE MESSY_H2OISO_CONVECT_TIEDTKE


   IMPLICIT NONE
   !PRIVATE
   SAVE
 
   PUBLIC :: DP


!   CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'h2oiso_convect'
!   CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '0.1'


   ! SUBROUTINES
!   PUBLIC :: h2oiso_cucall



 CONTAINS

! =======================================================================
! **********************************************************************
! =======================================================================


SUBROUTINE h2oiso_cucall(   kproma, kbdim, klev, klevp1, klevm1, ilab,        &
!re                     ktrac,                                            &
!re                     pxtm1,    pxtte,                                  &
                     ptm1,     pqm1,     pum1,     pvm1,               &
                     pxlm1,    pxim1,                                  &
                     ptte,     pqte,     pvom,     pvol,               &
                     pxlte,    pxite,                                  &
                     pverv,    pxtec,    pqtec,    pqhfla,             &
                     papp1,    paphp1,   pgeo,                         &
                     prsfc,    pssfc,    paprc,    paprs,              &
                     ktype,    ldland,                                 &
!re                     ptopmax,                                          &
!---wiso-code
                     kwiso,                                            &
                     pwisoqm1,                                         &
                     pwisoxlm1,pwisoxim1,                              &
                     pwisoqte,                                         &
                     pwisoxlte,pwisoxite,                              &
                     pwisoxtec,                                        &
                     pwisoqtec,                                        &
                     pwisorsfc,pwisossfc,                              &
                     pwisoaprc, pwisoaprs,                             &
                     time_step_len, lookupoverflow,                    &
                     nn, delta_time )
!---wiso-code-end
                     
!
!
!          *CUCALL* - MASTER ROUTINE - PROVIDES INTERFACE FOR:
!                     *CUMASTR* (CUMULUS PARAMETERIZATION)
!
!           M.TIEDTKE      E.C.M.W.F.     12/1989
!
!**   PURPOSE.
!     --------
!
!          *CUCALL* - INTERFACE FOR *CUMASTR*:
!                     PROVIDES INPUT FOR CUMASTR
!                     RECEIVES UPDATED TENDENCIES, PRECIPITATION.
!
!**   INTERFACE.
!     ----------
!
!          *CUCALL* IS CALLED FROM *PHYSC*
!
!     EXTERNALS.
!     ----------
!
!          CUMASTR, CUMASTRT OR CUMASTRH
!
  USE MESSY_convect,               ONLY: convect_param
  USE messy_main_constants_mem,    ONLY: M_air, M_H2O ! vtmpc1=rv/rd-1
  USE messy_main_tools,            ONLY: tlucua, jptlucu1, jptlucu2


  IMPLICIT NONE

  REAL(DP), PARAMETER :: vtmpc1 = M_air / M_H2O - 1.0_dp

  INTEGER, INTENT (IN) :: nn
  INTEGER, INTENT (IN) :: klev, klevm1, klevp1, kproma, kbdim!re, ktrac

  REAL(DP), INTENT(IN) :: time_step_len
  LOGICAL,  INTENT(INOUT) :: lookupoverflow

  REAL(DP), INTENT(IN) :: delta_time

!---wiso-code

  INTEGER, INTENT (IN) :: kwiso

!---wiso-code-end

  REAL(dp) ::  ptm1(kbdim,klev),         pqm1(kbdim,klev),              &
               pum1(kbdim,klev),         pvm1(kbdim,klev),              &
               ptte(kbdim,klev),         pqte(kbdim,klev),              &
               pvom(kbdim,klev),         pvol(kbdim,klev),              &
               pverv(kbdim,klev),        pgeo(kbdim,klev),              &
               papp1(kbdim,klev),        paphp1(kbdim,klevp1)
  REAL(dp) ::  paprc(kbdim),             paprs(kbdim),                  &
               prsfc(kbdim),             pssfc(kbdim)
  INTEGER  ::  ktype(kbdim)
  REAL(dp) ::  pqhfla(kbdim)
!re  REAL(dp)::  ptopmax(kbdim)
  INTEGER  ::  ilab(kbdim,klev)
  REAL(dp) ::  pxtec(kbdim,klev),        pqtec(kbdim,klev)
  REAL(dp) ::  pxlm1(kbdim,klev),        pxim1(kbdim,klev),             &
               pxlte(kbdim,klev),        pxite(kbdim,klev)

!---wiso-code

  REAL(dp) ::  pwisoqm1(kbdim,klev,kwiso), pwisoqte(kbdim,klev,kwiso)
  REAL(dp) ::  pwisoaprc(kbdim,kwiso),     pwisoaprs(kbdim,kwiso),      &
               pwisorsfc(kbdim,kwiso),     pwisossfc(kbdim,kwiso)
  REAL(dp) ::  pwisoxtec(kbdim,klev,kwiso),pwisoqtec(kbdim,klev,kwiso)
  REAL(dp) ::  pwisoxlm1(kbdim,klev,kwiso),pwisoxim1(kbdim,klev,kwiso), &
               pwisoxlte(kbdim,klev,kwiso),pwisoxite(kbdim,klev,kwiso)
              
!---wiso-code-end

  REAL(dp) ::  ztp1(kbdim,klev),         zqp1(kbdim,klev),              &
               zxp1(kbdim,klev),         ztvp1(kbdim,klev),             &
               zup1(kbdim,klev),         zvp1(kbdim,klev),              &
               ztu(kbdim,klev),          zqu(kbdim,klev),               &
               zlu(kbdim,klev),          zlude(kbdim,klev),             &
               zqude(kbdim,klev),                                       &
               zmfu(kbdim,klev),         zmfd(kbdim,klev),              &
               zqsat(kbdim,klev),        zrain(kbdim)
!re  INTEGER ::  itopec2(kbdim)
  INTEGER  ::  icbot(kbdim),             ictop(kbdim)
!re  REAL(dp)::  zxtp1(kbdim,klev,ktrac),  zxtu(kbdim,klev,ktrac),        &
!re              pxtm1(kbdim,klev,ktrac),  pxtte(kbdim,klev,ktrac)
!re  REAL(dp)::  ztopmax(kbdim)
  LOGICAL  ::  locum(kbdim),             ldland(kbdim)

  !  Local scalars: 
  REAL(dp) :: ztmst, zxlp1, zxip1
!re  INTEGER :: ilevmin
  INTEGER  :: jk, jl, jt, it

!---wiso-code

  REAL(dp) ::  zwisoqp1(kbdim,klev,kwiso),                              &
               zwisoxp1(kbdim,klev,kwiso),                              &
               zwisoqu(kbdim,klev,kwiso),                               &
               zwisolu(kbdim,klev,kwiso),  zwisolude(kbdim,klev,kwiso), &
               zwisoqude(kbdim,klev,kwiso),                             &
               zwisorain(kbdim,kwiso)

  !  Local scalars: 
  REAL(dp) ::  zwisoxlp1, zwisoxip1

!---wiso-code-end

  !  Executable statements 

  lookupoverflow = .FALSE.

!
!-----------------------------------------------------------------------
!*    1.           CALCULATE T,Q AND QS AT MAIN LEVELS
!*                 -----------------------------------
!
!
100 CONTINUE
  ztmst=time_step_len
  DO 120 jk=1,klev
     DO 110 jl=1,kproma
        ztp1(jl,jk)=ptm1(jl,jk)+ptte(jl,jk)*ztmst
        zqp1(jl,jk)=MAX(0._dp,pqm1(jl,jk)+pqte(jl,jk)*ztmst)
        zxlp1=pxlm1(jl,jk)+pxlte(jl,jk)*ztmst
        zxip1=pxim1(jl,jk)+pxite(jl,jk)*ztmst
        zxp1(jl,jk)=MAX(0._dp,zxlp1+zxip1)
        ztvp1(jl,jk)=ztp1(jl,jk)*(1._dp+vtmpc1*zqp1(jl,jk)-zxp1(jl,jk))
        zup1(jl,jk)=pum1(jl,jk)+pvom(jl,jk)*ztmst
        zvp1(jl,jk)=pvm1(jl,jk)+pvol(jl,jk)*ztmst
        it = INT(ztp1(jl,jk)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zqsat(jl,jk)=tlucua(it)/papp1(jl,jk)
        zqsat(jl,jk)=MIN(0.5_dp,zqsat(jl,jk))
        zqsat(jl,jk)=zqsat(jl,jk)/(1._dp-vtmpc1*zqsat(jl,jk))
110  END DO

!re     IF (lookupoverflow) CALL lookuperror ('cucall      ')
     IF (lookupoverflow) THEN 
        do jl=1,kproma
           if ( INT(ztp1(jl,jk)*1000.) <jptlucu1 .OR. &
                INT(ztp1(jl,jk)*1000.) >jptlucu2)     &
                print*, jk, jl,ztp1(jl,jk)*1000., ptm1(jl,jk), &
                ptte(jl,jk)*ztmst
        enddo
        RETURN!CALL FINISH ('convect_convec - lookuperror')
     ENDIF

!re     DO 1104 jt=1,ktrac
!re        DO 1102 jl=1,kproma
!re           zxtp1(jl,jk,jt)=pxtm1(jl,jk,jt)+pxtte(jl,jk,jt)*ztmst
!re1102    END DO
!re1104 END DO

120 END DO

!---wiso-code

  DO jt=1,kwiso
    DO jk=1,klev
      DO jl=1,kproma
        zwisoqp1(jl,jk,jt)=MAX(0._dp,pwisoqm1(jl,jk,jt)+pwisoqte(jl,jk,jt)*ztmst)
        zwisoxlp1=pwisoxlm1(jl,jk,jt)+pwisoxlte(jl,jk,jt)*ztmst
        zwisoxip1=pwisoxim1(jl,jk,jt)+pwisoxite(jl,jk,jt)*ztmst
        zwisoxp1(jl,jk,jt)=MAX(0._dp,zwisoxlp1+zwisoxip1)
      END DO
    END DO
  END DO

!---wiso-code-end

  DO 130 jl=1,kproma
     zrain(jl)=0.0_dp
     locum(jl)=.FALSE.
130 END DO

!---wiso-code

  DO jt=1,kwiso
    DO jl=1,kproma
      zwisorain(jl,jt)=0.0_dp
    END DO
  END DO
!---wiso-code-end
!
!
!-----------------------------------------------------------------------
!
!*    2.     CALL 'CUMASTR'(MASTER-ROUTINE FOR CUMULUS PARAMETERIZATION)
!*           -----------------------------------------------------------
!
!
200 CONTINUE

  SELECT CASE (convect_param)
  CASE(1)         ! Tiedtke - Nordeng

     CALL h2oiso_cumastr(kproma, kbdim, klev, klevp1, klevm1, ilab,           &
                  ztp1,     zqp1,     zxp1,     zup1,   zvp1,          &
                  ztvp1,    ldland, &!re,    ktrac,                          &
!re                  zxtp1,    zxtu,     pxtte,                           &
                  pverv,    zqsat,    pqhfla,                          &
                  paphp1,   pgeo,                                      &
                  ptte,     pqte,     pvom,     pvol,                  &
                  prsfc,    pssfc,    paprc,    paprs,  pxtec,         &
                  pqtec,    zqude,                                     &
                  locum,    ktype,    icbot,    ictop,                 &
                  ztu,      zqu,      zlu,      zlude,                 &
                  zmfu,     zmfd,     zrain,                           &
!---wiso-code
                  kwiso,                                               &
                  zwisoqp1, zwisoxp1,                                  &
                  pwisoqte,                                            &
                  pwisorsfc,pwisossfc,pwisoaprc,pwisoaprs,pwisoxtec,   &
                  pwisoqtec,zwisoqude,                                 &
                  zwisoqu,  zwisolu,  zwisolude,                       &
                  zwisorain, &
                  nn, time_step_len,  lookupoverflow, delta_time)
!---wiso-code-end
  CASE(2)         ! Tiedtke
     CALL h2oiso_cumastrt(kproma, kbdim, klev, klevp1, klevm1, ilab,          &
                  ztp1,     zqp1,     zxp1,     zup1,   zvp1,          &
                  ztvp1,    ldland, &!re,    ktrac,                          &
!re                  zxtp1,    zxtu,     pxtte,                           &
                  pverv,    zqsat,    pqhfla,                          &
                  paphp1,   pgeo,                                      &
                  ptte,     pqte,     pvom,     pvol,                  &
                  prsfc,    pssfc,    paprc,    paprs,  pxtec,         &
                  pqtec,    zqude,                                     &
                  locum,    ktype,    icbot,    ictop,                 &
                  ztu,      zqu,      zlu,      zlude,                 &
                  zmfu,     zmfd,     zrain,                           &
!---wiso-code
                  kwiso,                                               &
                  zwisoqp1, zwisoxp1,                                  &
                  pwisoqte,                                            &
                  pwisorsfc,pwisossfc,pwisoaprc,pwisoaprs,pwisoxtec,   &
                  pwisoqtec,zwisoqude,                                 &
                  zwisoqu,  zwisolu,  zwisolude,                       &
                  zwisorain, &
                  nn, time_step_len,  lookupoverflow, delta_time)
!---wiso-code-end
  CASE(3)         ! Tiedtke - Hybrid
     CALL h2oiso_cumastrh(kproma, kbdim, klev, klevp1, klevm1, ilab,          &
                  ztp1,     zqp1,     zxp1,     zup1,   zvp1,          &
                  ztvp1,    ldland, &!re,    ktrac,                          &
!re                  zxtp1,    zxtu,     pxtte,                           &
                  pverv,    zqsat,    pqhfla,                          &
                  paphp1,   pgeo,                                      &
                  ptte,     pqte,     pvom,     pvol,                  &
                  prsfc,    pssfc,    paprc,    paprs,  pxtec,         &
                  pqtec,    zqude,                                     &
                  locum,    ktype,    icbot,    ictop,                 &
                  ztu,      zqu,      zlu,      zlude,                 &
                  zmfu,     zmfd,     zrain,                           &
!---wiso-code
                  kwiso,                                               &
                  zwisoqp1, zwisoxp1,                                  &
                  pwisoqte,                                            &
                  pwisorsfc,pwisossfc,pwisoaprc,pwisoaprs,pwisoxtec,   &
                  pwisoqtec,zwisoqude,                                 &
                  zwisoqu,  zwisolu,  zwisolude,                       &
                  zwisorain, &
                  nn, time_step_len,  lookupoverflow, delta_time)
!---wiso-code-end

  END SELECT
    IF (lookupoverflow) RETURN
!
!
! ------------------------------------------------------------------
!
!*     3.     PRESSURE ALTITUDE OF CONVECTIVE CLOUD TOPS.
!             -------- -------- -- ---------- ----- -----
!
300 CONTINUE
!
!re  ilevmin=klev-4
!
!re  DO 301 jl=1,kproma
!re     itopec2(jl)=klevp1
!re301 END DO
!
!re  DO 303 jk=1,ilevmin
!re     DO 302 jl=1,kproma
!re        IF(ilab(jl,jk).EQ.2 .AND. itopec2(jl).EQ.klevp1) THEN
!re           itopec2(jl)=jk
!re        END IF
!re302  END DO
!re303 END DO
!
!re   ztopmax(1:kproma) = ptopmax(1:kproma)

!re   DO 304 jl=1,kproma
!re      IF(itopec2(jl).EQ.1) THEN
!re         ptopmax(jl)=papp1(jl,1)
!re      ELSE IF(itopec2(jl).NE.klevp1) THEN
!re         ptopmax(jl)=paphp1(jl,itopec2(jl))
!re      ELSE
!re         ptopmax(jl)=99999._dp
!re      END IF
!re      ptopmax(jl)=MIN(ptopmax(jl),ztopmax(jl))
!re 304 END DO
!
!
!---------------------------------------------------------------------
!
  RETURN
END SUBROUTINE h2oiso_cucall


!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
!


 SUBROUTINE h2oiso_cumastr(  kproma, kbdim, klev, klevp1, klevm1, ilab,       &
           pten,     pqen,     pxen,     puen,     pven,               &
           ptven,    ldland, &!re,    ktrac,                                 &
!re           pxten,    pxtu,     pxtte,                                  &
           pverv,    pqsen,    pqhfla,                                 &
           paphp1,   pgeo,                                             &
           ptte,     pqte,     pvom,     pvol,                         &
           prsfc,    pssfc,    paprc,    paprs,    pxtec,              &
           pqtec,    pqude,                                            &
           ldcum,    ktype,    kcbot,    kctop,                        &
           ptu,      pqu,      plu,      plude,                        &
           pmfu,     pmfd,     prain,                                  &
!---wiso-code
           kwiso,                                                      &
           pwisoqen, pwisoxen,                                         &
           pwisoqte,                                                   &
           pwisorsfc,pwisossfc,pwisoaprc,pwisoaprs,pwisoxtec,          &
           pwisoqtec,pwisoqude,                                        &
           pwisoqu,  pwisolu,  pwisolude,                              &
           pwisorain, &
           nn, time_step_len, lookupoverflow, delta_time)
!---wiso-code-end
!
!**** *CUMASTR*  MASTER ROUTINE FOR CUMULUS MASSFLUX-SCHEME
!
!     M.TIEDTKE      E.C.M.W.F.     1986/1987/1989
!
!     PURPOSE
!     -------
!
!          THIS ROUTINE COMPUTES THE PHYSICAL TENDENCIES OF THE
!     PROGNOSTIC VARIABLES T,Q,U AND V DUE TO CONVECTIVE PROCESSES.
!     PROCESSES CONSIDERED ARE: CONVECTIVE FLUXES, FORMATION OF
!     PRECIPITATION, EVAPORATION OF FALLING RAIN BELOW CLOUD BASE,
!     SATURATED CUMULUS DOWNDRAFTS.
!
!**   INTERFACE.
!     ----------
!
!          *CUMASTR* IS CALLED FROM *CUCALL*
!     THE ROUTINE TAKES ITS INPUT FROM THE LONG-TERM STORAGE
!     T,Q,U,V,PHI AND P AND MOISTURE TENDENCIES.
!     IT RETURNS ITS OUTPUT TO THE SAME SPACE
!      1.MODIFIED TENDENCIES OF MODEL VARIABLES
!      2.RATES OF CONVECTIVE PRECIPITATION
!        (USED IN SUBROUTINE SURF)
!
!     METHOD
!     -------
!
!     PARAMETERIZATION IS DONE USING A MASSFLUX-SCHEME.
!        (1) DEFINE CONSTANTS AND PARAMETERS
!        (2) SPECIFY VALUES (T,Q,QS...) AT HALF LEVELS AND
!            INITIALIZE UPDRAFT- AND DOWNDRAFT-VALUES IN 'CUINI'
!        (3) CALCULATE CLOUD BASE IN 'CUBASE'
!            AND SPECIFY CLOUD BASE MASSFLUX FROM PBL MOISTURE BUDGET
!        (4) DO CLOUD ASCENT IN 'CUASC' IN ABSENCE OF DOWNDRAFTS
!        (5) DO DOWNDRAFT CALCULATIONS:
!              (A) DETERMINE VALUES AT LFS IN 'CUDLFS'
!              (B) DETERMINE MOIST DESCENT IN 'CUDDRAF'
!              (C) RECALCULATE CLOUD BASE MASSFLUX CONSIDERING THE
!                  EFFECT OF CU-DOWNDRAFTS
!        (6) DO FINAL CLOUD ASCENT IN 'CUASC'
!        (7) DO FINAL ADJUSMENTS TO CONVECTIVE FLUXES IN 'CUFLX',
!            DO EVAPORATION IN SUBCLOUD LAYER
!        (8) CALCULATE INCREMENTS OF T AND Q IN 'CUDTDQ'
!        (9) CALCULATE INCREMENTS OF U AND V IN 'CUDUDV'
!
!     EXTERNALS.
!     ----------
!
!       CUINI:  INITIALIZES VALUES AT VERTICAL GRID USED IN CU-PARAMETR.
!       CUBASE: CLOUD BASE CALCULATION FOR PENETR.AND SHALLOW CONVECTION
!       CUASC:  CLOUD ASCENT FOR ENTRAINING PLUME
!       CUDLFS: DETERMINES VALUES AT LFS FOR DOWNDRAFTS
!       CUDDRAF:DOES MOIST DESCENT FOR CUMULUS DOWNDRAFTS
!       CUFLX:  FINAL ADJUSTMENTS TO CONVECTIVE FLUXES (ALSO IN PBL)
!       CUDQDT: UPDATES TENDENCIES FOR T AND Q
!       CUDUDV: UPDATES TENDENCIES FOR U AND V
!
!     SWITCHES.
!     --------
!
!          LMFPEN=.T.   PENETRATIVE CONVECTION IS SWITCHED ON
!          LMFSCV=.T.   SHALLOW CONVECTION IS SWITCHED ON
!          LMFMID=.T.   MIDLEVEL CONVECTION IS SWITCHED ON
!          LMFDD=.T.    CUMULUS DOWNDRAFTS SWITCHED ON
!          LMFDUDV=.T.  CUMULUS FRICTION SWITCHED ON
!
!     MODEL PARAMETERS (DEFINED IN SUBROUTINE CUPARAM)
!     ------------------------------------------------
!     ENTRPEN    ENTRAINMENT RATE FOR PENETRATIVE CONVECTION
!     ENTRSCV    ENTRAINMENT RATE FOR SHALLOW CONVECTION
!     ENTRMID    ENTRAINMENT RATE FOR MIDLEVEL CONVECTION
!     ENTRDD     ENTRAINMENT RATE FOR CUMULUS DOWNDRAFTS
!     CMFCTOP    RELATIVE CLOUD MASSFLUX AT LEVEL ABOVE NONBUOYANCY LEVEL
!     CMFCMAX    MAXIMUM MASSFLUX VALUE ALLOWED FOR
!     CMFCMIN    MINIMUM MASSFLUX VALUE (FOR SAFETY)
!     CMFDEPS    FRACTIONAL MASSFLUX FOR DOWNDRAFTS AT LFS
!     CPRCON     COEFFICIENT FOR CONVERSION FROM CLOUD WATER TO RAIN
!
!     REFERENCE.
!     ----------
!
!          PAPER ON MASSFLUX SCHEME (TIEDTKE,1989)
!

   USE messy_main_constants_mem,     ONLY: g, alv, als, tmelt, vtmpc1, rd
   USE MESSY_CONVECT_TIEDTKE_PARAM,  ONLY: entrpen, entrscv, lmfdd, cmfdeps, lmfdudv
   USE messy_main_tools,             ONLY: tlucua, tlucub, jptlucu1, jptlucu2

   IMPLICIT NONE

   !INTEGER, INTENT(IN) :: nn
   REAL(DP), INTENT(IN) :: time_step_len
   LOGICAL,  INTENT(INOUT) :: lookupoverflow

   INTEGER,  INTENT(IN) :: kbdim, klev, kproma, klevp1, klevm1, nn !re, ktrac

   REAL(dp), INTENT(IN) :: delta_time

   !---wiso-code

   INTEGER,  INTENT(IN) :: kwiso

   !---wiso-code-end

   REAL(dp) :: pten(kbdim,klev),        pqen(kbdim,klev),                  &
               pxen(kbdim,klev),        ptven(kbdim,klev),                 &
               puen(kbdim,klev),        pven(kbdim,klev),                  &
               ptte(kbdim,klev),        pqte(kbdim,klev),                  &
               pvom(kbdim,klev),        pvol(kbdim,klev),                  &
               pqsen(kbdim,klev),       pgeo(kbdim,klev),                  &
               paphp1(kbdim,klevp1),                                       &
               pverv(kbdim,klev),       pqude(kbdim,klev)                                            
   REAL(dp) :: ptu(kbdim,klev),         pqu(kbdim,klev),                   &
               plu(kbdim,klev),         plude(kbdim,klev),                 &
               pmfu(kbdim,klev),        pmfd(kbdim,klev),                  &
               paprc(kbdim),            paprs(kbdim),                      &
               prsfc(kbdim),            pssfc(kbdim),                      &
               prain(kbdim),            pqhfla(kbdim)
   ! convective gust?
   INTEGER  :: kcbot(kbdim),            kctop(kbdim),                      &
               ktype(kbdim)
   REAL(dp) :: pxtec(kbdim,klev),       pqtec(kbdim,klev)

   !---wiso-code
        
   REAL(dp) :: pwisoqen(kbdim,klev,kwiso),                                 &
               pwisoxen(kbdim,klev,kwiso),                                 &
               pwisoqte(kbdim,klev,kwiso)
   REAL(dp) :: pwisoqu(kbdim,klev,kwiso),                                  &
               pwisolu(kbdim,klev,kwiso), pwisolude(kbdim,klev,kwiso),     &
               pwisoaprc(kbdim,kwiso),    pwisoaprs(kbdim,kwiso),          &
               pwisorsfc(kbdim,kwiso),    pwisossfc(kbdim,kwiso),          &
               pwisorain(kbdim,kwiso)
   REAL(dp) :: pwisoxtec(kbdim,klev,kwiso),pwisoqtec(kbdim,klev,kwiso),    &
               pwisoqude(kbdim,klev,kwiso)

   !---wiso-code-end

   REAL(dp) :: ztenh(kbdim,klev),       zqenh(kbdim,klev),                 &
               zxenh(kbdim,klev),                                          &
               zgeoh(kbdim,klev),       zqsenh(kbdim,klev),                &
               ztd(kbdim,klev),         zqd(kbdim,klev),                   &
               zmfus(kbdim,klev),       zmfds(kbdim,klev),                 &
               zmfuq(kbdim,klev),       zmfdq(kbdim,klev),                 &
               zdmfup(kbdim,klev),      zdmfdp(kbdim,klev),                &
               zmful(kbdim,klev),       zrfl(kbdim),                       &
               zuu(kbdim,klev),         zvu(kbdim,klev),                   &
               zud(kbdim,klev),         zvd(kbdim,klev)
   REAL(dp) :: zcpen(kbdim,klev),       zcpcu(kbdim,klev)
   REAL(dp) :: zentr(kbdim),            zhcbase(kbdim),                    &
               zmfub(kbdim),            zmfub1(kbdim),                     &
               zdqpbl(kbdim),           zdqcv(kbdim)
   REAL(dp) :: zsfl(kbdim),             zdpmel(kbdim,klev)
   REAL(dp) :: zcape(kbdim),            zheat(kbdim)
   REAL(dp) :: zhmin(kbdim)
   REAL(dp) :: zhhatt(kbdim,klev)
   INTEGER  :: ihmin(kbdim)
   INTEGER  :: ilab(kbdim,klev),        idtop(kbdim),                      &
               ictop0(kbdim),           ilwmin(kbdim)
   !re REAL(dp) :: pxten(kbdim,klev,ktrac), pxtte(kbdim,klev,ktrac),           &
   !re             pxtu(kbdim,klev,ktrac),                                     &
   !re             zxtenh(kbdim,klev,ktrac),zxtd(kbdim,klev,ktrac),            &
   !re             zmfuxt(kbdim,klev,ktrac),zmfdxt(kbdim,klev,ktrac)
   LOGICAL  :: loddraf(kbdim),          ldland(kbdim)
   LOGICAL  :: ldcum(kbdim)

   LOGICAL  :: llo1, lo

! local variables

   REAL(dp) :: zcons2, ztau,  zqumqe, zdqmin, zmfmax, zalvs,  zalvdcp, zqalv, &
               zhsat,  zqsat, zes,    zcor,   zqst1,  zdqsdt, zgam,    zzz,   &
               zhhat,  zb,    zbi,    zro,    zdz,    zdhdz,  zdepth,  zfac,  &
               zrh, zeps, zpbmpt
   INTEGER  :: jl, jk, ikb, it, it1, jt, itopm2

   !---wiso-code

   REAL(dp) :: zwisoqenh(kbdim,klev,kwiso),                                &
               zwisoxenh(kbdim,klev,kwiso),                                &
               zwisoqd(kbdim,klev,kwiso),                                  &
               zwisomfuq(kbdim,klev,kwiso),zwisomfdq(kbdim,klev,kwiso),    &
               zwisodmfup(kbdim,klev,kwiso),zwisodmfdp(kbdim,klev,kwiso),  &
               zwisomful(kbdim,klev,kwiso)

   REAL(dp) :: zmelt_tmp(kbdim,klev),        &  ! temporary value
               zprec_frac(kbdim,klev),       &  ! temporary value
               zrain_tmp(kbdim,klev),        &  ! temporary value 
               zrfl_tmp(kbdim),              &  ! temporary value
               zsfl_tmp(kbdim),              &  ! temporary value
               zqsto(kbdim,klev),            &  ! old moisture value
               ztsto(kbdim,klev),            &  ! old temperature value
               zwisoqsto(kbdim,klev,kwiso),  &  ! old water isotope moisture value
               zwisorfl(kbdim,kwiso),        &  ! temporary value
               zwisosfl(kbdim,kwiso)            ! temporary value
           
   !---wiso-code-end
   !
   !  INTRINSIC FUNCTIONS
   INTRINSIC MIN, MAX
!
!  Executable statements

  lookupoverflow = .FALSE.
!-----------------------------------------------------------------------
!
!     1.           SPECIFY CONSTANTS AND PARAMETERS
!                  --------------------------------
!
100 CONTINUE
!
  zcons2=1._dp/(g*time_step_len)
  ztau=MIN(3._dp*3600._dp,7200._dp*63._dp/REAL(nn,dp))  ! op_ff_20161020 REAL(nn,dp)

!---wiso-code

  zprec_frac(:,:)=0._dp
  zrain_tmp(:,:)=0._dp

!---wiso-code-end
!
!----------------------------------------------------------------------
!
!*    2.           INITIALIZE VALUES AT VERTICAL GRID POINTS IN 'CUINI'
!                  ---------------------------------------------------
!
200 CONTINUE
  CALL h2oiso_cuini(kproma, kbdim, klev, klevp1, klevm1,                      &
             pten,     pqen,     pqsen,    pxen,     puen,     pven,   &
             ptven, &!re ,    ktrac,                                          &
!re             pxten,    zxtenh,   pxtu,     zxtd,     zmfuxt,   zmfdxt, &
             pverv,    pgeo,     paphp1,   zgeoh,                      &
             ztenh,    zqenh,    zqsenh,   zxenh,    ilwmin,           &
             ptu,      pqu,      ztd,      zqd,                        &
             zuu,      zvu,      zud,      zvd,                        &
             pmfu,     pmfd,     zmfus,    zmfds,                      &
             zmfuq,    zmfdq,    zdmfup,   zdmfdp,                     &
             zcpen,    zcpcu,                                          &
             zdpmel,   plu,      plude,    pqude,    ilab,             &
!---wiso-code
             kwiso,                                                    &
             pwisoqen, pwisoxen,                                       &
             zwisoqenh,zwisoxenh,                                      &
             pwisoqu,  zwisoqd,                                        &
             zwisomfuq,zwisomfdq,zwisodmfup,zwisodmfdp,                &
             pwisolu,  pwisolude,pwisoqude                             &
             )         
!---wiso-code-end
    IF (lookupoverflow) RETURN!re
!
!-----------------------------------------------------------------------
!
!*    3.0          CLOUD BASE CALCULATIONS
!                  -----------------------
!
300 CONTINUE
!
!*             (A) DETERMINE CLOUD BASE VALUES IN 'CUBASE'
!                  ---------------------------------------
!
  CALL h2oiso_cubase(kproma, kbdim, klev, klevp1, klevm1,                     &
              ztenh,    zqenh,    zgeoh,    paphp1,                    &
              ptu,      pqu,      plu,                                 &
              puen,     pven,     zuu,      zvu,                       &
              zcpcu,                                                   &
              ldcum,    kcbot,    ilab,                                &
!---wiso-code
              kwiso,                                                   &
              pwisoqu,  pwisolu                                        &
              )
!---wiso-code-end
    IF (lookupoverflow) RETURN!re
!

!*             (B) DETERMINE TOTAL MOISTURE CONVERGENCE AND
!*                 THEN DECIDE ON TYPE OF CUMULUS CONVECTION
!                  -----------------------------------------
!
  jk=1
  DO 310 jl=1,kproma
     zdqpbl(jl)=0.0_dp
     zdqcv(jl)=pqte(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))
     idtop(jl)=0
310 END DO
  DO 320 jk=2,klev
     DO 315 jl=1,kproma
        zdqcv(jl)=zdqcv(jl)+pqte(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))
        IF(jk.GE.kcbot(jl)) zdqpbl(jl)=zdqpbl(jl)+pqte(jl,jk)          &
                                       *(paphp1(jl,jk+1)-paphp1(jl,jk))
315  END DO
320 END DO
!
!*             (C) DETERMINE MOISTURE SUPPLY FOR BOUNDARY LAYER
!*                 AND DETERMINE CLOUD BASE MASSFLUX IGNORING
!*                 THE EFFECTS OF DOWNDRAFTS AT THIS STAGE
!                  ------------------------------------------
!
!DIR$ IVDEP
!DIR$ CONCURRENT
  DO 340 jl=1,kproma
     ikb=kcbot(jl)
     zqumqe=pqu(jl,ikb)+plu(jl,ikb)-zqenh(jl,ikb)
     zdqmin=MAX(0.01_dp*zqenh(jl,ikb),1.e-10_dp)
     llo1=zdqpbl(jl).GT.0._dp.AND.zqumqe.GT.zdqmin.AND.ldcum(jl)
     zmfub(jl)=MERGE(zdqpbl(jl)/(g*MAX(zqumqe,zdqmin)),0.01_dp,llo1)
     zmfmax=(paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
     zmfub(jl)=MIN(zmfub(jl),zmfmax)
     IF(.NOT.llo1) ldcum(jl)=.FALSE.
     ktype(jl)=MERGE(1,2,zdqcv(jl).GT.MAX(0._dp,-1.1_dp*pqhfla(jl)*g))
     zentr(jl)=MERGE(entrpen,entrscv,ktype(jl).EQ.1)
340 END DO
!
!-----------------------------------------------------------------------
!*    4.0          DETERMINE CLOUD ASCENT FOR ENTRAINING PLUME
!                  -------------------------------------------
!
400 CONTINUE
!
!*             (A) ESTIMATE CLOUD HEIGHT FOR ENTRAINMENT/DETRAINMENT
!*                 CALCULATIONS IN CUASC (MAX.POSSIBLE CLOUD HEIGHT
!*                 FOR NON-ENTRAINING PLUME, FOLLOWING A.-S.,1974)
!                  -------------------------------------------------
!
!DIR$ IVDEP
  DO 410 jl=1,kproma
     ikb=kcbot(jl)
     zalvs=MERGE(alv,als,ptu(jl,ikb)>tmelt)
     zhcbase(jl)=zcpcu(jl,ikb)*ptu(jl,ikb)+zgeoh(jl,ikb)+              &
                                           zalvs*pqu(jl,ikb)
     ictop0(jl)=kcbot(jl)-1
410 END DO
  DO 430 jk=klevm1,3,-1
!DIR$ IVDEP
     DO 420 jl=1,kproma
        zalvs=MERGE(alv,als,ztenh(jl,jk)>tmelt)
        zalvdcp=zalvs/zcpcu(jl,jk)
        zqalv=1._dp/zalvs
        zhsat=zcpcu(jl,jk)*ztenh(jl,jk)+zgeoh(jl,jk)+zalvs*zqsenh(jl,jk)
        it = NINT(ztenh(jl,jk)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zes=tlucua(it)/paphp1(jl,jk)
        zes=MIN(0.5_dp,zes)
        LO=zes<0.4_dp
        zcor=1._dp/(1._dp-vtmpc1*zes)
        zqsat=zes*zcor
        it1=it+1
        it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
        zqst1=tlucua(it1)/paphp1(jl,jk)
        zqst1=MIN(0.5_dp,zqst1)
        zqst1=zqst1/(1._dp-vtmpc1*zqst1)
        zdqsdt=(zqst1-zqsat)*1000._dp
        zgam=MERGE(zalvdcp*zdqsdt,zqsat*zcor*tlucub(it),LO)
        zzz=zcpcu(jl,jk)*ztenh(jl,jk)*vtmpc1
        zhhat=zhsat-(zzz+zgam*zzz)/(1._dp+zgam*zzz*zqalv)*             &
                          MAX(zqsenh(jl,jk)-zqenh(jl,jk),0._dp)
        zhhatt(jl,jk)=zhhat
        IF(jk.LT.ictop0(jl).AND.zhcbase(jl).GT.zhhat) ictop0(jl)=jk
420  END DO
430 END DO
!!
!!     DEEP CONVECTION IF CLOUD DEPTH > 200 HPA, ELSE SHALLOW
!!     (CLOUD DEPTH FROM NON-ENTRAINIG PLUME)
!!
!  DO jl=1,kproma
!     ktype(jl)=MERGE(1,2,                                              &
!                paphp1(jl,kcbot(jl))-paphp1(jl,ictop0(jl)).gt.2.E4_dp)
!     zentr(jl)=MERGE(entrpen,entrscv,ktype(jl).eq.1)
!  ENDDO
!!
     IF (lookupoverflow) RETURN !re CALL lookuperror ('cumastr ')
!
!                  FIND LOWEST POSSIBLE ORG. DETRAINMENT LEVEL
!                  -------------------------------------------
!
  DO jl=1,kproma
     zhmin(jl)=0._dp
     ihmin(jl)=0
     llo1=ldcum(jl).AND.ktype(jl).EQ.1
     IF(llo1) THEN
        ikb=kcbot(jl)
        ihmin(jl)=ikb
     ENDIF
  ENDDO
!
  zb=25._dp
  zbi=1._dp/(zb*g) 
  DO jk=klev,1,-1
     DO jl=1,kproma
        llo1=ldcum(jl).AND.ktype(jl).EQ.1.AND.ihmin(jl).EQ.kcbot(jl)
        IF(llo1.AND.jk.LT.kcbot(jl).AND.jk.GE.ictop0(jl)) THEN
           zalvs=MERGE(alv,als,ztenh(jl,jk)>tmelt)
           ikb=kcbot(jl)
           zro=paphp1(jl,jk)/(rd*ztenh(jl,jk)*                         &
                                           (1._dp+vtmpc1*zqenh(jl,jk)))
           zdz=(paphp1(jl,jk)-paphp1(jl,jk-1))/(g*zro)
           zdhdz=( zcpen(jl,jk-1)*pten(jl,jk-1)-                        &
                        zcpen(jl,jk)*pten(jl,jk)+                       &
                          zalvs*(pqen(jl,jk-1)-pqen(jl,jk))+            &
                                (pgeo(jl,jk-1)-pgeo(jl,jk)) )*g/        &
                                (pgeo(jl,jk-1)-pgeo(jl,jk))
           zdepth=zgeoh(jl,jk)-zgeoh(jl,ikb)
           zfac=SQRT(1._dp+zdepth*zbi)
           zhmin(jl)=zhmin(jl) + zdhdz*zfac*zdz
           zrh=-zalvs*(zqsenh(jl,jk)-zqenh(jl,jk))*zfac
           IF(zhmin(jl).GT.zrh) ihmin(jl)=jk
        ENDIF
     ENDDO
  ENDDO
!
  DO jl=1,kproma
     IF(ldcum(jl).AND.ktype(jl).EQ.1) THEN
        IF(ihmin(jl).LT.ictop0(jl)) ihmin(jl)=ictop0(jl)
     ENDIF
  ENDDO
!
!*             (B) DO ASCENT IN 'CUASC'IN ABSENCE OF DOWNDRAFTS
!                  --------------------------------------------
!
  CALL h2oiso_cuasc(kproma, kbdim, klev, klevp1, klevm1,                      &
             ztenh,    zqenh,    puen,     pven,                       &
!re             ktrac,                                                    &
!re             zxtenh,   pxten,    pxtu,     zmfuxt,                     &
             pten,     pqen,     pqsen,                                &
             pgeo,     zgeoh,    paphp1,                               &
             pqte,     pverv,    ilwmin,                               &
             ldcum,    ldland,   ktype,    ilab,                       &
             ptu,      pqu,      plu,      zuu,      zvu,              &
             pmfu,     zmfub,    zentr,                                &
             zmfus,    zmfuq,                                          &
             zmful,    plude,    pqude,    zdmfup,                     &
             ihmin,    zhhatt,   zhcbase,  zqsenh,                     &
             zcpen,    zcpcu,                                          &
             kcbot,    kctop,    ictop0,                               &
!---wiso-code
             kwiso,                                                    &
             zwisoqenh,                                                &
             pwisoqen,                                                 &
             pwisoqu,  pwisolu,                                        &
             zwisomfuq,                                                &
             zwisomful,pwisolude,pwisoqude,zwisodmfup                  &
             , nn, time_step_len)
!---wiso-code-end
!
!*         (C) CHECK CLOUD DEPTH AND CHANGE ENTRAINMENT RATE ACCORDINGLY
!              CALCULATE PRECIPITATION RATE (FOR DOWNDRAFT CALCULATION)
!              ---------------------------------------------------------
!
!DIR$ IVDEP
  DO 440 jl=1,kproma
     zpbmpt=paphp1(jl,kcbot(jl))-paphp1(jl,kctop(jl))
     IF(ldcum(jl).AND.ktype(jl).EQ.1.AND.zpbmpt.LT.2.e4_dp) ktype(jl)=2
     IF(ldcum(jl)) ictop0(jl)=kctop(jl)
     IF(ktype(jl).EQ.2) zentr(jl)=entrscv
     zrfl(jl)=zdmfup(jl,1)
440 END DO
  DO 460 jk=2,klev
     DO 450 jl=1,kproma
        zrfl(jl)=zrfl(jl)+zdmfup(jl,jk)
450  END DO
460 END DO
!
!-----------------------------------------------------------------------
!*    5.0          CUMULUS DOWNDRAFT CALCULATIONS
!                  ------------------------------
!
500 CONTINUE
!
  IF(lmfdd) THEN
!
!*             (A) DETERMINE LFS IN 'CUDLFS'
!                  -------------------------
!
     CALL h2oiso_cudlfs(kproma,   kbdim,    klev,     klevp1,                 &
                 ztenh,    zqenh,    puen,     pven,                   &
!re                 ktrac,                                                &
!re                 zxtenh,   pxtu,     zxtd,     zmfdxt,                 &
                 zgeoh,    paphp1,                                     &
                 ptu,      pqu,      zuu,      zvu,                    &
                 ldcum,    kcbot,    kctop,    zmfub,    zrfl,         &
                 ztd,      zqd,      zud,      zvd,                    &
                 pmfd,     zmfds,    zmfdq,    zdmfdp,                 &
                 zcpcu,                                                &
                 idtop,    loddraf,                                    &
!---wiso-code
                 kwiso,                                                &
                 zwisoqenh,                                            &
                 pwisoqu,                                              &
                 zwisoqd,                                              &
                 zwisomfdq,zwisodmfdp,                                 &
                 zdmfup,   zwisodmfup                                  &
                 )
!---wiso-code-end
!
!*            (B)  DETERMINE DOWNDRAFT T,Q AND FLUXES IN 'CUDDRAF'
!                  -----------------------------------------------
!
     CALL h2oiso_cuddraf(kproma,   kbdim,    klev,     klevp1,                &
                  ztenh,    zqenh,    puen,     pven,                  &
!re                  ktrac,                                               &
!re                  zxtenh,   zxtd,     zmfdxt,                          &
                  zgeoh,    paphp1,   zrfl,                            &
                  ztd,      zqd,      zud,      zvd,                   &
                  pmfd,     zmfds,    zmfdq,    zdmfdp,                &
                  zcpcu,                                               &
                  loddraf,                                             &
!---wiso-code
                  kwiso,                                               &
                  zwisoqenh,                                           &
                  zwisoqd,                                             &
                  zwisomfdq,zwisodmfdp,                                &
                  zdmfup,   zwisodmfup                                 &
                  )
!---wiso-code-end
!
  END IF
!
!*    5.5          RECALCULATE CLOUD BASE MASSFLUX FROM A
!*                 CAPE CLOSURE FOR DEEP CONVECTION (KTYPE=1)
!*                 AND BY PBL EQUILIBRUM TAKING DOWNDRAFTS INTO
!*                 ACCOUNT FOR SHALLOW CONVECTION (KTYPE=2)
!                  -------------------------------------------
!
  DO jl=1,kproma
     zheat(jl)=0._dp
     zcape(jl)=0._dp
     zmfub1(jl)=zmfub(jl)
  ENDDO
!
  DO jk=1,klev
     DO jl=1,kproma
        llo1=ldcum(jl).AND.ktype(jl).EQ.1
        IF(llo1.AND.jk.LE.kcbot(jl).AND.jk.GT.kctop(jl)) THEN
           ikb=kcbot(jl)
           zro=paphp1(jl,jk)/(rd*ztenh(jl,jk)*                         &
                                           (1._dp+vtmpc1*zqenh(jl,jk)))
           zdz=(paphp1(jl,jk)-paphp1(jl,jk-1))/(g*zro)
           zheat(jl)=zheat(jl) +                                        &
                (  (pten(jl,jk-1)-pten(jl,jk) + g*zdz/zcpcu(jl,jk))     &
                     /ztenh(jl,jk)                                      &
                    +  vtmpc1*(pqen(jl,jk-1)-pqen(jl,jk))  ) *          &
                       (g*(pmfu(jl,jk)+pmfd(jl,jk)))/zro
           zcape(jl)=zcape(jl) +                                        &
                         (g*(ptu(jl,jk)-ztenh(jl,jk))/ztenh(jl,jk)      &
                              +g*vtmpc1*(pqu(jl,jk)-zqenh(jl,jk))       &
                              -g*plu(jl,jk) ) * zdz
        ENDIF
     ENDDO
  ENDDO
!
  DO jl=1,kproma
     IF(ldcum(jl).AND.ktype(jl).EQ.1) THEN
        ikb=kcbot(jl)
        zmfub1(jl)=(zcape(jl)*zmfub(jl))/(zheat(jl)*ztau)
        zmfub1(jl)=MAX(zmfub1(jl),0.001_dp)
        zmfmax=(paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
        zmfub1(jl)=MIN(zmfub1(jl),zmfmax)
     ENDIF
  ENDDO
!
!*                 RECALCULATE CONVECTIVE FLUXES DUE TO EFFECT OF
!*                 DOWNDRAFTS ON BOUNDARY LAYER MOISTURE BUDGET
!*                 FOR SHALLOW CONVECTION (KTYPE=2)
!                  --------------------------------------------
!
!DIR$ IVDEP
  DO 520 jl=1,kproma
     IF(ktype(jl).EQ.2) THEN
        ikb=kcbot(jl)
        llo1=pmfd(jl,ikb).LT.0._dp.AND.loddraf(jl)
        zeps=MERGE(cmfdeps,0._dp,llo1)
        zqumqe=pqu(jl,ikb)+plu(jl,ikb)-                                 &
                    zeps*zqd(jl,ikb)-(1._dp-zeps)*zqenh(jl,ikb)
        zdqmin=MAX(0.01_dp*zqenh(jl,ikb),1.e-10_dp)
        zmfmax=(paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
        llo1=zdqpbl(jl).GT.0._dp.AND.zqumqe.GT.zdqmin.AND.ldcum(jl)    &
                    .AND.zmfub(jl).LT.zmfmax
        zmfub1(jl)=MERGE(zdqpbl(jl)/(g*MAX(zqumqe,zdqmin)),             &
                               zmfub(jl),llo1)
        zmfub1(jl)=MERGE(zmfub1(jl),zmfub(jl),                         &
                         ABS(zmfub1(jl)-zmfub(jl)).LT.0.2_dp*zmfub(jl))
     END IF
520 END DO
  DO 540 jk=1,klev
     DO 530 jl=1,kproma
        IF(ldcum(jl)) THEN
           zfac=zmfub1(jl)/MAX(zmfub(jl),1.e-10_dp)
           pmfd(jl,jk)=pmfd(jl,jk)*zfac
           zmfds(jl,jk)=zmfds(jl,jk)*zfac
           zmfdq(jl,jk)=zmfdq(jl,jk)*zfac
           zdmfdp(jl,jk)=zdmfdp(jl,jk)*zfac
        END IF
530  END DO
!
!---wiso-code

! Change water isotope fluxes after downdraft calculation and new cloud base mass flux
     DO jt=1,kwiso
       DO jl=1,kproma
         IF (ldcum(jl)) THEN
           zfac = zmfub1(jl)/MAX(zmfub(jl),1.E-10_dp)
           zwisomfdq(jl,jk,jt) = zwisomfdq(jl,jk,jt)*zfac
           zwisodmfdp(jl,jk,jt) = zwisodmfdp(jl,jk,jt)*zfac
         END IF
       END DO
     END DO

!---wiso-code-end

!re     DO 5304 jt=1,ktrac
!re        DO 5302 jl=1,kproma
!re           IF(ldcum(jl)) THEN
!re              zfac=zmfub1(jl)/MAX(zmfub(jl),1.e-10_dp)
!re              zmfdxt(jl,jk,jt)=zmfdxt(jl,jk,jt)*zfac
!re           ENDIF
!re5302    END DO
!re5304 END DO
!
540 END DO
!
!*                 NEW VALUES OF CLOUD BASE MASS FLUX
!                  ----------------------------------
!
  DO 550 jl=1,kproma
     IF(ldcum(jl)) zmfub(jl)=zmfub1(jl)
550 END DO
!
!-----------------------------------------------------------------------
!*    6.0          DETERMINE FINAL CLOUD ASCENT FOR ENTRAINING PLUME
!*                 FOR PENETRATIVE CONVECTION (TYPE=1),
!*                 FOR SHALLOW TO MEDIUM CONVECTION (TYPE=2)
!*                 AND FOR MID-LEVEL CONVECTION (TYPE=3).
!                  --------------------------------------------------
!
600 CONTINUE
  CALL h2oiso_cuasc(kproma, kbdim, klev, klevp1, klevm1,                      &
             ztenh,    zqenh,    puen,     pven,                       &
!re             ktrac,                                                    &
!re             zxtenh,   pxten,    pxtu,     zmfuxt,                     &
             pten,     pqen,     pqsen,                                &
             pgeo,     zgeoh,    paphp1,                               &
             pqte,     pverv,    ilwmin,                               &
             ldcum,    ldland,   ktype,    ilab,                       &
             ptu,      pqu,      plu,      zuu,      zvu,              &
             pmfu,     zmfub,    zentr,                                &
             zmfus,    zmfuq,                                          &
             zmful,    plude,    pqude,    zdmfup,                     &
             ihmin,    zhhatt,   zhcbase,  zqsenh,                     &
             zcpen,    zcpcu,                                          &
             kcbot,    kctop,    ictop0,                               &
!---wiso-code
             kwiso,                                                    &
             zwisoqenh,                                                &
             pwisoqen,                                                 &
             pwisoqu,  pwisolu,                                        &
             zwisomfuq,                                                &
             zwisomful,pwisolude,pwisoqude,zwisodmfup                  &
             , nn, time_step_len)
!---wiso-code-end
!
!-----------------------------------------------------------------------
!*    7.0          DETERMINE FINAL CONVECTIVE FLUXES IN 'CUFLX'
!                  ------------------------------------------
!
700 CONTINUE
  CALL h2oiso_cuflx(kproma,   kbdim,    klev,     klevp1,                     &
             pqen,     pqsen,    ztenh,    zqenh,                      &
!re             ktrac,                                                    &
!re             zxtenh,   zmfuxt,   zmfdxt,                               &
             paphp1,   zgeoh,                                          &
             kcbot,    kctop,    idtop,                                &
             ktype,    loddraf,  ldcum,                                &
             pmfu,     pmfd,     zmfus,    zmfds,                      &
             zmfuq,    zmfdq,    zmful,                                &
             zdmfup,   zdmfdp,   zrfl,     prain,                      &
             zcpcu,                                                    &
             pten,     zsfl,     zdpmel,   itopm2,                     &
!---wiso-code
             kwiso,                                                    &
             zwisoqenh,                                                &
             zwisomfuq,zwisomfdq,zwisomful,                            &
             zwisodmfdp,                                               &
             zmelt_tmp,zprec_frac,zrain_tmp,                           &
             zrfl_tmp, zsfl_tmp                                        &
             , time_step_len )
!---wiso-code-end
!
!-----------------------------------------------------------------------
!
!*    8.0          UPDATE TENDENCIES FOR T AND Q IN SUBROUTINE CUDTDQ
!                  --------------------------------------------------
!
!---wiso-code

!  Calculate old values (T-1)
  DO jk=1,klev
    DO jl=1,kproma
      zqsto(jl,jk)=pqen(jl,jk)-pqte(jl,jk)*time_step_len
      ztsto(jl,jk)=pten(jl,jk)-ptte(jl,jk)*time_step_len
    END DO
  END DO
  DO jt=1,kwiso
    DO jk=1,klev
      DO jl=1,kproma
        zwisoqsto(jl,jk,jt)=pwisoqen(jl,jk,jt)-pwisoqte(jl,jk,jt)*time_step_len
      END DO
    END DO
  END DO

!---wiso-code-end

800 CONTINUE
  CALL h2oiso_cudtdq(kproma, kbdim, klev, klevp1, itopm2, ldcum, &!re, ktrac,       &
              paphp1,   pten,     ptte,     pqte,                      &
              pxtec, &!re,    pxtte,    zmfuxt,   zmfdxt,                    &
              zmfus,    zmfds,    zmfuq,    zmfdq,                     &
              zmful,    zdmfup,   zdmfdp,   plude,                     &
              zdpmel,   zrfl,     zsfl,                                &
              zcpen,    pqtec,    pqude,                               &
              prsfc,    pssfc,    paprc,    paprs, delta_time)
!
!---wiso-code

!  update water isotope tendencies
  CALL h2oiso_cuwisodq(kproma, kbdim, klev, klevp1, itopm2, ldcum, kwiso,     &
                paphp1,   pwisoqte,                                    &
                pwisoxtec,                                             &
                zwisomfuq,zwisomfdq,                                   &
                zwisomful,zwisodmfup,zwisodmfdp,pwisolude,             &
                pwisoqtec,pwisoqude, delta_time)

! calculate precipitation of water isotopes and
! get precipitation into equilibrium with surrounding vapour 
  CALL h2oiso_cuwisoequ(kproma, kbdim, klev, klevp1, itopm2, ldcum, kwiso,    &
                         kctop,          kcbot,       paphp1,                  &
                         ztsto,          ptte,                                 & 
                         zqsto,          pqte,                                 &
                         zwisoqsto,      pwisoqte,                             &
                         ptu,            pqu,         pwisoqu,                 &
                           pten,                                                 &
                                 zrfl_tmp,       zsfl_tmp,                             &
                         zwisorfl,       zwisosfl,                             &
                         zwisodmfup,     zwisodmfdp,                           &
                         zmelt_tmp,                                            &
                         zprec_frac,     zrain_tmp,                            &
                         pwisorsfc,      pwisossfc,                            &
                         pwisoaprc,      pwisoaprs,                            &
                                 time_step_len, delta_time )
                     
!---wiso-code-end


!-----------------------------------------------------------------------
!
!*    9.0          UPDATE TENDENCIES FOR U AND U IN SUBROUTINE CUDUDV
!                  --------------------------------------------------
!
!!$900 CONTINUE
!!$  IF(lmfdudv) THEN
!!$     CALL h2oiso_cududv(kproma,   kbdim,    klev,     klevp1,                 &
!!$                 itopm2,   ktype,    kcbot,    paphp1,   ldcum,        &
!!$                 puen,     pven,     pvom,     pvol,                   &
!!$                 zuu,      zud,      zvu,      zvd,                    &
!!$                 pmfu,     pmfd)
!!$!
!!$  END IF
!!$!
!!$1000 CONTINUE
!

  RETURN
END SUBROUTINE h2oiso_cumastr

! ---------------------------------------------------------------------------------


!============================================================================
!============================================================================


! ---------------------------------------------------------------------------------

SUBROUTINE h2oiso_cumastrt( kproma, kbdim, klev, klevp1, klevm1, ilab,        &
           pten,     pqen,     pxen,     puen,     pven,               &
           ptven,    ldland, &!re,    ktrac,                                 &
!re           pxten,    pxtu,     pxtte,                                  &
           pverv,    pqsen,    pqhfla,                                 &
           paphp1,   pgeo,                                             &
           ptte,     pqte,     pvom,     pvol,                         &
           prsfc,    pssfc,    paprc,    paprs,    pxtec,              &
           pqtec,    pqude,                                            &
           ldcum,    ktype,    kcbot,    kctop,                        &
           ptu,      pqu,      plu,      plude,                        &
           pmfu,     pmfd,     prain,                                  &
!---wiso-code
           kwiso,                                                      &
           pwisoqen, pwisoxen,                                         &
           pwisoqte,                                                   &
           pwisorsfc,pwisossfc,pwisoaprc,pwisoaprs,pwisoxtec,          &
           pwisoqtec,pwisoqude,                                        &
           pwisoqu,  pwisolu,  pwisolude,                              &
           pwisorain, &
           nn, time_step_len, lookupoverflow, delta_time)
!---wiso-code-end
!
!**** *CUMASTRT*  MASTER ROUTINE FOR CUMULUS MASSFLUX-SCHEME
!
!     M.TIEDTKE      E.C.M.W.F.     1986/1987/1989
!
!     PURPOSE
!     -------
!
!          THIS ROUTINE COMPUTES THE PHYSICAL TENDENCIES OF THE
!     PROGNOSTIC VARIABLES T,Q,U AND V DUE TO CONVECTIVE PROCESSES.
!     PROCESSES CONSIDERED ARE: CONVECTIVE FLUXES, FORMATION OF
!     PRECIPITATION, EVAPORATION OF FALLING RAIN BELOW CLOUD BASE,
!     SATURATED CUMULUS DOWNDRAFTS.
!
!**   INTERFACE.
!     ----------
!
!          *CUMASTRT* IS CALLED FROM *CUCALL* IN CASE ICONV .EQ. 2
!     THE ROUTINE TAKES ITS INPUT FROM THE LONG-TERM STORAGE
!     T,Q,U,V,PHI AND P AND MOISTURE TENDENCIES.
!     IT RETURNS ITS OUTPUT TO THE SAME SPACE
!      1.MODIFIED TENDENCIES OF MODEL VARIABLES
!      2.RATES OF CONVECTIVE PRECIPITATION
!        (USED IN SUBROUTINE SURF)
!
!     METHOD
!     -------
!
!     PARAMETERIZATION IS DONE USING A MASSFLUX-SCHEME.
!        (1) DEFINE CONSTANTS AND PARAMETERS
!        (2) SPECIFY VALUES (T,Q,QS...) AT HALF LEVELS AND
!            INITIALIZE UPDRAFT- AND DOWNDRAFT-VALUES IN 'CUINI'
!        (3) CALCULATE CLOUD BASE IN 'CUBASE'
!            AND SPECIFY CLOUD BASE MASSFLUX FROM PBL MOISTURE BUDGET
!        (4) DO CLOUD ASCENT IN 'CUASC' IN ABSENCE OF DOWNDRAFTS
!        (5) DO DOWNDRAFT CALCULATIONS:
!              (A) DETERMINE VALUES AT LFS IN 'CUDLFS'
!              (B) DETERMINE MOIST DESCENT IN 'CUDDRAF'
!              (C) RECALCULATE CLOUD BASE MASSFLUX CONSIDERING THE
!                  EFFECT OF CU-DOWNDRAFTS
!        (6) DO FINAL CLOUD ASCENT IN 'CUASC'
!        (7) DO FINAL ADJUSMENTS TO CONVECTIVE FLUXES IN 'CUFLX',
!            DO EVAPORATION IN SUBCLOUD LAYER
!        (8) CALCULATE INCREMENTS OF T AND Q IN 'CUDTDQ'
!        (9) CALCULATE INCREMENTS OF U AND V IN 'CUDUDV'
!
!     EXTERNALS.
!     ----------
!
!       CUINI:  INITIALIZES VALUES AT VERTICAL GRID USED IN CU-PARAMETR.
!       CUBASE: CLOUD BASE CALCULATION FOR PENETR.AND SHALLOW CONVECTION
!       CUASCT: CLOUD ASCENT FOR ENTRAINING PLUME
!       CUDLFS: DETERMINES VALUES AT LFS FOR DOWNDRAFTS
!       CUDDRAF:DOES MOIST DESCENT FOR CUMULUS DOWNDRAFTS
!       CUFLX:  FINAL ADJUSTMENTS TO CONVECTIVE FLUXES (ALSO IN PBL)
!       CUDQDT: UPDATES TENDENCIES FOR T AND Q
!       CUDUDV: UPDATES TENDENCIES FOR U AND V
!
!     SWITCHES.
!     --------
!
!          LMFPEN=.T.   PENETRATIVE CONVECTION IS SWITCHED ON
!          LMFSCV=.T.   SHALLOW CONVECTION IS SWITCHED ON
!          LMFMID=.T.   MIDLEVEL CONVECTION IS SWITCHED ON
!          LMFDD=.T.    CUMULUS DOWNDRAFTS SWITCHED ON
!          LMFDUDV=.T.  CUMULUS FRICTION SWITCHED ON
!
!     MODEL PARAMETERS (DEFINED IN MODULE mo_cumulus_flux)
!     ----------------------------------------------------
!     ENTRPEN    ENTRAINMENT RATE FOR PENETRATIVE CONVECTION
!     ENTRSCV    ENTRAINMENT RATE FOR SHALLOW CONVECTION
!     ENTRMID    ENTRAINMENT RATE FOR MIDLEVEL CONVECTION
!     ENTRDD     ENTRAINMENT RATE FOR CUMULUS DOWNDRAFTS
!     CMFCTOP    RELATIVE CLOUD MASSFLUX AT LEVEL ABOVE NONBUOYANCY LEVEL
!     CMFCMAX    MAXIMUM MASSFLUX VALUE ALLOWED FOR
!     CMFCMIN    MINIMUM MASSFLUX VALUE (FOR SAFETY)
!     CMFDEPS    FRACTIONAL MASSFLUX FOR DOWNDRAFTS AT LFS
!     CPRCON     COEFFICIENT FOR CONVERSION FROM CLOUD WATER TO RAIN
!
!     REFERENCE.
!     ----------
!
!          PAPER ON MASSFLUX SCHEME (TIEDTKE,1989)
!
!re USE mo_kind,           ONLY: dp
!re USE mo_control,        ONLY: nn
!re USE mo_constants,      ONLY: g, alv, als, tmelt, vtmpc1, rd
!re USE mo_cumulus_flux,   ONLY: entrpen, entrscv, lmfdd, cmfdeps, lmfdudv
!re USE mo_convect_tables, ONLY: tlucua,                         & ! table a
!re                              tlucub,                         & ! table b
!re                              jptlucu1, jptlucu2,                       &
!re                              lookuperror, lookupoverflow
!re USE mo_time_control,   ONLY: time_step_len


   USE messy_main_constants_mem,     ONLY: g, alv, als, tmelt, vtmpc1, rd

   USE MESSY_CONVECT_TIEDTKE_PARAM,  ONLY: entrpen, entrscv, lmfdd, cmfdeps, lmfdudv

   USE messy_main_tools,             ONLY: tlucua, tlucub, jptlucu1, jptlucu2


!
IMPLICIT NONE
!
!
INTEGER, INTENT (IN) :: nn 
REAL(DP),INTENT(IN):: time_step_len!re 
LOGICAL, INTENT(INOUT):: lookupoverflow!re 

INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, klevm1!re, ktrac

REAL(DP), INTENT(IN):: delta_time!re
!
!---wiso-code

  INTEGER, INTENT (IN) :: kwiso

!---wiso-code-end

REAL(dp):: pten(kbdim,klev),        pqen(kbdim,klev),                  &
           pxen(kbdim,klev),        ptven(kbdim,klev),                 &
           puen(kbdim,klev),        pven(kbdim,klev),                  &
           ptte(kbdim,klev),        pqte(kbdim,klev),                  &
           pvom(kbdim,klev),        pvol(kbdim,klev),                  &
           pqsen(kbdim,klev),       pgeo(kbdim,klev),                  &
           paphp1(kbdim,klevp1),                                       &
           pverv(kbdim,klev)
REAL(dp):: ptu(kbdim,klev),         pqu(kbdim,klev),                   &
           plu(kbdim,klev),         plude(kbdim,klev),                 &
           pmfu(kbdim,klev),        pmfd(kbdim,klev),                  &
           paprc(kbdim),            paprs(kbdim),                      &
           prsfc(kbdim),            pssfc(kbdim),                      &
           prain(kbdim),            pqhfla(kbdim)
INTEGER :: kcbot(kbdim),            kctop(kbdim),                      &
           ktype(kbdim)
REAL(dp):: pxtec(kbdim,klev),       pqtec(kbdim,klev),                 &
           pqude(kbdim,klev)

!---wiso-code

REAL(dp):: pwisoqen(kbdim,klev,kwiso),                                 &
           pwisoxen(kbdim,klev,kwiso),                                 &
           pwisoqte(kbdim,klev,kwiso)
REAL(dp):: pwisoqu(kbdim,klev,kwiso),                                  &
           pwisolu(kbdim,klev,kwiso), pwisolude(kbdim,klev,kwiso),     &
           pwisoaprc(kbdim,kwiso),    pwisoaprs(kbdim,kwiso),          &
           pwisorsfc(kbdim,kwiso),    pwisossfc(kbdim,kwiso),          &
           pwisorain(kbdim,kwiso)
REAL(dp):: pwisoxtec(kbdim,klev,kwiso),pwisoqtec(kbdim,klev,kwiso),    &
           pwisoqude(kbdim,klev,kwiso)

!---wiso-code-end

REAL(dp):: ztenh(kbdim,klev),       zqenh(kbdim,klev),                 &
           zxenh(kbdim,klev),                                          &
           zgeoh(kbdim,klev),       zqsenh(kbdim,klev),                &
           ztd(kbdim,klev),         zqd(kbdim,klev),                   &
           zmfus(kbdim,klev),       zmfds(kbdim,klev),                 &
           zmfuq(kbdim,klev),       zmfdq(kbdim,klev),                 &
           zdmfup(kbdim,klev),      zdmfdp(kbdim,klev),                &
           zmful(kbdim,klev),       zrfl(kbdim),                       &
           zuu(kbdim,klev),         zvu(kbdim,klev),                   &
           zud(kbdim,klev),         zvd(kbdim,klev)
REAL(dp):: zcpen(kbdim,klev),       zcpcu(kbdim,klev)
REAL(dp):: zentr(kbdim),            zhcbase(kbdim),                    &
           zmfub(kbdim),            zmfub1(kbdim),                     &
           zdqpbl(kbdim),           zdqcv(kbdim)
REAL(dp):: zsfl(kbdim),             zdpmel(kbdim,klev)
INTEGER :: ilab(kbdim,klev),        idtop(kbdim),                      &
           ictop0(kbdim),           ilwmin(kbdim)
!re REAL(dp):: pxten(kbdim,klev,ktrac), pxtte(kbdim,klev,ktrac),           &
!re           pxtu(kbdim,klev,ktrac),                                     &
!re           zxtenh(kbdim,klev,ktrac),zxtd(kbdim,klev,ktrac),            &
!re           zmfuxt(kbdim,klev,ktrac),zmfdxt(kbdim,klev,ktrac)
LOGICAL :: loddraf(kbdim),          ldland(kbdim)
LOGICAL :: ldcum(kbdim)
LOGICAL :: llo1, lo
!
INTEGER :: jl, jk, ikb, it, it1, jt, itopm2
REAL(dp):: zcons2, zqumqe, zdqmin, zmfmax, zalvs, zalvdcp, zqalv       &
         , zhsat, zes, zcor, zqsat, zqst1, zdqsdt, zgam, zzz, zhhat    &
         , zfac, zpbmpt, zeps

!---wiso-code

REAL(dp):: zwisoqenh(kbdim,klev,kwiso),                                &
           zwisoxenh(kbdim,klev,kwiso),                                &
           zwisoqd(kbdim,klev,kwiso),                                  &
           zwisomfuq(kbdim,klev,kwiso),zwisomfdq(kbdim,klev,kwiso),    &
           zwisodmfup(kbdim,klev,kwiso),zwisodmfdp(kbdim,klev,kwiso),  &
           zwisomful(kbdim,klev,kwiso)

REAL(dp):: zmelt_tmp(kbdim,klev),        &  ! temporary value
           zprec_frac(kbdim,klev),       &  ! temporary value
           zrain_tmp(kbdim,klev),        &  ! temporary value 
           zrfl_tmp(kbdim),              &  ! temporary value
           zsfl_tmp(kbdim),              &  ! temporary value
           zqsto(kbdim,klev),            &  ! old moisture value
           ztsto(kbdim,klev),            &  ! old temperature value
           zwisoqsto(kbdim,klev,kwiso),  &  ! old water isotope moisture value
           zwisorfl(kbdim,kwiso),        &  ! temporary value
           zwisosfl(kbdim,kwiso)            ! temporary value

!---wiso-code-end
!
!  INTRINSIC FUNCTIONS
INTRINSIC MIN, MAX
!
!  Executable statements

  lookupoverflow = .FALSE.
!---------------------------------------------------------------------
!
!     1.           SPECIFY CONSTANTS AND PARAMETERS
!                  --------------------------------
!
100 CONTINUE
!
  zcons2=1._dp/(g*time_step_len)

!---wiso-code

  zprec_frac(:,:)=0._dp
  zrain_tmp(:,:)=0._dp

!---wiso-code-end
!
!---------------------------------------------------------------------
!
!*    2.           INITIALIZE VALUES AT VERTICAL GRID POINTS IN 'CUINI'
!                  ---------------------------------------------------
!
200 CONTINUE
  CALL h2oiso_cuini(kproma,   kbdim,    klev,     klevp1,   klevm1,           &
             pten,     pqen,     pqsen,    pxen,     puen,     pven,   &
             ptven, &!re,    ktrac,                                          &
!re             pxten,    zxtenh,   pxtu,     zxtd,     zmfuxt,   zmfdxt, &
             pverv,    pgeo,     paphp1,   zgeoh,                      &
             ztenh,    zqenh,    zqsenh,   zxenh,    ilwmin,           &
             ptu,      pqu,      ztd,      zqd,                        &
             zuu,      zvu,      zud,      zvd,                        &
             pmfu,     pmfd,     zmfus,    zmfds,                      &
             zmfuq,    zmfdq,    zdmfup,   zdmfdp,                     &
             zcpen,    zcpcu,                                          &
             zdpmel,   plu,      plude,    pqude,    ilab,             &
!---wiso-code
             kwiso,                                                    &
             pwisoqen, pwisoxen,                                       &
             zwisoqenh,zwisoxenh,                                      &
             pwisoqu,  zwisoqd,                                        &
             zwisomfuq,zwisomfdq,zwisodmfup,zwisodmfdp,                &
             pwisolu,  pwisolude,pwisoqude                             &
             )         
!---wiso-code-end
    IF (lookupoverflow) RETURN!re
!
!---------------------------------------------------------------------
!
!*    3.0          CLOUD BASE CALCULATIONS
!                  -----------------------
!
300 CONTINUE
!
!*             (A) DETERMINE CLOUD BASE VALUES IN 'CUBASE'
!                  ---------------------------------------
!
  CALL h2oiso_cubase(kproma, kbdim, klev, klevp1, klevm1,                     &
              ztenh,    zqenh,    zgeoh,    paphp1,                    &
              ptu,      pqu,      plu,                                 &
              puen,     pven,     zuu,      zvu,                       &
              zcpcu,                                                   &
              ldcum,    kcbot,    ilab,                                &
!---wiso-code
              kwiso,                                                   &
              pwisoqu,  pwisolu                                        &
              )
!---wiso-code-end
    IF (lookupoverflow) RETURN!re
!
!*             (B) DETERMINE TOTAL MOISTURE CONVERGENCE AND
!*                 THEN DECIDE ON TYPE OF CUMULUS CONVECTION
!                  -----------------------------------------
!
  jk=1
  DO 310 jl=1,kproma
     zdqpbl(jl)=0.0_dp
     zdqcv(jl)=pqte(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))
     idtop(jl)=0
310 END DO
  DO 320 jk=2,klev
     DO 315 jl=1,kproma
        zdqcv(jl)=zdqcv(jl)+pqte(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))
        IF(jk.GE.kcbot(jl)) zdqpbl(jl)=zdqpbl(jl)+pqte(jl,jk)          &
                             *(paphp1(jl,jk+1)-paphp1(jl,jk))
315  END DO
320 END DO
!
!*             (C) DETERMINE MOISTURE SUPPLY FOR BOUNDARY LAYER
!*                 AND DETERMINE CLOUD BASE MASSFLUX IGNORING
!*                 THE EFFECTS OF DOWNDRAFTS AT THIS STAGE
!                  ------------------------------------------
!
!DIR$ IVDEP
  DO 340 jl=1,kproma
     ikb=kcbot(jl)
     zqumqe=pqu(jl,ikb)+plu(jl,ikb)-zqenh(jl,ikb)
     zdqmin=MAX(0.01_dp*zqenh(jl,ikb),1.e-10_dp)
     llo1=zdqpbl(jl).GT.0._dp.AND.zqumqe.GT.zdqmin.AND.ldcum(jl)
     zmfub(jl)=MERGE(zdqpbl(jl)/(g*MAX(zqumqe,zdqmin)),0.01_dp,llo1)
     zmfmax=(paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
     zmfub(jl)=MIN(zmfub(jl),zmfmax)
     IF(.NOT.llo1) ldcum(jl)=.FALSE.
     ktype(jl)=MERGE(1,2,zdqcv(jl).GT.MAX(0._dp,-1.1_dp*pqhfla(jl)*g))
     zentr(jl)=MERGE(entrpen,entrscv,ktype(jl).EQ.1)
340 END DO
!
!---------------------------------------------------------------------
!*    4.0          DETERMINE CLOUD ASCENT FOR ENTRAINING PLUME
!                  -------------------------------------------
!
400 CONTINUE
!
!*             (A) ESTIMATE CLOUD HEIGHT FOR ENTRAINMENT/DETRAINMENT
!*                 CALCULATIONS IN CUASC (MAX.POSSIBLE CLOUD HEIGHT
!*                 FOR NON-ENTRAINING PLUME, FOLLOWING A.-S.,1974)
!                  -------------------------------------------------
!
!DIR$ IVDEP
  DO 410 jl=1,kproma
     ikb=kcbot(jl)
     zalvs=MERGE(alv,als,ptu(jl,ikb)>tmelt)
     zhcbase(jl)=zcpcu(jl,ikb)*ptu(jl,ikb)+zgeoh(jl,ikb)+zalvs*pqu(jl,ikb)
     ictop0(jl)=kcbot(jl)-1
410 END DO
  DO 430 jk=klevm1,3,-1
!DIR$ IVDEP
     DO 420 jl=1,kproma
        zalvs=MERGE(alv,als,ztenh(jl,jk)>tmelt)
        zalvdcp=zalvs/zcpcu(jl,jk)
        zqalv=1._dp/zalvs
        zhsat=zcpcu(jl,jk)*ztenh(jl,jk)+zgeoh(jl,jk)+zalvs*zqsenh(jl,jk)
        it = NINT(ztenh(jl,jk)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zes=tlucua(it)/paphp1(jl,jk)
        zes=MIN(0.5_dp,zes)
        LO=zes<0.4_dp
        zcor=1._dp/(1._dp-vtmpc1*zes)
        zqsat=zes*zcor
        it1=it+1
        it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
        zqst1=tlucua(it1)/paphp1(jl,jk)
        zqst1=MIN(0.5_dp,zqst1)
        zqst1=zqst1/(1._dp-vtmpc1*zqst1)
        zdqsdt=(zqst1-zqsat)*1000._dp
        zgam=MERGE(zalvdcp*zdqsdt,zqsat*zcor*tlucub(it),LO)
        zzz=zcpcu(jl,jk)*ztenh(jl,jk)*vtmpc1
        zhhat=zhsat-(zzz+zgam*zzz)/(1._dp+zgam*zzz*zqalv)*             &
                       MAX(zqsenh(jl,jk)-zqenh(jl,jk),0._dp)
        IF(jk.LT.ictop0(jl).AND.zhcbase(jl).GT.zhhat) ictop0(jl)=jk
420  END DO
430 END DO
!
     IF (lookupoverflow) RETURN !reCALL lookuperror ('h2oiso_cumastrt')
!!
!!     DEEP CONVECTION IF CLOUD DEPTH > 200 HPA, ELSE SHALLOW
!!     (CLOUD DEPTH FROM NON-ENTRAINIG PLUME)
!!
!  DO jl=1,kproma
!     ktype(jl)=MERGE(1,2,                                             &
!                paphp1(jl,kcbot(jl))-paphp1(jl,ictop0(jl)).gt.2.E4_dp)
!     zentr(jl)=MERGE(entrpen,entrscv,ktype(jl).eq.1)
!  ENDDO
!!
!*             (B) DO ASCENT IN 'CUASCT' IN ABSENCE OF DOWNDRAFTS
!                  ----------------------------------------------
!
  CALL h2oiso_cuasct(kproma, kbdim, klev, klevp1, klevm1,                     &
             ztenh,    zqenh,    puen,     pven,                       &
!re             ktrac,                                                    &
!re             zxtenh,   pxten,    pxtu,     zmfuxt,                     &
             pten,     pqen,     pqsen,                                &
             pgeo,     zgeoh,    paphp1,                               &
             pqte,     pverv,    ilwmin,                               &
             ldcum,    ldland,   ktype,    ilab,                       &
             ptu,      pqu,      plu,      zuu,      zvu,              &
             pmfu,     zmfub,    zentr,                                &
             zmfus,    zmfuq,                                          &
             zmful,    plude,    pqude,    zdmfup,                     &
             zcpen,    zcpcu,                                          &
             kcbot,    kctop,    ictop0,                               &
!---wiso-code
             kwiso,                                                    &
             zwisoqenh,                                                &
             pwisoqen,                                                 &
             pwisoqu,  pwisolu,                                        &
             zwisomfuq,                                                &
             zwisomful,pwisolude,pwisoqude,zwisodmfup                  &
             , nn, time_step_len)
!---wiso-code-end
!
!*        (C) CHECK CLOUD DEPTH AND CHANGE ENTRAINMENT RATE ACCORDINGLY
!             CALCULATE PRECIPITATION RATE (FOR DOWNDRAFT CALCULATION)
!             ---------------------------------------------------------
!
!DIR$ IVDEP
  DO 440 jl=1,kproma
     zpbmpt=paphp1(jl,kcbot(jl))-paphp1(jl,kctop(jl))
     IF(ldcum(jl).AND.ktype(jl).EQ.1.AND.zpbmpt.LT.2.e4_dp) ktype(jl)=2
     IF(ldcum(jl)) ictop0(jl)=kctop(jl)
     IF(ktype(jl).EQ.2) zentr(jl)=entrscv
     zrfl(jl)=zdmfup(jl,1)
440 END DO
  DO 460 jk=2,klev
     DO 450 jl=1,kproma
        zrfl(jl)=zrfl(jl)+zdmfup(jl,jk)
450  END DO
460 END DO
!
!---------------------------------------------------------------------
!*    5.0          CUMULUS DOWNDRAFT CALCULATIONS
!                  ------------------------------
!
500 CONTINUE
!
  IF(lmfdd) THEN
!
!*             (A) DETERMINE LFS IN 'CUDLFS'
!                  -------------------------
!
     CALL h2oiso_cudlfs(kproma,   kbdim,    klev,     klevp1,                 &
                 ztenh,    zqenh,    puen,     pven,                   &
!re                 ktrac,                                                &
!re                 zxtenh,   pxtu,     zxtd,     zmfdxt,                 &
                 zgeoh,    paphp1,                                     &
                 ptu,      pqu,      zuu,      zvu,                    &
                 ldcum,    kcbot,    kctop,    zmfub,    zrfl,         &
                 ztd,      zqd,      zud,      zvd,                    &
                 pmfd,     zmfds,    zmfdq,    zdmfdp,                 &
                 zcpcu,                                                &
                 idtop,    loddraf,                                    &
!---wiso-code
                 kwiso,                                                &
                 zwisoqenh,                                            &
                 pwisoqu,                                              &
                 zwisoqd,                                              &
                 zwisomfdq,zwisodmfdp,                                 &
                 zdmfup,   zwisodmfup                                  &
                 )
!---wiso-code-end
!
!*            (B)  DETERMINE DOWNDRAFT T,Q AND FLUXES IN 'CUDDRAF'
!                 -----------------------------------------------
!
     CALL h2oiso_cuddraf(kproma,   kbdim,    klev,     klevp1,                &
                  ztenh,    zqenh,    puen,     pven,                  &
!re                  ktrac,                                               &
!re                  zxtenh,   zxtd,     zmfdxt,                          &
                  zgeoh,    paphp1,   zrfl,                            &
                  ztd,      zqd,      zud,      zvd,                   &
                  pmfd,     zmfds,    zmfdq,    zdmfdp,                &
                  zcpcu,                                               &
                  loddraf,                                             &
!---wiso-code
                  kwiso,                                               &
                  zwisoqenh,                                           &
                  zwisoqd,                                             &
                  zwisomfdq,zwisodmfdp,                                &
                  zdmfup,   zwisodmfup                                 &
                  )
!---wiso-code-end
!
!*            (C)  RECALCULATE CONVECTIVE FLUXES DUE TO EFFECT OF
!*                 DOWNDRAFTS ON BOUNDARY LAYER MOISTURE BUDGET
!                  --------------------------------------------
!
!DIR$ IVDEP
  DO 520 jl=1,kproma
     IF(loddraf(jl)) THEN
        ikb=kcbot(jl)
        llo1=pmfd(jl,ikb).LT.0._dp
        zeps=MERGE(cmfdeps,0._dp,llo1)
        zqumqe=pqu(jl,ikb)+plu(jl,ikb)-                                &
                    zeps*zqd(jl,ikb)-(1._dp-zeps)*zqenh(jl,ikb)
        zdqmin=MAX(0.01_dp*zqenh(jl,ikb),1.e-10_dp)
        zmfmax=(paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
        llo1=zdqpbl(jl).GT.0._dp.AND.zqumqe.GT.zdqmin.AND.ldcum(jl)    &
                    .AND.zmfub(jl).LT.zmfmax
        zmfub1(jl)=MERGE(zdqpbl(jl)/(g*MAX(zqumqe,zdqmin)),            &
                               zmfub(jl),llo1)
        zmfub1(jl)=MERGE(zmfub1(jl),zmfub(jl),                         &
                           (ktype(jl).EQ.1.OR.ktype(jl).EQ.2) .AND.    &
                         ABS(zmfub1(jl)-zmfub(jl)).LT.0.2_dp*zmfub(jl))
     END IF
520 END DO
  DO 540 jk=1,klev
     DO 530 jl=1,kproma
        IF(loddraf(jl)) THEN
           zfac=zmfub1(jl)/MAX(zmfub(jl),1.e-10_dp)
           pmfd(jl,jk)=pmfd(jl,jk)*zfac
           zmfds(jl,jk)=zmfds(jl,jk)*zfac
           zmfdq(jl,jk)=zmfdq(jl,jk)*zfac
           zdmfdp(jl,jk)=zdmfdp(jl,jk)*zfac
        END IF
530  END DO
!
!---wiso-code

! Change water isotope fluxes after downdraft calculation and new cloud base mass flux
     DO jt=1,kwiso
       DO jl=1,kproma
         IF (loddraf(jl)) THEN
           zfac = zmfub1(jl)/MAX(zmfub(jl),1.e-10_dp)
           zwisomfdq(jl,jk,jt) = zwisomfdq(jl,jk,jt)*zfac
           zwisodmfdp(jl,jk,jt) = zwisodmfdp(jl,jk,jt)*zfac
         END IF
       END DO
     END DO

!---wiso-code-end

!re     DO 5304 jt=1,ktrac
!re        DO 5302 jl=1,kproma
!re           IF(loddraf(jl)) THEN
!re             zfac=zmfub1(jl)/MAX(zmfub(jl),1.e-10_dp)
!re              zmfdxt(jl,jk,jt)=zmfdxt(jl,jk,jt)*zfac
!re           ENDIF
!re5302    END DO
!re5304 END DO
!
540 END DO
!
!*                 NEW VALUES OF CLOUD BASE MASS FLUX
!                  ----------------------------------
  DO 550 jl=1,kproma
     IF(loddraf(jl)) zmfub(jl)=zmfub1(jl)
550 END DO
!
  END IF
!
!---------------------------------------------------------------------
!*    6.0          DETERMINE FINAL CLOUD ASCENT FOR ENTRAINING PLUME
!*                 FOR PENETRATIVE CONVECTION (TYPE=1),
!*                 FOR SHALLOW TO MEDIUM CONVECTION (TYPE=2)
!*                 AND FOR MID-LEVEL CONVECTION (TYPE=3).
!                  ----------------------------------------------------
!
600 CONTINUE
!
  CALL h2oiso_cuasct(kproma,  kbdim,    klev,     klevp1,   klevm1,           &
             ztenh,    zqenh,    puen,     pven,                       &
!re             ktrac,                                                    &
!re             zxtenh,   pxten,    pxtu,     zmfuxt,                     &
             pten,     pqen,     pqsen,                                &
             pgeo,     zgeoh,    paphp1,                               &
             pqte,     pverv,    ilwmin,                               &
             ldcum,    ldland,   ktype,    ilab,                       &
             ptu,      pqu,      plu,      zuu,      zvu,              &
             pmfu,     zmfub,    zentr,                                &
             zmfus,    zmfuq,                                          &
             zmful,    plude,    pqude,    zdmfup,                     &
             zcpen,    zcpcu,                                          &
             kcbot,    kctop,    ictop0,                               &
!---wiso-code
             kwiso,                                                    &
             zwisoqenh,                                                &
             pwisoqen,                                                 &
             pwisoqu,  pwisolu,                                        &
             zwisomfuq,                                                &
             zwisomful,pwisolude,pwisoqude,zwisodmfup                  &
             , nn, time_step_len)
!---wiso-code-end
!
!---------------------------------------------------------------------
!*    7.0          DETERMINE FINAL CONVECTIVE FLUXES IN 'CUFLX'
!                  ------------------------------------------
!
700 CONTINUE
  CALL h2oiso_cuflx(kproma,   kbdim,    klev,     klevp1,                     &
             pqen,     pqsen,    ztenh,    zqenh,                      &
!re             ktrac,                                                    &
!re             zxtenh,   zmfuxt,   zmfdxt,                               &
             paphp1,   zgeoh,                                          &
             kcbot,    kctop,    idtop,                                &
             ktype,    loddraf,  ldcum,                                &
             pmfu,     pmfd,     zmfus,    zmfds,                      &
             zmfuq,    zmfdq,    zmful,                                &
             zdmfup,   zdmfdp,   zrfl,     prain,                      &
             zcpcu,                                                    &
             pten,     zsfl,     zdpmel,   itopm2,                     &
!---wiso-code
             kwiso,                                                    &
             zwisoqenh,                                                &
             zwisomfuq,zwisomfdq,zwisomful,                            &
             zwisodmfdp,                                               &
             zmelt_tmp,zprec_frac,zrain_tmp,                           &
             zrfl_tmp, zsfl_tmp                                        &
             , time_step_len)
!---wiso-code-end
!
!---------------------------------------------------------------------
!
!*    8.0          UPDATE TENDENCIES FOR T AND Q IN SUBROUTINE CUDTDQ
!                  --------------------------------------------------
!
!---wiso-code

!  Calculate old values (T-1)
  DO jk=1,klev
    DO jl=1,kproma
      zqsto(jl,jk)=pqen(jl,jk)-pqte(jl,jk)*time_step_len
      ztsto(jl,jk)=pten(jl,jk)-ptte(jl,jk)*time_step_len
    END DO
  END DO
  DO jt=1,kwiso
    DO jk=1,klev
      DO jl=1,kproma
        zwisoqsto(jl,jk,jt)=pwisoqen(jl,jk,jt)-pwisoqte(jl,jk,jt)*time_step_len
      END DO
    END DO
  END DO

!---wiso-code-end

800 CONTINUE
  CALL h2oiso_cudtdq(kproma, kbdim, klev, klevp1, itopm2, ldcum, &!re, ktrac,       &
              paphp1,   pten,     ptte,     pqte,                      &
              pxtec, & !re,    pxtte,    zmfuxt,   zmfdxt,                    &
              zmfus,    zmfds,    zmfuq,    zmfdq,                     &
              zmful,    zdmfup,   zdmfdp,   plude,                     &
              zdpmel,   zrfl,     zsfl,                                &
              zcpen,    pqtec,    pqude,                               &
              prsfc,    pssfc,    paprc,    paprs, delta_time)
!
!---wiso-code

!  update water isotope tendencies
  CALL h2oiso_cuwisodq(kproma, kbdim, klev, klevp1, itopm2, ldcum, kwiso,     &
                paphp1,   pwisoqte,                                    &
                pwisoxtec,                                             &
                zwisomfuq,zwisomfdq,                                   &
                zwisomful,zwisodmfup,zwisodmfdp,pwisolude,             &
                pwisoqtec,pwisoqude, delta_time)

! calculate precipitation of water isotopes and
! get precipitation into equilibrium with surrounding vapour 
  CALL h2oiso_cuwisoequ(kproma, kbdim, klev, klevp1, itopm2, ldcum, kwiso,    &
                         kctop,          kcbot,       paphp1,                  &
                         ztsto,          ptte,                                 & 
                         zqsto,          pqte,                                 &
                         zwisoqsto,      pwisoqte,                             &
                         ptu,            pqu,         pwisoqu,                 &
                         pten,                                                 &
                                 zrfl_tmp,       zsfl_tmp,                             &
                         zwisorfl,       zwisosfl,                             &
                         zwisodmfup,     zwisodmfdp,                           &
                         zmelt_tmp,                                            &
                         zprec_frac,     zrain_tmp,                            &
                         pwisorsfc,      pwisossfc,                            &
                         pwisoaprc,      pwisoaprs,                            &
                                 time_step_len, delta_time )
                     
!---wiso-code-end

!---------------------------------------------------------------------
!
!*    9.0          UPDATE TENDENCIES FOR U AND U IN SUBROUTINE CUDUDV
!                  --------------------------------------------------
!
!re900 CONTINUE
!re  IF(lmfdudv) THEN
!re     CALL h2oiso_cududv(kproma,   kbdim,    klev,     klevp1,                 &
!re                itopm2,   ktype,    kcbot,    paphp1,   ldcum,        &
!re                 puen,     pven,     pvom,     pvol,                   &
!re                 zuu,      zud,      zvu,      zvd,                    &
!re                 pmfu,     pmfd)
!re!
!re  END IF
!
!re1000 CONTINUE
!
!
  RETURN
END SUBROUTINE h2oiso_cumastrt


! ---------------------------------------------------------------------------------




! ---------------------------------------------------------------------------------


SUBROUTINE h2oiso_cumastrh( kproma, kbdim, klev, klevp1, klevm1, ilab,        &
           pten,     pqen,     pxen,     puen,     pven,               &
           ptven,    ldland, & !re,    ktrac,                                 &
 !re          pxten,    pxtu,     pxtte,                                  &
           pverv,    pqsen,    pqhfla,                                 &
           paphp1,   pgeo,                                             &
           ptte,     pqte,     pvom,     pvol,                         &
           prsfc,    pssfc,    paprc,    paprs,    pxtec,              &
           pqtec,    pqude,                                            &
           ldcum,    ktype,    kcbot,    kctop,                        &
           ptu,      pqu,      plu,      plude,                        &
           pmfu,     pmfd,     prain,                                  &
!---wiso-code
           kwiso,                                                      &
           pwisoqen, pwisoxen,                                         &
           pwisoqte,                                                   &
           pwisorsfc,pwisossfc,pwisoaprc,pwisoaprs,pwisoxtec,          &
           pwisoqtec,pwisoqude,                                        &
           pwisoqu,  pwisolu,  pwisolude,                              &
           pwisorain, &
           nn, time_step_len,  lookupoverflow, delta_time)
!---wiso-code-end
!
!**** *CUMASTRH*  MASTER ROUTINE FOR CUMULUS MASSFLUX-SCHEME
!
!     M.TIEDTKE      E.C.M.W.F.     1986/1987/1989
!
!     PURPOSE
!     -------
!
!          THIS ROUTINE COMPUTES THE PHYSICAL TENDENCIES OF THE
!     PROGNOSTIC VARIABLES T,Q,U AND V DUE TO CONVECTIVE PROCESSES.
!     PROCESSES CONSIDERED ARE: CONVECTIVE FLUXES, FORMATION OF
!     PRECIPITATION, EVAPORATION OF FALLING RAIN BELOW CLOUD BASE,
!     SATURATED CUMULUS DOWNDRAFTS.
!
!**   INTERFACE.
!     ----------
!
!          *CUMASTRH* IS CALLED FROM *CUCALL* IN CASE ICONV .EQ. 3
!     THE ROUTINE TAKES ITS INPUT FROM THE LONG-TERM STORAGE
!     T,Q,U,V,PHI AND P AND MOISTURE TENDENCIES.
!     IT RETURNS ITS OUTPUT TO THE SAME SPACE
!      1.MODIFIED TENDENCIES OF MODEL VARIABLES
!      2.RATES OF CONVECTIVE PRECIPITATION
!        (USED IN SUBROUTINE SURF)
!
!     METHOD
!     -------
!
!     PARAMETERIZATION IS DONE USING A MASSFLUX-SCHEME.
!        (1) DEFINE CONSTANTS AND PARAMETERS
!        (2) SPECIFY VALUES (T,Q,QS...) AT HALF LEVELS AND
!            INITIALIZE UPDRAFT- AND DOWNDRAFT-VALUES IN 'CUINI'
!        (3) CALCULATE CLOUD BASE IN 'CUBASE'
!            AND SPECIFY CLOUD BASE MASSFLUX FROM PBL MOISTURE BUDGET
!        (4) DO CLOUD ASCENT IN 'CUASC' IN ABSENCE OF DOWNDRAFTS
!        (5) DO DOWNDRAFT CALCULATIONS:
!              (A) DETERMINE VALUES AT LFS IN 'CUDLFS'
!              (B) DETERMINE MOIST DESCENT IN 'CUDDRAF'
!              (C) RECALCULATE CLOUD BASE MASSFLUX CONSIDERING THE
!                  EFFECT OF CU-DOWNDRAFTS
!        (6) DO FINAL CLOUD ASCENT IN 'CUASC'
!        (7) DO FINAL ADJUSMENTS TO CONVECTIVE FLUXES IN 'CUFLX',
!            DO EVAPORATION IN SUBCLOUD LAYER
!        (8) CALCULATE INCREMENTS OF T AND Q IN 'CUDTDQ'
!        (9) CALCULATE INCREMENTS OF U AND V IN 'CUDUDV'
!
!     EXTERNALS.
!     ----------
!
!       CUINI:  INITIALIZES VALUES AT VERTICAL GRID USED IN CU-PARAMETR.
!       CUBASE: CLOUD BASE CALCULATION FOR PENETR.AND SHALLOW CONVECTION
!       CUASCT: CLOUD ASCENT FOR ENTRAINING PLUME
!       CUDLFS: DETERMINES VALUES AT LFS FOR DOWNDRAFTS
!       CUDDRAF:DOES MOIST DESCENT FOR CUMULUS DOWNDRAFTS
!       CUFLX:  FINAL ADJUSTMENTS TO CONVECTIVE FLUXES (ALSO IN PBL)
!       CUDQDT: UPDATES TENDENCIES FOR T AND Q
!       CUDUDV: UPDATES TENDENCIES FOR U AND V
!
!     SWITCHES.
!     --------
!
!          LMFPEN=.T.   PENETRATIVE CONVECTION IS SWITCHED ON
!          LMFSCV=.T.   SHALLOW CONVECTION IS SWITCHED ON
!          LMFMID=.T.   MIDLEVEL CONVECTION IS SWITCHED ON
!          LMFDD=.T.    CUMULUS DOWNDRAFTS SWITCHED ON
!          LMFDUDV=.T.  CUMULUS FRICTION SWITCHED ON
!
!     MODEL PARAMETERS (DEFINED IN MODULE mo_cumulus_flux)
!     ----------------------------------------------------
!     ENTRPEN    ENTRAINMENT RATE FOR PENETRATIVE CONVECTION
!     ENTRSCV    ENTRAINMENT RATE FOR SHALLOW CONVECTION
!     ENTRMID    ENTRAINMENT RATE FOR MIDLEVEL CONVECTION
!     ENTRDD     ENTRAINMENT RATE FOR CUMULUS DOWNDRAFTS
!     CMFCTOP    RELATIVE CLOUD MASSFLUX AT LEVEL ABOVE NONBUOYANCY LEVEL
!     CMFCMAX    MAXIMUM MASSFLUX VALUE ALLOWED FOR
!     CMFCMIN    MINIMUM MASSFLUX VALUE (FOR SAFETY)
!     CMFDEPS    FRACTIONAL MASSFLUX FOR DOWNDRAFTS AT LFS
!     CPRCON     COEFFICIENT FOR CONVERSION FROM CLOUD WATER TO RAIN
!
!     REFERENCE.
!     ----------
!
!          PAPER ON MASSFLUX SCHEME (TIEDTKE,1989)
!
!reUSE mo_kind,           ONLY: dp
!reUSE mo_control,        ONLY: nn
!reUSE mo_constants,      ONLY: g, alv, als, tmelt, vtmpc1, rd
!reUSE mo_cumulus_flux,   ONLY: entrpen, entrscv, lmfdd, cmfdeps, lmfdudv
!reUSE mo_convect_tables, ONLY: tlucua,                         & ! table a
!re                             tlucub,                         & ! table b
!re                             jptlucu1, jptlucu2,                       &
!re                             lookuperror, lookupoverflow
!reUSE mo_time_control,   ONLY: time_step_len
!

   USE messy_main_constants_mem,     ONLY: g, alv, als, tmelt, vtmpc1, rd

   USE MESSY_CONVECT_TIEDTKE_PARAM,  ONLY: entrpen, entrscv, lmfdd, cmfdeps, lmfdudv

   USE messy_main_tools,             ONLY: tlucua, tlucub, jptlucu1, jptlucu2


!
IMPLICIT NONE
!
!
INTEGER, INTENT (IN) :: nn 
REAL(DP),INTENT(IN):: time_step_len!re 
LOGICAL, INTENT(INOUT):: lookupoverflow!re 

INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, klevm1!re, ktrac

REAL(DP), INTENT(IN):: delta_time!re

!
!---wiso-code

  INTEGER, INTENT (IN) :: kwiso

!---wiso-code-end

REAL(dp):: pten(kbdim,klev),        pqen(kbdim,klev),                  &
           pxen(kbdim,klev),        ptven(kbdim,klev),                 &
           puen(kbdim,klev),        pven(kbdim,klev),                  &
           ptte(kbdim,klev),        pqte(kbdim,klev),                  &
           pvom(kbdim,klev),        pvol(kbdim,klev),                  &
           pqsen(kbdim,klev),       pgeo(kbdim,klev),                  &
           paphp1(kbdim,klevp1),                                       &
           pverv(kbdim,klev)
REAL(dp):: ptu(kbdim,klev),         pqu(kbdim,klev),                   &
           plu(kbdim,klev),         plude(kbdim,klev),                 &
           pmfu(kbdim,klev),        pmfd(kbdim,klev),                  &
           paprc(kbdim),            paprs(kbdim),                      &
           prsfc(kbdim),            pssfc(kbdim),                      &
           prain(kbdim),            pqhfla(kbdim)
INTEGER :: kcbot(kbdim),            kctop(kbdim),                      &
           ktype(kbdim)
REAL(dp):: pxtec(kbdim,klev),       pqtec(kbdim,klev),                 &
           pqude(kbdim,klev)

!---wiso-code

REAL(dp):: pwisoqen(kbdim,klev,kwiso),                                 &
           pwisoxen(kbdim,klev,kwiso),                                 &
           pwisoqte(kbdim,klev,kwiso)
REAL(dp):: pwisoqu(kbdim,klev,kwiso),                                  &
           pwisolu(kbdim,klev,kwiso), pwisolude(kbdim,klev,kwiso),     &
           pwisoaprc(kbdim,kwiso),    pwisoaprs(kbdim,kwiso),          &
           pwisorsfc(kbdim,kwiso),    pwisossfc(kbdim,kwiso),          &
           pwisorain(kbdim,kwiso)
REAL(dp):: pwisoxtec(kbdim,klev,kwiso),pwisoqtec(kbdim,klev,kwiso),    &
           pwisoqude(kbdim,klev,kwiso)

!---wiso-code-end

REAL(dp):: ztenh(kbdim,klev),       zqenh(kbdim,klev),                 &
           zxenh(kbdim,klev),                                          &
           zgeoh(kbdim,klev),       zqsenh(kbdim,klev),                &
           ztd(kbdim,klev),         zqd(kbdim,klev),                   &
           zmfus(kbdim,klev),       zmfds(kbdim,klev),                 &
           zmfuq(kbdim,klev),       zmfdq(kbdim,klev),                 &
           zdmfup(kbdim,klev),      zdmfdp(kbdim,klev),                &
           zmful(kbdim,klev),       zrfl(kbdim),                       &
           zuu(kbdim,klev),         zvu(kbdim,klev),                   &
           zud(kbdim,klev),         zvd(kbdim,klev)
REAL(dp):: zcpen(kbdim,klev),       zcpcu(kbdim,klev)
REAL(dp):: zentr(kbdim),            zhcbase(kbdim),                    &
           zmfub(kbdim),            zmfub1(kbdim),                     &
           zdqpbl(kbdim),           zdqcv(kbdim)
REAL(dp):: zsfl(kbdim),             zdpmel(kbdim,klev)
REAL(dp):: zcape(kbdim),            zheat(kbdim)
INTEGER :: ilab(kbdim,klev),        idtop(kbdim),                      &
           ictop0(kbdim),           ilwmin(kbdim)
!re REAL(dp):: pxten(kbdim,klev,ktrac), pxtte(kbdim,klev,ktrac),           &
!re            pxtu(kbdim,klev,ktrac),                                     &
!re            zxtenh(kbdim,klev,ktrac),zxtd(kbdim,klev,ktrac),            &
!re            zmfuxt(kbdim,klev,ktrac),zmfdxt(kbdim,klev,ktrac)
LOGICAL :: loddraf(kbdim),          ldland(kbdim)
LOGICAL :: ldcum(kbdim)
LOGICAL :: llo1, lo
!
INTEGER :: jl, jk, ikb, it, it1, jt, itopm2
REAL(dp):: zcons2, ztau, zqumqe, zdqmin, zmfmax, zalvs, zalvdcp, zqalv &
         , zhsat, zes, zcor, zqsat, zqst1, zdqsdt, zgam, zzz, zhhat    &
         , zro, zdz, zfac, zpbmpt, zeps

!---wiso-code

REAL(dp):: zwisoqenh(kbdim,klev,kwiso),                                &
           zwisoxenh(kbdim,klev,kwiso),                                &
           zwisoqd(kbdim,klev,kwiso),                                  &
           zwisomfuq(kbdim,klev,kwiso),zwisomfdq(kbdim,klev,kwiso),    &
           zwisodmfup(kbdim,klev,kwiso),zwisodmfdp(kbdim,klev,kwiso),  &
           zwisomful(kbdim,klev,kwiso)

REAL(dp):: zmelt_tmp(kbdim,klev),        &  ! temporary value
           zprec_frac(kbdim,klev),       &  ! temporary value
           zrain_tmp(kbdim,klev),        &  ! temporary value 
           zrfl_tmp(kbdim),              &  ! temporary value
           zsfl_tmp(kbdim),              &  ! temporary value
           zqsto(kbdim,klev),            &  ! old moisture value
           ztsto(kbdim,klev),            &  ! old temperature value
           zwisoqsto(kbdim,klev,kwiso),  &  ! old water isotope moisture value
           zwisorfl(kbdim,kwiso),        &  ! temporary value
           zwisosfl(kbdim,kwiso)            ! temporary value

!---wiso-code-end
!
!  INTRINSIC FUNCTIONS
INTRINSIC MIN, MAX
!
!  Executable statements

  lookupoverflow = .FALSE.
!-----------------------------------------------------------------------
!
!     1.           SPECIFY CONSTANTS AND PARAMETERS
!                  --------------------------------
!
100 CONTINUE
!
  zcons2=1._dp/(g*time_step_len)
  ztau=MIN(3._dp*3600._dp,7200._dp*63._dp/nn)

!---wiso-code

  zprec_frac(:,:)=0._dp
  zrain_tmp(:,:)=0._dp

!---wiso-code-end
!
!-----------------------------------------------------------------------
!
!*    2.           INITIALIZE VALUES AT VERTICAL GRID POINTS IN 'CUINI'
!                  ---------------------------------------------------
!
200 CONTINUE
  CALL h2oiso_cuini(kproma,   kbdim,    klev,     klevp1,   klevm1,           &
             pten,     pqen,     pqsen,    pxen,     puen,     pven,   &
             ptven, &!re,    ktrac,                                          &
!re             pxten,    zxtenh,   pxtu,     zxtd,     zmfuxt,   zmfdxt, &
             pverv,    pgeo,     paphp1,   zgeoh,                      &
             ztenh,    zqenh,    zqsenh,   zxenh,    ilwmin,           &
             ptu,      pqu,      ztd,      zqd,                        &
             zuu,      zvu,      zud,      zvd,                        &
             pmfu,     pmfd,     zmfus,    zmfds,                      &
             zmfuq,    zmfdq,    zdmfup,   zdmfdp,                     &
             zcpen,    zcpcu,                                          &
             zdpmel,   plu,      plude,    pqude,    ilab,             &
!---wiso-code
             kwiso,                                                    &
             pwisoqen, pwisoxen,                                       &
             zwisoqenh,zwisoxenh,                                      &
             pwisoqu,  zwisoqd,                                        &
             zwisomfuq,zwisomfdq,zwisodmfup,zwisodmfdp,                &
             pwisolu,  pwisolude,pwisoqude                             &
             )         
!---wiso-code-end
    IF (lookupoverflow) RETURN!re
!
!-----------------------------------------------------------------------
!
!*    3.0          CLOUD BASE CALCULATIONS
!                  -----------------------
!
300 CONTINUE
!
!*             (A) DETERMINE CLOUD BASE VALUES IN 'CUBASE'
!                  ---------------------------------------
!
  CALL h2oiso_cubase(kproma,   kbdim,    klev,     klevp1,   klevm1,          &
              ztenh,    zqenh,    zgeoh,    paphp1,                    &
              ptu,      pqu,      plu,                                 &
              puen,     pven,     zuu,      zvu,                       &
              zcpcu,                                                   &
              ldcum,    kcbot,    ilab,                                &
!---wiso-code
              kwiso,                                                   &
              pwisoqu,  pwisolu                                        &
              )
!---wiso-code-end
    IF (lookupoverflow) RETURN!re
!
!*             (B) DETERMINE TOTAL MOISTURE CONVERGENCE AND
!*                 THEN DECIDE ON TYPE OF CUMULUS CONVECTION
!                  -----------------------------------------
!
  jk=1
  DO 310 jl=1,kproma
     zdqpbl(jl)=0.0_dp
     zdqcv(jl)=pqte(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))
     idtop(jl)=0
310 END DO
  DO 320 jk=2,klev
     DO 315 jl=1,kproma
        zdqcv(jl)=zdqcv(jl)+pqte(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))
        IF(jk.GE.kcbot(jl)) zdqpbl(jl)=zdqpbl(jl)+pqte(jl,jk)          &
                                       *(paphp1(jl,jk+1)-paphp1(jl,jk))
315  END DO
320 END DO
!
!*             (C) DETERMINE MOISTURE SUPPLY FOR BOUNDARY LAYER
!*                 AND DETERMINE CLOUD BASE MASSFLUX IGNORING
!*                 THE EFFECTS OF DOWNDRAFTS AT THIS STAGE
!                  ------------------------------------------
!
!DIR$ IVDEP
  DO 340 jl=1,kproma
     ikb=kcbot(jl)
     zqumqe=pqu(jl,ikb)+plu(jl,ikb)-zqenh(jl,ikb)
     zdqmin=MAX(0.01_dp*zqenh(jl,ikb),1.e-10_dp)
     llo1=zdqpbl(jl).GT.0._dp.AND.zqumqe.GT.zdqmin.AND.ldcum(jl)
     zmfub(jl)=MERGE(zdqpbl(jl)/(g*MAX(zqumqe,zdqmin)),0.01_dp,llo1)
     zmfmax=(paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
     zmfub(jl)=MIN(zmfub(jl),zmfmax)
     IF(.NOT.llo1) ldcum(jl)=.FALSE.
     ktype(jl)=MERGE(1,2,zdqcv(jl).GT.MAX(0._dp,-1.1_dp*pqhfla(jl)*g))
     zentr(jl)=MERGE(entrpen,entrscv,ktype(jl).EQ.1)
340 END DO
!
!-----------------------------------------------------------------------
!*    4.0          DETERMINE CLOUD ASCENT FOR ENTRAINING PLUME
!                  -------------------------------------------
!
400 CONTINUE
!
!*             (A) ESTIMATE CLOUD HEIGHT FOR ENTRAINMENT/DETRAINMENT
!*                 CALCULATIONS IN CUASC (MAX.POSSIBLE CLOUD HEIGHT
!*                 FOR NON-ENTRAINING PLUME, FOLLOWING A.-S.,1974)
!                  -------------------------------------------------
!
!DIR$ IVDEP
  DO 410 jl=1,kproma
     ikb=kcbot(jl)
     zalvs=MERGE(alv,als,ptu(jl,ikb)>tmelt)
     zhcbase(jl)=zcpcu(jl,ikb)*ptu(jl,ikb)+zgeoh(jl,ikb)+              &
                                                      zalvs*pqu(jl,ikb)
     ictop0(jl)=kcbot(jl)-1
410 END DO
  DO 430 jk=klevm1,3,-1
!DIR$ IVDEP
     DO 420 jl=1,kproma
        zalvs=MERGE(alv,als,ztenh(jl,jk)>tmelt)
        zalvdcp=zalvs/zcpcu(jl,jk)
        zqalv=1._dp/zalvs
        zhsat=zcpcu(jl,jk)*ztenh(jl,jk)+zgeoh(jl,jk)+zalvs*zqsenh(jl,jk)
        it = NINT(ztenh(jl,jk)*1000._dp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zes=tlucua(it)/paphp1(jl,jk)
        zes=MIN(0.5_dp,zes)
        LO=zes<0.4_dp
        zcor=1._dp/(1._dp-vtmpc1*zes)
        zqsat=zes*zcor
        it1=it+1
        it1 = MAX(MIN(it1,jptlucu2),jptlucu1)
        zqst1=tlucua(it1)/paphp1(jl,jk)
        zqst1=MIN(0.5_dp,zqst1)
        zqst1=zqst1/(1._dp-vtmpc1*zqst1)
        zdqsdt=(zqst1-zqsat)*1000._dp
        zgam=MERGE(zalvdcp*zdqsdt,zqsat*zcor*tlucub(it),LO)
        zzz=zcpcu(jl,jk)*ztenh(jl,jk)*vtmpc1
        zhhat=zhsat-(zzz+zgam*zzz)/(1._dp+zgam*zzz*zqalv)*             &
                          MAX(zqsenh(jl,jk)-zqenh(jl,jk),0._dp)
        IF(jk.LT.ictop0(jl).AND.zhcbase(jl).GT.zhhat) ictop0(jl)=jk
420  END DO
430 END DO
!
     IF (lookupoverflow) RETURN !CALL lookuperror ('cumastrh')
!!
!!     DEEP CONVECTION IF CLOUD DEPTH > 200 HPA, ELSE SHALLOW
!!     (CLOUD DEPTH FROM NON-ENTRAINIG PLUME)
!!
!  DO jl=1,kproma
!     ktype(jl)=MERGE(1,2,                                               &
!                paphp1(jl,kcbot(jl))-paphp1(jl,ictop0(jl)).gt.2.E4_dp)
!     zentr(jl)=MERGE(entrpen,entrscv,ktype(jl).eq.1)
!  ENDDO
!!
!*             (B) DO ASCENT IN 'CUASCT' IN ABSENCE OF DOWNDRAFTS
!                  ----------------------------------------------
!
  CALL h2oiso_cuasct(kproma,  kbdim,    klev,     klevp1,   klevm1,           &
             ztenh,    zqenh,    puen,     pven,                       &
!re             ktrac,                                                    &
!re             zxtenh,   pxten,    pxtu,     zmfuxt,                     &
             pten,     pqen,     pqsen,                                &
             pgeo,     zgeoh,    paphp1,                               &
             pqte,     pverv,    ilwmin,                               &
             ldcum,    ldland,   ktype,    ilab,                       &
             ptu,      pqu,      plu,      zuu,      zvu,              &
             pmfu,     zmfub,    zentr,                                &
             zmfus,    zmfuq,                                          &
             zmful,    plude,    pqude,    zdmfup,                     &
             zcpen,    zcpcu,                                          &
             kcbot,    kctop,    ictop0,                               &
!---wiso-code
             kwiso,                                                    &
             zwisoqenh,                                                &
             pwisoqen,                                                 &
             pwisoqu,  pwisolu,                                        &
             zwisomfuq,                                                &
             zwisomful,pwisolude,pwisoqude,zwisodmfup                  &
             , nn, time_step_len)
!---wiso-code-end
!
!*         (C) CHECK CLOUD DEPTH AND CHANGE ENTRAINMENT RATE ACCORDINGLY
!              CALCULATE PRECIPITATION RATE (FOR DOWNDRAFT CALCULATION)
!              ---------------------------------------------------------
!
!DIR$ IVDEP
  DO 440 jl=1,kproma
     zpbmpt=paphp1(jl,kcbot(jl))-paphp1(jl,kctop(jl))
     IF(ldcum(jl).AND.ktype(jl).EQ.1.AND.zpbmpt.LT.2.e4_dp) ktype(jl)=2
     IF(ldcum(jl)) ictop0(jl)=kctop(jl)
     IF(ktype(jl).EQ.2) zentr(jl)=entrscv
     zrfl(jl)=zdmfup(jl,1)
440 END DO
  DO 460 jk=2,klev
     DO 450 jl=1,kproma
        zrfl(jl)=zrfl(jl)+zdmfup(jl,jk)
450  END DO
460 END DO
!
!-----------------------------------------------------------------------
!*    5.0          CUMULUS DOWNDRAFT CALCULATIONS
!                  ------------------------------
!
500 CONTINUE
!
  IF(lmfdd) THEN
!
!*             (A) DETERMINE LFS IN 'CUDLFS'
!                  -------------------------
!
     CALL h2oiso_cudlfs(kproma,   kbdim,    klev,     klevp1,                 &
                 ztenh,    zqenh,    puen,     pven,                   &
!re                 ktrac,                                                &
!re                 zxtenh,   pxtu,     zxtd,     zmfdxt,                 &
                 zgeoh,    paphp1,                                     &
                 ptu,      pqu,      zuu,      zvu,                    &
                 ldcum,    kcbot,    kctop,    zmfub,    zrfl,         &
                 ztd,      zqd,      zud,      zvd,                    &
                 pmfd,     zmfds,    zmfdq,    zdmfdp,                 &
                 zcpcu,                                                &
                 idtop,    loddraf,                                    &
!---wiso-code
                 kwiso,                                                &
                 zwisoqenh,                                            &
                 pwisoqu,                                              &
                 zwisoqd,                                              &
                 zwisomfdq,zwisodmfdp,                                 &
                 zdmfup,   zwisodmfup                                  &
                 )
!---wiso-code-end
!
!*            (B)  DETERMINE DOWNDRAFT T,Q AND FLUXES IN 'CUDDRAF'
!                  -----------------------------------------------
!
     CALL h2oiso_cuddraf(kproma,   kbdim,    klev,     klevp1,                &
                  ztenh,    zqenh,    puen,     pven,                  &
!re                  ktrac,                                               &
!re                  zxtenh,   zxtd,     zmfdxt,                          &
                  zgeoh,    paphp1,   zrfl,                            &
                  ztd,      zqd,      zud,      zvd,                   &
                  pmfd,     zmfds,    zmfdq,    zdmfdp,                &
                  zcpcu,                                               &
                  loddraf,                                             &
!---wiso-code
                  kwiso,                                               &
                  zwisoqenh,                                           &
                  zwisoqd,                                             &
                  zwisomfdq,zwisodmfdp,                                &
                  zdmfup,   zwisodmfup                                 &
                  )
!---wiso-code-end
!
  END IF
!
!*    5.5          RECALCULATE CLOUD BASE MASSFLUX FROM A
!*                 CAPE CLOSURE FOR DEEP CONVECTION (KTYPE=1)
!*                 AND BY PBL EQUILIBRUM TAKING DOWNDRAFTS INTO
!*                 ACCOUNT FOR SHALLOW CONVECTION (KTYPE=2)
!                  -------------------------------------------
!
  DO jl=1,kproma
     zheat(jl)=0._dp
     zcape(jl)=0._dp
     zmfub1(jl)=zmfub(jl)
  ENDDO
!
  DO jk=1,klev
     DO jl=1,kproma
        llo1=ldcum(jl).AND.ktype(jl).EQ.1
        IF(llo1.AND.jk.LE.kcbot(jl).AND.jk.GT.kctop(jl)) THEN
           ikb=kcbot(jl)
           zro=paphp1(jl,jk)/(rd*ztenh(jl,jk)*                         &
                                           (1._dp+vtmpc1*zqenh(jl,jk)))
           zdz=(paphp1(jl,jk)-paphp1(jl,jk-1))/(g*zro)
           zheat(jl)=zheat(jl) +                                       &
                (  (pten(jl,jk-1)-pten(jl,jk) + g*zdz/zcpcu(jl,jk))    &
                     /ztenh(jl,jk)                                     &
                    +  vtmpc1*(pqen(jl,jk-1)-pqen(jl,jk))  ) *         &
                       (g*(pmfu(jl,jk)+pmfd(jl,jk)))/zro
           zcape(jl)=zcape(jl) +                                       &
                         (g*(ptu(jl,jk)-ztenh(jl,jk))/ztenh(jl,jk)     &
                              +g*vtmpc1*(pqu(jl,jk)-zqenh(jl,jk))      &
                              -g*plu(jl,jk) ) * zdz
        ENDIF
     ENDDO
  ENDDO
!
  DO jl=1,kproma
     IF(ldcum(jl).AND.ktype(jl).EQ.1) THEN
        ikb=kcbot(jl)
        zmfub1(jl)=(zcape(jl)*zmfub(jl))/(zheat(jl)*ztau)
        zmfub1(jl)=MAX(zmfub1(jl),0.001_dp)
        zmfmax=(paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
        zmfub1(jl)=MIN(zmfub1(jl),zmfmax)
     ENDIF
  ENDDO
!
!*                 RECALCULATE CONVECTIVE FLUXES DUE TO EFFECT OF
!*                 DOWNDRAFTS ON BOUNDARY LAYER MOISTURE BUDGET
!*                 FOR SHALLOW CONVECTION (KTYPE=2)
!                  --------------------------------------------
!
!DIR$ IVDEP
  DO 520 jl=1,kproma
     IF(ktype(jl).EQ.2) THEN
        ikb=kcbot(jl)
        llo1=pmfd(jl,ikb).LT.0._dp.AND.loddraf(jl)
        zeps=MERGE(cmfdeps,0._dp,llo1)
        zqumqe=pqu(jl,ikb)+plu(jl,ikb)-                                &
                    zeps*zqd(jl,ikb)-(1._dp-zeps)*zqenh(jl,ikb)
        zdqmin=MAX(0.01_dp*zqenh(jl,ikb),1.e-10_dp)
        zmfmax=(paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
        llo1=zdqpbl(jl).GT.0._dp.AND.zqumqe.GT.zdqmin.AND.ldcum(jl)    &
                    .AND.zmfub(jl).LT.zmfmax
        zmfub1(jl)=MERGE(zdqpbl(jl)/(g*MAX(zqumqe,zdqmin)),            &
                               zmfub(jl),llo1)
        zmfub1(jl)=MERGE(zmfub1(jl),zmfub(jl),                         &
                         ABS(zmfub1(jl)-zmfub(jl)).LT.0.2_dp*zmfub(jl))
     END IF
520 END DO
  DO 540 jk=1,klev
     DO 530 jl=1,kproma
        IF(ldcum(jl)) THEN
           zfac=zmfub1(jl)/MAX(zmfub(jl),1.e-10_dp)
           pmfd(jl,jk)=pmfd(jl,jk)*zfac
           zmfds(jl,jk)=zmfds(jl,jk)*zfac
           zmfdq(jl,jk)=zmfdq(jl,jk)*zfac
           zdmfdp(jl,jk)=zdmfdp(jl,jk)*zfac
        END IF
530  END DO
!
!---wiso-code

! Change water isotope fluxes after downdraft calculation and new cloud base mass flux
     DO jt=1,kwiso
       DO jl=1,kproma
         IF (ldcum(jl)) THEN
           zfac = zmfub1(jl)/MAX(zmfub(jl),1.E-10_dp)
           zwisomfdq(jl,jk,jt) = zwisomfdq(jl,jk,jt)*zfac
           zwisodmfdp(jl,jk,jt) = zwisodmfdp(jl,jk,jt)*zfac
         END IF
       END DO
     END DO

!---wiso-code-end

!re     DO 5304 jt=1,ktrac
!re        DO 5302 jl=1,kproma
!re           IF(ldcum(jl)) THEN
!re              zfac=zmfub1(jl)/MAX(zmfub(jl),1.e-10_dp)
!re              zmfdxt(jl,jk,jt)=zmfdxt(jl,jk,jt)*zfac
!re           ENDIF
!re5302    END DO
!re5304 END DO
!
540 END DO
!
!*                 NEW VALUES OF CLOUD BASE MASS FLUX
!                  ----------------------------------
!
  DO 550 jl=1,kproma
     IF(ldcum(jl)) zmfub(jl)=zmfub1(jl)
550 END DO
!
!-----------------------------------------------------------------------
!*    6.0          DETERMINE FINAL CLOUD ASCENT FOR ENTRAINING PLUME
!*                 FOR PENETRATIVE CONVECTION (TYPE=1),
!*                 FOR SHALLOW TO MEDIUM CONVECTION (TYPE=2)
!*                 AND FOR MID-LEVEL CONVECTION (TYPE=3).
!                  --------------------------------------------------
!
600 CONTINUE
  CALL h2oiso_cuasct(kproma,  kbdim,    klev,     klevp1,   klevm1,           &
             ztenh,    zqenh,    puen,     pven,                       &
!re             ktrac,                                                    &
!re             zxtenh,   pxten,    pxtu,     zmfuxt,                     &
             pten,     pqen,     pqsen,                                &
             pgeo,     zgeoh,    paphp1,                               &
             pqte,     pverv,    ilwmin,                               &
             ldcum,    ldland,   ktype,    ilab,                       &
             ptu,      pqu,      plu,      zuu,      zvu,              &
             pmfu,     zmfub,    zentr,                                &
             zmfus,    zmfuq,                                          &
             zmful,    plude,    pqude,    zdmfup,                     &
             zcpen,    zcpcu,                                          &
             kcbot,    kctop,    ictop0,                               &
!---wiso-code
             kwiso,                                                    &
             zwisoqenh,                                                &
             pwisoqen,                                                 &
             pwisoqu,  pwisolu,                                        &
             zwisomfuq,                                                &
             zwisomful,pwisolude,pwisoqude,zwisodmfup                  &
             , nn, time_step_len)
!---wiso-code-end
!
!-----------------------------------------------------------------------
!*    7.0          DETERMINE FINAL CONVECTIVE FLUXES IN 'CUFLX'
!                  ------------------------------------------
!
700 CONTINUE
  CALL h2oiso_cuflx(kproma,   kbdim,    klev,     klevp1,                     &
             pqen,     pqsen,    ztenh,    zqenh,                      &
!re             ktrac,                                                    &
!re             zxtenh,   zmfuxt,   zmfdxt,                               &
             paphp1,   zgeoh,                                          &
             kcbot,    kctop,    idtop,                                &
             ktype,    loddraf,  ldcum,                                &
             pmfu,     pmfd,     zmfus,    zmfds,                      &
             zmfuq,    zmfdq,    zmful,                                &
             zdmfup,   zdmfdp,   zrfl,     prain,                      &
             zcpcu,                                                    &
             pten,     zsfl,     zdpmel,   itopm2,                     &
!---wiso-code
             kwiso,                                                    &
             zwisoqenh,                                                &
             zwisomfuq,zwisomfdq,zwisomful,                            &
             zwisodmfdp,                                               &
             zmelt_tmp,zprec_frac,zrain_tmp,                           &
             zrfl_tmp, zsfl_tmp                                        &
             , time_step_len)
!---wiso-code-end
!
!-----------------------------------------------------------------------
!
!*    8.0          UPDATE TENDENCIES FOR T AND Q IN SUBROUTINE CUDTDQ
!                  --------------------------------------------------
!
!---wiso-code

!  Calculate old values (T-1)
  DO jk=1,klev
    DO jl=1,kproma
      zqsto(jl,jk)=pqen(jl,jk)-pqte(jl,jk)*time_step_len
      ztsto(jl,jk)=pten(jl,jk)-ptte(jl,jk)*time_step_len
    END DO
  END DO
  DO jt=1,kwiso
    DO jk=1,klev
      DO jl=1,kproma
        zwisoqsto(jl,jk,jt)=pwisoqen(jl,jk,jt)-pwisoqte(jl,jk,jt)*time_step_len
      END DO
    END DO
  END DO

!---wiso-code-end

800 CONTINUE
  CALL h2oiso_cudtdq(kproma, kbdim, klev, klevp1, itopm2, ldcum, & !re, ktrac,       &
              paphp1,   pten,     ptte,     pqte,                      &
              pxtec, &!re ,    pxtte,    zmfuxt,   zmfdxt,                    &
              zmfus,    zmfds,    zmfuq,    zmfdq,                     &
              zmful,    zdmfup,   zdmfdp,   plude,                     &
              zdpmel,   zrfl,     zsfl,                                &
              zcpen,    pqtec,    pqude,                               &
              prsfc,    pssfc,    paprc,    paprs, delta_time)
!
!---wiso-code

!  update water isotope tendencies
  CALL h2oiso_cuwisodq(kproma, kbdim, klev, klevp1, itopm2, ldcum, kwiso,     &
                paphp1,   pwisoqte,                                    &
                pwisoxtec,                                             &
                zwisomfuq,zwisomfdq,                                   &
                zwisomful,zwisodmfup,zwisodmfdp,pwisolude,             &
                pwisoqtec,pwisoqude, delta_time)

! calculate precipitation of water isotopes and
! get precipitation into equilibrium with surrounding vapour 
  CALL h2oiso_cuwisoequ(kproma, kbdim, klev, klevp1, itopm2, ldcum, kwiso,    &
                         kctop,          kcbot,       paphp1,                  &
                         ztsto,          ptte,                                 & 
                         zqsto,          pqte,                                 &
                         zwisoqsto,      pwisoqte,                             &
                         ptu,            pqu,         pwisoqu,                 &
                         pten,                                                 &
                                 zrfl_tmp,       zsfl_tmp,                             &
                         zwisorfl,       zwisosfl,                             &
                         zwisodmfup,     zwisodmfdp,                           &
                         zmelt_tmp,                                            &
                         zprec_frac,     zrain_tmp,                            &
                         pwisorsfc,      pwisossfc,                            &
                         pwisoaprc,      pwisoaprs,                            &
                                 time_step_len, delta_time )
                     
!---wiso-code-end

!-----------------------------------------------------------------------
!
!*    9.0          UPDATE TENDENCIES FOR U AND U IN SUBROUTINE CUDUDV
!                  --------------------------------------------------
!
!re 900 CONTINUE
!re   IF(lmfdudv) THEN
!re      CALL h2oiso_cududv(kproma,   kbdim,    klev,     klevp1,                 &
!re                  itopm2,   ktype,    kcbot,    paphp1,   ldcum,        &
!re                  puen,     pven,     pvom,     pvol,                   &
!re                  zuu,      zud,      zvu,      zvd,                    &
!re                  pmfu,     pmfd)
!re !
!re   END IF
!re !
!re 1000 CONTINUE
!
  RETURN
END SUBROUTINE h2oiso_cumastrh




! **********************************************************************
END MODULE messy_h2oiso_convect
! **********************************************************************
