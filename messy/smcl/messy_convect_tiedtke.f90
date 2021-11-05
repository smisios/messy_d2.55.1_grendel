MODULE MESSY_CONVECT_TIEDTKE

! This module contains the subroutines for the Tiedtke Convection
! Scheme. They are called from the cumastr/cumastrh/cumastrt from 
! MESSY_convection or other subroutines within this module

! Author:  H.Tost,   MPICH, March 2004

! changes are mostly to eliminate USE-Statements and 
! defining undefined variables and adding IMPLICIT NONE


! original code from TIEDTKE (see below)
! op_ck_20031001+
!     UPDATE:
!     -------
!     
!          Michael Ponater, DLR OP, Aug. 2003
!          Update from Sabine Brinkop to ensure positive tracer concentration
!          during the model integration. 
!          Update in subroutines:  CUINI, CUASC, CUFLX, CUDLFS, CUDDRAF, 
!          CUBASMC, CUMASTR, CUDTDQ, VDIFF
!          Look for string "!pa31, SB" to identify the modifications
!
!          Reference: Brinkop, S., Sausen, R.: A Finite Difference 
!           Approximation for Convective Transports which Maintains Positive 
!           Tracer Concentrations. 
!           Beiträge zur Physik der Atmosphäre, 3, (1997), S. 245-248, 
! op_ck_20031001-

  USE messy_main_constants_mem,     ONLY: g, R_gas, M_air, T0, wp, &
       rd, rv, cpd=>cp_air, cpv, vtmpc1, vtmpc2, tmelt
  USE messy_convect_tiedtke_param
  USE messy_convect_mem

  USE messy_main_tools,             ONLY: jptlucu1, jptlucu2, tlucua, tlucub, tlucuc   
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
  REAL(dp), ALLOCATABLE, DIMENSION(:):: cevapcu

  LOGICAL :: lookupoverflow = .FALSE.

!======================================================================
CONTAINS
!======================================================================

SUBROUTINE tiedtke_cuadjtq(kproma,       kbdim,    klev,    & 
           kk,                                              & 
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


  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: kcall, kk, klev, kproma, kbdim

  !  Array arguments with intent(In):
  REAL(dp), INTENT (IN) :: pp(kbdim)
  LOGICAL, INTENT (IN) :: ldflag(kbdim)
 
  !  Array arguments with intent(InOut):
  REAL(dp), INTENT (INOUT) :: pq(kbdim,klev), pt(kbdim,klev)



  !  Local scalars: 
  REAL(dp)    :: zcond1, zqst1, zdqsdt, zqsat, zes, zcor, zlcdqsdt
  INTEGER :: isum, jl, it, it1
  LOGICAL :: LO

  !  Local arrays: 
  REAL(dp) ::    zcond(kbdim)


  !  Executable statements 

  lookupoverflow = .FALSE.
  zcond(:)=0._dp !mz_ht_20040121
!
!----------------------------------------------------------------------
!
!     2.           CALCULATE CONDENSATION AND ADJUST T AND Q ACCORDINGLY
!                  -----------------------------------------------------
!
!!200 CONTINUE
  IF (kcall.EQ.1 ) THEN

     isum=0
!DIR$ IVDEP
!OCL NOVREC
     DO 210 jl=1,kproma
        IF(ldflag(jl)) THEN
           it = NINT(pt(jl,kk)*1000.)
! ka_sv_20170410+
! limiting to below thermosphere
!!$        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
           IF ((it<jptlucu1 .OR. it>jptlucu2) .AND. &
                (pp(jl) >= 1.0_dp)) lookupoverflow = .TRUE.
! ka_sv_20170410-
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
           IF(ABS(zcond(jl)) > 0.) isum=isum+1
        END IF
210  END DO

     IF (lookupoverflow) RETURN

!     IF(isum.EQ.0) go to 230 ! mz_ht_20071217

!DIR$ IVDEP
!OCL NOVREC
     DO 220 jl=1,kproma
        IF(ldflag(jl)) THEN
           IF(ABS(zcond(jl)) > 0._dp) THEN
           it = NINT(pt(jl,kk)*1000.)
! ka_sv_20170410+
! limiting to below thermosphere
!!$        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
           IF ((it<jptlucu1 .OR. it>jptlucu2) .AND. &
                (pp(jl) >= 1.0_dp)) lookupoverflow = .TRUE.
! ka_sv_20170410-
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

     IF (lookupoverflow) RETURN

230  CONTINUE

  END IF

  IF(kcall.EQ.2) THEN

     isum=0
!DIR$ IVDEP
!OCL NOVREC
     DO 310 jl=1,kproma
        IF(ldflag(jl)) THEN
           it = NINT(pt(jl,kk)*1000._dp)
! ka_sv_20170410+
! limiting to below thermosphere
!!$        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
           IF ((it<jptlucu1 .OR. it>jptlucu2) .AND. &
                (pp(jl) >= 1.0_dp)) lookupoverflow = .TRUE.
! ka_sv_20170410-
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

     IF (lookupoverflow) RETURN

!     IF(isum.EQ.0) go to 330 ! mz_ht_20071217

!DIR$ IVDEP
!OCL NOVREC
     DO 320 jl=1,kproma
        IF(ldflag(jl) .AND. ABS(zcond(jl)) > 0._dp) THEN
           it = NINT(pt(jl,kk)*1000.)
! ka_sv_20170410+
! limiting to below thermosphere
!!$        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
           IF ((it<jptlucu1 .OR. it>jptlucu2) .AND. &
                (pp(jl) >= 1.0_dp)) lookupoverflow = .TRUE.
! ka_sv_20170410-
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

     IF (lookupoverflow) RETURN
330  CONTINUE

  END IF

  IF(kcall.EQ.0) THEN

     isum=0
!DIR$ IVDEP
!OCL NOVREC
     DO 410 jl=1,kproma
        it = NINT(pt(jl,kk)*1000.)
! ka_sv_20170410+
! limiting to below thermosphere
!!$        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
           IF ((it<jptlucu1 .OR. it>jptlucu2) .AND. &
                (pp(jl) >= 1.0_dp)) lookupoverflow = .TRUE.
! ka_sv_20170410-
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

     IF (lookupoverflow) RETURN

!     IF(isum.EQ.0) go to 430 ! mz_ht_20071217

!DIR$ IVDEP
!OCL NOVREC
     DO 420 jl=1,kproma
        it = NINT(pt(jl,kk)*1000.)
! ka_sv_20170410+
! limiting to below thermosphere
!!$        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
           IF ((it<jptlucu1 .OR. it>jptlucu2) .AND. &
                (pp(jl) >= 1.0_dp)) lookupoverflow = .TRUE.
! ka_sv_20170410-
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

     IF (lookupoverflow) RETURN

430  CONTINUE

  END IF

  IF(kcall.EQ.4) THEN

!DIR$ IVDEP
!OCL NOVREC
     DO 510 jl=1,kproma
        it = NINT(pt(jl,kk)*1000.)
! ka_sv_20170410+
! limiting to below thermosphere
!!$        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
           IF ((it<jptlucu1 .OR. it>jptlucu2) .AND. &
                (pp(jl) >= 1.0_dp)) lookupoverflow = .TRUE.
! ka_sv_20170410-
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

     IF (lookupoverflow) RETURN

!DIR$ IVDEP
!OCL NOVREC
     DO 520 jl=1,kproma
        it = NINT(pt(jl,kk)*1000.)
! ka_sv_20170410+
! limiting to below thermosphere
!!$        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
           IF ((it<jptlucu1 .OR. it>jptlucu2) .AND. &
                (pp(jl) >= 1.0_dp)) lookupoverflow = .TRUE.
! ka_sv_20170410-
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

     IF (lookupoverflow) RETURN

  END IF

  RETURN
END SUBROUTINE tiedtke_cuadjtq

!===========================================================================

SUBROUTINE tiedtke_cuasc(kproma, kbdim, klev, klevp1, klevm1,           &
           jrow, ztmst,                                                 &  !mz_ht_20040317+
           ptenh,    pqenh,    puen,     pven,                          &
           ktrac,                                                       &
           pxtenh,   pxten,    pxtu,     pmfuxt,                        &
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
           kcbot,    kctop,    kctop0,   kcum,                          &
           zpmfun)   ! op_ck_20031001/mz_ht_20040318 for posdef update
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
!          *CUBASMC* CALCULATE CLOUD BASE VALUES FOR MIDLEVEL CONVECTIONC
!          REFERENCE
!          ---------
!          (TIEDTKE,1989)
!


IMPLICIT NONE

!added variables  mz_ht_20040317+
REAL(dp), INTENT(IN) :: ztmst
INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, klevm1, ktrac, jrow
! mz_ht_20040317-

REAL(dp) :: ptenh(kbdim,klev),       pqenh(kbdim,klev),                &
            puen(kbdim,klev),        pven(kbdim,klev),                 &
            pten(kbdim,klev),        pqen(kbdim,klev),                 &
            pgeo(kbdim,klev),        pgeoh(kbdim,klev),                &
            paphp1(kbdim,klevp1),                                      &
            pqsen(kbdim,klev),       pqte(kbdim,klev),                 &
            pverv(kbdim,klev)
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
REAL(dp) :: pxtenh(kbdim,klev,ktrac),pxten(kbdim,klev,ktrac),          &
            pxtu(kbdim,klev,ktrac),  pmfuxt(kbdim,klev,ktrac)

! for brinkop_update
Real(dp) :: zpmfun(kbdim,klev) ! op_ck_20031001/ mz_ht_20040318


!mz_ht_20040317+ added variables

REAL(dp) :: zcons2,  zdz,    zdrodz,  zmfmax, zfac,                  &
            zmftest, zdprho,   zalvs,  zmse,    znevn,  zodmax,      &
            zqeen,   zseen,    zscde,  zga,     zdt,    zscod,       &
            zqude,   zqcod,    zmfusk, zmfuqk,  zmfulk, zxteen,      &
            zxtude,  zmfuxtk,  zbuo,   zdnoprc, zprcon, zlnew,       &
            zz,      zdmfeu,   zdmfdu, zbuoyz,  zzdmf,               &
            zzp        ! for brinkop update

INTEGER  :: kcum, is, ik, icall, ikb, ikt,   &
            jl,   jk, jt
! mz_ht_20040317-
!
!
!
!----------------------------------------------------------------------
!
!*    1.           SPECIFY PARAMETERS
!                  ------------------
!
!! 100 CONTINUE
  zcons2=1._dp/(g*ztmst)                ! mz_ht_20040317+
  zqold(:)= 0._dp   ! mz_ht_20031512
!
!
!----------------------------------------------------------------------
!
!     2.           SET DEFAULT VALUES
!                  ------------------
!
!! 200 CONTINUE
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
        zpmfun(jl,jk)=0._dp ! op_ck_20031001 for brinkop update mz_ht_20040318
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
!! 300 CONTINUE
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

  if (lpos_def) then    ! mz_ht_20040318+

! update for positive definite tracers by S.Brinkop
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
3102    END DO
3101 END DO
     !
     DO 311 jl=1,kproma
        pmfu(jl,klev)=pmfub(jl)
        pmfus(jl,klev)=pmfub(jl)*(pcpcu(jl,klev)*ptu(jl,klev)          &
                       + pgeoh(jl,klev))
        pmfuq(jl,klev)=pmfub(jl)*pqu(jl,klev)
        IF(lmfdudv) THEN
           zmfuu(jl)=pmfub(jl)*puu(jl,klev)
           zmfuv(jl)=pmfub(jl)*pvu(jl,klev)
        END IF
311  END DO
! op_ck_20031001-
  END if! mz_ht_20040318+ end of update part

!
  DO 3112 jt=1,ktrac
     DO 3110 jl=1,kproma
        IF(.NOT.ldcum(jl)) THEN
           pxtu(jl,klev,jt)=0._dp
        ENDIF
        pmfuxt(jl,klev,jt)=pmfub(jl)*pxtu(jl,klev,jt)
! update for pos_def
        if (lpos_def) &              ! mz_ht_20040318+ update for positive definite tracers by S.Brinkop
             pmfuxt(jl,klev,jt)=zpmfun(jl,klev)*pxtu(jl,klev,jt) ! op_ck_20031001
! update for pos_def
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
!! 350 CONTINUE
  DO jl=1,kproma
     IF(ktype(jl).EQ.1) THEN
        ikb=kcbot(jl)
        zbuoy(jl)=g*(ptu(jl,ikb)-ptenh(jl,ikb))/ptenh(jl,ikb) +         &
             &            g*vtmpc1*(pqu(jl,ikb)-pqenh(jl,ikb))
        IF(zbuoy(jl).GT.0._dp) THEN
           zdz=(pgeo(jl,ikb-1)-pgeo(jl,ikb))/g
           zdrodz=-LOG(pten(jl,ikb-1)/pten(jl,ikb))/zdz                 &
                         -g/(rd*ptenh(jl,ikb)*(1._dp+vtmpc1*pqenh(jl,ikb)))
! nb zoentr is here a fractional value
           zoentr(jl,ikb-1)=zbuoy(jl)*0.5_dp/(1._dp+zbuoy(jl)*zdz) + zdrodz
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
!! 400 CONTINUE
  DO 480 jk=klevm1,2,-1
!
!                  SPECIFY CLOUD BASE VALUES FOR MIDLEVEL CONVECTION
!                  IN *CUBASMC* IN CASE THERE IS NOT ALREADY CONVECTION
!                  ----------------------------------------------------
!
     ik=jk
     IF(lmfmid.AND.ik.LT.klevm1.AND.ik.GT.nmctop) THEN
        CALL tiedtke_cubasmc(kproma, kbdim, klev, ik, klab,                   &
                     pten,     pqen,     pqsen,    puen,     pven,            &
                     ktrac,                                                   &
                     pxten,    pxtu,     pmfuxt,                              &
                     pverv,    pgeo,     pgeoh,    ldcum,    ktype,           &
                     pmfu,     pmfub,    pentr,    kcbot,                     &
                     ptu,      pqu,      plu,      puu,      pvu,             &
                     pmfus,    pmfuq,    pmful,    pdmfup,   zmfuu,           &
                     pcpen,                                                   &
                     zmfuv,    paphp1,  zpmfun)  ! op_ck_20031001/mz_ht_20040318 for posdef tracers
     ENDIF
!
     is=0
     DO 410 jl=1,kproma
        is=is+klab(jl,jk+1)
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
!     IF(is.EQ.0) go to 480 ! mz_ht_20071217
!
!
!*                 SPECIFY TURBULENT ENTRAINMENT AND DETRAINMENTS
!                  RATES PLUS ORGANIZED DETRAINMENT RATES IN *CUENTR*
!                   -------------------------------------
!
     ik=jk
     CALL tiedtke_cuentr(    kproma, kbdim, klev, klevp1, ik,         &
           ptenh,    pqenh,    pqte,     paphp1,                      &
           klwmin,   ldcum,    ktype,    kcbot,    kctop0,            &
           zpbase,   pmfu,     pentr,    zodetr,                      &
           khmin,    pgeoh,                                           &
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
              zoentr(jl,jk)=MAX(zoentr(jl,jk)-MAX(zmftest-zmfmax,0._dp),0._dp)
           ELSE
              zoentr(jl,jk)=0._dp
           ENDIF
           IF(ktype(jl).EQ.1.AND.jk.LT.kcbot(jl).AND.jk.LE.khmin(jl)) THEN
!          limit organized detrainment to not allowing for too
!          deep clouds
              zalvs=MERGE(alv,als,ptu(jl,jk+1)>tmelt)
              zmse=pcpcu(jl,jk+1)*ptu(jl,jk+1)+zalvs*pqu(jl,jk+1)      &
                                                 +pgeoh(jl,jk+1)
              ikt=kctop0(jl)
              znevn=(pgeoh(jl,ikt)-pgeoh(jl,jk+1))*(zmse-phhatt(jl,jk+1))/g
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
                &         (1._dp/ptenh(jl,jk+1) + vtmpc1*zga)
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
           ptu(jl,jk)=(zmfusk*(1._dp/MAX(cmfcmin,pmfu(jl,jk)))-           &
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

              if (.not.lpos_def) then     ! mz_ht_20040318+ changed 
!                                         for positiv definite tracers update 
                 zxteen=pxtenh(jl,jk+1,jt)*(zdmfen(jl)+zoentr(jl,jk))
                 zxtude=pxtu(jl,jk+1,jt)*(zdmfde(jl)+zodetr(jl,jk))
                 zmfuxtk=pmfuxt(jl,jk+1,jt)+zxteen-zxtude
                 pxtu(jl,jk,jt)=zmfuxtk*(1._dp/MAX(cmfcmin,pmfu(jl,jk)))
              ELSE ! mz_ht_20040318+

! update for pos_def
! op_ck_20031001+   
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
                    pmfuxt(jl,jk,jt)=zmfuxtk/(1._dp+(zdmfde(jl)+zodetr(jl,jk)) &
                         /max(cmfcmin,pmfu(jl,jk)))
                    pxtu(jl,jk,jt)=pmfuxt(jl,jk,jt)/max(cmfcmin,pmfu(jl,jk))
                 ENDIF
! op_ck_20031001-
              ENDIF ! lpos_def     !mz_ht_20040318
! end of update for pos_def

           ENDIF
4202    END DO
4204 END DO
!
!! mz_ht_20040318+   for channel output
!     do jl=1,kproma
! channel objects set
!        conv_top(jl,jrow)      = REAL(kctop(jl),dp) 
!        massfu_asc(jl,jk,jrow) = pmfu(jl,jk)
!        u_entr(jl,jk,jrow)     = zoentr(jl,jk)+zdmfen(jl)
!        u_detr(jl,jk,jrow)     = zodetr(jl,jk)+zdmfde(jl)
!     enddo

! op_mm_20140131+
! use 2d fields instead of 3d fields 
     do jl=1,kproma
! channel objects set
!!$        conv_top(jl,jrow)      = REAL(kctop(jl),dp) ! op_mm_20140520
        conv_top_1d(jl)      = REAL(kctop(jl),dp)
        massfu_asc_2d(jl,jk) = pmfu(jl,jk)
        u_entr_2d(jl,jk)     = zoentr(jl,jk)+zdmfen(jl)
        u_detr_2d(jl,jk)     = zodetr(jl,jk)+zdmfde(jl)
     enddo
     ! op_mm_20140131-
! mz_ht_20031712-
!
!                  DO CORRECTIONS FOR MOIST ASCENT
!                  BY ADJUSTING T,Q AND L IN *CUADJTQ*
!                  -----------------------------------
!
     ik=jk
     icall=1
     CALL tiedtke_cuadjtq(kproma, kbdim, klev, ik,    &
          &       zph,      ptu,      pqu,      loflag,   icall)
     if (lookupoverflow) RETURN

     ! mz_ak_20051221+
     ! store change in liquid water in channel object
 !    del_liqwat(1:kproma,jk,jrow) = zqold(1:kproma)-pqu(1:kproma,jk) 
     ! op_ mm_20140131
     ! use 2d fields 
      del_liqwat_2d(1:kproma,jk) = zqold(1:kproma)-pqu(1:kproma,jk) 
     ! mz_ak_20051221-

!
!DIR$ IVDEP
!OCL NOVREC
     DO 440 jl=1,kproma
        IF(loflag(jl).AND.pqu(jl,jk).LT.zqold(jl)) THEN
           klab(jl,jk)=2
           plu(jl,jk)=plu(jl,jk)+zqold(jl)-pqu(jl,jk)
           
           ! mz_ht_20070325+
           IF (pten(jl,jk) > tmelt) THEN
           !!  cv_lwc(jl,jk,jrow) = plu(jl,jk)  ! op_mm_20140327 use 2d_fields 
              cv_lwc_2d(jl,jk) = plu(jl,jk) 
           ELSE
 !!            cv_iwc(jl,jk,jrow) = plu(jl,jk) 
              cv_iwc_2d(jl,jk) = plu(jl,jk) 
           ENDIF
           ! mz_ht_20070325-
           zbuo=ptu(jl,jk)*(1._dp+vtmpc1*pqu(jl,jk)-plu(jl,jk))-             &
                &   ptenh(jl,jk)*(1._dp+vtmpc1*pqenh(jl,jk))
           IF(klab(jl,jk+1).EQ.1) zbuo=zbuo+0.5_dp
           IF(zbuo.GT.0._dp.AND.pmfu(jl,jk).GE.0.01_dp*pmfub(jl).AND.           &
                &      jk.GE.kctop0(jl)) THEN
              kctop(jl)=jk
              ldcum(jl)=.TRUE.
              zdnoprc=MERGE(3.e4_dp,1.5e4_dp,ldland(jl))
              zprcon=MERGE(0._dp,cprcon,zpbase(jl)-paphp1(jl,jk).LT.zdnoprc)
              zlnew=plu(jl,jk)/(1._dp+zprcon*(pgeoh(jl,jk)-pgeoh(jl,jk+1)))
              pdmfup(jl,jk)=MAX(0._dp,(plu(jl,jk)-zlnew)*pmfu(jl,jk))
              ! mz_ht_20070325+
              IF (pten(jl,jk) > tmelt) THEN
 !!               cv_rform(jl,jk,jrow) = MAX(0._dp,(plu(jl,jk)-zlnew)) ! op_mm_20140327
                cv_rform_2d(jl,jk) = MAX(0._dp,(plu(jl,jk)-zlnew))     ! use 2d fields
              ELSE
                 cv_sform_2d(jl,jk) = MAX(0._dp,(plu(jl,jk)-zlnew))
!!                cv_sform(jl,jk,jrow) = MAX(0._dp,(plu(jl,jk)-zlnew))
              ENDIF
              ! mz_ht_20070325-
              plu(jl,jk)=zlnew
           ELSE 
              if (lpos_def) zpmfun(jl,jk)=0._dp ! op_ck_20031001/mz_ht_20040318
              klab(jl,jk)=0
              pmfu(jl,jk)=0._dp
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
              pmfuxt(jl,jk,jt)=pxtu(jl,jk,jt)*pmfu(jl,jk)
              if (lpos_def) pmfuxt(jl,jk,jt)=pxtu(jl,jk,jt)*zpmfun(jl,jk) ! op_ck_20031001/mz_ht_20040318+
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
!!                 zz=MERGE(3._dp,2._dp,zdmfen(jl).EQ.0._dp)
                 zz=MERGE(3._dp,2._dp,ABS(zdmfen(jl)).lt.TINY(0._dp))
              ELSE
!!                 zz=MERGE(1._dp,0._dp,zdmfen(jl).EQ.0._dp)
                 zz=MERGE(1._dp,0._dp,ABS(zdmfen(jl)).lt.TINY(0._dp))
              END IF
              zdmfeu=zdmfen(jl)+zz*zdmfde(jl)
              zdmfdu=zdmfde(jl)+zz*zdmfde(jl)
              zdmfdu=MIN(zdmfdu,0.75_dp*pmfu(jl,jk+1))
              zmfuu(jl)=zmfuu(jl)+                                      &
                   &         zdmfeu*puen(jl,jk)-zdmfdu*puu(jl,jk+1)
              zmfuv(jl)=zmfuv(jl)+                                      &
                   &         zdmfeu*pven(jl,jk)-zdmfdu*pvu(jl,jk+1)
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
           zbuoyz=g*(ptu(jl,jk)-ptenh(jl,jk))/ptenh(jl,jk) +               &
                &         g*vtmpc1*(pqu(jl,jk)-pqenh(jl,jk))-g*plu(jl,jk)
           zbuoyz=MAX(zbuoyz,0._dp)
           zdz=(pgeo(jl,jk-1)-pgeo(jl,jk))/g
           zdrodz=-LOG(pten(jl,jk-1)/pten(jl,jk))/zdz                      &
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
!! 500 CONTINUE
  DO 510 jl=1,kproma
     IF(kctop(jl).EQ.klevm1) ldcum(jl)=.FALSE.
     kcbot(jl)=MAX(kcbot(jl),kctop(jl))
510 END DO
  is=0
  DO 520 jl=1,kproma
     is=is+MERGE(1,0,ldcum(jl))
520 END DO
  kcum=is
!  IF(is.EQ.0) go to 800 ! mz_ht_20071217
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

! update for pos_def
        if (lpos_def) zpmfun(jl,jk)=pmfu(jl,jk) ! op_ck_20031001/mz_ht_20040318

     END IF
530 END DO
  DO 5312 jt=1,ktrac
     DO 5310 jl=1,kproma
        IF(ldcum(jl)) THEN
           jk=kctop(jl)-1
           pmfuxt(jl,jk,jt)=pxtu(jl,jk,jt)*pmfu(jl,jk)

! update for pos_def
           if (lpos_def) pmfuxt(jl,jk,jt)=zpmfun(jl,jk)*pxtu(jl,jk,jt) ! op_ck_20031001
  
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
800 CONTINUE
!
  RETURN
END SUBROUTINE tiedtke_cuasc

!=============================================================================

SUBROUTINE tiedtke_cuasct(   kproma, kbdim, klev, klevp1, klevm1,      &
           jrow, ztmst,                                                & !mz_ht_20040317+ added ztmst
           ptenh,    pqenh,    puen,     pven,                         &
           ktrac,                                                      &
           pxtenh,   pxten,    pxtu,     pmfuxt,                       &
           pten,     pqen,     pqsen,                                  &
           pgeo,     pgeoh,    paphp1,                                 &
           pqte,     pverv,    klwmin,                                 &
           ldcum,    ldland,   ktype,    klab,                         &
           ptu,      pqu,      plu,      puu,      pvu,                &
           pmfu,     pmfub,    pentr,                                  &
           pmfus,    pmfuq,                                            &
           pmful,    plude,    pqude,    pdmfup,                       &
           pcpen,    pcpcu,                                            &
           kcbot,    kctop,    kctop0,   kcum,                         &
           zpmfun)   ! op_ck_20031001/mz_ht_20040318 for posdef update
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
!          *CUBASMC* CALCULATE CLOUD BASE VALUES FOR MIDLEVEL CONVECTIONC
!
!          REFERENCE
!          ---------
!          (TIEDTKE,1989)
!


!added variables  mz_ht_20040317+
REAL(dp), INTENT(IN) :: ztmst
INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, klevm1, ktrac, jrow
! mz_ht_20040317-
!
REAL(dp) :: ptenh(kbdim,klev),       pqenh(kbdim,klev),                   &
            puen(kbdim,klev),        pven(kbdim,klev),                    &
            pten(kbdim,klev),        pqen(kbdim,klev),                    &
            pgeo(kbdim,klev),        pgeoh(kbdim,klev),                   &
            paphp1(kbdim,klevp1),                                         &
            pqsen(kbdim,klev),       pqte(kbdim,klev),                    &
            pverv(kbdim,klev)
!
REAL(dp) :: ptu(kbdim,klev),         pqu(kbdim,klev),                     &
            puu(kbdim,klev),         pvu(kbdim,klev),                     &
            pmfu(kbdim,klev),                                             &
            pmfub(kbdim),            pentr(kbdim),                        &
            pmfus(kbdim,klev),       pmfuq(kbdim,klev),                   &
            plu(kbdim,klev),         plude(kbdim,klev),                   &
            pqude(kbdim,klev),                                            &
            pmful(kbdim,klev),       pdmfup(kbdim,klev)
REAL(dp) :: pcpen(kbdim,klev),       pcpcu(kbdim,klev)
INTEGER  :: klwmin(kbdim),           ktype(kbdim),                        &
            klab(kbdim,klev),        kcbot(kbdim),                        &
            kctop(kbdim),            kctop0(kbdim)
LOGICAL  :: ldcum(kbdim),            ldland(kbdim)
!
REAL(dp) :: zdmfen(kbdim),           zdmfde(kbdim),                       &
            zmfuu(kbdim),            zmfuv(kbdim),                        &
            zpbase(kbdim),           zqold(kbdim)
REAL(dp) :: zph(kbdim)
LOGICAL  :: loflag(kbdim)
REAL(dp) :: pxtenh(kbdim,klev,ktrac),pxten(kbdim,klev,ktrac),             &
            pxtu(kbdim,klev,ktrac),  pmfuxt(kbdim,klev,ktrac)

! for brinkop_update
Real(dp) :: zpmfun(kbdim,klev) ! op_ck_20031001/ mz_ht_20040318
! mz_ht_20040317+

REAL(dp) :: zcons2,  zmfmax, zfac,                         &
            zmftest, zqeen,    zseen,  zscde,                        &
            zqude,   zmfusk,   zmfuqk, zmfulk,  zxteen,              &
            zxtude,  zmfuxtk,  zbuo,   zdnoprc, zprcon, zlnew,       &
            zz,      zdmfeu,   zdmfdu, zzdmf,                        &
            zzp          ! for brinkop update

INTEGER  :: kcum, is, ik, ikb, icall,    &
            jl,   jk, jt
! mz_ht_20040317-
!
!
!
!----------------------------------------------------------------------
!
!*    1.           SPECIFY PARAMETERS
!                  ------------------
!
!! 100 CONTINUE
  zcons2=1._dp/(g*ztmst)        !mz_ht_20040317 changed time_step_len to ztmst
  zqold(:)= 0._dp   ! mz_ht_20031512
!
!----------------------------------------------------------------------
!
!     2.           SET DEFAULT VALUES
!                  ------------------
!
!! 200 CONTINUE
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
        zpmfun(jl,jk)=0._dp ! op_ck_20031001 for brinkop update mz_ht_20040318
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
 
!
!----------------------------------------------------------------------
!
!     3.0          INITIALIZE VALUES AT LIFTING LEVEL
!                  ----------------------------------
!
!! 300 CONTINUE
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

! mz_ht_20040318+ update for posdef
  if (lpos_def) then
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
3102    END DO
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
311  END DO
! op_ck_20031001-
  ENDIF   !lpos_def 
!mz_ht_20040318- end of update for posdef

  DO 3112 jt=1,ktrac
     DO 3110 jl=1,kproma
        IF(.NOT.ldcum(jl)) THEN
           pxtu(jl,klev,jt)=0._dp
        ENDIF
        pmfuxt(jl,klev,jt)=pmfub(jl)*pxtu(jl,klev,jt)
    
        ! mz_ht_20040318+ update for posdef
        if (lpos_def) pmfuxt(jl,klev,jt)=zpmfun(jl,klev)*pxtu(jl,klev,jt) ! op_ck_20031001
        !mz_ht_20040318- end of update for posdef

3110 END DO
3112 END DO
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
!! 400 CONTINUE
  DO 480 jk=klevm1,2,-1
!
!                  SPECIFY CLOUD BASE VALUES FOR MIDLEVEL CONVECTION
!                  IN *CUBASMC* IN CASE THERE IS NOT ALREADY CONVECTION
!                  ----------------------------------------------------
!
     ik=jk
     IF(lmfmid.AND.ik.LT.klevm1.AND.ik.GT.nmctop) THEN
        CALL tiedtke_cubasmc(kproma,   kbdim,    klev,    ik,    klab, &
                     pten,     pqen,     pqsen,    puen,     pven,     &
                     ktrac,                                            &
                     pxten,    pxtu,     pmfuxt,                       &
                     pverv,    pgeo,     pgeoh,    ldcum,    ktype,    &
                     pmfu,     pmfub,    pentr,    kcbot,              &
                     ptu,      pqu,      plu,      puu,      pvu,      &
                     pmfus,    pmfuq,    pmful,    pdmfup,   zmfuu,    &
                     pcpen,                                            &
                     zmfuv,    paphp1,  zpmfun)  ! op_ck_20031001/mz_ht_20040318 for posdef tracers)
     ENDIF
!
     is=0
     DO 410 jl=1,kproma
        is=is+klab(jl,jk+1)
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
!     IF(is.EQ.0) go to 480  ! mz_ht_20071217
!
!
!*                 SPECIFY ENTRAINMENT RATES IN *CUENTRT*
!                  --------------------------------------
!
     ik=jk
     CALL tiedtke_cuentrt(kproma,  kbdim,    klev,     klevp1,   ik,   &
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
           ptu(jl,jk)=(zmfusk*(1._dp/MAX(cmfcmin,pmfu(jl,jk)))-           &
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
              
!mz_ht_20040318+ update for posdef
              if (lpos_def) then

! op_ck_20031001+       
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
                    zxteen=pxten(jl,jk,jt)*zdmfen(jl)
                    zmfuxtk=pmfuxt(jl,jk+1,jt)+zxteen
                    pmfuxt(jl,jk,jt)=zmfuxtk/(1._dp+(zdmfde(jl) &
                         /max(cmfcmin,pmfu(jl,jk))))
                    pxtu(jl,jk,jt)=pmfuxt(jl,jk,jt)/max(cmfcmin,pmfu(jl,jk))
                 ENDIF
! op_ck_20031001-

              ELSE
!mz_ht_20040318- end of update for posdef 

                 zxteen=pxtenh(jl,jk+1,jt)*zdmfen(jl)
                 zxtude=pxtu(jl,jk+1,jt)*zdmfde(jl)
                 zmfuxtk=pmfuxt(jl,jk+1,jt)+zxteen-zxtude
                 pxtu(jl,jk,jt)=zmfuxtk*(1._dp/MAX(cmfcmin,pmfu(jl,jk)))
              ENDIF   !lpos_def mz_ht_20040318+
           ENDIF
4202    END DO
4204 END DO
!
!! mz_ht_20040318+   for channel output
     do jl=1,kproma
! channel objects set
!!$        conv_top(jl,jrow)      = REAL(kctop(jl),dp)
!!$        massfu_asc(jl,jk,jrow) = pmfu(jl,jk)
!!$        u_entr(jl,jk,jrow)     = zdmfen(jl)
!!$        u_detr(jl,jk,jrow)     = zdmfde(jl)
! op_mm_20140327   for 1d/2d fields 
       conv_top_1d(jl)      = REAL(kctop(jl),dp)
       massfu_asc_2d(jl,jk) = pmfu(jl,jk)
       u_entr_2d(jl,jk)     = zdmfen(jl)
       u_detr_2d(jl,jk)     = zdmfde(jl)
     enddo
! mz_ht_20031712-

!
!                  DO CORRECTIONS FOR MOIST ASCENT
!                  BY ADJUSTING T,Q AND L IN *CUADJTQ*
!                  -----------------------------------
!
     ik=jk
     icall=1
     CALL tiedtke_cuadjtq(kproma,   kbdim,    klev,     ik,            &
                  zph,      ptu,      pqu,      loflag,   icall)
     IF (lookupoverflow) RETURN
!
!DIR$ IVDEP
!OCL NOVREC
     DO 440 jl=1,kproma
        IF(loflag(jl).AND.pqu(jl,jk).LT.zqold(jl)) THEN
           klab(jl,jk)=2
           plu(jl,jk)=plu(jl,jk)+zqold(jl)-pqu(jl,jk)
           ! mz_ht_20070325+
           IF (pten(jl,jk) > tmelt) THEN
       !!      cv_lwc(jl,jk,jrow) = plu(jl,jk)  ! op_mm_20140327 use 2d fields 
               cv_lwc_2d(jl,jk) = plu(jl,jk)  
           ELSE
              !!        cv_iwc(jl,jk,jrow) = plu(jl,jk) 
              cv_iwc_2d(jl,jk) = plu(jl,jk) 
           ENDIF
           ! mz_ht_20070325-
           zbuo=ptu(jl,jk)*(1._dp+vtmpc1*pqu(jl,jk)-plu(jl,jk))-          &
                    ptenh(jl,jk)*(1._dp+vtmpc1*pqenh(jl,jk))
           IF(klab(jl,jk+1).EQ.1) zbuo=zbuo+0.5_dp
           IF(zbuo.GT.0._dp.AND.pmfu(jl,jk).GE.0.1_dp*pmfub(jl)) THEN
              kctop(jl)=jk
              ldcum(jl)=.TRUE.
              zdnoprc=MERGE(3.e4_dp,1.5e4_dp,ldland(jl))
              zprcon=MERGE(0._dp,cprcon,                                  &
                               zpbase(jl)-paphp1(jl,jk).LT.zdnoprc)
              zlnew=plu(jl,jk)/(1._dp+zprcon*(pgeoh(jl,jk)-pgeoh(jl,jk+1)))
              pdmfup(jl,jk)=MAX(0._dp,(plu(jl,jk)-zlnew)*pmfu(jl,jk))
              ! mz_ht_20070325+
              IF (pten(jl,jk) > tmelt) THEN
!!                cv_rform(jl,jk,jrow) = MAX(0._dp,(plu(jl,jk)-zlnew)) ! op_mm_20140327
                 cv_rform_2d(jl,jk) = MAX(0._dp,(plu(jl,jk)-zlnew))    ! use 2d fields 
              ELSE
 !!               cv_sform(jl,jk,jrow) = MAX(0._dp,(plu(jl,jk)-zlnew))
                 cv_sform_2d(jl,jk) = MAX(0._dp,(plu(jl,jk)-zlnew))
              ENDIF
              ! mz_ht_20070325-
              plu(jl,jk)=zlnew
           ELSE
              if (lpos_def) zpmfun(jl,jk)=0._dp ! op_ck_20031001/mz_ht_20040318 for posdef update
              klab(jl,jk)=0
              pmfu(jl,jk)=0._dp
           END IF
        END IF
440  END DO
     DO 455 jl=1,kproma
        IF(loflag(jl)) THEN
           pmful(jl,jk)=plu(jl,jk)*pmfu(jl,jk)
           pmfus(jl,jk)=(pcpcu(jl,jk)*ptu(jl,jk)+pgeoh(jl,jk))         &
                                 *pmfu(jl,jk)
           pmfuq(jl,jk)=pqu(jl,jk)*pmfu(jl,jk)
        END IF
455  END DO
     DO 4554 jt=1,ktrac
        DO 4552 jl=1,kproma
           IF(loflag(jl)) THEN
              pmfuxt(jl,jk,jt)=pxtu(jl,jk,jt)*pmfu(jl,jk)

              ! mz_ht_20040318+ update for posdef
              if (lpos_def) pmfuxt(jl,klev,jt)=zpmfun(jl,klev)*pxtu(jl,klev,jt) ! op_ck_20031001
              !mz_ht_20040318- end of update for posdef

           ENDIF
4552    END DO
4554 END DO
!
     IF(lmfdudv) THEN
        DO 460 jl=1,kproma
           IF(loflag(jl)) THEN
              IF(ktype(jl).EQ.1.OR.ktype(jl).EQ.3) THEN
!!                zz=MERGE(3._dp,2._dp,zdmfen(jl).EQ.0._dp)
                zz=MERGE(3._dp,2._dp,ABS(zdmfen(jl)).lt.TINY(0._dp))
              ELSE
!!                zz=MERGE(1._dp,0._dp,zdmfen(jl).EQ.0._dp)
                zz=MERGE(1._dp,0._dp,ABS(zdmfen(jl)).lt.TINY(0._dp))
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
!! 500 CONTINUE
  DO 510 jl=1,kproma
     IF(kctop(jl).EQ.klevm1) ldcum(jl)=.FALSE.
     kcbot(jl)=MAX(kcbot(jl),kctop(jl))
510 END DO
  is=0
  DO 520 jl=1,kproma
     is=is+MERGE(1,0,ldcum(jl))
520 END DO
  kcum=is
!  IF(is.EQ.0) go to 800 ! mz_ht_20071217
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
! update for pos_def
        if (lpos_def) zpmfun(jl,jk)=pmfu(jl,jk) ! op_ck_20031001/mz_ht_20040318

530  END DO
  DO 5312 jt=1,ktrac
     DO 5310 jl=1,kproma
        IF(ldcum(jl)) THEN
           jk=kctop(jl)-1
           pmfuxt(jl,jk,jt)=pxtu(jl,jk,jt)*pmfu(jl,jk)
        
! update for pos_def
           if (lpos_def) pmfuxt(jl,jk,jt)=zpmfun(jl,jk)*pxtu(jl,jk,jt) ! op_ck_20031001
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
800 CONTINUE
!
  RETURN
END SUBROUTINE tiedtke_cuasct

!==============================================================================

SUBROUTINE tiedtke_cubase(   kproma, kbdim, klev, klevp1, klevm1,      &
           ptenh,    pqenh,    pgeoh,    paph,                         &
           ptu,      pqu,      plu,                                    &
           puen,     pven,     puu,      pvu,                          &
           pcpcu,                                                      &
           ldcum,    kcbot,    klab)
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

REAL(dp) ::    ptenh(kbdim,klev),       pqenh(kbdim,klev),                   &
           pgeoh(kbdim,klev),       paph(kbdim,klevp1)
!
REAL(dp) ::    ptu(kbdim,klev),         pqu(kbdim,klev),                     &
           plu(kbdim,klev)
REAL(dp) ::    puen(kbdim,klev),        pven(kbdim,klev),                    &
           puu(kbdim,klev),         pvu(kbdim,klev)
REAL(dp) ::    pcpcu(kbdim,klev)   
INTEGER :: klab(kbdim,klev),        kcbot(kbdim)
LOGICAL :: ldcum(kbdim)
!
REAL(dp) ::    zqold(kbdim)
REAL(dp) ::    zph(kbdim)
LOGICAL :: loflag(kbdim)

INTEGER :: jl, jk, is, ik, ikb, icall
REAL(dp)    :: zbuo, zz
!
!
!----------------------------------------------------------------------
!
!     1.           INITIALIZE VALUES AT LIFTING LEVEL
!                  ----------------------------------
!
!! 100 CONTINUE
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
!! 200 CONTINUE
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
           zbuo=ptu(jl,jk)*(1._dp+vtmpc1*pqu(jl,jk))-                     &
                          ptenh(jl,jk)*(1._dp+vtmpc1*pqenh(jl,jk))+0.5_dp
           IF(zbuo.GT.0._dp) klab(jl,jk)=1
           zqold(jl)=pqu(jl,jk)
        END IF
220  END DO
!
     ik=jk
     icall=1
     CALL tiedtke_cuadjtq(kproma, kbdim, klev, ik,                             &
                  zph,      ptu,      pqu,      loflag,   icall)
!
     IF (lookupoverflow) RETURN
!DIR$ IVDEP
!OCL NOVREC
     DO 240 jl=1,kproma
        IF(loflag(jl).AND.pqu(jl,jk).LT.zqold(jl)) THEN
           klab(jl,jk)=2
           plu(jl,jk)=plu(jl,jk)+zqold(jl)-pqu(jl,jk)
           zbuo=ptu(jl,jk)*(1._dp+vtmpc1*pqu(jl,jk)-plu(jl,jk))-          &
                         ptenh(jl,jk)*(1._dp+vtmpc1*pqenh(jl,jk))+0.5_dp
           IF(zbuo.GT.0._dp) THEN
              kcbot(jl)=jk
              ldcum(jl)=.TRUE.
           END IF
        END IF
240  END DO
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
END SUBROUTINE tiedtke_cubase

!==============================================================================
SUBROUTINE tiedtke_cubasmc(  kproma, kbdim, klev, kk, klab,            &
           pten,     pqen,     pqsen,    puen,     pven,               &
           ktrac,                                                      &
           pxten,    pxtu,     pmfuxt,                                 &
           pverv,    pgeo,     pgeoh,    ldcum,    ktype,              &
           pmfu,     pmfub,    pentr,    kcbot,                        &
           ptu,      pqu,      plu,      puu,      pvu,                &
           pmfus,    pmfuq,    pmful,    pdmfup,   pmfuu,              &
           pcpen,                                                      &
           pmfuv,    paphp1,  zpmfun) ! op_ck_20031001/mz_ht_20040318 for brinkop update
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


! mz_ht_20040317+ 
INTEGER, INTENT(IN) :: kproma, klev, ktrac, kbdim, kk
INTEGER ::             jl, jt, jk, ikb

REAL(dp) :: zzzmb, zzp
REAL(dp) :: zpmfun(kbdim,klev)   ! op_ck_20031001
REAL(dp) :: paphp1(kbdim,klev+1) ! op_ck_20031001
! mz_ht_20040317-

!
REAL(dp) :: pten(kbdim,klev),        pqen(kbdim,klev),                    &
            puen(kbdim,klev),        pven(kbdim,klev),                    &
            pqsen(kbdim,klev),       pverv(kbdim,klev),                   &
            pgeo(kbdim,klev),        pgeoh(kbdim,klev)
!
REAL(dp) :: ptu(kbdim,klev),         pqu(kbdim,klev),                     &
            puu(kbdim,klev),         pvu(kbdim,klev),                     &
            plu(kbdim,klev),         pmfu(kbdim,klev),                    &
            pmfub(kbdim),            pentr(kbdim),                        &
            pmfus(kbdim,klev),       pmfuq(kbdim,klev),                   &
            pmful(kbdim,klev),       pdmfup(kbdim,klev),                  &
            pmfuu(kbdim),            pmfuv(kbdim)
REAL(dp) :: pcpen(kbdim,klev)
INTEGER  :: ktype(kbdim),            kcbot(kbdim),                        &
            klab(kbdim,klev)
LOGICAL  :: ldcum(kbdim)
!
REAL(dp) :: pxten(kbdim,klev,ktrac), pxtu(kbdim,klev,ktrac),              &
            pmfuxt(kbdim,klev,ktrac)
LOGICAL  :: llo3(kbdim)
!
!
!----------------------------------------------------------------------
!
!*    1.           CALCULATE ENTRAINMENT AND DETRAINMENT RATES
!                  -------------------------------------------
!
!! 100 CONTINUE
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

!       mz_ht_20040318+
! update for pos_def
        if (lpos_def) zpmfun(jl,kk+1)=pmfub(jl) ! op_ck_20031001
!       mz_ht_20040318-

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
!DIR$ IVDEP
!OCL NOVREC

!       mz_ht_20040318+
! update for pos_def

  if (lpos_def) then     !mz_ht_20040318+

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
! end uf update for pos_def

  ELSE   ! mz_ht_20040318+ 

     DO 1506 jt=1,ktrac
        DO 1502 jl=1,kproma
           IF (llo3(jl)) THEN
              pxtu(jl,kk+1,jt)=pxten(jl,kk,jt)
              pmfuxt(jl,kk+1,jt)=pmfub(jl)*pxtu(jl,kk+1,jt)
           ENDIF
1502    END DO
1506 END DO

  ENDIF            ! mz_ht_20040318+  end of lpos_def
!
!
  RETURN
END SUBROUTINE tiedtke_cubasmc

!==============================================================================
SUBROUTINE tiedtke_cuddraf(  kproma, kbdim, klev, klevp1, jrow,        &   ! mz_ht_20040318+ jrow added
           ptenh,    pqenh,    puen,     pven,                         &
           ktrac,                                                      &
           pxtenh,   pxtd,     pmfdxt,                                 &
           pgeoh,    paphp1,   prfl,                                   &
           ptd,      pqd,      pud,      pvd,                          &
           pmfd,     pmfds,    pmfdq,    pdmfdp,                       &
           pcpcu,                                                      &
           lddraf,                                                     &
           pxten, & ! op_ck_20031001/mz_ht_20040318   for pos_def update
           pvddraf) ! op_mm_20140226 for convective gust
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

!
!mz_ht_20040317+ added variables
INTEGER, INTENT(IN):: kproma, kbdim, klev, ktrac, klevp1, jrow
!mz_ht_20040317-
INTEGER  :: jl, jk, jt, is, ik, itopde, icall

REAL(dp) :: ptenh(kbdim,klev),       pqenh(kbdim,klev),                   &
            puen(kbdim,klev),        pven(kbdim,klev),                    &
            pgeoh(kbdim,klev),       paphp1(kbdim,klevp1)
!
REAL(dp) :: ptd(kbdim,klev),         pqd(kbdim,klev),                     &
            pud(kbdim,klev),         pvd(kbdim,klev),                     &
            pmfd(kbdim,klev),        pmfds(kbdim,klev),                   &
            pmfdq(kbdim,klev),       pdmfdp(kbdim,klev),                  &
            prfl(kbdim)      ,       pvddraf(kbdim)  ! op_mm_20140226 for convective gust
REAL(dp) :: pcpcu(kbdim,klev)
LOGICAL  :: lddraf(kbdim)
!
REAL(dp) :: zdmfen(kbdim),           zdmfde(kbdim),                       &
            zcond(kbdim)
REAL(dp) :: zph(kbdim) 
! op_mm_20140327+
! for tke and convective gust 
REAL(dp) :: zvbuo(kbdim)
REAL(dp) :: ztddraf, zqprec  
! op_mm_20140327-

LOGICAL  :: llo2(kbdim)
REAL(dp) :: pxtenh(kbdim,klev,ktrac),pxtd(kbdim,klev,ktrac),              &
            pmfdxt(kbdim,klev,ktrac)
REAL(dp) :: pxten(kbdim,klev,ktrac) ! op_ck_20031001/mz_ht_20040318 for posdef update
LOGICAL  :: llo1

!mz_ht_20040317+ added local variables
REAL(dp) :: zentr,  zseen,  zqeen,   zsdde, zqdde,  zmfdsk, zmfdqk,         &
            zxteen, zxtdde, zmfdxtk, zbuo,  zdmfdp, zmfduk, zmfdvk
!mz_ht_20040317-
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
    zdmfde(:) = 0._dp
    zdmfen(:) = 0._dp
! op_mm_20140327+
    zvbuo(:)=0._dp
    pvddraf(:)=0._dp
! op_mm_20140327-
!! 100 CONTINUE
  DO 180 jk=3,klev
     is=0
     DO 110 jl=1,kproma
        zph(jl)=paphp1(jl,jk)
        llo2(jl)=lddraf(jl).AND.pmfd(jl,jk-1).LT.0._dp
        is=is+MERGE(1,0,llo2(jl))
110  END DO
!     IF(is.EQ.0) go to 180 ! mz_ht_20071217
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
           ptd(jl,jk)=(zmfdsk*(1._dp/MIN(-cmfcmin,pmfd(jl,jk)))-          &
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

! mz_ht_20040318+ update for posdef
              if (lpos_def) then
  ! op_ck_20031001+
                 zxteen=pxten(jl,jk-1,jt)*zdmfen(jl)
                 zmfdxtk=pmfdxt(jl,jk-1,jt)+zxteen
                 pmfdxt(jl,jk,jt)=zmfdxtk/(1._dp+zdmfde(JL) /                 &
                                  (min(pmfd(jl,jk),-cmfcmin)))
                 pxtd(jl,jk-1,jt)=pmfdxt(jl,jk,jt)       /                 &
                                  (min(pmfd(JL,JK),-cmfcmin))
  ! op_ck_20031001-
              ELSE
                 zxteen=pxtenh(jl,jk-1,jt)*zdmfen(jl)
                 zxtdde=pxtd(jl,jk-1,jt)*zdmfde(jl)
                 zmfdxtk=pmfdxt(jl,jk-1,jt)+zxteen-zxtdde
                 pxtd(jl,jk,jt)=zmfdxtk*(1._dp/MIN(-cmfcmin,pmfd(jl,jk)))
              ENDIF
! mz_ht_20040318- end of update for posdef
           ENDIF
1262    END DO
1264 END DO
!
!
     ik=jk
     icall=2
     CALL tiedtke_cuadjtq(kproma,   kbdim,    klev,     ik,                    &
                  zph,      ptd,      pqd,      llo2,     icall)
!
     IF (lookupoverflow) RETURN
     DO 150 jl=1,kproma
        IF(llo2(jl)) THEN
           zcond(jl)=zcond(jl)-pqd(jl,jk)
           zbuo=ptd(jl,jk)*(1._dp+vtmpc1*pqd(jl,jk))-                     &
                       ptenh(jl,jk)*(1._dp+vtmpc1*pqenh(jl,jk))
           llo1=zbuo.LT.0._dp.AND.(prfl(jl)-pmfd(jl,jk)*zcond(jl).GT.0._dp)
           pmfd(jl,jk)=MERGE(pmfd(jl,jk),0._dp,llo1)
           pmfds(jl,jk)=(pcpcu(jl,jk)*ptd(jl,jk)+pgeoh(jl,jk))         &
                                                    *pmfd(jl,jk)
           pmfdq(jl,jk)=pqd(jl,jk)*pmfd(jl,jk)
           zdmfdp=-pmfd(jl,jk)*zcond(jl)
           pdmfdp(jl,jk-1)=zdmfdp
           prfl(jl)=prfl(jl)+zdmfdp
           ! op_mm_20140226+
           ! copied from scr_conv_tiedtke to calculate gusts
           ! check if half levels are correct 
           ! pt -> ptenh
           ! pt_d -> ptd 
           IF (llo1) THEN
              ztddraf    = 0.5_dp*(ptd(jl,jk-1) + ptd(jl,jk)) 
              zqprec     = 0.0_dp              ! Neglect rain water, it gives
                                               ! unrealistic high values
              zvbuo(jl)   = zvbuo(jl) + 2.0_dp                         &
                      *((ptenh(jl,jk-1)-ztddraf  )/ptenh(jl,jk-1) + zqprec   ) &
                      *rd*ptenh(jl,jk-1)/                          &
                      (g*paphp1(jl,jk-1))*                         &
                      (paphp1(jl,jk)-paphp1(jl,jk-1))
           ENDIF
           ! op_m_20140226-


        END IF
150  END DO
!
     DO 1504 jt=1,ktrac
        DO 1502 jl=1,kproma
           IF(llo2(jl)) THEN
              pmfdxt(jl,jk,jt)=pxtd(jl,jk,jt)*pmfd(jl,jk)
! mz_ht_20040318+  update for pos_def              
              if (lpos_def) pmfdxt(jl,jk,jt)=pmfd(jl,jk)*pxtd(jl,jk-1,jt) ! op_ck_20031001
! mz_ht_20040318- end of update for pos_def

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
! mz_ht_20040318+   for channel output
     do jl=1,kproma
! op_mm_20140131
! for 2d fields 
!!        d_entr(jl,jk,jrow)      = zdmfen(jl)
!!        d_detr(jl,jk,jrow)      = zdmfde(jl)
!!        massfd_draf(jl,jk,jrow) = pmfd(jl,jk)
        d_entr_2d(jl,jk)      = zdmfen(jl)
        d_detr_2d(jl,jk)      = zdmfde(jl)
        massfd_draf_2d(jl,jk) = pmfd(jl,jk)
     enddo
! mz_ht_20031712-
180 END DO
!

! op_m_20140226+  
! Calculate the maximum possible convective gust
!
  DO jl=1,kproma
    pvddraf(jl) = SQRT(0.2_dp*MAX(zvbuo(jl),0.0_dp) + pvddraf(jl))
    pvddraf(jl) = MIN(pvddraf(jl),30.0_dp)   ! But do not allow convective
                                               ! gusts higher than 30 m/s
  ENDDO

! op_m_20140226-

!
  RETURN
END SUBROUTINE tiedtke_cuddraf
!=============================================================================

SUBROUTINE tiedtke_cudlfs(   kproma, kbdim, klev, klevp1,              &
           ptenh,    pqenh,    puen,     pven,                         &
           ktrac,                                                      &
           pxtenh,   pxtu,     pxtd,     pmfdxt,                       &
           pgeoh,    paphp1,                                           &
           ptu,      pqu,      puu,      pvu,                          &
           ldcum,    kcbot,    kctop,    pmfub,    prfl,               &
           ptd,      pqd,      pud,      pvd,                          &
           pmfd,     pmfds,    pmfdq,    pdmfdp,                       &
           pcpcu,                                                      &
           kdtop,    lddraf)
!
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


!mz_ht_20040317+ added variables
INTEGER, INTENT(IN):: kproma, kbdim, klev, ktrac, klevp1
!mz_ht_20040317-
!
REAL(dp) :: ptenh(kbdim,klev),       pqenh(kbdim,klev),                   &
            puen(kbdim,klev),        pven(kbdim,klev),                    &
            pgeoh(kbdim,klev),       paphp1(kbdim,klevp1),                &
            ptu(kbdim,klev),         pqu(kbdim,klev),                     &
            puu(kbdim,klev),         pvu(kbdim,klev),                     &
            pmfub(kbdim),            prfl(kbdim)
!
REAL(dp) :: ptd(kbdim,klev),         pqd(kbdim,klev),                     &
            pud(kbdim,klev),         pvd(kbdim,klev),                     &
            pmfd(kbdim,klev),        pmfds(kbdim,klev),                   &
            pmfdq(kbdim,klev),       pdmfdp(kbdim,klev)
REAL(dp) :: pcpcu(kbdim,klev)
INTEGER ::  kcbot(kbdim),            kctop(kbdim),                        &
            kdtop(kbdim)
LOGICAL ::  ldcum(kbdim),            lddraf(kbdim)
!
REAL(dp) :: ztenwb(kbdim,klev),      zqenwb(kbdim,klev),                  &
            zcond(kbdim)
REAL(dp) :: zph(kbdim)
LOGICAL  :: llo2(kbdim)
REAL(dp) ::  pxtenh(kbdim,klev,ktrac),pxtu(kbdim,klev,ktrac),              &
             pxtd(kbdim,klev,ktrac),  pmfdxt(kbdim,klev,ktrac)
LOGICAL  :: llo3(kbdim)
!
!mz_ht_20040317+ local variables added
INTEGER  :: jl, jk, jt, is, ik, icall, ke
REAL(dp) :: zttest, zqtest, zbuo, zmftop
!mz_ht_20040317-
!----------------------------------------------------------------------
!
!     1.           SET DEFAULT VALUES FOR DOWNDRAFTS
!                  ---------------------------------
!
!! 100 CONTINUE
  DO 110 jl=1,kproma
     lddraf(jl)=.FALSE.
     kdtop(jl)=klevp1
110 END DO
!
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
!! 200 CONTINUE
!
  ke=klev-3
  DO 290 jk=3,ke
!
!
!     2.1          CALCULATE WET-BULB TEMPERATURE AND MOISTURE
!                  FOR ENVIRONMENTAL AIR IN *CUADJTQ*
!                  -------------------------------------------
!
!! 210  CONTINUE
     is=0
     DO 212 jl=1,kproma
        ztenwb(jl,jk)=ptenh(jl,jk)
        zqenwb(jl,jk)=pqenh(jl,jk)
        zph(jl)=paphp1(jl,jk)
        llo2(jl)=ldcum(jl).AND.prfl(jl).GT.0._dp.AND..NOT.lddraf(jl)      &
                          .AND.(jk.LT.kcbot(jl).AND.jk.GT.kctop(jl))
        is=is+MERGE(1,0,llo2(jl))
212  END DO
!     IF(is.EQ.0) go to 290 ! mz_ht_20071217
!
     ik=jk
     icall=2
     CALL tiedtke_cuadjtq(kproma, kbdim, klev, ik,                             &
                  zph,      ztenwb,   zqenwb,   llo2,     icall)
     IF (lookupoverflow) RETURN
!
!
!     2.2          DO MIXING OF CUMULUS AND ENVIRONMENTAL AIR
!                  AND CHECK FOR NEGATIVE BUOYANCY.
!                  THEN SET VALUES FOR DOWNDRAFT AT LFS.
!                  ----------------------------------------
!
!! 220  CONTINUE
!DIR$ IVDEP
!OCL NOVREC
     DO 222 jl=1,kproma
        llo3(jl)=.FALSE.
        IF(llo2(jl)) THEN
           zttest=0.5_dp*(ptu(jl,jk)+ztenwb(jl,jk))
           zqtest=0.5_dp*(pqu(jl,jk)+zqenwb(jl,jk))
           zbuo=zttest*(1._dp+vtmpc1*zqtest)-                             &
                         ptenh(jl,jk)*(1._dp+vtmpc1*pqenh(jl,jk))
           zcond(jl)=pqenh(jl,jk)-zqenwb(jl,jk)
           zmftop=-cmfdeps*pmfub(jl)
           IF(zbuo.LT.0._dp.AND.prfl(jl).GT.10._dp*zmftop*zcond(jl)) THEN
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
     DO 2224 jt=1,ktrac
        DO 2222 jl=1,kproma
           IF(llo3(jl)) THEN
              pxtd(jl,jk,jt)=0.5_dp*(pxtu(jl,jk,jt)+pxtenh(jl,jk,jt))
              pmfdxt(jl,jk,jt)=pmfd(jl,jk)*pxtd(jl,jk,jt)
           ENDIF
2222    END DO
2224 END DO
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
END SUBROUTINE tiedtke_cudlfs
!==============================================================================

SUBROUTINE tiedtke_cudtdq(kproma, kbdim, klev, klevp1, ktopm2, ldcum, ktrac,  &
                          jrow,     paphp1,   pten,     ptte,     pqte,       &
                          pxtte,    pxtec,    pmfuxt,   pmfdxt,               &
                          pmfus,    pmfds,    pmfuq,    pmfdq,                &
                          pmful,    pdmfup,   pdmfdp,   plude,                &
                          pdpmel,   prfl,     psfl,                           &
                          pcpen,    pqtec,    pqude,                          &
                          prsfc,    pssfc,    paprc,    paprs,                &
                          delta_time, ptracconv, ztmst, pxten,  pdtke_con,    &
                          pqenh, ptenh) ! op_mm_20140327 added pdtke_con, pqenh and ptenh
   USE messy_main_constants_mem,        ONLY:rd, cp_air, rv  ! op_mm_20140227



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

!mz_ht_20040317+ added variables
INTEGER, INTENT(IN):: kproma, kbdim, klev, ktrac, klevp1, jrow
INTEGER  :: ptracconv(ktrac), jl, jk, jt, ktopm2
REAL(dp), INTENT(IN) :: delta_time, ztmst
!mz_ht_20040317-!
LOGICAL  llo1
!
REAL(dp) ::  ptte(kbdim,klev),        pqte(kbdim,klev),                 &
             pten(kbdim,klev),        paphp1(kbdim,klevp1),             &
             paprc(kbdim),            paprs(kbdim),                     &
             prsfc(kbdim),            pssfc(kbdim)
REAL(dp) ::  pmfus(kbdim,klev),       pmfds(kbdim,klev),                &
             pmfuq(kbdim,klev),       pmfdq(kbdim,klev),                &
             pmful(kbdim,klev),       plude(kbdim,klev),                &
             pdmfup(kbdim,klev),      pdmfdp(kbdim,klev),               &
             pqtec(kbdim,klev),       pqude(kbdim,klev),                &
             pxtec(kbdim,klev),       prfl(kbdim)
REAL(dp) ::    pdpmel(kbdim,klev),      psfl(kbdim)
REAL(dp) ::    pcpen(kbdim,klev)
! op_mm_20140327+
! for tke calculation
REAL(dp) ::    pdtke_con(kbdim,klev),pqenh(kbdim, klev),ptenh(kbdim, klev)
REAL(dp) ::    pcvfl_s, pcvfl_q 
! op_mm_20140327-
LOGICAL  :: ldcum(kbdim)
!
REAL(dp) ::    zmelt(kbdim)
REAL(dp) ::    zsheat(kbdim)
REAL(dp) ::    pxtte(kbdim,klev,ktrac), pmfuxt(kbdim,klev,ktrac),         &
               pmfdxt(kbdim,klev,ktrac)
!
REAL(dp) ::    zrcpm ! reciprocal value of specific heat of moist air

REAL(dp) ::    zalpha(kbdim,klev,ktrac),zgamma(kbdim,ktrac),     & ! op_ck_20031001
               ztend(kbdim,klev,ktrac), ztenm(kbdim,klev,ktrac), & ! op_ck_20031001
               pxten(kbdim,klev,ktrac)                             ! op_ck_2003100
! mz_ht_20040317+    added variables
REAL(dp) :: zdiagt, zalv, zdtdt, zdqdt, zdxtdt, zeps1 

! mz_ht_20040317-
!----------------------------------------------------------------------
!
!*    1.0          SPECIFY PARAMETERS
!                  ------------------
!
!! 100 CONTINUE
  zeps1=(1._dp-1.E-10_dp)   ! op_ck_20031001
  zdiagt=delta_time
!
!
!----------------------------------------------------------------------
!
!*    2.0          INCREMENTATION OF T AND Q TENDENCIES
!                  ------------------------------------
!
!! 200 CONTINUE
  DO 210 jl=1,kproma
     zmelt(jl)=0._dp
     zsheat(jl)=0._dp
210 END DO
!

if (lpos_def) then
! op_ck_20031001+
  DO 2100 jt=1,ktrac
     DO 2101 jl=1,kproma
        ztend(jl,1,jt)=(pmfuxt(jl,1,jt)+pmfdxt(jl,1,jt))             &
                      *ztmst*g/(paphp1(jl,2)-paphp1(jl,1))
        ztenm(jl,1,jt)=max(ztend(jl,1,jt),pxten(jl,1,jt)*zeps1)
        ztend(jl,klev,jt)=(pmfuxt(jl,klev,jt)+pmfdxt(jl,klev,jt))    &
                         *ztmst*g/(paphp1(jl,klev+1)-paphp1(jl,klev))
        ztenm(jl,klev,jt)=max(ztend(jl,klev,jt)                      &
                             ,pxten(jl,klev,jt)*zeps1)
        zgamma(jl,jt)=999._dp
2101    ENDDO
2100 ENDDO
!
  DO 2110 jt=1,ktrac
     DO 2111 jk=2,klev-1
        DO 2112 jl=1,kproma
           ztend(jl,jk,jt)=(pmfuxt(jl,jk,jt)-pmfuxt(jl,jk+1,jt)      &
                          +pmfdxt(jl,jk,jt)-pmfdxt(jl,jk+1,jt))      &
                          *ztmst*g/(paphp1(jl,jk+1)-paphp1(jl,jk))
           ztenm(jl,jk,jt)=max(ztend(jl,jk,jt),pxten(jl,jk,jt)*zeps1)
2112       ENDDO
2111    ENDDO
2110 ENDDO
!
  DO 2120 jt=1,ktrac
     DO 2121 jk=1,klev
        DO 2122 jl=1,kproma
!!          IF (pxten(jl,jk,jt).le.0._dp.and.ztenm(jl,jk,jt).eq.0._dp) THEN
          IF (pxten(jl,jk,jt).le.0._dp.and.ABS(ztenm(jl,jk,jt)).lt.TINY(0._dp)) THEN
              zalpha(jl,jk,jt)=1._dp
           ELSE
              zalpha(jl,jk,jt)=(pxten(jl,jk,jt)*zeps1)/(ztenm(jl,jk,jt))
           ENDIF
           zgamma(jl,jt)=min(zgamma(jl,jt),zalpha(jl,jk,jt))
           zgamma(jl,jt)=min(zgamma(jl,jt),1._dp)
           zgamma(jl,jt)=max(zgamma(jl,jt),0._dp)
2122       ENDDO 
2121    ENDDO  
2120 ENDDO
! op_ck_20031001-
ENDIF
  DO 250 jk=ktopm2,klev
!
     IF(jk.LT.klev) THEN
        DO 220 jl=1,kproma
           IF(ldcum(jl)) THEN
              llo1=(pten(jl,jk)-tmelt).GT.0._dp
              zalv=MERGE(alv,als,llo1)
              zrcpm=1._dp/pcpen(jl,jk)
              zdtdt=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*zrcpm*        &
                                  (pmfus(jl,jk+1)-pmfus(jl,jk)+       &
                                   pmfds(jl,jk+1)-pmfds(jl,jk)-       &
                                   alf*pdpmel(jl,jk)-                 &
                                   zalv*(pmful(jl,jk+1)-pmful(jl,jk)- &
                                   plude(jl,jk)-                      &
                                  (pdmfup(jl,jk)+pdmfdp(jl,jk))))
              ptte(jl,jk)=ptte(jl,jk)+zdtdt
     !!         conv_tte(jl,jk,jrow) = zdtdt ! op_mm_20140327
              conv_tte_2d(jl,jk) = zdtdt
! mim_sb_20090917+
              IF (l_lgmc_diag) THEN
! op_mm_20140327+
! 2D fields  
 !!                conv_tte_up_cond(jl,jk,jrow) =                          &
 !!                     -zalv*(pmful(jl,jk+1)-pmful(jl,jk) &
 !!                     -plude(jl,jk)-                     &
 !!                     pdmfup(jl,jk))*g/(paphp1(jl,jk+1)-paphp1(jl,jk))*zrcpm
 
                 conv_tte_up_cond_2d(jl,jk) =                          &
                      -zalv*(pmful(jl,jk+1)-pmful(jl,jk) &
                      -plude(jl,jk)-                     &
                      pdmfup(jl,jk))*g/(paphp1(jl,jk+1)-paphp1(jl,jk))*zrcpm
! op_mm_20140327-

!!                 conv_tte_do_melt(jl,jk,jrow) = &
!!                      (g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*zrcpm*alf*(-pdpmel(jl,jk))

                 conv_tte_do_melt_2d(jl,jk) = &
                      (g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*zrcpm*alf*(-pdpmel(jl,jk))
! op_mm_20140327+
! for 2d fields                 


!                 conv_tte_do_verd(jl,jk,jrow) = &
!                      (g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*zrcpm*zalv*(pdmfdp(jl,jk))

                 conv_tte_do_verd_2d(jl,jk) = &
                      (g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*zrcpm*zalv*(pdmfdp(jl,jk))

!!                 conv_tte_up_freeze(jl,jk,jrow) =                        &  
!!                      conv_tte_do_melt(jl,jk,jrow)+&
!!                      conv_tte_do_verd(jl,jk,jrow)+conv_tte_up_cond(jl,jk,jrow)
 !!                conv_tte_ev(jl,jk,jrow) = g/(paphp1(jl,jk+1)-paphp1(jl,jk))*zrcpm*        &
!!                      (pmfus(jl,jk+1)-pmfus(jl,jk)+                           &
!!                      pmfds(jl,jk+1)-pmfds(jl,jk))

                 conv_tte_up_freeze_2d(jl,jk) =                        &  
                      conv_tte_do_melt_2d(jl,jk)+&
                      conv_tte_do_verd_2d(jl,jk)+conv_tte_up_cond_2d(jl,jk)

                 conv_tte_ev_2d(jl,jk) = g/(paphp1(jl,jk+1)-paphp1(jl,jk))*zrcpm*        &
                      (pmfus(jl,jk+1)-pmfus(jl,jk)+                           &
                      pmfds(jl,jk+1)-pmfds(jl,jk))
! op_mm_20140327-

              END IF
! mim_sb_20090917-
              zdqdt=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*              &
                                  (pmfuq(jl,jk+1)-pmfuq(jl,jk)+       &
                                   pmfdq(jl,jk+1)-pmfdq(jl,jk)+       &
                                   pmful(jl,jk+1)-pmful(jl,jk)-       &
                                   plude(jl,jk)-                      &
                                  (pdmfup(jl,jk)+pdmfdp(jl,jk)))
!!              conv_qte(jl,jk,jrow) = zdqdt      ! op_mm_20140327 2d fields 
              conv_qte_2d(jl,jk) = zdqdt
              pqte(jl,jk)=pqte(jl,jk)+zdqdt
              pxtec(jl,jk)=(g/(paphp1(jl,jk+1)-                       &
                            paphp1(jl,jk)))*plude(jl,jk)
              pqtec(jl,jk)=(g/(paphp1(jl,jk+1)-                       &
                            paphp1(jl,jk)))*pqude(jl,jk)
              zsheat(jl)=zsheat(jl)+zalv*(pdmfup(jl,jk)+pdmfdp(jl,jk))
              zmelt(jl)=zmelt(jl)+pdpmel(jl,jk)
           END IF
220     END DO
!mz_ht_20040317+
!        if (MAXVAL(PTRACCONV).gt.0) THEN
!mz_ht_20040317-
           DO 2204 jt=1,ktrac
!mz_ht_20040317+
              IF (PTRACCONV(jt).eq.1) THEN
!              IF (ti_gp(jt)%tp%proc%nconvect == ON) THEN
!mz_ht_20040317-
                DO 2202 jl=1,kproma
                   IF(ldcum(jl)) THEN
                      if (lpos_def) then
                         zdxtdt=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))        &
                           *(pmfuxt(jl,jk+1,jt)-pmfuxt(jl,jk,jt) &
                           +pmfdxt(jl,jk+1,jt)-pmfdxt(jl,jk,jt))&
                           *zgamma(jl,jt) ! op_ck_20031001
                      ELSE
                         zdxtdt=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))        &
                           *(pmfuxt(jl,jk+1,jt)-pmfuxt(jl,jk,jt) &
                           +pmfdxt(jl,jk+1,jt)-pmfdxt(jl,jk,jt))
                      ENDIF
                      IF(.NOT.CVTRANS) &          ! mz_ht_20031512
                           pxtte(jl,jk,jt)=pxtte(jl,jk,jt)+zdxtdt
                   ENDIF
2202            END DO
             ENDIF
2204       END DO
!        ENDIF
!
!
     ELSE
        DO 230 jl=1,kproma
           IF(ldcum(jl)) THEN
              llo1=(pten(jl,jk)-tmelt).GT.0._dp
              zalv=MERGE(alv,als,llo1)
              zrcpm=1._dp/pcpen(jl,jk)
              zdtdt=-(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*zrcpm*       &
                     (pmfus(jl,jk)+pmfds(jl,jk)+alf*pdpmel(jl,jk)-    &
                                 zalv*(pmful(jl,jk)+pdmfup(jl,jk)     &
                                +pdmfdp(jl,jk)+plude(jl,jk)))
              ptte(jl,jk)=ptte(jl,jk)+zdtdt
!!              conv_tte(jl,jk,jrow) = zdtdt     ! op_mm_20140327 for 2d fields 
              conv_tte_2d(jl,jk) = zdtdt
! mim_sb_20090917+
              IF (l_lgmc_diag) THEN
! op_mm_20140327+
! for 2d fields                  
                
!!                 conv_tte_up_cond(jl,jk,jrow) =                          &
!!                      -zalv*(pmful(jl,jk)-               &
!!                      plude(jl,jk)-                      &
!!                      pdmfup(jl,jk))*g/(paphp1(jl,jk+1)-paphp1(jl,jk))*zrcpm
 
                conv_tte_up_cond_2d(jl,jk) =                          &
                      -zalv*(pmful(jl,jk)-               &
                      plude(jl,jk)-                      &
                      pdmfup(jl,jk))*g/(paphp1(jl,jk+1)-paphp1(jl,jk))*zrcpm


!!                 conv_tte_do_melt(jl,jk,jrow) = &
!!                      (g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*zrcpm*alf*(-pdpmel(jl,jk))

            conv_tte_do_melt_2d(jl,jk) = &
                      (g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*zrcpm*alf*(-pdpmel(jl,jk))

!!                 conv_tte_do_verd(jl,jk,jrow) = &
!!                      (g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*zrcpm*zalv*(pdmfdp(jl,jk))

            conv_tte_do_verd_2d(jl,jk) = &
                 (g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*zrcpm*zalv*(pdmfdp(jl,jk))


 !!                conv_tte_up_freeze(jl,jk,jrow) =                        & 
 !!                     conv_tte_do_melt(jl,jk,jrow)+conv_tte_do_verd(jl,jk,jrow)+ &
 !!                     conv_tte_up_cond(jl,jk,jrow)
 !!                conv_tte_ev(jl,jk,jrow) =-(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*zrcpm*       &
  !!                    (pmfus(jl,jk)+pmfds(jl,jk))

                 conv_tte_up_freeze_2d(jl,jk) =                        & 
                      conv_tte_do_melt_2d(jl,jk)+conv_tte_do_verd_2d(jl,jk)+ &
                      conv_tte_up_cond_2d(jl,jk)
                 conv_tte_ev_2d(jl,jk) =-(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*zrcpm*       &
                      (pmfus(jl,jk)+pmfds(jl,jk))
! op_mm_20140327-

              ENDIF
! mim_sb_20090917-
              zdqdt=-(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*             &
                        (pmfuq(jl,jk)+pmfdq(jl,jk)+plude(jl,jk)+     &
                        (pmful(jl,jk)+pdmfup(jl,jk)+pdmfdp(jl,jk)))
  !!            conv_qte(jl,jk,jrow) = zdqdt     ! op_mm_20140327 for 2d fields 
              conv_qte_2d(jl,jk) = zdqdt
              pqte(jl,jk)=pqte(jl,jk)+zdqdt
              pxtec(jl,jk)=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))        &
                           *plude(jl,jk)
              pqtec(jl,jk)=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))        &
                           *pqude(jl,jk)
              zsheat(jl)=zsheat(jl)+zalv*(pdmfup(jl,jk)+pdmfdp(jl,jk))
              zmelt(jl)=zmelt(jl)+pdpmel(jl,jk)
           END IF
230     END DO
!
!mz_ht_20040317+
 !       if (MAXVAL(PTRACCONV).gt.0) THEN
!mz_ht_20040317-
           DO 2304 jt=1,ktrac
!mz_ht_20040317+
              IF (PTRACCONV(jt).eq.1) THEN
 !             IF (ti_gp(jt)%tp%proc%nconvect == ON) THEN
!mz_ht_20040317-
                DO 2302 jl=1,kproma
                   IF(ldcum(jl)) THEN
                      if (lpos_def) then
                         zdxtdt=-(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))        &
                           *(pmfuxt(jl,jk,jt)+pmfdxt(jl,jk,jt))&
                           *zgamma(jl,jt) ! op_ck_20031001
                      ELSE
                         zdxtdt=-(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))        &
                           *(pmfuxt(jl,jk,jt)+pmfdxt(jl,jk,jt))
                      ENDIF
                      IF(.NOT.CVTRANS) &          ! mz_ht_20031512     
                           pxtte(jl,jk,jt)=pxtte(jl,jk,jt)+zdxtdt
                   ENDIF
2302            END DO
              END IF
2304       END DO
      !  ENDIF
!
     END IF
     
! mz_ht_20031030+
     WHERE (pten(1:kproma,jk).gt.tmelt)
! op_mm_20140327+
! for 2d fields 
  !!     cv_precnew(1:kproma,jk,jrow) = pdmfup(1:kproma,jk)+pdmfdp(1:kproma,jk) 
        cv_precnew_2d(1:kproma,jk) = pdmfup(1:kproma,jk)+pdmfdp(1:kproma,jk)
     ELSEWHERE
!!       cv_snownew(1:kproma,jk,jrow) = pdmfup(1:kproma,jk)+pdmfdp(1:kproma,jk)
        cv_snownew_2d(1:kproma,jk) = pdmfup(1:kproma,jk)+pdmfdp(1:kproma,jk)
     ENDWHERE
 ! op_mm_20140327-
! mz_ht_20031030-

   
250 END DO
  ! op_mm_20140227+ 
  ! calculation of tke
  ! copied from src_conv_tiedtke (COSMO) 
  ! zmfus -> pmfus ; zmfds -> pmfds
  ! zmfuq -> pmfuq ; zmfdq -> pmfdq 
  ! zqenh -> pqenh ; ztenh -> ptenh
  ! pph -> paphp1
  !set whole field to zero (as tke only computed when tdcum(jl))
  pdtke_con(1:kproma,:)=0.0_dp 
  
  ! for all model half levels (except the lower model boundary)
  DO jk = ktopm2,klev 
     DO jl = 1, kproma
        IF (ldcum(jl)) THEN ! for convective grid points only
           pcvfl_s  =  pmfus (jl,jk) + pmfds (jl,jk)
           pcvfl_q  =  pmfuq (jl,jk) + pmfdq (jl,jk)
           ! um_ak_20140604+ workaround in ECHAM paphp1 is 0. for jk = 1
           IF (jk == 1) THEN 
              pdtke_con(jl,jk) = 0.0_dp
           ELSE
           ! um_ak_20140604+   
           pdtke_con(jl,jk) = MAX( 0.0_dp, g*rd/paphp1(jl,jk) * &
                ( pcvfl_s/cp_air + &
                ptenh(jl,jk)*pcvfl_q*(rv/rd-1.0_dp)/ &
                (1.0_dp-(rv/rd-1.0_dp)*pqenh(jl,jk)) ) )
           END IF  ! um_ak_20140604
        END IF
     END DO
  END DO
  ! op_mm_20140227-

!
!
!---------------------------------------------------------------------
!
!      3.          UPDATE SURFACE FIELDS
!                  ---------------------
!
!! 300 CONTINUE
  DO 310 jl=1,kproma
     prsfc(jl)=prfl(jl)
     pssfc(jl)=psfl(jl)
     paprc(jl)=paprc(jl)+zdiagt*(prfl(jl)+psfl(jl))
     paprs(jl)=paprs(jl)+zdiagt*psfl(jl)
310 END DO
!
  RETURN
END SUBROUTINE tiedtke_cudtdq
!==============================================================================

SUBROUTINE tiedtke_cududv(   kproma,   kbdim,    klev,     klevp1,     &
           ktopm2,   ktype,    kcbot,    paphp1,   ldcum,              &
           puen,     pven,     pvom,     pvol,                         &
           puu,      pud,      pvu,      pvd,                          &
           pmfu,     pmfd)
!
!
!**** *CUDUDV* - UPDATES U AND V TENDENCIES,
!                DOES GLOBAL DIAGNOSTIC OF DISSIPATION
!
!          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
!
!**   INTERFACE.
!     ----------
!
!          *CUDUDV* IS CALLED FROM *CUMASTR*
!
!mz_ht_20040317+ added variables
INTEGER, INTENT(IN):: kproma, kbdim, klev, klevp1
!mz_ht_20040317-
REAL(dp) ::    puen(kbdim,klev),        pven(kbdim,klev),                    &
           pvol(kbdim,klev),        pvom(kbdim,klev),                    &
           paphp1(kbdim,klevp1)
REAL(dp) ::    puu(kbdim,klev),         pud(kbdim,klev),                     &
           pvu(kbdim,klev),         pvd(kbdim,klev),                     &
           pmfu(kbdim,klev),        pmfd(kbdim,klev)
INTEGER :: ktype(kbdim),            kcbot(kbdim)
LOGICAL :: ldcum(kbdim)
!
REAL(dp) ::    zmfuu(kbdim,klev),       zmfdu(kbdim,klev),                   &
           zmfuv(kbdim,klev),       zmfdv(kbdim,klev)
!mz_ht_20040317+ added variables

REAL(dp) :: zzp, zdudt,zdvdt 
INTEGER  :: jl, jk, ik, ikb, ktopm2
!
!----------------------------------------------------------------------
!
!*    1.0          CALCULATE FLUXES AND UPDATE U AND V TENDENCIES
!                  ----------------------------------------------
!

!! 100 CONTINUE
  IF(ktopm2.EQ.1) THEN
    DO 107 jk=2,klev
       ik=jk-1
       DO 106 jl=1,kproma
          IF(ldcum(jl)) THEN
             zmfuu(jl,jk)=pmfu(jl,jk)*(puu(jl,jk)-puen(jl,ik))
             zmfuv(jl,jk)=pmfu(jl,jk)*(pvu(jl,jk)-pven(jl,ik))
             zmfdu(jl,jk)=pmfd(jl,jk)*(pud(jl,jk)-puen(jl,ik))
             zmfdv(jl,jk)=pmfd(jl,jk)*(pvd(jl,jk)-pven(jl,ik))
          END IF
106    END DO
107  END DO
    DO 105 jl=1,kproma
      IF(ldcum(jl)) THEN
        zmfuu(jl,1)=zmfuu(jl,2)
        zmfuv(jl,1)=zmfuv(jl,2)
        zmfdu(jl,1)=zmfdu(jl,2)
        zmfdv(jl,1)=zmfdv(jl,2)
      END IF
105 END DO
  ELSE
    DO 120 jk=ktopm2,klev
       ik=jk-1
       DO 110 jl=1,kproma
          IF(ldcum(jl)) THEN
             zmfuu(jl,jk)=pmfu(jl,jk)*(puu(jl,jk)-puen(jl,ik))
             zmfuv(jl,jk)=pmfu(jl,jk)*(pvu(jl,jk)-pven(jl,ik))
             zmfdu(jl,jk)=pmfd(jl,jk)*(pud(jl,jk)-puen(jl,ik))
             zmfdv(jl,jk)=pmfd(jl,jk)*(pvd(jl,jk)-pven(jl,ik))
          END IF
110    END DO
120  END DO
  END IF

  DO 140 jk=ktopm2,klev
!DIR$ IVDEP
!OCL NOVREC
     DO 130 jl=1,kproma
        IF(ldcum(jl).AND.jk.GT.kcbot(jl)) THEN
           ikb=kcbot(jl)
           zzp=((paphp1(jl,klevp1)-paphp1(jl,jk))/                     &
                         (paphp1(jl,klevp1)-paphp1(jl,ikb)))
           zzp=MERGE(zzp**2,zzp,ktype(jl).EQ.3)
           zmfuu(jl,jk)=zmfuu(jl,ikb)*zzp
           zmfuv(jl,jk)=zmfuv(jl,ikb)*zzp
           zmfdu(jl,jk)=zmfdu(jl,ikb)*zzp
           zmfdv(jl,jk)=zmfdv(jl,ikb)*zzp
        END IF
130  END DO
140 END DO
!
  DO 190 jk=ktopm2,klev
!
     IF(jk.LT.klev) THEN
        DO 160 jl=1,kproma
           IF(ldcum(jl)) THEN
              zdudt=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*               &
                          (zmfuu(jl,jk+1)-zmfuu(jl,jk)+                &
                           zmfdu(jl,jk+1)-zmfdu(jl,jk))
              zdvdt=(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*               &
                          (zmfuv(jl,jk+1)-zmfuv(jl,jk)+                &
                           zmfdv(jl,jk+1)-zmfdv(jl,jk))
              pvom(jl,jk)=pvom(jl,jk)+zdudt
              pvol(jl,jk)=pvol(jl,jk)+zdvdt
           END IF
160     END DO
!
     ELSE
        DO 170 jl=1,kproma
           IF(ldcum(jl)) THEN
              zdudt=-(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*              &
                           (zmfuu(jl,jk)+zmfdu(jl,jk))
              zdvdt=-(g/(paphp1(jl,jk+1)-paphp1(jl,jk)))*              &
                           (zmfuv(jl,jk)+zmfdv(jl,jk))
              pvom(jl,jk)=pvom(jl,jk)+zdudt
              pvol(jl,jk)=pvol(jl,jk)+zdvdt
           END IF
170     END DO
     END IF
!
190 END DO
!
!
  RETURN
END SUBROUTINE tiedtke_cududv
!==============================================================================
SUBROUTINE tiedtke_cuentr(kproma, kbdim, klev, klevp1, kk,             &
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

!
!mz_ht_20040317+ added variables
INTEGER, INTENT(IN) :: kbdim, kproma, klev, kk, klevp1
!mz_ht_20040317-
REAL(dp) :: ptenh(kbdim,klev),       pqenh(kbdim,klev),                  &
            paphp1(kbdim,klevp1),                                        &
            pmfu(kbdim,klev),        pqte(kbdim,klev),                   &
            pentr(kbdim),            ppbase(kbdim)
REAL(dp) :: podetr(kbdim,klev)
REAL(dp) :: pgeoh (kbdim,klev)
INTEGER  :: khmin (kbdim)
INTEGER  :: klwmin(kbdim),           ktype(kbdim),                       &
            kcbot(kbdim),            kctop0(kbdim)
LOGICAL  :: ldcum(kbdim)
!
REAL(dp) :: pdmfen(kbdim),           pdmfde(kbdim)
!
LOGICAL  ::  llo1,llo2
!
!mz_ht_20040317+ local variables added
INTEGER  :: ikt, ikh, jl, iklwmin
REAL(dp) :: arg, zrg, zrrho, zdprho, zpmid, zentr, zentest, zzmzk, ztmzk, zorgde
!mz_ht_20040317-
!----------------------------------------------------------------------
!
!*    1.           CALCULATE ENTRAINMENT AND DETRAINMENT RATES
!                  -------------------------------------------
!
!! 100 CONTINUE
!
!
!*    1.1          SPECIFY ENTRAINMENT RATES FOR SHALLOW CLOUDS
!                  --------------------------------------------
!
!! 110 CONTINUE
!
!
!*    1.2          SPECIFY ENTRAINMENT RATES FOR DEEP CLOUDS
!                  -----------------------------------------
!
!! 120 CONTINUE
  zrg=1._dp/g
  DO 125 jl=1,kproma
     ppbase(jl)=paphp1(jl,kcbot(jl))
     zrrho=(rd*ptenh(jl,kk+1)*(1._dp+vtmpc1*pqenh(jl,kk+1)))              &
                    /paphp1(jl,kk+1)
     zdprho=(paphp1(jl,kk+1)-paphp1(jl,kk))*zrg
     zpmid=0.5_dp*(ppbase(jl)+paphp1(jl,kctop0(jl)))
     zentr=pentr(jl)*pmfu(jl,kk+1)*zdprho*zrrho
     llo1=kk.LT.kcbot(jl).AND.ldcum(jl)
     pdmfde(jl)=MERGE(zentr,0._dp,llo1)
     llo2=llo1.AND.ktype(jl).EQ.2.AND.                                 &
                   (ppbase(jl)-paphp1(jl,kk).LT.0.2e5_dp.OR.              &
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
     podetr(jl,kk)=0._dp
     IF(llo2.AND.kk.LE.khmin(jl).AND.kk.GE.kctop0(jl)) THEN
        ikt=kctop0(jl)
        ikh=khmin(jl)
        IF(ikh.GT.ikt) THEN
           zzmzk  =-(pgeoh(jl,ikh)-pgeoh(jl,kk))*zrg
           ztmzk  =-(pgeoh(jl,ikh)-pgeoh(jl,ikt))*zrg
           arg  =3.1415_dp*(zzmzk/ztmzk)*0.5_dp
           zorgde=TAN(arg)*3.1415_dp*0.5_dp/ztmzk
           zdprho=(paphp1(jl,kk+1)-paphp1(jl,kk))*(zrg*zrrho)
           podetr(jl,kk)=MIN(zorgde,centrmax)*pmfu(jl,kk+1)*zdprho
        ENDIF
     ENDIF
125 END DO
!
  RETURN
END SUBROUTINE tiedtke_cuentr

!==============================================================================

SUBROUTINE tiedtke_cuentrt(  kproma, kbdim, klev, klevp1, kk,          &
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


!mz_ht_20040317+ added variables
INTEGER, INTENT(IN) :: kbdim, kproma, klev, kk, klevp1
!mz_ht_20040317-
!
REAL(dp) :: ptenh(kbdim,klev),       pqenh(kbdim,klev),                  &
            paphp1(kbdim,klevp1),                                        &
            pmfu(kbdim,klev),        pqte(kbdim,klev),                   &
            pentr(kbdim),            ppbase(kbdim)
INTEGER  :: klwmin(kbdim),           ktype(kbdim),                       &
            kcbot(kbdim),            kctop0(kbdim)
LOGICAL  :: ldcum(kbdim)
!
REAL(dp) :: pdmfen(kbdim),           pdmfde(kbdim)
!
LOGICAL  ::  llo1,llo2

!mz_ht_20040317+ local variables added
INTEGER  :: jl, iklwmin
REAL(dp) :: zrg, zrrho, zdprho, zpmid, zentr, zentest
!mz_ht_20040317-
!
!----------------------------------------------------------------------
!
!*    1.           CALCULATE ENTRAINMENT AND DETRAINMENT RATES
!                  -------------------------------------------
!
!! 100 CONTINUE
!
!
!*    1.1          SPECIFY ENTRAINMENT RATES FOR SHALLOW CLOUDS
!                  --------------------------------------------
!
!! 110 CONTINUE
!
!
!*    1.2          SPECIFY ENTRAINMENT RATES FOR DEEP CLOUDS
!                  -----------------------------------------
!
!! 120 CONTINUE
  zrg=1._dp/g
  DO 125 jl=1,kproma
     ppbase(jl)=paphp1(jl,kcbot(jl))
     zrrho=(rd*ptenh(jl,kk+1)*(1._dp+vtmpc1*pqenh(jl,kk+1)))              &
                    /paphp1(jl,kk+1)
     zdprho=(paphp1(jl,kk+1)-paphp1(jl,kk))*zrg
     zpmid=0.5_dp*(ppbase(jl)+paphp1(jl,kctop0(jl)))
     zentr=pentr(jl)*pmfu(jl,kk+1)*zdprho*zrrho
     llo1=kk.LT.kcbot(jl).AND.ldcum(jl)
     pdmfde(jl)=MERGE(zentr,0._dp,llo1)
     llo2=llo1.AND.ktype(jl).EQ.2.AND.                                 &
                   (ppbase(jl)-paphp1(jl,kk).LT.0.2e5_dp.OR.              &
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
END SUBROUTINE tiedtke_cuentrt
!==============================================================================
SUBROUTINE tiedtke_cuflx(kproma, kbdim, klev, klevp1, jrow,            &
           pqen,     pqsen,    ptenh,    pqenh,                        &
           ktrac,                                                      &
           pxtenh,   pmfuxt,   pmfdxt,                                 &
           paphp1,   pgeoh,                                            &
           kcbot,    kctop,    kdtop,                                  &
           ktype,    lddraf,   ldcum,                                  &
           pmfu,     pmfd,     pmfus,    pmfds,                        &
           pmfuq,    pmfdq,    pmful,                                  &
           pdmfup,   pdmfdp,   prfl,     prain,                        &
           pcpcu,                                                      &
           pcpen,                                                   &   ! mim_sb_20090917
           pten,     psfl,     pdpmel,   ktopm2,   ztmst,              &
           zpmfun,   pxten)  ! mz_ht_20040318+ needed for posdef update
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
!          scav_cvdep for convective scavenging

USE messy_scav_inter,     ONLY: scav_cvdep
USE messy_scav_mem,       ONLY: lscav_cv, kpp_rain_cv, aero_flx
          !
!mz_ht_20040317+ added variables
INTEGER,  INTENT(IN) :: kbdim, kproma, klev, klevp1, ktrac, jrow
REAL(dp), INTENT(IN) :: ztmst
!mz_ht_20040317-
!
REAL(dp) :: pqen(kbdim,klev),        pqsen(kbdim,klev),                     &
            ptenh(kbdim,klev),       pqenh(kbdim,klev),                     &
            paphp1(kbdim,klevp1),    pgeoh(kbdim,klev)
!
REAL(dp) :: pmfu(kbdim,klev),        pmfd(kbdim,klev),                      &
            pmfus(kbdim,klev),       pmfds(kbdim,klev),                     &
            pmfuq(kbdim,klev),       pmfdq(kbdim,klev),                     &
            pdmfup(kbdim,klev),      pdmfdp(kbdim,klev),                    &
            pmful(kbdim,klev),                                              &
            prfl(kbdim),             prain(kbdim)
REAL(dp) :: pten(kbdim,klev),        pdpmel(kbdim,klev),                    &
            psfl(kbdim)
REAL(dp) :: pcpcu(kbdim,klev)                                               
REAL(dp) :: pcpen(kbdim,klev)                                         ! mim_sb_20090917
REAL(dp) :: zpsubcl(kbdim)
INTEGER  :: kcbot(kbdim),            kctop(kbdim),                          &
            kdtop(kbdim),            ktype(kbdim)
LOGICAL  :: lddraf(kbdim),           ldcum(kbdim)
REAL(dp) :: pxtenh(kbdim,klev,ktrac),                                      &
            pmfuxt(kbdim,klev,ktrac),pmfdxt(kbdim,klev,ktrac)

!mz_ht_20040317+ added local variables
INTEGER  :: jl, jk, jt, ikb, ktopm2
REAL(dp) :: zcons1, zcons2, zcucov, ztmelp2, zzp,   zfac,  zsnmlt, zrfl,   &
            zrnew,  zrmin,  zrfln,  zdrfl,   zrsum, zdpevap
REAL(dp) :: zrhoa, zrainh(kbdim), zcover_conv(kbdim)
! variables for posdef update
LOGICAL  :: ldp(kbdim,klev) ! op_ck_20031001
REAL(dp) :: zmass(kbdim,klev), zpmfun(kbdim,klev),pxten(kbdim,klev,ktrac) ! op_ck_20031001
REAL(dp) :: zmfu(kbdim)
!mz_ht_20040317-
REAL(dp) ::    zrcpm ! reciprocal value of specific heat of moist air ! mim_sb_20090917
!
!*             SPECIFY CONSTANTS
!

  zcons1=cpd/(alf*g*ztmst)
  zcons2=1._dp/(g*ztmst)
  zcucov=0.05_dp
  ztmelp2=tmelt+2._dp

! mz_ht_20040318+
! values for lnox set to zero in each time step
!!$ cu_bot(:,jrow)      = 0._dp
!!$ cu_top(:,jrow)      = 0._dp
  cu_bot_1d(:)      = 0._dp          ! op_mm_20140521
  cu_top_1d(:)      = 0._dp

!!$  cu_freeze(:,jrow)   = 0._dp
  cu_freeze_1d(:)   = 0._dp          ! op_mm_20140521
!  cu_uvelo(:,:,jrow)  = 0._dp       ! op_mm_20140327 for 2d fields 
  cu_uvelo_2d(:,:)  = 0._dp

  ! mz_pj_20050615+
!!$  cu_bot_mid(:,jrow)      = 0._dp
!!$  cu_top_mid(:,jrow)      = 0._dp
!!$  cu_freeze_mid(:,jrow)   = 0._dp  
  cu_bot_mid_1d(:)      = 0._dp
  cu_top_mid_1d(:)      = 0._dp
  cu_freeze_mid_1d(:)   = 0._dp  

  ! mz_pj_20050615-
! mz_ht_20040318-
!
!
!
!
!*    1.0          DETERMINE FINAL CONVECTIVE FLUXES
!                  ---------------------------------
!
!! 100 CONTINUE
!  itop=klev
  DO 110 jl=1,kproma
!     itop=MIN(itop,kctop(jl))
     IF(.NOT.ldcum(jl).OR.kdtop(jl).LT.kctop(jl)) lddraf(jl)=.FALSE.
     IF(.NOT.ldcum(jl)) ktype(jl)=0
110 END DO
  ktopm2=1

! mz_ht_20040318+ update for posdef 
  if (lpos_def) then 
! op_ck_20031001+
     DO 1100 jk=1,klev
        DO 1101 jl=1,kproma
           zmass(jl,jk)=0._dp
           ldp(jl,jk)=.FALSE.
1101    ENDDO
1100 ENDDO
     !
     DO 1110 jt=1,ktrac
        DO 1111 jk=ktopm2,klev
           DO 1112 jl=1,kproma
              IF (ldcum(jl).and.jk.ge.kctop(jl)-1) THEN
                 zmass(jl,jk)=zpmfun(jl,jk)
                 IF (lddraf(jl).and.jk.ge.kdtop(jl)) THEN
                    zmass(jl,jk)=zpmfun(jl,jk)+pmfd(jl,jk)
                 ELSE
                    pmfdxt(jl,jk,jt)=0._dp
                 ENDIF
              ELSE
                 pmfuxt(jl,jk,jt)=0._dp
                 pmfdxt(jl,jk,jt)=0._dp
              ENDIF
1112       ENDDO
1111    ENDDO
1110 ENDDO
!
     DO 1120 jk=ktopm2,klev
        DO 1121 jl=1,kproma
           IF (zmass(jl,jk).GE.0._dp) ldp(jl,jk)=.TRUE.
1121    ENDDO
1120 ENDDO    
! op_ck_20031001-
  ENDIF    ! of lpos_def
! mz_ht_20040318- end of update for posdef

           ! mz_ht_20031124+
           ! the wet deposition flux should be zero at the
           ! upper boundary and not
           ! taking any values from the timestep before
  do jl=1,kproma
    IF (lscav_cv)   kpp_rain_cv(jl,1,:,jrow)    = 0.0_dp
    IF (lscav_cv)   aero_flx(1)%aer_flx_cv(jl,1,:,jrow)     = 0.0_dp
    zrainh(jl)      = 0._dp
    zcover_conv(jl) = 0._dp
  enddo
           ! mz_ht_20031124-        

  DO 120 jk=ktopm2,klev
!DIR$ IVDEP
!OCL NOVREC
     DO 115 jl=1,kproma
        IF(ldcum(jl).AND.jk.GE.kctop(jl)-1) THEN      
! mim_sb_20090917+
! op_mm_20140327+
! for 2d fields 

           IF (l_lgmc_diag) THEN
              zrcpm=1._dp/pcpen(jl,jk)
              if (jk < klev) then
!!                 conv_tte_up(jl,jk,jrow) = &
!!                      (pmfus(jl,jk)-pmfus(jl,jk+1))*g/(paphp1(jl,jk)-paphp1(jl,jk+1))*zrcpm

                 conv_tte_up_2d(jl,jk) = &
                      (pmfus(jl,jk)-pmfus(jl,jk+1))*g/(paphp1(jl,jk)-paphp1(jl,jk+1))*zrcpm

 !!                conv_tte_su(jl,jk,jrow)   = &
 !!                     (-(pmfu(jl,jk)+pmfd(jl,jk))*(pcpcu(jl,jk)*ptenh(jl,jk)+pgeoh(jl,jk))         &
 !!                     + (pmfu(jl,jk+1)+ &
 !!                     pmfd(jl,jk+1))*(pcpcu(jl,jk+1)*ptenh(jl,jk+1)+pgeoh(jl,jk+1)))*  &
 !!!                     g/(paphp1(jl,jk)-paphp1(jl,jk+1))*zrcpm


                 conv_tte_su_2d(jl,jk)   = &
                      (-(pmfu(jl,jk)+pmfd(jl,jk))*(pcpcu(jl,jk)*ptenh(jl,jk)+pgeoh(jl,jk))         &
                      + (pmfu(jl,jk+1)+ &
                      pmfd(jl,jk+1))*(pcpcu(jl,jk+1)*ptenh(jl,jk+1)+pgeoh(jl,jk+1)))*  &
                      g/(paphp1(jl,jk+1)-paphp1(jl,jk))*zrcpm

              else
 !!                conv_tte_up(jl,jk,jrow) = pmfus(jl,jk)* g/(paphp1(jl,jk)-paphp1(jl,jk+1))*zrcpm
                 conv_tte_up_2d(jl,jk) = pmfus(jl,jk)* g/(paphp1(jl,jk)-paphp1(jl,jk+1))*zrcpm

  !!               conv_tte_su(jl,jk,jrow) = &
  !!                    -(pmfu(jl,jk)+pmfd(jl,jk))*(pcpcu(jl,jk)*ptenh(jl,jk)+pgeoh(jl,jk))*       &
  !!                    g/(paphp1(jl,jk)-paphp1(jl,jk+1))*zrcpm 
                 
                 conv_tte_su_2d(jl,jk) = &
                      -(pmfu(jl,jk)+pmfd(jl,jk))*(pcpcu(jl,jk)*ptenh(jl,jk)+pgeoh(jl,jk))*       &
                      g/(paphp1(jl,jk+1)-paphp1(jl,jk))*zrcpm 
              endif
           END IF
! op_mm_20140327-
! mim_sb_20090917-

         
           pmfus(jl,jk)=pmfus(jl,jk)-pmfu(jl,jk)*                      &
                               (pcpcu(jl,jk)*ptenh(jl,jk)+pgeoh(jl,jk))
           pmfuq(jl,jk)=pmfuq(jl,jk)-pmfu(jl,jk)*pqenh(jl,jk)
           IF(lddraf(jl).AND.jk.GE.kdtop(jl)) THEN
! mim_sb_20090917+
! op_mm_20140327+
! for 2d fields 

              IF (l_lgmc_diag) THEN
                 if (jk < klev) then
 !!                   conv_tte_do(jl,jk,jrow) = &
 !!                        (pmfds(jl,jk)-pmfds(jl,jk+1))*g/(paphp1(jl,jk)-paphp1(jl,jk+1))*zrcpm

                      conv_tte_do_2d(jl,jk) = &
                           - (pmfds(jl,jk)-pmfds(jl,jk+1))*g/(paphp1(jl,jk)-paphp1(jl,jk+1))*zrcpm  
                 else
                !!    conv_tte_do(jl,jk,jrow) = &
                !!         pmfds(jl,jk)*g/(paphp1(jl,jk)-paphp1(jl,jk+1))*zrcpm
                    conv_tte_do_2d(jl,jk) = &
                         - pmfds(jl,jk)*g/(paphp1(jl,jk)-paphp1(jl,jk+1))*zrcpm
                 endif
              END IF
! op_mm_20140327-
! mim_sb_20090917-
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

      ! mz_ht_20030630+
      ! Introduction of the calculation of the in-cloud 
      !    convective wet deposition fluxes also considering some of the 
      !    aqueous-phase chemistry. 
      IF (lpos_def) THEN
        zmfu(1:kproma) = zpmfun(1:kproma,jk)
      ELSE
        zmfu(1:kproma) = pmfu(1:kproma,jk)
      ENDIF
      IF (lscav_cv.and..not.CVTRANS) &
        CALL scav_cvdep (kproma, kbdim, klev,   jk,  jrow, ztmst,    &
          ktrac,  ktopm2, ldcum,  kctop,                             &
          kcbot,  kdtop,  lddraf, pmfu,   pmfd,   pdmfup, pdmfdp,    &
          pmfuxt, pmfdxt, zrainh, zcover_conv, pten, zmfu)
      ! mz_ht_20030630-

     DO 1154 jt=1,ktrac
        DO 1152 jl=1,kproma

! mz_ht_20040318+ update for posdef
           if (lpos_def) then 
              IF(ldcum(jl).AND.jk.GE.kctop(jl)-1) THEN             
 ! op_ck_20031001+
                 IF (ldp(jl,jk)) THEN
                    pmfuxt(jl,jk,jt)=pmfuxt(jl,jk,jt)                     &
                                    -zmass(jl,jk)*pxten(jl,jk-1,jt)
                 ELSE
                    pmfdxt(jl,jk,jt)=pmfdxt(jl,jk,jt)                     &
                                     -zmass(jl,jk)*pxten(jl,jk,jt)
                 ENDIF
              ELSE
                 pmfuxt(jl,jk,jt)=0._dp
                 pmfdxt(jl,jk,jt)=0._dp
              ENDIF
! op_ck_20031001-
           ELSE

              IF(ldcum(jl).AND.jk.GE.kctop(jl)-1) THEN
                 pmfuxt(jl,jk,jt)=pmfuxt(jl,jk,jt)                        &
                                 -pmfu(jl,jk)*pxtenh(jl,jk,jt)
                 IF(lddraf(jl).AND.jk.GE.kdtop(jl)) THEN
                    pmfdxt(jl,jk,jt)=pmfdxt(jl,jk,jt)                     &
                                    -pmfd(jl,jk)*pxtenh(jl,jk,jt)
                 ELSE
                    pmfdxt(jl,jk,jt)=0._dp
                 ENDIF
              ELSE
                 pmfuxt(jl,jk,jt)=0._dp
                 pmfdxt(jl,jk,jt)=0._dp
              ENDIF
           ENDIF   !lpos_def
! mz_ht_20040318- end of update for posdef

1152    END DO
1154 END DO
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
! op_ck_20031001+/ mz_ht_20040504 update for lpos_def
     if (.not.lpos_def) then
     
        DO 1254 jt=1,ktrac
!DIR$ IVDEP
!OCL NOVREC
           DO 1252 jl=1,kproma
              IF(ldcum(jl).AND.jk.GT.kcbot(jl)) THEN
                 ikb=kcbot(jl)
                 zzp=(paphp1(jl,klevp1)-paphp1(jl,jk))/                   &
                      (paphp1(jl,klevp1)-paphp1(jl,ikb))
                 zzp=MERGE(zzp**2,zzp,ktype(jl).EQ.3)
                 pmfuxt(jl,jk,jt)=pmfuxt(jl,ikb,jt)*zzp
              ENDIF
1252       END DO
1254    END DO

     ENDIF    ! lpos_def
! op_ck_20031001-/ mz_ht_20040504 update for lpos_def
!
130 END DO

! mim_sb_20090917+
  IF (l_lgmc_diag) THEN
     DO jk=ktopm2,klev
!DIR$ IVDEP
!OCL NOVREC
!op_mm_20140327+
! for 2d fields 
       do  jl=1,kproma
        IF(ldcum(jl).AND.jk.GT.kcbot(jl)) THEN
           zrcpm=1._dp/pcpen(jl,jk)
           if (jk < klev) then
!!            conv_tte_up(jl,jk,jrow) = (pmfus(jl,jk)-pmfus(jl,jk+1))*g/(paphp1(jl,jk)-paphp1(jl,jk+1))*zrcpm
              conv_tte_up_2d(jl,jk) = pmfus(jl,jk)* g/(paphp1(jl,jk)-paphp1(jl,jk+1))*zrcpm 
           else
 !!           conv_tte_up(jl,jk,jrow) = pmfus(jl,jk)* g/(paphp1(jl,jk)-paphp1(jl,jk+1))*zrcpm
              conv_tte_up_2d(jl,jk) = pmfus(jl,jk)* g/(paphp1(jl,jk)-paphp1(jl,jk+1))*zrcpm
           endif
        ENDIF
       enddo
     ENDDO
  END IF
! op_mm_20140327-
! mim_sb_20090917-  
  
!
!
!*    2.            CALCULATE RAIN/SNOW FALL RATES
!*                  CALCULATE MELTING OF SNOW
!*                  CALCULATE EVAPORATION OF PRECIP
!                   -------------------------------
!

! mz_ht_20040318+ changes in the original/updated to avoid unneccessary 
!                 calculations, explicit type casts. Pointers that give 
!                 values to messy_main_data

! mz_lg_20030827+ added to assign the parameters needed for the NOx
!     lighting calculations

      do jl=1,kproma
        if (ktype(jl).eq.1) then
!!$          cu_top(jl,jrow)=REAL(kctop(jl),dp)
!!$          cu_bot(jl,jrow)=REAL(kcbot(jl),dp)
           cu_top_1d(jl)=REAL(kctop(jl),dp)    ! op_mm_20140521
           cu_bot_1d(jl)=REAL(kcbot(jl),dp)

!  Determine 0C cloud level for calculating fraction 
!  of InterCloud lightning
!!$          cu_freeze(jl,jrow)=REAL(kctop(jl),dp)
           cu_freeze_1d(jl)=REAL(kctop(jl),dp)
          do jk=kctop(jl),kcbot(jl)
!!$            if (pten(jl,jk).le.273.15_dp) cu_freeze(jl,jrow)=REAL(jk,dp)
             if (pten(jl,jk).le.273.15_dp) cu_freeze_1d(jl)=REAL(jk,dp)
          enddo          
        endif
        ! mz_pj_20050615+
        if (ktype(jl).eq.3) then
!!$          cu_top_mid(jl,jrow)=REAL(kctop(jl),dp)
!!$          cu_bot_mid(jl,jrow)=REAL(kcbot(jl),dp)
           cu_top_mid_1d(jl)=REAL(kctop(jl),dp)
           cu_bot_mid_1d(jl)=REAL(kcbot(jl),dp)
!  Determine 0C cloud level for calculating fraction 
!  of InterCloud lightning
!!$          cu_freeze_mid(jl,jrow)=REAL(kctop(jl),dp)
             cu_freeze_mid_1d(jl)=REAL(kctop(jl),dp)   ! op_mm_20140521
          do jk=kctop(jl),kcbot(jl)
!!$            if (pten(jl,jk).le.273.15_dp) cu_freeze_mid(jl,jrow)=REAL(jk,dp)
             if (pten(jl,jk).le.273.15_dp) cu_freeze_mid_1d(jl)=REAL(jk,dp)
          enddo          
        endif
        ! mz_pj_20050615-
      enddo

! Calculate 'updraft' =  mass flux / density 
! xupdr 
!
      do jl=1,kproma
         if (kctop(jl).ne.0) then
            do jk=kctop(jl),kcbot(jl)
               zrhoa=paphp1(jl,jk)/(ptenh(jl,jk)*rd)
      !!         cu_uvelo(jl,jk,jrow)=pmfu(jl,jk)/zrhoa    ! op_mm_20140327 for 2d fields 
               cu_uvelo_2d(jl,jk)=pmfu(jl,jk)/zrhoa
            enddo
         endif
      enddo
! mz_lg_20030827- end

! mz_ht_20040318- end of update

!! 200   CONTINUE
  DO 210 jl=1,kproma
     prfl(jl)=0._dp
     psfl(jl)=0._dp
     prain(jl)=0._dp
210 END DO
  DO 220 jk=ktopm2,klev
     DO 215 jl=1,kproma
        IF(ldcum(jl)) THEN
           prain(jl)=prain(jl)+pdmfup(jl,jk)
           IF(pten(jl,jk).GT.tmelt) THEN
              prfl(jl)=prfl(jl)+pdmfup(jl,jk)+pdmfdp(jl,jk)
              IF(psfl(jl).GT.0._dp.AND.pten(jl,jk).GT.ztmelp2) THEN
                 zfac=zcons1*(1._dp+vtmpc2*pqen(jl,jk))                   &
                             *(paphp1(jl,jk+1)-paphp1(jl,jk))
                 zsnmlt=MIN(psfl(jl),zfac*(pten(jl,jk)-ztmelp2))
                 pdpmel(jl,jk)=zsnmlt
                 psfl(jl)=psfl(jl)-zsnmlt
                 prfl(jl)=prfl(jl)+zsnmlt
              END IF
           ELSE
              psfl(jl)=psfl(jl)+pdmfup(jl,jk)+pdmfdp(jl,jk)
           END IF
        END IF
     !mz_ht_20040212+
!!        cv_precflx(jl,jk,jrow) = prfl(jl)  ! op_mm_20140327
!!       cv_snowflx(jl,jk,jrow) = psfl(jl)   ! for 2d fields 
        cv_precflx_2d(jl,jk) = prfl(jl) ! for convec channel: rain flux
        cv_snowflx_2d(jl,jk) = psfl(jl) !                     snow flux
   !mz_ht_20040212-
215  END DO
220 END DO
  DO 230 jl=1,kproma
     prfl(jl)=MAX(prfl(jl),0._dp)
     psfl(jl)=MAX(psfl(jl),0._dp)
     zpsubcl(jl)=prfl(jl)+psfl(jl)
230 END DO
  DO 240 jk=ktopm2,klev
     DO 235 jl=1,kproma
       IF(ldcum(jl).AND.jk.GE.kcbot(jl).AND.zpsubcl(jl).GT.1.e-20_dp) THEN
           zrfl=zpsubcl(jl)
           zrnew=(MAX(0._dp,SQRT(zrfl/zcucov)-                            &
                        cevapcu(jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))*      &
                        MAX(0._dp,pqsen(jl,jk)-pqen(jl,jk))))             &
                        **2*zcucov
           zrmin=zrfl-zcucov*MAX(0._dp,0.8_dp*pqsen(jl,jk)-pqen(jl,jk))      &
                        *zcons2*(paphp1(jl,jk+1)-paphp1(jl,jk))
           zrnew=MAX(zrnew,zrmin)
           zrfln=MAX(zrnew,0._dp)
           zdrfl=MIN(0._dp,zrfln-zrfl)
           pdmfup(jl,jk)=pdmfup(jl,jk)+zdrfl
           zpsubcl(jl)=zrfln

           !mz_ht_20040212+   for convect channel
           zrsum=prfl(jl)+psfl(jl)
           zdpevap=zpsubcl(jl)-zrsum
!!           cv_precflx(jl,jk,jrow) = cv_precflx(jl,jk,jrow) + &     ! op_mm_20140327
!!                zdpevap * prfl(jl) * (1._dp/MAX(1.e-20_dp,zrsum))  ! for 2d fields 
 !!          cv_snowflx(jl,jk,jrow) = cv_snowflx(jl,jk,jrow) + &
 !!               zdpevap * psfl(jl) * (1._dp/MAX(1.e-20_dp,zrsum))
           cv_precflx_2d(jl,jk) = cv_precflx_2d(jl,jk) + &
                zdpevap * prfl(jl) * (1._dp/MAX(1.e-20_dp,zrsum))
           cv_snowflx_2d(jl,jk) = cv_snowflx_2d(jl,jk) + &
                zdpevap * psfl(jl) * (1._dp/MAX(1.e-20_dp,zrsum))
           !mz_ht_20040212-

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
END SUBROUTINE tiedtke_cuflx

!==============================================================================

SUBROUTINE tiedtke_cuini(kproma, kbdim, klev, klevp1, klevm1,          &
           pten,     pqen,     pqsen,    pxen,     puen,     pven,     &
           ptven,    ktrac,                                            &
           pxten,    pxtenh,   pxtu,     pxtd,     pmfuxt,   pmfdxt,   &
           pverv,    pgeo,     paphp1,   pgeoh,                        &
           ptenh,    pqenh,    pqsenh,   pxenh,    klwmin,             &
           ptu,      pqu,      ptd,      pqd,                          &
           puu,      pvu,      pud,      pvd,                          &
           pmfu,     pmfd,     pmfus,    pmfds,                        &
           pmfuq,    pmfdq,    pdmfup,   pdmfdp,                       &
           pcpen,    pcpcu,                                            &
           pdpmel,   plu,      plude,    pqude,    klab,               &
           zpmfun)   ! op_ck_20031001/mz_ht_20040318 for posdef update

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

INTEGER, INTENT (IN) :: kproma, kbdim, klev, klevp1, klevm1, ktrac

REAL(dp) :: pten(kbdim,klev),          pqen(kbdim,klev),                  &
            puen(kbdim,klev),          pven(kbdim,klev),                  &
            pqsen(kbdim,klev),         pverv(kbdim,klev),                 &
            pgeo(kbdim,klev),          pgeoh(kbdim,klev),                 &
            paphp1(kbdim,klevp1),      ptenh(kbdim,klev),                 &
            pxenh(kbdim,klev),         pxen(kbdim,klev),                  &
            ptven(kbdim,klev),                                           &
            pqenh(kbdim,klev),         pqsenh(kbdim,klev)
REAL(dp) :: pcpen(kbdim,klev),         pcpcu(kbdim,klev)
!
REAL(dp) :: ptu(kbdim,klev),           pqu(kbdim,klev),                   &
            ptd(kbdim,klev),           pqd(kbdim,klev),                   &
            puu(kbdim,klev),           pud(kbdim,klev),                   &
            pvu(kbdim,klev),           pvd(kbdim,klev),                   &
            pmfu(kbdim,klev),          pmfd(kbdim,klev),                  &
            pmfus(kbdim,klev),         pmfds(kbdim,klev),                 &
            pmfuq(kbdim,klev),         pmfdq(kbdim,klev),                 &
            pdmfup(kbdim,klev),        pdmfdp(kbdim,klev),                &
            plu(kbdim,klev),           plude(kbdim,klev),                 &
            pqude(kbdim,klev)
REAL(dp) :: pdpmel(kbdim,klev)
INTEGER  :: klab(kbdim,klev),          klwmin(kbdim)
!
REAL(dp) :: zwmax(kbdim)
REAL(dp) :: zph(kbdim)
LOGICAL  :: loflag(kbdim)
REAL(dp) :: pxten(kbdim,klev,ktrac),   pxtenh(kbdim,klev,ktrac),          &
            pxtu(kbdim,klev,ktrac),    pxtd(kbdim,klev,ktrac),            &
            pmfuxt(kbdim,klev,ktrac),  pmfdxt(kbdim,klev,ktrac)
REAL(dp) :: zpmfun(kbdim,klev) ! op_ck_20031001/mz_ht_20040318 for posdef update

INTEGER  :: jk, jl, jt, ik, icall
REAL(dp) :: zarg, zcpm, zzs
!
!
!----------------------------------------------------------------------
!
!*    1.           SPECIFY LARGE SCALE PARAMETERS AT HALF LEVELS
!*                 ADJUST TEMPERATURE FIELDS IF STATICLY UNSTABLE
!*                 FIND LEVEL OF MAXIMUM VERTICAL VELOCITY
!                  ----------------------------------------------
!
!! 100 CONTINUE
    DO 101 jk=1,klev
       DO 102 jl=1,kproma
!          pcpen(jl,jk)=cpd*(1.+vtmpc2*MAX(pqen(jl,jk),0.0))
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
     CALL tiedtke_cuadjtq(kproma, kbdim, klev, ik,                             &
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
        zzs=MAX(pcpcu(jl,jk)*ptenh(jl,jk)+pgeoh(jl,jk),                         &
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
!! 200 CONTINUE
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
        zpmfun(jl,jk)=0._dp ! op_ck_20031001/mz_ht_20040318  for posdef update
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
END SUBROUTINE tiedtke_cuini


!==============================================================================
!==============================================================================
END MODULE MESSY_CONVECT_TIEDTKE
