MODULE MESSY_CONVECT


! This module contains the organizing routines for each convection parametrizations.
! They are called from MESSY_CONVECTION_E5.

! The subroutines are in the following files:
! TIEDTKE :   MESSY_convection_tiedtke


! By now there are included:

! 1) TIEDTKE                    cumastr
! 2) TIEDTKE                    cumastrh
! 3) TIEDTKE                    cumastrt
! 4) ECMWF                      cumastr 
! 5) ZHANG / HACK / McFRALANE   cumastr
! 6) BECHTOLD                   cumastr
! 7) EMANUEL                    cumastr
! 8) DONNER (2006/2007)         cumastr

! Author : H.Tost, MPICH, March 2004

! original code from Tiedtke (see below)
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
 
 
!!$#ifdef __ibm__
!!$  USE mo_specfun, ONLY: merge
!!$#endif

  IMPLICIT NONE
  
  SAVE

  CHARACTER(len=*), PARAMETER :: MODSTR='convect'
  CHARACTER(LEN=*), PARAMETER :: modver='2.0'
  LOGICAL :: lconvection       ! global switch
  INTEGER :: convect_param     ! global switch

  LOGICAL :: altconv           ! switch for the two Zhang/Hack parametrisations
  LOGICAL :: evap_sub          ! switch for the subcloud evaporation in Zhang

  LOGICAL :: ltransport        ! switch for tracer transport in Emanuel scheme
 
  LOGICAL :: L_LGMC = .FALSE.  ! submodel LGMC is ON ! op_sb_20171018

CONTAINS

!=============================================================================
  SUBROUTINE init_convection(iou, fstat)

    ! SCAV MODULE ROUTINE (CORE)
    !
    ! READ CONVECTION NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! Author: H.Tost, MPICH, March 2004

    USE MESSY_CONVECT_TIEDTKE
    USE MESSY_CONVECT_MEM,            ONLY: ODEEP, OSHAL, ODOWN, OREFRESH_ALL, &
                                            OSETTADJ, OUVTRANS, OCHTRANS,      &
                                            KENSM, KICE, PTADJD, PTADJS
    USE MESSY_CONVECT_ECMWF_PARAM,    ONLY: PEN => LMFPEN, SHAL => LMFSCV,     &
                                            MID => LMFMID, DOWN => LMFDD,      &
                                            FRIC => LMFDUDV,                   &
                                            LMFTRAC, LEPCLD, LMFSCL_WSTAR

    IMPLICIT NONE

    INTRINSIC :: TRIM
    ! I/O
    INTEGER, INTENT(IN)  :: iou   ! logical I/O unit

    ! (LOCAL) NAMELIST VARIABLES FOR convection MODULE CONTROL

    INTEGER :: convparam, ensemble, ice
    REAL    :: adjtimed, adjtimes
    LOGICAL :: posdef, penetrative, shallow, midlevel, downdrafts, friction, &
               alterconv,                                                    &
               deep, tracertrans, refresh, adjustment

    NAMELIST /CTRL/            convparam   
    NAMELIST /CTRL_TIEDTKE/    penetrative, shallow, midlevel, &
                               downdrafts, friction, &
                               posdef, &
                               ! fb_mk_20120116+
                               rset_cmfctop , rset_cprcon, rset_entrscv, &
                               ! fb_mk_20120116-
                               ! fb_mk_20140210+
                               rset_entrmid, rset_entrpen
                               ! fb_mk_20140210-
                               
    NAMELIST /CTRL_ECMWF/      penetrative, shallow, midlevel, &
                               downdrafts, friction,           &
                               tracertrans, LEPCLD, LMFSCL_WSTAR
    
    NAMELIST /CTRL_ZHANG/      evap_sub, alterconv
    NAMELIST /CTRL_BECHTOLD/   deep, shallow, downdrafts, ice, friction, &
                               tracertrans, ensemble, refresh,           &
                               adjustment, adjtimed, adjtimes
    NAMELIST /CTRL_EMANUEL/    ltransport

    ! LOCAL
    LOGICAL                           :: lex          ! file exists ?
    INTEGER                           :: fstat        ! file status

    ! INITIALIZE GLOBAL CONTROL VARIABLES
    fstat = 0
    lconvection =.false.
    convect_param = 0
    lpos_def   = .false.

    lmfpen       = .false.
    lmfscv       = .false. 
    lmfmid       = .false.  
    lmfdd        = .false.
    lmfdudv      = .false.

    LMFTRAC      = .false.
    LEPCLD       = .false.
    LMFSCL_WSTAR = .false.
    
    altconv      = .false.
    EVAP_SUB     = .false.

    ODEEP        = .false.
    OSHAL        = .false.
    ODOWN        = .false.
    OREFRESH_ALL = .false.
    OSETTADJ     = .false.
    OUVTRANS     = .false.
    OCHTRANS     = .false.
    KENSM        = 0
    KICE         = 0 
    PTADJD       = 0.0_dp
    PTADJS       = 0.0_dp

    ltransport   = .false.

    ! INITIALIZE NAMELIST VARIABLES

    convparam = 0
    ensemble  = 0
    ice       = 0

    posdef      = .false.
! fb_mk_20120116+
    rset_cmfctop%l = .FALSE.
    rset_cmfctop%v = 0.0_dp
    rset_cprcon%l  = .FALSE.
    rset_cprcon%v  = 0.0_dp
    rset_entrscv%l = .FALSE.
    rset_entrscv%v = 0.0_dp
! fb_mk_20120116-
! fb_mk_20140210+
    rset_entrpen%l = .FALSE.
    rset_entrpen%v = 0.0_dp
    rset_entrmid%l = .FALSE.
    rset_entrmid%v = 0.0_dp
! fb_mk_20140210-

    penetrative = .false.
    shallow     = .false.
    midlevel    = .false.   
    downdrafts  = .false.   
    friction    = .false.

    alterconv   = .false.

    deep        = .false.
    tracertrans = .false.
    refresh     = .false.
    adjustment  = .false.
    
    adjtimed    = 0.0
    adjtimes    = 0.0
 
    ! INPUT NAMELIST
    ! CHECK IF FILE EXISTS, YES: SWITCH convect ON
    !                       NO : KEEP convect SWITCHED OFF
    INQUIRE(file=TRIM(modstr)//'.nml', exist=lex)
    IF (.NOT.lex) THEN
       WRITE(*,*) 'WARNING *** FILE '//TRIM(modstr)//'.nml'//'  NOT FOUND !'
       WRITE(*,*) ' CONVECTION SWITCHED OFF !'
       WRITE(*,*) '******************************************************'
       fstat=1
       RETURN
    END IF
    ! SET GLOBAL SWITCH
    lconvection = .true.

    ! READ NAMELIST
    OPEN(iou,file=TRIM(modstr)//'.nml')
    WRITE(*,*) 'Reading namelist from '//TRIM(modstr)//'.nml', &
         ' (unit ',iou,') ...'
    READ(iou, NML=CTRL, IOSTAT=fstat)
    IF (fstat /= 0) THEN
      WRITE(*,*) 'ERROR *** READ ERROR in NAMELIST ', &
        TRIM(modstr)//'.nml'//' !'
      WRITE(*,*) '******************************************************'
      RETURN
    END IF
    convect_param = convparam
    print*, "convect_param: ", convect_param
    IF (convect_param.gt.0.and.convect_param.le.3) then
      READ(iou, NML=CTRL_TIEDTKE, IOSTAT=fstat)
      IF (fstat /= 0) THEN
        WRITE(*,*) 'ERROR *** READ ERROR in NAMELIST ', &
          TRIM(modstr)//'.nml'//' !'
        WRITE(*,*) '******************************************************'
        RETURN
      END IF
    ENDIF

    IF (convect_param.eq.4) then
      READ(iou, NML=CTRL_ECMWF, IOSTAT=fstat)
      IF (fstat /= 0) THEN
        WRITE(*,*) 'ERROR *** READ ERROR in NAMELIST ', &
          TRIM(modstr)//'.nml'//' !'
        WRITE(*,*) '******************************************************'
        RETURN
      END IF
    ENDIF

    IF (convect_param.eq.5) then
      READ(iou, NML=CTRL_ZHANG, IOSTAT=fstat)
      IF (fstat /= 0) THEN
        WRITE(*,*) 'ERROR *** READ ERROR in NAMELIST ', &
          TRIM(modstr)//'.nml'//' !'
        WRITE(*,*) '******************************************************'
        RETURN
      END IF
    ENDIF

    IF (convect_param.eq.6) then
      READ(iou, NML=CTRL_BECHTOLD, IOSTAT=fstat)
      IF (fstat /= 0) THEN
        WRITE(*,*) 'ERROR *** READ ERROR in NAMELIST ', &
          TRIM(modstr)//'.nml'//' !'
        WRITE(*,*) '******************************************************'
        RETURN
      END IF
    ENDIF

    IF (convect_param.eq.7) then
      READ(iou, NML=CTRL_EMANUEL, IOSTAT=fstat)
      IF (fstat /= 0) THEN
        WRITE(*,*) 'ERROR *** READ ERROR in NAMELIST ', &
          TRIM(modstr)//'.nml'//' !'
        WRITE(*,*) '******************************************************'
        RETURN
      END IF
    ENDIF
    
    if (convect_param  < 1 .or. convect_param > 8) lconvection = .false.

    lmfpen   = penetrative
    lmfscv   = shallow 
    lmfmid   = midlevel  
    lmfdd    = downdrafts
    lmfdudv  = friction 
    lpos_def = posdef

  
    PEN     = penetrative
    SHAL    = shallow 
    MID     = midlevel  
    DOWN    = downdrafts
    FRIC    = friction 
    LMFTRAC = tracertrans
    
    altconv = alterconv

    ODEEP        = deep
    OSHAL        = shallow
    ODOWN        = downdrafts
    OREFRESH_ALL = refresh
    OSETTADJ     = adjustment
    OUVTRANS     = friction
    OCHTRANS     = tracertrans
    KENSM        = ensemble
    KICE         = ice
    PTADJD       = REAL(adjtimed,dp)
    PTADJS       = REAL(adjtimes,dp)

    CLOSE(iou)

    convect_param = convparam
    if (.not.lconvection) then
       WRITE(*,*) '******************************************************'
       WRITE(*,*) '       WARNING: NO CONVECTION is SELECTED !!!!!'
       WRITE(*,*) '******************************************************'
    ENDIF

  END SUBROUTINE init_convection
! =====================================================================================

! -------------------------------------------------------------------------
!mz_ht_20040205+
! subroutine that calculates the estimate of a convective cloud cover

SUBROUTINE calc_conv_cover(cover, pmassfu, zrhoa, ztmst)

  USE MESSY_MAIN_CONSTANTS_MEM,    ONLY: dp

  IMPLICIT NONE
  
  REAL(dp), INTENT(OUT) :: cover
  REAL(dp), INTENT(IN)  :: pmassfu, zrhoa, ztmst

!local
  REAL(dp) :: vup_min, covadj
!      estimated convective cloud cover adapted from MATCH from Mark Lawrence

!     first determine the fraction swept out by precip from the core cloud;
!        the precipitating part is the cloud width needed to support the 
!        updraft mass flux with a nominal updraft velocity (of 1 m/s);
!        this might miss out on some diffs betw land and sea, but should be 
!        good for a first attempt; 

!     ATTENTION: Do not sum up pclcover and conv_cover for a total cover!
!                They might get higher values than 1 !!!!!!

        cover = 0.
          
        vup_min = 1._dp   ! [m/s] necessary updraft velocity to support convective rainfall
                          ! still to be tested or modified

! mz_ht_20111130+      
!!$        cover = pmassfu  / (vup_min * zrhoa) 
        cover = MIN(1._dp, MAX(pmassfu  / (vup_min * zrhoa) , 0.0001_dp) )
! mz_ht_20111130-

!     next compute a "coverage adjustment"; this is to help reduce the
!     time-step dependence of the process, since it is expected that for longer
!     time-steps the convective precip will tend to strip out the entire
!     precipitating column.  (Could collect statistics in a future
!     study to see if the threshold chosen here is really appropriate...)
!     Assume:
!     1) the characteristic horiz dimension of any individual tower is 10 km;
!     2) the cells propogate horizontally with a characteristic speed of 10 m/s

!      covadj increases the conv_cover a lot and it is getting values higher than 1, 
!      has to be checked if it is appropiate here......

!      covadj = max( 1., 10.*ztmst/1.e4 ) 
      
!      cover =  cover * covadj

END SUBROUTINE calc_conv_cover
!mz_ht_20040205-
! -------------------------------------------------------------------------

! mz_jd_20161011+
! subroutine that calculates the area fraction of cloud cores according to Sikma 
! and Ouwersloot (2015). Next to this, standard deviation and skewness of the
! beta distribution of specific humidity are calculated, according to Tompkins
! (2002).
SUBROUTINE calc_conv_cover_sikma(cover_sikma, qt,    qs,   &
                                 qstdev,      Q2,    qskew,&
                                 xvar,        xskew, aclc  )

USE MESSY_MAIN_CONSTANTS_MEM,     ONLY: dp
USE MESSY_CLOUD_ORI,              ONLY: cbeta_pq ! beta distribution parameter p

IMPLICIT NONE

REAL(dp), INTENT(OUT)    :: cover_sikma,                 & ! area fraction of cloud cores
                            qskew,                       & ! skewness of total water amount
                            qstdev,                      & ! standard deviation of total water amount
                            Q2                             ! proxy for area fraction of cloud cores
REAL(dp), INTENT(IN)     :: qt,                          & ! total specific humidity
                            qs,                          & ! saturation specific humidity
                            xvar,                        & ! proxy for standard deviation of total water amount
                                                           ! (in fact: beta distribution width b-a)
                            xskew,                       & ! beta distribution parameter q
                            aclc                           ! area fraction of clouds

qstdev      = (xvar/(cbeta_pq + xskew)) * SQRT((cbeta_pq * xskew) / &
              (cbeta_pq + xskew + 1.0_dp))
qskew       = ((2.0_dp * (xskew - cbeta_pq)) / (cbeta_pq + xskew + 2.0_dp)) * &
              SQRT((cbeta_pq + xskew + 1.0_dp) / (cbeta_pq * xskew))

Q2          = MIN((qt - qs) / MAX(qstdev, 1.e-10_dp), -1.e-10_dp)

cover_sikma = 0.292_dp * Q2**(-2)

! Make sure that the area fraction of cloud cores is not larger than the area
! fraction of clouds.
cover_sikma = MAX(MIN(1._dp,cover_sikma),0.01_dp)

END SUBROUTINE
! mz_jd_20161011-

!======================================================================================

!--------------- Organizing routines for convection schemes ---------------------------

!======================================================================================
  

  SUBROUTINE tiedtke_cumastr(kproma, kbdim, klev, klevp1, klevm1,         &
             jrow,     ztmst,    nn,       ktrac,    ilab,                &
             pten,     pqen,     pxen,     puen,     pven,                &
             ptven,    ldland,                                            &
             pxten,    pxtu,     pxtte,                                   &
             pverv,    pqsen,    pqhfla,                                  &
             paphp1,   pgeo,                                              & 
             ptte,     pqte,     pvom,     pvol,                          &
             prsfc,    pssfc,    paprc,    paprs,    pxtec,               &
             pqtec,                                                       &
             ldcum,    ktype,    kcbot,    kctop,                         &
             pmfu,     pmfd,     prain,    ptracconv, delta_time,         &
             pvgustcon, pdtke_con) ! op_mm_20140327 added pvgustcon, pdtke_con
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
!          *CUMASTR* IS CALLED FROM *CONVECTION_E5*
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

    USE MESSY_CONVECT_TIEDTKE

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: kbdim, klev, ktrac, kproma, klevp1, klevm1, jrow, nn
    
    INTEGER  :: jl, jk, jt, ptracconv(ktrac)

    REAL(dp), INTENT(IN) :: ztmst, delta_time

    REAL(dp) :: pten(kbdim,klev),        pqen(kbdim,klev),                     &
                pxen(kbdim,klev),        ptven(kbdim,klev),                    &
                puen(kbdim,klev),        pven(kbdim,klev),                     &
                ptte(kbdim,klev),        pqte(kbdim,klev),                     &
                pvom(kbdim,klev),        pvol(kbdim,klev),                     &
                pqsen(kbdim,klev),       pgeo(kbdim,klev),                     &
                paphp1(kbdim,klevp1),                                          &
                pverv(kbdim,klev),       pqude(kbdim,klev)                                            
    REAL(dp) :: ptu(kbdim,klev),         pqu(kbdim,klev),                      &
                plu(kbdim,klev),         plude(kbdim,klev),                    &
                pmfu(kbdim,klev),        pmfd(kbdim,klev),                     &
                paprc(kbdim),            paprs(kbdim),                         &
                prsfc(kbdim),            pssfc(kbdim),                         &
                prain(kbdim),            pqhfla(kbdim)
    ! op_mm_20140226 calculation of convective gust:
    REAL(dp) :: pvgustcon(kbdim)      
    INTEGER  :: kcbot(kbdim),            kctop(kbdim),                         &
                ktype(kbdim)
    REAL(dp) :: pxtec(kbdim,klev),       pqtec(kbdim,klev)
    REAL(dp) :: ztenh(kbdim,klev),       zqenh(kbdim,klev),                    &
                zxenh(kbdim,klev),                                             &
                zgeoh(kbdim,klev),       zqsenh(kbdim,klev),                   &
                ztd(kbdim,klev),         zqd(kbdim,klev),                      &
                zmfus(kbdim,klev),       zmfds(kbdim,klev),                    &
                zmfuq(kbdim,klev),       zmfdq(kbdim,klev),                    &
                zdmfup(kbdim,klev),      zdmfdp(kbdim,klev),                   &
                zmful(kbdim,klev),       zrfl(kbdim),                          &
                zuu(kbdim,klev),         zvu(kbdim,klev),                      &
                zud(kbdim,klev),         zvd(kbdim,klev)
    REAL(dp) :: zcpen(kbdim,klev),       zcpcu(kbdim,klev)
    REAL(dp) :: zentr(kbdim),            zhcbase(kbdim),                       &
                zmfub(kbdim),            zmfub1(kbdim),                        &
                zdqpbl(kbdim),           zdqcv(kbdim),                         &
                zvddraf(kbdim)          ! op_mm_20140226
    REAL(dp) :: zsfl(kbdim),             zdpmel(kbdim,klev)
    REAL(dp) :: zcape(kbdim),            zheat(kbdim)
    REAL(dp) :: zhmin(kbdim)
    REAL(dp) :: zhhatt(kbdim,klev)
    INTEGER  :: ihmin(kbdim)
    INTEGER  :: ilab(kbdim,klev),        idtop(kbdim),                         &
                ictop0(kbdim),           ilwmin(kbdim)
    REAL(dp) :: pxten(kbdim,klev,ktrac), pxtte(kbdim,klev,ktrac),              &
                pxtu(kbdim,klev,ktrac),                                        &
                zxtenh(kbdim,klev,ktrac),zxtd(kbdim,klev,ktrac),               &
                zmfuxt(kbdim,klev,ktrac),zmfdxt(kbdim,klev,ktrac)
    ! op_mm_20140227 calculation of tke due to convection:
    REAL(dp) :: pdtke_con(kbdim,klev)   
    REAL(dp) :: zpmfun(kbdim,klev) ! op_ck_20031001
    LOGICAL  :: loddraf(kbdim),          ldland(kbdim)
    LOGICAL  :: ldcum(kbdim)

    LOGICAL  :: llo1, lo

! local variables

    REAL(dp) :: zcons2, ztau,  zqumqe, zdqmin, zmfmax, zalvs,  zalvdcp, zqalv, &
                zhsat,  zqsat, zes,    zcor,   zqst1,  zdqsdt, zgam,    zzz,   &
                zhhat,  zb,    zbi,    zro,    zdz,    zdhdz,  zdepth,  zfac,  &
                zrh,    zeps,  zpbmpt
    INTEGER  :: ikb, it, it1, icum, itopm2
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
    zcons2=1._dp/(g*ztmst)
    ztau=MIN(3._dp*3600._dp,7200._dp*63._dp/REAL(nn,dp)) ! op_ff_20161020 added _dp

!
!----------------------------------------------------------------------
!
!*    2.           INITIALIZE VALUES AT VERTICAL GRID POINTS IN 'CUINI'
!                  ---------------------------------------------------
!
200 CONTINUE
    call tiedtke_cuini(kproma, kbdim, klev, klevp1, klevm1,               &
               pten,     pqen,     pqsen,    pxen,     puen,     pven,    &
               ptven,    ktrac,                                           &
               pxten,    zxtenh,   pxtu,     zxtd,     zmfuxt,   zmfdxt,  &
               pverv,    pgeo,     paphp1,   zgeoh,                       &
               ztenh,    zqenh,    zqsenh,   zxenh,    ilwmin,            &
               ptu,      pqu,      ztd,      zqd,                         &
               zuu,      zvu,      zud,      zvd,                         &
               pmfu,     pmfd,     zmfus,    zmfds,                       &
               zmfuq,    zmfdq,    zdmfup,   zdmfdp,                      &
               zcpen,    zcpcu,                                           &
               zdpmel,   plu,      plude,    pqude,    ilab,              &
               zpmfun)   ! op_ck_20031001/mz_ht_20040318 for posdef update
    IF (lookupoverflow) RETURN
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
    CALL tiedtke_cubase(kproma, kbdim, klev, klevp1, klevm1,              &
                ztenh,    zqenh,    zgeoh,    paphp1,                     &
                ptu,      pqu,      plu,                                  &
                puen,     pven,     zuu,      zvu,                        &
                zcpcu,                                                    &
                ldcum,    kcbot,    ilab)
!
    IF (lookupoverflow) RETURN
!*             (B) DETERMINE TOTAL MOISTURE CONVERGENCE AND
!*                 THEN DECIDE ON TYPE OF CUMULUS CONVECTION
!                  -----------------------------------------
!
    jk=1
    DO 310 jl=1,kproma
       zdqpbl(jl)=0.0_dp ! op_ff_20161020 add _dp
       zdqcv(jl)=pqte(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))
       idtop(jl)=0
310 END DO
    DO 320 jk=2,klev
       DO 315 jl=1,kproma
          zdqcv(jl)=zdqcv(jl)+pqte(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))
          IF(jk.GE.kcbot(jl)) zdqpbl(jl)=zdqpbl(jl)+pqte(jl,jk)           &
                                        *(paphp1(jl,jk+1)-paphp1(jl,jk))
315    END DO
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
       zdqmin=MAX(0.01_dp*zqenh(jl,ikb),1.e-10_dp) ! op_ff_20161020
       llo1=zdqpbl(jl).GT.0._dp.AND.zqumqe.GT.zdqmin.AND.ldcum(jl) ! op_ff_20161020
       zmfub(jl)=MERGE(zdqpbl(jl)/(g*MAX(zqumqe,zdqmin)),0.01_dp,llo1)
       zmfmax=(paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
       zmfub(jl)=MIN(zmfub(jl),zmfmax)
       IF(.NOT.llo1) ldcum(jl)=.FALSE.
       ktype(jl)=MERGE(1,2,zdqcv(jl).GT.MAX(0._dp,-1.1_dp*pqhfla(jl)*g)) ! op_ff_20161020
       zentr(jl)=MERGE(entrpen,entrscv,ktype(jl).EQ.1)
340 END DO
!
!!$  base_f1(1:kproma,jrow) = zmfub(1:kproma)  
    base_f1_1d(1:kproma) = zmfub(1:kproma) ! op_mm_20140521  

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
       zhcbase(jl)=zcpcu(jl,ikb)*ptu(jl,ikb)+zgeoh(jl,ikb)+zalvs*pqu(jl,ikb)
       ictop0(jl)=kcbot(jl)-1
410 END DO
    DO 430 jk=klevm1,3,-1
       !DIR$ IVDEP
       DO 420 jl=1,kproma
          zalvs=MERGE(alv,als,ztenh(jl,jk)>tmelt)
          zalvdcp=zalvs/zcpcu(jl,jk)
          zqalv=1._dp/zalvs ! op_ff_20161020
          zhsat=zcpcu(jl,jk)*ztenh(jl,jk)+zgeoh(jl,jk)+zalvs*zqsenh(jl,jk)
          it = NINT(ztenh(jl,jk)*1000._dp) ! op_ff_20161020
! ka_sv_20170406+
! ! limiting to below thermosphere
!!$       IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
          IF ((it<jptlucu1 .OR. it>jptlucu2) .AND. &
               (paphp1(jl,jk+1) >= 1.0_dp)) lookupoverflow = .TRUE.
! ka_sv_20170406-
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
          zdqsdt=(zqst1-zqsat)*1000._dp ! op_ff_20161020
          zgam=MERGE(zalvdcp*zdqsdt,zqsat*zcor*tlucub(it),LO)
          zzz=zcpcu(jl,jk)*ztenh(jl,jk)*vtmpc1
          zhhat=zhsat-(zzz+zgam*zzz)/(1._dp+zgam*zzz*zqalv)*                 &
               MAX(zqsenh(jl,jk)-zqenh(jl,jk),0._dp)
          zhhatt(jl,jk)=zhhat
          IF(jk.LT.ictop0(jl).AND.zhcbase(jl).GT.zhhat) ictop0(jl)=jk
420    END DO
430 END DO
!!
!!     DEEP CONVECTION IF CLOUD DEPTH > 200 HPA, ELSE SHALLOW
!!     (CLOUD DEPTH FROM NON-ENTRAINIG PLUME)
!  DO jl=1,kproma
!     ktype(jl)=MERGE(1,2,                                               &
!                paphp1(jl,kcbot(jl))-paphp1(jl,ictop0(jl)).gt.2.E4)
!     zentr(jl)=MERGE(entrpen,entrscv,ktype(jl).eq.1)
!  ENDDO
!!
     IF (lookupoverflow) RETURN
!
!                  FIND LOWEST POSSIBLE ORG. DETRAINMENT LEVEL
!                  -------------------------------------------
!
  DO jl=1,kproma
     zhmin(jl)=0._dp ! op_ff_20161020
     ihmin(jl)=0
     llo1=ldcum(jl).AND.ktype(jl).EQ.1
     IF(llo1) THEN
        ikb=kcbot(jl)
        ihmin(jl)=ikb
     ENDIF
  ENDDO
!
  zb=25._dp
  zbi=1._dp/(zb*g) ! op_ff_20161020
  DO jk=klev,1,-1
     DO jl=1,kproma
        llo1=ldcum(jl).AND.ktype(jl).EQ.1.AND.ihmin(jl).EQ.kcbot(jl)
        IF(llo1.AND.jk.LT.kcbot(jl).AND.jk.GE.ictop0(jl)) THEN
           zalvs=MERGE(alv,als,ztenh(jl,jk)>tmelt)
           ikb=kcbot(jl)
           zro=paphp1(jl,jk)/(rd*ztenh(jl,jk)*(1._dp+vtmpc1*zqenh(jl,jk)))
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
  CALL tiedtke_cuasc(kproma, kbdim, klev, klevp1, klevm1, jrow, ztmst,  & !mz_ht_20040317+ added ztmst, jrow
             ztenh,    zqenh,    puen,     pven,                        &
             ktrac,                                                     &
             zxtenh,   pxten,    pxtu,     zmfuxt,                      &
             pten,     pqen,     pqsen,                                 &
             pgeo,     zgeoh,    paphp1,                                &
             pqte,     pverv,    ilwmin,                                &
             ldcum,    ldland,   ktype,    ilab,                        &
             ptu,      pqu,      plu,      zuu,      zvu,               &
             pmfu,     zmfub,    zentr,                                 &
             zmfus,    zmfuq,                                           &
             zmful,    plude,    pqude,    zdmfup,                      &
             ihmin,    zhhatt,   zhcbase,  zqsenh,                      &
             zcpen,    zcpcu,                                           &
             kcbot,    kctop,    ictop0,   icum,                        &
             zpmfun)   ! op_ck_20031001/mz_ht_20040318 for posdef update
  IF (lookupoverflow) RETURN
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
!!$  base_f2(1:kproma,jrow) = zmfub(1:kproma)  
  base_f2_1d(1:kproma) = zmfub(1:kproma)   ! op_mm_20140521

  IF(lmfdd) THEN
!
!*             (A) DETERMINE LFS IN 'CUDLFS'
!                  -------------------------
!
     call tiedtke_cudlfs(kproma,   kbdim,    klev,     klevp1,          &
                 ztenh,    zqenh,    puen,     pven,                    &
                 ktrac,                                                 &
                 zxtenh,   pxtu,     zxtd,     zmfdxt,                  &
                 zgeoh,    paphp1,                                      &
                 ptu,      pqu,      zuu,      zvu,                     &
                 ldcum,    kcbot,    kctop,    zmfub,    zrfl,          &
                 ztd,      zqd,      zud,      zvd,                     &
                 pmfd,     zmfds,    zmfdq,    zdmfdp,                  &
                 zcpcu,                                                 &
                 idtop,    loddraf)
!
!*            (B)  DETERMINE DOWNDRAFT T,Q AND FLUXES IN 'CUDDRAF'
!                  -----------------------------------------------
!
     call tiedtke_cuddraf(kproma,   kbdim,    klev,     klevp1, jrow,   &
                  ztenh,    zqenh,    puen,     pven,                   &
                  ktrac,                                                &
                  zxtenh,   zxtd,     zmfdxt,                           &
                  zgeoh,    paphp1,   zrfl,                             &
                  ztd,      zqd,      zud,      zvd,                    &
                  pmfd,     zmfds,    zmfdq,    zdmfdp,                 &
                  zcpcu,                                                &
                  loddraf,                                              &
                  pxten,  & ! op_ck_20031001/mz_ht_20040318   for pos_def update
                  zvddraf ) ! op_mm_20140226 convective gust
                       
!
  END IF
  ! op_sb_20171018+
  IF (L_LGMC) THEN
     ptd_gp(1:kproma,:,jrow) = ztd(1:kproma,:)
  ENDIF
  ! op_sb_20171018-

!
!*    5.5          RECALCULATE CLOUD BASE MASSFLUX FROM A
!*                 CAPE CLOSURE FOR DEEP CONVECTION (KTYPE=1)
!*                 AND BY PBL EQUILIBRUM TAKING DOWNDRAFTS INTO
!*                 ACCOUNT FOR SHALLOW CONVECTION (KTYPE=2)
!                  -------------------------------------------
!
  DO jl=1,kproma
     zheat(jl)=0._dp ! op_ff_20161020
     zcape(jl)=0._dp ! op_ff_20161020
     zmfub1(jl)=zmfub(jl)
 !!$    base_f3(jl,jrow) = zmfub(jl)
     base_f3_1d(jl) = zmfub(jl) ! op_mm_20140521
  ENDDO
!
  DO jk=1,klev
     DO jl=1,kproma
        llo1=ldcum(jl).AND.ktype(jl).EQ.1
        IF(llo1.AND.jk.LE.kcbot(jl).AND.jk.GT.kctop(jl)) THEN
           ikb=kcbot(jl)
           zro=paphp1(jl,jk)/(rd*ztenh(jl,jk)*(1._dp+vtmpc1*zqenh(jl,jk)))
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
  CAPE(1:kproma,jrow) = zcape(1:kproma)
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
        llo1=pmfd(jl,ikb).LT.0._dp.AND.loddraf(jl) ! op_ff_20161020
        zeps=MERGE(cmfdeps,0._dp,llo1)
        zqumqe=pqu(jl,ikb)+plu(jl,ikb)-                                 &
                    zeps*zqd(jl,ikb)-(1._dp-zeps)*zqenh(jl,ikb)
        zdqmin=MAX(0.01_dp*zqenh(jl,ikb),1.e-10_dp) ! op_ff_20161020
        zmfmax=(paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
        llo1=zdqpbl(jl).GT.0._dp.AND.zqumqe.GT.zdqmin.AND.ldcum(jl)        & ! op_ff_20161020
                    .AND.zmfub(jl).LT.zmfmax
        zmfub1(jl)=MERGE(zdqpbl(jl)/(g*MAX(zqumqe,zdqmin)),             &
                               zmfub(jl),llo1)
        zmfub1(jl)=MERGE(zmfub1(jl),zmfub(jl),                          &
                             ABS(zmfub1(jl)-zmfub(jl)).LT.0.2_dp*zmfub(jl)) ! op_ff_20161020
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
     DO 5304 jt=1,ktrac
        DO 5302 jl=1,kproma
           IF(ldcum(jl)) THEN
              zfac=zmfub1(jl)/MAX(zmfub(jl),1.e-10_dp)
              zmfdxt(jl,jk,jt)=zmfdxt(jl,jk,jt)*zfac
           ENDIF
5302    END DO
5304 END DO
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
  CALL tiedtke_cuasc(kproma, kbdim, klev, klevp1, klevm1, jrow, ztmst,  & !mz_ht_20040317+ ztmst added, jrow
             ztenh,    zqenh,    puen,     pven,                        &
             ktrac,                                                     &
             zxtenh,   pxten,    pxtu,     zmfuxt,                      &
             pten,     pqen,     pqsen,                                 &
             pgeo,     zgeoh,    paphp1,                                &
             pqte,     pverv,    ilwmin,                                &
             ldcum,    ldland,   ktype,    ilab,                        &
             ptu,      pqu,      plu,      zuu,      zvu,               &
             pmfu,     zmfub,    zentr,                                 &
             zmfus,    zmfuq,                                           &
             zmful,    plude,    pqude,    zdmfup,                      &
             ihmin,    zhhatt,   zhcbase,  zqsenh,                      &
             zcpen,    zcpcu,                                           &
             kcbot,    kctop,    ictop0,   icum,                        &
             zpmfun)   ! op_ck_20031001/mz_ht_20040318 for posdef update
  IF (lookupoverflow) RETURN
  ! op_sb_20171018+
  IF (L_LGMC) THEN
     ptu_gp(1:kproma,:,jrow) = ptu(1:kproma,:)
  END IF
  ! op_sb_20171018-

!
!-----------------------------------------------------------------------
!*    7.0          DETERMINE FINAL CONVECTIVE FLUXES IN 'CUFLX'
!                  ------------------------------------------
!
700 CONTINUE
  call tiedtke_cuflx(kproma,   kbdim,    klev,     klevp1, jrow,        &
             pqen,     pqsen,    ztenh,    zqenh,                       &
             ktrac,                                                     &
             zxtenh,   zmfuxt,   zmfdxt,                                &
             paphp1,   zgeoh,                                           &
             kcbot,    kctop,    idtop,                                 &
             ktype,    loddraf,  ldcum,                                 &
             pmfu,     pmfd,     zmfus,    zmfds,                       &
             zmfuq,    zmfdq,    zmful,                                 &
             zdmfup,   zdmfdp,   zrfl,     prain,                       &
             zcpcu,                                                     &
             zcpen,                                                   & ! mim_sb_20090917
             pten,     zsfl,     zdpmel,   itopm2, ztmst,               &
             zpmfun,   pxten)             ! mz_ht_20040318+ needed for posdef update
!
!-----------------------------------------------------------------------
!
!*    8.0          UPDATE TENDENCIES FOR T AND Q IN SUBROUTINE CUDTDQ
!                  --------------------------------------------------
!
800 CONTINUE
  call tiedtke_cudtdq(kproma, kbdim, klev, klevp1, itopm2, ldcum, ktrac,&
              jrow,     paphp1,   pten,     ptte,     pqte,             &
              pxtte,    pxtec,    zmfuxt,   zmfdxt,                     &
              zmfus,    zmfds,    zmfuq,    zmfdq,                      &
              zmful,    zdmfup,   zdmfdp,   plude,                      &
              zdpmel,   zrfl,     zsfl,                                 &
              zcpen,    pqtec,    pqude,                                &
              prsfc,    pssfc,    paprc,    paprs,                      &
              delta_time, ptracconv, ztmst, pxten, pdtke_con,           &
              zqenh, ztenh) ! op_mm_20140327 added ptke_con, zqenh and ztenh
!
!-----------------------------------------------------------------------
!
!*    9.0          UPDATE TENDENCIES FOR U AND U IN SUBROUTINE CUDUDV
!                  --------------------------------------------------
!
900 CONTINUE
  IF(lmfdudv) THEN
     call tiedtke_cududv(kproma,   kbdim,    klev,     klevp1,          &
                 itopm2,   ktype,    kcbot,    paphp1,   ldcum,         &
                 puen,     pven,     pvom,     pvol,                    &
                 zuu,      zud,      zvu,      zvd,                     &
                 pmfu,     pmfd)
!
  END IF
!
1000 CONTINUE
!
!!$  base_f4(1:kproma,jrow) = zmfub(1:kproma) 
  base_f4_1d(1:kproma) = zmfub(1:kproma)   ! op_mm_20140521
  ! mz_ak_20051221+
  ! cv_cldwater(1:kproma,1:klev,jrow) = plu(1:kproma,1:klev)
  ! op_mm_20140327+
  ! replace with 2d object 
  cv_cldwater_2d(1:kproma,1:klev) = plu(1:kproma,1:klev)
  ! op_mm_20140327-
  ! mz_ak_20051221-
  
  ! op_mm_20140226+
  ! add calculation of convective gust
  ! copied from src_conv_tiedtke (COSMO)
  ! set whole field to zero (as only convective parts get values
  pvgustcon(:)=0.0_dp
  DO jl = 1,kproma
     IF (ldcum(jl)) THEN
        ! correction of zvddraf (max. possible convective gust)
        IF (3600.0_dp*(prsfc(jl)+pssfc(jl)) <= 0.015_dp) zvddraf(jl)=0.0_dp ! op_ff_20161020
        pvgustcon(jl)= MAX(pvgustcon(jl),zvddraf(jl))
     END IF
  END DO
  ! op_mm_20140226-
     
  RETURN
END SUBROUTINE tiedtke_cumastr


!============================================================================
!============================================================================

SUBROUTINE tiedtke_cumastrt( kproma, kbdim, klev, klevp1, klevm1,       &
           jrow,     ztmst,    nn,       ktrac,    ilab,                &
           pten,     pqen,     pxen,     puen,     pven,                &
           ptven,    ldland,                                            &
           pxten,    pxtu,     pxtte,                                   &
           pverv,    pqsen,    pqhfla,                                  &
           paphp1,   pgeo,                                              &
           ptte,     pqte,     pvom,     pvol,                          &
           prsfc,    pssfc,    paprc,    paprs,    pxtec,               &
           pqtec,                                                       &
           ldcum,    ktype,    kcbot,    kctop,                         &
           pmfu,     pmfd,     prain,    ptracconv, delta_time,         &
           pvgustcon, pdtke_con )  ! op_mm_20140327 add pvgustcon and pdtke_con
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
!     MODEL PARAMETERS (DEFINED IN MODULE messy_convection_tiedtke_cumulus)
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

!
    USE MESSY_CONVECT_TIEDTKE

    IMPLICIT NONE

INTEGER, INTENT(IN) :: kbdim, klev, ktrac, kproma, klevp1, klevm1, jrow, nn

INTEGER  :: jl, jk, jt, ptracconv(ktrac)

REAL(dp), INTENT(IN) :: ztmst, delta_time

REAL(dp) :: pten(kbdim,klev),        pqen(kbdim,klev),                   &
            pxen(kbdim,klev),        ptven(kbdim,klev),                  &
            puen(kbdim,klev),        pven(kbdim,klev),                   &
            ptte(kbdim,klev),        pqte(kbdim,klev),                   &
            pvom(kbdim,klev),        pvol(kbdim,klev),                   &
            pqsen(kbdim,klev),       pgeo(kbdim,klev),                   &
            paphp1(kbdim,klevp1),                                        &
            pverv(kbdim,klev)
REAL(dp) :: ptu(kbdim,klev),         pqu(kbdim,klev),                    &
            plu(kbdim,klev),         plude(kbdim,klev),                  &
            pmfu(kbdim,klev),        pmfd(kbdim,klev),                   &
            paprc(kbdim),            paprs(kbdim),                       &
            prsfc(kbdim),            pssfc(kbdim),                       &
            prain(kbdim),            pqhfla(kbdim),                      &
            pvgustcon(kbdim) ! op_mm_20140226 calculation of convective gust
INTEGER  :: kcbot(kbdim),            kctop(kbdim),                       &
            ktype(kbdim)
REAL(dp) :: pxtec(kbdim,klev),       pqtec(kbdim,klev),                  &
            pqude(kbdim,klev)
REAL(dp) :: ztenh(kbdim,klev),       zqenh(kbdim,klev),                  &
            zxenh(kbdim,klev),                                           &
            zgeoh(kbdim,klev),       zqsenh(kbdim,klev),                 &
            ztd(kbdim,klev),         zqd(kbdim,klev),                    &
            zmfus(kbdim,klev),       zmfds(kbdim,klev),                  &
            zmfuq(kbdim,klev),       zmfdq(kbdim,klev),                  &
            zdmfup(kbdim,klev),      zdmfdp(kbdim,klev),                 &
            zmful(kbdim,klev),       zrfl(kbdim),                        &
            zuu(kbdim,klev),         zvu(kbdim,klev),                    &
            zud(kbdim,klev),         zvd(kbdim,klev)
REAL(dp) :: zcpen(kbdim,klev),       zcpcu(kbdim,klev)
REAL(dp) :: zentr(kbdim),            zhcbase(kbdim),                     &
            zmfub(kbdim),            zmfub1(kbdim),                      &
            zdqpbl(kbdim),           zdqcv(kbdim),                       &
            zvddraf(kbdim)   ! op_mm_20140226 for calculation of convective gust
REAL(dp) :: zsfl(kbdim),             zdpmel(kbdim,klev)
INTEGER  :: ilab(kbdim,klev),        idtop(kbdim),                       &
            ictop0(kbdim),           ilwmin(kbdim)
REAL(dp) :: pxten(kbdim,klev,ktrac), pxtte(kbdim,klev,ktrac),            &
            pxtu(kbdim,klev,ktrac),                                      &
            zxtenh(kbdim,klev,ktrac),zxtd(kbdim,klev,ktrac),             &
            zmfuxt(kbdim,klev,ktrac),zmfdxt(kbdim,klev,ktrac)

! op_mm_20140227 calculation of tke due to convection
REAL(dp) :: pdtke_con(kbdim,klev)  

LOGICAL  :: loddraf(kbdim),          ldland(kbdim)
LOGICAL  :: ldcum(kbdim)
LOGICAL  :: llo1, lo
!
REAL(dp) :: zcons2, ztau,  zqumqe, zdqmin, zmfmax, zalvs,  zalvdcp, zqalv, &
            zhsat,  zqsat, zes,    zcor,   zqst1,  zdqsdt, zgam,    zzz,   &
            zhhat,  zfac,  zeps,   zpbmpt
REAL(dp) :: zpmfun(kbdim,klev) ! op_ck_20031001
INTEGER  :: ikb, it, it1, icum, itopm2
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
  zcons2=1./(g*ztmst)
!
!---------------------------------------------------------------------
!
!*    2.           INITIALIZE VALUES AT VERTICAL GRID POINTS IN 'CUINI'
!                  ---------------------------------------------------
!
200 CONTINUE
  call tiedtke_cuini(kproma,   kbdim,    klev,     klevp1,   klevm1,  &
             pten,     pqen,     pqsen,    pxen,     puen,     pven,  &
             ptven,    ktrac,                                         &
             pxten,    zxtenh,   pxtu,     zxtd,     zmfuxt,   zmfdxt,&
             pverv,    pgeo,     paphp1,   zgeoh,                     &
             ztenh,    zqenh,    zqsenh,   zxenh,    ilwmin,          &
             ptu,      pqu,      ztd,      zqd,                       &
             zuu,      zvu,      zud,      zvd,                       &
             pmfu,     pmfd,     zmfus,    zmfds,                     &
             zmfuq,    zmfdq,    zdmfup,   zdmfdp,                    &
             zcpen,    zcpcu,                                         &
             zdpmel,   plu,      plude,    pqude,    ilab,            &
             zpmfun)   ! op_ck_20031001/mz_ht_20040318 for posdef update
  IF (lookupoverflow) RETURN
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
  CALL tiedtke_cubase(kproma, kbdim, klev, klevp1, klevm1,            &
              ztenh,    zqenh,    zgeoh,    paphp1,                   &
              ptu,      pqu,      plu,                                &
              puen,     pven,     zuu,      zvu,                      &
              zcpcu,                                                  &
              ldcum,    kcbot,    ilab)
  IF (lookupoverflow) RETURN
!
!*             (B) DETERMINE TOTAL MOISTURE CONVERGENCE AND
!*                 THEN DECIDE ON TYPE OF CUMULUS CONVECTION
!                  -----------------------------------------
!
  jk=1
  DO 310 jl=1,kproma
     zdqpbl(jl)=0.0
     zdqcv(jl)=pqte(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))
     idtop(jl)=0
310 END DO
  DO 320 jk=2,klev
     DO 315 jl=1,kproma
        zdqcv(jl)=zdqcv(jl)+pqte(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))
        IF(jk.GE.kcbot(jl)) zdqpbl(jl)=zdqpbl(jl)+pqte(jl,jk)         &
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
     zdqmin=MAX(0.01*zqenh(jl,ikb),1.e-10_dp)
     llo1=zdqpbl(jl).GT.0..AND.zqumqe.GT.zdqmin.AND.ldcum(jl)
     zmfub(jl)=MERGE(zdqpbl(jl)/(g*MAX(zqumqe,zdqmin)),0.01_dp,llo1)
     zmfmax=(paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
     zmfub(jl)=MIN(zmfub(jl),zmfmax)
     IF(.NOT.llo1) ldcum(jl)=.FALSE.
     ktype(jl)=MERGE(1,2,zdqcv(jl).GT.MAX(0._dp,-1.1*pqhfla(jl)*g))
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
        zqalv=1./zalvs
        zhsat=zcpcu(jl,jk)*ztenh(jl,jk)+zgeoh(jl,jk)+zalvs*zqsenh(jl,jk)
        it = NINT(ztenh(jl,jk)*1000.)
! ka_sv_20170406+
! ! limiting to below thermosphere
!!$       IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
          IF ((it<jptlucu1 .OR. it>jptlucu2) .AND. &
               (paphp1(jl,jk+1) >= 1.0_dp)) lookupoverflow = .TRUE.
! ka_sv_20170406-
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
        zdqsdt=(zqst1-zqsat)*1000.
        zgam=MERGE(zalvdcp*zdqsdt,zqsat*zcor*tlucub(it),LO)
        zzz=zcpcu(jl,jk)*ztenh(jl,jk)*vtmpc1
        zhhat=zhsat-(zzz+zgam*zzz)/(1._dp+zgam*zzz*zqalv)*               &
                       MAX(zqsenh(jl,jk)-zqenh(jl,jk),0._dp)
        IF(jk.LT.ictop0(jl).AND.zhcbase(jl).GT.zhhat) ictop0(jl)=jk
420  END DO
430 END DO
!
     IF (lookupoverflow) RETURN
!!
!!     DEEP CONVECTION IF CLOUD DEPTH > 200 HPA, ELSE SHALLOW
!!     (CLOUD DEPTH FROM NON-ENTRAINIG PLUME)
!!
!  DO jl=1,kproma
!     ktype(jl)=MERGE(1,2,                                             &
!                paphp1(jl,kcbot(jl))-paphp1(jl,ictop0(jl)).gt.2.E4)
!     zentr(jl)=MERGE(entrpen,entrscv,ktype(jl).eq.1)
!  ENDDO
!!
!*             (B) DO ASCENT IN 'CUASCT' IN ABSENCE OF DOWNDRAFTS
!                  ----------------------------------------------
!
  CALL tiedtke_cuasct(kproma, kbdim, klev, klevp1, klevm1, jrow,ztmst,& !mz_ht_20040317+ added ztmst, jrow
             ztenh,    zqenh,    puen,     pven,                      &
             ktrac,                                                   &
             zxtenh,   pxten,    pxtu,     zmfuxt,                    &
             pten,     pqen,     pqsen,                               &
             pgeo,     zgeoh,    paphp1,                              &
             pqte,     pverv,    ilwmin,                              &
             ldcum,    ldland,   ktype,    ilab,                      &
             ptu,      pqu,      plu,      zuu,      zvu,             &
             pmfu,     zmfub,    zentr,                               &
             zmfus,    zmfuq,                                         &
             zmful,    plude,    pqude,    zdmfup,                    &
             zcpen,    zcpcu,                                         &
             kcbot,    kctop,    ictop0,   icum,                      &
             zpmfun)   ! op_ck_20031001/mz_ht_20040318 for posdef update
  IF (lookupoverflow) RETURN
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
     call tiedtke_cudlfs(kproma,   kbdim,    klev,     klevp1,        &
                 ztenh,    zqenh,    puen,     pven,                  &
                 ktrac,                                               &
                 zxtenh,   pxtu,     zxtd,     zmfdxt,                &
                 zgeoh,    paphp1,                                    &
                 ptu,      pqu,      zuu,      zvu,                   &
                 ldcum,    kcbot,    kctop,    zmfub,    zrfl,        &
                 ztd,      zqd,      zud,      zvd,                   &
                 pmfd,     zmfds,    zmfdq,    zdmfdp,                &
                 zcpcu,                                               &
                 idtop,    loddraf)
!
!*            (B)  DETERMINE DOWNDRAFT T,Q AND FLUXES IN 'CUDDRAF'
!                 -----------------------------------------------
!
     call tiedtke_cuddraf(kproma,   kbdim,    klev,     klevp1, jrow, &
                  ztenh,    zqenh,    puen,     pven,                 &
                  ktrac,                                              &
                  zxtenh,   zxtd,     zmfdxt,                         &
                  zgeoh,    paphp1,   zrfl,                           &
                  ztd,      zqd,      zud,      zvd,                  &
                  pmfd,     zmfds,    zmfdq,    zdmfdp,               &
                  zcpcu,                                              &
                  loddraf,                                            &
                  pxten,  &  ! op_ck_20031001/mz_ht_20040318   for pos_def update
                  zvddraf )  ! op_mm_20140226 conective gust 
!
!*            (C)  RECALCULATE CONVECTIVE FLUXES DUE TO EFFECT OF
!*                 DOWNDRAFTS ON BOUNDARY LAYER MOISTURE BUDGET
!                  --------------------------------------------
!
!DIR$ IVDEP
  DO 520 jl=1,kproma
     IF(loddraf(jl)) THEN
        ikb=kcbot(jl)
        llo1=pmfd(jl,ikb).LT.0.
        zeps=MERGE(cmfdeps,0._dp,llo1)
        zqumqe=pqu(jl,ikb)+plu(jl,ikb)-                               &
                    zeps*zqd(jl,ikb)-(1._dp-zeps)*zqenh(jl,ikb)
        zdqmin=MAX(0.01*zqenh(jl,ikb),1.e-10_dp)
        zmfmax=(paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
        llo1=zdqpbl(jl).GT.0..AND.zqumqe.GT.zdqmin.AND.ldcum(jl)      &
                    .AND.zmfub(jl).LT.zmfmax
        zmfub1(jl)=MERGE(zdqpbl(jl)/(g*MAX(zqumqe,zdqmin)),           &
                               zmfub(jl),llo1)
        zmfub1(jl)=MERGE(zmfub1(jl),zmfub(jl),                        &
                           (ktype(jl).EQ.1.OR.ktype(jl).EQ.2) .AND.   &
                           ABS(zmfub1(jl)-zmfub(jl)).LT.0.2*zmfub(jl))
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
     DO 5304 jt=1,ktrac
        DO 5302 jl=1,kproma
           IF(loddraf(jl)) THEN
              zfac=zmfub1(jl)/MAX(zmfub(jl),1.e-10_dp)
              zmfdxt(jl,jk,jt)=zmfdxt(jl,jk,jt)*zfac
           ENDIF
5302    END DO
5304 END DO
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
  CALL tiedtke_cuasct(kproma, kbdim, klev, klevp1, klevm1, jrow,ztmst,& !mz_ht_20040317+ added ztmst, jrow
             ztenh,    zqenh,    puen,     pven,                      &
             ktrac,                                                   &
             zxtenh,   pxten,    pxtu,     zmfuxt,                    &
             pten,     pqen,     pqsen,                               &
             pgeo,     zgeoh,    paphp1,                              &
             pqte,     pverv,    ilwmin,                              &
             ldcum,    ldland,   ktype,    ilab,                      &
             ptu,      pqu,      plu,      zuu,      zvu,             &
             pmfu,     zmfub,    zentr,                               &
             zmfus,    zmfuq,                                         &
             zmful,    plude,    pqude,    zdmfup,                    &
             zcpen,    zcpcu,                                         &
             kcbot,    kctop,    ictop0,   icum,                      &
             zpmfun)   ! op_ck_20031001/mz_ht_20040318 for posdef update
  IF (lookupoverflow) RETURN
!
!---------------------------------------------------------------------
!*    7.0          DETERMINE FINAL CONVECTIVE FLUXES IN 'CUFLX'
!                  ------------------------------------------
!
700 CONTINUE
  call tiedtke_cuflx(kproma,   kbdim,    klev,     klevp1, jrow,      &
             pqen,     pqsen,    ztenh,    zqenh,                     &
             ktrac,                                                   &
             zxtenh,   zmfuxt,   zmfdxt,                              &
             paphp1,   zgeoh,                                         &
             kcbot,    kctop,    idtop,                               &
             ktype,    loddraf,  ldcum,                               &
             pmfu,     pmfd,     zmfus,    zmfds,                     &
             zmfuq,    zmfdq,    zmful,                               &
             zdmfup,   zdmfdp,   zrfl,     prain,                     &
             zcpcu,                                                   &
             zcpen,                                                   &  ! mim_sb_20090917
             pten,     zsfl,     zdpmel,   itopm2,  ztmst,            &
             zpmfun,   pxten)            ! mz_ht_20040318+ needed for posdef update
!
!
!---------------------------------------------------------------------
!
!*    8.0          UPDATE TENDENCIES FOR T AND Q IN SUBROUTINE CUDTDQ
!                  --------------------------------------------------
!
800 CONTINUE
  call tiedtke_cudtdq(kproma, kbdim, klev, klevp1, itopm2, ldcum, ktrac, &
              jrow,     paphp1,   pten,     ptte,     pqte,              &
              pxtte,    pxtec,    zmfuxt,   zmfdxt,                      &
              zmfus,    zmfds,    zmfuq,    zmfdq,                       &
              zmful,    zdmfup,   zdmfdp,   plude,                       &
              zdpmel,   zrfl,     zsfl,                                  &
              zcpen,    pqtec,    pqude,                                 &
              prsfc,    pssfc,    paprc,    paprs,                       &
              delta_time, ptracconv, ztmst, pxten, pdtke_con,            &
              zqenh, ztenh)
!
!---------------------------------------------------------------------
!
!*    9.0          UPDATE TENDENCIES FOR U AND U IN SUBROUTINE CUDUDV
!                  --------------------------------------------------
!
900 CONTINUE
  IF(lmfdudv) THEN
     call tiedtke_cududv(kproma,   kbdim,    klev,     klevp1,        &
                 itopm2,   ktype,    kcbot,    paphp1,   ldcum,       &
                 puen,     pven,     pvom,     pvol,                  &
                 zuu,      zud,      zvu,      zvd,                   &
                 pmfu,     pmfd)
!
  END IF
!
1000 CONTINUE
!
!
 ! op_mm_20140226+
 ! calculation of convective gust
 ! copied from src_conv_tiedtke (COSMO)
 ! set whole field to zero (as only convective parts get values)
 pvgustcon(:)=0.0_dp
 DO jl = 1,kproma
      IF (ldcum(jl)) THEN
        ! correction of zvddraf (max. possible convective gust)
        IF (3600.0_dp*(prsfc(jl)+pssfc(jl)) <= 0.015) zvddraf(jl)=0.0
        pvgustcon(jl)= MAX(pvgustcon(jl),zvddraf(jl))
     END IF
  END DO
 ! op_mm_20140226-

  RETURN
END SUBROUTINE tiedtke_cumastrt


!=========================================================================
!=========================================================================


SUBROUTINE tiedtke_cumastrh(kproma, kbdim, klev, klevp1, klevm1,        &
           jrow,     ztmst,    nn,       ktrac,    ilab,                &
           pten,     pqen,     pxen,     puen,     pven,                &
           ptven,    ldland,                                            &
           pxten,    pxtu,     pxtte,                                   &
           pverv,    pqsen,    pqhfla,                                  &
           paphp1,   pgeo,                                              &
           ptte,     pqte,     pvom,     pvol,                          &
           prsfc,    pssfc,    paprc,    paprs,    pxtec,               &
           pqtec,                                                       &
           ldcum,    ktype,    kcbot,    kctop,                         &
           pmfu,     pmfd,     prain,    ptracconv, delta_time,         &
           pvgustcon, pdtke_con) ! op_mm_20140327 added pvgustcon and pdtke_con
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
!       CUASCT:  CLOUD ASCENT FOR ENTRAINING PLUME
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
!     MODEL PARAMETERS (DEFINED IN MODULE messy_convection_tiedtke_cumulus)
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
    USE MESSY_CONVECT_TIEDTKE

    IMPLICIT NONE


INTEGER, INTENT(IN) :: kbdim, klev, ktrac, kproma, klevp1, klevm1, jrow, nn

INTEGER  :: jl, jk, jt, ptracconv(ktrac)

REAL(dp), INTENT(IN) :: ztmst, delta_time

REAL(dp) :: pten(kbdim,klev),        pqen(kbdim,klev),                     &
            pxen(kbdim,klev),        ptven(kbdim,klev),                    &
            puen(kbdim,klev),        pven(kbdim,klev),                     &
            ptte(kbdim,klev),        pqte(kbdim,klev),                     &
            pvom(kbdim,klev),        pvol(kbdim,klev),                     &
            pqsen(kbdim,klev),       pgeo(kbdim,klev),                     &
            paphp1(kbdim,klevp1),                                          &
            pverv(kbdim,klev)
REAL(dp) :: ptu(kbdim,klev),         pqu(kbdim,klev),                      &
            plu(kbdim,klev),         plude(kbdim,klev),                    &
            pmfu(kbdim,klev),        pmfd(kbdim,klev),                     &
            paprc(kbdim),            paprs(kbdim),                         &
            prsfc(kbdim),            pssfc(kbdim),                         &
            prain(kbdim),            pqhfla(kbdim)
INTEGER  :: kcbot(kbdim),            kctop(kbdim),                         &
            ktype(kbdim)
REAL(dp) :: pxtec(kbdim,klev),       pqtec(kbdim,klev),                    &
            pqude(kbdim,klev)
REAL(dp) :: ztenh(kbdim,klev),       zqenh(kbdim,klev),                    &
            zxenh(kbdim,klev),                                             &
            zgeoh(kbdim,klev),       zqsenh(kbdim,klev),                   &
            ztd(kbdim,klev),         zqd(kbdim,klev),                      &
            zmfus(kbdim,klev),       zmfds(kbdim,klev),                    &
            zmfuq(kbdim,klev),       zmfdq(kbdim,klev),                    &
            zdmfup(kbdim,klev),      zdmfdp(kbdim,klev),                   &
            zmful(kbdim,klev),       zrfl(kbdim),                          &
            zuu(kbdim,klev),         zvu(kbdim,klev),                      &
            zud(kbdim,klev),         zvd(kbdim,klev)
REAL(dp) :: zcpen(kbdim,klev),       zcpcu(kbdim,klev)
! op_mm_20140226+
! for calculation of convective gust
REAL(dp) :: zvddraf(kbdim),          pvgustcon(kbdim)
! op_mm_20140226-

REAL(dp) :: zentr(kbdim),            zhcbase(kbdim),                       &
            zmfub(kbdim),            zmfub1(kbdim),                        &
            zdqpbl(kbdim),           zdqcv(kbdim)
REAL(dp) :: zsfl(kbdim),             zdpmel(kbdim,klev)
REAL(dp) :: zcape(kbdim),            zheat(kbdim)
INTEGER  :: ilab(kbdim,klev),        idtop(kbdim),                         &
            ictop0(kbdim),           ilwmin(kbdim)
REAL(dp) :: pxten(kbdim,klev,ktrac), pxtte(kbdim,klev,ktrac),              &
            pxtu(kbdim,klev,ktrac),                                        &
            zxtenh(kbdim,klev,ktrac),zxtd(kbdim,klev,ktrac),               &
            zmfuxt(kbdim,klev,ktrac),zmfdxt(kbdim,klev,ktrac)
REAL(dp) :: pdtke_con(kbdim,klev)  ! op_m_20140227 tke due to convection
LOGICAL  :: loddraf(kbdim),          ldland(kbdim)
LOGICAL  :: ldcum(kbdim)
LOGICAL  :: llo1, lo
!
REAL(dp) :: zcons2, ztau,  zqumqe, zdqmin, zmfmax, zalvs,  zalvdcp, zqalv, &
            zhsat,  zqsat, zes,    zcor,   zqst1,  zdqsdt, zgam,    zzz,   &
            zhhat,  zro,   zdz,    zfac,  zeps,   zpbmpt
REAL(dp) :: zpmfun(kbdim,klev) ! op_ck_20031001
INTEGER  :: ikb, it, it1, icum, itopm2

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
  zcons2=1./(g*ztmst)
  ztau=MIN(3.*3600._dp,7200._dp*63./REAL(nn,dp))
!
!-----------------------------------------------------------------------
!
!*    2.           INITIALIZE VALUES AT VERTICAL GRID POINTS IN 'CUINI'
!                  ---------------------------------------------------
!
200 CONTINUE
  call tiedtke_cuini(kproma,   kbdim,    klev,     klevp1,   klevm1,    &
             pten,     pqen,     pqsen,    pxen,     puen,     pven,    &
             ptven,    ktrac,                                           &
             pxten,    zxtenh,   pxtu,     zxtd,     zmfuxt,   zmfdxt,  &
             pverv,    pgeo,     paphp1,   zgeoh,                       &
             ztenh,    zqenh,    zqsenh,   zxenh,    ilwmin,            &
             ptu,      pqu,      ztd,      zqd,                         &
             zuu,      zvu,      zud,      zvd,                         &
             pmfu,     pmfd,     zmfus,    zmfds,                       &
             zmfuq,    zmfdq,    zdmfup,   zdmfdp,                      &
             zcpen,    zcpcu,                                           &
             zdpmel,   plu,      plude,    pqude,    ilab,              &
             zpmfun)   ! op_ck_20031001/mz_ht_20040318 for posdef update
  IF (lookupoverflow) RETURN
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
  CALL tiedtke_cubase(kproma,   kbdim,    klev,     klevp1,   klevm1,   &
              ztenh,    zqenh,    zgeoh,    paphp1,                     &
              ptu,      pqu,      plu,                                  &
              puen,     pven,     zuu,      zvu,                        &
              zcpcu,                                                    &
              ldcum,    kcbot,    ilab)
  IF (lookupoverflow) RETURN
!
!*             (B) DETERMINE TOTAL MOISTURE CONVERGENCE AND
!*                 THEN DECIDE ON TYPE OF CUMULUS CONVECTION
!                  -----------------------------------------
!
  jk=1
  DO 310 jl=1,kproma
     zdqpbl(jl)=0.0
     zdqcv(jl)=pqte(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))
     idtop(jl)=0
310 END DO
  DO 320 jk=2,klev
     DO 315 jl=1,kproma
        zdqcv(jl)=zdqcv(jl)+pqte(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))
        IF(jk.GE.kcbot(jl)) zdqpbl(jl)=zdqpbl(jl)+pqte(jl,jk)           &
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
     zdqmin=MAX(0.01*zqenh(jl,ikb),1.e-10_dp)
     llo1=zdqpbl(jl).GT.0..AND.zqumqe.GT.zdqmin.AND.ldcum(jl)
     zmfub(jl)=MERGE(zdqpbl(jl)/(g*MAX(zqumqe,zdqmin)),0.01_dp,llo1)
     zmfmax=(paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
     zmfub(jl)=MIN(zmfub(jl),zmfmax)
     IF(.NOT.llo1) ldcum(jl)=.FALSE.
     ktype(jl)=MERGE(1,2,zdqcv(jl).GT.MAX(0._dp,-1.1*pqhfla(jl)*g))
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
     zhcbase(jl)=zcpcu(jl,ikb)*ptu(jl,ikb)+zgeoh(jl,ikb)+zalvs*pqu(jl,ikb)
     ictop0(jl)=kcbot(jl)-1
410 END DO
  DO 430 jk=klevm1,3,-1
!DIR$ IVDEP
     DO 420 jl=1,kproma
        zalvs=MERGE(alv,als,ztenh(jl,jk)>tmelt)
        zalvdcp=zalvs/zcpcu(jl,jk)
        zqalv=1./zalvs
        zhsat=zcpcu(jl,jk)*ztenh(jl,jk)+zgeoh(jl,jk)+zalvs*zqsenh(jl,jk)
        it = NINT(ztenh(jl,jk)*1000.)
! ka_sv_20170406+
! ! limiting to below thermosphere
!!$       IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
          IF ((it<jptlucu1 .OR. it>jptlucu2) .AND. &
               (paphp1(jl,jk+1) >= 1.0_dp)) lookupoverflow = .TRUE.
! ka_sv_20170406-
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
        zdqsdt=(zqst1-zqsat)*1000.
        zgam=MERGE(zalvdcp*zdqsdt,zqsat*zcor*tlucub(it),LO)
        zzz=zcpcu(jl,jk)*ztenh(jl,jk)*vtmpc1
        zhhat=zhsat-(zzz+zgam*zzz)/(1._dp+zgam*zzz*zqalv)*                 &
                          MAX(zqsenh(jl,jk)-zqenh(jl,jk),0._dp)
        IF(jk.LT.ictop0(jl).AND.zhcbase(jl).GT.zhhat) ictop0(jl)=jk
420  END DO
430 END DO
!
     IF (lookupoverflow) RETURN
!!
!!     DEEP CONVECTION IF CLOUD DEPTH > 200 HPA, ELSE SHALLOW
!!     (CLOUD DEPTH FROM NON-ENTRAINIG PLUME)
!!
!  DO jl=1,kproma
!     ktype(jl)=MERGE(1,2,                                               &
!                paphp1(jl,kcbot(jl))-paphp1(jl,ictop0(jl)).gt.2.E4)
!     zentr(jl)=MERGE(entrpen,entrscv,ktype(jl).eq.1)
!  ENDDO
!!
!*             (B) DO ASCENT IN 'CUASCT' IN ABSENCE OF DOWNDRAFTS
!                  ----------------------------------------------
!
  CALL tiedtke_cuasct(kproma, kbdim, klev, klevp1, klevm1,jrow, ztmst,  & !mz_ht_20040317+ added ztmst, jrow
             ztenh,    zqenh,    puen,     pven,                        &
             ktrac,                                                     &
             zxtenh,   pxten,    pxtu,     zmfuxt,                      &
             pten,     pqen,     pqsen,                                 &
             pgeo,     zgeoh,    paphp1,                                &
             pqte,     pverv,    ilwmin,                                &
             ldcum,    ldland,   ktype,    ilab,                        &
             ptu,      pqu,      plu,      zuu,      zvu,               &
             pmfu,     zmfub,    zentr,                                 &
             zmfus,    zmfuq,                                           &
             zmful,    plude,    pqude,    zdmfup,                      &
             zcpen,    zcpcu,                                           &
             kcbot,    kctop,    ictop0,   icum,                        &
             zpmfun)   ! op_ck_20031001/mz_ht_20040318 for posdef update
  IF (lookupoverflow) RETURN
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
     call tiedtke_cudlfs(kproma,   kbdim,    klev,     klevp1,          &
                 ztenh,    zqenh,    puen,     pven,                    &
                 ktrac,                                                 &
                 zxtenh,   pxtu,     zxtd,     zmfdxt,                  &
                 zgeoh,    paphp1,                                      &
                 ptu,      pqu,      zuu,      zvu,                     &
                 ldcum,    kcbot,    kctop,    zmfub,    zrfl,          &
                 ztd,      zqd,      zud,      zvd,                     &
                 pmfd,     zmfds,    zmfdq,    zdmfdp,                  &
                 zcpcu,                                                 &
                 idtop,    loddraf)
!
!*            (B)  DETERMINE DOWNDRAFT T,Q AND FLUXES IN 'CUDDRAF'
!                  -----------------------------------------------
!
     call tiedtke_cuddraf(kproma,   kbdim,    klev, klevp1, jrow,       &
                  ztenh,    zqenh,    puen,     pven,                   &
                  ktrac,                                                &
                  zxtenh,   zxtd,     zmfdxt,                           &
                  zgeoh,    paphp1,   zrfl,                             &
                  ztd,      zqd,      zud,      zvd,                    &
                  pmfd,     zmfds,    zmfdq,    zdmfdp,                 &
                  zcpcu,                                                &
                  loddraf,                                              &
                  pxten,  & ! op_ck_20031001/mz_ht_20040318   for pos_def update
                  zvddraf ) ! op_m_20140226 added zvddraf (convective gust) 
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
     zheat(jl)=0.
     zcape(jl)=0.
     zmfub1(jl)=zmfub(jl)
  ENDDO
!
  DO jk=1,klev
     DO jl=1,kproma
        llo1=ldcum(jl).AND.ktype(jl).EQ.1
        IF(llo1.AND.jk.LE.kcbot(jl).AND.jk.GT.kctop(jl)) THEN
           ikb=kcbot(jl)
           zro=paphp1(jl,jk)/(rd*ztenh(jl,jk)*(1._dp+vtmpc1*zqenh(jl,jk)))
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
  CAPE(1:kproma,jrow) = zcape(1:kproma)
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
        llo1=pmfd(jl,ikb).LT.0..AND.loddraf(jl)
        zeps=MERGE(cmfdeps,0._dp,llo1)
        zqumqe=pqu(jl,ikb)+plu(jl,ikb)-                                 &
                    zeps*zqd(jl,ikb)-(1._dp-zeps)*zqenh(jl,ikb)
        zdqmin=MAX(0.01*zqenh(jl,ikb),1.e-10_dp)
        zmfmax=(paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
        llo1=zdqpbl(jl).GT.0..AND.zqumqe.GT.zdqmin.AND.ldcum(jl)        &
                    .AND.zmfub(jl).LT.zmfmax
        zmfub1(jl)=MERGE(zdqpbl(jl)/(g*MAX(zqumqe,zdqmin)),             &
                               zmfub(jl),llo1)
        zmfub1(jl)=MERGE(zmfub1(jl),zmfub(jl),                          &
                             ABS(zmfub1(jl)-zmfub(jl)).LT.0.2*zmfub(jl))
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
     DO 5304 jt=1,ktrac
        DO 5302 jl=1,kproma
           IF(ldcum(jl)) THEN
              zfac=zmfub1(jl)/MAX(zmfub(jl),1.e-10_dp)
              zmfdxt(jl,jk,jt)=zmfdxt(jl,jk,jt)*zfac
           ENDIF
5302    END DO
5304 END DO
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
  CALL tiedtke_cuasct(kproma, kbdim, klev, klevp1, klevm1, jrow, ztmst, & !mz_ht_20040317+ added ztmst, jrow
             ztenh,    zqenh,    puen,     pven,                        &
             ktrac,                                                     &
             zxtenh,   pxten,    pxtu,     zmfuxt,                      &
             pten,     pqen,     pqsen,                                 &
             pgeo,     zgeoh,    paphp1,                                &
             pqte,     pverv,    ilwmin,                                &
             ldcum,    ldland,   ktype,    ilab,                        &
             ptu,      pqu,      plu,      zuu,      zvu,               &
             pmfu,     zmfub,    zentr,                                 &
             zmfus,    zmfuq,                                           &
             zmful,    plude,    pqude,    zdmfup,                      &
             zcpen,    zcpcu,                                           &
             kcbot,    kctop,    ictop0,   icum,                        &
             zpmfun)   ! op_ck_20031001/mz_ht_20040318 for posdef update
  IF (lookupoverflow) RETURN
!
!-----------------------------------------------------------------------
!*    7.0          DETERMINE FINAL CONVECTIVE FLUXES IN 'CUFLX'
!                  ------------------------------------------
!
700 CONTINUE
  call tiedtke_cuflx(kproma,   kbdim,    klev,     klevp1, jrow,        &
             pqen,     pqsen,    ztenh,    zqenh,                       &
             ktrac,                                                     &
             zxtenh,   zmfuxt,   zmfdxt,                                &
             paphp1,   zgeoh,                                           &
             kcbot,    kctop,    idtop,                                 &
             ktype,    loddraf,  ldcum,                                 &
             pmfu,     pmfd,     zmfus,    zmfds,                       &
             zmfuq,    zmfdq,    zmful,                                 &
             zdmfup,   zdmfdp,   zrfl,     prain,                       &
             zcpcu,                                                     &
             zcpen,                                                   & ! mim_sb_20090917
             pten,     zsfl,     zdpmel,   itopm2, ztmst,               &
             zpmfun,   pxten)            ! mz_ht_20040318+ needed for posdef update
!
!-----------------------------------------------------------------------
!
!*    8.0          UPDATE TENDENCIES FOR T AND Q IN SUBROUTINE CUDTDQ
!                  --------------------------------------------------
!
800 CONTINUE
  call tiedtke_cudtdq(kproma, kbdim, klev, klevp1, itopm2, ldcum, ktrac,&
              jrow,     paphp1,   pten,     ptte,     pqte,             &
              pxtte,    pxtec,    zmfuxt,   zmfdxt,                     &
              zmfus,    zmfds,    zmfuq,    zmfdq,                      &
              zmful,    zdmfup,   zdmfdp,   plude,                      &
              zdpmel,   zrfl,     zsfl,                                 &
              zcpen,    pqtec,    pqude,                                &
              prsfc,    pssfc,    paprc,    paprs,                      &
              delta_time, ptracconv, ztmst, pxten, pdtke_con,           &
              zqenh, ztenh) ! op_mm_20140327 added pdtke_con, zqenh, ztenh
!
!-----------------------------------------------------------------------
!
!*    9.0          UPDATE TENDENCIES FOR U AND U IN SUBROUTINE CUDUDV
!                  --------------------------------------------------
!
900 CONTINUE
  IF(lmfdudv) THEN
     call tiedtke_cududv(kproma,   kbdim,    klev,     klevp1,          &
                 itopm2,   ktype,    kcbot,    paphp1,   ldcum,         &
                 puen,     pven,     pvom,     pvol,                    &
                 zuu,      zud,      zvu,      zvd,                     &
                 pmfu,     pmfd)
!
  END IF
!
1000 CONTINUE
!

  ! op_mm_20140226+
  ! calculation of convective gust
  ! copied from src_conv_tiedtke (COSMO)
  ! set whole field to zero (as only convective parts get values
  pvgustcon(:)=0.0_dp
  DO jl = 1,kproma
     IF (ldcum(jl)) THEN
        ! correction of zvddraf (max. possible convective gust)
        IF (3600.0_dp*(prsfc(jl)+pssfc(jl)) <= 0.015) zvddraf(jl)=0.0
        pvgustcon(jl)= MAX(pvgustcon(jl),zvddraf(jl))
     END IF
  END DO
 ! op_mm_20140226-

  RETURN
END SUBROUTINE tiedtke_cumastrh
!=========================================================================

SUBROUTINE ECMWF_cumastr( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,   &
                          LDLAND,   PTSPHY,                     &
                          PTEN,     PQEN,     PUEN,     PVEN,     PLITOT, &
                          PVERVEL,  PQHFL,    PAHFS,    psstru,   psstrv, &
                          PAP,      PAPH,     PGEO,     PGEOH,            &
                          PTENT,    PTENQ,    PTENU,    PTENV,            &
                          PTENL,    PTENI,                                &
                          LDCUM,    KTYPE,    KCBOT,    KCTOP,            &
                          KBOTSC,   LDSC,                                 &
                          PTU,      PQU,      PLU,                        &
                          PMFLXR,   PMFLXS,                               &
                          PRAIN,                                          &
                          PMFU,     PMFD,                                 &
                          PMFUDE_RATE,        PMFDDE_RATE,                &
                          PMFUEN_RATE,        PMFDEN_RATE,                &
                          PCAPE,                                          &
                          KTRAC,    PCEN,     PTENC,                      &
                          lwc,      iwc,      rform,    sform)

!**** *CUMASTR*  MASTER ROUTINE FOR CUMULUS MASSFLUX-SCHEME

!     M.TIEDTKE      E.C.M.W.F.     1986/1987/1989
!     D.GREGORY      E.C.M.W.F.     1996
!     P.BECHTOLD     E.C.M.W.F.     2004

!     PURPOSE
!     -------

!          THIS ROUTINE COMPUTES THE PHYSICAL TENDENCIES OF THE
!     PROGNOSTIC VARIABLES T,Q,U, V AND TRACERS DUE TO CONVECTIVE PROCESSES.
!     PROCESSES CONSIDERED ARE: CONVECTIVE FLUXES, FORMATION OF
!     PRECIPITATION, EVAPORATION OF FALLING RAIN BELOW CLOUD BASE,
!     SATURATED CUMULUS DOWNDRAFTS.

!**   INTERFACE.
!     ----------

!          *CUMASTR* IS CALLED FROM *CUCALL*
!     THE ROUTINE TAKES ITS INPUT FROM THE LONG-TERM STORAGE
!     T,Q,U,V,PHI AND P AND MOISTURE TENDENCIES.
!     IT RETURNS ITS OUTPUT TO THE SAME SPACE
!      1.MODIFIED TENDENCIES OF MODEL VARIABLES
!      2.RATES OF CONVECTIVE PRECIPITATION
!        (USED IN SUBROUTINE SURF)
!      3.CLOUD BASE, CLOUD TOP AND PRECIP FOR RADIATION
!        (USED IN SUBROUTINE CLOUD)

!     METHOD
!     -------

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

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KTDIA*        START OF THE VERTICAL LOOP
!    *KLEV*         NUMBER OF LEVELS
!    *KTRAC*        NUMBER OF CHEMICAL TRACERS
!    *KSTEP*        CURRENT TIME STEP INDEX
!    *KSTART*       FIRST STEP OF MODEL

!     INPUT PARAMETERS (LOGICAL)

!    *LDLAND*       LAND SEA MASK (.TRUE. FOR LAND)

!     INPUT PARAMETERS (REAL)

!    *PTSPHY*       TIME STEP FOR THE PHYSICS                       S
!    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
!    *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1)  KG/KG
!    *PUEN*         PROVISIONAL ENVIRONMENT U-VELOCITY (T+1)       M/S
!    *PVEN*         PROVISIONAL ENVIRONMENT V-VELOCITY (T+1)       M/S
!    *PCEN*         PROVISIONAL ENVIRONMENT TRACER CONCENTRATIONS KG/KG 
!    *PLITOT*       GRID MEAN LIQUID WATER+ICE CONTENT            KG/KG
!    *PVERVEL*      VERTICAL VELOCITY                             PA/S
!    *PQSEN*        ENVIRONMENT SPEC. SATURATION HUMIDITY (T+1)   KG/KG
!    *PQHFL*        MOISTURE FLUX (EXCEPT FROM SNOW EVAP.)        KG/(SM2)
!    *PAHFS*        SENSIBLE HEAT FLUX                            W/M2
!    *PAP*          PROVISIONAL PRESSURE ON FULL LEVELS             PA
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS             PA
!    *PGEO*         GEOPOTENTIAL                                  M2/S2
!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2
!    *PSSTRU*       surface momentum flux u                - not used presently
!    *PSSTRV*       surface momentum flux v                - not used presently 

!    UPDATED PARAMETERS (REAL):

!    *PTENT*        TEMPERATURE TENDENCY                           K/S
!    *PTENQ*        MOISTURE TENDENCY                             KG/(KG S)
!    *PTENL*        LIQUID WATER TENDENCY                         KG/(KG S)
!    *PTENI*        ICE CONDENSATE TENDENCY                       KG/(KG S)
!    *PTENU*        TENDENCY OF U-COMP. OF WIND                    M/S2
!    *PTENV*        TENDENCY OF V-COMP. OF WIND                    M/S2
!    *PTENC*        TENDENCY OF CHEMICAL TRACERS                   1/S

!    OUTPUT PARAMETERS (LOGICAL):

!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS 
!    *LDSC*         FLAG: .TRUE. FOR SC-POINTS

!    OUTPUT PARAMETERS (INTEGER):

!    *KTYPE*        TYPE OF CONVECTION
!                       1 = PENETRATIVE CONVECTION
!                       2 = SHALLOW CONVECTION
!                       3 = MIDLEVEL CONVECTION
!    *KCBOT*        CLOUD BASE LEVEL
!    *KCTOP*        CLOUD TOP LEVEL
!    *KBOTSC*       CLOUD BASE LEVEL FOR SC-CLOUDS

!    OUTPUT PARAMETERS (REAL):

!    *PTU*          TEMPERATURE IN UPDRAFTS                         K
!    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                    KG/KG
!    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS              KG/KG
!    *PLUDE*        DETRAINED LIQUID WATER                        KG/(M2*S)
!    *PENTH*        INCREMENT OF DRY STATIC ENERGY                 J/(KG*S)
!    *PMFLXR*       CONVECTIVE RAIN FLUX                          KG/(M2*S)
!    *PMFLXS*       CONVECTIVE SNOW FLUX                          KG/(M2*S)
!    *PRAIN*        TOTAL PRECIP. PRODUCED IN CONV. UPDRAFTS      KG/(M2*S)
!                   (NO EVAPORATION IN DOWNDRAFTS)
!    *PMFU*         MASSFLUX UPDRAFTS                             KG/(M2*S)
!    *PMFD*         MASSFLUX DOWNDRAFTS                           KG/(M2*S)
!    *PMFUDE_RATE*  UPDRAFT DETRAINMENT RATE                      KG/(M3*S) modified to KG/(M2*S)
!    *PMFDDE_RATE*  DOWNDRAFT DETRAINMENT RATE                    KG/(M3*S) modified to KG/(M2*S)
!    *PMFUEN_RATE*  UPDRAFT ENTRAINMENT RATE                      KG/(M3*S) modified to KG/(M2*S)
!    *PMFUEN_RATE*  DOWNDRAFT ENTRAINMENT RATE                    KG/(M3*S) modified to KG/(M2*S)
!    *CAPE*         CONVECTVE AVAILABLE POTENTIAL ENERGY           J/KG

!     EXTERNALS.
!     ----------

!       CUINI:  INITIALIZES VALUES AT VERTICAL GRID USED IN CU-PARAMETR.
!       CUBASE: CLOUD BASE CALCULATION FOR PENETR.AND SHALLOW CONVECTION
!       CUASC:  CLOUD ASCENT FOR ENTRAINING PLUME
!       CUDLFS: DETERMINES VALUES AT LFS FOR DOWNDRAFTS
!       CUDDRAF:DOES MOIST DESCENT FOR CUMULUS DOWNDRAFTS
!       CUFLX:  FINAL ADJUSTMENTS TO CONVECTIVE FLUXES (ALSO IN PBL)
!       CUDQDT: UPDATES TENDENCIES FOR T AND Q
!       CUDUDV: UPDATES TENDENCIES FOR U AND V

!     SWITCHES.
!     --------

!          LMFPEN=.TRUE.   PENETRATIVE CONVECTION IS SWITCHED ON
!          LMFSCV=.TRUE.   SHALLOW CONVECTION IS SWITCHED ON
!          LMFMID=.TRUE.   MIDLEVEL CONVECTION IS SWITCHED ON
!          LMFDD=.TRUE.    CUMULUS DOWNDRAFTS SWITCHED ON
!          LMFDUDV=.TRUE.  CUMULUS FRICTION SWITCHED ON
!          LMFTRAC=.false. TRACER TRANSPORT (yet in user test phase)

!     MODEL PARAMETERS (DEFINED IN SUBROUTINE CUPARAM)
!     ------------------------------------------------
!     ENTRPEN    ENTRAINMENT RATE FOR PENETRATIVE CONVECTION
!     ENTRSCV    ENTRAINMENT RATE FOR SHALLOW CONVECTION
!     ENTRMID    ENTRAINMENT RATE FOR MIDLEVEL CONVECTION
!     ENTRDD     ENTRAINMENT RATE FOR CUMULUS DOWNDRAFTS
!     RMFCTOP    RELATIVE CLOUD MASSFLUX AT LEVEL ABOVE NONBUOYANCY LEVE
!     RMFCMAX    MAXIMUM MASSFLUX VALUE ALLOWED FOR
!     RMFCMIN    MINIMUM MASSFLUX VALUE (FOR SAFETY)
!     RMFDEPS    FRACTIONAL MASSFLUX FOR DOWNDRAFTS AT LFS
!     RPRCON     COEFFICIENT FOR CONVERSION FROM CLOUD WATER TO RAIN

!     REFERENCE.
!     ----------

!          PAPER ON MASSFLUX SCHEME (TIEDTKE,1989)
!          DRAFT PAPER ON MASSFLUX SCHEME (NORDENG, 1995)

!          MODIFICATIONS
!          -------------
!             92-09-21 : Update to Cy44      J.-J. MORCRETTE
!             96-03-12 : Introduce CAPE closure for deep convection
!                        (based upon previous work by Nordeng, 1995)
!             99-06-21 : Optimisation   D.Salmond
!             03-08-29 : Clean-up, deep/shallow switches  P.Bechtold
!             04-02-11 : Add tracer transport             P.Bechtold
!             05-02-11 : Positive scaling of total Mflux  P.Bechtold
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!----------------------------------------------------------------------

  USE MESSY_MAIN_CONSTANTS_MEM,   ONLY:  dp, RA => radius_earth
  USE MESSY_CONVECT_ECMWF_PARAM,  ONLY:  LHOOK,   DR_HOOK,                       &      !from YOMHOOK
                                         RG, RD, RCPD, RETV, RPI, RLVTT,         &
                                         RMFCFL, RMFSHCFL, RDEPTHS, LMFSCL_WSTAR,&
                                         RTAU, ENTRPEN, ENTRSCV,  NSMAX,         &  
                                         LMFPEN, LMFSCV, LMFDD, LMFDUDV, LMFTRAC,& 
                                         SATUR,                                  &      !subroutines
                                         FOELHMCU, FOEALFA                              !functions

  USE MESSY_CONVECT_ECMWF,        ONLY:  CUININ, CUBASEN, CUASCN, CUDLFSN,       &
                                         CUDDRAFN, CUFLXN, CUDTDQN, CUDUDV,      &  
                                         CUCTRACER
IMPLICIT NONE

INTEGER,INTENT(IN)    :: KLON
INTEGER,INTENT(IN)    :: KLEV
INTEGER,INTENT(IN)    :: KIDIA
INTEGER,INTENT(IN)    :: KFDIA
INTEGER,INTENT(IN)    :: KTRAC
INTEGER               :: KTDIA ! Argument NOT used
INTEGER               :: KSTEP ! Argument NOT used
INTEGER               :: KSTART ! Argument NOT used
LOGICAL           ,INTENT(IN)    :: LDLAND(KLON) 
REAL(dp)   ,INTENT(IN)    :: PTSPHY
REAL(dp)   ,INTENT(INOUT) :: PTEN(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PQEN(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PUEN(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PVEN(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PLITOT(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PVERVEL(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PQHFL(KLON,KLEV+1) 
REAL(dp)   ,INTENT(IN)    :: PAHFS(KLON,KLEV+1) 
REAL(dp)                  :: PSSTRU(KLON) ! Argument NOT used
REAL(dp)                  :: PSSTRV(KLON) ! Argument NOT used
REAL(dp)   ,INTENT(IN)    :: PAP(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1)
REAL(dp)   ,INTENT(IN)    :: PGEO(KLON,KLEV) 
REAL(dp)   ,INTENT(IN)    :: PGEOH(KLON,KLEV+1) 
REAL(dp)   ,INTENT(IN)    :: PCEN(KLON,KLEV,KTRAC) 
REAL(dp)   ,INTENT(INOUT) :: PTENT(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PTENQ(KLON,KLEV) 
REAL(dp)   ,INTENT(OUT)   :: PTENL(KLON,KLEV) 
REAL(dp)   ,INTENT(OUT)   :: PTENI(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PTENU(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PTENV(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PTENC(KLON,KLEV,KTRAC) 
LOGICAL           ,INTENT(INOUT) :: LDCUM(KLON)
INTEGER,INTENT(INOUT) :: KTYPE(KLON)
INTEGER,INTENT(INOUT) :: KCBOT(KLON)
INTEGER,INTENT(INOUT) :: KCTOP(KLON)
INTEGER,INTENT(OUT)   :: KBOTSC(KLON) 
LOGICAL           ,INTENT(OUT)   :: LDSC(KLON) 
REAL(dp)   ,INTENT(INOUT) :: PTU(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PQU(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PLU(KLON,KLEV) 
REAL(dp)   ,INTENT(INOUT) :: PMFLXR(KLON,KLEV+1) 
REAL(dp)   ,INTENT(INOUT) :: PMFLXS(KLON,KLEV+1) 
REAL(dp)   ,INTENT(OUT)   :: PRAIN(KLON) 
REAL(dp)   ,INTENT(INOUT) :: PMFU(KLON,KLEV)
REAL(dp)   ,INTENT(INOUT) :: PMFD(KLON,KLEV)
REAL(dp)   ,INTENT(INOUT) :: PMFUDE_RATE(KLON,KLEV)
REAL(dp)   ,INTENT(INOUT) :: PMFDDE_RATE(KLON,KLEV)
REAL(dp)   ,INTENT(INOUT) :: PMFUEN_RATE(KLON,KLEV)
REAL(dp)   ,INTENT(INOUT) :: PMFDEN_RATE(KLON,KLEV)
REAL(dp)   ,INTENT(OUT)   :: PCAPE(KLON) 
! mz_ht_20070622+
REAL(dp)   ,INTENT(INOUT) :: lwc(KLON,KLEV)
REAL(dp)   ,INTENT(INOUT) :: iwc(KLON,KLEV)
REAL(dp)   ,INTENT(INOUT) :: rform(KLON,KLEV)
REAL(dp)   ,INTENT(INOUT) :: sform(KLON,KLEV)
! mz_ht_20070622-

!*UPG change to operations
REAL(dp) :: PLUDE(KLON,KLEV) ! only local variable
REAL(dp) :: PENTH(KLON,KLEV) ! only local variable
REAL(dp) :: PQSEN(KLON,KLEV) ! only local variable
!*UPG
!     CHEMICAL TRACERS
!     3D DIAGNOSTICS FOR ERA40
!     OTHER CONVECTION DIAGNOSTICS
REAL(dp) ::     ZTENH(KLON,KLEV),       ZQENH(KLON,KLEV),&
 & ZQSENH(KLON,KLEV),&
 & ZTD(KLON,KLEV),         ZQD(KLON,KLEV),&
 & ZMFUS(KLON,KLEV),       ZMFDS(KLON,KLEV),&
 & ZMFUQ(KLON,KLEV),       ZMFDQ(KLON,KLEV),&
 & ZDMFUP(KLON,KLEV),      ZDMFDP(KLON,KLEV),&
 & ZMFUL(KLON,KLEV),       ZRFL(KLON),&
 & ZUU(KLON,KLEV),         ZVU(KLON,KLEV),&
 & ZUD(KLON,KLEV),         ZVD(KLON,KLEV)  
REAL(dp) ::     ZENTR(KLON),            ZHCBASE(KLON),&
 & ZMFUB(KLON),            ZMFUB1(KLON),&
 & ZDQPBL(KLON),           ZDQCV(KLON)  
REAL(dp) ::                  ZDPMEL(KLON,KLEV),ZLGLAC(KLON,KLEV)
REAL(dp) ::     ZDHPBL(KLON),       ZWUBASE(KLON)

REAL(dp) ::     ZDMFEN(KLON,KLEV),      ZDMFDE(KLON,KLEV)

INTEGER ::  ILAB(KLON,KLEV),        IDTOP(KLON),ICTOP0(KLON),ILWMIN(KLON)
INTEGER ::  IDPL(KLON) ! departure level for convection
REAL(dp) ::     ZCAPE(KLON), ZHEAT(KLON)
LOGICAL ::  LLDDRAF(KLON),  LLDDRAF2(KLON), LLDDRAF3(KLON), LLDCUM(KLON)
LOGICAL ::  LLO1, LLO2(KLON)

INTEGER :: IKB, ITOPM2, JK, JL, JN 

!*UPG change to operations
REAL(dp) ::   ZCONS2, ZDH,&
 & ZDQMIN, ZDZ, ZEPS, ZFAC, &
 & ZMFMAX, ZPBMPT, ZQUMQE, ZRO, zzz, zalv, zsfl(klon) 
!*UPG change to operations
          
REAL(dp) :: ZWN(KLON,KLEV)  ! normalized large-scale vertical velocity
REAL(dp) :: ZRHO            ! air density
REAL(dp) :: ZDXDY           ! horizontal grid surface

! scaling factor and momentum massflux
REAL(dp) :: ZMFS(KLON), ZMFUUS(KLON,KLEV), ZMFDUS(KLON,KLEV),& 
                            &  ZMFUDR(KLON,KLEV), ZMFDDR(KLON,KLEV)
REAL(dp) :: ZHOOK_HANDLE

!*UPG change to operations
!    LOCALS FOR CONSERVATION CHECK
LOGICAL :: LLCONSCHECK=.FALSE.
!!$LOGICAL :: LLCONSCHECK=.TRUE.
REAL(dp), DIMENSION(:,:), ALLOCATABLE   :: ZTENT, ZTENQ, ZTENU, ZTENV, ZSUMC
REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: ZTENC

  ! for Grant w* shallow closure
REAL(dp)  :: zmf_shal(klon)
!*UPG change to operations

!---------------------------------------------------------------------
!*UPG Change to operations call SATUR routine here

! INITITALIZATION OF UNDEFINED VARIABLES
LLDDRAF = .FALSE.
PMFLXR(:,:) = 0.0_dp
PMFLXS(:,:) = 0.0_dp


!     0.           Compute Saturation specific humidity
!                  ------------------------------------

    CALL SATUR (KIDIA , KFDIA , KLON  , KTDIA , KLEV,&
               &PAP,    PTEN  , PQSEN , 1  )


!*UPG
!---------------------------------------------------------------------

!     1.           SPECIFY CONSTANTS AND PARAMETERS
!                  --------------------------------

IF (LHOOK) CALL DR_HOOK('CUMASTRN',0,ZHOOK_HANDLE)
ZCONS2=RMFCFL/(RG*PTSPHY)

!----------------------------------------------------------------------

!*    2.           INITIALIZE VALUES AT VERTICAL GRID POINTS IN 'CUINI'
!                  ---------------------------------------------------

CALL CUININ &
 & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
 & PTEN,     PQEN,     PQSEN,    PUEN,     PVEN,&
 & PVERVEL,  PGEO,     PAPH,     PAP,&
 & ILWMIN,   ILAB,&
 & ZTENH,    ZQENH,    ZQSENH,   PGEOH,&
 & PTU,      PQU,      ZTD,      ZQD,&
 & ZUU,      ZVU,      ZUD,      ZVD,&
 & PLU     )  

!---------------------------------------------------------------------

!*    3.0          CLOUD BASE CALCULATIONS
!                  -----------------------

!*             (A) DETERMINE CLOUD BASE VALUES IN 'CUBASE'
!                  ---------------------------------------

! normalized vertical velocity for Trigger
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZRHO=PAP(JL,JK)/(RD*PTEN(JL,JK)*(1.0_dp+RETV*PQEN(JL,JK)))
    ZDXDY=RPI*RA/REAL(NSMAX)
    ZDXDY=MIN(ZDXDY,200.E3_dp)
    ZWN(JL,JK)=-PVERVEL(JL,JK)/(RG*ZRHO)*ZDXDY/25.E3_dp
  ENDDO
ENDDO

CALL CUBASEN &
 & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
 & ZTENH,    ZQENH,    PGEOH,    PAPH,&
 & PQHFL,    PAHFS,    PSSTRU,   PSSTRV,   ZWN,&
 & PTEN,     PQEN,     PGEO,&
 & PUEN,     PVEN,&
 & PTU,      PQU,      PLU,      ZUU,      ZVU,    ZWUBASE,&
 & ILAB,     LDCUM,    LDSC,     KCBOT,    KBOTSC,&
 & ICTOP0,   IDPL,     PCAPE  )   

!*             (B) DETERMINE TOTAL MOISTURE CONVERGENCE AND
!*                 DECIDE ON TYPE OF CUMULUS CONVECTION 
!*                 ONE THE BASIS OF THE DEPTH OF THE CONVECTION
!*                 DEEP IF CLOUD DEPTH > 200MB
!*                 SHALLOW IF CLOUD DEPTH <200MB
!                  -----------------------------------------

! CALCULATE COLUMN AND SUB CLOUD LAYER MOISTURE CONVERGENCE
! AND SUB CLOUD LAYER MOIST STATIC ENERGY CONVERGENCE

JK=1
DO JL=KIDIA,KFDIA
  ZDQCV(JL) =PTENQ(JL,JK)*(PAPH(JL,JK+1)-PAPH(JL,JK))
  ZDQPBL(JL)=0.0_dp
  ZDHPBL(JL)=0.0_dp
  IDTOP(JL)=0
ENDDO
DO JK=2,KLEV
  DO JL=KIDIA,KFDIA
    ZDQCV(JL)=ZDQCV(JL)+PTENQ(JL,JK)*(PAPH(JL,JK+1)-PAPH(JL,JK))
    IF(LDCUM(JL).and.JK >= KCBOT(JL)) THEN
      ZDQPBL(JL)=ZDQPBL(JL)+PTENQ(JL,JK)*(PAPH(JL,JK+1)-PAPH(JL,JK))
      ZDHPBL(JL)=ZDHPBL(JL)+(RLVTT*PTENQ(JL,JK)+RCPD*PTENT(JL,JK))&
       & *(PAPH(JL,JK+1)-PAPH(JL,JK))  
    ENDIF
  ENDDO
ENDDO

!*                 ESTIMATE CLOUD HEIGHT FOR ENTRAINMENT/DETRAINMENT
!*                 CALCULATIONS IN CUASC AND INITIAL DETERMINATION OF 
!*                 CLOUD TYPE
!*                 (MAX.POSSIBLE CLOUD HEIGHT
!*                 FOR NON-ENTRAINING PLUME, FOLLOWING A.-S.,1974)
!                  -------------------------------------------------

!DIR$ IVDEP
!OCL NOVREC      
DO JL=KIDIA,KFDIA
  IF (LDCUM(JL)) THEN
    IKB=KCBOT(JL)
    ZHCBASE(JL)=RCPD*PTU(JL,IKB)+PGEOH(JL,IKB)+RLVTT*PQU(JL,IKB)
    IF(.NOT.LDCUM(JL)) ICTOP0(JL)=-1
  ELSE
    ZHCBASE(JL)=0.0_dp
  ENDIF
ENDDO

!*                 SPECIFY INITIAL CLOUD TYPE
!*

DO JL=KIDIA,KFDIA
  IF (LDCUM(JL)) THEN
    ZPBMPT=PAPH(JL,KCBOT(JL))-PAPH(JL,ICTOP0(JL))
    IF (ZPBMPT >= RDEPTHS) THEN
      KTYPE(JL)=1
    ELSE
      KTYPE(JL)=2
    ENDIF
  ELSE
    KTYPE(JL)=0
  ENDIF
ENDDO

!*             (C) calculate initial updraught mass flux
!*                 and set lateral mixing rates
!*
!*                 for deep convection assume it is 10% of 
!*                 maximum value which is determined by the 
!*                 thickness of the layer and timestep
!*
!*                 for shallow convection calculated assuming
!*                 a balance of moist static energy in the 
!*                 sub-cloud layer (ignores present of downdraughts)
!                  ------------------------------------------

!DIR$ IVDEP
!OCL NOVREC
DO JL=KIDIA,KFDIA
  IF (LDCUM(JL)) THEN
    IKB=KCBOT(JL)
    ZMFMAX=(PAPH(JL,IKB)-PAPH(JL,IKB-1))*ZCONS2

!*UPG Grant
    zzz=MAX(0._dp,MIN(1.E3_dp, PGEO(JL,IKB)/RG ) )
    zmf_shal(JL)=.1_dp*(rg/pten(jl,klev)*zzz*&
                   &MAX(0._dp,-pahfs(jl,klev+1))/rcpd)**.3333
    zmf_shal(jl)=min(zmf_shal(jl),zmfmax)
!*UPG Grant

! deep convection

    IF (KTYPE(JL) == 1) THEN
      ZMFUB(JL)=zmfmax*0.1_dp
!      ZMFUB(JL)=zmfmax*0.1_dp*0.5_dp
      ZENTR(JL)=ENTRPEN

    ELSEIF (KTYPE(JL) == 2) THEN

! shallow convection

      ZQUMQE=PQU(JL,IKB)+PLU(JL,IKB)-ZQENH(JL,IKB)
      ZDQMIN=MAX(0.01_dp*ZQENH(JL,IKB),1.E-10_dp)
      ZDH=RCPD*(PTU(JL,IKB)-ZTENH(JL,IKB))+RLVTT*ZQUMQE
      ZDH=RG*MAX(ZDH,1.E5_dp*ZDQMIN)
      IF (ZDHPBL(JL) > 0.0_dp) THEN
        ZMFUB(JL)=ZDHPBL(JL)/ZDH
      !EPS: temporary solution
        if(ptsphy>1800._dp) then
          ZMFUB(JL)=MIN(ZMFUB(JL),RMFSHCFL*ZMFMAX)
        else
          ZMFUB(JL)=MIN(ZMFUB(JL),ZMFMAX)
        endif
      ELSE
        ZMFUB(JL)=ZMFMAX*0.1_dp
        LDCUM(JL)=.FALSE.
      ENDIF
      if(lmfscl_wstar) zmfub(jl)=zmf_shal(jl)
      ZENTR(JL)=ENTRSCV
    ENDIF

  ELSE

! no buoyancy cloud base from surface
! set cloud base mass flux and mixing rate
! to default value for safety

    ZMFUB(JL)=0.0_dp
    ZENTR(JL)=ENTRSCV
  ENDIF
ENDDO

!-----------------------------------------------------------------------

!*    4.0          DETERMINE CLOUD ASCENT FOR ENTRAINING PLUME
!                  -------------------------------------------

!*             (A) ESTIMATE CLOUD HEIGHT FOR ENTRAINMENT/DETRAINMENT
!*                 CALCULATIONS IN CUASC (MAX.POSSIBLE CLOUD HEIGHT
!*                 FOR NON-ENTRAINING PLUME, FOLLOWING A.-S.,1974)
!                  -------------------------------------------------

! CALCULATIONS NOW DONE IS SECTION 3 ABOVE SO THAT
! INITIAL CLOUD DEPTH CAN BE USED TO SPECIFY
! THE TYPE OF CONVECTION

!*             (B) DO ASCENT IN 'CUASC'IN ABSENCE OF DOWNDRAFTS
!                  --------------------------------------------

CALL CUASCN &
 & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
 & PTSPHY,&
 & ZTENH,    ZQENH,    PUEN,     PVEN,&
 & PTEN,     PQEN,     PQSEN,    PLITOT,&
 & PGEO,     PGEOH,    PAP,      PAPH,&
 & PTENQ,    PVERVEL,  zwubase,  ILWMIN,&
 & LDLAND,   LDCUM,    KTYPE,    ILAB,&
 & PTU,      PQU,      PLU,      ZUU,      ZVU,&
 & PMFU,     ZMFUB,    ZENTR,    ZLGLAC,&
 & ZMFUS,    ZMFUQ,    ZMFUL,    PLUDE,    ZDMFUP,&
 & ZDMFEN,&
 & KCBOT,    KCTOP,    ICTOP0,   IDPL,     PMFUDE_RATE,  PMFUEN_RATE, &
 & lwc,      iwc,      rform,    sform)   

!*         (C) CHECK CLOUD DEPTH AND CHANGE ENTRAINMENT RATE ACCORDINGLY
!              CALCULATE PRECIPITATION RATE (FOR DOWNDRAFT CALCULATION)
!              -----------------------------------------------------

!DIR$ IVDEP
!OCL NOVREC
DO JL=KIDIA,KFDIA
  IF (LDCUM(JL)) THEN
    ZPBMPT=PAPH(JL,KCBOT(JL))-PAPH(JL,KCTOP(JL))
    IF(KTYPE(JL) == 1.AND.ZPBMPT < RDEPTHS) KTYPE(JL)=2
    IF(KTYPE(JL) == 2.AND.ZPBMPT >= RDEPTHS) KTYPE(JL)=1
    ICTOP0(JL)=KCTOP(JL)
    IF(KTYPE(JL) == 1) ZENTR(JL)=ENTRPEN
    IF(KTYPE(JL) == 2) ZENTR(JL)=ENTRSCV
  ENDIF
  ZRFL(JL)=ZDMFUP(JL,1)
ENDDO
DO JK=2,KLEV
  DO JL=KIDIA,KFDIA
    ZRFL(JL)=ZRFL(JL)+ZDMFUP(JL,JK)
  ENDDO
ENDDO
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    PMFD(JL,JK)=0.0_dp
    ZMFDS(JL,JK)=0.0_dp
    ZMFDQ(JL,JK)=0.0_dp
    ZDMFDP(JL,JK)=0.0_dp
    ZDPMEL(JL,JK)=0.0_dp
  ENDDO
ENDDO

!-----------------------------------------------------------------------

!*    5.0          CUMULUS DOWNDRAFT CALCULATIONS
!                  ------------------------------

IF(LMFDD) THEN

!*             (A) DETERMINE LFS IN 'CUDLFS'
!                  -------------------------

  CALL CUDLFSN &
   & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
   & KCBOT,    KCTOP,    LDLAND,   LDCUM,&
   & ZTENH,    ZQENH,    PUEN,     PVEN,&
   & PTEN,     PQSEN,    PGEO,&
   & PGEOH,    PAPH,     PTU,      PQU,      PLU,&
   & ZUU,      ZVU,      ZMFUB,    ZRFL,&
   & ZTD,      ZQD,      ZUD,      ZVD,&
   & PMFD,     ZMFDS,    ZMFDQ,    ZDMFDP,&
   & IDTOP,    LLDDRAF )  

!*            (B)  DETERMINE DOWNDRAFT T,Q AND FLUXES IN 'CUDDRAF'
!                  -----------------------------------------------

  CALL CUDDRAFN &
   & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
   & LLDDRAF,&
   & ZTENH,    ZQENH,    PUEN,     PVEN,&
   & PGEO,     PGEOH,    PAPH,     ZRFL,&
   & ZTD,      ZQD,      ZUD,      ZVD,      PMFU,&
   & PMFD,     ZMFDS,    ZMFDQ,    ZDMFDP,&
   & ZDMFDE,   PMFDDE_RATE,        PMFDEN_RATE )  

ENDIF

!*                 (C)  RECALCULATE CLOUD BASE MASSFLUX FROM A
!*                 CAPE CLOSURE FOR DEEP CONVECTION (KTYPE=1)
!*                 AND BY PBL EQUILIBRUM TAKING DOWNDRAFTS INTO
!*                 ACCOUNT FOR SHALLOW CONVECTION (KTYPE=2)          
!                  --------------------------------------------

!   DEEP CONVECTION

DO JL=KIDIA,KFDIA
  ZHEAT(JL)=0.0_dp
  ZCAPE(JL)=0.0_dp
  ZMFUB1(JL)=ZMFUB(JL)
ENDDO

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    LLO1=LDCUM(JL).AND.KTYPE(JL) == 1
    IF(LLO1.AND.JK <= KCBOT(JL).AND.JK > KCTOP(JL)) THEN
      IKB=KCBOT(JL)
      ZRO=PAPH(JL,JK)/(RD*ZTENH(JL,JK))
      ZDZ=(PAPH(JL,JK)-PAPH(JL,JK-1))/(RG*ZRO)
      ZHEAT(JL)=ZHEAT(JL) +&
       & (  (PTEN(JL,JK-1)-PTEN(JL,JK) + RG*ZDZ/RCPD)/ZTENH(JL,JK)&
       & +  0.608_dp*(PQEN(JL,JK-1)-PQEN(JL,JK))  ) *&
       & (RG*(PMFU(JL,JK)+PMFD(JL,JK)))/ZRO  
      ZCAPE(JL)=ZCAPE(JL) +&
       & (RG*(PTU(JL,JK)-ZTENH(JL,JK))/ZTENH(JL,JK)&
       & +RG*0.608_dp*(PQU(JL,JK)-ZQENH(JL,JK))&
       & -RG*PLU(JL,JK) ) * ZDZ  
    ENDIF
  ENDDO
ENDDO

DO JL=KIDIA,KFDIA
  IF(LDCUM(JL).AND.KTYPE(JL) == 1) THEN
    IKB=KCBOT(JL)
    ZMFUB1(JL)=(ZCAPE(JL)*ZMFUB(JL))/(ZHEAT(JL)*RTAU)
    ZMFUB1(JL)=MAX(ZMFUB1(JL),0.001_dp)
    ZMFMAX=(PAPH(JL,IKB)-PAPH(JL,IKB-1))*ZCONS2
    ZMFUB1(JL)=MIN(ZMFUB1(JL),ZMFMAX)
  ENDIF
ENDDO

!  SHALLOW CONVECTION AND MID_LEVEL

!DIR$ IVDEP
!OCL NOVREC
DO JL=KIDIA,KFDIA
  IF ( LDCUM(JL) .AND. (KTYPE(JL) == 2.OR. KTYPE(JL) == 3) ) THEN
    IKB=KCBOT(JL)
    IF(PMFD(JL,IKB) < 0.0_dp) THEN
      ZEPS=-PMFD(JL,IKB)/MAX(ZMFUB(JL),1.E-10_dp)
    ELSE
      ZEPS=0.0_dp
    ENDIF
    ZQUMQE=PQU(JL,IKB)+PLU(JL,IKB)-&
     & ZEPS*ZQD(JL,IKB)-(1.0_dp-ZEPS)*ZQENH(JL,IKB)  
    ZDQMIN=MAX(0.01_dp*ZQENH(JL,IKB),1.E-10_dp)
! maximum permisable value of ud base mass flux 
    ZMFMAX=(PAPH(JL,IKB)-PAPH(JL,IKB-1))*ZCONS2

! shallow convection

    IF(KTYPE(JL) == 2) THEN
      ZDH=RCPD*(PTU(JL,IKB)-ZEPS*ZTD(JL,IKB)-&
       & (1.0_dp-ZEPS)*ZTENH(JL,IKB))+RLVTT*ZQUMQE  
      ZDH=RG*MAX(ZDH,1.E5_dp*ZDQMIN)
      IF(ZDHPBL(JL) > 0.0_dp) THEN
        ZMFUB1(JL)=ZDHPBL(JL)/ZDH
      ELSE
        ZMFUB1(JL)=ZMFUB(JL)
      ENDIF
   !EPS: temporary solution for EPS
      if(ptsphy>1800._dp) then
        zmfub1(jl)=min(zmfub1(jl),RMFSHCFL*zmfmax)
      else
        zmfub1(jl)=min(zmfub1(jl),zmfmax)
      endif
      if(lmfscl_wstar) zmfub1(jl)=zmf_shal(jl)
    ENDIF

! mid-level convection

    IF(KTYPE(JL) == 3)THEN
      ZMFUB1(JL)=ZMFUB(JL)*(1.0_dp+ZEPS)
      ZMFUB1(JL)=MIN(ZMFUB1(JL),ZMFMAX)
    ENDIF

  ENDIF
ENDDO

! rescale DD fluxes if deep and shallow convection

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    IF ( LLDDRAF(JL) .AND.( KTYPE(JL) == 1.OR. KTYPE(JL) == 2 ) ) THEN
      ZFAC=ZMFUB1(JL)/MAX(ZMFUB(JL),1.E-10_dp)
      PMFD(JL,JK)=PMFD(JL,JK)*ZFAC
      ZMFDS(JL,JK)=ZMFDS(JL,JK)*ZFAC
      ZMFDQ(JL,JK)=ZMFDQ(JL,JK)*ZFAC
      ZDMFDP(JL,JK)=ZDMFDP(JL,JK)*ZFAC
!  also rescale detrainment flux for ERA pp
      PMFDDE_RATE(JL,JK)=PMFDDE_RATE(JL,JK)*ZFAC
      PMFDEN_RATE(JL,JK)=PMFDEN_RATE(JL,JK)*ZFAC
    ENDIF
  ENDDO
ENDDO

DO JK=2,KLEV-1
  DO JL=KIDIA,KFDIA
    ZUU(JL,JK)=PUEN(JL,JK-1)
    ZVU(JL,JK)=PVEN(JL,JK-1)
  ENDDO
ENDDO

! reset updraught mass flux at cloud base

DO JL=KIDIA,KFDIA
  ZMFUB(JL)=ZMFUB1(JL)
ENDDO

!-----------------------------------------------------------------------

!*    6.0          DETERMINE FINAL CLOUD ASCENT FOR ENTRAINING PLUME
!*                 FOR PENETRATIVE CONVECTION (TYPE=1),
!*                 FOR SHALLOW TO MEDIUM CONVECTION (TYPE=2)
!*                 AND FOR MID-LEVEL CONVECTION (TYPE=3).
!                  -------------------------------------------------

CALL CUASCN &
 & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
 & PTSPHY,&
 & ZTENH,    ZQENH,    PUEN,     PVEN,&
 & PTEN,     PQEN,     PQSEN,    PLITOT,&
 & PGEO,     PGEOH,    PAP,      PAPH,&
 & PTENQ,    PVERVEL,  zwubase,  ILWMIN,&
 & LDLAND,   LDCUM,    KTYPE,    ILAB,&
 & PTU,      PQU,      PLU,      ZUU,      ZVU,&
 & PMFU,     ZMFUB,    ZENTR,    ZLGLAC,&
 & ZMFUS,    ZMFUQ,    ZMFUL,    PLUDE,    ZDMFUP,&
 & ZDMFEN,&
 & KCBOT,    KCTOP,    ICTOP0,   IDPL,     PMFUDE_RATE,   PMFUEN_RATE,&
 & lwc,      iwc,      rform,    sform   )   

!-----------------------------------------------------------------------

!*    6.5          IN CASE THAT EITHER DEEP OR SHALLOW IS SWITCHED OFF
!                  RESET LDCUM TO FALSE-> FLUXES SET TO ZERO IN CUFLXN
!                  ---------------------------------------------------

IF (.NOT.LMFSCV .OR. .NOT.LMFPEN) THEN
  DO JL=KIDIA,KFDIA
    LLO2(JL)=.FALSE.
    IF((.NOT.LMFSCV .AND. KTYPE(JL)==2).OR.(.NOT.LMFPEN .AND. KTYPE(JL)==1)) THEN
      LLO2(JL)=.TRUE.
      LDCUM(JL)=.FALSE.
    ENDIF
  ENDDO
ENDIF
!-----------------------------------------------------------------------

!*    7.0          DETERMINE FINAL CONVECTIVE FLUXES IN 'CUFLX'
!                  ------------------------------------------

CALL CUFLXN &
 & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
 & PTSPHY,&
 & PTEN,     PQEN,     PQSEN,    ZTENH,    ZQENH,&
 & PAPH,     PAP,      PGEOH,    LDLAND,   LDCUM,&
 & KCBOT,    KCTOP,    IDTOP,    ITOPM2,&
 & KTYPE,    LLDDRAF,&
 & PMFU,     PMFD,     ZMFUS,    ZMFDS,&
 & ZMFUQ,    ZMFDQ,    ZMFUL,    PLUDE,&
 & ZDMFUP,   ZDMFDP,   ZDPMEL,   ZLGLAC,&
 & PMFLXR,   PMFLXS,   PRAIN,    PMFDDE_RATE )  

! rescale DD fluxes if total mass flux becomes negative

ZMFS(:)=1._dp
DO JK=1,KLEV
   DO JL=KIDIA,KFDIA
      IF ( LLDDRAF(JL) .AND. JK>=IDTOP(JL)-1 ) THEN
         ZMFMAX=PMFU(JL,JK)*0.98_dp
         IF(PMFD(JL,JK)+ZMFMAX+1.E-10_dp<0._dp) THEN
            ZMFS(JL)=MIN(ZMFS(JL),-ZMFMAX/PMFD(JL,JK))
         ENDIF
      ENDIF
   ENDDO
ENDDO

DO JK=2,KLEV
  DO JL=KIDIA,KFDIA
    IF ( LLDDRAF(JL) .AND. ZMFS(JL)<1._dp ) THEN
      PMFD(JL,JK)=PMFD(JL,JK)*ZMFS(JL)
      ZMFDS(JL,JK)=ZMFDS(JL,JK)*ZMFS(JL)
      ZMFDQ(JL,JK)=ZMFDQ(JL,JK)*ZMFS(JL)
      PMFDDE_RATE(JL,JK)=PMFDDE_RATE(JL,JK)*ZMFS(JL)
      ZDMFDP(JL,JK)=ZDMFDP(JL,JK)*ZMFS(JL)
    ENDIF
  ENDDO
ENDDO

!*UPG change to operations
IF ( LLCONSCHECK ) THEN
    ALLOCATE(ZTENT(KLON,KLEV))
    ALLOCATE(ZTENQ(KLON,KLEV))
    ALLOCATE(ZTENU(KLON,KLEV))
    ALLOCATE(ZTENV(KLON,KLEV))
    DO JK=2,KLEV
    DO JL=KIDIA,KFDIA
       IF ( LDCUM(JL) ) THEN
            ZTENT(JL,JK)=PTENT(JL,JK)
            ZTENQ(JL,JK)=PTENQ(JL,JK)
            ZTENU(JL,JK)=PTENU(JL,JK)
            ZTENV(JL,JK)=PTENV(JL,JK)
       ENDIF
    ENDDO
    ENDDO
    IF ( LMFTRAC .AND. KTRAC>0 ) THEN
       ALLOCATE(ZTENC(KLON,KLEV,KTRAC))
       ALLOCATE(ZSUMC(KLON,4+KTRAC))
       DO JN=1,KTRAC
          DO JK=2,KLEV
          DO JL=KIDIA,KFDIA
             IF ( LDCUM(JL) ) THEN
                ZTENC(JL,JK,JN)=PTENC(JL,JK,JN)
             ENDIF
          ENDDO
          ENDDO
       ENDDO
    ELSE
       ALLOCATE(ZSUMC(KLON,4))
    ENDIF
ENDIF
!*UPG change to operations

!----------------------------------------------------------------------

!*    8.0          UPDATE TENDENCIES FOR T AND Q IN SUBROUTINE CUDTDQ
!                  --------------------------------------------------

CALL CUDTDQN &
 & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
 & ITOPM2,   KTYPE,    LDLAND,   LDCUM,    PTSPHY,&
 & PAPH,     PTEN,     ZLGLAC,   PLUDE,&
 & ZMFUS,    ZMFDS,    ZMFUQ,    ZMFDQ,&
 & ZMFUL,    ZDMFUP,   ZDMFDP,   ZDPMEL,&
 & PTENT,    PTENQ,    PENTH )  

!----------------------------------------------------------------------

!*    9.0          UPDATE TENDENCIES FOR U AND V IN SUBROUTINE CUDUDV
!                  --------------------------------------------------

IF(LMFDUDV) THEN

! EPS: Begin
! Intermediate Solution for stability in EPS: 
! rescale massfluxes of  shallow convection for stability in Momentum
!-------------------------------------------------------------------

  DO JL=KIDIA,KFDIA
    ZMFS(JL)=1.0_dp
    LLDDRAF2(JL)=LDCUM(JL).AND.KTYPE(JL) == 2 .AND. PTSPHY>1800._dp
  ENDDO
! IF(RMFSOLUV<=0.5_dp) THEN
    DO JK=2,KLEV
      DO JL=KIDIA,KFDIA
        IF(LLDDRAF2(JL)) THEN
          ZMFMAX=(PAPH(JL,JK)-PAPH(JL,JK-1))*ZCONS2
          IF(PMFU(JL,JK)>ZMFMAX.AND.JK>=KCTOP(JL)) &
           & ZMFS(JL)=MIN(ZMFS(JL),ZMFMAX/PMFU(JL,JK))  
        ENDIF
      ENDDO
    ENDDO
! ENDIF
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZMFUUS(JL,JK)=PMFU(JL,JK)
      ZMFDUS(JL,JK)=PMFD(JL,JK)
      IF(LLDDRAF2(JL).AND.JK>=KCTOP(JL)-1) THEN
        ZMFUUS(JL,JK)=PMFU(JL,JK)*ZMFS(JL)
        ZMFDUS(JL,JK)=PMFD(JL,JK)*ZMFS(JL)
      ENDIF
    ENDDO
  ENDDO

!-------------------------------------------------------------------
! End
! Intermediate Solution for stability in EPS: 
! For original code replace line
!  &, PUEN,     PVEN,     ZMFUUS,   ZMFDUS &
!by
!  &, PUEN,     PVEN,     PMFU,     PMFD

  CALL CUDUDV &
   & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,&
   & ITOPM2,   KTYPE,    KCBOT,    KCTOP,    LDCUM,    PTSPHY,&
   & PAPH,     PUEN,     PVEN,     ZMFUUS,   ZMFDUS,&
   & ZUU,      ZUD,      ZVU,      ZVD,&
   & PTENU,    PTENV     )  

ENDIF

!----------------------------------------------------------------------

!*   10.           IN CASE THAT EITHER DEEP OR SHALLOW IS SWITCHED OFF
!                  NEED TO SET SOME VARIABLES A POSTERIORI TO ZERO
!                  ---------------------------------------------------

IF (.NOT.LMFSCV .OR. .NOT.LMFPEN) THEN
  DO JK=2,KLEV
    DO JL=KIDIA,KFDIA
      IF(LLO2(JL).AND.JK>=KCTOP(JL)-1) THEN
        PTU(JL,JK)  =PTEN(JL,JK)
        PQU(JL,JK)  =PQEN(JL,JK)
        PLU(JL,JK)  =0.0_dp
        PENTH(JL,JK) =0.0_dp
        PMFUDE_RATE(JL,JK) =0.0_dp
        PMFUEN_RATE(JL,JK) =0.0_dp
        PMFDDE_RATE(JL,JK) =0.0_dp
        PMFDEN_RATE(JL,JK) =0.0_dp
      ENDIF
    ENDDO
  ENDDO
  DO JL=KIDIA,KFDIA
    IF(LLO2(JL)) THEN
      KCTOP(JL)=KLEV-1
      KCBOT(JL)=KLEV-1
    ENDIF
  ENDDO
ENDIF

!----------------------------------------------------------------------

!*   11.0          CHEMICAL TRACER TRANSPORT
!                  -------------------------

IF ( LMFTRAC .AND. KTRAC>0 ) THEN

 DO JL=KIDIA,KFDIA
    !IF( LDCUM(JL).AND.KTYPE(JL)/=3 ) THEN
     IF( LDCUM(JL).AND.KTYPE(JL)/=3.and.kcbot(jl)-kctop(jl)>=1 ) THEN
         LLDCUM(JL)=.TRUE.
         LLDDRAF3(JL)=LLDDRAF(JL)
     ELSE
         LLDCUM(JL)=.FALSE.
         LLDDRAF3(JL)=.FALSE.
     ENDIF
  ENDDO

! check and correct mass fluxes for CFL criterium

  ZMFS(:)=1.0_dp
! IF(RMFSOLCT<0.8_dp) THEN
    DO JK=2,KLEV
      DO JL=KIDIA,KFDIA
        IF(LLDCUM(JL).AND.JK>=KCTOP(JL)) THEN
          ZMFMAX=(PAPH(JL,JK)-PAPH(JL,JK-1))*0.8*ZCONS2
          IF(PMFU(JL,JK)>ZMFMAX) &
           & ZMFS(JL)=MIN(ZMFS(JL),ZMFMAX/PMFU(JL,JK))
        ENDIF
      ENDDO
    ENDDO
! ENDIF
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      IF(LLDCUM(JL).AND.JK>=KCTOP(JL)-1) THEN
        ZMFUUS(JL,JK)=PMFU(JL,JK)*ZMFS(JL)
        ZMFUDR(JL,JK)=PMFUDE_RATE(JL,JK)*ZMFS(JL)
      ELSE
        ZMFUUS(JL,JK)=0._dp
        ZMFUDR(JL,JK)=0._dp
      ENDIF
      IF ( LLDDRAF3(JL) .AND. JK>=IDTOP(JL)-1) THEN
        ZMFDUS(JL,JK)=PMFD(JL,JK)*ZMFS(JL)
        ZMFDDR(JL,JK)=PMFDDE_RATE(JL,JK)*ZMFS(JL)
      ELSE
        ZMFDUS(JL,JK)=0._dp
        ZMFDDR(JL,JK)=0._dp
      ENDIF
    ENDDO
  ENDDO

  CALL CUCTRACER &
   & ( KIDIA,    KFDIA,    KLON,     KTDIA,    KLEV,     KTRAC,&
   & KTYPE,    KCTOP,    KCBOT,    IDPL,     IDTOP,&
   & LLDCUM,   LLDDRAF3, PTSPHY,   PAPH,&
   & ZMFUUS,   ZMFDUS,   ZMFUDR,   ZMFDDR,&
   & PCEN,     PTENC     )  

ENDIF

!----------------------------------------------------------------------

!*   12.           PUT DETRAINMENT RATES FROM MFLX UNITS IN UNITS MFLX/M 
!                  FOR ERA40
!                  ---------------------------------------------------

! commented because MESSy / CVTRANS needs the rates in kg/(s m**2)

!!$DO JK=2,KLEV
!!$  DO JL=KIDIA,KFDIA
!!$    IF ( LDCUM(JL) ) THEN
!!$      ZRHO=(PGEOH(JL,JK)-PGEOH(JL,JK+1))/RG  ! dz
!!$      PMFUDE_RATE(JL,JK)=PMFUDE_RATE(JL,JK)/ZRHO
!!$      PMFDDE_RATE(JL,JK)=PMFDDE_RATE(JL,JK)/ZRHO
!!$    ENDIF
!!$  ENDDO
!!$ENDDO

!----------------------------------------------------------------------
!*UPG change to operations

IF ( LLCONSCHECK ) THEN

!*   13.0          CONSERVATION CHECK and CORRECTION
!                  ---------------------------------

    DO JL=KIDIA,KFDIA
      ZSUMC(JL,:)=0.
    ENDDO
    DO JK=KLEV,2,-1
    DO JL=KIDIA,KFDIA
     IF ( LDCUM(JL) .AND. JK>=KCTOP(JL)-1) THEN
       ZDZ=(PAPH(JL,JK+1)-PAPH(JL,JK))/RG
       ZSUMC(JL,1)=ZSUMC(JL,1)+(PTENQ(JL,JK)-ZTENQ(JL,JK))*ZDZ+PLUDE(JL,JK)
       
       ZALV=FOELHMCU(PTEN(JL,JK))
       ZSUMC(JL,2)=ZSUMC(JL,2)+RCPD*(PTENT(JL,JK)-ZTENT(JL,JK))*ZDZ-ZALV*PLUDE(JL,JK)
       ZSUMC(JL,3)=ZSUMC(JL,3)+(PTENU(JL,JK)-ZTENU(JL,JK))*ZDZ
       ZSUMC(JL,4)=ZSUMC(JL,4)+(PTENV(JL,JK)-ZTENV(JL,JK))*ZDZ
     ENDIF
    ENDDO
    ENDDO

    IF ( LMFTRAC .AND. KTRAC>0 ) THEN
      DO JN=1,KTRAC
         DO JK=KLEV,2,-1
         DO JL=KIDIA,KFDIA
            IF ( LDCUM(JL) .AND. JK>=KCTOP(JL)-1) THEN
               ZDZ=(PAPH(JL,JK+1)-PAPH(JL,JK))/RG
               ZSUMC(JL,4+JN)=ZSUMC(JL,4+JN)+(PTENC(JL,JK,JN)-ZTENC(JL,JK,JN))*ZDZ
            ENDIF
          ENDDO
          ENDDO
       ENDDO
     ENDIF

    DO JL=KIDIA,KFDIA
     IF ( LDCUM(JL) ) THEN
       ZALV=FOELHMCU(PTEN(JL,KLEV))
       ZSFL(JL)=PMFLXR(JL,KLEV+1)+PMFLXS(JL,KLEV+1)

       write(61,'(i4,a9,2f15.8,i4,a9,f15.8,a10,2f15.8)')jl,' CONS q: ',&
          &-zsumc(jl,1)*zalv,zsfl(jl)*zalv,ktype(jl),&
          &' CONS h: ',zsumc(jl,2),' CONS uv: ',zsumc(jl,3),zsumc(jl,4)
       if ( lmftrac .and. ktrac>0 ) then
          write(61,*)' Conserv Error Tracers 1-',ktrac,' :'
         do jn=1,ktrac
          write(61,'(i4,e12.4)')jn,zsumc(jl,4+jn)
         enddo
       endif

       IKB=KCTOP(JL)
       ZDZ=(PAPH(JL,KLEV+1)-PAPH(JL,IKB-1))/RG
!!$       if (ABS((ZSUMC(JL,1)+ZSFL(JL))/ZDZ) > 1.e-10_dp) THEN
!!$         print*, "in convect_ecmwf",&
!!$         (ZSUMC(JL,1)+ZSFL(JL))/ZDZ, zsumc(jl,1), zsfl(jl), zdz
!!$         print*, "ptenq", ptenq(jl,:)
!!$         print*, "ztenq", ztenq(jl,:)
!!$         print*, "plude", plude(jl,:)
!!$       ENDIF
       ! mz_ht_20070910+
!       ZSUMC(JL,1)=(ZSUMC(JL,1)+ZSFL(JL))/ZDZ
       ZSUMC(JL,1)=(ZSUMC(JL,1)+ZSFL(JL))
       PMFLXR(JL,KLEV+1) = PMFLXR(JL,KLEV+1) - ZSUMC(JL,1)
       ! mz_ht_20070910-
       ZSUMC(JL,2)=(ZSUMC(JL,2)-ZALV*ZSFL(JL))/(ZDZ*RCPD)
     ENDIF
    END DO
    
  !  DO JK=KLEV,2,-1
  !    DO JL=KIDIA,KFDIA
  !      IKB=KCTOP(JL)
  !      IF ( LDCUM(JL) .AND. JK >= IKB-1) THEN
  !!    PTENQ(JL,JK)=PTENQ(JL,JK)-ZSUMC(JL,1)
  !        PTENT(JL,JK)=PTENT(JL,JK)-ZSUMC(JL,2)
  !        PTENU(JL,JK)=PTENU(JL,JK)-ZSUMC(JL,3)
  !        PTENV(JL,JK)=PTENV(JL,JK)-ZSUMC(JL,4)
  !      ENDIF
  !    ENDDO
  !  ENDDO

    DEALLOCATE(ZSUMC)
    IF ( LMFTRAC .AND. KTRAC>0 ) THEN
       DEALLOCATE(ZTENC)
    ENDIF
    DEALLOCATE(ZTENV)
    DEALLOCATE(ZTENU)
    DEALLOCATE(ZTENQ)
    DEALLOCATE(ZTENT)

ENDIF
!----------------------------------------------------------------------
!*UPG Change to operations

!*    14.0         COMPUTE CONVECTIVE TENDENCIES FOR LIQUID AND SOLID
!                  CLOUD CONDENSATE, DO NOT CHANGE PRECIP UNITS IN M/S
!                  --------------------------------------------------

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    PTENL(JL,JK)=PLUDE(JL,JK)*RG/(PAPH(JL,JK+1)-PAPH(JL,JK))
    PTENI(JL,JK)=(1.0_dp-FOEALFA(PTEN(JL,JK)))*PTENL(JL,JK)
    PTENL(JL,JK)=PTENL(JL,JK)-PTENI(JL,JK)
    PMFLXR(JL,JK)=PMFLXR(JL,JK)!*1.E-3
    PMFLXS(JL,JK)=PMFLXS(JL,JK)!*1.E-3
  ENDDO
ENDDO
!!$  DO JL=KIDIA,KFDIA
!!$    PMFLXR(JL,KLEV+1)=PMFLXR(JL,KLEV+1)*1.E-3
!!$    PMFLXS(JL,KLEV+1)=PMFLXS(JL,KLEV+1)*1.E-3
!!$  ENDDO
!*UPG Change to operations


IF (LHOOK) CALL DR_HOOK('CUMASTRN',1,ZHOOK_HANDLE)
END SUBROUTINE ECMWF_CUMASTR

!=========================================================================

SUBROUTINE zhang_cumastr(plon    ,plond   ,plev    ,plevp   ,pcnst   , &  
! kproma, kbdim, nlev, nlevp1, 2 (only two species: water and ice)
! other tracers transport in cvtrans
                         lat     ,nstep   ,tdt     ,pmid    ,pint    , &
                         pdel    ,rpdel   ,zm      ,zi      ,tpert   , &
                         qpert   ,phis    ,ts      ,pblh    ,          & 
                         t       ,q1      ,q       ,          &
!                         cmfmc   , &
!                         zdu     ,cmfdqr  ,conicw  ,precc   ,&
                         cnt     ,cnb     ,&
                         dlf     ,&
                         pflx    ,sflx,                                &
                         omga    ,wpert3d ,                            &
!+mgl
!                         zmu     , chembgt ,                           &
!-mgl
!                         zmug    ,zmdg     ,zdug      ,zeug    ,       &  
!                         zedg    ,zdpg     ,dsubcld   ,zjtg    ,zjbg , &
!                         ideep   &
                          kdeep   &
!                         ,lengath&
                         )

!Convective adjustments using CCM3 parameterizations (Zhang + Hack).

  USE messy_main_constants_mem,    ONLY: dp, g, cpres => cp_air, Tmelt
  USE messy_convect_ZHANG
  USE messy_convect_zhang_param,   ONLY: aqsatd
 ! USE messy_convect_mem,           ONLY: massfu, massfd, u_entr,        &
 !                                        u_detr, d_entr, udetr_h, cape, &
 !                                        cv_lwc, cv_iwc, cv_rform, cv_sform
  ! op_mm_20140327+
  ! use 2d fields instead of 3d fields 
  USE messy_convect_mem,           ONLY: massfu, massfd, u_entr_2d,        &
                                         u_detr_2d, d_entr_2d, udetr_h_2d, &
                                         cape,cv_lwc_2d, cv_iwc_2d,        &
                                         cv_rform_2d, cv_sform_2d
  ! op_mm_20140327-

  implicit none
!----------------------------------------------------------------------
!----------------------------------------------------------------------

!Input arguments
  integer, INTENT(IN)  ::         &
              plev,               &   ! number of levels, nlev
              plevp,              &   ! number of levels + 1, nlevp
              plon,               &   ! "longitude index", kproma
              plond,              &   ! number of longitudes, kbdim
              pcnst,              &   ! number of tracers, ntrac
              lat,                &   ! "latitude index", jrow 
              nstep                   ! current time step index
  real(dp), INTENT(IN) :: tdt,   &    ! 2 delta-t (seconds)
              pmid(plond,plev),   &    ! pressure at midpoints
              pint(plond,plevp),  &    ! pressure at interfaces
              pdel(plond,plev),   &    ! delta-p across layers
              rpdel(plond,plev),  &    ! 1./pdel
              zm(plond,plev),     &    ! height above sfc at layer midpoints
              zi(plond,plevp)          ! interface height above sfc
  real(dp) :: tpert(plond),       &    ! PBL perturbation theta
              qpert(plond),       &    ! PBL perturbation specific humidity 
              phis(plond),        &    ! surface geopotential
              ts(plond),          &    ! surface temperature
              pblh(plond)              ! PBL height (provided by PBL routine)
  real(dp) :: omga(plond,plev)         ! vertical p velocity
  real(dp) :: wpert3d(plond,plev)      ! perturbation vertical velocities
!+mgl
  logical  :: chembgt                  ! t => perform budget computations
!-mgl

  !Input/output arguments
  real(dp), INTENT(INOUT) ::      &
              t(plond,plev),      &    ! temperature (t bar)
              q1(plond,plev),     &    ! specific humidity (sh bar)
              q(plond,plev,pcnst)      ! passive tracers

!Output arguments
  real(dp) ::                     &
              cmfmc(plond,plev),  &    ! moist convection cloud mass flux
              zdu(plond,plev),    &    ! du2 from conv_ccm, 1/s
              cmfdqr(plond,plev), &    ! dq/dt due to convective rainout 
              conicw(plond,plev), &    ! in cloud water mixing ratio (convective)
              precc(plond),       &    ! convective precipitation rate
              cnt(plond),         &    ! top level of convective activity   
              cnb(plond),         &    ! bottom level of convective activity
              dlf(plond,plev),    &    ! detraining cloud water from convection
              ! conv rain flux thru interface of that lev
              pflx(plond,plevp),  &        
              ! conv rain flux thru interface of that lev
              sflx(plond,plevp)    


!    gathered arrays of driving data for zhang convective transport
  real(dp) ::                             & 
              zmug(plond,plev),           &
              zmdg(plond,plev),           &
              zdug(plond,plev),           &
              zeug(plond,plev),           &
              zedg(plond,plev),           &
              zdpg(plond,plev),           &
              dsubcld(plond)
  integer  :: zjtg(plond),                &
              zjbg(plond),                & 
              ideep(plond),               &
              lengath
  integer ::  kdeep(plond)

!Local variables:

  integer  ::  i, k
    
  real(dp) ::                     &    
              zmdt(plond,plev),   &    ! zhang convective temperature tendency
              zmdq(plond,plev),   &    ! zhang convective moisture tendency
              cme(plond,plev),    &    ! zhang condensation - evaporation
              zmu(plond,plev),    &    ! mu2 from conv_ccm, kg/m2/s
              zmd(plond,plev),    &    ! md2 from conv_ccm, kg/m2/s
              zeu(plond,plev),    &    ! eu2 from conv_ccm, 1/s
              zed(plond,plev),    &    ! ed2 from conv_ccm, 1/s
              fracis(plond,plev,pcnst) ! fraction of constituent in interstitial gas

  real(dp) ::                     &
              precc2(plond),      &    ! convective-scale preciptn rate
              cnt2(plond),        &    ! top level of convective activity
              cnb2(plond),        &    ! bottom level of convective activity
              hketa(plond,plev),  &    ! eta from cmfmca
              hkbeta(plond,plev), &    ! beta from cmfmca
              zmu2(plon,plev),    &    ! convective cloud mass flux from cmfmca
!              udetr_h(plon,plev), &    ! convective updr. detrainemt from cmfmca
              cmfdt(plond,plev),  &    ! dt/dt due to moist convection
              cmfdq(plond,plev),  &    ! dq/dt due to moist convection
              qc2(plond,plev),    &    ! dq/dt due to rainout terms
              cmfdqr2(plond,plev),&    ! dq/dt due to moist convective rainout 
              conicw2(plond,plev),&    ! in cloud water mixing ratio (convective)
              cmfmc2(plond,plev), &    ! moist convection cloud mass flux
              cmfsl2(plond,plev), &    ! moist convection lw stat energy flux
              cmflq2(plond,plev), &    ! moist convection total water flux
              tpert2(plond),      &    ! perturbation T
              qpert2(plond)            ! perturbation q
  REAL(dp) :: rform2(plond,plev), &
              sform2(plond,plev), &
              lwc2(plond,plev),   & 
              iwc2(plond,plev)
! conv rain flux thru interface of that lev for Hack
  REAL(dp) :: rflx2(plond,plevp)   
! conv snow flux thru interface of that lev for Hack
  REAL(dp) :: sflx2(plond,plevp)       

  real(dp) :: mup(plond,plevp)
  real(dp) :: mdn(plond,plevp)
  real(dp) :: wpls(plond,plev)
  real(dp) :: wprime3d(plond,plev)     ! perturbation vertical velocities
  real(dp) :: capblfc(plond)           ! negative cape below level of free conv
  real(dp) :: pblt(plond)
  integer  :: ireset(plond)

  real(dp) :: told(plond, plev)

!+mgl to convert the mass fluxes properly
  logical  :: progCwat

!++mgl for adding the evaporation of convective precipitation
  real(dp) :: rain(plond,plev)
  real(dp) :: esat(plond,plev)
  real(dp) :: qsat(plond,plev), dummy(plond,plev)
  real(dp) :: rpdeli(plond,plev)
  real(dp) :: rl
  real(dp) :: ke
  real(dp) :: evap
  real(dp) :: tmpevp
  real(dp) :: convcnd(plond,plev)
  real(dp) :: convevp(plond,plev)
  real(dp) :: totevp(plond)
  real(dp) :: totcnd(plond)
  real(dp) :: pasthru
  real(dp) :: rh

!----------------------------------------------------------------------
  progcwat=.true.
    
!    Zhang convection.
  ideep(:)   = 0
  zmdt(:,:)  = 0._dp
  zmdq(:,:)  = 0._dp
  cmfmc(:,:) = 0._dp
  cme(:,:)   = 0._dp
  zmu(:,:)   = 0._dp
  zmd(:,:)   = 0._dp
  zeu(:,:)   = 0._dp
  zed(:,:)   = 0._dp
  zdu(:,:)   = 0._dp
  cmfdqr(:,:)= 0._dp
  fracis(:,:,:)= 1._dp
  precc(:)   = 0._dp
  SFLX(:,:)  = 0._dp
!  call qneg3( 'APHYS_bz', lat, q )

  told(:,:) = t(:,:)

  if ( altConv ) then
!       add the large scale vert vel (in m/s) to the
!       pert velocities from flow over orography and pbl processes
    do k = 2,plev
      do i = 1,plon      
        wpls(i,k) = - omga(i,k)*(zi(i,k)-zi(i,k+1))/pdel(i,k)
        wprime3d(i,k) = wpert3d(i,k) + wpls(i,k)
      end do
    end do
    call conv_ccm_pjr(plev,plevp ,plon    ,plond            , & 
                t       ,q1      ,precc   ,cnt     ,cnb     , &
                pblh    ,zm      ,phis    ,zi      ,zmdq    , &
                zmdt    ,pmid    ,pint    ,pdel    ,ts      , &
                .5*tdt  ,cmfmc   ,cme     ,lat              , &
                tpert   ,dlf     ,pflx    ,zdu              , &
                cmfdqr  ,mup     ,mdn     ,wprime3d,capblfc , &
                pblt    ,ireset  ,                            &
                zmug    ,zmdg    ,zdug    ,zeug    ,zedg    , &
                zdpg    ,dsubcld ,zjtg    ,zjbg    ,ideep   , &
                lengath ,conicw  ,                            &
!!                cape(:,lat), cv_lwc(:,:,lat), cv_iwc(:,:,lat),&
!!                cv_rform(:,:,lat), cv_sform(:,:,lat))
! op_mm_20140327+
! use 2d fields instead of 3d fields                 
                cape(:,lat), cv_lwc_2d(:,:), cv_iwc_2d(:,:),&
                cv_rform_2d(:,:), cv_sform_2d(:,:) )
! op_m_20140327-

  else
    call conv_ccm(plev    ,plevp   ,plon    ,plond            , &
                  t       ,q1      ,precc   ,cnt     ,cnb     , &
                  pblh    ,zm      ,phis    ,zi      ,zmdq    , &
                  zmdt    ,pmid    ,pint    ,pdel    ,ts      , &
                  .5*tdt  ,cmfmc   ,cme                       , &
                  tpert   ,dlf     ,pflx    ,zdu              , &
                  cmfdqr  ,                                     &
                  zmug    ,zmdg    ,zdug    ,zeug    ,zedg    , &
                  zdpg    ,dsubcld ,zjtg    ,zjbg    ,ideep   , &
                  lengath ,conicw  ,                            &
!                  cape(:,lat), cv_lwc(:,:,lat), cv_iwc(:,:,lat),&
!                 cv_rform(:,:,lat), cv_sform(:,:,lat) )
! op_mm_20140327+
! use 2d fields instead of 3d fields 
                  cape(:,lat), cv_lwc_2d(:,:), cv_iwc_2d(:,:),  &
                  cv_rform_2d(:,:), cv_sform_2d(:,:) )
! op_mm_20140327-
  end if

 
  if ( progCwat ) then
!       transport cloud water only
    call convtran                              &
         (q       ,2       ,zmug    ,zmdg    , & 
          zdug    ,zeug    ,zedg    ,zdpg    , &
          zjtg    ,zjbg    ,ideep   ,          &
          1       ,lengath ,                   & 
          .5*tdt  ,fracis  ,plond   ,plev)
  end if

!   determine mask whether deep convection occurs or not

  kdeep = 0
  do i=1,lengath
    if ( maxval(zmug(i,:)) .gt. 1.e-8_dp ) then
      kdeep(ideep(i)) = 1
    else
      kdeep(ideep(i)) = 0
    endif
  enddo


! mz_ht_20051201 
!############ insert subcloud evaporation according to Wilcox
  IF (EVAP_SUB) THEN

!++mgl - added Eric Wilcox's sub-cloud evaporation evaporation scheme
!     I've done this ahead of Hack, he did it in CCM after Hack.
!!$      if (pflx(1,plevp).gt.0.) then
!!$         write(*,*) 'CHECK WHY PFLX AT SURFACE IS NOT ALWAYS RAIN RATE!'
!!$         write(*,*) 'Precip at surface before correc: ', pflx(1,plevp)/1000._dp,precc(1)
!!$      end if
      rl = 2.5104E6_dp
!     cpres = 1004.64_dp used from messy_main_constants_mem
      ke = 1.0e-5_dp

! for the first round, skip the cloud condition and only use RH condition
!      do i = 1, plon
!        do k = 1, plev-1
!          rpdeli(i,k) = 1./(pmid(i,k+1) - pmid(i,k))
!        end do
!      end do
!C run cloud fraction scheme to find where the large-scale
!C clouds are.  only allow evaporation below the upper-trop
!C extended cloud deck
!      call cldfrc(pmid    ,rpdeli  ,t       ,q       ,omga    ,
!     $            cnt     ,cnb     ,cld     ,clc     ,pdel    ,
!     $            cmfmc   ,oro     ,snowh   )
      call aqsatd(t       ,pmid    ,esat    ,qsat    ,dummy, &
                  plond   ,plon    ,plev    ,1       ,plev    )

! find the precip evaporated in convective downdrafts (see
! comment below)
      do i = 1,plon
        convcnd(i,1) = 0._dp
        convevp(i,1) = 0._dp
        totevp(i) = 0._dp
        totcnd(i) = 0._dp
        do k = 2,plev
!          tmpevp = max(pflx(i,k),0.) - max(pflx(i,k-1),0.)
          tmpevp = max(pflx(i,k+1),0._dp) - max(pflx(i,k),0._dp)
          if (tmpevp < 0._dp) then
            convevp(i,k) = -tmpevp
            convcnd(i,k) = 0._dp
            totevp(i) = totevp(i) - tmpevp
          else
            if (tmpevp > 0._dp) then
              convcnd(i,k) = tmpevp
              convevp(i,k) = 0._dp
              totcnd(i) = totcnd(i) + tmpevp
            else
              convevp(i,k) = 0._dp
              convcnd(i,k) = 0._dp
            end if
          end if
        end do
      end do
      do i = 1, plon
!!$! rain in units of kg/m^2/s
!!$!        if (pflx(i,18) .gt. 0.) then
!!$        if (pflx(i,plev) .gt. 0.) then
!!$          rain(i) = max(pflx(i,1),0.)
!!$          do k = 2, plev
!!$! if evaporation in convective downdrafts, remove that first
!!$! because the only rain I want to make available for evaporation
!!$! is that rain which will not be evaporated in a convective downdraft
!!$            rain(i) = rain(i) - convevp(i,k)
!!$            totevp(i) = totevp(i) - convevp(i,k)
!!$            pasthru = max((totevp(i)-totcnd(i)),0.)
!!$            rain(i) = rain(i) - pasthru
!!$            rh = q1(i,k)/qsat(i,k)
!!$! physicist's way
!!$!              rh = es(i,k)/esat(i,k)
!!$!              if (cld(i,k) .lt. 0.5 .and. rh .lt. 1
!!$!++mgl first round, only use an RH threshold
!!$            if ( (rh .lt. 0.9) .and. (rain(i) .gt. 0) ) then
!!$!C compute the below-cloud evaporation. evap in units of kg/kg/s
!!$              evap = max(ke*(1.0-rh)*sqrt(rain(i)),0.)
!!$              evap = min(evap, (qsat(i,k)-q1(i,k))/tdt)
!!$              evap = min((rain(i)*g/pdel(i,k)), evap)
!!$! adjust the moisture and temperature fields accordingly
!!$              q1(i,k)  = q1(i,k) + evap*tdt
!!$              t(i,k)  = t(i,k) - evap*tdt*rl/cpres
!!$              zmdq(i,k) = zmdq(i,k) + evap
!!$              zmdt(i,k) = zmdt(i,k) - evap*rl/cpres
!!$            else
!!$              evap = 0.
!!$            end if
!!$            totcnd(i) = totcnd(i) - convcnd(i,k)
!!$! if condensation in convective updrafts, add it here
!!$            rain(i) = rain(i) + convcnd(i,k) -        &
!!$                       (evap*pdel(i,k)/g) + pasthru
!!$          end do
!!$        else
!!$          rain(i) = 0
!!$        end if
!!$      end do
!!$! convert units of rain from kg/m^2/s to m/s
!!$      do i = 1, plon
!!$        rain(i) = rain(i)/1000.
!!$!         precc(i) = rain(i) + precc2(i)
!!$!  the addition of precc2 is only needed if done after Hack
!!$        if (ABS(precc(i) - rain(i))/rain(i) > 1.e-3_dp) print*, "WARNING, Significant change", Precc(i), rain(i)
!!$        precc(i) = rain(i) 
!!$      end do
!!$!
!!$! end added block

! rain in units of kg/m^2/s
!        if (pflx(i,18) .gt. 0.) then
        if (pflx(i,plev) .gt. 0._dp) then
          rain(i,1) = max(pflx(i,1),0._dp)
          do k = 2, plev
! if evaporation in convective downdrafts, remove that first
! because the only rain I want to make available for evaporation
! is that rain which will not be evaporated in a convective downdraft
            rain(i,k) = rain(i,k-1) - convevp(i,k)
!!$            if ( (rain(i,k) == 0._dp) .and.(rain(i,k-1)>0._dp) )&
!!$              print*, "1. Fehler: ", rain(i,k), rain(i,k-1), convevp(i,k),i,k 
            totevp(i) = totevp(i) - convevp(i,k)
            pasthru = max((totevp(i)-totcnd(i)),0._dp)
            rain(i,k) = rain(i,k) - pasthru
!!$            if ( (rain(i,k) == 0._dp) .and.((rain(i,k)+pasthru)>0._dp)) &
!!$              print*, "2. Fehler: ", rain(i,k), rain(i,k-1), convevp(i,k), pasthru,i,k 
            rh = q1(i,k)/qsat(i,k)
! physicist's way
!              rh = es(i,k)/esat(i,k)
!              if (cld(i,k) .lt. 0.5 .and. rh .lt. 1
!++mgl first round, only use an RH threshold
            if ( (rh .lt. 0.9_dp) .and. (rain(i,k) .gt. 0._dp) ) then
!C compute the below-cloud evaporation. evap in units of kg/kg/s
              evap = max(ke*(1.0_dp-rh)*sqrt(rain(i,k)),0._dp)
              evap = min(evap, (qsat(i,k)-q1(i,k))/tdt)
              evap = min((rain(i,k)*g/pdel(i,k)), evap)
! adjust the moisture and temperature fields accordingly
              q1(i,k)  = q1(i,k) + evap*tdt
              t(i,k)  = t(i,k) - evap*tdt*rl/cpres
              zmdq(i,k) = zmdq(i,k) + evap
              zmdt(i,k) = zmdt(i,k) - evap*rl/cpres
            else
              evap = 0._dp
            end if
            totcnd(i) = totcnd(i) - convcnd(i,k)
! if condensation in convective updrafts, add it here
            rain(i,k) = rain(i,k) + convcnd(i,k) -        &
                       (evap*pdel(i,k)/g) + pasthru
!!$            if ( (rain(i,k) == 0._dp) .and. (pflx(i,k+1).ne.0._dp).and.(pflx(i,plev) >0._dp) ) then
!!$              print*, 'ERROR: ',"rain",rain(i,:)
!!$              print*, '       ',"pflx", pflx(i,:)
!!$              print*, '       ',"convcnd", convcnd(i,k)
!!$              print*, '       ', (evap*pdel(i,k)/g), pasthru,i,k
!!$            endif
          end do
        else
          rain(i,plev) = 0._dp
        end if
      end do
! convert units of rain from kg/m^2/s to m/s
      do i = 1, plon
!         precc(i) = rain(i) + precc2(i)
!  the addition of precc2 is only needed if done after Hack
    
!!$        if ( (rain(i,plev)>0._dp) .and. (ABS(precc(i) - rain(i,plev)/1000._dp)/rain(i,plev)/1000._dp > 1.e-3_dp) )&
!!$          print*, "WARNING, Significant change", Precc(i), rain(i,plev)
!!$        if (rain(i,plev) == 0._dp .and. precc(i) > 0._dp) print*, 'TEST: ',precc(i), rain(i,plev), pflx(i,plev)
!!$        precc(i) = rain(i,plev) /1000._dp
      end do

      do i=1,plon
        do k=1,plev
          if (pflx(i,plev) .gt. 0._dp) pflx(i,k+1) = rain(i,k)
        enddo
      enddo

!!$      print*, "In new Block"
    END IF
!############ end of subcloud evaporation according to Wilcox
! mz_ht_20051201 


!    Hack convection.

  cmfdt(:,:)   = 0._dp
  cmfdq(:,:)   = 0._dp
  cmfmc2(:,:)  = 0._dp
  cmfdqr2(:,:) = 0._dp
  cmfsl2(:,:)  = 0._dp 
  cmflq2(:,:)  = 0._dp
  rflx2(:,:)   = 0._dp
  sflx2(:,:)   = 0._dp
  lwc2(:,:)    = 0._dp
  iwc2(:,:)    = 0._dp
  rform2(:,:)  = 0._dp
  sform2(:,:)  = 0._dp
  
  do i=1,plon
    precc2(i) = 0._dp
    tpert2(i) = 0.
    qpert2(i) = qpert(i)
  end do

  call cmfmca(plev    ,plevp   ,plon    ,plond   ,pcnst   ,   &
              tdt     ,pmid    ,pdel    ,                     &
              rpdel   ,zm      ,tpert2  ,qpert2  ,phis    ,   &
              pblh    ,t       ,q1      ,cmfdt   ,cmfdq   ,   &
              cmfmc2  ,cmfdqr2 ,cmfsl2  ,cmflq2  ,precc2  ,   &
              qc2     ,cnt2    ,cnb2    ,                     &
              conicw2 ,                                       &
              q, pcnst, hketa, hkbeta   ,                     &
              rflx2, sflx2, lwc2, iwc2, rform2, sform2   )


!    Merge Zhang and Hack diagnostics.

  do i=1,plon
    precc(i) = precc(i) + precc2(i)
    if (cnt2(i) .lt. cnt(i)) cnt(i) = cnt2(i)
    if (cnb2(i) .gt. cnb(i)) cnb(i) = cnb2(i)
  end do

  do k=1,plev
    do i=1,plon
      cmfmc(i,k)  = cmfmc(i,k)  + cmfmc2(i,k)
      cmfdqr(i,k) = cmfdqr(i,k) + cmfdqr2(i,k)
      conicw(i,k)  = conicw(i,k) + conicw2(i,k)
    end do
  end do

  do k=2,plevp
    do i=1,plon
      IF (T(i,k-1) < Tmelt) THEN
        SFLX(i,k)   = PFLX(i,k)
        PFLX(i,k)   = 0._dp
      ENDIF
      pflx(i,k)   = pflx(i,k) + rflx2(i,k) 
      sflx(i,k)   = sflx(i,k) + sflx2(i,k) 
    enddo
  enddo
  do k=1,plev
    do i=1,plon
   !!   cv_lwc(i,k,lat)   = cv_lwc(i,k,lat)   + lwc2(i,k)
   !!   cv_iwc(i,k,lat)   = cv_iwc(i,k,lat)   + iwc2(i,k)
   !!   cv_rform(i,k,lat) = cv_rform(i,k,lat) + rform2(i,k)
   !!   cv_sform(i,k,lat) = cv_sform(i,k,lat) + sform2(i,k)
! op_mm_20140131
! use 2d fields instead of 3d fields        
       cv_lwc_2d(i,k)   = cv_lwc_2d(i,k)   + lwc2(i,k)
       cv_iwc_2d(i,k)   = cv_iwc_2d(i,k)   + iwc2(i,k)
       cv_rform_2d(i,k) = cv_rform_2d(i,k) + rform2(i,k)
       cv_sform_2d(i,k) = cv_sform_2d(i,k) + sform2(i,k)
    end do
  end do


  do k = 1, plev
    do i = 1, lengath
!          Convert mass flux from mb/s to kg/m^2/s
      zmu(ideep(i),k) = zmug(i,k) * 100. / g
      zmd(ideep(i),k) = zmdg(i,k) * 100. / g
      zeu(ideep(i),k) = zeug(i,k)
      zdu(ideep(i),k) = zdug(i,k)
      zed(ideep(i),k) = zedg(i,k)
    end do
  end do

 ! Hack cloud updraft massflux

!  udetr_h(:,:,lat) = 0._dp
  udetr_h_2d(:,:) = 0._dp   ! op_mm_20140327 use 2d fields 
  zmu2(:,:) = 0._dp
  do k=plev-1,1,-1
    do i=1,plon
      zmu2(i,k) = hketa(i,k) + hkbeta(i,k+1) * hketa(i,k+1)
!      udetr_h(i,k,lat) = (1._dp - hkbeta(i,k)) * hketa(i,k)
      udetr_h_2d(i,k) = (1._dp - hkbeta(i,k)) * hketa(i,k)  ! op_mm_20140327 use 2d field
    end do
  end do
  do i=1,plon
    zmu2(i,plev) = hketa(i,plev)
!    udetr_h(i,plev,lat) = (1._dp - hkbeta(i,plev)) * hketa(i,plev)
    udetr_h_2d(i,plev) = (1._dp - hkbeta(i,plev)) * hketa(i,plev) ! op_mm_20140327 use 2d field
  enddo
  
  massfu(1:plon,:,lat) = zmu(1:plon,:) + zmu2(1:plon,:)
  massfd(1:plon,:,lat) = zmd(1:plon,:)
  do k=1,plev      
    do i=1,plon
!!      u_entr(i,k,lat) = zeu(i,k) * pdel(i,k) / g
!!      u_detr(i,k,lat) = zdu(i,k) * pdel(i,k) / g
!!      d_entr(i,k,lat) = zed(i,k) * pdel(i,k) / g
! op_mm_20140131
! use 2d fields instead of 3d fields 
       u_entr_2d(i,k) = zeu(i,k) * pdel(i,k) / g
       u_detr_2d(i,k) = zdu(i,k) * pdel(i,k) / g
       d_entr_2d(i,k) = zed(i,k) * pdel(i,k) / g
    enddo
  enddo
  !!u_detr(1:plon,:,lat) =  u_detr(1:plon,:,lat) + udetr_h(1:plon,:,lat)
  u_detr_2d(1:plon,:) =  u_detr_2d(1:plon,:) + udetr_h_2d(1:plon,:)   ! op_mm_20140327 use 2d field

  return
end SUBROUTINE zhang_cumastr

!=========================================================================

SUBROUTINE bechtold_cumastr( &

!-------------------------------------------------------------------------------
!   ############################################################################
!    SUBROUTINE CONVECTION( 
                           KLON, KLEV, KIDIA, KFDIA, KBDIA, KTDIA,             &
                           PDTCONV,                                            &
                           PPABS, PPAH, PZZ, PDXDY, PHSFLX,                    &
                           PT, PRV, PRC, PRI, PU, PV, PW,                      &
                           KCOUNT, PTTEN, PRVTEN, PRCTEN, PRITEN,              &
!                          PPRTEN, PPRSTEN,                                    &
                           PUMF, PDMF,                                         &
                           PPRLFLX, PPRSFLX,                                   &
                           PCAPE,                                              &
                           KCLTOP, KCLBAS, LDEEP,                              &
!                          PURV, PURCI,                                        &
                           PUTEN, PVTEN,                                       &
                           KCH1, PCH1, PCH1TEN,                                &
                           ! for ERA40 / MESSy
                           PUDR, PDDR, PUER, PDER, WAT_DIAG,                   &
                           PURLIQ, PURICE, PURR_P, PURS_P )
!   ############################################################################
!
!!**** Interface routine to the fast Meso-NH convection code developed for ECMWF/ARPEGE
!!     having a structure typical for operational routines
!!     
!!     Transformations necessary to call deep+ shallow code
!!     - skip input vertical arrays/levels : bottom=1, top=KLEV
!!     - transform specific humidities in mixing ratio
!!
!!
!!    PURPOSE
!!    -------
!!      The routine interfaces the MNH convection code as developed for operational
!!      forecast models like ECMWF/ARPEGE or HIRLAM with the typical Meso-NH array structure
!!      Calls the deep and/or shallow convection routine
!!
!!
!!**  METHOD
!!    ------
!!     Returns one tendency for shallow+deep convection but each part can
!!     be activated/desactivated separately
!!     For deep convection one can enable up to 3 additional ensemble members
!!     - this substantially improves the smoothness of the scheme and reduces
!!       allows for runs with different cloud radii (entrainment rates) and
!!       reduces the arbitrariness inherent to convective trigger condition
!!      
!!     
!!
!!    EXTERNAL
!!    --------
!!    CONVECT_DEEP
!!    CONVECT_SHALLOW
!!    SU_CONVPAR, SU_CONVPAR1
!!    SUCST:   ECMWF/ARPEGE routine
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    11/12/98
!!      modified    20/03/2002 by P. Marquet : transformed for ARPEGE/Climat
!!                             (tsmbkind.h, REAL_B, INTEGER_M, _JPRB, "& &",
!!                              _ZERO_, _ONE_, _HALF_)
!!      modified    11/04/O2 allow for ensemble of deep updrafts/downdrafts
!!
!!    REFERENCE
!!    ---------
!!    Bechtold et al., 2001, Quart. J. Roy. Meteor. Soc., Vol 127, pp 869-886: 
!!           A mass flux convection scheme for regional and global models.
!!
!-------------------------------------------------------------------------------

!#include "tsmbkind.h"                                                                   
  USE MESSY_CONVECT_MEM,            ONLY: ODEEP, OSHAL, ODOWN, OREFRESH_ALL, &
                                          OSETTADJ, OUVTRANS, OCHTRANS,      &
                                          KENSM, KICE, PTADJD, PTADJS
  USE MESSY_MAIN_CONSTANTS_MEM,     ONLY: dp
  USE MESSY_CONVECT_BECHTOLD,       ONLY: convect_deep, convect_shallow
  
  USE messy_convect_bechtold_param, ONLY: su_convpar, su_convpar1, su_convpar_shal, &
                                          initialize_satw, RG


!
!*       0.    DECLARATIONS
!              ------------
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
                                                      
INTEGER,                    INTENT(IN)   :: KLON   ! horizontal dimension
INTEGER,                    INTENT(IN)   :: KLEV   ! vertical dimension
INTEGER,                    INTENT(IN)   :: KIDIA  ! value of the first point in x
INTEGER,                    INTENT(IN)   :: KFDIA  ! value of the last point in x
INTEGER,                    INTENT(IN)   :: KBDIA  ! vertical  computations start at
                                                   ! KBDIA that is at least 1
INTEGER,                    INTENT(IN)   :: KTDIA  ! vertical computations can be
                                                   ! limited to KLEV + 1 - KTDIA
                                                   ! default=1
REAL(dp),                   INTENT(IN)   :: PDTCONV! Interval of time between two
                                                   ! calls of the deep convection
                                                   ! scheme
!LOGICAL,                     INTENT(IN)   :: ODEEP  ! switch for deep convection
!LOGICAL,                     INTENT(IN)   :: OSHAL  ! switch for shallow convection
!LOGICAL,                     INTENT(IN)   :: OREFRESH_ALL ! refresh or not all 
                                                           ! tendencies  at every call
!LOGICAL,                     INTENT(IN)   :: ODOWN  ! take or not convective
                                                     ! downdrafts into account
!INTEGER,                     INTENT(IN)   :: KICE   ! flag for ice ( 1 = yes, 
                                                     !                0 = no ice )
!INTEGER,                     INTENT(IN)   :: KENSM  ! number of additional deep convection calls
                                                     ! for ensemble (presently limited to 3)
                                                     ! KENSM=0 corresponds to base run with
                                                     ! 1 deep and 1 shallow call
!LOGICAL,                      INTENT(IN)   :: OSETTADJ ! logical to set convective
                                                        ! adjustment time by user
!REAL(dp),                     INTENT(IN)   :: PTADJD ! user defined deep adjustment time (s)
!REAL(dp),                     INTENT(IN)   :: PTADJS ! user defined shal. adjustment time (s)
!
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN)   :: PT     ! grid scale T at time t  (K)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN)   :: PRV    ! grid scale water vapor  (kg/kg)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN)   :: PRC    ! grid scale r_c (kg/kg)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN)   :: PRI    ! grid scale r_i (kg/kg)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN)   :: PU     ! grid scale horiz. wind u (m/s) 
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN)   :: PV     ! grid scale horiz. wind v (m/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN)   :: PW     ! grid scale vertical velocity (m/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN)   :: PPABS  ! grid scale pressure (Pa)
REAL(dp), DIMENSION(KLON,KLEV+1),INTENT(IN)  :: PPAH   ! half level pressure (Pa)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN)   :: PZZ    ! geopotential (m2/s2) 
REAL(dp), DIMENSION(KLON),      INTENT(IN)   :: PDXDY  ! grid area (m2)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN)   :: PHSFLX ! turbulent sensible heat flux (W/m^2)

INTEGER, DIMENSION(KLON), INTENT(INOUT)  :: KCOUNT ! convective counter(recompute
                                                   ! tendency or keep it
!   
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT):: PTTEN  ! convective temperat. tendency (K/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT):: PRVTEN ! convective r_v tendency (1/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT):: PRCTEN ! convective r_c tendency (1/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT):: PRITEN ! convective r_i tendency (1/s)

!!$REAL(dp), DIMENSION(KLON),      INTENT(INOUT):: PPRTEN ! total surf precipitation tendency (m/s)
!!$REAL(dp), DIMENSION(KLON),      INTENT(INOUT):: PPRSTEN! solid surf precipitation tendency (m/s)
REAL(dp), DIMENSION(KLON)  :: PPRTEN ! total surf precipitation tendency (m/s)
REAL(dp), DIMENSION(KLON)  :: PPRSTEN! solid surf precipitation tendency (m/s)
!
! convective U,V transport (Cu friction)
!LOGICAL,                      INTENT(IN)        :: OUVTRANS ! flag to compute convective
                                                            ! transport for horizonal wind
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT)     :: PUTEN    ! convecctive u tendency (m/s^2)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT)     :: PVTEN    ! convecctive v tendency (m/s^2)

! convective Chemical Tracers:
!LOGICAL,                      INTENT(IN)        :: OCHTRANS ! flag to compute convective
                                                            ! transport for chemical tracer
INTEGER,                    INTENT(IN)            :: KCH1     ! number of species
REAL(dp), DIMENSION(KLON,KLEV,KCH1), INTENT(IN)   :: PCH1     ! grid scale chemical species
REAL(dp), DIMENSION(KLON,KLEV,KCH1), INTENT(OUT):: PCH1TEN  ! chemical convective tendency
                                                              ! (1/s)
!                                        
! Diagnostic variables:
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PPRLFLX! liquid precip flux  (mm/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PPRSFLX! solid precip flux   (mm/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUMF   ! updraft mass flux   (kg/s m2)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDMF   ! downdraft mass flux (kg/s m2)
!!$REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PURV   ! water vapor in updraft (kg/kg)
!!$REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PURCI  ! total condensate in updraft (kg/kg)
REAL(dp), DIMENSION(KLON,KLEV) :: PURV   ! water vapor in updraft (kg/kg)
REAL(dp), DIMENSION(KLON,KLEV) :: PURCI  ! total condensate in updraft (kg/kg)
INTEGER, DIMENSION(KLON),   INTENT(INOUT) :: KCLTOP ! cloud top level (number of model level)
INTEGER, DIMENSION(KLON),   INTENT(INOUT) :: KCLBAS ! cloud base level(number of model level)
                                                      ! they are given a value of
                                                      ! 0 if no convection
REAL(dp), DIMENSION(KLON),      INTENT(OUT)   :: PCAPE  ! CAPE (J/kg)

!
! special for ERA40/MESSy
! updraft detrainment rate   (kg/s m2)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUDR   
! downdraft detrainment rate (kg/s m2)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDDR   
! updraft entrainment rate   (kg/s m2)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUER   
! downdraft entrainment rate (kg/s m2)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDER   
! cv liq. water content (kg/kg)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PURLIQ 
! cv ice water content (kg/kg)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PURICE 
! cv liq prec production (kg/kg)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PURR_P 
! cv ice prec production (kg/kg)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PURS_P 
! 1 for deep convection, else 0
INTEGER, DIMENSION(KLON), INTENT(INOUT) :: ldeep        
REAL(dp), DIMENSION(KLON), INTENT(INOUT)      :: WAT_DIAG
!                                                
!*       0.2   Declarations of local variables :
!
INTEGER  :: JI, JK, JKP, JN  ! loop index
!
REAL(dp), DIMENSION(KLON)               :: ZTIMEC, ZPRLTEN
!
! Local arrays (upside/down) necessary for change of ECMWF arrays to convection arrays
REAL(dp), DIMENSION(KLON,KLEV) :: ZT     ! grid scale T at time t  (K)
REAL(dp), DIMENSION(KLON,KLEV) :: ZRV    ! grid scale water vapor  (kg/kg)
REAL(dp), DIMENSION(KLON,KLEV) :: ZRC    ! grid scale r_c mixing ratio (kg/kg)
REAL(dp), DIMENSION(KLON,KLEV) :: ZRI    ! grid scale r_i mixing ratio (kg/kg)
REAL(dp), DIMENSION(KLON,KLEV) :: ZU     ! grid scale horiz. wind u (m/s) 
REAL(dp), DIMENSION(KLON,KLEV) :: ZV     ! grid scale horiz. wind v (m/s)
REAL(dp), DIMENSION(KLON,KLEV) :: ZW     ! grid scale vertical velocity (m/s)
REAL(dp), DIMENSION(KLON,KLEV) :: ZW1    ! perturbed vertical velocity for ensemble (m/s)
REAL(dp), DIMENSION(KLON,KLEV) :: ZPABS  ! grid scale pressure (Pa)
REAL(dp), DIMENSION(KLON,0:KLEV)::ZPAH   ! grid scale half-level pressure (Pa)
REAL(dp), DIMENSION(KLON,KLEV) :: ZZZ    ! height of model layer (m) 
REAL(dp), DIMENSION(KLON,KLEV) :: ZHSFLX ! turbulent sensible heat flux (W/m^2)
!
REAL(dp), DIMENSION(KLON,KLEV) :: ZTTEN  ! convective temperat. tendency (K/s)
REAL(dp), DIMENSION(KLON,KLEV) :: ZRVTEN ! convective r_v tendency (1/s)
REAL(dp), DIMENSION(KLON,KLEV) :: ZRCTEN ! convective r_c tendency (1/s)
REAL(dp), DIMENSION(KLON,KLEV) :: ZRITEN ! convective r_i tendency (1/s)
REAL(dp), DIMENSION(KLON,KLEV) :: ZUTEN  ! convective u tendency (m/s^2)
REAL(dp), DIMENSION(KLON,KLEV) :: ZVTEN  ! convective m tendency (m/s^2)
REAL(dp), DIMENSION(KLON,KLEV) :: ZUMF   ! updraft mass flux   (kg/s m2)
REAL(dp), DIMENSION(KLON,KLEV) :: ZDMF   ! downdraft mass flux (kg/s m2)
REAL(dp), DIMENSION(KLON,KLEV) :: ZURV   ! water vapor in updrafts (kg/kg)
REAL(dp), DIMENSION(KLON,KLEV) :: ZURCI  ! total condensate in updrafts (kg/kg)
REAL(dp), DIMENSION(KLON,KLEV) :: ZPRLFLX! liquid precip flux  (m/s)
REAL(dp), DIMENSION(KLON,KLEV) :: ZPRSFLX! solid precip flux   (m/s)
INTEGER, DIMENSION(KLON)   :: ICLTOP ! cloud top level (number of model level)
INTEGER, DIMENSION(KLON)   :: ICLBAS ! cloud base level(number of model level)
REAL(dp), DIMENSION(KLON)      :: ZUMFSBAS ! norm. shallow cloud base mass flux (m/s)
REAL(dp), DIMENSION(KLON,KLEV,KCH1):: ZCH1     ! grid scale chemical species
REAL(dp), DIMENSION(KLON,KLEV,KCH1):: ZCH1TEN  ! chemical convective tendency
!!$REAL(dp), DIMENSION(KLON,KLEV,2):: ZCH1     ! grid scale chemical species
!!$REAL(dp), DIMENSION(KLON,KLEV,2):: ZCH1TEN  ! chemical convective tendency
!
! mz_ht_20071219+
REAL(dp), DIMENSION(KLON,KLEV) :: ZURCICE   !updraft ice   (kg/kg)
REAL(dp), DIMENSION(KLON,KLEV) :: ZURCLIQ   !updraft water (kg/kg)
REAL(dp), DIMENSION(KLON,KLEV) :: ZURCRPRO  !updraft rain prod. (kg/kg)
REAL(dp), DIMENSION(KLON,KLEV) :: ZURCSPRO  !updraft snow prod. (kg/kg)
! mz_ht_20071219-
! special for shallow convection
REAL(dp), DIMENSION(:,:), ALLOCATABLE   :: ZTTENS, ZRVTENS, ZRCTENS, ZRITENS, &
                                         & ZUTENS, ZVTENS, ZUMFS, ZURVS, ZURCIS
REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: ZCH1TENS
REAL(dp), DIMENSION(KLON,KLEV) :: ZURCICES   !updraft ice   (kg/kg)
REAL(dp), DIMENSION(KLON,KLEV) :: ZURCLIQS   !updraft water (kg/kg)
INTEGER, DIMENSION(:), ALLOCATABLE  :: ICLBASS, ICLTOPS
!
!*       0.5   Declarations of additional Ensemble fields:
!
INTEGER                             :: KENS     ! number of allowed 
                                                  ! additional deep convection calls
REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: ZTTENE   ! convective temperat. tendency (K/s)
REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: ZRVTENE  ! convective r_v tendency (1/s)
REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: ZRCTENE  ! convective r_c tendency (1/s)
REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: ZRITENE  ! convective r_i tendency (1/s)
REAL(dp), DIMENSION(:,:),   ALLOCATABLE :: ZPRLTENE ! liquid surf precipitation tendency (m/s)
REAL(dp), DIMENSION(:,:),   ALLOCATABLE :: ZPRSTENE ! solid surf precipitation tendency (m/s)
REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: ZUMFE    ! updraft mass flux   (kg/s m2)
REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: ZDMFE    ! downdraft mass flux (kg/s m2)
REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: ZURVE    ! updraft water vapor  (kg/kg)
REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: ZURCIE   ! updraft condensate   (kg/kg)
! mz_ht_20071219+
REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: ZURCICEE  !updraft ice   (kg/kg)
REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: ZURCLIQE  !updraft water (kg/kg)
REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: ZURCRPROE !updraft rain prod. (kg/kg)
REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: ZURCSPROE !updraft snow prod. (kg/kg)
! mz_ht_20071219-

REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: ZPRLFLXE ! liquid precip flux  (m/s)
REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: ZPRSFLXE ! solid precip flux   (m/s)
INTEGER, DIMENSION(:,:),ALLOCATABLE :: ICLTOPE  ! cloud top level (number of model level)
INTEGER, DIMENSION(:,:),ALLOCATABLE :: ICLBASE  ! cloud base level(number of model level)
REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: ZUTENE   ! convective u tendency (m/s^2)
REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: ZVTENE   ! convective u tendency (m/s^2)
REAL(dp), DIMENSION(:,:,:,:),ALLOCATABLE:: ZCH1TENE ! chemical convective tendency
REAL(dp), DIMENSION(:),     ALLOCATABLE :: ZEDUMMY  ! field not to be recomputed by ensemble
INTEGER, DIMENSION(:),  ALLOCATABLE :: IEDUMMY  ! field not to be recomputed by ensemble
REAL(dp), DIMENSION(:),     ALLOCATABLE :: ZWEIGHT  ! weighting factor for ensemble members
REAL(dp)                                :: ZSUM     ! sum of weighting factors
!
! special for ERA40 / MESSy output
REAL(dp), DIMENSION(KLON,KLEV)          :: ZUDR, ZDDR  ! updraft/downdraft detrainment rates (kg/s m3)
REAL(dp), DIMENSION(:,:),   ALLOCATABLE :: ZUDRS
REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: ZUDRE, ZDDRE   

REAL(dp), DIMENSION(KLON,KLEV)          :: ZUER, ZDER  ! updraft/downdraft entrainment rates (kg/s m3)
REAL(dp), DIMENSION(:,:),   ALLOCATABLE :: ZUERS
REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: ZUERE, ZDERE    
!-------------------------------------------------------------------------------
!
!
!*       .9   Setup fundamental thermodunamical/physical constants using ECMWF/ARPEGE routine
!             ------------------------------------------------------------------------------
!
!  CALL SUCST(54,20020211,0,0)
!  already done as parameters in messy_convect_bechtold_param ! mz_ht_16082004
    CALL INITIALIZE_SATW
!
!
!*       1.   Allocate 2D (horizontal, vertical) arrays and additional ensemble arrays
!             ------------------------------------------------------------------------
!
    ALLOCATE( ZTTENS(KLON,KLEV) ) 
    ALLOCATE( ZRVTENS(KLON,KLEV) ) 
    ALLOCATE( ZRCTENS(KLON,KLEV) )
    ALLOCATE( ZRITENS(KLON,KLEV) ) 
    ALLOCATE( ZUTENS(KLON,KLEV) )
    ALLOCATE( ZVTENS(KLON,KLEV) )
    ALLOCATE( ZCH1TENS(KLON,KLEV,KCH1) ) 
    ALLOCATE( ZUMFS(KLON,KLEV) )
    ALLOCATE( ZURCIS(KLON,KLEV) )
    ALLOCATE( ZURVS(KLON,KLEV) )
    ALLOCATE( ICLBASS(KLON) )
    ALLOCATE( ICLTOPS(KLON) )

    ALLOCATE( ZUDRS(KLON,KLEV) )
 
    ALLOCATE( ZUERS(KLON,KLEV) )

    KCLTOP(:)  = 0 ! set default value when no convection
    KCLBAS(:)  = 0 ! can be changed  depending on user
    ICLTOP(:)  = 0 
    ICLBAS(:)  = 0 
    ICLTOPS(:) = 0 
    ICLBASS(:) = 0 

    LDEEP(:)   = 0
    ZCH1TEN(:,:,:) = 0._dp
    PCH1TEN(:,:,:) = 0._dp
!
!
KENS = MIN( KENSM, 3 )
IF ( KENS > 0 ) THEN
    ALLOCATE( ZTTENE(KLON,KLEV,KENS) )
    ALLOCATE( ZRVTENE(KLON,KLEV,KENS) )
    ALLOCATE( ZRCTENE(KLON,KLEV,KENS) )
    ALLOCATE( ZRITENE(KLON,KLEV,KENS) )
    ALLOCATE( ZUMFE(KLON,KLEV,KENS) )
    ALLOCATE( ZDMFE(KLON,KLEV,KENS) )
    ALLOCATE( ZURVE(KLON,KLEV,KENS) )
    ALLOCATE( ZURCIE(KLON,KLEV,KENS) )
    ALLOCATE( ZUTENE(KLON,KLEV,KENS) )
    ALLOCATE( ZVTENE(KLON,KLEV,KENS) )
    ALLOCATE( ZCH1TENE(KLON,KLEV,KCH1,KENS) )
    ALLOCATE( ZPRLFLXE(KLON,KLEV,KENS) )
    ALLOCATE( ZPRSFLXE(KLON,KLEV,KENS) )
    ALLOCATE( ZPRLTENE(KLON,KENS) )
    ALLOCATE( ZPRSTENE(KLON,KENS) )
    ALLOCATE( ICLTOPE(KLON,KENS) )
    ALLOCATE( ICLBASE(KLON,KENS) )
    ALLOCATE( ZEDUMMY(KLON) )
    ALLOCATE( IEDUMMY(KLON) )
    ALLOCATE( ZWEIGHT(KENS) )
    IEDUMMY(:)   = 0
    ICLTOPE(:,:) = 1 ! set default value when no convection
    ICLBASE(:,:) = 1 

  
    ALLOCATE( ZUDRE(KLON,KLEV,KENS) ) 
    ALLOCATE( ZDDRE(KLON,KLEV,KENS) )

    ALLOCATE( ZUERE(KLON,KLEV,KENS) ) 
    ALLOCATE( ZDERE(KLON,KLEV,KENS) )

    ALLOCATE( ZURCICEE(KLON,KLEV,KENS) )
    ALLOCATE( ZURCLIQE(KLON,KLEV,KENS) )
    ALLOCATE( ZURCRPROE(KLON,KLEV,KENS) )
    ALLOCATE( ZURCSPROE(KLON,KLEV,KENS) )
END IF
!
!
!*       2.   Flip arrays upside-down as  first vertical level in convection is 1
!             --------------------------------------------------------------------
!
DO JK = 1, KLEV
   JKP = KLEV - JK + 1
   DO JI = KIDIA, KFDIA
      ZPABS(JI,JKP) = PPABS(JI,JK)
      ZPAH(JI,JKP)  = PPAH(JI,JK)
      ZZZ(JI,JKP)   = PZZ(JI,JK) / RG
      ZT(JI,JKP)    = PT(JI,JK)
!     input for PRV, PRC, PRI is a specific values
      ZRV(JI,JKP)   = PRV(JI,JK) / ( 1.0_dp - PRV(JI,JK) ) ! transform specific humidity
      ZRC(JI,JKP)   = PRC(JI,JK) / ( 1.0_dp - PRC(JI,JK) ) ! in mixing ratio
      ZRI(JI,JKP)   = PRI(JI,JK) / ( 1.0_dp - PRI(JI,JK) ) 
      ZU(JI,JKP)    = PU(JI,JK)
      ZV(JI,JKP)    = PV(JI,JK)
      ZW(JI,JKP)    = PW(JI,JK) 
      ZHSFLX(JI,JKP)= PHSFLX(JI,JK) 
   END DO
END DO
DO JI = KIDIA, KFDIA
      ZPAH(JI,0)   = PPAH(JI,KLEV+1)
END DO

IF ( OCHTRANS ) THEN
   DO JK = 1, KLEV
      JKP = KLEV - JK + 1
      DO JN = 1, KCH1
         DO JI = KIDIA, KFDIA
            ZCH1(JI,JKP,JN) = PCH1(JI,JK,JN)
         END DO
      END DO
   END DO
END IF
! 
!*       4.a  Call deep convection routine
!             ----------------------------
!
IF ( ODEEP ) THEN
!
! 1. Base version
!
    CALL SU_CONVPAR
!
    IF ( OSETTADJ ) ZTIMEC(:) = PTADJD

!
    CALL CONVECT_DEEP( KLON, KLEV, KIDIA, KFDIA, KBDIA, KTDIA,            &
                          PDTCONV, KICE, OREFRESH_ALL, ODOWN, OSETTADJ,   & 
                          ZPABS, ZPAH, ZZZ, PDXDY, ZTIMEC, ZHSFLX,        &
                          ZT, ZRV, ZRC, ZRI, ZU, ZV, ZW,                  &
                          KCOUNT, ZTTEN, ZRVTEN, ZRCTEN, ZRITEN,          &
                          ZPRLTEN, PPRSTEN,                               &
                          ICLTOP, ICLBAS, ZPRLFLX, ZPRSFLX,               & 
                          ZUMF, ZDMF, ZURV, ZURCI, PCAPE,                 &
                          OUVTRANS, ZUTEN, ZVTEN,                         &
                          OCHTRANS, KCH1, ZCH1, ZCH1TEN,                  &
                          ! for ERA40 /MESSy
                          ZUDR, ZDDR, ZUER, ZDER, WAT_DIAG,               &
                          ZURCICE, ZURCLIQ, ZURCRPRO, ZURCSPRO     ) 
!
!  2. Additional Ensemble members
!
    IF ( KENS > 0 ) THEN
!
    CALL SU_CONVPAR1
!
!* first member - changes in YOE_CONVPAR (cloud radius of 500 m)
!                                          specified in SU_CONVPAR1  
!
    CALL CONVECT_DEEP( KLON, KLEV, KIDIA, KFDIA, KBDIA, KTDIA,               &
                          PDTCONV, KICE, OREFRESH_ALL, ODOWN, OSETTADJ,      & 
                          ZPABS, ZPAH, ZZZ, PDXDY, ZTIMEC, ZHSFLX,           &
                          ZT, ZRV, ZRC, ZRI, ZU, ZV, ZW,                     &
                          IEDUMMY, ZTTENE(:,:,1), ZRVTENE(:,:,1),            &
                          ZRCTENE(:,:,1), ZRITENE(:,:,1),                    &
                          ZPRLTENE(:,1), ZPRSTENE(:,1),                      &
                          ICLTOPE(:,1), ICLBASE(:,1), ZPRLFLXE(:,:,1),       &
                          ZPRSFLXE(:,:,1), ZUMFE(:,:,1), ZDMFE(:,:,1),       &
                          ZURVE(:,:,1), ZURCIE(:,:,1), ZEDUMMY,              &
                          OUVTRANS, ZUTENE(:,:,1), ZVTENE(:,:,1),            &
                          OCHTRANS, KCH1, ZCH1, ZCH1TENE(:,:,:,1),           &
                          ! for ERA40 /MESSy
                          ZUDRE(:,:,1), ZDDRE(:,:,1), ZUERE(:,:,1),          &
                          ZDERE(:,:,1), WAT_DIAG,                            &
                          ZURCICEE(:,:,1), ZURCLIQE(:,:,1), ZURCRPROE(:,:,1),&
                          ZURCSPROE(:,:,1)     )  
                         
    END IF
!
    IF ( KENS > 1 ) THEN
!
    CALL SU_CONVPAR
!
!* second member (positive vertical velocity perturb for Trigger)
!
    ZW1=ZW*1.5_dp+1.E-4_dp
    CALL CONVECT_DEEP( KLON, KLEV, KIDIA, KFDIA, KBDIA, KTDIA,               &
                          PDTCONV, KICE, OREFRESH_ALL, ODOWN, OSETTADJ,      & 
                          ZPABS, ZPAH, ZZZ, PDXDY, ZTIMEC, ZHSFLX,           &
                          ZT, ZRV, ZRC, ZRI, ZU, ZV, ZW1,                    &
                          IEDUMMY, ZTTENE(:,:,2), ZRVTENE(:,:,2),            &
                          ZRCTENE(:,:,2), ZRITENE(:,:,2),                    &
                          ZPRLTENE(:,2), ZPRSTENE(:,2),                      &
                          ICLTOPE(:,2), ICLBASE(:,2), ZPRLFLXE(:,:,2),       &
                          ZPRSFLXE(:,:,2), ZUMFE(:,:,2), ZDMFE(:,:,2),       &
                          ZURVE(:,:,2), ZURCIE(:,:,2), ZEDUMMY,              &
                          OUVTRANS, ZUTENE(:,:,2), ZVTENE(:,:,2),            &
                          OCHTRANS, KCH1, ZCH1, ZCH1TENE(:,:,:,2),           &
                          ! for ERA40 /MESSy
                          ZUDRE(:,:,2), ZDDRE(:,:,2), ZUERE(:,:,2),          &
                          ZDERE(:,:,2), WAT_DIAG,                            &
                          ZURCICEE(:,:,2), ZURCLIQE(:,:,2), ZURCRPROE(:,:,2),&
                          ZURCSPROE(:,:,2)     )  
    END IF
!
    IF ( KENS > 2 ) THEN
!
!* third member (negative vertical velocity perturb for Trigger)
!
      ZW1=ZW*0.5_dp-1.E-4_dp
      CALL CONVECT_DEEP( KLON, KLEV, KIDIA, KFDIA, KBDIA, KTDIA,             &
                          PDTCONV, KICE, OREFRESH_ALL, ODOWN, OSETTADJ,      & 
                          ZPABS, ZPAH, ZZZ, PDXDY, ZTIMEC, ZHSFLX,           &
                          ZT, ZRV, ZRC, ZRI, ZU, ZV, ZW1,                    &
                          IEDUMMY, ZTTENE(:,:,3), ZRVTENE(:,:,3),            &
                          ZRCTENE(:,:,3), ZRITENE(:,:,3),                    &
                          ZPRLTENE(:,3), ZPRSTENE(:,3),                      &
                          ICLTOPE(:,3), ICLBASE(:,3), ZPRLFLXE(:,:,3),       &
                          ZPRSFLXE(:,:,3), ZUMFE(:,:,3), ZDMFE(:,:,3),       &
                          ZURVE(:,:,3), ZURCIE(:,:,3), ZEDUMMY,              &
                          OUVTRANS, ZUTENE(:,:,3), ZVTENE(:,:,3),            &
                          OCHTRANS, KCH1, ZCH1, ZCH1TENE(:,:,:,3),           &
                          ! for ERA40 /MESSy
                          ZUDRE(:,:,3), ZDDRE(:,:,3), ZUERE(:,:,3),          &
                          ZDERE(:,:,3), WAT_DIAG,                            &
                          ZURCICEE(:,:,3), ZURCLIQE(:,:,3), ZURCRPROE(:,:,3),&
                          ZURCSPROE(:,:,3)     )  
    END IF
!
  ENDIF
  IF ( .NOT. ODEEP ) THEN
  KCOUNT(:)     =0
  ZTTEN(:,:)    =0.0_dp
  ZRVTEN(:,:)   =0.0_dp
  ZRCTEN(:,:)   =0.0_dp
  ZRITEN(:,:)   =0.0_dp
  ZUTEN(:,:)    =0.0_dp
  ZVTEN(:,:)    =0.0_dp
  ZUMF(:,:)     =0.0_dp
  ZDMF(:,:)     =0.0_dp
  ZURV(:,:)     =0.0_dp
  ZURCI(:,:)    =0.0_dp
  ZCH1TEN(:,:,:)=0.0_dp
  ZPRLTEN(:)    =0.0_dp
  PPRSTEN(:)    =0.0_dp
  ZPRLFLX(:,:)  =0.0_dp
  ZPRSFLX(:,:)  =0.0_dp
  PCAPE(:)      =0.0_dp

  ZUDR(:,:)     =0.0_dp
  ZDDR(:,:)     =0.0_dp

  ZUER(:,:)     =0.0_dp
  ZDER(:,:)     =0.0_dp

  ZURCICE(:,:)  = 0.0_dp
  ZURCLIQ(:,:)  = 0.0_dp
  ZURCRPRO(:,:) = 0.0_dp
  ZURCSPRO(:,:) = 0.0_dp

END IF


do ji=1,klon
  if (maxval(ZUMF(ji,:)).gt.0.0_dp) ldeep(ji) = 1
enddo
!
!*       4.b  Call shallow convection routine
!             -------------------------------
!
IF ( OSHAL ) THEN
!
    IF ( .NOT. ODEEP ) CALL SU_CONVPAR 
    CALL SU_CONVPAR_SHAL
!
    CALL CONVECT_SHALLOW( KLON, KLEV, KIDIA, KFDIA, KBDIA, KTDIA,        &
                           & PDTCONV, KICE, OSETTADJ, PTADJS,            &
                           & ZPABS, ZPAH, ZZZ, PHSFLX(:,KLEV),           &
                           & ZT, ZRV, ZRC, ZRI, ZU, ZV, ZW,              &
                           & KCOUNT, ZTTENS, ZRVTENS, ZRCTENS, ZRITENS,  &
                           & ICLTOPS, ICLBASS, ZUMFS, ZURVS, ZURCIS,     &
                           & OUVTRANS, ZUTENS, ZVTENS,                   &
                           & OCHTRANS, KCH1, ZCH1, ZCH1TENS              &
                           ! for ERA40 /MESSy
                           &,ZUDRS, ZUERS, ZURCICES, ZURCLIQS           )  


ENDIF
IF ( .NOT. OSHAL ) THEN
  ZTTENS(:,:)     = 0.0_dp
  ZRVTENS(:,:)    = 0.0_dp
  ZRCTENS(:,:)    = 0.0_dp
  ZRITENS(:,:)    = 0.0_dp
  ZUTENS(:,:)     = 0.0_dp
  ZVTENS(:,:)     = 0.0_dp
  ZURVS(:,:)      = 0.0_dp
  ZURCIS(:,:)     = 0.0_dp
  ZUMFS(:,:)      = 0.0_dp
  ZCH1TENS(:,:,:) = 0.0_dp
  ZUDRS(:,:)      = 0.0_dp
  ZUERS(:,:)      = 0.0_dp
  ZURCICES(:,:)   = 0.0_dp
  ZURCLIQS(:,:)   = 0.0_dp
END IF
!
!*       5.  Add  - if activated - ensemble average values for deep 
!            and then shallow convective tendencies
!            ---------------------------------------------------------
!
ZSUM = 1.0_dp 
IF ( KENS > 0 ) THEN
    IF ( KENS == 1 ) ZWEIGHT(:) = 0.5_dp
    IF ( KENS >  1 ) ZWEIGHT(:) = 1.0_dp
    DO JN = 1, KENS
    DO JK = 1, KLEV
    DO JI = KIDIA, KFDIA
       ZTTEN(JI,JK)  = ZTTEN(JI,JK)  + ZWEIGHT(JN) * ZTTENE(JI,JK,JN)
       ZRVTEN(JI,JK) = ZRVTEN(JI,JK) + ZWEIGHT(JN) * ZRVTENE(JI,JK,JN)
       ZRCTEN(JI,JK) = ZRCTEN(JI,JK) + ZWEIGHT(JN) * ZRCTENE(JI,JK,JN) 
       ZRITEN(JI,JK) = ZRITEN(JI,JK) + ZWEIGHT(JN) * ZRITENE(JI,JK,JN)
       ZPRLFLX(JI,JK)= ZPRLFLX(JI,JK)+ ZWEIGHT(JN) * ZPRLFLXE(JI,JK,JN)
       ZPRSFLX(JI,JK)= ZPRSFLX(JI,JK)+ ZWEIGHT(JN) * ZPRSFLXE(JI,JK,JN)
       ZUMF(JI,JK)   = ZUMF(JI,JK)   + ZWEIGHT(JN) * ZUMFE(JI,JK,JN)
       ZDMF(JI,JK)   = ZDMF(JI,JK)   + ZWEIGHT(JN) * ZDMFE(JI,JK,JN)
       ZURV(JI,JK)   = ZURV(JI,JK)   + ZWEIGHT(JN) * ZURVE(JI,JK,JN)
       ZURCI(JI,JK)  = ZURCI(JI,JK)  + ZWEIGHT(JN) * ZURCIE(JI,JK,JN)

       ZUDR(JI,JK)   = ZUDR(JI,JK)   + ZWEIGHT(JN) * ZUDRE(JI,JK,JN) ! for ERA40
       ZDDR(JI,JK)   = ZDDR(JI,JK)   + ZWEIGHT(JN) * ZDDRE(JI,JK,JN)
       ZUER(JI,JK)   = ZUER(JI,JK)   + ZWEIGHT(JN) * ZUERE(JI,JK,JN) ! for ERA40
       ZDER(JI,JK)   = ZDER(JI,JK)   + ZWEIGHT(JN) * ZDERE(JI,JK,JN)

       ZURCICE(JI,JK)  = ZURCICE(JI,JK)  + ZWEIGHT(JN) * ZURCICEE(JI,JK,JN)
       ZURCLIQ(JI,JK)  = ZURCLIQ(JI,JK)  + ZWEIGHT(JN) * ZURCLIQE(JI,JK,JN)
       ZURCRPRO(JI,JK) = ZURCRPRO(JI,JK) + ZWEIGHT(JN) * ZURCRPROE(JI,JK,JN)
       ZURCSPRO(JI,JK) = ZURCSPRO(JI,JK) + ZWEIGHT(JN) * ZURCSPROE(JI,JK,JN)
    END DO
    END DO
    DO JI = KIDIA, KFDIA
       ZPRLTEN(JI)  = ZPRLTEN(JI)  + ZWEIGHT(JN) * ZPRLTENE(JI,JN)
       PPRSTEN(JI)  = PPRSTEN(JI)  + ZWEIGHT(JN) * ZPRSTENE(JI,JN)
       ICLTOP(JI)   = MAX(ICLTOP(JI), ICLTOPE(JI,JN))
       ICLBAS(JI)   = MAX(ICLBAS(JI), ICLBASE(JI,JN))
    END DO
       IF ( OUVTRANS ) THEN
         DO JK = 1, KLEV
         DO JI = KIDIA, KFDIA
            ZUTEN(JI,JK) = ZUTEN(JI,JK) + ZWEIGHT(JN) * ZUTENE(JI,JK,JN) 
            ZVTEN(JI,JK) = ZVTEN(JI,JK) + ZWEIGHT(JN) * ZVTENE(JI,JK,JN) 
         END DO
         END DO
       END IF
       IF ( OCHTRANS )  THEN
         DO JK = 1, KLEV
         DO JI = KIDIA, KFDIA
          ZCH1TEN(JI,JK,:) = ZCH1TEN(JI,JK,:) + ZWEIGHT(JN) * ZCH1TENE(JI,JK,:,JN)
         END DO
         END DO
       END IF
    END DO
!
    ZSUM = 1.0_dp / ( 1.0_dp + SUM( ZWEIGHT(:) ) )
END IF
!
    DO JK = 1, KLEV
    DO JI = KIDIA, KFDIA
     IF( ZPRLTEN(JI) > 0.0_dp ) THEN
       ZTTENS(JI,JK) = 0.0_dp
       ZRVTENS(JI,JK)= 0.0_dp
       ZRCTENS(JI,JK)= 0.0_dp
       ZRITENS(JI,JK)= 0.0_dp
       ZUMFS(JI,JK)  = 0.0_dp
       ZURVS(JI,JK)  = 0.0_dp
       ZURCIS(JI,JK) = 0.0_dp
       ZUDRS(JI,JK)  = 0.0_dp
       ZUERS(JI,JK)  = 0.0_dp
       ZURCICES(JI,JK)  = 0.0_dp
       ZURCLIQS(JI,JK)  = 0.0_dp
     END IF
       ZTTEN(JI,JK)  = ZTTEN(JI,JK)  * ZSUM + ZTTENS(JI,JK) 
       ZRVTEN(JI,JK) = ZRVTEN(JI,JK) * ZSUM + ZRVTENS(JI,JK)
       ZRCTEN(JI,JK) = ZRCTEN(JI,JK) * ZSUM + ZRCTENS(JI,JK)
       ZRITEN(JI,JK) = ZRITEN(JI,JK) * ZSUM + ZRITENS(JI,JK)
       ZPRLFLX(JI,JK)= ZPRLFLX(JI,JK)* ZSUM
       ZPRSFLX(JI,JK)= ZPRSFLX(JI,JK)* ZSUM
       ZUMF(JI,JK)   = ZUMF(JI,JK)   * ZSUM + ZUMFS(JI,JK)
       ZDMF(JI,JK)   = ZDMF(JI,JK)   * ZSUM
       ZURV(JI,JK)   = ZURV(JI,JK)   * ZSUM + ZURVS(JI,JK)
       ZURCI(JI,JK)  = ZURCI(JI,JK)  * ZSUM + ZURCIS(JI,JK)

       ZUDR(JI,JK)   = ZUDR(JI,JK)   * ZSUM + ZUDRS(JI,JK)
       ZUER(JI,JK)   = ZUER(JI,JK)   * ZSUM + ZUERS(JI,JK)

       ZURCICE(JI,JK)  = ZURCICE(JI,JK)  * ZSUM + ZURCICES(JI,JK)
       ZURCLIQ(JI,JK)  = ZURCLIQ(JI,JK)  * ZSUM + ZURCLIQS(JI,JK)
       ZURCRPRO(JI,JK) = ZURCRPRO(JI,JK) * ZSUM 
       ZURCSPRO(JI,JK) = ZURCSPRO(JI,JK) * ZSUM 
    END DO
    END DO
    DO JI = KIDIA, KFDIA
       PPRTEN(JI)   = ( ZPRLTEN(JI) + PPRSTEN(JI) ) * ZSUM 
       PPRSTEN(JI)  = PPRSTEN(JI)  * ZSUM
       ICLTOP(JI)   = MAX(ICLTOP(JI), ICLTOPS(JI))
       ICLBAS(JI)   = MAX(ICLBAS(JI), ICLBASS(JI))
    END DO
       IF ( OUVTRANS ) THEN
         DO JK = 1, KLEV
         DO JI = KIDIA, KFDIA
            ZUTEN(JI,JK) = ZUTEN(JI,JK) * ZSUM + ZUTENS(JI,JK)
            ZVTEN(JI,JK) = ZVTEN(JI,JK) * ZSUM + ZVTENS(JI,JK)
         END DO
         END DO
       END IF
       IF ( OCHTRANS ) THEN
         DO JK = 1, KLEV
         DO JI = KIDIA, KFDIA
            ZCH1TEN(JI,JK,:) = ZCH1TEN(JI,JK,:) * ZSUM + ZCH1TENS(JI,JK,:)
         END DO
         END DO
       END IF
!
!
!*       6.  Reflip arrays to ECMWF/ARPEGE vertical structure
!            change mixing ratios to sepcific humidity
!
   DO JK = 1, KLEV
   JKP = KLEV - JK + 1 
   DO JI = KIDIA, KFDIA
      PTTEN(JI,JK)  = ZTTEN(JI,JKP)
!      PRVTEN(JI,JK) = ZRVTEN(JI,JKP) / ( 1.0_dp + ZRV(JI,JKP) ) ** 2
!      PRCTEN(JI,JK) = ZRCTEN(JI,JKP) / ( 1.0_dp + ZRC(JI,JKP) ) ** 2
!      PRITEN(JI,JK) = ZRITEN(JI,JKP) / ( 1.0_dp + ZRI(JI,JKP) ) ** 2
      PRVTEN(JI,JK) = ZRVTEN(JI,JKP) * ( 1.0_dp - PRV(JI,JK) ) ** 2
      PRCTEN(JI,JK) = ZRCTEN(JI,JKP) * ( 1.0_dp - PRC(JI,JK) ) ** 2
      PRITEN(JI,JK) = ZRITEN(JI,JKP) * ( 1.0_dp - PRI(JI,JK) ) ** 2

      PUTEN(JI,JK)  = ZUTEN(JI,JKP)
      PVTEN(JI,JK)  = ZVTEN(JI,JKP)
      PUMF(JI,JK)   = ZUMF(JI,JKP)
      PDMF(JI,JK)   = ZDMF(JI,JKP)
      PURV(JI,JK)   = ZURV(JI,JKP) / ( 1.0_dp + ZURV(JI,JKP) )
      PURCI(JI,JK)  = ZURCI(JI,JKP)/ ( 1.0_dp + ZURCI(JI,JKP) )
      PPRLFLX(JI,JK)= ZPRLFLX(JI,JKP)*1000._dp
      PPRSFLX(JI,JK)= ZPRSFLX(JI,JKP)*1000._dp
!      PPRLFLX(JI,JK)= ZPRLFLX(JI,JKP)
!      PPRSFLX(JI,JK)= ZPRSFLX(JI,JKP)
      PUDR(JI,JK)   = ZUDR(JI,JKP)
      PUER(JI,JK)   = ZUER(JI,JKP)
      PDDR(JI,JK)   = -1._dp * ZDDR(JI,JKP)
      PDER(JI,JK)   = -1._dp * ZDER(JI,JKP)

      PURLIQ(JI,JK) = ZURCLIQ(JI,JKP)
      PURICE(JI,JK) = ZURCICE(JI,JKP)
      PURR_P(JI,JK) = ZURCRPRO(JI,JKP)
      PURS_P(JI,JK) = ZURCSPRO(JI,JKP)
   END DO
END DO

DO JI = KIDIA, KFDIA
   JK = ICLTOP(JI)
   KCLTOP(JI) = KLEV - JK + 1
   JK = ICLBAS(JI)
   KCLBAS(JI) = KLEV - JK + 1
!   IF ( ICLTOP(JI) == 1 ) KCLTOP(JI) = KLEV+1 ! 1 original ! mz_ht_20050914
!   IF ( ICLBAS(JI) == 1 ) KCLBAS(JI) = KLEV+1 ! 1 original ! mz_ht_20050914
   ! op_mm_20140324+
   ! change for use in COSMO
   ! If no convection occures klev=0
   ! THIS SHOULD NOT BE CHANGED. OTHERWISE COSMO GETS IN TROUBLE 
   IF ( ICLTOP(JI) == 1 ) KCLTOP(JI) = 0 
   IF ( ICLBAS(JI) == 1 ) KCLBAS(JI) = 0 
   IF ( ICLTOP(JI) == 0 ) KCLTOP(JI) = 0 
   IF ( ICLBAS(JI) == 0 ) KCLBAS(JI) = 0 
   ! op_mm_20140324-

END DO

IF ( OCHTRANS ) THEN
   DO JK = 1, KLEV
      JKP = KLEV - JK + 1
      DO JN = 1, KCH1
         DO JI = KIDIA, KFDIA
            PCH1TEN(JI,JK,JN) = ZCH1TEN(JI,JKP,JN)
         END DO
      END DO
   END DO
END IF

   
!*       7.  Deallocate local arrays
!
       DEALLOCATE( ICLBASS )
       DEALLOCATE( ICLTOPS )
       DEALLOCATE( ZUMFS )
       DEALLOCATE( ZURVS )
       DEALLOCATE( ZURCIS )
       DEALLOCATE( ZCH1TENS ) 
       DEALLOCATE( ZRCTENS )
       DEALLOCATE( ZRITENS ) 
       DEALLOCATE( ZTTENS )
       DEALLOCATE( ZRVTENS ) 
       DEALLOCATE( ZUDRS )
       DEALLOCATE( ZUERS )
       DEALLOCATE( ZUTENS ) 
       DEALLOCATE( ZVTENS ) 

IF ( KENS > 0 ) THEN
       DEALLOCATE( ZTTENE )
       DEALLOCATE( ZRVTENE )
       DEALLOCATE( ZRCTENE )
       DEALLOCATE( ZRITENE )
       DEALLOCATE( ZUTENE )
       DEALLOCATE( ZVTENE )
       DEALLOCATE( ZUMFE )
       DEALLOCATE( ZDMFE )
       DEALLOCATE( ZURVE )
       DEALLOCATE( ZURCIE )
       DEALLOCATE( ZCH1TENE )
       DEALLOCATE( ZPRLFLXE )
       DEALLOCATE( ZPRSFLXE )
       DEALLOCATE( ZPRLTENE )
       DEALLOCATE( ZPRSTENE )
       DEALLOCATE( ZEDUMMY )
       DEALLOCATE( IEDUMMY )
       DEALLOCATE( ICLTOPE )
       DEALLOCATE( ICLBASE )
       DEALLOCATE( ZWEIGHT )

       DEALLOCATE( ZUDRE )
       DEALLOCATE( ZDDRE )
       DEALLOCATE( ZUERE )
       DEALLOCATE( ZDERE )
END IF

!
!
!END SUBROUTINE CONVECTION
!
!----------------------------------------------------------------------------

END SUBROUTINE bechtold_cumastr

!=========================================================================

SUBROUTINE EMANUEL_CUMASTR(KPROMA,  NLEV, NTRAC, ZTMST,                    &
                           PPABS,   PPAH,                                  &
                           PT,      PQ,   PU,    PV,                       &
                           PXTP1,   PQSAT,                                 &
                           ! OUTPUT
                           itype,   PRECIP,                                &
                           PTTE,    PQTE, PUTE,  PVTE,                     &
                           PXTTE,                                          &
                           ZWD,     TPRIME,                                &
                           QPRIME,  CBMF,                                  &
                           kbot,    ktop, PUMF,  PDMF,                     &
                           PREC,    PRECFORM,                              &
                           PUER,    PUDR, PDER,  PDDR,  PCAPE,             &
                           LWC,     IWC,  RFORM, SFORM )


  USE MESSY_CONVECT_EMANUEL

  IMPLICIT NONE
! IN - OUT - Variables 

  INTEGER,  INTENT(IN) :: KPROMA, NLEV, NTRAC
  REAL(dp), INTENT(IN) :: PPABS(kproma,nlev), PPAH(kproma, nlev+1), ztmst

  REAL(dp) :: PT(kproma,nlev), PQ(kproma,nlev),    &
              PU(kproma,nlev), PV(kproma,nlev), PQSAT(kproma,nlev)
  REAL(dp) :: PXTP1(kproma,nlev,ntrac)

  INTEGER  :: ITYPE(kproma)
  REAL(dp) :: PRECIP(kproma), CBMF(kproma)
  REAL(dp) :: PTTE(kproma,nlev), PQTE(kproma,nlev),    &
              PUTE(kproma,nlev), PVTE(kproma,nlev)

  REAL(dp) :: PXTTE(kproma,nlev,ntrac)
  REAL(dp) :: ZWD(kproma), TPRIME(kproma), QPRIME(kproma)
  INTEGER  :: kbot(kproma), ktop(kproma)
  REAL(dp) :: PUMF(kproma,nlev), PDMF(kproma,nlev)
  REAL(dp) :: PUER(kproma,nlev), PUDR(kproma,nlev)
  REAL(dp) :: PDER(kproma,nlev), PDDR(kproma,nlev)
  REAL(dp) :: PREC(kproma,nlev), PRECFORM(kproma,nlev)
  REAL(dp) :: PCAPE(kproma)
  REAL(dp) :: LWC(kproma,nlev),  IWC(kproma,nlev)
  REAL(dp) :: RFORM(kproma,nlev),SFORM(kproma,nlev)
! Local variables

  REAL(dp) :: ZPABS(kproma,nlev), ZPAH(kproma, nlev+1)

  REAL(dp) :: ZT(kproma,nlev),    ZQ(kproma,nlev), QSAT(kproma,nlev),   &
              ZU(kproma,nlev),    ZV(kproma,nlev)
  REAL(dp) :: ZTTE(kproma,nlev),  ZQTE(kproma,nlev),    &
              ZUTE(kproma,nlev),  ZVTE(kproma,nlev)
  REAL(dp) :: ZUMF(kproma,nlev),  ZDMF(kproma,nlev)
  REAL(dp) :: ZUER(kproma,nlev),  ZUDR(kproma,nlev)
  REAL(dp) :: ZDER(kproma,nlev),  ZDDR(kproma,nlev)
  REAL(dp) :: ZPREC(kproma,nlev), ZFORM(kproma,nlev)
  REAL(dp) :: ZLWC(kproma,nlev),  ZIWC(kproma,nlev)
  REAL(dp) :: ZRFORM(kproma,nlev),ZSFORM(kproma,nlev)
  REAL(dp), ALLOCATABLE :: ZXTP1(:,:,:)
  REAL(dp), ALLOCATABLE :: ZXTTE(:,:,:)

  INTEGER  :: JL, JK, JKP, JI
  INTEGER  :: MINORIG

! the arrays have to be flipped, since for the convection calculation 
! the lowest model level is number one, in contrast to ECHAM5 where 
! the lowest model level is number nlev
  DO JK = 1, NLEV
    JKP = NLEV - JK + 1
    DO JL = 1, KPROMA
      ZPABS(JL,JKP)         = PPABS(JL,JK)/100._dp
      ZT(JL,JKP)            = PT(JL,JK)
      ZQ(JL,JKP)            = PQ(JL,JK)
      ZU(JL,JKP)            = PU(JL,JK)
      ZV(JL,JKP)            = PV(JL,JK)
      QSAT(JL,JKP)          = PQSAT(JL,JK)
    END DO
  END DO
  DO JK = 1, NLEV+1
    JKP = NLEV+1 - JK + 1
    DO JL = 1, KPROMA
      ZPAH(JL,JKP)          = PPAH(JL,JK)/100._dp
    END DO
  END DO

  ZTTE(:,:)    = 0._dp
  ZQTE(:,:)    = 0._dp
  ZUTE(:,:)    = 0._dp
  ZVTE(:,:)    = 0._dp
  ZUMF(:,:)    = 0._dp
  ZDMF(:,:)    = 0._dp
  ZPREC(:,:)   = 0._dp
  ZFORM(:,:)   = 0._dp
  ZUER(:,:)    = 0._dp
  ZUDR(:,:)    = 0._dp
  ZDER(:,:)    = 0._dp
  ZDDR(:,:)    = 0._dp
  ZRFORM(:,:)  = 0._dp
  ZSFORM(:,:)  = 0._dp
  ZLWC(:,:)    = 0._dp
  ZIWC(:,:)    = 0._dp

  ALLOCATE(ZXTP1(kproma,nlev,ntrans))
  !ALLOCATE(ZXTTE(kproma,nlev,ntrans))
  ! op_mm_20140327+
  !compiler workaround (G95 0.92, 0.93)
  ALLOCATE(ZXTTE(kproma,nlev,0:ntrans))
  ! op_mm_20140327-

  ZXTTE(:,:,:) = 0._dp
  ZXTP1(:,:,:) = 0._dp

  IF (ltransport) THEN
    DO JK = 1, NLEV
      JKP = NLEV - JK + 1
      DO JL = 1, KPROMA
        ZXTP1(JL,JKP,1:NTRANS) = PXTP1(JL,JK,1:NTRANS)
      END DO
    END DO
  ENDIF
!
! the convection routine is coded as a single column model
! therefore it must be called sequentially for all columns of the model
! WARNING: very bad performance on vector machines is to be expected !!!!

  DO JL = 1,kproma
    DO JI = 1, nlev
      MINORIG = JI
!      print*, "begin_cv", jl, JI
      CALL CONVECT_EMANUEL(NLEV+1,                                 &
                           ZT(JL,:),    ZQ(JL,:),   QSAT(JL,:),    &
                           ZU(JL,:),    ZV(JL,:),   ZXTP1(JL,:,:), &
                           ZPABS(JL,:), ZPAH(JL,:),                &
                           NLEV,        NLEV-1,     NTRANS,        &
                           ztmst,       ITYPE(JL),                 &
                           ZTTE(JL,:),  ZQTE(JL,:),                &
                           ZUTE(JL,:),  ZVTE(JL,:), ZXTTE(JL,:,:), &
                           PRECIP(JL),                             &
                           ZWD(JL),     TPRIME(JL), QPRIME(JL),    &
                           CBMF(JL),                               &
! extra diagnostics ! mz_ht_20070214
                           MINORIG,     kbot(jl),   ktop(jl),      &
                           ZUMF(JL,:),  ZDMF(JL,:),                &
                           ZPREC(JL,:), ZFORM(JL,:),               &
                           ZUER(JL,:),  ZUDR(JL,:),                &
                           ZDER(JL,:),  ZDDR(JL,:), PCAPE(JL),     &
                           ZLWC(JL,:),  ZIWC(JL,:),                &
                           ZRFORM(JL,:),ZSFORM(JL,:) )
      IF (ITYPE(jl) /= 4) exit
    enddo
!    print*, "end cv", kbot(jl), ktop(jl), precip(jl), ztte(jl,:)*ztmst
  ENDDO

! for the output the levels have to be flipped again to be consistent 
! with the notation of ECHAM5, and for a correct assignment of tendencies, etc.
  DO JK = 1, NLEV
    JKP = NLEV - JK + 1 
    DO JL = 1, KPROMA
      PTTE(JL,JK)          = ZTTE(JL,JKP)
      PQTE(JL,JK)          = ZQTE(JL,JKP)
      PUTE(JL,JK)          = ZUTE(JL,JKP)
      PVTE(JL,JK)          = ZVTE(JL,JKP)
      PUMF(JL,JK)          = ZUMF(JL,JKP)
      PDMF(JL,JK)          = -ZDMF(JL,JKP)
      PREC(JL,JK)          = ZPREC(JL,JKP)/86400._dp
      PRECFORM(JL,JK)      = ZFORM(JL,JKP)/86400._dp
      PUDR(JL,JK)          = ZUDR(JL,JKP)
      PUER(JL,JK)          = ZUER(JL,JKP)
      PDDR(JL,JK)          = -1._dp * ZDDR(JL,JKP)
      PDER(JL,JK)          = -1._dp * ZDER(JL,JKP)
      LWC(JL,JK)           = ZLWC(JL,JKP)
      IWC(JL,JK)           = ZIWC(JL,JKP)
      RFORM(JL,JK)         = ZRFORM(JL,JKP)
      SFORM(JL,JK)         = ZSFORM(JL,JKP)
    END DO
  END DO
  
  IF (ltransport) THEN
    DO JK = 1, NLEV
      JKP = NLEV - JK + 1
      DO JL = 1, KPROMA
        PXTTE(JL,JK,1:NTRANS) = ZXTTE(JL,JKP,1:NTRANS)
      END DO
    END DO
  ENDIF
  DEALLOCATE(ZXTP1)
  DEALLOCATE(ZXTTE)

  ! transform precip from mm/day to kg/m^2s = mm/s
  PRECIP(1:KPROMA) = PRECIP(1:KPROMA) / 86400._dp
 
  ! flip index of cloud base and cloud top
  DO jl=1,kproma
    kbot(jl) = NLEV + 1 - kbot(jl)
    ktop(jl) = NLEV + 1 - ktop(jl)
 
    ! op_mm_20140324+
   ! change for use in COSMO
   ! If no convection occures klev=0
   ! THIS IS SHOULD NOT BE CHANGED. OTHERWISE COSMO GETS IN TROUBLE    
   IF (  kbot(jl) == nlev)     kbot(jl)= 0 
   IF (  ktop(jl) == nlev)     ktop(jl)= 0 
   IF (  kbot(jl) == nlev+1)   kbot(jl)= 0 
   IF (  ktop(jl) == nlev+1)   ktop(jl)= 0 
   ! op_mm_20140324-

  ENDDO

END SUBROUTINE EMANUEL_CUMASTR

!=============================================================================
SUBROUTINE DONNER_CUMASTR(is, ie, js, je, dt, temp, mixing_ratio, pfull, &
                        phalf, omega, land, sfc_sh_flux, sfc_vapor_flux,&
                        tr_flux, tracers, Time, cell_cld_frac,  &
                        cell_liq_amt, cell_liq_size, cell_ice_amt,   &
                        cell_ice_size, meso_cld_frac, meso_liq_amt, &
                        meso_liq_size, meso_ice_amt, meso_ice_size,  &
                        nsum, precip, delta_temp, delta_vapor, detf, &
                        uceml_inter, mtot, donner_humidity_area,    &
                        donner_humidity_ratio, qtrtnd, &
                        ! mz_ht_20070421+                        
                        p_pe, exit_flag,               &
                        ! mz_ht_20070421-
                        qlin, qiin, qain,              &      ! optional
                        delta_ql, delta_qi, delta_qa)         ! optional

                        
!-------------------------------------------------------------------
!    donner_deep is the prognostic driver subroutine of donner_deep_mod.
!    it takes as input the temperature (temp), vapor mixing ratio 
!    (mixing_ratio), pressure at full and half-levels (pfull, phalf),
!    vertical velocity at full levels (omega), the large scale cloud 
!    variables (qlin, qiin, qain), the land fraction (land),  the heat 
!    (sfc_sh_flux) , moisture (sfc_vapor_flux) and tracer (tr_flux) 
!    fluxes across the surface that are to be seen by this parameter-
!    ization, the tracers to be transported by the donner convection
!    parameterization (tracers), and the current time (as time_type 
!    variable Time). the routine returns the precipitation (precip),
!    increments to the temperature (delta_temp) and mixing ratio 
!    (delta_vapor), the detrained mass flux (detf), upward cell mass 
!    flux at interface levels  (uceml_inter) and total mass flux at full
!    levels (mtot), two arrays needed to connect the donner convection 
!    and strat cloud parameterizations (donner_humidity_area, 
!    donner_humidity_ratio), increments to the cloudwater (delta_ql), 
!    cloudice (delta_qi) and cloud area (delta_qa) fields and tendencies
!    for those tracers that are to be transported by the donner convect-
!    ion parameterization (qtrtnd). there are an additional eleven arrays
!    defining the donner scheme cloud characteristics needed by the rad-
!    iation package, which are passed in and updated on donner calcul-
!    ation steps.
!-------------------------------------------------------------------
  ! mz_ht_20070421+
  USE MESSY_CONVECT_DONNER_TYPES_MOD,       ONLY: donner_initialized_type, &
                                  donner_save_type, donner_rad_type, &
                                  donner_nml_type, donner_param_type, &
                                  donner_column_diag_type, &
                                  donner_conv_type, donner_cape_type
  USE MESSY_CONVECT_DONNER_ADDITIONS,       ONLY: MODULE_IS_INITIALIZED,       &
                                  PARAM, COL_DIAG, NML, DON_SAVE, INITIALIZED, &
                                  NLEV_HIRES
  USE MESSY_CONVECT_DONNER_DEEP_K,          ONLY: don_d_donner_deep_k,         &
                                  don_d_dealloc_loc_vars_k
  ! mz_ht_20070421-

!--------------------------------------------------------------------
integer,                      intent(in)    :: is, ie, js, je
real,                         intent(in)    :: dt
real, dimension(:,:,:),       intent(in)    :: temp, mixing_ratio, &
                                               pfull, phalf, omega
real, dimension(:,:),         intent(in)    :: land
real, dimension(:,:),         intent(in)    :: sfc_sh_flux, &
                                               sfc_vapor_flux
real, dimension(:,:,:),       intent(in)    :: tr_flux 
real, dimension(:,:,:,:),     intent(in)    :: tracers 
! mz_ht_20070421+
!type(time_type),              intent(in)    :: Time
REAL, DIMENSION(:),           intent(in)    :: Time
! mz_ht_20070421-
real, dimension(:,:,:),       intent(inout) :: cell_cld_frac,  &
                                               cell_liq_amt,  &
                                               cell_liq_size, &
                                               cell_ice_amt,  &
                                               cell_ice_size, &
                                               meso_cld_frac,  &
                                               meso_liq_amt, &
                                               meso_liq_size, &
                                               meso_ice_amt,   &
                                               meso_ice_size
integer, dimension(:,:),      intent(inout) :: nsum
real, dimension(:,:),         intent(out)   :: precip      
real, dimension(:,:,:),       intent(out)   :: delta_temp, delta_vapor,&
                                               detf, uceml_inter, mtot, &
                                               donner_humidity_area,&
                                               donner_humidity_ratio
real, dimension(:,:,:,:),     intent(out)   :: qtrtnd 
integer,                      intent(in)    :: p_pe
real, dimension(:,:,:),       intent(in),                &
                                   optional :: qlin, qiin, qain
real, dimension(:,:,:),       intent(out),               &
                                   optional :: delta_ql, delta_qi, &
                                               delta_qa
LOGICAL, dimension(:,:),      intent(inout) :: exit_flag
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!   intent(in) variables:
!
!     is, ie         first and last values of i index values of points 
!                    in this physics window (processor coordinates)
!     js, je         first and last values of j index values of points 
!                    in this physics window (processor coordinates)
!     dt             physics time step [ sec ]
!     temp           temperature field at model levels [ deg K ]
!     mixing_ratio   vapor mixing ratio field at model levels 
!                    [ kg(h20) / kg(dry air) ]
!     pfull          pressure field on model full levels [ Pa ]
!     phalf          pressure field at half-levels 1:nlev+1  [ Pa ]
!     omega          model omega field at model full levels [ Pa / sec ]
!     qlin           large-scale cloud liquid specific humidity 
!                    [ kg(h2o) / kg (moist air) ]
!     qiin           large-scale cloud ice specific humidity 
!                    [ kg(h2o) / kg (moist air) ]
!     qain           large-scale cloud fraction  
!                    [ fraction ]
!     land           fraction of grid box covered by land
!                    [ fraction ]
!     sfc_sh_flux    sensible heat flux across the surface
!                    [ watts / m**2 ]
!     sfc_vapor_flux water vapor flux across the surface
!                    [ kg(h2o) / (m**2 sec) ]
!     tr_flux        surface flux of tracers transported by
!                    donner_deep_mod [ kg(tracer) / (m**2 sec) ]
!     tracers        tracer mixing ratios
!                    [ kg(tracer) / kg (dry air) ]
!     Time           current time (time_type)
!
!   intent(out) variables:
!
!     precip         precipitation generated by deep convection
!                    [ kg(h2o) / m**2 ]
!     delta_temp     temperature increment due to deep convection 
!                    [ deg K ]
!     delta_vapor    water vapor mixing ratio increment due to deep 
!                    convection [ kg(h2o) / kg (dry air) ]
!     detf           detrained cell mass flux at model levels 
!                    [ (kg / (m**2 sec) ) ]
!     uceml_inter    upward cell mass flux at interface levels 
!                    [ (kg / (m**2 sec) ) ]
!     mtot           mass flux at model full levels, convective plus 
!                    mesoscale, due to donner_deep_mod 
!                    [ (kg / (m**2 sec) ) ]
!     donner_humidity_area
!                    fraction of grid box in which humidity is affected
!                    by the deep convection, defined as 0.0 below cloud
!                    base and above the mesoscale updraft, and as the
!                    sum of the cell and mesoscale cloud areas in 
!                    between. it is used in strat_cloud_mod to determine
!                    the large-scale specific humidity field for the
!                    grid box. DO NOT use for radiation calculation,
!                    since not all of this area includes condensate.
!                    [ fraction ]
!     donner_humidity_ratio
!                    ratio of large-scale specific humidity to specific 
!                    humidity in environment outside convective system
!                    [ dimensionless ]
!     delta_ql       cloud water specific humidity increment due to 
!                    deep convection over the timestep
!                    [ kg (h2o) / kg (moist air) ]
!     delta_qi       cloud ice specific humidity increment due to deep 
!                    convection over the timestep 
!                    [ kg (h2o) / kg (moist air) ]
!     delta_qa       cloud area increment due to deep convection
!                    over the time step [ fraction ]
!     qtrtnd         tracer time tendencies due to deep convection
!                    during the time step
!                    [ kg(tracer) / (kg (dry air) sec) ]
!
!   intent(inout) variables:
!
!     cell_cld_frac  fractional coverage of convective cells in
!                    grid box [ dimensionless ]
!     cell_liq_amt   liquid water content of convective cells
!                    [ kg(h2o) / kg(air) ]
!     cell_liq_size  assumed effective size of cell liquid drops
!                    [ microns ]
!     cell_ice_amt   ice water content of cells
!                    [ kg(h2o) / kg(air) ]
!     cell_ice_size  generalized effective diameter for ice in
!                    convective cells [ microns ]
!     meso_cld_frac  fractional area of mesoscale clouds in grid
!                    box [ dimensionless ]
!     meso_liq_amt   liquid water content in mesoscale clouds
!                    [ kg(h2o) / kg(air) ]
!     meso_liq_size  assumed effective size of mesoscale drops
!                    [ microns ]
!     meso_ice_amt   ice water content of mesoscale elements
!                    [ kg(h2o) / kg(air) ]
!     meso_ice_size  generalized ice effective size for anvil ice
!                    [ microns ]
!     nsum           number of time levels over which the above variables
!                    have so far been summed
!
!--------------------------------------------------------------------


!--------------------------------------------------------------------
!    local variables:

      real,    dimension (size(temp,1), size(temp,2), size(temp,3)) :: &
                       temperature_forcing, moisture_forcing, pmass, &
                       qlin_arg, qiin_arg, qain_arg, delta_ql_arg, & 
                       delta_qi_arg, delta_qa_arg

      real,    dimension (size(temp,1), size(temp,2)) ::                &
                       parcel_rise, total_precip

      type(donner_conv_type)            :: Don_conv
      type(donner_cape_type)            :: Don_cape
      type(donner_rad_type)             :: Don_rad
      character(len=128)                :: ermesg
      integer                           :: isize, jsize, nlev_lsm
      integer                           :: ntr, me
      logical                           :: calc_conv_on_this_step 
      logical                           :: cloud_tracers_present
      integer                           :: num_cld_tracers
      integer                           :: k   

!--------------------------------------------------------------------
!   local variables:
!
!     temperature_forcing  temperature tendency due to donner convection
!                          [ deg K / sec ]
!     moisture_forcing     vapor mixing ratio tendency due to donner 
!                          convection [ kg(h2o) / (kg(dry air) sec ) ]
!     pmass                mass per unit area within the grid box
!                          [ kg (air) / (m**2) ]
!     parcel_rise          accumulated vertical displacement of a 
!                          near-surface parcel as a result of the lowest
!                          model level omega field [ Pa ]
!     total_precip         total precipitation rate produced by the
!                          donner parameterization [ mm / day ]
!     exit_flag            logical array indicating whether deep conv-
!                          ection exists in a column
!     Don_conv             donner_convection_type derived type variable 
!                          containing diagnostics and intermediate
!                          results describing the nature of the convec-
!                          tion produced by the donner parameterization
!     Don_cape             donner_cape type derived type variable con-
!                          taining diagnostics and intermediate results
!                          related to the cape calculation associated 
!                          with the donner convection parameterization
!     Don_rad              donner_rad_type derived type variable used
!                          to hold those fields needed to connect the
!                          donner deep convection parameterization and
!                          the model radiation package
!     ermesg               character string containing any error message
!                          that is returned from a kernel subroutine
!     isize                x-direction size of the current physics window
!     isize, jsize         y-direction size of the current physics window
!     nlev_lsm             number of model layers in large-scale model
!     ntr                  number of tracers to be transported by donner
!                          convection 
!     me                   local pe number
!     calc_conv_on_this_step 
!                          is this a step on which to calculate 
!                          convection ?
!     k                    do-loop index
!
!---------------------------------------------------------------------

! mz_ht_20070421+
!!$!---------------------------------------------------------------------
!!$!    check that the module has been initialized.
!!$!---------------------------------------------------------------------
!!$      if (.not. module_is_initialized) then
!!$        call error_mesg ('donner_deep_mod', 'donner_deep: &
!!$             &donner_deep_init was not called before subroutine   &
!!$                                                  &donner_deep', FATAL)
!!$      endif
!!$
!!$!----------------------------------------------------------------------
!!$!    determine if the arguments needed when run with the strat_cloud_mod 
!!$!    are present; set cloud_tracers_present appropriately.
!!$!----------------------------------------------------------------------
      num_cld_tracers = count( (/present(qlin), present(qiin),   &
                                 present(qain), present(delta_ql), &
                                 present(delta_qi),present(delta_qa)/) )
      if (num_cld_tracers == 0) then
        cloud_tracers_present = .false.
        qlin_arg = 0.
        qiin_arg = 0.
        qain_arg = 0.
      else if (num_cld_tracers == 6) then
        cloud_tracers_present = .true.
        qlin_arg = qlin 
        qiin_arg = qiin
        qain_arg = qain
!!$      else
!!$        call error_mesg ('donner_deep_mod','donner_deep: &
!!$                        &Either none or all of the cloud tracers '// &
!!$                         'and their tendencies must be present',FATAL)
      endif
      

!!$!--------------------------------------------------------------------
!!$!    if column diagnostics have been requested for any column, call 
!!$!    donner_column_control to define the components of the 
!!$!    donner_column_diag_type variable for the diagnostic columns in this 
!!$!    window. if column diagnostics have not been requested, the needed
!!$!    variables so indicating have already been set.
!!$!--------------------------------------------------------------------
!!$      if (Col_diag%num_diag_pts > 0) then
!!$        call donner_column_control (is, ie, js, je, Time)
!!$      endif
! mz_ht_20070421-

!-------------------------------------------------------------------
!    define the dimensions for the variables in this physics window.
!    define the pe number of the current pe.
!-------------------------------------------------------------------
      isize     = ie - is + 1
      jsize     = je - js + 1
      nlev_lsm  = size(temp,3)
      ntr       = size(tracers,4) 
! mz_ht_20070421+
!      me        = mpp_pe()
      me        = p_pe
! mz_ht_20070421-
!-----------------------------------------------------------------------
!    call the kernel subroutine don_d_donner_deep_k to obtain the
!    output fields resulting from the donner deep convection parameter-
!    ization.
!-----------------------------------------------------------------------
!!$      print*, is, ie, js, je, isize, jsize, nlev_lsm, NLEV_HIRES, ntr, me
!!$      print*, cloud_tracers_present, dt
!!$      print*, "t: ", temp
!!$      print*, "mr: ",mixing_ratio
!!$      print*, "pfull: ",pfull
!!$      print*, "phalf: ",phalf
!!$      print*, "omega: ", omega
!!$      print*, "qlin_arg: ", qlin_arg
!!$      print*, "qiin_arg: ", qiin_arg
!!$      print*, "qain_arg: ", qain_arg
!!$      print*, "land: ", land
!!$      print*, "sfc_sh_flux: ", sfc_sh_flux
!!$      print*, "sfc_vapor_flux: ", sfc_vapor_flux
!!$      print*, "tr_flux: ", tr_flux
!!$      print*, "tracers: ", tracers
!!$      print*, "cell_cld_frac: ", cell_cld_frac
!!$      print*, "cell_liq_amt: ", cell_liq_amt
!!$      print*, "cell_liq_size: ", cell_liq_size
!!$      print*, "cell_ice_amt: ", cell_ice_amt
!!$      print*, "cell_ice_size: ", cell_ice_size
!!$      print*, "meso_cld_frac: ", meso_cld_frac
!!$      print*, "meso_liq_amt: ", meso_liq_amt
!!$      print*, "meso_liq_size: ", meso_liq_size
!!$      print*, "meso_ice_amt: ", meso_ice_amt
!!$      print*, "meso_ice_size: ", meso_ice_size
!!$      print*, "nsum: ", nsum
!!$      print*, "examples from Param", PARAM%grav, Param%hls
!!$      print*, "examples from NML",  Nml%parcel_launch_level, &
!!$                                    Nml%allow_mesoscale_circulation
!!$      print*, "examples from Initialized", Initialized%do_donner_tracer

      call don_d_donner_deep_k(   &
            is, ie, js, je, isize, jsize, nlev_lsm, NLEV_HIRES, ntr, me,&
            cloud_tracers_present,    &
            dt, Param, Nml, temp, mixing_ratio, pfull, phalf, omega,   &
!           qlin, qiin, qain, land, sfc_sh_flux, sfc_vapor_flux,    &
            qlin_arg, qiin_arg, qain_arg, land, sfc_sh_flux,  &
            sfc_vapor_flux,    &
            tr_flux, tracers, cell_cld_frac, cell_liq_amt,      &
            cell_liq_size, cell_ice_amt, cell_ice_size, meso_cld_frac,  &
            meso_liq_amt, meso_liq_size, meso_ice_amt, meso_ice_size,  &
            nsum, precip, delta_temp, delta_vapor, detf, uceml_inter,  &
            mtot, donner_humidity_area, donner_humidity_ratio, &
            total_precip, temperature_forcing, moisture_forcing,    &
!           parcel_rise, delta_ql, delta_qi, delta_qa, qtrtnd,         &
            parcel_rise, delta_ql_arg, delta_qi_arg, delta_qa_arg,   &
            qtrtnd,         &
            calc_conv_on_this_step, ermesg, Initialized, Col_diag,   &
            Don_rad, Don_conv, Don_cape, Don_save, exit_flag)         

!----------------------------------------------------------------------
!    if strat_cloud is active, move the output arguments into the proper
!    locations.
!----------------------------------------------------------------------
      if (cloud_tracers_present) then
        delta_ql = delta_ql_arg
        delta_qi = delta_qi_arg
        delta_qa = delta_qa_arg
      endif
!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, process the error message.
!!!  HOW TO DISTINGUISH FATAL, WARNING, NOTE ??
!    FOR NOW, ALL messages considered FATAL.
!----------------------------------------------------------------------
! mz_ht_20070421+
!!$      if (trim(ermesg) /= ' ') then
!!$        print *, 'ermesg', ermesg, me
!!$        call error_mesg ('donner_deep_mod', ermesg, FATAL)
!!$      endif
! mz_ht_20070421-
!---------------------------------------------------------------------
!    if this is a calculation step for donner_deep, define a mass
!    weighting factor (mass per unit area) needed for some of the netcdf
!    diagnostics (pmass). call donner_deep_netcdf to send the requested 
!    diagnostic data to the diag_manager for output.
!---------------------------------------------------------------------
      if (calc_conv_on_this_step) then
        do k=1,nlev_lsm
          pmass(:,:,k) = (phalf(:,:,k+1) - phalf(:,:,k))/Param%GRAV   
        end do
!!$        call donner_deep_netcdf (is, ie, js, je, Time, Don_conv,  &
!!$                                 Don_cape, parcel_rise, pmass, &
!!$                                 total_precip, temperature_forcing, &
!!$                                 moisture_forcing)

!----------------------------------------------------------------------
!    on calculation steps, update the values of the cell and
!    mesoscale cloud variables to be returned to moist_processes_mod. 
!    (on non-calculation steps, the values that were passed in are 
!    simply passed back.)
!----------------------------------------------------------------------
        cell_cld_frac = Don_rad%cell_cloud_frac
        cell_liq_amt  = Don_rad%cell_liquid_amt
        cell_liq_size = Don_rad%cell_liquid_size
        cell_ice_amt  = Don_rad%cell_ice_amt
        cell_ice_size = Don_rad%cell_ice_size
        meso_cld_frac = Don_rad%meso_cloud_frac
        meso_liq_amt  = Don_rad%meso_liquid_amt
        meso_liq_size = Don_rad%meso_liquid_size
        meso_ice_amt  = Don_rad%meso_ice_amt
        meso_ice_size = Don_rad%meso_ice_size
        nsum          = Don_rad%nsum

!--------------------------------------------------------------------
!    call deallocate_local_variables to deallocate space used by the
!    local derived-type variables.
!--------------------------------------------------------------------
        call don_d_dealloc_loc_vars_k   &
               (Don_conv, Don_cape, Don_rad, ermesg)

!----------------------------------------------------------------------
!    determine if an error message was returned from the kernel routine.
!    if so, process the error message.
!----------------------------------------------------------------------
! mz_ht_20070421
!!$        if (trim(ermesg) /= ' ') then
!!$          call error_mesg ('donner_deep_mod', ermesg, FATAL)
!!$        endif
! mz_ht_20070421
      endif  ! (calc_conv_on_this_step)

!--------------------------------------------------------------------

END SUBROUTINE DONNER_CUMASTR

!================================================================

END MODULE MESSY_CONVECT








