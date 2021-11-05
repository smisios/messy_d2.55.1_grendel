MODULE MESSY_CONVECT_BECHTOLD

USE MESSY_CONVECT_BECHTOLD_PARAM, ONLY: dp,                     &
                                        RG, RV, RD, RCPD, RCPV, &
                                        RALPW, RTT, RCW, RLVTT, &
                                        RCS, RLSTT,             &
                                        RATM, RGAMW, RBETW,     &
                                        JCVEXB, JCVEXT

USE MESSY_CONVECT_BECHTOLD_DEEP,  ONLY: CONVECT_TRIGGER_FUNCT, CONVECT_SATMIXRATIO, &
                                        CONVECT_UPDRAFT, CONVECT_CONDENS,           &
                                        CONVECT_MIXING_FUNCT, CONVECT_TSTEP_PREF,   &
                                        CONVECT_DOWNDRAFT, CONVECT_PRECIP_ADJUST,   &
                                        CONVECT_CLOSURE,                            &
                                        CONVECT_UV_TRANSPORT, CONVECT_CHEM_TRANSPORT

USE MESSY_CONVECT_BECHTOLD_SHAL,  ONLY: CONVECT_TRIGGER_SHAL, CONVECT_UPDRAFT_SHAL, &
                                        CONVECT_CLOSURE_SHAL
         
IMPLICIT NONE
SAVE


INTRINSIC :: COUNT, EXP, INT, LOG, MAX, MAXVAL, MIN, REAL, PACK, SPREAD, UNPACK

CONTAINS
!############################################################################
  SUBROUTINE CONVECT_DEEP( KLON, KLEV, KIDIA, KFDIA, KBDIA, KTDIA,          &
                           & PDTCONV, KICE, OREFRESH, ODOWN, OSETTADJ,      &
                           & PPABST, PPAH, PZZ, PDXDY, PTIMEC, PHSFLX,      &
                           & PTT, PRVT, PRCT, PRIT, PUT, PVT, PWT,          &
                           & KCOUNT, PTTEN, PRVTEN, PRCTEN, PRITEN,         &
                           & PPRLTEN, PPRSTEN,                              &
                           & KCLTOP, KCLBAS, PPRLFLX, PPRSFLX,              &
                           & PUMF, PDMF, PURV, PURCI, PCAPE,                &
                           & OUVCONV, PUTEN, PVTEN,                         &
                           & OCH1CONV, KCH1, PCH1, PCH1TEN,                 &
                           ! for ERA40/MESSy
                           & PUDR, PDDR, PUER, PDER, WAT_DIAG,              &
                           & PURCICE, PURCLIQ, PURCRPRO, PURCSPRO ) 
!############################################################################

!!**** Monitor routine to compute all convective tendencies by calls
!!     of several subroutines.
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine the convective
!!      tendencies. The routine first prepares all necessary grid-scale
!!      variables. The final convective tendencies are then computed by
!!      calls of different subroutines.
!!
!!
!!**  METHOD
!!    ------
!!      We start by selecting convective columns in the model domain through
!!      the call of routine TRIGGER_FUNCT. Then, we allocate memory for the
!!      convection updraft and downdraft variables and gather the grid scale
!!      variables in convective arrays.
!!      The updraft and downdraft computations are done level by level starting
!!      at the  bottom and top of the domain, respectively.
!!      All computations are done on MNH thermodynamic levels. The depth
!!      of the current model layer k is defined by DP(k)=P(k-1)-P(k)
!!
!!
!!
!!    EXTERNAL
!!    --------
!!    CONVECT_TRIGGER_FUNCT
!!    CONVECT_SATMIXRATIO
!!    CONVECT_UPDRAFT
!!        CONVECT_CONDENS
!!        CONVECT_MIXING_FUNCT
!!    CONVECT_TSTEP_PREF
!!    CONVECT_DOWNDRAFT
!!    CONVECT_PRECIP_ADJUST
!!    CONVECT_CLOSURE
!!        CONVECT_CLOSURE_THRVLCL
!!        CONVECT_CLOSURE_ADJUST
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module YOMCST
!!          RG                   ! gravity constant
!!          RPI                  ! number Pi
!!          RATM                 ! reference pressure
!!          RD, RV               ! gaz  constants for dry air and water vapor
!!          RCPD, RCPV           ! specific heat for dry air and water vapor
!!          RALPW, RBETW, RGAMW  ! constants for water saturation pressure
!!          RTT                  ! triple point temperature
!!          RLVTT, RLSTT         ! vaporization, sublimation heat constant
!!          RCW, RCS             ! specific heat for liquid water and ice
!!          RHOH2O               ! density of liquid water (normally defined in Module YOETHF )
!!                               ! but redefine dhere as RATM/100
!!
!!      Module YOE_CONVPAREXT
!!          JCVEXB, JCVEXT       ! extra levels on the vertical boundaries
!!
!!      Module YOE_CONVPAR
!!          XA25                 ! reference grid area
!!          XCRAD                ! cloud radius
!!
!!
!!
!!    REFERENCE
!!    ---------
!!
!!      Bechtold et al., 2001, Quart. J. Roy. Meteor. Soc. :
!!           A mass flux convection scheme for regional and global models.
!!      Kain and Fritsch, 1990, J. Atmos. Sci., Vol. 47, 2784-2801.
!!      Kain and Fritsch, 1993, Meteor. Monographs, Vol. 24, 165-170.
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96
!!   Peter Bechtold 04/10/97 replace theta_il by enthalpy
!!         "        10/12/98 changes for ARPEGE
!!
!!     Holger Tost  13/08/04 changes for MESSy
!-------------------------------------------------------------------------------



!*       0.    DECLARATIONS
!              ------------

USE messy_convect_bechtold_param,    ONLY: XA25, XCRAD,            &
                                           RPI
IMPLICIT NONE

!*       0.1   Declarations of dummy arguments :


INTEGER,                        INTENT(IN) :: KLON     ! horizontal dimension
INTEGER,                        INTENT(IN) :: KLEV     ! vertical dimension
INTEGER,                        INTENT(IN) :: KIDIA    ! value of the first point in x
INTEGER,                        INTENT(IN) :: KFDIA    ! value of the last point in x
INTEGER,                        INTENT(IN) :: KBDIA    ! vertical  computations start at
!                                                      ! KBDIA that is at least 1
INTEGER,                        INTENT(IN) :: KTDIA    ! vertical computations can be
                                                       ! limited to KLEV + 1 - KTDIA
                                                       ! default=1
REAL(dp),                       INTENT(IN) :: PDTCONV  ! Interval of time between two
                                                       ! calls of the deep convection
                                                       ! scheme
INTEGER,                        INTENT(IN) :: KICE     ! flag for ice ( 1 = yes,
                                                        !                0 = no ice )
LOGICAL,                        INTENT(IN) :: OREFRESH ! refresh or not tendencies
                                                       ! at every call
LOGICAL,                        INTENT(IN) :: ODOWN    ! take or not convective
                                                        ! downdrafts into account
LOGICAL,                        INTENT(IN) :: OSETTADJ ! logical to set convective
                                                        ! adjustment time by user
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PTT      ! grid scale temperature (K)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PRVT     ! grid scale water vapor(kg/kg) 
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PRCT     ! grid scale r_c (kg/kg) 
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PRIT     ! grid scale r_i (kg/kg) 
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PUT      ! grid scale horiz. wind u (m/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PVT      ! grid scale horiz. wind v (m/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PWT      ! grid scale vertical
                                                       ! velocity (m/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PPABST   ! grid scale pressure at (Pa)
REAL(dp), DIMENSION(KLON,0:KLEV),INTENT(IN):: PPAH     ! grid scale half-level pressure at (Pa)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PZZ      ! height of model layer (m)
REAL(dp), DIMENSION(KLON),      INTENT(IN) :: PDXDY    ! horizontal grid area (m2)
REAL(dp), DIMENSION(KLON),      INTENT(IN) :: PTIMEC   ! value of convective adjustment
                                                       ! time if OSETTADJ=.TRUE.
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PHSFLX   ! turbulent heat flux (W/m2)

INTEGER,  DIMENSION(KLON),      INTENT(INOUT):: KCOUNT ! convective counter (recompute
                                                       ! tendency or keep it)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT):: PTTEN  ! convective temperature
                                                       ! tendency (K/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT):: PRVTEN ! convective r_v tendency (1/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT):: PRCTEN ! convective r_c tendency (1/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT):: PRITEN ! convective r_i tendency (1/s)
REAL(dp), DIMENSION(KLON),      INTENT(OUT):: PPRLTEN! liquid surf. precipitation
                                                       ! tendency (m/s)
REAL(dp), DIMENSION(KLON),      INTENT(OUT):: PPRSTEN! solid surf. precipitation
                                                       ! tendency (m/s)
INTEGER,  DIMENSION(KLON),      INTENT(INOUT):: KCLTOP ! cloud top level
INTEGER,  DIMENSION(KLON),      INTENT(INOUT):: KCLBAS ! cloud base level
                                                       ! they are given a value of
                                                       ! 0 if no convection
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT):: PPRLFLX! liquid precip flux (m/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT):: PPRSFLX! solid  precip flux (m/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT):: PUMF   ! updraft mass flux (kg/s m2)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT):: PDMF   ! downdraft mass flux (kg/s m2)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT):: PURV   ! updraft water vapor (kg/kg)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT):: PURCI  ! updraft liquid+ice condensate (kg/kg)
REAL(dp), DIMENSION(KLON),      INTENT(OUT):: PCAPE  ! maximum CAPE (J/kg)

LOGICAL,                        INTENT(IN) :: OUVCONV! include wind transport (Cu friction)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT):: PUTEN  ! convective u tendency (m/s^2)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT):: PVTEN  ! convective v tendency (m/s^2)

LOGICAL,                        INTENT(IN) :: OCH1CONV ! include tracer transport
INTEGER,                        INTENT(IN) :: KCH1     ! number of species

REAL(dp), DIMENSION(KLON,KLEV,KCH1), INTENT(IN) :: PCH1    ! grid scale chemical species
REAL(dp), DIMENSION(KLON,KLEV,KCH1), INTENT(OUT):: PCH1TEN ! species conv. tendency (1/s)

! for ERA40
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT)  :: PUDR   ! updraft detrainment rate (kg/s m2)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT)  :: PDDR   ! downdraft detrainment rate (kg/s m2)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT)  :: PUER   ! updraft entrainment rate (kg/s m2)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT)  :: PDER   ! downdraft entrainment rate (kg/s m2)

REAL(DP), DIMENSION(KLON),      INTENT(INOUT):: WAT_DIAG   ! correction for column water
                                                           ! relative to total precip
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PURCICE
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PURCLIQ
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PURCRPRO
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PURCSPRO
!*       0.2   Declarations of local fixed memory variables :

INTEGER  :: ITEST, ICONV, ICONV1    ! number of convective columns
INTEGER  :: IIB, IIE                ! horizontal loop bounds
INTEGER  :: IKB, IKE                ! vertical loop bounds
INTEGER  :: IKS                     ! vertical dimension
INTEGER  :: JI, JL                  ! horizontal loop index
INTEGER  :: JN                      ! number of tracers
INTEGER  :: JK, JKP, JKM            ! vertical loop index
INTEGER  :: IFTSTEPS                ! only used for chemical tracers
REAL(dp) :: RHOH2O                  ! density of liquid water (normally defined in Module YOETHF )
REAL(dp) :: ZEPS, ZEPSA, ZEPSB      ! R_d / R_v, R_v / R_d, RCPV / RCPD - ZEPSA
REAL(dp) :: ZCPORD, ZRDOCP          ! C_p/R_d,  R_d/C_p

LOGICAL,   DIMENSION(KLON, KLEV)        :: GTRIG3 ! 3D logical mask for convection
LOGICAL,   DIMENSION(KLON)              :: GTRIG  ! 2D logical mask for trigger test
REAL(dp),  DIMENSION(KLON,KLEV)         :: ZTHT, ZSTHV, ZSTHES  ! grid scale theta,
                                                              ! theta_v, theta_es
REAL(dp),  DIMENSION(KLON)              :: ZTIME  ! convective time period
REAL(dp),  DIMENSION(KLON)              :: ZWORK2, ZWORK2B ! work array
REAL(dp)                                :: ZW1 ! work


!*       0.2   Declarations of local allocatable  variables :

INTEGER, DIMENSION(:),ALLOCATABLE  :: IDPL    ! index for parcel departure level
INTEGER, DIMENSION(:),ALLOCATABLE  :: IPBL    ! index for source layer top
INTEGER, DIMENSION(:),ALLOCATABLE  :: ILCL    ! index for lifting condensation level
INTEGER, DIMENSION(:),ALLOCATABLE  :: IETL    ! index for zero buoyancy level
INTEGER, DIMENSION(:),ALLOCATABLE  :: ICTL    ! index for cloud top level
INTEGER, DIMENSION(:),ALLOCATABLE  :: ILFS    ! index for level of free sink
INTEGER, DIMENSION(:),ALLOCATABLE  :: IDBL    ! index for downdraft base level
INTEGER, DIMENSION(:),ALLOCATABLE  :: IML     ! melting level

INTEGER, DIMENSION(:), ALLOCATABLE :: ISDPL   ! index for parcel departure level
INTEGER, DIMENSION(:), ALLOCATABLE :: ISPBL   ! index for source layer top
INTEGER, DIMENSION(:), ALLOCATABLE :: ISLCL   ! index for lifting condensation level

REAL(dp), DIMENSION(:), ALLOCATABLE    :: ZSTHLCL ! updraft theta at LCL
REAL(dp), DIMENSION(:), ALLOCATABLE    :: ZSTLCL  ! updraft temp. at LCL
REAL(dp), DIMENSION(:), ALLOCATABLE    :: ZSRVLCL ! updraft rv at LCL
REAL(dp), DIMENSION(:), ALLOCATABLE    :: ZSWLCL  ! updraft w at LCL
REAL(dp), DIMENSION(:), ALLOCATABLE    :: ZSZLCL  ! LCL height
REAL(dp), DIMENSION(:), ALLOCATABLE    :: ZSTHVELCL! envir. theta_v at LCL
REAL(dp), DIMENSION(:), ALLOCATABLE    :: ZSDXDY  ! grid area (m^2)

! grid scale variables
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZZ      ! height of model layer (m)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZPRES   ! grid scale pressure
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZPAH    ! grid scale half-level pressure at (Pa)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZDPRES  ! pressure difference between
                                                  ! bottom and top of layer (Pa)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZU      ! grid scale horiz. u component on theta grid
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZV      ! grid scale horiz. v component on theta grid
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZW      ! grid scale vertical velocity on theta grid
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZHSFLX  ! turbulent sensible heat flux (W/m2)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZTT     ! temperature
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZTH     ! grid scale theta
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZTHV    ! grid scale theta_v
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZTHL    ! grid scale enthalpy (J/kg)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZTHES, ZTHEST ! grid scale saturated theta_e
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZRW     ! grid scale total water (kg/kg)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZRV     ! grid scale water vapor (kg/kg)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZRC     ! grid scale cloud water (kg/kg)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZRI     ! grid scale cloud ice (kg/kg)
REAL(dp), DIMENSION(:),   ALLOCATABLE  :: ZDXDY   ! grid area (m^2)

! updraft variables
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZUMF    ! updraft mass flux (kg/s)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZUER    ! updraft entrainment (kg/s)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZUDR    ! updraft detrainment (kg/s)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZUPR    ! updraft precipitation in
                                                  ! flux units (kg water / s)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZUTHL   ! updraft enthalpy (J/kg)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZUTHV   ! updraft theta_v (K)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZURW    ! updraft total water (kg/kg)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZURC    ! updraft cloud water (kg/kg)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZURI    ! updraft cloud ice   (kg/kg)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZURR    ! liquid precipit. (kg/kg)
                                                  ! produced in  model layer
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZURS    ! solid precipit. (kg/kg)
                                                  ! produced in  model layer
REAL(dp), DIMENSION(:),   ALLOCATABLE  :: ZUTPR   ! total updraft precipitation (kg/s)
REAL(dp), DIMENSION(:),   ALLOCATABLE  :: ZMFLCL  ! cloud base unit mass flux(kg/s)
REAL(dp), DIMENSION(:),   ALLOCATABLE  :: ZCAPE   ! available potent. energy
REAL(dp), DIMENSION(:),   ALLOCATABLE  :: ZTHLCL  ! updraft theta at LCL
REAL(dp), DIMENSION(:),   ALLOCATABLE  :: ZTLCL   ! updraft temp. at LCL
REAL(dp), DIMENSION(:),   ALLOCATABLE  :: ZRVLCL  ! updraft rv at LCL
REAL(dp), DIMENSION(:),   ALLOCATABLE  :: ZWLCL   ! updraft w at LCL
REAL(dp), DIMENSION(:),   ALLOCATABLE  :: ZZLCL   ! LCL height
REAL(dp), DIMENSION(:),   ALLOCATABLE  :: ZTHVELCL! envir. theta_v at LCL

! downdraft variables
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZDMF    ! downdraft mass flux (kg/s)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZDER    ! downdraft entrainment (kg/s)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZDDR    ! downdraft detrainment (kg/s)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZDTHL   ! downdraft enthalpy (J/kg)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZDRW    ! downdraft total water (kg/kg)
REAL(dp), DIMENSION(:),   ALLOCATABLE  :: ZMIXF   ! mixed fraction at LFS
REAL(dp), DIMENSION(:),   ALLOCATABLE  :: ZTPR    ! total surf precipitation (kg/s)
REAL(dp), DIMENSION(:),   ALLOCATABLE  :: ZSPR    ! solid surf precipitation (kg/s)
REAL(dp), DIMENSION(:),   ALLOCATABLE  :: ZDTEVR  ! donwndraft evapor. (kg/s)
REAL(dp), DIMENSION(:),   ALLOCATABLE  :: ZPREF   ! precipitation efficiency
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZDTEVRF ! donwndraft evapor. (kg/s)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZPRLFLX ! liquid precip flux
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZPRSFLX ! solid precip flux

! closure variables
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZLMASS  ! mass of model layer (kg)
REAL(dp), DIMENSION(:),   ALLOCATABLE  :: ZTIMEA  ! advective time period
REAL(dp), DIMENSION(:),   ALLOCATABLE  :: ZTIMEC, ZTIMED! time during which convection is
                                                  ! active at grid point (as ZTIME)

REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZTHC    ! conv. adj. grid scale theta
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZRVC    ! conv. adj. grid scale r_w
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZRCC    ! conv. adj. grid scale r_c
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZRIC    ! conv. adj. grid scale r_i
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZWSUB   ! envir. compensating subsidence (Pa/s)

LOGICAL,   DIMENSION(:),ALLOCATABLE  :: GTRIG1  ! logical mask for convection
LOGICAL,   DIMENSION(:),ALLOCATABLE  :: GWORK   ! logical work array
INTEGER,   DIMENSION(:),ALLOCATABLE  :: IINDEX, IJINDEX, IJSINDEX, IJPINDEX!hor.index

REAL(dp), DIMENSION(:),   ALLOCATABLE  :: ZCPH    ! specific heat C_ph
REAL(dp), DIMENSION(:),   ALLOCATABLE  :: ZLV, ZLS! latent heat of vaporis., sublim.
REAL(dp)                               :: ZES     ! saturation vapor mixng ratio

! for U, V transport:
REAL(dp),  DIMENSION(:,:), ALLOCATABLE :: ZUC     ! horizontal wind u (m/s)
REAL(dp),  DIMENSION(:,:), ALLOCATABLE :: ZVC     ! horizontal wind v (m/s)

! for Chemical Tracer transport:
REAL(dp),  DIMENSION(:,:,:), ALLOCATABLE:: ZCH1    ! grid scale chemical specy (kg/kg)
REAL(dp),  DIMENSION(:,:,:), ALLOCATABLE:: ZCH1C   ! conv. adjust. chemical specy 1
REAL(dp),  DIMENSION(:,:),   ALLOCATABLE:: ZWORK3  ! work arrays
LOGICAL,   DIMENSION(:,:,:), ALLOCATABLE:: GTRIG4  ! logical mask

! for water correction
REAL(dp), DIMENSION(:),   ALLOCATABLE :: FRAC_L    ! fraction of liquid precip per box to total
REAL(dp), DIMENSION(:),   ALLOCATABLE :: FRAC_S    ! fraction of solid precip per box to total
REAL(dp), DIMENSION(:),   ALLOCATABLE :: SUM_PREC  ! sum of all precip in this column
LOGICAL                               :: SW        ! flag for wind correction
!-------------------------------------------------------------------------------
RHOH2O = 1000._dp       ! in kg/m^3,  original: RATM *1.E-2_dp

!*       0.3    Compute loop bounds
!               -------------------

IIB    = KIDIA
IIE    = KFDIA
JCVEXB = MAX( 0, KBDIA - 1 )
IKB    = 1 + JCVEXB
IKS    = KLEV
JCVEXT = MAX( 0, KTDIA - 1 )
IKE    = IKS - JCVEXT


!*       0.5    Update convective counter ( where KCOUNT > 0
!               convection is still active ).
!               ---------------------------------------------

KCOUNT(IIB:IIE) = KCOUNT(IIB:IIE) - 1

IF ( OREFRESH ) THEN
  KCOUNT(:)       = 1
  KCOUNT(IIB:IIE) = 0 ! refresh or not at every call
ENDIF

GTRIG(:)  = KCOUNT(:) <= 0
ITEST     = COUNT( GTRIG(:) )

IF ( ITEST == 0 ) RETURN  ! if convection is already active at every grid point
                          ! exit CONVECT_DEEP


!*       0.7    Reset convective tendencies to zero if convective
!               counter becomes negative
!               -------------------------------------------------

GTRIG3(:,:) = SPREAD( GTRIG(:), DIM=2, NCOPIES=IKS )
WHERE ( GTRIG3(:,:) )
    PTTEN(:,:)    = 0.0_dp
    PRVTEN(:,:)   = 0.0_dp
    PRCTEN(:,:)   = 0.0_dp
    PRITEN(:,:)   = 0.0_dp
    PPRLFLX(:,:)  = 0.0_dp
    PPRSFLX(:,:)  = 0.0_dp
    PUTEN(:,:)    = 0.0_dp
    PVTEN(:,:)    = 0.0_dp
    PUMF(:,:)     = 0.0_dp
    PDMF(:,:)     = 0.0_dp
    PURV(:,:)     = 0.0_dp
    PURCI(:,:)    = 0.0_dp
    PUDR(:,:)     = 0.0_dp
    PDDR(:,:)     = 0.0_dp
    PUER(:,:)     = 0.0_dp
    PDER(:,:)     = 0.0_dp
    PURCICE(:,:)  = 0.0_dp
    PURCLIQ(:,:)  = 0.0_dp
    PURCRPRO(:,:) = 0.0_dp
    PURCSPRO(:,:) = 0.0_dp
END WHERE

WHERE ( GTRIG(:) )
   PPRLTEN(:) = 0.0_dp
   PPRSTEN(:) = 0.0_dp
!  KCLTOP(:)  = 1 ! already initialized in CONVECTION
!  KCLBAS(:)  = 1
   PCAPE(:)   = 0.0_dp
END WHERE

IF ( OCH1CONV ) THEN
   ALLOCATE( GTRIG4(KLON,KLEV,KCH1) )
   GTRIG4(:,:,:) = SPREAD( GTRIG3(:,:), DIM=3, NCOPIES=KCH1 )
   WHERE( GTRIG4(:,:,:) ) PCH1TEN(:,:,:) = 0.0_dp
   DEALLOCATE( GTRIG4 )
ENDIF


!*       1.     Initialize  local variables
!               ----------------------------

ZEPS   = RD / RV
ZEPSA  = RV / RD
ZEPSB  = RCPV / RCPD - ZEPSA
ZCPORD = RCPD / RD
ZRDOCP = RD / RCPD


!*       1.1    Set up grid scale theta, theta_v, theta_es
!               ------------------------------------------

ZTHT(:,:)   = 300._dp
ZSTHV(:,:)  = 300._dp
ZSTHES(:,:) = 400._dp

DO JK = IKB, IKE
DO JI = IIB, IIE
   IF ( PPABST(JI,JK) > 40.E2_dp ) THEN
      ZTHT(JI,JK)  = PTT(JI,JK) * ( RATM / PPABST(JI,JK) ) ** ZRDOCP
      ZSTHV(JI,JK) = ZTHT(JI,JK) * ( 1.0_dp + ZEPSA * PRVT(JI,JK) ) /           &
                   & ( 1.0_dp + PRVT(JI,JK) + PRCT(JI,JK) + PRIT(JI,JK) )

          ! use conservative Bolton (1980) formula for theta_e
          ! it is used to compute CAPE for undilute parcel ascent
          ! For economical reasons we do not use routine CONVECT_SATMIXRATIO here

      ZES = EXP( RALPW - RBETW / PTT(JI,JK) - RGAMW * LOG( PTT(JI,JK) ) )
      ZES = MIN( 1.0_dp, ZEPS * ZES / ( PPABST(JI,JK) - ZES ) )
      ZSTHES(JI,JK) = PTT(JI,JK) * ( ZTHT(JI,JK) / PTT(JI,JK) ) **         &
                    & ( 1.0_dp - 0.28_dp * ZES ) *                        &
                    &   EXP( ( 3374.6525_dp / PTT(JI,JK) - 2.5403_dp ) &
                    &   * ZES * ( 1.0_dp + 0.81_dp * ZES ) )
   ENDIF
ENDDO
ENDDO



!*       2.     Test for convective columns and determine properties at the LCL
!               --------------------------------------------------------------

!*       2.1    Allocate arrays depending on number of model columns that need
!               to be tested for convection (i.e. where no convection is present
!               at the moment.
!               --------------------------------------------------------------

     ALLOCATE( ZPRES(ITEST,IKS) )
     ALLOCATE( ZZ(ITEST,IKS) )
     ALLOCATE( ZW(ITEST,IKS) )
     ALLOCATE( ZTH(ITEST,IKS) )
     ALLOCATE( ZTHV(ITEST,IKS) )
     ALLOCATE( ZTHEST(ITEST,IKS) )
     ALLOCATE( ZRV(ITEST,IKS) )
     ALLOCATE( ZSTHLCL(ITEST) )
     ALLOCATE( ZSTLCL(ITEST) )
     ALLOCATE( ZSRVLCL(ITEST) )
     ALLOCATE( ZSWLCL(ITEST) )
     ALLOCATE( ZSZLCL(ITEST) )
     ALLOCATE( ZSTHVELCL(ITEST) )
     ALLOCATE( ISDPL(ITEST) )
     ALLOCATE( ISPBL(ITEST) )
     ALLOCATE( ISLCL(ITEST) )
     ALLOCATE( ZSDXDY(ITEST) )
     ALLOCATE( GTRIG1(ITEST) )
     ALLOCATE( ZCAPE(ITEST) )
     ALLOCATE( IINDEX(KLON) )
     ALLOCATE( IJSINDEX(ITEST) )
     ALLOCATE( ZHSFLX(ITEST,IKS) )
     DO JI = 1, KLON
        IINDEX(JI) = JI
     ENDDO
     IJSINDEX(:) = PACK( IINDEX(:), MASK=GTRIG(:) )

  DO JK = IKB, IKE
  DO JI = 1, ITEST
     JL = IJSINDEX(JI)
     ZPRES(JI,JK)  = PPABST(JL,JK)
     ZZ(JI,JK)     = PZZ(JL,JK)
     ZTH(JI,JK)    = ZTHT(JL,JK)
     ZTHV(JI,JK)   = ZSTHV(JL,JK)
     ZTHEST(JI,JK) = ZSTHES(JL,JK)
     ZRV(JI,JK)    = MAX( 0.0_dp, PRVT(JL,JK) )
     ZW(JI,JK)     = PWT(JL,JK)
     ZHSFLX(JI,JK) = PHSFLX(JL,JK)
  ENDDO
  ENDDO
  DO JI = 1, ITEST
     JL = IJSINDEX(JI)
     ZSDXDY(JI)    = PDXDY(JL)
  ENDDO

!*       2.2    Compute environm. enthalpy and total water = r_v + r_i + r_c
!               and envir. saturation theta_e
!               ------------------------------------------------------------


!*       2.3    Test for convective columns and determine properties at the LCL
!               --------------------------------------------------------------

     ISLCL(:) = MAX( IKB, 2 )   ! initialize DPL PBL and LCL
     ISDPL(:) = IKB
     ISPBL(:) = IKB


     CALL CONVECT_TRIGGER_FUNCT( ITEST, KLEV,                              &
                               & ZPRES, ZTH, ZTHV, ZTHEST,                 &
                               & ZRV, ZW, ZZ, ZSDXDY, ZHSFLX,              &
                               & ZSTHLCL, ZSTLCL, ZSRVLCL, ZSWLCL, ZSZLCL, &
                               & ZSTHVELCL, ISLCL, ISDPL, ISPBL, GTRIG1,   &
                               & ZCAPE )

     DO JI = 1, ITEST
        JL = IJSINDEX(JI)
        PCAPE(JL) = ZCAPE(JI)
     ENDDO

     DEALLOCATE( ZPRES )
     DEALLOCATE( ZZ )
     DEALLOCATE( ZTH )
     DEALLOCATE( ZTHV )
     DEALLOCATE( ZTHEST )
     DEALLOCATE( ZRV )
     DEALLOCATE( ZW )
     DEALLOCATE( ZCAPE )
     DEALLOCATE( ZHSFLX )


!*       3.     After the call of TRIGGER_FUNCT we allocate all the dynamic
!               arrays used in the convection scheme using the mask GTRIG, i.e.
!               we do calculus only in convective columns. This corresponds to
!               a GATHER operation.
!               --------------------------------------------------------------

     ICONV = COUNT( GTRIG1(:) )
     IF ( ICONV == 0 )  THEN
         DEALLOCATE( ZSTHLCL )
         DEALLOCATE( ZSTLCL )
         DEALLOCATE( ZSRVLCL )
         DEALLOCATE( ZSWLCL )
         DEALLOCATE( ZSZLCL )
         DEALLOCATE( ZSTHVELCL )
         DEALLOCATE( ZSDXDY )
         DEALLOCATE( ISLCL )
         DEALLOCATE( ISDPL )
         DEALLOCATE( ISPBL )
         DEALLOCATE( GTRIG1 )
         DEALLOCATE( IINDEX )
         DEALLOCATE( IJSINDEX )
         RETURN   ! no convective column has been found, exit CONVECT_DEEP
     ENDIF

     ! vertical index variables

         ALLOCATE( IDPL(ICONV) )
         ALLOCATE( IPBL(ICONV) )
         ALLOCATE( ILCL(ICONV) )
         ALLOCATE( ICTL(ICONV) )
         ALLOCATE( IETL(ICONV) )

         ! grid scale variables

         ALLOCATE( ZZ(ICONV,IKS) )
         ALLOCATE( ZPRES(ICONV,IKS) )
         ALLOCATE( ZPAH(ICONV,0:IKS) )
         ALLOCATE( ZDPRES(ICONV,IKS) )
         ALLOCATE( ZU(ICONV,IKS) )
         ALLOCATE( ZV(ICONV,IKS) )
         ALLOCATE( ZTT(ICONV, IKS) )
         ALLOCATE( ZTH(ICONV,IKS) )
         ALLOCATE( ZTHV(ICONV,IKS) )
         ALLOCATE( ZTHL(ICONV,IKS) )
         ALLOCATE( ZTHES(ICONV,IKS) )
         ALLOCATE( ZRV(ICONV,IKS) )
         ALLOCATE( ZRC(ICONV,IKS) )
         ALLOCATE( ZRI(ICONV,IKS) )
         ALLOCATE( ZRW(ICONV,IKS) )
         ALLOCATE( ZDXDY(ICONV) )

         ! updraft variables

         ALLOCATE( ZUMF(ICONV,IKS) )
         ALLOCATE( ZUER(ICONV,IKS) )
         ALLOCATE( ZUDR(ICONV,IKS) )
         ALLOCATE( ZUPR(ICONV,IKS) )
         ALLOCATE( ZUTHL(ICONV,IKS) )
         ALLOCATE( ZUTHV(ICONV,IKS) )
         ALLOCATE( ZURW(ICONV,IKS) )
         ALLOCATE( ZURC(ICONV,IKS) )
         ALLOCATE( ZURI(ICONV,IKS) )
         ALLOCATE( ZURR(ICONV,IKS) )
         ALLOCATE( ZURS(ICONV,IKS) )
         ALLOCATE( ZUTPR(ICONV) )
         ALLOCATE( ZTHLCL(ICONV) )
         ALLOCATE( ZTLCL(ICONV) )
         ALLOCATE( ZRVLCL(ICONV) )
         ALLOCATE( ZWLCL(ICONV) )
         ALLOCATE( ZMFLCL(ICONV) )
         ALLOCATE( ZZLCL(ICONV) )
         ALLOCATE( ZTHVELCL(ICONV) )
         ALLOCATE( ZCAPE(ICONV) )

         ! work variables

         ALLOCATE( IJINDEX(ICONV) )
         ALLOCATE( IJPINDEX(ICONV) )
         ALLOCATE( ZCPH(ICONV) )
         ALLOCATE( ZLV(ICONV) )
         ALLOCATE( ZLS(ICONV) )


!*           3.1    Gather grid scale and updraft base variables in
!                   arrays using mask GTRIG
!                   ---------------------------------------------------

         GTRIG(:)      = UNPACK( GTRIG1(:), MASK=GTRIG(:), FIELD=.FALSE. )
         IJINDEX(:)    = PACK( IINDEX(:), MASK=GTRIG(:) )

    DO JK = IKB, IKE
    DO JI = 1, ICONV
         JL = IJINDEX(JI)
         ZZ(JI,JK)     = PZZ(JL,JK)
         ZPRES(JI,JK)  = PPABST(JL,JK)
         ZPAH(JI,JK)   = PPAH(JL,JK)
         ZTT(JI,JK)    = PTT(JL,JK)
         ZTH(JI,JK)    = ZTHT(JL,JK)
         ZTHES(JI,JK)  = ZSTHES(JL,JK)
         ZRV(JI,JK)    = MAX( 0.0_dp, PRVT(JL,JK) )
         ZRC(JI,JK)    = MAX( 0.0_dp, PRCT(JL,JK) )
         ZRI(JI,JK)    = MAX( 0.0_dp, PRIT(JL,JK) )
         ZTHV(JI,JK)   = ZSTHV(JL,JK)
         ZU(JI,JK)     = PUT(JL,JK)
         ZV(JI,JK)     = PVT(JL,JK)
    ENDDO
    ENDDO
    ! mz_ht_20050913+   
    DO  JI = 1, ICONV
         JL = IJINDEX(JI)
         ZPAH(JI,0)    = PPAH(JL,0)
    ENDDO
    ! mz_ht_20050913-

    IF ( OSETTADJ ) THEN
         ALLOCATE( ZTIMED(ICONV) )
         DO JI = 1, ICONV
            JL = IJINDEX(JI)
            ZTIMED(JI) = PTIMEC(JL)
         ENDDO
    ENDIF

    DO JI = 1, ITEST
       IJSINDEX(JI) = JI
    ENDDO
    IJPINDEX(:) = PACK( IJSINDEX(:), MASK=GTRIG1(:) )
    DO JI = 1, ICONV
         JL = IJPINDEX(JI)
         IDPL(JI)      = ISDPL(JL)
         IPBL(JI)      = ISPBL(JL)
         ILCL(JI)      = ISLCL(JL)
         ZTHLCL(JI)    = ZSTHLCL(JL)
         ZTLCL(JI)     = ZSTLCL(JL)
         ZRVLCL(JI)    = ZSRVLCL(JL)
         ZWLCL(JI)     = ZSWLCL(JL)
         ZZLCL(JI)     = ZSZLCL(JL)
         ZTHVELCL(JI)  = ZSTHVELCL(JL)
         ZDXDY(JI)     = ZSDXDY(JL)
    ENDDO
         ALLOCATE( GWORK(ICONV) )
         GWORK(:)      = PACK( GTRIG1(:),  MASK=GTRIG1(:) )
         DEALLOCATE( GTRIG1 )
         ALLOCATE( GTRIG1(ICONV) )
         GTRIG1(:)     = GWORK(:)

         DEALLOCATE( GWORK )
         DEALLOCATE( IJPINDEX )
         DEALLOCATE( ISDPL )
         DEALLOCATE( ISPBL )
         DEALLOCATE( ISLCL )
         DEALLOCATE( ZSTHLCL )
         DEALLOCATE( ZSTLCL )
         DEALLOCATE( ZSRVLCL )
         DEALLOCATE( ZSWLCL )
         DEALLOCATE( ZSZLCL )
         DEALLOCATE( ZSTHVELCL )
         DEALLOCATE( ZSDXDY )


!*           3.2    Compute pressure difference
!                   ---------------------------------------------------

        ! mz_ht_20050913+       
        ZDPRES(:,IKB) = 0.0_dp
        DO JK = IKB + 1, IKE
            ZDPRES(:,JK)  = ZPRES(:,JK-1) - ZPRES(:,JK) 
        ENDDO
!!$         ZDPRES(:,:) = 0.0_dp
!!$         DO JK = IKB, IKE
!!$           ZDPRES(:,jk) = ZPAH(:,jk-1) - ZPAH(:,jk)
!!$         END DO


        ! mz_ht_20050913-

!*           3.3   Compute environm. enthalpy and total water = r_v + r_i + r_c
!                  ----------------------------------------------------------

        DO JK = IKB, IKE, 1
            ZRW(:,JK)  = ZRV(:,JK) + ZRC(:,JK) + ZRI(:,JK)
            ZCPH(:)    = RCPD + RCPV * ZRW(:,JK)
            ZLV(:)     = RLVTT + ( RCPV - RCW ) * ( ZTT(:,JK) - RTT ) ! compute L_v
            ZLS(:)     = RLSTT + ( RCPV - RCS ) * ( ZTT(:,JK) - RTT ) ! compute L_i
            ZTHL(:,JK) = ZCPH(:) * ZTT(:,JK) + ( 1.0_dp + ZRW(:,JK) ) * RG * ZZ(:,JK) &
                       & - ZLV(:) * ZRC(:,JK) - ZLS(:) * ZRI(:,JK)
        ENDDO


!*           4.     Compute updraft properties
!                   ----------------------------

!*           4.1    Set mass flux at LCL ( here a unit mass flux with w = 1 m/s )
!                   -------------------------------------------------------------

         DO JI = 1, ICONV
               JK = ILCL(JI) - 1
               ZMFLCL(JI) = ZPRES(JI,JK) / ( RD * ZTT(JI,JK) *                    &
                          & ( 1.0_dp + ZEPS * ZRVLCL(JI) ) ) * RPI * XCRAD * XCRAD &
                          & * MAX( 1.0_dp, ZDXDY(JI)/XA25 )
         ENDDO

         DEALLOCATE( ZCPH )
         DEALLOCATE( ZLV )
         DEALLOCATE( ZLS )


     CALL CONVECT_UPDRAFT( ICONV, KLEV,                                     &
                         & KICE, ZPRES, ZDPRES, ZZ, ZTHL, ZTHV, ZTHES, ZRW, &
                         & ZTHLCL, ZTLCL, ZRVLCL, ZWLCL, ZZLCL, ZTHVELCL,   &
                         & ZMFLCL, GTRIG1, ILCL, IDPL, IPBL,                &
                         & ZUMF, ZUER, ZUDR, ZUTHL, ZUTHV, ZURW,            &
                         & ZURC, ZURI, ZURR, ZURS, ZUPR,                    &
                         & ZUTPR, ZCAPE, ICTL, IETL, zpah                   )   

!*           4.2    In routine UPDRAFT GTRIG1 has been set to false when cloud
!                   thickness is smaller than 3 km
!                   -----------------------------------------------------------


     ICONV1 = COUNT(GTRIG1)

     IF ( ICONV1 > 0 )  THEN

!*       4.3    Allocate memory for downdraft variables
!               ---------------------------------------

! downdraft variables

        ALLOCATE( ILFS(ICONV) )
        ALLOCATE( IDBL(ICONV) )
        ALLOCATE( IML(ICONV) )
        ALLOCATE( ZDMF(ICONV,IKS) )
        ALLOCATE( ZDER(ICONV,IKS) )
        ALLOCATE( ZDDR(ICONV,IKS) )
        ALLOCATE( ZDTHL(ICONV,IKS) )
        ALLOCATE( ZDRW(ICONV,IKS) )
        ALLOCATE( ZLMASS(ICONV,IKS) )
        DO JK = IKB, IKE
           ZLMASS(:,JK)  = ZDXDY(:) * ZDPRES(:,JK) / RG  ! mass of model layer
        ENDDO
        ZLMASS(:,IKB) = ZLMASS(:,IKB+1)
        ALLOCATE( ZMIXF(ICONV) )
        ALLOCATE( ZTPR(ICONV) )
        ALLOCATE( ZSPR(ICONV) )
        ALLOCATE( ZDTEVR(ICONV) )
        ALLOCATE( ZPREF(ICONV) )
        ALLOCATE( ZDTEVRF(ICONV,IKS) )
        ALLOCATE( ZPRLFLX(ICONV,IKS) )
        ALLOCATE( ZPRSFLX(ICONV,IKS) )

! closure variables

        ALLOCATE( ZTIMEA(ICONV) )
        ALLOCATE( ZTIMEC(ICONV) )
        ALLOCATE( ZTHC(ICONV,IKS) )
        ALLOCATE( ZRVC(ICONV,IKS) )
        ALLOCATE( ZRCC(ICONV,IKS) )
        ALLOCATE( ZRIC(ICONV,IKS) )
        ALLOCATE( ZWSUB(ICONV,IKS) )

! correction variables
        ALLOCATE( FRAC_L(ICONV) )
        ALLOCATE( FRAC_S(ICONV) )
        ALLOCATE( SUM_PREC(ICONV) )
        FRAC_L(:)   = 0.0_dp
        FRAC_S(:)   = 0.0_dp
        SUM_PREC(:) = 0.0_dp
        

!*           5.     Compute downdraft properties
!                   ----------------------------

!*           5.1    Compute advective time period and precipitation
!                   efficiency as a function of mean ambient wind (shear)
!                   --------------------------------------------------------

        CALL CONVECT_TSTEP_PREF( ICONV, KLEV,                          &
                               & ZU, ZV, ZPRES, ZZ, ZDXDY, ILCL, ICTL, &
                               & ZTIMEA, ZPREF )

          ! exclude convective downdrafts if desired
        IF ( .NOT. ODOWN ) ZPREF(:) = 1.0_dp

          ! Compute the period during which convection is active
        ZTIMEC(:) = MAX( 1800._dp, MIN( 3600._dp, ZTIMEA(:) ) )
        ZTIMEC(:) = REAL( INT( ZTIMEC(:) / PDTCONV ),dp ) * PDTCONV
        ZTIMEC(:) = MAX( PDTCONV, ZTIMEC(:) ) ! necessary if PDTCONV > 1800
        IF ( OSETTADJ ) THEN
             ZTIMEC(:) = MAX( PDTCONV, ZTIMED(:) )
        ENDIF


!*           5.2    Compute melting level
!                   ----------------------

        IML(:) = IKB
        DO JK = IKE, IKB, -1
          WHERE( ZTT(:,JK) <= RTT )  IML(:) = JK
        ENDDO

        CALL CONVECT_DOWNDRAFT( ICONV, KLEV,                               &
                              & KICE, ZPRES, ZDPRES, ZZ, ZTH, ZTHES,       &
                              & ZRW, ZRC, ZRI,                             &
                              & ZPREF, ILCL, ICTL, IETL,                   &
                              & ZUTHL, ZURW, ZURC, ZURI,                   &
                              & ZDMF, ZDER, ZDDR, ZDTHL, ZDRW,             &
                              & ZMIXF, ZDTEVR, ILFS, IDBL, IML,            &
                              & ZDTEVRF                                    )

!*           6.     Adjust up and downdraft mass flux to be consistent
!                   with precipitation efficiency relation.
!                   ---------------------------------------------------

       CALL CONVECT_PRECIP_ADJUST( ICONV, KLEV,                              &
                                 & ZPRES,ZUMF, ZUER, ZUDR, ZUPR, ZUTPR, ZURW,&
                                 & ZDMF, ZDER, ZDDR, ZDTHL, ZDRW,            &
                                 & ZPREF, ZTPR, ZMIXF, ZDTEVR,               &
                                 & ILFS, IDBL, ILCL, ICTL, IETL,             &
                                 & ZDTEVRF                                   )


!*           7.     Determine adjusted environmental values assuming
!                   that all available buoyant energy must be removed
!                   within an advective time step ZTIMEC.
!                   ---------------------------------------------------
       CALL CONVECT_CLOSURE( ICONV, KLEV,                                &
                           & ZPRES, ZDPRES, ZZ, ZDXDY, ZLMASS,           &
                           & ZTHL, ZTH, ZRW, ZRC, ZRI, GTRIG1,           &
                           & ZTHC, ZRVC, ZRCC, ZRIC, ZWSUB,              &
                           & ILCL, IDPL, IPBL, ILFS, ICTL, IML,          &
                           & ZUMF, ZUER, ZUDR, ZUTHL, ZURW,              &
                           & ZURC, ZURI, ZUPR,                           &
                           & ZDMF, ZDER, ZDDR, ZDTHL, ZDRW,              &
                           & ZTPR, ZSPR, ZDTEVR,                         &
                           & ZCAPE, ZTIMEC,                              &
                           & IFTSTEPS,                                   &
                           & ZDTEVRF, ZPRLFLX, ZPRSFLX )



!*           8.     Determine the final grid-scale (environmental) convective
!                   tendencies and set convective counter
!                   --------------------------------------------------------


!*           8.1    Grid scale tendencies
!                   ---------------------

          ! in order to save memory, the tendencies are temporarily stored
          ! in the tables for the adjusted grid-scale values

      DO JK = IKB, IKE
         ZTHC(:,JK) = ( ZTHC(:,JK) - ZTH(:,JK) ) / ZTIMEC(:)             &
                    & * ( ZPRES(:,JK) / RATM ) ** ZRDOCP  ! change theta in temperature
         ZRVC(:,JK) = ( ZRVC(:,JK) - ZRW(:,JK) + ZRC(:,JK) + ZRI(:,JK) ) &
                    &                            / ZTIMEC(:)

         ZRCC(:,JK) = ( ZRCC(:,JK) - ZRC(:,JK) ) / ZTIMEC(:)
         ZRIC(:,JK) = ( ZRIC(:,JK) - ZRI(:,JK) ) / ZTIMEC(:)

         ZPRLFLX(:,JK) = ZPRLFLX(:,JK) / ( RHOH2O * ZDXDY(:) )
         ZPRSFLX(:,JK) = ZPRSFLX(:,JK) / ( RHOH2O * ZDXDY(:) )
      ENDDO
    
      ZPRLFLX(:,IKB) = ZPRLFLX(:,IKB+1)
      ZPRSFLX(:,IKB) = ZPRSFLX(:,IKB+1)
      ZTPR(:) = (ZPRLFLX(:,IKB+1) + ZPRSFLX(:,IKB+1))*RHOH2O                                 


!*           8.2    Apply conservation correction
!                   -----------------------------

          ! Compute vertical integrals - fluxes

       JKM = MAXVAL( ICTL(:) )
       ZWORK2(:) = 0.0_dp
       ZWORK2B(:)= 0.0_dp
       DO JK = IKB+1, JKM
         JKP = JK + 1
         DO JI = 1, ICONV
         ! ZW1 = _HALF_ *  (ZPRES(JI,JK-1) - ZPRES(JI,JKP)) / RG
           ZW1 = (ZPAH(JI,JK-1) - ZPAH(JI,JK)) / RG
           ZWORK2(JI) = ZWORK2(JI) + ( ZRVC(JI,JK) + ZRCC(JI,JK) + ZRIC(JI,JK) ) *   & ! moisture
              &                   ZW1
         ! ZWORK2B(JI) = ZWORK2B(JI) + (                                             & ! enthalpy
         !    &                        ( RCPD + RCPV * ZRW(JI,JK) )* ZTHC(JI,JK)   - &
         !    &  ( RLVTT + ( RCPV - RCW ) * ( ZTT(JI,JK) - RTT ) ) * ZRCC(JI,JK)   - &
         !    &  ( RLSTT + ( RCPV - RCW ) * ( ZTT(JI,JK) - RTT ) ) * ZRIC(JI,JK) ) * &
         !    &                   ZW1
           ZWORK2B(JI) = ZWORK2B(JI) + (                                             & ! enthalpy
              &                          RCPD * ZTHC(JI,JK)                        - &
              &  ( RLVTT + ( RCPV - RCW ) * ( ZTT(JI,JK) - RTT ) ) * ZRCC(JI,JK)   - &
              &  ( RLSTT + ( RCPV - RCW ) * ( ZTT(JI,JK) - RTT ) ) * ZRIC(JI,JK) ) * &
              &                   ZW1
         ENDDO
       ENDDO

          ! Budget error (compare integral to surface precip.)

       DO JI = 1, ICONV
         IF ( ZTPR(JI) > 0.0_dp) THEN
           JKP = ICTL(JI) 
         ! ZW1 = RG / ( ZPRES(JI,IKB) - ZPRES(JI,JKP) - &
         !       & _HALF_*(ZDPRES(JI,IKB+1) - ZDPRES(JI,JKP+1)) )
           ZW1 = RG / ( ZPAH(JI,IKB) - ZPAH(JI,JKP) )
           ZWORK2(JI) = ( ZTPR(JI)  + ZWORK2(JI) ) !* ZW1
           ZWORK2B(JI) = ( ZTPR(JI) *                                                &
           &  ( RLVTT + ( RCPV - RCW ) * ( ZTT(JI,IKB) - RTT ) ) - ZWORK2B(JI) )     &
           &                           * ZW1
         ENDIF
       ENDDO

          ! Apply uniform correction

       DO JK = JKM, IKB+1, -1
         DO JI = 1, ICONV
           IF ( ZTPR(JI) > 0.0_dp .AND. JK <= ICTL(JI) ) THEN
               ZW1 = ABS(ZRVC(JI,JK)) +  ABS(ZRCC(JI,JK)) +  ABS(ZRIC(JI,JK)) + 1.E-12_dp
               ! moisture
!               print*, zrvc(ji,jk), zrcc(ji,jk), zric(ji,jk), zw1, &
!                       ZWORK2(JI)/(ICTL(JI)-1)/((ZPAH(JI,JK-1) - ZPAH(JI,JK)) / RG) 
!!$               ZRVC(JI,JK) = ZRVC(JI,JK) - ABS(ZRVC(JI,JK))/ZW1 * &
!!$                             ZWORK2(JI)/(ICTL(JI)-1) /            &
!!$                             ((ZPAH(JI,JK-1) - ZPAH(JI,JK)) / RG)
!!$               ZRCC(JI,JK) = ZRCC(JI,JK) - ABS(ZRCC(JI,JK))/ZW1 * &
!!$                             ZWORK2(JI)/(ICTL(JI)-1) /            &
!!$                             ((ZPAH(JI,JK-1) - ZPAH(JI,JK)) / RG)
!!$               ZRIC(JI,JK) = ZRIC(JI,JK) - ABS(ZRIC(JI,JK))/ZW1 * &
!!$                             ZWORK2(JI)/(ICTL(JI)-1) /            &
!!$                             ((ZPAH(JI,JK-1) - ZPAH(JI,JK)) / RG)
              ! ZTHC(JI,JK) = ZTHC(JI,JK) + ZWORK2B(JI) / ( RCPD + RCPV * ZRW(JI,JK) )! enthalpy
                ZTHC(JI,JK) = ZTHC(JI,JK) + ZWORK2B(JI) / RCPD                        ! enthalpy
 
                ! mz_ht_20050915+
           ENDIF
         ENDDO
       ENDDO
       
       DO JI = 1, ICONV
         sum_prec(ji) = ZPRLFLX(JI,1) + ZPRSFLX(JI,1)

         IF ( SUM_PREC(JI) > 0.0_dp .AND. JK <= ICTL(JI) ) THEN

  !             Do not correct water vapour tendency, but precipitation because it is not a 
  !             prognostic quantity!
  !             Not uniform correction but weighted with fraction of snow/rain and box/column      
                ! liquid fraction
           FRAC_L(JI) = ZPRLFLX(JI,1) / SUM_PREC(JI)
                ! snow fraction
           FRAC_S(JI) = ZPRSFLX(JI,1) / SUM_PREC(JI)
                ! correction
           ZPRLFLX(JI,1) = ZPRLFLX(JI,1) + ZWORK2(JI) * FRAC_L(JI) / RHOH2O
                ! correction
           ZPRSFLX(JI,1) = ZPRSFLX(JI,1) + ZWORK2(JI) * FRAC_S(JI) / RHOH2O
         END IF
       END DO
       
       WAT_DIAG(:) = 0.0_dp
       DO JI=1,ICONV
         IF ( ZTPR(JI) > 0.0_dp ) THEN
           JL = IJINDEX(JI)
!           WAT_DIAG(JL)   = ZWORK2(JI) / ZTPR(JI)
           WAT_DIAG(JL)   = ZWORK2(JI) 
         ENDIF
       ENDDO
       
       DEALLOCATE(FRAC_L)
       DEALLOCATE(FRAC_S)
       DEALLOCATE(SUM_PREC)
          ! extend tendencies to first model level

    ! DO JI = 1, ICONV
    !    ZWORK2(JI) = ZDPRES(JI,IKB+1) + ZDPRES(JI,IKB+2)
    !    ZTHC(JI,IKB)  = ZTHC(JI,IKB+1) * ZDPRES(JI,IKB+2)/ZWORK2(JI)
    !    ZTHC(JI,IKB+1)= ZTHC(JI,IKB+1) * ZDPRES(JI,IKB+1)/ZWORK2(JI)
    !    ZRVC(JI,IKB)  = ZRVC(JI,IKB+1) * ZDPRES(JI,IKB+2)/ZWORK2(JI)
    !    ZRVC(JI,IKB+1)= ZRVC(JI,IKB+1) * ZDPRES(JI,IKB+1)/ZWORK2(JI)
    ! ENDDO

              ! execute a "scatter"= pack command to store the tendencies in
              ! the final 2D tables
      DO JK = IKB, IKE
      DO JI = 1, ICONV
         JL = IJINDEX(JI)
         PTTEN(JL,JK)   = ZTHC(JI,JK)
         PRVTEN(JL,JK)  = ZRVC(JI,JK)
         PRCTEN(JL,JK)  = ZRCC(JI,JK)
         PRITEN(JL,JK)  = ZRIC(JI,JK)

         PURV(JL,JK)    = ZURW(JI,JK) - ZURC(JI,JK) - ZURI(JI,JK)
         PURCI(JL,JK)   = ZURC(JI,JK) + ZURI(JI,JK)
         PURCICE(JL,JK) = ZURI(JI,JK)
         PURCLIQ(JL,JK) = ZURC(JI,JK)
         PURCRPRO(JL,JK)= ZURR(JI,JK)
         PURCSPRO(JL,JK)= ZURS(JI,JK)
         PPRLFLX(JL,JK) = ZPRLFLX(JI,JK)
         PPRSFLX(JL,JK) = ZPRSFLX(JI,JK)

      ENDDO
      ENDDO




!*           8.3    Convective rainfall tendency
!                   ----------------------------

     DO JI = 1, ICONV
        JL = IJINDEX(JI)
        PPRLTEN(JL) = PPRLFLX(JL,IKB+1)
        PPRSTEN(JL) = PPRSFLX(JL,IKB+1)
     ENDDO


!                   Cloud base and top levels
!                   -------------------------

     ILCL(:) = MIN( ILCL(:), ICTL(:) )
     DO JI = 1, ICONV
        JL = IJINDEX(JI)
        KCLTOP(JL) = ICTL(JI)
        KCLBAS(JL) = ILCL(JI)
     ENDDO


!*           8.4    Set convective counter
!                   ----------------------

         ! compute convective counter for just activated convective
         ! grid points
         ! If the advective time period is less than specified
         ! minimum for convective period, allow feedback to occur only
         ! during advective time

     ZTIME(:) = 1.0_dp
     ZWORK2(:) = 0.0_dp
     DO JI = 1, ICONV
       JL = IJINDEX(JI)
       ZTIME(JL)  =  ZTIMEC(JI)
       ZWORK2(JL) =  ZTIMEA(JI)
       ZWORK2(JL) =  MIN( ZWORK2(JL), ZTIME(JL) )
       ZWORK2(JL) =  MAX( ZWORK2(JL), PDTCONV )
       IF ( GTRIG(JL) )  KCOUNT(JL) = INT( ZWORK2(JL) / PDTCONV )
       IF ( GTRIG(JL) .AND. PPRLTEN(JL)<1.E-14_dp ) KCOUNT(JL) = 0
     ENDDO



!*           8.6    Compute convective tendencies for Tracers
!                   ------------------------------------------

  IF ( OCH1CONV ) THEN

       ALLOCATE( ZCH1(ICONV,IKS,KCH1) )
       ALLOCATE( ZCH1C(ICONV,IKS,KCH1) )
       ALLOCATE( ZWORK3(ICONV,KCH1) )

       DO JK = IKB, IKE
       DO JI = 1, ICONV
          JL = IJINDEX(JI)
          ZCH1(JI,JK,:) = PCH1(JL,JK,:)
       ENDDO
       ENDDO

      CALL CONVECT_CHEM_TRANSPORT( ICONV, KLEV, KCH1, ZCH1, ZCH1C,      &
                                 & IDPL, IPBL, ILCL, ICTL, ILFS, IDBL,  &
                                 & ZUMF, ZUER, ZUDR, ZDMF, ZDER, ZDDR,  &
                                 & ZTIMEC, ZDXDY, ZMIXF, ZLMASS, ZWSUB, &
                                 & IFTSTEPS )

       DO JK = IKB, IKE
       DO JN = 1, KCH1
          ZCH1C(:,JK,JN) = ( ZCH1C(:,JK,JN)- ZCH1(:,JK,JN) ) / ZTIMEC(:)
       ENDDO
       ENDDO


!*           8.7    Apply conservation correction
!                   -----------------------------

          ! Compute vertical integrals

       JKM = MAXVAL( ICTL(:) )
       ZWORK3(:,:) = 0.0_dp
       DO JK = IKB+1, JKM
         JKP = JK + 1
         DO JI = 1, ICONV
           ZW1 = (ZPAH(JI,JK-1) - ZPAH(JI,JK)) / RG
           ZWORK3(JI,:) = ZWORK3(JI,:) + ZCH1C(JI,JK,:) * ZW1 
         ENDDO
       ENDDO

          ! Mass error (integral must be zero)

       DO JI = 1, ICONV
         IF ( ZTPR(JI) > 0.0_dp) THEN
           JKP = ICTL(JI) 
         ! ZW1 = RG / ( ZPRES(JI,IKB) - ZPRES(JI,JKP) - &
         !      & _HALF_*(ZDPRES(JI,IKB+1) - ZDPRES(JI,JKP+1)) )
           ZW1 = RG / ( ZPAH(JI,IKB) - ZPAH(JI,JKP) )
           ZWORK3(JI,:) = ZWORK3(JI,:) * ZW1
         ENDIF
       ENDDO

          ! Apply uniform correction but assure positive mass at each level

       DO JK = JKM, IKB+1, -1
         DO JI = 1, ICONV
           IF ( ZTPR(JI) > 0.0_dp .AND. JK <= ICTL(JI) ) THEN
                ZCH1C(JI,JK,:) = ZCH1C(JI,JK,:) - ZWORK3(JI,:)
              ! ZCH1C(JI,JK,:) = MAX( ZCH1C(JI,JK,:), -ZCH1(JI,JK,:)/ZTIMEC(JI) )
           ENDIF
         ENDDO
       ENDDO
!
          ! extend tendencies to first model level

   !   DO JI = 1, ICONV
   !      ZWORK2(JI) = ZDPRES(JI,IKB+1) + ZDPRES(JI,IKB+2)
   !   ENDDO
   !   DO JN = 1, KCH1
   !   DO JI = 1, ICONV
   !     ZCH1(JI,IKB,JN)  = ZCH1(JI,IKB+1,JN) * ZDPRES(JI,IKB+2)/ZWORK2(JI)
   !     ZCH1(JI,IKB+1,JN)= ZCH1(JI,IKB+1,JN) * ZDPRES(JI,IKB+1)/ZWORK2(JI)
   !   ENDDO
   !   ENDDO

       DO JK = IKB, IKE
       DO JI = 1, ICONV
          JL = IJINDEX(JI)
          PCH1TEN(JL,JK,:) = ZCH1C(JI,JK,:)
       ENDDO
       ENDDO
  ENDIF
!
!*           8.8    Compute convective tendencies for wind
!                   --------------------------------------

  IF ( OUVCONV ) THEN

       ALLOCATE( ZUC(ICONV,IKS) )
       ALLOCATE( ZVC(ICONV,IKS) )

       CALL CONVECT_UV_TRANSPORT( ICONV, KLEV, ZU, ZV, ZUC, ZVC,        &
                                & IDPL, IPBL, ILCL, ICTL, ILFS, IDBL,   &
                                & ZUMF, ZUER, ZUDR, ZDMF, ZDER, ZDDR,   &
                                & ZTIMEC, ZDXDY, ZMIXF, ZLMASS, ZWSUB,  &
                                & IFTSTEPS )

       DO JK = IKB, IKE
!         do ji=1,iconv
!            if ( abs (ZUC(ji,jk)) > 100._dp) print*, "WARNING high windspeed: ",zuc(ji,jk), zu(ji,jk) 
!         enddo
          ZUC(:,JK) = ( ZUC(:,JK)- ZU(:,JK) ) / ZTIMEC(:)
          ZVC(:,JK) = ( ZVC(:,JK)- ZV(:,JK) ) / ZTIMEC(:)
       ENDDO


       do ji=1,iconv
          SW = .FALSE.
          do jk=ikb,ike
             IF ( zumf(ji,jk)*PDTCONV > 2._dp*zlmass(ji,jk) ) THEN
                print*, ji,jk, zumf(ji,jk), zlmass(ji,jk), &
                     zumf(ji,jk)*pdtconv, pdtconv, zuc(ji,jk), zu(ji,jk),&
                     zuc(ji,jk)*pdtconv, zvc(ji,jk), zv(ji,jk), zvc(ji,jk)*pdtconv 
                SW = .TRUE.
             END IF
          enddo
          IF (SW) THEN
          zuc(ji,:) = 0.0_dp
          zvc(ji,:) = 0.0_dp
             WRITE(*,*) ji, "Wind correction switched of because of too high massfluxes! "
          ENDIF
       enddo
!*           8.9    Apply conservation correction
!                   -----------------------------

          ! Compute vertical integrals

       JKM = MAXVAL( ICTL(:) )
       ZWORK2(:) = 0.0_dp
       ZWORK2B(:)= 0.0_dp
       DO JK = IKB+1, JKM
         JKP = JK + 1
         DO JI = 1, ICONV
         ! ZW1 = _HALF_ *  (ZPRES(JI,JK-1) - ZPRES(JI,JKP)) / RG
           ZW1 = (ZPAH(JI,JK-1) - ZPAH(JI,JK)) / RG
           ZWORK2(JI) = ZWORK2(JI) + ZUC(JI,JK) * ZW1
           ZWORK2B(JI)= ZWORK2B(JI)+ ZVC(JI,JK) * ZW1
         ENDDO
       ENDDO

          !  error (integral must be zero)
       DO JI = 1, ICONV
         IF ( ZTPR(JI) > 0.0_dp) THEN
           JKP = ICTL(JI)  
         ! ZW1 = RG / ( ZPRES(JI,IKB) - ZPRES(JI,JKP) - _HALF_*(ZDPRES(JI,IKB+1) - ZDPRES(JI,JKP+1)) )
           ZW1 = RG / ( ZPAH(JI,IKB) - ZPAH(JI,JKP) )
           ZWORK2(JI) = ZWORK2(JI) * ZW1
           ZWORK2B(JI)= ZWORK2B(JI)* ZW1
         ENDIF
       ENDDO

          ! Apply uniform correction 

       DO JK = JKM, IKB+1, -1
         DO JI = 1, ICONV
           IF ( ZTPR(JI) > 0.0_dp .AND. JK <= ICTL(JI) ) THEN
       !       if (ABS (zwork2b(ji)) > abs(zvc(ji,jk)))  print*,zvc(ji,jk), zwork2b(ji)
                ZUC(JI,JK) = ZUC(JI,JK) - ZWORK2(JI)
                ZVC(JI,JK) = ZVC(JI,JK) - ZWORK2B(JI)
       !         if ( abs(zuc(ji,jk))*2400. > 10._dp ) print*, zuc(ji,jk), zwork2(ji), zu(ji,jk)
       ! if ( abs(zvc(ji,jk))*2400. > 10._dp ) print*, zvc(ji,jk), zwork2b(ji), zv(ji,jk)
           ENDIF
         ENDDO
       ENDDO

          ! extend tendencies to first model level
 
   !  DO JI = 1, ICONV
   !     ZWORK2(JI) = ZDPRES(JI,IKB+1) + ZDPRES(JI,IKB+2)
   !     ZUC(JI,IKB)  = ZUC(JI,IKB+1) * ZDPRES(JI,IKB+2)/ZWORK2(JI)
   !     ZUC(JI,IKB+1)= ZUC(JI,IKB+1) * ZDPRES(JI,IKB+1)/ZWORK2(JI)
   !     ZVC(JI,IKB)  = ZVC(JI,IKB+1) * ZDPRES(JI,IKB+2)/ZWORK2(JI)
   !     ZVC(JI,IKB+1)= ZVC(JI,IKB+1) * ZDPRES(JI,IKB+1)/ZWORK2(JI)
   !  ENDDO

!
       DO JK = IKB, IKE
       DO JI = 1, ICONV
          JL = IJINDEX(JI)
          PUTEN(JL,JK)   = ZUC(JI,JK)
          PVTEN(JL,JK)   = ZVC(JI,JK)
       ENDDO
       ENDDO

       DEALLOCATE( ZUC )
       DEALLOCATE( ZVC )

  ENDIF

!
!*           9.     Write up- and downdraft mass fluxes
!                   ------------------------------------
!
    DO JK = IKB, IKE
       ZUMF(:,JK)  = ZUMF(:,JK) / ZDXDY(:) ! Mass flux per unit area
       ZDMF(:,JK)  = ZDMF(:,JK) / ZDXDY(:)
    ENDDO
    DO JK = IKB, IKE
!!$       ZUDR(:,JK)  = ZUDR(:,JK) / ( ZDXDY(:) * ZDPRES(:,JK) )! detrainment for ERA40
!!$       ZDDR(:,JK)  = ZDDR(:,JK) / ( ZDXDY(:) * ZDPRES(:,JK) )
      ZUDR(:,JK)  = ZUDR(:,JK) / ZDXDY(:)! detrainment for ERA40 / MESSy
      ZDDR(:,JK)  = ZDDR(:,JK) / ZDXDY(:)
      ZUER(:,JK)  = ZUER(:,JK) / ZDXDY(:)! entrainment for ERA40 / MESSy
      ZDER(:,JK)  = ZDER(:,JK) / ZDXDY(:)
    ENDDO
    ZWORK2(:) = 1.0_dp
    DO JK = IKB, IKE
    DO JI = 1, ICONV
       JL = IJINDEX(JI)
       IF ( PPRLTEN(JL) < 1.E-14_dp ) ZWORK2(JL) = 0.0_dp
       PUMF(JL,JK)     = ZUMF(JI,JK)     * ZWORK2(JL)
       PDMF(JL,JK)     = ZDMF(JI,JK)     * ZWORK2(JL)
       PUDR(JL,JK)     = ZUDR(JI,JK)     * ZWORK2(JL)
       PDDR(JL,JK)     = ZDDR(JI,JK)     * ZWORK2(JL)
       PUER(JL,JK)     = ZUER(JI,JK)     * ZWORK2(JL)
       PDER(JL,JK)     = ZDER(JI,JK)     * ZWORK2(JL)
! remapping already done above !!!!
       PURV(JL,JK)     = PURV(JL,JK)     * ZWORK2(JL)
       PURCI(JL,JK)    = PURCI(JL,JK)    * ZWORK2(JL)
       PURCICE(JL,JK)  = PURCICE(JL,JK)  * ZWORK2(JL)
       PURCLIQ(JL,JK)  = PURCLIQ(JL,JK)  * ZWORK2(JL)
       PURCRPRO(JL,JK) = PURCRPRO(JL,JK) * ZWORK2(JL)
       PURCSPRO(JL,JK) = PURCSPRO(JL,JK) * ZWORK2(JL)
    ENDDO
    ENDDO
!
!
!*           10.    Deallocate all local arrays
!                   ---------------------------
!
! downdraft variables
!
      DEALLOCATE( ZDMF )
      DEALLOCATE( ZDER )
      DEALLOCATE( ZDDR )
      DEALLOCATE( ZDTHL )
      DEALLOCATE( ZDRW )
      DEALLOCATE( ZLMASS )
      DEALLOCATE( ZMIXF )
      DEALLOCATE( ZTPR )
      DEALLOCATE( ZSPR )
      DEALLOCATE( ZDTEVR )
      DEALLOCATE( ZPREF )
      DEALLOCATE( IML )
      DEALLOCATE( ILFS )
      DEALLOCATE( IDBL )
      DEALLOCATE( ZDTEVRF )
      DEALLOCATE( ZPRLFLX )
      DEALLOCATE( ZPRSFLX )
!
!   closure variables
!
      DEALLOCATE( ZTIMEA )
      DEALLOCATE( ZTIMEC )
      DEALLOCATE( ZTHC )
      DEALLOCATE( ZRVC )
      DEALLOCATE( ZRCC )
      DEALLOCATE( ZRIC )
      DEALLOCATE( ZWSUB )
!
       IF ( OCH1CONV ) THEN
           DEALLOCATE( ZCH1 )
           DEALLOCATE( ZCH1C )
           DEALLOCATE( ZWORK3 )
       ENDIF
!
    ENDIF
!
!    vertical index
!
    DEALLOCATE( IDPL )
    DEALLOCATE( IPBL )
    DEALLOCATE( ILCL )
    DEALLOCATE( ICTL )
    DEALLOCATE( IETL )
!
! grid scale variables
!
    DEALLOCATE( ZZ )
    DEALLOCATE( ZPRES )
    DEALLOCATE( ZPAH )
    DEALLOCATE( ZDPRES )
    DEALLOCATE( ZU )
    DEALLOCATE( ZV )
    DEALLOCATE( ZTT )
    DEALLOCATE( ZTH )
    DEALLOCATE( ZTHV )
    DEALLOCATE( ZTHL )
    DEALLOCATE( ZTHES )
    DEALLOCATE( ZRW )
    DEALLOCATE( ZRV )
    DEALLOCATE( ZRC )
    DEALLOCATE( ZRI )
    DEALLOCATE( ZDXDY )
!
! updraft variables
!
    DEALLOCATE( ZUMF )
    DEALLOCATE( ZUER )
    DEALLOCATE( ZUDR )
    DEALLOCATE( ZUTHL )
    DEALLOCATE( ZUTHV )
    DEALLOCATE( ZURW )
    DEALLOCATE( ZURC )
    DEALLOCATE( ZURI )
    DEALLOCATE( ZURR )
    DEALLOCATE( ZURS )
    DEALLOCATE( ZUPR )
    DEALLOCATE( ZUTPR )
    DEALLOCATE( ZTHLCL )
    DEALLOCATE( ZTLCL )
    DEALLOCATE( ZRVLCL )
    DEALLOCATE( ZWLCL )
    DEALLOCATE( ZZLCL )
    DEALLOCATE( ZTHVELCL )
    DEALLOCATE( ZMFLCL )
    DEALLOCATE( ZCAPE )
    IF ( OSETTADJ ) DEALLOCATE( ZTIMED )
!
! work arrays
!
    DEALLOCATE( IINDEX )
    DEALLOCATE( IJINDEX )
    DEALLOCATE( IJSINDEX )
    DEALLOCATE( GTRIG1 )
!
!
END SUBROUTINE CONVECT_DEEP

!===============================================================================
!===============================================================================

!############################################################################
 SUBROUTINE CONVECT_SHALLOW( KLON, KLEV, KIDIA, KFDIA, KBDIA, KTDIA,        &
                           & PDTCONV, KICE, OSETTADJ, PTADJS,               &
                           & PPABST, PPAH, PZZ, PHSFLX,                     &
                           & PTT, PRVT, PRCT, PRIT, PUT, PVT, PWT,          &
                           & KCOUNT, PTTEN, PRVTEN, PRCTEN, PRITEN,         &
                           & KCLTOP, KCLBAS, PUMF, PURV, PURCI,             &
                           & OUVCONV, PUTEN, PVTEN,                         &
                           & OCH1CONV, KCH1, PCH1, PCH1TEN,                 &
                           ! for ERA40 / MESSY
                           & PUDR, PUER, PURCICE, PURCLIQ          ) 
!############################################################################

!!**** Monitor routine to compute all convective tendencies by calls
!!     of several subroutines.
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine the convective
!!      tendencies. The routine first prepares all necessary grid-scale
!!      variables. The final convective tendencies are then computed by
!!      calls of different subroutines.
!!
!!
!!**  METHOD
!!    ------
!!      We start by selecting convective columns in the model domain through
!!      the call of routine TRIGGER_FUNCT. Then, we allocate memory for the
!!      convection updraft and downdraft variables and gather the grid scale
!!      variables in convective arrays.
!!      The updraft and downdraft computations are done level by level starting
!!      at the  bottom and top of the domain, respectively.
!!      All computations are done on MNH thermodynamic levels. The depth
!!      of the current model layer k is defined by DP(k)=P(k-1)-P(k)
!!
!!
!!
!!    EXTERNAL
!!    --------
!!    CONVECT_TRIGGER_SHAL
!!    CONVECT_SATMIXRATIO
!!    CONVECT_UPDRAFT_SHAL
!!        CONVECT_CONDENS
!!        CONVECT_MIXING_FUNCT
!!    CONVECT_CLOSURE_SHAL
!!        CONVECT_CLOSURE_THRVLCL
!!        CONVECT_CLOSURE_ADJUST_SHAL
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module YOMCST
!!          RG                   ! gravity constant
!!          RPI                  ! number Pi
!!          RATM                 ! reference pressure
!!          RD, RV               ! gaz  constants for dry air and water vapor
!!          RCPD, RCPV           ! specific heat for dry air and water vapor
!!          RALPW, RBETW, RGAMW  ! constants for water saturation pressure
!!          RTT                  ! triple point temperature
!!          RLVTT, RLSTT         ! vaporization, sublimation heat constant
!!          RCW, RCS             ! specific heat for liquid water and ice
!!
!!      Module YOE_CONVPAREXT
!!          JCVEXB, JCVEXT       ! extra levels on the vertical boundaries
!!
!!      Module YOE_CONVPAR
!!          XA25                 ! reference grid area
!!          XCRAD                ! cloud radius
!!
!!
!!    REFERENCE
!!    ---------
!!
!!      Bechtold, 1997 : Meso-NH scientific  documentation (31 pp)
!!      Fritsch and Chappell, 1980, J. Atmos. Sci., Vol. 37, 1722-1761.
!!      Kain and Fritsch, 1990, J. Atmos. Sci., Vol. 47, 2784-2801.
!!      Kain and Fritsch, 1993, Meteor. Monographs, Vol. 24, 165-170.
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96
!!   Peter Bechtold 15/11/96 replace theta_il by enthalpy
!!         "        10/12/98 changes for ARPEGE
!!
!!     Holger Tost  13/08/04 changes for MESSy
!-------------------------------------------------------------------------------



!*       0.    DECLARATIONS
!              ------------

USE MESSY_CONVECT_BECHTOLD_PARAM,        ONLY: XCTIME_SHAL, XA25_S

IMPLICIT NONE

!*       0.1   Declarations of dummy arguments :


INTEGER,                      INTENT(IN) :: KLON     ! horizontal dimension
INTEGER,                      INTENT(IN) :: KLEV     ! vertical dimension
INTEGER,                      INTENT(IN) :: KIDIA    ! value of the first point in x
INTEGER,                      INTENT(IN) :: KFDIA    ! value of the last point in x
INTEGER,                      INTENT(IN) :: KBDIA    ! vertical  computations start at
!                                                    ! KBDIA that is at least 1
INTEGER,                      INTENT(IN) :: KTDIA    ! vertical computations can be
                                                     ! limited to KLEV + 1 - KTDIA
                                                     ! default=1
REAL(dp),                     INTENT(IN) :: PDTCONV  ! Interval of time between two
                                                     ! calls of the deep convection
                                                     ! scheme
INTEGER,                      INTENT(IN) :: KICE     ! flag for ice ( 1 = yes,
                                                     !                0 = no ice )
LOGICAL,                      INTENT(IN) :: OSETTADJ ! logical to set convective
                                                     ! adjustment time by user
REAL(dp),                       INTENT(IN) :: PTADJS   ! user defined adjustment time (s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PTT      ! grid scale temperature at (K)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PRVT     ! grid scale water vapor (kg/kg)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PRCT     ! grid scale r_c  (kg/kg)"
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PRIT     ! grid scale r_i (kg/kg)"
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PUT      ! grid scale horiz. wind u (m/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PVT      ! grid scale horiz. wind v (m/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PWT      ! grid scale vertical
                                                       ! velocity (m/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PPABST   ! grid scale pressure (Pa)
REAL(dp), DIMENSION(KLON,0:KLEV),INTENT(IN):: PPAH     ! grid scale half-level pressure at (Pa)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PZZ      ! height of model layer (m)
REAL(dp), DIMENSION(KLON),      INTENT(IN) :: PHSFLX   ! surface turb. sensible heat flux (W/m2)
INTEGER,DIMENSION(KLON),        INTENT(IN) :: KCOUNT   ! counter for already active
                                                       ! deep convection points

REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT):: PTTEN  ! convective temperature
                                                       ! tendency (K/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT):: PRVTEN ! convective r_v tendency (1/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT):: PRCTEN ! convective r_c tendency (1/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT):: PRITEN ! convective r_i tendency (1/s)
INTEGER,  DIMENSION(KLON),      INTENT(INOUT):: KCLTOP ! cloud top level
INTEGER,  DIMENSION(KLON),      INTENT(INOUT):: KCLBAS ! cloud base level
                                                       ! they are given a value of
                                                       ! 0 if no convection
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT):: PUMF   ! updraft mass flux (kg/s m2)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT)  :: PURV   ! updraft water vapor (kg/kg)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT)  :: PURCI  ! updraft liquid+ice condensate (kg/kg)

LOGICAL,                      INTENT(IN)   :: OUVCONV! include wind transport (Cu friction)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT):: PUTEN  ! convective u tendency (m/s^2)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT):: PVTEN  ! convective v tendency (m/s^2)


LOGICAL,                      INTENT(IN) :: OCH1CONV ! include tracer transport
INTEGER,                      INTENT(IN) :: KCH1     ! number of species
REAL(dp), DIMENSION(KLON,KLEV,KCH1), INTENT(IN) :: PCH1! grid scale chemical species
REAL(dp), DIMENSION(KLON,KLEV,KCH1), INTENT(OUT):: PCH1TEN! species conv. tendency (1/s)

! for ERA40 / MESSy
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT)  :: PUDR   ! updraft detrainment rate (kg/s m2)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT)  :: PUER   ! updraft entrainment rate (kg/s m2)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PURCICE ! cv ice water (kg/kg)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PURCLIQ ! cv liq water (kg/kg)

!*       0.2   Declarations of local fixed memory variables :

INTEGER    :: ITEST, ICONV, ICONV1    ! number of convective columns
INTEGER    :: IIB, IIE                ! horizontal loop bounds
INTEGER    :: IKB, IKE                ! vertical loop bounds
INTEGER    :: IKS                     ! vertical dimension
INTEGER    :: JI, JL                  ! horizontal loop index
INTEGER    :: JN                      ! number of tracers
INTEGER    :: JK, JKP, JKM            ! vertical loop index
INTEGER    :: IFTSTEPS                ! only used for chemical tracers
REAL(dp)   :: ZEPS, ZEPSA, ZEPSB      ! R_d / R_v, R_v / R_d, RCPV / RCPD - ZEPSA
REAL(dp)   :: ZCPORD, ZRDOCP          ! C_p/R_d,  R_d/C_p

LOGICAL,   DIMENSION(KLON,KLEV)        :: GTRIG3 ! 3D logical mask for convection
LOGICAL,   DIMENSION(KLON)             :: GTRIG  ! 2D logical mask for trigger test
REAL(dp),  DIMENSION(KLON,KLEV)        :: ZTHT, ZSTHV, ZSTHES  ! grid scale theta, theta_v
REAL(dp),  DIMENSION(KLON)             :: ZWORK2, ZWORK2B ! work array


!*       0.2   Declarations of local allocatable  variables :

INTEGER,  DIMENSION(:), ALLOCATABLE    :: IDPL    ! index for parcel departure level
INTEGER,  DIMENSION(:), ALLOCATABLE    :: IPBL    ! index for source layer top
INTEGER,  DIMENSION(:), ALLOCATABLE    :: ILCL    ! index for lifting condensation level
INTEGER,  DIMENSION(:), ALLOCATABLE    :: IETL    ! index for zero buoyancy level
INTEGER,  DIMENSION(:), ALLOCATABLE    :: ICTL    ! index for cloud top level
INTEGER,  DIMENSION(:), ALLOCATABLE    :: ILFS    ! index for level of free sink

INTEGER,  DIMENSION(:), ALLOCATABLE    :: ISDPL   ! index for parcel departure level
INTEGER,  DIMENSION(:), ALLOCATABLE    :: ISPBL   ! index for source layer top
INTEGER,  DIMENSION(:), ALLOCATABLE    :: ISLCL   ! index for lifting condensation level
REAL(dp), DIMENSION(:), ALLOCATABLE    :: ZSTHLCL ! updraft theta at LCL
REAL(dp), DIMENSION(:), ALLOCATABLE    :: ZSTLCL  ! updraft temp. at LCL
REAL(dp), DIMENSION(:), ALLOCATABLE    :: ZSRVLCL ! updraft rv at LCL
REAL(dp), DIMENSION(:), ALLOCATABLE    :: ZSWLCL  ! updraft w at LCL
REAL(dp), DIMENSION(:), ALLOCATABLE    :: ZSZLCL  ! LCL height
REAL(dp), DIMENSION(:), ALLOCATABLE    :: ZSTHVELCL! envir. theta_v at LCL
REAL(dp), DIMENSION(:), ALLOCATABLE    :: ZSDXDY  ! grid area (m^2)

! grid scale variables
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZZ      ! height of model layer (m)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZPRES   ! grid scale pressure
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZPAH    ! grid scale half-level pressure at (Pa)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZDPRES  ! pressure difference between
                                                ! bottom and top of layer (Pa)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZW      ! grid scale vertical velocity
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZTT     ! temperature
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZTH     ! grid scale theta
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZTHV    ! grid scale theta_v
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZTHL    ! grid scale enthalpy (J/kg)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZTHES, ZTHEST ! grid scale saturated theta_e
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZRW     ! grid scale total water (kg/kg)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZRV     ! grid scale water vapor (kg/kg)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZRC     ! grid scale cloud water (kg/kg)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZRI     ! grid scale cloud ice (kg/kg)
REAL(dp), DIMENSION(:),   ALLOCATABLE  :: ZDXDY   ! grid area (m^2)

! updraft variables
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZUMF    ! updraft mass flux (kg/s)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZUER    ! updraft entrainment (kg/s)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZUDR    ! updraft detrainment (kg/s)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZUTHL   ! updraft enthalpy (J/kg)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZUTHV   ! updraft theta_v (K)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZURW    ! updraft total water (kg/kg)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZURC    ! updraft cloud water (kg/kg)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZURI    ! updraft cloud ice   (kg/kg)
REAL(dp), DIMENSION(:),   ALLOCATABLE  :: ZMFLCL  ! cloud base unit mass flux(kg/s)
REAL(dp), DIMENSION(:),   ALLOCATABLE  :: ZCAPE   ! available potent. energy
REAL(dp), DIMENSION(:),   ALLOCATABLE  :: ZTHLCL  ! updraft theta at LCL
REAL(dp), DIMENSION(:),   ALLOCATABLE  :: ZTLCL   ! updraft temp. at LCL
REAL(dp), DIMENSION(:),   ALLOCATABLE  :: ZRVLCL  ! updraft rv at LCL
REAL(dp), DIMENSION(:),   ALLOCATABLE  :: ZWLCL   ! updraft w at LCL
REAL(dp), DIMENSION(:),   ALLOCATABLE  :: ZZLCL   ! LCL height
REAL(dp), DIMENSION(:),   ALLOCATABLE  :: ZTHVELCL! envir. theta_v at LCL

! downdraft variables
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZDMF    ! downdraft mass flux (kg/s)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZDER    ! downdraft entrainment (kg/s)
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZDDR    ! downdraft detrainment (kg/s)

! closure variables
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZLMASS  ! mass of model layer (kg)
REAL(dp), DIMENSION(:),   ALLOCATABLE  :: ZTIMEC  ! advective time period

REAL(dp),  DIMENSION(:,:), ALLOCATABLE  :: ZTHC    ! conv. adj. grid scale theta
REAL(dp),  DIMENSION(:,:), ALLOCATABLE  :: ZRVC    ! conv. adj. grid scale r_w
REAL(dp),  DIMENSION(:,:), ALLOCATABLE  :: ZRCC    ! conv. adj. grid scale r_c
REAL(dp),  DIMENSION(:,:), ALLOCATABLE  :: ZRIC    ! conv. adj. grid scale r_i
REAL(dp),  DIMENSION(:,:), ALLOCATABLE  :: ZWSUB   ! envir. compensating subsidence (Pa/s)

LOGICAL,     DIMENSION(:), ALLOCATABLE  :: GTRIG1  ! logical mask for convection
LOGICAL,     DIMENSION(:), ALLOCATABLE  :: GWORK   ! logical work array
INTEGER,     DIMENSION(:), ALLOCATABLE  :: IINDEX, IJINDEX, IJSINDEX, IJPINDEX!hor.index
INTEGER,     DIMENSION(:), ALLOCATABLE  :: ICOUNT  ! index for deep convective points
REAL(dp),    DIMENSION(:), ALLOCATABLE  :: ZCPH    ! specific heat C_ph
REAL(dp),    DIMENSION(:), ALLOCATABLE  :: ZLV, ZLS! latent heat of vaporis., sublim.
REAL(dp)                                :: ZES     ! saturation vapor mixng ratio
REAL(dp)                                :: ZW1, ZW2! work variables

! for U, V transport:
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZU      ! grid scale horiz. u component 
REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: ZV      ! grid scale horiz. v component
REAL(dp),  DIMENSION(:,:), ALLOCATABLE :: ZUC     ! horizontal wind u (m/s)
REAL(dp),  DIMENSION(:,:), ALLOCATABLE :: ZVC     ! horizontal wind v (m/s)

! Chemical Tracers:
REAL(dp),  DIMENSION(:,:,:), ALLOCATABLE:: ZCH1    ! grid scale chemical specy (kg/kg)
REAL(dp),  DIMENSION(:,:,:), ALLOCATABLE:: ZCH1C   ! conv. adjust. chemical specy 1
REAL(dp),  DIMENSION(:,:),   ALLOCATABLE:: ZWORK3  ! conv. adjust. chemical specy 1
LOGICAL,   DIMENSION(:,:,:), ALLOCATABLE:: GTRIG4  ! logical mask

!-------------------------------------------------------------------------------


!*       0.3    Compute loop bounds
!               -------------------

IIB    = KIDIA
IIE    = KFDIA
JCVEXB = MAX( 0, KBDIA - 1 )
IKB    = 1 + JCVEXB
IKS    = KLEV
JCVEXT = MAX( 0, KTDIA - 1)
IKE    = IKS - JCVEXT


!*       0.5    Update convective counter ( where KCOUNT > 0
!               convection is still active ).
!               ---------------------------------------------

GTRIG(:)       = .FALSE.
GTRIG(IIB:IIE) = .TRUE.
ITEST          = COUNT( GTRIG(:) )

IF ( ITEST == 0 ) RETURN



!*       0.7    Reset convective tendencies to zero if convective
!               counter becomes negative
!               -------------------------------------------------

GTRIG3(:,:) = SPREAD( GTRIG(:), DIM=2, NCOPIES=IKS )

WHERE ( GTRIG3(:,:) )
    PTTEN(:,:)   = 0.0_dp
    PRVTEN(:,:)  = 0.0_dp
    PRCTEN(:,:)  = 0.0_dp
    PRITEN(:,:)  = 0.0_dp
    PUTEN(:,:)   = 0.0_dp
    PVTEN(:,:)   = 0.0_dp
    PUMF(:,:)    = 0.0_dp
    PURV(:,:)    = 0.0_dp
    PURCI(:,:)   = 0.0_dp
    PUDR(:,:)    = 0.0_dp
    PUER(:,:)    = 0.0_dp
    PURCLIQ(:,:) = 0.0_dp
    PURCICE(:,:) = 0.0_dp
END WHERE

WHERE ( GTRIG(:) )
 ! KCLTOP(:)  = 1 ! already initialized in CONVECTION
 ! KCLBAS(:)  = 1
END WHERE

IF ( OCH1CONV ) THEN
    ALLOCATE( GTRIG4(KLON,KLEV,KCH1) )
    GTRIG4(:,:,:) = SPREAD( GTRIG3(:,:), DIM=3, NCOPIES=KCH1 )
    WHERE( GTRIG4(:,:,:) ) PCH1TEN(:,:,:) = 0.0_dp
    DEALLOCATE( GTRIG4 )
ENDIF


!*       1.     Initialize  local variables
!               ----------------------------

ZEPS   = RD  / RV
ZEPSA  = RV  / RD
ZEPSB  = RCPV / RCPD - ZEPSA
ZCPORD = RCPD / RD
ZRDOCP = RD  / RCPD


!*       1.1    Set up grid scale theta, theta_v, theta_es
!               ------------------------------------------

ZTHT  (:,:) = 300._dp
ZSTHV (:,:) = 300._dp
ZSTHES(:,:) = 400._dp

DO JK = IKB, IKE
DO JI = IIB, IIE
   IF ( PPABST(JI,JK) > 40.E2_dp ) THEN
      ZTHT(JI,JK)  = PTT(JI,JK) * ( RATM / PPABST(JI,JK) ) ** ZRDOCP
      ZSTHV(JI,JK) = ZTHT(JI,JK) * ( 1.0_dp + ZEPSA * PRVT(JI,JK) ) /         &
                   & ( 1.0_dp + PRVT(JI,JK) + PRCT(JI,JK) + PRIT(JI,JK) )

          ! use conservative Bolton (1980) formula for theta_e
          ! it is used to compute CAPE for undilute parcel ascent
          ! For economical reasons we do not use routine CONVECT_SATMIXRATIO here

      ZES = EXP( RALPW - RBETW / PTT(JI,JK) - RGAMW * LOG( PTT(JI,JK) ) )
      ZES = MIN( 1.0_dp, ZEPS * ZES / ( PPABST(JI,JK) - ZES ) )
      ZSTHES(JI,JK) = PTT(JI,JK) * ( ZTHT(JI,JK) / PTT(JI,JK) ) **           &
                   &  ( 1.0_dp - 0.28_dp * ZES )                            &
                   &    * EXP( ( 3374.6525_dp / PTT(JI,JK) - 2.5403_dp ) &
                   &    * ZES * ( 1.0_dp + 0.81_dp * ZES )  )
   ENDIF
ENDDO
ENDDO



!*       2.     Test for convective columns and determine properties at the LCL
!               --------------------------------------------------------------

!*       2.1    Allocate arrays depending on number of model columns that need
!               to be tested for convection (i.e. where no convection is present
!               at the moment.
!               --------------------------------------------------------------

  ALLOCATE( ZPRES(ITEST,IKS) )
  ALLOCATE( ZZ(ITEST,IKS) )
  ALLOCATE( ZW(ITEST,IKS) )
  ALLOCATE( ZTH(ITEST,IKS) )
  ALLOCATE( ZTHV(ITEST,IKS) )
  ALLOCATE( ZTHEST(ITEST,IKS) )
  ALLOCATE( ZRV(ITEST,IKS) )
  ALLOCATE( ZSTHLCL(ITEST) )
  ALLOCATE( ZSTLCL(ITEST) )
  ALLOCATE( ZSRVLCL(ITEST) )
  ALLOCATE( ZSWLCL(ITEST) )
  ALLOCATE( ZSZLCL(ITEST) )
  ALLOCATE( ZSTHVELCL(ITEST) )
  ALLOCATE( ISDPL(ITEST) )
  ALLOCATE( ISPBL(ITEST) )
  ALLOCATE( ISLCL(ITEST) )
  ALLOCATE( ICOUNT(ITEST) )
  ALLOCATE( ZSDXDY(ITEST) )
  ALLOCATE( GTRIG1(ITEST) )
  ALLOCATE( IINDEX(KLON) )
  ALLOCATE( IJSINDEX(ITEST) )

  DO JI = 1, KLON
    IINDEX(JI) = JI
  ENDDO

  IJSINDEX(:) = PACK( IINDEX(:), MASK=GTRIG(:) )

  DO JK = IKB, IKE
  DO JI = 1, ITEST
    JL = IJSINDEX(JI)
    ZPRES(JI,JK)  = PPABST(JL,JK)
    ZZ(JI,JK)     = PZZ(JL,JK)
    ZTH(JI,JK)    = ZTHT(JL,JK)
    ZTHV(JI,JK)   = ZSTHV(JL,JK)
    ZTHEST(JI,JK) = ZSTHES(JL,JK)
    ZRV(JI,JK)    = MAX( 0.0_dp, PRVT(JL,JK) )
    ZW(JI,JK)     = PWT(JL,JK)
  ENDDO
  ENDDO
  DO JI = 1, ITEST
    JL = IJSINDEX(JI)
    ZSDXDY(JI)    = XA25_S
    ICOUNT(JI)    = KCOUNT(JL)
  ENDDO

!*       2.2    Compute environm. enthalpy and total water = r_v + r_i + r_c
!               and envir. saturation theta_e
!               ------------------------------------------------------------


!*       2.3    Test for convective columns and determine properties at the LCL
!               --------------------------------------------------------------

  ISLCL(:) = MAX( IKB, 2 )   ! initialize DPL PBL and LCL
  ISDPL(:) = IKB
  ISPBL(:) = IKB

  CALL CONVECT_TRIGGER_SHAL(  ITEST, KLEV,                              &
                           &  ZPRES, ZTH, ZTHV, ZTHEST,                 &
                           &  ZRV, ZW, ZZ, ZSDXDY, ICOUNT,              &
                           &  ZSTHLCL, ZSTLCL, ZSRVLCL, ZSWLCL, ZSZLCL, &
                           &  ZSTHVELCL, ISLCL, ISDPL, ISPBL, GTRIG1    )

  DEALLOCATE( ZPRES )
  DEALLOCATE( ZZ )
  DEALLOCATE( ZTH )
  DEALLOCATE( ZTHV )
  DEALLOCATE( ZTHEST )
  DEALLOCATE( ZRV )
  DEALLOCATE( ZW )
  DEALLOCATE( ICOUNT )


!*       3.     After the call of TRIGGER_FUNCT we allocate all the dynamic
!               arrays used in the convection scheme using the mask GTRIG, i.e.
!               we do calculus only in convective columns. This corresponds to
!               a GATHER operation.
!               --------------------------------------------------------------

  ICONV = COUNT( GTRIG1(:) )
  IF ( ICONV == 0 )  THEN
      DEALLOCATE( ZSTHLCL )
      DEALLOCATE( ZSTLCL )
      DEALLOCATE( ZSRVLCL )
      DEALLOCATE( ZSWLCL )
      DEALLOCATE( ZSZLCL )
      DEALLOCATE( ZSTHVELCL )
      DEALLOCATE( ZSDXDY )
      DEALLOCATE( ISLCL )
      DEALLOCATE( ISDPL )
      DEALLOCATE( ISPBL )
      DEALLOCATE( GTRIG1 )
      DEALLOCATE( IINDEX )
      DEALLOCATE( IJSINDEX )
      RETURN   ! no convective column has been found, exit CONVECT_SHALLOW 
  ENDIF

   ! vertical index variables

   ALLOCATE( IDPL(ICONV) )
   ALLOCATE( IPBL(ICONV) )
   ALLOCATE( ILCL(ICONV) )
   ALLOCATE( ICTL(ICONV) )
   ALLOCATE( IETL(ICONV) )

   ! grid scale variables

   ALLOCATE( ZZ(ICONV,IKS) )
   ALLOCATE( ZPRES(ICONV,IKS) )
   ALLOCATE( ZPAH(ICONV,0:IKS) )
   ALLOCATE( ZDPRES(ICONV,IKS) )
   ALLOCATE( ZTT(ICONV, IKS) )
   ALLOCATE( ZTH(ICONV,IKS) )
   ALLOCATE( ZTHV(ICONV,IKS) )
   ALLOCATE( ZTHL(ICONV,IKS) )
   ALLOCATE( ZTHES(ICONV,IKS) )
   ALLOCATE( ZRV(ICONV,IKS) )
   ALLOCATE( ZRC(ICONV,IKS) )
   ALLOCATE( ZRI(ICONV,IKS) )
   ALLOCATE( ZRW(ICONV,IKS) )
   ALLOCATE( ZDXDY(ICONV) )

   ! updraft variables

   ALLOCATE( ZUMF(ICONV,IKS) )
   ALLOCATE( ZUER(ICONV,IKS) )
   ALLOCATE( ZUDR(ICONV,IKS) )
   ALLOCATE( ZUTHL(ICONV,IKS) )
   ALLOCATE( ZUTHV(ICONV,IKS) )
   ALLOCATE( ZURW(ICONV,IKS) )
   ALLOCATE( ZURC(ICONV,IKS) )
   ALLOCATE( ZURI(ICONV,IKS) )
   ALLOCATE( ZTHLCL(ICONV) )
   ALLOCATE( ZTLCL(ICONV) )
   ALLOCATE( ZRVLCL(ICONV) )
   ALLOCATE( ZWLCL(ICONV) )
   ALLOCATE( ZMFLCL(ICONV) )
   ALLOCATE( ZZLCL(ICONV) )
   ALLOCATE( ZTHVELCL(ICONV) )
   ALLOCATE( ZCAPE(ICONV) )

   ! work variables

   ALLOCATE( IJINDEX(ICONV) )
   ALLOCATE( IJPINDEX(ICONV) )
   ALLOCATE( ZCPH(ICONV) )
   ALLOCATE( ZLV(ICONV) )
   ALLOCATE( ZLS(ICONV) )


!*           3.1    Gather grid scale and updraft base variables in
!                   arrays using mask GTRIG
!                   ---------------------------------------------------

    GTRIG(:)      = UNPACK( GTRIG1(:), MASK=GTRIG(:), FIELD=.FALSE. )
    IJINDEX(:)    = PACK( IINDEX(:), MASK=GTRIG(:) )

    DO JK = IKB, IKE
    DO JI = 1, ICONV
         JL = IJINDEX(JI)
         ZZ(JI,JK)     = PZZ(JL,JK)
         ZPRES(JI,JK)  = PPABST(JL,JK)
         ZPAH(JI,JK)   = PPAH(JL,JK)
         ZTT(JI,JK)    = PTT(JL,JK)
         ZTH(JI,JK)    = ZTHT(JL,JK)
         ZTHES(JI,JK)  = ZSTHES(JL,JK)
         ZRV(JI,JK)    = MAX( 0.0_dp, PRVT(JL,JK) )
         ZRC(JI,JK)    = MAX( 0.0_dp, PRCT(JL,JK) )
         ZRI(JI,JK)    = MAX( 0.0_dp, PRIT(JL,JK) )
         ZTHV(JI,JK)   = ZSTHV(JL,JK)
    ENDDO
    ENDDO

    ! mz_ht_20050913+   
    DO  JI = 1, ICONV
         JL = IJINDEX(JI)
         ZPAH(JI,0)    = PPAH(JL,0)
    ENDDO
    ! mz_ht_20050913-

    DO JI = 1, ITEST
       IJSINDEX(JI) = JI
    ENDDO
    IJPINDEX(:) = PACK( IJSINDEX(:), MASK=GTRIG1(:) )
    DO JI = 1, ICONV
         JL = IJPINDEX(JI)
         IDPL(JI)      = ISDPL(JL)
         IPBL(JI)      = ISPBL(JL)
         ILCL(JI)      = ISLCL(JL)
         ZTHLCL(JI)    = ZSTHLCL(JL)
         ZTLCL(JI)     = ZSTLCL(JL)
         ZRVLCL(JI)    = ZSRVLCL(JL)
         ZWLCL(JI)     = ZSWLCL(JL)
         ZZLCL(JI)     = ZSZLCL(JL)
         ZTHVELCL(JI)  = ZSTHVELCL(JL)
         ZDXDY(JI)     = ZSDXDY(JL)
    ENDDO

    ALLOCATE( GWORK(ICONV) )
    GWORK(:)  = PACK( GTRIG1(:),  MASK=GTRIG1(:) )
    DEALLOCATE( GTRIG1 )
    ALLOCATE( GTRIG1(ICONV) )
    GTRIG1(:) = GWORK(:)

    DEALLOCATE( GWORK )
    DEALLOCATE( IJPINDEX )
    DEALLOCATE( ISDPL )
    DEALLOCATE( ISPBL )
    DEALLOCATE( ISLCL )
    DEALLOCATE( ZSTHLCL )
    DEALLOCATE( ZSTLCL )
    DEALLOCATE( ZSRVLCL )
    DEALLOCATE( ZSWLCL )
    DEALLOCATE( ZSZLCL )
    DEALLOCATE( ZSTHVELCL )
    DEALLOCATE( ZSDXDY )


!*           3.2    Compute pressure difference
!                   ---------------------------------------------------
    ! mz_ht_20050913+
!        ZDPRES(:,IKB) = 0.0_dp
!        DO JK = IKB + 1, IKE
!            ZDPRES(:,JK)  = ZPRES(:,JK-1) - ZPRES(:,JK)
!        ENDDO
         ZDPRES(:,:) = 0.0_dp
         DO JK = IKB, IKE
           ZDPRES(:,jk) = ZPAH(:,jk-1) - ZPAH(:,jk)
         END DO
    ! mz_ht_20050913-

!*           3.3   Compute environm. enthalpy and total water = r_v + r_i + r_c
!                  ----------------------------------------------------------

        DO JK = IKB, IKE, 1
            ZRW(:,JK)  = ZRV(:,JK) + ZRC(:,JK) + ZRI(:,JK)
            ZCPH(:)    = RCPD + RCPV * ZRW(:,JK)
            ZLV(:)     = RLVTT + ( RCPV - RCW ) * ( ZTT(:,JK) - RTT ) ! compute L_v
            ZLS(:)     = RLSTT + ( RCPV - RCS ) * ( ZTT(:,JK) - RTT ) ! compute L_i
            ZTHL(:,JK) = ZCPH(:) * ZTT(:,JK) + ( 1.0_dp + ZRW(:,JK) ) * RG * ZZ(:,JK) &
                       & - ZLV(:) * ZRC(:,JK) - ZLS(:) * ZRI(:,JK)
        ENDDO

        DEALLOCATE( ZCPH )
        DEALLOCATE( ZLV )
        DEALLOCATE( ZLS )


!*           4.     Compute updraft properties
!                   ----------------------------

!*           4.1    Set mass flux at LCL ( here a unit mass flux with w = 1 m/s )
!                   -------------------------------------------------------------

         DO JI = 1, ICONV
            ZDXDY(JI)  = XA25_S
            ZMFLCL(JI) = XA25_S * 1.e-3_dp
            JL = IJINDEX(JI)
             ! AGrant formulation
            JK = ILCL(JI)
            ZW1=MAX(0.0_dp,MIN(1.E3_dp,ZZ(JI,JK)) )
            ZMFLCL(JI) = .04_dp*(RG/PTT(JL,IKB)*&
                      &ZW1*MAX(0.0_dp,-PHSFLX(JL))/RCPD)**.3333_dp * XA25_S
         END DO


     CALL CONVECT_UPDRAFT_SHAL( ICONV, KLEV,                                     &
                              & KICE, ZPRES, ZDPRES, ZZ, ZTHL, ZTHV, ZTHES, ZRW, &
                              & ZTHLCL, ZTLCL, ZRVLCL, ZWLCL, ZZLCL, ZTHVELCL,   &
                              & ZMFLCL, GTRIG1, ILCL, IDPL, IPBL,                &
                              & ZUMF, ZUER, ZUDR, ZUTHL, ZUTHV, ZURW,            &
                              & ZURC, ZURI, ZCAPE, ICTL, IETL                    )



!*           4.2    In routine UPDRAFT GTRIG1 has been set to false when cloud
!                   thickness is smaller than 3 km
!                   -----------------------------------------------------------
    

     ICONV1 = COUNT(GTRIG1)

     IF ( ICONV1 > 0 )  THEN

!*       4.3    Allocate memory for downdraft variables
!               ---------------------------------------

! downdraft variables

        ALLOCATE( ZDMF(ICONV,IKS) )
        ALLOCATE( ZDER(ICONV,IKS) )
        ALLOCATE( ZDDR(ICONV,IKS) )
        ALLOCATE( ILFS(ICONV) )
        ALLOCATE( ZLMASS(ICONV,IKS) )
        ZDMF(:,:) = 0.0_dp
        ZDER(:,:) = 0.0_dp
        ZDDR(:,:) = 0.0_dp
        ILFS(:)   = IKB
        DO JK = IKB, IKE
           ZLMASS(:,JK)  = ZDXDY(:) * ZDPRES(:,JK) / RG  ! mass of model layer
        ENDDO
        ZLMASS(:,IKB) = ZLMASS(:,IKB+1)

! closure variables

        ALLOCATE( ZTIMEC(ICONV) )
        ALLOCATE( ZTHC(ICONV,IKS) )
        ALLOCATE( ZRVC(ICONV,IKS) )
        ALLOCATE( ZRCC(ICONV,IKS) )
        ALLOCATE( ZRIC(ICONV,IKS) )
        ALLOCATE( ZWSUB(ICONV,IKS) )


!*           5.     Compute downdraft properties
!                   ----------------------------

        ZTIMEC(:) = XCTIME_SHAL
        IF ( OSETTADJ ) ZTIMEC(:) = PTADJS

!*           7.     Determine adjusted environmental values assuming
!                   that all available buoyant energy must be removed
!                   within an advective time step ZTIMEC.
!                   ---------------------------------------------------

       CALL CONVECT_CLOSURE_SHAL( ICONV, KLEV,                         &
                                & ZPRES, ZDPRES, ZZ, ZDXDY, ZLMASS,    &
                                & ZTHL, ZTH, ZRW, ZRC, ZRI, GTRIG1,    &
                                & ZTHC, ZRVC, ZRCC, ZRIC, ZWSUB,       &
                                & ILCL, IDPL, IPBL, ICTL,              &
                                & ZUMF, ZUER, ZUDR, ZUTHL, ZURW,       &
                                & ZURC, ZURI, ZCAPE, ZTIMEC, IFTSTEPS  )




!*           8.     Determine the final grid-scale (environmental) convective
!                   tendencies and set convective counter
!                   --------------------------------------------------------


!*           8.1    Grid scale tendencies
!                   ---------------------

          ! in order to save memory, the tendencies are temporarily stored
          ! in the tables for the adjusted grid-scale values

      DO JK = IKB, IKE
         ZTHC(:,JK) = ( ZTHC(:,JK) - ZTH(:,JK) ) / ZTIMEC(:)             &
         & * ( ZPRES(:,JK) / RATM ) ** ZRDOCP ! change theta in temperature
         ZRVC(:,JK) = ( ZRVC(:,JK) - ZRW(:,JK) + ZRC(:,JK) + ZRI(:,JK) ) &
         &                                        / ZTIMEC(:)

         ZRCC(:,JK) = ( ZRCC(:,JK) - ZRC(:,JK) ) / ZTIMEC(:)
         ZRIC(:,JK) = ( ZRIC(:,JK) - ZRI(:,JK) ) / ZTIMEC(:)
      ENDDO


!*           8.2    Apply conservation correction
!                   -----------------------------

          ! adjustment at cloud top to smooth discontinuous profiles at PBL inversions
          ! (+ - - tendencies for moisture )


       DO JI = 1, ICONV
          JK = ICTL(JI)
          JKM= MAX(1,ICTL(JI)-1)
          JKP= MAX(1,ICTL(JI)-2)
          ZRVC(JI,JKM) = ZRVC(JI,JKM) + .5_dp * ZRVC(JI,JK)
          ZRCC(JI,JKM) = ZRCC(JI,JKM) + .5_dp * ZRCC(JI,JK)
          ZRIC(JI,JKM) = ZRIC(JI,JKM) + .5_dp * ZRIC(JI,JK)
          ZTHC(JI,JKM) = ZTHC(JI,JKM) + .5_dp * ZTHC(JI,JK)
          ZRVC(JI,JKP) = ZRVC(JI,JKP) + .3_dp * ZRVC(JI,JK)
          ZRCC(JI,JKP) = ZRCC(JI,JKP) + .3_dp * ZRCC(JI,JK)
          ZRIC(JI,JKP) = ZRIC(JI,JKP) + .3_dp * ZRIC(JI,JK)
          ZTHC(JI,JKP) = ZTHC(JI,JKP) + .3_dp * ZTHC(JI,JK)
          ZRVC(JI,JK)  = .2_dp * ZRVC(JI,JK)  
          ZRCC(JI,JK)  = .2_dp * ZRCC(JI,JK) 
          ZRIC(JI,JK)  = .2_dp * ZRIC(JI,JK)
          ZTHC(JI,JK)  = .2_dp * ZTHC(JI,JK)
       END DO
         

          ! compute vertical integrals - fluxes
       
       JKM = MAXVAL( ICTL(:) )
       ZWORK2(:) = 0.0_dp
       ZWORK2B(:)= 0.0_dp
       DO JK = IKB+1, JKM
         JKP = JK + 1
         DO JI = 1, ICONV
         ! ZW1 = _HALF_ *  (ZPRES(JI,JK-1) - ZPRES(JI,JKP)) / RG
           ZW1 = (ZPAH(JI,JK-1) - ZPAH(JI,JK)) / RG
           ZWORK2(JI) = ZWORK2(JI) + ( ZRVC(JI,JK) + ZRCC(JI,JK) + ZRIC(JI,JK) ) *   & ! moisture
              &                   ZW1
         ! ZWORK2B(JI) = ZWORK2B(JI) + (                                             & ! enthalpy
         !    &                        ( RCPD + RCPV * ZRW(JI,JK) )* ZTHC(JI,JK)   - &
         !    &  ( RLVTT + ( RCPV - RCW ) * ( ZTT(JI,JK) - RTT ) ) * ZRCC(JI,JK)   - &
         !    &  ( RLSTT + ( RCPV - RCW ) * ( ZTT(JI,JK) - RTT ) ) * ZRIC(JI,JK) ) * &
         !    &                   ZW1
           ZWORK2B(JI) = ZWORK2B(JI) + (                                             & ! enthalpy
              &                          RCPD * ZTHC(JI,JK)                        - &
              &  ( RLVTT + ( RCPV - RCW ) * ( ZTT(JI,JK) - RTT ) ) * ZRCC(JI,JK)   - &
              &  ( RLSTT + ( RCPV - RCW ) * ( ZTT(JI,JK) - RTT ) ) * ZRIC(JI,JK) ) * &
              &                   ZW1
         ENDDO
       ENDDO

          ! Budget error (integral must be zero)

       DO JI = 1, ICONV
           JKP = ICTL(JI) 
           IF ( ICTL(JI) > IKB+1 ) THEN
            ! ZW1 = RG / ( ZPRES(JI,IKB) - ZPRES(JI,JKP) - &
            !       & _HALF_*(ZDPRES(JI,IKB+1) - ZDPRES(JI,JKP+1)) )
              ZW1 = RG / ( ZPAH(JI,IKB) - ZPAH(JI,JKP) )
              ZWORK2(JI) =  ZWORK2(JI) * ZW1
              ZWORK2B(JI)= ZWORK2B(JI) * ZW1
           END IF
       ENDDO

          ! Apply uniform correction

       DO JK = JKM, IKB+1, -1
         DO JI = 1, ICONV
           IF ( ICTL(JI) > IKB+1 .AND. JK <= ICTL(JI) ) THEN
              ! ZW1 = ABS(ZRVC(JI,JK)) +  ABS(ZRCC(JI,JK)) +  ABS(ZRIC(JI,JK)) + 1.E-12_dp
              ! ZRVC(JI,JK) = ZRVC(JI,JK) - ABS(ZRVC(JI,JK))/ZW1*ZWORK2(JI)           ! moisture
                ZRVC(JI,JK) = ZRVC(JI,JK) - ZWORK2(JI)                                ! moisture
              ! ZRCC(JI,JK) = ZRCC(JI,JK) - ABS(ZRCC(JI,JK))/ZW1*ZWORK2(JI)
              ! ZRIC(JI,JK) = ZRIC(JI,JK) - ABS(ZRIC(JI,JK))/ZW1*ZWORK2(JI)
              ! ZTHC(JI,JK) = ZTHC(JI,JK) + ZWORK2B(JI) / ( RCPD + RCPV * ZRW(JI,JK) )! energy
                ZTHC(JI,JK) = ZTHC(JI,JK) + ZWORK2B(JI) / RCPD                        ! energy
           END IF
         ENDDO
       ENDDO

          ! extend tendencies to first model level

    ! DO JI = 1, ICONV
    !    ZWORK2(JI)    = ZDPRES(JI,IKB+1) + ZDPRES(JI,IKB+2)
    !    ZTHC(JI,IKB)  = ZTHC(JI,IKB+1) * ZDPRES(JI,IKB+2)/ZWORK2(JI)
    !    ZTHC(JI,IKB+1)= ZTHC(JI,IKB+1) * ZDPRES(JI,IKB+1)/ZWORK2(JI)
    !    ZRVC(JI,IKB)  = ZRVC(JI,IKB+1) * ZDPRES(JI,IKB+2)/ZWORK2(JI)
    !    ZRVC(JI,IKB+1)= ZRVC(JI,IKB+1) * ZDPRES(JI,IKB+1)/ZWORK2(JI)
    ! END DO


              ! execute a "scatter"= pack command to store the tendencies in
              ! the final 2D tables

      DO JK = IKB, IKE
      DO JI = 1, ICONV
         JL = IJINDEX(JI)
         PTTEN(JL,JK)   = ZTHC(JI,JK)
         PRVTEN(JL,JK)  = ZRVC(JI,JK)
         PRCTEN(JL,JK)  = ZRCC(JI,JK)
         PRITEN(JL,JK)  = ZRIC(JI,JK)

         PURV(JL,JK)    = ZURW(JI,JK) - ZURC(JI,JK) - ZURI(JI,JK)
         PURCI(JL,JK)   = ZURC(JI,JK) + ZURI(JI,JK)
         PURCLIQ(JL,JK) = ZURC(JI,JK)
         PURCICE(JL,JK) = ZURI(JI,JK)
    
      ENDDO
      ENDDO


!                   Cloud base and top levels
!                   -------------------------

     ILCL(:) = MIN( ILCL(:), ICTL(:) )
     DO JI = 1, ICONV
        JL = IJINDEX(JI)
        KCLTOP(JL) = ICTL(JI)
        KCLBAS(JL) = ILCL(JI)
     ENDDO


!*           8.6    Compute convective tendencies for Tracers
!                   ------------------------------------------

  IF ( OCH1CONV ) THEN

       ALLOCATE( ZCH1(ICONV,IKS,KCH1) )
       ALLOCATE( ZCH1C(ICONV,IKS,KCH1) )
       ALLOCATE( ZWORK3(ICONV,KCH1) )

       DO JK = IKB, IKE
       DO JI = 1, ICONV
          JL = IJINDEX(JI)
          ZCH1(JI,JK,:) = PCH1(JL,JK,:)
       ENDDO
       ENDDO

      CALL CONVECT_CHEM_TRANSPORT( ICONV, KLEV, KCH1, ZCH1, ZCH1C,          &
                                &  IDPL, IPBL, ILCL, ICTL, ILFS, ILFS,      &
                                &  ZUMF, ZUER, ZUDR, ZDMF, ZDER, ZDDR,      &
                                &  ZTIMEC, ZDXDY, ZDMF(:,1), ZLMASS, ZWSUB, &
                                &  IFTSTEPS )

       DO JK = IKB, IKE
       DO JN = 1, KCH1
          ZCH1C(:,JK,JN) = ( ZCH1C(:,JK,JN)- ZCH1(:,JK,JN) ) / ZTIMEC(:)
       ENDDO
       ENDDO

          ! extend tendencies to first model level

     ! DO JI = 1, ICONV
     !    ZWORK2(JI) = ZDPRES(JI,IKB+1) + ZDPRES(JI,IKB+2)
     ! ENDDO
     ! DO JN = 1, KCH1
     ! DO JI = 1, ICONV
     !   ZCH1(JI,IKB,JN)  = ZCH1(JI,IKB+1,JN) * ZDPRES(JI,IKB+2)/ZWORK2(JI)
     !   ZCH1(JI,IKB+1,JN)= ZCH1(JI,IKB+1,JN) * ZDPRES(JI,IKB+1)/ZWORK2(JI)
     ! ENDDO
     ! ENDDO


!*           8.7    Apply conservation correction
!                   -----------------------------

          ! Compute vertical integrals

       JKM = MAXVAL( ICTL(:) )
       ZWORK3(:,:) = 0.0_dp
       DO JK = IKB, JKM+1
         JKP = JK + 1
         DO JI = 1, ICONV
           ZWORK3(JI,:) = ZWORK3(JI,:) + ZCH1C(JI,JK,:) *                    &
                        &   (ZPAH(JI,JK-1) - ZPAH(JI,JK)) / RG
         ENDDO
       ENDDO

          ! Mass error (integral must be zero)

       DO JI = 1, ICONV
           JKP = ICTL(JI) + 1
           IF ( ICTL(JI) > IKB+1 ) THEN
            ! ZW1 = RG / ( ZPRES(JI,IKB) - ZPRES(JI,JKP) - &
            !       & _HALF_*(ZDPRES(JI,IKB+1) - ZDPRES(JI,JKP+1)) )
              ZW1 = RG / ( ZPAH(JI,IKB) - ZPAH(JI,JKP) )
              ZWORK3(JI,:) = ZWORK3(JI,:) * ZW1
           END IF
       ENDDO

          ! Apply uniform correction but assure positive mass at each level

       DO JK = JKM, IKB, -1
         DO JI = 1, ICONV
           IF ( ICTL(JI) > IKB+1 .AND. JK <= ICTL(JI) ) THEN
                ZCH1C(JI,JK,:) = ZCH1C(JI,JK,:) - ZWORK3(JI,:)
              ! ZCH1C(JI,JK,:) = MAX( ZCH1C(JI,JK,:), -ZCH1(JI,JK,:)/ZTIMEC(JI) )
           ENDIF
         ENDDO
       ENDDO

       DO JK = IKB, IKE
       DO JI = 1, ICONV
          JL = IJINDEX(JI)
          PCH1TEN(JL,JK,:) = ZCH1C(JI,JK,:)
       ENDDO
       ENDDO
  ENDIF

!*           8.8    Compute convective tendencies for wind
!                   --------------------------------------

  IF ( OUVCONV ) THEN

       ALLOCATE( ZU(ICONV,IKS) )
       ALLOCATE( ZV(ICONV,IKS) )
       DO JK = IKB, IKE
       DO JI = 1, ICONV
         JL = IJINDEX(JI)
         ZU(JI,JK)     = PUT(JL,JK)
         ZV(JI,JK)     = PVT(JL,JK)
       ENDDO
       ENDDO

       ALLOCATE( ZUC(ICONV,IKS) )
       ALLOCATE( ZVC(ICONV,IKS) )

       CALL CONVECT_UV_TRANSPORT( ICONV, KLEV, ZU, ZV, ZUC, ZVC,           &
                                & IDPL, IPBL, ILCL, ICTL, ILFS, ILFS,      &
                                & ZUMF, ZUER, ZUDR, ZDMF, ZDER, ZDDR,      &
                                & ZTIMEC, ZDXDY, ZDMF(:,1), ZLMASS, ZWSUB, &
                                & IFTSTEPS )

       DO JK = IKB, IKE
          ZUC(:,JK) = ( ZUC(:,JK)- ZU(:,JK) ) / ZTIMEC(:)
          ZVC(:,JK) = ( ZVC(:,JK)- ZV(:,JK) ) / ZTIMEC(:)
       ENDDO

!*           8.9    Apply conservation correction
!                   -----------------------------

          ! Compute vertical integrals

       JKM = MAXVAL( ICTL(:) )
       ZWORK2(:) = 0.0_dp
       ZWORK2B(:)= 0.0_dp
       DO JK = IKB+1, JKM
         JKP = JK + 1
         DO JI = 1, ICONV
           ZW1 = (ZPAH(JI,JK-1) - ZPAH(JI,JK)) / RG
           ZWORK2(JI) = ZWORK2(JI) + ZUC(JI,JK) * ZW1
           ZWORK2B(JI)= ZWORK2B(JI)+ ZVC(JI,JK) * ZW1
         ENDDO
       ENDDO

          !  error (integral must be zero)

       DO JI = 1, ICONV
         IF ( ICTL(JI) > IKB ) THEN
           JKP = ICTL(JI)
           ZW1 = RG / ( ZPAH(JI,IKB) - ZPAH(JI,JKP) )
           ZWORK2(JI) = ZWORK2(JI) * ZW1
           ZWORK2B(JI)= ZWORK2B(JI)* ZW1
         ENDIF
       ENDDO

          ! Apply uniform correction

       DO JK = JKM, IKB+1, -1
         DO JI = 1, ICONV
           IF ( ICTL(JI) > IKB .AND. JK <= ICTL(JI) ) THEN
              !  ZUC(JI,JK) = ZUC(JI,JK) - ZWORK2(JI)
              !  ZVC(JI,JK) = ZVC(JI,JK) - ZWORK2B(JI)
           ENDIF
         ENDDO
       ENDDO

          ! extend tendencies to first model level

    ! DO JI = 1, ICONV
    !    ZWORK2(JI) = ZDPRES(JI,IKB+1) + ZDPRES(JI,IKB+2)
    !    ZUC(JI,IKB)  = ZUC(JI,IKB+1) * ZDPRES(JI,IKB+2)/ZWORK2(JI)
    !    ZUC(JI,IKB+1)= ZUC(JI,IKB+1) * ZDPRES(JI,IKB+1)/ZWORK2(JI)
    !    ZVC(JI,IKB)  = ZVC(JI,IKB+1) * ZDPRES(JI,IKB+2)/ZWORK2(JI)
    !    ZVC(JI,IKB+1)= ZVC(JI,IKB+1) * ZDPRES(JI,IKB+1)/ZWORK2(JI)
    ! ENDDO


       DO JK = IKB, IKE
       DO JI = 1, ICONV
          JL = IJINDEX(JI)
          PUTEN(JL,JK)   = ZUC(JI,JK)
          PVTEN(JL,JK)   = ZVC(JI,JK)
       ENDDO
       ENDDO

       DEALLOCATE( ZU )
       DEALLOCATE( ZV )
       DEALLOCATE( ZUC )
       DEALLOCATE( ZVC )

  ENDIF

!*           9.     Write up- and downdraft mass fluxes
!                   -----------------------------------

    DO JK = IKB, IKE
       ZUMF(:,JK)  = ZUMF(:,JK) / ZDXDY(:) ! Mass flux per unit area
    ENDDO
    DO JK = IKB+1, IKE
!       ZUDR(:,JK)  = ZUDR(:,JK) / ( ZDXDY(:) * ZDPRES(:,JK) )! detrainment for ERA40
      ZUDR(:,JK)  = ZUDR(:,JK) / ZDXDY(:)
      ZUER(:,JK)  = ZUER(:,JK) / ZDXDY(:)
    ENDDO
    ZWORK2(:) = 1.0_dp
    DO JK = IKB, IKE
    DO JI = 1, ICONV
       JL = IJINDEX(JI)
       IF ( KCLTOP(JL) <= IKB+1 ) ZWORK2(JL) = 0.0_dp
       PUMF(JL,JK)    = ZUMF(JI,JK)    * ZWORK2(JL)
       PUDR(JL,JK)    = ZUDR(JI,JK)    * ZWORK2(JL)
       PUER(JL,JK)    = ZUER(JI,JK)    * ZWORK2(JL)
       PURCI(JL,JK)   = PURCI(JI,JK)   * ZWORK2(JL)
       PURV(JL,JK)    = PURV(JI,JK)    * ZWORK2(JL)
       PURCLIQ(JL,JK) = PURCLIQ(JI,JK) * ZWORK2(JL)
       PURCICE(JL,JK) = PURCICE(JI,JK) * ZWORK2(JL)
    ENDDO
    ENDDO


!*           10.    Deallocate all local arrays
!                   ---------------------------

! downdraft variables

      DEALLOCATE( ZDMF )
      DEALLOCATE( ZDER )
      DEALLOCATE( ZDDR )
      DEALLOCATE( ILFS )
      DEALLOCATE( ZLMASS )
!
!   closure variables
!
      DEALLOCATE( ZTIMEC )
      DEALLOCATE( ZTHC )
      DEALLOCATE( ZRVC )
      DEALLOCATE( ZRCC )
      DEALLOCATE( ZRIC )
      DEALLOCATE( ZWSUB )
!
       IF ( OCH1CONV ) THEN
           DEALLOCATE( ZCH1 )
           DEALLOCATE( ZCH1C )
           DEALLOCATE( ZWORK3 )
       ENDIF
!
    ENDIF
!
!    vertical index
!
    DEALLOCATE( IDPL )
    DEALLOCATE( IPBL )
    DEALLOCATE( ILCL )
    DEALLOCATE( ICTL )
    DEALLOCATE( IETL )
!
! grid scale variables
!
    DEALLOCATE( ZZ )
    DEALLOCATE( ZPRES )
    DEALLOCATE( ZPAH )
    DEALLOCATE( ZDPRES )
    DEALLOCATE( ZTT )
    DEALLOCATE( ZTH )
    DEALLOCATE( ZTHV )
    DEALLOCATE( ZTHL )
    DEALLOCATE( ZTHES )
    DEALLOCATE( ZRW )
    DEALLOCATE( ZRV )
    DEALLOCATE( ZRC )
    DEALLOCATE( ZRI )
    DEALLOCATE( ZDXDY )
!
! updraft variables
!
    DEALLOCATE( ZUMF )
    DEALLOCATE( ZUER )
    DEALLOCATE( ZUDR )
    DEALLOCATE( ZUTHL )
    DEALLOCATE( ZUTHV )
    DEALLOCATE( ZURW )
    DEALLOCATE( ZURC )
    DEALLOCATE( ZURI )
    DEALLOCATE( ZTHLCL )
    DEALLOCATE( ZTLCL )
    DEALLOCATE( ZRVLCL )
    DEALLOCATE( ZWLCL )
    DEALLOCATE( ZZLCL )
    DEALLOCATE( ZTHVELCL )
    DEALLOCATE( ZMFLCL )
    DEALLOCATE( ZCAPE )
!
! work arrays
!
    DEALLOCATE( IINDEX )
    DEALLOCATE( IJINDEX )
    DEALLOCATE( IJSINDEX )
    DEALLOCATE( GTRIG1 )
!
!
END SUBROUTINE CONVECT_SHALLOW

!========================================================================================
END MODULE MESSY_CONVECT_BECHTOLD
