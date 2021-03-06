! ******************************************************************
! ------------------------------------------------------------------
MODULE MESSY_MAIN_GRID_TRAFO_NRGD_BASE
! ------------------------------------------------------------------
! Author: Patrick Joeckel, MPICH, Mainz, June 2002
!         Astrid  Kerkweg, UniMz, Mainz, 2012-2013
! ******************************************************************

  USE MESSY_MAIN_TOOLS,  ONLY: STR
#if defined(MESSY)
  USE MESSY_MAIN_CONSTANTS_MEM,  ONLY: SP, DP, I4, I8, strlen_long
#else
  USE typeSizes,     ONLY:   SP => FourByteReal  &
                           , DP => EightByteReal &
                           , I4 => FourByteInt   &
                           , I8 => EightByteInt
#endif
  USE MESSY_MAIN_GRID_NETCDF,   ONLY: t_narray, INIT_NARRAY, COPY_NARRAY    & 
                                    , PRINT_NARRAY                          &
                                    , VTYPE_INT, VTYPE_UNDEF, VTYPE_DOUBLE  &
                                    , VTYPE_REAL, VTYPE_BYTE, VTYPE_CHAR    &
                                    , RGMLE, RGMLEC, RGMLIC, RGMLVM, RGMLVL &
                                    , RGMLI, POSITION, RGMSG, ERRMSG 
  USE MESSY_MAIN_GRID_TRAFO,    ONLY: RG_INT, RG_EXT, RG_IDX, RG_IXF

  IMPLICIT NONE

  PUBLIC

  INTRINSIC :: ASSOCIATED, PRESENT, REAL, SIZE, MAXVAL, MINVAL, NULL &
             , INT, ABS, DBLE, MAX, MIN, SIGN, TRIM

  PRIVATE   :: ASSOCIATED, PRESENT, REAL, SIZE, MAXVAL, MINVAL, NULL &
             , INT, ABS, DBLE, MAX, MIN, SIGN, TRIM


  TYPE t_axis
     ! hyper-axis (for curvilinear coordinates)
     TYPE (t_narray)                      :: dat  ! interface bounds
     LOGICAL                              :: lm = .false. ! modulo axis ?
     INTEGER (I4)                         :: ndp = 0  ! number of dependencies
     ! FIRST IN LIST MUST BE INDEPENDENT !!!
     ! E.G. IF DIM 3 DEPENDS ON DIM 1,2, and 5
     ! dep = (/ 3,1,2,5 /)
     INTEGER     , DIMENSION(:), POINTER  :: dep => NULL()  ! dependencies
  END TYPE t_axis

INTERFACE OVL
   ! CALCULATES THE OVERLAP BETWEEN TWO INTERVALS
   !  (AS FRACTION OF THE TWO INTERVALS)
   !  > ALSO APPLICABLE TO 'MODULO' INTERVALS
#if !(defined(__SX__))
   MODULE PROCEDURE OVL_RR      ! REAL - REAL
   MODULE PROCEDURE OVL_RD      ! REAL - DOUBLE PRECISION
   MODULE PROCEDURE OVL_DR      ! DOUBLE PRECISION - REAL
#endif
   MODULE PROCEDURE OVL_DD      ! DOUBLE PRECISION - DOUBLE PRECISION
END INTERFACE

INTERFACE OVL_1D
   ! CALCULATES THE OVERLAP MATRICES BETWEEN TWO CONTINUOUS
   ! INTERVAL SEQUENCES (IN UNITS OF INTERVAL FRACTIONS)
   !  > ALSO APPLICABLE TO 'MODULO' SEQUENCES
#if !(defined(__SX__))
   MODULE PROCEDURE OVL_1D_RR      ! REAL - REAL
   MODULE PROCEDURE OVL_1D_RD      ! REAL - DOUBLE PRECISION
   MODULE PROCEDURE OVL_1D_DR      ! DOUBLE PRECISION - REAL
#endif
   MODULE PROCEDURE OVL_1D_DD      ! DOUBLE PRECISION - DOUBLE PRECISION
END INTERFACE

CONTAINS


! ------------------------------------------------------------------
SUBROUTINE INIT_AXIS(a)

  IMPLICIT NONE

  ! I/O
  TYPE (t_axis), INTENT(INOUT) :: a
  
  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'INIT_AXIS'
  INTEGER :: status


  CALL INIT_NARRAY(a%dat)
  a%lm = .false.
  a%ndp = 0
  IF (ASSOCIATED(a%dep)) THEN
     IF (SIZE(a%dep) > 0) THEN
        DEALLOCATE(a%dep, STAT=status)
        CALL ERRMSG(substr,status,1)
     END IF
     NULLIFY(a%dep)
  END IF

END SUBROUTINE INIT_AXIS
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE PRINT_AXIS(a, str)

  IMPLICIT NONE

  ! I/O
  TYPE (t_axis),    INTENT(IN) :: a
  CHARACTER(LEN=*), INTENT(IN) :: str
  
  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'PRINT_AXIS'
  INTEGER                     :: i

  write (0,*) substr, ' ', str, ' LM / NDP ',a%lm, a%ndp
  DO i= 1, a%ndp
     write (0,*) substr, ' ', str, ' DEP ', i, a%dep(i)
  END DO
  CALL PRINT_NARRAY(a%dat, str)

END SUBROUTINE PRINT_AXIS
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE COPY_AXIS(d, s)

  IMPLICIT NONE

  ! I/O
  TYPE (t_axis), INTENT(OUT) :: d
  TYPE (t_axis), INTENT(IN)  :: s

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'COPY_AXIS'
  INTEGER :: status

  CALL COPY_NARRAY(d%dat, s%dat)
  d%lm   = s%lm
  d%ndp  = s%ndp
  IF (ASSOCIATED(s%dep)) THEN
     ALLOCATE(d%dep(d%ndp),STAT=status)
     CALL ERRMSG(substr,status,1)
     d%dep(:) = s%dep(:)
  END IF

END SUBROUTINE COPY_AXIS
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE OVL_RR (sl,sr,dl,dr,fs,fd,modulo)

  IMPLICIT NONE

  ! I/O
  REAL (SP), INTENT(IN)              :: sl, sr  ! 'box' edges 
  REAL (SP), INTENT(IN)              :: dl, dr  ! 'box' edges
  REAL (DP), INTENT(IN),  OPTIONAL   :: modulo  ! shift for 'modulo' axis
  REAL (DP), INTENT(OUT)             :: fs, fd  ! overlap fractions

  ! LOCAL
  REAL (DP) :: x1, x2
  REAL (DP) :: shift
  INTEGER   :: fx

  x1 = MAX(MIN(dl,dr),MIN(sl,sr))
  x2 = MIN(MAX(dl,dr),MAX(sl,sr))
  fx = INT(SIGN(REAL(-0.5, DP),x2-x1)+REAL(1., DP))
  fs = (x2-x1)*REAL(fx, DP)/ABS(sr-sl)
  fd = (x2-x1)*REAL(fx, DP)/ABS(dr-dl)

  IF (PRESENT(modulo)) THEN
     
     shift = REAL(INT((dl/modulo) + REAL(1.5, DP)), DP) * modulo
     x1 = MAX(DBLE(MIN(dl,dr))+shift,DBLE(MIN(sl,sr)))
     x2 = MIN(DBLE(MAX(dl,dr))+shift,DBLE(MAX(sl,sr)))
     fx = INT(SIGN(REAL(-0.5, DP),x2-x1)+REAL(1., DP))
     fs = fs + (x2-x1)*REAL(fx, DP)/ABS(sr-sl)
     fd = fd + (x2-x1)*REAL(fx, DP)/ABS(dr-dl)

     shift = REAL(INT((sl/modulo) + REAL(1.5, DP)), DP) * modulo
     x1 = MAX(DBLE(MIN(dl,dr)),DBLE(MIN(sl,sr))+shift)
     x2 = MIN(DBLE(MAX(dl,dr)),DBLE(MAX(sl,sr))+shift)
     fx = INT(SIGN(REAl(-0.5, DP),x2-x1)+REAL(1., DP))
     fs = fs + (x2-x1)*REAL(fx, DP)/ABS(sr-sl)
     fd = fd + (x2-x1)*REAL(fx, DP)/ABS(dr-dl)
  END IF

END SUBROUTINE OVL_RR
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE OVL_RD (sl,sr,dl,dr,fs,fd,modulo)

  IMPLICIT NONE

  ! I/O
  REAL (SP), INTENT(IN)             :: sl, sr  ! 'box' edges
  REAL (DP), INTENT(IN)             :: dl, dr  ! 'box' edges
  REAL (DP), INTENT(IN),  OPTIONAL  :: modulo  ! shift for 'modulo' axis
  REAL (DP), INTENT(OUT)            :: fs, fd  ! overlap fractions

  ! LOCAL
  REAL (DP) :: x1, x2
  REAL (DP) :: shift
  INTEGER   :: fx

  x1 = MAX(MIN(dl,dr),DBLE(MIN(sl,sr)))
  x2 = MIN(MAX(dl,dr),DBLE(MAX(sl,sr)))
  fx = INT(SIGN(REAL(-0.5, DP),x2-x1)+REAl(1., DP))
  fs = (x2-x1)*REAL(fx, DP)/ABS(sr-sl)
  fd = (x2-x1)*REAL(fx, DP)/ABS(dr-dl)

  IF (PRESENT(modulo)) THEN

     shift = REAL(INT((dl/modulo) + REAL(1.5, DP)), DP) * modulo
     x1 = MAX(MIN(dl,dr)+shift,DBLE(MIN(sl,sr)))
     x2 = MIN(MAX(dl,dr)+shift,DBLE(MAX(sl,sr)))
     fx = INT(SIGN(REAL(-0.5, DP),x2-x1)+REAL(1., DP))
     fs = fs + (x2-x1)*REAL(fx, DP)/ABS(sr-sl)
     fd = fd + (x2-x1)*REAL(fx, DP)/ABS(dr-dl)

     shift = REAL(INT((sl/modulo) + REAL(1.5, DP)), DP) * modulo
     x1 = MAX(MIN(dl,dr),DBLE(MIN(sl,sr))+shift)
     x2 = MIN(MAX(dl,dr),DBLE(MAX(sl,sr))+shift)
     fx = INT(SIGN(REAL(-0.5, DP),x2-x1)+REAL(1., DP))
     fs = fs + (x2-x1)*REAL(fx, DP)/ABS(sr-sl)
     fd = fd + (x2-x1)*REAL(fx, DP)/ABS(dr-dl)
  END IF

END SUBROUTINE OVL_RD
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE OVL_DR (sl,sr,dl,dr,fs,fd,modulo)

  IMPLICIT NONE

  ! I/O
  REAL (DP), INTENT(IN)             :: sl, sr  ! 'box' edges
  REAL (SP), INTENT(IN)             :: dl, dr  ! 'box' edges
  REAL (DP), INTENT(IN),  OPTIONAL  :: modulo  ! shift for 'modulo' axis
  REAL (DP), INTENT(OUT)            :: fs, fd  ! overlap fractions

  ! LOCAL
  REAL (DP) :: x1, x2
  REAL (DP) :: shift
  INTEGER   :: fx

  x1 = MAX(DBLE(MIN(dl,dr)),MIN(sl,sr))
  x2 = MIN(DBLE(MAX(dl,dr)),MAX(sl,sr))
  fx = INT(SIGN(REAL(-0.5, DP),x2-x1)+REAL(1., DP))
  fs = (x2-x1)*REAL(fx, DP)/ABS(sr-sl)
  fd = (x2-x1)*REAL(fx, DP)/ABS(dr-dl)

  IF (PRESENT(modulo)) THEN

     shift = REAL(INT((dl/modulo) + REAL(1.5, DP)), DP) * modulo
     x1 = MAX(DBLE(MIN(dl,dr))+shift,MIN(sl,sr))
     x2 = MIN(DBLE(MAX(dl,dr))+shift,MAX(sl,sr))
     fx = INT(SIGN(REAL(-0.5, DP),x2-x1)+REAL(1., DP))
     fs = fs + (x2-x1)*REAL(fx, DP)/ABS(sr-sl)
     fd = fd + (x2-x1)*REAL(fx, DP)/ABS(dr-dl)

     shift = REAL(INT((sl/modulo) + REAL(1.5, DP)), DP) * modulo
     x1 = MAX(DBLE(MIN(dl,dr)),MIN(sl,sr)+shift)
     x2 = MIN(DBLE(MAX(dl,dr)),MAX(sl,sr)+shift)
     fx = INT(SIGN(REAL(-0.5, DP),x2-x1)+REAL(1., DP))
     fs = fs + (x2-x1)*REAL(fx, DP)/ABS(sr-sl)
     fd = fd + (x2-x1)*REAL(fx, DP)/ABS(dr-dl)
  END IF

END SUBROUTINE OVL_DR
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE OVL_DD (sl,sr,dl,dr,fs,fd,modulo)

  IMPLICIT NONE

  ! I/O
  REAL (DP), INTENT(IN)             :: sl, sr  ! 'box' edges
  REAL (DP), INTENT(IN)             :: dl, dr  ! 'box' edges
  REAL (DP), INTENT(IN),  OPTIONAL  :: modulo  ! shift for 'modulo' axis
  REAL (DP), INTENT(OUT)            :: fs, fd  ! overlap fractions

  ! LOCAL
  REAL (DP) :: x1, x2
  REAL (DP) :: shift
  INTEGER   :: fx
  !CHARACTER(LEN=*), PARAMETER :: substr = 'OVL_DD'

  x1 = MAX(MIN(dl,dr),MIN(sl,sr))
  x2 = MIN(MAX(dl,dr),MAX(sl,sr))
  fx = INT(SIGN(REAL(-0.5, DP),x2-x1)+REAL(1., DP))
  fs = (x2-x1)*REAL(fx, DP)/ABS(sr-sl)
  fd = (x2-x1)*REAL(fx, DP)/ABS(dr-dl)

  IF (PRESENT(modulo)) THEN
     shift = REAL(INT((dl/modulo) + REAL(1.5, DP)), DP) * modulo
     x1 = MAX(MIN(dl,dr)+shift,MIN(sl,sr))
     x2 = MIN(MAX(dl,dr)+shift,MAX(sl,sr))
     fx = INT(SIGN(REAL(-0.5, DP),x2-x1)+REAL(1., DP))
     fs = fs + (x2-x1)*REAL(fx, DP)/ABS(sr-sl)
     fd = fd + (x2-x1)*REAL(fx, DP)/ABS(dr-dl)

     shift = REAL(INT((sl/modulo) + REAL(1.5, DP)), DP) * modulo
     x1 = MAX(MIN(dl,dr),MIN(sl,sr)+shift)
     x2 = MIN(MAX(dl,dr),MAX(sl,sr)+shift)
     fx = INT(SIGN(REAL(-0.5, DP),x2-x1)+REAL(1., DP))
     fs = fs + (x2-x1)*REAL(fx, DP)/ABS(sr-sl)
     fd = fd + (x2-x1)*REAL(fx, DP)/ABS(dr-dl)

  END IF

END SUBROUTINE OVL_DD
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE OVL_1D_RR(s, d, mfs, mfd, lmod)

  IMPLICIT NONE

  ! I/O
  REAL (SP), DIMENSION(:), INTENT(IN) :: s        ! source bounds
  REAL (SP), DIMENSION(:), INTENT(IN) :: d        ! dest. bounds
  REAL (DP), DIMENSION(:,:)           :: mfs, mfd ! overlap matrices
  LOGICAL,                 INTENT(IN) :: lmod     ! modulo axis ?

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'OVL_1D_RR'
  INTEGER   :: n,m    ! dimensions
  INTEGER   :: i,j    ! counter
  REAL (DP) :: vmod   ! modulo value

  ! INIT
  n = SIZE(s)-1
  m = SIZE(d)-1

  ! CHECK CONSISTENCY
  IF ((SIZE(mfs,1) /= n).OR.(SIZE(mfs,2) /= m)) THEN
     CALL RGMSG(substr, RGMLE, 'INVALID SIZE OF MATRIX !')
  END IF
  IF ((SIZE(mfd,2) /= n).OR.(SIZE(mfd,1) /= m)) THEN
     CALL RGMSG(substr, RGMLE, 'INVALID SIZE OF MATRIX !')
  END IF

  IF (lmod) THEN
     vmod = ABS(s(n+1))+ABS(s(1))
     DO i=1, n
        DO j=1, m
           CALL OVL_RR(s(i),s(i+1),d(j),d(j+1),mfs(i,j),mfd(j,i),vmod)
        END DO
     END DO
  ELSE
     DO i=1, n
        DO j=1, m
           CALL OVL_RR(s(i),s(i+1),d(j),d(j+1),mfs(i,j),mfd(j,i))
        END DO
     END DO
  END IF

END SUBROUTINE OVL_1D_RR
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE OVL_1D_RD(s, d, mfs, mfd, lmod)

  IMPLICIT NONE

  ! I/O
  REAL (SP), DIMENSION(:), INTENT(IN) :: s        ! source bounds
  REAL (DP), DIMENSION(:), INTENT(IN) :: d        ! dest. bounds
  REAL (DP), DIMENSION(:,:)           :: mfs, mfd ! overlap matrices
  LOGICAL,                 INTENT(IN) :: lmod     ! modulo axis ?

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'OVL_1D_RD'
  INTEGER   :: n,m    ! dimensions
  INTEGER   :: i,j    ! counter
!  INTEGER   :: status
  REAL (DP) :: vmod   ! modulo value

  ! INIT
  n = SIZE(s)-1
  m = SIZE(d)-1

  ! CHECK CONSISTENCY
  IF ((SIZE(mfs,1) /= n).OR.(SIZE(mfs,2) /= m)) THEN
     CALL RGMSG(substr, RGMLE, 'INVALID SIZE OF MATRIX !')
  END IF
  IF ((SIZE(mfd,2) /= n).OR.(SIZE(mfd,1) /= m)) THEN
     CALL RGMSG(substr, RGMLE, 'INVALID SIZE OF MATRIX !')
  END IF

  IF (lmod) THEN
     vmod = ABS(s(n+1))+ABS(s(1))
     DO i=1, n
        DO j=1, m
           CALL OVL_RD(s(i),s(i+1),d(j),d(j+1),mfs(i,j),mfd(j,i),vmod)
        END DO
     END DO
  ELSE
     DO i=1, n
        DO j=1, m
           CALL OVL_RD(s(i),s(i+1),d(j),d(j+1),mfs(i,j),mfd(j,i))
        END DO
     END DO
  END IF

END SUBROUTINE OVL_1D_RD
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE OVL_1D_DR(s, d, mfs, mfd, lmod)

  IMPLICIT NONE

  ! I/O
  REAL (DP), DIMENSION(:), INTENT(IN) :: s        ! source bounds
  REAL (SP), DIMENSION(:), INTENT(IN) :: d        ! dest. bounds
  REAL (DP), DIMENSION(:,:)           :: mfs, mfd ! overlap matrices
  LOGICAL,                 INTENT(IN) :: lmod     ! modulo axis ?

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'OVL_1D_DR'
  INTEGER   :: n,m    ! dimensions
  INTEGER   :: i,j    ! counter
!  INTEGER   :: status
  REAL (DP) :: vmod   ! modulo value

  ! INIT
  n = SIZE(s)-1
  m = SIZE(d)-1

  ! CHECK CONSISTENCY
  IF ((SIZE(mfs,1) /= n).OR.(SIZE(mfs,2) /= m)) THEN
     CALL RGMSG(substr, RGMLE, 'INVALID SIZE OF MATRIX !')
  END IF
  IF ((SIZE(mfd,2) /= n).OR.(SIZE(mfd,1) /= m)) THEN
     CALL RGMSG(substr, RGMLE, 'INVALID SIZE OF MATRIX !')
  END IF

  IF (lmod) THEN
     vmod = ABS(s(n+1))+ABS(s(1))
     DO i=1, n
        DO j=1, m
           CALL OVL_DR(s(i),s(i+1),d(j),d(j+1),mfs(i,j),mfd(j,i),vmod)
        END DO
     END DO
  ELSE
     DO i=1, n
        DO j=1, m
           CALL OVL_DR(s(i),s(i+1),d(j),d(j+1),mfs(i,j),mfd(j,i))
        END DO
     END DO
  END IF

END SUBROUTINE OVL_1D_DR
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE OVL_1D_DD(s, d, mfs, mfd, lmod)

  IMPLICIT NONE

  ! I/O
  REAL (DP), DIMENSION(:), INTENT(IN) :: s        ! source bounds
  REAL (DP), DIMENSION(:), INTENT(IN) :: d        ! dest. bounds
  REAL (DP), DIMENSION(:,:)           :: mfs, mfd ! overlap matrices
  LOGICAL,                 INTENT(IN) :: lmod     ! modulo axis ?

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'OVL_1D_DD'
  INTEGER   :: n,m    ! dimensions
  INTEGER   :: i,j    ! counter
!  INTEGER   :: status
  REAL (DP) :: vmod   ! modulo value

  ! INIT
  n = SIZE(s)-1
  m = SIZE(d)-1

  ! CHECK CONSISTENCY
  IF ((SIZE(mfs,1) /= n).OR.(SIZE(mfs,2) /= m)) THEN
     CALL RGMSG(substr, RGMLE, 'INVALID SIZE OF MATRIX !')
  END IF
  IF ((SIZE(mfd,2) /= n).OR.(SIZE(mfd,1) /= m)) THEN
     CALL RGMSG(substr, RGMLE, 'INVALID SIZE OF MATRIX !')
  END IF

  IF (lmod) THEN
     vmod = ABS(s(n+1))+ABS(s(1))
     DO i=1, n
        DO j=1, m
          CALL OVL_DD(s(i),s(i+1),d(j),d(j+1),mfs(i,j),mfd(j,i),vmod)
        END DO
     END DO
  ELSE
     DO i=1, n
        DO j=1, m
           CALL OVL_DD(s(i),s(i+1),d(j),d(j+1),mfs(i,j),mfd(j,i))
        END DO
     END DO
  END IF

END SUBROUTINE OVL_1D_DD
! ------------------------------------------------------------------

! ------------------------------------------------------------------
RECURSIVE SUBROUTINE NREGRID(s, sax, dax, d, RG_TYPE       &
                             ,sovl,  dovl, rcnt            &
                             ,gm, gnr, gsvec, gdvec, gdio  &
                             ,gsf, gdf                     &
                            )

  IMPLICIT NONE

  ! I/O
  TYPE (t_narray), DIMENSION(:), INTENT(IN)    :: s       ! source data fields
  TYPE (t_narray), DIMENSION(:), POINTER       :: d       ! dest. data fields
  TYPE (t_axis)  , DIMENSION(:), INTENT(IN)    :: sax     ! axes (source)
  TYPE (t_axis)  , DIMENSION(:), INTENT(INOUT) :: dax     ! axes (dest)
  INTEGER        , DIMENSION(:), INTENT(IN)    :: RG_TYPE ! regridding type
  REAL (DP),       DIMENSION(:), POINTER       :: sovl  ! glob. overlap fraction
  REAL (DP),       DIMENSION(:), POINTER       :: dovl  ! glob. overlap fraction
  INTEGER,         DIMENSION(:), POINTER       :: rcnt  ! recursion level counter

  ! FOR RECURSION
  INTEGER,      OPTIONAL, INTENT(IN)            :: gm    ! current dimension
  INTEGER,      OPTIONAL, INTENT(IN)            :: gnr   ! no of dims to regrid
  INTEGER,      OPTIONAL, DIMENSION(:)          :: gsvec ! cnt. vector (source)
  INTEGER,      OPTIONAL, DIMENSION(:)          :: gdvec ! cnt vector (dest.)
  INTEGER,      OPTIONAL, DIMENSION(:)          :: gdio  ! dimension order
  REAL (DP),    OPTIONAL, INTENT(IN)   :: gsf   ! current overlap fraction
  REAL (DP),    OPTIONAL, INTENT(IN)   :: gdf   ! current overlap fraction

  ! FOR 1st STEP -> RECURSION
  INTEGER, DIMENSION(:), ALLOCATABLE :: dio  ! dimension order
  INTEGER                        :: m     ! current dimenion
  INTEGER                        :: nr    ! no of dims to regrid
  INTEGER, DIMENSION(:), ALLOCATABLE :: svec ! cnt. vector (source)
  INTEGER, DIMENSION(:), ALLOCATABLE :: dvec ! cnt. vector (dest.)
  REAL (DP)                      :: sf    ! current overlap fraction
  REAL (DP)                      :: df    ! current overlap fraction
  REAL (DP)                      :: sf0   ! current overlap fraction
  REAL (DP)                      :: df0   ! current overlap fraction

  ! FOR INVARIANT DEPENDENT AXIS
  TYPE (t_axis)  , DIMENSION(:), ALLOCATABLE   :: psax   ! axes (source)
  TYPE (t_axis)  , DIMENSION(:), ALLOCATABLE   :: pdax   ! axes (dest)
  TYPE (t_narray), DIMENSION(:), POINTER       :: pd     ! dest. data
  REAL (DP),       DIMENSION(:), POINTER       :: psovl  ! glob. overlap fraction
  REAL (DP),       DIMENSION(:), POINTER       :: pdovl  ! glob. overlap fraction
  INTEGER,         DIMENSION(:), POINTER       :: prcnt ! recursion level counter

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'NREGRID'
  INTEGER                                 :: i,j, k ! counter
  INTEGER                                 :: nvar   ! number of fields
  INTEGER                                 :: id     ! dimenions loop ounter
  INTEGER                                 :: nid    ! count no. of dimensions
  INTEGER                                 :: ndim   ! number of dimensions
  REAL (DP) , DIMENSION(:), ALLOCATABLE   :: sb     ! source bounds
  REAL (DP) , DIMENSION(:), ALLOCATABLE   :: db     ! dest. bounds
  INTEGER                                 :: vo     ! offset for bound vec.
  INTEGER                                 :: vl     ! length of bound vec.
  INTEGER,    DIMENSION(:), ALLOCATABLE   :: vec    ! element vector
  REAL (DP),  DIMENSION(:,:), ALLOCATABLE :: ms ! overlap matrix
  REAL (DP),  DIMENSION(:,:), ALLOCATABLE :: md ! overlap matrix
  INTEGER                                 :: ns, nd ! dim. of ovl-matrices
  INTEGER                                 :: status ! error status
  INTEGER                                 :: vtype, svtype  ! variable type
  LOGICAL                                 :: lflag  ! switch
  LOGICAL                                 :: lsdef  ! flag for defined dim.
  LOGICAL                                 :: lddef  ! flag for defined dim.
  LOGICAL                                 :: lsdep  ! flag for dependent dim.
  INTEGER                                 :: novl   ! number of overlaps
  REAL (DP)                               :: zso, zdo ! local overlap
  INTEGER, DIMENSION(:), ALLOCATABLE      :: hio    ! for output
  CHARACTER(LEN=1000)                     :: hstr   ! for output

  ! INIT
  ! NUMBER OF DIMENSIONS
  ndim = SIZE(sax)

  IF (PRESENT(gm)) THEN
     m = gm
  ELSE
     m = 1
  END IF

  IF (PRESENT(gnr)) THEN
     nr = gnr
  ELSE
     nr = 0
  END IF

  IF (PRESENT(gsvec)) THEN
     ALLOCATE(svec(SIZE(gsvec)), STAT=status)
     CALL ERRMSG(substr,status,1)
     svec(:) = gsvec(:)
  ELSE
     ALLOCATE(svec(s(1)%n),STAT=status)
     CALL ERRMSG(substr,status,2)
     svec(:) = 1
  END IF

  IF (PRESENT(gdvec)) THEN
     ALLOCATE(dvec(SIZE(gdvec)), STAT=status)
     CALL ERRMSG(substr,status,3)
     dvec(:) = gdvec(:)
  ELSE
     ALLOCATE(dvec(s(1)%n),STAT=status)
     CALL ERRMSG(substr,status,4)
     dvec(:) = 1
  END IF

  IF (PRESENT(gdio)) THEN
     ALLOCATE(dio(SIZE(gdio)), STAT=status)
     CALL ERRMSG(substr,status,5)
     dio(:) = gdio(:)
  ELSE
     ALLOCATE(dio(ndim),STAT=status)
     CALL ERRMSG(substr,status,6)
     dio(:) = 0
  END IF

  IF (PRESENT(gsf)) THEN
     sf0 = gsf
  ELSE
     sf0 = REAL(1.0, DP)
  END IF

  IF (PRESENT(gdf)) THEN
     df0 = gdf
  ELSE
     df0 = REAL(1.0, DP)
  END IF

  nvar = SIZE(s)

  ! ....................................................................
  ! (A) ONLY TO BE DONE ONCE (FIRST STEP)
  IF (m == 1) THEN
     ! 1. CHECK SOME BASIC REQUIREMENTS
     IF (SIZE(RG_TYPE) /= nvar) THEN
       CALL RGMSG(substr,RGMLE, &
            'NUMBER OF REGRIDDING TYPE MISMATCH IN SOURCE !')
     END IF
     !
     DO k=1, nvar
           CALL RGMSG(substr,RGMLVL,'FOR VARIABLE ',k,':', .false.)
           CALL RGMSG(substr,RGMLVL,'VAR_DIMENSIONS: ',s(k)%n,' ', .false.)
           CALL RGMSG(substr,RGMLVL,'AXES          : ',SIZE(sax),' ')
        IF (s(k)%n /= SIZE(sax)) THEN
           CALL RGMSG(substr,RGMLE,'NUMBER OF DIMENSIONS MISMATCH', .false.)
           CALL RGMSG(substr,RGMLEC,'FOR VARIABLE ',k,':', .false.)
           CALL RGMSG(substr,RGMLEC,'VAR_DIMENSIONS: ',s(k)%n,' ', .false.)
           CALL RGMSG(substr,RGMLEC,'AXES          : ',SIZE(sax),' ')
        END IF
     END DO
     !
     DO k=1, nvar
        DO i=1, s(k)%n
           vtype = sax(i)%dat%type
           IF (vtype /= VTYPE_UNDEF) THEN
              IF (s(k)%dim(i) /= sax(i)%dat%dim(1)-1) THEN
                 CALL RGMSG(substr,RGMLE,'SOURCE DIMENSION MISMATCH', &
                      .false.)
                 CALL RGMSG(substr,RGMLEC,'FOR DIMENSION ',i,' ', .false.)
                 CALL RGMSG(substr,RGMLEC,'OF VARIABLE ',k,' ', .false.)
                 CALL RGMSG(substr,RGMLEC,'BETWEEN DATA WITH LENGTH ', &
                      s(k)%dim(i),' ', .false.)
                 CALL RGMSG(substr,RGMLEC,'AND AXIS WITH LENGTH ',     &
                      sax(i)%dat%dim(1),' ')
              END IF
           END IF
        END DO
     END DO
     !
     IF (SIZE(dax) /= SIZE(sax)) THEN
        CALL RGMSG(substr,RGMLE,'NUMBER OF DIMENSIONS OF SOURCE: ',  &
             SIZE(sax),' ',.false.)
        CALL RGMSG(substr,RGMLEC,'NUMBER OF DIMENSIONS OF DEST. : ',  &
             SIZE(dax),' ')
     END IF
     !

     ! 2. CHECK DIMENSIONS
     nid = 0
     CALL RGMSG(substr, RGMLI, 'SCANNING ',ndim, ' DIMENSIONS ...')
     ! 2.1: SOURCE DIMESNION TYPE
     DO id=1,ndim
        ! 2.1.1.: SOURCE DIM MUST BE REAL, DOUBLE, OR UNDEFINED
        svtype = sax(id)%dat%type
        SELECT CASE(svtype)
           CASE(VTYPE_REAL)
              CALL RGMSG(substr, RGMLIC, '... SOURCE DIMENSION',id,  &
                   ' IS OF TYPE REAL' )
              lsdef = .true.
           CASE(VTYPE_DOUBLE)
              CALL RGMSG(substr, RGMLIC, '... SOURCE DIMENSION',id,  &
                   ' IS OF TYPE DOUBLE PRECISION' )
              lsdef = .true.
           CASE(VTYPE_UNDEF)
              CALL RGMSG(substr, RGMLIC, '... SOURCE DIMENSION',id,  &
                   ' IS UNDEFINED' )
              lsdef = .false.
           CASE DEFAULT
              !CASE(VTYPE_INT)
              !CASE(VTYPE_CHAR)
              !CASE(VTYPE_BYTE)
              CALL RGMSG(substr, RGMLE, 'REGRIDDING IS ONLY POSSIBLE FOR', &
                   .false.)
              CALL RGMSG(substr, RGMLEC,                                  &
                   'DIMENSIONS OF TYPE REAL OR DOUBLE PRECISION !')
        END SELECT
        ! 2.1.2 CHECK DEPENDENCIES
        IF (sax(id)%ndp > 1) THEN
           lsdep = .true.
           DO i=1, sax(id)%ndp
              IF (sax(sax(id)%dep(i))%dat%type == VTYPE_UNDEF) THEN
                 CALL RGMSG(substr, RGMLE, 'SOURCE DIMENSION ', id,     &
                      ' IS DEPENDENT ON', .false.)
                 CALL RGMSG(substr, RGMLEC, 'DIMENSION ',sax(id)%dep(i), &
                      ' WHICH IS UNDEFINED !')
              END IF
           END DO
        ELSE
           lsdep = .false.
        END IF
        ! 2.2.1.: DEST. DIM MUST BE REAL, DOUBLE, OR UNDEFINED
        vtype = dax(id)%dat%type
        SELECT CASE(vtype)
           CASE(VTYPE_REAL)
              CALL RGMSG(substr, RGMLIC, '... DEST.  DIMENSION',id,  &
                   ' IS OF TYPE REAL' )
              lddef = .true.
           CASE(VTYPE_DOUBLE)
              CALL RGMSG(substr, RGMLIC, '... DEST.  DIMENSION',id,  &
                   ' IS OF TYPE DOUBLE PRECISION')
              lddef = .true.
           CASE(VTYPE_UNDEF)
              CALL RGMSG(substr, RGMLIC, '... DEST.  DIMENSION',id,  &
                   ' IS UNDEFINED')
              lddef = .false.
           CASE DEFAULT
              !CASE(VTYPE_INT)
              !CASE(VTYPE_CHAR)
              !CASE(VTYPE_BYTE)
              CALL RGMSG(substr, RGMLE, 'REGRIDDING IS ONLY POSSIBLE FOR', &
                   .false.)
              CALL RGMSG(substr, RGMLEC,                                 &
                   'DIMENSIONS OF TYPE REAL OR DOUBLE PRECISION !')
        END SELECT
        ! 2.2.2 CHECK DEPENDENCIES
        IF (dax(id)%ndp > 1) THEN
           DO i=1, dax(id)%ndp
              IF (dax(dax(id)%dep(i))%dat%type == VTYPE_UNDEF) THEN
                 CALL RGMSG(substr, RGMLE, 'DEST.  DIMENSION ', id,     &
                      ' IS DEPENDENT ON', .false.)
                 CALL RGMSG(substr, RGMLEC, 'DIMENSION ',dax(id)%dep(i), &
                      ' WHICH IS UNDEFINED !')
              END IF
           END DO
        END IF
        !
        ! 2.3 CHECK CONSISTENCY BETWEEN SOURCE AND DEST.
        ! 2.3.1 UNDEF dimensions cannot be regridded
        IF ((.NOT.lsdef).AND.(lddef)) THEN
           CALL RGMSG(substr, RGMLE,  &
                'DIMENSION OF TYPE UNDEFINED', .false.)
           CALL RGMSG(substr, RGMLEC, &
                'CANNOT BE REGRIDDED ! FOR INVARIANT', .false.)
           CALL RGMSG(substr, RGMLEC, &
                'DIMENSION, DEST. DIMENSION MUST ALSO', .false.)
           CALL RGMSG(substr, RGMLEC, 'BE UNDEFINED !')
        END IF
        ! 2.3.2 INVARIANT, DEPENDENT DIMENSIONS HAVE TO BE PRE-REGRIDDED
        IF (lsdep.AND.(.NOT.lddef)) THEN
           CALL RGMSG(substr, RGMLI, 'FOUND INVARIANT DEPENDENT DIMENSION')
           CALL RGMSG(substr, RGMLIC, &
                ' >>> START REGRIDDING OF DEPENDENT DIMENSION ...')
           ! REGRIDDING AXES
           dax(id)%lm = sax(id)%lm
           dax(id)%ndp = sax(id)%ndp
           ALLOCATE(dax(id)%dep(dax(id)%ndp),STAT=status)
           CALL ERRMSG(substr,status,7)
           dax(id)%dep(:) = sax(id)%dep(:)
           ALLOCATE(psax(sax(id)%ndp),STAT=status)
           CALL ERRMSG(substr,status,8)
           ALLOCATE(pdax(dax(id)%ndp),STAT=status)
           CALL ERRMSG(substr,status,9)
           CALL INIT_AXIS(psax(1)) ! FIRST DIMENSION IS INVARIANT ...
           CALL INIT_AXIS(pdax(1)) ! ...
           DO i=2, sax(id)%ndp
              CALL COPY_AXIS(psax(i), sax(sax(id)%dep(i)))
              CALL COPY_AXIS(pdax(i), dax(dax(id)%dep(i)))
           END DO
           !
           CALL NREGRID( (/sax(id)%dat/) ,psax, pdax, pd  &
                         ,(/ RG_INT/) , psovl, pdovl, prcnt)
           CALL NREGRID_STAT(psax, pdax, psovl, pdovl     &
                         ,(/sax(id)%dat/), pd, prcnt)
           !
           CALL COPY_NARRAY(dax(id)%dat, pd(1))
           DEALLOCATE(psax, pdax, psovl, pdovl, STAT=status)
           CALL ERRMSG(substr,status,10)
           NULLIFY(psovl, pdovl)
           DEALLOCATE(pd, STAT=status)
           CALL ERRMSG(substr,status,11)
           NULLIFY(pd)
           DEALLOCATE(prcnt, STAT=status)
           CALL ERRMSG(substr,status,12)
           NULLIFY(prcnt)
           CALL RGMSG(substr, RGMLIC, &
                ' <<< ... END REGRIDDING OF DEPENDENT DIMENSION !')
        END IF

     END DO
     CALL RGMSG(substr, RGMLIC, '... END SCANNING ',ndim, ' DIMENSIONS !')

     ! 3. CALCULATE REGRIDDING ORDER OF DIMENSIONS
     CALL RGMSG(substr, RGMLI, 'PROCESSING ORDER OF ',ndim, ' DIMENSIONS ...')
     DO id=1,s(1)%n
        ! 1st: ALL INDEPENDENT + TO BE REGRIDDED
        vtype = dax(id)%dat%type
        lflag = (sax(id)%ndp == 1).AND.   &    ! source dim. independent
                ! destination dimension available
                ((vtype == VTYPE_REAL).OR.(vtype == VTYPE_DOUBLE)).AND. &
                ! destination dimenions is independent
                (dax(id)%ndp == 1)
        IF (lflag) THEN
           nid = nid + 1
           dio(nid) = id
           CALL RGMSG(substr, RGMLIC, ' ... DIMENSION ',id, ' IS INDEPENDENT')
        END IF
     END DO
     DO id=1,s(1)%n
        ! 2nd: ALL DEPENDENT + TO BE REGRIDDED
        vtype = dax(id)%dat%type
        lflag = ((sax(id)%ndp > 1).OR.      &  !  source dim. indepndent
                 (dax(id)%ndp > 1)).AND.    &
                ! destination dimension available
                ((vtype == VTYPE_REAL).OR.(vtype == VTYPE_DOUBLE))
        IF (lflag) THEN
           nid = nid + 1
           dio(nid) = id
           CALL RGMSG(substr, RGMLIC, ' ... DIMENSION ',id, ' IS DEPENDENT')
        END IF
     END DO
     !
     nr = nid     ! THIS IS THE NUMBER OF DIMs TO REGRID
     !
     DO id=1,s(1)%n
        ! 3rd: REST IS INVARIANT
        vtype = dax(id)%dat%type
        lflag = (vtype /= VTYPE_REAL).AND.(vtype /= VTYPE_DOUBLE)
        IF (lflag) THEN
           nid = nid + 1
           dio(nid) = id
           CALL RGMSG(substr, RGMLIC, ' ... DIMENSION ',id, ' IS INVARIANT')
        END IF
     END DO
     CALL RGMSG(substr, RGMLIC, &
          '... END PROCESSING ORDER OF ',ndim, ' DIMENSIONS !')
     IF (nid < s(1)%n) THEN
        CALL RGMSG(substr, RGMLE, 'UNRECOGNIZED DIMENSION !', .false.)
        CALL RGMSG(substr, RGMLEC, 'CHECK TYPE OF BOUNDS', .false.)
        CALL RGMSG(substr, RGMLEC, '(MUST BE REAL OR DOUBLE PRECISION) !')
     END IF
     IF (nid > s(1)%n) THEN
        CALL RGMSG(substr, RGMLE, 'AMBIGIOUS DIMENSION !')
     END IF

     ! 4. ALLOCATE SPACE FOR REGRIDDED DATA
     ALLOCATE(d(nvar))

     DO k=1, nvar                  ! LOOP OVER VARIABLES
        nid = 1
        d(k)%n = s(k)%n            ! COPY NUMBER OF DIMENSIONS
        d(k)%type = s(k)%type
        ALLOCATE(d(k)%dim(d(k)%n),STAT=status)
        CALL ERRMSG(substr,status,13)
        DO id=1, s(k)%n
           vtype = dax(id)%dat%type
           IF ((vtype == VTYPE_REAL).OR.(vtype == VTYPE_DOUBLE)) THEN
              ! new destination dimension
              ! INTERFACES HAVE +1 ENTRY
              ! 1st DIMENSION MUST BE THE INDEPENDENT !!!
              d(k)%dim(id) = (dax(id)%dat%dim(1) -1)
              nid = nid * d(k)%dim(id)
              ! COPY MODULO INFO
              dax(id)%lm  = sax(id)%lm
           ELSE
              ! keep source dimension
              d(k)%dim(id) = s(k)%dim(id)
              nid = nid * d(k)%dim(id)
              CALL COPY_AXIS(dax(id), sax(id))
           END IF                               ! dest. dim or source dim.
        END DO

        vtype = s(k)%type
        SELECT CASE(vtype)
        CASE(VTYPE_REAL)
           ALLOCATE(d(k)%vr(nid),STAT=status)
           CALL ERRMSG(substr,status,14)
           d(k)%vr(:) = REAL(0.0, SP)
           d(k)%type = VTYPE_REAL
        CASE(VTYPE_DOUBLE)
           ALLOCATE(d(k)%vd(nid),STAT=status)
           CALL ERRMSG(substr,status,15)
           d(k)%vd(:) = REAL(0.0, DP)
           d(k)%type = VTYPE_DOUBLE
        CASE DEFAULT
           !CASE(VTYPE_INT)
           !CASE(VTYPE_CHAR)
           !CASE(VTYPE_BYTE)
           !CASE(VTYPE_UNDEF)
           CALL RGMSG(substr, RGMLE, 'REGRIDDING IS ONLY POSSIBLE FOR', &
                .false.)
           CALL RGMSG(substr, RGMLEC,                                 &
                'DATA OF TYPE REAL OR DOUBLE PRECISION !')
        END SELECT

     END DO  ! LOOP OVER VARIABLES

     ! 5. SOME DIAGNOSTIC OUTPUT
     CALL RGMSG(substr, RGMLVM, '    DIMENSIONS TO REGRID          : ', &
          dio(1:nr),' ')
     !
     ALLOCATE(hio(nr),STAT=status)
     CALL ERRMSG(substr,status,16)
     DO i=1, nr
        hio(i) = s(1)%dim(dio(i))
     END DO
     CALL RGMSG(substr, RGMLVM, '    - LENGTH(S) IN SOURCE         : ', &
          hio,' ')
     !
     DO i=1, nr
        hio(i) = d(1)%dim(dio(i))
     END DO
     CALL RGMSG(substr, RGMLVM, '    - LENGTH(S) IN DEST.          : ', &
          hio,' ')
     DEALLOCATE(hio,STAT=status)
     CALL ERRMSG(substr,status,17)
     !
     CALL RGMSG(substr, RGMLVM, '    INVARIANT DIMENSIONS          : ', &
          dio(nr+1:ndim), ' ')
     CALL RGMSG(substr, RGMLVM, '    NUMBER OF FIELDS              : ', &
          nvar,' ')
     !
     CALL RGMSG(substr, RGMLVM, '    INVARIANT DIMENSION LENGTH(S) : ')
     ALLOCATE(hio(ndim-nr),STAT=status)
     CALL ERRMSG(substr,status,18)
     DO k=1,nvar
        WRITE(hstr, *) '       FIELD ',k,': '
        DO i=nr+1,ndim
           hio(i-nr) = s(k)%dim(dio(i))
        END DO
        CALL RGMSG(substr, RGMLVM, TRIM(hstr), hio,' ')
     END DO
     DEALLOCATE(hio,STAT=status)
     CALL ERRMSG(substr,status,18)

     ! 6. ALLOCATE SPACE FOR DIAGNOSTIC OUTPUT AND INITIALIZE
     ALLOCATE(sovl(s(1)%n),STAT=status)
     CALL ERRMSG(substr,status,19)
     sovl(:) = REAL(0.0, DP)
     ALLOCATE(dovl(d(1)%n),STAT=status)
     CALL ERRMSG(substr,status,20)
     dovl(:) = REAL(0.0, DP)
     ALLOCATE(rcnt(ndim),STAT=status)
     CALL ERRMSG(substr,status,21)
     rcnt(:) = 0

     CALL RGMSG(substr, RGMLVL, '    STARTING NREGRID ... ')

  END IF ! (A) ONLY TO BE DONE ONCE AT FIRST STEP

  ! ....................................................................
  ! SET CURRENT DIMENSION
  nid = dio(m)
  rcnt(m) = rcnt(m) + 1

  ! ....................................................................
  ! (B) CONDITION FOR END OF RECURSION
  IF (m == SIZE(dio)) THEN   ! LAST (INVARIANT !) DIMENSION REACHED
     zso = REAL(0.0, DP)
     zdo = REAL(0.0, DP)
     DO k=1, nvar    ! LOOP OVER VARIABLES
        SELECT CASE(RG_TYPE(k))
        CASE(RG_EXT)
           DO j=1, s(k)%dim(nid)
              dvec(nid) = j
              svec(nid) = j     ! INVARIANT !!!
              ! CALCULATE DEST. DATA POINT
              vo = POSITION(d(k)%dim,dvec)   ! POSITION IN DEST.
              vl = POSITION(s(k)%dim,svec)   ! POSITION IN SOURCE
              vtype = s(k)%type
              SELECT CASE(vtype)
              CASE(VTYPE_REAL)
                 d(k)%vr(vo) = d(k)%vr(vo) + s(k)%vr(vl) * REAL(sf0, SP)
              CASE(VTYPE_DOUBLE)
                 d(k)%vd(vo) = d(k)%vd(vo) + s(k)%vd(vl) * sf0
              CASE DEFAULT
                 !CASE(VTYPE_INT)
                 !CASE(VTYPE_CHAR)
                 !CASE(VTYPE_BYTE)
                 !CASE(VTYPE_UNDEF)
                 CALL RGMSG(substr, RGMLE, &
                      'REGRIDDING IS ONLY POSSIBLE FOR DATA', .false.)
                 CALL RGMSG(substr, RGMLEC, &
                      'OF TYPE REAL OR DOUBLE PRECISION !')
              END SELECT
              ! OVERLAP IS 1 FOR INVARIANT DIMENSIONS
              ! INTERVAL LENGTH, NUMBER OF VARIABLES
              zso = zso + REAL(1.0, DP)/ &
                   (REAL(s(k)%dim(nid), DP)*REAL(nvar, DP))
              zdo = zdo + REAL(1.0, DP)/ &
                   (REAL(s(k)%dim(nid), DP)*REAL(nvar, DP))
           END DO
           ! RESET FOR NEXT/PREVIOUS RECURSION STEP
           dvec(nid) = 1
           svec(nid) = 1
        CASE(RG_INT, RG_IDX, RG_IXF)
           !CASE(RG_IDX, RG_IXF)
           ! NOTE: IDX VAR WAS CONVERTET TO INDEX-FRACTION
           DO j=1, s(k)%dim(nid)
              dvec(nid) = j
              svec(nid) = j     ! INVARIANT !!!
              ! CALCULATE DEST. DATA POINT
              vo = POSITION(d(k)%dim,dvec)   ! POSITION IN DEST.
              vl = POSITION(s(k)%dim,svec)   ! POSITION IN SOURCE
              vtype = s(k)%type
              SELECT CASE(vtype)
              CASE(VTYPE_REAL)
                d(k)%vr(vo) = d(k)%vr(vo) + s(k)%vr(vl) * REAL(df0, SP)
              CASE(VTYPE_DOUBLE)
                 d(k)%vd(vo) = d(k)%vd(vo) + s(k)%vd(vl) * df0
              CASE(VTYPE_INT)
                d(k)%vd(vo) = d(k)%vd(vo) + REAL(s(k)%vi(vl),DP) * df0
              CASE(VTYPE_BYTE)
                 d(k)%vd(vo) = d(k)%vd(vo) + REAL(s(k)%vb(vl),DP) * df0
              CASE DEFAULT
                 !CASE(VTYPE_CHAR)
                 !CASE(VTYPE_UNDEF)
                 CALL RGMSG(substr, RGMLE, &
                      'REGRIDDING IS ONLY POSSIBLE FOR DATA', .false.)
                 CALL RGMSG(substr, RGMLEC, &
                      'OF TYPE REAL OR DOUBLE PRECISION !')
              END SELECT
              ! OVERLAP IS 1 FOR INVARIANT DIMENSIONS:
              ! INTERVAL LENGTH, NUMBER OF VARIABLES
              zso = zso + REAL(1.0, DP)/ &
                   (REAL(s(k)%dim(nid), DP)*REAL(nvar, DP))
              zdo = zdo + REAL(1.0, DP)/ &
                   (REAL(s(k)%dim(nid), DP)*REAL(nvar, DP))
           END DO
           ! RESET FOR NEXT/PREVIOUS RECURSION STEP
           dvec(nid) = 1
           svec(nid) = 1
        !CASE DEFAULT
         END SELECT

     END DO  ! LOOP OVER VARIABLES

     ! SUM OVERLAP
     sovl(nid) = sovl(nid) + zso
     dovl(nid) = dovl(nid) + zdo

     ! CLEAN UP
     DEALLOCATE(dio, STAT=status)
     CALL ERRMSG(substr,status,22)
     DEALLOCATE(svec, dvec, STAT=status)
     CALL ERRMSG(substr,status,23)

     RETURN        ! DONE ... ... END RECURSION

  END IF ! (B) CONDITION FOR END OF RECURSION
  ! ....................................................................

  ! ....................................................................
  ! (C) RECURSION STEP
  IF (m <= nr) THEN           ! (C1): DIMENSION TO REGRID

     ! GET SOURCE BOUNDS
     ! 1st dimension is independent -> length of bounds vector
     vl = sax(nid)%dat%dim(1)        ! length is length of 1st (indep.) dim.
     ALLOCATE(sb(vl),STAT=status)
     CALL ERRMSG(substr,status,24)

     ALLOCATE(vec(sax(nid)%ndp),STAT=status)  ! dependency element vector
     CALL ERRMSG(substr,status,25)
     vec(1) = 1         ! start at 1st position of independent dim.
     DO i=2,sax(nid)%ndp        ! loop over dependent dims
        vec(i) = svec(sax(nid)%dep(i)) ! ... actual counter of dep. dims
     END DO

     ! CALCULATE OFFSET ...
     vo = POSITION(sax(nid)%dat%dim, vec)
     ! ... AND GET BOUNDS
     vtype = sax(nid)%dat%type

     SELECT CASE(vtype)
     CASE(VTYPE_REAL)
        sb(:) = sax(nid)%dat%vr(vo:(vo+vl-1))
     CASE(VTYPE_DOUBLE)
        sb(:) = sax(nid)%dat%vd(vo:(vo+vl-1))
     CASE DEFAULT
        !CASE(VTYPE_INT)
        !CASE(VTYPE_CHAR)
        !CASE(VTYPE_BYTE)
        !CASE(VTYPE_UNDEF)
        CALL RGMSG(substr, RGMLE, 'SOURCE BOUNDS MUST BE OF TYPE', .false.)
        CALL RGMSG(substr, RGMLEC, 'REAL OR DOUBLE PRECISION !')
     END SELECT

     !
     DEALLOCATE(vec,STAT=status)
     CALL ERRMSG(substr,status,26)

     ! GET DEST. BOUNDS
     ! 1st dimension is independent -> length of bounds vector
     vl = dax(nid)%dat%dim(1)        ! length is length of 1st (indep.) dim.
     ALLOCATE(db(vl),STAT=status)
     CALL ERRMSG(substr,status,27)

     ALLOCATE(vec(dax(nid)%ndp),STAT=status)  ! dependency element vector
     CALL ERRMSG(substr,status,28)
     vec(1) = 1         ! start at 1st position of independent dim.
     DO i=2,dax(nid)%ndp        ! loop over dependent dims
        vec(i) = dvec(dax(nid)%dep(i)) ! ... actual counter of dep. dims
     END DO

     ! CALCULATE OFFSET ...
     vo = POSITION(dax(nid)%dat%dim, vec)
     ! ... AND GET BOUNDS
     vtype = dax(nid)%dat%type
     SELECT CASE(vtype)
     CASE(VTYPE_REAL)
        db(:) = dax(nid)%dat%vr(vo:(vo+vl-1))
     CASE(VTYPE_DOUBLE)
        db(:) = dax(nid)%dat%vd(vo:(vo+vl-1))
     CASE DEFAULT
        !CASE(VTYPE_INT)
        !CASE(VTYPE_CHAR)
        !CASE(VTYPE_BYTE)
        !CASE(VTYPE_UNDEF)
        CALL RGMSG(substr, RGMLE, 'DEST. BOUNDS MUST BE OF TYPE', .false.)
        CALL RGMSG(substr, RGMLEC, 'REAL OR DOUBLE PRECISION !')
     END SELECT

     !
     DEALLOCATE(vec,STAT=status)
     CALL ERRMSG(substr,status,29)

    ! REGRID ALONG THIS DIMENSION
     ! ALLOCATE SPACE FOR MATRICES
     ns = SIZE(sb)-1
     nd = SIZE(db)-1
     ALLOCATE(ms(ns,nd), STAT=status)
     CALL ERRMSG(substr,status,30)
     ALLOCATE(md(nd,ns), STAT=status)
     CALL ERRMSG(substr,status,31)
     CALL OVL_1D(sb, db, ms, md, sax(nid)%lm)

     ! SAVE OVERLAP OF NEXT DIMENSION CALCULATED SO FAR ...
     zso = sovl(dio(m+1))
     zdo = dovl(dio(m+1))
     ! ... AND RESET VECTOR ELEMENT TO ZERO
     sovl(dio(m+1)) = REAL(0.0, DP)
     dovl(dio(m+1)) = REAL(0.0, DP)
     novl = 0
     ! LOOP OVER ALL MATRIX ELEMENTS AND REGRID ALONG NEXT DIMENSION ...
     ! ... IF MATRIX ELEMENT IS NON-ZERO
     !

     DO i=1, ns
        svec(nid) = i
        DO j=1, nd
           dvec(nid) = j
           sf = sf0*ms(svec(nid),dvec(nid))
           df = df0*md(dvec(nid),svec(nid))
           IF (sf > 0.0) THEN    ! NON-ZERO
              novl = novl + 1
              ! CALCULATE TOTAL OVERLAP FRACTION FOR THIS DIMENSION
              sovl(nid) = sovl(nid) +                 &
                   ms(svec(nid),dvec(nid))            &
                   *(sb(svec(nid)+1)-sb(svec(nid)))   &
                   /(sb(ns+1)-sb(1))
              dovl(nid) = dovl(nid) +                 &
                   md(dvec(nid),svec(nid))            &
                   *(db(dvec(nid)+1)-db(dvec(nid)))   &
                   /(db(nd+1)-db(1))
              ! REGRID NEXT DIMENSION
                 CALL NREGRID(s, sax, dax, d, RG_TYPE, sovl, dovl, rcnt  &
                      ,m+1, nr, svec, dvec                               &
                      ,dio, sf, df)
           END IF
        END DO
     END DO

    ! AVERAGE OVERLAP FRACTIONS:
     ! RESTORE OLD OVERLAP FRACTION AND ADD NEW
     IF (novl /= 0) THEN
        sovl(dio(m+1)) = zso + sovl(dio(m+1))/REAL(novl, DP)
        dovl(dio(m+1)) = zdo + dovl(dio(m+1))/REAL(novl, DP)
     ELSE
        sovl(dio(m+1)) = zso
        dovl(dio(m+1)) = zdo
     END IF

     ! RESET FOR NEXT/PREVIOUS RECURSION STEP
     svec(nid) = 1
     dvec(nid) = 1

     ! CLEAN
     DEALLOCATE(ms, md, stat=status)
     CALL ERRMSG(substr,status,32)
     DEALLOCATE(sb, db, stat=status)
     CALL ERRMSG(substr,status,33)

  ELSE  ! (C2): INVARIANT DIMENSION: ( m > nr )  .......................

     ! SAVE OVERLAP OF NEXT DIMENSION CALCULATED SO FAR ...
     zso = sovl(dio(m+1))
     zdo = dovl(dio(m+1))
     ! ... AND RESET VECTOR ELEMENT TO ZERO
     sovl(dio(m+1)) = REAL(0.0, DP)
     dovl(dio(m+1)) = REAL(0.0, DP)
     ! LOOP OVER ALL DIAGONAL MATRIX ELEMENTS ...
     ! ... AND REGRID NEXT DIMENSION
     DO j=1, s(1)%dim(nid)
        dvec(nid) = j
        svec(nid) = j     ! INVARIANT !!!
        ! REGRID NEXT DIMENSION
        CALL NREGRID(s, sax, dax, d, RG_TYPE, sovl, dovl, rcnt  &
             ,m+1, nr, svec, dvec                               &
             ,dio, sf0, df0 )
        ! OVERLAP IS 1 FOR INVARIANT DIMENSIONS
        sovl(nid) = sovl(nid)+REAL(1.0, DP)/REAL(s(1)%dim(nid), DP)
        dovl(nid) = dovl(nid)+REAL(1.0, DP)/REAL(s(1)%dim(nid), DP)
     END DO
     
     ! AVERAGE OVERLAP FRACTIONS:
     ! RESTORE OLD OVERLAP FRACTION AND ADD NEW
     sovl(dio(m+1)) = zso + sovl(dio(m+1))/REAL(s(1)%dim(nid), DP)
     dovl(dio(m+1)) = zdo + dovl(dio(m+1))/REAL(s(1)%dim(nid), DP)

     ! RESET FOR NEXT/PREVIOUS RECURSION STEP
     svec(nid) = 1
     dvec(nid) = 1


  END IF  ! (C) RECURSION STEP
  ! ....................................................................

  ! CLEAN UP
  DEALLOCATE(dio, STAT=status)
  CALL ERRMSG(substr,status,34)
  DEALLOCATE(svec, dvec, STAT=status)
  CALL ERRMSG(substr,status,35)

END SUBROUTINE NREGRID
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE NREGRID_STAT(sax, dax, sovl, dovl, nai, nao, rcnt)

  IMPLICIT NONE

  ! I/O
  TYPE (t_axis),   DIMENSION(:), INTENT(IN) :: sax, dax
  REAL (DP),       DIMENSION(:), INTENT(IN) :: sovl, dovl
  TYPE (t_narray), DIMENSION(:), INTENT(IN) :: nai, nao
  INTEGER,         DIMENSION(:), INTENT(IN) :: rcnt

  ! LOCAL
  INTEGER   :: i
  INTEGER   :: vtype
  REAL (DP) :: div

  WRITE(*,*) '    NREGRID STATISTICS:'
  WRITE(*,*) '    .......................................................'
  WRITE(*,*) '    NO. OF RECURSION LEVELS   : ',SIZE(rcnt)
  WRITE(*,*) '    RECURSION LEVELS PROCESSED: ',rcnt
  WRITE(*,*) ' '
  DO i=1, SIZE(sovl)
     !
     IF (i > 1) THEN
        IF (rcnt(i-1) /= 0) THEN
           div = REAL(rcnt(i-1), DP)
        ELSE
           div = REAL(1.0, DP)
        END IF
     ELSE
        div = REAL(1.0, DP)
     END IF
     !
     IF (sax(i)%lm.OR.dax(i)%lm) THEN
        WRITE(*,*) '    DIMENSION ',i,': (MODULO)'
     ELSE
        WRITE(*,*) '    DIMENSION ',i,':'
     END IF
     vtype = sax(i)%dat%type
     SELECT CASE(vtype)
     CASE(VTYPE_REAL)
        WRITE(*,*) '      SOURCE: ',MINVAL(sax(i)%dat%vr),&
             MAXVAL(sax(i)%dat%vr)
     CASE(VTYPE_DOUBLE)
        WRITE(*,*) '      SOURCE: ',MINVAL(sax(i)%dat%vd),&
             MAXVAL(sax(i)%dat%vd)
     CASE(VTYPE_INT)
        WRITE(*,*) '      SOURCE: ',MINVAL(sax(i)%dat%vi),&
             MAXVAL(sax(i)%dat%vi)
     CASE(VTYPE_BYTE)
        WRITE(*,*) '      SOURCE: ',MINVAL(sax(i)%dat%vb),&
             MAXVAL(sax(i)%dat%vb)
     CASE(VTYPE_CHAR)
        WRITE(*,*) '      SOURCE:  CHAR NOT SUPPORTED'
     CASE DEFAULT
        WRITE(*,*) '      SOURCE:  <UNDEFINED>'
     END SELECT
     WRITE(*,*) '      SOURCE REGION COVERED ON AVERAGE: ',sovl(i)/div
     !
     vtype = dax(i)%dat%type
     SELECT CASE(vtype)
     CASE(VTYPE_REAL)
        WRITE(*,*) '      DEST. : ',MINVAL(dax(i)%dat%vr),&
             MAXVAL(dax(i)%dat%vr)
     CASE(VTYPE_DOUBLE)
        WRITE(*,*) '      DEST. : ',MINVAL(dax(i)%dat%vd),&
             MAXVAL(dax(i)%dat%vd)
     CASE(VTYPE_INT)
        WRITE(*,*) '      DEST. : ',MINVAL(dax(i)%dat%vi),&
             MAXVAL(dax(i)%dat%vi)
     CASE(VTYPE_BYTE)
        WRITE(*,*) '      DEST. : ',MINVAL(dax(i)%dat%vb),&
             MAXVAL(dax(i)%dat%vb)
     CASE(VTYPE_CHAR)
        WRITE(*,*) '      DEST. :  CHAR NOT SUPPORTED'
     CASE DEFAULT
        WRITE(*,*) '      DEST. :  <UNDEFINED>'
     END SELECT
     !
     WRITE(*,*) '      DEST.  REGION COVERED ON AVERAGE: ',dovl(i)/div
  END DO
  !
  WRITE(*,*) ' '
  WRITE(*,*) '    VARIABLE RANGES: '
  DO i=1, SIZE(nai)
     vtype = nai(i)%type
     SELECT CASE(vtype)
     CASE(VTYPE_REAL)
        WRITE(*,*) '    (',i,'): <REAL>'
        WRITE(*,*) '      SOURCE: ',MINVAL(nai(i)%vr),MAXVAL(nai(i)%vr)
        WRITE(*,*) '      DEST. : ',MINVAL(nao(i)%vr),MAXVAL(nao(i)%vr)
     CASE(VTYPE_DOUBLE)
        WRITE(*,*) '    (',i,'): <DOUBLE PRECISION>'
        WRITE(*,*) '      SOURCE: ',MINVAL(nai(i)%vd),MAXVAL(nai(i)%vd)
        WRITE(*,*) '      DEST. : ',MINVAL(nao(i)%vd),MAXVAL(nao(i)%vd)
     CASE(VTYPE_INT)
        WRITE(*,*) '    (',i,'): <INTEGER>'
        WRITE(*,*) '      SOURCE: ',MINVAL(nai(i)%vi),MAXVAL(nai(i)%vi)
        WRITE(*,*) '      DEST. : ',MINVAL(nao(i)%vi),MAXVAL(nao(i)%vi)
     CASE(VTYPE_BYTE)
        WRITE(*,*) '    (',i,'): <BYTE>'
        WRITE(*,*) '      SOURCE: ',MINVAL(nai(i)%vb),MAXVAL(nai(i)%vb)
        WRITE(*,*) '      DEST. : ',MINVAL(nao(i)%vb),MAXVAL(nao(i)%vb)
     CASE(VTYPE_CHAR)
        WRITE(*,*) '    (',i,'): <CHAR>'
        WRITE(*,*) '      SOURCE:  CHAR NOT SUPPORTED'
        WRITE(*,*) '      DEST. :  CHAR NOT SUPPORTED'
     CASE DEFAULT
        WRITE(*,*) '    (',i,'): <UNDEFINED>'
        WRITE(*,*) '      SOURCE:  <UNDEFINED>'
        WRITE(*,*) '      DEST. :  <UNDEFINED>'
     END SELECT
  END DO
  !
  WRITE(*,*) '    .......................................................'

END SUBROUTINE NREGRID_STAT
! ------------------------------------------------------------------

SUBROUTINE SNREGRID(status, gin, gout, din, dout, lext)

  IMPLICIT NONE

  ! I/O
  INTEGER,                INTENT(OUT)   :: status
  REAL(DP), DIMENSION(:), INTENT(IN)    :: gin    ! bounds of input grid
  REAL(DP), DIMENSION(:), INTENT(IN)    :: gout   ! bounds of output grid
  REAL(DP), DIMENSION(:), INTENT(IN)    :: din    ! data on input grid
  REAL(DP), DIMENSION(:), INTENT(INOUT) :: dout   ! data on output grid
  LOGICAL,                INTENT(IN)    :: lext   ! extensive (T), intensive (F)

  ! LOCAL
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: mfs, mfd
  INTEGER                               :: n, m
  INTEGER                               :: i, j

  ! DATA sizes
  n = SIZE(gin)  - 1
  m = SIZE(gout) - 1

  ! CHECK FOR CONSISTENCY
  IF (SIZE(din) /= n) THEN
     status = 3
     RETURN
  ENDIF
  IF (SIZE(dout) /= m) THEN
     status = 4
     RETURN
  ENDIF

  ! ALLOCATE WORKSPACE
  ALLOCATE(mfs(n,m))
  ALLOCATE(mfd(m,n))

  CALL OVL_1D(gin, gout, mfs, mfd, .FALSE.)

  IF (lext) THEN
     ! extensive
     DO j=1,m
        dout(j) = 0.0_dp
        DO i=1, n
           dout(j) = dout(j) + din(i)*mfs(i,j)
        END DO
     END DO
  ELSE
     ! intensive
     DO j=1,m
        DO i=1,n
           dout(j) = dout(j) + din(i)*mfd(j,i)
        END DO
     END DO
  ENDIF

  ! CLEAN
  DEALLOCATE(mfs)
  DEALLOCATE(mfd)

  status = 0

END SUBROUTINE SNREGRID

! ******************************************************************
END MODULE MESSY_MAIN_GRID_TRAFO_NRGD_BASE
! ******************************************************************
