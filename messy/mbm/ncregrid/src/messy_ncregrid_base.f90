#ifndef _SX

! ******************************************************************
! ------------------------------------------------------------------
MODULE MESSY_NCREGRID_BASE
! ------------------------------------------------------------------
! Author: Patrick Joeckel, MPICH, Mainz, June 2002
! ******************************************************************

#if defined(MESSY)
  USE messy_main_constants_mem,  ONLY: SP, DP, I4, I8
#else
  USE typeSizes,     ONLY:   SP => FourByteReal  &
                           , DP => EightByteReal &
                           , I4 => FourByteInt   &
                           , I8 => EightByteInt
#endif

  IMPLICIT NONE

  PUBLIC

  INTRINSIC :: ASSOCIATED, PRESENT, PRODUCT, REAL, SIZE &
             , INT, ABS, DBLE, MAX, MIN, SIGN, TRIM     &
             , MAXVAL, MINVAL, IAND, NULL

  PRIVATE   :: ASSOCIATED, PRESENT, PRODUCT, REAL, SIZE &
             , INT, ABS, DBLE, MAX, MIN, SIGN, TRIM     &
             , MAXVAL, MINVAL, IAND, NULL

  ! REGRIDDING VARIABLE TYPES
  INTEGER, PARAMETER :: RG_INT = 1
  INTEGER, PARAMETER :: RG_EXT = 2
  INTEGER, PARAMETER :: RG_IDX = 3
  INTEGER, PARAMETER :: RG_IXF = 4

  ! REGRIDDING DATA TYPES
  INTEGER, PARAMETER :: VTYPE_UNDEF  = 0
  INTEGER, PARAMETER :: VTYPE_INT    = 1
  INTEGER, PARAMETER :: VTYPE_REAL   = 2
  INTEGER, PARAMETER :: VTYPE_DOUBLE = 3
  INTEGER, PARAMETER :: VTYPE_BYTE   = 4
  INTEGER, PARAMETER :: VTYPE_CHAR   = 5

  ! MESSAGE LEVEL TYPES
  INTEGER, PARAMETER :: RGMLE   = 0   ! ERROR
  INTEGER, PARAMETER :: RGMLEC  = 1   ! ERROR CONTINUED
  INTEGER, PARAMETER :: RGMLVL  = 2   ! LITTLE VERBOSE
  INTEGER, PARAMETER :: RGMLVLC = 3   ! LITTLE VERBOSE CONTINUED
  INTEGER, PARAMETER :: RGMLW   = 4   ! WARNING
  INTEGER, PARAMETER :: RGMLWC  = 5   ! WARNING CONTINUED
  INTEGER, PARAMETER :: RGMLVM  = 6   ! MEDIUM VERBOSE
  INTEGER, PARAMETER :: RGMLVMC = 7   ! MEDIUM VERBOSE CONTINUED
  INTEGER, PARAMETER :: RGMLI   = 8   ! INFO
  INTEGER, PARAMETER :: RGMLIC  = 9   ! INFO CONTINUED

  ! MESSAGE OUTPUT LEVEL
  INTEGER, PARAMETER :: MSGMODE_S  =  0  ! SILENT
  INTEGER, PARAMETER :: MSGMODE_E  =  1  ! ERROR MESSAGES
  INTEGER, PARAMETER :: MSGMODE_VL =  2  ! LITTLE VERBOSE
  INTEGER, PARAMETER :: MSGMODE_W  =  4  ! WARNING MESSAGES
  INTEGER, PARAMETER :: MSGMODE_VM =  8  ! MEDIUM VERBOSE
  INTEGER, PARAMETER :: MSGMODE_I  = 16  ! INFO MESSAGES
  INTEGER, SAVE      :: MSGMODE = MSGMODE_S + MSGMODE_E + MSGMODE_VL &
                                + MSGMODE_W + MSGMODE_VM + MSGMODE_I

  ! mz_pj_20080801+
  INTEGER, SAVE :: MY_RANK = -1
  ! mz_pj_20080801-

  TYPE narray
     ! n-dimenional array as 1D (LINEAR) array (REAL)
     INTEGER                              :: type = VTYPE_UNDEF
     INTEGER                              :: n = 0  ! number of dimensions
     INTEGER     , DIMENSION(:), POINTER  :: dim => NULL()! dim. vector
     REAL (SP)   , DIMENSION(:), POINTER  :: vr => NULL() ! real values
     REAL (DP)   , DIMENSION(:), POINTER  :: vd => NULL() ! double values
     INTEGER (I8), DIMENSION(:), POINTER  :: vi => NULL() ! integer values
     INTEGER (I4), DIMENSION(:), POINTER  :: vb => NULL() ! byte values
     CHARACTER,    DIMENSION(:), POINTER  :: vc => NULL() ! char. values
  END TYPE narray

  TYPE axis
     ! hyper-axis (for curvilinear coordinates)
     TYPE (narray)                        :: dat  ! interface bounds
     LOGICAL                              :: lm = .false. ! modulo axis ?
     INTEGER (I4)                         :: ndp = 0  ! number of dependencies
     ! FIRST IN LIST MUST BE INDEPENDENT !!!
     ! E.G. IF DIM 3 DEPENDS ON DIM 1,2, and 5
     ! dep = (/ 3,1,2,5 /)
     INTEGER     , DIMENSION(:), POINTER  :: dep => NULL()  ! dependencies
  END TYPE axis

!! NOTE: DOES NOT WORK PROPERLY FOR SOME COMPILERS ...
!  INTERFACE ASSIGNMENT (=)
!     MODULE PROCEDURE COPY_NARRAY
!     MODULE PROCEDURE COPY_AXIS
!  END INTERFACE

INTERFACE QSORT
   ! THE GOOD OLD (RECURSIVE) QUICKSORT ALGORITHM FOR LINEAR ARRAYS
   MODULE PROCEDURE QSORT_I    ! INTEGER
#if !(defined(__SX__))
   MODULE PROCEDURE QSORT_B    ! BYTE
   MODULE PROCEDURE QSORT_R    ! REAL
#endif
   MODULE PROCEDURE QSORT_D    ! DOUBLE PRECISION
END INTERFACE

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

INTERFACE RGMSG                    ! MESSAGE OUTPUT
   MODULE PROCEDURE RGMSG_C
   MODULE PROCEDURE RGMSG_I
   MODULE PROCEDURE RGMSG_IA
   MODULE PROCEDURE RGMSG_R
   MODULE PROCEDURE RGMSG_D ! op_mr_20170321
END INTERFACE

CONTAINS

! ------------------------------------------------------------------
SUBROUTINE INIT_NARRAY(na, n, dim, qtype)

  IMPLICIT NONE

  ! I/O
  TYPE (narray),         INTENT(INOUT)          :: na
  INTEGER,               INTENT(IN),   OPTIONAL :: n
  INTEGER, DIMENSION(:), INTENT(IN),   OPTIONAL :: dim
  INTEGER,                             OPTIONAL :: qtype

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'INIT_NARRAY'
  INTEGER                     :: status
  INTEGER                     :: len

  na%type = VTYPE_UNDEF

  na%n = 0
  IF (ASSOCIATED(na%dim)) THEN
     IF (SIZE(na%dim) > 0) THEN
        DEALLOCATE(na%dim, STAT=status)
        CALL ERRMSG(substr,status,1)
     END IF
     NULLIFY(na%dim)
  END IF

  IF (ASSOCIATED(na%vi)) THEN
     IF (SIZE(na%vi) > 0) THEN
        DEALLOCATE(na%vi, STAT=status)
        CALL ERRMSG(substr,status,2)
     END IF
     NULLIFY(na%vi)
  END IF

  IF (ASSOCIATED(na%vr)) THEN
     IF (SIZE(na%vr) > 0) THEN
        DEALLOCATE(na%vr, STAT=status)
        CALL ERRMSG(substr,status,3)
     END IF
     NULLIFY(na%vr)
  END IF

  IF (ASSOCIATED(na%vd)) THEN
     IF (SIZE(na%vd) > 0) THEN
        DEALLOCATE(na%vd, STAT=status)
        CALL ERRMSG(substr,status,4)
     END IF
     NULLIFY(na%vd)
  END IF

  IF (ASSOCIATED(na%vc)) THEN
     IF (SIZE(na%vc) > 0) THEN
        DEALLOCATE(na%vc, STAT=status)
        CALL ERRMSG(substr,status,5)
     END IF
     NULLIFY(na%vc)
  END IF

  IF (ASSOCIATED(na%vb)) THEN
     IF (SIZE(na%vb) > 0) THEN
        DEALLOCATE(na%vb, STAT=status)
        CALL ERRMSG(substr,status,6)
     END IF
     NULLIFY(na%vb)
  END IF  

  IF (PRESENT(dim) .AND. PRESENT(n)) THEN
     IF (SIZE(dim) /= n) THEN
        CALL RGMSG(substr, RGMLE,                &
             'NUMBER OF DIMENSIONS MISMATCH !' )
     END IF
  END IF

  len = 0
  IF (PRESENT(n)) THEN
     na%n = n
     ALLOCATE(na%dim(n), STAT=status)
     CALL ERRMSG(substr,status,7)
     na%dim(:) = 0
  END IF

  IF (PRESENT(dim)) THEN
     IF (.NOT. ASSOCIATED(na%dim)) THEN
        na%n = SIZE(dim)
        ALLOCATE(na%dim(na%n), STAT=status)
        CALL ERRMSG(substr,status,8)
     END IF
     na%dim(:) = dim(:)
     len = PRODUCT(na%dim(:))
  END IF

  IF (PRESENT(qtype).AND.(len > 0)) THEN
     SELECT CASE(qtype)
     CASE(VTYPE_INT)
        na%type = VTYPE_INT
        ALLOCATE(na%vi(len), STAT=status)
        CALL ERRMSG(substr,status,9)
        na%vi(:) = 0
     CASE(VTYPE_REAL)
        na%type = VTYPE_REAL
        ALLOCATE(na%vr(len), STAT=status)
        CALL ERRMSG(substr,status,10)
        na%vr(:) = 0.0
     CASE(VTYPE_DOUBLE)
        na%type = VTYPE_DOUBLE
        ALLOCATE(na%vd(len), STAT=status)
        CALL ERRMSG(substr,status,11)
        na%vd(:) = REAL(0.0, DP)
     CASE(VTYPE_CHAR)
        na%type = VTYPE_CHAR
        ALLOCATE(na%vc(len), STAT=status)
        CALL ERRMSG(substr,status,12)
        na%vc(:) = ' '
     CASE(VTYPE_BYTE)
        na%type = VTYPE_BYTE
        ALLOCATE(na%vb(len), STAT=status)
        CALL ERRMSG(substr,status,13)
        na%vi(:) = 0
     CASE(VTYPE_UNDEF)
        na%type = VTYPE_UNDEF
        ! DO NOTHING, KEEP UNDEFINED
     CASE DEFAULT
        CALL RGMSG(substr, RGMLE,                           &
             'REQUESTED TYPE FOR N-ARRAY IS UNRECOGNIZED !' )
     END SELECT
  END IF

END SUBROUTINE INIT_NARRAY
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE INIT_AXIS(a)

  IMPLICIT NONE

  ! I/O
  TYPE (axis), INTENT(INOUT) :: a
  
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
SUBROUTINE COPY_NARRAY(d, s)

  IMPLICIT NONE

  ! I/O
  TYPE (narray), INTENT(OUT) :: d
  TYPE (narray), INTENT(IN)  :: s
  
  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'COPY_NARRAY'
  INTEGER :: n
  INTEGER :: vtype
  INTEGER :: status

  ! INIT
  n = 0

  d%n = s%n
  IF (ASSOCIATED(s%dim)) THEN
     ALLOCATE(d%dim(d%n),STAT=status)
     CALL ERRMSG(substr,status,1)
     d%dim(:) = s%dim(:)
  END IF

  vtype = s%type
  SELECT CASE(vtype)
  CASE(VTYPE_INT)
     d%type = VTYPE_INT
     n = SIZE(s%vi)
     IF (n > 0) THEN
        ALLOCATE(d%vi(n),STAT=status)
        CALL ERRMSG(substr,status,2)
        d%vi(:) = s%vi(:) 
     END IF
  CASE(VTYPE_REAL)
     d%type = VTYPE_REAL
     n = SIZE(s%vr)
     IF (n > 0) THEN
        ALLOCATE(d%vr(n),STAT=status)
        CALL ERRMSG(substr,status,3)
        d%vr(:) = s%vr(:)
     END IF
  CASE(VTYPE_DOUBLE)
     d%type = VTYPE_DOUBLE
     n = SIZE(s%vd)
     IF (n > 0) THEN
        ALLOCATE(d%vd(n),STAT=status)
        CALL ERRMSG(substr,status,4)
        d%vd(:) = s%vd(:)
     END IF
  CASE(VTYPE_CHAR)
     d%type = VTYPE_CHAR
     n = SIZE(s%vc)
     IF (n > 0) THEN
        ALLOCATE(d%vc(n),STAT=status)
        CALL ERRMSG(substr,status,5)
        d%vc(:) = s%vc(:) 
     END IF
  CASE(VTYPE_BYTE)
     d%type = VTYPE_BYTE
     n = SIZE(s%vb)
     IF (n > 0) THEN
        ALLOCATE(d%vb(n),STAT=status)
        CALL ERRMSG(substr,status,6)
        d%vb(:) = s%vb(:)
     END IF
  CASE(VTYPE_UNDEF)
     d%type = VTYPE_UNDEF
     ! DO NOTHING, EMPTY N-ARRAY IS COPIED
  CASE DEFAULT
     CALL RGMSG(substr, RGMLE,                              &
          'N-ARRAY OF UNRECOGNIZED TYPE CANNOT BE COPIED !' )          
  END SELECT

END SUBROUTINE COPY_NARRAY
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE COPY_AXIS(d, s)

  IMPLICIT NONE

  ! I/O
  TYPE (axis), INTENT(OUT) :: d
  TYPE (axis), INTENT(IN)  :: s

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
SUBROUTINE SORT_NARRAY(na, nx, reverse)

  IMPLICIT NONE

  ! I/O
  TYPE (narray), INTENT(INOUT)          :: na
  TYPE (narray), INTENT(INOUT)          :: nx
  LOGICAL      , INTENT(IN)   ,OPTIONAL :: reverse

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'SORT_NARRAY'
  INTEGER       :: vtype
  LOGICAL       :: lrev
  INTEGER       :: n, i
  TYPE (narray) :: nh

  IF (PRESENT(reverse)) THEN
     lrev = reverse
  ELSE
     lrev = .false.  ! DEFAULT
  END IF

  IF (na%n == 0) THEN
     CALL RGMSG(substr, RGMLW, 'EMPTY ARRAY ! NOTHING TO DO !')
     RETURN
  END IF

  IF (na%n > 1) THEN
     CALL RGMSG(substr, RGMLW, &
          'SORTING A ',na%n,'-DIMENSIONAL N-ARRAY AS LINEAR ARRAY !')
  END IF

  IF (lrev) THEN

     vtype = nx%type
     IF (vtype /= VTYPE_INT) THEN
        CALL RGMSG(substr, RGMLE, 'INDEX N-ARRAY MUST BE OF TYPE INTEGER !')
     END IF
     n = SIZE(nx%vi)
     CALL COPY_NARRAY(nh, na)   ! COPY TO BE RE-ORDERED
     vtype = na%type
     SELECT CASE(vtype)
     CASE(VTYPE_REAL)
        DO i=1,n
           na%vr(nx%vi(i)) = nh%vr(i)
        END DO
     CASE(VTYPE_DOUBLE)
        DO i=1,n
           na%vd(nx%vi(i)) = nh%vd(i)
        END DO
     CASE(VTYPE_INT)
        DO i=1,n
           na%vi(nx%vi(i)) = nh%vi(i)
        END DO        
     CASE(VTYPE_BYTE)
        DO i=1,n
           na%vb(nx%vi(i)) = nh%vb(i)
        END DO        
     CASE(VTYPE_CHAR)
     CALL RGMSG(substr,RGMLE,'UN-SORTING OF TYPE CHAR IS NOT IMPLEMENTED !')
     CASE(VTYPE_UNDEF)
     CALL RGMSG(substr,RGMLE,'ARRAY OF UNDEFINED TYPE CANNOT BE UN-SORTED !')
     CASE DEFAULT
   CALL RGMSG(substr,RGMLE,'ARRAY OF UNRECOGNIZED TYPE CANNOT BE UN-SORTED !')
     END SELECT
     ! CLEAN UP
     CALL INIT_NARRAY(nh)

  ELSE

     CALL INIT_NARRAY(nx, na%n, na%dim, VTYPE_INT)
     DEALLOCATE(nx%vi)
     NULLIFY(nx%vi)
     !
     vtype = na%type
     SELECT CASE(vtype)
     CASE(VTYPE_REAL)
        CALL QSORT_R(na%vr, nx%vi)
     CASE(VTYPE_DOUBLE)
        CALL QSORT_D(na%vd, nx%vi)
     CASE(VTYPE_INT)
        CALL QSORT_I(na%vi, nx%vi)
     CASE(VTYPE_BYTE)
        CALL QSORT_B(na%vb, nx%vi)
     CASE(VTYPE_CHAR)
     CALL RGMSG(substr,RGMLE,'SORTING OF TYPE CHAR IS NOT IMPLEMENTED !')
     CASE(VTYPE_UNDEF)
     CALL RGMSG(substr,RGMLE,'ARRAY OF UNDEFINED TYPE CANNOT BE SORTED !')
     CASE DEFAULT
     CALL RGMSG(substr,RGMLE,'ARRAY OF UNRECOGNIZED TYPE CANNOT BE SORTED !')
     END SELECT

  END IF ! (reverse ?)

END SUBROUTINE SORT_NARRAY
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE REORDER_NARRAY(na, nx)

  IMPLICIT NONE

  ! I/O
  TYPE (narray), INTENT(INOUT) :: na   ! n-array to reorder
  TYPE (narray), INTENT(IN)    :: nx   ! index n-array

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'REORDER_NARRAY'
  TYPE (narray) :: nh   ! copy of na
  INTEGER       :: vtype
  INTEGER       :: i, n
  
  IF (na%n == 0) THEN
     CALL RGMSG(substr, RGMLW, 'EMPTY ARRAY ! NOTHING TO DO !')
     RETURN
  END IF

  vtype = nx%type
  IF (vtype /= VTYPE_INT) THEN
     CALL RGMSG(substr, RGMLE, 'INDEX N-ARRAY MUST BE OF TYPE INT !')
  END IF

  IF (na%n > 1) THEN
     CALL RGMSG(substr, RGMLW, 'REORDERING A ',na%n, &
          '-DIMENSIONAL N-ARRAY AS LINEAR ARRAY !')
  END IF

  IF (na%n /= nx%n) THEN
     CALL RGMSG(substr, RGMLE, 'DIMENSION MISMATCH BETWEEN N-ARRAY (',&
          na%n,')' , .false.)
     CALL RGMSG(substr, RGMLEC, 'AND INDEX N-ARRAY (',nx%n,') !')
  END IF

  IF (PRODUCT(na%dim) /= PRODUCT(nx%dim)) THEN
     CALL RGMSG(substr, RGMLE, 'LENGTH OF N-ARRAY MISMATCH BETWEEN', .false.)
     CALL RGMSG(substr, RGMLEC, 'N-ARRAY (',PRODUCT(na%dim),') AND', .false.)
     CALL RGMSG(substr, RGMLEC, 'INDEX N-ARRAY (',PRODUCT(nx%dim),') !')
  END IF

  vtype = na%type
  CALL COPY_NARRAY(nh, na)
  n = PRODUCT(na%dim)

  SELECT CASE(vtype)
  CASE(VTYPE_REAL)
     DO i=1,n
        na%vr(i) = nh%vr(nx%vi(i))
     END DO
  CASE(VTYPE_DOUBLE)
     DO i=1,n
        na%vd(i) = nh%vd(nx%vi(i))
     END DO
  CASE(VTYPE_INT)
     DO i=1,n
        na%vi(i) = nh%vi(nx%vi(i))
     END DO
  CASE(VTYPE_BYTE)
     DO i=1,n
        na%vb(i) = nh%vb(nx%vi(i))
     END DO
  CASE(VTYPE_CHAR)
     DO i=1,n
        na%vc(i) = nh%vc(nx%vi(i))
     END DO
  CASE(VTYPE_UNDEF)
  CALL RGMSG(substr, RGMLE, 'ARRAY OF UNDEFINED TYPE CANNOT BE UN-SORTED !')
  CASE DEFAULT
  CALL RGMSG(substr, RGMLE, 'ARRAY OF UNRECOGNIZED TYPE CANNOT BE UN-SORTED !')
  END SELECT
  
  ! CLEAN UP
  CALL INIT_NARRAY(nh)

END SUBROUTINE REORDER_NARRAY
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE DOUBLE_NARRAY(na)

  IMPLICIT NONE

  ! I/O
  TYPE (narray), INTENT(INOUT) :: na

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'DOUBLE_NARRAY'
  INTEGER :: vtype
  INTEGER :: status

  vtype = na%type
  SELECT CASE(vtype)
  CASE(VTYPE_REAL)
     ALLOCATE(na%vd(SIZE(na%vr)), STAT=status)
     CALL ERRMSG(substr,status,1)
     na%vd(:) = REAL(na%vr(:), DP)
     na%type = VTYPE_DOUBLE
     DEALLOCATE(na%vr, STAT=status)
     CALL ERRMSG(substr,status,2)
     NULLIFY(na%vr)
  CASE(VTYPE_DOUBLE)
     ! NOTHING TO DO
  CASE(VTYPE_INT)
     ALLOCATE(na%vd(SIZE(na%vi)), STAT=status)
     CALL ERRMSG(substr,status,3)
     na%vd(:) = REAL(na%vi(:), DP)
     na%type = VTYPE_DOUBLE
     DEALLOCATE(na%vi, STAT=status)
     CALL ERRMSG(substr,status,4)
     NULLIFY(na%vi)
  CASE(VTYPE_BYTE)
     ALLOCATE(na%vd(SIZE(na%vb)), STAT=status)
     CALL ERRMSG(substr,status,5)
     na%vd(:) = REAL(na%vb(:), DP)
     na%type = VTYPE_DOUBLE
     DEALLOCATE(na%vb, STAT=status)
     CALL ERRMSG(substr,status,6)
     NULLIFY(na%vb)
  CASE(VTYPE_CHAR)
     CALL RGMSG(substr, RGMLE, 'CHAR CANNOT BE CONVERTED TO DOUBLEPRECISION !')
  CASE(VTYPE_UNDEF)
     CALL RGMSG(substr, RGMLE, 'UNDEFINED N-ARRAY TYPE !')
  CASE DEFAULT
     CALL RGMSG(substr, RGMLE, 'UNRECOGNIZED N-ARRAY TYPE !')
  END SELECT

END SUBROUTINE DOUBLE_NARRAY
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE SCALE_NARRAY(na, sc)

  IMPLICIT NONE

  ! I/O
  TYPE (narray), INTENT(INOUT) :: na  ! N-array
  REAL(dp)     , INTENT(IN)    :: sc  ! scaling factor

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'SCALE_NARRAY'
  INTEGER :: vtype
  INTEGER :: status
  INTEGER :: i

  vtype = na%type

  SELECT CASE(vtype)
  CASE(VTYPE_REAL)
! op_mr_20170322+
!!$  na%vr(:) = na%vr(:) * REAL(sc, SP)
     na%vr(:) = REAL(REAL(na%vr(:), DP) * sc, SP)
! op_mr_20170322-
  CASE(VTYPE_DOUBLE)
!!$     na%vd(:) = na%vd(:) * REAL(sc, DP) 
     na%vd(:) = na%vd(:) * sc 
  CASE(VTYPE_INT)
     CALL RGMSG(substr, RGMLI, 'N-ARRAY OF TYPE INTEGER CONVERTED TO REAL !')
     ALLOCATE(na%vr(SIZE(na%vi)), STAT=status)
     na%type = VTYPE_REAL
     CALL ERRMSG(substr,status,1)
     DO i=1, SIZE(na%vi)
        na%vr(i) = REAL(na%vi(i), SP) * REAL(sc, SP)
     END DO
     DEALLOCATE(na%vi, STAT=status)
     CALL ERRMSG(substr,status,2)
     NULLIFY(na%vi)
  CASE(VTYPE_BYTE)
     CALL RGMSG(substr, RGMLI, 'N-ARRAY OF TYPE BYTE CONVERTED TO REAL !')
     ALLOCATE(na%vr(SIZE(na%vb)), STAT=status)
     na%type = VTYPE_REAL
     CALL ERRMSG(substr,status,3)
     DO i=1, SIZE(na%vb)
        na%vr(i) = REAL(INT(na%vb(i), I8), SP) * REAL(sc, SP)
     END DO
     DEALLOCATE(na%vb, STAT=status)
     CALL ERRMSG(substr,status,4)
     NULLIFY(na%vb)
  CASE(VTYPE_CHAR)
     CALL RGMSG(substr, RGMLE, 'N-ARRAY OF TYPE CHAR CANNOT BE SCALED !')
  CASE(VTYPE_UNDEF)
     CALL RGMSG(substr, RGMLE, 'UNDEFINED N-ARRAY TYPE !')
  CASE DEFAULT
     CALL RGMSG(substr, RGMLE, 'UNRECOGNIZED N-ARRAY TYPE !')
  END SELECT
  
END SUBROUTINE SCALE_NARRAY
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE CAT_NARRAY(na, nb)

  IMPLICIT NONE

  ! I/O
  TYPE(narray), INTENT(INOUT) :: na
  TYPE(narray), INTENT(in)    :: nb

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'CAT_NARRAY'
  TYPE(narray) :: nh
  INTEGER      :: vtype1, vtype2
  INTEGER      :: n,m,i

  vtype1 = na%type
  vtype2 = nb%type

  IF (vtype2 == VTYPE_UNDEF) THEN
     CALL RGMSG(substr, RGMLE, 'N-ARRAY TO BE APPENDED MUST BE DEFINED !')
  ELSE
     IF (vtype1 == VTYPE_UNDEF) THEN
        CALL COPY_NARRAY(na, nb)
        RETURN
     ELSE
        IF (vtype1 /= vtype2) THEN
           CALL RGMSG(substr, RGMLE, 'N-ARRAYS MUST BE OF SAME TYPE !')
        END IF
     END IF
     IF ((na%n /= 1).OR.(nb%n /= 1)) THEN
        CALL RGMSG(substr, RGMLE, 'N-ARRAYS MUST BE 1-DIMENSIONAL !')
     END IF
  END IF

  n = na%dim(1)
  m = nb%dim(1)
  CALL INIT_NARRAY(nh, 1, (/ n+m /), vtype2)

  SELECT CASE(vtype2)
  CASE(VTYPE_REAL)
     DO i=1, n
        nh%vr(i) = na%vr(i) 
     END DO
     DO i=1, m
        nh%vr(n+i) = nb%vr(i) 
     END DO
  CASE(VTYPE_DOUBLE)
     DO i=1, n
        nh%vd(i) = na%vd(i) 
     END DO
     DO i=1, m
        nh%vd(n+i) = nb%vd(i) 
     END DO
  CASE(VTYPE_INT)
     DO i=1, n
        nh%vi(i) = na%vi(i) 
     END DO
     DO i=1, m
        nh%vi(n+i) = nb%vi(i) 
     END DO
  CASE(VTYPE_BYTE)
     DO i=1, n
        nh%vb(i) = na%vb(i) 
     END DO
     DO i=1, m
        nh%vb(n+i) = nb%vb(i) 
     END DO
  CASE(VTYPE_CHAR)
     DO i=1, n
        nh%vc(i) = na%vc(i) 
     END DO
     DO i=1, m
        nh%vc(n+i) = nb%vc(i) 
     END DO
  CASE(VTYPE_UNDEF)
     ! NOTHING TO DO FOR UNDEFINED ARRAYS
  CASE DEFAULT
     ! NOTHING TO DO FOR UNRECOGNIZED ARRAYS
  END SELECT  

  CALL COPY_NARRAY(na, nh)
  ! CLEAN UP
  CALL INIT_NARRAY(nh)

END SUBROUTINE CAT_NARRAY
! ------------------------------------------------------------------

! ------------------------------------------------------------------
RECURSIVE SUBROUTINE QSORT_I(data,idx,ileft,iright)

  IMPLICIT NONE

  ! I/O
  INTEGER (I8), DIMENSION(:), INTENT(INOUT)        :: data   ! data to sort
  INTEGER (I8), DIMENSION(:), POINTER              :: idx    ! index list
  INTEGER (I8),               INTENT(IN), OPTIONAL :: ileft, iright

  ! LOCAL 
  CHARACTER(LEN=*), PARAMETER :: substr = 'QSORT_I'
  INTEGER (I8)  :: temp            ! temporal data
  INTEGER       :: n               ! LENGTH OF LIST
  INTEGER (I8)  :: left, right
  INTEGER (I8)  :: i, last, apu
  INTEGER (I8)  :: ti              ! temporal index
  INTEGER       :: status

  ! INIT
  n = SIZE(data)
  IF (.NOT.ASSOCIATED(idx)) THEN
     ALLOCATE(idx(n), STAT=status)
     CALL ERRMSG(substr,status,1)
     DO i=1, n
        idx(i) = i
     END DO
  END IF

  IF (PRESENT(ileft)) THEN
     left = ileft
  ELSE
     left = 1
  END IF
  IF (PRESENT(iright)) THEN
     right = iright
  ELSE
     right = n
  END IF

  IF(left >= right) RETURN

  apu = (left+right)/2

  temp = data(left)
  data(left) = data(apu)
  data(apu) = temp

  ti = idx(left)
  idx(left) = idx(apu)
  idx(apu) = ti

  last = left

  DO i=left+1,right
     IF(data(i) < data(left)) THEN
        last = last+1
        temp = data(last)
        data(last) = data(i)
        data(i) = temp

        ti = idx(last)
        idx(last) = idx(i)
        idx(i) = ti
     ENDIF
  END DO

  temp = data(left)
  data(left) = data(last)
  data(last) = temp

  ti = idx(left)
  idx(left) = idx(last)
  idx(last) = ti

  CALL QSORT_I(data,idx,left,last-1)
  CALL QSORT_I(data,idx,last+1,right)

  RETURN

END SUBROUTINE QSORT_I
! ------------------------------------------------------------------

! ------------------------------------------------------------------
RECURSIVE SUBROUTINE QSORT_B(data,idx,ileft,iright)

  IMPLICIT NONE

  ! I/O
  INTEGER (I4), DIMENSION(:), INTENT(INOUT)        :: data   ! data to sort
  INTEGER (I8), DIMENSION(:), POINTER              :: idx    ! index list
  INTEGER (I8),               INTENT(IN), OPTIONAL :: ileft, iright

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'QSORT_B'
  INTEGER (I4)  :: temp            ! temporal data
  INTEGER       :: n               ! LENGTH OF LIST
  INTEGER (I8)  :: left, right
  INTEGER (I8)  :: i, last, apu
  INTEGER (I8)  :: ti              ! temporal index
  INTEGER       :: status

  ! INIT
  n = SIZE(data)
  IF (.NOT.ASSOCIATED(idx)) THEN
     ALLOCATE(idx(n), STAT=status)
     CALL ERRMSG(substr,status,1)
     DO i=1, n
        idx(i) = i
     END DO
  END IF

  IF (PRESENT(ileft)) THEN
     left = ileft
  ELSE
     left = 1
  END IF
  IF (PRESENT(iright)) THEN
     right = iright
  ELSE
     right = n
  END IF

  IF(left >= right) RETURN

  apu = (left+right)/2

  temp = data(left)
  data(left) = data(apu)
  data(apu) = temp

  ti = idx(left)
  idx(left) = idx(apu)
  idx(apu) = ti

  last = left

  DO i=left+1,right
     IF(data(i) < data(left)) THEN
        last = last+1
        temp = data(last)
        data(last) = data(i)
        data(i) = temp

        ti = idx(last)
        idx(last) = idx(i)
        idx(i) = ti
     ENDIF
  END DO

  temp = data(left)
  data(left) = data(last)
  data(last) = temp

  ti = idx(left)
  idx(left) = idx(last)
  idx(last) = ti

  CALL QSORT_B(data,idx,left,last-1)
  CALL QSORT_B(data,idx,last+1,right)

  RETURN

END SUBROUTINE QSORT_B
! ------------------------------------------------------------------

! ------------------------------------------------------------------
RECURSIVE SUBROUTINE QSORT_R(data,idx,ileft,iright)

  IMPLICIT NONE

  ! I/O
  REAL (SP),    DIMENSION(:), INTENT(INOUT)        :: data   ! data to sort
  INTEGER (I8), DIMENSION(:), POINTER              :: idx    ! index list
  INTEGER (I8),               INTENT(IN), OPTIONAL :: ileft, iright

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'QSORT_R'
  REAL (SP)    :: temp            ! temporal data
  INTEGER      :: n               ! LENGTH OF LIST
  INTEGER (I8) :: left, right
  INTEGER (I8) :: i, last, apu
  INTEGER (I8) :: ti              ! temporal index
  INTEGER      :: status

  ! INIT
  n = SIZE(data)
  IF (.NOT.ASSOCIATED(idx)) THEN
     ALLOCATE(idx(n),STAT=status)
     CALL ERRMSG(substr,status,1)
     DO i=1, n
        idx(i) = i
     END DO
  END IF

  IF (PRESENT(ileft)) THEN
     left = ileft
  ELSE
     left = 1
  END IF
  IF (PRESENT(iright)) THEN
     right = iright
  ELSE
     right = n
  END IF

  IF(left >= right) RETURN

  apu = (left+right)/2

  temp = data(left)
  data(left) = data(apu)
  data(apu) = temp

  ti = idx(left)
  idx(left) = idx(apu)
  idx(apu) = ti

  last = left

  DO i=left+1,right
     IF(data(i) < data(left)) THEN
        last = last+1
        temp = data(last)
        data(last) = data(i)
        data(i) = temp

        ti = idx(last)
        idx(last) = idx(i)
        idx(i) = ti
     ENDIF
  END DO

  temp = data(left)
  data(left) = data(last)
  data(last) = temp

  ti = idx(left)
  idx(left) = idx(last)
  idx(last) = ti

  CALL QSORT_R(data,idx,left,last-1)
  CALL QSORT_R(data,idx,last+1,right)

  RETURN

END SUBROUTINE QSORT_R
! ------------------------------------------------------------------

! ------------------------------------------------------------------
RECURSIVE SUBROUTINE QSORT_D(data,idx,ileft,iright)

  IMPLICIT NONE

  ! I/O
  REAL (DP),    DIMENSION(:), INTENT(INOUT)        :: data   ! data to sort
  INTEGER (I8), DIMENSION(:), POINTER              :: idx    ! index list
  INTEGER (I8),               INTENT(IN), OPTIONAL :: ileft, iright

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'QSORT_D'
  REAL (DP)    :: temp            ! temporal data
  INTEGER      :: n               ! LENGTH OF LIST
  INTEGER (I8) :: left, right
  INTEGER (I8) :: i, last, apu
  INTEGER (I8) :: ti              ! temporal index
  INTEGER      :: status

  ! INIT
  n = SIZE(data)
  IF (.NOT.ASSOCIATED(idx)) THEN
     ALLOCATE(idx(n), STAT=status)
     CALL ERRMSG(substr,status,1)
     DO i=1, n
        idx(i) = i
     END DO
  END IF

  IF (PRESENT(ileft)) THEN
     left = ileft
  ELSE
     left = 1
  END IF
  IF (PRESENT(iright)) THEN
     right = iright
  ELSE
     right = n
  END IF

  IF(left >= right) RETURN

  apu = (left+right)/2

  temp = data(left)
  data(left) = data(apu)
  data(apu) = temp

  ti = idx(left)
  idx(left) = idx(apu)
  idx(apu) = ti

  last = left

  DO i=left+1,right
     IF(data(i) < data(left)) THEN
        last = last+1
        temp = data(last)
        data(last) = data(i)
        data(i) = temp

        ti = idx(last)
        idx(last) = idx(i)
        idx(i) = ti
     ENDIF
  END DO

  temp = data(left)
  data(left) = data(last)
  data(last) = temp

  ti = idx(left)
  idx(left) = idx(last)
  idx(last) = ti

  CALL QSORT_D(data,idx,left,last-1)
  CALL QSORT_D(data,idx,last+1,right)

  RETURN

END SUBROUTINE QSORT_D
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
!  REAL (DP), DIMENSION(:,:), POINTER  :: mfs, mfd ! overlap matrices
  REAL (DP), DIMENSION(:,:)           :: mfs, mfd ! overlap matrices
  LOGICAL,                 INTENT(IN) :: lmod     ! modulo axis ?

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'OVL_1D_RR'
  INTEGER   :: n,m    ! dimensions
  INTEGER   :: i,j    ! counter
!  INTEGER   :: status
  REAL (DP) :: vmod   ! modulo value

  ! INIT
  n = SIZE(s)-1
  m = SIZE(d)-1

  ! ALLOCATE MEMORY
!  ALLOCATE(mfs(n,m),STAT=status)
!  CALL ERRMSG('OVL_1D_RR',status,1)
!  ALLOCATE(mfd(m,n),STAT=status)
!  CALL ERRMSG('OVL_1D_RR',status,2)

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
!  REAL (DP), DIMENSION(:,:), POINTER  :: mfs, mfd ! overlap matrices
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

  ! ALLOCATE MEMORY
!  ALLOCATE(mfs(n,m),STAT=status)
!  CALL ERRMSG('OVL_1D_RD',status,1)
!  ALLOCATE(mfd(m,n),STAT=status)
!  CALL ERRMSG('OVL_1D_RD',status,2)

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
!  REAL (DP), DIMENSION(:,:), POINTER  :: mfs, mfd ! overlap matrices
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

  ! ALLOCATE MEMORY
!  ALLOCATE(mfs(n,m),STAT=status)
!  CALL ERRMSG('OVL_1D_DR',status,1)
!  ALLOCATE(mfd(m,n),STAT=status)
!  CALL ERRMSG('OVL_1D_DR',status,2)

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
!  REAL (DP), DIMENSION(:,:), POINTER  :: mfs, mfd ! overlap matrices
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

  ! ALLOCATE MEMORY
!  ALLOCATE(mfs(n,m),STAT=status)
!  CALL ERRMSG('OVL_1D_DD',status,1)
!  ALLOCATE(mfd(m,n),STAT=status)
!  CALL ERRMSG('OVL_1D_DD',status,2)

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
INTEGER FUNCTION POSITION(dim, vec)

! THIS FUNCTION CALCULATES THE POSITION NUMBER IN A
! 1-D (LINEAR) ARRAY, GIVEN THAT THE ARRAY SHOULD BE INTERPRETED
! AS N-D ARRAY WITH DIMENSIONS
! dim = (d1, d2, d3, ..., dN)
! OF THE ELEMENT
! vec = (v1, v2, v3, ..., vN)
!

  IMPLICIT NONE

  ! I/O
  INTEGER, DIMENSION(:), INTENT(IN) :: dim
  INTEGER, DIMENSION(:), INTENT(IN) :: vec

  ! LOCAL
  INTEGER :: i
  INTEGER :: n
  INTEGER :: dacc

  IF (SIZE(dim) /= SIZE(vec)) THEN
     POSITION = 0
     RETURN
  END IF

  DO i=1, SIZE(dim)
     IF (vec(i) > dim(i)) THEN
        POSITION = 0
        RETURN
     END IF
  END DO

  n = vec(1)
  dacc = 1
  DO i=2,SIZE(dim)
     dacc = dacc*dim(i-1)
     n = n + dacc*(vec(i)-1)
  END DO

  POSITION = n

END FUNCTION POSITION
! ------------------------------------------------------------------

#ifdef _XOLD
! ------------------------------------------------------------------
SUBROUTINE ELEMENT(dim, n, vec)

! THIS SUBROUTINE CALCULATES THE ELEMENT VECTOR
! vec = (v1, v2, v3, ..., vN)
! OF THE ELEMENT WITH POSITION n IN A
! 1-D (LINEAR) ARRAY, GIVEN THAT THE ARRAY SHOULD BE INTERPRETED
! AS N-D ARRAY WITH DIMENSIONS
! dim = (d1, d2, d3, ..., dN)
!

  IMPLICIT NONE

  ! I/O
  INTEGER, DIMENSION(:), INTENT(IN)  :: dim   ! dimension vector
  INTEGER,               INTENT(IN)  :: n     ! element in linear array
  INTEGER, DIMENSION(:), POINTER     :: vec   ! element vector

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'ELEMENT'
  INTEGER                             :: m    ! COPY of n
  INTEGER                             :: i    ! counter
  INTEGER , DIMENSION(:), ALLOCATABLE :: dacc
  INTEGER                             :: l    ! length of dim
  INTEGER                             :: status

  m = n
  l = SIZE(dim)
  IF (ASSOCIATED(vec)) THEN
     IF (SIZE(vec) > 0) THEN
        DEALLOCATE(vec, STAT=status)
        CALL ERRMSG(substr,status,1)
     END IF
     NULLIFY(vec)
  END IF
  ALLOCATE(vec(l), STAT=status)
  CALL ERRMSG(substr,status,2)
  vec(:) = 0

  ALLOCATE(dacc(l))
  dacc(1) = 1
  DO i=2, l
     dacc(i) = dacc(i-1)*dim(i-1)
  END DO

  IF (m > dacc(l)*dim(l)) RETURN

  DO i=l, 2, -1
     vec(i) = (m-1)/dacc(i)+1
     m = m - (vec(i)-1)*dacc(i)
  END DO
  vec(1) = m

  DEALLOCATE(dacc, stat=STATUS)
  CALL ERRMSG(substr,status,3)

END SUBROUTINE ELEMENT
! ------------------------------------------------------------------
#else
! ------------------------------------------------------------------
SUBROUTINE ELEMENT(dim, n, vec)

! THIS SUBROUTINE CALCULATES THE ELEMENT VECTOR
! vec = (v1, v2, v3, ..., vN)
! OF THE ELEMENT WITH POSITION n IN A
! 1-D (LINEAR) ARRAY, GIVEN THAT THE ARRAY SHOULD BE INTERPRETED
! AS N-D ARRAY WITH DIMENSIONS
! dim = (d1, d2, d3, ..., dN)
!

  IMPLICIT NONE
  INTRINSIC :: MOD

  ! I/O
  INTEGER, DIMENSION(:), INTENT(IN)  :: dim   ! dimension vector
  INTEGER,               INTENT(IN)  :: n     ! element in linear array
  INTEGER, DIMENSION(:), POINTER     :: vec   ! element vector

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'ELEMENT'
  INTEGER                             :: m    ! COPY of n
  INTEGER                             :: i    ! counter
  INTEGER                             :: l    ! length of dim
  INTEGER                             :: status

  m = n
  l = SIZE(dim)
  IF (ASSOCIATED(vec)) THEN
     IF (SIZE(vec) > 0) THEN
        DEALLOCATE(vec, STAT=status)
        CALL ERRMSG(substr,status,1)
     END IF
     NULLIFY(vec)
  END IF
  ALLOCATE(vec(l), STAT=status)
  CALL ERRMSG(substr,status,2)

  m=m-1
  DO i=1,l-1
     vec(i)=1+mod(m, dim(i))
     m=m/dim(i)
  ENDDO
  vec(l)=1+m
  if (vec(l)>dim(l)) then
     vec(:) = 0
     return
  endif

END SUBROUTINE ELEMENT
! ------------------------------------------------------------------
#endif

! ------------------------------------------------------------------
RECURSIVE SUBROUTINE NREGRID(s, sax, dax, d, RG_TYPE       &
                             ,sovl,  dovl, rcnt            &
                             ,gm, gnr, gsvec, gdvec, gdio  &
                             ,gsf, gdf                     &
                            )

  IMPLICIT NONE

  ! I/O
  TYPE (narray), DIMENSION(:), INTENT(IN)    :: s       ! source data fields
  TYPE (narray), DIMENSION(:), POINTER       :: d       ! dest. data fields
  TYPE (axis)  , DIMENSION(:), INTENT(IN)    :: sax     ! axes (source)
  TYPE (axis)  , DIMENSION(:), INTENT(INOUT) :: dax     ! axes (dest)
  INTEGER      , DIMENSION(:), INTENT(IN)    :: RG_TYPE ! regridding type
  REAL (DP),     DIMENSION(:), POINTER       :: sovl  ! glob. overlap fraction
  REAL (DP),     DIMENSION(:), POINTER       :: dovl  ! glob. overlap fraction
  INTEGER,       DIMENSION(:), POINTER       :: rcnt  ! recursion level counter

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
  TYPE (axis)  , DIMENSION(:), ALLOCATABLE   :: psax   ! axes (source)
  TYPE (axis)  , DIMENSION(:), ALLOCATABLE   :: pdax   ! axes (dest)
  TYPE (narray), DIMENSION(:), POINTER       :: pd     ! dest. data
  REAL (DP),     DIMENSION(:), POINTER       :: psovl  ! glob. overlap fraction
  REAL (DP),     DIMENSION(:), POINTER       :: pdovl  ! glob. overlap fraction
  INTEGER,       DIMENSION(:), POINTER       :: prcnt ! recursion level counter

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'NREGRID'
  INTEGER                                 :: i,j, k ! counter
  INTEGER                                 :: nvar   ! number of fields
  INTEGER                                 :: id     ! dimenions loop counter
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
           CALL RGMSG(substr, RGMLI, 'FOUND INVRIANT DEPENDENT DIMENSION')
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
        lflag = ((sax(id)%ndp > 1).OR.      &  !  source dim. independent
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
        CALL RGMSG(substr, RGMLE, 'AMBIGUOUS DIMENSION !')
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
           ! NOTE: IDX VAR WAS CONVERTED TO INDEX-FRACTION
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
  TYPE (axis),   DIMENSION(:), INTENT(IN) :: sax, dax
  REAL (DP),     DIMENSION(:), INTENT(IN) :: sovl, dovl
  TYPE (narray), DIMENSION(:), INTENT(IN) :: nai, nao
  INTEGER,       DIMENSION(:), INTENT(IN) :: rcnt

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

! ------------------------------------------------------------------
SUBROUTINE ERRMSG(routine, status, pos)

  IMPLICIT NONE

  ! I/O
  CHARACTER(LEN=*), INTENT(IN)  :: routine
  INTEGER,          INTENT(IN)  :: status
  INTEGER,          INTENT(IN)  :: pos

  IF (status == 0) THEN
     RETURN
  ELSE
     CALL RGMSG(routine, RGMLE,  'ERROR STATUS ',status,' ', .false.)
     CALL RGMSG(routine, RGMLEC, 'AT POSITION ',pos,' !')
  END IF

END SUBROUTINE ERRMSG
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE RGMSG_C(routine, level, c, lstop)

#if ! defined(NOMPI)
  USE messy_ncregrid_mpi,  ONLY: ncregrid_abort
#endif

  IMPLICIT NONE

  ! I/O
  CHARACTER(LEN=*), INTENT(IN)            :: routine
  INTEGER,          INTENT(IN)            :: level
  CHARACTER(LEN=*), INTENT(IN)            :: c
  LOGICAL,          INTENT(IN), OPTIONAL  :: lstop

  ! LOCAL
  LOGICAL :: llstop
  ! mz_pj_20080801+
  INTEGER           :: iou
  LOGICAL           :: opened
  CHARACTER(LEN=20) :: errfname = ''
  ! mz_pj_20080801-

  ! INIT
  IF (PRESENT(lstop)) THEN
     llstop = lstop
  ELSE
     IF ((level == RGMLE).OR.(level == RGMLEC)) THEN
        llstop = .true.      ! STOP ON ERROR
     ELSE
        llstop = .false.
     END IF
  END IF

  ! mz_pj_20080801+
  IF ((level == RGMLE).OR.(level == RGMLEC)) THEN
     DO iou=100,300
        INQUIRE(unit=iou,opened=opened)
        IF (.NOT.opened) EXIT
     END DO
     IF (MY_RANK < 0) THEN
        errfname = 'ERROR.ncregrid'
     ELSE
        WRITE(errfname,'(a14,a1,i4.4)') 'ERROR.ncregrid','.',MY_RANK
     ENDIF
  END IF
  ! mz_pj_20080801-

  SELECT CASE(level)
  CASE(RGMLE)   ! ERROR MESSAGE
     IF (IAND(MSGMODE, MSGMODE_E) == MSGMODE_E) THEN
        WRITE(*,*) '*** ',TRIM(routine),' ERROR: '
        WRITE(*,*) '    ',TRIM(c)
        ! mz_pj_20080801+
        OPEN(iou, FILE=TRIM(errfname), STATUS='UNKNOWN')
        WRITE(iou,*) '*** ',TRIM(routine),' ERROR: '
        WRITE(iou,*) '    ',TRIM(c)
        CLOSE(iou)
        ! mz_pj_20080801-
     END IF
  CASE(RGMLEC) ! ERROR MESSAGE CONTINUED
     IF (IAND(MSGMODE, MSGMODE_E) == MSGMODE_E) THEN
        WRITE(*,*) '    ',TRIM(c)
        ! mz_pj_20080801+
        OPEN(iou, FILE=TRIM(errfname), STATUS='OLD',POSITION='APPEND')
        WRITE(iou,*) '    ',TRIM(c)
        CLOSE(iou)
        ! mz_pj_20080801-
     END IF
  CASE(RGMLVL)  ! LITTLE VERBOSE
     IF (IAND(MSGMODE, MSGMODE_VL) == MSGMODE_VL) THEN
        WRITE(*,*) TRIM(c)
     END IF
  CASE(RGMLVLC) ! LITTLE VERBOSE CONTINUED
     IF (IAND(MSGMODE, MSGMODE_VL) == MSGMODE_VL) THEN
        WRITE(*,*) '    ',TRIM(c)
     END IF
  CASE(RGMLW) ! WARNING MESSAGE
     IF (IAND(MSGMODE, MSGMODE_W) == MSGMODE_W) THEN
        WRITE(*,*) '+++ ',TRIM(routine),' WARNING: '
        WRITE(*,*) '    ',TRIM(c)
     END IF
  CASE(RGMLWC) ! WARNING MESSAGE CONTINUED
     IF (IAND(MSGMODE, MSGMODE_W) == MSGMODE_W) THEN
        WRITE(*,*) '    ',TRIM(c)
     END IF
  CASE(RGMLVM)  ! MEDIUM VERBOSE
     IF (IAND(MSGMODE, MSGMODE_VM) == MSGMODE_VM) THEN
        WRITE(*,*) TRIM(c)
     END IF
  CASE(RGMLVMC) ! MEDIUM VERBOSE CONTINUED
     IF (IAND(MSGMODE, MSGMODE_VM) == MSGMODE_VM) THEN
        WRITE(*,*) '    ',TRIM(c)
     END IF
  CASE(RGMLI)  ! INFO MESSAGE
     IF (IAND(MSGMODE, MSGMODE_I) == MSGMODE_I) THEN
        WRITE(*,*) '=== ',TRIM(routine),' INFO: '
        WRITE(*,*) '    ',TRIM(c)
     END IF
  CASE(RGMLIC) ! INFO MESSAGE CONTINUED
     IF (IAND(MSGMODE, MSGMODE_I) == MSGMODE_I) THEN
        WRITE(*,*) '    ',TRIM(c)
     END IF
  CASE DEFAULT
  END SELECT

#if ! defined(NOMPI)
  IF (llstop) CALL ncregrid_abort('ncregrid')
#else
  IF (llstop) STOP
#endif

END SUBROUTINE RGMSG_C
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE RGMSG_I(routine, level, c1, i, c2, lstop)

  IMPLICIT NONE

  ! I/O
  CHARACTER(LEN=*), INTENT(IN)            :: routine
  INTEGER,          INTENT(IN)            :: level
  CHARACTER(LEN=*), INTENT(IN)            :: c1
  INTEGER,          INTENT(IN)            :: i
  CHARACTER(LEN=*), INTENT(IN)            :: c2
  LOGICAL,          INTENT(IN), OPTIONAL  :: lstop

  ! LOCAL
  CHARACTER(LEN=1000) :: istr = ''

  WRITE(istr,*) TRIM(c1),i,TRIM(c2)
  CALL RGMSG_C(routine, level, TRIM(istr), lstop)

END SUBROUTINE RGMSG_I
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE RGMSG_IA(routine, level, c1, i, c2, lstop)

  IMPLICIT NONE

  ! I/O
  CHARACTER(LEN=*),      INTENT(IN)            :: routine
  INTEGER,               INTENT(IN)            :: level
  CHARACTER(LEN=*),      INTENT(IN)            :: c1
  INTEGER, DIMENSION(:), INTENT(IN)            :: i
  CHARACTER(LEN=*),      INTENT(IN)            :: c2
  LOGICAL,               INTENT(IN), OPTIONAL  :: lstop

  ! LOCAL
  CHARACTER(LEN=1000) :: istr = ''

  WRITE(istr,*) TRIM(c1),i,TRIM(c2)
  CALL RGMSG_C(routine, level, TRIM(istr), lstop)

END SUBROUTINE RGMSG_IA
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE RGMSG_R(routine, level, c1, r, c2, lstop)

  IMPLICIT NONE

  ! I/O
  CHARACTER(LEN=*), INTENT(IN)            :: routine
  INTEGER,          INTENT(IN)            :: level
  CHARACTER(LEN=*), INTENT(IN)            :: c1
  REAL,             INTENT(IN)            :: r
  CHARACTER(LEN=*), INTENT(IN)            :: c2
  LOGICAL,          INTENT(IN), OPTIONAL  :: lstop

  ! LOCAL
  CHARACTER(LEN=1000) :: rstr = ''

  WRITE(rstr,*) TRIM(c1),r,TRIM(c2)
  CALL RGMSG_C(routine, level, TRIM(rstr), lstop)

END SUBROUTINE RGMSG_R
! ------------------------------------------------------------------

! op_mr_20170321+
SUBROUTINE RGMSG_D(routine, level, c1, d, c2, lstop)

  IMPLICIT NONE

  ! I/O
  CHARACTER(LEN=*), INTENT(IN)            :: routine
  INTEGER,          INTENT(IN)            :: level
  CHARACTER(LEN=*), INTENT(IN)            :: c1
  REAL(dp),         INTENT(IN)            :: d
  CHARACTER(LEN=*), INTENT(IN)            :: c2
  LOGICAL,          INTENT(IN), OPTIONAL  :: lstop

  ! LOCAL
  CHARACTER(LEN=1000) :: rstr = ''

  WRITE(rstr,*) TRIM(c1),d,TRIM(c2)
  CALL RGMSG_C(routine, level, TRIM(rstr), lstop)

END SUBROUTINE RGMSG_D
! ------------------------------------------------------------------
! op_mr_20170321-

! ******************************************************************
END MODULE MESSY_NCREGRID_BASE
! ******************************************************************

#else

! ******************************************************************
! ------------------------------------------------------------------
MODULE MESSY_NCREGRID_BASE
! ------------------------------------------------------------------
! Author: Patrick Joeckel, MPICH, Mainz, June 2002
! ******************************************************************

#if defined(MESSY)
  USE messy_main_constants_mem,  ONLY: SP, DP, I4, I8
#else
  USE typeSizes,     ONLY:   SP => FourByteReal  &
                           , DP => EightByteReal &
                           , I4 => FourByteInt   &
                           , I8 => EightByteInt
#endif

  IMPLICIT NONE

  INTRINSIC :: ASSOCIATED, PRESENT, PRODUCT, REAL, SIZE &
             , INT, ABS, DBLE, MAX, MIN, SIGN, TRIM     &
             , MAXVAL, MINVAL, IAND, NULL

  PRIVATE   :: ASSOCIATED, PRESENT, PRODUCT, REAL, SIZE &
             , INT, ABS, DBLE, MAX, MIN, SIGN, TRIM     &
             , MAXVAL, MINVAL, IAND, NULL

  ! REGRIDDING VARIABLE TYPES
  INTEGER, PARAMETER :: RG_INT = 1
  INTEGER, PARAMETER :: RG_EXT = 2
  INTEGER, PARAMETER :: RG_IDX = 3
  INTEGER, PARAMETER :: RG_IXF = 4

  ! REGRIDDING DATA TYPES
  INTEGER, PARAMETER :: VTYPE_UNDEF  = 0
  INTEGER, PARAMETER :: VTYPE_INT    = 1
  INTEGER, PARAMETER :: VTYPE_REAL   = 2
  INTEGER, PARAMETER :: VTYPE_DOUBLE = 3
  INTEGER, PARAMETER :: VTYPE_BYTE   = 4
  INTEGER, PARAMETER :: VTYPE_CHAR   = 5

  ! MESSAGE LEVEL TYPES
  INTEGER, PARAMETER :: RGMLE   = 0   ! ERROR
  INTEGER, PARAMETER :: RGMLEC  = 1   ! ERROR CONTINUED
  INTEGER, PARAMETER :: RGMLVL  = 2   ! LITTLE VERBOSE
  INTEGER, PARAMETER :: RGMLVLC = 3   ! LITTLE VERBOSE CONTINUED
  INTEGER, PARAMETER :: RGMLW   = 4   ! WARNING
  INTEGER, PARAMETER :: RGMLWC  = 5   ! WARNING CONTINUED
  INTEGER, PARAMETER :: RGMLVM  = 6   ! MEDIUM VERBOSE
  INTEGER, PARAMETER :: RGMLVMC = 7   ! MEDIUM VERBOSE CONTINUED
  INTEGER, PARAMETER :: RGMLI   = 8   ! INFO
  INTEGER, PARAMETER :: RGMLIC  = 9   ! INFO CONTINUED

  ! MESSAGE OUTPUT LEVEL
  INTEGER, PARAMETER :: MSGMODE_S  =  0  ! SILENT
  INTEGER, PARAMETER :: MSGMODE_E  =  1  ! ERROR MESSAGES
  INTEGER, PARAMETER :: MSGMODE_VL =  2  ! LITTLE VERBOSE
  INTEGER, PARAMETER :: MSGMODE_W  =  4  ! WARNING MESSAGES
  INTEGER, PARAMETER :: MSGMODE_VM =  8  ! MEDIUM VERBOSE
  INTEGER, PARAMETER :: MSGMODE_I  = 16  ! INFO MESSAGES
  INTEGER, SAVE      :: MSGMODE = MSGMODE_S + MSGMODE_E + MSGMODE_VL &
                                + MSGMODE_W + MSGMODE_VM + MSGMODE_I

  ! mz_pj_20080801+
  INTEGER, SAVE :: MY_RANK = -1
  ! mz_pj_20080801-

  TYPE narray
     ! n-dimenional array as 1D (LINEAR) array (REAL)
     INTEGER                              :: type = VTYPE_UNDEF
     INTEGER (I4)                         :: n = 0  ! number of dimensions
     INTEGER (I4), DIMENSION(:), POINTER  :: dim => NULL()! dim. vector
     REAL (SP)   , DIMENSION(:), POINTER  :: vr => NULL() ! real values
     REAL (DP)   , DIMENSION(:), POINTER  :: vd => NULL() ! double values
     INTEGER (I8), DIMENSION(:), POINTER  :: vi => NULL() ! integer values
     INTEGER (I4), DIMENSION(:), POINTER  :: vb => NULL() ! byte values
     CHARACTER,    DIMENSION(:), POINTER  :: vc => NULL() ! char. values
  END TYPE narray

  TYPE axis
     ! hyper-axis (for curvilinear coordinates)
     TYPE (narray)                        :: dat  ! interface bounds
     LOGICAL                              :: lm = .false. ! modulo axis ?
     INTEGER (I4)                         :: ndp = 0  ! number of dependencies
     ! FIRST IN LIST MUST BE INDEPENDENT !!!
     ! E.G. IF DIM 3 DEPENDS ON DIM 1,2, and 5
     ! dep = (/ 3,1,2,5 /)
     INTEGER (I4), DIMENSION(:), POINTER  :: dep => NULL()  ! dependencies
  END TYPE axis

!! NOTE: DOES NOT WORK PROPERLY FOR SOME COMPILERS ...
!  INTERFACE ASSIGNMENT (=)
!     MODULE PROCEDURE COPY_NARRAY
!     MODULE PROCEDURE COPY_AXIS
!  END INTERFACE

INTERFACE QSORT
   ! THE GOOD OLD (RECURSIVE) QUICKSORT ALGORITHM FOR LINEAR ARRAYS
   MODULE PROCEDURE QSORT_I    ! INTEGER
   MODULE PROCEDURE QSORT_D    ! DOUBLE PRECISION
END INTERFACE

INTERFACE OVL
   ! CALCULATES THE OVERLAP BETWEEN TWO INTERVALS
   !  (AS FRACTION OF THE TWO INTERVALS)
   !  > ALSO APPLICABLE TO 'MODULO' INTERVALS
   MODULE PROCEDURE OVL_DD      ! DOUBLE PRECISION - DOUBLE PRECISION
END INTERFACE

INTERFACE OVL_1D
   ! CALCULATES THE OVERLAP MATRICES BETWEEN TWO CONTINUOUS
   ! INTERVAL SEQUENCES (IN UNITS OF INTERVAL FRACTIONS)
   !  > ALSO APPLICABLE TO 'MODULO' SEQUENCES
   MODULE PROCEDURE OVL_1D_DD      ! DOUBLE PRECISION - DOUBLE PRECISION
END INTERFACE

INTERFACE RGMSG                    ! MESSAGE OUTPUT
   MODULE PROCEDURE RGMSG_C
   MODULE PROCEDURE RGMSG_I
   MODULE PROCEDURE RGMSG_IA
   MODULE PROCEDURE RGMSG_R
   MODULE PROCEDURE RGMSG_D ! op_mr_20170321
END INTERFACE

CONTAINS

! ------------------------------------------------------------------
SUBROUTINE INIT_NARRAY(na, n, dim, qtype)

  IMPLICIT NONE

  ! I/O
  TYPE (narray),         INTENT(INOUT)          :: na
  INTEGER,               INTENT(IN),   OPTIONAL :: n
  INTEGER, DIMENSION(:), INTENT(IN),   OPTIONAL :: dim
  INTEGER,                             OPTIONAL :: qtype

  ! LOCAL
!  CHARACTER(LEN=*), PARAMETER :: substr = 'INIT_NARRAY'
!  INTEGER                     :: status
  INTEGER                     :: len

  na%type = VTYPE_UNDEF

  na%n = 0
  IF (ASSOCIATED(na%dim)) THEN
     IF (SIZE(na%dim) > 0) THEN
        DEALLOCATE(na%dim)
!        DEALLOCATE(na%dim, STAT=status)
!        CALL ERRMSG(substr,status,1)
     END IF
     NULLIFY(na%dim)
  END IF

  IF (ASSOCIATED(na%vi)) THEN
     IF (SIZE(na%vi) > 0) THEN
        DEALLOCATE(na%vi)
!        DEALLOCATE(na%vi, STAT=status)
!        CALL ERRMSG(substr,status,2)
     END IF
     NULLIFY(na%vi)
  END IF

  IF (ASSOCIATED(na%vr)) THEN
     IF (SIZE(na%vr) > 0) THEN
        DEALLOCATE(na%vr)
!        DEALLOCATE(na%vr, STAT=status)
!        CALL ERRMSG(substr,status,3)
     END IF
     NULLIFY(na%vr)
  END IF

  IF (ASSOCIATED(na%vd)) THEN
     IF (SIZE(na%vd) > 0) THEN
        DEALLOCATE(na%vd)
!        DEALLOCATE(na%vd, STAT=status)
!        CALL ERRMSG(substr,status,4)
     END IF
     NULLIFY(na%vd)
  END IF

  IF (ASSOCIATED(na%vc)) THEN
     IF (SIZE(na%vc) > 0) THEN
        DEALLOCATE(na%vc)
!        DEALLOCATE(na%vc, STAT=status)
!        CALL ERRMSG(substr,status,5)
     END IF
     NULLIFY(na%vc)
  END IF

  IF (ASSOCIATED(na%vb)) THEN
     IF (SIZE(na%vb) > 0) THEN
        DEALLOCATE(na%vb)
!        DEALLOCATE(na%vb, STAT=status)
!        CALL ERRMSG(substr,status,6)
     END IF
     NULLIFY(na%vb)
  END IF  

!  IF (PRESENT(dim) .AND. PRESENT(n)) THEN
!     IF (SIZE(dim) /= n) THEN
!        CALL RGMSG(substr, RGMLE,                &
!             'NUMBER OF DIMENSIONS MISMATCH !' )
!     END IF
!  END IF

  len = 0
  IF (PRESENT(n)) THEN
     na%n = n
     ALLOCATE(na%dim(n))
!     ALLOCATE(na%dim(n), STAT=status)
!     CALL ERRMSG(substr,status,7)
     na%dim(:) = 0
  END IF

  IF (PRESENT(dim)) THEN
     IF (.NOT. ASSOCIATED(na%dim)) THEN
        na%n = SIZE(dim)
        ALLOCATE(na%dim(na%n))
!        ALLOCATE(na%dim(na%n), STAT=status)
!        CALL ERRMSG(substr,status,8)
     END IF
     na%dim(:) = dim(:)
     len = PRODUCT(na%dim(:))
  END IF

  IF (PRESENT(qtype).AND.(len > 0)) THEN
     SELECT CASE(qtype)
     CASE(VTYPE_INT)
        na%type = VTYPE_INT
        ALLOCATE(na%vi(len))
!        ALLOCATE(na%vi(len), STAT=status)
!        CALL ERRMSG(substr,status,9)
        na%vi(:) = 0
     CASE(VTYPE_REAL)
        na%type = VTYPE_REAL
        ALLOCATE(na%vr(len))
!        ALLOCATE(na%vr(len), STAT=status)
!        CALL ERRMSG(substr,status,10)
        na%vr(:) = 0.0
     CASE(VTYPE_DOUBLE)
        na%type = VTYPE_DOUBLE
        ALLOCATE(na%vd(len))
!        ALLOCATE(na%vd(len), STAT=status)
!        CALL ERRMSG(substr,status,11)
        na%vd(:) = REAL(0.0, DP)
     CASE(VTYPE_CHAR)
        na%type = VTYPE_CHAR
        ALLOCATE(na%vc(len))
!        ALLOCATE(na%vc(len), STAT=status)
!        CALL ERRMSG(substr,status,12)
        na%vc(:) = ' '
     CASE(VTYPE_BYTE)
        na%type = VTYPE_BYTE
        ALLOCATE(na%vb(len))
!        ALLOCATE(na%vb(len), STAT=status)
!        CALL ERRMSG(substr,status,13)
        na%vi(:) = 0
     CASE(VTYPE_UNDEF)
        na%type = VTYPE_UNDEF
        ! DO NOTHING, KEEP UNDEFINED
     CASE DEFAULT
!        CALL RGMSG(substr, RGMLE,                           &
!             'REQUESTED TYPE FOR N-ARRAY IS UNRECOGNIZED !' )
     END SELECT
  END IF

END SUBROUTINE INIT_NARRAY
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE INIT_AXIS(a)

  IMPLICIT NONE

  ! I/O
  TYPE (axis), INTENT(INOUT) :: a
  
  ! LOCAL
!  CHARACTER(LEN=*), PARAMETER :: substr = 'INIT_AXIS'
!  INTEGER :: status


  CALL INIT_NARRAY(a%dat)
  a%lm = .false.
  a%ndp = 0
  IF (ASSOCIATED(a%dep)) THEN
     IF (SIZE(a%dep) > 0) THEN
        DEALLOCATE(a%dep)
!        DEALLOCATE(a%dep, STAT=status)
!        CALL ERRMSG(substr,status,1)
     END IF
     NULLIFY(a%dep)
  END IF

END SUBROUTINE INIT_AXIS
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE COPY_NARRAY(d, s)

  IMPLICIT NONE

  ! I/O
  TYPE (narray), INTENT(OUT) :: d
  TYPE (narray), INTENT(IN)  :: s
  
  ! LOCAL
!  CHARACTER(LEN=*), PARAMETER :: substr = 'COPY_NARRAY'
  INTEGER :: n
  INTEGER :: vtype
!  INTEGER :: status

  ! INIT
  n = 0

  d%n = s%n
  IF (ASSOCIATED(s%dim)) THEN
     ALLOCATE(d%dim(d%n))
!     ALLOCATE(d%dim(d%n),STAT=status)
!     CALL ERRMSG(substr,status,1)
     d%dim(:) = s%dim(:)
  END IF

!  vtype = QTYPE_NARRAY(s)
  vtype = s%type

  SELECT CASE(vtype)
  CASE(VTYPE_INT)
     d%type = VTYPE_INT
     n = SIZE(s%vi)
     IF (n > 0) THEN
        ALLOCATE(d%vi(n))
!        ALLOCATE(d%vi(n),STAT=status)
!        CALL ERRMSG(substr,status,2)
        d%vi(:) = s%vi(:) 
     END IF
  CASE(VTYPE_REAL)
     d%type = VTYPE_REAL
     n = SIZE(s%vr)
     IF (n > 0) THEN
        ALLOCATE(d%vr(n))
!        ALLOCATE(d%vr(n),STAT=status)
!        CALL ERRMSG(substr,status,3)
        d%vr(:) = s%vr(:)
     END IF
  CASE(VTYPE_DOUBLE)
     d%type = VTYPE_DOUBLE
     n = SIZE(s%vd)
     IF (n > 0) THEN
        ALLOCATE(d%vd(n))
!        ALLOCATE(d%vd(n),STAT=status)
!        CALL ERRMSG(substr,status,4)
        d%vd(:) = s%vd(:)
     END IF
  CASE(VTYPE_CHAR)
     d%type = VTYPE_CHAR
     n = SIZE(s%vc)
     IF (n > 0) THEN
        ALLOCATE(d%vc(n))
!        ALLOCATE(d%vc(n),STAT=status)
!        CALL ERRMSG(substr,status,5)
        d%vc(:) = s%vc(:) 
     END IF
  CASE(VTYPE_BYTE)
     d%type = VTYPE_BYTE
     n = SIZE(s%vb)
     IF (n > 0) THEN
        ALLOCATE(d%vb(n))
!        ALLOCATE(d%vb(n),STAT=status)
!        CALL ERRMSG(substr,status,6)
        d%vb(:) = s%vb(:)
     END IF
  CASE(VTYPE_UNDEF)
     d%type = VTYPE_UNDEF
     ! DO NOTHING, EMPTY N-ARRAY IS COPIED
  CASE DEFAULT
!     CALL RGMSG(substr, RGMLE,                              &
!          'N-ARRAY OF UNRECOGNIZED TYPE CANNOT BE COPIED !' )          
  END SELECT

END SUBROUTINE COPY_NARRAY
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE COPY_AXIS(d, s)

  IMPLICIT NONE

  ! I/O
  TYPE (axis), INTENT(OUT) :: d
  TYPE (axis), INTENT(IN)  :: s

  ! LOCAL
!  CHARACTER(LEN=*), PARAMETER :: substr = 'COPY_AXIS'
!  INTEGER :: status

  CALL COPY_NARRAY(d%dat, s%dat)
  d%lm   = s%lm
  d%ndp  = s%ndp
  IF (ASSOCIATED(s%dep)) THEN
     ALLOCATE(d%dep(d%ndp))
!     ALLOCATE(d%dep(d%ndp),STAT=status)
!     CALL ERRMSG(substr,status,1)
     d%dep(:) = s%dep(:)
  END IF

END SUBROUTINE COPY_AXIS
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE SORT_NARRAY(na, nx, reverse)

  IMPLICIT NONE

  ! I/O
  TYPE (narray), INTENT(INOUT)          :: na
  TYPE (narray), INTENT(INOUT)          :: nx
  LOGICAL      , INTENT(IN)   ,OPTIONAL :: reverse

  ! LOCAL
!  CHARACTER(LEN=*), PARAMETER :: substr = 'SORT_NARRAY'
  INTEGER       :: vtype
  LOGICAL       :: lrev
  INTEGER       :: n, i
  TYPE (narray) :: nh

  IF (PRESENT(reverse)) THEN
     lrev = reverse
  ELSE
     lrev = .false.  ! DEFAULT
  END IF

  IF (na%n == 0) THEN
!     CALL RGMSG(substr, RGMLW, 'EMPTY ARRAY ! NOTHING TO DO !')
     RETURN
  END IF

!  IF (na%n > 1) THEN
!     CALL RGMSG(substr, RGMLW, &
!          'SORTING A ',na%n,'-DIMENSIONAL N-ARRAY AS LINEAR ARRAY !')
!  END IF

  IF (lrev) THEN

     vtype = nx%type
!     IF (vtype /= VTYPE_INT) THEN
!        CALL RGMSG(substr, RGMLE, 'INDEX N-ARRAY MUST BE OF TYPE INTEGER !')
!     END IF
     n = SIZE(nx%vi)
     CALL COPY_NARRAY(nh, na)   ! COPY TO BE RE-ORDERED
     vtype = na%type
     SELECT CASE(vtype)
     CASE(VTYPE_REAL)
        DO i=1,n
           na%vr(nx%vi(i)) = nh%vr(i)
        END DO
     CASE(VTYPE_DOUBLE)
        DO i=1,n
           na%vd(nx%vi(i)) = nh%vd(i)
        END DO
     CASE(VTYPE_INT)
        DO i=1,n
           na%vi(nx%vi(i)) = nh%vi(i)
        END DO        
     CASE(VTYPE_BYTE)
        DO i=1,n
           na%vb(nx%vi(i)) = nh%vb(i)
        END DO        
!     CASE(VTYPE_CHAR)
!     CALL RGMSG(substr,RGMLE,'UN-SORTING OF TYPE CHAR IS NOT IMPLEMENTED !')
!     CASE(VTYPE_UNDEF)
!     CALL RGMSG(substr,RGMLE,'ARRAY OF UNDEFINED TYPE CANNOT BE UN-SORTED !')
     CASE DEFAULT
!   CALL RGMSG(substr,RGMLE,'ARRAY OF UNRECOGNIZED TYPE CANNOT BE UN-SORTED !')
     END SELECT
     ! CLEAN UP
     CALL INIT_NARRAY(nh)

  ELSE

     CALL INIT_NARRAY(nx, na%n, na%dim, VTYPE_INT)
     DEALLOCATE(nx%vi)
     NULLIFY(nx%vi)
     !
     vtype = na%type
     SELECT CASE(vtype)
     CASE(VTYPE_REAL)
        CALL QSORT_R(na%vr, nx%vi)
     CASE(VTYPE_DOUBLE)
        CALL QSORT_D(na%vd, nx%vi)
     CASE(VTYPE_INT)
        CALL QSORT_I(na%vi, nx%vi)
     CASE(VTYPE_BYTE)
        CALL QSORT_B(na%vb, nx%vi)
!     CASE(VTYPE_CHAR)
!     CALL RGMSG(substr,RGMLE,'SORTING OF TYPE CHAR IS NOT IMPLEMENTED !')
!     CASE(VTYPE_UNDEF)
!     CALL RGMSG(substr,RGMLE,'ARRAY OF UNDEFINED TYPE CANNOT BE SORTED !')
     CASE DEFAULT
!     CALL RGMSG(substr,RGMLE,'ARRAY OF UNRECOGNIZED TYPE CANNOT BE SORTED !')
     END SELECT

  END IF ! (reverse ?)

END SUBROUTINE SORT_NARRAY
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE REORDER_NARRAY(na, nx)

  IMPLICIT NONE

  ! I/O
  TYPE (narray), INTENT(INOUT) :: na   ! n-array to reorder
  TYPE (narray), INTENT(IN)    :: nx   ! index n-array

  ! LOCAL
!  CHARACTER(LEN=*), PARAMETER :: substr = 'REORDER_NARRAY'
  TYPE (narray) :: nh   ! copy of na
  INTEGER       :: vtype
  INTEGER       :: i, n
  
  IF (na%n == 0) THEN
!     CALL RGMSG(substr, RGMLW, 'EMPTY ARRAY ! NOTHING TO DO !')
     RETURN
  END IF

  vtype = nx%type
!  IF (vtype /= VTYPE_INT) THEN
!     CALL RGMSG(substr, RGMLE, 'INDEX N-ARRAY MUST BE OF TYPE INT !')
!  END IF

!  IF (na%n > 1) THEN
!     CALL RGMSG(substr, RGMLW, 'REORDERING A ',na%n, &
!          '-DIMENSIONAL N-ARRAY AS LINEAR ARRAY !')
!  END IF

!  IF (na%n /= nx%n) THEN
!     CALL RGMSG(substr, RGMLE, 'DIMENSION MISMATCH BETWEEN N-ARRAY (',&
!          na%n,')' , .false.)
!     CALL RGMSG(substr, RGMLEC, 'AND INDEX N-ARRAY (',nx%n,') !')
!  END IF

!  IF (PRODUCT(na%dim) /= PRODUCT(nx%dim)) THEN
!     CALL RGMSG(substr, RGMLE, 'LENGTH OF N-ARRAY MISMATCH BETWEEN', .false.)
!     CALL RGMSG(substr, RGMLEC, 'N-ARRAY (',PRODUCT(na%dim),') AND', .false.)
!     CALL RGMSG(substr, RGMLEC, 'INDEX N-ARRAY (',PRODUCT(nx%dim),') !')
!  END IF

  vtype = na%type
  CALL COPY_NARRAY(nh, na)
  n = PRODUCT(na%dim)

  SELECT CASE(vtype)
  CASE(VTYPE_REAL)
     DO i=1,n
        na%vr(i) = nh%vr(nx%vi(i))
     END DO
  CASE(VTYPE_DOUBLE)
     DO i=1,n
        na%vd(i) = nh%vd(nx%vi(i))
     END DO
  CASE(VTYPE_INT)
     DO i=1,n
        na%vi(i) = nh%vi(nx%vi(i))
     END DO
  CASE(VTYPE_BYTE)
     DO i=1,n
        na%vb(i) = nh%vb(nx%vi(i))
     END DO
  CASE(VTYPE_CHAR)
     DO i=1,n
        na%vc(i) = nh%vc(nx%vi(i))
     END DO
!  CASE(VTYPE_UNDEF)
!  CALL RGMSG(substr, RGMLE, 'ARRAY OF UNDEFINED TYPE CANNOT BE UN-SORTED !')
  CASE DEFAULT
!  CALL RGMSG(substr, RGMLE, 'ARRAY OF UNRECOGNIZED TYPE CANNOT BE UN-SORTED !')
  END SELECT
  
  ! CLEAN UP
  CALL INIT_NARRAY(nh)

END SUBROUTINE REORDER_NARRAY
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE DOUBLE_NARRAY(na)

  IMPLICIT NONE

  ! I/O
  TYPE (narray), INTENT(INOUT) :: na

  ! LOCAL
!  CHARACTER(LEN=*), PARAMETER :: substr = 'DOUBLE_NARRAY'
  INTEGER :: vtype
!  INTEGER :: status

  vtype = na%type
  SELECT CASE(vtype)
  CASE(VTYPE_REAL)
     ALLOCATE(na%vd(SIZE(na%vr)))
!     ALLOCATE(na%vd(SIZE(na%vr)), STAT=status)
!     CALL ERRMSG(substr,status,1)
     na%vd(:) = REAL(na%vr(:), DP)
     na%type = VTYPE_DOUBLE
     DEALLOCATE(na%vr)
!     DEALLOCATE(na%vr, STAT=status)
!     CALL ERRMSG(substr,status,2)
     NULLIFY(na%vr)
  CASE(VTYPE_DOUBLE)
     ! NOTHING TO DO
  CASE(VTYPE_INT)
     ALLOCATE(na%vd(SIZE(na%vi)))
!     ALLOCATE(na%vd(SIZE(na%vi)), STAT=status)
!     CALL ERRMSG(substr,status,3)
     na%vd(:) = REAL(na%vi(:), DP)
     na%type = VTYPE_DOUBLE
     DEALLOCATE(na%vi)
!     DEALLOCATE(na%vi, STAT=status)
!     CALL ERRMSG(substr,status,4)
     NULLIFY(na%vi)
  CASE(VTYPE_BYTE)
     ALLOCATE(na%vd(SIZE(na%vi)))
!     ALLOCATE(na%vd(SIZE(na%vi)), STAT=status)
!     CALL ERRMSG(substr,status,5)
     na%vd(:) = REAL(na%vb(:), DP)
     na%type = VTYPE_DOUBLE
     DEALLOCATE(na%vb)
!     DEALLOCATE(na%vb, STAT=status)
!     CALL ERRMSG(substr,status,6)
     NULLIFY(na%vb)
!  CASE(VTYPE_CHAR)
!     CALL RGMSG(substr, RGMLE, 'CHAR CANNOT BE CONVERTED TO DOUBLEPRECISION !')
!  CASE(VTYPE_UNDEF)
!     CALL RGMSG(substr, RGMLE, 'UNDEFINED N-ARRAY TYPE !')
  CASE DEFAULT
!     CALL RGMSG(substr, RGMLE, 'UNRECOGNIZED N-ARRAY TYPE !')
  END SELECT

END SUBROUTINE DOUBLE_NARRAY
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE SCALE_NARRAY(na, sc)

  IMPLICIT NONE

  ! I/O
  TYPE (narray), INTENT(INOUT) :: na  ! N-array
  REAL(dp)     , INTENT(IN)    :: sc  ! scaling factor

  ! LOCAL
!  CHARACTER(LEN=*), PARAMETER :: substr = 'SCALE_NARRAY'
  INTEGER :: vtype
!  INTEGER :: status
  INTEGER :: i

  vtype = na%type

  SELECT CASE(vtype)
  CASE(VTYPE_REAL)
! op_mr_20170322+
!!$  na%vr(:) = na%vr(:) * REAL(sc, SP)
     na%vr(:) = REAL(REAL(na%vr(:), DP) * sc, SP)
! op_mr_20170322-
  CASE(VTYPE_DOUBLE)
!!$  na%vd(:) = na%vd(:) * REAL(sc, DP) 
     na%vd(:) = na%vd(:) * sc
  CASE(VTYPE_INT)

!     CALL RGMSG(SUBSTR, RGMLI, 'N-ARRAY OF TYPE INTEGER CONVERTED TO REAL !')
     ALLOCATE(NA%VR(SIZE(NA%VI)))
!     ALLOCATE(NA%VR(SIZE(NA%VI)), STAT=STATUS)
!     CALL ERRMSG(SUBSTR,STATUS,1)
     na%type = VTYPE_REAL
     DO I=1, SIZE(NA%VI)
        NA%VR(I) = REAL(NA%VI(I), SP) * REAL(SC, SP)
     END DO
     DEALLOCATE(NA%VI)
!     DEALLOCATE(NA%VI, STAT=STATUS)
!     CALL ERRMSG(SUBSTR,STATUS,2)
     NULLIFY(NA%VI)
  CASE(VTYPE_BYTE)
!     CALL RGMSG(SUBSTR, RGMLI, 'N-ARRAY OF TYPE BYTE CONVERTED TO REAL !')
     ALLOCATE(NA%VR(SIZE(NA%VB)))
!     ALLOCATE(NA%VR(SIZE(NA%VB)), STAT=STATUS)
!     CALL ERRMSG(SUBSTR,STATUS,3)
     na%type = VTYPE_REAL
     DO I=1, SIZE(NA%VB)
        NA%VR(I) = REAL(INT(NA%VB(I), I8), SP) * REAL(SC, SP)
     END DO
     DEALLOCATE(NA%VB)
!     DEALLOCATE(NA%VB, STAT=STATUS)
!     CALL ERRMSG(SUBSTR,STATUS,4)
     NULLIFY(NA%VB)
!  CASE(VTYPE_CHAR)
!     CALL RGMSG(SUBSTR, RGMLE, 'N-ARRAY OF TYPE CHAR CANNOT BE SCALED !')
!  CASE(VTYPE_UNDEF)
!     CALL RGMSG(SUBSTR, RGMLE, 'UNDEFINED N-ARRAY TYPE !')
  CASE DEFAULT
!     CALL RGMSG(SUBSTR, RGMLE, 'UNRECOGNIZED N-ARRAY TYPE !')
  END SELECT
  
END SUBROUTINE SCALE_NARRAY
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE CAT_NARRAY(NA, NB)

  IMPLICIT NONE

  ! I/O
  TYPE(NARRAY), INTENT(INOUT) :: NA
  TYPE(NARRAY), INTENT(IN)    :: NB

  ! LOCAL
!  CHARACTER(LEN=*), PARAMETER :: SUBSTR = 'CAT_NARRAY'
  TYPE(NARRAY) :: NH
  INTEGER      :: VTYPE1, VTYPE2
  INTEGER      :: N,M,I

  VTYPE1 = NA%TYPE
  VTYPE2 = NB%TYPE

  IF (VTYPE2 .NE. VTYPE_UNDEF) THEN
!  IF (VTYPE2 == VTYPE_UNDEF) THEN
!     CALL RGMSG(SUBSTR, RGMLE, 'N-ARRAY TO BE APPENDED MUST BE DEFINED !')
!  ELSE
     IF (VTYPE1 == VTYPE_UNDEF) THEN
        CALL COPY_NARRAY(NA, NB)
        RETURN
!     ELSE
!        IF (VTYPE1 /= VTYPE2) THEN
!           CALL RGMSG(SUBSTR, RGMLE, 'N-ARRAYS MUST BE OF SAME TYPE !')
!        END IF
     END IF
!     IF ((NA%N /= 1).OR.(NB%N /= 1)) THEN
!        CALL RGMSG(SUBSTR, RGMLE, 'N-ARRAYS MUST BE 1-DIMENSIONAL !')
!     END IF
  END IF

  N = NA%DIM(1)
  M = NB%DIM(1)
  CALL INIT_NARRAY(NH, 1, (/ N+M /), VTYPE2)

  SELECT CASE(VTYPE2)
  CASE(VTYPE_REAL)
     DO I=1, N
        NH%VR(I) = NA%VR(I) 
     END DO
     DO I=1, M
        NH%VR(N+I) = NB%VR(I) 
     END DO
  CASE(VTYPE_DOUBLE)
     DO I=1, N
        NH%VD(I) = NA%VD(I) 
     END DO
     DO I=1, M
        NH%VD(N+I) = NB%VD(I) 
     END DO
  CASE(VTYPE_INT)
     DO I=1, N
        NH%VI(I) = NA%VI(I) 
     END DO
     DO I=1, M
        NH%VI(N+I) = NB%VI(I) 
     END DO
  CASE(VTYPE_BYTE)
     DO I=1, N
        NH%VB(I) = NA%VB(I) 
     END DO
     DO I=1, M
        NH%VB(N+I) = NB%VB(I) 
     END DO
  CASE(VTYPE_CHAR)
     DO I=1, N
        NH%VC(I) = NA%VC(I) 
     END DO
     DO I=1, M
        NH%VC(N+I) = NB%VC(I) 
     END DO
!  CASE(VTYPE_UNDEF)
     ! NOTHING TO DO FOR UNDEFINED ARRAYS
  CASE DEFAULT
     ! NOTHING TO DO FOR UNRECOGNIZED ARRAYS
  END SELECT  

  CALL COPY_NARRAY(NA, NH)
  ! CLEAN UP
  CALL INIT_NARRAY(NH)

END SUBROUTINE CAT_NARRAY
! ------------------------------------------------------------------

! ------------------------------------------------------------------
RECURSIVE SUBROUTINE QSORT_I(DATA,IDX,ILEFT,IRIGHT)

  IMPLICIT NONE

  ! I/O
  INTEGER (I8), DIMENSION(:), INTENT(INOUT)        :: DATA   ! DATA TO SORT
  INTEGER (I8), DIMENSION(:), POINTER              :: IDX  ! INDEX LIST
  INTEGER (I8),               INTENT(IN), OPTIONAL :: ILEFT, IRIGHT

  ! LOCAL 
  CHARACTER(LEN=*), PARAMETER :: SUBSTR = 'QSORT_I'
  INTEGER (I8)  :: TEMP            ! TEMPORAL DATA
  INTEGER       :: N               ! LENGTH OF LIST
  INTEGER (I8)  :: LEFT, RIGHT
  INTEGER (I8)  :: I, LAST, APU
  INTEGER (I8)  :: TI              ! TEMPORAL INDEX
!  INTEGER       :: STATUS

  ! INIT
  N = SIZE(DATA)
  IF (.NOT.ASSOCIATED(IDX)) THEN
     ALLOCATE(IDX(N))
!     ALLOCATE(IDX(N), STAT=STATUS)
!     CALL ERRMSG(SUBSTR,STATUS,1)
     DO I=1, N
        IDX(I) = I
     END DO
  END IF

  IF (PRESENT(ILEFT)) THEN
     LEFT = ILEFT
  ELSE
     LEFT = 1
  END IF
  IF (PRESENT(IRIGHT)) THEN
     RIGHT = IRIGHT
  ELSE
     RIGHT = N
  END IF

  IF(LEFT >= RIGHT) RETURN

  APU = (LEFT+RIGHT)/2

  TEMP = DATA(LEFT)
  DATA(LEFT) = DATA(APU)
  DATA(APU) = TEMP

  TI = IDX(LEFT)
  IDX(LEFT) = IDX(APU)
  IDX(APU) = TI

  LAST = LEFT

  DO I=LEFT+1,RIGHT
     IF(DATA(I) < DATA(LEFT)) THEN
        LAST = LAST+1
        TEMP = DATA(LAST)
        DATA(LAST) = DATA(I)
        DATA(I) = TEMP

        TI = IDX(LAST)
        IDX(LAST) = IDX(I)
        IDX(I) = TI
     ENDIF
  END DO

  TEMP = DATA(LEFT)
  DATA(LEFT) = DATA(LAST)
  DATA(LAST) = TEMP

  TI = IDX(LEFT)
  IDX(LEFT) = IDX(LAST)
  IDX(LAST) = TI

  CALL QSORT_I(DATA,IDX,LEFT,LAST-1)
  CALL QSORT_I(DATA,IDX,LAST+1,RIGHT)

  RETURN

END SUBROUTINE QSORT_I
! ------------------------------------------------------------------

! ------------------------------------------------------------------
RECURSIVE SUBROUTINE QSORT_B(DATA,IDX,ILEFT,IRIGHT)

  IMPLICIT NONE

  ! I/O
  INTEGER (I4), DIMENSION(:), INTENT(INOUT)        :: DATA   ! DATA TO SORT
  INTEGER (I8), DIMENSION(:), POINTER              :: IDX  ! INDEX LIST
  INTEGER (I8),               INTENT(IN), OPTIONAL :: ILEFT, IRIGHT

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: SUBSTR = 'QSORT_B'
  INTEGER (I4)  :: TEMP            ! TEMPORAL DATA
  INTEGER       :: N               ! LENGTH OF LIST
  INTEGER (I8)  :: LEFT, RIGHT
  INTEGER (I8)  :: I, LAST, APU
  INTEGER (I8)  :: TI              ! TEMPORAL INDEX
!  INTEGER       :: STATUS

  ! INIT
  N = SIZE(DATA)
  IF (.NOT.ASSOCIATED(IDX)) THEN
     ALLOCATE(IDX(N))
!     ALLOCATE(IDX(N), STAT=STATUS)
!     CALL ERRMSG(SUBSTR,STATUS,1)
     DO I=1, N
        IDX(I) = I
     END DO
  END IF

  IF (PRESENT(ILEFT)) THEN
     LEFT = ILEFT
  ELSE
     LEFT = 1
  END IF
  IF (PRESENT(IRIGHT)) THEN
     RIGHT = IRIGHT
  ELSE
     RIGHT = N
  END IF

  IF(LEFT >= RIGHT) RETURN

  APU = (LEFT+RIGHT)/2

  TEMP = DATA(LEFT)
  DATA(LEFT) = DATA(APU)
  DATA(APU) = TEMP

  TI = IDX(LEFT)
  IDX(LEFT) = IDX(APU)
  IDX(APU) = TI

  LAST = LEFT

  DO I=LEFT+1,RIGHT
     IF(DATA(I) < DATA(LEFT)) THEN
        LAST = LAST+1
        TEMP = DATA(LAST)
        DATA(LAST) = DATA(I)
        DATA(I) = TEMP

        TI = IDX(LAST)
        IDX(LAST) = IDX(I)
        IDX(I) = TI
     ENDIF
  END DO

  TEMP = DATA(LEFT)
  DATA(LEFT) = DATA(LAST)
  DATA(LAST) = TEMP

  TI = IDX(LEFT)
  IDX(LEFT) = IDX(LAST)
  IDX(LAST) = TI

  CALL QSORT_B(DATA,IDX,LEFT,LAST-1)
  CALL QSORT_B(DATA,IDX,LAST+1,RIGHT)

  RETURN

END SUBROUTINE QSORT_B
! ------------------------------------------------------------------

! ------------------------------------------------------------------
RECURSIVE SUBROUTINE QSORT_R(DATA,IDX,ILEFT,IRIGHT)

  IMPLICIT NONE

  ! I/O
  REAL (SP),    DIMENSION(:), INTENT(INOUT)        :: DATA   ! DATA TO SORT
  INTEGER (I8), DIMENSION(:), POINTER              :: IDX  ! INDEX LIST
  INTEGER (I8),               INTENT(IN), OPTIONAL :: ILEFT, IRIGHT

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: SUBSTR = 'QSORT_R'
  REAL (SP)    :: TEMP            ! TEMPORAL DATA
  INTEGER      :: N               ! LENGTH OF LIST
  INTEGER (I8) :: LEFT, RIGHT
  INTEGER (I8) :: I, LAST, APU
  INTEGER (I8) :: TI              ! TEMPORAL INDEX
!  INTEGER      :: STATUS

  ! INIT
  N = SIZE(DATA)
  IF (.NOT.ASSOCIATED(IDX)) THEN
     ALLOCATE(IDX(N))
!     ALLOCATE(IDX(N),STAT=STATUS)
!     CALL ERRMSG(SUBSTR,STATUS,1)
     DO I=1, N
        IDX(I) = I
     END DO
  END IF

  IF (PRESENT(ILEFT)) THEN
     LEFT = ILEFT
  ELSE
     LEFT = 1
  END IF
  IF (PRESENT(IRIGHT)) THEN
     RIGHT = IRIGHT
  ELSE
     RIGHT = N
  END IF

  IF(LEFT >= RIGHT) RETURN

  APU = (LEFT+RIGHT)/2

  TEMP = DATA(LEFT)
  DATA(LEFT) = DATA(APU)
  DATA(APU) = TEMP

  TI = IDX(LEFT)
  IDX(LEFT) = IDX(APU)
  IDX(APU) = TI

  LAST = LEFT

  DO I=LEFT+1,RIGHT
     IF(DATA(I) < DATA(LEFT)) THEN
        LAST = LAST+1
        TEMP = DATA(LAST)
        DATA(LAST) = DATA(I)
        DATA(I) = TEMP

        TI = IDX(LAST)
        IDX(LAST) = IDX(I)
        IDX(I) = TI
     ENDIF
  END DO

  TEMP = DATA(LEFT)
  DATA(LEFT) = DATA(LAST)
  DATA(LAST) = TEMP

  TI = IDX(LEFT)
  IDX(LEFT) = IDX(LAST)
  IDX(LAST) = TI

  CALL QSORT_R(DATA,IDX,LEFT,LAST-1)
  CALL QSORT_R(DATA,IDX,LAST+1,RIGHT)

  RETURN

END SUBROUTINE QSORT_R
! ------------------------------------------------------------------

! ------------------------------------------------------------------
RECURSIVE SUBROUTINE QSORT_D(DATA,IDX,ILEFT,IRIGHT)

  IMPLICIT NONE

  ! I/O
  REAL (DP),    DIMENSION(:), INTENT(INOUT)        :: DATA   ! DATA TO SORT
  INTEGER (I8), DIMENSION(:), POINTER              :: IDX  ! INDEX LIST
  INTEGER (I8),               INTENT(IN), OPTIONAL :: ILEFT, IRIGHT

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: SUBSTR = 'QSORT_D'
  REAL (DP)    :: TEMP            ! TEMPORAL DATA
  INTEGER      :: N               ! LENGTH OF LIST
  INTEGER (I8) :: LEFT, RIGHT
  INTEGER (I8) :: I, LAST, APU
  INTEGER (I8) :: TI              ! TEMPORAL INDEX
!  INTEGER      :: STATUS

  ! INIT
  N = SIZE(DATA)
  IF (.NOT.ASSOCIATED(IDX)) THEN
     ALLOCATE(IDX(N))
!     ALLOCATE(IDX(N), STAT=STATUS)
!     CALL ERRMSG(SUBSTR,STATUS,1)
     DO I=1, N
        IDX(I) = I
     END DO
  END IF

  IF (PRESENT(ILEFT)) THEN
     LEFT = ILEFT
  ELSE
     LEFT = 1
  END IF
  IF (PRESENT(IRIGHT)) THEN
     RIGHT = IRIGHT
  ELSE
     RIGHT = N
  END IF

  IF(LEFT >= RIGHT) RETURN

  APU = (LEFT+RIGHT)/2

  TEMP = DATA(LEFT)
  DATA(LEFT) = DATA(APU)
  DATA(APU) = TEMP

  TI = IDX(LEFT)
  IDX(LEFT) = IDX(APU)
  IDX(APU) = TI

  LAST = LEFT

  DO I=LEFT+1,RIGHT
     IF(DATA(I) < DATA(LEFT)) THEN
        LAST = LAST+1
        TEMP = DATA(LAST)
        DATA(LAST) = DATA(I)
        DATA(I) = TEMP

        TI = IDX(LAST)
        IDX(LAST) = IDX(I)
        IDX(I) = TI
     ENDIF
  END DO

  TEMP = DATA(LEFT)
  DATA(LEFT) = DATA(LAST)
  DATA(LAST) = TEMP

  TI = IDX(LEFT)
  IDX(LEFT) = IDX(LAST)
  IDX(LAST) = TI

  CALL QSORT_D(DATA,IDX,LEFT,LAST-1)
  CALL QSORT_D(DATA,IDX,LAST+1,RIGHT)

  RETURN

END SUBROUTINE QSORT_D
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE OVL_RR (SL,SR,DL,DR,FS,FD,MODULO)

  IMPLICIT NONE

  ! I/O
  REAL (SP), INTENT(IN)              :: SL, SR  ! 'BOX' EDGES 
  REAL (SP), INTENT(IN)              :: DL, DR  ! 'BOX' EDGES
  REAL (DP), INTENT(IN),  OPTIONAL   :: MODULO  ! SHIFT FOR 'MODULO' AXIS
  REAL (DP), INTENT(OUT)             :: FS, FD  ! OVERLAP FRACTIONS

  ! LOCAL
  REAL (DP) :: X1, X2
  REAL (DP) :: SHIFT
  INTEGER   :: FX

  X1 = MAX(MIN(DL,DR),MIN(SL,SR))
  X2 = MIN(MAX(DL,DR),MAX(SL,SR))
  FX = INT(SIGN(REAL(-0.5, DP),X2-X1)+REAL(1., DP))
  FS = (X2-X1)*REAL(FX, DP)/ABS(SR-SL)
  FD = (X2-X1)*REAL(FX, DP)/ABS(DR-DL)

  IF (PRESENT(MODULO)) THEN
     
     SHIFT = REAL(INT((DL/MODULO) + REAL(1.5, DP)), DP) * MODULO
     X1 = MAX(DBLE(MIN(DL,DR))+SHIFT,DBLE(MIN(SL,SR)))
     X2 = MIN(DBLE(MAX(DL,DR))+SHIFT,DBLE(MAX(SL,SR)))
     FX = INT(SIGN(REAL(-0.5, DP),X2-X1)+REAL(1., DP))
     FS = FS + (X2-X1)*REAL(FX, DP)/ABS(SR-SL)
     FD = FD + (X2-X1)*REAL(FX, DP)/ABS(DR-DL)

     SHIFT = REAL(INT((SL/MODULO) + REAL(1.5, DP)), DP) * MODULO
     X1 = MAX(DBLE(MIN(DL,DR)),DBLE(MIN(SL,SR))+SHIFT)
     X2 = MIN(DBLE(MAX(DL,DR)),DBLE(MAX(SL,SR))+SHIFT)
     FX = INT(SIGN(REAL(-0.5, DP),X2-X1)+REAL(1., DP))
     FS = FS + (X2-X1)*REAL(FX, DP)/ABS(SR-SL)
     FD = FD + (X2-X1)*REAL(FX, DP)/ABS(DR-DL)
  END IF

END SUBROUTINE OVL_RR
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE OVL_RD (SL,SR,DL,DR,FS,FD,MODULO)

  IMPLICIT NONE

  ! I/O
  REAL (SP), INTENT(IN)             :: SL, SR  ! 'BOX' EDGES
  REAL (DP), INTENT(IN)             :: DL, DR  ! 'BOX' EDGES
  REAL (DP), INTENT(IN),  OPTIONAL  :: MODULO  ! SHIFT FOR 'MODULO' AXIS
  REAL (DP), INTENT(OUT)            :: FS, FD  ! OVERLAP FRACTIONS

  ! LOCAL
  REAL (DP) :: X1, X2
  REAL (DP) :: SHIFT
  INTEGER   :: FX

  X1 = MAX(MIN(DL,DR),DBLE(MIN(SL,SR)))
  X2 = MIN(MAX(DL,DR),DBLE(MAX(SL,SR)))
  FX = INT(SIGN(REAL(-0.5, DP),X2-X1)+REAL(1., DP))
  FS = (X2-X1)*REAL(FX, DP)/ABS(SR-SL)
  FD = (X2-X1)*REAL(FX, DP)/ABS(DR-DL)

  IF (PRESENT(MODULO)) THEN

     SHIFT = REAL(INT((DL/MODULO) + REAL(1.5, DP)), DP) * MODULO
     X1 = MAX(MIN(DL,DR)+SHIFT,DBLE(MIN(SL,SR)))
     X2 = MIN(MAX(DL,DR)+SHIFT,DBLE(MAX(SL,SR)))
     FX = INT(SIGN(REAL(-0.5, DP),X2-X1)+REAL(1., DP))
     FS = FS + (X2-X1)*REAL(FX, DP)/ABS(SR-SL)
     FD = FD + (X2-X1)*REAL(FX, DP)/ABS(DR-DL)

     SHIFT = REAL(INT((SL/MODULO) + REAL(1.5, DP)), DP) * MODULO
     X1 = MAX(MIN(DL,DR),DBLE(MIN(SL,SR))+SHIFT)
     X2 = MIN(MAX(DL,DR),DBLE(MAX(SL,SR))+SHIFT)
     FX = INT(SIGN(REAL(-0.5, DP),X2-X1)+REAL(1., DP))
     FS = FS + (X2-X1)*REAL(FX, DP)/ABS(SR-SL)
     FD = FD + (X2-X1)*REAL(FX, DP)/ABS(DR-DL)
  END IF

END SUBROUTINE OVL_RD
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE OVL_DR (SL,SR,DL,DR,FS,FD,MODULO)

  IMPLICIT NONE

  ! I/O
  REAL (DP), INTENT(IN)             :: SL, SR  ! 'BOX' EDGES
  REAL (SP), INTENT(IN)             :: DL, DR  ! 'BOX' EDGES
  REAL (DP), INTENT(IN),  OPTIONAL  :: MODULO  ! SHIFT FOR 'MODULO' AXIS
  REAL (DP), INTENT(OUT)            :: FS, FD  ! OVERLAP FRACTIONS

  ! LOCAL
  REAL (DP) :: X1, X2
  REAL (DP) :: SHIFT
  INTEGER   :: FX

  X1 = MAX(DBLE(MIN(DL,DR)),MIN(SL,SR))
  X2 = MIN(DBLE(MAX(DL,DR)),MAX(SL,SR))
  FX = INT(SIGN(REAL(-0.5, DP),X2-X1)+REAL(1., DP))
  FS = (X2-X1)*REAL(FX, DP)/ABS(SR-SL)
  FD = (X2-X1)*REAL(FX, DP)/ABS(DR-DL)

  IF (PRESENT(MODULO)) THEN

     SHIFT = REAL(INT((DL/MODULO) + REAL(1.5, DP)), DP) * MODULO
     X1 = MAX(DBLE(MIN(DL,DR))+SHIFT,MIN(SL,SR))
     X2 = MIN(DBLE(MAX(DL,DR))+SHIFT,MAX(SL,SR))
     FX = INT(SIGN(REAL(-0.5, DP),X2-X1)+REAL(1., DP))
     FS = FS + (X2-X1)*REAL(FX, DP)/ABS(SR-SL)
     FD = FD + (X2-X1)*REAL(FX, DP)/ABS(DR-DL)

     SHIFT = REAL(INT((SL/MODULO) + REAL(1.5, DP)), DP) * MODULO
     X1 = MAX(DBLE(MIN(DL,DR)),MIN(SL,SR)+SHIFT)
     X2 = MIN(DBLE(MAX(DL,DR)),MAX(SL,SR)+SHIFT)
     FX = INT(SIGN(REAL(-0.5, DP),X2-X1)+REAL(1., DP))
     FS = FS + (X2-X1)*REAL(FX, DP)/ABS(SR-SL)
     FD = FD + (X2-X1)*REAL(FX, DP)/ABS(DR-DL)
  END IF

END SUBROUTINE OVL_DR
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE OVL_DD (SL,SR,DL,DR,FS,FD,MODULO)

  IMPLICIT NONE

  ! I/O
  REAL (DP), INTENT(IN)             :: SL, SR  ! 'BOX' EDGES
  REAL (DP), INTENT(IN)             :: DL, DR  ! 'BOX' EDGES
  REAL (DP), INTENT(IN),  OPTIONAL  :: MODULO  ! SHIFT FOR 'MODULO' AXIS
  REAL (DP), INTENT(OUT)            :: FS, FD  ! OVERLAP FRACTIONS

  ! LOCAL
  REAL (DP) :: X1, X2
  REAL (DP) :: SHIFT
  INTEGER   :: FX

  X1 = MAX(MIN(DL,DR),MIN(SL,SR))
  X2 = MIN(MAX(DL,DR),MAX(SL,SR))
  FX = INT(SIGN(REAL(-0.5, DP),X2-X1)+REAL(1., DP))
  FS = (X2-X1)*REAL(FX, DP)/ABS(SR-SL)
  FD = (X2-X1)*REAL(FX, DP)/ABS(DR-DL)

  IF (PRESENT(MODULO)) THEN

     SHIFT = REAL(INT((DL/MODULO) + REAL(1.5, DP)), DP) * MODULO
     X1 = MAX(MIN(DL,DR)+SHIFT,MIN(SL,SR))
     X2 = MIN(MAX(DL,DR)+SHIFT,MAX(SL,SR))
     FX = INT(SIGN(REAL(-0.5, DP),X2-X1)+REAL(1., DP))
     FS = FS + (X2-X1)*REAL(FX, DP)/ABS(SR-SL)
     FD = FD + (X2-X1)*REAL(FX, DP)/ABS(DR-DL)

     SHIFT = REAL(INT((SL/MODULO) + REAL(1.5, DP)), DP) * MODULO
     X1 = MAX(MIN(DL,DR),MIN(SL,SR)+SHIFT)
     X2 = MIN(MAX(DL,DR),MAX(SL,SR)+SHIFT)
     FX = INT(SIGN(REAL(-0.5, DP),X2-X1)+REAL(1., DP))
     FS = FS + (X2-X1)*REAL(FX, DP)/ABS(SR-SL)
     FD = FD + (X2-X1)*REAL(FX, DP)/ABS(DR-DL)
  END IF

END SUBROUTINE OVL_DD
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE OVL_1D_RR(S, D, MFS, MFD, LMOD)

  IMPLICIT NONE

  ! I/O
  REAL (SP), DIMENSION(:), INTENT(IN) :: S        ! SOURCE BOUNDS
  REAL (SP), DIMENSION(:), INTENT(IN) :: D        ! DEST. BOUNDS
!  REAL (DP), DIMENSION(:,:), POINTER  :: MFS, MFD ! OVERLAP MATRICES
  REAL (DP), DIMENSION(:,:)           :: MFS, MFD ! OVERLAP MATRICES
  LOGICAL,                 INTENT(IN) :: LMOD     ! MODULO AXIS ?

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: SUBSTR = 'OVL_1D_RR'
  INTEGER   :: N,M    ! DIMENSIONS
  INTEGER   :: I,J    ! COUNTER
!  INTEGER   :: STATUS
  REAL (DP) :: VMOD   ! MODULO VALUE

  ! INIT
  N = SIZE(S)-1
  M = SIZE(D)-1

  ! ALLOCATE MEMORY
!  ALLOCATE(MFS(N,M),STAT=STATUS)
!  CALL ERRMSG('OVL_1D_RR',STATUS,1)
!  ALLOCATE(MFD(M,N),STAT=STATUS)
!  CALL ERRMSG('OVL_1D_RR',STATUS,2)

  ! CHECK CONSISTENCY
  IF ((SIZE(MFS,1) /= N).OR.(SIZE(MFS,2) /= M)) THEN
     CALL RGMSG(SUBSTR, RGMLE, 'INVALID SIZE OF MATRIX !')
  END IF
  IF ((SIZE(MFD,2) /= N).OR.(SIZE(MFD,1) /= M)) THEN
     CALL RGMSG(SUBSTR, RGMLE, 'INVALID SIZE OF MATRIX !')
  END IF

  IF (LMOD) THEN
     VMOD = ABS(S(N+1))+ABS(S(1))
     DO I=1, N
        DO J=1, M
           CALL OVL_RR(S(I),S(I+1),D(J),D(J+1),MFS(I,J),MFD(J,I),VMOD)
        END DO
     END DO
  ELSE
     DO I=1, N
        DO J=1, M
           CALL OVL_RR(S(I),S(I+1),D(J),D(J+1),MFS(I,J),MFD(J,I))
        END DO
     END DO
  END IF

END SUBROUTINE OVL_1D_RR
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE OVL_1D_RD(S, D, MFS, MFD, LMOD)

  IMPLICIT NONE

  ! I/O
  REAL (SP), DIMENSION(:), INTENT(IN) :: S        ! SOURCE BOUNDS
  REAL (DP), DIMENSION(:), INTENT(IN) :: D        ! DEST. BOUNDS
!  REAL (DP), DIMENSION(:,:), POINTER  :: MFS, MFD ! OVERLAP MATRICES
  REAL (DP), DIMENSION(:,:)           :: MFS, MFD ! OVERLAP MATRICES
  LOGICAL,                 INTENT(IN) :: LMOD     ! MODULO AXIS ?

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: SUBSTR = 'OVL_1D_RD'
  INTEGER   :: N,M    ! DIMENSIONS
  INTEGER   :: I,J    ! COUNTER
!  INTEGER   :: STATUS
  REAL (DP) :: VMOD   ! MODULO VALUE

  ! INIT
  N = SIZE(S)-1
  M = SIZE(D)-1

  ! ALLOCATE MEMORY
!  ALLOCATE(MFS(N,M),STAT=STATUS)
!  CALL ERRMSG('OVL_1D_RD',STATUS,1)
!  ALLOCATE(MFD(M,N),STAT=STATUS)
!  CALL ERRMSG('OVL_1D_RD',STATUS,2)

  ! CHECK CONSISTENCY
  IF ((SIZE(MFS,1) /= N).OR.(SIZE(MFS,2) /= M)) THEN
     CALL RGMSG(SUBSTR, RGMLE, 'INVALID SIZE OF MATRIX !')
  END IF
  IF ((SIZE(MFD,2) /= N).OR.(SIZE(MFD,1) /= M)) THEN
     CALL RGMSG(SUBSTR, RGMLE, 'INVALID SIZE OF MATRIX !')
  END IF

  IF (LMOD) THEN
     VMOD = ABS(S(N+1))+ABS(S(1))
     DO I=1, N
        DO J=1, M
           CALL OVL_RD(S(I),S(I+1),D(J),D(J+1),MFS(I,J),MFD(J,I),VMOD)
        END DO
     END DO
  ELSE
     DO I=1, N
        DO J=1, M
           CALL OVL_RD(S(I),S(I+1),D(J),D(J+1),MFS(I,J),MFD(J,I))
        END DO
     END DO
  END IF

END SUBROUTINE OVL_1D_RD
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE OVL_1D_DR(S, D, MFS, MFD, LMOD)

  IMPLICIT NONE

  ! I/O
  REAL (DP), DIMENSION(:), INTENT(IN) :: S        ! SOURCE BOUNDS
  REAL (SP), DIMENSION(:), INTENT(IN) :: D        ! DEST. BOUNDS
!  REAL (DP), DIMENSION(:,:), POINTER  :: MFS, MFD ! OVERLAP MATRICES
  REAL (DP), DIMENSION(:,:)           :: MFS, MFD ! OVERLAP MATRICES
  LOGICAL,                 INTENT(IN) :: LMOD     ! MODULO AXIS ?

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: SUBSTR = 'OVL_1D_DR'
  INTEGER   :: N,M    ! DIMENSIONS
  INTEGER   :: I,J    ! COUNTER
!  INTEGER   :: STATUS
  REAL (DP) :: VMOD   ! MODULO VALUE

  ! INIT
  N = SIZE(S)-1
  M = SIZE(D)-1

  ! ALLOCATE MEMORY
!  ALLOCATE(MFS(N,M),STAT=STATUS)
!  CALL ERRMSG('OVL_1D_DR',STATUS,1)
!  ALLOCATE(MFD(M,N),STAT=STATUS)
!  CALL ERRMSG('OVL_1D_DR',STATUS,2)

  ! CHECK CONSISTENCY
  IF ((SIZE(MFS,1) /= N).OR.(SIZE(MFS,2) /= M)) THEN
     CALL RGMSG(SUBSTR, RGMLE, 'INVALID SIZE OF MATRIX !')
  END IF
  IF ((SIZE(MFD,2) /= N).OR.(SIZE(MFD,1) /= M)) THEN
     CALL RGMSG(SUBSTR, RGMLE, 'INVALID SIZE OF MATRIX !')
  END IF

  IF (LMOD) THEN
     VMOD = ABS(S(N+1))+ABS(S(1))
     DO I=1, N
        DO J=1, M
           CALL OVL_DR(S(I),S(I+1),D(J),D(J+1),MFS(I,J),MFD(J,I),VMOD)
        END DO
     END DO
  ELSE
     DO I=1, N
        DO J=1, M
           CALL OVL_DR(S(I),S(I+1),D(J),D(J+1),MFS(I,J),MFD(J,I))
        END DO
     END DO
  END IF

END SUBROUTINE OVL_1D_DR
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE OVL_1D_DD(S, D, MFS, MFD, LMOD)

  IMPLICIT NONE

  ! I/O
  REAL (DP), DIMENSION(:), INTENT(IN) :: S        ! SOURCE BOUNDS
  REAL (DP), DIMENSION(:), INTENT(IN) :: D        ! DEST. BOUNDS
!  REAL (DP), DIMENSION(:,:), POINTER  :: MFS, MFD ! OVERLAP MATRICES
  REAL (DP), DIMENSION(:,:)           :: MFS, MFD ! OVERLAP MATRICES
  LOGICAL,                 INTENT(IN) :: LMOD     ! MODULO AXIS ?

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: SUBSTR = 'OVL_1D_DD'
  INTEGER   :: N,M    ! DIMENSIONS
  INTEGER   :: I,J    ! COUNTER
!  INTEGER   :: STATUS
  REAL (DP) :: VMOD   ! MODULO VALUE

  ! INIT
  N = SIZE(S)-1
  M = SIZE(D)-1

  ! ALLOCATE MEMORY
!  ALLOCATE(MFS(N,M),STAT=STATUS)
!  CALL ERRMSG('OVL_1D_DD',STATUS,1)
!  ALLOCATE(MFD(M,N),STAT=STATUS)
!  CALL ERRMSG('OVL_1D_DD',STATUS,2)

  ! CHECK CONSISTENCY
  IF ((SIZE(MFS,1) /= N).OR.(SIZE(MFS,2) /= M)) THEN
     CALL RGMSG(SUBSTR, RGMLE, 'INVALID SIZE OF MATRIX !')
  END IF
  IF ((SIZE(MFD,2) /= N).OR.(SIZE(MFD,1) /= M)) THEN
     CALL RGMSG(SUBSTR, RGMLE, 'INVALID SIZE OF MATRIX !')
  END IF

  IF (LMOD) THEN
     VMOD = ABS(S(N+1))+ABS(S(1))
     DO I=1, N
        DO J=1, M
           CALL OVL_DD(S(I),S(I+1),D(J),D(J+1),MFS(I,J),MFD(J,I),VMOD)
        END DO
     END DO
  ELSE
     DO I=1, N
        DO J=1, M
           CALL OVL_DD(S(I),S(I+1),D(J),D(J+1),MFS(I,J),MFD(J,I))
        END DO
     END DO
  END IF

END SUBROUTINE OVL_1D_DD
! ------------------------------------------------------------------

! ------------------------------------------------------------------
INTEGER FUNCTION POSITION(DIM, VEC)

! THIS FUNCTION CALCULATES THE POSITION NUMBER IN A
! 1-D (LINEAR) ARRAY, GIVEN THAT THE ARRAY SHOULD BE INTERPRETED
! AS N-D ARRAY WITH DIMENSIONS
! DIM = (D1, D2, D3, ..., DN)
! OF THE ELEMENT
! VEC = (V1, V2, V3, ..., VN)
!

  IMPLICIT NONE

  ! I/O
  INTEGER, DIMENSION(:), INTENT(IN) :: DIM
  INTEGER, DIMENSION(:), INTENT(IN) :: VEC

  ! LOCAL
  INTEGER :: I
  INTEGER :: N
  INTEGER :: DACC

! Inactive to allow SIZE(VEC) > SIZE(DIM) - NEC jh
! IF (SIZE(DIM) /= SIZE(VEC)) THEN
!    POSITION = 0
!    RETURN
! END IF

  DO I=1, SIZE(DIM)
     IF (VEC(I) > DIM(I)) THEN
        POSITION = 0
        RETURN
     END IF
  END DO

  N = VEC(1)
  DACC = 1
  DO I=2,SIZE(DIM)
     DACC = DACC*DIM(I-1)
     N = N + DACC*(VEC(I)-1)
  END DO

  POSITION = N

END FUNCTION POSITION
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE ELEMENT(DIM, N, VEC)

! THIS SUBROUTINE CALCULATES THE ELEMENT VECTOR
! VEC = (V1, V2, V3, ..., VN)
! OF THE ELEMENT WITH POSITION N IN A
! 1-D (LINEAR) ARRAY, GIVEN THAT THE ARRAY SHOULD BE INTERPRETED
! AS N-D ARRAY WITH DIMENSIONS
! DIM = (D1, D2, D3, ..., DN)
!

  IMPLICIT NONE

  ! I/O
  INTEGER, DIMENSION(:), INTENT(IN)  :: DIM   ! DIMENSION VECTOR
  INTEGER,               INTENT(IN)  :: N     ! ELEMENT IN LINEAR ARRAY
  INTEGER, DIMENSION(:), POINTER     :: VEC   ! ELEMENT VECTOR

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: SUBSTR = 'ELEMENT'
  INTEGER                             :: M    ! COPY OF N
  INTEGER                             :: I    ! COUNTER
  INTEGER , DIMENSION(:), ALLOCATABLE :: DACC
  INTEGER                             :: L    ! LENGTH OF DIM
!  INTEGER                             :: STATUS

  M = N
  L = SIZE(DIM)
  IF (ASSOCIATED(VEC)) THEN
     IF (SIZE(VEC) > 0) THEN
        DEALLOCATE(VEC)
!        DEALLOCATE(VEC, STAT=STATUS)
!        CALL ERRMSG(SUBSTR,STATUS,1)
     END IF
     NULLIFY(VEC)
  END IF
  ALLOCATE(VEC(L))
!  ALLOCATE(VEC(L), STAT=STATUS)
!  CALL ERRMSG(SUBSTR,STATUS,2)
  VEC(:) = 0

  ALLOCATE(DACC(L))
  DACC(1) = 1
  DO I=2, L
     DACC(I) = DACC(I-1)*DIM(I-1)
  END DO

  IF (M > DACC(L)*DIM(L)) RETURN

  DO I=L, 2, -1
     VEC(I) = (M-1)/DACC(I)+1
     M = M - (VEC(I)-1)*DACC(I)
  END DO
  VEC(1) = M

  DEALLOCATE(DACC)
!  DEALLOCATE(DACC, STAT=STATUS)
!  CALL ERRMSG(SUBSTR,STATUS,3)

END SUBROUTINE ELEMENT
! ------------------------------------------------------------------

! ------------------------------------------------------------------
RECURSIVE SUBROUTINE NREGRID(S, SAX, DAX, D, RG_TYPE       &
                             ,SOVL,  DOVL, RCNT            &
                             ,GM, GNR, GSVEC, GDVEC, GDIO  &
                             ,GSF, GDF                     &
                            )

  IMPLICIT NONE

  ! I/O
  TYPE (NARRAY), DIMENSION(:), INTENT(IN)    :: S       ! SOURCE DATA FIELDS
  TYPE (NARRAY), DIMENSION(:), POINTER       :: D       ! DEST. DATA FIELDS
  TYPE (AXIS)  , DIMENSION(:), INTENT(IN)    :: SAX     ! AXES (SOURCE)
  TYPE (AXIS)  , DIMENSION(:), INTENT(INOUT) :: DAX     ! AXES (DEST)
  INTEGER      , DIMENSION(:), INTENT(IN)    :: RG_TYPE ! REGRIDDING TYPE
  REAL (DP),     DIMENSION(:), POINTER       :: SOVL  ! GLOB. OVERLAP FRACTION
  REAL (DP),     DIMENSION(:), POINTER       :: DOVL  ! GLOB. OVERLAP FRACTION
  INTEGER,       DIMENSION(:), POINTER       :: RCNT  ! RECURSION LEVEL COUNTER

  ! FOR RECURSION
  INTEGER,      OPTIONAL, INTENT(IN)            :: GM    ! CURRENT DIMENSION
  INTEGER,      OPTIONAL, INTENT(IN)            :: GNR   ! NO OF DIMS TO REGRID
! INTEGER,      OPTIONAL, DIMENSION(:)          :: GSVEC ! CNT. VECTOR (SOURCE)
! INTEGER,      OPTIONAL, DIMENSION(:)          :: GDVEC ! CNT VECTOR (DEST.)
  INTEGER,      OPTIONAL, DIMENSION(:), TARGET  :: GSVEC ! CNT. VECTOR (SOURCE)
  INTEGER,      OPTIONAL, DIMENSION(:), TARGET  :: GDVEC ! CNT VECTOR (DEST.)
! INTEGER,      OPTIONAL, DIMENSION(:)          :: GDIO  ! DIMENSION ORDER
  INTEGER,      OPTIONAL, DIMENSION(:), TARGET  :: GDIO  ! DIMENSION ORDER
  REAL (DP),    OPTIONAL, INTENT(IN)   :: GSF   ! CURRENT OVERLAP FRACTION
  REAL (DP),    OPTIONAL, INTENT(IN)   :: GDF   ! CURRENT OVERLAP FRACTION

  ! FOR 1ST STEP -> RECURSION
! INTEGER, DIMENSION(:), ALLOCATABLE :: DIO  ! DIMENSION ORDER
  INTEGER, DIMENSION(:), POINTER :: DIO  ! DIMENSION ORDER
  INTEGER                        :: M     ! CURRENT DIMENION
  INTEGER                        :: NR    ! NO OF DIMS TO REGRID
! INTEGER, DIMENSION(:), ALLOCATABLE :: SVEC ! CNT. VECTOR (SOURCE)
! INTEGER, DIMENSION(:), ALLOCATABLE :: DVEC ! CNT. VECTOR (DEST.)
  INTEGER, DIMENSION(:), POINTER :: SVEC ! CNT. VECTOR (SOURCE)
  INTEGER, DIMENSION(:), POINTER :: DVEC ! CNT. VECTOR (DEST.)
  REAL (DP)                      :: SF    ! CURRENT OVERLAP FRACTION
  REAL (DP)                      :: DF    ! CURRENT OVERLAP FRACTION
  REAL (DP)                      :: SF0   ! CURRENT OVERLAP FRACTION
  REAL (DP)                      :: DF0   ! CURRENT OVERLAP FRACTION

  ! FOR INVARIANT DEPENDENT AXIS
  TYPE (AXIS)  , DIMENSION(:), ALLOCATABLE   :: PSAX   ! AXES (SOURCE)
  TYPE (AXIS)  , DIMENSION(:), ALLOCATABLE   :: PDAX   ! AXES (DEST)
  TYPE (NARRAY), DIMENSION(:), POINTER       :: PD     ! DEST. DATA
  REAL (DP),     DIMENSION(:), POINTER       :: PSOVL  ! GLOB. OVERLAP FRACTION
  REAL (DP),     DIMENSION(:), POINTER       :: PDOVL  ! GLOB. OVERLAP FRACTION
  INTEGER,       DIMENSION(:), POINTER       :: PRCNT ! RECURSION LEVEL COUNTER

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: SUBSTR = 'NREGRID'
  INTEGER                                 :: I,J, K ! COUNTER
  INTEGER                                 :: NVAR   ! NUMBER OF FIELDS
  INTEGER                                 :: ID     ! DIMENIONS LOOP OUNTER
  INTEGER                                 :: NID    ! COUNT NO. OF DIMENSIONS
  INTEGER                                 :: NDIM   ! NUMBER OF DIMENSIONS
  REAL (DP) , DIMENSION(:), ALLOCATABLE   :: SB     ! SOURCE BOUNDS
  REAL (DP) , DIMENSION(:), ALLOCATABLE   :: DB     ! DEST. BOUNDS
  INTEGER                                 :: VO     ! OFFSET FOR BOUND VEC.
  INTEGER                                 :: VL     ! LENGTH OF BOUND VEC.
! INTEGER,    DIMENSION(:), ALLOCATABLE   :: VEC    ! ELEMENT VECTOR
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: VEC    ! ELEMENT VECTOR
  REAL (DP),  DIMENSION(:,:), ALLOCATABLE :: MS ! OVERLAP MATRIX
  REAL (DP),  DIMENSION(:,:), ALLOCATABLE :: MD ! OVERLAP MATRIX
  INTEGER                                 :: NS, ND ! DIM. OF OVL-MATRICES
!  INTEGER                                 :: STATUS ! ERROR STATUS
  INTEGER                                 :: VTYPE, SVTYPE  ! VARIABLE TYPE
  LOGICAL                                 :: LFLAG  ! SWITCH
  LOGICAL                                 :: LSDEF  ! FLAG FOR DEFINED DIM.
  LOGICAL                                 :: LDDEF  ! FLAG FOR DEFINED DIM.
  LOGICAL                                 :: LSDEP  ! FLAG FOR DEPENDENT DIM.
  INTEGER                                 :: NOVL   ! NUMBER OF OVERLAPS
  REAL (DP)                               :: ZSO, ZDO ! LOCAL OVERLAP
! INTEGER, DIMENSION(:), ALLOCATABLE      :: HIO    ! FOR OUTPUT
! CHARACTER(LEN=1000)                     :: HSTR   ! FOR OUTPUT

! INTEGER                                 :: NPOS,DACCHELP
! INTEGER                                 :: DDIM,SDIM

  real (dp), allocatable, dimension(:) :: df_, sf_
  integer, allocatable, dimension(:)   :: i_, j_

  ! INIT
  ! NUMBER OF DIMENSIONS
  NDIM = SIZE(SAX)

  IF (PRESENT(GM)) THEN
     M = GM
  ELSE
     M = 1
  END IF

  IF (PRESENT(GNR)) THEN
     NR = GNR
  ELSE
     NR = 0
  END IF

  IF (PRESENT(GSVEC)) THEN
     SVEC => GSVEC
!     ALLOCATE(SVEC(SIZE(GSVEC)), STAT=STATUS)
!     CALL ERRMSG(SUBSTR,STATUS,1)
!     SVEC(:) = GSVEC(:)
  ELSE
     ALLOCATE(SVEC(S(1)%N))
!     ALLOCATE(SVEC(S(1)%N),STAT=STATUS)
!     CALL ERRMSG(SUBSTR,STATUS,2)
     SVEC(:) = 1
  END IF

  IF (PRESENT(GDVEC)) THEN
     DVEC => GDVEC
!     ALLOCATE(DVEC(SIZE(GDVEC)), STAT=STATUS)
!     CALL ERRMSG(SUBSTR,STATUS,3)
!     DVEC(:) = GDVEC(:)
  ELSE
     ALLOCATE(DVEC(S(1)%N))
!     ALLOCATE(DVEC(S(1)%N),STAT=STATUS)
!     CALL ERRMSG(SUBSTR,STATUS,4)
     DVEC(:) = 1
  END IF

  IF (PRESENT(GDIO)) THEN
     DIO => GDIO
!     ALLOCATE(DIO(SIZE(GDIO)), STAT=STATUS)
!     CALL ERRMSG(SUBSTR,STATUS,5)
!     DIO(:) = GDIO(:)
  ELSE
     ALLOCATE(DIO(NDIM))
!     ALLOCATE(DIO(NDIM),STAT=STATUS)
!     CALL ERRMSG(SUBSTR,STATUS,6)
     DIO(:) = 0
  END IF

  IF (PRESENT(GSF)) THEN
     SF0 = GSF
  ELSE
     SF0 = REAL(1.0, DP)
  END IF

  IF (PRESENT(GDF)) THEN
     DF0 = GDF
  ELSE
     DF0 = REAL(1.0, DP)
  END IF

  NVAR = SIZE(S)

  ! ....................................................................
  ! (A) ONLY TO BE DONE ONCE (FIRST STEP)
  IF (M == 1) THEN
     ! 1. CHECK SOME BASIC REQUIREMENTS
!     IF (SIZE(RG_TYPE) /= NVAR) THEN
!       CALL RGMSG(SUBSTR,RGMLE, &
!            'NUMBER OF REGRIDDING TYPE MISMATCH IN SOURCE !')
!     END IF
     !
!     DO K=1, NVAR
!        IF (S(K)%N /= SIZE(SAX)) THEN
!           CALL RGMSG(SUBSTR,RGMLE,'NUMBER OF DIMENSIONS MISMATCH', .FALSE.)
!           CALL RGMSG(SUBSTR,RGMLEC,'FOR VARIABLE ',K,':', .FALSE.)
!           CALL RGMSG(SUBSTR,RGMLEC,'VAR_DIMENSIONS: ',S(K)%N,' ', .FALSE.)
!           CALL RGMSG(SUBSTR,RGMLEC,'AXES          : ',SIZE(SAX),' ')
!        END IF
!     END DO
     !
!     DO K=1, NVAR
!        DO I=1, S(K)%N
!           VTYPE = SAX(I)%DAT%TYPE
!           IF (VTYPE /= VTYPE_UNDEF) THEN
!              IF (S(K)%DIM(I) /= SAX(I)%DAT%DIM(1)-1) THEN
!                 CALL RGMSG(SUBSTR,RGMLE,'SOURCE DIMENSION MISMATCH', &
!                      .FALSE.)
!                 CALL RGMSG(SUBSTR,RGMLEC,'FOR DIMENSION ',I,' ', .FALSE.)
!                 CALL RGMSG(SUBSTR,RGMLEC,'OF VARIABLE ',K,' ', .FALSE.)
!                 CALL RGMSG(SUBSTR,RGMLEC,'BETWEEN DATA WITH LENGTH ', &
!                      S(K)%DIM(I),' ', .FALSE.)
!                 CALL RGMSG(SUBSTR,RGMLEC,'AND AXIS WITH LENGTH ',     &
!                      SAX(I)%DAT%DIM(1),' ')
!              END IF
!           END IF
!        END DO
!     END DO
     !
!     IF (SIZE(DAX) /= SIZE(SAX)) THEN
!        CALL RGMSG(SUBSTR,RGMLE,'NUMBER OF DIMENSIONS OF SOURCE: ',  &
!             SIZE(SAX),' ',.FALSE.)
!        CALL RGMSG(SUBSTR,RGMLEC,'NUMBER OF DIMENSIONS OF DEST. : ',  &
!             SIZE(DAX),' ')
!     END IF
     !

     ! 2. CHECK DIMENSIONS
     NID = 0
!     CALL RGMSG(SUBSTR, RGMLI, 'SCANNING ',NDIM, ' DIMENSIONS ...')
     ! 2.1: SOURCE DIMESNION TYPE
     DO ID=1,NDIM
        ! 2.1.1.: SOURCE DIM MUST BE REAL, DOUBLE, OR UNDEFINED
        svtype = SAX(ID)%DAT%TYPE
        SELECT CASE(svtype)
           CASE(VTYPE_REAL)
!              CALL RGMSG(substr, RGMLIC, '... SOURCE DIMENSION',id,  &
!                   ' IS OF TYPE REAL' )
              lsdef = .true.
           CASE(VTYPE_DOUBLE)
!              CALL RGMSG(substr, RGMLIC, '... SOURCE DIMENSION',id,  &
!                   ' IS OF TYPE DOUBLE PRECISION' )
              lsdef = .true.
           CASE(VTYPE_UNDEF)
!              CALL RGMSG(substr, RGMLIC, '... SOURCE DIMENSION',id,  &
!                   ' IS UNDEFINED' )
              lsdef = .false.
           CASE DEFAULT
              !CASE(VTYPE_INT)
              !CASE(VTYPE_CHAR)
              !CASE(VTYPE_BYTE)
!              CALL RGMSG(substr, RGMLE, 'REGRIDDING IS ONLY POSSIBLE FOR', &
!                   .false.)
!              CALL RGMSG(substr, RGMLEC,                                  &
!                   'DIMENSIONS OF TYPE REAL OR DOUBLE PRECISION !')
        END SELECT
        ! 2.1.2 CHECK DEPENDENCIES
        IF (sax(id)%ndp > 1) THEN
           lsdep = .true.
!           DO i=1, sax(id)%ndp
!              IF (sax(sax(id)%dep(i))%dat%type == VTYPE_UNDEF) THEN
!                 CALL RGMSG(substr, RGMLE, 'SOURCE DIMENSION ', id,     &
!                      ' IS DEPENDENT ON', .false.)
!                 CALL RGMSG(substr, RGMLEC, 'DIMENSION ',sax(id)%dep(i), &
!                      ' WHICH IS UNDEFINED !')
!              END IF
!           END DO
        ELSE
           lsdep = .false.
        END IF
        ! 2.2.1.: DEST. DIM MUST BE REAL, DOUBLE, OR UNDEFINED
        vtype = dax(id)%dat%type
        SELECT CASE(vtype)
           CASE(VTYPE_REAL)
!              CALL RGMSG(substr, RGMLIC, '... DEST.  DIMENSION',id,  &
!                   ' IS OF TYPE REAL' )
              lddef = .true.
           CASE(VTYPE_DOUBLE)
!              CALL RGMSG(substr, RGMLIC, '... DEST.  DIMENSION',id,  &
!                   ' IS OF TYPE DOUBLE PRECISION')
              lddef = .true.
           CASE(VTYPE_UNDEF)
!              CALL RGMSG(substr, RGMLIC, '... DEST.  DIMENSION',id,  &
!                   ' IS UNDEFINED')
              lddef = .false.
           CASE DEFAULT
              !CASE(VTYPE_INT)
              !CASE(VTYPE_CHAR)
              !CASE(VTYPE_BYTE)
!              CALL RGMSG(substr, RGMLE, 'REGRIDDING IS ONLY POSSIBLE FOR', &
!                   .false.)
!              CALL RGMSG(substr, RGMLEC,                                 &
!                   'DIMENSIONS OF TYPE REAL OR DOUBLE PRECISION !')
        END SELECT
        ! 2.2.2 CHECK DEPENDENCIES
!        IF (dax(id)%ndp > 1) THEN
!           DO i=1, dax(id)%ndp
!              IF (dax(dax(id)%dep(i))%dat%type == VTYPE_UNDEF) THEN
!                 CALL RGMSG(substr, RGMLE, 'DEST.  DIMENSION ', id,     &
!                      ' IS DEPENDENT ON', .false.)
!                 CALL RGMSG(substr, RGMLEC, 'DIMENSION ',dax(id)%dep(i), &
!                      ' WHICH IS UNDEFINED !')
!              END IF
!           END DO
!        END IF
        !
        ! 2.3 CHECK CONSISTENCY BETWEEN SOURCE AND DEST.
        ! 2.3.1 UNDEF dimensions cannot be regridded
!        IF ((.NOT.lsdef).AND.(lddef)) THEN
!           CALL RGMSG(substr, RGMLE,  &
!                'DIMENSION OF TYPE UNDEFINED', .false.)
!           CALL RGMSG(substr, RGMLEC, &
!                'CANNOT BE REGRIDDED ! FOR INVARIANT', .false.)
!           CALL RGMSG(substr, RGMLEC, &
!                'DIMENSION, DEST. DIMENSION MUST ALSO', .false.)
!           CALL RGMSG(substr, RGMLEC, 'BE UNDEFINED !')
!        END IF
        ! 2.3.2 INVARIANT, DEPENDENT DIMENSIONS HAVE TO BE PRE-REGRIDDED
        IF (lsdep.AND.(.NOT.lddef)) THEN
!           CALL RGMSG(substr, RGMLI, 'FOUND INVRIANT DEPENDENT DIMENSION')
!           CALL RGMSG(substr, RGMLIC, &
!                ' >>> START REGRIDDING OF DEPENDENT DIMENSION ...')
           ! REGRIDDING AXES
           dax(id)%lm = sax(id)%lm
           dax(id)%ndp = sax(id)%ndp
           ALLOCATE(dax(id)%dep(dax(id)%ndp))
!           ALLOCATE(dax(id)%dep(dax(id)%ndp),STAT=status)
!           CALL ERRMSG(substr,status,7)
           dax(id)%dep(:) = sax(id)%dep(:)
           ALLOCATE(psax(sax(id)%ndp))
!           ALLOCATE(psax(sax(id)%ndp),STAT=status)
!           CALL ERRMSG(substr,status,8)
           ALLOCATE(pdax(dax(id)%ndp))
!           ALLOCATE(pdax(dax(id)%ndp),STAT=status)
!           CALL ERRMSG(substr,status,9)
           CALL INIT_AXIS(psax(1)) ! FIRST DIMENSION IS INVARIANT ...
           CALL INIT_AXIS(pdax(1)) ! ...
           DO i=2, sax(id)%ndp
              CALL COPY_AXIS(psax(i), sax(sax(id)%dep(i)))
              CALL COPY_AXIS(pdax(i), dax(dax(id)%dep(i)))
           END DO
           !
           CALL NREGRID( (/sax(id)%dat/) ,psax, pdax, pd  &
                         ,(/ RG_INT/) , psovl, pdovl, prcnt)
!           CALL NREGRID_STAT(psax, pdax, psovl, pdovl     &
!                         ,(/sax(id)%dat/), pd, prcnt)
           !
           CALL COPY_NARRAY(dax(id)%dat, pd(1))
           DEALLOCATE(psax, pdax, psovl, pdovl)
!           DEALLOCATE(psax, pdax, psovl, pdovl, STAT=status)
!           CALL ERRMSG(substr,status,10)
           NULLIFY(psovl, pdovl)
           DEALLOCATE(pd)
!           DEALLOCATE(pd, STAT=status)
!           CALL ERRMSG(substr,status,11)
           NULLIFY(pd)
           DEALLOCATE(prcnt)
!           DEALLOCATE(prcnt, STAT=status)
!           CALL ERRMSG(substr,status,12)
           NULLIFY(prcnt)
!           CALL RGMSG(substr, RGMLIC, &
!                ' <<< ... END REGRIDDING OF DEPENDENT DIMENSION !')
        END IF

     END DO
!     CALL RGMSG(substr, RGMLIC, '... END SCANNING ',ndim, ' DIMENSIONS !')

     ! 3. CALCULATE REGRIDDING ORDER OF DIMENSIONS
!     CALL RGMSG(substr, RGMLI, 'PROCESSING ORDER OF ',ndim, ' DIMENSIONS ...')
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
!           CALL RGMSG(substr, RGMLIC, ' ... DIMENSION ',id, ' IS INDEPENDENT')
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
!           CALL RGMSG(substr, RGMLIC, ' ... DIMENSION ',id, ' IS DEPENDENT')
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
!           CALL RGMSG(substr, RGMLIC, ' ... DIMENSION ',id, ' IS INVARIANT')
        END IF
     END DO
!     CALL RGMSG(substr, RGMLIC, &
!          '... END PROCESSING ORDER OF ',ndim, ' DIMENSIONS !')
!     IF (nid < s(1)%n) THEN
!        CALL RGMSG(substr, RGMLE, 'UNRECOGNIZED DIMENSION !', .false.)
!        CALL RGMSG(substr, RGMLEC, 'CHECK TYPE OF BOUNDS', .false.)
!        CALL RGMSG(substr, RGMLEC, '(MUST BE REAL OR DOUBLE PRECISION) !')
!     END IF
!     IF (nid > s(1)%n) THEN
!        CALL RGMSG(substr, RGMLE, 'AMBIGIOUS DIMENSION !')
!     END IF

     ! 4. ALLOCATE SPACE FOR REGRIDDED DATA
     ALLOCATE(d(nvar))
     DO k=1, nvar                  ! LOOP OVER VARIABLES
        nid = 1
        d(k)%n = s(k)%n            ! COPY NUMBER OF DIMENSIONS
        d(k)%type = s(k)%type
        ALLOCATE(d(k)%dim(d(k)%n))
!        ALLOCATE(d(k)%dim(d(k)%n),STAT=status)
!        CALL ERRMSG(substr,status,13)
! diese schleife verdoppeln
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
           END IF                               ! dest. dim or source dim.
        END DO

        DO id=1, s(k)%n
           vtype = dax(id)%dat%type
           IF ((vtype == VTYPE_REAL).OR.(vtype == VTYPE_DOUBLE)) THEN
              ! new destination dimension
              ! INTERFACES HAVE +1 ENTRY
              ! 1st DIMENSION MUST BE THE INDEPENDENT !!!
!              d(k)%dim(id) = (dax(id)%dat%dim(1) -1)
!              nid = nid * d(k)%dim(id)
              ! COPY MODULO INFO
!              dax(id)%lm  = sax(id)%lm
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
           ALLOCATE(d(k)%vr(nid))
!           ALLOCATE(d(k)%vr(nid),STAT=status)
!           CALL ERRMSG(substr,status,14)
           d(k)%vr(:) = REAL(0.0, SP)
           d(k)%type = VTYPE_REAL
        CASE(VTYPE_DOUBLE)
           ALLOCATE(d(k)%vd(nid))
!           ALLOCATE(d(k)%vd(nid),STAT=status)
!           CALL ERRMSG(substr,status,15)
           d(k)%vd(:) = REAL(0.0, DP)
           d(k)%type = VTYPE_DOUBLE
        CASE DEFAULT
           !CASE(VTYPE_INT)
           !CASE(VTYPE_CHAR)
           !CASE(VTYPE_BYTE)
           !CASE(VTYPE_UNDEF)
!           CALL RGMSG(substr, RGMLE, 'REGRIDDING IS ONLY POSSIBLE FOR', &
!                .false.)
!           CALL RGMSG(substr, RGMLEC,                                 &
!                'DATA OF TYPE REAL OR DOUBLE PRECISION !')
        END SELECT

     END DO  ! LOOP OVER VARIABLES

     ! 5. SOME DIAGNOSTIC OUTPUT
!     CALL RGMSG(substr, RGMLVM, '    DIMENSIONS TO REGRID          : ', &
!          dio(1:nr),' ')
     !
!     ALLOCATE(hio(nr),STAT=status)
!     CALL ERRMSG(substr,status,16)
!     DO i=1, nr
!        hio(i) = s(1)%dim(dio(i))
!     END DO
!     CALL RGMSG(substr, RGMLVM, '    - LENGTH(S) IN SOURCE         : ', &
!          hio,' ')
     !
!     DO i=1, nr
!        hio(i) = d(1)%dim(dio(i))
!     END DO
!     CALL RGMSG(substr, RGMLVM, '    - LENGTH(S) IN DEST.          : ', &
!          hio,' ')
!     DEALLOCATE(hio,STAT=status)
!     CALL ERRMSG(substr,status,17)
     !
!     CALL RGMSG(substr, RGMLVM, '    INVARIANT DIMENSIONS          : ', &
!          dio(nr+1:ndim), ' ')
!     CALL RGMSG(substr, RGMLVM, '    NUMBER OF FIELDS              : ', &
!          nvar,' ')
     !
!     CALL RGMSG(substr, RGMLVM, '    INVARIANT DIMENSION LENGTH(S) : ')
!     ALLOCATE(hio(ndim-nr),STAT=status)
!     CALL ERRMSG(substr,status,18)
!     DO k=1,nvar
!        WRITE(hstr, *) '       FIELD ',k,': '
!        DO i=nr+1,ndim
!           hio(i-nr) = s(k)%dim(dio(i))
!        END DO
!        CALL RGMSG(substr, RGMLVM, TRIM(hstr), hio,' ')
!     END DO
!     DEALLOCATE(hio,STAT=status)
!     CALL ERRMSG(substr,status,18)

     ! 6. ALLOCATE SPACE FOR DIAGNOSTIC OUTPUT AND INITIALIZE
     ALLOCATE(sovl(s(1)%n))
!     ALLOCATE(sovl(s(1)%n),STAT=status)
!     CALL ERRMSG(substr,status,19)
     sovl(:) = REAL(0.0, DP)
     ALLOCATE(dovl(d(1)%n))
!     ALLOCATE(dovl(d(1)%n),STAT=status)
!     CALL ERRMSG(substr,status,20)
     dovl(:) = REAL(0.0, DP)
     ALLOCATE(rcnt(ndim))
!     ALLOCATE(rcnt(ndim),STAT=status)
!     CALL ERRMSG(substr,status,21)
     rcnt(:) = 0

!     CALL RGMSG(substr, RGMLVL, '    STARTING NREGRID ... ')
     allocate(vec(16))  ! NEC jh

  END IF ! (A) ONLY TO BE DONE ONCE AT FIRST STEP
  ! ....................................................................

  ! SET CURRENT DIMENSION
  nid = dio(m)
  rcnt(m) = rcnt(m) + 1

  ! ....................................................................
  ! (B) CONDITION FOR END OF RECURSION
  IF (m == SIZE(dio)) THEN   ! LAST (INVARIANT !) DIMENSION REACHED
!!$     zso = REAL(0.0, DP)
!!$     zdo = REAL(0.0, DP)
!!$     DO k=1, nvar    ! LOOP OVER VARIABLES
!!$
!!$       vtype = s(k)%type
!!$
!!$       IF ( (SIZE(d(k)%dim).GT.4) .OR. (SIZE(s(k)%dim).GT.4) ) THEN
!!$
!!$         SELECT CASE(RG_TYPE(k))
!!$         CASE(RG_EXT)
!!$           DO j=1, s(k)%dim(nid)
!!$              dvec(nid) = j
!!$              svec(nid) = j     ! INVARIANT !!!
!!$              ! CALCULATE DEST. DATA POINT
!!$              vo = POSITION(d(k)%dim,dvec)   ! POSITION IN DEST.
!!$              vl = POSITION(s(k)%dim,svec)   ! POSITION IN SOURCE
!!$              IF(vtype.EQ.VTYPE_REAL) THEN
!!$                 d(k)%vr(vo) = d(k)%vr(vo) + s(k)%vr(vl) * REAL(sf0, SP)
!!$              ELSE
!!$                 d(k)%vd(vo) = d(k)%vd(vo) + s(k)%vd(vl) * sf0
!!$              END IF
!!$              ! OVERLAP IS 1 FOR INVARIANT DIMENSIONS
!!$              ! INTERVAL LENGTH, NUMBER OF VARIABLES
!!$              zso = zso + REAL(1.0, DP)/ &
!!$                   (REAL(s(k)%dim(nid), DP)*REAL(nvar, DP))
!!$              zdo = zdo + REAL(1.0, DP)/ &
!!$                   (REAL(s(k)%dim(nid), DP)*REAL(nvar, DP))
!!$           END DO
!!$           ! RESET FOR NEXT/PREVIOUS RECURSION STEP
!!$           dvec(nid) = 1
!!$           svec(nid) = 1
!!$         CASE(RG_INT, RG_IDX, RG_IXF)
!!$           !CASE(RG_IDX, RG_IXF)
!!$           ! NOTE: IDX VAR WAS CONVERTET TO INDEX-FRACTION
!!$           DO j=1, s(k)%dim(nid)
!!$              dvec(nid) = j
!!$              svec(nid) = j     ! INVARIANT !!!
!!$              ! CALCULATE DEST. DATA POINT
!!$              vo = POSITION(d(k)%dim,dvec)   ! POSITION IN DEST.
!!$              vl = POSITION(s(k)%dim,svec)   ! POSITION IN SOURCE
!!$              IF(vtype.EQ.VTYPE_REAL) THEN
!!$                 d(k)%vr(vo) = d(k)%vr(vo) + s(k)%vr(vl) * REAL(df0, SP)
!!$              ELSE
!!$                 d(k)%vd(vo) = d(k)%vd(vo) + s(k)%vd(vl) * df0
!!$              END IF
!!$              ! OVERLAP IS 1 FOR INVARIANT DIMENSIONS:
!!$              ! INTERVAL LENGTH, NUMBER OF VARIABLES
!!$              zso = zso + REAL(1.0, DP)/ &
!!$                   (REAL(s(k)%dim(nid), DP)*REAL(nvar, DP))
!!$              zdo = zdo + REAL(1.0, DP)/ &
!!$                   (REAL(s(k)%dim(nid), DP)*REAL(nvar, DP))
!!$           END DO
!!$           ! RESET FOR NEXT/PREVIOUS RECURSION STEP
!!$           dvec(nid) = 1
!!$           svec(nid) = 1
!!$         !CASE DEFAULT
!!$         END SELECT
!!$
!!$       ELSE
!!$
!!$         SELECT CASE(RG_TYPE(k))
!!$         CASE(RG_EXT)
!!$           ddim=SIZE(d(k)%dim)
!!$           sdim=SIZE(s(k)%dim)
!!$           DO j=1, s(k)%dim(nid)
!!$              dvec(nid) = j
!!$              svec(nid) = j     ! INVARIANT !!!
!!$              ! CALCULATE DEST. DATA POINT
!!$              npos = dvec(1)
!!$              SELECT CASE(ddim)
!!$              CASE(2) ! size(dim)=2
!!$                dacchelp = 1
!!$                dacchelp = dacchelp*d(k)%dim(1)
!!$                npos     = npos + dacchelp*(dvec(2)-1)
!!$              CASE(3) ! size(dim)=3
!!$                dacchelp = 1
!!$                dacchelp = dacchelp*d(k)%dim(1)
!!$                npos     = npos + dacchelp*(dvec(2)-1)
!!$                dacchelp = dacchelp*d(k)%dim(2)
!!$                npos     = npos + dacchelp*(dvec(3)-1)
!!$              CASE(4) ! size(dim)=4
!!$                dacchelp = 1
!!$                dacchelp = dacchelp*d(k)%dim(1)
!!$                npos     = npos + dacchelp*(dvec(2)-1)
!!$                dacchelp = dacchelp*d(k)%dim(2)
!!$                npos     = npos + dacchelp*(dvec(3)-1)
!!$                dacchelp = dacchelp*d(k)%dim(3)
!!$                npos     = npos + dacchelp*(dvec(4)-1)
!!$              END SELECT
!!$              vo = npos   ! POSITION IN DEST.
!!$
!!$              npos = svec(1)
!!$              SELECT CASE(sdim)
!!$              CASE(2) ! size(dim)=2
!!$                dacchelp = 1
!!$                dacchelp = dacchelp*s(k)%dim(1)
!!$                npos     = npos + dacchelp*(svec(2)-1)
!!$              CASE(3) ! size(dim)=3
!!$                dacchelp = 1
!!$                dacchelp = dacchelp*s(k)%dim(1)
!!$                npos     = npos + dacchelp*(svec(2)-1)
!!$                dacchelp = dacchelp*s(k)%dim(2)
!!$                npos     = npos + dacchelp*(svec(3)-1)
!!$              CASE(4) ! size(dim)=4
!!$                dacchelp = 1
!!$                dacchelp = dacchelp*s(k)%dim(1)
!!$                npos     = npos + dacchelp*(svec(2)-1)
!!$                dacchelp = dacchelp*s(k)%dim(2)
!!$                npos     = npos + dacchelp*(svec(3)-1)
!!$                dacchelp = dacchelp*s(k)%dim(3)
!!$                npos     = npos + dacchelp*(svec(4)-1)
!!$              END SELECT
!!$              vl = npos   ! POSITION IN SOURCE
!!$
!!$              IF(vtype.EQ.VTYPE_REAL) THEN
!!$                 d(k)%vr(vo) = d(k)%vr(vo) + s(k)%vr(vl) * REAL(sf0, SP)
!!$              ELSE
!!$                 d(k)%vd(vo) = d(k)%vd(vo) + s(k)%vd(vl) * sf0
!!$              END IF
!!$              ! OVERLAP IS 1 FOR INVARIANT DIMENSIONS
!!$              ! INTERVAL LENGTH, NUMBER OF VARIABLES
!!$              zso = zso + REAL(1.0, DP)/ &
!!$                   (REAL(s(k)%dim(nid), DP)*REAL(nvar, DP))
!!$              zdo = zdo + REAL(1.0, DP)/ &
!!$                   (REAL(s(k)%dim(nid), DP)*REAL(nvar, DP))
!!$           END DO
!!$           ! RESET FOR NEXT/PREVIOUS RECURSION STEP
!!$           dvec(nid) = 1
!!$           svec(nid) = 1
!!$         CASE(RG_INT, RG_IDX, RG_IXF)
!!$           !CASE(RG_IDX, RG_IXF)
!!$           ! NOTE: IDX VAR WAS CONVERTET TO INDEX-FRACTION
!!$           ddim=SIZE(d(k)%dim)
!!$           sdim=SIZE(s(k)%dim)
!!$           DO j=1, s(k)%dim(nid)
!!$              dvec(nid) = j
!!$              svec(nid) = j     ! INVARIANT !!!
!!$              ! CALCULATE DEST. DATA POINT
!!$
!!$              npos = dvec(1)
!!$              SELECT CASE(ddim)
!!$              CASE(2) ! size(dim)=2
!!$                dacchelp = 1
!!$                dacchelp = dacchelp*d(k)%dim(1)
!!$                npos     = npos + dacchelp*(dvec(2)-1)
!!$              CASE(3) ! size(dim)=3
!!$                dacchelp = 1
!!$                dacchelp = dacchelp*d(k)%dim(1)
!!$                npos     = npos + dacchelp*(dvec(2)-1)
!!$                dacchelp = dacchelp*d(k)%dim(2)
!!$                npos     = npos + dacchelp*(dvec(3)-1)
!!$              CASE(4) ! size(dim)=4
!!$                dacchelp = 1
!!$                dacchelp = dacchelp*d(k)%dim(1)
!!$                npos     = npos + dacchelp*(dvec(2)-1)
!!$                dacchelp = dacchelp*d(k)%dim(2)
!!$                npos     = npos + dacchelp*(dvec(3)-1)
!!$                dacchelp = dacchelp*d(k)%dim(3)
!!$                npos     = npos + dacchelp*(dvec(4)-1)
!!$              END SELECT
!!$              vo = npos   ! POSITION IN DEST.
!!$
!!$              npos = svec(1)
!!$              SELECT CASE(sdim)
!!$              CASE(2) ! size(dim)=2
!!$                dacchelp = 1
!!$                dacchelp = dacchelp*s(k)%dim(1)
!!$                npos     = npos + dacchelp*(svec(2)-1)
!!$              CASE(3) ! size(dim)=3
!!$                dacchelp = 1
!!$                dacchelp = dacchelp*s(k)%dim(1)
!!$                npos     = npos + dacchelp*(svec(2)-1)
!!$                dacchelp = dacchelp*s(k)%dim(2)
!!$                npos     = npos + dacchelp*(svec(3)-1)
!!$              CASE(4) ! size(dim)=4
!!$                dacchelp = 1
!!$                dacchelp = dacchelp*s(k)%dim(1)
!!$                npos     = npos + dacchelp*(svec(2)-1)
!!$                dacchelp = dacchelp*s(k)%dim(2)
!!$                npos     = npos + dacchelp*(svec(3)-1)
!!$                dacchelp = dacchelp*s(k)%dim(3)
!!$                npos     = npos + dacchelp*(svec(4)-1)
!!$              END SELECT
!!$              vl = npos   ! POSITION IN SOURCE
!!$
!!$              IF(vtype.EQ.VTYPE_REAL) THEN
!!$                 d(k)%vr(vo) = d(k)%vr(vo) + s(k)%vr(vl) * REAL(df0, SP)
!!$              ELSE
!!$                 d(k)%vd(vo) = d(k)%vd(vo) + s(k)%vd(vl) * df0
!!$              END IF
!!$              ! OVERLAP IS 1 FOR INVARIANT DIMENSIONS:
!!$              ! INTERVAL LENGTH, NUMBER OF VARIABLES
!!$              zso = zso + REAL(1.0, DP)/ &
!!$                   (REAL(s(k)%dim(nid), DP)*REAL(nvar, DP))
!!$              zdo = zdo + REAL(1.0, DP)/ &
!!$                   (REAL(s(k)%dim(nid), DP)*REAL(nvar, DP))
!!$           END DO
!!$           ! RESET FOR NEXT/PREVIOUS RECURSION STEP
!!$           dvec(nid) = 1
!!$           svec(nid) = 1
!!$         !CASE DEFAULT
!!$         END SELECT
!!$
!!$       END IF
!!$
!!$     END DO  ! LOOP OVER VARIABLES
!!$
!!$     ! SUM OVERLAP
!!$     sovl(nid) = sovl(nid) + zso
!!$     dovl(nid) = dovl(nid) + zdo

! Because of my introduction of NREGRID_last this sequence will be passed
! for m == 1 only - NEC jh

     allocate(i_(1),j_(1),df_(1),sf_(1))
     i_(1) = 1
     j_(1) = 1
     sf_(1) = sf0
     df_(1) = df0
     CALL NREGRID_last(s, sax, d, RG_TYPE, sovl, dovl        &
                  ,m, svec, dvec, dio, 1, i_, j_, sf_, df_)
     deallocate(i_,j_,df_,sf_,vec)

     ! CLEAN UP
!!$     DEALLOCATE(dio, STAT=status)
!!$     CALL ERRMSG(substr,status,22)
!!$     DEALLOCATE(svec, dvec, STAT=status)
!!$     CALL ERRMSG(substr,status,23)
     IF ( .not. PRESENT(GDIO)) THEN
       DEALLOCATE(dio)
     endif
     IF ( .not. PRESENT(GSVEC)) THEN
       DEALLOCATE(svec)
     endif
     IF ( .not. PRESENT(GDVEC)) THEN
       DEALLOCATE(dvec)
     endif

     RETURN        ! DONE ... ... END RECURSION

  END IF ! (B) CONDITION FOR END OF RECURSION
  ! ....................................................................

  ! ....................................................................
  ! (C) RECURSION STEP
  IF (m <= nr) THEN           ! (C1): DIMENSION TO REGRID

     ! GET SOURCE BOUNDS
     ! 1st dimension is independent -> length of bounds vector
     vl = sax(nid)%dat%dim(1)        ! length is length of 1st (indep.) dim.
     ALLOCATE(sb(vl))
!     ALLOCATE(sb(vl),STAT=status)
!     CALL ERRMSG(substr,status,24)

!     ALLOCATE(vec(sax(nid)%ndp),STAT=status)  ! dependency element vector
!     CALL ERRMSG(substr,status,25)
! NEC jh
     if ( size (vec) .lt. sax(nid)%ndp ) then
       deallocate(vec)
       allocate ( vec(int(real(sax(nid)%ndp)*1.5)))
     endif
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
!        CALL RGMSG(substr, RGMLE, 'SOURCE BOUNDS MUST BE OF TYPE', .false.)
!        CALL RGMSG(substr, RGMLEC, 'REAL OR DOUBLE PRECISION !')
     END SELECT
     !
!     DEALLOCATE(vec,STAT=status)
!     CALL ERRMSG(substr,status,26)

     ! GET DEST. BOUNDS
     ! 1st dimension is independent -> length of bounds vector
     vl = dax(nid)%dat%dim(1)        ! length is length of 1st (indep.) dim.
     ALLOCATE(db(vl))
!     ALLOCATE(db(vl),STAT=status)
!     CALL ERRMSG(substr,status,27)

!     ALLOCATE(vec(dax(nid)%ndp),STAT=status)  ! dependency element vector
!     CALL ERRMSG(substr,status,28)
! NEC jh
     if ( size (vec) .lt. dax(nid)%ndp ) then
       deallocate(vec)
       allocate ( vec(int(real(dax(nid)%ndp)*1.5)))
     endif
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
!        CALL RGMSG(substr, RGMLE, 'DEST. BOUNDS MUST BE OF TYPE', .false.)
!        CALL RGMSG(substr, RGMLEC, 'REAL OR DOUBLE PRECISION !')
     END SELECT
     !
!     DEALLOCATE(vec,STAT=status)
!     CALL ERRMSG(substr,status,29)

     ! REGRID ALONG THIS DIMENSION
     ! ALLOCATE SPACE FOR MATRICES
     ns = SIZE(sb)-1
     nd = SIZE(db)-1
     ALLOCATE(ms(ns,nd))
!     ALLOCATE(ms(ns,nd), STAT=status)
!     CALL ERRMSG(substr,status,30)
     ALLOCATE(md(nd,ns))
!     ALLOCATE(md(nd,ns), STAT=status)
!     CALL ERRMSG(substr,status,31)
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
!!$     DO i=1, ns
!!$        svec(nid) = i
!!$        DO j=1, nd
!!$           dvec(nid) = j
!!$           sf = sf0*ms(svec(nid),dvec(nid))
!!$           df = df0*md(dvec(nid),svec(nid))
!!$           IF (sf > 0.0) THEN    ! NON-ZERO
!!$              novl = novl + 1
!!$              ! CALCULATE TOTAL OVERLAP FRACTION FOR THIS DIMENSION
!!$              sovl(nid) = sovl(nid) +                 &
!!$                   ms(svec(nid),dvec(nid))            &
!!$                   *(sb(svec(nid)+1)-sb(svec(nid)))   &
!!$                   /(sb(ns+1)-sb(1))
!!$              dovl(nid) = dovl(nid) +                 &
!!$                   md(dvec(nid),svec(nid))            &
!!$                   *(db(dvec(nid)+1)-db(dvec(nid)))   &
!!$                   /(db(nd+1)-db(1))
!!$              ! REGRID NEXT DIMENSION
!!$              CALL NREGRID(s, sax, dax, d, RG_TYPE, sovl, dovl, rcnt  &
!!$                   ,m+1, nr, svec, dvec                               &
!!$                   ,dio, sf, df)
!!$           END IF
!!$        END DO
!!$     END DO

     allocate(i_(ns*nd),j_(ns*nd))

     novl = 0
     DO j=1, nd
       DO i=1, ns
         if ( sf0*ms(i,j) > 0.0) then
           novl = novl + 1
           i_(novl) = i
           j_(novl) = j
         endif
       enddo
     enddo

     allocate(df_(novl),sf_(novl))

     do k = 1, novl
       i = i_(k)
       j = j_(k)
       sf = sf0*ms(i,j)
       df = df0*md(j,i)
       sf_(k) = sf
       df_(k) = df
       ! CALCULATE TOTAL OVERLAP FRACTION FOR THIS DIMENSION
       sovl(nid) = sovl(nid) +                 &
                   ms(i,j)                            &
                   *(sb(        i+1)-sb(        i))   &
                   /(sb(ns+1)-sb(1))
       dovl(nid) = dovl(nid) +                 &
                   md(j,i)                            &
                   *(db(        j+1)-db(        j))   &
                   /(db(nd+1)-db(1))
     enddo

     IF (m == SIZE(dio)-1 ) then
       CALL NREGRID_last(s, sax, d, RG_TYPE, sovl, dovl        &
                  ,m+1, svec, dvec, dio, novl, i_, j_, sf_, df_ )
       rcnt(m+1) = rcnt(m+1) + novl
     else
       do k = 1, novl
         svec(nid) = i_(k)
         dvec(nid) = j_(k)
         sf = sf_(k)
         df = df_(k)
         ! REGRID NEXT DIMENSION
         CALL NREGRID(s, sax, dax, d, RG_TYPE, sovl, dovl, rcnt  &
                   ,m+1, nr, svec, dvec                               &
                   ,dio, sf, df)
       END DO
     endif

     deallocate(df_, sf_, i_,j_)

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
     DEALLOCATE(ms, md)
!     DEALLOCATE(ms, md, stat=status)
!     CALL ERRMSG(substr,status,32)
     DEALLOCATE(sb, db)
!     DEALLOCATE(sb, db, stat=status)
!     CALL ERRMSG(substr,status,33)

  ELSE  ! (C2): INVARIANT DIMENSION: ( m > nr )  .......................

     ! SAVE OVERLAP OF NEXT DIMENSION CALCULATED SO FAR ...
     zso = sovl(dio(m+1))
     zdo = dovl(dio(m+1))
     ! ... AND RESET VECTOR ELEMENT TO ZERO
     sovl(dio(m+1)) = REAL(0.0, DP)
     dovl(dio(m+1)) = REAL(0.0, DP)
     ! LOOP OVER ALL DIAGONAL MATRIX ELEMENTS ...
     ! ... AND REGRID NEXT DIMENSION
!!$     DO j=1, s(1)%dim(nid)
!!$        dvec(nid) = j
!!$        svec(nid) = j     ! INVARIANT !!!
!!$        ! REGRID NEXT DIMENSION
!!$        CALL NREGRID(s, sax, dax, d, RG_TYPE, sovl, dovl, rcnt  &
!!$             ,m+1, nr, svec, dvec                               &
!!$             ,dio, sf0, df0 )
!!$        ! OVERLAP IS 1 FOR INVARIANT DIMENSIONS
!!$        sovl(nid) = sovl(nid)+REAL(1.0, DP)/REAL(s(1)%dim(nid), DP)
!!$        dovl(nid) = dovl(nid)+REAL(1.0, DP)/REAL(s(1)%dim(nid), DP)
!!$     END DO

     IF (m == SIZE(dio)-1 ) then
       novl = s(1)%dim(nid)
       allocate(i_(novl),j_(novl),df_(novl),sf_(novl))
       DO j=1, novl
         i_(j) = j
         j_(j) = j
         sf_(j) = sf0
         df_(j) = df0
       enddo
       CALL NREGRID_last(s, sax, d, RG_TYPE, sovl, dovl        &
                  ,m+1, svec, dvec, dio, novl, i_, j_, sf_, df_ )
       rcnt(m+1) = rcnt(m+1) + novl
       deallocate(df_, sf_, i_,j_)
     else
       DO j=1, s(1)%dim(nid)
         dvec(nid) = j
         svec(nid) = j     ! INVARIANT !!!
         ! REGRID NEXT DIMENSION
         CALL NREGRID(s, sax, dax, d, RG_TYPE, sovl, dovl, rcnt  &
             ,m+1, nr, svec, dvec                               &
             ,dio, sf0, df0 )
         ! OVERLAP IS 1 FOR INVARIANT DIMENSIONS
!        sovl(nid) = sovl(nid)+REAL(1.0, DP)/REAL(s(1)%dim(nid), DP)
!        dovl(nid) = dovl(nid)+REAL(1.0, DP)/REAL(s(1)%dim(nid), DP)
       END DO
     endif
     sovl(nid) = sovl(nid)+1.0_dp
     dovl(nid) = dovl(nid)+1.0_dp

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
!  DEALLOCATE(dio, STAT=status)
!  CALL ERRMSG(substr,status,34)
!  DEALLOCATE(svec, dvec, STAT=status)
!  CALL ERRMSG(substr,status,35)
  IF ( .not. PRESENT(GDIO)) THEN
    DEALLOCATE(dio)
  endif
  IF ( .not. PRESENT(GSVEC)) THEN
    DEALLOCATE(svec)
  endif
  IF ( .not. PRESENT(GDVEC)) THEN
    DEALLOCATE(dvec)
  endif

  if ( m .eq. 1 ) deallocate (vec)

END SUBROUTINE NREGRID
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE NREGRID_last(s, sax, d, RG_TYPE, sovl, dovl        &
                  ,m, svec, dvec                               &
                  ,dio, novl, i_, j_, sf_, df_ )

  IMPLICIT NONE

  ! I/O
  TYPE (NARRAY), DIMENSION(:), INTENT(IN)    :: S       ! SOURCE DATA FIELDS
  TYPE (NARRAY), DIMENSION(:)                :: D       ! DEST. DATA FIELDS
  TYPE (AXIS)  , DIMENSION(:), INTENT(IN)    :: SAX     ! AXES (SOURCE)
  INTEGER      , DIMENSION(:), INTENT(IN)    :: RG_TYPE ! REGRIDDING TYPE
  REAL (DP),     DIMENSION(:)                :: SOVL  ! GLOB. OVERLAP FRACTION
  REAL (DP),     DIMENSION(:)                :: DOVL  ! GLOB. OVERLAP FRACTION
  integer, intent(in)                        :: novl
  real (dp), dimension(:), intent(in)        :: df_, sf_
  integer,   dimension(:), intent(in)        :: i_, j_

  INTEGER,                INTENT(IN)            :: M    ! CURRENT DIMENSION
! REAL (DP),              INTENT(IN)   :: SF0   ! CURRENT OVERLAP FRACTION
! REAL (DP),              INTENT(IN)   :: DF0   ! CURRENT OVERLAP FRACTION
  REAL (DP)                            :: SF0   ! CURRENT OVERLAP FRACTION
  REAL (DP)                            :: DF0   ! CURRENT OVERLAP FRACTION

  INTEGER, DIMENSION(:)          :: DIO  ! DIMENSION ORDER
  INTEGER, DIMENSION(:)          :: SVEC ! CNT. VECTOR (SOURCE)
  INTEGER, DIMENSION(:)          :: DVEC ! CNT. VECTOR (DEST.)

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: SUBSTR = 'NREGRID_last'
  INTEGER                                 :: I,J, K ! COUNTER
  INTEGER                                 :: NVAR   ! NUMBER OF FIELDS
  INTEGER                                 :: NID    ! COUNT NO. OF DIMENSIONS
  INTEGER                                 :: NDIM   ! NUMBER OF DIMENSIONS
  INTEGER                                 :: VO     ! OFFSET FOR BOUND VEC.
  INTEGER                                 :: VL     ! LENGTH OF BOUND VEC.
  INTEGER                                 :: VTYPE          ! VARIABLE TYPE

  INTEGER                                 :: NPOS,DACCHELP
  INTEGER                                 :: DDIM,SDIM

  integer                                 :: nid_1  ! dimension one level back

     ! INIT
     ! NUMBER OF DIMENSIONS
     NDIM = SIZE(SAX)
     NVAR = SIZE(S)


     ! SET CURRENT DIMENSION
     nid = dio(m)
     if ( m > 1 ) nid_1 = dio(m-1)


!    zso = REAL(0.0, DP)
!    zdo = REAL(0.0, DP)
     DO k=1, nvar    ! LOOP OVER VARIABLES

!       vtype = QTYPE_NARRAY(s(k))
        vtype = s(k)%type
       IF ( (SIZE(d(k)%dim).GT.4) .OR. (SIZE(s(k)%dim).GT.4) ) THEN

         do i = 1, novl
           if ( m > 1 ) then
             svec(nid_1) = i_(i)
             dvec(nid_1) = j_(i)
           endif
           sf0 = sf_(i)
           df0 = df_(i)

           SELECT CASE(RG_TYPE(k))
           CASE(RG_EXT)
             DO j=1, s(k)%dim(nid)
               dvec(nid) = j
               svec(nid) = j     ! INVARIANT !!!
               ! CALCULATE DEST. DATA POINT
               vo = POSITION(d(k)%dim,dvec)   ! POSITION IN DEST.
               vl = POSITION(s(k)%dim,svec)   ! POSITION IN SOURCE
               IF(vtype.EQ.VTYPE_REAL) THEN
                 d(k)%vr(vo) = d(k)%vr(vo) + s(k)%vr(vl) * REAL(sf0, SP)
               ELSE
                 d(k)%vd(vo) = d(k)%vd(vo) + s(k)%vd(vl) * sf0
               END IF
               ! OVERLAP IS 1 FOR INVARIANT DIMENSIONS
               ! INTERVAL LENGTH, NUMBER OF VARIABLES
!             zso = zso + REAL(1.0, DP)/ &
!                  (REAL(s(k)%dim(nid), DP)*REAL(nvar, DP))
!             zdo = zdo + REAL(1.0, DP)/ &
!                  (REAL(s(k)%dim(nid), DP)*REAL(nvar, DP))
             END DO
           ! RESET FOR NEXT/PREVIOUS RECURSION STEP
!          dvec(nid) = 1
!          svec(nid) = 1
           CASE(RG_INT, RG_IDX, RG_IXF)
           !CASE(RG_IDX, RG_IXF)
           ! NOTE: IDX VAR WAS CONVERTET TO INDEX-FRACTION
             DO j=1, s(k)%dim(nid)
               dvec(nid) = j
               svec(nid) = j     ! INVARIANT !!!
               ! CALCULATE DEST. DATA POINT
               vo = POSITION(d(k)%dim,dvec)   ! POSITION IN DEST.
               vl = POSITION(s(k)%dim,svec)   ! POSITION IN SOURCE
               IF(vtype.EQ.VTYPE_REAL) THEN
                 d(k)%vr(vo) = d(k)%vr(vo) + s(k)%vr(vl) * REAL(df0, SP)
               ELSE
                 d(k)%vd(vo) = d(k)%vd(vo) + s(k)%vd(vl) * df0
               END IF
               ! OVERLAP IS 1 FOR INVARIANT DIMENSIONS:
               ! INTERVAL LENGTH, NUMBER OF VARIABLES
!             zso = zso + REAL(1.0, DP)/ &
!                  (REAL(s(k)%dim(nid), DP)*REAL(nvar, DP))
!             zdo = zdo + REAL(1.0, DP)/ &
!                  (REAL(s(k)%dim(nid), DP)*REAL(nvar, DP))
            END DO
            ! RESET FOR NEXT/PREVIOUS RECURSION STEP
!          dvec(nid) = 1
!          svec(nid) = 1
           !CASE DEFAULT
           END SELECT
         enddo

       ELSE

         SELECT CASE(RG_TYPE(k))
         CASE(RG_EXT)
           ddim=SIZE(d(k)%dim)
           sdim=SIZE(s(k)%dim)
           DO j=1, s(k)%dim(nid)
             dvec(nid) = j
             svec(nid) = j     ! INVARIANT !!!
             do i = 1, novl
               if ( m > 1 ) then
                 svec(nid_1) = i_(i)
                 dvec(nid_1) = j_(i)
               endif
               sf0 = sf_(i)
               df0 = df_(i)
               ! CALCULATE DEST. DATA POINT
               npos = dvec(1)

               SELECT CASE(ddim)
               CASE(2) ! size(dim)=2
                 dacchelp = 1
                 dacchelp = dacchelp*d(k)%dim(1)
                 npos     = npos + dacchelp*(dvec(2)-1)
               CASE(3) ! size(dim)=3
                 dacchelp = 1
                 dacchelp = dacchelp*d(k)%dim(1)
                 npos     = npos + dacchelp*(dvec(2)-1)
                 dacchelp = dacchelp*d(k)%dim(2)
                 npos     = npos + dacchelp*(dvec(3)-1)
               CASE(4) ! size(dim)=4
                 dacchelp = 1
                 dacchelp = dacchelp*d(k)%dim(1)
                 npos     = npos + dacchelp*(dvec(2)-1)
                 dacchelp = dacchelp*d(k)%dim(2)
                 npos     = npos + dacchelp*(dvec(3)-1)
                 dacchelp = dacchelp*d(k)%dim(3)
                 npos     = npos + dacchelp*(dvec(4)-1)
               END SELECT
               vo = npos   ! POSITION IN DEST.

               npos = svec(1)
               SELECT CASE(sdim)
               CASE(2) ! size(dim)=2
                 dacchelp = 1
                 dacchelp = dacchelp*s(k)%dim(1)
                 npos     = npos + dacchelp*(svec(2)-1)
               CASE(3) ! size(dim)=3
                 dacchelp = 1
                 dacchelp = dacchelp*s(k)%dim(1)
                 npos     = npos + dacchelp*(svec(2)-1)
                 dacchelp = dacchelp*s(k)%dim(2)
                 npos     = npos + dacchelp*(svec(3)-1)
               CASE(4) ! size(dim)=4
                 dacchelp = 1
                 dacchelp = dacchelp*s(k)%dim(1)
                 npos     = npos + dacchelp*(svec(2)-1)
                 dacchelp = dacchelp*s(k)%dim(2)
                 npos     = npos + dacchelp*(svec(3)-1)
                 dacchelp = dacchelp*s(k)%dim(3)
                 npos     = npos + dacchelp*(svec(4)-1)
               END SELECT
               vl = npos   ! POSITION IN SOURCE

               IF(vtype.EQ.VTYPE_REAL) THEN
                 d(k)%vr(vo) = d(k)%vr(vo) + s(k)%vr(vl) * REAL(sf0, SP)
               ELSE
                 d(k)%vd(vo) = d(k)%vd(vo) + s(k)%vd(vl) * sf0
               END IF
               ! OVERLAP IS 1 FOR INVARIANT DIMENSIONS
               ! INTERVAL LENGTH, NUMBER OF VARIABLES
!             zso = zso + REAL(1.0, DP)/ &
!                  (REAL(s(k)%dim(nid), DP)*REAL(nvar, DP))
!             zdo = zdo + REAL(1.0, DP)/ &
!                  (REAL(s(k)%dim(nid), DP)*REAL(nvar, DP))
             enddo
           END DO
           ! RESET FOR NEXT/PREVIOUS RECURSION STEP
!          dvec(nid) = 1
!          svec(nid) = 1

           CASE(RG_INT, RG_IDX, RG_IXF)
           !CASE(RG_IDX, RG_IXF)
           ! NOTE: IDX VAR WAS CONVERTET TO INDEX-FRACTION
             ddim=SIZE(d(k)%dim)
             sdim=SIZE(s(k)%dim)
             DO j=1, s(k)%dim(nid)
               dvec(nid) = j
               svec(nid) = j     ! INVARIANT !!!
               do i = 1, novl
                 if ( m > 1 ) then
                   svec(nid_1) = i_(i)
                   dvec(nid_1) = j_(i)
                 endif
                 sf0 = sf_(i)
                 df0 = df_(i)
                ! CALCULATE DEST. DATA POINT

                npos = dvec(1)
                SELECT CASE(ddim)
                CASE(2) ! size(dim)=2
                  dacchelp = 1
                  dacchelp = dacchelp*d(k)%dim(1)
                  npos     = npos + dacchelp*(dvec(2)-1)
                CASE(3) ! size(dim)=3
                  dacchelp = 1
                  dacchelp = dacchelp*d(k)%dim(1)
                  npos     = npos + dacchelp*(dvec(2)-1)
                  dacchelp = dacchelp*d(k)%dim(2)
                  npos     = npos + dacchelp*(dvec(3)-1)
                CASE(4) ! size(dim)=4
                  dacchelp = 1
                  dacchelp = dacchelp*d(k)%dim(1)
                  npos     = npos + dacchelp*(dvec(2)-1)
                  dacchelp = dacchelp*d(k)%dim(2)
                  npos     = npos + dacchelp*(dvec(3)-1)
                  dacchelp = dacchelp*d(k)%dim(3)
                  npos     = npos + dacchelp*(dvec(4)-1)
                END SELECT
                vo = npos   ! POSITION IN DEST.

                npos = svec(1)
                SELECT CASE(sdim)
                CASE(2) ! size(dim)=2
                  dacchelp = 1
                  dacchelp = dacchelp*s(k)%dim(1)
                  npos     = npos + dacchelp*(svec(2)-1)
                CASE(3) ! size(dim)=3
                  dacchelp = 1
                  dacchelp = dacchelp*s(k)%dim(1)
                  npos     = npos + dacchelp*(svec(2)-1)
                  dacchelp = dacchelp*s(k)%dim(2)
                  npos     = npos + dacchelp*(svec(3)-1)
                CASE(4) ! size(dim)=4
                  dacchelp = 1
                  dacchelp = dacchelp*s(k)%dim(1)
                  npos     = npos + dacchelp*(svec(2)-1)
                  dacchelp = dacchelp*s(k)%dim(2)
                  npos     = npos + dacchelp*(svec(3)-1)
                  dacchelp = dacchelp*s(k)%dim(3)
                  npos     = npos + dacchelp*(svec(4)-1)
                END SELECT
                vl = npos   ! POSITION IN SOURCE

                IF(vtype.EQ.VTYPE_REAL) THEN
                  d(k)%vr(vo) = d(k)%vr(vo) + s(k)%vr(vl) * REAL(df0, SP)
                ELSE
                  d(k)%vd(vo) = d(k)%vd(vo) + s(k)%vd(vl) * df0
                END IF
                ! OVERLAP IS 1 FOR INVARIANT DIMENSIONS:
                ! INTERVAL LENGTH, NUMBER OF VARIABLES
!             zso = zso + REAL(1.0, DP)/ &
!                  (REAL(s(k)%dim(nid), DP)*REAL(nvar, DP))
!             zdo = zdo + REAL(1.0, DP)/ &
!                  (REAL(s(k)%dim(nid), DP)*REAL(nvar, DP))
              enddo
            END DO
            ! RESET FOR NEXT/PREVIOUS RECURSION STEP
!          dvec(nid) = 1
!          svec(nid) = 1
           !CASE DEFAULT
           END SELECT

         END IF

       END DO  ! LOOP OVER VARIABLES

       ! RESET FOR NEXT/PREVIOUS RECURSION STEP
       dvec(nid) = 1
       svec(nid) = 1

       ! SUM OVERLAP
!    sovl(nid) = sovl(nid) + zso
!    dovl(nid) = dovl(nid) + zdo
       sovl(nid) = sovl(nid) + real(novl,DP)
       dovl(nid) = dovl(nid) + real(novl,DP)

END SUBROUTINE NREGRID_last
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE NREGRID_STAT(sax, dax, sovl, dovl, nai, nao, rcnt)

  IMPLICIT NONE

  ! I/O
  TYPE (axis),   DIMENSION(:), INTENT(IN) :: sax, dax
  REAL (DP),     DIMENSION(:), INTENT(IN) :: sovl, dovl
  TYPE (narray), DIMENSION(:), INTENT(IN) :: nai, nao
  INTEGER,       DIMENSION(:), INTENT(IN) :: rcnt

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

! ------------------------------------------------------------------
SUBROUTINE ERRMSG(routine, status, pos)

  IMPLICIT NONE

  ! I/O
  CHARACTER(LEN=*), INTENT(IN)  :: routine
  INTEGER,          INTENT(IN)  :: status
  INTEGER,          INTENT(IN)  :: pos

  IF (status == 0) THEN
     RETURN
  ELSE
     CALL RGMSG(routine, RGMLE,  'ERROR STATUS ',status,' ', .false.)
     CALL RGMSG(routine, RGMLEC, 'AT POSITION ',pos,' !')
  END IF

END SUBROUTINE ERRMSG
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE RGMSG_C(routine, level, c, lstop)

#if ! defined(NOMPI)
  USE messy_ncregrid_mpi,  ONLY: ncregrid_abort
#endif

  IMPLICIT NONE

  ! I/O
  CHARACTER(LEN=*), INTENT(IN)            :: routine
  INTEGER,          INTENT(IN)            :: level
  CHARACTER(LEN=*), INTENT(IN)            :: c
  LOGICAL,          INTENT(IN), OPTIONAL  :: lstop

  ! LOCAL
  LOGICAL :: llstop
  ! mz_pj_20080801+
  INTEGER :: iou
  LOGICAL :: opened
  CHARACTER(LEN=20) :: errfname = ''
  ! mz_pj_20080801-

  ! INIT
  IF (PRESENT(lstop)) THEN
     llstop = lstop
  ELSE
     IF ((level == RGMLE).OR.(level == RGMLEC)) THEN
        llstop = .true.      ! STOP ON ERROR
     ELSE
        llstop = .false.
     END IF
  END IF

  ! mz_pj_20080801+
  IF ((level == RGMLE).OR.(level == RGMLEC)) THEN
     DO iou=100,300
        INQUIRE(unit=iou,opened=opened)
        IF (.NOT.opened) EXIT
     END DO
     IF (MY_RANK < 0) THEN
        errfname = 'ERROR.ncregrid'
     ELSE
        WRITE(errfname,'(a14,a1,i4.4)') 'ERROR.ncregrid','.',MY_RANK
     ENDIF
  END IF
  ! mz_pj_20080801-

  SELECT CASE(level)
  CASE(RGMLE)   ! ERROR MESSAGE
     IF (IAND(MSGMODE, MSGMODE_E) == MSGMODE_E) THEN
        WRITE(*,*) '*** ',TRIM(routine),' ERROR: '
        WRITE(*,*) '    ',TRIM(c)
        ! mz_pj_20080801+
        OPEN(iou, FILE=TRIM(errfname), STATUS='UNKNOWN')
        WRITE(iou,*) '*** ',TRIM(routine),' ERROR: '
        WRITE(iou,*) '    ',TRIM(c)
        CLOSE(iou)
        ! mz_pj_20080801-
     END IF
  CASE(RGMLEC) ! ERROR MESSAGE CONTINUED
     IF (IAND(MSGMODE, MSGMODE_E) == MSGMODE_E) THEN
        WRITE(*,*) '    ',TRIM(c)
        ! mz_pj_20080801+
        OPEN(iou, FILE=TRIM(errfname), STATUS='OLD',POSITION='APPEND')
        WRITE(iou,*) '    ',TRIM(c)
        CLOSE(iou)
        ! mz_pj_20080801-
     END IF
  CASE(RGMLVL)  ! LITTLE VERBOSE
     IF (IAND(MSGMODE, MSGMODE_VL) == MSGMODE_VL) THEN
        WRITE(*,*) TRIM(c)
     END IF
  CASE(RGMLVLC) ! LITTLE VERBOSE CONTINUED
     IF (IAND(MSGMODE, MSGMODE_VL) == MSGMODE_VL) THEN
        WRITE(*,*) '    ',TRIM(c)
     END IF
  CASE(RGMLW) ! WARNING MESSAGE
     IF (IAND(MSGMODE, MSGMODE_W) == MSGMODE_W) THEN
        WRITE(*,*) '+++ ',TRIM(routine),' WARNING: '
        WRITE(*,*) '    ',TRIM(c)
     END IF
  CASE(RGMLWC) ! WARNING MESSAGE CONTINUED
     IF (IAND(MSGMODE, MSGMODE_W) == MSGMODE_W) THEN
        WRITE(*,*) '    ',TRIM(c)
     END IF
  CASE(RGMLVM)  ! MEDIUM VERBOSE
     IF (IAND(MSGMODE, MSGMODE_VM) == MSGMODE_VM) THEN
        WRITE(*,*) TRIM(c)
     END IF
  CASE(RGMLVMC) ! MEDIUM VERBOSE CONTINUED
     IF (IAND(MSGMODE, MSGMODE_VM) == MSGMODE_VM) THEN
        WRITE(*,*) '    ',TRIM(c)
     END IF
  CASE(RGMLI)  ! INFO MESSAGE
     IF (IAND(MSGMODE, MSGMODE_I) == MSGMODE_I) THEN
        WRITE(*,*) '=== ',TRIM(routine),' INFO: '
        WRITE(*,*) '    ',TRIM(c)
     END IF
  CASE(RGMLIC) ! INFO MESSAGE CONTINUED
     IF (IAND(MSGMODE, MSGMODE_I) == MSGMODE_I) THEN
        WRITE(*,*) '    ',TRIM(c)
     END IF
  CASE DEFAULT
  END SELECT

#if ! defined(NOMPI)
  IF (llstop) CALL ncregrid_abort('ncregrid')
#else
  IF (llstop) STOP
#endif

END SUBROUTINE RGMSG_C
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE RGMSG_I(routine, level, c1, i, c2, lstop)

  IMPLICIT NONE

  ! I/O
  CHARACTER(LEN=*), INTENT(IN)            :: routine
  INTEGER,          INTENT(IN)            :: level
  CHARACTER(LEN=*), INTENT(IN)            :: c1
  INTEGER,          INTENT(IN)            :: i
  CHARACTER(LEN=*), INTENT(IN)            :: c2
  LOGICAL,          INTENT(IN), OPTIONAL  :: lstop

  ! LOCAL
  CHARACTER(LEN=1000) :: istr = ''

  WRITE(istr,*) TRIM(c1),i,TRIM(c2)
  CALL RGMSG_C(routine, level, TRIM(istr), lstop)

END SUBROUTINE RGMSG_I
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE RGMSG_IA(routine, level, c1, i, c2, lstop)

  IMPLICIT NONE

  ! I/O
  CHARACTER(LEN=*),      INTENT(IN)            :: routine
  INTEGER,               INTENT(IN)            :: level
  CHARACTER(LEN=*),      INTENT(IN)            :: c1
  INTEGER, DIMENSION(:), INTENT(IN)            :: i
  CHARACTER(LEN=*),      INTENT(IN)            :: c2
  LOGICAL,               INTENT(IN), OPTIONAL  :: lstop

  ! LOCAL
  CHARACTER(LEN=1000) :: istr = ''

  WRITE(istr,*) TRIM(c1),i,TRIM(c2)
  CALL RGMSG_C(routine, level, TRIM(istr), lstop)

END SUBROUTINE RGMSG_IA
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE RGMSG_R(routine, level, c1, r, c2, lstop)

  IMPLICIT NONE

  ! I/O
  CHARACTER(LEN=*), INTENT(IN)            :: routine
  INTEGER,          INTENT(IN)            :: level
  CHARACTER(LEN=*), INTENT(IN)            :: c1
  REAL,             INTENT(IN)            :: r
  CHARACTER(LEN=*), INTENT(IN)            :: c2
  LOGICAL,          INTENT(IN), OPTIONAL  :: lstop

  ! LOCAL
  CHARACTER(LEN=1000) :: rstr = ''

  WRITE(rstr,*) TRIM(c1),r,TRIM(c2)
  CALL RGMSG_C(routine, level, TRIM(rstr), lstop)

END SUBROUTINE RGMSG_R
! ------------------------------------------------------------------

! op_mr_20170321+
! ------------------------------------------------------------------
SUBROUTINE RGMSG_D(routine, level, c1, d, c2, lstop)

  IMPLICIT NONE

  ! I/O
  CHARACTER(LEN=*), INTENT(IN)            :: routine
  INTEGER,          INTENT(IN)            :: level
  CHARACTER(LEN=*), INTENT(IN)            :: c1
  REAL(dp),         INTENT(IN)            :: d
  CHARACTER(LEN=*), INTENT(IN)            :: c2
  LOGICAL,          INTENT(IN), OPTIONAL  :: lstop

  ! LOCAL
  CHARACTER(LEN=1000) :: rstr = ''

  WRITE(rstr,*) TRIM(c1),d,TRIM(c2)
  CALL RGMSG_C(routine, level, TRIM(rstr), lstop)

END SUBROUTINE RGMSG_D
! ------------------------------------------------------------------
! op_mr_20170321-

! ******************************************************************
END MODULE MESSY_NCREGRID_BASE
! ******************************************************************

#endif
