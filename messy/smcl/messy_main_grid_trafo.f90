! ******************************************************************
! ------------------------------------------------------------------
MODULE MESSY_MAIN_GRID_TRAFO
! ------------------------------------------------------------------
! Author: Patrick Joeckel, MPICH, Mainz, June 2002
!         Bastian Kern,    MPICH, Mainz, July 2009
!            extended for curvilinear ocean grids by
!         Astrid  Kerkweg, UniMz, Mainz, 2012-2013
!         -> re-structured/expanded for more general GRID application
! ******************************************************************

  ! MESSY/SMCL
  USE MESSY_MAIN_TOOLS,         ONLY: STR
  
  USE MESSY_MAIN_CONSTANTS_MEM, ONLY: dp, sp, i4, i8
  ! - GRID
  USE MESSY_MAIN_GRID,          ONLY: t_geohybgrid, INIT_GEOHYBGRID         &
                                    , COPY_GEOHYBGRID, RGEMPTY              &
                                    , PRINT_GEOHYBGRID, RGMAX
  USE MESSY_MAIN_GRID_NETCDF,   ONLY: t_narray,  t_ncvar, GRD_MAXSTRLEN     & 
                                    , NULL_DIMID, NULL_VARID, ERRMSG        &
                                    , RGMLE, RGMLEC, RGMLI, RGMLIC, RGMLW   &
                                    , RGMLVL, MSGMODE, MSGMODE_VM           &
                                    , VTYPE_DOUBLE, VTYPE_INT,  VTYPE_REAL  &
                                    , VTYPE_BYTE,   VTYPE_CHAR, VTYPE_UNDEF &
                                    ! FUNCTIONS
                                    , QDEF_NCVAR,  QCMP_NCDIM               &
                                    , POSITION,    ELEMENT                  & 
                                    ! SUBROUTINES
                                    , INIT_NCDIM,  COPY_NCDIM               &
                                    , INIT_NCVAR,  COPY_NCVAR               &
                                    , INIT_NCATT,  COPY_NCATT,  ADD_NCATT   &
                                    , INIT_NARRAY, COPY_NARRAY, SORT_NARRAY &
                                    , REORDER_NARRAY, DOUBLE_NARRAY, RGMSG  &
                                    , NARRAYMINVAL, NARRAYMAXVAL            &
                                    , PRINT_NCVAR
  USE netcdf                                    

  IMPLICIT NONE
  PUBLIC

  INTRINSIC :: TRIM, ASSOCIATED, ABS, TINY, PRESENT, PRODUCT   &
             , SIZE, MAXVAL, REAL, INT, MINVAL

  PRIVATE   :: TRIM, ASSOCIATED, ABS, TINY, PRESENT, PRODUCT   &
             , SIZE, MAXVAL, REAL, INT, MINVAL


  CHARACTER(*), PARAMETER, PUBLIC :: GRIDTRAFOVERS = '1.0'

  ! grid_trafo +
  ! INTERPOLATION TYPES
  INTEGER, PARAMETER, PUBLIC :: GTRF_NONE = 1
  INTEGER, PARAMETER, PUBLIC :: GTRF_SCRP = 2
  INTEGER, PARAMETER, PUBLIC :: GTRF_NRGD = 3
  INTEGER, PARAMETER, PUBLIC :: GTRF_MAX  = 3

  CHARACTER(LEN=4), DIMENSION(GTRF_MAX), PUBLIC :: GTRFCHAR = &
       (/'NONE', 'SCRP', 'NRGD'/)


  ! REGRIDDING VARIABLE TYPES
  INTEGER, PARAMETER, PUBLIC :: RG_INT = 1
  INTEGER, PARAMETER, PUBLIC :: RG_EXT = 2
  INTEGER, PARAMETER, PUBLIC :: RG_IDX = 3
  INTEGER, PARAMETER, PUBLIC :: RG_IXF = 4

  CHARACTER(LEN=3), DIMENSION(4), PUBLIC :: RGTstr = (/'INT','EXT','IDX','IXF'/)

  ! grid_trafo -
  ! INTERNAL ORDER OF LON, LAT, LEV
  ! NOTE: THIS HAS TO BE CONSISTENT WITH SEARCH ORDER IN
  !       SUBROUTINE GEOHYBGRID_AXES
  INTEGER,   PARAMETER, PUBLIC :: GORD_LON = 1
  INTEGER,   PARAMETER, PUBLIC :: GORD_LAT = 2
  INTEGER,   PARAMETER, PUBLIC :: GORD_LEV = 3

  INTEGER, PARAMETER          :: I_UNDEF     = -99
  INTEGER, PARAMETER          :: I_0_360     = 1 ![0,360]
  INTEGER, PARAMETER          :: I_m180_180  = 2 ![-180:180]
  INTEGER, PARAMETER          :: I_0_180     = 3 ![could be both]
  ! Interval in which the destination local grid is defined
  INTEGER, SAVE               :: ivaldefd = I_UNDEF
  ! Interval in which the source local grid is defined
  INTEGER, SAVE               :: ivaldefs = I_UNDEF
  LOGICAL, SAVE               :: lunstruct = .FALSE.

  REAL(dp), DIMENSION(2),SAVE :: minmaxdest

!! NOTE: DOES NOT WORK PROPERLY FOR SOME COMPILERS ...
!  INTERFACE ASSIGNMENT (=)
!     MODULE PROCEDURE COPY_GEOHYBGRID
!  END INTERFACE

  ! PUBLIC :: SWITCH_GEOHYBGRID
  ! PUBLIC :: CHECK_GEOHYBGRID
  ! PUBLIC :: SORT_GEOHYBGRID
  ! PUBLIC :: H2PSIG
  ! PUBLIC :: COMPLETE_GEOHYBGRID
  ! PUBLIC :: CHECK_NCVAR_ON_GEOHYBGRID
  ! PUBLIC :: SORT_GEOHYBGRID_NCVAR
  ! PUBLIC :: PACK_GEOHYBGRID_NCVAR
  ! PUBLIC :: BALANCE_GEOHYBGRID
  ! PUBLIC :: BALANCE_GEOHYBGRID_NCVAR
  ! PUBLIC :: BALANCE_GEOHYBGRID_TIME
  ! PUBLIC :: REDUCE_INPUT_GRID
  ! PUBLIC :: PHI2PHIROT
  ! PUBLIC :: PHIROT2PHI
  ! PUBLIC :: RLA2RLAROT
  ! PUBLIC :: RLAROT2RLA
  ! PUBLIC :: UV2UVROT_VEC
  ! PUBLIC :: UVROT2UV_VEC
  ! PUBLIC :: UVROT2UV_VEC_JLOOP
  ! PUBLIC :: EXPAND_INPUT_GRID_LON
  ! PUBLIC :: H2PSIG_3D
  ! PUBLIC :: RNGADJ_NARRAY
  ! PUBLIC :: RNGADJ_ARRAY
  ! PUBLIC :: IMMI_NARRAY
  ! PUBLIC :: IMMI_NCVAR
  
CONTAINS

! ------------------------------------------------------------------
SUBROUTINE SWITCH_GEOHYBGRID(g, lx, ly, lz)

  ! this subroutine initializes (=deletes) dimensions of a geohybgrid,
  ! if respective logicals are true:
  ! - lx = .FALSE: : re-initialise lonm and loni
  ! - ly = .FALSE: : re-initialise latm and lati
  ! - lz = .FALSE: : re-initialise vertical definition

  IMPLICIT NONE

  ! I/O
  TYPE(t_geohybgrid), INTENT(INOUT) :: g
  LOGICAL           , INTENT(IN)    :: lx, ly, lz

  IF (.NOT.lx) THEN
     CALL INIT_NCVAR(g%loni)
     CALL INIT_NCVAR(g%lonm)
  END IF

  IF (.NOT.ly) THEN
     CALL INIT_NCVAR(g%lati)
     CALL INIT_NCVAR(g%latm)
  END IF
  
  IF (.NOT. lx .OR. .NOT. ly) THEN
     CALL INIT_NCVAR(g%clonm)
     CALL INIT_NCVAR(g%cloni)
     CALL INIT_NCVAR(g%clati)
     CALL INIT_NCVAR(g%clatm)
     !
     CALL INIT_NCVAR(g%rlonm)
     CALL INIT_NCVAR(g%rloni)
     CALL INIT_NCVAR(g%rlati)
     CALL INIT_NCVAR(g%rlatm)
     !
     CALL INIT_NCVAR(g%ulonm)
     CALL INIT_NCVAR(g%uloni)
     CALL INIT_NCVAR(g%ulati)
     CALL INIT_NCVAR(g%ulatm)
     !
  ENDIF

  IF (.NOT.lz) THEN
     CALL INIT_NCVAR(g%hyai)
     CALL INIT_NCVAR(g%hyam)
     CALL INIT_NCVAR(g%hybi)
     CALL INIT_NCVAR(g%hybm)
     CALL INIT_NCVAR(g%p0)
     CALL INIT_NCVAR(g%ps)
     !
     CALL INIT_NCVAR(g%pressi)
     CALL INIT_NCVAR(g%pressm)
     !
  END IF

END SUBROUTINE SWITCH_GEOHYBGRID
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE CHECK_GEOHYBGRID(grid)

  IMPLICIT NONE

  ! I/O
  TYPE (t_geohybgrid), INTENT(INOUT) :: grid

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'CHECK_GEOHYBGRID'

  ! LONGITUDE
  IF (grid%lonm%ndims > 1) THEN
     CALL RGMSG(substr, RGMLE, &
          'LONGITUDE VARIABLE '''//TRIM(grid%lonm%name)//'''', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IN FILE '''//TRIM(grid%file)//'''', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IS ',grid%lonm%ndims,'-DIMENSIONAL !', .false.)
  END IF

  IF (grid%loni%ndims > 1) THEN
     CALL RGMSG(substr, RGMLE, &
          'LONGITUDE INTERFACES VARIABLE '''//TRIM(grid%loni%name)//'''',&
          .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IN FILE '''//TRIM(grid%file)//'''', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IS ',grid%loni%ndims,'-DIMENSIONAL !', .false.)
  END IF

  ! LATITUDE
  IF (grid%latm%ndims > 1) THEN
     CALL RGMSG(substr, RGMLE, &
          'LATITUDE VARIABLE '''//TRIM(grid%latm%name)//'''', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IN FILE '''//TRIM(grid%file)//'''', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IS ',grid%latm%ndims,'-DIMENSIONAL !', .false.)
  END IF

  IF (grid%lati%ndims > 1) THEN
     CALL RGMSG(substr, RGMLE, &
          'LATITUDE INTERFACES VARIABLE '''//TRIM(grid%lati%name)//'''',&
          .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IN FILE '''//TRIM(grid%file)//'''', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IS ',grid%lati%ndims,'-DIMENSIONAL !', .false.)
  END IF

  ! TIME
  IF (grid%timem%ndims > 1) THEN
     CALL RGMSG(substr, RGMLE, &
          'TIME VARIABLE '''//TRIM(grid%timem%name)//'''', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IN FILE '''//TRIM(grid%file)//'''', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IS ',grid%timem%ndims,'-DIMENSIONAL !')
  END IF

  IF (grid%timei%ndims > 1) THEN
     CALL RGMSG(substr, RGMLE, &
          'TIME INTERFACES VARIABLE '''//TRIM(grid%timei%name)//'''', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IN FILE '''//TRIM(grid%file)//'''', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IS ',grid%timei%ndims,'-DIMENSIONAL !')
  END IF

  ! HYA
  IF (grid%hyam%ndims > 1) THEN
     CALL RGMSG(substr, RGMLE, &
          'HYAM VARIABLE '''//TRIM(grid%hyam%name)//'''', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IN FILE '''//TRIM(grid%file)//'''', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IS ',grid%hyam%ndims,'-DIMENSIONAL !')
  END IF

  IF (grid%hyai%ndims > 1) THEN
     CALL RGMSG(substr, RGMLE, &
          'HYAI VARIABLE '''//TRIM(grid%hyai%name)//'''', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IN FILE '''//TRIM(grid%file)//'''', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IS ',grid%hyai%ndims,'-DIMENSIONAL !')
  END IF

  ! HYB
  IF (grid%hybm%ndims > 1) THEN
     CALL RGMSG(substr, RGMLE, &
          'HYBM VARIABLE '''//TRIM(grid%hybm%name)//'''', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IN FILE '''//TRIM(grid%file)//'''', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IS ',grid%hybm%ndims,'-DIMENSIONAL !')
  END IF

  IF (grid%hybi%ndims > 1) THEN
     CALL RGMSG(substr, RGMLE, &
          'HYBI VARIABLE '''//TRIM(grid%hybi%name)//'''', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IN FILE '''//TRIM(grid%file)//'''', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IS ',grid%hybi%ndims,'-DIMENSIONAL !')
  END IF

  ! CHECK COMPATIBILITY OF HYA AND HYB
  IF (QDEF_NCVAR(grid%hyai).AND.QDEF_NCVAR(grid%hybi)) THEN
     IF (grid%hyai%dim(1)%len /= grid%hybi%dim(1)%len) THEN
        CALL RGMSG(substr, RGMLE, &
             'DIMENSIONS OF HYAI AND HYBI', .false.)
        CALL RGMSG(substr, RGMLEC, &
             'IN FILE '''//TRIM(grid%file)//'''', .false.)
        CALL RGMSG(substr, RGMLEC, &
             'HAVE DIFFERENT LENGTH: ', .false.)
        CALL RGMSG(substr, RGMLEC, &
             'HYAI = '''//TRIM(grid%hyai%name)//''' :', &
             grid%hyai%dim(1)%len, ' ', .false.)
        CALL RGMSG(substr, RGMLEC, &
             'HYBI = '''//TRIM(grid%hybi%name)//''' :', &
             grid%hybi%dim(1)%len, ' ')
     END IF
  END IF

  IF (QDEF_NCVAR(grid%hyam).AND.QDEF_NCVAR(grid%hybm)) THEN
     IF (grid%hyam%dim(1)%len /= grid%hybm%dim(1)%len) THEN
        CALL RGMSG(substr, RGMLE, &
             'DIMENSIONS OF HYAM AND HYBM', .false.)
        CALL RGMSG(substr, RGMLEC, &
             'IN FILE '''//TRIM(grid%file)//'''', .false.)
        CALL RGMSG(substr, RGMLEC, &
             'HAVE DIFFERENT LENGTH: ', .false.)
        CALL RGMSG(substr, RGMLEC, &
             'HYAM = '''//TRIM(grid%hyam%name)//''' :', &
             grid%hyam%dim(1)%len, ' ', .false.)
        CALL RGMSG(substr, RGMLEC, &
             'HYBM = '''//TRIM(grid%hybm%name)//''' :', &
             grid%hybm%dim(1)%len, ' ')
     END IF
  END IF

  ! CHECK CONSISTENCY OF INTERFACES AND MIDs
  IF (QDEF_NCVAR(grid%hyai).AND.QDEF_NCVAR(grid%hyam)) THEN
     IF (grid%hyai%dim(1)%len /= (grid%hyam%dim(1)%len+1)) THEN
        CALL RGMSG(substr, RGMLE, &
             'DIMENSION LENGHTS OF HYAI AND HYAM', .false.)
        CALL RGMSG(substr, RGMLEC, &
             'IN FILE '''//TRIM(grid%file)//'''', .false.)
        CALL RGMSG(substr, RGMLEC, &
             'ARE NOT CONSISTENT: ', .false.)
        CALL RGMSG(substr, RGMLEC, &
             'HYAI = '''//TRIM(grid%hyai%name)//''' :', &
             grid%hyai%dim(1)%len, ' ', .false.)
        CALL RGMSG(substr, RGMLEC, &
             'HYAM = '''//TRIM(grid%hyam%name)//''' :', &
             grid%hyam%dim(1)%len, ' ')
     END IF
  END IF

  IF (QDEF_NCVAR(grid%hybi).AND.QDEF_NCVAR(grid%hybm)) THEN
     IF (grid%hybi%dim(1)%len /= (grid%hybm%dim(1)%len+1)) THEN
        CALL RGMSG(substr, RGMLE, &
             'DIMENSION LENGHTS OF HYBI AND HYBM', .false.)
        CALL RGMSG(substr, RGMLEC, &
             'IN FILE '''//TRIM(grid%file)//'''', .false.)
        CALL RGMSG(substr, RGMLEC, &
             'ARE NOT CONSISTENT: ', .false.)
        CALL RGMSG(substr, RGMLEC, &
             'HYBI = '''//TRIM(grid%hybi%name)//''' :', &
             grid%hybi%dim(1)%len, ' ', .false.)
        CALL RGMSG(substr, RGMLEC, &
             'HYBM = '''//TRIM(grid%hybm%name)//''' :', &
             grid%hybm%dim(1)%len, ' ')
     END IF
  END IF

  ! SURFACE PRESSURE (MUST BE ON LON, LAT, TIME)
  IF (QDEF_NCVAR(grid%ps)) THEN
     IF (grid%ps%ndims > 3) THEN
        CALL RGMSG(substr, RGMLE, &
             'SURFACE PRESSURE VARIABLE '''//TRIM(grid%ps%name)//'''', &
             .false.)
        CALL RGMSG(substr, RGMLEC, &
             'IN FILE '''//TRIM(grid%file)//'''', .false.)
        CALL RGMSG(substr, RGMLEC, &
             'IS ',grid%ps%ndims,'-DIMENSIONAL !')
     END IF
  END IF

  ! CHECK DIM's OF LENGTH 1
  IF (grid%lonm%ndims == 1) THEN
     IF (grid%lonm%dim(1)%len == 1) THEN
        IF (grid%loni%ndims /= 1) THEN
           IF ((ABS(grid%ranges(1,1) - RGEMPTY) <= TINY(RGEMPTY)).OR. &
                (ABS(grid%ranges(1,2) - RGEMPTY) <= TINY(RGEMPTY))) THEN
              CALL RGMSG(substr, RGMLE, &
                   'LENGTH OF LONGITUDE DIMENSION '''&
                   &//TRIM(grid%lonm%dim(1)%name)//''' IS 1 !',.false.)
              CALL RGMSG(substr, RGMLEC, &
                   'PLEASE SPECIFY ''?_lonr'' IN NAMELIST !')
           END IF
        END IF
     END IF
  END IF

  IF (grid%latm%ndims == 1) THEN
     IF (grid%latm%dim(1)%len == 1) THEN
        IF (grid%lati%ndims /= 1) THEN
           IF ((ABS(grid%ranges(2,1) - RGEMPTY) <= TINY(RGEMPTY)).OR. &
                (ABS(grid%ranges(2,2) - RGEMPTY) <= TINY(RGEMPTY))) THEN
              CALL RGMSG(substr, RGMLE, &
                   'LENGTH OF LATITUDE DIMENSION '''&
                   &//TRIM(grid%latm%dim(1)%name)//''' IS 1 !',.false.)
              CALL RGMSG(substr, RGMLEC, &
                   'PLEASE SPECIFY ''?_latr'' IN NAMELIST !')
           END IF
        END IF
     END IF
  END IF

  IF (grid%hyam%ndims == 1) THEN
     IF (grid%hyam%dim(1)%len == 1) THEN
        IF (grid%hyai%ndims /= 1) THEN
           IF ((ABS(grid%ranges(3,1) - RGEMPTY) <= TINY(RGEMPTY)).OR. &
                (ABS(grid%ranges(3,2) - RGEMPTY) <= TINY(RGEMPTY))) THEN
              CALL RGMSG(substr, RGMLE, &
                   'LENGTH OF LEVEL DIMENSION '''&
                   &//TRIM(grid%hyam%dim(1)%name)//''' IS 1 !',.false.)
              CALL RGMSG(substr, RGMLEC, &
                   'PLEASE SPECIFY ''?_hyar'' IN NAMELIST !')
           END IF
        END IF
     END IF
  END IF

  IF (grid%hybm%ndims == 1) THEN
     IF (grid%hybm%dim(1)%len == 1) THEN
        IF (grid%hybi%ndims /= 1) THEN
           IF ((ABS(grid%ranges(4,1) - RGEMPTY) <= TINY(RGEMPTY)).OR. &
                (ABS(grid%ranges(4,2) - RGEMPTY) <= TINY(RGEMPTY))) THEN
              CALL RGMSG(substr, RGMLE, &
                   'LENGTH OF LEVEL DIMENSION '''&
                   &//TRIM(grid%hybm%dim(1)%name)//''' IS 1 !',.false.)
              CALL RGMSG(substr, RGMLEC, &
                   'PLEASE SPECIFY ''?_hybr'' IN NAMELIST !')
           END IF
        END IF
     END IF
  END IF

  ! CONSISTENCY OF HYAI, HYBI, PS, P0 IS CHECKED IN
  ! -> GEOHYBGRID_AXES -> H2PSIG
  ! HYBRID LEVELS: HYAI, P0, HYBI, PS
  ! SIGMA LEVELS :           HYBI, (PS)
  ! CONST. PRESS.: HYAI, (P0)
  ! NOTE: PARAMETERS IN '()' ARE ONLY REQUIRED IF REGRIDDING IS
  !       PERFORMED IN PRESSURE COORDINATES

END SUBROUTINE CHECK_GEOHYBGRID
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE SORT_GEOHYBGRID(gi, go, gx, reverse)

  IMPLICIT NONE

  ! I/O
  TYPE (t_geohybgrid), INTENT(IN)    :: gi  ! INPUT GRID
  TYPE (t_geohybgrid), INTENT(OUT)   :: go  ! OUTPUT GRID
  TYPE (t_geohybgrid), INTENT(INOUT) :: gx  ! INDEX 'GRID'
  LOGICAL, OPTIONAL,   INTENT(IN)    :: reverse

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER        :: substr = 'SORT_GEOHYBGRID'
  LOGICAL                            :: lrev
  INTEGER                            :: i
  INTEGER                            :: olat, olon  ! order of dimensions in PS
  INTEGER, DIMENSION(:), ALLOCATABLE :: dil  ! dimension length vector of PS
  INTEGER, DIMENSION(:), POINTER     :: vec  ! element vector of PS
  INTEGER, DIMENSION(:), ALLOCATABLE :: svec ! element vector of sorted PS
  INTEGER                            :: vtype
  INTEGER                            :: status
  TYPE (t_narray)                    :: ts   ! temporal n-array for sigma level
  TYPE (t_narray)                    :: col  ! sigma level column
  TYPE (t_narray)                    :: cx   ! sigma level column sort indices

  IF (PRESENT(reverse)) THEN
     lrev = reverse
  ELSE
     lrev = .false.  ! DEFAULT
  END IF
  NULLIFY(vec)

  IF (lrev) THEN
     CALL COPY_GEOHYBGRID(go, gi)  ! INITIALIZE
     ! DO NOT OVERWRITE gx HERE !!!
     CALL RGMSG(substr, RGMLI, 'UN-SORTING GRID ...')
  ELSE
     CALL COPY_GEOHYBGRID(go, gi)
     CALL COPY_GEOHYBGRID(gx, gi)  ! INITIALIZE
     CALL RGMSG(substr, RGMLI, 'SORTING GRID ...')
  END IF

  IF (QDEF_NCVAR(go%lonm)) THEN
     CALL RGMSG(substr, RGMLIC, ' ... LONM : '''//TRIM(go%lonm%name)//'''')
     IF (.NOT.lrev) THEN
        CALL INIT_NARRAY(gx%lonm%dat)
        CALL IDX_NCVAR(gx%lonm)
     END IF
     CALL SORT_NARRAY(go%lonm%dat, gx%lonm%dat, lrev)
  END IF

  IF (QDEF_NCVAR(go%loni)) THEN
     CALL RGMSG(substr, RGMLIC, ' ... LONI : '''//TRIM(go%loni%name)//'''')
     IF (.NOT.lrev) THEN
        CALL INIT_NARRAY(gx%loni%dat)
        CALL IDX_NCVAR(gx%loni)
     END IF
     CALL SORT_NARRAY(go%loni%dat, gx%loni%dat, lrev)
  END IF

  IF (QDEF_NCVAR(go%latm)) THEN
     CALL RGMSG(substr, RGMLIC, ' ... LATM : '''//TRIM(go%latm%name)//'''')
     IF (.NOT.lrev) THEN
        CALL INIT_NARRAY(gx%latm%dat)
        CALL IDX_NCVAR(gx%latm)
     END IF
     CALL SORT_NARRAY(go%latm%dat, gx%latm%dat, lrev)
  END IF

  IF (QDEF_NCVAR(go%lati)) THEN
     CALL RGMSG(substr, RGMLIC, ' ... LATI : '''//TRIM(go%lati%name)//'''')
     IF (.NOT.lrev) THEN
        CALL INIT_NARRAY(gx%lati%dat)
        CALL IDX_NCVAR(gx%lati)
     END IF
     CALL SORT_NARRAY(go%lati%dat, gx%lati%dat, lrev)
  END IF

  ! SORT TIME HERE, OR DELETE IT IN gx
  IF (QDEF_NCVAR(go%timem)) THEN
     CALL RGMSG(substr, RGMLIC, ' ... TIMEM: '''//TRIM(go%timem%name)//'''')
     IF (.NOT.lrev) THEN
        CALL INIT_NARRAY(gx%timem%dat)
        CALL IDX_NCVAR(gx%timem)
     END IF
     CALL SORT_NARRAY(go%timem%dat, gx%timem%dat, lrev)
  END IF

  IF (QDEF_NCVAR(go%timei)) THEN
     CALL RGMSG(substr, RGMLIC, ' ... TIMEI: '''//TRIM(go%timei%name)//'''')
     IF (.NOT.lrev) THEN
        CALL INIT_NARRAY(gx%timei%dat)
        CALL IDX_NCVAR(gx%timei)
     END IF
     CALL SORT_NARRAY(go%timei%dat, gx%timei%dat, lrev)
  END IF

  ! P0
  IF (QDEF_NCVAR(go%p0)) THEN
     CALL RGMSG(substr, RGMLIC, ' ... P0   : '''//TRIM(go%p0%name)//'''')
     IF (.NOT.lrev) THEN
        CALL INIT_NARRAY(gx%p0%dat)
        CALL IDX_NCVAR(gx%p0)
     END IF
     CALL SORT_NARRAY(go%p0%dat, gx%p0%dat, lrev)
  END IF


  IF (QDEF_NCVAR(go%hyai).OR.QDEF_NCVAR(go%hybi)) THEN
     CALL RGMSG(substr, RGMLIC, ' ... HYBRID COEFFICIENTS (INTERFACES) :')
     CALL RGMSG(substr, RGMLIC, '     HYAI : '''//TRIM(go%hyai%name)//'''')
     CALL RGMSG(substr, RGMLIC, '     HYBI : '''//TRIM(go%hybi%name)//'''')
     IF (.NOT.lrev) THEN
        ! PRECALCULATE SIGMA LEVEL FIELD FOR ONE COLUMN
        CALL INIT_NARRAY(ts)
        CALL H2PSIG(ts,go%hyai%dat,go%hybi%dat,go%ps%dat,go%p0%dat,.false.)
        ! PICK OUT ONE COLUMN
        ! NOTE: FIRST DIMENSION IS ALONG HYBRID-COEFF. !
        !       (SEE H2PSIG)
        CALL INIT_NARRAY(col, 1, (/ ts%dim(1) /), VTYPE_DOUBLE)
        ALLOCATE(vec(SIZE(ts%dim)), STAT=status)
        CALL ERRMSG(substr,status,3)
        vec(:) = 1
        DO i=1, ts%dim(1)
           vec(1) = i
           col%vd(i) = ts%vd(POSITION(ts%dim,vec))
        END DO
        ! SORT THIS COLUMN
        CALL SORT_NARRAY(col, cx)
        ! SAVE SORT INDICES AND RE-ORDER DATA
        IF (QDEF_NCVAR(go%hyai)) THEN
           CALL INIT_NARRAY(gx%hyai%dat)
           CALL IDX_NCVAR(gx%hyai)
           CALL COPY_NARRAY(gx%hyai%dat, cx)
           CALL REORDER_NARRAY(go%hyai%dat, cx)
        END IF
        IF (QDEF_NCVAR(go%hybi)) THEN
           CALL INIT_NARRAY(gx%hybi%dat)
           CALL IDX_NCVAR(gx%hybi)
           CALL COPY_NARRAY(gx%hybi%dat, cx)
           CALL REORDER_NARRAY(go%hybi%dat, cx)
        END IF
        CALL INIT_NARRAY(ts)
        CALL INIT_NARRAY(col)
        CALL INIT_NARRAY(cx)
        DEALLOCATE(vec, STAT=status)
        CALL ERRMSG(substr,status,4)
        NULLIFY(vec)
     ELSE
        IF (QDEF_NCVAR(go%hyai)) THEN
           CALL SORT_NARRAY(go%hyai%dat, gx%hyai%dat, lrev)
        END IF
        IF (QDEF_NCVAR(go%hybi)) THEN
           CALL SORT_NARRAY(go%hybi%dat, gx%hybi%dat, lrev)
        END IF
     END IF
  END IF

  IF (QDEF_NCVAR(go%hyam).OR.QDEF_NCVAR(go%hybm)) THEN
     CALL RGMSG(substr, RGMLIC, ' ... HYBRID COEFFICIENTS (MID) :')
     CALL RGMSG(substr, RGMLIC, '     HYAM : '''//TRIM(go%hyam%name)//'''')
     CALL RGMSG(substr, RGMLIC, '     HYBM : '''//TRIM(go%hybm%name)//'''')
     !
     IF (.NOT.lrev) THEN
        ! PRECALCULATE SIGMA LEVEL FIELD FOR ONE COLUMN
        CALL INIT_NARRAY(ts)
        CALL H2PSIG(ts,go%hyam%dat,go%hybm%dat,go%ps%dat,go%p0%dat,.false.)

        ! PICK OUT ONE COLUMN
        ! NOTE: FIRST DIMENSION IS ALONG HYBRID-COEFF. !
        !       (SEE H2PSIG)
        CALL INIT_NARRAY(col, 1, (/ ts%dim(1) /), VTYPE_DOUBLE)
        ALLOCATE(vec(SIZE(ts%dim)), STAT=status)
        CALL ERRMSG(substr,status,1)
        vec(:) = 1
        DO i=1, ts%dim(1)
           vec(1) = i
           col%vd(i) = ts%vd(POSITION(ts%dim,vec))
        END DO
        ! SORT THIS COLUMN
        CALL SORT_NARRAY(col, cx)
        ! SAVE SORT INDICES AND RE-ORDER DATA
        IF (QDEF_NCVAR(go%hyam)) THEN
           CALL INIT_NARRAY(gx%hyam%dat)
           CALL IDX_NCVAR(gx%hyam)
           CALL COPY_NARRAY(gx%hyam%dat, cx)
           CALL REORDER_NARRAY(go%hyam%dat, cx)
        END IF
        IF (QDEF_NCVAR(go%hybm)) THEN
           CALL INIT_NARRAY(gx%hybm%dat)
           CALL IDX_NCVAR(gx%hybm)
           CALL COPY_NARRAY(gx%hybm%dat, cx)
           CALL REORDER_NARRAY(go%hybm%dat, cx)
        END IF
        CALL INIT_NARRAY(ts)
        CALL INIT_NARRAY(col)
        CALL INIT_NARRAY(cx)
        DEALLOCATE(vec, STAT=status)
        CALL ERRMSG(substr,status,2)
        NULLIFY(vec)
     ELSE
        IF (QDEF_NCVAR(go%hyam)) THEN
          CALL SORT_NARRAY(go%hyam%dat, gx%hyam%dat, lrev)
        END IF
        IF (QDEF_NCVAR(go%hybm)) THEN
           CALL SORT_NARRAY(go%hybm%dat, gx%hybm%dat, lrev)
        END IF
     END IF
  END IF

  IF (QDEF_NCVAR(go%ps)) THEN
     CALL RGMSG(substr, RGMLIC, ' ... PS   : '''//TRIM(go%ps%name)//'''')
     ! CHECK FOR LATM, LONM AND RESORT
     ALLOCATE(dil(go%ps%ndims), STAT=status)
     CALL ERRMSG(substr,status,5)
     ALLOCATE(svec(go%ps%ndims), STAT=status)
     CALL ERRMSG(substr,status,6)
     !
     olon = 0
     olat = 0
     dil(:) = 1
     svec(:) = 1
     DO i=1, go%ps%ndims
        ! CHECK DIM LENGTH AND NAME/ID
        IF (QDEF_NCVAR(go%lonm)) THEN
           IF (QCMP_NCDIM(go%ps%dim(i), go%lonm%dim(1)) > 1) THEN ! LONGITUDE
              olon = i
              dil(i) = go%ps%dim(i)%len
           END IF
        END IF
        IF (QDEF_NCVAR(go%latm)) THEN
           IF (QCMP_NCDIM(go%ps%dim(i), go%latm%dim(1)) > 1) THEN ! LATITUDE
              olat = i
              dil(i) = go%ps%dim(i)%len
           END IF
        END IF
     END DO

     IF (.NOT.lrev) THEN
        ! SPACE FOR INDICES
        CALL INIT_NARRAY(gx%ps%dat, go%ps%dat%n, go%ps%dat%dim, VTYPE_INT)
        CALL IDX_NCVAR(gx%ps)
     END IF

     vtype = gi%ps%dat%type

     IF (.NOT.lrev) THEN
        DO i=1, PRODUCT(dil)         ! LOOP OVER ALL FIELD VALUES
           CALL ELEMENT(dil,i,vec)   ! element vector in source array
           IF (olon > 0) THEN
              svec(olon) = INT(gx%lonm%dat%vi(vec(olon)))
           END IF
           IF (olat > 0) THEN
              svec(olat) = INT(gx%latm%dat%vi(vec(olat)))
           END IF
           !
           SELECT CASE (vtype)
           CASE (VTYPE_REAL)
              go%ps%dat%vr(i) = gi%ps%dat%vr(POSITION(dil, svec))
           CASE(VTYPE_DOUBLE)
              go%ps%dat%vd(i) = gi%ps%dat%vd(POSITION(dil, svec))
           CASE(VTYPE_INT)
              go%ps%dat%vi(i) = gi%ps%dat%vi(POSITION(dil, svec))
           CASE(VTYPE_BYTE)
              go%ps%dat%vb(i) = gi%ps%dat%vb(POSITION(dil, svec))
           CASE(VTYPE_CHAR)
              go%ps%dat%vc(i) = gi%ps%dat%vc(POSITION(dil, svec))
           CASE(VTYPE_UNDEF)
              CALL RGMSG(substr, RGMLE, &
                   'ARRAY OF UNDEFINED TYPE CANNOT BE SORTED !')
           CASE DEFAULT
              CALL RGMSG(substr, RGMLE, &
                   'ARRAY OF UNRECOGNIZED TYPE CANNOT BE SORTED !')
           END SELECT
           gx%ps%dat%vi(i) = POSITION(dil, svec)
           !
           DEALLOCATE(vec, STAT=status)
           CALL ERRMSG(substr,status,7)
           NULLIFY(vec)
        END DO   ! LOOP OVER ALL FIELD VALUES
     ELSE ! reverse
        DO i=1, PRODUCT(dil)         ! LOOP OVER ALL FIELD VALUES
           CALL ELEMENT(dil,i,vec)   ! element vector in source array
           IF (olon > 0) THEN
              svec(olon) = INT(gx%lonm%dat%vi(vec(olon)))
           END IF
           IF (olat > 0) THEN
                 svec(olat) = INT(gx%latm%dat%vi(vec(olat)))
           END IF
           !
           SELECT CASE (vtype)
           CASE (VTYPE_REAL)
              go%ps%dat%vr(POSITION(dil, svec)) = gi%ps%dat%vr(i)
           CASE(VTYPE_DOUBLE)
              go%ps%dat%vd(POSITION(dil, svec)) = gi%ps%dat%vd(i)
           CASE(VTYPE_INT)
              go%ps%dat%vi(POSITION(dil, svec)) = gi%ps%dat%vi(i)
           CASE(VTYPE_BYTE)
              go%ps%dat%vb(POSITION(dil, svec)) = gi%ps%dat%vb(i)
           CASE(VTYPE_CHAR)
              go%ps%dat%vc(POSITION(dil, svec)) = gi%ps%dat%vc(i)
           CASE(VTYPE_UNDEF)
              CALL RGMSG(substr, RGMLE, &
                   'ARRAY OF UNDEFINED TYPE CANNOT BE UN-SORTED !')
           CASE DEFAULT
              CALL RGMSG(substr, RGMLE, &
                   'ARRAY OF UNRECOGNIZED TYPE CANNOT BE UN-SORTED !')
           END SELECT
           !
           DEALLOCATE(vec, STAT=status)
           CALL ERRMSG(substr,status,8)
           NULLIFY(vec)
        END DO   ! LOOP OVER ALL FIELD VALUES
     END IF ! reverse ?

     DEALLOCATE(svec, dil, STAT=status)
     CALL ERRMSG(substr,status,9)
  END IF

  CALL RGMSG(substr, RGMLIC, '... END (UN-)SORTING GRID !')

CONTAINS

! --------------------------------------------
  SUBROUTINE IDX_NCVAR(var)

    IMPLICIT NONE

    ! I/O
    TYPE (t_ncvar), INTENT(INOUT) :: var

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr = 'IDX_NCVAR'
    INTEGER                      :: i
    INTEGER                      :: status
    CHARACTER(LEN=GRD_MAXSTRLEN) :: str

    str = TRIM(var%name)
    var%name  = TRIM(str)//'_IDX'
    var%xtype = NF90_INT
    var%id    = NULL_VARID
    IF (var%natts > 0) THEN
       DO i=1, var%natts
          CALL INIT_NCATT(var%att(i))
       END DO
       DEALLOCATE(var%att, STAT=status)
       CALL ERRMSG(substr,status,1)
       NULLIFY(var%att)
       var%natts = 0
    END IF
    CALL ADD_NCATT(var, 'RG_SORT_INDEX_OF', vs=TRIM(str))

   END SUBROUTINE IDX_NCVAR
! --------------------------------------------

END SUBROUTINE SORT_GEOHYBGRID
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE H2PSIG(psig, hya, hyb, ps, p0, lp)

  IMPLICIT NONE

  ! I/O
  TYPE (t_narray), INTENT(INOUT) :: psig
  TYPE (t_narray), INTENT(IN)    :: hya, hyb, ps, p0
  LOGICAL        , INTENT(IN)    :: lp   ! .true.  -> pressure axis
                                         ! .false. -> dimensionless axis
  ! LOCAL
  CHARACTER(LEN=*), PARAMETER    :: substr = 'H2PSIG'
  ! dimension length of linear array in psig
  INTEGER                        :: dim
  INTEGER                        :: status
  INTEGER                        :: i
  INTEGER                        :: i1,i2 ! little helpers
  TYPE (t_narray)                :: zhya, zhyb, zps, zp0     ! local arrays
  INTEGER, DIMENSION(:), POINTER :: vec                      ! element vector

  ! INIT
  NULLIFY(vec)
  CALL COPY_NARRAY(zhya, hya)
  CALL COPY_NARRAY(zhyb, hyb)
  CALL COPY_NARRAY(zps, ps)
  CALL COPY_NARRAY(zp0, p0)

  ! 1. CASE: HYBRID PRESSURE LEVELS
  IF ((hya%type /= VTYPE_UNDEF).AND.(hyb%type /= VTYPE_UNDEF)) THEN
     CALL RGMSG(substr, RGMLIC, '      -> HYBRID LEVELS ...')
     CALL RGMSG(substr, RGMLIC, '         ... HYA')
     CALL DOUBLE_NARRAY(zhya)
     CALL RGMSG(substr, RGMLIC, '         ... HYB')
     CALL DOUBLE_NARRAY(zhyb)
     CALL RGMSG(substr, RGMLIC, '         ... PS')
     CALL DOUBLE_NARRAY(zps)
     CALL RGMSG(substr, RGMLIC, '         ... P0')
     CALL DOUBLE_NARRAY(zp0)
     !
     ! ALLOCATE DIMENSION VECTOR
     dim = 1
     psig%n = hya%n + ps%n
     ALLOCATE(psig%dim(psig%n), STAT=status)
     CALL ERRMSG(substr,status,1)
     DO i=1, hya%n
        psig%dim(i) = zhya%dim(i)
        dim = dim * zhya%dim(i)
     END DO
     DO i=1, zps%n
        psig%dim(hya%n+i) = zps%dim(i)
        dim = dim * zps%dim(i)
     END DO
     !
     ! ALLOCATE DATA SPACE
     ALLOCATE(psig%vd(dim), STAT=status)
     CALL ERRMSG(substr,status,2)
     psig%type = VTYPE_DOUBLE
     !
     IF (lp) THEN  ! PRESSURE COORDINATES
        DO i=1, dim
           CALL ELEMENT(psig%dim, i, vec)  ! get element vector
           i1 = POSITION(zhya%dim,vec(1:zhya%n))
           i2 = POSITION(zps%dim, vec((zhya%n+1):))
           psig%vd(i) = (zhya%vd(i1) * &
                         zp0%vd(1) +   &
                         zhyb%vd(i1) * &
                         zps%vd(i2) )
        END DO
     ELSE          ! SIGMA COORDINATES
        DO i=1, dim
           CALL ELEMENT(psig%dim, i, vec)  ! get element vector
           i1 = POSITION(zhya%dim,vec(1:zhya%n))
           i2 = POSITION(zps%dim, vec((zhya%n+1):))
           psig%vd(i) = (zhya%vd(i1) * &
                         zp0%vd(1) +   &
                         zhyb%vd(i1) * &
                         zps%vd(i2) )/ &
                         zps%vd(i2)
        END DO
     END IF
  END IF

  ! 2. CASE: SIGMA LEVELS
  IF ((hya%type == VTYPE_UNDEF).AND.(hyb%type /= VTYPE_UNDEF)) THEN
     CALL RGMSG(substr, RGMLIC, '      -> SIGMA LEVELS ...')
     CALL RGMSG(substr, RGMLIC, '         ... HYB')
     CALL DOUBLE_NARRAY(zhyb)
     !
     IF (lp) THEN  ! PRESSURE COORDINATES
     CALL RGMSG(substr, RGMLIC, '         ... PS')
        CALL DOUBLE_NARRAY(zps)
        ! ALLOCATE DIMENSION VECTOR
        dim = 1
        psig%n = hyb%n + ps%n
        ALLOCATE(psig%dim(psig%n), STAT=status)
        CALL ERRMSG(substr,status,3)
        DO i=1, hyb%n
           psig%dim(i) = zhyb%dim(i)
           dim = dim * zhyb%dim(i)
        END DO
        DO i=1, zps%n
           psig%dim(hyb%n+i) = zps%dim(i)
           dim = dim * zps%dim(i)
        END DO
        !
        ! ALLOCATE DATA SPACE
        ALLOCATE(psig%vd(dim), STAT=status)
        CALL ERRMSG(substr,status,4)
        psig%type = VTYPE_DOUBLE
        !
        DO i=1, dim
           CALL ELEMENT(psig%dim, i, vec)  ! get element vector
           i1 = POSITION(zhyb%dim,vec(1:zhyb%n))
           i2 = POSITION(zps%dim, vec((zhyb%n+1):))
           psig%vd(i) = (zhyb%vd(i1) * zps%vd(i2))
        END DO
     ELSE          ! SIGMA COORDINATES
        ! ALLOCATE DIMENSION VECTOR
        dim = 1
        psig%n = hyb%n
        ALLOCATE(psig%dim(psig%n), STAT=status)
        CALL ERRMSG(substr,status,5)
        psig%dim(:) = zhyb%dim(:)
        DO i=1, psig%n
           dim = dim * psig%dim(i)
        END DO
        !
        ! ALLOCATE DATA SPACE
        ALLOCATE(psig%vd(dim), STAT=status)
        CALL ERRMSG(substr,status,6)
        psig%type = VTYPE_DOUBLE
        !
        psig%vd(:) = zhyb%vd(:)/MAXVAL(zhyb%vd)
     END IF
  END IF

  ! 3.CASE: CONST. PRESSURE LEVELS
  IF ((hya%type /= VTYPE_UNDEF).AND.(hyb%type == VTYPE_UNDEF)) THEN
     CALL RGMSG(substr, RGMLIC, '      -> CONSTANT PRESSURE LEVELS ...')
     CALL RGMSG(substr, RGMLIC, '         ... HYA')
     CALL DOUBLE_NARRAY(zhya)
     !
     ! ALLOCATE DIMENSION VECTOR
     dim = 1
     psig%n = hya%n
     ALLOCATE(psig%dim(psig%n), STAT=status)
     CALL ERRMSG(substr,status,7)
     psig%dim(:) = zhya%dim(:)
     DO i=1, psig%n
        dim = dim * psig%dim(i)
     END DO
     !
     ! ALLOCATE DATA SPACE
     ALLOCATE(psig%vd(dim), STAT=status)
     CALL ERRMSG(substr,status,8)
     psig%type = VTYPE_DOUBLE
     !
     IF (lp) THEN  ! PRESSURE COORDINATES
        CALL RGMSG(substr, RGMLIC, '         ... P0')
        CALL DOUBLE_NARRAY(zp0)
        psig%vd(:) = zhya%vd(:)*zp0%vd(1)
     ELSE          ! SIGMA COORDINATES
        psig%vd(:) = zhya%vd(:)/MAXVAL(zhya%vd)
     END IF
  END IF

  ! CLEAN UP
  IF (ASSOCIATED(vec)) THEN
     DEALLOCATE(vec, STAT=status)
     CALL ERRMSG(substr,status,9)
     NULLIFY(vec)
  END IF
  CALL INIT_NARRAY(zhya)
  CALL INIT_NARRAY(zhyb)
  CALL INIT_NARRAY(zps)
  CALL INIT_NARRAY(zp0)

END SUBROUTINE H2PSIG
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE H2PSIG_3d(psig, xs, ys, hya, hyb, ps, p0, lp)

  IMPLICIT NONE

  ! I/O
!!$  REAL(dp), DIMENSION(:,:,:), POINTER, INTENT(INOUT) :: psig
  REAL(dp), DIMENSION(:,:,:), POINTER                :: psig ! INOUT
  INTEGER,                             INTENT(IN)    :: xs,ys 
  TYPE (t_narray),                     INTENT(IN)    :: hya, hyb, p0
  REAL(dp), DIMENSION(:,:),            INTENT(IN)    :: ps
  LOGICAL,                             INTENT(IN)    :: lp   ! .true.  -> pressure axis
                                         ! .false. -> dimensionless axis
  ! LOCAL
  CHARACTER(LEN=*), PARAMETER    :: substr = 'H2PSIG_3d'
  INTEGER                        :: i1,i2,i3 ! little helpers
  TYPE (t_narray)                :: zhya, zhyb,  zp0     ! local arrays
  REAL(dp), DIMENSION(xs,ys)     :: zps

  ! INIT
  CALL COPY_NARRAY(zhya, hya)
  CALL COPY_NARRAY(zhyb, hyb)
  CALL COPY_NARRAY(zp0, p0)
  zps = ps

  ! 1. CASE: HYBRID PRESSURE LEVELS
  IF ((hya%type /= VTYPE_UNDEF).AND.(hyb%type /= VTYPE_UNDEF)) THEN
     CALL RGMSG(substr, RGMLIC, '      -> HYBRID LEVELS ...')
     CALL RGMSG(substr, RGMLIC, '         ... HYA')
     CALL DOUBLE_NARRAY(zhya)
     CALL RGMSG(substr, RGMLIC, '         ... HYB')
     CALL DOUBLE_NARRAY(zhyb)
     !CALL RGMSG(substr, RGMLIC, '         ... PS')
     !CALL DOUBLE_NARRAY(zps)
     CALL RGMSG(substr, RGMLIC, '         ... P0')
     IF (zp0%type == VTYPE_UNDEF) THEN
        CALL RGMSG(substr, RGMLE, '    ... P0 not available but required')
     END IF
     CALL DOUBLE_NARRAY(zp0)
     !
     ! ALLOCATE DATA SPACE
     ALLOCATE(psig(xs,ys,zhya%dim(1)))
     !
     IF (lp) THEN  ! PRESSURE COORDINATES
        DO i1 = 1, xs
           DO i2 = 1, ys
              DO i3 = 1, SIZE(psig,3)
                 psig(i1,i2,i3) = zhya%vd(i3) * zp0%vd(1)  &
                                + zhyb%vd(i3) * zps(i1,i2)
              END DO
           END DO
        END DO
     ELSE          ! SIGMA COORDINATES
        DO i1 = 1, xs
           DO i2 = 1, ys
              DO i3 = 1, SIZE(psig,3)
                 psig(i1,i2,i3) = ( zhya%vd(i3) * zp0%vd(1)  &
                                + zhyb%vd(i3) * zps(i1,i2) ) / zps(i1,i2)
              END DO
           END DO
        END DO
     END IF
  END IF

  ! 2. CASE: SIGMA LEVELS
  IF ((hya%type == VTYPE_UNDEF).AND.(hyb%type /= VTYPE_UNDEF)) THEN
     CALL RGMSG(substr, RGMLIC, '      -> SIGMA LEVELS ...')
     CALL RGMSG(substr, RGMLIC, '         ... HYB')
     CALL DOUBLE_NARRAY(zhyb)
     !
     ALLOCATE(psig(xs,ys,zhyb%dim(1)))
     IF (lp) THEN  ! PRESSURE COORDINATES
        !
        ! ALLOCATE DATA SPACE
        !
        DO i1 = 1, xs
           DO i2 = 1, ys
              DO i3 = 1, SIZE(psig,3)
                 psig(i1,i2,i3) = zhyb%vd(i3) * zps(i1,i2)
              END DO
           END DO
        END DO
     ELSE          ! SIGMA COORDINATES
        !
        DO i1 = 1, xs
           DO i2 = 1, ys
              DO i3 = 1, SIZE(psig,3)
                 psig(i1,i2,i3) = zhyb%vd(i3) /MAXVAL(zhyb%vd)
              END DO
           END DO
        END DO
     END IF
  END IF

  ! 3.CASE: CONST. PRESSURE LEVELS
  IF ((hya%type /= VTYPE_UNDEF).AND.(hyb%type == VTYPE_UNDEF)) THEN
     CALL RGMSG(substr, RGMLIC, '      -> CONSTANT PRESSURE LEVELS ...')
     CALL RGMSG(substr, RGMLIC, '         ... HYA')
     CALL DOUBLE_NARRAY(zhya)
     !
     !
     ALLOCATE(psig(xs,ys,zhya%dim(1)))
     IF (lp) THEN  ! PRESSURE COORDINATES
        !
        ! ALLOCATE DATA SPACE
        !
        CALL RGMSG(substr, RGMLIC, '         ... p0')
        CALL DOUBLE_NARRAY(zp0)
        DO i1 = 1, xs
           DO i2 = 1, ys
              DO i3 = 1, SIZE(psig,3)
                 psig(i1,i2,i3) = zhya%vd(i3) * zp0%vd(1)
              END DO
           END DO
        END DO
     ELSE          ! SIGMA COORDINATES
        !
        DO i1 = 1, xs
           DO i2 = 1, ys
              DO i3 = 1, SIZE(psig,3)
                 psig(i1,i2,i3) = zhya%vd(i3) /MAXVAL(zhya%vd)
              END DO
           END DO
        END DO
     END IF
  END IF

  CALL INIT_NARRAY(zhya)
  CALL INIT_NARRAY(zhyb)
  CALL INIT_NARRAY(zp0)

END SUBROUTINE H2PSIG_3D
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE COMPLETE_GEOHYBGRID(g, gx, oranges)

  ! Note: g must be already sorted

  IMPLICIT NONE

  ! I/O
  TYPE (t_geohybgrid),      INTENT(INOUT)           :: g
  TYPE (t_geohybgrid),      INTENT(INOUT), OPTIONAL :: gx   ! SORT INDICES
  REAL(DP), DIMENSION(4,2), INTENT(IN)   , OPTIONAL :: oranges

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'COMPLETE_GEOHYBGRID'
  REAL(DP), DIMENSION(4,2)    :: ranges
  LOGICAL                     :: llonm,  llatm,  lloni,  llati
  LOGICAL                     :: lrlonm, lrlatm, lrloni, lrlati
  LOGICAL                     :: lclonm, lclatm, lcloni, lclati

  IF (PRESENT(oranges)) THEN
     ranges = oranges
  ELSE
     ranges = g%ranges
  ENDIF

  CALL RGMSG(substr, RGMLI, 'FILE '''//TRIM(g%file)//'''')
  CALL RGMSG(substr, RGMLIC,'STEP ',g%t,' ')

  CALL RGMSG(substr, RGMLIC,'  loni <<--->> lonm ...')
  CALL IMMI_NARRAY(g%loni%dat, g%lonm%dat, .FALSE., lloni, llonm)
  CALL IMMI_NCVAR(g%loni, g%lonm)
  IF (PRESENT(gx)) THEN
     CALL IMMI_NARRAY_IDX(gx%loni%dat, gx%lonm%dat)
     CALL IMMI_NCVAR(gx%loni, gx%lonm)
  END IF
  ! CORRECT FOR RANGE
  CALL RNGADJ_NARRAY(g%loni%dat, ranges(1,:))

  CALL RGMSG(substr, RGMLIC,'  lati <<--->> latm ...')
  CALL IMMI_NARRAY(g%lati%dat, g%latm%dat, .FALSE., llati, llatm)
  CALL IMMI_NCVAR(g%lati, g%latm)
  IF (PRESENT(gx)) THEN
     CALL IMMI_NARRAY_IDX(gx%lati%dat, gx%latm%dat)
     CALL IMMI_NCVAR(gx%lati, gx%latm)
  END IF
  ! CORRECT FOR RANGE
  CALL RNGADJ_NARRAY(g%lati%dat, ranges(2,:))

  CALL RGMSG(substr, RGMLIC,'  cloni <<--->> clonm ...')
  CALL IMMI_NARRAY(g%cloni%dat, g%clonm%dat, .TRUE., lcloni, lclonm)

  CALL RGMSG(substr, RGMLIC,'  clati <<--->> clatm ...')
  CALL IMMI_NARRAY(g%clati%dat, g%clatm%dat, .TRUE., lclati, lclatm)

  CALL RGMSG(substr, RGMLIC,'  rloni <<--->> rlonm ...')
  CALL IMMI_NARRAY(g%rloni%dat, g%rlonm%dat, .FALSE., lrloni, lrlonm)
  CALL IMMI_NCVAR(g%rloni, g%rlonm)
  ! CORRECT FOR RANGE
  CALL RNGADJ_NARRAY(g%rloni%dat, ranges(1,:))

  CALL RGMSG(substr, RGMLIC,'  rlati <<--->> rlatm ...')
  CALL IMMI_NARRAY(g%rlati%dat, g%rlatm%dat, .FALSE., lrlati, lrlatm)
  CALL IMMI_NCVAR(g%rlati, g%rlatm)
  ! CORRECT FOR RANGE
  CALL RNGADJ_NARRAY(g%rlati%dat, ranges(2,:))

  IF (.NOT. llonm .AND. .NOT. lloni) THEN
     IF (.NOT. lcloni .AND. lrloni) THEN
        CALL RGMSG(substr, RGMLIC,'  rloni ---> cloni ...')
        CALL RGMSG(substr, RGMLIC,'and  rlati ---> clati ...')
        CALL IMMI_CLONI_CLATI(g)
     END IF
  END IF
  IF (.NOT. llonm .AND. .NOT. lloni) THEN
     IF (.NOT. lclonm .AND. lrlonm) THEN
        CALL RGMSG(substr, RGMLIC,'  rlonm ---> clonm ...')
        CALL RGMSG(substr, RGMLIC,'and  rlatm ---> clatm ...')
        CALL IMMI_CLONM_CLATM(g)
     END IF
  END IF
  CALL RGMSG(substr, RGMLIC,'  hyai <<--->> hyam ...')
  CALL IMMI_NARRAY(g%hyai%dat, g%hyam%dat)
  CALL IMMI_NCVAR(g%hyai, g%hyam)
  IF (PRESENT(gx)) THEN
     CALL IMMI_NARRAY_IDX(gx%hyai%dat, gx%hyam%dat)
     CALL IMMI_NCVAR(gx%hyai, gx%hyam)
  END IF
  ! CORRECT FOR RANGE
  CALL RNGADJ_NARRAY(g%hyai%dat, ranges(3,:), .false.)

  CALL RGMSG(substr, RGMLIC,'  hybi <<--->> hybm ...')
  CALL IMMI_NARRAY(g%hybi%dat, g%hybm%dat)
  CALL IMMI_NCVAR(g%hybi, g%hybm)
  IF (PRESENT(gx)) THEN
     CALL IMMI_NARRAY_IDX(gx%hybi%dat, gx%hybm%dat)
     CALL IMMI_NCVAR(gx%hybi, gx%hybm)
  END IF
  ! CORRECT FOR RANGE
  CALL RNGADJ_NARRAY(g%hybi%dat, ranges(4,:))

  CALL RGMSG(substr, RGMLIC,'  timei <<--->> timem ...')
  CALL IMMI_NARRAY(g%timei%dat, g%timem%dat)
  CALL IMMI_NCVAR(g%timei, g%timem)
  IF (PRESENT(gx)) THEN
     CALL IMMI_NARRAY_IDX(gx%timei%dat, gx%timem%dat)
     CALL IMMI_NCVAR(gx%timei, gx%timem)
  END IF

CONTAINS

!--------------------------------------------
SUBROUTINE IMMI_NARRAY_IDX(nai, nam)

  IMPLICIT NONE

  ! I/O
  TYPE (t_narray), INTENT(INOUT) :: nai, nam

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'IMMI_NARRAY_IDX'
  INTEGER                     :: itype, mtype
  INTEGER                     :: i
  INTEGER                     :: status

  itype = nai%type
  mtype = nam%type

  IF ((itype == VTYPE_INT).AND.(mtype == VTYPE_INT)) RETURN
  IF ((itype == VTYPE_UNDEF).AND.(mtype == VTYPE_UNDEF)) RETURN


  IF ((itype == VTYPE_INT).AND.(mtype == VTYPE_UNDEF)) THEN
     ! I -> M
     IF (nai%n /= 1) THEN
        CALL RGMSG(substr, RGMLE, &
             'INTERFACE ARRAY MUST BE 1-DIMENSIONAL !')
     END IF
     nam%n = nai%n
     ALLOCATE(nam%dim(nam%n), STAT=status)
     CALL ERRMSG(substr,status,1)
     nam%dim(1) = nai%dim(1) - 1
     ALLOCATE(nam%vi(nam%dim(1)), STAT=status)
     nam%type = VTYPE_INT
     CALL ERRMSG(substr,status,2)
     DO i=1, nai%dim(1) - 1
        IF (nai%vi(i) > nai%vi(nai%dim(1))) THEN
           nam%vi(i) = nai%vi(i) - 1
        ELSE
           nam%vi(i) = nai%vi(i)
        END IF
     END DO
     RETURN
  END IF

  IF ((mtype == VTYPE_INT).AND.(itype == VTYPE_UNDEF)) THEN
     ! M -> I
     IF (nam%n /= 1) THEN
        CALL RGMSG(substr, RGMLE, &
             'MID ARRAY MUST BE 1-DIMENSIONAL !')
     END IF
     nai%n = nam%n
     ALLOCATE(nai%dim(nam%n), STAT=status)
     CALL ERRMSG(substr,status,3)
     nai%dim(1) = nam%dim(1) + 1
     ALLOCATE(nai%vi(nai%dim(1)), STAT=status)
     nai%type = VTYPE_INT
     CALL ERRMSG(substr,status,4)
     nai%vi(nai%dim(1)) = nam%vi(nam%dim(1))
     DO i=1, nam%dim(1)
        IF (nam%vi(i) >= nam%vi(nam%dim(1))) THEN
           nai%vi(i) = nam%vi(i) + 1
        ELSE
           nai%vi(i) = nam%vi(i)
        END IF
     END DO
     IF (nai%vi(nai%dim(1)-1) == MAXVAL(nai%vi)) THEN
        i = INT(nai%vi(nai%dim(1)-1))
        nai%vi(nai%dim(1)-1) = nai%vi(nai%dim(1))
        nai%vi(nai%dim(1)) = INT(i, I8)
     END IF
     RETURN
  END IF

  CALL RGMSG(substr, RGMLE, &
       'N-ARRAYS MUST BE OF TYPE INTEGER OR UNDEFINED !')

END SUBROUTINE IMMI_NARRAY_IDX
! --------------------------------------------

  SUBROUTINE IMMI_CLONI_CLATI(g)

  ! in this subroutine cloni is calculated from rloni
  !
  IMPLICIT NONE

  ! I/O
  TYPE (t_geohybgrid), INTENT(INOUT) :: g

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'IMMI_CLONI_CLATI'
  INTEGER                     :: i,j,n

  IF (g%rloni%dat%type == VTYPE_UNDEF) RETURN
  IF (g%rlati%dat%type == VTYPE_UNDEF) RETURN
  IF (g%cloni%dat%type /= VTYPE_UNDEF) RETURN
  IF (g%clati%dat%type /= VTYPE_UNDEF) RETURN

  ! calculation only possible, of coordinates of rotated pole are known
  IF (g%pollon%dat%type == VTYPE_UNDEF .OR. g%pollat%dat%type == VTYPE_UNDEF &
       .OR. g%polgam%dat%type == VTYPE_UNDEF) THEN
     CALL RGMSG(substr, RGMLIC, '  ... rotated pole not defined')
     CALL RGMSG(substr, RGMLIC, '  ... calculation of cloni not possible')
     RETURN
  ENDIF

  IF (g%rlati%dat%type /= g%rloni%dat%type) RETURN
  
  CALL RGMSG(substr, RGMLIC, &
       '  ... rotated coordinates -> curvilinear interfaces')
  IF (g%rlati%dat%n /= 1) THEN
     CALL RGMSG(substr, RGMLE, &
          'rotated latitude IS ',g%rlati%dat%n,'-DIMENSIONAL !')
  END IF
  IF (g%rloni%dat%n /= 1) THEN
     CALL RGMSG(substr, RGMLE, &
          'rotated longitude IS ',g%rloni%dat%n,'-DIMENSIONAL !')
  END IF
  g%cloni%name  = 'lon_I'
  g%clati%name  = 'lon_I'
  g%cloni%id    = NULL_VARID
  g%clati%id    = NULL_VARID
  g%cloni%ndims = 2
  g%clati%ndims = 2
  ALLOCATE(g%cloni%dim(g%cloni%ndims))
  ALLOCATE(g%clati%dim(g%clati%ndims))
  g%cloni%dim(1) = g%rloni%dim(1)
  g%cloni%dim(2) = g%rlati%dim(1)
  g%clati%dim(1) = g%rloni%dim(1)
  g%clati%dim(2) = g%rlati%dim(1)
  g%cloni%dim(1)%name  = 'lon_I'
  g%cloni%dim(1)%id    = NULL_DIMID
  g%cloni%dim(1)%len   = g%rloni%dim(1)%len
  g%cloni%dim(1)%fuid  = .false.
  g%cloni%dim(1)%varid = NULL_VARID
  g%cloni%dim(2)%name  = 'lat_I'
  g%cloni%dim(2)%id    = NULL_DIMID
  g%cloni%dim(2)%len   = g%rlati%dim(1)%len
  g%cloni%dim(2)%fuid  = .false.
  g%cloni%dim(2)%varid = NULL_VARID
  g%clati%dim(1)%name  = 'lon_I'
  g%clati%dim(1)%id    = NULL_DIMID
  g%clati%dim(1)%len   = g%rloni%dim(1)%len
  g%clati%dim(1)%fuid  = .false.
  g%clati%dim(1)%varid = NULL_VARID
  g%clati%dim(2)%name  = 'lat_I'
  g%clati%dim(2)%id    = NULL_DIMID
  g%clati%dim(2)%len   = g%rlati%dim(1)%len
  g%clati%dim(2)%fuid  = .false.
  g%clati%dim(2)%varid = NULL_VARID
  
  SELECT CASE(g%rlati%dat%type)
  CASE(VTYPE_REAL)
     g%cloni%xtype = NF90_FLOAT
     CALL INIT_NARRAY(g%cloni%dat, g%cloni%ndims &
       , (/g%cloni%dim(1)%len,g%cloni%dim(2)%len/) &
       ,VTYPE_REAL)
     g%clati%xtype = NF90_FLOAT
     CALL INIT_NARRAY(g%clati%dat, g%clati%ndims &
       , (/g%clati%dim(1)%len,g%clati%dim(2)%len/) &
       ,VTYPE_REAL)

     DO i = 1, g%clati%dim(1)%len
        DO j=1, g%clati%dim(2)%len
           n = (j-1) * g%clati%dim(1)%len + i
           g%cloni%dat%vr(n) = REAL(rlarot2rla(REAL(g%rlati%dat%vr(j),dp) &
                , REAL(g%rloni%dat%vr(i),dp)                              &
                , REAL(g%pollat%dat%vr(1),dp), REAL(g%pollon%dat%vr(1),dp)&
                , REAL(g%polgam%dat%vr(1),dp)) ,sp)
           g%clati%dat%vr(n) = REAL(phirot2phi(REAL(g%rlati%dat%vr(j),dp) &
                , REAL(g%rloni%dat%vr(i),dp)                              &
                , REAL(g%pollat%dat%vr(1),dp), REAL(g%polgam%dat%vr(1),dp)),sp)
        END DO
     END DO

       ! DETERMINE MIN / MAX LONGITUDE
     g%minmaxlonlat(1,1) = MINVAL(g%cloni%dat%vr)
     g%minmaxlonlat(1,2) = MAXVAL(g%cloni%dat%vr)
     ! DETERMIN MIN / MAX LATITUDE
     g%minmaxlonlat(2,1) = MINVAL(g%clati%dat%vr)
     g%minmaxlonlat(2,2) = MAXVAL(g%clati%dat%vr)

     CALL RGMSG(substr, RGMLIC, '  ... O.K.')

     RETURN
  CASE(VTYPE_DOUBLE)
     g%cloni%xtype = NF90_DOUBLE
     CALL INIT_NARRAY(g%cloni%dat, g%cloni%ndims   &
       , (/g%cloni%dim(1)%len,g%cloni%dim(2)%len/) &
       ,VTYPE_DOUBLE)
     g%clati%xtype = NF90_DOUBLE
     CALL INIT_NARRAY(g%clati%dat, g%clati%ndims   &
       , (/g%clati%dim(1)%len,g%clati%dim(2)%len/) &
       ,VTYPE_DOUBLE)

     DO i = 1, g%clati%dim(1)%len*g%clati%dim(2)%len
        g%cloni%dat%vd(i) = rlarot2rla(g%rlati%dat%vd(i)      &
             , g%rloni%dat%vd(i), REAL(g%pollat%dat%vr(1),dp) &
             , REAL(g%pollon%dat%vr(1),dp), REAL(g%polgam%dat%vr(1),dp))
        g%clati%dat%vd(i) = phirot2phi(g%rlati%dat%vd(i)      &
             , g%rloni%dat%vd(i), REAL(g%pollat%dat%vr(1),dp) &
             , REAL(g%polgam%dat%vr(1),dp))
     END DO

       ! DETERMINE MIN / MAX LONGITUDE
     g%minmaxlonlat(1,1) = MINVAL(g%cloni%dat%vd)
     g%minmaxlonlat(1,2) = MAXVAL(g%cloni%dat%vd)
     ! DETERMIN MIN / MAX LATITUDE
     g%minmaxlonlat(2,1) = MINVAL(g%clati%dat%vd)
     g%minmaxlonlat(2,2) = MAXVAL(g%clati%dat%vd)

     CALL RGMSG(substr, RGMLIC, '  ... O.K.')

     RETURN
  CASE(VTYPE_INT)
     CALL RGMSG(substr, RGMLE, &
          'CURVILINEAR INTERFACES OF TYPE INTEGER ARE NOT SUPPORTED !')
     RETURN
  CASE(VTYPE_BYTE)
     CALL RGMSG(substr, RGMLE, &
          'CURVILINEAR INTERFACES OF TYPE BYTE ARE NOT SUPPORTED !')
  CASE(VTYPE_CHAR)
     CALL RGMSG(substr, RGMLE, &
          'CURVILINEAR INTERFACES OF TYPE CHAR ARE NOT SUPPORTED !')
  CASE(VTYPE_UNDEF)
     ! OK, nam -> nai
  CASE DEFAULT
     CALL RGMSG(substr, RGMLE, &
          'UNRECOGNIZED TYPE OF CURVILINEAR INTERFACE ARRAY !')
  END SELECT

END SUBROUTINE  IMMI_CLONI_CLATI

  SUBROUTINE IMMI_CLONM_CLATM(g)

  ! in this subroutine cloni is calculated from rloni
  !
  IMPLICIT NONE

  ! I/O
  TYPE (t_geohybgrid), INTENT(INOUT) :: g

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'IMMI_CLONM_CLATM'
  INTEGER                     :: i,j,n

  IF (g%rlonm%dat%type == VTYPE_UNDEF) RETURN
  IF (g%rlatm%dat%type == VTYPE_UNDEF) RETURN
  IF (g%clonm%dat%type /= VTYPE_UNDEF) RETURN
  IF (g%clatm%dat%type /= VTYPE_UNDEF) RETURN

  ! calculation only possible, of coordinates of rotated pole are known
  IF (g%pollon%dat%type == VTYPE_UNDEF .OR. g%pollat%dat%type == VTYPE_UNDEF &
       .OR. g%polgam%dat%type == VTYPE_UNDEF) THEN
     CALL RGMSG(substr, RGMLIC, '  ... rotated pole not defined')
     CALL RGMSG(substr, RGMLIC, '  ... calculation of cloni not possible')
     RETURN
  ENDIF

  IF (g%rlatm%dat%type /= g%rlonm%dat%type) RETURN
  
  CALL RGMSG(substr, RGMLIC, &
       '  ... rotated coordinates -> curvilinear interfaces')
  IF (g%rlatm%dat%n /= 1) THEN
     CALL RGMSG(substr, RGMLE, &
          'rotated latitude IS ',g%rlatm%dat%n,'-DIMENSIONAL !')
  END IF
  IF (g%rlonm%dat%n /= 1) THEN
     CALL RGMSG(substr, RGMLE, &
          'rotated longitude IS ',g%rlonm%dat%n,'-DIMENSIONAL !')
  END IF
  g%clonm%name  = 'lon_I'
  g%clatm%name  = 'lon_I'
  g%clonm%id    = NULL_VARID
  g%clatm%id    = NULL_VARID
  g%clonm%ndims = 2
  g%clatm%ndims = 2
  ALLOCATE(g%clonm%dim(g%clonm%ndims))
  ALLOCATE(g%clatm%dim(g%clatm%ndims))
  g%clonm%dim(1) = g%rlonm%dim(1)
  g%clonm%dim(2) = g%rlatm%dim(1)
  g%clatm%dim(1) = g%rlonm%dim(1)
  g%clatm%dim(2) = g%rlatm%dim(1)
  g%clonm%dim(1)%name  = 'lon_I'
  g%clonm%dim(1)%id    = NULL_DIMID
  g%clonm%dim(1)%len   = g%rlonm%dim(1)%len
  g%clonm%dim(1)%fuid  = .false.
  g%clonm%dim(1)%varid = NULL_VARID
  g%clonm%dim(2)%name  = 'lat_I'
  g%clonm%dim(2)%id    = NULL_DIMID
  g%clonm%dim(2)%len   = g%rlatm%dim(1)%len
  g%clonm%dim(2)%fuid  = .false.
  g%clonm%dim(2)%varid = NULL_VARID
  g%clatm%dim(1)%name  = 'lon_I'
  g%clatm%dim(1)%id    = NULL_DIMID
  g%clatm%dim(1)%len   = g%rlonm%dim(1)%len
  g%clatm%dim(1)%fuid  = .false.
  g%clatm%dim(1)%varid = NULL_VARID
  g%clatm%dim(2)%name  = 'lat_I'
  g%clatm%dim(2)%id    = NULL_DIMID
  g%clatm%dim(2)%len   = g%rlatm%dim(1)%len
  g%clatm%dim(2)%fuid  = .false.
  g%clatm%dim(2)%varid = NULL_VARID
  
  SELECT CASE(g%rlatm%dat%type)
  CASE(VTYPE_REAL)
     g%clonm%xtype = NF90_FLOAT
     CALL INIT_NARRAY(g%clonm%dat, g%clonm%ndims &
       , (/g%clonm%dim(1)%len,g%clonm%dim(2)%len/) &
       ,VTYPE_REAL)
     g%clatm%xtype = NF90_FLOAT
     CALL INIT_NARRAY(g%clatm%dat, g%clatm%ndims &
       , (/g%clatm%dim(1)%len,g%clatm%dim(2)%len/) &
       ,VTYPE_REAL)

     DO i = 1, g%clatm%dim(1)%len
        DO j=1, g%clatm%dim(2)%len
           n = (j-1) * g%clatm%dim(1)%len + i
           g%clonm%dat%vr(n) = REAL(rlarot2rla(REAL(g%rlatm%dat%vr(j),dp) &
                , REAL(g%rlonm%dat%vr(i),dp)                         &
                , REAL(g%pollat%dat%vr(1),dp), REAL(g%pollon%dat%vr(1),dp)&
                , REAL(g%polgam%dat%vr(1),dp)) ,sp)
           g%clatm%dat%vr(n) = REAL(phirot2phi(REAL(g%rlatm%dat%vr(j),dp) &
                , REAL(g%rlonm%dat%vr(i),dp), REAL(g%pollat%dat%vr(1),dp) &
                , REAL(g%polgam%dat%vr(1),dp)),sp)
        END DO
     END DO

       ! DETERMINE MIN / MAX LONGITUDE
     g%minmaxlonlat(1,1) = MINVAL(g%clonm%dat%vr)
     g%minmaxlonlat(1,2) = MAXVAL(g%clonm%dat%vr)
     ! DETERMIN MIN / MAX LATITUDE
     g%minmaxlonlat(2,1) = MINVAL(g%clatm%dat%vr)
     g%minmaxlonlat(2,2) = MAXVAL(g%clatm%dat%vr)

     CALL RGMSG(substr, RGMLIC, '  ... O.K.')

     RETURN
  CASE(VTYPE_DOUBLE)
     g%clonm%xtype = NF90_DOUBLE
     CALL INIT_NARRAY(g%clonm%dat, g%clonm%ndims &
       , (/g%clonm%dim(1)%len,g%clonm%dim(2)%len/) &
       ,VTYPE_DOUBLE)
     g%clatm%xtype = NF90_DOUBLE
     CALL INIT_NARRAY(g%clatm%dat, g%clatm%ndims &
       , (/g%clatm%dim(1)%len,g%clatm%dim(2)%len/) &
       ,VTYPE_DOUBLE)

     DO i = 1, g%clatm%dim(1)%len*g%clatm%dim(2)%len
        g%clonm%dat%vd(i) = rlarot2rla(g%rlatm%dat%vd(i)      &
             , g%rlonm%dat%vd(i), REAL(g%pollat%dat%vr(1),dp) &
             , REAL(g%pollon%dat%vr(1),dp), REAL(g%polgam%dat%vr(1),dp))
        g%clatm%dat%vd(i) = phirot2phi(g%rlatm%dat%vd(i)      &
             , g%rlonm%dat%vd(i), REAL(g%pollat%dat%vr(1),dp) &
             , REAL(g%polgam%dat%vr(1),dp))
     END DO

       ! DETERMINE MIN / MAX LONGITUDE
     g%minmaxlonlat(1,1) = MINVAL(g%clonm%dat%vd)
     g%minmaxlonlat(1,2) = MAXVAL(g%clonm%dat%vd)
     ! DETERMIN MIN / MAX LATITUDE
     g%minmaxlonlat(2,1) = MINVAL(g%clatm%dat%vd)
     g%minmaxlonlat(2,2) = MAXVAL(g%clatm%dat%vd)

     CALL RGMSG(substr, RGMLIC, '  ... O.K.')

     RETURN
  CASE(VTYPE_INT)
     CALL RGMSG(substr, RGMLE, &
          'CURVILINEAR INTERFACES OF TYPE INTEGER ARE NOT SUPPORTED !')
     RETURN
  CASE(VTYPE_BYTE)
     CALL RGMSG(substr, RGMLE, &
          'CURVILINEAR INTERFACES OF TYPE BYTE ARE NOT SUPPORTED !')
  CASE(VTYPE_CHAR)
     CALL RGMSG(substr, RGMLE, &
          'CURVILINEAR INTERFACES OF TYPE CHAR ARE NOT SUPPORTED !')
  CASE(VTYPE_UNDEF)
     ! OK, nam -> nai
  CASE DEFAULT
     CALL RGMSG(substr, RGMLE, &
          'UNRECOGNIZED TYPE OF CURVILINEAR INTERFACE ARRAY !')
  END SELECT

END SUBROUTINE  IMMI_CLONM_CLATM
! --------------------------------------------
END SUBROUTINE COMPLETE_GEOHYBGRID

! --------------------------------------------
  SUBROUTINE IMMI_NCVAR(vari, varm)

    IMPLICIT NONE

    ! I/O
    TYPE (t_ncvar), INTENT(INOUT) :: vari, varm

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr = 'IMMI_NCVAR'
    INTEGER                      :: i
    INTEGER                      :: status
    CHARACTER(LEN=GRD_MAXSTRLEN) :: str
    INTEGER                      :: vtype

    IF ((TRIM(vari%name) == '').AND.(TRIM(varm%name) == '')) RETURN
    IF ((TRIM(vari%name) /= '').AND.(TRIM(varm%name) /= '')) RETURN

    IF (TRIM(vari%name) /= '') THEN    ! i -> m
       str = TRIM(vari%name)
       varm%name  = TRIM(str)//'_M'
       varm%id    = NULL_VARID
       vtype = varm%dat%type
       SELECT CASE(vtype)
       CASE(VTYPE_REAL)
          varm%xtype = NF90_FLOAT
       CASE(VTYPE_DOUBLE)
          varm%xtype = NF90_DOUBLE
       CASE(VTYPE_INT)
          varm%xtype = NF90_INT
       CASE(VTYPE_BYTE)
          varm%xtype = NF90_BYTE
       CASE(VTYPE_CHAR)
          varm%xtype = NF90_CHAR
       CASE(VTYPE_UNDEF)
          CALL RGMSG(substr, RGMLE, &
               'TYPE OF VARIABLE '''//TRIM(varm%name)//''' IS UNDEFINED !')
       CASE DEFAULT
          CALL RGMSG(substr, RGMLE, &
               'TYPE OF VARIABLE '''//TRIM(varm%name)//''' IS UNRECOGNIZED !')
       END SELECT
       IF (vari%ndims > 0) THEN
          varm%ndims = vari%ndims
          ALLOCATE(varm%dim(varm%ndims), STAT=status)
          CALL ERRMSG(substr,status,1)
          DO i=1, varm%ndims
             varm%dim(i)%name  = TRIM(vari%dim(i)%name)//'_M'
             varm%dim(i)%id    = NULL_DIMID
             varm%dim(i)%len   = vari%dim(i)%len - 1
             varm%dim(i)%fuid  = vari%dim(i)%fuid
             varm%dim(i)%varid = NULL_VARID
          END DO
       END IF
       IF (varm%natts > 0) THEN
          DO i=1, varm%natts
             CALL INIT_NCATT(varm%att(i))
          END DO
          DEALLOCATE(varm%att, STAT=status)
          CALL ERRMSG(substr,status,2)
          NULLIFY(varm%att)
          varm%natts = 0
       END IF
       CALL ADD_NCATT(varm, 'RG_MID_VALUES_OF', replace=.true. &
                      ,vs=TRIM(str))
    ELSE    ! m -> i
       str = TRIM(varm%name)
       vari%name  = TRIM(str)//'_I'
       vari%id    = NULL_VARID
       vtype = vari%dat%type
       SELECT CASE(vtype)
       CASE(VTYPE_REAL)
          vari%xtype = NF90_FLOAT
       CASE(VTYPE_DOUBLE)
          vari%xtype = NF90_DOUBLE
       CASE(VTYPE_INT)
          vari%xtype = NF90_INT
       CASE(VTYPE_BYTE)
          vari%xtype = NF90_BYTE
       CASE(VTYPE_CHAR)
          vari%xtype = NF90_CHAR
       CASE(VTYPE_UNDEF)
          CALL RGMSG(substr, RGMLE, &
               'TYPE OF VARIABLE '''//TRIM(vari%name)//''' IS UNDEFINED !')
       CASE DEFAULT
          CALL RGMSG(substr, RGMLE, &
               'TYPE OF VARIABLE '''//TRIM(vari%name)//''' IS UNRECOGNIZED !')
       END SELECT
       IF (varm%ndims > 0) THEN
          vari%ndims = varm%ndims
          ALLOCATE(vari%dim(vari%ndims), STAT=status)
          CALL ERRMSG(substr,status,3)
          DO i=1, vari%ndims
             vari%dim(i)%name  = TRIM(varm%dim(i)%name)//'_I'
             vari%dim(i)%id    = NULL_DIMID
             vari%dim(i)%len   = varm%dim(i)%len + 1
             vari%dim(i)%fuid  = varm%dim(i)%fuid
             vari%dim(i)%varid = NULL_VARID
          END DO
       END IF
       IF (vari%natts > 0) THEN
          DO i=1, vari%natts
             CALL INIT_NCATT(vari%att(i))
          END DO
          DEALLOCATE(vari%att, STAT=status)
          CALL ERRMSG(substr,status,4)
          NULLIFY(vari%att)
          vari%natts = 0
       END IF
       CALL ADD_NCATT(vari, 'RG_INTERFACE_VALUES_OF', replace=.true. &
                      ,vs=TRIM(str))
    END IF  ! i <-> m

  END SUBROUTINE IMMI_NCVAR
! --------------------------------------------
! --------------------------------------------
SUBROUTINE IMMI_NARRAY(nai, nam, ltestonly, lai, lam)

  IMPLICIT NONE

  ! I/O
  TYPE (t_narray),   INTENT(INOUT) :: nai, nam
  LOGICAL, OPTIONAL, INTENT(IN)    :: ltestonly ! do not calculate 
                                      ! interfaces/mids from each other (does 
                                      ! not work for curvilinear grids)
  LOGICAL, OPTIONAL, INTENT(OUT)   :: lai, lam ! interface / mids defined
  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'IMMI_NARRAY'
  INTEGER                     :: itype, mtype
  INTEGER                     :: i
  LOGICAL                     :: lonlytest

  IF (PRESENT(ltestonly)) THEN
     lonlytest = ltestonly
  ELSE
     lonlytest = .FALSE.
  ENDIF
  
  itype = nai%type
  mtype = nam%type

  IF (PRESENT(lai)) lai = (itype /= VTYPE_UNDEF)
  IF (PRESENT(lam)) lam = (mtype /= VTYPE_UNDEF)

  ! BOTH UNDEFINED
  IF ((itype == VTYPE_UNDEF).AND.(mtype == VTYPE_UNDEF)) THEN
     CALL RGMSG(substr, RGMLIC, '  ... BOTH UNDEFINED ...')
     CALL RGMSG(substr, RGMLIC, '  ... O.K.')
     RETURN
  END IF

  ! BOTH DEFINED: CHECK FOR COMPATIBILITY
  IF ((itype /= VTYPE_UNDEF).AND.(mtype /= VTYPE_UNDEF)) THEN
     CALL RGMSG(substr, RGMLIC, '  ... BOTH DEFINED ... ')
     IF (itype /= mtype) THEN
        CALL RGMSG(substr, RGMLE, &
             'INTERFACES AND MIDS ARE OF DIFFERENT TYPE !')
     END IF
     IF (nai%n /= 1) THEN
        CALL RGMSG(substr, RGMLIC, &
             'INTERFACE ARRAY IS ',nai%n,'-DIMENSIONAL !')
     END IF
     IF (nam%n /= 1) THEN
        CALL RGMSG(substr, RGMLIC, &
             'MID ARRAY IS ',nam%n,'-DIMENSIONAL !')
     END IF
     IF (nai%n /= nam%n) THEN
        CALL RGMSG(substr, RGMLE, &
             'DIMENSIONS OF INTERFACE ARRAY AND MID ARRAY ARE INCOMPATIBLE !')
     END IF

     DO i=1, nai%n     
        IF (nai%dim(i) /= (nam%dim(i)+1)) THEN
           CALL RGMSG(substr, RGMLE, &
                'SIZES OF INTERFACE ARRAY AND MID ARRAY ARE INCOMPATIBLE !')
        END IF
     END DO
     CALL RGMSG(substr, RGMLIC, '  ... O.K.')
     RETURN
  END IF

  IF (lonlytest) RETURN

  SELECT CASE(itype)
  CASE(VTYPE_REAL)
     ! nai -> nam
     CALL RGMSG(substr, RGMLIC, '  ... INTERFACES -> MIDS')
     IF (nai%n /= 1) THEN
        CALL RGMSG(substr, RGMLE, &
             'INTERFACE ARRAY IS ',nai%n,'-DIMENSIONAL !')
     END IF
     nai%dim(1) = nai%dim(1)-1
     CALL INIT_NARRAY( nam, 1, (/ nai%dim(1) /), VTYPE_REAL)
     nai%dim(1) = nai%dim(1)+1
     DO i=1, nam%dim(1)
        nam%vr(i) = (nai%vr(i)+nai%vr(i+1))/2.
     END DO
     CALL RGMSG(substr, RGMLIC, '  ... O.K.')
     IF (PRESENT(lam)) lam = .TRUE.
     RETURN
  CASE(VTYPE_DOUBLE)
     ! nai -> nam
     CALL RGMSG(substr, RGMLIC, '  ... INTERFACES -> MIDS')
     IF (nai%n /= 1) THEN
        CALL RGMSG(substr, RGMLE, &
             'INTERFACE ARRAY IS ',nai%n,'-DIMENSIONAL !')
     END IF
     nai%dim(1) = nai%dim(1)-1
     CALL INIT_NARRAY( nam, 1, (/ nai%dim(1) /), VTYPE_DOUBLE )
     nai%dim(1) = nai%dim(1)+1
     DO i=1, nam%dim(1)
        nam%vd(i) = (nai%vd(i)+nai%vd(i+1))/2.
     END DO
     CALL RGMSG(substr, RGMLIC, '  ... O.K.')
     IF (PRESENT(lam)) lam = .TRUE.
     RETURN
  CASE(VTYPE_INT)
     CALL RGMSG(substr, RGMLIC, '  ... INTERFACES -> MIDS')
     CALL RGMSG(substr, RGMLW, &
          'INTERFACES OF TYPE INTEGER CONVERTED TO DOUBLE !')
     IF (nai%n /= 1) THEN
        CALL RGMSG(substr, RGMLE, &
             'INTERFACE ARRAY IS ',nai%n,'-DIMENSIONAL !')
     END IF
     nai%dim(1) = nai%dim(1)-1
     CALL INIT_NARRAY( nam, 1, (/ nai%dim(1) /), VTYPE_DOUBLE )
     nai%dim(1) = nai%dim(1)+1
     DO i=1, nam%dim(1)
        nam%vd(i) = REAL(nai%vi(i) + nai%vi(i+1), DP)/REAL(2., DP)
     END DO
     CALL RGMSG(substr, RGMLIC, '  ... O.K.')
     IF (PRESENT(lam)) lam = .TRUE.
     RETURN
  CASE(VTYPE_BYTE)
     CALL RGMSG(substr, RGMLE, &
          'INTERFACES OF TYPE BYTE ARE NOT SUPPORTED !')
  CASE(VTYPE_CHAR)
     CALL RGMSG(substr, RGMLE, &
          'INTERFACES OF TYPE CHAR ARE NOT SUPPORTED !')
  CASE(VTYPE_UNDEF)
     ! OK, nam -> nai
  CASE DEFAULT
     CALL RGMSG(substr, RGMLE, &
          'UNRECOGNIZED TYPE OF INTERFACE ARRAY !')
  END SELECT

  SELECT CASE(mtype)
  CASE(VTYPE_REAL)
     ! nam -> nai
     CALL RGMSG(substr, RGMLIC, '  ... MIDS -> INTERFACES')
     IF (nam%n /= 1) THEN
        CALL RGMSG(substr, RGMLE, &
             'INTERFACE ARRAY IS ',nam%n,'-DIMENSIONAL !')
     END IF
     nam%dim(1) = nam%dim(1)+1
     CALL INIT_NARRAY( nai, 1, (/ nam%dim(1) /), VTYPE_REAL )
     nam%dim(1) = nam%dim(1)-1
     DO i=2, nai%dim(1)-1
        nai%vr(i) = (nam%vr(i-1)+nam%vr(i))/2.
     END DO
     nai%vr(1)          = 2.0*nam%vr(1) - nai%vr(2)
     nai%vr(nai%dim(1)) = 2.0*nam%vr(nai%dim(1)-1)-nai%vr(nai%dim(1)-1)
     CALL RGMSG(substr, RGMLIC, '  ... O.K.')
     IF (PRESENT(lai)) lai = .TRUE.
     RETURN
  CASE(VTYPE_DOUBLE)
     ! nam -> nai
     CALL RGMSG(substr, RGMLIC, '  ... MIDS -> INTERFACES')
     IF (nam%n /= 1) THEN
        CALL RGMSG(substr, RGMLE, &
             'INTERFACE ARRAY IS ',nam%n,'-DIMENSIONAL !')
     END IF
     nam%dim(1) = nam%dim(1)+1
     CALL INIT_NARRAY( nai, 1, (/ nam%dim(1) /), VTYPE_DOUBLE )
     nam%dim(1) = nam%dim(1)-1
     DO i=2, nai%dim(1)-1
        nai%vd(i) = (nam%vd(i-1)+nam%vd(i))/2.
     END DO
     nai%vd(1)          = 2.0*nam%vd(1) - nai%vd(2)
     nai%vd(nai%dim(1)) = 2.0*nam%vd(nai%dim(1)-1)-nai%vd(nai%dim(1)-1)
     CALL RGMSG(substr, RGMLIC, '  ... O.K.')
     IF (PRESENT(lai)) lai = .TRUE.
     RETURN
  CASE(VTYPE_INT)
     CALL RGMSG(substr, RGMLW, &
          'MIDS OF TYPE INTEGER CONVERTED TO DOUBLE !')
     ! nam -> nai
     CALL RGMSG(substr, RGMLIC, '  ... MIDS -> INTERFACES')
     IF (nam%n /= 1) THEN
        CALL RGMSG(substr, RGMLE, &
             'INTERFACE ARRAY IS ',nam%n,'-DIMENSIONAL !')
     END IF
     nam%dim(1) = nam%dim(1)+1
     CALL INIT_NARRAY( nai, 1, (/ nam%dim(1) /), VTYPE_DOUBLE )
     nam%dim(1) = nam%dim(1)-1
     DO i=2, nai%dim(1)-1
        nai%vd(i) = REAL(nam%vi(i-1)+nam%vi(i), DP)/REAL(2., DP)
     END DO
     nai%vd(1)          = REAL(2.0, DP) * (REAL(nam%vi(1), DP) - nai%vd(2))
     nai%vd(nai%dim(1)) = REAL(2.0, DP) * &
          (REAL(nam%vi(nai%dim(1)-1), DP)-nai%vd(nai%dim(1)-1))
     CALL RGMSG(substr, RGMLIC, '  ... O.K.')
     IF (PRESENT(lai)) lai = .TRUE.
     RETURN
  CASE(VTYPE_BYTE)
     CALL RGMSG(substr, RGMLE, &
          'MIDS OF TYPE BYTE ARE NOT SUPPORTED !')
  CASE(VTYPE_CHAR)
     CALL RGMSG(substr, RGMLE, &
          'MIDS OF TYPE CHAR ARE NOT SUPPORTED !')
  CASE(VTYPE_UNDEF)
     ! OK, nai -> nam  ! THIS POSITION SHOULD NEVER BE REACHED
     CALL RGMSG(substr, RGMLE, &
          'TYPE OF INTERFACE ARRAY IS UNDEFINED !')
  CASE DEFAULT
     CALL RGMSG(substr, RGMLE, &
          'UNRECOGNIZED TYPE OF INTERFACE ARRAY !')
  END SELECT

END SUBROUTINE IMMI_NARRAY
! --------------------------------------------
! --------------------------------------------
  SUBROUTINE RNGADJ_NARRAY(na, mr, lmon)

    IMPLICIT NONE

    ! I/O
    TYPE (t_narray), INTENT(INOUT)        :: na
    REAL(DP)  ,      INTENT(IN)           :: mr(2)
    LOGICAL,         INTENT(IN), OPTIONAL :: lmon

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'RNGADJ_NARRAY'
    INTEGER                     :: n
    INTEGER                     :: vtype
    LOGICAL                     :: llmon

    IF (PRESENT(lmon)) THEN
       llmon = lmon
    ELSE
       llmon = .true.   ! DEFAULT: COORDINATE IS MONOTONIC
    END IF

    vtype = na%type

    IF (llmon) THEN      ! MONOTONIC COORDINATE

       SELECT CASE(vtype)
       CASE(VTYPE_REAL)
          n = SIZE(na%vr)
          IF (na%vr(1) <= na%vr(n)) THEN
             IF (ABS(MINVAL(mr) - RGEMPTY) >= TINY(RGEMPTY)) &
                  na%vr(1) = REAL(MINVAL(mr), SP)
             IF (ABS(MAXVAL(mr) - RGEMPTY) >= TINY(RGEMPTY)) &
                  na%vr(n) = REAL(MAXVAL(mr), SP)
          ELSE
             IF (ABS(MINVAL(mr) - RGEMPTY) >= TINY(RGEMPTY)) &
                  na%vr(n) = REAL(MINVAL(mr), SP)
             IF (ABS(MAXVAL(mr) - RGEMPTY) >= TINY(RGEMPTY)) &
                  na%vr(1) = REAL(MAXVAL(mr), SP)
          END IF
       CASE(VTYPE_DOUBLE)
          n = SIZE(na%vd)
          IF (na%vd(1) <= na%vd(n)) THEN
             IF (ABS(MINVAL(mr) - RGEMPTY) >= TINY(RGEMPTY)) &
                  na%vd(1) = REAL(MINVAL(mr), DP)
             IF (ABS(MAXVAL(mr) - RGEMPTY) >= TINY(RGEMPTY)) &
                  na%vd(n) = REAL(MAXVAL(mr), DP)
          ELSE
             IF (ABS(MINVAL(mr) - RGEMPTY) >= TINY(RGEMPTY)) &
                  na%vd(n) = REAL(MINVAL(mr), DP)
             IF (ABS(MAXVAL(mr) - RGEMPTY) >= TINY(RGEMPTY)) &
                  na%vd(1) = REAL(MAXVAL(mr), DP)
          END IF
       CASE(VTYPE_INT)
          n = SIZE(na%vi)
          IF (na%vi(1) <= na%vi(n)) THEN
             IF (ABS(MINVAL(mr) - RGEMPTY) >= TINY(RGEMPTY)) &
                  na%vi(1) = INT(MINVAL(mr), I8)
             IF (ABS(MAXVAL(mr) - RGEMPTY) >= TINY(RGEMPTY)) &
                  na%vi(n) = INT(MAXVAL(mr), I8)
          ELSE
             IF (ABS(MINVAL(mr) - RGEMPTY) >= TINY(RGEMPTY)) &
                  na%vi(n) = INT(MINVAL(mr), I8)
             IF (ABS(MAXVAL(mr) - RGEMPTY) >= TINY(RGEMPTY)) &
                  na%vi(1) = INT(MAXVAL(mr), I8)
          END IF
       CASE(VTYPE_BYTE)
          n = SIZE(na%vb)
          IF (na%vb(1) <= na%vb(n)) THEN
             IF (ABS(MINVAL(mr) - RGEMPTY) >= TINY(RGEMPTY)) &
                  na%vb(1) = INT(MINVAL(mr), I4)
             IF (ABS(MAXVAL(mr) - RGEMPTY) >= TINY(RGEMPTY)) &
                  na%vb(n) = INT(MAXVAL(mr), I4)
          ELSE
             IF (ABS(MINVAL(mr) - RGEMPTY) >= TINY(RGEMPTY)) &
                  na%vb(n) = INT(MINVAL(mr), I4)
             IF (ABS(MAXVAL(mr) - RGEMPTY) >= TINY(RGEMPTY)) &
                  na%vb(1) = INT(MAXVAL(mr), I4)
          END IF
       CASE(VTYPE_CHAR)
          CALL RGMSG(substr, RGMLE, 'N-ARRAY IS OF TYPE CHAR !')
       CASE(VTYPE_UNDEF)
          CALL RGMSG(substr, RGMLW, 'TYPE OF N-ARRAY IS UNDEFINED !')
       CASE DEFAULT
          CALL RGMSG(substr, RGMLE, 'TYPE OF N-ARRAY IS NOT RECOGNIZED !')
       END SELECT

    ELSE   ! NON-MONOTONIC COORDINATE

       SELECT CASE(vtype)
       CASE(VTYPE_REAL)
          n = SIZE(na%vr)
          IF (ABS(mr(1) - RGEMPTY) >= TINY(RGEMPTY)) &
               na%vr(1) = REAL(mr(1), SP)
          IF (ABS(mr(2) - RGEMPTY) >= TINY(RGEMPTY)) &
               na%vr(n) = REAL(mr(2), SP)
       CASE(VTYPE_DOUBLE)
          n = SIZE(na%vd)
          IF (ABS(mr(1) - RGEMPTY)>= TINY(RGEMPTY)) &
               na%vd(1) = REAL(mr(1), DP)
          IF (ABS(mr(2) - RGEMPTY)>= TINY(RGEMPTY)) &
               na%vd(n) = REAL(mr(2), DP)
       CASE(VTYPE_INT)
          n = SIZE(na%vi)
          IF (ABS(mr(1) - RGEMPTY)>= TINY(RGEMPTY)) &
               na%vi(1) = INT(mr(1), I8)
          IF (ABS(mr(2) - RGEMPTY)>= TINY(RGEMPTY)) &
               na%vi(n) = INT(mr(2), I8)
       CASE(VTYPE_BYTE)
          n = SIZE(na%vb)
          IF (ABS(mr(1) - RGEMPTY)>= TINY(RGEMPTY)) &
               na%vb(1) = INT(mr(1), I4)
          IF (ABS(mr(2) - RGEMPTY)>= TINY(RGEMPTY)) &
               na%vb(n) = INT(mr(2), I4)
       CASE(VTYPE_CHAR)
          CALL RGMSG(substr, RGMLE, 'N-ARRAY IS OF TYPE CHAR !')
       CASE(VTYPE_UNDEF)
          CALL RGMSG(substr, RGMLW, 'TYPE OF N-ARRAY IS UNDEFINED !')
       CASE DEFAULT
          CALL RGMSG(substr, RGMLE, 'TYPE OF N-ARRAY IS NOT RECOGNIZED !')
       END SELECT

    END IF ! MONOTONIC COORDINATE ?

  END SUBROUTINE RNGADJ_NARRAY
! --------------------------------------------

! --------------------------------------------
  SUBROUTINE RNGADJ_ARRAY(arr, mr)

    ! APPLY RANGES, ASSUME MONOTONIC AXIS

    IMPLICIT NONE

    ! I/O
    REAL(DP), DIMENSION(:), INTENT(INOUT)        :: arr
    REAL(DP)  ,             INTENT(IN)           :: mr(2)

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'RNGADJ_ARRAY'
    INTEGER                     :: n

    n = SIZE(arr)
    IF (arr(1) <= arr(n)) THEN
       IF (ABS(MINVAL(mr) - RGEMPTY) >= TINY(RGEMPTY)) &
            arr(1) = REAL(MINVAL(mr), DP)
       IF (ABS(MAXVAL(mr) - RGEMPTY) >= TINY(RGEMPTY)) THEN
          IF (ABS(MAXVAL(mr) - RGMAX) >= TINY(RGMAX)) THEN
             arr(n) = REAL(MAXVAL(mr), DP)
          ELSE
             ! ENSURE SURFACE PRESSURE
             !arr(n) = MAX(arr(n), 101325._dp)
          ENDIF
       END IF
    ELSE
       IF (ABS(MINVAL(mr) - RGEMPTY) >= TINY(RGEMPTY)) &
            arr(n) = REAL(MINVAL(mr), DP)
       IF (ABS(MAXVAL(mr) - RGEMPTY) >= TINY(RGEMPTY)) THEN
            IF (ABS(MAXVAL(mr) - RGMAX) >= TINY(RGMAX)) THEN
               arr(1) = REAL(MAXVAL(mr), DP)
            ELSE
             ! ENSURE SURFACE PRESSURE
             !arr(1) = MAX(arr(1), 101325._dp)
          ENDIF
               
         END IF
    END IF

  END SUBROUTINE RNGADJ_ARRAY
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE CHECK_NCVAR_ON_GEOHYBGRID(var, g, dims, axes, ok)

  IMPLICIT NONE

  ! I/O
  TYPE (t_ncvar)     , INTENT(IN)  :: var    ! variable
  TYPE (t_geohybgrid), INTENT(IN)  :: g      ! grid
  INTEGER, DIMENSION(:), POINTER   :: dims   ! order of g-dims in var
  INTEGER            , INTENT(OUT) :: axes(3)! dimension no. of lon->lat->lev
  LOGICAL            , INTENT(OUT) :: ok     ! conform ?

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'CHECK_NCVAR_ON_GEOHYBGRID'
  INTEGER                     :: i
  INTEGER                     :: ndimg            ! number of dimensions in g
  INTEGER                     :: ndimv            ! number of g-dims in var
  INTEGER                     :: status
  ! position of lon, lat, lev in AXES-LIST !!!
  INTEGER                     :: plon, plat, plev 
  LOGICAL                     :: lvert            ! vertical axis present ?
  LOGICAL                     :: llat, llon       ! lat, lon axes present ?

  ! INIT
  ok = .true.
  axes(:) = 0
  plon = 0
  plat = 0
  plev = 0

  ! ALLOCATE SPACE FOR 'ORDER OF DIMENSIONS'
  IF (ASSOCIATED(dims)) THEN
     DEALLOCATE(dims)
     NULLIFY(dims)
  END IF

  ALLOCATE(dims(var%ndims), STAT=status)
  CALL ERRMSG(substr,status,1)
  dims(:) = 0

  ! GET NUMBER OF HYBRID-DIMENSIONS IN GRID
  ndimg = 0

  IF (QDEF_NCVAR(g%lonm).OR.QDEF_NCVAR(g%clonm).OR.QDEF_NCVAR(g%ulonm)) THEN
     ndimg = ndimg + 1
     plon = ndimg
  END IF
  IF (QDEF_NCVAR(g%latm).OR.QDEF_NCVAR(g%clatm).OR.QDEF_NCVAR(g%ulatm)) THEN
     ndimg = ndimg + 1
     plat = ndimg
  END IF
  IF (QDEF_NCVAR(g%hyam).OR.QDEF_NCVAR(g%hybm)) THEN
     ndimg = ndimg + 1
     plev = ndimg
  END IF

  ! GET GRID-DIMENSION POSITIONS IN VARIABLE
  ndimv = 0
  lvert = .false.
  llat = .false.
  llon = .false.

  DO i=1, var%ndims  ! LOOP OVER VARIABLE DIMENSIONS

     ! CHECK DIM LENGTH AND NAME/ID
     IF (.NOT.llon) THEN  ! CHECK LON
        IF (ASSOCIATED(g%lonm%dim)) THEN
           llon = (QCMP_NCDIM(g%lonm%dim(1), var%dim(i)) > 1)
           IF (ASSOCIATED(g%loni%dim).AND.(.NOT.llon)) THEN
              llon = ((QCMP_NCDIM(g%lonm%dim(1), var%dim(i)) > 0).AND. &
                      (TRIM(g%lonm%dim(1)%name) ==                     &
                       TRIM(g%loni%dim(1)%name)//'_M') )
           END IF
        ENDIF
        IF (ASSOCIATED(g%clonm%dim) .AND. .NOT. llon) THEN
           llon = (QCMP_NCDIM(g%clonm%dim(1), var%dim(i)) > 1)
           IF (ASSOCIATED(g%cloni%dim).AND.(.NOT.llon)) THEN
              llon = ((QCMP_NCDIM(g%clonm%dim(1), var%dim(i)) > 0).AND. &
                   (TRIM(g%clonm%dim(1)%name) ==                     &
                   TRIM(g%cloni%dim(1)%name)//'_M') )
           END IF
        ENDIF

        IF (ASSOCIATED(g%ulonm%dim) .AND. .NOT. llon) THEN
           llon = (QCMP_NCDIM(g%ulonm%dim(1), var%dim(i)) > 1)
        ENDIF

        IF (llon) THEN
           ndimv = ndimv + 1
           dims(i) = plon
           axes(GORD_LON) = i
           CYCLE
        END IF
     END IF

     IF (.NOT.llat) THEN  ! CHECK LAT
        IF (ASSOCIATED(g%latm%dim)) THEN
           IF (g%latm%ndims == 1) THEN
              llat = (QCMP_NCDIM(g%latm%dim(1), var%dim(i)) > 1)
              IF (ASSOCIATED(g%lati%dim).AND.(.NOT.llat)) THEN
                 llat = ((QCMP_NCDIM(g%latm%dim(1), var%dim(i)) > 0).AND. &
                      (TRIM(g%latm%dim(1)%name) ==                     &
                      TRIM(g%lati%dim(1)%name)//'_M') )
              END IF
           ELSE
              llat = (QCMP_NCDIM(g%latm%dim(2), var%dim(i)) > 1)
              IF (ASSOCIATED(g%lati%dim).AND.(.NOT.llat)) THEN
                 llat = ((QCMP_NCDIM(g%latm%dim(2), var%dim(i)) > 0).AND. &
                      (TRIM(g%latm%dim(2)%name) ==                     &
                      TRIM(g%lati%dim(2)%name)//'_M') )
              END IF
           END IF
        END IF

        IF (ASSOCIATED(g%clatm%dim) .AND. .NOT. llat) THEN
           IF (g%clatm%ndims == 1) THEN
              llat = (QCMP_NCDIM(g%clatm%dim(1), var%dim(i)) > 1)
              IF (ASSOCIATED(g%clati%dim).AND.(.NOT.llat)) THEN
                 llat = ((QCMP_NCDIM(g%clatm%dim(1), var%dim(i)) > 0).AND. &
                      (TRIM(g%clatm%dim(1)%name) ==                     &
                      TRIM(g%clati%dim(1)%name)//'_M') )
              END IF
           ELSE
              llat = (QCMP_NCDIM(g%clatm%dim(2), var%dim(i)) > 1)
              IF (ASSOCIATED(g%clati%dim).AND.(.NOT.llat)) THEN
                 llat = ((QCMP_NCDIM(g%clatm%dim(2), var%dim(i)) > 0).AND. &
                      (TRIM(g%clatm%dim(2)%name) ==                     &
                      TRIM(g%clati%dim(2)%name)//'_M') )
              END IF
           END IF
        END IF
        IF (ASSOCIATED(g%ulatm%dim) .AND. .NOT. llat) THEN
           IF (g%ulatm%ndims == 1) THEN
              llat = (QCMP_NCDIM(g%ulatm%dim(1), var%dim(i)) > 1)
           ELSE
              llat = (QCMP_NCDIM(g%ulatm%dim(2), var%dim(i)) > 1)
           END IF
        END IF
        IF (llat) THEN
           ndimv = ndimv + 1
           dims(i) = plat
           axes(GORD_LAT) = i
           CYCLE
        END IF
     END IF ! llat

     IF (.NOT.lvert) THEN  ! CHECK HYA
        IF (ASSOCIATED(g%hyam%dim)) THEN
           lvert = (QCMP_NCDIM(g%hyam%dim(1), var%dim(i)) > 1)
           ! COMPLETE GEOHYBGRID WAS CALLED BEFORE, THUS ALLWAYS HYAM AND
           ! HYAI are associated
           IF (ASSOCIATED(g%hyai%dim).AND.(.NOT.lvert)) THEN
              lvert = (QCMP_NCDIM(g%hyai%dim(1), var%dim(i)) > 1) .OR. &
                  ((QCMP_NCDIM(g%hyam%dim(1), var%dim(i)) > 0).AND. &
                  (TRIM(g%hyam%dim(1)%name) == TRIM(g%hyai%dim(1)%name)//'_M'))
           END IF
           IF (lvert) THEN
              ndimv = ndimv + 1
              dims(i) = plev
              axes(GORD_LEV) = i
              CYCLE
           END IF
        END IF
     END IF

     IF (.NOT.lvert) THEN  ! CHECK HYB
        IF (ASSOCIATED(g%hybm%dim)) THEN
             lvert = (QCMP_NCDIM(g%hybm%dim(1), var%dim(i)) > 1)
             ! COMPLETE GEOHYBGRID WAS CALLED BEFORE, THUS ALLWAYS HYAM AND
             ! HYAI are associated
             IF (ASSOCIATED(g%hybi%dim).AND.(.NOT.lvert)) THEN
                lvert = (QCMP_NCDIM(g%hybi%dim(1), var%dim(i)) > 1) .OR. &
                   ((QCMP_NCDIM(g%hybm%dim(1), var%dim(i)) > 0).AND.     &
                   (TRIM(g%hybm%dim(1)%name) == TRIM(g%hybi%dim(1)%name)//'_M'))
             END IF
             IF (lvert) THEN
                ndimv = ndimv + 1
                dims(i) = plev
                axes(GORD_LEV) = i
                CYCLE
             END IF
        END IF
     END IF

  END DO  ! LOOP OVER DIMENSIONS

!  ok = (ndimv >= ndimg)  ! ALL GRID DIMS MUST BE RECOGNIZED !!!
                          ! REMAINING DIMS ARE 'FREE'
  ok = (ndimv >= ndimg) .OR. ((.NOT. lvert) .AND. (ndimv + 1 >= ndimg))

  ! ALL GRID DIMS MUST BE RECOGNIZED !!!
  ! REMAINING DIMS ARE 'FREE'
  ! 2D fields on 3D grid also possible

END SUBROUTINE CHECK_NCVAR_ON_GEOHYBGRID
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE SORT_GEOHYBGRID_NCVAR(var, gx, axes, svar, reverse)

  IMPLICIT NONE

  ! I/O
  TYPE (t_ncvar)     , INTENT(IN)           :: var     ! input variable
  TYPE (t_geohybgrid), INTENT(IN)           :: gx      ! hybrid grid with  ...
                                                       ! ... index information
  INTEGER            , INTENT(IN)           :: axes(3) ! lon,lat,lev dim. no.
  TYPE (t_ncvar)     , INTENT(INOUT)        :: svar    ! sorted variable
  LOGICAL            , INTENT(IN), OPTIONAL :: reverse

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER        :: substr = 'SORT_GEOHYBGRID_NCVAR'
  INTEGER                            :: i
  INTEGER                            :: vtype
  INTEGER, DIMENSION(:), ALLOCATABLE :: vdim   ! variable dimension vector
  INTEGER                            :: status
  INTEGER, DIMENSION(:), POINTER     :: vec    ! element vector
  LOGICAL                            :: lrev   ! local reverse flag
  
  ! INIT
  IF (PRESENT(reverse)) THEN
     lrev = reverse
  ELSE
     lrev = .false.   ! DEFAULT
  END IF
  !
  NULLIFY(vec)
  !
  CALL COPY_NCVAR(svar, var)
  !
  ALLOCATE(vdim(var%ndims), STAT=status)
  CALL ERRMSG(substr,status,1)

  DO i=1, var%ndims
     vdim(i) = var%dim(i)%len
  END DO

  vtype = var%dat%type

  DO i=1, PRODUCT(vdim)  ! LOOP OVER ALL ELEMENTS
     ! GET ELEMENT VECTOR
     CALL ELEMENT(vdim,i,vec)
     ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
     IF (axes(GORD_LON) > 0) THEN
        vec(axes(GORD_LON)) = INT(gx%lonm%dat%vi(vec(axes(GORD_LON))))
     END IF
     IF (axes(GORD_LAT) > 0) THEN
        vec(axes(GORD_LAT)) = INT(gx%latm%dat%vi(vec(axes(GORD_LAT))))
     END IF
     IF (axes(GORD_LEV) > 0) THEN
        IF (gx%hyam%dat%type /= VTYPE_UNDEF) THEN
           vec(axes(GORD_LEV)) = INT(gx%hyam%dat%vi(vec(axes(GORD_LEV))))
        ELSE
           vec(axes(GORD_LEV)) = INT(gx%hybm%dat%vi(vec(axes(GORD_LEV))))
        END IF
     END IF
     ! COPY DATA ELEMENT
     SELECT CASE(vtype)
     CASE(VTYPE_INT)
        IF (lrev) THEN   ! UNSORT
           svar%dat%vi(POSITION(vdim,vec)) = var%dat%vi(i)
        ELSE             ! SORT
           svar%dat%vi(i) = var%dat%vi(POSITION(vdim,vec))
        END IF
     CASE(VTYPE_REAL)
        IF (lrev) THEN
           svar%dat%vr(POSITION(vdim,vec)) = var%dat%vr(i)
        ELSE
           svar%dat%vr(i) = var%dat%vr(POSITION(vdim,vec))
        END IF
     CASE(VTYPE_DOUBLE)
        IF (lrev) THEN
           svar%dat%vd(POSITION(vdim,vec)) = var%dat%vd(i)
        ELSE
           svar%dat%vd(i) = var%dat%vd(POSITION(vdim,vec))
        END IF
     CASE(VTYPE_CHAR)
        IF (lrev) THEN
           svar%dat%vc(POSITION(vdim,vec)) = var%dat%vc(i)
        ELSE
           svar%dat%vc(i) = var%dat%vc(POSITION(vdim,vec))
        END IF
     CASE(VTYPE_BYTE)
        IF (lrev) THEN
           svar%dat%vb(POSITION(vdim,vec)) = var%dat%vb(i)
        ELSE
           svar%dat%vb(i) = var%dat%vb(POSITION(vdim,vec))
        END IF
     CASE(VTYPE_UNDEF)
        CALL RGMSG(substr, RGMLE, 'UNDEFINED VARIABLE CANNOT BE SORTED !')
     CASE DEFAULT
        CALL RGMSG(substr, RGMLE, 'UNRECOGNIZED TYPE OF VARIABLE !')
     END SELECT

     DEALLOCATE(vec, STAT=status)
     CALL ERRMSG(substr,status,2)
     NULLIFY(vec)
  END DO

  ! CLEAN UP
  DEALLOCATE(vdim, STAT=status)
  CALL ERRMSG(substr,status,3)

END SUBROUTINE SORT_GEOHYBGRID_NCVAR
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE PACK_GEOHYBGRID_NCVAR(vi, dims, axes ,vo, reverse)

  IMPLICIT NONE

  ! I/O
  TYPE (t_ncvar)                , INTENT(IN)    :: vi   ! input variable
  INTEGER, DIMENSION(:)         , INTENT(IN)    :: dims
  INTEGER                       , INTENT(IN)    :: axes(3)
  TYPE (t_ncvar)                , INTENT(INOUT) :: vo   ! output variable
  LOGICAL, OPTIONAL             , INTENT(IN)    :: reverse

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER        :: substr = 'PACK_GEOHYBGRID_NCVAR'
  TYPE (t_ncvar)                     :: vu     ! unpacked variable
  TYPE (t_ncvar)                     :: vp     ! packed variable
  INTEGER                            :: nfree  ! length of free dimension
  INTEGER                            :: npdim  ! number of dims in packed
  INTEGER, DIMENSION(:), ALLOCATABLE :: pdim   ! packed dimension vector
  INTEGER, DIMENSION(:), ALLOCATABLE :: pvec   ! packed position vector
  INTEGER, DIMENSION(:), ALLOCATABLE :: umdim  ! main unpacked dim. vector
  INTEGER, DIMENSION(:), ALLOCATABLE :: urdim  ! rest unpacked dim. vector
  INTEGER, DIMENSION(:), POINTER     :: umvec  ! unpacked main pos. vector
  INTEGER, DIMENSION(:), ALLOCATABLE :: urvec  ! unpacked rest pos. vector
  INTEGER, DIMENSION(:), ALLOCATABLE :: map    ! map indices u on p
  INTEGER, DIMENSION(:), ALLOCATABLE :: rmap   ! map invar. indices u on p
  INTEGER                            :: status
  INTEGER                            :: i, j
  INTEGER                            :: vtype
  LOGICAL                            :: lrev

  IF (PRESENT(reverse)) THEN
     lrev = reverse
  ELSE
     lrev = .false. ! DEFAULT
  END IF

  NULLIFY(umvec)

  ! INIT
  IF (lrev) THEN
     CALL COPY_NCVAR(vu, vo) ! output is unpacked (must be balanced on input!)
     CALL COPY_NCVAR(vp, vi) ! input is packed
  ELSE
     CALL COPY_NCVAR(vu, vi) ! input is unpacked
     CALL INIT_NCVAR(vp)     ! output is packed (new)
  END IF

  ! CHECK
  IF (SIZE(dims) /= vu%ndims) THEN
     CALL RGMSG(substr, RGMLE, &
          'DIMENSION MISMATCH BETWEEN NUMBER OF DIMENSIONS IN VARIABLE', &
          .false.)
     CALL RGMSG(substr, RGMLEC, &
          ''''//TRIM(vu%name)//''' (',vu%ndims,')', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'AND LENGTH OF DIMESNSION IN ORDER VECTOR (',SIZE(dims),') !')
  END IF

  ! GET NUMBER OF DIMS
  npdim = MAXVAL(dims) + 1      ! plus one 'free' dimension

  ! CALCULATE LENGTH OF FREE DIMENSION
  nfree = 1
  DO i=1, vu%ndims
     IF (dims(i) == 0) nfree = nfree * vu%dim(i)%len
  END DO
  
  IF (.NOT.lrev) THEN   ! packing mode
     ! SPACE FOR DIMENSIONS
     ALLOCATE(vp%dim(npdim),STAT=status)
     CALL ERRMSG(substr,status,1)
     ! COPY DIMS
     DO i=1,npdim-1
        vp%dim(i) = vu%dim(axes(i))
     END DO
     CALL INIT_NCDIM(vp%dim(npdim))
     vp%dim(npdim)%len  = nfree
     vp%dim(npdim)%name = TRIM(vu%name)//'_pd'
     !
     ! NCVAR
     vp%name  = TRIM(vu%name)//'_p'
     vp%xtype = vu%xtype
     vp%ndims = npdim
     vp%uid   = vu%uid
     vp%natts = vu%natts
     ALLOCATE(vp%att(vp%natts),STAT=status)
     CALL ERRMSG(substr,status,2)
     DO i=1, vp%natts
        CALL COPY_NCATT(vp%att(i), vu%att(i))
     END DO
     CALL ADD_NCATT(vp, 'RG_PACKED_VAR', vs=TRIM(vu%name))
  END IF  ! packing mode

  ! SPACE FOR DIMENSION/POSITION VECTORS
  ALLOCATE(pdim(npdim), STAT=status)
  CALL ERRMSG(substr,status,3)
  pdim(:) = 0
  ALLOCATE(pvec(npdim), STAT=status)
  CALL ERRMSG(substr,status,4)
  pvec(:) = 0
  ALLOCATE(map(npdim), STAT=status)
  CALL ERRMSG(substr,status,5)
  map(:) = 0
  !
  IF (vu%ndims >= npdim) THEN
     ALLOCATE(urdim(vu%ndims-npdim+1), STAT=status)
     CALL ERRMSG(substr,status,6)
     urdim(:) = 0
     !
     ALLOCATE(urvec(vu%ndims-npdim+1), STAT=status)
     CALL ERRMSG(substr,status,7)
     urvec(:) = 0
     !
     ALLOCATE(rmap(vu%ndims-npdim+1), STAT=status)
     CALL ERRMSG(substr,status,8)
     rmap(:) = 0
  END IF
  !
  ALLOCATE(umdim(vu%ndims), STAT=status)
  CALL ERRMSG(substr,status,9)
  umdim(:) = 0

  ! CALCULATE MAPPING OF INDICES, DIMENSION VECTORS
  j = 0
  DO i=1, vu%ndims  ! LOOP OVER VARIABLE DIMENSIONS
     umdim(i) = vu%dim(i)%len
     IF (dims(i) > 0) THEN
        pdim(dims(i)) = vu%dim(i)%len
        map(dims(i)) = i
     ELSE
        j = j+1
        urdim(j) = vu%dim(i)%len
        rmap(j)  = i
     END IF
  END DO
  pdim(npdim) = nfree

  ! ALLOCATE SPACE FOR PACKED DATA
  IF (lrev) THEN
     vtype = vp%dat%type
     ! NOTE: SPACE FOR 'UNPACKED' DATA
     !       WAS ALREADY ALLOCATED IN BALANCE_GEOHYBGRID_NCVAR
  ELSE
     vtype = vu%dat%type
     CALL INIT_NARRAY(vp%dat, npdim, pdim, vtype)
  END IF

  ! COPY DATA
  DO i=1, PRODUCT(pdim)  ! == SIZE(v?%dat%v?) ! LOOP OVER ALL ELEMENTS
     ! ELEMENT VECTOR IN UNPACKED SPACE
     CALL ELEMENT(umdim, i, umvec)
     ! ELEMENT VECTOR IN PACKED SPACE
     DO j=1, npdim-1
        pvec(j) = umvec(map(j))
     END DO
     IF (vu%ndims >= npdim) THEN
        DO j=1, vu%ndims-npdim+1
           urvec(j) = umvec(rmap(j))
        END DO
        pvec(npdim) = POSITION(urdim, urvec)
     ELSE
        pvec(npdim) = 1
     END IF
     j = POSITION(pdim, pvec)
     !
     SELECT CASE(vtype)
     CASE(VTYPE_INT)
        IF (lrev) THEN
           vu%dat%vi(i) = vp%dat%vi(j)
        ELSE
           vp%dat%vi(j) = vu%dat%vi(i)
        END IF
     CASE(VTYPE_REAL)
        IF (lrev) THEN
           vu%dat%vr(i) = vp%dat%vr(j)
        ELSE
           vp%dat%vr(j) = vu%dat%vr(i)
        END IF
     CASE(VTYPE_DOUBLE)
        IF (lrev) THEN
           vu%dat%vd(i) = vp%dat%vd(j)
        ELSE
           vp%dat%vd(j) = vu%dat%vd(i)
        END IF
     CASE(VTYPE_CHAR)
        IF (lrev) THEN
           vu%dat%vc(i) = vp%dat%vc(j)
        ELSE
           vp%dat%vc(j) = vu%dat%vc(i)
        END IF
     CASE(VTYPE_BYTE)
        IF (lrev) THEN
           vu%dat%vb(i) = vp%dat%vb(j)
        ELSE
           vp%dat%vb(j) = vu%dat%vb(i)
        END IF
     CASE(VTYPE_UNDEF)
        CALL RGMSG(substr, RGMLE, 'TYPE OF VARIABLE IS UNDEFINED !')
     CASE DEFAULT
        CALL RGMSG(substr, RGMLE, 'UNRECOGNIZED TYPE OF VARIABLE !')
     END SELECT
     !
     DEALLOCATE(umvec, STAT=status)
     CALL ERRMSG(substr,status,10)
     NULLIFY(umvec)
  END DO

  ! COPY RESULT TO OUTPUT
  CALL INIT_NCVAR(vo)
  IF (lrev) THEN
     CALL COPY_NCVAR(vo, vu)
  ELSE
     CALL COPY_NCVAR(vo, vp)
  END IF

  ! CLEAN UP
  CALL INIT_NCVAR(vu)
  CALL INIT_NCVAR(vp)
  !
  DEALLOCATE(pdim, STAT=status)
  CALL ERRMSG(substr,status,11)
  DEALLOCATE(pvec, STAT=status)
  CALL ERRMSG(substr,status,12)
  DEALLOCATE(map, STAT=status)
  CALL ERRMSG(substr,status,13)
  DEALLOCATE(umdim, STAT=status)
  CALL ERRMSG(substr,status,14)
  IF (vu%ndims >= npdim) THEN
     DEALLOCATE(rmap, STAT=status)
     CALL ERRMSG(substr,status,15)
     DEALLOCATE(urdim, STAT=status)
     CALL ERRMSG(substr,status,16)
     DEALLOCATE(urvec, STAT=status)
     CALL ERRMSG(substr,status,17)
  END IF

END SUBROUTINE PACK_GEOHYBGRID_NCVAR
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE BALANCE_GEOHYBGRID(gi, go)

  IMPLICIT NONE

  ! I/O
  TYPE (t_geohybgrid), INTENT(INOUT) :: gi      ! input grid
  TYPE (t_geohybgrid), INTENT(INOUT) :: go      ! output grid

  ! LOCAL
  !  CHARACTER(LEN=*), PARAMETER :: substr = 'BALANCE_GEOHYBGRID'

  IF (QDEF_NCVAR(gi%lonm).AND.(.NOT.QDEF_NCVAR(go%lonm))) THEN
     CALL COPY_NCVAR(go%lonm, gi%lonm)
  END IF

  IF (QDEF_NCVAR(gi%loni).AND.(.NOT.QDEF_NCVAR(go%loni))) THEN
     CALL COPY_NCVAR(go%loni, gi%loni)
  END IF

  IF (QDEF_NCVAR(gi%latm).AND.(.NOT.QDEF_NCVAR(go%latm))) THEN
     CALL COPY_NCVAR(go%latm, gi%latm)
  END IF

  IF (QDEF_NCVAR(gi%lati).AND.(.NOT.QDEF_NCVAR(go%lati))) THEN
     CALL COPY_NCVAR(go%lati, gi%lati)
  END IF

  IF (QDEF_NCVAR(gi%clonm).AND.(.NOT.QDEF_NCVAR(go%clonm))) THEN
     CALL COPY_NCVAR(go%clonm, gi%clonm)
  END IF

  IF (QDEF_NCVAR(gi%clatm).AND.(.NOT.QDEF_NCVAR(go%clatm))) THEN
     CALL COPY_NCVAR(go%clatm, gi%clatm)
  END IF

  IF (QDEF_NCVAR(gi%cloni).AND.(.NOT.QDEF_NCVAR(go%cloni))) THEN
     CALL COPY_NCVAR(go%cloni, gi%cloni)
  END IF

  IF (QDEF_NCVAR(gi%clati).AND.(.NOT.QDEF_NCVAR(go%clati))) THEN
     CALL COPY_NCVAR(go%clati, gi%clati)
  END IF
  IF (QDEF_NCVAR(gi%ulonm).AND.(.NOT.QDEF_NCVAR(go%ulonm))) THEN
     CALL COPY_NCVAR(go%ulonm, gi%ulonm)
  END IF

  IF (QDEF_NCVAR(gi%ulatm).AND.(.NOT.QDEF_NCVAR(go%ulatm))) THEN
     CALL COPY_NCVAR(go%ulatm, gi%ulatm)
  END IF

  IF (QDEF_NCVAR(gi%uloni).AND.(.NOT.QDEF_NCVAR(go%uloni))) THEN
     CALL COPY_NCVAR(go%uloni, gi%uloni)
  END IF

  IF (QDEF_NCVAR(gi%ulati).AND.(.NOT.QDEF_NCVAR(go%ulati))) THEN
     CALL COPY_NCVAR(go%ulati, gi%ulati)
  END IF

  IF (QDEF_NCVAR(gi%hyam).AND.(.NOT.QDEF_NCVAR(go%hyam))) THEN
     CALL COPY_NCVAR(go%hyam, gi%hyam)
  END IF

  IF (QDEF_NCVAR(gi%hyai).AND.(.NOT.QDEF_NCVAR(go%hyai))) THEN
     CALL COPY_NCVAR(go%hyai, gi%hyai)
  END IF

  IF (QDEF_NCVAR(gi%hybm).AND.(.NOT.QDEF_NCVAR(go%hybm))) THEN
     CALL COPY_NCVAR(go%hybm, gi%hybm)
  END IF

  IF (QDEF_NCVAR(gi%hybi).AND.(.NOT.QDEF_NCVAR(go%hybi))) THEN
     CALL COPY_NCVAR(go%hybi, gi%hybi)
  END IF

  ! SURFACE PRESSURE NEEDS TO BE PRE-REGRIDDED, IF NOT AVAILABLE
  ! THIS HAS TO BE DONE BEFORE (IN BALANCE_GEOHYBGRID_PS) TO
  ! ALLOW PROPER 'SORTING' OF THE GRID

END SUBROUTINE BALANCE_GEOHYBGRID
! ------------------------------------------------------------------


! ------------------------------------------------------------------
SUBROUTINE BALANCE_GEOHYBGRID_NCVAR(vari, axes, go, varo, lrgz)

  USE MESSY_MAIN_GRID_NETCDF,  ONLY: null_xtype

  IMPLICIT NONE

  ! Note: go must already be 'balanced'

  ! I/O
  TYPE(t_ncvar)     , INTENT(IN)  :: vari    ! input variable
  INTEGER           , INTENT(IN)  :: axes(3) ! dim.no of lon, lat, lev
  TYPE(t_geohybgrid), INTENT(IN)  :: go      ! output grid
  TYPE(t_ncvar)     , INTENT(INOUT) :: varo    ! output variable
  LOGICAL           , INTENT(IN), OPTIONAL :: lrgz

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'BALANCE_GEOHYBGRID_NCVAR'
  INTEGER                            :: i
  INTEGER                            :: vtype
  INTEGER, DIMENSION(:), ALLOCATABLE :: dimvec  ! dimension vector
  INTEGER                            :: status
  LOGICAL                            :: llrgz
  TYPE(t_ncvar)                      :: llatm
  TYPE(t_ncvar)                      :: llonm

  ! INIT
  !
  IF (PRESENT(lrgz)) THEN
     llrgz = lrgz
  ELSE
     llrgz = .TRUE.
  ENDIF

  CALL INIT_NCVAR(varo)

  CALL COPY_NCVAR(varo, vari)
  !
  ALLOCATE(dimvec(vari%ndims), STAT=status)
  CALL ERRMSG(substr,status,1)
  !
  vtype = vari%dat%type

  CALL INIT_NCVAR(llonm)
  CALL INIT_NCVAR(llatm)
  
   IF (go%lonm%xtype /= NULL_XTYPE) THEN
      CALL COPY_NCVAR(llonm, go%lonm)
      CALL COPY_NCVAR(llatm, go%latm)
   ELSEIF (go%clonm%xtype /= NULL_XTYPE) THEN
      CALL COPY_NCVAR(llonm, go%clonm)
      CALL COPY_NCVAR(llatm, go%clatm)
   ELSEIF (go%ulonm%xtype /= NULL_XTYPE) THEN
      CALL COPY_NCVAR(llonm, go%ulonm)
      CALL COPY_NCVAR(llatm, go%ulatm)
   ENDIF

  ! CHANGE DIMENSIONS
  DO i=1, vari%ndims
     dimvec(i) = vari%dim(i)%len
     IF (i == axes(GORD_LON)) THEN
        dimvec(i)   = llonm%dim(1)%len
        varo%dim(i) = llonm%dim(1)
     END IF
     IF (i == axes(GORD_LAT)) THEN
        IF (llatm%ndims == 2) THEN
           dimvec(i)   = llatm%dim(2)%len
           varo%dim(i) = llatm%dim(2)
        ELSE
           dimvec(i)   = llatm%dim(1)%len
           varo%dim(i) = llatm%dim(1)
        END IF
     END IF
     IF (i == axes(GORD_LEV) .AND. llrgz) THEN
        IF (go%hyam%dat%type /= VTYPE_UNDEF) THEN
           dimvec(i) = go%hyam%dim(1)%len
           varo%dim(i) = go%hyam%dim(1)
        ELSE
           dimvec(i) = go%hybm%dim(1)%len
           varo%dim(i) = go%hybm%dim(1)
        END IF
     END IF
  END DO

  CALL INIT_NARRAY(varo%dat,varo%ndims,dimvec,vtype)

  ! CLEAN UP
  DEALLOCATE(dimvec, STAT=status)
  CALL ERRMSG(substr,status,2)

  CALL INIT_NCVAR(llonm)
  CALL INIT_NCVAR(llatm)

END SUBROUTINE BALANCE_GEOHYBGRID_NCVAR
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE BALANCE_GEOHYBGRID_TIME(gi, go, lint)

  IMPLICIT NONE

  ! I/O
  TYPE (t_geohybgrid),    INTENT(INOUT) :: gi, go
  LOGICAL,                INTENT(IN)    :: lint

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'BALANCE_GEOHYBGRID_TIME'
  INTEGER                     :: i
  INTEGER                     :: uii, uio

  ! --- TIME BALANCING --- PART 1: SURFACE PRESSURE
  ! UNLIMITED ID OF PS
  !
  ! GET POS. OF TIME-DIMENSION-ID IN INPUT-GRID-PS
  uii = 0
  DO i=1, gi%ps%ndims
     IF (.NOT.QDEF_NCVAR(gi%timem)) EXIT
     IF (QCMP_NCDIM(gi%ps%dim(i), gi%timem%dim(1)) > 1) THEN
        uii = i
        EXIT
     END IF
  END DO
  !
  ! GET POS. OF TIME-DIMENSION-ID IN OUTPUT-GRID-PS
  uio = 0
  DO i=1, go%ps%ndims
     IF (.NOT.QDEF_NCVAR(go%timem)) EXIT
     IF (QCMP_NCDIM(go%ps%dim(i), go%timem%dim(1)) > 1) THEN
        uio = i
        EXIT
     END IF
  END DO
  !
  IF (lint) THEN
     IF (uio /= 0) THEN
        CALL RGMSG(substr, RGMLI, &
             '(PS): '''//TRIM(go%ps%dim(uio)%name)//''' -> '''&
             &//TRIM(gi%timem%dim(1)%name)//'''')
        CALL COPY_NCDIM(go%ps%dim(uio), gi%timem%dim(1))
     END IF
  ELSE
     IF (uii /= 0) THEN
        CALL RGMSG(substr, RGMLI, &
             '(PS): '''//TRIM(gi%ps%dim(uii)%name)//''' -> '''&
             &//TRIM(go%timem%dim(1)%name)//'''')
        CALL COPY_NCDIM(gi%ps%dim(uii), go%timem%dim(1))
     END IF
  END IF

  ! --- TIME BALANCING --- PART 2: GRID
  IF ((QDEF_NCVAR(gi%timem).AND.(.NOT.QDEF_NCVAR(go%timem))).OR. &
       lint) THEN
     CALL RGMSG(substr, RGMLI, &
          ''''//TRIM(go%timem%name)//''' -> '''//TRIM(gi%timem%name)//'''')
     CALL COPY_NCVAR(go%timem, gi%timem)
  ELSE
     CALL RGMSG(substr, RGMLI, &
          ''''//TRIM(gi%timem%name)//''' -> '''//TRIM(go%timem%name)//'''')
     CALL COPY_NCVAR(gi%timem, go%timem)
  END IF

  IF ((QDEF_NCVAR(gi%timei).AND.(.NOT.QDEF_NCVAR(go%timei))).OR. &
       lint) THEN
     CALL RGMSG(substr, RGMLI, &
          ''''//TRIM(go%timei%name)//''' -> '''//TRIM(gi%timei%name)//'''')
     CALL COPY_NCVAR(go%timei, gi%timei)
  ELSE
     CALL RGMSG(substr, RGMLI, &
          ''''//TRIM(gi%timei%name)//''' -> '''//TRIM(go%timei%name)//'''')
     CALL COPY_NCVAR(gi%timei, go%timei)
  END IF

  ! --- TIME BALANCING --- PART 3: TIME STEP
  IF (lint) THEN
     go%t = gi%t
  ELSE
     gi%t = go%t
  END IF

END SUBROUTINE BALANCE_GEOHYBGRID_TIME
! ------------------------------------------------------------------

SUBROUTINE REDUCE_INPUT_GRID(gs, gd)

  USE MESSY_MAIN_GRID_NETCDF,  ONLY: null_xtype
  USE MESSY_MAIN_GRID,         ONLY: IMPORT_GEOHYBGRID

  IMPLICIT NONE

  ! I/O
  TYPE(t_geohybgrid), INTENT(INOUT) :: gs   ! INPUT GRID (will be changed)
  TYPE(t_geohybgrid), INTENT(IN)    :: gd   ! OUTPUT (destination) GRID 

  ! LOCAL

  CHARACTER(LEN=*), PARAMETER :: substr='REDUCE_INPUT_GRID'
  TYPE(t_narray)              :: llons, llond
  REAL(dp)                    :: dlon, dlontmp
  TYPE(t_narray)              :: llats, llatd
  REAL(dp)                    :: dlat
  REAL(dp)                    :: lonmax, lonmin, latmax, latmin
  INTEGER                     :: locidim, locjdim
  INTEGER                     :: n, ix , jx
  INTEGER                     :: endi, endj
  TYPE(t_geohybgrid)          :: gi
  LOGICAL                     :: lwrap
  LOGICAL                     :: lusefull = .FALSE.
  INTEGER                     :: tag
  LOGICAL                     :: lonelong = .FALSE.
  INTEGER, DIMENSION(2)       :: hdids = -99

  IF (ALL(gd%minmaxlonlat < -360._dp)) RETURN

  ! INITIALISE MARKER
  ivaldefd = I_UNDEF
  ivaldefs = I_UNDEF
  lusefull = .FALSE.
  lwrap    = .FALSE.
  tag      = -99

  CALL INIT_NARRAY(llons)
  CALL INIT_NARRAY(llond)
  CALL COPY_GEOHYBGRID(gi, gs)
  CALL COMPLETE_GEOHYBGRID(gi)

  IF (gi%clonm%xtype == NULL_XTYPE .AND. gi%lonm%xtype == NULL_XTYPE &
       .AND. gi%loni%xtype == NULL_XTYPE &
       .AND. gi%ulonm%xtype == NULL_XTYPE) THEN
      CALL RGMSG(substr, RGMLI, &
          ' SKIPPING GRID REDUCTION ', .FALSE.)
     RETURN
  ENDIF
  ! I)   The source field should be greater than the destination field
  ! Ia) longitude field
  IF (gi%lonm%xtype /= NULL_XTYPE) THEN
     CALL COPY_NARRAY(llons, gi%lonm%dat)
     locidim = gi%lonm%dim(1)%len
     hdids(1) = gi%lonm%dim(1)%id
  ELSE IF (gi%clonm%xtype /= NULL_XTYPE) THEN
     CALL COPY_NARRAY(llons, gi%clonm%dat)
     locidim = gi%clonm%dim(1)%len
     hdids(1) = gi%clonm%dim(1)%id
  ELSE IF (gi%loni%xtype /= NULL_XTYPE) THEN
     CALL COPY_NARRAY(llons, gi%loni%dat)
     locidim = gi%loni%dim(1)%len
     hdids(1) = gi%loni%dim(1)%id
  ELSE IF (gi%ulonm%xtype /= NULL_XTYPE) THEN
     CALL COPY_NARRAY(llons, gi%ulonm%dat)
     locidim = gi%ulonm%dim(1)%len
     hdids(1) = gi%ulonm%dim(1)%id
  ELSE IF (gi%uloni%xtype /= NULL_XTYPE) THEN
     CALL COPY_NARRAY(llons, gi%uloni%dat)
     locidim = gi%uloni%dim(1)%len
     hdids(1) = gi%uloni%dim(1)%id
  ENDIF
  ! For easier handling double array
  IF (ASSOCIATED(llons%vr)) CALL DOUBLE_NARRAY(llons)
  IF (SIZE(llons%vd) == 1) THEN
     lonelong = .TRUE.
  ELSE
     lonelong = .FALSE.
  ENDIF

  IF  (gd%loni%xtype /= NULL_XTYPE) THEN
     CALL COPY_NARRAY(llond, gd%loni%dat)
  ELSE IF (gd%cloni%xtype /= NULL_XTYPE) THEN
     CALL COPY_NARRAY(llond, gd%cloni%dat)
  ELSE IF (gd%clonm%xtype /= NULL_XTYPE) THEN
     CALL COPY_NARRAY(llond, gd%clonm%dat)
  ELSE IF (gd%lonm%xtype /= NULL_XTYPE) THEN
     CALL COPY_NARRAY(llond, gd%lonm%dat)
  ELSE IF (gd%ulonm%xtype /= NULL_XTYPE) THEN
     CALL COPY_NARRAY(llond, gd%ulonm%dat)
  ENDIF
  ! For easier handly double array
  IF (ASSOCIATED(llond%vr)) CALL DOUBLE_NARRAY(llond)

  IF (lonelong) THEN
     dlon = 2._dp
  ELSE  IF (gs%ulonm%xtype /= NULL_XTYPE .OR. gs%uloni%xtype /= NULL_XTYPE) THEN
     ! dlon needs to be approximated otherwise
     ! TODO do it better
     dlon = 2._dp
  ELSE
     dlon = ABS(llons%vd(2) - llons%vd(1))
  ENDIF
  dlon = 2._dp * dlon

  ! TEST IF BOTH GRIDs are defined on the same region [-180,180[ or [0,360[
  gs%minmaxlonlat(1,1) = MINVAL(llons%vd)
  gs%minmaxlonlat(1,2) = MAXVAL(llons%vd)

  !!!! ***********************************************
  !!!! **** WORK AROUND for better results in ICON 
  !!!! **** or unstructured grids
  !!!! ***********************************************
  lunstruct = .FALSE.
  IF  (gd%ulonm%xtype /= NULL_XTYPE) lunstruct = .TRUE.  ! ICON OR CESM

  IF (.NOT. lunstruct) THEN
     lonmin = MINVAL(llond%vd, mask=llond%vd > -360._dp) - dlon
     lonmax = MAXVAL(llond%vd) + dlon
     dlontmp = 0._dp
  ELSE
     lonmin = gd%minmaxlonlat(1,1) 
     lonmax = gd%minmaxlonlat(1,2)
     dlontmp = dlon
  END IF

  IF (lonmin < 0._dp) THEN
     ivaldefd = I_m180_180  ! destination grid defined on [-180,180]
!???  ELSE IF(lonmax + dlontmp <= 180._dp) THEN
  ELSE IF(lonmax + dlontmp < 180._dp) THEN
     IF (lonmin - dlontmp <= 0._dp) THEN
        ivaldefd = I_m180_180  ! destination grid defined on [-180,180]
     ELSE
        ivaldefd = I_0_180
     END IF
  ELSE 
     ivaldefd = I_0_360
  END IF

  IF (lonelong) THEN
     ! calculation of ivaldefd required for EXPAND_INPUT_GRID_LON
     ! RETURN moved here
     CALL RGMSG(substr, RGMLI, &
          ' SKIPPING GRID REDUCTION: DIM(lon) = 1 ', .FALSE.)
     RETURN
  END IF

  ! SOURCE GRID INTERVAL
   IF (gs%minmaxlonlat(1,1) < 0._dp) THEN
     ivaldefs = I_m180_180  ! destination grid defined on [-180,180]
  ELSE IF (gs%minmaxlonlat(1,2) > 180._dp) THEN
     ivaldefs  = I_0_360
  ELSE
     ivaldefs = I_0_180
  END IF

  IF (.not. lunstruct) THEN
  IF ( (lonmax - lonmin) > 359._dp) THEN
     lusefull = .TRUE.
     tag = 25
  END IF
  IF (ivaldefd == ivaldefs) THEN
     ! EVERY THING OK destination grid smaller as source grid
      tag = 1

  ELSE  IF (ivaldefd == I_0_180) THEN
     ! FITS in BOTH INTERVALS OF INPUT DATA
     tag = 2
  ELSE  IF (ivaldefd == I_m180_180 .AND. ivaldefs == I_0_360) THEN
           IF (lonmin < 0._dp)  lonmin = lonmin + 360._dp
           IF (lonmax < 0._dp)  lonmax = lonmax + 360._dp
           tag = 8
           lwrap = (lonmin > lonmax)
     ELSE IF (ivaldefd == I_0_360 .AND. ivaldefs == I_m180_180 ) THEN
           IF (lonmin > 180._dp) lonmin = lonmin - 360._dp
           IF (lonmax > 180._dp) lonmax = lonmax - 360._dp

           lwrap = (lonmin > lonmax)
           tag = 10
     END IF

     IF (tag == -99) THEN
        lusefull = .TRUE.
        write (0,*) 'ERROR IMPORT_GRID USE FULL: ',tag
     END IF
  ELSE ! lunstruct
  IF (gd%minmaxlonlat(1,1) - dlon > gs%minmaxlonlat(1,1) .AND. &
       gd%minmaxlonlat(1,2) + dlon < gs%minmaxlonlat(1,2)) THEN
     ! EVERY THING OK destination grid smaller as source grid
     tag = 1
  ELSE IF ( (gd%minmaxlonlat(1,1) >= gs%minmaxlonlat(1,2)) .OR. &
          (gd%minmaxlonlat(1,2) <= gs%minmaxlonlat(1,1))) THEN
        lusefull = .TRUE.
        tag = 6
  ELSE IF (ivaldefd == ivaldefs .OR. ivaldefs == I_0_180) THEN
     ! grids are defined on overlapping intervals,  but min_dest < min_source
     ! or max_dest > max_source => destination grid larger as source grid
     ! use full source grid
     lusefull = .TRUE.
     tag = 7
  ELSE  IF (ivaldefd /= I_0_180) THEN

     IF (ivaldefd == I_m180_180 .AND. ivaldefs == I_0_360) THEN
        IF (gd%minmaxlonlat(1,2) + dlon < 0._dp) THEN
           CALL RGMSG(substr, RGMLI, &
                ' LONGITUDE SWITCH REQUIRED -360', .FALSE.)
           DO ix = 1, SIZE(llons%vd)
              IF ( llons%vd(ix) > 180.0_dp) THEN
                 llons%vd(ix) =  llons%vd(ix) - 360.0_dp
              ENDIF
           END DO
           tag = 8

        ELSE
           lusefull = .TRUE.
           tag = 9
        END IF
     ELSE IF (ivaldefd == I_0_360 .AND. ivaldefs == I_m180_180 ) THEN
        IF (gd%minmaxlonlat(1,1) - dlon > 180._dp) THEN
           CALL RGMSG(substr, RGMLI, &
                ' LONGITUDE SWITCH REQUIRED +360', .FALSE.)
           DO ix = 1, SIZE(llons%vd)
              IF ( llons%vd(ix) < 0.0_dp) THEN
                 llons%vd(ix) =  llons%vd(ix) + 360.0_dp
              ENDIF
           END DO
           tag = 10
        ELSE
           lusefull = .TRUE.
           tag = 11
        END IF
     END IF
     gs%minmaxlonlat(1,1) = MINVAL(llons%vd)
     gs%minmaxlonlat(1,2) = MAXVAL(llons%vd)
     IF (gd%minmaxlonlat(1,1) - dlon < gs%minmaxlonlat(1,1) .AND. &
          gd%minmaxlonlat(1,2) + dlon > gs%minmaxlonlat(1,2)) THEN
        lusefull = .TRUE.
        tag = tag + 20

     END IF
     IF (tag == -99) THEN
        lusefull = .TRUE.
        tag = tag + 30
     END IF
  ELSE
     ! OK
     tag = 40
  END IF

  END IF ! lunstruct
  ! ----------------------
  ! LATITUDE RANGE
  ! ----------------------

  CALL INIT_NARRAY(llats)
  CALL INIT_NARRAY(llatd)

  ! I)   The source field should be greater than the destination field
  ! Ia) longitude field

  IF (gi%latm%xtype /= NULL_XTYPE) THEN
     CALL COPY_NARRAY(llats, gi%latm%dat)
     locjdim = gi%latm%dim(1)%len
     hdids(2) = gi%latm%dim(1)%id
  ELSE IF (gi%clatm%xtype /= NULL_XTYPE) THEN
     CALL COPY_NARRAY(llats, gi%clatm%dat)
     locjdim = gi%clatm%dim(2)%len
     hdids(2) = gi%clatm%dim(2)%id
  ELSE IF (gi%lati%xtype /= NULL_XTYPE) THEN
     CALL COPY_NARRAY(llats, gi%lati%dat)
     locjdim = gi%lati%dim(1)%len
     hdids(2) = gi%lati%dim(1)%id
  ELSE IF (gi%ulatm%xtype /= NULL_XTYPE) THEN
     CALL COPY_NARRAY(llats, gi%ulatm%dat)
     locjdim = gi%ulatm%dim(2)%len
     hdids(2) = gi%ulatm%dim(2)%id
  ENDIF
  ! For easier handling double array
  IF (ASSOCIATED(llats%vr)) CALL DOUBLE_NARRAY(llats)

  IF (SIZE(llats%vd) == 1) THEN
      CALL RGMSG(substr, RGMLI, &
          ' SKIPPING GRID REDUCTION: DIM(lat) = 1 ', .FALSE.)
     RETURN
  ENDIF

  IF (gd%latm%xtype /= NULL_XTYPE) THEN
     CALL COPY_NARRAY(llatd, gd%latm%dat)
  ELSE IF (gd%clatm%xtype /= NULL_XTYPE) THEN
     CALL COPY_NARRAY(llatd, gd%clatm%dat)
  ELSE IF (gd%lati%xtype /= NULL_XTYPE) THEN
     CALL COPY_NARRAY(llatd, gd%lati%dat)
  ELSE IF (gd%ulatm%xtype /= NULL_XTYPE) THEN
     CALL COPY_NARRAY(llatd, gd%ulatm%dat)
  ENDIF
  ! For easier handly double array
  IF (ASSOCIATED(llatd%vr)) CALL DOUBLE_NARRAY(llatd)

  dlat = ABS(llats%vd(2) - llats%vd(1))
  IF (gs%ulatm%xtype /= NULL_XTYPE) THEN
     ! dlon needs to be approximated otherwise
     ! TODO do it better
     dlat = 2.0_dp
  ENDIF

  gs%minmaxlonlat(2,1) = MINVAL(llats%vd)
  gs%minmaxlonlat(2,2) = MAXVAL(llats%vd)

  ! TEST IF DESTINATION GRID is part of source grid 
  IF (gd%minmaxlonlat(2,1) > gs%minmaxlonlat(2,1) .AND. &
       gd%minmaxlonlat(2,2) < gs%minmaxlonlat(2,2)) THEN
     ! EVERYTHING is fine
  ELSE
     CALL RGMSG(substr, RGMLVL, &
          'DESTINATION LATITUDE RANGE LARGER THAN SOURCE GRID', .FALSE.)
     RETURN
  ENDIF

  IF (.NOT. lunstruct) THEN
     lonmin = lonmin - 3*dlon
     lonmax = lonmax + 3*dlon
     latmin = gd%minmaxlonlat(2,1) - 3*dlat
     latmax = gd%minmaxlonlat(2,2) + 3*dlat
  ELSE
     lonmin = gd%minmaxlonlat(1,1) - dlon
     lonmax = gd%minmaxlonlat(1,2) + dlon
     latmin = gd%minmaxlonlat(2,1) - dlat
     latmax = gd%minmaxlonlat(2,2) + dlat
  END IF

  ! latitudes are increasing
  gs%start(1) = locidim
  gs%start(2) = locjdim
  endi        = 1
  endj        = 1

  ! Test, if we have a destination grid which wraps around
  IF (lunstruct) THEN
     lwrap = .FALSE.
     IF (gd%minmaxlonlat(1,1) > gd%minmaxlonlat(1,2)) lwrap = .TRUE.
  END IF
  ! Distinguish between rectangular non rotated (lon/lat1d) and 
  ! curvilinear (lon/lat 2d) grids or unstructured grids
  IF (llons%n == 1 .AND. llats%n == 1) THEN
     ! rectangular non rotated grid 
     IF (lwrap) THEN
        DO ix = 1, locidim
           IF (llons%vd(ix) < lonmax) THEN
              endi = MAX(ix, endi)
           END IF
           IF (llons%vd(ix) > lonmin) THEN
              gs%start(1) = MIN(ix,gs%start(1))
           END IF
        END DO
     ELSE
        DO ix = 1, locidim
           IF (llons%vd(ix) > lonmin .AND. llons%vd(ix) < lonmax ) THEN
              gs%start(1) = MIN(ix,gs%start(1))
              endi        = MAX(ix,endi)
           END IF
        END DO
     END IF
     DO jx = 1, locjdim
        IF (llats%vd(jx) > latmin .AND. llats%vd(jx) < latmax ) THEN
           gs%start(2) = MIN(jx,gs%start(2))
           endj        = MAX(jx,endj)
        END IF
     END DO
  ELSE
     ! curvilinear or unstructured
     DO ix = 1, locidim
        DO jx = 1, locjdim
           n = (jx -1) * locidim + ix
           IF (lwrap) THEN
              IF (llats%vd(n) > latmin .AND. llats%vd(n) < latmax .AND. &
                 &(llons%vd(n) < lonmax .OR. llons%vd(n) > lonmin)) THEN
                 ! Point is required
                 gs%start(1) = MIN(ix,gs%start(1))
                 gs%start(2) = MIN(jx,gs%start(2))
                 endi        = MAX(ix,endi)
                 endj        = MAX(jx,endj)
              END IF
           ELSE
              IF (llats%vd(n) > latmin .AND. llats%vd(n) < latmax .AND. &
                   llons%vd(n) > lonmin .AND. llons%vd(n) < lonmax ) THEN
                 ! Point is required
                 gs%start(1) = MIN(ix,gs%start(1))
                 gs%start(2) = MIN(jx,gs%start(2))
                 endi        = MAX(ix,endi)
                 endj        = MAX(jx,endj)
              END IF
           END IF ! lwrap
        END DO
     END DO
  ENDIF
 
  IF (gd%minmaxlonlat(1,1) == RGEMPTY .OR. &
       gd%minmaxlonlat(1,2) == RGEMPTY) THEN
     gs%start(1) = 1
     endi        = locidim
  END IF
  IF (gd%minmaxlonlat(2,1) == RGEMPTY .OR. &
       gd%minmaxlonlat(2,2) == RGEMPTY) THEN 
     gs%start(2) = 1
     endj        = locjdim
  END IF

  IF (lwrap) THEN !! .AND. .NOT.(llons%n == 1 .AND. llats%n == 1)) THEN
     IF (.NOT. lunstruct) THEN
        gs%count(1) = locidim - gs%start(1) + 1 + endi
        gs%count(2) = endj - gs%start(2) + 1
     ELSE IF (.NOT.(llons%n == 1 .AND. llats%n == 1)) THEN
        gs%count(1) = locidim - gs%start(1) + 1 + endi
        gs%count(2) = endj - gs%start(2) + 1
     END IF
  ELSE
     gs%count(1) = endi - gs%start(1) + 1
     gs%count(2) = endj - gs%start(2) + 1
  END IF
  IF (lusefull) THEN
     IF (.NOT. lunstruct) THEN
        gs%start(1) = 1
        gs%count(1) = -99
     ELSE
        gs%start(:) = 1
        gs%count(:) = -99
     END IF
     write (0,*) 'FLIP TAKE ALL LONGITUDES (usefull)'
  END IF

  CALL INIT_NARRAY(llond)
  CALL INIT_NARRAY(llons)
  CALL INIT_NARRAY(llats)
  CALL INIT_NARRAY(llatd)

  CALL INIT_GEOHYBGRID(gi)

  ! ------------------------
  ! DETERIMINE REDUCED GRID
  ! ------------------------

  ! in case of reduced grid, no modulo axis exist
  IF (gs%count(1) /= locidim) THEN 
     gs%lonc  = .FALSE.
     gs%clonc = .FALSE.
     gs%rlonc = .FALSE.
  ENDIF
  ! set horizontal dimension ids for GEOHYBGRID import
  gs%hdimids = hdids

  CALL IMPORT_GEOHYBGRID(gs, gs%start, gs%count)

  ! re-set horizontal dimension ids for GEOHYBGRID because these
  ! were reset to empty strings during import.
  gs%hdimids = hdids

  ! SET RANGES
  IF (.NOT. lusefull)  gs%ranges(1,:)  = RGEMPTY
  gs%minmaxlonlat(:,:) = RGEMPTY

  IF (.NOT. lunstruct) THEN
     IF (lwrap .AND. .NOT. lusefull) THEN
        IF (ivaldefd == I_m180_180) THEN
           CALL FLIP_LONS(gs, -1)
        ELSE IF (ivaldefd == I_0_360) THEN
           CALL FLIP_LONS(gs, 1)
        END IF
     END IF
  ENDIF !
!  write (0,*) 'FLIP2 ', ivaldefd ,I_m180_180, I_0_360, gs%minmaxlonlat(1,:)

CONTAINS

    SUBROUTINE FLIP_LONS(g, flag)

      TYPE(t_geohybgrid), INTENT(INOUT) :: g
      INTEGER,            INTENT(IN)    :: flag

      !LOCAL
      TYPE(t_narray)                    :: llon

      SELECT CASE (flag)
      CASE (1)
         IF (ASSOCIATED(g%lonm%dat%vr)) THEN
            WHERE (g%lonm%dat%vr < 0._sp .AND. g%lonm%dat%vr > -360._sp) &
               g%lonm%dat%vr = g%lonm%dat%vr + 360._sp
            CALL COPY_NARRAY(llon, g%lonm%dat)
            CALL DOUBLE_NARRAY(llon)
         END IF
         IF (ASSOCIATED(g%loni%dat%vr)) THEN
            WHERE (g%loni%dat%vr < 0._sp.AND. g%loni%dat%vr > -360._sp) &
               g%loni%dat%vr = g%loni%dat%vr + 360._sp
            CALL COPY_NARRAY(llon, g%loni%dat)
            CALL DOUBLE_NARRAY(llon)
         END IF
         IF (ASSOCIATED(g%clonm%dat%vr)) THEN
            WHERE (g%clonm%dat%vr < 0._sp.AND. g%clonm%dat%vr > -360._sp) &
               g%clonm%dat%vr = g%clonm%dat%vr + 360._sp
            CALL COPY_NARRAY(llon, g%clonm%dat)
            CALL DOUBLE_NARRAY(llon)
         END IF
         IF (ASSOCIATED(g%loni%dat%vr)) THEN
            WHERE (g%cloni%dat%vr < 0._sp.AND. g%cloni%dat%vr > -360._sp) &
               g%loni%dat%vr = g%cloni%dat%vr + 360._sp
            CALL COPY_NARRAY(llon, g%cloni%dat)
            CALL DOUBLE_NARRAY(llon)
         END IF
         IF (ASSOCIATED(g%ulonm%dat%vr)) THEN
            WHERE (g%ulonm%dat%vr < 0._sp.AND. g%ulonm%dat%vr > -360._sp) &
               g%ulonm%dat%vr = g%ulonm%dat%vr + 360._sp
            CALL COPY_NARRAY(llon, g%ulonm%dat)
            CALL DOUBLE_NARRAY(llon)
         END IF
         IF (ASSOCIATED(g%uloni%dat%vr)) THEN
            WHERE (g%uloni%dat%vr < 0._sp.AND. g%uloni%dat%vr > -360._sp) &
               g%uloni%dat%vr = g%uloni%dat%vr + 360._sp
            CALL COPY_NARRAY(llon, g%uloni%dat)
            CALL DOUBLE_NARRAY(llon)
         END IF
         IF (ASSOCIATED(g%lonm%dat%vd)) THEN
            WHERE (g%lonm%dat%vd < 0._dp .AND. g%lonm%dat%vd > -360._dp) &
               g%lonm%dat%vd = g%lonm%dat%vd + 360._dp
            CALL COPY_NARRAY(llon, g%lonm%dat)
         END IF
         IF (ASSOCIATED(g%loni%dat%vd)) THEN
            WHERE (g%loni%dat%vd < 0._dp .AND. g%loni%dat%vd > -360._dp) &
               g%loni%dat%vd = g%loni%dat%vd + 360._dp
            CALL COPY_NARRAY(llon, g%loni%dat)
         END IF
         IF (ASSOCIATED(g%clonm%dat%vd)) THEN
            WHERE (g%clonm%dat%vd < 0._dp .AND. g%clonm%dat%vd > -360._dp) &
               g%clonm%dat%vd = g%clonm%dat%vd + 360._dp
            CALL COPY_NARRAY(llon, g%clonm%dat)
         END IF
         IF (ASSOCIATED(g%loni%dat%vd)) THEN
            WHERE (g%cloni%dat%vd < 0._dp .AND. g%cloni%dat%vd > -360._dp) &
               g%loni%dat%vd = g%cloni%dat%vd + 360._dp
            CALL COPY_NARRAY(llon, g%cloni%dat)
         END IF
         IF (ASSOCIATED(g%ulonm%dat%vd)) THEN
            WHERE (g%ulonm%dat%vd < 0._dp .AND. g%ulonm%dat%vd > -360._dp) &
               g%ulonm%dat%vd = g%ulonm%dat%vd + 360._dp
            CALL COPY_NARRAY(llon, g%ulonm%dat)
         END IF
         IF (ASSOCIATED(g%uloni%dat%vd)) THEN
            WHERE (g%uloni%dat%vd < 0._dp .AND. g%uloni%dat%vd > -360._dp) &
               g%uloni%dat%vd = g%uloni%dat%vd + 360._dp
            CALL COPY_NARRAY(llon, g%uloni%dat)
         END IF


      CASE (-1)

         IF (ASSOCIATED(g%lonm%dat%vr)) THEN
            WHERE (g%lonm%dat%vr > 180._sp) &
               g%lonm%dat%vr = g%lonm%dat%vr - 360._sp
            CALL COPY_NARRAY(llon, g%lonm%dat)
            CALL DOUBLE_NARRAY(llon)
         END IF
         IF (ASSOCIATED(g%loni%dat%vr)) THEN
            WHERE (g%loni%dat%vr > 180._sp) &
               g%loni%dat%vr = g%loni%dat%vr - 360._sp
            CALL COPY_NARRAY(llon, g%loni%dat)
            CALL DOUBLE_NARRAY(llon)
         END IF
         IF (ASSOCIATED(g%clonm%dat%vr)) THEN
            WHERE (g%clonm%dat%vr > 180._sp) &
               g%clonm%dat%vr = g%clonm%dat%vr - 360._sp
            CALL COPY_NARRAY(llon, g%clonm%dat)
            CALL DOUBLE_NARRAY(llon)
         END IF
         IF (ASSOCIATED(g%loni%dat%vr)) THEN
            WHERE (g%cloni%dat%vr > 180._sp) &
               g%loni%dat%vr = g%cloni%dat%vr - 360._sp
            CALL COPY_NARRAY(llon, g%cloni%dat)
            CALL DOUBLE_NARRAY(llon)
         END IF
         IF (ASSOCIATED(g%ulonm%dat%vr)) THEN
            WHERE (g%ulonm%dat%vr > 180._sp) &
               g%ulonm%dat%vr = g%ulonm%dat%vr - 360._sp
            CALL COPY_NARRAY(llon, g%ulonm%dat)
            CALL DOUBLE_NARRAY(llon)
         END IF
         IF (ASSOCIATED(g%uloni%dat%vr)) THEN
            WHERE (g%uloni%dat%vr > 180._sp) &
               g%uloni%dat%vr = g%uloni%dat%vr - 360._sp
            CALL COPY_NARRAY(llon, g%uloni%dat)
            CALL DOUBLE_NARRAY(llon)
         END IF
         IF (ASSOCIATED(g%lonm%dat%vd)) THEN
            WHERE (g%lonm%dat%vd > 180._dp) &
               g%lonm%dat%vd = g%lonm%dat%vd - 360._dp
            CALL COPY_NARRAY(llon, g%lonm%dat)
         END IF
         IF (ASSOCIATED(g%loni%dat%vd)) THEN
            WHERE (g%loni%dat%vd > 180._dp) &
               g%loni%dat%vd = g%loni%dat%vd - 360._dp
            CALL COPY_NARRAY(llon, g%loni%dat)
          END IF
         IF (ASSOCIATED(g%clonm%dat%vd)) THEN
            WHERE (g%clonm%dat%vd > 180._dp) &
               g%clonm%dat%vd = g%clonm%dat%vd - 360._dp
            CALL COPY_NARRAY(llon, g%clonm%dat)
         END IF
         IF (ASSOCIATED(g%loni%dat%vd)) THEN
            WHERE (g%cloni%dat%vd > 180._dp) &
               g%loni%dat%vd = g%cloni%dat%vd - 360._dp
           CALL COPY_NARRAY(llon, g%cloni%dat)
          END IF
         IF (ASSOCIATED(g%ulonm%dat%vd)) THEN
            WHERE (g%ulonm%dat%vd > 180._dp) &
               g%ulonm%dat%vd = g%ulonm%dat%vd - 360._dp
           CALL COPY_NARRAY(llon, g%ulonm%dat)
          END IF
         IF (ASSOCIATED(g%uloni%dat%vd)) THEN
            WHERE (g%uloni%dat%vd > 180._dp) &
                 g%uloni%dat%vd = g%uloni%dat%vd - 360._dp
            CALL COPY_NARRAY(llon, g%uloni%dat)
         END IF

      END SELECT

      g%minmaxlonlat(1,1) = MINVAL(llon%vd, mask= llon%vd > -360._dp)
      g%minmaxlonlat(1,2) = MAXVAL(llon%vd)

      CALL INIT_NARRAY(llon)

    END SUBROUTINE FLIP_LONS
END SUBROUTINE REDUCE_INPUT_GRID
! ------------------------------------------------------------------

!==================================================================
!==================================================================
! COPYied FROM COSMO MODEL (required here)
!==================================================================
!==================================================================

!==============================================================================
!==============================================================================
!+ Function for rotation of geographical coordinates
!------------------------------------------------------------------------------

FUNCTION  PHIROT2PHI ( phirot, rlarot, polphi, polgam )

!------------------------------------------------------------------------------
!
! Description:
!   This function converts phi from one rotated system to phi in another
!   system. If the optional argument polgam is present, the other system
!   can also be a rotated one, where polgam is the angle between the two
!   north poles.
!   If polgam is not present, the other system is the real geographical
!   system.
!
! Method:
!   Transformation formulas for converting between these two systems.
!
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!
! Declarations:
!
!------------------------------------------------------------------------------

! Parameter list:
REAL (KIND=dp), INTENT (IN)      ::        &
  polphi,   & ! latitude of the rotated north pole
  phirot,   & ! latitude in the rotated system
  rlarot      ! longitude in the rotated system

REAL (KIND=dp), INTENT (IN)      ::        &
  polgam      ! angle between the north poles of the systems

REAL (KIND=dp)                   ::        &
  phirot2phi  ! latitude in the geographical system

! Local variables
REAL (KIND=dp)                   ::        &
  zsinpol, zcospol, zphis, zrlas, zarg, zgam

REAL (KIND=dp), PARAMETER        ::        &
  zrpi18 = 57.2957795_dp,                  &
  zpir18 = 0.0174532925_dp

!------------------------------------------------------------------------------

! Begin function phirot2phi

  zsinpol     = SIN (zpir18 * polphi)
  zcospol     = COS (zpir18 * polphi)
 
  zphis       = zpir18 * phirot
  IF (rlarot > 180.0_dp) THEN
    zrlas = rlarot - 360.0_dp
  ELSE
    zrlas = rlarot
  ENDIF
  zrlas       = zpir18 * zrlas

  IF (polgam /= 0.0_dp) THEN
    zgam  = zpir18 * polgam
    zarg  = zsinpol*SIN (zphis) +                                           &
        zcospol*COS(zphis) * ( COS(zrlas)*COS(zgam) - SIN(zgam)*SIN(zrlas) )
  ELSE
    zarg  = zcospol * COS (zphis) * COS (zrlas) + zsinpol * SIN (zphis)
  ENDIF
 
  phirot2phi  = zrpi18 * ASIN (zarg)

END FUNCTION PHIROT2PHI

!==============================================================================
!==============================================================================

!------------------------------------------------------------------------------

FUNCTION  PHI2PHIROT ( phi, rla, polphi, pollam )

!------------------------------------------------------------------------------
! Description:
!   This routine converts phi from the real geographical system to phi
!   in the rotated system.
!
! Method:
!   Transformation formulas for converting between these two systems.
!
!------------------------------------------------------------------------------
! Parameter list:
REAL (KIND=dp), INTENT (IN)      ::        &
  polphi,  & ! latitude of the rotated north pole
  pollam,  & ! longitude of the rotated north pole
  phi,     & ! latitude in the geographical system
  rla        ! longitude in the geographical system

REAL (KIND=dp)                   ::        &
  phi2phirot ! longitude in the rotated system

! Local variables
REAL (KIND=dp)                       ::    &
  zsinpol, zcospol, zlampol, zphi, zrla, zarg1, zarg2, zrla1

REAL (KIND=dp), PARAMETER            ::    &
  zrpi18 = 57.2957795_dp,                  & !
  zpir18 = 0.0174532925_dp

!------------------------------------------------------------------------------

! Begin function phi2phirot

  zsinpol  = SIN (zpir18 * polphi)
  zcospol  = COS (zpir18 * polphi)
  zlampol  =      zpir18 * pollam
  zphi     =      zpir18 * phi
  IF (rla > 180.0_dp) THEN
    zrla1  = rla - 360.0_dp
  ELSE
    zrla1  = rla
  ENDIF
  zrla     = zpir18 * zrla1

  zarg1    = SIN (zphi) * zsinpol
  zarg2    = COS (zphi) * zcospol * COS (zrla - zlampol)

  phi2phirot = zrpi18 * ASIN (zarg1 + zarg2)

END FUNCTION PHI2PHIROT

!==============================================================================
!==============================================================================

!------------------------------------------------------------------------------

FUNCTION  RLAROT2RLA (phirot, rlarot, polphi, pollam, polgam)

!------------------------------------------------------------------------------
!
! Description:
!   This function converts lambda from one rotated system to lambda in another
!   system. If the optional argument polgam is present, the other system
!   can also be a rotated one, where polgam is the angle between the two
!   north poles.
!   If polgam is not present, the other system is the real geographical
!   system.
!
! Method:
!   Transformation formulas for converting between these two systems.
!
! Modules used:    NONE
!
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!
! Declarations:
!
!------------------------------------------------------------------------------

! Parameter list:
REAL (KIND=dp), INTENT (IN)      ::        &
  polphi,   & ! latitude of the rotated north pole
  pollam,   & ! longitude of the rotated north pole
  phirot,   & ! latitude in the rotated system
  rlarot      ! longitude in the rotated system

REAL (KIND=dp), INTENT (IN)      ::        &
  polgam      ! angle between the north poles of the systems

REAL (KIND=dp)                   ::        &
  rlarot2rla  ! longitude in the geographical system

! Local variables
REAL (KIND=dp)                   ::        &
  zsinpol, zcospol, zlampol, zphis, zrlas, zarg1, zarg2, zgam

REAL (KIND=dp), PARAMETER        ::        &
  zrpi18 = 57.2957795_dp,                  & !
  zpir18 = 0.0174532925_dp

!------------------------------------------------------------------------------

! Begin function rlarot2rla

  zsinpol = SIN (zpir18 * polphi)
  zcospol = COS (zpir18 * polphi)

  zlampol = zpir18 * pollam
  zphis   = zpir18 * phirot
  IF (rlarot > 180.0_dp) THEN
    zrlas = rlarot - 360.0_dp
  ELSE
    zrlas = rlarot
  ENDIF
  zrlas   = zpir18 * zrlas

  IF (polgam /= 0.0_dp) THEN
    zgam    = zpir18 * polgam
    zarg1   = SIN (zlampol) *                                                &
      (- zsinpol*COS(zphis) * (COS(zrlas)*COS(zgam) - SIN(zrlas)*SIN(zgam))  &
       + zcospol * SIN(zphis))                                               &
    - COS (zlampol)*COS(zphis) * (SIN(zrlas)*COS(zgam) + COS(zrlas)*SIN(zgam))

    zarg2   = COS (zlampol) *                                                &
      (- zsinpol*COS(zphis) * (COS(zrlas)*COS(zgam) - SIN(zrlas)*SIN(zgam))  &
       + zcospol * SIN(zphis))                                               &
    + SIN (zlampol)*COS(zphis) * (SIN(zrlas)*COS(zgam) + COS(zrlas)*SIN(zgam))
  ELSE
    zarg1   = SIN (zlampol) * (-zsinpol * COS(zrlas) * COS(zphis)  +    &
                                zcospol *              SIN(zphis)) -    &
              COS (zlampol) *             SIN(zrlas) * COS(zphis)
    zarg2   = COS (zlampol) * (-zsinpol * COS(zrlas) * COS(zphis)  +    &
                                zcospol *              SIN(zphis)) +   &
              SIN (zlampol) *             SIN(zrlas) * COS(zphis)
  ENDIF
 
  IF (zarg2 == 0.0) zarg2 = 1.0E-20_dp
 
  rlarot2rla = zrpi18 * ATAN2(zarg1,zarg2)
 
END FUNCTION RLAROT2RLA

!==============================================================================
!==============================================================================

!------------------------------------------------------------------------------

FUNCTION  RLA2RLAROT ( phi, rla, polphi, pollam, polgam )

!------------------------------------------------------------------------------
!
! Description:
!   This routine converts lambda from the real geographical system to lambda 
!   in the rotated system.
!
! Method:
!   Transformation formulas for converting between these two systems.
!
!------------------------------------------------------------------------------
!
! Parameter list:
REAL (KIND=dp), INTENT (IN)      ::        &
  polphi,  & ! latitude of the rotated north pole
  pollam,  & ! longitude of the rotated north pole
  phi,     & ! latitude in geographical system
  rla        ! longitude in geographical system

REAL (KIND=dp), INTENT (IN)      ::        &
  polgam      ! angle between the north poles of the systems

REAL (KIND=dp)                   ::        &
  rla2rlarot ! longitude in the the rotated system

! Local variables
REAL (KIND=dp)                       ::    &
  zsinpol, zcospol, zlampol, zphi, zrla, zarg1, zarg2, zrla1

REAL (KIND=dp), PARAMETER            ::    &
  zrpi18 = 57.2957795_dp,                  & !
  zpir18 = 0.0174532925_dp

!------------------------------------------------------------------------------

! Begin function rla2rlarot

  zsinpol  = SIN (zpir18 * polphi)
  zcospol  = COS (zpir18 * polphi)
  zlampol  =      zpir18 * pollam
  zphi     =      zpir18 * phi
  IF (rla > 180.0_dp) THEN
    zrla1  = rla - 360.0_dp
  ELSE
    zrla1  = rla
  ENDIF
  zrla     = zpir18 * zrla1

  zarg1    = - SIN (zrla-zlampol) * COS(zphi)
  zarg2    = - zsinpol * COS(zphi) * COS(zrla-zlampol) + zcospol * SIN(zphi)

  IF (zarg2 == 0.0) zarg2 = 1.0E-20_dp

  rla2rlarot = zrpi18 * ATAN2 (zarg1,zarg2)

  IF (polgam /= 0.0_dp ) THEN
    rla2rlarot = polgam + rla2rlarot
    IF (rla2rlarot > 180._dp) rla2rlarot = rla2rlarot -360._dp
  ENDIF

END FUNCTION RLA2RLAROT

!==============================================================================
!==============================================================================
! The following 2 subroutines have  been adopted from the COSMO model
!------------------------------------------------------------------------------

SUBROUTINE UVROT2UV_VEC(u, v, rlat, rlon, pollat, pollon, idim, jdim)

!------------------------------------------------------------------------------
!
! Description:
!   This routine converts the wind components u and v from the rotated
!   system to the real geographical system. This is the vectorized form
!   of the routine above, i.e. the computation is for a whole 2D field.
!
! Method:
!   Transformation formulas for converting between these two systems.
!
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Parameter list:
INTEGER, INTENT(IN)     ::    &
  idim, jdim        ! dimensions of the field

REAL (KIND=dp), INTENT (INOUT)       ::    &
  u  (idim,jdim), & ! wind components in the true geographical system
  v  (idim,jdim)    !

REAL (KIND=dp), INTENT (IN)          ::    &
  rlat(idim,jdim),& ! coordinates in the true geographical system
  rlon(idim,jdim),& !
  pollat, pollon    ! latitude and longitude of the north pole of the
                    ! rotated grid

! Local variables
REAL (KIND=dp)                       ::    &
  zsinpol, zcospol, zlonp, zlat, zarg1, zarg2, znorm, zugeo, zvgeo

INTEGER                              ::    i, j

! added PARAMETER attribute to enable inlining and vectorization 
! (otherwise zrpi18 / zpir18 get an implicit "SAVE" attribute which hinders inlining!)
REAL (KIND=dp), PARAMETER            ::    &
  zrpi18 = 57.2957795_dp,                  & !
  zpir18 = 0.0174532925_dp

!------------------------------------------------------------------------------
! Begin Subroutine uvrot2uv_vec
!------------------------------------------------------------------------------

! Converting from degree to radians
  zsinpol = SIN(pollat * zpir18)
  zcospol = COS(pollat * zpir18)

  DO j = 1, jdim
    DO i = 1, idim

      zlonp   = (pollon-rlon(i,j)) * zpir18
      zlat    =         rlat(i,j)  * zpir18

      zarg1   = zcospol*SIN(zlonp)
      zarg2   = zsinpol*COS(zlat) - zcospol*SIN(zlat)*COS(zlonp)
      znorm   = 1.0/SQRT(zarg1**2 + zarg2**2)

      ! Convert the u- and v-components
      zugeo   =  u(i,j)*zarg2*znorm + v(i,j)*zarg1*znorm
      zvgeo   = -u(i,j)*zarg1*znorm + v(i,j)*zarg2*znorm
      u(i,j) = zugeo
      v(i,j) = zvgeo

    ENDDO
  ENDDO

END SUBROUTINE UVROT2UV_VEC

!==============================================================================
SUBROUTINE UVROT2UV_VEC_JLOOP(u, v, rlat, rlon, pollat, pollon, idim, kdim)
! The subroutine was adopted from the subroutine utilities of the COSMO model. 
! The original code is (uvrot2uv_vec) has been adjusted to be called in j-loop.
!------------------------------------------------------------------------------
!
! Description:
!   This routine converts the wind components u and v from the rotated
!   system to the real geographical system. This is the vectorized form
!   of the routine above, i.e. the computation is for a whole 2D field.
!
! Method:
!   Transformation formulas for converting between these two systems.
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
  ! Parameter list:
!------------------------------------------------------------------------------
  ! Parameter list:
  INTEGER, INTENT(IN)     ::    &
       idim, kdim        ! dimensions of the field
  
  REAL(dp), INTENT (INOUT)       ::    &
       u  (idim, kdim), & ! wind components in the true geographical system
       v  (idim, kdim)    !
  
  REAL(dp), INTENT (IN)          ::    &
       rlat(idim),& ! coordinates in the true geographical system
       rlon(idim),& !
       pollat, pollon    ! latitude and longitude of the north pole of the
                         ! rotated grid
  
  ! Local variables
  REAL(dp)                       ::    &
       zsinpol, zcospol, zlonp, zlat, zarg1, zarg2, znorm, zugeo, zvgeo
  
  INTEGER                 ::    i, k

  ! added PARAMETER attribute to enable inlining and vectorization 
  ! (otherwise zrpi18 / zpir18 get an implicit "SAVE" attribute which hinders inlining!)
  REAL(dp), PARAMETER ::    zrpi18 = 57.2957795_dp
  REAL(dp), PARAMETER ::    zpir18 = 0.0174532925_dp
     

  !------------------------------------------------------------------------------
  ! Begin Subroutine uvrot2uv_vec_jloop
  !------------------------------------------------------------------------------
  
  ! Converting from degree to radians
  zsinpol = SIN(pollat * zpir18)
  zcospol = COS(pollat * zpir18)

  
  DO i = 1, idim
     zlonp   = (pollon-rlon(i)) * zpir18
     zlat    =         rlat(i)  * zpir18
     
     zarg1   = zcospol*SIN(zlonp)
     zarg2   = zsinpol*COS(zlat) - zcospol*SIN(zlat)*COS(zlonp)
     znorm   = 1.0/SQRT(zarg1**2 + zarg2**2)
     DO k = 1, kdim
        ! Convert the u- and v-components
        zugeo   =  u(i,k)*zarg2*znorm + v(i,k)*zarg1*znorm
        zvgeo   = -u(i,k)*zarg1*znorm + v(i,k)*zarg2*znorm
        u(i,k) = zugeo
        v(i,k) = zvgeo
     END DO
  ENDDO
     
END SUBROUTINE UVROT2UV_VEC_JLOOP
!==============================================================================

!==============================================================================
SUBROUTINE UV2UVROT_VEC(u, v, rlat, rlon, pollat, pollon, idim, jdim)

!------------------------------------------------------------------------------
!
! Description:
!   This routine converts the wind components u and v from the real
!   geographical system to the rotated system. This is the vectorized form
!   of the routine above, i.e. the computation is for a whole 2D field.
!
! Method:
!   Transformation formulas for converting between these two systems.
!
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Parameter list:
INTEGER, INTENT(IN)     ::    &
  idim, jdim        ! dimensions of the field

REAL (KIND=DP), INTENT (INOUT)       ::    &
  u  (idim,jdim), & ! wind components in the true geographical system
  v  (idim,jdim)    !

REAL (KIND=DP), INTENT (IN)          ::    &
  rlat(idim,jdim),& ! coordinates in the true geographical system
  rlon(idim,jdim),& !
  pollat, pollon    ! latitude and longitude of the north pole of the
                    ! rotated grid

! Local variables
REAL (KIND=DP)                       ::    &
  zsinpol, zcospol, zlonp, zlat, zarg1, zarg2, znorm, zurot, zvrot

INTEGER                 ::    i, j

! added PARAMETER attribute to enable inlining and vectorization 
! (otherwise zrpi18 / zpir18 get an implicit "SAVE" attribute which hinders inlining!)
REAL (KIND=DP), PARAMETER            ::    &
  zrpi18 = 57.2957795_dp,                  & !
  zpir18 = 0.0174532925_dp

!------------------------------------------------------------------------------
! Begin Subroutine uv2uvrot_vec
!------------------------------------------------------------------------------

  zsinpol = SIN ( pollat * zpir18 )
  zcospol = COS ( pollat * zpir18 )

  DO j = 1, jdim
    DO i = 1, idim

      zlonp   = ( pollon - rlon(i,j) ) * zpir18
      zlat    =            rlat(i,j)   * zpir18

      zarg1   = zcospol*SIN(zlonp)
      zarg2   = zsinpol*COS(zlat) - zcospol*SIN(zlat)*COS(zlonp)
      znorm   = 1.0_dp/SQRT( zarg1**2 + zarg2**2 )

      ! Transform the u and v wind components
      zurot =  u(i,j)*zarg2*znorm - v(i,j)*zarg1*znorm
      zvrot =  u(i,j)*zarg1*znorm + v(i,j)*zarg2*znorm
      u(i,j) = zurot
      v(i,j) = zvrot

    ENDDO
  ENDDO

END SUBROUTINE UV2UVROT_VEC
!==============================================================================

!==============================================================================
SUBROUTINE EXPAND_INPUT_GRID_LON(grid, var, lshort)

  IMPLICIT NONE
  
  ! IN-/OUTPUT
  TYPE(t_geohybgrid), INTENT(INOUT)     :: grid
  TYPE (t_ncvar), DIMENSION(:), POINTER :: var ! => NULL() ! list of variables
  LOGICAL, INTENT(IN), OPTIONAL         :: lshort
  ! LOCAL
  INTEGER                               :: status       ! error status
  CHARACTER(LEN=*), PARAMETER           :: substr = 'EXPAND_INPUT_GRID_LON'
  TYPE (t_ncvar)                        :: rvar2 
  INTEGER                               :: npdim  ! number of dims in packed
  INTEGER, DIMENSION(:), ALLOCATABLE    :: pdim   ! packed dimension vector
  INTEGER, DIMENSION(:), ALLOCATABLE    :: pvec   ! packed position vector
  INTEGER, DIMENSION(:), ALLOCATABLE    :: umdim  ! main unpacked dim. vector
  INTEGER, DIMENSION(:), POINTER        :: umvec  ! unpacked main pos. vector
  INTEGER, DIMENSION(:), ALLOCATABLE    :: shiftmap    ! map indices u on p
  INTEGER, DIMENSION(:), POINTER        :: dims     ! order of grid dimensions
  INTEGER                               :: axes(3)  ! grid axes of variable
  LOGICAL                               :: ok    = .FALSE.
  INTEGER                               :: ix, j, k
  LOGICAL                               :: llshort

  IF (PRESENT(lshort)) THEN
     llshort = lshort
  ELSE
     llshort = .FALSE.
  ENDIF

  NULLIFY(umvec)
  NULLIFY(dims)

  IF (.NOT. llshort) THEN
    CALL CHECK_NCVAR_ON_GEOHYBGRID(var(1), grid, dims, axes, ok)
    IF (.NOT.ok) CALL RGMSG(substr, RGMLE, 'VAR NOT GRID CONFORM !')
    DEALLOCATE(dims, STAT=status)
    CALL ERRMSG(substr,status,1)
    NULLIFY(dims)
 END IF

IF (.NOT. lunstruct) THEN
   grid%lonm%dim(1)%len=6
   CALL INIT_NARRAY(grid%lonm%dat,1,(/6/),VTYPE_DOUBLE)
   grid%lonm%dat%vd = (/-150._dp,-90._dp, -30._dp, 30._dp,90._dp,150._dp/) 
   grid%ranges(1,:) = (/-180._dp,180._dp/)
ELSE 
   grid%lonm%dim(1)%len=15
   CALL INIT_NARRAY(grid%lonm%dat,1,(/15/),VTYPE_DOUBLE)
   grid%lonm%dat%vd = (/-175._dp, -150._dp, -80.0_dp,-60._dp,-40._dp,-20._dp,-10._dp,-5._dp, 0._dp,5._dp,10._dp,20._dp, 40.0_dp, 75._dp, 150._dp/) 
   grid%ranges(1,:) = (/-180._dp,180._dp/)
ENDIF 

 grid%lonc = .TRUE.
 
 ! "delete" loni so that it can be recalculated with complete_geohybgrid
 CALL INIT_NCVAR(grid%loni) 
 
 IF (llshort) RETURN

 npdim   = var(1)%ndims
 ALLOCATE(pdim(npdim),shiftmap(npdim),pvec(npdim),umdim(npdim), STAT=status)
 CALL ERRMSG(substr,status,2)
 pdim = var(1)%dat%dim

 pdim(axes(GORD_LON)) = grid%lonm%dim(1)%len
 shiftmap(:) = 0
 
 DO ix=1, SIZE(var)
    CALL COPY_NCVAR(rvar2,var(ix))

    rvar2%dim(1)%len=grid%lonm%dim(1)%len
    CALL INIT_NARRAY(rvar2%dat, npdim, pdim, var(ix)%dat%type)
    
    ! CALCULATE MAPPING OF INDICES, DIMENSION VECTORS
    DO k=1, npdim ! LOOP OVER VARIABLE DIMENSIONS
       umdim(k) = var(ix)%dim(k)%len
       IF (pdim(k)>umdim(k)) shiftmap(k) = 1
    END DO
    
    ! COPY DATA
    DO k=1, PRODUCT(umdim)  ! LOOP OVER ALL ELEMENTS
       ! ELEMENT VECTOR IN UNPACKED SPACE
       CALL ELEMENT(umdim, k, umvec)
       SELECT CASE(var(ix)%dat%type)
       CASE(VTYPE_REAL)
          j = POSITION(pdim, umvec)
          rvar2%dat%vr(j) = var(ix)%dat%vr(k)
          j = POSITION(pdim, umvec+shiftmap*1)
          rvar2%dat%vr(j) = var(ix)%dat%vr(k)
          j = POSITION(pdim, umvec+shiftmap*2)
          rvar2%dat%vr(j) = var(ix)%dat%vr(k)
          j = POSITION(pdim, umvec+shiftmap*3)
          rvar2%dat%vr(j) = var(ix)%dat%vr(k)
          j = POSITION(pdim, umvec+shiftmap*4)
          rvar2%dat%vr(j) = var(ix)%dat%vr(k)
          j = POSITION(pdim, umvec+shiftmap*5)
          rvar2%dat%vr(j) = var(ix)%dat%vr(k)
          IF (grid%lonm%dim(1)%len > 6) THEN
             j = POSITION(pdim, umvec+shiftmap*6)
             rvar2%dat%vr(j) = var(ix)%dat%vr(k)
             j = POSITION(pdim, umvec+shiftmap*7)
             rvar2%dat%vr(j) = var(ix)%dat%vr(k)
             j = POSITION(pdim, umvec+shiftmap*8)
             rvar2%dat%vr(j) = var(ix)%dat%vr(k)
             j = POSITION(pdim, umvec+shiftmap*9)
             rvar2%dat%vr(j) = var(ix)%dat%vr(k)
             j = POSITION(pdim, umvec+shiftmap*10)
             rvar2%dat%vr(j) = var(ix)%dat%vr(k)
             j = POSITION(pdim, umvec+shiftmap*11)
             rvar2%dat%vr(j) = var(ix)%dat%vr(k)
             j = POSITION(pdim, umvec+shiftmap*12)
             rvar2%dat%vr(j) = var(ix)%dat%vr(k)
             j = POSITION(pdim, umvec+shiftmap*13)
             rvar2%dat%vr(j) = var(ix)%dat%vr(k)
             j = POSITION(pdim, umvec+shiftmap*14)
             rvar2%dat%vr(j) = var(ix)%dat%vr(k)
          END IF
      CASE(VTYPE_DOUBLE)
          j = POSITION(pdim, umvec)
          rvar2%dat%vd(j) = var(ix)%dat%vd(k)
          j = POSITION(pdim, umvec+shiftmap*1)
          rvar2%dat%vd(j) = var(ix)%dat%vd(k)
          j = POSITION(pdim, umvec+shiftmap*2)
          rvar2%dat%vd(j) = var(ix)%dat%vd(k)
          j = POSITION(pdim, umvec+shiftmap*3)
          rvar2%dat%vd(j) = var(ix)%dat%vd(k)
          j = POSITION(pdim, umvec+shiftmap*4)
          rvar2%dat%vd(j) = var(ix)%dat%vd(k)
          j = POSITION(pdim, umvec+shiftmap*5)
          rvar2%dat%vd(j) = var(ix)%dat%vd(k)
          IF (grid%lonm%dim(1)%len > 6) THEN
             j = POSITION(pdim, umvec+shiftmap*6)
             rvar2%dat%vd(j) = var(ix)%dat%vd(k)
             j = POSITION(pdim, umvec+shiftmap*7)
             rvar2%dat%vd(j) = var(ix)%dat%vd(k)
             j = POSITION(pdim, umvec+shiftmap*8)
             rvar2%dat%vd(j) = var(ix)%dat%vd(k)
             j = POSITION(pdim, umvec+shiftmap*9)
             rvar2%dat%vd(j) = var(ix)%dat%vd(k)
             j = POSITION(pdim, umvec+shiftmap*10)
             rvar2%dat%vd(j) = var(ix)%dat%vd(k)
             j = POSITION(pdim, umvec+shiftmap*11)
             rvar2%dat%vd(j) = var(ix)%dat%vd(k)
             j = POSITION(pdim, umvec+shiftmap*12)
             rvar2%dat%vd(j) = var(ix)%dat%vd(k)
             j = POSITION(pdim, umvec+shiftmap*13)
             rvar2%dat%vd(j) = var(ix)%dat%vd(k)
             j = POSITION(pdim, umvec+shiftmap*14)
             rvar2%dat%vd(j) = var(ix)%dat%vd(k)
          ENDIF
       CASE DEFAULT
          CALL RGMSG(substr, RGMLE, 'UNRECOGNIZED TYPE OF VARIABLE !')
       END SELECT
       !
       DEALLOCATE(umvec, STAT=status)
       CALL ERRMSG(substr,status,3)
       NULLIFY(umvec)
    END DO
    CALL COPY_NCVAR(var(ix),rvar2)
 END DO
 DEALLOCATE(pdim,shiftmap,pvec,umdim, STAT=status)
 CALL ERRMSG(substr,status,4)
 
END SUBROUTINE EXPAND_INPUT_GRID_LON
! ******************************************************************
! ******************************************************************
END MODULE MESSY_MAIN_GRID_TRAFO
! ******************************************************************
