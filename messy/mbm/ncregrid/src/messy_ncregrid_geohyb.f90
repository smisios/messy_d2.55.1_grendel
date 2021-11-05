#ifndef _SX

! ******************************************************************
! ------------------------------------------------------------------
MODULE MESSY_NCREGRID_GEOHYB
! ------------------------------------------------------------------
! Author: Patrick Joeckel, MPICH, Mainz, June 2002
! ******************************************************************
#ifdef _A2O
!------------------------------------------------------------------
! extended for curvilinear ocean grids by
! Bastian Kern, MPICH, Mainz, July 2009
!------------------------------------------------------------------
#endif
! extended for unstructured grids by
! Andreas Baumgaertner, CU Boulder / NCAR, August 2013

  USE MESSY_NCREGRID_BASE
  USE MESSY_NCREGRID_NETCDF

  IMPLICIT NONE

  INTRINSIC :: TRIM, ASSOCIATED, ABS, TINY, PRESENT, PRODUCT   &
             , SIZE, COS, MAXVAL, REAL, INT, MINVAL, IAND

  PRIVATE   :: TRIM, ASSOCIATED, ABS, TINY, PRESENT, PRODUCT   &
             , SIZE, COS, MAXVAL, REAL, INT, MINVAL, IAND

  ! INTERNAL ORDER OF LON, LAT, LEV
  ! NOTE: THIS HAS TO BE CONSISTENT WITH SEARCH ORDER IN
  !       SUBROUTINE GEOHYBGRID_AXES
  INTEGER,   PARAMETER :: GORD_LON = 1
  INTEGER,   PARAMETER :: GORD_LAT = 2
  INTEGER,   PARAMETER :: GORD_LEV = 3
  INTEGER,   PARAMETER :: GORD_COL = 1
  REAL (DP), PARAMETER :: PI = 3.141592653589_DP
  REAL,      PARAMETER :: RGEMPTY = -999.0

  TYPE geohybgrid
     CHARACTER(LEN=GRD_MAXSTRLEN) :: file    ! path/filename
     INTEGER                      :: t       ! time step
     TYPE (ncvar)                 :: lonm, latm !mid-layer
     TYPE (ncvar)                 :: loni, lati !interface
#ifdef _A2O
     ! curved grid variables
     TYPE(ncvar)                  :: clonm, clatm, cloni, clati 
#endif
     ! CESM HOMME-SE spectral element grid
     TYPE(ncvar)                  :: col ! columns
     TYPE (ncvar)                 :: hyam, hybm, timem !mid-layer
     TYPE (ncvar)                 :: hyai, hybi, timei !interface
     TYPE (ncvar)                 :: ps, p0
  END TYPE geohybgrid

!! NOTE: DOES NOT WORK PROPERLY FOR SOME COMPILERS ...
!  INTERFACE ASSIGNMENT (=)
!     MODULE PROCEDURE COPY_GEOHYBGRID
!  END INTERFACE

CONTAINS

! ------------------------------------------------------------------
SUBROUTINE INIT_GEOHYBGRID(grid)

  IMPLICIT NONE

  ! I/O
  TYPE (geohybgrid), INTENT(INOUT) :: grid

  grid%file = ''
  grid%t    = 1

  CALL INIT_NCVAR(grid%hyam)
  CALL INIT_NCVAR(grid%hybm)
  CALL INIT_NCVAR(grid%timem)
  CALL INIT_NCVAR(grid%hyai)
  CALL INIT_NCVAR(grid%hybi)
  CALL INIT_NCVAR(grid%timei)
  CALL INIT_NCVAR(grid%lonm)
  CALL INIT_NCVAR(grid%latm)
  CALL INIT_NCVAR(grid%loni)
  CALL INIT_NCVAR(grid%lati)
#ifdef _A2O
  CALL INIT_NCVAR(grid%clonm)
  CALL INIT_NCVAR(grid%clatm)
  CALL INIT_NCVAR(grid%cloni)
  CALL INIT_NCVAR(grid%clati)
#endif
  ! CESM HOMME-SE spectral element grid
  CALL INIT_NCVAR(grid%col)
  CALL INIT_NCVAR(grid%ps)
  CALL INIT_NCVAR(grid%p0)

END SUBROUTINE INIT_GEOHYBGRID
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE COPY_GEOHYBGRID(gd, gs)

  IMPLICIT NONE

  ! I/O
  TYPE (geohybgrid), INTENT(OUT) :: gd  ! destination
  TYPE (geohybgrid), INTENT(IN)  :: gs  ! source

  gd%file = TRIM(gs%file)
  gd%t    = gs%t

  CALL COPY_NCVAR(gd%hyam, gs%hyam)
  CALL COPY_NCVAR(gd%hybm, gs%hybm)
  CALL COPY_NCVAR(gd%timem, gs%timem)
  CALL COPY_NCVAR(gd%hyai, gs%hyai)
  CALL COPY_NCVAR(gd%hybi, gs%hybi)
  CALL COPY_NCVAR(gd%timei, gs%timei)
  CALL COPY_NCVAR(gd%lonm, gs%lonm)
  CALL COPY_NCVAR(gd%latm, gs%latm)
  CALL COPY_NCVAR(gd%loni, gs%loni)
  CALL COPY_NCVAR(gd%lati, gs%lati)
#ifdef _A2O
  CALL COPY_NCVAR(gd%clonm, gs%clonm)
  CALL COPY_NCVAR(gd%clatm, gs%clatm)
  CALL COPY_NCVAR(gd%cloni, gs%cloni)
  CALL COPY_NCVAR(gd%clati, gs%clati)
#endif
  CALL COPY_NCVAR(gd%col, gs%col)

  CALL COPY_NCVAR(gd%ps, gs%ps)
  CALL COPY_NCVAR(gd%p0, gs%p0)

END SUBROUTINE COPY_GEOHYBGRID
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE IMPORT_GEOHYBGRID(grid)

  ! NO CHECKING IN THIS ROUTINE !!!

  IMPLICIT NONE

  ! I/O
  TYPE (geohybgrid), INTENT(INOUT) :: grid

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'IMPORT_GEOHYBGRID'
  CHARACTER(LEN=GRD_MAXSTRLEN) :: name   ! variable name
  CHARACTER(LEN=GRD_MAXSTRLEN) :: file   ! netCDF filename
  INTEGER                      :: ustep  ! step along UNLIMITED DIM.
  INTEGER                      :: status
  INTEGER                      :: iostat
  INTEGER                      :: ndims
  INTEGER                      :: dimlen ! dimension length of ps
  INTEGER, DIMENSION(:), ALLOCATABLE :: dimvec ! dimension vector of ps
  REAL                         :: p      ! local pressure
  CHARACTER(LEN=GRD_MAXSTRLEN) :: unit   ! pressure unit
  INTEGER                      :: setuid ! force unlimited ID

  ustep = grid%t
  file = TRIM(grid%file)

  ! LONGITUDE
  IF (TRIM(grid%lonm%name) /= '') THEN
     name = TRIM(grid%lonm%name)
     CALL IMPORT_NCVAR(grid%lonm, varname=name, file=file)
  END IF
  IF (TRIM(grid%loni%name) /= '') THEN
     name = TRIM(grid%loni%name)
     CALL IMPORT_NCVAR(grid%loni, varname=name, file=file)
  END IF
#ifdef _A2O
  IF (TRIM(grid%clonm%name) /= '') THEN
     name = TRIM(grid%clonm%name)
     CALL IMPORT_NCVAR(grid%clonm, varname=name, file=file)
  END IF
  IF (TRIM(grid%cloni%name) /= '') THEN
     name = TRIM(grid%cloni%name)
     CALL IMPORT_NCVAR(grid%cloni, varname=name, file=file)
  END IF
#endif

  ! LATITUDE
  IF (TRIM(grid%latm%name) /= '') THEN
     name = TRIM(grid%latm%name)
     CALL IMPORT_NCVAR(grid%latm, varname=name, file=file)
  END IF
  IF (TRIM(grid%lati%name) /= '') THEN
     name = TRIM(grid%lati%name)
     CALL IMPORT_NCVAR(grid%lati, varname=name, file=file)
  END IF
#ifdef _A2O
  IF (TRIM(grid%clatm%name) /= '') THEN
     name = TRIM(grid%clatm%name)
     CALL IMPORT_NCVAR(grid%clatm, varname=name, file=file)
  END IF
  IF (TRIM(grid%clati%name) /= '') THEN
     name = TRIM(grid%clati%name)
     CALL IMPORT_NCVAR(grid%clati, varname=name, file=file)
  END IF
#endif
  IF (TRIM(grid%col%name) /= '') THEN
     name = TRIM(grid%col%name)
     CALL IMPORT_NCVAR(grid%col, varname=name, file=file)
  END IF

  ! TIME
  IF (TRIM(grid%timem%name) /= '') THEN
     name = TRIM(grid%timem%name)
     CALL IMPORT_NCVAR(grid%timem, ustep=ustep, varname=name, file=file)
     ! ASSOCIATE ALL TIME AXES WITH UNLIMITED ID
     grid%timem%dim(1)%fuid = .true.
     grid%timem%uid = grid%timem%dim(1)%id
  END IF
  IF (TRIM(grid%timei%name) /= '') THEN
     name = TRIM(grid%timei%name)
     CALL IMPORT_NCVAR(grid%timei, varname=name, file=file)
  END IF

  ! HYA
  IF (TRIM(grid%hyam%name) /= '') THEN
     name = TRIM(grid%hyam%name)
     CALL IMPORT_NCVAR(grid%hyam, varname=name, file=file)
  END IF
  IF (TRIM(grid%hyai%name) /= '') THEN
     name = TRIM(grid%hyai%name)
     CALL IMPORT_NCVAR(grid%hyai, varname=name, file=file)
  END IF

  ! HYB
  IF (TRIM(grid%hybm%name) /= '') THEN
     name = TRIM(grid%hybm%name)
     CALL IMPORT_NCVAR(grid%hybm, varname=name, file=file)
  END IF
  IF (TRIM(grid%hybi%name) /= '') THEN
     name = TRIM(grid%hybi%name)
     CALL IMPORT_NCVAR(grid%hybi, varname=name, file=file)
  END IF

! mz_ab_20130830
  ! SURFACE PRESSURE
  IF (TRIM(grid%ps%name) /= '') THEN
     ! CHECK FOR REAL VALUE
     READ(grid%ps%name,*,IOSTAT=iostat) p
     IF (iostat == 0) THEN    ! VALUE of PS
        CALL INIT_NCVAR(grid%ps)
        grid%ps%name  = 'ps'
        grid%ps%xtype = NF90_FLOAT
        !
        ! DIMS
        ndims = 0
        IF (grid%lonm%dat%type /= VTYPE_UNDEF) ndims=ndims+1
        IF (grid%latm%dat%type /= VTYPE_UNDEF) ndims=ndims+1
        grid%ps%ndims = ndims
        ALLOCATE(grid%ps%dim(grid%ps%ndims), STAT=status)
        CALL ERRMSG(substr,status,1)
        ALLOCATE(dimvec(grid%ps%ndims), STAT=status)
        CALL ERRMSG(substr,status,2)
        ndims = 0
        dimlen = 1
        IF (grid%lonm%dat%type /= VTYPE_UNDEF) THEN
           ndims = ndims + 1
           CALL COPY_NCDIM(grid%ps%dim(ndims), grid%lonm%dim(1))
           dimlen = dimlen * grid%lonm%dim(1)%len
           dimvec(ndims) = grid%lonm%dim(1)%len
        END IF
        IF (grid%latm%dat%type /= VTYPE_UNDEF) THEN
           ndims = ndims + 1
           CALL COPY_NCDIM(grid%ps%dim(ndims), grid%latm%dim(1))
           dimlen = dimlen * grid%latm%dim(1)%len
           dimvec(ndims) = grid%latm%dim(1)%len
        END IF
        !
        ! ATTRIBUTES ...
        ! ... LONGNAME ATTRIBUTE
        CALL ADD_NCATT(grid%ps, 'long_name'            &
                        ,vs='constant surface pressure')
        ! ... UNITS ATTRIBUT
        unit = ''
        READ(grid%ps%name,*,IOSTAT=iostat) p,unit
        IF (TRIM(unit) /= '') THEN
           CALL ADD_NCATT(grid%ps, 'units' ,vs=TRIM(unit))
        END IF
        !
        ! DATA
        CALL INIT_NARRAY(grid%ps%dat, ndims, dimvec, VTYPE_REAL)
        grid%ps%dat%vr(:) = p
        DEALLOCATE(dimvec, STAT=status)
        CALL ERRMSG(substr,status,3)
     ELSE  ! PS IS VARIABLE NAME
        name = TRIM(grid%ps%name)
        IF (ASSOCIATED(grid%timem%dim)) THEN
           setuid = grid%timem%dim(1)%id
        ELSE
           setuid = NULL_DIMID
        END IF
        CALL IMPORT_NCVAR(grid%ps, ustep=ustep, varname=name, file=file &
             , setuid=setuid)
     END IF ! PS IS VALUE OR NAME
  END IF

  IF (TRIM(grid%p0%name) /= '') THEN
     ! CHECK FOR REAL VALUE
     READ(grid%p0%name,*,IOSTAT=iostat) p
     IF (iostat == 0) THEN    ! VALUE of P0
        CALL INIT_NCVAR(grid%p0)
        grid%p0%name  = 'p0'
        grid%p0%xtype = NF90_FLOAT
        !
        ! NO DIMENSIONS !!!
        grid%p0%ndims = 1
        ALLOCATE(grid%p0%dim(1), STAT=status)
        CALL ERRMSG(substr,status,4)
        grid%p0%dim(1)%name = TRIM(grid%p0%name)//'_dim'
        grid%p0%dim(1)%id   = NULL_DIMID
        grid%p0%dim(1)%len  = 1
        grid%p0%dim(1)%fuid = .false.
        grid%p0%dim(1)%varid = NULL_VARID
        !
        ! ATTRIBUTES ...
        ! ... LONGNAME ATTRIBUTE
        CALL ADD_NCATT(grid%p0, 'long_name'            &
                        ,vs='reference pressure')
        ! ... UNITS ATTRIBUT
        unit = ''
        READ(grid%p0%name,*,IOSTAT=iostat) p,unit
        IF (TRIM(unit) /= '') THEN
           CALL ADD_NCATT(grid%p0, 'units' ,vs=TRIM(unit))
        END IF
        !
        ! DATA
        CALL INIT_NARRAY(grid%p0%dat, 1, (/ 1 /), VTYPE_REAL)
        grid%p0%dat%vr(1) = p
     ELSE  ! NAME of P0
        name = TRIM(grid%p0%name)
        IF (ASSOCIATED(grid%timem%dim)) THEN
           setuid = grid%timem%dim(1)%id
        ELSE
           setuid = NULL_DIMID
        END IF
        CALL IMPORT_NCVAR(grid%p0, varname=name, file=file, setuid=setuid)
     END IF
  END IF

END SUBROUTINE IMPORT_GEOHYBGRID
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE EXPORT_GEOHYBGRID(grid)

  ! NO CHECKING IN THIS ROUTINE !!!

  IMPLICIT NONE

  ! I/O
  TYPE (geohybgrid), INTENT(INOUT) :: grid

  ! LOCAL
  CHARACTER(LEN=GRD_MAXSTRLEN) :: file   ! netCDF filename
  INTEGER                      :: ustep  ! step along UNLIMITED DIM.

  ustep = grid%t
  file = TRIM(grid%file)

  ! LONGITUDE
  IF (TRIM(grid%lonm%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%lonm, file=file)
  END IF
  IF (TRIM(grid%loni%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%loni, file=file)
  END IF
#ifdef _A2O
  IF (TRIM(grid%clonm%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%clonm, file=file)
  END IF
  IF (TRIM(grid%cloni%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%cloni, file=file)
  END IF
#endif

  ! LATITUDE
  IF (TRIM(grid%latm%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%latm, file=file)
  END IF
  IF (TRIM(grid%lati%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%lati, file=file)
  END IF
#ifdef _A2O
  IF (TRIM(grid%clatm%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%clatm, file=file)
  END IF
  IF (TRIM(grid%clati%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%clati, file=file)
  END IF
#endif
  ! COLUMNS
  IF (TRIM(grid%col%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%col, file=file)
  END IF

  ! TIME
  IF (TRIM(grid%timem%name) /= '') THEN
     grid%timem%ustep = ustep
     CALL EXPORT_NCVAR(grid%timem, file=file)
  END IF
!
  ! IF 'COMPLETED', timei has also the UNLIMITED-ATTRIBUTE
  ! DO NOT EXPORT AT THE MOMENT
!  IF (TRIM(grid%timei%name) /= '') THEN
!     CALL EXPORT_NCVAR(grid%timei, file=file)
!  END IF

  ! HYA
  IF (TRIM(grid%hyam%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%hyam, file=file)
  END IF
  IF (TRIM(grid%hyai%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%hyai, file=file)
  END IF

  ! HYB
  IF (TRIM(grid%hybm%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%hybm, file=file)
  END IF
  IF (TRIM(grid%hybi%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%hybi, file=file)
  END IF

  ! SURFACE PRESSURE
  IF (TRIM(grid%ps%name) /= '') THEN
     grid%ps%ustep = ustep
     CALL EXPORT_NCVAR(grid%ps, file=file)
  END IF

  ! REFERENCE PRESSURE
  IF (TRIM(grid%p0%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%p0, file=file)
  END IF

END SUBROUTINE EXPORT_GEOHYBGRID
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE SWITCH_GEOHYBGRID(g, lx, ly, lz)

  IMPLICIT NONE

  ! I/O
  TYPE(geohybgrid), INTENT(INOUT) :: g
  LOGICAL         , INTENT(IN)    :: lx, ly, lz

  IF (.NOT.lx) THEN
     CALL INIT_NCVAR(g%loni)
     CALL INIT_NCVAR(g%lonm)
  END IF

  IF (.NOT.ly) THEN
     CALL INIT_NCVAR(g%lati)
     CALL INIT_NCVAR(g%latm)
  END IF

  IF (.NOT.lx) THEN
     CALL INIT_NCVAR(g%col)
  END IF

  IF (.NOT.lz) THEN
     CALL INIT_NCVAR(g%hyai)
     CALL INIT_NCVAR(g%hyam)
     CALL INIT_NCVAR(g%hybi)
     CALL INIT_NCVAR(g%hybm)
     CALL INIT_NCVAR(g%p0)
     CALL INIT_NCVAR(g%ps)
  END IF

END SUBROUTINE SWITCH_GEOHYBGRID
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE CHECK_GEOHYBGRID(grid, ranges)

  IMPLICIT NONE

  ! I/O
  TYPE (geohybgrid),        INTENT(INOUT) :: grid
  REAL(DP), DIMENSION(4,2), INTENT(IN)    :: ranges

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'CHECK_GEOHYBGRID'

  ! LONGITUDE
  IF (grid%lonm%ndims > 1) THEN
     CALL RGMSG(substr, RGMLE, &
          'LONGITUDE VARIABLE '''//TRIM(grid%lonm%name)//'''', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IN FILE '''//TRIM(grid%file)//'''', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IS ',grid%lonm%ndims,'-DIMENSIONAL !')
  END IF

  IF (grid%loni%ndims > 1) THEN
     CALL RGMSG(substr, RGMLE, &
          'LONGITUDE INTERFACES VARIABLE '''//TRIM(grid%loni%name)//'''',&
          .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IN FILE '''//TRIM(grid%file)//'''', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IS ',grid%loni%ndims,'-DIMENSIONAL !')
  END IF

  ! LATITUDE
  IF (grid%latm%ndims > 1) THEN
     CALL RGMSG(substr, RGMLE, &
          'LATITUDE VARIABLE '''//TRIM(grid%latm%name)//'''', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IN FILE '''//TRIM(grid%file)//'''', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IS ',grid%latm%ndims,'-DIMENSIONAL !')
  END IF

  IF (grid%lati%ndims > 1) THEN
     CALL RGMSG(substr, RGMLE, &
          'LATITUDE INTERFACES VARIABLE '''//TRIM(grid%lati%name)//'''',&
          .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IN FILE '''//TRIM(grid%file)//'''', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IS ',grid%lati%ndims,'-DIMENSIONAL !')
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
!        IF ((ranges(1,1) == RGEMPTY).OR.(ranges(1,2) == RGEMPTY)) THEN
        IF (grid%loni%ndims /= 1) THEN
           IF ((ABS(ranges(1,1) - RGEMPTY) <= TINY(RGEMPTY)).OR. &
                (ABS(ranges(1,2) - RGEMPTY) <= TINY(RGEMPTY))) THEN
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
!        IF ((ranges(2,1) == RGEMPTY).OR.(ranges(2,2) == RGEMPTY)) THEN
        IF (grid%lati%ndims /= 1) THEN
           IF ((ABS(ranges(2,1) - RGEMPTY) <= TINY(RGEMPTY)).OR. &
                (ABS(ranges(2,2) - RGEMPTY) <= TINY(RGEMPTY))) THEN
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
!        IF ((ranges(3,1) == RGEMPTY).OR.(ranges(3,2) == RGEMPTY)) THEN
        IF (grid%hyai%ndims /= 1) THEN
           IF ((ABS(ranges(3,1) - RGEMPTY) <= TINY(RGEMPTY)).OR. &
                (ABS(ranges(3,2) - RGEMPTY) <= TINY(RGEMPTY))) THEN
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
!        IF ((ranges(4,1) == RGEMPTY).OR.(ranges(4,2) == RGEMPTY)) THEN
        IF (grid%hybi%ndims /= 1) THEN
           IF ((ABS(ranges(4,1) - RGEMPTY) <= TINY(RGEMPTY)).OR. &
                (ABS(ranges(4,2) - RGEMPTY) <= TINY(RGEMPTY))) THEN
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
  TYPE (geohybgrid), INTENT(IN)    :: gi  ! INPUT GRID
  TYPE (geohybgrid), INTENT(OUT)   :: go  ! OUTPUT GRID
  TYPE (geohybgrid), INTENT(INOUT) :: gx  ! INDEX 'GRID'
  LOGICAL, OPTIONAL, INTENT(IN)    :: reverse

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'SORT_GEOHYBGRID'
  LOGICAL :: lrev
  INTEGER :: i
  INTEGER :: olat, olon  ! order of dimensions in PS
  INTEGER, DIMENSION(:), ALLOCATABLE :: dil  ! dimension length vector of PS
  INTEGER, DIMENSION(:), POINTER     :: vec  ! element vector of PS
  INTEGER, DIMENSION(:), ALLOCATABLE :: svec ! element vector of sorted PS
  INTEGER                            :: vtype
  INTEGER                            :: status
  TYPE(narray)                       :: ts   ! temporal n-array for sigma level
  TYPE(narray)                       :: col  ! sigma level column
  TYPE(narray)                       :: cx   ! sigma level column sort indices

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
              svec(olon) = gx%lonm%dat%vi(vec(olon))
           END IF
           IF (olat > 0) THEN
              svec(olat) = gx%latm%dat%vi(vec(olat))
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
              svec(olon) = gx%lonm%dat%vi(vec(olon))
           END IF
           IF (olat > 0) THEN
              svec(olat) = gx%latm%dat%vi(vec(olat))
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
    TYPE (ncvar), INTENT(INOUT) :: var

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'IDX_NCVAR'
    INTEGER :: i
    INTEGER :: status
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
SUBROUTINE GEOHYBGRID_AXES(g, a, g2, a2, pflag)

  ! Note: NO CHECKING
  !       g, g2 must be 'ordered', 'complete', and 'consistent'

  IMPLICIT NONE

  ! I/O
  TYPE (geohybgrid), INTENT(IN)      :: g  ! GEO-HYBRID-GRID
  TYPE (axis), DIMENSION(:), POINTER :: a  ! LIST OF AXES FOR REGRIDDING
  TYPE (geohybgrid), INTENT(IN)     , OPTIONAL :: g2  ! GEO-HYBRID-GRID
  TYPE (axis), DIMENSION(:), POINTER, OPTIONAL :: a2  ! LIST OF AXES
  LOGICAL, OPTIONAL, INTENT(IN)      :: pflag  ! .true.: pressure axis
                                               ! .false. sigma-axis (default)

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'GEOHYBGRID_AXES'
  INTEGER                        :: n      ! number of axes
  INTEGER                        :: status
  LOGICAL                        :: lp     ! local presure axis flag
  INTEGER                        :: i      ! counter
  INTEGER                        :: vtype
  INTEGER                        :: ndep_lat ! number of lat-dimension in list
  INTEGER                        :: ndep_lon ! number of lon-dimension in list
  INTEGER                        :: ndep_col ! number of col-dimension in list
  INTEGER                        :: dpc      ! dependency counter
  TYPE (narray)                  :: zps      ! local surface pressure

  ! CHECK SUBROUTINE CALL
  IF (PRESENT(g2).NEQV.PRESENT(a2)) THEN
     CALL RGMSG(substr, RGMLE, &
          '2nd GRID AND 2nd AXES LIST MUST BE PRESENT ON SUBROUTINE CALL !')
  END IF

  ! INIT
  IF (PRESENT(pflag)) THEN
     lp = pflag
  ELSE
     lp = .false.     ! DEFAULT: sigma coordinates
  END IF
  ndep_lon = 0
  ndep_lat = 0
  ndep_col = 0

  ! COUNT AXES
  n = 1             ! START WITH MINIMUM 1 AXIS (LAST AXIS IS 'FREE' DIM)

  IF (QDEF_NCVAR(g%loni)) n=n+1
  IF (QDEF_NCVAR(g%lati)) n=n+1
  IF (QDEF_NCVAR(g%col)) n=n+1
  IF (QDEF_NCVAR(g%hyai).OR.QDEF_NCVAR(g%hybi)) n=n+1

  ! ALLOCATE SPACE
  ALLOCATE(a(n), STAT=status)
  CALL ERRMSG(substr,status,1)
  DO i=1, n
     CALL INIT_AXIS(a(i))
  END DO

  IF (PRESENT(g2)) THEN
     ALLOCATE(a2(n), STAT=status)
     CALL ERRMSG(substr,status,2)
     DO i=1, n
        CALL INIT_AXIS(a2(i))
     END DO
  END IF

  ! ASSIGN DATA
  n = 0

  CALL RGMSG(substr, RGMLI, 'AXIS CONSTRUCTION ...')

  ! NOTE: THE ORDER OF TESTING THE FOLLOWING DIMENSIONS
  !       (HERE: LON, LAT, LEV) SHOULD NOT BE CHANGED !
  !       IT HAS TO BE CONSISTENT WITH OTHER ROUTINES

  ! LONGITUDE
  IF (QDEF_NCVAR(g%loni)) THEN
     CALL RGMSG(substr, RGMLIC, ' ... ADDING LONGITUDE AXIS ...')
     n=n+1
     a(n)%lm     = .true.     ! LONGITUDE IS MODULO AXIS
     a(n)%ndp    = 1          ! LONGITUDE IS ...
     ALLOCATE(a(n)%dep(1), STAT=status)
     CALL ERRMSG(substr,status,3)
     a(n)%dep(1) = n          ! ... INDEPENDENT
     CALL COPY_NARRAY(a(n)%dat, g%loni%dat)
     ndep_lon    = n
     IF (PRESENT(g2)) THEN
        IF (QDEF_NCVAR(g2%loni)) THEN
           CALL RGMSG(substr, RGMLIC, ' ... DIFFERENT FOR OUTPUT GRID ...')
           a2(n)%lm  = .true.
           a2(n)%ndp = 1
           ALLOCATE(a2(n)%dep(1), STAT=status)
           CALL ERRMSG(substr,status,4)
           a2(n)%dep(1) = n
           CALL COPY_NARRAY(a2(n)%dat, g2%loni%dat)
        ELSE !UNDEFINED => INVARIANT
           CALL RGMSG(substr, RGMLIC, ' ... INVARIANT ...')
        END IF
     END IF
  END IF

  ! LATITUDE
  IF (QDEF_NCVAR(g%lati)) THEN
     CALL RGMSG(substr, RGMLIC, ' ... ADDING LATITUDE AXIS ...')
     n=n+1
     a(n)%lm     = .false.    ! LATITUDE IS NON-MODULO AXIS
     a(n)%ndp    = 1          ! LATITUDE IS ...
     ALLOCATE(a(n)%dep(1), STAT=status)
     CALL ERRMSG(substr,status,5)
     a(n)%dep(1) = n          ! ... INDEPENDENT
     CALL COPY_NARRAY(a(n)%dat, g%lati%dat)
     ndep_lat    = n
     IF (PRESENT(g2)) THEN
        IF (QDEF_NCVAR(g2%lati)) THEN
           CALL RGMSG(substr, RGMLIC, ' ... DIFFERENT FOR OUTPUT GRID ...')
           a2(n)%lm  = .false.
           a2(n)%ndp = 1
           ALLOCATE(a2(n)%dep(1), STAT=status)
           CALL ERRMSG(substr,status,6)
           a2(n)%dep(1) = n
           CALL COPY_NARRAY(a2(n)%dat, g2%lati%dat)
        ELSE !UNDEFINED => INVARIANT
           CALL RGMSG(substr, RGMLIC, ' ... INVARIANT ...')
        END IF
     END IF
     !
     ! TAKE INTO ACCOUNT SPHERICAL GEOMETRY ...
     vtype = a(n)%dat%type
     SELECT CASE(vtype)
     CASE(VTYPE_REAL)
        a(n)%dat%vr = COS(((a(n)%dat%vr - REAL(90., SP))/ &
             REAL(180., SP))*REAL(PI,SP))
     CASE(VTYPE_DOUBLE)
        a(n)%dat%vd = COS(((a(n)%dat%vd - REAL(90., DP))/ &
             REAL(180., DP))*PI)
     CASE(VTYPE_INT)
        CALL RGMSG(substr, RGMLE, &
             'LATITUDE AXIS OF TYPE INTEGER IS NOT SUPPORTED !')
     CASE(VTYPE_BYTE)
        CALL RGMSG(substr, RGMLE, &
             'LATITUDE AXIS OF TYPE BYTE IS NOT SUPPORTED !')
     CASE(VTYPE_CHAR)
        CALL RGMSG(substr, RGMLE, &
             'LATITUDE AXIS OF TYPE CHAR IS NOT SUPPORTED !')
     CASE(VTYPE_UNDEF)
        CALL RGMSG(substr, RGMLE, &
             'LATITUDE AXIS IS UNDEFINED !')
     CASE DEFAULT
        CALL RGMSG(substr, RGMLE, &
             'TYPE OF LATITUDE AXIS IS NOT RECOGNIZED !')
     END SELECT
     !
     ! ... ALSO FOR 2nd GRID
     IF (PRESENT(g2)) THEN
        IF (QDEF_NCVAR(g2%lati)) THEN
           vtype = a2(n)%dat%type
           SELECT CASE(vtype)
           CASE(VTYPE_REAL)
              a2(n)%dat%vr = COS(((a2(n)%dat%vr - REAL(90., SP))/ &
                   REAL(180., SP))*REAL(PI,SP))
           CASE(VTYPE_DOUBLE)
              a2(n)%dat%vd = COS(((a2(n)%dat%vd - REAL(90., DP))/ &
                   REAL(180., DP))*PI)
           CASE(VTYPE_INT)
              CALL RGMSG(substr, RGMLE, &
                   'LATITUDE AXIS OF TYPE INTEGER IS NOT SUPPORTED !')
           CASE(VTYPE_BYTE)
              CALL RGMSG(substr, RGMLE, &
                   'LATITUDE AXIS OF TYPE BYTE IS NOT SUPPORTED !')
           CASE(VTYPE_CHAR)
              CALL RGMSG(substr, RGMLE, &
                   'LATITUDE AXIS OF TYPE CHAR IS NOT SUPPORTED !')
           CASE(VTYPE_UNDEF)
              CALL RGMSG(substr, RGMLE, &
                   'LATITUDE AXIS IS UNDEFINED !')
           CASE DEFAULT
              CALL RGMSG(substr, RGMLE, &
                   'TYPE OF LATITUDE AXIS IS NOT RECOGNIZED !')
           END SELECT
        END IF
     END IF
  END IF

  ! COLUMNS
  IF (QDEF_NCVAR(g%col)) THEN
     CALL RGMSG(substr, RGMLIC, ' ... ADDING COLUMNS PSEUDO-AXIS ...')
     n=n+1
     a(n)%lm     = .false.   
     a(n)%ndp    = 1         
     ALLOCATE(a(n)%dep(1), STAT=status)
     CALL ERRMSG(substr,status,3)
     a(n)%dep(1) = n          ! ... INDEPENDENT
     CALL COPY_NARRAY(a(n)%dat, g%col%dat)
     !ndep_lon    = n
  END IF

  ! LEVELS
  IF (QDEF_NCVAR(g%hyai).OR.QDEF_NCVAR(g%hybi)) THEN
     CALL RGMSG(substr, RGMLIC, ' ... ADDING VERTICAL AXIS ...')
     n=n+1
     a(n)%lm     = .false.     ! VERTICAL AXIS IS NON-MODULO AXIS

     ! CALCULATE AXIS DATA
     IF (lp) THEN
        CALL RGMSG(substr, RGMLIC, '     -> PRESSURE AXIS ...')
     ELSE
        CALL RGMSG(substr, RGMLIC, '     -> DIMENSIONLESS AXIS ...')
     END IF
     CALL PS2PS(g%ps, zps)
     CALL H2PSIG(a(n)%dat,g%hyai%dat,g%hybi%dat,zps,g%p0%dat,lp)
     CALL INIT_NARRAY(zps)

     ! SET DEPENDENCIES
     a(n)%ndp = a(n)%dat%n
     ALLOCATE(a(n)%dep(a(n)%ndp), STAT=status)
     CALL ERRMSG(substr,status,7)
     a(n)%dep(:) = 0
     dpc = 1
     a(n)%dep(dpc) = n
     !
     IF (a(n)%ndp > 1) THEN
        DO i=1, g%ps%ndims
           ! CHECK DIM LENGTH AND NAME/ID
           IF (QCMP_NCDIM(g%ps%dim(i), g%lonm%dim(1)) > 1) THEN ! LONGITUDE
              dpc = dpc + 1
              a(n)%dep(dpc) = ndep_lon
              CALL RGMSG(substr, RGMLIC, '     -> DEPENDING ON LONGITUDE ...')
           END IF
           IF (QCMP_NCDIM(g%ps%dim(i), g%latm%dim(1)) > 1) THEN ! LATITUDE
              dpc = dpc + 1
              a(n)%dep(dpc) = ndep_lat
              CALL RGMSG(substr, RGMLIC, '     -> DEPENDING ON LATITUDE ...')
           END IF
           IF (QDEF_NCVAR(g%col)) THEN ! op_pj_20170120
              IF (QCMP_NCDIM(g%ps%dim(i), g%col%dim(1)) > 1) THEN ! COLUMN
                 dpc = dpc + 1
                 a(n)%dep(dpc) = ndep_col
                 CALL RGMSG(substr, RGMLIC, '     -> DEPENDING ON COLUMN ...')
              END IF
           END IF                      ! op_pj_20170120
        END DO
        IF ((dpc > 1).AND.(dpc /= a(n)%ndp)) THEN
           CALL RGMSG(substr, RGMLE, 'DEPENDENCY MISMATCH !')
        END IF
     END IF

     ! 2nd GRID
     IF (PRESENT(g2)) THEN
        IF (QDEF_NCVAR(g2%hyai).OR.QDEF_NCVAR(g2%hybi)) THEN
           CALL RGMSG(substr, RGMLIC, ' ... DIFFERENT FOR OUTPUT GRID ...')
           a2(n)%lm     = .false.     ! VERTICAL AXIS IS NON-MODULO AXIS

           ! CALCULATE AXIS DATA
           CALL PS2PS(g2%ps, zps)
           CALL H2PSIG(a2(n)%dat,g2%hyai%dat,g2%hybi%dat, &
                zps,g2%p0%dat,lp)
           CALL INIT_NARRAY(zps)

           ! SET DEPENDENCIES
           a2(n)%ndp = a2(n)%dat%n
           ALLOCATE(a2(n)%dep(a2(n)%ndp), STAT=status)
           CALL ERRMSG(substr,status,8)
           a2(n)%dep(:) = 0
           dpc = 1
           a2(n)%dep(dpc) = n
           !
           IF (a2(n)%ndp > 1) THEN
              DO i=1, g2%ps%ndims
                 ! CHECK DIM LENGTH AND NAME/ID
                 IF (QDEF_NCVAR(g2%lonm)) THEN
                    IF (QCMP_NCDIM(g2%ps%dim(i), g2%lonm%dim(1)) > 1) THEN
                       ! LONGITUDE
                       dpc = dpc + 1
                       a2(n)%dep(dpc) = ndep_lon
                       CALL RGMSG(substr, RGMLIC, &
                            '     -> DEPENDING ON LONGITUDE ...')
                    END IF
                 END IF
                 IF (QDEF_NCVAR(g2%latm)) THEN
                    IF (QCMP_NCDIM(g2%ps%dim(i), g2%latm%dim(1)) > 1) THEN
                       ! LATITUDE
                       dpc = dpc + 1
                       a2(n)%dep(dpc) = ndep_lat
                       CALL RGMSG(substr, RGMLIC, &
                            '     -> DEPENDING ON LATITUDE ...')
                    END IF
                 END IF
                 IF (QDEF_NCVAR(g2%col)) THEN
                    IF (QCMP_NCDIM(g2%ps%dim(i), g2%col%dim(1)) > 1) THEN
                       ! COLUMN
                       dpc = dpc + 1
                       a2(n)%dep(dpc) = ndep_col
                       CALL RGMSG(substr, RGMLIC, &
                            '     -> DEPENDING ON COLUMN ...')
                    END IF
                 END IF
              END DO
              IF ((dpc > 1).AND.(dpc /= a2(n)%ndp)) THEN
                 CALL RGMSG(substr, RGMLE, 'DEPENDENCY MISMATCH !')
              END IF
           END IF
        ELSE ! UNDEFINED => INVARIANT
           CALL RGMSG(substr, RGMLIC, ' ... INVARIANT ...')
        END IF
     END IF ! 2nd GRID

  END IF ! LEVELS

  CALL RGMSG(substr, RGMLIC, '... END AXES CONSTRUCTION !')

END SUBROUTINE GEOHYBGRID_AXES
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE PS2PS(var, na)

  IMPLICIT NONE

  ! I/O
  TYPE (ncvar) , INTENT(IN)  :: var
  TYPE (narray), INTENT(OUT) :: na

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'PS2PS'
  INTEGER :: i
  INTEGER :: uidpos
  INTEGER, DIMENSION(:), POINTER :: vec
  INTEGER :: dc
  INTEGER :: status
  INTEGER :: vtype

  NULLIFY(vec)

  uidpos = 0
  DO i=1, var%ndims
     IF (var%dim(i)%fuid) THEN
        uidpos = i
        EXIT
     END IF
  END DO

  IF (uidpos == 0) THEN  ! no unlimited ID
     CALL COPY_NARRAY(na, var%dat)
  ELSE                   ! 'remove' unlimited ID
     IF ((var%dim(uidpos)%len /= 1).OR.(var%dat%dim(uidpos) /= 1)) THEN
        CALL RGMSG(substr, RGMLE, &
             'DIMENSION LENGTH OF UNLIMITED DIMENSION MUST BE 1 !')
     END IF
     ALLOCATE(vec(var%ndims-1), STAT=status)
     CALL ERRMSG(substr,status,1)
     dc = 0
     DO i=1, var%ndims
        IF (.NOT.var%dim(i)%fuid) THEN
           dc = dc + 1
           vec(dc) = var%dim(i)%len
        END IF
     END DO
     vtype = var%dat%type
     CALL INIT_NARRAY(na, var%ndims-1, vec, vtype)
     SELECT CASE(vtype)
     CASE(VTYPE_REAL)
        na%vr(:) = var%dat%vr(:)
     CASE(VTYPE_DOUBLE)
        na%vd(:) = var%dat%vd(:)
     CASE(VTYPE_INT)
        na%vi(:) = var%dat%vi(:)
     CASE(VTYPE_BYTE)
        na%vb(:) = var%dat%vb(:)
     CASE(VTYPE_CHAR)
        na%vc(:) = var%dat%vc(:)
     CASE(VTYPE_UNDEF)
        CALL RGMSG(substr, RGMLE, 'N-ARRAY IS UNDEFINED !')
     CASE DEFAULT
        CALL RGMSG(substr, RGMLE, 'N-ARRAY IS UNRECOGNIZED !')
     END SELECT
     DEALLOCATE(vec,STAT=status)
     CALL ERRMSG(substr,status,2)
     NULLIFY(vec)
  END IF

END SUBROUTINE PS2PS
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE H2PSIG(psig, hya, hyb, ps, p0, lp)

  IMPLICIT NONE

  ! I/O
  TYPE (narray), INTENT(INOUT) :: psig
  TYPE (narray), INTENT(IN)    :: hya, hyb, ps, p0
  LOGICAL      , INTENT(IN)    :: lp   ! .true.  -> pressure axis
                                       ! .false. -> dimensionless axis
  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'H2PSIG'
  INTEGER :: dim   ! dimension length of linear array in psig
  INTEGER :: status
  INTEGER :: i
  INTEGER :: i1,i2 ! little helpers
  TYPE (narray) :: zhya, zhyb, zps, zp0     ! local arrays
  INTEGER, DIMENSION(:), POINTER :: vec     ! element vector

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
SUBROUTINE COMPLETE_GEOHYBGRID(g, ranges, gx)

  ! Note: g must be already sorted

  IMPLICIT NONE

  ! I/O
  TYPE (geohybgrid),        INTENT(INOUT)           :: g
  REAL(DP), DIMENSION(4,2), INTENT(IN)              :: ranges
  TYPE (geohybgrid),        INTENT(INOUT), OPTIONAL :: gx   ! SORT INDICES

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'COMPLETE_GEOHYBGRID'

  CALL RGMSG(substr, RGMLI, 'FILE '''//TRIM(g%file)//'''')
  CALL RGMSG(substr, RGMLIC,'STEP ',g%t,' ')

  CALL RGMSG(substr, RGMLIC,'  loni <<--->> lonm ...')
  CALL IMMI_NARRAY(g%loni%dat, g%lonm%dat)
  CALL IMMI_NCVAR(g%loni, g%lonm)
  IF (PRESENT(gx)) THEN
     CALL IMMI_NARRAY_IDX(gx%loni%dat, gx%lonm%dat)
     CALL IMMI_NCVAR(gx%loni, gx%lonm)
  END IF
  ! CORRECT FOR RANGE
  CALL RNGADJ_NARRAY(g%loni%dat, ranges(1,:))

  CALL RGMSG(substr, RGMLIC,'  lati <<--->> latm ...')
  CALL IMMI_NARRAY(g%lati%dat, g%latm%dat)
  CALL IMMI_NCVAR(g%lati, g%latm)
  IF (PRESENT(gx)) THEN
     CALL IMMI_NARRAY_IDX(gx%lati%dat, gx%latm%dat)
     CALL IMMI_NCVAR(gx%lati, gx%latm)
  END IF
  ! CORRECT FOR RANGE
  CALL RNGADJ_NARRAY(g%lati%dat, ranges(2,:))

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

! --------------------------------------------
SUBROUTINE IMMI_NARRAY(nai, nam)

  IMPLICIT NONE

  ! I/O
  TYPE (narray), INTENT(INOUT) :: nai, nam

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'IMMI_NARRAY'
  INTEGER :: itype, mtype
  INTEGER :: i

  itype = nai%type
  mtype = nam%type

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
        CALL RGMSG(substr, RGMLE, &
             'INTERFACE ARRAY IS ',nai%n,'-DIMENSIONAL !')
     END IF
     IF (nam%n /= 1) THEN
        CALL RGMSG(substr, RGMLE, &
             'MID ARRAY IS ',nam%n,'-DIMENSIONAL !')
     END IF
     IF (nai%dim(1) /= (nam%dim(1)+1)) THEN
        CALL RGMSG(substr, RGMLE, &
             'SIZES OF INTERFACE ARRAY AND MID ARRAY ARE INCOMPATIBLE !')
     END IF
     CALL RGMSG(substr, RGMLIC, '  ... O.K.')
     RETURN
  END IF

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
SUBROUTINE IMMI_NARRAY_IDX(nai, nam)

  IMPLICIT NONE

  ! I/O
  TYPE (narray), INTENT(INOUT) :: nai, nam

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'IMMI_NARRAY_IDX'
  INTEGER :: itype, mtype
  INTEGER :: i
  INTEGER :: status

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
        i = nai%vi(nai%dim(1)-1)
        nai%vi(nai%dim(1)-1) = nai%vi(nai%dim(1))
        nai%vi(nai%dim(1)) = i
     END IF
     RETURN
  END IF

  CALL RGMSG(substr, RGMLE, &
       'N-ARRAYS MUST BE OF TYPE INTEGER OR UNDEFINED !')

END SUBROUTINE IMMI_NARRAY_IDX
! --------------------------------------------

! --------------------------------------------
  SUBROUTINE IMMI_NCVAR(vari, varm)

    IMPLICIT NONE

    ! I/O
    TYPE (ncvar), INTENT(INOUT) :: vari, varm

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'IMMI_NCVAR'
    INTEGER :: i
    INTEGER :: status
    CHARACTER(LEN=GRD_MAXSTRLEN) :: str
    INTEGER :: vtype

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
             !vari%dim(i)%fuid  = .false.  ! UNLIM-DIM IS ONLY ALONG MIDs !
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
  SUBROUTINE RNGADJ_NARRAY(na, mr, lmon)

    IMPLICIT NONE

    ! I/O
    TYPE (narray), INTENT(INOUT)        :: na
    REAL(DP),      INTENT(IN)           :: mr(2)
    LOGICAL,       INTENT(IN), OPTIONAL :: lmon

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'RNGADJ_NARRAY'
    INTEGER :: n
    INTEGER :: vtype
    LOGICAL :: llmon

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

END SUBROUTINE COMPLETE_GEOHYBGRID
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE CHECK_NCVAR_ON_GEOHYBGRID(var, g, dims, axes, ok)

  IMPLICIT NONE

  ! I/O
  TYPE (ncvar)     , INTENT(IN)  :: var    ! variable
  TYPE (geohybgrid), INTENT(IN)  :: g      ! grid
  INTEGER, DIMENSION(:), POINTER :: dims   ! order of g-dims in var
  INTEGER          , INTENT(OUT) :: axes(3)! dimension no. of lon -> lat -> lev
  LOGICAL          , INTENT(OUT) :: ok     ! conform ?

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'CHECK_NCVAR_ON_GEOHYBGRID'
  INTEGER :: i
  INTEGER :: ndimg            ! number of dimensions in g
  INTEGER :: ndimv            ! number of g-dims in var
  INTEGER :: status
  INTEGER :: plon, plat, plev ! position of lon, lat, lev in AXES-LIST !!!
  INTEGER :: pcol             ! mz_ab_20130823
  LOGICAL :: lvert            ! vertical axis present ?
  LOGICAL :: llat, llon       ! lat, lon axes present ?
  LOGICAL :: lcol ! mz_ab_20130823

  ! INIT
  ok = .true.
  axes(:) = 0
  plon = 0
  plat = 0
  plev = 0
  pcol = 0 ! mz_ab_20130823

  ! ALLOCATE SPACE FOR 'ORDER OF DIMENSIONS'
  ALLOCATE(dims(var%ndims), STAT=status)
  CALL ERRMSG(substr,status,1)
  dims(:) = 0

  ! GET NUMBER OF HYBRID-DIMENSIONS IN GRID
  ndimg = 0
  IF (QDEF_NCVAR(g%lonm)) THEN
     ndimg = ndimg + 1
     plon = ndimg
  END IF
  IF (QDEF_NCVAR(g%latm)) THEN
     ndimg = ndimg + 1
     plat = ndimg
  END IF
  IF (QDEF_NCVAR(g%col)) THEN
     ndimg = ndimg + 1
     pcol = ndimg
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
  lcol = .false. ! mz_ab_20130823

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
           IF (llon) THEN
              ndimv = ndimv + 1
              dims(i) = plon
              axes(GORD_LON) = i
              CYCLE
           END IF
        END IF
     END IF

     IF (.NOT.llat) THEN  ! CHECK LAT
        IF (ASSOCIATED(g%latm%dim)) THEN
#ifndef _A2O
           llat = (QCMP_NCDIM(g%latm%dim(1), var%dim(i)) > 1)
           IF (ASSOCIATED(g%lati%dim).AND.(.NOT.llat)) THEN
              llat = ((QCMP_NCDIM(g%latm%dim(1), var%dim(i)) > 0).AND. &
                      (TRIM(g%latm%dim(1)%name) ==                     &
                       TRIM(g%lati%dim(1)%name)//'_M') )
           END IF
#else
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
#endif
           IF (llat) THEN
              ndimv = ndimv + 1
              dims(i) = plat
              axes(GORD_LAT) = i
              CYCLE
           END IF
        END IF
     END IF

     IF (.NOT.lcol) THEN  ! CHECK COLUMNS
        IF (ASSOCIATED(g%col%dim)) THEN
           lcol = (QCMP_NCDIM(g%col%dim(1), var%dim(i)) > 1)
           IF (lcol) THEN
              ndimv = ndimv + 1
              dims(i) = pcol
              axes(GORD_COL) = i
              CYCLE
           END IF
        END IF
     END IF

     IF (.NOT.lvert) THEN  ! CHECK HYA
        IF (ASSOCIATED(g%hyam%dim)) THEN
           lvert = (QCMP_NCDIM(g%hyam%dim(1), var%dim(i)) > 1)
           IF (ASSOCIATED(g%hyai%dim).AND.(.NOT.lvert)) THEN
#ifndef _A2O
              lvert = ((QCMP_NCDIM(g%hyam%dim(1), var%dim(i)) > 0).AND. &
                       (TRIM(g%hyam%dim(1)%name) ==                     &
                        TRIM(g%hyai%dim(1)%name)//'_M') )
#else
              lvert = (QCMP_NCDIM(g%hyai%dim(1), var%dim(i)) > 1)
#endif
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
           IF (ASSOCIATED(g%hybi%dim).AND.(.NOT.lvert)) THEN
#ifndef _A2O
              lvert = ((QCMP_NCDIM(g%hybm%dim(1), var%dim(i)) > 0).AND. &
                       (TRIM(g%hybm%dim(1)%name) ==                     &
                        TRIM(g%hybi%dim(1)%name)//'_M') )
#else
              lvert = (QCMP_NCDIM(g%hybi%dim(1), var%dim(i)) > 1)
#endif
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

#ifndef _A2O
  ok = (ndimv >= ndimg)  ! ALL GRID DIMS MUST BE RECOGNIZED !!!
                         ! REMAINING DIMS ARE 'FREE'
#else
  ok = (ndimv >= ndimg) .OR. ((.NOT. lvert) .AND. (ndimv + 1 >= ndimg))
  ! ALL GRID DIMS MUST BE RECOGNIZED !!!
  ! REMAINING DIMS ARE 'FREE'
  ! 2D fields on 3D grid also possible
#endif

END SUBROUTINE CHECK_NCVAR_ON_GEOHYBGRID
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE SORT_GEOHYBGRID_NCVAR(var, gx, axes, svar, reverse)

  IMPLICIT NONE

  ! I/O
  TYPE (ncvar)     , INTENT(IN)           :: var     ! input variable
  TYPE (geohybgrid), INTENT(IN)           :: gx      ! hybrid grid with  ...
                                                     ! ... index information
  INTEGER          , INTENT(IN)           :: axes(3) ! lon,lat,lev dim. no.
  TYPE (ncvar)     , INTENT(OUT)          :: svar    ! sorted variable
  LOGICAL          , INTENT(IN), OPTIONAL :: reverse

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'SORT_GEOHYBGRID_NCVAR'
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
! op_pj_20130422+
!!$        vec(axes(GORD_LON)) = gx%lonm%dat%vi(vec(GORD_LON))
        vec(axes(GORD_LON)) = gx%lonm%dat%vi(vec(axes(GORD_LON)))
! op_pj_20130422-
     END IF
     IF (axes(GORD_LAT) > 0) THEN
! op_pj_20130422+
!!$        vec(axes(GORD_LAT)) = gx%latm%dat%vi(vec(GORD_LAT))
        vec(axes(GORD_LAT)) = gx%latm%dat%vi(vec(axes(GORD_LAT)))
! op_pj_20130422-
     END IF
     IF (axes(GORD_LEV) > 0) THEN
        IF (gx%hyam%dat%type /= VTYPE_UNDEF) THEN
! op_pj_20130422+
!!$           vec(axes(GORD_LEV)) = gx%hyam%dat%vi(vec(GORD_LEV))
           vec(axes(GORD_LEV)) = gx%hyam%dat%vi(vec(axes(GORD_LEV)))
! op_pj_20130422-
        ELSE
! op_pj_20130422+
!!$           vec(axes(GORD_LEV)) = gx%hybm%dat%vi(vec(GORD_LEV))
           vec(axes(GORD_LEV)) = gx%hybm%dat%vi(vec(axes(GORD_LEV)))
! op_pj_20130422-
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
  TYPE (ncvar)                  , INTENT(IN)    :: vi   ! input variable
  INTEGER, DIMENSION(:)         , INTENT(IN)    :: dims
  INTEGER                       , INTENT(IN)    :: axes(3)
  TYPE (ncvar)                  , INTENT(INOUT) :: vo   ! output variable
  LOGICAL, OPTIONAL             , INTENT(IN)    :: reverse

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'PACK_GEOHYBGRID_NCVAR'
  TYPE (ncvar)                       :: vu     ! unpacked variable
  TYPE (ncvar)                       :: vp     ! packed variable
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
  TYPE (geohybgrid), INTENT(INOUT) :: gi      ! input grid
  TYPE (geohybgrid), INTENT(INOUT) :: go      ! output grid

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

#ifdef _A2O
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
#endif

  IF (QDEF_NCVAR(gi%col).AND.(.NOT.QDEF_NCVAR(go%col))) THEN
     CALL COPY_NCVAR(go%col, gi%col)
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
SUBROUTINE BALANCE_GEOHYBGRID_TIME(gi, go, lint)

  IMPLICIT NONE

  ! I/O
  TYPE (geohybgrid),      INTENT(INOUT) :: gi, go
  LOGICAL,                INTENT(IN)    :: lint

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'BALANCE_GEOHYBGRID_TIME'
  INTEGER :: i
  INTEGER :: uii, uio

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

! ------------------------------------------------------------------
SUBROUTINE BALANCE_GEOHYBGRID_PS(gi, go, ranges)

  IMPLICIT NONE

  ! I/O
  TYPE (geohybgrid),          INTENT(INOUT) :: gi, go
  REAL(DP), DIMENSION(2,4,2), INTENT(IN)    :: ranges

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'BALANCE_GEOHYBGRID_PS'
  LOGICAL :: err
  LOGICAL :: linv   ! PRE-REGRID AVAILABLE go%ps

  err = (.NOT.QDEF_NCVAR(gi%ps)).AND.(.NOT.QDEF_NCVAR(go%ps))
  err = err .AND. ((QDEF_NCVAR(gi%hybi).OR.QDEF_NCVAR(gi%hybm)) .OR. &
                   (QDEF_NCVAR(go%hybi).OR.QDEF_NCVAR(go%hybm)))

  IF (err) THEN
     CALL RGMSG(substr, RGMLE, &
          'HYBRID-B-COEFFICIENTS NEED I_PS AND/OR G_PS IN NAMELIST !')
  END IF

  err = (.NOT.QDEF_NCVAR(gi%p0)).AND.(.NOT.QDEF_NCVAR(go%p0))
  err = err .AND. ((QDEF_NCVAR(gi%hyai).OR.QDEF_NCVAR(gi%hyam)) .OR. &
                   (QDEF_NCVAR(go%hyai).OR.QDEF_NCVAR(go%hyam)))

  IF (err) THEN
     CALL RGMSG(substr, RGMLE, &
          'HYBRID-A-COEFFICIENTS NEED I_P0 AND/OR G_P0 IN NAMELIST !')
  END IF

  ! RETURN, IF 2D
  err = (.NOT.QDEF_NCVAR(gi%p0)).AND.(.NOT.QDEF_NCVAR(gi%ps)).AND.     &
        (.NOT.QDEF_NCVAR(gi%hyai)).AND.(.NOT.QDEF_NCVAR(gi%hybi)).AND. &
        (.NOT.QDEF_NCVAR(gi%hyam)).AND.(.NOT.QDEF_NCVAR(gi%hybm))
  IF (err) THEN
     CALL RGMSG(substr, RGMLI, &
          'INPUT GRID IS 2-D! NO SURFACE PRESSURE REGRIDDING REQUIRED!')
     RETURN
  END IF

  ! BALANCE PO IN ANY CASE
  IF (QDEF_NCVAR(gi%p0).AND.(.NOT.QDEF_NCVAR(go%p0))) THEN
     CALL COPY_NCVAR(go%p0, gi%p0)
  ELSE
     IF (QDEF_NCVAR(go%p0).AND.(.NOT.QDEF_NCVAR(gi%p0))) THEN
        CALL COPY_NCVAR(gi%p0, go%p0)
     END IF
  END IF

  ! NOW PS ...

  ! CASE A: gi%ps AND go%ps BOTH AVAILABLE
  ! ADJUST TIME AXIS
  ! NOTE: THE TIME BALANCING (INCLUDING THAT FOR THE SURFACE PRESSURE
  !       TIME DIMENSION) IS PERFORMED IN
  !       SUBROUTINE BALANCE_GEOHYBGRID_TIME !

  ! CASE B: gi%ps AND go%ps BOTH UNAVAILABLE
  ! NOTHING TO DO !
  IF (.NOT.QDEF_NCVAR(gi%ps).AND.(.NOT.QDEF_NCVAR(go%ps))) THEN
     RETURN
  END IF

  ! SURFACE PRESSURE NEEDS TO BE PRE-REGRIDDED, IF NOT AVAILABLE
  ! CASE C: go%ps AVAILABLE, BUT ON WRONG HORIZONTAL GRID
  !
  linv = QDEF_NCVAR(go%ps).AND.( &
         ((.NOT.QDEF_NCVAR(go%lonm)).AND.(.NOT.QDEF_NCVAR(go%loni))).OR. &
         ((.NOT.QDEF_NCVAR(go%latm)).AND.(.NOT.QDEF_NCVAR(go%lati))) )
  !
  IF (QDEF_NCVAR(go%ps).AND.linv) THEN
     CALL RGMSG(substr, RGMLE,  &
          'REGRIDDING 3-D DISTRIBUTIONS', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'ONTO A DESTINATION SURFACE PRESSURE COORDINATE', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'USING (AN) INVARIANT HORIZONTAL DIMENSION(S)', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IS NOT POSSIBLE DUE TO A LACK OF INFORMATION', .false.)
     CALL RGMSG(substr, RGMLEC, &
          '(HORIZONTAL DESTINATION GRID)', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'FOR PRE-REGRIDDING THE DESTINATION SURFACE PRESSURE !', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'PLEASE PERFORM 2-D PRE-REGRIDDING OF SURFACE PRESSURE', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IN SEPARATE STEP!')
  END IF

  ! CASE D: gi%ps XOR go%ps NOT AVAILABLE
  IF (QDEF_NCVAR(gi%ps).AND.(.NOT.QDEF_NCVAR(go%ps))) THEN
     CALL REGRID_GEOHYBGRID_PS(gi, go, ranges(2,:,:))
  ELSE
     IF (QDEF_NCVAR(go%ps).AND.(.NOT.QDEF_NCVAR(gi%ps))) THEN
        CALL REGRID_GEOHYBGRID_PS(go, gi, ranges(1,:,:))
     END IF
  END IF

END SUBROUTINE BALANCE_GEOHYBGRID_PS
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE REGRID_GEOHYBGRID_PS(gi, go, ranges)

  IMPLICIT NONE

  ! I/O
  TYPE (geohybgrid),        INTENT(IN)    :: gi
  TYPE (geohybgrid),        INTENT(INOUT) :: go
  REAL(DP), DIMENSION(4,2), INTENT(IN)    :: ranges

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'REGRID_GEOHYBGRID_PS'
  TYPE (geohybgrid) :: gih, goh, gihs, gohs
  TYPE (geohybgrid) :: gix, gox
  INTEGER, DIMENSION(:),       POINTER :: dims      ! order of grid dimensions
  INTEGER                              :: axes(3)   ! grid axes of variable
  TYPE (axis)  , DIMENSION(:), POINTER :: ai, ao    ! source and dest. axes
  TYPE (ncvar)                         :: psi, pso
  TYPE (ncvar)                         :: psis, psos
  TYPE (ncvar)                         :: psip, psop
  TYPE (narray), DIMENSION(:), POINTER :: nao       ! re-gridder I/O
  REAL (DP),     DIMENSION(:), POINTER :: sovl, dovl
  INTEGER,       DIMENSION(:), POINTER :: rcnt
  LOGICAL                              :: ok
  INTEGER                              :: i
  INTEGER                              :: status

  ! INIT
  NULLIFY(dims)
  NULLIFY(ai)
  NULLIFY(ao)
  NULLIFY(nao)
  NULLIFY(sovl)
  NULLIFY(dovl)
  NULLIFY(rcnt)

  ! CREATE 2-D GRIDS
  CALL INIT_GEOHYBGRID(gih)
  CALL INIT_GEOHYBGRID(goh)
  !
  gih%file = TRIM(gi%file)
  gih%t    = gi%t
  !
  goh%file = TRIM(go%file)
  goh%t    = go%t
  !
  CALL COPY_NCVAR(gih%lonm, gi%lonm)
  CALL COPY_NCVAR(gih%loni, gi%loni)
  !
  CALL COPY_NCVAR(gih%latm, gi%latm)
  CALL COPY_NCVAR(gih%lati, gi%lati)
  CALL COPY_NCVAR(gih%col, gi%col)
  !
  CALL COPY_NCVAR(gih%timem, gi%timem)
  CALL COPY_NCVAR(gih%timei, gi%timei)
  !
  CALL COPY_NCVAR(goh%lonm, go%lonm)
  CALL COPY_NCVAR(goh%loni, go%loni)
  !
  CALL COPY_NCVAR(goh%latm, go%latm)
  CALL COPY_NCVAR(goh%lati, go%lati)
  CALL COPY_NCVAR(goh%col, go%col)
  !
  CALL COPY_NCVAR(goh%timem, go%timem)
  CALL COPY_NCVAR(goh%timei, go%timei)

  CALL SORT_GEOHYBGRID(gih, gihs, gix)
  CALL SORT_GEOHYBGRID(goh, gohs, gox)

  CALL COMPLETE_GEOHYBGRID(gihs, ranges, gix)
  CALL COMPLETE_GEOHYBGRID(gohs, ranges, gox)

  CALL BALANCE_GEOHYBGRID(gihs, gohs)
  CALL BALANCE_GEOHYBGRID(gix, gox)

  CALL GEOHYBGRID_AXES(gihs, ai, gohs, ao, .false.)

  CALL CHECK_NCVAR_ON_GEOHYBGRID(gi%ps, gihs, dims, axes, ok)
  IF (.NOT.ok) THEN
     CALL RGMSG(substr, RGMLE, 'PS NOT GRID CONFORM !')
  END IF

  CALL COPY_NCVAR(psi, gi%ps)
  CALL SORT_GEOHYBGRID_NCVAR(psi, gix, axes, psis)
  CALL BALANCE_GEOHYBGRID_NCVAR(psis, axes, gohs, psos)

  CALL PACK_GEOHYBGRID_NCVAR(psis, dims, axes, psip)
  CALL BALANCE_GEOHYBGRID_NCVAR(psip, axes, gohs, psop)

  DEALLOCATE(dims, STAT=status)
  CALL ERRMSG(substr,status,1)
  NULLIFY(dims)

  CALL NREGRID( (/ psip%dat /), ai, ao, nao, (/ RG_INT /), sovl, dovl, rcnt)
  IF (IAND(MSGMODE, MSGMODE_VM) == MSGMODE_VM) THEN
     CALL NREGRID_STAT(ai, ao, sovl, dovl, (/ psip%dat /), nao, rcnt)
  END IF

  CALL COPY_NARRAY(psop%dat, nao(1))
  CALL CHECK_NCVAR_ON_GEOHYBGRID(psos, gohs, dims, axes, ok)

  IF (.NOT.ok) THEN
     CALL RGMSG(substr, RGMLE, 'PS NOT GRID CONFORM !')
  END IF

  CALL PACK_GEOHYBGRID_NCVAR(psop, dims, axes, psos, .true.)

  CALL SORT_GEOHYBGRID_NCVAR(psos, gox, axes, pso, .true.)

  CALL COPY_NCVAR(go%ps, pso)

  ! CLEAN UP
  DEALLOCATE(dims, STAT=status)
  CALL ERRMSG(substr,status,2)
  NULLIFY(dims)
  DEALLOCATE(sovl, dovl, rcnt, STAT=status)
  CALL ERRMSG(substr,status,3)
  NULLIFY(sovl, dovl, rcnt)
  CALL INIT_NCVAR(psi)
  CALL INIT_NCVAR(pso)
  CALL INIT_NCVAR(psis)
  CALL INIT_NCVAR(psos)
  CALL INIT_NCVAR(psip)
  CALL INIT_NCVAR(psop)
  CALL INIT_GEOHYBGRID(gih)
  CALL INIT_GEOHYBGRID(goh)
  CALL INIT_GEOHYBGRID(gihs)
  CALL INIT_GEOHYBGRID(gohs)
  CALL INIT_GEOHYBGRID(gix)
  CALL INIT_GEOHYBGRID(gox)
  DO i=1, SIZE(nao)
     CALL INIT_NARRAY(nao(i))
  END DO
  DEALLOCATE(nao, STAT=status)
  CALL ERRMSG(substr,status,4)
  NULLIFY(nao)
  DO i=1, SIZE(ai)
     CALL INIT_AXIS(ai(i))
     CALL INIT_AXIS(ao(i))
  END DO
  DEALLOCATE(ai, ao, STAT=status)
  CALL ERRMSG(substr,status,5)
  NULLIFY(ai,ao)

END SUBROUTINE REGRID_GEOHYBGRID_PS
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE BALANCE_GEOHYBGRID_NCVAR(vari, axes, go, varo)

  IMPLICIT NONE

  ! Note: go must already be 'balanced'

  ! I/O
  TYPE (ncvar)     , INTENT(IN)  :: vari    ! input variable
  INTEGER          , INTENT(IN)  :: axes(3) ! dim.no of lon, lat, lev
  TYPE (geohybgrid), INTENT(IN)  :: go      ! output grid
  TYPE (ncvar)     , INTENT(OUT) :: varo    ! output variable

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'BALANCE_GEOHYBGRID_NCVAR'
  INTEGER                            :: i
  INTEGER                            :: vtype
  INTEGER, DIMENSION(:), ALLOCATABLE :: dimvec  ! dimension vector
  INTEGER                            :: status

  ! INIT
  !
  CALL COPY_NCVAR(varo, vari)
  !
  ALLOCATE(dimvec(vari%ndims), STAT=status)
  CALL ERRMSG(substr,status,1)
  !
  vtype = vari%dat%type

  ! CHANGE DIMENSIONS
  DO i=1, vari%ndims
     dimvec(i) = vari%dim(i)%len
     IF (i == axes(GORD_LON)) THEN
        dimvec(i) = go%lonm%dim(1)%len
        varo%dim(i) = go%lonm%dim(1)
     END IF
#ifndef _A2O
     IF (i == axes(GORD_LAT)) THEN
        dimvec(i) = go%latm%dim(1)%len
        varo%dim(i) = go%latm%dim(1)
     END IF
#else
     IF (i == axes(GORD_LAT)) THEN
        IF (go%latm%ndims == 2) THEN
           dimvec(i) = go%latm%dim(2)%len
           varo%dim(i) = go%latm%dim(2)
        ELSE
           dimvec(i) = go%latm%dim(1)%len
           varo%dim(i) = go%latm%dim(1)
        END IF
     END IF
#endif

     IF (i == axes(GORD_COL)) THEN
        IF (go%col%dat%type /= VTYPE_UNDEF) THEN ! op_pj_20170120
           dimvec(i) = go%col%dim(1)%len
           varo%dim(i) = go%col%dim(1)
        END IF                                   ! op_pj_20170120
     END IF

     IF (i == axes(GORD_LEV)) THEN
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

END SUBROUTINE BALANCE_GEOHYBGRID_NCVAR
! ------------------------------------------------------------------

! ******************************************************************
END MODULE MESSY_NCREGRID_GEOHYB
! ******************************************************************

#else

! ******************************************************************
! ------------------------------------------------------------------
MODULE MESSY_NCREGRID_GEOHYB
! ------------------------------------------------------------------
! Author: Patrick Joeckel, MPICH, Mainz, June 2002
! ******************************************************************

  USE MESSY_NCREGRID_BASE
  USE MESSY_NCREGRID_NETCDF

  IMPLICIT NONE

  INTRINSIC :: TRIM, ASSOCIATED, ABS, TINY, PRESENT, PRODUCT   &
             , SIZE, COS, MAXVAL, REAL, INT, MINVAL, IAND

  PRIVATE   :: TRIM, ASSOCIATED, ABS, TINY, PRESENT, PRODUCT   &
             , SIZE, COS, MAXVAL, REAL, INT, MINVAL, IAND

  ! INTERNAL ORDER OF LON, LAT, LEV
  ! NOTE: THIS HAS TO BE CONSISTENT WITH SEARCH ORDER IN
  !       SUBROUTINE GEOHYBGRID_AXES
  INTEGER,   PARAMETER :: GORD_LON = 1
  INTEGER,   PARAMETER :: GORD_LAT = 2
  INTEGER,   PARAMETER :: GORD_LEV = 3
  REAL (DP), PARAMETER :: PI = 3.141592653589_DP
  REAL,      PARAMETER :: RGEMPTY = -999.0

  TYPE geohybgrid
     CHARACTER(LEN=GRD_MAXSTRLEN)        :: file    ! path/filename
     INTEGER                             :: t       ! time step
     TYPE (ncvar)                        :: lonm, latm, hyam, hybm, timem
     TYPE (ncvar)                        :: loni, lati, hyai, hybi, timei
     TYPE (ncvar)                        :: ps, p0
  END TYPE geohybgrid

!! NOTE: DOES NOT WORK PROPERLY FOR SOME COMPILERS ...
!  INTERFACE ASSIGNMENT (=)
!     MODULE PROCEDURE COPY_GEOHYBGRID
!  END INTERFACE

CONTAINS

! ------------------------------------------------------------------
SUBROUTINE INIT_GEOHYBGRID(grid)

  IMPLICIT NONE

  ! I/O
  TYPE (geohybgrid), INTENT(INOUT) :: grid

  grid%file = ''
  grid%t    = 1

  CALL INIT_NCVAR(grid%lonm)
  CALL INIT_NCVAR(grid%latm)
  CALL INIT_NCVAR(grid%hyam)
  CALL INIT_NCVAR(grid%hybm)
  CALL INIT_NCVAR(grid%timem)
  CALL INIT_NCVAR(grid%loni)
  CALL INIT_NCVAR(grid%lati)
  CALL INIT_NCVAR(grid%hyai)
  CALL INIT_NCVAR(grid%hybi)
  CALL INIT_NCVAR(grid%timei)

  CALL INIT_NCVAR(grid%ps)
  CALL INIT_NCVAR(grid%p0)

END SUBROUTINE INIT_GEOHYBGRID
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE COPY_GEOHYBGRID(gd, gs)

  IMPLICIT NONE

  ! I/O
  TYPE (geohybgrid), INTENT(OUT) :: gd  ! destination
  TYPE (geohybgrid), INTENT(IN)  :: gs  ! source

  gd%file = TRIM(gs%file)
  gd%t    = gs%t

  CALL COPY_NCVAR(gd%lonm, gs%lonm)
  CALL COPY_NCVAR(gd%latm, gs%latm)
  CALL COPY_NCVAR(gd%hyam, gs%hyam)
  CALL COPY_NCVAR(gd%hybm, gs%hybm)
  CALL COPY_NCVAR(gd%timem, gs%timem)
  CALL COPY_NCVAR(gd%loni, gs%loni)
  CALL COPY_NCVAR(gd%lati, gs%lati)
  CALL COPY_NCVAR(gd%hyai, gs%hyai)
  CALL COPY_NCVAR(gd%hybi, gs%hybi)
  CALL COPY_NCVAR(gd%timei, gs%timei)

  CALL COPY_NCVAR(gd%ps, gs%ps)
  CALL COPY_NCVAR(gd%p0, gs%p0)

END SUBROUTINE COPY_GEOHYBGRID
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE IMPORT_GEOHYBGRID(grid)

  ! NO CHECKING IN THIS ROUTINE !!!

  IMPLICIT NONE

  ! I/O
  TYPE (geohybgrid), INTENT(INOUT) :: grid

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'IMPORT_GEOHYBGRID'
  CHARACTER(LEN=GRD_MAXSTRLEN) :: name   ! variable name
  CHARACTER(LEN=GRD_MAXSTRLEN) :: file   ! netCDF filename
  INTEGER                      :: ustep  ! step along UNLIMITED DIM.
  INTEGER                      :: status
  INTEGER                      :: iostat
  INTEGER                      :: ndims
  INTEGER                      :: dimlen ! dimension length of ps
  INTEGER, DIMENSION(:), ALLOCATABLE :: dimvec ! dimension vector of ps
  REAL                         :: p      ! local pressure
  CHARACTER(LEN=GRD_MAXSTRLEN) :: unit   ! pressure unit
  INTEGER                      :: setuid ! force unlimited ID

  ustep = grid%t
  file = TRIM(grid%file)

  ! LONGITUDE
  IF (TRIM(grid%lonm%name) /= '') THEN
     name = TRIM(grid%lonm%name)
     CALL IMPORT_NCVAR(grid%lonm, varname=name, file=file)
  END IF
  IF (TRIM(grid%loni%name) /= '') THEN
     name = TRIM(grid%loni%name)
     CALL IMPORT_NCVAR(grid%loni, varname=name, file=file)
  END IF

  ! LATITUDE
  IF (TRIM(grid%latm%name) /= '') THEN
     name = TRIM(grid%latm%name)
     CALL IMPORT_NCVAR(grid%latm, varname=name, file=file)
  END IF
  IF (TRIM(grid%lati%name) /= '') THEN
     name = TRIM(grid%lati%name)
     CALL IMPORT_NCVAR(grid%lati, varname=name, file=file)
  END IF

  ! TIME
  IF (TRIM(grid%timem%name) /= '') THEN
     name = TRIM(grid%timem%name)
     CALL IMPORT_NCVAR(grid%timem, ustep=ustep, varname=name, file=file)
     ! ASSOCIATE ALL TIME AXES WITH UNLIMITED ID
     grid%timem%dim(1)%fuid = .true.
     grid%timem%uid = grid%timem%dim(1)%id
  END IF
  IF (TRIM(grid%timei%name) /= '') THEN
     name = TRIM(grid%timei%name)
     CALL IMPORT_NCVAR(grid%timei, varname=name, file=file)
  END IF

  ! HYA
  IF (TRIM(grid%hyam%name) /= '') THEN
     name = TRIM(grid%hyam%name)
     CALL IMPORT_NCVAR(grid%hyam, varname=name, file=file)
  END IF
  IF (TRIM(grid%hyai%name) /= '') THEN
     name = TRIM(grid%hyai%name)
     CALL IMPORT_NCVAR(grid%hyai, varname=name, file=file)
  END IF

  ! HYB
  IF (TRIM(grid%hybm%name) /= '') THEN
     name = TRIM(grid%hybm%name)
     CALL IMPORT_NCVAR(grid%hybm, varname=name, file=file)
  END IF
  IF (TRIM(grid%hybi%name) /= '') THEN
     name = TRIM(grid%hybi%name)
     CALL IMPORT_NCVAR(grid%hybi, varname=name, file=file)
  END IF

  ! SURFACE PRESSURE
  IF (TRIM(grid%ps%name) /= '') THEN
     ! CHECK FOR REAL VALUE
     READ(grid%ps%name,*,IOSTAT=iostat) p
     IF (iostat == 0) THEN    ! VALUE of PS
        CALL INIT_NCVAR(grid%ps)
        grid%ps%name  = 'ps'
        grid%ps%xtype = NF90_FLOAT
        !
        ! DIMS
        ndims = 0
        IF (grid%lonm%dat%type /= VTYPE_UNDEF) ndims=ndims+1
        IF (grid%latm%dat%type /= VTYPE_UNDEF) ndims=ndims+1
        grid%ps%ndims = ndims
        ALLOCATE(grid%ps%dim(grid%ps%ndims), STAT=status)
        CALL ERRMSG(substr,status,1)
        ALLOCATE(dimvec(grid%ps%ndims), STAT=status)
        CALL ERRMSG(substr,status,2)
        ndims = 0
        dimlen = 1
        IF (grid%lonm%dat%type /= VTYPE_UNDEF) THEN
           ndims = ndims + 1
           CALL COPY_NCDIM(grid%ps%dim(ndims), grid%lonm%dim(1))
           dimlen = dimlen * grid%lonm%dim(1)%len
           dimvec(ndims) = grid%lonm%dim(1)%len
        END IF
        IF (grid%latm%dat%type /= VTYPE_UNDEF) THEN
           ndims = ndims + 1
           CALL COPY_NCDIM(grid%ps%dim(ndims), grid%latm%dim(1))
           dimlen = dimlen * grid%latm%dim(1)%len
           dimvec(ndims) = grid%latm%dim(1)%len
        END IF
        !
        ! ATTRIBUTES ...
        ! ... LONGNAME ATTRIBUTE
        CALL ADD_NCATT(grid%ps, 'long_name'            &
                        ,vs='constant surface pressure')
        ! ... UNITS ATTRIBUT
        unit = ''
        READ(grid%ps%name,*,IOSTAT=iostat) p,unit
        IF (TRIM(unit) /= '') THEN
           CALL ADD_NCATT(grid%ps, 'units' ,vs=TRIM(unit))
        END IF
        !
        ! DATA
        CALL INIT_NARRAY(grid%ps%dat, ndims, dimvec, VTYPE_REAL)
        grid%ps%dat%vr(:) = p
        DEALLOCATE(dimvec, STAT=status)
        CALL ERRMSG(substr,status,3)
     ELSE  ! PS IS VARIABLE NAME
        name = TRIM(grid%ps%name)
        IF (ASSOCIATED(grid%timem%dim)) THEN
           setuid = grid%timem%dim(1)%id
        ELSE
           setuid = NULL_DIMID
        END IF
        CALL IMPORT_NCVAR(grid%ps, ustep=ustep, varname=name, file=file &
             , setuid=setuid)
     END IF ! PS IS VALUE OR NAME
  END IF

  IF (TRIM(grid%p0%name) /= '') THEN
     ! CHECK FOR REAL VALUE
     READ(grid%p0%name,*,IOSTAT=iostat) p
     IF (iostat == 0) THEN    ! VALUE of P0
        CALL INIT_NCVAR(grid%p0)
        grid%p0%name  = 'p0'
        grid%p0%xtype = NF90_FLOAT
        !
        ! NO DIMENSIONS !!!
        grid%p0%ndims = 1
        ALLOCATE(grid%p0%dim(1), STAT=status)
        CALL ERRMSG(substr,status,4)
        grid%p0%dim(1)%name = TRIM(grid%p0%name)//'_dim'
        grid%p0%dim(1)%id   = NULL_DIMID
        grid%p0%dim(1)%len  = 1
        grid%p0%dim(1)%fuid = .false.
        grid%p0%dim(1)%varid = NULL_VARID
        !
        ! ATTRIBUTES ...
        ! ... LONGNAME ATTRIBUTE
        CALL ADD_NCATT(grid%p0, 'long_name'            &
                        ,vs='reference pressure')
        ! ... UNITS ATTRIBUT
        unit = ''
        READ(grid%p0%name,*,IOSTAT=iostat) p,unit
        IF (TRIM(unit) /= '') THEN
           CALL ADD_NCATT(grid%p0, 'units' ,vs=TRIM(unit))
        END IF
        !
        ! DATA
        CALL INIT_NARRAY(grid%p0%dat, 1, (/ 1 /), VTYPE_REAL)
        grid%p0%dat%vr(1) = p
     ELSE  ! NAME of P0
        name = TRIM(grid%p0%name)
        IF (ASSOCIATED(grid%timem%dim)) THEN
           setuid = grid%timem%dim(1)%id
        ELSE
           setuid = NULL_DIMID
        END IF
        CALL IMPORT_NCVAR(grid%p0, varname=name, file=file, setuid=setuid)
     END IF
  END IF

END SUBROUTINE IMPORT_GEOHYBGRID
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE EXPORT_GEOHYBGRID(grid)

  ! NO CHECKING IN THIS ROUTINE !!!

  IMPLICIT NONE

  ! I/O
  TYPE (geohybgrid), INTENT(INOUT) :: grid

  ! LOCAL
  CHARACTER(LEN=GRD_MAXSTRLEN) :: file   ! netCDF filename
  INTEGER                      :: ustep  ! step along UNLIMITED DIM.

  ustep = grid%t
  file = TRIM(grid%file)

  ! LONGITUDE
  IF (TRIM(grid%lonm%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%lonm, file=file)
  END IF
  IF (TRIM(grid%loni%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%loni, file=file)
  END IF

  ! LATITUDE
  IF (TRIM(grid%latm%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%latm, file=file)
  END IF
  IF (TRIM(grid%lati%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%lati, file=file)
  END IF

  ! TIME
  IF (TRIM(grid%timem%name) /= '') THEN
     grid%timem%ustep = ustep
     CALL EXPORT_NCVAR(grid%timem, file=file)
  END IF
!
  ! IF 'COMPLETED', timei has also the UNLIMITED-ATTRIBUTE
  ! DO NOT EXPORT AT THE MOMENT
!  IF (TRIM(grid%timei%name) /= '') THEN
!     CALL EXPORT_NCVAR(grid%timei, file=file)
!  END IF

  ! HYA
  IF (TRIM(grid%hyam%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%hyam, file=file)
  END IF
  IF (TRIM(grid%hyai%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%hyai, file=file)
  END IF

  ! HYB
  IF (TRIM(grid%hybm%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%hybm, file=file)
  END IF
  IF (TRIM(grid%hybi%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%hybi, file=file)
  END IF

  ! SURFACE PRESSURE
  IF (TRIM(grid%ps%name) /= '') THEN
     grid%ps%ustep = ustep
     CALL EXPORT_NCVAR(grid%ps, file=file)
  END IF

  ! REFERENCE PRESSURE
  IF (TRIM(grid%p0%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%p0, file=file)
  END IF

END SUBROUTINE EXPORT_GEOHYBGRID
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE SWITCH_GEOHYBGRID(g, lx, ly, lz)

  IMPLICIT NONE

  ! I/O
  TYPE(geohybgrid), INTENT(INOUT) :: g
  LOGICAL         , INTENT(IN)    :: lx, ly, lz

  IF (.NOT.lx) THEN
     CALL INIT_NCVAR(g%loni)
     CALL INIT_NCVAR(g%lonm)
  END IF

  IF (.NOT.ly) THEN
     CALL INIT_NCVAR(g%lati)
     CALL INIT_NCVAR(g%latm)
  END IF

  IF (.NOT.lz) THEN
     CALL INIT_NCVAR(g%hyai)
     CALL INIT_NCVAR(g%hyam)
     CALL INIT_NCVAR(g%hybi)
     CALL INIT_NCVAR(g%hybm)
     CALL INIT_NCVAR(g%p0)
     CALL INIT_NCVAR(g%ps)
  END IF

END SUBROUTINE SWITCH_GEOHYBGRID
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE CHECK_GEOHYBGRID(grid, ranges)

  IMPLICIT NONE

  ! I/O
  TYPE (geohybgrid),        INTENT(INOUT) :: grid
  REAL(DP), DIMENSION(4,2), INTENT(IN)    :: ranges

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'CHECK_GEOHYBGRID'

  ! LONGITUDE
  IF (grid%lonm%ndims > 1) THEN
     CALL RGMSG(substr, RGMLE, &
          'LONGITUDE VARIABLE '''//TRIM(grid%lonm%name)//'''', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IN FILE '''//TRIM(grid%file)//'''', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IS ',grid%lonm%ndims,'-DIMENSIONAL !')
  END IF

  IF (grid%loni%ndims > 1) THEN
     CALL RGMSG(substr, RGMLE, &
          'LONGITUDE INTERFACES VARIABLE '''//TRIM(grid%loni%name)//'''',&
          .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IN FILE '''//TRIM(grid%file)//'''', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IS ',grid%loni%ndims,'-DIMENSIONAL !')
  END IF

  ! LATITUDE
  IF (grid%latm%ndims > 1) THEN
     CALL RGMSG(substr, RGMLE, &
          'LATITUDE VARIABLE '''//TRIM(grid%latm%name)//'''', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IN FILE '''//TRIM(grid%file)//'''', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IS ',grid%latm%ndims,'-DIMENSIONAL !')
  END IF

  IF (grid%lati%ndims > 1) THEN
     CALL RGMSG(substr, RGMLE, &
          'LATITUDE INTERFACES VARIABLE '''//TRIM(grid%lati%name)//'''',&
          .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IN FILE '''//TRIM(grid%file)//'''', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IS ',grid%lati%ndims,'-DIMENSIONAL !')
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
!        IF ((ranges(1,1) == RGEMPTY).OR.(ranges(1,2) == RGEMPTY)) THEN
        IF (grid%loni%ndims /= 1) THEN
           IF ((ABS(ranges(1,1) - RGEMPTY) <= TINY(RGEMPTY)).OR. &
                (ABS(ranges(1,2) - RGEMPTY) <= TINY(RGEMPTY))) THEN
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
!        IF ((ranges(2,1) == RGEMPTY).OR.(ranges(2,2) == RGEMPTY)) THEN
        IF (grid%lati%ndims /= 1) THEN
           IF ((ABS(ranges(2,1) - RGEMPTY) <= TINY(RGEMPTY)).OR. &
                (ABS(ranges(2,2) - RGEMPTY) <= TINY(RGEMPTY))) THEN
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
!        IF ((ranges(3,1) == RGEMPTY).OR.(ranges(3,2) == RGEMPTY)) THEN
        IF (grid%hyai%ndims /= 1) THEN
           IF ((ABS(ranges(3,1) - RGEMPTY) <= TINY(RGEMPTY)).OR. &
                (ABS(ranges(3,2) - RGEMPTY) <= TINY(RGEMPTY))) THEN
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
!        IF ((ranges(4,1) == RGEMPTY).OR.(ranges(4,2) == RGEMPTY)) THEN
        IF (grid%hybi%ndims /= 1) THEN
           IF ((ABS(ranges(4,1) - RGEMPTY) <= TINY(RGEMPTY)).OR. &
                (ABS(ranges(4,2) - RGEMPTY) <= TINY(RGEMPTY))) THEN
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
  TYPE (geohybgrid), INTENT(IN)    :: gi  ! INPUT GRID
  TYPE (geohybgrid), INTENT(OUT)   :: go  ! OUTPUT GRID
  TYPE (geohybgrid), INTENT(INOUT) :: gx  ! INDEX 'GRID'
  LOGICAL, OPTIONAL, INTENT(IN)    :: reverse

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'SORT_GEOHYBGRID'
  LOGICAL :: lrev
  INTEGER :: i
  INTEGER :: olat, olon  ! order of dimensions in PS
  INTEGER, DIMENSION(:), ALLOCATABLE :: dil  ! dimension length vector of PS
  INTEGER, DIMENSION(:), POINTER     :: vec  ! element vector of PS
  INTEGER, DIMENSION(:), ALLOCATABLE :: svec ! element vector of sorted PS
  INTEGER                            :: vtype
  INTEGER                            :: status
  TYPE(narray)                       :: ts   ! temporal n-array for sigma level
  TYPE(narray)                       :: col  ! sigma level column
  TYPE(narray)                       :: cx   ! sigma level column sort indices

  INTEGER                             :: poshelp,dacchelp
  INTEGER                             :: ulimit
  INTEGER , DIMENSION(:), ALLOCATABLE :: daccelement
  INTEGER                             :: mhelp
  INTEGER, DIMENSION(:,:), POINTER    :: vechelp  ! element vector of PS
  INTEGER , DIMENSION(:,:), ALLOCATABLE :: daccelementhelp
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: svechelp ! element vect. of sorted PS

  NULLIFY(vechelp)

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

        SELECT CASE(SIZE(ts%dim))
        CASE(1) ! size(dim)=1
          DO i=1, ts%dim(1)
            vec(1) = i
            poshelp = vec(1)
            col%vd(i) = ts%vd(poshelp)
          END DO
        CASE(2) ! size(dim)=2
          DO i=1, ts%dim(1)
            vec(1) = i
            poshelp = vec(1)
            dacchelp = 1
            dacchelp = dacchelp*ts%dim(1)
            poshelp = poshelp + dacchelp*(vec(2)-1)
            col%vd(i) = ts%vd(poshelp)
          END DO
        CASE(3) ! size(dim)=3
          DO i=1, ts%dim(1)
            vec(1) = i
            poshelp = vec(1)
            dacchelp = 1
            dacchelp = dacchelp*ts%dim(1)
            poshelp = poshelp + dacchelp*(vec(2)-1)
            dacchelp = dacchelp*ts%dim(2)
            poshelp = poshelp + dacchelp*(vec(3)-1)
            col%vd(i) = ts%vd(poshelp)
          END DO
        CASE(4) ! size(dim)=4
          DO i=1, ts%dim(1)
            vec(1) = i
            poshelp = vec(1)
            dacchelp = 1
            dacchelp = dacchelp*ts%dim(1)
            poshelp = poshelp + dacchelp*(vec(2)-1)
            dacchelp = dacchelp*ts%dim(2)
            poshelp = poshelp + dacchelp*(vec(3)-1)
            dacchelp = dacchelp*ts%dim(3)
            poshelp = poshelp + dacchelp*(vec(4)-1)
            col%vd(i) = ts%vd(poshelp)
          END DO
        CASE DEFAULT ! size(dim).GT.4
          DO i=1, ts%dim(1)
            vec(1) = i
            col%vd(i) = ts%vd(POSITION(ts%dim,vec))
          END DO
        END SELECT

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

        SELECT CASE(SIZE(ts%dim))
        CASE(1) ! size(dim)=1
          DO i=1, ts%dim(1)
            vec(1) = i
            poshelp = vec(1)
            col%vd(i) = ts%vd(poshelp)
          END DO
        CASE(2) ! size(dim)=2
          DO i=1, ts%dim(1)
            vec(1) = i
            poshelp = vec(1)
            dacchelp = 1
            dacchelp = dacchelp*ts%dim(1)
            poshelp = poshelp + dacchelp*(vec(2)-1)
            col%vd(i) = ts%vd(poshelp)
          END DO
        CASE(3) ! size(dim)=3
          DO i=1, ts%dim(1)
            vec(1) = i
            poshelp = vec(1)
            dacchelp = 1
            dacchelp = dacchelp*ts%dim(1)
            poshelp = poshelp + dacchelp*(vec(2)-1)
            dacchelp = dacchelp*ts%dim(2)
            poshelp = poshelp + dacchelp*(vec(3)-1)
            col%vd(i) = ts%vd(poshelp)
          END DO
        CASE(4) ! size(dim)=4
          DO i=1, ts%dim(1)
            vec(1) = i
            poshelp = vec(1)
            dacchelp = 1
            dacchelp = dacchelp*ts%dim(1)
            poshelp = poshelp + dacchelp*(vec(2)-1)
            dacchelp = dacchelp*ts%dim(2)
            poshelp = poshelp + dacchelp*(vec(3)-1)
            dacchelp = dacchelp*ts%dim(3)
            poshelp = poshelp + dacchelp*(vec(4)-1)
            col%vd(i) = ts%vd(poshelp)
          END DO
        CASE DEFAULT ! size(dim).GT.4
          DO i=1, ts%dim(1)
            vec(1) = i
            col%vd(i) = ts%vd(POSITION(ts%dim,vec))
          END DO
        END SELECT

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
! bugfix 1
     svec(:) = 1
     dil(:) = 1

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
     ulimit = PRODUCT(dil)
     IF (ASSOCIATED(vec)) THEN
       IF (SIZE(vec) > 0) THEN
         DEALLOCATE(vec, STAT=status)
         CALL ERRMSG(substr,status,1)
       END IF
       NULLIFY(vec)
     END IF

     IF (.NOT.lrev) THEN

        SELECT CASE (vtype)
 
        CASE (VTYPE_REAL)

          SELECT CASE(SIZE(dil)) 
          CASE(1) !SIZE(dim)=1
            ALLOCATE(vechelp(ulimit,1), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              vechelp(i,1) = i
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              go%ps%dat%vr(i) = gi%ps%dat%vr(poshelp)
              gx%ps%dat%vi(i) = poshelp
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE(2) !SIZE(dim)=2
            ALLOCATE(daccelementhelp(ulimit,2))
            ALLOCATE(vechelp(ulimit,2), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              mhelp=i
              vechelp(i,1) = 0
              vechelp(i,2) = 0
              daccelementhelp(i,1) = 1
              daccelementhelp(i,2) = daccelementhelp(i,1)*dil(1)
              vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
              mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
              vechelp(i,1) = mhelp
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              dacchelp = 1
              dacchelp = dacchelp*dil(1)
              poshelp = poshelp + dacchelp*(svechelp(i,2)-1)
              go%ps%dat%vr(i) = gi%ps%dat%vr(poshelp)
              gx%ps%dat%vi(i) = poshelp
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(daccelementhelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE(3) !SIZE(dim)=3
            ALLOCATE(daccelementhelp(ulimit,3))
            ALLOCATE(vechelp(ulimit,3), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              mhelp=i
              vechelp(i,1) = 0
              vechelp(i,2) = 0
              vechelp(i,3) = 0
              daccelementhelp(i,1) = 1
              daccelementhelp(i,2) = daccelementhelp(i,1)*dil(1)
              daccelementhelp(i,3) = daccelementhelp(i,2)*dil(2)
              vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
              mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
              vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
              mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
              vechelp(i,1) = mhelp
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              dacchelp = 1
              dacchelp = dacchelp*dil(1)
              poshelp = poshelp + dacchelp*(svechelp(i,2)-1)
              dacchelp = dacchelp*dil(2)
              poshelp = poshelp + dacchelp*(svechelp(i,3)-1)
              go%ps%dat%vr(i) = gi%ps%dat%vr(poshelp)
              gx%ps%dat%vi(i) = poshelp
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(daccelementhelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE(4) !SIZE(dim)=4
            ALLOCATE(daccelementhelp(ulimit,4))
            ALLOCATE(vechelp(ulimit,4), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              mhelp=i
              vechelp(i,1) = 0
              vechelp(i,2) = 0
              vechelp(i,3) = 0
              vechelp(i,4) = 0
              daccelementhelp(i,1) = 1
              daccelementhelp(i,2) = daccelementhelp(i,1)*dil(1)
              daccelementhelp(i,3) = daccelementhelp(i,2)*dil(2)
              daccelementhelp(i,4) = daccelementhelp(i,3)*dil(3)
              vechelp(i,4) = (mhelp-1)/daccelementhelp(i,4)+1
              mhelp = mhelp - (vechelp(i,4)-1)*daccelementhelp(i,4)
              vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
              mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
              vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
              mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
              vechelp(i,1) = mhelp
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              dacchelp = 1
              dacchelp = dacchelp*dil(1)
              poshelp = poshelp + dacchelp*(svechelp(i,2)-1)
              dacchelp = dacchelp*dil(2)
              poshelp = poshelp + dacchelp*(svechelp(i,3)-1)
              dacchelp = dacchelp*dil(3)
              poshelp = poshelp + dacchelp*(svechelp(i,4)-1)
              go%ps%dat%vr(i) = gi%ps%dat%vr(poshelp)
              gx%ps%dat%vi(i) = poshelp
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(daccelementhelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE DEFAULT
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              CALL ELEMENT(dil,i,vec)   ! element vector in source array
              IF (olon > 0) THEN
                svec(olon) = gx%lonm%dat%vi(vec(olon))
              END IF
              IF (olat > 0) THEN
                svec(olat) = gx%latm%dat%vi(vec(olat))
              END IF
              !
              go%ps%dat%vr(i) = gi%ps%dat%vr(POSITION(dil, svec))
              gx%ps%dat%vi(i) = POSITION(dil, svec)
              !
              DEALLOCATE(vec, STAT=status)
              NULLIFY(vec)
            END DO   ! LOOP OVER ALL FIELD VALUES
          END SELECT

        CASE (VTYPE_DOUBLE)

          SELECT CASE(SIZE(dil)) 
          CASE(1) !SIZE(dim)=1
            ALLOCATE(vechelp(ulimit,1), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              vechelp(i,1) = i
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              go%ps%dat%vd(i) = gi%ps%dat%vd(poshelp)
              gx%ps%dat%vi(i) = poshelp
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE(2) !SIZE(dim)=2
            ALLOCATE(daccelementhelp(ulimit,2))
            ALLOCATE(vechelp(ulimit,2), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              mhelp=i
              vechelp(i,1) = 0
              vechelp(i,2) = 0
              daccelementhelp(i,1) = 1
              daccelementhelp(i,2) = daccelementhelp(i,1)*dil(1)
              vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
              mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
              vechelp(i,1) = mhelp
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              dacchelp = 1
              dacchelp = dacchelp*dil(1)
              poshelp = poshelp + dacchelp*(svechelp(i,2)-1)
              go%ps%dat%vd(i) = gi%ps%dat%vd(poshelp)
              gx%ps%dat%vi(i) = poshelp
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(daccelementhelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE(3) !SIZE(dim)=3
            ALLOCATE(daccelementhelp(ulimit,3))
            ALLOCATE(vechelp(ulimit,3), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              mhelp=i
              vechelp(i,1) = 0
              vechelp(i,2) = 0
              vechelp(i,3) = 0
              daccelementhelp(i,1) = 1
              daccelementhelp(i,2) = daccelementhelp(i,1)*dil(1)
              daccelementhelp(i,3) = daccelementhelp(i,2)*dil(2)
              vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
              mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
              vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
              mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
              vechelp(i,1) = mhelp
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              dacchelp = 1
              dacchelp = dacchelp*dil(1)
              poshelp = poshelp + dacchelp*(svechelp(i,2)-1)
              go%ps%dat%vd(i) = gi%ps%dat%vd(poshelp)
              gx%ps%dat%vi(i) = poshelp
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(daccelementhelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE(4) !SIZE(dim)=4
            ALLOCATE(daccelementhelp(ulimit,4))
            ALLOCATE(vechelp(ulimit,4), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              mhelp=i
              vechelp(i,1) = 0
              vechelp(i,2) = 0
              vechelp(i,3) = 0
              vechelp(i,4) = 0
              daccelementhelp(i,1) = 1
              daccelementhelp(i,2) = daccelementhelp(i,1)*dil(1)
              daccelementhelp(i,3) = daccelementhelp(i,2)*dil(2)
              daccelementhelp(i,4) = daccelementhelp(i,3)*dil(3)
              vechelp(i,4) = (mhelp-1)/daccelementhelp(i,4)+1
              mhelp = mhelp - (vechelp(i,4)-1)*daccelementhelp(i,4)
              vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
              mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
              vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
              mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
              vechelp(i,1) = mhelp
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              dacchelp = 1
              dacchelp = dacchelp*dil(1)
              poshelp = poshelp + dacchelp*(svechelp(i,2)-1)
              dacchelp = dacchelp*dil(2)
              poshelp = poshelp + dacchelp*(svechelp(i,3)-1)
              dacchelp = dacchelp*dil(3)
              poshelp = poshelp + dacchelp*(svechelp(i,4)-1)
              go%ps%dat%vd(i) = gi%ps%dat%vd(poshelp)
              gx%ps%dat%vi(i) = poshelp
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(daccelementhelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE DEFAULT
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              CALL ELEMENT(dil,i,vec)   ! element vector in source array
              IF (olon > 0) THEN
                svec(olon) = gx%lonm%dat%vi(vec(olon))
              END IF
              IF (olat > 0) THEN
                svec(olat) = gx%latm%dat%vi(vec(olat))
              END IF
              !
              go%ps%dat%vd(i) = gi%ps%dat%vd(POSITION(dil, svec))
              gx%ps%dat%vi(i) = POSITION(dil, svec)
              !
              DEALLOCATE(vec, STAT=status)
              NULLIFY(vec)
            END DO   ! LOOP OVER ALL FIELD VALUES
          END SELECT

        CASE (VTYPE_INT)

          SELECT CASE(SIZE(dil)) 
          CASE(1) !SIZE(dim)=1
            ALLOCATE(vechelp(ulimit,1), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              vechelp(i,1) = i
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              go%ps%dat%vi(i) = gi%ps%dat%vi(poshelp)
              gx%ps%dat%vi(i) = poshelp
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE(2) !SIZE(dim)=2
            ALLOCATE(daccelementhelp(ulimit,2))
            ALLOCATE(vechelp(ulimit,2), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              mhelp=i
              vechelp(i,1) = 0
              vechelp(i,2) = 0
              daccelementhelp(i,1) = 1
              daccelementhelp(i,2) = daccelementhelp(i,1)*dil(1)
              vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
              mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
              vechelp(i,1) = mhelp
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              dacchelp = 1
              dacchelp = dacchelp*dil(1)
              poshelp = poshelp + dacchelp*(svechelp(i,2)-1)
              go%ps%dat%vi(i) = gi%ps%dat%vi(poshelp)
              gx%ps%dat%vi(i) = poshelp
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(daccelementhelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE(3) !SIZE(dim)=3
            ALLOCATE(daccelementhelp(ulimit,3))
            ALLOCATE(vechelp(ulimit,3), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              mhelp=i
              vechelp(i,1) = 0
              vechelp(i,2) = 0
              vechelp(i,3) = 0
              daccelementhelp(i,1) = 1
              daccelementhelp(i,2) = daccelementhelp(i,1)*dil(1)
              daccelementhelp(i,3) = daccelementhelp(i,2)*dil(2)
              vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
              mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
              vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
              mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
              vechelp(i,1) = mhelp
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              dacchelp = 1
              dacchelp = dacchelp*dil(1)
              poshelp = poshelp + dacchelp*(svechelp(i,2)-1)
              dacchelp = dacchelp*dil(2)
              poshelp = poshelp + dacchelp*(svechelp(i,3)-1)
              go%ps%dat%vi(i) = gi%ps%dat%vi(poshelp)
              gx%ps%dat%vi(i) = poshelp
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(daccelementhelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE(4) !SIZE(dim)=4
            ALLOCATE(daccelementhelp(ulimit,4))
            ALLOCATE(vechelp(ulimit,4), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              mhelp=i
              vechelp(i,1) = 0
              vechelp(i,2) = 0
              vechelp(i,3) = 0
              vechelp(i,4) = 0
              daccelementhelp(i,1) = 1
              daccelementhelp(i,2) = daccelementhelp(i,1)*dil(1)
              daccelementhelp(i,3) = daccelementhelp(i,2)*dil(2)
              daccelementhelp(i,4) = daccelementhelp(i,3)*dil(3)
              vechelp(i,4) = (mhelp-1)/daccelementhelp(i,4)+1
              mhelp = mhelp - (vechelp(i,4)-1)*daccelementhelp(i,4)
              vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
              mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
              vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
              mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
              vechelp(i,1) = mhelp
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              dacchelp = 1
              dacchelp = dacchelp*dil(1)
              poshelp = poshelp + dacchelp*(svechelp(i,2)-1)
              dacchelp = dacchelp*dil(2)
              poshelp = poshelp + dacchelp*(svechelp(i,3)-1)
              dacchelp = dacchelp*dil(3)
              poshelp = poshelp + dacchelp*(svechelp(i,4)-1)
              go%ps%dat%vi(i) = gi%ps%dat%vi(poshelp)
              gx%ps%dat%vi(i) = poshelp
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(daccelementhelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE DEFAULT
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              CALL ELEMENT(dil,i,vec)   ! element vector in source array
              IF (olon > 0) THEN
                svec(olon) = gx%lonm%dat%vi(vec(olon))
              END IF
              IF (olat > 0) THEN
                svec(olat) = gx%latm%dat%vi(vec(olat))
              END IF
              !
              go%ps%dat%vi(i) = gi%ps%dat%vi(POSITION(dil, svec))
              gx%ps%dat%vi(i) = POSITION(dil, svec)
              !
              DEALLOCATE(vec, STAT=status)
              NULLIFY(vec)
            END DO   ! LOOP OVER ALL FIELD VALUES
          END SELECT

        CASE (VTYPE_BYTE)

          SELECT CASE(SIZE(dil)) 
          CASE(1) !SIZE(dim)=1
            ALLOCATE(vechelp(ulimit,1), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              vechelp(i,1) = i
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              go%ps%dat%vb(i) = gi%ps%dat%vb(poshelp)
              gx%ps%dat%vi(i) = poshelp
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE(2) !SIZE(dim)=2
            ALLOCATE(daccelementhelp(ulimit,2))
            ALLOCATE(vechelp(ulimit,2), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              mhelp=i
              vechelp(i,1) = 0
              vechelp(i,2) = 0
              daccelementhelp(i,1) = 1
              daccelementhelp(i,2) = daccelementhelp(i,1)*dil(1)
              vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
              mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
              vechelp(i,1) = mhelp
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              dacchelp = 1
              dacchelp = dacchelp*dil(1)
              poshelp = poshelp + dacchelp*(svechelp(i,2)-1)
              go%ps%dat%vb(i) = gi%ps%dat%vb(poshelp)
              gx%ps%dat%vi(i) = poshelp
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(daccelementhelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE(3) !SIZE(dim)=3
            ALLOCATE(daccelementhelp(ulimit,3))
            ALLOCATE(vechelp(ulimit,3), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              mhelp=i
              vechelp(i,1) = 0
              vechelp(i,2) = 0
              vechelp(i,3) = 0
              daccelementhelp(i,1) = 1
              daccelementhelp(i,2) = daccelementhelp(i,1)*dil(1)
              daccelementhelp(i,3) = daccelementhelp(i,2)*dil(2)
              vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
              mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
              vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
              mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
              vechelp(i,1) = mhelp
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              dacchelp = 1
              dacchelp = dacchelp*dil(1)
              poshelp = poshelp + dacchelp*(svechelp(i,2)-1)
              dacchelp = dacchelp*dil(2)
              poshelp = poshelp + dacchelp*(svechelp(i,3)-1)
              go%ps%dat%vb(i) = gi%ps%dat%vb(poshelp)
              gx%ps%dat%vi(i) = poshelp
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(daccelementhelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE(4) !SIZE(dim)=4
            ALLOCATE(daccelementhelp(ulimit,4))
            ALLOCATE(vechelp(ulimit,4), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              mhelp=i
              vechelp(i,1) = 0
              vechelp(i,2) = 0
              vechelp(i,3) = 0
              vechelp(i,4) = 0
              daccelementhelp(i,1) = 1
              daccelementhelp(i,2) = daccelementhelp(i,1)*dil(1)
              daccelementhelp(i,3) = daccelementhelp(i,2)*dil(2)
              daccelementhelp(i,4) = daccelementhelp(i,3)*dil(3)
              vechelp(i,4) = (mhelp-1)/daccelementhelp(i,4)+1
              mhelp = mhelp - (vechelp(i,4)-1)*daccelementhelp(i,4)
              vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
              mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
              vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
              mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
              vechelp(i,1) = mhelp
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              dacchelp = 1
              dacchelp = dacchelp*dil(1)
              poshelp = poshelp + dacchelp*(svechelp(i,2)-1)
              dacchelp = dacchelp*dil(2)
              poshelp = poshelp + dacchelp*(svechelp(i,3)-1)
              dacchelp = dacchelp*dil(3)
              poshelp = poshelp + dacchelp*(svechelp(i,4)-1)
              go%ps%dat%vb(i) = gi%ps%dat%vb(poshelp)
              gx%ps%dat%vi(i) = poshelp
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(daccelementhelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE DEFAULT
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              CALL ELEMENT(dil,i,vec)   ! element vector in source array
              IF (olon > 0) THEN
                svec(olon) = gx%lonm%dat%vi(vec(olon))
              END IF
              IF (olat > 0) THEN
                svec(olat) = gx%latm%dat%vi(vec(olat))
              END IF
              !
              go%ps%dat%vb(i) = gi%ps%dat%vb(POSITION(dil, svec))
              gx%ps%dat%vi(i) = POSITION(dil, svec)
              !
              DEALLOCATE(vec, STAT=status)
              NULLIFY(vec)
            END DO   ! LOOP OVER ALL FIELD VALUES
          END SELECT

        CASE (VTYPE_CHAR)

          SELECT CASE(SIZE(dil)) 
          CASE(1) !SIZE(dim)=1
            ALLOCATE(vechelp(ulimit,1), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              vechelp(i,1) = i
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              go%ps%dat%vc(i) = gi%ps%dat%vc(poshelp)
              gx%ps%dat%vi(i) = poshelp
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE(2) !SIZE(dim)=2
            ALLOCATE(daccelementhelp(ulimit,2))
            ALLOCATE(vechelp(ulimit,2), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              mhelp=i
              vechelp(i,1) = 0
              vechelp(i,2) = 0
              daccelementhelp(i,1) = 1
              daccelementhelp(i,2) = daccelementhelp(i,1)*dil(1)
              vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
              mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
              vechelp(i,1) = mhelp
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              dacchelp = 1
              dacchelp = dacchelp*dil(1)
              poshelp = poshelp + dacchelp*(svechelp(i,2)-1)
              go%ps%dat%vc(i) = gi%ps%dat%vc(poshelp)
              gx%ps%dat%vi(i) = poshelp
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(daccelementhelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE(3) !SIZE(dim)=3
            ALLOCATE(daccelementhelp(ulimit,3))
            ALLOCATE(vechelp(ulimit,3), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              mhelp=i
              vechelp(i,1) = 0
              vechelp(i,2) = 0
              vechelp(i,3) = 0
              daccelementhelp(i,1) = 1
              daccelementhelp(i,2) = daccelementhelp(i,1)*dil(1)
              daccelementhelp(i,3) = daccelementhelp(i,2)*dil(2)
              vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
              mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
              vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
              mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
              vechelp(i,1) = mhelp
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              dacchelp = 1
              dacchelp = dacchelp*dil(1)
              poshelp = poshelp + dacchelp*(svechelp(i,2)-1)
              dacchelp = dacchelp*dil(2)
              poshelp = poshelp + dacchelp*(svechelp(i,3)-1)
              go%ps%dat%vc(i) = gi%ps%dat%vc(poshelp)
              gx%ps%dat%vi(i) = poshelp
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(daccelementhelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE(4) !SIZE(dim)=4
            ALLOCATE(daccelementhelp(ulimit,4))
            ALLOCATE(vechelp(ulimit,4), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              mhelp=i
              vechelp(i,1) = 0
              vechelp(i,2) = 0
              vechelp(i,3) = 0
              vechelp(i,4) = 0
              daccelementhelp(i,1) = 1
              daccelementhelp(i,2) = daccelementhelp(i,1)*dil(1)
              daccelementhelp(i,3) = daccelementhelp(i,2)*dil(2)
              daccelementhelp(i,4) = daccelementhelp(i,3)*dil(3)
              vechelp(i,4) = (mhelp-1)/daccelementhelp(i,4)+1
              mhelp = mhelp - (vechelp(i,4)-1)*daccelementhelp(i,4)
              vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
              mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
              vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
              mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
              vechelp(i,1) = mhelp
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              dacchelp = 1
              dacchelp = dacchelp*dil(1)
              poshelp = poshelp + dacchelp*(svechelp(i,2)-1)
              dacchelp = dacchelp*dil(2)
              poshelp = poshelp + dacchelp*(svechelp(i,3)-1)
              dacchelp = dacchelp*dil(3)
              poshelp = poshelp + dacchelp*(svechelp(i,4)-1)
              go%ps%dat%vc(i) = gi%ps%dat%vc(poshelp)
              gx%ps%dat%vi(i) = poshelp
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(daccelementhelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE DEFAULT
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              CALL ELEMENT(dil,i,vec)   ! element vector in source array
              IF (olon > 0) THEN
                svec(olon) = gx%lonm%dat%vi(vec(olon))
              END IF
              IF (olat > 0) THEN
                svec(olat) = gx%latm%dat%vi(vec(olat))
              END IF
              !
              go%ps%dat%vc(i) = gi%ps%dat%vc(POSITION(dil, svec))
              gx%ps%dat%vi(i) = POSITION(dil, svec)
              !
              DEALLOCATE(vec, STAT=status)
              NULLIFY(vec)
            END DO   ! LOOP OVER ALL FIELD VALUES
          END SELECT

          CASE DEFAULT
             CALL RGMSG(substr, RGMLE, &
                   'ARRAY OF UNRECOGNIZED TYPE CANNOT BE SORTED !')
          END SELECT

     ELSE ! reverse

        SELECT CASE (vtype)
 
        CASE (VTYPE_REAL)

          SELECT CASE(SIZE(dil)) 
          CASE(1) !SIZE(dim)=1
            ALLOCATE(vechelp(ulimit,1), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              vechelp(i,1) = i
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              go%ps%dat%vr(poshelp) = gi%ps%dat%vr(i)
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE(2) !SIZE(dim)=2
            ALLOCATE(daccelementhelp(ulimit,2))
            ALLOCATE(vechelp(ulimit,2), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              mhelp=i
              vechelp(i,1) = 0
              vechelp(i,2) = 0
              daccelementhelp(i,1) = 1
              daccelementhelp(i,2) = daccelementhelp(i,1)*dil(1)
              vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
              mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
              vechelp(i,1) = mhelp
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              dacchelp = 1
              dacchelp = dacchelp*dil(1)
              poshelp = poshelp + dacchelp*(svechelp(i,2)-1)
              go%ps%dat%vr(poshelp) = gi%ps%dat%vr(i)
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(daccelementhelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE(3) !SIZE(dim)=3
            ALLOCATE(daccelementhelp(ulimit,3))
            ALLOCATE(vechelp(ulimit,3), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              mhelp=i
              vechelp(i,1) = 0
              vechelp(i,2) = 0
              vechelp(i,3) = 0
              daccelementhelp(i,1) = 1
              daccelementhelp(i,2) = daccelementhelp(i,1)*dil(1)
              daccelementhelp(i,3) = daccelementhelp(i,2)*dil(2)
              vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
              mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
              vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
              mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
              vechelp(i,1) = mhelp
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              dacchelp = 1
              dacchelp = dacchelp*dil(1)
              poshelp = poshelp + dacchelp*(svechelp(i,2)-1)
              dacchelp = dacchelp*dil(2)
              poshelp = poshelp + dacchelp*(svechelp(i,3)-1)
              go%ps%dat%vr(poshelp) = gi%ps%dat%vr(i)
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(daccelementhelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE(4) !SIZE(dim)=4
            ALLOCATE(daccelementhelp(ulimit,4))
            ALLOCATE(vechelp(ulimit,4), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              mhelp=i
              vechelp(i,1) = 0
              vechelp(i,2) = 0
              vechelp(i,3) = 0
              vechelp(i,4) = 0
              daccelementhelp(i,1) = 1
              daccelementhelp(i,2) = daccelementhelp(i,1)*dil(1)
              daccelementhelp(i,3) = daccelementhelp(i,2)*dil(2)
              daccelementhelp(i,4) = daccelementhelp(i,3)*dil(3)
              vechelp(i,4) = (mhelp-1)/daccelementhelp(i,4)+1
              mhelp = mhelp - (vechelp(i,4)-1)*daccelementhelp(i,4)
              vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
              mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
              vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
              mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
              vechelp(i,1) = mhelp
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              dacchelp = 1
              dacchelp = dacchelp*dil(1)
              poshelp = poshelp + dacchelp*(svechelp(i,2)-1)
              dacchelp = dacchelp*dil(2)
              poshelp = poshelp + dacchelp*(svechelp(i,3)-1)
              dacchelp = dacchelp*dil(3)
              poshelp = poshelp + dacchelp*(svechelp(i,4)-1)
              go%ps%dat%vr(poshelp) = gi%ps%dat%vr(i)
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(daccelementhelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE DEFAULT
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              CALL ELEMENT(dil,i,vec)   ! element vector in source array
              IF (olon > 0) THEN
                svec(olon) = gx%lonm%dat%vi(vec(olon))
              END IF
              IF (olat > 0) THEN
                svec(olat) = gx%latm%dat%vi(vec(olat))
              END IF
              !
              go%ps%dat%vr(POSITION(dil, svec)) = gi%ps%dat%vr(i)
              !
              DEALLOCATE(vec, STAT=status)
              NULLIFY(vec)
            END DO   ! LOOP OVER ALL FIELD VALUES
          END SELECT

        CASE (VTYPE_DOUBLE)

          SELECT CASE(SIZE(dil)) 
          CASE(1) !SIZE(dim)=1
            ALLOCATE(vechelp(ulimit,1), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              vechelp(i,1) = i
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              go%ps%dat%vd(poshelp) = gi%ps%dat%vd(i)
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE(2) !SIZE(dim)=2
            ALLOCATE(daccelementhelp(ulimit,2))
            ALLOCATE(vechelp(ulimit,2), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              mhelp=i
              vechelp(i,1) = 0
              vechelp(i,2) = 0
              daccelementhelp(i,1) = 1
              daccelementhelp(i,2) = daccelementhelp(i,1)*dil(1)
              vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
              mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
              vechelp(i,1) = mhelp
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              dacchelp = 1
              dacchelp = dacchelp*dil(1)
              poshelp = poshelp + dacchelp*(svechelp(i,2)-1)
              go%ps%dat%vd(poshelp) = gi%ps%dat%vd(i)
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(daccelementhelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE(3) !SIZE(dim)=3
            ALLOCATE(daccelementhelp(ulimit,3))
            ALLOCATE(vechelp(ulimit,3), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              mhelp=i
              vechelp(i,1) = 0
              vechelp(i,2) = 0
              vechelp(i,3) = 0
              daccelementhelp(i,1) = 1
              daccelementhelp(i,2) = daccelementhelp(i,1)*dil(1)
              daccelementhelp(i,3) = daccelementhelp(i,2)*dil(2)
              vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
              mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
              vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
              mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
              vechelp(i,1) = mhelp
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              dacchelp = 1
              dacchelp = dacchelp*dil(1)
              poshelp = poshelp + dacchelp*(svechelp(i,2)-1)
              go%ps%dat%vd(poshelp) = gi%ps%dat%vd(i)
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(daccelementhelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE(4) !SIZE(dim)=4
            ALLOCATE(daccelementhelp(ulimit,4))
            ALLOCATE(vechelp(ulimit,4), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              mhelp=i
              vechelp(i,1) = 0
              vechelp(i,2) = 0
              vechelp(i,3) = 0
              vechelp(i,4) = 0
              daccelementhelp(i,1) = 1
              daccelementhelp(i,2) = daccelementhelp(i,1)*dil(1)
              daccelementhelp(i,3) = daccelementhelp(i,2)*dil(2)
              daccelementhelp(i,4) = daccelementhelp(i,3)*dil(3)
              vechelp(i,4) = (mhelp-1)/daccelementhelp(i,4)+1
              mhelp = mhelp - (vechelp(i,4)-1)*daccelementhelp(i,4)
              vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
              mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
              vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
              mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
              vechelp(i,1) = mhelp
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              dacchelp = 1
              dacchelp = dacchelp*dil(1)
              poshelp = poshelp + dacchelp*(svechelp(i,2)-1)
              dacchelp = dacchelp*dil(2)
              poshelp = poshelp + dacchelp*(svechelp(i,3)-1)
              dacchelp = dacchelp*dil(3)
              poshelp = poshelp + dacchelp*(svechelp(i,4)-1)
              go%ps%dat%vd(poshelp) = gi%ps%dat%vd(i)
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(daccelementhelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE DEFAULT
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              CALL ELEMENT(dil,i,vec)   ! element vector in source array
              IF (olon > 0) THEN
                svec(olon) = gx%lonm%dat%vi(vec(olon))
              END IF
              IF (olat > 0) THEN
                svec(olat) = gx%latm%dat%vi(vec(olat))
              END IF
              !
              go%ps%dat%vd(POSITION(dil, svec)) = gi%ps%dat%vd(i)
              !
              DEALLOCATE(vec, STAT=status)
              NULLIFY(vec)
            END DO   ! LOOP OVER ALL FIELD VALUES
          END SELECT

        CASE (VTYPE_INT)

          SELECT CASE(SIZE(dil)) 
          CASE(1) !SIZE(dim)=1
            ALLOCATE(vechelp(ulimit,1), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              vechelp(i,1) = i
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              go%ps%dat%vi(poshelp) = gi%ps%dat%vi(i)
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE(2) !SIZE(dim)=2
            ALLOCATE(daccelementhelp(ulimit,2))
            ALLOCATE(vechelp(ulimit,2), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              mhelp=i
              vechelp(i,1) = 0
              vechelp(i,2) = 0
              daccelementhelp(i,1) = 1
              daccelementhelp(i,2) = daccelementhelp(i,1)*dil(1)
              vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
              mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
              vechelp(i,1) = mhelp
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              dacchelp = 1
              dacchelp = dacchelp*dil(1)
              poshelp = poshelp + dacchelp*(svechelp(i,2)-1)
              go%ps%dat%vi(poshelp) = gi%ps%dat%vi(i)
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(daccelementhelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE(3) !SIZE(dim)=3
            ALLOCATE(daccelementhelp(ulimit,3))
            ALLOCATE(vechelp(ulimit,3), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              mhelp=i
              vechelp(i,1) = 0
              vechelp(i,2) = 0
              vechelp(i,3) = 0
              daccelementhelp(i,1) = 1
              daccelementhelp(i,2) = daccelementhelp(i,1)*dil(1)
              daccelementhelp(i,3) = daccelementhelp(i,2)*dil(2)
              vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
              mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
              vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
              mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
              vechelp(i,1) = mhelp
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              dacchelp = 1
              dacchelp = dacchelp*dil(1)
              poshelp = poshelp + dacchelp*(svechelp(i,2)-1)
              dacchelp = dacchelp*dil(2)
              poshelp = poshelp + dacchelp*(svechelp(i,3)-1)
              go%ps%dat%vi(poshelp) = gi%ps%dat%vi(i)
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(daccelementhelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE(4) !SIZE(dim)=4
            ALLOCATE(daccelementhelp(ulimit,4))
            ALLOCATE(vechelp(ulimit,4), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              mhelp=i
              vechelp(i,1) = 0
              vechelp(i,2) = 0
              vechelp(i,3) = 0
              vechelp(i,4) = 0
              daccelementhelp(i,1) = 1
              daccelementhelp(i,2) = daccelementhelp(i,1)*dil(1)
              daccelementhelp(i,3) = daccelementhelp(i,2)*dil(2)
              daccelementhelp(i,4) = daccelementhelp(i,3)*dil(3)
              vechelp(i,4) = (mhelp-1)/daccelementhelp(i,4)+1
              mhelp = mhelp - (vechelp(i,4)-1)*daccelementhelp(i,4)
              vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
              mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
              vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
              mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
              vechelp(i,1) = mhelp
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              dacchelp = 1
              dacchelp = dacchelp*dil(1)
              poshelp = poshelp + dacchelp*(svechelp(i,2)-1)
              dacchelp = dacchelp*dil(2)
              poshelp = poshelp + dacchelp*(svechelp(i,3)-1)
              dacchelp = dacchelp*dil(3)
              poshelp = poshelp + dacchelp*(svechelp(i,4)-1)
              go%ps%dat%vi(poshelp) = gi%ps%dat%vi(i)
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(daccelementhelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE DEFAULT
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              CALL ELEMENT(dil,i,vec)   ! element vector in source array
              IF (olon > 0) THEN
                svec(olon) = gx%lonm%dat%vi(vec(olon))
              END IF
              IF (olat > 0) THEN
                svec(olat) = gx%latm%dat%vi(vec(olat))
              END IF
              !
              go%ps%dat%vi(POSITION(dil, svec)) = gi%ps%dat%vi(i)
              !
              DEALLOCATE(vec, STAT=status)
              NULLIFY(vec)
            END DO   ! LOOP OVER ALL FIELD VALUES
          END SELECT

        CASE (VTYPE_BYTE)

          SELECT CASE(SIZE(dil)) 
          CASE(1) !SIZE(dim)=1
            ALLOCATE(vechelp(ulimit,1), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              vechelp(i,1) = i
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              go%ps%dat%vb(poshelp) = gi%ps%dat%vb(i)
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE(2) !SIZE(dim)=2
            ALLOCATE(daccelementhelp(ulimit,2))
            ALLOCATE(vechelp(ulimit,2), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              mhelp=i
              vechelp(i,1) = 0
              vechelp(i,2) = 0
              daccelementhelp(i,1) = 1
              daccelementhelp(i,2) = daccelementhelp(i,1)*dil(1)
              vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
              mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
              vechelp(i,1) = mhelp
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              dacchelp = 1
              dacchelp = dacchelp*dil(1)
              poshelp = poshelp + dacchelp*(svechelp(i,2)-1)
              go%ps%dat%vb(poshelp) = gi%ps%dat%vb(i)
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(daccelementhelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE(3) !SIZE(dim)=3
            ALLOCATE(daccelementhelp(ulimit,3))
            ALLOCATE(vechelp(ulimit,3), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              mhelp=i
              vechelp(i,1) = 0
              vechelp(i,2) = 0
              vechelp(i,3) = 0
              daccelementhelp(i,1) = 1
              daccelementhelp(i,2) = daccelementhelp(i,1)*dil(1)
              daccelementhelp(i,3) = daccelementhelp(i,2)*dil(2)
              vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
              mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
              vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
              mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
              vechelp(i,1) = mhelp
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              dacchelp = 1
              dacchelp = dacchelp*dil(1)
              poshelp = poshelp + dacchelp*(svechelp(i,2)-1)
              dacchelp = dacchelp*dil(2)
              poshelp = poshelp + dacchelp*(svechelp(i,3)-1)
              go%ps%dat%vb(poshelp) = gi%ps%dat%vb(i)
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(daccelementhelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE(4) !SIZE(dim)=4
            ALLOCATE(daccelementhelp(ulimit,4))
            ALLOCATE(vechelp(ulimit,4), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              mhelp=i
              vechelp(i,1) = 0
              vechelp(i,2) = 0
              vechelp(i,3) = 0
              vechelp(i,4) = 0
              daccelementhelp(i,1) = 1
              daccelementhelp(i,2) = daccelementhelp(i,1)*dil(1)
              daccelementhelp(i,3) = daccelementhelp(i,2)*dil(2)
              daccelementhelp(i,4) = daccelementhelp(i,3)*dil(3)
              vechelp(i,4) = (mhelp-1)/daccelementhelp(i,4)+1
              mhelp = mhelp - (vechelp(i,4)-1)*daccelementhelp(i,4)
              vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
              mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
              vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
              mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
              vechelp(i,1) = mhelp
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              dacchelp = 1
              dacchelp = dacchelp*dil(1)
              poshelp = poshelp + dacchelp*(svechelp(i,2)-1)
              dacchelp = dacchelp*dil(2)
              poshelp = poshelp + dacchelp*(svechelp(i,3)-1)
              dacchelp = dacchelp*dil(3)
              poshelp = poshelp + dacchelp*(svechelp(i,4)-1)
              go%ps%dat%vb(poshelp) = gi%ps%dat%vb(i)
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(daccelementhelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE DEFAULT
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              CALL ELEMENT(dil,i,vec)   ! element vector in source array
              IF (olon > 0) THEN
                svec(olon) = gx%lonm%dat%vi(vec(olon))
              END IF
              IF (olat > 0) THEN
                svec(olat) = gx%latm%dat%vi(vec(olat))
              END IF
              !
              go%ps%dat%vb(POSITION(dil, svec)) = gi%ps%dat%vb(i)
              !
              DEALLOCATE(vec, STAT=status)
              NULLIFY(vec)
            END DO   ! LOOP OVER ALL FIELD VALUES
          END SELECT

        CASE (VTYPE_CHAR)

          SELECT CASE(SIZE(dil)) 
          CASE(1) !SIZE(dim)=1
            ALLOCATE(vechelp(ulimit,1), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              vechelp(i,1) = i
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              go%ps%dat%vc(poshelp) = gi%ps%dat%vc(i)
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE(2) !SIZE(dim)=2
            ALLOCATE(daccelementhelp(ulimit,2))
            ALLOCATE(vechelp(ulimit,2), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              mhelp=i
              vechelp(i,1) = 0
              vechelp(i,2) = 0
              daccelementhelp(i,1) = 1
              daccelementhelp(i,2) = daccelementhelp(i,1)*dil(1)
              vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
              mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
              vechelp(i,1) = mhelp
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              dacchelp = 1
              dacchelp = dacchelp*dil(1)
              poshelp = poshelp + dacchelp*(svechelp(i,2)-1)
              go%ps%dat%vc(poshelp) = gi%ps%dat%vc(i)
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(daccelementhelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE(3) !SIZE(dim)=3
            ALLOCATE(daccelementhelp(ulimit,3))
            ALLOCATE(vechelp(ulimit,3), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              mhelp=i
              vechelp(i,1) = 0
              vechelp(i,2) = 0
              vechelp(i,3) = 0
              daccelementhelp(i,1) = 1
              daccelementhelp(i,2) = daccelementhelp(i,1)*dil(1)
              daccelementhelp(i,3) = daccelementhelp(i,2)*dil(2)
              vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
              mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
              vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
              mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
              vechelp(i,1) = mhelp
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              dacchelp = 1
              dacchelp = dacchelp*dil(1)
              poshelp = poshelp + dacchelp*(svechelp(i,2)-1)
              dacchelp = dacchelp*dil(2)
              poshelp = poshelp + dacchelp*(svechelp(i,3)-1)
              go%ps%dat%vc(poshelp) = gi%ps%dat%vc(i)
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(daccelementhelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE(4) !SIZE(dim)=4
            ALLOCATE(daccelementhelp(ulimit,4))
            ALLOCATE(vechelp(ulimit,4), STAT=status)
            ALLOCATE(svechelp(ulimit,go%ps%ndims), STAT=status)
            svechelp(:,:) = 1
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              mhelp=i
              vechelp(i,1) = 0
              vechelp(i,2) = 0
              vechelp(i,3) = 0
              vechelp(i,4) = 0
              daccelementhelp(i,1) = 1
              daccelementhelp(i,2) = daccelementhelp(i,1)*dil(1)
              daccelementhelp(i,3) = daccelementhelp(i,2)*dil(2)
              daccelementhelp(i,4) = daccelementhelp(i,3)*dil(3)
              vechelp(i,4) = (mhelp-1)/daccelementhelp(i,4)+1
              mhelp = mhelp - (vechelp(i,4)-1)*daccelementhelp(i,4)
              vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
              mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
              vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
              mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
              vechelp(i,1) = mhelp
              IF (olon > 0) THEN
                svechelp(i,olon) = gx%lonm%dat%vi(vechelp(i,olon))
              END IF
              IF (olat > 0) THEN
                svechelp(i,olat) = gx%latm%dat%vi(vechelp(i,olat))
              END IF
              !
              poshelp = svechelp(i,1)
              dacchelp = 1
              dacchelp = dacchelp*dil(1)
              poshelp = poshelp + dacchelp*(svechelp(i,2)-1)
              dacchelp = dacchelp*dil(2)
              poshelp = poshelp + dacchelp*(svechelp(i,3)-1)
              dacchelp = dacchelp*dil(3)
              poshelp = poshelp + dacchelp*(svechelp(i,4)-1)
              go%ps%dat%vc(poshelp) = gi%ps%dat%vc(i)
              !
            END DO   ! LOOP OVER ALL FIELD VALUES
            DEALLOCATE(svechelp)
            DEALLOCATE(daccelementhelp)
            DEALLOCATE(vechelp, STAT=status)
            NULLIFY(vechelp)

          CASE DEFAULT
            DO i=1, ulimit         ! LOOP OVER ALL FIELD VALUES
              CALL ELEMENT(dil,i,vec)   ! element vector in source array
              IF (olon > 0) THEN
                svec(olon) = gx%lonm%dat%vi(vec(olon))
              END IF
              IF (olat > 0) THEN
                svec(olat) = gx%latm%dat%vi(vec(olat))
              END IF
              !
              go%ps%dat%vc(POSITION(dil, svec)) = gi%ps%dat%vc(i)
              !
              DEALLOCATE(vec, STAT=status)
              NULLIFY(vec)
            END DO   ! LOOP OVER ALL FIELD VALUES
          END SELECT

          CASE DEFAULT
             CALL RGMSG(substr, RGMLE, &
                   'ARRAY OF UNRECOGNIZED TYPE CANNOT BE SORTED !')
          END SELECT

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
    TYPE (ncvar), INTENT(INOUT) :: var

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'IDX_NCVAR'
    INTEGER :: i
    INTEGER :: status
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
SUBROUTINE GEOHYBGRID_AXES(g, a, g2, a2, pflag)

  ! Note: NO CHECKING
  !       g, g2 must be 'ordered', 'complete', and 'consistent'

  IMPLICIT NONE

  ! I/O
  TYPE (geohybgrid), INTENT(IN)      :: g  ! GEO-HYBRID-GRID
  TYPE (axis), DIMENSION(:), POINTER :: a  ! LIST OF AXES FOR REGRIDDING
  TYPE (geohybgrid), INTENT(IN)     , OPTIONAL :: g2  ! GEO-HYBRID-GRID
  TYPE (axis), DIMENSION(:), POINTER, OPTIONAL :: a2  ! LIST OF AXES
  LOGICAL, OPTIONAL, INTENT(IN)      :: pflag  ! .true.: pressure axis
                                               ! .false. sigma-axis (default)

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'GEOHYBGRID_AXES'
  INTEGER                        :: n      ! number of axes
  INTEGER                        :: status
  LOGICAL                        :: lp     ! local presure axis flag
  INTEGER                        :: i      ! counter
  INTEGER                        :: vtype
  INTEGER                        :: ndep_lat ! number of lat-dimension in list
  INTEGER                        :: ndep_lon ! number of lon-dimension in list
  INTEGER                        :: dpc      ! dependency counter
  TYPE (narray)                  :: zps      ! local surface pressure

  ! CHECK SUBROUTINE CALL
  IF (PRESENT(g2).NEQV.PRESENT(a2)) THEN
     CALL RGMSG(substr, RGMLE, &
          '2nd GRID AND 2nd AXES LIST MUST BE PRESENT ON SUBROUTINE CALL !')
  END IF

  ! INIT
  IF (PRESENT(pflag)) THEN
     lp = pflag
  ELSE
     lp = .false.     ! DEFAULT: sigma coordinates
  END IF
  ndep_lon = 0
  ndep_lat = 0

  ! COUNT AXES
  n = 1             ! START WITH MINIMUM 1 AXIS (LAST AXIS IS 'FREE' DIM)

  IF (QDEF_NCVAR(g%loni)) n=n+1
  IF (QDEF_NCVAR(g%lati)) n=n+1
  IF (QDEF_NCVAR(g%hyai).OR.QDEF_NCVAR(g%hybi)) n=n+1

  ! ALLOCATE SPACE
  ALLOCATE(a(n), STAT=status)
  CALL ERRMSG(substr,status,1)
  DO i=1, n
     CALL INIT_AXIS(a(i))
  END DO

  IF (PRESENT(g2)) THEN
     ALLOCATE(a2(n), STAT=status)
     CALL ERRMSG(substr,status,2)
     DO i=1, n
        CALL INIT_AXIS(a2(i))
     END DO
  END IF

  ! ASSIGN DATA
  n = 0

  CALL RGMSG(substr, RGMLI, 'AXIS CONSTRUCTION ...')

  ! NOTE: THE ORDER OF TESTING THE FOLLOWING DIMENSIONS
  !       (HERE: LON, LAT, LEV) SHOULD NOT BE CHANGED !
  !       IT HAS TO BE CONSISTENT WITH OTHER ROUTINES

  ! LONGITUDE
  IF (QDEF_NCVAR(g%loni)) THEN
     CALL RGMSG(substr, RGMLIC, ' ... ADDING LONGITUDE AXIS ...')
     n=n+1
     a(n)%lm     = .true.     ! LONGITUDE IS MODULO AXIS
     a(n)%ndp    = 1          ! LONGITUDE IS ...
     ALLOCATE(a(n)%dep(1), STAT=status)
     CALL ERRMSG(substr,status,3)
     a(n)%dep(1) = n          ! ... INDEPENDENT
     CALL COPY_NARRAY(a(n)%dat, g%loni%dat)
     ndep_lon    = n
     IF (PRESENT(g2)) THEN
        IF (QDEF_NCVAR(g2%loni)) THEN
           CALL RGMSG(substr, RGMLIC, ' ... DIFFERENT FOR OUTPUT GRID ...')
           a2(n)%lm  = .true.
           a2(n)%ndp = 1
           ALLOCATE(a2(n)%dep(1), STAT=status)
           CALL ERRMSG(substr,status,4)
           a2(n)%dep(1) = n
           CALL COPY_NARRAY(a2(n)%dat, g2%loni%dat)
        ELSE !UNDEFINED => INVARIANT
           CALL RGMSG(substr, RGMLIC, ' ... INVARIANT ...')
        END IF
     END IF
  END IF

  ! LATITUDE
  IF (QDEF_NCVAR(g%lati)) THEN
     CALL RGMSG(substr, RGMLIC, ' ... ADDING LATITUDE AXIS ...')
     n=n+1
     a(n)%lm     = .false.    ! LATITUDE IS NON-MODULO AXIS
     a(n)%ndp    = 1          ! LATITUDE IS ...
     ALLOCATE(a(n)%dep(1), STAT=status)
     CALL ERRMSG(substr,status,5)
     a(n)%dep(1) = n          ! ... INDEPENDENT
     CALL COPY_NARRAY(a(n)%dat, g%lati%dat)
     ndep_lat    = n
     IF (PRESENT(g2)) THEN
        IF (QDEF_NCVAR(g2%lati)) THEN
           CALL RGMSG(substr, RGMLIC, ' ... DIFFERENT FOR OUTPUT GRID ...')
           a2(n)%lm  = .false.
           a2(n)%ndp = 1
           ALLOCATE(a2(n)%dep(1), STAT=status)
           CALL ERRMSG(substr,status,6)
           a2(n)%dep(1) = n
           CALL COPY_NARRAY(a2(n)%dat, g2%lati%dat)
        ELSE !UNDEFINED => INVARIANT
           CALL RGMSG(substr, RGMLIC, ' ... INVARIANT ...')
        END IF
     END IF
     !
     ! TAKE INTO ACCOUNT SPHERICAL GEOMETRY ...
     vtype = a(n)%dat%type
     SELECT CASE(vtype)
     CASE(VTYPE_REAL)
        a(n)%dat%vr = COS(((a(n)%dat%vr - REAL(90., SP))/ &
             REAL(180., SP))*REAL(PI,SP))
     CASE(VTYPE_DOUBLE)
        a(n)%dat%vd = COS(((a(n)%dat%vd - REAL(90., DP))/ &
             REAL(180., DP))*PI)
     CASE(VTYPE_INT)
        CALL RGMSG(substr, RGMLE, &
             'LATITUDE AXIS OF TYPE INTEGER IS NOT SUPPORTED !')
     CASE(VTYPE_BYTE)
        CALL RGMSG(substr, RGMLE, &
             'LATITUDE AXIS OF TYPE BYTE IS NOT SUPPORTED !')
     CASE(VTYPE_CHAR)
        CALL RGMSG(substr, RGMLE, &
             'LATITUDE AXIS OF TYPE CHAR IS NOT SUPPORTED !')
     CASE(VTYPE_UNDEF)
        CALL RGMSG(substr, RGMLE, &
             'LATITUDE AXIS IS UNDEFINED !')
     CASE DEFAULT
        CALL RGMSG(substr, RGMLE, &
             'TYPE OF LATITUDE AXIS IS NOT RECOGNIZED !')
     END SELECT
     !
     ! ... ALSO FOR 2nd GRID
     IF (PRESENT(g2)) THEN
        IF (QDEF_NCVAR(g2%lati)) THEN
           vtype = a2(n)%dat%type
           SELECT CASE(vtype)
           CASE(VTYPE_REAL)
              a2(n)%dat%vr = COS(((a2(n)%dat%vr - REAL(90., SP))/ &
                   REAL(180., SP))*REAL(PI,SP))
           CASE(VTYPE_DOUBLE)
              a2(n)%dat%vd = COS(((a2(n)%dat%vd - REAL(90., DP))/ &
                   REAL(180., DP))*PI)
           CASE(VTYPE_INT)
              CALL RGMSG(substr, RGMLE, &
                   'LATITUDE AXIS OF TYPE INTEGER IS NOT SUPPORTED !')
           CASE(VTYPE_BYTE)
              CALL RGMSG(substr, RGMLE, &
                   'LATITUDE AXIS OF TYPE BYTE IS NOT SUPPORTED !')
           CASE(VTYPE_CHAR)
              CALL RGMSG(substr, RGMLE, &
                   'LATITUDE AXIS OF TYPE CHAR IS NOT SUPPORTED !')
           CASE(VTYPE_UNDEF)
              CALL RGMSG(substr, RGMLE, &
                   'LATITUDE AXIS IS UNDEFINED !')
           CASE DEFAULT
              CALL RGMSG(substr, RGMLE, &
                   'TYPE OF LATITUDE AXIS IS NOT RECOGNIZED !')
           END SELECT
        END IF
     END IF
  END IF

  ! LEVELS
  IF (QDEF_NCVAR(g%hyai).OR.QDEF_NCVAR(g%hybi)) THEN
     CALL RGMSG(substr, RGMLIC, ' ... ADDING VERTICAL AXIS ...')
     n=n+1
     a(n)%lm     = .false.     ! VERTICAL AXIS IS NON-MODULO AXIS

     ! CALCULATE AXIS DATA
     IF (lp) THEN
        CALL RGMSG(substr, RGMLIC, '     -> PRESSURE AXIS ...')
     ELSE
        CALL RGMSG(substr, RGMLIC, '     -> DIMENSIONLESS AXIS ...')
     END IF
     CALL PS2PS(g%ps, zps)
     CALL H2PSIG(a(n)%dat,g%hyai%dat,g%hybi%dat,zps,g%p0%dat,lp)
     CALL INIT_NARRAY(zps)

     ! SET DEPENDENCIES
     a(n)%ndp = a(n)%dat%n
     ALLOCATE(a(n)%dep(a(n)%ndp), STAT=status)
     CALL ERRMSG(substr,status,7)
     a(n)%dep(:) = 0
     dpc = 1
     a(n)%dep(dpc) = n
     !
     IF (a(n)%ndp > 1) THEN
        DO i=1, g%ps%ndims
           ! CHECK DIM LENGTH AND NAME/ID
           IF (QCMP_NCDIM(g%ps%dim(i), g%lonm%dim(1)) > 1) THEN ! LONGITUDE
              dpc = dpc + 1
              a(n)%dep(dpc) = ndep_lon
              CALL RGMSG(substr, RGMLIC, '     -> DEPENDING ON LONGITUDE ...')
           END IF
           IF (QCMP_NCDIM(g%ps%dim(i), g%latm%dim(1)) > 1) THEN ! LATITUDE
              dpc = dpc + 1
              a(n)%dep(dpc) = ndep_lat
              CALL RGMSG(substr, RGMLIC, '     -> DEPENDING ON LATITUDE ...')
           END IF
        END DO
        IF ((dpc > 1).AND.(dpc /= a(n)%ndp)) THEN
           CALL RGMSG(substr, RGMLE, 'DEPENDENCY MISMATCH !')
        END IF
     END IF

     ! 2nd GRID
     IF (PRESENT(g2)) THEN
        IF (QDEF_NCVAR(g2%hyai).OR.QDEF_NCVAR(g2%hybi)) THEN
           CALL RGMSG(substr, RGMLIC, ' ... DIFFERENT FOR OUTPUT GRID ...')
           a2(n)%lm     = .false.     ! VERTICAL AXIS IS NON-MODULO AXIS

           ! CALCULATE AXIS DATA
           CALL PS2PS(g2%ps, zps)
           CALL H2PSIG(a2(n)%dat,g2%hyai%dat,g2%hybi%dat, &
                zps,g2%p0%dat,lp)
           CALL INIT_NARRAY(zps)

           ! SET DEPENDENCIES
           a2(n)%ndp = a2(n)%dat%n
           ALLOCATE(a2(n)%dep(a2(n)%ndp), STAT=status)
           CALL ERRMSG(substr,status,8)
           a2(n)%dep(:) = 0
           dpc = 1
           a2(n)%dep(dpc) = n
           !
           IF (a2(n)%ndp > 1) THEN
              DO i=1, g2%ps%ndims
                 ! CHECK DIM LENGTH AND NAME/ID
                 IF (QDEF_NCVAR(g2%lonm)) THEN
                    IF (QCMP_NCDIM(g2%ps%dim(i), g2%lonm%dim(1)) > 1) THEN
                       ! LONGITUDE
                       dpc = dpc + 1
                       a2(n)%dep(dpc) = ndep_lon
                       CALL RGMSG(substr, RGMLIC, &
                            '     -> DEPENDING ON LONGITUDE ...')
                    END IF
                 END IF
                 IF (QDEF_NCVAR(g2%latm)) THEN
                    IF (QCMP_NCDIM(g2%ps%dim(i), g2%latm%dim(1)) > 1) THEN
                       ! LATITUDE
                       dpc = dpc + 1
                       a2(n)%dep(dpc) = ndep_lat
                       CALL RGMSG(substr, RGMLIC, &
                            '     -> DEPENDING ON LATITUDE ...')
                    END IF
                 END IF
              END DO
              IF ((dpc > 1).AND.(dpc /= a2(n)%ndp)) THEN
                 CALL RGMSG(substr, RGMLE, 'DEPENDENCY MISMATCH !')
              END IF
           END IF
        ELSE ! UNDEFINED => INVARIANT
           CALL RGMSG(substr, RGMLIC, ' ... INVARIANT ...')
        END IF
     END IF ! 2nd GRID

  END IF ! LEVELS

  CALL RGMSG(substr, RGMLIC, '... END AXES CONSTRUCTION !')

END SUBROUTINE GEOHYBGRID_AXES
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE PS2PS(var, na)

  IMPLICIT NONE

  ! I/O
  TYPE (ncvar) , INTENT(IN)  :: var
  TYPE (narray), INTENT(OUT) :: na

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'PS2PS'
  INTEGER :: i
  INTEGER :: uidpos
  INTEGER, DIMENSION(:), POINTER :: vec
  INTEGER :: dc
  INTEGER :: status
  INTEGER :: vtype

  NULLIFY(vec)

  uidpos = 0
  DO i=1, var%ndims
     IF (var%dim(i)%fuid) THEN
        uidpos = i
        EXIT
     END IF
  END DO

  IF (uidpos == 0) THEN  ! no unlimited ID
     CALL COPY_NARRAY(na, var%dat)
  ELSE                   ! 'remove' unlimited ID
     IF ((var%dim(uidpos)%len /= 1).OR.(var%dat%dim(uidpos) /= 1)) THEN
        CALL RGMSG(substr, RGMLE, &
             'DIMENSION LENGTH OF UNLIMITED DIMENSION MUST BE 1 !')
     END IF
     ALLOCATE(vec(var%ndims-1), STAT=status)
     CALL ERRMSG(substr,status,1)
     dc = 0
     DO i=1, var%ndims
        IF (.NOT.var%dim(i)%fuid) THEN
           dc = dc + 1
           vec(dc) = var%dim(i)%len
        END IF
     END DO
     vtype = var%dat%type
     CALL INIT_NARRAY(na, var%ndims-1, vec, vtype)
     SELECT CASE(vtype)
     CASE(VTYPE_REAL)
        na%vr(:) = var%dat%vr(:)
     CASE(VTYPE_DOUBLE)
        na%vd(:) = var%dat%vd(:)
     CASE(VTYPE_INT)
        na%vi(:) = var%dat%vi(:)
     CASE(VTYPE_BYTE)
        na%vb(:) = var%dat%vb(:)
     CASE(VTYPE_CHAR)
        na%vc(:) = var%dat%vc(:)
     CASE(VTYPE_UNDEF)
        CALL RGMSG(substr, RGMLE, 'N-ARRAY IS UNDEFINED !')
     CASE DEFAULT
        CALL RGMSG(substr, RGMLE, 'N-ARRAY IS UNRECOGNIZED !')
     END SELECT
     DEALLOCATE(vec,STAT=status)
     CALL ERRMSG(substr,status,2)
     NULLIFY(vec)
  END IF

END SUBROUTINE PS2PS
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE H2PSIG(psig, hya, hyb, ps, p0, lp)

  IMPLICIT NONE

  ! I/O
  TYPE (narray), INTENT(INOUT) :: psig
  TYPE (narray), INTENT(IN)    :: hya, hyb, ps, p0
  LOGICAL      , INTENT(IN)    :: lp   ! .true.  -> pressure axis
                                       ! .false. -> dimensionless axis
  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'H2PSIG'
  INTEGER :: dim   ! dimension length of linear array in psig
  INTEGER :: status
  INTEGER :: i
  INTEGER :: i1,i2 ! little helpers
  TYPE (narray) :: zhya, zhyb, zps, zp0     ! local arrays
  INTEGER, DIMENSION(:), POINTER :: vec     ! element vector

  INTEGER, DIMENSION(:,:), POINTER     :: vechelp  ! unpacked main pos. vector
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: daccelementhelp
  INTEGER, DIMENSION(:), ALLOCATABLE   :: i1help,i2help
  INTEGER                              :: mhelp,dacchelp,poshelp

  NULLIFY(vechelp)
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
     IF( (SIZE(psig%dim).GT.4) .OR. (SIZE(zhya%dim).GT.3) .OR. (SIZE(zps%dim).GT.3) ) THEN
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
     ELSE

       ALLOCATE(vechelp(dim,SIZE(psig%dim)), STAT=status)
       SELECT CASE(SIZE(psig%dim))
       CASE(1)
         DO i=1, dim
           vechelp(i,1) = i
         ENDDO
       CASE(2)
         ALLOCATE(daccelementhelp(dim,2))
         DO i=1, dim
           mhelp=i
           vechelp(i,1) = 0
           vechelp(i,2) = 0
           daccelementhelp(i,1) = 1
           daccelementhelp(i,2) = daccelementhelp(i,1)*psig%dim(1)
           vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
           mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
           vechelp(i,1) = mhelp
         ENDDO
       CASE(3)
         ALLOCATE(daccelementhelp(dim,3))
         DO i=1, dim
           mhelp=i
           vechelp(i,1) = 0
           vechelp(i,2) = 0
           vechelp(i,3) = 0
           daccelementhelp(i,1) = 1
           daccelementhelp(i,2) = daccelementhelp(i,1)*psig%dim(1)
           daccelementhelp(i,3) = daccelementhelp(i,2)*psig%dim(2)
           vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
           mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
           vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
           mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
           vechelp(i,1) = mhelp
         ENDDO
       CASE(4)
         ALLOCATE(daccelementhelp(dim,4))
         DO i=1, dim
           mhelp=i
           vechelp(i,1) = 0
           vechelp(i,2) = 0
           vechelp(i,3) = 0
           vechelp(i,4) = 0
           daccelementhelp(i,1) = 1
           daccelementhelp(i,2) = daccelementhelp(i,1)*psig%dim(1)
           daccelementhelp(i,3) = daccelementhelp(i,2)*psig%dim(2)
           daccelementhelp(i,4) = daccelementhelp(i,3)*psig%dim(3)
           vechelp(i,4) = (mhelp-1)/daccelementhelp(i,4)+1
           mhelp = mhelp - (vechelp(i,4)-1)*daccelementhelp(i,4)
           vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
           mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
           vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
           mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
           vechelp(i,1) = mhelp
         ENDDO
       END SELECT
       DEALLOCATE(daccelementhelp)

       ALLOCATE(i1help(dim))
       SELECT CASE(SIZE(zhya%dim))
       CASE(1)
         DO i=1, dim
           i1help(i) = vechelp(i,1)
         ENDDO
       CASE(2)
         DO i=1, dim
           poshelp = vechelp(i,1)
           dacchelp = 1
           dacchelp = dacchelp*zhya%dim(1)
           poshelp = poshelp + dacchelp*(vechelp(i,2)-1)
           i1help(i) = poshelp
         ENDDO
       CASE(3)
         DO i=1, dim
           poshelp = vechelp(i,1)
           dacchelp = 1
           dacchelp = dacchelp*zhya%dim(1)
           poshelp = poshelp + dacchelp*(vechelp(i,2)-1)
           dacchelp = dacchelp*zhya%dim(2)
           poshelp = poshelp + dacchelp*(vechelp(i,3)-1)
           i1help(i) = poshelp
         ENDDO
       END SELECT

       ALLOCATE(i2help(dim))
       SELECT CASE(SIZE(zps%dim))
       CASE(1)
         DO i=1, dim
           i2help(i) = vechelp(i,(zhya%n+1))
         ENDDO
       CASE(2)
         DO i=1, dim
           poshelp = vechelp(i,(zhya%n+1))
           dacchelp = 1
           dacchelp = dacchelp*zps%dim(1)
           poshelp = poshelp + dacchelp*(vechelp(i,(zhya%n+1)+1)-1)
           i2help(i) = poshelp
         ENDDO
       CASE(3)
         DO i=1, dim
           poshelp = vechelp(i,(zhya%n+1))
           dacchelp = 1
           dacchelp = dacchelp*zps%dim(1)
           poshelp = poshelp + dacchelp*(vechelp(i,(zhya%n+1)+1)-1)
           dacchelp = dacchelp*zps%dim(2)
           poshelp = poshelp + dacchelp*(vechelp(i,(zhya%n+1)+2)-1)
           i2help(i) = poshelp
         ENDDO
       END SELECT

       IF (lp) THEN  ! PRESSURE COORDINATES
         DO i=1, dim
           psig%vd(i) = (zhya%vd(i1help(i)) * &
                         zp0%vd(1) +   &
                         zhyb%vd(i1help(i)) * &
                         zps%vd(i2help(i)) )
         END DO
       ELSE          ! SIGMA COORDINATES
         DO i=1, dim
           psig%vd(i) = (zhya%vd(i1help(i)) * &
                         zp0%vd(1) +   &
                         zhyb%vd(i1help(i)) * &
                         zps%vd(i2help(i)) )/ &
                         zps%vd(i2help(i))
         END DO
       ENDIF

       DEALLOCATE(i1help)
       DEALLOCATE(i2help)
       DEALLOCATE(vechelp, STAT=status)
       NULLIFY(vechelp)


     ENDIF
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
SUBROUTINE COMPLETE_GEOHYBGRID(g, ranges, gx)

  ! Note: g must be already sorted

  IMPLICIT NONE

  ! I/O
  TYPE (geohybgrid),        INTENT(INOUT)           :: g
  REAL(DP), DIMENSION(4,2), INTENT(IN)              :: ranges
  TYPE (geohybgrid),        INTENT(INOUT), OPTIONAL :: gx   ! SORT INDICES

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'COMPLETE_GEOHYBGRID'

  CALL RGMSG(substr, RGMLI, 'FILE '''//TRIM(g%file)//'''')
  CALL RGMSG(substr, RGMLIC,'STEP ',g%t,' ')

  CALL RGMSG(substr, RGMLIC,'  loni <<--->> lonm ...')
  CALL IMMI_NARRAY(g%loni%dat, g%lonm%dat)
  CALL IMMI_NCVAR(g%loni, g%lonm)
  IF (PRESENT(gx)) THEN
     CALL IMMI_NARRAY_IDX(gx%loni%dat, gx%lonm%dat)
     CALL IMMI_NCVAR(gx%loni, gx%lonm)
  END IF
  ! CORRECT FOR RANGE
  CALL RNGADJ_NARRAY(g%loni%dat, ranges(1,:))

  CALL RGMSG(substr, RGMLIC,'  lati <<--->> latm ...')
  CALL IMMI_NARRAY(g%lati%dat, g%latm%dat)
  CALL IMMI_NCVAR(g%lati, g%latm)
  IF (PRESENT(gx)) THEN
     CALL IMMI_NARRAY_IDX(gx%lati%dat, gx%latm%dat)
     CALL IMMI_NCVAR(gx%lati, gx%latm)
  END IF
  ! CORRECT FOR RANGE
  CALL RNGADJ_NARRAY(g%lati%dat, ranges(2,:))

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

! --------------------------------------------
SUBROUTINE IMMI_NARRAY(nai, nam)

  IMPLICIT NONE

  ! I/O
  TYPE (narray), INTENT(INOUT) :: nai, nam

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'IMMI_NARRAY'
  INTEGER :: itype, mtype
  INTEGER :: i

  itype = nai%type
  mtype = nam%type

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
        CALL RGMSG(substr, RGMLE, &
             'INTERFACE ARRAY IS ',nai%n,'-DIMENSIONAL !')
     END IF
     IF (nam%n /= 1) THEN
        CALL RGMSG(substr, RGMLE, &
             'MID ARRAY IS ',nam%n,'-DIMENSIONAL !')
     END IF
     IF (nai%dim(1) /= (nam%dim(1)+1)) THEN
        CALL RGMSG(substr, RGMLE, &
             'SIZES OF INTERFACE ARRAY AND MID ARRAY ARE INCOMPATIBLE !')
     END IF
     CALL RGMSG(substr, RGMLIC, '  ... O.K.')
     RETURN
  END IF

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
SUBROUTINE IMMI_NARRAY_IDX(nai, nam)

  IMPLICIT NONE

  ! I/O
  TYPE (narray), INTENT(INOUT) :: nai, nam

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'IMMI_NARRAY_IDX'
  INTEGER :: itype, mtype
  INTEGER :: i
  INTEGER :: status

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
        i = nai%vi(nai%dim(1)-1)
        nai%vi(nai%dim(1)-1) = nai%vi(nai%dim(1))
        nai%vi(nai%dim(1)) = i
     END IF
     RETURN
  END IF

  CALL RGMSG(substr, RGMLE, &
       'N-ARRAYS MUST BE OF TYPE INTEGER OR UNDEFINED !')

END SUBROUTINE IMMI_NARRAY_IDX
! --------------------------------------------

! --------------------------------------------
  SUBROUTINE IMMI_NCVAR(vari, varm)

    IMPLICIT NONE

    ! I/O
    TYPE (ncvar), INTENT(INOUT) :: vari, varm

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'IMMI_NCVAR'
    INTEGER :: i
    INTEGER :: status
    CHARACTER(LEN=GRD_MAXSTRLEN) :: str
    INTEGER :: vtype

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
             !vari%dim(i)%fuid  = .false.  ! UNLIM-DIM IS ONLY ALONG MIDs !
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
  SUBROUTINE RNGADJ_NARRAY(na, mr, lmon)

    IMPLICIT NONE

    ! I/O
    TYPE (narray), INTENT(INOUT)        :: na
    REAL(DP),      INTENT(IN)           :: mr(2)
    LOGICAL,       INTENT(IN), OPTIONAL :: lmon

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'RNGADJ_NARRAY'
    INTEGER :: n
    INTEGER :: vtype
    LOGICAL :: llmon

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

END SUBROUTINE COMPLETE_GEOHYBGRID
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE CHECK_NCVAR_ON_GEOHYBGRID(var, g, dims, axes, ok)

  IMPLICIT NONE

  ! I/O
  TYPE (ncvar)     , INTENT(IN)  :: var    ! variable
  TYPE (geohybgrid), INTENT(IN)  :: g      ! grid
  INTEGER, DIMENSION(:), POINTER :: dims   ! order of g-dims in var
  INTEGER          , INTENT(OUT) :: axes(3)! dimension no. of lon -> lat -> lev
  LOGICAL          , INTENT(OUT) :: ok     ! conform ?

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'CHECK_NCVAR_ON_GEOHYBGRID'
  INTEGER :: i
  INTEGER :: ndimg            ! number of dimensions in g
  INTEGER :: ndimv            ! number of g-dims in var
  INTEGER :: status
  INTEGER :: plon, plat, plev ! position of lon, lat, lev in AXES-LIST !!!
  LOGICAL :: lvert            ! vertical axis present ?
  LOGICAL :: llat, llon       ! lat, lon axes present ?

  ! INIT
  ok = .true.
  axes(:) = 0
  plon = 0
  plat = 0
  plev = 0

  ! ALLOCATE SPACE FOR 'ORDER OF DIMENSIONS'
  ALLOCATE(dims(var%ndims), STAT=status)
  CALL ERRMSG(substr,status,1)
  dims(:) = 0

  ! GET NUMBER OF HYBRID-DIMENSIONS IN GRID
  ndimg = 0
  IF (QDEF_NCVAR(g%lonm)) THEN
     ndimg = ndimg + 1
     plon = ndimg
  END IF
  IF (QDEF_NCVAR(g%latm)) THEN
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
           IF (llon) THEN
              ndimv = ndimv + 1
              dims(i) = plon
              axes(GORD_LON) = i
              CYCLE
           END IF
        END IF
     END IF

     IF (.NOT.llat) THEN  ! CHECK LAT
        IF (ASSOCIATED(g%latm%dim)) THEN
           llat = (QCMP_NCDIM(g%latm%dim(1), var%dim(i)) > 1)
           IF (ASSOCIATED(g%lati%dim).AND.(.NOT.llat)) THEN
              llat = ((QCMP_NCDIM(g%latm%dim(1), var%dim(i)) > 0).AND. &
                      (TRIM(g%latm%dim(1)%name) ==                     &
                       TRIM(g%lati%dim(1)%name)//'_M') )
           END IF
           IF (llat) THEN
              ndimv = ndimv + 1
              dims(i) = plat
              axes(GORD_LAT) = i
              CYCLE
           END IF
        END IF
     END IF

     IF (.NOT.lvert) THEN  ! CHECK HYA
        IF (ASSOCIATED(g%hyam%dim)) THEN
           lvert = (QCMP_NCDIM(g%hyam%dim(1), var%dim(i)) > 1)
           IF (ASSOCIATED(g%hyai%dim).AND.(.NOT.lvert)) THEN
              lvert = ((QCMP_NCDIM(g%hyam%dim(1), var%dim(i)) > 0).AND. &
                       (TRIM(g%hyam%dim(1)%name) ==                     &
                        TRIM(g%hyai%dim(1)%name)//'_M') )
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
           IF (ASSOCIATED(g%hybi%dim).AND.(.NOT.lvert)) THEN
              lvert = ((QCMP_NCDIM(g%hybm%dim(1), var%dim(i)) > 0).AND. &
                       (TRIM(g%hybm%dim(1)%name) ==                     &
                        TRIM(g%hybi%dim(1)%name)//'_M') )
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

  ok = (ndimv >= ndimg)  ! ALL GRID DIMS MUST BE RECOGNIZED !!!
                         ! REMAINING DIMS ARE 'FREE'

END SUBROUTINE CHECK_NCVAR_ON_GEOHYBGRID
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE SORT_GEOHYBGRID_NCVAR(var, gx, axes, svar, reverse)

  IMPLICIT NONE

  ! I/O
  TYPE (ncvar)     , INTENT(IN)           :: var     ! input variable
  TYPE (geohybgrid), INTENT(IN)           :: gx      ! hybrid grid with  ...
                                                     ! ... index information
  INTEGER          , INTENT(IN)           :: axes(3) ! lon,lat,lev dim. no.
  TYPE (ncvar)     , INTENT(OUT)          :: svar    ! sorted variable
  LOGICAL          , INTENT(IN), OPTIONAL :: reverse

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'SORT_GEOHYBGRID_NCVAR'
  INTEGER                            :: i
  INTEGER                            :: vtype
  INTEGER, DIMENSION(:), ALLOCATABLE :: vdim   ! variable dimension vector
  INTEGER                            :: status
  INTEGER, DIMENSION(:), POINTER     :: vec    ! element vector
  LOGICAL                            :: lrev   ! local reverse flag

  LOGICAL                            ::  sw_hyam
  INTEGER                            ::  ulimit,poshelp,dacchelp,mhelp
  INTEGER, DIMENSION(:,:), POINTER   :: vechelp  ! element vector of PS
  INTEGER , DIMENSION(:,:), ALLOCATABLE :: daccelementhelp
 
  NULLIFY(vechelp) 
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

  sw_hyam = .FALSE.
  IF (axes(GORD_LEV) > 0) THEN
    IF (gx%hyam%dat%type /= VTYPE_UNDEF) THEN
      sw_hyam = .TRUE.
    ELSE
      sw_hyam = .FALSE.
    END IF
  END IF

  ulimit = PRODUCT(vdim)

  SELECT CASE(SIZE(vdim))

  CASE(1) ! size(vdim)=1
     ALLOCATE(vechelp(ulimit,1), STAT=status)
     
     SELECT CASE(vtype)
     CASE(VTYPE_INT)

       IF (lrev) THEN   ! UNSORT
         DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
            ! GET ELEMENT VECTOR
            vechelp(i,1) = i
            ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
            IF (axes(GORD_LON) > 0) THEN
               vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
            END IF
            IF (axes(GORD_LAT) > 0) THEN
               vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
            END IF
            IF (axes(GORD_LEV) > 0) THEN
               IF (sw_hyam) THEN
                  vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
               ELSE
                  vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
               END IF
            END IF
            ! COPY DATA ELEMENT
            poshelp = vechelp(i,1)
            svar%dat%vi(poshelp) = var%dat%vi(i)
         END DO

       ELSE             ! SORT

         DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
            ! GET ELEMENT VECTOR
            vechelp(i,1) = i
            ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
            IF (axes(GORD_LON) > 0) THEN
               vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
            END IF
            IF (axes(GORD_LAT) > 0) THEN
               vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
            END IF
            IF (axes(GORD_LEV) > 0) THEN
               IF (sw_hyam) THEN
                  vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
               ELSE
                  vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
               END IF
            END IF
            ! COPY DATA ELEMENT
            poshelp = vechelp(i,1)
            svar%dat%vi(i) = var%dat%vi(poshelp)
         END DO
       END IF

     CASE(VTYPE_REAL)

       IF (lrev) THEN   ! UNSORT
         DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
            ! GET ELEMENT VECTOR
            vechelp(i,1) = i
            ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
            IF (axes(GORD_LON) > 0) THEN
               vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
            END IF
            IF (axes(GORD_LAT) > 0) THEN
               vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
            END IF
            IF (axes(GORD_LEV) > 0) THEN
               IF (sw_hyam) THEN
                  vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
               ELSE
                  vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
               END IF
            END IF
            ! COPY DATA ELEMENT
            poshelp = vechelp(i,1)
            svar%dat%vr(poshelp) = var%dat%vr(i)
         END DO

       ELSE             ! SORT

         DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
            ! GET ELEMENT VECTOR
            vechelp(i,1) = i
            ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
            IF (axes(GORD_LON) > 0) THEN
               vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
            END IF
            IF (axes(GORD_LAT) > 0) THEN
               vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
            END IF
            IF (axes(GORD_LEV) > 0) THEN
               IF (sw_hyam) THEN
                  vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
               ELSE
                  vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
               END IF
            END IF
            ! COPY DATA ELEMENT
            poshelp = vechelp(i,1)
            svar%dat%vr(i) = var%dat%vr(poshelp)
         END DO
       END IF

     CASE(VTYPE_DOUBLE)

       IF (lrev) THEN   ! UNSORT
         DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
            ! GET ELEMENT VECTOR
            vechelp(i,1) = i
            ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
            IF (axes(GORD_LON) > 0) THEN
               vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
            END IF
            IF (axes(GORD_LAT) > 0) THEN
               vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
            END IF
            IF (axes(GORD_LEV) > 0) THEN
               IF (sw_hyam) THEN
                  vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
               ELSE
                  vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
               END IF
            END IF
            ! COPY DATA ELEMENT
            poshelp = vechelp(i,1)
            svar%dat%vd(poshelp) = var%dat%vd(i)
         END DO

       ELSE             ! SORT

         DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
            ! GET ELEMENT VECTOR
            vechelp(i,1) = i
            ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
            IF (axes(GORD_LON) > 0) THEN
               vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
            END IF
            IF (axes(GORD_LAT) > 0) THEN
               vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
            END IF
            IF (axes(GORD_LEV) > 0) THEN
               IF (sw_hyam) THEN
                  vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
               ELSE
                  vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
               END IF
            END IF
            ! COPY DATA ELEMENT
            poshelp = vechelp(i,1)
            svar%dat%vd(i) = var%dat%vd(poshelp)
         END DO
       END IF

     CASE(VTYPE_CHAR)

       IF (lrev) THEN   ! UNSORT
         DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
            ! GET ELEMENT VECTOR
            vechelp(i,1) = i
            ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
            IF (axes(GORD_LON) > 0) THEN
               vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
            END IF
            IF (axes(GORD_LAT) > 0) THEN
               vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
            END IF
            IF (axes(GORD_LEV) > 0) THEN
               IF (sw_hyam) THEN
                  vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
               ELSE
                  vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
               END IF
            END IF
            ! COPY DATA ELEMENT
            poshelp = vechelp(i,1)
            svar%dat%vc(poshelp) = var%dat%vc(i)
         END DO

       ELSE             ! SORT

         DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
            ! GET ELEMENT VECTOR
            vechelp(i,1) = i
            ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
            IF (axes(GORD_LON) > 0) THEN
               vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
            END IF
            IF (axes(GORD_LAT) > 0) THEN
               vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
            END IF
            IF (axes(GORD_LEV) > 0) THEN
               IF (sw_hyam) THEN
                  vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
               ELSE
                  vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
               END IF
            END IF
            ! COPY DATA ELEMENT
            poshelp = vechelp(i,1)
            svar%dat%vc(i) = var%dat%vc(poshelp)
         END DO
       END IF

     CASE(VTYPE_BYTE)

       IF (lrev) THEN   ! UNSORT
         DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
            ! GET ELEMENT VECTOR
            vechelp(i,1) = i
            ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
            IF (axes(GORD_LON) > 0) THEN
               vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
            END IF
            IF (axes(GORD_LAT) > 0) THEN
               vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
            END IF
            IF (axes(GORD_LEV) > 0) THEN
               IF (sw_hyam) THEN
                  vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
               ELSE
                  vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
               END IF
            END IF
            ! COPY DATA ELEMENT
            poshelp = vechelp(i,1)
            svar%dat%vb(poshelp) = var%dat%vb(i)
         END DO

       ELSE             ! SORT

         DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
            ! GET ELEMENT VECTOR
            vechelp(i,1) = i
            ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
            IF (axes(GORD_LON) > 0) THEN
               vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
            END IF
            IF (axes(GORD_LAT) > 0) THEN
               vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
            END IF
            IF (axes(GORD_LEV) > 0) THEN
               IF (sw_hyam) THEN
                  vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
               ELSE
                  vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
               END IF
            END IF
            ! COPY DATA ELEMENT
            poshelp = vechelp(i,1)
            svar%dat%vb(i) = var%dat%vb(poshelp)
         END DO
       END IF

     END SELECT
     DEALLOCATE(vechelp, STAT=status)
     NULLIFY(vechelp)

  CASE(2) ! size(vdim)=2
    ALLOCATE(daccelementhelp(ulimit,2))
    ALLOCATE(vechelp(ulimit,2), STAT=status)
    SELECT CASE(vtype)
    CASE(VTYPE_INT)

       IF (lrev) THEN   ! UNSORT
         DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
            ! GET ELEMENT VECTOR
            mhelp=i
            vechelp(i,1) = 0
            vechelp(i,2) = 0
            daccelementhelp(i,1) = 1
            daccelementhelp(i,2) = daccelementhelp(i,1)*vdim(1)
            vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
            mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
            vechelp(i,1) = mhelp
            ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
            IF (axes(GORD_LON) > 0) THEN
               vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
            END IF
            IF (axes(GORD_LAT) > 0) THEN
               vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
            END IF
            IF (axes(GORD_LEV) > 0) THEN
               IF (sw_hyam) THEN
                  vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
               ELSE
                  vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
               END IF
            END IF
            ! COPY DATA ELEMENT
            poshelp = vechelp(i,1)
            dacchelp = 1
            dacchelp = dacchelp*vdim(1)
            poshelp = poshelp + dacchelp*(vechelp(i,2)-1)
            svar%dat%vi(poshelp) = var%dat%vi(i)
         END DO

       ELSE             ! SORT

         DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
            ! GET ELEMENT VECTOR
            mhelp=i
            vechelp(i,1) = 0
            vechelp(i,2) = 0
            daccelementhelp(i,1) = 1
            daccelementhelp(i,2) = daccelementhelp(i,1)*vdim(1)
            vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
            mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
            vechelp(i,1) = mhelp
            ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
            IF (axes(GORD_LON) > 0) THEN
               vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
            END IF
            IF (axes(GORD_LAT) > 0) THEN
               vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
            END IF
            IF (axes(GORD_LEV) > 0) THEN
               IF (sw_hyam) THEN
                  vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
               ELSE
                  vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
               END IF
            END IF
            ! COPY DATA ELEMENT
            poshelp = vechelp(i,1)
            dacchelp = 1
            dacchelp = dacchelp*vdim(1)
            poshelp = poshelp + dacchelp*(vechelp(i,2)-1)
            svar%dat%vi(i) = var%dat%vi(poshelp)
         END DO
       END IF

    CASE(VTYPE_REAL)

       IF (lrev) THEN   ! UNSORT
         DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
            ! GET ELEMENT VECTOR
            mhelp=i
            vechelp(i,1) = 0
            vechelp(i,2) = 0
            daccelementhelp(i,1) = 1
            daccelementhelp(i,2) = daccelementhelp(i,1)*vdim(1)
            vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
            mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
            vechelp(i,1) = mhelp
            ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
            IF (axes(GORD_LON) > 0) THEN
               vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
            END IF
            IF (axes(GORD_LAT) > 0) THEN
               vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
            END IF
            IF (axes(GORD_LEV) > 0) THEN
               IF (sw_hyam) THEN
                  vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
               ELSE
                  vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
               END IF
            END IF
            ! COPY DATA ELEMENT
            poshelp = vechelp(i,1)
            dacchelp = 1
            dacchelp = dacchelp*vdim(1)
            poshelp = poshelp + dacchelp*(vechelp(i,2)-1)
            svar%dat%vr(poshelp) = var%dat%vr(i)
         END DO

       ELSE             ! SORT

         DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
            ! GET ELEMENT VECTOR
            mhelp=i
            vechelp(i,1) = 0
            vechelp(i,2) = 0
            daccelementhelp(i,1) = 1
            daccelementhelp(i,2) = daccelementhelp(i,1)*vdim(1)
            vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
            mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
            vechelp(i,1) = mhelp
            ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
            IF (axes(GORD_LON) > 0) THEN
               vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
            END IF
            IF (axes(GORD_LAT) > 0) THEN
               vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
            END IF
            IF (axes(GORD_LEV) > 0) THEN
               IF (sw_hyam) THEN
                  vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
               ELSE
                  vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
               END IF
            END IF
            ! COPY DATA ELEMENT
            poshelp = vechelp(i,1)
            dacchelp = 1
            dacchelp = dacchelp*vdim(1)
            poshelp = poshelp + dacchelp*(vechelp(i,2)-1)
            svar%dat%vr(i) = var%dat%vr(poshelp)
         END DO
       END IF

    CASE(VTYPE_DOUBLE)

       IF (lrev) THEN   ! UNSORT
         DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
            ! GET ELEMENT VECTOR
            mhelp=i
            vechelp(i,1) = 0
            vechelp(i,2) = 0
            daccelementhelp(i,1) = 1
            daccelementhelp(i,2) = daccelementhelp(i,1)*vdim(1)
            vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
            mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
            vechelp(i,1) = mhelp
            ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
            IF (axes(GORD_LON) > 0) THEN
               vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
            END IF
            IF (axes(GORD_LAT) > 0) THEN
               vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
            END IF
            IF (axes(GORD_LEV) > 0) THEN
               IF (sw_hyam) THEN
                  vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
               ELSE
                  vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
               END IF
            END IF
            ! COPY DATA ELEMENT
            poshelp = vechelp(i,1)
            dacchelp = 1
            dacchelp = dacchelp*vdim(1)
            poshelp = poshelp + dacchelp*(vechelp(i,2)-1)
            svar%dat%vd(poshelp) = var%dat%vd(i)
         END DO

       ELSE             ! SORT

         DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
            ! GET ELEMENT VECTOR
            mhelp=i
            vechelp(i,1) = 0
            vechelp(i,2) = 0
            daccelementhelp(i,1) = 1
            daccelementhelp(i,2) = daccelementhelp(i,1)*vdim(1)
            vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
            mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
            vechelp(i,1) = mhelp
            ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
            IF (axes(GORD_LON) > 0) THEN
               vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
            END IF
            IF (axes(GORD_LAT) > 0) THEN
               vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
            END IF
            IF (axes(GORD_LEV) > 0) THEN
               IF (sw_hyam) THEN
                  vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
               ELSE
                  vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
               END IF
            END IF
            ! COPY DATA ELEMENT
            poshelp = vechelp(i,1)
            dacchelp = 1
            dacchelp = dacchelp*vdim(1)
            poshelp = poshelp + dacchelp*(vechelp(i,2)-1)
            svar%dat%vd(i) = var%dat%vd(poshelp)
         END DO
       END IF

    CASE(VTYPE_CHAR)

       IF (lrev) THEN   ! UNSORT
         DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
            ! GET ELEMENT VECTOR
            mhelp=i
            vechelp(i,1) = 0
            vechelp(i,2) = 0
            daccelementhelp(i,1) = 1
            daccelementhelp(i,2) = daccelementhelp(i,1)*vdim(1)
            vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
            mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
            vechelp(i,1) = mhelp
            ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
            IF (axes(GORD_LON) > 0) THEN
               vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
            END IF
            IF (axes(GORD_LAT) > 0) THEN
               vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
            END IF
            IF (axes(GORD_LEV) > 0) THEN
               IF (sw_hyam) THEN
                  vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
               ELSE
                  vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
               END IF
            END IF
            ! COPY DATA ELEMENT
            poshelp = vechelp(i,1)
            dacchelp = 1
            dacchelp = dacchelp*vdim(1)
            poshelp = poshelp + dacchelp*(vechelp(i,2)-1)
            svar%dat%vc(poshelp) = var%dat%vc(i)
         END DO

       ELSE             ! SORT

         DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
            ! GET ELEMENT VECTOR
            mhelp=i
            vechelp(i,1) = 0
            vechelp(i,2) = 0
            daccelementhelp(i,1) = 1
            daccelementhelp(i,2) = daccelementhelp(i,1)*vdim(1)
            vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
            mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
            vechelp(i,1) = mhelp
            ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
            IF (axes(GORD_LON) > 0) THEN
               vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
            END IF
            IF (axes(GORD_LAT) > 0) THEN
               vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
            END IF
            IF (axes(GORD_LEV) > 0) THEN
               IF (sw_hyam) THEN
                  vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
               ELSE
                  vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
               END IF
            END IF
            ! COPY DATA ELEMENT
            poshelp = vechelp(i,1)
            dacchelp = 1
            dacchelp = dacchelp*vdim(1)
            poshelp = poshelp + dacchelp*(vechelp(i,2)-1)
            svar%dat%vc(i) = var%dat%vc(poshelp)
         END DO
       END IF

    CASE(VTYPE_BYTE)

       IF (lrev) THEN   ! UNSORT
         DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
            ! GET ELEMENT VECTOR
            mhelp=i
            vechelp(i,1) = 0
            vechelp(i,2) = 0
            daccelementhelp(i,1) = 1
            daccelementhelp(i,2) = daccelementhelp(i,1)*vdim(1)
            vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
            mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
            vechelp(i,1) = mhelp
            ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
            IF (axes(GORD_LON) > 0) THEN
               vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
            END IF
            IF (axes(GORD_LAT) > 0) THEN
               vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
            END IF
            IF (axes(GORD_LEV) > 0) THEN
               IF (sw_hyam) THEN
                  vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
               ELSE
                  vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
               END IF
            END IF
            ! COPY DATA ELEMENT
            poshelp = vechelp(i,1)
            dacchelp = 1
            dacchelp = dacchelp*vdim(1)
            poshelp = poshelp + dacchelp*(vechelp(i,2)-1)
            svar%dat%vb(poshelp) = var%dat%vb(i)
         END DO

       ELSE             ! SORT

         DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
            ! GET ELEMENT VECTOR
            mhelp=i
            vechelp(i,1) = 0
            vechelp(i,2) = 0
            daccelementhelp(i,1) = 1
            daccelementhelp(i,2) = daccelementhelp(i,1)*vdim(1)
            vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
            mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
            vechelp(i,1) = mhelp
            ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
            IF (axes(GORD_LON) > 0) THEN
               vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
            END IF
            IF (axes(GORD_LAT) > 0) THEN
               vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
            END IF
            IF (axes(GORD_LEV) > 0) THEN
               IF (sw_hyam) THEN
                  vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
               ELSE
                  vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
               END IF
            END IF
            ! COPY DATA ELEMENT
            poshelp = vechelp(i,1)
            dacchelp = 1
            dacchelp = dacchelp*vdim(1)
            poshelp = poshelp + dacchelp*(vechelp(i,2)-1)
            svar%dat%vb(i) = var%dat%vb(poshelp)
         END DO
       END IF

    END SELECT
    DEALLOCATE(daccelementhelp)
    DEALLOCATE(vechelp, STAT=status)
    NULLIFY(vechelp)

  CASE(3) ! size(vdim)=3
   ALLOCATE(daccelementhelp(ulimit,3))
   ALLOCATE(vechelp(ulimit,3), STAT=status)
   SELECT CASE(vtype)
   CASE(VTYPE_INT)

      IF (lrev) THEN   ! UNSORT
        DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
           ! GET ELEMENT VECTOR
           mhelp=i
           vechelp(i,1) = 0
           vechelp(i,2) = 0
           vechelp(i,3) = 0
           daccelementhelp(i,1) = 1
           daccelementhelp(i,2) = daccelementhelp(i,1)*vdim(1)
           daccelementhelp(i,3) = daccelementhelp(i,2)*vdim(2)
           vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
           mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
           vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
           mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
           vechelp(i,1) = mhelp

           ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
           IF (axes(GORD_LON) > 0) THEN
              vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
           END IF
           IF (axes(GORD_LAT) > 0) THEN
              vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
           END IF
           IF (axes(GORD_LEV) > 0) THEN
              IF (sw_hyam) THEN
                 vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
              ELSE
                 vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
              END IF
           END IF
           ! COPY DATA ELEMENT
           poshelp = vechelp(i,1)
           dacchelp = 1
           dacchelp = dacchelp*vdim(1)
           poshelp = poshelp + dacchelp*(vechelp(i,2)-1)
           dacchelp = dacchelp*vdim(2)
           poshelp = poshelp + dacchelp*(vechelp(i,3)-1)

           svar%dat%vi(poshelp) = var%dat%vi(i)
        END DO

      ELSE             ! SORT

        DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
           ! GET ELEMENT VECTOR
           mhelp=i
           vechelp(i,1) = 0
           vechelp(i,2) = 0
           vechelp(i,3) = 0
           daccelementhelp(i,1) = 1
           daccelementhelp(i,2) = daccelementhelp(i,1)*vdim(1)
           daccelementhelp(i,3) = daccelementhelp(i,2)*vdim(2)
           vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
           mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
           vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
           mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
           vechelp(i,1) = mhelp

           ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
           IF (axes(GORD_LON) > 0) THEN
              vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
           END IF
           IF (axes(GORD_LAT) > 0) THEN
              vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
           END IF
           IF (axes(GORD_LEV) > 0) THEN
              IF (sw_hyam) THEN
                 vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
              ELSE
                 vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
              END IF
           END IF
           ! COPY DATA ELEMENT
           poshelp = vechelp(i,1)
           dacchelp = 1
           dacchelp = dacchelp*vdim(1)
           poshelp = poshelp + dacchelp*(vechelp(i,2)-1)
           dacchelp = dacchelp*vdim(2)
           poshelp = poshelp + dacchelp*(vechelp(i,3)-1)

           svar%dat%vi(i) = var%dat%vi(poshelp)
        END DO
      END IF

   CASE(VTYPE_REAL)

      IF (lrev) THEN   ! UNSORT
        DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
           ! GET ELEMENT VECTOR
           mhelp=i
           vechelp(i,1) = 0
           vechelp(i,2) = 0
           vechelp(i,3) = 0
           daccelementhelp(i,1) = 1
           daccelementhelp(i,2) = daccelementhelp(i,1)*vdim(1)
           daccelementhelp(i,3) = daccelementhelp(i,2)*vdim(2)
           vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
           mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
           vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
           mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
           vechelp(i,1) = mhelp

!           ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
           IF (axes(GORD_LON) > 0) THEN
              vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
           END IF
           IF (axes(GORD_LAT) > 0) THEN
              vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
           END IF
           IF (axes(GORD_LEV) > 0) THEN
              IF (sw_hyam) THEN
                 vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
              ELSE
                 vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
              END IF
           END IF
           ! COPY DATA ELEMENT
           poshelp = vechelp(i,1)
           dacchelp = 1
           dacchelp = dacchelp*vdim(1)
           poshelp = poshelp + dacchelp*(vechelp(i,2)-1)
           dacchelp = dacchelp*vdim(2)
           poshelp = poshelp + dacchelp*(vechelp(i,3)-1)

           svar%dat%vr(poshelp) = var%dat%vr(i)
        END DO

      ELSE             ! SORT

        DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
           ! GET ELEMENT VECTOR
           mhelp=i
           vechelp(i,1) = 0
           vechelp(i,2) = 0
           vechelp(i,3) = 0
           daccelementhelp(i,1) = 1
           daccelementhelp(i,2) = daccelementhelp(i,1)*vdim(1)
           daccelementhelp(i,3) = daccelementhelp(i,2)*vdim(2)
           vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
           mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
           vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
           mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
           vechelp(i,1) = mhelp

           ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
           IF (axes(GORD_LON) > 0) THEN
              vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
           END IF
           IF (axes(GORD_LAT) > 0) THEN
              vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
           END IF
           IF (axes(GORD_LEV) > 0) THEN
              IF (sw_hyam) THEN
                 vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
              ELSE
                 vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
              END IF
           END IF
           ! COPY DATA ELEMENT
           poshelp = vechelp(i,1)
           dacchelp = 1
           dacchelp = dacchelp*vdim(1)
           poshelp = poshelp + dacchelp*(vechelp(i,2)-1)
           dacchelp = dacchelp*vdim(2)
           poshelp = poshelp + dacchelp*(vechelp(i,3)-1)

           svar%dat%vr(i) = var%dat%vr(poshelp)
        END DO
      END IF

   CASE(VTYPE_DOUBLE)

      IF (lrev) THEN   ! UNSORT
        DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
           ! GET ELEMENT VECTOR
           mhelp=i
           vechelp(i,1) = 0
           vechelp(i,2) = 0
           vechelp(i,3) = 0
           daccelementhelp(i,1) = 1
           daccelementhelp(i,2) = daccelementhelp(i,1)*vdim(1)
           daccelementhelp(i,3) = daccelementhelp(i,2)*vdim(2)
           vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
           mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
           vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
           mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
           vechelp(i,1) = mhelp

           ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
           IF (axes(GORD_LON) > 0) THEN
              vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
           END IF
           IF (axes(GORD_LAT) > 0) THEN
              vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
           END IF
           IF (axes(GORD_LEV) > 0) THEN
              IF (sw_hyam) THEN
                 vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
              ELSE
                 vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
              END IF
           END IF
           ! COPY DATA ELEMENT
           poshelp = vechelp(i,1)
           dacchelp = 1
           dacchelp = dacchelp*vdim(1)
           poshelp = poshelp + dacchelp*(vechelp(i,2)-1)
           dacchelp = dacchelp*vdim(2)
           poshelp = poshelp + dacchelp*(vechelp(i,3)-1)

           svar%dat%vd(poshelp) = var%dat%vd(i)
        END DO

      ELSE             ! SORT

        DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
           ! GET ELEMENT VECTOR
           mhelp=i
           vechelp(i,1) = 0
           vechelp(i,2) = 0
           vechelp(i,3) = 0
           daccelementhelp(i,1) = 1
           daccelementhelp(i,2) = daccelementhelp(i,1)*vdim(1)
           daccelementhelp(i,3) = daccelementhelp(i,2)*vdim(2)
           vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
           mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
           vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
           mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
           vechelp(i,1) = mhelp

           ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
           IF (axes(GORD_LON) > 0) THEN
              vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
           END IF
           IF (axes(GORD_LAT) > 0) THEN
              vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
           END IF
           IF (axes(GORD_LEV) > 0) THEN
              IF (sw_hyam) THEN
                 vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
              ELSE
                 vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
              END IF
           END IF
           ! COPY DATA ELEMENT
           poshelp = vechelp(i,1)
           dacchelp = 1
           dacchelp = dacchelp*vdim(1)
           poshelp = poshelp + dacchelp*(vechelp(i,2)-1)
           dacchelp = dacchelp*vdim(2)
           poshelp = poshelp + dacchelp*(vechelp(i,3)-1)

           svar%dat%vd(i) = var%dat%vd(poshelp)
        END DO
      END IF

   CASE(VTYPE_CHAR)

      IF (lrev) THEN   ! UNSORT
        DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
           ! GET ELEMENT VECTOR
           mhelp=i
           vechelp(i,1) = 0
           vechelp(i,2) = 0
           vechelp(i,3) = 0
           daccelementhelp(i,1) = 1
           daccelementhelp(i,2) = daccelementhelp(i,1)*vdim(1)
           daccelementhelp(i,3) = daccelementhelp(i,2)*vdim(2)
           vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
           mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
           vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
           mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
           vechelp(i,1) = mhelp

           ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
           IF (axes(GORD_LON) > 0) THEN
              vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
           END IF
           IF (axes(GORD_LAT) > 0) THEN
              vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
           END IF
           IF (axes(GORD_LEV) > 0) THEN
              IF (sw_hyam) THEN
                 vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
              ELSE
                 vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
              END IF
           END IF
           ! COPY DATA ELEMENT
           poshelp = vechelp(i,1)
           dacchelp = 1
           dacchelp = dacchelp*vdim(1)
           poshelp = poshelp + dacchelp*(vechelp(i,2)-1)
           dacchelp = dacchelp*vdim(2)
           poshelp = poshelp + dacchelp*(vechelp(i,3)-1)

           svar%dat%vc(poshelp) = var%dat%vc(i)
        END DO

      ELSE             ! SORT

        DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
           ! GET ELEMENT VECTOR
           mhelp=i
           vechelp(i,1) = 0
           vechelp(i,2) = 0
           vechelp(i,3) = 0
           daccelementhelp(i,1) = 1
           daccelementhelp(i,2) = daccelementhelp(i,1)*vdim(1)
           daccelementhelp(i,3) = daccelementhelp(i,2)*vdim(2)
           vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
           mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
           vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
           mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
           vechelp(i,1) = mhelp

           ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
           IF (axes(GORD_LON) > 0) THEN
              vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
           END IF
           IF (axes(GORD_LAT) > 0) THEN
              vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
           END IF
           IF (axes(GORD_LEV) > 0) THEN
              IF (sw_hyam) THEN
                 vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
              ELSE
                 vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
              END IF
           END IF
           ! COPY DATA ELEMENT
           poshelp = vechelp(i,1)
           dacchelp = 1
           dacchelp = dacchelp*vdim(1)
           poshelp = poshelp + dacchelp*(vechelp(i,2)-1)
           dacchelp = dacchelp*vdim(2)
           poshelp = poshelp + dacchelp*(vechelp(i,3)-1)

           svar%dat%vc(i) = var%dat%vc(poshelp)
        END DO
      END IF

   CASE(VTYPE_BYTE)

      IF (lrev) THEN   ! UNSORT
        DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
           ! GET ELEMENT VECTOR
           mhelp=i
           vechelp(i,1) = 0
           vechelp(i,2) = 0
           vechelp(i,3) = 0
           daccelementhelp(i,1) = 1
           daccelementhelp(i,2) = daccelementhelp(i,1)*vdim(1)
           daccelementhelp(i,3) = daccelementhelp(i,2)*vdim(2)
           vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
           mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
           vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
           mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
           vechelp(i,1) = mhelp

           ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
           IF (axes(GORD_LON) > 0) THEN
              vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
           END IF
           IF (axes(GORD_LAT) > 0) THEN
              vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
           END IF
           IF (axes(GORD_LEV) > 0) THEN
              IF (sw_hyam) THEN
                 vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
              ELSE
                 vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
              END IF
           END IF
           ! COPY DATA ELEMENT
           poshelp = vechelp(i,1)
           dacchelp = 1
           dacchelp = dacchelp*vdim(1)
           poshelp = poshelp + dacchelp*(vechelp(i,2)-1)
           dacchelp = dacchelp*vdim(2)
           poshelp = poshelp + dacchelp*(vechelp(i,3)-1)

           svar%dat%vb(poshelp) = var%dat%vb(i)
        END DO

      ELSE             ! SORT

        DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
           ! GET ELEMENT VECTOR
           mhelp=i
           vechelp(i,1) = 0
           vechelp(i,2) = 0
           vechelp(i,3) = 0
           daccelementhelp(i,1) = 1
           daccelementhelp(i,2) = daccelementhelp(i,1)*vdim(1)
           daccelementhelp(i,3) = daccelementhelp(i,2)*vdim(2)
           vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
           mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
           vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
           mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
           vechelp(i,1) = mhelp

           ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
           IF (axes(GORD_LON) > 0) THEN
              vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
           END IF
           IF (axes(GORD_LAT) > 0) THEN
              vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
           END IF
           IF (axes(GORD_LEV) > 0) THEN
              IF (sw_hyam) THEN
                 vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
              ELSE
                 vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
              END IF
           END IF
           ! COPY DATA ELEMENT
           poshelp = vechelp(i,1)
           dacchelp = 1
           dacchelp = dacchelp*vdim(1)
           poshelp = poshelp + dacchelp*(vechelp(i,2)-1)
           dacchelp = dacchelp*vdim(2)
           poshelp = poshelp + dacchelp*(vechelp(i,3)-1)

           svar%dat%vb(i) = var%dat%vb(poshelp)
        END DO
      END IF

   END SELECT
   DEALLOCATE(daccelementhelp)
   DEALLOCATE(vechelp, STAT=status)
   NULLIFY(vechelp)

  CASE(4) ! size(vdim)=4
    ALLOCATE(daccelementhelp(ulimit,4))
    ALLOCATE(vechelp(ulimit,4), STAT=status)
    SELECT CASE(vtype)
    CASE(VTYPE_INT)

       IF (lrev) THEN   ! UNSORT
         DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
            ! GET ELEMENT VECTOR
            mhelp=i
            vechelp(i,1) = 0
            vechelp(i,2) = 0
            vechelp(i,3) = 0
            vechelp(i,4) = 0
            daccelementhelp(i,1) = 1
            daccelementhelp(i,2) = daccelementhelp(i,1)*vdim(1)
            daccelementhelp(i,3) = daccelementhelp(i,2)*vdim(2)
            daccelementhelp(i,4) = daccelementhelp(i,3)*vdim(3)
            vechelp(i,4) = (mhelp-1)/daccelementhelp(i,4)+1
            mhelp = mhelp - (vechelp(i,4)-1)*daccelementhelp(i,4)
            vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
            mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
            vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
            mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
            vechelp(i,1) = mhelp

            ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
            IF (axes(GORD_LON) > 0) THEN
               vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
            END IF
            IF (axes(GORD_LAT) > 0) THEN
               vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
            END IF
            IF (axes(GORD_LEV) > 0) THEN
               IF (sw_hyam) THEN
                  vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
               ELSE
                  vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
               END IF
            END IF
            ! COPY DATA ELEMENT
            poshelp = vechelp(i,1)
            dacchelp = 1
            dacchelp = dacchelp*vdim(1)
            poshelp = poshelp + dacchelp*(vechelp(i,2)-1)
            dacchelp = dacchelp*vdim(2)
            poshelp = poshelp + dacchelp*(vechelp(i,3)-1)
            dacchelp = dacchelp*vdim(3)
            poshelp = poshelp + dacchelp*(vechelp(i,4)-1)

            svar%dat%vi(poshelp) = var%dat%vi(i)
         END DO

       ELSE             ! SORT

         DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
            ! GET ELEMENT VECTOR
            mhelp=i
            vechelp(i,1) = 0
            vechelp(i,2) = 0
            vechelp(i,3) = 0
            vechelp(i,4) = 0
            daccelementhelp(i,1) = 1
            daccelementhelp(i,2) = daccelementhelp(i,1)*vdim(1)
            daccelementhelp(i,3) = daccelementhelp(i,2)*vdim(2)
            daccelementhelp(i,4) = daccelementhelp(i,3)*vdim(3)
            vechelp(i,4) = (mhelp-1)/daccelementhelp(i,4)+1
            mhelp = mhelp - (vechelp(i,4)-1)*daccelementhelp(i,4)
            vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
            mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
            vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
            mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
            vechelp(i,1) = mhelp

            ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
            IF (axes(GORD_LON) > 0) THEN
               vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
            END IF
            IF (axes(GORD_LAT) > 0) THEN
               vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
            END IF
            IF (axes(GORD_LEV) > 0) THEN
               IF (sw_hyam) THEN
                  vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
               ELSE
                  vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
               END IF
            END IF
            ! COPY DATA ELEMENT
            poshelp = vechelp(i,1)
            dacchelp = 1
            dacchelp = dacchelp*vdim(1)
            poshelp = poshelp + dacchelp*(vechelp(i,2)-1)
            dacchelp = dacchelp*vdim(2)
            poshelp = poshelp + dacchelp*(vechelp(i,3)-1)
            dacchelp = dacchelp*vdim(3)
            poshelp = poshelp + dacchelp*(vechelp(i,4)-1)

            svar%dat%vi(i) = var%dat%vi(poshelp)
         END DO
       END IF

    CASE(VTYPE_REAL)

       IF (lrev) THEN   ! UNSORT
         DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
            ! GET ELEMENT VECTOR
            mhelp=i
            vechelp(i,1) = 0
            vechelp(i,2) = 0
            vechelp(i,3) = 0
            vechelp(i,4) = 0
            daccelementhelp(i,1) = 1
            daccelementhelp(i,2) = daccelementhelp(i,1)*vdim(1)
            daccelementhelp(i,3) = daccelementhelp(i,2)*vdim(2)
            daccelementhelp(i,4) = daccelementhelp(i,3)*vdim(3)
            vechelp(i,4) = (mhelp-1)/daccelementhelp(i,4)+1
            mhelp = mhelp - (vechelp(i,4)-1)*daccelementhelp(i,4)
            vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
            mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
            vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
            mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
            vechelp(i,1) = mhelp

            ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
            IF (axes(GORD_LON) > 0) THEN
               vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
            END IF
            IF (axes(GORD_LAT) > 0) THEN
               vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
            END IF
            IF (axes(GORD_LEV) > 0) THEN
               IF (sw_hyam) THEN
                  vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
               ELSE
                  vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
               END IF
            END IF
            ! COPY DATA ELEMENT
            poshelp = vechelp(i,1)
            dacchelp = 1
            dacchelp = dacchelp*vdim(1)
            poshelp = poshelp + dacchelp*(vechelp(i,2)-1)
            dacchelp = dacchelp*vdim(2)
            poshelp = poshelp + dacchelp*(vechelp(i,3)-1)
            dacchelp = dacchelp*vdim(3)
            poshelp = poshelp + dacchelp*(vechelp(i,4)-1)

            svar%dat%vr(poshelp) = var%dat%vr(i)
         END DO

       ELSE             ! SORT

         DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
            ! GET ELEMENT VECTOR
            mhelp=i
            vechelp(i,1) = 0
            vechelp(i,2) = 0
            vechelp(i,3) = 0
            vechelp(i,4) = 0
            daccelementhelp(i,1) = 1
            daccelementhelp(i,2) = daccelementhelp(i,1)*vdim(1)
            daccelementhelp(i,3) = daccelementhelp(i,2)*vdim(2)
            daccelementhelp(i,4) = daccelementhelp(i,3)*vdim(3)
            vechelp(i,4) = (mhelp-1)/daccelementhelp(i,4)+1
            mhelp = mhelp - (vechelp(i,4)-1)*daccelementhelp(i,4)
            vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
            mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
            vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
            mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
            vechelp(i,1) = mhelp

            ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
            IF (axes(GORD_LON) > 0) THEN
               vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
            END IF
            IF (axes(GORD_LAT) > 0) THEN
               vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
            END IF
            IF (axes(GORD_LEV) > 0) THEN
               IF (sw_hyam) THEN
                  vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
               ELSE
                  vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
               END IF
            END IF
            ! COPY DATA ELEMENT
            poshelp = vechelp(i,1)
            dacchelp = 1
            dacchelp = dacchelp*vdim(1)
            poshelp = poshelp + dacchelp*(vechelp(i,2)-1)
            dacchelp = dacchelp*vdim(2)
            poshelp = poshelp + dacchelp*(vechelp(i,3)-1)
            dacchelp = dacchelp*vdim(3)
            poshelp = poshelp + dacchelp*(vechelp(i,4)-1)

            svar%dat%vr(i) = var%dat%vr(poshelp)
         END DO
       END IF

    CASE(VTYPE_DOUBLE)

       IF (lrev) THEN   ! UNSORT
         DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
            ! GET ELEMENT VECTOR
            mhelp=i
            vechelp(i,1) = 0
            vechelp(i,2) = 0
            vechelp(i,3) = 0
            vechelp(i,4) = 0
            daccelementhelp(i,1) = 1
            daccelementhelp(i,2) = daccelementhelp(i,1)*vdim(1)
            daccelementhelp(i,3) = daccelementhelp(i,2)*vdim(2)
            daccelementhelp(i,4) = daccelementhelp(i,3)*vdim(3)
            vechelp(i,4) = (mhelp-1)/daccelementhelp(i,4)+1
            mhelp = mhelp - (vechelp(i,4)-1)*daccelementhelp(i,4)
            vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
            mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
            vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
            mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
            vechelp(i,1) = mhelp

            ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
            IF (axes(GORD_LON) > 0) THEN
               vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
            END IF
            IF (axes(GORD_LAT) > 0) THEN
               vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
            END IF
            IF (axes(GORD_LEV) > 0) THEN
               IF (sw_hyam) THEN
                  vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
               ELSE
                  vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
               END IF
            END IF
            ! COPY DATA ELEMENT
            poshelp = vechelp(i,1)
            dacchelp = 1
            dacchelp = dacchelp*vdim(1)
            poshelp = poshelp + dacchelp*(vechelp(i,2)-1)
            dacchelp = dacchelp*vdim(2)
            poshelp = poshelp + dacchelp*(vechelp(i,3)-1)
            dacchelp = dacchelp*vdim(3)
            poshelp = poshelp + dacchelp*(vechelp(i,4)-1)

            svar%dat%vd(poshelp) = var%dat%vd(i)
         END DO

       ELSE             ! SORT

         DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
            ! GET ELEMENT VECTOR
            mhelp=i
            vechelp(i,1) = 0
            vechelp(i,2) = 0
            vechelp(i,3) = 0
            vechelp(i,4) = 0
            daccelementhelp(i,1) = 1
            daccelementhelp(i,2) = daccelementhelp(i,1)*vdim(1)
            daccelementhelp(i,3) = daccelementhelp(i,2)*vdim(2)
            daccelementhelp(i,4) = daccelementhelp(i,3)*vdim(3)
            vechelp(i,4) = (mhelp-1)/daccelementhelp(i,4)+1
            mhelp = mhelp - (vechelp(i,4)-1)*daccelementhelp(i,4)
            vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
            mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
            vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
            mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
            vechelp(i,1) = mhelp

            ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
            IF (axes(GORD_LON) > 0) THEN
               vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
            END IF
            IF (axes(GORD_LAT) > 0) THEN
               vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
            END IF
            IF (axes(GORD_LEV) > 0) THEN
               IF (sw_hyam) THEN
                  vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
               ELSE
                  vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
               END IF
            END IF
            ! COPY DATA ELEMENT
            poshelp = vechelp(i,1)
            dacchelp = 1
            dacchelp = dacchelp*vdim(1)
            poshelp = poshelp + dacchelp*(vechelp(i,2)-1)
            dacchelp = dacchelp*vdim(2)
            poshelp = poshelp + dacchelp*(vechelp(i,3)-1)
            dacchelp = dacchelp*vdim(3)
            poshelp = poshelp + dacchelp*(vechelp(i,4)-1)

            svar%dat%vd(i) = var%dat%vd(poshelp)
         END DO
       END IF

    CASE(VTYPE_CHAR)

       IF (lrev) THEN   ! UNSORT
         DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
            ! GET ELEMENT VECTOR
            mhelp=i
            vechelp(i,1) = 0
            vechelp(i,2) = 0
            vechelp(i,3) = 0
            vechelp(i,4) = 0
            daccelementhelp(i,1) = 1
            daccelementhelp(i,2) = daccelementhelp(i,1)*vdim(1)
            daccelementhelp(i,3) = daccelementhelp(i,2)*vdim(2)
            daccelementhelp(i,4) = daccelementhelp(i,3)*vdim(3)
            vechelp(i,4) = (mhelp-1)/daccelementhelp(i,4)+1
            mhelp = mhelp - (vechelp(i,4)-1)*daccelementhelp(i,4)
            vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
            mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
            vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
            mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
            vechelp(i,1) = mhelp

            ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
            IF (axes(GORD_LON) > 0) THEN
               vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
            END IF
            IF (axes(GORD_LAT) > 0) THEN
               vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
            END IF
            IF (axes(GORD_LEV) > 0) THEN
               IF (sw_hyam) THEN
                  vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
               ELSE
                  vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
               END IF
            END IF
            ! COPY DATA ELEMENT
            poshelp = vechelp(i,1)
            dacchelp = 1
            dacchelp = dacchelp*vdim(1)
            poshelp = poshelp + dacchelp*(vechelp(i,2)-1)
            dacchelp = dacchelp*vdim(2)
            poshelp = poshelp + dacchelp*(vechelp(i,3)-1)
            dacchelp = dacchelp*vdim(3)
            poshelp = poshelp + dacchelp*(vechelp(i,4)-1)

            svar%dat%vc(poshelp) = var%dat%vc(i)
         END DO

       ELSE             ! SORT

         DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
            ! GET ELEMENT VECTOR
            mhelp=i
            vechelp(i,1) = 0
            vechelp(i,2) = 0
            vechelp(i,3) = 0
            vechelp(i,4) = 0
            daccelementhelp(i,1) = 1
            daccelementhelp(i,2) = daccelementhelp(i,1)*vdim(1)
            daccelementhelp(i,3) = daccelementhelp(i,2)*vdim(2)
            daccelementhelp(i,4) = daccelementhelp(i,3)*vdim(3)
            vechelp(i,4) = (mhelp-1)/daccelementhelp(i,4)+1
            mhelp = mhelp - (vechelp(i,4)-1)*daccelementhelp(i,4)
            vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
            mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
            vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
            mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
            vechelp(i,1) = mhelp

            ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
            IF (axes(GORD_LON) > 0) THEN
               vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
            END IF
            IF (axes(GORD_LAT) > 0) THEN
               vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
            END IF
            IF (axes(GORD_LEV) > 0) THEN
               IF (sw_hyam) THEN
                  vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
               ELSE
                  vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
               END IF
            END IF
            ! COPY DATA ELEMENT
            poshelp = vechelp(i,1)
            dacchelp = 1
            dacchelp = dacchelp*vdim(1)
            poshelp = poshelp + dacchelp*(vechelp(i,2)-1)
            dacchelp = dacchelp*vdim(2)
            poshelp = poshelp + dacchelp*(vechelp(i,3)-1)
            dacchelp = dacchelp*vdim(3)
            poshelp = poshelp + dacchelp*(vechelp(i,4)-1)

            svar%dat%vc(i) = var%dat%vc(poshelp)
         END DO
       END IF

    CASE(VTYPE_BYTE)

       IF (lrev) THEN   ! UNSORT
         DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
            ! GET ELEMENT VECTOR
            mhelp=i
            vechelp(i,1) = 0
            vechelp(i,2) = 0
            vechelp(i,3) = 0
            vechelp(i,4) = 0
            daccelementhelp(i,1) = 1
            daccelementhelp(i,2) = daccelementhelp(i,1)*vdim(1)
            daccelementhelp(i,3) = daccelementhelp(i,2)*vdim(2)
            daccelementhelp(i,4) = daccelementhelp(i,3)*vdim(3)
            vechelp(i,4) = (mhelp-1)/daccelementhelp(i,4)+1
            mhelp = mhelp - (vechelp(i,4)-1)*daccelementhelp(i,4)
            vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
            mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
            vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
            mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
            vechelp(i,1) = mhelp

            ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
            IF (axes(GORD_LON) > 0) THEN
               vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
            END IF
            IF (axes(GORD_LAT) > 0) THEN
               vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
            END IF
            IF (axes(GORD_LEV) > 0) THEN
               IF (sw_hyam) THEN
                  vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
               ELSE
                  vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
               END IF
            END IF
            ! COPY DATA ELEMENT
            poshelp = vechelp(i,1)
            dacchelp = 1
            dacchelp = dacchelp*vdim(1)
            poshelp = poshelp + dacchelp*(vechelp(i,2)-1)
            dacchelp = dacchelp*vdim(2)
            poshelp = poshelp + dacchelp*(vechelp(i,3)-1)
            dacchelp = dacchelp*vdim(3)
            poshelp = poshelp + dacchelp*(vechelp(i,4)-1)

            svar%dat%vb(poshelp) = var%dat%vb(i)
         END DO

       ELSE             ! SORT

         DO i=1, ulimit  ! LOOP OVER ALL ELEMENTS
            ! GET ELEMENT VECTOR
            mhelp=i
            vechelp(i,1) = 0
            vechelp(i,2) = 0
            vechelp(i,3) = 0
            vechelp(i,4) = 0
            daccelementhelp(i,1) = 1
            daccelementhelp(i,2) = daccelementhelp(i,1)*vdim(1)
            daccelementhelp(i,3) = daccelementhelp(i,2)*vdim(2)
            daccelementhelp(i,4) = daccelementhelp(i,3)*vdim(3)
            vechelp(i,4) = (mhelp-1)/daccelementhelp(i,4)+1
            mhelp = mhelp - (vechelp(i,4)-1)*daccelementhelp(i,4)
            vechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
            mhelp = mhelp - (vechelp(i,3)-1)*daccelementhelp(i,3)
            vechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
            mhelp = mhelp - (vechelp(i,2)-1)*daccelementhelp(i,2)
            vechelp(i,1) = mhelp

            ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
            IF (axes(GORD_LON) > 0) THEN
               vechelp(i,axes(GORD_LON)) = gx%lonm%dat%vi(vechelp(i,axes(GORD_LON)))
            END IF
            IF (axes(GORD_LAT) > 0) THEN
               vechelp(i,axes(GORD_LAT)) = gx%latm%dat%vi(vechelp(i,axes(GORD_LAT)))
            END IF
            IF (axes(GORD_LEV) > 0) THEN
               IF (sw_hyam) THEN
                  vechelp(i,axes(GORD_LEV)) = gx%hyam%dat%vi(vechelp(i,axes(GORD_LEV)))
               ELSE
                  vechelp(i,axes(GORD_LEV)) = gx%hybm%dat%vi(vechelp(i,axes(GORD_LEV)))
               END IF
            END IF
            ! COPY DATA ELEMENT
            poshelp = vechelp(i,1)
            dacchelp = 1
            dacchelp = dacchelp*vdim(1)
            poshelp = poshelp + dacchelp*(vechelp(i,2)-1)
            dacchelp = dacchelp*vdim(2)
            poshelp = poshelp + dacchelp*(vechelp(i,3)-1)
            dacchelp = dacchelp*vdim(3)
            poshelp = poshelp + dacchelp*(vechelp(i,4)-1)

            svar%dat%vb(i) = var%dat%vb(poshelp)
         END DO
       END IF

    END SELECT
    DEALLOCATE(daccelementhelp)
    DEALLOCATE(vechelp, STAT=status)
    NULLIFY(vechelp)

  CASE DEFAULT ! size(vdim).GT.4

    DO i=1, PRODUCT(vdim)  ! LOOP OVER ALL ELEMENTS
       ! GET ELEMENT VECTOR
       CALL ELEMENT(vdim,i,vec)
       ! CHANGE ELEMENT VECTOR ACCORDING TO SORT ORDER
       IF (axes(GORD_LON) > 0) THEN
          vec(axes(GORD_LON)) = gx%lonm%dat%vi(vec(axes(GORD_LON)))
       END IF
       IF (axes(GORD_LAT) > 0) THEN
          vec(axes(GORD_LAT)) = gx%latm%dat%vi(vec(axes(GORD_LAT)))
       END IF
       IF (axes(GORD_LEV) > 0) THEN
          IF (gx%hyam%dat%type /= VTYPE_UNDEF) THEN
             vec(axes(GORD_LEV)) = gx%hyam%dat%vi(vec(axes(GORD_LEV)))
          ELSE
             vec(axes(GORD_LEV)) = gx%hybm%dat%vi(vec(axes(GORD_LEV)))
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

  END SELECT

  ! CLEAN UP
  DEALLOCATE(vdim, STAT=status)
  CALL ERRMSG(substr,status,3)

END SUBROUTINE SORT_GEOHYBGRID_NCVAR
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE PACK_GEOHYBGRID_NCVAR(vi, dims, axes ,vo, reverse)

  IMPLICIT NONE

  ! I/O
  TYPE (ncvar)                  , INTENT(IN)    :: vi   ! input variable
  INTEGER, DIMENSION(:)         , INTENT(IN)    :: dims
  INTEGER                       , INTENT(IN)    :: axes(3)
  TYPE (ncvar)                  , INTENT(INOUT) :: vo   ! output variable
  LOGICAL, OPTIONAL             , INTENT(IN)    :: reverse

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'PACK_GEOHYBGRID_NCVAR'
  TYPE (ncvar)                       :: vu     ! unpacked variable
  TYPE (ncvar)                       :: vp     ! packed variable
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

  INTEGER, DIMENSION(:,:), POINTER     :: umvechelp  ! unpacked main pos. vector
  INTEGER                              :: ulimit,mhelp
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: daccelementhelp
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: pvechelp   ! packed position vector
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: urvechelp ! unpacked rest pos. vector
  INTEGER                              :: poshelp,dacchelp
  INTEGER, DIMENSION(:), ALLOCATABLE   :: jhelp
  LOGICAL                              :: lnotvec

  NULLIFY(umvec)

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
!qqqqq sollte vektorisierbar sein, sofern man fuer npdim-1 und vu%ndims-npdim+1 ein loop unrolling durchfuehren
!                                  kann. Vielleicht laeuft das programm auch nie in IF (vu%ndims >= npdim) THEN
!                                  rein.
  lnotvec = .FALSE.

  IF (vu%ndims >= npdim) THEN
    IF( (SIZE(umdim).GT.4) .OR. (SIZE(urdim).GT.2) .OR. (SIZE(pdim).GT.4) ) lnotvec = .TRUE.
  ELSE
    IF( (SIZE(umdim).GT.4) .OR. (SIZE(pdim).GT.4) ) lnotvec = .TRUE.
  ENDIF

  IF(lnotvec) THEN
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

  ELSE


    ulimit=PRODUCT(pdim)
    ALLOCATE(umvechelp(ulimit,SIZE(umdim)), STAT=status)
    SELECT CASE(SIZE(umdim))
    CASE(1)
      DO i=1, ulimit  ! == SIZE(v?%dat%v?) ! LOOP OVER ALL ELEMENTS
        umvechelp(i,1) = i
      ENDDO
    CASE(2)
      ALLOCATE(daccelementhelp(ulimit,2))
      DO i=1, ulimit  ! == SIZE(v?%dat%v?) ! LOOP OVER ALL ELEMENTS
        mhelp=i
        umvechelp(i,1) = 0
        umvechelp(i,2) = 0
        daccelementhelp(i,1) = 1
        daccelementhelp(i,2) = daccelementhelp(i,1)*umdim(1)
        umvechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
        mhelp = mhelp - (umvechelp(i,2)-1)*daccelementhelp(i,2)
        umvechelp(i,1) = mhelp
      ENDDO
    CASE(3)
      ALLOCATE(daccelementhelp(ulimit,3))
      DO i=1, ulimit  ! == SIZE(v?%dat%v?) ! LOOP OVER ALL ELEMENTS
        mhelp=i
        umvechelp(i,1) = 0
        umvechelp(i,2) = 0
        umvechelp(i,3) = 0
        daccelementhelp(i,1) = 1
        daccelementhelp(i,2) = daccelementhelp(i,1)*umdim(1)
        daccelementhelp(i,3) = daccelementhelp(i,2)*umdim(2)
        umvechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
        mhelp = mhelp - (umvechelp(i,3)-1)*daccelementhelp(i,3)
        umvechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
        mhelp = mhelp - (umvechelp(i,2)-1)*daccelementhelp(i,2)
        umvechelp(i,1) = mhelp
      ENDDO
    CASE(4)
      ALLOCATE(daccelementhelp(ulimit,4))
      DO i=1, ulimit  ! == SIZE(v?%dat%v?) ! LOOP OVER ALL ELEMENTS
        mhelp=i
        umvechelp(i,1) = 0
        umvechelp(i,2) = 0
        umvechelp(i,3) = 0
        umvechelp(i,4) = 0
        daccelementhelp(i,1) = 1
        daccelementhelp(i,2) = daccelementhelp(i,1)*umdim(1)
        daccelementhelp(i,3) = daccelementhelp(i,2)*umdim(2)
        daccelementhelp(i,4) = daccelementhelp(i,3)*umdim(3)
        umvechelp(i,4) = (mhelp-1)/daccelementhelp(i,4)+1
        mhelp = mhelp - (umvechelp(i,4)-1)*daccelementhelp(i,4)
        umvechelp(i,3) = (mhelp-1)/daccelementhelp(i,3)+1
        mhelp = mhelp - (umvechelp(i,3)-1)*daccelementhelp(i,3)
        umvechelp(i,2) = (mhelp-1)/daccelementhelp(i,2)+1
        mhelp = mhelp - (umvechelp(i,2)-1)*daccelementhelp(i,2)
        umvechelp(i,1) = mhelp
      ENDDO
    END SELECT
    DEALLOCATE(daccelementhelp)

    ALLOCATE(pvechelp(ulimit,npdim), STAT=status)
    DO j=1, npdim-1
      DO i=1, ulimit  ! == SIZE(v?%dat%v?) ! LOOP OVER ALL ELEMENTS
        pvechelp(i,j) = umvechelp(i,map(j))
      ENDDO
    ENDDO

    ALLOCATE(urvechelp(ulimit,vu%ndims-npdim+1), STAT=status)
    urvechelp(:,:) = 0
    IF (vu%ndims >= npdim) THEN
      DO j=1, vu%ndims-npdim+1
        DO i=1, ulimit  ! == SIZE(v?%dat%v?) ! LOOP OVER ALL ELEMENTS
          urvechelp(i,j) = umvechelp(i,rmap(j))
        ENDDO
      ENDDO

      SELECT CASE(SIZE(urdim))
      CASE(1)
        DO i=1, ulimit  ! == SIZE(v?%dat%v?) ! LOOP OVER ALL ELEMENTS
          pvechelp(i,npdim) = urvechelp(i,1)
        ENDDO
      CASE(2)
        DO i=1, ulimit  ! == SIZE(v?%dat%v?) ! LOOP OVER ALL ELEMENTS
          poshelp = urvechelp(i,1)
          dacchelp = 1
          dacchelp = dacchelp*urdim(1)
          poshelp = poshelp + dacchelp*(urvechelp(i,2)-1)
          pvechelp(i,npdim) = poshelp
        ENDDO
      END SELECT
    ELSE
      DO i=1, ulimit  ! == SIZE(v?%dat%v?) ! LOOP OVER ALL ELEMENTS
        pvechelp(i,npdim) = 1
      ENDDO
    ENDIF

    ALLOCATE(jhelp(ulimit), STAT=status)
    SELECT CASE(SIZE(pdim))
    CASE(1)
      DO i=1, ulimit  ! == SIZE(v?%dat%v?) ! LOOP OVER ALL ELEMENTS
        jhelp(i) = pvechelp(i,1)
      ENDDO
    CASE(2)
      DO i=1, ulimit  ! == SIZE(v?%dat%v?) ! LOOP OVER ALL ELEMENTS
        poshelp = pvechelp(i,1)
        dacchelp = 1
        dacchelp = dacchelp*pdim(1)
        poshelp = poshelp + dacchelp*(pvechelp(i,2)-1)
        jhelp(i) = poshelp
      ENDDO
    CASE(3)
      DO i=1, ulimit  ! == SIZE(v?%dat%v?) ! LOOP OVER ALL ELEMENTS
        poshelp = pvechelp(i,1)
        dacchelp = 1
        dacchelp = dacchelp*pdim(1)
        poshelp = poshelp + dacchelp*(pvechelp(i,2)-1)
        dacchelp = dacchelp*pdim(2)
        poshelp = poshelp + dacchelp*(pvechelp(i,3)-1)
        jhelp(i) = poshelp
      ENDDO
    CASE(4)
      DO i=1, ulimit  ! == SIZE(v?%dat%v?) ! LOOP OVER ALL ELEMENTS
        poshelp = pvechelp(i,1)
        dacchelp = 1
        dacchelp = dacchelp*pdim(1)
        poshelp = poshelp + dacchelp*(pvechelp(i,2)-1)
        dacchelp = dacchelp*pdim(2)
        poshelp = poshelp + dacchelp*(pvechelp(i,3)-1)
        dacchelp = dacchelp*pdim(3)
        poshelp = poshelp + dacchelp*(pvechelp(i,4)-1)
        jhelp(i) = poshelp
      ENDDO
    END SELECT

    SELECT CASE(vtype)
    CASE(VTYPE_INT)
      IF (lrev) THEN
        DO i=1, ulimit  ! == SIZE(v?%dat%v?) ! LOOP OVER ALL ELEMENTS
          vu%dat%vi(i) = vp%dat%vi(jhelp(i))
        ENDDO
      ELSE
        DO i=1, ulimit  ! == SIZE(v?%dat%v?) ! LOOP OVER ALL ELEMENTS
          vp%dat%vi(jhelp(i)) = vu%dat%vi(i)
        ENDDO
      END IF
    CASE(VTYPE_REAL)
      IF (lrev) THEN
        DO i=1, ulimit  ! == SIZE(v?%dat%v?) ! LOOP OVER ALL ELEMENTS
          vu%dat%vr(i) = vp%dat%vr(jhelp(i))
        ENDDO
      ELSE
        DO i=1, ulimit  ! == SIZE(v?%dat%v?) ! LOOP OVER ALL ELEMENTS
          vp%dat%vr(jhelp(i)) = vu%dat%vr(i)
        ENDDO
      END IF
    CASE(VTYPE_DOUBLE)
      IF (lrev) THEN
        DO i=1, ulimit  ! == SIZE(v?%dat%v?) ! LOOP OVER ALL ELEMENTS
          vu%dat%vd(i) = vp%dat%vd(jhelp(i))
        ENDDO
      ELSE
        DO i=1, ulimit  ! == SIZE(v?%dat%v?) ! LOOP OVER ALL ELEMENTS
          vp%dat%vd(jhelp(i)) = vu%dat%vd(i)
        ENDDO
      END IF
    CASE(VTYPE_CHAR)
      IF (lrev) THEN
        DO i=1, ulimit  ! == SIZE(v?%dat%v?) ! LOOP OVER ALL ELEMENTS
          vu%dat%vc(i) = vp%dat%vc(jhelp(i))
        ENDDO
      ELSE
        DO i=1, ulimit  ! == SIZE(v?%dat%v?) ! LOOP OVER ALL ELEMENTS
          vp%dat%vc(jhelp(i)) = vu%dat%vc(i)
        ENDDO
      END IF
    CASE(VTYPE_BYTE)
      IF (lrev) THEN
        DO i=1, ulimit  ! == SIZE(v?%dat%v?) ! LOOP OVER ALL ELEMENTS
          vu%dat%vb(i) = vp%dat%vb(jhelp(i))
        ENDDO
      ELSE
        DO i=1, ulimit  ! == SIZE(v?%dat%v?) ! LOOP OVER ALL ELEMENTS
          vp%dat%vb(jhelp(i)) = vu%dat%vb(i)
        ENDDO
      END IF
    END SELECT

    DEALLOCATE(jhelp)
    DEALLOCATE(pvechelp)
    DEALLOCATE(urvechelp)
    DEALLOCATE(umvechelp, STAT=status)
    NULLIFY(umvechelp)

  ENDIF

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
  TYPE (geohybgrid), INTENT(INOUT) :: gi      ! input grid
  TYPE (geohybgrid), INTENT(INOUT) :: go      ! output grid

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
SUBROUTINE BALANCE_GEOHYBGRID_TIME(gi, go, lint)

  IMPLICIT NONE

  ! I/O
  TYPE (geohybgrid),      INTENT(INOUT) :: gi, go
  LOGICAL,                INTENT(IN)    :: lint

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'BALANCE_GEOHYBGRID_TIME'
  INTEGER :: i
  INTEGER :: uii, uio

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

! ------------------------------------------------------------------
SUBROUTINE BALANCE_GEOHYBGRID_PS(gi, go, ranges)

  IMPLICIT NONE

  ! I/O
  TYPE (geohybgrid),          INTENT(INOUT) :: gi, go
  REAL(DP), DIMENSION(2,4,2), INTENT(IN)    :: ranges

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'BALANCE_GEOHYBGRID_PS'
  LOGICAL :: err
  LOGICAL :: linv   ! PRE-REGRID AVAILABLE go%ps

  err = (.NOT.QDEF_NCVAR(gi%ps)).AND.(.NOT.QDEF_NCVAR(go%ps))
  err = err .AND. ((QDEF_NCVAR(gi%hybi).OR.QDEF_NCVAR(gi%hybm)) .OR. &
                   (QDEF_NCVAR(go%hybi).OR.QDEF_NCVAR(go%hybm)))

  IF (err) THEN
     CALL RGMSG(substr, RGMLE, &
          'HYBRID-B-COEFFICIENTS NEED I_PS AND/OR G_PS IN NAMELIST !')
  END IF

  err = (.NOT.QDEF_NCVAR(gi%p0)).AND.(.NOT.QDEF_NCVAR(go%p0))
  err = err .AND. ((QDEF_NCVAR(gi%hyai).OR.QDEF_NCVAR(gi%hyam)) .OR. &
                   (QDEF_NCVAR(go%hyai).OR.QDEF_NCVAR(go%hyam)))

  IF (err) THEN
     CALL RGMSG(substr, RGMLE, &
          'HYBRID-A-COEFFICIENTS NEED I_P0 AND/OR G_P0 IN NAMELIST !')
  END IF

  ! RETURN, IF 2D
  err = (.NOT.QDEF_NCVAR(gi%p0)).AND.(.NOT.QDEF_NCVAR(gi%ps)).AND.     &
        (.NOT.QDEF_NCVAR(gi%hyai)).AND.(.NOT.QDEF_NCVAR(gi%hybi)).AND. &
        (.NOT.QDEF_NCVAR(gi%hyam)).AND.(.NOT.QDEF_NCVAR(gi%hybm))
  IF (err) THEN
     CALL RGMSG(substr, RGMLI, &
          'INPUT GRID IS 2-D! NO SURFACE PRESSURE REGRIDDING REQUIRED!')
     RETURN
  END IF

  ! BALANCE PO IN ANY CASE
  IF (QDEF_NCVAR(gi%p0).AND.(.NOT.QDEF_NCVAR(go%p0))) THEN
     CALL COPY_NCVAR(go%p0, gi%p0)
  ELSE
     IF (QDEF_NCVAR(go%p0).AND.(.NOT.QDEF_NCVAR(gi%p0))) THEN
        CALL COPY_NCVAR(gi%p0, go%p0)
     END IF
  END IF

  ! NOW PS ...

  ! CASE A: gi%ps AND go%ps BOTH AVAILABLE
  ! ADJUST TIME AXIS
  ! NOTE: THE TIME BALANCING (INCLUDING THAT FOR THE SURFACE PRESSURE
  !       TIME DIMENSION) IS PERFORMED IN
  !       SUBROUTINE BALANCE_GEOHYBGRID_TIME !

  ! CASE B: gi%ps AND go%ps BOTH UNAVAILABLE
  ! NOTHING TO DO !
  IF (.NOT.QDEF_NCVAR(gi%ps).AND.(.NOT.QDEF_NCVAR(go%ps))) THEN
     RETURN
  END IF

  ! SURFACE PRESSURE NEEDS TO BE PRE-REGRIDDED, IF NOT AVAILABLE
  ! CASE C: go%ps AVAILABLE, BUT ON WRONG HORIZONTAL GRID
  !
  linv = QDEF_NCVAR(go%ps).AND.( &
         ((.NOT.QDEF_NCVAR(go%lonm)).AND.(.NOT.QDEF_NCVAR(go%loni))).OR. &
         ((.NOT.QDEF_NCVAR(go%latm)).AND.(.NOT.QDEF_NCVAR(go%lati))) )
  !
  IF (QDEF_NCVAR(go%ps).AND.linv) THEN
     CALL RGMSG(substr, RGMLE,  &
          'REGRIDDING 3-D DISTRIBUTIONS', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'ONTO A DESTINATION SURFACE PRESSURE COORDINATE', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'USING (AN) INVARIANT HORIZONTAL DIMENSION(S)', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IS NOT POSSIBLE DUE TO A LACK OF INFORMATION', .false.)
     CALL RGMSG(substr, RGMLEC, &
          '(HORIZONTAL DESTINATION GRID)', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'FOR PRE-REGRIDDING THE DESTINATION SURFACE PRESSURE !', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'PLEASE PERFORM 2-D PRE-REGRIDDING OF SURFACE PRESSURE', .false.)
     CALL RGMSG(substr, RGMLEC, &
          'IN SEPARATE STEP!')
  END IF

  ! CASE D: gi%ps XOR go%ps NOT AVAILABLE
  IF (QDEF_NCVAR(gi%ps).AND.(.NOT.QDEF_NCVAR(go%ps))) THEN
     CALL REGRID_GEOHYBGRID_PS(gi, go, ranges(2,:,:))
  ELSE
     IF (QDEF_NCVAR(go%ps).AND.(.NOT.QDEF_NCVAR(gi%ps))) THEN
        CALL REGRID_GEOHYBGRID_PS(go, gi, ranges(1,:,:))
     END IF
  END IF

END SUBROUTINE BALANCE_GEOHYBGRID_PS
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE REGRID_GEOHYBGRID_PS(gi, go, ranges)

  IMPLICIT NONE

  ! I/O
  TYPE (geohybgrid),        INTENT(IN)    :: gi
  TYPE (geohybgrid),        INTENT(INOUT) :: go
  REAL(DP), DIMENSION(4,2), INTENT(IN)    :: ranges

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'REGRID_GEOHYBGRID_PS'
  TYPE (geohybgrid) :: gih, goh, gihs, gohs
  TYPE (geohybgrid) :: gix, gox
  INTEGER, DIMENSION(:),       POINTER :: dims      ! order of grid dimensions
  INTEGER                              :: axes(3)   ! grid axes of variable
  TYPE (axis)  , DIMENSION(:), POINTER :: ai, ao    ! source and dest. axes
  TYPE (ncvar)                         :: psi, pso
  TYPE (ncvar)                         :: psis, psos
  TYPE (ncvar)                         :: psip, psop
  TYPE (narray), DIMENSION(:), POINTER :: nao       ! re-gridder I/O
  REAL (DP),     DIMENSION(:), POINTER :: sovl, dovl
  INTEGER,       DIMENSION(:), POINTER :: rcnt
  LOGICAL                              :: ok
  INTEGER                              :: i
  INTEGER                              :: status

  ! INIT
  NULLIFY(dims)
  NULLIFY(ai)
  NULLIFY(ao)
  NULLIFY(nao)
  NULLIFY(sovl)
  NULLIFY(dovl)
  NULLIFY(rcnt)

  ! CREATE 2-D GRIDS
  CALL INIT_GEOHYBGRID(gih)
  CALL INIT_GEOHYBGRID(goh)
  !
  gih%file = TRIM(gi%file)
  gih%t    = gi%t
  !
  goh%file = TRIM(go%file)
  goh%t    = go%t
  !
  CALL COPY_NCVAR(gih%lonm, gi%lonm)
  CALL COPY_NCVAR(gih%loni, gi%loni)
  !
  CALL COPY_NCVAR(gih%latm, gi%latm)
  CALL COPY_NCVAR(gih%lati, gi%lati)
  !
  CALL COPY_NCVAR(gih%timem, gi%timem)
  CALL COPY_NCVAR(gih%timei, gi%timei)
  !
  CALL COPY_NCVAR(goh%lonm, go%lonm)
  CALL COPY_NCVAR(goh%loni, go%loni)
  !
  CALL COPY_NCVAR(goh%latm, go%latm)
  CALL COPY_NCVAR(goh%lati, go%lati)
  !
  CALL COPY_NCVAR(goh%timem, go%timem)
  CALL COPY_NCVAR(goh%timei, go%timei)

  CALL SORT_GEOHYBGRID(gih, gihs, gix)
  CALL SORT_GEOHYBGRID(goh, gohs, gox)

  CALL COMPLETE_GEOHYBGRID(gihs, ranges, gix)
  CALL COMPLETE_GEOHYBGRID(gohs, ranges, gox)

  CALL BALANCE_GEOHYBGRID(gihs, gohs)
  CALL BALANCE_GEOHYBGRID(gix, gox)

  CALL GEOHYBGRID_AXES(gihs, ai, gohs, ao, .false.)

  CALL CHECK_NCVAR_ON_GEOHYBGRID(gi%ps, gihs, dims, axes, ok)
  IF (.NOT.ok) THEN
     CALL RGMSG(substr, RGMLE, 'PS NOT GRID CONFORM !')
  END IF

  CALL COPY_NCVAR(psi, gi%ps)
  CALL SORT_GEOHYBGRID_NCVAR(psi, gix, axes, psis)
  CALL BALANCE_GEOHYBGRID_NCVAR(psis, axes, gohs, psos)

  CALL PACK_GEOHYBGRID_NCVAR(psis, dims, axes, psip)
  CALL BALANCE_GEOHYBGRID_NCVAR(psip, axes, gohs, psop)

  DEALLOCATE(dims, STAT=status)
  CALL ERRMSG(substr,status,1)
  NULLIFY(dims)

  CALL NREGRID( (/ psip%dat /), ai, ao, nao, (/ RG_INT /), sovl, dovl, rcnt)
  IF (IAND(MSGMODE, MSGMODE_VM) == MSGMODE_VM) THEN
     CALL NREGRID_STAT(ai, ao, sovl, dovl, (/ psip%dat /), nao, rcnt)
  END IF

  CALL COPY_NARRAY(psop%dat, nao(1))
  CALL CHECK_NCVAR_ON_GEOHYBGRID(psos, gohs, dims, axes, ok)

  IF (.NOT.ok) THEN
     CALL RGMSG(substr, RGMLE, 'PS NOT GRID CONFORM !')
  END IF

  CALL PACK_GEOHYBGRID_NCVAR(psop, dims, axes, psos, .true.)

  CALL SORT_GEOHYBGRID_NCVAR(psos, gox, axes, pso, .true.)

  CALL COPY_NCVAR(go%ps, pso)

  ! CLEAN UP
  DEALLOCATE(dims, STAT=status)
  CALL ERRMSG(substr,status,2)
  NULLIFY(dims)
  DEALLOCATE(sovl, dovl, rcnt, STAT=status)
  CALL ERRMSG(substr,status,3)
  NULLIFY(sovl, dovl, rcnt)
  CALL INIT_NCVAR(psi)
  CALL INIT_NCVAR(pso)
  CALL INIT_NCVAR(psis)
  CALL INIT_NCVAR(psos)
  CALL INIT_NCVAR(psip)
  CALL INIT_NCVAR(psop)
  CALL INIT_GEOHYBGRID(gih)
  CALL INIT_GEOHYBGRID(goh)
  CALL INIT_GEOHYBGRID(gihs)
  CALL INIT_GEOHYBGRID(gohs)
  CALL INIT_GEOHYBGRID(gix)
  CALL INIT_GEOHYBGRID(gox)
  DO i=1, SIZE(nao)
     CALL INIT_NARRAY(nao(i))
  END DO
  DEALLOCATE(nao, STAT=status)
  CALL ERRMSG(substr,status,4)
  NULLIFY(nao)
  DO i=1, SIZE(ai)
     CALL INIT_AXIS(ai(i))
     CALL INIT_AXIS(ao(i))
  END DO
  DEALLOCATE(ai, ao, STAT=status)
  CALL ERRMSG(substr,status,5)
  NULLIFY(ai,ao)

END SUBROUTINE REGRID_GEOHYBGRID_PS
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE BALANCE_GEOHYBGRID_NCVAR(vari, axes, go, varo)

  IMPLICIT NONE

  ! Note: go must already be 'balanced'

  ! I/O
  TYPE (ncvar)     , INTENT(IN)  :: vari    ! input variable
  INTEGER          , INTENT(IN)  :: axes(3) ! dim.no of lon, lat, lev
  TYPE (geohybgrid), INTENT(IN)  :: go      ! output grid
  TYPE (ncvar)     , INTENT(OUT) :: varo    ! output variable

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'BALANCE_GEOHYBGRID_NCVAR'
  INTEGER                            :: i
  INTEGER                            :: vtype
  INTEGER, DIMENSION(:), ALLOCATABLE :: dimvec  ! dimension vector
  INTEGER                            :: status

  ! INIT
  !
  CALL COPY_NCVAR(varo, vari)
  !
  ALLOCATE(dimvec(vari%ndims), STAT=status)
  CALL ERRMSG(substr,status,1)
  !
  vtype = vari%dat%type

  ! CHANGE DIMENSIONS
  DO i=1, vari%ndims
     dimvec(i) = vari%dim(i)%len
     IF (i == axes(GORD_LON)) THEN
        dimvec(i) = go%lonm%dim(1)%len
        varo%dim(i) = go%lonm%dim(1)
     END IF
     IF (i == axes(GORD_LAT)) THEN
        dimvec(i) = go%latm%dim(1)%len
        varo%dim(i) = go%latm%dim(1)
     END IF
     IF (i == axes(GORD_LEV)) THEN
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

END SUBROUTINE BALANCE_GEOHYBGRID_NCVAR
! ------------------------------------------------------------------

! ******************************************************************
END MODULE MESSY_NCREGRID_GEOHYB
! ******************************************************************

#endif
