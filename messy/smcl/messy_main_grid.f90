! ********************************************************************
! --------------------------------------------------------------------
MODULE MESSY_MAIN_GRID
! --------------------------------------------------------------------
! Author: Patrick Joeckel, MPICH, Mainz, June 2002
!         Bastian Kern,    MPICH, Mainz, July 2009
!         extended for curvilinear ocean grids by
!         Astrid  Kerkweg, UniMz, Mainz, 2012-2013
!         -> re-structured/expanded for more general GRID application
!         extended for horizontally unstructured grids by
!         Andreas Baumgaertner, DLR, 2015-2016
! ********************************************************************
!---------------------------------------------------------------------
!---------------------------------------------------------------------

  USE MESSY_MAIN_CONSTANTS_MEM, ONLY: dp
  USE MESSY_MAIN_GRID_NETCDF,   ONLY: t_ncvar, GRD_MAXSTRLEN                 &
                                    , NULL_DIMID, NULL_VARID                 &
                                    , NF90_FLOAT, VTYPE_UNDEF, VTYPE_REAL    &
                                    , INIT_NARRAY                            &
                                    ,               COPY_NCDIM               &
                                    , INIT_NCVAR,   COPY_NCVAR,   QCMP_NCVAR &
                                    , IMPORT_NCVAR, EXPORT_NCVAR, QDEF_NCVAR &
                                    , ADD_NCATT,    COPY_NCATT,   INIT_NCATT &
                                    , ERRMSG, PRINT_NCVAR, PRINT_NCATT       &
                                    , t_ncatt, NARRAYMINVAL, NARRAYMAXVAL

  IMPLICIT NONE
  !PRIVATE
  PUBLIC
  SAVE

  INTRINSIC :: TRIM, ASSOCIATED, PRESENT, NULL, ADJUSTL

  PRIVATE   :: TRIM, ASSOCIATED, PRESENT, NULL, ADJUSTL

  REAL(dp), PARAMETER             :: RGEMPTY = -999.0_dp

  ! Indicator to use the provided maximum pressure
  ! (required for regridding on height based grids with associated and thus
  !  changing pressure levels)
  REAL(dp), PARAMETER             :: RGMAX   = 999999.0_dp

  TYPE t_geohybgrid
     CHARACTER(LEN=GRD_MAXSTRLEN) :: name   = ''  ! grid name
     INTEGER                      :: ID     = -99 ! GRID ID
     CHARACTER(LEN=GRD_MAXSTRLEN) :: file   = ''  ! path/filename
     INTEGER                      :: t      =  1  ! time step
     TYPE(t_ncatt)                :: att          ! list of attributes
     REAL(dp), DIMENSION(4,2)     :: ranges = RGEMPTY
     INTEGER                      :: corners = 4  ! number of corners 
     ! 
     ! MINIUM / MAXIMUM GEOGRAPHICAL EXTENSION
     REAL(dp), DIMENSION(2,2)     :: minmaxlonlat = RGEMPTY
     INTEGER,  DIMENSION(2)       :: start        = -99
     INTEGER,  DIMENSION(2)       :: count        = -99
     INTEGER,  DIMENSION(2)       :: hdimids      = -99
     LOGICAL                      :: lonc         = .TRUE.
     TYPE (t_ncvar)               :: lonm, latm, hyam, hybm, timem !mid-layer
     TYPE (t_ncvar)               :: loni, lati, hyai, hybi, timei !interface
     TYPE (t_ncvar)               :: ps, p0
     ! curved grid variables
     ! geographical coordinates
     LOGICAL                      :: clonc  = .TRUE.
     TYPE (t_ncvar)               :: clonm, clatm, cloni, clati 
     ! rotated coordinates
     LOGICAL                      :: rlonc  = .TRUE.
     TYPE (t_ncvar)               :: rlonm, rlatm, rloni, rlati 
     ! rotated pole for rotated coordinates
     TYPE (t_ncvar)               :: pollon, pollat, polgam
     ! unstructured grid (CESM HOMME-SE spectral element grid)
     TYPE (t_ncvar)               :: ulonm, ulatm, uloni, ulati
     ! mask: only used in SCRIP   0: skipp / 1 : use
     TYPE (t_ncvar)               :: imask
     ! 3D pressure field:
     ! this is required to map pressure input fields to height grids
     ! or height input grids
     TYPE(t_ncvar)                :: pressm
     TYPE(t_ncvar)                :: pressi
   END TYPE t_geohybgrid

  TYPE t_geohybgrid_list
     TYPE(t_geohybgrid)               :: this
     TYPE(t_geohybgrid_list), POINTER :: next => NULL()
  END type t_geohybgrid_list

  TYPE(t_geohybgrid_list), POINTER    :: GEOHYBGRIDLIST => NULL()
  ! number of defined grids 
  INTEGER                             :: NGEOHYBGRID = 0 

  INTERFACE LOCATE_GEOHYBGRID
     MODULE PROCEDURE LOCATE_GEOHYBGRID_BY_ID
     MODULE PROCEDURE LOCATE_GEOHYBGRID_BY_NAME
  END INTERFACE
  !PUBLIC :: LOCATE_GEOHYBGRID

  INTERFACE NEW_GEOHYBGRID
     MODULE PROCEDURE NEW_GEOHYBGRID_BY_GRID
     MODULE PROCEDURE NEW_GEOHYBGRID_BY_COMPONENTS
  END INTERFACE
  !PUBLIC :: NEW_GEOHYBGRID

  PUBLIC :: COMPARE_TWO_GRIDS
  PUBLIC :: GEOHYBGRID_UPDATE_PRESS
  PUBLIC :: GEOHYBGRID_UPDATE_SURFPRESS

CONTAINS

! ------------------------------------------------------------------
SUBROUTINE INIT_GEOHYBGRID(grid)

  ! initialize geohybgrid grid

  IMPLICIT NONE

  ! I/O
  TYPE (t_geohybgrid), INTENT(INOUT)          :: grid

  grid%name         = ''
  grid%id           = -99
  grid%file         = ''
  grid%t            = 1

  grid%ranges       = RGEMPTY
  grid%corners      = 4
  grid%minmaxlonlat = RGEMPTY

  grid%start        = -99
  grid%count        = -99
  grid%hdimids      = -99

  grid%lonc         = .TRUE.
  grid%clonc        = .TRUE.
  grid%rlonc        = .TRUE.
  
  CALL INIT_NCATT(grid%att)
     
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

  CALL INIT_NCVAR(grid%clonm)
  CALL INIT_NCVAR(grid%clatm)
  CALL INIT_NCVAR(grid%cloni)
  CALL INIT_NCVAR(grid%clati)

  CALL INIT_NCVAR(grid%rlonm)
  CALL INIT_NCVAR(grid%rlatm)
  CALL INIT_NCVAR(grid%rloni)
  CALL INIT_NCVAR(grid%rlati)

  CALL INIT_NCVAR(grid%pollon)
  CALL INIT_NCVAR(grid%pollat)
  CALL INIT_NCVAR(grid%polgam)

  CALL INIT_NCVAR(grid%ulonm)
  CALL INIT_NCVAR(grid%ulatm)
  CALL INIT_NCVAR(grid%uloni)
  CALL INIT_NCVAR(grid%ulati)
  CALL INIT_NCVAR(grid%imask)
  CALL INIT_NCVAR(grid%pressm)
  CALL INIT_NCVAR(grid%pressi)

END SUBROUTINE INIT_GEOHYBGRID
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE COPY_GEOHYBGRID(gd, gs)

  ! COPY source grid (gs) to destination grid (gd)

  IMPLICIT NONE

  ! I/O
  TYPE (t_geohybgrid), INTENT(INOUT) :: gd  ! destination 
  TYPE (t_geohybgrid), INTENT(IN)    :: gs  ! source

  CALL INIT_GEOHYBGRID(gd)

  gd%name         = TRIM(gs%name)
  gd%id           = gs%id
  gd%file         = TRIM(gs%file)
  gd%t            = gs%t

  gd%ranges       = gs%ranges
  gd%corners      = gs%corners
  gd%minmaxlonlat = gs%minmaxlonlat

  gd%start        = gs%start
  gd%count        = gs%count
  gd%hdimids      = gs%hdimids

  gd%lonc         = gs%lonc 
  gd%clonc        = gs%clonc 
  gd%rlonc        = gs%rlonc 

  CALL COPY_NCATT(gd%att, gs%att)

  CALL COPY_NCVAR(gd%lonm,  gs%lonm)
  CALL COPY_NCVAR(gd%latm,  gs%latm)
  CALL COPY_NCVAR(gd%hyam,  gs%hyam)
  CALL COPY_NCVAR(gd%hybm,  gs%hybm)
  CALL COPY_NCVAR(gd%timem, gs%timem)

  CALL COPY_NCVAR(gd%loni,  gs%loni)
  CALL COPY_NCVAR(gd%lati,  gs%lati)
  CALL COPY_NCVAR(gd%hyai,  gs%hyai)
  CALL COPY_NCVAR(gd%hybi,  gs%hybi)
  CALL COPY_NCVAR(gd%timei, gs%timei)

  CALL COPY_NCVAR(gd%ps,    gs%ps)
  CALL COPY_NCVAR(gd%p0,    gs%p0)

  CALL COPY_NCVAR(gd%clonm, gs%clonm)
  CALL COPY_NCVAR(gd%clatm, gs%clatm)
  CALL COPY_NCVAR(gd%cloni, gs%cloni)
  CALL COPY_NCVAR(gd%clati, gs%clati)

  CALL COPY_NCVAR(gd%rlonm, gs%rlonm)
  CALL COPY_NCVAR(gd%rlatm, gs%rlatm)
  CALL COPY_NCVAR(gd%rloni, gs%rloni)
  CALL COPY_NCVAR(gd%rlati, gs%rlati)

  CALL COPY_NCVAR(gd%pollon,gs%pollon)
  CALL COPY_NCVAR(gd%pollat,gs%pollat)
  CALL COPY_NCVAR(gd%polgam,gs%polgam)

  CALL COPY_NCVAR(gd%ulonm,  gs%ulonm)
  CALL COPY_NCVAR(gd%ulatm,  gs%ulatm)
  CALL COPY_NCVAR(gd%uloni,  gs%uloni)
  CALL COPY_NCVAR(gd%ulati,  gs%ulati)
  CALL COPY_NCVAR(gd%imask,  gs%imask)
  CALL COPY_NCVAR(gd%pressm, gs%pressm)
  CALL COPY_NCVAR(gd%pressi, gs%pressi)

END SUBROUTINE COPY_GEOHYBGRID
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE IMPORT_GEOHYBGRID(grid, pstart, pcount)

  ! read geohybgrid from netcdf file
  ! NO CHECKING IN THIS ROUTINE !!!

  IMPLICIT NONE

  ! I/O
  TYPE (t_geohybgrid), INTENT(INOUT) :: grid
  INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL :: pstart
  INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL :: pcount

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER        :: substr = 'IMPORT_GEOHYBGRID'
  CHARACTER(LEN=GRD_MAXSTRLEN)       :: name   ! variable name
  CHARACTER(LEN=GRD_MAXSTRLEN)       :: file   ! netCDF filename
  INTEGER                            :: ustep  ! step along UNLIMITED DIM.
  INTEGER                            :: status
  INTEGER                            :: iostat
  INTEGER                            :: ndims
  INTEGER                            :: dimlen ! dimension length of ps
  INTEGER, DIMENSION(:), ALLOCATABLE :: dimvec ! dimension vector of ps
  REAL                               :: p      ! local pressure
  CHARACTER(LEN=GRD_MAXSTRLEN)       :: unit   ! pressure unit
  INTEGER                            :: setuid ! force unlimited ID

  
  INTEGER, DIMENSION(2)       :: zstart, zcount
  INTEGER, DIMENSION(2)       :: locstart, loccount
  REAL                        :: pollon, pollat, polgam
  INTEGER                     :: i, j
  CHARACTER(LEN=26)           :: press_string = 'constant surface pressure'
  INTEGER, DIMENSION(2)       :: hdids = -99

  IF (PRESENT(pstart)) THEN
     locstart = pstart
  ELSE
     locstart = 1
  ENDIF
  IF (PRESENT(pcount)) THEN
     loccount = pcount
  ELSE
     loccount = -99
  ENDIF

  ustep = grid%t
  file  = TRIM(grid%file)

  ! LONGITUDE
  IF (TRIM(grid%lonm%name) /= '') THEN
     zstart(1) = locstart(1)
     zstart(2) = 1 ! some values which is not used as lonm is 1-dim.
     zcount(1) = loccount(1)
     zcount(2) = -99 ! some values which is not used as lonm is 1-dim.
     hdids(1)  = grid%hdimids(1)
     hdids(2)  = -99
     name = TRIM(grid%lonm%name)
     CALL INIT_NCVAR(grid%lonm)
     IF (zcount(1) > 0 ) THEN
        CALL IMPORT_NCVAR(grid%lonm,  varname=name, file=file &
             , pstart=zstart, pcount=zcount, hdimids=hdids)
     ELSE
        CALL IMPORT_NCVAR(grid%lonm,  varname=name, file=file)
     END IF
  END IF
  IF (TRIM(grid%loni%name) /= '') THEN
     zstart(1) = locstart(1)
     zstart(2) = 1 ! some values which is not used as lonm is 1-dim.
     zcount(1) = loccount(1) + 1
     zcount(2) = -99 ! some values which is not used as lonm is 1-dim.
     hdids(1)  = grid%hdimids(1)
     hdids(2)  = -99
     name = TRIM(grid%loni%name)
     CALL INIT_NCVAR(grid%loni)
     IF (zcount(1) > 0 ) THEN
        CALL IMPORT_NCVAR(grid%loni,  varname=name, file=file &
             , pstart=zstart, pcount=zcount, hdimids=hdids)
     ELSE
        CALL IMPORT_NCVAR(grid%loni,  varname=name, file=file)
     END IF
  END IF

  IF (TRIM(grid%clonm%name) /= '') THEN
     zstart(:) = locstart(:)
     zcount(:) = loccount(:)
     hdids(:)  = grid%hdimids(:)
     name = TRIM(grid%clonm%name)
     CALL INIT_NCVAR(grid%clonm)
     IF (ANY(zcount < 0)) THEN
        CALL IMPORT_NCVAR(grid%clonm, varname=name, file=file )
     ELSE
        CALL IMPORT_NCVAR(grid%clonm, varname=name, file=file &
             , pstart=zstart, pcount=zcount, hdimids=hdids)
     END IF
  END IF
  IF (TRIM(grid%cloni%name) /= '') THEN
     zstart(:) = locstart(:)
     zcount(:) = loccount(:) + 1
     name = TRIM(grid%cloni%name)
     hdids(:)  = grid%hdimids(:)
     CALL INIT_NCVAR(grid%cloni)
     IF (ANY(zcount < 0 )) THEN
        CALL IMPORT_NCVAR(grid%cloni, varname=name, file=file)
     ELSE
        CALL IMPORT_NCVAR(grid%cloni, varname=name, file=file &
             , pstart=zstart, pcount=zcount, hdimids=hdids)
     END IF
  END IF

  IF (TRIM(grid%rlonm%name) /= '') THEN
     zstart(1) = locstart(1)
     zstart(2) = 1 ! some values which is not used as lonm is 1-dim.
     zcount(1) = loccount(1)
     zcount(2) = -99 ! some values which is not used as lonm is 1-dim.
     hdids(1)  = grid%hdimids(1)
     hdids(2)  = -99
     name = TRIM(grid%rlonm%name)
     CALL INIT_NCVAR(grid%rlonm)
     IF (zcount(1) > 0) THEN
        CALL IMPORT_NCVAR(grid%rlonm, varname=name, file=file &
             , pstart=zstart, pcount=zcount, hdimids=hdids)
     ELSE
        CALL IMPORT_NCVAR(grid%rlonm, varname=name, file=file)
     END IF
  END IF
  IF (TRIM(grid%rloni%name) /= '') THEN
     zstart(1) = locstart(1)
     zstart(2) = 1 ! some values which is not used as lonm is 1-dim.
     zcount(1) = loccount(1) + 1
     zcount(2) = -99 ! some values which is not used as lonm is 1-dim.
     hdids(1)  = grid%hdimids(1)
     hdids(2)  = -99
     name = TRIM(grid%rloni%name)
     CALL INIT_NCVAR(grid%rloni)
     IF (zcount(1) > 0) THEN
        CALL IMPORT_NCVAR(grid%rloni, varname=name, file=file &
             , pstart=zstart, pcount=zcount, hdimids=hdids)
     ELSE
        CALL IMPORT_NCVAR(grid%rloni, varname=name, file=file)
     END IF
  END IF

  ! LATITUDE
  IF (TRIM(grid%latm%name) /= '') THEN
     ! NOTE: latm is 1-dim. thus start and count are moved to index 1
     zstart(1) = locstart(2)
     zstart(2) = 1 ! some values which is not used as latm is 1-dim.
     zcount(1) = loccount(2)
     zcount(2) = -99 ! some values which is not used as latm is 1-dim.
     hdids(1)  = grid%hdimids(2)
     hdids(2)  = -99
     name = TRIM(grid%latm%name)
     CALL INIT_NCVAR(grid%latm)
     CALL IMPORT_NCVAR(grid%latm,  varname=name, file=file &
          , pstart=zstart, pcount=zcount, hdimids=hdids)
  END IF
  IF (TRIM(grid%lati%name) /= '') THEN
     zstart(1) = locstart(2)
     zstart(2) = 1 ! some values which is not used as lonm is 1-dim.
     zcount(1) = loccount(2) + 1
     zcount(2) = -99 ! some values which is not used as lonm is 1-dim.
     hdids(1)  = grid%hdimids(2)
     hdids(2)  = -99
     name = TRIM(grid%lati%name)
     CALL INIT_NCVAR(grid%lati)
     CALL IMPORT_NCVAR(grid%lati,  varname=name, file=file &
          , pstart=zstart, pcount=zcount, hdimids=hdids)
  END IF
  IF (TRIM(grid%clatm%name) /= '') THEN
     zstart(:) = locstart(:)
     zcount(:) = loccount(:)
     hdids(:)  = grid%hdimids(:)
     name = TRIM(grid%clatm%name)
     CALL INIT_NCVAR(grid%clatm)
     CALL IMPORT_NCVAR(grid%clatm, varname=name, file=file &
          , pstart=zstart, pcount=zcount, hdimids=hdids)
  END IF
  IF (TRIM(grid%clati%name) /= '') THEN
     zstart(:) = locstart(:)
     zcount(:) = loccount(:) + 1
     hdids(:)  = grid%hdimids(:)
     name = TRIM(grid%clati%name)
     CALL INIT_NCVAR(grid%clati)
     CALL IMPORT_NCVAR(grid%clati, varname=name, file=file &
          , pstart=zstart, pcount=zcount, hdimids=hdids)
  END IF
  IF (TRIM(grid%rlatm%name) /= '') THEN
     zstart(1) = locstart(2)
     zstart(2) = 1 ! some values which is not used as lonm is 1-dim.
     zcount(1) = loccount(2)
     zcount(2) = -99 ! some values which is not used as lonm is 1-dim.
     hdids(1)  = grid%hdimids(2)
     hdids(2)  = -99
     name = TRIM(grid%rlatm%name)
     CALL INIT_NCVAR(grid%rlatm)
     CALL IMPORT_NCVAR(grid%rlatm, varname=name, file=file &
          , pstart=zstart, pcount=zcount, hdimids=hdids)
  END IF
  IF (TRIM(grid%rlati%name) /= '') THEN
     zstart(1) = locstart(2)
     zstart(2) = 1 ! some values which is not used as lonm is 1-dim.
     zcount(1) = loccount(2) + 1
     zcount(2) = -99 ! some values which is not used as lonm is 1-dim.
     hdids(1)  = grid%hdimids(2)
     hdids(2)  = -99
      name = TRIM(grid%rlati%name)
     CALL INIT_NCVAR(grid%rlati)
     CALL IMPORT_NCVAR(grid%rlati, varname=name, file=file &
          , pstart=zstart, pcount=zcount, hdimids=hdids)
  END IF

  ! TIME
  IF (TRIM(grid%timem%name) /= '') THEN
     name = TRIM(grid%timem%name)
     CALL INIT_NCVAR(grid%timem)
     CALL IMPORT_NCVAR(grid%timem, ustep=ustep, varname=name, file=file)
     ! ASSOCIATE ALL TIME AXES WITH UNLIMITED ID
     grid%timem%dim(1)%fuid = .true.
     grid%timem%uid = grid%timem%dim(1)%id
  END IF
  IF (TRIM(grid%timei%name) /= '') THEN
     name = TRIM(grid%timei%name)
     CALL INIT_NCVAR(grid%timei)
     CALL IMPORT_NCVAR(grid%timei, varname=name, file=file)
  END IF

  ! HYA
  IF (TRIM(grid%hyam%name) /= '') THEN
     name = TRIM(grid%hyam%name)
     CALL INIT_NCVAR(grid%hyam)
     CALL IMPORT_NCVAR(grid%hyam, varname=name, file=file)
  END IF
  IF (TRIM(grid%hyai%name) /= '') THEN
     name = TRIM(grid%hyai%name)
     CALL INIT_NCVAR(grid%hyai)
     CALL IMPORT_NCVAR(grid%hyai, varname=name, file=file)
  END IF

  ! HYB
  IF (TRIM(grid%hybm%name) /= '') THEN
     name = TRIM(grid%hybm%name)
     CALL INIT_NCVAR(grid%hybm)
     CALL IMPORT_NCVAR(grid%hybm, varname=name, file=file)
  END IF
  IF (TRIM(grid%hybi%name) /= '') THEN
     name = TRIM(grid%hybi%name)
     CALL INIT_NCVAR(grid%hybi)
     CALL IMPORT_NCVAR(grid%hybi, varname=name, file=file)
  END IF

  ! SURFACE PRESSURE
  IF (TRIM(grid%ps%name) /= '') THEN
     ! CHECK FOR REAL VALUE
     READ(grid%ps%name,*,IOSTAT=iostat) p
     IF (iostat == 0) THEN    ! VALUE of PS
        ! below the important information is lost
        READ(grid%ps%name,*,IOSTAT=iostat) p, unit
        CALL INIT_NCVAR(grid%ps)
        grid%ps%name  = TRIM('ps_const')
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
        CALL ADD_NCATT(grid%ps, 'long_name',vs=press_string)
        ! ... UNITS ATTRIBUT
        unit = ''

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
        ! test for second call from reduce grid
        IF (.NOT. QDEF_NCVAR(grid%ps) &
             .OR.( QDEF_NCVAR(grid%ps) .AND. TRIM(grid%ps%name)/='ps_const') &
             ) THEN
           name = TRIM(grid%ps%name)
           IF (ASSOCIATED(grid%timem%dim)) THEN
              setuid = grid%timem%dim(1)%id
           ELSE
              setuid = NULL_DIMID
           END IF
           zstart(:) = locstart(:)
           zcount(:) = loccount(:)
           CALL INIT_NCVAR(grid%ps)
           CALL IMPORT_NCVAR(grid%ps, ustep=ustep, varname=name, file=file &
                , setuid=setuid, pstart=zstart, pcount=zcount              &
                , hdimids=grid%hdimids)
        ELSE
           ! CASE CONSTANT PS DEFINED IN FIRST READ => CHANGE DIMENSIONS
           ! 1. CHECK attributes for constant surface pressure
           ! 2. memorize const pressur
           ! 3. initialise and re-define constant pressure field
           att_loop: DO i = 1 , grid%ps%natts        
              IF (TRIM(grid%ps%att(i)%name)=='long_name') THEN
                 DO j = 1,grid%ps%att(i)%len
                    IF (grid%ps%att(i)%dat%vc(j) /= press_string(j:j)) &
                    CYCLE att_loop
                 END DO
                 EXIT att_loop
              END IF
           END DO att_loop
           IF (i <= grid%ps%natts) THEN
              ! constant pressure
              p = grid%ps%dat%vr(1)
              !CALL INIT_NARRAY(grid%ps%dat) 
              ndims = grid%ps%ndims
              ALLOCATE(dimvec(grid%ps%ndims))
              IF (loccount(1)>0) THEN
                 dimvec(1) = loccount(1)
                 grid%ps%dim(1)%len = loccount(1)
              ELSE
                 dimvec(1) = grid%ps%dim(1)%len
              ENDIF
              IF (loccount(2) > 0) THEN
                 dimvec(2) = loccount(2)
                 grid%ps%dim(2)%len = loccount(2)
              ELSE
                 dimvec(2) = grid%ps%dim(2)%len
              ENDIF
              CALL INIT_NARRAY(grid%ps%dat,ndims,dimvec,VTYPE_REAL )
              grid%ps%dat%vr(:) = p
              DEALLOCATE(dimvec)
           ENDIF
        ENDIF 
     END IF ! PS IS VALUE OR NAME
  END IF

  IF (TRIM(grid%p0%name) /= '') THEN
     ! CHECK FOR REAL VALUE
     READ(grid%p0%name,*,IOSTAT=iostat) p
     IF (iostat == 0) THEN    ! VALUE of P0
        ! below the important information is lost
        unit = ''
        READ(grid%p0%name,*,IOSTAT=iostat) p, unit 
        CALL INIT_NCVAR(grid%p0)
        grid%p0%name  = TRIM('p0')
        grid%p0%xtype = NF90_FLOAT
        !
        ! NO DIMENSIONS !!!
        grid%p0%ndims = 1
        ALLOCATE(grid%p0%dim(1), STAT=status)
        CALL ERRMSG(substr,status,4)
        grid%p0%dim(1)%name  = TRIM(grid%p0%name)//'_dim'
        grid%p0%dim(1)%id    = NULL_DIMID
        grid%p0%dim(1)%len   = 1
        grid%p0%dim(1)%fuid  = .false.
        grid%p0%dim(1)%varid = NULL_VARID
        !
        ! ATTRIBUTES ...
        ! ... LONGNAME ATTRIBUTE
        CALL ADD_NCATT(grid%p0, 'long_name', vs='reference pressure')
        ! ... UNITS ATTRIBUT
        IF (TRIM(unit) /= '') THEN
           CALL ADD_NCATT(grid%p0, 'units' ,vs=TRIM(unit))
        END IF
        !
        ! DATA
        CALL INIT_NARRAY(grid%p0%dat, 1, (/ 1 /), VTYPE_REAL)
        grid%p0%dat%vr(1) = p
     ELSE  ! NAME of P0
        ! test for second call from reduce grid
        IF (.NOT. QDEF_NCVAR(grid%p0)) THEN
           name = TRIM(grid%p0%name)
           IF (ASSOCIATED(grid%timem%dim)) THEN
              setuid = grid%timem%dim(1)%id
           ELSE
              setuid = NULL_DIMID
           END IF
           zstart(:) = locstart(:)
           zcount(:) = loccount(:)
           CALL INIT_NCVAR(grid%p0)
           CALL IMPORT_NCVAR(grid%p0, varname=name, file=file, setuid=setuid&
                , pstart=zstart, pcount=zcount, hdimids=hdids)
        ENDIF 
     END IF
  END IF

  ! 3D pressure
  ! a) interfaces 3D pressure
  IF (TRIM(grid%pressi%name) /= '') THEN
     zstart(:) = locstart(:)
     zcount(:) = loccount(:)
     name = TRIM(grid%pressi%name)
     CALL INIT_NCVAR(grid%pressi)
     CALL IMPORT_NCVAR(grid%pressi, varname=name, file=file &
          , pstart=zstart, pcount=zcount, hdimids=hdids)
  END IF

  ! b) mids 3D pressure
   IF (TRIM(grid%pressm%name) /= '') THEN
      zstart(:) = locstart(:)
      zcount(:) = loccount(:)
      name = TRIM(grid%pressm%name)
      CALL INIT_NCVAR(grid%pressm)
      CALL IMPORT_NCVAR(grid%pressm, varname=name, file=file &
           , pstart=zstart, pcount=zcount, hdimids=hdids)
  END IF

  ! READ ROTATED POLE COORDINATES
  IF (TRIM(grid%pollon%name) /= '') THEN
     ! CHECK FOR REAL VALUE
     READ(grid%pollon%name,*,IOSTAT=iostat) pollon
     IF (iostat == 0) THEN    ! VALUE of pollon
        CALL INIT_NCVAR(grid%pollon)
        grid%pollon%name  = TRIM('pollon')
        grid%pollon%xtype = NF90_FLOAT
        !
        ! NO DIMENSIONS !!!
        grid%pollon%ndims = 1
        ALLOCATE(grid%pollon%dim(1), STAT=status)
        CALL ERRMSG(substr,status,4)
        grid%pollon%dim(1)%name  = TRIM(grid%pollon%name)//'_dim'
        grid%pollon%dim(1)%id    = NULL_DIMID
        grid%pollon%dim(1)%len   = 1
        grid%pollon%dim(1)%fuid  = .false.
        grid%pollon%dim(1)%varid = NULL_VARID
        !
        ! ATTRIBUTES ...
        ! ... LONGNAME ATTRIBUTE
        CALL ADD_NCATT(grid%pollon, 'long_name', vs='grid_north_pole_longitude')
        ! ... UNITS ATTRIBUT
        unit = ''
        READ(grid%pollon%name,*,IOSTAT=iostat) pollon,unit
        IF (TRIM(unit) /= '') THEN
           CALL ADD_NCATT(grid%pollon, 'units' ,vs=TRIM(unit))
        END IF
        !
        ! DATA
        CALL INIT_NARRAY(grid%pollon%dat, 1, (/ 1 /), VTYPE_REAL)
        grid%pollon%dat%vr(1) = pollon
     ELSE  ! NAME of pollon
        ! test for second call from reduce grid
        IF (.NOT. QDEF_NCVAR(grid%pollon)) THEN
           name = TRIM(grid%pollon%name)
           IF (ASSOCIATED(grid%timem%dim)) THEN
              setuid = grid%timem%dim(1)%id
           ELSE
              setuid = NULL_DIMID
           END IF
           zstart(:) = locstart(:)
           zcount(:) = loccount(:)
           hdids = grid%hdimids
           CALL INIT_NCVAR(grid%pollon)
           CALL IMPORT_NCVAR(grid%pollon, varname=name, file=file &
                , setuid=setuid , pstart=zstart, pcount=zcount    &
                , hdimids=hdids)
        ENDIF 
     END IF
  END IF

  IF (TRIM(grid%pollat%name) /= '') THEN
     ! CHECK FOR REAL VALUE
     READ(grid%pollat%name,*,IOSTAT=iostat) pollat
     IF (iostat == 0) THEN    ! VALUE of pollat
        CALL INIT_NCVAR(grid%pollat)
        grid%pollat%name  = TRIM('pollat')
        grid%pollat%xtype = NF90_FLOAT
        !
        ! NO DIMENSIONS !!!
        grid%pollat%ndims = 1
        ALLOCATE(grid%pollat%dim(1), STAT=status)
        CALL ERRMSG(substr,status,4)
        grid%pollat%dim(1)%name  = TRIM(grid%pollat%name)//'_dim'
        grid%pollat%dim(1)%id    = NULL_DIMID
        grid%pollat%dim(1)%len   = 1
        grid%pollat%dim(1)%fuid  = .false.
        grid%pollat%dim(1)%varid = NULL_VARID
        !
        ! ATTRIBUTES ...
        ! ... LONGNAME ATTRIBUTE
        CALL ADD_NCATT(grid%pollat, 'long_name', vs='grid_north_pole_latitude')
        ! ... UNITS ATTRIBUT
        unit = ''
        READ(grid%pollat%name,*,IOSTAT=iostat) pollat,unit
        IF (TRIM(unit) /= '') THEN
           CALL ADD_NCATT(grid%pollat, 'units' ,vs=TRIM(unit))
        END IF
        !
        ! DATA
        CALL INIT_NARRAY(grid%pollat%dat, 1, (/ 1 /), VTYPE_REAL)
        grid%pollat%dat%vr(1) = pollat
     ELSE  ! NAME of pollat
        ! test for second call from reduce grid
        IF (.NOT. QDEF_NCVAR(grid%pollat)) THEN
           name = TRIM(grid%pollat%name)
           IF (ASSOCIATED(grid%timem%dim)) THEN
              setuid = grid%timem%dim(1)%id
           ELSE
              setuid = NULL_DIMID
           END IF
           zstart(:) = locstart(:)
           zcount(:) = loccount(:)
           hdids = grid%hdimids
           CALL INIT_NCVAR(grid%pollat)
           CALL IMPORT_NCVAR(grid%pollat, varname=name, file=file &
                , setuid=setuid , pstart=zstart, pcount=zcount    &
                , hdimids=hdids)
        ENDIF 
     END IF
  END IF

  IF (TRIM(grid%polgam%name) /= '') THEN
     ! CHECK FOR REAL VALUE
     READ(grid%polgam%name,*,IOSTAT=iostat) polgam
     IF (iostat == 0) THEN    ! VALUE of polgam
        CALL INIT_NCVAR(grid%polgam)
        grid%polgam%name  = TRIM('polgam')
        grid%polgam%xtype = NF90_FLOAT
        !
        ! NO DIMENSIONS !!!
        grid%polgam%ndims = 1
        ALLOCATE(grid%polgam%dim(1), STAT=status)
        CALL ERRMSG(substr,status,4)
        grid%polgam%dim(1)%name  = TRIM(grid%polgam%name)//'_dim'
        grid%polgam%dim(1)%id    = NULL_DIMID
        grid%polgam%dim(1)%len   = 1
        grid%polgam%dim(1)%fuid  = .false.
        grid%polgam%dim(1)%varid = NULL_VARID
        !
        ! ATTRIBUTES ...
        ! ... LONGNAME ATTRIBUTE
        CALL ADD_NCATT(grid%polgam, 'long_name', vs='grid_north_pole_angle')
        ! ... UNITS ATTRIBUT
        unit = ''
        READ(grid%polgam%name,*,IOSTAT=iostat) polgam,unit
        IF (TRIM(unit) /= '') THEN
           CALL ADD_NCATT(grid%polgam, 'units' ,vs=TRIM(unit))
        END IF
        !
        ! DATA
        CALL INIT_NARRAY(grid%polgam%dat, 1, (/ 1 /), VTYPE_REAL)
        grid%polgam%dat%vr(1) = polgam
     ELSE  ! NAME of polgam
        ! test for second call from reduce grid
        IF (.NOT. QDEF_NCVAR(grid%polgam)) THEN
           name = TRIM(grid%polgam%name)
           IF (ASSOCIATED(grid%timem%dim)) THEN
              setuid = grid%timem%dim(1)%id
           ELSE
              setuid = NULL_DIMID
           END IF
           zstart(:) = locstart(:)
           zcount(:) = loccount(:)
           hdids = grid%hdimids
           CALL INIT_NCVAR(grid%polgam)
           CALL IMPORT_NCVAR(grid%polgam, varname=name, file=file &
                , setuid=setuid , pstart=zstart, pcount=zcount    &
                , hdimids=hdids)
        ENDIF 
     END IF
  END IF

  ! READ UNSTRUCTURED GRID COORDINATES
  IF (TRIM(grid%ulatm%name) /= '') THEN
     zstart(1) = locstart(1)
     zstart(2) = 1
     zcount(1) = loccount(1)
     zcount(2) = 1
     hdids = grid%hdimids
     name = TRIM(grid%ulatm%name)
     CALL INIT_NCVAR(grid%ulatm)
     CALL IMPORT_NCVAR(grid%ulatm,  varname=name, file=file &
          , pstart=zstart, pcount=zcount, hdimids=hdids)
  END IF
  IF (TRIM(grid%ulonm%name) /= '') THEN
     zstart(1) = locstart(1)
     zstart(2) = 1
     zcount(1) = loccount(1)
     zcount(2) = 1
     hdids = grid%hdimids
     name = TRIM(grid%ulonm%name)
     CALL INIT_NCVAR(grid%ulonm)
     CALL IMPORT_NCVAR(grid%ulonm,  varname=name, file=file &
          , pstart=zstart, pcount=zcount, hdimids=hdids)
  END IF
  IF (TRIM(grid%ulati%name) /= '') THEN
     zstart(1) = locstart(1)
     zstart(2) = 1
     zcount(1) = loccount(1) + 1
     zcount(2) = loccount(2) + 1
     hdids = grid%hdimids
     name = TRIM(grid%ulati%name)
     CALL INIT_NCVAR(grid%ulati)
     CALL IMPORT_NCVAR(grid%ulati,  varname=name, file=file &
          , pstart=zstart, pcount=zcount, hdimids=hdids)
  END IF
  IF (TRIM(grid%uloni%name) /= '') THEN
     zstart(1) = locstart(1)
     zstart(2) = 1
     zcount(1) = loccount(1) + 1
     zcount(2) = loccount(2) + 1
     hdids = grid%hdimids
     name = TRIM(grid%uloni%name)
     CALL INIT_NCVAR(grid%uloni)
     CALL IMPORT_NCVAR(grid%uloni,  varname=name, file=file &
          , pstart=zstart, pcount=zcount, hdimids=hdids)
     grid%corners = grid%uloni%dim(2)%len
  END IF
  IF (TRIM(grid%imask%name) /= '') THEN
     IF ( QDEF_NCVAR(grid%lonm) .OR. QDEF_NCVAR(grid%rlonm))  THEN
        zstart(1) = locstart(1)
        zstart(2) = 1 ! some values which is not used as lonm is 1-dim.
        zcount(1) = loccount(1)
        zcount(2) = -99 ! some values which is not used as lonm is 1-dim.
        hdids(1) = grid%hdimids(1)
        hdids(2) = -99
     ELSE IF (QDEF_NCVAR(grid%ulonm)) THEN
        zstart(1) = locstart(1)
        zstart(2) = 1
        zcount(1) = loccount(1)
        zcount(2) = loccount(2)
        hdids = grid%hdimids
     END IF
     name = TRIM(grid%imask%name)
     CALL INIT_NCVAR(grid%imask)
     IF (zcount(1) > 0 ) THEN
        CALL IMPORT_NCVAR(grid%imask,  varname=name, file=file &
             , pstart=zstart, pcount=zcount, hdimids=hdids)
     ELSE
        CALL IMPORT_NCVAR(grid%imask,  varname=name, file=file)
     END IF
  END IF

END SUBROUTINE IMPORT_GEOHYBGRID
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE EXPORT_GEOHYBGRID(grid)

  ! write geohybgrid to netcdf file
  ! NO CHECKING IN THIS ROUTINE !!!

  IMPLICIT NONE

  ! I/O
  TYPE (t_geohybgrid), INTENT(INOUT) :: grid

  ! LOCAL
  CHARACTER(LEN=GRD_MAXSTRLEN) :: file     ! netCDF filename
  INTEGER                      :: ustep    ! step along UNLIMITED DIM.
  CHARACTER(LEN=GRD_MAXSTRLEN) :: timename ! name of unlimited variable

  timename = ''
  ustep = grid%t
  file = TRIM(grid%file)

  ! LONGITUDE
  IF (TRIM(grid%lonm%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%lonm, file=file)
  END IF
  IF (TRIM(grid%loni%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%loni, file=file)
  END IF
  IF (TRIM(grid%clonm%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%clonm, file=file)
  END IF
  IF (TRIM(grid%cloni%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%cloni, file=file)
  END IF
  IF (TRIM(grid%ulonm%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%ulonm, file=file)
  END IF
!!$  IF (TRIM(grid%uloni%name) /= '') THEN
!!$     CALL EXPORT_NCVAR(grid%uloni, file=file)
!!$  END IF

!!$  IF (TRIM(grid%rlonm%name) /= '') THEN
!!$     CALL EXPORT_NCVAR(grid%rlonm, file=file)
!!$  END IF
!!$  IF (TRIM(grid%rloni%name) /= '') THEN
!!$     CALL EXPORT_NCVAR(grid%rloni, file=file)
!!$  END IF

  ! LATITUDE
  IF (TRIM(grid%latm%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%latm, file=file)
  END IF
  IF (TRIM(grid%lati%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%lati, file=file)
  END IF

  IF (TRIM(grid%clatm%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%clatm, file=file)
  END IF
  IF (TRIM(grid%clati%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%clati, file=file)
  END IF
  IF (TRIM(grid%ulatm%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%ulatm, file=file)
  END IF
!!$  IF (TRIM(grid%ulati%name) /= '') THEN
!!$     CALL EXPORT_NCVAR(grid%ulati, file=file)
!!$  END IF

!!$  IF (TRIM(grid%rlatm%name) /= '') THEN
!!$     CALL EXPORT_NCVAR(grid%rlatm, file=file)
!!$  END IF
!!$  IF (TRIM(grid%rlati%name) /= '') THEN
!!$     CALL EXPORT_NCVAR(grid%rlati, file=file)
!!$  END IF

  ! TIME
  IF (TRIM(grid%timem%name) /= '') THEN
     grid%timem%ustep = ustep
     CALL unify_ulimit_name(timename, grid%timem)
     CALL EXPORT_NCVAR(grid%timem, file=file)
  END IF

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

  IF (TRIM(grid%pressi%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%pressi, file=file)
  END IF
  IF (TRIM(grid%pressm%name) /= '') THEN
     CALL EXPORT_NCVAR(grid%pressm, file=file)
  END IF
  ! SURFACE PRESSURE
  IF (TRIM(grid%ps%name) /= '') THEN
     grid%ps%ustep = ustep
     CALL unify_ulimit_name(timename, grid%ps)
     CALL EXPORT_NCVAR(grid%ps, file=file)
  END IF

  ! REFERENCE PRESSURE
  IF (TRIM(grid%p0%name) /= '') THEN
     CALL unify_ulimit_name(timename, grid%p0)
     CALL EXPORT_NCVAR(grid%p0, file=file)
  END IF

CONTAINS

SUBROUTINE UNIFY_ULIMIT_NAME(tname, var)

  IMPLICIT NONE

  ! I/O
  CHARACTER(LEN=GRD_MAXSTRLEN), INTENT(INOUT) :: tname
  TYPE(t_ncvar),                INTENT(INOUT) :: var
  ! LOCAL
  INTEGER   :: ix

  DO ix = 1, var%ndims
     IF (var%dim(ix)%fuid) THEN
        IF (TRIM(tname) /= '') THEN
           var%dim(ix)%name = TRIM(tname)
        ELSE
           tname = TRIM(var%dim(ix)%name)
        END IF
        EXIT
     END IF
  END DO

END SUBROUTINE UNIFY_ULIMIT_NAME

END SUBROUTINE EXPORT_GEOHYBGRID
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE PRINT_GEOHYBGRID(grid, str, verbose_level)

  ! print geohybgrid content

  TYPE(t_geohybgrid), INTENT(IN) :: grid
  CHARACTER(LEN=*),   INTENT(IN) :: str
  INTEGER, INTENT(IN), OPTIONAL  :: verbose_level

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'PRINT_GEOHYBGRID'
  INTEGER                     :: iverb 

  ! verbose_level
  ! 100: print all
  !  80: all but no attributes
  iverb = 100
  IF (PRESENT(verbose_level) ) iverb = verbose_level


  write (0,*) substr, ' ', str, ' NAME / ID: ', grid%name, grid%id
  write (0,*) substr, ' ', str, ' file: ', grid%file, 'TIMESTEP', grid%t
  write (0,*) substr, ' ', str, ' RANGES / CORNERS: ', grid%ranges, grid%corners
  write (0,*) substr, ' ', str, ' MIN/MAX LON/LAT: ', grid%minmaxlonlat

  IF (iverb > 80) CALL PRINT_NCATT(grid%att,   str)

  IF (QDEF_NCVAR(grid%lonm)) THEN
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  write (0,*) substr, ' lonm ', str, ''
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  CALL PRINT_NCVAR(grid%lonm,  str, verbose_level=iverb)
  ENDIF
  IF (QDEF_NCVAR(grid%latm)) THEN
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  write (0,*) substr, ' latm ', str, ''
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  CALL PRINT_NCVAR(grid%latm,  str, verbose_level=iverb)
  ENDIF
  IF (QDEF_NCVAR(grid%hyam)) THEN
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  write (0,*) substr, ' hyam ', str, ''
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  CALL PRINT_NCVAR(grid%hyam,  str, verbose_level=iverb)
  ENDIF
  IF (QDEF_NCVAR(grid%hybm)) THEN
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  write (0,*) substr, ' hybm ', str, ''
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  CALL PRINT_NCVAR(grid%hybm,  str, verbose_level=iverb)
  ENDIF
  IF (QDEF_NCVAR(grid%timem)) THEN
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  write (0,*) substr, ' timem ', str, ''
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  CALL PRINT_NCVAR(grid%timem, str, verbose_level=iverb)
  ENDIF

  IF (QDEF_NCVAR(grid%loni)) THEN
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  write (0,*) substr, ' loni ', str, ''
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  CALL PRINT_NCVAR(grid%loni,  str, verbose_level=iverb)
  ENDIF
  IF (QDEF_NCVAR(grid%lati)) THEN
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  write (0,*) substr, ' lati ', str, ''
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  CALL PRINT_NCVAR(grid%lati,  str, verbose_level=iverb)
  ENDIF
  IF (QDEF_NCVAR(grid%hyai)) THEN
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  write (0,*) substr, ' hyai ', str, ''
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  CALL PRINT_NCVAR(grid%hyai,  str, verbose_level=iverb)
  ENDIF
  IF (QDEF_NCVAR(grid%hybi)) THEN
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  write (0,*) substr, ' hybi ', str, ''
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  CALL PRINT_NCVAR(grid%hybi,  str, verbose_level=iverb)
  ENDIF
  IF (QDEF_NCVAR(grid%timei)) THEN
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  write (0,*) substr, ' timei ', str, ''
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  CALL PRINT_NCVAR(grid%timei, str, verbose_level=iverb)
  ENDIF

  IF (QDEF_NCVAR(grid%ps)) THEN
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  write (0,*) substr, ' ps ', str, ''
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  CALL PRINT_NCVAR(grid%ps,    str, verbose_level=iverb)
  ENDIF
  IF (QDEF_NCVAR(grid%p0)) THEN
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  write (0,*) substr, ' p0  ', str, ''
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  CALL PRINT_NCVAR(grid%p0,    str, verbose_level=iverb)
  ENDIF

  IF (QDEF_NCVAR(grid%pressi)) THEN
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  write (0,*) substr, ' pressi  ', str, ''
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  CALL PRINT_NCVAR(grid%pressi,    str, verbose_level=iverb)
  ENDIF
  IF (QDEF_NCVAR(grid%pressm)) THEN
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  write (0,*) substr, ' pressm  ', str, ''
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  CALL PRINT_NCVAR(grid%pressm,    str, verbose_level=iverb)
  ENDIF

  IF (QDEF_NCVAR(grid%clonm)) THEN
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  write (0,*) substr, ' clonm', str, ''
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  CALL PRINT_NCVAR(grid%clonm, str, verbose_level=iverb)
  ENDIF
  IF (QDEF_NCVAR(grid%clatm)) THEN
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  write (0,*) substr, ' clatm ', str, ''
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  CALL PRINT_NCVAR(grid%clatm, str, verbose_level=iverb)
  ENDIF
  IF (QDEF_NCVAR(grid%cloni)) THEN
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  write (0,*) substr, ' cloni ', str, ''
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  CALL PRINT_NCVAR(grid%cloni, str, verbose_level=iverb)
  ENDIF
  IF (QDEF_NCVAR(grid%clati)) THEN
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  write (0,*) substr, ' clati ', str, ''
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  CALL PRINT_NCVAR(grid%clati, str, verbose_level=iverb)
  ENDIF

  IF (QDEF_NCVAR(grid%rlonm)) THEN
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  write (0,*) substr, ' rlonm ', str, ''
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  CALL PRINT_NCVAR(grid%rlonm, str, verbose_level=iverb)
  ENDIF
  IF (QDEF_NCVAR(grid%rlatm)) THEN
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  write (0,*) substr, ' rlatm ', str, ''
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  CALL PRINT_NCVAR(grid%rlatm, str, verbose_level=iverb)
  ENDIF
  IF (QDEF_NCVAR(grid%rloni)) THEN
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  write (0,*) substr, ' rloni ', str, ''
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  CALL PRINT_NCVAR(grid%rloni, str, verbose_level=iverb)
  ENDIF
  IF (QDEF_NCVAR(grid%rlati)) THEN
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  write (0,*) substr, ' rlati ', str, ''
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  CALL PRINT_NCVAR(grid%rlati, str, verbose_level=iverb)
  ENDIF

  IF (QDEF_NCVAR(grid%ulonm)) THEN
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  write (0,*) substr, ' ulonm ', str, ''
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  CALL PRINT_NCVAR(grid%ulonm, str, verbose_level=iverb)
  ENDIF
  IF (QDEF_NCVAR(grid%ulatm)) THEN
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  write (0,*) substr, ' ulatm ', str, ''
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  CALL PRINT_NCVAR(grid%ulatm, str, verbose_level=iverb)
  ENDIF
  IF (QDEF_NCVAR(grid%uloni)) THEN
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  write (0,*) substr, ' uloni ', str, ''
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  CALL PRINT_NCVAR(grid%uloni, str, verbose_level=iverb)
  ENDIF
  IF (QDEF_NCVAR(grid%ulati)) THEN
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  write (0,*) substr, ' ulati ', str, ''
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  CALL PRINT_NCVAR(grid%ulati, str, verbose_level=iverb)
  ENDIF

  IF (QDEF_NCVAR(grid%pollon)) THEN
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  write (0,*) substr, ' pollon ', str, ''
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  CALL PRINT_NCVAR(grid%pollon, str, verbose_level=iverb)
  ENDIF
  IF (QDEF_NCVAR(grid%pollat)) THEN
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  write (0,*) substr, ' pollat ', str, ''
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  CALL PRINT_NCVAR(grid%pollat, str, verbose_level=iverb)
  ENDIF
  IF (QDEF_NCVAR(grid%polgam)) THEN
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  write (0,*) substr, ' polgam ', str, ''
  write (0,*) 'PPPPPPPPPPPPPPPPPPPPPPP'
  CALL PRINT_NCVAR(grid%polgam, str, verbose_level=iverb)
  ENDIF

END SUBROUTINE PRINT_GEOHYBGRID
! ------------------------------------------------------------------

! ===========================================================================
! ===========================================================================
! ===========================================================================
! ===========================================================================
SUBROUTINE LOCATE_GEOHYBGRID_BY_ID(status, ID, pgrid, grid)

  ! search grid according to the given ID
  ! OUTPUT: 
  !   - grid:  the grid itself
  !   - pgrid: a pointer to the grid

  IMPLICIT NONE

  ! I/O
  INTEGER, INTENT(OUT)                    :: status
  INTEGER, INTENT(IN)                     :: ID
  TYPE(t_geohybgrid),   POINTER, OPTIONAL :: pgrid      
  TYPE(t_geohybgrid),            OPTIONAL :: grid      

  ! LOCAL
  TYPE(t_geohybgrid_list), POINTER :: gi => NULL()
  TYPE(t_geohybgrid_list), POINTER :: ge => NULL()
  LOGICAL                          :: lexists
  !CHARACTER(LEN=*), PARAMETER     :: substr = 'locate_geohybgrid_id'

  lexists = .FALSE.
  status  = 1

  ! CHECKS
  IF (id <= 0) THEN
     status = 2007  ! INVALID GRID ID
     RETURN
  END IF
  !
  IF (.NOT. ASSOCIATED(GEOHYBGRIDLIST)) THEN
     status = 2004  ! GRID (ID) DOES NOT EXIST
     RETURN
  END IF

  gi => GEOHYBGRIDLIST
  DO
     IF (.NOT. ASSOCIATED(gi)) EXIT
     IF (id == gi%this%id) THEN
        lexists = .TRUE.
        EXIT
     END IF
     ge => gi
     gi => ge%next
  END DO

  IF (lexists) THEN
     IF (PRESENT(pgrid)) THEN
        pgrid => gi%this
     ENDIF
     IF (PRESENT(grid)) THEN
        CALL COPY_GEOHYBGRID(grid, gi%this)
     ENDIF
  ELSE
     status = 2003  ! GRID DOES NOT EXIST
     RETURN
  END IF
  NULLIFY(gi)
  NULLIFY(ge)

  status = 0
  
END SUBROUTINE LOCATE_GEOHYBGRID_BY_ID
! ===========================================================================

! ===========================================================================
SUBROUTINE LOCATE_GEOHYBGRID_BY_NAME(status, name, pgrid, grid, ID)

  ! search grid according to the given name
  ! OUTPUT: 
  !   - grid:  the grid itself
  !   - pgrid: a pointer to the grid
  !   - ID:    the grid ID 

  USE MESSY_MAIN_CONSTANTS_MEM, ONLY: strlen_medium

  IMPLICIT NONE

  INTRINSIC :: ADJUSTL, ASSOCIATED, LEN_TRIM, TRIM

  ! I/O
  INTEGER,            INTENT(OUT)             :: status
  CHARACTER(LEN=*),   INTENT(IN)              :: name
  TYPE(t_geohybgrid), POINTER,       OPTIONAL :: pgrid
  TYPE(t_geohybgrid), INTENT(INOUT), OPTIONAL :: grid
  INTEGER,            INTENT(OUT), OPTIONAL :: ID    

  ! LOCAL
  TYPE(t_geohybgrid_list), POINTER :: ai  => NULL()
  TYPE(t_geohybgrid_list), POINTER :: ae  => NULL()
  LOGICAL                          :: lexists

  ! INIT
  lexists = .FALSE.
  NULLIFY(pgrid)
  CALL INIT_GEOHYBGRID(grid)

  ! CHECKS
  IF (TRIM(name) == '') THEN
     status = 2010 ! GRID NAME IS EMPTY
     RETURN
  END IF
  IF (LEN_TRIM(ADJUSTL(name)) > STRLEN_MEDIUM) THEN
     status = 2001  ! GRID NAME TOO LONG
     RETURN
  END IF
  !
  IF (.NOT. ASSOCIATED(GEOHYBGRIDLIST)) THEN
     status = 2003  ! GRID (NAME) DOES NOT EXIST
     RETURN
  END IF

  ! CHECK, IF IT EXISTS
  ai => GEOHYBGRIDLIST
  DO
     IF (.NOT. ASSOCIATED(ai)) EXIT
     IF (TRIM(ADJUSTL(name)) == TRIM(ai%this%name)) THEN
        lexists = .TRUE.
        EXIT
     END IF
     ae => ai
     ai => ae%next
  END DO

  IF (lexists) THEN
     IF (PRESENT(pgrid)) THEN
        pgrid => ai%this
     ENDIF
     IF (PRESENT(grid)) THEN
        CALL COPY_GEOHYBGRID(grid, ai%this)
     ENDIF
     IF (PRESENT(ID)) THEN
        ID = ai%this%ID
     ENDIF
  ELSE
     status = 2003  ! GRID DOES NOT EXIST
     RETURN
  END IF

  NULLIFY(ai)
  NULLIFY(ae)
  status = 0

END SUBROUTINE LOCATE_GEOHYBGRID_BY_NAME
! ===========================================================================

! ===========================================================================
SUBROUTINE NEW_GEOHYBGRID_BY_GRID( status, id, grid)

  ! define new geohybgrid
  ! INPUT: 
  !  -  grid: the new grid define as structure
  ! OUTPUT:
  !  - id:    the id of the newly defined grid

  IMPLICIT NONE

  ! I/O
  INTEGER,            INTENT(OUT)        :: status
  INTEGER,            INTENT(OUT)        :: id
  TYPE(t_geohybgrid), INTENT(IN)         :: grid

  CALL new_geohybgrid_by_components(status, id, TRIM(grid%name), grid%corners &
               , grid%lonm,  grid%latm,  grid%hyam,  grid%hybm, grid%timem    &
               , grid%loni,  grid%lati,  grid%hyai,  grid%hybi, grid%timei    &
               , grid%clonm, grid%clatm, grid%cloni, grid%clati               & 
               , grid%rlonm, grid%rlatm, grid%rloni, grid%rlati               &
               , grid%ps,    grid%p0,    grid%file,  grid%t,    grid%ranges   &
               , grid%minmaxlonlat,    grid%pollon, grid%pollat,grid%polgam   &
               , grid%ulonm, grid%ulatm, grid%uloni, grid%ulati               & 
               , grid%lonc,  grid%clonc, grid%rlonc , grid%imask              &
               , grid%pressi, grid%pressm                                     &
               )

END SUBROUTINE NEW_GEOHYBGRID_BY_GRID
! ===========================================================================

! ===========================================================================
SUBROUTINE NEW_GEOHYBGRID_BY_COMPONENTS( status, id, name, corners     &
                                       , lonm, latm, hyam, hybm, timem &
                                       , loni, lati, hyai, hybi, timei &
                                       , clonm, clatm, cloni, clati    & 
                                       , rlonm, rlatm, rloni, rlati    &
                                       , ps,    p0,    file,  t        &
                                       , ranges, minmaxlonlat          &
                                       , pollon, pollat, polgam        &     
                                       , ulonm, ulatm, uloni, ulati    &
                                       , lonc, clonc, rlonc, imask     &
                                       , pressi, pressm                &
                                       )

  ! define new geohybgrid
  ! INPUT: 
  !  -  components of geohybgrid
  ! OUTPUT:
  !  - id:    the id of the newly defined grid

  USE MESSY_MAIN_TOOLS,     ONLY: int2str

  IMPLICIT NONE

  ! I/O
  INTEGER,                  INTENT(OUT)          :: status
  INTEGER,                  INTENT(OUT)          :: id
  CHARACTER(LEN=*),         INTENT(IN), OPTIONAL :: name 
  ! mid-layer
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: lonm
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: latm
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: hyam
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: hybm
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: timem
  ! interface layer
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: loni
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: lati
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: hyai
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: hybi
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: timei
  ! curvi-linear
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: clonm
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: clatm
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: cloni
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: clati
  ! rotated curvi-linear
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: rlonm
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: rlatm
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: rloni
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: rlati
  ! pressure
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: ps
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: p0
  ! rotated north pole
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: pollon
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: pollat
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: polgam
  ! unstructured grid
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: ulonm
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: ulatm
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: uloni
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: ulati

  LOGICAL,                  INTENT(IN), OPTIONAL :: lonc
  LOGICAL,                  INTENT(IN), OPTIONAL :: clonc
  LOGICAL,                  INTENT(IN), OPTIONAL :: rlonc
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: imask
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: pressi 
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: pressm

  ! if required file name for imported grid
  CHARACTER(LEN=*),         INTENT(IN), OPTIONAL :: file
  ! if required time step in file of imported grid
  INTEGER,                  INTENT(IN), OPTIONAL :: t
  INTEGER,                  INTENT(IN), OPTIONAL :: corners
  REAL(dp), DIMENSION(4,2), INTENT(IN), OPTIONAL :: ranges
  REAL(dp), DIMENSION(2,2), INTENT(IN), OPTIONAL :: minmaxlonlat

  ! LOCAL
  TYPE(t_geohybgrid_list), POINTER :: gi => NULL()
  TYPE(t_geohybgrid_list), POINTER :: ge => NULL()
  TYPE(t_geohybgrid),      POINTER :: lgrid
  INTEGER                          :: cid ! grid id returned by compare
  CHARACTER(LEN=GRD_MAXSTRLEN)     :: inpname = ' '     
  LOGICAL                          :: lname   = .FALSE.
  INTEGER, SAVE                    :: COUNT   = 0       
  CHARACTER(LEN=5)                 :: cstr = ''

  status = 2009 

  IF (PRESENT(name)) THEN
     lname   = .TRUE.
     inpname = name
  ELSE
     lname   = .FALSE.
     inpname = ''
  ENDIF

  ! GOTO END OF LIST
  gi => GEOHYBGRIDLIST
  DO
     IF (.NOT. ASSOCIATED(gi)) EXIT

     lgrid => gi%this
     ! do not compare pressi3d, as it is subject to change
     CALL COMPARE_TO_GRID(cid, lgrid, lonm,  latm,  hyam,  hybm, timem  &
                                    , loni,  lati,  hyai,  hybi, timei  &
                                    , clonm, clatm, cloni, clati        & 
                                    , rlonm, rlatm, rloni, rlati        & 
                                    , ps,    p0,    file,  t,    ranges &
                                    , minmaxlonlat, pollon, pollat, polgam &
                                    , ulonm, ulatm, uloni, ulati        &
                                    , corners, lonc, clonc, rlonc, imask)
     IF (cid > 0) THEN
        IF (.NOT. lname .OR. TRIM(ADJUSTL(inpname)) == '') THEN
           ID   =  cid
           status = 01
           RETURN
        END IF
     ENDIF
     IF (TRIM(gi%this%name) == TRIM(ADJUSTL(name))) THEN
        ! grid are not equal, but name exists
        status = 2009
     ENDIF
     ge => gi
     gi => gi%next
  END DO

!***************

  ! ADD NEW
  ALLOCATE(gi)
  NULLIFY(gi%next)
  IF (.NOT. ASSOCIATED(GEOHYBGRIDLIST)) THEN
     GEOHYBGRIDLIST => gi         ! SET POINTER TO FIRST GRID
  ELSE
     ge%next => gi                ! SET NEXT POINTER OF LAST GRID
     !                            ! TO NEW GRID
  END IF

  CALL INIT_GEOHYBGRID(gi%this)
  ! SET VALUES
  IF (lname .AND. TRIM(ADJUSTL(name)) /= '') THEN
     gi%this%name = TRIM(ADJUSTL(name))
  ELSE
     ! GENERATE GENERIC NAME
     COUNT = COUNT + 1
     CALL int2str(cstr, COUNT, '0')
     gi%this%name = 'GENERICGIRD'//cstr
  ENDIF
  ! COUNT AND SET ID
  NGEOHYBGRID = NGEOHYBGRID + 1
  gi%this%id  = NGEOHYBGRID

  ! DIMENSION GRID ... 
  IF (PRESENT(lonm))  CALL COPY_NCVAR(gi%this%lonm , lonm)
  IF (PRESENT(latm))  CALL COPY_NCVAR(gi%this%latm , latm)
  IF (PRESENT(hyam))  CALL COPY_NCVAR(gi%this%hyam , hyam)
  IF (PRESENT(hybm))  CALL COPY_NCVAR(gi%this%hybm , hybm)
  IF (PRESENT(timem)) CALL COPY_NCVAR(gi%this%timem, timem)

  IF (PRESENT(loni))  CALL COPY_NCVAR(gi%this%loni , loni)
  IF (PRESENT(lati))  CALL COPY_NCVAR(gi%this%lati , lati)
  IF (PRESENT(hyai))  CALL COPY_NCVAR(gi%this%hyai , hyai)
  IF (PRESENT(hybi))  CALL COPY_NCVAR(gi%this%hybi , hybi)
  IF (PRESENT(timei)) CALL COPY_NCVAR(gi%this%timei, timei)

  IF (PRESENT(clonm)) CALL COPY_NCVAR(gi%this%clonm, clonm)
  IF (PRESENT(clatm)) CALL COPY_NCVAR(gi%this%clatm, clatm)
  IF (PRESENT(cloni)) CALL COPY_NCVAR(gi%this%cloni, cloni)
  IF (PRESENT(clati)) CALL COPY_NCVAR(gi%this%clati, clati)

  IF (PRESENT(rlonm)) CALL COPY_NCVAR(gi%this%rlonm, rlonm)
  IF (PRESENT(rlatm)) CALL COPY_NCVAR(gi%this%rlatm, rlatm)
  IF (PRESENT(rloni)) CALL COPY_NCVAR(gi%this%rloni, rloni)
  IF (PRESENT(rlati)) CALL COPY_NCVAR(gi%this%rlati, rlati)

  IF (PRESENT(ps))    CALL COPY_NCVAR(gi%this%ps   , ps)
  IF (PRESENT(p0))    CALL COPY_NCVAR(gi%this%p0   , p0)

  IF (PRESENT(pollon))CALL COPY_NCVAR(gi%this%pollon, pollon)
  IF (PRESENT(pollat))CALL COPY_NCVAR(gi%this%pollat, pollat)
  IF (PRESENT(polgam))CALL COPY_NCVAR(gi%this%polgam, polgam)

  IF (PRESENT(ulonm)) CALL COPY_NCVAR(gi%this%ulonm, ulonm)
  IF (PRESENT(ulatm)) CALL COPY_NCVAR(gi%this%ulatm, ulatm)
  IF (PRESENT(uloni)) CALL COPY_NCVAR(gi%this%uloni, uloni)
  IF (PRESENT(ulati)) CALL COPY_NCVAR(gi%this%ulati, ulati)

  IF (PRESENT(lonc))     gi%this%lonc             = lonc
  IF (PRESENT(clonc))    gi%this%clonc            = clonc
  IF (PRESENT(rlonc))    gi%this%rlonc            = rlonc
  IF (PRESENT(imask))   CALL COPY_NCVAR(gi%this%imask,  imask)
  IF (PRESENT(pressi))  CALL COPY_NCVAR(gi%this%pressi, pressi)
  IF (PRESENT(pressm))  CALL COPY_NCVAR(gi%this%pressm, pressm)

  IF (PRESENT(file))         gi%this%file         = file
  IF (PRESENT(t))            gi%this%t            = t
  IF (PRESENT(corners))      gi%this%corners      = corners
  IF (PRESENT(ranges))       gi%this%ranges       = ranges
  IF (PRESENT(minmaxlonlat)) gi%this%minmaxlonlat = minmaxlonlat

  ID     = gi%this%id 

  NULLIFY(gi)
  NULLIFY(ge)
  NULLIFY(lgrid)

  status = 0

END SUBROUTINE NEW_GEOHYBGRID_BY_COMPONENTS
! ===========================================================================

! ===========================================================================
SUBROUTINE COMPARE_TO_GRID(id, grid, lonm,  latm,  hyam,  hybm, timem     &
                                   , loni,  lati,  hyai,  hybi, timei     &
                                   , clonm, clatm, cloni, clati           & 
                                   , rlonm, rlatm, rloni, rlati           & 
                                   , ps,    p0,    file,  t,    ranges    &
                                   , minmaxlonlat, pollon, pollat, polgam &
                                   , ulonm, ulatm, uloni, ulati           &
                                   , corners, lonc, clonc, rlonc, imask)

  ! compare grid components to given grid
  ! if all input grid components are equal to the grids components
  ! ID is set to the grid id.
  ! OUTPUT: ID 
  !      id = -99      if components are not equal
  !      id = grid%id  otherwise
  ! NOTE: pressi3d is a helper field,  required to enable regridding from 
  !       pressure based 
  !       coordinated to height coordinates. Thereofre pressi3d is subject to 
  !       change in each grid trafo time step. As it is not part of the grid
  !       definition itself, it is not compared !

  IMPLICIT NONE

  ! I/O
  INTEGER,            INTENT(OUT)                :: id
  TYPE(t_geohybgrid), POINTER                    :: grid

  ! mid-layer   
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: lonm
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: latm
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: hyam
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: hybm
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: timem
  ! interface layer
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: loni
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: lati
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: hyai
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: hybi
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: timei
  ! curvi-linear
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: clonm
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: clatm
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: cloni
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: clati
  ! rotated curvi-linear
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: rlonm
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: rlatm
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: rloni
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: rlati
  ! pressure
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: ps
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: p0
  ! rotated north pole
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: pollon
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: pollat
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: polgam
  ! if required file name for imported grid
  CHARACTER(LEN=*),         INTENT(IN), OPTIONAL :: file
  ! if required time step in file of imported grid
  INTEGER,                  INTENT(IN), OPTIONAL :: t
  INTEGER,                  INTENT(IN), OPTIONAL :: corners
  REAL(dp), DIMENSION(4,2), INTENT(IN), OPTIONAL :: ranges
  REAL(dp), DIMENSION(2,2), INTENT(IN), OPTIONAL :: minmaxlonlat
  ! unstructured 
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: ulonm
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: ulatm
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: uloni
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: ulati
  LOGICAL,                  INTENT(IN), OPTIONAL :: lonc
  LOGICAL,                  INTENT(IN), OPTIONAL :: clonc
  LOGICAL,                  INTENT(IN), OPTIONAL :: rlonc
  TYPE(t_ncvar),            INTENT(IN), OPTIONAL :: imask

  ! LOCAL
  INTEGER :: i,j

  ! INITIALIZE
  id = -99

  IF (PRESENT(lonm)) THEN
     IF ((.NOT. QDEF_NCVAR(grid%lonm) .AND. QDEF_NCVAR(lonm)) &
         .OR. (QDEF_NCVAR(grid%lonm) .AND. .NOT. QDEF_NCVAR(lonm))) RETURN
     IF  (QDEF_NCVAR(grid%lonm) .AND. QDEF_NCVAR(lonm)) THEN
        IF (.NOT. QCMP_NCVAR(grid%lonm, lonm)) RETURN
     END IF
  ELSE
     IF (QDEF_NCVAR(grid%lonm)) RETURN
  ENDIF

  IF (PRESENT(latm)) THEN
     IF ((.NOT. QDEF_NCVAR(grid%latm) .AND. QDEF_NCVAR(latm)) &
         .OR. (QDEF_NCVAR(grid%latm) .AND. .NOT. QDEF_NCVAR(latm))) RETURN
     IF (QDEF_NCVAR(grid%latm) .AND.  QDEF_NCVAR(latm)) THEN
        IF (.NOT. QCMP_NCVAR(grid%latm, latm)) RETURN
     END IF
  ELSE
     IF (QDEF_NCVAR(grid%latm)) RETURN
  ENDIF

  IF (PRESENT(hyam)) THEN
     IF ((.NOT. QDEF_NCVAR(grid%hyam) .AND. QDEF_NCVAR(hyam)) &
         .OR. (QDEF_NCVAR(grid%hyam) .AND. .NOT. QDEF_NCVAR(hyam))) RETURN
     IF (QDEF_NCVAR(grid%hyam) .AND. QDEF_NCVAR(hyam)) THEN
        IF (.NOT. QCMP_NCVAR(grid%hyam, hyam)) RETURN
     ENDIF
  ELSE
     IF (QDEF_NCVAR(grid%hyam)) RETURN
  ENDIF

  IF (PRESENT(hybm)) THEN
     IF ((.NOT. QDEF_NCVAR(grid%hybm) .AND. QDEF_NCVAR(hybm)) &
         .OR. (QDEF_NCVAR(grid%hybm) .AND. .NOT. QDEF_NCVAR(hybm))) RETURN
     IF (QDEF_NCVAR(grid%hybm) .AND. QDEF_NCVAR(hybm)) THEN
        IF (.NOT. QCMP_NCVAR(grid%hybm, hybm)) RETURN
     ENDIF
  ELSE
     IF (QDEF_NCVAR(grid%hybm)) RETURN
  ENDIF

  IF (PRESENT(timem)) THEN
     IF ((.NOT. QDEF_NCVAR(grid%timem) .AND. QDEF_NCVAR(timem)) &
         .OR. (QDEF_NCVAR(grid%timem) .AND. .NOT. QDEF_NCVAR(timem))) RETURN
     IF (QDEF_NCVAR(grid%timem) .AND. QDEF_NCVAR(timem)) THEN
        IF (.NOT. QCMP_NCVAR(grid%timem, timem)) RETURN
     ENDIF
  ELSE
     IF (QDEF_NCVAR(grid%timem)) RETURN
  ENDIF

  ! INTERFACE
  IF (PRESENT(loni)) THEN
     IF ((.NOT. QDEF_NCVAR(grid%loni) .AND. QDEF_NCVAR(loni)) &
         .OR. (QDEF_NCVAR(grid%loni) .AND. .NOT. QDEF_NCVAR(loni))) RETURN
     IF (QDEF_NCVAR(grid%loni) .AND. QDEF_NCVAR(loni)) THEN
        IF (.NOT. QCMP_NCVAR(grid%loni, loni)) RETURN
     ENDIF
  ELSE
     IF (QDEF_NCVAR(grid%loni)) RETURN
  ENDIF

  IF (PRESENT(lati)) THEN
     IF ((.NOT. QDEF_NCVAR(grid%lati) .AND. QDEF_NCVAR(lati)) &
         .OR. (QDEF_NCVAR(grid%lati) .AND. .NOT. QDEF_NCVAR(lati))) RETURN
     IF (QDEF_NCVAR(grid%lati) .AND. QDEF_NCVAR(lati)) THEN
        IF (.NOT. QCMP_NCVAR(grid%lati, lati)) RETURN
     ENDIF
  ELSE
     IF (QDEF_NCVAR(grid%lati)) RETURN
  ENDIF

  IF (PRESENT(hyai)) THEN
     IF ((.NOT. QDEF_NCVAR(grid%hyai) .AND. QDEF_NCVAR(hyai)) &
         .OR. (QDEF_NCVAR(grid%hyai) .AND. .NOT. QDEF_NCVAR(hyai))) RETURN
     IF (QDEF_NCVAR(grid%hyai) .AND. QDEF_NCVAR(hyai)) THEN
        IF (.NOT. QCMP_NCVAR(grid%hyai, hyai)) RETURN
     ENDIF
  ELSE
     IF (QDEF_NCVAR(grid%hyai)) RETURN
  ENDIF

  IF (PRESENT(hybi)) THEN
     IF ((.NOT. QDEF_NCVAR(grid%hybi) .AND. QDEF_NCVAR(hybi)) &
         .OR. (QDEF_NCVAR(grid%hybi) .AND. .NOT. QDEF_NCVAR(hybi))) RETURN
     IF (QDEF_NCVAR(grid%hybi) .AND. QDEF_NCVAR(hybi)) THEN
        IF (.NOT. QCMP_NCVAR(grid%hybi, hybi)) RETURN
     ENDIF
  ELSE
     IF (QDEF_NCVAR(grid%hybi)) RETURN
  ENDIF

  IF (PRESENT(timei)) THEN
     IF ((.NOT. QDEF_NCVAR(grid%timei) .AND. QDEF_NCVAR(timei)) &
         .OR. (QDEF_NCVAR(grid%timei) .AND. .NOT. QDEF_NCVAR(timei))) RETURN
     IF (QDEF_NCVAR(grid%timei) .AND. QDEF_NCVAR(timei)) THEN
        IF (.NOT. QCMP_NCVAR(grid%timei, timei)) RETURN
     ENDIF
  ELSE
     IF (QDEF_NCVAR(grid%timei)) RETURN
  ENDIF

  ! curvi-linear
  IF (PRESENT(clonm)) THEN
     IF ((.NOT. QDEF_NCVAR(grid%clonm) .AND. QDEF_NCVAR(clonm)) &
         .OR. (QDEF_NCVAR(grid%clonm) .AND. .NOT. QDEF_NCVAR(clonm))) RETURN
     IF (QDEF_NCVAR(grid%clonm) .AND. QDEF_NCVAR(clonm)) THEN
        IF (.NOT. QCMP_NCVAR(grid%clonm, clonm)) RETURN
     ENDIF
  ELSE
     IF (QDEF_NCVAR(grid%clonm)) RETURN
  ENDIF

  IF (PRESENT(clatm)) THEN
     IF ((.NOT. QDEF_NCVAR(grid%clatm) .AND. QDEF_NCVAR(clatm)) &
         .OR. (QDEF_NCVAR(grid%clatm) .AND. .NOT. QDEF_NCVAR(clatm))) RETURN
     IF (QDEF_NCVAR(grid%clatm) .AND. QDEF_NCVAR(clatm)) THEN
        IF (.NOT. QCMP_NCVAR(grid%clatm, clatm)) RETURN
     ENDIF
  ELSE
     IF (QDEF_NCVAR(grid%clatm)) RETURN
  ENDIF

  IF (PRESENT(cloni)) THEN
     IF ((.NOT. QDEF_NCVAR(grid%cloni) .AND. QDEF_NCVAR(cloni)) &
         .OR. (QDEF_NCVAR(grid%cloni) .AND. .NOT. QDEF_NCVAR(cloni))) RETURN
     IF (QDEF_NCVAR(grid%cloni) .AND. QDEF_NCVAR(cloni)) THEN
        IF (.NOT. QCMP_NCVAR(grid%cloni, cloni)) RETURN
     ENDIF
  ELSE
     IF (QDEF_NCVAR(grid%cloni)) RETURN
  ENDIF

  IF (PRESENT(clati)) THEN
     IF ((.NOT. QDEF_NCVAR(grid%clati) .AND. QDEF_NCVAR(clati)) &
         .OR. (QDEF_NCVAR(grid%clati) .AND. .NOT. QDEF_NCVAR(clati))) RETURN
     IF (QDEF_NCVAR(grid%clati) .AND. QDEF_NCVAR(clati)) THEN
        IF (.NOT. QCMP_NCVAR(grid%clati, clati)) RETURN
     ENDIF
  ELSE
     IF (QDEF_NCVAR(grid%clati)) RETURN
  ENDIF
  ! rotated curvi-linear
  IF (PRESENT(rlonm)) THEN

     IF ((.NOT. QDEF_NCVAR(grid%rlonm) .AND. QDEF_NCVAR(rlonm)) &
         .OR. (QDEF_NCVAR(grid%rlonm) .AND. .NOT. QDEF_NCVAR(rlonm))) RETURN
     IF (QDEF_NCVAR(grid%rlonm) .AND. QDEF_NCVAR(rlonm)) THEN
        IF (.NOT. QCMP_NCVAR(grid%rlonm, rlonm)) RETURN
     ENDIF
  ELSE
     IF (QDEF_NCVAR(grid%rlonm)) RETURN
  ENDIF

  IF (PRESENT(rlatm)) THEN
     IF ((.NOT. QDEF_NCVAR(grid%rlatm) .AND. QDEF_NCVAR(rlatm)) &
         .OR. (QDEF_NCVAR(grid%rlatm) .AND. .NOT. QDEF_NCVAR(rlatm))) RETURN
     IF (QDEF_NCVAR(grid%rlatm) .AND. QDEF_NCVAR(rlatm)) THEN
        IF (.NOT. QCMP_NCVAR(grid%rlatm, rlatm)) RETURN
     ENDIF
  ELSE
     IF (QDEF_NCVAR(grid%rlatm)) RETURN
  ENDIF

  IF (PRESENT(rloni)) THEN
     IF ((.NOT. QDEF_NCVAR(grid%rloni) .AND. QDEF_NCVAR(rloni)) &
         .OR. (QDEF_NCVAR(grid%rloni) .AND. .NOT. QDEF_NCVAR(rloni))) RETURN
     IF (QDEF_NCVAR(grid%rloni) .AND. QDEF_NCVAR(rloni)) THEN
        IF (.NOT. QCMP_NCVAR(grid%rloni, rloni)) RETURN
     ENDIF
  ELSE
     IF (QDEF_NCVAR(grid%rloni)) RETURN
  ENDIF

  IF (PRESENT(rlati)) THEN
     IF ((.NOT. QDEF_NCVAR(grid%rlati) .AND. QDEF_NCVAR(rlati)) &
         .OR. (QDEF_NCVAR(grid%rlati) .AND. .NOT. QDEF_NCVAR(rlati))) RETURN
     IF (QDEF_NCVAR(grid%rlati) .AND. QDEF_NCVAR(rlati)) THEN
        IF (.NOT. QCMP_NCVAR(grid%rlati, rlati)) RETURN
     ENDIF
  ELSE
     IF (QDEF_NCVAR(grid%rlati)) RETURN
  ENDIF

  ! PRESSURE
  IF (PRESENT(ps)) THEN
     IF ((.NOT. QDEF_NCVAR(grid%ps) .AND. QDEF_NCVAR(ps)) &
         .OR. (QDEF_NCVAR(grid%ps) .AND. .NOT. QDEF_NCVAR(ps))) RETURN
     IF (QDEF_NCVAR(grid%ps) .AND. QDEF_NCVAR(ps)) THEN
        IF (.NOT. QCMP_NCVAR(grid%ps, ps)) RETURN
     ENDIF
  ELSE
     IF (QDEF_NCVAR(grid%ps)) RETURN
  ENDIF

  IF (PRESENT(p0)) THEN
     IF ((.NOT. QDEF_NCVAR(grid%p0) .AND. QDEF_NCVAR(p0)) &
         .OR. (QDEF_NCVAR(grid%p0) .AND. .NOT. QDEF_NCVAR(p0))) RETURN
     IF (QDEF_NCVAR(grid%p0) .AND. QDEF_NCVAR(p0)) THEN
        IF (.NOT. QCMP_NCVAR(grid%p0, p0)) RETURN
     ENDIF
  ELSE
     IF (QDEF_NCVAR(grid%p0)) RETURN
  ENDIF

  IF (PRESENT(pollon)) THEN
     IF ((.NOT. QDEF_NCVAR(grid%pollon) .AND. QDEF_NCVAR(pollon)) &
         .OR. (QDEF_NCVAR(grid%pollon) .AND. .NOT. QDEF_NCVAR(pollon))) RETURN
     IF (QDEF_NCVAR(grid%pollon) .AND. QDEF_NCVAR(pollon)) THEN
        IF (.NOT. QCMP_NCVAR(grid%pollon, pollon)) RETURN
     ENDIF
  ELSE
     IF (QDEF_NCVAR(grid%pollon)) RETURN
 ENDIF
  IF (PRESENT(pollat)) THEN
     IF ((.NOT. QDEF_NCVAR(grid%pollat) .AND. QDEF_NCVAR(pollat)) &
         .OR. (QDEF_NCVAR(grid%pollat) .AND. .NOT. QDEF_NCVAR(pollat))) RETURN
     IF (QDEF_NCVAR(grid%pollat) .AND. QDEF_NCVAR(pollat)) THEN
        IF (.NOT. QCMP_NCVAR(grid%pollat, pollat)) RETURN
     ENDIF
  ELSE
     IF (QDEF_NCVAR(grid%pollat)) RETURN
  ENDIF
  IF (PRESENT(polgam)) THEN
     IF ((.NOT. QDEF_NCVAR(grid%polgam) .AND. QDEF_NCVAR(polgam)) &
         .OR. (QDEF_NCVAR(grid%polgam) .AND. .NOT. QDEF_NCVAR(polgam))) RETURN
     IF (QDEF_NCVAR(grid%polgam) .AND. QDEF_NCVAR(polgam)) THEN
        IF (.NOT. QCMP_NCVAR(grid%polgam, polgam)) RETURN
     ENDIF
  ELSE
     IF (QDEF_NCVAR(grid%polgam)) RETURN
  ENDIF

  ! unstructured
  IF (PRESENT(ulonm)) THEN
     IF ((.NOT. QDEF_NCVAR(grid%ulonm) .AND. QDEF_NCVAR(ulonm)) &
         .OR. (QDEF_NCVAR(grid%ulonm) .AND. .NOT. QDEF_NCVAR(ulonm))) RETURN
     IF (QDEF_NCVAR(grid%ulonm) .AND. QDEF_NCVAR(ulonm)) THEN
        IF (.NOT. QCMP_NCVAR(grid%ulonm, ulonm)) RETURN
     ENDIF
  ELSE
     IF (QDEF_NCVAR(grid%ulonm)) RETURN
  ENDIF

  IF (PRESENT(ulatm)) THEN
     IF ((.NOT. QDEF_NCVAR(grid%ulatm) .AND. QDEF_NCVAR(ulatm)) &
         .OR. (QDEF_NCVAR(grid%ulatm) .AND. .NOT. QDEF_NCVAR(ulatm))) RETURN
     IF (QDEF_NCVAR(grid%ulatm) .AND. QDEF_NCVAR(ulatm)) THEN
        IF (.NOT. QCMP_NCVAR(grid%ulatm, ulatm)) RETURN
     ENDIF
  ELSE
     IF (QDEF_NCVAR(grid%ulatm)) RETURN
  ENDIF

  IF (PRESENT(uloni)) THEN
     IF ((.NOT. QDEF_NCVAR(grid%uloni) .AND. QDEF_NCVAR(uloni)) &
         .OR. (QDEF_NCVAR(grid%uloni) .AND. .NOT. QDEF_NCVAR(uloni))) RETURN
     IF (QDEF_NCVAR(grid%uloni) .AND. QDEF_NCVAR(uloni)) THEN
        IF (.NOT. QCMP_NCVAR(grid%uloni, uloni)) RETURN
     ENDIF
  ELSE
     IF (QDEF_NCVAR(grid%uloni)) RETURN
  ENDIF

  IF (PRESENT(ulati)) THEN
     IF ((.NOT. QDEF_NCVAR(grid%ulati) .AND. QDEF_NCVAR(ulati)) &
         .OR. (QDEF_NCVAR(grid%ulati) .AND. .NOT. QDEF_NCVAR(ulati))) RETURN
     IF (QDEF_NCVAR(grid%ulati) .AND. QDEF_NCVAR(ulati)) THEN
        IF (.NOT. QCMP_NCVAR(grid%ulati, ulati)) RETURN
     ENDIF
  ELSE
     IF (QDEF_NCVAR(grid%ulati)) RETURN
  ENDIF

  IF (PRESENT(lonc)) THEN 
     IF (grid%lonc .AND. .NOT. lonc) RETURN
     IF (.NOT. grid%lonc .AND. lonc) RETURN
  ELSE
     IF (.NOT. grid%lonc) RETURN
  END IF
  IF (PRESENT(clonc)) THEN 
     IF (grid%clonc .AND. .NOT. clonc) RETURN
     IF (.NOT. grid%clonc .AND. clonc) RETURN
  ELSE
     IF (.NOT. grid%clonc) RETURN
  END IF
  IF (PRESENT(rlonc)) THEN 
     IF (grid%rlonc .AND. .NOT. rlonc) RETURN
     IF (.NOT. grid%rlonc .AND. rlonc) RETURN
  ELSE
     IF (.NOT. grid%rlonc) RETURN
  END IF
  IF (PRESENT(imask)) THEN 
     IF ((.NOT. QDEF_NCVAR(grid%imask) .AND. QDEF_NCVAR(imask)) &
         .OR. (QDEF_NCVAR(grid%imask) .AND. .NOT. QDEF_NCVAR(imask))) RETURN
     IF (QDEF_NCVAR(grid%imask) .AND. QDEF_NCVAR(imask)) THEN
        IF (.NOT. QCMP_NCVAR(grid%imask, imask)) RETURN
     ENDIF
  ELSE
     IF (QDEF_NCVAR(grid%imask)) RETURN
  END IF
  ! pressi3d is intentionally not compared (see above)
  ! FOR IMPORTED GRID
  IF (PRESENT(file)) THEN
     IF (TRIM(grid%file) /= TRIM(file)) RETURN
  ELSE
     IF (TRIM(grid%file) /= '' ) RETURN 
  ENDIF

  IF (PRESENT(t)) THEN
     IF (grid%t /= t) RETURN
  ELSE
     IF (grid%t /= 1) RETURN 
  ENDIF

  IF (PRESENT(corners)) THEN
     IF (grid%corners /= corners) RETURN
  ELSE
     IF (grid%corners /= 4) RETURN 
  ENDIF
 
  IF (PRESENT(ranges)) THEN
     DO i = 1, 4
        DO j= 1,2
           IF (grid%ranges(i,j) /=  ranges(i,j)) RETURN
        END DO
     END DO
  ELSE
     DO i = 1, 4
        DO j= 1,2
           IF (grid%ranges(i,j) /=  RGEMPTY) RETURN
        END DO
     END DO
  ENDIF
  IF (PRESENT(minmaxlonlat)) THEN
     DO i = 1, 2
        DO j= 1,2
           IF (grid%minmaxlonlat(i,j) /=  minmaxlonlat(i,j)) RETURN
        END DO
     END DO
  ELSE
     DO i = 1, 2
        DO j= 1,2
           IF (grid%minmaxlonlat(i,j) /=  RGEMPTY) RETURN
        END DO
     END DO
  ENDIF

  ID = grid%id

END SUBROUTINE COMPARE_TO_GRID
! ===========================================================================

SUBROUTINE COMPARE_TWO_GRIDS(agrid, bgrid, status)

  IMPLICIT NONE

  ! I/O
  TYPE(t_geohybgrid), POINTER     :: agrid
  TYPE(t_geohybgrid), POINTER     :: bgrid
  INTEGER,            INTENT(OUT) :: status
  ! LOCAL
  INTEGER                         :: ID

  status = 0
  ID = -99

  CALL COMPARE_TO_GRID(ID, agrid &
       , bgrid%lonm,  bgrid%latm , bgrid%hyam,   bgrid%hybm,   bgrid%timem  &
       , bgrid%loni,  bgrid%lati,  bgrid%hyai,   bgrid%hybi,   bgrid%timei  &
       , bgrid%clonm, bgrid%clatm, bgrid%cloni,  bgrid%clati                & 
       , bgrid%rlonm, bgrid%rlatm, bgrid%rloni,  bgrid%rlati                & 
       , bgrid%ps,    bgrid%p0,    bgrid%file,   bgrid%t,      bgrid%ranges &
       , bgrid%minmaxlonlat,       bgrid%pollon, bgrid%pollat, bgrid%polgam &
       , bgrid%ulonm, bgrid%ulatm, bgrid%uloni,  bgrid%ulati                &
       , bgrid%corners, bgrid%lonc, bgrid%clonc, bgrid%rlonc, bgrid%imask )

  IF (ID == -99) status = 2020 ! 'GRIDS ARE NOT EQUAL'

END SUBROUTINE COMPARE_TWO_GRIDS

! ====================================================================
SUBROUTINE GEOHYBGRID_UPDATE_PRESS(status, ID, pressi, pressm)

  ! UPDATE pressi field of grid
  ! 1. search grid according to the given ID
  ! 2. copy ncvar pressure to grid

  IMPLICIT NONE

  ! I/O
  INTEGER,            INTENT(OUT)            :: status
  INTEGER,            INTENT(IN)             :: ID
  TYPE(t_ncvar),      INTENT(IN), OPTIONAL   :: pressi
  TYPE(t_ncvar),      INTENT(IN), OPTIONAL   :: pressm

  ! LOCAL
  TYPE(t_geohybgrid_list), POINTER :: gi => NULL()
  TYPE(t_geohybgrid_list), POINTER :: ge => NULL()
  LOGICAL                          :: lexists
  !CHARACTER(LEN=*), PARAMETER     :: substr = 'GEOHYBGRID_UPDATE_PRESS'

  lexists = .FALSE.
  status  = 1

  ! CHECKS
  IF (id <= 0) THEN
     status = 2007  ! INVALID GRID ID
     RETURN
  END IF
  !
  IF (.NOT. ASSOCIATED(GEOHYBGRIDLIST)) THEN
     status = 2004  ! GRID (ID) DOES NOT EXIST
     RETURN
  END IF

  gi => GEOHYBGRIDLIST
  DO
     IF (.NOT. ASSOCIATED(gi)) EXIT
     IF (id == gi%this%id) THEN
        lexists = .TRUE.
        EXIT
     END IF
     ge => gi
     gi => ge%next
  END DO

  IF (lexists) THEN
     IF (PRESENT(pressi)) CALL COPY_NCVAR(gi%this%pressi,pressi)
     IF (PRESENT(pressm)) CALL COPY_NCVAR(gi%this%pressm,pressm)
  ELSE
     status = 2003  ! GRID DOES NOT EXIST
     RETURN
  END IF
  NULLIFY(gi)
  NULLIFY(ge)

  status = 0

END SUBROUTINE GEOHYBGRID_UPDATE_PRESS

! ====================================================================
SUBROUTINE GEOHYBGRID_UPDATE_SURFPRESS(status, ID, ps)

  ! UPDATE surface pressure of grid
  ! 1. search grid according to the given ID
  ! 2. copy ncvar pressure to grid

  IMPLICIT NONE

  ! I/O
  INTEGER,            INTENT(OUT)  :: status
  INTEGER,            INTENT(IN)   :: ID
  TYPE(t_ncvar),      INTENT(IN)   :: ps


  ! LOCAL
  TYPE(t_geohybgrid_list), POINTER :: gi => NULL()
  TYPE(t_geohybgrid_list), POINTER :: ge => NULL()
  LOGICAL                          :: lexists
  !CHARACTER(LEN=*), PARAMETER     :: substr = 'GEOHYBGRID_UPDATE_SURFPRESS'

  lexists = .FALSE.
  status  = 1

  ! CHECKS
  IF (id <= 0) THEN
     status = 2007  ! INVALID GRID ID
     RETURN
  END IF
  !
  IF (.NOT. ASSOCIATED(GEOHYBGRIDLIST)) THEN
     status = 2004  ! GRID (ID) DOES NOT EXIST
     RETURN
  END IF

  gi => GEOHYBGRIDLIST
  DO
     IF (.NOT. ASSOCIATED(gi)) EXIT
     IF (id == gi%this%id) THEN
        lexists = .TRUE.
        EXIT
     END IF
     ge => gi
     gi => ge%next
  END DO

  IF (lexists) THEN
     CALL COPY_NCVAR(gi%this%ps,ps)
  ELSE
     status = 2003  ! GRID DOES NOT EXIST
     RETURN
  END IF
  NULLIFY(gi)
  NULLIFY(ge)

  status = 0


END SUBROUTINE GEOHYBGRID_UPDATE_SURFPRESS

! ====================================================================
  SUBROUTINE CLEAN_GEOHYBGRID_LIST(status)

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,  INTENT(OUT) :: status

    ! LOCAL
    TYPE(t_geohybgrid_list),      POINTER :: psgi, psge
    TYPE(t_geohybgrid),           POINTER :: grid

    status = 0

    psgi => GEOHYBGRIDLIST
    channel_loop: DO
       IF (.NOT. ASSOCIATED(psgi)) EXIT

       grid => psgi%this

       psge => psgi
       psgi => psgi%next

       ! -------------------------------------------
       ! MEMORY: DATA
       CALL INIT_GEOHYBGRID(grid)
       DEALLOCATE(psge)
       NULLIFY(psge)
       !
       ! COUNT
       NGEOHYBGRID = NGEOHYBGRID - 1
       ! -------------------------------------------

    END DO channel_loop

    NULLIFY(GEOHYBGRIDLIST)

    IF (NGEOHYBGRID /= 0) THEN
       status = 3045 ! INTERNAL SCRIP GRID COUNT ERROR
       RETURN
    END IF

    status = 0

  END SUBROUTINE CLEAN_GEOHYBGRID_LIST
! ====================================================================

! ===========================================================================
CHARACTER(LEN=256) FUNCTION GRID_ERROR(status)

  ! This subroutine provides an error string for a given status

  IMPLICIT NONE
     
  INTEGER, INTENT(IN) :: status
    
  SELECT CASE(status)
     !
     ! GRID ERRORS
        !
  CASE(2001)
     grid_error = TRIM('GRID NAME TOO LONG')
  CASE(2002)
     grid_error = TRIM('GRID EXISTS ALREADY')
  CASE(2003)
     grid_error = TRIM('GRID DOES NOT EXIST')
  CASE(2004)
     grid_error = TRIM('GRID (ID) DOES NOT EXIST')
  CASE(2007)
     grid_error = TRIM('INVALID GRID ID')
  CASE(2008)
     grid_error = TRIM('GRID EXITS AND IS EQUAL')
  CASE(2009)
     grid_error = TRIM('GRID NAME EXISTS ALREADY')
  CASE(2010)
     grid_error = TRIM('GRID NAME IS EMPTY')
  CASE(2020)
     grid_error = TRIM('GRIDS ARE NOT EQUAL')
  CASE(2050)
     grid_error = TRIM('ERROR IN NEW_GEOHYBGRID')
  CASE(2100)
     grid_error = TRIM('ps grid dimensions do not match')
     !
     ! SCRIP
     ! 
  CASE(3001)
     grid_error = TRIM('SCRIP: ARRAY DIMENSIONS FOR LON and LAT not confrom')
  CASE(3002)
     grid_error = TRIM('SCRIP: DIMENSION SIZES FOR LON and LAT not confrom')
  CASE(3003)
     grid_error = TRIM('SCRIP: NEITHER vr nor vd are associated')
  CASE(3004)
     grid_error = TRIM('SCRIP: ARRAY DIMENSIONS FOR REGULAR GRID SHOULD BE 1')
  CASE(3005)
     grid_error = TRIM('SCRIP: NEITHER MIDPOINTs NOR BOX INTERFACES ARE DEFINED')
  CASE(3006)
     grid_error = &
          TRIM('SCRIP: ARRAY DIMENSIONS FOR CURVILINEAR GRID SHOULD BE 2')
  CASE(3007)
     grid_error = &
          TRIM('SCRIP: SCRIP can only handle equal regrid types per data set')
  CASE(3008)
     grid_error = TRIM('SCRIP: INDEX REGRIDDING NOT YET POSSIBLE IN SCRIP')
  CASE(3009)
     grid_error = TRIM('SCRIP: DATA (ID) DOES NOT EXIST')
  CASE(3010)
     grid_error = TRIM('SCRIP: COPY OF SCRIP DATA NOT IMPLEMENTED')
  CASE(3011)
     grid_error = TRIM('SCRIP: INDEX FRACTIONS NOT YET IMPLEMENTED')
  CASE(3012)
     grid_error = TRIM('SCRIP: map_type not known')
  CASE(3013)
     grid_error = TRIM('SCRIP: DATA DOES NOT EXIST')
  CASE(3014)
     grid_error = TRIM('SCRIP: ! VTYPE NOT IMPLEMENTED')
  CASE(3015)
     grid_error = TRIM('SCRIP: INVALID DATA ID')
  CASE(3016)
     grid_error = TRIM('SCRIP:  BALANCE_GEOHYBGRID_PS, PSD NOT PRESENT')
  CASE(3017)
     grid_error = TRIM('SCRIP: ')
  CASE(3018)
     grid_error = TRIM('SCRIP: longitude range too small < -180.')
  CASE(3019)
     grid_error = TRIM('SCRIP: longitude range too large > 360.')
  CASE(3020)
     grid_error = TRIM('SCRIP: longitude range too wide > 360.')
  CASE(3021)
     grid_error = TRIM('SCRIP: for curvi-linear grids interfaces are required for corner calculation ')
  CASE(3022)
     grid_error = TRIM('SCRIP: SIZES INCOMPATIBLE')
  CASE(3023)
     grid_error = TRIM('SCRIP: for curvi-linear grids mids are required for center calculation ')
  CASE(3024) 
     grid_error = TRIM('SCRIP: set_ranges: neither lloni nor llonm associated')
  CASE(3025) 
     grid_error = TRIM('SCRIP: mask and lon / lat arrays differ in longitude')
  CASE(3030)
     grid_error = TRIM('SCRIP: NORMALIZE OPTION NOT IMPLEMENTED')

  CASE(3040)
     grid_error = TRIM('SCRIP: INTERNAL SCRIP DATA COUNT ERROR')
  CASE(3045)
     grid_error = TRIM('SCRIP: INTERNAL SCRIP GRID COUNT ERROR')

  CASE(3090)
     grid_error = TRIM(&
          'SCRIP: curvilinear and only interfaces defined not yet implemented')
  CASE(3091)
     grid_error = TRIM(&
          'SCRIP: curvilinear and no interfaces defined not yet implemented')
  CASE(3998)
     grid_error = TRIM('SCRIP: integration stalled: num_subseg exceeded limit')
  CASE(3999)
     grid_error = TRIM('SCRIP: UNKOWN ERROR')
     ! 
  CASE(4000)
     grid_error = TRIM('TOOLS: DIMENSIONS OF VAR MUST BE PERDEFINED')
     ! FORCED EXIT FOR TESTS
  CASE(5000) 
     grid_error = TRIM('FORCED EXIT FOR TESTS')
  CASE(6000)
     grid_error = TRIM('IMPORT_GRID: namelist scan: no suitable namelist found')

  END SELECT
     
END FUNCTION GRID_ERROR
 ! ===========================================================================

! ******************************************************************
END MODULE MESSY_MAIN_GRID
! ******************************************************************
