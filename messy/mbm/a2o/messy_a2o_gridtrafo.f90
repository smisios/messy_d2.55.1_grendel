!********************************************************************
!--------------------------------------------------------------------
MODULE MESSY_A2O_GRIDTRAFO
  !------------------------------------------------------------------
  ! Author: Bastian Kern, MPICH, Mainz, July 2009 / September 2010
  ! Great parts of this code are from ncregrid
  ! (messy_ncregrid_control.f90) by Patrick Joeckel, MPICH, Mainz,
  ! with modifications.
  ! Parts are from messy_a2o_e5.f90 by Andrea Pozzer, MPICH, Mainz.
  ! Parts are from the mpi-om, MPI-M, Hamburg.
  !******************************************************************

  USE MESSY_NCREGRID_NETCDF, ONLY : NCVAR, NCATT, NARRAY                      &
       , GRD_MAXSTRLEN                                                        &
       , NF90_inq_libvers, NF90_GLOBAL, NF90_CHAR                             &
       , NULL_DIMID                                                           &
       , INIT_NCVAR, RENAME_NCVAR, COPY_NCVAR                                 &
       , IMPORT_NCVAR, SCAN_NCVAR, EXPORT_NCVAR                               &
       , INIT_NCATT, IMPORT_NCATT, ADD_NCATT, EXPORT_NCATT                    &
       , INIT_NARRAY, CAT_NARRAY, COPY_NARRAY                                 &
       , COPY_NCDIM, QCMP_NCDIM

  USE MESSY_NCREGRID_GEOHYB, ONLY : GORD_LON, GORD_LAT, GORD_LEV              &
       , GEOHYBGRID                                                           &
       , INIT_GEOHYBGRID, COPY_GEOHYBGRID, IMPORT_GEOHYBGRID                  &
       , CHECK_NCVAR_ON_GEOHYBGRID, EXPORT_GEOHYBGRID                         &
       , BALANCE_GEOHYBGRID_NCVAR, BALANCE_GEOHYBGRID                         &
       , BALANCE_GEOHYBGRID_TIME, PACK_GEOHYBGRID_NCVAR
       
  USE MESSY_MAIN_TOOLS, ONLY : READ_NML_OPEN, READ_NML_CHECK                  &
       , READ_NML_CLOSE, FIND_NEXT_FREE_UNIT

! um_ak_20130524 file name changed
!  USE MESSY_MAIN_GRIDTRAFO_SCRIP, ONLY:                                       &
  USE MESSY_MAIN_GRID_TRAFO_SCRP_BASE, ONLY:                                  &
         grid1_center_lat, grid1_center_lon                                   &
       , grid1_corner_lat, grid1_corner_lon                                   &
       , grid1_mask, grid1_area_in                                            &
       , grid2_center_lat, grid2_center_lon                                   &
       , grid2_corner_lat, grid2_corner_lon                                   &
       , grid2_mask, grid2_area_in                                            &
       , deg2rad                                                              &
       , ATM, OCES, OCEU, OCEV                                                &
       , NORM_OPT_NONE, NORM_OPT_DSTAREA, NORM_OPT_FRCAREA                    &
       , MAP_TYPE_CONSERV, MAP_TYPE_BILINEAR                                  &
       , MAP_TYPE_BICUBIC, MAP_TYPE_DISTWGT                                   &
       , MAP_TYPE, NORM_OPT                                                   &
       , num_links_map1, num_wts                                              &
       , grid1_add_map1, grid1_area                                           &
       , grid2_add_map1, grid2_area                                           &
       , LUSE_GRID_CENTERS                                                    &
       , wts_map1, grid2_frac                                                 &
       , PI2, PIH, ZERO, TINY                                                 &
       , REMAP_INIT, BOUNDS_CALC, REMAP_VARS, REMAP_CONSERV, REMAP_DEALLOC

  USE MESSY_MAIN_CONSTANTS_MEM, ONLY : RADIUS_EARTH, PI, DP, SP

  USE MESSY_NCREGRID_DIAG, ONLY : WRITE_NCVAR, WRITE_NARRAY

  USE MESSY_NCREGRID_BASE, ONLY : ERRMSG, RGMSG                               &
       , RGMLE, RGMLEC, RGMLVL, RGMLVLC, RGMLVM, RGMLVMC                      &
       , RGMLI, RGMLIC, RGMLW, RGMLWC                                         &
       , VTYPE_REAL, VTYPE_DOUBLE, VTYPE_CHAR, VTYPE_UNDEF                    &
       , POSITION, ELEMENT, DOUBLE_NARRAY

  IMPLICIT NONE

  INTRINSIC :: TRIM, DATE_AND_TIME, ANY, ASSOCIATED, PRESENT, SIZE, ABS, ACOS &
       , COS, SIN, MIN, MAX, NULL, PACK, UNPACK, SQRT, ACHAR, LEN_TRIM        &
       , ALLOCATED, INDEX, PRODUCT, MAXVAL, MINVAL

  PRIVATE   :: TRIM, DATE_AND_TIME, ANY, ASSOCIATED, PRESENT, SIZE, ABS, ACOS &
       , COS, SIN, MIN, MAX, NULL, PACK, UNPACK, SQRT, ACHAR, LEN_TRIM        &
       , ALLOCATED, INDEX, PRODUCT, MAXVAL, MINVAL

  CHARACTER(*), PARAMETER :: A2OVERS = 'v0.5 21 September 2010'

  ! control parameter
  INTEGER, SAVE           :: GT_CTRL, GT_NML, GT_STATUS

  INTEGER, PARAMETER      :: GT_SCAN = 1        ! control
  INTEGER, PARAMETER      :: GT_PROC = 2
  INTEGER, PARAMETER      :: GT_STOP = 3

  INTEGER, PARAMETER      :: NML_NEXT = 1       ! namelist
  INTEGER, PARAMETER      :: NML_STAY = 2

  INTEGER, PARAMETER      :: GTSTAT_START = 1   ! status
  INTEGER, PARAMETER      :: GTSTAT_CONT = 2
  INTEGER, PARAMETER      :: GTSTAT_STOP = 3

    TYPE INT_INTERP
       INTEGER :: source_grid
       INTEGER :: dest_grid
       INTEGER :: map_method
       INTEGER :: norm_option
    END TYPE INT_INTERP

    ! Calculated fields for the grid transformation
    TYPE GRIDTRAFO
       INTEGER                           :: nn
       INTEGER,  DIMENSION(:),   POINTER :: srcadd => NULL()   ! src addresses
       INTEGER,  DIMENSION(:),   POINTER :: destadd => NULL()  ! des addresses
       REAL(DP), DIMENSION(:),   POINTER :: srcarea => NULL()  ! src area
       REAL(DP), DIMENSION(:),   POINTER :: destarea => NULL() ! des area
       REAL(DP), DIMENSION(:,:), POINTER :: weight => NULL()   ! trf. weigths
       REAL(DP), DIMENSION(:),   POINTER :: frac => NULL()     ! grid fraction
    END TYPE GRIDTRAFO

    ! Grid transformation
    TYPE TRAFO
       TYPE (INT_INTERP) :: info
       TYPE (GRIDTRAFO), POINTER  :: transform
    END TYPE TRAFO

    TYPE CORNERARRAY
       REAL(DP), DIMENSION(:,:), POINTER :: corners
    END TYPE CORNERARRAY

    TYPE MASKARRAY
       LOGICAL, DIMENSION(:), POINTER :: mask
    END TYPE MASKARRAY

    ! if you change anything here, make sure, all the strings in the array    &
    ! have the same length (24)
    CHARACTER(LEN=24), DIMENSION(4), PARAMETER :: gridstring                  &
         = (/"ATMOSPHERE              "                                       &
         , "OCEAN SCALAR            "                                         &
         , "OCEAN U-VECTOR COMPONENT"                                         &
         , "OCEAN V-VECTOR COMPONENT"/)


CONTAINS

  !------------------------------------------------------------------
  SUBROUTINE GRIDTRAFO_CONTROL(TCTRL, TNML, TSTAT                             &
       , nmlfile                                                              &
       , iounit                                                               &
       , var                                                                  &
       , grid                                                                 &
       )

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN)          :: TCTRL   ! CTRL: GT_SCAN, GT_PROC, GT_STOP
    INTEGER, INTENT(IN)          :: TNML    ! NML CONTROL: NML_NEXT, NML_STAY
    INTEGER, INTENT(OUT)         :: TSTAT   ! STATUS: GTSTAT_START,
                                            ! GTSTAT_CONT, GTSTAT_STOP
    CHARACTER(LEN=*), INTENT(IN) :: nmlfile ! namelist file

    ! Optional I/O
    INTEGER, INTENT(IN), OPTIONAL                 :: iounit  ! i/o unit nmlfile
    TYPE (NCVAR), DIMENSION(:), POINTER, OPTIONAL :: var     ! variable list
    TYPE (GEOHYBGRID), INTENT(OUT), OPTIONAL      :: grid    ! out grid info

    ! Local
    CHARACTER(LEN=*), PARAMETER   :: substr = 'GRIDTRAFO_CONTROL'

    INTEGER :: iou        ! i/o unit namelist file
    INTEGER :: iout       ! i/o unit
    INTEGER :: i,j,k,l    ! counter variables
    INTEGER :: ii,jj      ! more counters
    INTEGER :: nvars      ! number of variables
    INTEGER :: nvecs      ! number of vectors
    INTEGER :: mvars      ! effective number of variables
    INTEGER :: status     ! error return status

    INTEGER :: setuid     ! unlimited id

    LOGICAL :: lex        ! file exists?
    LOGICAL :: lopn       ! file open?

    LOGICAL :: ldegrees

    INTEGER :: pos, posjm, posjp, poskm, poskp   ! position k/j:+/-

    TYPE (NCVAR), DIMENSION(:), POINTER :: xivar ! list of variables (input)
    TYPE (NCVAR), DIMENSION(:), POINTER :: xitvar! temporary variable
    TYPE (NCVAR), DIMENSION(:), POINTER :: xovar ! list of variables (output)
    TYPE (NCVAR)                        :: tvar  ! temporary variable

    ! Namelist/grids
    INTEGER, SAVE                    :: fstat       ! status of nmlfile
    TYPE (GEOHYBGRID), SAVE          :: gi0         ! raw input grid
    TYPE (GEOHYBGRID), SAVE          :: ggatm0      ! raw atm grdfile grid
    TYPE (GEOHYBGRID), SAVE          :: ggoces0     ! raw oces grdfile grid
    TYPE (GEOHYBGRID), SAVE          :: ggoceu0     ! raw oceu grdfile grid
    TYPE (GEOHYBGRID), SAVE          :: ggocev0     ! raw ocev grdfile grid
    TYPE (GEOHYBGRID)                :: gi     ! input grid
    TYPE (GEOHYBGRID)                :: go     ! output grid

    TYPE (GEOHYBGRID), DIMENSION(4)  :: ggrid  ! gridfile grid

    LOGICAL :: lover,lpress                  ! output ocean grid overlapping? &
                                             ! pressure var in nml varstr?

    CHARACTER (LEN=100*GRD_MAXSTRLEN), SAVE :: varstr  ! variable string
    CHARACTER (LEN=100*GRD_MAXSTRLEN), SAVE :: vecstr  ! variable string
    CHARACTER (LEN=GRD_MAXSTRLEN), SAVE :: outfile     ! NetCDF outfile
    INTEGER, SAVE                       :: i_t(4)      ! input time control
    INTEGER, SAVE                       :: g_t(3)      ! grid time control
    INTEGER, SAVE                       :: o_t(3)      ! output time control
    LOGICAL, SAVE                       :: lint        ! use input time?
    TYPE (NCATT), SAVE                  :: gtgatt      ! global ncdf attributes
    TYPE (NARRAY), SAVE                 :: nmlstr      ! namelist as string

    ! Time control
    INTEGER                             :: ttmin, ttmax, ttstep, ttret
    INTEGER                             :: tg     ! grid time step
    INTEGER                             :: to     ! output time step
    INTEGER                             :: tc     ! time step counter

    LOGICAL :: ok

    ! Transformation
    TYPE(TRAFO), DIMENSION(4,4,4,3) :: transformation
    LOGICAL, DIMENSION(4,4,4,3)     :: INTCALCTRAFO = .false. ! which         &
                                             ! transformations have to be done?

    ! ncvars for regridding
    TYPE(NCVAR), DIMENSION(:), ALLOCATABLE   :: nai, nao         ! in and out

    ! Temporary arrays for vector rotation
    REAL(DP), DIMENSION(:,:), ALLOCATABLE         :: xvecdat
    REAL(DP), DIMENSION(:,:), ALLOCATABLE         :: yvecdat
    REAL(DP), DIMENSION(:,:), ALLOCATABLE         :: vec_tmp_x
    REAL(DP), DIMENSION(:,:), ALLOCATABLE         :: vec_tmp_y

    ! Variables
    INTEGER, POINTER                      :: SOURCEGRID(:)    ! src-grid
    INTEGER, POINTER                      :: DESTGRID(:)      ! dest-grid
    INTEGER, POINTER                      :: MAPPING(:)       ! mapping method
    INTEGER, POINTER                      :: NORMALIZING(:)   ! norm. option
    INTEGER, POINTER                      :: TEMPSOURCE(:)    ! temporary field
    INTEGER, POINTER                      :: TEMPDEST(:)      ! temporary field
    INTEGER, POINTER                      :: TEMPMAP(:)       ! temporary field
    INTEGER, POINTER                      :: TEMPNORM(:)      ! temporary field
    CHARACTER(LEN=GRD_MAXSTRLEN), POINTER :: SOURCEGRIDstr(:) ! src-grid string
    CHARACTER(LEN=GRD_MAXSTRLEN), POINTER :: DESTGRIDstr(:)   ! dst-grid string
    CHARACTER(LEN=GRD_MAXSTRLEN), POINTER :: MAPstr(:)        ! rgr.type string
    CHARACTER(LEN=GRD_MAXSTRLEN), POINTER :: NORMstr(:)       ! rgr.type string

    CHARACTER(LEN=GRD_MAXSTRLEN), POINTER :: ivar(:)    ! input variable names
    CHARACTER(LEN=GRD_MAXSTRLEN), POINTER :: ovar(:)    ! output variable names
    CHARACTER(LEN=GRD_MAXSTRLEN), POINTER :: xvar(:)    ! vec x variable names
    CHARACTER(LEN=GRD_MAXSTRLEN), POINTER :: yvar(:)    ! vec y variable names
    INTEGER, POINTER                      :: xvaridx(:) ! vec x variable idx
    INTEGER, POINTER                      :: yvaridx(:) ! vec y variable idx

    ! Temporary variables
    INTEGER                               :: destadd,srcadd ! array addresses
    REAL(DP)                              :: factor         ! norm factor
    REAL(DP)                              :: value          ! temp value

    CHARACTER(LEN=100*GRD_MAXSTRLEN)      :: attstr  ! attribute string

    INTEGER, DIMENSION(:)       , POINTER :: dims    ! order of grid dimensions
    INTEGER                               :: axes(3) ! grid axes of variable

    ! Number of x/y points
    INTEGER, DIMENSION(4)                 :: nlon, nlat
    INTEGER, DIMENSION(:), POINTER        :: in_nlon, in_lev, in_nz

    ! Number of corners (should be 4)
    ! NEVER USED
    !INTEGER                               :: nc=4 ! (do not change!)

    ! Grid fields
    TYPE(NARRAY), DIMENSION(4), TARGET           :: longitude
    TYPE(NARRAY), DIMENSION(4), TARGET           :: latitude
    TYPE(MASKARRAY), DIMENSION(4), TARGET        :: gridmask
    TYPE(NARRAY), DIMENSION(4), TARGET           :: gridcellarea
    TYPE(NCVAR), DIMENSION(2)                    :: corners
    TYPE(CORNERARRAY), DIMENSION(4), TARGET      :: corner_lon
    TYPE(CORNERARRAY), DIMENSION(4), TARGET      :: corner_lat

    INTEGER, DIMENSION(4) :: gg_nlon, gg_nlat, gg_nc

    INTEGER, DIMENSION(3) :: dim_in, dim_out
    INTEGER, DIMENSION(:), POINTER :: srcvec, destvec

    ! Temporary variables
    REAL(DP), DIMENSION(:,:), ALLOCATABLE      :: outlat
    REAL(DP), DIMENSION(:,:), ALLOCATABLE      :: outlon

    INTEGER :: natype ! type of data array
    INTEGER :: nidx   ! number of associated vectors

    ! vectors for calculation of position in 1-d array
    INTEGER, DIMENSION(3) :: cornervec, dimvec
    INTEGER, DIMENSION(:), ALLOCATABLE :: ovec

    ! Variables for area calculations
    REAL(DP) :: sinlam, sinlap, lonm, lonp

    ! Temporary variables for overlapping ocean grids
    REAL(DP), DIMENSION(:), ALLOCATABLE :: out
    REAL(DP), DIMENSION(:), ALLOCATABLE :: stripes
    LOGICAL,  DIMENSION(:), ALLOCATABLE :: mask

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    IF ((TNML /= NML_NEXT).AND.(TNML /= NML_STAY)) THEN
       CALL RGMSG(substr, RGMLE, 'UNRECOGNIZED VALUE OF GT_NML!')
    END IF

    ! Default values
    IF (PRESENT(iounit)) THEN
       iou = iounit
    ELSE
       iou = 17
    END IF

    ! Intitialize output parameter
    IF (PRESENT(var)) THEN
       IF (ASSOCIATED(var)) THEN
          DO i=1, SIZE(var)
             CALL INIT_NCVAR(var(i))
          END DO
          DEALLOCATE(var, STAT=status)
          CALL ERRMSG(substr, status, 1)
       END IF
       NULLIFY(var)
    END IF

    IF (PRESENT(grid)) CALL INIT_GEOHYBGRID(grid)

    ! Check namelist file; open file
    INQUIRE(file=TRIM(nmlfile), exist=lex, opened=lopn, number=iout)
    IF (.NOT.lex) THEN
       CALL RGMSG(substr, RGMLE, &
            'NAMELIST-FILE '''//TRIM(nmlfile)//''' NOT FOUND !')
    END IF
    IF (lopn) THEN  ! File open
       IF (iout /= iou) THEN   ! Check unit
          CALL RGMSG(substr, RGMLE, &
               'NAMELIST-FILE '''//TRIM(nmlfile)//''' ALREADY OPEN ON UNIT'   &
               , iout,' !')
       ELSE
          TSTAT = GTSTAT_CONT
       END IF
    ELSE            ! File not open
       CALL RGMSG(substr, RGMLVL, &
            '-------------------------------------------------------')
       CALL RGMSG(substr, RGMLVL, &
            'START REGRIDDING PROCEDURE (A2O VERSION '//&
            TRIM(A2OVERS)//')' )
       CALL RGMSG(substr, RGMLVL, &
            ' <Author: Bastian Kern, MPICH, July 2009>' )
       CALL RGMSG(substr, RGMLVL, &
            ' ( COMPILED WITH netCDF Fortran90 LIBRARY' )
       CALL RGMSG(substr, RGMLVL, &
            '   VERSION '//TRIM(NF90_inq_libvers())//' )' )
       CALL RGMSG(substr, RGMLVL, &
            '-------------------------------------------------------')
       CALL RGMSG(substr, RGMLVL, &
            'OPENING NAMELIST FILE '''//TRIM(nmlfile)//''' ON UNIT'           &
            , iou,' !')
       OPEN(iou,file=TRIM(nmlfile))
       TSTAT = GTSTAT_START
    END IF

    ! Read next namelist
    IF ((TSTAT == GTSTAT_START).OR.(TNML == NML_NEXT)) THEN
       CALL RGMSG(substr, RGMLI, &
            'READING NAMELIST FROM FILE '''//TRIM(nmlfile)//''' ...')
       CALL READ_A2O_NML(iou, fstat                                           &
            , gi0, ggatm0, ggoces0, ggoceu0, ggocev0                          &
            , outfile                                                         &
            , varstr                                                          &
            , vecstr                                                          &
            , i_t, g_t, o_t, lint, lover                                      &
            ,nmlstr)
       ! Check end of file
       IF (fstat == 0) THEN    ! Go on
          CALL RGMSG(substr, RGMLIC, '... OK !')
          TSTAT = GTSTAT_CONT
       ELSE                    ! Error or end of file
          IF (fstat > 0) THEN  ! Error
             CLOSE(iou)
             CALL RGMSG(substr, RGMLE                                         &
                  , 'READ ERROR READING NAMELIST FROM FILE '''                &
                  //TRIM(nmlfile)//''' !' , .false.)
             CALL RGMSG(substr, RGMLEC                                        &
                  , 'STATUS: ',fstat,' ')
          ELSE                 ! End of file
             CALL RGMSG(substr, RGMLVL, &
                  'END OF FILE '''//TRIM(nmlfile)//''' REACHED !' )
          END IF               ! Error or end of file
          TSTAT = GTSTAT_STOP
       END IF
    END IF

    IF (TCTRL == GT_STOP) TSTAT = GTSTAT_STOP

    !------------------------------------------------------------------
    IF (TSTAT == GTSTAT_STOP) THEN

       CLOSE(iou)

       ! Clean up
       CALL RGMSG(substr, RGMLVM, 'CLEANING MEMORY ...')
       IF (PRESENT(grid)) THEN
          CALL INIT_GEOHYBGRID(grid)
       END IF
       !
       IF (PRESENT(var)) THEN
          IF (ASSOCIATED(var)) THEN
             DO i=1, SIZE(var)
                CALL INIT_NCVAR(var(i))
             END DO
             DEALLOCATE(var, STAT=status)
             CALL ERRMSG(substr, status, 2)
          END IF
          NULLIFY(var)
       END IF
       !

       CALL INIT_GEOHYBGRID(gi0)
       CALL INIT_GEOHYBGRID(ggatm0)
       CALL INIT_GEOHYBGRID(ggoces0)
       CALL INIT_GEOHYBGRID(ggoceu0)
       CALL INIT_GEOHYBGRID(ggocev0)
       CALL INIT_GEOHYBGRID(go)
       CALL INIT_NARRAY(nmlstr)
       CALL INIT_NCATT(gtgatt)

       CALL RGMSG(substr, RGMLVM, '... DONE')
       RETURN
    END IF

    ! --- TIME CONTROL ---------------------------------
    ! update time control from namelist
    ttmin = i_t(1)
    ttstep = i_t(2)
    ttmax = i_t(3)
    ttret = i_t(4)


    ! overwrite time control from subroutine call
    !    IF (PRESENT(tmin)) THEN
    !       ttmin = tmin
    !    END IF
    !    IF (PRESENT(tmax)) THEN
    !       ttmax = tmax
    !    END IF
    !    IF (PRESENT(tstep)) THEN
    !       ttstep = tstep
    !    END IF
    !    IF (PRESENT(tret)) THEN
    !       ttret = tret
    !    END IF

    ! check time control settings
    ! check infile time stepping / interface mode time stepping

    IF (ttmin <= 0) THEN
       CALL RGMSG(substr, RGMLE,  'i_t(1) IN NAMELIST OR tmin AT ', .false.)
       CALL RGMSG(substr, RGMLEC, 'SUBROUTINECALL MUST BE >0 !')
    END IF

    IF (ttstep <= 0) THEN
       CALL RGMSG(substr, RGMLE,  'i_t(2) IN NAMELIST OR tstep AT ', .false.)
       CALL RGMSG(substr, RGMLEC, 'SUBROUTINECALL MUST BE >0 !')
    END IF

    IF (ttmax < 0) THEN
       CALL RGMSG(substr, RGMLE,  'i_t(3) IN NAMELIST OR tmax AT ', .false.)
       CALL RGMSG(substr, RGMLEC, 'SUBROUTINECALL MUST BE >=0 !')
    END IF

    IF (ttret < 0) THEN
       CALL RGMSG(substr, RGMLE,  'i_t(4) IN NAMELIST OR tret AT ', .false.)
       CALL RGMSG(substr, RGMLEC, 'SUBROUTINECALL MUST BE >=0 !')
    END IF

    IF (ttmax == 0) THEN  ! special case
       CALL RGMSG(substr, RGMLW,  'tmax (i_t(3)) RESET TO tmin (i_t(1)) !')
       ttmax = ttmin
    END IF

    IF (ttret == 0) THEN  ! special case
       CALL RGMSG(substr, RGMLW,  'tret (i_t(4)) RESET TO tmin (i_t(1)) !')
       ttret = ttmin
    END IF

    IF (ttmax < ttmin) THEN
       CALL RGMSG(substr, RGMLE,  'i_t(3) IN NAMELIST OR tmax AT ', .false.)
       CALL RGMSG(substr, RGMLEC, 'SUBROUTINECALL MUST BE >= i_t(1) OR tmin!')
    END IF
    IF ( .NOT.( (ttret >= ttmin).AND.(ttret <= ttmax) ) ) THEN
       CALL RGMSG(substr, RGMLE,                                              &
            'NCREGRID IN INTERFACE MODE ONLY RETURNS DATA, IF', .false.)
       CALL RGMSG(substr, RGMLEC,                                             &
            '   i_t(1) <= i_t(4) <= i_t(3)  IN THE NAMELIST, OR', .false.)
       CALL RGMSG(substr, RGMLEC,                                             &
            '   tmin   <= tret   <= tmax    AT SUBROUTINECALL!', .false.)
       CALL RGMSG(substr, RGMLEC,                                             &
            'PLEASE CORRECT NAMELIST (MOST PROBABLY i_t(4)), OR', .false.)
       CALL RGMSG(substr, RGMLEC,                                             &
            'PARAMETERS IN SUBROUTINECALL (tret) ACCORDINGLY!')
    END IF

    ! check gridfile time stepping
    ! note: in principle it is possible to specify a grdfile (in the namelist)
    !       in interface mode; this overwrites the grid specified by the
    !       interface.
    IF (g_t(1) <= 0) THEN
       CALL RGMSG(substr, RGMLE, 'g_t(1) IN NAMELIST MUST BE >0 !')
    END IF
    IF (g_t(2) <= 0) THEN
       CALL RGMSG(substr, RGMLE, 'g_t(2) IN NAMELIST MUST BE >0 !')
    END IF
    IF (g_t(3) <= g_t(1)) THEN
       CALL RGMSG(substr, RGMLW, &
            'GRDFILE TIME STEP ',g_t(1),' WILL BE USED')
       CALL RGMSG(substr, RGMLWC, 'FOR ALL INFILE TIME STEPS !')
    END IF
    ! check outfile time stepping
    IF (o_t(1) <= 0) THEN
       CALL RGMSG(substr, RGMLE, 'o_t(1) IN NAMELIST MUST BE >0 !')
    END IF
    IF (o_t(1) <= 0) THEN
       CALL RGMSG(substr, RGMLE, 'o_t(2) IN NAMELIST MUST BE >0 !')
    END IF
    ! o_t(3) is not used
    ! --- END TIME CONTROL ---------------------------------

    ! get variable infos
    CALL PARSE_VARSTR(varstr, ivar, ovar, SOURCEGRID, SOURCEGRIDstr, DESTGRID &
         , DESTGRIDstr, MAPPING, MAPstr, NORMALIZING, NORMstr                 &
         ,TRIM(gi0%file))
    nvars = SIZE(ivar)
    ALLOCATE(xivar(nvars), STAT=status)
    CALL ERRMSG(substr, status, 3)

    ! get infos of vector variables
    nvecs = 0
    IF (TRIM(vecstr) /= '') THEN
       CALL PARSE_VECSTR(vecstr,xvar,yvar)
       nvecs = SIZE(xvar)
       CALL RGMSG(substr, RGMLVL, &
            'PARSED ',nvecs,' VECTORS')
    END IF

    ! start time loop
    tg = g_t(1) - g_t(2)
    to = o_t(1) - o_t(2)
    !--------------------------------------------------------------------
    DO tc = ttmin, ttmax, ttstep

       ! copy grids from namelist input
       CALL COPY_GEOHYBGRID(gi, gi0)
       CALL COPY_GEOHYBGRID(ggrid(atm), ggatm0)
       CALL COPY_GEOHYBGRID(ggrid(oces), ggoces0)
       CALL COPY_GEOHYBGRID(ggrid(oceu), ggoceu0)
       CALL COPY_GEOHYBGRID(ggrid(ocev), ggocev0)

       tg = tg + g_t(2)         ! next gridfile timestep
       to = to + o_t(2)         ! next outputfile timestep
       IF ((tg > g_t(3)).AND.(g_t(3) > 0)) THEN
          tg = g_t(1)
       END IF
       gi%t = tc
       DO i=1,4
          ggrid(i)%t = tg
       END DO
       CALL RGMSG(substr, RGMLVL, 'CURRENT TIME STEP: ',tc,' ')
       CALL RGMSG(substr, RGMLVL, 'GRID TIME STEP: ',tg,' ')

       ! TODO: time information messages

       ! get input grid
       CALL RGMSG(substr, RGMLVL, 'IMPORTING INPUT GRID...')
       CALL RGMSG(substr, RGMLVL, '... IMPORTING FROM '''//TRIM(gi%file)//'''')
       CALL IMPORT_GEOHYBGRID(gi)

       ! TODO: check input grid ??

       ! get destination grids
       DO i=1,4
          CALL RGMSG(substr, RGMLVL, 'IMPORTING '''//TRIM(gridstring(i))//    &
               ''' OUTPUT GRID...')
          CALL RGMSG(substr, RGMLVL                                           &
               , '... IMPORTING FROM '''//TRIM(ggrid(i)%file)//'''')
          CALL IMPORT_GEOHYBGRID(ggrid(i))
       END DO

       IF (lint) THEN  ! adjust output time stepping for interface mode
       CALL RGMSG(substr, RGMLVL                                              &
            , 'USING INPUT TIME STEPS')
          to = (tc - ttmin)/ttstep + 1
!          go%t = gi%t
          go%t = to
       ELSE
          ! TODO: time info from grid file ... later? / search to
          to = tg
          go%t = tg
       END IF

       ! get input variabels
       CALL RGMSG(substr, RGMLVL                                              &
            , 'IMPORTING ',nvars,' VARIABLE(S) FROM FILE '''                  &
            //TRIM(gi%file)//'''')
       CALL RGMSG(substr, RGMLVL, &
            '('''//TRIM(gi%timem%name)//''' STEP: ',tc,')')
       mvars = 0          ! number of effective/grid-conform variables
       lpress = .false.   ! is pressure variable already in varlist?

       DO i=1, nvars      ! Loop over variable names
          CALL RGMSG(substr, RGMLVL,' ... '''//TRIM(ivar(i))//''' ...')
          IF (ASSOCIATED(gi%timem%dim)) THEN
             setuid = gi%timem%dim(1)%id
             CALL RGMSG(substr, RGMLVL                                         &
                  ,'TIME ID SET TO: ',setuid,' ')
          ELSE
             setuid = NULL_DIMID
             CALL RGMSG(substr, RGMLVL                                         &
                  ,'NO TIME ID SET !!!')
          END IF
          CALL IMPORT_NCVAR(tvar, tc, varname=TRIM(ivar(i))                   &
               , file=TRIM(gi%file)                                           &
               , setuid=setuid)
          CALL WRITE_NCVAR(tvar,.false.)
          ! check
          CALL RGMSG(substr, RGMLIC, '  ... CHECKING FOR GRID CONFORMITY ...')
          CALL CHECK_NCVAR_ON_GEOHYBGRID(tvar, gi, dims, axes, ok)
          ! not needed here
          axes(:) = 0
          DEALLOCATE(dims, STAT=status)
          CALL ERRMSG(substr, status, 4)
          NULLIFY(dims)

          IF (ok) THEN
             ! rename
             IF (TRIM(ivar(i)) /= TRIM(ovar(i))) THEN
                CALL RGMSG(substr, RGMLVLC                                    &
                     , '  ... RENAMING TO '''//TRIM(ovar(i))//'''... ')
                CALL RENAME_NCVAR(tvar, TRIM(ovar(i)))
                CALL ADD_NCATT(tvar, 'GT_ORIGINAL_NAME', vs=TRIM(ivar(i)))
             END IF
             ! adjust time axis
             ! if input_time (= lint) == .false.
             ! variable must get correct time-dimension-ID
             IF (.NOT.lint) THEN
                IF (.NOT.(setuid == NULL_DIMID)) THEN
                   DO j=1, tvar%ndims
                      IF (tvar%dim(j)%fuid) THEN
                         EXIT
                      END IF
                   END DO
                   CALL RGMSG(substr, RGMLVMC,'  ... ADJUSTING TIME-DIM.: ')
                   CALL RGMSG(substr, RGMLVMC,'      '                        &
                        //TRIM(tvar%dim(j)%name)//' (ID = ',j,') -> ')
                   CALL RGMSG(substr, RGMLVMC,'      '                        &
                        //TRIM(ggrid(DESTGRID(i))%timem%dim(1)%name)//        &
                        ' (ID = ',setuid,') ...')
                   CALL COPY_NCDIM(tvar%dim(j), ggrid(DESTGRID(i))%timem%dim(1))
                END IF
             END IF
             ! save in list
             mvars = mvars + 1
             ! change indices to new (effective variable position)
             CALL COPY_NCVAR(xivar(mvars), tvar)
             SOURCEGRID(mvars)    = SOURCEGRID(i)          ! shift
             SOURCEGRIDstr(mvars) = SOURCEGRIDstr(i)       ! ./.
             DESTGRID(mvars)      = DESTGRID(i)
             DESTGRIDstr(mvars)   = DESTGRIDstr(i)
             MAPPING(mvars)       = MAPPING(i)
             MAPstr(mvars)        = MAPstr(i)
             NORMALIZING(mvars)   = NORMALIZING(i)
             NORMstr(mvars)       = NORMstr(i)
          ELSE
             CALL  RGMSG(substr, RGMLVMC                                      &
                  , '  ... NOT CONSISTENT WITH GRID! IGNORING !')
          END IF
          ! check if surface pressure will already be regridded if i_ps in nml
          ! -> ps is already in variable list?
          IF (TRIM(gi%ps%name) /= '') THEN
             lpress = .true.
             IF (TRIM(tvar%name) == TRIM(gi%ps%name)) THEN
                lpress = .false.
             END IF
          END IF

          ! clean
          CALL INIT_NCVAR(tvar)
       END DO        ! loop over variable names

       ! if surface pressure is needed but not in variable list,
       ! it has to be added
       IF (lpress) THEN
          ! allocate and initialise appropriate arrays
          ALLOCATE(xitvar(mvars+1), STAT=status)
          CALL ERRMSG(substr, status, 5)
          ALLOCATE(TEMPSOURCE(mvars+1), STAT=status)
          CALL ERRMSG(substr, status, 6)
          ALLOCATE(TEMPDEST(mvars+1), STAT=status)
          CALL ERRMSG(substr, status, 7)
          ALLOCATE(TEMPMAP(mvars+1), STAT=status)
          CALL ERRMSG(substr, status, 8)
          ALLOCATE(TEMPNORM(mvars+1), STAT=status)
          CALL ERRMSG(substr, status, 9)
          DO i=1, mvars
             CALL INIT_NCVAR(xitvar(i))
             CALL COPY_NCVAR(xitvar(i),xivar(i))
             CALL INIT_NCVAR(xivar(i))
             TEMPSOURCE(i)=SOURCEGRID(i)
             TEMPDEST(i)=DESTGRID(i)
             TEMPMAP(i)=MAPPING(i)
             TEMPNORM(i)=NORMALIZING(i)
          END DO
          ! find sourcegrid of pressure variable
          TEMPSOURCE(mvars+1)=SOURCEGRID(1)
          IF (SOURCEGRID(1) == atm) THEN
             TEMPDEST(mvars+1)= oces
          ELSE
             TEMPDEST(mvars+1)= atm
          END IF
          TEMPMAP(i)=map_type_conserv
          TEMPNORM(i)=norm_opt_none
          ! deallocate old arrays
          DEALLOCATE(xivar, STAT=status)
          CALL ERRMSG(substr,status,10)
          DEALLOCATE(SOURCEGRID, STAT=status)
          CALL ERRMSG(substr,status,11)
          DEALLOCATE(DESTGRID, STAT=status)
          CALL ERRMSG(substr,status,12)
          DEALLOCATE(MAPPING, STAT=status)
          CALL ERRMSG(substr,status,13)
          DEALLOCATE(NORMALIZING, STAT=status)
          CALL ERRMSG(substr,status,14)

          ! get pressure variable
          CALL INIT_NCVAR(xitvar(mvars+1))
          CALL IMPORT_NCVAR(xitvar(mvars+1), tc, varname=TRIM(gi%ps%name)     &
               ,file=TRIM(gi%file)                                            &
               ,setuid=setuid)
          mvars = mvars + 1
          ! set up new arrays
          ALLOCATE(xivar(mvars), STAT=status)
          CALL ERRMSG(substr, status, 15)
          ALLOCATE(SOURCEGRID(mvars), STAT=status)
          CALL ERRMSG(substr, status, 16)
          ALLOCATE(DESTGRID(mvars), STAT=status)
          CALL ERRMSG(substr, status, 17)
          ALLOCATE(MAPPING(mvars), STAT=status)
          CALL ERRMSG(substr, status, 18)
          ALLOCATE(NORMALIZING(mvars), STAT=status)
          CALL ERRMSG(substr, status, 19)
          ! copy varstr variables
          DO i=1,mvars
             CALL INIT_NCVAR(xivar(i))
             CALL COPY_NCVAR(xivar(i),xitvar(i))
             CALL INIT_NCVAR(xitvar(i))
             SOURCEGRID(i)=TEMPSOURCE(i)
             DESTGRID(i)=TEMPDEST(i)
             MAPPING(i)=TEMPMAP(i)
             NORMALIZING(i)=TEMPNORM(i)
          END DO
          ! clean up
          DEALLOCATE(xitvar, STAT=status)
          CALL ERRMSG(substr,status,20)
          DEALLOCATE(TEMPSOURCE, STAT=status)
          CALL ERRMSG(substr,status,21)
          DEALLOCATE(TEMPDEST, STAT=status)
          CALL ERRMSG(substr,status,22)
          DEALLOCATE(TEMPMAP, STAT=status)
          CALL ERRMSG(substr,status,23)
          DEALLOCATE(TEMPNORM, STAT=status)
          CALL ERRMSG(substr,status,24)
       END IF

       !-----------------------------------------------------------------------
       ! CONTROL
       SELECT CASE(TCTRL)
       !-----------------------------------------------------------------------
       ! SCAN AND RETURN DATA
       CASE(GT_SCAN)
          CALL RGMSG(substr, RGMLVL, 'SCAN MODE...')
          IF (PRESENT(var)) THEN
             CALL RGMSG(substr, RGMLVM, '... RETURNING RAW DATA')
             ALLOCATE(var(mvars), STAT=status)
             CALL ERRMSG(substr, status, 25)
             DO i=1, mvars
                CALL INIT_NCVAR(var(i))
                CALL COPY_NCVAR(var(i),xivar(i))
             END DO
          END IF
          IF (PRESENT(grid)) THEN
             CALL RGMSG(substr, RGMLVM, '...RETURNING INPUT GRID')
             CALL COPY_GEOHYBGRID(grid,gi)
          END IF
          CALL RGMSG(substr, RGMLVL, '...DONE (SCAN MODE)')

       !-----------------------------------------------------------------------
       ! REGRIDDING MODE
       CASE(GT_PROC)
          CALL RGMSG(substr, RGMLVL, 'REGRIDDING MODE...')

          ! first time step
          IF (tc == ttmin) THEN

             ! initialization of INTCALCTRAFO and vector indices

             ALLOCATE(xvaridx(nvecs), STAT=status)
             CALL ERRMSG(substr, status, 26)
             ALLOCATE(yvaridx(nvecs), STAT=status)
             CALL ERRMSG(substr, status, 27)
             ALLOCATE(in_nlon(mvars), STAT=status)
             CALL ERRMSG(substr, status, 28)
             ALLOCATE(in_nz(mvars), STAT=status)
             CALL ERRMSG(substr, status, 29)
             ALLOCATE(in_lev(mvars), STAT=status)
             CALL ERRMSG(substr, status, 30)

             ! initialize number of vector indices
             nidx=0

             DO i=1, mvars ! loop over variables

                ! set every transformation, that has to be done = true
                INTCALCTRAFO(SOURCEGRID(i), DESTGRID(i)                       &
                     , MAPPING(i), NORMALIZING(i)) = .true.

                DO j=1,nvecs ! loop over vectors

                   ! find index of vector components in input variables
                   IF (TRIM(xvar(j)) == TRIM(xivar(i)%name)) THEN
                      xvaridx(j)=i
                      nidx = nidx + 1
                      CALL RGMSG(substr, RGMLVL                               &
                           , 'PUTTING '''//TRIM(xivar(i)%name)//              &
                           ''' TO X-VECTOR LIST')
                   END IF
                   IF (TRIM(yvar(j)) == TRIM(xivar(i)%name)) THEN
                      yvaridx(j)=i
                      nidx = nidx + 1
                      CALL RGMSG(substr, RGMLVL                               &
                           , 'PUTTING '''//TRIM(xivar(i)%name)//              &
                           ''' TO Y-VECTOR LIST')
                   END IF

                END DO ! loop over vectors

             END DO ! loop over variables

             ! test if all vector variables in varstr
             IF (nidx /= (nvecs * 2)) THEN
                CALL RGMSG(substr, RGMLE                                      &
                     , 'ERROR DURING VECTOR INITIALISATION.', .false.)
                CALL RGMSG(substr, RGMLEC                                     &
                     , 'ALL VECTOR COMPONENTS HAVE TO BE INCLUDED')
                CALL RGMSG(substr, RGMLEC                                     &
                     , 'IN THE NAMELIST VAR LIST.')
             END IF


             DO i = 1, 4 ! loop over grids
                IF (TRIM(ggrid(i)%file) /= '') THEN
                   nlon(i) = ggrid(i)%lonm%dim(1)%len
                   IF (ggrid(i)%latm%ndims == 2) THEN
                      nlat(i) = ggrid(i)%latm%dim(2)%len
                   ELSE
                      nlat(i) = ggrid(i)%latm%dim(1)%len
                   END IF
                   IF ((ggrid(i)%lonm%ndims == 1)                             &
                        .AND. (ggrid(i)%latm%ndims == 1)) THEN

                      ! produce 2-d lon/lat grids for scrip remapping
                      CALL INIT_NARRAY(longitude(i), 2                        &
                           , (/nlon(i), nlat(i)/), VTYPE_DOUBLE)
                      CALL INIT_NARRAY(latitude(i), 2                         &
                           , (/nlon(i), nlat(i)/), VTYPE_DOUBLE)
                      
                      ! put everything in double precision narrays
                      natype = ggrid(i)%lonm%dat%type
                      SELECT CASE(natype)
                      CASE(VTYPE_REAL)
                         DO j = 1, nlat(i)
                            DO k = 1, nlon(i)
                               longitude(i)%vd(POSITION((/nlon(i), nlat(i)/)  &
                                    , (/k, j/))) = ggrid(i)%lonm%dat%vr(k)
                            END DO
                         END DO
                      CASE(VTYPE_DOUBLE)
                         DO j = 1, nlat(i)
                            DO k = 1, nlon(i)
                               longitude(i)%vd(POSITION((/nlon(i), nlat(i)/)  &
                                    , (/k, j/))) = ggrid(i)%lonm%dat%vd(k)
                            END DO
                         END DO
                      CASE(VTYPE_UNDEF)
                         CALL RGMSG(substr, RGMLE                             &
                              , 'ERROR! '''//TRIM(gridstring(i))//            &
                              ''' GRID DIMENSION UNDEFINED.'                  &
                              , .false.)
                      CASE DEFAULT
                         CALL RGMSG(substr, RGMLE                             &
                              , 'ERROR! REGRIDDING ONLY POSSIBLE FOR'         &
                              , .false.)
                         CALL RGMSG(substr, RGMLEC                            &
                              , 'DIMENSIONS OF TYPE REAL OR DOUBLE PRECISION !')
                      END SELECT
                      natype = ggrid(i)%latm%dat%type
                      SELECT CASE(natype)
                      CASE(VTYPE_REAL)
                         DO j = 1, nlon(i)
                            DO k = 1, nlat(i)
                               latitude(i)%vd(POSITION((/nlon(i), nlat(i)/)   &
                                    , (/j, k/))) = ggrid(i)%latm%dat%vr(k)
                            END DO
                         END DO
                      CASE(VTYPE_DOUBLE)
                         DO j = 1, nlon(i)
                            DO k = 1, nlat(i)
                               latitude(i)%vd(POSITION((/nlon(i), nlat(i)/)   &
                                    , (/j, k/))) = ggrid(i)%latm%dat%vd(k)
                            END DO
                         END DO
                      CASE(VTYPE_UNDEF)
                         CALL RGMSG(substr, RGMLE                             &
                              , 'ERROR! '''//TRIM(gridstring(i))//            &
                              ''' GRID DIMENSION UNDEFINED.'                  &
                              , .false.)
                      CASE DEFAULT
                         CALL RGMSG(substr, RGMLE                             &
                              , 'ERROR! REGRIDDING ONLY POSSIBLE FOR'         &
                              , .false.)
                         CALL RGMSG(substr, RGMLEC                            &
                              , 'DIMENSIONS OF TYPE REAL OR DOUBLE PRECISION !')
                      END SELECT
                   ELSE
                      ! already 2-d: copy...
                      CALL COPY_NARRAY(longitude(i), ggrid(i)%lonm%dat)
                      CALL COPY_NARRAY(latitude(i), ggrid(i)%latm%dat)
                      ! and put to double precision narrays...
                      CALL DOUBLE_NARRAY(longitude(i))
                      CALL DOUBLE_NARRAY(latitude(i))
                   END IF

                   ! for SCRIP library, fields have to be in radians...
                   ldegrees = .false.
                   DO j = 1, ggrid(i)%lonm%natts
                      ! test, if the fields are in degrees...
                      attstr = ''
                      IF (TRIM(ggrid(i)%lonm%att(j)%name) == 'units') THEN
                         DO k=1,ggrid(i)%lonm%att(j)%len
                            attstr(k:k) = ggrid(i)%lonm%att(j)%dat%vc(k)
                         END DO
                         IF (INDEX(TRIM(attstr),'degrees') /= 0)              &
                              ldegrees = .true.
                      END IF
                   END DO
                   IF (MAXVAL(longitude(i)%vd) > (2._dp * pi)                 &
                        .OR. (MINVAL(longitude(i)%vd) < (-2._dp * pi)))       &
                        ldegrees = .true.
                   ! put fields to radians
                   IF (ldegrees) THEN
                      CALL RADIANS(longitude(i),.true.)
                      CALL RADIANS(latitude(i))
                      CALL RESTRICT_RADIANS(longitude(i)%vd, latitude(i)%vd)
                   END IF

                   ! get grid corners...
                   IF (ggrid(i)%loni%ndims /= 0) THEN
                      ! check longitude corners...
                      ! if there are corners already in the grid...
                      axes(:)=0
                      ! find out where x,y axis are
                      CALL CHECK_NCVAR_ON_GEOHYBGRID(ggrid(i)%loni,           &
                           ggrid(i), dims, axes, ok)
                      gg_nlon(i) = ggrid(i)%loni%dat%dim(axes(GORD_LON))
                      gg_nlat(i) = ggrid(i)%loni%dat%dim(axes(GORD_LAT))
                      IF (axes(GORD_LEV) .ne. 0) THEN
                         dims(axes(GORD_LEV)) = 0
                      END IF
                      ! length of free dimensions (number of corners)
                      gg_nc(i) = 1
                      DO j = 1, ggrid(i)%loni%ndims
                         IF (dims(j) == 0)                                    &
                              gg_nc(i) = gg_nc(i) * ggrid(i)%loni%dim(j)%len
                      END DO
                      ! put corners in 1-d array
                      CALL PACK_GEOHYBGRID_NCVAR(ggrid(i)%loni, dims          &
                           , axes, corners(1))
                      ! make every data double precision
                      CALL DOUBLE_NARRAY(corners(1)%dat)

                      IF (ggrid(i)%lati%ndims /= 0) THEN
                         ! check latitude corners...
                         ! should be sorted like longitude corners,...
                         ! so we put latitude corners also in 1-d array
                         CALL PACK_GEOHYBGRID_NCVAR(ggrid(i)%lati, dims       &
                              , axes,corners(2))
                         ! and double precision
                         CALL DOUBLE_NARRAY(corners(2)%dat)
                      ELSE
                         ! lati has to be in grid!!
                         CALL RGMSG(substr, RGMLE                            &
                              , 'ERROR! lati has to be in grid!', .false.)
                      END IF
                      
                      ! for SCRIP library, fields have to be in radians...
                      ldegrees = .false.
                      DO j = 1, ggrid(i)%loni%natts
                         ! test, if the fields are in degrees...
                         attstr = ''
                         IF (TRIM(ggrid(i)%loni%att(j)%name) == 'units') THEN
                            DO k=1,ggrid(i)%loni%att(j)%len
                               attstr(k:k) = ggrid(i)%loni%att(j)%dat%vc(k)
                            END DO
                            IF (INDEX(TRIM(attstr),'degrees') /= 0)           &
                                 ldegrees = .true.
                         END IF
                      END DO
                      IF (MAXVAL(corners(1)%dat%vd) > (2._dp * pi)            &
                           .OR. (MINVAL(corners(1)%dat%vd) < (-2._dp * pi)))  &
                           ldegrees = .true.
                      ! put fields to radians
                      IF (ldegrees) THEN
                         CALL RADIANS(corners(1)%dat,.true.)
                         CALL RADIANS(corners(2)%dat)
                         CALL RESTRICT_RADIANS(corners(1)%dat%vd              &
                              , corners(2)%dat%vd)
                      END IF

                      ! new fields for input to SCRIP remapping...
                      ! they have to be of form array(#corners,data)...
                      ALLOCATE(corner_lon(i)%corners(gg_nc(i)                 &
                           , gg_nlon(i)*gg_nlat(i)), STAT=status)
                      CALL ERRMSG(substr, status, 31)
                      ALLOCATE(corner_lat(i)%corners(gg_nc(i)                 &
                           , gg_nlon(i)*gg_nlat(i)), STAT=status)
                      CALL ERRMSG(substr, status, 32)

                      DO j = 1, gg_nc(i)
                         corner_lon(i)%corners(j,:) =                         &
                              corners(1)%dat%vd((j-1)*gg_nlon(i)*gg_nlat(i)+1 &
                              : j*gg_nlon(i)*gg_nlat(i))
                         corner_lat(i)%corners(j,:) =                         &
                              corners(2)%dat%vd((j-1)*gg_nlon(i)*gg_nlat(i)+1 &
                              : j*gg_nlon(i)*gg_nlat(i))
                      END DO

                      ! clean up
                      CALL INIT_NARRAY(corners(1)%dat)
                      CALL INIT_NARRAY(corners(2)%dat)
                      CALL INIT_NCVAR(corners(1))
                      CALL INIT_NCVAR(corners(2))
                      DEALLOCATE(dims, STAT=status)
                      CALL ERRMSG(substr, status, 32)

                   ELSE
                      ! else calculate corners...
                      gg_nc(i) = 4
                      gg_nlon(i) = nlon(i)
                      gg_nlat(i) = nlat(i)
                      ALLOCATE(corner_lon(i)%corners(gg_nc(i)                 &
                           , gg_nlon(i)*gg_nlat(i)), STAT=status)
                      CALL ERRMSG(substr, status, 33)
                      ALLOCATE(corner_lat(i)%corners(gg_nc(i)                 &
                           , gg_nlon(i)*gg_nlat(i)), STAT=status)
                      CALL ERRMSG(substr, status, 34)
                      CALL INIT_NARRAY(gridcellarea(i), 2                     &
                           , (/nlon(i), nlat(i)/), VTYPE_DOUBLE)
                      ! calculation of grid cell areas...
                      ! only tested for regular atmospheric grid now!!
                      DO j = 1, gg_nlon(i)
                         DO k = 1, gg_nlat(i)
                            IF (j == 1) THEN
                               posjm = POSITION((/gg_nlon(i), gg_nlat(i)/)    &
                                    , (/gg_nlon(i), k/))
                            ELSE
                               posjm = POSITION((/gg_nlon(i), gg_nlat(i)/)    &
                                    , (/j-1, k/))
                            END IF
                            IF (j == gg_nlon(i)) THEN
                               posjp = POSITION((/gg_nlon(i), gg_nlat(i)/)    &
                                    , (/1, k/))
                            ELSE
                               posjp = POSITION((/gg_nlon(i), gg_nlat(i)/)    &
                                    ,(/j+1, k/))
                            END IF
                            IF (k == 1) THEN
                               poskm = POSITION((/gg_nlon(i), gg_nlat(i)/)    &
                                    , (/j, gg_nlat(i)/))
                            ELSE
                               poskm = POSITION((/gg_nlon(i), gg_nlat(i)/)    &
                                    , (/j, k-1/))
                            END IF
                            IF (k == gg_nlat(i)) THEN
                               poskp = POSITION((/gg_nlon(i), gg_nlat(i)/)    &
                                    , (/j, 1/))
                            ELSE
                               poskp = POSITION((/gg_nlon(i), gg_nlat(i)/)    &
                                    ,(/j, k+1/))
                            END IF
                            pos = POSITION((/gg_nlon(i), gg_nlat(i)/)         &
                                 , (/j, k/))

                            IF (ABS(longitude(i)%vd(pos)                      &
                                 - longitude(i)%vd(posjp)) > pi) THEN
                               DO l = 1, 2
                                  corner_lon(i)%corners(l,pos) =              &
                                       0.5_dp * (longitude(i)%vd(pos)         &
                                       + longitude(i)%vd(posjp))              &
                                       + pi
                               END DO
                            ELSE
                               DO l = 1, 2
                                  corner_lon(i)%corners(l,pos) =              &
                                       0.5_dp * (longitude(i)%vd(pos)         &
                                       + longitude(i)%vd(posjp))
                               END DO
                            END IF
                            IF (ABS(longitude(i)%vd(posjm)                    &
                                 - longitude(i)%vd(pos)) > pi) THEN
                               DO l = 3, 4
                                  corner_lon(i)%corners(l,pos) =              &
                                       0.5_dp * (longitude(i)%vd(posjm)       &
                                       + longitude(i)%vd(pos))                &
                                       + pi
                               END DO
                            ELSE
                               DO l = 3, 4
                                  corner_lon(i)%corners(l,pos) =              &
                                       0.5_dp * (longitude(i)%vd(posjm)       &
                                       + longitude(i)%vd(pos))
                               END DO
                            END IF
                            DO l = 1, 4, 3
                               corner_lat(i)%corners(l,pos) =                 &
                                    0.5_dp * (latitude(i)%vd(pos)             &
                                    + latitude(i)%vd(poskp))
                            END DO
                            DO l = 2, 3
                               corner_lat(i)%corners(l,pos) =                 &
                                    0.5_dp * (latitude(i)%vd(pos)             &
                                    + latitude(i)%vd(poskm))
                            END DO
                            IF (k == gg_nlat(i)) THEN
                               corner_lat(i)%corners(1,pos) = -pih
                               corner_lat(i)%corners(4,pos) = -pih
                            END IF
                            IF (k == 1) THEN
                               corner_lat(i)%corners(2,pos) = pih
                               corner_lat(i)%corners(3,pos) = pih
                            END IF

                            sinlap = SIN(0.5_dp *                             &
                                 (corner_lat(i)%corners(1,pos)                &
                                 + corner_lat(i)%corners(4,pos)))
                            sinlam = SIN(0.5_dp *                             &
                                 (corner_lat(i)%corners(2,pos)                &
                                 + corner_lat(i)%corners(3,pos)))
                            lonp = 0.5_dp * (corner_lon(i)%corners(1,pos)     &
                                 + corner_lon(i)%corners(2,pos))
                            lonm = 0.5_dp * (corner_lon(i)%corners(3,pos)     &
                                 + corner_lon(i)%corners(4,pos))
                            gridcellarea(i)%vd(pos) = ABS(sinlam - sinlap)    &
                                 * ABS(lonm - lonp)
                         END DO
                      END DO
                   END IF
                END IF
             END DO ! loop over grids

             ! calculate grid cell areas...
             CALL INIT_NARRAY(gridcellarea(2), 2, (/gg_nlon(2), gg_nlat(2)/)  &
                  , VTYPE_DOUBLE)
             CALL GRIDAREA(gridcellarea(2), longitude(3)%vd, latitude(3)%vd   &
                  , longitude(4)%vd, latitude(4)%vd, 2)
             CALL INIT_NARRAY(gridcellarea(3), 2, (/gg_nlon(3), gg_nlat(3)/)  &
                  , VTYPE_DOUBLE)
             CALL GRIDAREA(gridcellarea(3), longitude(2)%vd, latitude(2)%vd   &
                  , corner_lon(2)%corners(1,:), corner_lat(2)%corners(1,:), 3)
             CALL INIT_NARRAY(gridcellarea(4), 2, (/gg_nlon(4), gg_nlat(4)/)  &
                  , VTYPE_DOUBLE)
             CALL GRIDAREA(gridcellarea(4), corner_lon(2)%corners(1,:)        &
                  , corner_lat(2)%corners(1,:), longitude(2)%vd               &
                  , latitude(2)%vd, 4)

             DO i = 1, 4
                ALLOCATE(gridmask(i)%mask(gg_nlon(i)*gg_nlat(i)), STAT=status)
                CALL ERRMSG(substr, status, 35)
                gridmask(i)%mask(:) = .true.
             END DO

             CALL REGRID(INTCALCTRAFO, longitude, latitude                    &
                  , corner_lon, corner_lat, gridcellarea, gridmask            &
                  , gg_nlon, gg_nlat, transformation)

             DO i = 1, 4
                IF (ASSOCIATED(corner_lon(i)%corners)) THEN
                   DEALLOCATE(corner_lon(i)%corners, STAT=status)
                   CALL ERRMSG(substr, status, 36)
                END IF
                IF (ASSOCIATED(corner_lat(i)%corners)) THEN
                   DEALLOCATE(corner_lat(i)%corners, STAT=status)
                   CALL ERRMSG(substr, status, 37)
                END IF
                IF (ASSOCIATED(gridmask(i)%mask)) THEN
                   DEALLOCATE(gridmask(i)%mask, STAT=status)
                   CALL ERRMSG(substr, status, 38)
                END IF
                CALL INIT_NARRAY(gridcellarea(i))
             END DO

          END IF ! first time step

          ! mz_bk_20100924+
          DO i = 1, 4
             CALL INIT_NCVAR(ggrid(i)%loni)
             CALL INIT_NCVAR(ggrid(i)%lati)
          END DO
          ! mz_bk_20100924-


          !--------------------------------------------------------------------
          ! regridding of variables
          CALL RGMSG(substr, RGMLVL, &
               'INITIALIZING REGRIDDING OF VARIABLES')
          ALLOCATE(nai(mvars), STAT=status)
          CALL ERRMSG(substr, status, 39)
          ALLOCATE(nao(mvars), STAT=status)
          CALL ERRMSG(substr, status, 40)
          DO i=1, mvars
             CALL INIT_NARRAY(nai(i)%dat)
             CALL INIT_NARRAY(nao(i)%dat)
          END DO

          CALL RGMSG(substr, RGMLVL, &
               'PREPARING REGRIDDING OF ',mvars,' VARIABLE(S)...')

          DO i=1, mvars ! loop over variables
             axes(:)=0
             CALL CHECK_NCVAR_ON_GEOHYBGRID(xivar(i), gi, dims, axes, ok)

             in_nlon(i) = xivar(i)%dat%dim(axes(GORD_LON))
             in_lev(i)=0
             IF (axes(GORD_LEV) .ne. 0) THEN
                in_lev(i) = xivar(i)%dat%dim(axes(GORD_LEV))
                dims(axes(GORD_LEV)) = 0
             END IF
             !length of free dimensions
             in_nz(i) = 1
             DO j=1, xivar(i)%ndims
                IF (dims(j) == 0) in_nz(i) = in_nz(i) * xivar(i)%dim(j)%len
             END DO

             CALL PACK_GEOHYBGRID_NCVAR(xivar(i), dims          &
                  , axes, nai(i))
             ! make every data double precision
             CALL DOUBLE_NARRAY(nai(i)%dat)

             IF ((SOURCEGRID(i) == oces)                                      &
                  .AND. ((in_nlon(i) - 2) == gg_nlon(oces))) THEN
                ! temporaray out variable
                ALLOCATE(out(gg_nlon(oces)*gg_nlat(oces)*in_nz(i)), STAT=status)
                CALL ERRMSG(substr, status, 41)
                ! set mask to cut corners
                ALLOCATE(mask(PRODUCT(xivar(i)%dat%dim(:))), STAT=status)
                 CALL ERRMSG(substr, status, 42)
                 mask(:) = .true.
                 ! loop over all free dimensions (levs, n)
                 dimvec = (/in_nlon(i), gg_nlat(oces), in_nz(i)/)
                 DO j = 1, in_nz(i)
                    DO k = 1, gg_nlat(oces)
                       cornervec = (/1, k, j/)
                       mask(POSITION(dimvec, cornervec)) = .false.
                       cornervec = (/in_nlon(i), k, j/)
                       mask(POSITION(dimvec, cornervec)) = .false.
                    END DO
                 END DO
                ! cut corners
                CALL RGMSG(substr, RGMLVL, &
                     'CUTTING CORNERS!!')
                out = PACK(nai(i)%dat%vd,mask)
                ! copy back to variable
                nai(i)%dat%dim(dims(axes(GORD_LON))) = gg_nlon(oces)
                nai(i)%dim(dims(axes(GORD_LON)))%len = gg_nlon(oces)
                DEALLOCATE(nai(i)%dat%vd, STAT=status)
                CALL ERRMSG(substr, status, 43)
                ALLOCATE(nai(i)%dat%vd(gg_nlon(oces)*gg_nlat(oces)*in_nz(i))  &
                     , STAT=status)
                CALL ERRMSG(substr, status, 44)
                nai(i)%dat%vd = out
                ! deallocate temporary variables
                DEALLOCATE(out, STAT=status)
                CALL ERRMSG(substr, status, 45)
                DEALLOCATE(mask, STAT=status)
                CALL ERRMSG(substr, status, 46)
             END IF
             DEALLOCATE(dims, STAT=status)
             CALL ERRMSG(substr, status, 32)
          END DO
          CALL RGMSG(substr, RGMLVL, &
               'VARIABLE(S) PREPARED!')

          ! vector rotation for ocean -> atmosphere must be done before
          ! regridding. first rotate ocean vectors in north-east atmosphere
          ! grid orientation, then regrid.

          CALL RGMSG(substr, RGMLVL, &
               'PREPARING REGRIDDING OF ',nvecs,' VECTOR(S)...')
          DO i=1,nvecs ! loop over vectors
             ! ocean -> atmosphere
             IF ((DESTGRID(xvaridx(i))==atm)                                  &
                  .AND. (DESTGRID(yvaridx(i))==atm)) THEN
                ALLOCATE(xvecdat(gg_nlon(OCES), gg_nlat(OCES)), STAT=status)
                CALL ERRMSG(substr, status, 47)
                ALLOCATE(yvecdat(gg_nlon(OCES), gg_nlat(OCES)), STAT=status)
                CALL ERRMSG(substr, status, 48)
                ALLOCATE(outlon(gg_nlon(OCES), gg_nlat(OCES)), STAT=status)
                CALL ERRMSG(substr, status, 49)
                ALLOCATE(outlat(gg_nlon(OCES), gg_nlat(OCES)), STAT=status)
                CALL ERRMSG(substr, status, 50)
                DO ii = 1, in_nz(xvaridx(i))
                   DO j = 1, gg_nlon(OCES)
                      DO k = 1, gg_nlat(OCES)
                         xvecdat(j, k) = nai(xvaridx(i))%dat%vd((k-1)         &
                              * gg_nlon(OCES) + j + (ii-1) * gg_nlon(OCES)    &
                              * gg_nlat(OCES))
                         yvecdat(j, k) = nai(yvaridx(i))%dat%vd((k-1)         &
                              * gg_nlon(OCES) + j + (ii-1) * gg_nlon(OCES)    &
                              * gg_nlat(OCES))
                         outlon(j, k) = longitude(OCES)%vd(POSITION(          &
                              (/gg_nlon(OCES), gg_nlat(OCES)/), (/j, k/)))
                         outlat(j, k) = latitude(OCES)%vd(POSITION(           &
                              (/gg_nlon(OCES), gg_nlat(OCES)/), (/j, k/)))
                      END DO
                   END DO
                   CALL rotate_2_ne(xvecdat, yvecdat, outlon, outlat          &
                        , gg_nlon(OCES), gg_nlat(OCES))
                   CALL RGMSG(substr, RGMLVL, 'ROTATING VECTOR COMPONENTS')
                   CALL WRITE_NARRAY(nai(xvaridx(i))%dat,.false.)
                   DO j = 1, gg_nlon(OCES)
                      DO k=1 , gg_nlat(OCES)
                         nai(xvaridx(i))%dat%vd((k-1) * gg_nlon(OCES)         &
                              + j + (ii-1) * gg_nlon(OCES) * gg_nlat(OCES))   &
                              = xvecdat(j, k)
                         nai(yvaridx(i))%dat%vd((k-1) * gg_nlon(OCES)         &
                              + j + (ii-1) * gg_nlon(OCES) * gg_nlat(OCES))   &
                              = yvecdat(j, k)
                      END DO
                   END DO
                END DO
                DEALLOCATE(xvecdat, STAT=status)
                CALL ERRMSG(substr, status, 51)
                DEALLOCATE(yvecdat, STAT=status)
                CALL ERRMSG(substr, status, 52)
                DEALLOCATE(outlon, STAT=status)
                CALL ERRMSG(substr, status, 53)
                DEALLOCATE(outlat, STAT=status)
                CALL ERRMSG(substr, status, 54)
             END IF
          END DO ! vectors
          CALL RGMSG(substr, RGMLVL, &
               'VECTOR(S) PREPARED!')

          ! clean up
          DO i = 1, 4
             CALL INIT_NARRAY(longitude(i))
             CALL INIT_NARRAY(latitude(i))
          END DO
          
          ALLOCATE(srcvec(3))
          ALLOCATE(destvec(3))

          CALL RGMSG(substr, RGMLVL, &
               'REMAPPING OF ',mvars,' VARIABLE(S)...')
          DO i=1,mvars ! loop over variables
             ! setup array for output
             axes(:)=0
             CALL CHECK_NCVAR_ON_GEOHYBGRID(xivar(i), gi, dims, axes, ok)
             IF (axes(GORD_LEV) .ne. 0) THEN
                dims(axes(GORD_LEV)) = 0
             END IF
             !length of free dimensions
             in_nz(i) = 1
             DO j=1, xivar(i)%ndims
                IF (dims(j) == 0) in_nz(i) = in_nz(i) * xivar(i)%dim(j)%len
             END DO
             CALL COPY_NCVAR(nao(i),nai(i))
             CALL INIT_NARRAY(nao(i)%dat, 3,                                  &
                  (/nlon(DESTGRID(i)), nlat(DESTGRID(i)), in_nz(i)/)          &
                  , VTYPE_DOUBLE)
             CALL COPY_NCDIM(nao(i)%dim(axes(GORD_LON))                       &
                  , ggrid(DESTGRID(i))%lonm%dim(1))
             IF (ggrid(DESTGRID(i))%latm%ndims == 1) THEN
                CALL COPY_NCDIM(nao(i)%dim(axes(GORD_LAT))                    &
                     , ggrid(DESTGRID(i))%latm%dim(1))
             ELSE
                CALL COPY_NCDIM(nao(i)%dim(axes(GORD_LAT))                    &
                     , ggrid(DESTGRID(i))%latm%dim(2))
             END IF
                
             nao(i)%dat%vd(:) = 0.0
             ! regridding
             SELECT CASE (MAPPING(i))
             CASE (map_type_conserv)

                ! first order conservative transformation
                DO k = 1, in_nz(i) !loop over free levels
                   DO j = 1, transformation(SOURCEGRID(i), DESTGRID(i)        &
                        , MAPPING(i), NORMALIZING(i))%transform%nn
                      ! address in destination array
                      destadd = transformation(SOURCEGRID(i), DESTGRID(i)     &
                           , MAPPING(i), NORMALIZING(i))%transform%destadd(j)
                      ! address in source array
                      srcadd = transformation(SOURCEGRID(i), DESTGRID(i)      &
                           , MAPPING(i),NORMALIZING(i))%transform%srcadd(j)
                      ! normalizing factor
                      factor = 1.0_dp
                      ! factor for destarea
                      IF (NORMALIZING(i).eq.norm_opt_dstarea) THEN
                         IF (transformation(SOURCEGRID(i), DESTGRID(i)        &
                              , MAPPING(i), NORMALIZING(i))%transform         &
                              %frac(destadd).gt.tiny) THEN
                            factor = 1.0_dp/transformation(                   &
                                 SOURCEGRID(i), DESTGRID(i), MAPPING(i)       &
                                 , NORMALIZING(i))%transform%frac(destadd)
                         ELSE
                            factor = 0.0_dp
                         END IF
                      END IF
                      ! factor for no normalizing
                      IF (NORMALIZING(i).eq.norm_opt_none) THEN
                         IF (transformation(SOURCEGRID(i)                     &
                              , DESTGRID(i)                                   &
                              , MAPPING(i)                                    &
                              , NORMALIZING(i))%transform%frac(destadd)       &
                              .gt.tiny) THEN
                            factor = 1.0_dp / (transformation(SOURCEGRID(i)   &
                                 , DESTGRID(i)                                &
                                 , MAPPING(i)                                 &
                                 , NORMALIZING(i))%transform%frac(destadd)    &
                                 * transformation(SOURCEGRID(i)               &
                                 , DESTGRID(i)                                &
                                 , MAPPING(i)                                 &
                                 , NORMALIZING(i))%transform%destarea(destadd))
                         ELSE
                            factor = 0.0_dp
                         END IF
                      END IF
                      dim_in(1) = nai(i)%dat%dim(1)
                      dim_in(2) = nai(i)%dat%dim(2)
                      dim_in(3) = 1
                      CALL ELEMENT(dim_in, srcadd, srcvec)
                      dim_in(3) = in_nz(i)
                      srcvec(3) = k
                      dim_out(1) = nao(i)%dat%dim(1)
                      dim_out(2) = nao(i)%dat%dim(2)
                      dim_out(3) = 1
                      CALL ELEMENT(dim_out, destadd, destvec)
                      dim_out(3) = in_nz(i)
                      destvec(3) = k
                      value = nao(i)%dat%vd(POSITION(dim_out, destvec))
                      nao(i)%dat%vd(POSITION(dim_out, destvec))               &
                           = value                                            &
                           + transformation(SOURCEGRID(i)                     &
                           , DESTGRID(i)                                      &
                           , MAPPING(i)                                       &
                           , NORMALIZING(i)                                   &
                           )%transform%weight(1, j)                           &
                           * nai(i)%dat%vd(POSITION(dim_in, srcvec))          &
                           * factor
                   END DO
                END DO !loop over levels

             CASE (map_type_bilinear)
                factor = 1.0_dp
                IF (NORMALIZING(i).eq.norm_opt_dstarea) THEN
                END IF
                IF (NORMALIZING(i).eq.norm_opt_none) THEN
                END IF
             CASE (map_type_bicubic)
                factor = 1.0_dp
                IF (NORMALIZING(i).eq.norm_opt_dstarea) THEN
                END IF
                IF (NORMALIZING(i).eq.norm_opt_none) THEN
                END IF
             CASE (map_type_distwgt)
                factor = 1.0_dp
                IF (NORMALIZING(i).eq.norm_opt_dstarea) THEN
                END IF
                IF (NORMALIZING(i).eq.norm_opt_none) THEN
                END IF
             END SELECT

             ! clean up
             DEALLOCATE(dims, STAT=status)
             CALL ERRMSG(substr, status, 32)

          END DO ! variables

          DEALLOCATE(srcvec, destvec)

          CALL RGMSG(substr, RGMLVL, &
               'REMAPPING OF VARIABLE(S) DONE!')

          ! vector rotation
          CALL RGMSG(substr, RGMLVL, &
               'ROTATING ',nvecs,' VECTOR(S)...')
          DO i = 1, nvecs ! loop over vectors
             ! atmosphere -> ocean
             IF ((SOURCEGRID(xvaridx(i)) == ATM)                              &
                  .AND. (SOURCEGRID(yvaridx(i)) == ATM)) THEN
                ALLOCATE(xvecdat(gg_nlon(OCES), gg_nlat(OCES)), STAT=status)
                CALL ERRMSG(substr, status, 55)
                ALLOCATE(yvecdat(gg_nlon(OCES), gg_nlat(OCES)), STAT=status)
                CALL ERRMSG(substr, status, 56)
                ALLOCATE(vec_tmp_x(gg_nlon(OCES), gg_nlat(OCES)), STAT=status)
                CALL ERRMSG(substr, status, 57)
                ALLOCATE(vec_tmp_y(gg_nlon(OCES), gg_nlat(OCES)), STAT=status)
                CALL ERRMSG(substr, status, 58)

                DO ii = 1, in_nz(xvaridx(i))
                   DO j = 1, gg_nlon(OCES)
                      DO k = 1, gg_nlat(OCES)
                         xvecdat(j, k) = nao(xvaridx(i))%dat%vd((k-1)         &
                              * gg_nlon(OCES) + j + (ii-1) * gg_nlon(OCES)    &
                              * gg_nlat(OCES))
                         yvecdat(j, k) = nao(yvaridx(i))%dat%vd((k-1)         &
                              * gg_nlon(OCES) + j + (ii-1) * gg_nlon(OCES)    &
                              * gg_nlat(OCES))
                      END DO
                   END DO
                   vec_tmp_x = xvecdat
                   vec_tmp_y = yvecdat
                   SELECT CASE (DESTGRID(xvaridx(i)))
                   CASE (OCES)
                      ALLOCATE(outlon(gg_nlon(OCES), gg_nlat(OCES))           &
                           , STAT=status)
                      CALL ERRMSG(substr, status, 59)
                      DO j= 1, gg_nlon(OCES)
                         DO k = 1, gg_nlat(OCES)
                            outlon(j, k) = longitude(OCES)%vd(POSITION(       &
                           (/gg_nlon(OCES), gg_nlat(OCES)/), (/j, k/)))
                            outlat(j, k) = latitude(OCES)%vd(POSITION(        &
                           (/gg_nlon(OCES), gg_nlat(OCES)/), (/j, k/)))
                         END DO
                      END DO
                   CASE (OCEU)
                      ALLOCATE(outlon(gg_nlon(OCEU), gg_nlat(OCEU))           &
                           , STAT=status)
                      CALL ERRMSG(substr, status, 60)
                      DO j= 1, gg_nlon(OCEU)
                         DO k = 1, gg_nlat(OCEU)
                            outlon(j, k) = longitude(OCEU)%vd(POSITION(       &
                           (/gg_nlon(OCEU), gg_nlat(OCEU)/), (/j, k/)))
                            outlat(j, k) = latitude(OCEU)%vd(POSITION(        &
                           (/gg_nlon(OCEU), gg_nlat(OCEU)/), (/j, k/)))
                         END DO
                      END DO
                   CASE (OCEV)
                      ALLOCATE(outlon(gg_nlon(OCEV), gg_nlat(OCEV))           &
                           , STAT=status)
                      CALL ERRMSG(substr, status, 61)
                      DO j= 1, gg_nlon(OCEV)
                         DO k = 1, gg_nlat(OCEV)
                            outlon(j, k) = longitude(OCEV)%vd(POSITION(       &
                           (/gg_nlon(OCEV), gg_nlat(OCEV)/), (/j, k/)))
                            outlat(j, k) = latitude(OCEV)%vd(POSITION(        &
                           (/gg_nlon(OCEV), gg_nlat(OCEV)/), (/j, k/)))
                         END DO
                      END DO
                   END SELECT

                   CALL RGMSG(substr, RGMLVL, 'ROTATING VECTOR COMPONENTS')
                   CALL rotate_u(xvecdat, vec_tmp_y, outlon, outlat           &
                        , gg_nlon(OCES), gg_nlat(OCES))

                   CALL rotate_v(vec_tmp_x, yvecdat, outlon, outlat           &
                        , gg_nlon(OCES), gg_nlat(OCES))

                   DO j = 1, gg_nlon(OCES)
                      DO k = 1, gg_nlat(OCES)
                         nao(xvaridx(i))%dat%vd((k-1) * gg_nlon(OCES)         &
                              + j + (ii-1) * gg_nlon(OCES) * gg_nlat(OCES))   &
                              = xvecdat(j, k)
                         nao(yvaridx(i))%dat%vd((k-1) * gg_nlon(OCES)         &
                              + j + (ii-1) * gg_nlon(OCES) * gg_nlat(OCES))   &
                              = yvecdat(j,k)
                      END DO
                   END DO
                END DO
                DEALLOCATE(xvecdat, STAT=status)
                CALL ERRMSG(substr, status, 62)
                DEALLOCATE(yvecdat, STAT=status)
                CALL ERRMSG(substr, status, 63)
                DEALLOCATE(vec_tmp_x, STAT=status)
                CALL ERRMSG(substr, status, 64)
                DEALLOCATE(vec_tmp_y, STAT=status)
                CALL ERRMSG(substr, status, 65)
                DEALLOCATE(outlon, STAT=status)
                CALL ERRMSG(substr, status, 66)
                DEALLOCATE(outlat, STAT=status)
                CALL ERRMSG(substr, status, 67)
             END IF
          END DO ! vectors
          CALL RGMSG(substr, RGMLVL, &
               'ROTATING OF VECTOR(S) DONE!')

          !--------------------------------------------------------------------
          ! first  time step
          IF (tc == ttmin) THEN

             CALL RGMSG(substr, RGMLVL, 'INITIALIZING OUTPUT GRID')
             CALL INIT_GEOHYBGRID (go)

             go%file = TRIM(outfile)

             ! balance grids
             DO j = 1, 4 ! loop over grids

                IF (ANY(MASK=INTCALCTRAFO(:,j,:,:))) THEN
                   IF (j == atm) THEN
                      CALL BALANCE_GEOHYBGRID(ggrid(j),go)
                   ELSE
                      CALL RGMSG(substr, RGMLVL, 'BALANCING '                 &
                           //TRIM(gridstring(j)))
                      CALL INIT_NCVAR(tvar)
                      CALL IMPORT_NCVAR(tvar                                  &
                           , varname=TRIM(ggrid(j)%lonm%name)                 &
                           , file=TRIM(ggrid(j)%file))

                      IF (lover) THEN
                         IF (tvar%name /= '')                                 &
                              CALL OVERLAP(tvar, nlon(j), nlat(j))
                      END IF
                      
                      IF (go%lonm%name /= '') THEN
                         IF (go%lonm%name == ggrid(j)%lonm%name)              &
                              CALL RENAME_NCVAR(tvar                          &
                              , TRIM(ggrid(j)%lonm%name)                      &
                              //'_'//TRIM(gridstring(j)))
                         IF (go%clonm%name /= '') THEN
                            IF (go%clonm%name == ggrid(j)%lonm%name)          &
                                 CALL RENAME_NCVAR(tvar                       &
                                 , TRIM(ggrid(j)%lonm%name)//'_'              &
                                 //TRIM(gridstring(j)))
                            CALL COPY_NCVAR(go%cloni, tvar)
                         ELSE
                            CALL COPY_NCVAR(go%clonm, tvar)
                         END IF
                      ELSE
                         CALL COPY_NCVAR(go%lonm, tvar)
                      END IF
                      CALL INIT_NCVAR(tvar)
                      CALL IMPORT_NCVAR(tvar                                  &
                           , varname=TRIM(ggrid(j)%latm%name)                 &
                           , file=TRIM(ggrid(j)%file))
                      IF (lover) THEN
                         IF (tvar%name /= '')                                 &
                              CALL OVERLAP(tvar, nlon(j), nlat(j))
                      END IF
                      
                      IF (go%latm%name /= '') THEN
                         IF (go%latm%name == ggrid(j)%latm%name)              &
                              CALL RENAME_NCVAR(tvar                          &
                              , TRIM(ggrid(j)%latm%name)                      &
                              //'_'//TRIM(gridstring(j)))
                         IF (go%clatm%name /= '') THEN
                            IF (go%clatm%name == ggrid(j)%latm%name)          &
                                 CALL RENAME_NCVAR(tvar                       &
                                 , TRIM(ggrid(j)%latm%name)//'_'              &
                                 //TRIM(gridstring(j)))
                            CALL COPY_NCVAR(go%cloni, tvar)
                         ELSE
                            CALL COPY_NCVAR(go%clatm, tvar)
                         END IF
                      ELSE
                         CALL COPY_NCVAR(go%latm, tvar)
                      END IF
                      CALL INIT_NCVAR(tvar)
                   END IF
                END IF
             END DO ! loop over grids

             DO i = 1, 4
                CALL INIT_GEOHYBGRID(ggrid(i))
             END DO

          END IF ! first time step
          !--------------------------------------------------------------------

          CALL RGMSG(substr, RGMLVL, 'BALANCING GI <-> GO')
          CALL BALANCE_GEOHYBGRID(gi,go)

          CALL BALANCE_GEOHYBGRID_TIME(gi,go,lint)
          go%t = (tc-ttmin)/ttstep + 1
          CALL EXPORT_GEOHYBGRID(go)

          !--------------------------------------------------------------------
          ! first time step
          IF (tc == ttmin) THEN

             ! export the global attributes
             ! get old attribute
             CALL IMPORT_NCATT(gtgatt, varid=NF90_GLOBAL                      &
                  ,attname = 'GT_HISTORY'                                     &
                  ,file = TRIM(go%file), lnostop=.true.)
             ! check if empty or char-att
             IF ((gtgatt%dat%type == VTYPE_UNDEF).OR.                         &
                  (gtgatt%dat%type == VTYPE_CHAR)) THEN
                ! append new data
                CALL CAT_NARRAY(gtgatt%dat, nmlstr)
                gtgatt%xtype = NF90_CHAR
                gtgatt%len   = gtgatt%dat%dim(1)
             ELSE  ! exists with non-char type
                CALL RGMSG(substr, RGMLW, &
                     'ATTRIBUTE '''//TRIM(gtgatt%name)//'''')
                CALL RGMSG(substr, RGMLWC, &
                     'EXISTS ALREADY, BUT IS NOT OF TYPE CHAR !')
             END IF
             ! write attribute
             CALL EXPORT_NCATT(gtgatt, file=TRIM(go%file)                     &
                  ,clobber=.true.)
             ! clean up
             CALL INIT_NARRAY(nmlstr)
             
          END IF
          !--------------------------------------------------------------------
          
          ! setup output variables
          ALLOCATE(xovar(mvars), STAT=status)
          CALL ERRMSG(substr, status, 68)
          DO i=1, mvars ! loop over variables

             ! mz_bk_20100924+
             ! IF (axes(GORD_LEV) .ne. 0) THEN
             !    dims(axes(GORD_LEV)) = 0
             ! END IF
             ! mz_bk_20100924-

             ! if output variables should be on overlapping ocean gird
             ! add the two overlapping longitude lines
             IF ((lover) .AND. (DESTGRID(i) == oces)) THEN

                ALLOCATE(out((gg_nlon(oces) + 2) * gg_nlat(oces) * in_nz(i))  &
                     , STAT=status)
                CALL ERRMSG(substr, status, 69)
                ALLOCATE(mask((gg_nlon(oces) + 2) * gg_nlat(oces) * in_nz(i)) &
                     , STAT=status)
                CALL ERRMSG(substr, status, 70)
                ALLOCATE(stripes((gg_nlon(oces) + 2)*gg_nlat(oces)*in_nz(i))  &
                     , STAT=status)
                CALL ERRMSG(substr, status, 71)
                mask(:) = .true.
                ! loop over all free dimensions (levs, n)
                dimvec = (/gg_nlon(oces) + 2, gg_nlat(oces), in_nz(i)/)
                DO j = 1, in_nz(i)
                   DO k = 1, gg_nlat(oces)
                      cornervec = (/1, k, j/)
                      mask(POSITION(dimvec, cornervec)) = .false.
                      cornervec = (/gg_nlon(oces) + 2, k, j/)
                      mask(POSITION(dimvec, cornervec)) = .false.
                   END DO
                END DO

                ! cut corners
                stripes(:)=0.0
                dimvec = (/gg_nlon(oces) + 2, gg_nlat(oces), in_nz(i)/)
                DO j = 1, in_nz(i)
                   DO k = 1, gg_nlat(oces)
                      cornervec = (/1, k, j/)
                      stripes(POSITION(dimvec, cornervec))                    &
                           = nao(i)%dat%vd(POSITION((/gg_nlon(oces)           &
                           , gg_nlat(oces), in_nz(i)/), (/gg_nlon(oces)       &
                           , k, j/)))
                      cornervec = (/gg_nlon(oces) + 2, k, j/)
                      stripes(POSITION(dimvec, cornervec))                    &
                           = nao(i)%dat%vd(POSITION((/gg_nlon(oces)           &
                           , gg_nlat(oces), in_nz(i)/), (/1, k, j/)))
                   END DO
                END DO
                out = UNPACK(nao(i)%dat%vd,mask,stripes)

                nao(i)%dat%dim(dims(axes(GORD_LON))) = gg_nlon(oces) + 2
                nao(i)%dim(dims(axes(GORD_LON)))%len = gg_nlon(oces) + 2
                DEALLOCATE(nao(i)%dat%vd, STAT=status)
                CALL ERRMSG(substr, status, 72)
                ALLOCATE(nao(i)%dat%vd((gg_nlon(oces) + 2) * gg_nlat(oces)    &
                     * in_nz(i)), STAT=status)
                CALL ERRMSG(substr, status, 73)

                nao(i)%dat%vd=out

                DEALLOCATE(out, STAT=status)
                CALL ERRMSG(substr, status, 74)
                DEALLOCATE(mask, STAT=status)
                CALL ERRMSG(substr, status, 75)
                DEALLOCATE(stripes, STAT=status)
                CALL ERRMSG(substr, status, 76)

             END IF

             CALL CHECK_NCVAR_ON_GEOHYBGRID(xivar(i), gi, dims, axes, ok)
             CALL BALANCE_GEOHYBGRID_NCVAR(xivar(i),axes,go,xovar(i))

             IF (in_lev(i) == 0) THEN
                DEALLOCATE(dims)
                NULLIFY(dims)
                ALLOCATE(dims(3))
                dims(1) = 1
                dims(2) = 2
                dims(3) = 3
             ELSE
                axes(3)=3
                dims(1) = 1
                dims(2) = 2
                dims(3) = 3
                dims(4) = 4
                ! TODO: find out, which axis is the lev axis:
                ! -> compare len with in_lev -> dims(3) -> # (now 3)
                !, dims(4) -> # (now 4)
             END IF

             CALL COPY_NCDIM(xovar(i)%dim(dims(axes(GORD_LON)))               &
                  , nao(i)%dim(dims(axes(GORD_LON))))
             CALL COPY_NCDIM(xovar(i)%dim(dims(axes(GORD_LAT)))               &
                  , nao(i)%dim(dims(axes(GORD_LAT))))

             IF (QCMP_NCDIM(xovar(i)%dim(axes(GORD_LEV))                      &
                  , xivar(i)%dim(axes(GORD_LEV))) == 0) THEN
                CALL COPY_NCDIM(xovar(i)%dim(dims(axes(GORD_LEV)))            &
                     , xivar(i)%dim(dims(axes(GORD_LEV))))
             END IF

             ! TODO   !!! NOT WORKING FOR EVERYTHING!!!
             ALLOCATE(ovec(SIZE(dims)))
             DO j = 1, SIZE(dims)
                ovec(j)=xovar(i)%dim(j)%len
             END DO

             IF (in_lev(i) .ne. 0) THEN
             ELSE
                
             END IF


             CALL INIT_NARRAY(xovar(i)%dat, SIZE(dims), ovec, VTYPE_DOUBLE)

             DEALLOCATE(ovec)

             CALL PACK_GEOHYBGRID_NCVAR(nao(i),dims,axes,xovar(i),.true.)

             !test types and fill accordingly
             natype = xivar(i)%dat%type
             SELECT CASE(natype)
             CASE (VTYPE_REAL)
                ALLOCATE(xovar(i)%dat%vr(PRODUCT(xovar(i)%dat%dim(:)))        &
                     , STAT=status)
                CALL ERRMSG(substr, status, 77)
                xovar(i)%dat%vr = xovar(i)%dat%vd
                DEALLOCATE(xovar(i)%dat%vd, STAT=status)
                CALL ERRMSG(substr, status, 78)
                xovar(i)%dat%type = VTYPE_REAL
             CASE (VTYPE_DOUBLE)
             CASE(VTYPE_UNDEF)
                CALL RGMSG(substr, RGMLE                                      &
                     , 'ERROR! VARIABLE DIMENSION UNDEFINED.'                 &
                     , .false.)
             CASE DEFAULT
                CALL RGMSG(substr, RGMLE                                      &
                     , 'ERROR! REGRIDDING ONLY POSSIBLE FOR'                  &
                     , .false.)
                CALL RGMSG(substr, RGMLEC                                     &
                     , 'DIMENSIONS OF TYPE REAL OR DOUBLE PRECISION !')
             END SELECT

             xovar(i)%ustep = go%t
             CALL RGMSG(substr, RGMLVL  &
                  , 'OUTPUT TIMESTEP SET TO:',xovar(i)%ustep,' ')
             CALL EXPORT_NCVAR(xovar(i), file=go%file)

             ! clean up
             CALL INIT_NCVAR(nao(i))
             CALL INIT_NCVAR(nai(i))
             CALL INIT_NCVAR(xovar(i))
             CALL INIT_NCVAR(xivar(i))
             DEALLOCATE(dims, STAT=status)
             CALL ERRMSG(substr, status, 96)
             NULLIFY(dims)

          END DO ! variables

          ! clean up

          DEALLOCATE(nao, STAT=status)
          CALL ERRMSG(substr, status, 79)
          DEALLOCATE(nai, STAT=status)
          CALL ERRMSG(substr, status, 80)
          DEALLOCATE(xovar, STAT=status)
          CALL ERRMSG(substr, status, 81)

          CALL INIT_GEOHYBGRID(gi)

       !-----------------------------------------------------------------------
       CASE(GT_STOP)
       ! do nothing, never reached!!
       !-----------------------------------------------------------------------
       CASE DEFAULT
       !-----------------------------------------------------------------------
       END SELECT
       !-----------------------------------------------------------------------
    END DO ! time loop
    !--------------------------------------------------------------------------

    ! clean up

    IF (ASSOCIATED(xvaridx)) DEALLOCATE(xvaridx)
    IF (ASSOCIATED(yvaridx)) DEALLOCATE(yvaridx)
    IF (ASSOCIATED(xvar)) DEALLOCATE(xvar)
    IF (ASSOCIATED(yvar)) DEALLOCATE(yvar)

    DEALLOCATE(xivar, STAT=status)
    CALL ERRMSG(substr, status, 82)
    DEALLOCATE(ivar, STAT=status)
    CALL ERRMSG(substr, status, 83)
    DEALLOCATE(ovar, STAT=status)
    CALL ERRMSG(substr, status, 84)
    DEALLOCATE(SOURCEGRID, STAT=status)
    CALL ERRMSG(substr, status, 85)
    DEALLOCATE(DESTGRID, STAT=status)
    CALL ERRMSG(substr, status, 86)
    DEALLOCATE(MAPPING, STAT=status)
    CALL ERRMSG(substr, status, 87)
    DEALLOCATE(NORMALIZING, STAT=status)
    CALL ERRMSG(substr, status, 88)
    DEALLOCATE(SOURCEGRIDstr, STAT=status)
    CALL ERRMSG(substr, status, 89)
    DEALLOCATE(DESTGRIDstr, STAT=status)
    CALL ERRMSG(substr, status, 90)
    DEALLOCATE(MAPstr, STAT=status)
    CALL ERRMSG(substr, status, 91)
    DEALLOCATE(NORMstr, STAT=status)
    CALL ERRMSG(substr, status, 92)
    DEALLOCATE(in_nlon, STAT=status)
    CALL ERRMSG(substr, status, 93)
    DEALLOCATE(in_nz, STAT=status)
    CALL ERRMSG(substr, status, 94)
    DEALLOCATE(in_lev, STAT=status)
    CALL ERRMSG(substr, status, 95)

    DO i=1,4
       DO j=1,4
          DO ii=1,4
             DO jj=1,3
                IF(ASSOCIATED(transformation(i,j,ii,jj)%transform)) THEN
                   IF(ASSOCIATED(transformation(i,j,ii,jj)%transform%srcadd)) &
                        DEALLOCATE(transformation(i,j,ii,jj)%transform%srcadd)
                   IF(ASSOCIATED(transformation(i,j,ii,jj)                    &
                        %transform%destadd))                                  &
                        DEALLOCATE(transformation(i,j,ii,jj)%transform%destadd)
                   IF(ASSOCIATED(transformation(i,j,ii,jj)                    &
                        %transform%srcarea))                                  &
                        DEALLOCATE(transformation(i,j,ii,jj)%transform%srcarea)
                   IF(ASSOCIATED(transformation(i,j,ii,jj)                    &
                        %transform%destarea))                                 &
                        DEALLOCATE(transformation(i,j,ii,jj)                  &
                        %transform%destarea)
                   IF(ASSOCIATED(transformation(i,j,ii,jj)%transform%weight)) &
                        DEALLOCATE(transformation(i,j,ii,jj)%transform%weight)
                   IF(ASSOCIATED(transformation(i,j,ii,jj)%transform%frac))   &
                        DEALLOCATE(transformation(i,j,ii,jj)%transform%frac)
                   DEALLOCATE(transformation(i,j,ii,jj)%transform, STAT=status)
                   CALL ERRMSG(substr, status, 97)
                END IF
             END DO
          END DO
       END DO
    END DO

    CALL INIT_GEOHYBGRID(gi)
    CALL INIT_GEOHYBGRID(go)

  END SUBROUTINE GRIDTRAFO_CONTROL
  !------------------------------------------------------------------

  !------------------------------------------------------------------
  SUBROUTINE OVERLAP(longitude, nlon, nlat)
    ! returns longitude on an overlapping grid.
    ! it adds a stripe in front and in back of longitude.

    ! I/O
    TYPE(NCVAR), INTENT(INOUT) :: longitude
    INTEGER, INTENT(IN)        :: nlon, nlat


    REAL(DP), DIMENSION(:), ALLOCATABLE :: out, stripes
    LOGICAL, DIMENSION(:), ALLOCATABLE  :: mask

    TYPE(NCVAR) :: unpacked

    INTEGER :: ndims, natype, nz
    INTEGER :: i, j
    INTEGER, DIMENSION(3) :: axes
    INTEGER, DIMENSION(:), ALLOCATABLE :: dims
    INTEGER, DIMENSION(3) :: dimvec, cornervec

    INTEGER :: status

    ! Local
    CHARACTER(LEN=*), PARAMETER   :: substr = 'OVERLAP'

    ndims = longitude%ndims
    natype = longitude%dat%type
    CALL DOUBLE_NARRAY(longitude%dat)

    CALL INIT_NCVAR(unpacked)

    CALL COPY_NCVAR(unpacked, longitude)

    axes(:)=0
    DO i = 1, ndims
       IF (longitude%dim(i)%len == nlon) axes(GORD_LON)=i
       IF (longitude%dim(i)%len == nlat) axes(GORD_LAT)=i
    END DO

    ALLOCATE (dims(longitude%ndims))
    dims(:)=0
    dims(axes(GORD_LON))=GORD_LON
    dims(axes(GORD_LAT))=GORD_LAT

    nz = 1
    DO i = 1, ndims
       IF (dims(i) == 0) nz = nz * longitude%dim(i)%len
    END DO

    CALL PACK_GEOHYBGRID_NCVAR(longitude, dims, axes, longitude)

    ALLOCATE(out((nlon + 2)*nlat*nz), STAT=status)
    CALL ERRMSG(substr, status, 1)
    ALLOCATE(mask((nlon + 2)*nlat*nz), STAT=status)
    CALL ERRMSG(substr, status, 2)
    ALLOCATE(stripes((nlon + 2)*nlat*nz), STAT=status)
    CALL ERRMSG(substr, status, 3)
    mask(:) = .true.
    stripes(:) = 0.0
    ! loop over all free dimensions (levs, n)
    dimvec = (/nlon + 2, nlat, nz/)
    DO i = 1, nz
       DO j = 1, nlat
          cornervec = (/1, j, i/)
          mask(POSITION(dimvec, cornervec)) = .false.
          stripes(POSITION(dimvec, cornervec)) =                              &
               longitude%dat%vd(POSITION((/nlon, nlat, nz/), (/nlon, j, i/)))
          cornervec = (/nlon + 2, j, i/)
          mask(POSITION(dimvec, cornervec)) = .false.
          stripes(POSITION(dimvec, cornervec)) =                              &
               longitude%dat%vd(POSITION((/nlon, nlat, nz/), (/1, j, i/)))
       END DO
    END DO
    out = UNPACK(longitude%dat%vd,mask,stripes)

    longitude%dim(GORD_LON)%len = nlon + 2

    CALL INIT_NARRAY(longitude%dat, 3, (/nlon + 2, nlat, nz/), VTYPE_DOUBLE)
    longitude%dat%vd = out

    unpacked%dim(1)%len = nlon + 2
    CALL INIT_NARRAY(unpacked%dat, 2, (/nlon + 2, nlat/), VTYPE_DOUBLE)
    CALL PACK_GEOHYBGRID_NCVAR(longitude, dims, axes, unpacked, .true.)

    CALL COPY_NCVAR(longitude,unpacked)

    SELECT CASE(natype)
    CASE (VTYPE_REAL)
       ALLOCATE(longitude%dat%vr(PRODUCT(longitude%dat%dim(:))), STAT=status)
       CALL ERRMSG(substr, status, 4)
       longitude%dat%vr = longitude%dat%vd
       DEALLOCATE(longitude%dat%vd, STAT=status)
       CALL ERRMSG(substr, status, 5)
       ! Important to change the type here !!
       longitude%dat%type = VTYPE_REAL
    CASE (VTYPE_DOUBLE)
    CASE(VTYPE_UNDEF)
       CALL RGMSG(substr, RGMLE                                      &
            , 'ERROR! VARIABLE DIMENSION UNDEFINED.'                 &
            , .false.)
    CASE DEFAULT
       CALL RGMSG(substr, RGMLE                                      &
            , 'ERROR! ONLY POSSIBLE FOR'                  &
            , .false.)
       CALL RGMSG(substr, RGMLEC                                     &
            , 'DIMENSIONS OF TYPE REAL OR DOUBLE PRECISION !')
    END SELECT
    
  END SUBROUTINE OVERLAP
  !------------------------------------------------------------------

  !------------------------------------------------------------------
  SUBROUTINE REGRID(INTCALCTRAFO, longitude, latitude, corner_lon, corner_lat &
       , gridcellarea, gridmask, nlon, nlat, transformation)

    IMPLICIT NONE

    LOGICAL, DIMENSION(4,4,4,3), INTENT(IN)             :: INTCALCTRAFO

    TYPE(NARRAY), DIMENSION(:), TARGET, INTENT(IN)      :: longitude, latitude

    TYPE(CORNERARRAY), DIMENSION(4), TARGET, INTENT(IN) :: corner_lon         &
         , corner_lat

    TYPE(NARRAY), DIMENSION(:), TARGET, INTENT(IN)      :: gridcellarea

    TYPE(MASKARRAY), DIMENSION(4), TARGET, INTENT(IN)   :: gridmask

    TYPE(TRAFO), DIMENSION(4,4,4,3), INTENT(INOUT)      :: transformation

    INTEGER, DIMENSION(4), INTENT(IN), OPTIONAL         :: nlon, nlat

    INTEGER, DIMENSION(4)                               :: nnlon, nnlat

    INTEGER :: i, j, jj, kk ! counter

    INTEGER :: status

    ! Local
    CHARACTER(LEN=*), PARAMETER   :: substr = 'REGRID'

    IF (.NOT. PRESENT(nlon)) THEN
       DO i = 1, 4
          nnlon(i) = longitude(i)%dim(1)
       END DO
    ELSE
       nnlon = nlon
    END IF

    IF (.NOT. PRESENT(nlat)) THEN
       DO i = 1, 4
          nnlat(i) = longitude(i)%dim(2)
       END DO
    ELSE
       nnlat = nlat
    END IF


    DO i = 1, 4
       IF (ANY(MASK=INTCALCTRAFO(i,:,:,:))) THEN
          CALL RGMSG(substr, RGMLVL, 'SOURCEGRID: '//TRIM(gridstring(i)))
          grid1_center_lon => longitude(i)%vd
          grid1_center_lat => latitude(i)%vd
          grid1_corner_lon => corner_lon(i)%corners
          grid1_corner_lat => corner_lat(i)%corners
          grid1_area_in => gridcellarea(i)%vd
          grid1_mask => gridmask(i)%mask
          DO j = 1, 4
             IF (ANY(MASK=INTCALCTRAFO(i,j,:,:))) THEN
                CALL RGMSG(substr, RGMLVL, 'DESTGRID: '//TRIM(gridstring(j)))
                grid2_center_lon => longitude(j)%vd
                grid2_center_lat => latitude(j)%vd
                grid2_corner_lon => corner_lon(j)%corners
                grid2_corner_lat => corner_lat(j)%corners
                grid2_area_in => gridcellarea(j)%vd
                grid2_mask => gridmask(j)%mask
                DO jj=1,4 ! loop over mapping types
                   map_type = jj
                   IF (jj.eq.map_type_conserv) THEN
                      luse_grid_centers = .false.
                   ELSE
                      luse_grid_centers = .true.
                   END IF
                   DO kk=1,3 ! loop over normalizing options
                      IF (INTCALCTRAFO(i,j,jj,kk)) THEN
                         norm_opt = kk
                         CALL REMAP_INIT(4, nnlon(i), nnlat(i)                &
                              , 4, nnlon(j), nnlat(j))
                         CALL BOUNDS_CALC
                         CALL REMAP_VARS(1)
                         SELECT CASE(jj)
                         CASE (map_type_conserv)
                            CALL RGMSG(substr, RGMLVL,                        &
                                 'SCRIPTRAFO: CONSERVATIVE REMAPPING')
                            CALL REMAP_CONSERV
                         CASE (map_type_bilinear)
                            CALL RGMSG(substr, RGMLE                          &
                                 , 'ERROR: BILINEAR REMAPPING', .false.)
                            CALL RGMSG(substr, RGMLEC                         &
                                 , 'NOT IMPLEMENTED (YET?).')
                            !CALL REMAP_BILIN
                         CASE (map_type_bicubic)
                            CALL RGMSG(substr, RGMLE                          &
                                 , 'ERROR: BICUBIC REMAPPING', .false.)
                            CALL RGMSG(substr, RGMLEC                         &
                                 , 'NOT IMPLEMENTED (YET?).')
                            !CALL REMAP_BICUB
                         CASE (map_type_distwgt)
                            CALL RGMSG(substr, RGMLE                          &
                                 , 'ERROR: DISTANCE WEIGHTED REMAPPING'       &
                                 , .false.)
                            CALL RGMSG(substr, RGMLEC                         &
                                 , 'NOT IMPLEMENTED (YET?).')
                            !CALL REMAP_DISTWGT
                         CASE DEFAULT
                            CALL RGMSG(substr, RGMLE                          &
                                 , 'ERROR: UNKNOWN REMAPPING TYPE.'           &
                                 , .false.)
                         END SELECT
                         CALL REMAP_VARS(2)
                         transformation(i,j,jj,kk)%info%source_grid           &
                              = i
                         transformation(i,j,jj,kk)%info%dest_grid=j
                         transformation(i,j,jj,kk)%info%map_method=jj
                         transformation(i,j,jj,kk)%info%norm_option=kk
                         ALLOCATE(transformation(i,j,jj,kk)                   &
                              %transform, STAT=status)
                         CALL ERRMSG(substr, status, 1)
                         transformation(i,j,jj,kk)%transform%nn               &
                              = num_links_map1
                         ALLOCATE(transformation(i,j,jj,kk)                   &
                              %transform%srcadd(num_links_map1)               &
                              , STAT=status)
                         CALL ERRMSG(substr, status, 2)
                         ALLOCATE(transformation(i,j,jj,kk)                   &
                              %transform%destadd(num_links_map1)              &
                              , STAT=status)
                         CALL ERRMSG(substr, status, 3)
                         ALLOCATE(transformation(i,j,jj,kk)                   &
                              %transform%srcarea(nnlon(i)*nnlat(i))           &
                              , STAT=status)
                         CALL ERRMSG(substr, status, 4)
                         ALLOCATE(transformation(i,j,jj,kk)                   &
                              %transform%destarea(nnlon(j)*nnlat(j))          &
                              , STAT=status)
                         CALL ERRMSG(substr, status, 5)
                         ALLOCATE(transformation(i,j,jj,kk)                   &
                              %transform%weight(num_wts,num_links_map1)       &
                              , STAT=status)
                         CALL ERRMSG(substr, status, 6)
                         ALLOCATE(transformation(i,j,jj,kk)                   &
                              %transform%frac(nnlon(j)*nnlat(j))              &
                              , STAT=status)
                         CALL ERRMSG(substr, status, 7)
                         transformation(i,j,jj,kk)%transform%srcadd           &
                              = grid1_add_map1
                         transformation(i,j,jj,kk)%transform%destadd          &
                              = grid2_add_map1
                         transformation(i,j,jj,kk)%transform%srcarea          &
                              = grid1_area
                         transformation(i,j,jj,kk)                            &
                              %transform%destarea = grid2_area
                         transformation(i,j,jj,kk)%transform%weight           &
                              = wts_map1
                         transformation(i,j,jj,kk)%transform%frac             &
                              = grid2_frac
                         CALL REMAP_DEALLOC
                      END IF
                   END DO ! normalizing options
                END DO    ! mapping types

             END IF
          END DO
       END IF
    END DO

  END SUBROUTINE REGRID
  !------------------------------------------------------------------

  !------------------------------------------------------------------
  SUBROUTINE GRIDAREA(gridcellarea, lon1, lat1, lon2, lat2, type)
    ! calculates gridcell area

    IMPLICIT NONE

    ! I/O
    TYPE(NARRAY), INTENT(INOUT) :: gridcellarea
    REAL(DP), DIMENSION(:), INTENT(IN) :: lon1, lat1, lon2, lat2
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: dlx, dly
    INTEGER, INTENT(IN) :: type ! type of grid

    INTEGER :: posn, posim, posip, posjm, posjp

    INTEGER, DIMENSION(4) :: pos

    INTEGER :: nlon, nlat

    INTEGER :: i,j
    
    REAL(DP) :: siphup, siphum, siphvp, siphvm, cophup, cophum                &
         , cophvp, cophvm, silaup, silaum, silavp, silavm, colaup, colaum     &
         , colavp, colavm                                                     &
         , xup, xum, xvp, xvm, yup, yum, yvp, yvm, zup, zum, zvp, zvm
    
    INTEGER :: status

    ! Local
    CHARACTER(LEN=*), PARAMETER   :: substr = 'GRIDAREA'

    nlon = gridcellarea%dim(1)
    nlat = gridcellarea%dim(2)

    ALLOCATE(dlx(nlon,nlat), STAT=status)
    CALL ERRMSG(substr, status, 1)
    ALLOCATE(dly(nlon,nlat), STAT=status)
    CALL ERRMSG(substr, status, 2)

    DO i = 1, nlon
       DO j = 1, nlat
          
          IF (i == 1) THEN
             posim = POSITION((/nlon, nlat/), (/nlon, j/))
          ELSE
             posim = POSITION((/nlon, nlat/), (/i-1, j/))
          END IF
          IF (i == nlon) THEN
             posip = POSITION((/nlon, nlat/), (/1, j/))
          ELSE
             posip = POSITION((/nlon, nlat/),(/i+1, j/))
          END IF
          IF (j == 1) THEN
             posjm = POSITION((/nlon, nlat/), (/i, nlat/))
          ELSE
             posjm = POSITION((/nlon, nlat/), (/i, j-1/))
          END IF
          IF (j == nlat) THEN
             posjp = POSITION((/nlon, nlat/), (/i, 1/))
          ELSE
             posjp = POSITION((/nlon, nlat/),(/i, j+1/))
          END IF
          posn = POSITION((/nlon, nlat/), (/i, j/))

          IF (type == 2) THEN
             pos(1) = posip
             pos(2) = posn
             pos(3) = posjp
             pos(4) = posn
          ELSE IF (type == 3) THEN
             pos(1) = posn
             pos(2) = posim
             pos(3) = posjp
             pos(4) = posn
          ELSE IF (type == 4) THEN
             pos(1) = posip
             pos(2) = posn
             pos(3) = posn
             pos(4) = posjm
          ELSE
             CALL RGMSG(substr, RGMLEC                         &
                  , 'TYPE NOT SUPPORTED.')
          END IF

          siphup = sin(lon1(pos(1)))
          siphum = sin(lon1(pos(2)))
          siphvp = sin(lon2(pos(3)))
          siphvm = sin(lon2(pos(4)))
          
          cophup = cos(lon1(pos(1)))
          cophum = cos(lon1(pos(2)))
          cophvp = cos(lon2(pos(3)))
          cophvm = cos(lon2(pos(4)))
          
          silaup = sin(lat1(pos(1)))
          silaum = sin(lat1(pos(2)))
          silavp = sin(lat2(pos(3)))
          silavm = sin(lat2(pos(4)))
          
          colaup = cos(lat1(pos(1)))
          colaum = cos(lat1(pos(2)))
          colavp = cos(lat2(pos(3)))
          colavm = cos(lat2(pos(4)))
          
          xup = silaup * cophup
          xum = silaum * cophum
          yup = colaup * cophup
          yum = colaum * cophum
          zup = siphup
          zum = siphum
          
          xvp = silavp * cophvp
          xvm = silavm * cophvm
          yvp = colavp * cophvp
          yvm = colavm * cophvm
          zvp = siphvp
          zvm = siphvm
          
          dlx(i,j)=max(1.0_dp, (radius_earth                                  &
               * acos(min((xup*xum+yup*yum+zup*zum),1.0_dp))))
          dly(i,j)=max(1.0_dp, (radius_earth                                  &
               * acos(min((xvp*xvm+yvp*yvm+zvp*zvm),1.0_dp))))
       END DO
    END DO

    DO i = 1, nlon
       DO j = 1, nlat
          pos = POSITION((/nlon, nlat/), (/i, j/))          
          gridcellarea%vd(pos) = dlx(i,j) * dly(i,j) / (radius_earth ** 2)
       END DO
    END DO
    
    DEALLOCATE(dlx, dly)

  END SUBROUTINE GRIDAREA

  !------------------------------------------------------------------
  SUBROUTINE radians(na, llon)

    IMPLICIT NONE

    ! I/O
    TYPE(NARRAY), INTENT(INOUT)   :: na
    LOGICAL, INTENT(IN), OPTIONAL :: llon

    LOGICAL                       :: ldeg

    ldeg = .false.

    IF (PRESENT(llon)) ldeg = llon

    IF (ldeg) THEN
       WHERE (na%vd < 0._dp)
          na%vd = na%vd + 360._dp
       END WHERE
       WHERE (na%vd >= 360._dp)
          na%vd = na%vd - 360._dp
       END WHERE
    END IF
    na%vd = na%vd * deg2rad

  END SUBROUTINE radians
  !------------------------------------------------------------------

  !------------------------------------------------------------------
  SUBROUTINE restrict_radians(arraylon,arraylat)

    implicit none

    ! i/o
    real(dp), dimension(:), intent(inout) :: arraylon     
    real(dp), dimension(:), intent(inout) :: arraylat     

    !-----------------------------------------------------------------------
    !
    !     convert longitudes to 0,2pi interval
    !
    !-----------------------------------------------------------------------

    where (arraylon .gt. pi2)  arraylon = arraylon - pi2
    where (arraylon .lt. zero) arraylon = arraylon + pi2

    !-----------------------------------------------------------------------
    !
    !     make sure input latitude range is within the machine values
    !     for +/- pi/2 
    !
    !-----------------------------------------------------------------------

    where (arraylat >  pih) arraylat =  pih
    where (arraylat < -pih) arraylat = -pih

  end subroutine restrict_radians
  !------------------------------------------------------------------

  !------------------------------------------------------------------
  subroutine read_a2o_nml(iou, status, gi, gatm, goces, goceu, gocev          &
       , outfile, var, vec, i_t, g_t, o_t, lint, lover, nmlstr)

    ! i/o
    integer                     , intent(in)  :: iou     ! i/o unit
    integer                     , intent(out) :: status  ! file status
    type (geohybgrid)        , intent(out) :: gi      ! input grid
    type (geohybgrid)        , intent(out) :: gatm    ! out grid (atm.)
    type (geohybgrid)        , intent(out) :: goces   ! out grid (oce. scal)
    type (geohybgrid)        , intent(out) :: goceu   ! out grid (oce.u-vec)
    type (geohybgrid)        , intent(out) :: gocev   ! out grid (oce.v-vec)
    character(len=grd_maxstrlen), intent(out) :: outfile ! output filename
    character(len=100*grd_maxstrlen), intent(out) :: var ! variable string
    character(len=100*grd_maxstrlen), intent(out) :: vec ! vector string
    integer                     , intent(out) :: i_t(4)  ! input time control
    integer                     , intent(out) :: g_t(3)  ! grid time control
    integer                     , intent(out) :: o_t(3)  ! output time control
    logical                     , intent(out) :: lint    ! use input time?
    logical                     , intent(out) :: lover   ! overlap the corners
                                                         ! (on oceangrid)?
    type (narray), intent(out)                :: nmlstr  ! namelist as string

    ! local
    character(len=grd_maxstrlen)  :: infile        ! input filename
    character(grd_maxstrlen) :: atm_grdfile        ! atm. grid file (netcdf)
    character(grd_maxstrlen) :: oces_grdfile       ! ocean grid file (netcdf)
    character(grd_maxstrlen) :: oceu_grdfile       ! ocean grid file (netcdf)
    character(grd_maxstrlen) :: ocev_grdfile       ! ocean grid file (netcdf)
    character(grd_maxstrlen) :: i_latm, i_lonm     ! input latitude/longitude
    character(grd_maxstrlen) :: i_lati, i_loni     ! input latitude/longitude
    character(grd_maxstrlen) :: c_latm, c_lonm     ! curve latitude/longitude
    character(grd_maxstrlen) :: c_lati, c_loni     ! curve latitude/longitude
    character(grd_maxstrlen) :: i_hyam, i_hybm     ! input hybrid coeff.(midl.)
    character(grd_maxstrlen) :: i_hyai, i_hybi     ! input hybrid coeff.(intf.)
    character(grd_maxstrlen) :: i_depthm, i_depthi ! input depth (midl./intf.)
    character(grd_maxstrlen) :: i_timei, i_timem   ! input time
    character(grd_maxstrlen) :: i_ps, i_p0         ! input sfc./ref. pressure

    character(grd_maxstrlen) :: g_atm_latm, g_atm_lonm     ! out lat/lon (mid)
    character(grd_maxstrlen) :: g_atm_lati, g_atm_loni     ! out lat/lon (int)
    character(grd_maxstrlen) :: g_atm_timei, g_atm_timem   ! out time
                                                                
    character(grd_maxstrlen) :: g_oces_latm, g_oces_lonm   ! out lat/lon (mid)
    character(grd_maxstrlen) :: g_oces_lati, g_oces_loni   ! out lat/lon (int)
    character(grd_maxstrlen) :: g_oces_timei, g_oces_timem ! out time
                                                                
    character(grd_maxstrlen) :: g_oceu_latm, g_oceu_lonm   ! out lat/lon (mid)
    character(grd_maxstrlen) :: g_oceu_lati, g_oceu_loni   ! out lat/lon (int)
    character(grd_maxstrlen) :: g_oceu_timei, g_oceu_timem ! out time
                                                                
    character(grd_maxstrlen) :: g_ocev_latm, g_ocev_lonm   ! out lat/lon (mid)
    character(grd_maxstrlen) :: g_ocev_lati, g_ocev_loni   ! out lat/lon (int)
    character(grd_maxstrlen) :: g_ocev_timei, g_ocev_timem ! out time
    logical                  :: input_time                 ! use in time info?
    logical                  :: ocean_over                 ! overlap oc. crnrs?

    character(len=*), parameter  :: substr='a2o_read_nml'

    ! helpers
    character(len=grd_maxstrlen) :: tstr
    character(len=8)             :: date
    character(len=10)            :: time
    character(len=5)             :: zone


    namelist /regrid/ infile, outfile                                         &
         ,atm_grdfile, oces_grdfile, oceu_grdfile, ocev_grdfile               &
         ,i_latm, i_lonm, i_lati, i_loni                                      &
         ,c_latm, c_lonm, c_lati, c_loni                                      &
         ,i_hyam, i_hybm, i_hyai, i_hybi                                      &
         ,i_depthm, i_depthi                                                  &
         ,i_ps, i_p0, i_timei, i_timem                                        &
         ,g_atm_latm, g_atm_lonm, g_atm_lati, g_atm_loni                      &
         ,g_atm_timei, g_atm_timem                                            &
         ,g_oces_latm, g_oces_lonm, g_oces_lati, g_oces_loni                  &
         ,g_oces_timei, g_oces_timem                                          &
         ,g_oceu_latm, g_oceu_lonm, g_oceu_lati, g_oceu_loni                  &
         ,g_oceu_timei, g_oceu_timem                                          &
         ,g_ocev_latm, g_ocev_lonm, g_ocev_lati, g_ocev_loni                  &
         ,g_ocev_timei, g_ocev_timem                                          &
         ,var, vec                                                            &
         ,i_t, g_t, o_t, input_time, ocean_over

    ! init (default values)
    infile       = ''
    outfile      = ''
    atm_grdfile  = ''
    oces_grdfile = ''
    oceu_grdfile = ''
    ocev_grdfile = ''
    i_latm   = ''
    i_lonm   = ''
    i_lati   = ''
    i_loni   = ''
    c_latm   = ''
    c_lonm   = ''
    c_lati   = ''
    c_loni   = ''
    i_hyam   = ''
    i_hybm   = ''
    i_hyai   = ''
    i_hybi   = ''
    i_depthm = ''
    i_depthi = ''
    i_timem  = ''
    i_timei  = ''
    i_ps     = ''
    i_p0     = ''
    g_atm_latm  = ''
    g_atm_lonm  = ''
    g_atm_lati  = ''
    g_atm_loni  = ''
    g_atm_timem = ''
    g_atm_timei = ''
    g_oces_latm = ''
    g_oces_lonm = ''
    g_oces_lati = ''
    g_oces_loni = ''
    g_oces_timem= ''
    g_oces_timei= ''
    g_oceu_latm = ''
    g_oceu_lonm = ''
    g_oceu_lati = ''
    g_oceu_loni = ''
    g_oceu_timem= ''
    g_oceu_timei= ''
    g_ocev_latm = ''
    g_ocev_lonm = ''
    g_ocev_lati = ''
    g_ocev_loni = ''
    g_ocev_timem= ''
    g_ocev_timei= ''
    var = ''
    vec = ''
    i_t(:) = (/ 1,1,0,0 /)
    g_t(:) = (/ 1,1,0 /)
    o_t(:) = (/ 1,1,0 /)
    input_time = .true.  ! default: use input time info for output
    ocean_over = .false. ! default: no overlapping ocean grid in long. dir.

    call init_geohybgrid(gi)
    call init_geohybgrid(gatm)
    call init_geohybgrid(goces)
    call init_geohybgrid(goceu)
    call init_geohybgrid(gocev)

    read(iou, nml=regrid, iostat=status)
    if (status /=0) return

    call init_narray(nmlstr)
    call addline_cnarray(nmlstr, 'a2o version '//trim(a2overs))
    call date_and_time(date, time, zone)
    call addline_cnarray(nmlstr, 'gt_date: '//trim(date))
    call addline_cnarray(nmlstr, 'gt_time: '//trim(time)//trim(zone))
    call addline_cnarray(nmlstr, 'namelist:')
    call addline_cnarray(nmlstr, '----------------------------------')
    call addline_cnarray(nmlstr, 'infile  = '''//trim(infile)//''',')
    call addline_cnarray(nmlstr, 'outfile = '''//trim(outfile)//''',')
    call addline_cnarray(nmlstr, 'atm_grdfile = '''//trim(atm_grdfile)//''',')
    call addline_cnarray(nmlstr, 'oces_grdfile = '''//trim(oces_grdfile)//    &
         ''',')
    call addline_cnarray(nmlstr, 'oceu_grdfile = '''//trim(oceu_grdfile)//    &
         ''',')
    call addline_cnarray(nmlstr, 'ocev_grdfile = '''//trim(ocev_grdfile)//    &
         ''',')
    call addline_cnarray(nmlstr, 'i_latm = '''//trim(i_latm)//''',')
    call addline_cnarray(nmlstr, 'i_lati = '''//trim(i_lati)//''',')
    call addline_cnarray(nmlstr, 'i_lonm = '''//trim(i_lonm)//''',')
    call addline_cnarray(nmlstr, 'i_loni = '''//trim(i_loni)//''',')
    call addline_cnarray(nmlstr, 'c_latm = '''//trim(c_latm)//''',')
    call addline_cnarray(nmlstr, 'c_lati = '''//trim(c_lati)//''',')
    call addline_cnarray(nmlstr, 'c_lonm = '''//trim(c_lonm)//''',')
    call addline_cnarray(nmlstr, 'c_loni = '''//trim(c_loni)//''',')
    call addline_cnarray(nmlstr, 'i_hyam = '''//trim(i_hyam)//''',')
    call addline_cnarray(nmlstr, 'i_hyai = '''//trim(i_hyai)//''',')
    call addline_cnarray(nmlstr, 'i_hybm = '''//trim(i_hybm)//''',')
    call addline_cnarray(nmlstr, 'i_hybi = '''//trim(i_hybi)//''',')
    call addline_cnarray(nmlstr, 'i_depthm = '''//trim(i_depthm)//''',')
    call addline_cnarray(nmlstr, 'i_depthi = '''//trim(i_depthi)//''',')
    call addline_cnarray(nmlstr, 'i_timem = '''//trim(i_timem)//''',')
    call addline_cnarray(nmlstr, 'i_timei = '''//trim(i_timei)//''',')
    call addline_cnarray(nmlstr, 'i_ps = '''//trim(i_ps)//''',')
    call addline_cnarray(nmlstr, 'i_p0 = '''//trim(i_p0)//''',')
    call addline_cnarray(nmlstr, 'g_atm_latm = '''//trim(g_atm_latm)//''',')
    call addline_cnarray(nmlstr, 'g_atm_lati = '''//trim(g_atm_lati)//''',')
    call addline_cnarray(nmlstr, 'g_atm_lonm = '''//trim(g_atm_lonm)//''',')
    call addline_cnarray(nmlstr, 'g_atm_loni = '''//trim(g_atm_loni)//''',')
    call addline_cnarray(nmlstr, 'g_atm_timem = '''//trim(g_atm_timem)//''',')
    call addline_cnarray(nmlstr, 'g_atm_timei = '''//trim(g_atm_timei)//''',')
    call addline_cnarray(nmlstr, 'g_oces_latm = '''//trim(g_oces_latm)//''',')
    call addline_cnarray(nmlstr, 'g_oces_lati = '''//trim(g_oces_lati)//''',')
    call addline_cnarray(nmlstr, 'g_oces_lonm = '''//trim(g_oces_lonm)//''',')
    call addline_cnarray(nmlstr, 'g_oces_loni = '''//trim(g_oces_loni)//''',')
    call addline_cnarray(nmlstr, 'g_oces_timem = '''//trim(g_oces_timem)//    &
         ''',')
    call addline_cnarray(nmlstr, 'g_oces_timei = '''//trim(g_oces_timei)//    &
         ''',')
    call addline_cnarray(nmlstr, 'g_oceu_latm = '''//trim(g_oceu_latm)//''',')
    call addline_cnarray(nmlstr, 'g_oceu_lati = '''//trim(g_oceu_lati)//''',')
    call addline_cnarray(nmlstr, 'g_oceu_lonm = '''//trim(g_oceu_lonm)//''',')
    call addline_cnarray(nmlstr, 'g_oceu_loni = '''//trim(g_oceu_loni)//''',')
    call addline_cnarray(nmlstr, 'g_oceu_timem = '''//trim(g_oceu_timem)//    &
         ''',')
    call addline_cnarray(nmlstr, 'g_oceu_timei = '''//trim(g_oceu_timei)//    &
         ''',')
    call addline_cnarray(nmlstr, 'g_ocev_latm = '''//trim(g_ocev_latm)//''',')
    call addline_cnarray(nmlstr, 'g_ocev_lati = '''//trim(g_ocev_lati)//''',')
    call addline_cnarray(nmlstr, 'g_ocev_lonm = '''//trim(g_ocev_lonm)//''',')
    call addline_cnarray(nmlstr, 'g_ocev_loni = '''//trim(g_ocev_loni)//''',')
    call addline_cnarray(nmlstr, 'g_ocev_timem = '''//trim(g_ocev_timem)//    &
         ''',')
    call addline_cnarray(nmlstr, 'g_ocev_timei = '''//trim(g_ocev_timei)//    &
         ''',')
    call addline_cnarray(nmlstr, 'var     = '''//trim(var)//''',')
    call addline_cnarray(nmlstr, 'vec     = '''//trim(vec)//''',')
    !
    write(tstr,'(4(i6,",",1x))') i_t
    call addline_cnarray(nmlstr, 'i_t     = '//trim(tstr))
    write(tstr,'(3(i6,",",1x))') g_t
    call addline_cnarray(nmlstr, 'g_t     = '//trim(tstr))
    write(tstr,'(3(i6,",",1x))') o_t
    call addline_cnarray(nmlstr, 'o_t     = '//trim(tstr))
    !
    if (input_time) then
       call addline_cnarray(nmlstr, 'input_time= .true.')
    else
       call addline_cnarray(nmlstr, 'input_time= .false.')
    end if
    !
    if (ocean_over) then
       call addline_cnarray(nmlstr, 'ocean_over= .true.')
    else
       call addline_cnarray(nmlstr, 'ocean_over= .false.')
    end if
    !
    call addline_cnarray(nmlstr, '----------------------------------')

    gi%file       = trim(infile)
    gi%latm%name  = trim(i_latm)
    gi%lati%name  = trim(i_lati)
    gi%lonm%name  = trim(i_lonm)
    gi%loni%name  = trim(i_loni)
    gi%clatm%name = trim(c_latm)
    gi%clati%name = trim(c_lati)
    gi%clonm%name = trim(c_lonm)
    gi%cloni%name = trim(c_loni)
    if (trim(i_hyam) /= '') then
       gi%hyam%name = trim(i_hyam)
    else
       gi%hyam%name = trim(i_depthm)
    end if
    if (trim(i_hyai) /= '') then
       gi%hyai%name = trim(i_hyai)
    else
       gi%hyai%name = trim(i_depthi)
    end if
    gi%hybm%name  = trim(i_hybm)
    gi%hybi%name  = trim(i_hybi)
    gi%timei%name = trim(i_timei)
    gi%timem%name = trim(i_timem)
    gi%ps%name    = trim(i_ps)
    gi%p0%name    = trim(i_p0)

    gatm%file       = trim(atm_grdfile)
    gatm%latm%name  = trim(g_atm_latm)
    gatm%lati%name  = trim(g_atm_lati)
    gatm%lonm%name  = trim(g_atm_lonm)
    gatm%loni%name  = trim(g_atm_loni)
    gatm%timei%name = trim(g_atm_timei)
    gatm%timem%name = trim(g_atm_timem)

    goces%file       = trim(oces_grdfile)
    goces%latm%name  = trim(g_oces_latm)
    goces%lati%name  = trim(g_oces_lati)
    goces%lonm%name  = trim(g_oces_lonm)
    goces%loni%name  = trim(g_oces_loni)
    goces%timei%name = trim(g_oces_timei)
    goces%timem%name = trim(g_oces_timem)

    goceu%file       = trim(oceu_grdfile)
    goceu%latm%name  = trim(g_oceu_latm)
    goceu%lati%name  = trim(g_oceu_lati)
    goceu%lonm%name  = trim(g_oceu_lonm)
    goceu%loni%name  = trim(g_oceu_loni)
    goceu%timei%name = trim(g_oceu_timei)
    goceu%timem%name = trim(g_oceu_timem)

    gocev%file       = trim(ocev_grdfile)
    gocev%latm%name  = trim(g_ocev_latm)
    gocev%lati%name  = trim(g_ocev_lati)
    gocev%lonm%name  = trim(g_ocev_lonm)
    gocev%loni%name  = trim(g_ocev_loni)
    gocev%timei%name = trim(g_ocev_timei)
    gocev%timem%name = trim(g_ocev_timem)

    lint = input_time
    lover= ocean_over

  end subroutine read_a2o_nml
  !------------------------------------------------------------------

  !------------------------------------------------------------------
  subroutine addline_cnarray(na, str)

    implicit none

    ! i/o
    type (narray),    intent(inout) :: na
    character(len=*), intent(in)    :: str

    ! local
    character(len=*), parameter :: substr = 'addline_cnarray'
    integer :: n, m ! string length
    integer :: i ! counter
    type (narray) :: nas

    ! local
    if (na%type == vtype_undef) then
       ! new
       n = len_trim(str)
       call init_narray(na, 1, (/ n+1 /), vtype_char)   ! +1 for 'newline'
       do i=1, n
          na%vc(i) = str(i:i)
       end do
       na%vc(n+1) = achar(10)
    else
       if (na%type == vtype_char) then
          ! save old
          call init_narray(nas)
          call copy_narray(nas, na)
          n = len_trim(str)
          m = size(nas%vc)
          call init_narray(na, 1, (/ m+n+1 /), vtype_char ) ! +1 for 'newline'
          do i=1, m
             na%vc(i) = nas%vc(i)
          end do
          do i=1, n
             na%vc(m+i) = str(i:i)
          end do
          na%vc(n+m+1) = achar(10)
          call init_narray(nas)
       else
          call rgmsg(substr, rgmle, 'n-array must be of type char !')
       end if
    end if

  end subroutine addline_cnarray
  !------------------------------------------------------------------

  !------------------------------------------------------------------
  subroutine parse_vecstr(vec,xvar,yvar)

    implicit none

    character(len=*), intent(in)           :: vec       ! vector string
    character(len=grd_maxstrlen), pointer  :: xvar(:)   ! x vector names
    character(len=grd_maxstrlen), pointer  :: yvar(:)   ! y vector names

    integer                                :: nvec      ! number of vectors

    character(grd_maxstrlen), allocatable  :: lvec(:)   ! list of vector infos

    integer                                :: status
    integer                                :: i
    integer                                :: idx1, idx2
    integer                                :: idxl, idxr
    integer                                :: len

    character(len=*), parameter :: substr = 'parse_vecstr'

    ! count variables
    nvec = 1
    len = len_trim(vec)
    do i=1, len-1        ! last ';' optional !
       if (vec(i:i) == ';') nvec = nvec + 1
    end do
    ! allocate space
    allocate(xvar(nvec), stat=status)
    call errmsg(substr,status,1)
    allocate(yvar(nvec), stat=status)
    call errmsg(substr,status,2)
    allocate(lvec(nvec), stat=status)
    call errmsg(substr,status,3)
    !
    !parsing
    idx1 = 1
    idx2 = index(trim(vec), ';')
    do i=1, nvec-1
       lvec(i) = vec(idx1:idx2-1)
       idx1 = idx2+1
       idx2 = idx2+index(vec(idx1:),';')
    end do
    if (idx2 == (idx1-1)) idx2 = len+1    ! last ';' optional !
    lvec(nvec) = vec(idx1:idx2-1)
    ! fill ovec, ivec, normalize_opt, rgt, rgtstr
    do i=1, nvec  ! loop over all vec variables
       ! get indices
       len  = len_trim(lvec(i))

       idx1 = index(lvec(i),',')

       ! check for x variable name
       idxl = 1
       idxr = idx1-1
       xvar(i) = lvec(i)(idxl:idxr)
       if (trim(xvar(i)) /= '') then
          call rgmsg(substr, rgmlvl, 'x-vectorcomponent ok')
       else
          call rgmsg(substr, rgmle, 'error: empty x-vectorcomponent!', .false.)
       end if

       ! check for y variable name
       idxl = idx1+1
       idxr = len
       yvar(i) = lvec(i)(idxl:idxr)
       if (trim(yvar(i)) /= '') then
          call rgmsg(substr, rgmlvl, 'y-vectorcomponent ok')
       else
          call rgmsg(substr, rgmle, 'error: empty y-vectorcomponent!', .false.)
       end if

    end do  ! loop over all vec variables

    ! clean
    if (allocated(lvec)) then
       deallocate(lvec, stat=status)
       call errmsg(substr,status,3)
    end if

  end subroutine parse_vecstr
  !------------------------------------------------------------------

  !------------------------------------------------------------------
  subroutine parse_varstr(var, ivar, ovar, sourcegrid, sourcegridstr,         &
       destgrid, destgridstr, mapping, mapstr , normalizing, normstr, file)

    implicit none

    ! i/o
    character(len=*), intent(in)           :: var       ! variable string
    character(len=grd_maxstrlen), pointer  :: ivar(:)   ! input variable names
    character(len=grd_maxstrlen), pointer  :: ovar(:)   ! output variabel names
    integer                     , pointer  :: sourcegrid(:)  ! sourcegrid
    integer                     , pointer  :: destgrid(:)    ! destgrid
    integer                     , pointer  :: mapping(:)     ! regridding opt.
    integer                     , pointer  :: normalizing(:) ! normalizing opt.
    character(len=*)            , pointer  :: sourcegridstr(:) ! source string
    character(len=*)            , pointer  :: destgridstr(:)   ! dest string
    character(len=*)            , pointer  :: mapstr(:)        ! regridding str
    character(len=*)            , pointer  :: normstr(:)       ! norm string
    character(len=*)            , optional :: file             ! filename

    ! local
    character(len=*), parameter :: substr = 'parse_varstr'
    type (ncvar), dimension(:), pointer   :: kvar      ! variables from file
    integer                               :: nvar      ! number of variables
    integer                               :: status
    integer                               :: i,j
    integer                               :: idx1, idx2, idx3
    integer                               :: idxl, idxr
    character(grd_maxstrlen), allocatable :: lvar(:)   ! list of variable infos
    integer                               :: len
    character(grd_maxstrlen)              :: infile    ! filename
    character(grd_maxstrlen), allocatable :: optstr(:,:)
    character(grd_maxstrlen)              :: partstr

    ! init
    if (present(file)) then
       infile = trim(file)
    else
       infile = ''
    end if
    nullify(kvar)

    !  count number of variables
    if ((trim(var) == '').or.(trim(var) == ';')) then   ! empty var list
       ! no variables given -> scan for all
       nvar = 0
       call rgmsg(substr, rgmli, 'no variables listed !')
       if (trim(infile) == '') then
          call rgmsg(substr, rgmle,'no input file given !')
       else
          call rgmsg(substr, rgmlvm, &
               'scanning input file '''//trim(infile)//''' ...')
          call scan_ncvar(kvar, file=trim(infile))
          nvar = size(kvar)
          allocate(ivar(nvar), stat=status)
          call errmsg(substr,status,1)
          allocate(ovar(nvar), stat=status)
          call errmsg(substr,status,2)
          allocate(sourcegrid(nvar), stat=status)
          call errmsg(substr,status,3)
          allocate(sourcegridstr(nvar), stat=status)
          call errmsg(substr,status,4)
          allocate(destgrid(nvar), stat=status)
          call errmsg(substr,status,3)
          allocate(destgridstr(nvar), stat=status)
          call errmsg(substr,status,4)
          allocate(mapping(nvar), stat=status)
          call errmsg(substr,status,5)
          allocate(mapstr(nvar), stat=status)
          call errmsg(substr,status,6)
          allocate(normalizing(nvar), stat=status)
          call errmsg(substr,status,7)
          allocate(normstr(nvar), stat=status)
          call errmsg(substr,status,8)
          allocate(optstr(nvar,4), stat=status)
          call errmsg(substr,status,9)
          !
          ! init
          sourcegridstr(:) = 'default: atm'
          sourcegrid(:) = atm
          destgridstr(:) = 'default: oces'
          destgrid(:) = oces
          mapstr(:) = 'default: conservative'
          mapping(:) = map_type_conserv
          normstr(:) = 'default: none'
          normalizing(:) = norm_opt_none
          optstr(:,:) = ''
          !
          do i=1, nvar
             ivar(i) = trim(kvar(i)%name)
             ovar(i) = trim(kvar(i)%name)
             !! try to get rgtstr, rgt from attributes
             !rgtstr(i) = '...'
             !do j=1, kvar(i)%natts
             !   if (trim(kvar(i)%att(j)%name) == 'rg_type') then
             !      rgtstr(i) = string(kvar(i)%att(j)%dat%vc)
             !   end if
             !end do
          end do
          ! clean up
          do i=1, nvar
             call init_ncvar(kvar(i))
          end do
          deallocate(kvar, stat=status)
          call errmsg(substr,status,10)
          nullify(kvar)
       end if

    else     ! var list not empty

       ! count variables
       nvar = 1
       len = len_trim(var)
       do i=1, len-1        ! last ';' optional !
          if (var(i:i) == ';') nvar = nvar + 1
       end do
       ! allocate space
       allocate(ivar(nvar), stat=status)
       call errmsg(substr,status,11)
       allocate(ovar(nvar), stat=status)
       call errmsg(substr,status,12)
       allocate(sourcegrid(nvar), stat=status)
       call errmsg(substr,status,13)
       allocate(sourcegridstr(nvar), stat=status)
       call errmsg(substr,status,14)
       allocate(destgrid(nvar), stat=status)
       call errmsg(substr,status,15)
       allocate(destgridstr(nvar), stat=status)
       call errmsg(substr,status,16)
       allocate(mapping(nvar), stat=status)
       call errmsg(substr,status,17)
       allocate(mapstr(nvar), stat=status)
       call errmsg(substr,status,18)
       allocate(normalizing(nvar), stat=status)
       call errmsg(substr,status,19)
       allocate(normstr(nvar), stat=status)
       call errmsg(substr,status,20)
       allocate(lvar(nvar), stat=status)
       call errmsg(substr,status,21)
       allocate(optstr(nvar,4), stat=status)
       call errmsg(substr,status,22)
       !
       ! init
       sourcegridstr(:) = 'default: atm'
       sourcegrid(:) = atm
       destgridstr(:) = 'default: oces'
       destgrid(:) = oces
       mapstr(:) = 'default: conservative'
       mapping(:) = map_type_conserv
       normstr(:) = 'default: none'
       normalizing(:) = norm_opt_none
       optstr(:,:) = ''
       !
       ! parsing
       idx1 = 1
       idx2 = index(trim(var), ';')
       do i=1, nvar-1
          lvar(i) = var(idx1:idx2-1)
          idx1 = idx2+1
          idx2 = idx2+index(var(idx1:),';')
       end do
       if (idx2 == (idx1-1)) idx2 = len+1    ! last ';' optional !
       lvar(nvar) = var(idx1:idx2-1)
       ! fill ovar, ivar, normalize_opt, rgt, rgtstr
       do i=1, nvar  ! loop over all variables
          ! get indices
          len  = len_trim(lvar(i))

          idx1 = index(lvar(i),'=')
          if (idx1 == 0) idx1 = len+1

          idx2 = index(lvar(i),':')
          if (idx2 == 0) idx2 = len+1

          idx3 = index(lvar(i),',')
          if (idx3 == 0) idx3 = len+1

          ! check for destination variable name
          idxl = 1
          idxr = min(idx1, idx2)-1
          ovar(i) = lvar(i)(idxl:idxr)

          ! check for source variable name
          idxl = idx1+1
          idxr = min(idx2, idx3)-1
          if (idxl > len) idxl = 1
          if (idxr > len) idxr = len
          ivar(i) = lvar(i)(idxl:idxr)

          ! get substring sourcegrid, destinationgrid, mapping method 
          ! and normalizing option
          idxl = idx2+1
          if (idxl < len) then
             partstr = lvar(i)
             do j=1,4
                partstr = partstr(idxl:len)
                len = len_trim(partstr)
                idx3 = index(partstr,',')
                if (idx3 == 0) idx3 = len
                idxr = idx3+1
                read(partstr(1:idxr-1),*) optstr(i,j)
                idxl = idxr
                if (idxl > len) exit
             end do

          end if
       end do  ! loop over all variables

       ! clean
       if (allocated(lvar)) then
          deallocate(lvar, stat=status)
          call errmsg(substr,status,23)
       end if
    end if   ! var list empty or not


    do i=1, nvar
       select case (optstr(i,1))
       case ('atm','ATM')
          sourcegridstr(i)='atm'
          sourcegrid(i)=atm
       case ('oces','OCES')
          sourcegridstr(i)='oces'
          sourcegrid(i)=oces
       case ('oceu','OCEU')
          sourcegridstr(i)='oceu'
          sourcegrid(i)=oceu
       case ('ocev','OCEV')
          sourcegridstr(i)='ocev'
          sourcegrid(i)=ocev
       case default
          call rgmsg(substr, rgmlw, &
               'unknown source grid of varaible '''//trim(ivar(i))//'''')
          call rgmsg(substr, rgmlwc, &
               'using default value: atm!')
       end select

       select case (optstr(i,2))
       case ('atm','ATM')
          destgridstr(i)='atm'
          destgrid(i)=atm
       case ('oces','OCES')
          destgridstr(i)='oces'
          destgrid(i)=oces
       case ('oceu','OCEU')
          destgridstr(i)='oceu'
          destgrid(i)=oceu
       case ('ocev','OCEV')
          destgridstr(i)='ocev'
          destgrid(i)=ocev
       case default
          call rgmsg(substr, rgmlw, &
               'unknown destination grid of varaible '''//trim(ivar(i))//'''')
          call rgmsg(substr, rgmlwc, &
               'using default value: oces!')
       end select

       do j=1,4
          select case (optstr(i,j))
          case ('conservative','CONSERVATIVE')
             mapstr(i)='conservative'
             mapping(i)=map_type_conserv
          case ('bilinear','BILINEAR')
             mapstr(i)='bilinear'
             mapping(i)=map_type_bilinear
          case ('bicubic','BICUBIC')
             mapstr(i)='bicubic'
             mapping(i)=map_type_bicubic
          case ('distwgt','DISTWGT')
             mapstr(i)='distwgt'
             mapping(i)=map_type_distwgt
          case ('none','NONE')
             normstr(i)='none'
             normalizing(i)=norm_opt_none
          case ('dest','DEST')
             normstr(i)='dest'
             normalizing(i)=norm_opt_dstarea
          case ('frac','FRAC')
             normstr(i)='frac'
             normalizing(i)=norm_opt_frcarea
          end select
       end do

       select case (mapstr(i))
       case ('conservative','bilinear','bicubic','distwgt')
       case default
          call rgmsg(substr, rgmlw                                            &
               , 'unknown mapping method option of varaible '''               &
               //trim(ivar(i))//'''')
          call rgmsg(substr, rgmlwc                                           &
               , 'using default value: conservative!')
       end select

       select case (normstr(i))
       case ('none','frac','dest')
       case default
          call rgmsg(substr, rgmlw                                            &
               , 'unknown normalizing method option of varaible '''           &
               //trim(ivar(i))//'''')
          call rgmsg(substr, rgmlwc                                           &
               , 'using default value: none!')
       end select

    end do

    ! clean
    if (allocated(optstr)) then
       deallocate(optstr, stat=status)
       call errmsg(substr,status,21)
    end if

  end subroutine parse_varstr
  !------------------------------------------------------------------

  !------------------------------------------------------------------
  subroutine rotate_2_ne(u_i,v_j,lon,lat,ix,iy)
    !----------------------------------------------------------------------
    ! this subroutine is from mpio-m (with additional changes)
    !----------------------------------------------------------------------
    !
    !     rotation of vectors: in ocean models with rotated grids velocity
    !     vectors are given in the direction of grid lines and rows. they 
    !     have to be rotated in latitudinal and longitudinal direction.
    !
    !     note: this routine assumes positive meridional flow for a flow
    !           from grid point(i,j) to grid point(i,j+1) and positive 
    !           zonal flow for a flow from grid point(i,j) to point(i+1,j).
    !           this is not the case for mpi-om!
    !
    !           if this routine is used to rotate data of mpi-om, the 
    !           logical change_sign_v needs to be true.
    !j. jungclaus: 22.01.04:
    !note here for the coupling fields u-i,v_j are on the non-verlapping
    ! (ie-2) grid, furthermore, the velocity fields were previously
    ! interpolated onto the scalar points !
    !
    !h.haak: 07.10.2005 vectorisation and omp directives      
    !----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !     local variables
    !-----------------------------------------------------------------------

    integer, intent(in) :: ix,iy
    real(dp), intent(in) :: lon(ix,iy),lat(ix,iy)
    real(dp), intent(inout) :: u_i(ix,iy),v_j(ix,iy) ! vec components
    real(dp) :: u_lon(ix,iy),v_lat(ix,iy)            ! vec components in n-e 

    real(dp) :: dlat_i, dlat_j,dlon_i,dlon_j,dist_i,dist_j
    real(dp) :: lat_factor

    integer :: i,j,ip1,im1,jp1,jm1

    logical :: change_sign_u,change_sign_v
    !-----------------------------------------------------------------------
    ! specification whether change in sign is needed for the input arrays
    change_sign_u=.false.
    change_sign_v=.true.

    ! initialization
    ! --------------
    do i = 1, ix
       do j = 1, iy
          v_lat(i,j) = 0.0_dp 
          u_lon(i,j) = 0.0_dp 
       end do
    end do

    if (change_sign_u) then
       do i = 1, ix
          do j = 1, iy
             u_i(i,j)=u_i(i,j)*(-1.0_dp)
          end do
       end do
    endif

    if (change_sign_v) then
       do i = 1, ix
          do j = 1, iy
             v_j(i,j)=v_j(i,j)*(-1.0_dp)
          end do
       end do
    endif


    ! rotation
    ! --------
    do i = 1, ix
       do j = 1, iy

          ip1 = i + 1
          im1 = i - 1
          jp1 = j + 1
          jm1 = j - 1

          if (ip1 > ix) ip1 = ip1 - ix ! the 0-meridian
          if (im1 < 1 ) im1 = ix
          if (jp1 > iy) then           ! treatment of the last..
             jp1 = j
          endif
          if (jm1 < 1 ) then ! .. and the fist grid-row
             jm1 = j
          endif

          ! difference in latitudes
          dlat_i = lat(ip1,j) - lat(im1,j)
          dlat_j = lat(i,jp1) - lat(i,jm1)

          ! difference in longitudes                  
          dlon_i = lon(ip1,j) - lon(im1,j)
          if (dlon_i >   pi)  dlon_i = dlon_i - (2._dp*pi)
          if (dlon_i < (-pi)) dlon_i = dlon_i + (2._dp*pi)
          dlon_j = lon(i,jp1) - lon(i,jm1)
          if (dlon_j >   pi)  dlon_j = dlon_j - (2._dp*pi)
          if (dlon_j < (-pi)) dlon_j = dlon_j + (2._dp*pi)

          lat_factor = cos(lat(i,j))
          dlon_i = dlon_i * lat_factor
          dlon_j = dlon_j * lat_factor

          ! projection by scalar product
          ! ----------------------------
          u_lon(i,j) = u_i(i,j)*dlon_i + v_j(i,j)*dlat_i
          v_lat(i,j) = u_i(i,j)*dlon_j + v_j(i,j)*dlat_j

          dist_i = sqrt(dlon_i**2+dlat_i**2)
          dist_j = sqrt(dlon_j**2+dlat_j**2)
          if (dist_i /= 0.0_dp .and. dist_j /= 0.0_dp) then
             u_lon(i,j) = u_lon(i,j)/dist_i
             v_lat(i,j) = v_lat(i,j)/dist_j
          else
             u_lon(i,j) = 0.0_dp 
             v_lat(i,j) = 0.0_dp  
          endif

       end do
    end do

    ! write back to input field!
    do i=1,ix
       do j=1,iy
          u_i(i,j)=u_lon(i,j)
          v_j(i,j)=v_lat(i,j)
          !u_i(i,j)=u_lon(i,j)*weto_g(i+1,j,1)
          !v_j(i,j)=v_lat(i,j)*weto_g(i+1,j,1)
       end do
    end do

  end subroutine rotate_2_ne
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  subroutine rotate_u(fieldu,fieldv,lon,lat,ii,jj)
    !----------------------------------------------------------------------
    ! this subroutine is from mpio-m (with additional changes)
    !----------------------------------------------------------------------

    integer, intent(in)     :: ii,jj
    real(dp), intent(inout) :: fieldu(ii,jj),fieldv(ii,jj)
    real(dp), intent(in)    :: lon(ii,jj),lat(ii,jj)

    integer  :: i,j
    real(dp) :: uin,vin,uout,vout
    real(dp) :: deltxx,deltyy,deltxy,deltyx
    real(dp) :: proms
    real(dp) :: absin,absout

    integer  :: ip,im,jp,jm

    do i=1,ii 
       do j=1,jj 
          uin=fieldu(i,j)
          vin=fieldv(i,j)
          uout=0.0_dp
          vout=0.0_dp

          ip = i+1
          im = i-1
          jp = j+1
          jm = j-1
          if (ip > ii) ip = 1
          if (im == 0) im = ii
          if (jp > jj) jp = 1
          if (jm == 0) jm = jj

          !deltxx=lon(ip,j) - lon(i,j)
          !deltyy=lat(ip,jm)-lat(ip,j)
          !deltxy=lat(ip,j) - lat(i,j)
          !deltyx=lon(ip,jm)-lon(ip,j)

          deltxx=lon(ip,j) - lon(i,j)
          deltyy=lat(i,jp)-lat(i,j)
          deltxy=lat(ip,j) - lat(i,j)
          deltyx=lon(i,jp)-lon(i,j)

          if (deltxx .gt. pi) deltxx=deltxx-2._dp*pi
          if (deltxx .lt. -pi) deltxx=deltxx+2._dp*pi
          if (deltyx .lt. -pi) deltyx=deltyx+2._dp*pi
          if (deltyx .gt. pi) deltyx=deltyx-2._dp*pi

          !proms=cos(lat(ip,j))
          proms=cos(lat(i,j))
          deltxx=proms*deltxx
          deltyx=proms*deltyx

          !:: projection by scalar product
          uout= uin*deltxx + vin*deltxy
          vout= uin*deltyx + vin*deltyy
          uout=uout/sqrt(deltxx**2+deltxy**2)
          vout=vout/sqrt(deltyy**2+deltyx**2)
          !:: correct skale
          absin=sqrt(uin**2+vin**2)
          absout=sqrt(uout**2+vout**2)
          if (absout.gt.1.e-10_dp) then
             uout=uout*absin/absout
          else
             uout=0.0_dp
          endif
          fieldu(i,j)=uout
       enddo
    enddo
  end subroutine rotate_u
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  subroutine rotate_v(fieldu,fieldv,lon,lat,ii,jj)
    !----------------------------------------------------------------------
    ! this subroutine is from mpio-m (with additional changes)
    !----------------------------------------------------------------------

    integer, intent(in) :: ii,jj
    real(dp), intent(inout) :: fieldu(ii,jj),fieldv(ii,jj)
    real(dp), intent(in) :: lon(ii,jj),lat(ii,jj)

    integer :: i,j
    real(dp) :: uin,vin,uout,vout
    real(dp) :: deltxx,deltyy,deltxy,deltyx
    real(dp) :: proms
    real(dp) :: absin,absout

    integer :: ip,im,jp,jm

    do j=1,jj
       do i=1,ii
          uin=fieldu(i,j)
          vin=fieldv(i,j)
          uout=0.0_dp
          vout=0.0_dp

          ip = i+1
          im = i-1
          jp = j+1
          jm = j-1
          if (ip > ii) ip = 1
          if (im == 0) im = ii
          if (jp > jj) jp = 1
          if (jm == 0) jm = jj
          !DELTXX=lon(I,jp) - lon(im,jp)
          !DELTYY=lat(I,J) - lat(I,jp)
          !DELTXY=lat(I,jp) - lat(im,jp)
          !DELTYX=lon(I,J)-lon(I,jp)

          DELTXX=lon(ip,j) - lon(i,j)
          DELTYY=lat(I,jp) - lat(I,j)
          DELTXY=lat(ip,j) - lat(i,j)
          DELTYX=lon(I,jp)-lon(I,j)

          IF (DELTXX .GT. PI) DELTXX=DELTXX-2._dp*PI
          IF (DELTXX .LT. -PI) DELTXX=DELTXX+2._dp*PI
          IF (DELTYX .LT. -PI) DELTYX=DELTYX+2._dp*PI
          IF (DELTYX .GT. PI) DELTYX=DELTYX-2._dp*PI

          !PROMS=COS( lat(I,jp) )
          PROMS=COS( lat(I,j) )
          DELTXX=PROMS*DELTXX
          DELTYX=PROMS*DELTYX
          !:: PROJECTION BY SCALAR PRODUCT
          UOUT= UIN*DELTXX + VIN*DELTXY
          VOUT= UIN*DELTYX + VIN*DELTYY
          UOUT=UOUT/SQRT(DELTXX**2+DELTXY**2)
          VOUT=VOUT/SQRT(DELTYY**2+DELTYX**2)
          !:: CORRECT SKALE
          ABSIN=SQRT(UIN**2+VIN**2)
          ABSOUT=SQRT(UOUT**2+VOUT**2)
          IF (ABSOUT .GT. 1.E-10_dp) THEN
             VOUT=VOUT*ABSIN/ABSOUT
          ELSE
             VOUT=0.0_dp
          ENDIF
          FIELDV(I,J)=VOUT
       ENDDO
    ENDDO
  END SUBROUTINE ROTATE_V
  !-------------------------------------------------------------------------

END MODULE MESSY_A2O_GRIDTRAFO
